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
  use randomNumber_mod
  use controlVector_mod
  use gridVariableTransforms_mod
  use bMatrix_mod
  use humidityLimits_mod
  use localizationFunction_mod
  use varNameList_mod
  use codePrecision_mod
  use fileNames_mod
  use codTyp_mod
  use clib_interfaces_mod
  implicit none
  save
  private

  ! public types
  public :: struct_enkfInterpInfo

  ! public procedures
  public :: enkf_setupInterpInfo, enkf_LETKFanalyses, enkf_modifyAMSUBobsError
  public :: enkf_rejectHighLatIR, enkf_postProcess

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
            if (vnl_varLevelFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex)) == 'SF') then
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
            if (vnl_varLevelFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex)) == 'SF') then
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
  ! enkf_postProcess
  !----------------------------------------------------------------------
  subroutine enkf_postProcess(ensembleAnl, ensembleTrl, stateVectorHeightSfc)
    !
    !:Purpose:  Perform numerous post-processing steps to the ensemble
    !           produced by the LETKF algorithm.
    !
    implicit none

    ! Arguments
    type(struct_ens), pointer :: ensembleTrl
    type(struct_ens)          :: ensembleAnl
    type(struct_gsv)          :: stateVectorHeightSfc

    ! Locals
    integer                   :: ierr, nEns, dateStamp, datePrint, timePrint, imode, randomSeedRandomPert
    integer                   :: stepIndex, middleStepIndex, nulnam
    type(struct_hco), pointer :: hco_ens
    type(struct_vco), pointer :: vco_ens
    type(struct_gsv)          :: stateVectorMeanAnl, stateVectorMeanTrl
    type(struct_gsv)          :: stateVectorMeanInc
    type(struct_gsv)          :: stateVectorStdDevAnl, stateVectorStdDevAnlPert, stateVectorStdDevTrl
    type(struct_gsv)          :: stateVectorMeanIncSubSample
    type(struct_gsv)          :: stateVectorMeanAnlSubSample
    type(struct_gsv)          :: stateVectorMeanAnlSfcPres
    type(struct_gsv)          :: stateVectorMeanAnlSfcPresMpiGlb
    type(struct_ens)          :: ensembleTrlSubSample
    type(struct_ens)          :: ensembleAnlSubSample
    character(len=12)         :: etiketMean='', etiketStd=''
    character(len=256)        :: outFileName
    character(len=4), pointer :: varNames(:)

    integer, external         :: fnom, fclos, newdate

    ! Namelist variables
    integer  :: randomSeed           ! seed used for random perturbation additive inflation
    logical  :: writeSubSample       ! write sub-sample members for initializing medium-range fcsts
    real(8)  :: alphaRTPS            ! RTPS coefficient (between 0 and 1; 0 means no relaxation)
    real(8)  :: alphaRandomPert      ! Random perturbation additive inflation coeff (0->1)
    real(8)  :: alphaRandomPertSubSample ! Random perturbation additive inflation coeff for medium-range fcsts
    logical  :: imposeSaturationLimit  ! switch for choosing to impose saturation limit of humidity
    logical  :: imposeRttovHuLimits    ! switch for choosing to impose the RTTOV limits on humidity
    real(8)  :: weightRecenter         ! weight applied to EnVar recentering increment
    integer  :: numMembersToRecenter ! number of members that get recentered on EnVar analysis
    logical  :: useOptionTableRecenter ! use values in the optiontable file
    character(len=12) :: etiket0

    NAMELIST /NAMENSPOSTPROC/randomSeed, writeSubSample,  &
                             alphaRTPS, alphaRandomPert, alphaRandomPertSubSample,  &
                             imposeSaturationLimit, imposeRttovHuLimits,  &
                             weightRecenter, numMembersToRecenter, useOptionTableRecenter,  &
                             etiket0

    hco_ens => ens_getHco(ensembleAnl)
    vco_ens => ens_getVco(ensembleAnl)
    nEns = ens_getNumMembers(ensembleAnl)

    !- Setting default namelist variable values
    randomSeed            =  -999
    writeSubSample        = .false.
    alphaRTPS             =  0.0D0
    alphaRandomPert       =  0.0D0
    alphaRandomPertSubSample =  -1.0D0
    imposeSaturationLimit = .false.
    imposeRttovHuLimits   = .false.
    weightRecenter        = 0.0D0 ! means no recentering applied
    numMembersToRecenter  = -1    ! means all members recentered by default
    useOptionTableRecenter = .false.
    etiket0               = 'E26_0_0P'

    !- Read the namelist
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=namenspostproc, iostat=ierr)
    if ( ierr /= 0) call utl_abort('enkf_postProc: Error reading namelist')
    if ( mpi_myid == 0 ) write(*,nml=namenspostproc)
    ierr = fclos(nulnam)

    if (alphaRTPS < 0.0D0) alphaRTPS = 0.0D0
    if (alphaRandomPert < 0.0D0) alphaRandomPert = 0.0D0
    if (alphaRandomPertSubSample < 0.0D0) alphaRandomPertSubSample = 0.0D0
    if (numMembersToRecenter == -1) numMembersToRecenter = nEns ! default behaviour

    !- Allocate and compute ensemble mean Trl and Anl
    call gsv_allocate( stateVectorMeanTrl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_allocate( stateVectorMeanAnl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorMeanTrl)
    call gsv_zero(stateVectorMeanAnl)
    call ens_computeMean(ensembleTrl)
    call ens_computeMean(ensembleAnl)
    call ens_copyEnsMean(ensembleTrl, stateVectorMeanTrl)
    call ens_copyEnsMean(ensembleAnl, stateVectorMeanAnl)

    !- Allocate and compute ensemble spread stddev Trl and Anl (AnlPert computed later)
    call gsv_allocate( stateVectorStdDevTrl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_allocate( stateVectorStdDevAnl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeight_opt=.false., allocPressure_opt=.false. )
    call ens_computeStdDev(ensembleTrl)
    call ens_computeStdDev(ensembleAnl)
    call ens_copyEnsStdDev(ensembleTrl, stateVectorStdDevTrl)
    call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnl)

    !- Allocate and compute mean increment
    call gsv_allocate( stateVectorMeanInc, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorMeanInc)
    call gsv_copy(stateVectorMeanAnl, stateVectorMeanInc)
    call gsv_add(stateVectorMeanTrl, stateVectorMeanInc, scaleFactor_opt=-1.0D0)

    !- Apply RTPS, if requested
    if (alphaRTPS > 0.0D0) then
      call enkf_RTPS(ensembleAnl, ensembleTrl, stateVectorStdDevAnl, stateVectorStdDevTrl, stateVectorMeanAnl, alphaRTPS)
      ! recompute the analysis spread stddev
      call ens_computeStdDev(ensembleAnl)
      call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnl)
    end if

    !- Impose limits on humidity *before* random perturbations, if requested
    if (imposeSaturationLimit .or. imposeRttovHuLimits) then
      call tmg_start(102,'LETKF-imposeHulimits')
      if (mpi_myid == 0) write(*,*) ''
      if (mpi_myid == 0) write(*,*) 'midas-letkf: limits will be imposed on the humidity of analysis ensemble'
      if (mpi_myid == 0 .and. imposeSaturationLimit ) write(*,*) '              -> Saturation Limit'
      if (mpi_myid == 0 .and. imposeRttovHuLimits   ) write(*,*) '              -> Rttov Limit'
      if ( imposeSaturationLimit ) call qlim_saturationLimit(ensembleAnl)
      if ( imposeRttovHuLimits   ) call qlim_rttovLimit     (ensembleAnl)
      ! And recompute analysis mean
      call ens_computeMean(ensembleAnl)
      call ens_copyEnsMean(ensembleAnl, stateVectorMeanAnl)
      ! And recompute mean increment
      call gsv_copy(stateVectorMeanAnl, stateVectorMeanInc)
      call gsv_add(stateVectorMeanTrl, stateVectorMeanInc, scaleFactor_opt=-1.0D0)
      call tmg_stop(102)
    end if

    !- Recenter analysis ensemble on EnVar analysis
    if (weightRecenter > 0.0D0 .or. useOptionTableRecenter) then
      write(*,*) 'midas-letkf: Recenter analyses on EnVar analysis'
      call enkf_hybridRecentering(ensembleAnl, weightRecenter, useOptionTableRecenter, numMembersToRecenter)
      ! And recompute analysis mean
      call ens_computeMean(ensembleAnl)
      call ens_copyEnsMean(ensembleAnl, stateVectorMeanAnl)
      ! And recompute mean increment
      call gsv_copy(stateVectorMeanAnl, stateVectorMeanInc)
      call gsv_add(stateVectorMeanTrl, stateVectorMeanInc, scaleFactor_opt=-1.0D0)
      ! And recompute the analysis spread stddev
      call ens_computeStdDev(ensembleAnl)
      call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnl)
    end if

    !- If SubSample requested, copy sub-sample of analysis and trial members
    if (writeSubSample) then
      ! Copy sub-sampled analysis and trial ensemble members
      call enkf_selectSubSample(ensembleAnl, ensembleTrl,  &
                                ensembleAnlSubSample, ensembleTrlSubSample)

      ! Create subdirectory for outputting sub sample increments
      ierr = clib_mkdir_r('subspace')

      ! Allocate stateVectors to store and output sub-sampled ensemble mean analysis and increment
      call gsv_allocate( stateVectorMeanAnlSubSample, tim_nstepobsinc,  &
                         hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                         mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                         dataKind_opt=4, allocHeightSfc_opt=.true., &
                         allocHeight_opt=.false., allocPressure_opt=.false. )
      call gsv_zero(stateVectorMeanAnlSubSample)
      call gsv_allocate( stateVectorMeanIncSubSample, tim_nstepobsinc,  &
                         hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                         mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                         dataKind_opt=4, allocHeightSfc_opt=.true., &
                         allocHeight_opt=.false., allocPressure_opt=.false. )
      call gsv_zero(stateVectorMeanIncSubSample)

    end if

    !- Apply random additive inflation, if requested
    if (alphaRandomPert > 0.0D0) then
      ! If namelist value is -999, set random seed using the date (as in standard EnKF)
      if (randomSeed == -999) then
        imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
        dateStamp = tim_getDateStamp()
        ierr = newdate(dateStamp, datePrint, timePrint, imode)
        timePrint = timePrint/1000000
        datePrint =  datePrint*100 + timePrint
        ! Remove the year and add 9
        randomSeedRandomPert = 9 + datePrint - 1000000*(datePrint/1000000)
        write(*,*) 'midas-letkf: randomSeed for additive inflation set to ', randomSeedRandomPert
      else
        randomSeedRandomPert = randomSeed
      end if
      call tmg_start(101,'LETKF-randomPert')
      call enkf_addRandomPert(ensembleAnl, stateVectorMeanTrl, alphaRandomPert, randomSeedRandomPert)
      call tmg_stop(101)
    end if

    !- Recompute the analysis spread stddev after inflation and humidity limits
    call gsv_allocate( stateVectorStdDevAnlPert, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeight_opt=.false., allocPressure_opt=.false. )
    call ens_computeStdDev(ensembleAnl)
    call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnlPert)

    !- If SubSample requested, do remaining processing and output of sub-sampled members
    if (writeSubSample) then

      ! Apply random additive inflation to sub-sampled ensemble, if requested
      if (alphaRandomPertSubSample > 0.0D0) then
        ! If namelist value is -999, set random seed using the date (as in standard EnKF)
        if (randomSeed == -999) then
          imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
          dateStamp = tim_getDateStamp()
          ierr = newdate(dateStamp, datePrint, timePrint, imode)
          timePrint = timePrint/1000000
          datePrint =  datePrint*100 + timePrint
          ! Remove the year and add 9
          randomSeedRandomPert = 9 + datePrint - 1000000*(datePrint/1000000)
          write(*,*) 'midas-letkf: randomSeed for additive inflation set to ', randomSeedRandomPert
        else
          randomSeedRandomPert = randomSeed
        end if
        call tmg_start(101,'LETKF-randomPert')
        call enkf_addRandomPert(ensembleAnlSubSample, stateVectorMeanTrl,  &
                                alphaRandomPertSubSample, randomSeedRandomPert)
        call tmg_stop(101)
      end if

      ! Compute analysis mean of sub-sampled ensemble
      call ens_computeMean(ensembleAnlSubSample)

      ! Shift members to have same mean as full ensemble and impose humidity limits, if requested
      call ens_recenter(ensembleAnlSubSample, stateVectorMeanAnl,  &
                        recenteringCoeff_opt=1.0D0)

      ! Re-compute analysis mean of sub-sampled ensemble
      call ens_computeMean(ensembleAnlSubSample)
      call ens_copyEnsMean(ensembleAnlSubSample, stateVectorMeanAnlSubSample)

      ! And compute mean increment with respect to mean of full trial ensemble
      call gsv_copy(stateVectorMeanAnlSubSample, stateVectorMeanIncSubSample)
      call gsv_add(stateVectorMeanTrl, stateVectorMeanIncSubSample, scaleFactor_opt=-1.0D0)

    end if

    !
    !- Output everything
    !
    call tmg_start(4,'LETKF-writeOutput')


    !- Prepare stateVector with only MeanAnl surface pressure and surface height
    call gsv_allocate( stateVectorMeanAnlSfcPres, tim_nstepobsinc, hco_ens, vco_ens,   &
                       dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0'/) )
    call gsv_zero(stateVectorMeanAnlSfcPres)
    if (mpi_myid <= (nEns-1)) then
      call gsv_allocate( stateVectorMeanAnlSfcPresMpiGlb, tim_nstepobsinc, hco_ens, vco_ens,   &
                         dateStamp_opt=tim_getDateStamp(),  &
                         mpi_local_opt=.false., &
                         dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0'/) )
      call gsv_zero(stateVectorMeanAnlSfcPresMpiGlb)
    end if
    call gsv_copy(stateVectorMeanAnl, stateVectorMeanAnlSfcPres, allowMismatch_opt=.true.)
    call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorMeanAnlSfcPres)
    call gsv_transposeTilesToMpiGlobal(stateVectorMeanAnlSfcPresMpiGlb, stateVectorMeanAnlSfcPres)

    !- Output ens stddev and mean in trialrms, analrms and analpertrms files

    ! determine middle timestep for output of these files
    middleStepIndex = (tim_nstepobsinc + 1) / 2

    ! output trialmean, trialrms
    call enkf_getRmsEtiket(etiketMean, etiketStd, 'F', etiket0, nEns)
    call fln_ensTrlFileName(outFileName, '.', tim_getDateStamp())
    outFileName = trim(outFileName) // '_trialmean'
    call gsv_writeToFile(stateVectorMeanTrl, outFileName, trim(etiketMean),  &
                         typvar_opt='P', writeHeightSfc_opt=.false., numBits_opt=16,  &
                         stepIndex_opt=middleStepIndex, containsFullField_opt=.true.)
    call fln_ensTrlFileName(outFileName, '.', tim_getDateStamp())
    outFileName = trim(outFileName) // '_trialrms'
    call gsv_writeToFile(stateVectorStdDevTrl, outFileName, trim(etiketStd),  &
                         typvar_opt='P', writeHeightSfc_opt=.false., numBits_opt=16, &
                         stepIndex_opt=middleStepIndex, containsFullField_opt=.false.)

    ! output analmean, analrms
    call enkf_getRmsEtiket(etiketMean, etiketStd, 'A', etiket0, nEns)
    call fln_ensAnlFileName(outFileName, '.', tim_getDateStamp())
    outFileName = trim(outFileName) // '_analmean'
    call gsv_writeToFile(stateVectorMeanAnl, outFileName, trim(etiketMean),  &
                         typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                         stepIndex_opt=middleStepIndex, containsFullField_opt=.true.)
    call fln_ensAnlFileName(outFileName, '.', tim_getDateStamp())
    outFileName = trim(outFileName) // '_analrms'
    call gsv_writeToFile(stateVectorStdDevAnl, outFileName, trim(etiketStd),  &
                         typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                         stepIndex_opt=middleStepIndex, containsFullField_opt=.false.)

    ! output analpertmean, analpertrms
    call enkf_getRmsEtiket(etiketMean, etiketStd, 'P', etiket0, nEns)
    call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp() )
    outFileName = trim(outFileName) // '_analpertmean'
    call gsv_writeToFile(stateVectorMeanAnl, outFileName, trim(etiketMean),  &
                         typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                         stepIndex_opt=middleStepIndex, containsFullField_opt=.true.)
    call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp() )
    outFileName = trim(outFileName) // '_analpertrms'
    call gsv_writeToFile(stateVectorStdDevAnlPert, outFileName, trim(etiketStd),  &
                         typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                         stepIndex_opt=middleStepIndex, containsFullField_opt=.false.)

    !- Output the ensemble mean increment (include MeanAnl Psfc) and analysis

    ! convert transformed to model variables for ensemble mean of analysis and trial
    call gvt_transform(stateVectorMeanAnl,'AllTransformedToModel',allowOverWrite_opt=.true.)
    call gvt_transform(stateVectorMeanTrl,'AllTransformedToModel',allowOverWrite_opt=.true.)
    ! and recompute mean increment for converted model variables (e.g. VIS and PR)
    nullify(varNames)
    call gsv_varNamesList(varNames, stateVectorMeanAnl)
    call gsv_deallocate( stateVectorMeanInc )
    call gsv_allocate( stateVectorMeanInc, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles',  &
                       dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=varNames )
    call gsv_copy(stateVectorMeanAnl, stateVectorMeanInc)
    call gsv_add(stateVectorMeanTrl, stateVectorMeanInc, scaleFactor_opt=-1.0D0)
    deallocate(varNames)

    ! output ensemble mean increment
    call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp(), 0, ensFileNameSuffix_opt='inc' )
    do stepIndex = 1, tim_nstepobsinc
      call gsv_writeToFile(stateVectorMeanInc, outFileName, 'ENSMEAN_INC',  &
                           typvar_opt='R', writeHeightSfc_opt=.false., numBits_opt=16, &
                           stepIndex_opt=stepIndex, containsFullField_opt=.false.)
      call gsv_writeToFile(stateVectorMeanAnlSfcPres, outFileName, 'ENSMEAN_INC',  &
                           typvar_opt='A', writeHeightSfc_opt=.true., &
                           stepIndex_opt=stepIndex, containsFullField_opt=.true.)
    end do

    ! output ensemble mean analysis state
    call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp(), 0 )
    do stepIndex = 1, tim_nstepobsinc
      call gsv_writeToFile(stateVectorMeanAnl, outFileName, 'ENSMEAN_ANL',  &
                           typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                           stepIndex_opt=stepIndex, containsFullField_opt=.true.)
    end do

    !- Output all ensemble member analyses and increments
    ! convert transformed to model variables for analysis and trial ensembles
    call gvt_transform(ensembleAnl,'AllTransformedToModel',allowOverWrite_opt=.true.)
    call gvt_transform(ensembleTrl,'AllTransformedToModel',allowOverWrite_opt=.true.)
    call tmg_start(104,'LETKF-writeEns')
    call ens_writeEnsemble(ensembleAnl, '.', '', ' ', 'ENS_ANL', 'A',  &
                           numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                           containsFullField_opt=.true.)
    call tmg_stop(104)

    ! WARNING: Increment put in ensembleTrl for output
    call ens_add(ensembleAnl, ensembleTrl, scaleFactorInOut_opt=-1.0D0)
    call tmg_start(104,'LETKF-writeEns')
    call ens_writeEnsemble(ensembleTrl, '.', '', ' ', 'ENS_INC', 'R',  &
                           numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                           containsFullField_opt=.false.)
    ! Also write the reference (analysis) surface pressure to increment files
    call enkf_writeToAllMembers(stateVectorMeanAnlSfcPresMpiGlb, nEns,  &
                                etiket='ENS_INC', typvar='A', fileNameSuffix='inc',  &
                                ensPath='.')
    call tmg_stop(104)

    !- Output the sub-sampled ensemble analyses and increments
    if (writeSubSample) then

      ! Output the ensemble mean increment (include MeanAnl Psfc)
      call fln_ensAnlFileName( outFileName, 'subspace', tim_getDateStamp(), 0, ensFileNameSuffix_opt='inc' )
      do stepIndex = 1, tim_nstepobsinc
        call gsv_writeToFile(stateVectorMeanIncSubSample, outFileName, 'ENSMEAN_INC',  &
                             typvar_opt='R', writeHeightSfc_opt=.false., numBits_opt=16, &
                             stepIndex_opt=stepIndex, containsFullField_opt=.false.)
        call gsv_writeToFile(stateVectorMeanAnlSfcPres, outFileName, 'ENSMEAN_INC',  &
                             typvar_opt='A', writeHeightSfc_opt=.true., &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true.)
      end do

      ! Output the ensemble mean analysis state
      call fln_ensAnlFileName( outFileName, 'subspace', tim_getDateStamp(), 0 )
      do stepIndex = 1, tim_nstepobsinc
        call gsv_writeToFile(stateVectorMeanAnlSubSample, outFileName, 'ENSMEAN_ANL',  &
                             typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true.)
      end do

      ! Output the sub-sampled analysis ensemble members
      call tmg_start(104,'LETKF-writeEns')
      call ens_writeEnsemble(ensembleAnlSubSample, 'subspace', '', ' ', 'ENS_ANL', 'A',  &
                             numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                             containsFullField_opt=.true.)
      call tmg_stop(104)

      ! Output the sub-sampled ensemble increments (include MeanAnl Psfc)
      ! WARNING: Increment put in ensembleTrlSubSample for output
      call ens_add(ensembleAnlSubSample, ensembleTrlSubSample, scaleFactorInOut_opt=-1.0D0)
      call tmg_start(104,'LETKF-writeEns')
      call ens_writeEnsemble(ensembleTrlSubSample, 'subspace', '', ' ', 'ENS_INC', 'R',  &
                             numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                             containsFullField_opt=.false.)
      ! Also write the reference (analysis) surface pressure to increment files
      call enkf_writeToAllMembers(stateVectorMeanAnlSfcPresMpiGlb,  &
                                  ens_getNumMembers(ensembleAnlSubSample),  &
                                  etiket='ENS_INC', typvar='A', fileNameSuffix='inc',  &
                                  ensPath='subspace')
      call tmg_stop(104)
    end if

    call tmg_stop(4)

  end subroutine enkf_postProcess

  !----------------------------------------------------------------------
  ! enkf_getRMSEtiket (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_getRmsEtiket(etiketMean, etiketStd, etiketType, etiket0, nEns)
    !
    !:Purpose:   Return the appropriate string to use for the etiket
    !            in standard files containing the ensemble mean and
    !            spread.
    !
    implicit none

    ! arguments:
    character(len=*) :: etiketMean
    character(len=*) :: etiketStd
    character(len=*) :: etiketType
    character(len=*) :: etiket0
    integer          :: nEns

    if (trim(etiketType) == 'F') then

      ! create trialrms etiket, e.g. E2AVGTRPALL E24_3GMP0256
      etiketStd(1:5) = etiket0(1:5)
      etiketStd(6:7) = 'GM'
      etiketStd(8:8) = etiket0(8:8)
      write(etiketStd(9:12),'(I4.4)') nEns
      etiketMean(1:2) = etiket0(1:2)
      etiketMean(3:7) = 'AVGTR'
      etiketMean(8:8) = etiket0(8:8)
      etiketMean(9:11) = 'ALL'

    else if (trim(etiketType) == 'A') then

      ! create analrms etiket, e.g. E2AVGANPALL, E24_3_0P0256
      etiketStd(1:8) = etiket0(1:8)
      write(etiketStd(9:12),'(I4.4)') nEns
      etiketMean(1:2) = etiket0(1:2)
      etiketMean(3:7)='AVGAN'
      etiketMean(8:8) = etiket0(8:8)
      etiketMean(9:11) = 'ALL'

    else if (trim(etiketType) == 'P') then

      ! create analpertrms etiket, e.g. E2AVGPTPALL, E24_3PTP0256
      etiketStd(1:5) = etiket0(1:5)
      etiketStd(6:7) = 'PT'
      etiketStd(8:8) = etiket0(8:8)
      write(etiketStd(9:12),'(I4.4)') nEns
      etiketMean(1:2) = etiket0(1:2)
      etiketMean(3:7) = 'AVGPT'
      etiketMean(8:8) = etiket0(8:8)
      etiketMean(9:11) = 'ALL'

    else
      call utl_abort('midas-letkf: unknown value of etiketType')
    end if

  end subroutine enkf_getRmsEtiket

  !----------------------------------------------------------------------
  ! enkf_writeToAllMembers (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_writeToAllMembers(stateVector, nEns, etiket, typvar,  &
                                    fileNameSuffix, ensPath)
    !
    !:Purpose:   Write the contents of the supplied stateVector to all
    !            ensemble member files in an efficient parallel way.
    !
    implicit none

    ! arguments:
    type(struct_gsv) :: stateVector
    integer          :: nEns
    character(len=*) :: etiket
    character(len=*) :: typvar
    character(len=*) :: fileNameSuffix
    character(len=*) :: ensPath

    ! locals:
    integer            :: memberIndex, stepIndex, writeFilePE(nEns)
    character(len=4)   :: memberIndexStr
    character(len=256) :: outFileName

    do memberIndex = 1, nEns
      writeFilePE(memberIndex) = mod(memberIndex-1, mpi_nprocs)
    end do

    do memberIndex = 1, nEns

      if (mpi_myid == writeFilePE(memberIndex)) then

        call fln_ensAnlFileName( outFileName, ensPath, tim_getDateStamp(),  &
                                 memberIndex, ensFileNameSuffix_opt=fileNameSuffix )
        write(memberIndexStr,'(I4.4)') memberIndex

        do stepIndex = 1, tim_nstepobsinc
          call gsv_writeToFile(stateVector, outFileName,  &
                               trim(etiket) // memberIndexStr,  &
                               typvar_opt=trim(typvar), writeHeightSfc_opt=.true., &
                               stepIndex_opt=stepIndex, containsFullField_opt=.true., &
                               numBits_opt=16)
        end do

      end if

    end do

  end subroutine enkf_writeToAllMembers

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
  ! enkf_RTPS
  !--------------------------------------------------------------------------
  subroutine enkf_RTPS(ensembleAnl, ensembleTrl, stateVectorStdDevAnl,  &
                       stateVectorStdDevTrl, stateVectorMeanAnl, alphaRTPS)
    ! :Purpose: Apply Relaxation To Prior Spread ensemble inflation according
    !           to the factor alphaRTPS (usually between 0 and 1).
    implicit none

    ! Arguments
    type(struct_ens) :: ensembleAnl
    type(struct_ens) :: ensembleTrl
    type(struct_gsv) :: stateVectorStdDevAnl
    type(struct_gsv) :: stateVectorStdDevTrl
    type(struct_gsv) :: stateVectorMeanAnl
    real(8)          :: alphaRTPS

    ! Locals
    integer :: varLevIndex, latIndex, lonIndex, stepIndex, memberIndex
    integer :: nEns, numVarLev, myLonBeg, myLonEnd, myLatBeg, myLatEnd
    real(8) :: factorRTPS
    real(4), pointer     :: stdDevTrl_ptr_r4(:,:,:,:), stdDevAnl_ptr_r4(:,:,:,:)
    real(4), pointer     :: meanAnl_ptr_r4(:,:,:,:), memberAnl_ptr_r4(:,:,:,:)

    write(*,*) 'enkf_RTPS: Starting'

    stdDevTrl_ptr_r4 => gsv_getField_r4(stateVectorStdDevTrl)
    stdDevAnl_ptr_r4 => gsv_getField_r4(stateVectorStdDevAnl)
    meanAnl_ptr_r4 => gsv_getField_r4(stateVectorMeanAnl)

    nEns = ens_getNumMembers(ensembleAnl)
    numVarLev = ens_getNumK(ensembleAnl)
    call ens_getLatLonBounds(ensembleAnl, myLonBeg, myLonEnd, myLatBeg, myLatEnd)
    do varLevIndex = 1, numVarLev
      memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,varLevIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = 1, tim_nstepobsinc
            ! compute the inflation factor for RTPS
            if ( stdDevAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) > 0.0 ) then
              factorRTPS = 1.0D0 + alphaRTPS *  &
                           ( stdDevTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) -  &
                             stdDevAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) ) /  &
                           stdDevAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)
            else
              factorRTPS = 0.0D0
            end if
            ! apply the inflation factor to all Anl members (in place)
            if (factorRTPS > 0.0D0) then
              do memberIndex = 1, nEns
                memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) =  &
                     factorRTPS *  &
                     ( memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) -  &
                       meanAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) ) +  &
                     meanAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)
              end do ! memberIndex
            end if ! factorRTPS > 0
          end do ! stepIndex
        end do ! lonIndex
      end do ! latIndex
    end do ! varLevIndex

    write(*,*) 'enkf_RTPS: Finished'

  end subroutine enkf_RTPS

  !--------------------------------------------------------------------------
  ! enkf_addRandomPert
  !--------------------------------------------------------------------------
  subroutine enkf_addRandomPert(ensembleAnl, stateVectorMeanTrl, alphaRandomPert, randomSeed)
    ! :Purpose: Apply additive inflation using random perturbations from sampling
    !           the B matrix as defined by the regular namelist block NAMBHI, NAMBEN, etc.
    !           The scale factor alphaRandomPert (usually between 0 and 1) defines is used
    !           to simply multiply the resulting perturbations before adding to the original
    !           members. The perturbations have zero ensemble mean.
    implicit none

    ! Arguments
    type(struct_ens) :: ensembleAnl
    type(struct_gsv) :: stateVectorMeanTrl
    real(8)          :: alphaRandomPert
    integer          :: randomSeed

    ! Locals
    type(struct_gsv)         :: stateVectorPerturbation
    type(struct_gsv)         :: stateVectorPerturbationInterp
    type(struct_gsv)         :: statevectorTrlHU
    type(struct_gsv)         :: stateVectorVtr
    type(struct_vco), pointer :: vco_randomPert, vco_ens
    type(struct_hco), pointer :: hco_randomPert, hco_ens
    character(len=12)  :: etiket
    real(8), allocatable :: controlVector_mpiglobal(:), controlVector(:)
    real(8), allocatable :: perturbationMean(:,:,:)
    real(8), allocatable :: PsfcReference(:,:,:)
    real(8), pointer     :: perturbation_ptr(:,:,:)
    real(4), pointer     :: memberAnl_ptr_r4(:,:,:,:)
    integer :: cvIndex, memberIndex, varLevIndex, lonIndex, latIndex, stepIndex
    integer :: nEns, numVarLev, myLonBeg, myLonEnd, myLatBeg, myLatEnd
    logical, save :: firstCall = .true.

    ! Get ensemble dimensions
    nEns = ens_getNumMembers(ensembleAnl)
    numVarLev = ens_getNumK(ensembleAnl)
    call ens_getLatLonBounds(ensembleAnl, myLonBeg, myLonEnd, myLatBeg, myLatEnd)
    vco_ens => ens_getVco(ensembleAnl)
    hco_ens => ens_getHco(ensembleAnl)

    ! Define the horiz/vertical coordinate for perturbation calculation
    nullify(vco_randomPert)
    nullify(hco_randomPert)
    call hco_setupFromFile(hco_randomPert, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN
    if ( hco_randomPert % global ) then
      etiket = 'BGCK_STDDEV'
    else
      etiket = 'STDDEV'
    end if
    call vco_setupFromFile(vco_randomPert, './bgcov', etiket)
    if (firstCall) then
      call bmat_setup(hco_randomPert, vco_randomPert)
      firstCall = .false.
    end if
    call gvt_setup(hco_randomPert, vco_randomPert)

    call rng_setup(abs(randomSeed))

    allocate(controlVector(cvm_nvadim))
    allocate(controlVector_mpiglobal(cvm_nvadim_mpiglobal))
    allocate(perturbationMean(myLonBeg:myLonEnd,myLatBeg:myLatEnd,numVarLev))
    perturbationMean(:,:,:) = 0.0d0

    call gsv_allocate(stateVectorPerturbation, 1, hco_randomPert, vco_randomPert, &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.,  &
                      hInterpolateDegree_opt='LINEAR')
    call gsv_allocate(stateVectorPerturbationInterp, 1, hco_ens, vco_ens, &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.,  &
                      hInterpolateDegree_opt='LINEAR')
    perturbation_ptr => gsv_getField3d_r8(stateVectorPerturbationInterp)
    allocate(PsfcReference(myLonBeg:myLonEnd,myLatBeg:myLatEnd,1))
    PsfcReference(:,:,:) = 100000.0D0

    ! prepare the ensemble mean HU field for transforming LQ to HU perturbations
    call gsv_allocate(statevectorTrlHU, 1, hco_ens, vco_ens,   &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocHeightSfc_opt=.true., hInterpolateDegree_opt='LINEAR', &
                      hExtrapolateDegree_opt='MINIMUM', &
                      varNames_opt=(/'HU','P0'/) )
    call gsv_copy(stateVectorMeanTrl, stateVectorTrlHU, allowMismatch_opt=.true.)
    call gsv_allocate(stateVectorVtr, 1, hco_randomPert, vco_randomPert,   &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocHeightSfc_opt=.true., hInterpolateDegree_opt='LINEAR', &
                      hExtrapolateDegree_opt='MINIMUM', &
                      varNames_opt=(/'HU','P0'/) )
    call gsv_interpolate(stateVectorTrlHU, stateVectorVtr)

    do memberIndex = 1, nEns

      if( mpi_myid == 0 ) then
        write(*,*) 
        write(*,*) 'Computing random perturbation number= ', memberIndex
      end if

      ! global vector random control vector (independent of mpi topology)
      do cvIndex = 1, cvm_nvadim_mpiglobal
        controlVector_mpiglobal(cvIndex) = rng_gaussian()
      end do
      call bmat_reduceToMPILocal( controlVector, controlVector_mpiglobal )

      call bmat_sqrtB(controlVector, cvm_nvadim, &       ! IN
                      stateVectorPerturbation,   &       ! OUT
                      stateVectorRef_opt=stateVectorVtr) ! IN

      call gsv_interpolate(stateVectorPerturbation, stateVectorPerturbationInterp, &
                           PsfcReference_opt=PsfcReference)

      ! scale the perturbation by the specified factor
      call gsv_scale(stateVectorPerturbationInterp, alphaRandomPert)

      write(*,*) 'enkf_addRandomPert: member ', memberIndex, ', perturbation min/maxval = ',  &
                 minval(perturbation_ptr), maxval(perturbation_ptr)

      do varLevIndex = 1, numVarLev
        memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,varLevIndex)
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            perturbationMean(lonIndex, latIndex, varLevIndex) =   &
                 perturbationMean(lonIndex, latIndex, varLevIndex) +  &
                 perturbation_ptr(lonIndex, latIndex, varLevIndex) / real(nEns, 8)
            do stepIndex = 1, tim_nstepobsinc
              memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) =  &
                   memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) + perturbation_ptr(lonIndex, latIndex, varLevIndex)
            end do ! stepIndex
          end do ! lonIndex
        end do ! latIndex
      end do ! varLevIndex

    end do ! memberIndex

    ! remove the ensemble mean of the perturbations
    do varLevIndex = 1, numVarLev
      memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,varLevIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = 1, tim_nstepobsinc
            do memberIndex = 1, nEns
              memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) =  &
                   memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) - perturbationMean(lonIndex, latIndex, varLevIndex)
            end do ! memberIndex
          end do ! stepIndex
        end do ! lonIndex
      end do ! latIndex
    end do ! varLevIndex

    deallocate(controlVector)
    deallocate(controlVector_mpiglobal)
    deallocate(perturbationMean)
    call gsv_deallocate(stateVectorPerturbation)
    call gsv_deallocate(stateVectorPerturbationInterp)
    deallocate(PsfcReference)
    call gsv_deallocate(statevectorTrlHU)
    call gsv_deallocate(stateVectorVtr)

  end subroutine enkf_addRandomPert

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

  !-----------------------------------------------------------------
  ! enkf_selectSubSample
  !-----------------------------------------------------------------
  subroutine enkf_selectSubSample(ensembleAnl, ensembleTrl,  &
                                  ensembleAnlSubSample, ensembleTrlSubSample)
    ! :Purpose: Create sub-sampled ensembles of analyses and trials based on
    !           the contents of the ascii files 'sampletable' which lists
    !           the member indices for the subsample.
    implicit none

    ! Arguments
    type(struct_ens) :: ensembleAnl
    type(struct_ens) :: ensembleTrl
    type(struct_ens) :: ensembleAnlSubSample
    type(struct_ens) :: ensembleTrlSubSample

    ! Locals
    type(struct_gsv) :: stateVectorMember
    integer :: nulFile, ierr, status, numSubSample
    integer :: memberIndex, memberIndexSubSample, memberIndexFull
    integer :: memberIndexesSubSample(1000), memberIndexesFull(1000)
    integer, allocatable :: dateStampListInc(:)
    integer, external :: fnom, fclos

    numSubSample = 0

    nulFile = 0
    ierr = fnom(nulFile, './sampletable', 'FTN+SEQ+R/O', 0)
    do
      read(nulFile,*, IOSTAT=status) memberIndexSubSample, memberIndexFull
      if (status < 0) exit

      numSubSample = numSubSample + 1
      write(*,*) 'enkf_selectSubSample: ', memberIndexSubSample, memberIndexFull
      memberIndexesSubSample(numSubSample) = memberIndexSubSample
      memberIndexesFull(numSubSample) = memberIndexFull
    end do
    ierr = fclos(nulFile)

    write(*,*) 'enkf_selectSubSample: number of subSample members = ', numSubSample

    allocate(dateStampListInc(tim_nstepobsinc))
    call tim_getstamplist(dateStampListInc,tim_nstepobsinc,tim_getDatestamp())

    call ens_allocate(ensembleAnlSubSample, numSubSample, tim_nstepobsinc,  &
                      ens_getHco(ensembleAnl), ens_getVco(ensembleAnl), dateStampListInc)
    call ens_allocate(ensembleTrlSubSample, numSubSample, tim_nstepobsinc,  &
                      ens_getHco(ensembleAnl), ens_getVco(ensembleAnl), dateStampListInc)

    call gsv_allocate(stateVectorMember, tim_nstepobsinc,  &
                      ens_getHco(ensembleAnl), ens_getVco(ensembleAnl),  &
                      dateStamp_opt=tim_getDateStamp(),  &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                      dataKind_opt=4, allocHeightSfc_opt=.false., &
                      allocHeight_opt=.false., allocPressure_opt=.false. )

    do memberIndex = 1, numSubSample
      ! copy analysis ensemble member
      call ens_copyMember(ensembleAnl, stateVectorMember,  &
                          memberIndexesFull(memberIndex))
      call ens_insertMember(ensembleAnlSubSample, stateVectorMember,  &
                            memberIndexesSubSample(memberIndex))

      ! copy trial ensemble member
      call ens_copyMember(ensembleTrl, stateVectorMember,  &
                          memberIndexesFull(memberIndex))
      call ens_insertMember(ensembleTrlSubSample, stateVectorMember,  &
                            memberIndexesSubSample(memberIndex))
    end do

    call gsv_deallocate(stateVectorMember)
    deallocate(dateStampListInc)

  end subroutine enkf_selectSubSample

  !-----------------------------------------------------------------
  ! enkf_hybridRecentering
  !-----------------------------------------------------------------
  subroutine enkf_hybridRecentering(ensembleAnl, weightRecenter,  &
                                    useOptionTableRecenter, numMembersToRecenter)
    ! :Purpose: Modify an ensemble by recentering the members on a state provided
    !           in the file "recentering_analysis".
    !           The "weightRecenter" and "numMembersToRecenter" are used in the calculation
    !           to determine the amount of recentering and how many members it is
    !           applied to. Alternatively the information in the "optiontable" file
    !           can be used to perform a different amount of recentering on each
    !           member.
    implicit none

    ! Arguments
    type(struct_ens) :: ensembleAnl
    real(8)          :: weightRecenter
    logical          :: useOptionTableRecenter
    integer          :: numMembersToRecenter

    ! Locals
    type(struct_gsv) :: stateVectorRecenterAnl
    type(struct_hco), pointer :: hco_ens => null()
    type(struct_vco), pointer :: vco_ens => null()
    character(len=30)    :: recenterAnlFileName = 'recentering_analysis'
    character(len=20)    :: stringArray(100)
    character(len=1000)  :: textLine
    integer              :: stepIndex, memberIndex, columnIndex
    integer              :: numMembers, numColumns, nulFile, status
    logical              :: recenterAnlFileExists
    real(8), allocatable :: weightArray(:)
    integer, external    :: fnom, fclos

    ! check if recentering analysis file exists
    inquire(file=recenterAnlFileName, exist=recenterAnlFileExists)
    if (.not. recenterAnlFileExists) then
      write(*,*) 'enkf_hybridRecentering: RecenterAnlFileName = ', recenterAnlFileName
      call utl_abort('enkf_hybridRecentering: The recentering analysis file does not exist')
    end if

    numMembers = ens_getNumMembers(ensembleAnl)

    ! read the optiontable file, if requested
    if (useOptionTableRecenter) then
      write(*,*) 'enkf_hybridRecentering: using optiontable file to specify recentering weights.'
      nulFile = 0
      status = fnom(nulFile, './optiontable', 'FMT+SEQ+R/O', 0)
      read(nulFile,'(a)', IOSTAT=status) textLine
      if (status /= 0) then
        call utl_abort('enkf_hybridRecentering: unable to read optiontable file')
      end if
      call utl_parseColumns(textLine, numColumns)
      if (mpi_myid==0) write(*,*) 'enkf_hybridRecentering: optiontable file has ', numColumns, ' columns.'
      allocate( weightArray(0:numMembers) )
      rewind(nulFile)
      do memberIndex = 0, numMembers
        read(nulFile,'(a)') textLine
        call utl_parseColumns(textLine, numColumns, stringArray_opt=stringArray)
        if (mpi_myid==0) write(*,*) memberIndex, (stringArray(columnIndex),columnIndex=1,numColumns)
        read(stringArray(numColumns),'(f6.3)') weightArray(memberIndex)
        if (mpi_myid==0) write(*,*) 'weightArray = ', weightArray(memberIndex)
      end do
      status = fclos(nulFile)
    end if

    ! allocate and read in recentering analysis state
    hco_ens => ens_getHco(ensembleAnl)
    vco_ens => ens_getVco(ensembleAnl)
    call gsv_allocate( stateVectorRecenterAnl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.false., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorRecenterAnl)

    do stepIndex = 1, tim_nstepobsinc
      call gsv_readFromFile( stateVectorRecenterAnl, recenterAnlFileName, ' ', ' ',  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true., &
                             readHeightSfc_opt=.false. )
    end do

    ! apply recentering
    if (useOptionTableRecenter) then
      call ens_recenter(ensembleAnl, stateVectorRecenterAnl,  &
                        recenteringCoeffArray_opt=weightArray(1:numMembers),  &
                        numMembersToRecenter_opt=numMembersToRecenter)
    else
      call ens_recenter(ensembleAnl, stateVectorRecenterAnl,  &
                        recenteringCoeff_opt=weightRecenter,  &
                        numMembersToRecenter_opt=numMembersToRecenter)
    end if

    call gsv_deallocate(stateVectorRecenterAnl)
    if ( allocated(weightArray) ) deallocate(weightArray)

  end subroutine enkf_hybridRecentering

end module enkf_mod

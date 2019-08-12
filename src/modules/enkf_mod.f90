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
  use columnData_mod
  use timeCoord_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use ensembleObservations_mod
  use randomNumber_mod
  use controlVector_mod
  use gridVariableTransforms_mod
  use bMatrix_mod
  use localizationFunction_mod
  use varNameList_mod
  implicit none
  save
  private

  ! public types
  public :: struct_enkfInterpInfo

  ! public procedures
  public :: enkf_setupInterpInfo, enkf_interpWeights, enkf_addRandomPert, enkf_RTPS
  public :: enkf_selectSubSample, enkf_LETKFanalyses

  public :: enkf_computeColumnsMean, enkf_computeColumnsPerturbations
  public :: enkf_gatherHX

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
  subroutine enkf_LETKFanalyses(ensembleAnl, ensembleTrl, ensObs_mpiglobal,  &
                                       stateVectorMeanInc, stateVectorMeanTrl, stateVectorMeanAnl, &
                                       stateVectorDeterInc, stateVectorDeterTrl, stateVectorDeterAnl, &
                                       wInterpInfo, maxNumLocalObs,  &
                                       hLocalize, vLocalize, updateMembers, alphaRTPP, mpiDistribution)
    ! :Purpose: Local subroutine containing the code for computing
    !           the LETKF analyses for all ensemble members, ensemble
    !           mean and (if present) a deterministic state.
    implicit none

    ! Arguments
    type(struct_ens), pointer   :: ensembleTrl
    type(struct_ens)            :: ensembleAnl
    type(struct_eob)            :: ensObs_mpiglobal
    type(struct_gsv)            :: stateVectorMeanInc
    type(struct_gsv)            :: stateVectorMeanTrl
    type(struct_gsv)            :: stateVectorMeanAnl
    type(struct_gsv)            :: stateVectorDeterInc
    type(struct_gsv)            :: stateVectorDeterTrl
    type(struct_gsv)            :: stateVectorDeterAnl
    type(struct_enkfInterpInfo) :: wInterpInfo
    integer :: maxNumLocalObs
    real(8) :: hLocalize        ! horizontal localization radius (in km)
    real(8) :: vLocalize        ! vertical localization radius (in units of ln(Pressure in Pa))
    logical :: updateMembers
    real(8) :: alphaRTPP
    character(len=*) :: mpiDistribution

    ! Locals
    integer :: nEns, nLev_M, ierr
    integer :: memberIndex, memberIndex1, memberIndex2, procIndex, procIndexSend
    integer :: latIndex, lonIndex, stepIndex, kIndex, levIndex, levIndex2
    integer :: bodyIndex, localObsIndex, numLocalObs, numLocalObsFound
    integer :: countMaxExceeded, maxCountMaxExceeded, numGridPointWeights
    integer :: myNumLatLonRecv, myNumLatLonSend, numLatLonRecvMax
    integer :: numLatLonTotalUnique, latLonIndex
    integer :: sendTag, recvTag, nsize, numRecv, numSend
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd, numK
    integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    real(8) :: anlLat, anlLon, anlLogPres, distance
    real(8) :: localization

    integer, allocatable :: localBodyIndices(:)
    integer, allocatable :: myLatIndexesRecv(:), myLonIndexesRecv(:)
    integer, allocatable :: myLatIndexesSend(:), myLonIndexesSend(:)
    integer, allocatable :: myNumProcIndexesSend(:)
    integer, allocatable :: myProcIndexesRecv(:), myProcIndexesSend(:,:)
    integer, allocatable :: requestIdRecv(:), requestIdSend(:)

    real(8), allocatable :: distances(:)
    real(8), allocatable :: PaInv(:,:), PaSqrt(:,:), Pa(:,:), YbTinvR(:,:)
    real(8), allocatable :: weightsTemp(:)
    real(8), allocatable :: weightsMembers(:,:,:,:), weightsMembersLatLon(:,:,:)
    real(8), allocatable :: weightsMean(:,:,:,:), weightsMeanLatLon(:,:,:)
    real(8), allocatable :: weightsDeter(:,:,:,:), weightsDeterLatLon(:,:,:)
    real(8), allocatable :: memberAnlPert(:)
    real(4), allocatable :: logPres_M_r4(:,:,:)

    real(4), pointer     :: meanTrl_ptr_r4(:,:,:,:), meanAnl_ptr_r4(:,:,:,:), meanInc_ptr_r4(:,:,:,:)
    real(4), pointer     :: deterTrl_ptr_r4(:,:,:,:), deterAnl_ptr_r4(:,:,:,:), deterInc_ptr_r4(:,:,:,:)
    real(4), pointer     :: memberTrl_ptr_r4(:,:,:,:), memberAnl_ptr_r4(:,:,:,:)

    type(struct_hco), pointer :: hco_ens => null()

    logical :: deterExists
    logical :: firstTime = .true.

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
    myLonBeg = stateVectorMeanInc%myLonBeg
    myLonEnd = stateVectorMeanInc%myLonEnd
    myLatBeg = stateVectorMeanInc%myLatBeg
    myLatEnd = stateVectorMeanInc%myLatEnd
    numK     = stateVectorMeanInc%nk
    myLonBegHalo = wInterpInfo%myLonBegHalo
    myLonEndHalo = wInterpInfo%myLonEndHalo
    myLatBegHalo = wInterpInfo%myLatBegHalo
    myLatEndHalo = wInterpInfo%myLatEndHalo
    deterExists = stateVectorDeterTrl%allocated

    !
    ! Compute gridded 3D ensemble weights
    !
    allocate(localBodyIndices(maxNumLocalObs))
    allocate(distances(maxNumLocalObs))
    allocate(YbTinvR(nEns,maxNumLocalObs))
    allocate(PaInv(nEns,nEns))
    allocate(PaSqrt(nEns,nEns))
    allocate(Pa(nEns,nEns))
    allocate(memberAnlPert(nEns))
    allocate(weightsTemp(nEns))
    weightsTemp(:) = 0.0d0
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
    ! Weights for deterministic analysis
    if (deterExists) then
      allocate(weightsDeter(nEns,1,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
      weightsDeter(:,:,:,:) = 0.0d0
      allocate(weightsDeterLatLon(nEns,1,myNumLatLonSend))
      weightsDeterLatLon(:,:,:) = 0.0d0
    end if
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
        recvTag = (latIndex-1)*stateVectorMeanInc%ni + lonIndex

        nsize = nEns
        numRecv = numRecv + 1
        call mpi_irecv( weightsMean(:,1,lonIndex,latIndex),  &
                        nsize, mpi_datyp_real8, procIndex-1, recvTag,  &
                        mpi_comm_grid, requestIdRecv(numRecv), ierr )
        if (deterExists) then
          numRecv = numRecv + 1
          recvTag = recvTag + stateVectorMeanInc%ni*stateVectorMeanInc%nj
          call mpi_irecv( weightsDeter(:,1,lonIndex,latIndex),  &
                          nsize, mpi_datyp_real8, procIndex-1, recvTag,  &
                          mpi_comm_grid, requestIdRecv(numRecv), ierr )
        end if
        if (updateMembers) then
          nsize = nEns*nEns
          numRecv = numRecv + 1
          recvTag = recvTag + stateVectorMeanInc%ni*stateVectorMeanInc%nj
          call mpi_irecv( weightsMembers(:,:,lonIndex,latIndex),  &
                          nsize, mpi_datyp_real8, procIndex-1, recvTag,  &
                          mpi_comm_grid, requestIdRecv(numRecv), ierr )
        end if
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

        ! Get list of nearby observations and distances to gridpoint
        call tmg_start(9,'LETKF-getLocalBodyIndices')
        numLocalObs = eob_getLocalBodyIndices(ensObs_mpiglobal, localBodyIndices,     &
                                              distances, anlLat, anlLon, anlLogPres,  &
                                              hLocalize, vLocalize, numLocalObsFound)
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
          localization = lfn_Response(distances(localObsIndex),hLocalize)
          ! Vertical - use pressures at the grid point (not obs) location
          if (vLocalize > 0) then
            distance = abs( anlLogPres - ensObs_mpiglobal%logPres(bodyIndex) )
            localization = localization * lfn_Response(distance,vLocalize)
          end if
          call tmg_stop(18)

          do memberIndex = 1, nEns
            YbTinvR(memberIndex,localObsIndex) = ensObs_mpiglobal%Yb(memberIndex, bodyIndex) * &
                                                 localization * ensObs_mpiglobal%varObsInv(bodyIndex)
          end do

          if (localObsIndex == 1) PaInv(:,:) = 0.0D0
          do memberIndex2 = 1, nEns
            do memberIndex1 = 1, nEns
              PaInv(memberIndex1,memberIndex2) = PaInv(memberIndex1,memberIndex2) +  &
                   YbTinvR(memberIndex1,localObsIndex) * ensObs_mpiglobal%Yb(memberIndex2, bodyIndex)
            end do
          end do

        end do ! localObsIndex

        ! Rest of the computation of local weights for this grid point
        if (numLocalObs > 0) then

          ! Add second term of PaInv
          do memberIndex = 1, nEns
            PaInv(memberIndex,memberIndex) = PaInv(memberIndex,memberIndex) + real(nEns - 1,8)
          end do

          ! Compute Pa and sqrt(Pa) matrices from PaInv
          Pa(:,:) = PaInv(:,:)
          call tmg_start(90,'LETKF-matInverse')
          if (updateMembers) then
            call utl_matInverse(Pa, nEns, inverseSqrt_opt=PaSqrt)
          else
            call utl_matInverse(Pa, nEns)
          end if
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
              weightsMeanLatLon(memberIndex1,1,latLonIndex) = weightsMeanLatLon(memberIndex1,1,latLonIndex) +  &
                   Pa(memberIndex1,memberIndex2)*weightsTemp(memberIndex2)
            end do
          end do

          if (deterExists) then
            ! Compute deterministic analysis local weights as Pa * YbTinvR * (obs - deterYb)
            weightsTemp(:) = 0.0d0
            do localObsIndex = 1, numLocalObs
              bodyIndex = localBodyIndices(localObsIndex)
              do memberIndex = 1, nEns
                weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                           YbTinvR(memberIndex,localObsIndex) *  &
                                           ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                             ensObs_mpiglobal%deterYb(bodyIndex) )
              end do
            end do

            weightsDeterLatLon(:,1,latLonIndex) = 0.0d0
            do memberIndex2 = 1, nEns
              do memberIndex1 = 1, nEns
                weightsDeterLatLon(memberIndex1,1,latLonIndex) = weightsDeterLatLon(memberIndex1,1,latLonIndex) +  &
                     Pa(memberIndex1,memberIndex2)*weightsTemp(memberIndex2)
              end do
            end do
          end if

          ! Compute ensemble perturbation weights: (1-alphaRTPP)*[(Nens-1)^1/2*PaSqrt]+alphaRTPP*I
          if (updateMembers) then
            weightsMembersLatLon(:,:,latLonIndex) = (1.0d0 - alphaRTPP) * sqrt(real(nEns - 1,8)) * PaSqrt(:,:)
            do memberIndex = 1, nEns
              weightsMembersLatLon(memberIndex,memberIndex,latLonIndex) = alphaRTPP +  &
                   weightsMembersLatLon(memberIndex,memberIndex,latLonIndex)
            end do
          end if

        else
          ! no observations near this grid point, set weights to zero
          weightsMeanLatLon(:,1,latLonIndex) = 0.0d0
          if (deterExists) weightsDeterLatLon(:,1,latLonIndex) = 0.0d0
        end if ! numLocalObs > 0

        call tmg_stop(91)

        !
        ! Now post all send instructions (each lat-lon may be sent to multiple tasks)
        !
        call tmg_start(103,'LETKF-commWeights')
        latIndex = myLatIndexesSend(latLonIndex)
        lonIndex = myLonIndexesSend(latLonIndex)
        do procIndex = 1, myNumProcIndexesSend(latLonIndex)
          sendTag = (latIndex-1)*stateVectorMeanInc%ni + lonIndex
          procIndexSend = myProcIndexesSend(latLonIndex, procIndex)

          nsize = nEns
          numSend = numSend + 1
          call mpi_isend( weightsMeanLatLon(:,1,latLonIndex),  &
                          nsize, mpi_datyp_real8, procIndexSend-1, sendTag,  &
                          mpi_comm_grid, requestIdSend(numSend), ierr )
          if (deterExists) then
            numSend = numSend + 1
            sendTag = sendTag + stateVectorMeanInc%ni*stateVectorMeanInc%nj
            call mpi_isend( weightsDeterLatLon(:,1,latLonIndex),  &
                            nsize, mpi_datyp_real8, procIndexSend-1, sendTag,  &
                            mpi_comm_grid, requestIdSend(numSend), ierr )
          end if
          if (updateMembers) then
            nsize = nEns*nEns
            numSend = numSend + 1
            sendTag = sendTag + stateVectorMeanInc%ni*stateVectorMeanInc%nj
            call mpi_isend( weightsMembersLatLon(:,:,latLonIndex),  &
                            nsize, mpi_datyp_real8, procIndexSend-1, sendTag,  &
                            mpi_comm_grid, requestIdSend(numSend), ierr )
          end if
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

        if (deterExists) call enkf_interpWeights(wInterpInfo, weightsDeter)

        if (updateMembers) call enkf_interpWeights(wInterpInfo, weightsMembers)

      end if
      call tmg_stop(92)

      call tmg_start(100,'LETKF-applyWeights')

      !
      ! Apply the weights to compute the ensemble mean and members
      !
      meanInc_ptr_r4 => gsv_getField_r4(stateVectorMeanInc)
      meanTrl_ptr_r4 => gsv_getField_r4(stateVectorMeanTrl)
      meanAnl_ptr_r4 => gsv_getField_r4(stateVectorMeanAnl)
      if (deterExists) then
        deterInc_ptr_r4 => gsv_getField_r4(stateVectorDeterInc)
        deterTrl_ptr_r4 => gsv_getField_r4(stateVectorDeterTrl)
        deterAnl_ptr_r4 => gsv_getField_r4(stateVectorDeterAnl)
      end if
      do latIndex = myLatBeg, myLatEnd
        LON_LOOP5: do lonIndex = myLonBeg, myLonEnd

          ! skip this grid point if all weights zero (no nearby obs)
          if (all(weightsMean(:,1,lonIndex,latIndex) == 0.0d0)) cycle LON_LOOP5

          ! Compute the ensemble mean increment and analysis
          do kIndex = 1, numK
            ! Only treat kIndex values that correspond with current levIndex
            if (vnl_varLevelFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,kIndex)) == 'SF') then
              levIndex2 = nLev_M
            else
              levIndex2 = gsv_getLevFromK(stateVectorMeanInc,kIndex)
            end if
            if (levIndex2 /= levIndex) cycle
            memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,kIndex)
            do stepIndex = 1, tim_nstepobsinc
              ! mean increment
              do memberIndex = 1, nEns
                meanInc_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) = meanInc_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) +  &
                     weightsMean(memberIndex,1,lonIndex,latIndex) *  &
                     (memberTrl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) - meanTrl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex))
              end do ! memberIndex
              ! mean analysis
              meanAnl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) = meanTrl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) +  &
                   meanInc_ptr_r4(lonIndex,latIndex,kIndex,stepIndex)
            end do ! stepIndex
          end do ! kIndex

          ! Compute the deterministic increment and analysis
          if (deterExists) then
            do kIndex = 1, numK
              ! Only treat kIndex values that correspond with current levIndex
              if (vnl_varLevelFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,kIndex)) == 'SF') then
                levIndex2 = nLev_M
              else
                levIndex2 = gsv_getLevFromK(stateVectorMeanInc,kIndex)
              end if
              if (levIndex2 /= levIndex) cycle
              memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,kIndex)
              do stepIndex = 1, tim_nstepobsinc
                ! deterministic increment
                do memberIndex = 1, nEns
                  deterInc_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) = deterInc_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) +  &
                       weightsDeter(memberIndex,1,lonIndex,latIndex) *  &
                       (memberTrl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) - deterTrl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex))
                end do ! memberIndex
                ! deterministic analysis
                deterAnl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) = deterTrl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) +  &
                     deterInc_ptr_r4(lonIndex,latIndex,kIndex,stepIndex)
              end do ! stepIndex
            end do ! kIndex
          end if

          ! Compute the ensemble member analyses
          if (updateMembers) then
            do kIndex = 1, numK
              ! Only treat kIndex values that correspond with current levIndex
              if (vnl_varLevelFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,kIndex)) == 'SF') then
                levIndex2 = nLev_M
              else
                levIndex2 = gsv_getLevFromK(stateVectorMeanInc,kIndex)
              end if
              if (levIndex2 /= levIndex) cycle
              memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,kIndex)
              memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,kIndex)
              do stepIndex = 1, tim_nstepobsinc
                ! Compute analysis member perturbation
                memberAnlPert(:) = 0.0d0
                do memberIndex2 = 1, nEns
                  do memberIndex1 = 1, nEns
                    memberAnlPert(memberIndex2) = memberAnlPert(memberIndex2) + &
                         weightsMembers(memberIndex1,memberIndex2,lonIndex,latIndex) *  &
                         (memberTrl_ptr_r4(memberIndex1,stepIndex,lonIndex,latIndex) -  &
                          meanTrl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex))
                  end do ! memberIndex1
                end do ! memberIndex2
                ! Add analysis member perturbation to mean analysis
                memberAnl_ptr_r4(:,stepIndex,lonIndex,latIndex) =  &
                     meanAnl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) + memberAnlPert(:)
              end do ! stepIndex
            end do ! kIndex
          end if ! updateMembers

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

  end subroutine enkf_LETKFanalyses

  !----------------------------------------------------------------------
  ! enkf_computeLogPresM
  !----------------------------------------------------------------------
  subroutine enkf_computeLogPresM(logPres_M_r4,stateVectorMeanTrl)
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
    call vtr_transform(stateVectorMeanTrlPressure,'PsfcToP_nl')
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
  ! enkf_setupMpiDistribution
  !----------------------------------------------------------------------
  subroutine enkf_LETKFsetupMpiDistribution(myNumLatLonRecv, myNumLatLonSend, &
                                            myLatIndexesRecv, myLonIndexesRecv, &
                                            myLatIndexesSend, myLonIndexesSend, &
                                            myProcIndexesRecv, myProcIndexesSend, &
                                            myNumProcIndexesSend, mpiDistribution, wInterpInfo)

    ! :Purpose: Setup for distribution of grid points over mpi tasks.

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

  end subroutine enkf_LETKFsetupMpiDistribution

  !----------------------------------------------------------------------
  ! enkf_latLonAlreadyFound
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
                                  myLonBegHalo,myLonEndHalo,myLatBegHalo,myLatEndHalo)
    ! :Purpose: Setup the weights and lat/lon indices needed to bilinearly
    !           interpolate the LETKF weights from a coarse grid to the full
    !           resolution grid. The coarseness of the grid is specified by
    !           the weightLatLonStep argument.
    implicit none

    ! Arguments
    type(struct_enkfInterpInfo) :: wInterpInfo
    type(struct_hco) :: hco
    integer :: weightLatLonStep
    integer :: myLonBegHalo
    integer :: myLonEndHalo
    integer :: myLatBegHalo
    integer :: myLatEndHalo

    ! Locals
    integer :: lonIndex, latIndex, ni, nj
    real(8) :: interpWeightLon, interpWeightLat

    if (hco%grtyp == 'U') then
      call utl_abort('enkf_setupInterpInfo: Yin-Yang grid (U) not yet supported.')
    end if

    ni = hco%ni
    nj = hco%nj

    wInterpInfo%latLonStep   = weightLatLonStep
    wInterpInfo%myLonBegHalo = myLonBegHalo
    wInterpInfo%myLonEndHalo = myLonEndHalo
    wInterpInfo%myLatBegHalo = myLatBegHalo
    wInterpInfo%myLatEndHalo = myLatEndHalo

    allocate(wInterpInfo%numIndexes(myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
    if (weightLatLonStep > 1) then
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
      end if
      if (myLatEndHalo == nj) then
        do lonIndex = myLonBegHalo, myLonEndHalo, weightLatLonStep
          wInterpInfo%numIndexes(lonIndex,myLatEndHalo) = 0
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
          ! Find nearest grid point with a value towards left 
          wInterpInfo%numIndexes(lonIndex,latIndex) = 2
          wInterpInfo%lonIndexes(1,lonIndex,latIndex) = lonIndex
          wInterpInfo%lonIndexes(2,lonIndex,latIndex) = lonIndex
          wInterpInfo%latIndexes(1,lonIndex,latIndex) = myLatBegHalo +  &
               weightLatLonStep * floor(real(latIndex - myLatBegHalo)/real(weightLatLonStep)) 
          wInterpInfo%latIndexes(2,lonIndex,latIndex) = min(nj,  &
               wInterpInfo%latIndexes(1,lonIndex,latIndex) + weightLatLonStep)
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

  end subroutine enkf_setupInterpInfo

  !--------------------------------------------------------------------------
  ! enkf_interpWeights
  !--------------------------------------------------------------------------
  subroutine enkf_interpWeights(wInterpInfo, weights)
    ! :Purpose: Perform the bilinear interpolation of the weights
    !           using the precalculated interpolation info.
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
    integer :: kIndex, latIndex, lonIndex, stepIndex, memberIndex
    integer :: nEns, numK, myLonBeg, myLonEnd, myLatBeg, myLatEnd
    real(8) :: factorRTPS
    real(4), pointer     :: stdDevTrl_ptr_r4(:,:,:,:), stdDevAnl_ptr_r4(:,:,:,:)
    real(4), pointer     :: meanAnl_ptr_r4(:,:,:,:), memberAnl_ptr_r4(:,:,:,:)

    write(*,*) 'enkf_RTPS: Starting'

    stdDevTrl_ptr_r4 => gsv_getField_r4(stateVectorStdDevTrl)
    stdDevAnl_ptr_r4 => gsv_getField_r4(stateVectorStdDevAnl)
    meanAnl_ptr_r4 => gsv_getField_r4(stateVectorMeanAnl)

    nEns = ens_getNumMembers(ensembleAnl)
    numK = ens_getNumK(ensembleAnl)
    call ens_getLatLonBounds(ensembleAnl, myLonBeg, myLonEnd, myLatBeg, myLatEnd)
    do kIndex = 1, numK
      memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,kIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = 1, tim_nstepobsinc
            ! compute the inflation factor for RTPS
            if ( stdDevAnl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) > 0.0 ) then
              factorRTPS = 1.0D0 + alphaRTPS *  &
                           ( stdDevTrl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) -  &
                             stdDevAnl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) ) /  &
                           stdDevAnl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex)
            else
              factorRTPS = 0.0D0
            end if
            ! apply the inflation factor to all Anl members (in place)
            if (factorRTPS > 0.0D0) then
              do memberIndex = 1, nEns
                memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) =  &
                     factorRTPS *  &
                     ( memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) -  &
                       meanAnl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) ) +  &
                     meanAnl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex)
              end do ! memberIndex
            end if ! factorRTPS > 0
          end do ! stepIndex
        end do ! lonIndex
      end do ! latIndex
    end do ! kIndex

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
    integer :: cvIndex, memberIndex, kIndex, lonIndex, latIndex, stepIndex
    integer :: nEns, numK, myLonBeg, myLonEnd, myLatBeg, myLatEnd
    logical, save :: firstCall = .true.

    ! Get ensemble dimensions
    nEns = ens_getNumMembers(ensembleAnl)
    numK = ens_getNumK(ensembleAnl)
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
    allocate(perturbationMean(myLonBeg:myLonEnd,myLatBeg:myLatEnd,numK))
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

      do kIndex = 1, numK
        memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,kIndex)
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            perturbationMean(lonIndex, latIndex, kIndex) =   &
                 perturbationMean(lonIndex, latIndex, kIndex) +  &
                 perturbation_ptr(lonIndex, latIndex, kIndex) / real(nEns, 8)
            do stepIndex = 1, tim_nstepobsinc
              memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) =  &
                   memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) + perturbation_ptr(lonIndex, latIndex, kIndex)
            end do ! stepIndex
          end do ! lonIndex
        end do ! latIndex
      end do ! kIndex

    end do ! memberIndex

    ! remove the ensemble mean of the perturbations
    do kIndex = 1, numK
      memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,kIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = 1, tim_nstepobsinc
            do memberIndex = 1, nEns
              memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) =  &
                   memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) - perturbationMean(lonIndex, latIndex, kIndex)
            end do ! memberIndex
          end do ! stepIndex
        end do ! lonIndex
      end do ! latIndex
    end do ! kIndex

    deallocate(controlVector)
    deallocate(controlVector_mpiglobal)
    deallocate(perturbationMean)
    call gsv_deallocate(stateVectorPerturbation)
    call gsv_deallocate(stateVectorPerturbationInterp)
    deallocate(PsfcReference)
    call gsv_deallocate(statevectorTrlHU)
    call gsv_deallocate(stateVectorVtr)

  end subroutine enkf_addRandomPert

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
    integer :: fnom, fclos
    integer, allocatable :: dateStampListInc(:)

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
! The following routines are used only by the old ensembleH program
! and will soon be eliminated
!-----------------------------------------------------------------
  subroutine enkf_computeColumnsMean(column_mean, columns)
    implicit none

    ! Arguments:
    type(struct_columnData) :: column_mean, columns(:)

    ! Locals:
    logical :: verbose = .true.
    integer :: memberIndex, nEns, levIndex
    real(8) :: multFactor
    real(8), pointer :: column_ptr(:)

    call tmg_start(146,'ENKF_COLSMEAN')

    nEns = size(columns)
    multFactor = 1.0d0 / real(nEns,8)
    write(*,*) 'enkf_computeColumnsMean: nEns =', nEns

    call col_zero(column_mean)

    do memberIndex = 1, nEns

        column_mean%all(:,:) = column_mean%all(:,:) +  &
                               multFactor * columns(memberIndex)%all(:,:)

    end do

    if ( verbose ) then
      write(*,*) '======================='
      write(*,*) 'Contents of column_mean:'
      column_ptr => col_getColumn(column_mean,1,'UU')
      write(*,*) 'column UU = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'VV')
      write(*,*) 'column VV = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'TT')
      write(*,*) 'column TT = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'HU')
      write(*,*) 'column LQ = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'P0')
      write(*,*) 'column P0 = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'TG')
      write(*,*) 'column TG = ', column_ptr(:)
      write(*,*) '======================='

      do levIndex = 1, col_getNumLev(column_mean,'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getPressure(column_mean,levIndex,1,MM) = ',  &
                   levIndex,col_getPressure(column_mean,levIndex,1,'MM')
      end do
      do levIndex = 1, col_getNumLev(column_mean,'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getHeight(column_mean,levIndex,1,MM) = ',  &
                   levIndex,col_getHeight(column_mean,levIndex,1,'MM')
      end do
    end if

    call tmg_stop(146)

  end subroutine enkf_computeColumnsMean


  subroutine enkf_computeColumnsPerturbations(columns, column_mean)
    implicit none

    ! Arguments:
    type(struct_columnData) :: columns(:), column_mean

    ! Locals:
    logical :: verbose = .true.
    integer :: memberIndex, nEns, levIndex
    real(8), pointer :: column_ptr(:)

    call tmg_start(147,'ENKF_COLSPERTS')

    nEns = size(columns)
    write(*,*) 'enkf_computeColumnsPerturbations: nEns =', nEns

    !
    ! Remove ensemble mean from all variables, except: HeightSfc, oltv
    !
    do memberIndex = 1, nEns

        columns(memberIndex)%all(:,:) = columns(memberIndex)%all(:,:) -  &
                                        column_mean%all(:,:)
    end do

    if ( verbose ) then
      write(*,*) '======================='
      write(*,*) 'Contents of columns(1):'
      column_ptr => col_getColumn(columns(1),1,'UU')
      write(*,*) 'column UU = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'VV')
      write(*,*) 'column VV = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'TT')
      write(*,*) 'column TT = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'HU')
      write(*,*) 'column LQ = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'P0')
      write(*,*) 'column P0 = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'TG')
      write(*,*) 'column TG = ', column_ptr(:)
      write(*,*) '======================='

      do levIndex = 1, col_getNumLev(columns(1),'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getPressure(columns(1),levIndex,1,MM) = ',  &
                   levIndex,col_getPressure(columns(1),levIndex,1,'MM')
      end do
      do levIndex = 1, col_getNumLev(columns(1),'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getHeight(columns(1),levIndex,1,MM) = ',  &
                   levIndex,col_getHeight(columns(1),levIndex,1,'MM')
      end do
    end if

    call tmg_stop(147)

  end subroutine enkf_computeColumnsPerturbations


  subroutine enkf_gatherHX(HXens,HXensT_mpiglobal)
    implicit none

    ! Arguments:
    real(8) :: HXens(:,:)
    real(8), allocatable :: HXensT_mpiglobal(:,:)

    ! Locals:
    integer :: ierr, nEns, numBody, procIndex, memberIndex, numBody_mpiglobal
    integer :: allNumBody(mpi_nprocs), displs(mpi_nprocs)

    numBody = size(HXens,1)
    nEns     = size(HXens,2)

    write(*,*) 'enkf_gatherHX: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call rpn_comm_gather( numBody, 1, 'mpi_integer', allNumBody, 1, 'mpi_integer', &
                          0, 'GRID', ierr )
    if ( mpi_myid == 0 ) then
      displs(1) = 0
      do procIndex = 2, mpi_nprocs
        displs(procIndex) = displs(procIndex-1) + allNumBody(procIndex-1)
      end do
    else
      displs(:) = 0
    end if

    numBody_mpiglobal = sum(allNumBody(:))
    if( mpi_myid == 0 ) then
      allocate(HXensT_mpiglobal(nEns,numBody_mpiglobal))
    else
      allocate(HXensT_mpiglobal(nEns,1))
    end if

    do memberIndex = 1, nEns
      call rpn_comm_gatherv( HXens(:,memberIndex), numBody, 'mpi_double_precision', &
                             HXensT_mpiglobal(memberIndex,:), allNumBody, displs, &
                             'mpi_double_precision', 0, 'GRID', ierr )
    end do

    write(*,*) 'enkf_gatherHX: finished'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine enkf_gatherHX

end module enkf_mod

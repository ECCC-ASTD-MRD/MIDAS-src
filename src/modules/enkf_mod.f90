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
  !           an EnKF in MIDAS.
  !
  use mpi_mod
  use gridStateVector_mod
  use mathPhysConstants_mod
  use utilities_mod
  use fileNames_mod
  use varNameList_mod
  use tt2phi_mod
  use obsSpaceData_mod
  use columnData_mod

  use timeCoord_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use randomNumber_mod
  use controlVector_mod
  use variableTransforms_mod
  use bMatrix_mod
  implicit none
  save
  private

  ! public types
  public :: struct_enkfInterpInfo

  ! public procedures
  public :: enkf_setupInterpInfo, enkf_interpWeights, enkf_addRandomPert, enkf_RTPS

  public :: enkf_computeColumnsMean, enkf_computeColumnsPerturbations
  public :: enkf_gatherHX

  ! for weight interpolation
  type struct_enkfInterpInfo
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

  !--------------------------------------------------------------------------
  ! enkf_setupInterpInfo
  !--------------------------------------------------------------------------
  subroutine enkf_setupInterpInfo(wInterpInfo, weightLatLonStep, ni, nj,  &
                                  myLonBegHalo,myLonEndHalo,myLatBegHalo,myLatEndHalo)
    ! :Purpose: Setup the weights and lat/lon indices needed to bilinearly
    !           interpolate the LETKF weights from a coarse grid to the full
    !           resolution grid. The coarseness of the grid is specified by
    !           the weightLatLonStep argument.
    implicit none

    ! Arguments
    type(struct_enkfInterpInfo) :: wInterpInfo
    integer :: weightLatLonStep
    integer :: ni
    integer :: nj
    integer :: myLonBegHalo
    integer :: myLonEndHalo
    integer :: myLatBegHalo
    integer :: myLatEndHalo

    ! Locals
    integer :: lonIndex, latIndex
    real(8) :: interpWeightLon, interpWeightLat

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
        do latIndex = myLatBegHalo, myLatEndHalo
          if (wInterpInfo%numIndexes(myLonBegHalo,latIndex) == 0) cycle
          wInterpInfo%numIndexes(myLonEndHalo,latIndex) = 0
        end do
      end if
      if (myLatEndHalo == nj) then
        do lonIndex = myLonBegHalo, myLonEndHalo
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

    ! Get ensemble dimensions
    nEns = ens_getNumMembers(ensembleAnl)
    numK = ens_getNumK(ensembleAnl)
    call ens_getLatLonBounds(ensembleAnl, myLonBeg, myLonEnd, myLatBeg, myLatEnd)
    vco_ens => ens_getVco(ensembleAnl)
    hco_ens => ens_getHco(ensembleAnl)

    ! Define the horiz/vertical coordinate for perturbation calculation
    nullify(vco_randomPert)
    nullify(hco_randomPert)
    call hco_setupFromFile(hco_randomPert, './analysisgrid_forRandPert', 'ANALYSIS', 'Analysis' ) ! IN
    if ( hco_randomPert % global ) then
      etiket = 'BGCK_STDDEV'
    else
      etiket = 'STDDEV'
    end if
    call vco_setupFromFile(vco_randomPert, './bgcov', etiket)
    call bmat_setup(hco_randomPert, vco_randomPert)
    call vtr_setup(hco_randomPert, vco_randomPert)

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

      if( mpi_myid == 0 ) write(*,*) 'Computing random perturbation number= ', memberIndex

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

      write(*,*) 'enkf_addRandomPert: perturbation min/maxval = ',  &
                 minval(perturbation_ptr), maxval(perturbation_ptr)

      do kIndex = 1, numK
        memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,kIndex)
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            !ensemble_r4(lonIndex, latIndex, kIndex, memberIndex) = real(perturbation_ptr(lonIndex, latIndex, kIndex), 4)
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

  end subroutine enkf_addRandomPert


!-----------------------------------------------------------------
! The following routines are used only by the ensembleH program
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

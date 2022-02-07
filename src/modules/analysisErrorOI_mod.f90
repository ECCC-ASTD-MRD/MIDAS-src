!-------------------------------------- LICENCE BEGIN ------------------------------------
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

module analysisErrorOI_mod
  ! MODULE analysisErrorOI (prefix='aer' category='3. High-level transformations')
  !
  ! :Purpose: Calculate the analysis-error standard deviation.
  !           The method used is Optimal Interpolation,
  !           where it is assumed that only a subset of the
  !           total number of observations influence the analysis at a given grid point.
  !           By default, everything in the module is private.
  !           The data is accessed by external subroutines through public subroutines
  !           and functions calls.
  !

  use columnData_mod
  use gridStateVector_mod
  use kdtree2_mod
  use MathPhysConstants_mod
  use mpi_mod
  use obserrors_mod
  use obsSpaceData_mod
  use oceanMask_mod
  use physicsFunctions_mod
  use stateToColumn_mod
  use varNamelist_mod
  use utilities_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use obsOperators_mod

  implicit none

  private

  ! public subroutines and functions
  public :: aer_analysisError, aer_daysSinceLastObs

  type struct_neighborhood
   integer          :: numObs
   integer, pointer :: headerIndex(:)
   integer, pointer :: bodyIndex(:)
  end type struct_neighborhood

  type(struct_neighborhood), pointer :: influentObs(:,:)

  real(8), allocatable :: Lcorr(:,:)

  type(kdtree2), pointer :: tree

  integer :: ni, nj

  integer, external :: get_max_rss

  integer, parameter :: maxNumLocalGridptsSearch = 1000

  real(8) :: interpWeight(maxNumLocalGridptsSearch)
  integer :: obsLatIndex(maxNumLocalGridptsSearch), obsLonIndex(maxNumLocalGridptsSearch)

contains

  !--------------------------------------------------------------------------
  ! aer_analysisError
  !--------------------------------------------------------------------------
  subroutine aer_analysisError(obsSpaceData, hco_ptr, vco_ptr, trlmFileName)
    !
    ! :Purpose: Calculate analysis-error variance.
    !
    implicit none

    ! Arguments
    type(struct_obs), intent(in) :: obsSpaceData
    type(struct_hco), pointer    :: hco_ptr
    type(struct_vco), pointer    :: vco_ptr
    character(len=*), intent(in) :: trlmFileName

    ! Local Variables
    integer :: fnom, fclos, nulnam, ierr

    type(struct_gsv) :: stateVectorBkGnd
    type(struct_gsv) :: stateVectorAnal
    integer, allocatable :: numObs(:,:)
    real(8), pointer :: bkGndErrorStdDev_ptr(:,:,:,:), analysisErrorStdDev_ptr(:,:,:,:)

    real(8), allocatable :: obsOperator(:,:), Bmatrix(:,:), PHiA(:), &
                            innovCovariance(:,:), obsErrorVariance(:), &
                            PH(:,:), KH(:), IKH(:), innovCovarianceInverse(:,:)
    integer, parameter :: maxvar = 15000

    integer :: statei(maxvar)    ! Model grid coordinate i of neighbors
    integer :: statej(maxvar)    ! Model grid coordinate j of neighbors

    integer :: numVariables, varIndex1, varIndex2, currentAnalVarIndex

    logical :: found

    integer :: stepIndex, levIndex, lonIndex, latIndex, gridIndex, headerIndex, &
               bodyIndex, kIndex, procIndex
    integer :: influentObsIndex2, influentObsIndex

    real(8) :: scaling
    integer :: xStateIndex, yStateIndex

    integer :: numInfluentObs, xIndex1, yIndex1, xIndex2, yIndex2
    integer :: gridptCount, gridpt

    real(8) :: distance

    real(kdkind), allocatable :: positionArray(:,:)
    real(8),      allocatable :: latInRad(:,:), lonInRad(:,:)

    character(len=*), parameter :: myName = 'aer_analysisError'
    character(len=*), parameter :: correlationLengthFileName = './bgstddev'
    type(struct_gsv)            :: statevector
    real(4), pointer            :: field3D_r4_ptr(:,:,:)

    character(len=2 ) :: typvar
    character(len=12) :: etiket

    type(struct_columnData) :: column
    type(struct_columnData) :: columng

    real(8) :: leadTimeInHours

    ! namelist variables:
    real(8) :: maxAnalysisErrorStdDev
    namelist /namaer/ maxAnalysisErrorStdDev

    if( mpi_nprocs > 1 ) then
      write(*,*) 'mpi_nprocs = ',mpi_nprocs
      call utl_abort( myName// &
                      ': this version of the code should only be used with one mpi task.')
    end if
    if( mpi_myid > 0 ) return

    write(*,*) '**********************************************************'
    write(*,*) '** '//myName//': Calculate analysis-error variance **'
    write(*,*) '**********************************************************'

    ! reading namelist variables
    ! default values
    maxAnalysisErrorStdDev = 1.0

    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam, nml = namaer, iostat = ierr)
    if ( ierr /= 0 ) call utl_abort( myName//': Error reading namelist')
    ierr = fclos(nulnam)

    call gsv_allocate( stateVectorBkGnd, 1, hco_ptr, vco_ptr, dateStamp_opt=-1, &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       varNames_opt=(/'GLE'/), dataKind_opt=8 )
    call gsv_allocate( stateVectorAnal,  1, hco_ptr, vco_ptr, dateStamp_opt=-1, &
                       varNames_opt=(/'GLE'/), dataKind_opt=8 )

    etiket = '            '
    typvar = 'P@'

    call gsv_readFromFile( stateVectorBkGnd, trlmFileName, etiket, typvar )

    leadTimeInHours = real(stateVectorBkGnd%deet*stateVectorBkGnd%npasList(1),8)/3600.0d0
    call incdatr(stateVectorAnal%dateOriginList(1), stateVectorBkGnd%dateOriginList(1), &
                 leadTimeInHours)

    call ocm_copyMask(stateVectorBkGnd%oceanMask, stateVectorAnal%oceanMask)
    stateVectorAnal%etiket = stateVectorBkGnd%etiket

    call col_setVco(column, vco_ptr)
    call col_allocate(column,  obs_numHeader(obsSpaceData), varNames_opt=(/'GLE'/))
    call col_setVco(columng, vco_ptr)
    call col_allocate(columng, obs_numHeader(obsSpaceData), varNames_opt=(/'GLE'/))
    call s2c_tl(statevectorBkGnd, column, columng, obsSpaceData)

    ni = stateVectorBkGnd%hco%ni
    nj = stateVectorBkGnd%hco%nj

    allocate( lonInRad( ni, nj ), latInRad( ni, nj) )
    allocate( Lcorr( ni, nj ) )

    write(*,*) myName// &
               ': Correlation length scale 2D field will be read from the file: ', &
               correlationLengthFileName
    call gsv_allocate( statevector, 1, hco_ptr, vco_ptr, dateStamp_opt=-1, &
                       dataKind_opt=4, hInterpolateDegree_opt='LINEAR', &
                       varNames_opt=(/'GL'/) )
    call gsv_zero( statevector )
    call gsv_readFromFile( statevector, correlationLengthFileName, 'CORRLEN', ' ', &
                           unitConversion_opt = .false. )

    call gsv_getField( statevector, field3D_r4_ptr, 'GL' )
    ! Convert from km to meters
    Lcorr(:,:) = 1000.0d0*real(field3D_r4_ptr( :, :, 1 ), 8)
    write(*,*) myName//': correlation length scale 2D field for variable GL min/max: ', &
               minval( Lcorr(:,:) ), maxval( Lcorr(:,:) )
    call gsv_deallocate( statevector )

    ! create kdtree
    write(*,*) myName//': start creating kdtree for stateVectorBkGnd'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(positionArray(3,ni*nj))

    gridIndex = 0
    do latIndex = 1, nj
      do lonIndex = 1, ni

        gridIndex = gridIndex + 1
        latInRad(lonIndex,latIndex) = real( &
                                      stateVectorBkGnd%hco%lat2d_4(lonIndex,latIndex), 8)
        lonInRad(lonIndex,latIndex) = real( &
                                      stateVectorBkGnd%hco%lon2d_4(lonIndex,latIndex), 8)

        positionArray(:,gridIndex) = kdtree2_3dPosition(lonInRad(lonIndex,latIndex), &
                                                        latInRad(lonIndex,latIndex))

      end do
    end do

    write(*,*) myName//': latInRad min/max: ', minval(latInRad(:,:) ), &
                                               maxval( latInRad(:,:) )
    write(*,*) myName//': lonInRad min/max: ', minval(lonInRad(:,:) ), &
                                               maxval( lonInRad(:,:) )

    nullify(tree)
    tree => kdtree2_create(positionArray, sort=.false., rearrange=.true.) 

    deallocate(positionArray)

    write(*,*) myName//': done creating kdtree for stateVectorBkGnd'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! Go through all observations a first time to get
    ! the number of influential observations
    ! in order to allocate memory appropriately.

    allocate(influentObs(ni, nj))
    allocate(     numObs(ni, nj))

    call findObs(obsSpaceData, stateVectorBkGnd, numObs)

    ! Memory allocation

    do latIndex = 1, nj
      do lonIndex = 1, ni
        influentObs(lonIndex,latIndex)%numObs = numObs(lonIndex,latIndex)
        allocate(influentObs(lonIndex,latIndex)%headerIndex(numObs(lonIndex,latIndex)))
        allocate(influentObs(lonIndex,latIndex)%bodyIndex(numObs(lonIndex,latIndex)))
      end do
    end do

    ! Go through all observations a second time to get
    ! the indexes of the influential observations.
    ! This is not efficient to go through all observations twice
    ! but it saves lot of memory space.

    call findObs(obsSpaceData, stateVectorBkGnd, numObs)

    deallocate(numObs)

    write(*,*) myName//': obs found'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! Calculate analysis-error one analysis variable (grid point) at a time

    call gsv_getField(stateVectorBkGnd, bkGndErrorStdDev_ptr, 'GLE')
    call gsv_getField(stateVectorAnal,  analysisErrorStdDev_ptr, 'GLE')

    ! Initialisation
    do stepIndex = 1, stateVectorBkGnd%numStep
      do levIndex = 1, gsv_getNumLev(stateVectorBkGnd,vnl_varLevelFromVarname('GLE'))
        do latIndex = 1, nj
          do lonIndex = 1, ni
            analysisErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
               bkGndErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex)
          end do
        end do
      end do
    end do

    ! Only variables assigned within or by the loop can be private.

    STEP: do stepIndex = 1, stateVectorBkGnd%numStep
      LEVEL: do levIndex = 1, &
                           gsv_getNumLev(stateVectorBkGnd,vnl_varLevelFromVarname('GLE'))

        !$omp parallel do default(shared) schedule(dynamic) private(obsOperator, Bmatrix, &
        !$omp        PHiA, innovCovariance, innovCovarianceInverse, obsErrorVariance, &
        !$omp        PH, KH, IKH, &
        !$omp        statei, statej, &
        !$omp        numVariables, varIndex1, varIndex2, currentAnalVarIndex, found, &
        !$omp        influentObsIndex2, influentObsIndex, &
        !$omp        scaling, xStateIndex, yStateIndex, lonIndex, latIndex, headerIndex, &
        !$omp        bodyIndex, kIndex, procIndex, &
        !$omp        numInfluentObs, xIndex1, yIndex1, xIndex2, yIndex2, distance, &
        !$omp        interpWeight, obsLatIndex, obsLonIndex, gridptCount, gridpt)
        YINDEX: do latIndex = 1, nj
          XINDEX: do lonIndex = 1, ni

            numInfluentObs = influentObs(lonIndex,latIndex)%numObs

            if ( numInfluentObs == 0 .or. &
                 .not. stateVectorBkGnd%oceanMask%mask(lonIndex,latIndex,levIndex) ) cycle

            ! form the observation-error covariance (diagonal) matrix

            allocate(obsErrorVariance(numInfluentObs))

            do influentObsIndex = 1, numInfluentObs
              bodyIndex = influentObs(lonIndex,latIndex)%bodyIndex(influentObsIndex)
              obsErrorVariance(influentObsIndex) = ( obs_bodyElem_r(obsSpaceData,OBS_OER, &
                                                                    bodyIndex) )**2
            end do

            ! find all model variables involved here

            numVariables = 0
            statei(:) = 0
            statej(:) = 0
            do influentObsIndex = 1, numInfluentObs

              headerIndex = influentObs(lonIndex,latIndex)%headerIndex(influentObsIndex)
              bodyIndex   = influentObs(lonIndex,latIndex)%bodyIndex(influentObsIndex)

              scaling = oop_iceScaling(obsSpaceData, bodyIndex)

              if ( scaling /= 0.0d0 ) then

                do kIndex = stateVectorBkGnd%mykBeg, stateVectorBkGnd%mykEnd
                  do procIndex = 1, mpi_nprocs

                    call s2c_getWeightsAndGridPointIndexes(headerIndex, kIndex, &
                         stepIndex, procIndex, interpWeight, obsLatIndex, obsLonIndex, &
                         gridptCount)

                    do gridpt = 1, gridptCount

                      if ( interpWeight(gridpt) /= 0.0d0 ) then

                        xStateIndex = obsLonIndex(gridpt)
                        yStateIndex = obsLatIndex(gridpt)

                        found = .false.
                        do varIndex2=1,numVariables
                          if(xStateIndex == statei(varIndex2) .and. &
                             yStateIndex == statej(varIndex2)) then
                            found = .true.
                            exit
                          end if
                        end do
                        if(.not. found) then
                          numVariables = numVariables + 1
                          if(numVariables > maxvar) then
                            call utl_abort( myName//': Structure state'// &
                                 ' too small in subroutine'// &
                                 ' analysis_error_mod'// &
                                 ' Increase maxvar')
                          end if
                          statei(numVariables) = xStateIndex
                          statej(numVariables) = yStateIndex
                        end if

                      end if

                    end do

                  end do
                end do

              end if

            end do

            ! make sure that current analysis variable is part of the state
            ! vector even if it does not participate in obs calculation

            found = .false.
            do varIndex2=1,numVariables
              if(lonIndex == statei(varIndex2) .and. &
                 latIndex == statej(varIndex2)) then
                currentAnalVarIndex = varIndex2
                found = .true.
                exit
              end if
            end do
            if(.not. found) then
              numVariables = numVariables + 1
              if(numVariables > maxvar) then
                call utl_abort(myName//': Structure state too small in'// &
                     ' subroutine analysis_error_mod'// &
                     ' increase maxvar')
              end if
              currentAnalVarIndex = numVariables
              statei(numVariables) = lonIndex
              statej(numVariables) = latIndex
            end if

            ! form the observation operator matrix (obsOperator)

            allocate(obsOperator(numVariables,numInfluentObs))

            obsOperator(:,:) = 0.0d0

            varIndex2 = 0

            do influentObsIndex = 1, numInfluentObs

              headerIndex = influentObs(lonIndex,latIndex)%headerIndex(influentObsIndex)
              bodyIndex   = influentObs(lonIndex,latIndex)%bodyIndex(influentObsIndex)

              scaling = oop_iceScaling(obsSpaceData, bodyIndex)

              if ( scaling /= 0.0d0 ) then

                do kIndex = stateVectorBkGnd%mykBeg, stateVectorBkGnd%mykEnd
                  do procIndex = 1, mpi_nprocs

                    call s2c_getWeightsAndGridPointIndexes(headerIndex, kIndex, &
                         stepIndex, procIndex, interpWeight, obsLatIndex, obsLonIndex, &
                         gridptCount)

                    do gridpt = 1, gridptCount

                      if ( interpWeight(gridpt) /= 0.0d0 ) then

                        xStateIndex = obsLonIndex(gridpt)
                        yStateIndex= obsLatIndex(gridpt)

                        found = .false.
                        do while (.not. found)
                          if(varIndex2 < numVariables) then
                            varIndex2 = varIndex2 + 1
                          else
                            varIndex2 = 1
                          end if
                          if(xStateIndex == statei(varIndex2) .and. &
                             yStateIndex == statej(varIndex2)) then
                            found = .true.
                            obsOperator(varIndex2,influentObsIndex) = scaling* &
                                                                      interpWeight(gridpt)
                          end if
                        end do
                        if(.not. found) then
                          write(*,*) 'xStateIndex = ',xStateIndex
                          write(*,*) 'yStateIndex = ',yStateIndex
                          write(*,*) 'lonIndex = ',lonIndex
                          write(*,*) 'latIndex = ',latIndex
                          write(*,*) 'gridptCount = ',gridptCount
                          write(*,*) 'gridpt = ',gridpt
                          write(*,*) 'numVariables = ',numVariables
                          do varIndex2=1,numVariables
                            write(*,*) 'varIndex2 statei statej = ',statei(varIndex2), &
                                                                    statej(varIndex2)
                          end do
                          call utl_abort( myName//': not found in state vect.' )
                        end if

                      end if

                    end do

                  end do
                end do

              end if

            end do

            ! form the background-error covariance matrix

            allocate(Bmatrix(numVariables,numVariables))

            Bmatrix(:,:) = 0.0d0

            do varIndex2 = 1, numVariables

              xIndex2 = statei(varIndex2)
              yIndex2 = statej(varIndex2)

              do varIndex1 = varIndex2, numVariables

                xIndex1 = statei(varIndex1)
                yIndex1 = statej(varIndex1)

                if(xIndex2 == xIndex1 .and. yIndex2 == yIndex1) then
                  Bmatrix(varIndex1,varIndex2) = bkGndErrorStdDev_ptr(xIndex1,yIndex1,levIndex,stepIndex)**2
                else if(Lcorr(xIndex2,yIndex2) > 0.0d0) then
                  distance = phf_calcDistance(latInRad(xIndex2,yIndex2), &
                                              lonInRad(xIndex2,yIndex2), &
                                              latInRad(xIndex1,yIndex1), &
                                              lonInRad(xIndex1,yIndex1))
                  Bmatrix(varIndex1,varIndex2) = bkGndErrorStdDev_ptr(xIndex2,yIndex2,levIndex,stepIndex)* &
                                     bkGndErrorStdDev_ptr(xIndex1,yIndex1,levIndex,stepIndex) &
                                     *exp(-0.5*(distance/Lcorr(xIndex2,yIndex2))**2)
                end if

                ! symmetric matrix !
                Bmatrix(varIndex2,varIndex1) = Bmatrix(varIndex1,varIndex2)

              end do

            end do

            ! form the observation background error covariance matrix (PHT)

            allocate(PH(numInfluentObs,numVariables))

            ! PH = matmul (Bmatrix, transpose(obsOperator))
            do varIndex1 = 1, numVariables
              do influentObsIndex = 1, numInfluentObs
                PH(influentObsIndex,varIndex1) = dot_product(Bmatrix(:,varIndex1), &
                                                 obsOperator(:,influentObsIndex))
              end do
            end do

            !  covariance matrix of the innovation (HPHT + R)

            allocate(innovCovariance(numInfluentObs,numInfluentObs), &
                     innovCovarianceInverse(numInfluentObs,numInfluentObs))

            do influentObsIndex = 1, numInfluentObs

              do influentObsIndex2 = 1, numInfluentObs
                if ( influentObsIndex == influentObsIndex2 ) then
                  innovCovariance(influentObsIndex2,influentObsIndex) = &
                       dot_product(obsOperator(:,influentObsIndex), &
                       PH(influentObsIndex2,:)) + &
                       obsErrorVariance(influentObsIndex)
                else
                  innovCovariance(influentObsIndex2,influentObsIndex) = &
                       dot_product(obsOperator(:,influentObsIndex), &
                       PH(influentObsIndex2,:))
                end if
              end do

            end do

            ! Inverse of the covariance matrix of the innovation

            innovCovarianceInverse(:,:) = innovCovariance(:,:)
            call utl_matInverse(innovCovarianceInverse, numInfluentObs, printInformation_opt=.false.)

            ! Kalman gain; this is the row corresponding to the analysis variable

            allocate(PHiA(numInfluentObs))

            do influentObsIndex = 1, numInfluentObs
              PHiA(influentObsIndex) = dot_product(PH(:,currentAnalVarIndex), &
                                       innovCovarianceInverse(influentObsIndex,:))
            end do

            ! compute the error variance of the analysis

            allocate(KH(numVariables))

            ! KH = matmul (PHiA, obsOperator)
            do varIndex1 = 1, numVariables
              KH(varIndex1) = dot_product(PHiA, obsOperator(varIndex1,:))
            end do

            ! IKH = I - KH

            allocate(IKH(numVariables))

            IKH = -KH

            IKH(currentAnalVarIndex) = 1.0d0 - KH(currentAnalVarIndex)

            analysisErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
                 sqrt( dot_product (IKH, Bmatrix(:,currentAnalVarIndex)) )

            if(analysisErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex) < 0.0) then
              write(*,*) myName//'negative analysis-error Std dev. = ', &
                   analysisErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex), &
                   ' reset to zero at grid point (',lonIndex,latIndex,')'
              analysisErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
                   max(analysisErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex), 0.0)
            end if

            if(analysisErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex) > &
                  bkGndErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex)) then
              write(*,*) myName//'analysis-error Std dev. = ', &
                   analysisErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex), &
                   ' is larger than background-error ', &
                   'Std dev. is kept at = ', &
                   bkGndErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex)
              analysisErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
                   min(bkGndErrorStdDev_ptr(lonIndex,latIndex,levIndex,stepIndex), &
                   maxAnalysisErrorStdDev)
            end if

            deallocate(Bmatrix, obsErrorVariance, innovCovariance, &
                       innovCovarianceInverse, obsOperator, PH, PHiA, KH, IKH)
            deallocate(influentObs(lonIndex,latIndex)%headerIndex)
            deallocate(influentObs(lonIndex,latIndex)%bodyIndex)

          end do XINDEX
        end do YINDEX
        !$omp end parallel do

      end do LEVEL
    end do STEP

    deallocate(influentObs)
    deallocate( Lcorr )
    deallocate( lonInRad, latInRad )

    call gsv_writeToFile( stateVectorAnal, './anlm_000m', '', typvar_opt='A@', &
                          containsFullField_opt=.true. )

    call gsv_deallocate( stateVectorBkGnd )
    call gsv_deallocate( stateVectorAnal )

  end subroutine aer_analysisError

  !---------------------------------------------------------
  ! findObs
  !---------------------------------------------------------
  subroutine findObs(obsSpaceData, stateVectorBkGnd, numObs)
    !
    !:Purpose: Find all observations used for the local analysis.
    !
    implicit none

    ! Arguments
    type(struct_obs), intent(in)  :: obsSpaceData
    type(struct_gsv), intent(in)  :: stateVectorBkGnd
    integer,          intent(out) :: numObs(ni,nj)

    ! Local variables
    integer :: headerIndex, bodyIndexBeg, bodyIndexEnd, bodyIndex, kIndex, stepIndex, &
               procIndex
    integer :: gridptCount, gridpt, numLocalGridptsFoundSearch

    integer :: lonIndex, latIndex, resultsIndex, gridIndex
    type(kdtree2_result) :: searchResults(maxNumLocalGridptsSearch)
    real(kdkind)         :: refPosition(3), maxRadiusSquared
    real(4) :: footprintRadius_r4 ! (metres)
    real(4) :: influenceRadius_r4 ! (metres)

    real(8) :: obsLonInRad, obsLatInRad, maxLcorr

    call tmg_start(189,'AER_findObs')

    maxLcorr = maxval( Lcorr(:,:) )

    numObs(:,:) = 0

    HEADER_LOOP: do headerIndex = 1, obs_numHeader(obsSpaceData)

      bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1

      BODY_LOOP: do bodyIndex = bodyIndexBeg, bodyIndexEnd

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) &
             cycle BODY_LOOP

        footprintRadius_r4 = s2c_getFootprintRadius(obsSpaceData, stateVectorBkGnd, &
                                                    headerIndex)
        influenceRadius_r4 = max(0.0, footprintRadius_r4) + maxLcorr

        if ( maxLcorr == 0.0d0 ) then

          do kIndex = stateVectorBkGnd%mykBeg, stateVectorBkGnd%mykEnd
            do stepIndex = 1, stateVectorBkGnd%numStep
              do procIndex = 1, mpi_nprocs

                call s2c_getWeightsAndGridPointIndexes(headerIndex, &
                     kIndex, stepIndex, procIndex, interpWeight, obsLatIndex, &
                     obsLonIndex, gridptCount)

                GRIDPT_LOOP: do gridpt = 1, gridptCount

                  if ( interpWeight(gridpt) == 0.0d0 ) cycle GRIDPT_LOOP

                  lonIndex = obsLonIndex(gridpt)
                  latIndex = obsLatIndex(gridpt)

                  numObs(lonIndex,latIndex) = numObs(lonIndex,latIndex) + 1
                  if(associated(influentObs(lonIndex,latIndex)%bodyIndex)) then
                    if(numObs(lonIndex,latIndex) > influentObs(lonIndex,latIndex)%numObs) then
                      call utl_abort('findObs: Array too small')
                    end if
                    influentObs(lonIndex,latIndex)%headerIndex(numObs(lonIndex,latIndex)) = headerIndex
                    influentObs(lonIndex,latIndex)%bodyIndex(numObs(lonIndex,latIndex)) = bodyIndex
                  end if

                end do GRIDPT_LOOP

              end do
            end do
          end do

        else

          ! Determine the grid point nearest the observation.

          obsLonInRad = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex)
          obsLatInRad = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex)

          ! do the search
          maxRadiusSquared = real(influenceRadius_r4, 8) ** 2
          refPosition(:) = kdtree2_3dPosition(obsLonInRad, obsLatInRad)
          call kdtree2_r_nearest(tp=tree, qv=refPosition, r2=maxRadiusSquared, &
               nfound=numLocalGridptsFoundSearch, &
               nalloc=maxNumLocalGridptsSearch, &
               results=searchResults)

          if (numLocalGridptsFoundSearch > maxNumLocalGridptsSearch ) then
            call utl_abort('findObs: parameter maxNumLocalGridptsSearch must be increased')
          end if

          do resultsIndex = 1, numLocalGridptsFoundSearch

            gridIndex = searchResults(resultsIndex)%idx
            if ( gridIndex < 1 .or. &
                 gridIndex > stateVectorBkGnd%hco%ni * stateVectorBkGnd%hco%nj ) then
              write(*,*) 'analysisErrorStdDev: gridIndex=', gridIndex
              call utl_abort('findObs: gridIndex out of bound.')
            end if

            latIndex = (gridIndex - 1) / stateVectorBkGnd%hco%ni + 1
            lonIndex = gridIndex - (latIndex - 1) * stateVectorBkGnd%hco%ni
            if ( lonIndex < 1 .or. lonIndex > stateVectorBkGnd%hco%ni .or. &
                 latIndex < 1 .or. latIndex > stateVectorBkGnd%hco%nj ) then
              write(*,*) 'analysisErrorStdDev: lonIndex=', lonIndex, &
                                            ', latIndex=', latIndex
              call utl_abort('findObs: lonIndex/latIndex out of bound.')
            end if

            numObs(lonIndex,latIndex) = numObs(lonIndex,latIndex) + 1
            if(associated(influentObs(lonIndex,latIndex)%bodyIndex)) then
              if(numObs(lonIndex,latIndex) > influentObs(lonIndex,latIndex)%numObs) then
                call utl_abort('findObs: Array too small in subroutine findObs')
              end if
              influentObs(lonIndex,latIndex)%headerIndex(numObs(lonIndex,latIndex)) = &
                          headerIndex
              influentObs(lonIndex,latIndex)%bodyIndex(numObs(lonIndex,latIndex)) = &
                          bodyIndex
            end if

          end do

        end if

      end do BODY_LOOP

    end do HEADER_LOOP

    call tmg_stop(189)

  end subroutine findObs

  !---------------------------------------------------------
  ! aer_daysSinceLastObs
  !---------------------------------------------------------
  subroutine aer_daysSinceLastObs(obsSpaceData, hco_ptr, vco_ptr, trlmFileName)
    !
    !:Purpose: Update the field "days since last obs" with the newly assimilated obs.
    !
    implicit none

    ! Arguments
    type(struct_obs), intent(in)  :: obsSpaceData
    type(struct_hco), pointer    :: hco_ptr
    type(struct_vco), pointer    :: vco_ptr
    character(len=*), intent(in) :: trlmFileName

    ! Local Variables
    integer :: fnom, fclos, nulnam, ierr

    type(struct_gsv) :: stateVectorBkGnd
    type(struct_gsv) :: stateVectorAnal
    real(8), pointer :: bkGndDaysSinceLastObs_ptr(:,:,:,:), analysisDaysSinceLastObs_ptr(:,:,:,:)

    integer :: stepIndex, levIndex, lonIndex, latIndex, headerIndex
    integer :: bodyIndexBeg, bodyIndexEnd, bodyIndex, kIndex, procIndex

    integer :: gridptCount, gridpt

    character(len=*), parameter :: myName = 'aer_daysSinceLastObs'

    character(len=2 ) :: typvar
    character(len=12) :: etiket

    type(struct_columnData) :: column
    type(struct_columnData) :: columng

    real(8) :: leadTimeInHours

    if( mpi_nprocs > 1 ) then
      write(*,*) 'mpi_nprocs = ',mpi_nprocs
      call utl_abort( myName// &
                      ': this version of the code should only be used with one mpi task.')
    end if
    if( mpi_myid > 0 ) return

    write(*,*) '**********************************************************'
    write(*,*) '** '//myName//': Update the days since last obs **'
    write(*,*) '**********************************************************'

    call gsv_allocate( stateVectorBkGnd, 1, hco_ptr, vco_ptr, dateStamp_opt=-1, &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       varNames_opt=(/'DSLO'/), dataKind_opt=8 )
    call gsv_allocate( stateVectorAnal,  1, hco_ptr, vco_ptr, dateStamp_opt=-1, &
                       varNames_opt=(/'DSLO'/), dataKind_opt=8 )

    etiket = '            '
    typvar = 'P@'

    call gsv_readFromFile( stateVectorBkGnd, trlmFileName, etiket, typvar )

    leadTimeInHours = real(stateVectorBkGnd%deet*stateVectorBkGnd%npasList(1),8)/3600.0d0
    call incdatr(stateVectorAnal%dateOriginList(1), stateVectorBkGnd%dateOriginList(1), &
                 leadTimeInHours)

    call ocm_copyMask(stateVectorBkGnd%oceanMask, stateVectorAnal%oceanMask)
    stateVectorAnal%etiket = stateVectorBkGnd%etiket

    call col_setVco(column, vco_ptr)
    call col_allocate(column,  obs_numHeader(obsSpaceData), varNames_opt=(/'DSLO'/))
    call col_setVco(columng, vco_ptr)
    call col_allocate(columng, obs_numHeader(obsSpaceData), varNames_opt=(/'DSLO'/))
    call s2c_tl(statevectorBkGnd, column, columng, obsSpaceData)

    ni = stateVectorBkGnd%hco%ni
    nj = stateVectorBkGnd%hco%nj

    call gsv_getField(stateVectorBkGnd, bkGndDaysSinceLastObs_ptr, 'DSLO')
    call gsv_getField(stateVectorAnal,  analysisDaysSinceLastObs_ptr, 'DSLO')

    ! Initialisation
    do stepIndex = 1, stateVectorBkGnd%numStep
      do levIndex = 1, gsv_getNumLev(stateVectorBkGnd,vnl_varLevelFromVarname('DSLO'))
        do latIndex = 1, nj
          do lonIndex = 1, ni
            analysisDaysSinceLastObs_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
               bkGndDaysSinceLastObs_ptr(lonIndex,latIndex,levIndex,stepIndex)
          end do
        end do
      end do
    end do

    HEADER_LOOP: do headerIndex = 1, obs_numHeader(obsSpaceData)

      bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1

      BODY_LOOP: do bodyIndex = bodyIndexBeg, bodyIndexEnd

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) then
          cycle BODY_LOOP
        end if

        do kIndex = stateVectorBkGnd%mykBeg, stateVectorBkGnd%mykEnd
          do stepIndex = 1, stateVectorBkGnd%numStep
            do procIndex = 1, mpi_nprocs

              call s2c_getWeightsAndGridPointIndexes(headerIndex, &
                   kIndex, stepIndex, procIndex, interpWeight, obsLatIndex, &
                   obsLonIndex, gridptCount)

              GRIDPT_LOOP: do gridpt = 1, gridptCount

                if ( interpWeight(gridpt) == 0.0d0 ) cycle GRIDPT_LOOP

                lonIndex = obsLonIndex(gridpt)
                latIndex = obsLatIndex(gridpt)

                do levIndex = 1, gsv_getNumLev(stateVectorBkGnd,vnl_varLevelFromVarname('DSLO'))
                  analysisDaysSinceLastObs_ptr(lonIndex,latIndex,levIndex,stepIndex) = 0.0
                end do

              end do GRIDPT_LOOP

            end do
          end do
        end do

      end do BODY_LOOP

    end do HEADER_LOOP

    call gsv_writeToFile( stateVectorAnal, './anlm_000m', '', typvar_opt='A@', &
                          containsFullField_opt=.true. )

    call gsv_deallocate( stateVectorBkGnd )
    call gsv_deallocate( stateVectorAnal )

  end subroutine aer_daysSinceLastObs

end module analysisErrorOI_mod

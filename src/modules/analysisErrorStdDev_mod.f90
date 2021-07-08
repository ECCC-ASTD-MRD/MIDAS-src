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

module analysisError_mod
  ! MODULE analysisError (prefix='aer' category='3. High-level transformations')
  !
  ! :Purpose: Calculate the analysis-error standard deviation.
  !           The method used is Optimal Interpolation,
  !           where it is assumed that only a subset of the
  !           total number of observations influence the analysis at a given grid point.
  !           By default, everything in the module is private.
  !           The data is accessed by external subroutines through public subroutines
  !           and functions calls.
  !

  use bufr_mod
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

  implicit none

  private

  ! public subroutines and functions
  public :: aer_analysisError

  type struct_neighborhood
   integer          :: numObs
   integer, POINTER :: headerIndex(:)
   integer, POINTER :: bodyIndex(:)
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
  subroutine aer_analysisError(obsSpaceData, hco_ptr, vco_ptr)
    !
    ! :Purpose: Calculate analysis-error variance.
    !
    implicit none

    ! Arguments
    type(struct_obs), intent(in)  :: obsSpaceData
    type(struct_hco), pointer  :: hco_ptr
    type(struct_vco), pointer  :: vco_ptr

    ! Local Variables
    type(struct_gsv) :: stateVectorBkGnd
    type(struct_gsv) :: stateVectorAnal
    integer, allocatable :: numObs(:,:)
    real(8), pointer :: GLBE_ptr(:,:,:,:), GLAE_ptr(:,:,:,:)

    real(8), allocatable :: H(:,:), B(:,:), PHiA(:), A(:,:), &
                            R(:), PH(:,:), KH(:), IKH(:), Ainv(:,:)
    integer, parameter :: maxvar = 15000

    integer :: statei(maxvar)    ! Model grid coordinate i of neighbors
    integer :: statej(maxvar)    ! Model grid coordinate j of neighbors

    integer :: nmodel, imodel, jmodel, varind, errorflag

    logical :: found

    integer :: stepIndex, levIndex, lonIndex, latIndex, gridIndex, headerIndex, bodyIndex, bufrCode, kIndex, procIndex
    integer :: influentObsIndex2, influentObsIndex

    real(8) :: scaling
    integer :: iObs, jObs

    integer :: numInfluentObs, xIndex1, yIndex1, xIndex2, yIndex2
    integer :: gridptCount, gridpt

    ! max_avar is the maximum analysis-error variance
    real(8) :: distance, max_avar

    real(kdkind), allocatable :: positionArray(:,:)
    real(8),      allocatable :: latInRad(:,:), lonInRad(:,:)

    character(len=*), parameter :: myName = 'aer_analysisError'
    character(len=*), parameter :: correlationLengthFileName = './bgstddev'
    type(struct_gsv)            :: statevector
    real(4), pointer            :: field3D_r4_ptr(:,:,:)

    character(len=2 ) :: typvar
    character(len=12) :: etiket
    logical :: containsFullField

    type(struct_columnData) :: column
    type(struct_columnData) :: columng

    integer          :: idate, imonth
    integer          :: trackCellNum
    character(len=8) :: ccyymmdd

    real(8) :: leadTimeInHours

    if( mpi_myid > 0 ) return

    write(*,*) '**********************************************************'
    write(*,*) '** '//myName//': Calculate analysis-error variance **'
    write(*,*) '**********************************************************'

    call gsv_allocate( stateVectorBkGnd, 1, hco_ptr, vco_ptr, dateStamp_opt=-1, mpi_local_opt=.true., mpi_distribution_opt='Tiles', varNames_opt=(/'GLE'/), dataKind_opt=8 )
    call gsv_allocate( stateVectorAnal,  1, hco_ptr, vco_ptr, dateStamp_opt=-1, varNames_opt=(/'GLE'/), dataKind_opt=8 )

    etiket = '            '
    typvar = 'P@'
    containsFullField = .true.

    call gsv_readFromFile( stateVectorBkGnd, './trlm_02', etiket, typvar )

    leadTimeInHours = real(stateVectorBkGnd%deet*stateVectorBkGnd%npasList(1),8)/3600.0d0
    call incdatr(stateVectorAnal%dateOriginList(1), stateVectorBkGnd%dateOriginList(1), leadTimeInHours)

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

    write(*,*) myName//': Correlation length scale 2D field will be read from the file: ', correlationLengthFileName
    call gsv_allocate( statevector, 1, hco_ptr, vco_ptr, dateStamp_opt=-1, dataKind_opt=4, &
                       hInterpolateDegree_opt='LINEAR', varNames_opt=(/'GL'/) )
    call gsv_zero( statevector )
    call gsv_readFromFile(statevector, correlationLengthFileName, 'CORRLEN', ' ', unitConversion_opt = .false. )

    call gsv_getField( statevector, field3D_r4_ptr, 'GL' )
    Lcorr(:,:) = real(field3D_r4_ptr( :, :, 1 ), 8)
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
        latInRad(lonIndex,latIndex) = real(stateVectorBkGnd % hco % lat2d_4(lonIndex,latIndex), 8)
        lonInRad(lonIndex,latIndex) = real(stateVectorBkGnd % hco % lon2d_4(lonIndex,latIndex), 8)

        positionArray(:,gridIndex) = kdtree2_3dPosition(lonInRad(lonIndex,latIndex), latInRad(lonIndex,latIndex))

      end do
    end do

    write(*,*) myName//': latInRad min/max: ', minval(latInRad(:,:) ), maxval( latInRad(:,:) )
    write(*,*) myName//': lonInRad min/max: ', minval(lonInRad(:,:) ), maxval( lonInRad(:,:) )

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

    max_avar = 1.0

    call gsv_getField(stateVectorBkGnd, GLBE_ptr, 'GLE')
    call gsv_getField(stateVectorAnal,  GLAE_ptr, 'GLE')

    ! Initialisation
    do stepIndex = 1, stateVectorBkGnd%numStep
      do levIndex = 1, gsv_getNumLev(stateVectorBkGnd,vnl_varLevelFromVarname('GLE'))
        do latIndex = 1, nj
          do lonIndex = 1, ni
            GLAE_ptr(lonIndex,latIndex,levIndex,stepIndex) = GLBE_ptr(lonIndex,latIndex,levIndex,stepIndex)
          end do
        end do
      end do
    end do

    ! Only variables assigned within or by the loop can be private.

    STEP: do stepIndex = 1, stateVectorBkGnd%numStep
      LEVEL: do levIndex = 1, gsv_getNumLev(stateVectorBkGnd,vnl_varLevelFromVarname('GLE'))

!$omp parallel do default(shared) schedule(dynamic) private(H,B,PHiA,A,Ainv,R,PH,KH,IKH, &
!$omp        statei, statej, errorflag, &
!$omp        nmodel, imodel, jmodel, varind, found, influentObsIndex2, influentObsIndex, &
!$omp        scaling, iObs, jObs, lonIndex, latIndex, headerIndex, bodyIndex, bufrCode, kIndex, procIndex, &
!$omp        numInfluentObs, xIndex1, yIndex1, xIndex2, yIndex2, distance, &
!$omp        interpWeight, obsLatIndex, obsLonIndex, gridptCount, gridpt, idate, imonth, trackCellNum, ccyymmdd)
        YINDEX: do latIndex = 1, nj
          XINDEX: do lonIndex = 1, ni

            numInfluentObs = influentObs(lonIndex,latIndex)%numObs

            if ( numInfluentObs == 0 .or. .not. stateVectorBkGnd%oceanMask%mask(lonIndex,latIndex,levIndex) ) cycle

            ! form the observation-error covariance (diagonal) matrix

            allocate(R(numInfluentObs))

            do influentObsIndex = 1, numInfluentObs
              R(influentObsIndex) = obs_bodyElem_r(obsSpaceData,OBS_OER,influentObs(lonIndex,latIndex)%bodyIndex(influentObsIndex))
            end do

            ! From standard deviation to variance
            R = R*R

            ! find all model variables involved here

            nmodel = 0
            statei(:) = 0
            statej(:) = 0
            do influentObsIndex = 1, numInfluentObs

              headerIndex = influentObs(lonIndex,latIndex)%headerIndex(influentObsIndex)
              bodyIndex   = influentObs(lonIndex,latIndex)%bodyIndex(influentObsIndex)

              bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )

              scaling = 0.0d0

              select case (bufrCode)
              case(BUFR_ICEC, BUFR_ICEP)
                scaling = 100.0d0
              case(BUFR_ICEV)
                scaling = 1.0d0
              case(BUFR_ICES)
                idate = obs_headElem_i( obsSpaceData, OBS_DAT, headerIndex ) 
                write(ccyymmdd, FMT='(i8.8)') idate
                read(ccyymmdd(5:6), FMT='(i2)') imonth
                trackCellNum = obs_headElem_i( obsSpaceData, OBS_FOV, headerIndex )
                scaling = oer_ascatAnisIce(trackCellNum,imonth) - oer_ascatAnisOpenWater(trackCellNum,imonth)
              case default
                write(*,*) 'bufrCode = ',bufrCode
                call utl_abort( myName//': bufrCode not an ice code')
              end select

              if ( scaling /= 0.0d0 ) then

                do kIndex = stateVectorBkGnd%mykBeg, stateVectorBkGnd%mykEnd
                  do procIndex = 1, mpi_nprocs

                    call s2c_getWeightsAndGridPointIndexes(headerIndex, &
                         kIndex, stepIndex, procIndex, interpWeight, obsLatIndex, obsLonIndex, gridptCount)

                    do gridpt = 1, gridptCount

                      if ( interpWeight(gridpt) /= 0.0d0 ) then

                        iObs = obsLonIndex(gridpt)
                        jObs = obsLatIndex(gridpt)

                        found = .false.
                        do jmodel=1,nmodel
                          if(iObs == statei(jmodel) .and. &
                             jObs == statej(jmodel)) then
                            found = .true.
                            exit
                          end if
                        end do
                        if(.not. found) then
                          nmodel = nmodel + 1
                          if(nmodel > maxvar) then
                            call utl_abort( myName//': Structure state'// &
                                 ' too small in subroutine'// &
                                 ' analysis_error_mod'// &
                                 ' Increase maxvar')
                          end if
                          statei(nmodel) = iObs
                          statej(nmodel) = jObs
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
            do jmodel=1,nmodel
              if(lonIndex == statei(jmodel) .and. &
                 latIndex == statej(jmodel)) then
                varind = jmodel
                found = .true.
                exit
              end if
            end do
            if(.not. found) then
              nmodel = nmodel + 1
              if(nmodel > maxvar) then
                call utl_abort(myName//': Structure state too small in'// &
                     ' subroutine analysis_error_mod'// &
                     ' increase maxvar')
              end if
              varind = nmodel
              statei(nmodel) = lonIndex
              statej(nmodel) = latIndex
            end if

            ! form the observation operator matrix (H)

            allocate(H(nmodel,numInfluentObs))

            H(:,:) = 0.0d0

            jmodel = 0

            do influentObsIndex = 1, numInfluentObs

              headerIndex = influentObs(lonIndex,latIndex)%headerIndex(influentObsIndex)
              bodyIndex   = influentObs(lonIndex,latIndex)%bodyIndex(influentObsIndex)

              bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )

              scaling = 0.0d0

              select case (bufrCode)
              case(BUFR_ICEC, BUFR_ICEP)
                scaling = 100.0d0
              case(BUFR_ICEV)
                scaling = 1.0d0
              case(BUFR_ICES)
                idate = obs_headElem_i( obsSpaceData, OBS_DAT, headerIndex ) 
                write(ccyymmdd, FMT='(i8.8)') idate
                read(ccyymmdd(5:6), FMT='(i2)') imonth
                trackCellNum = obs_headElem_i( obsSpaceData, OBS_FOV, headerIndex )
                scaling = oer_ascatAnisIce(trackCellNum,imonth) - oer_ascatAnisOpenWater(trackCellNum,imonth)
              case default
                 write(*,*) 'bufrCode = ',bufrCode
                 call utl_abort( myName//': bufrCode not an ice code')
              end select

              if ( scaling /= 0.0d0 ) then

                do kIndex = stateVectorBkGnd%mykBeg, stateVectorBkGnd%mykEnd
                  do procIndex = 1, mpi_nprocs

                    call s2c_getWeightsAndGridPointIndexes(headerIndex, &
                         kIndex, stepIndex, procIndex, interpWeight, obsLatIndex, obsLonIndex, gridptCount)

                    do gridpt = 1, gridptCount

                      if ( interpWeight(gridpt) /= 0.0d0 ) then

                        iObs = obsLonIndex(gridpt)
                        jObs= obsLatIndex(gridpt)

                        found = .false.
                        do while (.not. found)
                          if(jmodel < nmodel) then
                            jmodel = jmodel + 1
                          else
                            jmodel = 1
                          end if
                          if(iObs == statei(jmodel) .and. &
                             jObs == statej(jmodel)) then
                            found = .true.
                            H(jmodel,influentObsIndex) = scaling*interpWeight(gridpt)
                          end if
                        end do
                        if(.not. found) then
                          write(*,*) 'iObs = ',iObs
                          write(*,*) 'jObs = ',jObs
                          write(*,*) 'lonIndex = ',lonIndex
                          write(*,*) 'latIndex = ',latIndex
                          write(*,*) 'gridptCount = ',gridptCount
                          write(*,*) 'gridpt = ',gridpt
                          write(*,*) 'nmodel = ',nmodel
                          do jmodel=1,nmodel
                            write(*,*) 'jmodel statei statej = ',statei(jmodel),statej(jmodel)
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

            allocate(B(nmodel,nmodel))

            B = 0.0d0

            do jmodel = 1, nmodel

              xIndex2 = statei(jmodel)
              yIndex2 = statej(jmodel)

              do imodel = jmodel, nmodel

                xIndex1 = statei(imodel)
                yIndex1 = statej(imodel)

                if(xIndex2 == xIndex1 .and. yIndex2 == yIndex1) then
                  B(imodel,jmodel) = GLBE_ptr(xIndex1,yIndex1,levIndex,stepIndex)**2
                else if(Lcorr(xIndex2,yIndex2) > 0.0d0) then
                  distance = 1.0d-3 * phf_calcDistance(latInRad(xIndex2,yIndex2), &
                                                       lonInRad(xIndex2,yIndex2), &
                                                       latInRad(xIndex1,yIndex1), &
                                                       lonInRad(xIndex1,yIndex1))
                  B(imodel,jmodel) = GLBE_ptr(xIndex2,yIndex2,levIndex,stepIndex)* &
                                     GLBE_ptr(xIndex1,yIndex1,levIndex,stepIndex) &
                                     *exp(-0.5*(distance/Lcorr(xIndex2,yIndex2))**2)
                end if

                ! symmetric matrix !
                B(jmodel,imodel) = B(imodel,jmodel)

              end do

            end do

            ! form the observation background error covariance matrix (PHT)

            allocate(PH(numInfluentObs,nmodel))

            ! PH = matmul (B, transpose(H))
            do imodel = 1, nmodel
              do influentObsIndex = 1, numInfluentObs
                PH(influentObsIndex,imodel) = dot_product(B(:,imodel), H(:,influentObsIndex))
              end do
            end do

            ! form the error covariance matrix background field (HPHT)

            allocate(A(numInfluentObs,numInfluentObs),Ainv(numInfluentObs,numInfluentObs))

            do influentObsIndex = 1, numInfluentObs

              do influentObsIndex2 = 1, numInfluentObs
                A(influentObsIndex2,influentObsIndex) = dot_product(H(:,influentObsIndex), PH(influentObsIndex2,:))
              end do

            end do

            !  covariance matrix of the innovation (HPHT + R)

            do influentObsIndex = 1, numInfluentObs
              A(influentObsIndex,influentObsIndex) = A(influentObsIndex,influentObsIndex) + R(influentObsIndex)
            end do

            ! Inverse of the covariance matrix of the innovation

            call FINDInv(A, Ainv, numInfluentObs, errorflag)

            ! Kalman gain; this is the row corresponding to the analysis variable

            allocate(PHiA(numInfluentObs))

            do influentObsIndex = 1, numInfluentObs
              PHiA(influentObsIndex) = dot_product(PH(:,varind), Ainv(influentObsIndex,:))
            end do

            ! compute the error variance of the analysis

            allocate(KH(nmodel))

            ! KH = matmul (PHiA, H)
            do imodel = 1, nmodel
              KH(imodel) = dot_product(PHiA, H(imodel,:))
            end do

            ! IKH = I - KH

            allocate(IKH(nmodel))

            IKH = -KH

            IKH(varind) = 1.0d0 - KH(varind)

            GLAE_ptr(lonIndex,latIndex,levIndex,stepIndex) = sqrt( dot_product (IKH, B(:,varind)) )

            if(GLAE_ptr(lonIndex,latIndex,levIndex,stepIndex) < 0.0) then
              write(*,*) myName//'negative analysis-error variance = ', &
                   GLAE_ptr(lonIndex,latIndex,levIndex,stepIndex), ' reset to zero at grid point (',lonIndex,latIndex,')'
              GLAE_ptr(lonIndex,latIndex,levIndex,stepIndex) = max(GLAE_ptr(lonIndex,latIndex,levIndex,stepIndex), 0.0)
            end if

            if(GLAE_ptr(lonIndex,latIndex,levIndex,stepIndex) > GLBE_ptr(lonIndex,latIndex,levIndex,stepIndex)) then
              write(*,*) myName//'analysis-error variance = ', &
                   GLAE_ptr(lonIndex,latIndex,levIndex,stepIndex), ' is larger than background-error ', &
                   'variance and kept at = ',GLBE_ptr(lonIndex,latIndex,levIndex,stepIndex)
              GLAE_ptr(lonIndex,latIndex,levIndex,stepIndex) = min(GLBE_ptr(lonIndex,latIndex,levIndex,stepIndex), max_avar)
            end if

            deallocate(B, R, A, Ainv, H, PH, PHiA, KH, IKH)
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

    call gsv_writeToFile( stateVectorAnal, './anlm_000m', '', typvar_opt='A@', containsFullField_opt=.true. )

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
    integer :: headerIndex, bodyIndexBeg, bodyIndexEnd, bodyIndex, kIndex, stepIndex, procIndex
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

    do headerIndex = 1, obs_numHeader(obsSpaceData)
      bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1          
      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) then

          footprintRadius_r4 = s2c_getFootprintRadius(obsSpaceData, stateVectorBkGnd, headerIndex)
          influenceRadius_r4 = max(0.0, footprintRadius_r4) + 1000.0*maxLcorr

          if ( maxLcorr == 0.0d0 ) then

            do kIndex = stateVectorBkGnd%mykBeg, stateVectorBkGnd%mykEnd
              do stepIndex = 1, stateVectorBkGnd%numStep
                do procIndex = 1, mpi_nprocs

                  call s2c_getWeightsAndGridPointIndexes(headerIndex, &
                       kIndex, stepIndex, procIndex, interpWeight, obsLatIndex, obsLonIndex, gridptCount)

                  do gridpt = 1, gridptCount

                    if ( interpWeight(gridpt) /= 0.0d0 ) then

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

                    end if

                  end do

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
              call utl_abort('findObs: the parameter maxNumLocalGridptsSearch must be increased')
            end if

            do resultsIndex = 1, numLocalGridptsFoundSearch

              gridIndex = searchResults(resultsIndex)%idx
              if ( gridIndex < 1 .or. gridIndex > stateVectorBkGnd%hco%ni * stateVectorBkGnd%hco%nj ) then
                write(*,*) 'analysisErrorStdDev: gridIndex=', gridIndex
                call utl_abort('findObs: gridIndex out of bound.')
              end if

              latIndex = (gridIndex - 1) / stateVectorBkGnd%hco%ni + 1
              lonIndex = gridIndex - (latIndex - 1) * stateVectorBkGnd%hco%ni
              if ( lonIndex < 1 .or. lonIndex > stateVectorBkGnd%hco%ni .or. &
                   latIndex < 1 .or. latIndex > stateVectorBkGnd%hco%nj ) then
                write(*,*) 'analysisErrorStdDev: lonIndex=', lonIndex, ',latIndex=', latIndex
                call utl_abort('findObs: lonIndex/latIndex out of bound.')
              end if

              numObs(lonIndex,latIndex) = numObs(lonIndex,latIndex) + 1
              if(associated(influentObs(lonIndex,latIndex)%bodyIndex)) then
                if(numObs(lonIndex,latIndex) > influentObs(lonIndex,latIndex)%numObs) then
                  call utl_abort('findObs: Array too small in subroutine findObs')
                end if
                influentObs(lonIndex,latIndex)%headerIndex(numObs(lonIndex,latIndex)) = headerIndex
                influentObs(lonIndex,latIndex)%bodyIndex(numObs(lonIndex,latIndex)) = bodyIndex
              end if

            end do

          end if

        end if

      end do

    end do

    call tmg_stop(189)

  end subroutine findObs

  !--------------------------------------------------------------------------
  ! FINDInv
  !--------------------------------------------------------------------------
  SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
    !
    ! :Purpose: Subroutine to find the inverse of a square matrix
    !           Author : Louisda16th a.k.a Ashwith J. Rego
    !           Reference : Algorithm has been well explained in:
    !           http://math.uww.edu/~mcfarlat/inverse.htm           
    !           http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
    !           Source code found at https://github.com/ashwith/workspace/blob/master/fortran/matrixinverse.f90
    !           Source code found in http://www.dreamincode.net/code/snippet1272.htm
    !
    implicit none

    ! arguments:
    integer, INTENT(IN) :: n
    integer, INTENT(OUT) :: errorflag ! Return error status. -1 for error, 0 for normal
    real(8), INTENT(IN), DIMENSION(n,n) :: matrix ! Input matrix
    real(8), INTENT(OUT), DIMENSION(n,n) :: inverse ! Inverted matrix

    ! locals:
    LOGICAL :: FLAG = .TRUE.

    integer :: i, j, k

    ! Double precision is required here to avoid round-off errors.

    real(8) :: m

    real(8), DIMENSION(n,2*n) :: augmatrix !augmented matrix

    ! Augment input matrix with an identity matrix

    DO i = 1, n

      DO j = 1, 2*n

        IF (j <= n ) THEN

          augmatrix(i,j) = matrix(i,j)

        ELSE IF ((i+n) == j) THEN

          augmatrix(i,j) = 1

        Else

          augmatrix(i,j) = 0

        END IF

      END DO

    END DO

    ! Reduce augmented matrix to upper triangular form

    DO k =1, n-1

      IF (augmatrix(k,k) == 0) THEN

        FLAG = .FALSE.

        DO i = k+1, n

          IF (augmatrix(i,k) /= 0) THEN

            DO j = 1,2*n

              augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)

            END DO

            FLAG = .TRUE.

            EXIT

          END IF

          IF (FLAG .EQV. .FALSE.) THEN

            write(*,*) 'matrix(:,',k,') = ',matrix(:,k)
            write(*,*) 'augmatrix(:,',k,') = ',augmatrix(:,k)
            write(*,*) 'n = ',n
            call utl_abort("FINDinv: Matrix is non - invertible test 1")

            inverse = 0

            errorflag = -1

          END IF

        END DO

      END IF

      DO j = k+1, n                       

        m = augmatrix(j,k)/augmatrix(k,k)

        DO i = k, 2*n

          augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)

        END DO

      END DO

    END DO

    ! Test for invertibility

    DO i = 1, n

      IF (augmatrix(i,i) == 0) THEN

        write(*,*) 'augmatrix(',i,',',i,') = ',augmatrix(i,i)
        write(*,*) 'n = ',n
        call utl_abort("FINDinv: Matrix is non - invertible test 2")

        inverse = 0

        errorflag = -1

      END IF

    END DO

    ! Make diagonal elements as 1

    DO i = 1 , n

      m = augmatrix(i,i)

      DO j = i , (2 * n)                               

        augmatrix(i,j) = (augmatrix(i,j) / m)

      END DO

    END DO

    ! Reduced right side half of augmented matrix to identity matrix

    DO k = n-1, 1, -1

      DO i =1, k

        m = augmatrix(i,k+1)

        DO j = k, (2*n)

          augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m

        END DO

      END DO

    END DO

    ! store answer

    DO i =1, n

      DO j = 1, n

        inverse(i,j) = augmatrix(i,j+n)

      END DO

    END DO

    errorflag = 0

  END SUBROUTINE FINDinv

end module analysisError_mod

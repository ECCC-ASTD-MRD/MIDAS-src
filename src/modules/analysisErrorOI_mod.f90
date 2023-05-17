
module analysisErrorOI_mod
  ! MODULE analysisErrorOI (prefix='aer' category='1. High-level functionality')
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
  use gridStateVectorFileIO_mod
  use kdtree2_mod
  use midasMpi_mod
  use obsSpaceData_mod
  use physicsFunctions_mod
  use stateToColumn_mod
  use varNamelist_mod
  use utilities_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use obsOperators_mod
  use message_mod
  use oceanMask_mod
  use timeCoord_mod

  implicit none

  private

  ! public subroutines and functions
  public :: aer_analysisError

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

  integer, parameter :: maxNumLocalGridptsSearch = 1500

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
    type(struct_obs), intent(in) :: obsSpaceData
    type(struct_hco), pointer    :: hco_ptr
    type(struct_vco), pointer    :: vco_ptr

    ! Local Variables
    integer :: fnom, fclos, nulnam, ierr

    type(struct_gsv) :: stateVectorAnlErrorStd    ! state vector for analysis error std deviation
    type(struct_gsv) :: stateVectorTrlErrorStd    ! state vector for background error std deviation
    real(8), pointer :: anlErrorStdDev_ptr(:,:,:,:) ! pointer for analysis error std deviation
    real(8), pointer :: trlErrorStdDev_ptr(:,:,:,:) ! pointer for background error std deviation
    integer, allocatable :: numObs(:,:)

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

    character(len=4), pointer :: analysisVariable(:)
    type(struct_gsv)          :: statevectorLcorr
    real(4), pointer          :: field3D_r4_ptr(:,:,:)
    type(struct_columnData)   :: column
    type(struct_columnData)   :: columng
    type(struct_ocm)          :: oceanMask
    real(8) :: leadTimeInHours

    character(len=20), parameter :: errorStddev_input  = 'errorstdev_in'        ! input  filename for anl ot trl error standard deviation
    character(len=20), parameter :: anlErrorStddev_output = 'anlerrorstdev_out' ! output filename for anl error std deviation
    character(len=20), parameter :: trlErrorStddev_output = 'trlerrorstdev_out' ! output filename for trl (background) error std deviation

    ! namelist variables:
    real(8)           :: maxAnalysisErrorStdDev ! maximum limit imposed on analysis error stddev
    logical           :: propagateAnalysisError ! propagate analysis error or not
    logical           :: propagateDSLO          ! propagate Days Since Last Obs field or not
    real(4)           :: errorGrowth            ! seaice: fraction of ice per hour, SST: estimated growth
    character(len=12) :: analysisEtiket         ! analysis field etiket in a standard file
    character(len=12) :: analErrorStdEtiket     ! analysis error standard deviation field etiket in the input/output standard files
    character(len=12) :: bckgErrorStdEtiket     ! background error standard deviation field etiket in the input/output standard files
    integer           :: hoursSinceLastAnalysis ! number of hours since the last analysis
    logical           :: saveTrlStdField        ! to save trial standard deviation field or not
    character(len=2)  :: inputTypeVar           ! typvar of the analysis error field in the input file 
    namelist /namaer/ maxAnalysisErrorStdDev, propagateAnalysisError, propagateDSLO, &
                      errorGrowth, analysisEtiket, analErrorStdEtiket, &
                      bckgErrorStdEtiket, hoursSinceLastAnalysis, saveTrlStdField, inputTypeVar

    if(mmpi_nprocs > 1) then
      write(*,*) 'mmpi_nprocs = ', mmpi_nprocs
      call utl_abort('aer_analysisError: this version of the code should only be used with one mpi task.')
    end if
    if(mmpi_myid > 0) return

    write(*,*) '**********************************************************'
    write(*,*) '** aer_analysisError: Calculate analysis-error variance **'
    write(*,*) '**********************************************************'

    ! default namelist variable values
    maxAnalysisErrorStdDev = 1.0d0
    propagateAnalysisError = .false.
    propagateDSLO = .false.
    errorGrowth = 1.0
    analysisEtiket = ''
    analErrorStdEtiket = 'A-ER STD DEV'
    bckgErrorStdEtiket = 'B-ER STD DEV'
    hoursSinceLastAnalysis = 6
    saveTrlStdField = .false.
    inputTypeVar = 'A@'
    
    ! read the namelist
    if (.not. utl_isNamelistPresent('namaer','./flnml')) then
      if (mmpi_myid == 0) then
        call msg('aer_analysisError:', ' namaer is missing in the namelist.')
        call msg('aer_analysisError:', ' the default values will be taken.')
      end if
    else
      ! reading namelist variables
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam, nml = namaer, iostat = ierr)
      if (ierr /= 0) call utl_abort('aer_analysisError:: Error reading namelist')
      ierr = fclos(nulnam)
    end if
    write(*, nml = namaer)

    nullify(analysisVariable)
    call gsv_varNamesList(analysisVariable)
    
    if (size(analysisVariable) > 1) then
      call utl_abort('aer_analysisError: Check namelist NAMSTATE. analysisVariable is greater than 1.')
    end if

    if (analysisVariable(1) == 'GL') then
      call msg('aer_analysisError:', ' computing seaice analysis error...')
    else if(analysisVariable(1) == 'TM') then
      call msg('aer_analysisError:', ' computing SST analysis error...')
    else
      call utl_abort('aer_analysisError:: The current code does not work with '&
                     //trim(analysisVariable(1))//' analysis variable.')
    end if

    call gsv_allocate(stateVectorAnlErrorStd, 1, hco_ptr, vco_ptr, dateStamp_opt = -1, &
                      mpi_local_opt = .true., mpi_distribution_opt = 'Tiles', &
                      varNames_opt = (/analysisVariable(1)/), dataKind_opt = 8)
    call gsv_getField(stateVectorAnlErrorStd, anlErrorStdDev_ptr)
    call gsv_allocate(stateVectorTrlErrorStd, 1, hco_ptr, vco_ptr, dateStamp_opt = -1, &
                      mpi_local_opt = .true., mpi_distribution_opt = 'Tiles', &
                      varNames_opt = (/analysisVariable(1)/), dataKind_opt = 8)
    call gsv_getField(stateVectorTrlErrorStd, trlErrorStdDev_ptr)

    if (propagateAnalysisError) then     
      call msg('aer_analysisError:', &
               ' analysis error std field is read from: '//trim(errorStddev_input))    
      call gio_readFromFile(stateVectorAnlErrorStd, errorStddev_input, ' ', &
                            inputTypeVar, containsFullField_opt = .false.)
      
      ! initialize trl error std deviation field:
      trlErrorStdDev_ptr(:,:,:,:) = anlErrorStdDev_ptr(:,:,:,:)

      call ocm_readMaskFromFile (oceanMask, hco_ptr, vco_ptr, errorStddev_input)
      call aer_propagateAnalysisError (stateVectorTrlErrorStd, oceanMask, &
                                       analysisVariable(1), &
                                       analysisEtiket, errorGrowth, &
                                       hco_ptr, vco_ptr)

      if (saveTrlStdField) then
        ! zap analysis error etiket with background error etiket
        stateVectorTrlErrorStd%etiket = bckgErrorStdEtiket
        
        ! copy mask from analysis error std deviation field to trl error std field
        call gsv_copyMask(stateVectorAnlErrorStd, stateVectorTrlErrorStd)
        
        ! update dateStamp from env variable
        call gsv_modifyDate(stateVectorTrlErrorStd, tim_getDateStamp(), &
                            modifyDateOrigin_opt = .true.)

        ! save background error (increased analysis error) standard deviation field
        call gio_writeToFile(stateVectorTrlErrorStd, trlErrorStddev_output, &
                             bckgErrorStdEtiket, typvar_opt = inputTypeVar, &
                             containsFullField_opt = .false.)
      end if

    else
      call msg('aer_analysisError:', &
               ' trial error std field is read from: '//trim(errorStddev_input))
      call gio_readFromFile(stateVectorTrlErrorStd, errorStddev_input, ' ', &
                            inputTypeVar, containsFullField_opt = .false.)
      call gsv_copyMask(stateVectorTrlErrorStd, stateVectorAnlErrorStd)
    end if

    leadTimeInHours = real(stateVectorTrlErrorStd%deet * stateVectorTrlErrorStd%npasList(1), 8) / 3600.0d0
    call incdatr(stateVectorAnlErrorStd%dateOriginList(1), &
                 stateVectorTrlErrorStd%dateOriginList(1), leadTimeInHours)

    stateVectorAnlErrorStd%etiket = analErrorStdEtiket

    call col_setVco(column, vco_ptr)
    call col_allocate(column,  obs_numHeader(obsSpaceData), varNames_opt=(/analysisVariable(1)/))
    call col_setVco(columng, vco_ptr)
    call col_allocate(columng, obs_numHeader(obsSpaceData), varNames_opt=(/analysisVariable(1)/))
    call s2c_tl(stateVectorTrlErrorStd, column, columng, obsSpaceData)

    ni = stateVectorTrlErrorStd%hco%ni
    nj = stateVectorTrlErrorStd%hco%nj

    allocate(lonInRad(ni, nj))
    allocate(latInRad(ni, nj))
    allocate(Lcorr(ni, nj))

    ! get correlation length scale field
    write(*,*) 'aer_analysisError: get correlation length scale field...'
    call gsv_allocate(statevectorLcorr, 1, hco_ptr, vco_ptr, dateStamp_opt = -1, &
                      dataKind_opt = 4, hInterpolateDegree_opt = 'LINEAR', &
                      varNames_opt = (/analysisVariable(1)/))
    call gsv_zero(statevectorLcorr)
    call gio_readFromFile(statevectorLcorr, './bgstddev', 'CORRLEN', ' ', &
                          unitConversion_opt = .false.)
    call gsv_getField(statevectorLcorr, field3D_r4_ptr, analysisVariable(1))
    ! Convert from km to meters
    Lcorr(:,:) = 1000.0d0 * real(field3D_r4_ptr(:, :, 1), 8)

    write(*,*) 'aer_analysisError: min/max correlation length scale 2D field: ', &
               minval(Lcorr(:,:) ), maxval(Lcorr(:,:))
    call gsv_deallocate(statevectorLcorr)

    ! create kdtree
    write(*,*) 'aer_analysisError: start creating kdtree for stateVectorTrlErrorStd'
    write(*,*) 'Memory Used: ', get_max_rss() / 1024, 'Mb'

    allocate(positionArray(3, ni * nj))

    gridIndex = 0
    do latIndex = 1, nj
      do lonIndex = 1, ni
        gridIndex = gridIndex + 1
        latInRad(lonIndex, latIndex) = real(stateVectorTrlErrorStd%hco%lat2d_4(lonIndex, latIndex), 8)
        lonInRad(lonIndex, latIndex) = real(stateVectorTrlErrorStd%hco%lon2d_4(lonIndex, latIndex), 8)
        positionArray(:, gridIndex) = kdtree2_3dPosition(lonInRad(lonIndex, latIndex), &
                                                        latInRad(lonIndex, latIndex))
      end do
    end do

    write(*,*) 'aer_analysisError: latInRad min/max: ', minval(latInRad(:,:)), maxval(latInRad(:,:))
    write(*,*) 'aer_analysisError: lonInRad min/max: ', minval(lonInRad(:,:)), maxval(lonInRad(:,:))

    nullify(tree)
    tree => kdtree2_create(positionArray, sort = .false., rearrange = .true.) 

    deallocate(positionArray)

    write(*,*) 'aer_analysisError: done creating kdtree for stateVectorTrlErrorStd'
    write(*,*) 'Memory Used: ', get_max_rss() / 1024, 'Mb'

    ! Go through all observations a first time to get
    ! the number of influential observations
    ! in order to allocate memory appropriately.

    allocate(influentObs(ni, nj))
    allocate(numObs(ni, nj))

    call msg('aer_analysisError:', ' looking for observations...')
    call findObs(obsSpaceData, stateVectorTrlErrorStd, numObs, analysisVariable(1))

    ! Memory allocation

    do latIndex = 1, nj
      do lonIndex = 1, ni
        influentObs(lonIndex, latIndex)%numObs = numObs(lonIndex,latIndex)
        allocate(influentObs(lonIndex, latIndex)%headerIndex(numObs(lonIndex, latIndex)))
        allocate(influentObs(lonIndex, latIndex)%bodyIndex(numObs(lonIndex, latIndex)))
      end do
    end do

    ! Go through all observations a second time to get
    ! the indexes of the influential observations.
    ! This is not efficient to go through all observations twice
    ! but it saves lot of memory space.

    call msg('aer_analysisError:', ' go through all observations a second time...')
    call findObs(obsSpaceData, stateVectorTrlErrorStd, numObs, analysisVariable(1))

    deallocate(numObs)

    write(*,*) 'aer_analysisError:: obs found'
    write(*,*) 'Memory Used: ',get_max_rss() / 1024, 'Mb'

    ! Calculate analysis-error one analysis variable (grid point) at a time

    ! Initialisation
    anlErrorStdDev_ptr(:,:,:,:) = trlErrorStdDev_ptr(:,:,:,:)
    call msg('aer_analysisError:', ' analysis error std field initialization...OK')

    ! Only variables assigned within or by the loop can be private.
    STEP: do stepIndex = 1, stateVectorTrlErrorStd%numStep
      LEVEL: do levIndex = 1, gsv_getNumLev(stateVectorTrlErrorStd, &
                                            vnl_varLevelFromVarname(analysisVariable(1)))

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

            numInfluentObs = influentObs(lonIndex, latIndex)%numObs

            if (numInfluentObs == 0 .or. &
                .not. stateVectorTrlErrorStd%oceanMask%mask(lonIndex, latIndex, levIndex)) cycle XINDEX

            ! form the observation-error covariance (diagonal) matrix

            allocate(obsErrorVariance(numInfluentObs))

            do influentObsIndex = 1, numInfluentObs
              bodyIndex = influentObs(lonIndex, latIndex)%bodyIndex(influentObsIndex)
              obsErrorVariance(influentObsIndex) = (obs_bodyElem_r(obsSpaceData, OBS_OER, &
                                                                   bodyIndex))**2
            end do

            ! find all model variables involved here

            numVariables = 0
            statei(:) = 0
            statej(:) = 0

            INFLUENTOBSCYCLE: do influentObsIndex = 1, numInfluentObs

              headerIndex = influentObs(lonIndex, latIndex)%headerIndex(influentObsIndex)
              bodyIndex   = influentObs(lonIndex, latIndex)%bodyIndex(influentObsIndex)

              if (analysisVariable(1) == 'GL') then
                scaling = oop_iceScaling(obsSpaceData, bodyIndex)
              else if (analysisVariable(1) == 'TM') then
                scaling = 1.0d0
              end if	

              if (scaling == 0.0d0) cycle INFLUENTOBSCYCLE
              KINDEXCYCLE: do kIndex = stateVectorTrlErrorStd%mykBeg, stateVectorTrlErrorStd%mykEnd
                PROCINDEXCYCLE: do procIndex = 1, mmpi_nprocs

                  call s2c_getWeightsAndGridPointIndexes(headerIndex, kIndex, stepIndex, procIndex, &
                                                         interpWeight, obsLatIndex, obsLonIndex, &
                                                         gridptCount)

                  GRIDPTCYCLE: do gridpt = 1, gridptCount

                    !if (interpWeight(gridpt) /= 0.0d0) then
                    if (interpWeight(gridpt) == 0.0d0) cycle GRIDPTCYCLE

                    xStateIndex = obsLonIndex(gridpt)
                    yStateIndex = obsLatIndex(gridpt)

                    found = .false.
                    do varIndex2=1, numVariables
                      if(xStateIndex == statei(varIndex2) .and. &
                         yStateIndex == statej(varIndex2)) then
                        found = .true.
                        exit
                      end if
                    end do

                    if (found) cycle GRIDPTCYCLE

                    numVariables = numVariables + 1
                    if(numVariables > maxvar) then
                      call utl_abort('aer_analysisError: Structure state'// &
                    	             ' too small in subroutine'// &
                    	             ' analysis_error_mod. Increase maxvar')
                    end if
                    statei(numVariables) = xStateIndex
                    statej(numVariables) = yStateIndex

                  end do GRIDPTCYCLE
                end do PROCINDEXCYCLE
              end do KINDEXCYCLE
            end do INFLUENTOBSCYCLE

            ! make sure that current analysis variable is part of the state
            ! vector even if it does not participate in obs calculation

            found = .false.
            do varIndex2 = 1, numVariables
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
                call utl_abort('aer_analysisError: Structure state too small in'// &
                     ' subroutine analysis_error_mod'// &
                     ' increase maxvar')
              end if
              currentAnalVarIndex = numVariables
              statei(numVariables) = lonIndex
              statej(numVariables) = latIndex
            end if

            ! form the observation operator matrix (obsOperator)

            allocate(obsOperator(numVariables, numInfluentObs))

            obsOperator(:,:) = 0.0d0

            varIndex2 = 0

            INFLUENTOBSCYCLE2: do influentObsIndex = 1, numInfluentObs

              headerIndex = influentObs(lonIndex, latIndex)%headerIndex(influentObsIndex)
              bodyIndex   = influentObs(lonIndex, latIndex)%bodyIndex(influentObsIndex)

              if (analysisVariable(1) == 'GL') then
                scaling = oop_iceScaling(obsSpaceData, bodyIndex)
              else if (analysisVariable(1) == 'TM') then
                scaling = 1.0d0
              end if

              if (scaling == 0.0d0) cycle INFLUENTOBSCYCLE2

              KINDEXCYCLE2: do kIndex = stateVectorTrlErrorStd%mykBeg, stateVectorTrlErrorStd%mykEnd
                PROCINDEXCYCLE2: do procIndex = 1, mmpi_nprocs

                  call s2c_getWeightsAndGridPointIndexes(headerIndex, kIndex, stepIndex, procIndex, &
                                                         interpWeight, obsLatIndex, obsLonIndex, &
                                                         gridptCount)

                  GRIDPTCYCLE2: do gridpt = 1, gridptCount

                    if (interpWeight(gridpt) == 0.0d0) cycle GRIDPTCYCLE2

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
                        obsOperator(varIndex2, influentObsIndex) = scaling * &
                        					   interpWeight(gridpt)
                      end if
                    end do

                    if(found) cycle GRIDPTCYCLE2

                    write(*,*) 'xStateIndex = ', xStateIndex
                    write(*,*) 'yStateIndex = ', yStateIndex
                    write(*,*) 'lonIndex = ', lonIndex
                    write(*,*) 'latIndex = ', latIndex
                    write(*,*) 'gridptCount = ', gridptCount
                    write(*,*) 'gridpt = ', gridpt
                    write(*,*) 'numVariables = ', numVariables
                    do varIndex2=1, numVariables
                      write(*,*) 'varIndex2 statei statej = ',statei(varIndex2), &
                    					      statej(varIndex2)
                    end do
                    call utl_abort('aer_analysisError: not found in state vect.')

                  end do GRIDPTCYCLE2
                end do PROCINDEXCYCLE2
              end do KINDEXCYCLE2
            end do INFLUENTOBSCYCLE2

            ! form the background-error covariance matrix

            allocate(Bmatrix(numVariables, numVariables))

            Bmatrix(:,:) = 0.0d0

            do varIndex2 = 1, numVariables

              xIndex2 = statei(varIndex2)
              yIndex2 = statej(varIndex2)

              do varIndex1 = varIndex2, numVariables

                xIndex1 = statei(varIndex1)
                yIndex1 = statej(varIndex1)

                if(xIndex2 == xIndex1 .and. yIndex2 == yIndex1) then
                  Bmatrix(varIndex1,varIndex2) = trlErrorStdDev_ptr(xIndex1, yIndex1, levIndex, stepIndex)**2
                else if(Lcorr(xIndex2,yIndex2) > 0.0d0) then
                  distance = phf_calcDistance(latInRad(xIndex2, yIndex2), lonInRad(xIndex2, yIndex2), &
                                              latInRad(xIndex1, yIndex1), lonInRad(xIndex1, yIndex1))
                  Bmatrix(varIndex1,varIndex2) = trlErrorStdDev_ptr(xIndex2, yIndex2, levIndex, stepIndex) * &
                                                 trlErrorStdDev_ptr(xIndex1, yIndex1, levIndex, stepIndex) * &
                                                 exp(-0.5 * (distance / Lcorr(xIndex2, yIndex2))**2)
                end if

                ! symmetric matrix !
                Bmatrix(varIndex2, varIndex1) = Bmatrix(varIndex1, varIndex2)

              end do

            end do

            ! form the observation background error covariance matrix (PHT)

            allocate(PH(numInfluentObs, numVariables))

            ! PH = matmul (Bmatrix, transpose(obsOperator))
            do varIndex1 = 1, numVariables
              do influentObsIndex = 1, numInfluentObs
                PH(influentObsIndex,varIndex1) = dot_product(Bmatrix(:, varIndex1), &
                                                 obsOperator(:, influentObsIndex))
              end do
            end do

            !  covariance matrix of the innovation (HPHT + R)

            allocate(innovCovariance(numInfluentObs, numInfluentObs), &
                     innovCovarianceInverse(numInfluentObs, numInfluentObs))

            do influentObsIndex = 1, numInfluentObs

              do influentObsIndex2 = 1, numInfluentObs
                if (influentObsIndex == influentObsIndex2) then
                  innovCovariance(influentObsIndex2, influentObsIndex) = &
                       dot_product(obsOperator(:, influentObsIndex), &
                       PH(influentObsIndex2,:)) + &
                       obsErrorVariance(influentObsIndex)
                else
                  innovCovariance(influentObsIndex2, influentObsIndex) = &
                       dot_product(obsOperator(:, influentObsIndex), &
                       PH(influentObsIndex2,:))
                end if
              end do

            end do

            ! Inverse of the covariance matrix of the innovation

            innovCovarianceInverse(:,:) = innovCovariance(:,:)
            call utl_matInverse(innovCovarianceInverse, numInfluentObs, &
                                printInformation_opt = .false.)

            ! Kalman gain; this is the row corresponding to the analysis variable

            allocate(PHiA(numInfluentObs))

            do influentObsIndex = 1, numInfluentObs
              PHiA(influentObsIndex) = dot_product(PH(:, currentAnalVarIndex), &
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

            anlErrorStdDev_ptr(lonIndex, latIndex, levIndex, stepIndex) = &
                 sqrt(dot_product (IKH, Bmatrix(:, currentAnalVarIndex)))

            if(anlErrorStdDev_ptr(lonIndex, latIndex, levIndex, stepIndex) < 0.0) then
              write(*,*) 'aer_analysisError: negative analysis-error Std dev. = ', &
                         anlErrorStdDev_ptr(lonIndex, latIndex, levIndex, stepIndex), &
                         ' reset to zero at grid point (',lonIndex, latIndex,')'
              anlErrorStdDev_ptr(lonIndex, latIndex, levIndex, stepIndex) = &
                   max(anlErrorStdDev_ptr(lonIndex, latIndex, levIndex, stepIndex), 0.0)
            end if

            if(anlErrorStdDev_ptr(lonIndex, latIndex, levIndex, stepIndex) > &
               trlErrorStdDev_ptr(lonIndex, latIndex, levIndex, stepIndex)) then
              write(*,*) 'aer_analysisError: analysis-error Std dev. = ', &
                         anlErrorStdDev_ptr(lonIndex, latIndex, levIndex, stepIndex), &
                         ' is larger than background-error ', &
                         'Std dev. is kept at = ', &
                         trlErrorStdDev_ptr(lonIndex, latIndex, levIndex, stepIndex)
              anlErrorStdDev_ptr(lonIndex, latIndex, levIndex, stepIndex) = &
                   min(trlErrorStdDev_ptr(lonIndex, latIndex, levIndex, stepIndex), &
                   maxAnalysisErrorStdDev)
            end if

            deallocate(Bmatrix)
            deallocate(obsErrorVariance)
            deallocate(innovCovariance)
            deallocate(innovCovarianceInverse)
            deallocate(obsOperator)
            deallocate(PH)
            deallocate(PHiA)
            deallocate(KH)
            deallocate(IKH)
            deallocate(influentObs(lonIndex, latIndex)%headerIndex)
            deallocate(influentObs(lonIndex, latIndex)%bodyIndex)

          end do XINDEX
        end do YINDEX
        !$omp end parallel do

      end do LEVEL
    end do STEP

    deallocate(influentObs)
    deallocate(Lcorr)
    deallocate(lonInRad)
    deallocate(latInRad)

    ! update dateStamp from env variable
    call gsv_modifyDate(stateVectorAnlErrorStd, tim_getDateStamp(), &
                        modifyDateOrigin_opt = .true.)

    ! save analysis error
    call msg('aer_analysisError:', ' writing analysis error std field to output file...')
    call gio_writeToFile(stateVectorAnlErrorStd, anlErrorStddev_output, &
                         analErrorStdEtiket, typvar_opt = inputTypeVar, &
                         containsFullField_opt = .false.)

    call col_deallocate(columng)
    call col_deallocate(column)
    call gsv_deallocate(stateVectorAnlErrorStd)
    call gsv_deallocate(stateVectorTrlErrorStd)

    if (analysisVariable(1) == 'GL') then
      ! Update the Days Since Last Obs
      call aer_daysSinceLastObs(obsSpaceData, hco_ptr, vco_ptr, &
                                errorStddev_input, anlErrorStddev_output, &
                                analysisVariable(1), propagateDSLO, &
                                hoursSinceLastAnalysis)
    end if
    
    call msg('aer_analysisError:', ' finished.')

  end subroutine aer_analysisError

  !---------------------------------------------------------
  ! findObs
  !---------------------------------------------------------
  subroutine findObs(obsSpaceData, stateVectorBkGnd, numObs, variableName)
    !
    !:Purpose: Find all observations used for the analysis.
    !
    implicit none

    ! Arguments
    type(struct_obs), intent(in)  :: obsSpaceData
    type(struct_gsv), intent(in)  :: stateVectorBkGnd
    integer,          intent(out) :: numObs(ni,nj)    ! number of observations found
    character(len=*), intent(in)  :: variableName     ! 'GL' for seaice or 'TM' for SST

    ! Local variables
    integer :: headerIndex, bodyIndexBeg, bodyIndexEnd, bodyIndex
    integer :: procIndex, kIndex, stepIndex
    integer :: gridptCount, gridpt, numLocalGridptsFoundSearch
    integer :: lonIndex, latIndex, resultsIndex, gridIndex
    type(kdtree2_result) :: searchResults(maxNumLocalGridptsSearch)
    real(kdkind)         :: refPosition(3), maxRadiusSquared
    real(4) :: footprintRadius_r4 ! (metres) used for seaice observations only
    real(4) :: influenceRadius_r4 ! (metres)
    real(8) :: obsLonInRad, obsLatInRad, maxLcorr

    call utl_tmg_start(122,'--AnalErrOI_FindObs')

    maxLcorr = maxval(Lcorr(:,:))
    numObs(:,:) = 0

    HEADER_LOOP: do headerIndex = 1, obs_numHeader(obsSpaceData)

      bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1

      BODY_LOOP: do bodyIndex = bodyIndexBeg, bodyIndexEnd

        if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY_LOOP

        if (trim(variableName) == 'GL') then
          footprintRadius_r4 = s2c_getFootprintRadius(obsSpaceData, stateVectorBkGnd, headerIndex)
          influenceRadius_r4 = max(0.0, footprintRadius_r4) + maxLcorr
        else if (trim(variableName) == 'TM') then 
          influenceRadius_r4 = maxLcorr
        else
          call utl_abort('findObs: The current code does not work with '&
                         //trim(variableName)//' analysis variable.')
        end if      

        if (maxLcorr == 0.0d0) then

          do kIndex = stateVectorBkGnd%mykBeg, stateVectorBkGnd%mykEnd
            do stepIndex = 1, stateVectorBkGnd%numStep
              do procIndex = 1, mmpi_nprocs

                call s2c_getWeightsAndGridPointIndexes(headerIndex, kIndex, stepIndex, &
                                                       procIndex, interpWeight, obsLatIndex, &
                                                       obsLonIndex, gridptCount)

                GRIDPT_LOOP: do gridpt = 1, gridptCount

                  if (interpWeight(gridpt) == 0.0d0) cycle GRIDPT_LOOP

                  lonIndex = obsLonIndex(gridpt)
                  latIndex = obsLatIndex(gridpt)

                  numObs(lonIndex, latIndex) = numObs(lonIndex, latIndex) + 1
                  if(associated(influentObs(lonIndex, latIndex)%bodyIndex)) then
                    if(numObs(lonIndex, latIndex) > influentObs(lonIndex, latIndex)%numObs) then
                      call utl_abort('findObs: Array too small')
                    end if
                    influentObs(lonIndex, latIndex)%headerIndex(numObs(lonIndex, latIndex)) = headerIndex
                    influentObs(lonIndex, latIndex)%bodyIndex(numObs(lonIndex, latIndex)) = bodyIndex
                  end if

                end do GRIDPT_LOOP

              end do
            end do
          end do

        else

          ! Determine the grid point nearest the observation.

          obsLonInRad = obs_headElem_r(obsSpaceData, obs_lon, headerIndex)
          obsLatInRad = obs_headElem_r(obsSpaceData, obs_lat, headerIndex)

          ! do the search
          maxRadiusSquared = real(influenceRadius_r4, 8) ** 2
          refPosition(:) = kdtree2_3dPosition(obsLonInRad, obsLatInRad)
          call kdtree2_r_nearest(tp = tree, qv = refPosition, r2 = maxRadiusSquared, &
                                 nfound = numLocalGridptsFoundSearch, &
                                 nalloc = maxNumLocalGridptsSearch, &
                                 results = searchResults)

          if (numLocalGridptsFoundSearch > maxNumLocalGridptsSearch) then
            call utl_abort('findObs: parameter maxNumLocalGridptsSearch must be increased')
          end if

          do resultsIndex = 1, numLocalGridptsFoundSearch

            gridIndex = searchResults(resultsIndex)%idx
            if (gridIndex < 1 .or. &
                gridIndex > stateVectorBkGnd%hco%ni * stateVectorBkGnd%hco%nj) then
              write(*,*) 'analysisErrorStdDev: gridIndex = ', gridIndex
              call utl_abort('findObs: gridIndex out of bound.')
            end if

            latIndex = (gridIndex - 1) / stateVectorBkGnd%hco%ni + 1
            lonIndex = gridIndex - (latIndex - 1) * stateVectorBkGnd%hco%ni
            if (lonIndex < 1 .or. lonIndex > stateVectorBkGnd%hco%ni .or. &
                latIndex < 1 .or. latIndex > stateVectorBkGnd%hco%nj) then
              write(*,*) 'analysisErrorStdDev: lonIndex = ', lonIndex, &
                                            ', latIndex = ', latIndex
              call utl_abort('findObs: lonIndex/latIndex out of bound.')
            end if

            numObs(lonIndex, latIndex) = numObs(lonIndex, latIndex) + 1
            if (associated(influentObs(lonIndex, latIndex)%bodyIndex)) then
              if (numObs(lonIndex, latIndex) > influentObs(lonIndex, latIndex)%numObs) then
                call utl_abort('findObs: Array too small in subroutine findObs')
              end if
              influentObs(lonIndex, latIndex)%headerIndex(numObs(lonIndex, latIndex)) = headerIndex
              influentObs(lonIndex, latIndex)%bodyIndex(numObs(lonIndex, latIndex)) = bodyIndex
            end if

          end do

        end if

      end do BODY_LOOP

    end do HEADER_LOOP

    call utl_tmg_stop(122)

  end subroutine findObs

  !---------------------------------------------------------
  ! aer_daysSinceLastObs
  !---------------------------------------------------------
  subroutine aer_daysSinceLastObs(obsSpaceData, hco_ptr, vco_ptr, &
                                  inputFileName, outputFileName, &
                                  variableName, propagateDSLO, &
                                  hoursSinceLastAnalysis)
    !
    !:Purpose: Update the field "days since last obs" with the newly assimilated obs.
    !
    implicit none

    ! Arguments
    type(struct_obs), intent(in) :: obsSpaceData
    type(struct_hco), pointer    :: hco_ptr
    type(struct_vco), pointer    :: vco_ptr
    character(len=*), intent(in) :: inputFileName  ! input file name
    character(len=*), intent(in) :: outputFileName ! output file name
    character(len=*), intent(in) :: variableName
    logical         , intent(in) :: propagateDSLO  ! propagate (increase) Days Since Last Obs in time
    integer         , intent(in) :: hoursSinceLastAnalysis

    ! Local Variables
    type(struct_gsv) :: stateVectorBkGnd
    type(struct_gsv) :: stateVectorAnal
    real(8), pointer :: bkGndDaysSinceLastObs_ptr(:,:,:,:)
    real(8), pointer :: analysisDaysSinceLastObs_ptr(:,:,:,:)
    integer :: stepIndex, levIndex, lonIndex, latIndex, headerIndex
    integer :: bodyIndexBeg, bodyIndexEnd, bodyIndex, kIndex, procIndex
    integer :: gridptCount, gridpt

    type(struct_columnData) :: column
    type(struct_columnData) :: columng

    real(8) :: leadTimeInHours

    if(mmpi_nprocs > 1) then
      write(*,*) 'mmpi_nprocs = ', mmpi_nprocs
      call utl_abort('aer_daysSinceLastObs: this version of the code should only be used with one mpi task.')
    end if
    if(mmpi_myid > 0) return

    write(*,*) '**********************************************************'
    write(*,*) '** aer_daysSinceLastObs: Update the days since last obs **'
    write(*,*) '**********************************************************'

    write(*,*) 'aer_daysSinceLastObs: input variable: ', trim(variableName)

    call gsv_allocate(stateVectorBkGnd, 1, hco_ptr, vco_ptr, dateStamp_opt = -1, &
                      mpi_local_opt = .true., mpi_distribution_opt = 'Tiles', &
                      varNames_opt=(/'DSLO'/), dataKind_opt = 8)
    call gsv_allocate(stateVectorAnal,  1, hco_ptr, vco_ptr, dateStamp_opt = -1, &
                      varNames_opt = (/'DSLO'/), dataKind_opt = 8)

    if (propagateDSLO) then
      call aer_propagateDSLO(stateVectorBkGnd, trim(variableName), &
                             inputFileName, outputFileName, &
                             hoursSinceLastAnalysis, hco_ptr, vco_ptr)
    else
      call gio_readFromFile(stateVectorBkGnd, inputFileName, '', 'P@')
    end if

    leadTimeInHours = real(stateVectorBkGnd%deet * stateVectorBkGnd%npasList(1), 8) / 3600.0d0
    call incdatr(stateVectorAnal%dateOriginList(1), stateVectorBkGnd%dateOriginList(1), &
                 leadTimeInHours)

    call gsv_copyMask(stateVectorBkGnd, stateVectorAnal)
    stateVectorAnal%etiket = stateVectorBkGnd%etiket

    call col_setVco(column, vco_ptr)
    call col_allocate(column,  obs_numHeader(obsSpaceData), varNames_opt=(/'DSLO'/))
    call col_setVco(columng, vco_ptr)
    call col_allocate(columng, obs_numHeader(obsSpaceData), varNames_opt=(/'DSLO'/))
    call s2c_tl(statevectorBkGnd, column, columng, obsSpaceData)

    ni = stateVectorBkGnd%hco%ni
    nj = stateVectorBkGnd%hco%nj

    call gsv_getField(stateVectorBkGnd,    bkGndDaysSinceLastObs_ptr, 'DSLO')
    call gsv_getField(stateVectorAnal , analysisDaysSinceLastObs_ptr, 'DSLO')

    ! Initialisation
    do stepIndex = 1, stateVectorBkGnd%numStep
      do levIndex = 1, gsv_getNumLev(stateVectorBkGnd,vnl_varLevelFromVarname('DSLO'))
        do latIndex = 1, nj
          do lonIndex = 1, ni
            analysisDaysSinceLastObs_ptr(lonIndex, latIndex, levIndex, stepIndex) = &
               bkGndDaysSinceLastObs_ptr(lonIndex, latIndex, levIndex, stepIndex)
          end do
        end do
      end do
    end do

    HEADER_LOOP: do headerIndex = 1, obs_numHeader(obsSpaceData)

      bodyIndexBeg = obs_headElem_i(obsSpaceData, obs_rln, headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData, obs_nlv, headerIndex) + bodyIndexBeg - 1

      BODY_LOOP: do bodyIndex = bodyIndexBeg, bodyIndexEnd

        if (obs_bodyElem_i(obsSpaceData, obs_ass, bodyIndex) /= obs_assimilated) then
          cycle BODY_LOOP
        end if

        do kIndex = stateVectorBkGnd%mykBeg, stateVectorBkGnd%mykEnd
          do stepIndex = 1, stateVectorBkGnd%numStep
            do procIndex = 1, mmpi_nprocs

              call s2c_getWeightsAndGridPointIndexes(headerIndex, kIndex, stepIndex, &
                                                     procIndex, interpWeight, obsLatIndex, &
                                                     obsLonIndex, gridptCount)

              GRIDPT_LOOP: do gridpt = 1, gridptCount

                if (interpWeight(gridpt) == 0.0d0) cycle GRIDPT_LOOP

                lonIndex = obsLonIndex(gridpt)
                latIndex = obsLatIndex(gridpt)

                do levIndex = 1, gsv_getNumLev(stateVectorBkGnd, vnl_varLevelFromVarname('DSLO'))
                  analysisDaysSinceLastObs_ptr(lonIndex, latIndex, levIndex, stepIndex) = 0.0
                end do

              end do GRIDPT_LOOP

            end do
          end do
        end do

      end do BODY_LOOP

    end do HEADER_LOOP

    call gio_writeToFile(stateVectorAnal, outputFileName, '', typvar_opt = 'A@', &
                         containsFullField_opt = .false. )

    call col_deallocate(columng)
    call col_deallocate(column)
    call gsv_deallocate(stateVectorBkGnd)
    call gsv_deallocate(stateVectorAnal)

  end subroutine aer_daysSinceLastObs

  !---------------------------------------------------------
  ! aer_propagateAnalysisError
  !---------------------------------------------------------
  subroutine aer_propagateAnalysisError(stateVectorErrorStd, oceanMask, variableName, &
                                        analysisEtiket, errorGrowth, hco_ptr, vco_ptr)
    !
    !:Purpose: read analysis error standard deviation field and propagate it forward in time
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: stateVectorErrorStd ! Input: analysis std error; Output: background std error. 
    type(struct_ocm), intent(in)    :: oceanMask           ! ocean-land mask (1=water, 0=land)
    character(len=*), intent(in)    :: variableName        ! variable name
    character(len=*), intent(in)    :: analysisEtiket      ! analysis etiket in the input std file 
    real(4)         , intent(in)    :: errorGrowth         ! seaice: fraction of ice per hour, SST: estimated growth
    type(struct_hco), pointer       :: hco_ptr             ! horizontal coordinates structure, pointer
    type(struct_vco), pointer       :: vco_ptr             ! vertical coordinates structure, pointer

    ! locals:
    integer :: latIndex, lonIndex, localLatIndex, localLonIndex 
    real(8), pointer :: stateVectorStdError_ptr(:,:,:)
    type(struct_gsv) :: stateVectorAnalysis
    real(8), pointer :: stateVectorAnalysis_ptr(:,:,:)
    real(4) :: totalLocalVariance
    integer :: pointCount

    write(*,*) ''
    write(*,*) 'aer_propagateAnalysisError: propagate analysis error forward in time for: ', &
               trim(variableName)
    
    ! read analysis itself (seaice concentration or SST analysis)
    call msg('aer_propagateAnalysisError:', ' reading analysis field...')
    call gsv_allocate(stateVectorAnalysis, 1, hco_ptr, vco_ptr, dateStamp_opt = -1, &
                      mpi_local_opt = .true., mpi_distribution_opt = 'Tiles', &
                      varNames_opt = (/trim(variableName)/), dataKind_opt = 8)
    call gio_readFromFile(stateVectorAnalysis, './anlm_000m', analysisEtiket, &
                          'A@', containsFullField_opt = .true.)
    call gsv_getField(stateVectorAnalysis, stateVectorAnalysis_ptr)

    ! initialize pointer for the error standard deviation field    
    call gsv_getField(stateVectorErrorStd, stateVectorStdError_ptr)
    
    ! calculation in variance unit
    stateVectorStdError_ptr(:,:,1) = stateVectorStdError_ptr(:,:,1)**2 

    pointCount = 0
    totalLocalVariance = 0.0d0

    do latIndex = 1, hco_ptr%nj
      do lonIndex = 1, hco_ptr%ni

        OCEANPOINTS: if (oceanMask%mask(lonIndex, latIndex, 1)) then
          do localLatIndex = latIndex - 1, latIndex + 1
            if (localLatIndex >= 1 .and. localLatIndex <= hco_ptr%nj) then
              do localLonIndex = lonIndex - 1,  lonIndex + 1
                if (localLonIndex >= 1 .and. localLonIndex <= hco_ptr%ni) then
                  if (oceanMask%mask(localLonIndex, localLatIndex, 1) .and. &
                      (localLonIndex /= lonIndex .or. localLatIndex /= latIndex)) then

                    pointCount = pointCount + 1
                    totalLocalVariance = totalLocalVariance + &
                                         (stateVectorAnalysis_ptr(localLonIndex, localLatIndex, 1) - &
                                          stateVectorAnalysis_ptr(lonIndex, latIndex, 1))**2

                  end if
                end if
              end do
            end if
          end do
        end if OCEANPOINTS

        if (pointCount > 0) then
          totalLocalVariance = totalLocalVariance / real(pointCount)
          if (stateVectorStdError_ptr(lonIndex, latIndex,1) < totalLocalVariance) then
            stateVectorStdError_ptr(lonIndex, latIndex, 1) = totalLocalVariance
          end if
        end if

      end do
    end do

    call gsv_deallocate(stateVectorAnalysis)

    ! Back to standard deviation
    stateVectorStdError_ptr(:, :, 1) = sqrt(stateVectorStdError_ptr(:, :, 1))

    ! Add the standard deviation increment
    do latIndex = 1, hco_ptr%nj
      do lonIndex = 1, hco_ptr%ni
        if (trim(variableName) == 'GL') then
          stateVectorStdError_ptr(lonIndex, latIndex, 1) = &
                     min(stateVectorStdError_ptr(lonIndex, latIndex, 1) + &
                         errorGrowth * tim_dstepobs, 1.0)
        else if (trim(variableName) == 'TM') then
          stateVectorStdError_ptr(lonIndex, latIndex, 1) = &
                         stateVectorStdError_ptr(lonIndex, latIndex, 1) + &
                         errorGrowth * tim_dstepobs / 2.0
        end if
      end do
    end do

  end subroutine aer_propagateAnalysisError

  !---------------------------------------------------------
  ! aer_propagateDSLO
  !---------------------------------------------------------
  subroutine aer_propagateDSLO(stateVectorErrorStd, variableName, &
                               inputFileName, outputFileName, &
                               hoursSinceLastAnalysis, hco_ptr, vco_ptr)
    !
    !:Purpose: propagate the field "days since last obs" in time.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: stateVectorErrorStd    ! read "analysed" state into it and increase it using 
                                                              ! hoursSinceLastAnalysis to make a background field 
    character(len=*), intent(in)    :: variableName           ! variable name ('GL' or 'TM') 
    character(len=*), intent(in)    :: inputFileName          ! input  file name
    character(len=*), intent(in)    :: outputFileName         ! output file name
    integer         , intent(in)    :: hoursSinceLastAnalysis ! hours since last analysis (namelist variable)
    type(struct_hco), pointer       :: hco_ptr
    type(struct_vco), pointer       :: vco_ptr

    ! locals
    real(8), pointer :: analysisDaysSinceLastObs_ptr(:,:,:)
    real(8) :: daysSinceLastAnalysis
    integer :: latIndex, lonIndex

    write(*,*) 'aer_propagateDSLO: propagating in time the Days Since Last Obs field...'
    write(*,*) 'aer_propagateDSLO: hours since last analysis: ', hoursSinceLastAnalysis

    daysSinceLastAnalysis = real(hoursSinceLastAnalysis, 8) / 24.0d0

    ! get DSLO field
    call gio_readFromFile(stateVectorErrorStd, inputFileName, '', 'A@')
    call gsv_getField(stateVectorErrorStd, analysisDaysSinceLastObs_ptr)

    ! Add number of days since last obs
    do latIndex = 1, hco_ptr%nj
      do lonIndex = 1, hco_ptr%ni
        analysisDaysSinceLastObs_ptr(lonIndex, latIndex, 1) = analysisDaysSinceLastObs_ptr(lonIndex, latIndex, 1) + &
                                                              daysSinceLastAnalysis
      end do
    end do

    ! update dateStamp from env variable
    call gsv_modifyDate(stateVectorErrorStd, tim_getDateStamp(), &
                        modifyDateOrigin_opt = .true.)

    ! save increased DSLO field into an fst-file
    call gio_writeToFile(stateVectorErrorStd, outputFileName, &
                         '', typvar_opt = 'P@', &
                         containsFullField_opt = .false.)

  end subroutine aer_propagateDSLO

end module analysisErrorOI_mod

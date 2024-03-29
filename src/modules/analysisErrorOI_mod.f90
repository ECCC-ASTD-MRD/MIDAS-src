
module analysisErrorOI_mod
  ! MODULE analysisErrorOI_mod (prefix='aer' category='1. High-level functionality')
  !
  !:Purpose:  Calculate the analysis-error standard deviation.
  !           The method used is Optimal Interpolation,
  !           where it is assumed that only a subset of the
  !           total number of observations influence the analysis at a given grid point.
  !

  use columnData_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use kdTree2_mod
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

  integer, external :: get_max_rss
  integer, parameter :: maxNumLocalGridptsSearch = 1500

contains

  !--------------------------------------------------------------------------
  ! aer_analysisError
  !--------------------------------------------------------------------------
  subroutine aer_analysisError(obsSpaceData, hco_ptr, vco_ptr)
    !
    ! :Purpose: Calculate analysis-error variance.
    !
    implicit none

    ! Arguments:
    type(struct_obs),          intent(inout) :: obsSpaceData ! observation data structure
    type(struct_hco), pointer, intent(in)    :: hco_ptr      ! horizontal grid definition
    type(struct_vco), pointer, intent(in)    :: vco_ptr      ! vertical grid definition

    ! Locals:
    integer :: fnom, fclos, nulnam, ierr
    type(struct_gsv) :: stateVectorAnlErrorStd      ! state vector for analysis error std deviation
    type(struct_gsv) :: stateVectorTrlErrorStd      ! state vector for background error std deviation
    real(8), pointer :: anlErrorStdDev_ptr(:,:,:,:) ! pointer for analysis error std deviation
    real(8), pointer :: trlErrorStdDev_ptr(:,:,:,:) ! pointer for background error std deviation
    integer, allocatable :: numObs(:,:)
    integer :: lonIndex, latIndex
    character(len=4), pointer :: analysisVariable(:)
    type(struct_gsv)          :: statevectorLcorr
    real(4), pointer          :: field3D_r4_ptr(:,:,:)
    type(struct_columnData)   :: column
    type(struct_columnData)   :: columng
    type(struct_ocm)          :: oceanMask
    real(8) :: leadTimeInHours
    character(len=20), parameter :: errStddevFileName_in  = 'errorstdev_in'        ! input  filename for anl ot trl error standard deviation
    character(len=20), parameter :: anlErrStddevFileName_out = 'anlerrorstdev_out' ! output filename for anl error std deviation
    character(len=20), parameter :: trlErrStddevFileName_out = 'trlerrorstdev_out' ! output filename for trl (background) error std deviation
    type(struct_neighborhood), pointer :: influentObs(:,:)
    real(8), allocatable :: Lcorr(:,:)
    
    ! namelist variables:
    real(8)           :: maxAnalysisErrorStdDev ! maximum limit imposed on analysis error stddev
    logical           :: propagateAnalysisError ! choose to propagate analysis error
    logical           :: propagateDSLO          ! choose to propagate Days Since Last Obs field
    real(4)           :: errorGrowth            ! seaice: fraction of ice per hour, SST: estimated growth
    character(len=12) :: analysisEtiket         ! analysis field etiket in a standard file
    character(len=12) :: anlErrorStdEtiket      ! analysis error standard deviation field etiket in the input/output standard files
    character(len=12) :: trlErrorStdEtiket      ! background error standard deviation field etiket in the input/output standard files
    integer           :: hoursSinceLastAnalysis ! number of hours since the last analysis
    logical           :: saveTrlStdField        ! choose to save trial standard deviation field
    character(len=2)  :: inputTypeVar           ! typvar of the analysis error field in the input file 
    character(len=2)  :: outputTypeVar          ! typvar of the analysis error field for the output file 
    real(4)           :: multFactorLcorr        ! multiplication scaling factor to increase the correlation length scale field
    namelist /namaer/ maxAnalysisErrorStdDev, propagateAnalysisError, propagateDSLO, &
                      errorGrowth, analysisEtiket, anlErrorStdEtiket, trlErrorStdEtiket, &
                      hoursSinceLastAnalysis, saveTrlStdField, inputTypeVar, outputTypeVar, &
                      multFactorLcorr

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
    anlErrorStdEtiket = 'A-ER STD DEV'
    trlErrorStdEtiket = 'B-ER STD DEV'
    hoursSinceLastAnalysis = 6
    saveTrlStdField = .false.
    inputTypeVar = 'P@'
    outputTypeVar = 'A@'
    multFactorLcorr = 1.0
    
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
               ' analysis error std field is read from: '//trim(errStddevFileName_in))    
      call gio_readFromFile(stateVectorAnlErrorStd, errStddevFileName_in, &
                            etiket_in = anlErrorStdEtiket, typvar_in = inputTypeVar, &
                            containsFullField_opt = .false.)
      
      ! initialize trl error std deviation field:
      trlErrorStdDev_ptr(:,:,:,:) = anlErrorStdDev_ptr(:,:,:,:)

      call ocm_readMaskFromFile (oceanMask, hco_ptr, vco_ptr, errStddevFileName_in)
      call aer_propagateAnalysisError (stateVectorTrlErrorStd, oceanMask, &
                                       analysisVariable(1), &
                                       analysisEtiket, errorGrowth, &
                                       hco_ptr, vco_ptr)

      ! impose maximum value on trial error standard deviation field
      trlErrorStdDev_ptr(:,:,:,:) = min(trlErrorStdDev_ptr(:,:,:,:), &
                                        maxAnalysisErrorStdDev)

      if (saveTrlStdField) then
        ! zap analysis error etiket with background error etiket
        stateVectorTrlErrorStd%etiket = trlErrorStdEtiket
        
        ! copy mask from analysis error std deviation field to trl error std field
        call gsv_copyMask(stateVectorAnlErrorStd, stateVectorTrlErrorStd)
        
        ! update dateStamp from env variable
        call gsv_modifyDate(stateVectorTrlErrorStd, tim_getDateStamp(), &
                            modifyDateOrigin_opt = .true.)

        ! save background error (increased analysis error) standard deviation field
        call gio_writeToFile(stateVectorTrlErrorStd, trlErrStddevFileName_out, &
                             trlErrorStdEtiket, typvar_opt = outputTypeVar, &
                             containsFullField_opt = .false.)
      end if

    else
      call msg('aer_analysisError:', &
               ' trial error std field is read from: '//trim(errStddevFileName_in))
      call gio_readFromFile(stateVectorTrlErrorStd, errStddevFileName_in, &
                            etiket_in = trlErrorStdEtiket, typvar_in = inputTypeVar, &
                            containsFullField_opt = .false.)
      call gsv_copyMask(stateVectorTrlErrorStd, stateVectorAnlErrorStd)
    end if

    leadTimeInHours = real(stateVectorTrlErrorStd%deet * stateVectorTrlErrorStd%npasList(1), 8) / 3600.0d0
    call incdatr(stateVectorAnlErrorStd%dateOriginList(1), &
                 stateVectorTrlErrorStd%dateOriginList(1), leadTimeInHours)

    stateVectorAnlErrorStd%etiket = anlErrorStdEtiket

    call col_setVco(column, vco_ptr)
    call col_allocate(column,  obs_numHeader(obsSpaceData), varNames_opt=(/analysisVariable(1)/))
    call col_setVco(columng, vco_ptr)
    call col_allocate(columng, obs_numHeader(obsSpaceData), varNames_opt=(/analysisVariable(1)/))
    call s2c_tl(stateVectorTrlErrorStd, column, columng, obsSpaceData)

    allocate(Lcorr(hco_ptr%ni, hco_ptr%nj))

    ! get correlation length scale field
    write(*,*) 'aer_analysisError: get correlation length scale field...'
    call gsv_allocate(statevectorLcorr, 1, hco_ptr, vco_ptr, dateStamp_opt = -1, &
                      dataKind_opt = 4, hInterpolateDegree_opt = 'LINEAR', &
                      varNames_opt = (/analysisVariable(1)/))
    call gsv_zero(statevectorLcorr)
    call gio_readFromFile(statevectorLcorr, './bgstddev', 'CORRLEN', ' ', &
                          unitConversion_opt = .false.)
    call gsv_getField(statevectorLcorr, field3D_r4_ptr, analysisVariable(1))

    ! apply multiplication scaling factor
    field3D_r4_ptr(:, :, 1) = field3D_r4_ptr(:, :, 1) * multFactorLcorr

    ! Convert from km to meters
    Lcorr(:,:) = 1000.0d0 * real(field3D_r4_ptr(:, :, 1), 8)

    write(*,*) 'aer_analysisError: min/max correlation length scale 2D field: ', &
               minval(Lcorr(:,:) ), maxval(Lcorr(:,:))
    call gsv_deallocate(statevectorLcorr)

    write(*,*) 'aer_analysisError: done creating kdtree for stateVectorTrlErrorStd'
    write(*,*) 'Memory Used: ', get_max_rss() / 1024, 'Mb'

    ! Go through all observations a first time to get
    ! the number of influential observations
    ! in order to allocate memory appropriately.

    allocate(influentObs(hco_ptr%ni, hco_ptr%nj))
    allocate(numObs(hco_ptr%ni, hco_ptr%nj))

    call msg('aer_analysisError:', ' looking for observations...')
    call findObs(obsSpaceData, stateVectorTrlErrorStd, numObs, &
                 analysisVariable(1), Lcorr, influentObs)

    ! Memory allocation

    do latIndex = 1, hco_ptr%nj
      do lonIndex = 1, hco_ptr%ni
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
    call findObs(obsSpaceData, stateVectorTrlErrorStd, numObs, &
                 analysisVariable(1), Lcorr, influentObs)

    deallocate(numObs)

    write(*,*) 'aer_analysisError:: obs found'
    write(*,*) 'Memory Used: ',get_max_rss() / 1024, 'Mb'

    ! compute analysis error standard deviation
    call aer_computeAnlErrorStd(obsSpaceData, stateVectorAnlErrorStd, stateVectorTrlErrorStd, &
                                analysisVariable(1), maxAnalysisErrorStdDev, influentObs, Lcorr)

    deallocate(influentObs)
    deallocate(Lcorr)

    ! update dateStamp from env variable
    call gsv_modifyDate(stateVectorAnlErrorStd, tim_getDateStamp(), &
                        modifyDateOrigin_opt = .true.)

    ! save analysis error
    call msg('aer_analysisError:', ' writing analysis error std field to output file...')
    call gio_writeToFile(stateVectorAnlErrorStd, anlErrStddevFileName_out, &
                         anlErrorStdEtiket, typvar_opt = outputTypeVar, &
                         containsFullField_opt = .false.)

    call col_deallocate(columng)
    call col_deallocate(column)
    call gsv_deallocate(stateVectorAnlErrorStd)
    call gsv_deallocate(stateVectorTrlErrorStd)

    if (analysisVariable(1) == 'GL') then
      ! Update the Days Since Last Obs
      call aer_daysSinceLastObs(obsSpaceData, hco_ptr, vco_ptr, &
                                errStddevFileName_in, anlErrStddevFileName_out, &
                                analysisVariable(1), propagateDSLO, &
                                hoursSinceLastAnalysis)
    end if
    
    call msg('aer_analysisError:', ' finished.')

  end subroutine aer_analysisError

  !---------------------------------------------------------
  ! findObs
  !---------------------------------------------------------
  subroutine findObs(obsSpaceData, stateVectorTrlErrorStd, numObs, variableName, Lcorr, influentObs)
    !
    !:Purpose: Find all observations used for the analysis.
    !
    implicit none

    ! Arguments:
    type(struct_obs),                   intent(in)    :: obsSpaceData     ! observation data structure
    type(struct_gsv),                   intent(in)    :: stateVectorTrlErrorStd ! state containing background error stddev
    integer,                            intent(out)   :: numObs(:,:)      ! number of observations found
    character(len=*),                   intent(in)    :: variableName     ! 'GL' for seaice or 'TM' for SST
    real(8)         ,                   intent(in)    :: Lcorr(:,:)       ! horizontal background-error correlation length scale
    type(struct_neighborhood), pointer, intent(inout) :: influentObs(:,:) ! details about observations to use in update

    ! Locals:
    integer :: headerIndex, bodyIndexBeg, bodyIndexEnd, bodyIndex
    integer :: procIndex, kIndex, stepIndex
    integer :: gridptCount, gridpt, numLocalGridptsFoundSearch
    integer :: lonIndex, latIndex, resultsIndex, gridIndex
    type(kdtree2_result) :: searchResults(maxNumLocalGridptsSearch)
    real(kdkind)         :: refPosition(3), maxRadiusSquared
    real(4) :: footprintRadius_r4 ! (metres) used for seaice observations only
    real(4) :: influenceRadius_r4 ! (metres)
    real(8) :: obsLonInRad, obsLatInRad, maxLcorr
    type(kdtree2), pointer :: tree
    real(kdkind), allocatable :: positionArray(:,:)
    real(8),      allocatable :: latInRad(:,:), lonInRad(:,:)
    real(8) :: interpWeight(maxNumLocalGridptsSearch)
    integer :: obsLatIndex(maxNumLocalGridptsSearch), obsLonIndex(maxNumLocalGridptsSearch)

    call utl_tmg_start(122,'--AnalErrOI_FindObs')

    ! create kdtree
    write(*,*) 'findObs: start creating kdtree for stateVectorTrlErrorStd'
    write(*,*) 'Memory Used: ', get_max_rss() / 1024, 'Mb'

    allocate(positionArray(3, stateVectorTrlErrorStd%hco%ni * stateVectorTrlErrorStd%hco%nj))
    allocate(lonInRad(stateVectorTrlErrorStd%hco%ni, stateVectorTrlErrorStd%hco%nj))
    allocate(latInRad(stateVectorTrlErrorStd%hco%ni, stateVectorTrlErrorStd%hco%nj))

    gridIndex = 0
    do latIndex = 1, stateVectorTrlErrorStd%hco%nj
      do lonIndex = 1, stateVectorTrlErrorStd%hco%ni
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
    deallocate(lonInRad)
    deallocate(latInRad)

    maxLcorr = maxval(Lcorr(:,:))
    numObs(:,:) = 0

    HEADER_LOOP: do headerIndex = 1, obs_numHeader(obsSpaceData)

      bodyIndexBeg = obs_headElem_i(obsSpaceData, obs_rln, headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData, obs_nlv, headerIndex) + bodyIndexBeg - 1

      BODY_LOOP: do bodyIndex = bodyIndexBeg, bodyIndexEnd

        if (obs_bodyElem_i(obsSpaceData, obs_ass, bodyIndex) /= obs_assimilated) cycle BODY_LOOP

        if (trim(variableName) == 'GL') then
          footprintRadius_r4 = s2c_getFootprintRadius(obsSpaceData, stateVectorTrlErrorStd, headerIndex)
          influenceRadius_r4 = max(0.0, footprintRadius_r4) + maxLcorr
        else if (trim(variableName) == 'TM') then 
          influenceRadius_r4 = maxLcorr
        else
          call utl_abort('findObs: The current code does not work with '&
                         //trim(variableName)//' analysis variable.')
        end if      

        if (maxLcorr == 0.0d0) then

          do kIndex = stateVectorTrlErrorStd%mykBeg, stateVectorTrlErrorStd%mykEnd
            do stepIndex = 1, stateVectorTrlErrorStd%numStep
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
                gridIndex > stateVectorTrlErrorStd%hco%ni * stateVectorTrlErrorStd%hco%nj) then
              write(*,*) 'analysisErrorStdDev: gridIndex = ', gridIndex
              call utl_abort('findObs: gridIndex out of bound.')
            end if

            latIndex = (gridIndex - 1) / stateVectorTrlErrorStd%hco%ni + 1
            lonIndex = gridIndex - (latIndex - 1) * stateVectorTrlErrorStd%hco%ni
            if (lonIndex < 1 .or. lonIndex > stateVectorTrlErrorStd%hco%ni .or. &
                latIndex < 1 .or. latIndex > stateVectorTrlErrorStd%hco%nj) then
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

    ! Arguments:
    type(struct_obs),          intent(inout) :: obsSpaceData   ! observation data structure
    type(struct_hco), pointer, intent(in)    :: hco_ptr        ! horizontal grid definition
    type(struct_vco), pointer, intent(in)    :: vco_ptr        ! vertical grid definition
    character(len=*),          intent(in)    :: inputFileName  ! input file name
    character(len=*),          intent(in)    :: outputFileName ! output file name
    character(len=*),          intent(in)    :: variableName   ! name of variable being treated
    logical         ,          intent(in)    :: propagateDSLO  ! propagate (increase) Days Since Last Obs in time
    integer         ,          intent(in)    :: hoursSinceLastAnalysis ! number of hours between analysis times

    ! Locals:
    type(struct_gsv) :: stateVectorTrlDSLO, stateVectorAnlDSLO
    real(8), pointer :: trlDSLO_ptr(:,:,:,:), anlDSLO_ptr(:,:,:,:)
    integer :: stepIndex, levIndex, lonIndex, latIndex, headerIndex
    integer :: bodyIndexBeg, bodyIndexEnd, bodyIndex, kIndex, procIndex
    integer :: gridptCount, gridpt
    type(struct_columnData) :: column, columng
    real(8) :: leadTimeInHours, interpWeight(maxNumLocalGridptsSearch)
    integer :: obsLatIndex(maxNumLocalGridptsSearch), obsLonIndex(maxNumLocalGridptsSearch)

    if(mmpi_nprocs > 1) then
      write(*,*) 'mmpi_nprocs = ', mmpi_nprocs
      call utl_abort('aer_daysSinceLastObs: this version of the code should only be used with one mpi task.')
    end if
    if(mmpi_myid > 0) return

    write(*,*) '**********************************************************'
    write(*,*) '** aer_daysSinceLastObs: Update the days since last obs **'
    write(*,*) '**********************************************************'

    write(*,*) 'aer_daysSinceLastObs: input variable: ', trim(variableName)

    call gsv_allocate(stateVectorTrlDSLO, 1, hco_ptr, vco_ptr, dateStamp_opt = -1, &
                      mpi_local_opt = .true., mpi_distribution_opt = 'Tiles', &
                      varNames_opt=(/'DSLO'/), dataKind_opt = 8)
    call gsv_allocate(stateVectorAnlDSLO,  1, hco_ptr, vco_ptr, dateStamp_opt = -1, &
                      varNames_opt = (/'DSLO'/), dataKind_opt = 8)

    if (propagateDSLO) then
      call aer_propagateDSLO(stateVectorTrlDSLO, inputFileName, outputFileName, &
                             hoursSinceLastAnalysis, hco_ptr)
    else
      write(*,*) 'aer_daysSinceLastObs: DSLO trial field is read from: ', trim(inputFileName)
      call gio_readFromFile(stateVectorTrlDSLO, inputFileName, &
                            etiket_in = ' ', typvar_in = 'P@')
    end if

    leadTimeInHours = real(stateVectorTrlDSLO%deet * stateVectorTrlDSLO%npasList(1), 8) / 3600.0d0
    call incdatr(stateVectorAnlDSLO%dateOriginList(1), stateVectorTrlDSLO%dateOriginList(1), &
                 leadTimeInHours)

    call gsv_copyMask(stateVectorTrlDSLO, stateVectorAnlDSLO)
    stateVectorAnlDSLO%etiket = stateVectorTrlDSLO%etiket

    call col_setVco(column, vco_ptr)
    call col_allocate(column,  obs_numHeader(obsSpaceData), varNames_opt=(/'DSLO'/))
    call col_setVco(columng, vco_ptr)
    call col_allocate(columng, obs_numHeader(obsSpaceData), varNames_opt=(/'DSLO'/))
    call s2c_tl(stateVectorTrlDSLO, column, columng, obsSpaceData)

    call gsv_getField(stateVectorTrlDSLO, trlDSLO_ptr, 'DSLO')
    call gsv_getField(stateVectorAnlDSLO, anlDSLO_ptr, 'DSLO')

    ! Initialisation
    do stepIndex = 1, stateVectorTrlDSLO%numStep
      do levIndex = 1, gsv_getNumLev(stateVectorTrlDSLO, vnl_varLevelFromVarname('DSLO'))
        do latIndex = 1, hco_ptr%nj
          do lonIndex = 1, hco_ptr%ni
            anlDSLO_ptr(lonIndex, latIndex, levIndex, stepIndex) = &
               trlDSLO_ptr(lonIndex, latIndex, levIndex, stepIndex)
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

        do kIndex = stateVectorTrlDSLO%mykBeg, stateVectorTrlDSLO%mykEnd
          do stepIndex = 1, stateVectorTrlDSLO%numStep
            do procIndex = 1, mmpi_nprocs

              call s2c_getWeightsAndGridPointIndexes(headerIndex, kIndex, stepIndex, &
                                                     procIndex, interpWeight, obsLatIndex, &
                                                     obsLonIndex, gridptCount)

              GRIDPT_LOOP: do gridpt = 1, gridptCount

                if (interpWeight(gridpt) == 0.0d0) cycle GRIDPT_LOOP

                lonIndex = obsLonIndex(gridpt)
                latIndex = obsLatIndex(gridpt)

                do levIndex = 1, gsv_getNumLev(stateVectorTrlDSLO, vnl_varLevelFromVarname('DSLO'))
                  anlDSLO_ptr(lonIndex, latIndex, levIndex, stepIndex) = 0.0
                end do

              end do GRIDPT_LOOP

            end do
          end do
        end do

      end do BODY_LOOP

    end do HEADER_LOOP

    call gio_writeToFile(stateVectorAnlDSLO, outputFileName, '', typvar_opt = 'A@', &
                         containsFullField_opt = .false. )

    call col_deallocate(columng)
    call col_deallocate(column)
    call gsv_deallocate(stateVectorTrlDSLO)
    call gsv_deallocate(stateVectorAnlDSLO)

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

    ! Arguments:
    type(struct_gsv),          intent(inout) :: stateVectorErrorStd ! Input: analysis std error; Output: background std error. 
    type(struct_ocm),          intent(in)    :: oceanMask           ! ocean-land mask (1=water, 0=land)
    character(len=*),          intent(in)    :: variableName        ! variable name
    character(len=*),          intent(in)    :: analysisEtiket      ! analysis etiket in the input std file 
    real(4)         ,          intent(in)    :: errorGrowth         ! seaice: fraction of ice per hour, SST: estimated growth
    type(struct_hco), pointer, intent(in)    :: hco_ptr             ! horizontal coordinates structure, pointer
    type(struct_vco), pointer, intent(in)    :: vco_ptr             ! vertical coordinates structure, pointer

    ! Locals:
    type(struct_gsv) :: stateVectorAnalysis
    integer :: latIndex, lonIndex, localLatIndex, localLonIndex, pointCount 
    real(8), pointer :: stateVectorStdError_ptr(:,:,:), stateVectorAnalysis_ptr(:,:,:)
    real(4) :: totalLocalVariance

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
  subroutine aer_propagateDSLO(stateVectorErrorStd, inputFileName, outputFileName, &
                               hoursSinceLastAnalysis, hco_ptr)
    !
    !:Purpose: propagate the field "days since last obs" in time.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),          intent(inout) :: stateVectorErrorStd    ! read "analysis" into it and increase by hoursSinceLastAnalysis
    character(len=*),          intent(in)    :: inputFileName          ! input file name
    character(len=*),          intent(in)    :: outputFileName         ! output file name
    integer         ,          intent(in)    :: hoursSinceLastAnalysis ! hours since last analysis (namelist variable)
    type(struct_hco), pointer, intent(in)    :: hco_ptr                ! horizontal grid definition

    ! Locals:
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

  !---------------------------------------------------------
  ! aer_computeAnlErrorStd
  !---------------------------------------------------------
  subroutine aer_computeAnlErrorStd(obsSpaceData, stateVectorAnlErrorStd, stateVectorTrlErrorStd, &
                                    analysisVariable, maxAnalysisErrorStdDev, influentObs, Lcorr)
    !
    !:Purpose: compute analysis error standard deviation for
    !          one analysis variable (grid point) at a time
    !
    implicit none

    ! Arguments:
    type(struct_obs)         ,          intent(in)    :: obsSpaceData           ! obsSpaceData structure
    type(struct_gsv)         ,          intent(inout) :: stateVectorAnlErrorStd ! state vector for analysis error std deviation
    type(struct_gsv)         ,          intent(in)    :: stateVectorTrlErrorStd ! state vector for background error std deviation
    character(len=*)         ,          intent(in)    :: analysisVariable       ! variable name ('GL' or 'TM') 
    real(8)                  ,          intent(in)    :: maxAnalysisErrorStdDev ! maximum limit imposed on analysis error stddev
    type(struct_neighborhood), pointer, intent(in)    :: influentObs(:,:)       ! details about observations to use in update
    real(8)                  ,          intent(in)    :: Lcorr(:,:)             ! horizontal background-error length scale

    ! Locals:
    integer :: latIndex, lonIndex, stepIndex, levIndex, kIndex
    integer :: numInfluentObs, bodyIndex, numVariables
    integer :: influentObsIndex2, influentObsIndex
    integer :: xStateIndex, yStateIndex, procIndex
    integer :: xIndex1, yIndex1, xIndex2, yIndex2
    integer :: gridptCount, gridpt, headerIndex
    integer :: varIndex1, varIndex2, currentAnalVarIndex, ni, nj
    real(8), allocatable :: latInRad(:,:), lonInRad(:,:)
    real(8), pointer     :: anlErrorStdDev_ptr(:,:,:,:)  ! pointer for analysis error std deviation
    real(8), pointer     :: trlErrorStdDev_ptr(:,:,:,:)  ! pointer for background error std deviation
    real(8), allocatable :: obsOperator(:,:), Bmatrix(:,:), PHiA(:)
    real(8), allocatable :: innovCovariance(:,:), obsErrorVariance(:)
    real(8), allocatable :: PH(:,:), KH(:), IKH(:), innovCovarianceInverse(:,:)
    integer, parameter :: maxvar = 15000
    integer :: statei(maxvar), statej(maxvar) ! Model grid coordinate i,j of neighbors
    real(8) :: scaling, distance
    logical :: found
    real(8) :: interpWeight(maxNumLocalGridptsSearch)
    integer :: obsLatIndex(maxNumLocalGridptsSearch), obsLonIndex(maxNumLocalGridptsSearch)

    call msg('aer_computeAnlErrorStd:', ' computing analysis error std field initialization...')

    ni = stateVectorTrlErrorStd%hco%ni
    nj = stateVectorTrlErrorStd%hco%nj

    ! compute lon/lat in radians
    allocate(lonInRad(ni, nj))
    allocate(latInRad(ni, nj))

    do latIndex = 1, nj
      do lonIndex = 1, ni
        latInRad(lonIndex, latIndex) = real(stateVectorTrlErrorStd%hco%lat2d_4(lonIndex, latIndex), 8)
        lonInRad(lonIndex, latIndex) = real(stateVectorTrlErrorStd%hco%lon2d_4(lonIndex, latIndex), 8)
      end do
    end do
    
    ! Initialisation of pointers
    call gsv_getField(stateVectorAnlErrorStd, anlErrorStdDev_ptr)
    call gsv_getField(stateVectorTrlErrorStd, trlErrorStdDev_ptr)
    anlErrorStdDev_ptr(:,:,:,:) = trlErrorStdDev_ptr(:,:,:,:)
    call msg('aer_computeAnlErrorStd:', ' analysis error std field initialization...OK')

    ! Only variables assigned within or by the loop can be private.
    STEP: do stepIndex = 1, stateVectorTrlErrorStd%numStep
      LEVEL: do levIndex = 1, gsv_getNumLev(stateVectorTrlErrorStd, &
                                            vnl_varLevelFromVarname(analysisVariable))

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

              if (analysisVariable == 'GL') then
                scaling = oop_iceScaling(obsSpaceData, bodyIndex)
              else if (analysisVariable == 'TM') then
                scaling = 1.0d0
              end if	

              if (scaling == 0.0d0) cycle INFLUENTOBSCYCLE
              KINDEXCYCLE: do kIndex = stateVectorTrlErrorStd%mykBeg, stateVectorTrlErrorStd%mykEnd
                PROCINDEXCYCLE: do procIndex = 1, mmpi_nprocs

                  call s2c_getWeightsAndGridPointIndexes(headerIndex, kIndex, stepIndex, procIndex, &
                                                         interpWeight, obsLatIndex, obsLonIndex, &
                                                         gridptCount)

                  GRIDPTCYCLE: do gridpt = 1, gridptCount

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

              if (analysisVariable == 'GL') then
                scaling = oop_iceScaling(obsSpaceData, bodyIndex)
              else if (analysisVariable == 'TM') then
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

    deallocate(lonInRad)
    deallocate(latInRad)

  end subroutine aer_computeAnlErrorStd

end module analysisErrorOI_mod

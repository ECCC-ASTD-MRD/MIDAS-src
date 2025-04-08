
module oceanBackground_mod
  ! MODULE oceanBackground_mod (prefix='obgd' category='1. High-level functionality')
  !
  !:Purpose: storage for ocean background related subroutines
  !
  use codtyp_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use utilities_mod
  
  implicit none
  save
  private

  ! Public functions/subroutines
  public :: obgd_computeSSTrial
 
  ! External functions
  integer, external :: fnom, fclos  

  contains

  !----------------------------------------------------------------------------------------
  ! obgd_computeSSTrial
  !----------------------------------------------------------------------------------------
  subroutine obgd_computeSSTrial(hco, vco, trialDateStamp, analysisDateStamp, &
                                 nmonthsClim, datestampClim, alphaClim, etiket_in)
    !
    !: Purpose: 1) to compute SST analysis anomaly w.r.t. climatology
    !              x_anomaly(t-1) = (xa(t-1) - xclim(t-1))
    !           2) to compute SST background   
    !              xb(t) = x_anomaly(t-1) * alpha + xclim(t)             
    implicit none

    ! Arguments:
    type(struct_hco), pointer, intent(in) :: hco               ! horizontal grid
    type(struct_vco), pointer, intent(in) :: vco               ! vertical grid
    integer                  , intent(in) :: trialDateStamp    ! trial datestamp
    integer                  , intent(in) :: analysisDateStamp ! datestamp for last analysis 
    integer                  , intent(in) :: nmonthsClim       ! number of climatological fields (= 12)
    integer                  , intent(in) :: datestampClim(:)  ! datestamps of input climatology fields
    real(8)                  , intent(in) :: alphaClim         ! scalling factor to relax towards climatology
    character(len=10)        , intent(in) :: etiket_in         ! etiket from namelist and for trial
    
    ! Locals:
    type(struct_gsv) :: stateVector, stateVectorAnomaly 
    real(8), pointer :: stateVector_ptr(:, :, :), stateVectorAnomaly_ptr(:, :, :)
    integer          :: lonIndex, latIndex, status
    real(8)          :: climatology_m1(hco % ni, hco % nj), climatology(hco % ni, hco % nj)
    character(len=12) :: etiket
    
    write(*,*) 'obgd_computeSSTrial: starting...'  
      
    ! make a copy of analysis file
    status = utl_copyFile('./analysis', './analysisAndAnomaly')

    ! get SST analysis
    call gsv_allocate(stateVector, 1, hco, vco, dataKind_opt = 8, &
                      dateStampList_opt = (/analysisDateStamp/), &
                      mpi_local_opt = .true., varNames_opt = (/'TM'/))
    call gio_readFromFile(stateVector, './analysisAndAnomaly', ' ','A', &
                          unitConversion_opt=.false., &
                          containsFullField_opt=.true.)
    call gsv_getField(stateVector, stateVector_ptr)

    ! initialize state vector for analysis anomaly
    call gsv_allocate(stateVectorAnomaly, 1, hco, vco, dataKind_opt = 8, &
                      dateStampList_opt = (/analysisDateStamp/), &
                      mpi_local_opt = .true., varNames_opt = (/'TM'/))
    call gsv_copy(stateVector, stateVectorAnomaly)
    call gsv_getField(stateVectorAnomaly, stateVectorAnomaly_ptr)
    
    call obgd_getClimatology(analysisDateStamp, hco, vco, nmonthsClim, datestampClim, climatology_m1)
    call obgd_getClimatology(trialDateStamp   , hco, vco, nmonthsClim, datestampClim, climatology)
    
    ! compute background state
    ! xb(t) = (xa(t-1) - xclim(t-1))*alpha + xclim(t)
    do lonIndex = stateVector%myLonBeg, stateVector%myLonEnd 
      do latIndex = stateVector%myLatBeg, stateVector%myLatEnd
        stateVectorAnomaly_ptr(lonIndex, latIndex, 1) = stateVector_ptr(lonIndex, latIndex, 1) - &
                                                        climatology_m1(lonIndex, latIndex)
        stateVector_ptr(lonIndex, latIndex, 1) = stateVectorAnomaly_ptr(lonIndex, latIndex, 1) * &
                                                 alphaClim + climatology(lonIndex, latIndex)
      end do
    end do

    ! save analysis anomaly in RPN standard file
    etiket = 'ANOMALY'
    call gio_writeToFile(stateVectorAnomaly, './analysisAndAnomaly', etiket, typvar_opt = 'A@')

    ! modify dateStamp (from analysis) with trial dateStamp
    call gsv_modifyDate(stateVector, trialDateStamp, modifyDateOrigin_opt = .true.)
    
    ! save trial field in RPN standard file
    call gio_writeToFile(stateVector, './trial', etiket_in, typvar_opt = 'P@')

    ! save climatology corresponding to the analysisDateStamps
    do lonIndex = stateVector%myLonBeg, stateVector%myLonEnd
      do latIndex = stateVector%myLatBeg, stateVector%myLatEnd
        stateVector_ptr(lonIndex, latIndex, 1) = climatology_m1(lonIndex, latIndex)
      end do
    end do
    etiket = 'CLIMATO'
    call gio_writeToFile(stateVector, './analysisAndAnomaly', etiket, typvar_opt = 'C@')

    call gsv_deallocate(stateVector)
    call gsv_deallocate(stateVectorAnomaly)

  end subroutine obgd_computeSSTrial
  
  !----------------------------------------------------------------------------------------
  ! obgd_getClimatology
  !----------------------------------------------------------------------------------------
  subroutine obgd_getClimatology(dateStamp, hco, vco, nmonthsClim, datestampClim, output)
    !
    !: Purpose: 1) to read SST climatological fields from a std file
    !           2) to interpolate the field in time using the current day (t) in current month (m)    
    !           SST(t) = SST_neighbourMonth * weight + SST_currentMonth * (1.0 - weight),
    !
    implicit none

    ! Arguments:
    integer                  , intent(in)    :: dateStamp        ! date stamp for the current day
    type(struct_hco), pointer, intent(in)    :: hco              ! horizontal grid
    type(struct_vco), pointer, intent(in)    :: vco              ! vertical grid
    integer                  , intent(in)    :: nmonthsClim      ! number of records in the climatology file 
    integer                  , intent(in)    :: datestampClim(:) ! datestamps in the climatology file
    real(8)                  , intent(inout) :: output(:,:)      ! interpolated SST field from climatology
  
    ! Locals:
    integer          :: hour, day, month, yyyy  ! these variables correspond to the current time (dateStamp)
    integer          :: ndays                   ! number of days in the current month
    integer          :: neighbourMonth          ! neighbour month to the current month, can be previous or next
    integer          :: neighbourYear
    type(struct_gsv) :: stateVector, stateVector_neighbourMonth
    real(8), pointer :: clim_ptr(:, :, :), clim_neighbourMonth_ptr(:, :, :)
    integer          :: lonIndex, latIndex
    real(8)          :: weight                  ! weight for linear interpolation of climatology in time
    real(8)          :: numberHours             ! number of hours between the 15th of the current and neighbour months
    integer          :: dataStampMonth          ! dataStamp of the 15th of the current month
    integer          :: dataStampNeighbourMonth ! dataStamp of the 15th of the neighbour month

    call tim_dateStampToYYYYMMDDHH(dateStamp, hour, day, month, ndays, yyyy)
    write(*,'(a,3i5,a,i12,a)') 'obgd_getClimatology: interpolating climatology for day/month/year (datestamp): ', &
    day, month, yyyy, '(', datestamp, ')'

    neighbourYear = yyyy

    if (day < 15) then
      if (month == 1) then
        neighbourMonth = 12
        neighbourYear  = yyyy - 1
      else
        neighbourMonth = month - 1
      end if
    else ! day >=15
      if (month == nmonthsClim) then
        neighbourMonth = 1
        neighbourYear  = yyyy + 1
      else
        neighbourMonth = month + 1
      end if
    end if

    ! compute datestamp of the 15th of the current month/year
    dataStampMonth = tim_yyyymmddhhToDatestamp(yyyy, month, 15, hour)
    ! compute datestamps of the 15th of the neighbour month/year
    dataStampNeighbourMonth = tim_yyyymmddhhToDatestamp(neighbourYear, neighbourMonth, 15, hour)
    write(*,*) 'obgd_getClimatology: the 15th of the current   month dateStamp: ', dataStampMonth
    write(*,*) 'obgd_getClimatology: the 15th of the neighbour month dateStamp: ', dataStampNeighbourMonth

    ! Difference (in hours) between the 15th of the current and neighbour months
    call difdatr(dataStampMonth, dataStampNeighbourMonth, numberHours)
    if (numberHours == 2.d0**30) then
      call utl_abort('obgd_getClimatology: difdatr received invalid arguments: '//&
                                           char(dataStampMonth)//' and '//&
                                           char(dataStampNeighbourMonth))
    end if 
    ! safe check: the number of hours cannot be greater than 31 days * 24h = 744h
    if (abs(numberHours) > 744.d0) then
      call utl_abort('obgd_getClimatology: number of hours between two months exceeds 744h!')
    end if
    ! safe check: the number of hours cannot be equal to zero
    if (numberHours == 0.d0) then
      call utl_abort('obgd_getClimatology: number of hours between two neighbour months is zero!')
    end if
 
    ! computing weight for linear interpolation of climatology in time
    ! The weight depends on the distance between the current day and the 15th of the month
    ! and the number of days between the 15th of the month and the neighbour month  
    weight = abs(real(day - 15, 8) * 24.d0 / numberHours)
    write(*,*) 'obgd_getClimatology: weight for the current date: ', weight

    ! get climatology, current month
    write(*,*) 'obgd_getClimatology: reading climatology, current month: ', &
               month, ', datestamp: ', datestampClim(month) 
    call gsv_allocate(stateVector, 1, hco, vco, dataKind_opt = 8, &
                      dateStampList_opt = datestampClim(month:month), mpi_local_opt = .true., &
                      varNames_opt = (/'TM'/), hInterpolateDegree_opt ='LINEAR')
    call gio_readFromFile(stateVector, './climatology', ' ',' ', &
                          unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVector, clim_ptr)
    
    ! get climatology, neighbour month
    write(*,*) 'obgd_getClimatology: reading climatology, neighbour month: ', &
               neighbourMonth, ', datestamp: ', datestampClim(neighbourMonth) 
    call gsv_allocate(stateVector_neighbourMonth, 1, hco, vco, dataKind_opt = 8, &
                      dateStampList_opt = datestampClim(neighbourMonth:neighbourMonth), mpi_local_opt = .true., &
                      varNames_opt = (/'TM'/), hInterpolateDegree_opt ='LINEAR')
    call gio_readFromFile(stateVector_neighbourMonth, './climatology', ' ',' ', &
                          unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVector_neighbourMonth, clim_neighbourMonth_ptr)

    do lonIndex = stateVector%myLonBeg, stateVector%myLonEnd 
      do latIndex = stateVector%myLatBeg, stateVector%myLatEnd
        output(lonIndex, latIndex) = clim_ptr(lonIndex, latIndex, 1) * (1.0d0 - weight) + &
                                     clim_neighbourMonth_ptr(lonIndex, latIndex, 1) * weight
      end do
    end do
    
    call gsv_deallocate(stateVector)
    call gsv_deallocate(stateVector_neighbourMonth)
  
  end subroutine obgd_getClimatology 
  
end module oceanBackground_mod  

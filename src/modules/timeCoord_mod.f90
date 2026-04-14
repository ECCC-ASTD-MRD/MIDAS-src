
module timeCoord_mod
  ! MODULE timeCoord_mod (prefix='tim' category='7. Low-level data objects')
  !
  !:Purpose:  To store public variables and procedures related to the time
  !           coordinate.
  !
  use rmn_fst98
  use rmn_date
  use midasMpi_mod
  use varNameList_mod
  use utilities_mod
  use runtimeInfo_mod
  use message_mod
  use utilities_mod
  use runtimeInfo_mod

  implicit none
  save
  private

  ! Public module variables
  real(8), public, protected :: tim_dstepobs, tim_dstepobsinc, tim_windowsize
  integer, public, protected :: tim_nstepobs, tim_nstepobsinc
  logical, public, protected :: tim_fullyUseExtremeTimeBins
  character(len=6), public, protected :: tim_referenceTime

  ! Public module procedures
  public :: tim_setup, tim_initialized
  public :: tim_getDateStamp, tim_setDateStamp, tim_getStampList, tim_getStepObsIndex
  public :: tim_getDateStampFromFile, tim_dateStampToYYYYMMDDHH, tim_getValidDateTimeFromList
  public :: tim_yyyymmddhhToDatestamp, tim_getHoursBetweenDates, tim_getSecondsBetweenDates

  ! Private variables
  character(len=4) :: varNameForDate
  integer :: datestamp = 0  ! datestamp is usually the centre of time window
  logical :: initialized = .false.

  ! module interfaces
  ! -----------------

  ! general interface for tim_dateStampToYYYYMMDDHH
  interface tim_dateStampToYYYYMMDDHH
    module procedure tim_dateStampToYYYYMMDDHHWithDaysInMonth
    module procedure tim_dateStampToYYYYMMDDHHOnly
    module procedure tim_dateStampToYYYYMMDDHHPrintable
  end interface tim_dateStampToYYYYMMDDHH

  ! general interface for tim_yyyymmddhhToDatestamp
  interface tim_yyyymmddhhToDatestamp
    module procedure tim_yyyymmddhhToDatestampWithyyyymmddhh
    module procedure tim_yyyymmddhhToDatestampWithDateTime
  end interface tim_yyyymmddhhToDatestamp

contains


  !----------------------------------------------------------------------------------------
  ! tim_setup
  !----------------------------------------------------------------------------------------
  subroutine tim_setup(fileNameForDate_opt)
    !
    ! :Purpose: Setup of obs time window size and related trial field
    !           time step for OmP determination.
    !
    implicit none

    ! Arguments:
    character(len=*), optional, intent(in) :: fileNameForDate_opt ! Optional Input file name to the the datetimestamp from

    ! Locals:
    integer :: prntdate, prnttime
    integer :: dateStampEnvVar

    call tim_readNml()

    if (initialized) then
      write(*,*) 'tim_setup: already initialized, just return'
      return
    end if

    ! First try to set dateStamp from MIDAS_DATE
    dateStampEnvVar = tim_getDateStampFromEnvVar()

    ! Possibly set the datestamp (except when set later from burp files)
    if (dateStampEnvVar /= 0) then
      write(*,*) 'tim_setup: ====================================================='
      write(*,*) 'tim_setup: DATESTAMP set by value in supplied MIDAS_DATE'
      write(*,*) 'tim_setup: ====================================================='
      dateStamp = dateStampEnvVar
      call tim_dateStampToYYYYMMDDHH(datestamp, prntdate, prnttime)
      write(*,*) 'tim_setup: printdate = ', prntdate
      write(*,*) 'tim_setup: printtime = ', prnttime
      write(*,*) 'tim_setup: datestamp = ', datestamp
    else if (present(fileNameForDate_opt)) then
      write(*,*) 'tim_setup: ====================================================='
      write(*,*) 'tim_setup: DATESTAMP set by value in supplied file'
      write(*,*) 'tim_setup: ====================================================='
      datestamp = tim_getDatestampFromFile(fileNameForDate_opt)
      call tim_dateStampToYYYYMMDDHH(datestamp, prntdate, prnttime)
      write(*,*) 'tim_setup: printdate = ', prntdate
      write(*,*) 'tim_setup: printtime = ', prnttime
      write(*,*) 'tim_setup: datestamp = ', datestamp
    else
      write(*,*) 'tim_setup: =========================================================='
      write(*,*) 'tim_setup: DATESTAMP not set in this subroutine, use tim_setDateStamp'
      write(*,*) 'tim_setup: =========================================================='
    end if

    if (mmpi_myid == 0) write(*,*) 'tim_setup: dobs_windowsize   = ', tim_windowsize
    if (mmpi_myid == 0) write(*,*) 'tim_setup: dstepobs          = ', tim_dstepobs
    if (mmpi_myid == 0) write(*,*) 'tim_setup: nstepobs          = ', tim_nstepobs
    if (mmpi_myid == 0) write(*,*) 'tim_setup: dstepobsinc       = ', tim_dstepobsinc
    if (mmpi_myid == 0) write(*,*) 'tim_setup: nstepobsinc       = ', tim_nstepobsinc
    if (mmpi_myid == 0) write(*,*) 'tim_setup: tim_referenceTime = ', tim_referenceTime

    if (tim_nstepobs == 0 .or. tim_nstepobsinc == 0) then
      call rti_abort('tim_setup: Wrong configuration, nstepobs/nstepobsinc can not be 0.')
    end if

    initialized = .true.

  end subroutine tim_setup


  !----------------------------------------------------------------------------------------
  ! tim_initialized
  !----------------------------------------------------------------------------------------
  function tim_initialized() result(initialized_out)
    implicit none

    ! Result:
    logical initialized_out ! Return if the module has been initialized or not

    initialized_out = initialized

  end function tim_initialized


  !----------------------------------------------------------------------------------------
  ! tim_readNml (private)
  !----------------------------------------------------------------------------------------
  subroutine tim_readNml()
    !
    ! :Purpose: Read the namelist block NAMTIME.
    !
    ! :Namelist parameters:
    !         :dstepobs:    time step (hrs) between successive trial fields
    !                       for use in OmP determation. Set to dwindowsize for
    !                       single trial field, i.e. use of 3dvar instead of 3dvar-FGAT.
    !                       nstepobs = number of trial fields
    !         :dstepobsinc: time step (hrs) between obs groupings in time. Set to
    !                       dwindowsize for use of a single obs group.
    !                       nstepobsinc = number of obs time intervals
    !         :dwindowsize: Time window size (hrs).
    !
    ! :Comment:
    !
    !   Provided dates and number of provided trial field files must be
    !   consistent with nstepobs, dstepobs and dwindowsize with reference datestamp
    !   corresponding to the date of the middle trial field file.
    !
    implicit none

    ! Locals:
    integer :: ierr
    logical, save :: firstCall = .true.
    integer :: nulnam

    ! Namelist variables:
    real(8) :: dstepobs      ! time step length for background state (in hours)
    real(8) :: dstepobsinc   ! time step length for increment and/or B matrix (in hours)
    real(8) :: dwindowsize   ! length of assimilation window (in hours)
    character(len=6) :: referenceTime  ! location of 'date' within the window: 'middle' or 'start'
    logical :: fullyUseExtremeTimeBins ! choose to use full-size bins at both ends of window (usually only half size)

    NAMELIST /NAMTIME/dstepobs, dstepobsinc, dwindowsize, referenceTime, fullyUseExtremeTimeBins

    if (.not.firstCall) then
      write(*,*) 'tim_readNml: already initialized, just return'
      return
    else
      firstCall = .false.
    end if

    ! Set default values for namelist variables
    dstepobs       = 6.0d0
    dstepobsinc    = 6.0d0
    dwindowsize    = 6.0d0
    referenceTime = 'middle'
    fullyUseExtremeTimeBins = .false.

    ! Read the namelist
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml = namtime, iostat = ierr)
    if (ierr /= 0) call rti_abort('tim_readNml: Error reading namelist')
    if (mmpi_myid == 0) write(*, nml = namtime)
    ierr = fclos(nulnam)

    ! Set the module variables for timestep length, number of timesteps and window length
    tim_dstepobs      = dstepobs
    tim_dstepobsinc   = dstepobsinc
    tim_windowsize    = dwindowsize
    tim_referenceTime = referenceTime
    tim_fullyUseExtremeTimeBins = fullyUseExtremeTimeBins

    if (tim_fullyUseExtremeTimeBins .and. tim_referencetime=='middle') then
      write(*,*) 'Warning: tim_fullyUseExtremeTimeBins = .true. and tim_referencetime = "middle" is a non-standard combination'
      write(*,*) 'Is it really what you want ?'
    end if

    if (dstepobs > dwindowsize) then
      if (mmpi_myid == 0) write(*,*) 'tim_readNml: dstepobs>dwindowsize. Reset to dwindowsize value.'
      tim_dstepobs = tim_windowsize
    end if
    if (dstepobsinc > dwindowsize) then
      if (mmpi_myid == 0) write(*,*) 'tim_readNml: dstepobsinc>dwindowsize. Reset to dwindowsize value.'
      tim_dstepobsinc = tim_windowsize
    end if

    if (tim_referenceTime == 'middle' .or. trim(tim_referenceTime) == 'end') then

      tim_nstepobs    = 2 * nint(((tim_windowsize - tim_dstepobs) / 2.d0) / tim_dstepobs) + 1
      tim_nstepobsinc = 2 * nint(((tim_windowsize - tim_dstepobsinc) / 2.d0) / tim_dstepobsinc) + 1

    else if (trim(tim_referenceTime) == 'start') then

      tim_nstepobs = max(nint(tim_windowsize / tim_dstepobs), 1)
      tim_nstepobsinc = max(nint(tim_windowsize / tim_dstepobsinc), 1)

    end if

  end subroutine tim_readNml


  !----------------------------------------------------------------------------------------
  ! tim_getDateStampFromEnvVar (private)
  !----------------------------------------------------------------------------------------
  function tim_getDateStampFromEnvVar() result(dateStamp)
    !
    !:Purpose: Determine the date from the environment variable 'MIDAS_DATE'.
    !
    implicit none

    ! Result:
    integer :: dateStamp ! datetime stamp extracted from environment variable 'MIDAS_DATE'

    ! Locals:
    integer    :: lengthValidDateStr, status
    integer(8) :: dateTimePrint, datePrint, timePrint
    character(len=256) :: validDateStr

    status = 0
    call get_environment_variable('MIDAS_DATE', validDateStr, lengthValidDateStr, status, .true.)

    if (status > 1) then
      call rti_abort('tim_getDateStampFromEnvVar: Problem when getting the environment variable MIDAS_DATE')
    end if

    if (status == 1) then
      write(*,*) 'tim_getDateStampFromEnvVar: WARNING: The environment variable MIDAS_DATE has not been detected!'
      dateStamp = 0
      return
    end if

    write(*,*)
    write(*,*) 'tim_getDateStampFromEnvVar: The environment variable MIDAS_DATE has correctly been detected'

    ! convert string to long integer
    read (validDateStr,*) dateTimePrint

    ! split dateTime long integer into separate date and time values
    if (lengthValidDateStr == 10) then
      datePrint = dateTimePrint/100
      timePrint = (dateTimePrint - datePrint*100) * 1000000
    else if (lengthValidDateStr == 12) then
      datePrint = dateTimePrint/10000
      timePrint = (dateTimePrint - datePrint*10000) * 10000
    else if (lengthValidDateStr == 14) then
      datePrint = dateTimePrint/1000000
      timePrint = (dateTimePrint - datePrint*1000000) * 100
    else if (lengthValidDateStr == 16) then
      datePrint = dateTimePrint/100000000
      timePrint = (dateTimePrint - datePrint*100000000)
    else
      write(*,*) 'length of MIDAS_DATE = ', lengthValidDateStr
      call rti_abort('tim_getDateStampFromEnvVar: Unexpected length of variable MIDAS_DATE')
    end if

    ! convert to CMC dateStamp
    datestamp = tim_yyyymmddhhToDatestamp(int(datePrint,4), int(timePrint,4))

    write(*,*) 'tim_getDateStampFromEnvVar: envVar, validDate, dateStamp = ', trim(validDateStr), dateTimePrint, dateStamp

  end function tim_getDateStampFromEnvVar


  !----------------------------------------------------------------------------------------
  ! tim_getDatestampFromFile
  !----------------------------------------------------------------------------------------
  function tim_getDatestampFromFile(fileName, varNameForDate_opt) result(dateStamp_out)
    !
    ! :Purpose: to extract the dateStamp from the supplied file.
    !
    implicit none

    ! Arguments:
    character(len=*),           intent(in) :: fileName           ! File used to estimate the valid date
    character(len=*), optional, intent(in) :: varNameForDate_opt ! Optional 'nomvar' to use to find the valid date
    ! Result:
    integer :: dateStamp_out

    ! Locals:
    integer :: nulFile, ierr
    integer, parameter :: maxNumDates = 2000
    integer :: numDates, ikeys(maxNumDates), varIndex
    integer :: prntdate, prnttime, windowIndex, windowsPerDay, dateStamp_tmp
    logical :: fileExists, foundWindow, foundVarNameInFile
    real(8) :: leadTimeInHours, windowBegHour, windowEndHour, fileHour, middleHour
    integer :: ideet, inpas, dateStamp_origin, ini, inj, ink, inbits, idatyp
    integer :: ip1, ip2, ip3, ig1, ig2, ig3, ig4, iswa, ilng, idltf, iubc
    integer :: iextra1, iextra2, iextra3
    character(len=2)  :: typvar
    character(len=4)  :: nomvar
    character(len=12) :: etiket
    character(len=1)  :: grtyp

    if (mmpi_myid == 0) then

      ! Check if file for any date within the analysis window (except the last) exists
      inquire(file=trim(fileName), exist=fileExists)
      if (.not.fileExists) then
        call rti_abort('tim_getDateStampFromFile: file not found '//trim(fileName))
      end if

      ! Determine variable to use for the date (default is P0)
      varNameForDate = 'P0'
      if (present(varNameForDate_opt)) then

        varNameForDate = trim(varNameForDate_opt)
        write(*,*) 'tim_getDateStampFromFile: defining dateStamp from the variable = ', varNameForDate

      ! If P0 not present, look for another suitable variable in the file
      else if (.not. vnl_varNamePresentInFile(varNameForDate, fileName = trim(fileName))) then

        foundVarNameInFile = .false.
        do varIndex = 1, vnl_numvarmax
          varNameForDate = vnl_varNameList(varIndex)
          ! check if variable is in the file
          if (.not. vnl_varNamePresentInFile(varNameForDate, fileName = trim(fileName))) cycle
          foundVarNameInFile = .true.
          exit
        end do

        if (.not. foundVarNameInFile) then
          call rti_abort('tim_getDateStampFromFile: NO variables found in the file!!!')
        end if

      end if

      ! Extract the datestamp from the file
      nulFile = 0
      ierr = fnom(nulFile, trim(fileName), 'RND+OLD+R/O', 0)
      ierr = fstouv(nulFile,'RND+OLD')
      ierr = fstinl(nulFile,ini,inj,ink,-1,' ',-1,-1,-1,' ', &
                    trim(varNameForDate),ikeys,numdates,maxnumdates)
      if (ikeys(1) <= 0) then
        call rti_abort('tim_getDateStampFromFile: Could not find variable ' //  &
                       trim(varNameForDate) // ' in the supplied file')
      end if
      write(*,*) 'tim_getDateStampFromFile: number of dates found = ', numDates
      ierr = fstprm(ikeys(1), dateStamp_origin, ideet, inpas, ini, inj, &
                    ink, inbits, idatyp, ip1, ip2, ip3, &
                    typvar, nomvar, etiket, grtyp, ig1, ig2, ig3, ig4, &
                    iswa, ilng, idltf, iubc, iextra1, iextra2, iextra3)
      leadTimeInHours = real(ideet*inpas,8)/3600.0d0
      call incdatr(dateStamp_out, dateStamp_origin, leadTimeInHours)

      ierr = fstfrm(nulFile)
      ierr = fclos(nulFile)

    end if

    call mmpi_bcast(dateStamp_out)

    if (tim_referenceTime == 'middle') then
      ! Modify date to ensure that it corresponds to the middle of the window
      ! Note: For this, we have to assume that the date in the file
      !       does NOT correspond to the final time of the window
      call tim_dateStampToYYYYMMDDHH(datestamp_out, prntdate, prnttime)
      fileHour = real(prnttime,8)/1000000.0d0
      windowsPerDay = nint(24.0d0 / tim_windowsize)
      foundWindow = .false.
      window_loop: do windowIndex = 0, windowsPerDay
        windowBegHour = (real(windowIndex,8) * tim_windowsize) - (tim_windowsize/2.0d0)
        windowEndHour = (real(windowIndex,8) * tim_windowsize) + (tim_windowsize/2.0d0)
        if (fileHour >= windowBegHour .and. fileHour < windowEndHour) then
          foundWindow = .true.
          middleHour = real(windowIndex,8) * tim_windowsize
          exit window_loop
        end if
      end do window_loop

      if (.not. foundWindow) then
        write(*,*) 'windowsPerDay, fileHour = ', windowsPerDay, fileHour
        call rti_abort('tim_getDateStampFromFile: could not determine assimilation window position')
      end if

      ! handle special case when window centered on hour 24
      if (nint(middleHour) == 24) then
        ! add 24h to dateStamp and recompute prntdate
        dateStamp_tmp = dateStamp_out
        call incdatr(dateStamp_out, dateStamp_tmp, 24.0d0)
        call tim_dateStampToYYYYMMDDHH(datestamp_out, prntdate, prnttime)

        ! subtract 24h from middleHour
        middleHour = 0.0d0
      end if

      prnttime = nint(middleHour) * 1000000
      datestamp_out = tim_yyyymmddhhToDatestamp(prntdate, prnttime)
    end if

  end function tim_getDateStampFromFile


  !----------------------------------------------------------------------------------------
  ! tim_setDatestamp
  !----------------------------------------------------------------------------------------
  subroutine tim_setDatestamp(datestamp_in)
    !
    ! :Purpose: to control access to the minimization object.  Sets the date
    !           of the window centre of analysis validity to the indicated value.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: datestamp_in ! set the internal 'datestamp' value in 'timeCoord_mod' module

    if (.not.initialized) call rti_abort('tim_setDateStamp: module not initialized')

    datestamp = datestamp_in

  end subroutine tim_setDatestamp


  !----------------------------------------------------------------------------------------
  ! tim_getDatestamp
  !----------------------------------------------------------------------------------------
  function tim_getDatestamp() result(datestamp_out)
    !
    ! :Purpose: to control access to the minimization object.  Returns the date
    !           of the window centre of analysis validity.
    !
    implicit none

    ! Result:
    integer :: datestamp_out ! get the internal 'datestamp' value in 'timeCoord_mod' module

    if (.not.initialized) call rti_abort('tim_getDateStamp: module not initialized')

    datestamp_out = datestamp

  end function tim_getDatestamp


  !----------------------------------------------------------------------------------------
  ! tim_getStampList
  !----------------------------------------------------------------------------------------
  subroutine tim_getStampList(dateStampList, numStep, referenceDateStamp)
    !
    ! :Purpose: Compute a list of STAMPS corresponding to stepobs time
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: numStep                ! number of step obs
    integer, intent(in)  :: referenceDateStamp     ! Synoptic time
    integer, intent(out) :: dateStampList(numStep) ! datestamp list

    ! Locals:
    integer :: stepIndex
    integer :: prntdate, prnttime
    real(8) :: dldelt ! delta time in hours between middle time and each step
    real(8) :: dtstep ! delta time in hours between step obs

    if (.not. initialized) call rti_abort('tim_getStampList: module not initialized')

    if (referenceDateStamp == -1) then

      if (mmpi_myid == 0) write(*,*) 'tim_getStampList: datestamp is not specified, keep as -1'
      dateStampList(:) = -1

    else

      if (tim_referenceTime == 'middle') then

        if (numStep > 1) then
          dtstep = tim_windowsize / (real(numStep - 1, 8))
        else
          dtstep = tim_windowsize
        end if

        do stepIndex = 1, numStep
          dldelt = (stepIndex - ((numStep - 1) / 2 + 1)) * dtstep
          call incdatr(dateStampList(stepIndex), referenceDateStamp, dldelt)
        end do

      else if (trim(tim_referencetime) == 'start') then

        dtstep = tim_windowsize / (real(numStep, 8))

        do stepIndex = 1, numStep
          dldelt = (stepIndex - 1) * dtstep
          call incdatr(dateStampList(stepIndex), referenceDateStamp, dldelt)
        end do

      else if (trim(tim_referencetime) == 'end') then

        if (numStep > 1) then

          call incdatr(dateStampList(1), referenceDateStamp, -tim_windowsize)

          dtstep = tim_windowsize / (real(numStep - 1, 8))

          do stepIndex = 2, numStep
            call incdatr(dateStampList(stepIndex), dateStampList(stepIndex - 1), dtstep)
          end do

        else
          call incdatr(dateStampList(1), referenceDateStamp, -tim_windowsize/2)
        end if

      end if

    end if ! datestamp is specified or not

    call msg('tim_getStampList', 'datestamp list of '//str(numStep)//' (numStep) states:')
    do stepIndex = 1, numStep
      call tim_dateStampToYYYYMMDDHH(dateStampList(stepIndex), prntdate, prnttime)
      write(*,*) stepIndex, dateStampList(stepIndex), prntdate, prnttime / 1000000, 'h'
    end do
    call msg('tim_getStampList', 'Completed')

  end subroutine tim_getStampList


  !----------------------------------------------------------------------------------------
  ! tim_getStepObsIndex
  !----------------------------------------------------------------------------------------
  function tim_getStepObsIndex(referenceDateStamp, obsYYYMMDD, obsHHMM, numStep) result(dnstepobs)
    !
    ! :Purpose: Return the stepobs index as a real number (-1.0 if out of range)
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: referenceDateStamp ! Synop CMC date-time stamp
    integer, intent(in)  :: obsYYYMMDD         ! Obs date YYYYMMDD
    integer, intent(in)  :: obsHHMM            ! Obs time HHMM
    integer, intent(in)  :: numStep            ! number of stepobs in assimilation window

    ! Result:
    real(8) :: dnstepobs ! number of stepobs from reference time

    ! Locals:
    real(8) :: dddt    ! delta time in hours
    integer :: istobs  ! obs CMC date-time stamp
    real(8) :: dlhours ! delta time from synop time

    if (.not. initialized) call rti_abort('tim_getStepObsIndex: module not initialized')

    ! Building observation stamp
    istobs = tim_yyyymmddhhToDatestamp(obsYYYMMDD, obsHHMM*10000)

    ! Difference (in hours) between obs time and reference time
    call difdatr(istobs, referenceDateStamp, dlhours)

    if (numStep > 1) then

      ! FGAT: more than 1 trial field in time window
      if (tim_referenceTime == 'middle') then
        dddt = tim_windowsize / (real(numStep - 1, 8))
      else if (tim_referenceTime == 'start') then
        dddt = tim_windowsize / (real(numStep, 8))
      else if (tim_referenceTime == 'end') then
        dddt = tim_windowsize / (real(numStep - 1, 8))
      end if

      dnstepobs = dlhours / dddt ! number of step obs from reference (e.g. synoptic)

      if (tim_referenceTime == 'middle') then
        dnstepobs = dnstepobs + real((numStep + 1) / 2, 8)
      else if (trim(tim_referencetime) == 'start') then
        dnstepobs = dnstepobs + 1.d0
      else if (trim(tim_referencetime) == 'end') then
        dnstepobs = dnstepobs + real(numStep, 8)
      end if

      if (dnstepobs < 0.5d0 .or. dnstepobs > (0.5d0 + real(numStep, 8))) dnstepobs = -1.0d0

    else
      ! 3D: only 1 trial field in time window
      if (tim_referenceTime == 'middle') then
        if (dlhours < -tim_windowsize / 2.0D0 .or. dlhours > tim_windowsize / 2.0D0) then
          ! outside time window
          dnstepobs = -1.0d0
        else
          ! inside time window
          dnstepobs = 1.0d0
        end if
      else if (tim_referenceTime == 'start') then
        dddt = tim_windowsize
        if (dlhours < -dddt / 2.0d0 .or. dlhours > tim_windowsize + dddt / 2.d0) then
          ! outside time window
          dnstepobs = -1.0d0
        else
          ! inside time window
          dnstepobs = 1.0d0
        end if
      else if (tim_referenceTime == 'end') then
        if (dlhours < -tim_windowsize .or. dlhours > 0.0D0) then
          ! outside time window
          dnstepobs = -1.0d0
        else
          ! inside time window
          dnstepobs = 1.0d0
        end if
      endif
    end if

  end function tim_getStepObsIndex

  !----------------------------------------------------------------------------------------
  ! tim_dateStampToYYYYMMDDHHPrintable
  !----------------------------------------------------------------------------------------
  subroutine tim_dateStampToYYYYMMDDHHPrintable(dateStamp, printableDate, printableTime, verbose_opt)
    !
    ! :Purpose: to get day (DD), month (MM), number of days in this month
    !           and year (YYYY) from dateStamp
    !
    implicit none

    ! Arguments:
    integer,           intent(in)  :: dateStamp     ! datetimestamp in CMC format
    integer,           intent(out) :: printableDate ! date extracted from 'datestamp' in format 'YYYYMMDD'
    integer,           intent(out) :: printableTime ! time extracted from 'datestamp' in format 'HHMMSS'
    logical, optional, intent(in)  :: verbose_opt   ! allows to write the output information

    ! Locals:
    character(len=8) :: yyyymmdd
    integer          :: mode, ierr, prntdate(2)
    logical          :: verbose

    if (present(verbose_opt)) then
      verbose = verbose_opt
    else
      verbose = .true.
    end if

    mode = -3 ! stamp to printable
    ierr = newdate(dateStamp, prntdate, printableTime, mode)
    if ( ierr < 0 ) then
      call rti_abort('tim_dateStampToYYYYMMDDHHPrintable: Invalid datestamp when calling ''newdate'' with mode=-3: ' // trim(utl_str(datestamp)))
    end if

    printableDate = prntdate(1)
    write(yyyymmdd,'(i8)') printableDate

    if(verbose) then
      write(*,*) 'tim_dateStampToYYYYMMDDHH: date = ', printableDate
      write(*,*) 'tim_dateStampToYYYYMMDDHH: time = ', printableTime
    end if

  end subroutine tim_dateStampToYYYYMMDDHHPrintable

  !----------------------------------------------------------------------------------------
  ! tim_dateStampToYYYYMMDDHHOnly
  !----------------------------------------------------------------------------------------
  subroutine tim_dateStampToYYYYMMDDHHOnly(dateStamp, time, dd, mm, yyyy, verbose_opt)
    !
    ! :Purpose: to get day (DD), month (MM), number of days in this month
    !           and year (YYYY) from dateStamp
    !
    implicit none

    ! Arguments:
    integer,           intent(in)  :: dateStamp   ! datetimestamp in CMC format
    integer,           intent(out) :: time        ! time extracted from 'datestamp' in format 'HHMM'
    integer,           intent(out) :: dd          ! day of month extracted from 'datestamp'
    integer,           intent(out) :: mm          ! month extracted from 'datestamp'
    integer,           intent(out) :: yyyy        ! year extracted from 'datestamp'
    logical, optional, intent(in)  :: verbose_opt ! allows to write the output information

    ! Locals:
    character(len=8) :: yyyymmdd
    integer          :: date
    logical          :: verbose

    if (present(verbose_opt)) then
      verbose = verbose_opt
    else
      verbose = .true.
    end if

    call tim_dateStampToYYYYMMDDHHPrintable(dateStamp, date, time, verbose_opt = .false.)
    write(yyyymmdd,'(i8)') date
    read (yyyymmdd(1:4), '(i4)') yyyy
    read (yyyymmdd(5:6), '(i2)') mm
    read (yyyymmdd(7:8), '(i2)') dd

    if(verbose) then
      write(*,*) 'tim_dateStampToYYYYMMDDHH: date  = ', yyyymmdd
      write(*,*) 'tim_dateStampToYYYYMMDDHH: year  = ', yyyy
      write(*,*) 'tim_dateStampToYYYYMMDDHH: month = ', mm
      write(*,*) 'tim_dateStampToYYYYMMDDHH: day   = ', dd
      write(*,*) 'tim_dateStampToYYYYMMDDHH: time  = ', time
    end if

  end subroutine tim_dateStampToYYYYMMDDHHOnly

  !----------------------------------------------------------------------------------------
  ! tim_dateStampToYYYYMMDDHHWithDaysInMonth
  !----------------------------------------------------------------------------------------
  subroutine tim_dateStampToYYYYMMDDHHWithDaysInMonth(dateStamp, prnttime, dd, mm, ndays, yyyy, verbose_opt)
    !
    ! :Purpose: to get day (DD), month (MM), number of days in this month
    !           and year (YYYY) from dateStamp
    !
    implicit none

    ! Arguments:
    integer,           intent(in)  :: dateStamp   ! datetimestamp in CMC format
    integer,           intent(out) :: prnttime    ! time extracted from 'datestamp' in format 'HHMM'
    integer,           intent(out) :: dd          ! day of month extracted from 'datestamp'
    integer,           intent(out) :: mm          ! month extracted from 'datestamp'
    integer,           intent(out) :: ndays       ! number of days in that month
    integer,           intent(out) :: yyyy        ! year extracted from 'datestamp'
    logical, optional, intent(in)  :: verbose_opt ! allows to write the output information

    ! Constants:
    character(len=3), parameter :: months(12) = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    ! Locals:
    logical :: verbose

    if (present(verbose_opt)) then
      verbose = verbose_opt
    else
      verbose = .true.
    end if

    call tim_dateStampToYYYYMMDDHHOnly(dateStamp, prnttime, dd, mm, yyyy, verbose_opt = .false.)

    ndays = tim_numberOfDaysInMonth(dateStamp,yyyy,mm,dd,prnttime)

    if(verbose) then
      write(*,*) 'tim_dateStampToYYYYMMDDHH: date = ', dateStamp
      write(*,*) 'tim_dateStampToYYYYMMDDHH: year = ', yyyy
      write(*,'(a,i5,a,i5,a)') ' tim_dateStampToYYYYMMDDHH: month = ', mm, ' ( '// months(mm)//' where there are ', ndays, ' days)'
      write(*,*) 'tim_dateStampToYYYYMMDDHH: day = ', dd
      write(*,*) 'tim_dateStampToYYYYMMDDHH: time = ', prnttime
    end if

  end subroutine tim_dateStampToYYYYMMDDHHWithDaysInMonth

  !----------------------------------------------------------------------------------------
  ! tim_yyyymmddhhToDatestampWithDateTime
  !----------------------------------------------------------------------------------------
  function tim_yyyymmddhhToDatestampWithDateTime(date, time) result(currentDateStamp)
    !
    !:Purpose: compute datestamp from year, month day and hour.
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: date ! date as an integer in format yyyymmdd
    integer, intent(in)  :: time ! time as an integer in format HHMMSShh

    ! Result:
    integer :: currentDateStamp ! datetimestamp in CMC format computed from 'date' and 'time'

    ! Locals:
    integer :: mode, ierr
    integer :: date_array(2)

    date_array(:) = date

    mode = 3
    ierr = newdate(currentDateStamp, date_array, time, mode)
    if ( ierr < 0 ) then
      call rti_abort('tim_yyyymmddhhToDatestampWithDateTime: Invalid datestamp when calling ''newdate'' with mode=3: ' // trim(utl_str(datestamp)))
    end if

  end function tim_yyyymmddhhToDatestampWithDateTime

  !----------------------------------------------------------------------------------------
  ! tim_yyyymmddhhToDatestampWithyyyymmddhh
  !----------------------------------------------------------------------------------------
  function tim_yyyymmddhhToDatestampWithyyyymmddhh(year, month, day, time) result(currentDateStamp)
    !
    !:Purpose: compute datestamp from year, month day and hour.
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: year  ! year in format yyyy to compute dateStamp
    integer, intent(in)  :: month ! month number in [1,12]
    integer, intent(in)  :: day   ! day number to compute dateStamp
    integer, intent(in)  :: time  ! time as an integer in format HHMMSShh

    ! Result:
    integer :: currentDateStamp

    ! Locals:
    integer :: date

    date = tim_dateFromYYYYMMDD(year, month, day)

    currentDateStamp = tim_yyyymmddhhToDatestampWithDateTime(date, time)

  end function tim_yyyymmddhhToDatestampWithyyyymmddhh

  !----------------------------------------------------------------------------------------
  ! tim_getValidDateTimeFromList
  !----------------------------------------------------------------------------------------
  subroutine tim_getValidDateTimeFromList(headDateValues, headTimeValues, validDate, validtime)
    !
    !:Purpose: From a collection of datetimes, find the assimilation window valid date and time
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: headDateValues(:) ! Array of dates (format 'YYYYDDMM') coherent with 'headTimeValues'
    integer, intent(in)  :: headTimeValues(:) ! Array of time (format 'HHMMSS') coherent with 'headDateValues'
    integer, intent(out) :: validDate         ! Assimilation window valid date
    integer, intent(out) :: validTime         ! Assimilation window valid time

    ! Locals:
    integer                 :: numDates, numWindowsPerDay, windowIndex
    integer                 :: windowBoundaryMin, windowBoundaryMax, validTimeMin, validTimeMax, validDateMin, validDateMax
    integer(8)              :: dateTimeMin, dateTimeMax, timeMin, timeMax, dateMin, dateMax
    integer(8), allocatable :: dateTimeValues(:), windowBoundaries(:)
    integer                 :: dateStampIn, dateStampOut

    call tim_readNml()

    numDates = size(headDateValues)
    checkNumDates: if (numDates > 0) then
      write(*,*) 'tim_getValidDateTimeFromList: check inputs: time min/max = ', minval(headTimeValues), maxval(headTimeValues)
      write(*,*) 'tim_getValidDateTimeFromList: check inputs: date min/max = ', minval(headDateValues), maxval(headDateValues)
      allocate(dateTimeValues(size(headDateValues)))
      dateTimeValues(:) = headDateValues(:)
      dateTimeValues(:) = 10000*dateTimeValues(:) + headTimeValues(:)

      dateTimeMin = minval(dateTimeValues(:))
      dateTimeMax = maxval(dateTimeValues(:))
      deallocate(dateTimeValues)
      dateMin = dateTimeMin/10000
      dateMax = dateTimeMax/10000
      timeMin = dateTimeMin - 10000*(dateTimeMin/10000)
      timeMax = dateTimeMax - 10000*(dateTimeMax/10000)
      ! convert from hhmm to just minutes: hhmm - 100*(hhmm/100) + 60*(hhmm/100)
      timeMin = timeMin - 100*(timeMin/100) + 60*(timeMin/100)
      timeMax = timeMax - 100*(timeMax/100) + 60*(timeMax/100)
      write(*,*) 'tim_getValidDateTimeFromList: min/max DateTime             = ', dateTimeMin, dateTimeMax
      write(*,*) 'tim_getValidDateTimeFromList: min/max time (in minutes)    = ', timeMin, timeMax
      if (tim_windowSize < 24.0d0) then
        numWindowsPerDay = nint(24.0/tim_windowSize)

        ! define boundaries between assimilation windows in minutes relative to 0UTC
        allocate(windowBoundaries(0:numWindowsPerDay))
        do windowIndex = 0, numWindowsPerDay
          ! example for windowSize=6h, boundaries = -3h,+3h,+9h,+15h,+21h
          windowBoundaries(windowIndex) = nint(windowIndex*60.0*tim_windowSize - 60.0*tim_windowSize/2.0)
        end do
        write(*,*) 'tim_getValidDateTimeFromList: boundaries (in minutes) = ', windowBoundaries(:)

        ! find left boundary of window where timeMin/Max are located
        windowBoundaryMin = -1
        windowBoundaryMax = -1
        do windowIndex = 0, numWindowsPerDay
          if (timeMin >= windowBoundaries(windowIndex)) then
            windowBoundaryMin = windowIndex
          end if
          if (timeMax >= windowBoundaries(windowIndex)) then
            windowBoundaryMax = windowIndex
          end if
        end do

        ! find validTimeMin/Max from left boundary
        validTimeMin = nint((windowBoundaries(windowBoundaryMin) + 60.0*tim_windowSize/2.0)/60.0)
        if (validTimeMin >= 24) then
          validTimeMin = 0
          dateStampIn = tim_yyyymmddhhToDatestamp(int(dateMin,4), validTimeMin)
          call incdat(dateStampOut, dateStampIn, 24) ! add 1 day to get validDate
          call tim_dateStampToYYYYMMDDHH(dateStampOut, validDateMin, validTimeMin)
          validTimeMin = 0
        else
          validDateMin = int(dateMin,4)
        end if
        validTimeMax = nint((windowBoundaries(windowBoundaryMax) + 60.0*tim_windowSize/2.0)/60.0)
        if (validTimeMax >= 24) then
          validTimeMax = 0
          dateStampIn = tim_yyyymmddhhToDatestamp(int(dateMax,4), validTimeMax)
          call incdat(dateStampOut, dateStampIn, 24) ! add 1 day to get validDate
          call tim_dateStampToYYYYMMDDHH(dateStampOut, validDateMax, validTimeMax)
          validTimeMax = 0
        else
          validDateMax = int(dateMax,4)
        end if
        write(*,*) 'tim_getValidDateTimeFromList: date from Min/Max = ', validDateMin, validDateMax
        write(*,*) 'tim_getValidDateTimeFromList: hour from Min/Max = ', validTimeMin, validTimeMax
        if (validTimeMin /= validTimeMax) call rti_abort('validTimeMin/Max not equal')
        validTime = validTimeMin
        validDate = validDateMin
        deallocate(windowBoundaries)
      else
        write(*,*) 'tim_getValidDateTimeFromList: WARNING: window size equal or greater than 1 day, cannot get dateStamp'
        validTime = 0
        validDate = 0
      end if
    end if checkNumDates

  end subroutine tim_getValidDateTimeFromList

  !----------------------------------------------------------------------------------------
  ! tim_getHoursBetweenDates(currentDate, referenceDate)
  !----------------------------------------------------------------------------------------
  function tim_getHoursBetweenDates(currentDateStamp, referenceDate) result(numberHours)
    !
    ! :Purpose: to compute number of hours between current and reference date
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: currentDateStamp ! current datestamp
    integer, intent(in)  :: referenceDate    ! date in print format yyyyddmm defined in 'namEnsPostProcModule' namelist

    ! Results:
    integer :: numberHours      ! number of hours between current and reference date

    ! Locals:
    integer :: refDateStamp

    write(*,*) 'tim_getHoursSinceReferenceDate: reference date: ', referenceDate
    refDateStamp = tim_yyyymmddhhToDatestamp(referenceDate, 0)
    write(*,*) 'tim_getHoursSinceReferenceDate: reference datestamp: ', refDateStamp

    ! Difference (in hours) between current date and reference date
    call difdat(currentDateStamp, refDateStamp, numberHours)
    write(*,*) 'tim_getHoursSinceReferenceDate: difference in hours: ', numberHours

  end function tim_getHoursBetweenDates

  !----------------------------------------------------------------------------------------
  ! tim_getSecondsBetweenDates(currentDate, referenceDate)
  !----------------------------------------------------------------------------------------
  subroutine tim_getSecondsBetweenDates(currentDateStamp, referenceDate, numberSeconds)
    !
    ! :Purpose: to compute number of seconds in integer(8) between current and reference date.
    !
    implicit none

    integer,    intent(in)  :: currentDateStamp ! current datestamp
    integer,    intent(in)  :: referenceDate    ! date in print format yyyyddmm defined in 'namEnsPostProcModule' namelist
    integer(8), intent(out) :: numberSeconds    ! number of seconds between current and reference date

    ! Locals:
    integer :: numberHours

    numberHours = tim_getHoursBetweenDates(currentDateStamp, referenceDate)
    numberSeconds = int(numberHours * 3600.0d0 , 8)

  end subroutine tim_getSecondsBetweenDates


  !----------------------------------------------------------------------------------------
  ! tim_numberOfDaysInMonth (private)
  !----------------------------------------------------------------------------------------
  function tim_numberOfDaysInMonth(dateStamp,yyyy,mm,dd,time) result(ndays)
    !
    ! :Purpose: Compute the number of days in the month given by 'dateStamp'
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: dateStamp ! datetimestamp in CMC format
    integer, intent(in) :: yyyy      ! year extracted from 'datestamp'
    integer, intent(in) :: mm        ! month extracted from 'datestamp'
    integer, intent(in) :: dd        ! day of month extracted from 'datestamp'
    integer, intent(in) :: time      ! time extracted from 'datestamp' in format 'HHMM'
    ! Results:
    integer :: ndays

    ! Locals:
    integer :: nextdate_year, nextdate_month, nextdate
    integer :: newdate, hours

    ! To compute the number of days in the month, we will compute the
    ! time difference between 'datestamp' and the date in one month from 'datestamp'

    if (mm == 12) then
      nextdate_year = yyyy+1
      nextdate_month = 1
    else
      nextdate_year = yyyy
      nextdate_month = mm + 1
    end if

    !! Build the new date from 'nextdate_year', 'nextdate_month' and 'dd'
    nextdate = tim_dateFromYYYYMMDD(nextdate_year, nextdate_month, dd)
    !! Compute a new datetimestamp from 'nextdate' and original 'time'
    newdate = tim_yyyymmddhhToDatestampWithDateTime(nextdate, time)
    !! Get the numbers of hours between 'dateStamp' and 'newdate'
    hours = tim_getHoursBetweenDates(newdate, dateStamp)

    !! From the number of hours, transform into the number of days
    ndays = hours/24

  end function tim_numberOfDaysInMonth

  !----------------------------------------------------------------------------------------
  ! tim_dateFromYYYYMMDD (private)
  !----------------------------------------------------------------------------------------
  function tim_dateFromYYYYMMDD(year,month,day) result(date)
    !
    ! :Purpose: Transform the separate 'yyyy', 'mm' and 'dd' to an
    !           integer representing this information
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: year  ! year extracted from 'datestamp'
    integer, intent(in) :: month ! month extracted from 'datestamp'
    integer, intent(in) :: day   ! day of month extracted from 'datestamp'

    ! Results:
    integer :: date

    date = year * 10000 + month * 100 + day

  end function tim_dateFromYYYYMMDD

end module timeCoord_mod

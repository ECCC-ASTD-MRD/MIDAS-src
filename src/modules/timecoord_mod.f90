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

module timeCoord_mod
  ! MODULE timeCoord (prefix='tim' category='7. Low-level data objects')
  !
  ! :Purpose: To store public variables and procedures related to the time
  !           coordinate.
  !
  use midasMpi_mod
  use varNameList_mod
  use utilities_mod
  implicit none
  save
  private
  
  ! public variables
  public :: tim_dstepobs, tim_dstepobsinc, tim_windowsize
  public :: tim_nstepobs, tim_nstepobsinc, tim_referencetime
  public :: tim_fullyUseExtremeTimeBins
  ! public procedures
  public :: tim_setup, tim_initialized
  public :: tim_getDateStamp, tim_setDateStamp, tim_getStampList, tim_getStepObsIndex
  public :: tim_getDateStampFromFile, tim_dateStampToYYYYMMDDHH, tim_getValidDateTimeFromList

  character(len=4) :: varNameForDate
  character(len=6) :: tim_referencetime
  real(8)   :: tim_dstepobs
  real(8)   :: tim_dstepobsinc
  real(8)   :: tim_windowsize
  integer   :: tim_nstepobs
  integer   :: tim_nstepobsinc
  logical   :: tim_fullyUseExtremeTimeBins
  integer   :: datestamp = 0  ! datestamp is usually the centre of time window
  logical   :: initialized = .false.

  integer, external :: get_max_rss

contains

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

    ! locals
    integer :: nulnam, ierr, fnom, fclos
    logical, save :: firstCall = .true.

    ! namelist variables:
    real(8) :: dstepobs      ! time step length for background state (in hours)
    real(8) :: dstepobsinc   ! time step length for increment and/or B matrix (in hours)
    real(8) :: dwindowsize   ! length of assimilation window (in hours)
    character(len=6) :: referencetime  ! location of 'date' within the window: 'middle' or 'start'
    logical :: fullyUseExtremeTimeBins ! choose to use full-size bins at both ends of window (usually only half size)

    NAMELIST /NAMTIME/dstepobs, dstepobsinc, dwindowsize, referencetime, fullyUseExtremeTimeBins

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
    read(nulnam, nml=namtime, iostat=ierr)
    if (ierr /= 0) call utl_abort('tim_readNml: Error reading namelist')
    if (mmpi_myid == 0) write(*,nml=namtime)
    ierr = fclos(nulnam)

    ! Set the module variables for timestep length, number of timesteps and window length
    tim_dstepobs      = dstepobs
    tim_dstepobsinc   = dstepobsinc
    tim_windowsize    = dwindowsize
    tim_referencetime = referenceTime
    tim_fullyUseExtremeTimeBins = fullyUseExtremeTimeBins

    if (tim_fullyUseExtremeTimeBins .and. tim_referencetime=="middle") then
      write(*,*) "Warning: tim_fullyUseExtremeTimeBins==.true. and tim_referencetime=='middle' is a non-standard combination"
      write(*,*) "Is it really what you want ?"
    end if

    if (dstepobs > dwindowsize) then
      if (mmpi_myid == 0) write(*,*) 'tim_readNml: dstepobs>dwindowsize. Reset to dwindowsize value.'
      tim_dstepobs = tim_windowsize 
    end if
    if (dstepobsinc > dwindowsize) then
      if (mmpi_myid == 0) write(*,*) 'tim_readNml: dstepobsinc>dwindowsize. Reset to dwindowsize value.'
      tim_dstepobsinc = tim_windowsize 
    end if 
    if (tim_referencetime == "middle") then
      tim_nstepobs    = 2*nint(((tim_windowsize - tim_dstepobs)/2.d0)/tim_dstepobs) + 1
      tim_nstepobsinc = 2*nint(((tim_windowsize - tim_dstepobsinc)/2.d0)/tim_dstepobsinc) + 1
    end if

    if (trim(tim_referencetime) == "start") then
      tim_nstepobs = max(nint(tim_windowsize/tim_dstepobs), 1)
      tim_nstepobsinc = max(nint(tim_windowsize/tim_dstepobsinc), 1)
    end if

  end subroutine tim_readNml


  subroutine tim_getDateStampFromEnvVar(dateStamp)
    !
    !:Purpose: Determine the date from the environment variable MIDAS_DATE.
    !
    implicit none

    ! arguments:
    integer, intent(inout) :: dateStamp

    ! locals:
    integer    :: lengthValidDateStr, status, datePrint, timePrint, imode, ierr
    integer(8) :: dateTimePrint
    character(len=256) :: validDateStr
    integer    :: newdate

    status = 0
    call get_environment_variable('MIDAS_DATE',validDateStr,lengthValidDateStr,status,.true.)

    if (status > 1) then
      call utl_abort('tim_getDateStampFromEnvVar: Problem when getting the environment variable MIDAS_DATE')
    end if
    if (status == 1) then
      write(*,*) 'tim_getDateStampFromEnvVar: WARNING: The environment variable MIDAS_DATE has not been detected!'
      !call utl_abort('tim_getDateStampFromEnvVar: The environment variable MIDAS_DATE has not been detected!')
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
      call utl_abort('tim_getDateStampFromEnvVar: Unexpected length of variable MIDAS_DATE')
    end if

    ! convert to CMC dateStamp
    imode = 3 ! printable to stamp
    ierr = newdate(datestamp, datePrint, timePrint, imode)

    write(*,*) 'tim_getDateStampFromEnvVar: envVar, validDate, dateStamp = ', trim(validDateStr), dateTimePrint, dateStamp

  end subroutine tim_getDateStampFromEnvVar


  subroutine tim_setup(fileNameForDate_opt)
    !
    ! :Purpose: Setup of obs time window size and related trial field 
    !           time step for OmP determination. 
    !
    implicit none

    ! arguments
    character(len=*), optional :: fileNameForDate_opt

    ! locals
    integer :: ierr, newdate, imode, prntdate, prnttime
    integer :: dateStampEnvVar

    call tim_readNml()

    if (initialized) then
      write(*,*) 'tim_setup: already initialized, just return'
      return
    end if

    ! First try to set dateStamp from MIDAS_DATE
    dateStampEnvVar = 0
    call tim_getDateStampFromEnvVar(dateStampEnvVar)

    ! Possibly set the datestamp (except when set later from burp files)
    if (dateStampEnvVar /= 0) then
      write(*,*) 'tim_setup: ====================================================='
      write(*,*) 'tim_setup: DATESTAMP set by value in supplied MIDAS_DATE'
      write(*,*) 'tim_setup: ====================================================='
      dateStamp = dateStampEnvVar
      imode = -3 ! stamp to printable
      ierr = newdate(datestamp, prntdate, prnttime, imode)
      write(*,*) 'tim_setup: printdate = ',prntdate
      write(*,*) 'tim_setup: printtime = ',prnttime
      write(*,*) 'tim_setup: datestamp = ',datestamp      
    else if (present(fileNameForDate_opt)) then
      write(*,*) 'tim_setup: ====================================================='
      write(*,*) 'tim_setup: DATESTAMP set by value in supplied file'
      write(*,*) 'tim_setup: ====================================================='
      datestamp = tim_getDatestampFromFile(fileNameForDate_opt)
      imode = -3 ! stamp to printable
      ierr = newdate(datestamp, prntdate, prnttime, imode)
      write(*,*) 'tim_setup: printdate = ',prntdate
      write(*,*) 'tim_setup: printtime = ',prnttime
      write(*,*) 'tim_setup: datestamp = ',datestamp
    else
      write(*,*) 'tim_setup: =========================================================='
      write(*,*) 'tim_setup: DATESTAMP not set in this subroutine, use tim_setDateStamp'
      write(*,*) 'tim_setup: =========================================================='      
    end if

    if (mmpi_myid == 0) write(*,*) 'tim_setup: dobs_windowsize=',tim_windowsize
    if (mmpi_myid == 0) write(*,*) 'tim_setup: dstepobs       =',tim_dstepobs
    if (mmpi_myid == 0) write(*,*) 'tim_setup: nstepobs       =',tim_nstepobs
    if (mmpi_myid == 0) write(*,*) 'tim_setup: dstepobsinc    =',tim_dstepobsinc
    if (mmpi_myid == 0) write(*,*) 'tim_setup: nstepobsinc    =',tim_nstepobsinc
    if (mmpi_myid == 0) write(*,*) 'tim_setup: tim_referencetime   =',tim_referencetime

    initialized = .true.

  end subroutine tim_setup


  function tim_initialized() result(initialized_out)
    implicit none
    logical initialized_out

    initialized_out = initialized

  end function tim_initialized


  function tim_getDatestampFromFile(fileName, varNameForDate_opt) result(dateStamp_out)
    !
    ! :Purpose: to extract the dateStamp from the supplied file.
    !
    implicit none

    ! arguments
    character(len=*)           :: fileName
    character(len=*), optional :: varNameForDate_opt
    integer :: dateStamp_out

    ! locals
    integer :: nulFile, ierr
    integer, parameter :: maxNumDates = 2000
    integer :: numDates, ikeys(maxNumDates), varIndex
    integer :: fnom, fstouv, fstinl, fstprm, fstfrm, fclos, newdate
    integer :: prntdate, prnttime, imode, windowIndex, windowsPerDay, dateStamp_tmp
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
        call utl_abort('tim_getDateStampFromFile: file not found '//trim(fileName))
      end if

      ! Determine variable to use for the date (default is P0)
      varNameForDate = 'P0'
      if (present(varNameForDate_opt)) then
      
	varNameForDate = trim(varNameForDate_opt)
        write(*,*) 'tim_getDateStampFromFile: defining dateStamp from the variable = ', varNameForDate
      
      ! If P0 not present, look for another suitable variable in the file
      else if (.not. utl_varNamePresentInFile(varNameForDate,fileName_opt=trim(fileName))) then
      
        foundVarNameInFile = .false.
        do varIndex = 1, vnl_numvarmax
          varNameForDate = vnl_varNameList(varIndex)
          ! check if variable is in the file
          if (.not. utl_varNamePresentInFile(varNameForDate,fileName_opt=trim(fileName))) cycle
          foundVarNameInFile = .true.
          exit      
        end do

        if (.not. foundVarNameInFile) then
          call utl_abort('tim_getDateStampFromFile: NO variables found in the file!!!')
        end if
	
      end if

      ! Extract the datestamp from the file
      nulFile = 0
      ierr = fnom(nulFile,trim(fileName),'RND+OLD+R/O',0)
      ierr = fstouv(nulFile,'RND+OLD')
      ierr = fstinl(nulFile,ini,inj,ink,-1,' ',-1,-1,-1,' ', &
                    trim(varNameForDate),ikeys,numdates,maxnumdates)
      if (ikeys(1) <= 0) then
        call utl_abort('tim_getDateStampFromFile: Could not find variable ' //  &
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

    call rpn_comm_bcast(dateStamp_out, 1, 'MPI_INTEGER', 0, 'GRID', ierr)

    if (tim_referenceTime == 'middle') then
      ! Modify date to ensure that it corresponds to the middle of the window
      ! Note: For this, we have to assume that the date in the file
      !       does NOT correspond to the final time of the window
      imode = -3 ! stamp to printable, time is HHMMSShh
      ierr = newdate(datestamp_out, prntdate, prnttime, imode)
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
        write(*,*) 'windowsPerDay, fileHour=', windowsPerDay, fileHour
        call utl_abort('tim_getDateStampFromFile: could not determine assimilation window position')
      end if

      ! handle special case when window centered on hour 24
      if (nint(middleHour) == 24) then
        ! add 24h to dateStamp and recompute prntdate
        dateStamp_tmp = dateStamp_out
        call incdatr(dateStamp_out, dateStamp_tmp, 24.0d0)
        imode = -3 ! stamp to printable, only prntdate will be used
        ierr = newdate(datestamp_out, prntdate, prnttime, imode)

        ! subtract 24h from middleHour
        middleHour = 0.0d0
      end if

      prnttime = nint(middleHour) * 1000000
      imode = 3 ! printable to stamp
      ierr = newdate(datestamp_out, prntdate, prnttime, imode)

    end if

  end function tim_getDateStampFromFile


  subroutine tim_setDatestamp(datestamp_in)
    !
    ! :Purpose: to control access to the minimization object.  Sets the date
    !           of the window centre of analysis validity to the indicated value.
    !
    implicit none
    integer, intent(in) :: datestamp_in

    if (.not.initialized) call utl_abort('tim_setDateStamp: module not initialized')

    datestamp = datestamp_in

  end subroutine tim_setDatestamp


  function tim_getDatestamp() result(datestamp_out)
    !
    ! :Purpose: to control access to the minimization object.  Returns the date
    !           of the window centre of analysis validity.
    !
    implicit none
    integer :: datestamp_out

    if (.not.initialized) call utl_abort('tim_getDateStamp: module not initialized')

    datestamp_out = datestamp

  end function tim_getDatestamp


  subroutine tim_getStampList(datestamplist, numStep, referenceDateStamp)
    !
    ! :Purpose: Compute a list of STAMPS corresponding to stepobs time
    !
    implicit none

    ! arguments
    integer, intent(in)  :: numStep ! number of step obs
    integer, intent(in)  :: referenceDateStamp ! Synoptic time
    integer, intent(out) :: datestamplist(numStep) 

    ! locals
    integer :: stepIndex
    real(8) :: dldelt ! delta time in hours between middle time and each step
    real(8) :: dtstep ! delta time in hours between step obs

    if (.not. initialized) call utl_abort('tim_getStampList: module not initialized')

    if (referenceDateStamp == -1) then

      if (mmpi_myid == 0) write(*,*) 'tim_getStampList: datestamp is not specified, keep as -1'
      datestamplist(:) = -1

    else

      if (tim_referencetime == "middle") then
        if (numStep > 1) then
          dtstep = tim_windowsize/(real(numStep-1,8))
        else
          dtstep = tim_windowsize
        end if

        do stepIndex = 1, numStep
          dldelt = (stepIndex - ((numStep-1)/2 + 1)) * dtstep
          call incdatr(datestamplist(stepIndex), referenceDateStamp, dldelt)
        end do
      end if

      if (trim(tim_referencetime) == "start") then
     
        dtstep = tim_windowsize/(real(numStep,8))
     
        do stepIndex = 1, numStep
          dldelt = (stepIndex - 1) * dtstep
          call incdatr(datestamplist(stepIndex), referenceDateStamp, dldelt)
        end do

      end if

    end if

  end subroutine tim_getStampList


  subroutine tim_getStepObsIndex(dnstepobs, referenceDateStamp, obsYYYMMDD, obsHHMM, numStep)
    !
    ! :Purpose: Return the stepobs index as a real number (-1.0 if out of range)
    !

    implicit none

    ! arguments
    real(8), intent(out) :: dnstepobs          ! number of stepobs from reference time
    integer, intent(in)  :: referenceDateStamp    ! Synop CMC date-time stamp
    integer, intent(in)  :: obsYYYMMDD         ! Obs date YYYYMMDD
    integer, intent(in)  :: obsHHMM            ! Obs time HHMM
    integer, intent(in)  :: numStep            ! number of stepobs in assimilation window

    ! locals
    integer :: newdate, istat, imode
    real(8) :: dddt      ! delta time in hours
    integer :: istobs    ! obs CMC date-time stamp
    integer :: itobs     ! obs time HHMMSShh
    real(8) :: dlhours   ! delta time from synop time

    if (.not. initialized) call utl_abort('tim_getStepObsIndex: module not initialized')

    ! Building observation stamp
    imode = 3 ! printable to stamp
    itobs = obsHHMM * 10000
    istat = newdate(istobs, obsYYYMMDD, itobs, imode)

    ! Difference (in hours) between obs time and reference time
    call difdatr(istobs, referenceDateStamp, dlhours)

    if (numStep > 1) then
      ! FGAT: more than 1 trial field in time window
      if (tim_referencetime == "middle") then
        dddt = tim_windowsize / (real(numStep-1,8))
      else
        dddt = tim_windowsize / (real(numStep,8))
      end if
      dnstepobs = dlhours / dddt ! number of step obs from reference (e.g. synoptic)
      if (tim_referencetime == "middle") then
        dnstepobs = dnstepobs + real((numStep+1)/2,8)
      end if
      if (trim(tim_referencetime) == "start") then
        dnstepobs = dnstepobs + 1.d0
      end if
      if (dnstepobs < 0.5d0.or.dnstepobs > (0.5d0+real(numStep,8))) dnstepobs = -1.0d0

    else
      ! 3D: only 1 trial field in time window
      if (tim_referencetime == "middle") then
        if (dlhours < -tim_windowsize/2.0D0 .or. dlhours > tim_windowsize/2.0D0) then
          ! outside time window
          dnstepobs = -1.0d0
        else
          ! inside time window
          dnstepobs = 1.0d0
        end if
      else
        dddt = tim_windowsize
        if (dlhours < -dddt/2.0d0 .or. dlhours > tim_windowsize + dddt/2.d0) then
          ! outside time window
          dnstepobs = -1.0d0
        else
          ! inside time window
          dnstepobs = 1.0d0
        end if
      endif
      
    end if

  end subroutine tim_getStepObsIndex
  
  !----------------------------------------------------------------------------------------
  ! tim_dateStampToDDMMYYYY
  !----------------------------------------------------------------------------------------
  subroutine tim_dateStampToYYYYMMDDHH(dateStamp, prnttime, dd, mm, ndays, yyyy, verbose_opt)
    !
    ! :Purpose: to get day (DD), month (MM), number of days in this month 
    !           and year (YYYY) from dateStamp
    !  
    
    implicit none
  
    ! arguments
    integer, intent(in)           :: dateStamp
    integer, intent(inout)        :: prnttime, dd, mm, ndays, yyyy
    logical, intent(in), optional :: verbose_opt
    
    ! locals
    character(len=8)            :: yyyymmdd
    character(len=3), parameter :: months(12) = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    integer                     :: ndaysM(12)    
    integer                     :: imode, ierr, newdate, prntdate
    logical                     :: verbose = .True.

    ndaysM(:) = [   31,    28,    31,    30,    31,    30,    31,    31,    30,    31,    30,    31]
    
    imode = -3 ! stamp to printable
    ierr = newdate(dateStamp, prntdate, prnttime, imode)
    write(yyyymmdd,'(i8)') prntdate
    read (yyyymmdd(7:8), '(i)') dd
    read (yyyymmdd(5:6), '(i)') mm
    read (yyyymmdd(1:4), '(i)') yyyy

    ! leap year    
    if (mm == 2 .and. mod(yyyy,4)==0) ndaysM(mm) = 29
    ndays = ndaysM(mm)
    
    if (present(verbose_opt)) verbose = verbose_opt
    
    if(verbose) then
      write(*,*) 'tim_dateStampToYYYYMMDDHH:  date = ', prntdate
      write(*,*) 'tim_dateStampToYYYYMMDDHH:  year = ', yyyy
      write(*,'(a,i5,a,i5,a)') 'tim_dateStampToYYYYMMDDHH: month = ', mm, ' ( '// months(mm)//' where there are ', ndays, ' days)' 
      write(*,*) 'tim_dateStampToYYYYMMDDHH:   day = ', dd
      write(*,*) 'tim_dateStampToYYYYMMDDHH:  time = ', prnttime
    end if     
  
  end subroutine tim_dateStampToYYYYMMDDHH

  !----------------------------------------------------------------------------------------
  ! tim_getValidDateTimeFromList
  !----------------------------------------------------------------------------------------
  subroutine tim_getValidDateTimeFromList(headDateValues, headTimeValues, validDate, validtime)
    implicit none

    ! arguments:
    integer, intent(in)  :: headDateValues(:)
    integer, intent(in)  :: headTimeValues(:)
    integer, intent(out) :: validDate
    integer, intent(out) :: validTime

    ! locals:
    integer                 :: numDates, numWindowsPerDay, windowIndex, timeMin, timeMax, dateMin, dateMax
    integer                 :: windowBoundaryMin, windowBoundaryMax, validTimeMin, validTimeMax, validDateMin, validDateMax 
    integer(8)              :: dateTimeMin, dateTimeMax
    integer(8), allocatable :: dateTimeValues(:), windowBoundaries(:)
    integer                 :: ier, imode
    integer                 :: newdate, dateStampIn, dateStampOut

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
          imode = 3
          ier = newdate(dateStampIn, dateMin, validTimeMin, imode)
          call incdat(dateStampOut, dateStampIn, 24) ! add 1 day to get validDate
          imode = -3
          ier = newdate(dateStampOut, validDateMin, validTimeMin, imode)
          validTimeMin = 0
        else
          validDateMin = dateMin
        end if
        validTimeMax = nint((windowBoundaries(windowBoundaryMax) + 60.0*tim_windowSize/2.0)/60.0)
        if (validTimeMax >= 24) then
          validTimeMax = 0
          imode = 3
          ier = newdate(dateStampIn, dateMax, validTimeMax, imode)
          call incdat(dateStampOut, dateStampIn, 24) ! add 1 day to get validDate
          imode = -3
          ier = newdate(dateStampOut, validDateMax, validTimeMax, imode)
          validTimeMax = 0
        else
          validDateMax = dateMax
        end if
        write(*,*) 'tim_getValidDateTimeFromList: date from Min/Max = ', validDateMin, validDateMax
        write(*,*) 'tim_getValidDateTimeFromList: hour from Min/Max = ', validTimeMin, validTimeMax
        if (validTimeMin /= validTimeMax) call utl_abort('validTimeMin/Max not equal')
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

end module timeCoord_mod

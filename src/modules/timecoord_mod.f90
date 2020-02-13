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
  ! MODULE timeCoord (prefix='tim' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: To store public variables and procedures related to the time
  !           coordinate.
  !
  use mpi_mod
  use mpivar_mod
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
  public :: tim_getDateStampFromFile

  character(len=4) :: varNameForDate = 'P0'
  character(len=6) ::  tim_referencetime
  real(8)   :: tim_dstepobs
  real(8)   :: tim_dstepobsinc
  real(8)   :: tim_windowsize
  integer   :: tim_nstepobs
  integer   :: tim_nstepobsinc
  logical   :: tim_fullyUseExtremeTimeBins
  integer   :: datestamp = 0      ! window centre of analysis validity
  logical   :: initialized = .false.

  integer, external :: get_max_rss

contains

  subroutine tim_setup(fileNameForDate_opt)
    !
    ! :Purpose: Setup of obs time window size and related trial field 
    !           time step for OmP determination. 
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

    ! arguments
    character(len=*), optional :: fileNameForDate_opt

    ! locals
    integer :: nulnam, ierr, fnom, fclos, newdate, imode, prntdate, prnttime, date
    real(8) :: dstepobs,dstepobsinc,dwindowsize
    character(len=6) :: referencetime
    logical :: fullyUseExtremeTimeBins
    NAMELIST /NAMTIME/dstepobs, dstepobsinc, dwindowsize, date, referencetime, fullyUseExtremeTimeBins

    if ( initialized ) then
      write(*,*) 'tim_setup: already initialized, just return'
      return
    end if

    ! Set default values for namelist variables
    dstepobs       = 6.0d0
    dstepobsinc    = 6.0d0      
    dwindowsize    = 6.0d0     
    date           = 0
    referenceTime = 'middle'
    fullyUseExtremeTimeBins = .false.

    ! Read the namelist
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=namtime, iostat=ierr)
    if (ierr /= 0) call utl_abort('tim_setup: Error reading namelist')
    if (mpi_myid == 0) write(*,nml=namtime)
    ierr = fclos(nulnam)

    ! Set the module variables for timestep length, number of timesteps and window length
    tim_dstepobs    = dstepobs
    tim_dstepobsinc = dstepobsinc
    tim_windowsize  = dwindowsize
    tim_referencetime = referenceTime
    tim_fullyUseExtremeTimeBins = fullyUseExtremeTimeBins

    if ( tim_fullyUseExtremeTimeBins .and. tim_referencetime=="middle"  ) then
      write(*,*) "Warning: tim_fullyUseExtremeTimeBins==.true. and tim_referencetime=='middle' is a non-standard combination"
      write(*,*) "Is is really what you want ?"
    end if

    if (dstepobs > dwindowsize) then
      if (mpi_myid == 0) write(*,*) 'tim_setup: dstepobs>dwindowsize. Reset to dwindowsize value.'
      tim_dstepobs = tim_windowsize 
    end if
    if (dstepobsinc > dwindowsize) then
      if (mpi_myid == 0) write(*,*) 'tim_setup: dstepobsinc>dwindowsize. Reset to dwindowsize value.'
      tim_dstepobsinc = tim_windowsize 
    end if 
    if ( tim_referencetime == "middle" ) then
      tim_nstepobs    = 2*nint(((tim_windowsize - tim_dstepobs)/2.d0)/tim_dstepobs) + 1
      tim_nstepobsinc = 2*nint(((tim_windowsize - tim_dstepobsinc)/2.d0)/tim_dstepobsinc) + 1
    end if

    if ( trim(tim_referencetime) == "start" ) then
      tim_nstepobs = max(nint(tim_windowsize/tim_dstepobs), 1)
      tim_nstepobsinc = max(nint(tim_windowsize/tim_dstepobsinc), 1)
    end if


    ! Possibly set the datestamp (except when set later from burp files)
    if ( present(fileNameForDate_opt) ) then
      write(*,*) 'tim_setup: ====================================================='
      write(*,*) 'tim_setup: DATESTAMP set by value in supplied file'
      write(*,*) 'tim_setup: ====================================================='
      datestamp = tim_getDatestampFromFile(fileNameForDate_opt)
      imode = -3 ! stamp to printable
      ierr = newdate(datestamp, prntdate, prnttime, imode)
      write(*,*) 'tim_setup: printdate = ',prntdate
      write(*,*) 'tim_setup: printtime = ',prnttime
      write(*,*) 'tim_setup: datestamp = ',datestamp
    else if ( date /= 0 ) then
      write(*,*) 'tim_setup: ====================================================='
      write(*,*) 'tim_setup: DATESTAMP set by value in namelist'
      write(*,*) 'tim_setup: ====================================================='
      prntdate = date/100
      prnttime = (date - prntdate*100) * 1000000
      imode = 3 ! printable to stamp
      ierr = newdate(datestamp, prntdate, prnttime, imode)
      write(*,*) 'tim_setup: printdate = ',prntdate
      write(*,*) 'tim_setup: printtime = ',prnttime
      write(*,*) 'tim_setup: datestamp = ',datestamp
    end if

    if (mpi_myid == 0) write(*,*) 'tim_setup: dobs_windowsize=',tim_windowsize
    if (mpi_myid == 0) write(*,*) 'tim_setup: dstepobs       =',tim_dstepobs
    if (mpi_myid == 0) write(*,*) 'tim_setup: nstepobs       =',tim_nstepobs
    if (mpi_myid == 0) write(*,*) 'tim_setup: dstepobsinc    =',tim_dstepobsinc
    if (mpi_myid == 0) write(*,*) 'tim_setup: nstepobsinc    =',tim_nstepobsinc
    if (mpi_myid == 0) write(*,*) 'tim_setup: tim_referencetime   =',tim_referencetime

    initialized = .true.

  end subroutine tim_setup


  function tim_initialized() result(initialized_out)
    implicit none
    logical initialized_out

    initialized_out = initialized

  end function tim_initialized


  function tim_getDatestampFromFile(fileName) result(dateStamp_out)
    !
    ! :Purpose: to extract the dateStamp from the supplied file.
    !
    implicit none

    ! arguments
    character(len=*) :: fileName
    integer :: dateStamp_out

    ! locals
    integer :: nulFile, ierr
    integer, parameter :: maxNumDates = 100
    integer :: numDates, ikeys(maxNumDates)
    integer :: fnom, fstouv, fstinl, fstprm, fstfrm, fclos, newdate
    integer :: prntdate, prnttime, imode, windowIndex, windowsPerDay, dateStamp_tmp
    logical :: fileExists, foundWindow
    real(8) :: leadTimeInHours, windowBegHour, windowEndHour, fileHour, middleHour
    integer :: ideet, inpas, dateStamp_origin, ini, inj, ink, inbits, idatyp, &
               ip1, ip2, ip3, ig1, ig2, ig3, ig4, iswa, ilng, idltf, iubc,   &
               iextra1, iextra2, iextra3
    character(len=2)  :: typvar
    character(len=4)  :: nomvar
    character(len=12) :: etiket
    character(len=1)  :: grtyp

    if ( mpi_myid == 0 ) then

      ! Check if file for any date within the analysis window (except the last) exists
      inquire(file=trim(fileName), exist=fileExists)
      if ( .not.fileExists ) then
        call utl_abort('tim_getDateStampFromFile: File not found')
      end if

      ! Extract the datestamp from the file
      nulFile = 0
      ierr = fnom(nulFile,trim(fileName),'RND+OLD+R/O',0)
      ierr = fstouv(nulFile,'RND+OLD')
      ierr = fstinl(nulFile,ini,inj,ink,-1,' ',-1,-1,-1,' ', &
                    trim(varNameForDate),ikeys,numdates,maxnumdates)
      if ( ikeys(1) <= 0 ) then
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
      if ( fileHour >= windowBegHour .and. fileHour < windowEndHour ) then
        foundWindow = .true.
        middleHour = real(windowIndex,8) * tim_windowsize
        exit window_loop
      end if
    end do window_loop

    if ( .not. foundWindow ) then
      write(*,*) 'windowsPerDay, fileHour=', windowsPerDay, fileHour
      call utl_abort('tim_getDateStampFromFile: could not determine assimilation window position')
    end if

    ! handle special case when window centered on hour 24
    if ( nint(middleHour) == 24 ) then
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

  end function tim_getDateStampFromFile


  subroutine tim_setDatestamp(datestamp_in)
    !
    ! :Purpose: to control access to the minimization object.  Sets the date
    !           of the window centre of analysis validity to the indicated value.
    !
    implicit none
    integer, intent(in) :: datestamp_in

    if ( .not.initialized ) call utl_abort('tim_setDateStamp: module not initialized')

    datestamp = datestamp_in

  end subroutine tim_setDatestamp


  function tim_getDatestamp() result(datestamp_out)
    !
    ! :Purpose: to control access to the minimization object.  Returns the date
    !           of the window centre of analysis validity.
    !
    implicit none
    integer :: datestamp_out

    if ( .not.initialized ) call utl_abort('tim_getDateStamp: module not initialized')

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
     
      dtstep = tim_windowsize/( real(numStep,8) )
     
      do stepIndex = 1, numStep
        dldelt = (stepIndex - 1) * dtstep
        call incdatr(datestamplist(stepIndex), referenceDateStamp, dldelt)
      end do

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

end module timeCoord_mod

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

!--------------------------------------------------------------------------
!! MODULE timeCoord (prefix="tim")
!!
!! *Purpose*: To store public variables and procedures related to the time coordinate.
!!
!--------------------------------------------------------------------------
MODULE timeCoord_mod
  use mpivar_mod
  use obsSpaceData_mod
  use utilities_mod
  implicit none
  save
  private
  
  ! public variables
  public :: tim_dstepobs, tim_dstepobsinc, tim_windowsize
  public :: tim_nstepobs, tim_nstepobsinc
  ! public procedures
  public :: tim_setup, tim_timeBinning, tim_initialized
  public :: tim_sutimeinterp, tim_setTimeInterpWeight, tim_getTimeInterpWeight
  public :: tim_getDateStamp, tim_setDateStamp, tim_getStampList, tim_getStepObsIndex
  public :: tim_getDateStampFromFile

  character(len=4) :: varNameForDate = 'P0'
  real(8)   :: tim_dstepobs
  real(8)   :: tim_dstepobsinc, tim_windowsize
  integer   :: tim_nstepobs
  integer   :: tim_nstepobsinc
  integer, parameter :: mxstepobs = 9 
  real(8), pointer   :: timeInterpWeight(:,:) => NULL() ! weights for linear temporal interpolation of increment to obs times
  integer            :: datestamp = 0    ! window centre of analysis validity
  logical            :: initialized = .false.

contains

  subroutine tim_setup(fileNameForDate_opt)
    !
    ! Purpose: Setup of obs time window size and related trial field 
    !          time step for OmP determination. 
    !  
    ! Revisions:
    !           Y. Rochon ARQI/AQRD July 2016
    !           - Allowance for obs time windows of size different from 6 hours.
    !             Related addition of dwindowsize and tim_windowsize 
    !
    ! Namelist parameters:
    !
    !   dstepobs          time step (hrs) between successive trial fields 
    !                     for use in OmP determation. Set to dwindowsize for
    !                     single trial field, i.e. use of 3dvar instead of 3dvar-FGAT.
    !                     nstepobs = number of trial fields 
    !   dstepobsinc       time step (hrs) between obs groupings in time. Set to
    !                     dwindowsize for use of a single obs group.
    !                     nstepobsinc = number of obs time intervals 
    !   dwindowsize       Time window size (hrs).
    !
    ! Comment:
    !
    ! - Provided dates and number of provided trial field files must be
    !   consistent with nstepobs, dstepobs and dwindowsize with reference datestamp
    !   corresponding to the date of the middle trial field file.
    !
    implicit none
    character(len=*), optional :: fileNameForDate_opt

    integer :: nulnam,ierr,fnom,fclos,newdate,imode,prntdate,prnttime,date
    real(8) :: dstepobs,dstepobsinc,dwindowsize
    NAMELIST /NAMTIME/dstepobs,dstepobsinc,dwindowsize,date

    if ( initialized ) then
      write(*,*) 'tim_setup: already initialized, just return'
      return
    end if

    ! Set default values for namelist variables
    dstepobs       = 6.0d0
    dstepobsinc    = 6.0d0      
    dwindowsize    = 6.0d0     
    date           = 0

    ! Read the namelist
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=namtime, iostat=ierr)
    if(ierr /= 0) call utl_abort('tim_setup: Error reading namelist')
    if(mpi_myid == 0) write(*,nml=namtime)
    ierr = fclos(nulnam)

    ! Set the module variables for timestep length, number of timesteps and window length
    tim_dstepobs    = dstepobs
    tim_dstepobsinc = dstepobsinc
    tim_windowsize  = dwindowsize
    if (dstepobs > dwindowsize) then
        if(mpi_myid == 0) write(*,*) 'tim_setup: dstepobs>dwindowsize. Reset to dwindowsize value.'
        tim_dstepobs = tim_windowsize 
    end if
    if (dstepobsinc > dwindowsize) then
        if(mpi_myid == 0) write(*,*) 'tim_setup: dstepobsinc>dwindowsize. Reset to dwindowsize value.'
        tim_dstepobsinc = tim_windowsize 
    end if      
    tim_nstepobs    = 2*nint(((tim_windowsize - tim_dstepobs)/2.d0)/tim_dstepobs) + 1
    tim_nstepobsinc = 2*nint(((tim_windowsize - tim_dstepobsinc)/2.d0)/tim_dstepobsinc) + 1

    ! Possibly set the datestamp (except when set later from burp files)
    if( present(fileNameForDate_opt) ) then
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
      prnttime = (date - prntdate*100) * 10000
      imode = 3 ! printable to stamp
      ierr = newdate(datestamp, prntdate, prnttime, imode)
      write(*,*) 'tim_setup: printdate = ',prntdate
      write(*,*) 'tim_setup: printtime = ',prnttime
      write(*,*) 'tim_setup: datestamp = ',datestamp
    end if

    if(mpi_myid == 0) write(*,*) 'tim_setup: dobs_windowsize=',tim_windowsize
    if(mpi_myid == 0) write(*,*) 'tim_setup: dstepobs       =',tim_dstepobs
    if(mpi_myid == 0) write(*,*) 'tim_setup: nstepobs       =',tim_nstepobs
    if(mpi_myid == 0) write(*,*) 'tim_setup: dstepobsinc    =',tim_dstepobsinc
    if(mpi_myid == 0) write(*,*) 'tim_setup: nstepobsinc    =',tim_nstepobsinc

    initialized = .true.

  end subroutine tim_setup

  function tim_initialized() result(initialized_out)
    implicit none
    logical initialized_out

    initialized_out = initialized

  end function tim_initialized

  function tim_getDatestampFromFile(fileName) result(dateStamp_out)
    !
    ! object: to extract the dateStamp from the supplied file.
    implicit none

    character(len=*) :: fileName
    integer :: dateStamp_out

    integer :: middleStep, nulFile, ierr
    integer, parameter :: maxNumDates = 100
    integer :: numDates, ikeys(maxNumDates)
    integer :: fnom, fstouv, fstinf, fstinl, fstprm, fstfrm, fclos
    character(len=2)   :: fileNumber
    logical :: fileExists
    real(8) :: leadTimeInHours
    integer :: ideet, inpas, dateStamp_origin, ini, inj, ink, inbits, idatyp, &
               ip1, ip2, ip3, ig1, ig2, ig3, ig4, iswa, ilng, idltf, iubc,   &
               iextra1, iextra2, iextra3
    character(len=2)  :: typvar
    character(len=4)  :: nomvar
    character(len=12) :: etiket
    character(len=1)  :: grtyp

    if ( mpi_myid == 0 ) then

      ! Check if file for middle of analysis window exists
      inquire(file=trim(fileName), exist=fileExists)
      if (.not. fileExists) then
        call utl_abort('tim_getDateStampFromFile: File not found')
      end if

      ! Extract the datestamp from the file
      nulFile = 0
      ierr = fnom(nulFile,trim(fileName),'RND+OLD+R/O',0)
      ierr = fstouv(nulFile,'RND+OLD')
      ierr = fstinl(nulFile,ini,inj,ink,-1,' ',-1,-1,-1,' ', &
                    trim(varNameForDate),ikeys,numdates,maxnumdates)
      if( ikeys(1) <= 0 ) then
        call utl_abort('tim_getDateStampFromFile: Could not find variable ' //  &
                       trim(varNameForDate) // ' in the supplied file')
      end if
      write(*,*) 'tim_getDateStampFromFile: number of dates found = ', numDates
      ierr = fstprm(ikeys((numDates+1)/2), dateStamp_origin, ideet, inpas, ini, inj, &
                    ink, inbits, idatyp, ip1, ip2, ip3, &
                    typvar, nomvar, etiket, grtyp, ig1, ig2, ig3, ig4, &
                    iswa, ilng, idltf, iubc, iextra1, iextra2, iextra3)
      leadTimeInHours = real(ideet*inpas,8)/3600.0d0
      call incdatr(dateStamp_out, dateStamp_origin, leadTimeInHours)

      ierr = fstfrm(nulFile)
      ierr = fclos(nulFile)

    end if

    call rpn_comm_bcast(dateStamp_out, 1, 'MPI_INTEGER', 0, 'GRID', ierr)

  end function tim_getDateStampFromFile


  subroutine tim_timeBinning(obsSpaceData,nstepobs)
    implicit none

    type(struct_obs) :: obsSpaceData
    integer :: nstepobs

    integer :: jstep,jobs,jfamily,jdata,iobs,idatabeg,idatend,nsize,ierr
    integer, allocatable :: idataass(:,:), inumheader(:,:)
    integer, allocatable :: my_idataass(:,:), my_inumheader(:,:)
    integer, parameter   :: nfamily = 10
    character(len=2)     :: familylist(nfamily)
    character(len=256)   :: formatspec,formatspec2
    real(8)              :: stepObsIndex

    if (.not. initialized) call utl_abort('tim_timeBinning: module not initialized')

    allocate(idataass(nfamily,nStepObs+1))
    allocate(my_idataass(nfamily,nStepObs+1))
    my_idataass(:,:) = 0
    allocate(inumheader(nfamily,nStepObs+1))
    allocate(my_inumheader(nfamily,nStepObs+1))
    my_inumheader(:,:) = 0

    familylist(1) = 'UA'
    familylist(2) = 'AI'
    familylist(3) = 'SF'
    familylist(4) = 'TO'
    familylist(5) = 'SW'
    familylist(6) = 'SC'
    familylist(7) = 'PR'
    familylist(8) = 'RO'
    familylist(9) = 'GP'
    familylist(10) = 'CH'

    do jobs = 1, obs_numheader(obsSpaceData)
       call tim_getStepObsIndex(stepObsIndex,tim_getDatestamp(), &
            obs_headElem_i(obsSpaceData,OBS_DAT,jobs), &
            obs_headElem_i(obsSpaceData,OBS_ETM,jobs),nstepobs)
       if(stepObsIndex > 0.0d0) then
          jstep = nint(stepObsIndex)
          idatabeg = obs_headElem_i(obsSpaceData,OBS_RLN,jobs)
          idatend = obs_headElem_i(obsSpaceData,OBS_NLV,jobs) + idatabeg - 1          
          do jfamily = 1, nfamily
             if(obs_getfamily(obsSpaceData,jobs) == familylist(jfamily)) then
                my_inumheader(jfamily,jstep) = my_inumheader(jfamily,jstep)+1
                my_inumheader(jfamily,nStepObs+1) = my_inumheader(jfamily,nStepObs+1)+1
                do jdata = idatabeg, idatend
                   if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,jdata) == 1) then
                      my_idataass(jfamily,jstep) = my_idataass(jfamily,jstep) + 1
                      my_idataass(jfamily,nStepObs+1) = &
                           my_idataass(jfamily,nStepObs+1) + 1
                   end if
                end do
             end if
          end do
       else
          write(*,*) 'TIM_TIMEBINNING: observation outside time window:',jobs,stepObsIndex
       end if
    end do

    formatspec ='(1X,A6,":"'
    do jstep = 1,nStepObs
       formatspec = trim(formatspec)//',1X,I7' ! this is for each time bin
    end do
    formatspec = trim(formatspec)//',1X,I9' ! this is for the total
    formatspec = trim(formatspec)//')'

    formatspec2 = '(1X,A6,":"'
    do jstep = 1,nStepObs
       formatspec2 = trim(formatspec2)//',1X,I7'
    end do
    formatspec2 = trim(formatspec2)//',1X,A9)'

    write(*,*) '-----------------------------------------------------------------'
    write(*,*) 'Distribution of number of headers over stepobs ON LOCAL PROCESSOR'
    write(*,trim(formatspec2)) 'Bin#',(jstep, jstep = 1, nStepObs),'Total'
    do jfamily = 1, nfamily
       write(*,trim(formatspec)) familylist(jfamily),(my_inumheader(jfamily,jstep), &
             jstep = 1, nStepObs+1)
    end do
    write(*,trim(formatspec)) 'ALL',(sum(my_inumheader(:,jstep)), jstep = 1, nStepObs+1)
    write(*,*) '----------------------------------------------------------------'
    write(*,*) 'Distribution of assimilated data over stepobs ON LOCAL PROCESSOR'
    write(*,trim(formatspec2)) 'Bin#', (jstep, jstep = 1, nStepObs),'Total'
    do jfamily = 1, nfamily
       write(*,trim(formatspec)) familylist(jfamily), (my_idataass(jfamily,jstep), &
             jstep = 1, nStepObs+1)
    end do
    write(*,trim(formatspec)) 'ALL',(sum(my_idataass(:,jstep)),jstep=1,nStepObs+1)
    write(*,*) '----------------------------------------------------------------'

    nsize = size(inumheader)
    call rpn_comm_allreduce(my_inumheader, inumheader, nsize, &
         "mpi_integer", "mpi_sum", "GRID", ierr)
    deallocate(my_inumheader) 
    nsize = size(idataass)
    call rpn_comm_allreduce(my_idataass, idataass, nsize, &
         "mpi_integer", "mpi_sum", "GRID", ierr)
    deallocate(my_idataass) 
    if(mpi_myid == 0) then
       write(*,*) '----------------------------------------------------------------'
       write(*,*) 'Distribution of number of headers over stepobs ON ALL PROCESSORS'
       write(*,trim(formatspec2)) 'Bin#', (jstep, jstep = 1, nStepObs), 'Total'
       do jfamily = 1, nfamily
          write(*,trim(formatspec)) familylist(jfamily), (inumheader(jfamily,jstep), &
                jstep = 1, nStepObs+1)
       end do
       write(*,trim(formatspec)) 'ALL', (sum(inumheader(:,jstep)), jstep = 1, nStepObs+1)
       write(*,*) '---------------------------------------------------------------'
       write(*,*) 'Distribution of assimilated data over stepobs ON ALL PROCESSORS'
       write(*,trim(formatspec2)) 'Bin#', (jstep, jstep = 1, nStepObs), 'Total'
       do jfamily = 1, nfamily
          write(*,trim(formatspec)) familylist(jfamily), (idataass(jfamily,jstep), &
                jstep = 1, nStepObs+1)
       end do
       write(*,trim(formatspec)) 'ALL', (sum(idataass(:,jstep)), jstep = 1, nStepObs+1)
       write(*,*) '---------------------------------------------------------------'
    end if

    deallocate(idataass)
    deallocate(inumheader)

  end subroutine tim_timeBinning


  subroutine tim_sutimeinterp(obsSpaceData)
    implicit none

    type(struct_obs) :: obsSpaceData

    integer :: jobs,jstep
    real(8) :: stepObsIndex

    if (.not. initialized) call utl_abort('tim_suTimeInterp: module not initialized')

    if(mpi_myid == 0) write(*,*) ' '
    if(mpi_myid == 0) write(*,*) '-------- ENTERING TIM_SUTIMEINTERP ---------'
    if(mpi_myid == 0) write(*,*) ' '

    ! Compute the number of step obs over the observation time window 
    if(mpi_myid == 0) write(*,*) 'TIM_SUTIMEINTERP: Number of step obs inc : ',tim_nstepobsinc

    if(associated(timeInterpWeight)) deallocate(timeInterpWeight)
    allocate(timeInterpWeight(obs_numHeader(obsSpaceData),mxstepobs))
    timeInterpWeight(:,:) = 0.0d0

    do jobs = 1, obs_numHeader(obsSpaceData)
      ! return the step stamp associated with date and time of the observation

      ! building the list of step stamp and counting number of obs in each step
      if(tim_nstepobsinc == 1) then
        call tim_setTimeInterpWeight(1.0d0,jobs,1)
      else
        call tim_getStepObsIndex(stepObsIndex,tim_getDatestamp(),  &
                             obs_headElem_i(obsSpaceData,OBS_DAT,jobs),  &
                             obs_headElem_i(obsSpaceData,OBS_ETM,jobs),tim_nstepobsinc)
        if(floor(stepObsIndex) >= tim_nstepobsinc) then
          write(*,*) 'tim_sutimeinterp: stepObsIndex too big=', jobs, stepObsIndex
          call tim_setTimeInterpWeight(1.0d0, jobs, tim_nstepobsinc)
        else if(floor(stepObsIndex) <= 0) then
          write(*,*) 'tim_sutimeinterp: stepObsIndex too small=',jobs, stepObsIndex
          call tim_setTimeInterpWeight(1.0d0, jobs, 1)
        else
          call tim_setTimeInterpWeight(1.0d0-(stepObsIndex-floor(stepObsIndex)), jobs, floor(stepObsIndex))
          call tim_setTimeInterpWeight(stepObsIndex-floor(stepObsIndex), jobs, floor(stepObsIndex)+1)
        end if
      end if

    end do

    if(mpi_myid == 0) write(*,*) ' '
    if(mpi_myid == 0) write(*,*) '-------- END OF TIM_SUTIMEINTERP ---------'
    if(mpi_myid == 0) write(*,*) ' '

  end subroutine tim_sutimeinterp


  subroutine tim_setTimeInterpWeight(weight_in,headerIndex,stepObs)
    implicit none
    integer, intent(in)      :: headerIndex, stepObs
    real(kind=8), intent(in) :: weight_in

    timeInterpWeight(headerIndex, stepObs) = weight_in

  end SUBROUTINE tim_setTimeInterpWeight


  function tim_getTimeInterpWeight(headerIndex,stepObs) result(weight_out)
    implicit none
    real(kind=8)        :: weight_out
    integer, intent(in) :: headerIndex, stepObs

    weight_out = timeInterpWeight(headerIndex, stepObs)

  end function tim_getTimeInterpWeight


  subroutine tim_setDatestamp(datestamp_in)
    !
    ! object: to control access to the minimization object.  Sets the date
    !         of the window centre of analysis validity to the indicated value.
    implicit none
    integer, intent(in) :: datestamp_in

    if (.not. initialized) call utl_abort('tim_setDateStamp: module not initialized')

    datestamp = datestamp_in

  end subroutine tim_setDatestamp


  function tim_getDatestamp() result(datestamp_out)
    !
    ! object: to control access to the minimization object.  Returns the date
    !         of the window centre of analysis validity.
    implicit none
    integer :: datestamp_out

    if (.not. initialized) call utl_abort('tim_getDateStamp: module not initialized')

    datestamp_out = datestamp

  end function tim_getDatestamp


  subroutine tim_getStampList(datestamplist,kstep,ksystamp)
    implicit none
    !
    ! Author: Simon Pellerin *ARMA/SMC Nov. 2001
    ! Purpose: Compute a list of STAMPS corresponding to stepobs time
    !
    ! Revisions:
    !           Y. Rochon ARQI/AQRD Julu 2016
    !           - Allowance for obs time windows of size different from 6 hours
    !             through use of tim_windowsize.
    !
    ! Dummies
    integer, intent(in) :: kstep ! number of step obs
    integer, intent(in) :: ksystamp ! Synoptic time
    integer, intent(out), dimension(kstep) :: datestamplist

    ! Locals
    integer :: jstep
    real(8) :: dldelt
    real(8) :: dtstep ! Del Time in hours between step obs

    if (.not. initialized) call utl_abort('tim_getStampList: module not initialized')

    if(kstep > 1) then
      dtstep = tim_windowsize/(real(kstep-1,8))
    else
      dtstep = tim_windowsize
    end if

    do jstep = 1,kstep
      dldelt = (jstep - ((kstep-1)/2+1)) * dtstep
      call incdatr(datestamplist(jstep),ksystamp, dldelt)
    end do

  end subroutine tim_getStampList


  subroutine tim_getStepObsIndex(dnstepobs,kstsyn,kdobs,ktobs,knstepobs)
    implicit none
    ! Author : Mark Buehner (based on stepobs.ftn90)
    !
    ! Purpose: Return the stepobs index as a real number (-1.0 if out of range)
    !
    ! Revisions:
    !           Y. Rochon ARQI/AQRD Julu 2016
    !           - Allowance for obs time windows of size different from 6 hours
    !             through use of tim_windowsize.
    !
    integer, intent(in) :: kstsyn    ! Synop CMC date-time stamp
    integer, intent(in) :: kdobs     ! Obs date YYYYMMDD
    integer, intent(in) :: ktobs     ! Obs time HHMM
    integer, intent(in) :: knstepobs ! number of stepobs in assimilation window

    real(8), intent(out):: dnstepobs ! number of stepobs from reference time

    ! Locals
    real(8) :: dddt      ! dstepobs
    integer :: newdate, istat, imode
    integer :: istobs    ! Obs CMC date-time stamp
    integer :: itobs     ! Obs time HHMMSShh
    real(8) :: dlhours   ! Del time from synop time

    if (.not. initialized) call utl_abort('tim_getStepObsIndex: module not initialized')

    ! Building observation stamp
    imode = 3 ! printable to stamp
    itobs = ktobs * 10000
    istat = newdate(istobs, kdobs, itobs, imode)

    ! Difference (in hours) between obs time and reference time
    call difdatr(istobs, kstsyn, dlhours)

    if(knstepobs > 1) then
      ! FGAT: more than 1 trial field in time window
      dddt = tim_windowsize/(real(knstepobs-1,8))
      dnstepobs = dlhours/dddt ! number of step obs from reference (e.g. synoptic)
      dnstepobs = dnstepobs + real((knstepobs+1)/2,8)
      if(dnstepobs < 0.5d0.or.dnstepobs > (0.5d0+real(knstepobs,8))) dnstepobs = -1.0d0

    else
      ! 3D-Var: only 1 trial field in time window
      if(dlhours < -tim_windowsize/2.0D0.or.dlhours > tim_windowsize/2.0D0) then
        ! outside time window
        dnstepobs = -1.0d0
      else
        ! inside time window
        dnstepobs = 1.0d0
      end if

    end if

  end subroutine tim_getStepObsIndex

end MODULE timeCoord_mod

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

!--------------------------------------------------------------------------
!!
!! *Purpose*: Main program for replacing functionality of Jeff's original 
!!            AddAnalInc (AAI) program
!!
!--------------------------------------------------------------------------
program main_aai
  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use gridStateVector_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use utilities_mod
  use ramDisk_mod
  implicit none

  type(struct_gsv) :: statevector_trial, statevector_increment,  &
                      statevector_Psfc, statevector_analysis

  type(struct_vco), pointer :: vco_trl => null()
  type(struct_hco), pointer :: hco_trl => null()

  integer              :: fclos, fnom, fstopc, ierr
  integer              :: stepIndex, numStep, middleStep
  integer              :: nulnam, dateStamp
  integer, allocatable :: dateStampList(:)
  integer              :: get_max_rss

  character(len=256)  :: trialFileName, incFileName, anlFileName
  character(len=4)    :: coffset

  real(8)             :: deltaHours
  real(8), pointer    :: PsfcTrial(:,:,:,:), PsfcAnalysis(:,:,:,:),  &
                         PsfcIncrement(:,:,:,:)
  real(8), pointer    :: GZsfc_increment(:,:), GZsfc_trial(:,:)

  ! namelist variables
  integer  :: writeNumBits
  logical  :: writeHiresIncrement
  character(len=12) :: etiket_rehm
  character(len=12) :: etiket_anlm
  NAMELIST /NAMAAI/writeHiresIncrement, etiket_rehm, etiket_anlm, writeNumBits

  write(*,'(/,' //  &
        '3(" *****************"),/,' //                   &
        '14x,"-- START OF MAIN_AAI             --",/,' //   &
        '14x,"-- Program for interpolating and adding analysis increment to background state --",/, ' //  &
        '14x,"-- Revision number ",a," --",/,' //  &
        '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  !
  !- MPI, TMG initialization
  !
  call mpi_initialize
  call tmg_init(mpi_myid, 'TMG_AAI' )

  call tmg_start(1,'MAIN')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  !
  !- Set/Read values for the namelist NAMAAI
  !
  
  !- Setting default values
  writeHiresIncrement = .true.
  etiket_rehm = 'INCREMENT'
  etiket_anlm = 'ANALYSIS'
  writeNumBits = 16

  !- Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namaai, iostat=ierr)
  if ( ierr /= 0) call utl_abort('main_aai: Error reading namelist')
  if ( mpi_myid == 0 ) write(*,nml=namaai)
  ierr = fclos(nulnam)

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  ! Setup timeCoord module (date read in NAMTIME namelist)
  call tim_setup
  numStep = tim_nstepobsinc
  allocate(dateStampList(numStep))
  call tim_getstamplist(dateStampList,numStep,tim_getDatestamp())

  !- Initialize variables of the model states
  call gsv_setup

  !- Initialize the trial state grid
  if (mpi_myid == 0) write(*,*) ''
  if (mpi_myid == 0) write(*,*) ' main_aai: Set hco parameters for trial grid'
  trialFileName = './trlm_01'
  call hco_SetupFromFile( hco_trl, trim(trialFileName), ' ')
  call vco_setupFromFile( vco_trl, trim(trialFileName) )

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  ! Read trial files, keeping HU as is (no conversion to LQ)
  !
  if(mpi_myid == 0) write(*,*) ''
  if(mpi_myid == 0) write(*,*) 'main_aai: reading background state for all time steps'
  call gsv_readTrials(hco_trl,vco_trl,statevector_trial,HUcontainsLQ_opt=.false.)

  !
  ! Read increment files only for Psfc
  !
  call gsv_allocate(statevector_Psfc, numStep, hco_trl, vco_trl, &
                    dateStamp=tim_getDateStamp(), mpi_local=.true.,  &
                    varName='P0', allocGZsfc=.true.)
  do stepIndex = 1, numStep
    dateStamp = datestamplist(stepIndex)
    if(mpi_myid == 0) write(*,*) ''
    if(mpi_myid == 0) write(*,*) 'main_aai: reading Psfc increment for time step: ',stepIndex, dateStamp

    call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
    if(nint(deltaHours*60.0d0).lt.0) then
      write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
    else
      write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
    endif
    incFileName = './rebm_' // trim(coffset) // 'm'

    PsfcTrial => gsv_getField_r8(statevector_trial,'P0')
    call gsv_readFromFile(statevector_Psfc, trim(incFileName), ' ', ' ', stepIndex,  &
                          HUcontainsLQ=.false., PsfcReference_opt=PsfcTrial(:,:,1,stepIndex))
  end do

  !
  ! Compute analysis Psfc to use for interpolation of increment
  !
  PsfcTrial => gsv_getField_r8(statevector_trial,'P0')
  PsfcIncrement => gsv_getField_r8(statevector_Psfc,'P0')
  PsfcAnalysis => gsv_getField_r8(statevector_Psfc,'P0')
  PsfcAnalysis(:,:,1,:) = PsfcTrial(:,:,1,:) + PsfcIncrement(:,:,1,:)

  !
  ! Copy the surface GZ from trial into statevector_Psfc
  !
  GZsfc_increment => gsv_getGZsfc(statevector_Psfc)
  GZsfc_trial => gsv_getGZsfc(statevector_trial)
  GZsfc_increment(:,:) = GZsfc_trial(:,:)

  !
  ! Read increment files (interpolating to trial grid/levels), also keeping HU as is
  !
  call gsv_allocate(statevector_increment, numStep, hco_trl, vco_trl, &
                    dateStamp=tim_getDateStamp(), mpi_local=.true., allocGZsfc=.true.)
  do stepIndex = 1, numStep
    dateStamp = datestamplist(stepIndex)
    if(mpi_myid == 0) write(*,*) ''
    if(mpi_myid == 0) write(*,*) 'main_aai: reading increment for time step: ',stepIndex, dateStamp

    call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
    if(nint(deltaHours*60.0d0).lt.0) then
      write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
    else
      write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
    endif
    incFileName = './rebm_' // trim(coffset) // 'm'

    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    call gsv_readFromFile(statevector_increment, trim(incFileName), ' ', ' ', stepIndex,  &
                          HUcontainsLQ=.false., PsfcReference_opt=PsfcAnalysis(:,:,1,stepIndex))
  end do

  !
  ! Write interpolated increment files
  !
  if( writeHiresIncrement ) then
    do stepIndex = 1, numStep
      dateStamp = datestamplist(stepIndex)
      if(mpi_myid == 0) write(*,*) ''
      if(mpi_myid == 0) write(*,*) 'main_aai: writing interpolated increment for time step: ',stepIndex, dateStamp

      call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
      if(nint(deltaHours*60.0d0).lt.0) then
        write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
      else
        write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
      endif
      incFileName = './rehm_' // trim(coffset) // 'm'

      write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
      call gsv_writeToFile(statevector_increment, trim(incFileName), etiket_rehm,  &
                           indexStep_in=stepIndex, typvar_in='R', numBits_opt=writeNumBits )

      ! Also write analysis value of Psfc and surface GZ to increment file
      call gsv_writeToFile(statevector_Psfc, trim(incFileName), etiket_rehm,  &
                           indexStep_in=stepIndex, typvar_in='A', writeGZsfc_opt=.true., &
                           numBits_opt=writeNumBits )
    end do
  end if
  
  !
  ! Calculate analysis state
  !
  call gsv_allocate(statevector_analysis, numStep, hco_trl, vco_trl, &
                    dateStamp=tim_getDateStamp(), mpi_local=.true., allocGZsfc=.true.)
  call gsv_copy(statevector_trial, statevector_analysis)
  call gsv_add(statevector_increment, statevector_analysis)

  !
  ! Write analysis state to file only at the central time
  !
  do stepIndex = 1, numStep
    dateStamp = datestamplist(stepIndex)
    if(mpi_myid == 0) write(*,*) ''
    if(mpi_myid == 0) write(*,*) 'main_aai: writing analysis for time step: ',stepIndex, dateStamp

    call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
    if(nint(deltaHours*60.0d0).lt.0) then
      write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
    else
      write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
    endif
    anlFileName = './anlm_' // trim(coffset) // 'm'

    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    call gsv_writeToFile(statevector_analysis, trim(anlFileName), etiket_anlm, indexStep_in=stepIndex, &
                         typvar_in='A', writeGZsfc_opt=.true., numBits_opt=writeNumBits )
  end do

  !
  !- MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_AAI' )
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- Ending
  !
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'
  if ( mpi_myid == 0 ) write(*,*) ' MAIN_AAI ENDS'
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'

end program main_aai

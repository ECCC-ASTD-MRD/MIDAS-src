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
program midas_addIncrement
  !
  ! **Purpose**: Main program for replacing functionality of Jeff's original 
  ! AddAnalInc (AAI) program
  !
  use mpi_mod
  use utilities_mod
  use ramDisk_mod
  use timeCoord_mod
  use gridStateVector_mod
  use addIncrement_mod
  implicit none

  integer :: get_max_rss
  integer :: ierr, fstopc
  character(len=256) :: trialFileName

  write(*,'(/,' //  &
        '3(" *****************"),/,' //                   &
        '14x,"-- START OF MIDAS-ADDINCREMENT         --",/,' //   &
        '14x,"-- Program for interpolating and adding analysis increment to background state --",/, ' //  &
        '14x,"-- Revision number ",a," --",/,' //  &
        '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  !
  !- Setup
  !

  !- MPI, TMG initialization
  call mpi_initialize
  call tmg_init(mpi_myid, 'TMG_ADDINC' )
  call tmg_start(1,'MAIN')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !- Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)

  !- Setup the ramdisk directory (if supplied)
  call ram_setup

  !- Setup timeCoord module (date read from trial file)
  trialFileName = './trlm_01'
  call tim_setup( fileNameForDate_opt=trim(trialFileName) )

  !- Initialize variables of the model states
  call gsv_setup

  !
  !- Interpolate the increment, add the increment and output the analysis
  !
  call adx_computeAndWriteAnalysis()

  !
  !- MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_ADDINC' )
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- Ending
  !
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'
  if ( mpi_myid == 0 ) write(*,*) ' MIDAS-ADDINCREMENT ENDS'
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'

end program midas_addIncrement

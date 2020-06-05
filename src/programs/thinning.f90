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

program midas_thinning
  !
  ! :Purpose: Thinning1 program to reduce the number of observation data
  !
  ! :Method:  Set bit 11 of *flag* according to an observation-type-specific
  !           algorithm.  Then delete all SQL records where bit 11 is set.
  !           So far, only aladin wind data are treated.
  !
  ! :Note:    In order for the SQL file size to be reduced, a script must
  !           subsequently execute the SQL command, *vacuum*.
  !
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use mpivar_mod
  use obsSpaceData_mod
  use obsFiles_mod
  use innovation_mod
  use thinning_mod
  use fSQLite

  implicit none

  integer :: istamp,exdb,exfin
  integer :: ierr, dateStamp
  type(struct_obs)    :: obsSpaceData
  character(len=48) :: obsMpiStrategy, varMode
  character(len=3)  :: obsColumnMode
  integer :: itemThinningFlagBitsList(15), numberThinningFlagBitsItems
  integer :: get_max_rss

  istamp = exdb('THINNING','DEBUT','NON')

  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM MIDAS-THINNING: --",/,' //   &
            '14x,"-- OBSERVATION THINNING          --",/, ' //&
            '14x,"-- VAR Revision number   ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  ! These values will be obtained from a thn_thin* method:
  numberThinningFlagBitsItems=0
  itemThinningFlagBitsList   =0

  ! MPI initilization
  call mpi_initialize

  call tmg_init(mpi_myid, 'TMG_THIN' )
  call tmg_start(1,'MAIN')

  varMode='bgck'  ! a necessary argument for obsf_setup


  ! 1. Top level setup

  call ram_setup


  ! 2. configuration the job

  ! Do initial set up
  call tmg_start(2,'PREMIN')

  obsMpiStrategy = 'LIKESPLITFILES'

  !     
  !- Initialize observation file names
  !
  call obsf_setup( dateStamp, varMode )

  !
  !- Setup and read observations
  !
  obsColumnMode='ALL'
  call inn_setupObs(obsSpaceData, obsColumnMode, obsMpiStrategy, varMode) ! IN

  write(*,*) 'Memory Used during set-up: ',get_max_rss()/1024,'Mb'


  ! 3. Do the Thinning

  ! Set bit 11 of flag, one observation type at a time
  call thn_thinAladin(obsSpaceData)

  ! Write bit11 to the sql files (FLG must be set in namSQLUpdate)
  call obsf_writeFiles(obsSpaceData)

  ! Delete the flagged observations, and make the files smaller
  call obsf_thinFiles(obsSpaceData)

  ! 4. Job termination

  istamp = exfin('THINNING','FIN','NON')

  ! deallocate obsSpaceData
  call obs_finalize(obsSpaceData)


  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_THIN' )

  call rpn_comm_finalize(ierr)

end program midas_thinning

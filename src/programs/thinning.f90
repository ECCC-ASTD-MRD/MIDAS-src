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
  !           algorithm.  Then remove all observations from SQL and/or burp files
  !           for which bit 11 is set. So far, only hyper-spectral IR, TOVS (amsua,
  !           amsub/mhs and atms), and aladin wind data are treated.
  !
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use obsFiles_mod
  use tovs_nl_mod
  use obsFilter_mod
  use thinning_mod
  use fSQLite

  implicit none

  type(struct_obs) :: obsSpaceData
  integer :: istamp,exdb,exfin
  integer :: ierr, dateStamp
  integer :: get_max_rss

  istamp = exdb('THINNING','DEBUT','NON')

  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM MIDAS-THINNING: --",/,' //   &
            '14x,"-- OBSERVATION THINNING          --",/, ' //&
            '14x,"-- VAR Revision number   ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  ! MPI initilization
  call mpi_initialize

  call tmg_init(mpi_myid, 'TMG_THIN' )
  call tmg_start(1,'MAIN')

  ! 1. Top level setup

  call ram_setup

  ! 2. configuration the job

  !- Initialize observation file names
  call obsf_setup( dateStamp, 'thinning' )

  !- Allocate obsSpaceData
  call obs_class_initialize('ALL')
  call obs_initialize( obsSpaceData, mpi_local=.true. )

  !- Read observations
  call tmg_start(2,'READ_OBS')
  call obsf_readFiles( obsSpaceData )
  call tmg_stop(2)

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !- Initialize TOVS processing
  if (obs_famExist( obsSpaceData, 'TO' )) call tvs_setup

  !- Select the elements to "assimilate" and apply rejection flags
  call filt_suprep( obsSpaceData )

  !- Setup timeCoord module
  call tim_setup()
  call tim_setDateStamp(dateStamp)

  ! 3. Do the Thinning

  !- Set bit 11 of flag, process one observation type at a time
  call thn_thinHyper(obsSpaceData)
  call thn_thinTovs(obsSpaceData)
  call thn_thinAircraft(obsSpaceData)
  call thn_thinAladin(obsSpaceData)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !- Update obs files and remove rejected obs (bit 11) from file (obsFileClean)
  call obsf_writeFiles(obsSpaceData, obsFileClean_opt=.true.)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! 4. Job termination

  istamp = exfin('THINNING','FIN','NON')

  !- deallocate obsSpaceData
  call obs_finalize(obsSpaceData)

  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_THIN' )

  call rpn_comm_finalize(ierr)

end program midas_thinning

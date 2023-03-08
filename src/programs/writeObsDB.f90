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

program midas_writeObsDB
  !
  ! :Purpose: Small program to create and write obsDB file.
  !
  use version_mod
  use codePrecision_mod
  use ramDisk_mod
  use obsSpaceData_mod
  use obsFiles_mod
  use utilities_mod
  use mathPhysConstants_mod
  use midasMpi_mod
  use fileNames_mod

  implicit none

  integer :: dateStampFromObs, headerIndex, ierr, nulnam
  type(struct_obs),       target :: obsSpaceData
  
  integer :: get_max_rss, fnom, fclos

  call ver_printNameAndVersion('writeObsDB','write obsDB file from scratch')

  !- 1.0 mpi
  call mmpi_initialize
  if (mmpi_nprocs /= 1) then
    write(*,*) 'mmpi_nprocs=', mmpi_nprocs
    call utl_abort('midas_writeObsDB: mmpi_nprocs should be equal to 1')
  end if

  !- 1.1 timings
  call tmg_init(mmpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'Main')

  !
  !- Initialize the ram disk
  !
  call ram_setup
  
  !
  !- Initialize observation file names, but don't use datestamp
  !
  call obsf_setup(dateStampFromObs, 'bgck')

  !
  !- Initialize constants
  !
  if (mmpi_myid == 0) then
    call mpc_printConstants(6)
    call pre_printPrecisions
  end if

  !
  !- Setup and read observations
  !
  !- Specify the active observation-array columns
  call obs_class_initialize('ALL')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Allocate memory for observation arrays
  call obs_initialize(obsSpaceData, mpi_local=obsf_filesSplit())
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Read the observations from files
  !
  call utl_tmg_start(11,'----ReadObsFiles')
  call obsf_readFiles( obsSpaceData )
  call utl_tmg_stop(11)

  ! 3 Write into the observation files
  write(*,*)
  write(*,*) '> midas-writeOBsDB: writing to file'
  call obsf_writeFiles(obsSpaceData)

  !
  ! 4.  Ending
  !
  write(*,*)
  write(*,*) '> midas-writeObsDB: Ending'
  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData

  call utl_tmg_stop(0)
  call tmg_terminate(mmpi_myid, 'TMG_INFO')

end program midas_writeObsDB
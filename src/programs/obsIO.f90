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
!! *Purpose*: This program is :
!!      - reading observation files
!!      - loading the information in 'struct_obs',
!!      - multiply the observation value by 'obsIOFactor' (specified in 'NAMOBSIO' namelist)
!!      - update the flag to value specified by 'obsIOFlag' (specified in 'NAMOBSIO' namelist)
!!      - update the observation files with those new values.
!!
!--------------------------------------------------------------------------
program midas_obsio
  use ramDisk_mod
  use utilities_mod
  use mpiVar_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use obsFiles_mod
  implicit none

  type(struct_obs), target  :: obsSpaceData

  integer :: datestamp, get_max_rss, ierr
  integer :: fnom, fclos, nulnam, headerIndex

  ! Namelist
  character(len=48) :: obsIOMpiStrategy
  real(kind=8) :: obsIOFactor
  integer :: obsIOFlag
  NAMELIST /NAMOBSIO/ obsIOFactor, obsIOFlag, obsIOMpiStrategy

  write(*,*) " ------------------------------------------"
  write(*,*) " ---  START OF MAIN PROGRAM midas-obsIO ---"
  write(*,*) " ---  Read and update observation files ---"
  write(*,*) " ------------------------------------------"

  !
  !- 1.  Settings and module initializations
  !
  write(*,*)
  write(*,*) '> midas-obsIO: setup - START'

  !- 1.0 Namelist default values
  obsIOFactor = 1.0
  obsIOFlag   = 1
  obsIOMpiStrategy = 'LIKESPLITFILES'

  nulnam = 0
  ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam,nml=namobsio,iostat=ierr)
  if (ierr /= 0) call utl_abort('midas-obsIO: Error reading namelist')
  if (mpi_myid == 0) write(*,nml=namobsio)
  ierr = fclos(nulnam)

  !- mpi
  call mpi_initialize

  !- initialize timings
  call tmg_init(mpi_myid, 'TMG_OBSIO' )
  call tmg_start(1,'MAIN')

  !- RAM disk usage
  call ram_setup

  !- Observation file names and get datestamp
  call obsf_setup(dateStamp, 'analysis' )

  write(*,*)
  write(*,*) '> midas-obsIO: dateStamp is ', dateStamp
  write(*,*) '> midas-obsIO: setup - END'
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  call obs_class_initialize('VAR')
  call obs_initialize( obsSpaceData, mpi_local=obsf_filesSplit() )

  call tmg_start(11,'READ_OBS')
  call obsf_readFiles( obsSpaceData )
  call tmg_stop(11)

  !- Update values in 'obsSpaceData'
  write(*,*)
  write(*,*) '> midas-obsIO: multiply obs value by ', obsIOFactor
  call obsIO_applyFactor(obsSpaceData,obsIOFactor)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  write(*,*)
  write(*,*) '> midas-obsIO: update flag ', obsIOFlag
  call obsIO_updateFlag(obsSpaceData,obsIOFlag)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! 2.4 Write the results

  ! 2.4.1 Into the listings
  write(*,*)
  write(*,*) '> midas-obsIO: printing the header and body for each obs'
  do headerIndex = 1, obs_numHeader(obsSpaceData)
    call obs_prnthdr(obsSpaceData,headerIndex)
    call obs_prntbdy(obsSpaceData,headerIndex)
  end do

  ! 2.4.2 Into the observation files
  write(*,*)
  write(*,*) '> midas-obsIO: writing to file'
  call obsf_writeFiles(obsSpaceData)

  !
  !- 3.  Ending
  !
  write(*,*)
  write(*,*) '> midas-obsIO: Ending'
  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData

  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_OBSIO' )

  call rpn_comm_finalize(ierr)

contains

  subroutine obsIO_updateFlag(obsSpaceData,flag)
    implicit none
    ! arguments
    type(struct_obs) :: obsSpaceData
    integer :: flag

    ! locals
    integer :: bodyIndex

    do bodyIndex=1,obs_numbody(obsSpaceData)
      call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex,flag)
      call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,1)
    end do

  end subroutine obsIO_updateFlag

  subroutine obsIO_applyFactor(obsSpaceData,factor)
    implicit none
    ! arguments
    type(struct_obs) :: obsSpaceData
    real(kind=8) :: factor

    ! locals
    integer :: bodyIndex
    real(kind=8) :: obsValue

    do bodyIndex=1,obs_numbody(obsSpaceData)
      obsValue = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
      obsValue = obsValue*factor
      call obs_bodySet_r(obsSpaceData,OBS_VAR,bodyIndex,obsValue)
    end do
  end subroutine obsIO_applyFactor


end program midas_obsio

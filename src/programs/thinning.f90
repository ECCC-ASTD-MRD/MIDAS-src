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
!! *Purpose*: Thinning1 program to reduce the number of observation data
!!
!--------------------------------------------------------------------------
program midas_thinning
  use ramDisk_mod
  use utilities_mod
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

  varMode='bgckConv'  ! a necessary argument for obsf_setup


  ! 1. Top level setup

  call ram_setup


  ! 2. configuration the job

  ! Do initial set up
  call tmg_start(2,'PREMIN')

  obsMpiStrategy = 'LIKESPLITFILES'

  call thin_setup('ALL') ! obsColumnMode 


  ! 3. Do the Thinning - set bit9 of iflag
  call obsf_thinFiles(obsSpaceData)

  ! 4. Job termination

  istamp = exfin('THINNING','FIN','NON')

  ! deallocate obsSpaceData
  call obs_finalize(obsSpaceData)


  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_THIN' )

  call rpn_comm_finalize(ierr)

contains

  !--------------------------------------------------------------------------
  !! *Purpose*: Control of the preprocessing of the variational assimilation
  !!
  !! Revisions:
  !!           Y.J. Rochon, Jan 2016
  !!           - Addition of test on availability of input trial fields according
  !!             to related observation families.
  !--------------------------------------------------------------------------
  subroutine thin_setup(obsColumnMode)
    implicit none

    character (len=*) :: obsColumnMode
    integer :: datestamp

    integer :: get_max_rss

    write(*,*) ''
    write(*,*) '-----------------------------------'
    write(*,*) '-- Starting subroutine thin_setup --'
    write(*,*) '-----------------------------------'

    !     
    !- Initialize observation file names
    !
    call obsf_setup( dateStamp, varMode )

    !
    !- Setup and read observations
    !
    call inn_setupObs(obsSpaceData, obsColumnMode, obsMpiStrategy, varMode) ! IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine thin_setup

end program midas_thinning

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
!! *Purpose*: Main program for O-F computations, background check, and thinning
!!            (O-F => Observation minus Forecast)
!!
!--------------------------------------------------------------------------
program midas_obsSelection
  !
  ! *Purpose*: Main program for O-F computations, background check, and thinning
  !            (O-F => Observation minus Forecast)
  !
  use oMinusF_mod
  use backgroundCheck_mod
  use thinning_mod
  use obsSpaceData_mod
  use columnData_mod
  use obsFiles_mod
  use utilities_mod
  use mpi_mod
  implicit none

  integer :: fnom, fclos, nulnam, ierr, headerIndex
  type(struct_columnData),target :: trlColumnOnAnlLev
  type(struct_columnData),target :: trlColumnOnTrlLev
  type(struct_obs),       target :: obsSpaceData

  write(*,*) " -------------------------------------------------"
  write(*,*) " ---  START OF MAIN PROGRAM midas-obsSelection ---"
  write(*,*) " ---  Computation of the innovation            ---"
  write(*,*) " -------------------------------------------------"

  !- 1.0 mpi
  call mpi_initialize

  !- 1.1 timings
  call tmg_init(mpi_myid, 'TMG_OBSSELECTION' )
  call tmg_start(1,'MAIN')

  if ( mpi_myid == 0 ) then
    call utl_writeStatus('VAR3D_BEG')
  endif


  ! 2.1 Calculate the Observation - Forecast difference
  call omf_oMinusF(trlColumnOnAnlLev, trlColumnOnTrlLev, obsSpaceData, &
                   'bgckConv', addHBHT=.true., addSigmaO=.true.)


  ! 2.2 Perform the background check
  !     The routine also calls compute_HBHT and writes to listings & obsSpaceData
  call bgck_bgcheck_conv(trlColumnOnAnlLev, trlColumnOnTrlLev, obsSpaceData)

  !
  ! 2.3 Thinning1:  Set bit 11 of flag, one observation type at a time
  !
  call thn_thinAladin(obsSpaceData)

  ! 3 Write the results

  ! 3.1 Into the listings
  write(*,*)
  write(*,*) '> midas-obsSelection: printing the FIRST header and body'
  do headerIndex = 1, min(1,obs_numHeader(obsSpaceData))
    call obs_prnthdr(obsSpaceData,headerIndex)
    call obs_prntbdy(obsSpaceData,headerIndex)
  end do
  ! 3.2 Into the observation files
  write(*,*)
  write(*,*) '> midas-obsSelection: writing to file'
  call obsf_writeFiles(obsSpaceData)

  ! Delete the flagged observations, and make the files smaller
  call obsf_thinFiles(obsSpaceData)

  !
  ! 4.  Ending
  !
  write(*,*)
  write(*,*) '> midas-obsSelection: Ending'
  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData

  call rpn_comm_finalize(ierr)

  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_OMINUSF' )

  if ( mpi_myid == 0 ) then
    call utl_writeStatus('VAR3D_END')
  endif

end program midas_obsSelection

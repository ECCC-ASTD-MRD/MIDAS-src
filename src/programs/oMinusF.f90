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

program midas_ominusf
  !
  ! :Purpose: Main program for Observation minus Forecast (O-F) computation
  !
  use version_mod
  use oMinusF_mod
  use obsSpaceData_mod
  use columnData_mod
  use obsFiles_mod
  use utilities_mod
  use mpi_mod
  use biasCorrectionSat_mod
  implicit none

  ! Namelist
  logical :: addHBHT
  logical :: addSigmaO
  NAMELIST /NAMOMF/addHBHT, addSigmaO
  
  integer :: fnom, fclos, nulnam, ierr, headerIndex
  type(struct_columnData),target  :: columnTrlOnAnlIncLev
  type(struct_columnData),target  :: columnTrlOnTrlLev
  type(struct_obs),       target  :: obsSpaceData

  call ver_printNameAndVersion('oMinusF','Computation of the innovation')

  !- 1.0 mpi
  call mpi_initialize

  !- 1.1 timings
  call tmg_init(mpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'Main')

  if ( mpi_myid == 0 ) then
    call utl_writeStatus('VAR3D_BEG')
  endif

  !- 1.2 Namelist
  addHBHT   = .false. ! default value
  addSigmaO = .false.

  nulnam = 0
  ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam,nml=namomf,iostat=ierr)
  if (ierr /= 0) call utl_abort('midas-OminusF: Error reading namelist')
  if (mpi_myid == 0) write(*,nml=namomf)
  ierr = fclos(nulnam)

  ! 2.1 Calculate the Observation - Forecast difference
  call omf_oMinusF(columnTrlOnAnlIncLev, columnTrlOnTrlLev, obsSpaceData, &
                   'OminusF', addHBHT, addSigmaO)

  call bcs_calcBias(obsSpaceData,columnTrlOnTrlLev) ! Fill in OBS_BCOR obsSpaceData column with computed bias correction

  call bcs_applyBiasCorrection(obsSpaceData,OBS_VAR,"TO") ! Apply bias correction to OBS

  call bcs_applyBiasCorrection(obsSpaceData,OBS_OMP,"TO") ! Apply bias correction to O-F

  ! 2.3 Write the results

  ! 2.3.1 Into the listings
  write(*,*)
  write(*,*) '> midas-OminusF: printing the FIRST header and body'
  do headerIndex = 1, min(1,obs_numHeader(obsSpaceData))
    call obs_prnthdr(obsSpaceData,headerIndex)
    call obs_prntbdy(obsSpaceData,headerIndex)
  end do
  ! 2.3.2 Into the observation files
  write(*,*)
  write(*,*) '> midas-OminusF: writing to file'
  call obsf_writeFiles(obsSpaceData)

  !
  !- 3.  Ending
  !
  write(*,*)
  write(*,*) '> midas-OminusF: Ending'
  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData

  call tmg_stop(0)
  call tmg_terminate(mpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr)

  if ( mpi_myid == 0 ) then
    call utl_writeStatus('VAR3D_END')
  endif

end program midas_ominusf

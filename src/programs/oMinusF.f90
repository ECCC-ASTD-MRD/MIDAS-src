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
  use midasMpi_mod
  use biasCorrectionSat_mod
  use ensembleObservations_mod
  use fileNames_mod
  implicit none

  ! Namelist
  integer :: nEns       ! ensemble size
  logical :: addHBHT    ! choose to add the value of HBHT to obsSpaceData so it can be output
  logical :: addSigmaO  ! choose to add the value of sigma_obs to obsSpaceData so it can be output
  NAMELIST /NAMOMF/addHBHT, addSigmaO, nEns
  
  integer :: fnom, fclos, nulnam, ierr, headerIndex
  
  type(struct_columnData),target  :: columnTrlOnAnlIncLev
  type(struct_columnData),target  :: columnTrlOnTrlLev
  type(struct_obs),       target  :: obsSpaceData
  type(struct_eob),       target  :: ensObs
  
  character(len=20)  :: oMinusFmode
  character(len=256) :: ensPathName = 'ensemble'
  character(len=256) :: ensFileName
  character(len=256) :: trlFileName = './trlm_01'

  logical :: ensFileExists, trlFileExists
  
  call ver_printNameAndVersion('oMinusF','Computation of the innovation')

  !
  !- 1.  Initialization
  !
  
  !- 1.1 mpi
  call mmpi_initialize

  !- 1.2 timings
  call tmg_init(mmpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'Main')

  if ( mmpi_myid == 0 ) then
    call utl_writeStatus('VAR3D_BEG')
  endif

  !- 1.3 Namelist
  addHBHT   = .false. ! default value
  addSigmaO = .false.
  nEns      = 20

  nulnam = 0
  ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam,nml=namomf,iostat=ierr)
  if (ierr /= 0) call utl_abort('midas-OminusF: Error reading namelist')
  if (mmpi_myid == 0) write(*,nml=namomf)
  ierr = fclos(nulnam)

  !- 1.4 Set mode
  inquire(file=trim(trlFileName),exist=trlFileExists)
  call fln_ensFileName(ensFileName, ensPathName, memberIndex_opt=1, &
                       shouldExist_opt=.false.)
  inquire(file=trim(ensFileName),exist=ensFileExists)

  if      (trlFileExists) then
    write(*,*)
    write(*,*) 'Trial/Prog file found'
    write(*,*) 'Setting mode to DETERMINISTIC'
    oMinusFmode = 'deterministic'
  else if (ensFileExists) then
    write(*,*)
    write(*,*) 'Ensemble file found'
    write(*,*) 'Setting mode to ENSEMBLE'
    oMinusFmode = 'ensemble'
  else
    write(*,*)
    write(*,*) 'trlFileName = ', trim(trlFileName)
    write(*,*) 'ensFileName = ', trim(ensFileName)
    call utl_abort('oMinusF : did not find a trial/prog or ensemble file')
  end if

  !
  !- 2.  Calculate the Observation - Forecast difference
  !
  if (trim(oMinusFmode) == 'deterministic') then

    !- 2.1 Compute O-F and store in obsSpaceDate
    call omf_oMinusF(columnTrlOnAnlIncLev, columnTrlOnTrlLev, obsSpaceData, &
                     'OminusF', addHBHT, addSigmaO)

    !- 2.2 Remove biases
    call bcs_calcBias(obsSpaceData,columnTrlOnTrlLev) ! Fill in OBS_BCOR obsSpaceData column with computed bias correction

    call bcs_applyBiasCorrection(obsSpaceData,OBS_VAR,"TO") ! Apply bias correction to OBS
    call bcs_applyBiasCorrection(obsSpaceData,OBS_OMP,"TO") ! Apply bias correction to O-F

    !- 2.3 Write the results

    !  2.3.1 Into the listings
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

  else ! ensemble

    !- 2.1 Compute O-F and store in ensObs
    call omf_oMinusFens( ensObs, obsSpaceData, nEns, ensPathName, &
                         'OminusF', addHBHT, addSigmaO)

    !- 2.3 Write the results
    call obsf_writeFiles(obsSpaceData, ensObs_opt=ensObs)
    
  end if  
    
  !
  !- 3.  Ending
  !
  write(*,*)
  write(*,*) '> midas-OminusF: Ending'

  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData
  
  call utl_tmg_stop(0)
  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr)

  if ( mmpi_myid == 0 ) then
    call utl_writeStatus('VAR3D_END')
  endif

end program midas_ominusf

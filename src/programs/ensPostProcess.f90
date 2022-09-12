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

program midas_ensPostProcess
  ! :Purpose: Post-processing program for the local ensemble transform Kalman filter (LETKF).
  !           Many aspects of this program are controlled throught the namelist
  !           block namEnsPostProc defined in epp_postProcess.
  use version_mod
  use midasMpi_mod
  use fileNames_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use utilities_mod
  use ramDisk_mod
  use ensPostProcess_mod
  implicit none

  type(struct_ens)          :: ensembleTrl
  type(struct_ens)          :: ensembleAnl
  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()
  type(struct_gsv)          :: stateVectorHeightSfc, stateVectorCtrlTrl

  character(len=256) :: ensPathNameAnl = 'ensemble_anal'
  character(len=256) :: ensPathNameTrl = 'ensemble_trial'
  character(len=256) :: ensFileName, gridFileName, ctrlFileName
  integer, allocatable :: dateStampList(:)
  integer :: ierr, nulnam, stepIndex, dateStamp
  logical :: targetGridFileExists
  integer, external :: get_max_rss, fstopc, fnom, fclos

  integer :: nEns             ! ensemble size
  logical :: readTrlEnsemble  ! activate reading of trial ensemble
  logical :: readAnlEnsemble  ! activate reading of analysis ensemble
  logical :: writeTrlEnsemble ! activate writing of the trial ensemble (useful when it's interpolated)
  character(len=12) :: hInterpolationDegree ! select degree of horizontal interpolation (if needed)
  NAMELIST /namEnsPostProc/nEns, readTrlEnsemble, readAnlEnsemble, &
                           writeTrlEnsemble, hInterpolationDegree

  call ver_printNameAndVersion('ensPostProcess','Program for post-processing of LETKF analysis ensemble')

  !
  !- 0. MPI, TMG and misc. initialization
  !
  call mmpi_initialize
  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  !- Setting default namelist variable values
  nEns = 256
  readTrlEnsemble  = .true.
  readAnlEnsemble  = .true.
  writeTrlEnsemble = .false.
  hInterpolationDegree = 'LINEAR' ! or 'CUBIC' or 'NEAREST'

  !- Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namEnsPostProc, iostat=ierr)
  if ( ierr /= 0) call utl_abort('midas-ensPostProcess: Error reading namelist')
  if ( mmpi_myid == 0 ) write(*,nml=namEnsPostProc)
  ierr = fclos(nulnam)

  if (.not.readTrlEnsemble .and. .not.readAnlEnsemble) then
    call utl_abort('midas-ensPostProcess: must read either Trial or Analysis ensemble')
  end if

  if (writeTrlEnsemble .and. .not.readTrlEnsemble) then
    call utl_abort('midas-ensPostProcess: cannot write Trial ensemble if it is not read')
  end if

  if (writeTrlEnsemble .and. readAnlEnsemble) then
    call utl_abort('midas-ensPostProcess: cannot write Trial ensemble when Analysis ensemble is read')
  end if

  !- 1. Initialize date/time-related info

  ! Setup timeCoord module, set dateStamp with value from trial or analysis ensemble member 1
  if (readTrlEnsemble) then
    call fln_ensFileName(ensFileName, ensPathNameTrl, memberIndex_opt=1, &
                         copyToRamDisk_opt=.false.)
  else
    call fln_ensFileName(ensFileName, ensPathNameAnl, memberIndex_opt=1, &
                         copyToRamDisk_opt=.false.)
  end if
  call tim_setup(fileNameForDate_opt=ensFileName)
  allocate(dateStampList(tim_nstepobsinc))
  call tim_getstamplist(dateStampList,tim_nstepobsinc,tim_getDatestamp())

  dateStamp = tim_getDatestamp()
  write(*,*) 'midas-ensPostProcess: analysis dateStamp = ', dateStamp

  !- 2. Initialize variables and grids

  !- Initialize variables of the model states
  call gsv_setup

  !- Initialize the grid from targetgrid file or from trial or analysis ensemble member 1
  if (mmpi_myid == 0) write(*,*) ''
  if (mmpi_myid == 0) write(*,*) 'midas-ensPostProcess: Set hco and vco parameters for ensemble grid'
  inquire(file='targetgrid', exist=targetGridFileExists)
  if (targetGridFileExists) then
    gridFileName = 'targetgrid'
  else if (readTrlEnsemble) then
    call fln_ensFileName(gridFileName, ensPathNameTrl, memberIndex_opt=1, &
                         copyToRamDisk_opt=.false.)
  else
    call fln_ensFileName(gridFileName, ensPathNameAnl, memberIndex_opt=1, &
                         copyToRamDisk_opt=.false.)
  end if
  if (mmpi_myid == 0) then
    write(*,*) 'midas-ensPostProcess: file use to define grid = ', trim(gridFileName)
  end if
  call hco_SetupFromFile(hco_ens, gridFileName, ' ', 'ENSFILEGRID')
  call vco_setupFromFile(vco_ens, gridFileName)

  !- Read the sfc height from trial ensemble member 1 - only if we are doing NWP!
  if (vco_getNumLev(vco_ens,'TH') > 0 .or. vco_getNumLev(vco_ens,'MM') > 0) then
    if (readTrlEnsemble) then
      call fln_ensFileName(ensFileName, ensPathNameTrl, memberIndex_opt=1, &
                           copyToRamDisk_opt=.false.)
    else
      call fln_ensFileName(ensFileName, ensPathNameAnl, memberIndex_opt=1, &
                           copyToRamDisk_opt=.false.)
    end if
    call gsv_allocate(stateVectorHeightSfc, 1, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                      hInterpolateDegree_opt=hInterpolationDegree, &
                      dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','TT'/))
    call gio_readFromFile(stateVectorHeightSfc, ensFileName, ' ', ' ',  &
                          containsFullField_opt=.true., readHeightSfc_opt=.true.)
  end if

  !- 3. Allocate and read ensembles

  call utl_tmg_start(2,'--ReadEnsemble')

  !- Allocate ensembles, read the Anl ensemble
  if (readAnlEnsemble) then
    call fln_ensFileName(ensFileName, ensPathNameAnl, resetFileInfo_opt=.true.)
    call ens_allocate(ensembleAnl, nEns, tim_nstepobsinc, hco_ens, vco_ens, &
                      dateStampList, hInterpolateDegree_opt=hInterpolationDegree)
    call ens_readEnsemble(ensembleAnl, ensPathNameAnl, biPeriodic=.false.)
  end if

  !- Allocate ensembles, read the Trl ensemble
  if (readTrlEnsemble) then
    call fln_ensFileName(ensFileName, ensPathNameAnl, resetFileInfo_opt=.true.)
    call ens_allocate(ensembleTrl, nEns, tim_nstepobsinc, hco_ens, vco_ens, &
                      dateStampList, hInterpolateDegree_opt=hInterpolationDegree)
    call ens_readEnsemble(ensembleTrl, ensPathNameTrl, biPeriodic=.false.)

  end if

  call utl_tmg_stop(2)

  !- Allocate and read the Trl control member
  if (readTrlEnsemble .and. readAnlEnsemble) then
    !- Allocate and read control member Trl
    call gsv_allocate( stateVectorCtrlTrl, tim_nstepobsinc, hco_ens, vco_ens, &
                       dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       hInterpolateDegree_opt=hInterpolationDegree, &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call fln_ensFileName(ctrlFileName, ensPathNameTrl, memberIndex_opt=0, &
                         copyToRamDisk_opt=.false.)
    do stepIndex = 1, tim_nstepobsinc
      call gio_readFromFile( stateVectorCtrlTrl, ctrlFileName, ' ', ' ',  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true., &
                             readHeightSfc_opt=.false. )
    end do
  end if

  !- 4. Post processing of the analysis results (if desired) and write everything to files
  call epp_postProcess(ensembleTrl, ensembleAnl, &
                       stateVectorHeightSfc, stateVectorCtrlTrl, &
                       writeTrlEnsemble)

  !
  !- 5. MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  if ( mmpi_myid == 0 ) write(*,*) ' --------------------------------'
  if ( mmpi_myid == 0 ) write(*,*) ' MIDAS-ensPostProcess ENDS'
  if ( mmpi_myid == 0 ) write(*,*) ' --------------------------------'

end program midas_ensPostProcess

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

program midas_ensembleH
  !
  ! :Purpose: Main program for applying the observation operator to an ensemble
  !           of states as the first step for most EnKF algorithms.
  !
  use version_mod
  use midasMpi_mod
  use fileNames_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use columnData_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use utilities_mod
  use ramDisk_mod
  use statetocolumn_mod
  use obsFiles_mod
  use obsSpaceData_mod
  use obsErrors_mod
  use innovation_mod
  use ensembleObservations_mod
  use ensembleStateVector_mod
  use enkf_mod
  implicit none

  type(struct_obs), target  :: obsSpaceData
  type(struct_ens)          :: ensembleTrl4D
  type(struct_gsv)          :: stateVectorMeanTrl4D
  type(struct_gsv)          :: stateVector4D
  type(struct_gsv)          :: stateVector4Dmod
  type(struct_gsv)          :: stateVectorWithZandP4D
  type(struct_gsv)          :: stateVectorHeightSfc
  type(struct_columnData)   :: column

  type(struct_eob), target  :: ensObs, ensObs_mpiglobal
  type(struct_eob), pointer :: ensObsGain, ensObsGain_mpiglobal

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()

  integer :: get_max_rss, fclos, fnom, fstopc, ierr
  integer :: memberIndex, nulnam, dateStamp
  integer :: nEnsGain, eigenVectorIndex, memberIndexInEnsObs
  integer, allocatable :: dateStampList(:)

  logical  :: useModulatedEns

  character(len=256)  :: ensFileName
  character(len=9)    :: obsColumnMode
  character(len=48)   :: obsMpiStrategy
  character(len=48)   :: midasMode
  character(len=10)   :: obsFileType

  ! namelist variables
  character(len=256) :: ensPathName
  character(len=20)  :: obsTimeInterpType ! type of time interpolation to obs time
  integer  :: nEns
  integer  :: numRetainedEigen ! number of retained eigenValues/Vectors of vertical localization matrix
                               !   used only when generating modulated ensembles.
  integer  :: fileMemberIndex1 ! first member number in ensemble set.
  real(8)  :: vLocalize        ! vertical localization radius (units: ln(Pressure in Pa) or meters)
                               !   used only when generating modulated ensembles.
  logical  :: writeEnsObsToFile
  NAMELIST /NAMENSEMBLEH/nEns, ensPathName, obsTimeInterpType, numRetainedEigen, &
                         vLocalize, writeEnsObsToFile, fileMemberIndex1

  midasMode = 'analysis'
  obsColumnMode = 'ENKFMIDAS'
  obsMpiStrategy = 'LIKESPLITFILES'

  call ver_printNameAndVersion('ensembleH','Program for applying H to ensemble')

  !
  !- 0. MPI, TMG initialization
  !
  call mmpi_initialize
  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  ! Setting default namelist variable values
  nEns                  = 10
  ensPathName           = 'ensemble'
  obsTimeInterpType     = 'LINEAR'
  numRetainedEigen      = 0
  vLocalize             = -1.0D0
  writeEnsObsToFile     = .false.
  fileMemberIndex1      = 1

  ! Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namensembleh, iostat=ierr)
  if (ierr /= 0) call utl_abort('midas-ensembleH: Error reading namelist')
  if (mmpi_myid == 0) write(*,nml=namensembleh)
  ierr = fclos(nulnam)
  
  if (numRetainedEigen < 0) call utl_abort('midas-ensembleH: numRetainedEigen should be ' // &
                                             'equal or greater than zero')

  useModulatedEns = (numRetainedEigen > 0)
  if (useModulatedEns .and. vLocalize <= 0) then
    call utl_abort('midas-ensembleH: vLocalize should be greater than zero for modulated ens')
  end if

  ! Read the observations
  call obsf_setup(dateStamp, midasMode, obsFileType_opt = obsFileType)
  if (obsFileType /= 'BURP' .and. obsFileType /= 'SQLITE') then
    call utl_abort('midas-ensembleH: only BURP and SQLITE are valid obs file formats')
  end if

  ! Use the first ensemble member to initialize datestamp and grid
  call fln_ensFileName(ensFileName, ensPathName, memberIndex_opt=1, &
                       fileMemberIndex1_opt=fileMemberIndex1)

  ! Setup timeCoord module, get datestamp from ensemble member
  call tim_setup(fileNameForDate_opt = ensFileName)
  allocate(dateStampList(tim_nstepobs))
  call tim_getstamplist(dateStampList,tim_nstepobs,tim_getDatestamp())
  
  write(*,*) 'midas-ensembleH: dateStamp = ',dateStamp

  !- Initialize variables of the model states
  call gsv_setup

  !- Initialize the Ensemble grid
  if (mmpi_myid == 0) write(*,*) ''
  if (mmpi_myid == 0) write(*,*) 'midas-ensembleH: Set hco and vco parameters for ensemble grid'
  call hco_SetupFromFile(hco_ens, ensFileName, ' ', 'ENSFILEGRID')
  call vco_setupFromFile(vco_ens, ensFileName)

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  ! read in the observations
  call inn_setupObs(obsSpaceData, hco_ens, obsColumnMode, obsMpiStrategy, midasMode)

  ! Initialize obs error covariances
  call oer_setObsErrors(obsSpaceData, midasMode)

  call col_setup
  call col_setVco(column, vco_ens)
  call col_allocate(column, obs_numheader(obsSpaceData))
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  ! Initialize the observation error covariances
  call oer_setObsErrors(obsSpaceData, midasMode) ! IN

  ! Allocate and initialize eob object for storing HX values
  call eob_allocate(ensObs, nEns, obs_numBody(obsSpaceData), obsSpaceData, &
                    fileMemberIndex1_opt=fileMemberIndex1)
  call eob_zero(ensObs)
  if (useModulatedEns) then
    nEnsGain = nEns * numRetainedEigen
    allocate(ensObsGain)
    call eob_allocate(ensObsGain, nEnsGain, obs_numBody(obsSpaceData), obsSpaceData, &
                      fileMemberIndex1_opt=fileMemberIndex1)
    call eob_zero(ensObsGain)
  else
    ensObsGain => ensObs
  end if
  
  ! Set lat, lon, obs values in ensObs
  call eob_setLatLonObs(ensObs)
  if (useModulatedEns) call eob_setLatLonObs(ensObsGain)


  ! Read the sfc height from ensemble member 1
  call gsv_allocate(stateVectorHeightSfc, 1, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                    mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                    dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','TT'/))
  call gio_readFromFile(stateVectorHeightSfc, ensFileName, ' ', ' ',  &
                        containsFullField_opt=.true., readHeightSfc_opt=.true.)

  ! Allocate statevector related to ensemble mean
  call gsv_allocate(stateVectorMeanTrl4D, tim_nstepobs, hco_ens, vco_ens, &
                    dateStamp_opt=tim_getDateStamp(),  &
                    mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                    dataKind_opt=4, allocHeightSfc_opt=.true., &
                    allocHeight_opt=.false., allocPressure_opt=.false.)
  call gsv_zero(stateVectorMeanTrl4D)
  call gsv_allocate(stateVector4D, tim_nstepobs, hco_ens, vco_ens, &
                    dateStamp_opt=tim_getDateStamp(),  &
                    mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                    dataKind_opt=4, allocHeightSfc_opt=.true., &
                    allocHeight_opt=.false., allocPressure_opt=.false.)
  call gsv_zero(stateVector4D)
  if (useModulatedEns) then
    ! same as stateVector4D
    call gsv_allocate(stateVector4Dmod, tim_nstepobs, hco_ens, vco_ens, &
                      dateStamp_opt=tim_getDateStamp(),  &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                      dataKind_opt=4, allocHeightSfc_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_zero(stateVector4Dmod)
  end if
  
  ! Allocate statevector for storing state with heights and pressures allocated (for s2c_nl)
  call gsv_allocate(stateVectorWithZandP4D, tim_nstepobs, hco_ens, vco_ens, &
                    dateStamp_opt=tim_getDateStamp(),  &
                    mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                    dataKind_opt=4, allocHeightSfc_opt=.true.)
  call gsv_zero(stateVectorWithZandP4D)

  ! Allocate ensembles, read the Trl ensemble
  call utl_tmg_start(2,'--ReadEnsemble')
  call ens_allocate(ensembleTrl4D, nEns, tim_nstepobs, hco_ens, vco_ens, dateStampList, &
                    fileMemberIndex1_opt=fileMemberIndex1)
  call ens_readEnsemble(ensembleTrl4D, ensPathName, biPeriodic=.false.)
  call utl_tmg_stop(2)

  ! Compute ensemble mean and copy to meanTrl stateVectors
  call ens_computeMean(ensembleTrl4D)
  call ens_copyEnsMean(ensembleTrl4D, stateVectorMeanTrl4D)

  do memberIndex = 1, nEns

    write(*,*) ''
    write(*,*) 'midas-letkf: apply nonlinear H to ensemble member ', memberIndex
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    ! copy 1 member to a stateVector
    call ens_copyMember(ensembleTrl4D, stateVector4D, memberIndex)
    call gsv_copy(stateVector4D, stateVectorWithZandP4D, allowVarMismatch_opt=.true., &
                  beSilent_opt=.true.)
    call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)

    ! Compute and set Yb in ensObs
    call s2c_nl(stateVectorWithZandP4D, obsSpaceData, column, hco_ens, &
                timeInterpType=obsTimeInterpType, dealloc_opt=.false., &
                beSilent_opt=.true.)

    ! Compute Y-H(X) in OBS_OMP
    call inn_computeInnovation(column, obsSpaceData, beSilent_opt=.true.)

    ! Copy to ensObs: Y-HX for this member
    call eob_setYb(ensObs, memberIndex)

    ! Compute and set Yb in ensObsGain
    do eigenVectorIndex = 1, numRetainedEigen
      if (mmpi_myid == 0) write(*,*) 'midas-ensembleH: apply nonlinear H to modulated member ', &
                                        eigenVectorIndex, '/', numRetainedEigen

      ! modulate the member with eigenvectors of vertical localization matrix
      call enkf_getModulatedState(stateVector4D, stateVectorMeanTrl4D, &
                                  vLocalize, numRetainedEigen, nEns, &
                                  eigenVectorIndex, stateVector4Dmod, &
                                  beSilent=.true.)

      call gsv_copy(stateVector4Dmod, stateVectorWithZandP4D, allowVarMismatch_opt=.true., &
                    beSilent_opt=.true.)
      call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)

      call s2c_nl(stateVectorWithZandP4D, obsSpaceData, column, hco_ens, &
                  timeInterpType=obsTimeInterpType, dealloc_opt=.false., &
                  beSilent_opt=.true.)

      ! Compute Y-H(X) in OBS_OMP
      call inn_computeInnovation(column, obsSpaceData, filterObsAndInitOer_opt=.false., &
                                 beSilent_opt=.true.)

      ! Copy to ensObsGain: Y-HX for this member
      memberIndexInEnsObs = (eigenVectorIndex - 1) * nEns + memberIndex
      call eob_setYb(ensObsGain, memberIndexInEnsObs)
    end do ! eigenVectorIndex
    
  end do
  call gsv_deallocate(stateVectorWithZandP4D)
  if (gsv_isAllocated(stateVector4Dmod)) call gsv_deallocate(stateVector4Dmod)
  call gsv_deallocate(stateVector4D)
  call gsv_deallocate(stateVectorMeanTrl4D)
  call gsv_deallocate(stateVectorHeightSfc)

  ! write local ensObs to file
  if ( writeEnsObsToFile ) call eob_writeToFilesMpiLocal(ensObs)

  ! Clean and globally communicate obs-related data, then write to files
  call eob_allGather(ensObs,ensObs_mpiglobal)
  call eob_writeToFiles(ensObs_mpiglobal, outputFilenamePrefix='eob_HX', writeObsInfo=.true.)
  if (useModulatedEns) then
    allocate(ensObsGain_mpiglobal)
    call eob_allGather(ensObsGain, ensObsGain_mpiglobal)
    call eob_writeToFiles(ensObsGain_mpiglobal, outputFilenamePrefix='eobGain_HX', writeObsInfo=.false.)
  end if

  !
  !- MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 7.  Ending
  !
  if (mmpi_myid == 0) write(*,*) ' --------------------------------'
  if (mmpi_myid == 0) write(*,*) ' MIDAS-ENSEMBLEH ENDS'
  if (mmpi_myid == 0) write(*,*) ' --------------------------------'

end program midas_ensembleH

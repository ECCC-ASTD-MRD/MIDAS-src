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
  !           block namEnsPostProc defined in enkf_postProcess.
  use mpi_mod
  use fileNames_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use analysisGrid_mod
  use timeCoord_mod
  use utilities_mod
  use ramDisk_mod
  use enkf_mod
  implicit none

  type(struct_ens), pointer :: ensembleTrl
  type(struct_ens) :: ensembleAnl
  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()
  type(struct_hco), pointer :: hco_ens_core => null()
  type(struct_gsv)          :: stateVectorHeightSfc

  character(len=256) :: ensPathNameAnl = 'ensemble_anal'
  character(len=256) :: ensPathNameTrl = 'ensemble_trial'
  character(len=256) :: ensFileName, ensFileBaseName
  integer, allocatable :: dateStampList(:)
  integer :: ierr, nulnam
  integer, external :: get_max_rss, fstopc, fnom, fclos

  integer :: nEns ! ensemble size
  NAMELIST /NAMLETKF/nEns

  write(*,'(/,' //  &
        '3(" *****************"),/,' //                   &
        '14x,"-- START OF MIDAS-ensPostProcess               --",/,' //   &
        '14x,"-- Program for post-processing of LETKF analysis ensemble --",/, ' //  &
        '14x,"-- Revision number ",a," --",/,' //  &
        '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  !
  !- 0. MPI, TMG and misc. initialization
  !
  call mpi_initialize
  call tmg_init(mpi_myid, 'TMG_LETKF')

  call tmg_start(1,'MAIN')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  !- Setting default namelist variable values
  nEns = 256
  !- Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namletkf, iostat=ierr)
  if ( ierr /= 0) call utl_abort('midas-ensPostProcess: Error reading namelist')
  if ( mpi_myid == 0 ) write(*,nml=namletkf)
  ierr = fclos(nulnam)

  !- 1. Initialize date/time-related info

  ! Setup timeCoord module, set dateStamp with value from trial ensemble member 1
  call fln_ensFileName(ensFileName, ensPathNameTrl, memberIndex_opt=1, ensFileBaseName_opt=ensFileBaseName)
  call tim_setup(fileNameForDate_opt=ensFileName)
  allocate(dateStampList(tim_nstepobsinc))
  call tim_getstamplist(dateStampList,tim_nstepobsinc,tim_getDatestamp())

  write(*,*) 'midas-ensPostProcess: analysis dateStamp = ',tim_getDatestamp()

  !- 2. Initialize variables and grids

  !- Initialize variables of the model states
  call gsv_setup

  !- Initialize the grid from trial ensemble member 1
  if (mpi_myid == 0) write(*,*) ''
  if (mpi_myid == 0) write(*,*) 'midas-ensPostProcess: Set hco and vco parameters for ensemble grid'
  call fln_ensFileName(ensFileName, ensPathNameTrl, memberIndex_opt=1)
  call hco_SetupFromFile(hco_ens, ensFileName, ' ', 'ENSFILEGRID')
  call vco_setupFromFile(vco_ens, ensFileName)
  if (vco_getNumLev(vco_ens, 'MM') /= vco_getNumLev(vco_ens, 'TH')) then
    call utl_abort('midas-letkf: nLev_M /= nLev_T - currently not supported')
  end if

  if ( hco_ens % global ) then
    call agd_SetupFromHCO( hco_ens ) ! IN
  else
    ! Setup the LAM analysis grid metrics
    hco_ens_core => hco_ens
    call agd_SetupFromHCO( hco_ens, hco_ens_core ) ! IN
  end if

  !- Read the sfc height from trial ensemble member 1
  call gsv_allocate(stateVectorHeightSfc, 1, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                    mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                    dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','TT'/))
  call gsv_readFromFile(stateVectorHeightSfc, ensFileName, ' ', ' ',  &
                        containsFullField_opt=.true., readHeightSfc_opt=.true.)

  !- 3. Allocate and read ensembles

  !- Allocate ensembles, read the Anl ensemble
  call fln_ensFileName(ensFileName, ensPathNameAnl, resetFileInfo_opt=.true.)
  call ens_allocate(ensembleAnl, nEns, tim_nstepobsinc, hco_ens, vco_ens, dateStampList)
  call ens_readEnsemble(ensembleAnl, ensPathNameAnl, biPeriodic=.false.)

  !- Allocate ensembles, read the Trl ensemble
  call fln_ensFileName(ensFileName, ensPathNameAnl, resetFileInfo_opt=.true.)
  allocate(ensembleTrl)
  call ens_allocate(ensembleTrl, nEns, tim_nstepobsinc, hco_ens, vco_ens, dateStampList)
  call ens_readEnsemble(ensembleTrl, ensPathNameTrl, biPeriodic=.false.)

  !- 4. Post processing of the analysis results (if desired) and write everything to files
  call tmg_start(8,'LETKF-postProcess')
  call enkf_postProcess(ensembleAnl, ensembleTrl, stateVectorHeightSfc)
  call tmg_stop(8)

  !
  !- 5. MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_LETKF')
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'
  if ( mpi_myid == 0 ) write(*,*) ' MIDAS-ensPostProcess ENDS'
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'

end program midas_ensPostProcess

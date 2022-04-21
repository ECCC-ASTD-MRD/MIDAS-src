!--------------------------------------- LICENCE BEGIN -----------------------------------
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

program midas_calcstats
  !
  ! :Purpose: Main program for computing background-error covariance file
  !           based on homogeneous and isotropic correlations.
  !
  use version_mod
  use mpi_mod
  use mpivar_mod
  use fileNames_mod
  use HorizontalCoord_mod
  use VerticalCoord_mod
  use calcstatsglb_mod
  use calcstatslam_mod
  use utilities_mod
  use ramDisk_mod
  use gridStateVector_mod
  use timeCoord_mod
  implicit none

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()

  character(len=256), parameter :: enspathname = 'ensemble'

  integer           :: fstopc
  integer           :: nulnam, ierr, fnom, fclos

  integer           :: nens   ! Ensemble size
  integer           :: ip2    ! Ensemble lead time (hour) selected within the file
  character(len=256) :: ensFileName
  character(len=60) :: mode
  character(len=4), pointer :: anlVar(:)

  NAMELIST /NAMCONF/mode
  NAMELIST /NAMENS/nens,ip2

  call ver_printNameAndVersion('calcStats','Compute the homogeneous-isotropic stats')

  !
  !- 1.  Initilization
  !
  ierr = fstopc('MSGLVL','ERRORS',0)

  !- 1.1 MPI and TMG
  call mpi_initialize
  call tmg_init(mpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  ! Setup time
  call tim_setup

  !- 1.2 Read NAMENS namelist
  nens              = 96                ! default value
  ip2               = -1                ! default value

  nulnam = 0
  ierr   = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read (nulnam,nml=namens)
  write(*     ,nml=namens)
  ierr=fclos(nulnam)

  !- 2.3 Initialize variables of the model states
  call gsv_setup

  call fln_ensfileName(ensFileName, ensPathName, memberIndex_opt=1)

  !- 1.3 Initialize the horizontal grid
  nullify(anlVar)
  call gsv_varNamesList(anlVar)
  call hco_SetupFromFile(hco_ens, ensFileName, ' ', 'Ensemble', varName_opt=anlVar(1)) ! IN

  !- 1.4 Initialize the vertical grid
  call vco_SetupFromFile( vco_ens,        & ! OUT
                          ensFileName, ' ') ! IN

  !- 1.5 Read NAMCONF namelist to find the mode
  mode  = 'BHI'  ! default value

  nulnam = 0
  ierr   = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read (nulnam,nml=namconf)
  write(*     ,nml=namconf)
  ierr   = fclos(nulnam)

  !
  !- 2. Select and launch the appropriate mode
  !

  !- 2.1 Module initialization
  if (hco_ens % global) then
     call csg_setup( nens, hco_ens, vco_ens) ! IN
  else
     call csl_setup( nens, hco_ens, vco_ens, ip2) ! IN
  end if

  !- 2.2 Mode selection
  select case(trim(mode))
  case ('BHI')
     if (hco_ens % global) then
        call csg_computeStats
     else
        call csl_computeBhi
     end if
  case ('BHI2')
     if (hco_ens % global) then
        call csg_computeStatsLatBands
     else
        call utl_abort('BHI2 mode is not available for LAM')
     end if
  case ('TOOLBOX')
     if (hco_ens % global) then
        call csg_toolbox
     else
        call csl_toolbox
     end if
  case ('STDDEV')
     if (hco_ens % global) then
        call csg_stddev
     else
        write(*,*)
        write(*,*) 'STDDEV mode is not availbale in LAM mode'
        call utl_abort('midas-calcstats')
     end if
  case ('POWERSPEC')
     if (hco_ens % global) then
        call csg_powerspec
     else
        write(*,*)
        write(*,*) 'POWERSPEC mode is not availbale in LAM mode'
        call utl_abort('midas-calcstats')
     end if
  case default
     write(*,*)
     write(*,*) 'Unknown value of MODE in global mode: ',mode
     call utl_abort('midas-calcstats')
  end select

  !
  !- 3.  Ending...
  !

  write(*,*)
  write(*,*) '---------------------'
  write(*,*) '> ENDING CALCBMATRIX '
  write(*,*) '---------------------'

  !
  !- 4.  MPI, tmg finalize
  !  
  call utl_tmg_stop(0)

  call tmg_terminate(mpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr) 

end program midas_calcstats

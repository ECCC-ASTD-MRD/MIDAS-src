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

!--------------------------------------------------------------------------
!!
!! *Purpose*: Main program for computing background-error covariance file
!!            based on homogeneous and isotropic correlations.
!!
!--------------------------------------------------------------------------
program main_calcstats
  use mpivar_mod
  use HorizontalCoord_mod
  use VerticalCoord_mod
  use calcstatsglb_mod
  use calcstatslam_mod
  use utilities_mod
  implicit none

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()

  character(len=256), parameter :: enspathname = './ensemble'
  character(len=4)   :: censnumber
  character(len=256), allocatable :: cflensin(:)

  integer           :: ens, fstopc
  integer           :: nulnam, ierr, fnom, fclos

  integer           :: nens   ! Ensemble size
  integer           :: ip1(1)
  integer           :: ip2    ! Ensemble lead time (hour) selected within the file
  character(len=256):: ensfilebasename
  character(len=60) :: mode
  character(len=2)  :: spatialDimensions

  NAMELIST /NAMCONF/mode
  NAMELIST /NAMENS/nens,ensfilebasename,spatialDimensions,ip1,ip2

  !
  !- 1.  Initilization
  !
  write(*,*)
  write(*,*) '-----------------------'
  write(*,*) '> STARTING CALCSTATS '
  write(*,*) '-----------------------'

  ierr = fstopc('MSGLVL','ERRORS',0)

  !- 1.1 MPI
  call mpi_initialize

  !- 1.2 Read NAMENS namelist
  nens              = 96                ! default value
  ensfilebasename   = '2011011918_006_' ! default value
  spatialDimensions = '3D'              ! default value
  ip1(:)            = -1                ! default value
  ip2               = -1                ! default value

  nulnam = 0
  ierr   = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read (nulnam,nml=namens)
  write(*     ,nml=namens)
  ierr=fclos(nulnam)

  allocate(cflensin(nens))
  do ens = 1,nens
    write(censnumber,'(i4.4)') ens
    cflensin(ens)= trim(enspathname) // '/' // trim(ensfilebasename) // censnumber
    write(*,*) 'ensemble file: ',ens, trim(cflensin(ens))
  end do

  !- 1.3 Initialize the horizontal grid
  call hco_SetupFromFile(hco_ens, trim(cflensin(1)), ' ', 'Ensemble' ) ! IN

  !- 1.4 Initialize the vertical grid
  select case(trim(spatialDimensions))
  case ('3D')
     call vco_SetupFromFile( vco_ens,        & ! OUT
                             cflensin(1), ' ') ! IN
  case ('2D')
     call vco_SetupManual( vco_ens,  & ! OUT
                           ip1, 1 )    ! IN
  case default
     write(*,*)
     write(*,*) 'Unknown value for spatialDimensions', spatialDimensions
     write(*,*) 'use 3D or 2D'
     call utl_abort('main_calcstats')
  end select

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
     call csg_setup( nens, cflensin, hco_ens, vco_ens) ! IN
  else
     call csl_setup( nens, cflensin, hco_ens, vco_ens, ip2) ! IN
  end if

  !- 2.2 Mode selection
  select case(trim(mode))
  case ('BHI')
     if (hco_ens % global) then
        call csg_computeStats
     else
        call csl_computeStats
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
        call utl_abort('main_calcstats')
     end if
  case ('POWERSPEC')
     if (hco_ens % global) then
        call csg_powerspec
     else
        write(*,*)
        write(*,*) 'POWERSPEC mode is not availbale in LAM mode'
        call utl_abort('main_calcstats')
     end if
  case default
     write(*,*)
     write(*,*) 'Unknown value of MODE in global mode: ',mode
     call utl_abort('main_calcstats')
  end select

  !
  !- 3.  Ending...
  !
  deallocate(cflensin)

  write(*,*)
  write(*,*) '---------------------'
  write(*,*) '> ENDING CALCBMATRIX '
  write(*,*) '---------------------'

  call rpn_comm_finalize(ierr) 

end program main_calcstats

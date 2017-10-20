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
!! *Purpose*: Main program for manipulating ensembles of states. Possible operations
!!            are (generally, all start with reading supplied ensemble):
!!            1. output_ensemble_mean
!!               - compute mean of input ensemble and write it out
!!            2. output_ensemble_perturbations
!!               - compute mean of input ensemble, subtract from
!!                 ensemble, and write out resulting ensemble perturbations
!!            3. output_recentered_ensemble
!!               - compute mean of input ensemble, subtract from 
!!                 ensemble, add supplied ensemble mean, and write out resulting ensemble
!!            4. output_ssensrf_iau_deterministic
!!               - compute mean of input background ensemble, subtract from 
!!                 ensemble, multiply by (1-tilde(K)), add random perturbations,
!!                 add supplied ensemble mean increment, and output resulting
!!                 analysis ensemble perturbations
!!               - result used as IAU increments added to deterministic background state
!!               - requires specification of sigma_obs for computing tilde(K)
!!            5. output_ssensrf_iau_ensemble
!!               - compute mean of input background ensemble, subtract from 
!!                 ensemble, multiply by (-tilde(K)), add random perturbations,
!!                 add ensemble mean increment with recentering (supplied mean 
!!                 analysis minus background ensemble mean)
!!               - result used as IAU increments added to ensemble of background states
!!               - requires specification of sigma_obs for computing tilde(K)
!!
!!            Should also add capability of adjusting analyses or analysis perturbations 
!!            to remove super-saturation, to avoid need of running AddAnalInc on each member
!!
!--------------------------------------------------------------------------
program main_ensManip
  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use gridStateVector_mod
  use ensembleStateVector_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use utilities_mod
  use ramDisk_mod
  implicit none

  type(struct_gsv) :: statevector_mean, statevector_member

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()
  type(struct_hco), pointer :: hco_ens_core => null()
  type(struct_ens)          :: ensemble

  integer              :: fclos, fnom, fstopc, newdate, ierr
  integer              :: memberIndex, stepIndex, numStep
  integer              :: idate, itime, nulnam
  integer              :: dateStamp, dateStamp_last 
  integer, allocatable :: dateStampList(:)
  integer              :: get_max_rss

  character(len=2)    :: hourstr, hourstr_last
  character(len=8)    :: datestr, datestr_last
  character(len=256)  :: ensFileName

  logical             :: makeBiPeriodic

  real(8)             :: delhh
  real(4), pointer    :: ensOneLevel(:,:,:,:)

  ! namelist variables
  character(len=2)   :: ctrlVarHumidity
  character(len=256) :: ensPathName, ensFileBaseName
  logical  :: write_mpi, output_ensemble_mean, output_ensemble_perturbations
  integer  :: nEns, date
  NAMELIST /NAMENSMANIP/nEns, date, ensPathName, ensFileBaseName, ctrlVarHumidity, write_mpi,  &
                        output_ensemble_mean, output_ensemble_perturbations

  write(*,'(/,' //  &
        '3(" *****************"),/,' //                   &
        '14x,"-- START OF MAIN_ENSMANIP             --",/,' //   &
        '14x,"-- Program for general manipulation of ensembles --",/, ' //  &
        '14x,"-- Revision number ",a," --",/,' //  &
        '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  !
  !- 0. MPI, TMG initialization
  !
  call mpi_initialize
  call tmg_init(mpi_myid, 'TMG_ENSMANIP' )

  call tmg_start(1,'MAIN')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  !
  !- 1. Set/Read values for the namelist NAMENSMANIP
  !
  
  !- 1.1 Setting default values
  nEns                          = 10
  date                          = 1900120100
  write_mpi                     = .false.
  ensPathName                   = 'ensemble'
  ensFileBaseName               = ''
  ctrlVarHumidity               = 'HU'
  output_ensemble_mean          = .false.
  output_ensemble_perturbations = .false.

  !- 1.2 Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namensmanip, iostat=ierr)
  if ( ierr /= 0) call utl_abort('main_ensManip: Error reading namelist')
  if ( mpi_myid == 0 ) write(*,nml=namensmanip)
  ierr = fclos(nulnam)

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 2.  Initialization
  !

  !- 2.1 Initialize the dates

  ! Analysis time
  idate   = date/100
  itime   = (date-idate*100)*1000000
  ierr    = newdate(dateStamp, idate, itime, 3)
  write(datestr,'(i8.8)') idate
  write(hourstr,'(i2.2)') itime/1000000
  if ( mpi_myid == 0 ) write(*,*)' datestr= ', datestr, ' hourstr= ', hourstr
  if ( mpi_myid == 0 ) write(*,*)' dateStamp= ', dateStamp

  ! Setup timeCoord module
  call tim_setup
  call tim_setDatestamp(dateStamp)
  numStep = tim_nstepobsinc
  allocate(dateStampList(numStep))
  call tim_getstamplist(dateStampList,numStep,tim_getDatestamp())

  ! Previous analysis time
  delhh = -tim_windowsize
  call incdatr( dateStamp_last, tim_getDatestamp(), delhh )
  ierr = newdate( dateStamp_last, idate, itime, -3 )
  write(datestr_last,'(i8.8)') idate
  write(hourstr_last,'(i2.2)') itime/1000000
  if ( mpi_myid == 0 ) write(*,*)' datestr_last= ', datestr_last, ' hourstr_last= ', hourstr_last
  if ( mpi_myid == 0 ) write(*,*)' dateStamp_last= ', dateStamp_last

  !- 2.3 Initialize variables of the model states
  call gsv_setup

  !- 2.4 Initialize the Ensemble grid
  if (mpi_myid == 0) write(*,*)''
  if (mpi_myid == 0) write(*,*)' Set hco parameters for ensemble grid'
  ! Use the first ensemble member to initialize the grid
  call ens_fileName( ensFileName, ensPathName, ensFileBaseName, 1 )
  call hco_SetupFromFile( hco_ens, ensFileName, ' ', 'ENSFILEGRID')
  call vco_setupFromFile( vco_ens, ensFileName )

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !- 2.5 Setup and read the background ensemble
  call tmg_start(2,'READ_ENSEMBLE')
  call ens_allocate(ensemble, nEns, numStep, hco_ens, vco_ens, dateStampList)
  makeBiPeriodic = .false.
  call ens_readEnsemble( ensemble, ensPathName, ensFileBaseName, makeBiPeriodic, ctrlVarHumidity )
  call tmg_stop(2)

  !- 2.6 Compute ensemble mean and subtract it from ensemble
  call tmg_start(3,'COMPUTE_MEAN')
  call ens_computeMean( ensemble )
  call ens_removeMean( ensemble )
  call tmg_stop(3)

  !- 3.0 Output the background ensemble mean, if requested
  if ( output_ensemble_mean ) then
    call tmg_start(4,'OUTPUT_MEAN')
    call ens_copyEnsMean( ensemble, statevector_mean )

    ! Filename for background ensemble mean
    ensFileName = './' // trim(ensfilebasename) // &
                  trim(datestr_last) // trim(hourstr_last) // '_006_ensmean'

    do stepIndex = 1, numStep
      if ( mpi_myid == 0 ) write(*,*) 'main_ensManip: writing time step ', stepIndex
      if ( write_mpi ) then
        call gsv_writeToFileMPI( statevector_mean, ensFileName, 'ENSMEAN', indexStep_in = stepIndex, typvar_in = 'P' )
      else
        call gsv_writeToFile( statevector_mean, ensFileName, 'ENSMEAN', indexStep_in = stepIndex, typvar_in = 'P' )
      end if
    end do

    call tmg_stop(4)
  end if

  !- 4.0 Output the background ensemble perturbations, if requested
  if ( output_ensemble_perturbations ) then
    call tmg_start(8,'OUTPUT_PERTURBATIONS')
    call ens_writeEnsemble( ensemble, './', 'pert_', ctrlVarHumidity, 'ENSPERT', 'P')
    call tmg_stop(8)
  end if

  !
  !- 6.  MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_ENSMANIP' )
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 7.  Ending
  !
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'
  if ( mpi_myid == 0 ) write(*,*) ' MAIN_ENSMANIP ENDS'
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'

end program main_ensManip
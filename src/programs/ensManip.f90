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
!!            2. output_ensemble_stddev
!!               - compute the input ensemble stddev and write it out
!!            3. output_recentered_ensemble
!!               - compute mean of input ensemble, subtract from 
!!                 ensemble, add supplied ensemble mean, and write out resulting ensemble
!!            4. output_ensemble_perturbations
!!               - compute mean of input ensemble, subtract from
!!                 ensemble, and write out resulting ensemble perturbations
!!
!!
!--------------------------------------------------------------------------
program midas_ensManip
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

  type(struct_gsv) :: statevector_mean, statevector_stddev, statevector_recenteringMean

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
  character(len=256)  :: ensFileName, recenteringMeanFileName

  logical             :: makeBiPeriodic

  real(8)             :: delhh
  real(4), pointer    :: ensOneLevel(:,:,:,:)

  ! namelist variables
  character(len=2)   :: ctrlVarHumidity
  character(len=256) :: ensPathName, ensFileBaseName
  logical  :: write_mpi, output_ensemble_mean, output_ensemble_stddev, output_ensemble_perturbations, recenter
  real(8)  :: recentering_coeff
  integer  :: nEns, date, numBits
  NAMELIST /NAMENSMANIP/nEns, date, ensPathName, ensFileBaseName, ctrlVarHumidity, write_mpi,  &
                        output_ensemble_mean, output_ensemble_stddev, output_ensemble_perturbations, &
                        recenter, recentering_coeff, numBits

  write(*,'(/,' //  &
        '3(" *****************"),/,' //                   &
        '14x,"-- START OF MIDAS-ENSMANIP             --",/,' //   &
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
  output_ensemble_stddev        = .false.
  output_ensemble_perturbations = .false.
  numBits                       = 16

  !- 1.2 Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namensmanip, iostat=ierr)
  if ( ierr /= 0) call utl_abort('midas-ensManip: Error reading namelist')
  if ( mpi_myid == 0 ) write(*,nml=namensmanip)
  ierr = fclos(nulnam)

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 2.  Initialization
  !

  !- 2.1 Initialize the dates

  ! Setup timeCoord module
  call ens_fileName( ensFileName, ensPathName, ensFileBaseName, 1 )
  call tim_setup( fileNameForDate_opt=ensFileName )
  numStep = tim_nstepobsinc
  allocate(dateStampList(numStep))
  call tim_getstamplist(dateStampList,numStep,tim_getDatestamp())

  ! Previous analysis time
  delhh = -tim_windowsize
  call incdatr( dateStamp_last, tim_getDatestamp(), delhh )
  ierr = newdate( dateStamp_last, idate, itime, -3 )
  write(datestr_last,'(i8.8)') idate
  write(hourstr_last,'(i2.2)') itime/1000000
  if ( mpi_myid == 0 ) write(*,*) 'midas-ensManip: datestr_last= ', datestr_last, ' hourstr_last= ', hourstr_last
  if ( mpi_myid == 0 ) write(*,*) 'midas-ensManip: dateStamp_last= ', dateStamp_last

  !- 2.3 Initialize variables of the model states
  call gsv_setup

  !- 2.4 Initialize the Ensemble grid
  if (mpi_myid == 0) write(*,*) ''
  if (mpi_myid == 0) write(*,*) 'midas-ensManip: Set hco parameters for ensemble grid'
  ! Use the first ensemble member to initialize the grid
  call ens_fileName( ensFileName, ensPathName, ensFileBaseName, 1 )
  call hco_SetupFromFile( hco_ens, ensFileName, ' ', 'ENSFILEGRID')
  call vco_setupFromFile( vco_ens, ensFileName )

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !- 2.5 Setup and read the ensemble
  call tmg_start(2,'READ_ENSEMBLE')
  call ens_allocate(ensemble, nEns, numStep, hco_ens, vco_ens, dateStampList)
  makeBiPeriodic = .false.
  call ens_readEnsemble( ensemble, ensPathName, ensFileBaseName, makeBiPeriodic, ctrlVarHumidity )
  call tmg_stop(2)

  !- 2.6 Compute ensemble mean
  call tmg_start(3,'COMPUTE_MEAN')
  call ens_computeMean( ensemble )
  call tmg_stop(3)

  !- 3.0 Output the ensemble mean, if requested
  if ( output_ensemble_mean ) then
    call tmg_start(4,'OUTPUT_MEAN')
    call ens_copyEnsMean( ensemble, statevector_mean )

    ! Filename for ensemble mean
    ensFileName = './' // trim(ensfilebasename) // &
                  trim(datestr_last) // trim(hourstr_last) // '_006_ensmean'

    do stepIndex = 1, numStep
      if ( mpi_myid == 0 ) write(*,*) 'midas-ensManip: writing time step ', stepIndex
      if ( write_mpi ) then
        call gsv_writeToFileMPI( statevector_mean, ensFileName, 'ENSMEAN', indexStep_in = stepIndex, typvar_in = 'P')
      else
        call gsv_writeToFile( statevector_mean, ensFileName, 'ENSMEAN', indexStep_in = stepIndex, typvar_in = 'P', numBits_opt = numBits)
      end if
    end do

    call tmg_stop(4)
  end if

  !- 4.0 Compute and output the ensemble spread stddev, if requested
  if ( output_ensemble_stddev ) then
    ! Compute the ensemble stddev and put in statevector_stddev
    call tmg_start(6,'COMPUTE_STDDEV')
    call ens_computeStdDev( ensemble )
    call tmg_stop(6)

    call tmg_start(7,'OUTPUT_STDDEV')
    call ens_copyEnsStdDev( ensemble, statevector_stddev )

    ! Filename for ensemble stddev
    ensFileName = './' // trim(ensfilebasename) // &
                  trim(datestr_last) // trim(hourstr_last) // '_006_ensstddev'

    ! Output the ensemble stddev
    do stepIndex = 1, numStep
      if ( mpi_myid == 0 ) write(*,*) 'midas-ensManip: writing time step ', stepIndex
      if ( write_mpi ) then
        call gsv_writeToFileMPI( statevector_stddev, ensFileName, 'ENSMEAN', indexStep_in = stepIndex, typvar_in = 'P')
      else
        call gsv_writeToFile( statevector_stddev, ensFileName, 'ENSMEAN', indexStep_in = stepIndex, typvar_in = 'P' , numBits_opt = numBits)
      end if
    end do

    call tmg_stop(7)
  end if

  if ( output_ensemble_perturbations .and. recenter ) then
     call utl_abort('midas-ensManip: You must choose between computing ensemble perturbations and recenter.')
  end if

  !- 5.0 Output the ensemble perturbations, if requested
  if ( output_ensemble_perturbations ) then
    call tmg_start(8,'OUTPUT_PERTURBATIONS')
    call ens_removeMean( ensemble )
    call ens_writeEnsemble( ensemble, './', 'pert_', ctrlVarHumidity, 'ENSPERT', 'P', numBits_opt = numBits)
    call tmg_stop(8)
  end if

  if (recenter) then
     ! read recentering mean in file '${DATE}_recenteringmean'
     call tmg_start(10,'READ_RECENTERINGMEAN')

     ! Filename for recentering mean
     recenteringMeanFileName = './' // trim(datestr_last) // trim(hourstr_last) // '_recenteringmean'

     call gsv_allocate(statevector_recenteringMean, numStep, hco_ens, vco_ens, &
          dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.)

     do stepIndex = 1, numStep
        dateStamp = datestamplist(stepIndex)
        if(mpi_myid == 0) write(*,*) ''
        if(mpi_myid == 0) write(*,*) 'midas-ensManip: reading recentering mean for time step: ',stepIndex, dateStamp
        call gsv_readFromFile(statevector_recenteringMean, trim(recenteringMeanFileName), ' ', ' ', stepIndex)
     end do

     call tmg_stop(10)

     ! Compute 'x_recenteringMean - x_ensembleMean' which is stored in 'statevector_recenteringMean'
     call gsv_add(statevector_mean, statevector_recenteringMean, -1.0d0)

     call tmg_start(11,'RECENTER_ENSEMBLE_MEMBERS')
     call ens_recenter(ensemble,statevector_recenteringMean,recentering_coeff)
     call tmg_stop(11)

     call tmg_start(12,'OUTPUT_RECENTER_MEMBERS')
     call ens_writeEnsemble( ensemble, './', 'recentered_', ctrlVarHumidity, 'ENSRECENTER', 'P', numBits_opt = numBits)
     call tmg_stop(12)
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
  if ( mpi_myid == 0 ) write(*,*) ' MIDAS-ENSMANIP ENDS'
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'

end program midas_ensManip

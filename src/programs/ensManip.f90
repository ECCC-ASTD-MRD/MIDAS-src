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
  !
  ! **Purpose**: Main program for manipulating ensembles of states. Possible operations
  ! are (generally, all start with reading supplied ensemble):
  !
  ! 1. output_ensemble_mean
  !    - compute mean of input ensemble and write it out
  ! 2. output_ensemble_stddev
  !    - compute the input ensemble stddev and write it out
  ! 3. output_recentered_ensemble
  !    - compute mean of input ensemble, subtract from 
  !    ensemble, add supplied ensemble mean, and write out resulting ensemble
  ! 4. output_ensemble_perturbations
  !    - compute mean of input ensemble, subtract from
  !    ensemble, and write out resulting ensemble perturbations
  !
  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use gridStateVector_mod
  use fileNames_mod
  use ensembleStateVector_mod
  use verticalCoord_mod
  use analysisGrid_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use utilities_mod
  use ramDisk_mod
  implicit none

  type(struct_gsv) :: statevector_mean, statevector_stddev
  type(struct_gsv) :: statevector_recenteringMean, statevector_alternativeEnsembleMean

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()
  type(struct_hco), pointer :: hco_ens_core => null()
  type(struct_ens)          :: ensemble

  integer              :: fclos, fnom, fstopc, ierr
  integer              :: memberIndex, stepIndex, numStep
  integer              :: nulnam, dateStamp
  integer, allocatable :: dateStampList(:)
  integer              :: get_max_rss

  character(len=256)  :: ensFileName, ensFileBaseName, recenteringMeanFileName, alternativeEnsembleMeanFileName

  logical             :: makeBiPeriodic

  real(4), pointer    :: ensOneLevel(:,:,:,:)

  ! namelist variables
  character(len=1)   :: ensembleTypVarOutput
  character(len=12)  :: hInterpolationDegree
  character(len=14)  :: ensembleEtiketOutput, ensembleControlMemberEtiket
  character(len=2)   :: ctrlVarHumidity
  character(len=256) :: ensPathName, alternativeEnsembleMean
  logical  :: output_ensemble_mean, output_ensemble_stddev, output_ensemble_perturbations
  logical  :: recenter, ensembleEtiketOutputAppendMemberNumber, recenterEnsembleControlMember
  real(8)  :: recentering_coeff
  integer  :: nEns, numBits
  NAMELIST /NAMENSMANIP/nEns, ensPathName, ctrlVarHumidity, alternativeEnsembleMean, ensembleEtiketOutput, ensembleTypVarOutput, &
                        output_ensemble_mean, output_ensemble_stddev, output_ensemble_perturbations, &
                        recenter, recentering_coeff, numBits, ensembleEtiketOutputAppendMemberNumber, &
                        recenterEnsembleControlMember, ensembleControlMemberEtiket, hInterpolationDegree

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
  ensPathName                   = 'ensemble'
  alternativeEnsembleMean       = ''
  ensembleTypVarOutput          = 'P'
  ensembleEtiketOutput          = 'ENSRECENTER'
  ensembleEtiketOutputAppendMemberNumber = .false.
  ensembleControlMemberEtiket   = 'ENSCTLMEM'
  recenterEnsembleControlMember    = .false.
  ctrlVarHumidity               = 'HU'
  output_ensemble_mean          = .false.
  output_ensemble_stddev        = .false.
  output_ensemble_perturbations = .false.
  numBits                       = 16
  recenter                      = .false.
  recentering_coeff             = 1.0
  hInterpolationDegree             = 'CUBIC' ! or "LINEAR" or "NEAREST"

  !- 1.2 Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namensmanip, iostat=ierr)
  if ( mpi_myid == 0 ) write(*,nml=namensmanip)
  if ( ierr /= 0) call utl_abort('midas-ensManip: Error reading namelist')
  ierr = fclos(nulnam)

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 2.  Initialization
  !

  !- 2.1 Initialize the dates

  ! Setup timeCoord module
  call fln_ensFileName( ensFileName, ensPathName, 1, ensFileBaseName_opt=ensFileBaseName)
  call tim_setup( fileNameForDate_opt=ensFileName )
  numStep = tim_nstepobsinc
  allocate(dateStampList(numStep))
  call tim_getstamplist(dateStampList,numStep,tim_getDatestamp())

  !- 2.3 Initialize variables of the model states
  call gsv_setup

  !- 2.4 Initialize the Ensemble grid
  if (mpi_myid == 0) write(*,*) ''
  if (mpi_myid == 0) write(*,*) 'midas-ensManip: Set hco parameters for ensemble grid'

  ! Use the first ensemble member to initialize the grid
  call hco_SetupFromFile( hco_ens, ensFileName, ' ', 'ENSFILEGRID')
  call vco_setupFromFile( vco_ens, ensFileName )

  if ( hco_ens % global ) then
    call agd_SetupFromHCO( hco_ens ) ! IN
  else
    !- Setup the LAM analysis grid metrics
    hco_ens_core => hco_ens
    call agd_SetupFromHCO( hco_ens, hco_ens_core ) ! IN
  end if

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !- 2.5 Setup and read the ensemble
  call tmg_start(2,'READ_ENSEMBLE')
  call ens_allocate(ensemble, nEns, numStep, hco_ens, vco_ens, dateStampList)
  makeBiPeriodic = .false.
  call ens_readEnsemble( ensemble, ensPathName, makeBiPeriodic, ctrlVarHumidity )
  call tmg_stop(2)

  if (trim(alternativeEnsembleMean) == '' .or. output_ensemble_mean) then
    !- 2.6 Compute ensemble mean
    call tmg_start(3,'COMPUTE_MEAN')
    call ens_computeMean( ensemble )
    call tmg_stop(3)
  end if

  !- 3.0 Output the ensemble mean, if requested
  if ( output_ensemble_mean ) then
    call tmg_start(4,'OUTPUT_MEAN')
    call ens_copyEnsMean( ensemble, statevector_mean )

    ! Filename for ensemble mean
    ensFileName = './' // trim(ensFileBaseName) // '_ensmean'

    do stepIndex = 1, numStep
      if ( mpi_myid == 0 ) write(*,*) 'midas-ensManip: writing time step ', stepIndex
      call gsv_writeToFile( statevector_mean, ensFileName, 'ENSMEAN',  &
           stepIndex_opt = stepIndex, typvar_opt = 'P', numBits_opt = numBits)
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
    ensFileName = './' // trim(ensFileBaseName) // '_ensstddev'

    ! Output the ensemble stddev
    do stepIndex = 1, numStep
      if ( mpi_myid == 0 ) write(*,*) 'midas-ensManip: writing time step ', stepIndex
      call gsv_writeToFile( statevector_stddev, ensFileName, 'ENSSTDDEV',   &
           stepIndex_opt = stepIndex, typvar_opt = 'P' , numBits_opt = numBits)
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
    call ens_writeEnsemble( ensemble, '.', 'pert_', ctrlVarHumidity, 'ENSPERT', 'P', numBits_opt = numBits)
    call tmg_stop(8)
  end if

  if (recenter) then
    ! read recentering mean in file '${DATE}_recenteringmean'
    call tmg_start(10,'READ_RECENTERINGMEAN')

    ! Filename for recentering mean
    recenteringMeanFileName = './' // trim(ensFileBaseName) // '_recenteringmean'

    call gsv_allocate(statevector_recenteringMean, numStep, hco_ens, vco_ens, &
         dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
         hInterpolateDegree_opt = hInterpolationDegree)

    do stepIndex = 1, numStep
      dateStamp = datestamplist(stepIndex)
      if(mpi_myid == 0) write(*,*) ''
      if(mpi_myid == 0) write(*,*) 'midas-ensManip: reading recentering mean for time step: ',stepIndex, dateStamp
      call gsv_readFromFile(statevector_recenteringMean, trim(recenteringMeanFileName), ' ', ' ',  &
           stepIndex_opt=stepIndex, unitConversion_opt=.true.)
    end do

    call tmg_stop(10)

    if (trim(alternativeEnsembleMean) /= '') then
      ! read ensemble center in file '${ensFileBaseName}_${alternativeEnsembleMean}'
      call tmg_start(11,'READ_ALTERNATIVEENSEMBLEMEAN')

      ! Filename for ensemble center
      alternativeEnsembleMeanFileName = './' // trim(ensFileBaseName) // trim(alternativeEnsembleMean)

      call gsv_allocate(statevector_alternativeEnsembleMean, numStep, hco_ens, vco_ens, &
           dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
           hInterpolateDegree_opt = hInterpolationDegree)

      do stepIndex = 1, numStep
        dateStamp = datestamplist(stepIndex)
        if(mpi_myid == 0) write(*,*) ''
        if(mpi_myid == 0) write(*,*) 'midas-ensManip: reading ensemble center for time step: ',stepIndex, dateStamp, trim(alternativeEnsembleMeanFileName)
        call gsv_readFromFile(statevector_alternativeEnsembleMean, trim(alternativeEnsembleMeanFileName), ' ', ' ',  &
             stepIndex_opt=stepIndex, unitConversion_opt=.true.)
      end do

      call tmg_stop(11)

      call tmg_start(12,'RECENTER_ENSEMBLE_MEMBERS')
      call ens_recenter(ensemble,statevector_recenteringMean,recentering_coeff,alternativeEnsembleMean_opt=statevector_alternativeEnsembleMean)
      call tmg_stop(12)

      if (recenterEnsembleControlMember) then
        call ens_recenterControlMember(ensemble,hco_ens,vco_ens,ensFileBaseName,'.', 'recentered_', &
             statevector_recenteringMean, recentering_coeff, ensembleControlMemberEtiket, &
             ensembleTypVarOutput, alternativeEnsembleMean_opt=statevector_alternativeEnsembleMean, numBits_opt = numBits)
      end if
    else
      call tmg_start(12,'RECENTER_ENSEMBLE_MEMBERS')
      call ens_recenter(ensemble,statevector_recenteringMean,recentering_coeff)
      call tmg_stop(12)

      if (recenterEnsembleControlMember) then
        call ens_recenterControlMember(ensemble,hco_ens,vco_ens,ensFileBaseName,ensPathName, 'recentered_', &
             statevector_recenteringMean, recentering_coeff, ensembleControlMemberEtiket, &
             ensembleTypVarOutput, numBits_opt = numBits)
      end if
    end if ! end of 'else' related to 'if (trim(alternativeEnsembleMean) /= '')'

    call tmg_start(130,'OUTPUT_RECENTER_MEMBERS')
    call ens_writeEnsemble( ensemble, '.', 'recentered_', ctrlVarHumidity, ensembleEtiketOutput, ensembleTypVarOutput, &
         numBits_opt = numBits, etiketAppendMemberNumber_opt = ensembleEtiketOutputAppendMemberNumber)
    call tmg_stop(130)
  end if ! end of 'if (recenter)'

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

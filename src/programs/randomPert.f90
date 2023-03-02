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

program midas_randomPert
  !
  ! :Purpose: Main program for generating an ensemble of random perturbations
  !           based on the B matrix (can be homogeneous/isotropic or
  !           ensemble-based).
  !
  use version_mod
  use midasMpi_mod
  use ramDisk_mod
  use controlVector_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use bmatrix_mod
  use interpolation_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use randomNumber_mod
  use utilities_mod
  use gridVariableTransforms_mod
  implicit none

  type(struct_gsv) :: stateVectorPert, stateVectorPertInterp
  type(struct_gsv) :: stateVectorEnsMean, stateVectorIce

  type(struct_vco), pointer :: vco_anl => null()
  type(struct_hco), pointer :: hco_anl => null()
  type(struct_hco), pointer :: hco_anlcore => null()
  type(struct_hco), pointer :: hco_target => null()
  type(struct_hco), pointer :: hco_targetcore => null()

  real(4), pointer :: seaice_ptr(:,:,:)
  real(8), pointer :: field(:,:,:), fieldInterp(:,:,:)

  integer :: fclos, fnom, fstopc, newdate, dateStamp, datePrevious, dateStampPrevious
  integer :: imode, ierr
  integer :: memberIndex, lonIndex, latIndex, cvIndex, levIndex, nkgdim
  integer :: datePrint, timePrint, nulnam, randomSeed
  integer :: get_max_rss, n_grid_point, n_grid_point_glb

  integer :: latPerPEa, latPerPEmaxa, myLatBega, myLatEnda
  integer :: lonPerPEa, lonPerPEmaxa, myLonBega, myLonEnda
  integer :: latPerPEt, latPerPEmaxt, myLatBegt, myLatEndt
  integer :: lonPerPEt, lonPerPEmaxt, myLonBegt, myLonEndt

  real(8) :: previousDateFactPrev, previousDateFactNoise
  real(8), allocatable :: controlVector(:), controlVector_mpiglobal(:)
  real(8), allocatable :: gdmean(:,:,:)
  real(4), allocatable :: ensemble_r4(:,:,:,:)
  real(4), allocatable :: ensemblePreviousDate_r4(:,:,:,:)
  real(8), allocatable :: avg_pturb_var(:)
  real(8), allocatable :: avg_pturb_var_glb(:)
  real(8), allocatable :: pturb_var(:,:,:)

  logical :: targetGridExists
  character(len=10) :: dateString, datePreviousString
  character(len=4)  :: memberString
  character(len=25) :: outFileName, inFileName
  character(len=64) :: ensMeanFileName = 'ensMeanState'
  character(len=2)  :: typvarOut

  ! Namelist variables
  logical :: remove_mean          ! choose to remove mean from perturbations
  logical :: smoothVariances      ! choose to impose horizontally constant perturbation variances
  logical :: mpiTopoIndependent   ! choose to compute random numbers with mpi-topology-independent method 
  logical :: readEnsMean          ! choose to read ens mean and add this to the perturbations
  logical :: setPertZeroUnderIce  ! choose to set perturbation to zero under sea ice (for SST)
  integer :: nens                 ! number of perturbations to compute
  integer :: seed                 ! initial value of the random seed
  integer :: numBits              ! number of bits to use when writing to standard files
  character(len=12) :: out_etiket ! the 'etiket' to write to standard files
  real(4) :: iceFractionThreshold ! ice fraction threshold to use in combination with 'setPertZeroUnderIce'
  real(4) :: previousDateFraction ! relative amount of previous date perturbations to include in current perturbations
  NAMELIST /NAMENKF/nens, seed, out_etiket, remove_mean,  &
                    smoothVariances, mpiTopoIndependent, numBits, &
                    readEnsMean, setPertZeroUnderIce, iceFractionThreshold, &
                    previousDateFraction

  call ver_printNameAndVersion('randomPert','Generation of random perturbations')

  !
  !- 0. MPI, tmg initialization
  !
  call mmpi_initialize
  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ierr = fstopc('MSGLVL','ERRORS',0)

  !
  !- 1. Set/Read values for the namelist NAMENKF
  !
  
  !- 1.1 Setting default values
  nens  = 10
  seed  = -999        ! If -999, set random seed using the date
  remove_mean = .true.
  out_etiket='RANDOM_PERT' 
  smoothVariances = .false.
  mpiTopoIndependent = .false.
  numBits = 32
  readEnsMean = .false.
  setPertZeroUnderIce = .false.
  iceFractionThreshold = 0.2
  previousDateFraction = -1.0

  !- 1.2 Read the namelist
  nulnam=0
  ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namenkf, iostat=ierr)
  if(ierr.ne.0) call utl_abort('midas-randomPert: Error reading namelist')
  if( mmpi_myid == 0 ) write(*,nml=namenkf)
  ierr=fclos(nulnam)

  if (readEnsMean) then
    typvarOut = 'A'
  else
    typvarOut = 'R'
  end if

  if (previousDateFraction > 1.0) then
    call utl_abort('midas-randomPert: previousDateFraction must be less than 1.0 to give stable results')
  end if

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 2.  Initialization
  !

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  !- 2.1 Set the dateStamp, either from env variable or ensMeanState file

  !- Initialize the Temporal grid and the dateStamp from env variable
  call tim_setup()

  ! If dateStamp not set by env variable, use date from ensMeanState, if available
  if (tim_getDateStamp() == 0) then
    if (readEnsMean) then
      dateStamp = tim_getDatestampFromFile(ensMeanFileName)
      call tim_setDatestamp(dateStamp)
    else
      call utl_abort('midas-randomPert: DateStamp must be set through env variable')
    end if
  end if
  dateStamp = tim_getDateStamp()
  imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
  ierr = newdate(dateStamp, datePrint, timePrint, imode)
  write(dateString, '(I10)') datePrint*100 + timePrint/1000000
  if( mmpi_myid == 0 ) then
    write(*,*) ' date= ', datePrint, ' time= ', timePrint, ' stamp= ', dateStamp
    write(*,*) ' dateString = ', dateString
  end if

  !- 2.2 Initialize variables of the model states
  call gsv_setup

  !- 2.3 Initialize the horizontal grids from analysisgrid and targetgrid files
  if (mmpi_myid == 0) write(*,*)
  if (mmpi_myid == 0) write(*,*) 'Set hco parameters for analysis grid'
  call hco_setupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    hco_anlcore => hco_anl
  else
    !- Iniatilized the core (Non-Exteded) analysis grid
    call hco_setupFromFile( hco_anlcore, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
  end if

  inquire(file='./targetgrid',exist=targetGridExists)
  if (targetGridExists) then
    if (mmpi_myid == 0) write(*,*)
    if (mmpi_myid == 0) write(*,*) 'Set hco parameters for target grid'

    call hco_setupFromFile(hco_target, './targetgrid', 'ANALYSIS', 'Target' ) ! IN

    if ( hco_target % global ) then
      hco_targetcore => hco_target
    else
      !- Iniatilized the core (Non-Exteded) analysis grid
      call hco_setupFromFile( hco_targetcore, './targetgrid', 'COREGRID', 'TargetCore' ) ! IN
    end if
  else
    ! no targetgrid file, so use analysisgrid grid instead
    hco_target => hco_anl
    hco_targetcore => hco_anlcore
  end if

  call mmpi_setup_latbands(hco_anl % nj,                & ! IN
                             latPerPEa, latPerPEmaxa, myLatBega, myLatEnda ) ! OUT
  call mmpi_setup_lonbands(hco_anl % ni,                & ! IN
                             lonPerPEa, lonPerPEmaxa, myLonBega, myLonEnda ) ! OUT

  call mmpi_setup_latbands(hco_target % nj,                & ! IN
                             latPerPEt, latPerPEmaxt, myLatBegt, myLatEndt ) ! OUT
  call mmpi_setup_lonbands(hco_target % ni,                & ! IN
                             lonPerPEt, lonPerPEmaxt, myLonBegt, myLonEndt ) ! OUT

  !- 2.4 Initialize the vertical coordinate from the analysisgrid file
  call vco_setupFromFile(vco_anl, './analysisgrid', ' ')
 
  !- 2.5 Initialize the B_hi matrix
  call bmat_setup(hco_anl, hco_anlcore, vco_anl)

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !- 2.6 Initialize the gridded variable transform module
  call gvt_setup(hco_anl,hco_anlcore,vco_anl)
  if ( gsv_varExist(varName='HU') ) call gvt_setupRefFromTrialFiles('HU')

  !- 2.7 Set randomSeed
  if (seed == -999) then
    ! If "seed" namelist value is -999, set random seed using the date
    dateStamp = tim_getDateStamp()
    if (dateStamp == -1) then
      call utl_abort('midas-randomPert: dateStamp is not set, cannot be used to set random seed')
    end if
    imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
    ierr = newdate(dateStamp, datePrint, timePrint, imode)
    timePrint = timePrint/1000000
    datePrint =  datePrint*100 + timePrint
    ! Remove the century, keeping 2 digits of the year
    randomSeed = datePrint - 100000000*(datePrint/100000000)
  else
    ! Otherwise, use value from namelist
    randomSeed = seed
  end if
  write(*,*) 'midas-randomPert: randomSeed for set to ', randomSeed

  !
  !- 3. Memory allocations
  !

  !- 3.1 Allocate the stateVectorPert
  call gsv_allocate(stateVectorPert, 1, hco_anl, vco_anl, &
                    dateStamp_opt=dateStamp, mpi_local_opt=.true., &
                    allocHeight_opt=.false., allocPressure_opt=.false., &
                    hInterpolateDegree_opt='LINEAR')
  nkgdim = stateVectorPert%nk
  allocate(ensemble_r4(myLonBega:myLonEnda, myLatBega:myLatEnda, nkgdim, nEns))

  !- 3.2 Allocate auxillary variables
  allocate(gdmean(myLonBega:myLonEnda, myLatBega:myLatEnda, nkgdim))
  allocate(controlVector(cvm_nvadim))
  allocate(controlVector_mpiglobal(cvm_nvadim_mpiglobal))

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 4. Compute an ensemble of random perturbations
  !

  if( mmpi_myid == 0 ) write(*,*) '******************'
  if( mmpi_myid == 0 ) write(*,*) 'midas-randomPert: COMPUTE the mean of the random ', &
                                  'perturbations of all the members'

  if( mpiTopoIndependent ) then
    call rng_setup(abs(randomSeed))
  else
    call rng_setup(abs(randomSeed+mmpi_myid))
  end if

  gdmean(:,:,:) = 0.0D0
  call gsv_getField(stateVectorPert,field)

  !- 4.1 Generate a (potentially) biased ensemble
  do memberIndex = 1, NENS

    if( mmpi_myid == 0 ) write(*,*) '...computing member number= ', memberIndex

    !- 4.1.1 Create a random control vector in spectral space

    if( mpiTopoIndependent ) then
      !- Global vector (for testing different mpi topology, less efficient)
      do cvIndex = 1, cvm_nvadim_mpiglobal
        controlVector_mpiglobal(cvIndex) = rng_gaussian()
      end do
      call bmat_reduceToMPILocal( controlVector, controlVector_mpiglobal )
    else
      !- Local vector (different randomSeed for each processor, more efficient)
      do cvIndex = 1, cvm_nvadim
        controlVector(cvIndex) = rng_gaussian()
      end do
    end if

    !- 4.1.2 Transform to control variables in physical space
    call bmat_sqrtB(controlVector, cvm_nvadim, & ! IN
                    stateVectorPert           )  ! OUT

    !- 4.1.3 Copy perturbations to big array and update ensemble sum
    !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex)    
    do levIndex = 1, nkgdim
      do latIndex = myLatBega, myLatEnda
        do lonIndex = myLonBega, myLonEnda
          ensemble_r4(lonIndex, latIndex, levIndex, memberIndex) = real(field(lonIndex, latIndex, levIndex), 4)
          gdmean(lonIndex, latIndex, levIndex) = gdmean(lonIndex, latIndex, levIndex) + field(lonIndex, latIndex, levIndex)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end do

  call gsv_deallocate(stateVectorPert)
  
  !- 4.2 Remove the ensemble mean
  if ( REMOVE_MEAN ) then

    !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex)    
    do levIndex = 1, nkgdim
      do latIndex = myLatBega, myLatEnda
        do lonIndex = myLonBega, myLonEnda
          gdmean(lonIndex, latIndex, levIndex) = gdmean(lonIndex, latIndex, levIndex) / real(NENS, 8)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE (memberIndex, levIndex, latIndex, lonIndex)    
    do memberIndex = 1, NENS
      do levIndex = 1, nkgdim
        do latIndex = myLatBega, myLatEnda
          do lonIndex = myLonBega, myLonEnda
            ensemble_r4(lonIndex, latIndex, levIndex, memberIndex) =  &
                ensemble_r4(lonIndex, latIndex, levIndex, memberIndex) -  &
                real(gdmean(lonIndex, latIndex, levIndex), 4)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end if
  
  !- 4.3 Smooth variances to horizontally constant values
  if ( smoothVariances ) then
  
    allocate(pturb_var(myLonBega:myLonEnda, myLatBega:myLatEnda, stateVectorPert%nk))
    allocate(avg_pturb_var(stateVectorPert%nk), avg_pturb_var_glb(stateVectorPert%nk))
    pturb_var(:,:,:) = 0.0D0
    avg_pturb_var(:) = 0.0D0
    avg_pturb_var_glb(:) = 0.0D0
  
    do memberIndex = 1, NENS  
      !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex)  
      do levIndex = 1, nkgdim
        do latIndex = myLatBega, myLatEnda
          do lonIndex = myLonBega, myLonEnda
               pturb_var(lonIndex, latIndex, levIndex) = pturb_var(lonIndex, latIndex, levIndex) +  &
                   real(ensemble_r4(lonIndex, latIndex, levIndex, memberIndex), 8)**2
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end do
  
    !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex)  
    do levIndex = 1, nkgdim
      do latIndex = myLatBega, myLatEnda
        do lonIndex = myLonBega, myLonEnda
          pturb_var(lonIndex, latIndex, levIndex) = pturb_var(lonIndex, latIndex, levIndex) / real(NENS, 8)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  
    do latIndex = myLatBega, myLatEnda
      do lonIndex = myLonBega, myLonEnda
        do levIndex = 1, nkgdim
          avg_pturb_var(levIndex) = avg_pturb_var(levIndex) + pturb_var(lonIndex, latIndex, levIndex)
        end do
      end do
    end do
  
    n_grid_point=(lonPerPEa)*(latPerPEa)
    call rpn_comm_allreduce(n_grid_point, n_grid_point_glb, 1,  &
                            "mpi_double_precision", "mpi_sum", "GRID", ierr)
    call rpn_comm_allreduce(avg_pturb_var, avg_pturb_var_glb, nkgdim,  &
                            "mpi_double_precision", "mpi_sum", "GRID", ierr)
  
    !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex)  
    do levIndex = 1, nkgdim
      do latIndex = myLatBega, myLatEnda
        do lonIndex = myLonBega, myLonEnda
          if( pturb_var(lonIndex, latIndex, levIndex) > 0.0d0 ) then
            pturb_var(lonIndex, latIndex, levIndex) = sqrt( pturb_var(lonIndex, latIndex, levIndex))
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  
    do levIndex = 1, nkgdim
      if( avg_pturb_var_glb(levIndex) > 0.0d0 ) then
        avg_pturb_var_glb(levIndex) = sqrt( avg_pturb_var_glb(levIndex) / real(n_grid_point_glb, 8) )
      end if
    end do

    !$OMP PARALLEL DO PRIVATE (memberIndex, levIndex, latIndex, lonIndex)    
    do memberIndex = 1, NENS
      do levIndex = 1, nkgdim
        if( avg_pturb_var_glb(levIndex) > 0.0d0 ) then
          do latIndex = myLatBega, myLatEnda
            do lonIndex = myLonBega, myLonEnda
              if( pturb_var(lonIndex, latIndex, levIndex) > 0.0d0 ) then
                ensemble_r4(lonIndex, latIndex, levIndex, memberIndex) =   &
                    ensemble_r4(lonIndex, latIndex, levIndex, memberIndex) *  &
                    real(avg_pturb_var_glb(levIndex)/pturb_var(lonIndex, latIndex, levIndex), 4)
              end if
            end do
          end do
        end if
      end do
    end do
    !$OMP END PARALLEL DO

  end if

  !- 4.4 Read ensemble mean state
  if ( readEnsMean ) then
    if (mmpi_myid == 0) write(*,*) 'midas-randomPert: reading ensemble mean state'

    call gsv_allocate(stateVectorEnsMean, 1, hco_target, vco_anl, &
                      dateStamp_opt=dateStamp, mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false., &
                      hInterpolateDegree_opt='LINEAR')
    call gio_readFromFile(stateVectorEnsMean, ensMeanFileName, ' ', ' ',  &
                          containsFullField_opt=.true.)

  end if

  !- 4.5 Set perturbations to zero under ice
  if ( setPertZeroUnderIce ) then

    ! read sea-ice analysis
    call gsv_allocate(stateVectorIce, 1, hco_anl, vco_anl, dataKind_opt = 4, &
                      datestamp_opt = -1, mpi_local_opt = .true.,    &
                      varNames_opt = (/'LG'/), hInterpolateDegree_opt ='LINEAR')
    call gio_readFromFile(stateVectorIce, './seaIceAnalysis', ' ','A', &
                          unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVectorIce,seaice_ptr)

    ! set perturbations to zero
    !$OMP PARALLEL DO PRIVATE (memberIndex, levIndex, latIndex, lonIndex)    
    do memberIndex = 1, NENS
      do levIndex = 1, nkgdim
        do latIndex = myLatBega, myLatEnda
          do lonIndex = myLonBega, myLonEnda
            if (seaice_ptr(lonIndex, latIndex, 1) >= iceFractionThreshold) then
              ensemble_r4(lonIndex, latIndex, levIndex, memberIndex) = 0.0
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call gsv_deallocate(stateVectorIce)

  end if

  !- 4.6 Add a fraction of previous date perturbations
  if (previousDateFraction > 0.0) then

    previousDateFactPrev  = previousDateFraction
    ! reduce the amount of random perturbation to maintain steady-state variance
    previousDateFactNoise = sqrt(1.0d0 - previousDateFraction**2)
    if (mmpi_myid == 0) then
      write(*,*) 'midas-randomPert: Scale factor applied to previous date       = ', &
           previousDateFactPrev
      write(*,*) 'midas-randomPert: Scale factor applied to random perturbation = ', &
           previousDateFactNoise
    end if

    ! determine dateStamp of previous date
    call incdatr(dateStampPrevious, dateStamp, -tim_windowsize)
    imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
    ierr    = newdate(dateStampPrevious, datePrint, timePrint, imode)
    datePrevious =  datePrint*100 + timePrint/1000000
    write(datePreviousString, '(I10)') datePrevious
    write(*,*) 'midas-randomPert: previous date, stamp = ', datePrevious, dateStampPrevious

    call gsv_allocate(stateVectorPert, 1, hco_target, vco_anl, &
                      dateStamp_opt=dateStampPrevious, mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false., &
                      hInterpolateDegree_opt='LINEAR')
    call gsv_getField(stateVectorPert,field)

    deallocate(gdmean)
    allocate(gdmean(myLonBegt:myLonEndt, myLatBegt:myLatEndt, nkgdim))
    gdmean(:,:,:) = 0.0D0
    allocate(ensemblePreviousDate_r4(myLonBegt:myLonEndt, myLatBegt:myLatEndt, nkgdim, nEns))

    do memberIndex = 1, NENS

      ! Read previous date perturbations (or perturbed analyses)
      write(memberString, '(I4.4)') memberIndex
      if (readEnsMean) then
        inFileName = './'//trim(datePreviousString)//'_000_'//trim(memberString)
      else
        inFileName = './pert_'//trim(datePreviousString)//'_'//trim(memberString)
      end if
      if( mmpi_myid == 0 ) write(*,*) 'midas-randomPert: reading previous date file = ', inFileName

      call gio_readFromFile(stateVectorPert, inFileName, ' ', ' ',  &
                            containsFullField_opt=.true.)

      ! Copy to big array and accumulate sum for computing mean
      !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex)    
      do levIndex = 1, nkgdim
        do latIndex = myLatBegt, myLatEndt
          do lonIndex = myLonBegt, myLonEndt
            ensemblePreviousDate_r4(lonIndex, latIndex, levIndex, memberIndex) = &
                 real(field(lonIndex, latIndex, levIndex), 4)
            gdmean(lonIndex, latIndex, levIndex) = gdmean(lonIndex, latIndex, levIndex) + &
                                                   field(lonIndex, latIndex, levIndex)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end do

    ! Finish computing mean and remove it
    !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex)    
    do levIndex = 1, nkgdim
      do latIndex = myLatBegt, myLatEndt
        do lonIndex = myLonBegt, myLonEndt
          gdmean(lonIndex, latIndex, levIndex) = gdmean(lonIndex, latIndex, levIndex) / real(NENS, 8)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE (memberIndex, levIndex, latIndex, lonIndex)    
    do memberIndex = 1, NENS
      do levIndex = 1, nkgdim
        do latIndex = myLatBegt, myLatEndt
          do lonIndex = myLonBegt, myLonEndt
            ensemblePreviousDate_r4(lonIndex, latIndex, levIndex, memberIndex) =  &
                ensemblePreviousDate_r4(lonIndex, latIndex, levIndex, memberIndex) -  &
                real(gdmean(lonIndex, latIndex, levIndex), 4)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call gsv_deallocate(stateVectorPert)

  end if

  !- 4.7 Write the perturbations
  do memberIndex = 1, NENS
    if( mmpi_myid == 0 ) write(*,*)
    if( mmpi_myid == 0 ) write(*,*) 'midas-randomPert: pre-processing for writing member number= ', memberIndex

    call gsv_allocate(stateVectorPert, 1, hco_anl, vco_anl, &
                      dateStamp_opt=dateStamp, mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false., &
                      hInterpolateDegree_opt='LINEAR')
    call gsv_getField(stateVectorPert,field)
    call gsv_allocate(stateVectorPertInterp, 1, hco_target, vco_anl, &
                      dateStamp_opt=dateStamp, mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false., &
                      hInterpolateDegree_opt='LINEAR')

    ! Copy mask if it exists
    call gsv_copyMask(stateVectorEnsMean, stateVectorPertInterp)

    !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex)    
    do levIndex = 1, nkgdim
      do latIndex = myLatBega, myLatEnda
        do lonIndex = myLonBega, myLonEnda
          field(lonIndex, latIndex, levIndex) = ensemble_r4(lonIndex, latIndex, levIndex, memberIndex)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    ! interpolate perturbation to the target grid
    call int_interp_gsv(stateVectorPert, stateVectorPertInterp)

    call gsv_getField(stateVectorPertInterp,fieldInterp)

    ! Average current and previous perturbations
    if (previousDateFraction > 0.0) then
      !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex)    
      do levIndex = 1, nkgdim
        do latIndex = myLatBegt, myLatEndt
          do lonIndex = myLonBegt, myLonEndt
            fieldInterp(lonIndex, latIndex, levIndex) =  &
                 previousDateFactPrev * &
                 ensemblePreviousDate_r4(lonIndex, latIndex, levIndex, memberIndex) + &
                 previousDateFactNoise * &
                 fieldInterp(lonIndex, latIndex, levIndex)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end if

    ! add perturbation to supplied ensemble mean
    if ( readEnsMean ) then
      if (mmpi_myid == 0) write(*,*) 'midas-randomPert: adding the ensemble mean and perturbation'
      call gsv_add(stateVectorEnsMean, stateVectorPertInterp)
    end if

    ! determine file name and write to file
    write(memberString, '(I4.4)') memberIndex
    if (readEnsMean) then
      outFileName = './'//trim(dateString)//'_000_'//trim(memberString)
    else
      outFileName = './pert_'//trim(dateString)//'_'//trim(memberString)
    end if

    if( mmpi_myid == 0 ) write(*,*) 'midas-randomPert: processing file = ', outFileName
    stateVectorPertInterp%etiket = 'UNDEFINED'
    call gio_writeToFile(stateVectorPertInterp, outFileName, out_etiket,              & ! IN
                         numBits_opt=numBits, unitConversion_opt=.true., &  ! IN
                         containsFullField_opt=readEnsMean, typvar_opt=typvarOut) 

    call gsv_deallocate(stateVectorPertInterp)
    call gsv_deallocate(stateVectorPert)

  end do

  ! Write out ensemble mean state, if it exists
  if (readEnsMean) then
    ! determine file name and write to file
    memberIndex = 0
    write(memberString, '(I4.4)') memberIndex
    outFileName = './'//trim(dateString)//'_000_'//trim(memberString)

    if( mmpi_myid == 0 ) write(*,*) 'midas-randomPert: processing file = ', outFileName
    stateVectorEnsMean%etiket = 'UNDEFINED'
    call gio_writeToFile(stateVectorEnsMean, outFileName, out_etiket,              & ! IN
                         numBits_opt=numBits, unitConversion_opt=.true., &  ! IN
                         containsFullField_opt=.true., typvar_opt=typvarOut)
  end if

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 5.  Memory deallocations
  !
  deallocate(gdmean)
  deallocate(ensemble_r4)
  if (previousDateFraction > 0.0) then
    deallocate(ensemblePreviousDate_r4)
  end if
  deallocate(controlVector)  
  deallocate(controlVector_mpiglobal)  

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call utl_tmg_stop(0)

  !
  !- 6.  MPI, tmg finalize
  !  
  call tmg_terminate(mmpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 7.  Ending
  !
  if( mmpi_myid == 0 ) write(*,*) '---------------------------------'
  if( mmpi_myid == 0 ) write(*,*) ' MIDAS-RANDOMPERT ENDS'
  if( mmpi_myid == 0 ) write(*,*) '---------------------------------'

end program midas_randomPert

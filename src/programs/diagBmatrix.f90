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

program midas_diagBmatrix
  !
  ! :Purpose: Main program for computing diagnostics of the B and L matrices
  !
  use version_mod
  use midasMpi_mod
  use controlVector_mod
  use gridVariableTransforms_mod
  use varNameList_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use ensembleStateVector_mod
  use bmatrix_mod
  use bmatrixEnsemble_mod
  use localization_mod
  use horizontalCoord_mod
  use advection_mod
  use verticalCoord_mod
  use timeCoord_mod
  use randomNumber_mod
  use utilities_mod
  use ramDisk_mod
  IMPLICIT NONE

  type(struct_gsv) :: statevector, statevectorEnsAmplitude
  type(struct_ens) :: ensAmplitude
  type(struct_hco), pointer :: hco_anl  => null()
  type(struct_hco), pointer :: hco_core => null()
  type(struct_vco), pointer :: vco_anl  => null()
  type(struct_loc), pointer :: loc      => null()
  type(struct_adv), pointer :: adv_amplitudeAssimWindow

  real(8), pointer :: field4d(:,:,:,:)
  real(8), pointer :: field3d(:,:,:)
  real(8), pointer :: ensAmplitude_oneLev(:,:,:,:)
  real(4), allocatable :: randomEns(:,:,:,:)
  real(8), allocatable :: mean(:,:,:)
  real(8), allocatable :: stddev(:,:,:)
  real(8), allocatable :: stddev_zm(:,:),stddev_zm2(:,:)
  real(8), allocatable :: stddev_dm(:,:),stddev_dm2(:,:)
  real(8), allocatable :: zonalMeanStddev(:)
  real(8), allocatable :: controlVector(:), controlVector_global(:)

  real(8) :: centralValue, centralValueLocal

  integer :: fclos, fnom, fstopc, newdate, get_max_rss
  integer :: ierr, nsize, iseed, nulnam, nultxt
  integer :: ensIndex, index, kIndex, nkgdim, levIndex, lonIndex, latIndex
  integer :: dateTime, datePrint, timePrint, dateStamp, numLoc, numStepAmplitude
  integer :: nlevs, nlevs2, varIndex, ip3
  integer :: locIndex, stepIndexInc, nEns, numBensInstance, instanceIndex
  integer :: amp3dStepIndex, nLonLatPos, lonLatPosIndex
  integer :: oneobs_timeStepIndex

  integer, parameter :: lonPosIndex = 1
  integer, parameter :: latPosIndex = 2

  integer :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax
  integer :: myLatBeg, myLatEnd
  integer :: myLonBeg, myLonEnd
  integer, allocatable :: dateStampList(:)

  character(len=128) :: filename, filenameInc, filenameIncNorm, filenameEnsAmp
  character(len=10)  :: datestr
  character(len=4)   :: varName
  character(len=1)   :: locIndexString
  character(len=2)   :: instanceIndexString

  character(len=4), parameter  :: varNameALFAatm(1) = (/ 'ALFA' /)
  character(len=4), parameter  :: varNameALFAsfc(1) = (/ 'ALFS' /)
  character(len=4)             :: varNameALFA(1)

  ! namelist variables
  integer :: numperturbations        ! number of perturbations for randomization estimate of stddev
  integer :: nrandseed               ! initial random seed value
  integer :: oneobs_levs(100)        ! list of level indexes where B matrix columns are computed 
  integer :: oneobs_lonlat(100,2)    ! list of lon,lat index pairs where B matrix columns are computed
  character(len=128) :: oneobs_timeStep ! can be 'first', 'last' or 'middle'
  character(len=4) :: oneobs_varName ! can be 'all' or a specific variable name (default='all')
  logical :: writeEnsAmplitude       ! choose to write ensemble amplitude fields (for ensemble B)
  logical :: writeTextStddev         ! choose to write stddev to text file in addition to standard file
  logical :: writePsiChiStddev       ! choose to also write stddev of Psi/Chi in addition to UU/VV

  namelist /namdiag/numperturbations, nrandseed, oneobs_levs, oneobs_lonlat, &
                    oneobs_varName, oneobs_timeStep, writeEnsAmplitude, writeTextStddev, writePsiChiStddev

  call ver_printNameAndVersion('diagBmatrix','Diagnositcs of the B matrix')

  ! MPI, tmg initialization
  call mmpi_initialize
  call tmg_init(mmpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'Main')
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Set default values for namelist NAMDIAG parameters
  numperturbations  = -1
  nrandseed         =  1
  oneobs_varName    = 'all'
  oneobs_timeStep   = 'middle'
  oneobs_levs(:)    = -1
  oneobs_lonlat(:,:)= -1
  writeEnsAmplitude = .false.
  writeTextStddev   = .false.
  writePsiChiStddev = .false.

  ! Read the parameters from NAMDIAG
  nulnam=0
  ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam,nml=namdiag,iostat=ierr)
  if(ierr.ne.0) call utl_abort('midas-diagBmatrix: Error reading namelist')
  write(*,nml=namdiag)
  ierr=fclos(nulnam)

  nlevs=0
  do index = 1, size(oneobs_levs)
    if (oneobs_levs(index) >= 1) nlevs=nlevs+1
  end do
  nLonLatPos=0
  do index = 1, size(oneobs_lonlat(:,lonPosIndex))
    if (oneobs_lonlat(index,lonPosIndex) >= 1 .and. oneobs_lonlat(index,latPosIndex) >= 1) nLonLatPos=nLonLatPos+1  
  end do

  ! Top Level Control setup
  call ram_setup

  !- Initialize the Temporal grid and set dateStamp from env variable
  call tim_setup()
  if (tim_getDateStamp() == 0) then
    call utl_abort('midas-diagBmatrix: date must be set by env variable MIDAS_DATE')
  end if

  ! Build date-time string from dateStamp
  dateStamp = tim_getDateStamp()
  ierr = newdate(dateStamp,datePrint,timePrint,-3)
  dateTime = datePrint*100 + timePrint/1000000
  write(datestr,'(i10.10)') dateTime
  write(*,*)' datePrint= ',datePrint,' timePrint= ',timePrint
  write(*,*)' date= ',dateTime,' stamp= ',dateStamp

  ! Initialize variables of the model states
  call gsv_setup

  ! Initialize the Analysis horizontal grid
  call hco_SetupFromFile( hco_anl,'./analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    hco_core => hco_anl
  else
    !- Iniatilized the core (Non-Exteded) analysis grid
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
  end if

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Initialize the vertical coordinate from the statistics file
  call vco_SetupFromFile( vco_anl,        & ! OUT
                          './analysisgrid') ! IN

  ! Allocate the statevector
  call gsv_allocate(statevector, tim_nstepobsinc, hco_anl, vco_anl, &
                    datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
                    allocHeight_opt=.false., allocPressure_opt=.false.)
  call gsv_zero(statevector)
  nkgdim = statevector%nk

  ! Setup the B matrix
  call bmat_setup(hco_anl,hco_core,vco_anl)

  !- Initialize the gridded variable transform module
  call gvt_setup(hco_anl,hco_core,vco_anl)
  if ( gsv_varExist(varName='HU') ) call gvt_setupRefFromTrialFiles('HU')
  
  ! Setup of the L matrix done in bmat_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !==============================================
  !- Compute columns of B and L matrices
  !==============================================
  !
  if ( nLevs >= 1 .and. nLonLatPos >= 1 ) then

    !
    !- Compute columns of the B matrix
    !
    select case(trim(oneobs_timeStep))
    case ('first')
      oneobs_timeStepIndex = 1
    case ('middle')
      if (mod(tim_nstepobsinc,2) == 0) then
        write(*,*)
        write(*,*) 'odd number of nstepobsinc a required for obs place in the middle of the analysis window '
        write(*,*) 'tim_nstepobsinc = ', tim_nstepobsinc
        call utl_abort('midas-diagBmatrix')
      end if
      oneobs_timeStepIndex = (tim_nstepobsinc+1)/2
    case ('last')
      oneobs_timeStepIndex = tim_nstepobsinc
    case default
      write(*,*)
      write(*,*) 'Unsupported oneobs_timeStep : ', trim(oneobs_timeStep)
      call utl_abort('midas-diagBmatrix')
    end select

    allocate(controlVector(cvm_nvadim))

    write(*,*) '********************************************'
    write(*,*) 'midas-diagBmatrix: Compute columns of B matrix'
    write(*,*) '********************************************'
    write(*,*)
    write(*,*) ' temporal location          = ',trim(oneobs_timeStep), oneobs_timeStepIndex
    write(*,*) ' number of levels            = ',nLevs
    write(*,*) ' number of lon-lat positions = ', nLonLatPos

    do varIndex = 1, vnl_numvarmax

      if ( .not. gsv_varExist(varName=vnl_varNameList(varIndex)) ) cycle
      if ( trim(oneobs_varName)  /= 'all' .and. (trim(oneobs_varName) /= trim(vnl_varNameList(varIndex))) ) cycle

      filenameInc    = 'columnB_'      // trim(vnl_varNameList(varIndex)) // '_' // datestr // '.fst'
      filenameIncNorm= 'columnBnorm_'  // trim(vnl_varNameList(varIndex)) // '_' // datestr // '.fst'
      filenameEnsAmp = 'ensAmplitude_' // trim(vnl_varNameList(varIndex)) // '_'

      write(*,*)
      write(*,*) 'midas-diagBmatrix: simulating a pseudo-observation of ', trim(vnl_varNameList(varIndex))

      if(vnl_varLevelFromVarname(vnl_varNameList(varIndex)).eq.'SF') then
        nlevs2 = 1
      else
        nlevs2 = nlevs
      end if

      ip3 = 0
      do levIndex = 1, nlevs2
        do lonLatPosIndex = 1, nLonLatPos

          latIndex = oneobs_lonlat(lonLatPosIndex,latPosIndex)
          lonIndex = oneobs_lonlat(lonLatPosIndex,lonPosIndex)

          call gsv_zero(statevector)                   
          call gsv_getField(statevector,field4d,vnl_varNameList(varIndex))

          if ( latIndex >= statevector%myLatBeg .and. latIndex <= statevector%myLatEnd .and. &
               lonIndex >= statevector%myLonBeg .and. lonIndex <= statevector%myLonEnd ) then
            if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'SF') then
              field4d(lonIndex,latIndex,1                    ,oneobs_timeStepIndex) = 1.0D0
            else
              field4d(lonIndex,latIndex,oneobs_levs(levIndex),oneobs_timeStepIndex) = 1.0D0
            end if
          end if

          controlVector(:)=0.0d0
          call bmat_sqrtBT(controlVector,cvm_nvadim,statevector)
          call bmat_sqrtB (controlVector,cvm_nvadim,statevector)
          
          write(*,*)'midas-diagBmatrix: writing out the column of B, levIndex,lonIndex,latIndex=',levIndex,lonIndex,latIndex

          ip3 = ip3 + 1
          
          do stepIndexInc = 1, tim_nstepobsinc
            call gio_writeToFile(statevector,filenameInc,'1OBS_'//trim(vnl_varNameList(varIndex)),  &
                 stepIndex_opt=stepIndexInc, ip3_opt=ip3, unitConversion_opt=.true.)
          end do

          ! Normalized the result to get correlation-like pattern
          centralValueLocal = 0.d0
          if ( latIndex >= statevector%myLatBeg .and. latIndex <= statevector%myLatEnd .and. &
               lonIndex >= statevector%myLonBeg .and. lonIndex <= statevector%myLonEnd ) then
            if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)).eq.'SF') then
              centralValueLocal = field4d(lonIndex,latIndex,1                    ,oneobs_timeStepIndex)
            else
              centralValueLocal = field4d(lonIndex,latIndex,oneobs_levs(levIndex),oneobs_timeStepIndex)
            end if
          end if
          call rpn_comm_allreduce(centralValueLocal, centralValue, 1,  &
                                  "MPI_DOUBLE_PRECISION", "MPI_SUM", "GRID", ierr)
          
          write(*,*) 'midas-diagBmatrix: centralValue found = ', centralValue
          
          if (centralValue /= 0.d0) then
            call gsv_scale(statevector,1.d0/centralValue)
          else
            call utl_abort('midas-diagBmatrix: central value equals 0!')
          end if
          
          do stepIndexInc = 1, tim_nstepobsinc
            call gio_writeToFile(statevector,filenameIncNorm,'1OBSNRM_'//trim(vnl_varNameList(varIndex)), &
                                 stepIndex_opt=stepIndexInc, ip3_opt=ip3,  &
                                 unitConversion_opt=.false.)
          end do

          ! Write the ensemble amplitude fields (i.e., the alphas) when Bens is active
          if (writeEnsAmplitude) call ben_writeAmplitude('./',filenameEnsAmp, ip3) ! IN

        end do
      end do
    end do

    deallocate(controlVector)

    !
    !- Compute columns of the L matrix
    !
    numLoc = ben_getNumLoc()
    numBensInstance = ben_getNumInstance()
    numStepAmplitude = ben_getNumStepAmplitudeAssimWindow()
    amp3dStepIndex   = ben_getAmp3dStepIndexAssimWindow()

    if (numStepAmplitude > 1) then
      allocate(datestampList(numStepAmplitude))
      call tim_getstamplist(dateStampList,numStepAmplitude,tim_getDatestamp())
      nEns = ben_getnEns()
      adv_amplitudeAssimWindow => ben_getAmplitudeAssimWindow()
    else
      allocate(datestampList(1))
      datestampList(1) = dateStamp
    end if

    do instanceIndex = 1, numBensInstance
      do locIndex = 1, numLoc ! (this loop will be done only when localization is used in B)
        loc => ben_getLoc(locIndex,instanceIndex_opt=instanceIndex)

        if (loc%vco%Vcode == 5002 .or. loc%vco%Vcode == 5005) then
          varNameALFA(:) = varNameALFAatm(:)
        else ! vco_anl%Vcode == -1
          varNameALFA(:) = varNameALFAsfc(:)
        end if

        call mmpi_setup_latbands(loc%hco%nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
        call mmpi_setup_lonbands(loc%hco%ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)

        call ens_allocate(ensAmplitude, loc%nEnsOverDimension, numStepAmplitude, loc%hco, loc%vco, &
                          datestampList=dateStampList, varNames_opt=varNameALFA, dataKind_opt=8)

        call gsv_allocate(statevectorEnsAmplitude, numStepAmplitude, loc%hco, loc%vco, &
                          dateStampList_opt=dateStampList, varNames_opt=varNameALFA, dataKind_opt=8, &
                          mpi_local_opt=.true.)

        allocate(controlVector(loc%cvDim))

        write(*,*) '********************************************'
        write(*,*) 'midas-diagBmatrix: Compute columns of L matrix'
        write(*,*) '********************************************'
      
        write(*,*) ' number of levels            = ', nlevs
        write(*,*) ' number of lon-lat positions = ', nLonLatPos

        if (numLoc > 1) then
          write(locIndexString,'(i1)') locIndex
          filename = 'columnL_' // trim(loc%locType) // '_' // locIndexString // '_' // datestr // '.fst'
        else
          if (numBensInstance > 1) then
            write(instanceIndexString,'(I2.2)') instanceIndex
            filename = 'columnL_i' // trim(instanceIndexString) // '_' // trim(loc%locType) // '_' // datestr // '.fst'
          else
            filename = 'columnL_' // trim(loc%locType) // '_' // datestr // '.fst'
          end if
        end if

        ip3 = 0
        do levIndex = 1, nlevs2
          do lonLatPosIndex = 1, nLonLatPos
            
            latIndex = oneobs_lonlat(lonLatPosIndex,latPosIndex)
            lonIndex = oneobs_lonlat(lonLatPosIndex,lonPosIndex)
            
            call ens_zero(ensAmplitude)
            call gsv_zero(statevectorEnsAmplitude)
            if ( latIndex >= statevector%myLatBeg .and. latIndex <= statevector%myLatEnd .and. &
                 lonIndex >= statevector%myLonBeg .and. lonIndex <= statevector%myLonEnd ) then
              ensAmplitude_oneLev => ens_getOneLev_r8(ensAmplitude,oneobs_levs(levIndex))
              ensAmplitude_oneLev(:,oneobs_timeStepIndex,lonIndex,latIndex) = 1.d0
            end if
            controlVector(:)=0.0d0

            if (numStepAmplitude > 1) then
              call adv_ensemble_ad(ensAmplitude,                 & ! INOUT
                                   adv_amplitudeAssimWindow, nEns)  ! IN
            end if

            call loc_LsqrtAd(loc,           & ! IN
                             ensAmplitude,  & ! IN
                             controlVector, & ! OUT
                             amp3dStepIndex)  ! IN
            call loc_Lsqrt  (loc,           & ! IN
                             controlVector, & ! IN
                             ensAmplitude,  & ! OUT
                             amp3dStepIndex)  ! IN

            if (numStepAmplitude > 1) then
              call adv_ensemble_tl(ensAmplitude,                 & ! INOUT
                                   adv_amplitudeAssimWindow, nEns) ! IN
            end if

            write(*,*)'midas-diagBmatrix: writing out the column of L, levIndex,lonIndex,latIndex=',levIndex,lonIndex,latIndex
            call flush(6)

            ip3 = ip3 + 1
            call ens_copyMember(ensAmplitude, statevectorEnsAmplitude, 1)
            do stepIndexInc = 1, numStepAmplitude
              call gio_writeToFile(statevectorEnsAmplitude,filename,'ONEOBS',  &
                                   stepIndex_opt=stepIndexInc, ip3_opt=ip3,unitConversion_opt=.false.)
            end do

          end do
        end do

        deallocate(controlVector)
        call ens_deallocate(ensAmplitude)
        call gsv_deallocate(statevectorEnsAmplitude)

      end do ! localization index (if active in B)
    end do ! Bens instance

  end if ! if any oneobs selected

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !==============================================
  !- Compute the stddev from random perturbations
  !==============================================
  !
  if ( numperturbations > 1 ) then

    write(*,*) '********************************************'
    write(*,*) 'midas-diagBmatrix: Compute the stddev from random perturbations'
    write(*,*) '********************************************'

    allocate(controlVector(cvm_nvadim))
    allocate(controlVector_global(cvm_nvadim_mpiglobal))

    call gsv_zero(statevector)

    ! Allocate the randomEns, mean and stddev
    allocate(randomEns(statevector%myLonBeg:statevector%myLonEnd,statevector%myLatBeg:statevector%myLatEnd,nkgdim,numperturbations))
    allocate(mean(statevector%myLonBeg:statevector%myLonEnd,statevector%myLatBeg:statevector%myLatEnd,nkgdim))
    allocate(stddev(statevector%myLonBeg:statevector%myLonEnd,statevector%myLatBeg:statevector%myLatEnd,nkgdim))

    iseed = abs(nrandseed)
    call rng_setup(iseed)

    call gsv_getField(statevector,field3d)

    !
    !- Compute the ensemble of random perturbations
    !
    do ensIndex = 1, numperturbations
      write(*,*) 'midas-diagBmatrix: computing member number= ',ensIndex
      call flush(6)

      !- Global vector (same for each processors)
      do index = 1, cvm_nvadim_mpiglobal
        controlVector_global(index) = rng_gaussian()
      end do

      !- Extract only the subvector for this processor
      call bmat_reduceToMPILocal(controlVector,       & ! OUT
                                 controlVector_global ) ! IN

      !- Transform to control variables in physical space
      call bmat_sqrtB(controlVector,cvm_nvadim,statevector)

      if ( writePsiChiStddev ) call gvt_transform(statevector,'UVtoPsiChi')

      !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)    
      do kIndex = 1, nkgdim
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            randomEns(lonIndex,latIndex,kIndex,ensIndex) = field3d(lonIndex,latIndex,kIndex)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end do ! Loop on ens member

    !
    !- Compute the ensemble mean
    !
    mean(:,:,:) = 0.0d0
    do ensIndex = 1, numperturbations
      !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)
      do kIndex = 1, nkgdim
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            mean(lonIndex,latIndex,kIndex) = mean(lonIndex,latIndex,kIndex) + randomEns(lonIndex,latIndex,kIndex,ensIndex)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end do

    !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)    
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          mean(lonIndex,latIndex,kIndex) = mean(lonIndex,latIndex,kIndex)/real(numperturbations,8)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !
    !- Remove the ensemble mean from the ensemble
    !
    !$OMP PARALLEL DO PRIVATE (lonIndex,ensIndex,latIndex,kIndex)    
    do ensIndex = 1, numperturbations
      do kIndex = 1, nkgdim
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            randomEns(lonIndex,latIndex,kIndex,ensIndex) = randomEns(lonIndex,latIndex,kIndex,ensIndex) - mean(lonIndex,latIndex,kIndex)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    deallocate(mean)

    !
    !- Compute the ensemble stddev
    !
    stddev(:,:,:) = 0.0d0

    do ensIndex = 1, numperturbations
      !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)
      do lonIndex = statevector%myLonBeg, statevector%myLonEnd
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do kIndex = 1, nkgdim
            stddev(lonIndex,latIndex,kIndex) = stddev(lonIndex,latIndex,kIndex) + &
                 (randomEns(lonIndex,latIndex,kIndex,ensIndex)**2)/real(numperturbations,8)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end do
    deallocate(randomEns)

    !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          stddev(lonIndex,latIndex,kIndex) = sqrt(stddev(lonIndex,latIndex,kIndex))
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !- Insert results in statevector
    !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)    
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          field3d(lonIndex,latIndex,kIndex) = stddev(lonIndex,latIndex,kIndex)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    deallocate(stddev)

    !- Write to file
    call gio_writeToFile(statevector,'stddev_' // datestr // '.fst','GD_STDDEV',  &
                         unitConversion_opt=.true.)

    !
    !- Compute the zonal mean std dev
    !
    write(*,*) 'midas-diagBmatrix: Compute the zonal mean stddev'
    call flush(6)

    allocate(stddev_zm(hco_anl%nj,nkgdim))
    allocate(stddev_zm2(hco_anl%nj,nkgdim))
    stddev_zm(:,:) = 0.d0
    stddev_zm2(:,:) = 0.d0

    !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          stddev_zm(latIndex,kIndex) = stddev_zm(latIndex,kIndex) + (field3d(lonIndex,latIndex,kIndex)**2)/real(hco_anl%ni,8)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    nsize = statevector%nj*nkgdim
    call rpn_comm_allreduce(stddev_zm,stddev_zm2,nsize,  &
         "MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)

    !- Insert results in statevector
    !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)    
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          if(stddev_zm2(latIndex,kIndex)> 0.0d0) then
            field3d(lonIndex,latIndex,kIndex) = sqrt(stddev_zm2(latIndex,kIndex))
          else
            field3d(lonIndex,latIndex,kIndex) = 0.0d0
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    deallocate(stddev_zm)
    deallocate(stddev_zm2)

    call gio_writeToFile(statevector,'stddev_' // datestr // '.fst','ZM_STDDEV',  &
                         unitConversion_opt=.true.)

    ! Write the zonal mean stddev to a text file, if requested
    if ( writeTextStddev ) then
      allocate(zonalMeanStddev(statevector%latPerPE * mmpi_npey))
      do varIndex = 1, vnl_numvarmax
        if (.not. gsv_varExist(varName=vnl_varNameList(varIndex)) ) cycle

        write(*,*) 'midas-diagBmatrix: writing zonal mean stddev to text file for variable: ', vnl_varNameList(varIndex)
        call gsv_getField(statevector,field3d,vnl_varNameList(varIndex))

        varName = vnl_varNameList(varIndex)
        if ( varName == 'HU' ) varName = 'LQ'
        filename = 'stddev_ZM_' // trim(varName) // '_' // datestr // '.txt'
        nultxt = 0
        if ( mmpi_myid == 0 ) ierr = fnom(nultxt,trim(filename),'FTN',0)

        do levIndex = 1, gsv_getNumLevFromVarName(statevector,vnl_varNameList(varIndex))
          nsize = statevector%latPerPE
          call rpn_comm_gather(field3d(1,:,levIndex), nsize, 'mpi_double_precision',  &
                               zonalMeanStddev,     nsize, 'mpi_double_precision', 0, 'NS', ierr )
          if ( mmpi_myid == 0 ) then
            do latIndex = 1, statevector%nj
              write(nultxt,*) field3d(1,latIndex,levIndex)
            end do
          end if
        end do

        if ( mmpi_myid == 0 ) ierr = fclos(nulnam)

      end do
      deallocate(zonalMeanStddev)

    end if

    !
    !- Compute the domain mean std dev
    !
    write(*,*) 'midas-diagBmatrix: Compute the domain mean stddev'
    call flush(6)

    allocate(stddev_dm(hco_anl%ni,nkgdim))
    allocate(stddev_dm2(hco_anl%ni,nkgdim))
    stddev_dm(:,:) = 0.d0
    stddev_dm2(:,:) = 0.d0

    !$OMP PARALLEL DO PRIVATE (latIndex,kIndex)
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          stddev_dm(lonIndex,kIndex) = stddev_dm(lonIndex,kIndex) + (field3d(lonIndex,latIndex,kIndex)**2)/real(hco_anl%nj,8)
        end do
      end do
    end do
    !$OMP END PARALLEL DO 

    nsize = statevector%ni*nkgdim
    call rpn_comm_allreduce(stddev_dm,stddev_dm2,nsize,  &
         "MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)

    !- Insert results in statevector
    !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          if(stddev_dm2(lonIndex,kIndex)> 0.0d0) then
            field3d(lonIndex,latIndex,kIndex) = sqrt(stddev_dm2(lonIndex,kIndex))
          else
            field3d(lonIndex,latIndex,kIndex) = 0.0d0
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    deallocate(stddev_dm)
    deallocate(stddev_dm2)
    deallocate(controlVector)
    deallocate(controlVector_global)

    call gio_writeToFile(statevector,'stddev_' // datestr // '.fst','DM_STDDEV',  &
                         unitConversion_opt=.true.)

  end if ! if numperturbations.gt.1

  call gsv_deallocate(statevector)

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! MPI, tmg finalize
  call utl_tmg_stop(0)
  call tmg_terminate(mmpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr) 

  write(*,*) ' --------------------------------'
  write(*,*) ' midas-diagBmatrix ENDS'
  write(*,*) ' --------------------------------'

end program midas_diagBmatrix

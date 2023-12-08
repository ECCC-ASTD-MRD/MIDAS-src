
module biasCorrectionSat_mod
  ! MODULE biasCorrectionSat_mod (prefix='bcs' category='1. High-level functionality')
  !
  !:Purpose: Performs the bias correction for satellite radiance
  !          data. This includes both the traditional approach based
  !          on regression and the variational bias correction approach
  !          for estimating the bias. Existing bias correction estimates
  !          can also be applied to observations.
  !
  use utilities_mod
  use ramDisk_mod
  use MathPhysConstants_mod
  use obsSpaceData_mod
  use controlVector_mod
  use midasMpi_mod
  use rttov_const, only : ninst
  use tovsNL_mod
  use timeCoord_mod
  use columnData_mod
  use codePrecision_mod
  use localizationFunction_mod
  use HorizontalCoord_mod
  use verticalCoord_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use stateToColumn_mod
  use codtyp_mod
  use timeCoord_mod
  use clibInterfaces_mod
  use obserrors_mod
  use fSQLite
  use obsfiles_mod

  implicit none
  save
  private

  public :: bcs_setup, bcs_calcBias_tl, bcs_calcBias_ad, bcs_writeBias, bcs_finalize
  public :: bcs_removeBiasCorrection, bcs_refreshBiasCorrection
  public :: bcs_do_regression, bcs_filterObs, bcs_computeResidualsStatistics, bcs_calcBias
  public :: bcs_removeOutliers, bcs_applyBiasCorrection
  public :: bcs_mimicSatbcor
  public :: bcs_readConfig
  public :: bcs_dumpBiasToSqliteAfterThinning

  type  :: struct_chaninfo
    integer :: numActivePredictors
    logical :: isDynamic
    integer :: channelNum
    character(len=1) :: bcmode
    character(len=1) :: bctype
    integer, allocatable  :: predictorIndex(:)
    real(8), allocatable  :: coeff(:)
    real(8), allocatable  :: coeffIncr(:)
    real(8), allocatable  :: coeff_fov(:)
    real(8), allocatable  :: coeff_offset(:)
    integer               :: coeff_nobs
    real(8), allocatable  :: coeffIncr_fov(:)
    real(8), allocatable  :: stddev(:)
    real(8), allocatable  :: coeffCov(:,:)
  end type struct_chaninfo

  type  :: struct_bias
    type(struct_chaninfo), allocatable :: chans(:)
    integer :: numscan
    integer :: numChannels
    real(8), allocatable  :: BHalfScanBias(:,:)
    real(8), allocatable  :: BMinusHalfScanBias(:,:)
  end type struct_bias

  type(struct_bias), allocatable  :: bias(:)
  type(struct_vco),       pointer :: vco_mask => null()
  type(struct_hco),       pointer :: hco_mask => null()
  type(struct_columnData) :: column_mask
  logical               :: initialized = .false.
  logical               :: bcs_mimicSatbcor
  logical               :: doRegression
  integer, parameter    :: NumPredictors = 7
  integer, parameter    :: NumPredictorsBcif = 6
  integer, parameter    :: maxfov = 120
  integer, parameter    :: maxNumInst = 25
  integer, parameter    :: maxPassiveChannels = 15
  
  real(8), allocatable  :: trialHeight300m1000(:)
  real(8), allocatable  :: trialHeight50m200(:)
  real(8), allocatable  :: trialHeight1m10(:)
  real(8), allocatable  :: trialHeight5m50(:)
  real(8), allocatable  :: RadiosondeWeight(:)
  real(8), allocatable  :: trialTG(:)
  integer               :: nobs
  integer, external     :: fnom, fclos 
  character(len=2), parameter  :: predTab(0:7) = [ "SB", "KK","T1", "T2", "T3", "T4", "SV", "TG"]
  integer               :: passiveChannelNumber(maxNumInst)
  ! Namelist variables
  character(len=5) :: biasMode  ! "varbc" for varbc, "reg" to compute bias correction coefficients by regression, "apply" to compute and apply bias correction
  logical  :: biasActive        ! logical variable to activate the module
  logical  :: outstats          ! flag to activate output of residual statistics in "reg" mode 
  logical  :: mimicSatbcor      ! in "reg" mode compute regression coefficients the same way as the original satbcor program
  logical  :: weightedestimate  ! flag to activate radiosonde weighting for bias correction computation in "reg" mode
  logical  :: filterObs         ! flag to activate additional observation filtering in "reg" mode. If it is .false. only observations selected for assimilation will be used in the linear regression
  logical  :: removeBiasCorrection  ! flag to activate removal of an already present bias correction
  logical  :: refreshBiasCorrection !flag to replace an existing bias correction with a new one
  logical  :: centerPredictors      ! flag to transparently remove predictor mean in "reg" mode (more stable problem; very little impact on the result)
  logical  :: outCoeffCov           ! flag to activate output of coefficients error covariance (useful for EnKF system)
  logical  :: outOmFPredCov         ! flag to activate output of O-F/predictors coefficients covariances and correlations
  real(8)  :: scanBiasCorLength     ! if positive and .not. mimicSatBcor use error correlation between scan positions with the given correlation length
  real(8)  :: bg_stddev(NumPredictors) ! background error for predictors ("varbc" mode)
  character(len=7) :: cinst(maxNumInst)   ! to read the bcif file for each instrument in cinst
  character(len=3) :: cglobal(maxNumInst) ! a "global" parameter and
  integer          :: nbscan(maxNumInst)  ! the number of scan positions are necessary
  integer          :: passiveChannelList(maxNumInst, maxPassiveChannels)
  ! To understand the meaning of the following parameters controling filtering,
  ! please see  https://wiki.cmc.ec.gc.ca/images/f/f6/Unified_SatRad_Dyn_bcor_v19.pdf pages 20-22
  logical  :: offlineMode   ! flag to select offline mode for bias correction computation
  logical  :: allModeSsmis  ! flag to select "ALL" mode for SSMIS
  logical  :: allModeTovs   ! flag to select "ALL" mode for TOVS (AMSU-A, AMSU-B, MHS, ATMS, MWHS-2)
  logical  :: allModeCsr    ! flag to select "ALL" mode for CSR (GOES, SEVIRI, MVIRI, ABI, etc..)
  logical  :: allModeHyperIr! flag to select "ALL" mode for hyperSpectral Infrared (AIRS, IASI, CrIS)
  logical  :: dumpToSqliteAfterThinning  ! option to output all usefull parameters to sqlite files after thinning
  namelist /nambiassat/ biasActive, biasMode, bg_stddev, removeBiasCorrection, refreshBiasCorrection
  namelist /nambiassat/ centerPredictors, scanBiasCorLength, mimicSatbcor, weightedEstimate
  namelist /nambiassat/ cglobal, cinst, nbscan, passiveChannelList, filterObs, outstats, outCoeffCov
  namelist /nambiassat/ offlineMode, allModeSsmis, allModeTovs, allModeCsr, allModeHyperIr
  namelist /nambiassat/ dumpToSqliteAfterThinning, outOmFPredCov
contains
 
  !-----------------------------------------------------------------------
  ! bcs_readConfig
  !-----------------------------------------------------------------------
  subroutine bcs_readConfig()
    !
    ! :Purpose: Read nambiassat namelist section
    !
    implicit none

    ! Locals:
    integer  :: ierr, nulnam
    logical, save :: firstCall = .true.
    integer :: instrumentIndex, channelIndex
    
    if (.not. firstCall) return
    firstCall = .false.

    ! set default values for namelist variables
    biasActive = .false.
    biasMode = "varbc"
    bg_stddev(:) = 0.0d0
    removeBiasCorrection = .false.
    filterObs = .false.
    refreshBiasCorrection = .false.
    centerPredictors = .false.
    mimicSatbcor = .true.
    scanBiasCorLength = -1.d0
    weightedEstimate = .false.
    outCoeffCov = .false.
    outOmFPredCov = .false.
    nbscan(:) = -1
    cinst(:) = "XXXXXXX"
    cglobal(:) = "XXX"
    outstats = .false.
    offlineMode = .false.
    allModeSsmis = .true.
    allModeTovs = .true.
    allModeCsr = .true.
    allModeHyperIr = .false.
    dumpToSqliteAfterThinning = .false.
    passiveChannelNumber(:) = 0
    passiveChannelList(:,:) = -1
    !
    ! read in the namelist NAMBIASSAT
    if (utl_isNamelistPresent('nambiassat', './flnml')) then
      nulnam = 0
      ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
      read(nulnam,nml=nambiassat,iostat=ierr)
      if (ierr /= 0) call utl_abort('bcs_readConfig: Error reading namelist section nambiassat')
      if (mmpi_myid == 0) write(*,nml=nambiassat)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'bcs_readconfig: nambiassat is missing in the namelist. The default value will be taken.'
    end if

    bcs_mimicSatbcor = mimicSatbcor
    doRegression = (trim(biasMode) == "reg")
    
    do instrumentIndex = 1, maxNumInst
      do channelIndex = 1, maxPassiveChannels
        if (passiveChannelList(instrumentIndex,channelIndex) > 0) then
          passiveChannelNumber(instrumentIndex) = passiveChannelNumber(instrumentIndex) + 1
        end if
      end do
    end do

  end subroutine bcs_readConfig

  !-----------------------------------------------------------------------
  ! bcs_setup
  !-----------------------------------------------------------------------
  subroutine bcs_setup()
    implicit none

    ! Locals:
    integer  :: cvdim
    integer  :: iSensor, iPredictor, instIndex
    integer  :: iChan
    integer  :: iPred, jPred, kPred, iScan1, iScan2
    character(len=85)  :: bcifFile
    character(len=10)  :: instrName, instrNamecoeff, satNamecoeff 
    logical            :: bcifExists
    ! variables from background coeff file
    integer            :: nfov, exitCode
    character(len=2)   :: predBCIF(tvs_maxchannelnumber,numpredictorsBCIF)
    integer            :: canBCIF(tvs_maxchannelnumber), npredBCIF(tvs_maxchannelnumber), ncanBcif, npredictors
    character(len=1)   :: bcmodeBCIF(tvs_maxchannelnumber), bctypeBCIF(tvs_maxchannelnumber)
    character(len=3)   :: global
    character(len=128) :: errorMessage
    real(8), allocatable :: Bmatrix(:,:)

    call bcs_readConfig()

    cvdim = 0

    if (biasActive) then

      if (scanBiasCorLength > 0.d0) call lfn_Setup('FifthOrder')

      allocate(bias(tvs_nSensors))

      do iSensor = 1, tvs_nSensors

        write(*,*) "bcs_setup: iSensor = ", iSensor
       
        instrName = InstrNametoCoeffFileName(tvs_instrumentName(iSensor))
        instrNamecoeff = InstrNameinCoeffFile(tvs_instrumentName(iSensor))
        satNamecoeff = SatNameinCoeffFile(tvs_satelliteName(iSensor)) 

        bcifFile = 'bcif_' // trim(instrName)

        global = "XXX"
        nfov = -1
        do instIndex = 1, size(cinst)
          if (trim(instrNamecoeff) == trim(cinst(instIndex))) then
            global = cglobal(instIndex)
            nfov = nbscan(instIndex)
          end if
        end do
        if (nfov == -1) then
          write(*,*) "bcs_setup: Problem with instrName ", instrNamecoeff
          write(*,'(15(A10,1x))')  cinst(:)
          call utl_abort('bcs_setup: check nambiassat namelist')
        end if

        inquire(file=trim(bcifFile), exist=bcifExists)
        if (bcifExists) then
          !read bcif (Bias Correction Information File)
          call read_bcif(bcifFile, tvs_isNameHyperSpectral(instrName), ncanBcif, &
               canBCIF, bcmodeBCIF, bctypeBCIF, npredBCIF, predBCIF, global, exitcode)
          if (exitcode /= 0) then
            call utl_abort("bcs_setup: Problem in read_bcif while reading " // trim(bcifFile))
          end if

          bias(iSensor)%numChannels = ncanBcif

          allocate(bias(iSensor)%chans(ncanBcif))
          
          do ichan=1, ncanBcif 
            bias(iSensor) % chans(ichan) % channelNum = canBCIF(ichan + 1) 
            bias(iSensor) % chans(ichan) % coeff_nobs = 0
            bias(iSensor) % chans(ichan) % bcmode =  bcmodeBCIF(ichan + 1)
            bias(iSensor) % chans(ichan) % bctype =  bctypeBCIF(ichan + 1)
            bias(iSensor) % chans(ichan) % isDynamic = ( (biasmode == "varbc" .and. bcmodeBCIF(ichan + 1) == "D") .or. biasmode /= "varbc")
            npredictors =  1 + npredBCIF(ichan + 1)
            bias(iSensor) % chans(ichan) % numActivePredictors = npredictors

            allocate( bias(iSensor) % chans(ichan) % stddev( npredictors ) )
            allocate( bias(iSensor) % chans(ichan) % coeffIncr( npredictors ) )
            allocate( bias(iSensor) % chans(ichan) % coeff( npredictors ) )
            allocate( bias(iSensor) % chans(ichan) % coeff_offset( npredictors ) )
            allocate( bias(iSensor) % chans(ichan) % predictorIndex( npredictors ) )
            bias(iSensor) % chans(ichan) % stddev(:) = 0.d0
            bias(iSensor) % chans(ichan) % coeffIncr(:) = 0.d0
            bias(iSensor) % chans(ichan) % coeff(:) = MPC_missingValue_R8
            bias(iSensor) % chans(ichan) % coeff_offset(:) = 0.d0

            bias(iSensor)%chans(ichan)% predictorIndex(1) = 1 !the constant term is always included
            jPred = 1
            do ipred = 1, npredBCIF(ichan + 1)
              jPred =  jPred + 1
              select case(predBCIF(ichan + 1, ipred))
              case('T1')
                kpred = 2
              case('T2')
                kpred = 3
              case('T3')
                kpred = 4
              case('T4')
                kpred = 5
              case('SV')
                kpred = 6
              case('TG')
                kpred = 7
              case default
                write(errorMessage,*) "bcs_setup: Unknown predictor ", predBCIF(ichan+1, ipred), ichan, ipred
                call utl_abort(errorMessage)
              end select
              bias(iSensor)%chans(ichan)%predictorIndex(jPred) = kpred
            end do
          end do
        else
          call utl_abort("bcs_setup: Error : " // trim(bcifFile) // " not present !")
        end if

        bias(iSensor)%numscan = nfov 

        allocate( bias(iSensor) % BHalfScanBias (nfov,nfov))
        if (doRegression) allocate( bias(iSensor) % BMinusHalfScanBias (nfov,nfov))
        allocate( Bmatrix(nfov,nfov))
        do ichan=1, ncanBcif
          allocate( bias(iSensor) % chans(ichan) % coeffIncr_fov(nfov))
          allocate( bias(iSensor) % chans(ichan) % coeff_fov(nfov))
          bias(iSensor) % chans(ichan) % coeffIncr_fov(:) = 0.d0
          bias(iSensor) % chans(ichan) % coeff_fov(:) = MPC_missingValue_R8
        end do

        do ichan = 1, ncanBcif
          if (bias(iSensor)%chans(ichan)%isDynamic) then
            do iPredictor = 1, bias(iSensor)%chans(ichan)%numActivePredictors
              bias(iSensor)%chans(ichan)%stddev(iPredictor) = bg_stddev(bias(iSensor)%chans(ichan)%PredictorIndex(iPredictor))
            end do
          end if
        end do

        if (trim(biasMode) == "varbc") then
          !change dimension of control vector
          do  ichan = 1, ncanBcif
            if (bias(iSensor)%chans(ichan)%isDynamic) &
                 cvdim = cvdim + bias(iSensor)%chans(ichan)%numActivePredictors - 1 + bias(iSensor)%numScan
          end do
        end if

        if (allocated(Bmatrix)) then
          if (scanBiasCorLength > 0.d0) then
            do iScan2 = 1, nfov
              do iScan1 = 1, nfov
                Bmatrix(iScan1,iScan2) = bg_stddev(1) * bg_stddev(1) * lfn_Response(1.d0*abs(iScan1-iScan2), scanBiasCorLength)
              end do
            end do
          else
            Bmatrix(:,:) = 0.d0
            do iScan1 = 1, nfov
              Bmatrix(iScan1,iScan1) = bg_stddev(1) * bg_stddev(1)
            end do
          end if
          bias(iSensor)%BHalfScanBias(:,:) =  Bmatrix(:,:)
          call utl_matsqrt(bias(iSensor)%BHalfScanBias, nfov, 1.d0, printInformation_opt=.true.)
          if (doRegression) then
            bias(iSensor)%BMinusHalfScanBias(:,:) = Bmatrix(:,:)
            call utl_matsqrt(bias(iSensor)%BMinusHalfScanBias, nfov, -1.d0, printInformation_opt=.true.)
          end if
          deallocate(Bmatrix)
        end if

      end do

    end if

    if  (trim(biasMode) == "varbc" .and. cvdim > 0) then
      if (mmpi_myid > 0) cvdim = 0 ! for minimization, all coefficients only on task 0
      call cvm_setupSubVector('BIAS', 'BIAS', cvdim)
    end if

    call  bcs_readCoeffs() ! Read coefficient files in the case of bias correction application (biasMode=="apply")

  end subroutine bcs_setup

  !-----------------------------------------------------------------------
  ! bcs_readCoeffs
  !-----------------------------------------------------------------------
  subroutine bcs_readCoeffs()
    !
    ! :Purpose: Fill the bias structure with read static and dynamic bias correction coefficient files
    !
    implicit none

    ! Locals:
    integer :: iSensor, iSat, jchannel, jChan
    integer :: satIndexDynamic, satIndexStatic
    integer :: chanindexDynamic, chanindexStatic
    character(len=10) :: instrName, instrNamecoeff, satNamecoeff
    character(len=64) :: dynamicCoeffFile, staticCoeffFile
    logical           :: corrected
    integer           :: nfov, npredictors
    character(len=10) :: satsDynamic(tvs_nsensors)       ! satellite names
    integer           :: chansDynamic(tvs_nsensors,tvs_maxchannelnumber)    ! channel numbers
    real(8)           :: fovbiasDynamic(tvs_nsensors,tvs_maxchannelnumber,maxfov)! bias as F(fov)
    real(8)           :: coeffDynamic(tvs_nsensors,tvs_maxchannelnumber,NumPredictors + 1)
    integer           :: nsatDynamic
    integer           :: nchanDynamic(tvs_nsensors)      !number of channels
    integer           :: nfovDynamic
    integer           :: npredDynamic(tvs_nsensors,tvs_maxchannelnumber)    !number of predictors
    character(len=7)  :: cinstrumDynamic      ! string: instrument (e.g. AMSUB) 
    character(len=2)  :: ptypesDynamic(tvs_nsensors,tvs_maxchannelnumber,NumPredictors)
    integer           :: ndataDynamic(tvs_nsensors,tvs_maxchannelnumber)    !number of channels
    character(len=10) :: satsStatic(tvs_nsensors)       ! satellite names
    integer           :: chansStatic(tvs_nsensors,tvs_maxchannelnumber)    !channel numbers 
    real(8)           :: fovbiasStatic(tvs_nsensors,tvs_maxchannelnumber,maxfov)!bias as F(fov)
    real(8)           :: coeffStatic(tvs_nsensors, tvs_maxchannelnumber,NumPredictors + 1)
    integer           :: nsatStatic
    integer           :: nchanStatic(tvs_nsensors)      !number of channels
    integer           :: nfovStatic
    integer           :: npredStatic(tvs_nsensors,tvs_maxchannelnumber)    ! number of predictors
    character(len=7)  :: cinstrumStatic      ! string: instrument (e.g. AMSUB) 9
    character(len=2)  :: ptypesStatic(tvs_nsensors,tvs_maxchannelnumber,NumPredictors)
    integer           :: ndataStatic(tvs_nsensors,tvs_maxchannelnumber)    ! number of channels

    if (biasActive .and. biasMode == "apply") then

      ! 1 fichier de coefficient par intrument avec les differentes plateformes
      ! Cas particulier GEORAD (CSR) 
      
      do iSensor = 1, tvs_nSensors

        write(*,*) "bcs_readCoeffs: iSensor = ", iSensor
       
        instrName = InstrNametoCoeffFileName(tvs_instrumentName(iSensor))
        instrNamecoeff = InstrNameinCoeffFile(tvs_instrumentName(iSensor))
        satNamecoeff = SatNameinCoeffFile(tvs_satelliteName(iSensor)) 

        dynamicCoeffFile = "coeff_file_" // trim(instrName)
        staticCoeffFile = "static_coeff_file_" // trim(instrName)

        if ( tvs_isNameGeostationary(instrName)) then
          dynamicCoeffFile = trim(dynamicCoeffFile) // "." // trim(satNamecoeff)
          staticCoeffFile = trim(staticCoeffFile) // "." // trim(satNamecoeff) 
        end if

        call read_coeff(satsDynamic, chansDynamic, fovbiasDynamic, coeffDynamic, nsatDynamic, nchanDynamic, nfovDynamic, &
             npredDynamic, cinstrumDynamic, dynamicCoeffFile, ptypesDynamic, ndataDynamic)

        call read_coeff(satsStatic, chansStatic, fovbiasStatic, coeffStatic, nsatStatic, nchanStatic, nfovStatic, &
             npredStatic, cinstrumStatic, staticCoeffFile, ptypesStatic, ndataStatic)
        write(*,*) "bcs_readCoeffs: cinstrumDynamic = ", cinstrumDynamic
        write(*,*) "bcs_readCoeffs: cinstrumStatic = ", cinstrumStatic

        satIndexDynamic = -1
        do iSat = 1, nsatDynamic
          if (trim(satNameCoeff) /= trim(satsDynamic(iSat)) .or. trim(instrNamecoeff) /= trim(cinstrumDynamic)) cycle
          satIndexDynamic = iSat
        end do

        satIndexStatic = -1
        do iSat = 1, nsatStatic
          if (trim(satNameCoeff) /= trim(satsStatic(iSat)) .or. trim(instrNamecoeff) /= trim(cinstrumStatic)) cycle
          satIndexStatic = iSat
        end do

        nfov =  bias(iSensor)%numscan

        do jChannel = 1, bias(iSensor)%numChannels

          chanindexDynamic = -1
          if (satIndexDynamic > 0) then
            do jChan = 1, nchanDynamic(satIndexDynamic)
              if (chansDynamic(satIndexDynamic,jChan) == bias(iSensor)%chans(jChannel)%channelNum) then
                chanindexDynamic = jChan
                exit
              end if
            end do
          end if

          chanindexStatic = -1
          if (satIndexStatic > 0) then
            do jChan = 1, nchanStatic(satIndexStatic)
              if (chansStatic(satIndexStatic,jChan) == bias(iSensor)%chans(jChannel)%channelNum) then
                chanindexStatic = jChan
                exit
              end if
            end do
          end if
          npredictors = bias(iSensor)%chans(jChannel)%numActivePredictors
       
          corrected = .true.
          select case(bias(iSensor)%chans(jChannel)%bcmode)
          case("D")
            if (chanindexDynamic > 0) then
              if (bias(iSensor)%chans(jChannel)%bctype=="C") &
                   bias(iSensor)%chans(jChannel)%coeff(1:npredictors) = coeffDynamic(satIndexDynamic,chanindexDynamic,1:npredictors)
              if (bias(iSensor)%chans(jChannel)%bctype=="C" .or. bias(iSensor)%chans(jChannel)%bctype=="F") &
                   bias(iSensor)%chans(jChannel)%coeff_fov(1:nfov) = fovbiasDynamic(satIndexDynamic,chanindexDynamic,1:nfov)
            else if (chanindexStatic > 0) then
              if (bias(iSensor)%chans(jChannel)%bctype=="C") &
                   bias(iSensor)%chans(jChannel)%coeff(1:npredictors) = coeffStatic(satIndexStatic,chanindexStatic,1:npredictors)
              if (bias(iSensor)%chans(jChannel)%bctype=="C" .or. bias(iSensor)%chans(jChannel)%bctype=="F") &
                   bias(iSensor)%chans(jChannel)%coeff_fov(1:nfov) = fovbiasStatic(satIndexStatic,chanindexStatic,1:nfov)
            else
              corrected = .false.
            end if
          case("S")
            if (chanindexStatic > 0) then
              if (bias(iSensor)%chans(jChannel)%bctype=="C") &
                   bias(iSensor)%chans(jChannel)%coeff(1:npredictors) = coeffStatic(satIndexStatic,chanindexStatic,1:npredictors)
              if (bias(iSensor)%chans(jChannel)%bctype=="C" .or. bias(iSensor)%chans(jChannel)%bctype=="F") & 
                   bias(iSensor)%chans(jChannel)%coeff_fov(1:nfov) = fovbiasStatic(satIndexStatic,chanindexStatic,1:nfov)
            else
              corrected = .false.
            end if
          end select
          
          if (.not. corrected) then
            write(*,*) "bcs_readCoeffs: Warning: channel ", bias(iSensor)%chans(jChannel)%channelNum, " of ", &
                 trim(instrName), " ", trim(satNamecoeff), " not corrected!"
          end if

        end do

      end do
        
    end if

  end subroutine bcs_readCoeffs

  !-----------------------------------------------------------------------
  ! bcs_computePredictorBiases
  !-----------------------------------------------------------------------
  subroutine bcs_computePredictorBiases(obsSpaceData)
    !
    ! :Purpose: to compute predictor average
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData

    ! Locals:
    real(8) :: predictor(NumPredictors)
    integer :: iobs, nsize, i, j, npred
    integer :: headerIndex, idatyp, indxtovs
    integer :: iSensor, iFov, iPredictor, ierr
    integer :: bodyIndex, jpred, chanIndx
    real(8), allocatable ::  temp_offset(:,:)
    integer, allocatable ::  temp_nobs(:)
    real(8), allocatable ::  temp_offset2(:,:,:)
    integer, allocatable ::  temp_nobs2(:,:)

    if (centerPredictors) then
      write(*,*) "bcs_computePredictorBiases: start"
      npred = 0
      do iSensor = 1, tvs_nSensors
        npred = max(npred, maxval(bias(iSensor)%chans(:)%numActivePredictors))
      end do
      allocate(temp_offset2(tvs_nsensors,maxval(bias(:)%numChannels),2:npred))
      allocate(temp_nobs2(tvs_nsensors,maxval(bias(:)%numChannels)))
      temp_offset2(:,:,:) = 0.d0
      temp_nobs2(:,:) = 0

      call obs_set_current_header_list(obsSpaceData, 'TO')
      iobs = 0
      HEADER: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER

        ! process only radiance data to be assimilated?
        idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
        if (.not. tvs_isIdBurpTovs(idatyp)) then
          write(*,*) 'bcs_computePredictorBiases: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
          cycle HEADER
        end if
        iobs = iobs + 1
        
        indxTovs = tvs_tovsIndex(headerIndex)
        if (indxTovs < 0) cycle HEADER

        iSensor = tvs_lsensor(indxTovs)

        call obs_set_current_body_list(obsSpaceData, headerIndex)
        iFov = obs_headElem_i(obsSpaceData, OBS_FOV, headerIndex)

        BODY: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY

          if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY   

          call bcs_getChannelIndex(obsSpaceData, iSensor, chanIndx, bodyIndex)
          if (chanindx > 0) then
            call bcs_getPredictors(predictor, headerIndex, iobs, chanIndx, obsSpaceData)
            do iPredictor = 2, bias(iSensor)%chans(chanIndx)%NumActivePredictors
              jPred = bias(iSensor)%chans(chanIndx)%PredictorIndex(iPredictor)
              temp_offset2(iSensor,chanIndx,iPredictor) = temp_offset2(iSensor,chanIndx,iPredictor) + predictor(jPred)
            end do
            temp_nobs2(iSensor,chanIndx) = temp_nobs2(iSensor,chanIndx) + 1
          end if
        end do BODY
      end do HEADER

      allocate(temp_offset(maxval(bias(:)%numChannels),2:npred))
      allocate(temp_nobs(maxval(bias(:)%numChannels)))

      do iSensor = 1, tvs_nSensors 
        temp_offset(:,:) = 0.0d0
        temp_offset(:,:) = temp_offset2(iSensor,:,:)
        call mmpi_allreduce_sumR8_2d( temp_offset, "GRID" )

        do i = 1, bias(iSensor)%numChannels
          do j = 2, bias(iSensor)%chans(i)%numActivePredictors
            bias(iSensor)%chans(i)%coeff_offset(j) = temp_offset(i,j)
          end do
        end do

        temp_nobs(:) = 0
        nsize = size(temp_nobs)
        call rpn_comm_allreduce(temp_nobs2(iSensor,:), temp_nobs, nsize, "mpi_integer", &
             "mpi_sum", "GRID", ierr)
        if (ierr /= 0) then
          call utl_abort("bcs_computePredictorBiases: Erreur de communication MPI 2")
        end if
       
        do i = 1, bias(iSensor)%numChannels
          bias(iSensor)%chans(i)%coeff_nobs = temp_nobs(i)
        end do
        do i = 1, bias(iSensor)%numChannels
          if (bias(iSensor)%chans(i)%coeff_nobs > 0) then
            bias(iSensor)%chans(i)%coeff_offset(:) = bias(iSensor)%chans(i)%coeff_offset / bias(iSensor)%chans(i)%coeff_nobs
          end if
        end do
      end do

      deallocate(temp_offset)
      deallocate(temp_nobs)
      deallocate(temp_offset2)
      deallocate(temp_nobs2)

      write(*,*) "bcs_computePredictorBiases: end"
    end if
   
  end subroutine bcs_computePredictorBiases

  !-----------------------------------------------------------------------
  ! bcs_calcBias
  !-----------------------------------------------------------------------
  subroutine bcs_calcBias(obsSpaceData, columnTrlOnTrlLev)
    !
    ! :Purpose:  to fill OBS_BCOR column of ObsSpaceData body with bias correction computed from read coefficient file
    !
    implicit none

    ! Arguments:
    type(struct_obs),        intent(inout) :: obsSpaceData
    type(struct_columnData), intent(inout) :: columnTrlOnTrlLev

    ! Locals:
    integer  :: headerIndex, bodyIndex, iobs, indxtovs, idatyp
    integer  :: iSensor, iPredictor, chanIndx
    integer  :: iScan, iFov, jPred
    real(8)  :: predictor(NumPredictors)
    real(8)  :: biasCor

    if (.not. biasActive) return

    write(*,*) "bcs_calcBias: start"

    if (.not. allocated(trialHeight300m1000)) then
      call bcs_getTrialPredictors(obsSpaceData, columnTrlOnTrlLev)
    end if

    iobs = 0
    call obs_set_current_header_list(obsSpaceData, 'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if (.not. tvs_isIdBurpTovs(idatyp)) then
        write(*,*) 'bcs_calBias: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER
      end if

      iobs = iobs + 1
      
      indxtovs = tvs_tovsIndex(headerIndex)
      if (indxtovs < 0) cycle HEADER

      iSensor = tvs_lsensor(indxTovs)

      call obs_set_current_body_list(obsSpaceData, headerIndex)
      iFov = obs_headElem_i(obsSpaceData, OBS_FOV, headerIndex)

      if (bias(iSensor)%numScan > 1) then
        iScan = iFov
      else
        iScan = 1
      end if

      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY

        if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY   
        if (obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex) == MPC_missingValue_R8) cycle BODY

        call bcs_getChannelIndex(obsSpaceData, iSensor, chanIndx, bodyIndex)
        if (chanindx > 0) then
          if (bias(iSensor)%chans(chanIndx)%isDynamic .and. bias(iSensor)%numScan >0) then
            if ( bias(iSensor)%chans(chanIndx)%coeff_fov(iScan)/= MPC_missingValue_R8 .and. &
                 all( bias(iSensor)%chans(chanIndx)%coeff(1:bias(iSensor)%chans(chanIndx)%NumActivePredictors)/= MPC_missingValue_R8) ) then
              biasCor = bias(iSensor)%chans(chanIndx)%coeff_fov(iScan) + &
                   bias(iSensor)%chans(chanIndx)%coeff(1)
              call bcs_getPredictors(predictor, headerIndex, iobs, chanIndx, obsSpaceData)
              do iPredictor = 2, bias(iSensor)%chans(chanIndx)%NumActivePredictors
                jPred = bias(iSensor)%chans(chanIndx)%PredictorIndex(iPredictor)
                biasCor = biasCor + predictor(jPred) * bias(iSensor)%chans(chanIndx)%coeff(iPredictor)
              end do
              biasCor = -1.d0 * biascor
              call obs_bodySet_r( obsSpaceData, OBS_BCOR, bodyIndex, biasCor)
            end if
          end if
        end if
      end do BODY
    end do HEADER

    write(*,*) "bcs_calcBias: end"

  end subroutine bcs_calcBias

  !-----------------------------------------------------------------------
  ! bcs_dumpBiasToSqliteAfterThinning
  !-----------------------------------------------------------------------
  subroutine bcs_dumpBiasToSqliteAfterThinning(obsSpaceData)
    !
    ! :Purpose:  to dump bias correction coefficients and predictors in dedicated sqlite files 
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout)     :: obsSpaceData

    ! Locals:
    integer  :: headerIndex, bodyIndex,iobs, indxtovs, idatyp
    integer  :: sensorIndex, iPredictor, chanIndx, codeTypeIndex, fileIndex, searchIndex
    integer  :: iScan, iFov, jPred, burpChanIndex
    real(8)  :: predictor(NumPredictors)
    real(8)  :: biasCor
    real(8)  :: sunzen, sunaz, satzen, sataz
    type(fSQL_DATABASE), allocatable  :: db(:) ! SQLIte  file handle
    type(fSQL_STATEMENT), allocatable :: stmtPreds(:), stmtCoeffs(:) ! precompiled SQLite statements
    type(fSQL_STATUS)    :: stat                 ! SQLiteerror status
    character(len=512)   :: queryCreate, queryPreds, queryCoeffs, queryTrim
    logical, allocatable :: first(:)
    integer, allocatable :: fileIndexes(:), obsOffset(:),  dataOffset(:)
    character(len=*), parameter :: myName = 'bcs_dumpBiasToSqliteAfterThinning::'
    character(len=30)  :: fileNameExtension
    character(len=4)   :: cmyidx, cmyidy
    integer            :: tovsCodeTypeListSize, tovsCodeTypeList(ninst)
    integer            :: tovsFileNameListSize
    character(len=20)  :: tovsFileNameList(30)
    character(len=256) :: fileName
    integer :: tovsAllCodeTypeListSize, tovsAllCodeTypeList(ninst)

    if (.not. biasActive) return
    if (.not. dumpToSqliteAfterThinning) return

    write(*,*) "bcs_dumpBiasToSqliteAfterThinning: start"

    ! get list of all possible tovs codetype values and unique list of corresponding filenames
    call tvs_getAllIdBurpTovs(tovsAllCodeTypeListSize, tovsAllCodeTypeList)
    write(*,*) 'bcs_dumpBiasToSqliteAfterThinning: tovsAllCodeTypeListSize = ', tovsAllCodeTypeListSize
    write(*,*) 'bcs_dumpBiasToSqliteAfterThinning: tovsAllCodeTypeList = ', tovsAllCodeTypeList(1:tovsAllCodeTypeListSize)
    
    tovsFileNameListSize = 0
    tovsFileNameList(:) = 'XXXXX'
    do codeTypeIndex = 1, tovsAllCodeTypeListSize
      fileName = getObsFileName(tovsAllCodeTypeList(codeTypeIndex))
      if (all(tovsFileNameList(:) /= fileName)) then
        tovsFileNameListSize = tovsFileNameListSize + 1
        tovsFileNameList(tovsFileNameListSize) = fileName
      end if
    end do
    write(*,*) 'bcs_dumpBiasToSqliteAfterThinning: tovsFileNameListSize = ', tovsFileNameListSize
    write(*,*) 'bcs_dumpBiasToSqliteAfterThinning: tovsFileNameList = ', tovsFileNameList(1:tovsFileNameListSize)
    
    allocate(db(tovsFileNameListSize))
    allocate(stmtPreds(tovsFileNameListSize))
    allocate(stmtCoeffs(tovsFileNameListSize))
    allocate(first(tovsFileNameListSize))
    first(:) = .true.
    allocate(fileIndexes(size(obsf_fileName)))
    fileIndexes(:) = -1
    do fileIndex = 1, tovsFileNameListSize
      do searchIndex = 1, size(obsf_fileName)
        if (index(trim(obsf_fileName(searchIndex)), trim(tovsFileNameList(fileIndex))) >0) then
          fileIndexes(searchIndex) = fileIndex
        end if
      end do
    end do
    write(*,*) 'bcs_dumpBiasToSqliteAfterThinning: fileIndexes', fileIndexes(1:tovsFileNameListSize)
    allocate(obsOffset(tovsFileNameListSize))
    allocate(dataOffset(tovsFileNameListSize))
    do fileIndex = 1, tovsFileNameListSize
      fileName = tovsFileNameList(fileIndex)
      write(*,*) 'tovs filename = ', fileName
      ! get list of codetypes associated with this filename
      tovsCodeTypeListSize = 0
      tovsCodeTypeList(:) = MPC_missingValue_INT
      do codeTypeIndex = 1, tovsAllCodeTypeListSize
        if (fileName == getObsFileName(tovsAllCodeTypeList(codeTypeIndex))) then
          tovsCodeTypeListSize = tovsCodeTypeListSize + 1
          tovsCodeTypeList(tovsCodeTypeListSize) = tovsAllCodeTypeList(codeTypeIndex)
        end if
      end do
      
      write(*,*) 'bcs_dumpBiasToSqliteAfterThinning: tovsCodeTypeListSize = ', tovsCodeTypeListSize
      write(*,*) 'bcs_dumpBiasToSqliteAfterThinning: tovsCodeTypeList = ', tovsCodeTypeList(1:tovsCodeTypeListSize) 
      call getInitialIdObsData(obsSpaceData, 'TO', obsOffset(fileIndex), dataOffset(fileIndex), &
           codeTypeList_opt=tovsCodeTypeList(1:tovsCodeTypeListSize))
      write(*,*) 'bcs_dumpBiasToSqliteAfterThinning: obsOffset(fileIndex), dataOffset(fileIndex)', fileIndex, obsOffset(fileIndex), dataOffset(fileIndex) 
    end do
    
    iobs = 0
    call obs_set_current_header_list(obsSpaceData, 'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if (.not. tvs_isIdBurpTovs(idatyp)) then
        write(*,*) 'bcs_dumpBiasToSqliteAfterThinning: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER
      end if
      iobs = iobs + 1
      fileIndex = fileIndexes(obs_headElem_i(obsSpaceData, OBS_IDF, headerIndex))
      indxtovs = tvs_tovsIndex(headerIndex)
      if (indxtovs < 0) cycle HEADER
      sensorIndex = tvs_lsensor(indxTovs)
      if (first(fileIndex)) then
        if (obs_mpiLocal(obsSpaceData)) then
          write(cmyidy,'(i4.4)') (mmpi_myidy + 1)
          write(cmyidx,'(i4.4)') (mmpi_myidx + 1)
          fileNameExtension = trim(cmyidx) // '_' // trim(cmyidy)
        else
          fileNameExtension = ' '
        end if

        fileName = 'obs/bcr' // trim(tovsFileNameList(fileIndex)) &
             // '_' // trim(filenameExtension)

        call fSQL_open(db(fileIndex), fileName, stat)
        write(*,*) 'bcs_dumpBiasToSqliteAfterThinning: Open ', trim(fileName)
        if (fSQL_error(stat) /= FSQL_OK) call handleError(stat, 'fSQL_open: ')

        ! Create the tables
        queryCreate = 'CREATE TABLE predictors(id_data integer, id_obs integer, predIndex integer, PredictorValue real,' // &
                      'PredictorType varchar(2), bcor real, fov integer, sunzen real, sunaz real, satzen real, sataz real);' // &
                      'CREATE TABLE  coeffs2(predIndex integer, coeff real, instrument varchar(10), platform varchar(16), vcoord integer, fov integer);'
        call fSQL_do_many(db(fileIndex), queryCreate, stat)
        if (fSQL_error(stat) /= FSQL_OK) call handleError(stat, 'fSQL_do_many: ')

        queryPreds = 'insert into predictors (id_data, id_obs, predIndex, predictorValue,' // &
                     ' predictorType, bcor, fov, sunzen, sunaz, satzen, sataz) values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);'
        queryCoeffs = 'insert into coeffs2 (predIndex, coeff, instrument, platform, vcoord, fov) values(?, ?, ?, ?, ?, ?); '
        write(*,*) myName // ' Insert query Predictors   = ', trim(queryPreds)
        write(*,*) myName // ' Insert query Coeffs = ', trim(queryCoeffs)
        call fSQL_begin(db(fileIndex))
        call fSQL_prepare(db(fileIndex), queryPreds, stmtPreds(fileIndex), stat)
        if (fSQL_error(stat) /= FSQL_OK) call handleError(stat, 'fSQL_prepare 1: ')
        call fSQL_prepare(db(fileIndex), queryCoeffs, stmtCoeffs(fileIndex), stat)
        if (fSQL_error(stat) /= FSQL_OK) call handleError(stat, 'fSQL_prepare 2: ')
        first(fileIndex) = .false.
      end if
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      iFov = obs_headElem_i(obsSpaceData, OBS_FOV, headerIndex)
      sunzen = obs_headElem_r(obsSpaceData, OBS_SUN, headerIndex)
      sunaz = obs_headElem_r(obsSpaceData, OBS_SAZ, headerIndex)
      satzen = obs_headElem_r(obsSpaceData, OBS_SZA, headerIndex)
      sataz = obs_headElem_r(obsSpaceData, OBS_AZA, headerIndex)
      if (bias(sensorIndex)%numScan > 1) then
        iScan = iFov
      else
        iScan = 1
      end if
      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY
        if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY   
        if (obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex) == MPC_missingValue_R8) cycle BODY
        if (btest(obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex), 11)) cycle BODY 
        call bcs_getChannelIndex(obsSpaceData, sensorIndex, chanIndx, bodyIndex)
        if (chanindx > 0) then
          biasCor = 0.0d0
          if (bias(sensorIndex)%chans(chanIndx)%isDynamic .and. bias(sensorIndex)%numScan > 0) then
            call bcs_getPredictors(predictor, headerIndex, iobs, chanIndx, obsSpaceData)
            biasCor = bias(sensorIndex)%chans(chanIndx)%coeff_fov(iScan) + &
                 bias(sensorIndex)%chans(chanIndx)%coeff(1) 
            burpChanIndex = nint(obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex))                   
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=1, int_var=0)
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=2, real8_var=bias(sensorIndex)%chans(chanIndx)%coeff_fov(iScan))
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=3, char_var=trim(tvs_instrumentName(sensorIndex)))
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=4, char_var=trim(tvs_satelliteName(sensorIndex)))
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=5, int_var=burpChanIndex)
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=6, int_var=iScan)
            call fSQL_exec_stmt(stmtCoeffs(fileIndex))
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=1, int_var=1)
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=2, real8_var=bias(sensorIndex)%chans(chanIndx)%coeff(1))
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=3, char_var=trim(tvs_instrumentName(sensorIndex)))
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=4, char_var=trim(tvs_satelliteName(sensorIndex)))
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=5, int_var=burpChanIndex)
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=6, int_var=iScan)
            call fSQL_exec_stmt(stmtCoeffs(fileIndex))
            do iPredictor = 2, bias(sensorIndex)%chans(chanIndx)%NumActivePredictors
              jPred = bias(sensorIndex)%chans(chanIndx)%PredictorIndex(iPredictor)
              biasCor = biasCor + predictor(jPred) * bias(sensorIndex)%chans(chanIndx)%coeff(iPredictor)
            end do
          end if
          biasCor = -1.d0 * biascor
          do iPredictor = 2, bias(sensorIndex)%chans(chanIndx)%NumActivePredictors
            jPred = bias(sensorIndex)%chans(chanIndx)%PredictorIndex(iPredictor)
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=1, int_var=jPred)
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=2, real8_var=bias(sensorIndex)%chans(chanIndx)%coeff(iPredictor))
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=3, char_var=trim(tvs_instrumentName(sensorIndex)))
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=4, char_var=trim(tvs_satelliteName(sensorIndex)))
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=5, int_var=burpChanIndex)
            call fSQL_bind_param(stmtCoeffs(fileIndex), param_index=6, int_var=iScan)
            call fSQL_exec_stmt (stmtCoeffs(fileIndex))
            call fSQL_bind_param(stmtPreds(fileIndex), param_index=1, int_var=bodyIndex + dataOffset(fileIndex))
            call fSQL_bind_param(stmtPreds(fileIndex), param_index=2, int_var=headerIndex + obsOffset(fileIndex))
            call fSQL_bind_param(stmtPreds(fileIndex), param_index=3, int_var=jPred)
            call fSQL_bind_param(stmtPreds(fileIndex), param_index=4, real8_var=predictor(jPred))
            call fSQL_bind_param(stmtPreds(fileIndex), param_index=5, char_var=predtab(jPred))
            call fSQL_bind_param(stmtPreds(fileIndex), param_index=6, real8_var=biascor)
            call fSQL_bind_param(stmtPreds(fileIndex), param_index=7, int_var=iScan)
            call fSQL_bind_param(stmtPreds(fileIndex), param_index=8, real8_var=sunzen)
            call fSQL_bind_param(stmtPreds(fileIndex), param_index=9, real8_var=sunaz )
            call fSQL_bind_param(stmtPreds(fileIndex), param_index=10, real8_var=satzen )
            call fSQL_bind_param(stmtPreds(fileIndex), param_index=11, real8_var=sataz)
            call fSQL_exec_stmt (stmtPreds(fileIndex))
          end do
        end if
      end do BODY
    end do HEADER
    do fileIndex = 1, tovsFileNameListSize
      if (.not. first(fileIndex)) then
        call fSQL_finalize(stmtCoeffs(fileIndex))
        call fSQL_finalize(stmtPreds(fileIndex))
        queryTrim = 'create table coeffs as select distinct * from coeffs2; drop table coeffs2;'
        call fSQL_do_many(db(fileIndex), queryTrim, stat)
        if (fSQL_error(stat) /= FSQL_OK) call handleError(stat, 'fSQL_do_many: ')
        call fSQL_commit(db(fileIndex), stat)
        if (fSQL_error(stat) /= FSQL_OK) call handleError(stat, 'fSQL_commit: ')
        call fSQL_close(db(fileIndex), stat)
        if (fSQL_error(stat) /= FSQL_OK) call handleError(stat, 'fSQL_close: ')
      end if
    end do
    deallocate(dataOffset)
    deallocate(obsOffset)
    deallocate(fileIndexes)
    deallocate(first)
    deallocate(stmtCoeffs)
    deallocate(stmtPreds)
    deallocate(db)
    write(*,*) "bcs_dumpBiasToSqliteAfterThinning: end"
  contains

    subroutine handleError(stat, message)
      implicit none

      ! Arguments:
      type(FSQL_STATUS), intent(in) :: stat
      character(len=*),  intent(in) :: message

      write(*,*) message, fSQL_errmsg(stat)
      call utl_abort(trim(message))
    end subroutine handleError

  end subroutine bcs_dumpBiasToSqliteAfterThinning

  !---------------------------------------
  ! bcs_computeResidualsStatistics
  !----------------------------------------
  subroutine bcs_computeResidualsStatistics(obsSpaceData, prefix)
    !
    ! :Purpose: to compute residuals mean and standard deviation by intrument, channel and scan position
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout)  :: obsSpaceData
    character(len=*), intent(in)     :: prefix

    ! Locals:
    real(8), allocatable :: tbias(:,:), tstd(:,:)
    integer, allocatable :: tcount(:,:)
    real(8), allocatable :: biasMpiGlobal(:,:), stdMpiGLobal(:,:)
    integer, allocatable :: countMpiGlobal(:,:)
    integer :: sensorIndex, headerIndex, bodyIndex
    integer :: nchans, nscan
    integer :: iSensor, iScan, chanIndx, iFov
    real(8) :: OmF, bcor
    integer :: ierr, nulfile1, nulfile2
    character(len=10) :: instrName, satNamecoeff
    character(len=72) :: errorMessage

    if (.not. biasActive) return

    if (.not. outstats) return

    write(*,*) "bcs_computeResidualsStatistics: start"

    SENSORS:do sensorIndex = 1, tvs_nsensors

      if  (.not. tvs_isReallyPresentMpiGLobal(sensorIndex)) cycle SENSORS

      write(*,*) "bcs_computeResidualsStatistics: sensorIndex ", sensorIndex

      nchans = bias(sensorIndex)%numChannels
      nscan = bias(sensorIndex)%numscan
     
      allocate(tbias(nchans,nscan))
      tbias(:,:) = 0.d0
      allocate(tstd(nchans,nscan))
      tstd(:,:) = 0.d0
      allocate(tcount(nchans,nscan))
      tcount(:,:) = 0

      call obs_set_current_header_list(obsSpaceData, 'TO')

      HEADER: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER
        if (tvs_tovsIndex(headerIndex) < 0) cycle HEADER
        iSensor = tvs_lsensor(tvs_tovsIndex(headerIndex))
        if (iSensor /= sensorIndex) cycle HEADER
          
        iFov = obs_headElem_i(obsSpaceData, OBS_FOV, headerIndex)
        if (nscan > 1) then
          iScan = iFov
        else
          iScan = 1
        end if

        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY
          
          if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY 
          call bcs_getChannelIndex(obsSpaceData, iSensor, chanIndx, bodyIndex)
          if (chanindx > 0) then
            OmF = obs_bodyElem_r(obsSpaceData, OBS_OMP, bodyIndex)
            bcor =  obs_bodyElem_r(obsSpaceData, OBS_BCOR, bodyIndex)
            if (OmF /= MPC_missingValue_R8 .and. bcor /= MPC_missingValue_R8) then
              tbias(chanIndx,iScan) = tbias(chanIndx,iScan) + (OmF + bcor)
              tstd(chanIndx,iScan) = tstd(chanIndx,iScan) + (OmF + bcor) ** 2
              tcount(chanIndx,iScan) = tcount(chanIndx,iScan) + 1
            end if
          end if
        end do BODY
      end do HEADER

      allocate( biasMpiGlobal(nchans,nscan) )
      allocate( stdMpiGlobal(nchans,nscan) )
      allocate( countMpiGlobal(nchans,nscan) )

      call mmpi_reduce_sumR8_2d( tbias, biasMpiGlobal, 0, "GRID" )
      call mmpi_reduce_sumR8_2d( tstd, stdMpiGlobal, 0, "GRID" )
      call rpn_comm_reduce(tcount, countMpiGlobal, size(countMpiGlobal), "mpi_integer", "MPI_SUM", 0, "GRID", ierr)
      if (ierr /=0) then
        write(errorMessage,*) "bcs_computeResidualsStatistics: MPI communication error 3", ierr 
        call utl_abort(errorMessage)
      end if

      if (mmpi_myId == 0) then
        where(countMpiGlobal > 0) 
          biasMpiGlobal = biasMpiGlobal / countMpiGlobal
          stdMpiGlobal = sqrt(stdMpiGlobal/ countMpiGlobal  - biasMpiGlobal**2)
        end where
      
        instrName = InstrNametoCoeffFileName(tvs_instrumentName(sensorIndex))
        satNamecoeff = SatNameinCoeffFile(tvs_satelliteName(sensorIndex)) 

        nulfile1 = 0
        ierr = fnom(nulfile1, './std_' // trim(instrName) // '_' // trim(satNamecoeff) // trim(prefix) // '.dat', 'FTN+FMT', 0)
        nulfile2 = 0
        ierr = fnom(nulfile2, './mean_' // trim(instrName) // '_' // trim(satNamecoeff) // trim(prefix) // '.dat', 'FTN+FMT', 0)

        do chanIndx = 1, nchans
          if (sum(countMpiGlobal(chanIndx, :)) > 0) then
            write(nulfile2,'(i4,1x,100e14.6)') chanIndx, biasMpiGlobal(chanindx,:)
            write(nulfile1,'(i4,1x,100e14.6)') chanIndx, stdMpiGlobal(chanindx,:)
          end if
        end do
        ierr = fclos(nulfile1)
        ierr = fclos(nulfile2)
       
      end if

      deallocate(biasMpiGlobal)
      deallocate(stdMpiGlobal)
      deallocate(countMpiGlobal)
      deallocate(tbias)
      deallocate(tstd)
      deallocate(tcount)
 
    end do SENSORS

    write(*,*) "bcs_computeResidualsStatistics: end"
    
  end subroutine bcs_computeResidualsStatistics

  !---------------------------------------
  !  bcs_removeOutliers
  !---------------------------------------- 
  subroutine  bcs_removeOutliers(obsSpaceData)
    !
    ! :Purpose: to remove outliers (too large OmF) from linear regression
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout)  :: obsSpaceData

    ! Locals:
    real(8), allocatable :: tbias(:,:), tstd(:,:)
    integer, allocatable :: tcount(:,:)
    real(8), allocatable :: biasMpiGlobal(:,:), stdMpiGLobal(:,:)
    integer, allocatable :: countMpiGlobal(:,:)
    integer :: sensorIndex, headerIndex, bodyIndex
    integer :: nchans, nfiles
    integer :: iSensor, chanIndx, timeIndex
    real(8) :: OmF
    real(8), parameter :: alpha = 5.d0
    real(8) :: stepObsIndex
    integer :: ierr
    character(len=72) :: errorMessage

    if (.not. biasActive) return

    write(*,*) "bcs_removeOutliers: start"

    SENSORS:do sensorIndex = 1, tvs_nsensors

      if  (.not. tvs_isReallyPresentMpiGLobal(sensorIndex)) cycle SENSORS

      write(*,*) "bcs_removeOutliers: sensorIndex ", sensorIndex

      nchans = bias(sensorIndex)%numChannels
      nfiles = tim_nstepobs
     
      allocate(tbias(nchans,nfiles))
      tbias(:,:) = 0.d0
      allocate(tstd(nchans,nfiles))
      tstd(:,:) = 0.d0
      allocate(tcount(nchans,nfiles))
      tcount(:,:) = 0

      call obs_set_current_header_list(obsSpaceData, 'TO')

      HEADER1: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER1
        if (tvs_tovsIndex(headerIndex) < 0) cycle HEADER1
        iSensor = tvs_lsensor(tvs_tovsIndex(headerIndex))
        if (iSensor /= sensorIndex) cycle HEADER1

        call tim_getStepObsIndex(stepObsIndex, tim_getDatestamp(), &
             obs_headElem_i(obsSpaceData, OBS_DAT, headerIndex), &
             obs_headElem_i(obsSpaceData, OBS_ETM, headerIndex), tim_nstepobs)

        timeIndex = nint(stepObsIndex)
        if  (timeIndex < 0) cycle HEADER1
        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY1: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY1
          
          if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated .and. &
               .not. btest(obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex), 6)) then 
            call bcs_getChannelIndex(obsSpaceData, iSensor, chanIndx, bodyIndex)
            if (chanindx > 0) then
              OmF = obs_bodyElem_r(obsSpaceData, OBS_OMP, bodyIndex)
              if (OmF /= MPC_missingValue_R8) then
                tbias(chanIndx,timeindex) = tbias(chanIndx,timeindex) + OmF
                tstd(chanIndx,timeindex) = tstd(chanIndx,timeindex) + (OmF ** 2)
                tcount(chanIndx,timeindex) =  tcount(chanIndx,timeindex) + 1
              end if
            end if
          end if
        end do BODY1
      end do HEADER1

      allocate(biasMpiGlobal(nchans,nfiles))
      allocate(stdMpiGlobal(nchans,nfiles))
      allocate(countMpiGlobal(nchans,nfiles))

      call mmpi_reduce_sumR8_2d( tbias, biasMpiGlobal, 0, "GRID" )
      call mmpi_reduce_sumR8_2d( tstd, stdMpiGlobal, 0, "GRID" )
      call rpn_comm_reduce(tcount, countMpiGlobal, size(countMpiGlobal), "mpi_integer", "MPI_SUM", 0, "GRID", ierr)
      if (ierr /=0) then
        write(errorMessage,*) "bcs_removeOutliers: MPI communication error 3", ierr 
        call utl_abort(errorMessage)
      end if

      if (mmpi_myId == 0) then
        where(countMpiGlobal > 0) 
          biasMpiGlobal = biasMpiGlobal / countMpiGlobal
          stdMpiGlobal = sqrt(stdMpiGlobal/ countMpiGlobal  - biasMpiGlobal**2)
        end where
      end if

      call rpn_comm_bcast(countMpiGlobal, nchans*nfiles, "mpi_integer", 0, "GRID", ierr)
      if (ierr /=0) then
        write(errorMessage,*) "bcs_removeOutliers: MPI communication error 4", ierr 
        call utl_abort(errorMessage)
      end if
      call rpn_comm_bcast(stdMpiGlobal, nchans*nfiles, "mpi_double_precision", 0, "GRID", ierr)
      if (ierr /=0) then
        write(errorMessage,*) "bcs_removeOutliers: MPI communication error 5", ierr 
        call utl_abort(errorMessage)
      end if

      if (sum(countMpiGlobal) /= 0) then

        tcount(:,:) = 0

        call obs_set_current_header_list(obsSpaceData, 'TO')

        HEADER2: do
          headerIndex = obs_getHeaderIndex(obsSpaceData)
          if (headerIndex < 0) exit HEADER2
          if (tvs_tovsIndex(headerIndex) < 0) cycle HEADER2
          iSensor = tvs_lsensor(tvs_tovsIndex(headerIndex))
          if (iSensor /= sensorIndex) cycle HEADER2

          call tim_getStepObsIndex(stepObsIndex, tim_getDatestamp(), &
               obs_headElem_i(obsSpaceData, OBS_DAT, headerIndex), &
               obs_headElem_i(obsSpaceData, OBS_ETM, headerIndex), tim_nstepobs)

          timeIndex = nint(stepObsIndex)
          if  (timeIndex < 0) cycle HEADER2
          
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY2: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY2
          
            if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated .and. &
                 .not. btest(obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex), 6)) then 

              call bcs_getChannelIndex(obsSpaceData, iSensor, chanIndx, bodyIndex)

              if (chanindx > 0) then
             
                OmF = obs_bodyElem_r(obsSpaceData, OBS_OMP, bodyIndex)
                if (countMpiGlobal(chanIndx, timeindex) > 2 .and.  &
                     abs(OmF) > alpha * stdMpiGlobal(chanIndx, timeindex)) then
                  call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
                  tcount(chanIndx,timeindex) = tcount(chanIndx,timeindex) + 1
                end if
             
              end if
            end if
          end do BODY2
        end do HEADER2

        do chanIndx = 1, nchans
          if (sum(tcount(chanIndx,:)) > 0) then
            write(*,'(A,1x,i2,1x,i4,1x,i6,1x,30e14.6)') "bcs_removeOutliers:", sensorindex, chanIndx, sum(tcount(chanIndx,:)), stdMpiGlobal(chanIndx,:)
          end if
        end do

      end if

      deallocate(biasMpiGlobal)
      deallocate(stdMpiGlobal)
      deallocate(countMpiGlobal)
      deallocate(tbias)
      deallocate(tstd)
      deallocate(tcount)

    end do SENSORS

    write(*,*) "bcs_removeOutliers: end"

  end subroutine bcs_removeOutliers

  !---------------------------------------
  ! bcs_calcBias_tl
  !---------------------------------------- 
  subroutine bcs_calcBias_tl(cv_in, obsColumnIndex, obsSpaceData, columnTrlOnTrlLev)
    !
    ! :Purpose: tl of bias computation (for varBC)
    !
    implicit none

    ! Arguments:
    real(8),                 intent(in)    :: cv_in(:)
    integer,                 intent(in)    :: obsColumnIndex
    type(struct_obs),        intent(inout) :: obsSpaceData
    type(struct_columnData), intent(inout) :: columnTrlOnTrlLev

    ! Locals:
    integer  :: headerIndex, bodyIndex, iobs, indxtovs, idatyp
    integer  :: iSensor, iPredictor, chanIndx
    integer  :: iScan, iFov, jPred
    real(8)  :: predictor(NumPredictors)
    real(8), pointer :: cv_bias(:)
    real(8), target  :: dummy4Pointer(1)
    real(8)  :: biasCor

    if (.not. biasActive) return

    if (.not. allocated(trialHeight300m1000)) then
      call bcs_getTrialPredictors(obsSpaceData, columnTrlOnTrlLev)
      call bcs_computePredictorBiases(obsSpaceData)
      call bcs_getRadiosondeWeight(obsSpaceData, lmodify_obserror_opt=.true.)
    end if

    nullify(cv_bias)
    if (mmpi_myid == 0) then
      if (cvm_subVectorExists('BIAS')) then
        cv_Bias => cvm_getSubVector(cv_in, 'BIAS')
        write(*,*) 'bcs_calcBias_tl: maxval(cv_bias) = ', maxval(cv_bias(:))
      else
        write(*,*) 'bcs_calcBias_tl: control vector does not include bias coefficients'
        return
      end if
   else
      cv_bias => dummy4Pointer
   end if

    ! get bias coefficients
    call bcs_cvToCoeff(cv_bias)

    ! apply bias increment to specified obs column
    iobs = 0
    call obs_set_current_header_list(obsSpaceData, 'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if (.not. tvs_isIdBurpTovs(idatyp)) then
        write(*,*) 'bcs_calBias_tl: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER
      end if
      
      iobs = iobs + 1

      indxtovs = tvs_tovsIndex(headerIndex)
      if (indxtovs < 0) cycle HEADER

      iSensor = tvs_lsensor(indxTovs)

      call obs_set_current_body_list(obsSpaceData, headerIndex)
      iFov = obs_headElem_i(obsSpaceData, OBS_FOV, headerIndex)

      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY

        if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY   

        call bcs_getChannelIndex(obsSpaceData, iSensor, chanIndx, bodyIndex)

        if (chanindx > 0) then
          biasCor = 0.0d0
          if (bias(iSensor)%chans(chanIndx)%isDynamic .and. bias(iSensor)%numScan >0) then
            call bcs_getPredictors(predictor, headerIndex, iobs, chanindx, obsSpaceData)
            do iPredictor = 1, bias(iSensor)%chans(chanIndx)%NumActivePredictors
              jPred = bias(iSensor)%chans(chanIndx)%PredictorIndex(iPredictor)
              if (iPredictor == 1) then
                if (bias(iSensor)%numScan > 1) then
                  iScan = iFov
                else
                  iScan = 1
                end if
                biasCor = biasCor + predictor(jPred) * bias(iSensor)%chans(chanIndx)%coeffIncr_fov(iScan) 
              else
                biasCor = biasCor + predictor(jPred) * bias(iSensor)%chans(chanIndx)%coeffIncr(iPredictor) 
              end if
            end do
          end if
          call obs_bodySet_r(obsSpaceData, obsColumnIndex, bodyIndex, &
               obs_bodyElem_r(obsSpaceData, obsColumnIndex, bodyIndex) - biasCor)
        end if
      end do BODY
    end do HEADER


  end subroutine bcs_calcBias_tl
 
  !----------------------
  ! bcs_getTrialPredictors
  !----------------------
  subroutine bcs_getTrialPredictors(obsSpaceData, columnTrlOnTrlLev)
    !
    ! :Purpose: get predictors from trial fields
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: columnTrlOnTrlLev
    type(struct_obs),        intent(inout) :: obsSpaceData

    ! Locals:
    integer  :: headerIndex, idatyp, iobs
    real(8)  :: height1, height2

    if (tvs_nobtov > 0) then
      allocate(trialHeight300m1000(tvs_nobtov))
      allocate(trialHeight50m200(tvs_nobtov))
      allocate(trialHeight5m50(tvs_nobtov))
      allocate(trialHeight1m10(tvs_nobtov))
      allocate(trialTG(tvs_nobtov))
      allocate(RadiosondeWeight(tvs_nobtov))
    else
      write(*,*) 'bcs_getTrialPredictors: No radiance OBS found'
      return
    end if
    
    iobs = 0

    call obs_set_current_header_list(obsSpaceData, 'TO')

    HEADER2: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER2
      idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if (.not.  tvs_isIdBurpTovs(idatyp)) then
        write(*,*) 'bcs_getTrialPredictors: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER2
      end if
      iobs = iobs + 1

      height1 = logInterpHeight(columnTrlOnTrlLev, headerIndex, 1000.d0)
      height2 = logInterpHeight(columnTrlOnTrlLev, headerIndex, 300.d0)
      
      trialHeight300m1000(iobs) = height2 - height1

      height1 = logInterpHeight(columnTrlOnTrlLev, headerIndex, 200.d0)
      height2 = logInterpHeight(columnTrlOnTrlLev, headerIndex, 50.d0)

      trialHeight50m200(iobs) = height2 - height1

      height1 = height2
      height2 = logInterpHeight(columnTrlOnTrlLev, headerIndex, 5.d0)

      trialHeight5m50(iobs) = height2 - height1

      height1 = logInterpHeight(columnTrlOnTrlLev, headerIndex, 10.d0)
      height2 = logInterpHeight(columnTrlOnTrlLev, headerIndex, 1.d0)

      trialHeight1m10(iobs) = height2 - height1

      trialTG(iobs) = col_getElem(columnTrlOnTrlLev, 1, headerIndex, 'TG')

    end do HEADER2

    if (trialTG(1) > 150.0d0) then
      write(*,*) 'bcs_getTrialPredictors: converting TG from Kelvin to deg_C'
      trialTG(:) = trialTG(:) - MPC_K_C_DEGREE_OFFSET_R8
    end if

    trialHeight300m1000(:) = 0.1d0 * trialHeight300m1000(:) ! conversion factor
    trialHeight50m200(:) = 0.1d0 * trialHeight50m200(:)
    trialHeight5m50(:) = 0.1d0 * trialHeight5m50(:)
    trialHeight1m10(:) =  0.1d0 *  trialHeight1m10(:)

    write(*,*) 'bcs_getTrialPredictors: end'

  contains

    function logInterpHeight(columnTrlOnTrlLev, headerIndex, p) result(height)
      implicit none

      ! Arguments:
      type(struct_columnData), intent(inout) :: columnTrlOnTrlLev
      integer,                 intent(in)    :: headerIndex
      real(8),                 intent(in)    :: p
      ! Result:
      real(8) :: height

      ! Locals:
      integer :: jk, nlev, ik
      real(8) :: zpt, zpb, zwt, zwb
      real(8), pointer :: col_ptr(:)

      ik = 1
      nlev = col_getNumLev(columnTrlOnTrlLev, 'TH')
      do jk = 2, nlev - 1
        zpt = col_getPressure(columnTrlOnTrlLev, jk, headerIndex, 'TH') * MPC_MBAR_PER_PA_R8
        if(p > zpt) ik = jk
      end do
      zpt = col_getPressure(columnTrlOnTrlLev, ik, headerIndex, 'TH') * MPC_MBAR_PER_PA_R8
      zpb = col_getPressure(columnTrlOnTrlLev, ik + 1, headerIndex, 'TH') * MPC_MBAR_PER_PA_R8

      zwb = log(p/zpt) / log(zpb/zpt)
      zwt = 1.d0 - zwb
      col_ptr => col_getColumn(columnTrlOnTrlLev, headerIndex, 'Z_T')

      height = zwb * col_ptr(ik+1) + zwt * col_ptr(ik)
   
    end function logInterpHeight

  end subroutine bcs_getTrialPredictors

  !----------------------
  ! bcs_cvToCoeff
  !----------------------
  subroutine bcs_cvToCoeff(cv_bias)
    !
    ! :Purpose: get coefficient increment from control vector
    !
    implicit none

    ! Arguments:
    real(8), intent(in) :: cv_bias(:)

    ! Locals:
    integer  :: index_cv, iSensor, iChannel, iPredictor, iScan
    integer  :: nsize, ierr
 
    if (mmpi_myid == 0) then
      write(*,*) 'bcs_cvToCoeff: start'
      index_cv = 0
      ! initialize of coeffIncr
      do iSensor = 1, tvs_nSensors
        if (bias(iSensor)%numScan > 0) then
          do iChannel = 1, bias(iSensor)%numChannels
            if (bias(iSensor)%chans(iChannel)%isDynamic) then
              bias(iSensor)%chans(iChannel)%coeffIncr(:) = 0.0d0
              bias(iSensor)%chans(iChannel)%coeffIncr_fov(:) = 0.0d0
            end if
          end do
        end if
      end do

      do iSensor = 1, tvs_nSensors
        if (bias(iSensor)%numScan > 0) then
          do iChannel = 1, bias(iSensor)%numChannels
            if (bias(iSensor)%chans(iChannel)%isDynamic) then
              do iPredictor = 1, bias(iSensor)%chans(iChannel)%numActivePredictors
                if (iPredictor == 1) then
                  do iScan = 1, bias(iSensor)%numScan
                    index_cv = index_cv + 1
                    bias(iSensor)%chans(iChannel)%coeffIncr_fov(iScan) = bias(iSensor)%chans(iChannel)%stddev(iPredictor) * cv_bias(index_cv)
                  end do
                else
                  index_cv = index_cv + 1
                  bias(iSensor)%chans(iChannel)%coeffIncr(iPredictor) = bias(iSensor)%chans(iChannel)%stddev(iPredictor) * cv_bias(index_cv)
                end if
              end do !iPredictor
            end if ! isDynamic
          end do !iChannel
        end if
      end do !iSensor

    end if
    
    ! for constant part
    do iSensor = 1, tvs_nSensors
      if (bias(iSensor)%numScan > 0) then
        nsize = bias(iSensor)%numScan 
        do iChannel = 1, bias(iSensor)%numChannels
          if (bias(iSensor)%chans(iChannel)%isDynamic) &
               call rpn_comm_bcast(bias(iSensor)%chans(iChannel)%coeffIncr_fov, nsize, "mpi_double_precision", 0, "GRID", ierr)
        end do
      end if
    end do

    ! for predictor part
    do iSensor = 1, tvs_nSensors
      if (bias(iSensor)%numScan > 0) then
        do iChannel = 1, bias(iSensor)%numChannels
          if (bias(iSensor)%chans(iChannel)%isDynamic) then
            nsize = bias(iSensor)%chans(iChannel)%numActivePredictors 
            call rpn_comm_bcast(bias(iSensor)%chans(iChannel)%coeffIncr, nsize, "mpi_double_precision", 0, "GRID", ierr)
          end if
        end do
      end if
    end do

  end subroutine bcs_cvToCoeff

  !-----------------------------------
  ! bcs_getPredictors
  !---------------------------------- 
  subroutine bcs_getPredictors(predictor, headerIndex, obsIndex, chanindx, obsSpaceData)
    !
    ! :Purpose: get predictors
    !
    implicit none

    ! Arguments:
    real(8),          intent(out)   :: predictor(NumPredictors)
    integer,          intent(in)    :: headerIndex
    integer,          intent(in)    :: obsIndex
    integer,          intent(in)    :: chanindx
    type(struct_obs), intent(inout) :: obsSpaceData

    ! Locals:
    integer  :: iSensor, iPredictor, jPredictor
    real(8)  :: zenithAngle

    predictor(:) = 0.0d0
    
    do iPredictor = 1, NumPredictors

      if (iPredictor == 1) then
        ! constant
        predictor(iPredictor) = 1.0d0
      else if (iPredictor == 2) then
        ! Height300-Height1000 (dam) /1000 T1
        predictor(iPredictor) = trialHeight300m1000(obsIndex) / 1000.0d0
      else if (iPredictor == 3) then
        ! Height50-Height200 (dam) /1000   T2
        predictor(iPredictor) = trialHeight50m200(obsIndex) / 1000.0d0
      else if (iPredictor == 4) then
        ! Height5-Height50 (dam) /1000    T3
        predictor(iPredictor) = trialHeight5m50(obsIndex) / 1000.0d0
      else if (iPredictor == 5) then
        ! Height1-Height10 (dam) /1000    T4
        predictor(iPredictor) = trialHeight1m10(obsIndex) / 1000.0d0        
      else if (iPredictor == 6) then
        ! SV secant of satellite zenith angle minus one
        zenithAngle = obs_headElem_r(obsSpaceData, OBS_SZA, headerIndex) 
        if (zenithAngle < 75.) predictor(iPredictor) = (1.d0 / cos(zenithAngle * MPC_RADIANS_PER_DEGREE_R8)) - 1.d0
      else if (iPredictor == 7) then
        ! skin temperature (C) /10
        predictor(iPredictor) = trialTG(obsIndex)
      end if

    end do

    iSensor = tvs_lsensor(tvs_tovsIndex(headerIndex))
    do  iPredictor = 1, bias(iSensor)%chans(chanIndx)%numActivePredictors
      jPredictor = bias(iSensor)%chans(chanIndx)%predictorIndex(iPredictor)
      predictor(jPredictor) = predictor(jPredictor) - bias(iSensor)%chans(chanindx)%coeff_offset(iPredictor)
    end do
   
  end subroutine bcs_getPredictors

  !---------------------------------------------
  ! bcs_calcBias_ad
  !---------------------------------------------
  subroutine bcs_calcBias_ad(cv_out, obsColumnIndex, obsSpaceData)
    !
    ! :Purpose: bias computation adjoint (for varBC)
    !
    implicit none

    ! Arguments:
    real(8),          intent(in)    :: cv_out(:)
    integer,          intent(in)    :: obsColumnIndex
    type(struct_obs), intent(inout) :: obsSpaceData

    ! Locals:
    integer  :: headerIndex, bodyIndex, iobs, idatyp
    integer  :: iSensor, iChannel, iPredictor, chanIndx
    integer  :: iScan, iFOV, jPred
    real(8)  :: predictor(NumPredictors)
    real(8), pointer  :: cv_bias(:)
    real(8), target  :: dummy4Pointer(1)
    real(8)  :: biasCor

    if (.not. biasActive) return

    if (mmpi_myid == 0) write(*,*) 'bcs_calcBias_ad: start'

    nullify(cv_bias)
    if (mmpi_myid == 0) then
      if (cvm_subVectorExists('BIAS')) then
        cv_bias => cvm_getSubVector(cv_out, 'BIAS')
      else
        write(*,*) 'bcs_calcBias_ad: control vector does not include bias coefficients'
        return
      end if
    else
      cv_bias => dummy4Pointer
    end if

    ! adjoint of applying bias increment to specified obs column
    do iSensor = 1, tvs_nSensors
      if (bias(iSensor)%numScan > 0) then
        do iChannel = 1, bias(iSensor)%numChannels
          if (bias(iSensor)%chans(iChannel)%isDynamic) then
            bias(iSensor)%chans(iChannel)%coeffIncr(:) = 0.0d0
            bias(iSensor)%chans(iChannel)%coeffIncr_fov(:) = 0.0d0
          end if
        end do
      end if
    end do

    iobs = 0
    call obs_set_current_header_list(obsSpaceData, 'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if (.not.  tvs_isIdBurpTovs(idatyp)) then
        write(*,*) 'bcs_calcBias_ad: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER
      end if
      iobs = iobs + 1
      if (tvs_tovsIndex(headerIndex) < 0) cycle HEADER

      iSensor = tvs_lsensor(tvs_tovsIndex(headerIndex))

      call obs_set_current_body_list(obsSpaceData, headerIndex)
      iFov = obs_headElem_i(obsSpaceData, OBS_FOV, headerIndex)

      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY

        if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY

        call bcs_getChannelIndex(obsSpaceData, iSensor, chanIndx, bodyIndex)

        if (chanindx > 0) then
          if (bias(iSensor)%chans(chanIndx)%isDynamic) then
            call bcs_getPredictors(predictor, headerIndex, iobs, chanIndx, obsSpaceData)
            biasCor = obs_bodyElem_r(obsSpaceData, obsColumnIndex, bodyIndex)
            do iPredictor = 1, bias(iSensor)%chans(chanIndx)%numActivePredictors
              jPred = bias(iSensor)%chans(chanIndx)%PredictorIndex(iPredictor)
              if (jPred == 1) then
                if (bias(iSensor)%numScan > 1) then
                  iScan = iFov
                else
                  iScan = 1
                end if
                bias(iSensor)%chans(chanIndx)%coeffIncr_fov(iScan) = bias(iSensor)%chans(chanIndx)%coeffIncr_fov(iScan) &
                     + predictor(jPred) * biasCor
              else
                bias(iSensor)%chans(chanIndx)%coeffIncr(iPredictor) = bias(iSensor)%chans(chanIndx)%coeffIncr(iPredictor) & 
                     + predictor(jPred) * biasCor
              end if
            end do !iPredictor
          end if
        end if

      end do BODY

    end do HEADER

    ! put the coefficients into the control vector
    call bcs_cvToCoeff_ad(cv_bias)

    if (mmpi_myid == 0) then
      write(*,*) 'bcs_calcBias_ad: maxval(cv_bias) = ', maxval(cv_bias(:))
    end if

  end subroutine bcs_calcBias_ad

  !----------------------------------------------------
  ! bcs_cvToCoeff_ad
  !----------------------------------------------------
  subroutine bcs_cvToCoeff_ad(cv_bias)
    !
    ! :Purpose: adjoint of control vector to coeff transfer (for varBC)
    !
    implicit none

    ! Arguments:
    real(8), intent(inout)  :: cv_bias(:)

    ! Locals:
    integer  :: index_cv, iSensor, iChannel, iPredictor, iScan
    integer  :: nChan, nScan
    integer  :: nsize, iChan
    real(8), allocatable  :: temp_coeffIncr(:), temp_coeffIncr_fov(:)

    write(*,*) "bcs_cvToCoeff_ad: start"

    do iSensor = 1, tvs_nSensors
      if (bias(iSensor)%numScan > 0) then
        nChan = bias(iSensor)%numChannels
        do ichan = 1, nChan
          if (bias(iSensor)%chans(ichan)%isDynamic) then
            nSize = bias(iSensor)%chans(iChan)%numActivePredictors
            allocate(temp_coeffIncr(nSize))
            temp_coeffIncr(:) = 0.0d0
            call mmpi_reduce_sumR8_1d( bias(iSensor)%chans(ichan)%coeffIncr(:), temp_coeffIncr, 0, "GRID" )
            bias(iSensor)%chans(ichan)%coeffIncr(:) = temp_coeffIncr(:)
            deallocate(temp_coeffIncr)
          end if
        end do
      end if
    end do

    do iSensor = 1, tvs_nSensors
      if (bias(iSensor)%numScan > 0) then
        nChan = bias(iSensor)%numChannels
        nScan = bias(iSensor)%numScan
        nsize = nChan * nScan
        if (nsize > 0) then
          allocate(temp_coeffIncr_fov(1:nScan))
          do ichan = 1, nChan
            if (bias(iSensor)%chans(ichan)%isDynamic) then
              temp_coeffIncr_fov(:) = 0.0d0
              call mmpi_reduce_sumR8_1d( bias(iSensor)%chans(ichan)%coeffIncr_fov, temp_coeffIncr_fov, 0, "GRID" )
              bias(iSensor)%chans(iChan)%coeffIncr_fov(:) = temp_coeffIncr_fov(:)
            end if
          end do
          deallocate(temp_coeffIncr_fov)
        end if
      end if
    end do

    if (mmpi_myid == 0) then
      cv_bias(:) = 0.d0
      index_cv = 0
      do iSensor = 1, tvs_nSensors
        if (bias(iSensor)%numScan > 0) then
          do iChannel = 1, bias(iSensor)%numChannels
            if (bias(iSensor)%chans(iChannel)%isDynamic) then
              do iPredictor = 1, bias(iSensor)%chans(iChannel)%numActivePredictors
                if (iPredictor == 1) then
                  do iScan = 1, bias(iSensor)%numScan
                    index_cv = index_cv + 1
                    cv_bias(index_cv) = bias(iSensor)%chans(iChannel)%stddev(iPredictor) * bias(iSensor)%chans(iChannel)%coeffIncr_fov(iScan)
                  end do
                else
                  index_cv = index_cv + 1
                  cv_bias(index_cv) = bias(iSensor)%chans(iChannel)%stddev(iPredictor) * bias(iSensor)%chans(iChannel)%coeffIncr(iPredictor)
                end if
              end do
            end if
          end do
        end if
      end do
    end if
    
  end subroutine bcs_cvToCoeff_ad

  !-----------------------------------------
  ! bcs_writeBias
  !-----------------------------------------
  subroutine bcs_writeBias(cv_in_opt)
    !
    ! :Purpose: to write bias increments and coefficients (varBC)
    !
    implicit none

    ! Arguments:
    real(8), optional, intent(in)  :: cv_in_opt(:)

    ! Locals:
    integer  :: iSensor, iChannel, iPredictor, iScan
    integer  :: jSensor, iChannel2
    integer  :: nulfile_inc, nulfile_fov, ierr
    real(8), pointer :: cv_bias(:)
    real(8), target  :: dummy4Pointer(1)
    character(len=80) :: BgFileName
    !for background coeff and write out
    integer             :: iInstr
    integer             :: numCoefFile, jCoef, kCoef
    character(len=10)   :: coefInstrName(tvs_nSensors), temp_instrName
    character(len=25)   :: filecoeff
    logical             :: coeffExists
    ! these variables are not used but need to be present to satisfy bcs_updateCoeff interface
    ! some bcs_updateCoeff arguments could be made optional (todo)
    integer            :: chans(tvs_nSensors, tvs_maxChannelNumber), nsat, nfov
    integer            :: nchan(tvs_nSensors)
    character(len=10)  :: sats(tvs_nsensors) ! satellite names
    character(len=7)   :: cinstrum           ! string: instrument (e.g. AMSUB)

    if (.not. biasActive) return

    if (present(cv_in_opt)) then
      nullify(cv_bias)
      if (mmpi_myid == 0) then
        if (cvm_subVectorExists('BIAS')) then
          cv_bias => cvm_getSubVector(cv_in_opt, 'BIAS')
          write(*,*) 'bcs_writeBias: maxval(cv_bias) = ', maxval(cv_bias(:))
        else
          write(*,*) 'bcs_writeBias: control vector does not include bias coefficients'
          return
        end if
      else
        cv_bias => dummy4Pointer
      end if
      call bcs_cvToCoeff(cv_bias)
    end if

    ! apply transformation to account for predictor offset

    do iSensor = 1, tvs_nSensors
      if (bias(iSensor)%numScan > 0) then
        do iChannel = 1, bias(iSensor)%numChannels
          do iPredictor = 2, bias(iSensor)%chans(iChannel)%numActivePredictors
            if (mimicSatbcor) then
              bias(iSensor)%chans(iChannel)%coeffIncr(1) = bias(iSensor)%chans(iChannel)%coeffIncr(1) - &
                   bias(iSensor)%chans(iChannel)%coeffIncr(iPredictor) * bias(iSensor)%chans(iChannel)%coeff_offset(iPredictor)
            else
              do iScan = 1, bias(iSensor)%numScan
                bias(iSensor)%chans(iChannel)%coeffIncr_fov(iScan) = bias(iSensor)%chans(iChannel)%coeffIncr_fov(iScan) - &
                     bias(iSensor)%chans(iChannel)%coeffIncr(iPredictor) * bias(iSensor)%chans(iChannel)%coeff_offset(iPredictor)
              end do
            end if
          end do
          if (bias(iSensor)%chans(iChannel)%numActivePredictors > 0 .and. mmpi_myId ==0 .and. (.not. doRegression)) &
               write(*,*) "bcs_writeBias: bias(iSensor)%chans(iChannel)%coeffIncr(:) = ",  bias(iSensor)%chans(iChannel)%coeffIncr(:)
        end do
      end if
    end do

    if (doRegression) then
      call bcs_writeCoeff()
      return
    end if

    if (mmpi_myId == 0) then

      ! write out bias coefficient increments in ascii file
      nulfile_inc = 0
      ierr = fnom(nulfile_inc, './satbias_increment.dat', 'FTN+FMT', 0)

      do iSensor = 1, tvs_nSensors
        write(nulfile_inc,'(/,1X,"Sensor Index=",I3,", Satellite Name=",A15,", Instrument Name=",A15)') &
             iSensor, tvs_satelliteName(iSensor), tvs_instrumentName(iSensor)
        if (bias(iSensor)%numScan > 0) then
          do iChannel = 1, bias(iSensor)%numChannels
            if (bias(iSensor)%chans(iChannel)%isDynamic) then
              iChannel2 = bias(iSensor)%chans(iChannel)%channelNum
              if (sum(bias(iSensor)%chans(iChannel)%coeffIncr(:)) /= 0.0d0) &
                   write(nulfile_inc,'(3X,"Channel number=",I4)') iChannel2
              do iPredictor = 2, bias(iSensor)%chans(iChannel)%numActivePredictors
                if (bias(iSensor)%chans(iChannel)%coeffIncr(iPredictor) /= 0.0d0) &
                     write(nulfile_inc,'(5X,"Predictor number=",I4,", Coefficient=",e12.4)') &
                     iPredictor, bias(iSensor)%chans(iChannel)%coeffIncr(iPredictor)
              end do
            end if
          end do
        end if
      end do

      ierr = fclos(nulfile_inc)

      ! write out fovbias coefficient increments in ascii file
      nulfile_fov = 0
      ierr = fnom(nulfile_fov, './fovbias_incre.dat', 'FTN+FMT', 0)
      do iSensor = 1, tvs_nSensors
        write(nulfile_fov,'(/,1X,"Sensor Index=",I3,", Satellite Name=",A15,", Instrument Name=",A15)') &
             iSensor, tvs_satelliteName(iSensor), tvs_instrumentName(iSensor)
        if (bias(iSensor)%numScan > 0) then
          do iChannel = 1, bias(iSensor)%numChannels
            if (bias(iSensor)%chans(iChannel)%isDynamic) then
              iChannel2 = bias(iSensor)%chans(iChannel)%channelNum
              if (sum(bias(iSensor)%chans(iChannel)%coeffIncr_fov(:)) /= 0.0d0) then
                write(nulfile_fov,'(3X,"Channel number=",I4)') iChannel2 
                write(nulfile_fov,*) bias(iSensor)%chans(iChannel)%coeffIncr_fov(:)
              end if
            end if
          end do
        end if
      end do
      ierr = fclos(nulfile_fov)

    end if

    ! Find the background coeff_file number and name
    do iSensor = 1, tvs_nSensors
      numCoefFile = 0
      jCoef = 0
      do jSensor = 1, tvs_nSensors
        temp_instrName = InstrNametoCoeffFileName(tvs_instrumentName(jSensor))
        filecoeff = 'coeff_file_' // trim(temp_instrName) // ''
        inquire(file=trim(filecoeff), exist=coeffExists)

        if (coeffExists) then
          numCoefFile = numCoefFile + 1
          jCoef = jCoef + 1
          coefInstrName(jCoef) = temp_instrName
        end if
        if (jSensor > 1) then
          do kCoef = 1, jCoef - 1
            if (temp_instrName == coefInstrName(kCoef)) then
              numCoefFile = numCoefFile - 1
              jCoef = jCoef - 1
            end if
          end do
        end if 
      end do
    end do

    ! update coeff_file_instrument and write out
    do iInstr=1, numCoefFile 
      BgFileName = './coeff_file_' // coefInstrName(iInstr)
      call bcs_updateCoeff(tvs_nSensors, NumPredictors, BgFileName, sats, chans, nsat, nchan, nfov, cinstrum)
    end do

  end subroutine bcs_writeBias

  !-----------------------------------------
  ! bcs_updateCoeff
  !-----------------------------------------
  subroutine bcs_updateCoeff(maxsat, maxpred, coeff_file, sats, chans, nsat, nchan, nfov, cinstrum, updateCoeff_opt)
    !
    ! :Purpose: to read, and optionaly update and write out, the coeff files (varBC).
    !
    implicit none

    ! Arguments:
    integer,           intent(in)  :: maxsat, maxpred
    character(len=*),  intent(in)  :: coeff_file
    logical, optional, intent(in)  :: updateCoeff_opt
    integer,           intent(out) :: chans(maxsat,tvs_maxChannelNumber) ! channel numbers
    integer,           intent(out) :: nsat
    integer,           intent(out) :: nfov
    integer,           intent(out) :: nchan(maxsat) ! number of channels
    character(len=10), intent(out) :: sats(maxsat)  ! Satellite names
    character(len=*),  intent(out) :: cinstrum      ! instrument (e.g. AMSUB)

    ! Locals:
    real(8)            :: fovbias(maxsat,tvs_maxChannelNumber,maxfov)
    real(8)            :: coeff(maxsat,tvs_maxChannelNumber,maxpred)
    character(len=2) :: ptypes(maxsat,tvs_maxChannelNumber,maxpred) 
    integer            :: npred(maxsat,tvs_maxChannelNumber)           ! number of predictors
    integer            :: ndata(maxsat,tvs_maxChannelNumber)
    ! LOCAL for reading background coeff file
    integer            :: iSat, jChan, kPred, kFov
    logical            :: verbose
    ! update coeff files
    real               :: fovbias_an(maxsat,tvs_maxChannelNumber,maxfov)
    real               :: coeff_an(maxsat,tvs_maxChannelNumber,maxpred) 
    integer            :: iSensor, jChannel, iFov, iPred, totPred
    character(len=10)  :: tmp_SatName, tmp_InstName 
    ! write out files 
    integer            :: iuncoef2, ierr, numPred
    character(len=80):: filename2
    logical            :: updateCoeff_opt2
    !   sats(nsat)            = satellite names
    !   chans(nsat, nchan(i))  = channel numbers of each channel of each satellite i
    !   npred(nsat, nchan(i))  = number of predictors for each channel of each satellite i
    !   fovbias(i, j, k)        = bias for satellite i, channel j, FOV k   k=1,nfov
    !     if FOV not considered for instrument, nfov = 1 and fovbias is global bias for channel
    !   coeff(i, j, 1)          = regression constant
    !   coeff(i, j, 2), ..., coeff(i, j, npred(i, j)) = predictor coefficients
    !   nsat, nchan, nfov, cinstrum (output) are determined from file
    !   if returned nsat = 0, coeff_file was empty
    !   maxpred (input) is max number of predictors
    !   maxsat (input)  is max number of satellites

    ! There are three parts in this subroutine, read, update and write out the coeff files
    ! 
    !- 1. read in the background coeff files, this program is read_coeff from genbiascorr
    ! 
    if (present(updateCoeff_opt)) then
      updateCoeff_opt2 = updateCoeff_opt
    else
      updateCoeff_opt2 = .true.
    end if

    verbose = .false.
   
    call read_coeff(sats, chans, fovbias, coeff, nsat, nchan, nfov, &
         npred, cinstrum, coeff_file, ptypes, ndata)

    ! Transfer of coefficients read from file to the bias structure
    satloop :do iSat = 1, nsat  !backgroud sat
      instloop:do iSensor = 1, tvs_nSensors
        ! for Satellite Name
        tmp_SatName = SatNameinCoeffFile(tvs_satelliteName(iSensor))
        ! for Instrument Name
        tmp_InstName = InstrNameinCoeffFile(tvs_instrumentName(iSensor))
        
        if (trim(tmp_SatName) /= trim(sats(iSat)) .or. trim(tmp_InstName) /= trim(cinstrum)) cycle instloop
        write(*,*) "bcs_updateCoeff: " // tmp_SatName // " " // tmp_InstName

        if (.not. allocated(bias(iSensor)%chans)) cycle instloop

        chan1loop:do jChan = 1, nchan(iSat)
          chan2loop:do jChannel = 1, bias(iSensor)%numChannels  

            if (chans(iSat, jChan) /= bias(iSensor)%chans(jChannel)%channelNum) cycle chan2loop

            ! part 1 for coeffIncr
            do iFov = 1, nfov
              bias(iSensor)%chans(jchannel)%coeff_fov(iFov) = fovbias(iSat, jChan, iFov)
            end do ! iFov
            
            ! part 2 for coeffIncr_fov
            totPred  = bias(iSensor)%chans(jchannel)%NumActivePredictors 
            do iPred = 1, totPred
              bias(iSensor)%chans(jchannel)%coeff(iPred) = coeff(iSat, jChan, iPred)
            end do ! iPred
            
          end do chan2loop ! jChannel
        end do chan1loop !jChan
        
      end do instloop
    end do satloop
    
    if (.not. updateCoeff_opt2) return 
    
    !
    !- 2.update coeff and fovbias  
    !
    coeff_an(:,:,:) = coeff(:,:,:)
    fovbias_an(:,:,:) = fovbias(:,:,:)

    do iSat = 1, nsat  !backgroud sat
      do iSensor = 1, tvs_nSensors
        ! for Satellite Name
        tmp_SatName = SatNameinCoeffFile(tvs_satelliteName(iSensor))
        ! for Instrument Name
        tmp_InstName = InstrNameinCoeffFile(tvs_instrumentName(iSensor))
        if (trim(tmp_SatName) /= trim(sats(iSat)) .or. trim(tmp_InstName) /= trim(cinstrum)) cycle 
        do jChan = 1, nchan(iSat)
          do jChannel = 1, bias(iSensor)%numChannels  

            if (chans(iSat, jChan) /= bias(iSensor)%chans(jChannel)%channelNum) cycle

            ! part 1 for coeffIncr
            do iFov = 1, nfov
              fovbias_an(iSat, jChan, iFov) = fovbias(iSat, jChan, iFov) + bias(iSensor)%chans(jchannel)%coeffIncr_fov(iFov)
            end do ! iFov

            ! part 2 for coeffIncr_fov
            totPred  = bias(iSensor)%chans(jchannel)%NumActivePredictors 
            do iPred = 1, totPred
              coeff_an(iSat, jChan, iPred) = coeff(iSat, jChan, iPred) + bias(iSensor)%chans(jchannel)%coeffIncr(iPred)
            end do ! iPred

          end do ! jChannel
        end do !jChan
      end do !iSensor
    end do ! iSat

    !
    !- 3. Write out updated_coeff
    ! 

    if (mmpi_myId == 0) then

      iuncoef2 = 0
      filename2 = './anlcoeffs_' // cinstrum 
      ierr = fnom(iuncoef2, filename2, 'FTN+FMT', 0)

      write(*,*) 'bcs_updateCoeff: write in bcs_updateCoeff'
   
      do iSat = 1, nsat
        do jChan = 1, nchan(iSat)
          numPred = npred(iSat, jChan)
          if (sum(abs(coeff_an(iSat, jchan, 1:numpred+1))) /= 0.d0 .and. sum(abs(fovbias_an(iSat, jChan, 1:nfov))) /= 0.d0) then 
            write(iuncoef2,'(A52,A8,1X,A7,1X,I6,1X,I8,1X,I2,1X,I3)') &
                 'SATELLITE, INSTRUMENT, CHANNEL, NOBS, NPRED, NSCAN: ', sats(iSat), cinstrum, chans(iSat,jChan), ndata(isat,jchan), numPred, nfov
            write(iuncoef2,'(A7,6(1X,A2))') 'PTYPES:', (ptypes(iSat, jChan, kPred), kPred = 1, numPred)
            write(iuncoef2,'(120(1x,ES17.10))') (fovbias_an(iSat, jChan, kFov), kFov = 1, nfov)
            write(iuncoef2,*) (coeff_an(iSat, jChan, kPred), kPred = 1, numPred + 1)
          end if
        end do
      end do

      ierr = fclos(iuncoef2) 
   
      write(*,*) 'bcs_updateCoeff: finish writing coeffient file', filename2
    
    end if

  end subroutine bcs_updateCoeff

  !-----------------------------------------
  ! bcs_writeCoeff
  !-----------------------------------------
  subroutine bcs_writeCoeff()
    !
    ! :Purpose:  write out  the coeff files (regression case).
    !
    implicit none

    ! Locals:
    integer            :: iuncoef, numPred, ierr
    character(len=80)  :: filename
    character(len=80)  :: instrName, satNamecoeff
    integer :: sensorIndex, nchans, nscan, nfov, kpred, kFov, jChan

    if (mmpi_myId == 0) then

      SENSORS:do sensorIndex = 1, tvs_nsensors

        if  (.not. tvs_isReallyPresentMpiGLobal(sensorIndex)) cycle SENSORS

        write(*,*) "bcs_writeCoeff: sensorIndex ", sensorIndex

        nchans = bias(sensorIndex)%numChannels
        nscan = bias(sensorIndex)%numscan

        instrName = InstrNameinCoeffFile(tvs_instrumentName(sensorIndex))
        satNamecoeff = SatNameinCoeffFile(tvs_satelliteName(sensorIndex)) 

        filename = './anlcoeffs_' // trim(instrName)  
        call utl_open_asciifile(filename, iuncoef)
        nfov = bias(sensorIndex)%numScan
        do jChan = 1, nchans
          if (bias(sensorIndex)%chans(jChan)%coeff_nobs > 0) then
            numPred = bias(sensorIndex)%chans(jChan)%numActivePredictors 
          
            write(iuncoef,'(A52,A8,1X,A7,1X,I6,1X,I8,1X,I2,1X,I3)') 'SATELLITE, INSTRUMENT, CHANNEL, NOBS, NPRED, NSCAN: ',  &
                 satNameCoeff, instrName, bias(sensorIndex)%chans(jChan)%channelNum, bias(sensorIndex)%chans(jChan)%coeff_nobs, numPred - 1, nfov
            write(iuncoef,'(A7,6(1X,A2))') 'PTYPES:',  (predtab(bias(sensorIndex)%chans(jChan)%predictorIndex(kPred)), kPred = 2, numPred)
            write(iuncoef,'(120(1x,ES17.10))') (bias(sensorIndex)%chans(jChan)%coeff_fov(kFov), kFov = 1, nfov)
            write(iuncoef,'(12(1x,ES17.10))') (bias(sensorIndex)%chans(jChan)%coeff(kPred), kPred = 1, numPred)
          end if
        end do

        ierr = fclos(iuncoef) 

        if (outCoeffCov) then
          filename = './anlcoeffsCov_' // trim(instrName)  
          call utl_open_asciifile(filename, iuncoef)
          do jChan = 1, nchans
            if (bias(sensorIndex)%chans(jChan)%coeff_nobs > 0) then
              numPred = bias(sensorIndex)%chans(jChan)%numActivePredictors 
          
              write(iuncoef,'(A38,A8,1X,A7,1X,I6,1X,I2)') 'SATELLITE, INSTRUMENT, CHANNEL, NPRED: ',  &
                   satNameCoeff, instrName, bias(sensorIndex)%chans(jChan)%channelNum, numPred
              do kpred =1, numPred
                write(iuncoef,'(10e14.6)') bias(sensorIndex)%chans(jChan)%coeffCov(kpred, :)
              end do
            end if
          end do
          ierr = fclos(iuncoef) 
        end if

      end do SENSORS

    end if

 end subroutine bcs_writeCoeff

 !-----------------------------------------
 ! bcs_removeBiasCorrection
 !-----------------------------------------
 subroutine bcs_removeBiasCorrection(obsSpaceData, family_opt)
    !
    ! :Purpose: to remove bias correction from OBS_VAR
    !           After the call OBS_VAR contains the uncorrected
    !           observation and OBS_BCOR is set to zero
    !
    implicit none

    ! Arguments:
    type(struct_obs),           intent(inout) :: obsSpaceData
    character(len=2), optional, intent(in)    :: family_opt

    ! Locals:
    integer :: nbcor
    integer :: bodyIndex
    real(8) :: biascor, Obs

    if (.not. removeBiasCorrection) return

    if (mmpi_myid == 0) write(*,*) 'bcs_removeBiasCorrection: start'

    nbcor = 0
    call obs_set_current_body_list(obsSpaceData, family_opt)
    
    BODY: do
      bodyIndex = obs_getBodyIndex(obsSpaceData)
      if (bodyIndex < 0) exit BODY
      if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY  
      biasCor = obs_bodyElem_r(obsSpaceData, OBS_BCOR, bodyIndex)
      Obs = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex)
      if (biasCor /= MPC_missingValue_R8 .and. Obs /= MPC_missingValue_R8) then
        call obs_bodySet_r(obsSpaceData, OBS_VAR, bodyIndex, real(Obs - biasCor, pre_obsReal))
        call obs_bodySet_r(obsSpaceData, OBS_BCOR, bodyIndex, real(0.d0, pre_obsReal))
        nbcor = nbcor + 1
      end if
    end do BODY

    if (mmpi_myid == 0) then
      write(*,*) 'bcs_removeBiasCorrection: removed bias correction for ', nbcor, ' observations'
      write(*,*) 'bcs_removeBiasCorrection exiting'
    end if

  end subroutine bcs_removeBiasCorrection

  !-----------------------------------------
  ! bcs_filterObs
  !-----------------------------------------
  subroutine bcs_filterObs(obsSpaceData)
    !
    ! :Purpose: to filter radiance observations to include into
    !           bias correction offline computation
    !           (same rules as in bgck.gen_table)
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData

    ! Locals:
    integer :: bodyIndex, headerIndex, instrumentIndex, sensorIndex
    integer, allocatable :: instrumentList(:)
    integer :: assim, flag, codtyp, channelNumber
    integer :: isatBufr, instBufr, iplatform, isat, inst, idsat, chanIndx
    logical :: lHyperIr, lGeo, lSsmis, lTovs
    logical :: condition, condition1, condition2, channelIsAllsky
    logical :: channelIsPassive
    character(len=10)  :: instrName
    
    if (.not. filterObs) return

    if (mmpi_myid == 0) write(*,*) 'bcs_filterObs: start'

    allocate(instrumentList(tvs_nsensors))
    instrumentList(:) = -1
    do sensorIndex = 1, tvs_nsensors
      instrName = InstrNameinCoeffFile(tvs_instrumentName(sensorIndex))
      do instrumentIndex = 1,  maxNumInst
        if (trim(instrName) == trim(cinst(instrumentIndex)))    then
          instrumentList(sensorIndex) = instrumentIndex
        end if
      end do
      if (instrumentList(sensorIndex) == -1) then
        call utl_abort('bcs_filterObs: Instrument ' // trim(instrName) // &
            'missing in CINST table fron NAMBIASSAT namelist section')
      end if
    end do
    
    call obs_set_current_header_list(obsSpaceData, 'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)

      lHyperIr = .false.
      lGeo = .false.
      lSsmis = .false.
      lTovs = .false.
 
      select case (codtyp_get_name(codtyp))
      case("ssmis")
        lSsmis = .true.
      case("csr")
        lGeo = .true.
      case("airs","iasi","cris","crisfsr")
        lHyperIr = .true.
      case default
        lTovs = .true.
      end select

      isatBufr = obs_headElem_i(obsSpaceData, OBS_SAT, headerIndex) !BUFR element 1007
      instBufr = obs_headElem_i(obsSpaceData, OBS_INS, headerIndex) !BUFR element 2019

      call tvs_mapSat(isatBufr, iplatform, isat)
      call tvs_mapInstrum(instBufr, inst)

      idsat = -1
      do sensorIndex = 1, tvs_nsensors
        if (tvs_platforms(sensorIndex) == iplatform .and.  &
            tvs_satellites(sensorIndex) == isat     .and.  &
            tvs_instruments(sensorIndex) == inst)   then
          idsat = sensorIndex
          exit
        end if
      end do
      if (idsat == -1) cycle HEADER

      call obs_set_current_body_list(obsSpaceData, headerIndex)

      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY

        ! determine if instrument/channel function in all-sky mode
        channelNumber = nint(obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex))
        channelIsAllsky = .false.
        if ((tvs_isInstrumUsingCLW(tvs_instruments(idsat)) .or. &
             tvs_isInstrumUsingHydrometeors(tvs_instruments(idsat))) .and. &
            oer_useStateDepSigmaObs(channelNumber, idsat)) then
          channelIsAllsky = .true.
        end if
        channelIsPassive = .false.
        if (passiveChannelNumber(instrumentList(idsat)) > 0) then
          channelIsPassive = (utl_findloc(passiveChannelList(instrumentList(idsat),1:passiveChannelNumber(instrumentList(idsat))), channelNumber) > 0)
        end if
        assim = obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex)
        if (assim == obs_notAssimilated) then
          call bcs_getChannelIndex(obsSpaceData, idsat, chanIndx, bodyIndex)
          if (chanIndx > 0) then
            flag = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex)

            if (lSsmis) then   ! SSM/I and SSMIS
              !   Here "good" means data that are not rejected by first QC program satqc_ssmi(s) (bit 7 ON).
              !   Bit 11 is ON for data that are unselected by UTIL or for uncorrected data (or both).
              !   Data rejected by first QC program satqc_ssmi(s) have bit 7 switched ON only (in addition to bit 9) as
              !   rogue/topo checks are skipped. So if bit 16 (rogue) is ON, bit 7 must be off.
              if (.not. offlineMode .and. .not. channelIsPassive) then
                if (allModeSsmis) then
                  !  FLAG test: all good data (corrected/selected or not) that have passed all QC (bit 9 OFF)
                  condition1 = .not. btest(flag, 9) !' AND (FLAG & 512 = 0)'
                  !  FLAG test: uncorrected good data that failed rogue check only ([bit 9 ON] + bit 6 OFF + bit 16 ON + bit 18 OFF + [bit 7 OFF])
                  condition2 = .not. btest(flag, 6) .and. btest(flag, 16) .and. .not. btest(flag, 18) !' AND (FLAG & 64 = 0) AND (FLAG &  65536 = 65536) AND (FLAG & 262144 = 0)'
                  condition = condition1 .or. condition2
                else
                  !  FLAG test: corrected/selected good data that have passed QC (bits 9,11 OFF) --> data to be assimilated
                  condition = .not. btest(flag, 9) .and. .not. btest(flag, 11)       !' AND (FLAG & 512 = 0) AND (FLAG & 2048 = 0)'
                end if
              else
                ! OFFLINE MODE --> want all observations except data rejected for any reason other than rogue innovation check
                condition1 = .not. btest(flag, 9) !' AND (FLAG & 512 = 0)'  
                ! all good data that passed all QC    
                ! "good" data that failed rogue check [bit 9 ON, bit 7 OFF, bit 18 OFF]
                condition2 = btest(flag, 9) .and. .not. btest(flag, 7) .and. .not. btest(flag, 18) !' AND (FLAG & 512 = 512) AND (FLAG & 128 = 0) AND (FLAG & 262144 = 0)'
                condition = condition1 .or. condition2
              end if
            else if(lTovs) then
              ! AMSU-A, AMSU-B/MHS, ATMS, MWHS-2
              !  In AMSU case, bit 11 is set for data that are not bias corrected or for unselected channels.
              !    BUT unlike other instruments, all AMSU data are bias corrected, whether selected or not
              !    so bit 11 = unselected channel (like bit 8 for AIRS/IASI)
              !  Bit 9 is set for all other rejections including rogue (9+16) and topography (9+18).
              !  In addition, bit 7 is set for channels with bad data or data that should not be assimilated.
              if (.not. offlineMode .and. .not. channelIsPassive) then
                if (allModeTovs) then
                  !  FLAG test: all data (selected or not) that have passed QC (bit 9 OFF)
                  condition1 = .not. btest(flag, 9) !' AND (FLAG & 512 = 0)'
                  !  FLAG test: uncorrected (bit 6 OFF) data that failed rogue check only (bit (9)/16 ON, 18,7 OFF)
                  !             NOTE: As all AMSU data are normally bias corrected, query2 will return nothing
                  condition2 = btest(flag, 16) .and. .not. btest(flag, 6) .and. .not. btest(flag, 18) .and. .not. btest(flag, 7)!' AND (FLAG & 64 = 0) AND (FLAG &  65536 = 65536) AND (FLAG & 262144 = 0) AND (FLAG & 128 = 0)'
                  condition = condition1 .or. condition2
                else
                  !  FLAG test: selected data (bit 11 OFF) that have passed QC (bit 9 OFF)
                  condition = .not. btest(flag, 9) .and. .not. btest(flag, 11) !' AND (FLAG & 512 = 0) AND (FLAG & 2048 = 0)'
                end if
              else    ! OFFLINE MODE --> want all observations except data rejected for any reason other than rogue check
                condition1 = .not. btest(flag, 9) !' AND (FLAG & 512 = 0)'  
                ! all good data that passed all QC    
                ! "good" data that failed rogue check [bit 9 ON, bit 7 OFF, bit 18 OFF]
                condition2 =  btest(flag, 9) .and. .not. btest(flag, 7) .and. .not. btest(flag, 18)  !' AND (FLAG & 512 = 512) AND (FLAG & 128 = 0) AND (FLAG & 262144 = 0)'
                condition = condition1 .or. condition2
              end if
              if (channelIsAllsky) condition = condition .and. .not. btest(flag, 23)
            else if(lGeo) then  ! CSR case
              !    No flag check        =                all data that have passed QC/filtering
              !  (FLAG & 2048 = 0)      = bit 11 OFF --> corrected/selected data that have passed QC/filtering
              if (allModeCsr .or. offlineMode .or. channelIsPassive) then
                condition = .true.
              else        
                condition = .not. btest(flag, 18) ! ' AND (FLAG & 2048 = 0)' 
              endif
            else if (lHyperIr) then ! AIRS, IASI and CRIS
              !  (FLAG & 2560 = 0)     = bits 9, 11 OFF       --> data that passed QC (rogue and other)
              !  (FLAG & 11010176 = 0) = bits 7,19,21,23 OFF  --> "good" data (corrected/selected or not)!  (FLAG & 64 = 64)  = bit 6 ON        --> bias corrected data
              !  (FLAG & 256 = 0)  = bit 8 OFF       --> passed selction (not blacklisted, UTIL=1)
              !  (FLAG & 2048 = 2048)   = bit 11 ON
              !  (FLAG & 65536 = 65536) = bit 16 ON  --> rogue check failure
              !  (FLAG & 524288 = 0)    = bit 19 OFF --> not surface affected [experimental, bit 11 may not be on if data
              !                                          are to be assimilated]
              !  (FLAG & 2097152 = 0)   = bit 21 OFF --> not rejected due to model top transmittance
              !  (FLAG & 8388608 = 0)   = bit 23 OFF --> "clear sky" radiance [experimental, bit 11 may not be on if cloudy
              !                                          data are assimilated]!  (FLAG & 128 = 0)      = bit 7 OFF  --> not shortwave channel during day
              !  (FLAG & 512 = 0)      = bit 9 OFF  --> non-erroneous data that passed O-P rogue check
              !! AIRS, IASI:!    bit  8 ON: blacklisted/unselected channel (UTIL=0)
              !    bit  9 ON: erroneous/suspect data (9), data failed O-P check (9+16)
              !    bit 11 ON: cloud (11+23), surface (11+19), model top transmittance (11+21), shortwave channel+daytime (11+7)
              !               not bias corrected (11) (with bit 6 OFF)
              if (.not. offlineMode .and. .not. channelIsPassive) then
                if (allModeHyperIr) then        
                  ! good data that have passed all QC (bits 9 and 7,19,21,23 OFF), corrected/selected or not
                  condition1  = .not. btest(flag, 9) .and. .not. btest(flag, 7) .and. .not. btest(flag, 19) .and. .not. btest(flag, 21) .and. .not. btest(flag, 23) !' AND (FLAG & 512 = 0) AND (FLAG & 11010176 = 0)'
                  ! uncorrected (6 OFF, [11 ON]) good data (7,19,21,23 OFF) that failed QC rogue check only (bits [9],16 ON), selected or not
                  condition2  = .not. btest(flag, 6) .and. btest(flag, 11) .and.  .not. btest(flag, 17) .and. .not. btest(flag, 19) .and. .not. btest(flag, 21) .and. .not. btest(flag, 23) 
                  !' AND (FLAG & 64 = 0) AND (FLAG & 65536 = 65536) AND (FLAG & 11010176 = 0)'
                  condition = condition1 .or. condition2
                else 
                  ! corrected data that passed all QC and selection excluding cloud/sfc affected obs
                  condition =  .not. btest(flag, 9) .and. .not. btest(flag, 11) .and.  .not. btest(flag, 8) .and. .not. btest(flag, 23) .and. .not. btest(flag, 19) 
                  !' AND (FLAG & 2560 = 0) AND (FLAG & 256 = 0) AND (FLAG & 8388608 = 0) AND (FLAG & 524288 = 0)'
                end if
              else! OFFLINE MODE --> Want all observations except data rejected for any reason other than innovation rogue check
                !   Assumes that type S or N correction has been applied to all data/channels (all data "corrected")
                ! data that passed all QC 
                condition1 =  .not. btest(flag, 9) .and. .not. btest(flag, 7) .and. .not. btest(flag, 19) .and. .not. btest(flag, 21) .and. .not. btest(flag, 23) 
                !' AND (FLAG & 512 = 0) AND (FLAG & 11010176 = 0)'
                ! good data (7,19,21,23 OFF) that failed QC rogue check only (bits [9],16 ON)
                condition2 = btest(flag, 9) .and. btest(flag, 16) .and. .not. btest(flag, 7) .and. .not. btest(flag, 19) .and. .not. btest(flag, 21) .and. .not. btest(flag, 23) !' AND (FLAG & 65536 = 65536) AND (FLAG & 11010176 = 0)'
                condition = condition1 .or. condition2
              end if
            end if

            if (condition) assim = obs_Assimilated

            call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, assim)
          end if
        end if
      end do BODY
    end do HEADER

    deallocate(instrumentList)
 
  end subroutine bcs_filterObs

  !-----------------------------------------
  ! bcs_applyBiasCorrection
  !-----------------------------------------
  subroutine bcs_applyBiasCorrection(obsSpaceData, column, family_opt)
    !
    ! :Purpose: to apply bias correction from OBS_BCOR to
    !           obsSpaceData column column.
    !           After the call obsSpaceData body column contains the corrected
    !           observation or O-F and OBS_BCOR is not modified.
    implicit none

    ! Arguments:
    type(struct_obs),           intent(inout) :: obsSpaceData
    integer,                    intent(in)    :: column !obsSpaceData column
    character(len=2), optional, intent(in)    :: family_opt

    ! Locals:
    integer :: nbcor
    integer :: bodyIndex
    real(8) :: biascor, Obs
    integer :: flag

    if (.not. biasActive) return
    if (trim(biasMode) /= "apply") return

    if (mmpi_myid == 0) write(*,*) 'bcs_applyBiasCorrection: start'

    nbcor = 0
    call obs_set_current_body_list(obsSpaceData, family_opt)
    
    BODY: do
      bodyIndex = obs_getBodyIndex(obsSpaceData)
      if (bodyIndex < 0) exit BODY
      if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY  
      biasCor = obs_bodyElem_r(obsSpaceData, OBS_BCOR, bodyIndex)
      if (biasCor /= MPC_missingValue_R8) then
        Obs = obs_bodyElem_r(obsSpaceData, column, bodyIndex)
        if (Obs /= MPC_missingValue_R8) then
          call obs_bodySet_r(obsSpaceData, column, bodyIndex, real(Obs + biasCor, pre_obsReal))
          flag = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex)
          flag = ibset(flag, 6)
          call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, flag)
          nbcor = nbcor + 1
        end if
      end if
    end do BODY

    if (mmpi_myid == 0) then
      write(*,*) 'bcs_applyBiasCorrection: apply bias correction for ', nbcor, ' observations'
      write(*,*) 'bcs_applyBiasCorrection exiting'
    end if

  end subroutine bcs_applyBiasCorrection

  !-----------------------------------------
  ! bcs_refreshBiasCorrection
  !-----------------------------------------
  subroutine bcs_refreshBiasCorrection(obsSpaceData, columnTrlOnTrlLev)
    !
    ! :Purpose: to apply bias correction from read coefficient file to OBS_VAR
    !           After the call OBS_VAR contains the corrected observation
    !           and OBS_BCOR is set to applied bias correction
    !
    implicit none

    ! Arguments:
    type(struct_obs),        intent(inout) :: obsSpaceData
    type(struct_columnData), intent(inout) :: columnTrlOnTrlLev

    if (.not.biasActive) return
    if (.not. refreshBiasCorrection) return

    if (mmpi_myid == 0) write(*,*) 'bcs_refreshBiasCorrection: start'
    call bcs_calcBias(obsSpaceData, columnTrlOnTrlLev)
    call bcs_applyBiasCorrection(obsSpaceData, OBS_VAR, "TO")
    if (mmpi_myid == 0) write(*,*) 'bcs_refreshBiasCorrection: exit'

  end subroutine bcs_refreshBiasCorrection

  !-----------------------------------------
  ! bcs_getRadiosondeWeight
  !-----------------------------------------
  subroutine bcs_getRadiosondeWeight(obsSpaceData, lmodify_obserror_opt)
    !
    ! :Purpose: initialize the weights to give more importance to data near radiosonde stations
    !
    implicit none

    ! Arguments:
    type(struct_obs),  intent(inout) :: obsSpaceData
    logical, optional, intent(in)    :: lmodify_obserror_opt

    ! Locals:
    integer :: iobs, headerIndex, idatyp, nobs, bodyIndex, stepIndex
    logical :: lmodify_obserror
    real(8) :: sigmaObs
    type(struct_gsv)      :: stateVector_mask, stateVector_mask_4d

    lmodify_obserror =.false.

    if (present(lmodify_obserror_opt)) lmodify_obserror = lmodify_obserror_opt


    if (weightedEstimate) then
      call hco_SetupFromFile(hco_mask, './raob_masque.std', 'WEIGHT', GridName_opt='RadiosondeWeight', varName_opt='WT')
      call vco_SetupFromFile(vco_mask, './raob_masque.std')   ! IN

      call gsv_allocate(stateVector_mask, 1, hco_mask, vco_mask, dateStampList_opt=[-1], varNames_opt=["WT"], &
           dataKind_opt=4, mpi_local_opt=.true., mpi_distribution_opt="Tiles")
      call gsv_allocate(stateVector_mask_4d, tim_nstepobs, hco_mask, vco_mask, dateStampList_opt=[ (-1, stepIndex = 1, tim_nstepobs) ], varNames_opt=["WT"], &
           dataKind_opt=4, mpi_local_opt=.true., mpi_distribution_opt="Tiles")

      call gio_readFromFile(stateVector_mask, './raob_masque.std', 'WEIGHT', 'O', unitConversion_opt=.false., &
           containsFullField_opt=.false.)

      do stepIndex = 1, tim_nstepobs
        call gsv_copy(stateVector_mask, stateVector_mask_4d, stepIndexOut_opt=stepIndex)
      end do
      call gsv_deallocate(stateVector_mask)
      call col_setVco(column_mask, vco_mask)
      nobs = obs_numHeader(obsSpaceData)
      call col_allocate(column_mask, nobs, beSilent_opt=.false., varNames_opt=["WT"])

      call s2c_nl(stateVector_mask_4d, obsSpaceData, column_mask, hco_mask, &
                   'NEAREST', varName_opt="WT", moveObsAtPole_opt=.true.)

      call obs_set_current_header_list(obsSpaceData, 'TO')
      iobs = 0
      HEADER: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER
          
        ! process only radiance data to be assimilated?
        idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
        if (.not. tvs_isIdBurpTovs(idatyp)) then
          write(*,*) 'bcs_getRadiosondeWeight: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
          cycle HEADER
        end if
        iobs = iobs + 1

        RadiosondeWeight(iobs) = col_getElem(column_mask, 1, headerIndex) 

        if (lmodify_obserror) then
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY
            
            if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY

            sigmaObs = obs_bodyElem_r(obsSpaceData, OBS_OER, bodyIndex)

            sigmaObs = sigmaObs / sqrt(RadiosondeWeight(iobs))
            call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, sigmaObs) 

          end do BODY

        end if

      end do HEADER
      call gsv_deallocate(stateVector_mask_4d)
    else
      RadiosondeWeight(:) = 1.d0
    end if

  end subroutine bcs_getRadiosondeWeight

  !-----------------------------------------
  ! bcs_do_regression
  !-----------------------------------------
  subroutine bcs_do_regression(columnTrlOnTrlLev, obsSpaceData)
    !
    ! :Purpose: compute the bias correction coefficients by linear regresion
    !
    implicit none

    ! Arguments:
    type(struct_obs),        intent(inout) :: obsSpaceData
    type(struct_columnData), intent(inout) :: columnTrlOnTrlLev

    ! Locals:
    integer    :: iSensor, iChannel, npred, nchans, nscan, ndim, ndimmax
    integer    :: sensorIndex, iPred1, jPred1, iobs
    integer    :: headerIndex, idatyp, nPredMax, ierr, iFov, iScan, idim
    integer    :: indxtovs, bodyIndex, chanIndx, predstart, ntot
    real(8)    :: OmF, sigmaObs, lambda, norm
    real(8)    :: predictor(NumPredictors)
    real(8), allocatable :: Matrix(:,:,:), Vector(:,:)
    real(8), allocatable :: matrixMpiGlobal(:,:,:), vectorMpiGlobal(:,:)
    real(8), allocatable :: pIMatrix(:,:), OmFBias(:,:), omfBiasMpiGlobal(:,:)
    real(8), allocatable :: BMatrixMinusOne(:,:), LineVec(:,:)
    integer, allocatable :: OmFCount(:,:), omfCountMpiGlobal(:,:)
    character(len=80)    :: errorMessage

    write(*,*) "bcs_do_regression: start"
    if (.not. allocated(trialHeight300m1000)) then
      call bcs_getTrialPredictors(obsSpaceData, columnTrlOnTrlLev)
      call bcs_computePredictorBiases(obsSpaceData)
    end if

    call  bcs_getRadiosondeWeight(obsSpaceData)

    if (outOmFPredCov) call bcs_outputCvOmPPred(obsSpaceData)

    SENSORS:do sensorIndex = 1, tvs_nsensors

      if  (.not. tvs_isReallyPresentMpiGLobal(sensorIndex)) cycle SENSORS

      nchans = bias(sensorIndex)%numChannels
      nscan = bias(sensorIndex)%numscan
      npredMax = maxval(bias(sensorIndex)%chans(:)%numActivePredictors)

      if (mimicSatbcor) then
        ndimmax = npredMax
        allocate(OmFBias(nchans,nscan))
        OmFBias(:,:) = 0.d0
      else
        ndimmax = npredMax + nscan - 1
      end if

      allocate(OmFCount(nchans,nscan))
      OmFCount(:,:) = 0

      allocate(Matrix(nchans,ndimmax,ndimmax))
      Matrix(:,:,:) = 0.d0

      allocate(Vector(nchans,ndimmax))
      Vector(:,:) = 0.d0

      allocate(LineVec(1,ndimmax))

      allocate(pIMatrix(ndimmax,ndimmax))


      ! First pass throught ObsSpaceData to estimate scan biases and count data

      call obs_set_current_header_list(obsSpaceData, 'TO')
      HEADER1: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER1
          
        ! process only radiance data to be assimilated?
        idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
        if (.not. tvs_isIdBurpTovs(idatyp)) then
          write(*,*) 'bcs_do_regression: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
          cycle HEADER1
        end if
        indxtovs = tvs_tovsIndex(headerIndex)
        if (indxtovs < 0) cycle HEADER1

        iSensor = tvs_lsensor(indxTovs)
        if (iSensor /= sensorIndex) cycle HEADER1
        iFov = obs_headElem_i(obsSpaceData, OBS_FOV, headerIndex)
        if (nscan > 1) then
          iScan = iFov
        else
          iScan = 1
        end if

        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY1: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY1
            
          if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY1 
          call bcs_getChannelIndex(obsSpaceData, iSensor, chanIndx, bodyIndex)
          if (chanindx > 0) then
            if  (mimicSatBcor) then
              OmF = obs_bodyElem_r(obsSpaceData, OBS_OMP, bodyIndex)
              OmFBias(chanIndx,iScan) = OmFBias(chanIndx,iScan) + OmF
            end if
            OmFCount(chanIndx,iScan) =  OmFCount(chanIndx,iScan) + 1
          end if
        end do BODY1
      end do HEADER1

      if (mimicSatbcor) allocate(omfBiasMpiGlobal(nchans,nscan))

      allocate(omfCountMpiGlobal(nchans,nscan))

      if (mimicSatbcor) then
        call mmpi_reduce_sumR8_2d( OmFBias, omfBiasMpiGlobal, 0, "GRID" )
      end if
      call rpn_comm_reduce(OmFCount, omfCountMpiGlobal, size(omfCountMpiGlobal), "mpi_integer", "MPI_SUM", 0, "GRID", ierr)

      if (ierr /= 0) then
        write(errorMessage,*) "bcs_do_regression: MPI communication error 2", ierr 
        call utl_abort(errorMessage)
      end if
      if (mimicSatbcor)  then
        if (mmpi_myId == 0) then
          where(omfCountMpiGlobal == 0) omfBiasMpiGlobal = 0.d0
          where(omfCountMpiGlobal > 0) omfBiasMpiGlobal = omfBiasMpiGlobal / omfCountMpiGlobal
        end if
        call rpn_comm_bcast(omfBiasMpiGlobal, size(omfBiasMpiGlobal), "mpi_double_precision", 0, "GRID", ierr)
        if (ierr /= 0) then
          write(errorMessage,*) "bcs_do_regression: MPI communication error 3", ierr 
          call utl_abort(errorMessage)
        end if
        do iChannel = 1, nchans
          bias(sensorIndex)%chans(iChannel)%coeff_fov(:) = omfBiasMpiGlobal(iChannel, :)
        end do
        deallocate(omfBiasMpiGlobal)
      end if
     
      ! Second pass to fill matrices and vectors
      call obs_set_current_header_list(obsSpaceData, 'TO')
      iobs = 0
      HEADER2: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER2

        ! process only radiance data to be assimilated?
        idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
        if (.not. tvs_isIdBurpTovs(idatyp)) then
          write(*,*) 'bcs_do_regression: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
          cycle HEADER2
        end if
        iobs = iobs + 1
        indxtovs = tvs_tovsIndex(headerIndex)
        if (indxtovs < 0) cycle HEADER2

        iSensor = tvs_lsensor(indxTovs)
        if (iSensor /= sensorIndex) cycle HEADER2
          
        call obs_set_current_body_list(obsSpaceData, headerIndex)
        iFov = obs_headElem_i(obsSpaceData, OBS_FOV, headerIndex)
        if (nscan > 1) then
          iScan = iFov
        else
          iScan = 1
        end if
        BODY2: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY2
          if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY2 
          call bcs_getChannelIndex(obsSpaceData, iSensor, chanIndx, bodyIndex)
          if (chanIndx > 0) then
            call bcs_getPredictors(predictor, headerIndex, iobs, chanIndx, obsSpaceData)
            OmF = obs_bodyElem_r(obsSpaceData, OBS_OMP, bodyIndex)

            if (mimicSatbcor) OmF = OmF - bias(sensorIndex)%chans(chanIndx)%coeff_fov(iScan)
            
            LineVec(:,:) = 0.d0

            if (mimicSatbcor) then
              idim = 0
              predstart = 1
              lambda = 1.d0
            else
              LineVec(1,iScan) = 1.d0
              idim = nscan
              predstart = 2
              sigmaObs = obs_bodyElem_r(obsSpaceData, OBS_OER, bodyIndex)
              lambda = 1.d0 / (sigmaObs ** 2)
            end if
          
            lambda = lambda * RadiosondeWeight(iobs)

            do iPred1 = predstart, bias(iSensor)%chans(chanIndx)%NumActivePredictors
              jPred1 = bias(iSensor)%chans(chanIndx)%PredictorIndex(iPred1)
              idim = idim + 1
              LineVec(1,idim) = predictor(jPred1)
            end do
            Matrix(chanindx,:,:) = Matrix(chanindx,:,:) + matmul(transpose(LineVec),LineVec) * lambda
            Vector(chanIndx,:) =  Vector(chanIndx,:) + LineVec(1,:) * OmF  * lambda
          end if
        end do BODY2
      end do HEADER2

      if (mmpi_myId == 0) then
        allocate(matrixMpiGlobal(nchans,ndimmax,ndimmax))
        allocate(vectorMpiGlobal(nchans,ndimmax))
      else
        allocate(matrixMpiGlobal(1,1,1))
        allocate(vectorMpiGlobal(1,1))
      end if

      ! communication MPI pour tout avoir sur tache 0
      call mmpi_reduce_sumR8_3d( Matrix, matrixMpiGlobal, 0, "GRID" )
      call mmpi_reduce_sumR8_2d( Vector, vectorMpiGlobal, 0, "GRID" )

      do iChannel = 1, nchans

        if (mmpi_myId == 0) then
          ntot = sum(omfCountMpiGlobal(iChannel, :))
          bias(sensorIndex)%chans(iChannel)%coeff_nobs = ntot
          if (ntot > 0 .and. .not. mimicSatbcor) then
            norm = 1.d0 / (ntot) 
            matrixMpiGlobal(iChannel,:,:) =  matrixMpiGlobal(iChannel,:,:) * norm
            vectorMpiGlobal(iChannel,:) = vectorMpiGlobal(iChannel,:) * norm
          end if

          nPred =  bias(sensorIndex)%chans(iChannel)%numActivePredictors
          if (mimicSatbcor) then
            ndim = npred
          else
            ndim = npred + nscan -1 
            allocate(BMatrixMinusOne(ndim,ndim))
            BMatrixMinusOne(:,:) = 0.d0 
            BMatrixMinusOne(1:nscan,1:nscan) = matmul(bias(sensorIndex)%BMinusHalfScanBias, bias(sensorIndex)%BMinusHalfScanBias)
            do iPred1 = 2, bias(sensorIndex)%chans(iChannel)%numActivePredictors
              BMatrixMinusOne(nscan - 1 + iPred1,nscan - 1 + iPred1) = (1.d0 / (bias(sensorIndex)%chans(iChannel)%stddev(iPred1)) ** 2)
            end do
            matrixMpiGlobal(iChannel,1:ndim,1:ndim) = matrixMpiGlobal(iChannel,1:ndim,1:ndim) + BMatrixMinusOne(:,:)
            deallocate(BMatrixMinusOne)
          end if

          pIMatrix(:,:) = 0.d0
          call utl_pseudo_inverse(matrixMpiGlobal(iChannel, 1:ndim, 1:ndim), pIMatrix(1:ndim, 1:ndim))
          LineVec(1,1:ndim) = matmul(pIMatrix(1:ndim,1:ndim), vectorMpiGlobal(iChannel,1:ndim))
          !call dsymv("L", ndim, 1.d0, pIMatrix, ndim,vectorMpiGlobal(iChannel,:), 1, 0.d0, LineVec(1,1:ndim), 1)
        end if

        call rpn_comm_bcast(ndim, 1, "mpi_integer", 0, "GRID", ierr)
        call rpn_comm_bcast(npred, 1, "mpi_integer", 0, "GRID", ierr)

        call rpn_comm_bcast(LineVec(1,1:ndim), ndim, "mpi_double_precision", 0, "GRID", ierr)
        if (ierr /= 0) then
          write(errorMessage,*) "bcs_do_regression: MPI communication error 6", ierr 
          call utl_abort(errorMessage)
        end if

        if (outCoeffCov) then
          allocate (bias(sensorIndex)%chans(iChannel)%coeffCov(ndim,ndim)) 
          call rpn_comm_bcast(pIMatrix(1:ndim, 1:ndim), ndim * ndim, "mpi_double_precision", 0, "GRID", ierr)
          bias(sensorIndex)%chans(iChannel)%coeffCov(:,:) = pIMatrix(1:ndim,1:ndim)
        end if

        if (mimicSatbcor) then
          bias(sensorIndex)%chans(iChannel)%coeff(:) = LineVec(1,1:npred)
        else
          bias(sensorIndex)%chans(iChannel)%coeff_fov(:) = LineVec(1,1:nscan)
          bias(sensorIndex)%chans(iChannel)%coeff(1) = 0.d0
          bias(sensorIndex)%chans(iChannel)%coeff(2:) = LineVec(1,nscan+1:ndim)
        end if

      end do

      deallocate(LineVec)
      deallocate(Matrix)
      deallocate(Vector  )
      deallocate(omfCountMpiGlobal)
      deallocate(matrixMpiGlobal)
      deallocate(vectorMpiGlobal)
      deallocate(pIMatrix)
      if (allocated(OmFBias)) deallocate(OmFBias)
      deallocate(OmFCount)
       
    end do SENSORS

  end subroutine bcs_do_regression

  !-----------------------------------------
  ! bcs_outputCvOmPPred
  !-----------------------------------------
  subroutine bcs_outputCvOmPPred(obsSpaceData)
    !
    ! :Purpose: compute and output OmF-predictors covariance and correlation matrices
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout)        :: obsSpaceData

    ! Locals:
    integer :: sensorIndex, headerIndex, bodyIndex, channelIndex, predictorIndex,predictorIndex2,nchans
    integer :: idatyp, indxtovs, iSensor, chanIndx, Iobs, ierr
    Real(8):: OmF
    real(8), allocatable :: OmFBias(:), Matrix(:,:,:), PredBias(:,:)
    integer, allocatable :: Count(:), CountMpiGlobal(:)
    real(8), allocatable :: OmFBiasMpiGlobal(:), predBiasMpiGlobal(:,:), MatrixMpiGLobal(:,:,:)
    character(len=128)   :: errorMessage
    real(8) :: vector(1,numPredictors), predictor(numPredictors),correlation(numPredictors,numPredictors)
    real(8) :: sigma(numPredictors)
    integer :: iuncov, iuncorr

    write(*,*) "bcs_outputCvOmPPred: Starting"

    SENSORS:do sensorIndex = 1, tvs_nsensors
  
      if  (.not. tvs_isReallyPresentMpiGLobal(sensorIndex)) cycle SENSORS
      write(*,*) "sensor ", sensorIndex

      nchans = bias(sensorIndex)%numChannels ! from bcif

      allocate(OmFBias(nchans))
      OmFBias(:) = 0.d0

      allocate(predBias(nchans,numPredictors))
      predBias(:,:) = 0.d0

      allocate(Count(nchans))
      Count(:) = 0

      iobs = 0
      ! First pass throught ObsSpaceData to estimate biases and count data
      call obs_set_current_header_list(obsSpaceData, 'TO')
      HEADER1: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER1
        ! process only radiance data to be assimilated?
        idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
        if (.not. tvs_isIdBurpTovs(idatyp)) then
          write(*,*) 'bcs_outputCvOmPPred: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
          cycle HEADER1
        end if
        iobs = iobs + 1
        indxtovs = tvs_tovsIndex(headerIndex)
        if (indxtovs < 0) cycle HEADER1
        iSensor = tvs_lsensor(indxTovs)
        if (iSensor /= sensorIndex) cycle HEADER1
        
        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY1: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY1
          if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY1 
          call bcs_getChannelIndex(obsSpaceData, iSensor, chanIndx, bodyIndex)
          if (chanindx > 0) then
            OmF = obs_bodyElem_r(obsSpaceData, OBS_OMP, bodyIndex)
            OmFBias(chanIndx) = OmFBias(chanIndx) + OmF
            count(chanIndx) = count(chanIndx) + 1
            call bcs_getPredictors(predictor, headerIndex, iobs, chanIndx, obsSpaceData)
            predBias(chanIndx,:) = predBias(chanIndx,:) + predictor(:)
          end if
        end do BODY1
      end do HEADER1

      allocate(omfBiasMpiGlobal(nchans))
      allocate(countMpiGlobal(nchans))
      allocate(predBiasMpiGlobal(nchans,numPredictors))

      call mmpi_reduce_sumR8_1d(OmFBias, omfBiasMpiGlobal, 0, "GRID" )
      call mmpi_reduce_sumR8_2d(predBias, predBiasMpiGlobal, 0, "GRID" )

      call rpn_comm_reduce(count, countMpiGlobal, size(countMpiGlobal), "mpi_integer", "MPI_SUM", 0, "GRID", ierr)
      if (ierr /= 0) then
        write(errorMessage,*) "bcs_outputCvOmPPred: MPI communication error 1", ierr 
        call utl_abort(errorMessage)
      end if

      if (mmpi_myId == 0) then
        where(countMpiGlobal == 0) omfBiasMpiGlobal = 0.d0
        where(countMpiGlobal > 0) omfBiasMpiGlobal = omfBiasMpiGlobal / countMpiGlobal
        do channelIndex =1, nchans
          if (countMpiGlobal(channelIndex) == 0) predBiasMpiGlobal(channelIndex,:) = 0.d0
          if (countMpiGlobal(channelIndex) > 0) predBiasMpiGlobal(channelIndex,:) = predBiasMpiGlobal(channelIndex,:) / countMpiGlobal(channelIndex)
        end do
      end if
      call rpn_comm_bcast(omfBiasMpiGlobal, size(omfBiasMpiGlobal), "mpi_double_precision", 0, "GRID", ierr)
      if (ierr /= 0) then
        write(errorMessage,*) "bcs_outputCvOmPPred: MPI communication error 4", ierr 
        call utl_abort(errorMessage)
      end if
      call rpn_comm_bcast(predBiasMpiGlobal, size(predBiasMpiGlobal), "mpi_double_precision", 0, "GRID", ierr)
      if (ierr /= 0) then
        write(errorMessage,*) "bcs_outputCvOmPPred: MPI communication error 3", ierr 
        call utl_abort(errorMessage)
      end if
      call rpn_comm_bcast(countMpiGlobal, size(countMpiGlobal), "mpi_integer", 0, "GRID", ierr)
      if (ierr /= 0) then
        write(errorMessage,*) "bcs_outputCvOmPPred: MPI communication error 4", ierr 
        call utl_abort(errorMessage)
      end if

      deallocate(OmFBias)
      deallocate(predBias)
      deallocate(Count)
      allocate(matrix(nchans,numPredictors,numPredictors))
      matrix(:,:,:) = 0.d0

      ! Second pass to fill covariance Matrix
      call obs_set_current_header_list(obsSpaceData, 'TO')
      iobs = 0
      HEADER2: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER2

        ! process only radiance data to be assimilated?
        idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
        if (.not. tvs_isIdBurpTovs(idatyp)) then
          write(*,*) 'bcs_outputCvOmPPred: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
          cycle HEADER2
        end if
        iobs = iobs + 1
        indxtovs = tvs_tovsIndex(headerIndex)
        if (indxtovs < 0) cycle HEADER2

        iSensor = tvs_lsensor(indxTovs)
        if (iSensor /= sensorIndex) cycle HEADER2
          
        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY2: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY2
          if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY2 
          call bcs_getChannelIndex(obsSpaceData, iSensor, chanIndx, bodyIndex)
          if (chanIndx > 0) then
            call bcs_getPredictors(predictor, headerIndex, iobs, chanIndx, obsSpaceData)
            predictor(:) = predictor(:) - predBiasMpiGLobal(chanIndx,:)
            OmF = obs_bodyElem_r(obsSpaceData, OBS_OMP, bodyIndex)
            OmF = OmF - omfBiasMpiGlobal(chanIndx)
            vector(1,1) = OmF
            vector(1,2:numPredictors) = predictor(2:numPredictors)
            Matrix(chanindx,:,:) = Matrix(chanindx,:,:) + matmul(transpose(vector),vector)
          end if
        end do BODY2
      end do HEADER2

      if (mmpi_myId == 0) then
        allocate(matrixMpiGlobal(nchans,numPredictors,numPredictors))
      else
        allocate(matrixMpiGlobal(1,1,1))
      end if

      ! communication MPI pour tout avoir sur tache 0
      call mmpi_reduce_sumR8_3d(matrix, matrixMpiGlobal, 0, "GRID" )
      deallocate(matrix)
      deallocate(OmFBiasMpiGLobal)
      deallocate(predBiasMpiGLobal)
      if (mmpi_myId == 0) then
        call utl_open_asciifile("covarianceMatrix_" // trim(tvs_instrumentName(sensorIndex)) //  &
             "_" // trim(tvs_satelliteName(sensorIndex)) // ".dat", iuncov)
        call utl_open_asciifile("correlationMatrix_" // trim(tvs_instrumentName(sensorIndex)) //  &
             "_" // trim(tvs_satelliteName(sensorIndex)) // ".dat", iuncorr)
        do channelIndex = 1, nchans
          if (countMpiGlobal(channelIndex) > 1) then
            matrixMpiGlobal(channelIndex,:,:) = matrixMpiGlobal(channelIndex,:,:) / countMpiGlobal(channelIndex)
           
            write(iuncov,*) "OmF Pred covariance Matrix for channel ", &
                 bias(sensorIndex)%chans(channelIndex)%channelNum, "instrument ", &
                 trim(tvs_instrumentName(sensorIndex))," ", &
                 trim(tvs_satelliteName(sensorIndex))
            write(iuncov,'(10x,A6)',advance="no") "OmF"
            do predictorIndex = 2, numPredictors
              write(iuncov,'(T6,A6,1x)',advance="no") predTab(predictorIndex) 
            end do
            write(iuncov,*)
            write(iuncov,'(A6)',advance="no") "Omf"
            write(iuncov,'(100f12.6)') matrixMpiGlobal(channelIndex,1,:)
            do predictorIndex = 2, numPredictors
              write(iuncov,'(A6)',advance="no") predTab(predictorIndex)
              write(iuncov,'(100f12.6)') matrixMpiGlobal(channelIndex,predictorIndex,:)
            end do
            sigma(:) = 0
            do predictorIndex = 1,numPredictors
              if ( matrixMpiGlobal(channelIndex,predictorIndex,predictorIndex) >0.d0) &
                   sigma(predictorIndex) = sqrt(matrixMpiGlobal(channelIndex,predictorIndex,predictorIndex))
            end do
            do predictorIndex = 1, numPredictors
              do predictorIndex2 =1, numPredictors
                correlation(predictorIndex, predictorIndex2) =  &
                     matrixMpiGlobal(channelIndex,predictorIndex,predictorIndex2) / &
                     (sigma(predictorIndex) * sigma(predictorIndex2) )
              end do
            end do
            write(iuncorr,*) "OmF Pred correlation Matrix for channel ", &
                 bias(sensorIndex)%chans(channelIndex)%channelNum,"instrument ", &
                 trim(tvs_instrumentName(sensorIndex))," ", &
                 trim(tvs_satelliteName(sensorIndex))
            write(iuncorr,'(10x,A6)',advance="no") "OmF"
            do predictorIndex = 2, numPredictors
              write(iuncorr,'(T6,A6,1x)',advance="no") predTab(predictorIndex) 
            end do
            write(iuncorr,*)
            write(iuncorr,'(A6)',advance="no") "Omf"
            write(iuncorr,'(100f12.6)') correlation(1,:)
            do predictorIndex = 2, numPredictors
              write(iuncorr,'(A6)',advance="no") predTab(predictorIndex)
              write(iuncorr,'(100f12.6)') correlation(predictorIndex,:)
            end do
          end if
        end do
        ierr = fclos(iuncov)
        ierr = fclos(iuncorr)
      end if
          
      deallocate(countMpiGlobal)
      deallocate(matrixMpiGlobal)
       
    end do SENSORS

  end subroutine bcs_outputCvOmPPred

  !----------------------
  ! bcs_Finalize
  !----------------------
  subroutine bcs_Finalize
    !
    ! :Purpose: release allocated memory for the module
    !
    implicit none

    ! Locals:
    integer    :: iSensor, iChannel

    if (.not. biasActive) return

    if (allocated(trialHeight300m1000)) deallocate(trialHeight300m1000)
    if (allocated(trialHeight50m200)) deallocate(trialHeight50m200)
    if (allocated(trialHeight1m10)) deallocate(trialHeight1m10)
    if (allocated(trialHeight5m50)) deallocate(trialHeight5m50)
    if (allocated(trialTG)) deallocate(trialTG)
    if (allocated(RadiosondeWeight)) deallocate(RadiosondeWeight)

    do iSensor = 1, tvs_nSensors
      if (allocated(bias(iSensor)%BHalfScanBias)) &
           deallocate(bias(iSensor)%BHalfScanBias)
      if (allocated(bias(iSensor)%BMinusHalfScanBias)) &
           deallocate(bias(iSensor)%BMinusHalfScanBias)
      do iChannel =1, bias(iSensor)%numChannels
        deallocate(bias(iSensor)%chans(iChannel)%stddev)
        deallocate(bias(iSensor)%chans(iChannel)%coeffIncr)
        deallocate(bias(iSensor)%chans(iChannel)%predictorIndex)
        if (allocated(bias(iSensor)%chans(iChannel)%coeffIncr_fov)) deallocate(bias(iSensor)%chans(iChannel)%coeffIncr_fov)
        deallocate(bias(iSensor)%chans(iChannel)%coeff_offset)
        if (allocated(bias(iSensor)%chans(iChannel)%coeff)) deallocate(bias(iSensor)%chans(iChannel)%coeff)
        if (allocated(bias(iSensor)%chans(iChannel)%coeff_fov))  deallocate(bias(iSensor)%chans(iChannel)%coeff_fov)
        if (allocated(bias(iSensor)%chans(iChannel)%coeffCov)) deallocate(bias(iSensor)%chans(iChannel)%coeffCov)
      end do
      deallocate(bias(iSensor)%chans)
    end do

  end subroutine bcs_Finalize 
 
  !-----------------------------
  ! InstrNametoCoeffFileName 
  !-----------------------------
  function InstrNametoCoeffFileName(nameIn) result(nameOut)
    implicit none

    ! Arguments:
    character(len=10), intent(in) :: nameIn
    ! Result:
    character(len=10)             :: nameOut

    ! Locals:
    character(len=10)  :: temp_instrName
    integer            :: ierr

    temp_instrName = nameIn
    ierr = clib_tolower(temp_instrName) 
    if (trim(temp_instrName) == 'mhs') then
      nameOut = 'amsub'
    else if (trim(temp_instrName) == 'goesimager') then
      nameOut = 'cgoes'
    else if (trim(temp_instrName) == 'gmsmtsat') then
      nameOut = 'mtsat'
    else if (trim(temp_instrName) == 'mviri') then
      nameOut = 'mets7'
    else
      nameOut = temp_instrName
    end if

  end function InstrNametoCoeffFileName 

  !-----------------------------
  ! InstrNameinCoeffFile
  !-----------------------------
  function InstrNameinCoeffFile(nameIn) result(nameOut)
    implicit none
    
    ! Arguments:
    character(len=10), intent(in) :: nameIn
    ! Result:
    character(len=10)             :: nameOut

    if (trim(nameIn) == 'MHS') then
      nameOut = 'AMSUB'
    else if (trim(nameIn) == 'GOESIMAGER') then
      nameOut = 'CGOES' 
    else if (trim(nameIn) == 'GMSMTSAT') then
      nameOut = 'MTSAT' 
    else if (trim(nameIn) == 'MVIRI') then
      nameOut = 'METS7' 
    else 
      nameOut = nameIn
    end if

  end function InstrNameinCoeffFile

  !-----------------------------
  ! SatNameinCoeffFile
  !-----------------------------
  function SatNameinCoeffFile(nameIn) result(nameOut)
    implicit none

    ! Arguments:
    character(len=10), intent(in) :: nameIn
    ! Result:
    character(len=10)             :: nameOut

    if (trim(nameIn) == 'MSG1') then
      nameOut = 'METSAT8'
    else if (trim(nameIn) == 'MSG2') then
      nameOut = 'METSAT9'
    else if (trim(nameIn) == 'MSG3') then
      nameOut = 'METSAT10'
    else if (trim(nameIn) == 'MSG4') then
      nameOut = 'METSAT11'
    else if (trim(nameIn) == 'METEOSAT7') then
      nameOut = 'METSAT7' 
    else 
      nameOut = nameIn
    end if

  end function SatNameinCoeffFile
  
  !-----------------------------------------
  ! read_bcif
  !-----------------------------------------
  subroutine read_bcif(bcifFile, hspec, ncan, can, bcmode, bctype, npred, pred, global_opt, exitcode)
    !
    ! :Purpose: to read channel-specific bias correction (BC) information (predictors) for instrument from BCIF.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in)  :: bcifFile
    logical,          intent(in)  :: hspec
    integer,          intent(out) :: exitcode
    integer,          intent(out) :: ncan
    integer,          intent(out) :: can(tvs_maxchannelnumber)
    integer,          intent(out) :: npred(tvs_maxchannelnumber)
    character(len=3), intent(in)  :: global_opt
    character(len=1), intent(out) :: bcmode(tvs_maxchannelnumber)
    character(len=1), intent(out) :: bctype(tvs_maxchannelnumber)
    character(len=2), intent(out) :: pred(tvs_maxchannelnumber,numpredictorsBCIF)

    ! Locals:
    character(len=7)             :: instrum
    integer                      :: i, j, ier, ii, iun
    character(len=64)            :: line
    integer                      :: xcan, xnpred, chknp
    character(len=1)             :: xbcmode, xbctype
    character(len=2)             :: xpred(numpredictors)

    ! Reads channel-specific bias correction (BC) information (predictors) for instrument from BCIF.
    ! Channel 0 values are global or default values (optionally applied to all channels).
    ! Returns BC information for all channels to calling routine.
    !
    !              Sample BCIF
    ! AMSUA  15                                                   <--- instrument and number of channels
    !CHAN  MODE  TYPE  NPRED PRED1 PRED2 PRED3 PRED4 PRED5 PRED6
    !   0     D     C      4    T1    T2    T3    T4    XX    XX  <--- channel "0" global or default values
    !   1     D     C      2    T1    T2    XX    XX    XX    XX
    !   2     D     C      3    BT    T1    T2    XX    XX    XX
    !   3     S     C      2    T3    T4    XX    XX    XX    XX
    ! ....
    ! ....
    ! ===================  24 APRIL 2014    LIST-DIRECTED I/O VERSION ==============================================
    !   CALL read_bcif(iunbc, bc_instrum, bc_ncan, bc_can, bc_mode, bc_type, bc_npred, bc_pred, global_opt, exitcode)
    !
    !  global_opt = NON    Read channel-specific data for ALL ncan channels from BCIF (channel 0 ignored)
    !               OUI    Read channel 0 data and apply to all ncan channels (global values)
    !               DEF    Read channel 0 data and apply as default values for all ncan channels;
    !                      then scan the rest of the BCIF for any channel-specific data that will
    !                      override the default values.
    !
    !  NOTE: For hyperspectral instruments (e.g. AIRS, IASI, CrIS) the BCIF must always contain records for ALL ncan channels.
    !          If global_opt = OUI, only the channel numbers are needed in column 1 (CHAN) to get the list
    !            of channel numbers.
    !          If global_opt = DEF, other column data (MODE, TYPE, NPRED, PRED1,...) are entered only for
    !            those channels for which the default (channel 0) values are to be overridden.
    !
    !        For standard instruments (AMSU, SSM/I, SSMIS), with consecutive channels 1,2,3,...,ncan:
    !           If global_opt = OUI, only the channel 0 record is needed in the BCIF (any other records after
    !             channel 0 are ignored (not read).
    !           If global_opt = DEF, the channel 0 record (default values) and only records for those channels
    !             for which values are different from defaults are needed in the BCIF.
    
    exitcode = -1

    iun = 0
    ier = fnom(iun, bcifFile, 'FMT', 0)
    if (ier /= 0) then
      call utl_abort('read_bcif: ERROR - Problem opening the bcif file!' // trim(bcifFile)) 
    end if
    
    read(iun,*) instrum, ncan
    read(iun,'(A64)') line

    ! For GLOBAL option, read global values from first line (channel 0) and clone to all channels
    if (global_opt == 'OUI' .or. global_opt == 'DEF') then 
      ! Read channel 0 information
      read(iun, *, iostat=ier) can(1), bcmode(1), bctype(1), npred(1), (pred(1,j), j = 1, numpredictorsbcif)
      if (ier /= 0) then
        write(*,*) 'read_BCIF: Error reading channel 0 data!'
        exitcode = ier
        return
      end if
      if (can(1) /= 0) then
        write(*,*) 'read_BCIF: Channel 0 global values not found!'
        exitcode = -1
        return
      end if
      ! Clone channel 0 information to all ncan channels
      if (.not. hspec) then
        ! For instruments with consecutive channels 1,2,3,...ncan (e.g. AMSU, SSM/I)
        !  -- no need to read the channel numbers from the BCIF
        do i = 2, ncan+1
          can(i) = i - 1
          bcmode(i) = bcmode(1)
          bctype(i) = bctype(1)
          npred(i) = npred(1)
          do j = 1, numpredictorsbcif
            pred(i,j) = pred(1,j)
          end do
        end do
      else
        ! For hyperspectral instruments (channel subsets), read the channel numbers from the BCIF
        do i = 2, ncan + 1
          read(iun, *, iostat=ier) xcan
          if (ier /= 0) then
            write(*,*) 'read_BCIF: Error reading channel numbers!'
            exitcode = ier
            return
          end if
          can(i) = xcan
          bcmode(i) = bcmode(1)
          bctype(i) = bctype(1)
          npred(i) = npred(1)
          do j = 1, numpredictorsbcif
            pred(i,j) = pred(1,j)
          end do
        end do
        ! Reposition the file to just after channel 0 record
        rewind (iun)
        read(iun,*) instrum, ncan
        read(iun,'(A64)') line
        read(iun,*,iostat=ier) xcan, xbcmode, xbctype, xnpred, (xpred(j), j = 1, numpredictorsbcif)
      end if
      ! For global_opt == 'DEF' check for channel-specific information and overwrite the default (channel 0) values
      ! for the channel with the values from the file
      if (global_opt == 'DEF') then
        if (.not. hspec) then
          do
            read(iun,*,iostat=ier) xcan, xbcmode, xbctype, xnpred, (xpred(j), j = 1, numpredictorsbcif)
            if (ier < 0) exit  
            if (ier > 0) then
              write(*,*) 'read_BCIF: Error reading file!'
              exitcode = ier
              return
            end if
            ii = xcan + 1
            if (ii > ncan + 1) then
              write(*,*) 'read_BCIF: Channel number in BCIF exceeds number of channels!'
              write(*,'(A,1X,I4,1X,I4)') '           Channel, ncan = ', xcan, ncan
              exitcode = -1
              return
            end if
            bcmode(ii) = xbcmode
            bctype(ii) = xbctype
            npred(ii) = xnpred
            do j = 1, numpredictorsbcif
              pred(ii,j) = xpred(j)
            end do
          end do
        else
          ! For hyperspectral instruments
          do i = 2, ncan + 1
            read(iun, *,iostat=ier) xcan, xbcmode, xbctype, xnpred, (xpred(j), j = 1, numpredictorsbcif)
            if (ier /= 0) cycle  
            bcmode(i) = xbcmode
            bctype(i) = xbctype
            npred(i) = xnpred
            do j = 1, numpredictorsbcif
              pred(i,j) = xpred(j)
            end do
          end do
        end if
      end if
    end if
    ! Non-GLOBAL: Read the entire file for channel specific values (all channels)
    ! ---------------------------------------------------------------------------------------------------------------------
    if (global_opt == 'NON') then
      ii = 1
      do
        read(iun,*,iostat=ier) can(ii), bcmode(ii), bctype(ii), npred(ii), (pred(ii,j), j = 1, numpredictorsbcif)
        if (ier < 0) exit  
        if (ier > 0) then
          write(*,*) 'read_BCIF: Error reading file!'
          exitcode = ier
          return
        end if
        if (ii == 1) then
          if (can(ii) /= 0) then
            write(*,*) 'read_BCIF: Channel 0 global/default values not found!'
            exitcode = -1
            return
          end if
        end if
        ii = ii + 1
      end do
      if (ii - 2 < ncan) then
        write(*,*) 'read_BCIF: Number of channels in file is less than specified value (NCAN). Changing value of NCAN.'
        ncan = ii - 2
      end if
      if (ii > ncan + 2) then
        write(*,*) 'read_BCIF: ERROR -- Number of channels in file is greater than specified value (NCAN)!'
        exitcode = -1
        return
      end if
    end if
    ! ---------------------------------------------------------------------------------------------------------------------
    write(*,*) ' '
    write(*,*) 'read_BCIF: Bias correction information for each channel (from BCIF):'
    write(*,'(1X,A7,1X,I4)') instrum, ncan
    write(*,*) line
    do i = 1, ncan + 1
      chknp = count(pred(i,:) /= 'XX')
      if (chknp /= npred(i)) npred(i) = chknp
      if (npred(i) == 0 .and. bctype(i) == 'C') bctype(i) = 'F'
      write(*,'(I4,2(5X,A1),5X,I2,6(4X,A2))') can(i), bcmode(i), bctype(i), npred(i), (pred(i,j), j = 1, numpredictorsbcif)
    end do
    write(*,*) 'read_BCIF: exit '

    ier = fclos(iun)
    exitcode = 0

  end subroutine read_bcif
  
  !-----------------------------------------
  ! read_coeff
  !-----------------------------------------
  subroutine read_coeff(sats, chans, fovbias, coeff, nsat, nchan, nfov, npred, cinstrum, coeff_file, ptypes, ndata)
    !
    ! :Purpose: to read radiance bias correction coefficients file
    !
    implicit none

    ! Arguments:
    character(len=10), intent(out) :: sats(:)       ! dim(maxsat), satellite names 1
    integer,           intent(out) :: chans(:,:)    ! dim(maxsat, maxchan), channel numbers 2
    real(8),           intent(out) :: fovbias(:,:,:)! dim(maxsat,maxchan,maxfov), bias as F(fov) 3
    real(8),           intent(out) :: coeff(:,:,:)  ! dim(maxsat,maxchan,maxpred+1) 4
    integer,           intent(out) :: nsat          !5
    integer,           intent(out) :: nchan(:)      ! dim(maxsat), number of channels 6
    integer,           intent(out) :: nfov          !7
    integer,           intent(out) :: npred(:,:)    ! dim(maxsat, maxchan), number of predictors !8
    character(len=7),  intent(out) :: cinstrum      ! string: instrument (e.g. AMSUB) 9
    character(len=*),  intent(in)  :: coeff_file    ! 10
    character(len=2),  intent(out) :: ptypes(:,:,:) ! dim(maxsat,maxchan,maxpred) 11
    integer,           intent(out) :: ndata(:,:)    ! dim(maxsat, maxchan), number of channels 12

    ! Locals:
    character(len=8)               :: sat
    character(len=120)             :: line
    integer                        :: chan
    integer                        :: nbfov, nbpred, i, j, k, ier, istat, ii, nobs
    logical                        :: newsat, fileExists
    real                           :: dummy
    integer                        :: iun
    integer                        :: maxsat    
    integer                        :: maxpred 

    !   sats(nsat)            = satellite names
    !   chans(nsat,nchan(i))  = channel numbers of each channel of each satellite i
    !   npred(nsat,nchan(i))  = number of predictors for each channel of each satellite i
    !   fovbias(i,j,k)        = bias for satellite i, channel j, FOV k   k=1,nfov
    !     if FOV not considered for instrument, nfov = 1 and fovbias is global bias for channel
    !   coeff(i,j,1)          = regression constant
    !   coeff(i,j,2), ..., coeff(i,j,npred(i,j)) = predictor coefficients
    !   nsat, nchan, nfov, cinstrum (output) are determined from file
    !   if returned nsat = 0, coeff_file was empty

    inquire(file=trim(coeff_file), exist=fileExists)
    if (fileExists) then
      iun = 0
      ier = fnom(iun, coeff_file, 'FMT', 0)
      if (ier /= 0) then
        call utl_abort('read_coeff: ERROR - Problem opening the coefficient file! ' // trim(coeff_file))
      end if

      write(*,*)
      write(*,*) 'read_coeff: Bias correction coefficient file open = ', coeff_file

    end if

    maxsat =  size(sats)
    maxpred = size(ptypes, dim = 3)

    coeff(:,:,:) = 0.0
    fovbias(:,:,:) = 0.0
    sats(:) = 'XXXXXXXX'
    cinstrum = 'XXXXXXX'
    chans(:,:) = 0
    npred(:,:) = 0
    nsat = 0
    nchan(:) = 0
    nfov = 0
    ptypes(:,:,:) = 'XX'

    if (.not. fileExists) then
      write(*,*) 'read_coeff: Warning- File ' // trim(coeff_file) //'not there.'
      return
    end if

    read(iun,*,iostat=istat)

    if (istat /= 0) then
      write(*,*) 'read_coeff: Warning- File appears empty.'
      ier = fclos(iun)
      return
    end if

    rewind(iun)

    ii = 0

    ! Loop over the satellites/channels in the file

    do
      read(iun,'(A)',iostat=istat) line
      if (istat < 0) exit
      if (line(1:3) == 'SAT') then
        newsat = .true.
        read(line,'(T53,A8,1X,A7,1X,I6,1X,I8,1X,I2,1X,I3)',iostat=istat) sat, cinstrum, chan, nobs, nbpred, nbfov
        if (istat /= 0) then
          call utl_abort('read_coeff: ERROR - reading data from SATELLITE line in coeff file!')
        end if
        do i = 1, maxsat
          if (trim(sats(i)) == trim(sat)) then
            newsat = .false.
            ii = i
          end if
        end do
        if (newsat) then
          ii = ii + 1
          if (ii > maxsat) then
            call utl_abort('read_coeff: ERROR - max number of satellites exceeded in coeff file!')
          end if
          sats(ii) = sat
          if (ii > 1) nchan(ii - 1) = j
          j = 1
        else
          j = j + 1
        end if
        chans(ii,j) = chan
        npred(ii,j) = nbpred
        ndata(ii,j) = nobs
        if (nbpred > maxpred) then
          call utl_abort('read_coeff: ERROR - max number of predictors exceeded in coeff file!')
        end if
        read(iun,'(A)',iostat=istat) line
        if (line(1:3) /= 'PTY') then
          call utl_abort('read_coeff: ERROR - list of predictors is missing in coeff file!')
        end if
        if (nbpred > 0) then
          read(line,'(T8,6(1X,A2))', iostat = istat) (ptypes(ii,j,k), k = 1, nbpred)
          if (istat /= 0) then
            call utl_abort('read_coeff: ERROR - reading predictor types from PTYPES line in coeff file!')
          end if
        end if
        read(iun,*,iostat=istat) (fovbias(ii,j,k), k= 1, nbfov)
        if (istat /= 0) then
          call utl_abort('read_coeff: ERROR - reading fovbias in coeff file!')
        end if
        if (nbpred > 0) then
          read(iun,*,iostat=istat) (coeff(ii,j,k), k = 1, nbpred + 1)
        else
          read(iun,*,iostat=istat) dummy
        end if
        if (istat /= 0) then
          call utl_abort('read_coeff: ERROR - reading coeff in coeff file!')
        end if
      else
        exit
      end if

    end do

    if (ii == 0) then
      call utl_abort('read_coeff: ERROR - No data read from coeff file!')
    end if

    nsat = ii
    nfov = nbfov
    nchan(ii) = j

    write(*,*) ' '
    write(*,*) 'read_coeff: ------------- BIAS CORRECTION COEFFICIENT FILE ------------------ '
    write(*,*) ' '
    write(*,*) 'read_coeff: Number of satellites =     ', nsat
    write(*,*) 'read_coeff: Number of FOV =            ', nfov
    write(*,*) 'read_coeff: Max number of predictors = ', maxval(npred)
    write(*,*) ' '
    do i = 1, nsat
      write(*,*) 'read_coeff:  Satellite = ' // sats(i)
      write(*,*) 'read_coeff:     Number of channels = ', nchan(i)
      write(*,*) 'read_coeff:     predictors, fovbias, coeff for each channel: '
      do j = 1, nchan(i)
        write(*,*) i, chans(i, j)
        if (npred(i, j) > 0) then
          write(*,'(6(1X,A2))') (ptypes(i,j,k), k = 1, npred(i,j))
        else
          write(*,'(A)') 'read_coeff: No predictors'
        end if
        write(*,*) (fovbias(i,j,k), k = 1, nfov)
        write(*,*) (coeff(i,j,k), k = 1, npred(i,j) + 1)
      end do
    end do
    write(*,*) 'read_coeff: exit'

    ier = fclos(iun)

  end subroutine read_coeff

  !-----------------------------------------
  ! bcs_getChannelIndex
  !-----------------------------------------
  subroutine bcs_getChannelIndex(obsSpaceData, idsat, chanIndx,indexBody)
    !
    ! :Purpose: to get the channel index (wrt bcif channels)
    !
    implicit none

    ! Arguments:
    integer,          intent(in)    :: idsat
    integer,          intent(in)    :: indexBody
    integer,          intent(out)   :: chanIndx
    type(struct_obs), intent(inout) :: obsSpaceData

    ! Locals:
    logical, save :: first =.true.
    integer :: ichan, isensor, indx 
    integer, allocatable, save :: Index(:,:)
    
    if (first) then
      allocate(Index(tvs_nsensors,tvs_maxChannelNumber))
      Index(:,:) = -1
      do isensor = 1, tvs_nsensors
        channels:do ichan = 1,  tvs_maxChannelNumber
          indexes: do indx =1, bias(isensor)%numChannels
            if (ichan == bias(isensor)%chans(indx)%channelNum) then
              Index(isensor,ichan) = indx
              exit indexes
            end if
          end do indexes
        end do channels
      end do
      first = .false.
    end if
    
    ichan = nint(obs_bodyElem_r(obsSpaceData, OBS_PPP, indexBody))
    ichan = max(0, min(ichan, tvs_maxChannelNumber + 1))
    ichan = ichan - tvs_channelOffset(idsat)
    
    chanIndx = Index(idsat,ichan)

  end subroutine bcs_getChannelIndex

  !-----------------------------------------
  ! getObsFileName
  !-----------------------------------------
  function getObsFileName(codetype) result(fileName)
    !
    ! :Purpose: Return the part of the observation file name associated
    !           with the type of observation it contains.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: codeType
    ! Result:
    character(len=20)   :: fileName
 
    if (codtyp_get_name(codeType) == 'radianceclear') then
      fileName  = 'csr'
    else if (codtyp_get_name(codeType) == 'mhs' .or. codtyp_get_name(codeType) == 'amsub') then
      fileName = 'to_amsub'
    else if (codtyp_get_name(codeType) == 'amsua') then
      fileName = 'to_amsua'
    else if (codtyp_get_name(codeType) == 'ssmi') then
      fileName = 'ssmis'
    else if (codtyp_get_name(codeType) == 'crisfsr') then
      fileName = 'cris'
    else
      fileName = codtyp_get_name(codeType)
    end if
    
  end function getObsFileName

  !-----------------------------------------
  ! getInitialIdObsData
  !-----------------------------------------
  subroutine getInitialIdObsData(obsSpaceData, obsFamily, idObs, idData, codeTypeList_opt)
    !
    ! :Purpose: Compute initial value for idObs and idData that will ensure
    !           unique values over all mpi tasks
    !
    implicit none

    ! Arguments:
    type(struct_obs),  intent(inout) :: obsSpaceData
    character(len=*),  intent(in)    :: obsFamily    
    integer,           intent(out)   :: idObs
    integer,           intent(out)   :: idData
    integer, optional, intent(in)    :: codeTypeList_opt(:)

    ! Locals:
    integer                :: headerIndex, numHeader, numBody, codeType, ierr
    integer, allocatable   :: allNumHeader(:), allNumBody(:)

    numHeader = 0
    numBody = 0
    call obs_set_current_header_list(obsSpaceData, obsFamily)
    HEADERCOUNT: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADERCOUNT
      if (present(codeTypeList_opt)) then
        codeType  = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
        if (all(codeTypeList_opt(:) /= codeType)) cycle HEADERCOUNT
      end if
      numHeader = numHeader + 1
      numBody = numBody + obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex)
    end do HEADERCOUNT
    allocate(allNumHeader(mmpi_nprocs))
    allocate(allNumBody(mmpi_nprocs))
    call rpn_comm_allgather(numHeader, 1, 'mpi_integer',       &
                            allNumHeader, 1, 'mpi_integer', 'GRID', ierr)
    call rpn_comm_allgather(numBody, 1, 'mpi_integer',       &
                            allNumBody, 1, 'mpi_integer', 'GRID', ierr)
    if (mmpi_myid > 0) then
      idObs = sum(allNumHeader(1:mmpi_myid))
      idData = sum(allNumBody(1:mmpi_myid))
    else
      idObs = 0
      idData = 0
    end if
    deallocate(allNumHeader)
    deallocate(allNumBody)

  end subroutine getInitialIdObsData

end module biasCorrectionSat_mod


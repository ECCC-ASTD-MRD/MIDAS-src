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

MODULE biasCorrection_mod
  ! MODULE biasCorrection_mod (prefix="bias" category='1. High-level functionality')
  !
  ! :Purpose: Performs the variational bias correction for satellite radiance
  !           data
  !
  use utilities_mod
  use ramDisk_mod
  use MathPhysConstants_mod
  use obsSpaceData_mod
  use controlVector_mod
  use mpi_mod
  use mpivar_mod
  use tovs_nl_mod
  use timeCoord_mod
  use columnData_mod
  use codePrecision_mod
  use localizationFunction_mod
  use HorizontalCoord_mod
  use verticalCoord_mod
  use gridStateVector_mod
  use stateToColumn_mod
  use codtyp_mod
  use timeCoord_mod
  use clib_interfaces_mod

  implicit none
  save
  private

  public               :: bias_setup,bias_calcBias_tl,bias_calcBias_ad, bias_writeBias, bias_finalize
  public               :: bias_removeBiasCorrection, bias_refreshBiasCorrection, bias_readConfig
  public               :: bias_do_regression, bias_filterObs, bias_computeResidualsStatistics, bias_calcBias
  public               :: bias_removeOutliers, bias_applyBiasCorrection
  public               :: bias_mimicSatbcor

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
  type(struct_gsv)      :: statevector_mask
  type(struct_columnData) :: column_mask
  logical               :: initialized = .false.
  integer, parameter    :: NumPredictors = 6
  integer, parameter    :: maxfov = 120

  real(8), allocatable  :: trialHeight300m1000(:)
  real(8), allocatable  :: trialHeight50m200(:)
  real(8), allocatable  :: trialHeight1m10(:)
  real(8), allocatable  :: trialHeight5m50(:)
  real(8), allocatable  :: RadiosondeWeight(:)
  real(8), allocatable  :: trialTG(:)
  integer               :: nobs
  character(len=5)      :: biasMode
  logical  :: biasActive, loutstats
  logical  :: doRegression
  logical  :: lMimicSatbcor, lweightedEstimate, filterObs
  logical  :: bias_mimicSatbcor
  real(8)  :: bg_stddev(NumPredictors),predScalingFactor(NumPredictors),predOffset(NumPredictors)
  real(8)  :: scanBiasCorLength 
  logical  :: removeBiasCorrection, refreshBiasCorrection, centerPredictors,loutCoeffCov
  character(len=3) :: cglobal(25)
  character(len=7) :: cinst(25)
  integer :: nbscan(25)
  logical :: offline_opt, all_opt_ssmis, all_opt_tovs, all_opt_csr, all_opt_hyperIr
  integer, external            :: fnom, fclos 
  namelist /nambias/ biasActive,biasMode,bg_stddev,removeBiasCorrection,refreshBiasCorrection
  namelist /nambias/ centerPredictors,doRegression, scanBiasCorLength,  lMimicSatbcor, lweightedEstimate
  namelist /nambias/ cglobal, cinst, nbscan,filterObs ,loutstats,loutCoeffCov
  namelist /nambias/ offline_opt, all_opt_ssmis, all_opt_tovs, all_opt_csr, all_opt_hyperIr

CONTAINS
 
  !-----------------------------------------------------------------------
  ! bias_readConfig
  !-----------------------------------------------------------------------
  subroutine bias_readConfig()
    !
    ! :Purpose: Read nambias namelist section
    !
    implicit none
    !Locals:
    integer  :: ierr,nulnam
  
    ! set default values for namelist variables
    biasActive   = .false.
    biasMode = "varbc"
    bg_stddev(:) = 0.0d0
    removeBiasCorrection = .false.
    filterObs = .false.
    refreshBiasCorrection = .false.
    centerPredictors = .false.
    doRegression = .false.
    lMimicSatbcor = .true.
    scanBiasCorLength = -1.d0
    lweightedEstimate = .false.
    loutCoeffCov = .false.
    nbscan(:) = -1
    cinst(:) = "XXXXXXX"
    cglobal(:) = "XXX"
    loutstats = .false.
    offline_opt=.false.
    all_opt_ssmis=.true.
    all_opt_tovs=.true.
    all_opt_csr=.true.
    all_opt_hyperIr=.false.

    !
    ! read in the namelist NAMBIAS
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=nambias,iostat=ierr)
    if ( ierr /= 0 .and. mpi_myid == 0 )  &
         write(*,*) 'bias_readConfig: WARNING: Error reading namelist, assume it will not be used!'
    if ( mpi_myid == 0 ) write(*,nml=nambias)
    ierr = fclos(nulnam)
    bias_mimicSatbcor = lMimicSatbcor

 end subroutine bias_readConfig


  !-----------------------------------------------------------------------
  ! bias_setup
  !-----------------------------------------------------------------------
  subroutine bias_setup()
    implicit none

    integer  :: cvdim
    integer  :: iSensor,iPredictor, instIndex
    integer  :: iChan
    integer  :: iPred,jPred, kPred,iScan1, iScan2
    character(len=85)  :: bcifFile
    character(len=10)  :: instrName, instrNamecoeff, satNamecoeff 
    logical            :: bcifExists
   
    !variables from background coeff file
    integer            :: nfov, exitCode
    integer            :: nchan(tvs_nSensors)
    character(len=2)   :: predBCIF(tvs_maxchannelnumber,numpredictors)
    integer            :: canBCIF(tvs_maxchannelnumber), npredBCIF(tvs_maxchannelnumber), ncanBcif, npredictors
    character(len=1)   :: bcmodeBCIF(tvs_maxchannelnumber), bctypeBCIF(tvs_maxchannelnumber)
    character(len=3)   :: global
    character(len=128) :: errorMessage
    real(8), allocatable :: Bmatrix(:,:)

    cvdim = 0

    if ( biasActive ) then

      if (scanBiasCorLength > 0.d0) call lfn_Setup('FifthOrder')


      allocate(bias(tvs_nSensors))

      do iSensor = 1, tvs_nSensors

        write(*,*) "bias_setup: iSensor = ",iSensor
       
        instrName = InstrNametoCoeffFileName(tvs_instrumentName(iSensor))
        instrNamecoeff = InstrNameinCoeffFile(tvs_instrumentName(iSensor))
        satNamecoeff = SatNameinCoeffFile(tvs_satelliteName(iSensor)) 

        bcifFile = 'bcif_' // trim(instrName)

        global = "XXX"
        nfov = -1
        do instIndex =1,size(cinst)
          if (trim(instrNamecoeff) == trim(cinst(instIndex))) then
            global = cglobal(instIndex)
            nfov = nbscan(instIndex)
          end if
        end do
        if ( nfov == -1) then
          write(*,*) "bias_setup: Problem with instrName ",instrNamecoeff
          write(*,'(15(A10,1x))')  cinst(:)
          call utl_abort('bias_setup: check nambias namelist')
        end if

        inquire(file=trim(bcifFile),exist = bcifExists)
        if ( bcifExists ) then
          !read bcif (Bias Correction Information File)
         
          call read_bcif(bcifFile, tvs_isNameHyperSpectral( instrName ),ncanBcif, &
               canBCIF, bcmodeBCIF, bctypeBCIF, npredBCIF, predBCIF, global, exitcode)

          if (exitcode /= 0) then
            call utl_abort("bias_setup: Problem in read_bcif while reading " // trim(bcifFile) )
          end if

          bias(iSensor)%numChannels = ncanBcif

          allocate( bias(iSensor)%chans(ncanBcif) )
          
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
            allocate(  bias(iSensor)%chans(ichan)% predictorIndex( npredictors ) )
            bias(iSensor) % chans(ichan) % stddev(:) = 0.d0
            bias(iSensor) % chans(ichan) % coeffIncr(:) = 0.d0
            bias(iSensor) % chans(ichan) % coeff(:) = 0.d0
            bias(iSensor) % chans(ichan) % coeff_offset(:) = 0.d0

            bias(iSensor)%chans(ichan)% predictorIndex(1) = 1 !the constant term is always included
            jPred = 1
            do ipred = 1, npredBCIF(ichan + 1)
              jPred =  jPred + 1
              select case(predBCIF(ichan+1,ipred))
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
              case default
                write(errorMessage,*) "bias_setup: Unknown predictor ", predBCIF(ichan+1,ipred), ichan, ipred
                call utl_abort(errorMessage)
              end select
              bias(iSensor)%chans(ichan)% predictorIndex(jPred) = kpred
            end do
          end do
        else
          call utl_abort( "bias_setup: Error : " // trim(bcifFile) // " not present !" )
        end if

        bias(iSensor) % numscan = nfov 

        allocate( bias(iSensor) % BHalfScanBias (nfov,nfov))
        if (doRegression) allocate( bias(iSensor) % BMinusHalfScanBias (nfov,nfov))
        allocate( Bmatrix(nfov,nfov))
        do ichan=1, ncanBcif
          allocate( bias(iSensor) % chans(ichan) % coeffIncr_fov(nfov))
          allocate( bias(iSensor) % chans(ichan) % coeff_fov(nfov))
          bias(iSensor) % chans(ichan) % coeffIncr_fov(:) = 0.d0
          bias(iSensor) % chans(ichan) % coeff_fov(:) = 0.d0
        end do

        do ichan =1, ncanBcif
          if ( bias(iSensor) % chans(ichan) %isDynamic ) then
            do iPredictor = 1, bias(iSensor)% chans(ichan) % numActivePredictors
              bias(iSensor)% chans(ichan) % stddev(iPredictor) = bg_stddev( bias(iSensor)% chans(ichan) % PredictorIndex(iPredictor) )
            end do
          end if
        end do

        if ( trim(biasMode) == "varbc" ) then
          !change dimension of control vector
          do  ichan=1, ncanBcif
            if  (bias(iSensor)% chans(ichan) % isDynamic) &
                 cvdim = cvdim + bias(iSensor)% chans(ichan) % numActivePredictors - 1 + bias(iSensor)%numScan
          end do
        end if

        if (allocated(Bmatrix)) then
          if (scanBiasCorLength > 0.d0) then
            do iScan2=1,nfov
              do iScan1=1,nfov
                Bmatrix(iScan1,iScan2) =   bg_stddev(1) * bg_stddev(1) * lfn_Response(1.d0*abs(iScan1-iScan2),scanBiasCorLength)
              end do
            end do
          else
            Bmatrix(:,:)=0.d0
            do iScan1=1,nfov
              Bmatrix(iScan1,iScan1) =   bg_stddev(1) * bg_stddev(1)
            end do
          end if
          bias(iSensor) %BHalfScanBias(:,:) =  Bmatrix(:,:)
          call utl_matsqrt(bias(iSensor) %BHalfScanBias,nfov,1.d0,printInformation_opt=.true.)
          if (doRegression) then
            bias(iSensor) %BMinusHalfScanBias(:,:) =  Bmatrix(:,:)
            call utl_matsqrt(bias(iSensor) %BMinusHalfScanBias,nfov,-1.d0,printInformation_opt=.true.)
          end if
          deallocate(Bmatrix)
        end if

      end do

    end if

    if  ( trim(biasMode) == "varbc" .and. cvdim > 0 ) then
      if ( mpi_myid > 0 ) cvdim = 0 ! for minimization, all coefficients only on task 0
      call cvm_setupSubVector('BIAS', 'BIAS', cvdim)
    end if

    call  bias_readCoeffs() ! Read coefficient files in the case of bias correction application (biasMode=="apply")

  end subroutine bias_setup

  !-----------------------------------------------------------------------
  ! bias_readCoeffs
  !-----------------------------------------------------------------------
  subroutine bias_readCoeffs()
    !
    ! :Purpose: Fill the bias structure with read static and dynamic bias correction coefficient files
    !
    implicit none
    !Locals:
    integer :: iSensor, iSat, jchannel, jChan
    integer :: satIndexDynamic, satIndexStatic
    integer :: chanindexDynamic, chanindexStatic
    character(len=10) :: instrName, instrNamecoeff, satNamecoeff
    character(len=64) :: dynamicCoeffFile, staticCoeffFile
    logical           :: corrected
    integer           :: nfov, npredictors

    character(len=10) :: satsDynamic(tvs_nsensors)       ! dim(maxsat), satellite names 1
    integer           :: chansDynamic(tvs_nsensors,tvs_maxchannelnumber)    ! dim(maxsat, maxchan), channel numbers 2
    real(8)           :: fovbiasDynamic(tvs_nsensors,tvs_maxchannelnumber,maxfov)! dim(maxsat,maxchan,maxfov), bias as F(fov) 3
    real(8)           :: coeffDynamic(tvs_nsensors,tvs_maxchannelnumber,NumPredictors+1)  ! dim(maxsat,maxchan,maxpred+1) 4
    integer           :: nsatDynamic          !5
    integer           :: nchanDynamic(tvs_nsensors)      ! dim(maxsat), number of channels 6
    integer           :: nfovDynamic          !7
    integer           :: npredDynamic(tvs_nsensors,tvs_maxchannelnumber)    ! dim(maxsat, maxchan), number of predictors !8
    character(len=7)  :: cinstrumDynamic      ! string: instrument (e.g. AMSUB) 9
    character(len=2)  :: ptypesDynamic(tvs_nsensors,tvs_maxchannelnumber,NumPredictors) ! dim(maxsat,maxchan,maxpred) 11
    integer           :: ndataDynamic(tvs_nsensors,tvs_maxchannelnumber)    ! dim(maxsat, maxchan), number of channels 12

    character(len=10) :: satsStatic(tvs_nsensors)       ! dim(maxsat), satellite names 1
    integer           :: chansStatic(tvs_nsensors,tvs_maxchannelnumber)    ! dim(maxsat, maxchan), channel numbers 2
    real(8)           :: fovbiasStatic(tvs_nsensors,tvs_maxchannelnumber,maxfov)! dim(maxsat,maxchan,maxfov), bias as F(fov) 3
    real(8)           :: coeffStatic(tvs_nsensors,tvs_maxchannelnumber,NumPredictors+1)  ! dim(maxsat,maxchan,maxpred+1) 4
    integer           :: nsatStatic          !5
    integer           :: nchanStatic(tvs_nsensors)      ! dim(maxsat), number of channels 6
    integer           :: nfovStatic          !7
    integer           :: npredStatic(tvs_nsensors,tvs_maxchannelnumber)    ! dim(maxsat, maxchan), number of predictors !8
    character(len=7)  :: cinstrumStatic      ! string: instrument (e.g. AMSUB) 9
    character(len=2)  :: ptypesStatic(tvs_nsensors,tvs_maxchannelnumber,NumPredictors) ! dim(maxsat,maxchan,maxpred) 11
    integer           :: ndataStatic(tvs_nsensors,tvs_maxchannelnumber)    ! dim(maxsat, maxchan), number of channels 12


    if (biasActive .and. biasMode == "apply") then

      ! 1 fichier de coefficient par intrument avec les differentes plateformes
      ! Cas particulier GEORAD (CSR) 
      
      do iSensor = 1, tvs_nSensors

        write(*,*) "bias_readCoeffs: iSensor = ",iSensor
       
        instrName = InstrNametoCoeffFileName(tvs_instrumentName(iSensor))
        instrNamecoeff = InstrNameinCoeffFile(tvs_instrumentName(iSensor))
        satNamecoeff = SatNameinCoeffFile(tvs_satelliteName(iSensor)) 

        dynamicCoeffFile = "coeff_file_" // trim(instrName)
        staticCoeffFile = "static_coeff_file_" // trim(instrName)

        if (  tvs_isNameGeostationary(instrName) ) then
          dynamicCoeffFile = trim(dynamicCoeffFile) // "." // trim( satNamecoeff )
          staticCoeffFile = trim(staticCoeffFile) // "." // trim( satNamecoeff ) 
        end if

        call read_coeff(satsDynamic, chansDynamic, fovbiasDynamic, coeffDynamic, nsatDynamic, nchanDynamic, nfovDynamic, &
             npredDynamic, cinstrumDynamic, dynamicCoeffFile, ptypesDynamic,ndataDynamic)

        call read_coeff(satsStatic, chansStatic, fovbiasStatic, coeffStatic, nsatStatic, nchanStatic, nfovStatic, &
             npredStatic, cinstrumStatic, staticCoeffFile, ptypesStatic,ndataStatic)
        write(*,*) "bias_readCoeffs: cinstrumDynamic= ", cinstrumDynamic
        write(*,*) "bias_readCoeffs: cinstrumStatic= ", cinstrumStatic

        satIndexDynamic = -1
        do iSat = 1, nsatDynamic
          if ( trim(satNameCoeff) /= trim(satsDynamic(iSat)) .or. trim(instrNamecoeff) /= trim(cinstrumDynamic) ) cycle
          satIndexDynamic = iSat
        end do

        satIndexStatic = -1
        do iSat = 1, nsatStatic
          if ( trim(satNameCoeff) /= trim(satsStatic(iSat)) .or. trim(instrNamecoeff) /= trim(cinstrumStatic) ) cycle
          satIndexStatic = iSat
        end do

        nfov =  bias(iSensor) % numscan

        do jChannel = 1, bias(iSensor)%numChannels

          chanindexDynamic = -1
          if ( satIndexDynamic > 0) then
            do jChan = 1, nchanDynamic(satIndexDynamic)
              if ( chansDynamic(satIndexDynamic,jChan) == bias(iSensor)%chans(jChannel)%channelNum ) then
                chanindexDynamic = jChan
                exit
              end if
            end do
          end if

          chanindexStatic = -1
          if ( satIndexStatic > 0) then
            do jChan = 1, nchanStatic(satIndexStatic)
              if ( chansStatic(satIndexStatic,jChan) == bias(iSensor)%chans(jChannel)%channelNum ) then
                chanindexStatic = jChan
                exit
              end if
            end do
          end if
          npredictors = bias(iSensor) % chans(jChannel) % numActivePredictors
       
          corrected =.true.
          select case(bias(iSensor) % chans(jChannel) % bcmode)
          case("D")
            if ( chanindexDynamic > 0) then
              if (bias(iSensor) % chans(jChannel) % bctype=="C") &
                   bias(iSensor) % chans(jChannel) % coeff(1:npredictors) = coeffDynamic(satIndexDynamic,chanindexDynamic,1:npredictors)
              if (bias(iSensor) % chans(jChannel) % bctype=="C" .or. bias(iSensor) % chans(jChannel) % bctype=="F") &
                   bias(iSensor) % chans(jChannel) % coeff_fov(1:nfov) = fovbiasDynamic(satIndexDynamic,chanindexDynamic,1:nfov)
            else if ( chanindexStatic > 0) then
              if (bias(iSensor) % chans(jChannel) % bctype=="C") &
                   bias(iSensor) % chans(jChannel) % coeff(1:npredictors) = coeffStatic(satIndexStatic,chanindexStatic,1:npredictors)
              if (bias(iSensor) % chans(jChannel) % bctype=="C" .or. bias(iSensor) % chans(jChannel) % bctype=="F") &
                   bias(iSensor) % chans(jChannel) % coeff_fov(1:nfov) = fovbiasStatic(satIndexStatic,chanindexStatic,1:nfov)
            else
              corrected = .false.
            end if
          case("S")
            if ( chanindexStatic > 0) then
              if (bias(iSensor) % chans(jChannel) % bctype=="C") &
                   bias(iSensor) % chans(jChannel) % coeff(1:npredictors) = coeffStatic(satIndexStatic,chanindexStatic,1:npredictors)
              if (bias(iSensor) % chans(jChannel) % bctype=="C" .or. bias(iSensor) % chans(jChannel) % bctype=="F") & 
                   bias(iSensor) % chans(jChannel) % coeff_fov(1:nfov) = fovbiasStatic(satIndexStatic,chanindexStatic,1:nfov)
            else
              corrected = .false.
            end if
          end select
          
          if (.not. corrected) then
            Write(*,*) "bias_readCoeffs: Warning: channel ",  bias(iSensor) % chans(jChannel) % channelNum, " of ", &
                 trim( instrName )," ",trim( satNamecoeff )," not corrected!"
          end if

        end do

      end do
        
    end if

  end subroutine bias_readCoeffs


  !-----------------------------------------------------------------------
  ! bias_computePredictorBiases
  !-----------------------------------------------------------------------
  subroutine bias_computePredictorBiases(obsSpaceData)
    !
    ! :Purpose: to compute predictor average
    !
    implicit none
    !Arguments:
    type(struct_obs), intent(inout)  :: obsSpaceData
    !Locals:
    real(8) :: predictor(NumPredictors)
    integer :: iobs,iChannel,nsize,i,j,npred
    integer :: headerIndex, idatyp, indxtovs
    integer :: iSensor,iFov,iPredictor,ierr
    integer :: bodyIndex, jpred, chanIndx
    real(8), allocatable ::  temp_offset(:,:)
    integer, allocatable ::  temp_nobs(:)
    real(8), allocatable ::  temp_offset2(:,:,:)
    integer, allocatable ::  temp_nobs2(:,:)

    if (centerPredictors) then
      write(*,*) "bias_computePredictorBiases: start"
      npred = 0
      do iSensor=1, tvs_nSensors
        npred = max(npred,maxval(bias(iSensor)%chans(:)%numActivePredictors))
      end do
      allocate( temp_offset2( tvs_nsensors, maxval(bias(:)%numChannels), 2:npred ) )
      allocate( temp_nobs2( tvs_nsensors, maxval(bias(:)%numChannels)) )
      temp_offset2(:,:,:) = 0.d0
      temp_nobs2(:,:) = 0

      call obs_set_current_header_list(obsSpaceData,'TO')
      iobs = 0
      HEADER: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if ( headerIndex < 0 ) exit HEADER

        ! process only radiance data to be assimilated?
        idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
        if ( .not. tvs_isIdBurpTovs(idatyp) ) cycle HEADER

        iobs = iobs + 1
        
        indxTovs = tvs_tovsIndex(headerIndex)
        if ( indxTovs < 0 ) cycle HEADER

        iSensor = tvs_lsensor( indxTovs )

        call obs_set_current_body_list(obsSpaceData, headerIndex)
        iFov = obs_headElem_i(obsSpaceData,OBS_FOV,headerIndex)

        BODY: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if ( bodyIndex < 0 ) exit BODY

          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY   

          call bias_getChannelIndex(obsSpaceData,iSensor,chanIndx,bodyIndex)
          if (chanindx > 0) then
            call bias_getPredictors(predictor,headerIndex,iobs,chanIndx,obsSpaceData)
            do iPredictor = 2, bias(iSensor)%chans(chanIndx)%NumActivePredictors
              jPred = bias(iSensor)%chans(chanIndx)%PredictorIndex(iPredictor)
              temp_offset2(iSensor,chanIndx,iPredictor) =  temp_offset2(iSensor,chanIndx,iPredictor) + predictor(jPred)
            end do
            temp_nobs2(iSensor,chanIndx) =  temp_nobs2(iSensor,chanIndx) + 1
          end if
        end do BODY
      end do HEADER

      allocate( temp_offset(maxval(bias(:)%numChannels), 2:npred))
      allocate( temp_nobs( maxval(bias(:)%numChannels) ) )

      do iSensor =1,tvs_nSensors 
        temp_offset(:,:) = 0.0d0
        nsize=size( temp_offset )
        call rpn_comm_allreduce(temp_offset2(iSensor,:,:),temp_offset(:,:),nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
        if ( ierr /= 0) then
          call utl_abort("bias_computePredictorBiases: Erreur de communication MPI 1")
        end if
       
        do i=1, bias(iSensor)%numChannels
          do j=2,bias(iSensor)%chans(i)%numActivePredictors
            bias(iSensor)%chans(i)%coeff_offset(j) = temp_offset(i,j)
          end do
        end do

        temp_nobs(:) = 0
        nsize=size( temp_nobs )
        call rpn_comm_allreduce(temp_nobs2(iSensor,:),temp_nobs,nsize,"mpi_integer","mpi_sum","GRID",ierr)
        if ( ierr /= 0) then
          call utl_abort("bias_computePredictorBiases: Erreur de communication MPI 2")
        end if

       
        do i=1, bias(iSensor)%numChannels
          bias(iSensor)%chans(i)%coeff_nobs = temp_nobs(i)
        end do
        
        do i=1, bias(iSensor)%numChannels
          if (bias(iSensor)%chans(i)%coeff_nobs > 0) then
            bias(iSensor)%chans(i)%coeff_offset(:) = bias(iSensor)%chans(i)%coeff_offset / bias(iSensor)%chans(i)%coeff_nobs
          end if
        end do

      end do

      deallocate( temp_offset )
      deallocate( temp_nobs )
      deallocate( temp_offset2 )
      deallocate( temp_nobs2 )

      write(*,*) "bias_computePredictorBiases: end"
    end if
   

  end subroutine bias_computePredictorBiases

  !-----------------------------------------------------------------------
  ! bias_calcBias
  !-----------------------------------------------------------------------
  subroutine bias_calcBias(obsSpaceData,columnhr)
    !
    ! :Purpose:  to fill OBS_BCOR column of ObsSpaceData body with bias correction computed from read coefficient file
    !
    implicit none
    !Arguments:
    type(struct_obs)        :: obsSpaceData
    type(struct_columnData) :: columnhr
    !Locals:
    integer  :: headerIndex,bodyIndex,iobs, indxtovs, idatyp
    integer  :: iSensor,iPredictor,index_cv,chanIndx
    integer  :: iScan, iFov, jPred
    real(8)  :: predictor(NumPredictors)
    real(8)  :: biasCor

    if ( .not. biasActive ) return

    write(*,*) "bias_calcBias: start"

    if ( .not. allocated(trialHeight300m1000) ) then
      call bias_getTrialPredictors(obsSpaceData,columnhr)
    end if

    iobs = 0
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if ( .not. tvs_isIdBurpTovs(idatyp) ) cycle HEADER

      iobs = iobs + 1
      
      indxtovs = tvs_tovsIndex(headerIndex)
      if ( indxtovs < 0 ) cycle HEADER

      iSensor = tvs_lsensor( indxTovs )

      call obs_set_current_body_list(obsSpaceData, headerIndex)
      iFov = obs_headElem_i(obsSpaceData,OBS_FOV,headerIndex)

      if ( bias(iSensor)%numScan > 1 ) then
        iScan = iFov
      else
        iScan = 1
      end if


      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY   

        call bias_getChannelIndex(obsSpaceData,iSensor,chanIndx,bodyIndex)
        if (chanindx > 0) then
          biasCor = 0.0d0
          if (bias(iSensor)%chans(chanIndx)%isDynamic .and. bias(iSensor)%numScan >0) then
            call bias_getPredictors(predictor, headerIndex, iobs, chanIndx, obsSpaceData)
            biasCor = bias(iSensor)%chans(chanIndx)%coeff_fov(iScan) + &
                 bias(iSensor)%chans(chanIndx)%coeff(1) 
            do iPredictor = 2, bias(iSensor)%chans(chanIndx)%NumActivePredictors
              jPred = bias(iSensor)%chans(chanIndx)%PredictorIndex(iPredictor)
              biasCor = biasCor + predictor(jPred) * bias(iSensor)%chans(chanIndx)%coeff(iPredictor)
            end do
          end if

          biasCor = -1.d0 * biascor

          call obs_bodySet_r( obsSpaceData, OBS_BCOR, bodyIndex, biasCor)

        end if
      end do BODY
    end do HEADER

    write(*,*) "bias_calcBias: end"

  end subroutine bias_calcBias


  !---------------------------------------
  ! bias_computeResidualsStatistics
  !----------------------------------------
  subroutine bias_computeResidualsStatistics(obsSpaceData,prefix)
    !
    ! :Purpose: to compute residuals mean and standard deviation by intrument, channel and scan position
    !
    implicit none
    !Arguments:
    type(struct_obs), intent(inout)  :: obsSpaceData
    character(len=*), intent(in)     :: prefix
    !Locals:
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

    if ( .not. biasActive ) return

    if (.not. loutstats) return

    write(*,*) "bias_computeResidualsStatistics: start"


    SENSORS:do sensorIndex = 1, tvs_nsensors

      if  (.not. tvs_isReallyPresentMpiGLobal(sensorIndex)) cycle SENSORS

      write(*,*) "bias_computeResidualsStatistics: sensorIndex ",  sensorIndex

      nchans = bias(sensorIndex)%numChannels
      nscan   = bias(sensorIndex)% numscan
     
      allocate(tbias(nchans,nscan) )
      tbias(:,:) = 0.d0
      allocate( tstd(nchans,nscan) )
      tstd(:,:) = 0.d0
      allocate( tcount(nchans,nscan) )
      tcount(:,:) = 0

      call obs_set_current_header_list(obsSpaceData,'TO')

      HEADER: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if ( headerIndex < 0 ) exit HEADER
        if ( tvs_tovsIndex(headerIndex) < 0) cycle HEADER
        iSensor = tvs_lsensor(tvs_tovsIndex(headerIndex))
        if (iSensor /= sensorIndex) cycle HEADER
          
        iFov = obs_headElem_i(obsSpaceData,OBS_FOV,headerIndex)
        if ( nscan > 1 ) then
          iScan = iFov
        else
          iScan = 1
        end if

        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if ( bodyIndex < 0 ) exit BODY
          
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY 
          call bias_getChannelIndex(obsSpaceData,iSensor,chanIndx,bodyIndex)
          if (chanindx > 0) then
            OmF = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
            bcor =  obs_bodyElem_r(obsSpaceData,OBS_BCOR,bodyIndex)
            if (OmF /= MPC_missingValue_R8 .and. bcor /= MPC_missingValue_R8) then
              tbias(chanIndx,iScan) = tbias(chanIndx,iScan) + ( OmF + bcor )
              tstd(chanIndx,iScan) = tstd(chanIndx,iScan) + ( OmF + bcor ) ** 2
              tcount(chanIndx,iScan) =  tcount(chanIndx,iScan) + 1
            end if
          end if
        end do BODY
      end do HEADER

      allocate( biasMpiGlobal(nchans,nscan) )
      allocate( stdMpiGlobal(nchans,nscan) )
      allocate( countMpiGlobal(nchans,nscan) )

      call rpn_comm_reduce(tbias , biasMpiGlobal, size(biasMpiGlobal), "MPI_DOUBLE_PRECISION" , "MPI_SUM", 0, "GRID", ierr )
      if (ierr /= 0) then
        Write(errorMessage,*) "bias_computeResidualsStatistics: MPI communication error 1",  ierr 
        call utl_abort(errorMessage)
      end if

      call rpn_comm_reduce(tstd , stdMpiGlobal, size(stdMpiGlobal), "MPI_DOUBLE_PRECISION" , "MPI_SUM", 0, "GRID", ierr )
      if (ierr /=0 ) then
        Write(errorMessage,*) "bias_computeResidualsStatistics: MPI communication error 2",  ierr 
        call utl_abort(errorMessage)
      end if

      call rpn_comm_reduce(tcount , countMpiGlobal, size(countMpiGlobal), "MPI_INTEGER" , "MPI_SUM", 0, "GRID", ierr )
      if (ierr /=0 ) then
        Write(errorMessage,*) "bias_computeResidualsStatistics: MPI communication error 3",  ierr 
        call utl_abort(errorMessage)
      end if

      if ( mpi_myId == 0 ) then
        where( countMpiGlobal > 0 ) 
          biasMpiGlobal = biasMpiGlobal / countMpiGlobal
          stdMpiGlobal = sqrt( stdMpiGlobal/ countMpiGlobal  - biasMpiGlobal**2)
        end where
      
        instrName = InstrNametoCoeffFileName(tvs_instrumentName(sensorIndex))
        satNamecoeff = SatNameinCoeffFile(tvs_satelliteName(sensorIndex)) 

        nulfile1 = 0
        ierr = fnom(nulfile1, './std_' // trim(instrName) // '_' // trim(satNamecoeff) // trim(prefix) // '.dat', 'FTN+FMT',0)
        nulfile2 = 0
        ierr = fnom(nulfile2, './mean_' // trim(instrName) // '_' // trim(satNamecoeff) // trim(prefix) // '.dat', 'FTN+FMT',0)

        do chanIndx=1, nchans
          if ( sum( countMpiGlobal(chanIndx,:)) > 0 ) then
            write(nulfile2,'(i4,1x,100e14.6)') chanIndx,  biasMpiGlobal(chanindx,:)
            write(nulfile1,'(i4,1x,100e14.6)') chanIndx,  stdMpiGlobal(chanindx,:)
          end if
        end do
        ierr = fclos(nulfile1)
        ierr = fclos(nulfile2)
       
      end if

      deallocate( biasMpiGlobal )
      deallocate( stdMpiGlobal )
      deallocate( countMpiGlobal )

      deallocate(tbias)
      deallocate( tstd )
      deallocate( tcount )
 
    end do SENSORS

    write(*,*) "bias_computeResidualsStatistics: end"
    
  end subroutine bias_computeResidualsStatistics

  !---------------------------------------
  !  bias_removeOutliers
  !---------------------------------------- 
  subroutine  bias_removeOutliers(obsSpaceData)
    !
    ! :Purpose: to remove outliers (too large OmF) from linear regression
    !
    implicit none
    !Arguments:
    type(struct_obs)  :: obsSpaceData
    !Locals:
    real(8), allocatable :: tbias(:,:), tstd(:,:)
    integer, allocatable :: tcount(:,:)
    real(8), allocatable :: biasMpiGlobal(:,:), stdMpiGLobal(:,:)
    integer, allocatable :: countMpiGlobal(:,:)
    integer :: sensorIndex, headerIndex, bodyIndex
    integer :: nchans, nfiles
    integer :: iSensor, chanIndx, timeIndex
    real(8) ::  OmF
    real(8), parameter :: alpha = 5.d0
    real(8) :: stepObsIndex
    integer :: ierr
    character(len=72) :: errorMessage

    if ( .not. biasActive ) return

    write(*,*) "bias_removeOutliers: start"

    SENSORS:do sensorIndex = 1, tvs_nsensors

      if  (.not. tvs_isReallyPresentMpiGLobal(sensorIndex)) cycle SENSORS

      write(*,*) "bias_removeOutliers: sensorIndex ",  sensorIndex

      nchans = bias(sensorIndex)%numChannels
      nfiles = tim_nstepobs
     
      allocate(tbias(nchans,nfiles) )
      tbias(:,:) = 0.d0
      allocate( tstd(nchans,nfiles) )
      tstd(:,:) = 0.d0
      allocate( tcount(nchans,nfiles) )
      tcount(:,:) = 0

      call obs_set_current_header_list(obsSpaceData,'TO')

      HEADER1: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if ( headerIndex < 0 ) exit HEADER1
        if ( tvs_tovsIndex(headerIndex) < 0 ) cycle HEADER1
        iSensor = tvs_lsensor(tvs_tovsIndex(headerIndex))
        if (iSensor /= sensorIndex) cycle HEADER1

        call tim_getStepObsIndex(stepObsIndex, tim_getDatestamp(), &
             obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex), &
             obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex),tim_nstepobs)

        timeIndex = nint( stepObsIndex )
        if  ( timeIndex < 0 ) cycle HEADER1
        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY1: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if ( bodyIndex < 0 ) exit BODY1
          
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated .and. &
               .not. btest( obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),6) ) then 
            call bias_getChannelIndex(obsSpaceData,iSensor,chanIndx,bodyIndex)
            if (chanindx > 0) then
              OmF = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
              if ( OmF /= MPC_missingValue_R8) then
                tbias(chanIndx,timeindex) = tbias(chanIndx,timeindex) +  OmF
                tstd(chanIndx,timeindex) = tstd(chanIndx,timeindex) + ( OmF ) ** 2
                tcount(chanIndx,timeindex) =  tcount(chanIndx,timeindex) + 1
              end if
            end if
          end if
        end do BODY1
      end do HEADER1

      allocate( biasMpiGlobal(nchans,nfiles) )
      allocate( stdMpiGlobal(nchans,nfiles) )
      allocate( countMpiGlobal(nchans,nfiles) )

      call rpn_comm_reduce(tbias , biasMpiGlobal, size(biasMpiGlobal), "MPI_DOUBLE_PRECISION" , "MPI_SUM", 0, "GRID", ierr )
      if (ierr /=0 ) then
        Write(errorMessage,*) "bias_removeOutliers: MPI communication error 1",  ierr 
        call utl_abort(errorMessage)
      end if

      call rpn_comm_reduce(tstd , stdMpiGlobal, size(stdMpiGlobal), "MPI_DOUBLE_PRECISION" , "MPI_SUM", 0, "GRID", ierr )
      if (ierr /=0 ) then
        Write(errorMessage,*) "bias_removeOutliers: MPI communication error 2",  ierr 
        call utl_abort(errorMessage)
      end if

      call rpn_comm_reduce(tcount , countMpiGlobal, size(countMpiGlobal), "MPI_INTEGER" , "MPI_SUM", 0, "GRID", ierr )
      if (ierr /=0 ) then
        Write(errorMessage,*) "bias_removeOutliers: MPI communication error 3",  ierr 
        call utl_abort(errorMessage)
      end if

      if ( mpi_myId == 0 ) then
        where( countMpiGlobal > 0 ) 
          biasMpiGlobal = biasMpiGlobal / countMpiGlobal
          stdMpiGlobal = sqrt( stdMpiGlobal/ countMpiGlobal  - biasMpiGlobal**2)
        end where
      end if

      call rpn_comm_bcast(countMpiGlobal,nchans*nfiles,"mpi_integer",0,"GRID",ierr)
      if (ierr /=0) then
        Write(errorMessage,*) "bias_removeOutliers: MPI communication error 4",  ierr 
        call utl_abort(errorMessage)
      end if
      call rpn_comm_bcast(stdMpiGlobal,nchans*nfiles,"mpi_double_precision",0,"GRID",ierr)
      if (ierr /=0) then
        Write(errorMessage,*) "bias_removeOutliers: MPI communication error 5",  ierr 
        call utl_abort(errorMessage)
      end if

      if ( sum(countMpiGlobal)/=0 ) then

        tcount(:,:) = 0

        call obs_set_current_header_list(obsSpaceData,'TO')

        HEADER2: do
          headerIndex = obs_getHeaderIndex(obsSpaceData)
          if ( headerIndex < 0 ) exit HEADER2
          if (tvs_tovsIndex(headerIndex) < 0) cycle HEADER2
          iSensor = tvs_lsensor(tvs_tovsIndex(headerIndex))
          if (iSensor /= sensorIndex) cycle HEADER2

          call tim_getStepObsIndex(stepObsIndex, tim_getDatestamp(), &
               obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex), &
               obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex),tim_nstepobs)

          timeIndex = nint( stepObsIndex )
          if  ( timeIndex <0 ) cycle HEADER2
          
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY2: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if ( bodyIndex < 0 ) exit BODY2
          
            if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated .and. &
                 .not. btest( obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),6) ) then 

              call bias_getChannelIndex(obsSpaceData,iSensor,chanIndx,bodyIndex)

              if (chanindx > 0) then
             
                OmF = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
                if (countMpiGlobal(chanIndx,timeindex) >2 .and.  &
                     abs(OmF) > alpha * stdMpiGlobal(chanIndx,timeindex) ) then
                  call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
                  tcount(chanIndx,timeindex) = tcount(chanIndx,timeindex) +1
                end if
             
              end if
            end if
          end do BODY2
        end do HEADER2

        do chanIndx =1, nchans
          if ( sum( tcount(chanIndx,:)) > 0) then
            Write(*,'(A,1x,i2,1x,i4,1x,i6,1x,30e14.6)') "bias_removeOutliers:", sensorindex, chanIndx, sum( tcount(chanIndx,:)), stdMpiGlobal(chanIndx,:)
          end if
        end do

      end if

      deallocate( biasMpiGlobal )
      deallocate( stdMpiGlobal )
      deallocate( countMpiGlobal )

      deallocate( tbias )
      deallocate( tstd )
      deallocate( tcount )

    end do SENSORS

    write(*,*) "bias_removeOutliers: end"

  end subroutine bias_removeOutliers



  !---------------------------------------
  ! bias_calcBias_tl
  !---------------------------------------- 
  subroutine bias_calcBias_tl(cv_in,obsColumnIndex,obsSpaceData,columnhr)
    !
    ! :Purpose: tl of bias computation (for varBC)
    !
    implicit none
    !Arguments:
    real(8)  :: cv_in(:)
    integer  :: obsColumnIndex
    type(struct_obs)  :: obsSpaceData
    type(struct_columnData) :: columnhr
    !Locals:
    integer  :: headerIndex,bodyIndex,iobs, indxtovs, idatyp
    integer  :: iSensor,iPredictor,index_cv,chanIndx
    integer  :: iScan, iFov, jPred
    real(8)  :: predictor(NumPredictors)
    real(8), pointer :: cv_bias(:)
    real(8), target  :: dummy4Pointer(1)
    real(8)  :: biasCor

    if ( .not. biasActive ) return

    if ( .not. allocated(trialHeight300m1000) ) then
      call bias_getTrialPredictors(obsSpaceData,columnhr)
      call bias_computePredictorBiases(obsSpaceData)
      call bias_getRadiosondeWeight(obsSpaceData, lmodify_obserror_opt=.true.)
    end if

    nullify(cv_bias)
    if ( mpi_myid == 0) then
      if ( cvm_subVectorExists('BIAS') ) then
        cv_Bias => cvm_getSubVector(cv_in,'BIAS')
        write(*,*) 'bias_calcBias_tl: maxval(cv_bias)=',maxval(cv_bias(:))
      else
        write(*,*) 'bias_calcBias_tl: control vector does not include bias coefficients'
        return
      end if
   else
      cv_bias => dummy4Pointer
   end if

    ! get bias coefficients
    call bias_cvToCoeff(cv_bias)

    ! apply bias increment to specified obs column
    iobs = 0
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if ( .not. tvs_isIdBurpTovs(idatyp) ) cycle HEADER
      
      iobs = iobs + 1

      indxtovs = tvs_tovsIndex(headerIndex)
      if ( indxtovs < 0 ) cycle HEADER

      iSensor = tvs_lsensor( indxTovs )

      call obs_set_current_body_list(obsSpaceData, headerIndex)
      iFov = obs_headElem_i(obsSpaceData,OBS_FOV,headerIndex)

      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY   

        call bias_getChannelIndex(obsSpaceData,iSensor,chanIndx,bodyIndex)

        if (chanindx > 0) then
          biasCor = 0.0d0
          if (bias(iSensor)%chans(chanIndx)%isDynamic .and. bias(iSensor)%numScan >0) then
            call bias_getPredictors(predictor,headerIndex,iobs,chanindx,obsSpaceData)
            do iPredictor = 1, bias(iSensor)%chans(chanIndx)%NumActivePredictors
              jPred = bias(iSensor)%chans(chanIndx)%PredictorIndex(iPredictor)
              if ( iPredictor == 1 ) then
                if ( bias(iSensor)%numScan > 1 ) then
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
          call obs_bodySet_r( obsSpaceData, obsColumnIndex, bodyIndex, &
               obs_bodyElem_r(obsSpaceData,obsColumnIndex,bodyIndex) - biasCor)
        end if
      end do BODY
    end do HEADER


  end subroutine bias_calcBias_tl

 
  !----------------------
  ! bias_getTrialPredictors
  !----------------------
  subroutine bias_getTrialPredictors(obsSpaceData,columnhr)
    !
    ! :Purpose: get predictors from trial fields
    !
    implicit none
    !Arguments:
    type(struct_columnData) :: columnhr
    type(struct_obs)        :: obsSpaceData
    !Locals:
    integer  :: headerIndex, idatyp, nobs
    real(8)  :: height1, height2

    ! count number of tovs locations
    nobs = 0
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if ( .not.  tvs_isIdBurpTovs(idatyp) ) cycle HEADER 
      nobs = nobs + 1
    end do HEADER
    
    if ( nobs > 0 ) then
      allocate(trialHeight300m1000(nobs))
      allocate(trialHeight50m200(nobs))
      allocate(trialHeight5m50(nobs))
      allocate(trialHeight1m10(nobs))
      allocate(trialTG(nobs))
      allocate(RadiosondeWeight(nobs))
    else
      write(*,*) 'bias_getTrialPredictors: NO OBS found'
      return
    end if

    nobs = 0

    call obs_set_current_header_list(obsSpaceData,'TO')

    HEADER2: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER2
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if ( .not.  tvs_isIdBurpTovs(idatyp) ) cycle HEADER2 
      nobs = nobs + 1

      height1 = logInterpHeight(columnhr,headerIndex,1000.d0)
      height2 = logInterpHeight(columnhr,headerIndex,300.d0)
      
      trialHeight300m1000(nobs) = height2 - height1

      height1 = logInterpHeight(columnhr,headerIndex,200.d0)
      height2 = logInterpHeight(columnhr,headerIndex,50.d0)

      trialHeight50m200(nobs) = height2 - height1

      height1 = height2
      height2 = logInterpHeight(columnhr,headerIndex,5.d0)

      trialHeight5m50(nobs) = height2 - height1

      height1 = logInterpHeight(columnhr,headerIndex,10.d0)
      height2 = logInterpHeight(columnhr,headerIndex,1.d0)

      trialHeight1m10(nobs) = height2 - height1

      trialTG(nobs) = col_getElem(columnhr,1,headerIndex,'TG')

    end do HEADER2

    if ( trialTG(1) > 150.0d0) then
      write(*,*) 'bias_getTrialPredictors: converting TG from Kelvin to deg_C'
      trialTG(:) = trialTG(:) - MPC_K_C_DEGREE_OFFSET_R8
    end if

    trialHeight300m1000(:) = 0.1d0 * trialHeight300m1000(:) ! conversion factor
    trialHeight50m200(:) = 0.1d0 * trialHeight50m200(:)
    trialHeight5m50(:) = 0.1d0 * trialHeight5m50(:)
    trialHeight1m10(:) =  0.1d0 *  trialHeight1m10(:)

    write(*,*) 'bias_getTrialPredictors: end'

  contains

    function logInterpHeight(columnhr,headerIndex,P) result(height)
      !Arguments:
      type(struct_columnData), intent(inout) :: columnhr
      integer, intent(in) :: headerIndex
      real(8), intent(in) :: P
      real(8) :: height
      !Locals
      integer :: jk, nlev, ik
      real(8) :: zpt, zpb, zwt, zwb
      real(8), pointer :: col_ptr(:)

      ik = 1
      nlev = col_getNumLev(COLUMNHR,'TH')
      do jk = 2,NLEV - 1
        ZPT = col_getPressure(COLUMNHR,jk,headerIndex,'TH')* MPC_MBAR_PER_PA_R8
        if( P > ZPT ) ik = jk
      end do
      ZPT = col_getPressure(COLUMNHR,ik,headerIndex,'TH') * MPC_MBAR_PER_PA_R8
      ZPB = col_getPressure(COLUMNHR,ik+1,headerIndex,'TH') * MPC_MBAR_PER_PA_R8

      zwb  = log(p/zpt) / log(zpb/zpt)
      zwt  = 1.d0 - zwb
      col_ptr=>col_getColumn(columnhr,headerIndex,'Z_T')

      height = zwb * col_ptr(ik+1) + zwt * col_ptr(ik)
   
    end function logInterpHeight

  end subroutine bias_getTrialPredictors


  !----------------------
  ! bias_cvToCoeff
  !----------------------
  subroutine bias_cvToCoeff(cv_bias)
    !
    ! :Purpose: get coefficient increment from control vector
    !
    implicit none
    !Parameter:
    real(8)  :: cv_bias(:)
    !Locals:
    integer  :: index_cv,iSensor,iChannel,iPredictor,iScan
    integer  :: nsize,ierr
 
    if ( mpi_myid == 0 ) then
      write(*,*) 'bias_cvToCoeff: start'
      index_cv = 0
      ! initialize of coeffIncr
      do iSensor = 1, tvs_nSensors
        if (bias(iSensor)%numScan > 0) then
          do iChannel=1, bias(iSensor)%numChannels
            if ( bias(iSensor)%chans(iChannel)%isDynamic) then
              bias(iSensor)%chans(iChannel)%coeffIncr(:) = 0.0d0
              bias(iSensor)%chans(iChannel)%coeffIncr_fov(:) = 0.0d0
            end if
          end do
        end if
      end do

      do iSensor = 1, tvs_nSensors
        if (bias(iSensor)%numScan > 0) then
          do iChannel = 1, bias(iSensor)%numChannels
            if ( bias(iSensor)%chans(iChannel)%isDynamic ) then
              do iPredictor = 1,  bias(iSensor)%chans(iChannel)%numActivePredictors
                if ( iPredictor == 1 ) then
                  do iScan = 1, bias(iSensor)%numScan
                    index_cv = index_cv + 1
                    bias(iSensor)%chans(iChannel)%coeffIncr_fov(iScan) = bias(iSensor)%chans(iChannel)%stddev(iPredictor)*cv_bias(index_cv)
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
          if ( bias(iSensor)%chans(iChannel)%isDynamic ) &
               call rpn_comm_bcast(bias(iSensor)%chans(iChannel)%coeffIncr_fov,nsize,"mpi_double_precision",0,"GRID",ierr)
        end do
      end if
    end do

    ! for predictor part
    do iSensor = 1, tvs_nSensors
      if (bias(iSensor)%numScan > 0) then
        do iChannel = 1, bias(iSensor)%numChannels
          if ( bias(iSensor)%chans(iChannel)%isDynamic ) then
            nsize = bias(iSensor)%chans(iChannel)%numActivePredictors 
            call rpn_comm_bcast(bias(iSensor)%chans(iChannel)%coeffIncr,nsize,"mpi_double_precision",0,"GRID",ierr)
          end if
        end do
      end if
    end do

  end subroutine bias_cvToCoeff

  !-----------------------------------
  ! bias_getPredictors
  !---------------------------------- 
   subroutine bias_getPredictors(predictor,headerIndex,obsIndex,chanindx,obsSpaceData)
     !
     ! :Purpose: get predictors
     !
    implicit none
    !Arguments:
    real(8), intent(out)            :: predictor(NumPredictors)
    integer, intent(in)             :: headerIndex, obsIndex, chanindx
    type(struct_obs), intent(inout) :: obsSpaceData
    !Locals:
    integer  :: iSensor, iPredictor, jPredictor
    real(8)  :: zenithAngle

    predictor(:) = 0.0d0
    
    do iPredictor = 1, NumPredictors

      if ( iPredictor == 1 ) then
        ! constant
        predictor(iPredictor) = 1.0d0
      else if ( iPredictor == 2 ) then
        ! Height300-Height1000 (dam) /1000 T1
        predictor(iPredictor) = trialHeight300m1000(obsIndex)/1000.0d0
      else if ( iPredictor == 3 ) then
        ! Height50-Height200 (dam) /1000   T2
        predictor(iPredictor) = trialHeight50m200(obsIndex)/1000.0d0
      else if ( iPredictor == 4 ) then
        ! Height5-Height50 (dam) /1000    T3
        predictor(iPredictor) = trialHeight5m50(obsIndex)/1000.0d0
      else if ( iPredictor == 5 ) then
        ! Height1-Height10 (dam) /1000    T4
        predictor(iPredictor) = trialHeight1m10(obsIndex)/1000.0d0
        !      else if ( iPredictor == 6 ) then
        !        ! skin temperature (C) /10
        !        predictor(iPredictor) = trialTG(obsIndex)
      else if ( iPredictor == 6 ) then
        ! SV secant of satellite zenith angle minus one
        zenithAngle = obs_headElem_r(obsSpaceData,OBS_SZA,headerIndex) 
        if (zenithAngle < 75. ) predictor(iPredictor) = (1.d0 /cos( zenithAngle * MPC_RADIANS_PER_DEGREE_R8 ) ) - 1.d0
        !if (headerIndex == 4731) then
        ! print *,"ZXC", obs_headElem_r(obsSpaceData,OBS_SZA,headerIndex), predictor(iPredictor)
        ! call flush(6)
        !end if
      end if

    end do

    iSensor = tvs_lsensor(tvs_tovsIndex(headerIndex))
    do  iPredictor = 1,  bias(iSensor)%chans(chanIndx)%numActivePredictors
      jPredictor = bias(iSensor)%chans(chanIndx)%predictorIndex(iPredictor)
      predictor(jPredictor) =  predictor(jPredictor) - bias(iSensor)%chans(chanindx)%coeff_offset(iPredictor)
    end do
   
  end subroutine bias_getPredictors


  !---------------------------------------------
  ! bias_calcBias_ad
  !---------------------------------------------
  subroutine bias_calcBias_ad(cv_out,obsColumnIndex,obsSpaceData)
    !
    ! :Purpose: bias computation adjoint (for varBC)
    !
    implicit none
    !Arguments:
    real(8), intent(in)  :: cv_out(:)
    integer, intent(in)  :: obsColumnIndex
    type(struct_obs)     :: obsSpaceData
    !Locals:
    integer  :: headerIndex, bodyIndex, iobs, idatyp
    integer  :: iSensor, iChannel, iPredictor, index_cv, nsize, ierr, chanIndx
    integer  :: iScan, iFOV, jPred
    real(8)  :: predictor(NumPredictors)
    real(8), pointer  :: cv_bias(:)
    real(8), target  :: dummy4Pointer(1)
    real(8)  :: biasCor

    if ( .not. biasActive ) return

    if ( mpi_myid == 0 ) write(*,*) 'bias_calcBias_ad: start'

    nullify(cv_bias)
    if ( mpi_myid == 0 ) then
      if ( cvm_subVectorExists('BIAS') ) then
        cv_bias => cvm_getSubVector(cv_out,'BIAS')
      else
        write(*,*) 'bias_calcBias_ad: control vector does not include bias coefficients'
        return
      end if
    else
      cv_bias => dummy4Pointer
    end if

    ! adjoint of applying bias increment to specified obs column
    do iSensor = 1, tvs_nSensors
      if (bias(iSensor)%numScan > 0) then
        do iChannel = 1, bias(iSensor)%numChannels
          if ( bias(iSensor)%chans(iChannel)%isDynamic) then
            bias(iSensor)%chans(iChannel)%coeffIncr(:) = 0.0d0
            bias(iSensor)%chans(iChannel)%coeffIncr_fov(:) = 0.0d0
          end if
        end do
      end if
    end do

    iobs = 0
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if ( .not.  tvs_isIdBurpTovs(idatyp) ) cycle HEADER
      iobs = iobs + 1
      if ( tvs_tovsIndex(headerIndex) < 0 ) cycle HEADER

      iSensor = tvs_lsensor(tvs_tovsIndex(headerIndex))

      call obs_set_current_body_list(obsSpaceData, headerIndex)
      iFov = obs_headElem_i(obsSpaceData,OBS_FOV,headerIndex)

      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY

        call bias_getChannelIndex(obsSpaceData,iSensor,chanIndx,bodyIndex)

        if (chanindx > 0) then
          if (bias(iSensor)%chans(chanIndx)%isDynamic) then
            call bias_getPredictors(predictor,headerIndex,iobs,chanIndx,obsSpaceData)
            biasCor = obs_bodyElem_r(obsSpaceData,obsColumnIndex,bodyIndex)
            do iPredictor = 1, bias(iSensor)%chans(chanIndx)%numActivePredictors
              jPred = bias(iSensor)%chans(chanIndx)%PredictorIndex(iPredictor)
              if ( jPred == 1 ) then
                if ( bias(iSensor)%numScan > 1) then
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
    call bias_cvToCoeff_ad(cv_bias)

    if ( mpi_myid == 0 ) then
      write(*,*) 'bias_calcBias_ad: maxval(cv_bias)=',maxval(cv_bias(:))
    end if

  end subroutine bias_calcBias_ad

  !----------------------------------------------------
  ! bias_cvToCoeff_ad
  !----------------------------------------------------
  subroutine bias_cvToCoeff_ad(cv_bias)
    !
    ! :Purpose: adjoint of control vector to coeff transfer (for varBC)
    !
    implicit none
    !Parameter:
    real(8)  :: cv_bias(:)
    !Locals:
    integer  :: index_cv, iSensor, iChannel, iPredictor, iScan
    integer  :: nChan, nScan
    integer  :: nsize, ierr, iChan
    real(8), allocatable  :: temp_coeffIncr(:), temp_coeffIncr_fov(:)


    Write(*,*) "bias_cvToCoeff_ad: start"

    do iSensor = 1, tvs_nSensors
      if (bias(iSensor)%numScan > 0) then
        nChan = bias(iSensor)%numChannels
        do ichan = 1, nChan
          if ( bias(iSensor) % chans(ichan) %isDynamic ) then
            nSize  = bias(iSensor)%chans(iChan)%numActivePredictors
            allocate(temp_coeffIncr(nSize))
            temp_coeffIncr(:) = 0.0d0
            call rpn_comm_reduce(bias(iSensor)%chans(ichan)%coeffIncr(:),temp_coeffIncr,nsize,"mpi_double_precision","mpi_sum",0,"GRID",ierr)
            bias(iSensor)%chans(ichan)%coeffIncr(:) = temp_coeffIncr(:)
            deallocate(temp_coeffIncr)
          end if
        end do
      end if
    end do

    do iSensor = 1, tvs_nSensors
      if (bias(iSensor)%numScan > 0) then
        nChan = bias(iSensor)%numChannels
        nScan  = bias(iSensor)%numScan
        nsize = nChan * nScan
        if (nsize > 0) then
          allocate(temp_coeffIncr_fov(1:nScan))
          do ichan = 1, nChan
            if ( bias(iSensor) % chans(ichan) %isDynamic) then
              temp_coeffIncr_fov(:) = 0.0d0
              call rpn_comm_reduce(bias(iSensor)%chans(ichan)%coeffIncr_fov,temp_coeffIncr_fov,nscan,"mpi_double_precision","mpi_sum",0,"GRID",ierr)
              bias(iSensor)%chans(iChan)%coeffIncr_fov(:) = temp_coeffIncr_fov(:)
            end if
          end do
          deallocate(temp_coeffIncr_fov)
        end if
      end if
    end do

    if ( mpi_myid == 0 ) then
      cv_bias(:) = 0.d0
      index_cv = 0
      do iSensor = 1, tvs_nSensors
        if (bias(iSensor)%numScan > 0) then
          do iChannel = 1, bias(iSensor)%numChannels
            if (bias(iSensor) % chans(iChannel) %isDynamic) then
              do iPredictor = 1, bias(iSensor)%chans(iChannel)%numActivePredictors
                if ( iPredictor == 1 ) then
                  do iScan = 1, bias(iSensor)%numScan
                    index_cv  =  index_cv  +  1
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
    
  end subroutine bias_cvToCoeff_ad

  !-----------------------------------------
  ! bias_writeBias
  !-----------------------------------------
  subroutine bias_writeBias(cv_in)
    !
    ! :Purpose: to write bias increments and coefficients (varBC)
    !
    implicit none
    !Parameter:
    real(8), optional  :: cv_in(:)
    !Locals:
    integer  :: iSensor, iChannel, iPredictor, iScan, instIndex
    integer  :: jSensor, jChannel, iChannel2, myNbScan
    integer  :: nulfile_inc, nulfile_fov, ierr
    integer, external  :: fnom, fclos
    real(8), pointer :: cv_bias(:)
    real(8), target  :: dummy4Pointer(1)
    character(len=80) :: BgFileName
    !for background coeff and write out
    integer             :: iInstr
    integer             :: numCoefFile,jCoef,kCoef
    character(len=10)   :: coefInstrName(tvs_nSensors), temp_instrName
    character(len=25)   :: filecoeff
    logical             :: coeffExists
    ! these variables are not used but need to be present to satisfy bias_updateCoeff interface
    ! some bias_updateCoeff arguments could be made optional (todo)
    integer            :: chans(tvs_nSensors, tvs_maxChannelNumber), nsat, nfov
    integer            :: nchan(tvs_nSensors)
    character(len=10)  :: sats(tvs_nsensors) ! satellite names
    character(len=7)   :: cinstrum           ! string: instrument (e.g. AMSUB)

    if ( .not. biasActive ) return

    if ( present(cv_in) ) then
      nullify(cv_bias)
      if ( mpi_myid == 0 ) then
        if ( cvm_subVectorExists('BIAS') ) then
          cv_bias => cvm_getSubVector(cv_in,'BIAS')
          write(*,*) 'bias_writeBias: maxval(cv_bias)=',maxval(cv_bias(:))
        else
          write(*,*) 'bias_writeBias: control vector does not include bias coefficients'
          return
        end if
      else
        cv_bias => dummy4Pointer
      end if
      call bias_cvToCoeff(cv_bias)
    end if

    ! apply transformation to account for predictor offset

    do iSensor = 1, tvs_nSensors
      if (bias(iSensor)%numScan > 0) then
        do iChannel = 1, bias(iSensor)%numChannels
          do iPredictor = 2, bias(iSensor)%chans(iChannel)%numActivePredictors
            if (lMimicSatbcor) then
              bias(iSensor)%chans(iChannel)%coeffIncr(1) = bias(iSensor)%chans(iChannel)%coeffIncr(1) - &
                   bias(iSensor)%chans(iChannel)%coeffIncr(iPredictor) *  bias(iSensor)%chans(iChannel)%coeff_offset(iPredictor)
            else
              do iScan = 1, bias(iSensor)%numScan
                bias(iSensor)%chans(iChannel)%coeffIncr_fov(iScan) = bias(iSensor)%chans(iChannel)%coeffIncr_fov(iScan) - &
                     bias(iSensor)%chans(iChannel)%coeffIncr(iPredictor) *  bias(iSensor)%chans(iChannel)%coeff_offset(iPredictor)
              end do
            end if
          end do
          if (bias(iSensor)%chans(iChannel)%numActivePredictors > 0 .and. mpi_myId ==0 .and. (.not. doRegression) ) &
               Write(*,*) "bias_writeBias: bias(iSensor)%chans(iChannel)%coeffIncr(:) =",  bias(iSensor)%chans(iChannel)%coeffIncr(:)
        end do
      end if
    end do

    if (  doRegression ) then
      call bias_writeCoeff()
      return
    end if

    if (mpi_myId==0) then

      ! write out bias coefficient increments in ascii file
      nulfile_inc = 0
      ierr = fnom(nulfile_inc,'./satbias_increment.dat','FTN+FMT',0)

      do iSensor = 1, tvs_nSensors
        write(nulfile_inc,'(/,1X,"Sensor Index=",I3,", Satellite Name=",A15,", Instrument Name=",A15)') &
             iSensor,tvs_satelliteName(iSensor),tvs_instrumentName(iSensor)
        if (bias(iSensor)%numScan > 0) then
          do iChannel = 1, bias(iSensor)%numChannels
            if (bias(iSensor)%chans(iChannel)%isDynamic) then
              iChannel2 = bias(iSensor)%chans(iChannel)%channelNum
              if ( sum(bias(iSensor)%chans(iChannel)%coeffIncr(:)) /= 0.0d0 ) &
                   write(nulfile_inc,'(3X,"Channel number=",I4)') iChannel2
              do iPredictor = 2, bias(iSensor)%chans(iChannel)%numActivePredictors
                if ( bias(iSensor)%chans(iChannel)%coeffIncr(iPredictor) /= 0.0d0 ) &
                     write(nulfile_inc,'(5X,"Predictor number=",I4,", Coefficient=",e12.4)') &
                     iPredictor,bias(iSensor)%chans(iChannel)%coeffIncr(iPredictor)
              end do
            end if
          end do
        end if
      end do

      ierr = fclos(nulfile_inc)

      ! write out fovbias coefficient increments in ascii file
      nulfile_fov = 0
      ierr = fnom(nulfile_fov,'./fovbias_incre.dat','FTN+FMT',0)
      do iSensor = 1, tvs_nSensors
        write(nulfile_fov,'(/,1X,"Sensor Index=",I3,", Satellite Name=",A15,", Instrument Name=",A15)') &
             iSensor,tvs_satelliteName(iSensor),tvs_instrumentName(iSensor)
        if (bias(iSensor)%numScan > 0) then
          do iChannel = 1, bias(iSensor)%numChannels
            if (bias(iSensor)%chans(iChannel)%isDynamic) then
              iChannel2 =  bias(iSensor)%chans(iChannel)%channelNum
              if ( sum(bias(iSensor)%chans(iChannel)%coeffIncr_fov(:)) /= 0.0d0 ) then
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
        filecoeff = 'coeff_file_'//trim(temp_instrName)//''
        inquire(file=trim(filecoeff),exist = coeffExists)

        if ( coeffExists ) then
          numCoefFile = numCoefFile + 1
          jCoef = jCoef + 1
          coefInstrName(jCoef) = temp_instrName
        end if
        if ( jSensor > 1 ) then
          do kCoef = 1, jCoef-1
            if ( temp_instrName == coefInstrName(kCoef) ) then
              numCoefFile = numCoefFile - 1
              jCoef = jCoef -1
            end if
          end do
        end if 
      end do
    end do

    ! update coeff_file_instrument and write out
    do iInstr=1, numCoefFile 
      BgFileName = './coeff_file_' // coefInstrName(iInstr)
      call bias_updateCoeff(tvs_nSensors,NumPredictors,BgFileName,sats,chans,nsat,nchan,nfov,cinstrum)
    end do


  end subroutine bias_writeBias

  
  !-----------------------------------------
  ! bias_updateCoeff
  !-----------------------------------------
  subroutine bias_updateCoeff(maxsat,maxpred,coeff_file,sats,chans,nsat,nchan,nfov,cinstrum,updateCoeff_opt)
    !
    ! :Purpose: to read, and optionaly update and write out, the coeff files (varBC).
    !
    implicit none
    !Arguments:
    integer,intent(in)            :: maxsat, maxpred
    character(len=*), intent(in)  :: coeff_file
    logical, optional, intent(in) :: updateCoeff_opt
    integer, intent(out)          :: chans(maxsat, tvs_maxChannelNumber)       ! channel numbers
    integer, intent(out)          :: nsat, nfov
    integer, intent(out)          :: nchan(maxsat)  ! number of channels
    character(len=10), intent(out):: sats(maxsat)   ! dim(maxsat), satellite names
    character(len=*), intent(out) :: cinstrum       ! instrument (e.g. AMSUB)
    !Locals:
    real(8)            :: fovbias(maxsat,tvs_maxChannelNumber,maxfov)     ! dim(maxsat,tvs_maxchannelnumber,maxfov), bias as F(fov)
    real(8)            :: coeff(maxsat,tvs_maxChannelNumber,maxpred)       ! dim(maxsat,tvs_maxchannelnumber,maxpred)
    character(len=2)   :: ptypes(maxsat,tvs_maxChannelNumber,maxpred) ! dim(maxsat,tvs_maxChannelNumbermaxchan,maxpred)
    integer            :: npred(maxsat, tvs_maxChannelNumber)       ! dim(maxsat, tvs_maxchannelnumber), number of predictors
    integer            :: ndata(maxsat, tvs_maxChannelNumber)
    ! LOCAL for reading background coeff file
    character(len=10)  :: sat
    character(len=120) :: line
    integer            :: chan, nbfov, nbpred
    integer            :: mincol, maxcol, nbcol
    integer            :: ier, istat 
    integer            :: iSat, jChan, kPred, kFov, totSat
    logical            :: newsat, verbose
    real               :: dummy
    integer            :: iun
    ! update coeff files
    real               :: fovbias_an(maxsat,tvs_maxChannelNumber,maxfov)
    real               :: coeff_an(maxsat,tvs_maxChannelNumber,maxpred) 
    integer            :: iSensor, jChannel , iFov, iPred, totPred
    character(len=10)  :: tmp_SatName, tmp_InstName 
    ! write out files 
    integer            :: iuncoef2, ierr, numPred
    character(len=80)  :: filename2
    logical            :: updateCoeff_opt2
    !   sats(nsat)            = satellite names
    !   chans(nsat,nchan(i))  = channel numbers of each channel of each satellite i
    !   npred(nsat,nchan(i))  = number of predictors for each channel of each satellite i
    !   fovbias(i,j,k)        = bias for satellite i, channel j, FOV k   k=1,nfov
    !     if FOV not considered for instrument, nfov = 1 and fovbias is global bias for channel
    !   coeff(i,j,1)          = regression constant
    !   coeff(i,j,2), ..., coeff(i,j,npred(i,j)) = predictor coefficients
    !   nsat, nchan, nfov, cinstrum (output) are determined from file
    !   if returned nsat = 0, coeff_file was empty
    !   maxpred (input) is max number of predictors
    !   maxsat (input)  is max number of satellites
    ! There are three parts in this subroutine, read, update and write out the coeff files
    ! 
    !- 1. read in the background coeff files, this program is read_coeff from genbiascorr
    ! 
    if ( present(updateCoeff_opt) ) then
      updateCoeff_opt2 = updateCoeff_opt
    else
      updateCoeff_opt2 = .true.
    end if

    verbose = .false.
   
    call read_coeff(sats, chans, fovbias, coeff, nsat, nchan, nfov, &
         npred, cinstrum, coeff_file, ptypes,ndata)

    ! Transfer of coefficients read from file to the bias structure
    satloop :do iSat = 1, nsat  !backgroud sat
      instloop:do iSensor = 1, tvs_nSensors
        ! for Satellite Name
        tmp_SatName = SatNameinCoeffFile(tvs_satelliteName(iSensor))
        ! for Instrument Name
        tmp_InstName = InstrNameinCoeffFile(tvs_instrumentName(iSensor))
        
        if ( trim(tmp_SatName) /= trim(sats(iSat)) .or. trim(tmp_InstName) /= trim(cinstrum) ) cycle instloop
        write(*,*) "bias_updateCoeff: " //tmp_SatName//" "//tmp_InstName

        if ( .not. allocated(bias(iSensor)%chans ) ) cycle instloop

        chan1loop:do jChan = 1, nchan(iSat)
          chan2loop:do jChannel = 1, bias(iSensor)%numChannels  

            if ( chans(iSat,jChan) /= bias(iSensor)%chans(jChannel)%channelNum ) cycle chan2loop

            ! part 1 for coeffIncr
            do iFov = 1, nfov
              bias(iSensor)%chans(jchannel)%coeff_fov(iFov) = fovbias(iSat,jChan,iFov)
            end do ! iFov
            
            ! part 2 for coeffIncr_fov
            totPred  = bias(iSensor)%chans(jchannel)%NumActivePredictors 
            do iPred = 1, totPred
              bias(iSensor)%chans(jchannel)%coeff(iPred) = coeff(iSat,jChan,iPred)
            end do ! iPred
            
            
          end do chan2loop ! jChannel
        end do chan1loop !jChan
        
      end do instloop
    end do satloop
    
    if ( .not. updateCoeff_opt2 ) return 
    
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
        if ( trim(tmp_SatName) /= trim(sats(iSat)) .or. trim(tmp_InstName) /= trim(cinstrum) ) cycle 
        do jChan = 1, nchan(iSat)
          do jChannel = 1, bias(iSensor)%numChannels  

            if ( chans(iSat,jChan) /= bias(iSensor)%chans(jChannel)%channelNum ) cycle

            ! part 1 for coeffIncr
            do iFov = 1, nfov
              fovbias_an(iSat,jChan,iFov) = fovbias(iSat,jChan,iFov) + bias(iSensor)%chans(jchannel)%coeffIncr_fov(iFov)
            end do ! iFov

            ! part 2 for coeffIncr_fov
            totPred  = bias(iSensor)%chans(jchannel)%NumActivePredictors 
            do iPred = 1, totPred
              coeff_an(iSat,jChan,iPred) = coeff(iSat,jChan,iPred) + bias(iSensor)%chans(jchannel)%coeffIncr(iPred)
            end do ! iPred

          end do ! jChannel
        end do !jChan
      end do !iSensor
    end do ! iSat

    !
    !- 3. Write out updated_coeff
    ! 

    if (mpi_myId == 0 ) then

      iuncoef2 = 0
      filename2 ='./anlcoeffs_'//cinstrum 
      ierr = fnom(iuncoef2, filename2,'FTN+FMT',0)

      write(*,*) 'bias_updateCoeff: write in bias_updateCoeff'
   
      do iSat = 1, nsat
        do jChan = 1, nchan(iSat)      
          numPred = npred(iSat,jChan)
          if (sum(abs(coeff_an(iSat,jchan,1:numpred+1)))/=0.d0 .and. sum(abs(fovbias_an(iSat,jChan,1:nfov)))/=0.d0) then 
            write(iuncoef2, '(A52,A8,1X,A7,1X,I6,1X,I8,1X,I2,1X,I3)') &
                 'SATELLITE, INSTRUMENT, CHANNEL, NOBS, NPRED, NSCAN: ', sats(iSat), cinstrum, chans(iSat,jChan), ndata(isat,jchan), numPred, nfov
            write(iuncoef2, '(A7,6(1X,A2))') 'PTYPES:', (ptypes(iSat,jChan,kPred) , kPred=1,numPred)
            write(iuncoef2,'(120(1x,ES17.10))') (fovbias_an(iSat,jChan,kFov),kFov=1,nfov)
            write(iuncoef2,*) (coeff_an(iSat,jChan,kPred),kPred=1,numPred+1)
          end if
        end do
      end do

      ierr = fclos(iuncoef2) 
   
      write(*,*) 'bias_updateCoeff: finish writing coeffient file', filename2
    
    end if

  end subroutine bias_updateCoeff


  !-----------------------------------------
  ! bias_writeCoeff
  !-----------------------------------------
  subroutine bias_writeCoeff()
    !
    ! :Purpose:  write out  the coeff files (regression case).
    !
    implicit none
    !Locals:
    integer            :: iuncoef, ierr, numPred
    character(len=80)  :: filename
    character(len=80)  :: instrName, satNamecoeff
    character(len=2), parameter  :: predTab(2:6) = (/ "T1", "T2", "T3", "T4", "SV"/)
    integer :: sensorIndex, nchans, nscan, nfov, kpred, kFov, jChan

    if (mpi_myId == 0 ) then

      SENSORS:do sensorIndex = 1, tvs_nsensors

        if  (.not. tvs_isReallyPresentMpiGLobal(sensorIndex)) cycle SENSORS

        write(*,*) "bias_writeCoeff: sensorIndex ",  sensorIndex

        nchans = bias(sensorIndex)%numChannels
        nscan   = bias(sensorIndex)% numscan

        instrName = InstrNameinCoeffFile(tvs_instrumentName(sensorIndex))
        satNamecoeff = SatNameinCoeffFile(tvs_satelliteName(sensorIndex)) 

        iuncoef = 0
        filename = './anlcoeffs_' // trim(instrName)  
        call utl_open_asciifile(filename,iuncoef)
        nfov = bias(sensorIndex) % numScan
        do jChan = 1, nchans
          if ( bias(sensorIndex) % chans(jChan) %coeff_nobs > 0) then
            numPred = bias(sensorIndex) % chans(jChan) % numActivePredictors 
          
            write(iuncoef, '(A52,A8,1X,A7,1X,I6,1X,I8,1X,I2,1X,I3)') 'SATELLITE, INSTRUMENT, CHANNEL, NOBS, NPRED, NSCAN: ',  &
                 satNameCoeff, instrName, bias(sensorIndex) % chans(jChan) %channelNum, bias(sensorIndex) % chans(jChan) %coeff_nobs, numPred -1, nfov
            write(iuncoef, '(A7,6(1X,A2))') 'PTYPES:',  ( predtab(bias(sensorIndex) % chans(jChan) %predictorIndex(kPred)) , kPred=2, numPred )
            write(iuncoef,'(120(1x,ES17.10))') (bias(sensorIndex) % chans(jChan) %coeff_fov(kFov),kFov=1,nfov)
            write(iuncoef,'(12(1x,ES17.10))') (bias(sensorIndex) % chans(jChan) %coeff(kPred),kPred=1,numPred)
          end if
        end do

        close(iuncoef) 

        if (loutCoeffCov) then
          iuncoef = 0
          filename ='./anlcoeffsCov_'// trim(instrName)  
          call utl_open_asciifile(filename,iuncoef)
          do jChan = 1, nchans
            if ( bias(sensorIndex) % chans(jChan) %coeff_nobs > 0) then
              numPred = bias(sensorIndex) % chans(jChan) % numActivePredictors 
          
              write(iuncoef, '(A38,A8,1X,A7,1X,I6,1X,I2)') 'SATELLITE, INSTRUMENT, CHANNEL, NPRED: ',  &
                   satNameCoeff, instrName, bias(sensorIndex) % chans(jChan) %channelNum, numPred
              do kpred =1, numPred
                write(iuncoef, '(10e14.6)')  bias(sensorIndex)%chans(jChan)%coeffCov(kpred,:)
              end do
            end if
          end do
          close(iuncoef) 
        end if

      end do SENSORS

    end if

 end subroutine bias_writeCoeff

 !-----------------------------------------
 ! bias_removeBiasCorrection
 !-----------------------------------------
 subroutine bias_removeBiasCorrection(obsSpaceData,family_opt)
    !
    ! :Purpose: to remove bias correction from OBS_VAR
    !           After the call OBS_VAR contains the uncorrected
    !           observation and OBS_BCOR is set to zero
    !
    implicit none
    !Arguments:
    type(struct_obs)                       :: obsSpaceData
    character(len=2), intent(in), optional :: family_opt
    !Locals:
    integer :: nbcor
    integer :: bodyIndex
    real(8) :: biascor, Obs

    if (.not. removeBiasCorrection) return

    if ( mpi_myid == 0 ) write(*,*) 'bias_removeBiasCorrection: start'

    nbcor = 0
    call obs_set_current_body_list(obsSpaceData,family_opt)
    
    BODY: do
      bodyIndex = obs_getBodyIndex(obsSpaceData)
      if ( bodyIndex < 0 ) exit BODY
      if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY  
      biasCor = obs_bodyElem_r(obsSpaceData,OBS_BCOR,bodyIndex)
      Obs =  obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
      if (biasCor /= MPC_missingValue_R8 .and. Obs /= MPC_missingValue_R8) then
        call obs_bodySet_r(obsSpaceData, OBS_VAR, bodyIndex, real(Obs - biasCor,OBS_REAL))
        call obs_bodySet_r(obsSpaceData, OBS_BCOR, bodyIndex, real(0.d0,OBS_REAL))
        nbcor = nbcor + 1
      end if
    end do BODY

    if ( mpi_myid == 0 ) then
      write(*,*) 'bias_removeBiasCorrection: removed bias correction for ', nbcor, ' observations'
      write(*,*) 'bias_removeBiasCorrection exiting'
    end if

  end subroutine bias_removeBiasCorrection

  !-----------------------------------------
  ! bias_filterObs
  !-----------------------------------------
  subroutine bias_filterObs(obsSpaceData)
    !
    ! :Purpose: to filter radiance observations to include into
    !           bias correction offline computation
    !           (same rules as in bgck.gen_table)
    !
    implicit none
    !Parameter:
    type(struct_obs), intent(inout) :: obsSpaceData
    !Locals:
    integer :: bodyIndex, headerIndex
    integer :: assim,flag,codtyp
    integer :: isatBufr, instBufr, iplatform, isat, inst ,idsat,i,chanIndx
    logical :: lHyperIr, lGeo, lSsmis, lTovs
    logical :: condition, condition1, condition2

    if (.not. filterObs ) return

    if ( mpi_myid == 0 ) write(*,*) 'bias_filterObs: start'

    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER

      codtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)

      lHyperIr = .false.
      lGeo     = .false.
      lSsmis   = .false.
      lTovs    = .false.

      select case ( codtyp_get_name(codtyp) )
      case("ssmis")
        lSsmis = .true.
      case("csr")
        lGeo = .true.
      case("airs","iasi","cris","crisfsr")
        lHyperIr = .true.
      case default
        lTovs = .true.
      end select

      isatBufr = obs_headElem_i(obsSpaceData,OBS_SAT,headerIndex) !BUFR element 1007
      instBufr = obs_headElem_i(obsSpaceData,OBS_INS,headerIndex) !BUFR element 2019

      call tvs_mapSat(isatBufr,iplatform,isat)
      call tvs_mapInstrum(instBufr,inst)

      idsat = -1
      do i =1, tvs_nsensors
         if ( tvs_platforms(i)  == iplatform .and.  &
              tvs_satellites(i) == isat      .and.  &
              tvs_instruments(i)== inst )    then
          idsat = i
          exit
        end if
      end do
      if (idsat == -1) cycle HEADER

      call obs_set_current_body_list(obsSpaceData, headerIndex)

      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        assim =  obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex)

        if ( assim == obs_notAssimilated ) then
          call bias_getChannelIndex(obsSpaceData, idsat, chanIndx,bodyIndex)
          if (chanIndx > 0) then
            flag = obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex)

            if ( lSsmis ) then   ! SSM/I and SSMIS
              !   Here "good" means data that are not rejected by first QC program satqc_ssmi(s) (bit 7 ON).
              !   Bit 11 is ON for data that are unselected by UTIL or for uncorrected data (or both).
              !   Data rejected by first QC program satqc_ssmi(s) have bit 7 switched ON only (in addition to bit 9) as
              !   rogue/topo checks are skipped. So if bit 16 (rogue) is ON, bit 7 must be off.
              if ( .not. offline_opt ) then
                if ( all_opt_ssmis ) then
                  !  FLAG test: all good data (corrected/selected or not) that have passed all QC (bit 9 OFF)
                  condition1  = .not. btest(flag,9) !' AND (FLAG & 512 = 0)'
                  !  FLAG test: uncorrected good data that failed rogue check only ([bit 9 ON] + bit 6 OFF + bit 16 ON + bit 18 OFF + [bit 7 OFF])
                  condition2  = .not. btest(flag,6) .and. btest(flag,16) .and. .not. btest(flag,18) !' AND (FLAG & 64 = 0) AND (FLAG &  65536 = 65536) AND (FLAG & 262144 = 0)'
                  condition = condition1 .or. condition2
                else
                  !  FLAG test: corrected/selected good data that have passed QC (bits 9,11 OFF) --> data to be assimilated
                  condition =   .not. btest(flag,9) .and. .not. btest(flag,11)       !' AND (FLAG & 512 = 0) AND (FLAG & 2048 = 0)'
                end if
              else   ! OFFLINE MODE --> want all observations except data rejected for any reason other than rogue innovation check
                condition1  = .not. btest(flag,9) !' AND (FLAG & 512 = 0)'  
                ! all good data that passed all QC    
                ! "good" data that failed rogue check [bit 9 ON, bit 7 OFF, bit 18 OFF]
                condition2  = btest(flag,9) .and. .not. btest(flag,7) .and. .not. btest(flag,18) !' AND (FLAG & 512 = 512) AND (FLAG & 128 = 0) AND (FLAG & 262144 = 0)'
                condition = condition1 .or. condition2
              end if
            else if( lTovs ) then
              ! AMSU-A, AMSU-B/MHS, ATMS, MWHS-2
              !  In AMSU case, bit 11 is set for data that are not bias corrected or for unselected channels.
              !    BUT unlike other instruments, all AMSU data are bias corrected, whether selected or not
              !    so bit 11 = unselected channel (like bit 8 for AIRS/IASI)
              !  Bit 9 is set for all other rejections including rogue (9+16) and topography (9+18).
              !  In addition, bit 7 is set for channels with bad data or data that should not be assimilated.
              if ( .not. offline_opt ) then
                if ( all_opt_tovs ) then
                  !  FLAG test: all data (selected or not) that have passed QC (bit 9 OFF)
                  condition1  = .not. btest(flag,9) !' AND (FLAG & 512 = 0)'
                  !  FLAG test: uncorrected (bit 6 OFF) data that failed rogue check only (bit (9)/16 ON, 18,7 OFF)
                  !             NOTE: As all AMSU data are normally bias corrected, query2 will return nothing
                  condition2  =  btest(flag,16) .and. .not. btest(flag,6) .and. .not. btest(flag,18) .and. .not. btest(flag,7)!' AND (FLAG & 64 = 0) AND (FLAG &  65536 = 65536) AND (FLAG & 262144 = 0) AND (FLAG & 128 = 0)'
                  condition = condition1 .or. condition2
                else
                  !  FLAG test: selected data (bit 11 OFF) that have passed QC (bit 9 OFF)
                  condition  =  .not. btest(flag,9) .and.  .not. btest(flag,11) !' AND (FLAG & 512 = 0) AND (FLAG & 2048 = 0)'
                end if
              else    ! OFFLINE MODE --> want all observations except data rejected for any reason other than rogue check
                condition1  = .not. btest(flag,9) !' AND (FLAG & 512 = 0)'  
                ! all good data that passed all QC    
                ! "good" data that failed rogue check [bit 9 ON, bit 7 OFF, bit 18 OFF]
                condition2  =  btest(flag,9) .and. .not. btest(flag,7) .and. .not. btest(flag,18)  !' AND (FLAG & 512 = 512) AND (FLAG & 128 = 0) AND (FLAG & 262144 = 0)'
                condition = condition1 .or. condition2
              end if
            else if( lGeo ) then   !  AIRS, IASI, CSR, CRIS  
              ! CSR case!    No flag check        =                all data that have passed QC/filtering
              !  (FLAG & 2048 = 0)      = bit 11 OFF --> corrected/selected data that have passed QC/filtering
              if ( all_opt_csr .or. offline_opt ) then
                condition  = .true. !
              else        
                condition  = .not. btest(flag,18) ! ' AND (FLAG & 2048 = 0)' 
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
              !                                          data are assimilated]!  (FLAG & 128 = 0 )      = bit 7 OFF  --> not shortwave channel during day
              !  (FLAG & 512 = 0 )      = bit 9 OFF  --> non-erroneous data that passed O-P rogue check
              !! AIRS, IASI:!    bit  8 ON: blacklisted/unselected channel (UTIL=0)
              !    bit  9 ON: erroneous/suspect data (9), data failed O-P check (9+16)
              !    bit 11 ON: cloud (11+23), surface (11+19), model top transmittance (11+21), shortwave channel+daytime (11+7)
              !               not bias corrected (11) (with bit 6 OFF)
              if ( .not. offline_opt ) then
                if ( all_opt_hyperIr ) then        
                  ! good data that have passed all QC (bits 9 and 7,19,21,23 OFF), corrected/selected or not
                  condition1  = .not. btest(flag,9) .and. .not. btest(flag,7) .and. .not. btest(flag,19) .and. .not. btest(flag,21) .and. .not. btest(flag,23) !' AND (FLAG & 512 = 0) AND (FLAG & 11010176 = 0)'
                  ! uncorrected (6 OFF, [11 ON]) good data (7,19,21,23 OFF) that failed QC rogue check only (bits [9],16 ON), selected or not
                  condition2  = .not. btest(flag,6) .and. btest(flag,11) .and.  .not. btest(flag,17) .and. .not. btest(flag,19) .and. .not. btest(flag,21) .and. .not. btest(flag,23) 
                  !' AND (FLAG & 64 = 0) AND (FLAG & 65536 = 65536) AND (FLAG & 11010176 = 0)'
                  condition = condition1 .or. condition2
                else 
                  ! corrected data that passed all QC and selection excluding cloud/sfc affected obs
                  condition =  .not. btest(flag,9) .and. .not. btest(flag,11) .and.  .not. btest(flag,8) .and. .not. btest(flag,23) .and. .not. btest(flag,19) 
                  !' AND (FLAG & 2560 = 0) AND (FLAG & 256 = 0) AND (FLAG & 8388608 = 0) AND (FLAG & 524288 = 0)'
                end if
              else! OFFLINE MODE --> Want all observations except data rejected for any reason other than innovation rogue check
                !   Assumes that type S or N correction has been applied to all data/channels (all data "corrected")
                ! data that passed all QC 
                condition1 =  .not. btest(flag,9) .and. .not. btest(flag,7) .and. .not. btest(flag,19) .and. .not. btest(flag,21) .and. .not. btest(flag,23) 
                !' AND (FLAG & 512 = 0) AND (FLAG & 11010176 = 0)'
                ! good data (7,19,21,23 OFF) that failed QC rogue check only (bits [9],16 ON)
                condition2 = btest(flag,9) .and. btest(flag,16) .and. .not. btest(flag,7) .and. .not. btest(flag,19) .and. .not. btest(flag,21) .and. .not. btest(flag,23) !' AND (FLAG & 65536 = 65536) AND (FLAG & 11010176 = 0)'
                condition = condition1 .or. condition2
              end if
            end if

            if (condition) assim = obs_Assimilated

            call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, assim)
          end if
        end if
      end do BODY
    end do HEADER
 
  end subroutine bias_filterObs

  !-----------------------------------------
  ! bias_applyBiasCorrection
  !-----------------------------------------
  subroutine bias_applyBiasCorrection(obsSpaceData,column,family_opt)
    !
    ! :Purpose: to apply bias correction from OBS_BCOR to
    !           obsSpaceData column column.
    !           After the call obsSpaceData body column contains the corrected
    !           observation or O-F and OBS_BCOR is not modified.
    implicit none
    !Arguments:
    type(struct_obs)                       :: obsSpaceData
    integer, intent(in)                    :: column !obsSpaceData column
    character(len=2), intent(in), optional :: family_opt
    !Locals:
    integer :: nbcor
    integer :: bodyIndex
    real(8) :: biascor, Obs
    integer :: flag

    if ( .not. biasActive ) return
    if ( trim(biasMode) /= "apply") return

    if ( mpi_myid == 0 ) write(*,*) 'bias_applyBiasCorrection: start'

    nbcor = 0
    call obs_set_current_body_list(obsSpaceData,family_opt)
    
    BODY: do
      bodyIndex = obs_getBodyIndex(obsSpaceData)
      if ( bodyIndex < 0 ) exit BODY
      if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY  
      biasCor = obs_bodyElem_r(obsSpaceData,OBS_BCOR,bodyIndex)
      if (biasCor /= MPC_missingValue_R8) then
        Obs =  obs_bodyElem_r(obsSpaceData,column,bodyIndex)
        if (Obs /= MPC_missingValue_R8) then
          call obs_bodySet_r(obsSpaceData, column, bodyIndex, real(Obs + biasCor,OBS_REAL))
          flag = obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex)
          flag = ibset(flag, 6)
          call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, flag)
          nbcor = nbcor + 1
        end if
      end if
    end do BODY

    if ( mpi_myid == 0 ) then
      write(*,*) 'bias_applyBiasCorrection: apply bias correction for ', nbcor, ' observations'
      write(*,*) 'bias_applyBiasCorrection exiting'
    end if

  end subroutine bias_applyBiasCorrection


  !-----------------------------------------
  ! bias_refreshBiasCorrection
  !-----------------------------------------
  subroutine bias_refreshBiasCorrection(obsSpaceData,columnhr)
    !
    ! :Purpose: to apply bias correction from read coefficient file to OBS_VAR
    !           After the call OBS_VAR contains the corrected observation
    !           and OBS_BCOR is set to applied bias correction
    !
    implicit none
    !Arguments:
    type(struct_obs)        :: obsSpaceData
    type(struct_columnData) :: columnhr

    if ( .not.biasActive ) return
    if ( .not. refreshBiasCorrection) return

    if ( mpi_myid == 0 ) write(*,*) 'bias_refreshBiasCorrection: start'
    call bias_calcBias(obsSpaceData,columnhr)
    call bias_applyBiasCorrection(obsSpaceData,OBS_VAR,"TO")
    if ( mpi_myid == 0 ) write(*,*) 'bias_refreshBiasCorrection: exit'

  end subroutine bias_refreshBiasCorrection


  !-----------------------------------------
  ! bias_applyBiasCorrection
  !-----------------------------------------
  subroutine bias_getRadiosondeWeight(obsSpaceData,lmodify_obserror_opt)
    !
    ! :Purpose: initialize the weights to give more importance to data near radiosonde stations
    !
    implicit none
    !Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData
    logical, intent(in), optional   :: lmodify_obserror_opt
    !Locals:
    integer :: iobs, headerIndex, idatyp, nobs, bodyIndex
    logical :: lmodify_obserror
    real(8) :: sigmaObs

    lmodify_obserror =.false.

    if ( present(lmodify_obserror_opt) ) lmodify_obserror = lmodify_obserror_opt


    if (lweightedEstimate) then
      call hco_SetupFromFile(hco_mask, './raob_masque.std', 'WEIGHT', GridName_opt='RadiosondeWeight',varName_opt='WT' )
      call vco_SetupFromFile(vco_mask, './raob_masque.std')   ! IN

      call gsv_allocate(statevector_mask, 1, hco_mask, vco_mask, dateStampList_opt=(/-1/), varNames_opt=(/"WT"/), &
           dataKind_opt=4, mpi_local_opt=.true.,mpi_distribution_opt="Tiles")
     
      call gsv_readFromFile(statevector_mask,'./raob_masque.std' , 'WEIGHT', 'O', unitConversion_opt=.false., &
           containsFullField_opt=.false.)

      call col_setVco(column_mask,vco_mask)
      nobs = obs_numHeader(obsSpaceData)
      call col_allocate(column_mask, nobs,beSilent_opt=.false., varNames_opt=(/"WT"/) )
     
      call s2c_nl( stateVector_mask, obsSpaceData, column_mask, 'NEAREST', varName_opt="WT", moveObsAtPole_opt=.true.)

      call obs_set_current_header_list(obsSpaceData,'TO')
      iobs = 0
      HEADER: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if ( headerIndex < 0 ) exit HEADER
          
        ! process only radiance data to be assimilated?
        idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
        if ( .not. tvs_isIdBurpTovs(idatyp) ) cycle HEADER
      
        iobs = iobs + 1

        RadiosondeWeight(iobs) = col_getElem(column_mask,1,headerIndex) 

        if (lmodify_obserror) then
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if ( bodyIndex < 0 ) exit BODY
            
            if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY

            sigmaObs = obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)

            sigmaObs = sigmaObs / sqrt( RadiosondeWeight(iobs) )
            call obs_bodySet_r( obsSpaceData, OBS_OER, bodyIndex, sigmaObs) 

          end do BODY

        end if

      end do HEADER

    else
      RadiosondeWeight(:) = 1.d0
    end if

  end subroutine bias_getRadiosondeWeight


  !-----------------------------------------
  ! bias_do_regression
  !-----------------------------------------
  subroutine bias_do_regression(columnhr,obsSpaceData)
    !
    ! :Purpose: compute the bias correction coefficients by linear regresion
    !
    implicit none
    !Arguments:
    type(struct_obs), intent(inout)        :: obsSpaceData
    type(struct_columnData), intent(inout) :: columnhr
    !Locals:
    integer    :: iSensor, iChannel, nobs, ifail, npred,nchans, nscan, ndim, ndimmax
    integer    :: sensorIndex,iPred1,jPred1,iobs
    integer    :: headerIndex, idatyp,nPredMax,ierr,iFov, iScan, idim
    integer    :: indxtovs, bodyIndex, chanIndx, predstart, ntot
    real(8)    :: OmF, sigmaObs, lambda,norm
    real(8)    :: predictor(NumPredictors)
    real(8), allocatable :: Matrix(:,:,:), Vector(:,:)
    real(8), allocatable :: matrixMpiGlobal(:,:,:), vectorMpiGlobal(:,:)
    real(8), allocatable :: pIMatrix(:,:), OmFBias(:,:), omfBiasMpiGlobal(:,:)
    real(8), allocatable :: BMatrixMinusOne(:,:), LineVec(:,:)
    integer, allocatable :: OmFCount(:,:), omfCountMpiGlobal(:,:)
    character(len=80)    :: errorMessage

    Write(*,*) "bias_do_regression: start"
    if ( .not. allocated(trialHeight300m1000) ) then
      call bias_getTrialPredictors(obsSpaceData,columnhr)
      call bias_computePredictorBiases(obsSpaceData)
    end if

    call  bias_getRadiosondeWeight(obsSpaceData)

    SENSORS:do sensorIndex = 1, tvs_nsensors

      if  (.not. tvs_isReallyPresentMpiGLobal(sensorIndex)) cycle SENSORS

      nchans = bias(sensorIndex)%numChannels
      nscan = bias(sensorIndex)% numscan
      npredMax = maxval( bias(sensorIndex)%chans(:)%numActivePredictors )

      if ( lMimicSatbcor ) then
        ndimmax = npredMax
        allocate( OmFBias(nchans,nscan) )
        OmFBias(:,:) = 0.d0
      else
        ndimmax = npredMax + nscan - 1
      end if

      allocate( OmFCount(nchans,nscan) )
      OmFCount(:,:) = 0

      allocate( Matrix( nchans,ndimmax,ndimmax) )
      Matrix(:,:,:) = 0.d0

      allocate( Vector(nchans,ndimmax ) )
      Vector(:,:) = 0.d0

      allocate( LineVec(1,ndimmax ) )

      allocate(pIMatrix(ndimmax,ndimmax))


      ! First pass throught ObsSpaceData to estimate scan biases and count data

      call obs_set_current_header_list(obsSpaceData,'TO')
      HEADER1: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if ( headerIndex < 0 ) exit HEADER1
          
        ! process only radiance data to be assimilated?
        idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
        if ( .not. tvs_isIdBurpTovs(idatyp) ) cycle HEADER1
      
        indxtovs = tvs_tovsIndex(headerIndex)
        if ( indxtovs < 0 ) cycle HEADER1

        iSensor = tvs_lsensor( indxTovs )
        if (iSensor /= sensorIndex) cycle HEADER1
        iFov = obs_headElem_i(obsSpaceData,OBS_FOV,headerIndex)
        if ( nscan > 1 ) then
          iScan = iFov
        else
          iScan = 1
        end if

        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY1: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if ( bodyIndex < 0 ) exit BODY1
            
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY1 
          call bias_getChannelIndex(obsSpaceData,iSensor,chanIndx,bodyIndex)
          if (chanindx > 0) then
            if  (lMimicSatBcor) then
              OmF = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
              OmFBias(chanIndx,iScan) = OmFBias(chanIndx,iScan) + OmF
            end if
            OmFCount(chanIndx,iScan) =  OmFCount(chanIndx,iScan) + 1
          end if
        end do BODY1
      end do HEADER1

      if ( lMimicSatbcor ) allocate( omfBiasMpiGlobal(nchans,nscan) )

      allocate( omfCountMpiGlobal(nchans,nscan) )

      if ( lMimicSatbcor ) then
        call rpn_comm_reduce(OmFBias , omfBiasMpiGlobal, size(omfBiasMpiGlobal), "MPI_DOUBLE_PRECISION" , "MPI_SUM", 0, "GRID", ierr )
        if (ierr /=0) then
          Write(errorMessage,*) "bias_do_regression: MPI communication error 1",  ierr 
          call utl_abort(errorMessage)
        end if
      end if
      call rpn_comm_reduce(OmFCount , omfCountMpiGlobal, size(omfCountMpiGlobal), "MPI_INTEGER" , "MPI_SUM", 0, "GRID", ierr )

      if (ierr /=0) then
        Write(errorMessage,*) "bias_do_regression: MPI communication error 2",  ierr 
        call utl_abort(errorMessage)
      end if
      if ( lMimicSatbcor)  then
        if (mpi_myId == 0) then
          where( omfCountMpiGlobal == 0 ) omfBiasMpiGlobal = 0.d0
          where( omfCountMpiGlobal > 0 ) omfBiasMpiGlobal = omfBiasMpiGlobal / omfCountMpiGlobal
        end if
        call rpn_comm_bcast(omfBiasMpiGlobal, size(omfBiasMpiGlobal), "MPI_DOUBLE_PRECISION" , 0, "GRID",ierr )
        if (ierr /=0) then
          Write(errorMessage,*) "bias_do_regression: MPI communication error 3",  ierr 
          call utl_abort(errorMessage)
        end if
        do iChannel = 1, nchans
          bias(sensorIndex)%chans(iChannel)%coeff_fov(:) = omfBiasMpiGlobal(iChannel,:)
        end do
        deallocate( omfBiasMpiGlobal )
      end if
     
      ! Second pass to fill matrices and vectors
      call obs_set_current_header_list(obsSpaceData,'TO')
      iobs = 0
      HEADER2: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if ( headerIndex < 0 ) exit HEADER2

        ! process only radiance data to be assimilated?
        idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
        if ( .not. tvs_isIdBurpTovs(idatyp) ) cycle HEADER2

        iobs = iobs + 1
        indxtovs = tvs_tovsIndex(headerIndex)
        if ( indxtovs < 0 ) cycle HEADER2

        iSensor = tvs_lsensor( indxTovs )
        if (iSensor /= sensorIndex) cycle HEADER2
          
        call obs_set_current_body_list(obsSpaceData, headerIndex)
        iFov = obs_headElem_i(obsSpaceData,OBS_FOV,headerIndex)
        if ( nscan > 1 ) then
          iScan = iFov
        else
          iScan = 1
        end if
        BODY2: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if ( bodyIndex < 0 ) exit BODY2
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY2 

          call bias_getChannelIndex(obsSpaceData,iSensor,chanIndx,bodyIndex)
          if (chanIndx > 0) then
            call bias_getPredictors(predictor,headerIndex,iobs,chanIndx,obsSpaceData)
            OmF = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)

            if (lMimicSatbcor) OmF = OmF - bias(sensorIndex)%chans(chanIndx)%coeff_fov(iScan)
            
            LineVec(:,:) = 0.d0

            if (lMimicSatbcor) then
              idim = 0
              predstart = 1
              lambda = 1.d0
            else
              LineVec(1,iScan) = 1.d0
              idim = nscan
              predstart = 2
              sigmaObs = obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)
              lambda = 1.d0 / ( sigmaObs **2 )
            end if
          
            lambda = lambda * RadiosondeWeight(iobs)

            do iPred1 = predstart, bias(iSensor)%chans(chanIndx)%NumActivePredictors
              jPred1 = bias(iSensor)%chans(chanIndx)%PredictorIndex(iPred1)
              idim = idim + 1
              LineVec(1,idim) =  predictor(jPred1)
            end do

            Matrix(chanindx,:,:) = Matrix(chanindx,:,:) + matmul( transpose(LineVec),LineVec ) * lambda

            Vector(chanIndx,:) =  Vector(chanIndx,:) + LineVec(1,:) * OmF  * lambda
          end if
        end do BODY2
      end do HEADER2

      if (mpi_myId == 0) then
        allocate( matrixMpiGlobal(nchans,ndimmax,ndimmax) )
        allocate( vectorMpiGlobal(nchans,ndimmax ) )
      else
        allocate( matrixMpiGlobal(1,1,1) )
        allocate( vectorMpiGlobal(1,1)   )
      end if

      ! communication MPI pour tout avoir sur tache 0
      call rpn_comm_reduce(Matrix , matrixMpiGlobal, size(Matrix), "MPI_DOUBLE_PRECISION" , "MPI_SUM", 0, "GRID", ierr )
      if (ierr /=0) then
        Write(errorMessage,*) "bias_do_regression: MPI communication error 4",  ierr 
        call utl_abort(errorMessage)
      end if
      call rpn_comm_reduce(Vector , vectorMpiGlobal, size(Vector), "MPI_DOUBLE_PRECISION" , "MPI_SUM", 0, "GRID", ierr )
      if (ierr /=0) then
        Write(errorMessage,*) "bias_do_regression: MPI communication error 5",  ierr 
        call utl_abort(errorMessage)
      end if

      do iChannel = 1, nchans

        if (mpi_myId == 0) then
          ntot = sum(omfCountMpiGlobal(iChannel,:) )
          bias(sensorIndex)%chans(iChannel)%coeff_nobs = ntot
          if (ntot > 0 .and. .not. lMimicSatbcor ) then
            norm = 1.d0 / ( ntot ) 
            matrixMpiGlobal(iChannel,:,:) =  matrixMpiGlobal(iChannel,:,:) * norm
            vectorMpiGlobal(iChannel,:) = vectorMpiGlobal(iChannel,:) * norm
          end if

          nPred =  bias(sensorIndex)%chans(iChannel)%numActivePredictors
          if ( lMimicSatbcor ) then
            ndim = npred
          else
            ndim = npred + nscan -1 
            allocate( BMatrixMinusOne(ndim,ndim) )
            BMatrixMinusOne(:,:) = 0.d0 
            BMatrixMinusOne(1:nscan,1:nscan) =  matmul(bias(sensorIndex)%BMinusHalfScanBias, bias(sensorIndex)%BMinusHalfScanBias)
            do iPred1 = 2, bias(sensorIndex)%chans(iChannel)%numActivePredictors
              BMatrixMinusOne(nscan-1+iPred1,nscan-1+iPred1) =  ( 1.d0 / ( bias(sensorIndex)% chans(iChannel) % stddev(iPred1) )**2 )
            end do
            matrixMpiGlobal(iChannel,1:ndim,1:ndim) =   matrixMpiGlobal(iChannel,1:ndim,1:ndim) + BMatrixMinusOne(:,:)
            deallocate( BMatrixMinusOne)
          end if

          pIMatrix(:,:) = 0.d0
          call utl_pseudo_inverse(matrixMpiGlobal(iChannel,1:ndim,1:ndim), pIMatrix(1:ndim,1:ndim) )
          LineVec(1,1:ndim) = matmul(pIMatrix(1:ndim,1:ndim), vectorMpiGlobal(iChannel,1:ndim))
          !call dsymv("L", ndim, 1.d0, pIMatrix, ndim,vectorMpiGlobal(iChannel,:), 1, 0.d0, LineVec(1,1:ndim) , 1)
        end if

        call rpn_comm_bcast(ndim,  1, "MPI_INTEGER", 0, "GRID", ierr )
        call rpn_comm_bcast(npred, 1, "MPI_INTEGER", 0, "GRID", ierr )

        call rpn_comm_bcast(LineVec(1,1:ndim), ndim, "MPI_DOUBLE_PRECISION", 0, "GRID", ierr )
        if (ierr /=0) then
          Write(errorMessage,*) "bias_do_regression: MPI communication error 6",  ierr 
          call utl_abort(errorMessage)
        end if

        if (loutCoeffCov) then
          allocate ( bias(sensorIndex)%chans(iChannel)%coeffCov(ndim,ndim) ) 
          call rpn_comm_bcast(pIMatrix(1:ndim,1:ndim), ndim*ndim, "MPI_DOUBLE_PRECISION", 0, "GRID", ierr )
          bias(sensorIndex)%chans(iChannel)%coeffCov(:,:) = pIMatrix(1:ndim,1:ndim)
        end if

        if ( lMimicSatbcor ) then
          bias(sensorIndex)%chans(iChannel)%coeff(:) =  LineVec(1,1:npred)
        else
          bias(sensorIndex)%chans(iChannel)%coeff_fov(:) = LineVec(1,1:nscan)
          bias(sensorIndex)%chans(iChannel)%coeff(1) = 0.d0
          bias(sensorIndex)%chans(iChannel)%coeff(2:) =  LineVec(1,nscan+1:ndim)
        end if

      end do

      deallocate( LineVec )
      deallocate( Matrix )
      deallocate( Vector   )
      deallocate( omfCountMpiGlobal )
      deallocate( matrixMpiGlobal )
      deallocate( vectorMpiGlobal )
      deallocate( pIMatrix )
      if ( allocated(OmFBias) ) deallocate( OmFBias )
      deallocate( OmFCount )
       
    end do SENSORS

  end subroutine bias_do_regression


  !----------------------
  ! bias_Finalize
  !----------------------
  subroutine bias_Finalize
    !
    ! :Purpose: release allocated memory for the module
    !
    implicit none
    !Locals:
    integer    :: iSensor, iChannel

    if ( .not. biasActive ) return

    deallocate(trialHeight300m1000)
    deallocate(trialHeight50m200)
    deallocate(trialHeight1m10)
    deallocate(trialHeight5m50)
    deallocate(trialTG)
    deallocate(RadiosondeWeight)

    do iSensor = 1, tvs_nSensors
      if ( allocated( bias(iSensor)%BHalfScanBias) ) &
           deallocate( bias(iSensor)%BHalfScanBias )
      if ( allocated( bias(iSensor)%BMinusHalfScanBias) ) &
           deallocate( bias(iSensor)%BMinusHalfScanBias )
      do iChannel =1, bias(iSensor)% numChannels
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

  end subroutine bias_Finalize 
 
  !-----------------------------
  ! InstrNametoCoeffFileName 
  !-----------------------------
  function InstrNametoCoeffFileName(nameIn) result(nameOut)
    implicit none
    !Arguments:
    character(len=10), intent(in) :: nameIn
    character(len=10)             :: nameOut
    !Locals
    character(len=10)  :: temp_instrName
    integer            :: ierr

    temp_instrName = nameIn
    ierr = clib_tolower(temp_instrName) 
    if ( trim(temp_instrName) == 'mhs' ) then
      nameOut = 'amsub'
    else if ( trim(temp_instrName) == 'goesimager' ) then
      nameOut = 'cgoes'
    else if ( trim(temp_instrName) == 'gmsmtsat' ) then
      nameOut = 'mtsat'
    else if ( trim(temp_instrName) == 'mviri' ) then
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
    
    character(len=10), intent(in) :: nameIn
    character(len=10)             :: nameOut

    if ( trim(nameIn) == 'MHS' ) then
      nameOut = 'AMSUB'
    else if ( trim(nameIn) == 'GOESIMAGER' ) then
      nameOut = 'CGOES' 
    else if ( trim(nameIn) == 'GMSMTSAT' ) then
      nameOut = 'MTSAT' 
    else if ( trim(nameIn) == 'MVIRI' ) then
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
    !Arguments:
    character(len=10), intent(in) :: nameIn
    character(len=10)             :: nameOut

    if ( trim(nameIn) == 'MSG1' ) then
      nameOut = 'METSAT8'
    else if ( trim(nameIn) == 'MSG2' ) then
      nameOut = 'METSAT9'
    else if ( trim(nameIn) == 'MSG3' ) then
      nameOut = 'METSAT10'
    else if ( trim(nameIn) == 'MSG4' ) then
      nameOut = 'METSAT11'
    else if ( trim(nameIn) == 'METEOSAT7' ) then
      nameOut = 'METSAT7' 
    else 
      nameOut = nameIn
    end if

  end function SatNameinCoeffFile

  
  !-----------------------------------------
  ! bias_read_bcif
  !-----------------------------------------
  subroutine read_bcif(bcifFile, hspec, ncan, can, bcmode, bctype, npred, pred, global_opt, exitcode)
    !
    ! :Purpose: to read channel-specific bias correction (BC) information (predictors) for instrument from BCIF.
    !
    implicit none
    !Arguments:
    character(len=*), intent(in)  :: bcifFile
    logical, intent(in)           :: hspec
    integer, intent(out)          :: exitcode, ncan
    integer, intent(out)          :: can(tvs_maxchannelnumber), npred(tvs_maxchannelnumber)
    character(len=3), intent(in)  :: global_opt
    character(len=1), intent(out) :: bcmode(tvs_maxchannelnumber), bctype(tvs_maxchannelnumber)
    character(len=2), intent(out) :: pred(tvs_maxchannelnumber,numpredictors)
    !Locals:
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
    ier = fnom(iun,bcifFile,'FMT',0)
    if ( ier /= 0 ) then
      call utl_abort('read_bcif: ERROR - Problem opening the bcif file!' // trim(bcifFile) ) 
    end if
    
    read(iun, *) instrum, ncan
    read(iun, '(A64)') line

    ! For GLOBAL option, read global values from first line (channel 0) and clone to all channels
    if ( global_opt == 'OUI' .or. global_opt == 'DEF' ) then 
      ! Read channel 0 information
      read(iun,*,IOSTAT=ier) can(1), bcmode(1), bctype(1), npred(1), (pred(1,j), j=1,numpredictors)
      if ( ier /= 0 ) then
        write(*,*) 'read_BCIF: Error reading channel 0 data!'
        exitcode = ier
        return
      end if
      if ( can(1) /= 0 ) then
        write(*,*) 'read_BCIF: Channel 0 global values not found!'
        exitcode = -1
        return
      end if
      ! Clone channel 0 information to all ncan channels
      if ( .not. hspec ) then
        ! For instruments with consecutive channels 1,2,3,...ncan (e.g. AMSU, SSM/I)
        !  -- no need to read the channel numbers from the BCIF
        do i = 2, ncan+1
          can(i)    = i-1
          bcmode(i) = bcmode(1)
          bctype(i) = bctype(1)
          npred(i)  = npred(1)
          do j = 1, numpredictors
            pred(i,j) = pred(1,j)
          end do
        end do
      else
        ! For hyperspectral instruments (channel subsets), read the channel numbers from the BCIF
        do i = 2, ncan+1
          read(iun,*,IOSTAT=ier) xcan
          if ( ier /= 0 ) then
            write(*,*) 'read_BCIF: Error reading channel numbers!'
            exitcode = ier
            return
          end if
          can(i)    = xcan
          bcmode(i) = bcmode(1)
          bctype(i) = bctype(1)
          npred(i)  = npred(1)
          do j = 1, numpredictors
            pred(i,j) = pred(1,j)
          end do
        end do
        ! Reposition the file to just after channel 0 record
        REWIND (UNIT=iun)
        read(iun, *) instrum, ncan
        read(iun, '(A64)') line
        read(iun,*,IOSTAT=ier) xcan, xbcmode, xbctype, xnpred, (xpred(j), j=1,numpredictors)
      end if
      ! For global_opt == 'DEF' check for channel-specific information and overwrite the default (channel 0) values
      ! for the channel with the values from the file
      if ( global_opt == 'DEF' ) then
        if ( .not. hspec ) then
          do
            read(iun,*,IOSTAT=ier) xcan, xbcmode, xbctype, xnpred, (xpred(j), j=1,numpredictors)
            if ( ier < 0 ) exit  
            if ( ier > 0 ) then
              write(*,*) 'read_BCIF: Error reading file!'
              exitcode = ier
              return
            end if
            ii = xcan+1
            if ( ii > ncan+1 ) then
              write(*,*) 'read_BCIF: Channel number in BCIF exceeds number of channels!'
              write(*,'(A,1X,I4,1X,I4)') '           Channel, ncan = ', xcan, ncan
              exitcode = -1
              return
            end if
            bcmode(ii) = xbcmode
            bctype(ii) = xbctype
            npred(ii)  = xnpred
            do j = 1, numpredictors
              pred(ii,j) = xpred(j)
            end do
          end do
        else
          ! For hyperspectral instruments
          do i = 2, ncan+1
            read(iun,*,IOSTAT=ier) xcan, xbcmode, xbctype, xnpred, (xpred(j), j=1,numpredictors)
            if ( ier /= 0 ) CYCLE  
            bcmode(i) = xbcmode
            bctype(i) = xbctype
            npred(i)  = xnpred
            do j = 1, numpredictors
              pred(i,j) = xpred(j)
            end do
          end do
        end if
      end if
    end if
    ! Non-GLOBAL: Read the entire file for channel specific values (all channels)
    ! ---------------------------------------------------------------------------------------------------------------------
    if ( global_opt == 'NON' ) then
      ii = 1
      do
        read(iun,*,IOSTAT=ier) can(ii), bcmode(ii), bctype(ii), npred(ii), (pred(ii,j), j=1,numpredictors)
        if ( ier < 0 ) exit  
        if ( ier > 0 ) then
          write(*,*) 'read_BCIF: Error reading file!'
          exitcode = ier
          return
        end if
        if ( ii == 1 ) then
          if ( can(ii) /= 0 ) then
            write(*,*) 'read_BCIF: Channel 0 global/default values not found!'
            exitcode = -1
            return
          end if
        end if
        ii = ii + 1
      end do
      if ( ii-2 < ncan ) then
        write(*,*) 'read_BCIF: Number of channels in file is less than specified value (NCAN). Changing value of NCAN.'
        ncan = ii - 2
      end if
      if ( ii > ncan+2 ) then
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
    do i = 1, ncan+1
      chknp = count(pred(i,:) /= 'XX')
      if ( chknp /= npred(i) ) npred(i) = chknp
      if ( npred(i) == 0 .and. bctype(i) == 'C' ) bctype(i) = 'F'
      write(*,'(I4,2(5X,A1),5X,I2,6(4X,A2))') can(i), bcmode(i), bctype(i), npred(i), (pred(i,j), j=1,numpredictors)
    end do
    write(*,*) 'read_BCIF: exit '

    ier = fclos(iun)
    exitcode = 0

  end subroutine read_bcif

  
  !-----------------------------------------
  ! bias_applyBiasCorrection
  !-----------------------------------------
  subroutine read_coeff(sats,chans,fovbias,coeff,nsat,nchan,nfov,npred,cinstrum,coeff_file,ptypes,ndata)
    !
    ! :Purpose: to read radiance bias correction coefficients file
    !
    implicit none
    !Arguments:
    character(len=10), intent(out) :: sats(:)       ! dim(maxsat), satellite names 1
    integer, intent(out)           :: chans(:,:)    ! dim(maxsat, maxchan), channel numbers 2
    real(8), intent(out)           :: fovbias(:,:,:)! dim(maxsat,maxchan,maxfov), bias as F(fov) 3
    real(8), intent(out)           :: coeff(:,:,:)  ! dim(maxsat,maxchan,maxpred+1) 4
    integer, intent(out)           :: nsat          !5
    integer, intent(out)           :: nchan(:)      ! dim(maxsat), number of channels 6
    integer, intent(out)           :: nfov          !7
    integer, intent(out)           :: npred(:,:)    ! dim(maxsat, maxchan), number of predictors !8
    character(len=7), intent(out)  :: cinstrum      ! string: instrument (e.g. AMSUB) 9
    character(len=*), intent(in)   :: coeff_file    ! 10
    character(len=2), intent(out)  :: ptypes(:,:,:) ! dim(maxsat,maxchan,maxpred) 11
    integer, intent(out)           :: ndata(:,:)    ! dim(maxsat, maxchan), number of channels 12
    !Locals:
    character(len=8)               :: sat
    character(len=120)             :: line
    integer                        :: chan
    integer                        :: nbfov, nbpred, i, j, k, ier, istat, ii, nobs
    logical                        :: newsat
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

    iun = 0
    ier = fnom(iun,coeff_file,'FMT',0)
    if ( ier /= 0 ) then
      call utl_abort('read_coeff: ERROR - Problem opening the coefficient file! ' // trim(coeff_file) )
    end if

    write(*,*)
    write(*,*) 'read_coeff: Bias correction coefficient file open = ', coeff_file

    maxsat =  size( sats )
    maxpred = size(ptypes, dim=3)

    coeff(:,:,:)    = 0.0
    fovbias(:,:,:)  = 0.0
    sats(:)         = 'XXXXXXXX'
    cinstrum        = 'XXXXXXX'
    chans(:,:)      = 0
    npred(:,:)      = 0
    nsat            = 0
    nchan(:)        = 0
    nfov            = 0
    ptypes(:,:,:)   = 'XX'
    
    read(iun,*,IOSTAT=istat)
    if ( istat < 0 ) then
      write(*,*) 'read_coeff: ERROR- File appears empty.'
      return
    end if
    rewind(iun)

    ii = 0

    ! Loop over the satellites/channels in the file

    do
      read(iun,'(A)',IOSTAT=istat) line
      if ( istat < 0 ) exit
      if ( line(1:3) == 'SAT' ) then
        newsat = .true.
        read(line,'(T53,A8,1X,A7,1X,I6,1X,I8,1X,I2,1X,I3)',IOSTAT=istat) sat, cinstrum, chan, nobs, nbpred, nbfov
        if ( istat /= 0 ) then
          call utl_abort('read_coeff: ERROR - reading data from SATELLITE line in coeff file!')
        end if
        do i = 1, maxsat
          if ( trim(sats(i)) == trim(sat) ) then
            newsat = .false.
            ii = i
          end if
        end do
        if ( newsat ) then
          ii = ii + 1
          if ( ii > maxsat ) then
            call utl_abort('read_coeff: ERROR - max number of satellites exceeded in coeff file!')
          end if
          sats(ii) = sat
          if (ii > 1) nchan(ii-1) = j
          j = 1
        else
          j = j + 1
        end if
        chans(ii, j) = chan
        npred(ii, j) = nbpred
        ndata(ii, j) = nobs
        if ( nbpred > maxpred ) then
          call utl_abort('read_coeff: ERROR - max number of predictors exceeded in coeff file!')
        end if
        read(iun,'(A)',IOSTAT=istat) line
        if ( line(1:3) /= 'PTY' ) then
          call utl_abort('read_coeff: ERROR - list of predictors is missing in coeff file!')
        end if
        if ( nbpred > 0 ) then
          read(line,'(T8,6(1X,A2))',IOSTAT=istat) (ptypes(ii,j,k),k=1,nbpred)
          if ( istat /= 0 ) then
            call utl_abort('read_coeff: ERROR - reading predictor types from PTYPES line in coeff file!')
          end if
        end if
        read(iun,*,IOSTAT=istat) (fovbias(ii,j,k),k=1,nbfov)
        if ( istat /= 0 ) then
          call utl_abort('read_coeff: ERROR - reading fovbias in coeff file!')
        end if
        if ( nbpred > 0 ) then
          read(iun,*,IOSTAT=istat) (coeff(ii,j,k),k=1,nbpred+1)
        else
          read(iun,*,IOSTAT=istat) dummy
        end if
        if ( istat /= 0 ) then
          call utl_abort('read_coeff: ERROR - reading coeff in coeff file!')
        end if
      else
        exit
      end if

    end do

    if ( ii == 0 ) then
      call utl_abort('read_coeff: ERROR - No data read from coeff file!')
    end if

    nsat      = ii
    nfov      = nbfov
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
        write(*,*) i, chans(i,j)
        if ( npred(i,j) > 0 ) then
          write(*,'(6(1X,A2))') (ptypes(i,j,k),k=1,npred(i,j))
        else
          write(*,'(A)') 'read_coeff: No predictors'
        end if
        write(*,*) (fovbias(i,j,k),k=1,nfov)
        write(*,*) (coeff(i,j,k),k=1,npred(i,j)+1)
      end do
    end do
    write(*,*) 'read_coeff: exit'

    ier = fclos(iun)

  end subroutine read_coeff

  !-----------------------------------------
  ! bias_getChannelIndex
  !-----------------------------------------
  subroutine bias_getChannelIndex(obsSpaceData, idsat, chanIndx,indexBody)
    !
    ! :Purpose: to get the channel index (wrt bcif channels)
    !
    implicit none
    !Arguments:
    integer, intent(in)  :: idsat, indexBody
    integer, intent(out) :: chanIndx
    type(struct_obs),intent(inout)  :: obsSpaceData
    !Locals:
    logical, save :: first =.true.
    integer :: ichan, isensor, indx 
    integer, allocatable, save :: Index(:,:)
    
    if (first) then
      allocate( Index(tvs_nsensors, tvs_maxChannelNumber ) )
      Index(:,:) = -1
      do isensor = 1, tvs_nsensors
        channels:do ichan = 1,  tvs_maxChannelNumber
          indexes: do indx =1, bias(isensor) % numChannels
            if ( ichan == bias(isensor) %chans(indx)%channelNum ) then
              Index(isensor,ichan) = indx
              exit indexes
            end if
          end do indexes
        end do channels
      end do
      first = .false.
    end if
    
    ichan = nint( obs_bodyElem_r(obsSpaceData,OBS_PPP,indexBody) )
    ichan = max(0,min(ichan,tvs_maxChannelNumber+1))
    ichan = ichan - tvs_channelOffset(idsat)
    
    chanIndx = Index(idsat,ichan)

  end subroutine bias_getChannelIndex


end MODULE biascorrection_mod


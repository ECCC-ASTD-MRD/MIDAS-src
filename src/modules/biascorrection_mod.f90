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

  implicit none
  save
  private

  public               :: bias_setup,bias_calcBias_tl,bias_calcBias_ad, bias_writeBias, bias_finalize

  type  :: struct_bias
    integer  :: numActivePredictors, numChannels, numScan
    real(8),allocatable  :: stddev(:,:)
    real(8),allocatable  :: coeffIncr(:,:)
    integer,allocatable  :: predictorIndex(:)
    real(8),allocatable  :: coeffIncr_fov(:,:)
    integer,allocatable  :: channelNum(:)
  end type struct_bias

  type(struct_bias),allocatable  :: bias(:)
   
  logical               :: initialized = .false.
  integer, parameter    :: maxNumChannels = tvs_maxChannelNumber
  integer, parameter    :: NumPredictors = 11
  integer, parameter    :: maxfov = 120
  real(8), allocatable  :: trialHeight300m1000(:)
  real(8), allocatable  :: trialHeight50m200(:)
  real(8), allocatable  :: trialHeight1m10(:)
  real(8), allocatable  :: trialHeight5m50(:)
  real(8), allocatable  :: trialTG(:)
  integer               :: nobs

  logical  :: lvarbc
  real(8)  :: bg_stddev(NumPredictors)
  logical  :: abortIfNoBackground
  namelist /nambias/lvarbc,bg_stddev,abortIfNoBackground

CONTAINS
 
  !-----------------------------------------------------------------------
  ! bias_setup
  !-----------------------------------------------------------------------
  SUBROUTINE bias_setup()
    implicit none

    integer  :: cvdim
    integer  :: iSensor,iPredictor
    integer  :: jPredictor
    integer  :: ierr,nulnam
    integer  :: fnom,fclos
    integer  :: jPred 
    integer  :: idNum(tvs_nSensors,NumPredictors)
    integer  :: activePredictors(tvs_nSensors)
    character(len=85)  :: filecoeff
    character(len=10)  :: instrName, instrNamecoeff, satNamecoeff 
    logical            :: coeffExists 

    !variables from background coeff file
    character(len=10)  :: sats(tvs_nSensors)
    integer            :: chans(tvs_nSensors, tvs_maxChannelNumber), nsat, nfov, iSat
    integer            :: nchan(tvs_nSensors)
    character(len=6)   :: cinstrum

    ! set default values for namelist variables
    lvarbc   = .false.
    bg_stddev(:) = 0.0d0
    abortIfNoBackground = .true.

    ! read in the namelist NAMBIAS
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=nambias,iostat=ierr)
    if ( ierr /= 0 .and. mpi_myid == 0 ) write(*,*) 'WARNING: bias_setup: Error reading namelist, ' //  &
                             'assume it will not be used!'
    if ( mpi_myid == 0 ) write(*,nml=nambias)
    ierr = fclos(nulnam)

    cvdim = 0

    if ( lvarbc ) then
      allocate(bias(tvs_nSensors))

      do iSensor = 1, tvs_nSensors

        bias(iSensor)%numActivePredictors = 0
        activePredictors(iSensor) = 0
        jPred = 0

        ! set predictors, numScan and numChannels
        iPredictor = 1  ! constant predictor for all sensors
        activePredictors(iSensor)=ibset(activePredictors(iSensor),iPredictor-1)
        iPredictor=2  ! height300-height1000
        activePredictors(iSensor)=ibset(activePredictors(iSensor),iPredictor-1)
        iPredictor=3  ! height50-height200
        activePredictors(iSensor)=ibset(activePredictors(iSensor),iPredictor-1)
        iPredictor=7  ! height5-height50
        activePredictors(iSensor)=ibset(activePredictors(iSensor),iPredictor-1)

        do iPredictor=1,NumPredictors
          if ( btest(activePredictors(iSensor),iPredictor-1) ) then
            bias(iSensor)%numActivePredictors = bias(iSensor)%numActivePredictors + 1
            jPred = bias(iSensor)%numActivePredictors
            idNum(iSensor,jPred) = iPredictor
          end if
        end do

        instrName = InstrNametoCoeffFileName(tvs_instrumentName(iSensor))
        instrNamecoeff = InstrNameinCoeffFile(tvs_instrumentName(iSensor))
        satNamecoeff = SatNameinCoeffFile(tvs_satelliteName(iSensor))

        filecoeff = 'coeff_file_'//trim(instrName)//''
        inquire(file=trim(filecoeff),exist = coeffExists)
         
        if ( coeffExists ) then
          call bias_updateCoeff(tvs_nSensors,NumPredictors,filecoeff,sats,chans,nsat,nchan,nfov,cinstrum, updateCoeff_opt = .false.)
          do iSat = 1, nsat
            if ( sats(iSat) /= trim(satNamecoeff) .or. cinstrum /= trim(instrNamecoeff) ) cycle
            bias(iSensor)%numScan = nfov
            allocate( bias(iSensor)%channelNum(nchan(iSat)) ) 
            bias(iSensor)%channelNum(:)  = chans(iSat,:)
            bias(iSensor)%numChannels = maxval(bias(iSensor)%channelNum(:))
          end do
        end if
        
        cvdim = cvdim + bias(iSensor)%numChannels*(bias(iSensor)%numActivePredictors-1) &
                      + bias(iSensor)%numChannels*bias(iSensor)%numScan
        write(*,*) 'bias_setup: Instrument=',iSensor,tvs_satelliteName(iSensor),tvs_instrumentName(iSensor), &
                   ', numActivePredictors=',bias(iSensor)%numActivePredictors, &
                   ', activePredictors=',(btest(activePredictors(iSensor),iPredictor-1),iPredictor=1,NumPredictors)
      end do

      do iSensor =  1, tvs_nSensors

        allocate( bias(iSensor)%stddev(1:bias(iSensor)%numChannels, 1:bias(iSensor)%numActivePredictors))
        allocate( bias(iSensor)%coeffIncr(1:bias(iSensor)%numChannels, 2:bias(iSensor)%numActivePredictors))
        allocate( bias(iSensor)%coeffIncr_fov(1:bias(iSensor)%numChannels, 1:bias(iSensor)%numScan))
        allocate( bias(iSensor)%predictorIndex (1:bias(iSensor)%numActivePredictors))

        do iPredictor = 1, bias(iSensor)%numActivePredictors
          bias(iSensor)%predictorIndex(iPredictor) = idNum(iSensor,iPredictor)
        end do   
       
        do iPredictor = 1, bias(iSensor)%numActivePredictors
          bias(iSensor)%stddev(:,iPredictor) = bg_stddev( bias(iSensor)%PredictorIndex(iPredictor) )
        end do

        ! Make AMSU-A ch 13-14 static
        if ( tvs_instruments(iSensor) == 3 ) then
          if ( bias(iSensor)%numChannels /= 0 ) then
            bias(iSensor)%stddev(13:14, :) = 0.0d0
          end if
        end if
      end do
    end if

    if ( cvdim > 0 ) then
      if ( mpi_myid > 0 ) cvdim = 0 ! for minimization, all coefficients only on task 0
      call cvm_setupSubVector('BIAS', 'BIAS', cvdim)
    end if

  END SUBROUTINE bias_setup

  !---------------------------------------
  ! bias_calcBias_tl
  !---------------------------------------- 
  SUBROUTINE bias_calcBias_tl(cv_in,cv_dim,obsColumnIndex,obsSpaceData,columnhr)
    implicit none

    integer  :: cv_dim
    real(8)  :: cv_in(cv_dim)
    integer  :: obsColumnIndex
    type(struct_obs)  :: obsSpaceData
    type(struct_columnData) :: columnhr

    integer  :: index_header,index_body,iobs, indxtovs, idatyp
    integer  :: iSensor,iChannel,iPredictor,index_cv
    integer  :: iScan, iFov, jPred
    real(8)  :: predictor(NumPredictors)
    real(8),pointer  :: cv_bias(:)
    real(8)  :: biasCor
    logical,save  :: firstTime=.true.

    if ( .not. lvarbc ) return

    if ( firstTime ) then
      call bias_getTrialPredictors(obsSpaceData,columnhr)
      call bias_calcMeanPredictors(obsSpaceData)
      firstTime = .false.
    end if

    if ( mpi_myid == 0) then
      if ( cvm_subVectorExists('BIAS') ) then
        cv_Bias => cvm_getSubVector(cv_in,'BIAS')
        write(*,*) 'bias_calcBias_tl: maxval(cv_bias)=',maxval(cv_bias(:))
      else
        write(*,*) 'bias_calcBias_tl: control vector does not include bias coefficients'
        return
      end if
    end if

    ! get bias coefficients
    call bias_cvToCoeff(cv_bias)

    ! apply bias increment to specified obs column
    iobs = 0
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      index_header = obs_getHeaderIndex(obsSpaceData)
      if ( index_header < 0 ) exit HEADER

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,index_header)
      if ( .not. tvs_isIdBurpTovs(idatyp) ) cycle HEADER
      
      indxtovs = tvs_ltovsno(index_header)
      if ( indxtovs == 0 ) then
        call utl_abort('bias_calcBias_tl')
      end if

      iobs = iobs + 1
      iSensor = tvs_lsensor(tvs_ltovsno(index_header))
      call bias_getPredictors(predictor,index_header,iobs,obsSpaceData)

      call obs_set_current_body_list(obsSpaceData, index_header)
      iFov = obs_headElem_i(obsSpaceData,OBS_FOV,index_header)

      BODY: do
        index_body = obs_getBodyIndex(obsSpaceData)
        if ( index_body < 0 ) exit BODY

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) /= obs_assimilated ) cycle BODY   

        iChannel = nint(obs_bodyElem_r(obsSpaceData,OBS_PPP,index_body))
        iChannel = max(0,min(iChannel,tvs_maxChannelNumber+1))
        iChannel = iChannel-tvs_channelOffset(iSensor)
        biasCor = 0.0d0

        do iPredictor = 1, bias(iSensor)%NumActivePredictors
          jPred = bias(iSensor)%PredictorIndex(iPredictor)
          if ( iPredictor == 1 ) then

            if ( bias(iSensor)%numScan > 1 ) then
              iScan = iFov
            else
              iScan = 1
            end if
            biasCor = biasCor + predictor(jPred) * bias(iSensor)%coeffIncr_fov(iChannel,iScan) 
          else
            biasCor = biasCor + predictor(jPred) * bias(iSensor)%coeffIncr(iChannel,iPredictor) 
          end if
          
        end do
        call obs_bodySet_r( obsSpaceData, obsColumnIndex, index_body, &
          obs_bodyElem_r(obsSpaceData,obsColumnIndex,index_body) - biasCor) 
      end do BODY
    end do HEADER

  END SUBROUTINE bias_calcBias_tl

 
  !----------------------
  ! bias_getTrialPredictors
  !----------------------
  SUBROUTINE bias_getTrialPredictors(obsSpaceData,columnhr)
    implicit none
    type(struct_columnData) :: columnhr
    type(struct_obs)  :: obsSpaceData
    integer  :: index_header, idatyp, nobs
    real(8)  :: height1, height2
 
    ! count number of tovs locations
    nobs = 0
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      index_header = obs_getHeaderIndex(obsSpaceData)
      if ( index_header < 0 ) exit HEADER
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,index_header)
      if ( .not.  tvs_isIdBurpTovs(idatyp) ) cycle HEADER 
      nobs = nobs + 1
    end do HEADER
    
    if ( nobs > 0 ) then
      allocate(trialHeight300m1000(nobs))
      allocate(trialHeight50m200(nobs))
      allocate(trialHeight5m50(nobs))
      allocate(trialHeight1m10(nobs))
      allocate(trialTG(nobs))
    else
      write(*,*) 'bias_getTrialPredictors: NO OBS found'
      return
    end if

    nobs = 0

    call obs_set_current_header_list(obsSpaceData,'TO')

    HEADER2: do
      index_header = obs_getHeaderIndex(obsSpaceData)
      if ( index_header < 0 ) exit HEADER2
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,index_header)
      if ( .not.  tvs_isIdBurpTovs(idatyp) ) cycle HEADER2 
      nobs = nobs + 1

      height1 = logInterpHeight(columnhr,index_header,1000.d0)
      height2 = logInterpHeight(columnhr,index_header,300.d0)
      
      trialHeight300m1000(nobs) = height2 - height1

      height1 = logInterpHeight(columnhr,index_header,200.d0)
      height2 = logInterpHeight(columnhr,index_header,50.d0)

      trialHeight50m200(nobs) = height2 - height1

      height1 = height2
      height2 = logInterpHeight(columnhr,index_header,5.d0)

      trialHeight5m50(nobs) = height2 - height1

      height1 = logInterpHeight(columnhr,index_header,10.d0)
      height2 = logInterpHeight(columnhr,index_header,1.d0)

      trialHeight1m10(nobs) = height2 - height1

      trialTG(nobs) = col_getElem(columnhr,1,index_header,'TG')

    end do HEADER2

    if ( trialTG(1) > 150.0d0) then
      write(*,*) 'bias_getTrialPredictors_forTG: converting TG from Kelvin to deg_C'
      trialTG(:) = trialTG(:) - MPC_K_C_DEGREE_OFFSET_R8
    end if
 
    trialHeight300m1000(:) = 0.1d0 * trialHeight300m1000(:) ! conversion from to
    trialHeight50m200(:) = 0.1d0 * trialHeight50m200(:)
    trialHeight5m50(:) = 0.1d0 * trialHeight5m50(:)
    trialHeight1m10(:) =  0.1d0 *  trialHeight1m10(:)

    write(*,*) 'bias_getTrialPredictors done'

  contains

    function logInterpHeight(columnhr,headerIndex,P) result(height)
      integer,intent(in) :: headerIndex
      Real(8), intent(in) :: P
      Real(8) :: height
      type(struct_columnData),intent(inout) :: columnhr

      integer :: jk, nlev, ik
      real(8) :: zpt, zpb, zwt, zwb
      real(8),pointer :: col_ptr(:)

      IK = 1
      nlev = COL_GETNUMLEV(COLUMNHR,'TH')
      do JK = 2,NLEV - 1
        ZPT = col_getPressure(COLUMNHR,JK,headerIndex,'TH')* MPC_MBAR_PER_PA_R8
        if( P > ZPT ) IK = JK
      end do
      ZPT = col_getPressure(COLUMNHR,IK,headerIndex,'TH') * MPC_MBAR_PER_PA_R8
      ZPB = col_getPressure(COLUMNHR,IK+1,headerIndex,'TH') * MPC_MBAR_PER_PA_R8

      zwb  = log(p/zpt) / log(zpb/zpt)
      zwt  = 1.d0 - zwb
      col_ptr=>col_getColumn(columnhr,headerIndex,'Z_T')

      height = zwb * col_ptr(ik+1) + zwt * col_ptr(ik)
   
    end function logInterpHeight

  END SUBROUTINE bias_getTrialPredictors


  !---------------------------------------
  ! bias_cvToCoeff
  !--------------------------------------
  SUBROUTINE bias_cvToCoeff(cv_bias)
    implicit none

    real(8)  :: cv_bias(:)
    integer  :: index_cv,iSensor,iChannel,iPredictor,iScan
    integer  :: nsize,ierr

    if ( mpi_myid == 0 ) write(*,*) 'bias_cvToCoeff starting'
 
    if ( mpi_myid == 0 ) then

      index_cv = 0
      ! initialize of coeffIncr
      do iSensor = 1, tvs_nSensors
        bias(iSensor)%coeffIncr(:,:) = 0.0d0
        bias(iSensor)%coeffIncr_fov(:,:) = 0.0d0
      end do

      do iSensor = 1, tvs_nSensors
        do iChannel = 1, bias(iSensor)%numChannels
          do iPredictor = 1, bias(iSensor)%numActivePredictors
            if ( iPredictor == 1 ) then
              do iScan = 1, bias(iSensor)%numScan
                index_cv = index_cv + 1
                bias(iSensor)%coeffIncr_fov(iChannel,iScan) = bias(iSensor)%stddev(iChannel,iPredictor)*cv_bias(index_cv)
              end do
            else
              index_cv = index_cv + 1
              bias(iSensor)%coeffIncr(iChannel,iPredictor) = bias(iSensor)%stddev(iChannel,iPredictor)*cv_bias(index_cv)
            end if
          end do !iPredictor
        end do !iChannel
      end do !iSensor
    end if
   
    
    ! for constant part
    do iSensor = 1, tvs_nSensors
      nsize = bias(iSensor)%numScan*bias(iSensor)%numChannels
      call rpn_comm_bcast(bias(iSensor)%coeffIncr_fov,nsize,"mpi_double_precision",0,"GRID",ierr)
    end do

    ! for predictor part
    do iSensor = 1, tvs_nSensors
      nsize = (bias(iSensor)%numActivePredictors-1)*bias(iSensor)%numChannels
      call rpn_comm_bcast(bias(iSensor)%coeffIncr,nsize,"mpi_double_precision",0,"GRID",ierr)
    end do

  END SUBROUTINE bias_cvToCoeff

  !-----------------------------------
  ! bias_getPredictors
  !---------------------------------- 
  SUBROUTINE bias_getPredictors(predictor,index_header,index_obs,obsSpaceData)
    implicit none

    real(8)  :: predictor(NumPredictors)
    integer  :: index_header,index_obs
    type(struct_obs)  :: obsSpaceData

    integer  :: iSensor,iPredictor
    integer  :: activePredictors(tvs_nSensors)

    predictor(:) = 0.0d0
    iSensor = tvs_lsensor(tvs_ltovsno(index_header))

    do iPredictor = 1, NumPredictors

      if ( iPredictor == 1 ) then
        ! constant
        predictor(iPredictor) = 1.0d0
      else if ( iPredictor == 2 ) then
        ! height300-height1000 (dam) /1000
        predictor(iPredictor) = trialHeight300m1000(index_obs)/1000.0d0  !-0.85d0 tried de-biasing
      else if ( iPredictor == 3 ) then
        ! height50-height200 (dam) /1000
        predictor(iPredictor) = trialHeight50m200(index_obs)/1000.0d0  !-0.85d0 tried de-biasing
      else if ( iPredictor == 4 ) then
        ! skin temperature (C) /10
        predictor(iPredictor) = trialTG(index_obs)
      else if ( iPredictor == 6 ) then
        ! height1-height10 (dam) /1000
        predictor(iPredictor) = trialHeight1m10(index_obs)/1000.0d0
      else if ( iPredictor == 7 ) then
        ! height5-height50 (dam) /1000
        predictor(iPredictor) = trialHeight5m50(index_obs)/1000.0d0  !-0.72d0 tried de-biasing
      else if ( iPredictor == 9 ) then
        ! satellite zenith angle (degrees/100)^1
        predictor(iPredictor) = obs_headElem_r(obsSpaceData,OBS_SZA,index_header)
      else if ( iPredictor == 10 ) then
        ! satellite zenith angle (degrees/100)^2
        predictor(iPredictor) = ( obs_headElem_r(obsSpaceData,OBS_SZA,index_header) )**2
      else if ( iPredictor == 11 ) then
        ! satellite zenith angle (degrees/100)^3
        predictor(iPredictor) = ( obs_headElem_r(obsSpaceData,OBS_SZA,index_header) )**3
      end if

    end do

  END SUBROUTINE bias_getPredictors

  !--------------------------------------------------
  ! bias_calcMeanPredictors
  !---------------------------------------------------
  SUBROUTINE bias_calcMeanPredictors(obsSpaceData)
    implicit none

    type(struct_obs)  :: obsSpaceData

    integer  :: index_header,iobs,idatyp
    integer  :: iSensor,iPredictor
    integer  :: countObs(NumPredictors)
    real(8)  :: meanAbsPredictors(NumPredictors)
    real(8)  :: meanPredictors(NumPredictors)
    real(8)  :: minPredictors(NumPredictors)
    real(8)  :: maxPredictors(NumPredictors)
    real(8)  :: predictor(NumPredictors)

    countObs(:) = 0
    meanAbsPredictors(:) = 0.0d0
    meanPredictors(:) = 0.0d0
    minPredictors(:) = 999999.0d0
    maxPredictors(:) = -999999.0d0
    iSensor = -1

    ! loop over all observations
    iobs = 0
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      index_header = obs_getHeaderIndex(obsSpaceData)
      if ( index_header < 0 ) exit HEADER
      
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,index_header)
      if ( .not.  tvs_isIdBurpTovs(idatyp) ) cycle HEADER       

      iobs = iobs + 1
      iSensor = tvs_lsensor(tvs_ltovsno(index_header))
      call bias_getPredictors(predictor,index_header,iobs,obsSpaceData)

      do iPredictor = 2, bias(iSensor)%NumActivePredictors
        minPredictors(iPredictor) = min(minPredictors(iPredictor),predictor(bias(iSensor)%predictorIndex(iPredictor)))
        maxPredictors(iPredictor) = max(maxPredictors(iPredictor),predictor(bias(iSensor)%predictorIndex(iPredictor)))
        meanPredictors(iPredictor) = meanPredictors(iPredictor) + predictor(bias(iSensor)%predictorIndex(iPredictor))
        meanAbsPredictors(iPredictor) = meanAbsPredictors(iPredictor) + abs(predictor(bias(iSensor)%predictorIndex(iPredictor)))
        countObs(iPredictor) = countObs(iPredictor) + 1
      end do

    end do HEADER

    if ( iSensor /= -1 ) then
      !  do iPredictor=1,bias(iSensor)%NumActivePredictors
      do iPredictor = 2,bias(iSensor)%NumActivePredictors
        if ( countObs(iPredictor) > 0 ) then
          meanPredictors(iPredictor) = meanPredictors(iPredictor)/real(countObs(iPredictor),8)
          meanAbsPredictors(iPredictor) = meanAbsPredictors(iPredictor)/real(countObs(iPredictor),8)
          write(*,*) 'bias_calcMeanPredictors: mean(abs),mean,min,max=',iPredictor, &
            meanAbsPredictors(iPredictor),meanPredictors(iPredictor),minPredictors(iPredictor),maxPredictors(iPredictor)
        end if
      end do
    end if

  END SUBROUTINE bias_calcMeanPredictors

  !---------------------------------------------
  ! bias_calcBias_ad
  !---------------------------------------------
  SUBROUTINE bias_calcBias_ad(cv_out,cv_dim,obsColumnIndex,obsSpaceData)
    implicit none

    integer  :: cv_dim
    real(8)  :: cv_out(cv_dim)
    integer  :: obsColumnIndex
    type(struct_obs)  :: obsSpaceData

    integer  :: index_header,index_body,iobs, idatyp
    integer  :: iSensor,iChannel,iPredictor,index_cv,nsize,ierr
    integer  :: iScan, iFOV, jPred
    real(8)  :: predictor(NumPredictors)
    real(8),pointer  :: cv_bias(:)
    real(8)  :: biasCor

    if ( .not. lvarbc ) return

    if ( mpi_myid == 0 ) write(*,*) 'Starting bias_calcBias_ad'

    if ( mpi_myid == 0 ) then
      if ( cvm_subVectorExists('BIAS') ) then
        cv_bias => cvm_getSubVector(cv_out,'BIAS')
      else
        write(*,*) 'bias_calcBias_ad: control vector does not include bias coefficients'
        return
      end if
    end if

    ! adjoint of applying bias increment to specified obs column
    do iSensor = 1, tvs_nSensors
      bias(iSensor)%coeffIncr(:,:) = 0.0d0
      bias(iSensor)%coeffIncr_fov(:,:) = 0.0d0
    end do

    iobs = 0
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      index_header = obs_getHeaderIndex(obsSpaceData)
      if ( index_header < 0 ) exit HEADER
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,index_header)
      if ( .not.  tvs_isIdBurpTovs(idatyp) ) cycle HEADER  

      iobs = iobs + 1
      iSensor = tvs_lsensor(tvs_ltovsno(index_header))
      call bias_getPredictors(predictor,index_header,iobs,obsSpaceData)
      call obs_set_current_body_list(obsSpaceData, index_header)
      iFov = obs_headElem_i(obsSpaceData,OBS_FOV,index_header)

      BODY: do
        index_body = obs_getBodyIndex(obsSpaceData)
        if ( index_body < 0 ) exit BODY

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) /= obs_assimilated ) cycle BODY  
        iChannel = nint(obs_bodyElem_r(obsSpaceData,OBS_PPP,index_body))
        iChannel = max(0,min(iChannel,tvs_maxChannelNumber+1))
        iChannel = iChannel-tvs_channelOffset(iSensor)

        biasCor = obs_bodyElem_r(obsSpaceData,obsColumnIndex,index_body)
       
        do iPredictor = 1, bias(iSensor)%numActivePredictors
          jPred = bias(iSensor)%PredictorIndex(iPredictor)

          if ( iPredictor == 1 ) then
            if ( bias(iSensor)%numScan > 1) then
              iScan = iFov
            else
              iScan = 1
            end if
            bias(iSensor)%coeffIncr_fov(iChannel,iScan) = bias(iSensor)%coeffIncr_fov(iChannel,iScan) &
                 + predictor(jPred)*biasCor
          else
            bias(iSensor)%coeffIncr(iChannel,iPredictor) = bias(iSensor)%coeffIncr(iChannel,iPredictor) & 
               + predictor(jPred)*biasCor
          end if

        end do !iPredictor
      end do BODY

    end do HEADER

    ! put the coefficients into the control vector
    call bias_cvToCoeff_ad(cv_bias)

    if ( mpi_myid == 0 ) then
      write(*,*) 'bias_calcBias_ad: maxval(cv_bias)=',maxval(cv_bias(:))
    end if

  END SUBROUTINE bias_calcBias_ad

  !----------------------------------------------------
  ! bias_cvToCoeff_ad
  !----------------------------------------------------
  SUBROUTINE bias_cvToCoeff_ad(cv_bias)
    implicit none

    real(8)  :: cv_bias(:)
    integer  :: index_cv,iSensor,iChannel,iPredictor,iScan
    integer  :: nChan,nPre,nScan
    integer  :: nsize,ierr
    real(8),allocatable  :: temp_coeffIncr(:,:), temp_coeffIncr_fov(:,:)

    do iSensor = 1, tvs_nSensors
      nChan = bias(iSensor)%numChannels
      nPre  = bias(iSensor)%numActivePredictors
      allocate(temp_coeffIncr(1:nChan, 2:nPre)) 
      
      temp_coeffIncr(:,:) = 0.0d0
      nsize = bias(iSensor)%numChannels * (bias(iSensor)%numActivePredictors-1)
      call rpn_comm_reduce(bias(iSensor)%coeffIncr,temp_coeffIncr,nsize,"mpi_double_precision","mpi_sum",0,"GRID",ierr)
      bias(iSensor)%coeffIncr(:,:) = temp_coeffIncr(:,:)
      deallocate(temp_coeffIncr)
    end do

    do iSensor = 1, tvs_nSensors
      nChan = bias(iSensor)%numChannels
      nScan  = bias(iSensor)%numScan
      allocate(temp_coeffIncr_fov(1:nChan, 1:nScan)) 
      temp_coeffIncr_fov(:,:) = 0.0d0
      nsize = bias(iSensor)%numChannels * bias(iSensor)%numScan
      call rpn_comm_reduce(bias(iSensor)%coeffIncr_fov,temp_coeffIncr_fov,nsize,"mpi_double_precision","mpi_sum",0,"GRID",ierr)
      bias(iSensor)%coeffIncr_fov(:,:) = temp_coeffIncr_fov(:,:)
      deallocate(temp_coeffIncr_fov)
    end do

    if ( mpi_myid == 0 ) then
      index_cv = 0
      do iSensor = 1, tvs_nSensors
        do iChannel = 1, bias(iSensor)%numChannels
          do iPredictor = 1, bias(iSensor)%numActivePredictors
            if ( iPredictor == 1 ) then
              do iScan = 1, bias(iSensor)%numScan
                index_cv  =  index_cv  +  1
                cv_bias(index_cv) = bias(iSensor)%stddev(iChannel,iPredictor) * bias(iSensor)%coeffIncr_fov(iChannel,iScan)
              end do
            else
              index_cv = index_cv + 1
              cv_bias(index_cv) = bias(iSensor)%stddev(iChannel,iPredictor) * bias(iSensor)%coeffIncr(iChannel,iPredictor)
            end if
          end do
        end do
      end do
    end if
    
  END SUBROUTINE bias_cvToCoeff_ad

  !-----------------------------------------
  ! bias_writeBias
  !-----------------------------------------
  SUBROUTINE bias_writeBias(cv_in,cv_dim)
    implicit none

    integer  :: cv_dim
    real(8)  :: cv_in(cv_dim)

    integer  :: iSensor,iChannel,iPredictor
    integer  :: jSensor,jChannel
    integer  :: fnom,fclos,nulfile_inc,nulfile_fov,ierr
    real(8),pointer   :: cv_bias(:)
    character(len=80) :: BgFileName
    real(8)           :: biasCoeff_bg(tvs_nSensors,maxNumChannels,NumPredictors)
    logical           :: fileExists

    !for background coeff and write out
    integer             :: iInstr, iuncoef
    real(8)             :: fovbias_bg(tvs_nSensors,maxNumChannels,maxfov)
    integer             :: numCoefFile,jCoef,kCoef
    character(len=10)   :: coefInstrName(tvs_nSensors), temp_instrName, instrName
    character(len=25)   :: filecoeff
    logical             :: coeffExists
    ! these variables are not used
    character(len=10)  :: sats(tvs_nSensors)
    integer            :: chans(tvs_nSensors, tvs_maxChannelNumber), nsat, nfov
    integer            :: nchan(tvs_nSensors)
    character(len=6)   :: cinstrum

    if ( .not. lvarbc ) return

    if ( mpi_myid == 0 ) then
      if ( cvm_subVectorExists('BIAS') ) then
        cv_bias => cvm_getSubVector(cv_in,'BIAS')
        write(*,*) 'bias_writeBias: maxval(cv_bias)=',maxval(cv_bias(:))
      else
        write(*,*) 'bias_writeBias: control vector does not include bias coefficients'
        return
      end if
    end if

    call bias_cvToCoeff(cv_bias)

    ! write out bias coefficient increments in ascii file
    nulfile_inc = 0
    ierr = fnom(nulfile_inc,'./satbias_increment.dat','FTN+FMT',0)

    do iSensor = 1, tvs_nSensors
      write(nulfile_inc,'(/,1X,"Sensor Index=",I3,", Satellite Name=",A15,", Instrument Name=",A15)') &
             iSensor,tvs_satelliteName(iSensor),tvs_instrumentName(iSensor)
      do iChannel = 1, bias(iSensor)%numChannels
        if ( sum(bias(iSensor)%coeffIncr(iChannel,:)) /= 0.0d0 ) &
          write(nulfile_inc,'(3X,"Channel number=",I4)') iChannel
        do iPredictor = 2, bias(iSensor)%numActivePredictors
          if ( bias(iSensor)%coeffIncr(iChannel,iPredictor) /= 0.0d0 ) &
            write(nulfile_inc,'(5X,"Predictor number=",I4,", Coefficient=",e12.4)') &
                 iPredictor,bias(iSensor)%coeffIncr(iChannel,iPredictor)
        end do
      end do
    end do

    ierr = fclos(nulfile_inc)

    ! write out fovbias coefficient increments in ascii file
    nulfile_fov = 0
    ierr = fnom(nulfile_fov,'./fovbias_incre.dat','FTN+FMT',0)
    do iSensor = 1, tvs_nSensors
      write(nulfile_fov,'(/,1X,"Sensor Index=",I3,", Satellite Name=",A15,", Instrument Name=",A15)') &
                       iSensor,tvs_satelliteName(iSensor),tvs_instrumentName(iSensor)
      do iChannel = 1, bias(iSensor)%numChannels
        if ( sum(bias(iSensor)%coeffIncr_fov(iChannel,:)) /= 0.0d0 ) &
          write(nulfile_fov,'(3X,"Channel number=",I4)') iChannel 
        if ( sum(bias(iSensor)%coeffIncr_fov(iChannel,:)) /= 0.0d0 ) & 
          write(nulfile_fov,*) bias(iSensor)%coeffIncr_fov(iChannel,:)
      end do
    end do
    ierr = fclos(nulfile_fov)

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
    if ( mpi_myid == 0 ) then
      do iInstr=1, numCoefFile 

        biasCoeff_bg(:,:,:) = 0.0
        fovbias_bg(:,:,:) = 0.0
        BgFileName ='./coeff_file_'//coefInstrName(iInstr)//''
        call bias_updateCoeff(tvs_nSensors,NumPredictors,BgFileName,sats,chans,nsat,nchan,nfov,cinstrum)

      end do
    end if

  END SUBROUTINE bias_writeBias

  !--------------------------------------
  ! bias_updateCoeff
  !--------------------------------------
  SUBROUTINE bias_updateCoeff(maxsat,maxpred,coeff_file,sats,chans,nsat,nchan,nfov,cinstrum,updateCoeff_opt)
    implicit none

    ! There are three parts in this subroutine, read, update and write out the coeff files
    ! IN
    integer            :: maxsat, maxpred
    character(len=80)  :: coeff_file
    logical,optional   :: updateCoeff_opt

    ! OUT 
    character(len=10)  :: sats(maxsat)        ! dim(maxsat), satellite names
    integer            :: chans(maxsat, maxNumChannels)       ! dim(maxsat, maxchan), channel numbers
    integer            :: nsat, nfov
    integer            :: nchan(maxsat)       ! dim(maxsat), number of channels
    character(len=6)   :: cinstrum    ! string: instrument (e.g. AMSUB)
 
    ! Local
    real(8)            :: fovbias(maxsat,maxNumChannels,maxfov)     ! dim(maxsat,maxchan,maxfov), bias as F(fov)
    real(8)            :: coeff(maxsat,maxNumChannels,maxpred)       ! dim(maxsat,maxchan,maxpred+1)
    character(len=2)   :: ptypes(maxsat,maxNumChannels,maxpred) ! dim(maxsat,maxNumChannelsmaxchan,maxpred)
    integer            :: npred(maxsat, maxNumChannels)       ! dim(maxsat, maxchan), number of predictors

    ! LOCAL for reading background coeff file
    character(len=10)  :: sat
    character(len=120) :: line
    integer            :: chan, ndata, nbfov, nbpred
    integer            :: ier, istat 
    integer            :: iSat, jChan, kPred, kFov, totSat
    logical            :: newsat, verbose
    real               :: dummy
    integer            :: iun
    integer            :: fnom, fclos

    ! update coeff files
    real               :: fovbias_an(maxsat,maxNumChannels,maxfov)
    real               :: coeff_an(maxsat,maxNumChannels,maxpred) 
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

    ! 
    !- 1. read in the background coeff files, this program is read_coeff from genbiascorr
    ! 
    if ( present(updateCoeff_opt) ) then
      updateCoeff_opt2 = updateCoeff_opt
    else
      updateCoeff_opt2 = .true.
    end if

    verbose = .false.
    coeff    = 0.0
    fovbias  = 0.0
    sats     = 'XXXXXXXX'
    cinstrum = 'XXXXXX'
    chans    = 0
    npred    = 0
    nsat     = 0
    nchan    = 0
    nfov     = 0
    ptypes   = 'XX'

    write(*,*) 'bias_updateCoeff: Reading coeff file starts'  
    iun = 0
    ier = fnom(iun,coeff_file,'FMT',0)
    if ( ier /= 0 ) then
      write(*,*) 'bias_updateCoeff: ERROR - Problem opening the coeff file!', coeff_file 
      call utl_abort('bias_updateCoeff for read_coeff')
    end if

    write(*,*) 'bias_updateCoeff: open background file open = ', coeff_file
    read(iun,*,iostat=istat)
    if ( istat < 0 ) then
      write(*,*) 'bias_updateCoeff: File appears empty', coeff_file
      return
    end if
    rewind(iun)

    totSat = 0
    ! Loop over the satellites/channels in the file
    do
      read(iun,'(A)',iostat=istat) line
      if ( istat < 0 ) exit
      if ( line(1:3) == 'SAT' ) then
        newsat = .true.
        read(line,'(T53,A8,1X,A6,1X,I6,1X,I8,1X,I2,1X,I3)',iostat=istat) sat, cinstrum, chan, ndata, nbpred, nbfov
        do iSat = 1, maxsat
          if ( trim(sats(iSat)) == trim(sat) ) then
            newsat = .false.
            totSat = iSat
          end if
        end do
        if ( newsat ) then
          totSat = totSat + 1
          sats(totSat) = sat
          if ( totSat > 1 ) nchan(totSat-1) = jChan
          jChan = 1
        else
          jChan = jChan + 1
        end if
        chans(totSat, jChan) = chan
        npred(totSat, jChan) = nbpred
        read(iun,'(A)',iostat=istat) line
        if ( nbpred > 0 ) then
          read(line,'(T8,6(1X,A2))',iostat=istat) (ptypes(totSat,jChan,kPred),kPred=1,nbpred)
        end if
        read(iun,*,iostat=istat) (fovbias(totSat,jChan,kFov),kFov=1,nbfov)
        if ( nbpred > 0 ) then
          read(iun,*,iostat=istat) (coeff(totSat,jChan,kPred),kPred=1,nbpred+1)
        else
          read(iun,*,iostat=istat) dummy
        end if
      end if
    end do

    if ( totSat == 0 ) then
      write(*,*) 'bias_updateCoeff: No data read from coeff file!', coeff_file
      call utl_abort('bias_updateCoeff')
    end if
    nsat      = totSat 
    nfov      = nbfov
    nchan(totSat) = jChan

    if ( verbose ) then
      write(*,*) ' '
      write(*,*) ' ------------- BIAS CORRECTION COEFFICIENT FILE ------------------ '
      write(*,*) ' '
      write(*,*) ' Number of satellites =     ', nsat
      write(*,*) ' Number of FOV =            ', nfov
      write(*,*) ' Max number of predictors = ', maxval(npred)
      write(*,*) ' '
      do iSat = 1, nsat
        write(*,*) '  Satellite = ' // sats(iSat)
        write(*,*) '     Number of channels = ', nchan(iSat)
        write(*,*) '     predictors, fovbias, coeff for each channel: '
        do jChan = 1, nchan(iSat)
          write(*,*) iSat, chans(iSat,jChan)
          if ( npred(iSat,jChan) > 0 ) then
            write(*,'(6(1X,A2))') (ptypes(iSat,jChan,kPred),kPred=1,npred(iSat,jChan))
          else
            write(*,'(A)') 'No predictors'
          end if
          write(*,*) (fovbias(iSat,jChan,kFov),kFov=1,nfov)
          write(*,*) (coeff(iSat,jChan,kPred),kPred=1,npred(iSat,jChan)+1)
        end do
      end do
        write(*,*) ' '
    end if
    ier = fclos(iun)

    write(*,*) 'bias_updateCoef: Reading coeff file done', coeff_file


    if ( updateCoeff_opt2 == .false. ) return 

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

            if ( chans(iSat,jChan) /= jChannel ) cycle 

            ! part 1 for coeffIncr
            do iFov = 1, nfov
              if ( bias(iSensor)%coeffIncr_fov(jChannel,iFov) /= 0.0d0 ) then
                fovbias_an(iSat,jChan,iFov) = fovbias(iSat,jChan,iFov) + bias(iSensor)%coeffIncr_fov(jChannel,iFov)
              end if
            end do ! iFov

            ! part 2 for coeffIncr_fov
            totPred  = bias(iSensor)%NumActivePredictors 
            do iPred = 1, totPred
              if ( iPred == 1 ) then
                coeff_an(iSat,jChan,iPred) = coeff(iSat,jChan,iPred)
              else
                if ( bias(iSensor)%coeffIncr(jChannel,iPred) /= 0.0d0 ) then
                  coeff_an(iSat,jChan,iPred) = coeff(iSat,jChan,iPred) + bias(iSensor)%coeffIncr(jChannel,iPred)
                end if
              end if
            end do ! iPred

          end do ! jChannel
        end do !jChan
      end do !iSensor
    end do ! iSat

    !
    !- 3. Write out updated_coeff
    ! 
    iuncoef2 = 0
    filename2 ='./anlcoeffs_'//cinstrum//'' 
    ierr = fnom(iuncoef2, filename2,'FTN+FMT',0)

    write(*,*) 'bias_updateCoeff: write in bias_updateCoeff'
   
    do iSat = 1, nsat
      do jChan = 1, nchan(iSat)      
        numPred = npred(iSat,jChan)
        write(iuncoef2, '(A52,A8,1X,A6,1X,I6,1X,I8,1X,I2,1X,I3)') &
          'SATELLITE, INSTRUMENT, CHANNEL, NOBS, NPRED, NSCAN: ', sats(iSat), cinstrum, chans(iSat,jChan), ndata, numPred, nfov
        write(iuncoef2, '(A7,6(1X,A2))') 'PTYPES:', (ptypes(iSat,jChan,kPred) , kPred=1,numPred)
        write(iuncoef2,'(120(1x,ES17.10))') (fovbias_an(iSat,jChan,kFov),kFov=1,nfov)
        write(iuncoef2,*) (coeff_an(iSat,jChan,kPred),kPred=1,numPred+1)
      end do
    end do

    ierr = fclos(iuncoef2) 
   
    write(*,*) 'bias_updateCoeff: finish writing coeffient file', filename2
    
  END SUBROUTINE bias_updateCoeff

  !----------------------
  ! bias_Finalize
  !----------------------
  SUBROUTINE bias_Finalize
    implicit none

    integer    :: iSensor

    if ( .not.lvarbc ) return

    deallocate(trialHeight300m1000)
    deallocate(trialHeight50m200)
    deallocate(trialHeight1m10)
    deallocate(trialHeight5m50)
    deallocate(trialTG)

    do iSensor = 1, tvs_nSensors
      deallocate(bias(iSensor)%stddev)
      deallocate(bias(iSensor)%coeffIncr)
      deallocate(bias(iSensor)%predictorIndex)
      deallocate(bias(iSensor)%coeffIncr_fov)
    end do

  END SUBROUTINE bias_Finalize 

  !-----------------------------
  ! Lower
  !-----------------------------
  function Lower(s1) result(s2) 
    implicit none
 
    character(*)        :: s1
    character(len(s1))  :: s2
    character           :: ch 
    integer, parameter  :: duc = ichar('A') -ichar('a')
    integer             :: iStr 

    do iStr = 1, len(s1)
      ch = s1(iStr:iStr)
      if (ch >= 'A' .and. ch <= 'Z') ch = char(ichar(ch)-duc)
      s2(iStr:iStr) =  ch
    end do
  
  end function Lower
 
  !-----------------------------
  ! InstrNametoCoeffFileName 
  !-----------------------------
  function InstrNametoCoeffFileName(nameIn) result(nameOut)
    implicit none

    character(len=10)  :: nameIn, temp_instrName
    character(len=10)  :: nameOut
  
    temp_instrName = Lower(nameIn)
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
    
    character(len=10)  :: nameIn
    character(len=10)  :: nameOut

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
    
    character(len=10)  :: nameIn
    character(len=10)  :: nameOut

    if ( trim(nameIn) == 'MSG2' ) then
      nameOut = 'METSAT9'
    else if ( trim(nameIn) == 'MSG3' ) then
      nameOut = 'METSAT10' 
    else if ( trim(nameIn) == 'METEOSAT7' ) then
      nameOut = 'METSAT7' 
    else 
      nameOut = nameIn
    end if

  end function SatNameinCoeffFile

END MODULE biasCorrection_mod

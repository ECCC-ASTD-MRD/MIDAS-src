!--------------------------------------------------------------------------
! MODULE biascorrection (Variational bias correction.  prefix="bias")
!
! Purpose: 
!
! Subroutines (public):
!    bias_setup
!
! Dependencies:
!
!--------------------------------------------------------------------------
module biascorrection_mod
  use utilities_mod
  use ramDisk_mod
  use MathPhysConstants_mod
  use obsSpaceData_mod
  use controlVector_mod
  use mpivar_mod
  use tovs_nl_mod
  use timeCoord_mod

  implicit none
  save
  private

  public               :: bias_setup,bias_calcBias_tl,bias_calcBias_ad, bias_writeBias, bias_Finalize

  type  :: struct_bias
   integer  :: numActivePredictors, numChannels, numScan
   real(8),allocatable  :: stddev(:,:)
   real(8),allocatable  :: coeffIncr(:,:)
   integer,allocatable  :: predictorIndex(:)
   real(8),allocatable  :: coeffIncr_fov(:,:)
   integer,allocatable  :: ChannelNum(:)
  end type struct_bias

  type(struct_bias),allocatable  :: bias(:)
   
  logical               :: initialized = .false.
  integer               :: cvdim
  integer, parameter    :: maxNumChannels = tvs_maxChannelNumber
  integer, parameter    :: NumPredictors = 11
  integer, parameter    :: maxfov = 120
  real(8), allocatable  :: trialGZ300m1000(:)
  real(8), allocatable  :: trialGZ50m200(:)
  real(8), allocatable  :: trialGZ1m10(:)
  real(8), allocatable  :: trialGZ5m50(:)
  real(8), allocatable  :: trialTG(:)
  integer               :: nobs

  logical  :: lvarbc
  real(8)  :: bg_stddev(NumPredictors)
  logical  :: abortIfNoBackground
  namelist /nambias/lvarbc,bg_stddev,abortIfNoBackground

CONTAINS

  subroutine bias_setup(cvdim_out)
    implicit none

    integer  :: cvdim_out
    integer  :: iSensor,iPredictor
    integer  :: jPredictor
    integer  :: ierr,nulnam
    integer  :: fnom,fclos
    integer  :: jj 
    integer  :: idNum(tvs_nSensors,NumPredictors)
    integer  :: activePredictors(tvs_nSensors)
    character(len=85)  :: filecoeff
    character(len=10)  :: instrName, instrNamecoeff, satNamecoeff 
    logical            :: coeffExists 

    !variables from background coeff file
    character(len=10)  :: sats(tvs_nSensors)
    integer            :: chans(tvs_nSensors, tvs_maxChannelNumber), nsat, nfov, ii
    integer            :: nchan(tvs_nSensors)
    character(len=6)   :: cinstrum

    ! set default values for namelist variables
    lvarbc   = .false.
    bg_stddev(:) = 0.0d0
    abortIfNoBackground = .true.

    ! read in the namelist NAMBIAS
    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=nambias,iostat=ierr)
    if(ierr.ne.0) call utl_abort('bias_setup: Error reading namelist')
    write(*,nml=nambias)
    ierr=fclos(nulnam)

    cvdim = 0

    if(lvarbc) then
      allocate(bias(tvs_nSensors))

      do iSensor = 1, tvs_nSensors

        bias(iSensor)%numActivePredictors = 0
        activePredictors(iSensor) = 0
        jj = 0

        ! set predictors, numScan and numChannels
        iPredictor = 1  ! constant predictor for all sensors
        activePredictors(iSensor)=ibset(activePredictors(iSensor),iPredictor-1)
        iPredictor=2  ! GZ300-GZ1000
        activePredictors(iSensor)=ibset(activePredictors(iSensor),iPredictor-1)
        iPredictor=3  ! GZ50-GZ200
        activePredictors(iSensor)=ibset(activePredictors(iSensor),iPredictor-1)
        iPredictor=7  ! GZ5-GZ50
        activePredictors(iSensor)=ibset(activePredictors(iSensor),iPredictor-1)

        do iPredictor=1,NumPredictors
          if( btest(activePredictors(iSensor),iPredictor-1) ) then
            bias(iSensor)%numActivePredictors = bias(iSensor)%numActivePredictors + 1
            jj = bias(iSensor)%numActivePredictors
            idNum(iSensor,jj) = iPredictor
          end if
        end do

        instrName = InstrNametoCoeffFileName(tvs_instrumentName(iSensor))
        instrNamecoeff = InstrNameinCoeffFile(tvs_instrumentName(iSensor))
        satNamecoeff = SatNameinCoeffFile(tvs_satelliteName(iSensor))

        filecoeff = 'coeff_file_'//trim(instrName)//''
        inquire(file=trim(filecoeff),exist = coeffExists)
         
        if( coeffExists ) then
          call read_coeff(tvs_nSensors,NumPredictors,filecoeff,sats,chans,nsat,nchan,nfov,cinstrum)
          do ii = 1, nsat
            if( sats(ii) == trim(satNamecoeff) .and. cinstrum == trim(instrNamecoeff) ) then
              bias(iSensor)%numScan = nfov
              allocate( bias(iSensor)%ChannelNum(nchan(ii)) ) 
              bias(iSensor)%ChannelNum(:)  = chans(ii,:)
              bias(iSensor)%numChannels = maxval(bias(iSensor)%ChannelNum(:))
            end if
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
        if( tvs_instruments(iSensor) == 3 ) then
          if( bias(iSensor)%numChannels /= 0 ) then
            bias(iSensor)%stddev(13:14, :) = 0.0d0
          end if
        end if
      end do
    end if

    if( mpi_myid == 0) then
      cvdim_out = cvdim
    else
      cvdim_out = 0
    end if

  end subroutine bias_setup

  !---------------------------------------
  ! calcBias_tl
  !---------------------------------------- 
  subroutine bias_calcBias_tl(cv_in,cv_dim,obsColumnIndex,multFactor,lobsSpaceData)
    implicit none

    integer  :: cv_dim
    real(8)  :: cv_in(cv_dim)
    integer  :: obsColumnIndex
    real(8)  :: multFactor
    type(struct_obs)  :: lobsSpaceData

    integer  :: index_header,index_body,iobs, indxtovs, idatyp
    integer  :: iSensor,iChannel,iPredictor,index_cv
    integer  :: iScan, iFov, jj, iiScan
    real(8)  :: predictor(NumPredictors)
    real(8),pointer  :: cv_bias(:)
    real(8)  :: biasCor
    logical,save  :: firstTime=.true.

    if(.not.lvarbc) return

    if( firstTime ) then
      call bias_getTrialPredictors(lobsSpaceData)
      call bias_calcMeanPredictors(lobsSpaceData)
      firstTime=.false.
    end if

    if( mpi_myid == 0) then
      if( cvm_subVectorExists(cvm_Bias) ) then
        cv_Bias => cvm_getSubVector(cv_in,cvm_Bias)
        write(*,*) 'bias_calcBias_tl: maxval(cv_bias)=',maxval(cv_bias(:))
        write(*,*) 'bias_calcBias_tl: cv_bias(69)=',cv_bias(69)
      else
        write(*,*) 'bias_calcBias_tl: control vector does not include bias coefficients'
        return
      end if
    end if

    ! get bias coefficients
    call bias_cvToCoeff(cv_bias)

    ! apply bias increment to specified obs column
    iobs = 0
    call obs_set_current_header_list(lobsSpaceData,'TO')
    HEADER: do
      index_header = obs_getHeaderIndex(lobsSpaceData)
      if (index_header < 0) exit HEADER

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(lobsSpaceData,OBS_ITY,index_header)
      if ( .not. tvs_isIdBurpTovs(idatyp) ) cycle HEADER
      
      indxtovs = tvs_ltovsno(index_header)
      if ( indxtovs == 0 ) then
        call utl_abort('bias_calcBias_tl')
      end if

      iobs = iobs + 1
      iSensor = tvs_lsensor(tvs_ltovsno(index_header))
      call bias_getPredictors(predictor,index_header,iobs,lobsSpaceData)

      call obs_set_current_body_list(lobsSpaceData, index_header)
      iFov = obs_headElem_i(lobsSpaceData,OBS_FOV,index_header)

      BODY: do
        index_body = obs_getBodyIndex(lobsSpaceData)
        if( index_body < 0 ) exit BODY

        if( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) .ne. 1 ) cycle BODY   

        iChannel = nint(obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body))
        iChannel = max(0,min(iChannel,tvs_maxChannelNumber+1))
        iChannel = iChannel-tvs_channelOffset(iSensor)
        biasCor = 0.0d0

        do iPredictor = 1, bias(iSensor)%NumActivePredictors
          jj = bias(iSensor)%PredictorIndex(iPredictor)
          if( iPredictor == 1 ) then

            if (bias(iSensor)%numScan .gt. 1) then
              iScan = iFov
            else
              iScan = 1
            end if
            biasCor = biasCor + predictor(jj) * bias(iSensor)%coeffIncr_fov(iChannel,iScan) 
          else
            biasCor = biasCor + predictor(jj) * bias(iSensor)%coeffIncr(iChannel,iPredictor) 
          end if
          
        end do
        call obs_bodySet_r( lobsSpaceData, obsColumnIndex, index_body, &
          obs_bodyElem_r(lobsSpaceData,obsColumnIndex,index_body) - &
          multFactor*biasCor )
      end do BODY
    end do HEADER

  end subroutine bias_calcBias_tl

  !----------------------
  ! getTrialPredictors
  !----------------------
  subroutine bias_getTrialPredictors(lobsSpaceData)
    implicit none

    type(struct_obs)  :: lobsSpaceData
    character(len=80)  :: trialfilename="./trlp"
    integer  :: iSensor,iPredictor
    integer  :: ierr,nulfst,iset
    integer  :: fnom,fclos,fstouv,fstfrm,vezgdef,ezdefset
    integer  :: fstinf,fstprm,fstlir,ezqkdef,ezsint
    integer  :: index_header, idatyp
    integer  :: obsgid,trlgid
    real(8)  :: lat_r8,lon_r8
    real(8), allocatable  ::  dlonfld(:), dlatfld(:)
    real(8)  :: zig1,zig2,zig3,zig4
    integer  :: ig1obs,ig2obs,ig3obs,ig4obs
    character(len=1)  :: clgrtyp
    character(len=2)  :: cltypvar
    character(len=4)  :: varName
    character(len=12) :: cletiket
    integer  :: key,dateo,ip1,ip2,ip3,ni,nj,nk
    integer  :: deet,npas,nbits,datyp,ig1,ig2,ig3,ig4,swa,lng,dltf,ubc,extra1,extra2,extra3
    real(8), allocatable  :: varTrial1(:,:),varTrial2(:,:)
    real(8), allocatable  :: varTrial3(:,:)

    integer           :: timeIndex
    character(len=80) :: trialfile
    character(len=2)  :: flnum
    integer  :: nultrlm
    logical  :: trialExists, FileExists

    ! count number of tovs locations
    nobs = 0
    call obs_set_current_header_list(lobsSpaceData,'TO')
    HEADER: do
      index_header = obs_getHeaderIndex(lobsSpaceData)
      if (index_header < 0) exit HEADER
      idatyp = obs_headElem_i(lobsSpaceData,OBS_ITY,index_header)
      if ( .not.  tvs_isIdBurpTovs(idatyp) ) cycle HEADER 
      nobs = nobs + 1
    end do HEADER
    allocate(dlonfld(nobs),dlatfld(nobs))

    ! define "Y" grid of observations for ezscint
    nobs = 0
    call obs_set_current_header_list(lobsSpaceData,'TO')
    HEADER2: do
      index_header = obs_getHeaderIndex(lobsSpaceData)
      if (index_header < 0) exit HEADER2
      idatyp = obs_headElem_i(lobsSpaceData,OBS_ITY,index_header)
      if ( .not.  tvs_isIdBurpTovs(idatyp) ) cycle HEADER2

      nobs = nobs + 1

      lat_r8 = obs_headElem_r(lobsSpaceData,OBS_LAT,index_header)
      lon_r8 = obs_headElem_r(lobsSpaceData,OBS_LON,index_header)
      if( lon_r8 < 0.0d0 ) lon_r8 = lon_r8 + 2*MPC_PI_R8
      if( lon_r8 >= 2.0d0*MPC_PI_R8) lon_r8 = lon_r8 - 2*MPC_PI_R8
      dlatfld(nobs) = lat_r8*MPC_DEGREES_PER_RADIAN_R8
      dlonfld(nobs) = lon_r8*MPC_DEGREES_PER_RADIAN_R8

    end do HEADER2

    zig1 = 0.0D0
    zig2 = 0.0D0
    zig3 = 1.0D0
    zig4 = 1.0D0
    call utl_cxgaig('L',ig1obs,ig2obs,ig3obs,ig4obs,zig1,zig2,zig3,zig4)

    if ( nobs > 0 ) then
      obsgid = utl_ezgdef(nobs,1,'Y','L',ig1obs,ig2obs,ig3obs,ig4obs,dlonfld(:),dlatfld(:))
    else
      write(*,*) 'bias_getTrialPredictors: NO OBS found'
      obsgid = -999
    end if

    ! read trial fields
    nulfst = 0
    trialfilename = './trlp'

    inquire(file=trim(trialfilename), exist=FileExists)
    write(*,*) 'bias_getTrialfile', FileExists
    ierr = fnom(nulfst,trialfilename,'RND+OLD+R/O',0)
    if ( ierr /= 0 ) then
      write(*,*) 'bias_getTrialPredictors: trialfilename=',trialfilename
      call utl_abort('bias_getTrialPredictors: Error opening trial file')
    end if
    ierr = fstouv(nulfst,'RND+OLD')
    write(*,*) 'bias_getTrialPredictors: FileName=', trialfilename 
    write(*,*) 'bias_getTrialPredictors: opened as unit # ',nulfst ! Determine grid size and EZSCINT ID for GZ
    dateo  = -1
    cletiket = ' '
    ip1    = -1
    ip2    = -1
    ip3    = -1
    cltypvar = ' '
    varName = 'GZ'
    key = fstinf( nulfst,                                    & ! IN
                  ni,nj,nk,                                  & ! OUT
                  dateo,cletiket,ip1,ip2,ip3,cltypvar,varName )! IN
    if (key < 0) then
      write(*,*) 'bias_getTrialPredictors: Unable to find trial field = ',varName
      call utl_abort('bias_getTrialPredictors')
    end if
    ierr = fstprm( key,                                               & ! IN
                   dateo, deet, npas, ni, nj, nk, nbits,              & ! OUT
                   datyp, ip1, ip2, ip3, cltypvar, varName, cletiket, & ! OUT
                   clgrtyp, ig1, ig2, ig3,                            & ! OUT
                   ig4, swa, lng, dltf, ubc, extra1, extra2, extra3 )   ! OUT
    trlgid = ezqkdef(ni,nj,clgrtyp,ig1,ig2,ig3,ig4,nulfst)
    write(*,*) 'bias_getTrialPredictors: trlgid=',trlgid,ni,nj

    iset = ezdefset(obsgid,trlgid)

    allocate(varTrial1(ni,nj))
    allocate(varTrial2(ni,nj))
    allocate(trialGZ300m1000(nobs))
    allocate(trialGZ50m200(nobs))
    allocate(trialGZ5m50(nobs))
    allocate(trialGZ1m10(nobs))
    allocate(trialTG(nobs))

    ! GZ300 - GZ1000
    varName= 'GZ'
    ip1 = 300
    key = utl_fstlir(varTrial1,nulfst,ni,nj,nk,-1,' ',ip1,-1,-1,' ',varName)
    if ( key < 0 ) then
      write(*,*) 'bias_getTrialPredictors: Unable to find trial field - varname,ip1= ',varName,ip1
      call utl_abort('bias_getTrialPredictors')
    end if
    ip1 = 1000
    key = utl_fstlir(varTrial2,nulfst,ni,nj,nk,-1,' ',ip1,-1,-1,' ',varName)
    if (key < 0) then
      write(*,*) 'bias_getTrialPredictors: Unable to find trial field - varname,ip1= ',varName,ip1
      call utl_abort('bias_getTrialPredictors')
    end if
    varTrial1(:,:) = varTrial1(:,:) - varTrial2(:,:)
    ierr = utl_ezsint(trialGZ300m1000,varTrial1,nobs,1,1,ni,nj,1)
    if( ierr < 0 ) call utl_abort('bias_getTrialPredictors: error interpolating GZ300-GZ1000')

    ! GZ50 - GZ200
    varName = 'GZ'
    ip1 = 50
    key = utl_fstlir(varTrial1,nulfst,ni,nj,nk,-1,' ',ip1,-1,-1,' ',varName)
    if (key < 0) then
      write(*,*) 'bias_getTrialPredictors: Unable to find trial field - varname,ip1= ',varName,ip1
      call utl_abort('bias_getTrialPredictors')
    end if
    ip1 = 200
    key = utl_fstlir(varTrial2,nulfst,ni,nj,nk,-1,' ',ip1,-1,-1,' ',varName)
    if (key < 0) then
      write(*,*) 'bias_getTrialPredictors: Unable to find trial field - varname,ip1= ',varName,ip1
      call utl_abort('bias_getTrialPredictors')
    end if
    varTrial1(:,:) = varTrial1(:,:) - varTrial2(:,:)
    ierr = utl_ezsint(trialGZ50m200,varTrial1,nobs,1,1,ni,nj,1)
    if( ierr < 0 ) call utl_abort('bias_getTrialPredictors: error interpolating GZ50-GZ200')

    ! GZ5 - GZ50
    varName = 'GZ'
    ip1 = 1900 ! 5hPa
    key = utl_fstlir(varTrial1,nulfst,ni,nj,nk,-1,' ',ip1,-1,-1,' ',varName)
    if (key < 0) then
      write(*,*) 'bias_getTrialPredictors: Unable to find trial field - varname,ip1= ',varName,ip1
      call utl_abort('bias_getTrialPredictors')
    end if
    ip1 = 50
    key = utl_fstlir(varTrial2,nulfst,ni,nj,nk,-1,' ',ip1,-1,-1,' ',varName)
    if (key < 0) then
      write(*,*) 'bias_getTrialPredictors: Unable to find trial field - varname,ip1= ',varName,ip1
      call utl_abort('bias_getTrialPredictors')
    end if
    varTrial1(:,:) = varTrial1(:,:) - varTrial2(:,:)
    ierr = utl_ezsint(trialGZ5m50,varTrial1,nobs,1,1,ni,nj,1)
    if( ierr < 0 ) call utl_abort('bias_getTrialPredictors: error interpolating GZ5-GZ50')

    ! GZ1 - GZ10
    varName = 'GZ'
    ip1 = 1820  ! 1hPa
    key = utl_fstlir(varTrial1,nulfst,ni,nj,nk,-1,' ',ip1,-1,-1,' ',varName)
    if ( key < 0 ) then
      write(*,*) 'bias_getTrialPredictors: Unable to find trial field - varname,ip1= ',varName,ip1
      call utl_abort('bias_getTrialPredictors')
    end if
    ip1 = 10
    key = utl_fstlir(varTrial2,nulfst,ni,nj,nk,-1,' ',ip1,-1,-1,' ',varName)
    if ( key < 0 ) then
      write(*,*) 'bias_getTrialPredictors: Unable to find trial field - varname,ip1= ',varName,ip1
      call utl_abort('bias_getTrialPredictors')
    end if
    varTrial1(:,:) = varTrial1(:,:) - varTrial2(:,:)
    ierr = utl_ezsint(trialGZ1m10,varTrial1,nobs,1,1,ni,nj,1)
    if( ierr < 0 ) call utl_abort('bias_getTrialPredictors: error interpolating GZ1-GZ10')

    ! Determine grid size and EZSCINT ID for TG

    timeIndex = nint((real(tim_nStepObs,8)+1.0d0)/2.0d0)
    write(flnum,'(I2.2)') timeIndex
    trialfile = './trlm_'//trim(flnum)

    inquire(file=trim(trialfile),exist=trialExists)

    if( .not.trialExists ) then
      write(*,*) 'File missing for bias_getTrialPredictors_forTG=',trialfile
      call utl_abort('bias_getTrialPredictors_forTG:DID NOT FIND THE MODEL-LEVEL TRIAL FIELD FILE')
    else
      nultrlm = 0 
      ierr = fnom(nultrlm,trim(trialfile),'RND+OLD+R/O',0)
        if ( ierr /= 0 ) then
          write(*,*) 'bias_getTrialPredictors_forTG: trialfilename=',trialfile
          call utl_abort('bias_getTrialPredictors_forTG: Error opening trial file')
        end if
      
      ierr = fstouv(nultrlm,'RND+OLD')

      write(*,*) 'bias_getTrialPredictors_forTG :', trialfile
      write(*,*) 'opened as unit file ', nultrlm

    end if

    dateo  = -1
    cletiket = ' '
    ip1    = -1
    ip2    = -1
    ip3    = -1
    cltypvar = ' '
    varName = 'TG'
    key = fstinf( nultrlm,                                    & ! IN
                  ni,nj,nk,                                  & ! OUT
                  dateo,cletiket,ip1,ip2,ip3,cltypvar,varName )! IN
    if (key < 0) then
      write(6,*)
      write(6,*) 'Bias_getTrialPredictors: Unable to find trial field = ',varName
      stop
    end if
    ierr = fstprm( key,                                               & ! IN
                   dateo, deet, npas, ni, nj, nk, nbits,              & ! OUT
                   datyp, ip1, ip2, ip3, cltypvar, varName, cletiket, & ! OUT
                   clgrtyp, ig1, ig2, ig3,                            & ! OUT
                   ig4, swa, lng, dltf, ubc, extra1, extra2, extra3 )   ! OUT
    trlgid = ezqkdef(ni,nj,clgrtyp,ig1,ig2,ig3,ig4,nultrlm)
    write(*,*) 'bias_getTrialPredictors_forTG: trlgid=',trlgid,ni,nj

    allocate(varTrial3(ni,nj))

    iset = ezdefset(obsgid,trlgid)

    ! TG (skin temperature)
    varName = 'TG'
    ip1 = -1
    key = utl_fstlir(varTrial3,nultrlm,ni,nj,nk,-1,' ',ip1,-1,-1,' ',varName)
    if (key < 0) then
      write(*,*) 'bias_getTrialPredictors_forTG: Unable to find trial field - varname,ip1= ',varName,ip1
      call utl_abort('bias_getTrialPredictors_forTG')
    end if
    ierr = utl_ezsint(trialTG,varTrial3,nobs,1,1,ni,nj,1)
    if( ierr < 0 ) call utl_abort('bias_getTrialPredictors_forTG: error interpolating TG')

    if( trialTG(1) > 150.0d0) then
      write(*,*) 'bias_getTrialPredictors_forTG: converting TG from Kelvin to deg_C'
      trialTG(:) = trialTG(:) - MPC_K_C_DEGREE_OFFSET_R8
    end if
 
    ierr = fstfrm(nulfst)
    ierr = fclos(nulfst)

    ierr = fstfrm(nultrlm)
    ierr = fclos(nultrlm)

    deallocate(dlonfld)
    deallocate(dlatfld)
    deallocate(varTrial1)
    deallocate(varTrial2)
    deallocate(varTrial3)

    write(*,*) 'bias_getTrialPredictors done'

  end subroutine bias_getTrialPredictors

  !---------------------------------------
  ! bias_cvToCoeff
  !--------------------------------------
  subroutine bias_cvToCoeff(cv_bias)
    implicit none

    real(8)  :: cv_bias(:)
    integer  :: index_cv,iSensor,iChannel,iPredictor,iScan
    integer  :: nsize,ierr

    if( mpi_myid == 0 ) write(*,*) 'bias_cvToCoeff starting'
 
    if( mpi_myid == 0 ) then

      index_cv = 0
      ! initialize of coeffIncr
      do iSensor = 1, tvs_nSensors
        bias(iSensor)%coeffIncr(:,:) = 0.0d0
        bias(iSensor)%coeffIncr_fov(:,:) = 0.0d0
      end do

      do iSensor = 1, tvs_nSensors
        do iChannel = 1, bias(iSensor)%numChannels
          do iPredictor = 1, bias(iSensor)%numActivePredictors
            if( iPredictor == 1) then
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

    call flush(6)

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

    if( mpi_myid == 0 ) write(*,*) 'bias_cvToCoeff finishes', maxval(bias(tvs_nSensors)%coeffIncr)
  end subroutine bias_cvToCoeff

  !-----------------------------------
  !getPredictors
  !---------------------------------- 
  subroutine bias_getPredictors(predictor,index_header,index_obs,lobsSpaceData)
    implicit none

    real(8)  :: predictor(NumPredictors)
    integer  :: index_header,index_obs
    type(struct_obs)  :: lobsSpaceData

    integer  :: iSensor,iPredictor
    integer  :: activePredictors(tvs_nSensors)

    predictor(:) = 0.0d0
    iSensor = tvs_lsensor(tvs_ltovsno(index_header))

    do iPredictor = 1, NumPredictors

        if( iPredictor == 1 ) then
          ! constant
          predictor(iPredictor) = 1.0d0
        else if( iPredictor == 2 ) then
          ! GZ300-GZ1000 (dam) /1000
          predictor(iPredictor) = (trialGZ300m1000(index_obs)/1000.0d0)  !-0.85d0 tried de-biasing
        else if( iPredictor == 3 ) then
          ! GZ50-GZ200 (dam) /1000
          predictor(iPredictor) = (trialGZ50m200(index_obs)/1000.0d0)  !-0.85d0 tried de-biasing
        else if( iPredictor == 4 ) then
          ! skin temperature (C) /10
          predictor(iPredictor) = trialTG(index_obs)
        else if( iPredictor == 6 ) then
          ! GZ1-GZ10 (dam) /1000
          predictor(iPredictor) = trialGZ1m10(index_obs)/1000.0d0
        else if( iPredictor == 7 ) then
          ! GZ5-GZ50 (dam) /1000
          predictor(iPredictor) = (trialGZ5m50(index_obs)/1000.0d0)  !-0.72d0 tried de-biasing
        else if( iPredictor == 9 ) then
          ! satellite zenith angle (degrees/100)^1
          predictor(iPredictor) = ((obs_headElem_i(lobsSpaceData,OBS_SZA,index_header)-9000)/100.0)
        else if( iPredictor == 10 ) then
          ! satellite zenith angle (degrees/100)^2
          predictor(iPredictor) = ((obs_headElem_i(lobsSpaceData,OBS_SZA,index_header)-9000)/100.0)**2
        else if( iPredictor == 11 ) then
          ! satellite zenith angle (degrees/100)^3
          predictor(iPredictor) = ((obs_headElem_i(lobsSpaceData,OBS_SZA,index_header)-9000)/100.0)**3
        end if

    end do

  end subroutine bias_getPredictors

  !--------------------------------------------------
  ! calcMeanPredictors
  !---------------------------------------------------
  subroutine bias_calcMeanPredictors(lobsSpaceData)
    implicit none

    type(struct_obs)  :: lobsSpaceData

    integer  :: index_header,iobs, idatyp
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

    ! loop over all observations
    iobs = 0
    call obs_set_current_header_list(lobsSpaceData,'TO')
    HEADER: do
      index_header = obs_getHeaderIndex(lobsSpaceData)
      if (index_header < 0) exit HEADER
      
      IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,index_header)

      IF ( .not.  tvs_isIdBurpTovs(IDATYP) ) cycle HEADER       

      iobs = iobs + 1
      iSensor = tvs_lsensor(tvs_ltovsno(index_header))
      call bias_getPredictors(predictor,index_header,iobs,lobsSpaceData)
      do iPredictor = 2, bias(iSensor)%NumActivePredictors
          minPredictors(iPredictor) = min(minPredictors(iPredictor),predictor(bias(iSensor)%predictorIndex(iPredictor)))
          maxPredictors(iPredictor) = max(maxPredictors(iPredictor),predictor(bias(iSensor)%predictorIndex(iPredictor)))
          meanPredictors(iPredictor) = meanPredictors(iPredictor) + predictor(bias(iSensor)%predictorIndex(iPredictor))
          meanAbsPredictors(iPredictor) = meanAbsPredictors(iPredictor) + abs(predictor(bias(iSensor)%predictorIndex(iPredictor)))
          countObs(iPredictor) = countObs(iPredictor) + 1
      end do
    end do HEADER

    !  do iPredictor=1,bias(iSensor)%NumActivePredictors
    do iPredictor=2,bias(iSensor)%NumActivePredictors
      if(countObs(iPredictor).gt.0) then
        meanPredictors(iPredictor)=meanPredictors(iPredictor)/real(countObs(iPredictor),8)
        meanAbsPredictors(iPredictor)=meanAbsPredictors(iPredictor)/real(countObs(iPredictor),8)
        write(*,*) 'bias_calcMeanPredictors: mean(abs),mean,min,max=',iPredictor, &
                   meanAbsPredictors(iPredictor),meanPredictors(iPredictor),minPredictors(iPredictor),maxPredictors(iPredictor)
      endif
    enddo

  end subroutine bias_calcMeanPredictors

  !---------------------------------------------
  ! calcBias_ad
  !---------------------------------------------
  subroutine bias_calcBias_ad(cv_out,cv_dim,obsColumnIndex,multFactor,lobsSpaceData)
    implicit none

    integer  :: cv_dim
    real(8)  :: cv_out(cv_dim)
    integer  :: obsColumnIndex
    real(8)  :: multFactor
    type(struct_obs)  :: lobsSpaceData

    integer  :: index_header,index_body,iobs, idatyp
    integer  :: iSensor,iChannel,iPredictor,index_cv,nsize,ierr
    integer  :: iScan, iFOV, jj
    real(8)  :: predictor(NumPredictors)
    real(8),pointer  :: cv_bias(:)
    real(8)  :: biasCor

    if( .not.lvarbc ) return

    if( mpi_myid == 0 ) write(*,*) 'Starting bias_calcBias_ad'

    if( mpi_myid == 0 ) then
      if( cvm_subVectorExists(cvm_Bias) ) then
        cv_bias => cvm_getSubVector(cv_out,cvm_Bias)
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
    call obs_set_current_header_list(lobsSpaceData,'TO')
    HEADER: do
      index_header = obs_getHeaderIndex(lobsSpaceData)
      if (index_header < 0) exit HEADER
      idatyp = obs_headElem_i(lobsSpaceData,OBS_ITY,index_header)
      IF ( .not.  tvs_isIdBurpTovs(idatyp) ) cycle HEADER  

      iobs = iobs + 1
      iSensor = tvs_lsensor(tvs_ltovsno(index_header))
      call bias_getPredictors(predictor,index_header,iobs,lobsSpaceData)
      call obs_set_current_body_list(lobsSpaceData, index_header)
      iFov = obs_headElem_i(lobsSpaceData,OBS_FOV,index_header)

      BODY: do
        index_body = obs_getBodyIndex(lobsSpaceData)
        if (index_body < 0) exit BODY

        if (obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) .ne. 1) cycle BODY  !PDD
        iChannel = nint(obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body))
        iChannel = max(0,min(iChannel,tvs_maxChannelNumber+1))
        iChannel = iChannel-tvs_channelOffset(iSensor)

        biasCor = multFactor*obs_bodyElem_r(lobsSpaceData,obsColumnIndex,index_body)
       
        do iPredictor = 1, bias(iSensor)%numActivePredictors
          jj = bias(iSensor)%PredictorIndex(iPredictor)

          if (iPredictor == 1) then
            if(bias(iSensor)%numScan .gt. 1) then
              iScan = iFov
            else
              iScan = 1
            endif
            bias(iSensor)%coeffIncr_fov(iChannel,iScan) = bias(iSensor)%coeffIncr_fov(iChannel,iScan) &
                 + predictor(jj)*biasCor
          else
            bias(iSensor)%coeffIncr(iChannel,iPredictor) = bias(iSensor)%coeffIncr(iChannel,iPredictor) & 
               + predictor(jj)*biasCor
            
          end if
        end do
      end do BODY

    end do HEADER

    ! put the coefficients into the control vector
    call bias_cvToCoeff_ad(cv_bias)

    if (mpi_myid == 0) then
      write(*,*) 'bias_calcBias_ad: maxval(cv_bias)=',maxval(cv_bias(:))
      write(*,*) 'bias_calcBias_ad: (cv_bias(69)=',cv_bias(69)
    end if

  end subroutine bias_calcBias_ad

  !----------------------------------------------------
  ! cvToCoeff_ad
  !----------------------------------------------------
  subroutine bias_cvToCoeff_ad(cv_bias)
    implicit none

    real(8)  :: cv_bias(:)
    integer  :: index_cv,iSensor,iChannel,iPredictor,iScan
    integer  :: nChan,nPre,nScan
    integer  :: nsize,ierr
    real(8)  :: dl_jbias
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

    if (mpi_myid == 0) then
      index_cv = 0
      dl_jbias = 0.0
      do iSensor = 1, tvs_nSensors
        do iChannel = 1, bias(iSensor)%numChannels
          do iPredictor = 1, bias(iSensor)%numActivePredictors
            if (iPredictor == 1) then
              do iScan = 1, bias(iSensor)%numScan
                index_cv  =  index_cv  +  1
                cv_bias(index_cv) = bias(iSensor)%stddev(iChannel,iPredictor) * bias(iSensor)%coeffIncr_fov(iChannel,iScan)
              end do
            else
              index_cv = index_cv + 1
              cv_bias(index_cv) = bias(iSensor)%stddev(iChannel,iPredictor) * bias(iSensor)%coeffIncr(iChannel,iPredictor)
              dl_jbias = dl_jbias+(cv_bias(index_cv)*cv_bias(index_cv))/2.0  
            end if
          end do
        end do
      end do
      write(*,'(6X,"SIMVAR:  Jbias = ",G23.16,6X)') dl_jbias
    end if
    
  end subroutine bias_cvToCoeff_ad

  !-----------------------------------------
  ! writeBias
  !-----------------------------------------
  subroutine bias_writeBias(cv_in,cv_dim)
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
    integer             :: numCoefFile,jj,kk
    character(len=10)   :: coefInstrName(tvs_nSensors), temp_instrName, instrName
    character(len=25)   :: filecoeff
    logical             :: coeffExists

    if( .not.lvarbc ) return

    if (mpi_myid == 0) then
      if ( cvm_subVectorExists(cvm_Bias) ) then
        cv_bias => cvm_getSubVector(cv_in,cvm_Bias)
        write(*,*) 'bias_writeBias: maxval(cv_bias)=',maxval(cv_bias(:))
      else
        write(*,*) 'bias_writeBias: control vector does not include bias coefficients'
        return
      endif
    end if

    call bias_cvToCoeff(cv_bias)

    ! write out bias coefficient increments in ascii file
    nulfile_inc = 0
    ierr = fnom(nulfile_inc,'./satbias_increment.dat','FTN+FMT',0)

    do iSensor = 1, tvs_nSensors
      write(nulfile_inc,'(/,1X,"Sensor Index=",I3,", Satellite Name=",A15,", Instrument Name=",A15)') &
             iSensor,tvs_satelliteName(iSensor),tvs_instrumentName(iSensor)
      do iChannel = 1, bias(iSensor)%numChannels
        if ( sum(bias(iSensor)%coeffIncr(iChannel,:)) /= 0.0d0) &
          write(nulfile_inc,'(3X,"Channel number=",I4)') iChannel
        do iPredictor = 2, bias(iSensor)%numActivePredictors
          if ( bias(iSensor)%coeffIncr(iChannel,iPredictor) /= 0.0d0) &
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
        if (sum(bias(iSensor)%coeffIncr_fov(iChannel,:)) /= 0.0d0) &
          write(nulfile_fov,'(3X,"Channel number=",I4)') iChannel 
        if (sum(bias(iSensor)%coeffIncr_fov(iChannel,:)) /= 0.0d0) & 
          write(nulfile_fov,*) bias(iSensor)%coeffIncr_fov(iChannel,:)
      end do
    end do
    ierr=fclos(nulfile_fov)

    ! Find the background coeff_file number and name
    do iSensor = 1, tvs_nSensors

        numCoefFile = 0
        jj = 0
        do jSensor = 1, tvs_nSensors
          temp_instrName = InstrNametoCoeffFileName(tvs_instrumentName(jSensor))
          filecoeff = 'coeff_file_'//trim(temp_instrName)//''
          inquire(file=trim(filecoeff),exist = coeffExists)

          if ( coeffExists ) then
            numCoefFile = numCoefFile + 1
            jj = jj + 1
            coefInstrName(jj) = temp_instrName
          end if
          if (jSensor > 1 ) then
            do kk = 1, jj-1
              if ( temp_instrName == coefInstrName(kk) ) then
                numCoefFile = numCoefFile - 1
                jj = jj -1
              end if
            end do
          end if 
        end do

    end do
    write(*,*) 'DPP-tst', numCoefFile
    write(*,*) 'DPP-tst2', coefInstrName(1:numCoefFile)

    ! update coeff_file_instrument and write out
    if( mpi_myid == 0) then
      do iInstr=1, numCoefFile 
        biasCoeff_bg(:,:,:) = 0.0
        fovbias_bg(:,:,:) = 0.0
        BgFileName ='./coeff_file_'//coefInstrName(iInstr)//''
        call bias_updateCoeff(tvs_nSensors,NumPredictors,BgFileName)
      end do
    end if

  end subroutine bias_writeBias

  !--------------------------------------
  ! updateCoeff
  !--------------------------------------
  subroutine bias_updateCoeff(maxsat,maxpred,coeff_file)
    implicit none

    ! There are three parts in this subroutine, read, update and write out the coeff files
    ! IN
    integer            :: maxsat, maxpred
    character(len=80)  :: coeff_file

    ! local 
    character(len=10)  :: sats(maxsat)        ! dim(maxsat), satellite names
    integer            :: chans(maxsat, maxNumChannels)       ! dim(maxsat, maxchan), channel numbers
    real(8)            :: fovbias(maxsat,maxNumChannels,maxfov)     ! dim(maxsat,maxchan,maxfov), bias as F(fov)
    real(8)            :: coeff(maxsat,maxNumChannels,maxpred)       ! dim(maxsat,maxchan,maxpred+1)
    integer            :: nsat, nfov
    integer            :: nchan(maxsat)       ! dim(maxsat), number of channels
    integer            :: npred(maxsat, maxNumChannels)       ! dim(maxsat, maxchan), number of predictors
    character(len=6)   :: cinstrum    ! string: instrument (e.g. AMSUB)
    character(len=2)   :: ptypes(maxsat,maxNumChannels,maxpred) ! dim(maxsat,maxNumChannelsmaxchan,maxpred)

    ! LOCAL for reading background coeff file
    character(len=10)  :: sat
    character(len=120) :: line
    integer            :: chan, ndata, nbfov, nbpred
    integer            :: i, j, k, ier, istat, ii, iSat
    logical            :: newsat, verbose
    real               :: dummy
    integer            :: iun
    integer            :: fnom, fclos
    external           :: fnom, fclos

    ! update coeff files
    real               :: fovbias_an(maxsat,maxNumChannels,maxfov)
    real               :: coeff_an(maxsat,maxNumChannels,maxpred) 
    integer            :: iSensor, jChannel , iFov, iPred, totPred, jChan
    character(len=10)  :: tmp_SatName, tmp_InstName 

    ! write out files 
    integer            :: iuncoef2, ierr, numPred
    character(len=80)  :: filename2

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

    write(*,*) 'Reading coeff file starts'  
    iun = 0
    ier = FNOM(iun,coeff_file,'FMT',0)
    if( ier == 0 ) then 
      write(*,*) 'Bias correction coefficient file open = ', coeff_file
      read(iun,*,IOSTAT=istat)
      if( istat < 0 ) then
        write(*,*) 'ERROR- File appears empty.'
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
          read(line,'(T53,A8,1X,A6,1X,I6,1X,I8,1X,I2,1X,I3)',IOSTAT=istat) sat, cinstrum, chan, ndata, nbpred, nbfov
          do i = 1, maxsat
            if ( trim(sats(i)) == trim(sat) ) then
              newsat = .false.
              ii = i
            end if
          end do
          if ( newsat ) then
            ii = ii + 1
            sats(ii) = sat
            if (ii > 1) nchan(ii-1) = j
            j = 1
          else
            j = j + 1
          end if
          chans(ii, j) = chan
          npred(ii, j) = nbpred
          read(iun,'(A)',IOSTAT=istat) line
          if( nbpred > 0 ) then
            read(line,'(T8,6(1X,A2))',IOSTAT=istat) (ptypes(ii,j,k),k=1,nbpred)
          end if
          read(iun,*,IOSTAT=istat) (fovbias(ii,j,k),k=1,nbfov)
          if( nbpred > 0 ) then
            read(iun,*,IOSTAT=istat) (coeff(ii,j,k),k=1,nbpred+1)
          else
            read(iun,*,IOSTAT=istat) dummy
          end if
        end if
      end do

      if( ii == 0 ) then
        write(*,*) ' ERROR - No data read from coeff file!'
        call utl_abort('bias_updateCoeff')
      end if
      nsat      = ii
      nfov      = nbfov
      nchan(ii) = j

      if( verbose ) then
        write(*,*) ' '
        write(*,*) ' ------------- BIAS CORRECTION COEFFICIENT FILE ------------------ '
        write(*,*) ' '
        write(*,*) ' Number of satellites =     ', nsat
        write(*,*) ' Number of FOV =            ', nfov
        write(*,*) ' Max number of predictors = ', maxval(npred)
        write(*,*) ' '
        do i = 1, nsat
          write(*,*) '  Satellite = ' // sats(i)
          write(*,*) '     Number of channels = ', nchan(i)
          write(*,*) '     predictors, fovbias, coeff for each channel: '
          do j = 1, nchan(i)
            write(*,*) i, chans(i,j)
            if( npred(i,j) > 0 ) then
              write(*,'(6(1X,A2))') (ptypes(i,j,k),k=1,npred(i,j))
            else
              write(*,'(A)') 'No predictors'
            end if
            write(*,*) (fovbias(i,j,k),k=1,nfov)
            write(*,*) (coeff(i,j,k),k=1,npred(i,j)+1)
          end do
        end do
          write(*,*) ' '
      end if
    else
      write(*,*) 'READ_COEFF: ERROR - Problem opening the coeff file!'
      call utl_abort('bias_updateCoeff')
    end if 
    ier = FCLOS(iun)

    write(*,*) 'Reading coeff file done', coeff_file

    !
    !- 2.update coeff and fovbias  
    !
    coeff_an(:,:,:) = coeff(:,:,:)
    fovbias_an(:,:,:) = fovbias(:,:,:)

    do i = 1, nsat
      do iSensor = 1, tvs_nSensors
        ! for Satellite Name
        tmp_SatName = SatNameinCoeffFile(tvs_satelliteName(iSensor))
        ! for Instrument Name
        tmp_InstName = InstrNameinCoeffFile(tvs_instrumentName(iSensor))
        if( trim(tmp_SatName) == trim(sats(i)) .and. trim(tmp_InstName) == trim(cinstrum) ) then
          do j = 1, nchan(i)
            do jChan = 1, bias(iSensor)%numChannels  
              if( chans(i,j) == jChan ) then
                ! part 1 for coeffIncr
                do iFov = 1, nfov
                  if( bias(iSensor)%coeffIncr_fov(jChan,iFov) /= 0.0d0 ) then
                    fovbias_an(i,j,iFov) = fovbias(i,j,iFov) + bias(iSensor)%coeffIncr_fov(jChan,iFov)
                  end if
                end do ! iFov
                ! part 2 for coeffIncr_fov
                totPred  = bias(iSensor)%NumActivePredictors 
                do iPred = 1, totPred
                  if( iPred == 1 ) then
                    coeff_an(i,j,iPred) = coeff(i,j,iPred)
                  else
                    if( bias(iSensor)%coeffIncr(jChan,iPred) /= 0.0d0 ) then
                      coeff_an(i,j,iPred) = coeff(i,j,iPred) + bias(iSensor)%coeffIncr(jChan,iPred)
                    end if
                  end if
                end do ! iPred
              end if
            end do ! jChan
          end do !j 
        end if ! sat and instr
      end do !iSendor
    end do ! i

    !
    !- 3. Write out updated_coeff
    ! 

    iuncoef2 = 0
    filename2 ='./Newcoeffs_'//cinstrum//'' 
    ierr = fnom(iuncoef2, filename2,'FTN+FMT',0)

    write(*,*) 'Write in bias_updateCoeff'
   
    do i = 1, nsat
      do j = 1, nchan(i)      
        numPred = npred(i,j)
        write(iuncoef2, '(A52,A8,1X,A6,1X,I6,1X,I8,1X,I2,1X,I3)') &
          'SATELLITE, INSTRUMENT, CHANNEL, NOBS, NPRED, NSCAN: ', sats(i), cinstrum, chans(i,j), ndata, numPred, nfov
        write(iuncoef2, '(A7,6(1X,A2))') 'PTYPES:', (ptypes(i,j,k), k=1,numPred)
        write(iuncoef2,'(120(1x,ES17.10))') (fovbias_an(i,j,k),k=1,nfov)
        write(iuncoef2,*) (coeff_an(i,j,k),k=1,numPred+1)
      end do
    end do

    ierr = fclos(iuncoef2) 
   
    write(*,*) 'Finish writing in bias_updateCoeff'
    
  end subroutine bias_updateCoeff

  !------------------
  ! read_coeff
  !------------------
  subroutine read_coeff(maxsat,maxpred,coeff_file,sats,chans,nsat,nchan,nfov,cinstrum)
    implicit none

    ! IN
    integer            :: maxsat, maxpred
    character(len=80)  :: coeff_file
    ! OUT
    character(len=10)  :: sats(maxsat)           ! dim(maxsat), satellite names
    integer            :: chans(maxsat, maxNumChannels)  ! dim(maxsat, maxchan), channel numbers
    integer            :: nsat, nfov
    integer            :: nchan(maxsat)       ! dim(maxsat), number of channels
    character(len=6)   :: cinstrum    ! string: instrument (e.g. AMSUB)

    !local 
    real(8)            :: fovbias(maxsat,maxNumChannels,maxfov)     ! dim(maxsat,maxchan,maxfov), bias as F(fov)
    real(8)            :: coeff(maxsat,maxNumChannels,maxfov)       ! dim(maxsat,maxchan,maxpred+1)
    integer            :: npred(maxsat, maxNumChannels)       ! dim(maxsat, maxchan), number of predictors
    character(len=2)   :: ptypes(maxsat,maxNumChannels,maxpred) ! dim(maxsat,maxchan,maxpred)
    character(len=10)  :: sat
    character(len=120) :: line
    integer            :: chan
    integer            :: ndata, nbfov, nbpred, i, j, k, ier, istat, ii
    logical            :: newsat, verbose
    real               :: dummy
    integer            :: iun

    integer            :: fnom, fclos
    external           :: fnom, fclos

    !   sats(nsat)            = satellite names
    !   chans(nsat,nchan(i))  = channel numbers of each channel of each satellite i
    !   npred(nsat,nchan(i))  = number of predictors for each channel of each satellite i
    !   nsat, nchan, nfov, cinstrum (output) are determined from file
    !   if returned nsat = 0, coeff_file was empty

    !   fovbias(i,j,k)        = bias for satellite i, channel j, FOV k   k=1,nfov
    !     if FOV not considered for instrument, nfov = 1 and fovbias is global bias for channel
    !   coeff(i,j,1)          = regression constant
    !   coeff(i,j,2), ..., coeff(i,j,npred(i,j)) = predictor coefficients
    !   maxpred (input) is max number of predictors
    !   maxsat (input)  is max number of satellites

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

    write(*,*) 'Reading coeff file starts'
    iun = 0
    ier = FNOM(iun,coeff_file,'FMT',0)
    if (ier == 0)  then 
      write(*,*) 'Bias correction coefficient file open = ', coeff_file
      read(iun,*,iostat=istat)
      if( istat < 0 ) then
        write(*,*) 'Error: File appears empty.'
        return
      end if 
      rewind(iun)
      ii = 0
      ! Loop over the satellites/channels in the file
      do
        read(iun,'(A)',iostat=istat) line
        if ( istat < 0 ) exit
          if ( line(1:3) == 'SAT' ) then
            newsat = .true.
            read(line,'(T53,A8,1X,A6,1X,I6,1X,I8,1X,I2,1X,I3)',iostat=istat) sat, cinstrum, chan, ndata, nbpred, nbfov
            do i = 1, maxsat
              if ( trim(sats(i)) == trim(sat) ) then
                newsat = .false.
                ii = i
              endif
            end do
            if ( newsat ) then
              ii = ii + 1
              sats(ii) = sat

              if (ii > 1) nchan(ii-1) = j
              j = 1
            else
              j = j + 1
            end if
            chans(ii, j) = chan
            npred(ii, j) = nbpred

            read(iun,'(A)',iostat=istat) line
            if ( nbpred > 0 ) then
              read(line,'(T8,6(1X,A2))',iostat=istat) (ptypes(ii,j,k),k=1,nbpred)
            end if
            read(iun,*,iostat=istat) (fovbias(ii,j,k),k=1,nbfov)
            if ( nbpred > 0 ) then
              read(iun,*,iostat=istat) (coeff(ii,j,k),k=1,nbpred+1)
            else
              read(iun,*,iostat=istat) dummy
            end if
        endif
      end do
      if ( ii == 0 ) then
        write(*,*) ' Error - No data read from coeff file!'
        call utl_abort('read_coeff')
      endif

      nsat      = ii
      nfov      = nbfov
      nchan(ii) = j

      if (verbose) then
        write(*,*) ' '
        write(*,*) ' ------------- BIAS CORRECTION COEFFICIENT FILE ------------------ '
        write(*,*) ' '
        write(*,*) ' Number of satellites =     ', nsat
        write(*,*) ' Number of FOV =            ', nfov
        write(*,*) ' Max number of predictors = ', maxval(npred)
        write(*,*) ' '
        do i = 1, nsat
          write(*,*) '  Satellite = ' // sats(i)
          write(*,*) '     Number of channels = ', nchan(i)
          write(*,*) '     predictors, fovbias, coeff for each channel: '
          do j = 1, nchan(i)
            write(*,*) i, chans(i,j)
            if ( npred(i,j) > 0 ) then
              write(*,'(6(1X,A2))') (ptypes(i,j,k),k=1,npred(i,j))
            else
              write(*,'(A)') 'No predictors'
            end if
            write(*,*) (fovbias(i,j,k),k=1,nfov)
            write(*,*) (coeff(i,j,k),k=1,npred(i,j)+1)
          end do
        end do
          write(*,*) ' '
      end if
    else
      write(*,*) 'Read_coeff:Error - Problem opening the coeff file!'
      call utl_abort('read_coeff')
    end if

    ier = FCLOS(iun)

    write(*,*) 'Reading coeff file done', coeff_file

  end subroutine read_coeff

  !----------------------
  ! Finalize
  !----------------------
  subroutine bias_Finalize
    implicit none

    integer    :: iSensor

    if (.not.lvarbc) return

    deallocate(trialGZ300m1000)
    deallocate(trialGZ50m200)
    deallocate(trialGZ1m10)
    deallocate(trialGZ5m50)
    deallocate(trialTG)

    do iSensor = 1, tvs_nSensors
      deallocate(bias(iSensor)%stddev)
      deallocate(bias(iSensor)%coeffIncr)
      deallocate(bias(iSensor)%predictorIndex)
      deallocate(bias(iSensor)%coeffIncr_fov)
    end do
  end subroutine bias_Finalize 

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
    if( trim(temp_instrName) == 'mhs' ) then
      nameOut = 'amsub'
    else if( trim(temp_instrName) == 'goesimager' ) then
      nameOut = 'cgoes'
    else if( trim(temp_instrName) == 'gmsmtsat' ) then
      nameOut = 'mtsat'
    else if( trim(temp_instrName) == 'mviri' ) then
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

    if( trim(nameIn) == 'MHS' ) then
      nameOut = 'AMSUB'
    else if( trim(nameIn) == 'GOESIMAGER' ) then
      nameOut = 'CGOES' 
    else if( trim(nameIn) == 'GMSMTSAT' ) then
      nameOut = 'MTSAT' 
    else if( trim(nameIn) == 'MVIRI' ) then
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

    if( trim(nameIn) == 'MSG2' ) then
      nameOut = 'METSAT9'
    else if( trim(nameIn) == 'MSG3' ) then
      nameOut = 'METSAT10' 
    else if( trim(nameIn) == 'METEOSAT7' ) then
      nameOut = 'METSAT7' 
    else 
      nameOut = nameIn
    end if

  end function SatNameinCoeffFile

end module biascorrection_mod

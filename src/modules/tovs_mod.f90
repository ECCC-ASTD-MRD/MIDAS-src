
module tovs_mod
  ! MODULE tovs_mod (prefix='tvs' category='5. Observation operators')
  !
  !:Purpose:  Derived types, public variables and procedures related to RTTOV
  !
  
  use rttovInterfaces_mod
  use rttov_types, only :              &
       rttov_coefs                    ,&
       rttov_fast_coef                ,&
       rttov_scatt_coef               ,&
       rttov_options                  ,&
       rttov_options_scatt            ,&
       rttov_profile                  ,&
       rttov_profile_cloud            ,&
       rttov_radiance                 ,&
       rttov_transmission             ,&
       rttov_chanprof                 ,&
       rttov_emissivity               ,&
       rttov_scatt_emis_retrieval_type,&
       rttov_geometry                 ,&
       rttov_profile_scatt_aux        ,&
       rttov_radiance2                ,&
       rttov_reflectivity             ,&
       eddington_sfc_type
  use rttov_const, only :           &
       platform_name               ,&
       nplatforms                  ,&
       inst_name                   ,&
       ninst                       ,&
       inst_id_goesim              ,&
       inst_id_gmsim               ,&
       inst_id_mtsatim             ,&
       inst_id_amsua               ,&
       inst_id_mhs                 ,&
       sensor_id_mw                ,&
       sensor_id_po                ,&
       platform_id_jpss            ,&
       platform_id_himawari        ,&
       platform_id_eos             ,&
       errorstatus_success         ,&
       surftype_land               ,&
       surftype_seaice             ,&
       surftype_sea                ,&
       watertype_ocean_water       ,&
       ngases_max                  ,&
       gas_id_mixed                ,&
       gas_unit_specconc           ,&
       interp_rochon_loglinear_wfn ,&
       zenmax                      ,&
       zenmaxv9                    ,&
       o3min                       ,&
       o3max                       ,&
       tmin                        ,&
       tmax                        ,&
       qmin                        ,&
       qmax                        ,&
       elevmax                     ,&
       wmax                        ,&
       pmin                        ,&
       pmax                        ,&
       errorstatus_fatal           ,&
       sensor_id_po                ,&
       min_reflectivity            ,&
       min_radiance_radar
  use parkind1, only : jpim, jplm, jprb
  use midasMpi_mod
  use message_mod
  use codtyp_mod
  use utilities_mod
  use obsSpaceData_mod
  use mathPhysConstants_mod
  use climatologies_mod
  use bufr_mod
  use columnData_mod 
  use mod_rttov_emis_atlas
  use verticalCoord_mod
  use codePrecision_mod
  use humidityLimits_mod
  use interpolation_mod
  use surfaceEmissivity_mod
  use clibInterfaces_mod
  use ramDisk_mod

  implicit none
  save
  private

  ! Public procedures
  public :: tvs_allocateSurfaceParameters, tvs_allocateEmissivity
  public :: tvs_fillProfiles, tvs_rttov, tvs_printDetailledOmfStatistics, tvs_allocTransmission
  public :: tvs_deallocateProfilesNlTlAd
  public :: tvs_setupAlloc,tvs_setup, tvs_isIdBurpTovs, tvs_isIdBurpHyperSpectral, tvs_isIdBurpInst, tvs_getAllIdBurpTovs
  public :: tvs_isNameHyperSpectral
  public :: tvs_isNameGeostationary
  public :: tvs_getInstrumentId, tvs_getPlatformId, tvs_mapSat, tvs_mapInstrum
  public :: tvs_ChangedStypValue
  public :: tvs_getLocalChannelIndexFromChannelNumber
  public :: tvs_getProfile, tvs_interp_sfc
  public :: tvs_getCorrectedSatelliteAzimuth
  public :: tvs_isInstrumUsingCLW, tvs_isInstrumUsingHydrometeors, tvs_getChannelNumIndexFromPPP
  public :: tvs_isInstrumAllskyTtAssim, tvs_isInstrumAllskyHuAssim
  public :: tvs_useSfcEmissObsSpace, tvs_emis_read_climatology, tvs_pcnt_box
  public :: tvs_rttov_tl, tvs_rttov_ad, tvs_rttov_k
  public :: tvs_checkAllskyChanNum, tvs_isChanNumInAllskyNamtovList
  
  type surface_params
    real(8)   :: albedo   ! surface albedo (0-1)
    real(8)   :: ice      ! ice cover (0-1)
    real(8)   :: snow     ! snow cover (0-1)
    real(8)   :: pcnt_wat ! water percentage in pixel containing profile (0-1)
    real(8)   :: pcnt_reg ! water percentage in an area around profile (0-1)
    integer   :: ltype    ! surface type (1,...,20)
  end type surface_params

  ! Public module parameters
  integer, public, parameter :: tvs_maxChannelNumber    = 8461 ! Max. value for channel number
  integer, public, parameter :: tvs_maxNumberOfChannels = 2211 ! Max. no. of channels (for one profile/spectra)
  integer, public, parameter :: tvs_maxNumberOfSensors  = 100  ! Max no sensors to be used
  real(pre_obsReal), public, parameter :: tvs_defaultEmissivity = 0.95

  ! Protected module variables
  logical, public, protected :: tvs_oldFashionIRSeaEmiss ! use of the old Masuda IR emissivity instead of built-in RTTOV IREMIS
  logical, public, protected :: tvs_oldFashionIRLandEmiss ! use of the CERES instead of built-in RTTOV
  logical, public, protected :: tvs_useWaterFraction ! use of water fraction to compute hyperspectral IR surface emissivity
  logical, public, protected :: tvs_debug ! Logical key controlling statements to be  executed while debugging TOVS only
  real(8), public, protected, allocatable :: tvs_emissivity(:,:) ! Surface emissivities organized by profiles and channels
  integer, public, protected :: tvs_headerEnd ! header index of the last radiance observation
  integer, public, protected :: tvs_nsensors  ! Number of individual sensors.
  logical, public, protected :: tvs_mwAllskyAssim
  logical, public, protected :: tvs_computeJacobian ! Compute Jacobian for brightness temperature
  integer, public, protected :: tvs_platforms(tvs_maxNumberOfSensors)    ! RTTOV platform ID's (e.g., 1=NOAA; 2=DMSP; ...)
  integer, public, protected :: tvs_satellites(tvs_maxNumberOfSensors)   ! RTTOV satellite ID's (e.g., 1 to 16 for NOAA; ...)
  integer, public, protected :: tvs_instruments(tvs_maxNumberOfSensors)  ! RTTOVinstrument ID's (e.g., 3=AMSU-A; 4=AMSU-B; 6=SSMI; ...)
  integer, public, protected :: tvs_channelOffset(tvs_maxNumberOfSensors)! BURP to RTTOV channel mapping offset
  character(len=15), public, protected :: tvs_satelliteName(tvs_maxNumberOfSensors)
  character(len=15), public, protected :: tvs_instrumentName(tvs_maxNumberOfSensors)
  integer, public, protected, allocatable :: tvs_listSensors(:,:)     ! Sensor list
  integer, public, protected, allocatable :: tvs_nchan(:)             ! Number of channels per instrument (local)
  integer, public, protected, allocatable :: tvs_ichan(:,:)           ! List of channels per instrument (local)
  integer, public, protected, allocatable :: tvs_lsensor(:)           ! Sensor number for each profile
  integer, public, protected, allocatable :: tvs_nchanMpiGlobal(:)    ! Number of channels per instrument (global)
  logical, public, protected, allocatable :: tvs_isReallyPresent(:)   ! Logical flag to identify instruments really assimilated (local)
  logical, public, protected, allocatable :: tvs_isReallyPresentMpiGLobal(:) ! Logical flag to identify instruments really assimilated (global)
  type(surface_params), public, protected, allocatable :: tvs_surfaceParameters(:) ! surface parameters
  type(rttov_transmission), public, protected , allocatable :: tvs_transmission(:) ! transmittances all profiles for HIR quality control
  type(rttov_coefs), public, protected, allocatable         :: tvs_coefs(:)      ! rttov coefficients
  type(rttov_options), public, protected, allocatable       :: tvs_opts(:)       ! rttov options  
  type(rttov_radiance), public, protected, allocatable      :: tvs_radiance(:)   ! radiances organized by profile

  ! Private module parameters
  integer, parameter :: kslon=2160, kslat=1080 ! CERES file dimension in grid points
  integer, parameter :: maxsize = 100 ! Max number of instruments
  integer, parameter :: tvs_nlevels = 101  ! Maximum No. of RTTOV pressure levels including 'rttov top' at 0.005 hPa
  integer, parameter :: atmsTtChanNum(15) = (/  1,  2,  3,  4,  5,  6, 7, 8, 9, 10, &
                                               11, 12, 13, 14, 15 /)                  ! ATMS temperature channel numbers
  integer, parameter :: atmsHuChanNum(7)  = (/ 16, 17, 18, 19, 20, 21, 22 /)          ! ATMS humidity channel numbers
  real(8), parameter :: microg2kg   = 1.0d-9 ! units conversion from micrograms/kg to kg/kg

  ! Private module variables
  logical :: tvs_useO3FromTrials ! Determine if ozone model field (.true.) or climatology (.false.) is used
  logical :: tvs_useO3FromTrials_tl ! Determine if ozone if part of assimilation controal variable 
  integer :: tvs_maxNumberOfRadiances ! Max no of computed radiances for one sensor
  integer :: tvs_numMWInstrumUsingCLW
  integer :: tvs_numMWInstrumUsingHydrometeors
  logical :: tvs_mwInstrumUsingCLW_tl
  logical :: tvs_mwInstrumUsingHydrometeors_tl
  integer :: tvs_channelsUsingHydrometeors(tvs_maxNumberOfSensors,tvs_maxNumberOfChannels) ! List of channels using full set of hydromet variable, used in all-sky HU
  integer :: tvs_channelsUsingClw(tvs_maxNumberOfSensors,tvs_maxNumberOfChannels) ! List of channels using CLW, used in all-sky TT
  type(rttov_scatt_coef), allocatable    :: tvs_coef_scatt(:) ! rttovscatt coefficients
  type(rttov_options_scatt), allocatable :: tvs_opts_scatt(:) ! rttovscatt options
  integer, allocatable :: tvs_bodyIndexFromBtIndex(:,:)      ! Provides RTTOV bodyIndex in ObsSpaceData based on btIndex for each sensor
  integer, allocatable :: tvs_bodyIndexFromBtIndexScatt(:,:) ! Provides RTTOVScatt bodyIndex in ObsSpaceData based on btIndex for each sensor
  type(rttov_emis_atlas_data), allocatable :: tvs_atlas(:)   ! Emissivity atlases
  integer :: instrumentIdsUsingCLW(tvs_maxNumberOfSensors)
  integer :: instrumentIdsUsingHydrometeors(tvs_maxNumberOfSensors)
  logical :: tvs_copyCoefficientFileToRamDisk
  real(8) :: tvs_cloudScaleFactor     ! cloud scale factor used for rttov NL in all-sky assimilation
  real(8) :: tvs_cloudScaleFactor_tl  ! cloud scale factor used for rttov TL/AD in all-sky assimilation
  logical :: tvs_regLimitExtrap                    ! use RTTOV reg_limit_extrap option
  logical :: tvs_doAzimuthCorrection(tvs_maxNumberOfSensors)
  logical :: tvs_isAzimuthValid(tvs_maxNumberOfSensors)
  logical :: tvs_userDefinedDoAzimuthCorrection
  logical :: tvs_userDefinedIsAzimuthValid
  logical :: tvs_useSfcEmissObsSpace                   ! Logical key to use MW surface emissivity from ObsSpaceData
  logical :: tvs_irEmissAngularCorrection
  character(len=8) :: radiativeTransferCode            ! RadiativeTransferCode : TOVS radiation model used
  real(8), allocatable :: tvs_emissivityFromTrl(:,:) ! Surface emissivities extracted from trial
  type(rttov_chanprof), allocatable :: tvs_chanProf(:,:)
  type(rttov_chanprof), allocatable :: tvs_chanProfScatt(:,:)
  ! High resolution surface fields
  integer :: surfaceType(kslon,kslat)  
  real(8) :: waterFraction(kslon,kslat) 
  type(rttov_profile), target, allocatable :: tvs_profiles_nl(:)    ! all profiles on trial vertical coordinate for nl obs operator
  type(rttov_profile), target, allocatable :: tvs_profiles_tlad(:)  ! all profiles on increments vertical coordinates for linearized obs. operator
  type(rttov_profile_cloud), target, allocatable :: tvs_cld_profiles_nl(:)! rttov scatt cloud profiles on trial vertical coordinate
  type(rttov_profile_cloud), target, allocatable :: tvs_cld_profiles_tlad(:) ! rttov scatt cloud profiles on increment vertical coordinates

  ! Namelist variables:
  logical :: useMWEmissivityAtlas ! Flag to activate use of RTTOV built-in MW emissivity Atlases
  integer :: mWAtlasId            ! MW Atlas Id used when useMWEmissivityAtlas == .true. ; 1 TELSEM2, 2 CNRM atlas

contains

  !--------------------------------------------------------------------------
  ! tvs_allocateSurfaceParameters
  !--------------------------------------------------------------------------
  subroutine tvs_allocateSurfaceParameters()
    !
    ! :Purpose:  Allocate the global array `tvs_surfaceParameters`.
    !
    implicit none

    allocate(tvs_surfaceParameters(tvs_headerEnd))

  end subroutine tvs_allocateSurfaceParameters

  !--------------------------------------------------------------------------
  ! tvs_allocateEmissivity
  !--------------------------------------------------------------------------
  subroutine tvs_allocateEmissivity(maxChannelNumber)
    !
    ! :Purpose:  Allocate the global array `tvs_emissivity`.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: maxChannelNumber  ! maximum number of channels for all sensors 

    allocate(tvs_emissivity(maxChannelNumber,tvs_headerEnd))

  end subroutine tvs_allocateEmissivity

  !--------------------------------------------------------------------------
  !  validateRttovVariable
  !--------------------------------------------------------------------------
  subroutine validateRttovVariable(value, variableName, obsSpaceData, headerIndex, minimum_opt, maximum_opt) 
    !
    ! :Purpose:  check variable validity for RTTOV, 
    !            if invalid replace by acceptable value and reject
    !
    implicit none

    ! Arguments:
    real(8),           intent(inout) :: value           ! value of the RTTOV input variable to check for validity
    character(len=*),  intent(in)    :: variableName    ! Name of the checked variable for output in listings
    type(struct_obs),  intent(inout) :: obsSpaceData    ! obsSpaceData object
    integer,           intent(in)    :: headerIndex     ! observation index in obsSpaceData header table 
    real(8), optional, intent(in)    :: minimum_opt     ! Minimum acceptable value
    real(8), optional, intent(in)    :: maximum_opt     ! Maximum acceptable value

    ! Locals:
    real(8)                          :: minimum, maximum

    if (present(minimum_opt)) then
      minimum = minimum_opt
    else
      minimum = - huge(0.d0)
    end if
    
    if (present(maximum_opt)) then
      maximum = maximum_opt
    else
      maximum = huge(0.d0)
    end if
    
    if (value < minimum .or. value > maximum) then
      write(*,*) 'validateRttovVariable: !!! WARNING !!!'
      write(*,*) 'validateRttovVariable: INVALID ' // trim(variableName)
      write(*,*) 'validateRttovVariable: headerIndex ', headerIndex, " !"
      write(*,*) 'validateRttovVariable: ', value, ' should be between ', minimum, ' and ', maximum
      value = max(minimum, min(value, maximum) )
      write(*,*) 'validateRttovVariable: replaced with ', value, ' !'
      call rejectObs(obsSpaceData, headerIndex)
    end if

  end subroutine validateRttovVariable

  !--------------------------------------------------------------------------
  !  validateRttovProfile
  !--------------------------------------------------------------------------
  subroutine validateRttovProfile(profile, profileName, minimum, maximum, obsSpaceData, headerIndex) 
    !
    ! :Purpose:  check profile validity for RTTOV, 
    !            if invalid replace by acceptable value(s) and reject
    !
    implicit none

    ! Arguments:
    real(8),          intent(inout) :: profile(:)    ! Vertical profile of input variables to check for validity
    character(len=*), intent(in)    :: profileName   ! Name of the checked profile for output in listings
    real(8),          intent(in)    :: minimum       ! Minimum acceptable value
    real(8),          intent(in)    :: maximum       ! Maximum acceptable value
    type(struct_obs), intent(inout) :: obsSpaceData  ! obsSpaceData structure
    integer,          intent(in)    :: headerIndex   ! observation index in obsSpaceData header table

    ! Locals:
    logical                         :: ltest(size(profile))
    integer                         :: levelIndex
    
    ltest(:) = (profile(:) > maximum .or. profile(:) < minimum)
    
    if (any(ltest)) then
      write(*,*) 'validateRttovProfile: !!! WARNING !!!'
      write(*,*) 'validateRttovProfile: some INVALID ' // trim(profileName)
      write(*,*) 'validateRttovProfile: headerIndex ', headerIndex, " !"
      do levelIndex = 1, size(profile)
        if (ltest(levelIndex)) then
          write(*,*) 'validateRttovProfile: ', profile(levelIndex), &
              'should be between ', minimum, ' and ', maximum
          profile(levelIndex) = max(minimum, min(profile(levelIndex), maximum) )
          write(*,*) 'validateRttovProfile: replaced with ', profile(levelIndex) 
        end if
      end do
      call rejectObs(obsSpaceData, headerIndex)
    end if
    
  end subroutine validateRttovProfile
  
  !--------------------------------------------------------------------------
  ! rejectObs
  !--------------------------------------------------------------------------
  subroutine rejectObs(obsSpaceData, headerIndex)
    !
    ! :Purpose: reject all data corresponding to headerIndex
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData ! obsSpaceData structure
    integer,          intent(in)    :: headerIndex  ! observation index in obsSpaceData header table

    ! Locals:
    integer :: bodyIndex
    
    call obs_set_current_body_list(obsSpaceData, headerIndex)
    BODY:do 
      bodyIndex = obs_getBodyIndex(obsSpaceData)
      if (bodyIndex < 0) exit BODY
      call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
    end do BODY
    
  end subroutine rejectObs
  
  !--------------------------------------------------------------------------
  ! tvs_setupAlloc
  !--------------------------------------------------------------------------
  subroutine tvs_setupAlloc(obsSpaceData)
    !
    ! :Purpose: Memory allocation for the non linear radiative transfer model variables.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData  ! obsSpaceData structure
    
    ! Locals:
    integer :: satelliteCode, instrumentCode, iplatform, isat, instrum
    integer :: idatyp, sensorIndex
    integer :: channelNumber, nosensor, channelIndex
    integer :: errorStatus
    integer :: headerIndex, bodyIndex, taskIndex
    logical :: runObsOperatorWithClw
    logical :: runObsOperatorWithHydrometeors
    logical, allocatable :: logicalBuffer(:)
    integer, allocatable :: sensorTotalNumberOfProfiles(:)
    integer, allocatable :: ichanMpiGlobal(:,:)  ! List of channels per instrument
    character(len=32) :: hydroTableFileName
    character(len=300) :: filePrefix
    character(len=400) :: fileName
    character(len=512) :: fullNameWithPath, path
    logical :: fileExists
    integer :: extensionIndex
    integer,parameter :: nExtensions = 4
    character(len=4), parameter :: extensionList(nExtensions) = ['.bin', '.h5 ', '.H5 ', '.dat'] 

    if (tvs_nsensors == 0) return

    !  1. Determine the number of radiances to be assimilated.
    !      Construct a list of channels for each sensor.
    !      Construct a list of sensor number for each profile

    write(*,*) 'tvs_setupAlloc: Starting' 

    allocate(tvs_nchan(tvs_nsensors))
    allocate(tvs_ichan(tvs_maxNumberOfChannels,tvs_nsensors))
    allocate(tvs_lsensor(obs_numheader(obsSpaceData)))
    allocate(tvs_isReallyPresent(tvs_nsensors))
    allocate(tvs_nchanMpiGlobal(tvs_nsensors))
    allocate(tvs_isReallyPresentMpiGlobal(tvs_nsensors))
    allocate(sensorTotalNumberOfProfiles(tvs_nsensors))

    tvs_nchan(:) = 0 
    tvs_ichan(:,:) = 0
    sensorTotalNumberOfProfiles(:) = 0
    tvs_isReallyPresent(:) = .true.
    tvs_lsensor(:) = -1
    tvs_headerEnd = -1

    ! Loop over all header indices of the 'TO' family
    ! Set the header list & start at the beginning of the list
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)

      if (.not. tvs_isIdBurpTovs(idatyp)) then
        write(*,*) 'tvs_setupAlloc: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        call rejectObs(obsSpaceData, headerIndex)
        cycle HEADER   ! Proceed to the next headerIndex
      end if
      tvs_headerEnd = headerIndex
     
      !    Construct list of channels for each sensor:
      !          map burp satellite info to RTTOV platform and satellite.
      satelliteCode = obs_headElem_i(obsSpaceData,OBS_SAT,headerIndex)
      call tvs_mapSat(satelliteCode,iplatform,isat)
      if (iplatform == -1) then
        write(*,*) 'Unknown OBS_SAT !', satelliteCode
        call utl_abort('tvs_setupAlloc')
      end if
      !    map burp instrument info to RTTOV instrument.
      instrumentCode = obs_headElem_i(obsSpaceData,OBS_INS,headerIndex)
      call tvs_mapInstrum(instrumentCode,instrum)
      if (instrum == -1) then
        write(*,*) 'Unknown OBS_INS !', instrumentCode
        call utl_abort('tvs_setupAlloc')
      end if
      !    find sensor number for this obs.
      nosensor = 0
      do sensorIndex = 1, tvs_nsensors
        if (iplatform == tvs_platforms  (sensorIndex) .and. &
            isat      == tvs_satellites (sensorIndex) .and. &
            instrum   == tvs_instruments(sensorIndex)      ) then
          nosensor = sensorIndex
          exit
        end if
      end do

      if (nosensor > 0) then
        tvs_lsensor(headerIndex) = nosensor
        sensorTotalNumberOfProfiles(nosensor) = sensorTotalNumberOfProfiles(nosensor) + 1
      else
        write(*,*) ' tvs_setupAlloc: Warning Invalid Sensor ', iplatform, isat, instrum, ' skipping ...'
      end if

      ! Loop over all body indices (still in the 'TO' family)
      ! Set the body list & start at the beginning of the list
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY
        if (nosensor > 0) then
          if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) then
            call tvs_getChannelNumIndexFromPPP( obsSpaceData, headerIndex, bodyIndex, &
                                                channelNumber, channelIndex )
            if (channelIndex == 0) then
              tvs_nchan(nosensor) = tvs_nchan(nosensor) + 1
              tvs_ichan(tvs_nchan(nosensor),nosensor) = channelNumber
            end if
            
            if (tvs_debug .and. mmpi_myid == 0 .and. &
                 trim(tvs_instrumentName(nosensor)) == 'AMSUA') then
              write(*,*) 'test channelNumber:', headerIndex, bodyIndex, nosensor, &
                   tvs_satelliteName(nosensor), channelNumber, channelIndex
            end if
          end if
        else           
          ! set to notAssimilated if instrument not in NAMTOV namelist
          call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
        end if
      end do BODY
    end do HEADER    
   
    tvs_maxNumberOfRadiances = maxval(tvs_nchan(:) * sensorTotalNumberOfProfiles(:))
    deallocate(sensorTotalNumberOfProfiles)

    if (.not. tvs_userDefinedDoAzimuthCorrection) then 
      ! tvs_doAzimuthCorrection user defined values will be overwriten by the old default values 
      do sensorIndex = 1, tvs_nsensors
        tvs_doAzimuthCorrection(sensorIndex) = ( tvs_platforms(sensorIndex) /= platform_id_eos .and. &
             ( tvs_instruments(sensorIndex) == inst_id_amsua .or. tvs_instruments(sensorIndex) == inst_id_mhs )  )     
      end do
      if (mmpi_myId == 0) write(*,*) ' tvs_setupAlloc: Warning tvs_doAzimuthCorrection user defined values overwriten by the old default values'
    end if

    if (.not. tvs_userDefinedIsAzimuthValid) then 
      ! tvs_isAzimuthValid  user defined values will be overwriten by the current default values 
      do sensorIndex = 1, tvs_nsensors
        tvs_isAzimuthValid(sensorIndex) = .not. ( tvs_isInstrumGeostationary(tvs_instruments(sensorIndex)) )
      end do
      if (mmpi_myId == 0) write(*,*) ' tvs_setupAlloc: Warning tvs_isAzimuthValid user defined values overwriten by the current default values'
    end if

    if (mmpi_myId == 0) then
      write(*,*) ' tvs_setupAlloc: platform satellite id tvs_doAzimuthCorrection tvs_isAzimuthValid'
      do sensorIndex = 1, tvs_nsensors
        write(*,'(18x,a,1x,a,1x,i2,1x,L1,10x,L1)') inst_name(tvs_instruments(sensorIndex)), &
             platform_name(tvs_platforms(sensorIndex)), tvs_satellites(sensorIndex), &
             tvs_doAzimuthCorrection(sensorIndex), tvs_isAzimuthValid(sensorIndex)
      end do
    end if
    
    ! Sort list of channels in ascending order.Also force at least one channel, if none are found.
    do sensorIndex = 1, tvs_nsensors
      call isort(tvs_ichan(:,sensorIndex),tvs_nchan(sensorIndex))
      if (tvs_nchan(sensorIndex) == 0) then
        tvs_isReallyPresent(sensorIndex) = .false.
        tvs_nchan(sensorIndex) = 1
        tvs_ichan(1,sensorIndex) = 1
      end if
    end do

    write(*,*) ' tvs_setupAlloc: tvs_headerEnd = ', tvs_headerEnd
    
    allocate(ichanMpiGlobal(tvs_maxNumberOfChannels,tvs_nsensors))
    do sensorIndex = 1, tvs_nsensors
      call tvs_getCommonChannelSet(tvs_ichan(:,sensorIndex),tvs_nchanMpiGlobal(sensorIndex), ichanMpiGlobal(:,sensorIndex))
      write(*,*) 'sensorIndex,tvs_nchan(sensorIndex),tvs_nchanMpiGlobal(sensorIndex)', sensorIndex, tvs_nchan(sensorIndex),tvs_nchanMpiGlobal(sensorIndex)
    end do
    deallocate(ichanMpiGlobal)
    
    if (mmpi_myid == 0) then
      allocate(logicalBuffer(mmpi_nprocs))
    else
      allocate(logicalBuffer(1))
    end if
    
    do sensorIndex = 1, tvs_nsensors
      call rpn_comm_gather(tvs_isReallyPresent ( sensorIndex ) , 1, 'MPI_LOGICAL', logicalBuffer, 1,'MPI_LOGICAL', 0, 'GRID', errorStatus )
      if (mmpi_myid == 0) then
        tvs_isReallyPresentMpiGlobal ( sensorIndex ) =.false.
        do taskIndex = 1, mmpi_nprocs
          tvs_isReallyPresentMpiGlobal ( sensorIndex ) =  tvs_isReallyPresentMpiGlobal ( sensorIndex ) .or. logicalBuffer(taskIndex)
        end do
      end if
      call rpn_comm_bcast(tvs_isReallyPresentMpiGlobal ( sensorIndex ), 1, 'MPI_LOGICAL', 0, 'GRID', errorStatus )
    end do
    
    deallocate(logicalBuffer)
    
    if (tvs_debug .and. mmpi_myid == 0) then
      do sensorIndex = 1, tvs_nsensors
        write(*,*) 'sensorIndex, tvs_instrumentName(sensorIndex), tvs_satelliteName(sensorIndex)'
        write(*,*) sensorIndex, tvs_instrumentName(sensorIndex), tvs_satelliteName(sensorIndex)
        write(*,*) 'tvs_channelOffset(sensorIndex), tvs_nchan(sensorIndex)'
        write(*,*) tvs_channelOffset(sensorIndex), tvs_nchan(sensorIndex)
        write(*,*) 'tvs_ichan(1:tvs_nchan(sensorIndex),sensorIndex)'
        write(*,*) tvs_ichan(1:tvs_nchan(sensorIndex),sensorIndex)
        write(*,*) 
      end do
    end if

    !  3. Initialize TOVS radiance transfer model

    if (radiativeTransferCode == 'RTTOV') then

      write(*,'(//,10x,A)') '-rttov_setup: initializing the TOVS radiative transfer model'

      allocate(tvs_coefs(tvs_nsensors))
      allocate(tvs_listSensors(3,tvs_nsensors))
      allocate(tvs_opts(tvs_nsensors))
      if (tvs_numMWInstrumUsingHydrometeors > 0) then
        allocate(tvs_opts_scatt(tvs_nsensors))
        allocate(tvs_coef_scatt(tvs_nsensors))
      end if

      do sensorIndex=1, tvs_nsensors
        tvs_listSensors(1,sensorIndex) = tvs_platforms  (sensorIndex)
        tvs_listSensors(2,sensorIndex) = tvs_satellites (sensorIndex)
        tvs_listSensors(3,sensorIndex) = tvs_instruments(sensorIndex)

        runObsOperatorWithClw = (tvs_numMWInstrumUsingCLW /= 0 .and. &
                                 tvs_isInstrumUsingCLW(tvs_instruments(sensorIndex)))
        runObsOperatorWithHydrometeors = (tvs_numMWInstrumUsingHydrometeors /= 0 .and. &
                                          tvs_isInstrumUsingHydrometeors(tvs_instruments(sensorIndex)))

        !< General configuration options
        tvs_opts(sensorIndex) % config % apply_reg_limits = .true. ! if true application of profiles limits
        tvs_opts(sensorIndex) % config % verbose = .false. ! verbose output
        tvs_opts(sensorIndex) % config % do_checkinput = .true. ! to check if input profiles are within absolute and regression limits
        tvs_opts(sensorIndex) % config % fix_hgpl = .false. ! for backward compatibility with RTTOV-12 should be changed later
        !< General RT options
        tvs_opts(sensorIndex) % rt_all % switchrad = .true.  ! to use brightness temperature (true) or radiance (false) units in AD routine
        tvs_opts(sensorIndex) % rt_all % use_q2m = .false.   ! if true use of surface humidity (false for compatibility with the way rttov 8.7 was compiled)
        tvs_opts(sensorIndex) % rt_all % addrefrac = .true.  ! to account for atmospheric refraction
        tvs_opts(sensorIndex) % rt_all % dtau_test = .true.  ! for backward compatibility with RTTOV-12 may be changed later
        tvs_opts(sensorIndex) % rt_all % use_t2m_opdep = .false. ! for backward compatibility with RTTOV-12 may be changed later
        !< VIS/IR RT options
        tvs_opts(sensorIndex) % rt_ir % addsolar = .false.  ! to model solar component in the near IR (2000 cm-1 et plus)
        tvs_opts(sensorIndex) % rt_ir % addaerosl = .false. ! to account for scattering due to aerosols
        tvs_opts(sensorIndex) % rt_ir % addclouds = .false. ! to account for scattering due to clouds
        tvs_opts(sensorIndex) % rt_ir % ir_sea_emis_model = 2 ! ISEM (ir_sea_emis_model 1) useful for GEORAD
                                                              ! 2 selects IREMIS which is more sophisticated 
        tvs_opts(sensorIndex) % rt_ir % pc % ipcreg = -1         ! index of the required PC predictors... to see later
        tvs_opts(sensorIndex) % rt_ir % pc % addpc = .false.     ! to carry out principal component calculations 
        tvs_opts(sensorIndex) % rt_ir % pc % addradrec = .false. ! to reconstruct radiances from principal components
        !< MW RT options
        tvs_opts(sensorIndex) % rt_mw % clw_data = tvs_isInstrumUsingCLW(tvs_instruments(sensorIndex)) ! disponibilite du profil d'eau liquide
        tvs_opts(sensorIndex) % rt_mw % fastem_version = 6  ! use fastem version 6 microwave sea surface emissivity model (1-6)
        tvs_opts(sensorIndex) % rt_mw % clw_scheme = 1 ! default and recommended is 2 just for backward compatibility
        !< Interpolation options
        tvs_opts(sensorIndex) % interpolation % addinterp = .true. ! use of internal profile interpolator (rt calculation on model levels)
        tvs_opts(sensorIndex) % interpolation % lgradp = .true.    ! allow tl/ad of user pressure levels
        tvs_opts(sensorIndex) % interpolation % interp_mode = interp_rochon_loglinear_wfn ! see table 9 page 37 of RTTOV 12.1 users guide
        tvs_opts(sensorIndex) % interpolation % reg_limit_extrap = tvs_regLimitExtrap 

        tvs_opts(sensorIndex) % rt_all % co2_data = .false.
        tvs_opts(sensorIndex) % rt_all % n2o_data = .false.
        tvs_opts(sensorIndex) % rt_all % co_data  = .false.
        tvs_opts(sensorIndex) % rt_all % ch4_data = .false.

        if (runObsOperatorWithHydrometeors) then
          tvs_opts_scatt(sensorIndex) % interp_mode =  tvs_opts(sensorIndex) % interpolation % interp_mode ! Set interpolation method
          tvs_opts_scatt(sensorIndex) % reg_limit_extrap = tvs_regLimitExtrap 
          tvs_opts_scatt(sensorIndex) % fastem_version = tvs_opts(sensorIndex) % rt_mw % fastem_version  
          tvs_opts_scatt(sensorIndex) % supply_foam_fraction = .false.
          tvs_opts_scatt(sensorIndex) % use_t2m_opdep = tvs_opts(sensorIndex) % rt_all % use_t2m_opdep
          tvs_opts_scatt(sensorIndex) % use_q2m = tvs_opts(sensorIndex) % rt_all % use_q2m
          tvs_opts_scatt(sensorIndex) % lgradp = .true.
          tvs_opts_scatt(sensorIndex) % lusercfrac = .false. !< Switch to enable user-specified effective cloud fraction ??
          tvs_opts_scatt(sensorIndex) % config % do_checkinput = tvs_opts(sensorIndex) % config % do_checkinput
          tvs_opts_scatt(sensorIndex) % config % apply_reg_limits = tvs_opts(sensorIndex) % config % apply_reg_limits
          tvs_opts_scatt(sensorIndex) % config % verbose = .true.
          tvs_opts_scatt(sensorIndex) % config % fix_hgpl= tvs_opts(sensorIndex) % config % fix_hgpl
          ! other option may be considered:
          !real(jprb)    :: cc_threshold          = 1.E-3_jprb    !< Minimum effective cloud fraction threshold to consider scattering
          !real(jprb)    :: ice_polarisation      = 1.40_jprb     !< Polarised scattering factor for ice hydrometeors (<0 = no polarisation)
          !logical(jplm) :: ozone_data            = .false.       !< Switch to enable input of O3 profile
                                                                  ! because standard RTTOV coefficients in the MW have no ozone sensitivity
          !logical(jplm) :: rad_down_lin_tau      = .true.        !< Linear-in-tau or layer-mean for downwelling radiances
          !logical(jplm) :: hydro_cfrac_tlad      = .true.        !< Switch for hydrometeor TL/AD sensitivity to effective cfrac
          !logical(jplm) :: zero_hydro_tlad       = .false.       !< Switch for hydrometeor TL/AD sensitivity in layers with zero
                                                                  !   hydrometeor concentration
        end if
        
        errorStatus = errorStatus_success
        call utl_tmg_start(16,'----RttovSetup')
        write(*,*) ' sensorIndex,tvs_nchan(sensorIndex)', sensorIndex, tvs_nchan(sensorIndex)
        
        path = './'
        fullNameWithPath = 'NOTFOUND'
        if (tvs_copyCoefficientFileToRamdisk .and. ram_getRamDiskDir() /= ' ') then
          call rttov_coeffname(errorStatus, tvs_listSensors(:,sensorIndex), 'rtcoef', filePrefix)
          do extensionIndex = 1, nExtensions
            fileName = './' // trim(filePrefix) // trim(extensionList(extensionIndex))
            inquire(file=fileName, exist=fileExists)
            if (fileExists) then
              fullNameWithPath = ram_fullWorkingPath(fileName) ! copy to ramdisk and return path on ramdisk
              path = trim(ram_getRamDiskDir())
              exit
            end if
          end do
          if (trim(fullNameWithPath) == 'NOTFOUND') then
            call utl_abort('tvs_setupAlloc: unable to find coefficient file starting with ' // trim(filePrefix))
          end if
        end if
        
        write(*,*) 'tvs_setupAlloc: calling rttov_read_coefs with path= ', trim(path)
        call rttov_read_coefs(                              &
            errorStatus,                                    &! out
            tvs_coefs(sensorIndex),                         &
            tvs_opts(sensorIndex),                          &
            instrument=tvs_listSensors(:,sensorIndex),      &! in
            path=path,                                      &
            channels=tvs_ichan(1:tvs_nchan(sensorIndex),sensorIndex) ) ! in option
        
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'rttov_read_coefs: fatal error reading coefficients',errorStatus,sensorIndex,tvs_listSensors(1:3,sensorIndex)
          call utl_abort('tvs_setupAlloc')
        end if

        if (tvs_copyCoefficientFileToRamdisk .and. ram_getRamDiskDir() /= ' ') then
          errorStatus = ram_remove(fullNameWithPath)
        end if
       
        if (runObsOperatorWithHydrometeors) then
          hydrotableFileName = 'hydrotable_' // trim(platform_name(tvs_platforms(sensorIndex))) // '_' // &
               trim(inst_name(tvs_instruments(sensorIndex))) // '.dat'
          call rttov_read_scattcoeffs(errorStatus, tvs_opts_scatt(sensorIndex), tvs_coefs(sensorIndex), &
               tvs_coef_scatt(sensorIndex), file_coef=hydrotableFilename)
          if (errorStatus /= errorStatus_success) then
            write(*,*) 'rttov_read_scattcoeffs: fatal error reading RTTOV-SCATT coefficients', hydrotableFileName
            call utl_abort('tvs_setupAlloc')
          end if
        end if
        call utl_tmg_stop(16)

        tvs_opts(sensorIndex) % rt_all % ozone_data = ( tvs_coefs(sensorIndex) % coef % nozone > 0 ) ! profil d'ozone disponible

        ! Ensure the options and coefficients are consistent
        call rttov_user_options_checkinput(errorStatus, tvs_opts(sensorIndex), tvs_coefs(sensorIndex))
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'error in rttov options', errorStatus
          call utl_abort('tvs_setupAlloc')
        end if
      end do
      
      !   4. Memory allocations for radiative tranfer model variables

      ! Radiance by profile
      allocate(tvs_radiance(tvs_headerEnd))
      do headerIndex = 1, tvs_headerEnd
        sensorIndex = tvs_lsensor(headerIndex)
        if (sensorIndex <= 0) cycle
        ! allocate BT equivalent to total direct, tl and ad radiance output
        allocate(tvs_radiance(headerIndex) % bt (tvs_nchan(sensorIndex)))
        tvs_radiance(headerIndex) % bt (:) = 0.d0
        nullify (tvs_radiance(headerIndex) % clear)
      end do

    end if !if (radiativeTransferCode == 'RTTOV')
  
    write(*,*) 'Leaving tvs_setupAlloc'

  end subroutine tvs_setupAlloc

  !--------------------------------------------------------------------------
  ! tvs_getProfile
  !--------------------------------------------------------------------------
  subroutine tvs_getProfile(profiles, profileType, cld_profiles_opt)
    !
    ! :Purpose: sets profiles as a pointer of type rttov_profile
    !           based on profileType equal to nl or tlad. 
    ! 
    implicit none

    ! Arguments:
    type(rttov_profile),       pointer,           intent(inout) :: profiles(:)         ! pointer to profile structures
    character(len=*),                             intent(in)    :: profileType         ! profile type ('nl' or 'tlad')
    type(rttov_profile_cloud), pointer, optional, intent(inout) :: cld_profiles_opt(:) ! Do we also need cloud profiles (for all sky assimilation) ?

    select case(trim( profileType))
      case('nl')
        profiles => tvs_profiles_nl
        if (present(cld_profiles_opt)) cld_profiles_opt => tvs_cld_profiles_nl
      case('tlad')
        profiles => tvs_profiles_tlad
        if (present(cld_profiles_opt)) cld_profiles_opt => tvs_cld_profiles_tlad
      case default
        call utl_abort('tvs_getProfile: invalid profileType ' // profileType)
    end select

  end subroutine tvs_getProfile

  !--------------------------------------------------------------------------
  ! tvs_allocTransmission
  !--------------------------------------------------------------------------
  subroutine tvs_allocTransmission(nlv_T)
    !
    ! :Purpose: Allocate the global rttov transmission structure used
    !           when this is needed for some purpose (e.g. used in 
    !           LETKF to determine peak pressure level of each radiance
    !           channel for vertical localization).
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: nlv_T ! number of model thermodynamical vertical levels

    ! Locals:
    integer :: headerIndex, sensorIndex, nchannels

    if (allocated(tvs_transmission)) deallocate(tvs_transmission)
    allocate(tvs_transmission(tvs_headerEnd))

    do headerIndex = 1, tvs_headerEnd
      sensorIndex = tvs_lsensor(headerIndex)
      if (sensorIndex <= 0) cycle
      nchannels = tvs_nchan(sensorIndex)
      ! allocate transmittance from surface and from pressure levels
      allocate(tvs_transmission(headerIndex) % tau_total(nchannels))
      allocate(tvs_transmission(headerIndex) % tau_levels(nlv_T,nchannels))
    end do

  end subroutine tvs_allocTransmission

  !--------------------------------------------------------------------------
  ! tvs_setup
  !--------------------------------------------------------------------------
  subroutine tvs_setup
    !
    ! :Purpose: to read namelist NAMTOV, initialize the observation error covariance and setup RTTOV.
    !
    implicit none

    ! Locals:
    integer :: sensorIndex, ierr
    integer :: instrumentIndex, numMWInstrumToUseCLW, numMWInstrumToUseHydrometeors

    ! Namelist variables: (local)
    character(len=8)  :: crtmodl ! For now, must equal 'RTTOV'
    integer :: nsensors ! MUST NOT BE INCLUDED IN NAMELIST!
    character(len=15) :: csatid(tvs_maxNumberOfSensors)        ! List of satellite names
    character(len=15) :: cinstrumentid(tvs_maxNumberOfSensors) ! List of incstrument names
    logical :: ldbgtov  ! Choose to print simulated and observed Tb to listing
    logical :: useO3FromTrials ! Choose to use ozone model fields (otherwise climatology)
    logical :: useO3FromTrials_tl ! Choose to get contributions of radiance to the ozone analysis
    logical :: regLimitExtrap ! Choose to use RTTOV reg_limit_extrap option
    logical :: doAzimuthCorrection(tvs_maxNumberOfSensors) ! Choose to apply correction to azimuth angle
    logical :: isAzimuthValid(tvs_maxNumberOfSensors) ! Indicate if azimuth angle is valid
    logical :: userDefinedDoAzimuthCorrection ! Indicate if user defined azimuth correction is to be used
    logical :: userDefinedIsAzimuthValid ! Indicate if user defined azimuth angle is valid
    logical :: copyCoefficientFileToRamDisk ! Copy RTTOV coefficient files to ramdisk
    logical :: mwInstrumUsingCLW_tl ! Choose to use CLW increment in TL/AD if exists as state variable
    logical :: mwInstrumUsingHydrometeors_tl ! Choose to all hydomet variables in TL/AD if exist as state variables
    real(8) :: cloudScaleFactor     ! Scale factor applied in rttov NL to model produced clouds to account for bias
    real(8) :: cloudScaleFactor_tl  ! Scale factor applied in rttov TL/AD to cloud increments
    character(len=15) :: instrumentNamesUsingCLW(tvs_maxNumberOfSensors) ! List of inst names using CLW
    character(len=15) :: instrumentNamesUsingHydrometeors(tvs_maxNumberOfSensors) ! List of inst name using full set of hydromet variables
    integer :: channelsUsingHydrometeors(tvs_maxNumberOfSensors,tvs_maxNumberOfChannels) ! List of channels using full set of hydromet variables, used in all-sky HU
    integer :: channelsUsingClw(tvs_maxNumberOfSensors,tvs_maxNumberOfChannels) ! List of channels using CLW, used in all-sky TT
    logical :: mwAllskyAssim ! High-level key to activate all-sky treatment of MW radiances
    logical :: computeJacobian !Choose to compute Jacobian for brightness temperature
    logical :: oldFashionIRSeaEmiss ! if .true. use of the old Masuda HIRS resolution IR emissivity instead of built-in RTTOV IREMIS
    logical :: oldFashionIRLandEmiss ! if .true. use of the old CERES based land emissivity instead of Borbas 2xxx
    logical :: irEmissAngularCorrection ! apply emissivity angular correction for hyperspectral IR when using RTTOV built-in land emissivity atlas
    logical :: useWaterFraction ! use of water fraction to compute hyperspectral IR surface emissivity when using RTTOV built-in land emissivity atlas
    
    namelist /NAMTOV/ nsensors, csatid, cinstrumentid
    namelist /NAMTOV/ ldbgtov,useO3FromTrials, useO3FromTrials_tl
    namelist /NAMTOV/ crtmodl
    namelist /NAMTOV/ useMWEmissivityAtlas, mWAtlasId
    namelist /NAMTOV/ mwInstrumUsingCLW_tl, instrumentNamesUsingCLW
    namelist /NAMTOV/ mwInstrumUsingHydrometeors_tl, instrumentNamesUsingHydrometeors
    namelist /NAMTOV/ channelsUsingHydrometeors, channelsUsingClw
    namelist /NAMTOV/ regLimitExtrap, doAzimuthCorrection, userDefinedDoAzimuthCorrection
    namelist /NAMTOV/ isAzimuthValid, userDefinedIsAzimuthValid
    namelist /NAMTOV/ cloudScaleFactor, cloudScaleFactor_tl 
    namelist /NAMTOV/ mwAllskyAssim, copyCoefficientFileToRamDisk, computeJacobian
    namelist /NAMTOV/ oldFashionIRSeaEmiss, oldFashionIRLandEmiss, irEmissAngularCorrection
    namelist /NAMTOV/ useWaterFraction
    
    ! Use MW surface emissivity from ObsSpaceData
    tvs_useSfcEmissObsSpace = .false.
    tvs_headerEnd = -1

    ! return if the NAMTOV does not exist
    if (.not. utl_isNamelistPresent('NAMTOV','./flnml')) then
      write(*,*)
      write(*,*) 'tvs_setup: Namelist block NAMTOV is missing in the namelist.'
      write(*,*) '           Skipping tvs_setup.'
      return
    end if
 
    !   1.1 Default values for namelist variables
    nsensors = MPC_missingValue_INT
    csatid(:) = '***UNDEFINED***'
    cinstrumentid(:) = '***UNDEFINED***'
    doAzimuthCorrection(:) = .false.
    isAzimuthValid(:) = .false.
    ldbgtov = .false.
    useO3FromTrials = .false.
    useO3FromTrials_tl = .false.
    userDefinedDoAzimuthCorrection = .false.
    userDefinedIsAzimuthValid = .false.
    crtmodl = 'RTTOV'
    useMWEmissivityAtlas = .false.
    mWAtlasId = 1 !Default to TELSEM-2
    mwInstrumUsingCLW_tl = .false.
    mwInstrumUsingHydrometeors_tl = .false.
    instrumentNamesUsingCLW(:) = '***UNDEFINED***'
    channelsUsingClw(:,:) = -1
    instrumentNamesUsingHydrometeors(:) = '***UNDEFINED***'
    channelsUsingHydrometeors(:,:) = -1
    regLimitExtrap = .true.
    cloudScaleFactor = 0.5D0
    cloudScaleFactor_tl = 1.0D0
    mwAllskyAssim = .false.
    copyCoefficientFileToRamDisk = .true.
    computeJacobian = .false.
    oldFashionIRSeaEmiss = .true.
    oldFashionIRLandEmiss = .true.
    irEmissAngularCorrection = .false.
    useWaterFraction = .true.
    
    !   1.2 Read the NAMELIST NAMTOV to modify them
    call utl_tmg_start(181,'low-level--readNML')
    read(utl_flnml, nml=namtov, iostat=ierr)
    if (ierr /= 0) call utl_abort('tvs_setup: Error reading namelist NAMTOV')

    if (mmpi_myid == 0) write(*,nml=namtov)
    call utl_tmg_stop(181)

    !  1.3 Transfer namelist variables to module variables
    if (nsensors /= MPC_missingValue_INT) then
      call utl_abort('tvs_setup: check namelist section NAMTOV; nsensors should be removed as it is' // &
          ' now computed by Midas from cinstrumentid and csatid arrays')
    end if
    
    tvs_nsensors = 0
    sensor_loop:do sensorIndex = 1, tvs_maxNumberOfSensors
      if (cinstrumentid(sensorIndex) /= "***UNDEFINED***" .and. &
          csatid(sensorIndex) /= "***UNDEFINED***" ) then
        tvs_nsensors = tvs_nsensors + 1
      else
        exit sensor_loop
      end if
    end do sensor_loop
    
    tvs_debug = ldbgtov
    radiativeTransferCode = crtmodl
    tvs_useO3FromTrials = useO3FromTrials
    tvs_useO3FromTrials_tl = useO3FromTrials_tl
    tvs_instrumentName(:) = cinstrumentid(:)
    tvs_satelliteName(:) = csatid(:)
    tvs_mwInstrumUsingCLW_tl = mwInstrumUsingCLW_tl
    tvs_mwInstrumUsingHydrometeors_tl = mwInstrumUsingHydrometeors_tl
    tvs_regLimitExtrap = regLimitExtrap
    tvs_userDefinedDoAzimuthCorrection = userDefinedDoAzimuthCorrection
    tvs_userDefinedIsAzimuthValid = userDefinedIsAzimuthValid
    tvs_doAzimuthCorrection(:) =  doAzimuthCorrection(:)
    tvs_isAzimuthValid(:) =  isAzimuthValid(:)
    tvs_cloudScaleFactor = cloudScaleFactor 
    tvs_cloudScaleFactor_tl = cloudScaleFactor_tl 
    tvs_mwAllskyAssim = mwAllskyAssim
    tvs_copyCoefficientFileToRamDisk = copyCoefficientFileToRamDisk
    tvs_computeJacobian = computeJacobian
    tvs_channelsUsingHydrometeors(:,:) = channelsUsingHydrometeors(:,:)
    tvs_channelsUsingClw(:,:) = channelsUsingClw(:,:)
    tvs_oldFashionIRSeaEmiss = oldFashionIRSeaEmiss
    tvs_oldFashionIRLandEmiss = oldFashionIRLandEmiss
    tvs_irEmissAngularCorrection = irEmissAngularCorrection
    tvs_useWaterFraction = useWaterFraction
    
    !  1.4 Validate namelist values
    if (tvs_nsensors == 0) then
      if (mmpi_myid == 0) then 
        write(*,*) ' ====================================================='
        write(*,*) ' tvs_setup: Number of sensors is zero, skipping setup'
        write(*,*) ' ====================================================='
      end if
      return
    end if

    if (radiativeTransferCode /= 'RTTOV') then
      write(*,'(A)') ' Invalid radiation model name'
      call utl_abort('tvs_setup')
    end if

    if (tvs_nsensors > tvs_maxNumberOfSensors) then
      write(*,'(A)') ' Number of sensors (tvs_nsensors) is greater than maximum allowed (tvs_maxNumberOfSensors)'
      call utl_abort('tvs_setup')
    end if

    !  1.5 Print the content of this NAMELIST
    if (mmpi_myid == 0) then
      write(*,'(A)') 
      write(*,'(3X,A)') '- Parameters used for TOVS processing (read in NAMTOV)'
      write(*,'(3X,A)') '  ----------------------------------------------------'
      write(*,'(6X,A,2X,L1)') 'TOVS debug                           : ', tvs_debug
      write(*,'(6X,A,2X,L1)') 'Use of UW IR land emissivity atlases : ', .not. tvs_oldFashionIRLandEmiss
      write(*,'(6X,A,2X,L1)') 'Use of MW land emissivity atlases    : ', useMWEmissivityAtlas
      if (useMWEmissivityAtlas) then
        write(*,'(6X,A,2X,I1)') 'MW atlas Id                          : ', mWAtlasId
      end if
      write(*,'(6X,A,2X,L1)') 'Use of reg_limit_extrap              : ', regLimitExtrap
      write(*,'(6X,A,2X,A)')  'Radiative transfer model             : ', radiativeTransferCode
      write(*,'(6X,A,2X,I3)') 'Number of sensors                    : ', tvs_nsensors
      write(*,"(6X,'Satellite ids          : ',10A10)") (tvs_satelliteName(sensorIndex), sensorIndex=1,tvs_nsensors)
      write(*,"(6X,'Instrument ids         : ',10A10)") (tvs_instrumentName(sensorIndex), sensorIndex=1,tvs_nsensors)
      write(*,'(A)') 
      write(*,'(A)') 
      write(*,'(A)') 
      write(*,'(3X,A)') '- Reading and initialization in preparation to the TOVS processing'
      write(*,'(5X,A)') '----------------------------------------------------------------'
    end if

    !  1.6 Set up platform, satellite, instrument and channel mapping

    call sensors()

    ! Get the name and number of instruments to use CLW
    instrumentIdsUsingCLW(:) = -1
    do instrumentIndex = 1, tvs_nsensors
      instrumentIdsUsingCLW(instrumentIndex) = tvs_getInstrumentId(instrumentNamesUsingCLW(instrumentIndex))
      if (instrumentNamesUsingCLW(instrumentIndex) /= '***UNDEFINED***') then
        if (instrumentIdsUsingCLW(instrumentIndex) == -1) then
          write(*,*) instrumentIndex, instrumentNamesUsingCLW(instrumentIndex)
          call utl_abort('tvs_setup: Unknown instrument name to use CLW')
        end if
      else
        numMWInstrumToUseCLW = instrumentIndex - 1
        exit
      end if
    end do

    ! Get the name and number of instruments to use hydrometeors
    instrumentIdsUsingHydrometeors(:) = -1
    do instrumentIndex = 1, tvs_nsensors
      instrumentIdsUsingHydrometeors(instrumentIndex) = &
                  tvs_getInstrumentId(instrumentNamesUsingHydrometeors(instrumentIndex))
      if (instrumentNamesUsingHydrometeors(instrumentIndex) /= '***UNDEFINED***') then
        if (instrumentIdsUsingHydrometeors(instrumentIndex) == -1) then
          write(*,*) instrumentIndex, instrumentNamesUsingHydrometeors(instrumentIndex)
          call utl_abort('tvs_setup: Unknown instrument name to use hydrometeors')
        end if
      else
        numMWInstrumToUseHydrometeors = instrumentIndex - 1
        exit
      end if
    end do

    ! check instrument is either using CLW or hydrometeors for non-ATMS instruments
    do instrumentIndex = 1, numMWInstrumToUseHydrometeors
      if (numMWInstrumToUseCLW == 0) exit

      if (any(instrumentIdsUsingCLW(1:numMWInstrumToUseCLW) == &
              instrumentIdsUsingHydrometeors(instrumentIndex))) then
        if (instrumentIdsUsingHydrometeors(instrumentIndex) /= tvs_getInstrumentId('atms')) then
          write(*,*) instrumentIndex, instrumentNamesUsingHydrometeors(instrumentIndex)
          call utl_abort('tvs_setup: all-sky TtHu for this intrument is not assimilated yet')
        end if
      end if
    end do

    tvs_numMWInstrumUsingCLW = numMWInstrumToUseCLW
    tvs_numMWInstrumUsingHydrometeors = numMWInstrumToUsehydrometeors

    if (mmpi_myid == 0) then
      write(*,*) 'tvs_setup: Instrument IDs to use CLW: ', instrumentIdsUsingCLW(1:numMWInstrumToUseCLW)
      write(*,*) 'tvs_setup: Number of instruments to use CLW: ', numMWInstrumToUseCLW

      write(*,*) 'tvs_setup: Instrument IDs to use hydrometeors: ', &
                             instrumentIdsUsingHydrometeors(1:numMWInstrumToUseHydrometeors)
      write(*,*) 'tvs_setup: Number of instruments to use hydrometeors: ', &
                             numMWInstrumToUseHydrometeors
    end if

  end subroutine tvs_setup

  !--------------------------------------------------------------------------
  ! tvs_cleanup
  !--------------------------------------------------------------------------
  subroutine tvs_cleanup
    !
    ! :Purpose: release memory used by RTTOV.
    !
    implicit none

    ! Locals:
    integer :: allocStatus
    integer :: sensorIndex,headerIndex,nlv_T

    write(*,*) 'tvs_cleanup: Starting'

    if (radiativeTransferCode == 'RTTOV') then

      !___ radiance by profile

      do headerIndex = 1, tvs_headerEnd
        ! deallocate BT equivalent to total direct, tl and ad radiance output
        deallocate(tvs_radiance(headerIndex) % bt)
      end do
      deallocate(tvs_radiance)
      do headerIndex = 1, tvs_headerEnd
        sensorIndex = tvs_lsensor(headerIndex)
        nlv_T = tvs_coefs(sensorIndex) % coef % nlevels
        ! deallocate model profiles atmospheric arrays with RTTOV levels dimension
        call rttov_alloc_prof(allocStatus,1,tvs_profiles_nl(headerIndex),nlv_T, &    ! 1 = nprofiles un profil a la fois
            tvs_opts(sensorIndex), asw=0,coefs=tvs_coefs(sensorIndex), init=.false. ) ! asw =0 deallocation
        if (allocStatus /= 0) call utl_abort('tvs_cleanup: memory deallocation error for tvs_profiles_nl')
        call rttov_alloc_prof(allocStatus,1,tvs_profiles_tlad(headerIndex),nlv_T, &    ! 1 = nprofiles un profil a la fois
             tvs_opts(sensorIndex),asw=0,coefs=tvs_coefs(sensorIndex),init=.false. ) ! asw =0 deallocation
        if (allocStatus /= 0) call utl_abort('tvs_cleanup: memory deallocation error for tvs_profiles_tlad')
      end do

      deallocate(tvs_profiles_nl)
      deallocate(tvs_profiles_tlad)

      do sensorIndex = tvs_nsensors, 1, -1
        call rttov_dealloc_coefs(allocStatus, tvs_coefs(sensorIndex))
        write(*,*) 'Deallocating coefficient structure for instrument', sensorIndex
        if (allocStatus /= 0) call utl_abort('tvs_cleanup: memory deallocation error in rttov_dealloc_coefs')
      end do

      deallocate(tvs_coefs)
      deallocate(tvs_listSensors)
      deallocate(tvs_opts)
    end if ! if (radiativeTransferCode == 'RTTOV')

    deallocate (tvs_nchan)
    deallocate (tvs_ichan)
    deallocate (tvs_lsensor)
    deallocate (tvs_isReallyPresent)
    deallocate (tvs_nchanMpiGlobal)
    deallocate (tvs_chanProf)
    deallocate (tvs_chanProfScatt)
    deallocate (tvs_bodyIndexFromBtIndex)
    deallocate (tvs_bodyIndexFromBtIndexScatt)

    write(*,*) 'tvs_cleanup: Finished'

  end subroutine tvs_cleanup

  !--------------------------------------------------------------------------
  ! tvs_deallocateProfilesNlTlAd
  !--------------------------------------------------------------------------
  subroutine tvs_deallocateProfilesNlTlAd
    !
    ! :Purpose: release memory used by RTTOV.
    !
    implicit none

    write(*,*) 'tvs_deallocateProfilesNlTlAd: Starting'

    if (radiativeTransferCode == 'RTTOV') then
      if (allocated(tvs_profiles_nl))   deallocate(tvs_profiles_nl)
      if (allocated(tvs_profiles_tlad)) deallocate(tvs_profiles_tlad)
      if (allocated(tvs_cld_profiles_nl)) deallocate(tvs_cld_profiles_nl)
      if (allocated(tvs_cld_profiles_tlad)) deallocate(tvs_cld_profiles_tlad)
    end if

    write(*,*) 'tvs_deallocateProfilesNlTlAd: Finished'

  end subroutine tvs_deallocateProfilesNlTlAd

  !--------------------------------------------------------------------------
  ! sensors
  !--------------------------------------------------------------------------
  subroutine sensors
    !
    !:Purpose: Initialisation of the RTTOV platform, satellite
    !          and instrument ID's. Also set burp to RTTOV channel
    !          mapping offset.
    !          To verify and transfom the sensor information contained in the
    !          NAMTOV namelist into the variables required by RTTTOV:
    !          platform, satellite and instrument ID's.
    !
    implicit none

    ! Locals:
    integer :: sensorIndex, instrumentIndex, platformIndex
    integer :: ipos1, ipos2
    integer :: numerosat, ierr, kindex
    character(len=15) :: tempocsatid
    logical, save :: firstCall=.true.
    integer, save :: ioffset1b(0:ninst-1)
    character(len=15) :: tempo_inst

    ! Namelist variables:
    character(len=8) :: listinstrum(0:ninst-1) ! List of instrument names
    integer          :: listoffset(0:ninst-1)  ! Corresponding list of channel offset values
    namelist /NAMCHANOFFSET/ listoffset, listinstrum

    !  1.0 Go through sensors and set RTTOV variables

    do sensorIndex=1, tvs_nsensors
      tvs_platforms  (sensorIndex) = -1
      tvs_satellites (sensorIndex) = -1
      tvs_instruments(sensorIndex) = -1
      tvs_channelOffset(sensorIndex) = -1
    end do

    if (firstCall) then
      if (utl_isNamelistPresent('NAMCHANOFFSET', './flnml')) then
        call utl_abort('sensors: NAMCHANOFFSET namelist section should be now in flnml_static !')
      end if
      ! read the namelist
      call utl_tmg_start(181,'low-level--readNML')
      listoffset(:) = 0
      listinstrum(:) = 'XXXXXXXX'
      read(utl_flnml_static,nml=NAMCHANOFFSET, iostat=ierr)
      if (ierr /= 0) then
        write(*,*) 'Error while reading NAMCHANOFFSET namelist section in flnml_static file !'
        call utl_abort('sensors')
      end if
      do instrumentIndex=0, ninst - 1
        if (listinstrum(instrumentIndex) /= 'XXXXXXXX') then
          ioffset1b( tvs_getInstrumentId( listinstrum(instrumentIndex) ) )  = listoffset(instrumentIndex)
        end if
      end do
      call utl_tmg_stop(181)
      firstCall = .false.
    end if

    !  1.1 Set platforms and satellites

    ! N.B.: Special cases for satellites TERRA and AQUA.
    !       For consistency with the RTTOV nomenclature, rename:
    !       TERRA  to  eos1
    !       AQUA   to  eos2
    !       NPP    to  jpss0
    !       JPSS    to  jpss0
    !       HMWARI    to  himawari
    !       FY-3C    to  FY3-3
    do sensorIndex = 1, tvs_nsensors
      if      (tvs_satelliteName(sensorIndex) == 'TERRA') then
        tempocsatid = 'eos1'
      else if (tvs_satelliteName(sensorIndex) == 'AQUA') then
        tempocsatid = 'eos2'
      else if (tvs_satelliteName(sensorIndex) == 'NPP') then
        tempocsatid = 'jpss0'
      else if (tvs_satelliteName(sensorIndex) == 'JPSS') then
        tempocsatid = 'jpss0'
      else if (tvs_satelliteName(sensorIndex)(1:6) == 'HMWARI') then
        tempocsatid = 'himawari' // trim(tvs_satelliteName(sensorIndex) (7:15))
      else if (tvs_satelliteName(sensorIndex) == 'FY-3C') then
        TEMPOCSATID = 'FY3-3'
      else
        call up2low(tvs_satelliteName(sensorIndex),tempocsatid)
      end if
      do platformIndex = 1, nplatforms
        ipos1 = len_trim(platform_name(platformIndex))
        ipos2 = index(tempocsatid,platform_name(platformIndex)(1:ipos1))
        if (ipos2 == 1) then
          tvs_platforms(sensorIndex) = platformIndex
          kindex = platformIndex
        end if
      end do
      if (tvs_platforms(sensorIndex) < 0) then
        write(*,'(A)') ' Satellite ' // trim(tempocsatid) // ' not supported.'
        call utl_abort('SENSORS')
      else
        ipos1 = len_trim(platform_name(kindex))
        ipos2 = len_trim(tempocsatid)
        read(tempocsatid(ipos1+1:ipos2),*,IOSTAT=ierr) numerosat
        numerosat = abs ( numerosat )
        if (ierr /= 0) then
          write(*,'(A,1x,i6,1x,i3,1x,i3,1x,A15)') 'Problem while reading satellite number', &
               ierr, ipos1, ipos2, tempocsatid
          call utl_abort('SENSORS')
        else
          tvs_satellites(sensorIndex) = numerosat
        end if
      end if
    end do

    !   1.2 Set instruments,
    !     also set channel offset, which is in fact a channel mapping between
    !     the channel number in BURP files and the channel number used in
    !     RTTOV.

    do sensorIndex = 1, tvs_nsensors
      if (tvs_instrumentName(sensorIndex)(1:10) == 'GOESIMAGER') then !cas particulier
        tvs_instruments(sensorIndex) = inst_id_goesim
      else if (tvs_satelliteName(sensorIndex)(1:5) == 'MTSAT') then   ! autre cas particulier
        tvs_instruments(sensorIndex) = inst_id_gmsim
      else                                                            ! cas general                 
        call up2low(tvs_instrumentName(sensorIndex),tempo_inst)
        do instrumentIndex = 0, ninst -1 
          if (trim(tempo_inst) == trim(inst_name(instrumentIndex))) then
            tvs_instruments(sensorIndex) = instrumentIndex
          end if
        end do
      end if
      if (tvs_instruments(sensorIndex) < 0) then
        write(*,'(A)') ' INSTRUMENT '// trim( tvs_instrumentName(sensorIndex)) // ' not supported.'
        call utl_abort('SENSORS')
      end if
      tvs_channelOffset(sensorIndex) = ioffset1b(tvs_instruments(sensorIndex))
    end do

    !    1.3 Print the RTTOV related variables

    if (mmpi_myid == 0) then
      write(*,*)
      write(*,'(3X,A)') '- SENSORS. Variables prepared for RTTOV:'
      write(*,'(3X,A)') '  ----------------------------------------'
      write(*,*)
      write(*,'(6X,A,I3)')   'Number of sensors       : ', tvs_nsensors
      write(*,"('Platform numbers        : ',6X,10I3)")  (tvs_platforms(sensorIndex), sensorIndex=1,tvs_nsensors)
      write(*,"('Satellite numbers       : ',6X,10I3)")  (tvs_satellites(sensorIndex), sensorIndex=1,tvs_nsensors)
      write(*,"('Instrument numbers      : ',6X,10I3)")  (tvs_instruments(sensorIndex), sensorIndex=1,tvs_nsensors)
      write(*,"('Channel mapping offsets : ',6X,10I3)")  (tvs_channelOffset(sensorIndex), sensorIndex=1,tvs_nsensors)
    end if

  end subroutine sensors

  !--------------------------------------------------------------------------
  !  tvs_getAllIdBurpTovs
  !--------------------------------------------------------------------------
  subroutine tvs_getAllIdBurpTovs(idatypListSize, idatypList)
    !
    ! :Purpose: Function to return a list of all idatyp (a.k.a. codtyp) values
    !           for all possible radiance observations (according to the namelist)
    !
    implicit none

    ! Argument:
    integer, intent(out) :: idatypListSize ! number of BUFR codtyp in the list
    integer, intent(out) :: idatypList(:)  ! list of radiance BUFR codtyp
    
    ! Locals:
    logical, save :: firstCall=.true.
    integer, save :: ninst_tovs
    integer :: ierr, instrumentIndex 
    integer, save :: list_inst(ninst)

    ! Namelist variables:
    character(len=22) :: inst_names(ninst) ! List of instrument names for all radiance types
    namelist /namtovsinst/ inst_names

    idatypList(:) = MPC_missingValue_int
    if (tvs_nsensors == 0) then
      ! no tovs data will be read, therefore false
      idatypListSize = 0      
      return
    end if

    if (firstCall) then
      if (utl_isNamelistPresent('NAMTOVSINST', './flnml')) then
        call utl_abort('tvs_getAllIdBurpTovs: NAMTOVSINST namelist section should be now in flnml_static !')
      end if
      ninst_tovs = 0
      list_inst(:) = -1
      inst_names(:) = 'XXXXXX'
      call utl_tmg_start(181,'low-level--readNML')
      read(utl_flnml_static, nml=namtovsinst, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_getAllIdBurpTovs: Error reading NAMTOVSINST namelist section in flnml_static file')
      if (mmpi_myid == 0) write(*,nml=namtovsinst)
      call utl_tmg_stop(181)
      do instrumentIndex=1, ninst
        if (inst_names(instrumentIndex) == 'XXXXXX') then
          ninst_tovs = instrumentIndex - 1
          exit
        else
          list_inst(instrumentIndex) = codtyp_get_codtyp( inst_names(instrumentIndex) )
          if (list_inst(instrumentIndex) < 0) then
            write(*,*) inst_names(instrumentIndex)
            call utl_abort('tvs_isIdBurpTovs: unknown instrument in namtovsinst namelist')
          end if
        end if
      end do
      if (ninst_tovs == 0) call utl_abort('tvs_getAllIdBurpTovs: Empty namtovsinst namelist')
      firstCall = .false.
    end if

    idatypList(1:ninst_tovs) = list_inst(1:ninst_tovs)
    idatypListSize = ninst_tovs

  end subroutine tvs_getAllIdBurpTovs

  !--------------------------------------------------------------------------
  !  tvs_countProfiles
  !--------------------------------------------------------------------------
  function tvs_countProfiles(sensorIndex) result(profileCount)
    !
    ! :Purpose: Function to count the number of radiances for a given sensor
    !
    implicit none
    ! Arguments:
    integer, intent(in) :: sensorIndex  ! RTTOV sensor index
    ! Result:
    integer             :: profileCount
    ! Local:
    integer :: headerIndex
    
    profileCount = 0
    do headerIndex = 1, tvs_headerEnd
      if (tvs_lsensor(headerIndex) == sensorIndex) then
        profileCount = profileCount + 1
      end if
    end do
    
  end function tvs_countProfiles

  !--------------------------------------------------------------------------
  !  tvs_isIdBurpTovs
  !--------------------------------------------------------------------------
  logical function tvs_isIdBurpTovs(idatyp)
    !
    ! :Purpose: Function to check if the given idatyp (a.k.a. codtyp) 
    !           corresponds to a radiance
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: idatyp  ! BUFR codtyp
    
    ! Locals:
    logical, save :: firstCall=.true.
    integer, save :: ninst_tovs
    integer :: ierr, instrumentIndex 
    integer, save :: list_inst(ninst)

    ! Namelist variables:
    character(len=22) :: inst_names(ninst) ! List of instrument names for all radiance types
    namelist /namtovsinst/ inst_names

    if (tvs_nsensors == 0) then
      ! no tovs data will be read, therefore false
      tvs_isIdBurpTovs = .false.
      return
    end if

    if (firstCall) then
      if (utl_isNamelistPresent('NAMTOVSINST', './flnml')) then
        call utl_abort('tvs_isIdBurpTovs: NAMTOVSINST namelist section should be now in flnml_static !')
      end if
      call utl_tmg_start(181,'low-level--readNML')
      ninst_tovs = 0
      list_inst(:) = -1
      inst_names(:) = 'XXXXXX'
      read(utl_flnml_static, nml=namtovsinst, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isIdBurpTovs: Error reading NAMTOVSINST namelist section in flnml_static file')
      if (mmpi_myid == 0) write(*,nml=namtovsinst)
      call utl_tmg_stop(181)
      do instrumentIndex=1, ninst
        if (inst_names(instrumentIndex) == 'XXXXXX') then
          ninst_tovs= instrumentIndex - 1
          exit
        else
          list_inst(instrumentIndex) = codtyp_get_codtyp( inst_names(instrumentIndex) )
          if (list_inst(instrumentIndex) < 0) then
            write(*,*) inst_names(instrumentIndex)
            call utl_abort('tvs_isIdBurpTovs: unknown instrument in namtovsinst namelist')
          end if
        end if
      end do
      if ( ninst_tovs == 0 ) call utl_abort('tvs_isIdBurpTovs: Empty namtovsinst namelist')
      firstCall = .false.
    end if
    
    tvs_isIdBurpTovs = .false.

    do instrumentIndex = 1, ninst_tovs
      if (idatyp == list_inst(instrumentIndex)) then
        tvs_isIdBurpTovs = .true.
        exit
      end if
    end do

  end function tvs_isIdBurpTovs

  !--------------------------------------------------------------------------
  !  tvs_isIdBurpHyperSpectral
  !--------------------------------------------------------------------------
  logical function tvs_isIdBurpHyperSpectral(idatyp)
    !
    ! :Purpose: Function to check if the given idatyp (a.k.a. codtyp) 
    !           corresponds to a hyper-spectral infrared radiance
    !
    implicit none

    ! Argument:
    integer, intent(in) :: idatyp ! BUFR codtyp
    
    ! Locals:
    logical, save :: firstCall=.true.
    integer, save :: ninst_hyper
    integer :: ierr, instrumentIndex 
    integer, save :: list_inst(ninst)

    ! Namelist variables:
    character(len=22) :: name_inst(ninst) ! List of instrument names for hyperspectral IR
    namelist /namhyper/ name_inst

    if (tvs_nsensors == 0) then
      ! no tovs data will be read, therefore false
      tvs_isIdBurpHyperSpectral = .false.
      return
    end if

    if (firstCall) then
      if (utl_isNamelistPresent('NAMHYPER', './flnml')) then
        call utl_abort('tvs_isIdBurpHyperSpectral: NAMHYPER namelist section should be now in flnml_static !')
      end if
      ninst_hyper = 0
      list_inst(:) = -1
      name_inst(:) = 'XXXXXX'
      call utl_tmg_start(181,'low-level--readNML')
      read(utl_flnml_static, nml=namhyper, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isIdBurpHyperSpectral: Error reading NAMHYPER namelist section in flnml_static file')
      if (mmpi_myid == 0) write(*,nml=namhyper)
      call utl_tmg_stop(181)
      do instrumentIndex=1, ninst
        if (name_inst(instrumentIndex) == 'XXXXXX' ) then
          ninst_hyper = instrumentIndex - 1
          exit
        else
          list_inst(instrumentIndex) = codtyp_get_codtyp( name_inst(instrumentIndex) )
          if (list_inst(instrumentIndex) < 0) then
            write(*,*) name_inst(instrumentIndex)
            call utl_abort('tvs_isIdBurpHyperSpectral: unknown instrument in namhyper namelist')
          end if
        end if
      end do
      if ( ninst_hyper == 0 ) call utl_abort('tvs_isIdBurpHyperSpectral: Empty namhyper namelist')
      firstCall = .false.
    end if

    tvs_isIdBurpHyperSpectral = .false.

    do instrumentIndex = 1, ninst_hyper
      if (idatyp == list_inst(instrumentIndex)) then
        tvs_isIdBurpHyperSpectral = .true.
        exit
      end if
    end do

  end function tvs_isIdBurpHyperSpectral

  !--------------------------------------------------------------------------
  !  tvs_isIdBurpInst
  !--------------------------------------------------------------------------
  logical function tvs_isIdBurpInst(idburp,cinst)
    !
    ! :Purpose: function to check if the provided idburp (a.k.a. codtyp) corresponds to instrument cinst
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: idburp ! Input BURP codtyp
    character(len=*), intent(in) :: cinst  ! Input instrument name

    if (tvs_nsensors == 0) then
      ! no tovs data will be read, therefore false
      tvs_isIdBurpInst = .false.
      return
    end if

    tvs_isIdBurpInst = ( idburp == codtyp_get_codtyp(cinst) )

  end function tvs_isIdBurpInst

  !--------------------------------------------------------------------------
  !  tvs_getPlatformId
  !--------------------------------------------------------------------------
  integer function tvs_getPlatformId(name)
    !
    ! :Purpose: return RTTOV platform id (>0) from platform name.
    !           -1 if not found
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: name ! Platform name

    ! Locals:
    integer           :: platformIndex, length, ipos
    character(len=64) :: tempo_name

    tvs_getPlatformId = -1
    length = len_trim(name)
    call up2low(name(1:length),tempo_name(1:length))

    if (index(tempo_name(1:length),'npp') /= 0) then
      tvs_getPlatformId = platform_id_jpss
    else if (index(tempo_name(1:length),'hmwari') /= 0) then
      tvs_getPlatformId = platform_id_himawari
    else
      do platformIndex = 1, nplatforms
        ipos = index(tempo_name(1:length),trim(platform_name(platformIndex)))
        if (ipos == 1) then
          tvs_getPlatformId = platformIndex
          exit
        end if
      end do
    end if

  end function tvs_getPlatformId

  !--------------------------------------------------------------------------
  !  tvs_getInstrumentId
  !--------------------------------------------------------------------------
  integer function tvs_getInstrumentId(name)
    !
    ! :Purpose: return RTTOV instrument id from intrument name. 0 is a valid answer.
    !           -1 if not found
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: name ! Instrument name

    ! Locals:
    integer           :: instrumentIndex, length
    character(len=64) :: tempo_name

    tvs_getInstrumentId = -1
    length = len_trim(name)
    call up2low(name(1:length),tempo_name(1:length))
    if (trim(tempo_name(1:length)) == 'goesim') then
      tvs_getInstrumentId = inst_id_goesim
    else if (trim(tempo_name(1:length)) == 'gmsim') then
      tvs_getInstrumentId = inst_id_gmsim
    else if (trim(tempo_name(1:length)) == 'mtsatim') then
      tvs_getInstrumentId = inst_id_mtsatim
    else
      do instrumentIndex = 0, ninst - 1
        if (trim(inst_name(instrumentIndex)) == tempo_name(1:length)) then
          tvs_getInstrumentId = instrumentIndex
          exit
        end if
      end do
    end if
    
  end function tvs_getInstrumentId

  !--------------------------------------------------------------------------
  !  tvs_isInstrumHyperSpectral
  !--------------------------------------------------------------------------
  logical function tvs_isInstrumHyperSpectral(instrum)
    !
    ! :Purpose: given an RTTOV instrument code return if it is an hyperspectral one
    !           information from namelist NAMHYPER
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: instrum     ! input RTTOV instrument code

    ! Locals:
    integer :: ierr, instrumentIndex 
    integer, save :: list_inst(maxsize), ninst_hir
    logical, save :: firstCall = .true.

    ! Namelist variables:
    character (len=8) :: name_inst(maxsize) ! List of instrument names for hyperspectral IR
    namelist /NAMHYPER/ name_inst
    
    if (firstCall) then
      if (utl_isNamelistPresent('NAMHYPER', './flnml')) then
        call utl_abort('tvs_isInstrumHyperSpectral: NAMHYPER namelist section should be now in flnml_static !')
      end if
      call utl_tmg_start(181,'low-level--readNML')
      ninst_hir = 0
      name_inst(:) = 'XXXXXXX'
      read(utl_flnml_static, nml=namhyper, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isInstrumHyperSpectral: Error reading namelist section NAMHYPER in flnm_static file')
      if (mmpi_myid == 0) write(*,nml=namhyper)
      call utl_tmg_stop(181)
      list_inst(:) = -1
      do instrumentIndex=1, maxsize
        list_inst(instrumentIndex) = tvs_getInstrumentId( name_inst(instrumentIndex) )
        if (name_inst(instrumentIndex) /= 'XXXXXXX') then
          if (list_inst(instrumentIndex) == -1) then
            write(*,*) instrumentIndex,name_inst(instrumentIndex)
            call utl_abort('tvs_isInstrumHyperSpectral: Unknown instrument name')
          end if
        else
          ninst_hir = instrumentIndex - 1
          exit
        end if
      end do
      firstCall = .false.
      if (ninst_hir == 0) then
        write(*,*) 'tvs_isInstrumHyperSpectral: Warning : empty namhyper namelist !'
      end if
    end if
    tvs_isInstrumHyperSpectral = .false.
    do instrumentIndex = 1, ninst_hir
      if (instrum == list_inst(instrumentIndex)) then
        tvs_isInstrumHyperSpectral = .true.
        exit
      end if
    end do

  end function tvs_isInstrumHyperSpectral

  !--------------------------------------------------------------------------
  !  tvs_isNameHyperSpectral
  !--------------------------------------------------------------------------
  logical function tvs_isNameHyperSpectral(instrumentName)
    !
    ! :Purpose: given an instrument name
    !           returns if it is an hyperspectral one
    !           (information from namelist NAMHYPER)
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: instrumentName ! Instrument name

    !Locals:
    integer :: ierr, instrumentIndex 
    integer, save :: ninst_hir
    logical, save :: firstCall = .true.
    character (len=8) :: name2

    ! Namelist variables:
    character (len=8),save  :: name_inst(maxsize) ! List of instrument names for hyperspectral IR
    namelist /NAMHYPER/ name_inst

    if (firstCall) then
      ninst_hir = 0
      name_inst(:) = 'XXXXXXX'
      call utl_tmg_start(181,'low-level--readNML')
      read(utl_flnml_static, nml=namhyper, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isNameHyperSpectral: Error reading NAMHYPER namelist section in flnml_static file')
      if (mmpi_myid == 0) write(*,nml=namhyper)
      call utl_tmg_stop(181)
      do instrumentIndex = 1, maxsize
        if (name_inst(instrumentIndex) == 'XXXXXXX') then
          ninst_hir = instrumentIndex -1
          exit
        end if
      end do
      firstCall = .false.
      if (ninst_hir == 0) then
        write(*,*) 'tvs_isNameHyperSpectral: Warning : empty namhyper namelist !'
      end if
    end if

    tvs_isNameHyperSpectral = .false.

    call up2low(instrumentName, name2)

    do instrumentIndex = 1, ninst_hir
      if (trim(name2) == name_inst(instrumentIndex)) then
        tvs_isNameHyperSpectral = .true.
        exit
      end if
    end do

  end function tvs_isNameHyperSpectral

  !--------------------------------------------------------------------------
  !  tvs_isInstrumGeostationary
  !--------------------------------------------------------------------------
  logical function tvs_isInstrumGeostationary(instrum)
    !
    ! :Purpose: given an RTTOV instrument code return if it is a Geostationnary Imager
    !           information from namelist NAMGEO
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: instrum ! input RTTOV instrument code

    ! Locals:
    integer :: ierr, instrumentIndex 
    integer, save :: list_inst(maxsize), ninst_geo
    logical, save :: firstCall = .true.

    ! Namelist variables:
    character(len=8) :: name_inst(maxsize) ! List of instrument names for geostationary
    namelist /NAMGEO/ name_inst

    if (firstCall) then
      if (utl_isNamelistPresent('NAMGEO', './flnml')) then
        call utl_abort('tvs_isInstrumGeostationary: NAMGEO namelist section should be now in flnml_static !')
      end if
      ninst_geo = 0
      name_inst(:) = 'XXXXXX'
      call utl_tmg_start(181,'low-level--readNML')
      read(utl_flnml_static, nml=namgeo, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isInstrumGeostationary: Error reading namelist section NAMGEO in flnml_static file')
      if (mmpi_myid == 0) write(*,nml=namgeo)
      call utl_tmg_stop(181)
      list_inst(:) = -1
      do instrumentIndex=1, maxsize
        list_inst(instrumentIndex) = tvs_getInstrumentId( name_inst(instrumentIndex) )
        if (name_inst(instrumentIndex) /= 'XXXXXX') then
          if (list_inst(instrumentIndex) == -1) then
            write(*,*) instrumentIndex,name_inst(instrumentIndex)
            call utl_abort('tvs_isInstrumGeostationary: Unknown instrument name')
          end if
        else
          ninst_geo = instrumentIndex - 1
          exit
        end if
      end do
      firstCall = .false.
      if (ninst_geo == 0) then
        write(*,*) 'tvs_isInstrumGeostationary: Warning : empty namgeo namelist !'
      end if
    end if
    tvs_isInstrumGeostationary = .false.
    do instrumentIndex = 1, ninst_geo
      if (instrum == list_inst(instrumentIndex)) then
        tvs_isInstrumGeostationary = .true.
        exit
      end if
    end do
    
  end function tvs_isInstrumGeostationary

  !--------------------------------------------------------------------------
  !  tvs_isInstrumUsingCLW
  !--------------------------------------------------------------------------
  function tvs_isInstrumUsingCLW(instrumId) result(idExist)
    !
    ! :Purpose: given an RTTOV instrument code return if it is in the list to use CLW
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: instrumId     ! input RTTOV instrument code
    ! Result:
    logical             :: idExist

    ! Locals:
    integer :: instrumentIndex 

    idExist = .false.
    do instrumentIndex = 1, tvs_numMWInstrumUsingCLW
      if (instrumId == instrumentIdsUsingCLW(instrumentIndex)) then
        idExist = .true.
        exit
      end if
    end do

  end function tvs_isInstrumUsingCLW

  !--------------------------------------------------------------------------
  !  tvs_getClwIndex
  !--------------------------------------------------------------------------
  function tvs_getClwIndex(instrumId) result(clwIndex)
    !
    ! :Purpose: given an RTTOV instrument code return if it is in the list to use CLW
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: instrumId     ! input RTTOV instrument code
    ! Result:
    integer             :: clwIndex

    ! Locals:
    integer :: instrumentIndex 

    clwIndex = -1
    do instrumentIndex = 1, tvs_numMWInstrumUsingCLW
      if (instrumId == instrumentIdsUsingCLW(instrumentIndex)) then
        clwIndex = instrumentIndex
        exit
      end if
    end do

  end function tvs_getClwIndex

  !--------------------------------------------------------------------------
  !  tvs_isInstrumUsingHydrometeors
  !--------------------------------------------------------------------------
  function tvs_isInstrumUsingHydrometeors(instrumId) result(idExist)
    !
    ! :Purpose: given an RTTOV instrument code return if it is in the list to use Hydrometeors
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: instrumId     ! input RTTOV instrument code
    ! Result:
    logical             :: idExist

    idExist = tvs_getHydrometeorsIndex(instrumId) > 0
    
  end function tvs_isInstrumUsingHydrometeors

  !--------------------------------------------------------------------------
  !  tvs_getHydrometeorsIndex
  !--------------------------------------------------------------------------
  function tvs_getHydrometeorsIndex(instrumId) result(hydrometeorsIndex)
    !
    ! :Purpose: given an RTTOV instrument code return if it is in the list to use Hydrometeors
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: instrumId     ! input RTTOV instrument code
    ! Result:
    integer             :: hydrometeorsIndex

    ! Locals:
    integer :: instrumentIndex 

    hydrometeorsIndex = -1
    do instrumentIndex = 1, tvs_numMWInstrumUsingHydrometeors
      if (instrumId == instrumentIdsUsingHydrometeors(instrumentIndex)) then
        hydrometeorsIndex = instrumentIndex
        exit
      end if
    end do

  end function tvs_getHydrometeorsIndex
  
  !--------------------------------------------------------------------------
  !  tvs_checkAllskyChanNum
  !--------------------------------------------------------------------------
  subroutine tvs_checkAllskyChanNum(useStateDepSigmaObs)
    !
    ! :Purpose: Check all-sky channel numbers in the array filled with symmetricObsErr ascii 
    !           file match the all-sky channel numbers from NAMTOV.
    !
    implicit none

    ! Argument:
    logical, intent(in) :: useStateDepSigmaObs(:,:) ! array filled with symmetricObsErr ascii
    
    ! Locals:
    integer :: sensorIndex, sensorIndex2, instrumId
    integer :: chanNum, chanNumWithOffset
    integer :: chanNumWithOffsetStart, chanNumWithOffsetEnd
    integer :: numChanInAsciiList, numChanInNamtovList
    logical :: isChannelInNamtovList, isChannelInAsciiList

    sensorLoop: do sensorIndex = 1, tvs_nsensors
      if (.not. tvs_mwAllskyAssim .or. .not. any(useStateDepSigmaObs(:,sensorIndex))) cycle sensorLoop

      instrumId = tvs_instruments(sensorIndex)
    
      ! all-sky HU
      if (tvs_isInstrumUsingHydrometeors(instrumId)) then
        sensorIndex2 = tvs_getHydrometeorsIndex(instrumId)

        ! check tvs_channelsUsingHydrometeors and useStateDepSigmaObs have same length for this instrument
        if (instrumId == tvs_getInstrumentId('atms')) then
          chanNumWithOffsetStart = atmsHuChanNum(1)
          chanNumWithOffsetEnd = atmsHuChanNum(size(atmsHuChanNum(:)))
        else
          chanNumWithOffsetStart = 1
          chanNumWithOffsetEnd = tvs_maxChannelNumber
        end if
        numChanInAsciiList = count(useStateDepSigmaObs(chanNumWithOffsetStart:chanNumWithOffsetEnd,sensorIndex))
        numChanInNamtovList = count(tvs_channelsUsingHydrometeors(sensorIndex2,:)>0)
        if (numChanInAsciiList /= numChanInNamtovList) then
          write(*,*) 'tvs_checkAllskyChanNum: sensorIndex=', sensorIndex, &
                    ', numChanInAsciiList=', numChanInAsciiList, &
                    ', numChanInNamtovList=', numChanInNamtovList, &
                    ', chanNumWithOffsetStart=', chanNumWithOffsetStart, &
                    ', chanNumWithOffsetEnd=', chanNumWithOffsetEnd, &
                    ', useStateDepSigmaObs=', useStateDepSigmaObs(chanNumWithOffsetStart:chanNumWithOffsetEnd,sensorIndex)
                    
          call utl_abort('tvs_checkAllskyChanNum: numChanInAsciiList /= numChan in tvs_channelsUsingHydrometeors')
        end if

        ! check channel list from NAMTOV match channel list from symmetricObsErr ascii file
        do chanNumWithOffset = chanNumWithOffsetStart, chanNumWithOffsetEnd
          chanNum = chanNumWithOffset - tvs_channelOffset(sensorIndex)

          isChannelInNamtovList = (utl_findloc(tvs_channelsUsingHydrometeors(sensorIndex2,:),chanNum) > 0)
          isChannelInAsciiList = useStateDepSigmaObs(chanNumWithOffset,sensorIndex)

          if (chanNum == -1 .and. .not. isChannelInAsciiList) cycle

          if ((      isChannelInNamtovList .and. .not. isChannelInAsciiList) .or. &
              (.not. isChannelInNamtovList .and.       isChannelInAsciiList)) then

            write(*,*) 'tvs_checkAllskyChanNum: sensorIndex=', sensorIndex, &
                      ', sensorIndex2=', sensorIndex2, &
                      ', inst=', instrumId, &
                      ', chanNum (or tvs_channelsUsingHydrometeors)=', chanNum, &
                      ', isChannelInNamtovList=', isChannelInNamtovList, &
                      ', chanNumWithOffset=', chanNumWithOffset, &
                      ', isChannelInAsciiList=', isChannelInAsciiList, &
                      ', useStateDepSigmaObs=', useStateDepSigmaObs(chanNumWithOffset,sensorIndex)

            call utl_abort('tvs_checkAllskyChanNum: useStateDepSigmaObs and tvs_channelsUsingHydrometeors not matching')
          end if
        end do ! do chanNum = 1, tvs_maxNumberOfChannels
      end if ! if (tvs_isInstrumUsingHydrometeors(instrumId)) then

      ! all-sky TT
      if (tvs_isInstrumUsingCLW(instrumId)) then
        sensorIndex2 = tvs_getClwIndex(instrumId)

        ! check tvs_channelsUsingClw and useStateDepSigmaObs have same length for this instrument
        if (instrumId == tvs_getInstrumentId('atms')) then
          chanNumWithOffsetStart = atmsTtChanNum(1)
          chanNumWithOffsetEnd = atmsTtChanNum(size(atmsTtChanNum(:)))
        else
          chanNumWithOffsetStart = 1
          chanNumWithOffsetEnd = tvs_maxChannelNumber
        end if
        numChanInAsciiList = count(useStateDepSigmaObs(chanNumWithOffsetStart:chanNumWithOffsetEnd,sensorIndex))
        numChanInNamtovList = count(tvs_channelsUsingClw(sensorIndex2,:)>0)
        if (numChanInAsciiList /= numChanInNamtovList) then
          write(*,*) 'tvs_checkAllskyChanNum: sensorIndex=', sensorIndex, &
                    ', numChanInAsciiList=', numChanInAsciiList, &
                    ', numChanInNamtovList=', numChanInNamtovList, &            
                    ', chanNumWithOffsetStart=', chanNumWithOffsetStart, &
                    ', chanNumWithOffsetEnd=', chanNumWithOffsetEnd, &
                    ', useStateDepSigmaObs=', useStateDepSigmaObs(chanNumWithOffsetStart:chanNumWithOffsetEnd,sensorIndex)
                    
          call utl_abort('tvs_checkAllskyChanNum: numChanInAsciiList /= numChan in tvs_channelsUsingClw')
        end if          

        ! check channel list from NAMTOV match channel list from symmetricObsErr ascii file
        do chanNumWithOffset = chanNumWithOffsetStart, chanNumWithOffsetEnd
          chanNum = chanNumWithOffset - tvs_channelOffset(sensorIndex)

          isChannelInNamtovList = (utl_findloc(tvs_channelsUsingClw(sensorIndex2,:),chanNum) > 0)
          isChannelInAsciiList = useStateDepSigmaObs(chanNumWithOffset,sensorIndex)

          if (chanNum == -1 .and. .not. isChannelInAsciiList) cycle

          if ((      isChannelInNamtovList .and. .not. isChannelInAsciiList) .or. &
              (.not. isChannelInNamtovList .and.       isChannelInAsciiList)) then

            write(*,*) 'tvs_checkAllskyChanNum: sensorIndex=', sensorIndex, &
                      ', sensorIndex2=', sensorIndex2, &
                      ', inst=', instrumId, &
                      ', chanNum (or tvs_channelsUsingClw)=', chanNum, &
                      ', isChannelInNamtovList=', isChannelInNamtovList, &
                      ', chanNumWithOffset=', chanNumWithOffset, &
                      ', isChannelInAsciiList=', isChannelInAsciiList, &
                      ', useStateDepSigmaObs=', useStateDepSigmaObs(chanNumWithOffset,sensorIndex)

            call utl_abort('tvs_checkAllskyChanNum: useStateDepSigmaObs and tvs_channelsUsingClw not matching')
          end if
        end do ! do chanNum = 1, tvs_maxNumberOfChannels
      end if ! if (tvs_isInstrumUsingCLW(instrumId)) then
    end do sensorLoop

  end subroutine tvs_checkAllskyChanNum

  !--------------------------------------------------------------------------
  !  tvs_isChanNumInAllskyNamtovList
  !--------------------------------------------------------------------------
  function tvs_isChanNumInAllskyNamtovList(instrumId,allskyTtHu,channelNumber) result(isChannelInNamtovList)
    !
    ! :Purpose: check if channel number is in NAMTOV list for all-sky TT/HU of the instrument.
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: instrumId     ! input RTTOV instrument code
    character(len=2), intent(in) :: allskyTtHu    ! 'TT' for all-sky temperature, HU for all-sky humidity
    integer,          intent(in) :: channelNumber ! channel number

    ! Result:
    logical :: isChannelInNamtovList

    ! Locals:
    integer :: sensorIndex2 
    
    isChannelInNamtovList = .false.

    ! all-sky HU
    if (allskyTtHu == 'HU') then
      if (tvs_isInstrumUsingHydrometeors(instrumId)) then
        sensorIndex2 = tvs_getHydrometeorsIndex(instrumId)
        isChannelInNamtovList = (utl_findloc(tvs_channelsUsingHydrometeors(sensorIndex2,:),channelNumber) > 0)
      end if
    end if

    ! all-sky TT
    if (allskyTtHu == 'TT') then
      if (tvs_isInstrumUsingCLW(instrumId)) then
        sensorIndex2 = tvs_getClwIndex(instrumId)
        isChannelInNamtovList = (utl_findloc(tvs_channelsUsingClw(sensorIndex2,:),channelNumber) > 0)
      end if
    end if

  end function tvs_isChanNumInAllskyNamtovList

  !--------------------------------------------------------------------------
  !  tvs_isInstrumAllskyTtAssim
  !--------------------------------------------------------------------------
  function tvs_isInstrumAllskyTtAssim(instrumId) result(allskyTtAssim)
    !
    ! :Purpose: determine if all-sky temperature-channel assimilation is asked for the instrument.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: instrumId     ! input RTTOV instrument code
    ! Result:
    logical             :: allskyTtAssim

    allskyTtAssim = (tvs_mwAllskyAssim .and. tvs_isInstrumUsingCLW(instrumId))

  end function tvs_isInstrumAllskyTtAssim

  !--------------------------------------------------------------------------
  !  tvs_isInstrumAllskyHuAssim
  !--------------------------------------------------------------------------
  function tvs_isInstrumAllskyHuAssim(instrumId) result(allskyHuAssim)
    !
    ! :Purpose: determine if all-sky humidity-channel assimilation is asked for the instrument.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: instrumId     ! input RTTOV instrument code
    ! Result:
    logical             :: allskyHuAssim

    allskyHuAssim = (tvs_mwAllskyAssim .and. tvs_isInstrumUsingHydrometeors(instrumId))

  end function tvs_isInstrumAllskyHuAssim

  !--------------------------------------------------------------------------
  !  tvs_mapInstrum
  !--------------------------------------------------------------------------
  subroutine tvs_mapInstrum(instrumburp,instrum)
    !
    ! :Purpose:  Map BUFR satellite instrument (element #2019) to RTTOV instrument.
    !            A negative value is returned, if no match in found.
    !
    ! :Table of  RTTOV instrument identifier:
    !
    ! ==================  =====================  ==================
    ! Instrument          Instrument identifier  Sensor type
    ! ==================  =====================  ==================
    !               HIRS               0                     ir
    !                MSU               1                     mw
    !                SSU               2                     ir
    !              AMSUA               3                     mw
    !              AMSUB               4                     mw
    !              AVHRR               5                     ir
    !               SSMI               6                     mw
    !              VTPR1               7                     ir
    !              VTPR2               8                     ir
    !                TMI               9                     mw
    !              SSMIS              10                     mw
    !               AIRS              11                     ir
    !              MODIS              13                     ir
    !               ATSR              14                     ir
    !                MHS              15                     mw
    !               ATMS              19                     mw
    !              MVIRI              20                     ir
    !             SEVIRI              21                     ir
    !         GOESIMAGER              22                     ir
    !        GOESSOUNDER              23                     ir
    !   GMS/MTSAT IMAGER              24                     ir
    !          FY2-VISSR              25                     ir
    !          FY1-MVISR              26                     ir
    !                AHI              56                     ir
    ! ==================  =====================  ==================
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: instrumburp  ! BUFR satellite instrument (element #2019)
    integer, intent(out) :: instrum      ! RTTOV instrument ID numbers (e.g. 3 for  AMSUA)
  
    ! Locals:  
    integer :: instrumentIndex, numinstburp
    integer, parameter :: mxinstrumburp = 100
    logical, save :: firstCall = .true.
    integer :: ierr

    ! Namelist variables:
    integer, save ::   listburp(mxinstrumburp)           ! List of instrument ID values from obs file
    character(len=8), save :: listinstrum(mxinstrumburp) ! List of instrument names
    namelist /NAMINST/ listburp, listinstrum

    !      Table of BURP satellite sensor identifier element #002019

    !   1.0 Find instrument

    if (firstCall) then
      if (utl_isNamelistPresent('NAMINST', './flnml')) then
        call utl_abort('tvs_mapInstrum: NAMINST namelist section should be now in flnml_static !')
      end if
      
      ! set the default values
      listburp(:) = -1
      listinstrum(:) = 'XXXXXXXX'

      ! read the namelist
      call utl_tmg_start(181,'low-level--readNML')
      read(utl_flnml_static, nml=NAMINST, iostat=ierr)
      if (ierr /= 0) then
        write(*,*) 'Error while reading NAMINST namelist section in flnml_static file !'
        call utl_abort('tvs_mapInstrum')
      end if
      call utl_tmg_stop(181)

      ! figure out how many valid elements in the lists
      do instrumentIndex=1, mxinstrumburp
        if (listburp(instrumentIndex) == -1) then
          numinstburp = instrumentIndex - 1
          exit
        end if
      end do
      if (numinstburp > mxinstrumburp) then
        call utl_abort('tvs_mapInstrum: exceeded maximum number of platforms')
      end if
      write(*,*) 'tvs_mapInstrum: number of satellites found in namelist = ',numinstburp
      write(*,*) 'tvs_mapInstrum: listburp   = ',listburp(1:numinstburp)
      write(*,*) 'tvs_mapInstrum: listinstrum    = ',listinstrum(1:numinstburp)
      firstCall = .false.
    end if

    instrum = -1
    do instrumentIndex = 1, mxinstrumburp
      if (instrumburp == listburp(instrumentIndex)) then
        instrum = tvs_getInstrumentId( listinstrum(instrumentIndex) )
        exit
      end if
    end do

  end subroutine tvs_mapInstrum

  !--------------------------------------------------------------------------
  !  tvs_isNameGeostationary
  !--------------------------------------------------------------------------
  logical function tvs_isNameGeostationary(instrumentName)
    !
    ! :Purpose: given an instrument name following BUFR convention
    !           returns if it is a Geostationnary Imager
    !           (information from namelist NAMGEOBUFR)
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: instrumentName ! Instrument name

    !Locals:
    integer :: ierr, instrumentIndex
    integer, save :: ninst_geo
    logical, save :: firstCall = .true.

    ! Namelist variables:
    character (len=8), save :: name_inst(maxsize) ! List of instrument names for geostationary
    namelist /NAMGEOBUFR/ name_inst

    if (firstCall) then
      if (utl_isNamelistPresent('NAMGEOBUFR', './flnml')) then
        call utl_abort('tvs_isNameGeostationary: NAMGEOBUFR namelist section should be now in flnml_static !')
      end if
      ninst_geo = 0
      name_inst(:) = 'XXXXXXXX'
      call utl_tmg_start(181,'low-level--readNML')
      read(utl_flnml_static, nml=namgeobufr, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isNameGeostationary: Error reading namelist section NAMGEOBUFR in flnml_static_file')
      if (mmpi_myid == 0) write(*,nml=namgeobufr)
      call utl_tmg_stop(181)
      do instrumentIndex = 1, maxsize
        if (name_inst(instrumentIndex) == 'XXXXXXXX') then
          ninst_geo = instrumentIndex - 1
          exit
        end if
      end do
      firstCall = .false.
      if (ninst_geo == 0) then
        write(*,*) 'tvs_isNameGeostationary: Warning : empty namgeobufr namelist !' 
      end if
    end if
    
    tvs_isNameGeostationary = .false.
    do instrumentIndex = 1, ninst_geo
      if (trim(instrumentName) == trim(name_inst(instrumentIndex))) then
        tvs_isNameGeostationary= .true.
        exit
      end if
    end do
    
  end function tvs_isNameGeostationary

  !--------------------------------------------------------------------------
  !  tvs_mapSat
  !--------------------------------------------------------------------------
  subroutine tvs_mapSat(isatBURP,iplatform,isat)
    !
    ! :Purpose:  Map BUFR satellite identifier (element #1007)
    !            to RTTOV platform and satellite.
    !            Negative values are returned, if no match in found.
    !
    ! :Table of  RTTOV platform identifier:
    !
    ! ========          =========================
    ! Platform          RTTOV platform identifier
    ! ========          =========================
    !     NOAA               1
    !     DMSP               2
    ! METEOSAT               3
    !     GOES               4
    !      GMS               5
    !      FY2               6
    !     TRMM               7
    !      ERS               8
    !      EOS               9
    !    METOP              10
    !  ENVISAT              11
    !      MSG              12
    !      FY1              13
    !    ADEOS              14
    !    MTSAT              15
    ! CORIOLIS              16
    !      NPP              17
    ! ========          ===========================
    !
    ! :Example: 
    !          NOAA15, which has a BUFR satellite identifier value of 206,
    !          is mapped into the following:
    !          RTTOV platform  =  1,
    !          RTTOV satellite = 15.
    !
    ! :Arguments:
    !     :isatBURP: BUFR satellite identifier
    !     :iplatform: RTTOV platform ID numbers (e.g. 1 for  NOAA)
    !     :isat: RTTOV satellite ID numbers (e.g. 15)
    !
    implicit none
    
    ! Arguments:
    integer, intent(in)  :: isatburp   ! BUFR satellite identifier
    integer, intent(out) :: iplatform  ! RTTOV platform ID numbers (e.g. 1 for  NOAA)
    integer, intent(out) :: isat       ! RTTOV satellite ID numbers (e.g. 15)

    ! Locals:
    integer           :: satelliteIndex, ierr
    logical, save     :: firstCall=.true.
    integer, parameter:: mxsatburp = 100
    integer, save     :: numsatburp

    ! Namelist variables:
    integer, save          :: listburp(mxsatburp) ! Table of BURP satellite identifier element #001007
    character(len=8), save :: listplat(mxsatburp) ! Table of RTTOV platform identifier
    integer, save          :: listsat (mxsatburp) ! Table of RTTOV satellite identifier

    namelist /NAMSAT/ listburp, listplat, listsat

    !     Fill tables from namelist at the first call 
    if (firstCall) then
      if (utl_isNamelistPresent('NAMSAT', './flnml')) then
        call utl_abort('tvs_mapSat: NAMSAT namelist section should be now in flnml_static !')
      end if
      ! set the default values
      listburp(:) = -1
      listsat(:) = -1
      listplat(:) = 'XXXXXXXX'
      ! read the namelist
      call utl_tmg_start(181,'low-level--readNML')
      read(utl_flnml_static, nml=NAMSAT, iostat = ierr)
      if (ierr /= 0) then
        write(*,*) 'Error while reading NAMSAT namelist section in flnml_static file !'
        call utl_abort('tvs_mapSat')
      end if
      call utl_tmg_stop(181)

      !  Figure out how many valid elements in the lists
      do satelliteIndex=1, mxsatburp
        if (listburp(satelliteIndex) == -1) then
          numsatburp = satelliteIndex - 1
          exit
        end if
      end do
      if (numsatburp >= mxsatburp) then
        call utl_abort('tvs_mapSat: exceeded maximum number of platforms')
      end if
      write(*,*) 'tvs_mapSat: number of satellites found in namelist = ',numsatburp
      write(*,*) 'tvs_mapSat: listburp   = ',listburp(1:numsatburp)
      write(*,*) 'tvs_mapSat: listsat    = ',listsat(1:numsatburp)
      write(*,*) 'tvs_mapSat: listplat   = ',listplat(1:numsatburp)
      firstCall = .false.
    end if

    !   Find platform and satellite

    iplatform = -1
    isat      = -1
    do satelliteIndex = 1, numsatburp
      if (isatburp == listburp(satelliteIndex)) then
        iplatform = tvs_getPlatformId( listplat(satelliteIndex) )
        isat = listsat (satelliteIndex)
        exit
      end if
    end do

  end subroutine tvs_mapSat

  !--------------------------------------------------------------------------
  !  tvs_getChanProf
  !--------------------------------------------------------------------------
  subroutine tvs_getChanprof(sensorHeaderIndexes, obsSpaceData, chanprof, lchannel_subset_opt, iptobs_cma_opt, channelList_opt, excludeChannelsFromList_opt)
    ! 
    ! :Purpose: subroutine to initialize the chanprof structure used by RTTOV
    !
    implicit none

    ! Arguments:
    integer,              intent(in)  :: sensorHeaderIndexes(:)     ! indexes of radiance observations in header table for the currently processed sensor
    type(struct_obs),     intent(in)  :: obsSpaceData               ! obsSpaceData structure
    type(rttov_chanprof), intent(out) :: chanprof(:)                ! chanprof RTTOV structure
    logical,    optional, intent(out) :: lchannel_subset_opt(:,:)   ! logical array for channel selection by profile for RttovScatt
    integer,    optional, intent(out) :: iptobs_cma_opt(:)          ! list of observation locations in obsSpace Data body table 
    integer,    optional, intent(in)  :: channelList_opt(:)         ! list of channel to select or exclude
    logical,    optional, intent(in)  :: excludeChannelsFromList_opt! .true. to exclude channels from list; .false. to select them  

    ! Locals:
    integer :: btCount, profileIndex, headerIndex, bodyIndexBeg, bodyIndexEnd, bodyIndex, channelNumber
    integer :: channelIndex
    logical :: isChannelInList, excludeChannelsFromList, selected

    ! Build the list of channels/profiles indices

    excludeChannelsFromList = .false.
    if (present(excludeChannelsFromList_opt)) then
      excludeChannelsFromList = excludeChannelsFromList_opt
    end if
    
    if (present(lchannel_subset_opt)) lchannel_subset_opt(:,:) = .false.
    
    btCount = 0
    
    do profileIndex = 1, size(sensorHeaderIndexes)
      headerIndex = sensorHeaderIndexes(profileIndex)
      if (headerIndex <= 0) cycle
      bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1
      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) then
          call tvs_getChannelNumIndexFromPPP( obsSpaceData, headerIndex, bodyIndex, &
              channelNumber, channelIndex )
          if (channelIndex > 0) then
            isChannelInList = .true.
            if (present(channelList_opt) ) then
              isChannelInList = (utl_findloc(channelList_opt, channelNumber) > 0 )
            end if
            if (excludeChannelsFromList) then
              selected = .not. isChannelInList
            else
              selected = isChannelInList
            end if
            if (selected) then
              btCount = btCount + 1
              chanprof(btCount) % prof = profileIndex
              chanprof(btCount) % chan = channelIndex
              if (present(iptobs_cma_opt)) iptobs_cma_opt(btCount) = bodyIndex
              if (present(lchannel_subset_opt)) lchannel_subset_opt(profileIndex,channelIndex) = .true.
            end if
          else
            write(*,*) 'tvs_getChanProf: strange channel number', channelNumber
          end if
        end if
      end do
    end do
  
  end subroutine tvs_getChanprof

  !--------------------------------------------------------------------------
  !  tvs_countRadiances
  !--------------------------------------------------------------------------
  integer function tvs_countRadiances(sensorHeaderIndexes, obsSpaceData)
    !
    ! :Purpose: to count all radiances selected for assimilation
    !
    implicit none

    ! Arguments:
    integer,           intent(in)    :: sensorHeaderIndexes(:) ! indexes in header table of radiance observations for the currently processed sensor
    type(struct_obs),  intent(inout) :: obsSpaceData           ! obsSpaceData structure
    
    ! Locals:
    integer :: profileIndex, headerIndex, bodyIndex

    tvs_countRadiances = 0
    do profileIndex = 1, size(sensorHeaderIndexes)
      headerIndex = sensorHeaderIndexes(profileIndex)
      if (headerIndex <= 0) cycle
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      BODY:do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY
        if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) tvs_countRadiances  = tvs_countRadiances + 1
      end do BODY
    end do

  end function tvs_countRadiances
  
  !--------------------------------------------------------------------------
  !  tvs_countRadiancesScatt
  !--------------------------------------------------------------------------
  integer function tvs_countRadiancesScatt(sensorHeaderIndexes, obsSpaceData, scattChannelList, sensorIndex)
    !
    ! :Purpose: to count radiances selected for assimilation that need to be simulated using RttovScatt
    !
    implicit none

    ! Arguments:
    integer,           intent(in)    :: sensorHeaderIndexes(:) ! indexes of radiance observations for the currently processed sensor
    type(struct_obs),  intent(inout) :: obsSpaceData         ! obsSpaceData structure
    integer,           intent(in)    :: scattChannelList(:)  ! list of channel numbers to process using RttovScatt
    integer,           intent(in)    :: sensorIndex          ! sensor index in NAMTOV namelist section
    
    ! Locals:
    integer :: profileIndex, headerIndex, bodyIndex, channelNumber

    tvs_countRadiancesScatt = 0
    do profileIndex = 1, size(sensorHeaderIndexes)
      headerIndex = sensorHeaderIndexes(profileIndex)
      if (headerIndex <= 0) cycle
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      BODY:do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY
        if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) then
          channelNumber = nint(obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex))
          channelNumber = max(0 , min(channelNumber, tvs_maxChannelNumber + 1))
          channelNumber = channelNumber - tvs_channelOffset(sensorIndex)
          if (utl_findloc(scattChannelList,channelNumber) > 0) then
            tvs_countRadiancesScatt  = tvs_countRadiancesScatt + 1
          end if
        end if
      end do BODY
    end do
  
  end function tvs_countRadiancesScatt

  !--------------------------------------------------------------------------
  !  tvs_ChangedStypValue(obsspacedata, headerIndex)
  !--------------------------------------------------------------------------
  integer function tvs_ChangedStypValue(obsSpaceData, headerIndex)
    !
    ! :Purpose: to obtain new STYP value given observed STYP and TTYP value
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: headerIndex  ! observation index in obsSpaceData header table
    type(struct_obs), intent(in) :: obsSpaceData ! obsSpaceData structure
    
    ! Locals:
    integer :: terrainType
    integer :: landSea 

    terrainType = obs_headElem_i(obsSpaceData,OBS_TTYP,headerIndex)
    landSea     = obs_headElem_i(obsSpaceData,OBS_STYP,headerIndex)

    if ( terrainType ==  0 ) then
      tvs_ChangedStypValue = surftype_seaice
    else
      tvs_ChangedStypValue = landSea
    end if

  end function tvs_ChangedStypValue

  !--------------------------------------------------------------------------
  !  tvs_getHIREmissivities
  !--------------------------------------------------------------------------
  subroutine tvs_getHIREmissivities(sensorHeaderIndexes, obsSpaceData, surfem)
    !
    ! :Purpose: to get emissivity for Hyperspectral Infrared Sounders (AIRS, IASI, CrIS, ...)
    !
    implicit none

    ! Arguments:
    integer,          intent(in)  :: sensorHeaderIndexes(:) ! header indexes of radiance observations for the currently processed sensor
    type(struct_obs), intent(in)  :: obsSpaceData           ! obsSpaceData structure
    real(8),          intent(out) :: surfem(:)              ! surface emissivity

    ! Locals:
    integer :: count, profileIndex, bodyIndexBeg, bodyIndexEnd, bodyIndex, headerIndex

    count = 0 
    surfem(:) = 0.98d0
    do profileIndex = 1, size(sensorHeaderIndexes)
      headerIndex = sensorHeaderIndexes(profileIndex)
      if (headerIndex <= 0 ) cycle
      bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1
      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) then
          count = count + 1
          surfem ( count ) = obs_bodyElem_r(obsSpaceData,OBS_SEM,bodyIndex)
        end if
      end do
    end do

  end subroutine tvs_getHIREmissivities

  !--------------------------------------------------------------------------
  !  tvs_getOtherEmissivities
  !--------------------------------------------------------------------------
  subroutine tvs_getOtherEmissivities(chanprof, sensorHeaderIndexes, sensorType, instrument, surfem, calcemis)
    !
    ! :Purpose: to get emissivity for microwave sounders and infrared geostationary imagers
    !
    implicit none

    ! Arguments:
    type(rttov_chanprof), intent(in)  :: chanprof(:)           ! chanprof RTTOV structure
    integer,              intent(in)  :: sensorHeaderIndexes(:)! indexes of radiance observations for the currently processed sensor
    integer,              intent(in)  :: sensorType            ! RTTOV sensor type
    integer,              intent(in)  :: instrument            ! RTTOV instrument code
    real(8),              intent(out) :: surfem(:)             ! surface emissivity
    logical,              intent(out) :: calcemis(:)           ! flag to request emissivity computation by RTTOV
    
    ! Locals:
    integer :: radianceIndex, profileIndex, headerIndex, surfaceType

    do radianceIndex = 1, size(chanprof)
      profileIndex = chanprof(radianceIndex) % prof
      headerIndex = sensorHeaderIndexes(profileIndex)
      surfaceType = tvs_profiles_nl(headerIndex) % skin % surftype
      if ( sensorType == sensor_id_mw ) then
        if ( surfaceType == surftype_land .or. &
             surfaceType == surftype_seaice     ) then
          calcemis(radianceIndex) = .false.
          surfem (radianceIndex) = 0.75d0
        else
          calcemis(radianceIndex) = .true.
          surfem (radianceIndex) = 0.d0
        end if
      else if ( tvs_isInstrumHyperSpectral(instrument) ) then
        calcemis(radianceIndex) = .false. 
      else if ( tvs_isInstrumGeostationary(instrument) ) then
        calcemis(radianceIndex) = .true.
        surfem (radianceIndex) = 0.d0
      else
        write(*,*) sensorType,instrument
        call utl_abort('tvs_getOtherEmissivities. invalid sensor type or unknown IR instrument')
      end if
    end do
   
  end subroutine tvs_getOtherEmissivities
  
  !--------------------------------------------------------------------------
  !  tvs_fillProfiles
  !--------------------------------------------------------------------------
  subroutine tvs_fillProfiles(columnTrl, obsSpaceData, datestamp, profileType, beSilent)
    !
    ! :Purpose:  to fill in tvs_profiles_nl structure before call to non-linear, 
    !            tangent-linear or adjoint of RTTOV
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrl    ! trial field column structure
    type(struct_obs),        intent(inout) :: obsSpaceData ! obsSpaceData structure
    integer,                 intent(in)    :: datestamp    ! CMC date stamp
    character(len=*),        intent(in)    :: profileType  ! profile type (could be 'nl' or 'tlad')
    logical,                 intent(in)    :: beSilent     ! verbosity flag

    ! Locals:
    integer :: instrum, iplatform
    integer :: sensorIndex, headmax
    integer :: profileCount, headerIndex
    integer :: profileIndex, levelIndex
    integer :: ilowlvl_M,ilowlvl_T,nlv_M,nlv_T
    integer :: Vcode
    integer :: ierr,day,month,year,ijour,itime
    integer :: allocStatus    
    integer, external ::  newdate
    integer, allocatable :: sensorHeaderIndexes(:)  
    type(struct_vco), pointer :: vco
    real(8), allocatable :: pressure (:,:)
    real(8), allocatable :: height(:,:)
    real(8), allocatable :: latitudes(:)
    real(8), allocatable :: ozone(:,:)
    character(len=4)     :: ozoneVarName
    real(8), allocatable :: clw   (:,:)
    real(8), allocatable :: ciw   (:,:)
    real(8), allocatable :: rainflux  (:,:)
    real(8), allocatable :: snowflux  (:,:)
    real(8), allocatable :: cloudFraction(:,:)
    logical, allocatable :: surfTypeIsWater(:)
    logical :: runObsOperatorWithClw
    logical :: runObsOperatorWithHydrometeors
    type(rttov_profile), pointer :: profiles(:)
    type(rttov_profile_cloud), pointer :: cld_profiles(:)
    real(8), pointer :: column_ptr(:),column_ptrHU(:)
    real(8) :: zmax, wind

    if ( .not. beSilent ) write(*,*) 'tvs_fillProfiles: Starting'
  
    if (tvs_headerEnd < 0) return    ! exit if there are no tovs data

    if ( tvs_numMWInstrumUsingCLW > 0 .and. .not. col_varExist(columnTrl,'LWCR') ) then
      call utl_abort('tvs_fillProfiles: if number of instrument to use CLW greater than zero, ' // &
                     'the LWCR variable must be included as an analysis variable in NAMSTATE. ')
    end if
    if (tvs_numMWInstrumUsingHydrometeors > 0) then
      if (.not. (col_varExist(columnTrl,'LWCR') .and. col_varExist(columnTrl,'IWCR'))) then
        call utl_abort('tvs_fillProfiles: if number of instrument to use hydrometeors greater than zero, ' // &
                       'the LWCR/IWCR variables must be included as an analysis variable in NAMSTATE. ')
      end if
      if (.not. beSilent .and. .not. (col_varExist(columnTrl,'RF') .and. col_varExist(columnTrl,'SF') .and. &
          col_varExist(columnTrl,'CLDR'))) then
        write(*,*) 'tvs_fillProfiles: one of RF/SF/CLDR does not exist in NAMSTATE'
      end if
    end if

    if ( (tvs_numMWInstrumUsingCLW == 0 .and. tvs_numMWInstrumUsingHydrometeors == 0 .and. &
            tvs_mwAllskyAssim) .or. &
         (tvs_numMWInstrumUsingCLW  > 0 .and. tvs_numMWInstrumUsingHydrometeors == 0 .and. &
            .not. tvs_mwAllskyAssim) .or. &
         (tvs_numMWInstrumUsingCLW == 0 .and. tvs_numMWInstrumUsingHydrometeors  > 0 .and. &
            .not. tvs_mwAllskyAssim) ) then
      call utl_abort('tvs_fillProfiles: number of instrument to use CLW/hydrometeors do not match ' // &
                     'all-sky namelist variable.')
    end if

    if (tvs_useO3FromTrials .and. .not. col_varExist(columnTrl,'TO3') .and. &
        .not. col_varExist(columnTrl,'O3L') ) then
      call utl_abort('tvs_fillProfiles: if tvs_useO3FromTrials is set to .true. the ozone variable ' // &
                     'must be included as an analysis variable in NAMSTATE. ')
    else if (tvs_useO3FromTrials) then 
      if (col_varExist(columnTrl,'TO3') ) then
        ozoneVarName = 'TO3'
      else
        ozoneVarName = 'O3L'
      end if
    end if

    if ( profileType == 'nl' ) then
      if ( .not. allocated( tvs_profiles_nl) ) then
        allocate(tvs_profiles_nl(tvs_headerEnd))
        if (tvs_numMWInstrumUsingHydrometeors > 0) then
          allocate(tvs_cld_profiles_nl(tvs_headerEnd))
        end if
      end if
    else if ( profileType == 'tlad' ) then
      if ( .not. allocated(tvs_profiles_tlad) ) then
        allocate(tvs_profiles_tlad(tvs_headerEnd))
        if (tvs_numMWInstrumUsingHydrometeors > 0) then
          allocate(tvs_cld_profiles_tlad(tvs_headerEnd))
        end if
      else
        return
      end if
    else
      write(*,*) 'Invalid  profileType ', profileType
      call utl_abort('tvs_fillProfiles')
    end if

    if ( .not. beSilent ) write(*,*) 'tvs_fillProfiles: profileType is ', profileType

    call tvs_getProfile(profiles, profileType, cld_profiles)

    !
    !     1.    Set index for model's lowest level and model top
    !     .     ------------------------------------------------
    
    nlv_M = col_getNumLev(columnTrl,'MM')
    nlv_T = col_getNumLev(columnTrl,'TH')

    if (  col_getPressure(columnTrl,1,1,'TH') < col_getPressure(columnTrl,nlv_T,1,'TH') ) then
      ilowlvl_M = nlv_M
      ilowlvl_T = nlv_T
    else
      ilowlvl_M = 1
      ilowlvl_T = 1
    end if

    vco => col_getVco(columnTrl)
    Vcode = vco % Vcode

    ierr = newdate(datestamp,ijour,itime,-3)
    if (ierr < 0) then
      write(*,*) 'Invalid datestamp ',datestamp,ijour,itime,ierr
      call utl_abort('tvs_fillProfiles')
    end if
    year= ijour / 10000
    month = mod(ijour / 100,100)
    day = mod(ijour,100)

    !  1.2   Read ozone climatology

    if (.not. tvs_useO3FromTrials) call clm_readFields()
   
    !     2.  Fill profiles structure
    
    ! loop over all instruments
    sensor_loop: do sensorIndex=1, tvs_nsensors

      runObsOperatorWithClw = (col_varExist(columnTrl,'LWCR') .and. tvs_numMWInstrumUsingCLW /= 0 .and. & 
                               tvs_isInstrumUsingCLW(tvs_instruments(sensorIndex)))

      runObsOperatorWithHydrometeors = (col_varExist(columnTrl,'LWCR') .and. col_varExist(columnTrl,'IWCR') .and. &
                                        tvs_isInstrumUsingHydrometeors(tvs_instruments(sensorIndex)))

      ! first loop over all obs.
      bobs1: do headerIndex = 1, tvs_headerEnd
        if (tvs_lsensor(headerIndex) == sensorIndex) then
          headmax = headerIndex
        end if
      end do bobs1

      profileCount = tvs_countProfiles(sensorIndex)

      if (profileCount == 0) cycle sensor_loop

      if (tvs_coefs(sensorIndex) % coef % fmv_model_ver >= 9) then
        zmax = zenmaxv9
      else
        zmax = zenmax
      end if
      
      allocate(sensorHeaderIndexes(profileCount))
      sensorHeaderIndexes(:) = -1
      allocate(latitudes(profileCount))
      allocate(ozone(nlv_T,profileCount)) 
      allocate(pressure(nlv_T,profileCount))

      if (runObsOperatorWithClw .or. runObsOperatorWithHydrometeors) then
        allocate(clw(nlv_T,profileCount))
        clw(:,:) = qlim_getMinValueCloud('LWCR')
      end if
      if (runObsOperatorWithHydrometeors) then
        allocate(ciw          (nlv_T,profileCount))
        allocate(rainFlux     (nlv_T,profileCount))
        allocate(snowFlux     (nlv_T,profileCount))
        allocate(cloudFraction(nlv_T,profileCount))
        ciw(:,:) = qlim_getMinValueCloud('IWCR')
        rainFlux(:,:) = qlim_getMinValueCloud('RF')
        snowFlux(:,:) = qlim_getMinValueCloud('SF')
        cloudFraction(:,:) = qlim_getMinValueCloud('CLDR')
      end if
      allocate(surfTypeIsWater(profileCount)) 
      surfTypeIsWater(:) = .false.
      
      allocate (height(nlv_T,profileCount))

      profileCount = 0
      ! second loop over all obs.
      bobs2: do headerIndex = 1, headmax
        if (tvs_lsensor(headerIndex) /= sensorIndex) cycle bobs2
        profileCount = profileCount + 1
        sensorHeaderIndexes(profileCount) = headerIndex

        call rttov_alloc_prof(                 &
             allocStatus,                      &
             1,                                & ! 1 = nprofiles un profil a la fois
             profiles(headerIndex:headerIndex),&
             nlv_T,                            & 
             tvs_opts(sensorIndex),            &
             asw=1,                            & ! asw =1 allocation
             coefs=tvs_coefs(sensorIndex),     &
             init=.true. )
        if (allocStatus /= 0) call utl_abort('tvs_fillProfiles: memory allocation error in rttov_alloc_prof')
        if (runObsOperatorWithHydrometeors) then
          call rttov_alloc_scatt_prof(                &   
               allocstatus,                           &
               1,                                     &
               cld_profiles(headerIndex:headerIndex), &
               nlv_T,                                 &
               nhydro=5,                              & ! depending on what is defined in the Mie tables
               nhydro_frac=1,                         & ! 1 cloud fraction for all variable or nhydro 1 cloud fraction for each variable
               asw=1_jpim,                            & ! 1 => allocate
               init=.true.,                           & ! initialize profiles to zero
               flux_conversion=[1,2,0,0,0] )            !flux_conversion  input units: 0 (default) => kg/kg,
                                                        ! 1,2 => kg/m2/s, optional for rain, snow
          if (allocStatus /= 0) call utl_abort('tvs_fillProfiles: memory allocation error in rttov_alloc_scatt_prof')
        end if

        !    extract land/sea/sea-ice flag (0=land, 1=sea, 2=sea-ice)
        profiles(headerIndex) % skin % surftype = tvs_ChangedStypValue(obsSpaceData,headerIndex)

        !    extract satellite zenith and azimuth angle, 
        !    sun zenith angle, cloud fraction, latitude and longitude
        profiles(headerIndex) % zenangle   = obs_headElem_r(obsSpaceData,OBS_SZA,headerIndex)

        call validateRttovVariable(profiles(headerIndex) % zenangle, "satellite zenith angle", &
            obsSpaceData, headerIndex, 0.d0, zmax)
 
        profiles(headerIndex) % azangle = tvs_getCorrectedSatelliteAzimuth(obsSpaceData, headerIndex)
        profiles(headerIndex) % sunazangle  = obs_headElem_r(obsSpaceData,OBS_SAZ,headerIndex) ! necessaire pour radiation solaire
        iplatform = tvs_coefs(sensorIndex) % coef % id_platform
        instrum = tvs_coefs(sensorIndex) % coef % id_inst
        profiles(headerIndex) % sunzenangle = obs_headElem_r(obsSpaceData,OBS_SUN,headerIndex)
        if (tvs_opts(sensorIndex) % rt_ir % addsolar) then
          call validateRttovVariable(profiles(headerIndex) % sunzenangle, "sun zenith angle", &
              obsSpaceData, headerIndex, 0.d0)
        end if
        latitudes(profileCount) = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex) * MPC_DEGREES_PER_RADIAN_R8
        !1d-5 was chosen as a threshold because it is the high precision for latitudes in BUFR/BURP files
        if (latitudes(profileCount) > 90.d0 .and. (latitudes(profileCount)-90.d0) < 1.d-5) latitudes(profileCount) = 90.d0
        if (latitudes(profileCount) < -90.d0 .and. (-latitudes(profileCount)-90.d0) < 1.d-5) latitudes(profileCount) = -90.d0
        call validateRttovVariable(latitudes(profileCount), 'latitude', obsSpaceData, headerIndex, -90.d0, 90.d0) 
        
        profiles(headerIndex) % longitude = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex) * MPC_DEGREES_PER_RADIAN_R8

        surfTypeIsWater(profileCount) = ( tvs_ChangedStypValue(obsSpaceData,headerIndex) == surftype_sea )

        do levelIndex = 1, nlv_T
          pressure(levelIndex,profileCount) = col_getPressure(columnTrl,levelIndex,headerIndex,'TH') * MPC_MBAR_PER_PA_R8
          height(levelIndex,profileCount) = col_getHeight(columnTrl,levelIndex,headerIndex,'TH') ! in meters
          if ((runObsOperatorWithClw .and. surfTypeIsWater(profileCount)) .or. &
              (runObsOperatorWithHydrometeors .and. surfTypeIsWater(profileCount))) then

            ! cloud liquid water content
            clw(levelIndex,profileCount) = col_getElem(columnTrl,levelIndex,headerIndex,'LWCR')
            if (clw(levelIndex,profileCount) < qlim_getMinValueCloud('LWCR') .or. &
                clw(levelIndex,profileCount) > qlim_getMaxValueCloud('LWCR')) then
              write(*,*) 'tvs_fillProfiles: clw=' , clw(:,profileCount) 
              call utl_abort('tvs_fillProfiles: columnTrl has clw outside RTTOV bounds')
            end if

            clw(levelIndex,profileCount) = clw(levelIndex,profileCount) * tvs_cloudScaleFactor
          end if
          if (runObsOperatorWithHydrometeors .and. surfTypeIsWater(profileCount)) then
            ! cloud ice water content
            ciw(levelIndex,profileCount) = col_getElem(columnTrl,levelIndex,headerIndex,'IWCR')
            if (ciw(levelIndex,profileCount) < qlim_getMinValueCloud('IWCR') .or. &
                ciw(levelIndex,profileCount) > qlim_getMaxValueCloud('IWCR')) then
              write(*,*) 'tvs_fillProfiles: ciw=' , ciw(:,profileCount) 
              call utl_abort('tvs_fillProfiles: columnTrl has ciw outside RTTOV bounds')
            end if

            ! rain flux (zero, if not part of control variables)
            if (col_varExist(columnTrl,'RF')) then
              rainFlux(levelIndex,profileCount) = col_getElem(columnTrl,levelIndex,headerIndex,'RF')
            end if
            if (rainFlux(levelIndex,profileCount) < qlim_getMinValueCloud('RF') .or. &
                rainFlux(levelIndex,profileCount) > qlim_getMaxValueCloud('RF')) then
              write(*,*) 'tvs_fillProfiles: rainFlux=' , rainFlux(:,profileCount) 
              call utl_abort('tvs_fillProfiles: columnTrl has rain flux outside RTTOV bounds')
            end if

            ! snow flux (zero, if not part of control variables)
            if (col_varExist(columnTrl,'SF')) then
              snowFlux(levelIndex,profileCount) = col_getElem(columnTrl,levelIndex,headerIndex,'SF')
            end if
            if (snowFlux(levelIndex,profileCount) < qlim_getMinValueCloud('SF') .or. &
                snowFlux(levelIndex,profileCount) > qlim_getMaxValueCloud('SF')) then
              write(*,*) 'tvs_fillProfiles: snowFlux=' , snowFlux(:,profileCount) 
              call utl_abort('tvs_fillProfiles: columnTrl has snow flux outside RTTOV bounds')
            end if

            ! cloud fraction (zero, if not part of control variables)
            if (col_varExist(columnTrl,'CLDR')) then
              cloudFraction(levelIndex,profileCount) = col_getElem(columnTrl,levelIndex,headerIndex,'CLDR')
            end if
            if (cloudFraction(levelIndex,profileCount) < qlim_getMinValueCloud('CLDR') .or. &
                cloudFraction(levelIndex,profileCount) > qlim_getMaxValueCloud('CLDR')) then
              write(*,*) 'tvs_fillProfiles: cloudFraction=' , cloudFraction(:,profileCount) 
              call utl_abort('tvs_fillProfiles: columnTrl has cloud fraction outside RTTOV bounds')
            end if

            ciw(levelIndex,profileCount) = ciw(levelIndex,profileCount) * tvs_cloudScaleFactor
            rainFlux(levelIndex,profileCount) = rainFlux(levelIndex,profileCount) * tvs_cloudScaleFactor
            snowFlux(levelIndex,profileCount) = snowFlux(levelIndex,profileCount) * tvs_cloudScaleFactor

          end if ! runObsOperatorWithHydrometeors .and. surfTypeIsWater
        end do ! levelIndex
        
        ! Constituents assumed to be in micrograms/kg
        if (tvs_coefs(sensorIndex) % coef % nozone > 0 .and. tvs_useO3FromTrials) then
          ! Get ozone from trial field
          column_ptr => col_getColumn(columnTrl, headerIndex, trim(ozoneVarName) )
          ozone(:,profileCount) = column_ptr(:)
        else if (tvs_coefs(sensorIndex) % coef % nozone > 0 .and. .not. tvs_useO3FromTrials) then
          ! Get ozone profiles (ug/kg) from climatology
          column_ptr => col_getColumn(columnTrl, headerIndex,'TT' )
          column_ptrHU => col_getColumn(columnTrl, headerIndex,'HU' )
          call clm_setColumn(nlv_T,pressure(:,profileCount)*MPC_PA_PER_MBAR_R8, &
                             height(:,profileCount),latitudes(profileCount), &
                             obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)*MPC_DEGREES_PER_RADIAN_R8, &
                             profileCount,maxsize,BUFR_NECH_O3,tt_opt=column_ptr,hu_opt=column_ptrHU, &
                             climatProfile_opt=ozone(:,profileCount))
        end if

      end do bobs2

      !   2.5  Fill profiles structure

      do  profileIndex = 1 , profileCount 
        headerIndex = sensorHeaderIndexes(profileIndex)
        profiles(headerIndex) % gas_units       = gas_unit_specconc ! all gas profiles are supposed to be provided in kg/kg (specific humidity, i.e. mass mixing ratio [kg/kg] over wet air)
        profiles(headerIndex) % id              = '' ! profile id, up to 128 characters, to consider for use
        profiles(headerIndex) % nlevels         = nlv_T
        profiles(headerIndex) % nlayers         = nlv_T - 1
        profiles(headerIndex) % date(1)         = year
        profiles(headerIndex) % date(2)         = month
        profiles(headerIndex) % date(3)         = day
        profiles(headerIndex) % latitude        = latitudes(profileIndex)
        profiles(headerIndex) % elevation       = 0.001d0 * col_getHeight(columnTrl,ilowlvl_T,headerIndex,'TH') ! unite km
        call validateRttovVariable(profiles(headerIndex) % elevation, "elevation", obsSpaceData, headerIndex, maximum_opt=elevmax)
        profiles(headerIndex) % skin % watertype = watertype_ocean_water !utilise pour calcul rayonnement solaire reflechi seulement
        profiles(headerIndex) % skin % t         = col_getElem(columnTrl,1,headerIndex,'TG')
        call validateRttovVariable( profiles(headerIndex) % skin % t, "skin temperature", obsSpaceData, headerIndex, tmin, tmax) 
        profiles(headerIndex) % skin % salinity  = 35.d0 ! for FASTEM-4 only to revise (practical salinity units)
        profiles(headerIndex) % skin % fastem(:) = 0.0d0
        profiles(headerIndex) % skin % snow_fraction  = 0.d0 ! Surface coverage snow fraction(0-1), used only by IR emissivity atlas
        profiles(headerIndex) % skin % soil_moisture  = 0.d0 ! soil moisure (m**3/m**3) not yet used
        profiles(headerIndex) % s2m % t          = col_getElem(columnTrl,ilowlvl_T,headerIndex,'TT')
        call validateRttovVariable(profiles(headerIndex) % s2m % t, '2m air temperature', &
            obsSpaceData, headerIndex, tmin, tmax) 
        profiles(headerIndex) % s2m % p         = col_getElem(columnTrl,1      ,headerIndex,'P0')*MPC_MBAR_PER_PA_R8
        call validateRttovVariable(profiles(headerIndex) % s2m % p, 'surface pressure', &
            obsSpaceData, headerIndex, pmin, pmax) 
        profiles(headerIndex) % s2m % u         = col_getElem(columnTrl,ilowlvl_M,headerIndex,'UU')
        profiles(headerIndex) % s2m % v         = col_getElem(columnTrl,ilowlvl_M,headerIndex,'VV')
        wind = sqrt( profiles(headerIndex) % s2m % u ** 2 + &
            profiles(headerIndex) % s2m % v ** 2 )
        if ( wind > wmax ) then
          write(*,*) 'tvs_fillProfiles: !!! WARNING !!!'
          write(*,*) 'tvs_fillProfiles: INVALID 10m wind speed'
          write(*,*) 'tvs_fillProfiles: headerIndex ', headerIndex, " !"
          write(*,*) 'tvs_fillProfiles: modulus ', wind, ' larger than ', wmax, 'set to zero !'
          profiles(headerIndex) % s2m % u = 0.d0
          profiles(headerIndex) % s2m % v = 0.d0
          call rejectObs(obsSpaceData, headerIndex)
        end if
        profiles(headerIndex) % s2m % o         = 0.0d0 !surface ozone never used
        profiles(headerIndex) % s2m % wfetc     = 100000.0d0 ! Wind fetch (in meter for rttov10 ?) used to calculate reflection of solar radiation by sea surface
        profiles(headerIndex) % icede_param     = 0
        profiles(headerIndex) % Be              = 0.4d0 ! earth magnetic field strength (gauss) (must be non zero)
        profiles(headerIndex) % cosbk           = 0.0d0 ! cosine of the angle between the earth magnetic field and wave propagation direction
        profiles(headerIndex) % p(:)            = pressure(:,profileIndex)
        call validateRttovVariable(profiles(headerIndex) % p(nlv_T), "pressure profile near surface", obsSpaceData, headerIndex, maximum_opt=2000.d0)
        call validateRttovVariable(profiles(headerIndex) % p(1), "pressure profile near top", obsSpaceData, headerIndex, 0.d0) 
        !RTTOV scatt needs half pressure levels (see figure 5 of RTTOV 12 User's Guide)
        if (runObsOperatorWithHydrometeors) then
          cld_profiles(headerIndex) % ph (1) = 0.d0
          cld_profiles(headerIndex) % cfrac = 0.d0
          do levelIndex = 1, nlv_T - 1
            cld_profiles(headerIndex) % ph (levelIndex+1) = 0.5d0 * (profiles(headerIndex) % p(levelIndex) + &
                                                                   profiles(headerIndex) % p(levelIndex+1))
          end do
          cld_profiles(headerIndex) % ph (nlv_T+1) = profiles(headerIndex) % s2m % p
        end if
        column_ptr => col_getColumn(columnTrl, headerIndex,'TT' )
        profiles(headerIndex) % t(:)   = column_ptr(:)
        call validateRttovProfile(profiles(headerIndex) % t, 'temperature', tmin, tmax, obsSpaceData, headerIndex) 
        if (tvs_coefs(sensorIndex) % coef % nozone > 0) then
          profiles(headerIndex) % o3(:) = ozone(:,profileIndex) * microg2kg ! micrograms/kg to kg/kg
          call validateRttovProfile(profiles(headerIndex) % o3, 'ozone', o3min, o3max, obsSpaceData, headerIndex)
        end if

        column_ptr => col_getColumn(columnTrl, headerIndex,'HU' )
        profiles(headerIndex) % q(:)            =  column_ptr(:)
        call validateRttovProfile(profiles(headerIndex) % q, 'water vapor', qmin, qmax, obsSpaceData, headerIndex)
        profiles(headerIndex) % ctp = 1013.25d0
        profiles(headerIndex) % cfraction = 0.d0
        if (runObsOperatorWithClw) profiles(headerIndex) % clw(:) = clw(:,profileIndex)
        if (runObsOperatorWithHydrometeors) then
          cld_profiles(headerIndex) % hydro(:,1) = rainFlux(:,profileIndex)
          cld_profiles(headerIndex) % hydro(:,2) = snowFlux(:,profileIndex)
          cld_profiles(headerIndex) % hydro(:,4) = clw(:,profileIndex)
          cld_profiles(headerIndex) % hydro(:,5) = ciw(:,profileIndex)

          do levelIndex = 1, nlv_T
            if (cld_profiles(headerIndex) % hydro(levelIndex,1) > qlim_getMinValueCloud('RF') .or. &
                cld_profiles(headerIndex) % hydro(levelIndex,2) > qlim_getMinValueCloud('SF') .or. &
                cld_profiles(headerIndex) % hydro(levelIndex,4) > qlim_getMinValueCloud('LWCR') .or. &
                cld_profiles(headerIndex) % hydro(levelIndex,5) > qlim_getMinValueCloud('IWCR')) then
              
              ! set to overcast cloud, if CLDR not part of control variables
              if (col_varExist(columnTrl,'CLDR')) then
                cld_profiles(headerIndex) % hydro_frac(levelIndex,1) = cloudFraction(levelIndex,profileIndex)
              else
                cld_profiles(headerIndex) % hydro_frac(levelIndex,1) = 1.0d0
              end if
            else
              cld_profiles(headerIndex) % hydro_frac(levelIndex,1) = qlim_getMinValueCloud('CLDR')  
            end if

          end do ! levelIndex
        end if  ! runObsOperatorWithHydrometeors
      end do ! profileIndex

      ! Extract emissivity from background column object to be used in the computation
      ! of non-linear RTTOV
      if (col_varExist(columnTrl, 'EMMW')) then
        call sse_extractEmissivityCol(columnTrl, tvs_emissivityFromTrl, profileCount, sensorHeaderIndexes, tvs_headerEnd)
      end if

      deallocate(pressure)
      if (allocated(ozone)) deallocate(ozone)
      deallocate(latitudes)
      deallocate(sensorHeaderIndexes)
      if (allocated(clw)) then
        deallocate(clw)
      end if
      if (allocated(ciw)) then
        deallocate(ciw)
        deallocate(rainFlux)
        deallocate(snowFlux)
        deallocate(cloudFraction)
      end if
      deallocate(surfTypeIsWater)
      deallocate (height)
     
    end do sensor_loop

  end subroutine tvs_fillProfiles

  !--------------------------------------------------------------------------
  !  tvs_getCorrectedSatelliteAzimuth
  !--------------------------------------------------------------------------
  function tvs_getCorrectedSatelliteAzimuth(obsSpaceData, headerIndex) result(correctedAzimuth)
    !
    ! :Purpose: get properly corrected satellite Azimuth Angle from obsSpaceData header
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(in) :: obsSpaceData     ! obsSpaceData structure
    integer,          intent(in) :: headerIndex      ! observation index in obsSpaceData header table
    ! Result:
    real(8)                      :: correctedAzimuth ! corrected azimuth (function result)

    ! Locals:
    integer :: sensorNo

    correctedAzimuth = obs_headElem_r(obsSpaceData,OBS_AZA,headerIndex)

    sensorNo = tvs_lsensor(headerIndex)

    if (.not. tvs_isAzimuthValid(sensorNo)) then
      correctedAzimuth = obs_missingValue_R
      return
    end if

    if ( tvs_doAzimuthCorrection(sensorNo) ) then
      ! Correction sur la definition de l'angle.
      correctedAzimuth = obs_headElem_r(obsSpaceData,OBS_SAZ,headerIndex) + correctedAzimuth
      if ( correctedAzimuth > 360.d0 ) correctedAzimuth = correctedAzimuth - 360.d0
    end if

  end function tvs_getCorrectedSatelliteAzimuth

  !--------------------------------------------------------------------------
  !  tvs_rttov
  !--------------------------------------------------------------------------
  subroutine tvs_rttov(obsSpaceData, bgckMode, beSilent)
    !
    ! :Purpose: Interface for RTTOV non linear operator
    !           tvs_fillProfiles should be called before
    !
    implicit none

    ! Arguments:
    type(struct_obs),  intent(inout) :: obsSpaceData          ! obsSpaceData structure
    logical,           intent(in)    :: bgckMode              ! flag to transfer transmittances and cloudy overcast radiances in bgck mode 
    logical,           intent(in)    :: beSilent              ! verbosity flag

    ! Locals:
    integer :: nlv_T
    integer :: btCount
    integer :: btCountScatt
    integer :: allocStatus
    integer :: rttov_err_stat ! rttov error return code
    integer :: nthreads,max_nthreads
    integer :: sensorIndex
    integer :: channelIndex
    integer :: hydroChannelsCount
    integer :: profileCount
    integer :: profileIndex, levelIndex, btIndex, headerIndex
    integer :: instrum
    integer :: sensorType        ! sensor type (1=infrared; 2=microwave; 3=high resolution; 4=polarimetric)
    integer, allocatable :: sensorHeaderIndexes(:)
    type (rttov_emissivity), pointer :: emissivity_local(:)      ! emissivity structure with input and output
    type (rttov_emissivity), pointer :: emissivity_localScatt(:) ! emissivity structure with input and output
    type (rttov_chanprof), allocatable :: chanprof1(:)
    type (rttov_radiance) :: radiancedata_d, radiancedata_d1, radiancedata_dScatt
    type (rttov_transmission) :: transmission, transmission1
    logical, pointer :: calcemis(:), calcemisScatt(:)
    real(8), allocatable  :: surfem1(:), surfem1Scatt(:)
    integer, allocatable  :: frequencies(:)
    logical, allocatable  :: lchannelSubset(:,:)
    integer               :: profileIndex2, tb1, tb2
    integer :: bodyIndex
    real(8) :: clearMwRadiance
    logical :: runObsOperatorWithClw
    logical :: runObsOperatorWithHydrometeors

    if (tvs_headerEnd < 0) return       ! exit if there are not tovs data
    
    if (.not. beSilent) write(*,*) 'tvs_rttov: Starting'
    if (.not. beSilent) call msg_memUsage('tvs_rttov')

    !   1.  Get number of threads available and allocate memory for some variables

    max_nthreads = mmpi_numThread

    !   1.1   Read surface information
    if (bgckMode) call tvs_emis_read_climatology

    !   2.  Computation of hx for tovs data only

    ! Loop over all sensors specified by user
    sensor_loop:do sensorIndex = 1, tvs_nsensors
      
      runObsOperatorWithClw = (tvs_numMWInstrumUsingCLW /= 0 .and. &
          tvs_isInstrumUsingCLW(tvs_instruments(sensorIndex)) )
          
      runObsOperatorWithHydrometeors = (tvs_numMWInstrumUsingHydrometeors /= 0 .and. &
          tvs_isInstrumUsingHydrometeors(tvs_instruments(sensorIndex)) ) 

      sensorType = tvs_coefs(sensorIndex) % coef % id_sensor
      instrum = tvs_coefs(sensorIndex) % coef % id_inst
      profileCount = tvs_countProfiles(sensorIndex)
      call utl_reallocate(sensorHeaderIndexes, profileCount)
      call tvs_setupPointers(runObsOperatorWithHydrometeors, sensorIndex, btCount, btCountScatt, &
          hydroChannelsCount, sensorHeaderIndexes, lChannelSubset, obsSpaceData, &
          irBgckMode_opt= (bgckMode .and. tvs_isInstrumHyperSpectral(instrum) ) )

      if (profileCount == 0) cycle sensor_loop
      if (btCount == 0 .and. btCountScatt==0) cycle sensor_loop
      
      !  loop over all obs.
      obs_loop: do headerIndex = 1, tvs_headerEnd 
        !    Currently processed sensor?
        if (tvs_lsensor(headerIndex) == sensorIndex) then
          nlv_T = tvs_profiles_nl(headerIndex) % nlevels
          exit obs_loop
        end if
      end do obs_loop
      
      !    2.1  Calculate the actual number of threads which will be used.
      
      nthreads = min(max_nthreads, profileCount )  
      
      !    2.2  Prepare all input variables required by rttov.

      if (btCountScatt > 0 .and. (tvs_useSfcEmissObsSpace .or. allocated(tvs_emissivityFromTrl))) then 
        call utl_abort('tvs_rttov_tl: RTTOV scatt does not support the inclusion of surface emissivity in the analysis variable or read from ObsSpaceData')
      end if

      if (btCount > 0) then
        call rttov_alloc_direct(         &
            allocStatus,                 &
            asw=1,                       &
            nprofiles=profileCount,      & ! (not used)
            nchanprof=btCount,           &
            nlevels=nlv_T,               &
            opts=tvs_opts(sensorIndex),  &
            coefs=tvs_coefs(sensorIndex),&
            transmission=transmission,   &
            radiance=radiancedata_d,     &
            calcemis=calcemis,           &
            emissivity=emissivity_local, &
            init=.true.)
        if (allocStatus /= 0) call utl_abort('tvs_rttov: memory allocation error 1 in rttov_alloc_direct')
        allocate(surfem1(btCount))
        
        if (tvs_isInstrumHyperSpectral(instrum)) then
          !     get Hyperspectral IR emissivities
          surfem1(:) = 0.
          if (bgckMode) then
            call emis_getIrEmissivity (surfem1,tvs_nchan(sensorIndex),sensorIndex,profileCount,btCount,sensorHeaderIndexes)
          else
            call tvs_getHIREmissivities(sensorHeaderIndexes, obsSpaceData, surfem1)
          end if
        end if
        
        call tvs_getOtherEmissivities(tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes, sensorType, instrum, surfem1, calcemis)
        
        if (sensorType == sensor_id_mw) then
          if (allocated(tvs_emissivityFromTrl)) then
            ! Read surface emissivity from column when it's included as an analysis variable

            ! Set the default surface emissivity values
            emissivity_local(:) % emis_in = surfem1(:)

            ! Setup the surface emissvity from column object to rttov emissivity_local
            call sse_setupEmissivityfromState(emissivity_local, obsSpaceData, tvs_bodyIndexFromBtIndex(:,sensorIndex), tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes, &
                                              tvs_nsensors, tvs_lsensor, tvs_instrumentName, &
                                              tvs_maxChannelNumber, tvs_channelOffset, tvs_ichan, tvs_profiles_nl(:) % skin % surftype, &
                                              emissivityProfDt_opt = tvs_emissivityFromTrl)
          else if (tvs_useSfcEmissObsSpace) then

            ! Set the default surface emissivity values
            emissivity_local(:) % emis_in = surfem1(:)

            ! Setup the surface emissvity from obsSpaceData Object 
            call sse_emissFromObsSpace(obsSpaceData, emissivity_local, tvs_bodyIndexFromBtIndex(:,sensorIndex), tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes)
          else
            ! Read surface emissivity from emissivity atlas
            call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorIndex, tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes)
          end if
        else
          emissivity_local(:) % emis_in = surfem1(:)
        end if
        
        !   2.3  Compute radiance with rttov_direct

        rttov_err_stat = 0 

        if (bgckMode .and. tvs_isInstrumHyperSpectral(instrum)) then
          write(*,*) 'for bgck IR: call rttov_parallel_direct for each profile...'

          ! allocate transmitance structure for 1 profile
          call rttov_alloc_transmission(allocStatus, transmission1, nlevels=nlv_T, &
              nchanprof=tvs_nchan(sensorIndex), asw=1, init=.true.)
          if (allocStatus /= 0) call utl_abort('tvs_rttov: memory allocation error in rttov_alloc_transmission')
          ! allocate radiance structure for 1 profile
          call rttov_alloc_rad (allocStatus, tvs_nchan(sensorIndex), radiancedata_d1, nlv_T, asw=1, init=.true.)
          if (allocStatus /= 0) call utl_abort('tvs_rttov: memory allocation error in rttov_alloc_rad')
          ! allocate chanprof for 1 profile
          allocate(chanprof1(tvs_nchan(sensorIndex)))
          do channelIndex = 1, tvs_nchan(sensorIndex)
            chanprof1(channelIndex) % prof = 1
            chanprof1(channelIndex) % chan = channelIndex
          end do
          
          do profileIndex2 = 1, profileCount
            tb1 = 1 + (profileIndex2-1) * tvs_nchan(sensorIndex) 
            tb2 = profileIndex2 * tvs_nchan(sensorIndex)

            call rttov_parallel_direct(                                                                &
                rttov_err_stat,                                                                        & ! out
                chanprof1,                                                                             & ! in
                tvs_opts(sensorIndex),                                                                 & ! in
                tvs_profiles_nl(sensorHeaderIndexes(profileIndex2):sensorHeaderIndexes(profileIndex2)),& ! in
                tvs_coefs(sensorIndex),                                                                & ! in
                transmission1,                                                                         & ! inout
                radiancedata_d1,                                                                       & ! inout
                calcemis=calcemis(tb1:tb2),                                                            & ! in
                emissivity=emissivity_local(tb1:tb2),                                                  & ! inout
                nthreads=nthreads )

            ! copy contents of single profile structures into complete structures
            transmission % tau_total(tb1:tb2)             = transmission1 % tau_total(:)
            transmission % tau_levels(:,tb1:tb2)          = transmission1 % tau_levels(:,:)
            transmission % tausun_levels_path1(:,tb1:tb2) = transmission1 % tausun_levels_path1(:,:)
            transmission % tausun_levels_path2(:,tb1:tb2) = transmission1 % tausun_levels_path2(:,:)
            transmission % tausun_total_path1(tb1:tb2)    = transmission1 % tausun_total_path1(:)
            transmission % tausun_total_path2(tb1:tb2)    = transmission1 % tausun_total_path2(:)
            radiancedata_d % clear(tb1:tb2)      = radiancedata_d1 % clear(:)
            radiancedata_d % total(tb1:tb2)      = radiancedata_d1 % total(:)
            radiancedata_d % bt_clear(tb1:tb2)   = radiancedata_d1 % bt_clear(:)
            radiancedata_d % bt(tb1:tb2)         = radiancedata_d1 % bt(:)
            radiancedata_d % refl_clear(tb1:tb2) = radiancedata_d1 % refl_clear(:)
            radiancedata_d % refl(tb1:tb2)       = radiancedata_d1 % refl(:)
            radiancedata_d % overcast(:,tb1:tb2) = radiancedata_d1 % overcast(:,:)
            radiancedata_d % cloudy(tb1:tb2)     = radiancedata_d1 % cloudy(:)
            
          end do

          ! transmittance deallocation for 1 profile
          deallocate(chanprof1)
          ! transmittance deallocation for 1 profile
          call rttov_alloc_transmission(allocStatus, transmission1, nlevels=nlv_T,  &
              nchanprof=tvs_nchan(sensorIndex), asw=0)
          if (allocStatus /= 0) call utl_abort('tvs_rttov: memory deallocation error in rttov_alloc_transmission')
          ! radiance deallocation for 1 profile
          call rttov_alloc_rad (allocStatus, tvs_nchan(sensorIndex), radiancedata_d1, nlv_T, asw=0)
          if (allocStatus /= 0) call utl_abort('tvs_rttov: memory deallocation error in rttov_alloc_rad')
        else

          ! run clear-sky RTTOV, save the radiances in OBS_BTCL of obsSpaceData 
          if (runObsOperatorWithClw  .and. obs_columnActive_RB(obsSpaceData, OBS_BTCL)) then
            
            ! set the cloud profile in tvs_profiles_nl to zero
            call updateCloudInTovsProfile(sensorHeaderIndexes,  &
                                           nlv_T,               &
                                           mode='save',         &
                                           beSilent=.true.)
            call rttov_parallel_direct(                         &
                rttov_err_stat,                                 & ! out
                tvs_chanProf(1:btCount,sensorIndex),            & ! in
                tvs_opts(sensorIndex),                          & ! in
                tvs_profiles_nl(sensorHeaderIndexes(:)),        & ! in
                tvs_coefs(sensorIndex),                         & ! in
                transmission,                                   & ! inout
                radiancedata_d,                                 & ! inout
                calcemis=calcemis,                              & ! in
                emissivity=emissivity_local,                    & ! inout
                nthreads=nthreads      )   

            ! save in obsSpaceData
            do btIndex = 1, btCount
              clearMwRadiance = radiancedata_d % bt(btIndex)
              bodyIndex = tvs_bodyIndexFromBtIndex(btIndex,sensorIndex)
              if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated) then
                call obs_bodySet_r(obsSpaceData, OBS_BTCL, bodyIndex, clearMwRadiance)
              end if
            end do

            ! tvs_profiles_nl
            call updateCloudInTovsProfile(sensorHeaderIndexes, &
                                          nlv_T,               &
                                          mode='restore',      &
                                          beSilent=.true.)
          end if

          if (.not. beSilent) write(*,*) 'before rttov_parallel_direct...', sensorIndex, profileCount, btCount
          call rttov_parallel_direct(                  &
              rttov_err_stat,                          & ! out
              tvs_chanProf(1:btCount,sensorIndex),     & ! in
              tvs_opts(sensorIndex),                   & ! in
              tvs_profiles_nl(sensorHeaderIndexes(:)), & ! in
              tvs_coefs(sensorIndex),                  & ! in
              transmission,                            & ! inout
              radiancedata_d,                          & ! inout
              calcemis=calcemis,                       & ! in
              emissivity=emissivity_local,             & ! inout
              nthreads=nthreads      )
          if (.not. beSilent) write(*,*) 'after rttov_parallel_direct...'
          if (rttov_err_stat /= 0) then
            write(*,*) 'Error in rttov_parallel_direct',rttov_err_stat
            call utl_abort('tvs_rttov')
          end if
        end if ! if (bgckMode .and. tvs_isInstrumHyperSpectral(instrum))
        
        !    2.4  Store hx in the structure tvs_radiance
        
        do btIndex = 1, btCount
          profileIndex = tvs_chanProf(btIndex,sensorIndex) % prof
          channelIndex = tvs_chanProf(btIndex,sensorIndex) % chan
          headerIndex = sensorHeaderIndexes(profileIndex)
          tvs_radiance(headerIndex) % bt(channelIndex) = radiancedata_d % bt(btIndex)

          if (bgckMode) then
            if (.not. associated(tvs_radiance(headerIndex) % clear)) then 
              allocate(tvs_radiance(headerIndex) % clear ( tvs_nchan(sensorIndex)))
              !  allocate overcast black cloud sky radiance output
              allocate(tvs_radiance(headerIndex) % overcast (nlv_T-1,tvs_nchan(sensorIndex)))
            end if
            tvs_radiance(headerIndex) % clear(channelIndex) =  &
                radiancedata_d % clear(btIndex)
            do levelIndex = 1, nlv_T - 1
              tvs_radiance(headerIndex) % overcast(levelIndex,channelIndex) =   &
                  radiancedata_d % overcast(levelIndex,btIndex)
            end do
          end if
          
          if (.not. allocated(tvs_transmission)) call tvs_allocTransmission(nlv_T)
   
          do levelIndex = 1, nlv_T
            tvs_transmission(headerIndex) % tau_levels(levelIndex,channelIndex) = &
                transmission % tau_levels(levelIndex,btIndex)
          end do
          
          tvs_transmission(headerIndex) % tau_total(channelIndex) = &
              transmission % tau_total(btIndex)

          if (allocated(tvs_emissivity)) then
            tvs_emissivity(channelIndex,headerIndex) = emissivity_local(btIndex) % emis_out
          end if
          
        end do

        !    Deallocate memory
        call rttov_alloc_direct(         &
            allocStatus,                 &
            asw=0,                       &
            nprofiles=profileCount,      & ! (not used)
            nchanprof=btCount,           &
            nlevels=nlv_T,               &
            opts=tvs_opts(sensorIndex),  &
            coefs=tvs_coefs(sensorIndex),&
            transmission=transmission,   &
            radiance=radiancedata_d,     &
            calcemis=calcemis,           &
            emissivity=emissivity_local, &
            init=.true.)
        if (allocStatus /= 0) call utl_abort('tvs_rttov: memory deallocation error 1 in rttov_alloc_direct')
        
        deallocate(surfem1)
      end if ! if (btCount > 0)
      
      if (btCountScatt > 0) then
         
        allocate(surfem1Scatt(btCountScatt))
        allocate(frequencies (btCountScatt))
        call rttov_alloc_direct(             &
            allocStatus,                     &
            asw=1,                           &
            nprofiles=profileCount,          & ! (not used)
            nchanprof=btCountScatt,          &
            nlevels=nlv_T,                   &
            opts=tvs_opts(sensorIndex),      &
            coefs=tvs_coefs(sensorIndex),    &
            radiance=radiancedata_dScatt,    &
            calcemis=calcemisScatt,          &
            emissivity=emissivity_localScatt,&
            init=.true.)
        if (allocStatus /= 0) call utl_abort('tvs_rttov: memory allocation error 2 in rttov_alloc_direct')
        
        call rttov_scatt_setupindex(                       &
            rttov_err_stat,                                &
            profileCount,                                  &  ! number of profiles
            tvs_nchan(sensorIndex),                        &  ! number of channels 
            tvs_coefs(sensorIndex),                        &  ! coef structure read in from rttov coef file
            tvs_coef_scatt(sensorIndex),                   &  ! 
            btCountScatt,                                  &  ! number of calculated channels
            tvs_chanProfScatt(1:btCountScatt,sensorIndex), &  ! channels and profile numbers
            frequencies,                                   &  ! array, frequency number for each channel
            lchannelSubset )                                  ! OPTIONAL array of logical flags to indicate a subset of channels
        deallocate(lchannelSubset)
        call tvs_getOtherEmissivities(tvs_chanProfScatt(1:btCountScatt,sensorIndex), sensorHeaderIndexes, sensorType, instrum, surfem1Scatt, calcemisScatt)
        
        call tvs_getMWemissivityFromAtlas(surfem1Scatt(1:btcountScatt), emissivity_localScatt, sensorIndex, tvs_chanProfScatt(1:btCountScatt,sensorIndex), &
                                          sensorHeaderIndexes)
      
        !   2.3  Compute radiance with rttov_direct

        rttov_err_stat = 0 
        ! run clear-sky RTTOV, save the radiances in OBS_BTCL of obsSpaceData
        if (runObsOperatorWithHydrometeors .and. obs_columnActive_RB(obsSpaceData, OBS_BTCL)) then
          ! run rttovScatt
          ! set the cloud profile in tvs_cld_profiles_nl to zero
          call updateCloudInTovsCloudProfile(sensorHeaderIndexes,      &
                                             nlv_T,                    &
                                             mode='save',              &
                                             beSilent=.true.)
          call tvs_rttovScatt(                                         &
              rttov_err_stat,                                          &! out
              tvs_opts_scatt(sensorIndex),                             &! in
              nlv_T,                                                   &! in
              tvs_chanProfScatt(1:btCountScatt,sensorIndex),           &! in
              frequencies,                                             &! in
              tvs_profiles_nl(sensorHeaderIndexes(:)),                 &! in
              tvs_cld_profiles_nl(sensorHeaderIndexes(:)),             &! in
              tvs_coefs(sensorIndex),                                  &! in
              tvs_coef_scatt(sensorIndex),                             &! in
              calcemisScatt,                                           &! in
              emissivity_localScatt,                                   &! inout
              transmission,                                            &
              radiancedata_dScatt) 

          ! save in obsSpaceData
          do btIndex = 1, btCountScatt
            clearMwRadiance = radiancedata_dScatt % bt(btIndex)
            bodyIndex = tvs_bodyIndexFromBtIndexScatt(btIndex,sensorIndex)
            profileIndex = tvs_chanprofScatt(btIndex,sensorIndex) % prof
            channelIndex = tvs_chanprofScatt(btIndex,sensorIndex) % chan
            headerIndex = sensorHeaderIndexes(profileIndex)
            if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated) then
              call obs_bodySet_r(obsSpaceData, OBS_BTCL, bodyIndex, clearMwRadiance)
            end if

            if (.not. allocated(tvs_transmission)) call tvs_allocTransmission(nlv_T)
   
            do levelIndex = 1, nlv_T
              tvs_transmission(headerIndex) % tau_levels(levelIndex,channelIndex) = &
                  transmission % tau_levels(levelIndex,btIndex)
            end do
            tvs_transmission(headerIndex) % tau_total(channelIndex) = &
                transmission % tau_total(btIndex)           
          end do
          
          ! restore the cloud profiles in ...
          call updateCloudInTovsCloudProfile(sensorHeaderIndexes,  &
                                             nlv_T,                &
                                             mode='restore',       &
                                             beSilent=.true.)
          
        end if ! run clear-sky RTTOV
        
        if (.not. beSilent) write(*,*) 'before rttov_scatt...', sensorIndex, profileCount, btCountScatt
        call tvs_rttovScatt(                                &
            rttov_err_stat,                                 &! out
            tvs_opts_scatt(sensorIndex),                    &! in
            nlv_T,                                          &! in
            tvs_chanprofScatt(1:btCountScatt,sensorIndex),  &! in
            frequencies,                                    &! in
            tvs_profiles_nl(sensorHeaderIndexes(:)),        &! in
            tvs_cld_profiles_nl(sensorHeaderIndexes(:)),    &! in
            tvs_coefs(sensorIndex),                         &! in
            tvs_coef_scatt(sensorIndex),                    &! in
            calcemisScatt,                                  &! in
            emissivity_localScatt,                          &! inout
            transmission,                                   &! allocated inside
            radiancedata_dScatt) 
        if (.not. beSilent) write(*,*) 'after rttov_scatt...'
        if (rttov_err_stat /= 0) then
          write(*,*) 'Error in rttov_scatt',rttov_err_stat
          call utl_abort('tvs_rttov')
        end if

        !    2.4  Store hx in the structure tvs_radiance        
        do btIndex = 1, btCountScatt
          profileIndex = tvs_chanProfScatt(btIndex,sensorIndex) % prof
          channelIndex = tvs_chanprofScatt(btIndex,sensorIndex) % chan
          headerIndex = sensorHeaderIndexes(profileIndex)
          tvs_radiance(headerIndex) % bt(channelIndex) = radiancedata_dScatt % bt(btIndex)

          if (bgckMode) then
            if (.not. associated(tvs_radiance(headerIndex) % clear)) then 
              allocate(tvs_radiance(headerIndex) % clear(tvs_nchan(sensorIndex)))
              !  allocate overcast black cloud sky radiance output
              allocate(tvs_radiance(headerIndex) % overcast(nlv_T-1,tvs_nchan(sensorIndex)))
            end if
            tvs_radiance(headerIndex) % clear(channelIndex) =  &
                radiancedata_dScatt % clear(btIndex)
            do levelIndex = 1, nlv_T - 1
              tvs_radiance(headerIndex) % overcast(levelIndex,channelIndex) = &
                  radiancedata_dScatt % overcast(levelIndex,btIndex)
            end do
            if (.not. allocated(tvs_transmission)) call tvs_allocTransmission(nlv_T)
          end if

          if (allocated(tvs_emissivity)) then
            tvs_emissivity(channelIndex,headerIndex) = emissivity_localScatt(btIndex) % emis_out
          end if

          if (allocated(tvs_transmission)) then
            do levelIndex = 1, nlv_T
              tvs_transmission(headerIndex) % tau_levels(levelIndex,channelIndex) = &
                  transmission % tau_levels(levelIndex,btIndex)
            end do
          
            tvs_transmission(headerIndex) % tau_total(channelIndex) = &
                transmission % tau_total(btIndex)
          end if
        end do

        !    Deallocate memory
        call rttov_alloc_direct(              &
             allocStatus,                     &
             asw=0,                           &
             nprofiles=profileCount,          & ! (not used)
             nchanprof=btCountScatt,          &
             nlevels=nlv_T,                   &
             opts=tvs_opts(sensorIndex),      &
             coefs=tvs_coefs(sensorIndex),    &
             radiance=radiancedata_dScatt,    &
             calcemis=calcemisScatt,          &
             emissivity=emissivity_localScatt,&
             init=.true.)
        if (allocStatus /= 0) call utl_abort('tvs_rttov: memory deallocation error 2 in rttov_alloc_direct')
        deallocate(surfem1Scatt)
        deallocate(frequencies)
        deallocate(transmission % tau_total)
        deallocate(transmission % tau_levels) 
      end if ! if (btCountScatt > 0)
      
      if (.not. beSilent) call msg_memUsage('tvs_rttov')
    end do sensor_loop
    
    if (allocated(tvs_bodyIndexFromBtIndex)) then
      deallocate( tvs_bodyIndexFromBtIndex )
      deallocate( tvs_bodyIndexFromBtIndexScatt )
      deallocate( tvs_chanProf )
      deallocate( tvs_chanProfScatt )
    end if
    write(*,*) 'tvs_rttov: Finished'

  end subroutine tvs_rttov

  !--------------------------------------------------------------------------
  !  tvs_getMWemissivityFromAtlas
  !--------------------------------------------------------------------------
  subroutine tvs_getMWemissivityFromAtlas(originalEmissivity, updatedEmissivity, sensorIndex, chanprof, sensorHeaderIndexes)
    implicit none

    ! Arguments:
    real(8),                 intent(in)  :: originalEmissivity(:) ! initial emissivity (typically 0.75 over land)
    type(rttov_emissivity),  intent(out) :: updatedEmissivity(:)  ! emissivity from atlas where possible
    integer,                 intent(in)  :: sensorIndex           ! sensor Index 
    type(rttov_chanprof),    intent(in)  :: chanprof(:)           ! chanprof RTTOV structure
    integer,                 intent(in)  :: sensorHeaderIndexes(:)! indexes of radiance observations for the currently processed sensor

    ! Locals:
    integer :: returnCode
    real(8) :: mWAtlasSurfaceEmissivity(size(originalEmissivity))
    integer :: btCount, profileCount
    integer :: profileIndex, btIndex, headerIndex
    
    btCount = size(originalEmissivity)
    if (useMWEmissivityAtlas) then

      if (.not. allocated (tvs_atlas)) allocate(tvs_atlas(tvs_nsensors))
      if (.not. tvs_atlas(sensorIndex) % init) then
        call rttov_setup_emis_atlas( returnCode, &! out
             tvs_opts(sensorIndex),              &! in
             tvs_profiles_nl(1) % date(2),       &! in
             atlas_type_mw,                      &! in
             tvs_atlas(sensorIndex),             &! inout
             atlas_id = mWAtlasId,               &! in ! 1 TELSEM2, 2 CNRM
             coefs = tvs_coefs(sensorIndex)  )
        if (returnCode /= 0) then
          write(*,*) 'Error in rttov_atlas_setup MW',returnCode
          call utl_abort('tvs_getMWemissivityFromAtlas')
        end if
      end if
   
      call rttov_get_emis( returnCode,             & ! out
           tvs_opts(sensorIndex),                  & ! in
           chanprof,                               & ! in
           tvs_profiles_nl(sensorHeaderIndexes(:)),& ! in
           tvs_coefs(sensorIndex),                 & ! in
           tvs_atlas(sensorIndex),                 & ! in
           mWAtlasSurfaceEmissivity)                ! out
    
      if (returnCode /= 0) then
        write(*,*) 'Error in rttov_get_emis MW', returnCode
        call utl_abort('tvs_getMWemissivityFromAtlas')
      end if

      profileCount = size( sensorHeaderIndexes )

      do profileIndex = 1, profileCount !loop on profiles
        headerIndex = sensorHeaderIndexes(profileIndex)
        do btIndex = 1, btCount !loop on channels
          if (chanprof(btIndex) % prof == profileIndex) then
            ! Now we have 0.75 in originalEmissivity(:) for land and sea ice
            ! and the MW atlas emissivity in mWAtlasSurfaceEmissivity(:)
            if (tvs_profiles_nl(headerIndex) % skin % surftype == surftype_land .and. &
                mWAtlasSurfaceEmissivity(btIndex) > 0.d0 .and. &
                mWAtlasSurfaceEmissivity(btIndex) <= 1.d0) then ! check for missing values
              updatedEmissivity(btIndex) % emis_in = mWAtlasSurfaceEmissivity(btIndex)
            else
              updatedEmissivity(btIndex) % emis_in = originalEmissivity(btIndex)
            end if
            ! Note that emissivity above sea-ice is not modified
          end if
        end do
      end do
    else
      updatedEmissivity(:) % emis_in = originalEmissivity(:)
    end if
    
  end subroutine tvs_getMWemissivityFromAtlas

  !--------------------------------------------------------------------------
  !  comp_ir_emiss
  !--------------------------------------------------------------------------
  subroutine comp_ir_emiss (emiss, wind, angle, nChannels, nProfiles, mchannel)
    !
    ! :Purpose: Computes water infrared emissivity for a specific set of
    !           channel indices, wind speed and zenith angle.
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: nChannels                  ! number of channels to process
    integer, intent(in)  :: nProfiles                  ! number of locations
    real(8), intent(out) :: emiss(nChannels,nProfiles) ! emissivities (0.-1.)
    real(8), intent(in)  :: wind(nProfiles)            ! surface wind speed (m/s)
    real(8), intent(in)  :: angle(nProfiles)           ! viewing angle (deg)
    integer, intent(in)  :: mchannel(nChannels)        ! vector of channel indices to process

    ! Locals:
    integer, parameter :: MaxWn = 19
    integer, parameter :: Nparm = 3
    integer, parameter :: MaxChan = 19

    real(8), parameter :: Theta(Nparm,MaxWn) = [  &
         1700.381d0, 25.28534d0, 144.1023d0,      &
         1738.149d0, 25.67787d0, 146.6139d0,      & 
         1769.553d0, 26.05250d0, 148.6586d0,      &
         1778.610d0, 26.12333d0, 149.5127d0,      &
         1794.245d0, 26.18523d0, 150.5874d0,      &
         1791.904d0, 26.19991d0, 150.7076d0,      &
         1806.872d0, 26.37132d0, 151.7191d0,      &
         1926.078d0, 27.63825d0, 160.7103d0,      &
         1969.155d0, 28.02767d0, 163.6069d0,      &
         1975.549d0, 27.86465d0, 164.6228d0,      &
         1991.288d0, 27.94312d0, 166.2924d0,      &
         2082.691d0, 28.93558d0, 172.4025d0,      &
         2182.872d0, 29.89974d0, 179.5839d0,      &
         2338.510d0, 31.27507d0, 191.2063d0,      &
         2164.615d0, 28.97152d0, 182.6279d0,      &
         2099.714d0, 29.91868d0, 178.4015d0,      &
         1857.644d0, 29.13640d0, 160.9822d0,      &
         1610.696d0, 26.48602d0, 142.2768d0,      &
         1503.969d0, 24.97931d0, 133.4392d0 ]

    real(8), parameter ::  C(Nparm,2,MaxWn) = [                                  &  
         0.9715104043561414d0,-1.2034233230944147D-06, -5.8742655960993913D-07,  &
         0.9263932848727608d0,-9.4908630939690859D-04, 2.2831134823358876D-05,   &
         0.9732503924722753d0,-1.2007007329295099D-06, -5.8767355551283423D-07,  &
         0.9290947860585505d0,-9.5233413988900161D-04, 2.2640835623043761D-05,   &
         0.9745005204317289d0, 1.2857517639804244D-06, -7.1356127087301190D-07,  &
         0.9310852809117095d0,-9.5453509182819095D-04, 2.2562638663187251D-05,   &
         0.9756204829761132d0, 1.2979181109743976D-06, -7.1406809362820139D-07,  &
         0.9329073568177888d0,-9.5627536945214183D-04, 2.2442358508999558D-05,   &
         0.9764012672766408d0,-2.0826654381361387D-06, -4.9103920569405721D-07,  &
         0.9341937281933334d0,-9.5764423928102976D-04, 2.2326701573603621D-05,   &
         0.9770513558720460d0, 4.1867599593267133D-07, -6.1768073971231931D-07,  &
         0.9352981872014672d0,-9.5833614545300181D-04, 2.2261996883974513D-05,   &
         0.9775970810179080d0,-1.2289690625562906D-06, -5.2953762169985775D-07,  &
         0.9362188153954743d0,-9.5950872922696905D-04, 2.2251301675423482D-05,   &
         0.9830610391451819d0, 2.7693589475690676D-07, -5.1580217018207558D-07,  &
         0.9461121192685766d0,-9.5718115604053031D-04, 2.1087308573177295D-05,   &
         0.9840097930773377d0,-1.4987900189155091D-06, -3.8281408128977038D-07,  &
         0.9479758694804105d0,-9.5451252460440695D-04, 2.0800627740862229D-05,   &
         0.9851056150728580d0,-6.5768237152417477D-07, -4.2053769829400935D-07,  &
         0.9502084544618980d0,-9.4965534997704157D-04,  2.0326602209199427D-05,  &
         0.9862706396188835d0,-2.3713068057993353D-06, -2.8671134918457728D-07,  &
         0.9526580467595886d0,-9.4614505430749598D-04,  2.0001856872595840D-05,  &
         0.9875307221489201d0, 1.3003462826947714D-07, -4.1335288320283954D-07,  &
         0.9554195617948236d0,-9.3806678196435643D-04,  1.9407754748128057D-05,  &
         0.9891153260567763d0,-8.0730206675976713D-07, -3.1811320626834656D-07,  &
         0.9590558393678170d0,-9.2716287670223167D-04, 1.8690586764925213D-05,   &
         0.9913304557147524d0,-2.1153391229093421D-08, -3.1094269595901165D-07,  &
         0.9644162604969492d0,-9.0342753739935612D-04, 1.7274993357160937D-05,   &
         0.9925188366950193d0,-4.6365959315123653D-07, -2.7020120347068712D-07,  &
         0.9667877170960085d0,-9.0665804037922043D-04, 1.7083616616646458D-05,   &
         0.9919408379810360d0,-2.0563508815953840D-06, -1.8066722718403761D-07,  &
         0.9627535343397309d0,-9.7537134133678965D-04,  1.9698263973541952D-05,  &
         0.9889406296815972d0,-2.3713068057993353D-06, -2.8671134918457728D-07,  &
         0.9506051906192242d0,-1.0642902225813857D-03,  2.4235485973033298D-05,  &
         0.9828819693848310d0,-7.4086701870172759D-07, -6.2949258820534062D-07,  &
         0.9329616683158125d0,-1.0728027288012200D-03, 2.7209071863380586D-05,   &
         0.9767410313266288d0,-9.1750038410238915D-07, -7.9177921107781349D-07,  &
         0.9192775350344998d0,-1.0369254272157462D-03, 2.8000594542037504D-05 ]

    real(8) :: a(MaxChan), b(MaxChan), cc(MaxChan)
    real(8) :: ww
    integer :: bandIndex, channelIndex, profileIndex

    do channelIndex = 1 , nChannels

      bandIndex = Mchannel(channelIndex)

      do profileIndex = 1, nProfiles

        ww = wind(profileIndex)
        a(channelIndex) = c(1,1,bandIndex) + c(2,1,bandIndex) * ww    &  
             + c(3,1,bandIndex) * ww * ww
        b(channelIndex) = c(1,2,bandIndex) + c(2,2,bandIndex) * ww    &
             + c(3,2,bandIndex)* ww * ww

        cc(channelIndex) = Theta(1,bandIndex) + Theta(2,bandIndex) * ww

        emiss(channelIndex,profileIndex) = a(channelIndex) + (b(channelIndex) - a(channelIndex)) *   & 
             exp(( (theta(3,bandIndex) - 60.d0)**2.d0              &
             - (angle(profileIndex) - theta(3,bandIndex))**2.d0 ) / CC(channelIndex))
       
      end do
      
    end do

  end subroutine comp_ir_emiss

  !--------------------------------------------------------------------------
  !  tvs_pcnt_box
  !--------------------------------------------------------------------------
  subroutine tvs_pcnt_box(f_low, nprf, ilat, ilon, ireduc)
    !
    ! :Purpose: Computes a low resolution feature form a high
    !           resolution one by averaging.
    !           example: use for percentage of water
    !
    implicit none

    ! Arguments:
    integer, intent(in)   :: nprf              ! Number of profiles
    real(8), intent(out)  :: f_low(nprf)       ! Low resolution field 
    integer, intent(in)   :: ilat(nprf)        ! Y-coordinate of profile
    integer, intent(in)   :: ilon(nprf)        ! X-coordinate of profile
    integer, intent(in)   :: ireduc            ! Means a 2xireduc+1 by 2xireduc+1 averaging

    ! Locals
    integer :: nplon, jdlo1, jdlo2, jlon1, jlon2
    integer :: nx, ilat1, ilat2, ilon1, ilon2, profileIndex, lonIndex, latIndex
   
    profiles : do profileIndex = 1, nprf

      nplon = 0

      ! normal limits

      ilat1=max(ilat(profileIndex)-ireduc,1)
      ilat2=min(ilat(profileIndex)+ireduc,kslat)
      ilon1=max(ilon(profileIndex)-ireduc,1)
      ilon2=min(ilon(profileIndex)+ireduc,kslon)

      if (ilon1 == 1 .or. ilon2 == kslon) then
        ! border cases for longitudes
        jdlo1 = ilon(profileIndex)-ireduc
        jdlo2 = ilon(profileIndex)+ireduc

        if (jdlo1 <= 0) then
          nplon = 1
          jlon1 = kslon + jdlo1
          jlon2 = kslon
        else if (jdlo2 > kslon) then
          nplon = 1
          jlon1 = 1
          jlon2 = jdlo2 - kslon
        end if
      end if

      nx = 0
      f_low(profileIndex) = 0.d0
     
      do latIndex = ilat1, ilat2

        do lonIndex = ilon1, ilon2
          nx = nx + 1
          f_low(profileIndex) = f_low(profileIndex) + waterFraction(lonIndex,latIndex)         
        end do
        
        if (nplon == 1) then
          ! additional cases at border 1-kslon
          do lonIndex = jlon1, jlon2
            nx = nx + 1
            f_low(profileIndex) = f_low(profileIndex) + waterFraction(lonIndex,latIndex)         
          end do
        end if

      end do
      
      f_low(profileIndex) = f_low(profileIndex) / dble(nx)

    end do profiles

  end subroutine tvs_pcnt_box

  !--------------------------------------------------------------------------
  !  tvs_emis_read_climatology
  !--------------------------------------------------------------------------
  subroutine tvs_emis_read_climatology
    !
    ! :Purpose: Read information about ceres surface type and water fraction.
    !
    implicit none
    
    ! Locals:
    integer            :: nisf,njsf,nksf
    integer            :: niwa,njwa,nkwa
    character(len=20)  :: ceresFile
    integer, external  :: fnom,fstouv,fstfrm,fclos,fstlir
    integer            :: isftest
    integer            :: iv1,iv2,iv3,iv4,iv5,iv6

    isftest = 0

    !* get surface type and water fraction
    ceresFile = 'ceres_global.std'
    iv1 = fnom(isftest,ceresFile,'RND+R/O',0)
    iv2 = fstouv(isftest,'RND')
    iv3 = fstlir(surfaceType,isftest,nisf,njsf,nksf,-1,'SFC-TYPE',-1,-1,-1,'','TY')
    iv4 = utl_fstlir(waterFraction,isftest,niwa,njwa,nkwa,-1,'WATER_FR',-1,-1,-1,'','W%')
    iv5 = fstfrm(isftest)
    iv6 = fclos(isftest)

    if (iv1 < 0 .or. iv2 < 0 .or. iv3 < 0 .or. iv4 < 0 .or. iv5 < 0 .or. iv6 < 0) then
      write(*,*) 'LES iv DE CERES ',iv1,iv2,iv3,iv4,iv5,iv6
      write(*,*) 'THESE NUMBER SHOULD NOT BE NEGATIVE WHEN DOING AIRS BACKGROUND CHECK'
      call utl_abort('Problem with file ceres_global.std in tvs_emis_read_climatology ')
    end if
   
  end subroutine tvs_emis_read_climatology

  !--------------------------------------------------------------------------
  !  emis_getIrEmissivity
  !--------------------------------------------------------------------------
  subroutine emis_getIrEmissivity (surfem1, nchn, sensorIndex, nprf, btCount, sensorHeaderIndexes)
    !
    ! :Purpose: Assign new ir surface emissivities based on
    !           cmc analysis surface albedo, sea ice fraction and snow mask
    !           in addition to ceres surface type and water fraction. 
    !           This is a subroutine that can apply to any IR instrument.
    !
    implicit none
   
    ! Arguments:
    integer, intent(in)  :: nprf                     ! Number of profiles
    integer, intent(in)  :: btCount                  ! Total number of observations treated
    real(8), intent(out) :: surfem1(btCount)         ! IR surface emissivity estimate (0-1)
    integer, intent(in)  :: nchn                     ! Number of channels
    integer, intent(in)  :: sensorindex              ! Sensor number
    integer, intent(in)  :: sensorHeaderIndexes(nprf)! header indexes of radiance observations for the currently processed sensor

    ! Locals:
    integer :: channelIndex, profileIndex, headerIndex, errorStatus
    integer :: btIndex
    integer :: ilat(nprf), ilon(nprf), surfType(nprf)
    real(8) :: latitudes(nprf), longitudes(nprf), satzang(nprf)
    real(8) :: wind_sfc(nprf), f_low(nprf), waven(nchn), em_oc(nchn,nprf), emi_mat(nchn,20)
    real(8) :: emissivity(btCount)
    logical :: thermal(btCount), calcemis(btCount)
    type(rttov_geometry) :: geometry(nprf) 
    real(8) :: uOfWLandWSurfaceEmissivity(btcount)
    
    ! Information to extract (transvidage)
    ! latitudes(nprf) -- latitude (-90 to 90)
    ! longitudes(nprf) -- longitude (0 to 360)
    ! satzang(nprf) -- satellite zenith angle (deg)

    do profileIndex = 1, nprf
      headerIndex = sensorHeaderIndexes(profileIndex)
      latitudes(profileIndex)  = tvs_profiles_nl(headerIndex) % latitude
      longitudes(profileIndex) = tvs_profiles_nl(headerIndex) % longitude
      satzang(profileIndex)    = tvs_profiles_nl(headerIndex) % zenangle
    end do

    !  Assign surface properties from grid to profiles
    call tvs_interp_sfc(ilat,ilon, nprf,latitudes,longitudes,sensorHeaderIndexes)

    !  Find the sensor bands (central) wavenumbers
    do channelIndex = 1, nchn      
      waven(channelIndex) = tvs_coefs(sensorIndex) % coef % ff_cwn(channelIndex)
    end do 

    if (tvs_oldFashionIRLandEmiss) then

      !  Get the CERES emissivity matrix for all sensor wavenumbers and surface types
      call ceres_ematrix(emi_mat, waven, nchn)
      
    else
      
      !  Get the Camel V2 emissivity climatology
      if (.not. allocated (tvs_atlas)) allocate(tvs_atlas(tvs_nsensors))
      if (.not. tvs_atlas(sensorIndex) % init) then
        call rttov_setup_emis_atlas(errorStatus, & ! out
            tvs_opts(sensorIndex),               & ! in
            tvs_profiles_nl(1) % date(2),        & ! in
            atlas_type_ir,                       & ! in
            tvs_atlas(sensorIndex),              & ! inout
            atlas_id = camel_clim_atlas_id,      & ! in
            coefs = tvs_coefs(sensorIndex),      & ! in
            ir_atlas_read_std = .false.,         & ! in
            ir_atlas_ang_corr = tvs_irEmissAngularCorrection ) ! in
            
        if (errorStatus /= 0) then
          write(*,*) 'Error in rttov_atlas_setup IR', errorStatus
          call utl_abort('emis_getIrEmissivity')
        end if
      end if
      
    end if

    if (tvs_oldFashionIRSeaEmiss) then
      do profileIndex = 1, nprf
        !       find surface wind
        headerIndex = sensorHeaderIndexes(profileIndex)
        wind_sfc(profileIndex) = min(sqrt(tvs_profiles_nl(headerIndex) % s2m %u**2 + tvs_profiles_nl(headerIndex) % s2m % v**2 + 1.d-12),15.d0)
      end do
      ! find ocean emissivities  
      call emi_sea (em_oc, waven,satzang,wind_sfc,nprf,nchn)
    else
      thermal(:) = .true.
      calcemis(:) = .true.
      do profileIndex = 1, nprf
        geometry(profileIndex) % normzen = tvs_profiles_nl(sensorHeaderIndexes(profileIndex))%zenangle / 60.0_jprb
        ! Save RTTOV surface type
        headerIndex = sensorHeaderIndexes(profileIndex)
        surfType(profileIndex) = tvs_profiles_nl(headerIndex) % skin % surftype
        tvs_profiles_nl(headerIndex) % skin % surftype = surftype_sea
      end do
      call rttov_calcemis_ir(                      &
          errorStatus,                             &
          tvs_opts(sensorIndex),                   &
          tvs_chanProf(1:btCount,sensorIndex),     &
          tvs_profiles_nl(sensorHeaderIndexes),    &
          geometry,                                &
          tvs_coefs(sensorIndex),                  &
          thermal,                                 &
          calcemis,                                &
          emissivity)
      if (errorStatus /= 0) then
        call utl_abort('emis_getIrEmissivity: problem in rttov_calcemis_ir')
      end if
      !Restore RTTOV surface type
      do profileIndex = 1, nprf
        headerIndex = sensorHeaderIndexes(profileIndex)
        tvs_profiles_nl(headerIndex) % skin % surftype = surfType(profileIndex)
      end do
      do btIndex = 1, btCount
        profileIndex = tvs_chanProf(btIndex,sensorIndex)%prof
        channelIndex = tvs_chanProf(btIndex,sensorIndex)%chan
        em_oc(channelIndex,profileIndex) = emissivity(btIndex)
      end do
    end if

    ! Get surface emissivities
    do profileIndex = 1, nprf
      headerIndex = sensorHeaderIndexes(profileIndex)
      !       set albedo to 0.6 where snow is present
      if (tvs_profiles_nl(headerIndex) % skin % surftype == surftype_land .and. tvs_surfaceParameters(headerIndex) % snow > 0.999) tvs_surfaceParameters(headerIndex) % albedo = 0.6
      !       if albedo too high no water
      if (tvs_surfaceParameters(headerIndex) % albedo >= 0.55) tvs_surfaceParameters(headerIndex) % pcnt_wat = 0.
     
      if (tvs_oldFashionIRLandEmiss) then
        !       if water and CMC ice present then sea ice
        if (tvs_profiles_nl(headerIndex) % skin % surftype == surftype_sea .and. tvs_surfaceParameters(headerIndex) % ice > 0.001) tvs_surfaceParameters(headerIndex) % ltype = 20
        !       if land and CMC snow present then snow
        if (tvs_profiles_nl(headerIndex) % skin % surftype == surftype_land .and. tvs_surfaceParameters(headerIndex) % snow > 0.999) tvs_surfaceParameters(headerIndex) % ltype = 15
        do channelIndex = 1, nchn
          surfem1((profileIndex-1)*nchn+channelIndex) =  tvs_surfaceParameters(headerIndex) % pcnt_wat * em_oc(channelIndex,profileIndex)  +  &
              ( 1.d0 - tvs_surfaceParameters(headerIndex) % pcnt_wat ) * emi_mat(channelIndex,tvs_surfaceParameters(headerIndex) % ltype)
        end do
      else

        !       if water and CMC ice present then sea ice
        if (tvs_profiles_nl(headerIndex) % skin % surftype == surftype_sea .and. tvs_surfaceParameters(headerIndex) % ice > 0.001) tvs_profiles_nl(headerIndex) % skin % surftype = surftype_seaice
        !       if land and CMC snow present then snow
        if (tvs_profiles_nl(headerIndex) % skin % surftype == surftype_land .and. tvs_surfaceParameters(headerIndex) % snow > 0.999) tvs_profiles_nl(headerIndex) % skin % snow_fraction = 1.d0
        call rttov_get_emis(errorStatus,          & ! out
            tvs_opts(sensorIndex),                & ! in
            tvs_chanprof(1:btCount,sensorIndex),  & ! in
            tvs_profiles_nl(sensorHeaderIndexes), & ! in
            tvs_coefs(sensorIndex),               & ! in
            tvs_atlas(sensorIndex),               & ! inout
            uOfWLandWSurfaceEmissivity )            ! out
        if (errorStatus /= 0) then
          write(*,*) 'Error in rttov_get_emis IR', errorStatus
          call utl_abort('emis_getIrEmissivity')
        end if
        if (uOfWLandWSurfaceEmissivity((profileIndex-1)*nchn + 1) < 0.d0) then
          ! land Emissivity atlas is not defined where there is too much water
          do channelIndex = 1, nchn
            surfem1((profileIndex-1)*nchn+channelIndex) = em_oc(channelIndex,profileIndex)
          end do
        else
          if (tvs_useWaterFraction) then
            do channelIndex = 1, nchn
              surfem1((profileIndex-1)*nchn+channelIndex) = tvs_surfaceParameters(headerIndex) % pcnt_wat * em_oc(channelIndex,profileIndex) +  &
                  ( 1.d0 - tvs_surfaceParameters(headerIndex) % pcnt_wat ) * uOfWLandWSurfaceEmissivity((profileIndex-1)*nchn+channelIndex)
            end do
          else
            do channelIndex = 1, nchn
              surfem1((profileIndex-1)*nchn+channelIndex) = uOfWLandWSurfaceEmissivity((profileIndex-1)*nchn+channelIndex)
            end do
          end if
        end if
        
      end if
      
    end do
    
    ! Find the regional water fraction (here in a 15x15 pixel box centered on profile)
    call tvs_pcnt_box (f_low, nprf, ilat, ilon, 7)

    do profileIndex = 1, nprf
      tvs_surfaceParameters(sensorHeaderIndexes(profileIndex)) % pcnt_reg = f_low(profileIndex)
    end do

  end subroutine emis_getIrEmissivity

  !--------------------------------------------------------------------------
  !  tvs_interp_sfc
  !--------------------------------------------------------------------------
  subroutine tvs_interp_sfc (ilat, ilon, nprf, latitudes, longitudes, sensorHeaderIndexes,skipAlbedo_opt)
    !
    ! :Purpose: Associate surface albedo, ice fraction, snow depth 
    !           and ceres surface type and water fraction to observations profiles.
    !
    implicit none

    ! Arguments:
    integer,           intent(in)  :: nprf                      ! number of profiles
    integer,           intent(out) :: ilat(nprf)                ! y-coordinate of profile
    integer,           intent(out) :: ilon(nprf)                ! x-coordinate of profile 
    real(8),           intent(in)  :: latitudes(nprf)           ! latitude (-90s to 90n)
    real(8),           intent(in)  :: longitudes(nprf)          ! longitude (0 to 360)
    integer,           intent(in)  :: sensorHeaderIndexes(nprf) ! header indexes of radiance observations for the currently processed sensor
    logical,optional,  intent(in)  :: skipAlbedo_opt            ! skip the section of code about Albedo  

    ! Locals:
    character(len=20)  :: cfile3,cfile5
    integer            :: iun3,iun5
    integer            ::                        iv7
    integer            :: ix1,ix2,ix3,ix4,ix5,        ix8,ix9,ix10,ix11,ix12
    integer            ::         iy3,iy4,iy5,        iy8,iy9,iy10
    integer            :: iz1,iz2,iz3,iz4,iz5,        iz8,iz9,iz10,iz11,iz12
    integer            :: ni3,nj3,nk3
    integer            :: ni4,nj4,nk4
    integer            :: ni5,nj5,nk5
    integer            :: dateo,deet,npas,nbits,datyp
    integer            :: ip1,ip2,ip3
    integer            :: ig13,ig23,ig33,ig43
    integer            :: ig14,ig24,ig34,ig44
    integer            :: ig15,ig25,ig35,ig45
    integer            :: swa,lng,dltf,ubc,ex1,ex2,ex3
    integer            :: profileIndex, headerIndex
    character(len=1)   :: typvar
    character(len=1)   :: grtyp3,grtyp4,grtyp5
    character(len=2)   :: nomvar, snowvar
    character(len=8)   :: etiket
    integer, external  :: fnom,fstouv,fstinf,fstprm,fstfrm,fclos
    integer, external  :: ezqkdef,ezdefset
    real(8)            :: zig1,zig2,zig3,zig4
    integer            :: ig1obs,ig2obs,ig3obs,ig4obs
    real (8)           :: alat, alon, zzlat, zzlon
    logical            :: skipAlbedo
    ! fields on input grid
    real(8), allocatable :: glace(:,:), neige(:,:), alb(:,:)
    ! fields on output grid
    real(8)              :: glace_intrpl(nprf,1), neige_intrpl(nprf,1), alb_intrpl(nprf,1)

    ! printout header
    write(*,*) 
    write(*,*) 'SUBROUTINE tvs_interp_sfc'
    write(*,*) '---------------------'
    write(*,*) ' called multiple time by bunch of ',nprf,' profiles'
    write(*,*) ' <RETURN CODES> SHOULD NOT BE NEGATIVE'
    write(*,*) '---------------------------------------------------'

    if (present(skipAlbedo_opt)) then
      skipAlbedo = skipAlbedo_opt
    else
      skipAlbedo = .false.
    end if

    ! --- FOR CERES VARIABLES -------------
    !  Get number of pixels per degree of lat or lon

    alat = dble(kslat)/180.d0
    alon = dble(kslon)/360.d0

    do profileIndex=1, nprf

      ! get lat and lon within limits if necessary
      zzlat = min(latitudes(profileIndex),89.999d0)
      zzlat = max(Zzlat,-89.999d0)
      
      zzlon = min(longitudes(profileIndex),359.999d0)
      zzlon = max(zzlon,0.d0)

      !  Find in which surface field pixel is located the observation profile

      ! Note : CERES grid at 1/6 resolution 
      !         N-S : starts at N pole and excludes S pole
      !         W-E : starts at longitude 0 and excludes longitude 360

      ilat(profileIndex) = max( nint((zzlat + 90.d0) * alat),1) 
      ilon(profileIndex) = nint(zzlon * alon) + 1
      if (ilon(profileIndex) > kslon) ilon(profileIndex) = 1

      !  Assign surface caracteristics to observation profiles
      headerIndex = sensorHeaderIndexes(profileIndex)
      tvs_surfaceParameters(headerIndex) % ltype    = surfaceType(ilon(profileIndex),ilat(profileIndex))
      tvs_surfaceParameters(headerIndex) % pcnt_wat = waterFraction(ilon(profileIndex),ilat(profileIndex))

    end do

    if (skipAlbedo) return
    
    !  For ice, snow and albedo variables -------------

    iun3 = 0
    iun5 = 0

    ! File names
    cfile3 = 'sfc4airs'          ! for ice fraction and snow cover
    cfile5 = 'sfc4airs_newalb'   ! for albedo

    ! fnom: make the connections with the external files name
    ! success = 0
    write(*,*) 
    ix1 = fnom(iun3,cfile3,'RND+R/O',0)
    write(*,*) 'file = sfc4airs         : fnom   : return = ', ix1

    iz1 = fnom(iun5,cfile5,'RND+R/O',0)
    write(*,*) 'file = sfc4airs_newalb  : fnom   : return = ', iz1

    ! fstouv: open the standard files
    ! success = number of records found in the file
    write(*,*) 
    ix2 = fstouv(iun3,'RND')
    write(*,*) 'file = sfc4airs         : fstouv : return = ', ix2
    iz2 = fstouv(iun5,'RND')
    write(*,*) 'file = sfc4airs_newalb  : fstouv : return = ', iz2

    ! fstinf: locate the records that matches the search keys
    ! success = handle of the record found after the search
    ! desired output = handle
    write(*,*) 
    ix3 = fstinf(iun3,ni3,nj3,nk3,-1,'',-1,-1,-1,'','LG')
    write(*,*) 'variable = LG           : fstinf : return = ', ix3

    snowvar='SD'
    iy3 = fstinf(iun3,ni4,nj4,nk4,-1,'',-1,-1,-1,'',snowvar)
    write(*,*) 'variable = ', snowvar, '           : fstinf : return = ', iy3
    if (iy3  <  0) then
      write(*,*) 'did not find ''SD'' so look for ''NE'''
      snowvar='NE'
      iy3 = fstinf(iun3,ni4,nj4,nk4,-1,'',-1,-1,-1,'',snowvar)
      write(*,*) 'variable = ', snowvar, '           : fstinf : return = ', iy3
    end if

    iz3 = fstinf(iun5,ni5,nj5,nk5,-1,'',-1,-1,-1,'','AL')
    write(*,*) 'variable = AL           : fstinf : return = ', iz3

    ! fstprm: get the description informations of the record given the key
    ! success = 0
    ! desired output = nix,njx,grtypx,igxx,ig1x,ig2x,ig3x,ig4x

    write(*,*) 
    ix4 = fstprm(ix3, dateo,deet,npas,ni3,nj3,nk3,nbits,datyp, &
         ip1,ip2,ip3,typvar,nomvar,etiket,grtyp3,  &
         ig13,ig23,ig33,ig43,swa,lng,dltf,ubc,ex1,ex2,ex3)
    write(*,*) 'variable = LG           : fstprm : return = ', ix4

    iy4 = fstprm(iy3, dateo,deet,npas,ni4,nj4,nk4,nbits,datyp, &
         ip1,ip2,ip3,typvar,nomvar,etiket,grtyp4,  &
         ig14,ig24,ig34,ig44,swa,lng,dltf,ubc,ex1,ex2,ex3)
    write(*,*) 'variable = ', snowvar, '           : fstprm : return = ', iy4

    iz4 = fstprm(iz3, dateo,deet,npas,ni5,nj5,nk5,nbits,datyp, &
         ip1,ip2,ip3,typvar,nomvar,etiket,grtyp5,  &
         ig15,ig25,ig35,ig45,swa,lng,dltf,ubc,ex1,ex2,ex3)
    write(*,*) 'variable = AL           : fstprm : return = ', iz4

    ! allocation of the field on the grid
    allocate (glace(ni3,nj3))
    allocate (neige(ni4,nj4))
    allocate (alb  (ni5,nj5))

    ! utl_fstlir: read records data (field on the grid) given the key
    ! success = handle of the record
    ! desired output = FIELD
    write(*,*) 

    ix5 = utl_fstlir(glace, iun3,ni3,nj3,nk3,-1,'',-1,-1,-1,'','LG')
    write(*,*) 'variable = LG           : utl_fstlir : return = ', ix5
    iy5 = utl_fstlir(neige, iun3,ni4,nj4,nk4,-1,'',-1,-1,-1,'',snowvar)
    write(*,*) 'variable = ', snowvar, '           : utl_fstlir : return = ', iy5
    iz5 = utl_fstlir(alb,   iun5,ni5,nj5,nk5,-1,'',-1,-1,-1,'','AL')
    write(*,*) 'variable = AL           : utl_fstlir : return = ', iz5

    ! int_CXGAIG: define the grid descriptors (integer form) of the
    !          observation profile output grid
    ! desired output = ig1OBS, ig2OBS, ig3OBS, ig4OBS
    zig1 = 0.0d0
    zig2 = 0.0d0
    zig3 = 1.0d0
    zig4 = 1.0d0

    call int_cxgaig('L',ig1OBS,ig2OBS,ig3OBS,ig4OBS,zig1,zig2,zig3,zig4)

    ! int_EZGDEF: define the grid of the observations profiles (output grid)
    ! of type Y containing the lat-lon of profiles
    ! success = token to identify the grid
    ! desired output = token
    write(*,*) 
    iv7 = int_ezgdef(nprf,1,'Y','L',ig1obs,ig2obs,ig3obs,ig4obs,longitudes,latitudes)
    write(*,*) 'apply to all variables  : int_EZGDEF : return = ', iv7
    
    ! EZQKDEF: define the grid of the records data (input grid)
    ! success = token to identify the grid
    ! desired output = token
    ! EZDEFSET: interpolate from input grids to output grid
    ! success = key
    ! int_hInterpScalar: interpolation of the field on the input grid to observation profiles
    ! success = 0
    ! desired output = FIELD_intrpl
    write(*,*) 
    ix8 = ezqkdef(ni3,nj3,grtyp3,ig13,ig23,ig33,ig43,iun3)
    write(*,*) 'variable = LG           : ezqkdef  : return = ', ix8
    
    ix9 = ezdefset(iv7,ix8)
    write(*,*) 'variable = LG           : ezdefset : return = ', ix9

    ix10 = int_hInterpScalar(glace_intrpl,glace,interpDegree='NEAREST')
    write(*,*) 'variable = LG           : int_hInterpScalar  : return = ', ix10

    write(*,*) 

    iy8 = ezqkdef(ni4,nj4,grtyp4,ig14,ig24,ig34,ig44,iun3)
    write(*,*) 'variable = ', snowvar, '           : ezqkdef  : return = ', iy8

    iy9 = ezdefset(iv7,iy8)
    write(*,*) 'variable = ', snowvar, '           : ezdefset : return = ', iy9

    iy10 = int_hInterpScalar(neige_intrpl,neige,interpDegree='NEAREST')
    write(*,*) 'variable = ', snowvar, '           : int_hInterpScalar  : return = ', iy10

    write(*,*) 

    iz8 = ezqkdef(ni5,nj5,grtyp5,ig15,ig25,ig35,ig45,iun5)
    write(*,*) 'variable = AL           : ezqkdef  : return = ', iz8

    iz9 = ezdefset(iv7,iz8)
    write(*,*) 'variable = AL           : ezdefset : return = ', iz9

    iz10 = int_hInterpScalar(alb_intrpl,alb,interpDegree='NEAREST')
    write(*,*) 'variable = AL           : int_hInterpScalar  : return = ', iz10

    ! fstfrm: close the standard files
    ! success = 0
    write(*,*) 
    ix11 = fstfrm(iun3)
    write(*,*) 'file = sfc4airs         : fstfrm : return = ', ix11
    
    iz11 = fstfrm(iun5)
    write(*,*) 'file = sfc4airs_newalb  : fstfrm : return = ', iz11
 
    ! fclos: release the connections with the external files name
    ! success = 0

    write(*,*) 

    ix12 = fclos(iun3)
    write(*,*) 'file = sfc4airs         : fclos  : return = ', ix12

    iz12 = fclos(iun5)
    write(*,*) 'file = sfc4airs_newalb  : fclos  : return = ', iz12

    ! assign surface caracteristics to observation profiles

    do profileIndex = 1, nprf
      headerIndex = sensorHeaderIndexes(profileIndex)
      tvs_surfaceParameters(headerIndex) % ice    = glace_intrpl(profileIndex,1)
      tvs_surfaceParameters(headerIndex) % snow   = neige_intrpl(profileIndex,1)
      tvs_surfaceParameters(headerIndex) % albedo = alb_intrpl(profileIndex,1)
    end do

    deallocate(glace,neige,alb)

  end subroutine tvs_interp_sfc

  !--------------------------------------------------------------------------
  !  ceres_ematrix
  !--------------------------------------------------------------------------
  subroutine ceres_ematrix(emi_mat, waven, nchn)
    !
    ! :Purpose: Set up emissivity versus fixed wavenumbers and surface types.
    !
    ! :CERES:
    ! Emissivity data available at low spectral resolution: only 14 values 
    ! to cover the entire spectrum. Thus, this can be used as a nominal value.
    ! The error associated with this emissivity can roughly be estimated to
    ! increase with lower emissivity as : (1-EMI)*0.5 
    ! (i.e. as large as 0.10 for EMI=0.80 but better than 0.01 for EMI > 0.98)
    ! -No dependence on viewing angle is assumed.
    ! -Not to be used for oceans uncovered by ice.
    !
    ! :Longwave Emmissivities in 12 original Fu bands + 2 extra to cover the range:
    !
    ! Longwave spectral intervals [cm-1] for the Fu & Liou code.
    !
    ! ====  ==========  ==========  ==========  ===========  ==========  ==========  =========  =========  =========  =========  =========  =============
    ! Band       1          2           3           4           5            6           7          8          9          10         11          12
    !       2200-1900   1900-1700   1700-1400   1400-1250    1250-1100   1100-980     980-800    800-670    670-540    540-400    400-280    280-0 
    ! ====  ==========  ==========  ==========  ===========  ==========  ==========  =========  =========  =========  =========  =========  =============
    !
    ! Two additional LW spectral intervals have been added in beyond 2200cm-1.
    !
    ! =====   ===========   ===========
    ! Band        13            14
    !          2500-2200     2850-2500
    ! =====   ===========   ===========
    !
    ! Emissivity ems(band(1))   from April data, Table2 of Chen et al
    ! 11th Conf Sat Met, Madison, WI, p 514
    ! here regoganized as 14 13 1 2 ... 12 above
    !
    ! :20 surface types:
    !
    ! ===================  ===================  ===================  =====================
    !  1= evergreen nleaf   2= evergreen bleaf   3= deciduous nleaf   4= deciduous bleaf
    !  5= mixed forests     6= closed shrubs     7= open shrubs       8= woody savanna
    !  9= savanna          10= grasslands       11= perma wet        12= croplands
    ! 13= urban            14= mosaic           15= snow             16= barren (deserts)
    ! 17= water            18= toundra          19= fresh snow       20= sea ice
    ! ===================  ===================  ===================  =====================
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: nchn              ! number of bands for which emissivity is needed
    real(8), intent(in)  :: waven(nchn)       ! wavenumbers (cm-1)
    real(8), intent(out) :: emi_mat(nchn, 20) ! emissivity (0.0-1.0)

    ! Locals:
    integer          :: bandIndex, channelIndex, typeIndex
    real(8)          :: dum

    ! CERES bands central wavenumber (covers 3.7 micron to 71.4 mic)
    integer, parameter :: nb=14
    real(8), parameter :: mid(nb) =(/                                             &
         2675.d0, 2350.d0, 2050.d0, 1800.d0, 1550.d0, 1325.d0, 1175.d0, 1040.d0,  &
         890.d0,  735.d0,  605.d0,  470.d0,  340.d0,  140.d0 /)

    ! CERES emissivity per wavenumber and surface types
    real(8), parameter ::  emi_tab(nb,20)=(/                                      &
         0.951d0, 0.989d0, 0.989d0, 0.989d0, 0.990d0, 0.991d0, 0.991d0, 0.990d0,  &
         0.990d0, 0.995d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                    &
         0.956d0, 0.989d0, 0.989d0, 0.989d0, 0.990d0, 0.991d0, 0.991d0, 0.990d0,  &
         0.990d0, 0.995d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                    &
         0.929d0, 0.985d0, 0.985d0, 0.986d0, 0.984d0, 0.983d0, 0.979d0, 0.980d0,  &
         0.973d0, 0.987d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                    &
         0.943d0, 0.985d0, 0.985d0, 0.986d0, 0.984d0, 0.983d0, 0.979d0, 0.980d0,  &
         0.973d0, 0.987d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                    &
         0.945d0, 0.987d0, 0.987d0, 0.987d0, 0.987d0, 0.987d0, 0.985d0, 0.985d0,  &
         0.982d0, 0.991d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                    &
         0.933d0, 0.949d0, 0.949d0, 0.970d0, 0.974d0, 0.971d0, 0.947d0, 0.958d0,  &
         0.966d0, 0.975d0, 0.984d0, 0.984d0, 0.984d0, 0.984d0,                    &
         0.873d0, 0.873d0, 0.873d0, 0.934d0, 0.944d0, 0.939d0, 0.873d0, 0.904d0,  &
         0.936d0, 0.942d0, 0.951d0, 0.951d0, 0.951d0, 0.951d0,                    &
         0.930d0, 0.987d0, 0.987d0, 0.990d0, 0.992d0, 0.993d0, 0.983d0, 0.975d0,  &
         0.985d0, 0.993d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                    &
         0.926d0, 0.987d0, 0.987d0, 0.990d0, 0.992d0, 0.993d0, 0.983d0, 0.975d0,  &
         0.985d0, 0.993d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                    &
         0.899d0, 0.987d0, 0.987d0, 0.990d0, 0.992d0, 0.993d0, 0.983d0, 0.975d0,  &
         0.985d0, 0.993d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                    &
         0.951d0, 0.983d0, 0.983d0, 0.987d0, 0.987d0, 0.988d0, 0.983d0, 0.981d0,  &
         0.987d0, 0.982d0, 0.986d0, 0.986d0, 0.986d0, 0.986d0,                    &
         0.924d0, 0.987d0, 0.987d0, 0.990d0, 0.992d0, 0.993d0, 0.983d0, 0.975d0,  &
         0.985d0, 0.993d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                    &
         0.929d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,  &
         1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                    &
         0.926d0, 0.987d0, 0.987d0, 0.989d0, 0.989d0, 0.990d0, 0.984d0, 0.980d0,  &
         0.983d0, 0.992d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                    &
         0.972d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,  &
         1.000d0, 0.999d0, 0.999d0, 0.999d0, 0.999d0, 0.999d0,                    &
         0.866d0, 0.835d0, 0.835d0, 0.916d0, 0.934d0, 0.923d0, 0.835d0, 0.877d0,  &
         0.921d0, 0.926d0, 0.934d0, 0.934d0, 0.934d0, 0.934d0,                    &
         0.973d0, 0.979d0, 0.979d0, 0.983d0, 0.982d0, 0.982d0, 0.984d0, 0.987d0,  &
         0.989d0, 0.972d0, 0.972d0, 0.972d0, 0.972d0, 0.972d0,                    &
         0.968d0, 0.947d0, 0.947d0, 0.967d0, 0.988d0, 0.979d0, 0.975d0, 0.977d0,  &
         0.992d0, 0.989d0, 0.989d0, 0.989d0, 0.989d0, 0.989d0,                    &
         0.984d0, 0.988d0, 0.988d0, 0.988d0, 0.988d0, 0.988d0, 0.988d0, 0.988d0,  &
         0.988d0, 0.988d0, 0.988d0, 0.988d0, 0.988d0, 0.988d0,                    &
         0.964d0, 0.979d0, 0.979d0, 0.979d0, 0.979d0, 0.979d0, 0.979d0, 0.979d0,  &
         0.979d0, 0.979d0, 0.979d0, 0.979d0, 0.979d0, 0.979d0  /)

    do typeIndex = 1, 20
      do channelIndex = 1, nchn
        if (waven(channelIndex) > mid(1)) then
          emi_mat(channelIndex,typeIndex) = emi_tab(1,typeIndex)
        else if (waven(channelIndex) < mid(nb)) then
          emi_mat(channelIndex,typeIndex) = emi_tab(nb,typeIndex)
        else
          do bandIndex = 1, nb - 1
            if (waven(channelIndex) <= mid(bandIndex) .and. waven(channelIndex) >= mid(bandIndex + 1)) then
              dum = ( waven(channelIndex) - mid(bandIndex) ) / ( mid(bandIndex + 1) - mid(bandIndex) )
              emi_mat(channelIndex,typeIndex) = emi_tab(bandIndex,typeIndex) + &
                  (emi_tab(bandIndex + 1,typeIndex) - emi_tab(bandIndex,typeIndex)) * dum
              exit
            end if
          end do
        end if
      end do
    end do

  end subroutine ceres_ematrix

  !--------------------------------------------------------------------------
  !  emi_sea
  !--------------------------------------------------------------------------
  subroutine emi_sea(em_oc, wnum, angle, wind, nProfiles, nChannels)
    !
    ! :Purpose: GET OCEAN SURFACE EMISSIVITY
    ! :Note:    IMEM(NCHANNELS), set to zero initially, on next call IMEM will have the
    !           right boundary channel to save search time in interpolation.
    !           IOPT=1 means activate IMEM option (all calls ask for same channels)
    !
    !           To get surface ocean emissivity for a group of channels with
    !           wavenumbers WNUM (cm-1) looking at one point with surface
    !           wind speed wind from angle angle.
    !           Based on Masuda,1988, Remote Sens. of Envir, 313-329.
    !           Coded emissivity routine based on Masuda's data by Tom Kleespies
    !           Covers 650-2857 cm-1 or 3.1-15.4 microns
    !
    ! :CAUTION: extrapolated values from 769-650 cm-1
    !           and interpolated values between 2439-1250 cm-1
    !
    implicit none

    ! Arguments:
    integer, intent(in)    :: nProfiles                  ! Number of profiles
    integer, intent(in)    :: nChannels                  ! Number of channels
    real(8), intent(out)   :: em_oc(nChannels,nProfiles) ! Ocean emissivities (0.-1.)
    real(8), intent(in)    :: wnum(nChannels)            ! Channel wavenumbers (cm-1)
    real(8), intent(in)    :: angle(nProfiles)           ! Viewing angle (deg)
    real(8), intent(in)    :: wind(nProfiles)            ! Surface wind speed (m/s)

    ! Locals:
    integer     :: channelIndex, bandIndex, profileIndex
    integer     :: imem(nChannels) 
    integer     :: mchan(2)
    real(8)     :: dum
    real(8)     :: emi2(2,nProfiles)

    ! Masuda's 19 wavelengths converted to wavenumber
    real(8), parameter :: refw(19)=[ 2857.1d0, 2777.7d0, 2702.7d0, 2631.6d0, 2564.1d0, &
         2500.0d0, 2439.0d0, 1250.0d0, 1190.5d0, 1136.3d0,                             &
         1087.0d0, 1041.7d0, 1000.0d0, 952.38d0, 909.09d0,                             &
         869.57d0, 833.33d0, 800.00d0, 769.23d0]

    ! imem options
    imem(:) = 0

    do channelIndex = 1, nChannels

      !  out of range
      if (wnum(channelIndex) < 645.d0 .or. wnum(channelIndex) > refw(1)) then
        write(*,'(A,1x,e12.4)') ' fatal: wavenumber out of range in emi_sea', wnum(channelIndex)
        stop
      else if (wnum(channelIndex) <= refw(19) .and. wnum(channelIndex) > 645.d0) then
        !  extrapolated from 769 cm-1 to 645 cm-1: NOT FROM REAL DATA
        !  nevertheless thought to be much better than unity
        !  this is a region of relatively rapid emissivity change
        !  worst estimates for 700-645 cm-1, but these channels do not
        !  see the surface (strong co2 absorption).
        imem(channelIndex) = 18
      else
        !  CAUTION interpolation on large interval 1250-2439 cm-1
        !  where no data is available except that of ASTER. ASTER
        !  shows a relatively smooth variation with wavelength except
        !  for a sharp drop at 1600 cm-1 with highs at 1550 and 1650 cm-1
        !  with peak-to-peak variation of 1.5% in that narrow range.
        !  Worst estimates would be between 1400-1800 cm-1 in HIRS ch 12
        !  which only in very cold atmospheres sees the surface.
        do bandIndex = 1, 18
          if (wnum(channelIndex) > refw(bandIndex + 1) .and. wnum(channelIndex) <= refw(bandIndex)) then
            imem(channelIndex) = bandIndex
          end if
        end do

      end if
   
      mchan(1) = imem(channelIndex)
      mchan(2) = imem(channelIndex) + 1
 
      dum = ( wnum(channelIndex) - refw(mchan(1)) ) / ( refw(mchan(2)) - refw(mchan(1)) )

      call comp_ir_emiss(emi2, wind, angle, 2, nProfiles, mchan)

      ! interpolation/extrapolation in wavenumber 

      do profileIndex = 1, nProfiles
        em_oc(channelIndex,profileIndex) = emi2(1,profileIndex) + ( emi2(2,profileIndex) - emi2(1,profileIndex) ) * dum
      end do

    end do

  end subroutine emi_sea

  !--------------------------------------------------------------------------
  !  tvs_getCommonChannelSet
  !--------------------------------------------------------------------------
  subroutine tvs_getCommonChannelSet(channels, countUniqueChannel, listAll)
    !
    !:Purpose: get common channels among all MPI tasks
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: channels(:)        ! input channel list
    integer, intent(out) :: countUniqueChannel ! number of unique channels
    integer, intent(out) :: listAll(:)         ! output list of unque channels

    ! Locals:
    integer :: channelsb(tvs_maxChannelNumber)
    integer :: ierr, allChannelIndex, channelIndex
    integer, allocatable :: listGlobal(:)
    logical :: found
     
    if (size(channels) > tvs_maxChannelNumber) then
      write(*,*) 'You need to increase tvs_maxChannelNumber in tovsNL_mod !',size(channels), tvs_maxChannelNumber
      call utl_abort('tvs_getCommonChannelSet')
    end if

    if (mmpi_myid == 0) then
      allocate(listGlobal(mmpi_nprocs*tvs_maxChannelNumber))
    else
      allocate(listGlobal(1))
    end if

    listAll(:) = 0
    listGlobal(:) = 0
    channelsb(:) = 0
    channelsb(1:size(channels)) = channels(:)

    call rpn_comm_barrier('GRID',ierr)

    call rpn_comm_gather(channelsb, tvs_maxChannelNumber, 'MPI_INTEGER', listGlobal, &
         tvs_maxChannelNumber, 'MPI_INTEGER', 0, 'GRID', ierr) 
    countUniqueChannel = 0
    if (mmpi_myid == 0) then
      call isort(listGlobal, mmpi_nprocs*tvs_maxChannelNumber)
      do allChannelIndex = 1, mmpi_nprocs * tvs_maxChannelNumber
        if (listGlobal(allChannelIndex) > 0) then
          found = .false.
          LOOPJ: do channelIndex = countUniqueChannel, 1, -1
            if (listGlobal(allChannelIndex) == listAll(channelIndex)) then
              found =.true.
              exit LOOPJ
            end if
          end do LOOPJ
          if (.not. found) then
            countUniqueChannel = countUniqueChannel + 1
            listAll(countUniqueChannel) = listGlobal(allChannelIndex)
          end if
        end if
      end do
    end if
    
    call rpn_comm_bcast(countUniqueChannel, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(listAll(1:countUniqueChannel), countUniqueChannel, 'MPI_INTEGER', 0, 'GRID', ierr)

    deallocate(listGlobal)

  end subroutine tvs_getCommonChannelSet

  !--------------------------------------------------------------------------
  !   tvs_printDetailledOmfStatistics
  !--------------------------------------------------------------------------
  subroutine tvs_printDetailledOmfStatistics(obsSpaceData)
    !
    ! :Purpose: Print channel by channnel O-F statistics fro radiances
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData! obsSpacaData structure

    ! Locals:
    integer :: sensorIndex, channelIndex
    real(8) :: zjoch(0:tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    real(8) :: zavgnrm(0:tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    real(pre_obsReal) :: zdtb, obsPRM
    integer :: nchanperline, startChannel, endChannel
    integer :: count, incanjo
    integer :: idatyp
    integer :: rttovChannelNumber, bufrChannelNumber
    integer :: inobsch(0:tvs_maxChannelNumber, tvs_maxNumberOfSensors)
    integer :: lcanjo(tvs_maxChannelNumber)
    integer :: headerIndex, bodyIndex
    real(8) :: sigmaObs

    write(*,*) 'tvs_printDetailledOmfStatistics: Starting'

    if (tvs_headerEnd < 0) return    ! exit if there are not tovs data

    ! 1.  Computation of (hx - z)/sigma for tovs data only

    count  = 0
    inobsch(:,:) = 0
    zjoch  (:,:) = 0.0d0
    zavgnrm(:,:) = 0.0d0

    ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! 1.1  Extract general information for this observation point
      !      ------------------------------------------------------

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if (.not. tvs_isIdBurpTovs(idatyp)) then
        write(*,*) 'tvs_printDetailledOmfStatistics: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER
      end if
       
      sensorIndex = tvs_lsensor(headerIndex)

      ! Set the body list
      ! (& start at the beginning of the list)
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      count = 0
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY
        
        ! Only consider if flagged for assimilation
        if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated) cycle BODY                

        call tvs_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                           rttovChannelNumber, channelIndex )
        bufrChannelNumber = rttovChannelNumber + tvs_channelOffset(sensorIndex)
        if (channelIndex == 0) then
          write(*,'(A)') '  tvs_printDetailledOmfStatistics: error with channel number'
          call utl_abort(' tvs_printDetailledOmfStatistics')
        end if

        zdtb = obs_bodyElem_r(obsSpaceData,OBS_PRM,bodyIndex) - &
             tvs_radiance(headerIndex) % bt(channelIndex)
        if (tvs_debug) then
          obsPRM = obs_bodyElem_r(obsSpaceData,OBS_PRM,bodyIndex)
          write(*,'(a,i4,2f8.2,f6.2)') ' rttovChannelNumber,sim,obs,diff= ', &
               rttovChannelNumber, tvs_radiance(headerIndex) % bt(channelIndex), &
               obsPRM, -zdtb
        end if

        sigmaObs = obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)

        if (sigmaObs == MPC_missingValue_R8) cycle body

        count = count + 1
        inobsch(bufrChannelNumber,sensorIndex) = inobsch(bufrChannelNumber,sensorIndex) + 1
        zjoch(bufrChannelNumber,sensorIndex)   = &
             zjoch(bufrChannelNumber,sensorIndex) &
             + zdtb * zdtb / (sigmaObs * sigmaObs)
        zavgnrm(bufrChannelNumber,sensorIndex)   = &
             zavgnrm(bufrChannelNumber,sensorIndex) - &
             zdtb / sigmaObs
      end do BODY

    end do HEADER

    !   2.  Close up, print summary
    !   .   -----------------------

    ! printout of mean jo and normalized average for each sensor.

    nchanperline = 18
    if (count > 0) then
      write(*,*)
      write(*,*)
      write(*,'(10x,A)') '- tvs_printDetailledOmfStatistics: computing jo and residuals to tovs  observations'

      do sensorIndex = 1, tvs_nsensors
        inobsch(0,sensorIndex) = sum ( inobsch(1:,sensorIndex) )
        zjoch(0,sensorIndex) = sum( zjoch(1:,sensorIndex) )
        zavgnrm(0,sensorIndex) = sum( zavgnrm(1:,sensorIndex) )
      end do

      do sensorIndex = 1, tvs_nsensors
        incanjo = 0
        do channelIndex = 0, tvs_maxChannelNumber
          if (inobsch(channelIndex, sensorIndex) /= 0) then
            incanjo = incanjo + 1
            lcanjo(incanjo) = channelIndex
          end if
        end do
        if (incanjo /= 0) then
          write(*,"(/1x,'sensor #',i2,'. platform: ',a, 'instrument: ',a)") &
               sensorIndex, tvs_satelliteName(sensorIndex), tvs_instrumentName(sensorIndex)
          do startChannel = 1, incanjo, nchanperline
            endChannel = min(startChannel + nchanperline - 1 , incanjo)
            if (startChannel == 1) then
              write(*,"(1x,'channel',t13,'   all',17i6)") (lcanjo(channelIndex), channelIndex=startChannel+1, endChannel)
            else
              write(*,"(1x,'channel',t13,18i6)") (lcanjo(channelIndex), channelIndex=startChannel, endChannel)
            end if
            write(*,"(1x,'no. obs.',t13,18i6)") (inobsch(lcanjo(channelIndex),sensorIndex), channelIndex=startChannel, endChannel)
            write(*,"(1x,'mean jo',t13,18f6.2)") &
                 (zjoch(lcanjo(channelIndex),sensorIndex)/max(1,inobsch(lcanjo(channelIndex),sensorIndex)), channelIndex=startChannel,endChannel)
            write(*,"(1x,'norm. bias',t13,18f6.2,/)") &
                 (zavgnrm(lcanjo(channelIndex),sensorIndex)/max(1,inobsch(lcanjo(channelIndex), sensorIndex)) , channelIndex=startChannel, endChannel)
          end do
        end if
      end do
    end if

  end subroutine  tvs_printDetailledOmfStatistics

  !--------------------------------------------------------------------------
  !  tvs_getLocalChannelIndexFromChannelNumber
  !--------------------------------------------------------------------------
  subroutine tvs_getLocalChannelIndexFromChannelNumber(idsat,channelIndex_out,channelNumber_in)
    !
    ! :Purpose: to get local channel index from channel number
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: idsat            ! Satellite index
    integer, intent(out) :: channelIndex_out ! Channel index
    integer, intent(in)  :: channelNumber_in ! Channel number

    ! Locals:
    logical, save              :: firstCall =.true.
    integer                    :: channelNumber, sensorIndex, channelIndex 
    integer, allocatable, save :: savedChannelIndexes(:,:)

    if (firstCall) then
      allocate(savedChannelIndexes(tvs_nsensors, tvs_maxChannelNumber ))
      savedChannelIndexes(:,:) = -1
      do sensorIndex = 1, tvs_nsensors
        channels:do channelNumber = 1, tvs_maxChannelNumber
          indexes: do channelIndex = 1, tvs_nchan(sensorIndex)
            if (channelNumber == tvs_ichan(channelIndex,sensorIndex)) then
              savedChannelIndexes(sensorIndex,channelNumber) = channelIndex
              exit indexes
            end if
          end do indexes
        end do channels
      end do
      firstCall = .false.
    end if

    channelIndex_out = savedChannelIndexes(idsat,channelNumber_in)

    if (channelIndex_out == -1) then
      write(*,*) 'channel number requested = ', channelNumber_in
      write(*,*) 'idsat = ', idsat
      write(*,*) 'tvs_getLocalChannelIndexFromChannelNumber: warning channel not found'  
    end if

  end subroutine tvs_getLocalChannelIndexFromChannelNumber

  !--------------------------------------------------------------------------
  !  updateCloudInTovsProfile
  !--------------------------------------------------------------------------
  subroutine updateCloudInTovsProfile(sensorHeaderIndexes, nlv_T, mode, beSilent)
    !
    ! :Purpose: Modify the cloud in tvs_profiles_nl structure of rttov.
    !
    implicit none
    
    ! Arguments:
    integer,      intent(in) :: sensorHeaderIndexes(:) ! indexes of radiance observations for the currently processed sensor
    integer,      intent(in) :: nlv_T                  ! number of model vertical thermodynamical levels
    character(*), intent(in) :: mode                   ! save or restore
    logical,      intent(in) :: beSilent               ! flag to control verbosity

    ! Locals:
    integer :: profileIndex, profileCount
    real(8), allocatable, save :: cloudProfileToStore(:,:)

    if (.not. beSilent) write(*,*) 'updateCloudInTovsProfile: Starting'
    if (.not. beSilent) call msg_memUsage('updateCloudInTovsProfile')

    profileCount = size(sensorHeaderIndexes)

    if (trim(mode) == 'save') then 
      if (allocated(cloudProfileToStore)) deallocate(cloudProfileToStore)
      allocate(cloudProfileToStore(nlv_T,profileCount))

      do profileIndex = 1, profileCount
        cloudProfileToStore(:,profileIndex) = tvs_profiles_nl(sensorHeaderIndexes(profileIndex)) % clw(:)
        tvs_profiles_nl(sensorHeaderIndexes(profileIndex)) % clw(:) = qlim_getMinValueCloud('LWCR') 
      end do

    else if (trim(mode) == 'restore') then 
      do profileIndex = 1, profileCount
        tvs_profiles_nl(sensorHeaderIndexes(profileIndex)) % clw(:) = cloudProfileToStore(:,profileIndex)
      end do

      deallocate(cloudProfileToStore)

    else
      call utl_abort('updateCloudInTovsProfile: mode should be either "save" or "restore"')
    end if

  end subroutine updateCloudInTovsProfile

  !--------------------------------------------------------------------------
  !  updateCloudInTovsCloudProfile
  !--------------------------------------------------------------------------
  subroutine updateCloudInTovsCloudProfile(sensorHeaderIndexes, nlv_T, mode, beSilent)
    !
    ! :Purpose: Modify the cloud in tvs_cld_profiles_nl structure of rttovScatt.
    !
    implicit none
    
    ! Arguments:
    integer,      intent(in) :: sensorHeaderIndexes(:) ! indexes of radiance observations for the currently processed sensor
    integer,      intent(in) :: nlv_T                  ! number of model vertical thermodynamical levels   
    character(*), intent(in) :: mode                   ! save or restore
    logical,      intent(in) :: beSilent               ! flag to control verbosity

    ! Locals:
    integer :: profileIndex, profileCount, headerIndex
    real(8), allocatable, save :: rainFluxProfileToStore(:,:)
    real(8), allocatable, save :: snowFluxProfileToStore(:,:)
    real(8), allocatable, save :: clwProfileToStore(:,:)
    real(8), allocatable, save :: ciwProfileToStore(:,:)
    real(8), allocatable, save :: cloudFractionProfileToStore(:,:)

    profileCount = size(sensorHeaderIndexes)

    if (.not. beSilent) write(*,*) 'updateCloudInTovsCloudProfile: Starting', profileCount
    if (.not. beSilent) call msg_memUsage('updateCloudInTovsCloudProfile')

    if (trim(mode) == 'save') then 
      if (allocated(rainFluxProfileToStore)) deallocate(rainFluxProfileToStore)
      if (allocated(snowFluxProfileToStore)) deallocate(snowFluxProfileToStore)
      if (allocated(clwProfileToStore)) deallocate(clwProfileToStore)
      if (allocated(ciwProfileToStore)) deallocate(ciwProfileToStore)
      if (allocated(cloudFractionProfileToStore)) deallocate(cloudFractionProfileToStore)
      allocate(rainFluxProfileToStore(nlv_T,profileCount))
      allocate(snowFluxProfileToStore(nlv_T,profileCount))
      allocate(clwProfileToStore(nlv_T,profileCount))
      allocate(ciwProfileToStore(nlv_T,profileCount))
      allocate(cloudFractionProfileToStore(nlv_T,profileCount))

      do profileIndex = 1, profileCount
        headerIndex = sensorHeaderIndexes(profileIndex)
        rainFluxProfileToStore(:,profileIndex) = tvs_cld_profiles_nl(headerIndex) % hydro(:,1)
        snowFluxProfileToStore(:,profileIndex) = tvs_cld_profiles_nl(headerIndex) % hydro(:,2)
        clwProfileToStore(:,profileIndex) = tvs_cld_profiles_nl(headerIndex) % hydro(:,4)
        ciwProfileToStore(:,profileIndex) = tvs_cld_profiles_nl(headerIndex) % hydro(:,5)
        cloudFractionProfileToStore(:,profileIndex) = tvs_cld_profiles_nl(headerIndex) % hydro_frac(:,1)
        tvs_cld_profiles_nl(headerIndex) % hydro(:,1) = qlim_getMinValueCloud('RF')
        tvs_cld_profiles_nl(headerIndex) % hydro(:,2) = qlim_getMinValueCloud('SF')
        tvs_cld_profiles_nl(headerIndex) % hydro(:,4) = qlim_getMinValueCloud('LWCR')
        tvs_cld_profiles_nl(headerIndex) % hydro(:,5) = qlim_getMinValueCloud('IWCR')
        tvs_cld_profiles_nl(headerIndex) % hydro_frac(:,1) = qlim_getMinValueCloud('CLDR')
      end do

    else if (trim(mode) == 'restore') then 
      do profileIndex = 1, profileCount
        headerIndex = sensorHeaderIndexes(profileIndex)
        tvs_cld_profiles_nl(headerIndex) % hydro(:,1) = rainFluxProfileToStore(:,profileIndex)
        tvs_cld_profiles_nl(headerIndex) % hydro(:,2) = snowFluxProfileToStore(:,profileIndex)
        tvs_cld_profiles_nl(headerIndex) % hydro(:,4) = clwProfileToStore(:,profileIndex)
        tvs_cld_profiles_nl(headerIndex) % hydro(:,5) = ciwProfileToStore(:,profileIndex)
        tvs_cld_profiles_nl(headerIndex) % hydro_frac(:,1) = cloudFractionProfileToStore(:,profileIndex)
      end do

      deallocate(rainFluxProfileToStore)
      deallocate(snowFluxProfileToStore)
      deallocate(clwProfileToStore)
      deallocate(ciwProfileToStore)
      deallocate(cloudFractionProfileToStore)

    else
      call utl_abort('updateCloudInTovsCloudProfile: mode should be either "save" or "restore"')
    end if

  end subroutine updateCloudInTovsCloudProfile

  !--------------------------------------------------------------------------
  !  tvs_getChannelNumIndexFromPPP
  !--------------------------------------------------------------------------
  subroutine tvs_getChannelNumIndexFromPPP( obsSpaceData, headerIndex, bodyIndex, &
                                            channelNumber, channelIndex )
    !
    ! :Purpose: Get channel number/index from obs_ppp for TO observations.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(in)  :: obsSpaceData  ! obsSpaceData structure
    integer,          intent(in)  :: headerIndex   ! observation index in obsSpaceData header table
    integer,          intent(in)  :: bodyIndex     ! observation index in obsSpaceData body table
    integer,          intent(out) :: channelNumber ! channel number
    integer,          intent(out) :: channelIndex  ! channel index in tvs_ichan

    ! Locals:
    integer :: sensorIndex

    sensorIndex = tvs_lsensor(headerIndex)

    channelNumber = nint(obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex))
    channelNumber = max( 0 , min( channelNumber , tvs_maxChannelNumber + 1))
    channelNumber = channelNumber - tvs_channelOffset(sensorIndex)
    channelIndex = utl_findloc(tvs_ichan(:,sensorIndex),channelNumber)

  end subroutine tvs_getChannelNumIndexFromPPP

  !--------------------------------------------------------------------------
  !  tvs_writeJacobianAscii
  !--------------------------------------------------------------------------
  subroutine tvs_writeJacobianAscii(jacobian, jacobian_emiss, profiles, chanprof, obsSpaceData, satelliteName, instrumentName, &
                                    bodyIndexFromBtIndex, sensorHeaderIndexes, btCount)
    !
    ! :Purpose: Write the computed Jacobian into ASCII files
    !
    implicit none

    ! Arguments:
    type(rttov_profile), pointer,    intent(in)    :: jacobian(:)              ! Jacobian
    type(rttov_emissivity), pointer, intent(in)    :: jacobian_emiss(:)        ! Surface Emissivity Jacobian
    integer,                         intent(in)    :: bodyIndexFromBtIndex(:)  ! Provides the bodyIndex in ObsSpaceData based on btIndex
    integer,                         intent(in)    :: btCount                  ! Total number of simulated radiances
    type(struct_obs),                intent(in)    :: obsSpaceData             ! ObsSpaceData structure
    character(len=15),               intent(in)    :: satelliteName            ! Satellite Name
    character(len=15),               intent(in)    :: instrumentName           ! Instrument Name
    type (rttov_profile), pointer,   intent(in)    :: profiles(:)              ! Input profiles from background state
    type (rttov_chanprof),           intent(in)    :: chanprof(:)              ! Chanprof structure    
    integer,                         intent(in)    :: sensorHeaderIndexes(:)   ! Sensor obsSpaceData header indexes
   
    ! Locals:
    character(len=4)               :: cmyidx, cmyidy, strNumLev
    character(len=9)               :: cmyid
    character(len=1024)            :: fileName
    integer                        :: btIndex, bodyIndex
    integer(8)                     :: obsIdd, obsIdo
    integer                        :: profileIndex, headerIndex
    integer                        :: err, iunit, numLev
    integer, external              :: fnom,fclos
    character(len = 12), parameter :: dirName = 'tvs_jacobian'

    err = clib_isdir(trim(dirName))

    ! Create directory if it doesn't exists
    if (err /= CLIB_OK ) then
      err = clib_mkdir_r(trim(dirName))
    end if

    write(cmyidy,'(I4.4)') (mmpi_myidy + 1)
    write(cmyidx,'(I4.4)') (mmpi_myidx + 1)
    cmyid = trim(cmyidx) // '_' // trim(cmyidy)

    fileName = 'tvs_jacobian_' // trim(satelliteName) //'_'// trim(instrumentName) //'_'// cmyid

    iunit = 0
    err = fnom(iunit, trim(dirName) // '/' // trim(fileName), 'FTN+SEQ+R/W', 0)
    if (err /= 0) then 
      call utl_abort('tvs_writeJacobianAscii: Error writing Jacobian files')
    end if 

    do btIndex = 1, btCount
      profileIndex = chanprof(btIndex) % prof
      headerIndex = sensorHeaderIndexes(profileIndex)
      obsIdo = obs_headPrimaryKey(obsSpaceData, headerIndex)
      bodyIndex = bodyIndexFromBtIndex(btIndex)

      if (bodyIndex > 0) then
        obsIdd = obs_bodyPrimaryKey(obsSpaceData, bodyIndex)

        if (size(profiles(headerIndex) % p(:)) /= size(jacobian(btIndex) % t(:)) .or. &
            size(profiles(headerIndex) % p(:)) /= size(jacobian(btIndex) % q(:)) .or. &
            size(profiles(headerIndex) % p(:)) /= size(jacobian(btIndex) % p(:))) then
          call utl_abort('tvs_writeJacobianAscii: Number of pressure levels does not match ' // &
                          'the number of model levels in Jacobian')
        end if

        numLev = size(profiles(headerIndex) % p(:))
        write (strNumLev,'(I4)') numLev

        write(iunit,'(I20, I20, F16.2, F16.2, I4, ' &
                      // trim(strNumLev) // 'E16.5E2, ' // &
                      'E16.5E2, E16.5E2, E16.5E2, E16.5E2, ' &
                      // trim(strNumLev) // 'E16.5E2,' &
                      // trim(strNumLev) // 'E16.5E2,' &
                      // trim(strNumLev) // 'E16.5E2)') &
              obsIdo, obsIdd, profiles(headerIndex) % latitude, profiles(headerIndex) % longitude, numLev, &
              profiles(headerIndex) % p(:), &
              jacobian_emiss(btIndex) % emis_out, jacobian(btIndex) % skin % t, jacobian(btIndex) % s2m % t, jacobian(btIndex) % s2m % p, &
              jacobian(btIndex) % t(:), jacobian(btIndex) % q(:), jacobian(btIndex) % p(:)
      end if
    end do

    err = fclos(iunit)
    
  end subroutine tvs_writeJacobianAscii

  !--------------------------------------------------------------------------
  !  tvs_rttovScatt
  !--------------------------------------------------------------------------
  subroutine tvs_rttovScatt(errorstatus, opts_scatt, nlevels, chanprof, frequencies, &
      profiles, cld_profiles, coef_rttov, coef_scatt, calcemis, emissivity, transmission, &
      radiance, cfrac, emis_retrieval_terms, reflectivity)
    !
    ! :Purpose: modified rttov_scatt subroutine with additional transmission argument
    !
    implicit none

    ! Arguments:
    integer(kind=jpim),                              intent(out)    :: errorstatus                   ! Error return flag
    type(rttov_options_scatt),                       intent(in)     :: opts_scatt                    ! RTTOV-SCATT options
    integer(kind=jpim),                              intent(in)     :: nlevels                       ! Number of levels
    type(rttov_chanprof),                            intent(in)     :: chanprof(:)                   ! Channel and profile indices
    integer(kind=jpim),                              intent(in)     :: frequencies (size(chanprof))  ! Frequency indices
    type(rttov_profile),                             intent(in)     :: profiles(:)                   ! thermodynamic profiles
    type(rttov_profile_cloud),                       intent(in)     :: cld_profiles (size(profiles)) ! Cloud profiles
    type(rttov_coefs),                               intent(in)     :: coef_rttov                    ! RTTOV Coefficients
    type(rttov_scatt_coef),                          intent(in)     :: coef_scatt                    ! RTTOV_SCATT Coefficients
    logical(kind=jplm),                              intent(in)     :: calcemis   (size(chanprof))   ! Switch for emmissivity calculation
    type(rttov_emissivity),                          intent(inout)  :: emissivity (size(chanprof))   ! Surface emissivity
    type(rttov_transmission),                        intent(inout)  :: transmission                  ! Transmittances
    type(rttov_radiance),                            intent(inout)  :: radiance                      ! Radiances
    real(kind=jprb), optional,                       intent(out)    :: cfrac (size(profiles))        ! Cloud fraction (diagnostic)
    type(rttov_scatt_emis_retrieval_type), optional, intent(inout)  :: emis_retrieval_terms          ! Optional for all-sky emis retrievals
    type(rttov_reflectivity), optional,              intent(inout)  :: reflectivity                  ! Optional for radar

    ! Locals:
    integer(kind=jpim), target :: sa__mclayer(size(chanprof))
    real(kind=jprb), target :: sa__cfrac   (size(profiles))
    real(kind=jprb), target :: sa__ems_bnd (size(chanprof))
    real(kind=jprb), target :: sa__ref_bnd (size(chanprof))
    real(kind=jprb), target :: sa__ems_cld (size(chanprof))
    real(kind=jprb), target :: sa__ref_cld (size(chanprof))
    real(kind=jprb), target :: sa__tbd (size(chanprof),nlevels+1)
    real(kind=jprb), target :: sa__tsfc (size(chanprof))
    real(kind=jprb), target :: sa__tcosmic (size(chanprof))
    real(kind=jprb), target :: sa__delta  (size(chanprof),nlevels)
    real(kind=jprb), target :: sa__tau    (size(chanprof),nlevels)
    real(kind=jprb), target :: sa__int_tau(size(chanprof),nlevels)
    real(kind=jprb), target :: sa__ext    (size(chanprof),nlevels)
    real(kind=jprb), target :: sa__ssa    (size(chanprof),nlevels)
    real(kind=jprb), target :: sa__asm    (size(chanprof),nlevels)
    real(kind=jprb), target :: sa__zef    (size(chanprof),nlevels)
    real(kind=jprb), target :: sa__lambda (size(chanprof),nlevels)
    real(kind=jprb), target :: sa__h      (size(chanprof),nlevels)
    real(kind=jprb), target :: sa__b0     (size(chanprof),nlevels)
    real(kind=jprb), target :: sa__b1     (size(chanprof),nlevels)
    real(kind=jprb), target :: sa__bn     (size(chanprof),nlevels)
    real(kind=jprb), target :: sa__btop   (size(chanprof))
    real(kind=jprb), target :: sa__bsfc   (size(chanprof))
    real(kind=jprb), target :: sa__dz     (size(profiles),nlevels)
    real(kind=jprb), target :: sa__hydro  (size(profiles),nlevels,cld_profiles(1)%nhydro)
    real(kind=jprb), target :: sf__t_down (size(chanprof))
    real(kind=jprb), target :: sf__t_up   (size(chanprof))
    real(kind=jprb), target :: sf__tau    (size(chanprof))
    real(kind=jprb), target :: sr__upclear(size(chanprof))
    real(kind=jprb), target :: sr__dnclear(size(chanprof))
    real(kind=jprb), target :: sr__refldnclear(size(chanprof))
    real(kind=jprb), target :: sr__up  (nlevels-1, size(chanprof))
    real(kind=jprb), target :: sr__down(nlevels-1, size(chanprof))
    real(kind=jprb), target :: sr__surf(nlevels-1, size(chanprof))
    integer(kind=jpim) :: nprofiles, nchannels
    logical(kind=jplm) :: lpolarimetric(size(chanprof)), lthermal(size(chanprof))
    integer(kind=jpim) :: iprof, ichan, ilayer
    integer(kind=jpim) :: ichan_act
    real(kind=jprb)    :: rad_cld(size(chanprof))
    real(kind=jprb)    :: zlayers(nlevels) !altitude from ground in km
    logical(kind=jplm) :: lreflectivity
    
    ! Variables for emissivity calculations
    type (eddington_sfc_type) :: sfc_terms ! Downward radiance source terms, Upward radiance source terms, Total transmittancs
    type (rttov_radiance2)    :: radiance2               
    type (rttov_geometry)     :: angles (size(profiles))
    type (rttov_profile_scatt_aux) :: scatt_aux  
    type (rttov_options)    :: opts
    character (len=80) :: errMessage
    character (len=15) :: NameOfRoutine = 'tvs_rttovScatt '

    lreflectivity = present(reflectivity)

    nprofiles = size(profiles)
    nchannels = size(chanprof)

    errorstatus = errorstatus_success
  
    !we want to be able to get trasmittance as an argument. It is no longer possible to use automatic allocation
    allocate(transmission % tau_total(size(chanprof)))
    allocate(transmission % tau_levels(nlevels,size(chanprof))) 
 
    scatt_aux % cfrac    => sa__cfrac
    scatt_aux % ems_bnd  => sa__ems_bnd
    scatt_aux % ref_bnd  => sa__ref_bnd
    scatt_aux % ems_cld  => sa__ems_cld
    scatt_aux % ref_cld  => sa__ref_cld
    scatt_aux % tbd      => sa__tbd
    scatt_aux % tsfc     => sa__tsfc
    scatt_aux % tcosmic  => sa__tcosmic
    scatt_aux % mclayer  => sa__mclayer
    scatt_aux % delta    => sa__delta
    scatt_aux % tau      => sa__tau
    scatt_aux % ext      => sa__ext
    scatt_aux % ssa      => sa__ssa
    scatt_aux % asm      => sa__asm
    scatt_aux % int_tau  => sa__int_tau
    scatt_aux % zef      => sa__zef
    scatt_aux % lambda   => sa__lambda
    scatt_aux % h        => sa__h
    scatt_aux % b0       => sa__b0
    scatt_aux % b1       => sa__b1
    scatt_aux % bn       => sa__bn
    scatt_aux % btop     => sa__btop
    scatt_aux % bsfc     => sa__bsfc
    scatt_aux % dz       => sa__dz
    scatt_aux % hydro    => sa__hydro
    
    sfc_terms % down   => sf__t_down
    sfc_terms % up     => sf__t_up
    sfc_terms % tau    => sf__tau
    
    radiance2 % upclear => sr__upclear
    radiance2 % dnclear => sr__dnclear
    radiance2 % refldnclear => sr__refldnclear
    radiance2 % up      => sr__up  
    radiance2 % down    => sr__down
    radiance2 % surf    => sr__surf
   
    ! Check inputs
    ! ------------
    do iprof = 1, nprofiles
      if (profiles(iprof) % s2m % p /= cld_profiles(iprof) % ph(nlevels+1)) then
        errorstatus = errorstatus_fatal
        write( errMessage, '( "Surface pressure and lowest half level should be identical")' )
      end if
      if (cld_profiles(iprof) % nhydro /= coef_scatt % nhydro) then
        errorstatus = errorstatus_fatal
        write( errMessage, '( "Number of hydrometeors differs between inputs and scattering coefficients ")' )
      end if
      if (cld_profiles(iprof) % nhydro_frac /= coef_scatt % nhydro .and. &
          cld_profiles(iprof) % nhydro_frac /= 1_JPIM) then
        errorstatus = errorstatus_fatal
        write( errMessage, '( "Number of hydrometeor fractions should be 1 or nhydro ")' )
      end if
      if (errorstatus == errorstatus_fatal) then
        call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
        return
      end if
    end do

    ! Identify polarimetric channels for fix to use only the clear sky part of RT calculation 
    do ichan = 1, nchannels
      ichan_act = chanprof(ichan)%chan
      lpolarimetric(ichan) = ( (coef_rttov % coef % id_sensor == sensor_id_po) &
          .and.   (coef_rttov % coef % fastem_polar(ichan_act) + 1_jpim >= 6_jpim) )
    end do
      
    !* 1.   Gas absorption
    
    ! Profiles will be interpolated from model/RTTOV-SCATT levels to 
    ! RTTOV coefficient levels within RTTOV itself.
    opts % interpolation % addinterp   = .true.
    opts % rt_ir % addclouds           = .false.
    opts % rt_ir % addsolar            = .false.
    opts % rt_ir % addaerosl           = .false.
    opts % rt_ir % pc % addpc          = .false.
    opts % rt_ir % pc % addradrec      = .false.
    opts % rt_mw % clw_data            = .false.
    
    opts % rt_mw % fastem_version           = opts_scatt % fastem_version
    opts % rt_mw % supply_foam_fraction     = opts_scatt % supply_foam_fraction
    opts % rt_all % ozone_data              = opts_scatt % ozone_data
    opts % rt_all % use_t2m_opdep           = opts_scatt % use_t2m_opdep
    opts % rt_all % use_q2m                 = opts_scatt % use_q2m
    opts % rt_all % addrefrac               = opts_scatt % addrefrac
    opts % rt_all % rad_down_lin_tau        = opts_scatt % rad_down_lin_tau
    opts % rt_all % dtau_test               = opts_scatt % dtau_test
    opts % config                           = opts_scatt % config
    opts % interpolation % interp_mode      = opts_scatt % interp_mode
    opts % interpolation % reg_limit_extrap = opts_scatt % reg_limit_extrap
    opts % interpolation % lgradp           = opts_scatt % lgradp

    call rttov_direct(          &!
        errorstatus,            &! out
        chanprof,               &! in
        opts,                   &! in
        profiles,               &! in
        coef_rttov,             &! in
        transmission,           &! inout
        radiance,               &! inout
        radiance2 = radiance2,  &! inout 
        calcemis   = calcemis,  &! in
        emissivity = emissivity) ! inout

    if (errorstatus == errorstatus_fatal) Then
      write(errMessage, '( "error in rttov_direct")')
      call rttov_errorreport(errorstatus_fatal, errMessage, NameOfRoutine)
      return
    end if
    
    scatt_aux % ems_cld (:) = emissivity (:) % emis_in
    scatt_aux % ref_cld (:) = 1.0_JPRB - emissivity (:) % emis_in

    !* 2.   Initialisations for Eddington
    call rttov_iniscatt(       &
        errorstatus,           &! out
        opts,                  &! in
        opts_scatt,            &! in
        lreflectivity,         &! in
        nlevels,               &! in
        nchannels,             &! in
        nprofiles,             &! in
        chanprof,              &! in
        frequencies,           &! in
        profiles,              &! in
        cld_profiles,          &! in
        coef_rttov%coef,       &! in
        coef_scatt,            &! in
        transmission,          &! in
        calcemis,              &! in
        opts_scatt%lusercfrac, &! in
        angles,                &! out
        scatt_aux)              ! inout
    
    if (errorstatus == errorstatus_fatal) Then
      write(errMessage, '( "error in rttov_iniscatt")' )
      call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
      return
    end if
    
    if (.not. lreflectivity) then

      !* 3. Eddington (in temperature space)
      call rttov_eddington(        &
          opts_scatt%cc_threshold, &! in
          nlevels,                 &! in
          nchannels,               &! in
          nprofiles,               &! in
          chanprof,                &! in
          angles,                  &! in
          scatt_aux,               &! in
          rad_cld,                 &! out
          sfc_terms = sfc_terms)    ! inout, optional, Upward and downward radiance source terms, Total transmittances

      ! Emissivity retrieval terms
      if (present(emis_retrieval_terms)) then
        call rttov_scatt_emis_terms( &
            opts_scatt%cc_threshold, &! in
            chanprof,                &! in
            coef_rttov,              &! in
            scatt_aux,               &! in
            emissivity,              &! in
            transmission,            &! in
            radiance2,               &! in
            sfc_terms,               &! in
            emis_retrieval_terms)     ! inout
      end if

      !* 4.  Combine clear and cloudy parts
      do ichan = 1, nchannels
        iprof = chanprof(ichan)%prof
        
        if (scatt_aux % cfrac (iprof) > opts_scatt % cc_threshold .and. .not. lpolarimetric(ichan)) then
          radiance % total (ichan) = rad_cld (ichan) * scatt_aux % cfrac (iprof) &
              + radiance % clear (ichan) * (1.0_JPRB - scatt_aux % cfrac (iprof))
        else
          radiance % total (ichan) = radiance % clear (ichan)
        end if
      end do
      
      ! Return the cloud fraction actually used - this is diagnostic output
      ! only provided by the forward model.
      if (present(cfrac)) then
        cfrac(:) = scatt_aux % cfrac (:)
      end if
      
    else
      
      reflectivity % zef(:,:)          = 0._JPRB
      reflectivity % azef(:,:)         = 0._JPRB
      radiance % geometric_height(:,:) = 0._JPRB
      radiance % total (:)             = min_radiance_radar
      radiance % clear (:)             = min_radiance_radar
      
      do ichan = 1, nchannels
        iprof            = chanprof(ichan) % prof
        zlayers(nlevels) = profiles(iprof) % elevation
        
        do ilayer = nlevels, 1, -1
          
          if (ilayer < nlevels) zlayers(ilayer) = zlayers(ilayer+1) + scatt_aux%dz(iprof,ilayer+1)
          
          if (scatt_aux % zef (ichan,ilayer) > min_reflectivity) then
            
            reflectivity % zef (ilayer,ichan)  = scatt_aux % zef (ichan,ilayer)
            
            reflectivity % azef (ilayer,ichan) = scatt_aux % zef (ichan,ilayer) + &
                2*10*log10( scatt_aux % int_tau (ichan,ilayer) )
            
          else
            reflectivity % zef (ilayer,ichan)  = min_reflectivity
            reflectivity % azef (ilayer,ichan) = min_reflectivity
          end if
          
          ! Approximate altitude of the middle of the RTTOV-SCATT layer (corresponding to the IFS full level), back from km to m
          radiance % geometric_height (ilayer,ichan) = 1000.0_JPRB * (zlayers(ilayer) + scatt_aux % dz(iprof,ilayer)/2.)
        end do
      end do

    end if

    lthermal = .true.
    call rttov_calcbt(chanprof, coef_rttov % coef, lthermal, radiance)
    
  end subroutine tvs_rttovScatt

  !--------------------------------------------------------------------------
  !  tvs_setupPointers
  !--------------------------------------------------------------------------
  subroutine tvs_setupPointers(runObsOperatorWithHydrometeors, sensorIndex, btCount, &
      btCountScatt, hydroChannelsCount, sensorHeaderIndexes, &
      lChannelSubset, obsSpaceData, irBgckMode_opt)
    !
    ! :Purpose: Allocate and initialize tvs_bodyIndexFromBtIndex*  tvs_chanProf*
    !           module variables plus some other local variables.
    !
    implicit none
    
    ! Arguments:
    logical,              intent(in)     :: runObsOperatorWithHydrometeors ! flag to control rttovScatt use in linearized RTTOV
    integer,              intent(in)     :: sensorIndex                    ! sensor Index in NAMTOV namelist section
    integer,              intent(out)    :: btCount                        ! number of BTs simulated using RTTOV
    integer,              intent(out)    :: btCountScatt                   ! number of BTs simulated using RttovScatt
    integer,              intent(out)    :: hydroChannelsCount             ! number of channels simulated using RttovScatt
    integer,              intent(out)    :: sensorHeaderIndexes(:)         ! indexes of radiance observations for the currently processed sensor
    logical, allocatable, intent(out)    :: lChannelSubset(:,:)            ! logical array to setup RttovScatt
    type(struct_obs),     intent(inout)  :: obsSpaceData                   ! obsSpaceData structure
    logical, optional,    intent(in)     :: irBgckMode_opt                 ! background check mode option

    ! Locals:
    integer :: hydroSensorIndex, channelIndex, profileCount
    integer :: btIndex, profileIndex, headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd
    integer :: channelNumber
    logical :: irBgckMode

    if (present(irBgckMode_opt)) then
      irBgckMode = irBgckMode_opt
    else
      irBgckMode = .false.
    end if
    
    if (.not. allocated(tvs_bodyIndexFromBtIndex)) then
      ! These arays are overdimensionned for convenience
      allocate(tvs_bodyIndexFromBtIndex(tvs_maxNumberOfRadiances,tvs_nsensors))
      allocate(tvs_bodyIndexFromBtIndexScatt(tvs_maxNumberOfRadiances,tvs_nsensors))
      allocate(tvs_chanProf(tvs_maxNumberOfRadiances,tvs_nsensors))
      allocate(tvs_chanProfScatt(tvs_maxNumberOfRadiances,tvs_nsensors))
      tvs_bodyIndexFromBtIndex(:,:) = -1
      tvs_bodyIndexFromBtIndexScatt(:,:) = -1
      tvs_chanProf(:,:) % chan = 0
      tvs_chanProf(:,:) % prof = 0
    end if

    hydroSensorIndex = tvs_getHydrometeorsIndex(tvs_instruments(sensorIndex))
    hydroChannelsCount = 0
    if (hydroSensorIndex > 0) then
      do channelIndex = 1, tvs_maxNumberOfChannels
        if (tvs_channelsUsingHydrometeors(hydroSensorIndex,channelIndex) > 0) then
          hydroChannelsCount = hydroChannelsCount + 1
        end if
      end do
      if (hydroChannelsCount == 0) then
        call utl_abort('tvs_setupPointers: you have to initialize channelsUsingHydrometeors(:,:) in NAMTOV namelist section')
      end if
    end if

    profileCount = 0
    do headerIndex = 1, tvs_headerEnd
      ! Currently processed sensor?
      if (tvs_lsensor(headerIndex) == sensorIndex) then
        profileCount = profileCount + 1
        sensorHeaderIndexes(profileCount) = headerIndex
      end if
    end do
    if (profileCount == 0) return

    if (irBgckMode) then
      btCount = profileCount * tvs_nchan(sensorIndex)
      btCountScatt = 0
    else
      btCount = tvs_countRadiances(sensorHeaderIndexes, obsSpaceData)
      if (runObsOperatorWithHydrometeors) then
        btCountScatt = tvs_countRadiancesScatt(sensorHeaderIndexes, obsSpaceData, &
            tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount), sensorIndex)
      else
        btCountScatt = 0
      end if
      btCount = btCount - btCountScatt
    end if
    
    if (tvs_bodyIndexFromBtIndex(1,sensorIndex) == -1) then
      if (irBgckMode) then
        btIndex = 0
        do profileIndex = 1, profileCount
          do channelIndex = 1, tvs_nchan(sensorIndex)
            btIndex = btIndex + 1
            tvs_chanprof(btIndex,sensorIndex) % prof = profileIndex
            tvs_chanprof(btIndex,sensorIndex) % chan = channelIndex
          end do
        end do
        
        do profileIndex = 1, profileCount
          headerIndex = sensorHeaderIndexes(profileIndex)
          if (headerIndex <= 0) cycle
          bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
          bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1
          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            call tvs_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                channelNumber, channelIndex )
            if (channelIndex > 0) then
              tvs_bodyIndexFromBtIndex((profileIndex-1)*tvs_nchan(sensorIndex)+channelIndex,sensorIndex) = bodyIndex
            else
              write(*,*) 'tvs_rttov: strange channel number',channelNumber
            end if
          end do
        end do
      else
        if (btCountScatt > 0) then
          call tvs_getChanprof(sensorHeaderIndexes, obsSpaceData, tvs_chanProf(1:btCount,sensorIndex), &
              iptobs_cma_opt = tvs_bodyIndexFromBtIndex(:,sensorIndex), &
              channelList_opt=tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount), &
              excludeChannelsFromList_opt=.true.)
        else
          call tvs_getChanprof(sensorHeaderIndexes, obsSpaceData, tvs_chanProf(1:btCount,sensorIndex), &
              iptobs_cma_opt = tvs_bodyIndexFromBtIndex(:,sensorIndex))
        end if
      end if
    end if
    if (tvs_bodyIndexFromBtIndexScatt(1,sensorIndex) == -1 .and. btCountScatt > 0) then
      if (allocated(lChannelSubset)) deallocate(lChannelSubset)
      allocate(lChannelSubset(profileCount,tvs_nchan(sensorIndex)))
      call tvs_getChanprof(sensorHeaderIndexes, obsSpaceData, tvs_chanProfScatt(1:btCountScatt,sensorIndex), &
          lchannel_subset_opt = lChannelSubset, iptobs_cma_opt = tvs_bodyIndexFromBtIndexScatt(:,sensorIndex), &
          channelList_opt=tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount))
    end if
    
  end subroutine tvs_setupPointers

  !--------------------------------------------------------------------------
  !  tvs_rttov_tl
  !--------------------------------------------------------------------------
  subroutine tvs_rttov_tl(columnAnlInc, columnTrlOnAnlIncLev, obsSpaceData)
    !
    ! :Purpose: Tangent linear of computation of radiance with rttov_tl
    !   
    implicit none

    ! Arguments:
    type(struct_obs),        intent(inout) :: obsSpaceData         ! obsSpaceData structure
    type(struct_columnData), intent(in)    :: columnAnlInc         ! column structure for pertubation profile
    type(struct_columnData), intent(in)    :: columnTrlOnAnlIncLev ! column structure for background profile

    ! Locals:
    type(struct_vco), pointer :: vco_anl
    integer, allocatable :: sensorHeaderIndexes(:) 
    integer :: allocStatus
    integer :: sensorIndex
    integer :: hydroChannelsCount
    integer :: ilowlvl_M,ilowlvl_T,profileCount,headerIndex,levelIndex,nlv_M,nlv_T
    integer :: profileIndex
    integer :: Vcode
    character(len=4) :: ozoneVarName
    logical, allocatable :: surfTypeIsWater(:)
    real(8), pointer :: delTT(:), delHU(:), delP(:)
    real(8), pointer :: delO3(:)
    real(8), pointer :: delCLW(:)
    real(8), pointer :: delCIW(:), delRF(:), delSF(:)
    integer :: btCount, btcountScatt
    integer :: nthreads
    integer :: btIndex, bodyIndex
    integer :: instrum
    integer :: sensorType   !sensor type(1=infrared; 2=microwave; 3=high resolution, 4=polarimetric)
    integer :: errorStatus
    logical, allocatable :: lChannelSubset(:,:)
    real(8), allocatable :: surfem1(:)
    real(8), allocatable :: surfem1Scatt(:)
    integer, allocatable  :: frequencies(:)
    type(rttov_emissivity), pointer :: emissivity_local(:)
    type(rttov_emissivity), pointer :: emissivity_localScatt(:)
    type(rttov_emissivity), pointer :: emissivity_tl(:)
    type(rttov_emissivity), pointer :: emissivity_tlScatt(:)
    
    type(rttov_radiance) :: radiancedata_d   ! radiances full structure buffer used in rttov calls
    type(rttov_radiance) :: radiancedata_tl  ! tl radiances full structure buffer used in rttov calls
    type(rttov_radiance) :: radiancedata_dScatt   ! radiances full structure buffer used in rttov calls
    type(rttov_radiance) :: radiancedata_tlScatt  ! tl radiances full structure buffer used in rttov calls
    type(rttov_transmission) :: transmission       ! transmission
    type(rttov_transmission) :: transmission_tl    ! transmission tl
    type(rttov_profile), pointer :: profilesdata_tl(:) ! tl profiles buffer used in rttov calls
    type(rttov_profile_cloud), pointer :: cld_profiles_tl(:) !tl profiles buffer used in RttovScatt calls
    logical, pointer :: calcemis(:)
    logical, pointer :: calcemisScatt(:)
    logical :: runObsOperatorWithClw_tl
    logical :: runObsOperatorWithHydrometeors_tl
    real(8) :: obsOMP
    type(rttov_profile), pointer :: profiles(:)
    type(rttov_profile_cloud), pointer :: cld_profiles(:)
         
    if (tvs_headerEnd < 0) return       ! exit if there are not tovs data

    write(*,*) 'tvs_rttov_tl: Starting'

    call tvs_getProfile(profiles, 'tlad', cld_profiles)

    if (tvs_useO3FromTrials .and. .not. col_varExist(columnTrlOnAnlIncLev,'TO3') .and. .not.  col_varExist(columnTrlOnAnlIncLev,'O3L') ) then
      call utl_abort('tvs_rttov_tl: if tvs_useO3FromTrials is set to .true. the ozone variable must be included as an analysis variable in NAMSTATE.')
    else if (tvs_useO3FromTrials) then 
      if (col_varExist(columnTrlOnAnlIncLev,'TO3')) then
        ozoneVarName = 'TO3'
      else
        ozoneVarName = 'O3L'
      end if 
    end if

    !  1.  Set index for model's lowest level and model top

    nlv_M = col_getNumLev(columnTrlOnAnlIncLev,'MM')
    nlv_T = col_getNumLev(columnTrlOnAnlIncLev,'TH')

    if (col_getPressure(columnTrlOnAnlIncLev,1,1,'TH') < col_getPressure(columnTrlOnAnlIncLev,nlv_T,1,'TH')) then
      ilowlvl_M = nlv_M
      ilowlvl_T = nlv_T
    else
      ilowlvl_M = 1
      ilowlvl_T = 1
    end if

    vco_anl => col_getVco(columnTrlOnAnlIncLev)
    Vcode = vco_anl % Vcode
      
    !     1.  Get number of threads available and allocate memory for some variables
    !     .   ---------------------------------------------------------------------- 
    
    ! 2.  Computation of hx for tovs data only
    
    ! Loop over all sensors specified by user

    sensor_loop:  do sensorIndex = 1, tvs_nsensors

      runObsOperatorWithClw_tl = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
                                 tvs_isInstrumUsingCLW(tvs_instruments(sensorIndex)) .and. &
                                 tvs_mwInstrumUsingCLW_tl
      runObsOperatorWithHydrometeors_tl = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
                                          col_varExist(columnTrlOnAnlIncLev,'IWCR') .and. &
                                          tvs_isInstrumUsingHydrometeors(tvs_instruments(sensorIndex)) .and. &
                                          tvs_mwInstrumUsingHydrometeors_tl
      profileCount = tvs_countProfiles(sensorIndex)
      call utl_reallocate(sensorHeaderIndexes, profileCount)
      call tvs_setupPointers(runObsOperatorWithHydrometeors_tl, sensorIndex, btCount, btCountScatt, &
          hydroChannelsCount, sensorHeaderIndexes, lChannelSubset, obsSpaceData)
      
      if (profileCount == 0) cycle sensor_loop
      if (btCount == 0 .and. btCountScatt == 0) cycle  sensor_loop
      
      if (runObsOperatorWithClw_tl) write(*,*) 'tvs_rttov_tl: using clw_data'
      if (runObsOperatorWithHydrometeors_tl) write(*,*) 'tvs_rttov_tl: using hydrometeor data'
      sensorType = tvs_coefs(sensorIndex) % coef % id_sensor
      instrum = tvs_coefs(sensorIndex) % coef % id_inst
      
      allocate(profilesdata_tl(profileCount))
      allocate(cld_profiles_tl(profileCount))
      allocate(surfTypeIsWater(profileCount))
      surfTypeIsWater(:) = .false.

      profileCount = 0
      obs_loop: do headerIndex = 1, tvs_headerEnd  
        if (tvs_lsensor(headerIndex) /= sensorIndex) cycle obs_loop
        profileCount = profileCount + 1
        surfTypeIsWater(profileCount) = ( tvs_ChangedStypValue(obsSpaceData,headerIndex) == surftype_sea )
        sensorHeaderIndexes(profileCount) = headerIndex
      end do obs_loop

      call rttov_alloc_prof(            &
          allocStatus,                  &
          nprofiles=profileCount,       &
          profiles=profilesdata_tl,     &
          nlevels=nlv_T,                &
          opts=tvs_opts(sensorIndex),   &
          asw=1,                        &
          coefs=tvs_coefs(sensorIndex), &
          init=.true.)
      if (allocStatus /= 0) call utl_abort('tvs_rttov_tl: memory allocation error in rttov_alloc_prof')
      call rttov_alloc_scatt_prof (allocStatus,      &
                                   profileCount,     &
                                   cld_profiles_tl,  &
                                   nlv_T,            &
                                   nhydro=5,         &
                                   nhydro_frac=1,    &
                                   asw=1,            &
                                   init=.true.,      &  
                                   flux_conversion=[1,2,0,0,0])
      if (allocStatus /= 0) call utl_abort('tvs_rttov_tl: memory allocation error in rttov_alloc_scatt_prof')
      
      profileLoop:do profileIndex = 1, profileCount
        profilesdata_tl(profileIndex) % gas_units = gas_unit_specconc ! all gas profiles should be provided in kg/kg
        profilesdata_tl(profileIndex) % nlevels   =  nlv_T
        profilesdata_tl(profileIndex) % nlayers   =  nlv_T - 1
        if (tvs_coefs(sensorIndex) % coef % nozone > 0) then
          if (tvs_useO3FromTrials_tl) then
            delO3 => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),trim(ozoneVarName))
            profilesdata_tl(profileIndex) % o3(1:nlv_T) =  delO3(1:nlv_T) * 1.0d-9 ! Assumes model ozone in ug/kg
            profilesdata_tl(profileIndex) % s2m % o  = col_getElem(columnAnlInc,ilowlvl_T,sensorHeaderIndexes(profileIndex),trim(ozoneVarName)) * 1.0d-9 ! Assumes model ozone in ug/kg
          else
            profilesdata_tl(profileIndex) % o3(:) =  0.0d0
          end if
        end if

        ! using the zero CLW value for land FOV
        if (runObsOperatorWithClw_tl) then 
          if (surfTypeIsWater(profileIndex)) then
            delCLW => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'LWCR')
            profilesdata_tl(profileIndex) % clw(1:nlv_T)  = delCLW(:) * tvs_cloudScaleFactor_tl
          else
            profilesdata_tl(profileIndex) % clw(1:nlv_T)  = 0.d0
          end if
        end if

        if (runObsOperatorWithHydrometeors_tl) then 
          if (surfTypeIsWater(profileIndex)) then
            ! rain flux
            if (col_varExist(columnAnlInc,'RF')) then
              delRF => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'RF')
              cld_profiles_tl(profileIndex) % hydro(1:nlv_T,1) = delRF(:) * tvs_cloudScaleFactor_tl
            else
              cld_profiles_tl(profileIndex) % hydro(1:nlv_T,1) = 0.0d0
            end if

            ! snow flux
            if (col_varExist(columnAnlInc,'SF')) then
              delSF => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'SF')
              cld_profiles_tl(profileIndex) % hydro(1:nlv_T,2)  = delSF(:) * tvs_cloudScaleFactor_tl
            else
              cld_profiles_tl(profileIndex) % hydro(1:nlv_T,2) = 0.0d0
            end if

            ! graupel
            cld_profiles_tl(profileIndex) % hydro(1:nlv_T,3)  = 0.d0 ! no information for graupel

            ! cloud liquid water content
            delCLW => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'LWCR')
            cld_profiles_tl(profileIndex) % hydro(1:nlv_T,4) = delCLW(:) * tvs_cloudScaleFactor_tl

            ! cloud ice water content
            delCIW => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'IWCR')
            cld_profiles_tl(profileIndex) % hydro(1:nlv_T,5)  = delCIW(:) * tvs_cloudScaleFactor_tl
          else
            cld_profiles_tl(profileIndex) % hydro(1:nlv_T,1:5)  = 0.d0
          end if ! surfTypeIsWater

          cld_profiles_tl(profileIndex) % hydro_frac(1:nlv_T,1) = 0.d0   ! no perturbation on cloud fraction as it is a binary variable (or or 1.0) in this implementation
        end if ! if (runObsOperatorWithHydrometeors_tl)
        
        profilesdata_tl(profileIndex) % ctp             = 0.0d0
        profilesdata_tl(profileIndex) % cfraction       = 0.0d0
        profilesdata_tl(profileIndex) % zenangle        = 0.0d0
        profilesdata_tl(profileIndex) % azangle         = 0.0d0
        profilesdata_tl(profileIndex) % skin % surftype = 0
        profilesdata_tl(profileIndex) % skin % t        = col_getElem(columnAnlInc,1,sensorHeaderIndexes(profileIndex),'TG')
        profilesdata_tl(profileIndex) % skin % fastem(:)= 0.0d0
        profilesdata_tl(profileIndex) % skin % salinity = 0.0d0
        profilesdata_tl(profileIndex) % s2m % t         = col_getElem(columnAnlInc,ilowlvl_T,sensorHeaderIndexes(profileIndex),'TT')        
        profilesdata_tl(profileIndex) % s2m % q         = 0.d0

        profilesdata_tl(profileIndex) % s2m % p         = col_getElem(columnAnlInc,1,sensorHeaderIndexes(profileIndex),'P0')*MPC_MBAR_PER_PA_R8
        profilesdata_tl(profileIndex) % s2m % u         = col_getElem(columnAnlInc,ilowlvl_M,sensorHeaderIndexes(profileIndex),'UU')
        profilesdata_tl(profileIndex) % s2m % v         = col_getElem(columnAnlInc,ilowlvl_M,sensorHeaderIndexes(profileIndex),'VV')

        delP => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'P_T')
        profilesdata_tl(profileIndex) % p(1:nlv_T)    = delP(:) * MPC_MBAR_PER_PA_R8
        delTT => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'TT')
        profilesdata_tl(profileIndex) % t(1:nlv_T)    = delTT(:)
        delHU => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'HU')
        profilesdata_tl(profileIndex) % q(1:nlv_T)    = delHU(:)
        if (runObsOperatorWithHydrometeors_tl) then
          cld_profiles_tl(profileIndex) % ph (1) = 0.d0
          cld_profiles_tl(profileIndex) % cfrac = 0.d0
          do levelIndex = 1, nlv_T - 1
            cld_profiles_tl(profileIndex) % ph (levelIndex+1) = 0.5d0 * (profilesdata_tl(profileIndex) % p(levelIndex) + profilesdata_tl(profileIndex) % p(levelIndex+1) )
          end do
          cld_profiles_tl(profileIndex) % ph (nlv_T+1) = profilesdata_tl(profileIndex) % s2m % p
        end if
      end do profileLoop

      deallocate(surfTypeIsWater) 

      ! allocate profiledata_tl structures
      if (btCount > 0) then
        call rttov_alloc_tl(                 &
            allocStatus,                     &
            asw=1,                           &
            nprofiles=profileCount,          &
            nchanprof=btCount,               &
            nlevels=nlv_T,                   &
            opts=tvs_opts(sensorIndex),      &
            coefs=tvs_coefs(sensorIndex),    &
            transmission=transmission,       &
            transmission_tl=transmission_tl, &
            radiance=radiancedata_d,         &
            radiance_tl=radiancedata_tl,     &
            calcemis=calcemis,               &
            emissivity=emissivity_local,     &
            emissivity_tl=emissivity_tl,     &
            init=.true.)
        !   Prepare all input variables required by rttov.
        if (allocStatus /= 0) call utl_abort('tvs_rttov_tl: memory allocation error 1 in rttov_alloc_tl')
        allocate(surfem1(btCount))
         !    get Hyperspecral IR emissivities
        if (tvs_isInstrumHyperSpectral(instrum)) call tvs_getHIREmissivities(sensorHeaderIndexes, &
            obsSpaceData, surfem1)
       
        call tvs_getOtherEmissivities(tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes, sensorType, instrum, surfem1, calcemis)
        
        if (sensorType == sensor_id_mw) then
          if (col_varExist(columnAnlInc, 'EMMW')) then
            ! Read surface emissivity from column when it's included as an analysis variable
            ! Set the default surface emissivity values
            emissivity_local(:) % emis_in = surfem1(:)

            ! Setup the emissivity_tl from column object
            call sse_setupEmissivityfromState(emissivity_local, obsSpaceData, tvs_bodyIndexFromBtIndex(:,sensorIndex), tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes, &
                                              tvs_nsensors, tvs_lsensor, tvs_instrumentName, &
                                              tvs_maxChannelNumber, tvs_channelOffset, tvs_ichan, profiles(:) % skin % surftype, &
                                              emissivityProfDt_opt = tvs_emissivityFromTrl)
          else if (tvs_useSfcEmissObsSpace) then
            ! Set the default surface emissivity values
            emissivity_local(:) % emis_in = surfem1(:)
  
            ! Setup the surface emissvity from obsSpaceData Object 
            call sse_emissFromObsSpace(obsSpaceData, emissivity_local, tvs_bodyIndexFromBtIndex(:,sensorIndex), tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes)    
          else
            ! Read surface emissivity from emissivity atlas
            call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorIndex, tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes)
          end if
        else
          emissivity_local(:) % emis_in = surfem1(:)
        end if
  
        !  2.3  Compute tl radiance with rttov_tl
        
        if (col_varExist(columnAnlInc, 'EMMW') .and. sensorType == sensor_id_mw) then
          call sse_setupEmissivityfromState(emissivity_tl, obsSpaceData, tvs_bodyIndexFromBtIndex(:,sensorIndex), tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes, &
                                  tvs_nsensors, tvs_lsensor, tvs_instrumentName, tvs_maxChannelNumber, tvs_channelOffset, tvs_ichan, profiles(:) % skin % surftype, columProfTl_opt  = columnAnlInc)
        else
          emissivity_tl(:) % emis_in = 0.0d0
        end if

        errorStatus = errorStatus_success

        !  set nthreads to actual number of threads which will be used.

        nthreads = min(mmpi_numThread, profileCount)  
        call rttov_parallel_tl(                             &
            errorStatus,                                    & ! out
            tvs_chanProf(1:btCount,sensorIndex),            & ! in
            tvs_opts(sensorIndex),                          & ! in
            profiles(sensorHeaderIndexes(:)),               & ! in
            profilesdata_tl,                                & ! inout
            tvs_coefs(sensorIndex),                         & ! in
            transmission,                                   & ! inout
            transmission_tl,                                & ! inout
            radiancedata_d,                                 & ! inout
            radiancedata_tl,                                & ! inout
            calcemis=calcemis,                              & ! in
            emissivity=emissivity_local,                    & ! in
            emissivity_tl=emissivity_tl,                    & ! inout
            nthreads=nthreads )                               ! in

        if (errorStatus /= errorStatus_success) then
          write(*,*) 'Error in rttov_parallel_tl', errorStatus
          write(*,*) 'temperature           profile=',profiles(sensorHeaderIndexes(1)) % t(:)
          write(*,*) 'temperature increment profile=',profilesdata_tl(1) % t(:)
          call utl_abort('tvs_rttov_tl')
        end if

        !  2.4  Store hx in obsSpaceData,OBS_WORK
      
        do btIndex = 1, btCount
          bodyIndex = tvs_bodyIndexFromBtIndex(btIndex,sensorIndex)
          call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex, &
              radiancedata_tl % bt(btIndex) )
          if (tvs_debug) then
            obsOMP = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
            write(*,'(a,i4,2f8.2)') ' ichn,sim,obs= ', &
                tvs_chanProf(btIndex,sensorIndex) % chan, radiancedata_tl % bt(btIndex), obsOMP
          end if
        end do
        call rttov_alloc_tl(                 &
           allocStatus,                      &
           asw=0,                            &
           nprofiles=profileCount,           &
           nchanprof=btCount,                &
           nlevels=nlv_T,                    &
           opts=tvs_opts(sensorIndex),       &
           coefs=tvs_coefs(sensorIndex),     &
           transmission=transmission,        &
           transmission_tl=transmission_tl,  &
           radiance=radiancedata_d,          &
           radiance_tl=radiancedata_tl,      &
           calcemis=calcemis,                &
           emissivity=emissivity_local,      &
           emissivity_tl=emissivity_tl )
        if (allocStatus /= 0) call utl_abort('tvs_rttov_tl: memory deallocation 1 error in rttov_alloc_tl')
        deallocate(surfem1)
 
      end if ! if (btCount > 0)
      
      if (btCountScatt > 0) then 
        call rttov_alloc_tl(                  &
            allocStatus,                      &  
            asw=1,                            &
            nprofiles=profileCount,           &
            nchanprof=btCountScatt,           &
            nlevels=nlv_T,                    &
            opts=tvs_opts(sensorIndex),       &
            coefs=tvs_coefs(sensorIndex),     &
            radiance=radiancedata_dScatt,     &
            radiance_tl=radiancedata_tlScatt, &
            calcemis=calcemisScatt,           &
            emissivity=emissivity_localScatt, &
            emissivity_tl=emissivity_tlScatt, &
            init=.true.)
        if (allocStatus /= 0) call utl_abort('tvs_rttov_tl: memory allocation error 2 in rttov_alloc_tl')
        
        ! Prepare all input variables required by rttovScatt.
     
        allocate(surfem1Scatt(btCountScatt))
        allocate(frequencies(btCountScatt))
        call rttov_scatt_setupindex(                          &
            errorStatus,                                      &
            profileCount,                                     & ! number of profiles
            tvs_nchan(sensorIndex),                           & ! number of channels 
            tvs_coefs(sensorIndex),                           & ! coef structure read in from rttov coef file
            tvs_coef_scatt(sensorIndex),                      & ! coef structure read in from rttov coef file
            btCountScatt,                                     & ! number of calculated channels
            tvs_chanProfScatt(1:btCountScatt,sensorIndex),    & ! channels and profile numbers
            frequencies,                                      & ! array, frequency number for each channel
            lChannelSubset )                                    ! OPTIONAL array of logical flags to indicate a subset of channels
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'tvs_rttov_tl: fatal error in rttov_scatt_setupindex ', errorStatus
          call utl_abort('tvs_rttov_tl')
        end if
      
        call tvs_getOtherEmissivities(tvs_chanProfScatt(1:btCountScatt,sensorIndex), sensorHeaderIndexes, sensorType, instrum, surfem1Scatt, calcemisScatt)

        call tvs_getMWemissivityFromAtlas(surfem1Scatt(1:btcountScatt), emissivity_localScatt, sensorIndex, &
            tvs_chanProfScatt(1:btCountScatt,sensorIndex), sensorHeaderIndexes)
        errorStatus = errorStatus_success
        emissivity_tlScatt(:) % emis_in = 0.0d0
        call rttov_scatt_tl(                                  &
            errorStatus,                                      & ! out
            tvs_opts_scatt(sensorIndex),                      & ! in
            nlv_T,                                            & ! in
            tvs_chanProfScatt(1:btCountScatt,sensorIndex),    & ! in
            frequencies,                                      & ! in
            profiles(sensorHeaderIndexes(:)),                 & ! in  
            cld_profiles(sensorHeaderIndexes(:)),             & ! in
            tvs_coefs(sensorIndex),                           & ! in
            tvs_coef_scatt(sensorIndex),                      & ! in
            calcemisScatt,                                    & ! in
            emissivity_localScatt,                            & ! inout
            profilesdata_tl,                                  & ! in
            cld_profiles_tl,                                  & ! in
            emissivity_tlScatt,                               & ! inout
            radiancedata_dScatt,                              & ! inout
            radiancedata_tlScatt)                               ! inout
        
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'Error in rttov_scatt_tl', errorStatus
          write(*,*) 'temperature           profile=',profiles(sensorHeaderIndexes(1)) % t(:)
          write(*,*) 'temperature increment profile=',profilesdata_tl(1) % t(:)
          call utl_abort('tvs_rttov_tl')
        end if

        !  2.4  Store hx in obsSpaceData,OBS_WORK
      
        do btIndex = 1, btCountScatt
          bodyIndex = tvs_bodyIndexFromBtIndexScatt(btIndex,sensorIndex)
          call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex, &
              radiancedata_tlScatt % bt(btIndex) )
          if (tvs_debug) then
            obsOMP = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
            write(*,'(a,i4,2f8.2)') ' ichn,sim,obs= ', &
                tvs_chanprofScatt(btIndex,sensorIndex) % chan, radiancedata_tlScatt % bt(btIndex), obsOMP
          end if
        end do
        deallocate(surfem1Scatt)
        deallocate(frequencies)
        call rttov_alloc_tl(                  &
            allocStatus,                      &
            asw=0,                            &
            nprofiles=profileCount,           &
            nchanprof=btCountScatt,           &
            nlevels=nlv_T,                    &
            opts=tvs_opts(sensorIndex),       &
            coefs=tvs_coefs(sensorIndex),     &
            radiance=radiancedata_dScatt,     &
            radiance_tl=radiancedata_tlScatt, &
            calcemis=calcemisScatt,           &
            emissivity=emissivity_localScatt, &
            emissivity_tl=emissivity_tlScatt )
        if (allocStatus /= 0) call utl_abort('tvs_rttov_tl: memory deallocation error 2 in rttov_alloc_tl')
      end if  ! if (btCountScatt > 0)
      
      call rttov_alloc_scatt_prof (allocStatus,                 &
                                   profileCount,                &
                                   cld_profiles_tl,             &
                                   nlv_T,                       &
                                   nhydro=5,                    &
                                   nhydro_frac=1,               &
                                   asw=0,                       &   
                                   flux_conversion=[1,2,0,0,0])
      if (allocStatus /= 0) call utl_abort('tvs_rttov_tl: memory deallocation error in rttov_alloc_scatt_prof')
      deallocate(cld_profiles_tl)
      
      call rttov_alloc_prof(            &
          allocStatus,                  &
          nprofiles=profileCount,       &
          profiles=profilesdata_tl,     &
          nlevels=nlv_T,                &
          opts=tvs_opts(sensorIndex),   &
          asw=0,                        &
          coefs=tvs_coefs(sensorIndex), &
          init=.true.)
      if (allocStatus /= 0) call utl_abort('tvs_rttov_tl: memory deallocation error in rttov_alloc_prof')
      deallocate (sensorHeaderIndexes)
    end do sensor_loop

    nullify( profiles )
    if (allocated(tvs_bodyIndexFromBtIndex)) then
      deallocate( tvs_bodyIndexFromBtIndex )
      deallocate( tvs_bodyIndexFromBtIndexScatt )
      deallocate( tvs_chanProf )
      deallocate( tvs_chanProfScatt )
    end if
    write(*,*) 'tvs_rttov_tl: Finished'

  end subroutine tvs_rttov_tl

  !--------------------------------------------------------------------------
  !  tvs_rttov_ad
  !--------------------------------------------------------------------------
  subroutine tvs_rttov_ad( columnAnlInc, columnTrlOnAnlIncLev, obsSpaceData )
    !
    ! :Purpose: Adjoint of computation of radiance with rttov_ad
    !

    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout)    :: columnAnlInc         ! column structure for pertubation profile 
    type(struct_columnData), intent(in)       :: columnTrlOnAnlIncLev ! column structure for background profile
    type(struct_obs),        intent(inout)    :: obsSpaceData         ! obsSpaceData structure

    ! Locals:
    type(struct_vco), pointer :: vco_anl
    integer, allocatable :: sensorHeaderIndexes(:) 
    integer :: allocStatus
    integer :: nthreads
    integer :: sensorIndex
    integer :: hydroChannelsCount
    integer :: ilowlvl_T,ilowlvl_M,profileCount,headerIndex,nlv_M,nlv_T
    integer :: profileIndex, levelIndex
    integer :: Vcode
    real(8), allocatable :: tt_ad(:,:)
    real(8), allocatable :: hu_ad(:,:)
    real(8), allocatable :: pressure_ad(:,:)  
    real(8), allocatable :: ozone_ad(:,:)
    character(len=4) :: ozoneVarName
    real(8), allocatable :: clw_ad(:,:),clwScatt_ad(:,:)
    real(8), allocatable :: ciw_ad(:,:), rf_ad(:,:), sf_ad(:,:)
    logical, allocatable :: surfTypeIsWater(:)
    logical, allocatable :: lChannelSubset(:,:)
    real(8), pointer :: uu_column(:),vv_column(:),tt_column(:),hu_column(:),ps_column(:)
    real(8), pointer :: tg_column(:),p_column(:),o3_column(:),clw_column(:)
    real(8), pointer :: ciw_column(:), rf_column(:),sf_column(:)
    integer :: btCount, btCountScatt
    integer :: instrum
    integer :: btIndex, bodyIndex
    integer :: sensorType   ! sensor type (1=infrared; 2=microwave; 3=high resolution, 4=polarimetric)  
    integer :: errorStatus
    real(8), allocatable :: surfem1(:)
    real(8), allocatable :: surfem1Scatt(:)
    integer, allocatable :: frequencies(:)
    type(rttov_emissivity), pointer :: emissivity_local(:)
    type(rttov_emissivity), pointer :: emissivity_ad(:)
    type(rttov_emissivity), pointer :: emissivity_localScatt(:)
    type(rttov_emissivity), pointer :: emissivity_adScatt(:)
    type(rttov_transmission) :: transmission,transmission_ad
    type(rttov_radiance) :: radiancedata_ad, radiancedata_d
    type(rttov_radiance) :: radiancedata_adScatt, radiancedata_dScatt
    type(rttov_profile), pointer  :: profilesdata_ad(:) ! ad profiles buffer used in rttov calls
    type(rttov_profile), pointer  :: profiles(:)
    type(rttov_profile_cloud), pointer  :: cld_profiles(:)
    type(rttov_profile_cloud), pointer  :: cld_profiles_ad(:)
    logical, pointer :: calcemis(:)
    logical, pointer :: calcemisScatt(:)
    logical :: runObsOperatorWithClw_ad
    logical :: runObsOperatorWithHydrometeors_ad
         
    if (tvs_headerEnd < 0) return      ! exit if there are not tovs data
    write(*,*) 'tvs_rttov_ad: Starting'

    call tvs_getProfile(profiles, 'tlad', cld_profiles)

    if (tvs_useO3FromTrials .and. .not. col_varExist(columnTrlOnAnlIncLev,'TO3') .and. .not.  col_varExist(columnTrlOnAnlIncLev,'O3L') ) then
      call utl_abort('tvs_rttov_ad: if tvs_useO3FromTrials is set to .true. the ozone variable must be included as an analysis variable in NAMSTATE.')
    else if (tvs_useO3FromTrials) then 
      if (col_varExist(columnTrlOnAnlIncLev,'TO3')) then
        ozoneVarName = 'TO3'
      else
        ozoneVarName = 'O3L'
      end if
    end if

    !     1.    Set index for model's lowest level and model top

    nlv_M = col_getNumLev(columnTrlOnAnlIncLev,'MM')
    nlv_T = col_getNumLev(columnTrlOnAnlIncLev,'TH')

    if (col_getPressure(columnTrlOnAnlIncLev,1,1,'TH') < col_getPressure(columnTrlOnAnlIncLev,nlv_T,1,'TH')) then
      ilowlvl_M = nlv_M
      ilowlvl_T = nlv_T
    else
      ilowlvl_M = 1
      ilowlvl_T = 1
    end if

    vco_anl => col_getVco(columnTrlOnAnlIncLev)
    Vcode = vco_anl % Vcode

    !     1.  Get number of threads available and allocate memory for some variables
 
  

    !     2.  Computation of adjoint hx for tovs data only

    ! Loop over all sensors specified by user

    sensor_loop:do sensorIndex = 1, tvs_nsensors
      
      runObsOperatorWithClw_ad = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
                                 tvs_isInstrumUsingCLW(tvs_instruments(sensorIndex)) .and. &
                                 tvs_mwInstrumUsingCLW_tl
      
      runObsOperatorWithHydrometeors_ad = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
                                          col_varExist(columnTrlOnAnlIncLev,'IWCR') .and. &
                                          tvs_isInstrumUsingHydrometeors(tvs_instruments(sensorIndex)) .and. &
                                          tvs_mwInstrumUsingHydrometeors_tl
      profileCount = tvs_countProfiles(sensorIndex)
      call utl_reallocate(sensorHeaderIndexes, profileCount)
      call tvs_setupPointers(runObsOperatorWithHydrometeors_ad, sensorIndex, btCount, btCountScatt, &
          hydroChannelsCount, sensorHeaderIndexes, lChannelSubset, obsSpaceData)
      
      if (profileCount == 0) cycle sensor_loop
      if (btCount == 0 .and. btCountScatt == 0) cycle sensor_loop
      
      sensorType = tvs_coefs(sensorIndex) % coef % id_sensor
      instrum = tvs_coefs(sensorIndex) % coef % id_inst
     
      allocate(tt_ad(nlv_T,profileCount))
      allocate(hu_ad(nlv_T,profileCount))
      allocate(pressure_ad(nlv_T,profileCount))
      if (tvs_useO3FromTrials_tl) then
        if (tvs_coefs(sensorIndex) % coef % nozone > 0) then
          allocate(ozone_ad(nlv_T,profileCount))
        end if
      end if
      if (runObsOperatorWithClw_ad) allocate(clw_ad(nlv_T,profileCount))
      if (runObsOperatorWithHydrometeors_ad) then
        allocate(clwScatt_ad(nlv_T,profileCount))
        allocate(ciw_ad(nlv_T,profileCount))
        allocate(rf_ad(nlv_T,profileCount))
        allocate(sf_ad(nlv_T,profileCount))
      end if
      allocate(surfTypeIsWater(profileCount))
      surfTypeIsWater(:) = .false.

      !  2.1  Calculate the actual number of threads which will be used.

      nthreads = min(mmpi_numThread, profileCount )  
      allocate(profilesdata_ad(profileCount))
      call rttov_alloc_prof(            &
          allocStatus,                  &
          nprofiles=profileCount,       &
          profiles=profilesdata_ad,     &
          nlevels=nlv_T,                &
          opts=tvs_opts(sensorIndex),   &
          asw=1,                        &
          coefs=tvs_coefs(sensorIndex), &
          init=.true.)
      if (allocStatus /= 0) call utl_abort('tvs_rttov_ad: memory allocation error in rttov_alloc_prof')
      !  2.2  Prepare all input variables required by rttov_ad.
      if (btCount > 0) then
        call rttov_alloc_ad(                 &
            allocStatus,                     &
            asw=1,                           &
            nprofiles=profileCount,          &
            nchanprof=btCount,               &
            nlevels=nlv_T,                   &
            opts=tvs_opts(sensorIndex),      &
            coefs=tvs_coefs(sensorIndex),    &
            transmission=transmission,       &
            transmission_ad=transmission_ad, &
            radiance=radiancedata_d,         &
            radiance_ad=radiancedata_ad,     &
            calcemis=calcemis,               &
            emissivity=emissivity_local,     &
            emissivity_ad=emissivity_ad,     &
            init=.true.)
        if (allocStatus /= 0) call utl_abort('tvs_rttov_ad: memory allocation error 1 in rttov_alloc_ad')
        allocate(surfem1(btCount))
      
        !  get Hyperspectral IR emissivities
        if (tvs_isInstrumHyperSpectral(instrum)) call tvs_getHIREmissivities(sensorHeaderIndexes, obsSpaceData, surfem1)
        
        !     get non Hyperspectral IR emissivities
        call tvs_getOtherEmissivities(tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes, sensorType, instrum, surfem1, calcemis)

        if (sensorType == sensor_id_mw) then
          if (col_varExist(columnAnlInc, 'EMMW')) then
            ! Read surface emissivity from column when it's included as an analysis variable

            ! Set the default surface emissivity values
            emissivity_local(:) % emis_in = surfem1(:)

            ! Setup the surface emissvity from column object to rttov emissivity_local
            call sse_setupEmissivityfromState(emissivity_local, obsSpaceData, tvs_bodyIndexFromBtIndex(:,sensorIndex), tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes, &
                                        tvs_nsensors, tvs_lsensor, tvs_instrumentName, &
                                        tvs_maxChannelNumber, tvs_channelOffset, tvs_ichan, profiles(:) % skin % surftype, &
                                        emissivityProfDt_opt = tvs_emissivityFromTrl)
          else if (tvs_useSfcEmissObsSpace) then
            ! Set the default surface emissivity values
            emissivity_local(:) % emis_in = surfem1(:)

            ! Setup the surface emissvity from obsSpaceData Object 
            call sse_emissFromObsSpace(obsSpaceData, emissivity_local, tvs_bodyIndexFromBtIndex(:,sensorIndex), tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes)    
          else
            ! Read surface emissivity from emissivity atlas
            call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorIndex, tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes)
          end if
        else
          emissivity_local(:) % emis_in = surfem1(:)
        end if
        
        do btIndex = 1, btCount
          bodyIndex = tvs_bodyIndexFromBtIndex(btIndex,sensorIndex)
          radiancedata_ad % bt( btIndex ) = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
        end do
      
        errorStatus = errorStatus_success
        emissivity_ad(:) % emis_in = 0.0d0
        emissivity_ad(:) % emis_out = 0.0d0
        call rttov_parallel_ad(                             &
            errorstatus,                                    &! out
            tvs_chanProf(1:btCount,sensorIndex),            &! in
            tvs_opts(sensorIndex),                          &! in
            profiles(sensorHeaderIndexes(:)),               &! in
            profilesdata_ad,                                &! in
            tvs_coefs(sensorIndex),                         &! in
            transmission,                                   &! inout
            transmission_ad,                                &! inout
            radiancedata_d,                                 &! inout
            radiancedata_ad,                                &! inout
            calcemis=calcemis,                              &! in
            emissivity=emissivity_local,                    &! inout
            emissivity_ad=emissivity_ad,                    &! inout
            nthreads = nthreads )
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'Error in rttov_parallel_ad', errorStatus
          call utl_abort('tvs_rttov_ad')
        end if

        call rttov_alloc_ad(                 &
            allocStatus,                     &
            asw=0,                           &
            nprofiles=profileCount,          &
            nchanprof=btCount,               &
            nlevels=nlv_T,                   &
            opts=tvs_opts(sensorIndex),      &
            coefs=tvs_coefs(sensorIndex),    &
            transmission=transmission,       &
            transmission_ad=transmission_ad, &
            radiance=radiancedata_d,         &
            radiance_ad=radiancedata_ad,     &
            calcemis=calcemis,               &
            emissivity=emissivity_local)
        if (allocStatus /= 0) call utl_abort('tvs_rttov_ad: memory deallocation error 1 in rttov_alloc_ad')
        deallocate(surfem1)
      end if ! if (btCount > 0)

      if (btCountScatt > 0) then
        call rttov_alloc_ad(                   &
            allocStatus,                       &
            asw=1,                             &
            nprofiles=profileCount,            &
            nchanprof=btCountScatt,            &
            nlevels=nlv_T,                     &
            opts=tvs_opts(sensorIndex),        &
            coefs=tvs_coefs(sensorIndex),      &
            radiance=radiancedata_dScatt,      &
            radiance_ad=radiancedata_adScatt,  &
            calcemis=calcemisScatt,            &
            emissivity=emissivity_localScatt,  &
            emissivity_ad=emissivity_adScatt,  &
            init=.true.)
        if (allocStatus /= 0) call utl_abort('tvs_rttov_ad: memory allocation error 2 in rttov_alloc_ad')
        allocate(surfem1Scatt(btCountScatt))
        allocate(frequencies(btCountScatt))
        allocate(cld_profiles_ad(profileCount))
        call rttov_alloc_scatt_prof(allocStatus,      &
                                    profileCount,     &
                                    cld_profiles_ad,  &
                                    nlv_T,            &
                                    nhydro=5,         &
                                    nhydro_frac=1,    &
                                    asw=1,            &
                                    init=.true.,      &
                                    flux_conversion=[1,2,0,0,0])
        if (allocStatus /= 0) call utl_abort('tvs_rttov_ad: memory allocation error in rttov_alloc_scatt_prof')
        ! Build the list of channels/profiles indices
        call rttov_scatt_setupindex(                          &
            errorStatus,                                      &
            profileCount,                                     &  ! number of profiles
            tvs_nchan(sensorIndex),                           &  ! number of channels 
            tvs_coefs(sensorIndex),                           &  ! coef structure read in from rttov coef file
            tvs_coef_scatt(sensorIndex),                      &  ! coef structure read in from rttov coef file
            btcountScatt,                                     &  ! number of calculated channels
            tvs_chanProfScatt(1:btCountScatt,sensorIndex),    &  ! channels and profile numbers
            frequencies,                                      &  ! array, frequency number for each channel
            lChannelSubset)                                      ! OPTIONAL array of logical flags to indicate a subset of channels
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'tvs_rttov_ad: fatal error in rttov_scatt_setupindex ', errorStatus
          call utl_abort('tvs_rttov_ad')
        end if
        !     get non Hyperspectral IR emissivities
        call tvs_getOtherEmissivities(tvs_chanProfScatt(1:btCountScatt,sensorIndex), sensorHeaderIndexes, &
            sensorType, instrum, surfem1Scatt, calcemisScatt)

        call tvs_getMWemissivityFromAtlas(surfem1Scatt(1:btcountScatt), emissivity_localScatt, sensorIndex, &
            tvs_chanProfScatt(1:btCountScatt,sensorIndex), sensorHeaderIndexes)
        
        do btIndex = 1, btCountScatt
          bodyIndex = tvs_bodyIndexFromBtIndexScatt(btIndex,sensorIndex)
          radiancedata_adScatt % bt( btIndex ) = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
        end do
      
        errorStatus = errorStatus_success
        emissivity_adScatt(:) % emis_in = 0.0d0
        emissivity_adScatt(:) % emis_out = 0.0d0
        !  2.3  Compute ad radiance with rttov_ad
        call rttov_scatt_ad(                                  & 
            errorStatus,                                      &! out
            tvs_opts_scatt(sensorIndex),                      &! in
            nlv_T,                                            &! in
            tvs_chanProfScatt(1:btCountScatt,sensorIndex),    &! in
            frequencies,                                      &! in
            profiles(sensorHeaderIndexes(:)),                 &! in
            cld_profiles(sensorHeaderIndexes(:)),             &! in
            tvs_coefs(sensorIndex),                           &! in
            tvs_coef_scatt(sensorIndex),                      &! in
            calcemisScatt,                                    &! in
            emissivity_localScatt,                            &! inout
            profilesdata_ad,                                  &! inout
            cld_profiles_ad,                                  &! inout
            emissivity_adScatt,                               &! inout
            radiancedata_dScatt,                              &! inout
            radiancedata_adScatt)                              ! inout
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'Error in rttov_scatt_ad', errorStatus
          call utl_abort('tvs_rttov_ad')
        end if
      
        call rttov_alloc_ad(                  &
            allocStatus,                      &
            asw=0,                            &
            nprofiles=profileCount,           &
            nchanprof=btCountScatt,           &
            nlevels=nlv_T,                    &
            opts=tvs_opts(sensorIndex),       &
            coefs=tvs_coefs(sensorIndex),     &
            radiance=radiancedata_dScatt,     &
            radiance_ad=radiancedata_adScatt, &
            calcemis=calcemisScatt,           &
            emissivity=emissivity_localScatt, &
            emissivity_ad=emissivity_adScatt )
        if (allocStatus /= 0) call utl_abort('tvs_rttov_ad: memory deallocation error 2 in rttov_alloc_ad')
        deallocate(surfem1Scatt)
        deallocate(frequencies)
      end if !if (btCountScatt > 0)

      !   2.0  Store adjoints in columnData object
      tt_ad(:,:) = 0.d0
      hu_ad(:,:) = 0.d0
      pressure_ad(:,:) = 0.d0
      if (tvs_useO3FromTrials_tl) then
        if (tvs_coefs(sensorIndex) % coef % nozone > 0) ozone_ad(:,:) = 0.d0
      end if
      if (runObsOperatorWithClw_ad) clw_ad(:,:) = 0.d0
      if (runObsOperatorWithHydrometeors_ad) then
        clwScatt_ad(:,:) = 0.d0
        ciw_ad(:,:) = 0.d0
        rf_ad(:,:) = 0.d0
        sf_ad(:,:) = 0.d0
      end if
      
      do profileIndex= 1, profileCount
        headerIndex = sensorHeaderIndexes(profileIndex)
        ps_column => col_getColumn(columnAnlInc,headerIndex,'P0')
        p_column  => col_getColumn(columnAnlInc,headerIndex,'P_T')
        tg_column => col_getColumn(columnAnlInc,headerIndex,'TG')
        tt_column => col_getColumn(columnAnlInc,headerIndex,'TT')
        hu_column => col_getColumn(columnAnlInc,headerIndex,'HU')
        uu_column => col_getColumn(columnAnlInc,headerIndex,'UU')
        vv_column => col_getColumn(columnAnlInc,headerIndex,'VV')
      
        tt_ad(:,profileIndex) = profilesdata_ad(profileIndex) % t(:)
        hu_ad(:,profileIndex) = profilesdata_ad(profileIndex) % q(:)
        pressure_ad(:,profileIndex) =  profilesdata_ad(profileIndex) % p(:)
        tg_column(1) = profilesdata_ad(profileIndex) % skin % t
        tt_column(ilowlvl_T) = profilesdata_ad(profileIndex) % s2m % t
        ps_column(1) = profilesdata_ad(profileIndex) % s2m % p * MPC_MBAR_PER_PA_R8
        hu_column(ilowlvl_T) = 0.d0 
        uu_column(ilowlvl_M) = profilesdata_ad(profileIndex) % s2m % u
        vv_column(ilowlvl_M) = profilesdata_ad(profileIndex) % s2m % v

        if (tvs_useO3FromTrials_tl) then
          if (tvs_coefs(sensorIndex) % coef % nozone > 0) then
            ozone_ad(:,profileIndex) = profilesdata_ad(profileIndex) % o3(:)
          end if
        end if
        !end of the block of code to be removed later
        
        if (runObsOperatorWithClw_ad) then
          clw_ad(:,profileIndex) = profilesdata_ad(profileIndex) % clw(:)
        end if
      
        if (runObsOperatorWithHydrometeors_ad .and. btCountScatt > 0) then
          rf_ad(:,profileIndex)  = cld_profiles_ad(profileIndex) % hydro(:,1)
          sf_ad(:,profileIndex)  = cld_profiles_ad(profileIndex) % hydro(:,2)
          clwScatt_ad(:,profileIndex) = cld_profiles_ad(profileIndex) % hydro(:,4)
          ciw_ad(:,profileIndex) = cld_profiles_ad(profileIndex) % hydro(:,5)
        end if
      end do

      ! Store surface emissivity adjoint into column object
      if (col_varExist(columnAnlInc, 'EMMW') .and. sensorType == sensor_id_mw) then
        ! Setup emissivity in column object from emissivity_ad
        call sse_setupEmissivityfromState(emissivity_ad, obsSpaceData, tvs_bodyIndexFromBtIndex(:,sensorIndex), & 
                                          tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes, &
                                          tvs_nsensors, tvs_lsensor, tvs_instrumentName, &
                                          tvs_maxChannelNumber, tvs_channelOffset, tvs_ichan, profiles(:) % skin % surftype, &
                                          columProfAd_opt = columnAnlInc)
      end if
    
      call rttov_alloc_ad(                   &
            allocStatus,                     &
            asw=0,                           &
            nprofiles=profileCount,          &
            nchanprof=btCount,               &
            nlevels=nlv_T,                   &
            opts=tvs_opts(sensorIndex),      &
            coefs=tvs_coefs(sensorIndex),    &
            emissivity_ad=emissivity_ad )
      if (allocStatus /= 0) call utl_abort('tvs_rttov_ad: memory deallocation error 1 in rttov_alloc_ad')

      call rttov_alloc_prof(            &
          allocStatus,                  &
          nprofiles=profileCount,       &
          profiles=profilesdata_ad,     &
          nlevels=nlv_T,                &
          opts=tvs_opts(sensorIndex),   &
          asw=0,                        &
          coefs=tvs_coefs(sensorIndex), &
          init=.true.)
      if (allocStatus /= 0) call utl_abort('tvs_rttov_ad: memory deallocation error in rttov_alloc_prof')
      deallocate(profilesdata_ad)
    
      if (btCountScatt > 0) then
        call rttov_alloc_scatt_prof(  &
            allocStatus,              &
            profileCount,             &
            cld_profiles_ad,          &
            nlv_T,                    &
            nhydro=5,                 &
            nhydro_frac=1,            &
            asw=0,                    &     
            flux_conversion=[1,2,0,0,0])
        if (allocStatus /= 0) call utl_abort('tvs_rttov_ad: memory deallocation error in rttov_alloc_scatt_prof')
        deallocate(cld_profiles_ad)
      end if

      !     .  2.1  Store adjoints in columnData object
      !     .       -----------------------------------

      do profileIndex = 1, profileCount 
        p_column  => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'P_T')
        tt_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'TT')
        hu_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'HU')
      
        do levelIndex = 1, nlv_T
          p_column(levelIndex) = p_column(levelIndex) + pressure_ad(levelIndex,profileIndex) * MPC_MBAR_PER_PA_R8
          tt_column(levelIndex) = tt_column(levelIndex) + tt_ad(levelIndex,profileIndex)
          hu_column(levelIndex) = hu_column(levelIndex) + hu_ad(levelIndex,profileIndex)
        end do
      end do

      if (tvs_useO3FromTrials_tl) then
        if (tvs_coefs(sensorIndex) % coef % nozone > 0) then
          do  profileIndex = 1, profileCount 
            o3_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex),trim(ozoneVarName))
            do levelIndex = 1, col_getNumLev(columnAnlInc,'TH')
              o3_column(levelIndex) = o3_column(levelIndex) + ozone_ad(levelIndex,profileIndex) * 1.0d-9
            end do
          end do
        end if
      end if

      if (runObsOperatorWithClw_ad) then
        do  profileIndex = 1 , profileCount 
          surfTypeIsWater(profileIndex) = (tvs_ChangedStypValue(obsSpaceData,sensorHeaderIndexes(profileIndex)) == surftype_sea)
          if (surfTypeIsWater(profileIndex)) then
            clw_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'LWCR')
            do levelIndex = 1, col_getNumLev(columnAnlInc, 'TH')
              clw_column(levelIndex) = clw_column(levelIndex) + &
                  clw_ad(levelIndex,profileIndex) * tvs_cloudScaleFactor_tl
            end do
          end if
        end do
      end if
    
      if (runObsOperatorWithHydrometeors_ad) then
        do  profileIndex = 1 , profileCount
          surfTypeIsWater(profileIndex) = (tvs_ChangedStypValue(obsSpaceData,sensorHeaderIndexes(profileIndex)) == surftype_sea)
          if (surfTypeIsWater(profileIndex)) then
            ! rain flux
            if (col_varExist(columnAnlInc,'RF')) then
              rf_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'RF')
              do levelIndex = 1, col_getNumLev(columnAnlInc, 'TH')
                rf_column(levelIndex) = rf_column(levelIndex) + &
                  rf_ad(levelIndex,profileIndex) * tvs_cloudScaleFactor_tl
              end do
            end if

            ! snow flux
            if (col_varExist(columnAnlInc,'SF')) then
              sf_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'SF')
              do levelIndex = 1, col_getNumLev(columnAnlInc, 'TH')
                sf_column(levelIndex) = sf_column(levelIndex) + &
                  sf_ad(levelIndex,profileIndex) * tvs_cloudScaleFactor_tl
              end do
            end if

            ! cloud liquid/ice water content
            clw_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'LWCR')
            ciw_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'IWCR')
            do levelIndex = 1, col_getNumLev(columnAnlInc,'TH')
              clw_column(levelIndex) = clw_column(levelIndex) + &
                clwScatt_ad(levelIndex,profileIndex) * tvs_cloudScaleFactor_tl
              ciw_column(levelIndex) = ciw_column(levelIndex) + &
                ciw_ad(levelIndex,profileIndex) * tvs_cloudScaleFactor_tl
            end do
          end if ! surfTypeIsWater
        end do ! profileIndex
      end if ! runObsOperatorWithHydrometeors_ad

      deallocate(tt_ad)
      deallocate(hu_ad)
      deallocate(pressure_ad)
      if (tvs_useO3FromTrials_tl) then
        if (tvs_coefs(sensorIndex) % coef % nozone > 0) then
          deallocate(ozone_ad)
        end if
      end if
      if (allocated(clw_ad)) deallocate(clw_ad)
      if (allocated(clwScatt_ad)) deallocate(clwScatt_ad)
      deallocate(surfTypeIsWater)
      if (allocated(ciw_ad)) then
        deallocate(ciw_ad)
        deallocate(rf_ad)
        deallocate(sf_ad)
      end if
    end do sensor_loop

    ! 3.  Close up

    nullify(profiles)
    if (allocated(tvs_bodyIndexFromBtIndex)) then
      deallocate( tvs_bodyIndexFromBtIndex )
      deallocate( tvs_bodyIndexFromBtIndexScatt )
      deallocate( tvs_chanProf )
      deallocate( tvs_chanProfScatt )
    end if
    write(*,*) 'tvs_rttov_ad: Finished'

  end subroutine tvs_rttov_ad

  !--------------------------------------------------------------------------
  !  tvs_rttov_k
  !--------------------------------------------------------------------------
  subroutine tvs_rttov_k(columnTrlOnAnlIncLev, obsSpaceData)
    !
    ! :Purpose: Compute the Jacobian of radiance using rttov_k
    !
    implicit none
  
    ! Arguments:
    type(struct_obs),        intent(inout) :: obsSpaceData         ! obsSpaceData structure
    type(struct_columnData), intent(in)    :: columnTrlOnAnlIncLev ! column structure for background profile

    ! Locals:
    type(rttov_emissivity), pointer    :: emissivity_local(:)
    type(rttov_emissivity), pointer    :: emissivity_k(:)
    type (rttov_profile), pointer      :: profiles(:)
    type (rttov_profile), pointer      :: profiles_k(:)
    type(rttov_profile_cloud), pointer :: cld_profiles(:)
    type(rttov_transmission)           :: transmission
    type(rttov_transmission)           :: transmission_k
    type(rttov_radiance)               :: radiancedata_d 
    type(rttov_radiance)               :: radiancedata_k 
    integer, allocatable               :: sensorHeaderIndexes(:) 
    integer                            :: allocStatus
    integer                            :: profileCount, btCount, btCountScatt
    integer                            :: sensorIndex
    integer                            :: nlv_T
    integer                            :: instrum
    integer                            :: nthreads
    integer                            :: sensorType
    integer                            :: hydroChannelsCount
    real(8), allocatable               :: surfem1(:)
    integer                            :: errorstatus
    logical                            :: runObsOperatorWithHydrometeors_k
    logical, pointer                   :: calcemis(:)
    logical, allocatable               :: lChannelSubset(:,:)

    if (tvs_headerEnd < 0) return ! exit if there are not tovs data
  
    call tvs_getProfile(profiles, 'nl', cld_profiles)
  
    ! Set index for model's lowest level and model top
  
    nlv_T = col_getNumLev(columnTrlOnAnlIncLev, 'TH')
        
    ! Loop over all sensors specified by user
  
    sensor_loop: do sensorIndex = 1, tvs_nsensors
      
      runObsOperatorWithHydrometeors_k = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
                                         col_varExist(columnTrlOnAnlIncLev,'IWCR') .and. &
                                         tvs_isInstrumUsingHydrometeors(tvs_instruments(sensorIndex))
      profileCount = tvs_countProfiles(sensorIndex)
      call utl_reallocate(sensorHeaderIndexes, profileCount)
      call tvs_setupPointers(runObsOperatorWithHydrometeors_k, sensorIndex, btCount, btCountScatt, &
          hydroChannelsCount, sensorHeaderIndexes, lChannelSubset, obsSpaceData)
      
      if (profileCount == 0) cycle sensor_loop
      if (btCount == 0 .and. btCountScatt == 0) cycle sensor_loop
      
      sensorType = tvs_coefs(sensorIndex) % coef % id_sensor
      instrum = tvs_coefs(sensorIndex) % coef % id_inst

      if (btCount > 0) then
        call rttov_alloc_k(                  &
            allocStatus,                     &
            asw=1,                           & ! to allocate
            nprofiles=profileCount,          &
            nchanprof=btCount,               &
            nlevels=nlv_T,                   &
            opts=tvs_opts(sensorIndex),      &
            profiles_k=profiles_k,           &
            coefs=tvs_coefs(sensorIndex),    &
            transmission=transmission,       &
            transmission_k=transmission_k,   &
            radiance=radiancedata_d,         &
            radiance_k=radiancedata_k,       &
            calcemis=calcemis,               &
            emissivity=emissivity_local,     &
            emissivity_k=emissivity_k,       &
            init=.true.)
        if (allocStatus /= 0) call utl_abort('tvs_rttov_k: memory allocation error in rttov_alloc_k')
          
        ! Set nthreads to actual number of threads which will be used.
        nthreads = min(mmpi_numThread, profileCount)  
  
        ! Prepare all input variables required by rttov.
    
        write(*,*) 'tvs_rttov_k: Get surface emissiviy'
        allocate(surfem1(btCount))
        !    get Hyperspecral IR emissivities
        if (tvs_isInstrumHyperSpectral(instrum)) call tvs_getHIREmissivities(sensorHeaderIndexes, &
                                                          obsSpaceData, surfem1)
  
        call tvs_getOtherEmissivities(tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes, sensorType, instrum, surfem1, calcemis)
  
        if (sensorType == sensor_id_mw) then
          call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorIndex, tvs_chanProf(1:btCount,sensorIndex), sensorHeaderIndexes)
        else
          emissivity_local(:) % emis_in = surfem1(:)
        end if
   
        ! Intialize variables'
        emissivity_k % emis_in = 0.0
        emissivity_k % emis_out = 0.0
  
        call rttov_init_prof(profiles_k)
        call rttov_init_transmission(transmission_k)
        call rttov_init_rad(radiancedata_k)

        radiancedata_k % bt = 1.0
        radiancedata_k % total(:) = 1.0

        call rttov_parallel_k(                               &
            errorstatus,                                     & ! out
            tvs_chanProf(1:btCount,sensorIndex),             & ! in
            tvs_opts(sensorIndex),                           & ! in
            profiles(sensorHeaderIndexes(:)),                & ! in
            profiles_k,                                      & ! inout
            tvs_coefs(sensorIndex),                          & ! in
            transmission,                                    & ! inout
            transmission_k,                                  & ! inout
            radiancedata_d,                                  & ! inout
            radiancedata_k,                                  & ! inout
            calcemis=calcemis,                               & ! in
            emissivity=emissivity_local,                     & ! in
            emissivity_k=emissivity_k,                       & ! inout
            nthreads=nthreads )                                ! in
                   
        if (errorstatus /= errorStatus_success) then
          write(*,*) "Error in rttov_parallel_k", errorstatus
          write(*,*) 'temperature profile=', profiles(sensorHeaderIndexes(1)) % t(:)
          write(*,*) 'temperature Jacobian profile=', profiles_k(1) % t(:)
          call utl_abort('tovs_rttov_k')
        end if
  
        ! Write Jacobian to ASCII files
        call tvs_writeJacobianAscii(profiles_k, emissivity_k, profiles, tvs_chanProf(1:btCount,sensorIndex), &
            obsSpaceData, tvs_satelliteName(sensorIndex), tvs_instrumentName(sensorIndex), &
            tvs_bodyIndexFromBtIndex(:,sensorIndex), sensorHeaderIndexes, btCount)
  
        ! deallocate profiledata structures
        call rttov_alloc_k(                  &
            allocStatus,                     &
            asw=0,                           &
            nprofiles=profileCount,          &
            nchanprof=btCount,               &
            nlevels=nlv_T,                   &
            opts=tvs_opts(sensorIndex),      &
            profiles_k=profiles_k,           &
            coefs=tvs_coefs(sensorIndex),    &
            transmission=transmission,       &
            transmission_k=transmission_k,   &
            radiance=radiancedata_d,         &
            radiance_k=radiancedata_k,       &
            calcemis=calcemis,               &
            emissivity=emissivity_local,     &
            emissivity_k=emissivity_k)
        if (allocStatus /= 0) call utl_abort('tvs_rttov_k: memory deallocation error in rttov_alloc_k')
        deallocate(surfem1)
      end if  !if (btCount > 0)
      
      if (btCountScatt > 0) then
        call utl_abort("tvs_rttov_k: jacobians not (yet) available when rttov_scatt is used !")
      end if
    end do sensor_loop
  
    nullify(profiles)
    if (allocated(tvs_bodyIndexFromBtIndex)) then
      deallocate( tvs_bodyIndexFromBtIndex )
      deallocate( tvs_bodyIndexFromBtIndexScatt )
      deallocate( tvs_chanProf )
      deallocate( tvs_chanProfScatt )
    end if
    write(*,*) 'tvs_rttov_k: finished'
    
  end subroutine tvs_rttov_k

end module tovs_mod

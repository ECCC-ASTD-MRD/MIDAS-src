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

module tovs_nl_mod
  ! MODULE tovs_nl (prefix='tvs' category='5. Observation operators')
  !
  ! :Purpose: Derived types, public variables and procedures related to the
  !           nonlinear version of RTTOV
  !

  use rttov_interfaces_mod
  use rttov_types, only :   &
       rttov_coefs         ,&
       rttov_fast_coef     ,&
       rttov_scatt_coef    ,&
       rttov_options       ,&
       rttov_options_scatt ,&
       rttov_profile       ,&
       rttov_profile_cloud ,&
       rttov_radiance      ,&
       rttov_transmission  ,&
       rttov_chanprof      ,&
       rttov_emissivity

  use rttov_const, only :   &
       platform_name       ,&
       nplatforms          ,&
       inst_name           ,&
       ninst               ,&
       inst_id_goesim      ,&
       inst_id_gmsim       ,&
       inst_id_mtsatim     ,&
       inst_id_amsua       ,&
       inst_id_mhs         ,&
       sensor_id_mw        ,&
       sensor_id_po        ,&
       platform_id_jpss    ,&
       platform_id_himawari,&
       platform_id_eos     ,&
       errorstatus_success ,&
       mair, mh2o, mo3     ,&
       surftype_land       ,&
       surftype_seaice     ,&
       surftype_sea        ,&
       ngases_max          ,&
       gas_id_mixed        ,&
       gas_unit_specconc   ,&
       interp_rochon_loglinear_wfn, &
       zenmax

  use parkind1, only : jpim, jplm

  use rttov_fast_coef_utils_mod, only: set_pointers, set_fastcoef_level_bounds
  use rttov_solar_refl_mod, only : rttov_refl_water_interp

  use midasMpi_mod
  use codtyp_mod
  use mpi
  use utilities_mod
  use obsSpaceData_mod
  use earthConstants_mod
  use MathPhysConstants_mod
  use ozoneclim_mod
  use columnData_mod 
  use mod_rttov_emis_atlas
  use verticalCoord_mod
  use codePrecision_mod
  use humidityLimits_mod
  use interpolation_mod

  implicit none
  save
  private

  type surface_params
    real(8)              :: albedo   ! surface albedo (0-1)
    real(8)              :: ice      ! ice cover (0-1) 
    real(8)              :: snow     ! snow cover (0-1)
    real(8)              :: pcnt_wat ! water percentage in pixel containing profile (0-1)
    real(8)              :: pcnt_reg ! water percentage in an area around profile (0-1)
    integer              :: ltype    ! surface type (1,...,20)
  end type surface_params

  ! public variables (parameters)
  public :: tvs_maxChannelNumber, tvs_maxNumberOfChannels, tvs_maxNumberOfSensors, tvs_defaultEmissivity
  ! public variables (non-parameters)
  public :: tvs_nchan, tvs_ichan, tvs_lsensor, tvs_headerIndex, tvs_tovsIndex, tvs_nobtov
  public :: tvs_nchanMpiGlobal, tvs_ichanMpiGlobal
  public :: tvs_isReallyPresent,tvs_listSensors
  public :: tvs_isReallyPresentMpiGlobal
  public :: tvs_nsensors, tvs_platforms, tvs_satellites, tvs_instruments, tvs_channelOffset
  public :: tvs_debug, tvs_satelliteName, tvs_instrumentName, tvs_useO3Climatology
  public :: tvs_coefs, tvs_opts, tvs_transmission,tvs_emissivity
  public :: tvs_coef_scatt, tvs_opts_scatt
  public :: tvs_radiance, tvs_surfaceParameters
  public :: tvs_numMWInstrumUsingCLW, tvs_numMWInstrumUsingHydrometeors
  public :: tvs_mwInstrumUsingCLW_tl, tvs_mwInstrumUsingHydrometeors_tl
  public :: tvs_mwAllskyAssim
  ! public procedures
  public :: tvs_fillProfiles, tvs_rttov, tvs_printDetailledOmfStatistics, tvs_allocTransmission, tvs_cleanup
  public :: tvs_deallocateProfilesNlTlAd
  public :: tvs_setupAlloc,tvs_setup, tvs_isIdBurpTovs, tvs_isIdBurpHyperSpectral, tvs_isIdBurpInst, tvs_getAllIdBurpTovs
  public :: tvs_isInstrumGeostationary,  tvs_isNameHyperSpectral
  public :: tvs_isNameGeostationary
  public :: tvs_getInstrumentId, tvs_getPlatformId, tvs_mapSat, tvs_mapInstrum
  public :: tvs_isInstrumHyperSpectral, tvs_getChanprof, tvs_countRadiances
  public :: tvs_ChangedStypValue
  public :: tvs_getHIREmissivities, tvs_getOtherEmissivities, tvs_rttov_read_coefs
  public :: tvs_getLocalChannelIndexFromChannelNumber
  public :: tvs_getMWemissivityFromAtlas, tvs_getProfile
  public :: tvs_getCorrectedSatelliteAzimuth
  public :: tvs_isInstrumUsingCLW, tvs_isInstrumUsingHydrometeors, tvs_getChannelNumIndexFromPPP
  public :: tvs_isInstrumAllskyTtAssim, tvs_isInstrumAllskyHuAssim, tvs_isInstrumAllskyTtHuAssim
  ! Module parameters
  ! units conversion from  mixing ratio to ppmv and vice versa
  real(8), parameter :: qMixratio2ppmv  = (1000000.0d0 * mair) / mh2o
  real(8), parameter :: qppmv2Mixratio  = mh2o / (1000000.0d0 * mair)
  real(8), parameter :: o3Mixratio2ppmv = (1000000.0d0 * mair) / mo3
  real(8), parameter :: o3ppmv2Mixratio = mo3 / (1000000.0d0 * mair)
  real(pre_obsReal), parameter :: tvs_defaultEmissivity = 0.95

  integer, parameter :: tvs_maxChannelNumber   = 8461   ! Max. value for channel number
  integer, parameter :: tvs_maxNumberOfChannels = 2211  ! Max. no. of channels (for one profile/spectra)
  integer, parameter :: tvs_maxNumberOfSensors  = 100   ! Max no sensors to be used
  integer, parameter :: tvs_nlevels     = 101           ! Maximum No. of RTTOV pressure levels including 'rttov top' at 0.005 hPa

  ! Module variables
  integer, allocatable :: tvs_bodyIndexFromBtIndex(:)  ! Provides the bodyIndex in ObsSpaceData based on btIndex
  integer, allocatable :: tvs_nchan(:)             ! Number of channels per instrument (local)
  integer, allocatable :: tvs_ichan(:,:)           ! List of channels per instrument (local)
  integer, allocatable :: tvs_nchanMpiGlobal(:)     ! Number of channels per instrument (global)
  integer, allocatable :: tvs_ichanMpiGlobal(:,:)  ! List of channels per instrument  (global)
  integer, allocatable :: tvs_lsensor(:)           ! Sensor number for each profile
  integer, allocatable :: tvs_headerIndex (:)      ! Observation position in obsSpaceData header for each profile
  integer, allocatable :: tvs_tovsIndex (:)        ! Index in TOVS structures for each observation header in obsSpaceData
  logical, allocatable :: tvs_isReallyPresent(:)   ! Logical flag to identify instruments really assimilated (local)
  logical, allocatable :: tvs_isReallyPresentMpiGLobal(:)   ! Logical flag to identify instruments really assimilated (global)
  integer, allocatable :: tvs_listSensors(:,:)     ! Sensor list
  type (rttov_emis_atlas_data), allocatable :: tvs_atlas(:)     ! Emissivity atlases
  type(surface_params), allocatable :: tvs_surfaceParameters(:) ! surface parameters 
  integer tvs_nobtov                               ! Number of tovs observations
  integer tvs_nsensors                             ! Number of individual sensors.
  integer tvs_platforms(tvs_maxNumberOfSensors)    ! RTTOV platform ID's (e.g., 1=NOAA; 2=DMSP; ...)
  integer tvs_satellites(tvs_maxNumberOfSensors)   ! RTTOV satellite ID's (e.g., 1 to 16 for NOAA; ...)
  integer tvs_instruments(tvs_maxNumberOfSensors)  ! RTTOVinstrument ID's (e.g., 3=AMSU-A; 4=AMSU-B; 6=SSMI; ...)
  integer tvs_channelOffset(tvs_maxNumberOfSensors)! BURP to RTTOV channel mapping offset
  integer instrumentIdsUsingCLW(tvs_maxNumberOfSensors)
  integer instrumentIdsUsingHydrometeors(tvs_maxNumberOfSensors)
  integer tvs_numMWInstrumUsingCLW 
  integer tvs_numMWInstrumUsingHydrometeors
  logical tvs_mwInstrumUsingCLW_tl
  logical tvs_mwInstrumUsingHydrometeors_tl
  logical tvs_mwAllskyAssim
  logical :: tvs_mpiTask0ReadCoeffs 
  real(8) :: tvs_cloudScaleFactor 
  logical tvs_debug                                ! Logical key controlling statements to be  executed while debugging TOVS only
  logical tvs_useO3Climatology                     ! Determine if ozone model field or climatology is used
                                                   ! If ozone model field is specified, related increments will be generated in assimilation
  logical tvs_regLimitExtrap                       ! use RTTOV reg_limit_extrap option
  logical tvs_doAzimuthCorrection(tvs_maxNumberOfSensors)
  logical tvs_isAzimuthValid(tvs_maxNumberOfSensors)
  logical tvs_userDefinedDoAzimuthCorrection
  logical tvs_userDefinedIsAzimuthValid

  character(len=15) tvs_satelliteName(tvs_maxNumberOfSensors)
  character(len=15) tvs_instrumentName(tvs_maxNumberOfSensors)
  character(len=8) radiativeTransferCode           ! RadiativeTransferCode : TOVS radiation model used
  real(8), allocatable :: tvs_emissivity(:,:)      ! Surface emissivities organized by profiles and channels
  integer, parameter :: kslon=2160, kslat=1080     ! CERES file dimension in grid points

  ! High resolution surface fields
  integer :: surfaceType(kslon,kslat)  
  real(8) :: waterFraction(kslon,kslat) 

  ! Derived typeso
  type(rttov_coefs), allocatable           :: tvs_coefs(:)          ! rttov coefficients
  type(rttov_options), allocatable         :: tvs_opts(:)           ! rttov options
  type(rttov_scatt_coef),allocatable       :: tvs_coef_scatt(:)     ! rttovscatt coefficients
  type(rttov_options_scatt), allocatable   :: tvs_opts_scatt(:)     ! rttovscatt options
  type(rttov_profile), target, allocatable :: tvs_profiles_nl(:)    ! all profiles on trial vertical coordinate for nl obs operator
  type(rttov_profile), target, allocatable :: tvs_profiles_tlad(:)  ! all profiles on increments vertical coordinates for linearized obs. operator
  type(rttov_radiance), allocatable        :: tvs_radiance(:)       ! radiances organized by profile
  type(rttov_transmission), allocatable    :: tvs_transmission(:)   ! transmittances all profiles for HIR quality control
  type(rttov_profile_cloud), target, allocatable :: tvs_cld_profiles_nl(:)! rttov scatt cloud profiles on trial vertical coordinate
  type(rttov_profile_cloud), target, allocatable :: tvs_cld_profiles_tlad(:) ! rttov scatt cloud profiles on increment vertical coordinates

  ! Namelist variables:
  logical useUofWIREmiss                           ! Flag to activate use of RTTOV U of W IR emissivity Atlases
  logical useMWEmissivityAtlas                     ! Flag to activate use of RTTOV built-in MW emissivity Atlases      
  integer mWAtlasId                                ! MW Atlas Id used when useMWEmissivityAtlas == .true. ; 1 TELSEM2, 2 CNRM atlas

  integer, external :: get_max_rss
 
contains

  !--------------------------------------------------------------------------
  ! tvs_setupAlloc
  !--------------------------------------------------------------------------
  subroutine tvs_setupAlloc(obsSpaceData)
    !
    ! :Purpose: Memory allocation for the non linear radiative transfer model variables.
    !
    implicit none

    !Arguments:
    type(struct_obs) :: obsSpaceData

    ! Locals:
    integer :: allocStatus(9)
    integer :: satelliteCode, instrumentCode, iplatform, isat, instrum
    integer :: tovsIndex, idatyp, sensorIndex
    integer :: channelNumber, nosensor, channelIndex
    integer :: errorStatus(1)
    integer :: headerIndex, bodyIndex, taskIndex
    logical :: runObsOperatorWithClw
    logical :: runObsOperatorWithHydrometeors
    logical, allocatable :: logicalBuffer(:)
    character(len=32) :: hydroTableFilename

    if (tvs_nsensors == 0) return

    !  1. Determine the number of radiances to be assimilated.
    !      Construct a list of channels for each sensor.
    !      Construct a list of sensor number for each profile

    write(*,*) 'tvs_setupAlloc: Starting' 

    allocStatus = 0
    allocate (tvs_nchan(tvs_nsensors),                         stat= allocStatus(1))
    allocate (tvs_ichan(tvs_maxNumberOfChannels,tvs_nsensors), stat= allocStatus(2))
    allocate (tvs_lsensor(obs_numheader(obsSpaceData)),        stat= allocStatus(3))
    allocate (tvs_headerIndex (obs_numheader(obsSpaceData)),   stat= allocStatus(4))
    allocate (tvs_tovsIndex(obs_numheader(obsSpaceData)),      stat= allocStatus(5))
    allocate (tvs_isReallyPresent(tvs_nsensors),               stat= allocStatus(6))
    allocate (tvs_nchanMpiGlobal(tvs_nsensors),                stat= allocStatus(7))
    allocate (tvs_ichanMpiGlobal(tvs_maxNumberOfChannels,tvs_nsensors), stat= allocStatus(8))
    allocate (tvs_isReallyPresentMpiGlobal(tvs_nsensors), stat= allocStatus(9))

    call utl_checkAllocationStatus(allocStatus, ' tvs_setupAlloc')

    tvs_nchan(:) = 0 
    tvs_ichan(:,:) = 0
    tvs_isReallyPresent(:) = .true.
    tvs_lsensor(:) = -1
    tvs_headerIndex(:) = -1
    tvs_tovsIndex (:) = -1

    tvs_nobtov = 0

    ! Loop over all header indices of the 'TO' family
    ! Set the header list & start at the beginning of the list
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)

      if ( .not. tvs_isIdBurpTovs(idatyp) ) then
        write(*,*) 'tvs_setupAlloc: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY2:do 
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY2
          call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
        end do BODY2
        cycle HEADER   ! Proceed to the next headerIndex
      end if
      tvs_nobtov = tvs_nobtov + 1
     
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
      nosensor =0
      do sensorIndex = 1, tvs_nsensors
        if ( iplatform == tvs_platforms  (sensorIndex) .and. &
             isat      == tvs_satellites (sensorIndex) .and. &
             instrum   == tvs_instruments(sensorIndex)      ) then
          nosensor = sensorIndex
          exit
        end if
      end do

      if (nosensor > 0) then
        tvs_lsensor(tvs_nobtov) = nosensor
        tvs_headerIndex (tvs_nobtov) = headerIndex
        tvs_tovsIndex (headerIndex) = tvs_nobtov
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
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) then
            call tvs_getChannelNumIndexFromPPP( obsSpaceData, headerIndex, bodyIndex, &
                                                channelNumber, channelIndex )
            if ( channelIndex == 0 ) then
              tvs_nchan(nosensor) = tvs_nchan(nosensor) + 1
              tvs_ichan(tvs_nchan(nosensor),nosensor) = channelNumber
            end if
            
            if ( tvs_debug .and. mmpi_myid == 0 .and. &
                 trim(tvs_instrumentName(nosensor)) == 'AMSUA' ) then
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

    if ( .not. tvs_userDefinedDoAzimuthCorrection) then 
      ! tvs_doAzimuthCorrection user defined values will be overwriten by the old default values 
      do sensorIndex = 1, tvs_nsensors
        tvs_doAzimuthCorrection(sensorIndex) = ( tvs_platforms(sensorIndex) /= platform_id_eos .and. &
             ( tvs_instruments(sensorIndex) == inst_id_amsua .or. tvs_instruments(sensorIndex) == inst_id_mhs )  )     
      end do
      if ( mmpi_myId == 0 ) write(*,*) ' tvs_setupAlloc: Warning tvs_doAzimuthCorrection user defined values overwriten by the old default values'
    end if

    if ( .not. tvs_userDefinedIsAzimuthValid ) then 
      ! tvs_isAzimuthValid  user defined values will be overwriten by the current default values 
      do sensorIndex = 1, tvs_nsensors
        tvs_isAzimuthValid(sensorIndex) = .not. ( tvs_isInstrumGeostationary(tvs_instruments(sensorIndex)) )
      end do
      if ( mmpi_myId == 0 ) write(*,*) ' tvs_setupAlloc: Warning tvs_isAzimuthValid user defined values overwriten by the current default values'
    end if

    if ( mmpi_myId == 0 ) then
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
      if ( tvs_nchan(sensorIndex) == 0 ) then
        tvs_isReallyPresent ( sensorIndex ) =.false.
        tvs_nchan(sensorIndex) = 1
        tvs_ichan(1,sensorIndex) = 1
      end if
    end do

    write(*,*) ' tvs_setupAlloc: tvs_nobtov = ', tvs_nobtov

    do sensorIndex = 1, tvs_nsensors
      call tvs_getCommonChannelSet(tvs_ichan(:,sensorIndex),tvs_nchanMpiGlobal(sensorIndex), tvs_ichanMpiGlobal(:,sensorIndex))
      print *, 'sensorIndex,tvs_nchan(sensorIndex),tvs_nchanMpiGlobal(sensorIndex)', sensorIndex, tvs_nchan(sensorIndex),tvs_nchanMpiGlobal(sensorIndex)
    end do

    if (mmpi_myid ==0) then
      allocate(logicalBuffer(mmpi_nprocs))
    else
      allocate(logicalBuffer(1))
    end if
    
    do sensorIndex = 1, tvs_nsensors
      call RPN_COMM_gather( tvs_isReallyPresent ( sensorIndex ) , 1, 'MPI_LOGICAL', logicalBuffer, 1,'MPI_LOGICAL', 0, 'GRID', errorStatus(1) )
      if (mmpi_myid ==0) then
        tvs_isReallyPresentMpiGlobal ( sensorIndex ) =.false.
        do taskIndex=1, mmpi_nprocs
          tvs_isReallyPresentMpiGlobal ( sensorIndex ) =  tvs_isReallyPresentMpiGlobal ( sensorIndex ) .or. logicalBuffer(taskIndex)
        end do
      end if
      call rpn_comm_bcast(tvs_isReallyPresentMpiGlobal ( sensorIndex ), 1, 'MPI_LOGICAL', 0, 'GRID', errorStatus(1) )
    end do
    
    deallocate(logicalBuffer)

    if ( tvs_debug .and. mmpi_myid == 0 ) then
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

    if ( radiativeTransferCode == 'RTTOV' ) then

      write(*,'(//,10x,A)') '-rttov_setup: initializing the TOVS radiative transfer model'

      allocate (tvs_coefs(tvs_nsensors)          ,stat= allocStatus(1))
      allocate (tvs_listSensors (3,tvs_nsensors) ,stat= allocStatus(2))
      allocate (tvs_opts (tvs_nsensors)          ,stat= allocStatus(3))
      if (tvs_numMWInstrumUsingHydrometeors  > 0) then
        allocate (tvs_opts_scatt (tvs_nsensors) ,stat= allocStatus(4))
        allocate (tvs_coef_scatt (tvs_nsensors) ,stat= allocStatus(5))
      end if
      call utl_checkAllocationStatus(allocStatus(1:5), ' tvs_setupAlloc before rttov initialization')

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
        write(*,*) ' sensorIndex,tvs_nchan(sensorIndex)',  sensorIndex,tvs_nchan(sensorIndex)
        if ( tvs_mpiTask0ReadCoeffs ) then
          call tvs_rttov_read_coefs(errorStatus(1), tvs_coefs(sensorIndex), tvs_opts(sensorIndex), & 
               tvs_ichan(1:tvs_nchan(sensorIndex),sensorIndex), tvs_listSensors(:,sensorIndex))
        else
          call rttov_read_coefs (                               &
               errorStatus(1),                                  &! out
               tvs_coefs(sensorIndex),                          &
               tvs_opts(sensorIndex),                           &
               instrument= tvs_listSensors(:,sensorIndex),      &! in
               channels=  tvs_ichan(1:tvs_nchan(sensorIndex),sensorIndex) )     ! in option
        end if
        if (errorStatus(1) /= errorStatus_success) then
          write(*,*) 'rttov_read_coefs: fatal error reading coefficients',errorStatus,sensorIndex,tvs_listSensors(1:3,sensorIndex)
          call utl_abort('tvs_setupAlloc')
        end if
       
        if (runObsOperatorWithHydrometeors) then
          hydrotableFilename = 'hydrotable_' // trim(platform_name(tvs_platforms(sensorIndex))) // '_' // &
               trim(inst_name(tvs_instruments(sensorIndex))) // '.dat'
          call rttov_read_scattcoeffs(errorstatus(1), tvs_opts_scatt(sensorIndex), tvs_coefs(sensorIndex), &
               tvs_coef_scatt(sensorIndex), file_coef=hydrotableFilename)
          if (errorstatus(1) /= errorstatus_success) then
            write(*,*) 'rttov_read_scattcoeffs: fatal error reading RTTOV-SCATT coefficients', hydrotableFilename
            call utl_abort('tvs_setupAlloc')
          end if
        end if
        call utl_tmg_stop(16)

        tvs_opts(sensorIndex) % rt_all % ozone_data = ( tvs_coefs(sensorIndex) % coef % nozone > 0 ) ! profil d'ozone disponible

        ! Ensure the options and coefficients are consistent
        call rttov_user_options_checkinput(errorStatus(1), tvs_opts(sensorIndex), tvs_coefs(sensorIndex))
        if (errorStatus(1) /= errorStatus_success) then
          write(*,*) 'error in rttov options',errorStatus
          call utl_abort('tvs_setupAlloc')
        end if
       
      end do


      !   4. Memory allocations for radiative tranfer model variables

      ! Radiance by profile

      allocate( tvs_radiance(tvs_nobtov) ,stat= allocStatus(1))

      call utl_checkAllocationStatus(allocStatus(1:1), ' tvs_setupAlloc radiances 1')
  
      do tovsIndex = 1, tvs_nobtov
        sensorIndex = tvs_lsensor(tovsIndex)
        if (sensorIndex > -1) then
          ! allocate BT equivalent to total direct, tl and ad radiance output
          allocate( tvs_radiance(tovsIndex)  % bt  ( tvs_nchan(sensorIndex) ) ,stat= allocStatus(1))
          tvs_radiance(tovsIndex)  % bt  ( : ) = 0.d0
          call utl_checkAllocationStatus(allocStatus(1:1), ' tvs_setupAlloc radiances 2')
          nullify (tvs_radiance(tovsIndex)  % clear )
        end if
      end do

    end if
  
    write(*,*) 'Leaving tvs_setupAlloc'

  end subroutine tvs_setupAlloc

  !--------------------------------------------------------------------------
  ! tvs_getProfile
  !--------------------------------------------------------------------------
  subroutine tvs_getProfile(profiles, profileType, cld_profiles_opt)
    ! :Purpose: sets profiles as a pointer of type rttov_profile
    !           based on profileType equal to nl or tlad. 
    ! 
    implicit none

    !Arguments:
    type(rttov_profile), pointer, intent(inout)       :: profiles(:)
    type(rttov_profile_cloud), pointer, intent(inout), optional :: cld_profiles_opt(:)
    character(len=*), intent(in)                      :: profileType

    select case( trim( profileType) )
      case('nl')
        profiles => tvs_profiles_nl
        if (present(cld_profiles_opt)) cld_profiles_opt => tvs_cld_profiles_nl
      case('tlad')
        profiles => tvs_profiles_tlad
        if (present(cld_profiles_opt)) cld_profiles_opt => tvs_cld_profiles_tlad
      case default
        call utl_abort('tvs_getProfile: invalid profileType ' // profileType )
    end select

  end subroutine tvs_getProfile

  !--------------------------------------------------------------------------
  ! tvs_allocTransmission
  !--------------------------------------------------------------------------
  subroutine tvs_allocTransmission(nlevels)

    ! :Purpose: Allocate the global rttov transmission structure used
    !           when this is needed for some purpose (e.g. used in 
    !           LETKF to determine peak pressure level of each radiance
    !           channel for vertical localization).
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: nlevels 
    ! Locals:
    integer :: allocStatus(2), jo, isens, nc

    allocStatus(:) = 0
    allocate( tvs_transmission(tvs_nobtov), stat=allocStatus(1))
    call utl_checkAllocationStatus(allocStatus(1:1), ' tvs_allocTransmission')

    do jo = 1, tvs_nobtov
      isens = tvs_lsensor(jo)
      nc = tvs_nchan(isens)
      ! allocate transmittance from surface and from pressure levels
      allocate( tvs_transmission(jo) % tau_total(nc),     stat= allocStatus(1))
      allocate( tvs_transmission(jo) % tau_levels(nlevels,nc), stat= allocStatus(2))
      call utl_checkAllocationStatus(allocStatus, ' tvs_allocTransmission')
    end do

  end subroutine tvs_allocTransmission



  !--------------------------------------------------------------------------
  ! tvs_setup
  !--------------------------------------------------------------------------
  subroutine tvs_setup
    !
    ! :Purpose: to read namelist NAMTOV, initialize the observation error covariance and setup RTTOV-12.
    !

    implicit none
    ! Locals:
    integer  sensorIndex, nsensors, ierr, nulnam
    integer, external :: fclos, fnom
    integer :: instrumentIndex, numMWInstrumToUseCLW, numMWInstrumToUseHydrometeors

    ! Namelist variables: (local)
    character(len=8)  :: crtmodl
    character(len=15) :: csatid(tvs_maxNumberOfSensors)
    character(len=15) :: cinstrumentid(tvs_maxNumberOfSensors)
    logical :: ldbgtov
    logical :: useO3Climatology
    logical :: regLimitExtrap
    logical :: doAzimuthCorrection(tvs_maxNumberOfSensors)
    logical :: isAzimuthValid(tvs_maxNumberOfSensors)
    logical :: userDefinedDoAzimuthCorrection
    logical :: userDefinedIsAzimuthValid
    logical :: mpiTask0ReadCoeffs
    logical :: mwInstrumUsingCLW_tl
    logical :: mwInstrumUsingHydrometeors_tl
    real(8) :: cloudScaleFactor 
    character(len=15) :: instrumentNamesUsingCLW(tvs_maxNumberOfSensors)
    character(len=15) :: instrumentNamesUsingHydrometeors(tvs_maxNumberOfSensors)
    logical :: mwAllskyAssim

    namelist /NAMTOV/ csatid, cinstrumentid
    namelist /NAMTOV/ ldbgtov,useO3Climatology
    namelist /NAMTOV/ useUofWIREmiss, crtmodl
    namelist /NAMTOV/ useMWEmissivityAtlas, mWAtlasId
    namelist /NAMTOV/ mwInstrumUsingCLW_tl, instrumentNamesUsingCLW
    namelist /NAMTOV/ mwInstrumUsingHydrometeors_tl, instrumentNamesUsingHydrometeors
    namelist /NAMTOV/ regLimitExtrap, doAzimuthCorrection, userDefinedDoAzimuthCorrection
    namelist /NAMTOV/ isAzimuthValid, userDefinedIsAzimuthValid, cloudScaleFactor 
    namelist /NAMTOV/ mwAllskyAssim, mpiTask0ReadCoeffs

    ! return if the NAMTOV does not exist
    if ( .not. utl_isNamelistPresent('NAMTOV','./flnml') ) then
      write(*,*)
      write(*,*) 'tvs_setup: Namelist block NAMTOV is missing in the namelist.'
      write(*,*) '           Skipping tvs_setup.'
      return
    end if
 
    !   1.1 Default values for namelist variables

    csatid(:) = '***UNDEFINED***'
    cinstrumentid(:) = '***UNDEFINED***'
    doAzimuthCorrection(:) = .false.
    isAzimuthValid(:) = .false.
    !csatid(1) = 'NOAA16'
    !cinstrumentid(1) = 'AMSUA'
    ldbgtov = .false.
    useO3Climatology = .true.
    userDefinedDoAzimuthCorrection = .false.
    userDefinedIsAzimuthValid = .false.
    crtmodl = 'RTTOV'
    useUofWIREmiss = .false.
    useMWEmissivityAtlas = .false.
    mWAtlasId = 1 !Default to TELSEM-2
    mwInstrumUsingCLW_tl = .false.
    mwInstrumUsingHydrometeors_tl = .false.
    instrumentNamesUsingCLW(:) = '***UNDEFINED***'
    instrumentNamesUsingHydrometeors(:) = '***UNDEFINED***'
    regLimitExtrap = .true.
    cloudScaleFactor = 0.5D0
    mwAllskyAssim = .false.
    mpiTask0ReadCoeffs = .true.

    !   1.2 Read the NAMELIST NAMTOV to modify them
 
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam, nml=namtov, iostat=ierr)
    if (ierr /= 0) call utl_abort('tvs_setup: Error reading namelist NAMTOV')

    if (mmpi_myid == 0) write(*,nml=namtov)
    ierr = fclos(nulnam)

    !  1.3 Transfer namelist variables to module variables
 
    tvs_nsensors = 0
    sensor_loop:do sensorIndex = 1, tvs_maxNumberOfSensors
      if (cinstrumentid(sensorIndex) /= "***UNDEFINED***" .and. &
          csatid(sensorIndex) /= "***UNDEFINED***" ) then
        tvs_nsensors = tvs_nsensors + 1
      else
        exit sensor_loop
      end if
    end do sensor_loop
    nsensors = 0
    do sensorIndex = 1, tvs_maxNumberOfSensors
      if (cinstrumentid(sensorIndex) /= "***UNDEFINED***" .and. &
          csatid(sensorIndex) /= "***UNDEFINED***" ) then
        nsensors = nsensors + 1
      end if
    end do
    if (nsensors /= tvs_nsensors) then
      call utl_abort('tvs_setup: check namelist cinstrumentid and csatid arrays indexes should be defined from 1 with increasing values (no hole)')
    end if
    
    tvs_debug = ldbgtov
    radiativeTransferCode = crtmodl
    tvs_useO3Climatology = useO3Climatology
    tvs_instrumentName(:) = cinstrumentid(:)
    tvs_satelliteName(:) = csatid(:)
    tvs_mwInstrumUsingCLW_tl = mwInstrumUsingCLW_tl
    tvs_regLimitExtrap = regLimitExtrap
    tvs_userDefinedDoAzimuthCorrection = userDefinedDoAzimuthCorrection
    tvs_userDefinedIsAzimuthValid = userDefinedIsAzimuthValid
    tvs_doAzimuthCorrection(:) =  doAzimuthCorrection(:)
    tvs_isAzimuthValid(:) =  isAzimuthValid(:)
    tvs_cloudScaleFactor = cloudScaleFactor 
    tvs_mwAllskyAssim = mwAllskyAssim
    tvs_mpiTask0ReadCoeffs = mpiTask0ReadCoeffs

    !  1.4 Validate namelist values
    
    if ( tvs_nsensors == 0 ) then
      if(mmpi_myid==0) then 
        write(*,*) ' ====================================================='
        write(*,*) ' tvs_setup: Number of sensors is zero, skipping setup'
        write(*,*) ' ====================================================='
      end if
      return
    end if

    if ( radiativeTransferCode /= 'RTTOV' ) then
      write(*,'(A)') ' Invalid radiation model name'
      call utl_abort('tvs_setup')
    end if

    if ( tvs_nsensors > tvs_maxNumberOfSensors ) then
      write(*,'(A)') ' Number of sensors (tvs_nsensors) is greater than maximum allowed (tvs_maxNumberOfSensors)'
      call utl_abort('tvs_setup')
    end if

    !  1.5 Print the content of this NAMELIST

    if(mmpi_myid == 0) then
      write(*,'(A)') 
      write(*,'(3X,A)') '- Parameters used for TOVS processing (read in NAMTOV)'
      write(*,'(3X,A)') '  ----------------------------------------------------'
      write(*,'(6X,A,2X,L1)') 'TOVS debug                           : ', tvs_debug
      write(*,'(6X,A,2X,L1)') 'Use of UW IR land emissivity atlases : ', useUofWIREmiss
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
      if ( instrumentNamesUsingCLW(instrumentIndex) /= '***UNDEFINED***' ) then
        if ( instrumentIdsUsingCLW(instrumentIndex) == -1 ) then
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
      if ( instrumentNamesUsingHydrometeors(instrumentIndex) /= '***UNDEFINED***' ) then
        if ( instrumentIdsUsingHydrometeors(instrumentIndex) == -1 ) then
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
      if (numMWInstrumToUseCLW == 0 .or. &
          instrumentIdsUsingHydrometeors(instrumentIndex) == tvs_getInstrumentId('atms')) then
        exit
      end if

      if (any(instrumentIdsUsingCLW(1:numMWInstrumToUseCLW) == &
              instrumentIdsUsingHydrometeors(instrumentIndex))) then
        write(*,*) instrumentIndex, instrumentNamesUsingHydrometeors(instrumentIndex)
        call utl_abort('tvs_setup: this instrument is mentioned in instrumentIdsUsingCLW namelist')
      end if
    end do

    tvs_numMWInstrumUsingCLW = numMWInstrumToUseCLW
    tvs_numMWInstrumUsingHydrometeors = numMWInstrumToUsehydrometeors

    if ( mmpi_myid == 0 ) then
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
    ! :Purpose: release memory used by RTTOV-12.
    !
    implicit none
    integer :: allocStatus(8)
    integer :: iSensor,iObs,nChans,nl

    Write(*,*) 'tvs_cleanup: Starting'

    allocStatus(:) = 0

    if ( radiativeTransferCode == 'RTTOV' ) then

      !___ radiance by profile

      do iObs = 1, tvs_nobtov
        iSensor = tvs_lsensor(iObs)
        nchans = tvs_nchan(isensor)
        ! deallocate BT equivalent to total direct, tl and ad radiance output
        deallocate( tvs_radiance(iObs)  % bt ,stat= allocStatus(1))
        call utl_checkAllocationStatus(allocStatus(1:1), ' tvs_cleanup radiances 1',.false.)
      end do

      deallocate( tvs_radiance ,stat= allocStatus(1))
      call utl_checkAllocationStatus(allocStatus(1:1), ' tvs_cleanup radiances 2')

      do iObs = 1, tvs_nobtov
        iSensor = tvs_lsensor(iObs)
        nl = tvs_coefs(iSensor) % coef % nlevels
        ! deallocate model profiles atmospheric arrays with RTTOV levels dimension
        call rttov_alloc_prof(allocStatus(1),1,tvs_profiles_nl(iObs),nl, &    ! 1 = nprofiles un profil a la fois
             tvs_opts(iSensor),asw=0,coefs=tvs_coefs(iSensor),init=.false. ) ! asw =0 deallocation
        call rttov_alloc_prof(allocStatus(2),1,tvs_profiles_tlad(iObs),nl, &    ! 1 = nprofiles un profil a la fois
             tvs_opts(iSensor),asw=0,coefs=tvs_coefs(iSensor),init=.false. ) ! asw =0 deallocation
        call utl_checkAllocationStatus(allocStatus(1:2), 'profiles deallocation in tvs_cleanup',.false.)
      end do

      deallocate(tvs_profiles_nl,   stat=allocStatus(1) )
      deallocate(tvs_profiles_tlad, stat=allocStatus(2) )
      call utl_checkAllocationStatus(allocStatus(1:2), ' tvs_setupAlloc tvs_profiles_nl/tlad')

      do iSensor = tvs_nsensors,1,-1
        call rttov_dealloc_coefs(allocStatus(1),  tvs_coefs(iSensor) )
        Write(*,*) 'Deallocating coefficient structure for instrument', iSensor
        call utl_checkAllocationStatus(allocStatus(1:1), ' rttov_dealloc_coefs in tvs_cleanup', .false.)
      end do

      deallocate (tvs_coefs       ,stat= allocStatus(1))
      deallocate (tvs_listSensors ,stat= allocStatus(2))
      deallocate (tvs_opts        ,stat= allocStatus(3))

      call utl_checkAllocationStatus(allocStatus(1:3), ' tvs_cleanup', .false.)

    end if

    deallocate (tvs_nchan,          stat= allocStatus(1))
    deallocate (tvs_ichan,          stat= allocStatus(2))
    deallocate (tvs_lsensor,        stat= allocStatus(3))
    deallocate (tvs_headerIndex,    stat= allocStatus(4))
    deallocate (tvs_tovsIndex,      stat= allocStatus(5))
    deallocate (tvs_isReallyPresent,stat= allocStatus(6))
    deallocate (tvs_nchanMpiGlobal, stat= allocStatus(7))
    deallocate (tvs_ichanMpiGlobal, stat= allocStatus(8))

    call utl_checkAllocationStatus(allocStatus, ' tvs_cleanup', .false.)

    Write(*,*) 'tvs_cleanup: Finished'

  end subroutine tvs_cleanup

  !--------------------------------------------------------------------------
  ! tvs_deallocateProfilesNlTlAd
  !--------------------------------------------------------------------------
  subroutine tvs_deallocateProfilesNlTlAd
    !
    ! :Purpose: release memory used by RTTOV-12.
    !
    implicit none
    integer :: allocStatus(8)

    Write(*,*) 'tvs_deallocateProfilesNlTlAd: Starting'

    allocStatus(:) = 0

    if ( radiativeTransferCode == 'RTTOV' ) then
      if ( allocated(tvs_profiles_nl) ) deallocate(tvs_profiles_nl, stat=allocStatus(1))
      if ( allocated(tvs_profiles_tlad) ) deallocate(tvs_profiles_tlad, stat=allocStatus(2))
      call utl_checkAllocationStatus(allocStatus(1:2), ' tvs_profiles_nl tvs_profiles_tlad', .false.)
    end if

    Write(*,*) 'tvs_deallocateProfilesNlTlAd: Finished'

  end subroutine tvs_deallocateProfilesNlTlAd

  !--------------------------------------------------------------------------
  ! sensors
  !--------------------------------------------------------------------------
  subroutine sensors
    !
    !:Purpose: Initialisation of the RTTOV-10 platform, satellite
    !          and instrument ID's. Also set burp to RTTOV channel
    !          mapping offset.
    !          To verify and transfom the sensor information contained in the
    !          NAMTOV namelist into the variables required by RTTTOV-7:
    !          platform, satellite and instrument ID's.
    !
    implicit none

    !Locals:
    integer sensorIndex, instrumentIndex, platformIndex
    integer ipos1, ipos2
    integer numerosat, ierr, kindex, nulnam
    character(len=15) :: tempocsatid
    logical, save :: first=.true.
    integer, save :: ioffset1b(0:ninst-1)
    character(len=8) :: listinstrum(0:ninst-1)
    character(len=15) :: tempo_inst
    integer:: listoffset(0:ninst-1)
    namelist /NAMCHANOFFSET/ listoffset, listinstrum
    integer, external :: fnom, fclos

    !  1.0 Go through sensors and set RTTOV-10 variables

    do sensorIndex=1, tvs_nsensors
      tvs_platforms  (sensorIndex) = -1
      tvs_satellites (sensorIndex) = -1
      tvs_instruments(sensorIndex) = -1
      tvs_channelOffset(sensorIndex) = -1
    end do

    if (first) then
      if ( utl_isNamelistPresent('NAMCHANOFFSET', './flnml') ) then
        call utl_abort('sensors: NAMCHANOFFSET namelist section should be now in flnml_static !')
      end if
      ! read the namelist
      nulnam = 0
      ierr = fnom(nulnam,'./flnml_static','FTN+SEQ+R/O',0)
      if (ierr /= 0) then
        write(*,*) 'Error while opening flnml_static namelist file !'
        call utl_abort('sensors')
      end if
      listoffset(:) = 0
      listinstrum(:) = 'XXXXXXXX'
      read(nulnam,NAMCHANOFFSET, iostat=ierr)
      if (ierr /= 0) then
        write(*,*) 'Error while reading NAMCHANOFFSET namelist section in flnml_static file !'
        call utl_abort('sensors')
      end if
      do instrumentIndex=0, ninst - 1
        if ( listinstrum(instrumentIndex) /= 'XXXXXXXX' ) then
          ioffset1b( tvs_getInstrumentId( listinstrum(instrumentIndex) ) )  = listoffset(instrumentIndex)
        end if
      end do
      ierr = fclos(nulnam)
      first = .false.
    end if

    !  1.1 Set platforms and satellites

    ! N.B.: Special cases for satellites TERRA and AQUA.
    !       For consistency with the RTTOV-10 nomenclature, rename:
    !       TERRA  to  eos1
    !       AQUA   to  eos2
    !       NPP    to  jpss0
    !       JPSS    to  jpss0
    !       HMWARI    to  himawari
    !       FY-3C    to  FY3-3
    do sensorIndex = 1, tvs_nsensors
      if    ( tvs_satelliteName(sensorIndex) == 'TERRA' ) then
        tempocsatid = 'eos1'
      else if ( tvs_satelliteName(sensorIndex) == 'AQUA'  ) then
        tempocsatid = 'eos2'
      else if ( tvs_satelliteName(sensorIndex) == 'NPP'  ) then
        tempocsatid = 'jpss0'
      else if ( tvs_satelliteName(sensorIndex) == 'JPSS'  ) then
        tempocsatid = 'jpss0'
      else if ( tvs_satelliteName(sensorIndex)(1:6) == 'HMWARI'  ) then
        tempocsatid = 'himawari' // trim(tvs_satelliteName(sensorIndex) (7:15))
      else if ( tvs_satelliteName(sensorIndex) == 'FY-3C'  ) then
        TEMPOCSATID = 'FY3-3'
      else
        call up2low(tvs_satelliteName(sensorIndex),tempocsatid)
      end if
      do platformIndex = 1, nplatforms
        ipos1 = len_trim(platform_name(platformIndex))
        ipos2 = index(tempocsatid,platform_name(platformIndex)(1:ipos1))
        if ( ipos2 == 1 ) then
          tvs_platforms(sensorIndex) = platformIndex
          kindex = platformIndex
        end if
      end do
      if ( tvs_platforms(sensorIndex) < 0 ) then
        write(*,'(A)') ' Satellite ' // trim(tempocsatid) // ' not supported.'
        call utl_abort('SENSORS')
      else
        ipos1 = len_trim(platform_name(kindex))
        ipos2 = len_trim(tempocsatid)
        read(tempocsatid(ipos1+1:ipos2),*,IOSTAT=ierr) numerosat
        numerosat = abs ( numerosat )
        if ( ierr /= 0) then
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
    !     RTTOV-10.

    do sensorIndex = 1, tvs_nsensors
      if ( tvs_instrumentName(sensorIndex)(1:10) == 'GOESIMAGER') then !cas particulier
        tvs_instruments(sensorIndex) = inst_id_goesim
      else if ( tvs_satelliteName(sensorIndex)(1:5) == 'MTSAT') then !autre cas particulier
        tvs_instruments(sensorIndex) = inst_id_gmsim
      else 
        call up2low(tvs_instrumentName(sensorIndex),tempo_inst)
        do instrumentIndex = 0, ninst -1 
          if ( trim(tempo_inst) == trim(inst_name(instrumentIndex))) then
            tvs_instruments(sensorIndex) = instrumentIndex
          end if
        end do
      end if
      if ( tvs_instruments(sensorIndex) < 0 ) then
        write(*,'(A)') ' INSTRUMENT '// trim( tvs_instrumentName(sensorIndex)) // ' not supported.'
        call utl_abort('SENSORS')
      end if
      tvs_channelOffset(sensorIndex) = ioffset1b(tvs_instruments(sensorIndex))
    end do

    !    1.3 Print the RTTOV-12 related variables

    if (mmpi_myid == 0) then
      write(*,*)
      write(*,'(3X,A)') '- SENSORS. Variables prepared for RTTOV-12:'
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
    integer :: idatypListSize
    integer :: idatypList(:)
    
    ! Locals:
    logical, save :: first=.true.
    integer, save :: ninst_tovs
    integer :: nulnam, ierr, instrumentIndex 
    integer, external :: fnom, fclos
    integer, save :: list_inst(ninst)
    character(len=22) :: inst_names(ninst)
    namelist /namtovsinst/ inst_names

    if (tvs_nsensors == 0) then
      ! no tovs data will be read, therefore false
      idatypList(:) = MPC_missingValue_int
      idatypListSize = 0      
      return
    end if

    if (first) then
      if ( utl_isNamelistPresent('NAMTOVSINST', './flnml') ) then
        call utl_abort('tvs_getAllIdBurpTovs: NAMTOVSINST namelist section should be now in flnml_static !')
      end if
      nulnam = 0
      ninst_tovs = 0
      list_inst(:) = -1
      inst_names(:) = 'XXXXXX'
      ierr = fnom(nulnam,'./flnml_static','FTN+SEQ+R/O',0)
      read(nulnam, nml=namtovsinst, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_getAllIdBurpTovs: Error reading NAMTOVSINST namelist section in flnml_static file')
      if (mmpi_myid == 0) write(*,nml=namtovsinst)
      ierr = fclos(nulnam)
      do instrumentIndex=1, ninst
        if (inst_names(instrumentIndex) == 'XXXXXX' ) then
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
      if ( ninst_tovs == 0 ) call utl_abort('tvs_getAllIdBurpTovs: Empty namtovsinst namelist')
      first = .false.
    end if

    idatypList(:) = MPC_missingValue_int
    idatypList(1:ninst_tovs) = list_inst(1:ninst_tovs)
    idatypListSize = ninst_tovs

  end subroutine tvs_getAllIdBurpTovs

  !--------------------------------------------------------------------------
  !  tvs_isIdBurpTovs
  !--------------------------------------------------------------------------
  logical function tvs_isIdBurpTovs(idatyp)
    !
    ! :Purpose: Function to check if the given idatyp (a.k.a. codtyp) 
    !           corresponds to a radiance
    !
    implicit none

    ! Argument:
    integer, intent(in) :: idatyp
    
    ! Locals:
    logical, save :: first=.true.
    integer, save :: ninst_tovs
    integer :: nulnam, ierr, instrumentIndex 
    integer, external :: fnom, fclos
    integer, save :: list_inst(ninst)
    character(len=22) :: inst_names(ninst)
    namelist /namtovsinst/ inst_names

    if (tvs_nsensors == 0) then
      ! no tovs data will be read, therefore false
      tvs_isIdBurpTovs = .false.
      return
    end if

    if (first) then
       if ( utl_isNamelistPresent('NAMTOVSINST', './flnml') ) then
        call utl_abort('tvs_isIdBurpTovs: NAMTOVSINST namelist section should be now in flnml_static !')
      end if
      nulnam = 0
      ninst_tovs = 0
      list_inst(:) = -1
      inst_names(:) = 'XXXXXX'
      ierr = fnom(nulnam,'./flnml_static','FTN+SEQ+R/O',0)
      read(nulnam, nml=namtovsinst, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isIdBurpTovs: Error reading NAMTOVSINST namelist section in flnml_static file')
      if (mmpi_myid == 0) write(*,nml=namtovsinst)
      ierr = fclos(nulnam)
      do instrumentIndex=1, ninst
        if (inst_names(instrumentIndex) == 'XXXXXX' ) then
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
      first = .false.
    end if

    tvs_isIdBurpTovs = .false.

    do instrumentIndex = 1, ninst_tovs
      if (idatyp == list_inst(instrumentIndex) ) then
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
    integer, intent(in) :: idatyp
    
    ! Locals:
    logical, save :: first=.true.
    integer, save :: ninst_hyper
    integer :: nulnam, ierr, instrumentIndex 
    integer, external :: fnom, fclos
    integer, save :: list_inst(ninst)
    character(len=22) :: name_inst(ninst)
    namelist /namhyper/ name_inst

    if (tvs_nsensors == 0) then
      ! no tovs data will be read, therefore false
      tvs_isIdBurpHyperSpectral = .false.
      return
    end if

    if (first) then
      if ( utl_isNamelistPresent('NAMHYPER', './flnml') ) then
        call utl_abort('tvs_isIdBurpHyperSpectral: NAMHYPER namelist section should be now in flnml_static !')
      end if
      nulnam = 0
      ninst_hyper = 0
      list_inst(:) = -1
      name_inst(:) = 'XXXXXX'
      ierr = fnom(nulnam,'./flnml_static','FTN+SEQ+R/O',0)
      read(nulnam, nml=namhyper, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isIdBurpHyperSpectral: Error reading NAMHYPER namelist section in flnml_static file')
      if (mmpi_myid == 0) write(*,nml=namhyper)
      ierr = fclos(nulnam)
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
      first = .false.
    end if

    tvs_isIdBurpHyperSpectral = .false.

    do instrumentIndex = 1, ninst_hyper
      if (idatyp == list_inst(instrumentIndex) ) then
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
    integer,          intent(in) :: idburp ! Input codtyp
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
    implicit none

    !Argument:
    character(len=*), intent(in) :: name !Platform name
    !Locals:
    integer           :: platformIndex, length, ipos
    character(len=64) :: tempo_name

    tvs_getPlatformId = -1
    length = len_trim(name)
    call up2low(name(1:length),tempo_name(1:length))

    if ( index(tempo_name(1:length),'npp') /= 0 ) then
      tvs_getPlatformId = platform_id_jpss
    else if ( index(tempo_name(1:length),'hmwari') /= 0 ) then
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

    !Argument:
    character(len=*), intent(in) :: name ! Instrument name

    !Locals:
    integer           :: instrumentIndex, length
    character(len=64) :: tempo_name

    tvs_getInstrumentId = -1
    length = len_trim(name)
    call up2low(name(1:length),tempo_name(1:length))
    if ( trim(tempo_name(1:length)) == 'goesim' ) then
      tvs_getInstrumentId = inst_id_goesim
    else if ( trim(tempo_name(1:length)) == 'gmsim' ) then
      tvs_getInstrumentId = inst_id_gmsim
    else if ( trim(tempo_name(1:length)) == 'mtsatim' ) then
      tvs_getInstrumentId = inst_id_mtsatim
    else
      do instrumentIndex = 0, ninst - 1
        if (trim(inst_name(instrumentIndex)) == tempo_name(1:length) ) then
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

    ! Argument:
    integer, intent(in) :: instrum     ! input Rttov instrument code

    ! Locals:
    integer, parameter :: maxsize = 100
    integer :: nulnam, ierr, instrumentIndex 
    integer, save :: list_inst(maxsize), ninst_hir
    logical, save :: first = .true.
    integer, external :: fclos, fnom
    character (len=8) :: name_inst(maxsize)
    namelist /NAMHYPER/ name_inst
    
    if (first) then
      if ( utl_isNamelistPresent('NAMHYPER', './flnml') ) then
        call utl_abort('tvs_isInstrumHyperSpectral: NAMHYPER namelist section should be now in flnml_static !')
      end if
      nulnam = 0
      ninst_hir = 0
      name_inst(:) = 'XXXXXXX'
      ierr = fnom(nulnam,'./flnml_static','FTN+SEQ+R/O',0)
      read(nulnam,nml=namhyper, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isInstrumHyperSpectral: Error reading namelist section NAMHYPER in flnm_static file')
      if (mmpi_myid == 0) write(*,nml=namhyper)
      ierr = fclos(nulnam)
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
      first = .false.
      if (ninst_hir == 0) then
        write(*,*) 'tvs_isInstrumHyperSpectral: Warning : empty namhyper namelist !'
      end if
    end if
    tvs_isInstrumHyperSpectral = .false.
    do instrumentIndex =1, ninst_hir
      if ( instrum == list_inst(instrumentIndex)) then
        tvs_isInstrumHyperSpectral = .true.
        exit
      end if
    end do

  end function tvs_isInstrumHyperSpectral

  !--------------------------------------------------------------------------
  !  tvs_isNameHyperSpectral
  !--------------------------------------------------------------------------
  logical function tvs_isNameHyperSpectral(cinstrum)
    !
    ! :Purpose: given an instrument name
    !           returns if it is an hyperspectral one
    !           (information from namelist NAMHYPER)
    implicit none
    !Arguments:
    character(len=*), intent(in) :: cinstrum
    !Locals:
    integer, parameter :: maxsize = 20
    integer :: nulnam, ierr, i 
    integer, save :: ninst_hir
    logical, save :: lfirst = .true.
    integer, external :: fclos, fnom
    character (len=8),save  :: name_inst(maxsize)
    character (len=8) :: name2
    namelist /NAMHYPER/ name_inst

    if (lfirst) then
      nulnam = 0
      ninst_hir = 0
      name_inst(:) = 'XXXXXXX'
      ierr = fnom(nulnam,'./flnml_static','FTN+SEQ+R/O',0)
      read(nulnam,nml=namhyper, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isNameHyperSpectral: Error reading NAMHYPER namelist section in flnml_static file')
      if (mmpi_myid == 0) write(*,nml=namhyper)
      ierr = fclos(nulnam)
      do i=1, maxsize
        if (name_inst(i) == 'XXXXXXX') then
          ninst_hir = i -1
          exit
        end if
      end do
      lfirst = .false.
      if (ninst_hir == 0) then
        write(*,*) 'tvs_isNameHyperSpectral: Warning : empty namhyper namelist !'
      end if
    end if

    tvs_isNameHyperSpectral = .false.

    call up2low(cinstrum, name2)

    do i=1, ninst_hir
      if ( trim(name2) == name_inst(i)) then
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

    ! Argument:
    integer, intent(in) :: instrum ! input Rttov instrument code

    ! Locals:
    integer, parameter :: maxsize = 100
    integer :: nulnam, ierr, instrumentIndex 
    integer, save :: list_inst(maxsize), ninst_geo
    logical, save :: first = .true.
    character(len=8) :: name_inst(maxsize)
    integer, external :: fnom, fclos

    namelist /NAMGEO/ name_inst
    if (first) then
      if ( utl_isNamelistPresent('NAMGEO', './flnml') ) then
        call utl_abort('tvs_isInstrumGeostationary: NAMGEO namelist section should be now in flnml_static !')
      end if
      nulnam = 0
      ninst_geo = 0
      name_inst(:) = 'XXXXXX'
      ierr = fnom(nulnam,'./flnml_static','FTN+SEQ+R/O',0)
      read(nulnam,nml=namgeo, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isInstrumGeostationary: Error reading namelist section NAMGEO in flnml_static file')
      if (mmpi_myid == 0) write(*,nml=namgeo)
      ierr = fclos(nulnam)
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
      first = .false.
      if (ninst_geo == 0) then
        write(*,*) 'tvs_isInstrumGeostationary: Warning : empty namgeo namelist !'
      end if
    end if
    tvs_isInstrumGeostationary = .false.
    do instrumentIndex=1, ninst_geo
      if ( instrum == list_inst(instrumentIndex)) then
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

    ! Argument:
    integer, intent(in) :: instrumId     ! input Rttov instrument code
    logical             :: idExist

    ! Locals:
    integer :: instrumentIndex 

    idExist = .false.
    do instrumentIndex = 1, tvs_numMWInstrumUsingCLW
      if ( instrumId == instrumentIdsUsingCLW(instrumentIndex) ) then
        idExist = .true.
        exit
      end if
    end do

  end function tvs_isInstrumUsingCLW

  !--------------------------------------------------------------------------
  !  tvs_isInstrumUsingHydrometeors
  !--------------------------------------------------------------------------
  function tvs_isInstrumUsingHydrometeors(instrumId) result(idExist)
    !
    ! :Purpose: given an RTTOV instrument code return if it is in the list to use Hydrometeors
    !
    implicit none

    ! Argument:
    integer, intent(in) :: instrumId     ! input Rttov instrument code
    logical             :: idExist

    ! Locals:
    integer :: instrumentIndex 

    idExist = .false.
    do instrumentIndex = 1, tvs_numMWInstrumUsingHydrometeors
      if ( instrumId == instrumentIdsUsingHydrometeors(instrumentIndex) ) then
        idExist = .true.
        exit
      end if
    end do

  end function tvs_isInstrumUsingHydrometeors

  !--------------------------------------------------------------------------
  !  tvs_isInstrumAllskyTtAssim
  !--------------------------------------------------------------------------
  function tvs_isInstrumAllskyTtAssim(instrumId) result(allskyTtAssim)
    !
    ! :Purpose: determine if all-sky temperature-channel assimilation is asked for the instrument.
    !
    implicit none

    ! Argument:
    integer, intent(in) :: instrumId     ! input Rttov instrument code
    logical             :: allskyTtAssim

    allskyTtAssim = (tvs_mwAllskyAssim .and. tvs_isInstrumUsingCLW(instrumId) .and. &
                     .not. tvs_isInstrumUsingHydrometeors(instrumId))

  end function tvs_isInstrumAllskyTtAssim

  !--------------------------------------------------------------------------
  !  tvs_isInstrumAllskyHuAssim
  !--------------------------------------------------------------------------
  function tvs_isInstrumAllskyHuAssim(instrumId) result(allskyHuAssim)
    !
    ! :Purpose: determine if all-sky humidity-channel assimilation is asked for the instrument.
    !
    implicit none

    ! Argument:
    integer, intent(in) :: instrumId     ! input Rttov instrument code
    logical             :: allskyHuAssim

    allskyHuAssim = (tvs_mwAllskyAssim .and. tvs_isInstrumUsingHydrometeors(instrumId) .and. &
                     .not. tvs_isInstrumUsingCLW(instrumId))

  end function tvs_isInstrumAllskyHuAssim

  !--------------------------------------------------------------------------
  !  tvs_isInstrumAllskyTtHuAssim
  !--------------------------------------------------------------------------
  function tvs_isInstrumAllskyTtHuAssim(instrumId) result(AllskyTtHuAssim)
    !
    ! :Purpose: determine if all-sky temperature- AND humidity-channel assimilation is asked for the instrument.
    !
    implicit none

    ! Argument:
    integer, intent(in) :: instrumId     ! input Rttov instrument code
    logical             :: AllskyTtHuAssim

    AllskyTtHuAssim = (tvs_mwAllskyAssim .and. tvs_isInstrumUsingCLW(instrumId) .and. &
                       tvs_isInstrumUsingHydrometeors(instrumId))

  end function tvs_isInstrumAllskyTtHuAssim

  !--------------------------------------------------------------------------
  !  tvs_mapInstrum
  !--------------------------------------------------------------------------
  subroutine tvs_mapInstrum(instrumburp,instrum)
    !
    ! :Purpose:  Map burp satellite instrument (element #2019) to RTTOV-7 instrument.
    !            A negative value is returned, if no match in found.
    !
    ! :Table of  RTTOV-7 instrument identifier:
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
    implicit none
    !Arguments:
    integer, intent(in)  :: instrumburp  ! burp satellite instrument (element #2019)
    integer, intent(out) :: instrum      ! RTTOV-7 instrument ID numbers (e.g. 3 for  AMSUA)
  
    !Locals:  
    integer instrumentIndex, numinstburp
    integer, parameter :: mxinstrumburp   = 100
    integer, save ::   listburp(mxinstrumburp)
    character(len=8), save :: listinstrum(mxinstrumburp)
    namelist /NAMINST/ listburp, listinstrum
    logical, save :: first = .true.
    integer :: nulnam, ier
    integer, external :: fnom, fclos

    !      Table of BURP satellite sensor identifier element #002019

    !   1.0 Find instrument

    if (first) then
      if ( utl_isNamelistPresent('NAMINST', './flnml') ) then
        call utl_abort('tvs_mapInstrum: NAMINST namelist section should be now in flnml_static !')
      end if
      
      ! set the default values
      listburp(:) = -1
      listinstrum(:) = 'XXXXXXXX'

      ! read the namelist
      nulnam = 0
      ier = fnom(nulnam,'./flnml_static','FTN+SEQ+R/O',0)
      if (ier /= 0) then
        write(*,*) 'Error while opening flnml_static namelist file !'
        call utl_abort('tvs_mapInstrum')
      end if
      read(nulnam,NAMINST,iostat=ier)
      if (ier /= 0) then
        write(*,*) 'Error while reading NAMINST namelist section in flnml_static file !'
        call utl_abort('tvs_mapInstrum')
      end if
      ier = fclos(nulnam)

      ! figure out how many valid elements in the lists
      do instrumentIndex=1, mxinstrumburp
        if(listburp(instrumentIndex) == -1) then
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
      first=.false.
    end if

    instrum = -1
    do instrumentIndex=1, mxinstrumburp
      if ( instrumburp == listburp(instrumentIndex) ) then
        instrum = tvs_getInstrumentId( listinstrum(instrumentIndex) )
        exit
      end if
    end do

  end subroutine tvs_mapInstrum

  !--------------------------------------------------------------------------
  !  tvs_isNameGeostationary
  !--------------------------------------------------------------------------
  logical function tvs_isNameGeostationary(cinstrum)
    ! :Purpose: given an instrument name following BUFR convention
    !           returns if it is a Geostationnary Imager
    !           (information from namelist NAMGEOBUFR)
    implicit none
    !Arguments:
    character(len=*), intent(in) :: cinstrum
    !Locals:
    integer, parameter :: maxsize = 100
    integer :: nulnam, ierr, i 
    integer, save :: ninst_geo
    logical, save :: lfirst = .true.
    character (len=8),save :: name_inst(maxsize)
    integer, external :: fnom, fclos

    namelist /NAMGEOBUFR/ name_inst
    if (lfirst) then
      if ( utl_isNamelistPresent('NAMGEOBUFR', './flnml') ) then
        call utl_abort('tvs_isNameGeostationary: NAMGEOBUFR namelist section should be now in flnml_static !')
      end if
      nulnam = 0
      ninst_geo = 0
      name_inst(:) = 'XXXXXXXX'
      ierr = fnom(nulnam,'./flnml_static','FTN+SEQ+R/O',0)
      read(nulnam,nml=namgeobufr, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isNameGeostationary: Error reading namelist section NAMGEOBUFR in flnml_static_file')
      if (mmpi_myid == 0) write(*,nml=namgeobufr)
      ierr = fclos(nulnam)
      do i=1, maxsize
        if (name_inst(i) == 'XXXXXXXX') then
          ninst_geo = i - 1
          exit
        end if
      end do
      lfirst = .false.
      if (ninst_geo == 0) then
        write(*,*) 'tvs_isNameGeostationary: Warning : empty namgeobufr namelist !' 
      end if
    end if
    
    tvs_isNameGeostationary = .false.
    do i=1, ninst_geo
      if ( trim(cinstrum) == trim(name_inst(i)) ) then
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
    ! :Purpose:  Map burp satellite identifier (element #1007)
    !            to RTTOV-7 platform and satellite.
    !            Negative values are returned, if no match in found.
    !
    ! :Table of  RTTOV-7 platform identifier:
    !
    ! ========          ===========================
    ! Platform          RTTOV-7 platform identifier
    ! ========          ===========================
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
    !          NOAA15, which has a burp satellite identifier value of 206,
    !          is mapped into the following:
    !          RTTOV-7 platform  =  1,
    !          RTTOV-7 satellite = 15.
    !
    ! :Arguments:
    !     :isatBURP: BURP satellite identifier
    !     :iplatform: RTTOV-7 platform ID numbers (e.g. 1 for  NOAA)
    !     :isat: RTTOV-7 satellite ID numbers (e.g. 15)
    !

    implicit none
    
    ! Arguments:
    integer, intent(in)  :: isatburp   ! BURP satellite identifier
    integer, intent(out) :: iplatform  ! RTTOV-7 platform ID numbers (e.g. 1 for  NOAA)
    integer, intent(out) :: isat       ! RTTOV-7 satellite ID numbers (e.g. 15)

    ! Locals:
    integer           :: satelliteIndex, ierr, nulnam
    logical, save     :: first=.true.
    integer, external :: fnom, fclos
    integer, parameter:: mxsatburp = 100
    integer, save     :: numsatburp
    integer, save     :: listburp(mxsatburp)         ! Table of BURP satellite identifier element #001007
    character(len=8), save :: listplat(mxsatburp)! Table of RTTOV platform identifier
    integer, save :: listsat (mxsatburp)         ! Table of RTTOV satellite identifier

    namelist /NAMSAT/ listburp, listplat, listsat

    !     Fill tables from namelist at the first call 
    if (first) then
      if ( utl_isNamelistPresent('NAMSAT', './flnml') ) then
        call utl_abort('tvs_mapSat: NAMSAT namelist section should be now in flnml_static !')
      end if
      ! set the default values
      listburp(:) = -1
      listsat(:) = -1
      listplat(:) = 'XXXXXXXX'
      ! read the namelist
      nulnam = 0
      ierr = fnom(nulnam,'./flnml_static','FTN+SEQ+R/O',0)
      if (ierr /= 0) then
        write(*,*) 'Error while opening namelist flnml_static file !'
        call utl_abort('tvs_mapSat')
      end if
      read(nulnam, NAMSAT, iostat = ierr)
      if (ierr /= 0) then
        write(*,*) 'Error while reading NAMSAT namelist section in flnml_static file !'
        call utl_abort('tvs_mapSat')
      end if
      ierr = fclos(nulnam)

      !  Figure out how many valid elements in the lists
      do satelliteIndex=1, mxsatburp
        if(listburp(satelliteIndex) == -1) then
          numsatburp = satelliteIndex - 1
          exit
        end if
      end do
      if(numsatburp >= mxsatburp) then
        call utl_abort('tvs_mapSat: exceeded maximum number of platforms')
      end if
      write(*,*) 'tvs_mapSat: number of satellites found in namelist = ',numsatburp
      write(*,*) 'tvs_mapSat: listburp   = ',listburp(1:numsatburp)
      write(*,*) 'tvs_mapSat: listsat    = ',listsat(1:numsatburp)
      write(*,*) 'tvs_mapSat: listplat   = ',listplat(1:numsatburp)
      first = .false.
    end if

    !   Find platform and satellite

    iplatform = -1
    isat      = -1
    do satelliteIndex=1, numsatburp
      if ( isatburp == listburp(satelliteIndex) ) then
        iplatform = tvs_getPlatformId( listplat(satelliteIndex) )
        isat = listsat (satelliteIndex)
        exit
      end if
    end do

  end subroutine tvs_mapSat

  !--------------------------------------------------------------------------
  !  tvs_getChanProf
  !--------------------------------------------------------------------------
  subroutine tvs_getChanprof(sensorTovsIndexes, obsSpaceData, chanprof, lchannel_subset_opt, iptobs_cma_opt)
    ! 
    ! :Purpose: subroutine to initialize the chanprof structure used by RTTOV
    !
    implicit none

    ! Arguments:
    integer, intent(in)              :: sensorTovsIndexes(:)
    type(struct_obs), intent(in)     :: obsSpaceData
    type(rttov_chanprof), intent(out):: chanprof(:)
    logical, intent(out), optional   :: lchannel_subset_opt(:,:)
    integer, intent(out), optional   :: iptobs_cma_opt(:)

    ! Locals:
    integer :: count, profileIndex, headerIndex, istart, iend, bodyIndex, channelNumber, iobs
    integer :: ChannelIndex

    ! Build the list of channels/profiles indices
    count = 0
    if (present( lchannel_subset_opt )) lchannel_subset_opt(:,:) = .false.
         
    do profileIndex = 1, size(sensorTovsIndexes)
      iobs = sensorTovsIndexes(profileIndex)
      headerIndex = tvs_headerIndex(iobs)
      if (headerIndex > 0) then
        istart = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
        iend= obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + istart - 1
        do bodyIndex = istart, iend
          if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) then
            call tvs_getChannelNumIndexFromPPP( obsSpaceData, headerIndex, bodyIndex, &
                                                channelNumber, channelIndex )
            if (channelIndex >0) then
              count =  count + 1
              chanprof(count)%prof = profileIndex
              chanprof(count)%chan = channelIndex
              if (present(iptobs_cma_opt)) iptobs_cma_opt(count) = bodyIndex
              if (present( lchannel_subset_opt )) lchannel_subset_opt(profileIndex,channelIndex) = .true.
            else
              write(*,*) 'tvs_getChanProf: strange channel number',channelNumber
            end if
          end if
        end do
      end if
    end do
  
  end subroutine tvs_getChanprof

  !--------------------------------------------------------------------------
  !  tvs_countRadiances
  !--------------------------------------------------------------------------
  integer function tvs_countRadiances(sensorTovsIndexes, obsSpaceData, assim_flag_val_opt)
    !
    ! :Purpose: to count radiances selected for assimilation
    !
    implicit none
    integer, intent(in)          :: sensorTovsIndexes(:)
    type(struct_obs)             :: obsSpaceData
    integer, intent(in),optional :: assim_flag_val_opt
    

    integer :: profileIndex, headerIndex, istart, iend, bodyIndex, iobs, assim_flag_val

    if (present(assim_flag_val_opt)) then
      assim_flag_val = assim_flag_val_opt
    else
      assim_flag_val = 1
    end if

    tvs_countRadiances = 0
    do profileIndex = 1, size(sensorTovsIndexes)
      iobs = sensorTovsIndexes(profileIndex)
      headerIndex = tvs_headerIndex(iobs)
      if (headerIndex > 0) then
        istart = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
        iend = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + istart - 1
        do bodyIndex = istart, iend
          if(obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) tvs_countRadiances  = tvs_countRadiances + 1
        end do
      end if
    end do

  end function tvs_countRadiances

  !--------------------------------------------------------------------------
  !  tvs_ChangedStypValue(obsspacedata, headerIndex)
  !--------------------------------------------------------------------------
  integer function tvs_ChangedStypValue(obsSpaceData, headerIndex)
    !
    ! :Purpose: to obtain new STYP value given observed STYP and TTYP value
    !
    implicit none
    ! Arguments:
    integer, intent(in)          :: headerIndex
    type(struct_obs)             :: obsSpaceData
    
    ! Locals:
    integer :: terrainType
    integer :: landSea 

    terrainType = obs_headElem_i(obsSpaceData,OBS_TTYP,headerIndex)
    landSea     = obs_headElem_i(obsSpaceData,OBS_STYP,headerIndex)

    if ( terrainType ==  0 ) then
      tvs_ChangedStypValue = 2
    else
      tvs_ChangedStypValue = landSea
    end if

  end function tvs_ChangedStypValue

  !--------------------------------------------------------------------------
  !  tvs_getHIREmissivities
  !--------------------------------------------------------------------------
  subroutine tvs_getHIREmissivities(sensorTovsIndexes, obsSpaceData, surfem)
    !
    ! :Purpose: to get emissivity for Hyperspectral Infrared Sounders (AIRS, IASI, CrIS, ...)
    !
    implicit none

    !Arguments:
    integer, intent(in)          :: sensorTovsIndexes(:)
    type(struct_obs), intent(in) :: obsSpaceData
    real(8), intent(out)         :: surfem(:)

    integer :: count, profileIndex, iobs, istart, iend, bodyIndex, headerIndex

    count = 0 
    surfem(:) = 0.98d0
    do profileIndex = 1, size(sensorTovsIndexes)
      iobs = sensorTovsIndexes(profileIndex)
      headerIndex = tvs_headerIndex(iobs)
      if (headerIndex > 0 ) then
        istart = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
        iend = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + istart - 1
        do bodyIndex = istart, iend
          if(obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) then
            count = count + 1
            surfem ( count ) = obs_bodyElem_r(obsSpaceData,OBS_SEM,bodyIndex)
          end if
        end do
      end if
    end do

  end subroutine tvs_getHIREmissivities

  !--------------------------------------------------------------------------
  !  tvs_getOtherEmissivities
  !--------------------------------------------------------------------------
  subroutine tvs_getOtherEmissivities(chanprof, sensorTovsIndexes, sensorType, instrument, surfem, calcemis)
    !
    ! :Purpose: to get emissivity for microwave sounders and infrared geostationary imagers
    !
    implicit none

    ! Arguments:
    type(rttov_chanprof), intent(in) :: chanprof(:)
    integer, intent(in)              :: sensorTovsIndexes(:)
    integer, intent(in)              :: sensorType
    integer, intent(in)              :: instrument
    real(8), intent(out)             :: surfem(:)
    logical, intent(out)             :: calcemis(:)
    
    ! Locals:
    integer :: radiance_index, profileIndex, iobs, surfaceType

    do radiance_index = 1, size(chanprof)
      profileIndex = chanprof(radiance_index)%prof
      iobs = sensorTovsIndexes(profileIndex)
      surfaceType = tvs_profiles_nl(iobs) % skin % surftype
      if ( sensorType == sensor_id_mw ) then
        if ( surfaceType == surftype_land .or. &
             surfaceType == surftype_seaice     ) then
          calcemis(radiance_index) = .false.
          surfem (radiance_index) = 0.75d0
        else
          calcemis(radiance_index) = .true.
          surfem (radiance_index) = 0.d0
        end if
      else if ( tvs_isInstrumHyperSpectral(instrument) ) then
        calcemis(radiance_index) = .false. 
      else if ( tvs_isInstrumGeostationary(instrument) ) then
        calcemis(radiance_index) = .true.
        surfem (radiance_index) = 0.d0
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
    type(struct_columnData), intent(in) :: columnTrl    ! Column structure
    type(struct_obs),        intent(in) :: obsSpaceData ! obsSpaceData structure
    integer,                 intent(in) :: datestamp    ! CMC date stamp
    character (len=*), intent(in) :: profileType
    logical,                 intent(in) :: beSilent     ! To control verbosity

    ! Locals:
    integer :: instrum, iplatform
    integer :: nobmax
    integer :: sensorIndex, tovsIndex
    integer :: profileCount, headerIndex
    integer :: profileIndex, levelIndex
    integer :: ilowlvl_M,ilowlvl_T,nlv_M,nlv_T
    integer :: Vcode
    integer :: ierr,day,month,year,ijour,itime
    integer :: allocStatus(12)
    
    integer,external ::  omp_get_num_threads
    integer,external ::  newdate

    integer, allocatable :: sensorTovsIndexes(:)
    integer, allocatable :: sensorHeaderIndexes(:)
  
    type(struct_vco), pointer :: vco

    real(8), allocatable :: pressure (:,:)
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
    real(8), pointer :: column_ptr(:)

    if ( .not. beSilent ) write(*,*) 'tvs_fillProfiles: Starting'
  
    if (tvs_nobtov == 0) return    ! exit if there are no tovs data

    if ( tvs_numMWInstrumUsingCLW > 0 .and. .not. col_varExist(columnTrl,'LWCR') ) then
      call utl_abort('tvs_fillProfiles: if number of instrument to use CLW greater than zero, ' // &
                     'the LWCR variable must be included as an analysis variable in NAMSTATE. ')
    end if
    if (tvs_numMWInstrumUsingHydrometeors > 0 .and. &
        .not. (col_varExist(columnTrl,'LWCR') .and. col_varExist(columnTrl,'IWCR') .and. &
               col_varExist(columnTrl,'RF')   .and. col_varExist(columnTrl,'SF')   .and. &
               col_varExist(columnTrl,'CLDR'))) then
      call utl_abort('tvs_fillProfiles: if number of instrument to use hydrometeors greater than zero, ' // &
                     'the LWCR/IWCR/RF/SF/CLDR variables must be included as an analysis variable in NAMSTATE. ')
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

    if (.not. tvs_useO3Climatology .and. .not. col_varExist(columnTrl,'TO3') .and. &
        .not. col_varExist(columnTrl,'O3L') ) then
      call utl_abort('tvs_fillProfiles: if tvs_useO3Climatology is set to .false. the ozone variable ' // &
                     'must be included as an analysis variable in NAMSTATE. ')
    else if (.not.tvs_useO3Climatology) then 
      if (col_varExist(columnTrl,'TO3') ) then
        ozoneVarName = 'TO3'
      else
        ozoneVarName = 'O3L'
      end if
    end if

    allocStatus(:) = 0

    if ( profileType == 'nl' ) then
      if ( .not. allocated( tvs_profiles_nl) ) then
        allocate(tvs_profiles_nl(tvs_nobtov) , stat=allocStatus(1) )
        if (tvs_numMWInstrumUsingHydrometeors  > 0) then
          allocate(tvs_cld_profiles_nl(tvs_nobtov) , stat=allocStatus(2))
        end if
        call utl_checkAllocationStatus(allocStatus(1:2), ' tvs_fillProfiles tvs_profiles_nl')
      end if
    else if ( profileType == 'tlad' ) then
      if ( .not. allocated( tvs_profiles_tlad) ) then
        allocate(tvs_profiles_tlad(tvs_nobtov) , stat=allocStatus(1) )
        if (tvs_numMWInstrumUsingHydrometeors  > 0) then
          allocate(tvs_cld_profiles_tlad(tvs_nobtov) , stat=allocStatus(2))
        end if
        call utl_checkAllocationStatus(allocStatus(1:1), ' tvs_fillProfiles tvs_profiles_tlad')
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
    Vcode = vco%Vcode

    ierr = newdate(datestamp,ijour,itime,-3)
    if (ierr < 0) then
      write(*,*) 'Invalid datestamp ',datestamp,ijour,itime,ierr
      call utl_abort('tvs_fillProfiles')
    end if
    year= ijour / 10000
    month = mod(ijour / 100,100)
    day = mod(ijour,100)

    !  1.2   Read ozone climatology

    if (tvs_useO3Climatology) call ozo_read_climatology(datestamp)

    !     2.  Fill profiles structure
    
    ! loop over all instruments
    sensor_loop: do sensorIndex=1, tvs_nsensors

      runObsOperatorWithClw = (col_varExist(columnTrl,'LWCR') .and. tvs_numMWInstrumUsingCLW /= 0 .and. & 
                               tvs_isInstrumUsingCLW(tvs_instruments(sensorIndex)))

      runObsOperatorWithHydrometeors = (col_varExist(columnTrl,'LWCR') .and. col_varExist(columnTrl,'IWCR') .and. &
                                        col_varExist(columnTrl,'RF')   .and. col_varExist(columnTrl,'SF')   .and. &
                                        col_varExist(columnTrl,'CLDR')                                      .and. &
                                        tvs_isInstrumUsingHydrometeors(tvs_instruments(sensorIndex)))

      if (runObsOperatorWithClw .and. runObsOperatorWithHydrometeors) then
        call utl_abort('tvs_fillProfiles: this instrument is mentioned in using CLW and hydrometeors.')
      end if

      ! first loop over all obs.
      profileCount = 0
      bobs1: do tovsIndex = 1, tvs_nobtov
        if (tvs_lsensor(tovsIndex) == sensorIndex) then
          profileCount = profileCount + 1
          NOBMAX = tovsIndex
        end if
      end do bobs1

      if (profileCount == 0) cycle sensor_loop

      allocStatus(:) = 0
      allocate (sensorTovsIndexes(profileCount),                     stat = allocStatus(1) )
      allocate (sensorHeaderIndexes(profileCount),                   stat = allocStatus(2) )
      allocate (latitudes(profileCount),                             stat = allocStatus(3) )
      allocate (ozone(nlv_T,profileCount),                           stat = allocStatus(4) ) 
      allocate (pressure(nlv_T,profileCount),                        stat = allocStatus(5) )
      if (runObsOperatorWithClw .or. runObsOperatorWithHydrometeors) then
        allocate (clw       (nlv_T,profileCount),stat= allocStatus(6))
        clw(:,:) = qlim_getMinValueCloud('LWCR')
      end if
      if (runObsOperatorWithHydrometeors) then
        allocate (ciw       (nlv_T,profileCount),stat= allocStatus(7))
        allocate (rainFlux  (nlv_T,profileCount),stat= allocStatus(8))
        allocate (snowFlux  (nlv_T,profileCount),stat= allocStatus(9))
        allocate (cloudFraction(nlv_T,profileCount),stat= allocStatus(10))
        ciw(:,:) = qlim_getMinValueCloud('IWCR')
        rainFlux(:,:) = qlim_getMinValueCloud('RF')
        snowFlux(:,:) = qlim_getMinValueCloud('SF')
        cloudFraction(:,:) = qlim_getMinValueCloud('CLDR')
      end if
      allocate (surfTypeIsWater(profileCount),stat= allocStatus(11)) 
      surfTypeIsWater(:) = .false.

      call utl_checkAllocationStatus(allocStatus, ' tvs_fillProfiles')

      profileCount = 0

      ! second loop over all obs.
      bobs2: do tovsIndex = 1, NOBMAX
        if (tvs_lsensor(tovsIndex) /= sensorIndex) cycle bobs2
        profileCount = profileCount + 1
        sensorTovsIndexes(profileCount) = tovsIndex
        headerIndex = tvs_headerIndex(tovsIndex)
        sensorHeaderIndexes(profileCount) = headerIndex

        call rttov_alloc_prof(                 &
             allocStatus(1),                   &
             1,                                & ! 1 = nprofiles un profil a la fois
             profiles(tovsIndex:tovsIndex),    &
             nlv_T,                            & 
             tvs_opts(sensorIndex),            &
             asw=1,                            & ! asw =1 allocation
             coefs=tvs_coefs(sensorIndex),     &
             init=.true. )                       
        if (runObsOperatorWithHydrometeors) then
          call rttov_alloc_scatt_prof(            &   
               allocstatus(2),                    &
               1,                                 &
               cld_profiles(tovsIndex:tovsIndex), &
               nlv_T,                             &
               nhydro=5,                          & ! depending on what is defined in the Mie tables
               nhydro_frac=1,                     & ! 1 cloud fraction for all variable or nhydro 1 cloud fraction for each variable
               asw=1_jpim,                        & ! 1 => allocate
               init=.true.,                       & ! initialize profiles to zero
               flux_conversion=[1,2,0,0,0] )        !flux_conversion  input units: 0 (default) => kg/kg,
                                                    ! 1,2 => kg/m2/s, optional for rain, snow
        end if
        call utl_checkAllocationStatus(allocStatus(1:2), ' tvs_setupAlloc tvs_fillProfiles')

        !    extract land/sea/sea-ice flag (0=land, 1=sea, 2=sea-ice)
        profiles(tovsIndex) % skin % surftype = tvs_ChangedStypValue(obsSpaceData,headerIndex)

        !    extract satellite zenith and azimuth angle, 
        !    sun zenith angle, cloud fraction, latitude and longitude
        profiles(tovsIndex) % zenangle   = obs_headElem_r(obsSpaceData,OBS_SZA,headerIndex)

        !pour ne pas faire planter RTTOV dans le cas (rare) ou l'angle zenithal n'est pas defini ou invalide         
        if (profiles(tovsIndex) % zenangle < 0.0d0 .or. &
             profiles(tovsIndex) % zenangle > zenmax ) then
          write(*,*) '!!! WARNING !!!'
          write(*,*) 'INVALID ZENITH ANGLE'
          write(*,*) 'angle, profile number, sensor', profiles(tovsIndex) % zenangle, tovsIndex, sensorIndex
          write(*,*) 'replaced by 0.0 !!!'
          profiles(tovsIndex) % zenangle = 0.d0
        end if
 
        profiles(tovsIndex) % azangle = tvs_getCorrectedSatelliteAzimuth(obsSpaceData, headerIndex)
        profiles(tovsIndex) % sunazangle  = obs_headElem_r(obsSpaceData,OBS_SAZ,headerIndex) ! necessaire pour radiation solaire
        iplatform = tvs_coefs(sensorIndex) % coef % id_platform
        instrum = tvs_coefs(sensorIndex) % coef % id_inst
        
        profiles(tovsIndex) % sunzenangle = obs_headElem_r(obsSpaceData,OBS_SUN,headerIndex)
        latitudes(profileCount) = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex) *MPC_DEGREES_PER_RADIAN_R8
        profiles(tovsIndex) % longitude =  obs_headElem_r(obsSpaceData,OBS_LON,headerIndex) *MPC_DEGREES_PER_RADIAN_R8

        surfTypeIsWater(profileCount) = ( tvs_ChangedStypValue(obsSpaceData,headerIndex) == surftype_sea )

        do levelIndex = 1, nlv_T
          pressure(levelIndex,profileCount) = col_getPressure(columnTrl,levelIndex,headerIndex,'TH') * MPC_MBAR_PER_PA_R8
          if ((runObsOperatorWithClw .and. surfTypeIsWater(profileCount)) .or. &
              (runObsOperatorWithHydrometeors .and. surfTypeIsWater(profileCount))) then
            clw(levelIndex,profileCount) = col_getElem(columnTrl,levelIndex,headerIndex,'LWCR')

            if (clw(levelIndex,profileCount) < qlim_getMinValueCloud('LWCR') .or. &
                clw(levelIndex,profileCount) > qlim_getMaxValueCloud('LWCR')) then
              write(*,*) 'tvs_fillProfiles: clw=' , clw(:,profileCount) 
              call utl_abort('tvs_fillProfiles: columnTrl has clw outside RTTOV bounds')
            end if

            clw(levelIndex,profileCount) = clw(levelIndex,profileCount) * tvs_cloudScaleFactor
          end if
          if (runObsOperatorWithHydrometeors .and. surfTypeIsWater(profileCount)) then
            ciw(levelIndex,profileCount) = col_getElem(columnTrl,levelIndex,headerIndex,'IWCR')
            rainFlux(levelIndex,profileCount) = col_getElem(columnTrl,levelIndex,headerIndex,'RF')
            snowFlux(levelIndex,profileCount) = col_getElem(columnTrl,levelIndex,headerIndex,'SF')
            cloudFraction(levelIndex,profileCount) = col_getElem(columnTrl,levelIndex,headerIndex,'CLDR')

            if (ciw(levelIndex,profileCount) < qlim_getMinValueCloud('IWCR') .or. &
                ciw(levelIndex,profileCount) > qlim_getMaxValueCloud('IWCR')) then
              write(*,*) 'tvs_fillProfiles: ciw=' , ciw(:,profileCount) 
              call utl_abort('tvs_fillProfiles: columnTrl has ciw outside RTTOV bounds')
            end if

            if (rainFlux(levelIndex,profileCount) < qlim_getMinValueCloud('RF') .or. &
                rainFlux(levelIndex,profileCount) > qlim_getMaxValueCloud('RF')) then
              write(*,*) 'tvs_fillProfiles: rainFlux=' , rainFlux(:,profileCount) 
              call utl_abort('tvs_fillProfiles: columnTrl has rain flux outside RTTOV bounds')
            end if

            if (snowFlux(levelIndex,profileCount) < qlim_getMinValueCloud('SF') .or. &
                snowFlux(levelIndex,profileCount) > qlim_getMaxValueCloud('SF')) then
              write(*,*) 'tvs_fillProfiles: snowFlux=' , snowFlux(:,profileCount) 
              call utl_abort('tvs_fillProfiles: columnTrl has snow flux outside RTTOV bounds')
            end if

            if (cloudFraction(levelIndex,profileCount) < qlim_getMinValueCloud('CLDR') .or. &
                cloudFraction(levelIndex,profileCount) > qlim_getMaxValueCloud('CLDR')) then
              write(*,*) 'tvs_fillProfiles: cloudFraction=' , cloudFraction(:,profileCount) 
              call utl_abort('tvs_fillProfiles: columnTrl has cloud fraction outside RTTOV bounds')
            end if

            ciw(levelIndex,profileCount) = ciw(levelIndex,profileCount) * tvs_cloudScaleFactor
            rainFlux(levelIndex,profileCount) = rainFlux(levelIndex,profileCount) * tvs_cloudScaleFactor
            snowFlux(levelIndex,profileCount) = snowFlux(levelIndex,profileCount) * tvs_cloudScaleFactor
            !should we introduce a separate scaling factor for various variables ? 
          end if
        end do
        
        if (tvs_coefs(sensorIndex) %coef % nozone > 0 .and. .not. tvs_useO3Climatology) then
          column_ptr => col_getColumn(columnTrl, headerIndex, trim(ozoneVarName) )
          ! Conversion from microgram/km to ppmv (to have the same units as climatology when tvs_useO3Climatology is .true.
          ! Conversion to kg/kg for use by RTTOV in done later
          ozone(:,profileCount) = column_ptr(:) * 1.0D-9 * o3Mixratio2ppmv
        end if

      end do bobs2

      !    2.5  Get ozone profiles (ppmv) from climatology if necessary
      if (tvs_coefs(sensorIndex) %coef % nozone > 0 .and. tvs_useO3Climatology) then
        call ozo_get_profile (ozone, latitudes, pressure, nlv_T, profileCount)
      end if

      !   2.5  Fill profiles structure

      do  profileIndex = 1 , profileCount 
        tovsIndex = sensorTovsIndexes(profileIndex)
        headerIndex = sensorHeaderIndexes(profileIndex)
        profiles(tovsIndex) % gas_units       = gas_unit_specconc ! all gas profiles are supposed to be provided in kg/kg (specific humidity, i.e. mass mixing ratio [kg/kg] over wet air)
        profiles(tovsIndex) % id              = '' ! profile id, up to 128 characters, to consider for use
        profiles(tovsIndex) % nlevels         = nlv_T
        profiles(tovsIndex) % nlayers         = nlv_T - 1
        profiles(tovsIndex) % date(1)         = year
        profiles(tovsIndex) % date(2)         = month
        profiles(tovsIndex) % date(3)         = day
        profiles(tovsIndex) % latitude        = latitudes(profileIndex)
        profiles(tovsIndex) % elevation       = 0.001d0 * col_getHeight(columnTrl,ilowlvl_T,headerIndex,'TH') ! unite km
        profiles(tovsIndex) % skin % watertype= 1 !utilise pour calcul rayonnement solaire reflechi seulement
        profiles(tovsIndex) % skin % t        = col_getElem(columnTrl,1,headerIndex,'TG')
        profiles(tovsIndex) % skin % salinity = 35.d0 ! for FASTEM-4 only to revise (practical salinity units)
        profiles(tovsIndex) % skin % fastem(:)= 0.0d0
        profiles(tovsIndex) % skin % snow_fraction  = 0.d0 ! Surface coverage snow fraction(0-1), used only by IR emissivity atlas
        profiles(tovsIndex) % skin % soil_moisture  = 0.d0 ! soil moisure (m**3/m**3) not yet used
        profiles(tovsIndex) % s2m % t         = col_getElem(columnTrl,ilowlvl_T,headerIndex,'TT')
        profiles(tovsIndex) % s2m % q         = 0.3D6  * qppmv2Mixratio ! a value between 0 and 0.6d6 so that RTTOV will not complain; not used
        profiles(tovsIndex) % s2m % p         = col_getElem(columnTrl,1      ,headerIndex,'P0')*MPC_MBAR_PER_PA_R8
        profiles(tovsIndex) % s2m % u         = col_getElem(columnTrl,ilowlvl_M,headerIndex,'UU')
        profiles(tovsIndex) % s2m % v         = col_getElem(columnTrl,ilowlvl_M,headerIndex,'VV')
        profiles(tovsIndex) % s2m % o         = 0.0d0 !surface ozone never used
        profiles(tovsIndex) % s2m % wfetc     = 100000.0d0 ! Wind fetch (in meter for rttov10 ?) used to calculate reflection of solar radiation by sea surface
        profiles(tovsIndex) % icede_param     = 0
        profiles(tovsIndex) % Be              = 0.4d0 ! earth magnetic field strength (gauss) (must be non zero)
        profiles(tovsIndex) % cosbk           = 0.0d0 ! cosine of the angle between the earth magnetic field and wave propagation direction
        profiles(tovsIndex) % p(:)            = pressure(:,profileIndex)
        !RTTOV scatt needs half pressure levels (see figure 5 of RTTOV 12 User's Guide)
        if (runObsOperatorWithHydrometeors) then
          cld_profiles(tovsIndex) % ph (1) = 0.d0
          cld_profiles(tovsIndex) % cfrac = 0.d0
          do levelIndex = 1, nlv_T - 1
            cld_profiles(tovsIndex) % ph (levelIndex+1) = 0.5d0 * (profiles(tovsIndex) % p(levelIndex) + &
                                                                   profiles(tovsIndex) % p(levelIndex+1))
          end do
          cld_profiles(tovsIndex) % ph (nlv_T+1) = profiles(tovsIndex) % s2m % p
        end if
        column_ptr => col_getColumn(columnTrl, headerIndex,'TT' )
        profiles(tovsIndex) % t(:)   = column_ptr(:)
        
        if (tvs_coefs(sensorIndex) %coef %nozone > 0) then
          profiles(tovsIndex) % o3(:) = ozone(:,profileIndex) * o3ppmv2Mixratio ! Climatology output is ppmv over dry air                                                                                                            ! because atmosphere is very dry where there is significant absorption by ozone)
          if (.not. tvs_useO3Climatology)  then
            profiles(tovsIndex) % s2m % o  = col_getElem(columnTrl,ilowlvl_T,headerIndex,trim(ozoneVarName)) * 1.0d-9 ! Assumes model ozone in ug/kg
          end if
        end if

        column_ptr => col_getColumn(columnTrl, headerIndex,'HU' )
        profiles(tovsIndex) % q(:)            =  column_ptr(:)

        profiles(tovsIndex) % ctp = 1013.25d0
        profiles(tovsIndex) % cfraction = 0.d0
        if (runObsOperatorWithClw) then
          profiles(tovsIndex) % clw(:) = clw(:,profileIndex)
        else if (runObsOperatorWithHydrometeors) then
          cld_profiles(tovsIndex) % hydro(:,1) = rainFlux(:,profileIndex)
          cld_profiles(tovsIndex) % hydro(:,2) = snowFlux(:,profileIndex)
          cld_profiles(tovsIndex) % hydro(:,4) = clw(:,profileIndex)
          cld_profiles(tovsIndex) % hydro(:,5) = ciw(:,profileIndex)
          
          where (cld_profiles(tovsIndex) % hydro(:,1) > qlim_getMinValueCloud('RF') .or. &
                 cld_profiles(tovsIndex) % hydro(:,2) > qlim_getMinValueCloud('SF') .or. &
                 cld_profiles(tovsIndex) % hydro(:,4) > qlim_getMinValueCloud('LWCR') .or. &
                 cld_profiles(tovsIndex) % hydro(:,5) > qlim_getMinValueCloud('IWCR'))
                 cld_profiles(tovsIndex) % hydro_frac(:,1) = cloudFraction(:,profileIndex)
          elsewhere
            cld_profiles(tovsIndex) % hydro_frac(:,1) = qlim_getMinValueCloud('CLDR')
          end where
        end if     
      end do

      deallocate (pressure,            stat = allocStatus(2))
      deallocate (ozone,               stat = allocStatus(3))
      deallocate (latitudes,           stat = allocStatus(4))
      deallocate (sensorHeaderIndexes, stat = allocStatus(5))
      deallocate (sensorTovsIndexes,   stat = allocStatus(6))
      if (tvs_coefs(sensorIndex) %coef %nozone > 0 .and. .not.tvs_useO3Climatology) then
        deallocate (ozone,             stat = allocStatus(7))
      end if
      if ( allocated(clw) ) then
        deallocate (clw,      stat= allocStatus(8))
      end if
      if ( allocated(ciw) ) then
        deallocate (ciw,      stat= allocStatus(9))
        deallocate (rainFlux, stat= allocStatus(10))
        deallocate (snowFlux, stat= allocStatus(11))
        deallocate (cloudFraction, stat= allocStatus(12))
      end if
      deallocate (surfTypeIsWater,stat= allocStatus(13))

      call utl_checkAllocationStatus(allocStatus, ' tvs_fillProfiles', .false.)
     
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
    ! Arguments
    type(struct_obs), intent(in) :: obsSpaceData     ! obsSpaceData structure
    integer, intent(in)          :: headerIndex      ! location in header
    real(8)                      :: correctedAzimuth ! corrected azimuth (function result)
    ! Locals
    integer :: sensorNo, tovsIndex

    correctedAzimuth = obs_headElem_r(obsSpaceData,OBS_AZA,headerIndex)

    tovsIndex = tvs_tovsIndex (headerIndex)
    if ( tovsIndex < 0) return

    sensorNo  = tvs_lsensor(tovsIndex)

    if ( .not. tvs_isAzimuthValid(sensorNo) ) then
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

    ! Arguments
    type(struct_obs), intent(inout) :: obsSpaceData  ! obsSpaceData structure
    logical, intent(in)             :: bgckMode       ! flag to transfer transmittances and cloudy overcast radiances in bgck mode 
    logical, intent(in)             :: beSilent       ! flag to control verbosity

    ! Locals:
    integer :: nlv_T
    integer :: btCount
    integer :: allocStatus(4)
    integer :: rttov_err_stat ! rttov error return code
    integer, external :: omp_get_num_threads
    integer :: nthreads,max_nthreads
    integer :: sensorId, tovsIndex
    integer :: channelIndex, channelIndexFound, channelNumber
    integer :: profileCount
    integer :: profileIndex, levelIndex, jj, btIndex
    integer :: instrum
    integer :: sensorType        ! sensor type (1=infrared; 2=microwave; 3=high resolution,4=polarimetric)
    integer, allocatable:: sensorTovsIndexes(:)  
    
    type (rttov_emissivity), pointer :: emissivity_local(:)    ! emissivity structure with input and output
    type (rttov_chanprof), pointer :: chanprof(:)
    type (rttov_chanprof), allocatable :: chanprof1(:)
    type (rttov_radiance) :: radiancedata_d, radiancedata_d1
    type (rttov_transmission) :: transmission, transmission1
    integer              :: asw
    logical, pointer :: calcemis(:)
    real(8), allocatable  :: surfem1(:)
    integer, allocatable  :: frequencies(:)
    real(8), allocatable  :: uOfWLandWSurfaceEmissivity(:)
    logical, allocatable  :: lchannel_subset(:,:)
    integer               :: profileIndex2, tb1, tb2
    integer :: istart, iend, bodyIndex, headerIndex
    real(8) :: clearMwRadiance
    logical :: ifBodyIndexFound
    logical :: runObsOperatorWithClw
    logical :: runObsOperatorWithHydrometeors

    if ( .not. beSilent ) write(*,*) 'tvs_rttov: Starting'
    if ( .not. beSilent ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if (tvs_nobtov == 0) return                  ! exit if there are not tovs data

    !   1.  Get number of threads available and allocate memory for some variables
    !$omp parallel
    max_nthreads = omp_get_num_threads()
    !$omp end parallel
    allocStatus(:) = 0
    allocate(sensorTovsIndexes(tvs_nobtov),stat=allocStatus(1))
    call utl_checkAllocationStatus(allocStatus(1:1), ' tvs_rttov sensorTovsIndexes')
    
    
    !   1.1   Read surface information
    if ( bgckMode ) call EMIS_READ_CLIMATOLOGY

    !   2.  Computation of hx for tovs data only

    ! Loop over all sensors specified by user
    sensor_loop:do sensorId = 1, tvs_nsensors
   
      sensorType = tvs_coefs(sensorId) % coef % id_sensor
      instrum = tvs_coefs(sensorId) % coef % id_inst

      runObsOperatorWithClw = (tvs_numMWInstrumUsingCLW /= 0 .and. &
                               tvs_isInstrumUsingCLW(tvs_instruments(sensorId)))
      runObsOperatorWithHydrometeors = (tvs_numMWInstrumUsingHydrometeors /= 0 .and. &
                                        tvs_isInstrumUsingHydrometeors(tvs_instruments(sensorId)))
                                        
      if (runObsOperatorWithClw .and. runObsOperatorWithHydrometeors) then
        call utl_abort('tvs_rttov: this instrument is mentioned in using CLW and hydrometeors.')
      end if
    
      !  loop over all obs.
      profileCount = 0
      obs_loop: do tovsIndex = 1, tvs_nobtov
        
        !    Currently processed sensor?
        if ( tvs_lsensor(tovsIndex) == sensorId ) then
          profileCount = profileCount + 1
          sensorTovsIndexes(profileCount) = tovsIndex
          nlv_T = tvs_profiles_nl(tovsIndex) % nlevels
        end if
      end do obs_loop
      
      if (profileCount == 0) cycle sensor_loop
      
      !    2.1  Calculate the actual number of threads which will be used.
      
      nthreads = min(max_nthreads, profileCount )  
      
      !    2.2  Prepare all input variables required by rttov.
      
      if ( bgckMode .and. tvs_isInstrumHyperSpectral(instrum) ) then
        btCount = profileCount * tvs_nchan(sensorId)
      else
        btCount = tvs_countRadiances(sensorTovsIndexes(1:profileCount), obsSpaceData)
      end if
      
      if (allocated(tvs_bodyIndexFromBtIndex)) deallocate(tvs_bodyIndexFromBtIndex)
      allocate(tvs_bodyIndexFromBtIndex(btCount))

      if ( btCount == 0 ) cycle sensor_loop
      tvs_bodyIndexFromBtIndex(:) = -1      
      allocate ( surfem1          (btCount) ,stat=allocStatus(1))

      asw = 1 ! Allocation
      call rttov_alloc_direct(            &
              allocStatus(2),             &
              asw,                        &
              nprofiles=profileCount,     & ! (not used)
              nchanprof=btCount,          &
              nlevels=nlv_T,              &
              chanprof=chanprof,          &
              opts=tvs_opts(sensorId),    &
              coefs=tvs_coefs(sensorId),  &
              transmission=transmission,  &
              radiance=radiancedata_d,    &
              calcemis=calcemis,          &
              emissivity=emissivity_local,&
              init=.true.)

      if (runObsOperatorWithHydrometeors) then
        allocate ( frequencies(btCount), stat=allocStatus(3))
      end if

      if (useUofWIREmiss) then
        allocate ( uOfWLandWSurfaceEmissivity(btCount), stat=allocStatus(4) )
      end if
      call utl_checkAllocationStatus(allocStatus, ' tvs_rttov')
      
      !     get Hyperspectral IR emissivities
      if ( tvs_isInstrumHyperSpectral(instrum) ) then
        surfem1(:) = 0.
        if ( bgckMode ) then
          call emis_getIrEmissivity (surfem1,tvs_nchan(sensorId),sensorId,profileCount,btCount,sensorTovsIndexes)
        else
          call tvs_getHIREmissivities(sensorTovsIndexes(1:profileCount), obsSpaceData, surfem1)
        end if
      end if
      
      if ( bgckMode .and. tvs_isInstrumHyperSpectral(instrum) ) then
        btIndex = 0
        do profileIndex = 1 , profileCount
          do  channelIndex = 1,tvs_nchan(sensorId)
            btIndex = btIndex + 1
            chanprof(btIndex)%prof = profileIndex
            chanprof(btIndex)%chan = channelIndex
          end do
        end do
        
        do profileIndex = 1, profileCount
          headerIndex = tvs_headerIndex(sensorTovsIndexes(profileIndex))
          if (headerIndex > 0) then
            istart = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
            iend = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + istart - 1
            do bodyIndex = istart, iend
              call tvs_getChannelNumIndexFromPPP( obsSpaceData, headerIndex, bodyIndex, &
                  channelNumber, channelIndex )
              if (channelIndex > 0) then
                tvs_bodyIndexFromBtIndex((profileIndex-1)*tvs_nchan(sensorId)+channelIndex) = bodyIndex
              else
                write(*,*) 'tvs_rttov: strange channel number',channelNumber
              end if
            end do
          end if
        end do
        
      else
        allocate( lchannel_subset(profileCount,tvs_nchan(sensorId)) )
        call tvs_getChanprof(sensorTovsIndexes(1:profileCount), obsSpaceData, chanprof, lchannel_subset_opt = lchannel_subset, iptobs_cma_opt = tvs_bodyIndexFromBtIndex)
        if (runObsOperatorWithHydrometeors) then
          call rttov_scatt_setupindex (  &
               rttov_err_stat,           &
               profileCount,             &  ! number of profiles
               tvs_nchan(sensorId),      &  ! number of channels 
               tvs_coefs(sensorId),      &  ! coef structure read in from rttov coef file
               tvs_coef_scatt(sensorId), &  ! 
               btcount,                  &  ! number of calculated channels
               chanprof,                 &  ! channels and profile numbers
               frequencies,              &  ! array, frequency number for each channel
               lchannel_subset )            ! OPTIONAL array of logical flags to indicate a subset of channels
        end if
        deallocate( lchannel_subset )
      end if
                 
      call tvs_getOtherEmissivities(chanprof, sensorTovsIndexes, sensorType, instrum, surfem1, calcemis)
      
      if (useUofWIREmiss .and. tvs_isInstrumHyperSpectral(instrum) .and. bgckMode) then
        if (.not. allocated (tvs_atlas)) allocate(tvs_atlas(tvs_nsensors))
        if ( .not. tvs_atlas(sensorId)%init) then
          call rttov_setup_emis_atlas( rttov_err_stat,  &! out
               tvs_opts(sensorId),                      &! in
               tvs_profiles_nl(1)%date(2) ,             &! in
               atlas_type_ir,                           &! in
               tvs_atlas(sensorId),                     &! in
               ir_atlas_ang_corr = .false.,             &! in
               ir_atlas_read_std = .false.,             &! in
               coefs = tvs_coefs(sensorId)  )
          if (rttov_err_stat/=0) then
            write(*,*) 'Error in rttov_atlas_setup IR',rttov_err_stat
            call utl_abort('tvs_rttov')
          end if
        end if
          
        call rttov_get_emis( rttov_err_stat,                  & ! out
             tvs_opts(sensorId),                              & ! in
             chanprof(1:btCount),                             & ! in
             tvs_profiles_nl(sensorTovsIndexes(1:profileCount)), & ! in
             tvs_coefs(sensorId),                             & ! in
             tvs_atlas(sensorId),                             & ! inout
             uOfWLandWSurfaceEmissivity(1:btCount) )            ! out

        if (rttov_err_stat /= 0) then
          write(*,*) 'Error in rttov_get_emis IR', rttov_err_stat
          call utl_abort('tvs_rttov')
        end if
              
        do profileIndex=1, profileCount !loop on profiles
          jj = sensorTovsIndexes(profileIndex)
          do btIndex=1, btCount !loop on channels
            if (chanprof(btIndex)%prof==profileIndex) then
              ! surftype: 0 land, 1 sea, 2 sea-ice
              ! this logic is primitive and could be improved for example using
              ! additional criteria based on emissivity_std and emissivity_flg
              !Definition of emis_flag:
              ! emis_flag:Flag_0 = '0 = sea, no MOD11 data' ;
              ! emis_flag:Flag_1 = '1 = land where BF method was applied' ;
              ! emis_flag:Flag_2 = '2 = land where data was filled with average (original UWiremis bfemis_flag=2 or 3 or 4' ;
              ! emis_flag:Flag_3 = '3 = contains inland water or coastline by the sea/land mask where the BF method was used' ;
              ! emis_flag:Flag_4 = '4 = contains inland water or coastline by the sea/land mask where data was filled with average original UWiremis bfemis_flag=2 or 3 or 4' ;
              ! emis_flag:Flag_5 = '5 = contains coastline by land fraction where the BF method was used' ;
              ! emis_flag:Flag_6 = '6 = contains coastline by land fraction where data was filled with average (original UWiremis bfemis_flag=2 or 3 or 4' ;
              ! other information that could be useful for quality control can be found in the in the profile_qc structure
              ! Now we have the 'traditionnal' emissivity in surfem1(:)
              ! and University of Wisconsin emissivity in uOfWLandWSurfaceEmissivity(:)
              if (tvs_profiles_nl(jj)% skin % surftype == surftype_land .and. &
                   uOfWLandWSurfaceEmissivity(btIndex) > 0.5 ) then
                emissivity_local(btIndex)%emis_in = uOfWLandWSurfaceEmissivity(btIndex)
              else
                emissivity_local(btIndex)%emis_in = surfem1(btIndex)
              end if
            end if
          end do
        end do

      else if (sensorType == sensor_id_mw) then
        call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorId, chanprof, &
             sensorTovsIndexes(1:profileCount))
      else
        emissivity_local(:)%emis_in = surfem1(:)
      end if

      !   2.3  Compute radiance with rttov_direct

      rttov_err_stat = 0 

      if( bgckMode .and. tvs_isInstrumHyperSpectral(instrum) ) then
        write(*,*) 'for bgck IR: call rttov_parallel_direct for each profile...'

        asw = 1 ! 1 to allocate,0 to deallocate
        ! allocate transmitance structure for 1 profile
        call rttov_alloc_transmission(allocStatus(1), transmission1, nlevels=nlv_T, &
             nchanprof=tvs_nchan(sensorId), asw=asw, init=.true.)
        ! allocate radiance structure for 1 profile
        call rttov_alloc_rad (allocStatus(2),tvs_nchan(sensorId), radiancedata_d1,nlv_T,asw,init=.true.)
        ! allocate chanprof for 1 profile
        allocate(chanprof1(tvs_nchan(sensorId)))
        do  channelIndex = 1,tvs_nchan(sensorId)
          chanprof1(channelIndex)%prof = 1
          chanprof1(channelIndex)%chan = channelIndex
        end do

        do profileIndex2 = 1, profileCount
          tb1 = 1 + (profileIndex2-1) * tvs_nchan(sensorId) 
          tb2 = profileIndex2 * tvs_nchan(sensorId)

          call rttov_parallel_direct(                                                              &
               rttov_err_stat,                                                                     & ! out
               chanprof1,                                                                          & ! in
               tvs_opts(sensorId),                                                                 & ! in
               tvs_profiles_nl(sensorTovsIndexes(profileIndex2):sensorTovsIndexes(profileIndex2)), & ! in
               tvs_coefs(sensorId),                                                                & ! in
               transmission1,                                                                      & ! inout
               radiancedata_d1,                                                                    & ! inout
               calcemis=calcemis(tb1:tb2),                                                         & ! in
               emissivity=emissivity_local(tb1:tb2),                                               & ! inout
               nthreads=nthreads )   

          ! copy contents of single profile structures into complete structures
          transmission%tau_total(tb1:tb2)             = transmission1%tau_total(:)
          transmission%tau_levels(:,tb1:tb2)          = transmission1%tau_levels(:,:)
          transmission%tausun_levels_path1(:,tb1:tb2) = transmission1%tausun_levels_path1(:,:)
          transmission%tausun_levels_path2(:,tb1:tb2) = transmission1%tausun_levels_path2(:,:)
          transmission%tausun_total_path1(tb1:tb2)    = transmission1%tausun_total_path1(:)
          transmission%tausun_total_path2(tb1:tb2)    = transmission1%tausun_total_path2(:)
          radiancedata_d%clear(tb1:tb2)      = radiancedata_d1%clear(:)
          radiancedata_d%total(tb1:tb2)      = radiancedata_d1%total(:)
          radiancedata_d%bt_clear(tb1:tb2)   = radiancedata_d1%bt_clear(:)
          radiancedata_d%bt(tb1:tb2)         = radiancedata_d1%bt(:)
          radiancedata_d%refl_clear(tb1:tb2) = radiancedata_d1%refl_clear(:)
          radiancedata_d%refl(tb1:tb2)       = radiancedata_d1%refl(:)
          radiancedata_d%overcast(:,tb1:tb2) = radiancedata_d1%overcast(:,:)
          radiancedata_d%cloudy(tb1:tb2)     = radiancedata_d1%cloudy(:)

        end do

        ! transmittance deallocation for 1 profile
        deallocate(chanprof1)
        asw = 0 ! 1 to allocate,0 to deallocate
        ! transmittance deallocation for 1 profile
        call rttov_alloc_transmission(allocStatus(1),transmission1,nlevels=nlv_T,  &
             nchanprof=tvs_nchan(sensorId), asw=asw )
        ! radiance deallocation for 1 profile
        call rttov_alloc_rad (allocStatus(2), tvs_nchan(sensorId), radiancedata_d1, nlv_T, asw)

      else

        ! run clear-sky RTTOV, save the radiances in OBS_BTCL of obsSpaceData 
        if ((runObsOperatorWithClw .or. runObsOperatorWithHydrometeors) .and. &
            obs_columnActive_RB(obsSpaceData, OBS_BTCL)) then

          ! run rttovScatt
          if (runObsOperatorWithHydrometeors) then
            ! set the cloud profile in tvs_cld_profiles_nl to zero
            call updateCloudInTovsCloudProfile(sensorTovsIndexes(1:profileCount), &
                                          nlv_T,                      &
                                          mode='save',                &
                                          beSilent=.true.)
            call rttov_scatt(                                         &
                 rttov_err_stat,                                      &! out
                 tvs_opts_scatt(sensorId),                            &! in
                 nlv_T,                                               &! in
                 chanprof,                                            &! in
                 frequencies,                                         &! in
                 tvs_profiles_nl(sensorTovsIndexes(1:profileCount)),  &! in
                 tvs_cld_profiles_nl(sensorTovsIndexes(1:profileCount)), &! in
                 tvs_coefs(sensorId),                                 &! in
                 tvs_coef_scatt(sensorId),                            &! in
                 calcemis,                                            &! in
                 emissivity_local,                                    &! inout
                 radiancedata_d) 
          else
            ! set the cloud profile in tvs_profiles_nl to zero
            call updateCloudInTovsProfile(sensorTovsIndexes(1:profileCount), &
                                          nlv_T,                      &
                                          mode='save',                &
                                          beSilent=.true.)
            call rttov_parallel_direct(                               &
                 rttov_err_stat,                                      & ! out
                 chanprof,                                            & ! in
                 tvs_opts(sensorId),                                  & ! in
                 tvs_profiles_nl(sensorTovsIndexes(1:profileCount)),  & ! in
                 tvs_coefs(sensorId),                                 & ! in
                 transmission,                                        & ! inout
                 radiancedata_d,                                      & ! inout
                 calcemis=calcemis,                                   & ! in
                 emissivity=emissivity_local,                         & ! inout
                 nthreads=nthreads      )   
          end if ! run rttovScatt

          ! save in obsSpaceData
          loopClearSky1: do btIndex = 1, btCount
            profileIndex = chanprof(btIndex)%prof
            channelIndex = chanprof(btIndex)%chan
            tovsIndex = sensorTovsIndexes(profileIndex)

            clearMwRadiance = radiancedata_d % bt(btIndex)

            headerIndex = tvs_headerIndex(tovsIndex)
            if ( headerIndex < 1 ) cycle loopClearSky1
            istart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
            iend = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + istart - 1

            ifBodyIndexFound = .false.
            loopClearSky2: do bodyIndex = istart, iend
              call tvs_getChannelNumIndexFromPPP( obsSpaceData, headerIndex, bodyIndex, &
                                                  channelNumber, channelIndexFound )
              if ( channelIndex == channelIndexFound ) then
                ifBodyIndexFound = .true.
                exit loopClearSky2
              end if
            end do loopClearSky2

            if ( .not. ifBodyIndexFound ) call utl_abort('tvs_rttov: bodyIndex not found.')

            if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated) then
              call obs_bodySet_r(obsSpaceData, OBS_BTCL, bodyIndex, clearMwRadiance)
            end if
          end do loopClearSky1

          ! restore the cloud profiles in ...
          if (runObsOperatorWithHydrometeors) then
            ! tvs_cld_profiles_nl
            call updateCloudInTovsCloudProfile(sensorTovsIndexes(1:profileCount), &
                                          nlv_T,                             &
                                          mode='restore',                    &
                                          beSilent=.true.)
          else
            ! tvs_profiles_nl
            call updateCloudInTovsProfile(sensorTovsIndexes(1:profileCount), &
                                          nlv_T,                             &
                                          mode='restore',                    &
                                          beSilent=.true.)
          end if
        end if ! run clear-sky RTTOV

        if (runObsOperatorWithHydrometeors) then
          if (.not. beSilent) write(*,*) 'before rttov_scatt...', sensorID, profileCount
          call rttov_scatt(                                         &
               rttov_err_stat,                                      &! out
               tvs_opts_scatt(sensorId),                            &! in
               nlv_T,                                               &! in
               chanprof,                                            &! in
               frequencies,                                         &! in
               tvs_profiles_nl(sensorTovsIndexes(1:profileCount)),  &! in
               tvs_cld_profiles_nl(sensorTovsIndexes(1:profileCount)), &! in
               tvs_coefs(sensorId),                                 &! in
               tvs_coef_scatt(sensorId),                            &! in
               calcemis,                                            &! in
               emissivity_local,                                    &! inout
               radiancedata_d) 
          if ( .not. beSilent ) write(*,*) 'after rttov_scatt...'
        else
          if (.not. beSilent) write(*,*) 'before rttov_parallel_direct...', sensorID, profileCount
          
          call rttov_parallel_direct(                               &
               rttov_err_stat,                                      & ! out
               chanprof,                                            & ! in
               tvs_opts(sensorId),                                  & ! in
               tvs_profiles_nl(sensorTovsIndexes(1:profileCount)),  & ! in
               tvs_coefs(sensorId),                                 & ! in
               transmission,                                        & ! inout
               radiancedata_d,                                      & ! inout
               calcemis=calcemis,                                   & ! in
               emissivity=emissivity_local,                         & ! inout
               nthreads=nthreads      )
          if ( .not. beSilent ) write(*,*) 'after rttov_parallel_direct...'
          
        end if

      end if ! if (bgckMode .and. tvs_isInstrumHyperSpectral(instrum))

      if ( .not. beSilent ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'     

      if (rttov_err_stat /= 0) then
        write(*,*) 'Error in rttov_parallel_direct',rttov_err_stat
        call utl_abort('tvs_rttov')
      end if
                                        
      !    2.4  Store hx in the structure tvs_radiance

      do btIndex = 1, btCount
        profileIndex = chanprof(btIndex)%prof
        channelIndex = chanprof(btIndex)%chan
        tovsIndex = sensorTovsIndexes(profileIndex)
        tvs_radiance(tovsIndex) % bt(channelIndex) = radiancedata_d % bt(btIndex)

        if ( bgckMode ) then
          if ( .not. associated(tvs_radiance(tovsIndex)  % clear)) then 
            allocStatus = 0
            allocate( tvs_radiance(tovsIndex)  % clear  ( tvs_nchan(sensorId)  ), stat= allocStatus(1) )
            !  allocate overcast black cloud sky radiance output
            allocate( tvs_radiance(tovsIndex)  % overcast  (nlv_T - 1, tvs_nchan(sensorId) ), stat=allocStatus(2))
            call utl_checkAllocationStatus(allocStatus(1:2), ' tvs_rttov')
          end if
          tvs_radiance(tovsIndex) % clear(channelIndex) =  &
               radiancedata_d %clear(btIndex)
          do levelIndex = 1, nlv_T - 1
            tvs_radiance(tovsIndex) % overcast(levelIndex,channelIndex) =   &
                 radiancedata_d % overcast(levelIndex,btIndex)
          end do
          if (.not. allocated(tvs_transmission)) call tvs_allocTransmission(nlv_T)
        end if

        if ( allocated( tvs_transmission) ) then
          do levelIndex = 1, nlv_T
            tvs_transmission(tovsIndex) % tau_levels(levelIndex,channelIndex) = &
                 transmission % tau_levels(levelIndex,btIndex)
          end do
          
          tvs_transmission(tovsIndex) % tau_total(channelIndex) = &
               transmission % tau_total(btIndex)
        end if

        if ( allocated(tvs_emissivity) ) then
          tvs_emissivity(channelIndex,tovsIndex) = emissivity_local(btIndex)%emis_out
        end if
          
      end do

      ! Append Surface Emissivity into ObsSpaceData
      do btIndex = 1, btCount
        bodyIndex = tvs_bodyIndexFromBtIndex(btIndex)
        if (bodyIndex > 0) call obs_bodySet_r(obsSpaceData, OBS_SEM, bodyIndex, emissivity_local(btIndex)%emis_out)
      end do

      !    Deallocate memory
      asw = 0 ! 0 to deallocate
      call rttov_alloc_direct(         &
           allocStatus(1),             &
           asw,                        &
           nprofiles=profileCount,     & ! (not used)
           nchanprof=btCount,          &
           nlevels=nlv_T,              &
           chanprof=chanprof,          &
           opts=tvs_opts(sensorId),    &
           coefs=tvs_coefs(sensorId),  &
           transmission=transmission,  &
           radiance=radiancedata_d,    &
           calcemis=calcemis,          &
           emissivity=emissivity_local,&
           init=.true.)

 
      if (useUofWIREmiss) then
        deallocate ( uOfWLandWSurfaceEmissivity  ,stat=allocStatus(2) )
      end if
      deallocate ( surfem1    ,stat=allocStatus(3) )
      if (allocated(frequencies)) deallocate (frequencies, stat=allocStatus(4))
      call utl_checkAllocationStatus(allocStatus, ' tvs_rttov', .false.)
      
    end do sensor_loop
    
    deallocate(sensorTovsIndexes)

  end subroutine tvs_rttov

  !--------------------------------------------------------------------------
  !  tvs_getMWemissivityFromAtlas
  !--------------------------------------------------------------------------
  subroutine tvs_getMWemissivityFromAtlas(originalEmissivity, updatedEmissivity, sensorId, chanprof, sensorTovsIndexes)
    implicit none

    real(8), intent(in)                  :: originalEmissivity(:)
    type (rttov_emissivity), intent(out) :: updatedEmissivity(:)
    integer, intent(in)                  :: sensorId
    type (rttov_chanprof), intent(in)    :: chanprof(:)
    integer, intent(in)                  :: sensorTovsIndexes(:)

    integer :: returnCode
    real(8) :: mWAtlasSurfaceEmissivity(size(originalEmissivity))
    integer :: btCount, profileCount
    integer :: profileIndex, btIndex, sensorIndex
    
    btCount = size( originalEmissivity )
    if (useMWEmissivityAtlas) then

      if (.not. allocated (tvs_atlas)) allocate(tvs_atlas(tvs_nsensors))
      if ( .not. tvs_atlas(sensorId)%init) then
        call rttov_setup_emis_atlas( returnCode, &! out
             tvs_opts(sensorId),                 &! in
             tvs_profiles_nl(1)%date(2),         &! in
             atlas_type_mw,                      &! in
             tvs_atlas(sensorId),                &! inout
             atlas_id = mWAtlasId,               &! in ! 1 TELSEM2, 2 CNRM
             coefs = tvs_coefs(sensorId)  )
        if (returnCode /= 0) then
          write(*,*) 'Error in rttov_atlas_setup MW',returnCode
          call utl_abort('tvs_getMWemissivityFromAtlas')
        end if
      end if
   
      call rttov_get_emis( returnCode,            & ! out
           tvs_opts(sensorId),                    & ! in
           chanprof,                              & ! in
           tvs_profiles_nl(sensorTovsIndexes(:)), & ! in
           tvs_coefs(sensorId),                   & ! in
           tvs_atlas(sensorId),                   & ! in
           mWAtlasSurfaceEmissivity)                ! out
    
      if (returnCode /= 0) then
        write(*,*) 'Error in rttov_get_emis MW', returnCode
        call utl_abort('tvs_getMWemissivityFromAtlas')
      end if

      profileCount = size( sensorTovsIndexes )

      do profileIndex=1, profileCount !loop on profiles
        sensorIndex = sensorTovsIndexes(profileIndex)
  
        do btIndex=1, btCount !loop on channels
          if (chanprof(btIndex)%prof==profileIndex) then
            ! Now we have 0.75 in originalEmissivity(:) for land and sea ice
            ! and the MW atlas emissivity in mWAtlasSurfaceEmissivity(:)

            if ( tvs_profiles_nl(sensorIndex)% skin % surftype == surftype_land .and. &
                 mWAtlasSurfaceEmissivity(btIndex) > 0.d0 .and. &
                 mWAtlasSurfaceEmissivity(btIndex) <= 1.d0 ) then ! check for missing values
              updatedEmissivity(btIndex)%emis_in = mWAtlasSurfaceEmissivity(btIndex)

            else
              updatedEmissivity(btIndex)%emis_in = originalEmissivity(btIndex)
            end if
            ! Note that emissivity above sea-ice is not modified
          end if
        end do
      end do
    else
      updatedEmissivity(:)%emis_in = originalEmissivity(:)
    end if
  end subroutine tvs_getMWemissivityFromAtlas

  !--------------------------------------------------------------------------
  !  comp_ir_emiss
  !--------------------------------------------------------------------------
  subroutine comp_ir_emiss (emiss, wind, angle, nchn, np, mchannel)
    !
    ! :Purpose: Computes water infrared emissivity for a specific set of
    !           channel indices, wind speed and zenith angle.
    !
  
    implicit none

    !Arguments
    real(8), intent(out) :: emiss(nchn,np) ! emissivities (0.-1.)
    real(8), intent(in)  :: wind(np) ! wind: surface wind speed (m/s)
    real(8), intent(in)  :: angle(np) ! angle: viewing angle (deg)
    integer, intent(in)  :: nchn ! number of channels to process
    integer, intent(in)  :: np     !number of locations
    integer, intent(in)  :: mchannel(nchn) ! vector of channel indices to process

    !Locals
    integer, parameter :: MaxWn = 19
    integer, parameter :: Nparm=3
    integer, parameter :: MaxChan=19

    real (8),parameter :: Theta(Nparm,MaxWn) = (/ &
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
         1503.969d0, 24.97931d0, 133.4392d0 /)

    real (8),parameter ::  C(Nparm,2,MaxWn) = (/                                 &  
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
         0.9192775350344998d0,-1.0369254272157462D-03, 2.8000594542037504D-05 /)

    real (8) a(MaxChan),b(MaxChan),cc(MaxChan)  ! local variable
    real (8) ww
    integer Index,Ichan,IP


    do Ichan = 1 , Nchn

      Index = Mchannel(Ichan)

      do IP=1,NP

        ww = wind(IP)
        a(Ichan) = c(1,1,Index) + c(2,1,Index) * ww    &  
             + c(3,1,Index) * ww * ww
        b(Ichan) = c(1,2,Index) + c(2,2,Index) * ww    &
             + c(3,2,Index)* ww * ww

        cc(Ichan) = Theta(1,Index) + Theta(2,Index) * ww

        emiss(Ichan,IP) = a(Ichan) + (b(Ichan) - a(Ichan)) *   & 
             exp(( (theta(3,Index) - 60.d0)**2.d0              &
             - (angle(IP) - theta(3,Index))**2.d0 ) / CC(Ichan))
       
      end do
      
    end do

  end subroutine comp_ir_emiss

  !--------------------------------------------------------------------------
  !  pcnt_box
  !--------------------------------------------------------------------------
  subroutine pcnt_box(f_low, f_high, nprf, ilat, ilon, klat, klon, ireduc)
    !
    ! :Purpose: Computes a low resolution feature form a high
    !           resolution one by averaging.
    !           example: use for percentage of water
    implicit none

    ! Arguments:
    real(8), intent(out)  :: f_low(nprf)       ! Low resolution field
    real(8), intent(in)   :: f_high(klon, klat)! High resolution field 
    integer, intent(in)   :: nprf              ! Number of profiles
    integer, intent(in)   :: ilat(nprf)        ! Y-coordinate of profile
    integer, intent(in)   :: ilon(nprf)        ! X-coordinate of profile
    integer, intent(in)   :: klon              ! Max value of latitude indices
    integer, intent(in)   :: klat              ! Max value of longitude indices
    integer, intent(in)   :: ireduc            ! Means a 2xireduc+1 by 2xireduc+1 averaging

    ! Locals
    integer :: nplon, jdlo1, jdlo2, jlon1, jlon2
    integer :: nx, ilat1, ilat2, ilon1, ilon2, jn, ii, jj
   
    profiles : do jn = 1,nprf

      nplon = 0

      ! normal limits

      ilat1=max(ilat(jn)-ireduc,1)
      ilat2=min(ilat(jn)+ireduc,klat)
      ilon1=max(ilon(jn)-ireduc,1)
      ilon2=min(ilon(jn)+ireduc,klon)

      if (ilon1 == 1 .or. ilon2 == klon) then
        ! border cases for longitudes
        jdlo1 = ilon(jn)-ireduc
        jdlo2 = ilon(jn)+ireduc

        if ( jdlo1 <= 0 ) then
          nplon = 1
          jlon1 = klon + jdlo1
          jlon2 = klon
        else if ( jdlo2 > klon ) then
          nplon = 1
          jlon1 = 1
          jlon2 = jdlo2 - klon
        end if
      end if

      nx = 0
      f_low(jn) = 0.d0
     
      do jj = ilat1, ilat2

        do ii = ilon1, ilon2
          nx = nx + 1
          f_low(jn) = f_low(jn) + f_high(ii,jj)         
        end do
        
        if (nplon == 1) then
          ! additional cases at border 1-klon
          do ii = jlon1, jlon2
            nx = nx + 1
            f_low(jn) = f_low(jn) + f_high(ii,jj)         
          end do
        end if

      end do
      
      f_low(jn) = f_low(jn) / dble(nx)

    end do profiles

  end subroutine pcnt_box

  !--------------------------------------------------------------------------
  !  emis_read_climatology
  !--------------------------------------------------------------------------
  subroutine emis_read_climatology
    !
    ! :Purpose: Read information about ceres surface type and water fraction.
    !
    ! :Arguments:
    !        :none:
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
      call utl_abort('Problem with file ceres_global.std in emis_read_climatology ')
    end if
   
  end subroutine emis_read_climatology

  !--------------------------------------------------------------------------
  !  emis_getIrEmissivity
  !--------------------------------------------------------------------------
  subroutine emis_getIrEmissivity (surfem1, nchn, sensorIndex, nprf, nchannels_max, sensorTovsIndexes)
    !
    ! :Purpose: Assign new ir surface emissivities based on
    !           cmc analysis surface albedo, sea ice fraction and snow mask
    !           in addition to ceres surface type and water fraction. 
    !           This is a subroutine that can apply to any instrument.
    !
    implicit none
   
    ! Arguments:
    real(8), intent(out) :: surfem1(nchannels_max) ! IR surface emissivity estimate (0-1)
    integer, intent(in)  :: nchn                   ! Number of channels
    integer, intent(in)  :: sensorindex            ! Sensor number
    integer, intent(in)  :: nprf                   ! Number of profiles
    integer, intent(in)  :: nchannels_max          ! Total number of observations treated
    integer, intent(in)  :: sensorTovsIndexes( nprf )         ! Profile position number

    !Locals:
    integer :: jc,jn
    integer :: ilat(nprf), ilon(nprf)
    real(8) :: latitudes(nprf), longitudes(nprf), satzang(nprf)
    real(8) :: wind_sfc(nprf), f_low(nprf), waven(nchn), em_oc(nchn,nprf), emi_mat(nchn,20)


    ! Information to extract (transvidage)
    ! latitudes(nprf) -- latitude (-90 to 90)
    ! longitudes(nprf) -- longitude (0 to 360)
    ! satzang(nprf) -- satellite zenith angle (deg)

    do jn = 1, nprf
      latitudes(jn)   = tvs_profiles_nl(sensorTovsIndexes(jn))% latitude
      longitudes(jn)  = tvs_profiles_nl(sensorTovsIndexes(jn))% longitude
      satzang(jn)     = tvs_profiles_nl(sensorTovsIndexes(jn))% zenangle
    end do

    !  Assign surface properties from grid to profiles
    call interp_sfc(ilat,ilon, nprf,latitudes,longitudes,sensorTovsIndexes)


    !  Find the sensor bands (central) wavenumbers
    do jc = 1, nchn      
      waven(jc) = tvs_coefs(sensorIndex) % coef % ff_cwn(jc)
    end do


    !  Get the CERES emissivity matrix for all sensor wavenumbers and surface types
    call ceres_ematrix(emi_mat, waven,nchn)


    ! Refine water emissivities

    do jn = 1, nprf
      !       find surface wind
      wind_sfc(jn) = min(sqrt(tvs_profiles_nl(sensorTovsIndexes(jn))%S2M%U**2 + tvs_profiles_nl(sensorTovsIndexes(jn))%S2M%V**2 + 1.d-12),15.d0)
    end do

    !     find new ocean emissivities     

    do jc = 1, nchn
      em_oc(jc,:)= emi_mat(jc,17)
    end do
    
    call emi_sea (em_oc, waven,satzang,wind_sfc,nprf,nchn)
    

    ! Get surface emissivities

    do jn = 1, nprf
      !       set albedo to 0.6 where snow is present
      if ( tvs_profiles_nl(sensorTovsIndexes(jn))%SKIN%SURFTYPE == 0 .and. tvs_surfaceParameters(sensorTovsIndexes(jn))%snow > 0.999 ) tvs_surfaceParameters(sensorTovsIndexes(jn))%albedo = 0.6
      !       if albedo too high no water
      if ( tvs_surfaceParameters(sensorTovsIndexes(jn))%albedo >= 0.55 ) tvs_surfaceParameters(sensorTovsIndexes(jn))%pcnt_wat = 0.
      !       if water and CMC ice present then sea ice
      if ( tvs_profiles_nl(sensorTovsIndexes(jn))%SKIN%SURFTYPE == 1 .and. tvs_surfaceParameters(sensorTovsIndexes(jn))%ice > 0.001 ) tvs_surfaceParameters(sensorTovsIndexes(jn))%ltype = 20
      !       if land and CMC snow present then snow
      if ( tvs_profiles_nl(sensorTovsIndexes(jn))%SKIN%SURFTYPE == 0 .and. tvs_surfaceParameters(sensorTovsIndexes(jn))%snow > 0.999 ) tvs_surfaceParameters(sensorTovsIndexes(jn))%ltype = 15
      do jc=1,nchn
        surfem1((jn-1)*nchn+jc) =  tvs_surfaceParameters(sensorTovsIndexes(jn))%pcnt_wat * em_oc(jc,jn)  +   &
             ( 1.d0 - tvs_surfaceParameters(sensorTovsIndexes(jn))%pcnt_wat ) * emi_mat(jc,tvs_surfaceParameters(sensorTovsIndexes(jn))%ltype)
      end do
    end do

    ! Find the regional water fraction (here in a 15x15 pixel box centered on profile)
    call pcnt_box (f_low, waterFraction,nprf,ilat,ilon,kslat,kslon,7)

    do jn = 1, nprf
      tvs_surfaceParameters(sensorTovsIndexes(jn))%pcnt_reg = f_low(jn)
    end do

  end subroutine emis_getIrEmissivity

  !--------------------------------------------------------------------------
  !  interp_sfc
  !--------------------------------------------------------------------------
  subroutine interp_sfc (ilat, ilon, nprf, latitudes, longitudes, sensorTovsIndexes)
    !
    ! :Purpose: Associate surface albedo, ice fraction, snow depth 
    !           and ceres surface type and water fraction to observations profiles.

    implicit none

    ! Arguments
    integer, intent(out) :: ilat(nprf)   ! y-coordinate of profile
    integer, intent(out) :: ilon(nprf)   ! x-coordinate of profile 
    integer, intent(in)  :: nprf         ! number of profiles
    real(8), intent(in)  :: latitudes(nprf)   ! latitude (-90s to 90n)
    real(8), intent(in)  :: longitudes(nprf)   ! longitude (0 to 360)
    integer, intent(in)  :: sensorTovsIndexes(nprf) ! observation index

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
    integer            :: jn
    character(len=1)   :: typvar
    character(len=1)   :: grtyp3,grtyp4,grtyp5
    character(len=2)   :: nomvar, snowvar
    character(len=8)   :: etiket
    integer, external  :: fnom,fstouv,fstinf,fstprm,fstfrm,fclos
    integer, external  :: ezqkdef,ezdefset
    real(8)            :: zig1,zig2,zig3,zig4
    integer            :: ig1obs,ig2obs,ig3obs,ig4obs
    real (8)           :: alat, alon, zzlat, zzlon
    ! fields on input grid
    real(8), allocatable :: glace(:,:), neige(:,:), alb(:,:)
    ! fields on output grid
    real(8)              :: glace_intrpl(nprf,1), neige_intrpl(nprf,1), alb_intrpl(nprf,1)


    ! printout header
    write(*,*) 
    write(*,*) 'SUBROUTINE interp_sfc'
    write(*,*) '---------------------'
    write(*,*) ' called multiple time by bunch of ',nprf,' profiles'
    write(*,*) ' <RETURN CODES> SHOULD NOT BE NEGATIVE'
    write(*,*) '---------------------------------------------------'


    ! --- FOR CERES VARIABLES -------------
    !  Get number of pixels per degree of lat or lon

    alat = dble(kslat)/180.d0
    alon = dble(kslon)/360.d0

    do jn=1, nprf

      ! get lat and lon within limits if necessary
      zzlat = min(latitudes(jn),89.999d0)
      zzlat = max(Zzlat,-89.999d0)
      
      zzlon = min(longitudes(jn),359.999d0)
      zzlon = max(zzlon,0.d0)

      !  Find in which surface field pixel is located the observation profile

      ! Note : CERES grid at 1/6 resolution 
      !         N-S : starts at N pole and excludes S pole
      !         W-E : starts at longitude 0 and excludes longitude 360

      ilat(jn) = max( nint((zzlat + 90.d0) * alat),1) 
      ilon(jn) = nint(zzlon * alon) + 1
      if (ilon(jn) > kslon) ilon(jn) = 1

      !  Assign surface caracteristics to observation profiles

      tvs_surfaceParameters(sensorTovsIndexes(jn)) % ltype    = surfaceType(ilon(jn),ilat(jn))
      tvs_surfaceParameters(sensorTovsIndexes(jn)) % pcnt_wat = waterFraction(ilon(jn),ilat(jn))

    end do

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
    if ( iy3  <  0 ) then
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
    allocate ( glace  (ni3,nj3) )
    allocate ( neige  (ni4,nj4) )
    allocate ( alb    (ni5,nj5) )


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

    do jn=1, nprf
      tvs_surfaceParameters(sensorTovsIndexes(jn))%ice      = glace_intrpl(jn,1)
      tvs_surfaceParameters(sensorTovsIndexes(jn))%snow     = neige_intrpl(jn,1)
      tvs_surfaceParameters(sensorTovsIndexes(jn))%albedo   = alb_intrpl(jn,1)
    end do

    deallocate(glace,neige,alb)

  end subroutine interp_sfc

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
    implicit none

    ! Arguments
    integer, intent(in) :: nchn              ! number of bands for which emissivity is needed
    real(8), intent(in) :: waven(nchn)       ! wavenumbers (cm-1)
    real(8), intent(out):: emi_mat(nchn, 20) ! emissivity (0.0-1.0)

    ! locals
    integer          :: i, nc, nt
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



    do nt = 1, 20
      do nc = 1, nchn
        if ( waven(nc) > mid(1) ) then
          emi_mat(nc,nt) = emi_tab(1,nt)
        else if ( waven(nc) < mid(nb) ) then
          emi_mat(nc,nt) = emi_tab(nb,nt)
        else
          do i = 1, nb - 1
            if ( waven(nc) <= mid(i) .and. waven(nc) >= mid(i + 1) ) then
              dum = ( waven(nc) - mid(i) ) / ( mid(i + 1) - mid(i) )
              emi_mat(nc,nt) = emi_tab(i,nt) + ( emi_tab(i + 1,nt) - emi_tab(i,nt) ) * dum
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
  subroutine emi_sea(em_oc, wnum, angle, wind, np, nc)
    !
    ! :Purpose: GET OCEAN SURFACE EMISSIVITY
    ! :Note:    IMEM(NC), set to zero initially, on next call IMEM will have the
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
    real(8), intent(out)   :: em_oc(nc,np) ! Ocean emissivities (0.-1.)
    real(8), intent(in)    :: wnum(nc)     ! Channel wavenumbers (cm-1)
    real(8), intent(in)    :: angle(np)    ! Viewing angle (deg)
    real(8), intent(in)    :: wind(np)     ! Surface wind speed (m/s)
    integer, intent(in)    :: np           ! Number of profiles
    integer, intent(in)    :: nc           ! Number of channels

    ! Locals
    integer     :: i, k, l
    integer     :: imem(nc) 
    integer     :: mchan(2)
    real(8)     :: dum
    real(8)     :: emi2(2,np)

    ! Masuda's 19 wavelengths converted to wavenumber
    real(8), parameter :: refw(19)=(/ 2857.1d0, 2777.7d0, 2702.7d0, 2631.6d0, 2564.1d0, &
         2500.0d0, 2439.0d0, 1250.0d0, 1190.5d0, 1136.3d0,                           &
         1087.0d0, 1041.7d0, 1000.0d0, 952.38d0, 909.09d0,                           &
         869.57d0, 833.33d0, 800.00d0, 769.23d0/)


    ! imem options

    imem(:) = 0

    do I = 1, nc

      !  out of range
      if ( wnum(I) < 645.d0 .or. wnum(I) > refw(1) ) then
        write(*,'(A,1x,e12.4)') ' fatal: wavenumber out of range in emi_sea', wnum(I)
        stop
      else if ( wnum(I) <= refw(19) .and. wnum(I) > 645.d0 ) then
        !  extrapolated from 769 cm-1 to 645 cm-1: NOT FROM REAL DATA
        !  nevertheless thought to be much better than unity
        !  this is a region of relatively rapid emissivity change
        !  worst estimates for 700-645 cm-1, but these channels do not
        !  see the surface (strong co2 absorption).
        imem(I) = 18
      else
        !  CAUTION interpolation on large interval 1250-2439 cm-1
        !  where no data is available except that of ASTER. ASTER
        !  shows a relatively smooth variation with wavelength except
        !  for a sharp drop at 1600 cm-1 with highs at 1550 and 1650 cm-1
        !  with peak-to-peak variation of 1.5% in that narrow range.
        !  Worst estimates would be between 1400-1800 cm-1 in HIRS ch 12
        !  which only in very cold atmospheres sees the surface.
        do k = 1, 18
          if ( wnum(I) > refw(k + 1) .and. wnum(I) <= refw(k) ) then
            imem(I) = k
          end if
        end do

      end if
   
      mchan(1)= imem(I)
      mchan(2)= imem(I) + 1

      dum = ( wnum(I) - refw(mchan(1)) ) / ( refw(mchan(2)) - refw(mchan(1)) )

      call COMP_IR_EMISS(emi2, wind,angle,2,np,mchan)

      ! interpolation/extrapolation in wavenumber 

      do L = 1, np
  
        em_oc(I,L) = emi2(1,L) + ( emi2(2,L) - emi2(1,L) ) * dum
          
      end do

    end do


  end subroutine emi_sea


  !--------------------------------------------------------------------------
  !  tvs_getCommonChannelSet
  !--------------------------------------------------------------------------
  subroutine tvs_getCommonChannelSet(channels,countUniqueChannel, listAll)
    !
    !:Purpose: get common channels among all MPI tasks
    !
    implicit none
    !Arguments:
    integer, intent(in) :: channels(:)
    integer, intent(out):: countUniqueChannel, listAll(:)
    !Locals:
    integer :: channelsb(tvs_maxChannelNumber)
    integer :: ierr, i, j
    integer, allocatable :: listGlobal(:)
    logical :: found
     
    if (size(channels) > tvs_maxChannelNumber) then
      write(*,*) 'You need to increase tvs_maxChannelNumber in tovs_nl_mod !',size(channels), tvs_maxChannelNumber
      call utl_abort('tvs_getCommonChannelSet')
    end if

    if (mmpi_myid ==0) then
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
    if ( mmpi_myid == 0 ) then
      call isort(listGlobal, mmpi_nprocs*tvs_maxChannelNumber)
      do i=1, mmpi_nprocs * tvs_maxChannelNumber
        if (listGlobal(i) > 0) then
          found = .false.
          LOOPJ: do j=countUniqueChannel,1,-1
            if (listGlobal(i) == listAll(j) ) then
              found =.true.
              exit LOOPJ
            end if
          end do LOOPJ
          if (.not.found) then
            countUniqueChannel = countUniqueChannel + 1
            listAll(countUniqueChannel) = listGlobal(i)
          end if
        end if
      end do
    end if
    
    call rpn_comm_bcast(countUniqueChannel, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(listAll(1:countUniqueChannel), countUniqueChannel, 'MPI_INTEGER', 0, 'GRID', ierr)

    deallocate(listGlobal)

  end subroutine tvs_getCommonChannelSet


  !--------------------------------------------------------------------------
  !  tvs_rttov_read_coefs
  !--------------------------------------------------------------------------
  subroutine tvs_rttov_read_coefs(err, coefs, opts, channels, instrument)
    !
    ! :Purpose: MPI wrapper for rttov_read_coefs
    !           the coefficient files are read by MPI task 0
    !           and then broadcasted to the other tasks according to the selected
    !           channels. Argument channels is mandatory (it is optional in rttov_setup)
    !           optional argument channels_rec was removed (it is useful only in principal component mode)
    !           other optionnal arguments were removed :
    !                          *  form_coef, to specify format
    !                          *  form_scaer,        
    !                          *  form_sccld,        
    !                          *  form_pccoef,       
    !                          *  file_coef, to specify filename
    !                          *  file_scaer,        
    !                          *  file_sccld,        
    !                          *  file_pccoef,       
    !                          *  file_id_coef, to specify fortran unit number
    !                          *  file_id_scaer,     
    !                          *  file_id_sccld,     
    !                          *  file_id_pccoef,
    !                          *  path,  to specify the path to look for coefficient files
    !
    !           if necessary these arguments could be  added (ask S. Heilliette)
    !           also this subroutine will work only for clear sky radiance computations
    !           if somebody wants to do realistic cloud or aerosol affected radiance simulations
    !           some changes are needed. Ask me in that case. (S. Heilliette) 
    !           It is implicitely assumed that the options are the same for all MPI tasks for a given instrument
    !           No check will be done (options for task 0 will be used for all tasks). 
    !           Only differences in channel lists are accounted for.
    !

    implicit none

    ! Arguments:
    integer(kind=jpim), intent(out) :: err          ! Error status
    type(rttov_coefs),  intent(out) :: coefs        ! Rttov coefficient structure
    type(rttov_options), intent(in) :: opts         ! Rttov option structure
    integer(kind=jpim), intent(in)  :: channels(:)  ! Channel list
    integer(kind=jpim), intent(in)  :: instrument(3)! Instrument vector

    ! Locals:
    real(8), allocatable :: bigArray(:,:,:,:)
    integer :: i, j, ichan, ierr, countUniqueChannel, indexchan(size(channels)), listAll(tvs_maxChannelNumber)
    logical :: associated0
    integer :: nlte_count, nlte_start,isol,isat,nlte_file_nchan
    integer, allocatable :: nlte_chans(:) 

    write(*,*) 'tvs_rttov_read_coefs: Starting'

    ! First step: we should determine a common set of channels among MPI tasks
    call tvs_getCommonChannelSet(channels,countUniqueChannel, listAll)

    ! Second step: mpi task 0 will do the job
    if ( mmpi_myid == 0 ) then
      call rttov_read_coefs ( &
           err,             &! out
           coefs,           &
           opts,            &
           instrument=instrument,      &! in
           channels=listAll(1:countUniqueChannel)  )     ! in option
    else
      call rttov_nullify_coef(coefs%coef)
    end if

    ! Third step: common (i.e. independent from the channel list) parameters are simply broadcasted to other processors
    ! Scalar and fixed size arrays and  strings first
    coefs%initialised = .true.                                                    ! Logical flag for initialization
    call rpn_comm_bcast(coefs%coef%id_platform, 1, 'MPI_INTEGER', 0, 'GRID', ierr)! OK
    call rpn_comm_bcast(coefs%coef%id_sat, 1, 'MPI_INTEGER', 0, 'GRID', ierr)     ! OK
    call rpn_comm_bcast(coefs%coef%id_inst, 1, 'MPI_INTEGER', 0, 'GRID', ierr)    ! OK
    call rpn_comm_bcast(coefs%coef%id_sensor, 1, 'MPI_INTEGER', 0, 'GRID', ierr)  ! OK 
    call rpn_comm_bcast(coefs%coef%id_comp_lvl, 1, 'MPI_INTEGER', 0, 'GRID', ierr)! OK
    call rpn_comm_bcast(coefs%coef%id_comp_pc, 1, 'MPI_INTEGER', 0, 'GRID', ierr) ! OK

    call rpn_comm_bcast(coefs%coef%fmv_model_ver, 1, 'MPI_INTEGER', 0, 'GRID', ierr) !ok
    call rpn_comm_bcast(coefs%coef%fmv_chn, 1, 'MPI_INTEGER', 0, 'GRID', ierr) !ok
    call rpn_comm_bcast(coefs%coef%fmv_gas, 1, 'MPI_INTEGER', 0, 'GRID', ierr) !ok
    call rpn_comm_bcast(coefs%coef%fmv_ori_nchn, 1, 'MPI_INTEGER', 0, 'GRID', ierr) !ok

    call rpn_comm_bcast(coefs%coef%nmixed, 1, 'MPI_INTEGER', 0, 'GRID', ierr)          ! number of variables/predictors for Mixed Gases
    call rpn_comm_bcast(coefs%coef%nwater, 1, 'MPI_INTEGER', 0, 'GRID', ierr)          ! number of variables/predictors for Water Vapour
    call rpn_comm_bcast(coefs%coef%nozone, 1, 'MPI_INTEGER', 0, 'GRID', ierr)          ! number of variables/predictors for Ozone
    call rpn_comm_bcast(coefs%coef%nwvcont, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of variables/predictors for WV continuum
    call rpn_comm_bcast(coefs%coef%nco2, 1, 'MPI_INTEGER', 0, 'GRID', ierr)            ! number of variables/predictors for CO2
    call rpn_comm_bcast(coefs%coef%nn2o , 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of variables/predictors for N2O
    call rpn_comm_bcast(coefs%coef%nco, 1, 'MPI_INTEGER', 0, 'GRID', ierr)             ! number of variables/predictors for CO
    call rpn_comm_bcast(coefs%coef%nch4, 1, 'MPI_INTEGER', 0, 'GRID', ierr)            ! number of variables/predictors for CH4
    call rpn_comm_bcast(coefs%coef%nso2, 1, 'MPI_INTEGER', 0, 'GRID', ierr)            ! number of variables/predictors for SO2

    call rpn_comm_bcast(coefs%coef%nlevels, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of levels(pres/absorber) same for all gases
    call rpn_comm_bcast(coefs%coef%nlayers, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of layers(pres/absorber) nlevels-1
    call rpn_comm_bcast(coefs%coef%pmc_nlay, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%pmc_nvar, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%IncZeeman, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)       ! Flag to include Zeeman effect for this sensor
    call rpn_comm_bcast(coefs%coef%solarcoef, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)       ! Flag to include solar reflection
    call rpn_comm_bcast(coefs%coef%nltecoef, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)        ! Flag to include nlte corrections
    call rpn_comm_bcast(coefs%coef%pmc_shift, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)

    call rpn_comm_bcast(coefs%coef%ncmixed, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of coefficients for Mixed Gases
    call rpn_comm_bcast(coefs%coef%ncwater, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of coefficients for Water Vapour
    call rpn_comm_bcast(coefs%coef%ncozone, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of coefficients for Ozone
    call rpn_comm_bcast(coefs%coef%ncwvcont, 1, 'MPI_INTEGER', 0, 'GRID', ierr)        ! number of coefficients for WV continuum
    call rpn_comm_bcast(coefs%coef%ncco2, 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for CO2
    call rpn_comm_bcast(coefs%coef%ncn2o, 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for N2O
    call rpn_comm_bcast(coefs%coef%ncco , 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for CO
    call rpn_comm_bcast(coefs%coef%ncch4, 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for CH4
    call rpn_comm_bcast(coefs%coef%ncso2, 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for SO2

    call rpn_comm_bcast(coefs%coef%nccmixed, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of coefficients for Mixed Gases
    call rpn_comm_bcast(coefs%coef%nccwater, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of coefficients for Water Vapour
    call rpn_comm_bcast(coefs%coef%nccozone, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of coefficients for Ozone
    call rpn_comm_bcast(coefs%coef%nccwvcont, 1, 'MPI_INTEGER', 0, 'GRID', ierr)        ! number of coefficients for WV continuum
    call rpn_comm_bcast(coefs%coef%nccco2, 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for CO2
    call rpn_comm_bcast(coefs%coef%nccn2o, 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for N2O
    call rpn_comm_bcast(coefs%coef%nccco , 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for CO
    call rpn_comm_bcast(coefs%coef%nccch4, 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for CH4
    call rpn_comm_bcast(coefs%coef%nccso2, 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for SO2

    call rpn_comm_bcast(coefs%coef%ws_nomega, 1, 'MPI_INTEGER', 0, 'GRID', ierr)

    call rpn_comm_bcast(coefs%coef%id_creation_date, 3, 'MPI_INTEGER', 0, 'GRID', ierr) ! OK 
    call rpn_comm_bcastc(coefs%coef%id_creation, 80, 'MPI_CHARACTER', 0, 'GRID', ierr)  ! OK 
    call rpn_comm_bcastc(coefs%coef%id_Common_name, 32, 'MPI_CHARACTER', 0, 'GRID', ierr) !OK
    do i=1, 100
      call rpn_comm_bcastc(coefs%coef%line_by_line(i), 132, 'MPI_CHARACTER', 0, 'GRID', ierr) !ok
      call rpn_comm_bcastc(coefs%coef%readme_srf(i), 132, 'MPI_CHARACTER', 0, 'GRID', ierr) !ok
    end do
    call rpn_comm_bcastc(coefs%coef%fmv_model_def, 32, 'MPI_CHARACTER', 0, 'GRID', ierr)  !OK
    call rpn_comm_bcast(coefs%coef%fc_planck_c1, 1, 'MPI_REAL8', 0, 'GRID', ierr) ! first radiation constant (mW/(m2*sr*cm-4)) !ok
    call rpn_comm_bcast(coefs%coef%fc_planck_c2, 1, 'MPI_REAL8', 0, 'GRID', ierr) !second radiation constant (mW/(m2*sr*cm-4)) !ok
    call rpn_comm_bcast(coefs%coef%fc_sat_height, 1, 'MPI_REAL8', 0, 'GRID', ierr)! satellite nominal altitude (km) !ok

    call rpn_comm_bcast(coefs%coef%pmc_lengthcell, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%pmc_tempcell, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%pmc_betaplus1, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    ! FASTEM section
    call rpn_comm_bcast(coefs%coef%ssirem_ver, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%iremis_version, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%iremis_ncoef, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%iremis_angle0, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%iremis_tskin0, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%ratoe, 1, 'MPI_REAL8', 0, 'GRID', ierr) 
    ! then variable size vectors
    ! this one must be done first because it is used to dimension other ones ....
    if (mmpi_myid > 0) allocate( coefs%coef%fmv_lvl(coefs%coef%fmv_gas))
    call rpn_comm_bcast(coefs%coef%fmv_lvl, coefs%coef%fmv_gas, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (mmpi_myid > 0) then
      if (coefs%coef%nltecoef) allocate( coefs%coef%nlte_coef)
      allocate( coefs%coef%fmv_gas_id(coefs%coef%fmv_gas))
      allocate( coefs%coef%fmv_gas_pos(ngases_max)) !size is from rttov consts
      allocate( coefs%coef%fmv_var(coefs%coef%fmv_gas))
      allocate( coefs%coef%fmv_coe(coefs%coef%fmv_gas)) ! number of coefficients by gas (fmv_gas)
      allocate( coefs%coef%fmv_ncorr(coefs%coef%fmv_gas)) ! number of coefs by gas for correction term (fmv_gas) (v13 only)
      allocate( coefs%coef%ws_npoint(coefs%coef%ws_nomega) )
      allocate( coefs%coef%ws_k_omega(coefs%coef%ws_nomega) )
      allocate( coefs%coef%ref_prfl_p(coefs%coef%nlevels ) )
      allocate( coefs%coef%ref_prfl_t(coefs%coef%nlevels, coefs%coef%fmv_gas ) )
      allocate( coefs%coef%ref_prfl_mr(coefs%coef%nlevels, coefs%coef%fmv_gas ) )
      allocate( coefs%coef%bkg_prfl_mr(coefs%coef%nlevels, coefs%coef%fmv_gas ) )
      allocate( coefs%coef%lim_prfl_p(coefs%coef%nlevels ) )
      allocate( coefs%coef%lim_prfl_tmax(coefs%coef%fmv_lvl(gas_id_mixed) ) )
      allocate( coefs%coef%lim_prfl_tmin(coefs%coef%fmv_lvl(gas_id_mixed) ) )
      allocate( coefs%coef%lim_prfl_gmax(coefs%coef%fmv_lvl(gas_id_mixed), coefs%coef%fmv_gas ) )
      allocate( coefs%coef%lim_prfl_gmin(coefs%coef%fmv_lvl(gas_id_mixed), coefs%coef%fmv_gas ) )
      allocate( coefs%coef%env_prfl_tmax(coefs%coef%nlevels ) )
      allocate( coefs%coef%env_prfl_tmin(coefs%coef%nlevels ) )
      allocate( coefs%coef%env_prfl_gmax(coefs%coef%nlevels, coefs%coef%fmv_gas ) )
      allocate( coefs%coef%env_prfl_gmin(coefs%coef%nlevels, coefs%coef%fmv_gas ) )
      allocate( coefs%coef%dpp(0:coefs%coef%nlayers) )
      allocate( coefs%coef%dp(coefs%coef%nlayers) )
      allocate( coefs%coef%tstar(coefs%coef%nlayers) )
      allocate( coefs%coef%tstar_r(coefs%coef%nlayers) )
      allocate( coefs%coef%tstar_wsum_r(0:coefs%coef%nlayers) )
      if (coefs%coef%fmv_model_ver == 8) &
           allocate( coefs%coef%tstarmod_wsum_r(coefs%coef%nlayers) )
      if (coefs%coef%fmv_model_ver <= 9) &
           allocate( coefs%coef%tstar_uwsum_r(0:coefs%coef%nlayers) )
      allocate( coefs%coef%wstar(coefs%coef%nlayers) )
      allocate( coefs%coef%wstar_r(coefs%coef%nlayers) )
      allocate( coefs%coef%wstar_wsum_r(0:coefs%coef%nlayers) )
      allocate( coefs%coef%wtstar_wsum_r(0:coefs%coef%nlayers) )
    end if

    call broadcastI41dArray( coefs%coef%fmv_gas_id )
    call broadcastI41dArray( coefs%coef%fmv_gas_pos )
    call broadcastI41dArray( coefs%coef%fmv_var )
    call broadcastI41dArray( coefs%coef%fmv_coe )
    call broadcastI41dArray( coefs%coef%fmv_ncorr )
    call broadcastR81dArray( coefs%coef%ws_npoint )
    call broadcastR81dArray( coefs%coef%ws_k_omega )
    call broadcastR81dArray( coefs%coef%ref_prfl_p )
    call broadcastR82dArray( coefs%coef%ref_prfl_t )
    call broadcastR82dArray( coefs%coef%ref_prfl_mr )
    call broadcastR82dArray( coefs%coef%bkg_prfl_mr )
    call broadcastR81dArray( coefs%coef%lim_prfl_p )
    call broadcastR81dArray( coefs%coef%lim_prfl_tmax )
    call broadcastR81dArray( coefs%coef%lim_prfl_tmin )
    call broadcastR82dArray( coefs%coef%lim_prfl_gmax )
    call broadcastR82dArray( coefs%coef%lim_prfl_gmin )
    call broadcastR81dArray( coefs%coef%env_prfl_tmax )
    call broadcastR81dArray( coefs%coef%env_prfl_tmin )
    call broadcastR82dArray( coefs%coef%env_prfl_gmax )
    call broadcastR82dArray( coefs%coef%env_prfl_gmin )
    call broadcastR81dArray( coefs%coef%dp )
    call broadcastR81dArray( coefs%coef%dpp )
    call broadcastR81dArray( coefs%coef%tstar ) 
    call broadcastR81dArray( coefs%coef%tstar_r )
    call broadcastR81dArray( coefs%coef%tstar_wsum_r )
    if (coefs%coef%fmv_model_ver == 8) &
         call broadcastR81dArray( coefs%coef%tstarmod_wsum_r )
    if (coefs%coef%fmv_model_ver <= 9) &
         call broadcastR81dArray( coefs%coef%tstar_uwsum_r )
    call broadcastR81dArray( coefs%coef%wstar )
    call broadcastR81dArray( coefs%coef%wstar_r )
    call broadcastR81dArray( coefs%coef%wstar_wsum_r )
    call broadcastR81dArray( coefs%coef%wtstar_wsum_r )
    if (coefs%coef%nozone > 0) then
      if (mmpi_myid > 0) then
        allocate( coefs%coef%to3star(coefs%coef%nlayers) )
        allocate( coefs%coef%ostar(coefs%coef%nlayers) )
        allocate( coefs%coef%to3star_r(coefs%coef%nlayers) )
        allocate( coefs%coef%ostar_r(coefs%coef%nlayers) )
        allocate( coefs%coef%ostar_wsum_r(0:coefs%coef%nlayers) )
      end if
      call rpn_comm_bcast(coefs%coef%to3star, size(coefs%coef%to3star) , 'MPI_REAL8', 0, 'GRID', ierr) 
      call rpn_comm_bcast(coefs%coef%ostar, size(coefs%coef%ostar) , 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%to3star_r, size(coefs%coef%to3star_r) , 'MPI_REAL8', 0, 'GRID', ierr) 
      call rpn_comm_bcast(coefs%coef%ostar_r, size(coefs%coef%ostar_r) , 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%ostar_wsum_r, size(coefs%coef%ostar_wsum_r) , 'MPI_REAL8', 0, 'GRID', ierr) 
    end if
    if ( coefs%coef%nco2 > 0) then
      if (mmpi_myid>0) then
        allocate( coefs%coef%co2star(coefs%coef%nlayers) )
        allocate( coefs%coef%co2star_r(coefs%coef%nlayers) )
        allocate( coefs%coef%co2star_wsum_r(0:coefs%coef%nlayers) )
      end if
      call rpn_comm_bcast(coefs%coef%co2star, size(coefs%coef%co2star) , 'MPI_REAL8', 0, 'GRID', ierr) 
      call rpn_comm_bcast(coefs%coef%co2star_r, size(coefs%coef%co2star_r) , 'MPI_REAL8', 0, 'GRID', ierr) 
      call rpn_comm_bcast(coefs%coef%co2star_wsum_r, size(coefs%coef%co2star_wsum_r) , 'MPI_REAL8', 0, 'GRID', ierr) 
    end if
    if ( coefs%coef%nn2o > 0) then
      if (mmpi_myid>0) then
        allocate( coefs%coef%n2ostar(coefs%coef%nlayers) )
        allocate( coefs%coef%n2ostar_r(coefs%coef%nlayers) )
        allocate( coefs%coef%n2ostar_wsum_r(0:coefs%coef%nlayers) )
        allocate( coefs%coef%n2otstar_wsum_r(0:coefs%coef%nlayers) )
      end if
      call rpn_comm_bcast(coefs%coef%n2ostar, size(coefs%coef%n2ostar) , 'MPI_REAL8', 0, 'GRID', ierr) 
      call rpn_comm_bcast(coefs%coef%n2ostar_r, size(coefs%coef%n2ostar_r) , 'MPI_REAL8', 0, 'GRID', ierr) 
      call rpn_comm_bcast(coefs%coef%n2ostar_wsum_r, size(coefs%coef%n2ostar_wsum_r) , 'MPI_REAL8', 0, 'GRID', ierr) 
      call rpn_comm_bcast(coefs%coef%n2otstar_wsum_r, size(coefs%coef%n2otstar_wsum_r) , 'MPI_REAL8', 0, 'GRID', ierr) 
    end if
    if ( coefs%coef%nco > 0) then
      if (mmpi_myid>0) then
        allocate( coefs%coef%costar(coefs%coef%nlayers) )
        allocate( coefs%coef%costar_r(coefs%coef%nlayers) )
        allocate( coefs%coef%costar_wsum_r(0:coefs%coef%nlayers) )
        allocate( coefs%coef%cotstar_wsum_r(0:coefs%coef%nlayers) )
      end if
      call rpn_comm_bcast(coefs%coef%costar, size(coefs%coef%costar) , 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%costar_r, size(coefs%coef%costar_r) , 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%costar_wsum_r, size(coefs%coef%costar_wsum_r) , 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%cotstar_wsum_r, size(coefs%coef%cotstar_wsum_r) , 'MPI_REAL8', 0, 'GRID', ierr) 
    end if
    if ( coefs%coef%nch4 > 0) then
      if (mmpi_myid>0) then
        allocate( coefs%coef%ch4star(coefs%coef%nlayers) )
        allocate( coefs%coef%ch4star_r(coefs%coef%nlayers) )
        allocate( coefs%coef%ch4star_wsum_r(0:coefs%coef%nlayers) )
        allocate( coefs%coef%ch4tstar_wsum_r(coefs%coef%nlayers) )
      end if
      call rpn_comm_bcast(coefs%coef%ch4star, size(coefs%coef%ch4star) , 'MPI_REAL8', 0, 'GRID', ierr) 
      call rpn_comm_bcast(coefs%coef%ch4star_r, size(coefs%coef%ch4star_r) , 'MPI_REAL8', 0, 'GRID', ierr) 
      call rpn_comm_bcast(coefs%coef%ch4star_wsum_r, size(coefs%coef%ch4star_wsum_r) , 'MPI_REAL8', 0, 'GRID', ierr) 
      call rpn_comm_bcast(coefs%coef%ch4tstar_wsum_r, size(coefs%coef%ch4tstar_wsum_r) , 'MPI_REAL8', 0, 'GRID', ierr) 
    end if
    if (coefs%coef%nso2 > 0) then
      if (mmpi_myid>0) then
        allocate( coefs%coef%so2star(coefs%coef%nlayers) )
        allocate( coefs%coef%so2star_r(coefs%coef%nlayers) )
        allocate( coefs%coef%so2star_wsum_r(0:coefs%coef%nlayers) )
        allocate( coefs%coef%so2tstar_wsum_r(coefs%coef%nlayers) )
      end if
      call rpn_comm_bcast(coefs%coef%so2star, size(coefs%coef%so2star) , 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%so2star_r, size(coefs%coef%so2star_r) , 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%so2star_wsum_r, size(coefs%coef%so2star_wsum_r) , 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%so2tstar_wsum_r, size(coefs%coef%so2tstar_wsum_r) , 'MPI_REAL8', 0, 'GRID', ierr)
    end if
    ! Fourth step: channel dependent parameters are extracted according to the channel list and sent to each MPI task
    coefs%coef%fmv_chn = size( channels )

    if (mmpi_myid == 0) deallocate ( coefs%coef%ff_ori_chn)
    allocate ( coefs%coef%ff_ori_chn(coefs%coef%fmv_chn) )
    coefs%coef%ff_ori_chn =  channels

    do i=1,  coefs%coef%fmv_chn
      loopj2:do j=1, countUniqueChannel
        if (listAll(j) == channels(i)) then
          indexchan(i) = j
          exit loopj2
        end if
      end do loopj2
    end do
    
    ! 1D arrays first
    call extractI41dArray(coefs%coef%pw_val_chn, countUniqueChannel,indexchan)
    call extractI41dArray(coefs%coef%ff_val_chn, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%ff_cwn, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%ff_bco, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%ff_bcs, countUniqueChannel,indexchan) 
    call extractR81dArray(coefs%coef%ff_gam, countUniqueChannel,indexchan)
    call extractI41dArray(coefs%coef%tt_val_chn, countUniqueChannel,indexchan) 
    call extractR81dArray(coefs%coef%tt_a0, countUniqueChannel,indexchan) 
    call extractR81dArray(coefs%coef%tt_a1, countUniqueChannel,indexchan) 
    call extractI41dArray(coefs%coef%ss_val_chn, countUniqueChannel,indexchan) 
    call extractR81dArray(coefs%coef%ss_solar_spectrum, countUniqueChannel,indexchan)
    if ( coefs%coef%fmv_model_ver > 9) then
      call extractR81dArray(coefs%coef%ss_rayleigh_ext, countUniqueChannel,indexchan)
    end if
    call extractCmplx81dArray(coefs%coef%woc_waopc_ow, countUniqueChannel, indexchan) 
    call extractCmplx81dArray(coefs%coef%woc_waopc_fw, countUniqueChannel, indexchan)
  
    call extractI41dArray(coefs%coef%fastem_polar, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%pol_phi, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%ssirem_a0, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%ssirem_a1, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%ssirem_a2, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%ssirem_xzn1, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%ssirem_xzn2, countUniqueChannel,indexchan)

    call extractR81dArray(coefs%coef%planck1, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%planck2, countUniqueChannel,indexchan)
    if (coefs%coef%id_sensor == sensor_id_mw .or. coefs%coef%id_sensor == sensor_id_po) &
         call extractR81dArray(coefs%coef%frequency_ghz, countUniqueChannel,indexchan)

    call extractR81dArray(coefs%coef%pmc_pnominal, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%pmc_ppmc, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%pol_fac_v, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%pol_fac_h, countUniqueChannel,indexchan)

    ! 2D arrays
    call extractR82dArray(coefs%coef%iremis_coef,coefs%coef%iremis_ncoef,countUniqueChannel,indexchan)
    ! 3D arrays
    call extractR83dArray(coefs%coef%pmc_coef,coefs%coef%pmc_nlay,countUniqueChannel,coefs%coef%pmc_nvar,indexchan)

    ! then coefficients. It is more complicated with RTTOV12
    call dispatch_fast_coef(err, coefs%coef%thermal, coefs%coef%fmv_gas_id, coefs%coef%fmv_coe, coefs%coef%fmv_model_ver, &
         coefs%coef%nlayers, coefs%coef%fmv_gas)
    if (coefs%coef%fmv_model_ver > 9) THEN
      call dispatch_fast_coef(err, coefs%coef%thermal_corr, coefs%coef%fmv_gas_id, coefs%coef%fmv_ncorr, coefs%coef%fmv_model_ver, &
           coefs%coef%nlayers, coefs%coef%fmv_gas)
    end if
    if (coefs%coef%solarcoef) then
      call dispatch_fast_coef(err, coefs%coef%solar, coefs%coef%fmv_gas_id, coefs%coef%fmv_coe, coefs%coef%fmv_model_ver, &
           coefs%coef%nlayers, coefs%coef%fmv_gas)
      if (coefs%coef%fmv_model_ver > 9) THEN
        call dispatch_fast_coef(err, coefs%coef%solar_corr, coefs%coef%fmv_gas_id, coefs%coef%fmv_ncorr, coefs%coef%fmv_model_ver, &
             coefs%coef%nlayers, coefs%coef%fmv_gas)
      end if
    end if

    if (coefs%coef%nltecoef) then

      call rpn_comm_bcast(coefs%coef%nlte_coef%ncoef, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%nlte_coef%nsol, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%nlte_coef%nsat, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%nlte_coef%nchan, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%nlte_coef%start_chan, 1, 'MPI_INTEGER', 0, 'GRID', ierr)

      allocate(nlte_chans(coefs%coef%fmv_chn)) ! Index of selected channels in nlte_coefs array in the file

      nlte_count = 0  ! Number of NLTE channels being read in
      nlte_start = 0  ! Index (in input channel list) of first NLTE channel being read in
      do i = 1, coefs%coef%fmv_chn
        if (channels(i) >= coefs%coef%nlte_coef%start_chan .and. &
             channels(i) < coefs%coef%nlte_coef%start_chan + coefs%coef%nlte_coef%nchan) then
          nlte_count = nlte_count + 1
          nlte_chans(nlte_count) = channels(i) - coefs%coef%nlte_coef%start_chan + 1
          if (nlte_count == 1) nlte_start = i
        end if
      end do

      coefs%coef%nltecoef = ( nlte_count > 0)

      nlte_file_nchan  = coefs%coef%nlte_coef%nchan

      ! Reset NLTE channel variables according to input channel limit
      coefs%coef%nlte_coef%start_chan  = nlte_start
      coefs%coef%nlte_coef%nchan      = nlte_count

      if (mmpi_myid > 0) then
         allocate (coefs%coef%nlte_coef%sec_sat(coefs%coef%nlte_coef%nsat) )
         allocate (coefs%coef%nlte_coef%sol_zen_angle(coefs%coef%nlte_coef%nsol) )
         allocate (coefs%coef%nlte_coef%sat_zen_angle(coefs%coef%nlte_coef%nsat) )
         allocate (coefs%coef%nlte_coef%cos_sol(coefs%coef%nlte_coef%nsol) )
      end if
      call rpn_comm_bcast(coefs%coef%nlte_coef%sec_sat, coefs%coef%nlte_coef%nsat, 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%nlte_coef%sol_zen_angle, coefs%coef%nlte_coef%nsol, 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%nlte_coef%sat_zen_angle, coefs%coef%nlte_coef%nsat, 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%nlte_coef%cos_sol, coefs%coef%nlte_coef%nsol, 'MPI_REAL8', 0, 'GRID', ierr)

      allocate(bigArray(coefs%coef%nlte_coef%ncoef, coefs%coef%nlte_coef%nsat, &
           coefs%coef%nlte_coef%nsol, nlte_file_nchan) )

      if (mmpi_myid == 0) then
         do ichan = 1, nlte_file_nchan
            do isol = 1, coefs%coef%nlte_coef%nsol
               do isat = 1, coefs%coef%nlte_coef%nsat
                  do I=1, coefs%coef%nlte_coef%ncoef
                     bigArray(i, isat, isol, ichan) = coefs%coef%nlte_coef%coef(i,isat,isol,ichan)
                  end do
               end do
            end do
         end do
         deallocate(  coefs%coef%nlte_coef%coef )
      end if

      call rpn_comm_bcast(bigArray, size(bigArray), 'MPI_REAL8', 0, 'GRID', ierr)
      allocate(coefs%coef%nlte_coef%coef(coefs%coef%nlte_coef%ncoef, coefs%coef%nlte_coef%nsat, &
           coefs%coef%nlte_coef%nsol, nlte_count))
      coefs%coef%nlte_coef%coef(:,:,:,:) = bigArray(:,:,:,nlte_chans(1:nlte_count))
      deallocate(nlte_chans, bigArray)
       
    end if

    if (mmpi_myid==0 .and. associated(coefs%coef%bounds) )  deallocate(coefs%coef%bounds)
  
    !allocate bounds array to store opdep calculation layer limits
    !1st dim: upper boundary layer [ub](above which coefs all zeros), lower boundary layer [lb]
    !4th dim: thermal layer limits, solar layer limits
    allocate(coefs%coef%bounds(2, coefs%coef%fmv_gas, coefs%coef%fmv_chn, 2))
    call set_fastcoef_level_bounds(coefs%coef, coefs%coef%thermal, thermal = .true._jplm)
    ! if the solar_fast_coefficients section is not present then point the solar coefs to the thermal coefs
    if (coefs%coef%solarcoef) then
      call set_fastcoef_level_bounds(coefs%coef, coefs%coef%solar, thermal = .false._jplm)
    else
      coefs%coef%solar => coefs%coef%thermal
      coefs%coef%solar_corr => coefs%coef%thermal_corr
      coefs%coef%bounds(:,:,:,2) = coefs%coef%bounds(:,:,:,1)
    end if

    coefs%coef%ff_val_bc = any(coefs%coef%ff_bco(:) /= 0.0d0) .or. any(coefs%coef%ff_bcs(:) /= 1.d0)
    coefs%coef%ff_val_gam = any(coefs%coef%ff_gam(:) /= 1.d0)

    ! surface water reflectance for visible/near-ir channels
    if (any(coefs%coef%ss_val_chn == 2)) then
      if ( mmpi_myid==0) deallocate(coefs%coef%refl_visnir_ow, &
           coefs%coef%refl_visnir_fw, stat = err)
      allocate(coefs%coef%refl_visnir_ow(coefs%coef%fmv_chn), &
           coefs%coef%refl_visnir_fw(coefs%coef%fmv_chn), stat = err)
      call rttov_refl_water_interp(coefs%coef%ff_cwn, coefs%coef%refl_visnir_ow, coefs%coef%refl_visnir_fw)
    end if

    if (coefs%coef%pmc_shift .and. mmpi_myid > 0) then
      allocate(coefs%coef%pmc_ppmc(coefs%coef%fmv_chn), stat = err)
    else
      nullify(coefs%coef%pmc_pnominal, coefs%coef%pmc_coef, coefs%coef%pmc_ppmc)
    end if

    contains

    subroutine nullify_gas_coef_pointers(fast_coef)
      type(rttov_fast_coef), intent(inout) :: fast_coef
      nullify (fast_coef%mixedgas,&
           fast_coef%watervapour, &
           fast_coef%ozone,       &
           fast_coef%wvcont,      &
           fast_coef%co2,         &
           fast_coef%n2o,         &
           fast_coef%co,          &
           fast_coef%ch4,         &
           fast_coef%so2)
    end subroutine nullify_gas_coef_pointers

    subroutine dispatch_fast_coef(err, fast_coef, gas_ids, ncoefs, version, nlayers, ngas)
      integer,                        intent(out)   :: err
      type(rttov_fast_coef), pointer, intent(inout) :: fast_coef(:)
      integer(jpim),         intent(in)           :: gas_ids(:)
      integer(jpim),         intent(in)           :: ncoefs(:)
      integer(jpim),         intent(in)           :: version
      integer(jpim),         intent(in)           :: nlayers
      integer(jpim),         intent(in)           :: ngas

      integer(jpim) :: channelIndex, gasIndex, layerIndex, coefIndex
      real(8), allocatable :: bigArray(:,:,:,:)
      logical :: allocated0

      allocate(bigArray(countUniqueChannel,maxval(ncoefs),ngas,nlayers), stat=err )
      bigArray(:,:,:,:) = 0.0d0

      if (mmpi_myid > 0) then
        allocate (fast_coef(countUniqueChannel) )
        do channelIndex = 1, countUniqueChannel
          allocate(fast_coef(channelIndex)%gasarray(ngas) )
          do gasIndex = 1, ngas
            allocate (fast_coef(channelIndex)%gasarray(gasIndex)%coef( ncoefs(gasIndex), nlayers) )
          end do
        end do
      end if

      do channelIndex = 1, countUniqueChannel
        do gasIndex = 1, ngas
          call broadcastR82dArray( fast_coef(channelIndex)%gasarray(gasIndex)%coef )
        end do
      end do

      do channelIndex = 1, countUniqueChannel  
        do gasIndex = 1, ngas
          associated0 = associated( fast_coef(channelIndex)%gasarray(gasIndex)%coef )
          call rpn_comm_bcast(associated0, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
          if (associated0) then
            do layerIndex=1, nlayers
              do coefIndex=1, ncoefs(gasIndex)
                bigArray(channelIndex,coefIndex,gasIndex,layerIndex) = fast_coef(channelIndex)%gasarray(gasIndex)%coef(coefIndex,layerIndex)
              end do
            end do
          end if
        end do
      end do

      do channelIndex = 1, countUniqueChannel
        do gasIndex = 1, ngas
          associated0 = associated(fast_coef(channelIndex)%gasarray(gasIndex)%coef)
          if (associated0) deallocate(fast_coef(channelIndex)%gasarray(gasIndex)%coef)
        end do
        deallocate(fast_coef(channelIndex)%gasarray)
      end do
      deallocate(fast_coef)
      allocate (fast_coef(coefs%coef%fmv_chn) )
      do channelIndex = 1, coefs%coef%fmv_chn
        allocate (fast_coef(channelIndex)%gasarray(ngas) )
        call nullify_gas_coef_pointers( fast_coef(channelIndex) )
        do gasIndex = 1, ngas
          if (any( bigArray(indexchan(channelIndex),:,gasIndex,:) /= 0.) ) then
            allocate (fast_coef(channelIndex)%gasarray(gasIndex)%coef( ncoefs(gasIndex), nlayers) )
            do layerIndex=1, nlayers
              do coefIndex=1, ncoefs(gasIndex)
                fast_coef(channelIndex)%gasarray(gasIndex)%coef(coefIndex,layerIndex)  = bigArray(indexchan(channelIndex),coefIndex,gasIndex,layerIndex)
              end do
            end do
          end if
          call set_pointers(fast_coef(channelIndex), gasIndex, gas_ids(gasIndex))
        end do
      end do

      deallocate(bigArray, stat=err )

    end subroutine dispatch_fast_coef

  end subroutine tvs_rttov_read_coefs

  subroutine extractI41dArray(array,oldSize,index)
    implicit none
    integer, pointer :: array(:)
    integer, intent(in) :: oldSize
    integer, intent(in) :: index(:)
    ! Locals
    integer :: newSize, tmpI41d(oldSize), ierr, trueSize

    if (mmpi_myid == 0) then
      if (associated(array)) then
        trueSize = size(array)
      else
        trueSize = 0
      end if
    end if
    ierr = 0
    call rpn_comm_bcast(trueSize, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) 'error 1 in rpn_comm_bcast', ierr, trueSize
      call utl_abort('extractI41dArray')
    end if
    if (trueSize < 1) return
    
    if (trueSize /= oldSize) then
      write(*,*) 'extractI41dArray: should not happen ', trueSize, oldSize
    end if
    
    newSize = size( index )

    if (mmpi_myid > 0) allocate( array(oldSize))
    ierr = 0
    call rpn_comm_bcast(array, oldSize, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) 'error 2 in rpn_comm_bcast', ierr, array(:)
      call utl_abort('extractI41dArray')
    end if
    tmpI41d = array
    deallocate( array )
    allocate( array(newSize))
    array( : ) =  tmpI41d( index(:) )
  end subroutine extractI41dArray
  
  subroutine extractR81dArray(array,oldSize,index)
    implicit none
    real(8), pointer :: array(:)
    integer, intent(in) :: oldSize, index(:)
    !Locals
    integer :: newSize, ierr, trueSize
    real(8) :: tmpR81d(oldSize)
    
    if (mmpi_myid == 0) then
      if (associated(array)) then
        trueSize = size(array)
      else
        trueSize = 0
      end if
    end if
    ierr = 0
    call rpn_comm_bcast(trueSize, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) 'error 1 in rpn_comm_bcast', ierr, trueSize
      call utl_abort('extractR81dArray')
    end if
    if (trueSize < 1) return
    
    if (trueSize /= oldSize) then
      write(*,*) 'extractR81dArray: should not happen ', trueSize, oldSize
    end if
    
    newSize = size( index )

    if (mmpi_myid > 0) allocate( array(oldSize))
    ierr = 0
    call rpn_comm_bcast(array, oldSize, 'MPI_REAL8', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) 'error 2 in rpn_comm_bcast', ierr, array(:)
      call utl_abort('extractR81dArray')
    end if
    tmpR81d = array
    deallocate( array )
    allocate( array(newSize))
    array( : ) =  tmpR81d( index(:) )
  end subroutine extractR81dArray

  subroutine extractR82dArray(array,oldSize1,oldSize2,index)
    !second dimension is for channels
    implicit none
    real(8), pointer :: array(:,:)
    integer, intent(in) :: oldSize1, oldSize2,index(:)
    !Locals
    integer :: newSize, ierr, trueSize,i
    real(8) :: tmpR82d(oldSize1,oldsize2)
    
    if (mmpi_myid == 0) then
      if (associated(array)) then
        trueSize = size(array)
      else
        trueSize = 0
      end if
    end if
    ierr = 0
    call rpn_comm_bcast(trueSize, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) 'error 1 in rpn_comm_bcast', ierr, trueSize
      call utl_abort('extractR82dArray')
    end if
    if (trueSize < 1) return
    
    if (trueSize /= oldSize1 * oldSize2) then
      write(*,*) 'extractR82dArray: should not happen ', trueSize, oldSize1, oldSize2
    end if

    newSize = size( index )

    if (mmpi_myid > 0) allocate( array(oldSize1,oldSize2) )
    ierr = 0
    call rpn_comm_bcast(array, oldSize1*oldSize2, 'MPI_REAL8', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) 'error 2 in rpn_comm_bcast', ierr, array(:,:)
      call utl_abort('extractR82dArray')
    end if
    tmpR82d = array
    deallocate( array )
    allocate( array(oldSize1,newSize))
    do i=1, newSize
      array( :,i) =  tmpR82d(:, index(i) )
    end do
  end subroutine extractR82dArray

  subroutine extractR83dArray(array,oldSize1,oldSize2,oldSize3,index)
    !second dimension is for channels
    implicit none
    real(8), pointer :: array(:,:,:)
    integer, intent(in) :: oldSize1, oldSize2,oldSize3,index(:)
    !Locals
    integer :: newSize, ierr, trueSize,i
    real(8) :: tmpR83d(oldSize1,oldSize2,oldSize3)
    
    if (mmpi_myid == 0) then
      if (associated(array)) then
        trueSize = size(array)
      else
        trueSize = 0
      end if
    end if
    ierr = 0
    call rpn_comm_bcast(trueSize, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) 'error 1 in rpn_comm_bcast', ierr, trueSize
      call utl_abort('extractR83dArray')
    end if
    if (trueSize < 1) return

    if (trueSize /= oldSize1 * oldSize2 * oldSize3) then
      write(*,*) 'extractR83dArray: should not happen ', trueSize, oldSize1, oldSize2, oldSize3
    end if
  
    newSize = size( index )
  
    if (mmpi_myid > 0) allocate( array(oldSize1,oldSize2, oldSIze3) )
    ierr = 0
    call rpn_comm_bcast(array, oldSize1*oldSize2*oldSize3, 'MPI_REAL8', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) 'error 2 in rpn_comm_bcast', ierr, array(:,:,:)
      call utl_abort('extractR83dArray')
    end if
    tmpR83d = array
    deallocate( array )
    allocate( array(oldSize1,newSize,oldSize3))
    do i=1, newSize
      array( :,i,:) =  tmpR83d(:, index(i),: )
    end do
  end subroutine extractR83dArray

  subroutine extractCmplx81dArray(array,oldSize,index)
    implicit none
    complex(kind=8), pointer :: array(:)
    integer, intent(in) :: oldSize, index(:)
    !Locals
    integer :: newSize, ierr, trueSize
    complex(kind=8) :: tmpCx81d(oldSize)

    if (mmpi_myid == 0) then
      if (associated(array)) then
        trueSize = size(array)
      else
        trueSize = 0
      end if
    end if
    ierr = 0
    call rpn_comm_bcast(trueSize, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) 'error 1 in rpn_comm_bcast', ierr, trueSize
      call utl_abort('extractCmplx81dArray')
    end if
    if (trueSize < 1) return
  
    if (trueSize /= oldSize) then
      write(*,*) 'extractCmplx81dArray: should not happen ', trueSize, oldSize
    end if

    newSize = size( index )
    
    if (mmpi_myid > 0) allocate( array(oldSize))
    ierr = 0
    call rpn_comm_bcast(array, oldSize, 'MPI_COMPLEX8', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) 'error 2 in rpn_comm_bcast', ierr, array(:)
      call utl_abort('extractCmplx81dArray')
    end if
    tmpCx81d = array
    deallocate( array )
    allocate( array(newSize))
    array( : ) =  tmpCx81d( index(:) )
  end subroutine extractCmplx81dArray
  

  subroutine broadcastR82dArray(array)
    implicit none
    real(kind=8), pointer :: array(:,:)
    !Locals
    logical :: associated0
    integer :: ierr

    associated0 = associated(array)
    call rpn_comm_bcast(associated0, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) 'error 1 in rpn_comm_bcast', ierr, associated0
      call utl_abort('broadcastR82dArray')
    end if
    ierr = 0
    if (associated0) call rpn_comm_bcast(array, size(array) , 'MPI_REAL8', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) 'error 2 in rpn_comm_bcast', ierr, size(array,dim=1), size(array,dim=2)
      call utl_abort('broadcastR82dArray')
    end if
   
  end subroutine broadcastR82dArray

  subroutine broadcastR81dArray(array)
    implicit none
    real(kind=8), pointer :: array(:)
    !Locals
    logical :: associated0
    integer :: ierr
    
    associated0 = associated(array)
    call rpn_comm_bcast(associated0, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    if (ierr/=0) then
      write(*,*) 'error 1 in rpn_comm_bcast', ierr, associated0
      call utl_abort('broadcastR81dArray')
    end if
    ierr = 0
    if (associated0) call rpn_comm_bcast(array, size(array) , 'MPI_REAL8', 0, 'GRID', ierr)
    if (ierr/=0) then
      write(*,*) 'error 2 in rpn_comm_bcast', ierr, size(array)
      call utl_abort('broadcastR81dArray')
    end if
    
  end subroutine broadcastR81dArray


  subroutine broadcastI41dArray(array)
    implicit none
    integer(kind=4), pointer :: array(:)
    !Locals
    logical :: associated0
    integer :: ierr

    associated0 = associated(array)
    call rpn_comm_bcast(associated0, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    if (ierr/=0) then
      write(*,*) 'error 1 in rpn_comm_bcast', ierr, associated0
      call utl_abort('broadcastI41dArray')
    end if
    ierr = 0
    if (associated0) call rpn_comm_bcast(array, size(array) , 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr/=0) then
      write(*,*) 'error 2 in rpn_comm_bcast', ierr, size(array)
      call utl_abort('broadcastI41dArray')
    end if

  end subroutine broadcastI41dArray

  !--------------------------------------------------------------------------
  !   tvs_printDetailledOmfStatistics
  !--------------------------------------------------------------------------
  subroutine tvs_printDetailledOmfStatistics(obsSpaceData)
    !
    ! :Purpose: Print channel by channnel O-F statistics fro radiances
    !
    implicit none

    !Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData! obsSpacaData structure

    ! Locals:
    integer :: sensorIndex, channelIndex, tovsIndex
    real(8) zjoch  (0:tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    real(8) zavgnrm(0:tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    real(pre_obsReal) :: zdtb, obsPRM
    integer nchanperline, startChannel, endChannel
    integer count, incanjo
    integer idatyp
    integer rttovChannelNumber, bufrChannelNumber
    integer inobsch(0:tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    integer lcanjo(tvs_maxChannelNumber)
    integer :: headerIndex, bodyIndex
    real(8) :: sigmaObs

    write(*,*) 'tvs_printDetailledOmfStatistics: Starting'

    if ( tvs_nobtov == 0) return    ! exit if there are not tovs data

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
      if ( .not. tvs_isIdBurpTovs(idatyp) ) then
        write(*,*) 'tvs_printDetailledOmfStatistics: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER
      end if
      tovsIndex = tvs_tovsIndex(headerIndex)
      if ( tovsIndex == -1 ) cycle HEADER
       
      sensorIndex = tvs_lsensor(tovsIndex)

      ! Set the body list
      ! (& start at the beginning of the list)
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      count = 0
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY
        
        ! Only consider if flagged for assimilation
        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY                

        call tvs_getChannelNumIndexFromPPP( obsSpaceData, headerIndex, bodyIndex, &
                                            rttovChannelNumber, channelIndex )
        bufrChannelNumber = rttovChannelNumber + tvs_channelOffset(sensorIndex)
        if ( channelIndex == 0 ) then
          write(*,'(A)') '  tvs_printDetailledOmfStatistics: error with channel number'
          call utl_abort(' tvs_printDetailledOmfStatistics')
        end if

        zdtb = obs_bodyElem_r(obsSpaceData,OBS_PRM,bodyIndex) - &
             tvs_radiance (tovsIndex) % bt(channelIndex)
        if ( tvs_debug ) then
          obsPRM = obs_bodyElem_r(obsSpaceData,OBS_PRM,bodyIndex)
          write(*,'(a,i4,2f8.2,f6.2)') ' rttovChannelNumber,sim,obs,diff= ', &
               rttovChannelNumber,  tvs_radiance (tovsIndex) % bt(channelIndex), &
               obsPRM, -zdtb
        end if

        sigmaObs = obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)

        if ( sigmaObs == MPC_missingValue_R8) cycle body

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
    if ( count > 0 ) then
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
          if ( inobsch(channelIndex, sensorIndex) /= 0 ) then
            incanjo = incanjo + 1
            lcanjo(incanjo) = channelIndex
          end if
        end do
        if ( incanjo /= 0 ) then
          write(*,"(/1x,'sensor #',i2,'. platform: ',a, 'instrument: ',a)") &
               sensorIndex, tvs_satelliteName(sensorIndex), tvs_instrumentName(sensorIndex)
          do startChannel = 1, incanjo, nchanperline
            endChannel = min(startChannel + nchanperline - 1 , incanjo)
            if ( startChannel == 1 ) then
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

    !Arguments:
    integer, intent(in)  :: idsat            ! Satellite index
    integer, intent(out) :: channelIndex_out ! Channel index
    integer, intent(in)  :: channelNumber_in ! Channel number

    ! Locals:
    logical, save              :: first =.true.
    integer                    :: channelNumber, sensorIndex, channelIndex 
    integer, allocatable, save :: savedChannelIndexes(:,:)

    if (first) then
      allocate( savedChannelIndexes(tvs_nsensors, tvs_maxChannelNumber ) )
      savedChannelIndexes(:,:) = -1
      do sensorIndex = 1, tvs_nsensors
        channels:do channelNumber = 1,  tvs_maxChannelNumber
          indexes: do channelIndex =1, tvs_nchan(sensorIndex)
            if ( channelNumber == tvs_ichan(channelIndex,sensorIndex) ) then
              savedChannelIndexes(sensorIndex,channelNumber) = channelIndex
              exit indexes
            end if
          end do indexes
        end do channels
      end do
      first = .false.
    end if

    channelIndex_out = savedChannelIndexes(idsat,channelNumber_in)

    if (channelIndex_out == -1) then
      write(*,*) 'channel number requested = ', channelNumber_in
      write(*,*) 'idsat = ', idsat
      write(*,*) 'tvs_getLocalChannelIndexFromChannelNumber: warning channel not found'  
    end if

  end subroutine tvs_getLocalChannelIndexFromChannelNumber


  subroutine updateCloudInTovsProfile(sensorTovsIndexes, nlv_T, mode, beSilent)
    !
    ! :Purpose: Modify the cloud in tvs_profiles_nl structure of rttov.
    !
    implicit none
    
    ! Arguments:
    integer,      intent(in) :: sensorTovsIndexes(:)
    integer,      intent(in) :: nlv_T
    character(*), intent(in) :: mode         ! save or restore
    logical,      intent(in) :: beSilent     ! flag to control verbosity

    ! Locals:
    integer :: profileIndex, profileCount
    real(8), allocatable, save :: cloudProfileToStore(:,:)

    if ( .not. beSilent ) write(*,*) 'updateCloudInTovsProfile: Starting'
    if ( .not. beSilent ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    profileCount = size(sensorTovsIndexes)

    if ( trim(mode) == 'save' ) then 
      if (allocated(cloudProfileToStore)) deallocate(cloudProfileToStore)
      allocate(cloudProfileToStore(nlv_T,profileCount))

      do profileIndex = 1, profileCount
        cloudProfileToStore(:,profileIndex) = tvs_profiles_nl(sensorTovsIndexes(profileIndex)) % clw(:)
        tvs_profiles_nl(sensorTovsIndexes(profileIndex)) % clw(:) = qlim_getMinValueCloud('LWCR') 
      end do

    else if ( trim(mode) == 'restore' ) then 
      do profileIndex = 1, profileCount
        tvs_profiles_nl(sensorTovsIndexes(profileIndex)) % clw(:) = cloudProfileToStore(:,profileIndex)
      end do

      deallocate(cloudProfileToStore)

    else
      call utl_abort('updateCloudInTovsProfile: mode should be either "save" or "restore"')

    end if

  end subroutine updateCloudInTovsProfile


  subroutine updateCloudInTovsCloudProfile(sensorTovsIndexes, nlv_T, mode, beSilent)
    !
    ! :Purpose: Modify the cloud in tvs_cld_profiles_nl structure of rttovScatt.
    !
    implicit none
    
    ! Arguments:
    integer,      intent(in) :: sensorTovsIndexes(:)
    integer,      intent(in) :: nlv_T
    character(*), intent(in) :: mode         ! save or restore
    logical,      intent(in) :: beSilent     ! flag to control verbosity

    ! Locals:
    integer :: profileIndex, profileCount
    real(8), allocatable, save :: rainFluxProfileToStore(:,:)
    real(8), allocatable, save :: snowFluxProfileToStore(:,:)
    real(8), allocatable, save :: clwProfileToStore(:,:)
    real(8), allocatable, save :: ciwProfileToStore(:,:)
    real(8), allocatable, save :: cloudFractionProfileToStore(:,:)

    if ( .not. beSilent ) write(*,*) 'updateCloudInTovsCloudProfile: Starting'
    if ( .not. beSilent ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    profileCount = size(sensorTovsIndexes)

    if ( trim(mode) == 'save' ) then 
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
        rainFluxProfileToStore(:,profileIndex) = tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro(:,1)
        snowFluxProfileToStore(:,profileIndex) = tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro(:,2)
        clwProfileToStore(:,profileIndex) = tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro(:,4)
        ciwProfileToStore(:,profileIndex) = tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro(:,5)
        cloudFractionProfileToStore(:,profileIndex) = tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro_frac(:,1)
        tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro(:,1) = qlim_getMinValueCloud('RF')
        tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro(:,2) = qlim_getMinValueCloud('SF')
        tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro(:,4) = qlim_getMinValueCloud('LWCR')
        tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro(:,5) = qlim_getMinValueCloud('IWCR')
        tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro_frac(:,1) = qlim_getMinValueCloud('CLDR')
      end do

    else if ( trim(mode) == 'restore' ) then 
      do profileIndex = 1, profileCount
        tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro(:,1) = rainFluxProfileToStore(:,profileIndex)
        tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro(:,2) = snowFluxProfileToStore(:,profileIndex)
        tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro(:,4) = clwProfileToStore(:,profileIndex)
        tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro(:,5) = ciwProfileToStore(:,profileIndex)
        tvs_cld_profiles_nl(sensorTovsIndexes(profileIndex)) % hydro_frac(:,1) = cloudFractionProfileToStore(:,profileIndex)
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
    type(struct_obs),   intent(in) :: obsSpaceData
    integer,            intent(in) :: headerIndex
    integer,            intent(in) :: bodyIndex
    integer,           intent(out) :: channelNumber
    integer,           intent(out) :: channelIndex

    ! Locals:
    integer :: tovsIndex, sensorIndex

    tovsIndex = tvs_tovsIndex(headerIndex)
    sensorIndex = tvs_lsensor(tovsIndex)

    channelNumber = nint(obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex))
    channelNumber = max( 0 , min( channelNumber , tvs_maxChannelNumber + 1))
    channelNumber = channelNumber - tvs_channelOffset(sensorIndex)
    channelIndex = utl_findArrayIndex(tvs_ichan(:,sensorIndex),tvs_nchan(sensorIndex),channelNumber)

  end subroutine tvs_getChannelNumIndexFromPPP

end module tovs_nl_mod

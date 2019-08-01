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
  ! MODULE tovs_nl (prefix='tvs' category='4. Observation operators')
  !
  ! :Purpose: Derived types, public variables and procedures related to the
  !           nonlinear version of RTTOV
  !

  use rttov_interfaces_mod
  use rttov_types, only :   &
       rttov_coefs         ,&
       rttov_options       ,&
       rttov_profile       ,&
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
       ngases_max          ,&
       gas_id_mixed        ,&
       gas_unit_specconc   ,&
       zenmax

  use parkind1, only : jpim, jplm

  use rttov_fast_coef_utils_mod, only: set_pointers, set_fastcoef_level_bounds

  use mpi_mod
  use codtyp_mod
  use mpi
  use utilities_mod
  use obsSpaceData_mod
  use EarthConstants_mod
  use MathPhysConstants_mod
  use ozoneclim_mod
  use columnData_mod 
  use presProfileOperators_mod
  use tovs_extrap_mod
  use mod_rttov_emis_atlas
  use rMatrix_mod 
  use verticalCoord_mod

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

  ! public derived type through inheritance (from module rttov_types)
  public :: rttov_radiance, rttov_profile, rttov_chanprof, rttov_coefs, rttov_transmission, rttov_options, rttov_emissivity

  ! public variables (parameters)
  public :: tvs_maxChannelNumber, tvs_maxNumberOfChannels, tvs_maxNumberOfSensors,  tvs_nlevels, tvs_mesosphereLapseRate
  ! public variables (non-parameters)
  public :: tvs_nchan, tvs_ichan, tvs_lsensor, tvs_headerIndex, tvs_tovsIndex, tvs_nobtov
  public :: tvs_isReallyPresent,tvs_listSensors
  public :: tvs_nsensors, tvs_platforms, tvs_satellites, tvs_instruments, tvs_channelOffset
  public :: tvs_debug, tvs_satelliteName, tvs_instrumentName
  public :: platform_name, inst_name ! (from rttov)
  public :: tvs_coefs, tvs_opts, tvs_profiles, tvs_transmission,tvs_emissivity
  public :: tvs_radiance, tvs_surfaceParameters

  ! public procedures
  public :: tvs_fillProfiles, tvs_rttov, tvs_calc_jo, tvs_allocTransmission
  public :: tvs_setupAlloc,tvs_setup, tvs_isIdBurpTovs, tvs_isIdBurpInst
  public :: tvs_getInstrumentId, tvs_getPlatformId, tvs_mapSat, tvs_mapInstrum
  public :: tvs_isInstrumHyperSpectral, tvs_getChanprof, tvs_countRadiances
  public :: tvs_getHIREmissivities, tvs_getOtherEmissivities, tvs_rttov_read_coefs, tvs_getChannelIndexFromChannelNumber
  ! Module parameters
  ! units conversion from  mixing ratio to ppmv and vice versa
  real(8), parameter :: qMixratio2ppmv  = (1000000.0d0 * mair) / mh2o
  real(8), parameter :: qppmv2Mixratio  = mh2o / (1000000.0d0 * mair)
  real(8), parameter :: o3Mixratio2ppmv = (1000000.0d0 * mair) / mo3
  real(8), parameter :: o3ppmv2Mixratio = mo3 / (1000000.0d0 * mair)

  integer, parameter :: tvs_maxChannelNumber   = 8461   ! Max. value for channel number
  integer, parameter :: tvs_maxNumberOfChannels = 2211  ! Max. no. of channels (for one profile/spectra)
  integer, parameter :: tvs_maxNumberOfSensors  = 40    ! Max no sensors to be used
  integer, parameter :: tvs_nlevels     = 101           ! Maximum No. of RTTOV pressure levels including "rttov top" at 0.005 hPa

  ! S. Heilliette this parameter was computed from the mean lapse rate between 50 km and 85 km
  ! of the US standard atmosphere from data contained in "AFGL Atmospheric Constituent Profiles (0-120km)"
  ! afgl 1986, g.p. anderson, j.h. chetwynd, s.a. clough, e. p. shettle and f.x. kneizys
  ! unit is K/log(P), a positive value corresponds to decrease of temperature with height
  real(8), parameter :: tvs_mesosphereLapseRate=16.2d0

  ! Module variables
  integer, allocatable :: tvs_nchan(:)             ! Number of channels per instrument
  integer, allocatable :: tvs_ichan(:,:)           ! List of channels per instrument
  integer, allocatable :: tvs_lsensor(:)           ! Sensor number for each profile
  integer, allocatable :: tvs_headerIndex (:)      ! Observation position in obsSpaceData header for each profile
  integer, allocatable :: tvs_tovsIndex (:)          ! Index in TOVS structures for each observation header in obsSpaceData
  logical, allocatable :: tvs_isReallyPresent(:)   ! Logical flag to identify instruments really assimilated
  integer, allocatable :: tvs_listSensors(:,:)     ! Sensor list
  type(surface_params), allocatable :: tvs_surfaceParameters(:) ! surface parameters 
  integer tvs_nobtov                               ! Number of tovs observations
  integer tvs_nsensors                             ! Number of individual sensors.
  integer tvs_platforms(tvs_maxNumberOfSensors)    ! RTTOV platform ID's (e.g., 1=NOAA; 2=DMSP; ...)
  integer tvs_satellites(tvs_maxNumberOfSensors)   ! RTTOV satellite ID's (e.g., 1 to 16 for NOAA; ...)
  integer tvs_instruments(tvs_maxNumberOfSensors)  ! RTTOVinstrument ID's (e.g., 3=AMSU-A; 4=AMSU-B; 6=SSMI; ...)
  integer tvs_channelOffset(tvs_maxNumberOfSensors)! BURP to RTTOV channel mapping offset
  logical tvs_debug                                ! Logical key controlling statements to be  executed while debugging TOVS only
  logical useUofWIREmiss                           ! Flag to activate use of RTTOV U of W emissivity Atlases
  character(len=15) tvs_satelliteName(tvs_maxNumberOfSensors)
  character(len=15) tvs_instrumentName(tvs_maxNumberOfSensors)
  character(len=8) radiativeTransferCode           ! RadiativeTransferCode : TOVS radiation model used
  real(8) , allocatable :: tvs_emissivity(:,:)     ! Surface emissivities organized by profiles and channels
  integer, parameter :: kslon=2160, kslat=1080     ! CERES file dimension in grid points

  ! High resolution surface fields
  integer :: surfaceType(kslon,kslat)  
  real(8) :: waterFraction(kslon,kslat) 

  ! Derived types
  type( rttov_coefs ),        allocatable :: tvs_coefs(:)          ! coefficients
  type( rttov_options ),      allocatable :: tvs_opts(:)           ! options
  type( rttov_profile ),      allocatable :: tvs_profiles(:)       ! profiles, all profiles
  type( rttov_radiance ),     allocatable :: tvs_radiance(:)       ! radiances organized by profile
  type( rttov_transmission ), allocatable :: tvs_transmission(:)   ! transmittances all profiles for HIR quality control

  integer, external :: get_max_rss
 
contains

  !--------------------------------------------------------------------------
  ! tvs_setupAlloc
  !--------------------------------------------------------------------------
  subroutine tvs_setupAlloc(lobsSpaceData)
    !
    ! :Purpose: Memory allocation for the non linear radiative transfer model variables.
    !
    implicit none

    !Parameters:
    type(struct_obs) :: lobsSpaceData

    ! Locals:
    integer :: allocStatus(6)
    integer :: satelliteCode, instrumentCode, iplatform, isat, instrum
    integer :: tovsIndex, idatyp, sensorIndex
    integer :: channelNumber, nosensor, channelIndex
    integer :: errorstatus(1)
    integer :: headerIndex, bodyIndex

    if (tvs_nsensors == 0) return

    !  1. Determine the number of radiances to be assimilated.
    !      Construct a list of channels for each sensor.
    !      Construct a list of sensor number for each profile

    write(*,*) "Entering tvs_setupAlloc" 

    allocStatus = 0
    allocate (tvs_nchan(tvs_nsensors),                         stat= allocStatus(1))
    allocate (tvs_ichan(tvs_maxNumberOfChannels,tvs_nsensors), stat= allocStatus(2))
    allocate (tvs_lsensor(obs_numheader(lobsSpaceData)),       stat= allocStatus(3))
    allocate (tvs_headerIndex (obs_numheader(lobsSpaceData)),       stat= allocStatus(4))
    allocate (tvs_tovsIndex(obs_numheader(lobsSpaceData)),       stat= allocStatus(5))
    allocate (tvs_isReallyPresent(tvs_nsensors),               stat= allocStatus(6))

    call utl_checkAllocationStatus(allocStatus, " tvs_setupAlloc")
  
    tvs_nchan(:) = 0 
    tvs_ichan(:,:) = 0
    tvs_isReallyPresent(:) = .true.
    tvs_lsensor(:) = -1
    tvs_headerIndex(:) = -1
    tvs_tovsIndex (:) = -1


    tvs_nobtov = 0

    ! Loop over all header indices of the 'TO' family
    ! Set the header list & start at the beginning of the list
    call obs_set_current_header_list(lobsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(lobsSpaceData)
      if (headerIndex < 0) exit HEADER

      idatyp = obs_headElem_i(lobsSpaceData,OBS_ITY,headerIndex)

      if ( .not. tvs_isIdBurpTovs(idatyp) ) cycle HEADER   ! Proceed to the next headerIndex

      tvs_nobtov = tvs_nobtov + 1
     
      !    Construct list of channels for each sensor:
      !          map burp satellite info to RTTOV platform and satellite.
      satelliteCode = obs_headElem_i(lobsSpaceData,OBS_SAT,headerIndex)
      call tvs_mapSat(satelliteCode,iplatform,isat)
      if (iplatform == -1) then
        write(*,*) "Unknown OBS_SAT !", satelliteCode
        call utl_abort('tvs_setupAlloc')
      end if
      !    map burp instrument info to RTTOV instrument.
      instrumentCode = obs_headElem_i(lobsSpaceData,OBS_INS,headerIndex)
      call tvs_mapInstrum(instrumentCode,instrum)
      if (instrum == -1) then
        write(*,*) "Unknown OBS_INS !", instrumentCode
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
        write(*,*) " tvs_setupAlloc: Warning Invalid Sensor ", iplatform, isat, instrum, " skipping ..."
        cycle HEADER
      end if

      ! Loop over all body indices (still in the 'TO' family)
      ! Set the body list & start at the beginning of the list
      call obs_set_current_body_list(lobsSpaceData, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(lobsSpaceData)
        if (bodyIndex < 0) exit BODY

        if ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) then
          channelNumber = nint(obs_bodyElem_r(lobsSpaceData,OBS_PPP,bodyIndex))
          channelNumber = max(0,min(channelNumber,tvs_maxChannelNumber+1))
          
          channelNumber = channelNumber - tvs_channelOffset(nosensor)

          channelIndex = utl_findArrayIndex(tvs_ichan(:,nosensor),tvs_nchan(nosensor),channelNumber)
          if ( channelIndex == 0 ) then
            tvs_nchan(nosensor) = tvs_nchan(nosensor) + 1
            tvs_ichan(tvs_nchan(nosensor),nosensor) = channelNumber
          end if
        end if
      end do BODY
    end do HEADER

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



    !  3. Initialize TOVS radiance transfer model

    if ( radiativeTransferCode == 'RTTOV' ) then

      write(*,'(//,10x,A)') "-rttov_setup: initializing the TOVS radiative transfer model"

      allocate (tvs_coefs(tvs_nsensors)          ,stat= allocStatus(1))
      allocate (tvs_listSensors (3,tvs_nsensors) ,stat= allocStatus(2))
      allocate (tvs_opts (tvs_nsensors)          ,stat= allocStatus(3))

      call utl_checkAllocationStatus(allocStatus(1:3), " tvs_setupAlloc before rttov initialization")

      do sensorIndex=1, tvs_nsensors
        tvs_listSensors(1,sensorIndex) = tvs_platforms  (sensorIndex)
        tvs_listSensors(2,sensorIndex) = tvs_satellites (sensorIndex)
        tvs_listSensors(3,sensorIndex) = tvs_instruments(sensorIndex)

        !< General configuration options
        tvs_opts(sensorIndex) % config % apply_reg_limits = .true. ! if true application of profiles limits
        tvs_opts(sensorIndex) % config % verbose = .false. ! verbose output
        tvs_opts(sensorIndex) % config % do_checkinput = .true. ! to check if input profiles are within absolute and regression limits
        !< General RT options
        tvs_opts(sensorIndex) % rt_all % switchrad = .true.  ! to use brightness temperature (true) or radiance (false) units in AD routine
        tvs_opts(sensorIndex) % rt_all % use_q2m = .false.   ! if true use of surface humidity (false for compatibility with the way rttov 8.7 was compiled)
        tvs_opts(sensorIndex) % rt_all % addrefrac = .true. ! to account for atmospheric refraction
        !< VIS/IR RT options
        tvs_opts(sensorIndex) % rt_ir % addsolar = .false.  ! to model solar component in the near IR (2000 cm-1 et plus)
        tvs_opts(sensorIndex) % rt_ir % addaerosl = .false. ! to account for scattering due to aerosols
        tvs_opts(sensorIndex) % rt_ir % addclouds = .false. ! to account for scattering due to clouds
        tvs_opts(sensorIndex) % rt_ir % ir_sea_emis_model = 2 ! ISEM (ir_sea_emis_model 1) useful for GEORAD
                                             ! 2 would select IREMIS which is more sophisticated (to try)
        tvs_opts(sensorIndex) % rt_ir % pc % ipcreg = -1         ! index of the required PC predictors... to see later
        tvs_opts(sensorIndex) % rt_ir % pc % addpc = .false.     ! to carry out principal component calculations 
        tvs_opts(sensorIndex) % rt_ir % pc % addradrec = .false. ! to reconstruct radiances from principal components
        !< MW RT options
        tvs_opts(sensorIndex) % rt_mw % clw_data = .false.  ! profil d'eau liquide pas disponible
        tvs_opts(sensorIndex) % rt_mw % fastem_version = 6  ! use fastem version 6 microwave sea surface emissivity model (1-6)
        !< Interpolation options
        tvs_opts(sensorIndex) % interpolation % addinterp = .false. ! use of internal profile interpolator (rt calculation on model levels)
        tvs_opts(sensorIndex) % interpolation % lgradp = .false.    ! allow tl/ad of user pressure levels

        tvs_opts(sensorIndex) % rt_ir % co2_data = .false.
        tvs_opts(sensorIndex) % rt_ir % n2o_data = .false.
        tvs_opts(sensorIndex) % rt_ir % co_data  = .false.
        tvs_opts(sensorIndex) % rt_ir % ch4_data = .false.

        errorstatus = errorstatus_success

        call tmg_start(83,'RTTOV_SETUP')
        call tvs_rttov_read_coefs(errorstatus(1), tvs_coefs(sensorIndex), tvs_opts(sensorIndex), tvs_ichan(1:tvs_nchan(sensorIndex),sensorIndex), tvs_listSensors(:,sensorIndex))
        if (errorstatus(1) /= errorstatus_success) then
          write(*,*) 'tvs_rttov_read_coefs: fatal error reading coefficients',errorstatus,sensorIndex,tvs_listSensors(1:3,sensorIndex)
          call utl_abort('tvs_setupAlloc')
        end if
        call tmg_stop(83)

        tvs_opts(sensorIndex) % rt_ir % ozone_data = ( tvs_coefs(sensorIndex) % coef % nozone > 0 ) ! profil d'ozone disponible

        ! Ensure the options and coefficients are consistent
        call rttov_user_options_checkinput(errorstatus(1), tvs_opts(sensorIndex), tvs_coefs(sensorIndex))
        if (errorstatus(1) /= errorstatus_success) then
          write(*,*) 'error in rttov options',errorstatus
          call utl_abort('tvs_setupAlloc')
        end if
       
      end do


      !   4. Memory allocations for radiative tranfer model variables

      !  Profiles

      allocate(tvs_profiles(tvs_nobtov) , stat=allocStatus(1) )
      call utl_checkAllocationStatus(allocStatus(1:1), " tvs_setupAlloc tvs_profiles 1")

      do tovsIndex = 1, tvs_nobtov
        sensorIndex = tvs_lsensor(tovsIndex)
        if (sensorIndex > -1) then
          ! allocate model profiles atmospheric arrays with RTTOV levels dimension
          call rttov_alloc_prof(errorstatus(1),1,tvs_profiles(tovsIndex),  tvs_coefs(sensorIndex) % coef % nlevels  , &    ! 1 = nprofiles un profil a la fois
               tvs_opts(sensorIndex),asw=1,coefs=tvs_coefs(sensorIndex),init=.false. ) ! asw =1 allocation
          call utl_checkAllocationStatus(errorstatus(1:1), " tvs_setupAlloc tvs_profiles 2")
        end if
      end do

      ! Radiance by profile

      allocate( tvs_radiance(tvs_nobtov) ,stat= allocStatus(1))

      call utl_checkAllocationStatus(allocStatus(1:1), " tvs_setupAlloc radiances 1")
  
      do tovsIndex = 1, tvs_nobtov
        sensorIndex = tvs_lsensor(tovsIndex)
        if (sensorIndex > -1) then
          ! allocate BT equivalent to total direct, tl and ad radiance output
          allocate( tvs_radiance(tovsIndex)  % bt  ( tvs_nchan(sensorIndex) ) ,stat= allocStatus(1))
          tvs_radiance(tovsIndex)  % bt  ( : ) = 0.d0
          call utl_checkAllocationStatus(allocStatus(1:1), " tvs_setupAlloc radiances 2")
        end if
      end do

    end if
  
    write(*,*) "Leaving tvs_setupAlloc"

  end subroutine tvs_setupAlloc



  !--------------------------------------------------------------------------
  ! tvs_allocTransmission
  !--------------------------------------------------------------------------
  subroutine tvs_allocTransmission
    ! :Purpose: Allocate the global rttov transmission structure used
    !           when this is needed for some purpose (e.g. used in 
    !           LETKF to determine peak pressure level of each radiance
    !           channel for vertical localization).
    !
    implicit none

    ! Locals:
    integer :: allocStatus(2), jo, isens, nc, nl

    allocStatus(:) = 0
    allocate( tvs_transmission(tvs_nobtov), stat=allocStatus(1))
    call utl_checkAllocationStatus(allocStatus(1:1), " tvs_allocTransmission")

    do jo = 1, tvs_nobtov
      isens = tvs_lsensor(jo)
      nc = tvs_nchan(isens)
      nl = tvs_coefs(isens) % coef % nlevels
      ! allocate transmittance from surface and from pressure levels
      allocate( tvs_transmission(jo) % tau_total ( nc ), stat= allocStatus(1))
      allocate( tvs_transmission(jo) % tau_levels(nl,nc), stat= allocStatus(2))
      call utl_checkAllocationStatus(allocStatus, " tvs_allocTransmission")
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
    !Locals:
    integer  sensorIndex, ierr, nulnam
    integer, external :: fclos, fnom
    integer :: nsensors
    character(len=15) :: csatid(tvs_maxNumberOfSensors), cinstrumentid(tvs_maxNumberOfSensors)
    character(len=8)  :: crtmodl
    logical :: ldbgtov

    namelist /NAMTOV/ nsensors, csatid, cinstrumentid
    namelist /NAMTOV/ ldbgtov
    namelist /NAMTOV/ useUofWIREmiss, crtmodl

 
    !   1.1 Default values for namelist variables

    nsensors = 0
    csatid(:) = '***UNDEFINED***'
    cinstrumentid(:) = '***UNDEFINED***'
    csatid(1) = 'NOAA16'
    cinstrumentid(1) = 'AMSUA'
    ldbgtov = .false.
    crtmodl = 'RTTOV'
    useUofWIREmiss = .false.

    !   1.2 Read the NAMELIST NAMTOV to modify them
 
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam, nml=namtov, iostat=ierr)
    if (ierr /= 0) call utl_abort('tvs_setup: Error reading namelist')
    if (mpi_myid == 0) write(*,nml=namtov)
    ierr = fclos(nulnam)


    !  1.3 Transfer namelist variables to module variables
 
    tvs_nsensors = nsensors
    tvs_debug = ldbgtov
    radiativeTransferCode = crtmodl
    tvs_instrumentName(:) = cinstrumentid(:)
    tvs_satelliteName(:) = csatid(:)

    !  1.4 Validate namelist values
    
    if ( tvs_nsensors == 0 ) then
      if(mpi_myid==0) then 
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

    if(mpi_myid == 0) then
      write(*,'(A)') 
      write(*,'(3X,A)') '- Parameters used for TOVS processing (read in NAMTOV)'
      write(*,'(3X,A)') '  ----------------------------------------------------'
      write(*,'(6X,A,2X,L1)') 'TOVS debug                           : ', tvs_debug
      write(*,'(6X,A,2X,L1)') 'Use of UW IR land emissivity atlases : ', useUofWIREmiss
      write(*,'(6X,A,2X,A)')  'Radiative transfer model             : ', radiativeTransferCode
      write(*,'(6X,A,2X,I3)') 'Number of sensors                    : ', tvs_nsensors
      write(*,'(6X,"Satellite ids          : ",10A10)') (tvs_satelliteName(sensorIndex), sensorIndex=1,tvs_nsensors)
      write(*,'(6X,"Instrument ids         : ",10A10)') (tvs_instrumentName(sensorIndex), sensorIndex=1,tvs_nsensors)
      write(*,'(A)') 
      write(*,'(A)') 
      write(*,'(A)') 
      write(*,'(3X,A)') '- Reading and initialization in preparation to the TOVS processing'
      write(*,'(5X,A)') '----------------------------------------------------------------'
    end if

    !  1.6 Set up platform, satellite, instrument and channel mapping

    call sensors()

  end subroutine tvs_setup

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
     ! read the namelist
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) then
        write(*,*) "Error while opening namelist file !"
        call utl_abort("sensors")
      end if
      listoffset(:) = 0
      listinstrum(:) = "XXXXXXXX"
      read(nulnam,NAMCHANOFFSET, iostat=ierr)
      if (ierr /= 0) then
        write(*,*) "Error while reading namelist file !"
        call utl_abort("sensors")
      end if
      do instrumentIndex=0, ninst - 1
        if ( listinstrum(instrumentIndex) /= "XXXXXXXX" ) then
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
          write(*,'(A,1x,i6,1x,i3,1x,i3,1x,A15)') "Problem while reading satellite number", &
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
      if ( tvs_instrumentName(sensorIndex)(1:10) == "GOESIMAGER") then !cas particulier
        tvs_instruments(sensorIndex) = inst_id_goesim
      else if ( tvs_satelliteName(sensorIndex)(1:5) == "MTSAT") then !autre cas particulier
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

    if (mpi_myid == 0) then
      write(*,*)
      write(*,'(3X,A)') '- SENSORS. Variables prepared for RTTOV-12:'
      write(*,'(3X,A)') '  ----------------------------------------'
      write(*,*)
      write(*,'(6X,A,I3)')   "Number of sensors       : ", tvs_nsensors
      write(*,'("Platform numbers        : ",6X,10I3)')  (tvs_platforms(sensorIndex), sensorIndex=1,tvs_nsensors)
      write(*,'("Satellite numbers       : ",6X,10I3)')  (tvs_satellites(sensorIndex), sensorIndex=1,tvs_nsensors)
      write(*,'("Instrument numbers      : ",6X,10I3)')  (tvs_instruments(sensorIndex), sensorIndex=1,tvs_nsensors)
      write(*,'("Channel mapping offsets : ",6X,10I3)')  (tvs_channelOffset(sensorIndex), sensorIndex=1,tvs_nsensors)
    end if

  end subroutine sensors

  !--------------------------------------------------------------------------
  !  tvs_isIdBurpTovs
  !--------------------------------------------------------------------------
  logical function tvs_isIdBurpTovs(idatyp)
    !
    ! :Purpose: Function to check if the given idatyp (a.k.a. codtyp) corresponds to a radiance
    !
    implicit none

    ! Argument:
    integer, intent(in) :: idatyp
    
    !Locals:
    logical, save :: first=.true.
    integer,save :: ninst_tovs
    integer :: nulnam, ierr, instrumentIndex 
    integer,external :: fnom, fclos
    integer, save :: list_inst(ninst)
    character (len=22) :: inst_names(ninst)
    namelist /namtovsinst/ inst_names

    if (tvs_nsensors == 0) then
      ! no tovs data will be read, therefore false
      tvs_isIdBurpTovs = .false.
      return
    end if

    if (first) then
      nulnam = 0
      ninst_tovs = 0
      list_inst(:) = -1
      inst_names(:) = "XXXXXX"
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam, nml=namtovsinst, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isIdBurpTovs: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=namtovsinst)
      ierr = fclos(nulnam)
      do instrumentIndex=1, ninst
        if (inst_names(instrumentIndex) == "XXXXXX" ) then
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

    do instrumentIndex=1, ninst_tovs
      if (idatyp == list_inst(instrumentIndex) ) then
        tvs_isIdBurpTovs = .true.
        exit
      end if
    end do

  end function tvs_isIdBurpTovs

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
    character(len=*),intent(in) :: name !Platform name

    !Locals:
    integer           :: platformIndex, length, ipos
    character(len=64) :: tempo_name

    tvs_getPlatformId = -1
    length = len_trim(name)
    call up2low(name(1:length),tempo_name(1:length))

    if ( index(tempo_name(1:length),"npp") /= 0 ) then
      tvs_getPlatformId = platform_id_jpss
    else if ( index(tempo_name(1:length),"hmwari") /= 0 ) then
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
    character(len=*),intent(in) :: name ! Instrument name

    !Locals:
    integer           :: instrumentIndex, length
    character(len=64) :: tempo_name

    tvs_getInstrumentId = -1
    length = len_trim(name)
    call up2low(name(1:length),tempo_name(1:length))
    if ( trim(tempo_name(1:length)) == "goesim" ) then
      tvs_getInstrumentId = inst_id_goesim
    else if ( trim(tempo_name(1:length)) == "gmsim" ) then
      tvs_getInstrumentId = inst_id_gmsim
    else if ( trim(tempo_name(1:length)) == "mtsatim" ) then
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
    integer,intent(in) :: instrum     ! input Rttov instrument code

    ! Locals:
    integer ,parameter :: maxsize = 100
    integer :: nulnam, ierr, instrumentIndex 
    integer, save :: list_inst(maxsize), ninst_hir
    logical, save :: first = .true.
    integer, external :: fclos, fnom
    character (len=7) :: name_inst(maxsize)
    namelist /NAMHYPER/ name_inst

    if (first) then
      nulnam = 0
      ninst_hir = 0
      name_inst(:) = "XXXXXXX"
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namhyper, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isInstrumHyperSpectral: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=namhyper)
      ierr = fclos(nulnam)
      list_inst(:) = -1
      do instrumentIndex=1, maxsize
        list_inst(instrumentIndex) = tvs_getInstrumentId( name_inst(instrumentIndex) )
        if (name_inst(instrumentIndex) /= "XXXXXXX") then
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
        write(*,*) "tvs_isInstrumHyperSpectral: Warning : empty namhyper namelist !"
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
  !  tvs_isInstrumGeostationary
  !--------------------------------------------------------------------------
  logical function tvs_isInstrumGeostationary(instrum)
    !
    ! :Purpose: given an RTTOV instrument code return if it is a Geostationnary Imager
    !           information from namelist NAMGEO
    !
    implicit none

    ! Argument:
    integer,intent(in) :: instrum ! input Rttov instrument code

    ! Locals:
    integer ,parameter :: maxsize = 100
    integer :: nulnam, ierr, instrumentIndex 
    integer, save :: list_inst(maxsize), ninst_geo
    logical, save :: first = .true.
    character (len=8) :: name_inst(maxsize)
    integer, external :: fnom, fclos

    namelist /NAMGEO/ name_inst
    if (first) then
      nulnam = 0
      ninst_geo = 0
      name_inst(:) = "XXXXXX"
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namgeo, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isInstrumGeostationary: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=namgeo)
      ierr = fclos(nulnam)
      list_inst(:) = -1
      do instrumentIndex=1, maxsize
        list_inst(instrumentIndex) = tvs_getInstrumentId( name_inst(instrumentIndex) )
        if (name_inst(instrumentIndex) /= "XXXXXX") then
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
        write(*,*) "tvs_isInstrumGeostationary: Warning : empty namgeo namelist !"
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

    ! Arguments
    integer,intent(in)  :: instrumburp  ! burp satellite instrument (element #2019)
    integer,intent(out) :: instrum      ! RTTOV-7 instrument ID numbers (e.g. 3 for  AMSUA)
    
    ! locals  
    integer instrumentIndex, numinstburp
    integer,parameter :: mxinstrumburp   = 100
    integer,save ::   listburp(mxinstrumburp)
    character(len=8),save :: listinstrum(mxinstrumburp)
    namelist /NAMINST/ listburp, listinstrum
    logical,save :: first = .true.
    integer :: nulnam, IER
    integer, external :: fnom, fclos

    !      Table of BURP satellite sensor identifier element #002019

    !   1.0 Find instrument

    if (first) then
      ! set the default values
      listburp(:) = -1
      listinstrum(:) = "XXXXXXXX"

      ! read the namelist
      nulnam = 0
      IER = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (IER /= 0) then
        write(*,*) "Error while opening namelist file !"
        call utl_abort("tvs_mapInstrum")
      end if
      read(nulnam,NAMINST,iostat=ier)
      if (IER /= 0) then
        write(*,*) "Error while reading namelist file !"
        call utl_abort("tvs_mapInstrum")
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
    integer,intent(in)  :: isatburp   ! BURP satellite identifier
    integer,intent(out) :: iplatform  ! RTTOV-7 platform ID numbers (e.g. 1 for  NOAA)
    integer,intent(out) :: isat       ! RTTOV-7 satellite ID numbers (e.g. 15)

    ! Locals:
    integer           :: satelliteIndex, ierr, nulnam
    logical, save     :: first=.true.
    integer, external :: fnom, fclos
    integer, parameter:: mxsatburp = 100
    integer,save      :: numsatburp
    integer,save :: listburp(mxsatburp)         ! Table of BURP satellite identifier element #001007
    character(len=8),save :: listplat(mxsatburp)! Table of RTTOV platform identifier
    integer,save :: listsat (mxsatburp)         ! Table of RTTOV satellite identifier

    namelist /NAMSAT/ listburp, listplat, listsat

    !     Fill tables from namelist at the first call 
    if (first) then
      ! set the default values
      listburp(:) = -1
      listsat(:) = -1
      listplat(:) = "XXXXXXXX"
      ! read the namelist
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) then
        write(*,*) "Error while opening namelist file !"
        call utl_abort("tvs_mapSat")
      end if
      read(nulnam, NAMSAT, iostat = ierr)
      if (ierr /= 0) then
        write(*,*) "Error while reading namelist file !"
        call utl_abort("tvs_mapSat")
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
  subroutine tvs_getChanprof(sensor_id, iptobs, ObsSpaceData, chanprof, iptobs_cma_opt)
    ! 
    ! :Purpose: subroutine to initialize the chanprof structure used by RTTOV
    !
    implicit none

    ! Arguments:
    integer, intent(in)              :: sensor_id
    integer, intent(in)              :: iptobs(:)
    type(struct_obs), intent(in)     :: ObsSpaceData
    type(rttov_chanprof), intent(out):: chanprof(:)
    integer, intent(out),optional    :: iptobs_cma_opt(:)

    ! Locals:
    integer :: count, profileIndex, headerIndex, istart, iend, bodyIndex, channelNumber, nrank, iobs

    ! Build the list of channels/profiles indices
    count = 0
         
    do profileIndex = 1, size(iptobs)
      iobs = iptobs(profileIndex)
      headerIndex = tvs_headerIndex(iobs)
      if (headerIndex > 0) then
        istart = obs_headElem_i(ObsSpaceData,OBS_RLN,headerIndex)
        iend= obs_headElem_i(ObsSpaceData,OBS_NLV,headerIndex) + istart - 1
        do bodyIndex = istart, iend
          if (obs_bodyElem_i(ObsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) then
            channelNumber = nint(obs_bodyElem_r(ObsSpaceData,OBS_PPP,bodyIndex))
            channelNumber = max(0, min(channelNumber,tvs_maxChannelNumber + 1))
            channelNumber = channelNumber - tvs_channelOffset(sensor_id)
            do nrank = 1, tvs_nchan(sensor_id)
              if ( channelNumber == tvs_ichan(nrank,sensor_id) ) exit
            end do
            if (nrank /= tvs_nchan(sensor_id)+1) then
              count =  count + 1
              chanprof(count)%prof = profileIndex
              chanprof(count)%chan = nrank
              if (present(iptobs_cma_opt)) iptobs_cma_opt(count) = bodyIndex
            else
              write(*,*) "strange channel number",channelNumber
            end if
          end if
        end do
      end if
    end do
  
  end subroutine tvs_getChanprof

  !--------------------------------------------------------------------------
  !  tvs_countRadiances
  !--------------------------------------------------------------------------
  integer function tvs_countRadiances(iptobs, ObsSpaceData, assim_flag_val_opt)
    !
    ! :Purpose: to count radiances selected for assimilation
    !
    implicit none
    integer, intent(in)          :: iptobs(:)
    type(struct_obs)             :: ObsSpaceData
    integer, intent(in),optional :: assim_flag_val_opt
    

    integer :: profileIndex, headerIndex, istart, iend, bodyIndex, iobs, assim_flag_val

    if (present(assim_flag_val_opt)) then
      assim_flag_val = assim_flag_val_opt
    else
      assim_flag_val = 1
    end if

    tvs_countRadiances = 0
    do profileIndex = 1, size(iptobs)
      iobs = iptobs(profileIndex)
      headerIndex = tvs_headerIndex(iobs)
      if (headerIndex > 0) then
        istart = obs_headElem_i(ObsSpaceData,OBS_RLN,headerIndex)
        iend = obs_headElem_i(ObsSpaceData,OBS_NLV,headerIndex) + istart - 1
        do bodyIndex = istart, iend
          if(obs_bodyElem_i(ObsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) tvs_countRadiances  = tvs_countRadiances + 1
        end do
      end if
    end do

  end function tvs_countRadiances

  !--------------------------------------------------------------------------
  !  tvs_getHIREmissivities
  !--------------------------------------------------------------------------
  subroutine tvs_getHIREmissivities(iptobs, ObsSpaceData, surfem)
    !
    ! :Purpose: to get emissivity for Hyperspectral Infrared Sounders (AIRS, IASI, CrIS, ...)
    !
    implicit none

    !Arguments:
    integer, intent(in)          :: iptobs(:)
    type(struct_obs), intent(in) :: ObsSpaceData
    real(8), intent(out)         :: surfem(:)

    integer :: count, profileIndex, iobs, istart, iend, bodyIndex, headerIndex

    count = 0 
    surfem(:) = 0.98d0
    do profileIndex = 1, size(iptobs)
      iobs = iptobs(profileIndex)
      headerIndex = tvs_headerIndex(iobs)
      if (headerIndex > 0 ) then
        istart = obs_headElem_i(ObsSpaceData,OBS_RLN,headerIndex)
        iend = obs_headElem_i(ObsSpaceData,OBS_NLV,headerIndex) + istart - 1
        do bodyIndex = istart, iend
          if(obs_bodyElem_i(ObsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) then
            count = count + 1
            surfem ( count ) = obs_bodyElem_r(ObsSpaceData,OBS_SEM,bodyIndex)
          end if
        end do
      end if
    end do

  end subroutine tvs_getHIREmissivities

  !--------------------------------------------------------------------------
  !  tvs_getOtherEmissivities
  !--------------------------------------------------------------------------
  subroutine tvs_getOtherEmissivities(chanprof, iptobs, sensor_type, instrument, surfem, calcemis)
    !
    ! :Purpose: to get emissivity for microwave sounders ans infrared geostationary imagers
    !
    implicit none

    ! Arguments:
    type(rttov_chanprof),intent(in) :: chanprof(:)
    integer, intent(in)             :: iptobs(:)
    integer, intent(in)             :: sensor_type
    integer, intent(in)             :: instrument
    real(8), intent(out)            :: surfem(:)
    logical, intent(out)            :: calcemis(:)
    
    ! Locals:
    integer :: radiance_index, profileIndex, iobs, surface_type

    do radiance_index = 1, size(chanprof)
      profileIndex = chanprof(radiance_index)%prof
      iobs = iptobs(profileIndex)
      surface_type = tvs_profiles(iobs) % skin % surftype
      if     (sensor_type == sensor_id_mw ) then
        if ( surface_type == surftype_land .or. &
             surface_type == surftype_seaice     ) then
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
        write(*,*) sensor_type,instrument
        call utl_abort('tvs_getOtherEmissivities. invalid sensor type or unknown IR instrument')
      end if
    end do
   
  end subroutine tvs_getOtherEmissivities

  !--------------------------------------------------------------------------
  !  tvs_fillProfiles
  !--------------------------------------------------------------------------
  subroutine tvs_fillProfiles(columnghr,lobsSpaceData,datestamp,limlvhu,beSilent)
    !
    ! :Purpose:  to fill in tvs_profiles structure before call to non-linear, 
    !            tangent-linear or adjoint of RTTOV
    !
    implicit none

    ! Arguments:
    type(struct_columnData),intent(in) :: columnghr    ! Column structure
    type(struct_obs),       intent(in) :: lobsSpaceData! obsSpaceData structure
    integer,                intent(in) :: datestamp    ! CMC date stamp
    real(8),                intent(in) :: limlvhu      ! humidity value in the stratosphere (if extrapolated)
    logical,                intent(in) :: beSilent     ! To control verbosity

    ! Locals:
    logical :: diagTtop,TopAt10hPa
    integer :: ksurf, jpmotop, jpmolev 
    integer :: instrum, iplatform
    integer :: nlevels,nobmax
    integer :: sensorIndex, tovsIndex
    integer :: profileCount, headerIndex
    integer :: profileIndex, levelIndex
    integer :: ilowlvl_M,ilowlvl_T,nlv_M,nlv_T
    integer :: status, Vcode
    integer :: ierr,day,month,year,ijour,itime
    integer :: allocStatus(15)
    
    integer,external ::  omp_get_num_threads
    integer,external ::  newdate

    integer, allocatable :: iptobs    (:) 
    integer, allocatable :: iptobscma (:) 
  
    type(struct_vco), pointer :: vco

    real(8), allocatable :: to    (:,:)
    real(8), allocatable :: huo   (:,:)
    real(8), allocatable :: loghuo(:,:)
    real(8), allocatable :: logzhu(:,:)
    real(8), allocatable :: toext (:,:)
    real(8), allocatable :: qoext (:,:)
    real(8), allocatable :: zvlev (:,:)
    real(8), allocatable :: zt    (:,:)
    real(8), allocatable :: zhu   (:,:)
    real(8), allocatable :: zht   (:,:)
    real(8), allocatable :: xpres (:)
    real(8), allocatable :: zlat  (:)
    real(8), allocatable :: toto3obs(:),PP(:,:) 
    real(8), allocatable :: ozo     (:,:)
    
    real(8) :: zlon
    real(8) :: zptop, zptopmbs
 

    if ( .not. beSilent ) write(*,*) "Entering tvs_fillProfiles subroutine"
  
    if (tvs_nobtov == 0) return    ! exit if there are no tovs data


    !  1.    Set index for model's lowest level and model top
    
    nlv_M = col_getNumLev(columnghr,'MM')
    nlv_T = col_getNumLev(columnghr,'TH')

    if (  col_getPressure(columnghr,1,1,'TH') < col_getPressure(columnghr,nlv_T,1,'TH') ) then
      ilowlvl_M = nlv_M
      ilowlvl_T = nlv_T
    else
      ilowlvl_M = 1
      ilowlvl_T = 1
    end if

    vco => col_getVco(columnghr)
    status = vgd_get(vco%vgrid,key='ig_1 - vertical coord code',value=Vcode)
    diagTtop = (Vcode.eq.5002)
    if ( .not. beSilent ) write(*,*) 'tvs_fillProfiles: diagTtop=', diagTtop 

    ! find model level top, within 0.000001 mbs.
    zptop    = col_getPressure(columnghr,1,1,'TH')
    zptopmbs = zptop / 100.d0
    zptopmbs = zptopmbs - 0.000001d0
    if ( .not. beSilent ) write(*,*) 'tvs_fillProfiles: zptopmbs=',zptopmbs

    TopAt10hPa = ( abs( zptopmbs - 10.0d0 ) <= .1d0 )

    ierr = newdate(datestamp,ijour,itime,-3)
    if (ierr < 0) then
      write(*,*) "Invalid datestamp ",datestamp,ijour,itime,ierr
      call utl_abort('tvs_fillProfiles')
    end if
    year= ijour / 10000
    month = mod(ijour / 100,100)
    day = mod(ijour,100)

    !  1.2   Read ozone climatology

    call ozo_read_climatology(datestamp)

    !     2.  Fill profiles structure
    
    ! loop over all instruments
    sensor_loop: do sensorIndex=1,tvs_nsensors

      ! first loop over all obs.
      profileCount = 0
      bobs1: do tovsIndex = 1, tvs_nobtov
        if (tvs_lsensor(tovsIndex) == sensorIndex) then
          profileCount = profileCount + 1
          NOBMAX = tovsIndex
        end if
      end do bobs1

      if (profileCount == 0) cycle sensor_loop

      nlevels = tvs_coefs(sensorIndex) %coef% nlevels
      allocate ( xpres (nlevels) )
      xpres = tvs_coefs(sensorIndex)% coef % ref_prfl_p
      jpmotop = 1
      do levelIndex = 2, nlevels
        if ( zptopmbs >= xpres(levelIndex - 1) .and. &
             zptopmbs < xpres(levelIndex)        ) then 
          jpmotop = levelIndex
          exit
        end if
      end do
      if ( .not. beSilent ) write(*,*) 'tvs_fillProfiles: jpmotop=', sensorIndex, jpmotop
      jpmolev = nlevels - jpmotop + 1

      allocStatus(:) = 0
      allocate (iptobs    (profileCount)        ,stat= allocStatus(1) )
      allocate (iptobscma (profileCount)        ,stat= allocStatus(2) )
      allocate (zlat      (profileCount)        ,stat= allocStatus(3) )
      allocate (ozo       (nlevels,profileCount),stat= allocStatus(4)) 
      allocate (to        (jpmolev,profileCount),stat= allocStatus(5))
      allocate (huo       (jpmolev,profileCount),stat= allocStatus(6))
      allocate (loghuo    (jpmolev,profileCount),stat= allocStatus(7))	
      allocate (toext     (nlevels  ,profileCount),stat= allocStatus(8))
      allocate (qoext     (nlevels  ,profileCount),stat= allocStatus(9))
      allocate (zvlev     (nlv_T,profileCount),stat= allocStatus(10))
      allocate (zt        (nlv_T,profileCount),stat= allocStatus(11))
      allocate (zhu       (nlv_T,profileCount),stat= allocStatus(12))
      allocate (logzhu    (nlv_T,profileCount),stat= allocStatus(13))
      allocate (zht       (nlv_T,profileCount),stat= allocStatus(14))
      call utl_checkAllocationStatus(allocStatus, " tvs_fillProfiles")
      
      profileCount = 0

      ! second loop over all obs.
      bobs2: do tovsIndex = 1, NOBMAX
        if (tvs_lsensor(tovsIndex) /= sensorIndex) cycle bobs2
        profileCount = profileCount + 1
        iptobs(profileCount) = tovsIndex
        headerIndex = tvs_headerIndex(tovsIndex)
        iptobscma(profileCount) = headerIndex

        !    extract land/sea/sea-ice flag (0=land, 1=sea, 2=sea-ice)
        ksurf = obs_headElem_i(lobsSpaceData,OBS_OFL,headerIndex)
        tvs_profiles(tovsIndex) % skin % surftype = ksurf

        !    extract satellite zenith and azimuth angle, 
        !    sun zenith angle, cloud fraction, latitude and longitude
        tvs_profiles(tovsIndex) % zenangle   = obs_headElem_r(lobsSpaceData,OBS_SZA,headerIndex)

        !pour ne pas faire planter RTTOV dans le cas (rare) ou l'angle zenithal n'est pas defini ou invalide         
        if (tvs_profiles(tovsIndex) % zenangle < 0.0d0 .or. &
             tvs_profiles(tovsIndex) % zenangle > zenmax ) then
          write(*,*) "!!! WARNING !!!"
          write(*,*) "INVALID ZENITH ANGLE"
          write(*,*) "angle, profile number, sensor", tvs_profiles(tovsIndex) % zenangle, tovsIndex, sensorIndex
          write(*,*) "replaced by 0.0 !!!"
          tvs_profiles(tovsIndex) % zenangle = 0.d0
        end if
 
        tvs_profiles(tovsIndex) % azangle  = obs_headElem_r(lobsSpaceData,OBS_AZA,headerIndex)
        tvs_profiles(tovsIndex) % sunazangle  = obs_headElem_r(lobsSpaceData,OBS_SAZ,headerIndex) ! necessaire pour radiation solaire
        iplatform = tvs_coefs(sensorIndex) % coef % id_platform
        instrum = tvs_coefs(sensorIndex) % coef % id_inst
        if ( (instrum == inst_id_amsua .or. instrum == inst_id_mhs) .and. iplatform /= platform_id_eos ) then
          !Correction sur la definition de l'angle. A ammeliorer. Ok pour l'instant.
          tvs_profiles(tovsIndex) % azangle   = tvs_profiles(tovsIndex) % sunazangle + tvs_profiles(tovsIndex) % azangle
          if ( tvs_profiles(tovsIndex) % azangle > 360.d0 ) tvs_profiles(tovsIndex) % azangle = tvs_profiles(tovsIndex) % azangle - 360.d0
        end if
        tvs_profiles(tovsIndex) % sunzenangle = obs_headElem_r(lobsSpaceData,OBS_SUN,headerIndex)
        zlat(profileCount) = obs_headElem_r(lobsSpaceData,OBS_LAT,headerIndex) *MPC_DEGREES_PER_RADIAN_R8
        zlon = obs_headElem_r(lobsSpaceData,OBS_LON,headerIndex) *MPC_DEGREES_PER_RADIAN_R8
        tvs_profiles(tovsIndex) % longitude = zlon
        do levelIndex = 1, nlv_T
          zt   (levelIndex,profileCount) = col_getElem(columnghr,levelIndex,headerIndex,'TT')
          zhu  (levelIndex,profileCount) = col_getElem(columnghr,levelIndex,headerIndex,'HU')
          zvlev(levelIndex,profileCount) = col_getPressure(columnghr,levelIndex,headerIndex,'TH') * MPC_MBAR_PER_PA_R8
          zht  (levelIndex,profileCount) = col_getHeight(columnghr,levelIndex,headerIndex,'TH')
        end do

        if (diagTtop) then
          ! Fix temporaire (?) pour eviter probleme au toit avec GEM 4: on ne veut pas utiliser
          ! le premier niveau de GEM qui est disgnostique (extrapole a partir des deux niveaux plus bas)
          ! (grosse varibilite de la temperature au dernier niveau thermo due a l'extrapolation utilisee)
          zt   (1,profileCount) =  zt   (2,profileCount) + tvs_mesosphereLapseRate *  &
               log( col_getPressure(columnghr,1,headerIndex,'TH') /  &
               col_getPressure(columnghr,2,headerIndex,'TH') )
          zhu  (1,profileCount) =  zhu  (2,profileCount)         ! extrapolation valeur constante pour H2O peu important a cette hauteur
        end if
        
      end do bobs2
 
      !   2.1  Vertical interpolation of model temperature, logarithm of
      !           specific humidity and height levels to pressure levels
      !           required by tovs rt model


      !$omp parallel do private(profileIndex)
      do profileIndex=1, profileCount
        call ppo_IntAvg (zvlev(:,profileIndex:profileIndex),zt(:,profileIndex:profileIndex),nlv_T,1, &
             jpmolev,xpres(jpmotop:nlevels),to(:,profileIndex:profileIndex))
        logzhu(:,profileIndex) = log( zhu(:,profileIndex) )
        call ppo_IntAvg (zvlev(:,profileIndex:profileIndex),logzhu(:,profileIndex:profileIndex),nlv_T,1, &
             jpmolev,xpres(jpmotop:nlevels),loghuo(:,profileIndex:profileIndex))
        huo(:,profileIndex) = exp ( loghuo(:,profileIndex) )
      end do
      !$omp end parallel do


      !    2.2  Extrapolation of temperature profile above model top

      toext(:,:) = 0.0d0
      if ( .not. TopAt10hPa) then ! si le toit n'est pas a 10. hPa 
        do profileIndex=1,profileCount
          toext(jpmotop:nlevels,profileIndex) = to(1:jpmolev,profileIndex)
          ! New approach based on a specified lapse rate
          do levelIndex=1,jpmotop-1
            toext(levelIndex,profileIndex) = toext(jpmotop,profileIndex)  + &
                 ( log(xpres(levelIndex)/xpres(jpmotop)) * tvs_mesosphereLapseRate )
          end do
        end do
      else
        ! old code for temperature profile extrapolation (only apropriate if model top at 10 hPa)
        call extrap (to,toext,jpmolev,nlevels,profileCount)
      end if

      !   2.4  Extrapolation of humidity profile (kg/kg)
      !           above rlimlvhu (normally 300mbs or 70mbs)

      qoext(:,:) = 0.0d0
      
      do profileIndex = 1, profileCount
        do levelIndex = 1, jpmolev
          qoext(nlevels - jpmolev + levelIndex,profileIndex) = huo(levelIndex,profileIndex) !exp(lqo(levelIndex,profileIndex)) 
        end do
      end do

      if ( .not. TopAt10hPa ) then ! if model top not at 10. hPa
        qoext(1:jpmotop,1:profileCount) = MPC_MINIMUM_HU_R8 ! to replace with limlvhu ?
      else                    
        if ( tvs_debug ) then
          do profileIndex = 1, profileCount
            write(*,*)'qoext*1000 avant exthum4    = '
            write(*,'(1x,10f8.4)')(qoext(levelIndex,profileIndex)*1000.d0,levelIndex=1,nlevels)
            write(*,*)' '
          end do
        end if
        call exthum4 (xpres(1:nlevels),qoext,limlvhu)
        if ( tvs_debug ) then
          do profileIndex = 1, profileCount
            write(*,*)'qoext*1000 apres exthum4    = '
            write(*,'(1x,10f8.4)')(qoext(levelIndex,profileIndex)*1000.d0,levelIndex=1,nlevels)
            write(*,*)' '
          end do
        end if
      end if


      !    2.5  Get ozone profiles (ppmv)
      if (tvs_coefs(sensorIndex) %coef % nozone > 0) then
        allocate ( toto3obs(profileCount) )     
        toto3obs(:) = 0.d0
        allocate( PP(nlevels,profileCount) )
        do profileIndex=1,profileCount
          PP(1:nlevels,profileIndex)=xpres(1:nlevels)
        end do
        call ozo_get_profile (ozo,toto3obs,zlat,pp,nlevels,profileCount,datestamp)
        deallocate( PP )
        deallocate ( toto3obs )
      end if


      !    2.6  Fill profiles structure

      do  profileIndex = 1 , profileCount 
        tovsIndex = iptobs(profileIndex)
        headerIndex = iptobscma(profileIndex)
        tvs_profiles(tovsIndex) % gas_units       = gas_unit_specconc ! all gas profiles are supposed to be provided in kg/kg (specific humidity, i.e. mass mixing ratio [kg/kg] over wet air)
        tvs_profiles(tovsIndex) % id              = "" ! profile id, up to 128 characters, to consider for use
        tvs_profiles(tovsIndex) % nlevels         = nlevels
        tvs_profiles(tovsIndex) % nlayers         = nlevels - 1
        tvs_profiles(tovsIndex) % date(1)         = year
        tvs_profiles(tovsIndex) % date(2)         = month
        tvs_profiles(tovsIndex) % date(3)         = day
        tvs_profiles(tovsIndex) % latitude        = zlat(profileIndex)
        tvs_profiles(tovsIndex) % elevation       = 0.001d0 * zht(ilowlvl_T,profileIndex) ! unite km
        tvs_profiles(tovsIndex) % skin % watertype= 1 !utilise pour calcul rayonnement solaire reflechi seulement
        tvs_profiles(tovsIndex) % skin % t        = col_getElem(columnghr,1,headerIndex,'TG')
        tvs_profiles(tovsIndex) % skin % salinity = 35.d0 ! for FASTEM-4 only to revise (practical salinity units)
        tvs_profiles(tovsIndex) % skin % fastem(:)= 0.0d0
        tvs_profiles(tovsIndex) % skin % snow_fraction  = 0.d0 ! Surface coverage snow fraction(0-1), used only by IR emissivity atlas
        tvs_profiles(tovsIndex) % skin % soil_moisture  = 0.d0 ! soil moisure (m**3/m**3) not yet used
        tvs_profiles(tovsIndex) % s2m % t         = col_getElem(columnghr,ilowlvl_T,headerIndex,'TT')
        tvs_profiles(tovsIndex) % s2m % q         = 0.3D6  * qppmv2Mixratio ! a value between 0 and 0.6d6 so that RTTOV will not complain; not used
        tvs_profiles(tovsIndex) % s2m % p         = col_getElem(columnghr,1      ,headerIndex,'P0')*MPC_MBAR_PER_PA_R8
        tvs_profiles(tovsIndex) % s2m % u         = col_getElem(columnghr,ilowlvl_M,headerIndex,'UU')
        tvs_profiles(tovsIndex) % s2m % v         = col_getElem(columnghr,ilowlvl_M,headerIndex,'VV')
        tvs_profiles(tovsIndex) % s2m % o         = 0.0d0 !surface ozone never used
        tvs_profiles(tovsIndex) % s2m % wfetc     = 100000.0d0 ! Wind fetch (in meter for rttov10 ?) used to calculate reflection of solar radiation by sea surface
        tvs_profiles(tovsIndex) % idg             = 0
        tvs_profiles(tovsIndex) % Be              = 0.4d0 ! earth magnetic field strength (gauss) (must be non zero)
        tvs_profiles(tovsIndex) % cosbk           = 0.0d0 ! cosine of the angle between the earth magnetic field and wave propagation direction
        tvs_profiles(tovsIndex) % p(:)            = tvs_coefs(sensorIndex) %coef% ref_prfl_p(:)
        tvs_profiles(tovsIndex) % t(:)            = toext(:,profileIndex)
        if (tvs_coefs(sensorIndex) %coef %nozone > 0) &
             tvs_profiles(tovsIndex) % o3(:) = ozo(:,profileIndex) * o3ppmv2Mixratio ! Climatology output is ppmv (over dry or wet air? not sure but this conversion is only approximate but it should not matter                                                                                                             ! because atmosphere is very dry where there is significant absorption by ozone)
        tvs_profiles(tovsIndex) % q(:)            = qoext(:,profileIndex)
        tvs_profiles(tovsIndex) % ctp = 1013.25d0
        tvs_profiles(tovsIndex) % cfraction = 0.d0
        
      end do

      deallocate (xpres     ,stat= allocStatus(1))
      deallocate (zht       ,stat= allocStatus(2))
      deallocate (logzhu    ,stat= allocStatus(3))
      deallocate (zhu       ,stat= allocStatus(4))
      deallocate (zt        ,stat= allocStatus(5))
      deallocate (zvlev     ,stat= allocStatus(6))
      deallocate (qoext     ,stat= allocStatus(7))
      deallocate (toext     ,stat= allocStatus(8))
      deallocate (loghuo    ,stat= allocStatus(9))
      deallocate (huo       ,stat= allocStatus(10))
      deallocate (to        ,stat= allocStatus(11))
      deallocate (ozo       ,stat= allocStatus(12))
      deallocate (zlat      ,stat= allocStatus(13))
      deallocate (iptobscma ,stat= allocStatus(14))
      deallocate (iptobs    ,stat= allocStatus(15))
    
      call utl_checkAllocationStatus(allocStatus, " tvs_fillProfiles", .false.)
     
    end do sensor_loop

  end subroutine tvs_fillProfiles

  !--------------------------------------------------------------------------
  !  exthum4
  !--------------------------------------------------------------------------
  subroutine exthum4( ppres, humidity, limlvhu )
    ! :Purpose: extrapolate upper level humidity profile
    !           (adapted from exthum by J. Eyre).
    !           To extend mixing ratio profile into stratosphere in a reasonable way.
    !
    !
    ! :Method:
    !     take top tropospheric mixing ratio (e.g. near 300 mb) and
    !     extrapolate with given fall off into lower stratosphere
    !     (e.g. to 70 mb).  constrain mixing ratio to be >= zwmin
    !     (e.g. 0.000003 kg/kg). In upper strat, mixing ratio = zwmin.
    !
    ! :Externals:
    !          none.
    !
    ! :Reference:
    !     ecmwf tech mem 176.
    !
    implicit none

    ! Arguments
    real(8), intent(in)    :: ppres(:)      ! Pressure levels
    real(8), intent(inout) :: humidity(:,:) ! Humidity profiles
    real(8), intent(in)    :: limlvhu       ! Humidity is extrapolated for pressures < limlvhu

    ! locals
    real(8),allocatable :: zpres3(:)
    real(8) zwb
    real(8),parameter :: zp1 = 70.0d0  !press limits (in hpa) of region to be extrapolated
    integer :: topIndex, profileIndex,levelIndex,nlevels,nprofiles

    nlevels = size( ppres )
    nprofiles = size( humidity, dim =2)
    
    !      Find top level of given profile
    topIndex = 0
    do levelIndex=nlevels,1,-1
      if (ppres(levelIndex) < limlvhu) then
        topIndex = levelIndex
        exit
      end if
    end do

    !  Null extrapolation case
    if (topIndex == 0) return

    allocate ( zpres3( nlevels ) )
      !   Constants defining p**3 fall off around tropopause
    do levelIndex=1,topIndex
      zpres3(levelIndex)=(ppres(levelIndex)/ppres(topIndex+1))**3
    end do

    do profileIndex=1,nlevels
      zwb=humidity(topIndex+1,profileIndex)
      do levelIndex=1,topIndex
        if (ppres(levelIndex)<zp1) then
          humidity(levelIndex,profileIndex)=MPC_MINIMUM_HU_R8
        else
          humidity(levelIndex,profileIndex)=max((zwb*zpres3(levelIndex)),MPC_MINIMUM_HU_R8)
        end if
      end do
    end do

    deallocate ( zpres3 )
       
  end subroutine exthum4

  !--------------------------------------------------------------------------
  !  htextrap
  !--------------------------------------------------------------------------
  subroutine htextrap( profout, profin, xpres, jplev, jpmolev, jpmotop, nprf )
    !
    ! :Purpose: extrapolate height profiles above 10mb model top
    !           on rttov levels up to 0.1mb (rttov levels 1 to 7)
    !           using 10 rttov height levels from 100mb to 10mb
    !           (rttov levels 8 to 17) for linear fit.
    !
    !                #. linear extrapolation following
    !                #. profout(m) = a * ln(xpres(mb)) + b
    !                #. and solve a and b by least square method 
    !
    implicit none

    !  Arguments:
    real(8), intent(out) ::  profout(jplev,nprf) ! Height profiles  -extrapolated- (m)
    real(8), intent(in)  ::  profin(jpmolev,nprf)! Height profiles  -to be extrapolated- (m)
    real(8), intent(in)  ::  xpres(jplev)        ! Pressure levels of rttov model (hpa)
    integer, intent(in)  ::  jplev               ! Number of pressure levels of rttov model
    integer, intent(in)  ::  jpmolev             ! Number of rttov model levels below nwp model top
    integer, intent(in)  :: jpmotop              ! First rttov model level under nwp model top
    integer, intent(in)  ::  nprf                ! Number of profiles

    ! Locals:
    integer     :: i, jk, jn
    real(8)      :: lnx_sum, lnx_avg, y_sum, y_avg, a_num, a_den, a, b
    integer, parameter :: nl = 10  ! number of points used in the extrapolation

    do jn = 1, nprf
      
      lnx_sum = 0.d0
      y_sum   = 0.d0
      a_num   = 0.d0
      a_den   = 0.d0


      ! Find averaged values of height and ln ( pressure )

      do i = 1, nl
        lnx_sum = lnx_sum + log(xpres(jpmotop + i - 1))
        y_sum   = y_sum   + profin(i,jn)
      end do

      lnx_avg = lnx_sum / nl
      y_avg   = y_sum   / nl

      !  Find constants a and b by least-square method
      
      do i = 1, nl
        a_num = a_num + ( log(xpres(jpmotop + i - 1)) - lnx_avg ) * ( profin(i,jn) - y_avg )
        a_den = a_den + ( log(xpres(jpmotop + i - 1)) - lnx_avg )**2
      end do

      a = a_num / a_den

      b = y_avg - A * lnx_avg


      !   Initialize height for rttov levels under nwp model top

      do jk = 1, jpmolev
        profout(jplev - jpmolev + jk,jn) = profin(jk,jn)
      end do


      !   Extrapolate height for rttov levels above nwp model top

      do jk = 1, jpmotop - 1
        profout(jk,jn) = a * log(xpres(jk)) + b
      end do

    end do

  end subroutine htextrap

  !--------------------------------------------------------------------------
  !  tvs_rttov
  !--------------------------------------------------------------------------
  subroutine tvs_rttov( lobsSpaceData, bgckMode, beSilent )
    !
    ! :Purpose: Interface for RTTOV non linear operator
    !           tvs_fillProfiles should be called before
    !
    implicit none

    ! Arguments
    type(struct_obs), intent(inout) :: lobsSpaceData  ! obsSpaceData structure
    logical, intent(in)             :: bgckMode       ! flag to transfer transmittances and cloudy overcast radiances in bgck mode 
    logical, intent(in)             :: beSilent       ! flag to control verbosity

    ! Locals:
    integer :: nlevels
    integer :: btCount
    integer :: allocStatus(3)
    integer :: rttov_err_stat ! rttov error return code
    integer, external :: omp_get_num_threads
    integer :: nthreads,max_nthreads
    integer :: sensor_id, tovsIndex
    integer :: channelIndex
    integer :: profileCount
    integer :: profileIndex, levelIndex, jj, btIndex
    integer :: instrum
    integer :: sensor_type        ! sensor type (1=infrared; 2=microwave; 3=high resolution,4=polarimetric)
    integer, allocatable:: iptobs  (:)  
    
    type (rttov_emissivity), pointer :: emissivity_local (:)    ! emissivity structure with input and output
    type (rttov_chanprof), pointer :: chanprof(:)
    type (rttov_chanprof), allocatable :: chanprof1(:)
    type (rttov_radiance) :: radiancedata_d, radiancedata_d1
    type (rttov_transmission) :: transmission, transmission1
    type (rttov_emis_atlas_data), allocatable, save :: Atlas(:)
    logical, save        :: first=.true.
    integer              :: asw
    logical, pointer :: calcemis  (:)
    real*8, allocatable  :: surfem1 (:)
    real*8, allocatable  :: surfem2 (:)
    integer              :: profileIndex2, tb1, tb2

    if ( .not. beSilent ) write(*,*) "Entering tvs_rttov subroutine"
    if ( .not. beSilent ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if (tvs_nobtov == 0) return                  ! exit if there are not tovs data

    !   1.  Get number of threads available and allocate memory for some variables
    !$omp parallel
    max_nthreads = omp_get_num_threads()
    !$omp end parallel
    allocStatus(:) = 0

    allocate(iptobs(tvs_nobtov),stat=allocStatus(1))
    call utl_checkAllocationStatus(allocStatus(1:1), " tvs_rttov iptobs")
    
    
    !   1.1   Read surface information
    if ( bgckMode ) call EMIS_READ_CLIMATOLOGY

    !   2.  Computation of hx for tovs data only

    ! Loop over all sensors specified by user
    sensor_loop:do sensor_id = 1, tvs_nsensors
   
      nlevels = tvs_coefs(sensor_id)% coef % nlevels
      sensor_type = tvs_coefs(sensor_id) % coef % id_sensor
      instrum = tvs_coefs(sensor_id) % coef % id_inst
    
      !  loop over all obs.
      profileCount = 0
      obs_loop: do tovsIndex = 1, tvs_nobtov
        
        !    Currently processed sensor?
        if ( tvs_lsensor(tovsIndex) == sensor_id ) then
          profileCount = profileCount + 1
          iptobs(profileCount) = tovsIndex
        end if
      end do obs_loop
      
      if (profileCount == 0) cycle sensor_loop
      
      !    2.1  Calculate the actual number of threads which will be used.
      
      nthreads = min(max_nthreads, profileCount )  
      
      !    2.2  Prepare all input variables required by rttov.
      
      if ( bgckMode .and. tvs_isInstrumHyperSpectral(instrum) ) then
        btCount = profileCount * tvs_nchan(sensor_id)
      else
        btCount = tvs_countRadiances(iptobs(1:profileCount), lobsSpaceData)
      end if
      
      if ( btCount == 0 ) cycle sensor_loop
      
      allocate ( surfem1          (btCount) ,stat=allocStatus(1))

      asw = 1 ! Allocation
      call rttov_alloc_direct(            &
              allocStatus(2),             &
              asw,                        &
              nprofiles=profileCount,     & ! (not used)
              nchanprof=btCount,          &
              nlevels=nlevels,            &
              chanprof=chanprof,          &
              opts=tvs_opts(sensor_id),   &
              coefs=tvs_coefs(sensor_id), &
              transmission=transmission,  &
              radiance=radiancedata_d,    &
              calcemis=calcemis,          &
              emissivity=emissivity_local,&
              init=.true.)


      if (useUofWIREmiss) then
        allocate ( surfem2(btCount)        ,stat=allocStatus(3) )
      end if
      call utl_checkAllocationStatus(allocStatus, " tvs_rttov")
      
      
      !     get Hyperspectral IR emissivities
      if ( tvs_isInstrumHyperSpectral(instrum) ) then
        surfem1(:) = 0.
        if ( bgckMode ) then
          call emis_getIrEmissivity (surfem1,tvs_nchan(sensor_id),sensor_id,profileCount,btCount,iptobs)
        else
          call tvs_getHIREmissivities(iptobs(1:profileCount), lobsSpaceData, surfem1)
        end if
      end if
      
      if ( bgckMode .and. tvs_isInstrumHyperSpectral(instrum) ) then
        btIndex = 0
        do profileIndex = 1 , profileCount
          do  channelIndex = 1,tvs_nchan(sensor_id)
            btIndex = btIndex + 1
            chanprof(btIndex)%prof = profileIndex
            chanprof(btIndex)%chan = channelIndex
          end do
        end do
      else
        call tvs_getChanprof(sensor_id, iptobs(1:profileCount), lobsSpaceData, chanprof)
      end if
                 
      call tvs_getOtherEmissivities(chanprof, iptobs, sensor_type, instrum, surfem1, calcemis)
      
      if (useUofWIREmiss .and. tvs_isInstrumHyperSpectral(instrum) .and. bgckMode) then
        if (first) then
          if (.not. allocated (Atlas)) allocate(Atlas(tvs_nsensors))
          call rttov_setup_emis_atlas( rttov_err_stat,       &! out
               tvs_opts(sensor_id),                    &! in
               tvs_profiles(1)%date(2) ,               &! in
               atlas_type_ir,                          &! in
               atlas(sensor_id),                       &! in
               ir_atlas_ang_corr = .false.,            &! in
               ir_atlas_read_std = .false.,            &! in
               coefs = tvs_coefs(sensor_id)  )
          if (rttov_err_stat/=0) then
            write(*,*) "Error in rttov_atlas_setup ",rttov_err_stat
            call utl_abort('tvs_rttov')
          end if
          first=.false.
        end if
        
        call rttov_get_emis( rttov_err_stat        , &   ! out
             tvs_opts(sensor_id)                   , &   ! in
             chanprof(1:btCount)                  , &   ! in
             tvs_profiles(iptobs(1:profileCount)) , &   ! in
             tvs_coefs(sensor_id)                  , &   ! in
             Atlas(sensor_id)                      , &   ! in
             surfem2(1:btCount)     ) ! out

        if (rttov_err_stat /= 0) then
          write(*,*) "Error in rttov_get_emis ", rttov_err_stat
          call utl_abort('tvs_rttov')
        end if
              
        do profileIndex=1, profileCount !loop on profiles
          jj = iptobs(profileIndex)
          do btIndex=1, btCount !loop on channels
            if (chanprof(btIndex)%prof==profileIndex) then
              ! surftype: 0 land, 1 sea, 2 sea-ice
              ! this logic is primitive and could be improved for example using
              ! additional criteria based on emissivity_std and emissivity_flg
              !Definition of emis_flag:
              ! emis_flag:Flag_0 = "0 = sea, no MOD11 data" ;
              ! emis_flag:Flag_1 = "1 = land where BF method was applied" ;
              ! emis_flag:Flag_2 = "2 = land where data was filled with average (original UWiremis bfemis_flag=2 or 3 or 4" ;
              ! emis_flag:Flag_3 = "3 = contains inland water or coastline by the sea/land mask where the BF method was used" ;
              ! emis_flag:Flag_4 = "4 = contains inland water or coastline by the sea/land mask where data was filled with average original UWiremis bfemis_flag=2 or 3 or 4" ;
              ! emis_flag:Flag_5 = "5 = contains coastline by land fraction where the BF method was used" ;
              ! emis_flag:Flag_6 = "6 = contains coastline by land fraction where data was filled with average (original UWiremis bfemis_flag=2 or 3 or 4" ;
              ! other information that could be useful for quality control can be found in the in the profile_qc structure
              ! Now we have the "traditionnal" emissivity in surfem1(:)
              ! and University of Wisconsin emissivity in surfem2(:)
              if (tvs_profiles(jj)% skin % surftype == 0 .and. &
                   surfem2(btIndex) > 0.5 ) then
                emissivity_local(btIndex)%emis_in = surfem2(btIndex)
              else
                emissivity_local(btIndex)%emis_in = surfem1(btIndex)
              end if
            end if
          end do
        end do
      else
        emissivity_local(:)%emis_in = surfem1(:)
      end if
        
      !   2.3  Compute radiance with rttov_direct

      rttov_err_stat = 0 

      if( bgckMode .and. tvs_isInstrumHyperSpectral(instrum) ) then
        write(*,*) 'for bgck IR: call rttov_parallel_direct for each profile...'

        asw = 1 ! 1 to allocate,0 to deallocate
        ! allocate transmitance structure for 1 profile
        call rttov_alloc_transmission(allocStatus(1), transmission1, nlevels=nlevels, &
             nchanprof=tvs_nchan(sensor_id), asw=asw, init=.true. )
        ! allocate radiance structure for 1 profile
        call rttov_alloc_rad (allocStatus(1),tvs_nchan(sensor_id), radiancedata_d1,nlevels,asw,init=.true.)
        ! allocate chanprof for 1 profile
        allocate(chanprof1(tvs_nchan(sensor_id)))
        do  channelIndex = 1,tvs_nchan(sensor_id)
          chanprof1(channelIndex)%prof = 1
          chanprof1(channelIndex)%chan = channelIndex
        end do

        do profileIndex2 = 1, profileCount
          tb1 = 1 + (profileIndex2-1) * tvs_nchan(sensor_id) 
          tb2 = profileIndex2 * tvs_nchan(sensor_id)
          call rttov_parallel_direct(       &
               rttov_err_stat,              & ! out
               chanprof1,                   & ! in
               tvs_opts(sensor_id),         & ! in
               tvs_profiles(iptobs(profileIndex2):iptobs(profileIndex2)),  & ! in
               tvs_coefs(sensor_id),                 & ! in
               transmission1,                        & ! inout
               radiancedata_d1,                      & ! inout
               calcemis=calcemis(tb1:tb2),           & ! in
               emissivity=emissivity_local(tb1:tb2), & ! inout
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
        call rttov_alloc_transmission(allocStatus(1),transmission1,nlevels=nlevels,  &
             nchanprof=tvs_nchan(sensor_id), asw=asw )
        ! radiance deallocation for 1 profile
        call rttov_alloc_rad (allocStatus(1), tvs_nchan(sensor_id), radiancedata_d1, nlevels, asw)

      else

        call rttov_parallel_direct(       &
             rttov_err_stat,              & ! out
             chanprof,                    & ! in
             tvs_opts(sensor_id),         & ! in
             tvs_profiles(iptobs(1:profileCount)),  & ! in
             tvs_coefs(sensor_id),        & ! in
             transmission,                & ! inout
             radiancedata_d,              & ! inout
             calcemis=calcemis,           & ! in
             emissivity=emissivity_local, & ! inout
             nthreads=nthreads      )   

      end if

      if ( .not. beSilent ) write(*,*) 'after rttov_parallel_direct...'
      if ( .not. beSilent ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'     

      if (rttov_err_stat /= 0) then
        write(*,*) "Error in rttov_parallel_direct",rttov_err_stat
        call utl_abort('tvs_rttov')
      end if
                                        
      !    2.4  Store hx in the structure tvs_radiance

      do btIndex = 1, btCount
        profileIndex = chanprof(btIndex)%prof
        channelIndex = chanprof(btIndex)%chan
        tovsIndex = iptobs(profileIndex)
        tvs_radiance(tovsIndex) % bt(channelIndex) =     &
             radiancedata_d % bt(btIndex)
        if ( bgckMode ) then
          tvs_radiance(tovsIndex) % clear(channelIndex) =  &
               radiancedata_d %clear(btIndex)
          do levelIndex = 1, nlevels - 1
            tvs_radiance(tovsIndex) % overcast(levelIndex,channelIndex) =   &
                 radiancedata_d % overcast(levelIndex,btIndex)
          end do
        end if

        if ( allocated( tvs_transmission) ) then
          do levelIndex = 1, nlevels
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

      !    Deallocate memory
      asw = 0 ! 0 to deallocate
      call rttov_alloc_direct(         &
           allocStatus(1),             &
           asw,                        &
           nprofiles=profileCount,     & ! (not used)
           nchanprof=btCount,          &
           nlevels=nlevels,            &
           chanprof=chanprof,          &
           opts=tvs_opts(sensor_id),   &
           coefs=tvs_coefs(sensor_id), &
           transmission=transmission,  &
           radiance=radiancedata_d,    &
           calcemis=calcemis,          &
           emissivity=emissivity_local,&
           init=.true.)

 
      if (useUofWIREmiss) then
        deallocate ( surfem2  ,stat=allocStatus(2) )
      end if
      deallocate ( surfem1    ,stat=allocStatus(3) )
      call utl_checkAllocationStatus(allocStatus, " tvs_rttov", .false.)
      
    end do sensor_loop
    
    deallocate(iptobs)

  end subroutine tvs_rttov

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
    integer ,parameter :: MaxWn = 19
    integer ,parameter :: Nparm=3
    integer ,parameter :: MaxChan=19

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
    real (8),intent(out) :: f_low(nprf)       ! Low resolution field
    real (8),intent(in)  :: f_high(klon, klat)! High resolution field 
    integer, intent(in)  :: nprf              ! Number of profiles
    integer, intent(in)  :: ilat(nprf)        ! Y-coordinate of profile
    integer, intent(in)  :: ilon(nprf)        ! X-coordinate of profile
    integer, intent(in)  :: klon              ! Max value of latitude indices
    integer, intent(in)  :: klat              ! Max value of longitude indices
    integer, intent(in)  :: ireduc            ! Means a 2xireduc+1 by 2xireduc+1 averaging

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
  subroutine emis_getIrEmissivity (surfem1, nchn, sensorIndex, nprf, nchannels_max, iptobs)
    !
    ! :Purpose: Assign new ir surface emissivities based on
    !           cmc analysis surface albedo, sea ice fraction and snow mask
    !           in addition to ceres surface type and water fraction. 
    !           This is a subroutine that can apply to any instrument.
    !           However, due to the necessity of specifying the instrument
    !           bands wavenumbers, the use of this subroutine for a new instrument
    !           would require the minor following changes.
    !
    implicit none
   
    ! Arguments:
    real(8),intent(out) :: surfem1(nchannels_max) ! IR surface emissivity estimate (0-1)
    integer,intent(in)  :: nchn                   ! Number of channels
    integer,intent(in)  :: sensorindex            ! Sensor number
    integer,intent(in)  :: nprf                   ! Number of profiles
    integer,intent(in)  :: nchannels_max          ! Total number of observations treated
    integer,intent(in)  :: iptobs( nprf )         ! Profile position number

    !Locals:
    integer :: jc,jn
    integer :: ilat(nprf), ilon(nprf)
    real(8) :: zlat(nprf), zlon(nprf), satzang(nprf)
    real (8) :: wind_sfc(nprf), f_low(nprf), waven(nchn), em_oc(nchn,nprf), emi_mat(nchn,20)


    ! Information to extract (transvidage)
    ! zlat(nprf) -- latitude (-90 to 90)
    ! zlon(nprf) -- longitude (0 to 360)
    ! satzang(nprf) -- satellite zenith angle (deg)

    do jn = 1, nprf
      zlat(jn)    = tvs_profiles(iptobs(jn))% latitude
      zlon(jn)    = tvs_profiles(iptobs(jn))% longitude
      satzang(jn) = tvs_profiles(iptobs(jn))% zenangle
    end do

    !  Assign surface properties from grid to profiles
    call interp_sfc(ilat,ilon, nprf,zlat,zlon,iptobs)


    !  Find the sensor bands (central) wavenumbers
    do jc = 1, nchn      
      waven(jc) = tvs_coefs(sensorIndex) % coef % ff_cwn(jc)
    end do


    !  Get the CERES emissivity matrix for all sensor wavenumbers and surface types
    call ceres_ematrix(emi_mat, waven,nchn)


    ! Refine water emissivities

    do jn = 1, nprf
      !       find surface wind
      wind_sfc(jn) = min(sqrt(tvs_profiles(iptobs(jn))%S2M%U**2 + tvs_profiles(iptobs(jn))%S2M%V**2 + 1.d-12),15.d0)
    end do

    !     find new ocean emissivities     

    do jc = 1, nchn
      em_oc(jc,:)= emi_mat(jc,17)
    end do
    
    call emi_sea (em_oc, waven,satzang,wind_sfc,nprf,nchn)
    

    ! Get surface emissivities

    do jn = 1, nprf
      !       set albedo to 0.6 where snow is present
      if ( tvs_profiles(iptobs(jn))%SKIN%SURFTYPE == 0 .and. tvs_surfaceParameters(iptobs(jn))%snow > 0.999 ) tvs_surfaceParameters(iptobs(jn))%albedo = 0.6
      !       if albedo too high no water
      if ( tvs_surfaceParameters(iptobs(jn))%albedo >= 0.55 ) tvs_surfaceParameters(iptobs(jn))%pcnt_wat = 0.
      !       if water and CMC ice present then sea ice
      if ( tvs_profiles(iptobs(jn))%SKIN%SURFTYPE == 1 .and. tvs_surfaceParameters(iptobs(jn))%ice > 0.001 ) tvs_surfaceParameters(iptobs(jn))%ltype = 20
      !       if land and CMC snow present then snow
      if ( tvs_profiles(iptobs(jn))%SKIN%SURFTYPE == 0 .and. tvs_surfaceParameters(iptobs(jn))%snow > 0.999 ) tvs_surfaceParameters(iptobs(jn))%ltype = 15
      do jc=1,nchn
        surfem1((jn-1)*nchn+jc) =  tvs_surfaceParameters(iptobs(jn))%pcnt_wat * em_oc(jc,jn)  +   &
             ( 1.d0 - tvs_surfaceParameters(iptobs(jn))%pcnt_wat ) * emi_mat(jc,tvs_surfaceParameters(iptobs(jn))%ltype)
      end do
    end do

    ! Find the regional water fraction (here in a 15x15 pixel box centered on profile)
    call pcnt_box (f_low, waterFraction,nprf,ilat,ilon,kslat,kslon,7)

    do jn = 1, nprf
      tvs_surfaceParameters(iptobs(jn))%pcnt_reg = f_low(jn)
    end do

  end subroutine emis_getIrEmissivity

  !--------------------------------------------------------------------------
  !  interp_sfc
  !--------------------------------------------------------------------------
  subroutine interp_sfc (ilat, ilon, nprf, zlat, zlon, iptobs)
    !
    ! :Purpose: Associate surface albedo, ice fraction, snow depth 
    !           and ceres surface type and water fraction to observations profiles.

    implicit none

    ! Arguments
    integer,intent(out) :: ilat(nprf)   ! y-coordinate of profile
    integer,intent(out) :: ilon(nprf)   ! x-coordinate of profile 
    integer,intent(in)  :: nprf         ! number of profiles
    real(8),intent(in)  :: zlat(nprf)   ! latitude (-90s to 90n)
    real(8),intent(in)  :: zlon(nprf)   ! longitude (0 to 360)
    integer,intent(in)  :: iptobs(nprf) ! observation index

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
    integer,external   :: fnom,fstouv,fstinf,fstprm,fstfrm,fclos
    integer,external   :: ezqkdef,ezdefset
    real(8)            :: zig1,zig2,zig3,zig4
    integer            :: ig1obs,ig2obs,ig3obs,ig4obs
    real (8)           :: alat, alon, zzlat, zzlon
    ! fields on input grid
    real(8),allocatable:: glace(:,:), neige(:,:), alb(:,:)
    ! fields on output grid
    real(8)            :: glace_intrpl(nprf), neige_intrpl(nprf), alb_intrpl(nprf)


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
      zzlat = min(zlat(jn),89.999d0)
      zzlat = max(Zzlat,-89.999d0)
      
      zzlon = min(zlon(jn),359.999d0)
      zzlon = max(zzlon,0.d0)

      !  Find in which surface field pixel is located the observation profile

      ! Note : CERES grid at 1/6 resolution 
      !         N-S : starts at N pole and excludes S pole
      !         W-E : starts at longitude 0 and excludes longitude 360

      ilat(jn) = max( nint((zzlat + 90.d0) * alat),1) 
      ilon(jn) = nint(zzlon * alon) + 1
      if (ilon(jn) > kslon) ilon(jn) = 1

      !  Assign surface caracteristics to observation profiles

      tvs_surfaceParameters(iptobs(jn)) % ltype    = surfaceType(ilon(jn),ilat(jn))
      tvs_surfaceParameters(iptobs(jn)) % pcnt_wat = waterFraction(ilon(jn),ilat(jn))

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

    ! utl_CXGAIG: define the grid descriptors (integer form) of the
    !          observation profile output grid
    ! desired output = ig1OBS, ig2OBS, ig3OBS, ig4OBS
    zig1 = 0.0d0
    zig2 = 0.0d0
    zig3 = 1.0d0
    zig4 = 1.0d0

    call utl_cxgaig('L',ig1OBS,ig2OBS,ig3OBS,ig4OBS,zig1,zig2,zig3,zig4)


    ! utl_EZGDEF: define the grid of the observations profiles (output grid)
    ! of type Y containing the lat-lon of profiles
    ! success = token to identify the grid
    ! desired output = token
    write(*,*) 
    iv7 = utl_ezgdef(nprf,1,'Y','L',ig1obs,ig2obs,ig3obs,ig4obs,zlon,zlat)
    write(*,*) 'apply to all variables  : utl_EZGDEF : return = ', iv7
    

    ! EZQKDEF: define the grid of the records data (input grid)
    ! success = token to identify the grid
    ! desired output = token
    ! EZDEFSET: interpolate from input grids to output grid
    ! success = key
    ! utl_EZSINT: interpolation of the field on the input grid to observation profiles
    ! success = 0
    ! desired output = FIELD_intrpl
    write(*,*) 
    ix8 = ezqkdef(ni3,nj3,grtyp3,ig13,ig23,ig33,ig43,iun3)
    write(*,*) 'variable = LG           : ezqkdef  : return = ', ix8
    
    ix9 = ezdefset(iv7,ix8)
    write(*,*) 'variable = LG           : ezdefset : return = ', ix9

    ix10 = utl_ezsint(glace_intrpl,glace,interpDegree='NEAREST')
    write(*,*) 'variable = LG           : utl_ezsint  : return = ', ix10

    write(*,*) 

    iy8 = ezqkdef(ni4,nj4,grtyp4,ig14,ig24,ig34,ig44,iun3)
    write(*,*) 'variable = ', snowvar, '           : ezqkdef  : return = ', iy8

    iy9 = ezdefset(iv7,iy8)
    write(*,*) 'variable = ', snowvar, '           : ezdefset : return = ', iy9

    iy10 = utl_ezsint(neige_intrpl,neige,interpDegree='NEAREST')
    write(*,*) 'variable = ', snowvar, '           : utl_ezsint  : return = ', iy10

    write(*,*) 

    iz8 = ezqkdef(ni5,nj5,grtyp5,ig15,ig25,ig35,ig45,iun5)
    write(*,*) 'variable = AL           : ezqkdef  : return = ', iz8

    iz9 = ezdefset(iv7,iz8)
    write(*,*) 'variable = AL           : ezdefset : return = ', iz9

    iz10 = utl_ezsint(alb_intrpl,alb,interpDegree='NEAREST')
    write(*,*) 'variable = AL           : utl_ezsint  : return = ', iz10


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
      tvs_surfaceParameters(iptobs(jn))%ice      = glace_intrpl(jn)
      tvs_surfaceParameters(iptobs(jn))%snow     = neige_intrpl(jn)
      tvs_surfaceParameters(iptobs(jn))%albedo   = alb_intrpl(jn)
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
    real (8),intent(in) :: waven(nchn)       ! wavenumbers (cm-1)
    real (8),intent(out):: emi_mat(nchn, 20) ! emissivity (0.0-1.0)

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
    real (8),intent(out)  :: em_oc(nc,np) ! Ocean emissivities (0.-1.)
    real (8),intent(in)   :: wnum(nc)     ! Channel wavenumbers (cm-1)
    real (8),intent(in)   :: angle(np)    ! Viewing angle (deg)
    real (8),intent(in)   :: wind(np)     ! Surface wind speed (m/s)
    integer,intent(in)    :: np           ! Number of profiles
    integer,intent(in)    :: nc           ! Number of channels

    ! Locals
    integer      :: i, k, l
    integer      :: imem(nc) 
    integer      :: mchan(2)
    real (8)     :: dum
    real (8)     :: emi2(2,np)

    ! Masuda's 19 wavelengths converted to wavenumber
    real (8), parameter :: refw(19)=(/ 2857.1d0, 2777.7d0, 2702.7d0, 2631.6d0, 2564.1d0, &
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
  !  tvs_rttov_read_coefs
  !--------------------------------------------------------------------------
  subroutine tvs_rttov_read_coefs(err, coefs, opts, channels, instrument)
    !
    !  :Purpose: MPI wrapper for rttov_read_coefs
    !            the coefficient files are read by MPI task 0
    !            and then broadcasted to the other tasks according to the selected
    !            channels. Argument channels is mandatory (it is optional in rttov_setup)
    !            optional argument channels_rec was removed (it is useful only in principal component mode)
    !            other optionnal arguments were removed :
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
    !            if necessary these arguments could be  added (ask S. Heilliette)
    !            also this subroutine will work only for clear sky radiance computations
    !            if somebody wants to do realistic cloud or aerosol affected radiance simulations
    !            some changes are needed. Ask me in that case. (S. Heilliette) 
    !            It is implicitely assumed that the options are the same for all MPI tasks for a given instrument
    !            No check will be done (options for task 0 will be used for all tasks). 
    !            Only differences in channel lists are accounted for.
    !

    implicit none

    ! Arguments:
    integer(kind=jpim), intent(out) :: err          ! Error status
    type(rttov_coefs),  intent(out) :: coefs        ! Rttov coefficient structure
    type(rttov_options),intent(in)  :: opts         ! Rttov option structure
    integer(kind=jpim), intent(in)  :: channels(:)  ! Channel list
    integer(kind=jpim), intent(in)  :: instrument(3)! Instrument vector

    ! Locals:
    integer, allocatable :: listGlobal(:)
    real(8), allocatable :: bigArray(:,:,:,:)
    integer :: i, j, k, l, ichan,igas,ierr, countUniqueChannel, indexchan(size(channels)),channelsb(tvs_maxChannelNumber), listAll(tvs_maxChannelNumber)
    logical :: found, associated0
    integer :: nlte_count, nlte_start,isol,isat,nlte_file_nchan
    integer, allocatable :: nlte_chans(:) 

    write(*,*) "Entering tvs_rttov_read_coefs"

    if (size(channels) > tvs_maxChannelNumber) then
      write(*,*) 'You need to increase tvs_maxChannelNumber in tovs_nl_mod !',size(channels), tvs_maxChannelNumber
      call utl_abort("tvs_rttov_setup")
    end if

    ! First step: we should determine a common set of channels among MPI tasks
    if (mpi_myid ==0) then
      allocate(listGlobal(mpi_nprocs*tvs_maxChannelNumber))
    else
      allocate(listGlobal(1))
    end if

    listAll(:) = 0
    listGlobal(:) = 0
    channelsb(:) = 0
    channelsb(1:size(channels)) = channels(:)

    call rpn_comm_barrier("GRID",ierr)

    call rpn_comm_gather(channelsb, tvs_maxChannelNumber, 'MPI_INTEGER', listGlobal, tvs_maxChannelNumber, 'MPI_INTEGER', &
         0, 'GRID', ierr) 
    countUniqueChannel = 0
    if ( mpi_myid == 0 ) then
      call isort(listGlobal, mpi_nprocs*tvs_maxChannelNumber)
      do i=1, mpi_nprocs * tvs_maxChannelNumber
        if (listGlobal(i) >0) then
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

    ! Second step: mpi task 0 will do the job
    if ( mpi_myid == 0 ) then
      call rttov_read_coefs ( &
           ERR,             &! out
           coefs,           &
           opts,            &
           instrument=instrument,      &! in
           channels=listAll(1:countUniqueChannel)  )     ! in option
    else
      call rttov_nullify_coef(coefs%coef)
    end if

    ! Third step: common (i.e. independent from the channel list) parameters are simply broadcasted to other processors
    ! Scalar and fixed size arrays and  strings first

    call rpn_comm_bcast(coefs%initialised, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)       ! Logical flag for initialization
    call rpn_comm_bcast(coefs%coef%id_platform, 1, 'MPI_INTEGER', 0, 'GRID', ierr) !ok
    call rpn_comm_bcast(coefs%coef%id_sat, 1, 'MPI_INTEGER', 0, 'GRID', ierr)  !ok
    call rpn_comm_bcast(coefs%coef%id_inst, 1, 'MPI_INTEGER', 0, 'GRID', ierr)  !ok
    call rpn_comm_bcast(coefs%coef%id_sensor, 1, 'MPI_INTEGER', 0, 'GRID', ierr) !ok
    call rpn_comm_bcast(coefs%coef%id_comp_lvl, 1, 'MPI_INTEGER', 0, 'GRID', ierr)!ok
    call rpn_comm_bcast(coefs%coef%id_comp_pc, 1, 'MPI_INTEGER', 0, 'GRID', ierr)!ok
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
    call rpn_comm_bcast(coefs%coef%nlevels, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of levels(pres/absorber) same for all gases
    call rpn_comm_bcast(coefs%coef%nlayers, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of layers(pres/absorber) nlevels-1
    call rpn_comm_bcast(coefs%coef%pmc_nlay, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%pmc_nvar, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%IncZeeman, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)       ! Flag to include Zeeman effect for this sensor
    call rpn_comm_bcast(coefs%coef%solarcoef, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)       ! Flag to include solar reflection
    call rpn_comm_bcast(coefs%coef%nltecoef, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)        ! Flag to include nlte corrections
    call rpn_comm_bcast(coefs%coef%pmc_shift, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%ff_val_bc, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    call rpn_comm_bcast(coefs%coef%ff_val_gam, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)

    call rpn_comm_bcast(coefs%coef%ncmixed, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of coefficients for Mixed Gases
    call rpn_comm_bcast(coefs%coef%ncwater, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of coefficients for Water Vapour
    call rpn_comm_bcast(coefs%coef%ncozone, 1, 'MPI_INTEGER', 0, 'GRID', ierr)         ! number of coefficients for Ozone
    call rpn_comm_bcast(coefs%coef%ncwvcont, 1, 'MPI_INTEGER', 0, 'GRID', ierr)        ! number of coefficients for WV continuum
    call rpn_comm_bcast(coefs%coef%ncco2, 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for CO2
    call rpn_comm_bcast(coefs%coef%ncn2o, 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for N2O
    call rpn_comm_bcast(coefs%coef%ncco , 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for CO
    call rpn_comm_bcast(coefs%coef%ncch4, 1, 'MPI_INTEGER', 0, 'GRID', ierr)           ! number of coefficients for CH4

    call rpn_comm_bcast(coefs%coef%ws_nomega, 1, 'MPI_INTEGER', 0, 'GRID', ierr)

    call rpn_comm_bcast(coefs%coef%id_creation_date, 3, 'MPI_INTEGER', 0, 'GRID', ierr) 
    call rpn_comm_bcastc(coefs%coef%id_creation, 80, 'MPI_CHARACTER', 0, 'GRID', ierr) 
    call rpn_comm_bcastc(coefs%coef%id_Common_name, 32, 'MPI_CHARACTER', 0, 'GRID', ierr) 
    do i=1, 100
      call rpn_comm_bcastc(coefs%coef%line_by_line(i), 132, 'MPI_CHARACTER', 0, 'GRID', ierr) !ok
      call rpn_comm_bcastc(coefs%coef%readme_srf(i), 132, 'MPI_CHARACTER', 0, 'GRID', ierr) !ok
    end do
    call rpn_comm_bcastc(coefs%coef%fmv_model_def, 32, 'MPI_CHARACTER', 0, 'GRID', ierr)  !ok
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
    call rpn_comm_bcast(coefs%coef%mwcldtop, 1, 'MPI_INTEGER', 0, 'GRID', ierr) 
    ! then variable size vectors
    ! this one must be done first because it is used to dimension other ones ....
    if (mpi_myid > 0) allocate( coefs%coef%fmv_lvl(coefs%coef%fmv_gas))
    call rpn_comm_bcast(coefs%coef%fmv_lvl, coefs%coef%fmv_gas, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (mpi_myid > 0) then
      if (coefs%coef%nltecoef) allocate( coefs%coef%nlte_coef)
      allocate( coefs%coef%fmv_gas_id(coefs%coef%fmv_gas))
      allocate( coefs%coef%fmv_gas_pos(ngases_max)) !size is from rttov consts
      allocate( coefs%coef%fmv_var(coefs%coef%fmv_gas))
      allocate( coefs%coef%fmv_coe(coefs%coef%fmv_gas)) ! number of coefficients by gas (fmv_gas)
      allocate( coefs%coef%ws_npoint(coefs%coef%ws_nomega) )
      allocate( coefs%coef%ws_k_omega(coefs%coef%ws_nomega) )
      allocate( coefs%coef%ref_prfl_p(coefs%coef%fmv_lvl(gas_id_mixed) ) )
      allocate( coefs% coef % ref_prfl_t ( coefs%coef % fmv_lvl(gas_id_mixed), coefs%coef % fmv_gas ) )
      allocate( coefs% coef % ref_prfl_mr ( coefs%coef % fmv_lvl(gas_id_mixed), coefs%coef % fmv_gas ) )
      allocate( coefs% coef % bkg_prfl_mr ( coefs%coef % fmv_lvl(gas_id_mixed), coefs%coef % fmv_gas ) )
      allocate( coefs%coef%lim_prfl_p( coefs%coef%fmv_lvl(gas_id_mixed) ) )
      allocate( coefs%coef%lim_prfl_tmax( coefs%coef%fmv_lvl(gas_id_mixed) ) )
      allocate( coefs%coef%lim_prfl_tmin( coefs%coef%fmv_lvl(gas_id_mixed) ) )
      allocate( coefs%coef%lim_prfl_gmax( coefs%coef%fmv_lvl(gas_id_mixed), coefs%coef % fmv_gas ) )
      allocate( coefs%coef%lim_prfl_gmin( coefs%coef%fmv_lvl(gas_id_mixed), coefs%coef % fmv_gas ) )
      allocate( coefs%coef%env_prfl_tmax( coefs%coef%fmv_lvl(gas_id_mixed) ) )
      allocate( coefs%coef%env_prfl_tmin( coefs%coef%fmv_lvl(gas_id_mixed) ) )
      allocate( coefs%coef%env_prfl_gmax( coefs%coef%fmv_lvl(gas_id_mixed), coefs%coef % fmv_gas ) )
      allocate( coefs%coef%env_prfl_gmin( coefs%coef%fmv_lvl(gas_id_mixed), coefs%coef % fmv_gas ) )
      allocate( coefs%coef%dpp(0:coefs%coef%nlayers) )
      allocate( coefs%coef%dp(coefs%coef%nlayers) )
      allocate( coefs%coef%tstar(coefs%coef%nlayers) )
      allocate( coefs%coef%wstar(coefs%coef%nlayers) )
    end if

    call broadcastI41dArray( coefs%coef%fmv_gas_id )
    call broadcastI41dArray( coefs%coef%fmv_gas_pos )
    call broadcastI41dArray( coefs%coef%fmv_var )
    call broadcastI41dArray( coefs%coef%fmv_coe )
    call broadcastR81dArray( coefs%coef%ws_npoint )
    call broadcastR81dArray( coefs%coef%ws_k_omega )
    call broadcastR81dArray( coefs%coef%ref_prfl_p )
    call broadcastR82dArray( coefs%coef%ref_prfl_t )
    call broadcastR82dArray( coefs%coef%ref_prfl_mr )
    call broadcastR82dArray( coefs%coef%bkg_prfl_mr )
    call broadcastR81dArray( coefs%coef%lim_prfl_p )
    call broadcastR81dArray( coefs%coef%lim_prfl_tmax )
    call broadcastR81dArray( coefs%coef%lim_prfl_tmin )
    call broadcastR82dArray( coefs%coef%lim_prfl_gmax)
    call broadcastR82dArray( coefs%coef%lim_prfl_gmin)
    call broadcastR81dArray( coefs%coef%env_prfl_tmax )
    call broadcastR81dArray( coefs%coef%env_prfl_tmin )
    call broadcastR82dArray( coefs%coef%env_prfl_gmax )
    call broadcastR82dArray( coefs%coef%env_prfl_gmin )
    call broadcastR81dArray( coefs%coef%dp )
    call broadcastR81dArray( coefs%coef%dpp )
    call broadcastR81dArray( coefs%coef%tstar ) 
    call broadcastR81dArray( coefs%coef%wstar )

    if (coefs%coef%nozone > 0) then
      if (mpi_myid > 0) then
        allocate( coefs%coef%to3star(coefs%coef%nlayers) )
        allocate( coefs%coef%ostar(coefs%coef%nlayers) )
      end if
      call rpn_comm_bcast(coefs%coef%to3star, size(coefs%coef%to3star) , 'MPI_REAL8', 0, 'GRID', ierr) 
      call rpn_comm_bcast(coefs%coef%ostar, size(coefs%coef%ostar) , 'MPI_REAL8', 0, 'GRID', ierr)
    end if
    if ( coefs%coef%nco2 > 0) then
      if (mpi_myid>0) allocate( coefs%coef%co2star(coefs%coef%nlayers) )
      call rpn_comm_bcast(coefs%coef%co2star, size(coefs%coef%co2star) , 'MPI_REAL8', 0, 'GRID', ierr) 
    end if
    if (coefs%coef%nso2 > 0) then
      if (mpi_myid>0) allocate( coefs%coef%so2star(coefs%coef%nlayers) )
      call rpn_comm_bcast(coefs%coef%so2star, size(coefs%coef%so2star) , 'MPI_REAL8', 0, 'GRID', ierr)
    end if
    
    if (coefs%coef%fmv_model_ver == 9) then
      if (mpi_myid > 0) then
        allocate( coefs%coef%n2ostar(coefs%coef%nlayers) )
        allocate( coefs%coef%costar(coefs%coef%nlayers) )
        allocate( coefs%coef%ch4star(coefs%coef%nlayers) )
      end if
      call rpn_comm_bcast(coefs%coef%n2ostar, size(coefs%coef%n2ostar) , 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%costar, size(coefs%coef%costar) , 'MPI_REAL8', 0, 'GRID', ierr) 
      call rpn_comm_bcast(coefs%coef%ch4star, size(coefs%coef%ch4star) , 'MPI_REAL8', 0, 'GRID', ierr) 
    end if

    ! Fourth step: channel dependent parameters are extracted according to the channel list and sent to each MPI task
    coefs%coef%fmv_chn = size( channels )

    if (mpi_myid == 0) deallocate ( coefs%coef%ff_ori_chn)
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
    call extractR81dArray(coefs%coef%refl_visnir_fw, countUniqueChannel,indexchan)
    call extractR81dArray(coefs%coef%refl_visnir_ow, countUniqueChannel,indexchan)
    call extractCmplx81dArray(coefs%coef%woc_waopc_ow, countUniqueChannel, indexchan) 
    call extractCmplx81dArray(coefs%coef%woc_waopc_fw, countUniqueChannel, indexchan)
  
    call extractI41dArray(coefs%coef%fastem_polar, countUniqueChannel,indexchan)
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

    ! 2D arrays
    call extractR82dArray(coefs%coef%iremis_coef,coefs%coef%iremis_ncoef,countUniqueChannel,indexchan)
    ! 3D arrays
    call extractR83dArray(coefs%coef%pmc_coef,coefs%coef%pmc_nlay,countUniqueChannel,coefs%coef%pmc_nvar,indexchan)

    ! then coefficients. It is more complicated with RTTOV12

    if (mpi_myid > 0) then
      allocate (coefs%coef%thermal(countUniqueChannel) )
      do ichan = 1, countUniqueChannel
        allocate(coefs%coef%thermal(ichan)%gasarray(coefs%coef%fmv_gas) )
        nullify (coefs%coef%thermal(ichan)%mixedgas)
        nullify (coefs%coef%thermal(ichan)%watervapour)
        nullify (coefs%coef%thermal(ichan)%ozone)
        nullify (coefs%coef%thermal(ichan)%wvcont)
        nullify (coefs%coef%thermal(ichan)%co2)
        nullify (coefs%coef%thermal(ichan)%n2o)
        nullify (coefs%coef%thermal(ichan)%co)
        nullify (coefs%coef%thermal(ichan)%ch4)
        nullify (coefs%coef%thermal(ichan)%so2)
        do igas = 1, coefs%coef%fmv_gas
          allocate (coefs%coef%thermal(ichan)%gasarray(igas)%coef( coefs%coef%fmv_coe(igas), coefs%coef%nlayers) )
        end do
      end do
    end if
    
    do ichan = 1, countUniqueChannel
      do igas = 1, coefs%coef%fmv_gas
        call broadcastR82dArray( coefs%coef%thermal(ichan)%gasarray(igas)%coef )
      end do
    end do

    if (coefs%coef%solarcoef) then
      if (mpi_myid>0) then
        allocate (coefs%coef%solar(countUniqueChannel) )
        do ichan = 1, countUniqueChannel
          allocate (coefs%coef%solar(ichan)%gasarray(coefs%coef%fmv_gas) )
          nullify (coefs%coef%solar(ichan)%mixedgas)
          nullify (coefs%coef%solar(ichan)%watervapour)
          nullify (coefs%coef%solar(ichan)%ozone)
          nullify (coefs%coef%solar(ichan)%wvcont)
          nullify (coefs%coef%solar(ichan)%co2)
          nullify (coefs%coef%solar(ichan)%n2o)
          nullify (coefs%coef%solar(ichan)%co)
          nullify (coefs%coef%solar(ichan)%ch4)
          nullify (coefs%coef%solar(ichan)%so2)
          do igas = 1, coefs%coef%fmv_gas
            allocate (coefs%coef%solar(ichan)%gasarray(igas)%coef( coefs%coef%fmv_coe(igas), coefs%coef%nlayers) )
          end do
        end do
      end if
      
      do ichan = 1, countUniqueChannel
        do igas = 1, coefs%coef%fmv_gas
          call broadcastR82dArray( coefs%coef%solar(ichan)%gasarray(igas)%coef )
        end do
      end do
    end if
  
    allocate(bigArray(countUniqueChannel,maxval(coefs%coef%fmv_coe),coefs%coef%fmv_gas,coefs%coef%nlayers) )
    bigArray(:,:,:,:) = 0.0d0

    do ichan = 1, countUniqueChannel  
      do igas = 1, coefs%coef%fmv_gas
        associated0 = associated( coefs%coef%thermal(ichan)%gasarray(igas)%coef )
        call rpn_comm_bcast(associated0, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
        if (associated0) then
          do l=1, coefs%coef%nlayers
            do k=1, coefs%coef%fmv_coe(igas)
              bigArray(ichan,k,igas,l) = coefs%coef%thermal(ichan)%gasarray(igas)%coef(k,l)
            end do
          end do
        end if
      end do
    end do
    
    do ichan = 1, countUniqueChannel
      do igas = 1, coefs%coef%fmv_gas
        associated0 = associated(coefs%coef%thermal(ichan)%gasarray(igas)%coef)
        if (associated0) deallocate(coefs%coef%thermal(ichan)%gasarray(igas)%coef)
      end do
      deallocate(coefs%coef%thermal(ichan)%gasarray)
    end do
    deallocate(coefs%coef%thermal)
    
    allocate (coefs%coef%thermal(coefs%coef%fmv_chn) )
    do ichan = 1, coefs%coef%fmv_chn
      allocate (coefs%coef%thermal(ichan)%gasarray(coefs%coef%fmv_gas) )
      
      nullify (coefs%coef%thermal(ichan)%mixedgas)
      nullify (coefs%coef%thermal(ichan)%watervapour)
      nullify (coefs%coef%thermal(ichan)%ozone)
      nullify (coefs%coef%thermal(ichan)%wvcont)
      nullify (coefs%coef%thermal(ichan)%co2)
      nullify (coefs%coef%thermal(ichan)%n2o)
      nullify (coefs%coef%thermal(ichan)%co)
      nullify (coefs%coef%thermal(ichan)%ch4)
      nullify (coefs%coef%thermal(ichan)%so2)

      do igas = 1, coefs%coef%fmv_gas
        if (any( bigArray(indexchan(ichan),:,igas,:) /= 0.) ) then
          allocate (coefs%coef%thermal(ichan)%gasarray(igas)%coef( coefs%coef%fmv_coe(igas), coefs%coef%nlayers) )
          do l=1, coefs%coef%nlayers
            do k=1, coefs%coef%fmv_coe(igas)
              coefs%coef%thermal(ichan)%gasarray(igas)%coef(k,l)  = bigArray(indexchan(ichan),k,igas,l)
            end do
          end do
        end if
        call set_pointers(coefs%coef%thermal(ichan), igas, coefs%coef%fmv_gas_id(igas))
      end do
    end do

    if (coefs % coef % solarcoef) then
      bigArray(:,:,:,:) = 0.0d0
      do ichan = 1, countUniqueChannel  
        do igas = 1, coefs%coef%fmv_gas
          associated0 = associated( coefs%coef%solar(ichan)%gasarray(igas)%coef )
          call rpn_comm_bcast(associated0, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
          if (associated0) then
            do l=1,coefs % coef % nlayers
              do k=1,coefs % coef % fmv_coe(igas)
                bigArray(ichan,k,igas,l) = coefs % coef % solar(ichan) % gasarray(igas) % coef(k,l)
              end do
            end do
          end if
        end do
      end do
  
      do ichan = 1, countUniqueChannel
        do igas = 1, coefs%coef%fmv_gas
          associated0 = associated(coefs%coef%solar(ichan)%gasarray(igas)%coef)
          if (associated0) deallocate(coefs%coef%solar(ichan)%gasarray(igas)%coef)
        end do
        deallocate(coefs%coef%solar(ichan)%gasarray)   
      end do
      deallocate(coefs%coef%solar)
     
      allocate (coefs%coef%solar(coefs%coef%fmv_chn) )
      do ichan = 1, coefs%coef%fmv_chn
        allocate (coefs%coef%solar(ichan)%gasarray(coefs%coef%fmv_gas) )
        nullify (coefs%coef%solar(ichan)%mixedgas)
        nullify (coefs%coef%solar(ichan)%watervapour)
        nullify (coefs%coef%solar(ichan)%ozone)
        nullify (coefs%coef%solar(ichan)%wvcont)
        nullify (coefs%coef%solar(ichan)%co2)
        nullify (coefs%coef%solar(ichan)%n2o)
        nullify (coefs%coef%solar(ichan)%co)
        nullify (coefs%coef%solar(ichan)%ch4)
        nullify (coefs%coef%solar(ichan)%so2)
        do igas = 1, coefs%coef%fmv_gas
          if (any( bigArray(indexchan(ichan),:,igas,:) /= 0.) ) then
            allocate (coefs%coef%solar(ichan)%gasarray(igas)%coef( coefs%coef%fmv_coe(igas), coefs%coef%nlayers) )
            do l=1, coefs%coef%nlayers
              do k=1, coefs%coef%fmv_coe(igas)
                coefs%coef%solar(ichan)%gasarray(igas)%coef(k,l)  = bigArray(indexchan(ichan),k,igas,l)
              end do
            end do
          end if
          call set_pointers(coefs%coef%solar(ichan), igas, coefs%coef%fmv_gas_id(igas))
        end do
      end do
    end if

    deallocate(bigArray)
 
    if (mpi_myid==0 .and. associated(coefs%coef%bounds) )  deallocate(coefs%coef%bounds)
  

    !Allocate bounds array to store opdep calculation layer limits
    !1st dim: upper boundary layer [ub](above which coefs all zeros), lower boundary layer [lb]
    !4th dim: thermal layer limits, solar layer limits
    allocate(coefs%coef%bounds(2, coefs%coef%fmv_gas, coefs%coef%fmv_chn, 2))
    call set_fastcoef_level_bounds(coefs%coef, coefs%coef%thermal, thermal = .true._jplm)
    ! If the SOLAR_FAST_COEFFICIENTS section is not present then point the solar coefs to the thermal coefs
    if (coefs%coef%solarcoef) then
      call set_fastcoef_level_bounds(coefs%coef, coefs%coef%solar, thermal = .false._jplm)
    else
      coefs%coef%solar => coefs%coef%thermal
      coefs%coef%bounds(:,:,:,2) = coefs%coef%bounds(:,:,:,1)
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

      ! Reset NLTE channel variables according to input channel list
      coefs%coef%nlte_coef%start_chan  = nlte_start
      coefs%coef%nlte_coef%nchan      = nlte_count

      if (mpi_myid > 0) then
         allocate (coefs%coef%nlte_coef%sec_sat(coefs%coef%nlte_coef%nsat) )
         allocate (coefs%coef%nlte_coef%sol_zen_angle(coefs%coef%nlte_coef%nsol) )
      end if
      call rpn_comm_bcast(coefs%coef%nlte_coef%sec_sat, coefs%coef%nlte_coef%nsat, 'MPI_REAL8', 0, 'GRID', ierr)
      call rpn_comm_bcast(coefs%coef%nlte_coef%sol_zen_angle, coefs%coef%nlte_coef%nsol, 'MPI_REAL8', 0, 'GRID', ierr)
     
      allocate(bigArray(coefs%coef%nlte_coef%ncoef, coefs%coef%nlte_coef%nsat, &
           coefs%coef%nlte_coef%nsol, nlte_file_nchan) )

      if (mpi_myid == 0) then
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

    deallocate(listGlobal)
 
  end subroutine tvs_rttov_read_coefs

  subroutine extractI41dArray(array,oldSize,index)
    implicit none
    integer, pointer :: array(:)
    integer, intent(in) :: oldSize
    integer, intent(in) :: index(:)
    ! Locals
    integer :: newSize, tmpI41d(oldSize), ierr, trueSize

    if (mpi_myid == 0) then
      if (associated(array)) then
        trueSize = size(array)
      else
        trueSize = 0
      end if
    end if
    ierr = 0
    call rpn_comm_bcast(trueSize, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) "error 1 in rpn_comm_bcast", ierr, trueSize
      call utl_abort("extractI41dArray")
    end if
    if (trueSize < 1) return
    
    if (trueSize /= oldSize) then
      write(*,*) "extractI41dArray: should not happen ", trueSize, oldSize
    end if
    
    newSize = size( index )

    if (mpi_myid > 0) allocate( array(oldSize))
    ierr = 0
    call rpn_comm_bcast(array, oldSize, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) "error 2 in rpn_comm_bcast", ierr, array(:)
      call utl_abort("extractI41dArray")
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
    
    if (mpi_myid == 0) then
      if (associated(array)) then
        trueSize = size(array)
      else
        trueSize = 0
      end if
    end if
    ierr = 0
    call rpn_comm_bcast(trueSize, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) "error 1 in rpn_comm_bcast", ierr, trueSize
      call utl_abort("extractR81dArray")
    end if
    if (trueSize < 1) return
    
    if (trueSize /= oldSize) then
      write(*,*) "extractR81dArray: should not happen ", trueSize, oldSize
    end if
    
    newSize = size( index )

    if (mpi_myid > 0) allocate( array(oldSize))
    ierr = 0
    call rpn_comm_bcast(array, oldSize, 'MPI_REAL8', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) "error 2 in rpn_comm_bcast", ierr, array(:)
      call utl_abort("extractR81dArray")
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
    
    if (mpi_myid == 0) then
      if (associated(array)) then
        trueSize = size(array)
      else
        trueSize = 0
      end if
    end if
    ierr = 0
    call rpn_comm_bcast(trueSize, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) "error 1 in rpn_comm_bcast", ierr, trueSize
      call utl_abort("extractR82dArray")
    end if
    if (trueSize < 1) return
    
    if (trueSize /= oldSize1 * oldSize2) then
      write(*,*) "extractR82dArray: should not happen ", trueSize, oldSize1, oldSize2
    end if

    newSize = size( index )

    if (mpi_myid > 0) allocate( array(oldSize1,oldSize2) )
    ierr = 0
    call rpn_comm_bcast(array, oldSize1*oldSize2, 'MPI_REAL8', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) "error 2 in rpn_comm_bcast", ierr, array(:,:)
      call utl_abort("extractR82dArray")
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
    
    if (mpi_myid == 0) then
      if (associated(array)) then
        trueSize = size(array)
      else
        trueSize = 0
      end if
    end if
    ierr = 0
    call rpn_comm_bcast(trueSize, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) "error 1 in rpn_comm_bcast", ierr, trueSize
      call utl_abort("extractR83dArray")
    end if
    if (trueSize < 1) return

    if (trueSize /= oldSize1 * oldSize2 * oldSize3) then
      write(*,*) "extractR83dArray: should not happen ", trueSize, oldSize1, oldSize2, oldSize3
    end if
  
    newSize = size( index )
  
    if (mpi_myid > 0) allocate( array(oldSize1,oldSize2, oldSIze3) )
    ierr = 0
    call rpn_comm_bcast(array, oldSize1*oldSize2*oldSize3, 'MPI_REAL8', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) "error 2 in rpn_comm_bcast", ierr, array(:,:,:)
      call utl_abort("extractR83dArray")
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

    if (mpi_myid == 0) then
      if (associated(array)) then
        trueSize = size(array)
      else
        trueSize = 0
      end if
    end if
    ierr = 0
    call rpn_comm_bcast(trueSize, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) "error 1 in rpn_comm_bcast", ierr, trueSize
      call utl_abort("extractCmplx81dArray")
    end if
    if (trueSize < 1) return
  
    if (trueSize /= oldSize) then
      write(*,*) "extractCmplx81dArray: should not happen ", trueSize, oldSize
    end if

    newSize = size( index )
    
    if (mpi_myid > 0) allocate( array(oldSize))
    ierr = 0
    call rpn_comm_bcast(array, oldSize, 'MPI_COMPLEX8', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) "error 2 in rpn_comm_bcast", ierr, array(:)
      call utl_abort("extractCmplx81dArray")
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
      write(*,*) "error 1 in rpn_comm_bcast", ierr, associated0
      call utl_abort("broadcastR82dArray")
    end if
    ierr = 0
    if (associated0) call rpn_comm_bcast(array, size(array) , 'MPI_REAL8', 0, 'GRID', ierr)
    if (ierr /= 0) then
      write(*,*) "error 2 in rpn_comm_bcast", ierr, size(array,dim=1), size(array,dim=2)
      call utl_abort("broadcastR82dArray")
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
      write(*,*) "error 1 in rpn_comm_bcast", ierr, associated0
      call utl_abort("broadcastR81dArray")
    end if
    ierr = 0
    if (associated0) call rpn_comm_bcast(array, size(array) , 'MPI_REAL8', 0, 'GRID', ierr)
    if (ierr/=0) then
      write(*,*) "error 2 in rpn_comm_bcast", ierr, size(array)
      call utl_abort("broadcastR81dArray")
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
      write(*,*) "error 1 in rpn_comm_bcast", ierr, associated0
      call utl_abort("broadcastI41dArray")
    end if
    ierr = 0
    if (associated0) call rpn_comm_bcast(array, size(array) , 'MPI_INTEGER', 0, 'GRID', ierr)
    if (ierr/=0) then
      write(*,*) "error 2 in rpn_comm_bcast", ierr, size(array)
      call utl_abort("broadcastI41dArray")
    end if

  end subroutine broadcastI41dArray

  !--------------------------------------------------------------------------
  !  tvs_calc_jo
  !--------------------------------------------------------------------------
  subroutine tvs_calc_jo(pjo,llprint,lobsSpaceData,dest_obs)
    !
    ! *Purpose*: Computation of Jo and the residuals to the tovs observations
    !
    implicit none

    !Arguments:
    real(8),          intent(out)   :: pjo          ! Computed Tovs observation cost fucntion
    logical,          intent(in)    :: llprint      ! Logical flag to control printout by channel
    type(struct_obs), intent(inout) :: lobsSpaceData! obsSpacaData structure
    integer,          intent(in)    :: dest_obs     ! obsSpaceData body destinationcolumn (ex: OBS_OMP or OBS_OMA)

    ! Locals:
    integer :: sensorIndex, channelIndex, tovsIndex
    real*8 dlsum, zdtb !, zjon,zgami,zqcarg
    real*8 zjoch  (0:tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    real*8 zavgnrm(0:tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    integer j, i, nchanperline, indxs, indxe
    integer inobsjo, incanjo
    integer idatyp
    integer channelNumber, ichobs_a
    integer inobsch(0:tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    integer lcanjo(tvs_maxChannelNumber)
    integer :: headerIndex, bodyIndex
    real (8) :: x(tvs_maxChannelNumber),y(tvs_maxChannelNumber)
    real (8) :: sigmaObs
    integer :: list_chan(tvs_maxChannelNumber)
    integer :: count

    if ( llprint ) write(*,*) "Entering tvs_calc_jo subroutine"

    pjo = 0.d0

    if ( tvs_nobtov == 0) return    ! exit if there are not tovs data

    ! 1.  Computation of (hx - z)/sigma for tovs data only

    dlsum    = 0.d0
    inobsjo  = 0
    inobsch(:,:) = 0
    zjoch  (:,:) = 0.0d0
    zavgnrm(:,:) = 0.0d0

    ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(lobsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(lobsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! 1.1  Extract general information for this observation point
      !      ------------------------------------------------------

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(lobsSpaceData,OBS_ITY,headerIndex)
      if ( .not. tvs_isIdBurpTovs(idatyp) ) cycle HEADER

      tovsIndex = tvs_tovsIndex(headerIndex)
      if ( tovsIndex == -1 ) cycle HEADER
       
      sensorIndex = tvs_lsensor(tovsIndex)

      ! Set the body list
      ! (& start at the beginning of the list)
      call obs_set_current_body_list(lobsSpaceData, headerIndex)
      count = 0
      BODY: do 
        bodyIndex = obs_getBodyIndex(lobsSpaceData)
        if (bodyIndex < 0 ) then
          if (count > 0 .and. rmat_lnondiagr) then
            call rmat_sqrtRm1(sensorIndex,count,x(1:count),y(1:count),list_chan(1:count),tovsIndex)
            dlsum =  dlsum + 0.5d0*dot_product(y(1:count),y(1:count))
          end if
          exit BODY
        end if

        ! Only consider if flagged for assimilation
        if ( obs_bodyElem_i(lobsSpaceData,obs_aSS,bodyIndex) /= obs_assimilated ) cycle BODY                

        channelNumber = nint(obs_bodyElem_r(lobsSpaceData,OBS_PPP,bodyIndex))
        channelNumber = max( 0 , min( channelNumber , tvs_maxChannelNumber + 1))
        ichobs_a = channelNumber
        channelNumber = channelNumber - tvs_channelOffset(sensorIndex)
        channelIndex = utl_findArrayIndex(tvs_ichan(:,sensorIndex),tvs_nchan(sensorIndex),channelNumber)
        if ( channelIndex == 0 ) then
          write(*,'(A)') ' tvs_calc_jo: error with channel number'
          call utl_abort('tvs_calc_jo')
        end if

        zdtb = obs_bodyElem_r(lobsSpaceData,OBS_PRM,bodyIndex) - &
             tvs_radiance (tovsIndex) % bt(channelIndex)
        if ( tvs_debug ) then
          write(*,'(a,i4,2f8.2,f6.2)') ' channelNumber,sim,obs,diff= ', &
               channelNumber,  tvs_radiance (tovsIndex) % bt(channelIndex), &
               obs_bodyElem_r(lobsSpaceData,OBS_PRM,bodyIndex), -zdtb
        end if
        call obs_bodySet_r(lobsSpaceData,dest_obs,bodyIndex, zdtb)

        sigmaObs = obs_bodyElem_r(lobsSpaceData,OBS_OER,bodyIndex)

        if ( sigmaObs == MPC_missingValue_R8) cycle body

        ! Comment out the modification of Jobs due to varqc for now, since this is probably
        ! only needed for use of nonlinear obs operator in minimization, which is not yet
        ! functional, but this interferes with doing ensemble of analyses (M. Buehner, Dec. 2013)
        !if (.not. min_lvarqc .or. obs_bodyElem_r(lobsSpaceData,OBS_POB,bodyIndex).eq.0.0d0) then
        dlsum =  dlsum &
             + (zdtb * zdtb) / (2.d0 * sigmaObs * sigmaObs)
        !else
        !  compute contribution of data with varqc
 
        !   zgami = obs_bodyElem_r(lobsSpaceData,OBS_POB,bodyIndex)
        !   zjon = (zdtb* &
        !           zdtb)/2.d0
        !   zqcarg = zgami + exp(-1.0d0*zjon)
        !   dlsum= dlsum - log(zqcarg/(zgami+1.d0))
        !end if
        count = count + 1
        x(count) = zdtb
        list_chan(count) = channelNumber

        inobsjo = inobsjo + 1
        inobsch(ichobs_a,SensorIndex) = inobsch(ichobs_a,SensorIndex) + 1
        zjoch(ichobs_a,SensorIndex)   = &
             zjoch(ichobs_a,SensorIndex) &
             + zdtb * zdtb / (sigmaObs * sigmaObs)
        zavgnrm(ichobs_a,SensorIndex)   = &
             zavgnrm(ichobs_a,SensorIndex) - &
             zdtb / sigmaObs
      end do BODY

    end do HEADER

    !   2.  Close up, print summary
    !   .   -----------------------

    pjo = dlsum

    if ( pjo == 0.d0) return


    ! printout of mean jo and normalized average for each sensor.

    nchanperline = 18
    if ( llprint .and. inobsjo > 0 ) then
      write(*,*)
      write(*,*)
      write(*,'(10x,A)') "-tvs_calc_jo: computing jo and residuals to tovs  observations"

      do sensorIndex = 1, tvs_nsensors
        do i = 1, tvs_maxChannelNumber
          inobsch(0,sensorIndex) = inobsch(0,sensorIndex) + &
               inobsch(i,sensorIndex)
          zjoch(0,sensorIndex)   = zjoch(0,sensorIndex) + &
               zjoch(i,sensorIndex)
          zavgnrm(0,sensorIndex) = zavgnrm(0,sensorIndex) + &
               zavgnrm(i,sensorIndex)
        end do
      end do

      do sensorIndex = 1, tvs_nsensors
        incanjo = 0
        do i = 0, tvs_maxChannelNumber
          if ( inobsch(i,sensorIndex) /= 0 ) then
            incanjo = incanjo + 1
            lcanjo(incanjo) = i
          end if
        end do
        if ( incanjo /= 0 ) then
          write(*,'(/1x,"sensor #",i2,". platform: ",a, "instrument: ",a)') &
               sensorIndex, tvs_satelliteName(sensorIndex), tvs_instrumentName(sensorIndex)
          do j = 1, incanjo, nchanperline
            indxs = j
            indxe = min(j + nchanperline - 1 , incanjo)
            if ( j == 1 ) then
              write(*,'(1x,"channel",t13,"   all",17i6)') (lcanjo(i),i=indxs+1,indxe)
            else
              write(*,'(1x,"channel",t13,18i6)') (lcanjo(i),i=indxs,indxe)
            end if
            write(*,'(1x,"no. obs.",t13,18i6)') (inobsch(lcanjo(i),sensorIndex),i=indxs,indxe)
            write(*,'(1x,"mean jo",t13,18f6.2)') &
                 (zjoch(lcanjo(i),sensorIndex)/max(1,inobsch(lcanjo(i),sensorIndex)) ,i=indxs,indxe)
            write(*,'(1x,"norm. bias",t13,18f6.2,/)') &
                 (zavgnrm(lcanjo(i),sensorIndex)/max(1,inobsch(lcanjo(i),sensorIndex)) ,i=indxs,indxe)
          end do
        end if
      end do
    end if

  end subroutine tvs_calc_jo

  !--------------------------------------------------------------------------
  !  tvs_getChannelIndexFromChannelNumber
  !--------------------------------------------------------------------------
  subroutine tvs_getChannelIndexFromChannelNumber(idsat,chanIndx,chanNum)
    !
    ! *Purpose*: Get channel Index from channel number for given intrument
    !
    implicit none

    !Arguments:
    integer, intent(in)  :: idsat   ! Satellite index
    integer, intent(out) :: chanIndx! Channel index
    integer, intent(in)  :: chanNum ! Channel number

    ! Locals:
    logical, save              :: first =.true.
    integer                    :: channelNumber, sensorIndex, channelIndex 
    integer, allocatable, save :: index(:,:)

    if (first) then
      allocate( index(tvs_nsensors, tvs_maxChannelNumber ) )
      index(:,:) = -1
      do sensorIndex = 1, tvs_nsensors
        channels:do channelNumber = 1,  tvs_maxChannelNumber
          indexes: do channelIndex =1, tvs_nchan(sensorIndex)
            if ( channelNumber == tvs_ichan(channelIndex,sensorIndex) ) then
              index(sensorIndex,channelNumber) = channelIndex
              exit indexes
            end if
          end do indexes
        end do channels
      end do
      first = .false.
    end if

    chanIndx = index(idsat,chanNum)

  end subroutine tvs_getChannelIndexFromChannelNumber

end module tovs_nl_mod

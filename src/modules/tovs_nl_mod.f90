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

!--------------------------------------------------------------------------
!! MODULE tovs_nl (prefix= no defined prefix)
!!
!! *Purpose*: Derived types, public variables and procedures related to the nonlinear
!!            version of RTTOV
!!
!--------------------------------------------------------------------------
module tovs_nl_mod

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
  use hirchannels_mod
  use obsSpaceData_mod
  use EarthConstants_mod
  use MathPhysConstants_mod
  use ozoneclim_mod
  use columnData_mod 
  use presProfileOperators_mod
  use tovs_extrap_mod
  use mod_rttov_emis_atlas
  use rMatrix_mod 

  implicit none
  save
  private

  type surface_params
    real(8)              :: ALBEDO   ! surface albedo (0-1)
    real(8)              :: ICE      ! ice cover (0-1) 
    real(8)              :: SNOW     ! snow cover (0-1)
    real(8)              :: PCNT_WAT ! water percentage in pixel containing profile (0-1)
    real(8)              :: PCNT_REG ! water percentage in an area around profile (0-1)
    integer              :: LTYPE    ! surface type (1,...,20)
  end type surface_params

  ! public derived type through inheritance (from module rttov_types)
  public :: rttov_radiance, rttov_profile, rttov_chanprof, rttov_coefs, rttov_transmission, rttov_options, rttov_emissivity

  ! public variables (parameters)
  public :: tvs_maxChannelNumber, tvs_maxNumberOfChannels, tvs_maxNumberOfSensors,  tvs_nlevels, tvs_mesosphereLapseRate
  ! public variables (non-parameters)
  public :: tvs_nchan, tvs_ichan, tvs_lsensor, tvs_lobsno, tvs_ltovsno, tvs_nobtov
  public :: tvs_isReallyPresent,tvs_listSensors
  public :: tvs_nsensors, tvs_platforms, tvs_satellites, tvs_instruments, tvs_channelOffset
  public :: tvs_debug, tvs_satelliteName, tvs_instrumentName
  public :: platform_name, inst_name ! (from rttov)
  public :: tvs_coefs, tvs_opts, tvs_profiles, tvs_transmission,tvs_emissivity
  public :: tvs_radiance, tvs_surfaceParameters

  ! public procedures
  public :: tvs_fillProfiles, tvs_rttov, tvs_calc_jo
  public :: tvs_setupAlloc,tvs_setup, tvs_isIdBurpTovs, tvs_isIdBurpInst
  public :: tvs_getInstrumentId, tvs_getPlatformId, tvs_mapSat, tvs_mapInstrum
  public :: tvs_isInstrumHyperSpectral, tvs_getChanprof, tvs_countRadiances
  public :: tvs_getHIREmissivities, tvs_getOtherEmissivities, tvs_rttov_read_coefs
  ! Module parameters
  ! units conversion from  mixing ratio to ppmv and vice versa
  real(8), Parameter :: qMixratio2ppmv  = (1000000.0d0 * mair) / mh2o
  real(8), Parameter :: qppmv2Mixratio  = mh2o / (1000000.0d0 * mair)
  real(8), Parameter :: o3Mixratio2ppmv = (1000000.0d0 * mair) / mo3
  real(8), Parameter :: o3ppmv2Mixratio = mo3 / (1000000.0d0 * mair)

  integer, parameter :: tvs_maxChannelNumber   = 8461   ! Max. value for channel number
  integer, parameter :: tvs_maxNumberOfChannels = 1305  ! Max. no. of channels (for one profile/spectra)
  integer, parameter :: tvs_maxNumberOfSensors  = 40    ! Max no sensors to be used
  integer, parameter :: tvs_nlevels     = 54            ! No. of RTTOV pressure levels including "rttov top" at 0.005 hPa
!**********************************************************
! S. Heilliette this parameter was computed from the mean lapse rate between 50 km and 85 km
! of the US standard atmosphere from data contained in "AFGL Atmospheric Constituent Profiles (0-120km)"
! AFGL 1986, G.P. ANDERSON, J.H. CHETWYND, S.A. CLOUGH, E. P. SHETTLE and F.X. KNEIZYS
  real(8), parameter :: tvs_mesosphereLapseRate=16.2d0
! unit is K/log(P), a positive value corresponds to decrease of temperature with height
!**********************************************************

  ! Module variables
  integer, allocatable :: tvs_nchan(:)              ! number of channels per instrument
  integer, allocatable :: tvs_ichan(:,:)            ! list of channels per instrument
  integer, allocatable :: tvs_lsensor(:)            ! sensor number for each profile
  integer, allocatable :: tvs_lobsno (:)            ! observation number in cma for each profile
  integer, allocatable :: tvs_ltovsno (:)           ! index in TOVS structures for each observation in cma
  logical, allocatable :: tvs_isReallyPresent(:)   ! logical flag to identify instruments really assimilated
  integer, allocatable :: tvs_listSensors(:,:)     ! sensor list
  type( surface_params ) , allocatable :: tvs_surfaceParameters(:)    ! surface parameters 
  integer tvs_nobtov

!   Variables from comtov.cdk
!     tvs_nsensors           : number of individual sensors.
!     tvs_platforms(MXPLATFORM)  : platform ID's (e.g., 1=NOAA; 2=DMSP; ...)
!     tvs_satellites(tvs_maxNumberOfSensors)  : satellite ID's (e.g., 1 to 16 for NOAA; ...)
!     tvs_instruments(tvs_maxNumberOfSensors) : instrument ID's (e.g., 3=AMSU-A; 4=AMSU-B; 6=SSMI; ...)
!     tvs_channelOffset(tvs_maxNumberOfSensors) : BURP to RTTOV channel mapping offset
!     tvs_debug            : logical key controlling statements to be
!     .                    executed while debugging TOVS only
!     radiativeTransferCode            : TOVS radiation model used:
!                             RTTOV, EUMETSAT NWP SAF radiation model
  integer tvs_nsensors
  integer tvs_platforms(tvs_maxNumberOfSensors), tvs_satellites(tvs_maxNumberOfSensors)
  integer tvs_instruments(tvs_maxNumberOfSensors), tvs_channelOffset(tvs_maxNumberOfSensors)
  logical tvs_debug, useUofWIREmiss
  character(len=15) tvs_satelliteName(tvs_maxNumberOfSensors), tvs_instrumentName(tvs_maxNumberOfSensors)
  character(len=8) radiativeTransferCode
  real(8) , allocatable :: tvs_emissivity(:,:)   ! surface emissivities organized by profiles and channels
! CERES file dimension in grid points
  integer, parameter :: KSLON=2160, KSLAT=1080

! Variables on standard files
  integer ::  JTYPE(KSLON,KSLAT)       ! surface type
  real(8) :: WATERF(KSLON,KSLAT)       ! water fraction


  ! Derived types

  type( rttov_coefs ),        allocatable :: tvs_coefs(:)          ! coefficients
  type( rttov_options ),      allocatable :: tvs_opts(:)           ! options
  type( rttov_profile ),      allocatable :: tvs_profiles(:)       ! profiles, all profiles
  type( rttov_radiance ),     allocatable :: tvs_radiance(:)       ! radiances organized by profile
  type( rttov_transmission ), allocatable :: tvs_transmission(:)   ! transmittances all profiles for HIR quality control

  integer, external :: get_max_rss
 
contains

!--------------------------------------------------------------------------
!!
!!  s/r tvs_setupAlloc : Memory allocation for the non linear radiative transfer model
!!                 variables.
!!          (original name of routine: sutovalo)
!!
!!Author  : J. Halle *CMDA/AES Oct 1999
!!
!!     Purpose: to allocate memory for the radiative transfer model variables.
!!
!! Revision:
!!           S. Pellerin *ARMA/SMC May 2000
!!            - Fix for F90 conversion
!!           C. Chouinard *ARMA/SMC Aug 2000
!!            - remove reference to nincrem in memory allocation
!!           JM Belanger *CMDA/SMC!  aug 2000
!!            - 32 bits conversion
!!           J. Halle *CMDA/AES  dec 2000
!!            - adapt to TOVS level 1b.
!!           J. Halle CMDA/SMC May 2002
!!            - adapt to RTTOV-7 code
!!           J. Halle CMDA/SMC Feb 2003
!!            - add codtyp for AMSUB (=181).
!!           J. Halle CMDA/SMC Nov 2004
!!            - adapt to RTTOV-8;
!!            - convert to Fortran 90.
!!           A. Beaulne CMDA/SMC June 2006
!!            - modifications for AIRS
!!            - allocation of ozone profiles
!!           R. Sarrazin  CMDA   April 2008
!!            - adapt to CSR
!!           S. Heilliette
!!            - adapt to IASI
!!            - adapt to rttov 10.0 (october 2010)
!!           S. Macpherson
!!            - adapt to ATMS (codtyp 192)
!!           S.  Heilliette
!!            - adapt to CrIS (codtyp 193)
!--------------------------------------------------------------------------

  subroutine tvs_setupAlloc(lobsSpaceData)
    implicit none

    type(struct_obs) :: lobsSpaceData

    integer :: alloc_status(6)

    integer ::  ival, IPLATFORM, ISAT, INSTRUM, KRTID

    integer ::  JO, IDATYP, J, JI, JK
    integer ::  ISENS, NC, NL
    integer ::  ICHN, nosensor, INDXCHN
    integer ::  ERRORSTATUS(1),ASW
    integer ::  index_header, index_body

    if (tvs_nsensors == 0) return

    !     1. Determine the number of radiances to be assimilated.
    !        Construct a list of channels for each sensor.
    !        Construct a list of sensor number for each profile
    !     .  ---------------------------------------------------

    write(*,*) "Entering tvs_setupAlloc" 

    alloc_status = 0
    allocate (tvs_nchan(tvs_nsensors),                         stat= alloc_status(1))
    allocate (tvs_ichan(tvs_maxNumberOfChannels,tvs_nsensors), stat= alloc_status(2))
    allocate (tvs_lsensor(obs_numheader(lobsSpaceData)),       stat= alloc_status(3))
    allocate (tvs_lobsno (obs_numheader(lobsSpaceData)),       stat= alloc_status(4))
    allocate (tvs_ltovsno(obs_numheader(lobsSpaceData)),       stat= alloc_status(5))
    allocate (tvs_isReallyPresent(tvs_nsensors),               stat= alloc_status(6))

    call utl_checkAllocationStatus(alloc_status, " tvs_setupAlloc")
  
    tvs_nchan(:) = 0 
    tvs_ichan(:,:) = 0
    tvs_ltovsno(:) = 0
    tvs_isReallyPresent(:) = .true.

    tvs_nobtov = 0

    ! loop over all header indices of the 'TO' family
    ! Set the header list & start at the beginning of the list
    call obs_set_current_header_list(lobsSpaceData,'TO')
    HEADER: do
      index_header = obs_getHeaderIndex(lobsSpaceData)
      if (index_header < 0) exit HEADER

      IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,index_header)

      if ( .not. tvs_isIdBurpTovs(IDATYP) ) cycle HEADER   ! Proceed to the next header_index

      tvs_nobtov = tvs_nobtov + 1
     
      !    Construct list of channels for each sensor:
      !          map burp satellite info to RTTOV platform and satellite.
      IVAL = obs_headElem_i(lobsSpaceData,OBS_SAT,index_header)
      call tvs_mapSat(IVAL,IPLATFORM,ISAT)
      if (IPLATFORM == -1) then
        write(*,*) "Unknown OBS_SAT !",IVAL
        call utl_abort('tvs_setupAlloc')
      end if
      !    map burp instrument info to RTTOV instrument.
      IVAL = obs_headElem_i(lobsSpaceData,OBS_INS,index_header)
      call tvs_mapInstrum(IVAL,INSTRUM)
      if (INSTRUM == -1) then
        write(*,*) "Unknown OBS_INS !",IVAL
        call utl_abort('tvs_setupAlloc')
      end if
      !    find sensor number for this obs.
      nosensor =0
      do KRTID = 1, tvs_nsensors
        if ( IPLATFORM == tvs_platforms  (KRTID) .AND. &
             ISAT      == tvs_satellites (KRTID) .AND. &
             INSTRUM   == tvs_instruments(KRTID)      ) THEN
          nosensor = KRTID
          exit
        end if
      end do

      if (nosensor > 0) then
        tvs_lsensor(tvs_nobtov) = nosensor
        tvs_lobsno (tvs_nobtov) = index_header
        tvs_ltovsno (index_header) = tvs_nobtov
      else
        write(*,*) "IPLATFORM,ISAT,INSTRUM ",IPLATFORM,ISAT,INSTRUM
        write(*,'(A)') ' tvs_setupAlloc: Invalid Sensor'
        do KRTID = 1, tvs_nsensors
          print *,krtid,tvs_platforms  (KRTID),tvs_satellites (KRTID),tvs_instruments(KRTID)
        end do
        call utl_Abort('tvs_setupAlloc')
      endif

      ! loop over all body indices (still in the 'TO' family)
      ! Set the body list & start at the beginning of the list
      call obs_set_current_body_list(lobsSpaceData, index_header)
      BODY: do 
        index_body = obs_getBodyIndex(lobsSpaceData)
        if (index_body < 0) exit BODY

        if ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) == 1 ) THEN
          ICHN = nint(obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body))
          ICHN = max(0,min(ICHN,tvs_maxChannelNumber+1))
          
          ICHN = ICHN - tvs_channelOffset(nosensor)

          INDXCHN = utl_findArrayIndex(tvs_ichan(:,nosensor),tvs_nchan(nosensor),ichn)
          if ( indxchn == 0 ) then
            tvs_nchan(nosensor) = tvs_nchan(nosensor) + 1
            tvs_ichan(tvs_nchan(nosensor),nosensor) = ichn
          end if
        end if
      end do BODY
    end do HEADER

    ! Sort list of channels in ascending order.Also force at least one channel, if none are found.
    do ji = 1, tvs_nsensors
      call isort(tvs_ichan(:,ji),tvs_nchan(ji))
      if ( tvs_nchan(ji) == 0 ) then
        tvs_isReallyPresent ( ji ) =.false.
        tvs_nchan(ji) = 1
        tvs_ichan(1,ji) = 1
      end if
    end do

    write(*,*) ' tvs_setupAlloc: tvs_nobtov = ', tvs_nobtov



    !     3. Initialize TOVS radiance transfer model
    !     .  ---------------------------------------

    if ( radiativeTransferCode == 'RTTOV' ) THEN

      write(*,'(//,10x,A)') "-rttov_setup: initializing the TOVS radiative transfer model"

      allocate (tvs_coefs(tvs_nsensors)          ,stat= alloc_status(1))
      allocate (tvs_listSensors (3,tvs_nsensors) ,stat= alloc_status(2))
      allocate (tvs_opts (tvs_nsensors)          ,stat= alloc_status(3))

      call utl_checkAllocationStatus(alloc_status(1:3), " tvs_setupAlloc before rttov initialization")

      do JK=1, tvs_nsensors
        tvs_listSensors(1,JK) = tvs_platforms  (JK)
        tvs_listSensors(2,JK) = tvs_satellites (JK)
        tvs_listSensors(3,JK) = tvs_instruments(JK)

        !< General configuration options
        tvs_opts(JK) % config % apply_reg_limits = .true. ! if true application of profiles limits
        tvs_opts(JK) % config % verbose = .false. ! verbose output
        tvs_opts(JK) % config % do_checkinput = .true. ! to check if input profiles are within absolute and regression limits
        !< General RT options
        tvs_opts(JK) % rt_all % switchrad = .true.  ! to use brightness temperature (true) or radiance (false) units in AD routine
        tvs_opts(JK) % rt_all % use_q2m = .false.   ! if true use of surface humidity (false for compatibility with the way rttov 8.7 was compiled)
        tvs_opts(JK) % rt_all % addrefrac = .true. ! to account for atmospheric refraction
        !< VIS/IR RT options
        tvs_opts(JK) % rt_ir % addsolar = .false.  ! to model solar component in the near IR (2000 cm-1 et plus)
        tvs_opts(JK) % rt_ir % addaerosl = .false. ! to account for scattering due to aerosols
        tvs_opts(JK) % rt_ir % addclouds = .false. ! to account for scattering due to clouds
        tvs_opts(JK) % rt_ir % ir_sea_emis_model = 2 ! ISEM (ir_sea_emis_model 1) useful for GEORAD
                                             ! 2 would select IREMIS which is more sophisticated (to try)
        tvs_opts(JK) % rt_ir % pc % ipcreg = -1         ! index of the required PC predictors... to see later
        tvs_opts(JK) % rt_ir % pc % addpc = .false.     ! to carry out principal component calculations 
        tvs_opts(JK) % rt_ir % pc % addradrec = .false. ! to reconstruct radiances from principal components
        !< MW RT options
        tvs_opts(JK) % rt_mw % clw_data = .false.  ! profil d'eau liquide pas disponible
        tvs_opts(JK) % rt_mw % fastem_version = 6  ! use fastem version 6 microwave sea surface emissivity model (1-6)
        !< Interpolation options
        tvs_opts(JK) % interpolation % addinterp = .false. ! use of internal profile interpolator (rt calculation on model levels)
        tvs_opts(JK) % interpolation % lgradp = .false.    ! allow tl/ad of user pressure levels

        tvs_opts(JK) % rt_ir % co2_data = .false.
        tvs_opts(JK) % rt_ir % n2o_data = .false.
        tvs_opts(JK) % rt_ir % co_data  = .false.
        tvs_opts(JK) % rt_ir % ch4_data = .false.

        errorstatus = errorstatus_success

        call tmg_start(83,'RTTOV_SETUP')
        call tvs_rttov_read_coefs(errorstatus(1), tvs_coefs(jk), tvs_opts(jk), tvs_ichan(1:tvs_nchan(JK),JK), tvs_listSensors(:,JK))
        if (errorstatus(1) /= errorstatus_success) THEN
          write(*,*) 'tvs_rttov_read_coefs: fatal error reading coefficients',errorstatus,JK,tvs_listSensors(1:3,JK)
          call utl_abort('tvs_setupAlloc')
        end if
        call tmg_stop(83)

        tvs_opts(JK) % rt_ir % ozone_data = ( tvs_coefs(jk) % coef % nozone > 0 ) ! profil d'ozone disponible

! Ensure the options and coefficients are consistent
        call rttov_user_options_checkinput(errorstatus(1), tvs_opts(jk), tvs_coefs(jk))
        if (errorstatus(1) /= errorstatus_success) THEN
          write(*,*) 'error in rttov options',errorstatus
          call utl_abort('tvs_setupAlloc')
        end if
       
      end do

!    .   3.1 Validate RTTOV dimensions
!     .       -------------------------

!   Verify that all coefficient files have the same number of levels, since
!   the rest of the processing assumes this!

      do jk = 2, tvs_nsensors
        if ( tvs_coefs(jk) % coef %nlevels /= tvs_coefs(1) % coef % nlevels ) then
          write(*,'(A)') ' Number of levels not identical in all coef files'
          call utl_abort('tvs_setupAlloc')
        end if
      end do

    end if

!-----------------------------------------------------------------------


!     2. Memory allocation for radiative tranfer model variables
!     .  -----------------------------------------------------

!___ profiles

    allocate(tvs_profiles(tvs_nobtov) , stat=alloc_status(1) )
    call utl_checkAllocationStatus(alloc_status(1:1), " tvs_setupAlloc tvs_profiles 1")
 

    asw = 1 ! to allocate
    do jo = 1, tvs_nobtov
      isens = tvs_lsensor(jo)
      nl = tvs_coefs(isens) % coef % nlevels
     ! allocate model profiles atmospheric arrays with RTTOV levels dimension
      call rttov_alloc_prof(errorstatus(1),1,tvs_profiles(jo),nl, &
           tvs_opts(isens),asw,coefs=tvs_coefs(isens),init=.false. )

      call utl_checkAllocationStatus(errorstatus(1:1), " tvs_setupAlloc tvs_profiles 2")
     
    end do

!___ radiance by profile

    allocate( tvs_radiance(tvs_nobtov) ,stat= alloc_status(1))

    call utl_checkAllocationStatus(alloc_status(1:1), " tvs_setupAlloc radiances 1")
  
    do jo = 1, tvs_nobtov
      isens = tvs_lsensor(jo)
      nc = tvs_nchan(isens)
      nl = tvs_coefs(isens) % coef % nlevels
      ! allocate BT equivalent to total direct, tl and ad radiance output
      allocate( tvs_radiance(jo)  % bt  ( nc ) ,stat= alloc_status(1))
     
      tvs_radiance(jo)  % bt  ( : ) = 0.d0
    
      ! allocate clear sky radiance/BT output
      allocate( tvs_radiance(jo)  % clear  ( nc ) ,stat= alloc_status(2) )
      tvs_radiance(jo)  % clear  ( : ) = 0.d0

      call utl_checkAllocationStatus(alloc_status(1:2), " tvs_setupAlloc radiances 2")
     
    end do


    write(*,*) "Leaving tvs_setupAlloc"

  end subroutine tvs_setupAlloc

!--------------------------------------------------------------------------
!!  s/r tvs_setup : Initialisation of the TOVS processing and radiative
!!     .        transfer model.
!!    -------------------
!!     Purpose: to read namelist NAMTOV, initialize the observation error covariance
!!              and setup RTTOV-12.
!--------------------------------------------------------------------------

  subroutine tvs_setup
    implicit none
    integer  JK, IERR, nulnam
    integer ,external :: fclos, fnom
    integer :: nsensors
    character(len=15) :: csatid(tvs_maxNumberOfSensors), cinstrumentid(tvs_maxNumberOfSensors)
    character(len=8)  :: crtmodl
    logical :: ldbgtov

    namelist /NAMTOV/ nsensors, csatid, cinstrumentid
    namelist /NAMTOV/ ldbgtov
    namelist /NAMTOV/ useUofWIREmiss, crtmodl

 
  !     .  1.1 Default values for namelist variables
  !     .      -------------------------------------

    nsensors = 0
    csatid(:) = '***UNDEFINED***'
    cinstrumentid(:) = '***UNDEFINED***'
    csatid(1) = 'NOAA16'
    cinstrumentid(1) = 'AMSUA'
    ldbgtov = .FALSE.
    crtmodl = 'RTTOV'
    useUofWIREmiss = .false.

  !     .  1.2 Read the NAMELIST NAMTOV to modify them
  !     .      ---------------------------------------
 
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam, nml=namtov, iostat=ierr)
    if (ierr /= 0) call utl_abort('tvs_setup: Error reading namelist')
    if (mpi_myid == 0) write(*,nml=namtov)
    ierr = fclos(nulnam)


    !   .  1.3 Transfer namelist variables to module variables
    !   .      -- ---------------------------------------------
 
    tvs_nsensors = nsensors
    tvs_debug = ldbgtov
    radiativeTransferCode = crtmodl
    tvs_instrumentName(:) = cinstrumentid(:)
    tvs_satelliteName(:) = csatid(:)

    !   .  1.4 Validate namelist values
    !   .      ------------------------
    
    if ( tvs_nsensors == 0 ) then
      if(mpi_myid==0) then 
        write(*,*) ' ====================================================='
        write(*,*) ' tvs_setup: Number of sensors is zero, skipping setup'
        write(*,*) ' ====================================================='
      end if
      return
    end if

    if ( radiativeTransferCode /= 'RTTOV' ) THEN
      write(*,'(A)') ' Invalid radiation model name'
      call utl_abort('tvs_setup')
    end if

    if ( tvs_nsensors > tvs_maxNumberOfSensors ) THEN
      write(*,'(A)') ' Number of sensors (tvs_nsensors) is greater than maximum allowed (tvs_maxNumberOfSensors)'
      call utl_abort('tvs_setup')
    end if

    !   .  1.5 Print the content of this NAMELIST
    !   .      ----------------------------------

    if(mpi_myid == 0) then
      write(*,'(A)') 
      write(*,'(3X,A)') '- Parameters used for TOVS processing (read in NAMTOV)'
      write(*,'(3X,A)') '  ----------------------------------------------------'
      write(*,'(6X,A,2X,L1)') 'TOVS debug                           : ', tvs_debug
      write(*,'(6X,A,2X,L1)') 'Use of UW IR land emissivity atlases : ', useUofWIREmiss
      write(*,'(6X,A,2X,A)')  'Radiative transfer model             : ', radiativeTransferCode
      write(*,'(6X,A,2X,I3)') 'Number of sensors                    : ', tvs_nsensors
      write(*,'(6X,"Satellite ids          : ",10A10)') (tvs_satelliteName(JK), JK=1,tvs_nsensors)
      write(*,'(6X,"Instrument ids         : ",10A10)') (tvs_instrumentName(JK), JK=1,tvs_nsensors)
      write(*,'(A)') 
      write(*,'(A)') 
      write(*,'(A)') 
      write(*,'(3X,A)') '- Reading and initialization in preparation to the TOVS processing'
      write(*,'(5X,A)') '----------------------------------------------------------------'
    end if

    !   .  1.6 Set up platform, satellite, instrument and channel mapping
    !   .      ----------------------------------------------------------

    call SENSORS

  end subroutine tvs_setup


!--------------------------------------------------------------------------
!!**s/r SENSORS : Initialisation of the RTTOV-10 platform, satellite
!!                and instrument ID's. Also set burp to RTTOV channel
!!                mapping offset.
!!
!!*    Purpose: to verify and transfom the sensor information contained in the
!!              NAMTOV namelist into the variables required by RTTTOV-7:
!!              platform, satellite and instrument ID's.
!!
!!Variables:
!!     i : tvs_nsensors      : number of sensors
!!     i : tvs_satelliteName        : satellite ID (e.g. 'NOAA15')
!!     i : tvs_instrumentName : instrument ID (e.g. 'AMSUA')
!!     o : tvs_platforms      : RTTOV platform ID numbers (e.g. 1 for  NOAA)
!!     o : SATELLITE     : RTTOV satellite ID numbers (e.g. 15)
!!     o : INSTRUMENT    : RTTOV instrument ID numbers (e.g. 3 for AMSUA)
!!     o : tvs_channelOffset    : BURP to RTTOV channel mapping offset
!--------------------------------------------------------------------------

  subroutine SENSORS
    implicit none

    integer J, K, IPOS1, IPOS2
    integer NUMEROSAT, IERR, KINDEX, nulnam
    character(len=15) :: TEMPOCSATID
    logical, save :: LFIRST=.true.
    integer, save :: IOFFSET1B(0:ninst-1)
    character(len=8) :: LISTINSTRUM(0:ninst-1)
    character(len=15) :: tempo_inst
    integer:: LISTOFFSET(0:ninst-1)
    namelist /NAMCHANOFFSET/ LISTOFFSET, LISTINSTRUM
    integer ,external :: fnom, fclos

!
!*    .  1.0 Go through sensors and set RTTOV-10 variables
!     .      --------------------------------------------

    do J=1, tvs_nsensors
      tvs_platforms  (J) = -1
      tvs_satellites (J) = -1
      tvs_instruments(J) = -1
      tvs_channelOffset(J) = -1
    end do

    if (LFIRST) then
     ! read the namelist
      NULNAM = 0
      IERR = FNOM(NULNAM,'./flnml','FTN+SEQ+R/O',0)
      if (IERR /= 0) then
        write(*,*) "Error while opening namelist file !"
        call utl_abort("sensors")
      end if
      LISTOFFSET(:) = 0
      LISTINSTRUM(:) = "XXXXXXXX"
      read(NULNAM,NAMCHANOFFSET, iostat=ierr)
      if (IERR /= 0) then
        write(*,*) "Error while reading namelist file !"
        call utl_abort("sensors")
      end if
      do j=0, ninst - 1
        if ( LISTINSTRUM(j) /= "XXXXXXXX" ) then
          IOFFSET1B( tvs_getInstrumentId( LISTINSTRUM(j) ) )  = LISTOFFSET(j)
        end if
      end do
      ierr = FCLOS(NULNAM)
      LFIRST = .false.
    end if

!*    .  1.1 Set platforms and satellites
!     .      ----------------------------
!
!** N.B.: Special cases for satellites TERRA and AQUA.
!**       For consistency with the RTTOV-10 nomenclature, rename:
!**       TERRA  to  EOS1
!**       AQUA   to  EOS2
!**       NPP    to  NPP0
    do J = 1, tvs_nsensors
      if    ( tvs_satelliteName(J) == 'TERRA' ) THEN
        TEMPOCSATID = 'eos1'
      else if ( tvs_satelliteName(J) == 'AQUA'  ) THEN
        TEMPOCSATID = 'eos2'
      else if ( tvs_satelliteName(J) == 'NPP'  ) THEN
        TEMPOCSATID = 'jpss0'
      else if ( tvs_satelliteName(J) == 'JPSS'  ) THEN
        TEMPOCSATID = 'jpss0'
      else if ( tvs_satelliteName(J)(1:6) == 'HMWARI'  ) THEN
        TEMPOCSATID = 'himawari' // trim(tvs_satelliteName(J) (7:15))
      else
        call up2low(tvs_satelliteName(J),TEMPOCSATID)
      end if
      do K = 1, nplatforms
        IPOS1 = LEN_TRIM(platform_name(K))
        IPOS2 = INDEX(TEMPOCSATID,platform_name(K)(1:IPOS1))
        if ( IPOS2 == 1 ) THEN
          tvs_platforms(J) = K
          KINDEX = K
        end if
      end do
      if ( tvs_platforms(J) < 0 ) THEN
        write(*,'(A)') ' Satellite ' // trim(TEMPOCSATID) // ' not supported.'
        call utl_abort('SENSORS')
      else
        IPOS1 = LEN_TRIM(platform_name(KINDEX))
        IPOS2 = LEN_TRIM(TEMPOCSATID)
        read(TEMPOCSATID(IPOS1+1:IPOS2),*,IOSTAT=IERR) NUMEROSAT
        NUMEROSAT = ABS ( NUMEROSAT )
        if ( IERR /= 0) THEN
          write(*,'(A,1x,i6,1x,i3,1x,i3,1x,A15)') "Problem while reading satellite number", &
               IERR, ipos1, ipos2, TEMPOCSATID
          call utl_abort('SENSORS')
        else
          tvs_satellites(J) = NUMEROSAT
        end if
      end if
    end do

!*    .  1.2 Set instruments,
!     .      also set channel offset, which is in fact a channel mapping between
!     .      the channel number in BURP files and the channel number used in
!     .      RTTOV-10.
!     .      --------------------------------------------------------------------

    do J = 1, tvs_nsensors
      if ( tvs_satelliteName(J)(1:4) == "GOES") then !cas particulier
        tvs_instruments(J) = inst_id_goesim
      else if ( tvs_satelliteName(J)(1:5) == "MTSAT") then !autre cas particulier
        tvs_instruments(J) = inst_id_gmsim
      else 
        call up2low(tvs_instrumentName(J),tempo_inst)
        do K = 0, ninst -1 
          if ( trim(tempo_inst) == trim(inst_name(K))) THEN
            tvs_instruments(J) = K
          end if
        end do
      end if
      if ( tvs_instruments(J) < 0 ) THEN
        write(*,'(A)') ' INSTRUMENT '// trim( tvs_instrumentName(J)) // ' not supported.'
        call utl_abort('SENSORS')
      end if
      tvs_channelOffset(J) = IOFFSET1B(tvs_instruments(J))
    end do

!*    .   1.3 Print the RTTOV-10 related variables
!     .       -----------------------------------

    if (mpi_myid == 0) then
      write(*,*)
      write(*,'(3X,A)') '- SENSORS. Variables prepared for RTTOV-12:'
      write(*,'(3X,A)') '  ----------------------------------------'
      write(*,*)
      write(*,'(6X,A,I3)')   "Number of sensors       : ", tvs_nsensors
      write(*,'("Platform numbers        : ",6X,10I3)')  (tvs_platforms(J), J=1,tvs_nsensors)
      write(*,'("Satellite numbers       : ",6X,10I3)')  (tvs_satellites(J), J=1,tvs_nsensors)
      write(*,'("Instrument numbers      : ",6X,10I3)')  (tvs_instruments(J), J=1,tvs_nsensors)
      write(*,'("Channel mapping offsets : ",6X,10I3)')  (tvs_channelOffset(J), J=1,tvs_nsensors)
    end if

  end subroutine sensors


!--------------------------------------------------------------------------
!! Function to check if the given idatyp (a.k.a. codtyp) corresponds
!!  to a radiance
!--------------------------------------------------------------------------

  logical function tvs_isIdBurpTovs(idatyp)
    implicit none
    integer ,intent(in) :: idatyp
!*********************************************
    logical ,save :: lfirst=.true.
    integer,save :: ninst_tovs
    integer :: nulnam, ierr, i 
    integer,external :: fnom, fclos
    integer, save :: list_inst(ninst)
    character (len=22) :: inst_names(ninst)
    namelist /namtovsinst/ inst_names

    if (tvs_nsensors == 0) then
    ! no tovs data will be read, therefore false
      tvs_isIdBurpTovs = .false.
      return
    end if

    if (lfirst) then
      nulnam = 0
      ninst_tovs = 0
      list_inst(:) = -1
      inst_names(:) = "XXXXXX"
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam, nml=namtovsinst, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isIdBurpTovs: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=namtovsinst)
      ierr = fclos(nulnam)
      do i=1, ninst
        if (inst_names(i) == "XXXXXX" ) then
          ninst_tovs= i - 1
          exit
        else
          list_inst(i) = codtyp_get_codtyp( inst_names(i) )
          if (list_inst(i) < 0) then
            write(*,*) inst_names(i)
            call utl_abort('tvs_isIdBurpTovs: unknown instrument in namtovsinst namelist')
          end if
        end if
      end do
      if ( ninst_tovs == 0 ) call utl_abort('tvs_isIdBurpTovs: Empty namtovsinst namelist')
      lfirst = .false.
    end if

    tvs_isIdBurpTovs = .false.

    do i=1, ninst_tovs
      if (idatyp == list_inst(i) ) then
        tvs_isIdBurpTovs = .true.
        exit
      end if
    end do

  end function tvs_isIdBurpTovs

!--------------------------------------------------------------------------
!! function to check if the provided idburp (a.k.a. codtyp)
!! corresponds to instrument cinst 
!--------------------------------------------------------------------------

  logical function tvs_isIdBurpInst(idburp,cinst)
    implicit none
    integer ,intent(in) :: idburp
    character (len=*) ,intent(in) :: cinst

    if (tvs_nsensors == 0) then
      ! no tovs data will be read, therefore false
      tvs_isIdBurpInst = .false.
      return
    end if

    tvs_isIdBurpInst = ( idburp == codtyp_get_codtyp(cinst) )

  end function tvs_isIdBurpInst

!--------------------------------------------------------------------------
!! S. Heilliette November 2015
!! return RTTOV platform id (>0) from
!! platform name.
!! -1 if not found
!--------------------------------------------------------------------------

  integer function tvs_getPlatformId(name)
    implicit none
    character(len=*),intent(in) :: name
!************************************
    integer :: i, ipos1, ipos2
    character(len=64) :: tempo_name
!*************************************

    tvs_getPlatformId = -1
    ipos1 = LEN_TRIM(name)
    call up2low(name(1:ipos1),tempo_name(1:ipos1))

    if ( INDEX(tempo_name(1:ipos1),"npp") /= 0 ) then
      tvs_getPlatformId = platform_id_jpss
    else if ( INDEX(tempo_name(1:ipos1),"hmwari") /= 0 ) then
      tvs_getPlatformId = platform_id_himawari
    else
      do i = 1, nplatforms
        ipos2 = INDEX(tempo_name(1:ipos1),trim(platform_name(i)))
        if ( ipos2 == 1) then
          tvs_getPlatformId = i
          exit
        end if
      end do
    end if

  end function tvs_getPlatformId

!--------------------------------------------------------------------------
!! S. Heilliette November 2015
!! return RTTOV instrument id from
!! intrument name. 0 is a valid answer.
!! -1 if not found
!--------------------------------------------------------------------------

  integer function tvs_getInstrumentId(name)
    implicit none
    character(len=*),intent(in) :: name
!************************************
    integer :: i, length
    character(len=64) :: tempo_name
!***********************************
    tvs_getInstrumentId = -1
    length = LEN_TRIM(name)
    call up2low(name(1:length),tempo_name(1:length))
    if ( trim(tempo_name(1:length)) == "goesim" ) then
      tvs_getInstrumentId = inst_id_goesim
    else if ( trim(tempo_name(1:length)) == "gmsim" ) then
      tvs_getInstrumentId = inst_id_gmsim
    else if ( trim(tempo_name(1:length)) == "mtsatim" ) then
      tvs_getInstrumentId = inst_id_mtsatim
    else
      do i = 0, ninst - 1
        if (trim(inst_name(i)) == tempo_name(1:length) ) then
          tvs_getInstrumentId = i
          exit
        end if
      end do
    end if
  end function tvs_getInstrumentId


!--------------------------------------------------------------------------
!! S. Heilliette November 2015
!! given an RTTOV instrument code
!! return if it is an hyperspectral one
!! information from namelist NAMHYPER
!--------------------------------------------------------------------------

  logical function tvs_isInstrumHyperSpectral(instrum)
    implicit none
    integer,intent(in) :: instrum
!******************************************
    integer ,parameter :: maxsize = 100
    integer :: nulnam, ierr, i 
    integer ,save :: list_inst(maxsize), ninst_hir
    logical, save :: lfirst = .true.
    integer ,external :: fclos, fnom
    character (len=6) :: name_inst(maxsize)
    namelist /NAMHYPER/ name_inst

    if (lfirst) then
      nulnam = 0
      ninst_hir = 0
      name_inst(:) = "XXXXXX"
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namhyper, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isInstrumHyperSpectral: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=namhyper)
      ierr = fclos(nulnam)
      list_inst(:) = -1
      do i=1, maxsize
        list_inst(i) = tvs_getInstrumentId( name_inst(i) )
        if (name_inst(i) /= "XXXXXX") then
          if (list_inst(i) == -1) then
            write(*,*) i,name_inst(i)
            call utl_abort('tvs_isInstrumHyperSpectral: Unknown instrument name')
          end if
        else
          ninst_hir = i -1
          exit
        end if
      end do
      lfirst = .false.
      if (ninst_hir == 0) then
        write(*,*) "tvs_isInstrumHyperSpectral: Warning : empty namhyper namelist !"
      end if
    end if
    tvs_isInstrumHyperSpectral = .false.
    do i=1, ninst_hir
      if ( instrum == list_inst(i)) then
        tvs_isInstrumHyperSpectral = .true.
        exit
      end if
    end do

  end function tvs_isInstrumHyperSpectral


!--------------------------------------------------------------------------
!! S. Heilliette November 2015
!! given an RTTOV instrument code
!! return if it is a Geostationnary Imager
!! information from namelist NAMGEO
!--------------------------------------------------------------------------

  logical function tvs_isInstrumGeostationary(instrum)
    implicit none
    integer,intent(in) :: instrum
!******************************************
    integer ,parameter :: maxsize = 100
    integer :: nulnam, ierr, i 
    integer ,save :: list_inst(maxsize), ninst_geo
    logical, save :: lfirst = .true.
    character (len=8) :: name_inst(maxsize)
    integer ,external :: fnom, fclos

    namelist /NAMGEO/ name_inst
    if (lfirst) then
      nulnam = 0
      ninst_geo = 0
      name_inst(:) = "XXXXXX"
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namgeo, iostat=ierr)
      if (ierr /= 0) call utl_abort('tvs_isInstrumGeostationary: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=namgeo)
      ierr = fclos(nulnam)
      list_inst(:) = -1
      do i=1, maxsize
        list_inst(i) = tvs_getInstrumentId( name_inst(i) )
        if (name_inst(i) /= "XXXXXX") then
          if (list_inst(i) == -1) then
            write(*,*) i,name_inst(i)
            call utl_abort('tvs_isInstrumGeostationary: Unknown instrument name')
          end if
        else
          ninst_geo = i - 1
          exit
        end if
      end do
      lfirst = .false.
      if (ninst_geo == 0) then
        write(*,*) "tvs_isInstrumGeostationary: Warning : empty namgeo namelist !"
      end if
    end if
    tvs_isInstrumGeostationary = .false.
    do i=1, ninst_geo
      if ( instrum == list_inst(i)) then
        tvs_isInstrumGeostationary = .true.
        exit
      end if
    end do
    
  end function tvs_isInstrumGeostationary

!--------------------------------------------------------------------------
!!**s/r tvs_mapInstrum : Map burp satellite instrument (element #2019)
!!                    to RTTOV-10 instrument.
!!
!!Author  : J. Halle *CMDA/SMC May 2002
!!
!!Revision 001: J. Halle *CMDA/SMC Sept. 2005
!!              - add MHS.
!!         002: R Sarrazin CMDA, April 2008
!!              - comment on MTSAT imager instrument number
!!         003: R. McTaggart-Cowan *RPN  Mar 2012
!!              - Use assumed-length declarations for string dummy args
!!         003: S. Macpherson ARMA, August 2011
!!              - add NPP/ATMS codtyp=192
!!         004: C. Cote Juillet 2015
!!              - Ajout Himawari-8/AHI codtyp=185
!!         005: S. Heilliette Novembre 2015
!!              - suppression des tables hardcodees
!!                lecture depuis une namelist
!!
!!
!!*    Purpose:  Map burp satellite instrument (element #2019)
!!               to RTTOV-7 instrument.
!!*              A negative value is returned, if no match in found.
!!
!!               Table of  RTTOV-7 instrument identifier
!                ---------------------------------------
!!               Instrument        Instrument identifier  Sensor type
!                ---------         ---------------------  -----------
!!               HIRS               0                     ir
!!               MSU                1                     mw
!!               SSU                2                     ir
!!               AMSUA              3                     mw
!!               AMSUB              4                     mw
!!               AVHRR              5                     ir
!!               SSMI               6                     mw
!!               VTPR1              7                     ir
!!               VTPR2              8                     ir
!!               TMI                9                     mw
!!               SSMIS             10                     mw
!!               AIRS              11                     ir
!!               MODIS             13                     ir
!!               ATSR              14                     ir
!!               MHS               15                     mw
!!               ATMS              19                     mw
!!               MVIRI             20                     ir
!!               SEVIRI            21                     ir
!!               GOESIMAGER        22                     ir
!!               GOESSOUNDER       23                     ir
!!               GMS/MTSAT IMAGER  24                     ir
!!               FY2-VISSR         25                     ir
!!               FY1-MVISR         26                     ir
!!               AHI               56                     ir
!!
!!Arguments:
!!     i : INSTRUMBURP : burp satellite instrument (element #2019)
!!     o : INSTRUM     : RTTOV-7 instrument ID numbers (e.g. 3 for  AMSUA)
!--------------------------------------------------------------------------

  subroutine tvs_mapInstrum(INSTRUMBURP,INSTRUM)
    implicit none
    integer J,INSTRUMBURP,INSTRUM,numinstburp
    integer,parameter :: MXINSTRUMBURP   = 100
    integer,save ::   LISTBURP(MXINSTRUMBURP)
    character(len=8),save :: LISTINSTRUM(MXINSTRUMBURP)
    namelist /NAMINST/ LISTBURP, LISTINSTRUM
    logical,save :: LFIRST = .true.
    integer :: NULNAM, IER
    integer ,external :: fnom, fclos
!*************
!      TESTING ONLY: burp 2047 is AIRS.
!*************
!
!*            Table of BURP satellite sensor identifier element #002019
!*            ----------------------------------------------------------

!
!
!*    .  1.0 Find instrument
!     .      ---------------------------

    if (LFIRST) THEN
! set the default values
      LISTBURP(:) = -1
      LISTINSTRUM(:) = "XXXXXXXX"

! read the namelist
      NULNAM = 0
      IER = FNOM(NULNAM,'./flnml','FTN+SEQ+R/O',0)
      if (IER /= 0) then
        write(*,*) "Error while opening namelist file !"
        call utl_abort("tvs_mapInstrum")
      end if
      read(NULNAM,NAMINST,iostat=ier)
      if (IER /= 0) then
        write(*,*) "Error while reading namelist file !"
        call utl_abort("tvs_mapInstrum")
      end if
      ier = FCLOS(NULNAM)

! figure out how many valid elements in the lists
      do J=1, MXINSTRUMBURP
        if(listburp(j) == -1) then
          numinstburp = j - 1
          exit
        end if
      end do
      if (numinstburp > MXINSTRUMBURP) then
        call utl_abort('tvs_mapInstrum: exceeded maximum number of platforms')
      end if
      write(*,*) 'tvs_mapInstrum: number of satellites found in namelist = ',numinstburp
      write(*,*) 'tvs_mapInstrum: listburp   = ',listburp(1:numinstburp)
      write(*,*) 'tvs_mapInstrum: listinstrum    = ',listinstrum(1:numinstburp)
      LFIRST=.FALSE.
    end if

    INSTRUM = -1
    do J=1, MXINSTRUMBURP
      if ( INSTRUMBURP == LISTBURP(J) ) THEN
        INSTRUM = tvs_getInstrumentId( LISTINSTRUM(J) )
        EXIT
      end if
    end do

  end subroutine tvs_mapInstrum

!--------------------------------------------------------------------------
!!**s/r tvs_mapSat : Map burp satellite identifier (element #1007)
!!                to RTTOV-10 platform and satellite.
!!
!!Author  :       J. Halle *CMDA/SMC May 2002
!!
!!Revision 001  : J. Halle *CMDA/AES Jul 2005
!!                . add NOAA-18.
!!
!!Revision 002  : J. Halle *CMDA/AES May 2007
!!                . add METOP 1,2,3.
!!
!!Revision 003  : R. Sarrazin CMDA   Apr 2008
!!                  add MTSAT1, GOES13 and MSG2, modif to MSG1
!!
!!Revision 004  : C. Cote  mars 2009
!!                . add NOAA-19.
!!
!!Revision 005  : S. Macpherson *ARMA  Jul 2010
!!                . add SSMIS satellites DMSP17-18
!!
!!Revision 006  : S. Macpherson *ARMA  Feb 2013
!!                . add NPP/ATMS codtyp=192
!!
!!Revision 007  : J. Morneau CMDA    FEb 2014
!!                . add GOES15 and MTSAT-2
!!                . add GOES14
!!                . add MSG3 and MSG4
!!Revision 008 :  S. Heilliette    March 2014
!!                . Major Modification:
!!                . information is now read from namelist
!!                . instead of being hardcoded
!!
!!    Purpose:  Map burp satellite identifier (element #1007)
!!               to RTTOV-7 platform and satellite.
!!               Negative values are returned, if no match in found.
!!
!                ---------------------------------------------
!!               Table of  RTTOV-7 platform identifier
!                ---------------------------------------------
!!               Platform          RTTOV-7 platform identifier
!                ---------         ---------------------------
!!               NOAA               1
!!               DMSP               2
!!               METEOSAT           3
!!               GOES               4
!!               GMS                5
!!               FY2                6
!!               TRMM               7
!!               ERS                8
!!               EOS                9
!!               METOP             10
!!               ENVISAT           11
!!               MSG               12
!!               FY1               13
!!               ADEOS             14
!!               MTSAT             15
!!               CORIOLIS          16
!!               NPP               17
!!               ---------------------------------------------
!!
!!               Example: NOAA15, which has a burp satellite identifier value of 206,
!!                        is mapped into the following:
!!                                RTTOV-7 platform  =  1,
!!                                RTTOV-7 satellite = 15.
!!
!!
!!
!!Arguments:
!!     i : ISATBURP      : BURP satellite identifier
!!     o : IPLATFORM     : RTTOV-7 platform ID numbers (e.g. 1 for  NOAA)
!!     o : ISAT          : RTTOV-7 satellite ID numbers (e.g. 15)
!--------------------------------------------------------------------------

  subroutine tvs_mapSat(ISATBURP,IPLATFORM,ISAT)
    implicit none
    
    integer J,ISATBURP,IPLATFORM,ISAT
    integer IER,NULNAM
    logical ,save :: LFIRST=.TRUE.
    integer ,external :: FNOM, FCLOS

    integer, parameter :: mxsatburp = 100
    integer,save :: NUMSATBURP

! Table of BURP satellite identifier element #001007
    integer,save :: LISTBURP(mxsatburp)
! Table of RTTOV platform identifier
    character(len=8),save :: LISTPLAT(mxsatburp)
! Table of RTTOV satellite identifier
    integer,save :: LISTSAT (mxsatburp)

    namelist /NAMSAT/ LISTBURP, LISTPLAT, LISTSAT

!     Fill tables from namelist at the first call 
!     -------------------------------------------
!
    if (LFIRST) THEN
! set the default values
      LISTBURP(:) = -1
      LISTSAT(:) = -1
      LISTPLAT(:) = "XXXXXXXX"

! read the namelist
      NULNAM = 0
      IER = FNOM(NULNAM,'./flnml','FTN+SEQ+R/O',0)
      if (IER /= 0) then
        write(*,*) "Error while opening namelist file !"
        call utl_abort("tvs_mapSat")
      end if
      read(NULNAM, NAMSAT, iostat = ier)
      if (IER /= 0) then
        write(*,*) "Error while reading namelist file !"
        call utl_abort("tvs_mapSat")
      end if
      ier = FCLOS(NULNAM)

! figure out how many valid elements in the lists
      do J=1, MXSATBURP
        if(listburp(j) == -1) then
          numsatburp = j - 1
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
      LFIRST = .FALSE.
    end if

!     .  Find platform and satellite
!     .      ---------------------------
!
    IPLATFORM = -1
    ISAT      = -1
    do J=1, NUMSATBURP
      if ( ISATBURP == LISTBURP(J) ) THEN
        IPLATFORM = tvs_getPlatformId( LISTPLAT(J) )
        ISAT = LISTSAT (J)
        EXIT
      end if
    end do

  end subroutine tvs_mapSat

!--------------------------------------------------------------------------
!! subroutine to initialize the chanprof structure used by RTTOV
!--------------------------------------------------------------------------

  subroutine TVS_getChanprof(sensor_id, iptobs, nprofiles, ObsSpaceData, chanprof, iptobs_cma,assim_flag_val)
    implicit none
    integer ,intent(in) :: nprofiles, sensor_id, iptobs(:)
    integer ,intent(out),optional :: iptobs_cma(:)
    integer ,intent(in) ,optional :: assim_flag_val
    type(struct_obs) :: ObsSpaceData
    type(rttov_chanprof) :: chanprof(:)
!******************************************************************************
    integer :: count, profile_index, header_index, istart, iend, body_index, ichn, nrank, iobs, assim_flag_val2
!******************************************************************************
! Build the list of channels/profiles indices

    if (present(assim_flag_val)) then
      assim_flag_val2 = assim_flag_val
    else
      assim_flag_val2 = 1
    end if

    count = 0
         
    do profile_index = 1,  nprofiles
      iobs = iptobs(profile_index)
      header_index = tvs_lobsno(iobs)
      istart = obs_headElem_i(ObsSpaceData,OBS_RLN,header_index)
      iend= obs_headElem_i(ObsSpaceData,OBS_NLV,header_index) + istart - 1
      do body_index = istart, iend
        if (obs_bodyElem_i(ObsSpaceData,OBS_ASS,body_index) == assim_flag_val2) then
          ichn = nint(obs_bodyElem_r(ObsSpaceData,OBS_PPP,body_index))
          ichn = max(0, min(ichn,tvs_maxChannelNumber + 1))
          ICHN = ICHN - tvs_channelOffset(sensor_id)
          do nrank = 1, tvs_nchan(sensor_id)
            if ( ichn == tvs_ichan(nrank,sensor_id) ) exit
          end do
          if (nrank /= tvs_nchan(sensor_id)+1) then
            count =  count + 1
            chanprof(count)%prof = profile_index
            chanprof(count)%chan = nrank
            if (present(iptobs_cma)) iptobs_cma(count) = body_index
          else
            write(*,*) "strange channel number",ichn
          end if
        end if
      end do
    end do
  
  end subroutine TVS_getChanprof


!--------------------------------------------------------------------------
!! subroutine to count radiances selected for assimilation
!--------------------------------------------------------------------------

  integer function tvs_countRadiances(iptobs, nprofiles, ObsSpaceData,assim_flag_val)
    implicit none
    integer ,intent(in) :: nprofiles, iptobs(:)
    integer ,intent(in),optional :: assim_flag_val
    type(struct_obs) :: ObsSpaceData
!******************************************************************************
    integer :: profile_index, header_index, istart, iend, body_index, iobs, assim_flag_val2

    if (present(assim_flag_val)) then
      assim_flag_val2 = assim_flag_val
    else
      assim_flag_val2 = 1
    end if

    tvs_countRadiances = 0
    do profile_index = 1, nprofiles
      iobs = iptobs(profile_index)
      header_index = tvs_lobsno(iobs)
      istart = obs_headElem_i(ObsSpaceData,OBS_RLN,header_index)
      iend = obs_headElem_i(ObsSpaceData,OBS_NLV,header_index) + istart - 1
      do body_index = istart, iend
        if(obs_bodyElem_i(ObsSpaceData,OBS_ASS,body_index) == assim_flag_val2) tvs_countRadiances  = tvs_countRadiances + 1
      end do
    end do

  end function tvs_countRadiances

!--------------------------------------------------------------------------
!! subroutine to get emissivity for Hyperspectral Infrared Sounders (AIRS, IASI, CrIS, ...)
!--------------------------------------------------------------------------
  subroutine tvs_getHIREmissivities(sensor_id, iptobs, nprofiles, ObsSpaceData, surfem,assim_flag_val)

    implicit none
    integer ,intent(in) :: nprofiles, sensor_id, iptobs(:)
    integer ,intent(in),optional :: assim_flag_val
    type(struct_obs) :: ObsSpaceData
    real(8), intent(out) :: surfem(:)
!***************************************************************************
    integer :: count, profile_index, iobs, istart, iend, index_body, index_header, assim_flag_val2
!***************************************************************************

    if (present(assim_flag_val)) then
      assim_flag_val2 = assim_flag_val
    else
      assim_flag_val2 = 1
    end if

    count = 0 
    surfem(:) = 0.98d0
    do profile_index = 1, nprofiles
      iobs = iptobs(profile_index)
      index_header = tvs_lobsno(iobs)
      istart = obs_headElem_i(ObsSpaceData,OBS_RLN,index_header)
      iend = obs_headElem_i(ObsSpaceData,OBS_NLV,index_header) + istart - 1
      do index_body = istart, iend
        if(obs_bodyElem_i(ObsSpaceData,OBS_ASS,index_body) == assim_flag_val2) then
          count = count + 1
          surfem ( count ) = obs_bodyElem_r(ObsSpaceData,OBS_SEM,index_body)
        end if
      end do
    end do

  end subroutine tvs_getHIREmissivities


!--------------------------------------------------------------------------
!! subroutine to get emissivity for microwave souders ans infrared geostationary imagers 
!--------------------------------------------------------------------------

  subroutine tvs_getOtherEmissivities(chanprof, iptobs, nradiances, sensor_type, instrument, surfem, calcemis)
    implicit none
    integer ,intent(in) :: nradiances, iptobs(:),sensor_type, instrument
    real(8), intent(out) :: surfem(:)
    logical, intent(out) :: calcemis(:)
    type(rttov_chanprof),intent(in) :: chanprof(:)
!*******************************************************************************************************
    integer :: radiance_index, profile_index, iobs, surface_type
!*******************************************************************************************************
    do radiance_index = 1, nradiances
      profile_index = chanprof(radiance_index)%prof
      iobs = iptobs(profile_index)
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
!! Subroutine to fill in tvs_profiles structure before call to non-linear, tangent-linear or
!! adjoint of RTTOV
!--------------------------------------------------------------------------

  subroutine tvs_fillProfiles(columnghr,lobsSpaceData,datestamp,LIMLVHU,bgckMode,beSilent)
    implicit none
    type(struct_obs),intent(in) :: lobsSpaceData
    type(struct_columnData),intent(in) :: columnghr
    integer,intent(in) :: datestamp
    real(8),intent(in)  :: LIMLVHU
    logical,intent(in) :: bgckMode
    logical,intent(in) :: beSilent

    logical :: diagTtop,TopAt10hPa

    integer :: ksurf, jpmotop, jpmolev
    integer :: isatzen 
    integer :: isatazim, instrum, iplatform
    integer :: isunazim
    integer :: isunza
    integer :: nlevels,nobmax
    integer :: j, i, sensor_id, iobs, jj
    integer :: count_profile, header_index
    integer :: jn, jl
    integer :: ilowlvl_M,ilowlvl_T,nlv_M,nlv_T
    integer :: status, Vcode
    integer :: ier,day,month,year,ijour,itime
    integer :: alloc_status(14)
    
    integer,external ::  omp_get_num_threads
    integer,external ::  newdate

    integer, allocatable :: iptobs    (:) 
    integer, allocatable :: iptobscma (:) 
  
    type(struct_vco), pointer :: vco

    real(8), allocatable :: to    (:,:)
    real(8), allocatable :: lqo   (:,:)
    real(8), allocatable :: zho   (:,:)
    real(8), allocatable :: toext (:,:)
    real(8), allocatable :: qoext (:,:)
    real(8), allocatable :: zvlev (:,:)
    real(8), allocatable :: zt    (:,:)
    real(8), allocatable :: zlq   (:,:)
    real(8), allocatable :: zht   (:,:)
    real(8), allocatable :: xpres (:)
    real(8), allocatable :: zlat  (:)
    real(8), allocatable :: toto3obs(:),PP(:,:) 
    real(8), allocatable :: ozo     (:,:)
 
    real(8) :: zlon
    real(8) :: zptop, zptopmbs
 

    if ( .not. beSilent ) write(*,*) "Entering tvs_fillProfiles subroutine"
  
    if (tvs_nobtov == 0) return    ! exit if there are no tovs data

!
!     1.    Set index for model's lowest level and model top
!     .     ------------------------------------------------
    
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

    ier = newdate(datestamp,ijour,itime,-3)
    if (ier < 0) then
      write(*,*) "Invalid datestamp ",datestamp,ijour,itime,ier
      call utl_abort('tvs_fillProfiles')
    end if
    year= ijour / 10000
    MONTH = MOD(ijour / 100,100)
    DAY = MOD(ijour,100)

!     1.2   Read ozone climatology
!     .     ----------------------

    call ozo_read_climatology(datestamp)

!
!     2.  Fill profiles structure
!     .   -----------------------

! loop over all instruments
    sensor_loop: do sensor_id=1,tvs_nsensors

! first loop over all obs.
      count_profile = 0
      bobs1: do iobs = 1, tvs_nobtov
        if (tvs_lsensor(iobs) == sensor_id) then
          count_profile = count_profile + 1
          NOBMAX = iobs
        end if
      end do bobs1

      if (count_profile == 0) cycle sensor_loop

      nlevels = tvs_coefs(sensor_id) %coef% nlevels
      allocate ( xpres (nlevels) )
      xpres = tvs_coefs(1)% coef % ref_prfl_p
      jpmotop = 1
      do jl = 2, nlevels
        if ( zptopmbs >= xpres(jl - 1) .and. &
             zptopmbs < xpres(jl)        ) then 
          jpmotop = jl
          exit
        end if
      end do
      if ( .not. beSilent ) write(*,*) 'tvs_fillProfiles: jpmotop=', sensor_id, jpmotop
      jpmolev = nlevels - jpmotop + 1

      alloc_status(:) = 0
      allocate (iptobs    (count_profile)        ,stat= alloc_status(1) )
      allocate (iptobscma (count_profile)        ,stat= alloc_status(2) )
      allocate (zlat      (count_profile)        ,stat= alloc_status(3) )
      allocate (ozo       (nlevels,count_profile),stat= alloc_status(4)) 
      allocate (to        (jpmolev,count_profile),stat= alloc_status(5))
      allocate (lqo       (jpmolev,count_profile),stat= alloc_status(6))
      allocate (zho       (jpmolev,count_profile),stat= alloc_status(7))
      allocate (toext     (nlevels  ,count_profile),stat= alloc_status(8))
      allocate (qoext     (nlevels  ,count_profile),stat= alloc_status(9))
      allocate (zvlev     (nlv_T,count_profile),stat= alloc_status(10))
      allocate (zt        (nlv_T,count_profile),stat= alloc_status(11))
      allocate (zlq       (nlv_T,count_profile),stat= alloc_status(12))
      allocate (zht       (nlv_T,count_profile),stat= alloc_status(13))
      call utl_checkAllocationStatus(alloc_status, " tvs_fillProfiles")
      
      count_profile = 0

! second loop over all obs.
      bobs2: do iobs = 1, NOBMAX
        if (tvs_lsensor(iobs) /= sensor_id) cycle bobs2
        count_profile = count_profile + 1
        iptobs(count_profile) = iobs
        header_index = tvs_lobsno(iobs)
        iptobscma(count_profile) = header_index

        !    extract land/sea/sea-ice flag (0=land, 1=sea, 2=sea-ice)
        ksurf = obs_headElem_i(lobsSpaceData,OBS_OFL,header_index)
        tvs_profiles(iobs) % skin % surftype = ksurf

        !    extract satellite zenith and azimuth angle, 
        !    sun zenith angle, cloud fraction, latitude and longitude
        isatzen = obs_headElem_i(lobsSpaceData,OBS_SZA,header_index)
        tvs_profiles(iobs) % zenangle   = (isatzen - 9000) / 100.0

        !pour ne pas faire planter RTTOV dans le cas (rare) ou isatzen n'est pas defini ou invalide         
        if (tvs_profiles(iobs) % zenangle < 0.0d0 .or. &
             tvs_profiles(iobs) % zenangle > zenmax ) then
          write(*,*) "!!! WARNING !!!"
          write(*,*) "INVALID ZENITH ANGLE"
          write(*,*) "angle, profile number, sensor", tvs_profiles(iobs) % zenangle, iobs, sensor_id
          write(*,*) "replaced by 0.0 !!!"
          tvs_profiles(iobs) % zenangle = 0.d0
        end if
 
        isatazim = obs_headElem_i(lobsSpaceData,OBS_AZA,header_index) ! Satellite Azimuth Angle
        isunazim = obs_headElem_i(lobsSpaceData,OBS_SAZ,header_index) ! Sun Azimuth Angle
        tvs_profiles(iobs) % azangle   = ( isatazim / 100.0d0 )
        tvs_profiles(iobs) % sunazangle  =  ( isunazim / 100.0d0 )! necessaire pour radiation solaire
        iplatform = tvs_coefs(sensor_id) % coef % id_platform
        instrum = tvs_coefs(sensor_id) % coef % id_inst
        if ( (instrum == inst_id_amsua .or. instrum == inst_id_mhs) .and. iplatform /= platform_id_eos ) then
          !Correction sur la definition de l'angle. A ammeliorer. Ok pour l'instant.
          tvs_profiles(iobs) % azangle   = (isatazim + isunazim) / 100.d0
          if ( tvs_profiles(iobs) % azangle > 360.d0 ) tvs_profiles(iobs) % azangle = tvs_profiles(iobs) % azangle - 360.d0
        end if
        isunza = obs_headElem_i(lobsSpaceData,OBS_SUN,header_index)
        tvs_profiles(iobs) % sunzenangle = (isunza - 9000) / 100.0d0
        zlat(count_profile) = obs_headElem_r(lobsSpaceData,OBS_LAT,header_index) *MPC_DEGREES_PER_RADIAN_R8
        zlon = obs_headElem_r(lobsSpaceData,OBS_LON,header_index) *MPC_DEGREES_PER_RADIAN_R8
        tvs_profiles(iobs) % longitude = zlon
        do jl = 1, nlv_T
          zt   (jl,count_profile) = col_getElem(columnghr,jl,header_index,'TT')
          zlq  (jl,count_profile) = col_getElem(columnghr,jl,header_index,'HU')
          zvlev(jl,count_profile) = col_getPressure(columnghr,jl,header_index,'TH') * MPC_MBAR_PER_PA_R8
          zht  (jl,count_profile) = col_getHeight(columnghr,jl,header_index,'TH') / rg
        end do

        if (diagTtop) then
          ! Fix temporaire (?) pour eviter probleme au toit avec GEM 4: on ne veut pas utiliser
          ! le premier niveau de GEM qui est disgnostique (extrapole a partir des deux niveaux plus bas)
          ! (grosse varibilite de la temperature au dernier niveau thermo due a l'extrapolation utilisee)
          zt   (1,count_profile) =  zt   (2,count_profile) + tvs_mesosphereLapseRate *  &
               log( col_getPressure(columnghr,1,header_index,'TH') /  &
               col_getPressure(columnghr,2,header_index,'TH') )
          zlq  (1,count_profile) =  zlq  (2,count_profile)         ! extrapolation valeur constante pour H2O peu important a cette hauteur
!!!!
        end if
        
      end do bobs2
 
      !     .  2.1  Vertical interpolation of model temperature, logarithm of
      !             specific humidity and height levels to pressure levels
      !             required by tovs rt model
      !     .       ------------------------------------------


!$omp parallel do private(jn)
      do jn=1, count_profile
        call ppo_IntAvg (zvlev(:,jn:jn),zt(:,jn:jn),nlv_T,nlv_T,1, &
             jpmolev,xpres(jpmotop:nlevels),to(:,jn:jn))
        call ppo_IntAvg (zvlev(:,jn:jn),zlq(:,jn:jn),nlv_T,nlv_T,1, &
             jpmolev,xpres(jpmotop:nlevels),lqo(:,jn:jn))
      end do
!$omp end parallel do


      !     .  2.2  Extrapolation of temperature profile above model top
      !     .       ----------------------------------------------------
      toext(:,:) = 0.0d0
      if ( .not. TopAt10hPa) then ! si le toit n'est pas a 10. hPa 
        do jn=1,count_profile
          toext(jpmotop:nlevels,jn) = to(1:jpmolev,jn)
          ! New approach based on a specified lapse rate
          do jl=1,jpmotop-1
            toext(jl,jn) = toext(jpmotop,jn)  + &
                 ( log(xpres(jl)/xpres(jpmotop)) * tvs_mesosphereLapseRate )
          end do
        end do
      else
        ! old code for temperature profile extrapolation (only apropriate if model top at 10 hPa)
        call extrap (to,toext,jpmolev,nlevels,count_profile)
      end if

      !     .  2.4  Extrapolation of humidity profile (kg/kg)
      !             above rlimlvhu (normally 300mbs or 70mbs)
      !     .       -----------------------------------------

      qoext(:,:) = 0.0d0
     
      do jn = 1, count_profile
        do jl = 1, jpmolev
          qoext(nlevels - jpmolev + jl,jn) = exp(lqo(jl,jn)) 
        end do
      end do

      if ( .not. TopAt10hPa ) then ! if model top not at 10. hPa
        qoext(1:jpmotop,1:count_profile) = MPC_MINIMUM_HU_R8 ! to replace with LIMLVHU ?
      else                    
        if ( tvs_debug ) then
          do jn = 1, count_profile
            write(*,*)'qoext*1000 avant exthum4    = '
            write(*,'(1x,10f8.4)')(qoext(i,jn)*1000.d0,i=1,nlevels)
            write(*,*)' '
          end do
        end if
        call exthum4 (count_profile,nlevels,xpres(1:nlevels),qoext,LIMLVHU)
        if ( tvs_debug ) then
          do jn = 1, count_profile
            write(*,*)'qoext*1000 apres exthum4    = '
            write(*,'(1x,10f8.4)')(qoext(i,jn)*1000.d0,i=1,nlevels)
            write(*,*)' '
          end do
        end if
      end if


      !     .  2.5  Get ozone profiles (ppmv)
!     .       -------------------------
      if (tvs_coefs(sensor_id) %coef % nozone > 0) then
        allocate ( toto3obs(count_profile) )     
        toto3obs(:) = 0.d0
        allocate( PP(nlevels,count_profile) )
        do J=1,count_profile
          PP(1:nlevels,J)=xpres(1:nlevels)
        end do
        call ozo_get_profile (ozo,toto3obs,zlat,pp,nlevels,count_profile,datestamp)
        deallocate( PP )
        deallocate ( toto3obs )
      end if


        !     .  2.6  Fill profiles structure
        !     .       -----------------------

      do  j = 1 , count_profile 
        jj = iptobs(j)
        tvs_profiles(jj) % gas_units       = gas_unit_specconc ! all gas profiles are supposed to be provided in kg/kg (specific humidity, i.e. mass mixing ratio [kg/kg] over wet air)
        tvs_profiles(jj) % id              = "" ! profile id, up to 128 characters, to consider for use
        tvs_profiles(jj) % nlevels         = nlevels
        tvs_profiles(jj) % nlayers         = nlevels - 1
        
        tvs_profiles(jj) % date(1)         = year
        tvs_profiles(jj) % date(2)         = month
        tvs_profiles(jj) % date(3)         = day
        
        
        tvs_profiles(jj) % latitude        = zlat(j)
        
        tvs_profiles(jj) % elevation       = 0.001d0 * zht(ilowlvl_T,j) ! unite km

        tvs_profiles(jj) % skin % watertype= 1 !utilise pour calcul rayonnement solaire reflechi seulement
        tvs_profiles(jj) % skin % t        = col_getElem(columnghr,1,iptobscma(j),'TG')
        tvs_profiles(jj) % skin % salinity = 35.d0 ! for FASTEM-4 only to revise (practical salinity units)
        tvs_profiles(jj) % skin % fastem(:)= 0.0d0
        tvs_profiles(jj) % skin % snow_fraction  = 0.d0 ! Surface coverage snow fraction(0-1), used only by IR emissivity atlas
        tvs_profiles(jj) % skin % soil_moisture  = 0.d0 ! soil moisure (m**3/m**3) not yet used
!
        tvs_profiles(jj) % s2m % t         = col_getElem(columnghr,ilowlvl_T,iptobscma(j),'TT')
!        tvs_profiles(jj) % s2m % q         = exp(col_getElem(columnghr,ilowlvl_T,iptobscma(j),'HU')) * qMixratio2ppmv
        tvs_profiles(jj) % s2m % q         = 0.3D6  * qppmv2Mixratio ! a value between 0 and 0.6d6 so that RTTOV will not complain; not used
        tvs_profiles(jj) % s2m % p         = col_getElem(columnghr,1      ,iptobscma(j),'P0')*MPC_MBAR_PER_PA_R8
        tvs_profiles(jj) % s2m % u         = col_getElem(columnghr,ilowlvl_M,iptobscma(j),'UU')
        tvs_profiles(jj) % s2m % v         = col_getElem(columnghr,ilowlvl_M,iptobscma(j),'VV')
        tvs_profiles(jj) % s2m % o         = 0.0d0 !surface ozone never used
        tvs_profiles(jj) % s2m % wfetc     = 100000.0d0 ! Wind fetch (in meter for rttov10 ?) used to calculate reflection of solar radiation by sea surface
!
        tvs_profiles(jj) % idg             = 0
        
        tvs_profiles(jj) % Be              = 0.4d0 ! earth magnetic field strength (gauss) (must be non zero)
        tvs_profiles(jj) % cosbk           = 0.0d0 ! cosine of the angle between the earth magnetic field and wave propagation direction
          
        tvs_profiles(jj) % p(:)            = tvs_coefs(sensor_id) %coef% ref_prfl_p(:)
        tvs_profiles(jj) % t(:)            = toext(:,j)
        if (tvs_coefs(sensor_id) %coef %nozone > 0) tvs_profiles(jj) % o3(:) = ozo(:,j) * o3ppmv2Mixratio ! Climatology output is ppmv (over dry or wet air? not sure but this conversion is only approximate but it should not matter                                                                                                             ! because atmosphere is very dry where there is significant absorption by ozone)
        tvs_profiles(jj) % q(:)            = qoext(:,j)
        tvs_profiles(jj) % ctp = 1013.25d0
        tvs_profiles(jj) % cfraction = 0.d0
        
      end do

      deallocate (xpres     ,stat= alloc_status(1))
      deallocate (zht       ,stat= alloc_status(2))
      deallocate (zlq       ,stat= alloc_status(3))
      deallocate (zt        ,stat= alloc_status(4))
      deallocate (zvlev     ,stat= alloc_status(5))
      deallocate (qoext     ,stat= alloc_status(6))
      deallocate (toext     ,stat= alloc_status(7))
      deallocate (zho       ,stat= alloc_status(8))
      deallocate (lqo       ,stat= alloc_status(9))
      deallocate (to        ,stat= alloc_status(10))
      deallocate (ozo       ,stat= alloc_status(11))
      deallocate (zlat      ,stat= alloc_status(12))
      deallocate (iptobscma ,stat= alloc_status(13))
      deallocate (iptobs    ,stat= alloc_status(14))
    
      call utl_checkAllocationStatus(alloc_status, " tvs_fillProfiles", .false.)
     
    end do sensor_loop

  end subroutine tvs_fillProfiles

!--------------------------------------------------------------------------
!!**** *exthum4* - extrapolate upper level humidity profile.
!!                 (adapted from exthum by J. Eyre)
!!
!!     purpose.
!      --------
!!          to extend mixing ratio profile into stratosphere in
!!          a reasonable way.
!!
!!**   interface.
!      ----------
!!          *call* *exthum4(knpf,klapf,ppres,pav)*
!!               *knpf*:  no. of profiles to be processed.
!!               *klapf*: length of atm. profiles.
!!               *ppres*: pressure levels of atm. profiles.
!!               *pav*:   humidity profiles.
!!
!!     method.
!      -------
!!          take top tropospheric mixing ratio (e.g. near 300 mb) and
!!          extrapolate with given fall off into lower stratosphere
!!          (e.g. to 70 mb).  constrain mixing ratio to be >= zwmin
!!          (e.g. 0.000003 kg/kg).   in upper strat, mixing ratio = zwmin.
!!
!!     externals.
!      ----------
!!          none.
!!
!!     reference.
!      ----------
!!          ecmwf tech mem 176.
!--------------------------------------------------------------------------

  subroutine EXTHUM4(KNPF,KLAPF,PPRES,PAV,LIMLVHU)
    implicit none
   
    integer klapf, knpf
    real(8) PPRES(*),PAV(KLAPF,*)
    real(8) :: LIMLVHU
    real(8) :: ZPRES3(KLAPF)
    real(8) zwb
    real(8),parameter :: ZP1 = 70.0D0  !PRESS LIMITS (IN HPA) OF REGION to be extrapolated
    integer :: j, inlvw, jnpf
    
    !
    !          find top level of given profile
    inlvw = 0
    do J=KLAPF,1,-1
      if (PPRES(J) < LIMLVHU) THEN
        inlvw = J
        exit
      end if
    end do
    !
    !** Null extrapolation case
    !
    if (inlvw == 0) then
      return
    else
      !
      !          constants defining p**3 fall off around tropopause
      do J=1,inlvw
        ZPRES3(J)=(PPRES(J)/PPRES(inlvw+1))**3
      end do

      do JNPF=1,KNPF
        ZWB=PAV(inlvw+1,JNPF)
        do J=1,inlvw
          if (PPRES(J)<ZP1) THEN
            PAV(J,JNPF)=MPC_MINIMUM_HU_R8
          else
            PAV(J,JNPF)=max((ZWB*ZPRES3(J)),MPC_MINIMUM_HU_R8)
          end if
        end do
      end do
    end if

  end subroutine EXTHUM4

!--------------------------------------------------------------------------
!!**ID HTEXTRAP -- EXTRAPOLATION OF HEIGHT PROFILES
!!
!!       AUTHOR:   A. BEAULNE (CMDA/SMC) March 2006
!!
!!       REVISION:
!!
!!       OBJECT:  EXTRAPOLATE HEIGHT PROFILES ABOVE 10MB MODEL TOP
!!                ON RTTOV LEVELS UP TO 0.1MB (RTTOV LEVELS 1 TO 7)
!!                USING 10 RTTOV HEIGHT LEVELS FROM 100MB TO 10MB
!!                (RTTOV LEVELS 8 TO 17) FOR LINEAR FIT.
!!
!!                -- LINEAR EXTRAPOLATION FOLLOWING
!!                --    PROFOUT(m) = A * ln(XPRES(mb)) + B
!!                -- AND SOLVE A AND B BY LEAST SQUARE METHOD 
!!
!!       ARGUMENTS:
!!          INPUT:
!!            -PROFIN(JPMOLEV,NPRF) :: HEIGHT PROFILES  -TO BE EXTRAPOLATED- (M)
!!            -XPRES(JPLEV)         :: PRESSURE LEVELS OF RTTOV MODEL (HPA)
!!            -JPLEV                :: NUMBER OF PRESSURE LEVELS OF RTTOV MODEL
!!            -JPMOLEV              :: NUMBER OF RTTOV MODEL LEVELS BELOW NWP MODEL TOP
!!            -JPMOTOP              :: FIRST RTTOV MODEL LEVEL UNDER NWP MODEL TOP
!!            -NPRF                 :: NUMBER OF PROFILES
!!
!!          OUTPUT:
!!            -PROFOUT(JPLEV,NPRF)  :: HEIGHT PROFILES  -EXTRAPOLATED- (M)
!--------------------------------------------------------------------------

  subroutine HTEXTRAP ( PROFOUT, profin,xpres,jplev,jpmolev,jpmotop,nprf )
    implicit none


    integer      :: I, JK, JN, NPRF, JPMOLEV, JPMOTOP, JPLEV
    real(8)      :: LNX_SUM, LNX_AVG, Y_SUM, Y_AVG, A_NUM, A_DEN, A, B
    real(8)      :: XPRES(JPLEV), PROFIN(JPMOLEV,NPRF), PROFOUT(JPLEV,NPRF)

    integer, parameter :: NL = 10  ! number of points used in the extrapolation


    do JN = 1, NPRF
      
      LNX_SUM = 0.d0
      Y_SUM   = 0.d0
      A_NUM   = 0.d0
      A_DEN   = 0.d0


      !*      FIND AVERAGED VALUES OF HEIGHT AND LN ( PRESSURE )

      do I = 1, NL
        LNX_SUM = LNX_SUM + LOG(XPRES(JPMOTOP + I - 1))
        Y_SUM   = Y_SUM   + PROFIN(I,JN)
      end do

      LNX_AVG = LNX_SUM / NL
      Y_AVG   = Y_SUM   / NL

      !*      FIND CONSTANTS A AND B BY LEAST-SQUARE METHOD
      
      do I = 1, NL
        A_NUM = A_NUM + ( LOG(XPRES(JPMOTOP + I - 1)) - LNX_AVG ) * ( PROFIN(I,JN) - Y_AVG )
        A_DEN = A_DEN + ( LOG(XPRES(JPMOTOP + I - 1)) - LNX_AVG )**2
      end do

      A = A_NUM / A_DEN

      B = Y_AVG - A * LNX_AVG


      !*      INITIALIZE HEIGHT FOR RTTOV LEVELS UNDER NWP MODEL TOP

      do JK = 1, JPMOLEV
        PROFOUT(JPLEV - JPMOLEV + JK,JN) = PROFIN(JK,JN)
      end do


      !*      EXTRAPOLATE HEIGHT FOR RTTOV LEVELS ABOVE NWP MODEL TOP

      do JK = 1, JPMOTOP - 1
        PROFOUT(JK,JN) = A * LOG(XPRES(JK)) + B
      end do

    end do

  end subroutine HTEXTRAP

!--------------------------------------------------------------------------
!! Interface for RTTOV non linear operator
!! tvs_fillProfiles should be called before
!--------------------------------------------------------------------------

  subroutine tvs_rttov(lcolumnhr,lobsSpaceData,bgckMode,beSilent)
    implicit none

    type(struct_obs) :: lobsSpaceData
    type(struct_columnData) :: lcolumnhr
    logical :: bgckMode, beSilent
  
    integer :: ichn
    integer :: isurface
    integer :: nlevels
    integer :: count_tb
    integer :: alloc_status(6)
    integer :: rttov_err_stat ! rttov error return code
    integer ,external :: omp_get_num_threads
    integer :: nthreads,max_nthreads
    integer :: sensor_id, obs_index
    integer :: channel_index
    integer :: count_profile
    integer :: profile_index, level_index, jj, tb_index
    integer :: instrum
    integer :: sensor_type        ! sensor type (1=infrared; 2=microwave; 3=high resolution,4=polarimetric)
    integer ,allocatable:: iptobs  (:)  
    
    type (rttov_emissivity), allocatable :: emissivity_local (:)    ! emissivity structure with input and output
    type (rttov_chanprof), allocatable :: chanprof(:), chanprof1(:)
    type (rttov_radiance) :: radiancedata_d, radiancedata_d1
    type (rttov_transmission) :: transmission, transmission1
    type (rttov_emis_atlas_data),allocatable,save :: Atlas(:)
    logical ,save        :: first=.true.
    logical              :: init
    integer              :: asw,imonth
    logical, allocatable :: calcemis  (:)
    real*8, allocatable  :: surfem1 (:)
    real*8, allocatable  :: surfem2 (:)
    integer              :: profileIndex2, tb1, tb2
!*****************************************************************

    if ( .not. beSilent ) write(*,*) "Entering tvs_rttov subroutine"
    if ( .not. beSilent ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if (tvs_nobtov == 0) return                  ! exit if there are not tovs data

    !     1.  Get number of threads available and allocate memory for some variables
    !     .   ---------------------------------------------------------------------- 
    !           
!$omp parallel
    max_nthreads = omp_get_num_threads()
!$omp end parallel
    alloc_status(1) = 0
    allocate(iptobs(tvs_nobtov),stat=alloc_status(1))
    call utl_checkAllocationStatus(alloc_status(1:1), " tvs_rttov iptobs")
    
    
    !     1.1   Read surface information
    !     .     ------------------------
    if ( bgckMode ) call EMIS_READ_CLIMATOLOGY

    !
    !     2.  Computation of hx for tovs data only
    !     .   ------------------------------------


    ! Loop over all sensors specified by user
    sensor_loop:do sensor_id = 1, tvs_nsensors
   
      nlevels = tvs_coefs(sensor_id)% coef % nlevels
      sensor_type = tvs_coefs(sensor_id) % coef % id_sensor
      instrum = tvs_coefs(sensor_id) % coef % id_inst
    
      !  loop over all obs.
      count_profile = 0
      obs_loop: do obs_index = 1, tvs_nobtov
        
        !    Currently processed sensor?
        if ( tvs_lsensor(obs_index) == sensor_id ) then
          count_profile = count_profile + 1
          iptobs(count_profile) = obs_index
        end if
      end do obs_loop
      
      if (count_profile == 0) cycle sensor_loop
      
      !     .  2.1  Calculate the actual number of threads which will be used.
      !     .       ----------------------------------------------------------
      
      nthreads = min(max_nthreads, count_profile )  
      
      !     .  2.2  Prepare all input variables required by rttov.
      !     .       ---------------------------------------------------------
      
      if ( bgckMode .and. tvs_isInstrumHyperSpectral(instrum) ) then
        count_tb = count_profile * tvs_nchan(sensor_id)
      else
        count_tb = tvs_countRadiances(iptobs, count_profile, lobsSpaceData)
      end if
      
      if ( count_tb == 0 ) cycle sensor_loop
      
      alloc_status(:) = 0
      allocate ( surfem1          (count_tb) ,stat=alloc_status(2))
      allocate ( chanprof         (count_tb) ,stat=alloc_status(3))
      allocate ( emissivity_local (count_tb) ,stat=alloc_status(4))
      allocate ( calcemis         (count_tb) ,stat=alloc_status(5))
      if (useUofWIREmiss) then
        allocate ( surfem2(count_tb)        ,stat=alloc_status(6) )
      end if
      call utl_checkAllocationStatus(alloc_status, " tvs_rttov")
      
      
      !     get Hyperspectral IR emissivities
      if ( tvs_isInstrumHyperSpectral(instrum) ) then
        surfem1(:) = 0.
        if ( bgckMode ) then
          call EMIS_GET_IR_EMISSIVITY (SURFEM1,tvs_nchan(sensor_id),sensor_id,count_profile,count_tb,iptobs)
        else
          call tvs_getHIREmissivities(sensor_id, iptobs, count_profile, lobsSpaceData, surfem1)
        end if
      end if
      
      if ( bgckMode .and. tvs_isInstrumHyperSpectral(instrum) ) then
        ichn = 0
        do profile_index = 1 , count_profile
          do  channel_index = 1,tvs_nchan(sensor_id)
            ichn = ichn + 1
            chanprof(ichn)%prof = profile_index
            chanprof(ichn)%chan = channel_index
          end do
        end do
      else
        call TVS_getChanprof(sensor_id, iptobs, count_profile, lobsSpaceData, chanprof)
      end if
      
      asw = 1 ! 1 to allocate,0 to deallocate
      ! allocate transmitance structure
      call rttov_alloc_transmission(alloc_status(1), transmission, nlevels=nlevels, &
           nchanprof=count_tb, asw=asw, init=.true. )
      call utl_checkAllocationStatus(alloc_status, " tvs_rttov transmittances")
      
      ! allocate radiance structure
      call rttov_alloc_rad (alloc_status(1),count_tb, radiancedata_d,nlevels,asw,init=.true.)
      call utl_checkAllocationStatus(alloc_status, " tvs_rttov radiances")
      
      
      call tvs_getOtherEmissivities(chanprof, iptobs, count_tb, sensor_type, instrum, surfem1, calcemis)
      
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
          if (rttov_err_stat/=0) THEN
            write(*,*) "Error in rttov_atlas_setup ",rttov_err_stat
            call utl_abort('tvs_rttov')
          end if
          first=.false.
        end if
        
        call rttov_get_emis( rttov_err_stat        , &   ! out
             tvs_opts(sensor_id)                   , &   ! in
             chanprof(1:count_tb)                  , &   ! in
             tvs_profiles(iptobs(1:count_profile)) , &   ! in
             tvs_coefs(sensor_id)                  , &   ! in
             Atlas(sensor_id)                      , &   ! in
             surfem2(1:count_tb)     ) ! out

        if (rttov_err_stat /= 0) THEN
          write(*,*) "Error in rttov_get_emis ", rttov_err_stat
          call utl_abort('tvs_rttov')
        end if
              
        do profile_index=1, count_profile !loop on profiles
          jj = iptobs(profile_index)
          do tb_index=1, count_tb !loop on channels
            if (chanprof(tb_index)%prof==profile_index) then
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
                   surfem2(tb_index) > 0.5 ) then
                emissivity_local(tb_index)%emis_in = surfem2(tb_index)
              else
                emissivity_local(tb_index)%emis_in = surfem1(tb_index)
              end if
            end if
          end do
        end do
      else
        emissivity_local(:)%emis_in = surfem1(:)
      end if
        
      !     .  2.3  Compute radiance with rttov direct
      !     .       ----------------------------------

      rttov_err_stat = 0 

      if( bgckMode .and. tvs_isInstrumHyperSpectral(instrum) ) then
        write(*,*) 'for bgck IR: call rttov_parallel_direct for each profile...'

        asw = 1 ! 1 to allocate,0 to deallocate
        ! allocate transmitance structure for 1 profile
        call rttov_alloc_transmission(alloc_status(1), transmission1, nlevels=nlevels, &
             nchanprof=tvs_nchan(sensor_id), asw=asw, init=.true. )
        ! allocate radiance structure for 1 profile
        call rttov_alloc_rad (alloc_status(1),tvs_nchan(sensor_id), radiancedata_d1,nlevels,asw,init=.true.)
        ! allocate chanprof for 1 profile
        allocate(chanprof1(tvs_nchan(sensor_id)))
        do  channel_index = 1,tvs_nchan(sensor_id)
          chanprof1(channel_index)%prof = 1
          chanprof1(channel_index)%chan = channel_index
        end do

        do profileIndex2 = 1, count_profile
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
        call rttov_alloc_transmission(alloc_status(1),transmission1,nlevels=nlevels,  &
             nchanprof=tvs_nchan(sensor_id), asw=asw )
        ! radiance deallocation for 1 profile
        call rttov_alloc_rad (alloc_status(1), tvs_nchan(sensor_id), radiancedata_d1, nlevels, asw)

      else

        call rttov_parallel_direct(       &
             rttov_err_stat,              & ! out
             chanprof,                    & ! in
             tvs_opts(sensor_id),         & ! in
             tvs_profiles(iptobs(1:count_profile)),  & ! in
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
                                        
      !     .  2.4  Store hx in the structure tvs_radiance
      !     .       ------------------------------------

      do tb_index = 1, count_tb
        profile_index = chanprof(tb_index)%prof
        ichn = chanprof(tb_index)%chan
        obs_index = iptobs(profile_index)
        tvs_radiance(obs_index) % bt(ichn) =     &
             radiancedata_d % bt(tb_index)
        tvs_radiance(obs_index) % clear(ichn) =  &
             radiancedata_d %clear(tb_index)
        if ( bgckMode ) then
          do level_index = 1, nlevels - 1
            tvs_radiance(obs_index) % overcast(level_index,ichn) =   &
                 radiancedata_d % overcast(level_index,tb_index)
          end do
          do level_index = 1, nlevels
            tvs_transmission(obs_index) % tau_levels(level_index,ichn) = &
                 transmission % tau_levels(level_index,tb_index)
          end do
          
          tvs_transmission(obs_index) % tau_total(ichn) = &
               transmission % tau_total(tb_index)
          tvs_emissivity(ichn,obs_index) = emissivity_local(tb_index)%emis_out
          
        end if
        
      end do

      !     de-allocate memory
      asw = 0 ! 1 to allocate,0 to deallocate
      ! transmittance deallocation
      call rttov_alloc_transmission(alloc_status(1),transmission,nlevels=nlevels,  &
           nchanprof=count_tb, asw=asw )
      call utl_checkAllocationStatus(alloc_status, " tvs_rttov transmittances", .false.)
      ! radiance deallocation      
      call rttov_alloc_rad (alloc_status(1), count_tb, radiancedata_d, nlevels, asw)
      call utl_checkAllocationStatus(alloc_status, " tvs_rttov radiances", .false.)
      
      if (useUofWIREmiss) then
        deallocate ( surfem2         ,stat=alloc_status(1) )
      end if
      deallocate ( emissivity_local  ,stat=alloc_status(2) )
      deallocate ( calcemis          ,stat=alloc_status(3) )
      deallocate ( chanprof          ,stat=alloc_status(4) )
      deallocate ( surfem1           ,stat=alloc_status(5) )
      call utl_checkAllocationStatus(alloc_status, " tvs_rttov", .false.)
      
    enddo sensor_loop
    
    deallocate(iptobs)

  end subroutine tvs_rttov

!--------------------------------------------------------------------------
!!**ID COMP_IR_EMISS -- INFRARED EMISSIVITY COMPUTATION
!!
!!       AUTHOR:   Thomas J. Kleespies               8 February 1998
!!                 Physics Branch
!!                 Satellite Research Laboratory
!!                 Office of Research and Applications
!!                 NOAA/NESDIS
!!                 301-763-8136 x126
!!                 301-763-8108 FAX
!!                 Mailing Address: 810 NSC E/RA-14
!!                                  NOAA/NESDIS
!!                                  Washington, D.C. 20233
!!                 Email: TKleespies@nesdis.noaa.gov
!!
!!                 L. GARAND     modified for NP points
!!                 A. BEAULNE (CMDA/SMC)       April 2006  (ADAPT TO 3DVAR)
!!
!!       REVISION:
!!
!!       OBJECT:   COMPUTES WATER INFRARED EMISSIVITY FOR A SPECIFIC SET OF
!!          CHANNEL INDICES, WIND SPEED AND ZENITH ANGLE.
!!
!!          Restrictions:  Must be compiled with /EXTend_SOURCE or it's equivalent
!!
!!       ARGUMENTS:
!!          INPUT:
!!            -WIND(NP)         : SURFACE WIND SPEED (M/S)
!!            -ANGLE(NP)        : VIEWING ANGLE (DEG)
!!            -NCHN             : NUMBER OF CHANNELS TO PROCESS
!!            -NP               : NUMBER OF LOCATIONS
!!            -MCHANNEL(NCHN)   : VECTOR OF CHANNEL INDICES TO PROCESS
!!
!!          OUTPUT:
!!            -EMISS(NCHN,NP)   : EMISSIVITIES (0.-1.)
!--------------------------------------------------------------------------

  subroutine COMP_IR_EMISS (EMISS, wind,angle,nchn,np,mchannel)
    implicit None
    integer ,intent(in) :: nchn,np
    real (8)    ,intent(out):: Emiss(Nchn,NP)
    real (8)    ,intent(in) :: Wind(NP),Angle(NP)
    integer ,intent(in) :: Mchannel(Nchn)
!***********************************************
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

    real (8) A(MaxChan),B(MaxChan),CC(MaxChan)  ! local variable
    real (8) WW
    integer Index,Ichan,IP


    do Ichan = 1 , Nchn

      Index = Mchannel(Ichan)

      do IP=1,NP

        WW = WIND(IP)
        A(Ichan) = C(1,1,Index) + C(2,1,Index) * WW    &  
             + C(3,1,Index) * WW * WW
        B(Ichan) = C(1,2,Index) + C(2,2,Index) * WW    &
             + C(3,2,Index)* WW * WW

        CC(Ichan) = Theta(1,Index) + Theta(2,Index) * WW

        Emiss(Ichan,IP) = A(Ichan) + (B(Ichan) - A(Ichan)) *   & 
             Exp(( (Theta(3,Index) - 60.d0)**2.d0              &
             - (Angle(IP) - Theta(3,Index))**2.d0 ) / CC(Ichan))
       
      end do
      
    end do

  end subroutine COMP_IR_EMISS

!--------------------------------------------------------------------------
!!**ID PCNT_BOX -- COMPUTES A LOW_RESOLUTION FEATURE FROM HIGH RESOLUTION
!!
!!       AUTHOR:   L. GARAND (ARMA) AND A. BEAULNE (CMDA/SMC) June 2006
!!
!!       REVISION:
!!
!!       OBJECT:   COMPUTES A LOW RESOLUTION FEATURE FORM A HIGH
!!                 RESOLUTION ONE BY AVERAGING.
!!                 EXAMPLE: USE FOR PERCENTAGE OF WATER
!!
!!
!!       ARGUMENTS:
!!          INPUT:
!!            -F_HIGH(KLON,KLAT)   : HIGH RESOLUTION FIELD 
!!            -NPRF                : NUMBER OF PROFILES
!!            -ILAT(NPRF)          : Y-COORDINATE OF PROFILE
!!            -ILON(NPRF)          : X-COORDINATE OF PROFILE
!!            -KLAT                : MAX VALUE OF LATITUDE INDICES
!!            -KLON                : MAX VALUE OF LONGITUDE INDICES
!!            -IREDUC              : MEANS A 2xIREDUC+1 BY 2xIREDUC+1 AVERAGING
!!
!!          OUTPUT:
!!            -FLOW(NPRF)          : LOW RESOLUTION FIELD
!--------------------------------------------------------------------------

  subroutine PCNT_BOX(F_LOW, f_high,nprf,ilat,ilon,klat,klon,ireduc)
    implicit none
    integer ,intent(in) :: NPRF,KLON,KLAT,ireduc
    integer ,intent(in) :: ILAT(NPRF), ILON(NPRF)
    real (8),intent(in)    :: F_HIGH(KLON,KLAT)
    real (8),intent(out)   :: F_LOW(NPRF)
!*************************************************************
    integer :: NPLON, JDLO1, JDLO2, JLON1, JLON2
    integer :: NX, ILAT1, ILAT2, ILON1, ILON2, JN, ii, jj
   
    profiles : do JN = 1,NPRF

      NPLON = 0

      ! normal limits

      ilat1=max(ilat(JN)-IREDUC,1)
      ilat2=min(ilat(JN)+IREDUC,KLAT)
      ilon1=max(ilon(JN)-IREDUC,1)
      ilon2=min(ilon(JN)+IREDUC,KLON)

      if (ilon1 == 1 .or. ilon2 == klon) then
        ! border cases for longitudes
        JDLO1 = ILON(JN)-IREDUC
        JDLO2 = ILON(JN)+IREDUC

        if ( JDLO1 <= 0 ) THEN
          NPLON = 1
          JLON1 = KLON + JDLO1
          JLON2 = KLON
        else if ( JDLO2 > KLON ) THEN
          NPLON = 1
          JLON1 = 1
          JLON2 = JDLO2 - KLON
        end if
      end if

      NX = 0
      F_LOW(JN) = 0.d0
     
      do JJ = ILAT1, ILAT2

        do II = ILON1, ILON2
          NX = NX + 1
          F_LOW(JN) = F_LOW(JN) + F_HIGH(II,JJ)         
        end do
        
        if (NPLON == 1) THEN
          ! additional cases at border 1-KLON
          do II = JLON1, JLON2
            NX = NX + 1
            F_LOW(JN) = F_LOW(JN) + F_HIGH(II,JJ)         
          end do
        end if

      end do
      
      F_LOW(JN) = F_LOW(JN) / dble(NX)

    end do profiles

  end subroutine PCNT_BOX

!--------------------------------------------------------------------------
!!**ID emis_read_climatology -- READ INFORMATION FOR IR SURFACE EMISSIVITIES COMPUTATION
!!
!!       AUTHOR:   A. BEAULNE (CMDA/SMC) March 2006
!!
!!       OBJECT:   READ INFORMATION ABOUT CERES SURFACE TYPE AND WATER FRACTION.
!!
!!       REVISION: A. Beaulne (CMDA/SMC) July 2013
!!                 Use only for reading Ceres information,
!!                 so albedo, ice and snow now done in interp_sfc.ftn90.
!!       ARGUMENTS:
!!          INPUT:  NONE
!!          OUTPUT: NONE
!--------------------------------------------------------------------------

  subroutine emis_read_climatology
    implicit none
    
    integer            :: NISF,NJSF,NKSF
    integer            :: NIWA,NJWA,NKWA
    character(len=20)  :: CFILE
    integer,external   :: FNOM,FSTOUV,FSTFRM,FCLOS,FSTLIR
    integer            :: isftest
    integer            :: iv1,iv2,iv3,iv4,iv5,iv6

    isftest = 0


    !* get surface type and water fraction
    CFILE = 'ceres_global.std'
    IV1 = FNOM(ISFTEST,CFILE,'RND+R/O',0)
    IV2 = FSTOUV(ISFTEST,'RND')
    IV3 = FSTLIR(JTYPE,ISFTEST,NISF,NJSF,NKSF,-1,'SFC-TYPE',-1,-1,-1,'','TY')
    IV4 = UTL_FSTLIR(WATERF,ISFTEST,NIWA,NJWA,NKWA,-1,'WATER_FR',-1,-1,-1,'','W%')
    IV5 = FSTFRM(ISFTEST)
    IV6 = FCLOS(ISFTEST)

    if (iv1 < 0 .or. iv2 < 0 .or. iv3 < 0 .or. iv4 < 0 .or. iv5 < 0 .or. iv6 < 0) then
      write(*,*) 'LES IV DE CERES ',iv1,iv2,iv3,iv4,iv5,iv6
      write(*,*) 'THESE NUMBER SHOULD NOT BE NEGATIVE WHEN DOING AIRS BACKGROUND CHECK'
      call utl_abort('Problem with file ceres_global.std in emis_read_climatology ')
    end if
   
  end subroutine emis_read_climatology

!--------------------------------------------------------------------------
!!* This is a subroutine that can apply to any instrument.
!!* However, due to the necessity of specifying the instrument
!!* bands wavenumbers, the use of this subroutine for a new instrument
!!* would require the minor following changes.
!!*
!!*     - Continue the "find the bands (central) wavenumber" if loop
!!*         for your specific instrument
!!
!!**ID EMIS_GET_IR_EMISSIVITY -- ASSIGN NEW IR SURFACE EMISSIVITIES 
!!
!!       SCIENCE:  L. GARAND
!!       AUTHOR:   A. BEAULNE (CMDA/SMC) June 2006
!!
!!       OBJECT:   ASSIGN NEW IR SURFACE EMISSIVITIES BASED ON
!!                 CMC ANALYSIS SURFACE ALBEDO, SEA ICE FRACTION AND SNOW MASK
!!                 IN ADDITION TO CERES SURFACE TYPE AND WATER FRACTION
!!
!!       ARGUMENTS:
!!          INPUT:
!!            -NCHN           : NUMBER OF CHANNELS
!!            -KRTID          : SENSOR NUMBER
!!            -NPRF           : NUMBER OF PROFILES
!!            -NCHANNELS_MAX  : TOTAL NUMBER OF OBSERVATIONS TREATED
!!            -IPTOBS(NPRF)   : PROFILE POSITION NUMBER
!!
!!          OUTPUT:
!!            -SURFEM1(NCHANNELS_MAX)  : IR SURFACE EMISSIVITY ESTIMATE (0-1)
!--------------------------------------------------------------------------

  subroutine EMIS_GET_IR_EMISSIVITY ( SURFEM1, nchn,krtid,nprf,nchannels_max,iptobs)
    implicit none
    integer,intent(in) :: NPRF,NCHANNELS_MAX
    integer,intent(in) :: NCHN,IPTOBS(NPRF),KRTID
    real(8),intent(out) :: SURFEM1(NCHANNELS_MAX)
!****************************************************

    integer :: JC,JN,ICHN
    integer :: ILAT(NPRF), ILON(NPRF)
    real(8) :: ZLAT(NPRF), ZLON(NPRF), SATZANG(NPRF)
      
    real (8) :: WIND_SFC(NPRF), F_LOW(NPRF), WAVEN(NCHN), EM_OC(NCHN,NPRF), EMI_MAT(NCHN,20)


    !* information to extract (transvidage)
    !--------------------------------------
    !
    ! ZLAT(NPRF) -- latitude (-90 to 90)
    ! ZLON(NPRF) -- longitude (0 to 360)
    ! SATZANG(NPRF) -- satellite zenith angle (deg)

    do JN = 1, NPRF
      ZLAT(JN)    = tvs_profiles(IPTOBS(JN))% latitude
      ZLON(JN)    = tvs_profiles(IPTOBS(JN))% longitude
      SATZANG(JN) = tvs_profiles(IPTOBS(JN))% ZENANGLE
    end do

    !     assign surface properties from grid to profiles
    call INTERP_SFC(ILAT,ILON, nprf,zlat,zlon,iptobs)


    !* find the sensor bands (central) wavenumbers
    if ( tvs_instruments(KRTID) == 11 ) THEN ! --AIRS--
      
      do JC = 1, NCHN
        ICHN = tvs_ichan(JC,KRTID)
        WAVEN(JC) = hir_get_wavn("AIRS",ICHN)
      end do

    else if ( tvs_instruments(KRTID) == 16 ) THEN ! --IASI--

      do JC = 1, NCHN
        ICHN = tvs_ichan(JC,KRTID)
        WAVEN(JC) =  hir_get_wavn("IASI",ICHN)
      end do

    else if ( tvs_instruments(KRTID) == 27 ) THEN ! --CrIS--
      
      do JC = 1, NCHN
        ICHN = tvs_ichan(JC,KRTID)
        WAVEN(JC) =  hir_get_wavn("CRIS",ICHN)
      end do


    end if


    !* get the CERES emissivity matrix for all sensor wavenumbers and surface types


    call CERES_EMATRIX(EMI_MAT, waven,nchn)


!* refine water emissivities

    

    do JN = 1, NPRF

!       find surface wind

      WIND_SFC(JN) = min(sqrt(tvs_profiles(IPTOBS(JN))%S2M%U**2 + tvs_profiles(IPTOBS(JN))%S2M%V**2 + 1.d-12),15.d0)
 

    end do

    !     find new ocean emissivities     

    do JC = 1, NCHN
      EM_OC(JC,:)= EMI_MAT(JC,17)
    end do
    
    call EMI_SEA (EM_OC, waven,satzang,wind_sfc,nprf,nchn)
    

!* get surface emissivities

    do JN = 1, NPRF

!       set albedo to 0.6 where snow is present

      if ( tvs_profiles(IPTOBS(JN))%SKIN%SURFTYPE == 0 .AND. tvs_surfaceParameters(IPTOBS(JN))%SNOW > 0.999 ) tvs_surfaceParameters(IPTOBS(JN))%ALBEDO = 0.6

!       if albedo too high no water

      if ( tvs_surfaceParameters(IPTOBS(JN))%ALBEDO >= 0.55 ) tvs_surfaceParameters(IPTOBS(JN))%PCNT_WAT = 0.

!       if water and CMC ice present then sea ice

      if ( tvs_profiles(IPTOBS(JN))%SKIN%SURFTYPE == 1 .and. tvs_surfaceParameters(IPTOBS(JN))%ICE > 0.001 ) tvs_surfaceParameters(IPTOBS(JN))%LTYPE = 20

!       if land and CMC snow present then snow

      if ( tvs_profiles(IPTOBS(JN))%SKIN%SURFTYPE == 0 .and. tvs_surfaceParameters(IPTOBS(JN))%SNOW > 0.999 ) tvs_surfaceParameters(IPTOBS(JN))%LTYPE = 15

      do JC=1,NCHN

        SURFEM1((JN-1)*NCHN+JC) =  tvs_surfaceParameters(IPTOBS(JN))%PCNT_WAT * EM_OC(JC,JN)  +   &
             ( 1.d0 - tvs_surfaceParameters(IPTOBS(JN))%PCNT_WAT ) * EMI_MAT(JC,tvs_surfaceParameters(IPTOBS(JN))%LTYPE)

      end do
      
    end do

!* find the regional water fraction (here in a 15x15 pixel box centered on profile)

    call PCNT_BOX (F_LOW, waterf,nprf,ilat,ilon,kslat,kslon,7)

    do JN = 1, NPRF
      tvs_surfaceParameters(IPTOBS(JN))%PCNT_REG = F_LOW(JN)
    end do


  end subroutine EMIS_GET_IR_EMISSIVITY

!--------------------------------------------------------------------------
!!**ID INTERP_SFC -- ASSOCIATE SURFACE FIELDS TO OBSERVATION PROFILES
!!
!!       AUTHOR:   L. GARAND
!!                 A. BEAULNE (CMDA/SMC) March 2006  (ADAPT TO 3DVAR)
!!
!!       OBJECT:   ASSOCIATE SURFACE ALBEDO, ICE FRACTION, SNOW DEPTH 
!!          AND CERES SURFACE TYPE AND WATER FRACTION TO OBSERVATIONS PROFILES.
!!
!!       REVISION: A. BEAULNE (CMDA/SMC) July 2013
!!                 - GRID COMPUTATION PREVIOUSLY DONE IN SFC_EMISS FOR
!!                   ICE,SNOW AND ALBEDO NOW DONE HERE FOR
!!                   GENERALIZATION IN ACCEPTING ANY KIND OF GRID
!!
!!       ARGUMENTS:
!!          INPUT:
!!            -NPRF           : NUMBER OF PROFILES
!!            -ZLAT(NPRF)     : LATITUDE (-90S TO 90N)
!!            -ZLON(NPRF)     : LONGITUDE (0 TO 360)
!!
!!          OUTPUT:
!!            -ILAT(NPRF)     : Y-COORDINATE OF PROFILE
!!            -ILON(NPRF)     : X-COORDINATE OF PROFILE 
!--------------------------------------------------------------------------

  subroutine INTERP_SFC (ILAT,ILON, nprf,zlat,zlon,iptobs)
    implicit none
    integer,intent(in) :: NPRF, IPTOBS(NPRF)
    real(8),intent(in) :: ZLAT(NPRF), ZLON(NPRF)
    integer,intent(out):: ILAT(NPRF), ILON(NPRF)
!**********************************************************
    character(len=20)  :: CFILE3,CFILE5
    integer            :: iun3,iun5
    integer            ::                     IV6,IV7
    integer            :: IX1,IX2,IX3,IX4,IX5,        IX8,IX9,IX10,IX11,IX12
    integer            ::         IY3,IY4,IY5,        IY8,IY9,IY10
    integer            :: IZ1,IZ2,IZ3,IZ4,IZ5,        IZ8,IZ9,IZ10,IZ11,IZ12
    integer            :: NI3,NJ3,NK3
    integer            :: NI4,NJ4,NK4
    integer            :: NI5,NJ5,NK5
    integer            :: DATEO,DEET,NPAS,NBITS,DATYP
    integer            :: IP1,IP2,IP3
    integer            :: IG13,IG23,IG33,IG43
    integer            :: IG14,IG24,IG34,IG44
    integer            :: IG15,IG25,IG35,IG45
    integer            :: SWA,LNG,DLTF,UBC,EX1,EX2,EX3
    integer            :: JN
    character(len=1)   :: TYPVAR
    character(len=1)   :: GRTYP3,GRTYP4,GRTYP5
    character(len=2)   :: NOMVAR, snowvar
    character(len=8)   :: ETIKET
    integer,external   :: FNOM,FSTOUV,FSTINF,FSTPRM,FSTFRM,FCLOS
    integer,external   :: ezsetopt,ezqkdef,ezdefset
    real(8)            :: zig1,zig2,zig3,zig4
    integer            :: ig1obs,ig2obs,ig3obs,ig4obs
    real (8)           :: ALAT, ALON, ZZLAT, ZZLON

! fields on input grid

    real(8),ALLOCATABLE :: GLACE(:,:), NEIGE(:,:), ALB(:,:)

! fields on output grid

    real(8) :: GLACE_INTRPL(NPRF), NEIGE_INTRPL(NPRF), ALB_INTRPL(NPRF)


    ! printout header
    write(*,*) 
    write(*,*) 'SUBROUTINE INTERP_SFC'
    write(*,*) '---------------------'
    write(*,*) ' called multiple time by bunch of ',nprf,' profiles'
    write(*,*) ' <RETURN CODES> SHOULD NOT BE NEGATIVE'
    write(*,*) '---------------------------------------------------'


    !* --- FOR CERES VARIABLES -------------
    !* get number of pixels per degree of lat or lon

    ALAT = DBLE(KSLAT)/180.d0
    ALON = DBLE(KSLON)/360.d0

    do JN=1, NPRF

      !* get lat and lon within limits if necessary
      ZZLAT = min(ZLAT(JN),89.999d0)
      ZZLAT = max(ZZLAT,-89.999d0)
      
      ZZLON = min(ZLON(JN),359.999d0)
      ZZLON = max(ZZLON,0.d0)

      !* find in which surface field pixel is located the observation profile

      !* Note : CERES grid at 1/6 resolution 
      !*         N-S : starts at N pole and excludes S pole
      !*         W-E : starts at longitude 0 and excludes longitude 360

      ILAT(JN) = max( nint((ZZLAT + 90.d0) * ALAT),1) 
      ILON(JN) = nint(ZZLON * ALON) + 1
      if (ILON(JN) > KSLON) ILON(JN) = 1

!* assign surface caracteristics to observation profiles

      tvs_surfaceParameters(IPTOBS(JN)) % LTYPE    = JTYPE(ILON(JN),ILAT(JN))
      tvs_surfaceParameters(IPTOBS(JN)) % PCNT_WAT = WATERF(ILON(JN),ILAT(JN))

    end do



    !* --- FOR ICE, SNOW AND ALBEDO VARIABLES -------------

    iun3 = 0
    iun5 = 0

    ! files name
    CFILE3 = 'sfc4airs'          ! for ice fraction and snow cover
    CFILE5 = 'sfc4airs_newalb'   ! for albedo


    ! FNOM: make the connections with the external files name
    ! success = 0
    write(*,*) 
    IX1 = FNOM(iun3,CFILE3,'RND+R/O',0)
    write(*,*) 'file = sfc4airs         : FNOM   : return = ', IX1

    IZ1 = FNOM(iun5,CFILE5,'RND+R/O',0)
    write(*,*) 'file = sfc4airs_newalb  : FNOM   : return = ', IZ1


    ! FSTOUV: open the standard files
    ! success = number of records found in the file
    write(*,*) 
    IX2 = FSTOUV(iun3,'RND')
    write(*,*) 'file = sfc4airs         : FSTOUV : return = ', IX2
    IZ2 = FSTOUV(iun5,'RND')
    write(*,*) 'file = sfc4airs_newalb  : FSTOUV : return = ', IZ2


    ! FSTINF: locate the records that matches the search keys
    ! success = handle of the record found after the search
    ! desired output = handle
    write(*,*) 
    IX3 = FSTINF(iun3,NI3,NJ3,NK3,-1,'',-1,-1,-1,'','LG')
    write(*,*) 'variable = LG           : FSTINF : return = ', IX3

    snowvar='SD'
    IY3 = FSTINF(iun3,NI4,NJ4,NK4,-1,'',-1,-1,-1,'',snowvar)
    write(*,*) 'variable = ', snowvar, '           : FSTINF : return = ', IY3
    if ( IY3  <  0 ) then
      write(*,*) 'did not find ''SD'' so look for ''NE'''
      snowvar='NE'
      IY3 = FSTINF(iun3,NI4,NJ4,NK4,-1,'',-1,-1,-1,'',snowvar)
      write(*,*) 'variable = ', snowvar, '           : FSTINF : return = ', IY3
    end if

    IZ3 = FSTINF(iun5,NI5,NJ5,NK5,-1,'',-1,-1,-1,'','AL')
    write(*,*) 'variable = AL           : FSTINF : return = ', IZ3


    ! FSTPRM: get the description informations of the record given the key
    ! success = 0
    ! desired output = NIx,NJx,GRTYPx,IGxx,IG1x,IG2x,IG3x,IG4x

    write(*,*) 
    IX4 = FSTPRM(ix3, DATEO,DEET,NPAS,NI3,NJ3,NK3,NBITS,DATYP, &
         IP1,IP2,IP3,TYPVAR,NOMVAR,ETIKET,GRTYP3,  &
         IG13,IG23,IG33,IG43,SWA,LNG,DLTF,UBC,EX1,EX2,EX3)
    write(*,*) 'variable = LG           : FSTPRM : return = ', IX4

    IY4 = FSTPRM(iy3, DATEO,DEET,NPAS,NI4,NJ4,NK4,NBITS,DATYP, &
         IP1,IP2,IP3,TYPVAR,NOMVAR,ETIKET,GRTYP4,  &
         IG14,IG24,IG34,IG44,SWA,LNG,DLTF,UBC,EX1,EX2,EX3)
    write(*,*) 'variable = ', snowvar, '           : FSTPRM : return = ', IY4

    IZ4 = FSTPRM(iz3, DATEO,DEET,NPAS,NI5,NJ5,NK5,NBITS,DATYP, &
         IP1,IP2,IP3,TYPVAR,NOMVAR,ETIKET,GRTYP5,  &
         IG15,IG25,IG35,IG45,SWA,LNG,DLTF,UBC,EX1,EX2,EX3)
    write(*,*) 'variable = AL           : FSTPRM : return = ', IZ4


    ! allocation of the field on the grid
    allocate ( GLACE  (NI3,NJ3) )
    allocate ( NEIGE  (NI4,NJ4) )
    allocate ( ALB    (NI5,NJ5) )


    ! UTL_FSTLIR: read records data (field on the grid) given the key
    ! success = handle of the record
    ! desired output = FIELD
    write(*,*) 

    IX5 = utl_fstlir(GLACE, iun3,NI3,NJ3,NK3,-1,'',-1,-1,-1,'','LG')
    write(*,*) 'variable = LG           : UTL_FSTLIR : return = ', IX5
    IY5 = utl_fstlir(NEIGE, iun3,NI4,NJ4,NK4,-1,'',-1,-1,-1,'',snowvar)
    write(*,*) 'variable = ', snowvar, '           : UTL_FSTLIR : return = ', IY5
    IZ5 = utl_fstlir(ALB,   iun5,NI5,NJ5,NK5,-1,'',-1,-1,-1,'','AL')
    write(*,*) 'variable = AL           : UTL_FSTLIR : return = ', IZ5


    ! EZSETOPT: set nearest neighbor interpolation option within EZSCINT package
    ! success = 0
    write(*,*) 
    IV6 = ezsetopt('INTERP_DEGREE','NEAREST')
    write(*,*) 'apply to all variables  : ezsetopt : return = ', IV6


    ! UTL_CXGAIG: define the grid descriptors (integer form) of the
    !          observation profile output grid
    ! desired output = IG1OBS, IG2OBS, IG3OBS, IG4OBS
    zig1 = 0.0D0
    zig2 = 0.0D0
    zig3 = 1.0D0
    zig4 = 1.0D0

    call utl_cxgaig('L',IG1OBS,IG2OBS,IG3OBS,IG4OBS,zig1,zig2,zig3,zig4)


    ! UTL_EZGDEF: define the grid of the observations profiles (output grid)
    ! of type Y containing the lat-lon of profiles
    ! success = token to identify the grid
    ! desired output = token
    write(*,*) 
    IV7 = utl_ezgdef(nprf,1,'Y','L',ig1obs,ig2obs,ig3obs,ig4obs,zlon,zlat)
    write(*,*) 'apply to all variables  : UTL_EZGDEF : return = ', IV7
    

    ! EZQKDEF: define the grid of the records data (input grid)
    ! success = token to identify the grid
    ! desired output = token
    ! EZDEFSET: interpolate from input grids to output grid
    ! success = key
    ! UTL_EZSINT: interpolation of the field on the input grid to observation profiles
    ! success = 0
    ! desired output = FIELD_INTRPL
    write(*,*) 
    IX8 = ezqkdef(ni3,nj3,grtyp3,ig13,ig23,ig33,ig43,iun3)
    write(*,*) 'variable = LG           : ezqkdef  : return = ', IX8
    
    IX9 = ezdefset(iv7,ix8)
    write(*,*) 'variable = LG           : ezdefset : return = ', IX9

    IX10 = utl_ezsint(GLACE_INTRPL,glace,nprf,1,1,ni3,nj3,1)
    write(*,*) 'variable = LG           : utl_ezsint  : return = ', IX10

    write(*,*) 

    IY8 = ezqkdef(ni4,nj4,grtyp4,ig14,ig24,ig34,ig44,iun3)
    write(*,*) 'variable = ', snowvar, '           : ezqkdef  : return = ', IY8

    IY9 = ezdefset(iv7,iy8)
    write(*,*) 'variable = ', snowvar, '           : ezdefset : return = ', IY9

    IY10 = utl_ezsint(NEIGE_INTRPL,neige,nprf,1,1,ni4,nj4,1)
    write(*,*) 'variable = ', snowvar, '           : utl_ezsint  : return = ', IY10

    write(*,*) 

    IZ8 = ezqkdef(ni5,nj5,grtyp5,ig15,ig25,ig35,ig45,iun5)
    write(*,*) 'variable = AL           : ezqkdef  : return = ', IZ8

    IZ9 = ezdefset(iv7,iz8)
    write(*,*) 'variable = AL           : ezdefset : return = ', IZ9

    IZ10 = utl_ezsint(ALB_INTRPL,alb,nprf,1,1,ni5,nj5,1)
    write(*,*) 'variable = AL           : utl_ezsint  : return = ', IZ10


    ! FSTFRM: close the standard files
    ! success = 0
    write(*,*) 
    IX11 = FSTFRM(iun3)
    write(*,*) 'file = sfc4airs         : FSTFRM : return = ', IX11
    
    IZ11 = FSTFRM(iun5)
    write(*,*) 'file = sfc4airs_newalb  : FSTFRM : return = ', IZ11
 

    ! FCLOS: release the connections with the external files name
    ! success = 0

    write(*,*) 

    IX12 = FCLOS(iun3)
    write(*,*) 'file = sfc4airs         : FCLOS  : return = ', IX12

    IZ12 = FCLOS(iun5)
    write(*,*) 'file = sfc4airs_newalb  : FCLOS  : return = ', IZ12

    ! assign surface caracteristics to observation profiles

    do JN=1, NPRF
      tvs_surfaceParameters(IPTOBS(JN))%ICE      = GLACE_INTRPL(JN)
      tvs_surfaceParameters(IPTOBS(JN))%SNOW     = NEIGE_INTRPL(JN)
      tvs_surfaceParameters(IPTOBS(JN))%ALBEDO   = ALB_INTRPL(JN)
    end do

    deallocate(GLACE,NEIGE,ALB)

  end subroutine INTERP_SFC

!--------------------------------------------------------------------------
!!**ID CERES_EMATRIX -- SET UP EMISSIVITIES
!!
!!       AUTHOR:   L. GARAND              Sept 2004
!!                 A. BEAULNE (CMDA/SMC) March 2006  (ADAPT TO 3DVAR)
!!
!!       REVISION:
!!
!!       OBJECT:   SET UP EMISSIVITY VERSUS FIXED WAVENUMBERS AND SURFACE TYPES
!!
!!         CERES
!!         -----
!!         Emissivity data available at low spectral resolution: only 14 values 
!!         to cover the entire spectrum. Thus, this can be used as a nominal value.
!!         The error associated with this emissivity can roughly be estimated to
!!         increase with lower emissivity as : (1-EMI)*0.5 
!!         (i.e. as large as 0.10 for EMI=0.80 but better than 0.01 for EMI > 0.98)
!!         -No dependence on viewing angle is assumed.
!!         -Not to be used for oceans uncovered by ice.
!!
!!         Longwave Emmissivities in 12 original Fu bands + 2 extra to cover the range
!!         ---------------------------------------------------------------------------
!!         Longwave spectral intervals [cm-1] for the Fu & Liou code:
!!
!!         Band  1          2          3          4          5          6
!!           2200-1900, 1900-1700, 1700-1400, 1400-1250, 1250-1100, 1100-980,
!!         Band  7          8          9         10         11         12
!!            980-800,   800-670,   670-540,  540-400,    400-280,   280-0 
!!
!!         Two additional LW spectral intervals have been added in beyond 2200cm-1.
!!         Band        13              14
!!                  2500-2200       2850-2500
!!
!!         Emissivity    ems(band(1))   from April data, Table2 of Chen et al
!!         11th Conf Sat Met, Madison, WI, p 514
!!          here regoganized as 14 13 1 2 ... 12 above
!!
!!         20 surface types
!!         ----------------
!!          1= evergreen nleaf  2= evergreen bleaf 3= deciduous nleaf  4= deciduous bleaf
!!          5= mixed forests    6= closed shrubs   7= open shrubs      8= woody savanna
!!          9= savanna         10= grasslands     11= perma wet       12= croplands
!!         13= urban           14= mosaic         15= snow            16= barren (deserts)
!!         17= water           18= toundra        19= fresh snow      20= sea ice
!!
!!
!!       ARGUMENTS:
!!          INPUT:
!!            -WAVEN(NCHN)   : WAVENUMBERS (CM-1)
!!            -NCHN          : NUMBER OF BANDS FOR WHICH EMISSIVITY IS NEEDED
!!
!!          OUTPUT:
!!            -EMI_MAT(NCHN,NTYPE) : EMISSIVITY (0.0-1.0)
!--------------------------------------------------------------------------

  subroutine CERES_EMATRIX(EMI_MAT, waven,nchn)
    implicit none
    integer ,intent(in) ::  NCHN
    real (8),intent(in) :: WAVEN(NCHN)
    real (8),intent(out):: EMI_MAT(NCHN,20)
    !*********************************************
    integer            :: I, NC, NT
    real  (8)          :: DUM

    ! CERES bands central wavenumber (covers 3.7 micron to 71.4 mic)
    integer ,parameter :: NB=14
    real  (8)          :: MID(NB)

    ! CERES emissivity per wavenumber and surface types

    real  (8)          :: EMI_TAB(NB,20)

    data MID /                                                   &
         2675.d0, 2350.d0, 2050.d0, 1800.d0, 1550.d0, 1325.d0, 1175.d0, 1040.d0,  &
         890.d0,  735.d0,  605.d0,  470.d0,  340.d0,  140.d0 /

    data EMI_TAB /                                               &
         0.951d0, 0.989d0, 0.989d0, 0.989d0, 0.990d0, 0.991d0, 0.991d0, 0.990d0,  &
         0.990d0, 0.995d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                &
         0.956d0, 0.989d0, 0.989d0, 0.989d0, 0.990d0, 0.991d0, 0.991d0, 0.990d0,  &
         0.990d0, 0.995d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                &
         0.929d0, 0.985d0, 0.985d0, 0.986d0, 0.984d0, 0.983d0, 0.979d0, 0.980d0,  &
         0.973d0, 0.987d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                &
         0.943d0, 0.985d0, 0.985d0, 0.986d0, 0.984d0, 0.983d0, 0.979d0, 0.980d0,  &
         0.973d0, 0.987d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                &
         0.945d0, 0.987d0, 0.987d0, 0.987d0, 0.987d0, 0.987d0, 0.985d0, 0.985d0,  &
         0.982d0, 0.991d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                &
         0.933d0, 0.949d0, 0.949d0, 0.970d0, 0.974d0, 0.971d0, 0.947d0, 0.958d0,  &
         0.966d0, 0.975d0, 0.984d0, 0.984d0, 0.984d0, 0.984d0,                &
         0.873d0, 0.873d0, 0.873d0, 0.934d0, 0.944d0, 0.939d0, 0.873d0, 0.904d0,  &
         0.936d0, 0.942d0, 0.951d0, 0.951d0, 0.951d0, 0.951d0,                &
         0.930d0, 0.987d0, 0.987d0, 0.990d0, 0.992d0, 0.993d0, 0.983d0, 0.975d0,  &
         0.985d0, 0.993d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                &
         0.926d0, 0.987d0, 0.987d0, 0.990d0, 0.992d0, 0.993d0, 0.983d0, 0.975d0,  &
         0.985d0, 0.993d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                &
         0.899d0, 0.987d0, 0.987d0, 0.990d0, 0.992d0, 0.993d0, 0.983d0, 0.975d0,  &
         0.985d0, 0.993d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                &
         0.951d0, 0.983d0, 0.983d0, 0.987d0, 0.987d0, 0.988d0, 0.983d0, 0.981d0,  &
         0.987d0, 0.982d0, 0.986d0, 0.986d0, 0.986d0, 0.986d0,                &
         0.924d0, 0.987d0, 0.987d0, 0.990d0, 0.992d0, 0.993d0, 0.983d0, 0.975d0,  &
         0.985d0, 0.993d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                &
         0.929d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,  &
         1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                &
         0.926d0, 0.987d0, 0.987d0, 0.989d0, 0.989d0, 0.990d0, 0.984d0, 0.980d0,  &
         0.983d0, 0.992d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,                &
         0.972d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0,  &
         1.000d0, 0.999d0, 0.999d0, 0.999d0, 0.999d0, 0.999d0,                &
         0.866d0, 0.835d0, 0.835d0, 0.916d0, 0.934d0, 0.923d0, 0.835d0, 0.877d0,  &
         0.921d0, 0.926d0, 0.934d0, 0.934d0, 0.934d0, 0.934d0,                &
         0.973d0, 0.979d0, 0.979d0, 0.983d0, 0.982d0, 0.982d0, 0.984d0, 0.987d0,  &
         0.989d0, 0.972d0, 0.972d0, 0.972d0, 0.972d0, 0.972d0,                &
         0.968d0, 0.947d0, 0.947d0, 0.967d0, 0.988d0, 0.979d0, 0.975d0, 0.977d0,  &
         0.992d0, 0.989d0, 0.989d0, 0.989d0, 0.989d0, 0.989d0,                &
         0.984d0, 0.988d0, 0.988d0, 0.988d0, 0.988d0, 0.988d0, 0.988d0, 0.988d0,  &
         0.988d0, 0.988d0, 0.988d0, 0.988d0, 0.988d0, 0.988d0,                &
         0.964d0, 0.979d0, 0.979d0, 0.979d0, 0.979d0, 0.979d0, 0.979d0, 0.979d0,  &
         0.979d0, 0.979d0, 0.979d0, 0.979d0, 0.979d0, 0.979d0  /



    do NT = 1, 20
      do NC = 1, NCHN
        if ( WAVEN(NC) > MID(1) ) THEN
          EMI_MAT(NC,NT) = EMI_TAB(1,NT)
        else if ( WAVEN(NC) < MID(NB) ) THEN
          EMI_MAT(NC,NT) = EMI_TAB(NB,NT)
        else
          do I = 1, NB - 1
            if ( WAVEN(NC) <= MID(I) .AND. WAVEN(NC) >= MID(I + 1) ) THEN
              DUM = ( WAVEN(NC) - MID(I) ) / ( MID(I + 1) - MID(I) )
              EMI_MAT(NC,NT) = EMI_TAB(I,NT) + ( EMI_TAB(I + 1,NT) - EMI_TAB(I,NT) ) * DUM
              EXIT
            end if
          end do
        end if
      end do
    end do


  end subroutine CERES_EMATRIX

!--------------------------------------------------------------------------
!!**ID EMI_SEA -- GET OCEAN SURFACE EMISSIVITY
!!
!!       AUTHOR:   L. GARAND                March 1999
!!                       improved with IMEM       2004
!!                 A. BEAULNE (CMDA/SMC)    April 2006  (ADAPT TO 3DVAR)
!!
!!       REVISION:
!!
!!       OBJECT:    GET OCEAN SURFACE EMISSIVITY
!!
!!         Note: 
!!         IMEM(NC), set to zero initially, on next call IMEM will have the
!!         right boundary channel to save search time in interpolation.
!!         IOPT=1 means activate IMEM option (all calls ask for same channels)
!!
!!         To get surface ocean emissivity for a group of channels with
!!         wavenumbers WNUM (cm-1) looking at one point with surface
!!         wind speed WIND from angle ANGLE.
!!         Based on Masuda,1988, Remote Sens. of Envir, 313-329.
!!         Coded emissivity routine based on Masuda's data by Tom Kleespies
!!         Covers 650-2857 cm-1 or 3.1-15.4 microns
!!
!!         CAUTION: extrapolated values from 769-650 cm-1
!!         and interpolated values between 2439-1250 cm-1
!!
!!       ARGUMENTS:
!!          INPUT:
!!            -WNUM(NC)       : CHANNEL WAVENUMBERS (CM-1)
!!            -ANGLE          : VIEWING ANGLE (DEG)
!!            -WIND           : SURFACE WIND SPEED (M/S)
!!            -NP             : NUMBER OF PROFILES
!!            -NC             : NUMBER OF CHANNELS
!!
!!          OUTPUT:
!!            -EM_OC(NC,NP)   : OCEAN EMISSIVITIES (0.-1.)
!--------------------------------------------------------------------------
  subroutine EMI_SEA(EM_OC, wnum,angle,wind,np,nc)
    implicit none
    integer,intent(in) :: NC,NP
    real (8),intent(in)   :: WNUM(NC),ANGLE(NP),WIND(NP)
    real (8),intent(out)  :: EM_OC(NC,NP)
!*******************************************************
    integer      :: I,K,L
    integer      :: IMEM(NC) !,IOPT
    integer      :: MCHAN(2)
    real (8)     :: DUM
    real (8)     :: REFW(19),EMI2(2,NP)


    !* Masuda's 19 wavelengths converted to wavenumber

    data REFW/ 2857.1d0, 2777.7d0, 2702.7d0, 2631.6d0, 2564.1d0,  &
         2500.0d0, 2439.0d0, 1250.0d0, 1190.5d0, 1136.3d0,  &
         1087.0d0, 1041.7d0, 1000.0d0, 952.38d0, 909.09d0,  &
         869.57d0, 833.33d0, 800.00d0, 769.23d0/


    !* IMEM options

!    IOPT = 1
    IMEM(:) = 0

    do I = 1, NC

!       if ( IMEM(I) > 0 .AND. IOPT == 1 ) GO TO 50

      !* out of range
      if ( WNUM(I) < 645.d0 .OR. WNUM(I) > REFW(1) ) THEN
        write(*,'(A,1x,e12.4)') ' fatal: wavenumber out of range in emi_sea', WNUM(I)
        stop
      else if ( WNUM(I) <= REFW(19) .AND. WNUM(I) > 645.d0 ) THEN
        !* extrapolated from 769 cm-1 to 645 cm-1: NOT FROM REAL DATA
        !* nevertheless thought to be much better than unity
        !* this is a region of relatively rapid emissivity change
        !* worst estimates for 700-645 cm-1, but these channels do not
        !* see the surface (strong co2 absorption).
        IMEM(I) = 18
      else
        !* CAUTION interpolation on large interval 1250-2439 cm-1
        !* where no data is available except that of ASTER. ASTER
        !* shows a relatively smooth variation with wavelength except
        !* for a sharp drop at 1600 cm-1 with highs at 1550 and 1650 cm-1
        !* with peak-to-peak variation of 1.5% in that narrow range.
        !* Worst estimates would be between 1400-1800 cm-1 in HIRS ch 12
        !* which only in very cold atmospheres sees the surface.
        do K = 1, 18
          if ( WNUM(I) > REFW(K + 1) .AND. WNUM(I) <= REFW(K) ) THEN
            IMEM(I) = K
          end if
        end do

      end if
   
      MCHAN(1)= IMEM(I)
      MCHAN(2)= IMEM(I) + 1

      DUM = ( WNUM(I) - REFW(MCHAN(1)) ) / ( REFW(MCHAN(2)) - REFW(MCHAN(1)) )

      call COMP_IR_EMISS(EMI2, wind,angle,2,np,mchan)

!* INTERPOLATION/EXTRAPOLATION in wavenumber 

      do L = 1, NP
  
        EM_OC(I,L) = EMI2(1,L) + ( EMI2(2,L) - EMI2(1,L) ) * DUM
          
      end do

    end do


  end subroutine EMI_SEA

!--------------------------------------------------------------------------
!! SUBROUTINE  tvs_rttov_read_coefs
!!
!! *Purpose*: MPI wrapper for rttov_read_coefs
!!            the coefficient files are read by MPI task 0
!!            and then broadcasted to the other tasks according to the selected
!!            channels. Argument channels is mandatory (it is optional in rttov_setup)
!!            optional argument channels_rec was removed (it is useful only in principal component mode)
!!            other optionnal arguments were removed :
!!                            form_coef,         & !to specify format
!!                            form_scaer,        &
!!                            form_sccld,        &
!!                            form_pccoef,       &
!!                            file_coef,         & !to specify filename
!!                            file_scaer,        &
!!                            file_sccld,        &
!!                            file_pccoef,       &
!!                            file_id_coef,      & !to specify fortran unit number
!!                            file_id_scaer,     &
!!                            file_id_sccld,     &
!!                            file_id_pccoef,    &
!!                            path                 ! to specify the path to look for coefficient files
!!            if necessary these arguments could be  added (ask S. Heilliette)
!!            also this subroutine will work only for clear sky radiance computations
!!            if somebody wants to do realistic cloud or aerosol affected radiance simulations
!!            some changes are needed. Ask me in that case. (S. Heilliette) 
!! It is implicitely assumed that the options are the same for all MPI tasks for a given instrument
!! No check will be done (options for task 0 will be used for all tasks). Only differences in channel lists are accounted for.
!! S. Heilliette May 2017
!--------------------------------------------------------------------------

  subroutine tvs_rttov_read_coefs(err, coefs, opts, channels, instrument)
    implicit none

    integer(kind=jpim), intent(out) :: err
    TYPE(rttov_coefs),  intent(out) :: coefs
    TYPE(rttov_options),intent(in)  :: opts
    integer(kind=jpim), intent(in)  :: channels(:)
    integer(kind=jpim), intent(in)  :: instrument(3)

    integer ,allocatable :: listGlobal(:)
    real(8) ,allocatable :: bigArray(:,:,:,:)
    integer :: i, j, k, l, ichan,igas,ierr, countUniqueChannel, indexchan(size(channels)),channelsb(tvs_maxChannelNumber), listAll(tvs_maxChannelNumber)
    logical :: found, associated0
    integer :: nlte_count, nlte_start,isol,isat,nlte_file_nchan
    integer ,allocatable :: nlte_chans(:) 

    Write(*,*) "Entering tvs_rttov_read_coefs"

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
    if (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) &
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

    if (coefs%coef%solarcoef) THEN
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
        if (ASSOCIATED0) deallocate(coefs%coef%thermal(ichan)%gasarray(igas)%coef)
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

    if (coefs % coef % solarcoef) THEN
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
          associated0 = ASSOCIATED(coefs%coef%solar(ichan)%gasarray(igas)%coef)
          if (ASSOCIATED0) deallocate(coefs%coef%solar(ichan)%gasarray(igas)%coef)
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
    call set_fastcoef_level_bounds(coefs%coef, coefs%coef%thermal, thermal = .TRUE._jplm)
    ! If the SOLAR_FAST_COEFFICIENTS section is not present then point the solar coefs to the thermal coefs
    if (coefs%coef%solarcoef) THEN
      call set_fastcoef_level_bounds(coefs%coef, coefs%coef%solar, thermal = .FALSE._jplm)
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
        if (channels(i) >= coefs%coef%nlte_coef%start_chan .AND. &
             channels(i) < coefs%coef%nlte_coef%start_chan + coefs%coef%nlte_coef%nchan) THEN
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
    integer ,Pointer :: array(:)
    integer ,intent(in) :: oldSize, index(:)
    
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
    real(8) ,Pointer :: array(:)
    integer ,intent(in) :: oldSize, index(:)

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
    real(8) ,Pointer :: array(:,:)
    integer ,intent(in) :: oldSize1, oldSize2,index(:)
    
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
    real(8) ,Pointer :: array(:,:,:)
    integer ,intent(in) :: oldSize1, oldSize2,oldSize3,index(:)

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
    complex(kind=8) ,pointer :: array(:)
    integer ,intent(in) :: oldSize, index(:)
    
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
    real(kind=8) ,pointer :: array(:,:)
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
    real(kind=8) ,pointer :: array(:)
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
    integer(kind=4) ,pointer :: array(:)
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
!! *Purpose*: Computation of Jo and the residuals to the tovs observations
!!
!! @author j. halle *cmda/aes  december 14, 2004
!!
!!revision 001  : a. beaulne *cmda/smc  june 2006
!!                  -modifications for AIRS (codtyp 183)
!!
!!revision 002  : r. sarrazin cmda  april 2008
!!                  -modifications for CSR (codtyp 185)
!!
!!revision 003  : s. heilliette
!!                  -modifications for IASI (codtyp 186)
!!
!!revision 004  : s. heilliette
!!                  -modifications for RTTOV-10 (December 2010)
!!
!!revision 005  : s. macpherson
!!                  -modifications for ATMS (codtyp 192)
!!                  -modifications for CrIS (codtyp 193)
!!revision 006  : s. macpherson  nov 2012
!!                  - remove #include "comtovst.cdk"
!!
!--------------------------------------------------------------------------
  subroutine tvs_calc_jo(pjo,llprint,lobsSpaceData,dest_obs)
    implicit none

    real(8) :: pjo
    logical :: llprint
    type(struct_obs) :: lobsSpaceData
    integer, intent(in) :: dest_obs       ! probably set to OBS_OMP or OBS_OMA

    integer :: isens, indxchn, indxtovs
    real*8 dlsum, zdtb, zjon,zgami,zqcarg
    real*8 zjoch  (0:tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    real*8 zavgnrm(0:tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    integer j, i, krtid, nchanperline, indxs, indxe
    integer inobsjo, incanjo
    integer idatyp
    integer ichn, ichOBS_A, jl
    integer inobsch(0:tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    integer lcanjo(tvs_maxChannelNumber)
    integer :: index_header, index_body
    
    real (8) :: x(tvs_maxChannelNumber),y(tvs_maxChannelNumber)
    integer :: list_chan(tvs_maxChannelNumber)
    integer :: count

    if ( llprint ) write(*,*) "Entering tvs_calc_jo subroutine"

    if ( tvs_nobtov == 0) return    ! exit if there are not tovs data

    ! 1.  Computation of (hx - z)/sigma for tovs data only
    !     ------------------------------------------------

    dlsum    = 0.D0
    inobsjo  = 0
    do j = 1, tvs_nsensors
      do i = 0, tvs_maxChannelNumber
        inobsch(i,j) = 0
        zjoch  (i,j) = 0.0D0
        zavgnrm(i,j) = 0.0D0
      end do
    end do

    ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(lobsSpaceData,'TO')
    HEADER: do
      index_header = obs_getHeaderIndex(lobsSpaceData)
      if (index_header < 0) exit HEADER

      ! 1.1  Extract general information for this observation point
      !      ------------------------------------------------------

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(lobsSpaceData,OBS_ITY,index_header)
      if ( .not. tvs_isIdBurpTovs(idatyp) ) cycle HEADER

      indxtovs = tvs_ltovsno(index_header)
      if ( indxtovs == 0 ) then
        write(*,'(A)') ' tvs_calc_jo: error with indxtovs'
        call utl_abort('tvs_calc_jo')
      endif

      isens = tvs_lsensor(indxtovs)

      ! Set the body list
      ! (& start at the beginning of the list)
      call obs_set_current_body_list(lobsSpaceData, index_header)
      count = 0
      BODY: do 
        index_body = obs_getBodyIndex(lobsSpaceData)
        if (index_body < 0 ) then
          if (count > 0 .and. rmat_lnondiagr) then
            call rmat_sqrtRm1(isens,count,x(1:count),y(1:count),list_chan(1:count),indxtovs)
            dlsum =  dlsum + 0.5d0*dot_product(y(1:count),y(1:count))
          end if
          exit BODY
        end if

        ! Only consider if flagged for assimilation
        if ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) /= 1 ) cycle BODY                

        ichn = nint(obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body))
        ichn = max( 0 , min( ichn , tvs_maxChannelNumber + 1))
        ichOBS_A = max(0 , min(ichn , tvs_maxChannelNumber + 1))
        ichn = ichn - tvs_channelOffset(isens)
        indxchn = utl_findArrayIndex(tvs_ichan(:,isens),tvs_nchan(isens),ichn)
        if ( indxchn == 0 ) then
          write(*,'(A)') ' tvs_calc_jo: error with channel number'
          call utl_abort('tvs_calc_jo')
        end if

        zdtb = obs_bodyElem_r(lobsSpaceData,OBS_PRM,index_body) - &
             tvs_radiance (indxtovs) % bt(indxchn)
        if ( tvs_debug ) then
          write(*,'(a,i4,2f8.2,f6.2)') ' ichn,sim,obs,diff= ', &
               ichn,  tvs_radiance (indxtovs) % bt(indxchn), &
               obs_bodyElem_r(lobsSpaceData,OBS_PRM,index_body), -zdtb
        end if
        call obs_bodySet_r(lobsSpaceData,dest_obs,index_body, zdtb)

        ! Comment out the modification of Jobs due to varqc for now, since this is probably
        ! only needed for use of nonlinear obs operator in minimization, which is not yet
        ! functional, but this interferes with doing ensemble of analyses (M. Buehner, Dec. 2013)
        !if (.not. min_lvarqc .or. obs_bodyElem_r(lobsSpaceData,OBS_POB,index_body).eq.0.0d0) then
        dlsum =  dlsum &
             + (obs_bodyElem_r(lobsSpaceData,dest_obs,index_body) * &
             obs_bodyElem_r(lobsSpaceData,dest_obs,index_body)) &
             / (2.D0 * obs_bodyElem_r(lobsSpaceData,OBS_OER,index_body) &
             * obs_bodyElem_r(lobsSpaceData,OBS_OER,index_body))
        !else
        !  compute contribution of data with varqc
 
        !   zgami = obs_bodyElem_r(lobsSpaceData,OBS_POB,index_body)
        !   zjon = (obs_bodyElem_r(lobsSpaceData,dest_obs,index_body)* &
        !           obs_bodyElem_r(lobsSpaceData,dest_obs,index_body))/2.D0
        !   zqcarg = zgami + exp(-1.0D0*zjon)
        !   dlsum= dlsum - log(zqcarg/(zgami+1.D0))
        !end if
        count = count + 1
        x(count) = zdtb
        list_chan(count) = ichn

        inobsjo = inobsjo + 1
        inobsch(ichOBS_A,Isens) = inobsch(ichOBS_A,Isens) + 1
        zjoch(ichOBS_A,Isens)   = &
             zjoch(ichOBS_A,Isens) &
             + obs_bodyElem_r(lobsSpaceData,dest_obs,index_body) * &
             obs_bodyElem_r(lobsSpaceData,dest_obs,index_body) &
             / (obs_bodyElem_r(lobsSpaceData,OBS_OER,index_body) * &
             obs_bodyElem_r(lobsSpaceData,OBS_OER,index_body))
        zavgnrm(ichOBS_A,Isens)   = &
             zavgnrm(ichOBS_A,Isens) - &
             obs_bodyElem_r(lobsSpaceData,dest_obs,index_body) / &
             obs_bodyElem_r(lobsSpaceData,OBS_OER,index_body)
      end do BODY

    end do HEADER

    !   2.  Close up, print summary
    !   .   -----------------------

    pjo = dlsum

    ! printout of mean jo and normalized average for each sensor.

    nchanperline = 18
    if ( llprint .and. inobsjo > 0 ) then
      write(*,*)
      write(*,*)
      write(*,'(10x,A)') "-tvs_calc_jo: computing jo and residuals to tovs  observations"

      do krtid = 1, tvs_nsensors
        do i = 1, tvs_maxChannelNumber
          inobsch(0,krtid) = inobsch(0,krtid) + &
               inobsch(i,krtid)
          zjoch(0,krtid)   = zjoch(0,krtid) + &
               zjoch(i,krtid)
          zavgnrm(0,krtid) = zavgnrm(0,krtid) + &
               zavgnrm(i,krtid)
        end do
      end do

      do jl = 1, tvs_nsensors
        incanjo = 0
        do i = 0, tvs_maxChannelNumber
          if ( inobsch(i,jl) /= 0 ) then
            incanjo = incanjo + 1
            lcanjo(incanjo) = i
          end if
        end do
        if ( incanjo /= 0 ) then
          write(*,'(/1x,"sensor #",i2,". platform: ",a, &
               &   "instrument: ",a)') &
               jl, tvs_satelliteName(jl), tvs_instrumentName(jl)
          do j = 1, incanjo, nchanperline
            indxs = j
            indxe = min(j + nchanperline - 1 , incanjo)
            if ( j == 1 ) then
              write(*,'(1x,"channel",t13,"   all",17i6)') &
                   (lcanjo(i),i=indxs+1,indxe)
            else
              write(*,'(1x,"channel",t13,18i6)') &
                   (lcanjo(i),i=indxs,indxe)
            end if
            write(*,'(1x,"no. obs.",t13,18i6)') &
                 (inobsch(lcanjo(i),jl),i=indxs,indxe)
            write(*,'(1x,"mean jo",t13,18f6.2)') &
                 (zjoch(lcanjo(i),jl)/max(1,inobsch(lcanjo(i),jl)) &
                 ,i=indxs,indxe)
            write(*,'(1x,"norm. bias",t13,18f6.2,/)') &
                 (zavgnrm(lcanjo(i),jl)/max(1,inobsch(lcanjo(i),jl)) &
                 ,i=indxs,indxe)
          end do
        end if
      end do
    end if

end subroutine tvs_calc_jo



end module tovs_nl_mod
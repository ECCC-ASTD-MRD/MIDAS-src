
module tovsLin_mod
  ! MODULE tovsLin_mod (prefix='tvslin' category='5. Observation operators')
  !
  !:Purpose:  Derived types, public variables and procedures related to the 
  !           tangent-linear and adjoint versions of RTTOV
  !
  use rttovInterfaces_mod
  use rttov_types, only :   &
       rttov_profile       ,&
       rttov_profile_cloud ,&
       rttov_radiance      ,&
       rttov_transmission  ,&
       rttov_chanprof      ,&
       rttov_emissivity    ,&
       rttov_options       ,&
       rttov_options_scatt ,&
       rttov_coefs         ,&
       rttov_scatt_coef
  use rttov_const, only : &
      gas_unit_specconc  ,&
      sensor_id_mw       ,&
      surftype_sea       ,&
      errorStatus_success
  use parkind1, only : jpim, jprb, jplm
  use verticalCoord_mod
  use tovsNL_mod
  use utilities_mod
  use MathPhysConstants_mod
  use obsSpaceData_mod
  use columnData_mod
  use midasMpi_mod
 
  implicit none
  save
  private

  public :: tvslin_rttov_tl, tvslin_rttov_ad, tvslin_rttov_k


contains

  !--------------------------------------------------------------------------
  !  tvslin_rttov_tl
  !--------------------------------------------------------------------------
  subroutine tvslin_rttov_tl(columnAnlInc, columnTrlOnAnlIncLev, obsSpaceData)
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
    integer, allocatable :: sensorTovsIndexes(:) 
    integer, allocatable :: sensorHeaderIndexes(:) 
    integer :: allocStatus(5)
    integer :: nobmax
    integer :: sensorIndex, tovsIndex, channelIndex
    integer :: hydroSensorIndex, hydroChannelsCount
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
    integer,external :: omp_get_num_threads
    integer :: nthreads, max_nthreads
    integer :: btIndex, bodyIndex
    integer :: instrum
    integer :: sensorType   !sensor type(1=infrared; 2=microwave; 3=high resolution, 4=polarimetric)
    integer :: errorStatus
    integer, allocatable :: sensorBodyIndexes(:)
    integer, allocatable :: sensorBodyIndexesScatt(:)
    logical, allocatable :: lchannel_subset(:,:)
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
    type(rttov_chanprof), pointer :: chanprof(:)
    type(rttov_chanprof), pointer :: chanprofScatt(:)
    logical, pointer :: calcemis(:)
    logical, pointer :: calcemisScatt(:)
    logical :: runObsOperatorWithClw_tl
    logical :: runObsOperatorWithHydrometeors_tl
    real(8) :: obsOMP
    type (rttov_profile), pointer :: profiles(:)
    type(rttov_profile_cloud), pointer :: cld_profiles(:)
         
    if (tvs_nobtov == 0) return       ! exit if there are not tovs data

    write(*,*) 'tvslin_rttov_tl: Starting'

    call tvs_getProfile(profiles, 'tlad', cld_profiles)

    if (.not. tvs_useO3Climatology .and. .not. col_varExist(columnTrlOnAnlIncLev,'TO3') .and. .not.  col_varExist(columnTrlOnAnlIncLev,'O3L') ) then
      call utl_abort('tvslin_rttov_tl: if tvs_useO3Climatology is set to .true. the ozone variable must be included as an analysis variable in NAMSTATE.')
    else if (.not.tvs_useO3Climatology) then 
      if (col_varExist(columnTrlOnAnlIncLev,'TO3')) then
        ozoneVarName = 'TO3'
      else
        ozoneVarName = 'O3L'
      end if 
    end if

    !  1.  Set index for model's lowest level and model top

    nlv_M = col_getNumLev(columnTrlOnAnlIncLev,'MM')
    nlv_T = col_getNumLev(columnTrlOnAnlIncLev,'TH')

    if ( col_getPressure(columnTrlOnAnlIncLev,1,1,'TH') < col_getPressure(columnTrlOnAnlIncLev,nlv_T,1,'TH') ) then
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

    max_nthreads = mmpi_numThread

    allocStatus(:) = 0
    allocate(sensorTovsIndexes(tvs_nobtov), stat=allocStatus(1))
    call utl_checkAllocationStatus(allocStatus(1:1), ' tvslin_rttov_tl sensorTovsIndexes')
    
    ! 2.  Computation of hx for tovs data only

    
    ! Loop over all sensors specified by user

    sensor_loop:  do sensorIndex = 1, tvs_nsensors

      hydroSensorIndex = tvs_getHydrometeorsIndex(tvs_instruments(sensorIndex))
      hydroChannelsCount = 0
      if ( hydroSensorIndex > 0 ) then
        do channelIndex = 1, tvs_maxNumberOfChannels
          if (tvs_channelsUsingHydrometeors(hydroSensorIndex,channelIndex) > 0) then
            hydroChannelsCount = hydroChannelsCount + 1
          end if
        end do
        if (hydroChannelsCount == 0) then
          call utl_abort('tvslin_rttov_tl: you have to initialize channelsUsingHydrometeors(:,:) in NAMTOV namelist section')
        end if
      end if
      
      runObsOperatorWithClw_tl = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
                                 tvs_isInstrumUsingCLW(tvs_instruments(sensorIndex)) .and. &
                                 tvs_mwInstrumUsingCLW_tl
      runObsOperatorWithHydrometeors_tl = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
                                          col_varExist(columnTrlOnAnlIncLev,'IWCR') .and. &
                                          tvs_isInstrumUsingHydrometeors(tvs_instruments(sensorIndex)) .and. &
                                          tvs_mwInstrumUsingHydrometeors_tl .and. &
                                          hydroChannelsCount > 0
      
      sensorType = tvs_coefs(sensorIndex) % coef % id_sensor
      instrum = tvs_coefs(sensorIndex) % coef % id_inst
      !  loop over all obs.
      profileCount = 0
      do tovsIndex = 1, tvs_nobtov
        !    Currently processed sensor?
        if ( tvs_lsensor(tovsIndex) == sensorIndex ) then
          profileCount = profileCount + 1
          sensorTovsIndexes(profileCount) = tovsIndex
        end if
      end do
     
      if (profileCount == 0) cycle sensor_loop
      nobmax = sensorTovsIndexes(profileCount)
      !     compute the number of calculated radiances for one call
      btCount = tvs_countRadiances(sensorTovsIndexes(1:profileCount), obsSpaceData)
      if (runObsOperatorWithHydrometeors_tl ) then
        btCountScatt = tvs_countRadiancesScatt(sensorTovsIndexes(1:profileCount), obsSpaceData, &
            tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount), sensorIndex)
      else
        btCountScatt = 0
      end if
      btCount = btCount - btCountScatt
      if ( btCount == 0 .and. btCountScatt == 0) cycle  sensor_loop
   
      allocate(sensorHeaderIndexes (profileCount), stat=allocStatus(1))
      allocate(profilesdata_tl(profileCount),      stat=allocStatus(2))
      allocate(cld_profiles_tl(profileCount),      stat=allocStatus(3))
      allocate(surfTypeIsWater(profileCount),      stat=allocStatus(4))
      call utl_checkAllocationStatus(allocStatus(1:4), ' tvslin_rttov_tl')
      if (runObsOperatorWithClw_tl) write(*,*) 'tvslin_rttov_tl: using clw_data'
      if (runObsOperatorWithHydrometeors_tl) write(*,*) 'tvslin_rttov_tl: using hydrometeor data'
 
      sensorHeaderIndexes(:) = 0 
     
      surfTypeIsWater(:) = .false.

      profileCount = 0

      obs_loop: do tovsIndex = 1, nobmax
        if (tvs_lsensor(tovsIndex) /= sensorIndex) cycle obs_loop
        headerIndex = tvs_headerIndex(tovsIndex)
        profileCount = profileCount + 1
        surfTypeIsWater(profileCount) = ( tvs_ChangedStypValue(obsSpaceData,headerIndex) == surftype_sea )
        sensorHeaderIndexes(profileCount) = headerIndex
      end do obs_loop

      call rttov_alloc_prof(            &
          allocStatus(1),               &
          nprofiles=profileCount,       &
          profiles=profilesdata_tl,     &
          nlevels=nlv_T,                &
          opts=tvs_opts(sensorIndex),   &
          asw=1,                        &
          coefs=tvs_coefs(sensorIndex), &
          init=.true.)
     
      call rttov_alloc_scatt_prof ( allocStatus(2),   &
                                    profileCount,     &
                                    cld_profiles_tl,  &
                                    nlv_T,            &
                                    nhydro=5,         &
                                    nhydro_frac=1,    &
                                    asw=1,            &
                                    init=.true.,      &  
                                    flux_conversion=[1,2,0,0,0])
      
      do profileIndex = 1, profileCount
        profilesdata_tl(profileIndex) % gas_units = gas_unit_specconc ! all gas profiles should be provided in kg/kg
        profilesdata_tl(profileIndex) % nlevels   =  nlv_T
        profilesdata_tl(profileIndex) % nlayers   =  nlv_T - 1
        if (tvs_coefs(sensorIndex) % coef % nozone > 0) then
          if (tvs_useO3Climatology) then
            profilesdata_tl(profileIndex) % o3(:) =  0.0d0
          else
            delO3 => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),trim(ozoneVarName))
            profilesdata_tl(profileIndex) % o3(1:nlv_T) =  delO3(1:nlv_T) * 1.0d-9 ! Assumes model ozone in ug/kg
            profilesdata_tl(profileIndex) % s2m % o  = col_getElem(columnAnlInc,ilowlvl_T,sensorHeaderIndexes(profileIndex),trim(ozoneVarName)) * 1.0d-9 ! Assumes model ozone in ug/kg
          end if
        end if

        ! using the zero CLW value for land FOV
        if (runObsOperatorWithClw_tl) then 
          if (surfTypeIsWater(profileIndex)) then
            delCLW => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'LWCR')
            profilesdata_tl(profileIndex) % clw(1:nlv_T)  = delCLW(:)
          else
            profilesdata_tl(profileIndex) % clw(1:nlv_T)  = 0.d0
          end if
        end if

        if (runObsOperatorWithHydrometeors_tl) then 
          if (surfTypeIsWater(profileIndex)) then
            ! rain flux
            if (col_varExist(columnAnlInc,'RF')) then
              delRF => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'RF')
              cld_profiles_tl(profileIndex) % hydro(1:nlv_T,1) = delRF(:)
            else
              cld_profiles_tl(profileIndex) % hydro(1:nlv_T,1) = 0.0d0
            end if

            ! snow flux
            if (col_varExist(columnAnlInc,'SF')) then
              delSF => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'SF')
              cld_profiles_tl(profileIndex) % hydro(1:nlv_T,2)  = delSF(:)
            else
              cld_profiles_tl(profileIndex) % hydro(1:nlv_T,2) = 0.0d0
            end if

            ! graupel
            cld_profiles_tl(profileIndex) % hydro(1:nlv_T,3)  = 0.d0 ! no information for graupel

            ! cloud liquid water content
            delCLW => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'LWCR')
            cld_profiles_tl(profileIndex) % hydro(1:nlv_T,4) = delCLW(:)

            ! cloud ice water content
            delCIW => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'IWCR')
            cld_profiles_tl(profileIndex) % hydro(1:nlv_T,5)  = delCIW(:)
          else
            cld_profiles_tl(profileIndex) % hydro(1:nlv_T,1:5)  = 0.d0
          end if ! surfTypeIsWater

          cld_profiles_tl(profileIndex) % hydro_frac(1:nlv_T,1) = 0.d0   ! no perturbation on cloud fraction as it is a binary variable (or or 1.0) in this implementation
        end if ! runObsOperatorWithHydrometeors_tl
        
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
      end do

      deallocate(sensorHeaderIndexes,  stat=allocStatus(2) )
      deallocate(surfTypeIsWater, stat=allocStatus(3)) 
      call utl_checkAllocationStatus(allocStatus(1:3), ' tvslin_rttov_tl', .false.)

      ! allocate profiledata_tl structures
      if (btCount > 0) then
        call rttov_alloc_tl(                 &
            allocStatus(1),                  &
            asw=1,                           &
            nprofiles=profileCount,          &
            nchanprof=btCount,               &
            nlevels=nlv_T,                   &
            chanprof=chanprof,               &
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
     
        allocate(surfem1(btCount)           ,stat=allocStatus(2))
        allocate(sensorBodyIndexes(btCount) ,stat=allocStatus(3))
        call utl_checkAllocationStatus(allocStatus(1:3), ' tvslin_rtttov_tl check 1')
         !    get Hyperspecral IR emissivities
        if ( tvs_isInstrumHyperSpectral(instrum) ) call tvs_getHIREmissivities(sensorTovsIndexes(1:profileCount), &
            obsSpaceData, surfem1)
        if (btCountScatt > 0) then
          call tvs_getChanprof(sensorTovsIndexes(1:profileCount), obsSpaceData, chanprof, &
              iptobs_cma_opt =sensorBodyIndexes , &
              channelList_opt=tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount), &
              excludeChannelsFromList_opt=.true.)
        else
          call tvs_getChanprof(sensorTovsIndexes(1:profileCount), obsSpaceData, chanprof, &
              iptobs_cma_opt =sensorBodyIndexes)
        end if
        call tvs_getOtherEmissivities(chanprof, sensorTovsIndexes, sensorType, instrum, surfem1, calcemis)
        
        if (sensorType == sensor_id_mw) then
          call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorIndex, chanprof, sensorTovsIndexes(1:profileCount))
        else
          emissivity_local(:) % emis_in = surfem1(:)
        end if
        errorStatus = errorStatus_success
        emissivity_tl(:) % emis_in = 0.0d0
        !  set nthreads to actual number of threads which will be used.

        nthreads = min(max_nthreads, profileCount)  
        call rttov_parallel_tl(                             &
            errorStatus,                                    & ! out
            chanprof,                                       & ! in
            tvs_opts(sensorIndex),                          & ! in
            profiles(sensorTovsIndexes(1:profileCount)),    & ! in
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
          write(*,*) 'temperature           profile=',profiles(sensorTovsIndexes(1)) % t(:)
          write(*,*) 'temperature increment profile=',profilesdata_tl(1) % t(:)
          call utl_abort('tvslin_rttov_tl')
        end if

        !  2.4  Store hx in obsSpaceData,OBS_WORK
      
        do btIndex = 1, btCount
          bodyIndex = sensorBodyIndexes(btIndex)
          call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex, &
              radiancedata_tl % bt(btIndex) )
          if ( tvs_debug ) then
            obsOMP = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
            write(*,'(a,i4,2f8.2)') ' ichn,sim,obs= ', &
                chanprof(btIndex) % chan, radiancedata_tl % bt(btIndex), obsOMP
          end if
        end do
        call rttov_alloc_tl(                 &
           allocStatus(1),                   &
           asw=0,                            &
           nprofiles=profileCount,           &
           nchanprof=btCount,                &
           nlevels=nlv_T,                    &
           chanprof=chanprof,                &
           opts=tvs_opts(sensorIndex),       &
           coefs=tvs_coefs(sensorIndex),     &
           transmission=transmission,        &
           transmission_tl=transmission_tl,  &
           radiance=radiancedata_d,          &
           radiance_tl=radiancedata_tl,      &
           calcemis=calcemis,                &
           emissivity=emissivity_local,      &
           emissivity_tl=emissivity_tl )

        deallocate(surfem1,          stat=allocStatus(2))
        deallocate(sensorBodyIndexes,stat=allocStatus(3))
        call utl_checkAllocationStatus(allocStatus(1:3), ' tvslin_rtttov_tl Z', .false.)
 
      end if
      
      if (btCountScatt >0) then 
        call rttov_alloc_tl(                  &
            allocStatus(1),                   &
            asw=1,                            &
            nprofiles=profileCount,           &
            nchanprof=btCountScatt,           &
            nlevels=nlv_T,                    &
            chanprof=chanprofScatt,           &
            opts=tvs_opts(sensorIndex),       &
            coefs=tvs_coefs(sensorIndex),     &
            radiance=radiancedata_dScatt,     &
            radiance_tl=radiancedata_tlScatt, &
            calcemis=calcemisScatt,           &
            emissivity=emissivity_localScatt, &
            emissivity_tl=emissivity_tlScatt, &
            init=.true.)
        
        
        ! Prepare all input variables required by rttovScatt.
     
        allocate(surfem1Scatt(btCountScatt)           ,stat=allocStatus(2))
        allocate(sensorBodyIndexesScatt(btCountScatt) ,stat=allocStatus(3))
        allocate(frequencies(btCountScatt)            ,stat=allocStatus(4))
        call utl_checkAllocationStatus(allocStatus(1:4), ' tvslin_rtttov_tl check 2')

        allocate(lchannel_subset(profileCount,tvs_nchan(sensorIndex)))
        
        call tvs_getChanprof(sensorTovsIndexes(1:profileCount), obsSpaceData, chanprofScatt, &
            lchannel_subset_opt = lchannel_subset, iptobs_cma_opt =sensorBodyIndexesScatt , &
            channelList_opt=tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount))
        call rttov_scatt_setupindex (     &
            errorStatus,                  &
            profileCount,                 & ! number of profiles
            tvs_nchan(sensorIndex),       & ! number of channels 
            tvs_coefs(sensorIndex),       & ! coef structure read in from rttov coef file
            tvs_coef_scatt(sensorIndex),  & ! coef structure read in from rttov coef file
            btcountScatt,                 & ! number of calculated channels
            chanprofScatt,                & ! channels and profile numbers
            frequencies,                  & ! array, frequency number for each channel
            lchannel_subset )               ! OPTIONAL array of logical flags to indicate a subset of channels
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'tvslin_rttov_tl: fatal error in rttov_scatt_setupindex ', errorStatus
          call utl_abort('tvslin_rttov_tl')
        end if
        deallocate(lchannel_subset)
      
        call tvs_getOtherEmissivities(chanprofScatt, sensorTovsIndexes, sensorType, instrum, surfem1Scatt, calcemisScatt)

        call tvs_getMWemissivityFromAtlas(surfem1Scatt(1:btcountScatt), emissivity_localScatt, sensorIndex, chanprofScatt, sensorTovsIndexes(1:profileCount))
        errorStatus = errorStatus_success
        emissivity_tlScatt(:) % emis_in = 0.0d0
        call rttov_scatt_tl(                                &
            errorStatus,                                    & ! out
            tvs_opts_scatt(sensorIndex),                    & ! in
            nlv_T,                                          & ! in
            chanprofScatt,                                  & ! in
            frequencies,                                    & ! in
            profiles(sensorTovsIndexes(1:profileCount)),    & ! in  
            cld_profiles(sensorTovsIndexes(1:profileCount)),& ! in
            tvs_coefs(sensorIndex),                         & ! in
            tvs_coef_scatt(sensorIndex),                    & ! in
            calcemisScatt,                                  & ! in
            emissivity_localScatt,                          & ! inout
            profilesdata_tl,                                & ! in
            cld_profiles_tl,                                & ! in
            emissivity_tlScatt,                             & ! inout
            radiancedata_dScatt,                            & ! inout
            radiancedata_tlScatt)                             ! inout
        
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'Error in rttov_scatt_tl', errorStatus
          write(*,*) 'temperature           profile=',profiles(sensorTovsIndexes(1)) % t(:)
          write(*,*) 'temperature increment profile=',profilesdata_tl(1) % t(:)
          call utl_abort('tvslin_rttov_tl')
        end if

        !  2.4  Store hx in obsSpaceData,OBS_WORK
      
        do btIndex = 1, btCountScatt
          bodyIndex = sensorBodyIndexesScatt(btIndex)
          call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex, &
              radiancedata_tlScatt % bt(btIndex) )
          if ( tvs_debug ) then
            obsOMP = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
            write(*,'(a,i4,2f8.2)') ' ichn,sim,obs= ', &
                chanprofScatt(btIndex) % chan, radiancedata_tlScatt % bt(btIndex), obsOMP
          end if
        end do
        deallocate (surfem1Scatt           ,stat=allocStatus(1))
        deallocate (sensorBodyIndexesScatt ,stat=allocStatus(2))
        deallocate (frequencies            ,stat=allocStatus(3))
        call rttov_alloc_tl(                  &
            allocStatus(4),                   &
            asw=0,                            &
            nprofiles=profileCount,           &
            nchanprof=btCountScatt,           &
            nlevels=nlv_T,                    &
            chanprof=chanprofScatt,           &
            opts=tvs_opts(sensorIndex),       &
            coefs=tvs_coefs(sensorIndex),     &
            radiance=radiancedata_dScatt,     &
            radiance_tl=radiancedata_tlScatt, &
            calcemis=calcemisScatt,           &
            emissivity=emissivity_localScatt, &
            emissivity_tl=emissivity_tlScatt )
        call utl_checkAllocationStatus(allocStatus(1:4), ' tvslin_rtttov_tl check 3')
      end if
      
      call rttov_alloc_scatt_prof (allocStatus(1),              &
                                   profileCount,                &
                                   cld_profiles_tl,             &
                                   nlv_T,                       &
                                   nhydro=5,                    &
                                   nhydro_frac=1,               &
                                   asw=0,                       &   
                                   flux_conversion=[1,2,0,0,0])
      deallocate(cld_profiles_tl, stat=allocStatus(2))
      
      call rttov_alloc_prof(            &
          allocStatus(3),               &
          nprofiles=profileCount,       &
          profiles=profilesdata_tl,     &
          nlevels=nlv_T,                &
          opts=tvs_opts(sensorIndex),   &
          asw=0,                        &
          coefs=tvs_coefs(sensorIndex), &
          init=.true.)
      
      call utl_checkAllocationStatus(allocStatus(1:3), ' tvslin_rtttov_tl check 4')
    end do sensor_loop

    deallocate (sensorTovsIndexes)
    nullify( profiles )
    write(*,*) 'tvslin_rttov_tl: Finished'

  end subroutine tvslin_rttov_tl

  !--------------------------------------------------------------------------
  !  tvslin_rttov_ad
  !--------------------------------------------------------------------------
  subroutine tvslin_rttov_ad( columnAnlInc, columnTrlOnAnlIncLev, obsSpaceData )
    !
    ! :Purpose: Adjoint of computation of radiance with rttov_ad
    !

    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout)    :: columnAnlInc
    type(struct_columnData), intent(in)       :: columnTrlOnAnlIncLev
    type(struct_obs),        intent(inout)    :: obsSpaceData

    ! Locals:
    type(struct_vco), pointer :: vco_anl
    integer, allocatable :: sensorTovsIndexes(:) 
    integer, allocatable :: sensorHeaderIndexes(:) 
    integer :: allocStatus(17)
    integer :: nthreads
    integer :: nobmax
    integer :: sensorIndex, tovsIndex, channelIndex
    integer :: hydroSensorIndex, hydroChannelsCount
    integer :: ilowlvl_T,ilowlvl_M,profileCount,headerIndex,nlv_M,nlv_T
    integer :: profileIndex, levelIndex
    integer :: Vcode
    real(8), allocatable :: tt_ad(:,:)
    real(8), allocatable :: hu_ad(:,:)
    real(8), allocatable :: pressure_ad(:,:)  
    real(8), allocatable :: ozone_ad(:,:)
    character(len=4) :: ozoneVarName
    real(8), allocatable :: clw_ad(:,:)
    real(8), allocatable :: ciw_ad(:,:), rf_ad(:,:), sf_ad(:,:)
    logical, allocatable :: surfTypeIsWater(:), lchannel_subset(:,:)
    real(8), pointer :: uu_column(:),vv_column(:),tt_column(:),hu_column(:),ps_column(:)
    real(8), pointer :: tg_column(:),p_column(:),o3_column(:),clw_column(:)
    real(8), pointer :: ciw_column(:), rf_column(:),sf_column(:)
    integer :: btCount, btCountScatt
    integer :: max_nthreads
    integer :: instrum
    integer :: btIndex, bodyIndex
    integer :: sensorType   ! sensor type (1=infrared; 2=microwave; 3=high resolution, 4=polarimetric)    
    integer, allocatable :: sensorBodyIndexes(:)
    integer, allocatable :: sensorBodyIndexesScatt(:)
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
    type(rttov_chanprof), pointer :: chanprof(:)
    type(rttov_chanprof), pointer :: chanprofScatt(:)
    logical, pointer :: calcemis  (:)
    logical, pointer :: calcemisScatt  (:)
    logical :: runObsOperatorWithClw_ad
    logical :: runObsOperatorWithHydrometeors_ad
         
    if (tvs_nobtov == 0) return      ! exit if there are not tovs data
    write(*,*) 'tvslin_rttov_ad: Starting'

    call tvs_getProfile(profiles, 'tlad', cld_profiles)

    if (.not. tvs_useO3Climatology .and. .not. col_varExist(columnTrlOnAnlIncLev,'TO3') .and. .not.  col_varExist(columnTrlOnAnlIncLev,'O3L') ) then
      call utl_abort('tvslin_rttov_ad: if tvs_useO3Climatology is set to .true. the ozone variable must be included as an analysis variable in NAMSTATE.')
    else if (.not.tvs_useO3Climatology) then 
      if (col_varExist(columnTrlOnAnlIncLev,'TO3')) then
        ozoneVarName = 'TO3'
      else
        ozoneVarName = 'O3L'
      end if 
    end if

    !     1.    Set index for model's lowest level and model top

    nlv_M = col_getNumLev(columnTrlOnAnlIncLev,'MM')
    nlv_T = col_getNumLev(columnTrlOnAnlIncLev,'TH')

    if (  col_getPressure(columnTrlOnAnlIncLev,1,1,'TH') < col_getPressure(columnTrlOnAnlIncLev,nlv_T,1,'TH') ) then
      ilowlvl_M = nlv_M
      ilowlvl_T = nlv_T
    else
      ilowlvl_M = 1
      ilowlvl_T = 1
    end if

    vco_anl => col_getVco(columnTrlOnAnlIncLev)
    Vcode = vco_anl % Vcode

    !     1.  Get number of threads available and allocate memory for some variables
 
    max_nthreads = mmpi_numThread

    allocate(sensorTovsIndexes(tvs_nobtov))

    !     2.  Computation of adjoint hx for tovs data only

    ! Loop over all sensors specified by user

    sensor_loop:do  sensorIndex = 1, tvs_nsensors

      hydroSensorIndex = tvs_getHydrometeorsIndex(tvs_instruments(sensorIndex))
      hydroChannelsCount = 0
      if ( hydroSensorIndex > 0 ) then
        do channelIndex = 1, tvs_maxNumberOfChannels
          if (tvs_channelsUsingHydrometeors(hydroSensorIndex,channelIndex) > 0) then
            hydroChannelsCount = hydroChannelsCount + 1
          end if
        end do
        if (hydroChannelsCount ==0) then
          call utl_abort('tvslin_rttov_ad: you have to initialize channelsUsingHydrometeors(:,:) in NAMTOV namelist section')
        end if
      end if

      runObsOperatorWithClw_ad = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
                                 tvs_isInstrumUsingCLW(tvs_instruments(sensorIndex)) .and. &
                                 tvs_mwInstrumUsingCLW_tl
      
      runObsOperatorWithHydrometeors_ad = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
                                          col_varExist(columnTrlOnAnlIncLev,'IWCR') .and. &
                                          tvs_isInstrumUsingHydrometeors(tvs_instruments(sensorIndex)) .and. &
                                          tvs_mwInstrumUsingHydrometeors_tl .and. &
                                          hydroChannelsCount > 0

      sensorType = tvs_coefs(sensorIndex) % coef % id_sensor
      instrum = tvs_coefs(sensorIndex) % coef % id_inst

      profileCount = 0
      do tovsIndex = 1, tvs_nobtov
        !    Currently processed sensor?
        if ( tvs_lsensor(tovsIndex) == sensorIndex ) then
          profileCount = profileCount + 1
          sensorTovsIndexes(profileCount) = tovsIndex
        end if
      end do

      if (profileCount == 0) cycle sensor_loop

      nobmax = sensorTovsIndexes(profileCount)

      !  compute the number of radiances/tbs to be calculated
      btCount = tvs_countRadiances(sensorTovsIndexes(1:profileCount), obsSpaceData)
      if (runObsOperatorWithHydrometeors_ad) then
        btCountScatt = tvs_countRadiancesScatt(sensorTovsIndexes(1:profileCount), obsSpaceData, &
            tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount), sensorIndex)
      else
        btCountScatt = 0
      end if
      btCount = btCount - btCountScatt
      
      if (btCount == 0 .and. btCountScatt == 0) cycle sensor_loop
     
      allocStatus(:) = 0
      allocate(sensorHeaderIndexes(profileCount),       stat=allocStatus(1))
      allocate(tt_ad              (nlv_T,profileCount), stat=allocStatus(2))
      allocate(hu_ad              (nlv_T,profileCount), stat=allocStatus(3))
      allocate(pressure_ad        (nlv_T,profileCount), stat=allocStatus(4))
      if (.not. tvs_useO3Climatology) then
        if (tvs_coefs(sensorIndex) % coef % nozone > 0) then
          allocate(ozone_ad(nlv_T,profileCount),   stat=allocStatus(5))
        end if
      end if
      if (runObsOperatorWithClw_ad .or. runObsOperatorWithHydrometeors_ad) then
        allocate (clw_ad(nlv_T,profileCount), stat=allocStatus(6))
      end if
      allocate (surfTypeIsWater(profileCount),stat=allocStatus(7))
      surfTypeIsWater(:) = .false.
      if (runObsOperatorWithHydrometeors_ad) then
        allocate (ciw_ad(nlv_T,profileCount), stat=allocStatus(8))
        allocate (rf_ad(nlv_T,profileCount),  stat=allocStatus(9))
        allocate (sf_ad(nlv_T,profileCount), stat=allocStatus(10))
      end if

      call utl_checkAllocationStatus(allocStatus, ' tvslin_rttov_ad')

      profileCount = 0       
      ! loop over all obs.
      obs_loop: do tovsIndex = 1, nobmax
        if (tvs_lsensor(tovsIndex)/=sensorIndex) cycle obs_loop
        headerIndex = tvs_headerIndex(tovsIndex)
        profileCount = profileCount + 1
        sensorHeaderIndexes(profileCount) = headerIndex
      end do obs_loop
     
      !  2.1  Calculate the actual number of threads which will be used.

      nthreads = min(max_nthreads, profileCount )  
      allocate(profilesdata_ad(profileCount))
      call rttov_alloc_prof(            &
          allocStatus(1),               &
          nprofiles=profileCount,       &
          profiles=profilesdata_ad,     &
          nlevels=nlv_T,                &
          opts=tvs_opts(sensorIndex),   &
          asw=1,                        &
          coefs=tvs_coefs(sensorIndex), &
          init=.true.)
      
      !  2.2  Prepare all input variables required by rttov_ad.
      if (btCount > 0) then
        call rttov_alloc_ad(                 &
            allocStatus(1),                  &
            asw=1,                           &
            nprofiles=profileCount,          &
            nchanprof=btCount,               &
            nlevels=nlv_T,                   &
            chanprof=chanprof,               &
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

        allocate(surfem1(btCount), stat=allocStatus(2))
        allocate(sensorBodyIndexes(btCount), stat=allocStatus(3))
        
        call utl_checkAllocationStatus(allocStatus(1:3), ' tvslin_rttov_ad')
      
        !  get Hyperspectral IR emissivities
      
        if ( tvs_isInstrumHyperSpectral(instrum) ) call tvs_getHIREmissivities(sensorTovsIndexes(1:profileCount), obsSpaceData, surfem1)

        ! Build the list of channels/profiles indices
        if (btCountScatt >0) then
          call tvs_getChanprof(sensorTovsIndexes(1:profileCount), obsSpaceData, chanprof, &
              iptobs_cma_opt =sensorBodyIndexes , &
              channelList_opt=tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount), &
              excludeChannelsFromList_opt=.true.)
        else
          call tvs_getChanprof(sensorTovsIndexes(1:profileCount), obsSpaceData, chanprof, &
              iptobs_cma_opt =sensorBodyIndexes)
        end if
        
        !     get non Hyperspectral IR emissivities
        call tvs_getOtherEmissivities(chanprof, sensorTovsIndexes, sensorType, instrum, surfem1, calcemis)

        if (sensorType == sensor_id_mw) then
          call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorIndex, chanprof, sensorTovsIndexes(1:profileCount))
        else
          emissivity_local(:) % emis_in = surfem1(:)
        end if
        
        do btIndex = 1, btCount
          bodyIndex = sensorBodyIndexes(btIndex)
          radiancedata_ad % bt( btIndex ) = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
        end do
      
        errorStatus = errorStatus_success
        emissivity_ad(:) % emis_in = 0.0d0
        emissivity_ad(:) % emis_out = 0.0d0
        call rttov_parallel_ad(                             &
            errorstatus,                                    &! out
            chanprof,                                       &! in
            tvs_opts(sensorIndex),                          &! in
            profiles(sensorTovsIndexes(1:profileCount)),    &! in
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
          call utl_abort('tvslin_rttov_ad')
        end if

        call rttov_alloc_ad(                 &
            allocStatus(3),                  &
            asw=0,                           &
            nprofiles=profileCount,          &
            nchanprof=btCount,               &
            nlevels=nlv_T,                   &
            chanprof=chanprof,               &
            opts=tvs_opts(sensorIndex),      &
            coefs=tvs_coefs(sensorIndex),    &
            transmission=transmission,       &
            transmission_ad=transmission_ad, &
            radiance=radiancedata_d,         &
            radiance_ad=radiancedata_ad,     &
            calcemis=calcemis,               &
            emissivity=emissivity_local,     &
            emissivity_ad=emissivity_ad )
     
      
        deallocate(surfem1,           stat=allocStatus(4))
        deallocate(sensorBodyIndexes, stat=allocStatus(5))
        call utl_checkAllocationStatus(allocStatus(1:6), ' tvslin_rttov_ad', .false.)
        
      end if

      if (btCountScatt > 0) then
        call rttov_alloc_ad(                   &
            allocStatus(1),                    &
            asw=1,                             &
            nprofiles=profileCount,            &
            nchanprof=btCountScatt,            &
            nlevels=nlv_T,                     &
            chanprof=chanprofScatt,            &
            opts=tvs_opts(sensorIndex),        &
            coefs=tvs_coefs(sensorIndex),      &
            radiance=radiancedata_dScatt,      &
            radiance_ad=radiancedata_adScatt,  &
            calcemis=calcemisScatt,            &
            emissivity=emissivity_localScatt,  &
            emissivity_ad=emissivity_adScatt,  &
            init=.true.)

        allocate(surfem1Scatt(btCountScatt),           stat=allocStatus(2))
        allocate(frequencies(btCountScatt),            stat=allocStatus(3))
        allocate(sensorBodyIndexesScatt(btCountScatt), stat=allocStatus(4))
        allocate(cld_profiles_ad(profileCount),        stat=allocStatus(5))
        call rttov_alloc_scatt_prof (allocStatus(6),   &
                                     profileCount,     &
                                     cld_profiles_ad,  &
                                     nlv_T,            &
                                     nhydro=5,         &
                                     nhydro_frac=1,    &
                                     asw=1,            &     
                                     flux_conversion=[1,2,0,0,0])
        call utl_checkAllocationStatus(allocStatus(1:6), ' tvslin_rttov_ad')

        ! Build the list of channels/profiles indices
        allocate(lchannel_subset(profileCount,tvs_nchan(sensorIndex)))

        call tvs_getChanprof(sensorTovsIndexes(1:profileCount), obsSpaceData, chanprofScatt, &
            lchannel_subset_opt = lchannel_subset, iptobs_cma_opt =sensorBodyIndexesScatt , &
            channelList_opt=tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount))
      
        call rttov_scatt_setupindex (     &
            errorStatus,                  &
            profileCount,                 &  ! number of profiles
            tvs_nchan(sensorIndex),       &  ! number of channels 
            tvs_coefs(sensorIndex),       &  ! coef structure read in from rttov coef file
            tvs_coef_scatt(sensorIndex),  &  ! coef structure read in from rttov coef file
            btcountScatt,                 &  ! number of calculated channels
            chanprofScatt,                &  ! channels and profile numbers
            frequencies,                  &  ! array, frequency number for each channel
            lchannel_subset )                ! OPTIONAL array of logical flags to indicate a subset of channels
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'tvslin_rttov_ad: fatal error in rttov_scatt_setupindex ', errorStatus
          call utl_abort('tvslin_rttov_ad')
        end if
        deallocate(lchannel_subset)
        !     get non Hyperspectral IR emissivities
        call tvs_getOtherEmissivities(chanprofScatt, sensorTovsIndexes, sensorType, instrum, surfem1Scatt, calcemisScatt)

        call tvs_getMWemissivityFromAtlas(surfem1Scatt(1:btcountScatt), emissivity_localScatt, sensorIndex, chanprofScatt, sensorTovsIndexes(1:profileCount))
        
        do btIndex = 1, btCountScatt
          bodyIndex = sensorBodyIndexesScatt(btIndex)
          radiancedata_adScatt % bt( btIndex ) = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
        end do
      
        errorStatus = errorStatus_success
        emissivity_adScatt(:) % emis_in = 0.0d0
        emissivity_adScatt(:) % emis_out = 0.0d0
        !  2.3  Compute ad radiance with rttov_ad
        call rttov_scatt_ad(                                & 
            errorStatus,                                    &! out
            tvs_opts_scatt(sensorIndex),                    &! in
            nlv_T,                                          &! in
            chanprofScatt,                                  &! in
            frequencies,                                    &! in
            profiles(sensorTovsIndexes(1:profileCount)),    &! in
            cld_profiles(sensorTovsIndexes(1:profileCount)),&! in
            tvs_coefs(sensorIndex),                         &! in
            tvs_coef_scatt(sensorIndex),                    &! in
            calcemisScatt,                                  &! in
            emissivity_localScatt,                          &! inout
            profilesdata_ad,                                &! inout
            cld_profiles_ad,                                &! inout
            emissivity_adScatt,                             &! inout
            radiancedata_dScatt,                            &! inout
            radiancedata_adScatt)                            ! inout
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'Error in rttov_scatt_ad', errorStatus
          call utl_abort('tvslin_rttov_ad')
        end if
      
        call rttov_alloc_ad(                  &
            allocStatus(3),                   &
            asw=0,                            &
            nprofiles=profileCount,           &
            nchanprof=btCountScatt,           &
            nlevels=nlv_T,                    &
            chanprof=chanprofScatt,           &
            opts=tvs_opts(sensorIndex),       &
            coefs=tvs_coefs(sensorIndex),     &
            radiance=radiancedata_dScatt,     &
            radiance_ad=radiancedata_adScatt, &
            calcemis=calcemisScatt,           &
            emissivity=emissivity_localScatt, &
            emissivity_ad=emissivity_adScatt )
      
        deallocate (surfem1Scatt,           stat=allocStatus(4))
        deallocate (sensorBodyIndexesScatt, stat=allocStatus(5))
        deallocate (frequencies,            stat=allocStatus(6))
        call utl_checkAllocationStatus(allocStatus(1:6), ' tvslin_rttov_ad', .false.)
      end if

      !   2.0  Store adjoints in columnData object
      tt_ad(:,:) = 0.d0
      hu_ad(:,:) = 0.d0
      pressure_ad(:,:) = 0.d0
      if (.not. tvs_useO3Climatology) then
        if (tvs_coefs(sensorIndex) % coef % nozone > 0) ozone_ad(:,:) = 0.d0
      end if
      if (runObsOperatorWithClw_ad .or. runObsOperatorWithHydrometeors_ad) clw_ad(:,:) = 0.d0
      if (runObsOperatorWithHydrometeors_ad) then
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
 
        if (.not. tvs_useO3Climatology) then
          if (tvs_coefs(sensorIndex) % coef % nozone > 0) then
            ! This step is just to transfer the value for ilowlvl_T to the memory space defined by 'col_getColumn(...trim(ozoneVarName))  
            o3_column => col_getColumn(columnAnlInc,headerIndex,trim(ozoneVarName))
            o3_column(ilowlvl_T) =  profilesdata_ad(profileIndex) % s2m % o * 1.0d-9
            ozone_ad(:,profileIndex) = profilesdata_ad(profileIndex) % o3(:)
          end if
        end if
      
        if (runObsOperatorWithClw_ad) then
          clw_ad(:,profileIndex) = profilesdata_ad(profileIndex) % clw(:)
        end if
      
        if (runObsOperatorWithHydrometeors_ad) then
          rf_ad(:,profileIndex)  = cld_profiles_ad(profileIndex) % hydro(:,1)
          sf_ad(:,profileIndex)  = cld_profiles_ad(profileIndex) % hydro(:,2)
          clw_ad(:,profileIndex) = cld_profiles_ad(profileIndex) % hydro(:,4)
          ciw_ad(:,profileIndex) = cld_profiles_ad(profileIndex) % hydro(:,5)
        end if
      end do
    
      call rttov_alloc_prof(            &
          allocStatus(1),               &
          nprofiles=profileCount,       &
          profiles=profilesdata_ad,     &
          nlevels=nlv_T,                &
          opts=tvs_opts(sensorIndex),   &
          asw=0,                        &
          coefs=tvs_coefs(sensorIndex), &
          init=.true.)
    
      deallocate(profilesdata_ad)
    
      if (btCountScatt > 0) then
        call rttov_alloc_scatt_prof(  &
            allocStatus(2),           &
            profileCount,             &
            cld_profiles_ad,          &
            nlv_T,                    &
            nhydro=5,                 &
            nhydro_frac=1,            &
            asw=0,                    &     
            flux_conversion=[1,2,0,0,0])
        deallocate(cld_profiles_ad, stat=allocStatus(2))
      end if

      !     .  2.1  Store adjoints in columnData object
      !     .       -----------------------------------

      do  profileIndex = 1 , profileCount 
        p_column  => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'P_T')
        tt_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'TT')
        hu_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'HU')
      
        do levelIndex = 1, nlv_T
          p_column(levelIndex) = p_column(levelIndex) + pressure_ad(levelIndex,profileIndex) * MPC_MBAR_PER_PA_R8
          tt_column(levelIndex) = tt_column(levelIndex) + tt_ad(levelIndex,profileIndex)
          hu_column(levelIndex) = hu_column(levelIndex) + hu_ad(levelIndex,profileIndex)
        end do
      end do

      if (.not. tvs_useO3Climatology) then
        if (tvs_coefs(sensorIndex) % coef % nozone > 0) then
          do  profileIndex = 1, profileCount 
            o3_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex),trim(ozoneVarName))
            do levelIndex = 1, col_getNumLev(columnAnlInc,'TH')
              o3_column(levelIndex) = o3_column(levelIndex) +  ozone_ad(levelIndex,profileIndex) * 1.0d-9
            end do
          end do
        end if
      end if

      if (runObsOperatorWithClw_ad) then
        do  profileIndex = 1 , profileCount 
          surfTypeIsWater(profileIndex) = ( tvs_ChangedStypValue(obsSpaceData,sensorHeaderIndexes(profileIndex)) == surftype_sea )
          if (surfTypeIsWater(profileIndex)) then
            clw_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex),'LWCR')
            do levelIndex = 1, col_getNumLev(columnAnlInc,'TH')
              clw_column(levelIndex) = clw_column(levelIndex) + &
                  clw_ad(levelIndex,profileIndex)
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
              rf_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex),'RF')
              do levelIndex = 1, col_getNumLev(columnAnlInc,'TH')
                rf_column(levelIndex) = rf_column(levelIndex) + rf_ad(levelIndex,profileIndex)
              end do
            end if

            ! snow flux
            if (col_varExist(columnAnlInc,'SF')) then
              sf_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex),'SF')
              do levelIndex = 1, col_getNumLev(columnAnlInc,'TH')
                sf_column(levelIndex) = sf_column(levelIndex) + sf_ad(levelIndex,profileIndex)
              end do
            end if

            ! cloud liquid/ice water content
            clw_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex),'LWCR')
            ciw_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex),'IWCR')
            do levelIndex = 1, col_getNumLev(columnAnlInc,'TH')
              clw_column(levelIndex) = clw_column(levelIndex) + clw_ad(levelIndex,profileIndex)
              ciw_column(levelIndex) = ciw_column(levelIndex) + ciw_ad(levelIndex,profileIndex)
            end do
          end if ! surfTypeIsWater
        end do ! profileIndex
      end if ! runObsOperatorWithHydrometeors_ad

      deallocate(sensorHeaderIndexes, stat=allocStatus(1))
      deallocate(tt_ad,               stat=allocStatus(2))
      deallocate(hu_ad,               stat=allocStatus(3))
      deallocate(pressure_ad,         stat=allocStatus(4))
      if (.not. tvs_useO3Climatology) then
        if (tvs_coefs(sensorIndex) % coef % nozone > 0) then
          deallocate(ozone_ad,        stat=allocStatus(5))
        end if
      end if
      if ( allocated(clw_ad) ) then
        deallocate(clw_ad,stat=allocStatus(6))
      end if
      deallocate(surfTypeIsWater,stat=allocStatus(7))
      if ( allocated(ciw_ad) ) then
        deallocate(ciw_ad, stat=allocStatus(8))
        deallocate(rf_ad,  stat=allocStatus(9))
        deallocate(sf_ad,  stat=allocStatus(10))
      end if
      
      call utl_checkAllocationStatus(allocStatus, ' tvslin_fill_profiles_ad', .false.)
      
    end do sensor_loop

    ! 3.  Close up

    deallocate(sensorTovsIndexes)
    nullify(profiles)
    write(*,*) 'tvslin_rttov_ad: Finished'

  end subroutine tvslin_rttov_ad

  !--------------------------------------------------------------------------
  !  tvslin_rttov_K
  !--------------------------------------------------------------------------
  subroutine tvslin_rttov_k(columnTrlOnAnlIncLev, obsSpaceData)
    !
    ! :Purpose: Compute the Jacobian of radiance using rttov_k
    !
    implicit none
  
    ! Arguments:
    type(struct_obs),        intent(inout) :: obsSpaceData         ! obsSpaceData structure
    type(struct_columnData), intent(in)    :: columnTrlOnAnlIncLev ! column structure for background profile

    ! Locals:
    type(rttov_emissivity), pointer :: emissivity_local(:)
    type(rttov_emissivity), pointer :: emissivity_localScatt(:)
    type(rttov_emissivity), pointer :: emissivity_k(:)
    type(rttov_emissivity), pointer :: emissivity_kScatt(:)
    type(rttov_chanprof), pointer   :: chanprof(:)
    type(rttov_chanprof), pointer   :: chanprofScatt(:)
    type (rttov_profile), pointer   :: profiles(:)
    type (rttov_profile), pointer   :: profiles_k(:)
    type(rttov_profile_cloud), pointer :: cld_profiles(:)
    type(rttov_transmission)        :: transmission
    type(rttov_transmission)        :: transmission_k
    type(rttov_radiance)            :: radiancedata_d 
    type(rttov_radiance)            :: radiancedata_k
    type(rttov_radiance)            :: radiancedata_dScatt 
    type(rttov_radiance)            :: radiancedata_kScatt
    integer, allocatable            :: sensorTovsIndexes(:) 
    integer                         :: allocStatus(4)
    integer                         :: nobmax, profileCount, btCount, btCountScatt
    integer                         :: sensorIndex, tovsIndex, channelIndex
    integer                         :: nlv_T
    integer                         :: instrum, asw
    integer                         :: nthreads, max_nthreads
    integer                         :: sensorType
    integer                         :: hydroSensorIndex, hydroChannelsCount
    real(8), allocatable            :: surfem1(:)
    real(8), allocatable            :: surfem1Scatt(:)
    integer, allocatable            :: sensorBodyIndexes(:)
    integer, allocatable            :: sensorBodyIndexesScatt(:)
    integer, allocatable            :: frequencies(:)
    logical, allocatable            :: lchannel_subset(:,:)
    integer                         :: errorstatus
    logical                         :: runObsOperatorWithHydrometeors_k
    logical, pointer                :: calcemis(:)
    logical, pointer                :: calcemisScatt(:)

    if (tvs_nobtov == 0) return ! exit if there are not tovs data
  
    call tvs_getProfile(profiles, 'nl', cld_profiles)
  
    ! Set index for model's lowest level and model top
  
    nlv_T = col_getNumLev(columnTrlOnAnlIncLev,'TH')
  
    max_nthreads = mmpi_numThread 

    allocStatus(:) = 0
    allocate(sensorTovsIndexes(tvs_nobtov), stat=allocStatus(1))
    call utl_checkAllocationStatus(allocStatus(1:1), " tvslin_rttov_tl sensorTovsIndexes")
        
    ! Loop over all sensors specified by user
  
    sensor_loop: do sensorIndex = 1, tvs_nsensors
      
      hydroSensorIndex = tvs_getHydrometeorsIndex(tvs_instruments(sensorIndex))
      hydroChannelsCount = 0
      if ( hydroSensorIndex > 0 ) then
        do channelIndex = 1, tvs_maxNumberOfChannels
          if (tvs_channelsUsingHydrometeors(hydroSensorIndex,channelIndex) > 0) then
            hydroChannelsCount = hydroChannelsCount + 1
          end if
        end do
        if (hydroChannelsCount ==0) then
          call utl_abort('tvslin_rttov_k: you have to initialize channelsUsingHydrometeors(:,:) in NAMTOV namelist section')
        end if
      end if
      
      runObsOperatorWithHydrometeors_k = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
                                         col_varExist(columnTrlOnAnlIncLev,'IWCR') .and. &
                                         tvs_isInstrumUsingHydrometeors(tvs_instruments(sensorIndex)) .and. &
                                         hydroChannelsCount > 0
 
      sensorType = tvs_coefs(sensorIndex) % coef % id_sensor
      instrum = tvs_coefs(sensorIndex) % coef % id_inst

      ! Loop over all obs.
      profileCount = 0
      do tovsIndex = 1, tvs_nobtov
        if (tvs_lsensor(tovsIndex) == sensorIndex) then
          profileCount = profileCount + 1
          sensorTovsIndexes(profileCount) = tovsIndex
        end if
      end do
      
      if (profileCount == 0) cycle sensor_loop
      nobmax = sensorTovsIndexes(profileCount)

      ! Compute the number of calculated radiances for one call
      btCount = tvs_countRadiances(sensorTovsIndexes(1:profileCount), obsSpaceData)
      if (runObsOperatorWithHydrometeors_k) then
        btCountScatt = tvs_countRadiancesScatt(sensorTovsIndexes(1:profileCount), obsSpaceData, &
            tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount), sensorIndex)
      else
        btCountScatt = 0
      end if
      btCount = btCount - btCountScatt
      
      if (btCount == 0 .and. btCountScatt == 0) cycle sensor_loop
      
      if (btCount > 0) then
        call rttov_alloc_k(                  &
            allocStatus(1),                  &
            asw=1,                           & ! to allocate
            nprofiles=profileCount,          &
            nchanprof=btCount,               &
            nlevels=nlv_T,                   &
            chanprof=chanprof,               &
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
  
        call utl_checkAllocationStatus(allocStatus(1:1), " tovs_rtttov_k rttov_alloc_k 1")
          
        ! Set nthreads to actual number of threads which will be used.
  
        nthreads = min(max_nthreads, profileCount)  
  
        ! Prepare all input variables required by rttov.
        
        allocate(surfem1(btCount),           stat=allocStatus(1))
        allocate(sensorBodyIndexes(btCount), stat=allocStatus(2))
        call utl_checkAllocationStatus(allocStatus(1:2), " tovs_rtttov_k")
        
        write(*,*) 'tvslin_rttov_k: Get surface emissiviy'

        !    get Hyperspecral IR emissivities
        if (tvs_isInstrumHyperSpectral(instrum)) call tvs_getHIREmissivities(sensorTovsIndexes(1:profileCount), &
                                                                           obsSpaceData, surfem1)
  
        call tvs_getChanprof(sensorTovsIndexes(1:profileCount), obsSpaceData, chanprof, &
            iptobs_cma_opt = sensorBodyIndexes)
  
        call tvs_getOtherEmissivities(chanprof, sensorTovsIndexes, sensorType, instrum, surfem1, calcemis)
  
        if (sensorType == sensor_id_mw) then
          call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorIndex, chanprof, sensorTovsIndexes(1:profileCount))
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
            chanprof,                                        & ! in
            tvs_opts(sensorIndex),                           & ! in
            profiles(sensorTovsIndexes(1:profileCount)),     & ! in
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
          write(*,*) 'temperature profile=', profiles(sensorTovsIndexes(1)) % t(:)
          write(*,*) 'temperature Jacobian profile=', profiles_k(1) % t(:)
          call utl_abort('tovs_rttov_k')
        end if
  
        ! Write Jacobian to ASCII files
        call tvs_writeJacobianAscii(profiles_k, emissivity_k, profiles, chanProf(1:btCount), &
            obsSpaceData, tvs_satelliteName(sensorIndex), tvs_instrumentName(sensorIndex), &
            sensorBodyIndexes, sensorTovsIndexes, btCount)
  
        ! deallocate profiledata structures
        asw = 0 ! 1 to allocate
        call rttov_alloc_k(                  &
            allocStatus(1),                  &
            asw,                             &
            nprofiles=profileCount,          &
            nchanprof=btCount,               &
            nlevels=nlv_T,                   &
            chanprof=chanprof,               &
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
        
        deallocate(surfem1, stat=allocStatus(2))
        deallocate(sensorBodyIndexes, stat=allocStatus(3))
        call utl_checkAllocationStatus(allocStatus(1:3), "tvslin_rtttov_k", .false.)
      end if
      
      if (btCountScatt > 0) then
        call rttov_alloc_k(                   &
            allocStatus(1),                   &
            asw=1,                            & ! to allocate
            nprofiles=profileCount,           &
            nchanprof=btCountScatt,           &
            nlevels=nlv_T,                    &
            chanprof=chanprofScatt,           &
            opts=tvs_opts(sensorIndex),       &
            profiles_k=profiles_k,            &
            coefs=tvs_coefs(sensorIndex),     &
            radiance=radiancedata_dScatt,     &
            radiance_k=radiancedata_kScatt,   &
            calcemis=calcemisScatt,           &
            emissivity=emissivity_localScatt, &
            emissivity_k=emissivity_kScatt,   &
            init=.true.)
        ! Prepare all input variables required by rttovScatt.
     
        allocate(surfem1Scatt(btCountScatt),           stat=allocStatus(2))
        allocate(sensorBodyIndexesScatt(btCountScatt), stat=allocStatus(3))
        allocate(frequencies(btCountScatt),            stat=allocStatus(4))
        call utl_checkAllocationStatus(allocStatus(1:4), ' tvslin_rtttov_tl check 2')

        allocate(lchannel_subset(profileCount,tvs_nchan(sensorIndex)))
        
        call tvs_getChanprof(sensorTovsIndexes(1:profileCount), obsSpaceData, chanprofScatt, &
            lchannel_subset_opt = lchannel_subset, iptobs_cma_opt =sensorBodyIndexesScatt , &
            channelList_opt=tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount))
        call rttov_scatt_setupindex (     &
            errorStatus,                  &
            profileCount,                 & ! number of profiles
            tvs_nchan(sensorIndex),       & ! number of channels 
            tvs_coefs(sensorIndex),       & ! coef structure read in from rttov coef file
            tvs_coef_scatt(sensorIndex),  & ! coef structure read in from rttov coef file
            btcountScatt,                 & ! number of calculated channels
            chanprofScatt,                & ! channels and profile numbers
            frequencies,                  & ! array, frequency number for each channel
            lchannel_subset )               ! OPTIONAL array of logical flags to indicate a subset of channels
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'tvslin_rttov_tl: fatal error in rttov_scatt_setupindex ', errorStatus
          call utl_abort('tvslin_rttov_tl')
        end if
        deallocate(lchannel_subset)
          
        ! Prepare all input variables required by rttov.
        
        allocate(surfem1Scatt(btCountScatt),           stat=allocStatus(1))
        allocate(sensorBodyIndexesScatt(btCountScatt), stat=allocStatus(2))
        call utl_checkAllocationStatus(allocStatus(1:2), " tovs_rtttov_k")
        
        write(*,*) 'tvslin_rttov_k: Get surface emissivity for rttovScatt_k'
        call tvs_getOtherEmissivities(chanprofScatt, sensorTovsIndexes, sensorType, instrum, surfem1Scatt, calcemisScatt)
  
        call tvs_getMWemissivityFromAtlas(surfem1Scatt(1:btcountScatt), emissivity_localScatt, sensorIndex, chanprofScatt, sensorTovsIndexes(1:profileCount))
      
        ! Initialize variables
        emissivity_k % emis_in = 0.0
        emissivity_k % emis_out = 0.0
  
        call rttov_init_prof(profiles_k)
        call rttov_init_rad(radiancedata_k)

        radiancedata_k % bt = 1.0
        radiancedata_k % total(:) = 1.0
        
        call  tvslin_rttovScattK(                              &
              errorStatus,                                     &
              chanprofScatt,                                   &
              frequencies,                                     &
              tvs_opts(sensorIndex),                           &
              tvs_opts_scatt(sensorIndex),                     &
              profiles(sensorTovsIndexes(1:profileCount)),     &
              profiles_k,                                      &
              cld_profiles(sensorTovsIndexes(1:profileCount)), &
              tvs_coefs(sensorIndex),                          &
              tvs_coef_scatt(sensorIndex),                     & 
              radiancedata_dScatt,                             &
              radiancedata_kScatt,                             &
              calcemisScatt,                                   &
              emissivity_localScatt,                           &
              emissivity_k       )
                   
        if (errorstatus /= errorStatus_success) then
          write(*,*) "Error in tvslin_rttovScattK", errorstatus
          write(*,*) 'temperature profile=', profiles(sensorTovsIndexes(1)) % t(:)
          write(*,*) 'temperature Jacobian profile=', profiles_k(1) % t(:)
          call utl_abort('tvslin_rttovScattK')
        end if
  
        ! Write Jacobian to ASCII files
        call tvs_writeJacobianAscii(profiles_k, emissivity_k, profiles, chanProfScatt(1:btCountScatt), &
            obsSpaceData, tvs_satelliteName(sensorIndex), tvs_instrumentName(sensorIndex), &
            sensorBodyIndexesScatt, sensorTovsIndexes, btCountScatt)
  
        ! deallocate RTTOV structures
        call rttov_alloc_k(                   &
            allocStatus(1),                   &
            asw=0,                            &
            nprofiles=profileCount,           &
            nchanprof=btCountScatt,           &
            nlevels=nlv_T,                    &
            chanprof=chanprofScatt,           &
            opts=tvs_opts(sensorIndex),       &
            profiles_k=profiles_k,            &
            coefs=tvs_coefs(sensorIndex),     &
            radiance=radiancedata_dScatt,     &
            radiance_k=radiancedata_kScatt,   &
            calcemis=calcemisScatt,           &
            emissivity=emissivity_localScatt, &
            emissivity_k=emissivity_kScatt)
        
        deallocate(surfem1Scatt, stat=allocStatus(2))
        deallocate(sensorBodyIndexesScatt, stat=allocStatus(3))
        call utl_checkAllocationStatus(allocStatus(1:3), "tvslin_rtttov_k", .false.)
      end if
      
    end do sensor_loop
  
    deallocate (sensorTovsIndexes)
    nullify(profiles)

    write(*,*) 'tvslin_rttov_k: finished'
  end subroutine tvslin_rttov_k


  !--------------------------------------------------------------------------
  !  tvslin_rttovScattK
  !--------------------------------------------------------------------------
  subroutine tvslin_rttovScattK(&
              errorstatus,      &
              chanprof,         &
              frequencies,      &
              opts,             &
              opts_scatt,       &
              profiles,         &
              profiles_k,       &
              cld_profiles,     &
              coefs,            &
              coef_scatt,       &
              radiance,         &
              radiance_k,       &
              calcemis,         &
              emissivity,       &
              emissivity_k       )
    !
    ! :Purpose: to compute RttovScatt Jacobians from repeated calls to rttov_scatt_tl
    !
    implicit none

    ! Arguments:
    integer(kind=jpim),        intent(out)        :: errorstatus
    type(rttov_chanprof),      intent(in)         :: chanprof(:)
    integer,                   intent(in)         :: frequencies(:)
    type(rttov_options),       intent(in)         :: opts
    type(rttov_options_scatt), intent(in)         :: opts_scatt
    type(rttov_profile),       intent(in)         :: profiles(:)
    type(rttov_profile),       intent(inout)      :: profiles_k(:)
    type(rttov_profile_cloud), intent(in)         :: cld_profiles(:)
    type(rttov_coefs),         intent(in), target :: coefs
    type(rttov_scatt_coef),    intent(in), target :: coef_scatt
    type(rttov_radiance),      intent(inout)      :: radiance
    type(rttov_radiance),      intent(inout)      :: radiance_k
    logical(kind=jplm),        intent(in)         :: calcemis(:)
    type(rttov_emissivity),    intent(inout)      :: emissivity(:)
    type(rttov_emissivity),    intent(inout)      :: emissivity_k(:)
    ! Locals:
    integer :: allocStatus(3)
    integer :: levelIndex, nlv_T, profileIndex, nProfiles
    integer :: btIndex, btCount
    type(rttov_profile), allocatable :: profiles_tl(:) ! tl profiles buffer used in rttov calls
    type(rttov_emissivity) :: emissivity_tl(size(chanprof))
    type(rttov_radiance) :: radiancedata_tl  ! tl radiances full structure buffer used in rttov calls
    type(rttov_profile_cloud), allocatable :: cld_profiles_tl(:) !tl profiles buffer used in RttovScatt calls
    
    nprofiles = size(profiles)
    btCount = size(chanProf)
    nlv_T = profiles(1) % nlevels
    
    call rttov_alloc_prof(    &
        allocStatus(1),       &
        nprofiles=nprofiles,  &
        profiles=profiles_tl, &
        nlevels=nlv_T,        &
        opts=opts,            &
        asw=1,                &
        coefs=coefs,          &
        init=.true.)

    call rttov_alloc_rad(         &
        allocStatus(2),           &
        nchanprof=btcount,        &
        radiance=radiancedata_tl, &
        nlevels=nlv_T,            &
        asw=1,                    &
        init=.true.)
     
    call rttov_alloc_scatt_prof(&
        allocStatus(3),         &
        nprofiles,              &
        cld_profiles_tl,        &
        nlv_T,                  &
        nhydro=5,               &
        nhydro_frac=1,          &
        asw=1,                  &
        init=.true.,            &  
        flux_conversion=[1,2,0,0,0])
    
    call utl_checkAllocationStatus(allocStatus, ' tvslin_rttovScattK')

    emissivity_tl(:) % emis_in = 0.d0
    emissivity_tl(:) % emis_out = 0.d0

    !Temperature Jacobian
    do levelIndex = 1, nlv_T
      do profileIndex = 1, nProfiles
        profiles_tl(profileIndex)%t(levelIndex) = 1.d0
      end do
      call rttov_scatt_tl(  &
            errorStatus,    & ! out
            opts_scatt,     & ! in
            nlv_T,          & ! in
            chanprof,       & ! in
            frequencies,    & ! in
            profiles,       & ! in  
            cld_profiles,   & ! in
            coefs,          & ! in
            coef_scatt,     & ! in
            calcemis,       & ! in
            emissivity,     & ! inout
            profiles_tl,    & ! in
            cld_profiles_tl,& ! in
            emissivity_tl,  & ! inout
            radiance,       & ! inout
            radiancedata_tl)  ! inout

      do btIndex = 1, btCount
        profiles_k(btIndex)%t(levelIndex) = radiancedata_tl%total(btIndex)
      end do
      do profileIndex = 1, nProfiles
        profiles_tl(profileIndex)%t(levelIndex) = 0d0
      end do
    end do
    
    !Then humidity Jacobians
    do levelIndex = 1, nlv_T
      do profileIndex = 1, nProfiles
        profiles_tl(profileIndex)%q(levelIndex) = 1.d0
      end do
      call rttov_scatt_tl(  &
            errorStatus,    & ! out
            opts_scatt,     & ! in
            nlv_T,          & ! in
            chanprof,       & ! in
            frequencies,    & ! in
            profiles,       & ! in  
            cld_profiles,   & ! in
            coefs,          & ! in
            coef_scatt,     & ! in
            calcemis,       & ! in
            emissivity,     & ! inout
            profiles_tl,    & ! in
            cld_profiles_tl,& ! in
            emissivity_tl,  & ! inout
            radiance,       & ! inout
            radiancedata_tl)  ! inout

      do btIndex = 1, btCount
        profiles_k(btIndex)%q(levelIndex) = radiancedata_tl%total(btIndex)
      end do
      do profileIndex = 1, nProfiles
        profiles_tl(profileIndex)%q(levelIndex) = 0d0
      end do
    end do
    
    !Pressure Jacobians
    do levelIndex = 1, nlv_T
      do profileIndex = 1, nProfiles
        profiles_tl(profileIndex)%p(levelIndex) = 1.d0
      end do
      call rttov_scatt_tl(  &
            errorStatus,    & ! out
            opts_scatt,     & ! in
            nlv_T,          & ! in
            chanprof,       & ! in
            frequencies,    & ! in
            profiles,       & ! in  
            cld_profiles,   & ! in
            coefs,          & ! in
            coef_scatt,     & ! in
            calcemis,       & ! in
            emissivity,     & ! inout
            profiles_tl,    & ! in
            cld_profiles_tl,& ! in
            emissivity_tl,  & ! inout
            radiance,       & ! inout
            radiancedata_tl)  ! inout

      do btIndex = 1, btCount
        profiles_k(btIndex) % p(levelIndex) = radiancedata_tl%total(btIndex)
      end do
      do profileIndex = 1, nProfiles
        profiles_tl(profileIndex) % p(levelIndex) = 0d0
      end do
    end do

    !emmissivity Jacobians
    
    do btIndex = 1, btCount
      emissivity_tl(btIndex) % emis_in = 1.d0
    end do
    call rttov_scatt_tl(  &
          errorStatus,    & ! out
          opts_scatt,     & ! in
          nlv_T,          & ! in
          chanprof,       & ! in
          frequencies,    & ! in
          profiles,       & ! in  
          cld_profiles,   & ! in
          coefs,          & ! in
          coef_scatt,     & ! in
          calcemis,       & ! in
          emissivity,     & ! inout
          profiles_tl,    & ! in
          cld_profiles_tl,& ! in
          emissivity_tl,  & ! inout
          radiance,       & ! inout
          radiancedata_tl)  ! inout
    do btIndex = 1, btCount
      emissivity_k(btIndex) % emis_out = radiancedata_tl%total(btIndex)
      emissivity_tl(btIndex) % emis_in = 0.d0
    end do

    !Skin temperature Jacobians
    do profileIndex = 1, nProfiles
      profiles_tl(profileIndex) % skin % t = 1.d0
    end do
    call rttov_scatt_tl(  &
          errorStatus,    & ! out
          opts_scatt,     & ! in
          nlv_T,          & ! in
          chanprof,       & ! in
          frequencies,    & ! in
          profiles,       & ! in  
          cld_profiles,   & ! in
          coefs,          & ! in
          coef_scatt,     & ! in
          calcemis,       & ! in
          emissivity,     & ! inout
          profiles_tl,    & ! in
          cld_profiles_tl,& ! in
          emissivity_tl,  & ! inout
          radiance,       & ! inout
          radiancedata_tl)  ! inout
    do btIndex = 1, btCount
      profiles_k(btIndex) % skin % t = radiancedata_tl%total(btIndex)
    end do
    do profileIndex = 1, nProfiles
      profiles_tl(profileIndex) % skin % t = 0.d0
    end do

    !2m air temperature Jacobians
    do profileIndex = 1, nProfiles
      profiles_tl(profileIndex) % s2m % t = 1.d0
    end do
    call rttov_scatt_tl(  &
          errorStatus,    & ! out
          opts_scatt,     & ! in
          nlv_T,          & ! in
          chanprof,       & ! in
          frequencies,    & ! in
          profiles,       & ! in  
          cld_profiles,   & ! in
          coefs,          & ! in
          coef_scatt,     & ! in
          calcemis,       & ! in
          emissivity,     & ! inout
          profiles_tl,    & ! in
          cld_profiles_tl,& ! in
          emissivity_tl,  & ! inout
          radiance,       & ! inout
          radiancedata_tl)  ! inout
    do btIndex = 1, btCount
      profiles_k(btIndex) % s2m % t = radiancedata_tl%total(btIndex)
    end do
    do profileIndex = 1, nProfiles
      profiles_tl(profileIndex) % s2m % t = 0.d0
    end do

    !surface pressure Jacobians
    do profileIndex = 1, nProfiles
      profiles_tl(profileIndex) % s2m % p = 1.d0
    end do
    call rttov_scatt_tl(  &
          errorStatus,    & ! out
          opts_scatt,     & ! in
          nlv_T,          & ! in
          chanprof,       & ! in
          frequencies,    & ! in
          profiles,       & ! in  
          cld_profiles,   & ! in
          coefs,          & ! in
          coef_scatt,     & ! in
          calcemis,       & ! in
          emissivity,     & ! inout
          profiles_tl,    & ! in
          cld_profiles_tl,& ! in
          emissivity_tl,  & ! inout
          radiance,       & ! inout
          radiancedata_tl)  ! inout
    do btIndex = 1, btCount
      profiles_k(btIndex) % s2m % p = radiancedata_tl%total(btIndex)
    end do
    do profileIndex = 1, nProfiles
      profiles_tl(profileIndex) % s2m % p = 0.d0
    end do
    
    ! mempory deallocation
    call rttov_alloc_prof(    &
        allocStatus(1),       &
        nprofiles=nprofiles,  &
        profiles=profiles_tl, &
        nlevels=nlv_T,        &
        opts=opts,            &
        asw=0,                &
        coefs=coefs)
    
    call rttov_alloc_rad(         &
        allocStatus(2),           &
        nchanprof=btcount,        &
        radiance=radiancedata_tl, &
        nlevels=nlv_T,            &
        asw=0                     )
     
    call rttov_alloc_scatt_prof(&
        allocStatus(3),         &
        nprofiles,              &
        cld_profiles_tl,        &
        nlv_T,                  &
        nhydro=5,               &
        nhydro_frac=1,          &
        asw=0                     )
    
    call utl_checkAllocationStatus(allocStatus, ' tvslin_rttovScattK',.false.)
    
  end subroutine tvslin_rttovScattK
  
end module tovsLin_mod

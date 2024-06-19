
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
  use varNameList_mod
  use surfaceEmissivity_mod
 
  implicit none
  save
  private

  ! Public procedures
  public :: tvslin_rttov_tl, tvslin_rttov_ad, tvslin_rttov_k
  
  integer, allocatable :: tvslin_bodyIndexFromBtIndex(:,:)
  integer, allocatable :: tvslin_bodyIndexFromBtIndexScatt(:,:)
  type(rttov_chanprof), allocatable :: tvslin_chanProf(:,:)
  type(rttov_chanprof), allocatable :: tvslin_chanProfScatt(:,:)
  
contains

  !--------------------------------------------------------------------------
  !  tvslin_setupPointers
  !--------------------------------------------------------------------------
  subroutine tvslin_setupPointers(runObsOperatorWithHydrometeors, sensorIndex, btCount, &
      btCountScatt, hydroChannelsCount, profileCount, sensorTovsIndexes, &
      lChannelSubset, obsSpaceData)
    !
    ! :Purpose: Allocate and initialize tvslin_bodyIndexFromBtIndex*  tvslin_chanProf*
    !           module variables plus some other local variables.
    !
    implicit none
    
    ! Arguments:
    logical, intent(in)             :: runObsOperatorWithHydrometeors ! flag to control rttovScatt use in linearized RTTOV
    integer, intent(in)             :: sensorIndex                    ! sensor Index in NAMTOV namelist section
    integer, intent(out)            :: btCount                        ! number of BTs simulated using Rttov
    integer, intent(out)            :: btCountScatt                   ! number of BTs simulated using RttovScatt
    integer, intent(out)            :: hydroChannelsCount             ! number of channels simulated using RttovScatt
    integer, intent(out)            :: profileCount                   ! number of profiles for the current sensor
    integer, intent(out)            :: sensorTovsIndexes(:)           ! 
    logical, pointer, intent(inout) :: lChannelSubset(:,:)            ! logical array to setup RttovScatt
    type(struct_obs), intent(inout) :: obsSpaceData                   ! obsSpaceData structure

    ! Locals:
    integer :: tovsIndex, hydroSensorIndex, channelIndex
    
    if (.not. allocated(tvslin_bodyIndexFromBtIndex)) then
      allocate(tvslin_bodyIndexFromBtIndex(tvs_nsensors,tvs_maxNumberOfRadiances))
      allocate(tvslin_bodyIndexFromBtIndexScatt(tvs_nsensors,tvs_maxNumberOfRadiances))
      allocate(tvslin_chanProf(tvs_nsensors,tvs_maxNumberOfRadiances))
      allocate(tvslin_chanProfScatt(tvs_nsensors,tvs_maxNumberOfRadiances))
      tvslin_bodyIndexFromBtIndex(:,:) = -1
      tvslin_bodyIndexFromBtIndexScatt(:,:) = -1
      tvslin_chanProf(:,:) % chan = 0
      tvslin_chanProf(:,:) % prof = 0
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
        call utl_abort('tvslin_setupPointers: you have to initialize channelsUsingHydrometeors(:,:) in NAMTOV namelist section')
      end if
    end if

    profileCount = 0
    do tovsIndex = 1, tvs_nobtov
      ! Currently processed sensor?
      if (tvs_lsensor(tovsIndex) == sensorIndex) then
        profileCount = profileCount + 1
        sensorTovsIndexes(profileCount) = tovsIndex
      end if
    end do

    if (profileCount == 0) return
    
    btCount = tvs_countRadiances(sensorTovsIndexes(1:profileCount), obsSpaceData)
    if (runObsOperatorWithHydrometeors) then
      btCountScatt = tvs_countRadiancesScatt(sensorTovsIndexes(1:profileCount), obsSpaceData, &
          tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount), sensorIndex)
    else
      btCountScatt = 0
    end if
    btCount = btCount - btCountScatt

    if (tvslin_bodyIndexFromBtIndex(sensorIndex,1) == -1) then
      if (btCountScatt > 0) then
        call tvs_getChanprof(sensorTovsIndexes(1:profileCount), obsSpaceData, tvslin_chanProf(sensorIndex,1:btCount), &
            iptobs_cma_opt = tvslin_bodyIndexFromBtIndex(sensorIndex,:), &
            channelList_opt=tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount), &
            excludeChannelsFromList_opt=.true.)
      else
        call tvs_getChanprof(sensorTovsIndexes(1:profileCount), obsSpaceData, tvslin_chanProf(sensorIndex,1:btCount), &
            iptobs_cma_opt = tvslin_bodyIndexFromBtIndex(sensorIndex,:))
      end if
    end if
    if (tvslin_bodyIndexFromBtIndexScatt(sensorIndex,1) == -1 .and. btCountScatt > 0) then
      if (associated(lChannelSubset)) deallocate(lChannelSubset)
      allocate(lChannelSubset(profileCount,tvs_nchan(sensorIndex)))
      call tvs_getChanprof(sensorTovsIndexes(1:profileCount), obsSpaceData, tvslin_chanProfScatt(sensorIndex,1:btCountScatt), &
          lchannel_subset_opt = lChannelSubset, iptobs_cma_opt = tvslin_bodyIndexFromBtIndexScatt(sensorIndex,:), &
          channelList_opt=tvs_channelsUsingHydrometeors(hydroSensorIndex,1:hydroChannelsCount))
    end if
    
  end subroutine tvslin_setupPointers

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
    integer :: allocStatus
    integer :: nobmax
    integer :: sensorIndex, tovsIndex
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
    logical, pointer     :: lChannelSubset(:,:)
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

    allocate(sensorTovsIndexes(tvs_nobtov))
    
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
      
      call tvslin_setupPointers(runObsOperatorWithHydrometeors_tl, sensorIndex, btCount, btCountScatt, &
          hydroChannelsCount, profileCount, sensorTovsIndexes, lChannelSubset, obsSpaceData)

      if (profileCount == 0) cycle sensor_loop
      if (btCount == 0 .and. btCountScatt == 0) cycle  sensor_loop
      
      if (runObsOperatorWithClw_tl) write(*,*) 'tvslin_rttov_tl: using clw_data'
      if (runObsOperatorWithHydrometeors_tl) write(*,*) 'tvslin_rttov_tl: using hydrometeor data'
      sensorType = tvs_coefs(sensorIndex) % coef % id_sensor
      instrum = tvs_coefs(sensorIndex) % coef % id_inst
      
      nobmax = sensorTovsIndexes(profileCount)      
      allocate(sensorHeaderIndexes(profileCount))
      allocate(profilesdata_tl(profileCount))
      allocate(cld_profiles_tl(profileCount))
      allocate(surfTypeIsWater(profileCount))
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
          allocStatus,                  &
          nprofiles=profileCount,       &
          profiles=profilesdata_tl,     &
          nlevels=nlv_T,                &
          opts=tvs_opts(sensorIndex),   &
          asw=1,                        &
          coefs=tvs_coefs(sensorIndex), &
          init=.true.)
      if (allocStatus /= 0) call utl_abort('tvslin_rttov_tl: memory allocation error in rttov_alloc_prof')
      call rttov_alloc_scatt_prof (allocStatus,      &
                                   profileCount,     &
                                   cld_profiles_tl,  &
                                   nlv_T,            &
                                   nhydro=5,         &
                                   nhydro_frac=1,    &
                                   asw=1,            &
                                   init=.true.,      &  
                                   flux_conversion=[1,2,0,0,0])
      if (allocStatus /= 0) call utl_abort('tvslin_rttov_tl: memory allocation error in rttov_alloc_scatt_prof')
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

      deallocate(sensorHeaderIndexes)
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
        if (allocStatus /= 0) call utl_abort('tvslin_rttov_tl: memory allocation error 1 in rttov_alloc_tl')
        allocate(surfem1(btCount))
         !    get Hyperspecral IR emissivities
        if ( tvs_isInstrumHyperSpectral(instrum) ) call tvs_getHIREmissivities(sensorTovsIndexes(1:profileCount), &
            obsSpaceData, surfem1)
       
        call tvs_getOtherEmissivities(tvslin_chanProf(sensorIndex,1:btCount), sensorTovsIndexes, sensorType, instrum, surfem1, calcemis)
        
        if (sensorType == sensor_id_mw) then
          if (col_varExist(columnAnlInc, 'EMMW')) then
            ! Read surface emissivity from column when it's included as an analysis variable
            ! Set the default surface emissivity values
            emissivity_local(:) % emis_in = surfem1(:)

            ! Setup the emissivity_tl from column object
            call sse_setupEmissivityfromState(emissivity_local, obsSpaceData, tvslin_bodyIndexFromBtIndex(sensorIndex,:), tvslin_chanProf(sensorIndex,1:btCount), sensorTovsIndexes, &
                                              tvs_tovsIndex, tvs_headerIndex, tvs_nsensors, tvs_lsensor, tvs_instrumentName, &
                                              tvs_maxChannelNumber, tvs_channelOffset, tvs_ichan, profiles(:) % skin % surftype, &
                                              emissivityProfDt_opt = tvs_emissivityFromTrl)
          else if (tvs_useSfcEmissObsSpace) then
            ! Set the default surface emissivity values
            emissivity_local(:) % emis_in = surfem1(:)
  
            ! Setup the surface emissvity from obsSpaceData Object 
            call sse_emissFromObsSpace(obsSpaceData, emissivity_local, tvslin_bodyIndexFromBtIndex(sensorIndex,:), tvslin_chanProf(sensorIndex,1:btCount), sensorTovsIndexes(1:profileCount), tvs_headerIndex)    
          else
            ! Read surface emissivity from emissivity atlas
            call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorIndex, tvslin_chanProf(sensorIndex,1:btCount), sensorTovsIndexes(1:profileCount))
          end if
        else
          emissivity_local(:) % emis_in = surfem1(:)
        end if
  
        !  2.3  Compute tl radiance with rttov_tl
        
        if (col_varExist(columnAnlInc, 'EMMW') .and. sensorType == sensor_id_mw) then
          call sse_setupEmissivityfromState(emissivity_tl, obsSpaceData, tvslin_bodyIndexFromBtIndex(sensorIndex,:), tvslin_chanProf(sensorIndex,1:btCount), sensorTovsIndexes, tvs_tovsIndex, tvs_headerIndex, &
                                  tvs_nsensors, tvs_lsensor, tvs_instrumentName, tvs_maxChannelNumber, tvs_channelOffset, tvs_ichan, profiles(:) % skin % surftype, columProfTl_opt  = columnAnlInc)
        else
          emissivity_tl(:) % emis_in = 0.0d0
        end if

        errorStatus = errorStatus_success

        !  set nthreads to actual number of threads which will be used.

        nthreads = min(mmpi_numThread, profileCount)  
        call rttov_parallel_tl(                             &
            errorStatus,                                    & ! out
            tvslin_chanProf(sensorIndex,1:btCount),         & ! in
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
          bodyIndex = tvslin_bodyIndexFromBtIndex(sensorIndex,btIndex)
          call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex, &
              radiancedata_tl % bt(btIndex) )
          if ( tvs_debug ) then
            obsOMP = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
            write(*,'(a,i4,2f8.2)') ' ichn,sim,obs= ', &
                tvslin_chanProf(sensorIndex,btIndex) % chan, radiancedata_tl % bt(btIndex), obsOMP
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
        if (allocStatus /= 0) call utl_abort('tvslin_rttov_tl: memory deallocation 1 error in rttov_alloc_tl')
        deallocate(surfem1)
 
      end if
      
      if (btCountScatt >0) then 
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
        if (allocStatus /= 0) call utl_abort('tvslin_rttov_tl: memory allocation error 2 in rttov_alloc_tl')
        
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
            tvslin_chanProfScatt(sensorIndex,1:btCountScatt), & ! channels and profile numbers
            frequencies,                                      & ! array, frequency number for each channel
            lChannelSubset )                                    ! OPTIONAL array of logical flags to indicate a subset of channels
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'tvslin_rttov_tl: fatal error in rttov_scatt_setupindex ', errorStatus
          call utl_abort('tvslin_rttov_tl')
        end if
      
        call tvs_getOtherEmissivities(tvslin_chanProfScatt(sensorIndex,1:btCountScatt), sensorTovsIndexes, sensorType, instrum, surfem1Scatt, calcemisScatt)

        call tvs_getMWemissivityFromAtlas(surfem1Scatt(1:btcountScatt), emissivity_localScatt, sensorIndex, &
            tvslin_chanProfScatt(sensorIndex,1:btCountScatt), sensorTovsIndexes(1:profileCount))
        errorStatus = errorStatus_success
        emissivity_tlScatt(:) % emis_in = 0.0d0
        call rttov_scatt_tl(                                  &
            errorStatus,                                      & ! out
            tvs_opts_scatt(sensorIndex),                      & ! in
            nlv_T,                                            & ! in
            tvslin_chanProfScatt(sensorIndex,1:btCountScatt), & ! in
            frequencies,                                      & ! in
            profiles(sensorTovsIndexes(1:profileCount)),      & ! in  
            cld_profiles(sensorTovsIndexes(1:profileCount)),  & ! in
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
          write(*,*) 'temperature           profile=',profiles(sensorTovsIndexes(1)) % t(:)
          write(*,*) 'temperature increment profile=',profilesdata_tl(1) % t(:)
          call utl_abort('tvslin_rttov_tl')
        end if

        !  2.4  Store hx in obsSpaceData,OBS_WORK
      
        do btIndex = 1, btCountScatt
          bodyIndex = tvslin_bodyIndexFromBtIndexScatt(sensorIndex,btIndex)
          call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex, &
              radiancedata_tlScatt % bt(btIndex) )
          if ( tvs_debug ) then
            obsOMP = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
            write(*,'(a,i4,2f8.2)') ' ichn,sim,obs= ', &
                tvslin_chanprofScatt(sensorIndex,btIndex) % chan, radiancedata_tlScatt % bt(btIndex), obsOMP
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
        if (allocStatus /= 0) call utl_abort('tvslin_rttov_tl: memory deallocation error 2 in rttov_alloc_tl')
      end if
      
      call rttov_alloc_scatt_prof (allocStatus,                 &
                                   profileCount,                &
                                   cld_profiles_tl,             &
                                   nlv_T,                       &
                                   nhydro=5,                    &
                                   nhydro_frac=1,               &
                                   asw=0,                       &   
                                   flux_conversion=[1,2,0,0,0])
      if (allocStatus /= 0) call utl_abort('tvslin_rttov_tl: memory deallocation error in rttov_alloc_scatt_prof')
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
      if (allocStatus /= 0) call utl_abort('tvslin_rttov_tl: memory deallocation error in rttov_alloc_prof')
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
    integer :: allocStatus
    integer :: nthreads
    integer :: nobmax
    integer :: sensorIndex, tovsIndex
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
    logical, pointer :: lChannelSubset(:,:)
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
 
    allocate(sensorTovsIndexes(tvs_nobtov))

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
      
      call tvslin_setupPointers(runObsOperatorWithHydrometeors_ad, sensorIndex, btCount, btCountScatt, &
          hydroChannelsCount, profileCount, sensorTovsIndexes, lChannelSubset, obsSpaceData)
      
      if (profileCount == 0) cycle sensor_loop
      if (btCount == 0 .and. btCountScatt == 0) cycle sensor_loop
      
      sensorType = tvs_coefs(sensorIndex) % coef % id_sensor
      instrum = tvs_coefs(sensorIndex) % coef % id_inst
      nobmax = sensorTovsIndexes(profileCount)
     
      allocate(sensorHeaderIndexes(profileCount))
      allocate(tt_ad(nlv_T,profileCount))
      allocate(hu_ad(nlv_T,profileCount))
      allocate(pressure_ad(nlv_T,profileCount))
      if (.not. tvs_useO3Climatology) then
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

      profileCount = 0       
      ! loop over all obs.
      obs_loop: do tovsIndex = 1, nobmax
        if (tvs_lsensor(tovsIndex)/=sensorIndex) cycle obs_loop
        headerIndex = tvs_headerIndex(tovsIndex)
        profileCount = profileCount + 1
        sensorHeaderIndexes(profileCount) = headerIndex
      end do obs_loop
     
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
      if (allocStatus /= 0) call utl_abort('tvslin_rttov_ad: memory allocation error in rttov_alloc_prof')
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
        if (allocStatus /= 0) call utl_abort('tvslin_rttov_ad: memory allocation error 1 in rttov_alloc_ad')
        allocate(surfem1(btCount))
      
        !  get Hyperspectral IR emissivities
        if ( tvs_isInstrumHyperSpectral(instrum) ) call tvs_getHIREmissivities(sensorTovsIndexes(1:profileCount), obsSpaceData, surfem1)
        
        !     get non Hyperspectral IR emissivities
        call tvs_getOtherEmissivities(tvslin_chanProf(sensorIndex,1:btCount), sensorTovsIndexes, sensorType, instrum, surfem1, calcemis)

        if (sensorType == sensor_id_mw) then
          if (col_varExist(columnAnlInc, 'EMMW')) then
            ! Read surface emissivity from column when it's included as an analysis variable

            ! Set the default surface emissivity values
            emissivity_local(:) % emis_in = surfem1(:)

            ! Setup the surface emissvity from column object to rttov emissivity_local
            call sse_setupEmissivityfromState(emissivity_local, obsSpaceData, tvslin_bodyIndexFromBtIndex(sensorIndex,:), tvslin_chanProf(sensorIndex,1:btCount), sensorTovsIndexes, &
                                        tvs_tovsIndex, tvs_headerIndex, tvs_nsensors, tvs_lsensor, tvs_instrumentName, &
                                        tvs_maxChannelNumber, tvs_channelOffset, tvs_ichan, profiles(:) % skin % surftype, &
                                        emissivityProfDt_opt = tvs_emissivityFromTrl)
          else if (tvs_useSfcEmissObsSpace) then
            ! Set the default surface emissivity values
            emissivity_local(:) % emis_in = surfem1(:)

            ! Setup the surface emissvity from obsSpaceData Object 
            call sse_emissFromObsSpace(obsSpaceData, emissivity_local, tvslin_bodyIndexFromBtIndex(sensorIndex,:), tvslin_chanProf(sensorIndex,1:btCount), sensorTovsIndexes(1:profileCount), tvs_headerIndex)    
          else
            ! Read surface emissivity from emissivity atlas
            call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorIndex, tvslin_chanProf(sensorIndex,1:btCount), sensorTovsIndexes(1:profileCount))
          end if
        else
          emissivity_local(:) % emis_in = surfem1(:)
        end if
        
        do btIndex = 1, btCount
          bodyIndex = tvslin_bodyIndexFromBtIndex(sensorIndex,btIndex)
          radiancedata_ad % bt( btIndex ) = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
        end do
      
        errorStatus = errorStatus_success
        emissivity_ad(:) % emis_in = 0.0d0
        emissivity_ad(:) % emis_out = 0.0d0
        call rttov_parallel_ad(                             &
            errorstatus,                                    &! out
            tvslin_chanProf(sensorIndex,1:btCount),         &! in
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
        if (allocStatus /= 0) call utl_abort('tvslin_rttov_ad: memory deallocation error 1 in rttov_alloc_ad')
        deallocate(surfem1)
      end if

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
        if (allocStatus /= 0) call utl_abort('tvslin_rttov_ad: memory allocation error 2 in rttov_alloc_ad')
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
        if (allocStatus /= 0) call utl_abort('tvslin_rttov_ad: memory allocation error in rttov_alloc_scatt_prof')
        ! Build the list of channels/profiles indices
        call rttov_scatt_setupindex(                          &
            errorStatus,                                      &
            profileCount,                                     &  ! number of profiles
            tvs_nchan(sensorIndex),                           &  ! number of channels 
            tvs_coefs(sensorIndex),                           &  ! coef structure read in from rttov coef file
            tvs_coef_scatt(sensorIndex),                      &  ! coef structure read in from rttov coef file
            btcountScatt,                                     &  ! number of calculated channels
            tvslin_chanProfScatt(sensorIndex,1:btCountScatt), &  ! channels and profile numbers
            frequencies,                                      &  ! array, frequency number for each channel
            lChannelSubset)                                      ! OPTIONAL array of logical flags to indicate a subset of channels
        if (errorStatus /= errorStatus_success) then
          write(*,*) 'tvslin_rttov_ad: fatal error in rttov_scatt_setupindex ', errorStatus
          call utl_abort('tvslin_rttov_ad')
        end if
        !     get non Hyperspectral IR emissivities
        call tvs_getOtherEmissivities(tvslin_chanProfScatt(sensorIndex,1:btCountScatt), sensorTovsIndexes, &
            sensorType, instrum, surfem1Scatt, calcemisScatt)

        call tvs_getMWemissivityFromAtlas(surfem1Scatt(1:btcountScatt), emissivity_localScatt, sensorIndex, &
            tvslin_chanProfScatt(sensorIndex,1:btCountScatt), sensorTovsIndexes(1:profileCount))
        
        do btIndex = 1, btCountScatt
          bodyIndex = tvslin_bodyIndexFromBtIndexScatt(sensorIndex,btIndex)
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
            tvslin_chanProfScatt(sensorIndex,1:btCountScatt), &! in
            frequencies,                                      &! in
            profiles(sensorTovsIndexes(1:profileCount)),      &! in
            cld_profiles(sensorTovsIndexes(1:profileCount)),  &! in
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
          call utl_abort('tvslin_rttov_ad')
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
        if (allocStatus /= 0) call utl_abort('tvslin_rttov_ad: memory deallocation error 2 in rttov_alloc_ad')
        deallocate(surfem1Scatt)
        deallocate(frequencies)
      end if

      !   2.0  Store adjoints in columnData object
      tt_ad(:,:) = 0.d0
      hu_ad(:,:) = 0.d0
      pressure_ad(:,:) = 0.d0
      if (.not. tvs_useO3Climatology) then
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
 
        if (.not. tvs_useO3Climatology) then
          if (tvs_coefs(sensorIndex) % coef % nozone > 0) then
            ! This step is just to transfer the value for ilowlvl_T to the memory space defined by 'col_getColumn(...trim(ozoneVarName))  
            o3_column => col_getColumn(columnAnlInc, headerIndex, trim(ozoneVarName))
            o3_column(ilowlvl_T) = profilesdata_ad(profileIndex) % s2m % o * 1.0d-9
            ozone_ad(:,profileIndex) = profilesdata_ad(profileIndex) % o3(:)
          end if
        end if
      
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
        call sse_setupEmissivityfromState(emissivity_ad, obsSpaceData, tvslin_bodyIndexFromBtIndex(sensorIndex,:), & 
                                          tvslin_chanProf(sensorIndex,1:btCount), sensorTovsIndexes, &
                                          tvs_tovsIndex, tvs_headerIndex, tvs_nsensors, tvs_lsensor, tvs_instrumentName, &
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
      if (allocStatus /= 0) call utl_abort('tvslin_rttov_ad: memory deallocation error 1 in rttov_alloc_ad')

      call rttov_alloc_prof(            &
          allocStatus,                  &
          nprofiles=profileCount,       &
          profiles=profilesdata_ad,     &
          nlevels=nlv_T,                &
          opts=tvs_opts(sensorIndex),   &
          asw=0,                        &
          coefs=tvs_coefs(sensorIndex), &
          init=.true.)
      if (allocStatus /= 0) call utl_abort('tvslin_rttov_ad: memory deallocation error in rttov_alloc_prof')
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
        if (allocStatus /= 0) call utl_abort('tvslin_rttov_ad: memory deallocation error in rttov_alloc_scatt_prof')
        deallocate(cld_profiles_ad)
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
              rf_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'RF')
              do levelIndex = 1, col_getNumLev(columnAnlInc, 'TH')
                rf_column(levelIndex) = rf_column(levelIndex) + rf_ad(levelIndex,profileIndex)
              end do
            end if

            ! snow flux
            if (col_varExist(columnAnlInc,'SF')) then
              sf_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'SF')
              do levelIndex = 1, col_getNumLev(columnAnlInc, 'TH')
                sf_column(levelIndex) = sf_column(levelIndex) + sf_ad(levelIndex,profileIndex)
              end do
            end if

            ! cloud liquid/ice water content
            clw_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'LWCR')
            ciw_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'IWCR')
            do levelIndex = 1, col_getNumLev(columnAnlInc,'TH')
              clw_column(levelIndex) = clw_column(levelIndex) + clwScatt_ad(levelIndex,profileIndex)
              ciw_column(levelIndex) = ciw_column(levelIndex) + ciw_ad(levelIndex,profileIndex)
            end do
          end if ! surfTypeIsWater
        end do ! profileIndex
      end if ! runObsOperatorWithHydrometeors_ad

      deallocate(sensorHeaderIndexes)
      deallocate(tt_ad)
      deallocate(hu_ad)
      deallocate(pressure_ad)
      if (.not. tvs_useO3Climatology) then
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
    type(rttov_emissivity), pointer    :: emissivity_local(:)
    type(rttov_emissivity), pointer    :: emissivity_k(:)
    type (rttov_profile), pointer      :: profiles(:)
    type (rttov_profile), pointer      :: profiles_k(:)
    type(rttov_profile_cloud), pointer :: cld_profiles(:)
    type(rttov_transmission)           :: transmission
    type(rttov_transmission)           :: transmission_k
    type(rttov_radiance)               :: radiancedata_d 
    type(rttov_radiance)               :: radiancedata_k 
    integer, allocatable               :: sensorTovsIndexes(:) 
    integer                            :: allocStatus
    integer                            :: nobmax, profileCount, btCount, btCountScatt
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
    logical, pointer                   :: lChannelSubset(:,:)

    if (tvs_nobtov == 0) return ! exit if there are not tovs data
  
    call tvs_getProfile(profiles, 'nl', cld_profiles)
  
    ! Set index for model's lowest level and model top
  
    nlv_T = col_getNumLev(columnTrlOnAnlIncLev, 'TH')

    allocate(sensorTovsIndexes(tvs_nobtov))
        
    ! Loop over all sensors specified by user
  
    sensor_loop: do sensorIndex = 1, tvs_nsensors
      
      runObsOperatorWithHydrometeors_k = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
                                         col_varExist(columnTrlOnAnlIncLev,'IWCR') .and. &
                                         tvs_isInstrumUsingHydrometeors(tvs_instruments(sensorIndex))
      
      call tvslin_setupPointers(runObsOperatorWithHydrometeors_k, sensorIndex, btCount, btCountScatt, &
          hydroChannelsCount, profileCount, sensorTovsIndexes, lChannelSubset, obsSpaceData)
      
      if (profileCount == 0) cycle sensor_loop
      if (btCount == 0 .and. btCountScatt == 0) cycle sensor_loop
      
      sensorType = tvs_coefs(sensorIndex) % coef % id_sensor
      instrum = tvs_coefs(sensorIndex) % coef % id_inst
      nobmax = sensorTovsIndexes(profileCount)

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
        if (allocStatus /= 0) call utl_abort('tvslin_rttov_k: memory allocation error in rttov_alloc_k')
          
        ! Set nthreads to actual number of threads which will be used.
        nthreads = min(mmpi_numThread, profileCount)  
  
        ! Prepare all input variables required by rttov.
    
        write(*,*) 'tvslin_rttov_k: Get surface emissiviy'
        allocate(surfem1(btCount))
        !    get Hyperspecral IR emissivities
        if (tvs_isInstrumHyperSpectral(instrum)) call tvs_getHIREmissivities(sensorTovsIndexes(1:profileCount), &
                                                          obsSpaceData, surfem1)
  
        call tvs_getOtherEmissivities(tvslin_chanProf(sensorIndex,1:btCount), sensorTovsIndexes, sensorType, instrum, surfem1, calcemis)
  
        if (sensorType == sensor_id_mw) then
          call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorIndex, tvslin_chanProf(sensorIndex,1:btCount), sensorTovsIndexes(1:profileCount))
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
            tvslin_chanProf(sensorIndex,1:btCount),          & ! in
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
        call tvs_writeJacobianAscii(profiles_k, emissivity_k, profiles, tvslin_chanProf(sensorIndex,1:btCount), &
            obsSpaceData, tvs_satelliteName(sensorIndex), tvs_instrumentName(sensorIndex), &
            tvslin_bodyIndexFromBtIndex(sensorIndex,:), sensorTovsIndexes, btCount)
  
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
        if (allocStatus /= 0) call utl_abort('tvslin_rttov_k: memory deallocation error in rttov_alloc_k')
        deallocate(surfem1)
      end if
      
      if (btCountScatt > 0) then
        call utl_abort("tvslin_rttov_k: jacobians not (yet) available when rttov_scatt is used !")
      end if
      
    end do sensor_loop
  
    deallocate (sensorTovsIndexes)
    nullify(profiles)

    write(*,*) 'tvslin_rttov_k: finished'
  end subroutine tvslin_rttov_k

end module tovsLin_mod

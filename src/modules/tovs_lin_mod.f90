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

module tovs_lin_mod
  ! MODULE tovs_lin (prefix='tvslin' category='4. Observation operators')
  !
  ! :Purpose: Derived types, public variables and procedures related to the 
  !           tangent-linear and adjoint versions of RTTOV
  !
  use rttov_interfaces_mod

  use rttov_types, only : rttov_profile, rttov_radiance
  use rttov_const, only : gas_unit_specconc, sensor_id_mw, surftype_sea
  use parkind1, only : jpim, jprb
  use verticalCoord_mod
  use tovs_nl_mod
  use utilities_mod
  use MathPhysConstants_mod
  use obsFilter_mod
  use obsSpaceData_mod
  use columnData_mod
  use tovs_extrap_mod
  use presProfileOperators_mod
 
  implicit none
  save
  private

  public :: tvslin_rttov_tl, tvslin_rttov_ad

  ! public derived types
  ! public derived type through inheritance (from module rttov_types)
  public :: rttov_radiance


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
    type(struct_obs)        :: obsSpaceData  ! obsSpaceData structure
    type(struct_columnData) :: columnAnlInc        ! column structure for pertubation profile
    type(struct_columnData) :: columnTrlOnAnlIncLev ! column structure for background profile

    ! Locals:
    type(struct_vco), pointer :: vco_anl
    integer, allocatable :: sensorTovsIndexes(:) 
    integer, allocatable :: sensorHeaderIndexes(:) 
    integer :: allocStatus(3)
    integer :: nobmax
    integer :: sensorIndex, tovsIndex
    integer :: ilowlvl_M,ilowlvl_T,profileCount,headerIndex,nlv_M,nlv_T
    integer :: profileIndex
    integer :: status, Vcode

    character(len=4) :: ozoneVarName
    logical, allocatable :: surfTypeIsWater(:)
    real(8), pointer :: delTT(:), delHU(:), delP(:)
    real(8), pointer :: delO3(:)
    real(8), pointer :: delCLW(:)
    integer :: btCount
    integer,external :: omp_get_num_threads
    integer :: nthreads, max_nthreads
    integer :: btIndex, bodyIndex
    integer :: instrum
    integer :: sensorType   !sensor type(1=infrared; 2=microwave; 3=high resolution, 4=polarimetric)
    integer :: errorstatus
    integer,allocatable :: sensorBodyIndexes(:)
    real(8), allocatable :: surfem1(:)
    type(rttov_emissivity), pointer :: emissivity_local(:)
    type(rttov_emissivity), pointer :: emissivity_tl(:)
    type(rttov_radiance) :: radiancedata_d   ! radiances full structure buffer used in rttov calls
    type(rttov_radiance) :: radiancedata_tl  ! tl radiances full structure buffer used in rttov calls
    type(rttov_transmission) :: transmission       ! transmission
    type(rttov_transmission) :: transmission_tl    ! transmission tl
    type(rttov_profile), pointer :: profilesdata_tl(:) ! tl profiles buffer used in rttov calls
    type(rttov_chanprof), pointer :: chanprof(:)
    logical, pointer :: calcemis(:)
    logical :: runObsOperatorWithClw_tl
    integer :: asw
    real(8) :: obsOMP
    type (rttov_profile), pointer :: profiles(:)
         
    if (tvs_nobtov == 0) return       ! exit if there are not tovs data

    call tvs_getProfile(profiles, 'tlad')

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

    status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode)
    
  
    !     1.  Get number of threads available and allocate memory for some variables
    !     .   ---------------------------------------------------------------------- 

    !$omp parallel 
    max_nthreads = omp_get_num_threads()
    !$omp end parallel

    allocStatus(:) = 0
    allocate ( sensorTovsIndexes(tvs_nobtov), stat = allocStatus(1) )
    call utl_checkAllocationStatus(allocStatus(1:1), " tvslin_rttov_tl sensorTovsIndexes")
    
    ! 2.  Computation of hx for tovs data only

    
    ! Loop over all sensors specified by user

    sensor_loop:  do sensorIndex = 1, tvs_nsensors

      runObsOperatorWithClw_tl = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
        tvs_opts(sensorIndex) % rt_mw % clw_data .and. &
        tvs_mwInstrumUsingCLW_tl
       
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
      btCount = tvs_countRadiances(sensorTovsIndexes(1:profileCount), obsSpaceData, &
           assim_flag_val_opt=obs_assimilated)
      if ( btCount == 0 ) cycle  sensor_loop
   
      allocate (sensorHeaderIndexes (profileCount), stat= allocStatus(1))
      allocate (profilesdata_tl(profileCount),      stat= allocStatus(2))
      if ( runObsOperatorWithClw_tl ) then
        write(*,*) 'tvslin_rttov_tl: using clw_data'
      end if
      allocate (surfTypeIsWater(profileCount),stat= allocStatus(3))
      call utl_checkAllocationStatus(allocStatus, " tvslin_rttov_tl")
 
      sensorHeaderIndexes(:) = 0 
      
      surfTypeIsWater(:) = .false.

      ! allocate profiledata_tl structures
      asw = 1 ! 1 to allocate
      call rttov_alloc_tl(                     &
              allocStatus(1),                  &
              asw,                             &
              nprofiles=profileCount,          &
              nchanprof=btCount,               &
              nlevels=nlv_T,                   &
              chanprof=chanprof,               &
              opts=tvs_opts(sensorIndex),      &
              profiles_tl=profilesdata_tl,     &
              coefs=tvs_coefs(sensorIndex),    &
              transmission=transmission,       &
              transmission_tl=transmission_tl, &
              radiance=radiancedata_d,         &
              radiance_tl=radiancedata_tl,     &
              calcemis=calcemis,               &
              emissivity=emissivity_local,     &
              emissivity_tl=emissivity_tl,     &
              init=.true.)

      call utl_checkAllocationStatus(allocStatus(1:1), " tovs_rtttov_tl rttov_alloc_tl 1")

      profileCount = 0

      obs_loop: do tovsIndex = 1, nobmax
        if (tvs_lsensor(tovsIndex) /= sensorIndex) cycle obs_loop
        headerIndex = tvs_headerIndex(tovsIndex)
        profileCount = profileCount + 1
        surfTypeIsWater(profileCount) = ( tvs_ChangedStypValue(obsSpaceData,headerIndex) == surftype_sea )
        sensorHeaderIndexes(profileCount) = headerIndex
      end do obs_loop

      do  profileIndex = 1 , profileCount
        profilesdata_tl(profileIndex) % gas_units       = gas_unit_specconc ! all gas profiles should be provided in kg/kg
        profilesdata_tl(profileIndex) % nlevels         =  nlv_T
        profilesdata_tl(profileIndex) % nlayers         =  nlv_T - 1
        if (tvs_coefs(sensorIndex)%coef%nozone > 0) then
          if (tvs_useO3Climatology) then
            profilesdata_tl(profileIndex) % o3(:) =  0.0d0
          else
            delO3 => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),trim(ozoneVarName))
            profilesdata_tl(profileIndex) % o3(1:nlv_T) =  delO3(1:nlv_T) * 1.0d-9 ! Assumes model ozone in ug/kg
            profilesdata_tl(profileIndex) % s2m % o  = col_getElem(columnAnlInc,ilowlvl_T,sensorHeaderIndexes(profileIndex),trim(ozoneVarName)) * 1.0d-9 ! Assumes model ozone in ug/kg
          end if
        end if

        ! using the zero CLW value for land FOV
        if ( runObsOperatorWithClw_tl ) then 
          if ( surfTypeIsWater(profileIndex) ) then
            delCLW => col_getColumn(columnAnlInc,sensorHeaderIndexes(profileIndex),'LWCR')
            profilesdata_tl(profileIndex) % clw(1:nlv_T)  = delCLW(:)
          else
            profilesdata_tl(profileIndex) % clw(1:nlv_T)  = 0.d0
          end if
        end if
         
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
      end do

      deallocate (sensorHeaderIndexes,  stat= allocStatus(1) )
      deallocate (surfTypeIsWater,stat= allocStatus(2)) 
      call utl_checkAllocationStatus(allocStatus, "tvslin_rttov_tl", .false.)

      !  set nthreads to actual number of threads which will be used.

      nthreads = min(max_nthreads, profileCount)  

      !   2.2  Prepare all input variables required by rttov.
     
      allocate ( surfem1(btCount)           ,stat=allocStatus(1))
      allocate ( sensorBodyIndexes(btCount) ,stat=allocStatus(2))
      call utl_checkAllocationStatus(allocStatus(1:2), " tovs_rtttov_tl")
    
      !    get Hyperspecral IR emissivities
      if ( tvs_isInstrumHyperSpectral(instrum) ) call tvs_getHIREmissivities(sensorTovsIndexes(1:profileCount), &
           obsSpaceData, surfem1)

      call tvs_getChanprof(sensorIndex, sensorTovsIndexes(1:profileCount), obsSpaceData, chanprof, &
           iptobs_cma_opt=sensorBodyIndexes)

      call tvs_getOtherEmissivities(chanprof, sensorTovsIndexes, sensorType, instrum, surfem1, calcemis)

      if (sensorType == sensor_id_mw) then
        call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorIndex, chanprof, sensorTovsIndexes(1:profileCount))
      else
        emissivity_local(:)%emis_in = surfem1(:)
      end if
 
      !  2.3  Compute tl radiance with rttov_tl
      
      errorstatus   = 0
      emissivity_tl(:)%emis_in = 0.0d0
      call rttov_parallel_tl(                                     &
           errorstatus,                                           & ! out
           chanprof,                                              & ! in
           tvs_opts(sensorIndex),                                 & ! in
           profiles(sensorTovsIndexes(1:profileCount)),           & ! in
           profilesdata_tl,                                       & ! inout
           tvs_coefs(sensorIndex),                                & ! in
           transmission,                                          & ! inout
           transmission_tl,                                       & ! inout
           radiancedata_d,                                        & ! inout
           radiancedata_tl,                                       & ! inout
           calcemis=calcemis,                                     & ! in
           emissivity=emissivity_local,                           & ! in
           emissivity_tl=emissivity_tl,                           & ! inout
           nthreads=nthreads )                                      ! in
               
      if (errorstatus /= 0) then
        Write(*,*) "Error in rttov_parallel_tl",errorstatus
        write(*,*) 'temperature           profile=',profiles(sensorTovsIndexes(1)) % t(:)
        write(*,*) 'temperature increment profile=',profilesdata_tl(1) % t(:)
        call utl_abort('tovs_rttov_tl')
      end if

      !  2.4  Store hx in obsSpaceData,OBS_WORK
      
      do btIndex = 1, btCount
        
        bodyIndex = sensorBodyIndexes(btIndex)
        call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex, &
             radiancedata_tl % bt(btIndex) )
        if ( tvs_debug ) then
          obsOMP = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
          write(*,'(a,i4,2f8.2)') ' ichn,sim,obs= ', &
               chanprof(btIndex)%chan, radiancedata_tl % bt(btIndex), obsOMP
        end if
        
      end do
 
      ! de-allocate memory
      
      asw = 0 ! 0 to deallocate
      call rttov_alloc_tl(                   &
           allocStatus(1),                   &
           asw,                              &
           nprofiles=profileCount,           &
           nchanprof=btCount,                &
           nlevels=nlv_T,                    &
           chanprof=chanprof,                &
           opts=tvs_opts(sensorIndex),       &
           profiles_tl=profilesdata_tl,      &
           coefs=tvs_coefs(sensorIndex),     &
           transmission=transmission,        &
           transmission_tl=transmission_tl,  &
           radiance=radiancedata_d,          &
           radiance_tl=radiancedata_tl,      &
           calcemis=calcemis,                &
           emissivity=emissivity_local,      &
           emissivity_tl=emissivity_tl )

      deallocate ( surfem1,          stat=allocStatus(2) )
      deallocate ( sensorBodyIndexes,stat=allocStatus(3) )
      call utl_checkAllocationStatus(allocStatus(1:3), " tvslin_rtttov_tl", .false.)
      
    end do sensor_loop

    deallocate ( sensorTovsIndexes )
    nullify( profiles )


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
    type(struct_columnData) :: columnAnlInc
    type(struct_columnData) :: columnTrlOnAnlIncLev
    type(struct_obs)        :: obsSpaceData

    ! locals
    type(struct_vco), pointer :: vco_anl
    integer, allocatable :: sensorTovsIndexes(:) 
    integer, allocatable :: sensorHeaderIndexes(:) 

    integer :: allocStatus(17)
    integer :: omp_get_num_threads, nthreads
    integer :: nobmax
    integer :: sensorIndex, tovsIndex
    integer :: ilowlvl_T,ilowlvl_M,profileCount,headerIndex,nlv_M,nlv_T
    integer :: profileIndex, levelIndex
    integer :: status, Vcode
    real(8), allocatable :: tt_ad(:,:)
    real(8), allocatable :: hu_ad(:,:)
    real(8), allocatable :: pressure_ad(:,:)  
    real(8), allocatable :: ozone_ad(:,:)
    character(len=4) :: ozoneVarName
    real(8), allocatable :: clw_ad(:,:)
    logical, allocatable :: surfTypeIsWater(:)

    real(8), pointer :: uu_column(:),vv_column(:),tt_column(:),hu_column(:),ps_column(:),  &
                        tg_column(:),p_column(:),o3_column(:),clw_column(:)

    integer :: btCount
    integer :: max_nthreads
    integer :: instrum
    integer :: btIndex, bodyIndex
    integer :: sensorType   ! sensor type (1=infrared; 2=microwave; 3=high resolution, 4=polarimetric)
    
    integer, allocatable :: sensorBodyIndexes(:)
    integer :: errorstatus 
    
    real(8), allocatable :: surfem1(:) 
    type(rttov_emissivity), pointer :: emissivity_local(:)
    type(rttov_emissivity), pointer :: emissivity_ad(:)
    type(rttov_transmission) :: transmission,transmission_ad
    type(rttov_radiance) :: radiancedata_ad, radiancedata_d
    
    type(rttov_profile), pointer  :: profilesdata_ad(:) ! ad profiles buffer used in rttov calls
    type(rttov_profile), pointer  :: profiles(:)
    type(rttov_chanprof), pointer :: chanprof(:)
    integer :: asw
    logical, pointer :: calcemis  (:)
    logical :: runObsOperatorWithClw_ad
         
    if (tvs_nobtov == 0) return      ! exit if there are not tovs data

    call tvs_getProfile(profiles, 'tlad')

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
    status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode)



    !     1.  Get number of threads available and allocate memory for some variables
 
    !$omp parallel 
    max_nthreads = omp_get_num_threads()
    !$omp end parallel

    allocate ( sensorTovsIndexes(tvs_nobtov) )

    !     2.  Computation of adjoint hx for tovs data only

    ! Loop over all sensors specified by user

    sensor_loop:do  sensorIndex = 1, tvs_nsensors

      runObsOperatorWithClw_ad = col_varExist(columnTrlOnAnlIncLev,'LWCR') .and. &
        tvs_opts(sensorIndex) % rt_mw % clw_data .and. &
        tvs_mwInstrumUsingCLW_tl
     
      sensorType = tvs_coefs(sensorIndex) % coef% id_sensor
      instrum = tvs_coefs(sensorIndex) % coef% id_inst

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

      if (btCount == 0) cycle sensor_loop
     
      allocStatus(:) = 0
      allocate (sensorHeaderIndexes(profileCount),       stat= allocStatus(1) )
      allocate (tt_ad              (nlv_T,profileCount), stat= allocStatus(2) )
      allocate (hu_ad              (nlv_T,profileCount), stat= allocStatus(3))
      allocate (pressure_ad        (nlv_T,profileCount), stat= allocStatus(4))
      if (.not. tvs_useO3Climatology) then
        if (tvs_coefs(sensorIndex) %coef %nozone > 0) then
          allocate (ozone_ad(nlv_T,profileCount),   stat= allocStatus(5) )
        end if
      end if
      if ( runObsOperatorWithClw_ad ) then
        allocate (clw_ad(nlv_T,profileCount), stat= allocStatus(6))
      end if
      allocate (surfTypeIsWater(profileCount),stat= allocStatus(7))
      surfTypeIsWater(:) = .false.

      call utl_checkAllocationStatus(allocStatus, " tvslin_fill_profiles_ad")

      !  loop over all obs.
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

      !  2.2  Prepare all input variables required by rttov_ad.

      asw = 1 ! 1 for allocation, 0 for deallocation

      call rttov_alloc_ad(                  &
           allocStatus(1),                  &
           asw,                             &
           profileCount,                    &
           btCount,                         &
           nlv_T,                           &
           chanprof,                        &
           opts=tvs_opts(sensorIndex),      &
           profiles_ad=profilesdata_ad,     &
           coefs=tvs_coefs(sensorIndex),    &
           transmission= transmission,      &
           transmission_ad= transmission_ad,&
           radiance=radiancedata_d,         &
           radiance_ad=radiancedata_ad,     &
           calcemis=calcemis,               &
           emissivity=emissivity_local,     &
           emissivity_ad=emissivity_ad,     &
           init=.true.)

      allocate ( surfem1(btCount), stat=allocStatus(2))
      allocate ( sensorBodyIndexes(btCount), stat=allocStatus(3))
      call utl_checkAllocationStatus(allocStatus(1:3), " tvslin_rttov_ad")
      
      !  get Hyperspectral IR emissivities
      
      if ( tvs_isInstrumHyperSpectral(instrum) ) call tvs_getHIREmissivities(sensorTovsIndexes(1:profileCount), obsSpaceData, surfem1)

      ! Build the list of channels/profiles indices

      call tvs_getChanprof(sensorIndex, sensorTovsIndexes(1:profileCount), obsSpaceData, chanprof, sensorBodyIndexes)
      !     get non Hyperspectral IR emissivities
      call tvs_getOtherEmissivities(chanprof, sensorTovsIndexes, sensorType, instrum, surfem1, calcemis)

      if (sensorType == sensor_id_mw) then
        call tvs_getMWemissivityFromAtlas(surfem1(1:btcount), emissivity_local, sensorIndex, chanprof, sensorTovsIndexes(1:profileCount))
      else
        emissivity_local(:)%emis_in = surfem1(:)
      end if
        
      do btIndex = 1, btCount
        bodyIndex = sensorBodyIndexes(btIndex)
        radiancedata_ad % bt( btIndex ) = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
      end do

      !  2.3  Compute ad radiance with rttov_ad

      errorstatus  = 0
      emissivity_ad(:) % emis_in = 0.0d0
      emissivity_ad(:) % emis_out = 0.0d0
  
      call rttov_parallel_ad(                               &
           errorstatus,                                     & ! out
           chanprof,                                        & ! in
           tvs_opts(sensorIndex),                           & ! in
           profiles(sensorTovsIndexes(1:profileCount)),     & ! in
           profilesdata_ad,                                 & ! in
           tvs_coefs(sensorIndex),                          & ! in
           transmission,                                    & ! inout
           transmission_ad,                                 & ! inout
           radiancedata_d,                                  & ! inout
           radiancedata_ad,                                 & ! inout
           calcemis=calcemis,                               & ! in
           emissivity=emissivity_local,                     & ! inout
           emissivity_ad=emissivity_ad,                     & ! inout
           nthreads = nthreads )

      if (errorstatus /= 0) then
        Write(*,*) "Error in rttov_parallel_ad", errorstatus
        call utl_abort('tvslin_rttov_ad')
      end if

      !   2.0  Store adjoints in columnData object
      tt_ad(:,:) = 0.d0
      hu_ad(:,:) = 0.d0
      pressure_ad(:,:) = 0.d0
      if (.not. tvs_useO3Climatology) then
        if (tvs_coefs(sensorIndex) %coef %nozone > 0) ozone_ad(:,:) = 0.d0
      endif
      if ( runObsOperatorWithClw_ad ) clw_ad(:,:) = 0.d0

      do btIndex = 1, btCount
        profileIndex = chanprof(btIndex)%prof
        headerIndex = sensorHeaderIndexes(profileIndex)

        ps_column => col_getColumn(columnAnlInc,headerIndex,'P0')
        p_column  => col_getColumn(columnAnlInc,headerIndex,'P_T')
        tg_column => col_getColumn(columnAnlInc,headerIndex,'TG')
        tt_column => col_getColumn(columnAnlInc,headerIndex,'TT')
        hu_column => col_getColumn(columnAnlInc,headerIndex,'HU')
        uu_column => col_getColumn(columnAnlInc,headerIndex,'UU')
        vv_column => col_getColumn(columnAnlInc,headerIndex,'VV')

        tt_ad(:,profileIndex) =  profilesdata_ad(profileIndex) % t(:)
        hu_ad(:,profileIndex) = profilesdata_ad(profileIndex) % q(:)
        pressure_ad(:,profileIndex)   =  profilesdata_ad(profileIndex) % p(:)
        tg_column(1) = profilesdata_ad(profileIndex) % skin % t 
        tt_column(ilowlvl_T) = profilesdata_ad(profileIndex) % s2m % t
        ps_column(1) = profilesdata_ad(profileIndex) % s2m % p * MPC_MBAR_PER_PA_R8
        hu_column(ilowlvl_T) = 0.d0 
        uu_column(ilowlvl_M) = profilesdata_ad(profileIndex) % s2m % u
        vv_column(ilowlvl_M) = profilesdata_ad(profileIndex) % s2m % v

        if (.not. tvs_useO3Climatology) then
          if (tvs_coefs(sensorIndex) %coef %nozone > 0) then
            ! This step is just to transfer the value for ilowlvl_T to the memory space defined by 'col_getColumn(...trim(ozoneVarName))  
            o3_column => col_getColumn(columnAnlInc,headerIndex,trim(ozoneVarName))
            o3_column(ilowlvl_T) =  profilesdata_ad(profileIndex) % s2m % o * 1.0d-9
            ozone_ad(:,profileIndex) = profilesdata_ad(profileIndex) % o3(:)
          end if
        end if

        if ( runObsOperatorWithClw_ad ) then
          clw_ad(:,profileIndex) = profilesdata_ad(profileIndex) % clw(:)
        end if
      end do

      !     .  2.1  Store adjoints in columnData object
      !     .       -----------------------------------

      do  profileIndex = 1 , profileCount 
        ps_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'P0')
        p_column  => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'P_T')
        tt_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'TT')
        hu_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex), 'HU')
        
        do levelIndex = 1, nlv_T
          p_column(levelIndex) = p_column(levelIndex)  + pressure_ad  (levelIndex,profileIndex) * MPC_MBAR_PER_PA_R8
          tt_column(levelIndex) = tt_column(levelIndex) + tt_ad  (levelIndex,profileIndex)
          hu_column(levelIndex) = hu_column(levelIndex) + hu_ad (levelIndex,profileIndex)
        end do
      end do

      if (.not. tvs_useO3Climatology) then
        if (tvs_coefs(sensorIndex) %coef %nozone > 0) then
          do  profileIndex = 1 , profileCount 
            o3_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex),trim(ozoneVarName))
            do levelIndex = 1, col_getNumLev(columnAnlInc,'TH')
              o3_column(levelIndex) = o3_column(levelIndex) +  ozone_ad(levelIndex,profileIndex) * 1.0d-9
            end do
          end do
        end if
      end if

      if ( runObsOperatorWithClw_ad ) then
        do  profileIndex = 1 , profileCount 
          surfTypeIsWater(profileIndex) = ( tvs_ChangedStypValue(obsSpaceData,sensorHeaderIndexes(profileIndex)) == surftype_sea )
          if ( surfTypeIsWater(profileIndex) ) then
            clw_column => col_getColumn(columnAnlInc, sensorHeaderIndexes(profileIndex),'LWCR')
            do levelIndex = 1, col_getNumLev(columnAnlInc,'TH')
              clw_column(levelIndex) = clw_column(levelIndex) + &
                                       clw_ad(levelIndex,profileIndex)
            end do
          end if
        end do
      end if

      deallocate (sensorHeaderIndexes, stat= allocStatus(1) )
      deallocate (tt_ad,               stat= allocStatus(2) )
      deallocate (hu_ad,               stat= allocStatus(3) )
      deallocate (pressure_ad,         stat= allocStatus(4) )
      if (.not. tvs_useO3Climatology) then
        if (tvs_coefs(sensorIndex) %coef %nozone > 0) then
          deallocate (ozone_ad,        stat= allocStatus(5) )
        end if 
      end if
      if ( runObsOperatorWithClw_ad ) then
        deallocate (clw_ad,stat=allocStatus(6))
      end if
      deallocate (surfTypeIsWater,stat=allocStatus(7))
      
      call utl_checkAllocationStatus(allocStatus, " tvslin_fill_profiles_ad", .false.)
    
      !     de-allocate memory

      asw = 0 ! 0 to deallocate
      call rttov_alloc_ad(                  &
           allocStatus(1),                  &
           asw,                             &
           profileCount,                    &
           btCount,                         &
           nlv_T,                           &
           chanprof,                        &
           opts=tvs_opts(sensorIndex),      &
           profiles_ad=profilesdata_ad,     &
           coefs=tvs_coefs(sensorIndex),    &
           transmission= transmission,      &
           transmission_ad= transmission_ad,&
           radiance=radiancedata_d,         &
           radiance_ad=radiancedata_ad,     &
           calcemis=calcemis,               &
           emissivity=emissivity_local,     &
           emissivity_ad=emissivity_ad )
     
      
      deallocate ( surfem1,           stat=allocStatus(2))
      deallocate ( sensorBodyIndexes, stat=allocStatus(3))
      call utl_checkAllocationStatus(allocStatus(1:3), " tvslin_rttov_ad", .false.)
     
    end do sensor_loop

    ! 3.  Close up

    deallocate ( sensorTovsIndexes )
    nullify( profiles )

  end subroutine tvslin_rttov_ad




 
end module tovs_lin_mod


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
  use rttov_const ,only : gas_unit_specconc
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
  !  lexthum4
  !--------------------------------------------------------------------------
  subroutine lexthum4(pressure, humidityTl, humidity)
    !
    ! :Purpose: tangent linear of extrapolaation of upper level humidity profile
    !           (adapted from exthumtl by J. Eyre).
    !           to extend mixing ratio profile into stratosphere in
    !           a reasonable way.
    !           Take top tropospheric mixing ratio (e.g. near 300 mb) and
    !           extrapolate with given fall off into lower stratosphere
    !           (e.g. to 70 mb).  constrain mixing ratio to be >= zwmin
    !           (e.g. 0.000003 kg/kg).   In upper strat, mixing ratio = zwmin.
    !           Reference: ecmwf tech mem 176.
    !

    implicit none

    !Arguments
    real(8),intent(in)    :: pressure(:)        ! Pressure levels of atm. profiles.
    real(8),intent(inout) :: humidityTl(:,:) ! Tl of humidity profiles.
    real(8),intent(in)    :: humidity(:,:)   ! Humidity profiles.

    !Locals:
    integer :: nlevels
    integer :: nprofiles

    real(8), allocatable :: zpres3(:)

    real(8) zwb5,zwb
    real(8), parameter :: pressureLimit = 70.0d0  ! press limits (in hpa) of region  to be extrapolated
    integer :: topIndex, profileIndex, levelIndex

    nlevels = size( pressure )
    nprofiles = size( humidity, dim =2)


    !    extrapolate humidity profile.

    !    find top level of given profile
    topIndex = -1
    do levelIndex = nlevels, 1, -1
      if (pressure(levelIndex) < FILT_RLIMLVHU) then
        topIndex = levelIndex
        exit
      end if
    end do

    ! Null extrapolation case
    if (topIndex == -1) return
    allocate ( zpres3( nlevels ) )

    ! Constants defining p**3 fall off around tropopause
    do levelIndex=1, topINdex
      zpres3(levelIndex) = ( pressure(levelIndex) / pressure(topIndex + 1) ) ** 3
    end do

    do profileIndex = 1, nprofiles
      zwb = humidityTl (topIndex + 1, profileIndex)
      zwb5 = humidity(topIndex + 1, profileIndex)
      do levelIndex = 1, topIndex
        if (pressure(levelIndex) < pressureLimit) then
          humidityTl(levelIndex,profileIndex) = 0.d0
        else
          if ( zwb5 * zpres3(levelIndex) <= MPC_MINIMUM_HU_R8 ) then
            humidityTl(levelIndex,profileIndex) = 0.d0
          else
            humidityTl(levelIndex,profileIndex) = zwb * zpres3(levelIndex)
          end if
        end if
      end do
    end do

    deallocate ( zpres3 )
    
  end subroutine lexthum4


  !--------------------------------------------------------------------------
  !  tvslin_rttov_tl
  !--------------------------------------------------------------------------
  subroutine tvslin_rttov_tl(column, columng, obsSpaceData)
    !
    ! :Purpose: Tangent linear of computation of radiance with rttov_tl
    !   
    implicit none

    ! Arguments:
    type(struct_obs)        :: obsSpaceData ! obsSpaceData structure
    type(struct_columnData) :: column        ! column structure for pertubation profile
    type(struct_columnData) :: columng       ! column structure for background profile

    ! Locals:
    type(struct_vco), pointer :: vco_anl
    integer, allocatable :: sensorTovsIndexes(:) 
    integer, allocatable :: sensorHeaderIndexes(:) 
    integer :: allocStatus(17)
    logical :: diagTtop,TopAt10hPa
    integer :: nobmax
    integer :: sensorIndex, tovsIndex
    integer :: ilowlvl_M,ilowlvl_T,profileCount,headerIndex,nlv_M,nlv_T
    integer :: levelIndex, profileIndex
    integer :: status, Vcode
    integer :: modelTopIndex, levelsBelowModelTop

    real(8), allocatable :: interpolatedTemperature(:,:)
    real(8), allocatable :: interpolatedHumidityTl(:,:)
    real(8), allocatable :: loginterpolatedHumidityTl(:,:)
    real(8), allocatable :: extrapolatedTemperatureTl(:,:)
    real(8), allocatable :: extrapolatedHumidityTl(:,:)
    real(8), allocatable :: modelPressureLevels(:,:)
    real(8), allocatable :: dPdPs(:,:)
    real(8), allocatable :: modelTemperatureTl(:,:)
    real(8), allocatable :: modelHumidityTl(:,:)
    real(8), allocatable :: logModelHumidityTl(:,:)
    real(8), allocatable :: modelTemperature(:,:)
    real(8), allocatable :: modelHumidity(:,:)
    real(8), allocatable :: logModelHumidity(:,:)
    real(8), allocatable :: extrapolatedHumidity(:,:)
    real(8), allocatable :: modelPressureTl(:,:)
    real(8), allocatable :: rttovPressureLevels(:)
    real(8) :: modelTopPressure
    real(8), pointer :: delTT(:), delHU(:), TTb(:), HUb(:), Pres(:), delP(:)
    integer :: nRttovLevels
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
    type(rttov_emissivity), pointer :: emissivityTl(:)
    type(rttov_radiance) :: radiancedata_d   ! radiances full structure buffer used in rttov calls
    type(rttov_radiance) :: radiancedataTl  ! tl radiances full structure buffer used in rttov calls
    type(rttov_transmission) :: transmission       ! transmission
    type(rttov_transmission) :: transmissionTl    ! transmission tl
    type(rttov_profile), pointer :: profilesdataTl(:) ! tl profiles buffer used in rttov calls
    type(rttov_chanprof), pointer :: chanprof(:)
    logical, pointer :: calcemis(:)
    integer ::  asw
         
    if (tvs_nobtov == 0) return       ! exit if there are not tovs data

    !  1.  Set index for model's lowest level and model top

    nlv_M = col_getNumLev(columng,'MM')
    nlv_T = col_getNumLev(columng,'TH')

    if ( col_getPressure(columng,1,1,'TH') < col_getPressure(columng,nlv_T,1,'TH') ) then
      ilowlvl_M = nlv_M
      ilowlvl_T = nlv_T
    else
      ilowlvl_M = 1
      ilowlvl_T = 1
    end if

    vco_anl => col_getVco(columng)
    status = vgd_get(vco_anl%vgrid, key='ig_1 - vertical coord code', value = Vcode)
    diagTtop = (Vcode==5002)
    
    ! find model level top, within 0.000001 mbs.
    modelTopPressure = (col_getPressure(columng, 1, 1, 'TH')  * MPC_MBAR_PER_PA_R8) - 0.000001d0

    TopAt10hPa = ( abs( modelTopPressure - 10.0d0 ) <= .1d0 )

    !  1.  Get number of threads available and allocate memory for some variables

    !$omp parallel 
    max_nthreads = omp_get_num_threads()
    !$omp end parallel

    allocStatus(:) = 0
    allocate ( sensorTovsIndexes(tvs_nobtov), stat = allocStatus(1) )
    call utl_checkAllocationStatus(allocStatus(1:1), " tvslin_rttov_tl sensorTovsIndexes")
    
    ! 2.  Computation of hx for tovs data only

    
    ! Loop over all sensors specified by user

    sensor_loop:  do sensorIndex = 1, tvs_nsensors
       
      nRttovLevels = tvs_coefs(sensorIndex) % coef % nlevels
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
   
      allocate (rttovPressureLevels(nRttovLevels))
      rttovPressureLevels = tvs_coefs(sensorIndex)% coef % ref_prfl_p
      modelTopIndex = 1
      do levelIndex = 2, nRttovLevels
        if ( modelTopPressure >= rttovPressureLevels(levelIndex - 1).and. &
             modelTopPressure < rttovPressureLevels(levelIndex)       ) then
          modelTopIndex = levelIndex
          exit
        end if
      end do
      levelsBelowModelTop = (nRttovLevels - modelTopIndex + 1)

      allocate (sensorHeaderIndexes (profileCount),      stat= allocStatus(1) )
      allocate (interpolatedTemperature     (levelsBelowModelTop,profileCount),  stat= allocStatus(2) )
      allocate (interpolatedHumidityTl    (levelsBelowModelTop,profileCount),  stat= allocStatus(3) )
      allocate (loginterpolatedHumidityTl (levelsBelowModelTop,profileCount),  stat= allocStatus(4) )
      allocate (extrapolatedTemperatureTl  (nRttovLevels  ,profileCount),stat= allocStatus(5) )
      allocate (extrapolatedHumidityTl  (nRttovLevels  ,profileCount),stat= allocStatus(6) )
      allocate (modelPressureLevels     (nlv_T,profileCount),    stat= allocStatus(7) )
      allocate (dPdPs     (nlv_T,profileCount),    stat= allocStatus(8) )
      allocate (modelTemperatureTl     (nlv_T,profileCount),    stat= allocStatus(9) )
      allocate (modelHumidityTl    (nlv_T,profileCount),    stat= allocStatus(10))
      allocate (logModelHumidityTl (nlv_T,profileCount),    stat= allocStatus(11))
      allocate (modelTemperature        (nlv_T,profileCount),    stat= allocStatus(12))
      allocate (modelHumidity       (nlv_T,profileCount),    stat= allocStatus(13))
      allocate (logModelHumidity    (nlv_T,profileCount),    stat= allocStatus(14))
      allocate (extrapolatedHumidity     (nRttovLevels,profileCount),  stat= allocStatus(15))
      allocate (modelPressureTl     (nlv_T,profileCount),    stat= allocStatus(16))
      allocate (profilesdataTl(profileCount),     stat= allocStatus(17))

      call utl_checkAllocationStatus(allocStatus, " tvslin_rttov_tl")
 
      sensorHeaderIndexes(:) = 0 
      extrapolatedTemperatureTl(:,:) = 0.0d0
      modelPressureLevels(:,:) = 0.0d0
      dPdPs(:,:) = 0.0d0
      modelTemperatureTl(:,:) = 0.0d0
      modelHumidityTl(:,:) = 0.0d0
      modelTemperature(:,:) = 0.0d0
      modelHumidity(:,:) = 0.0d0
      extrapolatedHumidity(:,:) = 0.0d0
      modelPressureTl(:,:) = 0.0d0
      interpolatedTemperature(:,:) = 0.0d0
      interpolatedHumidityTl(:,:) = 0.0d0

      ! allocate profiledataTl structures
      asw = 1 ! 1 to allocate
      call rttov_alloc_tl(                  &
              allocStatus(1),               &
              asw,                          &
              nprofiles=profileCount,       &
              nchanprof=btCount,            &
              nlevels=nRttovLevels,              &
              chanprof=chanprof,            &
              opts=tvs_opts(sensorIndex),   &
              profiles_tl=profilesdataTl,  &
              coefs=tvs_coefs(sensorIndex), &
              transmission=transmission,    &
              transmission_tl=transmissionTl, &
              radiance=radiancedata_d,      &
              radiance_tl=radiancedataTl,  &
              calcemis=calcemis,            &
              emissivity=emissivity_local,  &
              emissivity_tl=emissivityTl,  &
              init=.true.)

      call utl_checkAllocationStatus(allocStatus(1:1), " tovs_rtttov_tl rttov_alloc_tl 1")

      profileCount = 0

      obs_loop: do tovsIndex = 1, nobmax
        if (tvs_lsensor(tovsIndex) /= sensorIndex) cycle obs_loop

        headerIndex = tvs_headerIndex(tovsIndex)
        profileCount = profileCount + 1

        delP => col_getColumn(column,headerIndex,'P_T')
        delTT => col_getColumn(column,headerIndex,'TT')
        delHU => col_getColumn(column,headerIndex,'HU')
        TTb => col_getColumn(columng,headerIndex,'TT')
        HUb => col_getColumn(columng,headerIndex,'HU')
        Pres => col_getColumn(columng,headerIndex,'P_T')
        do levelIndex = 1, nlv_T
          modelPressureTl(levelIndex,profileCount) = delP(levelIndex) * MPC_MBAR_PER_PA_R8
          modelTemperatureTl(levelIndex,profileCount) = delTT(levelIndex)
          modelHumidityTl(levelIndex,profileCount) = delHU(levelIndex)
          modelTemperature(levelIndex,profileCount)  = TTb(levelIndex)
          modelHumidity(levelIndex,profileCount)  = HUb(levelIndex)
          modelPressureLevels(levelIndex,profileCount)  = Pres(levelIndex) *MPC_MBAR_PER_PA_R8
          dPdPs(levelIndex,profileCount)  = 1.0d0
        end do
        
        ! Fix pour eviter probleme au toit avec GEM 4
        ! (grosse varibilite temperature au dernier niveau thermo due 
        !  a l'extrapolation utilisee)
        if ( diagTtop ) then
          modelTemperatureTl(1,profileCount) =  0.d0
          modelHumidityTl(1,profileCount) =  0.d0
          modelTemperature(1,profileCount) =  modelTemperature(2,profileCount) + tvs_mesosphereLapseRate *  &
               log( col_getPressure(columng,1,headerIndex,'TH') /  &
               col_getPressure(columng,2,headerIndex,'TH') )
          modelHumidity(1,profileCount) =  modelHumidity(2,profileCount)
        end if
        
        sensorHeaderIndexes(profileCount) = headerIndex
        
      end do obs_loop

      !  2.1  Vertical interpolation of model temperature and logarithm of
      !       specific humidity to pressure levels required by tovs rt model

      do profileIndex = 1, profileCount
        extrapolatedHumidity(1:nRttovLevels,profileIndex) =  tvs_profiles(sensorTovsIndexes(profileIndex)) % q(1:nRttovLevels)
      end do
      
      interpolatedTemperature (:,:) = 0.0d0
      interpolatedHumidityTl(:,:) = 0.0d0
      
      !$omp parallel do private(profileIndex)
      do profileIndex=1, profileCount

        call ppo_IntAvgTl_v2(modelPressureLevels(:,profileIndex:profileIndex),dPdPs(:,profileIndex:profileIndex),modelTemperatureTl(:,profileIndex:profileIndex), &
             modelTemperature(:,profileIndex:profileIndex),modelPressureTl(:,profileIndex:profileIndex),nlv_T,1, &
             levelsBelowModelTop,rttovPressureLevels(modelTopIndex:nRttovLevels),interpolatedTemperature(:,profileIndex:profileIndex))

        logModelHumidity(:,profileIndex) = log( modelHumidity(:,profileIndex) )
        logModelHumidityTl(:,profileIndex) = modelHumidityTl(:,profileIndex) / modelHumidity(:,profileIndex)
        call ppo_IntAvgTl_v2(modelPressureLevels(:,profileIndex:profileIndex),dPdPs(:,profileIndex:profileIndex),logModelHumidityTl(:,profileIndex:profileIndex), &
             logModelHumidity(:,profileIndex:profileIndex),modelPressureTl(:,profileIndex:profileIndex),nlv_T,1, &
             levelsBelowModelTop,rttovPressureLevels(modelTopIndex:nRttovLevels),loginterpolatedHumidityTl(:,profileIndex:profileIndex))

        interpolatedHumidityTl(:,profileIndex) = loginterpolatedHumidityTl(:,profileIndex) * extrapolatedHumidity(modelTopIndex:nRttovLevels,profileIndex)


      end do
      !$omp end parallel do
      
      !  2.2  Extrapolation of temperature profile above 10mb
      extrapolatedTemperatureTl(:,:) = 0.0d0
      if ( .not. TopAt10hPa ) then
        do profileIndex = 1, profileCount
          extrapolatedTemperatureTl(modelTopIndex:nRttovLevels,profileIndex) = interpolatedTemperature(1:levelsBelowModelTop,profileIndex)
          extrapolatedTemperatureTl(1:modelTopIndex-1,profileIndex) = 0.d0
        end do
      else
        call lextrap (interpolatedTemperature,extrapolatedTemperatureTl,levelsBelowModelTop,nRttovLevels,profileCount)
      end if
      
      !   2.3  Extrapolation of humidity profile (kg/kg)
      !        above rlimlvhu (normally 300mbs or 70mbs)
      
      do profileIndex = 1, profileCount
        do levelIndex = 1, modelTopIndex - 1
          extrapolatedHumidityTl(levelIndex,profileIndex) = 0.d0
        end do
        do levelIndex = 1, levelsBelowModelTop
          extrapolatedHumidityTl(nRttovLevels - levelsBelowModelTop + levelIndex,profileIndex) = interpolatedHumidityTl(levelIndex,profileIndex)
        end do
      end do
      
      if ( TopAt10hPa ) then
        if ( tvs_debug ) then
          do profileIndex = 1, profileCount
            write(*,*)'extrapolatedHumidity_tl*1000 avant exthum4    = '
            write(*,'(1x,10f8.4)')(extrapolatedHumidityTl(levelIndex,profileIndex) * 1000.d0,levelIndex=1,nRttovLevels)
            write(*,*)' '
          end do
        end if
        call lexthum4 (rttovPressureLevels(1:nRttovLevels),extrapolatedHumidityTl,extrapolatedHumidity)
        if ( tvs_debug ) then
          do profileIndex = 1, profileCount
            write(*,*)'extrapolatedHumidityTl*1000 apres exthum4    = '
            write(*,'(1x,10f8.4)')(extrapolatedHumidityTl(levelIndex,profileIndex) * 1000.d0,levelIndex=1,nRttovLevels)
            write(*,*)' '
          end do
        end if
      end if

      do  profileIndex = 1 , profileCount
        profilesdataTl(profileIndex) % gas_units       = gas_unit_specconc ! all gas profiles should be provided in kg/kg
        profilesdataTl(profileIndex) % nlevels         =  nRttovLevels
        profilesdataTl(profileIndex) % nlayers         =  nRttovLevels - 1
        if (tvs_coefs(sensorIndex)%coef%nozone > 0) profilesdataTl(profileIndex) % o3(:) =  0.d0
        
        profilesdataTl(profileIndex) % ctp             = 0.0d0
        profilesdataTl(profileIndex) % cfraction       = 0.0d0
        profilesdataTl(profileIndex) % zenangle        = 0.0d0
        profilesdataTl(profileIndex) % azangle         = 0.0d0
        profilesdataTl(profileIndex) % skin % surftype = 0
        profilesdataTl(profileIndex) % skin % t        = col_getElem(column,1,sensorHeaderIndexes(profileIndex),'TG')
        profilesdataTl(profileIndex) % skin % fastem(:)= 0.0d0
        profilesdataTl(profileIndex) % skin % salinity = 0.0d0
        profilesdataTl(profileIndex) % s2m % t         = col_getElem(column,ilowlvl_T,sensorHeaderIndexes(profileIndex),'TT')

        
        profilesdataTl(profileIndex) % s2m % q         = 0.d0

        profilesdataTl(profileIndex) % s2m % p         = col_getElem(column,1,sensorHeaderIndexes(profileIndex),'P0')*MPC_MBAR_PER_PA_R8
        profilesdataTl(profileIndex) % s2m % u         = col_getElem(column,ilowlvl_M,sensorHeaderIndexes(profileIndex),'UU')
        profilesdataTl(profileIndex) % s2m % v         = col_getElem(column,ilowlvl_M,sensorHeaderIndexes(profileIndex),'VV')
        
        profilesdataTl(profileIndex) % p(1:nRttovLevels)    = 0.d0
        profilesdataTl(profileIndex) % t(1:nRttovLevels)    = extrapolatedTemperatureTl(1:nRttovLevels,profileIndex)
        profilesdataTl(profileIndex) % q(1:nRttovLevels)    = extrapolatedHumidityTl(1:nRttovLevels,profileIndex)
      end do

      deallocate (sensorHeaderIndexes, stat= allocStatus(1) )
      deallocate (interpolatedTemperature,         stat= allocStatus(2) )
      deallocate (interpolatedHumidityTl,        stat= allocStatus(3) )
      deallocate (loginterpolatedHumidityTl,     stat= allocStatus(4) )
      deallocate (extrapolatedTemperatureTl,      stat= allocStatus(5) )
      deallocate (extrapolatedHumidityTl,      stat= allocStatus(6) )
      deallocate (modelPressureLevels,         stat= allocStatus(7) )
      deallocate (dPdPs,         stat= allocStatus(8) )
      deallocate (modelTemperatureTl,         stat= allocStatus(9) )
      deallocate (modelHumidityTl,        stat= allocStatus(10))
      deallocate (logModelHumidityTl,     stat= allocStatus(11))
      deallocate (modelTemperature,            stat= allocStatus(12))
      deallocate (modelHumidity,           stat= allocStatus(13))
      deallocate (logModelHumidity,        stat= allocStatus(14))
      deallocate (extrapolatedHumidity,         stat= allocStatus(15))
      deallocate (modelPressureTl,         stat= allocStatus(16))
      deallocate (rttovPressureLevels,         stat= allocStatus(17))
      call utl_checkAllocationStatus(allocStatus, "tvslin_rttov_tl", .false.)

      !     set nthreads to actual number of threads which will be used.

      nthreads = min(max_nthreads, profileCount)  

      !   2.2  Prepare all input variables required by rttov.
     
      allocate ( surfem1      (btCount)       ,stat=allocStatus(1))
      allocate ( sensorBodyIndexes   (btCount)      ,stat=allocStatus(2))
      call utl_checkAllocationStatus(allocStatus(1:2), " tovs_rtttov_tl")
    
      !    get Hyperspecral IR emissivities
      if ( tvs_isInstrumHyperSpectral(instrum) ) call tvs_getHIREmissivities(sensorTovsIndexes(1:profileCount), &
           obsSpaceData, surfem1)

      call tvs_getChanprof(sensorIndex, sensorTovsIndexes(1:profileCount), obsSpaceData, chanprof, &
           iptobs_cma_opt=sensorBodyIndexes)

      call tvs_getOtherEmissivities(chanprof, sensorTovsIndexes, sensorType, instrum, surfem1, calcemis)
      
      emissivity_local(:)%emis_in = surfem1(:)

 
      !  2.3  Compute tl radiance with rttov_tl
      
      errorstatus   = 0
      emissivityTl(:)%emis_in = 0.0d0
      call tmg_start(87,'rttov_tl')
      call rttov_parallel_tl(         &
           errorstatus,               & ! out
           chanprof,                  & ! in
           tvs_opts(sensorIndex),     & ! in
           tvs_profiles(sensorTovsIndexes(1:profileCount)),  & ! in
           profilesdataTl,           & ! inout
           tvs_coefs(sensorIndex),    & ! in
           transmission,              & ! inout
           transmissionTl,           & ! inout
           radiancedata_d,            & ! inout
           radiancedataTl,           & ! inout
           calcemis,                  & ! in
           emissivity_local,          & ! in
           emissivityTl,             & ! inout
           nthreads=nthreads )          ! in
               
      if (errorstatus /= 0) then
        Write(*,*) "Error in rttov_parallel_tl",errorstatus
        write(*,*) 'temperature           profile=',tvs_profiles(sensorTovsIndexes(1)) % t(:)
        write(*,*) 'temperature increment profile=',profilesdataTl(1) % t(:)
        call utl_abort('tovs_rttov_tl')
      end if

      call tmg_stop(87)

      !  2.4  Store hx in obsSpaceData,OBS_WORK
      
      do btIndex = 1, btCount
        
        bodyIndex = sensorBodyIndexes(btIndex)
        call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex, &
             radiancedataTl % bt(btIndex) )
        if ( tvs_debug ) then
          write(*,'(a,i4,2f8.2)') ' ichn,sim,obs= ', &
               chanprof(btIndex)%chan,   radiancedataTl % bt(btIndex), &
               obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
        end if
        
      end do

 
      ! de-allocate memory
      
      asw = 0 ! 0 to deallocate
      call rttov_alloc_tl(                   &
           allocStatus(1),                   &
           asw,                              &
           nprofiles=profileCount,           &
           nchanprof=btCount,                &
           nlevels=nRttovLevels,                  &
           chanprof=chanprof,                &
           opts=tvs_opts(sensorIndex),       &
           profiles_tl=profilesdataTl,      &
           coefs=tvs_coefs(sensorIndex),     &
           transmission=transmission,        &
           transmission_tl=transmissionTl,  &
           radiance=radiancedata_d,          &
           radiance_tl=radiancedataTl,      &
           calcemis=calcemis,                &
           emissivity=emissivity_local,      &
           emissivity_tl=emissivityTl )

      deallocate ( surfem1,    stat=allocStatus(2) )
      deallocate ( sensorBodyIndexes,stat=allocStatus(3) )
      call utl_checkAllocationStatus(allocStatus(1:3), " tvslin_rtttov_tl", .false.)
      
    end do sensor_loop

    deallocate ( sensorTovsIndexes )


  end subroutine tvslin_rttov_tl

  !--------------------------------------------------------------------------
  !  tvslin_rttov_ad
  !--------------------------------------------------------------------------
  subroutine tvslin_rttov_ad( column, columng, obsSpaceData )
    !
    ! :Purpose: Adjoint of computation of radiance with rttov_ad
    !

    implicit none

    ! Arguments:
    type(struct_columnData) :: column
    type(struct_columnData) :: columng
    type(struct_obs)        :: obsSpaceData

    ! locals
    type(struct_vco), pointer :: vco_anl
    integer, allocatable :: sensorTovsIndexes(:) 
    integer, allocatable :: sensorHeaderIndexes(:) 
    
    integer :: allocStatus(17)
    logical :: diagTtop,TopAt10hPa
    integer :: omp_get_num_threads, nthreads
    integer :: nRttovLevels,nobmax
    integer :: sensorIndex, tovsIndex
    integer :: ilowlvl_T,ilowlvl_M,profileCount,headerIndex,nlv_M,nlv_T
    integer :: profileIndex, levelIndex
    integer :: status, Vcode
    integer :: modelTopIndex, levelsBelowModelTop
    
    real(8), allocatable :: interpolatedTemperatureAd    (:,:)
    real(8), allocatable :: interpolatedHumidityAd   (:,:)
    real(8), allocatable :: loginterpolatedHumidityAd(:,:)
    real(8), allocatable :: extrapolatedTemperatureAd (:,:)
    real(8), allocatable :: extrapolatedHumidityAd (:,:)
    real(8), allocatable :: modelPressureLevels    (:,:)
    real(8), allocatable :: dPdPs    (:,:)
    real(8), allocatable :: modelTemperatureAd    (:,:)
    real(8), allocatable :: modelHumidityAd   (:,:)
    real(8), allocatable :: logModelHumidityAd(:,:)
    real(8), allocatable :: modelTemperature       (:,:)
    real(8), allocatable :: modelHumidity      (:,:)
    real(8), allocatable :: logModelHumidity   (:,:)
    real(8), allocatable :: extrapolatedHumidity    (:,:)
    real(8), allocatable :: modelPressureAd    (:,:)
    real(8), allocatable :: rttovPressureLevels    (:)

    real(8) :: modelTopPressure
   
    real(8), pointer :: uu_column(:),vv_column(:),tt_column(:),hu_column(:),ps_column(:),tg_column(:),p_column(:)
    real(8), pointer :: TTb(:), HUb(:), Pres(:)

    integer :: btCount
    integer :: max_nthreads
    integer :: instrum
    integer :: btIndex, bodyIndex
    integer :: sensorType   ! sensor type (1=infrared; 2=microwave; 3=high resolution, 4=polarimetric)
    
    integer, allocatable :: sensorBodyIndexes(:)
    integer :: errorstatus 
    
    real(8), allocatable :: surfem1(:) 
    type(rttov_emissivity), pointer :: emissivity_local (:)
    type(rttov_emissivity), pointer :: emissivityAd (:)
    type(rttov_transmission) :: transmission,transmissionAd
    type(rttov_radiance) :: radiancedataAd, radiancedata_d
    
    type(rttov_profile), pointer  :: profilesdataAd (:) ! ad profiles buffer used in rttov calls
    type(rttov_chanprof), pointer :: chanprof(:)
    integer :: asw
    logical, pointer :: calcemis  (:)
         
    if (tvs_nobtov == 0) return      ! exit if there are not tovs data


    !     1.    Set index for model's lowest level and model top

    nlv_M = col_getNumLev(columng,'MM')
    nlv_T = col_getNumLev(columng,'TH')

    if (  col_getPressure(columng,1,1,'TH') < col_getPressure(columng,nlv_T,1,'TH') ) then
      ilowlvl_M = nlv_M
      ilowlvl_T = nlv_T
    else
      ilowlvl_M = 1
      ilowlvl_T = 1
    end if

    vco_anl => col_getVco(columng)
    status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode)
    diagTtop = (Vcode == 5002)

    ! find model level top, within 0.000001 mbs.
    modelTopPressure = ( col_getPressure(columng, 1, 1, 'TH')  * MPC_MBAR_PER_PA_R8) - 0.000001d0

    TopAt10hPa = ( abs( modelTopPressure - 10.0d0 ) <= .1d0 )


    !     1.  Get number of threads available and allocate memory for some variables
 
    !$omp parallel 
    max_nthreads = omp_get_num_threads()
    !$omp end parallel

    allocate ( sensorTovsIndexes(tvs_nobtov) )

    !     2.  Computation of adjoint hx for tovs data only

    ! Loop over all sensors specified by user

    sensor_loop:do  sensorIndex = 1, tvs_nsensors
     
      nRttovLevels = tvs_coefs(sensorIndex) %coef % nlevels
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

     
      allocate (rttovPressureLevels(nRttovLevels))
      rttovPressureLevels = tvs_coefs(sensorIndex)% coef % ref_prfl_p
      modelTopIndex = 1
      do levelIndex = 2, nRttovLevels
        if ( modelTopPressure >= rttovPressureLevels(levelIndex - 1) .and.    &
             modelTopPressure < rttovPressureLevels(levelIndex)        ) then
          modelTopIndex = levelIndex
          exit
        end if
      end do

      levelsBelowModelTop = (nRttovLevels - modelTopIndex + 1)
     
      allocStatus(:) = 0
      allocate (sensorHeaderIndexes(profileCount)    ,stat= allocStatus(1) )
      allocate (interpolatedTemperatureAd    (levelsBelowModelTop,profileCount),stat= allocStatus(2) )
      allocate (interpolatedHumidityAd   (levelsBelowModelTop,profileCount),stat= allocStatus(3) )
      allocate (loginterpolatedHumidityAd(levelsBelowModelTop,profileCount),stat= allocStatus(4) )
      allocate (extrapolatedTemperatureAd (nRttovLevels,profileCount),stat= allocStatus(5) )
      allocate (extrapolatedHumidityAd (nRttovLevels,profileCount),stat= allocStatus(6) )
      allocate (modelPressureLevels    (nlv_T,profileCount)  ,stat= allocStatus(7) )
      allocate (dPdPs    (nlv_T,profileCount)  ,stat= allocStatus(8) )
      allocate (modelTemperatureAd    (nlv_T,profileCount)  ,stat= allocStatus(9) )
      allocate (modelHumidityAd   (nlv_T,profileCount)  ,stat= allocStatus(10))
      allocate (logModelHumidityAd(nlv_T,profileCount)  ,stat= allocStatus(11))
      allocate (modelTemperature       (nlv_T,profileCount)  ,stat= allocStatus(12))
      allocate (modelHumidity      (nlv_T,profileCount)  ,stat= allocStatus(13))
      allocate (logModelHumidity   (nlv_T,profileCount)  ,stat= allocStatus(14))
      allocate (extrapolatedHumidity    (nRttovLevels,profileCount),stat= allocStatus(15))
      allocate (modelPressureAd    (nlv_T,profileCount)  ,stat= allocStatus(16))

      call utl_checkAllocationStatus(allocStatus, " tvslin_fill_profiles_ad")
      !  loop over all obs.
      profileCount = 0 
      
      ! loop over all obs.
      obs_loop: do tovsIndex = 1, nobmax
        if (tvs_lsensor(tovsIndex)/=sensorIndex) cycle obs_loop
        headerIndex = tvs_headerIndex(tovsIndex)
        profileCount = profileCount + 1
        
        TTb => col_getColumn(columng,headerIndex,'TT')
        HUb => col_getColumn(columng,headerIndex,'HU')
        Pres => col_getColumn(columng,headerIndex,'P_T')
        do levelIndex = 1, nlv_T
          modelTemperature(levelIndex,profileCount) = TTb(levelIndex)
          modelHumidity(levelIndex,profileCount) = HUb(levelIndex)
          modelPressureLevels(levelIndex,profileCount) = Pres(levelIndex) * MPC_MBAR_PER_PA_R8
          dPdPs(levelIndex,profileCount)  = 1.0d0
        end do
        
        ! Fix pour eviter probleme au toit avec GEM 4
        ! (grosse variabilite de la temperature au dernier niveau thermo due 
        !  a l'extrapolation utilisee)
        if (diagTtop) then
          modelTemperature(1,profileCount) =  modelTemperature(2,profileCount) + tvs_mesosphereLapseRate *  &
               log( col_getPressure(columng,1,headerIndex,'TH') /  &
               col_getPressure(columng,2,headerIndex,'TH') )
          modelHumidity(1,profileCount) =  modelHumidity(2,profileCount)
        end if

        sensorHeaderIndexes(profileCount) = headerIndex
      end do obs_loop
     
      !  2.1  Calculate the actual number of threads which will be used.

      nthreads = min(max_nthreads, profileCount )  

      !  2.2  Prepare all input variables required by rttov_ad.

      asw = 1

      call rttov_alloc_ad(                  &
           allocStatus(1),                  &
           asw,                             &
           profileCount,                    &
           btCount,                         &
           nRttovLevels,                         &
           chanprof,                        &
           opts=tvs_opts(sensorIndex),      &
           profiles_ad=profilesdataAd,     &
           coefs=tvs_coefs(sensorIndex),    &
           transmission= transmission,      &
           transmission_ad= transmissionAd,&
           radiance=radiancedata_d,         &
           radiance_ad=radiancedataAd,     &
           calcemis=calcemis,               &
           emissivity=emissivity_local,     &
           emissivity_ad=emissivityAd,     &
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

      emissivity_local(:)%emis_in = surfem1(:)
        
      do btIndex = 1, btCount
        bodyIndex = sensorBodyIndexes(btIndex)
        radiancedataAd % bt( btIndex ) = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
      end do


      !  2.3  Compute ad radiance with rttov_ad

      errorstatus  = 0
      emissivityAd(:) % emis_in = 0.0d0
      emissivityAd(:) % emis_out = 0.0d0
     
      call tmg_start(84,'rttov_ad')
      call rttov_parallel_ad(        &
           errorstatus,              & ! out
           chanprof,                 & ! in
           tvs_opts(sensorIndex),    & ! in
           tvs_profiles(sensorTovsIndexes(1:profileCount)), & ! in
           profilesdataAd,          & ! in
           tvs_coefs(sensorIndex),   & ! in
           transmission,             & ! inout
           transmissionAd,          & ! inout
           radiancedata_d,           & ! inout
           radiancedataAd,          & ! inout
           calcemis=calcemis,           & ! in
           emissivity=emissivity_local, & ! inout
           emissivity_ad=emissivityAd, & ! inout
           nthreads = nthreads )

      if (errorstatus /= 0) then
        Write(*,*) "Error in rttov_parallel_ad", errorstatus
        call utl_abort('tvslin_rttov_ad')
      end if

      call tmg_stop(84)

      !.. store results from rttov_ad into profilesAd
      extrapolatedTemperatureAd(:,:) = 0.d0
      extrapolatedHumidityAd(:,:) = 0.d0

      do btIndex = 1, btCount
        
        profileIndex = chanprof(btIndex)%prof
        headerIndex = sensorHeaderIndexes(profileIndex)

        ps_column => col_getColumn(column,headerIndex,'P0')
        p_column  => col_getColumn(column,headerIndex,'P_T')
        tg_column => col_getColumn(column,headerIndex,'TG')
        tt_column => col_getColumn(column,headerIndex,'TT')
        hu_column => col_getColumn(column,headerIndex,'HU')
        uu_column => col_getColumn(column,headerIndex,'UU')
        vv_column => col_getColumn(column,headerIndex,'VV')

        extrapolatedTemperatureAd(:,profileIndex) =  profilesdataAd(profileIndex) % t(:)
        extrapolatedHumidityAd(:,profileIndex) = profilesdataAd(profileIndex) % q(:)
        tg_column(1) = profilesdataAd(profileIndex) % skin % t 
        tt_column(ilowlvl_T) = profilesdataAd(profileIndex) % s2m % t
        ps_column(1) = profilesdataAd(profileIndex) % s2m % p * MPC_MBAR_PER_PA_R8
        hu_column(ilowlvl_T) = 0.d0 
        uu_column(ilowlvl_M) = profilesdataAd(profileIndex) % s2m % u
        vv_column(ilowlvl_M) = profilesdataAd(profileIndex) % s2m % v
      end do


      !  2.4  Adjoint of filling profilesAd structure
      do profileIndex = 1, profileCount
        extrapolatedHumidity(:,profileIndex) =  tvs_profiles(sensorTovsIndexes(profileIndex)) % q(:)
      end do

    
      !   2.3  Adjoint of extrapolation of humidity profile (kg/kg)
      !        above rlimlvhu (normally 300mbs or 70mbs)
      if ( TopAt10hPa ) then
        if ( tvs_debug ) then
          do profileIndex = 1, profileCount
            write(*,*)'extrapolatedHumidityAd*1000 avant aexthum4    = '
            write(*,'(1x,10f8.4)')(extrapolatedHumidityAd(levelIndex,profileIndex)*1000.d0,levelIndex=1,nRttovLevels)
            write(*,*)' '
          end do
        end if
        call aexthum4 (rttovPressureLevels(1:nRttovLevels),extrapolatedHumidityAd,extrapolatedHumidity)
        if ( tvs_debug ) then
          do profileIndex = 1, profileCount
            write(*,*)'extrapolatedHumidityAd*1000 apres aexthum4    = '
            write(*,'(1x,10f8.4)')(extrapolatedHumidityAd(levelIndex,profileIndex)*1000.d0,levelIndex=1,nRttovLevels)
            write(*,*)' '
          end do
        end if
      end if

      ! adjoint of conversion lnq --> q
      interpolatedHumidityAd(:,:) = 0.0d0
      do profileIndex = 1, profileCount
        do levelIndex = 1, levelsBelowModelTop
          interpolatedHumidityAd(levelIndex,profileIndex) = extrapolatedHumidityAd(nRttovLevels-levelsBelowModelTop + levelIndex,profileIndex)
        end do
      end do

      !   2.2  Adjoint of extrapolation of temperature profile above 10mb
      interpolatedTemperatureAd(:,:) = 0.0d0
      if ( .not. TopAt10hPa ) then
        do profileIndex = 1, profileCount
          interpolatedTemperatureAd(1:levelsBelowModelTop,profileIndex) = interpolatedTemperatureAd(1:levelsBelowModelTop,profileIndex) + &
               extrapolatedTemperatureAd(modelTopIndex:nRttovLevels,profileIndex)
        end do
      else
        call aextrap (interpolatedTemperatureAd,extrapolatedTemperatureAd,levelsBelowModelTop,nRttovLevels,profileCount)
      end if
 
      !   2.1  Adjoint of vertical interpolation of model temperature and logarithm of
      !        specific humidity to pressure levels required by tovs rt model
      
      modelTemperatureAd(:,:) = 0.0d0
      modelHumidityAd(:,:) = 0.0d0
      modelPressureAd(:,:) = 0.0d0

      call tmg_start(75,'intavgad')
      !$omp parallel do private(profileIndex)
      do profileIndex = 1, profileCount

        modelTemperatureAd(:, profileIndex) = 0.0d0
        call ppo_IntAvgAd_v2(modelPressureLevels(:,profileIndex:profileIndex), dPdPs(:,profileIndex:profileIndex), &
             modelTemperatureAd(:,profileIndex:profileIndex), modelTemperature(:,profileIndex:profileIndex), &
             modelPressureAd(:,profileIndex:profileIndex),nlv_T,1, &
             levelsBelowModelTop,rttovPressureLevels(modelTopIndex:nRttovLevels),interpolatedTemperatureAd(:,profileIndex:profileIndex))
        

        logModelHumidity(:,profileIndex) = log( modelHumidity(:,profileIndex) ) 
        logModelHumidityAd(:,profileIndex) = 0.d0
        loginterpolatedHumidityAd(:,profileIndex) = 0.d0
        loginterpolatedHumidityAd(:,profileIndex) = loginterpolatedHumidityAd(:,profileIndex) + interpolatedHumidityAd(:,profileIndex) * extrapolatedHumidity(modelTopIndex:nRttovLevels,profileIndex)
        call ppo_IntAvgAd_v2(modelPressureLevels(:,profileIndex:profileIndex),dPdPs(:,profileIndex:profileIndex), &
             logModelHumidityAd(:,profileIndex:profileIndex), logModelHumidity(:,profileIndex:profileIndex), &
             modelPressureAd(:,profileIndex:profileIndex),nlv_T,1, &
             levelsBelowModelTop,rttovPressureLevels(modelTopIndex:nRttovLevels), loginterpolatedHumidityAd(:,profileIndex:profileIndex))

        modelHumidityAd(:,profileIndex) = modelHumidityAd(:,profileIndex) + logModelHumidityAd(:,profileIndex) / modelHumidity(:,profileIndex)

      end do
      !$omp end parallel do
      call tmg_stop(75)
 
      ! Fix pour eviter probleme au toit avec GEM 4
      ! (grosse variabilite temperature au dernier niveau thermo due 
      !  a l'extrapolation utilisee)
      if ( diagTtop ) then
        do profileIndex = 1, profileCount
          modelTemperatureAd (1,profileIndex) = 0.d0
          modelHumidityAd(1,profileIndex) = 0.d0
        end do
      end if

      !   2.1  Store adjoints in columnData object

      do  profileIndex = 1 , profileCount 
        ps_column => col_getColumn(column, sensorHeaderIndexes(profileIndex), 'P0')
        p_column  => col_getColumn(column, sensorHeaderIndexes(profileIndex), 'P_T')
        tt_column => col_getColumn(column, sensorHeaderIndexes(profileIndex), 'TT')
        hu_column => col_getColumn(column, sensorHeaderIndexes(profileIndex), 'HU')
        
        do levelIndex = 1, col_getNumLev(column,'TH')
          p_column(levelIndex) = p_column(levelIndex)  + modelPressureAd  (levelIndex,profileIndex) * MPC_MBAR_PER_PA_R8
          tt_column(levelIndex) = tt_column(levelIndex) + modelTemperatureAd  (levelIndex,profileIndex)
          hu_column(levelIndex) = hu_column(levelIndex) + modelHumidityAd (levelIndex,profileIndex)
        end do
      end do

      deallocate (sensorHeaderIndexes,stat= allocStatus(1) )
      deallocate (interpolatedTemperatureAd,    stat= allocStatus(2) )
      deallocate (interpolatedHumidityAd,   stat= allocStatus(3) )
      deallocate (loginterpolatedHumidityAd,stat= allocStatus(4) )
      deallocate (extrapolatedTemperatureAd, stat= allocStatus(5) )
      deallocate (extrapolatedHumidityAd, stat= allocStatus(6) )
      deallocate (modelPressureLevels,    stat= allocStatus(7) )
      deallocate (dPdPs,    stat= allocStatus(8) )
      deallocate (modelTemperatureAd,    stat= allocStatus(9) )
      deallocate (modelHumidityAd,   stat= allocStatus(10))
      deallocate (logModelHumidityAd,stat= allocStatus(11))
      deallocate (modelTemperature,       stat= allocStatus(12))
      deallocate (modelHumidity,      stat= allocStatus(13))
      deallocate (logModelHumidity,   stat= allocStatus(14))
      deallocate (extrapolatedHumidity,    stat= allocStatus(15))
      deallocate (modelPressureAd,    stat= allocStatus(16))
      deallocate (rttovPressureLevels,    stat= allocStatus(17))
      
      call utl_checkAllocationStatus(allocStatus, " tvslin_fill_profiles_ad", .false.)
    
      !     de-allocate memory

      asw = 0 ! 0 to deallocate
      call rttov_alloc_ad(                  &
           allocStatus(1),                  &
           asw,                             &
           profileCount,                    &
           btCount,                         &
           nRttovLevels,                         &
           chanprof,                        &
           opts=tvs_opts(sensorIndex),      &
           profiles_ad=profilesdataAd,      &
           coefs=tvs_coefs(sensorIndex),    &
           transmission= transmission,      &
           transmission_ad= transmissionAd, &
           radiance=radiancedata_d,         &
           radiance_ad=radiancedataAd,      &
           calcemis=calcemis,               &
           emissivity=emissivity_local,     &
           emissivity_ad=emissivityAd )
     
      
      deallocate ( surfem1,           stat=allocStatus(2))
      deallocate ( sensorBodyIndexes, stat=allocStatus(3))
      call utl_checkAllocationStatus(allocStatus(1:3), " tvslin_rttov_ad", .false.)
     
    end do sensor_loop

    ! 3.  Close up

    deallocate ( sensorTovsIndexes )

  end subroutine tvslin_rttov_ad

  !--------------------------------------------------------------------------
  !  aexthum4
  !--------------------------------------------------------------------------
  subroutine aexthum4 (pressure, humidityAd, humidity )
    !
    ! :Purpose: adjoint of extrapolation of upper level humidity profile
    !           (adapted from exthumad by J. Eyre).
    !           To extend mixing ratio profile into stratosphere in a reasonable way.
    !           Take top tropospheric mixing ratio (e.g. near 300 mb) and
    !           extrapolate with given fall off into lower stratosphere
    !           (e.g. to 70 mb). Constrain mixing ratio to be >= zwmin
    !           (e.g. 0.000003 kg/kg). In upper strat, mixing ratio = zwmin.
    ! 
    implicit none
    
    ! Arguments:
    real(8),intent(in)    :: pressure(:)         ! Pressure levels of atm. profiles
    real(8),intent(inout) :: humidityAd(:,: ) ! Adjoint of humidity profiles.
    real(8),intent(in)    :: humidity(:,:)    ! Humidity profiles.

    ! locals:
    integer :: nlevels
    integer :: nprofiles
    real(8), allocatable :: zpres3(:)
    real(8)              :: zwmix, zwb
    real(8), parameter   :: pressureLimit = 70.0d0 ! press limits (in hpa) of region to be extrapolated
    integer              :: topIndex, profileIndex, levelIndex


    nlevels = size( pressure )
    nprofiles = size( humidity, dim =2)
    allocate ( zpres3( nlevels ) )
    
    !  find top level of given profile
    topIndex = -1
    do levelIndex = nlevels, 1, -1
      if (pressure(levelIndex) < FILT_RLIMLVHU) then
        topIndex = levelIndex
        exit
      end if
    end do
 
    !  Null extrapolation case
 
    if (topIndex == -1) return

    !    Constants defining p**3 fall off around tropopause
    do levelIndex = 1, topIndex
      zpres3(levelIndex) = ( pressure(levelIndex) / pressure(topIndex+1) )**3
    end do

    do profileIndex = 1, nprofiles
      zwb = 0.d0
      do levelIndex = 1, topIndex
        zwmix = humidityAd(levelIndex,profileIndex)
        humidityAd(levelIndex,profileIndex) = 0.d0
        if (pressure(levelIndex) >= pressureLimit) then
          if (humidity(levelIndex,profileIndex) > MPC_MINIMUM_HU_R8) then
            zwb = zwb + zwmix * zpres3(levelIndex)
          end if
        end if
      end do
      humidityAd(topIndex + 1,profileIndex) = humidityAd(topIndex + 1,profileIndex) + zwb
    end do

    deallocate ( zpres3 )

  end subroutine aexthum4

 
end module tovs_lin_mod


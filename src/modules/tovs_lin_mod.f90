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
  subroutine lexthum4(ppres, humidityTl, humidity)
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
    real(8),intent(in)    :: ppres(:)        ! Pressure levels of atm. profiles.
    real(8),intent(inout) :: humidityTl(:,:) ! Tl of humidity profiles.
    real(8),intent(in)    :: humidity(:,:)   ! Humidity profiles.

    !Locals:
    integer :: nlevels
    integer :: nprofiles

    REAL(8),allocatable :: zpres3( : )

    real(8) zwb5,zwb
    real(8) ,parameter :: zp1 = 70.0d0  ! press limits (in hpa) of region  to be extrapolated
    integer :: topIndex, profileIndex, levelIndex

    nlevels = size( ppres )
    nprofiles = size( humidity, dim =2)


    !    extrapolate humidity profile.

    !    find top level of given profile
    topIndex = -1
    do levelIndex=nlevels,1,-1
      if (ppres(levelIndex) < FILT_RLIMLVHU) then
        topIndex = levelIndex
        exit
      end if
    end do

    ! Null extrapolation case
    if (topIndex == -1) return
    allocate ( zpres3( nlevels ) )

    ! Constants defining p**3 fall off around tropopause
    do levelIndex=1, topINdex
      zpres3(levelIndex)=(ppres(levelIndex) / ppres(topIndex + 1)) ** 3
    end do

    do profileIndex=1,nprofiles
      zwb =humidityTl (topIndex + 1,profileIndex)
      zwb5=humidity(topIndex + 1,profileIndex)
      do levelIndex=1,topIndex
        if (ppres(levelIndex) < zp1) then
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
  subroutine tvslin_rttov_tl(column, columng, lobsSpaceData)
    !
    ! :Purpose: Tangent linear of computation of radiance with rttov_tl
    !   
    implicit none

    ! Arguments:
    type(struct_obs)        :: lobsSpaceData ! obsSpaceData structure
    type(struct_columnData) :: column        ! column structure for pertubation profile
    type(struct_columnData) :: columng       ! column structure for background profile

    ! Locals:
    type(struct_vco), pointer :: vco_anl
    integer, allocatable :: iptobs    (:) 
    integer, allocatable :: iptobs_header (:) 
    integer :: allocStatus(17)
    logical :: diagTtop,TopAt10hPa
    integer :: nobmax
    integer :: sensorIndex, tovsIndex
    integer :: ilowlvl_M,ilowlvl_T,profileCount,headerIndex,nlv_M,nlv_T
    integer :: levelIndex, profileIndex
    integer :: status, Vcode
    integer :: jpmotop, jpmolev

    real(8), allocatable :: to_tl    (:,:)
    real(8), allocatable :: huo_tl   (:,:)
    real(8), allocatable :: loghuo_tl(:,:)
    real(8), allocatable :: toext_tl (:,:)
    real(8), allocatable :: qoext_tl (:,:)
    real(8), allocatable :: zvlev    (:,:)
    real(8), allocatable :: dPdPs    (:,:)
    real(8), allocatable :: zt_tl    (:,:)
    real(8), allocatable :: zhu_tl   (:,:)
    real(8), allocatable :: logzhu_tl(:,:)
    real(8), allocatable :: zt       (:,:)
    real(8), allocatable :: zhu      (:,:)
    real(8), allocatable :: logzhu   (:,:)
    real(8), allocatable :: qoext    (:,:)
    real(8), allocatable :: zp_tl    (:,:)
    real(8), allocatable :: xpres    (:)
    real(8) :: zptop, zptopmbs
    real(8), pointer :: delTT(:), delHU(:), TTb(:), HUb(:), Pres(:), delP(:)
    integer :: nlevels
    integer :: btCount
    integer,external :: omp_get_num_threads
    integer :: nthreads, max_nthreads
    integer :: btIndex, bodyIndex
    integer :: instrum
    integer :: sensor_type   !sensor type(1=infrared; 2=microwave; 3=high resolution, 4=polarimetric)
    integer :: errorstatus
    integer,allocatable :: iptobs_body  (:)
    real(8), allocatable :: surfem1       (:) 
    type(rttov_emissivity), pointer :: emissivity_local (:)
    type(rttov_emissivity), pointer :: emissivity_tl (:)
    type(rttov_radiance) :: radiancedata_d   ! radiances full structure buffer used in rttov calls
    type(rttov_radiance) :: radiancedata_tl  ! tl radiances full structure buffer used in rttov calls
    type(rttov_transmission) :: transmission       ! transmission
    type(rttov_transmission) :: transmission_tl    ! transmission tl
    type(rttov_profile)  , pointer :: profilesdata_tl (:) ! tl profiles buffer used in rttov calls
    type(rttov_chanprof) , pointer :: chanprof (:)
    logical              :: init
    logical, pointer :: calcemis  (:)
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
    status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode)
    diagTtop = (Vcode==5002)
    
    ! find model level top, within 0.000001 mbs.
    zptop    = col_getPressure(columng,1,1,'TH')
    zptopmbs = zptop / 100.d0
    zptopmbs = zptopmbs - 0.000001d0

    TopAt10hPa = ( abs( zptopmbs - 10.0d0 ) <= .1d0 )

    !  1.  Get number of threads available and allocate memory for some variables

    !$omp parallel 
    max_nthreads = omp_get_num_threads()
    !$omp end parallel

    allocStatus(:) = 0
    allocate ( iptobs(tvs_nobtov),stat=allocStatus(1) )
    call utl_checkAllocationStatus(allocStatus(1:1), " tvslin_rttov_tl iptobs")
    
    ! 2.  Computation of hx for tovs data only

    
    ! Loop over all sensors specified by user

    sensor_loop:  do sensorIndex = 1, tvs_nsensors
       
      nlevels = tvs_coefs(sensorIndex) % coef % nlevels
      sensor_type = tvs_coefs(sensorIndex) % coef % id_sensor
      instrum = tvs_coefs(sensorIndex) % coef % id_inst
      !  loop over all obs.
      profileCount = 0
      do tovsIndex = 1, tvs_nobtov
        !    Currently processed sensor?
        if ( tvs_lsensor(tovsIndex) == sensorIndex ) then
          profileCount = profileCount + 1
          iptobs(profileCount) = tovsIndex
        end if
      end do
     
      if (profileCount == 0) cycle sensor_loop
      nobmax = iptobs(profileCount)
      !     compute the number of calculated radiances for one call
      btCount = tvs_countRadiances(iptobs(1:profileCount), lobsSpaceData, &
           assim_flag_val_opt=obs_assimilated)
      if ( btCount == 0 ) cycle  sensor_loop
   
      allocate (xpres(nlevels))
      xpres = tvs_coefs(sensorIndex)% coef % ref_prfl_p
      jpmotop = 1
      do levelIndex = 2, nlevels
        if ( zptopmbs >= xpres(levelIndex - 1).and. &
             zptopmbs < xpres(levelIndex)       ) then
          jpmotop = levelIndex
          exit
        end if
      end do
      jpmolev = (nlevels - jpmotop + 1)

      allocate (iptobs_header (profileCount),      stat= allocStatus(1) )
      allocate (to_tl     (jpmolev,profileCount),  stat= allocStatus(2) )
      allocate (huo_tl    (jpmolev,profileCount),  stat= allocStatus(3) )
      allocate (loghuo_tl (jpmolev,profileCount),  stat= allocStatus(4) )
      allocate (toext_tl  (nlevels  ,profileCount),stat= allocStatus(5) )
      allocate (qoext_tl  (nlevels  ,profileCount),stat= allocStatus(6) )
      allocate (zvlev     (nlv_T,profileCount),    stat= allocStatus(7) )
      allocate (dPdPs     (nlv_T,profileCount),    stat= allocStatus(8) )
      allocate (zt_tl     (nlv_T,profileCount),    stat= allocStatus(9) )
      allocate (zhu_tl    (nlv_T,profileCount),    stat= allocStatus(10))
      allocate (logzhu_tl (nlv_T,profileCount),    stat= allocStatus(11))
      allocate (zt        (nlv_T,profileCount),    stat= allocStatus(12))
      allocate (zhu       (nlv_T,profileCount),    stat= allocStatus(13))
      allocate (logzhu    (nlv_T,profileCount),    stat= allocStatus(14))
      allocate (qoext     (nlevels,profileCount),  stat= allocStatus(15))
      allocate (zp_tl     (nlv_T,profileCount),    stat= allocStatus(16))
      allocate (profilesdata_tl(profileCount),     stat= allocStatus(17))

      call utl_checkAllocationStatus(allocStatus, " tvslin_rttov_tl")
 
      iptobs_header(:) = 0 
      toext_tl (:,:) = 0.0d0
      zvlev    (:,:) = 0.0d0
      dPdPs    (:,:) = 0.0d0
      zt_tl    (:,:) = 0.0d0
      zhu_tl   (:,:) = 0.0d0
      zt       (:,:) = 0.0d0
      zhu      (:,:) = 0.0d0
      qoext    (:,:) = 0.0d0
      zp_tl    (:,:) = 0.0d0
      to_tl    (:,:) = 0.0d0
      huo_tl   (:,:) = 0.0d0

      ! allocate profiledata_tl structures
      init = .true.
      asw = 1 ! 1 to allocate
      call rttov_alloc_tl(                  &
              allocStatus(1),               &
              asw,                          &
              nprofiles=profileCount,       &
              nchanprof=btCount,            &
              nlevels=nlevels,              &
              chanprof=chanprof,            &
              opts=tvs_opts(sensorIndex),   &
              profiles_tl=profilesdata_tl,  &
              coefs=tvs_coefs(sensorIndex), &
              transmission=transmission,    &
              transmission_tl=transmission_tl, &
              radiance=radiancedata_d,      &
              radiance_tl=radiancedata_tl,  &
              calcemis=calcemis,            &
              emissivity=emissivity_local,  &
              emissivity_tl=emissivity_tl,  &
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
          zp_tl (levelIndex,profileCount) = delP(levelIndex) * MPC_MBAR_PER_PA_R8
          zt_tl (levelIndex,profileCount) = delTT(levelIndex)
          zhu_tl(levelIndex,profileCount) = delHU(levelIndex)
          zt   (levelIndex,profileCount)  = TTb(levelIndex)
          zhu  (levelIndex,profileCount)  = HUb(levelIndex)
          zvlev(levelIndex,profileCount)  = Pres(levelIndex) *MPC_MBAR_PER_PA_R8
          dPdPs(levelIndex,profileCount)  = 1.0d0
        end do
        
        ! Fix pour eviter probleme au toit avec GEM 4
        ! (grosse varibilite temperature au dernier niveau thermo due 
        !  a l'extrapolation utilisee)
        if ( diagTtop ) then
          zt_tl   (1,profileCount) =  0.d0
          zhu_tl  (1,profileCount) =  0.d0
          zt   (1,profileCount) =  zt   (2,profileCount) + tvs_mesosphereLapseRate *  &
               log( col_getPressure(columng,1,headerIndex,'TH') /  &
               col_getPressure(columng,2,headerIndex,'TH') )
          zhu  (1,profileCount) =  zhu  (2,profileCount)
        end if
        
        iptobs_header(profileCount) = headerIndex
        
      end do obs_loop

      !  2.1  Vertical interpolation of model temperature and logarithm of
      !       specific humidity to pressure levels required by tovs rt model

      do profileIndex = 1, profileCount
        qoext(1:nlevels,profileIndex) =  tvs_profiles(iptobs(profileIndex)) % q(1:nlevels)
      end do
      
      to_tl (:,:) = 0.0d0
      huo_tl(:,:) = 0.0d0
      
      !$omp parallel do private(profileIndex)
      do profileIndex=1, profileCount

        call ppo_IntAvgTl_v2(zvlev(:,profileIndex:profileIndex),dPdPs(:,profileIndex:profileIndex),zt_tl(:,profileIndex:profileIndex), &
             zt(:,profileIndex:profileIndex),zp_tl(:,profileIndex:profileIndex),nlv_T,1, &
             jpmolev,xpres(jpmotop:nlevels),to_tl(:,profileIndex:profileIndex))

        logzhu(:,profileIndex) = log( zhu(:,profileIndex) )
        logzhu_tl(:,profileIndex) = zhu_tl(:,profileIndex) / zhu(:,profileIndex)
        call ppo_IntAvgTl_v2(zvlev(:,profileIndex:profileIndex),dPdPs(:,profileIndex:profileIndex),logzhu_tl(:,profileIndex:profileIndex), &
             logzhu(:,profileIndex:profileIndex),zp_tl(:,profileIndex:profileIndex),nlv_T,1, &
             jpmolev,xpres(jpmotop:nlevels),loghuo_tl(:,profileIndex:profileIndex))

        huo_tl(:,profileIndex) = loghuo_tl(:,profileIndex) * qoext(jpmotop:nlevels,profileIndex)


      end do
      !$omp end parallel do
      
      !  2.2  Extrapolation of temperature profile above 10mb
      toext_tl(:,:) = 0.0d0
      if ( .not. TopAt10hPa ) then
        do profileIndex = 1, profileCount
          toext_tl(jpmotop:nlevels,profileIndex) = to_tl(1:jpmolev,profileIndex)
          toext_tl(1:jpmotop-1,profileIndex) = 0.d0
        end do
      else
        call lextrap (to_tl,toext_tl,jpmolev,nlevels,profileCount)
      end if
      
      !   2.3  Extrapolation of humidity profile (kg/kg)
      !        above rlimlvhu (normally 300mbs or 70mbs)
      
      do profileIndex = 1, profileCount
        do levelIndex = 1, jpmotop - 1
          qoext_tl(levelIndex,profileIndex) = 0.d0
        end do
        do levelIndex = 1, jpmolev
          qoext_tl(nlevels - jpmolev + levelIndex,profileIndex) = huo_tl(levelIndex,profileIndex)
        end do
      end do
      
      if ( TopAt10hPa ) then
        if ( tvs_debug ) then
          do profileIndex = 1, profileCount
            write(*,*)'qoext_tl*1000 avant exthum4    = '
            write(*,'(1x,10f8.4)')(qoext_tl(levelIndex,profileIndex) * 1000.d0,levelIndex=1,nlevels)
            write(*,*)' '
          end do
        end if
        call lexthum4 (xpres(1:nlevels),qoext_tl,qoext)
        if ( tvs_debug ) then
          do profileIndex = 1, profileCount
            write(*,*)'qoext_tl*1000 apres exthum4    = '
            write(*,'(1x,10f8.4)')(qoext_tl(levelIndex,profileIndex) * 1000.d0,levelIndex=1,nlevels)
            write(*,*)' '
          end do
        end if
      end if

      do  profileIndex = 1 , profileCount
        profilesdata_tl(profileIndex) % gas_units       = gas_unit_specconc ! all gas profiles should be provided in kg/kg
        profilesdata_tl(profileIndex) % nlevels         =  nlevels
        profilesdata_tl(profileIndex) % nlayers         =  nlevels - 1
        if (tvs_coefs(sensorIndex)%coef%nozone > 0) profilesdata_tl(profileIndex) % o3(:) =  0.d0
        
        profilesdata_tl(profileIndex) % ctp             = 0.0d0
        profilesdata_tl(profileIndex) % cfraction       = 0.0d0
        profilesdata_tl(profileIndex) % zenangle        = 0.0d0
        profilesdata_tl(profileIndex) % azangle         = 0.0d0
        profilesdata_tl(profileIndex) % skin % surftype = 0
        profilesdata_tl(profileIndex) % skin % t        = col_getElem(column,1,iptobs_header(profileIndex),'TG')
        profilesdata_tl(profileIndex) % skin % fastem(:)= 0.0d0
        profilesdata_tl(profileIndex) % skin % salinity = 0.0d0
        profilesdata_tl(profileIndex) % s2m % t         = col_getElem(column,ilowlvl_T,iptobs_header(profileIndex),'TT')

        
        profilesdata_tl(profileIndex) % s2m % q         = 0.d0

        profilesdata_tl(profileIndex) % s2m % p         = col_getElem(column,1,iptobs_header(profileIndex),'P0')*MPC_MBAR_PER_PA_R8
        profilesdata_tl(profileIndex) % s2m % u         = col_getElem(column,ilowlvl_M,iptobs_header(profileIndex),'UU')
        profilesdata_tl(profileIndex) % s2m % v         = col_getElem(column,ilowlvl_M,iptobs_header(profileIndex),'VV')
        
        profilesdata_tl(profileIndex) % p(1:nlevels)    = 0.d0
        profilesdata_tl(profileIndex) % t(1:nlevels)    = toext_tl(1:nlevels,profileIndex)
        profilesdata_tl(profileIndex) % q(1:nlevels)    = qoext_tl(1:nlevels,profileIndex)
      end do

      deallocate (iptobs_header, stat= allocStatus(1) )
      deallocate (to_tl,         stat= allocStatus(2) )
      deallocate (huo_tl,        stat= allocStatus(3) )
      deallocate (loghuo_tl,     stat= allocStatus(4) )
      deallocate (toext_tl,      stat= allocStatus(5) )
      deallocate (qoext_tl,      stat= allocStatus(6) )
      deallocate (zvlev,         stat= allocStatus(7) )
      deallocate (dPdPs,         stat= allocStatus(8) )
      deallocate (zt_tl,         stat= allocStatus(9) )
      deallocate (zhu_tl,        stat= allocStatus(10))
      deallocate (logzhu_tl,     stat= allocStatus(11))
      deallocate (zt,            stat= allocStatus(12))
      deallocate (zhu,           stat= allocStatus(13))
      deallocate (logzhu,        stat= allocStatus(14))
      deallocate (qoext,         stat= allocStatus(15))
      deallocate (zp_tl,         stat= allocStatus(16))
      deallocate (xpres,         stat= allocStatus(17))
      call utl_checkAllocationStatus(allocStatus, "tvslin_rttov_tl", .false.)

      !     set nthreads to actual number of threads which will be used.

      nthreads = min(max_nthreads, profileCount)  

      !   2.2  Prepare all input variables required by rttov.
     
      allocate ( surfem1      (btCount)       ,stat=allocStatus(1))
      allocate ( iptobs_body   (btCount)      ,stat=allocStatus(2))
      call utl_checkAllocationStatus(allocStatus(1:2), " tovs_rtttov_tl")
    
      !    get Hyperspecral IR emissivities
      if ( tvs_isInstrumHyperSpectral(instrum) ) call tvs_getHIREmissivities(iptobs(1:profileCount), &
           lobsSpaceData, surfem1)

      call tvs_getChanprof(sensorIndex, iptobs(1:profileCount), lobsSpaceData, chanprof, &
           iptobs_cma_opt=iptobs_body)

      call tvs_getOtherEmissivities(chanprof, iptobs, sensor_type, instrum, surfem1, calcemis)
      
      emissivity_local(:)%emis_in = surfem1(:)

 
      !  2.3  Compute tl radiance with rttov_tl
      
      errorstatus   = 0
      emissivity_tl(:)%emis_in = 0.0d0
      call tmg_start(87,'rttov_tl')
      call rttov_parallel_tl(         &
           errorstatus,               & ! out
           chanprof,                  & ! in
           tvs_opts(sensorIndex),     & ! in
           tvs_profiles(iptobs(1:profileCount)),  & ! in
           profilesdata_tl,           & ! inout
           tvs_coefs(sensorIndex),    & ! in
           transmission,              & ! inout
           transmission_tl,           & ! inout
           radiancedata_d,            & ! inout
           radiancedata_tl,           & ! inout
           calcemis,                  & ! in
           emissivity_local,          & ! in
           emissivity_tl,             & ! inout
           nthreads=nthreads )          ! in
               
      if (errorstatus /= 0) then
        Write(*,*) "Error in rttov_parallel_tl",errorstatus
        write(*,*) 'temperature           profile=',tvs_profiles(iptobs(1)) % t(:)
        write(*,*) 'temperature increment profile=',profilesdata_tl(1) % t(:)
        call utl_abort('tovs_rttov_tl')
      end if

      call tmg_stop(87)

      !  2.4  Store hx in obsSpaceData,OBS_WORK
      
      do btIndex = 1, btCount
        
        bodyIndex = iptobs_body(btIndex)
        call obs_bodySet_r(lobsSpaceData,OBS_WORK,bodyIndex, &
             radiancedata_tl % bt(btIndex) )
        if ( tvs_debug ) then
          write(*,'(a,i4,2f8.2)') ' ichn,sim,obs= ', &
               chanprof(btIndex)%chan,   radiancedata_tl % bt(btIndex), &
               obs_bodyElem_r(lobsSpaceData,OBS_OMP,bodyIndex)
        end if
        
      end do

 
      ! de-allocate memory
      
      asw = 0 ! 0 to deallocate
      call rttov_alloc_tl(                   &
           allocStatus(1),                   &
           asw,                              &
           nprofiles=profileCount,           &
           nchanprof=btCount,                &
           nlevels=nlevels,                  &
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

      deallocate ( surfem1,    stat=allocStatus(2) )
      deallocate ( iptobs_body,stat=allocStatus(3) )
      call utl_checkAllocationStatus(allocStatus(1:3), " tvslin_rtttov_tl", .false.)
      
    end do sensor_loop

    deallocate ( iptobs )


  end subroutine tvslin_rttov_tl

  !--------------------------------------------------------------------------
  !  tvslin_rttov_ad
  !--------------------------------------------------------------------------
  subroutine tvslin_rttov_ad( column, columng, lobsSpaceData )
    !
    ! :Purpose: Adjoint of computation of radiance with rttov_ad
    !

    implicit none

    ! Arguments:
    type(struct_columnData) :: column
    type(struct_columnData) :: columng
    type(struct_obs)        :: lobsSpaceData

    ! locals
    type(struct_vco), pointer :: vco_anl
    integer, allocatable :: iptobs(:) 
    integer, allocatable :: iptobs_header(:) 
    
    integer :: allocStatus(17)
    logical :: diagTtop,TopAt10hPa
    integer :: omp_get_num_threads, nthreads
    integer :: nlevels,nobmax
    integer :: sensorIndex, tovsIndex
    integer :: ilowlvl_T,ilowlvl_M,profileCount,headerIndex,nlv_M,nlv_T
    integer :: profileIndex, levelIndex
    integer :: status, Vcode
    integer :: jpmotop, jpmolev
    
    real(8), allocatable :: to_ad    (:,:)
    real(8), allocatable :: huo_ad   (:,:)
    real(8), allocatable :: loghuo_ad(:,:)
    real(8), allocatable :: toext_ad (:,:)
    real(8), allocatable :: qoext_ad (:,:)
    real(8), allocatable :: zvlev    (:,:)
    real(8), allocatable :: dPdPs    (:,:)
    real(8), allocatable :: zt_ad    (:,:)
    real(8), allocatable :: zhu_ad   (:,:)
    real(8), allocatable :: logzhu_ad(:,:)
    real(8), allocatable :: zt       (:,:)
    real(8), allocatable :: zhu      (:,:)
    real(8), allocatable :: logzhu   (:,:)
    real(8), allocatable :: qoext    (:,:)
    real(8), allocatable :: zp_ad    (:,:)
    real(8), allocatable :: xpres    (:)

    real(8) :: zptop, zptopmbs
   
    real(8), pointer :: uu_column(:),vv_column(:),tt_column(:),hu_column(:),ps_column(:),tg_column(:),p_column(:)
    real(8), pointer :: TTb(:), HUb(:), Pres(:)

    integer :: btCount
    integer :: max_nthreads
    integer :: instrum
    integer :: btIndex, bodyIndex
    integer :: sensor_type   ! sensor type (1=infrared; 2=microwave; 3=high resolution, 4=polarimetric)
    
    integer, allocatable :: iptobs_body(:)
    integer :: errorstatus 
    
    real(8), allocatable :: surfem1(:) 
    type(rttov_emissivity), pointer :: emissivity_local (:)
    type(rttov_emissivity), pointer :: emissivity_ad (:)
    type(rttov_transmission) :: transmission,transmission_ad
    type(rttov_radiance) :: radiancedata_ad, radiancedata_d
    
    type(rttov_profile), pointer  :: profilesdata_ad (:) ! ad profiles buffer used in rttov calls
    type(rttov_chanprof), pointer :: chanprof(:)
    logical :: init
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
    zptop    = col_getPressure(columng,1,1,'TH')
    zptopmbs = zptop / 100.d0
    zptopmbs = zptopmbs - 0.000001d0

    TopAt10hPa = ( abs( zptopmbs - 10.0d0 ) <= .1d0 )


    !     1.  Get number of threads available and allocate memory for some variables
 
    !$omp parallel 
    max_nthreads = omp_get_num_threads()
    !$omp end parallel

    allocate ( iptobs(tvs_nobtov) )

    !     2.  Computation of adjoint hx for tovs data only

    ! Loop over all sensors specified by user

    sensor_loop:do  sensorIndex = 1, tvs_nsensors
     
      nlevels = tvs_coefs(sensorIndex) %coef % nlevels
      sensor_type = tvs_coefs(sensorIndex) % coef% id_sensor
      instrum = tvs_coefs(sensorIndex) % coef% id_inst

      profileCount = 0
      do tovsIndex = 1, tvs_nobtov
        !    Currently processed sensor?
        if ( tvs_lsensor(tovsIndex) == sensorIndex ) then
          profileCount = profileCount + 1
          iptobs(profileCount) = tovsIndex
        end if
      end do

      if (profileCount == 0) cycle sensor_loop

      nobmax = iptobs(profileCount)

      !  compute the number of radiances/tbs to be calculated
      btCount = tvs_countRadiances(iptobs(1:profileCount), lobsSpaceData)

      if (btCount == 0) cycle sensor_loop

     
      allocate (xpres(nlevels))
      xpres = tvs_coefs(sensorIndex)% coef % ref_prfl_p
      jpmotop = 1
      do levelIndex = 2, nlevels
        if ( zptopmbs >= xpres(levelIndex - 1 ) .and.    &
             zptopmbs < xpres(levelIndex)        ) then
          jpmotop = levelIndex
          exit
        end if
      end do

      jpmolev = (nlevels - jpmotop + 1)
     
      allocStatus(:) = 0
      allocate (iptobs_header(profileCount)    ,stat= allocStatus(1) )
      allocate (to_ad    (jpmolev,profileCount),stat= allocStatus(2) )
      allocate (huo_ad   (jpmolev,profileCount),stat= allocStatus(3) )
      allocate (loghuo_ad(jpmolev,profileCount),stat= allocStatus(4) )
      allocate (toext_ad (nlevels,profileCount),stat= allocStatus(5) )
      allocate (qoext_ad (nlevels,profileCount),stat= allocStatus(6) )
      allocate (zvlev    (nlv_T,profileCount)  ,stat= allocStatus(7) )
      allocate (dPdPs    (nlv_T,profileCount)  ,stat= allocStatus(8) )
      allocate (zt_ad    (nlv_T,profileCount)  ,stat= allocStatus(9) )
      allocate (zhu_ad   (nlv_T,profileCount)  ,stat= allocStatus(10))
      allocate (logzhu_ad(nlv_T,profileCount)  ,stat= allocStatus(11))
      allocate (zt       (nlv_T,profileCount)  ,stat= allocStatus(12))
      allocate (zhu      (nlv_T,profileCount)  ,stat= allocStatus(13))
      allocate (logzhu   (nlv_T,profileCount)  ,stat= allocStatus(14))
      allocate (qoext    (nlevels,profileCount),stat= allocStatus(15))
      allocate (zp_ad    (nlv_T,profileCount)  ,stat= allocStatus(16))

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
          zt   (levelIndex,profileCount) = TTb(levelIndex)
          zhu  (levelIndex,profileCount) = HUb(levelIndex)
          zvlev(levelIndex,profileCount) = Pres(levelIndex) * MPC_MBAR_PER_PA_R8
          dPdPs(levelIndex,profileCount)  = 1.0d0
        end do
        
        ! Fix pour eviter probleme au toit avec GEM 4
        ! (grosse variabilite de la temperature au dernier niveau thermo due 
        !  a l'extrapolation utilisee)
        if (diagTtop) then
          zt   (1,profileCount) =  zt   (2,profileCount) + tvs_mesosphereLapseRate *  &
               log( col_getPressure(columng,1,headerIndex,'TH') /  &
               col_getPressure(columng,2,headerIndex,'TH') )
          zhu  (1,profileCount) =  zhu  (2,profileCount)
        end if

        iptobs_header(profileCount) = headerIndex
      end do obs_loop
     
      !  2.1  Calculate the actual number of threads which will be used.

      nthreads = min(max_nthreads, profileCount )  

      !  2.2  Prepare all input variables required by rttov_ad.

      asw = 1
      init = .true.

      call rttov_alloc_ad(                  &
           allocStatus(1),                  &
           asw,                             &
           profileCount,                    &
           btCount,                         &
           nlevels,                         &
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
           init=init)

      allocate ( surfem1      (btCount)        ,stat=allocStatus(2))
      allocate ( iptobs_body   (btCount)       ,stat=allocStatus(3))
      call utl_checkAllocationStatus(allocStatus(1:3), " tvslin_rttov_ad")
      
      !  get Hyperspectral IR emissivities
      
      if ( tvs_isInstrumHyperSpectral(instrum) ) call tvs_getHIREmissivities(iptobs(1:profileCount), lobsSpaceData, surfem1)

      ! Build the list of channels/profiles indices

      call tvs_getChanprof(sensorIndex, iptobs(1:profileCount), lobsSpaceData, chanprof, iptobs_body)
      !     get non Hyperspectral IR emissivities
      call tvs_getOtherEmissivities(chanprof, iptobs, sensor_type, instrum, surfem1, calcemis)

      emissivity_local(:)%emis_in = surfem1(:)
        
      do btIndex = 1, btCount
        bodyIndex = iptobs_body(btIndex)
        radiancedata_ad % bt( btIndex ) = obs_bodyElem_r(lobsSpaceData,OBS_WORK,bodyIndex)
      end do


      !  2.3  Compute ad radiance with rttov_ad

      errorstatus  = 0
      emissivity_ad(:) % emis_in = 0.0d0
      emissivity_ad(:) % emis_out = 0.0d0
     
      call tmg_start(84,'rttov_ad')
      call rttov_parallel_ad(        &
           errorstatus,              & ! out
           chanprof,                 & ! in
           tvs_opts(sensorIndex),    & ! in
           tvs_profiles(iptobs(1:profileCount)), & ! in
           profilesdata_ad,          & ! in
           tvs_coefs(sensorIndex),   & ! in
           transmission,             & ! inout
           transmission_ad,          & ! inout
           radiancedata_d,           & ! inout
           radiancedata_ad,          & ! inout
           calcemis=calcemis,           & ! in
           emissivity=emissivity_local, & ! inout
           emissivity_ad=emissivity_ad, & ! inout
           nthreads = nthreads )

      if (errorstatus /= 0) then
        Write(*,*) "Error in rttov_parallel_ad", errorstatus
        call utl_abort('tvslin_rttov_ad')
      end if

      call tmg_stop(84)

      !.. store results from rttov_ad into profiles_ad
      toext_ad(:,:) = 0.d0
      qoext_ad(:,:) = 0.d0

      do btIndex = 1, btCount
        
        profileIndex = chanprof(btIndex)%prof
        headerIndex = iptobs_header(profileIndex)

        ps_column => col_getColumn(column,headerIndex,'P0')
        p_column  => col_getColumn(column,headerIndex,'P_T')
        tg_column => col_getColumn(column,headerIndex,'TG')
        tt_column => col_getColumn(column,headerIndex,'TT')
        hu_column => col_getColumn(column,headerIndex,'HU')
        uu_column => col_getColumn(column,headerIndex,'UU')
        vv_column => col_getColumn(column,headerIndex,'VV')

        toext_ad(:,profileIndex)   =  profilesdata_ad(profileIndex) % t(:)
        qoext_ad(:,profileIndex)   =  profilesdata_ad(profileIndex) % q(:)
        tg_column(1)                =  profilesdata_ad(profileIndex) % skin % t 
        tt_column(ilowlvl_T)        =  profilesdata_ad(profileIndex) % s2m % t
        ps_column(1)                =  profilesdata_ad(profileIndex) % s2m % p * MPC_MBAR_PER_PA_R8
        hu_column(ilowlvl_T)        =  0.d0 
        uu_column(ilowlvl_M)        =  profilesdata_ad(profileIndex) % s2m % u
        vv_column(ilowlvl_M)        =  profilesdata_ad(profileIndex) % s2m % v
      end do


      !  2.4  Adjoint of filling profiles_ad structure
      do profileIndex = 1, profileCount
        qoext(:,profileIndex) =  tvs_profiles(iptobs(profileIndex)) % q(:)
      end do

    
      !   2.3  Adjoint of extrapolation of humidity profile (kg/kg)
      !        above rlimlvhu (normally 300mbs or 70mbs)
      if ( TopAt10hPa ) then
        if ( tvs_debug ) then
          do profileIndex = 1, profileCount
            write(*,*)'qoext_ad*1000 avant aexthum4    = '
            write(*,'(1x,10f8.4)')(qoext_ad(levelIndex,profileIndex)*1000.d0,levelIndex=1,nlevels)
            write(*,*)' '
          end do
        end if
        call aexthum4 (xpres(1:nlevels),qoext_ad,qoext)
        if ( tvs_debug ) then
          do profileIndex = 1, profileCount
            write(*,*)'qoext_ad*1000 apres aexthum4    = '
            write(*,'(1x,10f8.4)')(qoext_ad(levelIndex,profileIndex)*1000.d0,levelIndex=1,nlevels)
            write(*,*)' '
          end do
        end if
      end if

      ! adjoint of conversion lnq --> q
      huo_ad(:,:) = 0.0d0
      do profileIndex = 1, profileCount
        do levelIndex = 1, jpmolev
          huo_ad(levelIndex,profileIndex) = qoext_ad(nlevels-jpmolev + levelIndex,profileIndex)
        end do
      end do

      !   2.2  Adjoint of extrapolation of temperature profile above 10mb
      to_ad(:,:) = 0.0d0
      if ( .not. TopAt10hPa ) then
        do profileIndex = 1, profileCount
          to_ad(1:jpmolev,profileIndex) = to_ad(1:jpmolev,profileIndex) + &
               toext_ad(jpmotop:nlevels,profileIndex)
        end do
      else
        call aextrap (to_ad,toext_ad,jpmolev,nlevels,profileCount)
      end if
 
      !   2.1  Adjoint of vertical interpolation of model temperature and logarithm of
      !        specific humidity to pressure levels required by tovs rt model
      
      zt_ad (:,:) = 0.0d0
      zhu_ad(:,:) = 0.0d0
      zp_ad (:,:) = 0.0d0

      call tmg_start(75,'intavgad')
      !$omp parallel do private(profileIndex)
      do profileIndex = 1, profileCount

        zt_ad(:,profileIndex) = 0.0d0
        call ppo_IntAvgAd_v2(zvlev(:,profileIndex:profileIndex), dPdPs(:,profileIndex:profileIndex), &
             zt_ad(:,profileIndex:profileIndex), zt(:,profileIndex:profileIndex), &
             zp_ad(:,profileIndex:profileIndex),nlv_T,1, &
             jpmolev,xpres(jpmotop:nlevels),to_ad(:,profileIndex:profileIndex))
        

        logzhu(:,profileIndex) = log( zhu(:,profileIndex) ) 
        logzhu_ad(:,profileIndex) = 0.d0
        loghuo_ad(:,profileIndex) = 0.d0
        loghuo_ad(:,profileIndex) = loghuo_ad(:,profileIndex) + huo_ad(:,profileIndex) * qoext(jpmotop:nlevels,profileIndex)
        call ppo_IntAvgAd_v2(zvlev(:,profileIndex:profileIndex),dPdPs(:,profileIndex:profileIndex), &
             logzhu_ad(:,profileIndex:profileIndex), logzhu(:,profileIndex:profileIndex), &
             zp_ad(:,profileIndex:profileIndex),nlv_T,1, &
             jpmolev,xpres(jpmotop:nlevels), loghuo_ad(:,profileIndex:profileIndex))

        zhu_ad(:,profileIndex) = zhu_ad(:,profileIndex) + logzhu_ad(:,profileIndex) / zhu(:,profileIndex)

      end do
      !$omp end parallel do
      call tmg_stop(75)
 
      ! Fix pour eviter probleme au toit avec GEM 4
      ! (grosse variabilite temperature au dernier niveau thermo due 
      !  a l'extrapolation utilisee)
      if ( diagTtop ) then
        do profileIndex = 1, profileCount
          zt_ad (1,profileIndex) = 0.d0
          zhu_ad(1,profileIndex) = 0.d0
        end do
      end if

      !   2.1  Store adjoints in columnData object

      do  profileIndex = 1 , profileCount 
        ps_column => col_getColumn(column,iptobs_header(profileIndex),'P0')
        p_column  => col_getColumn(column,iptobs_header(profileIndex),'P_T')
        tt_column => col_getColumn(column,iptobs_header(profileIndex),'TT')
        hu_column => col_getColumn(column,iptobs_header(profileIndex),'HU')
        
        do levelIndex = 1, col_getNumLev(column,'TH')
          p_column(levelIndex) = p_column(levelIndex)  + zp_ad  (levelIndex,profileIndex) * MPC_MBAR_PER_PA_R8
          tt_column(levelIndex) = tt_column(levelIndex) + zt_ad  (levelIndex,profileIndex)
          hu_column(levelIndex) = hu_column(levelIndex) + zhu_ad (levelIndex,profileIndex)
        end do
      end do

      deallocate (iptobs_header,stat= allocStatus(1) )
      deallocate (to_ad    ,stat= allocStatus(2) )
      deallocate (huo_ad   ,stat= allocStatus(3) )
      deallocate (loghuo_ad,stat= allocStatus(4) )
      deallocate (toext_ad ,stat= allocStatus(5) )
      deallocate (qoext_ad ,stat= allocStatus(6) )
      deallocate (zvlev    ,stat= allocStatus(7) )
      deallocate (dPdPs    ,stat= allocStatus(8) )
      deallocate (zt_ad    ,stat= allocStatus(9) )
      deallocate (zhu_ad   ,stat= allocStatus(10))
      deallocate (logzhu_ad,stat= allocStatus(11))
      deallocate (zt       ,stat= allocStatus(12))
      deallocate (zhu      ,stat= allocStatus(13))
      deallocate (logzhu   ,stat= allocStatus(14))
      deallocate (qoext    ,stat= allocStatus(15))
      deallocate (zp_ad    ,stat= allocStatus(16))
      deallocate (xpres    ,stat= allocStatus(17))
      
      call utl_checkAllocationStatus(allocStatus, " tvslin_fill_profiles_ad", .false.)
    
      !     de-allocate memory

      asw = 0 ! 0 to deallocate
      call rttov_alloc_ad(                  &
           allocStatus(1),                  &
           asw,                             &
           profileCount,                    &
           btCount,                         &
           nlevels,                         &
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
     
      
      deallocate ( surfem1         ,stat=allocStatus(2))
      deallocate ( iptobs_body     ,stat=allocStatus(3))
      call utl_checkAllocationStatus(allocStatus(1:3), " tvslin_rttov_ad", .false.)
     
    end do sensor_loop

    ! 3.  Close up

    deallocate ( iptobs )

  end subroutine tvslin_rttov_ad

  !--------------------------------------------------------------------------
  !  aexthum4
  !--------------------------------------------------------------------------
  subroutine aexthum4 (ppres, humidityAd, humidity )
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
    real(8),intent(in)    :: ppres(:)         ! Pressure levels of atm. profiles
    real(8),intent(inout) :: humidityAd(:,: ) ! Adjoint of humidity profiles.
    real(8),intent(in)    :: humidity(:,:)    ! Humidity profiles.

    ! locals:
    integer :: nlevels
    integer :: nprofiles
    real(8),allocatable :: zpres3(:)
    real(8)             :: zwmix, zwb
    real(8), parameter  :: zp1 = 70.0d0 ! press limits (in hpa) of region to be extrapolated
    integer             :: topIndex, profileIndex, levelIndex


    nlevels = size( ppres )
    nprofiles = size( humidity, dim =2)
    allocate ( zpres3( nlevels ) )
    
    !  find top level of given profile
    topIndex = -1
    do levelIndex=nlevels,1,-1
      if (ppres(levelIndex) < FILT_RLIMLVHU) then
        topIndex = levelIndex
        exit
      end if
    end do
 
    !  Null extrapolation case
 
    if (topIndex == -1) return

    !    Constants defining p**3 fall off around tropopause
    do levelIndex=1,topIndex
      zpres3(levelIndex) = (ppres(levelIndex)/ppres(topIndex+1))**3
    end do

    do profileIndex=1,nprofiles
      zwb = 0.d0
      do levelIndex=1,topIndex
        zwmix = humidityAd(levelIndex,profileIndex)
        humidityAd(levelIndex,profileIndex) = 0.d0
        if (ppres(levelIndex) >= zp1) then
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


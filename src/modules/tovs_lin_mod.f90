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
!! MODULE tovs_lin (prefix= no defined prefix)
!!
!! *Purpose*: Derived types, public variables and procedures related to the 
!!            tangent-linear and adjoint versions of RTTOV
!!
!--------------------------------------------------------------------------
module tovs_lin_mod
  use rttov_interfaces_mod

  use rttov_types, only : rttov_profile, rttov_radiance
  use rttov_const ,only : gas_unit_specconc
  use parkind1, only : jpim, jprb

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

  subroutine LEXTHUM4(KNPF,KLAPF,PPRES,PAV,PAV5)
!
!**** *lexthum4* - tangent linear of extrapolaation of upper level humidity profile.
!                 (adapted from exthumtl by J. Eyre)
!
!     purpose.
!     --------
!          to extend mixing ratio profile into stratosphere in
!          a reasonable way.
!
!**   interface.
!     ----------
!          *call* *lexthum4(knpf,klapf,ppres,pav,pav5)*
!               *knpf*:  no. of profiles to be processed.
!               *klapf*: length of atm. profiles.
!               *ppres*: pressure levels of atm. profiles.
!               *pav*:   gradient humidity profiles.
!               *pav5*:  humidity profiles.
!
!     method.
!     -------
!          take top tropospheric mixing ratio (e.g. near 300 mb) and
!          extrapolate with given fall off into lower stratosphere
!          (e.g. to 70 mb).  constrain mixing ratio to be >= zwmin
!          (e.g. 0.000003 kg/kg).   in upper strat, mixing ratio = zwmin.
!
!     externals.
!     ----------
!          none.
!
!     reference.
!     ----------
!          ecmwf tech mem 176.
!

    implicit none

    integer :: klapf, knpf
    REAL(8) PPRES(*),PAV(KLAPF,*),PAV5(KLAPF,*)

    REAL(8) :: ZPRES3(KLAPF)

    real(8) zwb5,zwb
    real(8) ,parameter :: ZP1 = 70.0D0  ! PRESS LIMITS (IN HPA) OF REGION  to be extrapolated
    integer :: inlvw,j,jnpf,ierr

    !
    !
    !*         1.   extrapolate humidity profile.
    !               ----------- -------- -------

    !          find top level of given profile
    INLVW = -1
    do J=KLAPF,1,-1
      if (PPRES(J) < FILT_RLIMLVHU) then
        INLVW = J
        exit
      end if
    end do

    !** Null extrapolation case
    if (INLVW == -1) return

    !          constants defining p**3 fall off around tropopause
    do J=1,INLVW
      ZPRES3(J)=(PPRES(J) / PPRES(INLVW + 1)) ** 3
    end do

    do JNPF=1,KNPF
      ZWB =PAV (INLVW + 1,JNPF)
      ZWB5=PAV5(INLVW + 1,JNPF)
      do J=1,INLVW
        if (PPRES(J) < ZP1) then
          PAV(J,JNPF) = 0.D0
        else
          if ( ZWB5 * ZPRES3(J) <= MPC_MINIMUM_HU_R8 ) then
            PAV(J,JNPF) = 0.D0
          else
            PAV(J,JNPF) = ZWB * ZPRES3(J)
          end if
        end if
      end do
    end do
    
  end subroutine LEXTHUM4

  subroutine tvslin_rttov_tl(column, columng, lobsSpaceData, obs_ass_val)
!--------------------------------------------------------------------------
!! *Purpose*: Tangent linear of computation of radiance with rttov_tl
!!
!! @author j. halle *cmda/aes  april 19, 2005
!!
!
!revision 001  : a. beaulne *cmda/msc  june 2006
!                  - addition of ozone and IR surface emissivities
!revision 002  : r. sarrazin cmda   april 2008
!                  - adapt to CSR
!revision 003  : s. heilliette
!                  - adapt to IASI
!         S. heilliette:
!                  - adaptation to rttov 10.0 (october 2010)
!revision 004  : s. macpherson  nov 2012
!                  - remove #include "comtovst.cdk"
!
!--------------------------------------------------------------------------
   
    implicit none

    type(struct_obs) :: lobsSpaceData
    type(struct_columnData) :: column,columng
    integer ,intent(in) :: obs_ass_val
    type(struct_vco), pointer :: vco_anl
    integer, allocatable :: iptobs    (:) 
    integer, allocatable :: iptobs_header (:) 
    integer :: alloc_status(14)
    logical :: diagTtop,TopAt10hPa
    integer :: nobmax
    integer :: j, i, sensor_id, iobs, jj, stat
    integer :: ilowlvl_M,ilowlvl_T,count_profile,headerIndex,nlv_M,nlv_T
    integer :: jk, jn, jl
    integer :: status, Vcode
    integer :: jpmotop, jpmolev

    real(8), allocatable :: to_tl    (:,:)
    real(8), allocatable :: lqo_tl   (:,:)
    real(8), allocatable :: toext_tl (:,:)
    real(8), allocatable :: qoext_tl (:,:)
    real(8), allocatable :: zvlev    (:,:)
    real(8), allocatable :: dPdPs    (:,:)
    real(8), allocatable :: zt_tl    (:,:)
    real(8), allocatable :: zlq_tl   (:,:)
    real(8), allocatable :: zt       (:,:)
    real(8), allocatable :: zlq      (:,:)
    real(8), allocatable :: qoext    (:,:)
    real(8), allocatable :: zps_tl   (:)
    real(8), allocatable :: xpres    (:)
   
    real(8) :: zptop, zptopmbs
    real(8), pointer :: delTT(:), delLQ(:), TTb(:), HUb(:), Pres(:)

    integer :: isurface
    integer :: nlevels
    integer :: count_tb
    integer,external :: omp_get_num_threads
    integer :: nthreads, max_nthreads
    integer :: profile_index, index_tb, bodyIndex, ichn, istart, iend
    integer :: instrum
    integer :: sensor_type   !sensor type(1=infrared; 2=microwave; 3=high resolution, 4=polarimetric)
    integer :: channel_index

    integer :: errorstatus
    integer,allocatable :: iptobs_body  (:)

    real*8, allocatable :: surfem1       (:) 
    type(rttov_emissivity), pointer :: emissivity_local (:)
    type(rttov_emissivity), pointer :: emissivity_tl (:)
   
    type(rttov_radiance) :: radiancedata_d   ! radiances full structure buffer used in rttov calls
    type(rttov_radiance) :: radiancedata_tl  ! tl radiances full structure buffer used in rttov calls
    type(rttov_transmission) :: transmission       ! transmission
    type(rttov_transmission) :: transmission_tl    ! transmission tl
    type(rttov_profile)  , pointer :: profilesdata_tl (:) ! tl profiles buffer used in rttov calls
    type(rttov_chanprof) , pointer :: chanprof (:)
    logical ,save        :: first=.true.
    logical              :: init
    
    logical, pointer :: calcemis  (:)
    integer ::  asw
         
    if (tvs_nobtov == 0) return       ! exit if there are not tovs data

    !     1.    Set index for model's lowest level and model top
    !     .     ------------------------------------------------

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

    !     1.  Get number of threads available and allocate memory for some variables
    !     .   ---------------------------------------------------------------------- 

!$omp parallel 
    max_nthreads = omp_get_num_threads()
!$omp end parallel

    alloc_status = 0
    allocate ( iptobs(tvs_nobtov),stat=alloc_status(1) )
    call utl_checkAllocationStatus(alloc_status(1:1), " tvslin_rttov_tl iptobs")
    !
    !     2.  Computation of hx for tovs data only
    !     .   ------------------------------------

    
    ! Loop over all sensors specified by user

    sensor_loop:  do sensor_id = 1, tvs_nsensors
       
      nlevels = tvs_coefs(sensor_id)%coef % nlevels
      sensor_type = tvs_coefs(sensor_id) % coef % id_sensor
      instrum = tvs_coefs(sensor_id) % coef % id_inst
      !  loop over all obs.
      count_profile = 0
      do iobs = 1, tvs_nobtov
        !    Currently processed sensor?
        if ( tvs_lsensor(iobs) == sensor_id ) then
          count_profile = count_profile + 1
          iptobs(count_profile) = iobs
        end if
      end do
     
      if (count_profile == 0) cycle sensor_loop
      NOBMAX = iptobs(count_profile)
      !     compute the number of calculated radiances for one call
      count_tb = tvs_countRadiances(iptobs, count_profile, lobsSpaceData, &
           assim_flag_val_opt=obs_ass_val)
      if ( count_tb == 0 ) cycle  sensor_loop
   
      allocate (xpres(nlevels))
      xpres = tvs_coefs(sensor_id)% coef % ref_prfl_p
      jpmotop = 1
      do jl = 2, nlevels
        if ( zptopmbs >= xpres(jl - 1).and. &
             zptopmbs < xpres(jl)       ) then
          jpmotop = jl
          exit
        end if
      end do
      jpmolev = (nlevels - jpmotop + 1)

      allocate (iptobs_header (count_profile)      ,stat= alloc_status(1) )
      allocate (to_tl     (jpmolev,count_profile)  ,stat= alloc_status(2) )
      allocate (lqo_tl    (jpmolev,count_profile)  ,stat= alloc_status(3) )
      allocate (toext_tl  (nlevels  ,count_profile),stat= alloc_status(4) )
      allocate (qoext_tl  (nlevels  ,count_profile),stat= alloc_status(5) )
      allocate (zvlev     (nlv_T,count_profile)    ,stat= alloc_status(6) )
      allocate (dPdPs     (nlv_T,count_profile)    ,stat= alloc_status(7) )
      allocate (zt_tl     (nlv_T,count_profile)    ,stat= alloc_status(8) )
      allocate (zlq_tl    (nlv_T,count_profile)    ,stat= alloc_status(9) )
      allocate (zt        (nlv_T,count_profile)    ,stat= alloc_status(10))
      allocate (zlq       (nlv_T,count_profile)    ,stat= alloc_status(11))
      allocate (qoext     (nlevels,count_profile)  ,stat= alloc_status(12))
      allocate (zps_tl    (count_profile)          ,stat= alloc_status(13))
      call utl_checkAllocationStatus(alloc_status, "  tvslin_rttov_tl")
 
      iptobs_header(:) = 0 
      toext_tl (:,:) = 0.0d0
      zvlev    (:,:) = 0.0d0
      dPdPs    (:,:) = 0.0d0
      zt_tl    (:,:) = 0.0d0
      zlq_tl   (:,:) = 0.0d0
      zt       (:,:) = 0.0d0
      zlq      (:,:) = 0.0d0
      qoext    (:,:) = 0.0d0
      zps_tl   (:)   = 0.0d0
      to_tl    (:,:) = 0.0d0
      lqo_tl   (:,:) = 0.0d0

      ! allocate profiledata_tl structures
      init = .true.
      asw = 1 ! 1 to allocate
      call rttov_alloc_tl( &
              alloc_status(1),  &
              asw=1,              &
              nprofiles=count_profile,        &
              nchanprof=count_tb,        &
              nlevels=nlevels,          &
              chanprof=chanprof,         &
              opts=tvs_opts(sensor_id),             &
              profiles_tl=profilesdata_tl,      &
              coefs=tvs_coefs(sensor_id),            &
              transmission=transmission,     &
              transmission_tl=transmission_tl,  &
              radiance=radiancedata_d,         &
              radiance_tl=radiancedata_tl,      &
              calcemis=calcemis,         &
              emissivity=emissivity_local,       &
              emissivity_tl=emissivity_tl,    &
              init=.true.)

      call utl_checkAllocationStatus(alloc_status(1:1), " tovs_rtttov_tl rttov_alloc_tl 1")

      count_profile = 0

      obs_loop: do iobs = 1, NOBMAX
        if (tvs_lsensor(iobs) /= sensor_id) cycle obs_loop

        headerIndex = tvs_lobsno(iobs)
        count_profile = count_profile + 1
       
        zps_tl (count_profile) = col_getElem(column,1,headerIndex,'P0') * MPC_MBAR_PER_PA_R8
        delTT => col_getColumn(column,headerIndex,'TT')
        delLQ => col_getColumn(column,headerIndex,'HU')
        TTb => col_getColumn(columng,headerIndex,'TT')
        HUb => col_getColumn(columng,headerIndex,'HU')
        Pres => col_getColumn(columng,headerIndex,'PRES','TH')
        do jl = 1, nlv_T
          zt_tl (jl,count_profile) = delTT(jl)
          zlq_tl(jl,count_profile) = delLQ(jl)
          zt   (jl,count_profile)  = TTb(jl)
          zlq  (jl,count_profile)  = HUb(jl)
          zvlev(jl,count_profile)  = Pres(jl) *MPC_MBAR_PER_PA_R8
          dPdPs(jl,count_profile)  = col_getPressureDeriv(columng,jl,headerIndex,'TH')
        end do
        
        ! Fix pour eviter probleme au toit avec GEM 4
        ! (grosse varibilite temperature au dernier niveau thermo due 
        !  a l'extrapolation utilisee)
        if ( diagTtop ) then
          zt_tl   (1,count_profile) =  0.d0
          zlq_tl  (1,count_profile) =  0.d0
          zt   (1,count_profile) =  zt   (2,count_profile) + tvs_mesosphereLapseRate *  &
               log( col_getPressure(columng,1,headerIndex,'TH') /  &
               col_getPressure(columng,2,headerIndex,'TH') )
          zlq  (1,count_profile) =  zlq  (2,count_profile)
        end if
        
        iptobs_header(count_profile) = headerIndex
        
      end do obs_loop

      !     .  2.1  Vertical interpolation of model temperature and logarithm of
      !             specific humidity to pressure levels required by tovs rt model
      !     .       --------------------------------------------------------------

      do jn = 1, count_profile
        qoext(1:nlevels,jn) =  tvs_profiles(iptobs(jn)) % q(1:nlevels)
      end do
      
      to_tl (:,:) = 0.0d0
      lqo_tl(:,:) = 0.0d0
      
!$omp parallel do private(jn)
      do jn=1, count_profile

        call ppo_IntAvgTl(zvlev(:,jn:jn),dPdPs(:,jn:jn),zt_tl(:,jn:jn), &
             zt(:,jn:jn),zps_tl(jn:jn),nlv_T,nlv_T,1, &
             jpmolev,xpres(jpmotop:nlevels),to_tl(:,jn:jn))

        call ppo_IntAvgTl(zvlev(:,jn:jn),dPdPs(:,jn:jn),zlq_tl(:,jn:jn), &
             zlq(:,jn:jn),zps_tl(jn:jn),nlv_T,nlv_T,1, &
             jpmolev,xpres(jpmotop:nlevels),lqo_tl(:,jn:jn))

      end do
!$omp end parallel do
      
      !     .  2.2  Extrapolation of temperature profile above 10mb
      !     .       -----------------------------------------------
      toext_tl(:,:) = 0.0d0
      if ( .not. TopAt10hPa ) then
        do jn = 1, count_profile
          toext_tl(jpmotop:nlevels,jn) = to_tl(1:jpmolev,jn)
          toext_tl(1:jpmotop-1,jn) = 0.d0
        end do
      else
        call lextrap (to_tl,toext_tl,jpmolev,nlevels,count_profile)
      end if
      
      !     .  2.3  Extrapolation of humidity profile (kg/kg)
      !             above rlimlvhu (normally 300mbs or 70mbs)
      !     .       -----------------------------------------
      
      do jn = 1, count_profile
        do jk = 1, jpmotop - 1
          qoext_tl(jk,jn) = 0.d0
        end do
        do jk = 1, jpmolev
          qoext_tl(nlevels - jpmolev + jk,jn) = qoext(nlevels - jpmolev + jk,jn) * lqo_tl(jk,jn)
        end do
      end do
      
      if ( TopAt10hPa ) then
        if ( tvs_debug ) then
          do jn = 1, count_profile
            write(*,*)'qoext_tl*1000 avant exthum4    = '
            write(*,'(1x,10f8.4)')(qoext_tl(i,jn) * 1000.d0,i=1,nlevels)
            write(*,*)' '
          end do
        end if
        call lexthum4 (count_profile,nlevels,xpres(1:nlevels),qoext_tl,qoext)
        if ( tvs_debug ) then
          do jn = 1, count_profile
            write(*,*)'qoext_tl*1000 apres exthum4    = '
            write(*,'(1x,10f8.4)')(qoext_tl(i,jn) * 1000.d0,i=1,nlevels)
            write(*,*)' '
          end do
        end if
      end if

      do  j = 1 , count_profile
        profilesdata_tl(j) % gas_units       = gas_unit_specconc ! all gas profiles should be provided in kg/kg
        profilesdata_tl(j) % nlevels         =  nlevels
        profilesdata_tl(j) % nlayers         =  nlevels - 1
        if (tvs_coefs(sensor_id)%coef%nozone > 0) profilesdata_tl(j) % o3(:) =  0.d0
        
        profilesdata_tl(j) % ctp             = 0.0d0
        profilesdata_tl(j) % cfraction       = 0.0d0
        profilesdata_tl(j) % zenangle        = 0.0d0
        profilesdata_tl(j) % azangle         = 0.0d0
        profilesdata_tl(j) % skin % surftype = 0
        profilesdata_tl(j) % skin % t        = col_getElem(column,1,iptobs_header(j),'TG')
        profilesdata_tl(j) % skin % fastem(:)= 0.0d0
        profilesdata_tl(j) % skin % salinity = 0.0d0
        profilesdata_tl(j) % s2m % t         = col_getElem(column,ilowlvl_T,iptobs_header(j),'TT')

        
        profilesdata_tl(j) % s2m % q         = 0.d0

        profilesdata_tl(j) % s2m % p         = col_getElem(column,1,iptobs_header(j),'P0')*MPC_MBAR_PER_PA_R8
        profilesdata_tl(j) % s2m % u         = col_getElem(column,ilowlvl_M,iptobs_header(j),'UU')
        profilesdata_tl(j) % s2m % v         = col_getElem(column,ilowlvl_M,iptobs_header(j),'VV')
        
        profilesdata_tl(j) % p(1:nlevels)    = 0.d0
        profilesdata_tl(j) % t(1:nlevels)    = toext_tl(1:nlevels,j)
        profilesdata_tl(j) % q(1:nlevels)    = qoext_tl(1:nlevels,j)
      end do

      deallocate (iptobs_header ,stat= alloc_status(1) )
      deallocate (to_tl     ,stat= alloc_status(2) )
      deallocate (lqo_tl    ,stat= alloc_status(3) )
      deallocate (toext_tl  ,stat= alloc_status(4) )
      deallocate (qoext_tl  ,stat= alloc_status(5) )
      deallocate (zvlev     ,stat= alloc_status(6) )
      deallocate (dPdPs     ,stat= alloc_status(7) )
      deallocate (zt_tl     ,stat= alloc_status(8) )
      deallocate (zlq_tl    ,stat= alloc_status(9))
      deallocate (zt        ,stat= alloc_status(10))
      deallocate (zlq       ,stat= alloc_status(11))
      deallocate (qoext     ,stat= alloc_status(12))
      deallocate (zps_tl    ,stat= alloc_status(13))
      deallocate (xpres     ,stat= alloc_status(14))
      call utl_checkAllocationStatus(alloc_status, "tvslin_rttov_tl", .false.)

      !     set nthreads to actual number of threads which will be used.

      nthreads = min(max_nthreads, count_profile)  

      !     .  2.2  Prepare all input variables required by rttov.
      !     .       ---------------------------------------------------------
     
      allocate ( surfem1      (count_tb)       ,stat=alloc_status(1))
      allocate ( iptobs_body   (count_tb)      ,stat=alloc_status(2))
      call utl_checkAllocationStatus(alloc_status(1:2), " tovs_rtttov_tl")
    
      !     get Hyperspecral IR emissivities
      if ( tvs_isInstrumHyperSpectral(instrum) ) call tvs_getHIREmissivities(sensor_id, iptobs, &
           count_profile, lobsSpaceData, surfem1,assim_flag_val_opt=obs_ass_val)

      call TVS_getChanprof(sensor_id, iptobs, count_profile, lobsSpaceData, chanprof, &
           iptobs_cma_opt=iptobs_body,assim_flag_val_opt=obs_ass_val)

      call tvs_getOtherEmissivities(chanprof, iptobs, count_tb, sensor_type, instrum, surfem1, calcemis)
      
      emissivity_local(:)%emis_in = surfem1(:)

 
      !     .  2.3  Compute tl radiance with rttov_tl
      !     .       ---------------------------------
      
      errorstatus   = 0
      emissivity_tl(:)%emis_in = 0.0d0
      call tmg_start(87,'rttov_tl')
      call rttov_parallel_tl(         &
           errorstatus,               & ! out
           chanprof,                  & ! in
           tvs_opts(sensor_id),       & ! in
           tvs_profiles(iptobs(1:count_profile)),  & ! in
           profilesdata_tl,           & ! inout
           tvs_coefs(sensor_id),      & ! in
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

      !     .  2.4  Store hx in obsSpaceData,OBS_WORK
      !     .       ------------------------------------
      
      do index_tb = 1, count_tb
        
        bodyIndex = iptobs_body(index_tb)
        call obs_bodySet_r(lobsSpaceData,OBS_WORK,bodyIndex, &
             radiancedata_tl % bt(index_tb) )
        if ( tvs_debug ) then
          write(*,'(a,i4,2f8.2)') ' ichn,sim,obs= ', &
               chanprof(index_tb)%chan,   radiancedata_tl % bt(index_tb), &
               obs_bodyElem_r(lobsSpaceData,OBS_OMP,bodyIndex)
        end if
        
      end do

 
      !     de-allocate memory
      
      asw = 0 ! 0 to deallocate
      call rttov_alloc_tl( &
           alloc_status(1),              &
           asw=1,              &
           nprofiles=count_profile,        &
           nchanprof=count_tb,        &
           nlevels=nlevels,          &
           chanprof=chanprof,         &
           opts=tvs_opts(sensor_id),             &
           profiles_tl=profilesdata_tl,      &
           coefs=tvs_coefs(sensor_id),            &
           transmission=transmission,     &
           transmission_tl=transmission_tl,  &
           radiance=radiancedata_d,         &
           radiance_tl=radiancedata_tl,      &
           calcemis=calcemis,         &
           emissivity=emissivity_local,       &
           emissivity_tl=emissivity_tl )

      deallocate ( surfem1         ,stat=alloc_status(2) )
      deallocate ( iptobs_body     ,stat=alloc_status(3) )
      call utl_checkAllocationStatus(alloc_status(1:3), " tvslin_rtttov_tl", .false.)
      
    end do sensor_loop

    deallocate ( iptobs )


  end subroutine tvslin_rttov_tl


  subroutine tvslin_rttov_ad(column,columng,lobsSpaceData)
!--------------------------------------------------------------------------
!! *Purpose*: Adjoint of computation of radiance with rttov_ad
!!
!! @author j. halle *cmda/aes  april 19, 2005
!!
!revision 001  : a. beaulne *cmda/smc  june 2006
!                  -addition of ozone and IR surface emissivities
!revision 002  : r. sarrazin cmda  april 2008
!                  -adapt to CSR
!revision 003  : s. heilliette
!                  -adapt to IASI
!           S. Heilliette
!              - adaptation to rttov 10.0 (october 2010)
!revision 004  : s. macpherson  nov 2012
!                  - remove #include "comtovst.cdk"
!--------------------------------------------------------------------------
   

    implicit none

    type(struct_columnData) :: column,columng

    type(struct_vco), pointer :: vco_anl
    integer, allocatable :: iptobs    (:) 
    integer, allocatable :: iptobs_header (:) 
    
    integer :: alloc_status(14)
    logical :: diagTtop,TopAt10hPa
    integer :: omp_get_num_threads, nthreads
    integer :: nlevels,nobmax
    integer :: sensor_id, iobs, stat
    integer :: ilowlvl_T,ilowlvl_M,count_profile,headerIndex,nlv_M,nlv_T
    integer :: profile_index, level_index, rttov_index
    integer :: status, Vcode
    integer :: jpmotop, jpmolev
    
    real(8), allocatable :: to_ad    (:,:)
    real(8), allocatable :: lqo_ad   (:,:)
    real(8), allocatable :: toext_ad (:,:)
    real(8), allocatable :: qoext_ad (:,:)
    real(8), allocatable :: zvlev    (:,:)
    real(8), allocatable :: dPdPs    (:,:)
    real(8), allocatable :: zt_ad    (:,:)
    real(8), allocatable :: zlq_ad   (:,:)
    real(8), allocatable :: zt       (:,:)
    real(8), allocatable :: zlq      (:,:)
    real(8), allocatable :: qoext    (:,:)
    real(8), allocatable :: zps_ad   (:)
    real(8), allocatable :: xpres    (:)

    real(8) :: zptop, zptopmbs
   
    real(8), pointer :: uu_column(:),vv_column(:),tt_column(:),hu_column(:),ps_column(:),tg_column(:)
    real(8), pointer :: TTb(:), HUb(:), Pres(:)

    type(struct_obs) :: lobsSpaceData

    integer :: count_tb
    integer :: max_nthreads
    integer :: instrum
    integer :: tb_index, istart, iend, ichn, channel_index, bodyIndex
    integer :: sensor_type   ! sensor type (1=infrared; 2=microwave; 3=high resolution, 4=polarimetric)
    
    integer, allocatable :: iptobs_body      (:)
    integer :: errorstatus 
    
    real*8, allocatable :: surfem1       (:) 
    type (rttov_emissivity), pointer :: emissivity_local (:)
    type (rttov_emissivity), pointer :: emissivity_ad (:)
    type (rttov_transmission) :: transmission,transmission_ad
    type (rttov_radiance) :: radiancedata_ad, radiancedata_d
    
    type(rttov_profile)   , pointer :: profilesdata_ad (:) ! ad profiles buffer used in rttov calls
    type(rttov_chanprof)  , pointer :: chanprof(:)
    logical              :: init
    integer :: asw
    logical, pointer :: calcemis  (:)
         
    if (tvs_nobtov == 0) return      ! exit if there are not tovs data


    !     1.    Set index for model's lowest level and model top
    !     .     ------------------------------------------------

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
    !     .   ---------------------------------------------------------------------- 
 
!$omp parallel 
    max_nthreads = omp_get_num_threads()
!$omp end parallel

    allocate ( iptobs(tvs_nobtov) )

    !
    !     2.  Computation of adjoint hx for tovs data only
    !     .   --------------------------------------------

    ! Loop over all sensors specified by user

    sensor_loop:do  sensor_id = 1, tvs_nsensors
     
      nlevels = tvs_coefs(sensor_id) %coef % nlevels
      sensor_type = tvs_coefs(sensor_id) % coef% id_sensor
      instrum = tvs_coefs(sensor_id) % coef% id_inst

      count_profile = 0
      do iobs = 1, tvs_nobtov
        !    Currently processed sensor?
        if ( tvs_lsensor(iobs) == sensor_id ) then
          count_profile = count_profile + 1
          iptobs(count_profile) = iobs
        end if
      end do

      if (count_profile == 0) cycle sensor_loop

      NOBMAX = iptobs(count_profile)

      !     compute the number of radiances/tbs to be calculated
      count_tb = tvs_countRadiances(iptobs, count_profile, lobsSpaceData)

      if (count_tb == 0) cycle sensor_loop

     
      allocate (xpres(nlevels))
      xpres = tvs_coefs(sensor_id)% coef % ref_prfl_p
      jpmotop = 1
      do level_index = 2, nlevels
        if ( zptopmbs >= xpres(level_index - 1 ) .and.    &
             zptopmbs < xpres(level_index)        ) then
          jpmotop = level_index
          exit
        end if
      end do

      jpmolev = (nlevels - jpmotop + 1)
     
      alloc_status(:) = 0
      allocate (iptobs_header(count_profile)    ,stat= alloc_status(1) )
      allocate (to_ad    (jpmolev,count_profile),stat= alloc_status(2) )
      allocate (lqo_ad   (jpmolev,count_profile),stat= alloc_status(3) )
      allocate (toext_ad (nlevels,count_profile),stat= alloc_status(4) )
      allocate (qoext_ad (nlevels,count_profile),stat= alloc_status(5) )
      allocate (zvlev    (nlv_T,count_profile)  ,stat= alloc_status(6) )
      allocate (dPdPs    (nlv_T,count_profile)  ,stat= alloc_status(7) )
      allocate (zt_ad    (nlv_T,count_profile)  ,stat= alloc_status(8) )
      allocate (zlq_ad   (nlv_T,count_profile)  ,stat= alloc_status(9))
      allocate (zt       (nlv_T,count_profile)  ,stat= alloc_status(10))
      allocate (zlq      (nlv_T,count_profile)  ,stat= alloc_status(11))
      allocate (qoext    (nlevels,count_profile),stat= alloc_status(12))
      allocate (zps_ad   (count_profile)        ,stat= alloc_status(13))

      call utl_checkAllocationStatus(alloc_status, " tvslin_fill_profiles_ad")
      !  loop over all obs.
      count_profile = 0 
      
      ! loop over all obs.
      obs_loop: do iobs = 1, NOBMAX
        if (tvs_lsensor(iobs)/=sensor_id) cycle obs_loop
        headerIndex = tvs_lobsno(iobs)
        count_profile = count_profile + 1
        
        TTb => col_getColumn(columng,headerIndex,'TT')
        HUb => col_getColumn(columng,headerIndex,'HU')
        Pres => col_getColumn(columng,headerIndex,'PRES','TH')
        do level_index = 1, nlv_T
          zt   (level_index,count_profile) = TTb(level_index)
          zlq  (level_index,count_profile) = HUb(level_index)
          zvlev(level_index,count_profile) = Pres(level_index) * MPC_MBAR_PER_PA_R8
          dPdPs(level_index,count_profile) = col_getPressureDeriv(columng,level_index,headerIndex,'TH')
        end do
        
        ! Fix pour eviter probleme au toit avec GEM 4
        ! (grosse variabilite de la temperature au dernier niveau thermo due 
        !  a l'extrapolation utilisee)
        if (diagTtop) then
          zt   (1,count_profile) =  zt   (2,count_profile) + tvs_mesosphereLapseRate *  &
               log( col_getPressure(columng,1,headerIndex,'TH') /  &
               col_getPressure(columng,2,headerIndex,'TH') )
          zlq  (1,count_profile) =  zlq  (2,count_profile)
        end if

        iptobs_header(count_profile) = headerIndex
      end do obs_loop
     
      !     .  2.1  Calculate the actual number of threads which will be used.
      !     .       ----------------------------------------------------------

      nthreads = min(max_nthreads, count_profile )  

      !     .  2.2  Prepare all input variables required by rttov_ad.
      !     .       ---------------------------------------------------------

      asw = 1
      init = .true.

      call rttov_alloc_ad(                  &
           alloc_status(1),                 &
           asw,                             &
           count_profile,                   &
           count_tb,                        &
           nlevels,                         &
           chanprof,                        &
           opts=tvs_opts(sensor_id),        &
           profiles_ad=profilesdata_ad,     &
           coefs=tvs_coefs(sensor_id),       &
           transmission= transmission,      &
           transmission_ad= transmission_ad,&
           radiance=radiancedata_d,         &
           radiance_ad=radiancedata_ad,     &
           calcemis=calcemis,               &
           emissivity=emissivity_local,     &
           emissivity_ad=emissivity_ad,     &
           init=init)

      allocate ( surfem1      (count_tb)        ,stat=alloc_status(2))
      allocate ( iptobs_body   (count_tb)       ,stat=alloc_status(3))
      call utl_checkAllocationStatus(alloc_status(1:3), " tvslin_rttov_ad")
      
!     get Hyperspectral IR emissivities
      
      if ( tvs_isInstrumHyperSpectral(instrum) ) call tvs_getHIREmissivities(sensor_id, iptobs, count_profile, lobsSpaceData, surfem1)

      ! Build the list of channels/profiles indices

      call TVS_getChanprof(sensor_id, iptobs, count_profile, lobsSpaceData, chanprof, iptobs_body)
      !     get non Hyperspectral IR emissivities
      call tvs_getOtherEmissivities(chanprof, iptobs, count_tb, sensor_type, instrum, surfem1, calcemis)

      emissivity_local(:)%emis_in = surfem1(:)
        
      ! allocate transmittance structures
      
!      call rttov_init_rad( radiancedata_d )

      !.. fill radiancedata_ad arrays
!      call rttov_init_rad( radiancedata_ad ) ! set irrrelevant fields to zero

      do tb_index = 1, count_tb
        bodyIndex = iptobs_body(tb_index)
        radiancedata_ad % bt( tb_index ) = obs_bodyElem_r(lobsSpaceData,OBS_WORK,bodyIndex)
      end do


      !     .  2.3  Compute ad radiance with rttov_ad
      !     .       ---------------------------------

      errorstatus  = 0
      emissivity_ad(:) % emis_in = 0.0d0
      emissivity_ad(:) % emis_out = 0.0d0
     
      call tmg_start(84,'rttov_ad')
      call rttov_parallel_ad(        &
           errorstatus,              & ! out
           chanprof,                 & ! in
           tvs_opts(sensor_id),      & ! in
           tvs_profiles(iptobs(1:count_profile)), & ! in
           profilesdata_ad,          & ! in
           tvs_coefs(sensor_id),     & ! in
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
      toext_ad = 0.d0
      qoext_ad = 0.d0

      do tb_index = 1, count_tb
        
        profile_index = chanprof(tb_index)%prof
        headerIndex = iptobs_header(profile_index)

        ps_column => col_getColumn(column,headerIndex,'P0')
        tg_column => col_getColumn(column,headerIndex,'TG')
        tt_column => col_getColumn(column,headerIndex,'TT')
        hu_column => col_getColumn(column,headerIndex,'HU')
        uu_column => col_getColumn(column,headerIndex,'UU')
        vv_column => col_getColumn(column,headerIndex,'VV')

        toext_ad(:,profile_index)   =  profilesdata_ad(profile_index) % t(:)
        qoext_ad(:,profile_index)   =  profilesdata_ad(profile_index) % q(:)
        tg_column(1)                =  profilesdata_ad(profile_index) % skin % t 
        tt_column(ilowlvl_T)        =  profilesdata_ad(profile_index) % s2m % t
        ps_column(1)                =  profilesdata_ad(profile_index) % s2m % p * MPC_MBAR_PER_PA_R8
        hu_column(ilowlvl_T)        =  0.d0 
        uu_column(ilowlvl_M)        =  profilesdata_ad(profile_index) % s2m % u
        vv_column(ilowlvl_M)        =  profilesdata_ad(profile_index) % s2m % v
      end do


      !     .  2.4  Adjoint of filling profiles_ad structure
      !     .       ----------------------------------------
      do profile_index = 1, count_profile
        qoext(:,profile_index) =  tvs_profiles(iptobs(profile_index)) % q(:)
      end do

    
      !     .  2.3  Adjoint of extrapolation of humidity profile (kg/kg)
      !             above rlimlvhu (normally 300mbs or 70mbs)
      !     .  
      if ( TopAt10hPa ) then
        if ( tvs_debug ) then
          do profile_index = 1, count_profile
            write(*,*)'qoext_ad*1000 avant aexthum4    = '
            write(*,'(1x,10f8.4)')(qoext_ad(level_index,profile_index)*1000.d0,level_index=1,nlevels)
            write(*,*)' '
          end do
        end if
        call aexthum4 (count_profile,nlevels,xpres(1:nlevels),qoext_ad,qoext)
        if ( tvs_debug ) then
          do profile_index = 1, count_profile
            write(*,*)'qoext_ad*1000 apres aexthum4    = '
            write(*,'(1x,10f8.4)')(qoext_ad(level_index,profile_index)*1000.d0,level_index=1,nlevels)
            write(*,*)' '
          end do
        end if
      end if

      ! adjoint of conversion lnq --> q
      lqo_ad(:,:) = 0.0d0
      do profile_index = 1, count_profile
        do level_index = 1, jpmolev
          lqo_ad(level_index,profile_index) = qoext_ad(nlevels-jpmolev + &
               level_index,profile_index) * qoext(nlevels-jpmolev+level_index,profile_index)
        end do
      end do

      !     .  2.2  Adjoint of extrapolation of temperature profile above 10mb
      !     .       ----------------------------------------------------------
      to_ad(:,:) = 0.0d0
      if ( .not. TopAt10hPa ) then
        do profile_index = 1, count_profile
          to_ad(1:jpmolev,profile_index) = to_ad(1:jpmolev,profile_index) + &
               toext_ad(jpmotop:nlevels,profile_index)
        end do
      else
        call aextrap (to_ad,toext_ad,jpmolev,nlevels,count_profile)
      end if
 
      !     .  2.1  Adjoint of vertical interpolation of model temperature and logarithm of
      !             specific humidity to pressure levels required by tovs rt model
      !     .       -----------------------------------------------------------------------
      
      zt_ad (:,:) = 0.0d0
      zlq_ad(:,:) = 0.0d0
      zps_ad(:)   = 0.0d0

      call tmg_start(75,'intavgad')
!$omp parallel do private(profile_index)
      do profile_index = 1, count_profile
        
        call ppo_IntAvgAd(zvlev(:,profile_index:profile_index), dPdPs(:,profile_index:profile_index), &
             zt_ad(:,profile_index:profile_index), zt(:,profile_index:profile_index), &
             zps_ad(profile_index:profile_index), nlv_T,nlv_T,1, &
             jpmolev,xpres(jpmotop:nlevels),to_ad(:,profile_index:profile_index))
        
        call ppo_IntAvgAd(zvlev(:,profile_index:profile_index),dPdPs(:,profile_index:profile_index), &
             zlq_ad(:,profile_index:profile_index), zlq(:,profile_index:profile_index), &
             zps_ad(profile_index:profile_index), nlv_T,nlv_T,1, &
             jpmolev,xpres(jpmotop:nlevels), lqo_ad(:,profile_index:profile_index))
       
      end do
!$omp end parallel do
      call tmg_stop(75)
 
      ! Fix pour eviter probleme au toit avec GEM 4
      ! (grosse variabilite temperature au dernier niveau thermo due 
      !  a l'extrapolation utilisee)
      if ( diagTtop ) then
        do profile_index = 1, count_profile
          zt_ad (1,profile_index) = 0.d0
          zlq_ad(1,profile_index) = 0.d0
        end do
      end if

      !     .  2.1  Store adjoints in columnData object
      !     .       -----------------------------------

      do  profile_index = 1 , count_profile 
        ps_column => col_getColumn(column,iptobs_header(profile_index),'P0')
        tt_column => col_getColumn(column,iptobs_header(profile_index),'TT')
        hu_column => col_getColumn(column,iptobs_header(profile_index),'HU')
        
        ps_column(1) = ps_column(1) + zps_ad  (profile_index) * MPC_MBAR_PER_PA_R8
        do level_index = 1, col_getNumLev(column,'TH')
          tt_column(level_index) = tt_column(level_index) + zt_ad  (level_index,profile_index)
          hu_column(level_index) = hu_column(level_index) + zlq_ad (level_index,profile_index)
        end do
      end do

      deallocate (iptobs_header,stat= alloc_status(1) )
      deallocate (to_ad    ,stat= alloc_status(2) )
      deallocate (lqo_ad   ,stat= alloc_status(3) )
      deallocate (toext_ad ,stat= alloc_status(4) )
      deallocate (qoext_ad ,stat= alloc_status(5) )
      deallocate (zvlev    ,stat= alloc_status(6) )
      deallocate (dPdPs    ,stat= alloc_status(7) )
      deallocate (zt_ad    ,stat= alloc_status(8) )
      deallocate (zlq_ad   ,stat= alloc_status(9) )
      deallocate (zt       ,stat= alloc_status(10))
      deallocate (zlq      ,stat= alloc_status(11))
      deallocate (qoext    ,stat= alloc_status(12))
      deallocate (zps_ad   ,stat= alloc_status(13))
      deallocate (xpres    ,stat= alloc_status(14))
      
      call utl_checkAllocationStatus(alloc_status, " tvslin_fill_profiles_ad", .false.)
    
      !     de-allocate memory

      asw = 0 ! 0 to deallocate
      call rttov_alloc_ad(                  &
           alloc_status(1),                 &
           asw,                             &
           count_profile,                   &
           count_tb,                        &
           nlevels,                         &
           chanprof,                        &
           opts=tvs_opts(sensor_id),        &
           profiles_ad=profilesdata_ad,     &
           coefs=tvs_coefs(sensor_id),       &
           transmission= transmission,      &
           transmission_ad= transmission_ad,&
           radiance=radiancedata_d,         &
           radiance_ad=radiancedata_ad,     &
           calcemis=calcemis,               &
           emissivity=emissivity_local,     &
           emissivity_ad=emissivity_ad )
     

      
      deallocate ( surfem1         ,stat=alloc_status(2))
      deallocate ( iptobs_body     ,stat=alloc_status(3))
      call utl_checkAllocationStatus(alloc_status(1:3), " tvslin_rttov_ad", .false.)
     
    end do sensor_loop

    !     3.  Close up
    !     .   --------

    deallocate ( iptobs )

  end subroutine tvslin_rttov_ad


  subroutine AEXTHUM4(KNPF,KLAPF,PPRES,PAV,PAV5)
!
!*****aexthum4* - adjoint of extrapolation of upper level humidity profile.
!                (adapted from exthumad by J. Eyre)
!
!     purpose.
!     --------
!          ad of routine
!          to extend mixing ratio profile into stratosphere in
!          a reasonable way.
!
!**   interface.
!     ----------
!          *call* *aexthum4(knpf,klapf,ppres,pav,pav5)*
!               *knpf*:  no. of profiles to be processed.
!               *klapf*: length of atm. profiles.
!               *ppres*: pressure levels of atm. profiles.
!               *pav*:   ad of humidity profiles.
!               *pav5*:  humidity profiles.
!
!     method.
!     -------
!          take top tropospheric mixing ratio (e.g. near 300 mb) and
!          extrapolate with given fall off into lower stratosphere
!          (e.g. to 70 mb).  constrain mixing ratio to be >= zwmin
!          (e.g. 0.000003 kg/kg).   in upper strat, mixing ratio = zwmin.
!
!     externals.
!     ----------
!          none.
!

   
    implicit none

    integer :: klapf,knpf
    REAL(8) PPRES(*),PAV(KLAPF,*), PAV5(KLAPF,*)

    REAL(8) :: ZPRES3(KLAPF)

    real(8) zwmix,zwb
    real(8),parameter ::  ZP1 = 70.0D0 ! PRESS LIMITS (IN HPA) OF REGION to be extrapolated
    integer :: inlvw,j,jnpf
    

    !          find top level of given profile
    INLVW = -1
    do J=KLAPF,1,-1
      if (PPRES(J) < FILT_RLIMLVHU) then
        INLVW = J
        exit
      end if
    end do
!
    !** Null extrapolation case
!
    if (INLVW == -1) return

    !          constants defining p**3 fall off around tropopause
    do J=1,INLVW
      ZPRES3(J) = (PPRES(J)/PPRES(INLVW+1))**3
    end do

    do JNPF=1,KNPF
      ZWB = 0.D0
      do J=1,INLVW
        ZWMIX = PAV(J,JNPF)
        PAV(J,JNPF) = 0.D0
        if (PPRES(J) >= ZP1) then
          if (PAV5(J,JNPF) > MPC_MINIMUM_HU_R8) then
            ZWB = ZWB + ZWMIX * ZPRES3(J)
          end if
        end if
      end do
      PAV(INLVW + 1,JNPF) = PAV(INLVW + 1,JNPF) + ZWB
    end do
   

  end subroutine AEXTHUM4

 
End module tovs_lin_mod


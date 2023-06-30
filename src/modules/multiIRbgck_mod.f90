
module multiIRbgck_mod
  ! MODULE multiIRbgck_mod (prefix='irbg' category='1. High-level functionality')
  !
  !:Purpose:  Variables for multispectral infrared background check and quality
  !           control.
  !
  use rttovInterfaces_mod
  use tovsNL_mod
  use rttov_const, only : inst_id_iasi
  use rttov_types, only :   &
       rttov_coefs         ,&
       rttov_profile       ,&
       rttov_radiance      ,&
       rttov_transmission  ,&
       rttov_chanprof      ,&
       rttov_emissivity
  use utilities_mod
  use obsSpaceData_mod
  use midasMpi_mod
  use columnData_mod
  use MathPhysConstants_mod
  implicit none
  save
  private
  ! Public functions (methods)
  public :: irbg_setup, irbg_bgCheckIR

  integer, parameter :: nClassAVHRR = 7
  integer, parameter :: NIR = 3, NVIS = 3
  integer, parameter :: nChanAVHRR = NIR + NVIS
  integer, parameter :: nmaxinst = 10
 
  ! Cloud top units : (1) mb, (2) meters
  ! (subroutines cloud_height (iopt1) and cloud_top (iopt2))

  integer, parameter        :: iopt1 = 2   ! verify subr input if iopt1 changes
  integer, parameter        :: iopt2 = 1

  ! Cloud top based on which background profile matching (subroutine cloud_top)
  ! (0) brightness temperature, (1) radiance, (2) both

  integer, parameter        :: ihgt = 2

  ! Highest flag in post files (value of N in 2^N)

  integer, parameter :: bitflag = 29

  real(8),parameter ::  albedoThresholdQC= 0.25d0 ! max albedo allowed over the water 

  real(8),parameter :: seuilalb_static(NIR,0:2)= reshape( (/ 70.0,67.0,50.0, &
                                                             40.0,37.0,37.0, &
                                                             70.0,57.0,40. /),(/3,3/) ) 
  real(8),parameter :: seuilalb_homog(NIR,0:2)= reshape( (/ 15.0,18.0,13.0, &
                                                            9.0,10.0,10.0, &
                                                            18.0,16.0,10.0 /),(/3,3/) )
  real(8) :: seuilbt_homog(NVIS+1:nChanAVHRR,0:2,1:2)= reshape( (/5.d0, 4.d0, 4.d0, 4.d0, 3.d0, 3.d0, &
                                                                  5.d0, 4.d0, 4.d0, 5.d0, 5.d0, 5.d0, &
                                                                  4.d0, 3.d0, 3.d0, 5.d0, 5.d0, 5.d0/), (/3,3,2/) )

  ! Number of channels to use for cloud top height detection
  ! with the "background profile matching" method (subroutine cloud_top)
  integer, parameter        :: nch_he = 4

  ! Number of channels pairs to use for cloud top height detection
  ! with the CO2-slicing method. (subroutine co2_slicing)
  integer, parameter  :: nco2 = 13

  ! Namelist variables:
  integer :: ninst                      ! MUST NOT BE INCLUDED IN NAMELIST!
  character(len=7) :: inst(nmaxinst)    ! List of instrument names
  integer :: iwindow(nmaxinst)          ! Ref window channel for clear/cloudy profile detection
  integer :: iwindow_alt(nmaxinst)      ! Alternate window channel for clear/cloudy profile detection
  integer :: ilist1(nmaxinst,nch_he)    ! Chan numbers for cloud top height detection: background profile matching 
  integer :: ilist2(nmaxinst,nco2)      ! Chan numbers for cloud top height detection: CO2-slicing
  integer :: ilist2_pair(nmaxinst,nco2) ! Chan number pairs for cloud top height detection: CO2-slicing
  real(8) :: dtw                        ! Max delta allowed btwn guess and true skin temp over water
  real(8) :: dtl                        ! Max delta allowed btwn guess and true skin temp over land
  real(8) :: pco2min                    ! Min RTTOV level for lev_start variable entering CO2 slicing in mb
  real(8) :: pco2max                    ! Max RTTOV level for lev_start variable entering CO2 slicing in mb
  real(8) :: night_ang                  ! Min solar zenith angle for night (between 90 and 180)
  real(8) :: crisCloudFractionThreshold ! threshold for CrIS cloud detection from VIIRS cloud mask

  type(rttov_coefs) :: coefs_avhrr

  type avhrr_bgck_iasi
     real(8) :: radmoy(nClassAVHRR,nChanAVHRR)
     real(8) :: radstd(nClassAVHRR,nChanAVHRR)
     real(8) :: cfrac(nClassAVHRR)
     real(8) :: tbmoy(nClassAVHRR,NVIS+1:NVIS+NIR)
     real(8) :: tbstd(nClassAVHRR,NVIS+1:NVIS+NIR)
     real(8) :: albedmoy(nClassAVHRR,1:NVIS)
     real(8) :: albedstd(nClassAVHRR,1:NVIS)
     real(8) :: tbstd_pixelIASI(NVIS+1:NVIS+NIR)
     real(8) :: albstd_pixeliasi(1:NVIS)
     real(8) :: radclearcalc(NVIS+1:NVIS+NIR)
     real(8) :: tbclearcalc(NVIS+1:NVIS+NIR)
     real(8),allocatable :: radovcalc(:,:)
     real(8) :: transmsurf(NVIS+1:NVIS+NIR)
     real(8) :: emiss(NVIS+1:NVIS+NIR)
  end type avhrr_bgck_iasi

  type(avhrr_bgck_iasi), allocatable :: avhrr_bgck(:) ! avhrr parameters for IASI quality control

contains


  !--------------------------------------------------------------------------
  ! irbg_init
  !--------------------------------------------------------------------------
  subroutine  irbg_init()
    !
    !:Purpose: This subroutine reads the namelist section NAMBGCKIR
    !          for the module.
    implicit none

    ! Locals:
    integer :: nulnam, ierr
    logical, save :: first = .true.
    integer, external :: fnom, fclos
    integer :: instrumentIndex
    namelist /NAMBGCKIR/ ninst, inst, iwindow, iwindow_alt, ilist1, ilist2, ilist2_pair
    namelist /NAMBGCKIR/ dtw, dtl, pco2min, pco2max, night_ang, crisCloudFractionThreshold

    if (first) then

      ! set the default values for namelist variables
      ninst = MPC_missingValue_INT
      inst(:) = '*EMPTY*'
      iwindow(:) = 0
      iwindow_alt(:) = 0
      ilist1(:,:) = 0
      ilist2(:,:) = 0
      ilist2_pair(:,:) = 0
      dtw = 0.0d0
      dtl = 0.0d0
      pco2min = 0.0d0
      pco2max = 0.0d0
      night_ang = 0.0d0
      crisCloudFractionThreshold = -1.d0

      ! read the namelist
      nulnam = 0
      ierr = fnom(nulnam, './flnml','FTN+SEQ+R/O', 0)
      read(nulnam, nml=NAMBGCKIR, iostat=ierr)
      if (ierr /= 0) call utl_abort('irbg_init: Error reading namelist')
      if (mmpi_myid == 0) write(*, nml=NAMBGCKIR)
      if (ninst /= MPC_missingValue_INT) then
        call utl_abort('irbg_init: check NAMBGCKIR namelist section: ninst should be removed')
      end if
      ninst = 0
      do instrumentIndex = 1, nmaxinst
        if (inst(instrumentIndex) == '*EMPTY*') exit
        ninst = ninst + 1
      end do
      ierr = fclos(nulnam)
      first = .false.

    end if

  end subroutine irbg_init

  !--------------------------------------------------------------------------
  ! irbg_setup
  !--------------------------------------------------------------------------
  subroutine irbg_setup()
    !
    ! :Purpose: Memory allocation for the Hyperspectral Infrared
    !           background check variables.
    !
    implicit none

    ! Locals:
    integer :: allocStatus(2)
    integer :: tovsIndex, maxChannelNumber
    integer :: sensorIndex, channelNumber

    !     Memory allocation for background check related variables
    allocStatus(:) = 0
    allocate( tvs_surfaceParameters(tvs_nobtov), stat=allocStatus(1))
    call utl_checkAllocationStatus(allocStatus(1:1), " irbg_setup tvs_surfaceParameters")

    !___ emissivity by profile

    maxChannelNumber = 1
    do tovsIndex = 1, tvs_nobtov
      sensorIndex = tvs_lsensor(tovsIndex)
      channelNumber = tvs_nchan(sensorIndex)
      if (channelNumber > maxChannelNumber) maxChannelNumber=channelNumber
    end do

    allocate( tvs_emissivity (maxChannelNumber,tvs_nobtov), stat=allocStatus(1))
    call utl_checkAllocationStatus(allocStatus(1:1), " irbg_setup tvs_emissivity")
    
    do sensorIndex = 1, tvs_nsensors
      if ( tvs_instruments(sensorIndex) == inst_id_iasi ) then
        allocate ( avhrr_bgck(tvs_nobtov), stat=allocStatus(1))
        call utl_checkAllocationStatus(allocStatus(1:1), " irbg_setup avhrr_bgck")
        exit
      end if
    end do

  end subroutine irbg_setup

  !--------------------------------------------------------------------------
  ! irbg_bgCheckIR
  !--------------------------------------------------------------------------
  subroutine irbg_bgCheckIR(columnTrlOnTrlLev, obsSpaceData)
    !
    ! :Purpose: Do background check on all hyperspectral infrared observations
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs),        intent(inout) :: obsSpaceData
   
    ! Locals:
    integer,allocatable :: nobir(:)
    integer             :: headerIndex, idatyp, sensorIndex, instrumentIndex
    logical             :: irDataPresent

    call utl_tmg_start(115,'--BgckInfrared')

    irDataPresent = .false.
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER0
      idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if ( tvs_isIdBurpHyperSpectral(idatyp) ) then
        irDataPresent = .true.
      end if
    end do HEADER0

    if ( .not. irDataPresent ) then
      write(*,*) 'WARNING: WILL NOT RUN irbg_bgCheckIR since no IR'
      return
    end if

    write(*,'(A)') " ****************"
    write(*,'(A)') " BEGIN IR BACKGROUND CHECK"
    write(*,'(A)') " **************** **************** ****************"

    !  Preliminary initializations

    tvs_nobtov = 0
    call irbg_init()
    allocate (nobir(ninst))
    nobir = 0
  

    ! Loop over all header indices of the 'TO' family
    ! Set the header list (and start at the beginning of the list)
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)

      if ( .not. tvs_isIdBurpTovs(idatyp) ) then
        write(*,*) 'irbg_bgCheckIR: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER   ! Proceed to the next header_index
      end if
      tvs_nobtov = tvs_nobtov + 1
      do instrumentIndex=1, ninst
        if ( tvs_isIdBurpInst(idatyp,inst(instrumentIndex)) ) then
          nobir(instrumentIndex) = nobir(instrumentIndex) + 1
          exit
        end if
      end do
    end do HEADER
    
    do instrumentIndex=1, ninst
      if (nobir(instrumentIndex) > 0) then
        do sensorIndex=1,tvs_nsensors
          if (tvs_instruments(sensorIndex) ==  tvs_getInstrumentId(inst(instrumentIndex)) ) then
            call irbg_doQualityControl (columnTrlOnTrlLev, obsSpaceData, inst(instrumentIndex), sensorIndex)
          end if
        end do
      end if
    end do
    deallocate (nobir)

    call utl_tmg_stop(115)

  end subroutine irbg_bgCheckIR
 
  !--------------------------------------------------------------------------
  !  bgck_get_qcid
  !--------------------------------------------------------------------------
  subroutine bgck_get_qcid(instrumentName, qcid )
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in)  :: instrumentName
    integer,          intent(out) :: qcid

    ! Locals:
    integer :: i 

    qcid = -1

    do i=1, ninst
      if (trim(instrumentName) == trim(inst(i))) then
        qcid = i
        exit
      end if
    end do

    if (qcid == -1) then
      write(*,*) "Unknown instrument ", instrumentName
      call utl_abort('bgck_get_qcid')
    end if

  end subroutine bgck_get_qcid

  !--------------------------------------------------------------------------
  ! irbg_doQualityControl
  !--------------------------------------------------------------------------
  subroutine irbg_doQualityControl (columnTrlOnTrlLev, obsSpaceData, instrumentName, id_opt)
    !
    ! :Purpose: Quality control of hyperspectral infrared observations.
    !           assign assimilation flags to observations 
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs),        intent(inout) :: obsSpaceData
    character(len=*),        intent(in)    :: instrumentName
    integer,       optional, intent(in)    :: id_opt

    ! Locals:
    integer       :: jc, nchn, levelIndex, bitIndex, channelNumber, classIndex
    integer       :: columnIndex
    integer       :: iwindo, iwindo_alt
    integer       :: bodyIndex, bodyStart, bodyEnd, headerIndex
    integer       :: idatyp
    real(8)       :: difftop_min
    integer       :: modelTopIndex
    integer       :: count
    real(8)       :: t_effective
    integer       :: allocStatus(27)
    real(8)       :: tg, p0, tskinRetrieved, ptop_T, qs
    real(8), allocatable :: tt(:), height(:,:)
    real(8), allocatable :: pressure(:,:)
    real(8), allocatable :: btObsErr(:), btObs(:), btCalc(:), rcal_clr(:), sfctau(:)
    real(8), allocatable :: radObs(:), cloudyRadiance(:,:), transm(:,:), emi_sfc(:) 
    real(8), allocatable :: ptop_bt(:), ptop_rd(:)
    real(8), allocatable :: pmin(:), dtaudp1(:), maxwf(:)
    real(8), allocatable :: cloudyRadiance_avhrr(:,:)
    integer, allocatable :: rejflag(:,:) 
    integer, allocatable :: ntop_bt(:), ntop_rd(:)
    integer, allocatable :: minp(:), fate(:), channelIndexes(:)
    real(8) :: clfr, sunZenithAngle, satelliteAzimuthAngle, satelliteZenithAngle, sunAzimuthAngle
    real(8) :: albedo, ice, pcnt_wat, pcnt_reg
    real(8) :: ptop_eq,ptop_mb
    real(8) :: ptop_co2(nco2),fcloud_co2(nco2)
    real(8) :: etop,vtop,ecf,vcf,heff
    real(8) :: tampon,cfsub
    real(8) :: tskinRetrieved_avhrr(nClassAVHRR), sfctau_avhrr(NIR), emi_sfc_avhrr(NIR), rcal_clr_avhrr(NIR)
    real(8) :: ptop_bt_avhrr(NIR,nClassAVHRR), ptop_rd_avhrr(NIR,nClassAVHRR)
    real(8) :: btObs_avhrr(NIR,nClassAVHRR), radObs_avhrr(NIR,nClassAVHRR), ptop_eq_avhrr(nClassAVHRR)
    real(8) :: cfrac_avhrr
    real(8) :: avhrr_surfem1(NIR)
    real(8) :: albedoThreshold(NIR)
    integer :: ksurf,ltype
    integer :: cldflag, lev_start   
    integer :: gncldflag
    integer :: ichref
    integer :: ntop_eq, ntop_mb
    integer :: ngood
    integer :: ntop_co2(nco2)
    integer :: cldflag_avhrr(nClassAVHRR),lev_start_avhrr(nClassAVHRR),ichref_avhrr(nClassAVHRR),ntop_rd_avhrr(NIR,nClassAVHRR)
    integer :: ntop_bt_avhrr(NIR,nClassAVHRR),ntop_eq_avhrr(nClassAVHRR)
    logical :: assim_all
    logical :: sunZenithAnglePresent
    logical :: satelliteAzimuthAnglePresent, satelliteZenithAnglePresent, sunAzimuthAnglePresent
    integer, parameter :: nn=2
    integer, parameter :: ilist_avhrr(nn)=(/ 2 ,3 /)
    integer :: countInvalidChannels
    logical :: bad
    real(8), parameter :: sunzenmax=87.12d0
    real(8) :: minpavhrr(2:3)
    real(8) :: anisot,zlamb,zcloud,scos,del,deltaphi
    integer :: ier,ijour,iloc(2:3),co2min(1),co2max(1)
    integer :: channelIndex,ilist_co2(nco2),ilist_co2_pair(nco2),ilist_he(nch_he)
    integer :: nlv_T,id,sensorIndex, qcid, nchannels
    logical :: liasi,lairs,lcris
    type (rttov_profile), pointer :: profiles(:)

    write (*,*) "Entering irbg_doQualityControl"

    call tvs_getProfile(profiles, 'nl')

    liasi= ( trim(instrumentName) == "IASI" .or.  trim(instrumentName) == "iasi")
    lairs= ( trim(instrumentName) == "AIRS" .or.  trim(instrumentName) == "airs")
    lcris= ( trim(instrumentName) == "CRIS" .or.  trim(instrumentName) == "cris" .or. trim(instrumentName) == "CRISFSR" .or. &
             trim(instrumentName) == "crisfsr" )

    call bgck_get_qcid(instrumentName,qcid)
    
    if (present(id_opt)) then
      id = id_opt
    else
      !  Find sensor number corresponding to the desired instrument
      id = -1
      do sensorIndex = 1, tvs_nsensors
        if ( trim(tvs_instrumentName(sensorIndex)) == TRIM(instrumentName)) then
          id = sensorIndex
          exit
        end if
      end do
      if (id < 0) call utl_abort("irbg_doQualityControl: should not happen !")
    end if

    !  Find number of profiles 
    count = 0

    ! Loop over all header indices of the 'TO' family
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      if ( tvs_tovsIndex( headerIndex ) < 0) cycle HEADER
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if ( tvs_isIdBurpInst(idatyp,instrumentName) .and. tvs_lsensor(tvs_tovsIndex (headerIndex)) == id ) then
        count = count + 1
      end if
    end do HEADER

    if ( count == 0 ) return
    ! Find number of channels and RTTOV levels

    nchn = tvs_coefs(id) % coef % fmv_chn

    nlv_T = col_getNumLev(columnTrlOnTrlLev,'TH')

    write(*,*) ' irbg_doQualityControl - nchn ', nchn
   

    ! information to extract (transvidage)
    !
    ! tg -- guess skin temperatures (deg K)
    ! p0 -- surface pressure (hPa)
    ! tt(nlv_T-1) -- temperature profiles on NWP model levels (deg K)
    ! height(nlv_T-1,1) -- height profiles on NWP model levels (m)
    ! qs -- surface specific humidity (kg/kg)
    ! btObsErr(nchn) -- observation error standard deviation
    ! btObs(nchn) -- observed brightness temperatures (deg K)
    ! btCalc(nchn) -- computed brightness temperatures (deg K)
    ! rcal_clr(nchn) -- computed clear radiances (mw/m2/sr/cm-1)
    ! sfctau(nchn) -- surface to space transmittances (0-1)
    ! cloudyRadiance(nchn,nlv_T-1) -- overcast cloudy radiances (mw/m2/sr/cm-1)
    ! transm(nchn,nlv_T) -- layer to space transmittances (0-1)
    ! emi_sfc(nchn) -- surface emissivities (0-1)
    ! ksurf -- surface type in obs file (0, 1)
    ! clfr -- cloud fraction (%)
    ! sunZenithAngle -- sun zenith angle (deg)
    ! satelliteAzimuthAngle -- satellite azimuth angle (deg)
    ! satelliteZenithAngle -- satellite zenith angle (deg)
    ! albedo -- surface albedo (0-1)
    ! ice -- ice fraction (0-1)
    ! ltype -- surface type (1,...,20)
    ! pcnt_wat -- water fraction (0-1)
    ! pcnt_reg -- water fraction in the area (0-1)
    ! radObs(nchn) -- observed radiances (mW/m2/sr/cm-1)

    allocStatus(:) = 0
 
    allocate ( btObsErr(nchn),                           stat= allocStatus(1))
    allocate ( btObs(nchn),                              stat= allocStatus(2))
    allocate ( btCalc(nchn),                             stat= allocStatus(3))
    allocate ( rcal_clr(nchn),                           stat= allocStatus(4))
    allocate ( sfctau(nchn),                             stat= allocStatus(5))
    allocate ( cloudyRadiance(nchn,nlv_T-1),             stat= allocStatus(6))
    allocate ( transm(nchn,nlv_T),                       stat= allocStatus(7))
    allocate ( emi_sfc(nchn),                            stat= allocStatus(8))
    allocate ( radObs(nchn),                             stat= allocStatus(11))
    allocate ( rejflag(nchn,0:bitflag),                  stat= allocStatus(12))
    allocate ( ntop_bt(nchn),                            stat= allocStatus(13))
    allocate ( ntop_rd(nchn),                            stat= allocStatus(14))
    allocate ( ptop_bt(nchn),                            stat= allocStatus(15))
    allocate ( ptop_rd(nchn),                            stat= allocStatus(16))
    allocate ( minp(nchn),                               stat= allocStatus(17))
    allocate ( pmin(nchn),                               stat= allocStatus(18))
    allocate ( dtaudp1(nchn),                            stat= allocStatus(19))
    allocate ( fate(nchn),                               stat= allocStatus(20))
    if (liasi) allocate ( cloudyRadiance_avhrr(NIR,nlv_T-1), stat= allocStatus(21))
    allocate ( maxwf(nchn),                              stat= allocStatus(22))
    allocate ( pressure(nlv_T-1,1),                      stat= allocStatus(24))
    allocate ( tt(nlv_T-1),                              stat= allocStatus(25))
    allocate ( height(nlv_T-1,1),                        stat= allocStatus(26))
    allocate ( channelIndexes(nchn),                     stat= allocStatus(27))
    call utl_checkAllocationStatus(allocStatus, " irbg_doQualityControl 1")
  
    difftop_min = 100000.d0
    modelTopIndex = 1

    ptop_T = col_getPressure(columnTrlOnTrlLev,1,1,'TH')

    write(*,*) 'TOIT DU MODELE (MB)'
    write(*,*) 0.01d0 * ptop_T
    write(*,*) 'NIVEAU DU MODELE DE TRANSFERT RADIATIF LE PLUS PRES DU TOIT DU MODELE'
    write(*,*) modelTopIndex

    tvs_nobtov = 0

    ! Loop over all header indices of the 'TO' family
    call obs_set_current_header_list(obsSpaceData, 'TO')
    HEADER_2: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER_2

      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)

      if ( tvs_isIdBurpTovs(idatyp) ) tvs_nobtov = tvs_nobtov + 1
      if ( tvs_tovsIndex(headerIndex) < 0) cycle HEADER_2
      if ( tvs_isIdBurpInst(idatyp,instrumentName) .and. tvs_lsensor(tvs_tovsIndex (headerIndex)) == id) then
        btObs(:)    = -1.d0
        btCalc(:)   = -1.d0
        btObsErr(:) = -1.d0
        rcal_clr(:) = -1.d0
        sfctau(:)   = -1.d0
        cloudyRadiance(:,:)   = -1.d0
        transm(:,:) = -1.d0
        emi_sfc(:)  = -1.d0
        rejflag(:,:) = 0
        channelIndexes(:) = -1

        if (liasi) then
          classIndex = 1
          do columnIndex = OBS_CF1, OBS_CF7
            avhrr_bgck(headerIndex)%cfrac(classIndex) = obs_headElem_r(obsSpaceData,columnIndex,headerIndex)
            classIndex = classIndex + 1
          end do
          classIndex = 1
          channelNumber = 1
          do columnIndex = OBS_M1C1, OBS_M7C6
            avhrr_bgck(headerIndex)%radmoy(classIndex,channelNumber) = obs_headElem_r(obsSpaceData,columnIndex,headerIndex)
            channelNumber = channelNumber + 1
            if (channelNumber > nChanAVHRR) then
              channelNumber = 1
              classIndex = classIndex + 1
            end if
          end do
          classIndex = 1
          channelNumber = 1
          do columnIndex = OBS_S1C1, OBS_S7C6
            avhrr_bgck(headerIndex)%radstd(classIndex,channelNumber) = obs_headElem_r(obsSpaceData,columnIndex,headerIndex)
            channelNumber =channelNumber  + 1
            if (channelNumber> nChanAVHRR) then
              channelNumber = 1
              classIndex = classIndex + 1
            end if
          end do
          sunAzimuthAngle = obs_headElem_r(obsSpaceData,OBS_SAZ,headerIndex)
          sunAzimuthAnglePresent = ( abs(sunAzimuthAngle - MPC_missingValue_R8) > 0.01 ) 
        end if

        tg = col_getElem(columnTrlOnTrlLev, 1, headerIndex, 'TG')
        p0 = col_getElem(columnTrlOnTrlLev, 1, headerIndex, 'P0') * MPC_MBAR_PER_PA_R8

        do levelIndex = 1, nlv_T - 1
          tt(levelIndex) = col_getElem(columnTrlOnTrlLev, levelIndex+1, headerIndex, 'TT')
          height(levelIndex,1) = col_getHeight(columnTrlOnTrlLev, levelIndex+1, headerIndex, 'TH')
          pressure(levelIndex,1)= col_getPressure(columnTrlOnTrlLev, levelIndex+1, headerIndex, 'TH') * MPC_MBAR_PER_PA_R8
        end do
        qs = col_getElem(columnTrlOnTrlLev, nlv_T, headerIndex, 'HU')

        bodyStart   = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
        bodyEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyStart - 1
        bad = .false.
        if (lcris) bad=( obs_headElem_i(obsSpaceData, OBS_GQF, headerIndex)/=0 .or. &
             obs_headElem_i(obsSpaceData, OBS_GQL, headerIndex) /=0)
        if (liasi) bad=( obs_headElem_i(obsSpaceData, OBS_GQF, headerIndex)/=0 .or. &
             obs_headElem_i(obsSpaceData, OBS_GQL, headerIndex) >1) 

        nchannels = 0 ! number of channels available at that observation point
        do bodyIndex= bodyStart, bodyEnd
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) then
            channelNumber = nint(obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex))
            channelNumber = max( 0,min( channelNumber,tvs_maxChannelNumber + 1 ) )
            call tvs_getLocalChannelIndexFromChannelNumber(id,channelIndex,channelNumber)
            if (channelIndex == -1) cycle
            nchannels = nchannels + 1
            channelIndexes(nchannels) = channelIndex
            btObsErr(channelIndex) = obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)
            btObs(channelIndex) = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)

            !  Flag check on observed BTs ***
            if (.not.liasi .and. btest(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),2)) rejflag(channelIndex,9) = 1
            if (bad) rejflag(channelIndex,9) = 1

            !  Gross check on observed BTs ***
            if (btObs(channelIndex)<150.d0) rejflag(channelIndex,9) = 1
            if (btObs(channelIndex)>350.d0) rejflag(channelIndex,9) = 1
          end if
        end do

        if (nchannels==0) cycle HEADER_2

        do jc = 1, nchannels
          channelIndex = channelIndexes(jc)
          if ( channelIndex ==-1) cycle
          btCalc(channelIndex) = tvs_radiance(tvs_nobtov) % bt(channelIndex)
          rcal_clr(channelIndex) = tvs_radiance(tvs_nobtov) % clear(channelIndex)
          sfctau(channelIndex) = tvs_transmission(tvs_nobtov) % tau_total(channelIndex)
          do levelIndex = 1, nlv_T 
            transm(channelIndex,levelIndex) = tvs_transmission(tvs_nobtov) % tau_levels(levelIndex, channelIndex)
          end do          
          do levelIndex = 1, nlv_T - 1
            cloudyRadiance(channelIndex,levelIndex) = tvs_radiance(tvs_nobtov) % overcast(levelIndex, channelIndex)
          end do
          emi_sfc(channelIndex) = tvs_emissivity(channelIndex,tvs_nobtov)
          !  Gross check on computed BTs ***
          if (btCalc(channelIndex) < 150.d0) rejflag(channelIndex,9) = 1
          if (btCalc(channelIndex) > 350.d0) rejflag(channelIndex,9) = 1
        end do

        ksurf = profiles(tvs_nobtov) % skin % surftype

        !Test pour detecter l angle zenithal  manquant (-1) ou anormal
        ! (angle negatif ou superieur a 75 degres )
        satelliteZenithAngle= obs_headElem_r(obsSpaceData, OBS_SZA, headerIndex)
        if ( satelliteZenithAngle < 0 .or. satelliteZenithAngle > 75. ) then
          rejflag(:,9) = 1
          satelliteZenithAnglePresent = .false.
        else
          satelliteZenithAnglePresent = .true.
        end if
        clfr = 0.
        if (lairs .or. lcris) clfr = obs_headElem_r(obsSpaceData, OBS_CLF, headerIndex)

        sunZenithAngle = profiles(tvs_nobtov) % sunzenangle
        sunZenithAnglePresent = ( abs(sunZenithAngle - MPC_missingValue_R8) > 0.01 ) 

        if (liasi) then
          satelliteAzimuthAngle = profiles(tvs_nobtov) % azangle
          satelliteAzimuthAnglePresent = ( abs(satelliteAzimuthAngle - MPC_missingValue_R8) > 0.01 ) 
        end if
        albedo =  tvs_surfaceParameters(tvs_nobtov) % albedo
        ice =  tvs_surfaceParameters(tvs_nobtov) % ice
        ltype =  tvs_surfaceParameters(tvs_nobtov) % ltype
        if (ltype == 20) ksurf = 2
        pcnt_wat =  tvs_surfaceParameters(tvs_nobtov) % pcnt_wat
        pcnt_reg =  tvs_surfaceParameters(tvs_nobtov) % pcnt_reg
           
        !  Find TOA radiances converted from observed BT's

        radObs(:) = -1.d0
        
        channels: do jc = 1, nchannels
          channelIndex = channelIndexes(jc)
          if (channelIndex == -1) cycle channels
          if ( rejflag(channelIndex,9) == 1 ) cycle channels
          t_effective =  tvs_coefs(id) % coef % ff_bco(channelIndex) &
               + tvs_coefs(id) % coef % ff_bcs(channelIndex) * btObs(channelIndex)
          radObs(channelIndex) =  tvs_coefs(id) % coef % planck1(channelIndex) / &
               ( exp( tvs_coefs(id) % coef % planck2(channelIndex) / t_effective ) - 1.d0 )
        end do channels

        !  Set height fields to 'height above ground' fields
        do levelIndex = 1, nlv_T - 1
          height(levelIndex,1) = height(levelIndex,1) - height(nlv_T-1,1)
        end do

        ! ///// ---------------------------------------------------- /////
        ! ///// Determination of the clear/cloudy profiles (cldflag) /////
        ! ///// ---------------------------------------------------- /////           
        cldflag = 0
        
        ! Reference for window channel
        call tvs_getLocalChannelIndexFromChannelNumber(id, iwindo, iwindow(qcid) )
        call tvs_getLocalChannelIndexFromChannelNumber(id, iwindo_alt, iwindow_alt(qcid) )
        ichref = -1
        if (iwindo /= -1) then
          if ( rejflag(iwindo,9) == 0) then
            ichref = iwindo
          end if
        end if
        if (ichref == -1 .and. iwindo_alt /= -1) then
          if ( rejflag(iwindo_alt,9) == 0) then
            ichref = iwindo_alt
          end if
        end if

        if (ichref == -1) then
          cldflag = -1
          rejflag(:,9) = 1
          write(*,*) 'WARNING'
          write(*,*) 'WINDOW AND ALTERNATE WINDOW CHANNEL OBSERVATIONS'
          write(*,*) 'HAVE BEEN REJECTED.                             '
          write(*,*) 'ALL '//instrumentName//' OBSERVATIONS FROM THIS PROFILE REJECTED'
        end if

        !  -- Cloud top based on matching observed brightness temperature 
        !  -- at a reference surface channel with background temperature profile (ptop_eq)
        !  -- on guess vertical levels.

        lev_start = 0

        !iopt2=1 : calcul de la hauteur en hPa ptop_mb et du ntop_mb correspondant
        call cloud_height ( ptop_mb, ntop_mb, btObs, cldflag, tt, &
             height(:,1), p0, pressure(:,1), ichref, lev_start, iopt2 )

        !iopt1=2 : calcul de la hauteur em metres ptop_eq et du ntop_eq correspondant
        call cloud_height ( ptop_eq, ntop_eq, btObs, cldflag, tt, &
             height(:,1), p0, pressure(:,1), ichref, lev_start, iopt1 )

        if (liasi) then
          ! appel de RTTOV pour calculer les radiances des 3 canaux IR (3b, 4 et 5) de AVHRR 3
           
          call get_avhrr_emiss(emi_sfc(channelIndexes(1:nchannels)),tvs_coefs(id) % coef % ff_cwn(channelIndexes(1:nchannels)), &
               nchannels,avhrr_surfem1)

          call tovs_rttov_avhrr_for_IASI(headerIndex,avhrr_surfem1,tvs_satellites(id))
                 
          !The value computed will be used only if sunZenithAnglePresent is true
          call convert_avhrr(sunZenithAngle, avhrr_bgck(headerIndex) )
          call stat_avhrr(avhrr_bgck(headerIndex))
          
          lev_start_avhrr(:) = 0
          cldflag_avhrr(:) = 0
          do classIndex = 1, nClassAVHRR
            btObs_avhrr(:,classIndex) = avhrr_bgck(headerIndex) % tbmoy(classIndex,:)
            radObs_avhrr(1:NIR,classIndex) = avhrr_bgck(headerIndex) % radmoy(classIndex,NVIS + 1:NIR + NVIS)
            rcal_clr_avhrr(:) = avhrr_bgck(headerIndex) % radclearcalc(:)
            emi_sfc_avhrr(:) = avhrr_bgck(headerIndex) % emiss(:)
            sfctau_avhrr(:) = avhrr_bgck(headerIndex) % transmsurf(:)
            
            do levelIndex = 1, nlv_T - 1
              cloudyRadiance_avhrr(:,levelIndex) = avhrr_bgck(headerIndex) % radovcalc(levelIndex,:)
            end do
           
            if (btObs_avhrr(2,classIndex) > 100.d0 ) then
              ichref_avhrr(classIndex) = 2
            else if (btObs_avhrr(3,classIndex) > 100.d0 ) then
              ichref_avhrr(classIndex) = 3
            else
              ichref_avhrr(classIndex) = -1
              cldflag_avhrr(classIndex) = -1
            end if
            
            call cloud_height (ptop_eq_avhrr(classIndex),ntop_eq_avhrr(classIndex), btObs_avhrr(:,classIndex),cldflag_avhrr(classIndex),tt, &
                 height(:,1),p0,pressure,ichref_avhrr(classIndex),lev_start_avhrr(classIndex),iopt1)
          end do
          
        end if

        !  -- Clear/cloudy profile detection using the garand & nadon algorithm

        call garand1998nadon (cldflag, btObs,tg,tt, &
             height(:,1),ptop_eq,ntop_eq,ichref)

        if (liasi) then
          do classIndex=1,nClassAVHRR
            call garand1998nadon (cldflag_avhrr(classIndex), btObs_avhrr(:,classIndex),tg,tt, &
                 height(:,1),ptop_eq_avhrr(classIndex),ntop_eq_avhrr(classIndex),ichref_avhrr(classIndex))
          end do
        end if
        
        ! Further tests to remove potential cloudy profiles
        !  Test # A 
        !  In daytime, set cloudy if cloud fraction over 5% 
        cfsub = -1.d0
        if (lairs) then
          if ( cldflag == 0 .and. clfr > 5.d0 .and. sunZenithAngle < 90.d0 .and. sunZenithAnglePresent) then
            cldflag = 1
            cfsub = 0.01d0 * clfr !conversion % -> 0-1
          end if
        end if
        if (lcris) then
          if ( cldflag == 0 .and. clfr >= 0.d0 .and. clfr > crisCloudFractionThreshold &
               .and. crisCloudFractionThreshold >= 0.d0 ) then
            cldflag = 1
            cfsub = 0.01d0 * clfr !conversion % -> 0-1
          end if
        end if
        !  Test # B 
        !  Set cloudy if temperature difference between guess (tg)     
        !  and estimated true (tskinRetrieved) skin temperatures is over threshold 
        
        call estim_ts(tskinRetrieved, tg, emi_sfc, rcal_clr, radObs, &
             sfctau, cldflag, ichref, tvs_coefs(id) )

        if ( cldflag == 0 .and. ksurf == 1 &
             .and. abs(tskinRetrieved-tg) > dtw ) cldflag = 1 

        if ( cldflag == 0 .and. ksurf /= 1 &
             .and. abs(tskinRetrieved-tg) > dtl ) cldflag = 1

        if (liasi) then

          do classIndex = 1, nClassAVHRR
            call estim_ts(tskinRetrieved_avhrr(classIndex), tg,emi_sfc_avhrr,rcal_clr_avhrr,radObs_avhrr(:,classIndex), &
                 sfctau_avhrr,cldflag_avhrr(classIndex),ichref_avhrr(classIndex), coefs_avhrr)
          end do

          do classIndex = 1, nClassAVHRR
            if ( cldflag_avhrr(classIndex) == 0 .and. ksurf == 1 &
                 .and. abs(tskinRetrieved_avhrr(classIndex)-tg) > dtw ) cldflag_avhrr(classIndex) = 1
              
            if ( cldflag_avhrr(classIndex) == 0 .and. ksurf /= 1 &
                 .and. abs(tskinRetrieved_avhrr(classIndex)-tg) > dtl ) cldflag_avhrr(classIndex) = 1
              
          end do

          !criteres AVHRR utilisant les canaux visibles (de jour seulement)
          if ( sunZenithAngle < sunzenmax .and. sunZenithAnglePresent .and. satelliteAzimuthAnglePresent .and. &
               sunAzimuthAnglePresent .and. satelliteZenithAnglePresent) then 
            anisot = 1.d0
            
            if (albedo < albedoThresholdQC) then
              deltaphi = abs(satelliteAzimuthAngle - sunAzimuthAngle )
              if (deltaphi > 180.d0) deltaphi = 360.d0 - deltaphi
              call visocn(sunZenithAngle,satelliteZenithAngle,deltaphi,anisot,zlamb,zcloud,IER)
              albedoThreshold = 10.d0 * max(1.d0,anisot) 
            else
              albedoThreshold = 100.d0 * albedo + 10.d0
            end if
              
            if (anisot < 1.5d0) then !to avoid sun glint
              scos = cos ( sunZenithAngle * MPC_DEGREES_PER_RADIAN_R8 )
              call  cor_albedo (del, scos )
              albedoThreshold = albedoThreshold * del
              do classIndex = 1, nClassAVHRR
                if (avhrr_bgck(headerIndex)%albedmoy(classIndex,1) > albedoThreshold(1) ) then
                  cldflag_avhrr(classIndex) = 1
                end if
                !static AVHRR thresholds v3
                do channelIndex = 1, NVIS
                  if (avhrr_bgck(headerIndex)%albedmoy(classIndex,channelIndex) > seuilalb_static(channelIndex,ksurf) ) then
                    cldflag_avhrr(classIndex) = 1
                  end if
                end do
              end do
             
            end if
          end if

          ! Calcul de la pseudo fraction nuageuse AVHRR

          cfrac_avhrr = 0.d0
          do classIndex = 1, nClassAVHRR
            if (cldflag_avhrr(classIndex) == 1) cfrac_avhrr = cfrac_avhrr + avhrr_bgck(headerIndex) % cfrac(classIndex)
          end do

          cfsub = -1.0d0
          if ( cldflag == 0 .and. cfrac_avhrr > 5.d0 ) then
            cldflag = 1
            cfsub = 0.01d0 * min(cfrac_avhrr, 100.d0) !conversion % -> 0-1 avec seuil car parfois cfrac_avhrr=101
          end if

          !AVHRR Homogeneity criteria
          if (cldflag == 0 .and. sunZenithAnglePresent) then
            ijour = 1
            if (sunZenithAngle < 90.d0) ijour=2
            ! 1 NUIT
            ! 2 JOUR
            if (ijour == 2) then
              do channelIndex=1,NVIS
                if (avhrr_bgck(headerIndex)%albstd_pixeliasi(channelIndex) > seuilalb_homog(channelIndex,ksurf) ) cldflag = 1
              end do
            end if
            do channelIndex=NVIS+1,NVIS+NIR
              if (avhrr_bgck(headerIndex)%tbstd_pixelIASI(channelIndex) > seuilbt_homog(channelIndex,ksurf,ijour)) cldflag = 1
            end do
          end if
        end if

        gncldflag = cldflag

        ! ///// ------------------------------------------------------- /////
        ! ///// DETERMINATION OF THE ASSIMILABLE OBSERVATIONS (rejflag) /////
        ! ///// ------------------------------------------------------- /////


        !  -- FIRST TestS TO REJECT OBSERVATIONS


        ! *** Test # 1 ***
        ! *** Do not assimilate where cloudy ***

        if ( cldflag == 1 ) then
          rejflag(:,11) = 1
          rejflag(:,23) = 1
        end if

        ! *** Test # 2 ***
        ! *** Gross check on valid BTs ***

        !     already done


        ! -- Cloud top based on matching 
        ! -- observed brightness temperature with background temperature profiles (ptop_bt)
        ! -- or computed observed radiances with background radiance profiles (ptop_rd)
        ! -- on rttov vertical levels

        lev_start = 0

        do channelIndex = 1, nch_he
          call tvs_getLocalChannelIndexFromChannelNumber(id,ilist_HE(channelIndex),ilist1(qcid,channelIndex))
        end do
        !here the case for which one of the channelindex is -1 is treated in cloud_top
        call cloud_top ( ptop_bt,ptop_rd,ntop_bt,ntop_rd, &
             btObs,tt,height(:,1),rcal_clr,p0,radObs,cloudyRadiance,pressure(:,1), &
             cldflag, lev_start, iopt2, ihgt, ilist_he,rejflag_opt=rejflag,ichref_opt=ichref)

        if (liasi) then
          lev_start_avhrr(:) = 0
          do classIndex = 1, nClassAVHRR
            call cloud_top( ptop_bt_avhrr(:,classIndex),ptop_rd_avhrr(:,classIndex),ntop_bt_avhrr(:,classIndex),ntop_rd_avhrr(:,classIndex), &
                 btObs_avhrr(:,classIndex),tt,height(:,1),rcal_clr_avhrr,p0,radObs_avhrr(:,classIndex),cloudyRadiance_avhrr,pressure(:,1), &
                 cldflag_avhrr(classIndex),lev_start_avhrr(classIndex),iopt2,ihgt,ilist_avhrr)
          end do
        end if

        !  -- reference channel for co2-slicing

        do channelIndex = 1, nco2
          call tvs_getLocalChannelIndexFromChannelNumber(id, ilist_co2(channelIndex), ilist2(qcid,channelIndex)  )
          call tvs_getLocalChannelIndexFromChannelNumber(id, ilist_co2_pair(channelIndex), ilist2_pair(qcid,channelIndex)  )
        end do

        countInvalidChannels = 0
        do channelIndex = 1, nco2
          if (  ilist_co2(channelIndex) == -1 .or. ilist_co2_pair(channelIndex) == -1 ) then
            countInvalidChannels = countInvalidChannels + 1
          else if ( rejflag(ilist_co2(channelIndex),9) == 1 .or. &
               rejflag(ilist_co2_pair(channelIndex),9) == 1 ) then
            countInvalidChannels = countInvalidChannels + 1
          end if
        end do
         
        if (countInvalidChannels == nco2) then
          cldflag = -1
          rejflag(:,9) = 1
          write(*,*) 'WARNING'
          write(*,*) 'CO2 REFERENCE AND ALTERNATE CHANNEL OBSERVATIONS'
          write(*,*) 'HAVE BEEN REJECTED.                             '
          write(*,*) 'ALL '//instrumentName//' OBSERVATIONS FROM THIS PROFILE REJECTED'
        end if

        !   Equivalent height of selected window channel
        heff = MPC_missingValue_R8
        call tvs_getLocalChannelIndexFromChannelNumber(id,channelIndex,ilist1(qcid,2))
        if (channelIndex /= -1) heff = ptop_rd( channelIndex )
              
        if (ichref == iwindo_alt) then
          call tvs_getLocalChannelIndexFromChannelNumber(id,channelIndex,ilist1(qcid,3))
          if (channelIndex /= -1) heff = ptop_rd( channelIndex )
        end if
        !  Cloud top based on co2 slicing 

        co2min = minloc( abs( pressure(:,1) - pco2min ) )
        co2max = minloc( abs( pressure(:,1) - pco2max ) )
        
        lev_start = max( min(lev_start,co2max(1)), co2min(1) )

        call co2_slicing ( ptop_co2, ntop_co2, fcloud_co2, &
             rcal_clr, cloudyRadiance, radObs, p0, pressure(:,1), cldflag, rejflag, &
             lev_start, ichref, ilist_co2, ilist_co2_pair)

        !  -- Find consensus cloud top and fraction
 
        call seltop ( etop,vtop,ecf,vcf,ngood, heff,ptop_co2,fcloud_co2, &
             cfsub,ptop_mb,p0,cldflag,gncldflag )

        if (liasi) then
          ! Correction pour les nuages trop bas:
          ! en principe Pco2 < Heff.
          ! on cherche les cas pathologiques avec Pco2>Min(Heff(AVHRR))
          minpavhrr(2:3) = 12200
          iloc(2:3) = -1      ! pour eviter les catastrophes...
          do classIndex=1,nClassAVHRR
            if (avhrr_bgck(headerIndex)%cfrac(classIndex) > 0.d0) then
              if (ptop_rd_avhrr(2,classIndex) < minpavhrr(2)) then
                iloc(2) = classIndex
                minpavhrr(2) = ptop_rd_avhrr(2,classIndex)
              end if
              if (ptop_rd_avhrr(3,classIndex) < minpavhrr(3)) then
                iloc(3) = classIndex
                minpavhrr(3) = ptop_rd_avhrr(3,classIndex)
              end if
            end if
          end do
          if ( iloc(2) /= -1 .and. iloc(3) /= -1) then ! pour eviter les catastrophes...
            ! on se limite aux cas "surs" ou les deux hauteurs effectives sont > a Pco2
            ! et ou un accord raisonnable existe entre les deux hauteurs effectives
            if ( iloc(2) == iloc(3) .and. &
                 minpavhrr(2) < etop .and. &
                 minpavhrr(3) < etop .and. &
                 abs(minpavhrr(2)- minpavhrr(3)) < 25.d0 .and. &
                 cldflag_avhrr(iloc(2)) /= -1 .and. cldflag_avhrr(iloc(3)) /= -1) then
              
              if (ecf == 0.d0 .and. cldflag == 1) then
                ! cas predetermine nuageux mais ramene a clair 
                ecf = 0.01d0 * min(100.d0,cfrac_avhrr)
                ! cette ligne peut generer des fractions nuageuses inferieures a 20 %.
                etop = 0.5d0 * (minpavhrr(2) + minpavhrr(3))
              end if

              if (ecf > 0.d0 .and. cldflag == 1) then
                !cas predetermine nuageux pas ramene clair (==normal)
                etop = 0.5d0 * ( minpavhrr(2) + minpavhrr(3))
              end if

              if (cldflag == 0) then
                !cas predetermine clair ... que faire
                cldflag = 1
                etop = 0.5d0 * (minpavhrr(2) + minpavhrr(3))
                ecf = 0.01d0 * min(100.d0,cfrac_avhrr)
              end if
            end if
          end if
        end if

        !  -- Find minimum level of sensitivity for channel assimilation not sensible to clouds        
        call min_pres_new (maxwf, minp, pmin, dtaudp1, p0, transm(:,2:nlv_T), pressure(:,1), cldflag, modelTopIndex)
        !  -- ASSIMILATION OF OBSERVATIONS WHEN CLOUDY PROFILES

        ! *** Test # 3 ***
        ! *** Assimilation above clouds (refinement of test 1)             ***
        ! *** Set security margin to 2x the std on height from CO2-slicing *** 

        tampon = max(50.d0, 2.d0*vtop)                                                          

        do channelIndex = 1, nchn        
          if ( rejflag(channelIndex,11) == 1 .and. rejflag(channelIndex,23) == 1 .and. etop - tampon > pmin(channelIndex) ) then
            rejflag(channelIndex,11) = 0
            rejflag(channelIndex,23) = 0
          end if
        end do

        !     Look at the fate of the observations
        fate(:) = sum(rejflag(:,:), DIM=2)            

        !     Further reasons to reject observations
        do channelIndex = 1, nchn

          if ( fate(channelIndex) == 0 ) then

            ! *** Test # 4 ***
            ! *** Background check, do not assimilate if O-P > 3sigma ***

            if ( abs(btObs(channelIndex) - btCalc(channelIndex)) > 3.d0 * btObsErr(channelIndex) ) then
              rejflag(channelIndex,9) = 1
              rejflag(channelIndex,16) = 1
            end if

            ! *** Test # 5 ***
            ! *** Do not assimilate shortwave channels during the day ***

            if ( tvs_coefs(id) % coef % ff_cwn(channelIndex) >= 2000.d0 .and. sunZenithAngle < night_ang .and. sunZenithAnglePresent) then
              rejflag(channelIndex,11) = 1
              rejflag(channelIndex,7)  = 1
            end if

            ! *** Test # 6 ***
            ! *** Do not assimilate surface channels over land ***

            if ( minp(channelIndex) == nlv_T .or. p0-pmin(channelIndex) < 100.d0 ) then
              if ( ksurf == 0 ) then
                rejflag(channelIndex,11) = 1    !!! comment this line if assimilation under conditions
                rejflag(channelIndex,19) = 1    !!! comment this line if assimilation under conditions
                if ( pcnt_wat > 0.01d0 .or. pcnt_reg > 0.1d0 .or. emi_sfc(channelIndex) < 0.97d0 ) then
                  rejflag(channelIndex,11) = 1
                  rejflag(channelIndex,19) = 1
                end if

                ! *** Test # 7 ***
                ! *** Do not assimilate surface channels over water under conditions ***

              else if ( ksurf == 1 ) then
                if ( pcnt_wat < 0.99d0 .or. pcnt_reg < 0.97d0 .or. &
                     ice > 0.001d0 .or. albedo >= albedoThresholdQC .or. emi_sfc(channelIndex) < 0.9d0 ) then
                  rejflag(channelIndex,11) = 1   
                  rejflag(channelIndex,19) = 1   
                end if
                
                ! *** Test # 8 ***
                ! *** Do not assimilate surface channels over sea ice ***
                          
              else if ( ksurf == 2 ) then
                rejflag(channelIndex,11) = 1
                rejflag(channelIndex,19) = 1   
              end if
            end if
            
          end if

          ! *** Test # 9 ***
          ! *** Do not assimilate if jacobian has a significant contribution over model top ***

          ! Condition valid if model top at 10mb or lower only
          if ( nint(ptop_T) >= 1000 ) then
            if ( rejflag(channelIndex,9) /= 1 .and. dtaudp1(channelIndex) > 0.50d0 ) then
              rejflag(channelIndex,11) = 1
              rejflag(channelIndex,21) = 1
            end if
          end if
        
          ! Condition valid if model top at 10mb or lower only
          if ( nint(ptop_T) >= 1000 ) then
            if ( rejflag(channelIndex,9) /= 1 .and. transm(channelIndex,1) < 0.99d0 ) then
              rejflag(channelIndex,11) = 1
              rejflag(channelIndex,21) = 1 
            end if
          end if

          ! Condition valid if model top is higher than 10 mb
          ! with computation made on model levels instead of RTTOV levels
          ! this test is a litle bit more strict than before
          ! decrease of the order of  1-2% in the number of high peaking radiance IR assimilated
          ! not dramatic so could be kept as is
          if ( nint(ptop_T) < 1000 ) then
            if ( rejflag(channelIndex,9) /= 1 .and. transm(channelIndex,1) < 0.95d0 ) then
              rejflag(channelIndex,11) = 1
              rejflag(channelIndex,21) = 1 
            end if
          end if

        end do

        nchannels =0 
        do bodyIndex= bodyStart, bodyEnd
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) then
            nchannels =  nchannels + 1
            if (btest(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),8)) rejflag(channelIndexes(nchannels),8) = 1
          end if
        end do

        !  For each profile, are all non-blacklisted channels assimilated

        assim_all = .true.
        fate(:) = sum(rejflag(:,:),DIM=2)            
        
        chn: do channelIndex = 1, nchn
          if ( rejflag(channelIndex,8) == 0 ) then
            if ( fate(channelIndex) /= 0 ) then
              assim_all = .false.
              exit chn
            end if
          end if
        end do chn

        if  (.not.assim_all) then
          call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex,ibset(obs_headElem_i(obsSpaceData,OBS_ST1,headerIndex),6) )
        end if
        !  -- Addition of background check parameters to burp file

        call obs_headSet_r(obsSpaceData, OBS_ETOP, headerIndex, etop )
        call obs_headSet_r(obsSpaceData, OBS_VTOP, headerIndex, vtop )
        call obs_headSet_r(obsSpaceData, OBS_ECF,  headerIndex, 100.d0 * ecf )
        call obs_headSet_r(obsSpaceData, OBS_VCF,  headerIndex, 100.d0 * vcf )
        call obs_headSet_r(obsSpaceData, OBS_HE,   headerIndex, heff )
        call obs_headSet_r(obsSpaceData, OBS_ZTSR, headerIndex, tskinRetrieved )
        call obs_headSet_i(obsSpaceData, OBS_NCO2, headerIndex, ngood)
        call obs_headSet_r(obsSpaceData, OBS_ZTM,  headerIndex, tt(nlv_T-1) )
        call obs_headSet_r(obsSpaceData, OBS_ZTGM, headerIndex, tg )
        call obs_headSet_r(obsSpaceData, OBS_ZLQM, headerIndex, qs )
        call obs_headSet_r(obsSpaceData, OBS_ZPS,  headerIndex, 100.d0 * p0 )
        call obs_headSet_i(obsSpaceData, OBS_STYP, headerIndex, ksurf )

        do bodyIndex = bodyStart, bodyEnd
          channelNumber = nint(obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex))
          channelNumber = max(0, min(channelNumber, tvs_maxChannelNumber + 1))
          call tvs_getLocalChannelIndexFromChannelNumber(id,channelIndex,channelNumber)
          if (channelIndex == -1) cycle
          call obs_bodySet_r(obsSpaceData,OBS_SEM,bodyIndex,emi_sfc(channelIndex))
          do bitIndex = 0, bitflag
            if ( rejflag(channelIndex,bitIndex) == 1 ) &
                 call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex,IBSET(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),bitIndex))
          end do
        end do
          
      end if

    end do HEADER_2

    deallocate ( channelIndexes,           stat= allocStatus(1))
    deallocate ( height,                   stat= allocStatus(2))
    deallocate ( tt,                       stat= allocStatus(3))
    deallocate ( pressure,                 stat= allocStatus(4))
    deallocate ( maxwf,                    stat= allocStatus(5))
    if (liasi) deallocate ( cloudyRadiance_avhrr , stat= allocStatus(6))
    deallocate ( fate,                     stat= allocStatus(7))
    deallocate ( dtaudp1,                  stat= allocStatus(8))
    deallocate ( pmin,                     stat= allocStatus(9))
    deallocate ( minp,                     stat= allocStatus(10))
    deallocate ( ptop_rd,                  stat= allocStatus(11))
    deallocate ( ptop_bt,                  stat= allocStatus(12))
    deallocate ( ntop_rd,                  stat= allocStatus(13))
    deallocate ( ntop_bt,                  stat= allocStatus(14))
    deallocate ( rejflag,                  stat= allocStatus(15))
    deallocate ( radObs,                   stat= allocStatus(16))
    deallocate ( emi_sfc,                  stat= allocStatus(17))
    deallocate ( transm,                   stat= allocStatus(18))
    deallocate ( cloudyRadiance,           stat= allocStatus(19))
    deallocate ( sfctau,                   stat= allocStatus(20))
    deallocate ( rcal_clr,                 stat= allocStatus(21))
    deallocate ( btCalc,                   stat= allocStatus(22))
    deallocate ( btObs,                    stat= allocStatus(23))
    deallocate ( btObsErr,                 stat= allocStatus(24))
    call utl_checkAllocationStatus(allocStatus, " irbg_doQualityControl", .false.)
    nullify(profiles)

  end subroutine irbg_doQualityControl


  !--------------------------------------------------------------------------
  ! convert_avhrr
  !--------------------------------------------------------------------------
  subroutine convert_avhrr(sunzen, avhrr)
    !
    ! :Purpose: conversion des radiance IR en temperatures de brillance
    !           et des radiances visibles en "albedo"
  
    implicit none

    ! Arguments:
    real(8),                intent(in)    :: sunzen ! Solar zenith angle
    type (avhrr_bgck_iasi), intent(inout) :: avhrr  ! Structure containing AVHRR observations

    ! Locals:
    integer :: classIndex
    real(8) :: bt(NIR), dbtsdrad(NIR)
    real(8) :: freq(NIR), offset(NIR), slope(NIR)

    freq = coefs_avhrr%coef%ff_cwn (:)
    offset = coefs_avhrr%coef%ff_bco(:)
    slope = coefs_avhrr%coef%ff_bcs(:)

    do classIndex=1, nClassAVHRR
      call calcbt(avhrr % radmoy(classIndex,4:6), bt, dbtsdrad, freq, offset, slope)
      avhrr % tbmoy(classIndex,4:6) = bt(1:3)
      avhrr % tbstd(classIndex,4:6) = avhrr % radstd(classIndex,4:6) * dbtsdrad(1:3)
      call calcreflect(avhrr % radmoy(classIndex,1:3), sunzen, avhrr % albedmoy(classIndex,1:3) )
      call calcreflect(avhrr % radstd(classIndex,1:3), sunzen, avhrr % albedstd(classIndex,1:3) )
    end do

  end subroutine convert_avhrr

  !--------------------------------------------------------------------------
  ! calcreflect
  !--------------------------------------------------------------------------
  subroutine calcreflect(rad, sunzen, reflect)
    !
    ! :Purpose: Computes Top of Atmosphere Albedo as defined by equation (4)
    !           of Rao et al. Int. J. of Remote Sensing, 2003, vol 24, no 9, 1913-1924
    !           
    implicit none

    ! Arguments:
    real(8), intent(in) :: rad(nvis)     ! radiances array
    real(8), intent(in) :: sunzen        ! Sun zenith angle
    real(8), intent(out):: reflect(nvis) ! TOA albedo en %

    ! Locals:
    real(8),parameter :: solar_filtered_irradiance(nvis) = (/139.873215d0,232.919556d0,14.016470d0/)
    ! equivalent widths, integrated solar irradiance,  effective central wavelength
    !0.084877,139.873215,0.632815
    !0.229421,232.919556,0.841679
    !0.056998,14.016470,1.606119
    real(8) :: radb ! radiance en W/m2/str
    integer :: i

    do i = 1, nvis
      if (rad(i) >= 0.0d0 ) then
        radb = rad(i) / 1000.0d0
        reflect(i) = (MPC_PI_R8 * radb) / solar_filtered_irradiance(i)
        if (sunzen < 90.0d0 ) reflect(i) = reflect(i) / cos(sunzen * MPC_RADIANS_PER_DEGREE_R8)
      else
        reflect(i) = -1
      end if
    end do
  
  end subroutine calcreflect

  !--------------------------------------------------------------------------
  ! calcbt
  !--------------------------------------------------------------------------
  subroutine calcbt(rad, tb, dtbsdrad, freq, offset, slope)
    !
    ! :Purpose: Computes brightness temperature (bt) and the first derivative of
    !            bt with respect to radiance, from radiance, frequencies
    !           
    !
    implicit none

    ! Arguments:
    real(8), intent(in)  :: rad(:)              ! Radiance
    real(8), intent(in)  :: freq(size(rad))     ! Channel wavenumber (cm-1)
    real(8), intent(in)  :: offset(size(rad))   ! 
    real(8), intent(in)  :: slope(size(rad))    !
    real(8), intent(out) :: tb(size(rad))       ! Brightness Temperature
    real(8), intent(out) :: dtbsdrad(size(rad)) ! Derivative of tb wrt radiance

    ! Locals:
    integer  :: channelIndex, nchan
    real(8)  :: radtotal, tstore, planck1, planck2
    real(8), parameter :: c1= 1.19106590d-05   ! First Planck constant
    real(8), parameter :: c2= 1.438833d0       ! Second Planck constant 

    nchan = size(rad)

    do channelIndex = 1, nchan
      if (rad(channelIndex) > 1.d-20) then
        planck2 = c2 * freq(channelIndex)
        planck1 = c1 * ( freq(channelIndex) ** 3 ) 
        tstore = planck2 / log( 1.0d0 + planck1 / rad(channelIndex) )
        tb(channelIndex) = ( tstore - offset(channelIndex) ) / slope(channelIndex)
        
        radtotal = rad(channelIndex)
        
        dtbsdrad(channelIndex) = planck1 * tstore ** 2 / ( planck2 * radtotal * ( radtotal + planck1 ) )
        
        dtbsdrad(channelIndex) = dtbsdrad(channelIndex) / slope(channelIndex)
        
      else
        tb(channelIndex) = 0.d0
        dtbsdrad(channelIndex) = 0.d0
      end if
      
    end do

  end subroutine calcbt

  !--------------------------------------------------------------------------
  ! stat_avhrr
  !--------------------------------------------------------------------------
  subroutine stat_avhrr( avhrr )
    !
    ! :Purpose: calcul de statistiques
    !           sur l'information sous-pixel AVHRR
    !
    implicit none

    ! Arguments:
    type (avhrr_bgck_iasi), intent(inout) :: avhrr

    ! Locals:
    integer :: classIndex, channelIndex
    real(8) :: sumFrac(nvis+nir), sumBt(nvis+1:nvis+nir),sumBt2(nvis+1:nvis+nir)
    real(8) :: sumAlb(1:nvis),sumAlb2(1:nvis)
    
    sumFrac(:) = 0.d0
    sumBt(:) = 0.d0
    sumBt2(:) = 0.d0
    sumAlb(:) = 0.d0
    sumAlb2(:) = 0.d0
    
    do classIndex = 1, nClassAVHRR
      if (avhrr%cfrac(classIndex) > 0.d0 ) then
        do channelIndex = 1, NVIS
          if (avhrr%albedmoy(classIndex,channelIndex) >= 0.d0 ) then
            sumFrac(channelIndex) = sumFrac(channelIndex) + avhrr%cfrac(classIndex)
            sumAlb(channelIndex) = sumAlb(channelIndex) + avhrr%cfrac(classIndex) * avhrr%albedmoy(classIndex,channelIndex)
            sumAlb2(channelIndex) = sumAlb2(channelIndex) + avhrr%cfrac(classIndex) * ( avhrr%albedmoy(classIndex,channelIndex)**2 + avhrr%albedstd(classIndex,channelIndex)**2)
          end if
        end do
        do channelIndex = 1+NVIS, NVIS+NIR
          if (avhrr%tbmoy(classIndex,channelIndex) > 0.d0 ) then
            sumFrac(channelIndex) = sumFrac(channelIndex) + avhrr%cfrac(classIndex)
            sumBt(channelIndex) = sumBt(channelIndex) + avhrr%cfrac(classIndex) * avhrr%tbmoy(classIndex,channelIndex)
            sumBt2(channelIndex) = sumBt2(channelIndex) + avhrr%cfrac(classIndex) * (avhrr%tbmoy(classIndex,channelIndex)**2 + avhrr%tbstd(classIndex,channelIndex)**2 )
          end if
        end do
      end if
    end do
      
    do channelIndex = 1, NVIS
      if (sumFrac(channelIndex) > 0.d0 ) then
        sumAlb(channelIndex) = sumAlb(channelIndex) / sumFrac(channelIndex)
        sumAlb2(channelIndex) = sumAlb2(channelIndex)/sumFrac(channelIndex) - sumAlb(channelIndex)**2
        if (sumAlb2(channelIndex) > 0.d0) then
          sumAlb2(channelIndex) = sqrt( sumAlb2(channelIndex) )
        else
          sumAlb2(channelIndex) = 0.d0
        end if
      end if
    end do
      
    do channelIndex = NVIS+1, NVIS+NIR
      if (sumFrac(channelIndex) > 0.d0 ) then
        sumBt(channelIndex) = sumBt(channelIndex) / sumFrac(channelIndex)
        sumBt2(channelIndex) = sumBt2(channelIndex)/sumFrac(channelIndex) - sumBt(channelIndex)**2
        if (sumBt2(channelIndex) > 0.d0) then
          sumBt2(channelIndex)= sqrt ( sumBt2(channelIndex) )
        else
          sumBt2(channelIndex) = 0.d0
        end if
      end if
    end do
      
    avhrr % tbstd_pixelIASI = sumBt2
    avhrr % albstd_pixeliasi = sumAlb2
      
  end subroutine stat_avhrr

  !--------------------------------------------------------------------------
  ! co2_slicing
  !--------------------------------------------------------------------------
  subroutine co2_slicing (ptop, ntop, fcloud,                              &
       rcal, cloudyRadiance, radObs, p0, plev, cldflag, rejflag, &
       lev_start, ichref, ilist, ilist_pair)
    !
    ! :Purpose: cloud top height computation.
    !           cloud top from co2 slicing and cloud fraction estimate
    !
    implicit none

    ! Arguments:
    real(8), intent(in)    :: rcal(:)                ! computed clear radiances (mW/m2/sr/cm-1)
    real(8), intent(in)    :: plev(:)                ! pressure levels (hPa)
    real(8), intent(in)    :: cloudyRadiance(size(rcal),size(plev)) ! computed cloud radiances from each level (mW/m2/sr/cm-1)
    real(8), intent(in)    :: radObs(size(rcal))              ! Observed radiances (mW/m2/sr/cm-1)
    real(8), intent(in)    :: p0                        ! surface pressure (hPa)
    integer, intent(in)    :: ichref                    ! window channel to predetermine clear
    integer, intent(in)    :: cldflag                   ! Cloudy flag (0 Clear, 1 Cloudy, -1 undefined)
    integer, intent(in)    :: rejflag(1:,0:)            ! flags for rejected observations
    integer, intent(in)    :: ilist(nco2)               ! first list of channels
    integer, intent(in)    :: ilist_pair(nco2)          ! second list channe
    integer, intent(inout) :: lev_start                 ! Level to start iteration (ideally tropopause)
    real(8), intent(out)   :: ptop(nco2)                ! Cloud top (hPa)
    real(8), intent(out)   :: fcloud(nco2)              ! Cloud fraction
    integer, intent(out)   :: ntop(nco2)                ! Nearest pressure level corresponding to ptop (ptop <= p0)

    ! Locals
    integer             :: j,jch,jc,jpmax,jmax, nlev,nchn
    integer             :: sumrej
    real(8)             :: rapg
    real(8),allocatable :: drap(:,:), a_drap(:), fc(:,:)
    real(8)             :: val,val1,val2,val3,fcint
    real(8)             :: emi_ratio
    integer             :: jc_pair
    integer             :: iter,niter
    real(8), parameter  :: eps = 1.D-12

    ptop(:) = -1.d0
    ntop(:) = -1
    fcloud(:) = -1.d0

    !  Profile not assimilated if data from 2 windows channels bad
    !  and/or if data from 2 reference co2 channels bad

    if ( cldflag == -1 ) return

    nlev = size(plev)
    nchn = size(rcal)

    !  Define closest level jpmax to surface pressure p0

    jpmax = nlev
    
    do J = lev_start, nlev
      if ( plev(J) > p0 ) then
        jpmax = J
        exit
      end if
    end do
    
    !     define jmax as last level for co2-slicing calculations
  
    jmax = jpmax - 1
    
    !     predetermined clear window channel, all nco2 estimates clear

    sumrej = sum(rejflag(ichref,:))

    if ( sumrej == 0 ) then
      ptop(:) = p0
      ntop(:) = jpmax
      fcloud(:) = 0.d0
      return
    end if

    allocate(fc(nchn,nlev), drap(nco2,nlev), a_drap(nlev) )

    channels: do jch = 1, nco2
      
      jc = ilist(jch)
      jc_pair = ilist_pair(jch)
      fc(jc_pair,:) = rcal(jc_pair) - cloudyRadiance(jc_pair,:)
      niter = 1
      if ( jch > 13) niter = 2 
     
      iteration: do iter = 1, niter
        drap(jch,:)   = 9999.d0
        ntop(jch) = -1
        !         calcul emi_ratio
        if (jch > 13) then       
          if ( iter == 1 ) then
            emi_ratio = 1.0376d0
          else
            emi_ratio = 1.09961d0 - 0.09082d0 * fcloud(jch)
          end if
        else
          emi_ratio = 1.0d0
        end if

        fc(jc,:) = rcal(jc) - cloudyRadiance(jc,:)

        !      Gross check failure

        if ( rejflag(jc,9) == 1 ) cycle channels
        if ( rejflag(jc_pair,9) == 1 ) cycle channels

        if ( abs( rcal(jc_pair) - radObs(jc_pair) ) > eps ) then
          rapg = (rcal(jc) - radObs(jc)) / (rcal(jc_pair) - radObs(jc_pair))
        else
          rapg = 0.0d0
        end if

        do J = lev_start, jpmax
          if ( fc(jc,J) > 0.d0 .and. fc(jc_pair,J) > 0.d0 )  &
               drap(jch,J) = rapg - (fc(jc,J) / fc(jc_pair,J)) * emi_ratio
        end do

        a_drap(:) = abs( drap(jch,:) )

        levels: do J = lev_start + 1, jmax

          ! *         do not allow fc negative (i.e. drap(jch,j) = 9999.)

          if ( drap(jch,J) > 9000.d0 .and. &
               a_drap(J-1) < eps .and. &
               a_drap(J+1) < eps ) cycle channels

          val = drap(jch,J) / drap(jch,J - 1)

          ! *         find first, hopefully unique, zero crossing

          if ( val < 0.d0 ) then

            ! *         conditions near zero crossing of isolated minimum need monotonically
            ! *         decreasing drap from j-3 to j-1 as well increasing from j to j+1

            val1 = drap(jch,J - 2) / drap(jch,J - 1)
            val2 = drap(jch,J - 3) / drap(jch,J - 1)
            val3 = drap(jch,J) / drap(jch,J + 1)

            if ( val1 > 0.d0 .and.  & 
                 val2 > 0.d0 .and.  & 
                 val3 > 0.d0 .and.  &
                 a_drap(J-2) > a_drap(J-1) .and.  &
                 a_drap(J-3) > a_drap(J-2) .and.  &
                 a_drap(J)   < 9000.d0     .and.  &
                 a_drap(J+1) > a_drap(J) )        &
                 then
              ptop(jch) = plev(J)
              ntop(jch) = J
            end if
            
            exit levels
                      
          end if
              
        end do levels

        j = ntop(jch)

        ! *       special cases of no determination

        if ( j < 1) then
          ptop(jch)   = -1.d0
          ntop(jch)   = -1
          fcloud(jch) = -1.d0
          cycle channels
        end if
      
        if ( J <= lev_start .or. drap(jch,J) > 9000.d0 ) then
          !if ( iter == 1) then
          ptop(jch) = -1.d0
          ntop(jch) = -1
          fcloud(jch) = -1.d0
          !end if
          cycle channels
        end if
        
        if ( abs( cloudyRadiance(jc,J) - rcal(jc) ) > 0.d0 )  &
             fcloud(jch) = (radObs(jc) - rcal(jc)) /  &
             (cloudyRadiance(jc,J) - rcal(jc))

        ! *       find passage to zero if it exists and interpolate to exact pressure

        ptop(jch) = plev(J - 1) - drap(jch,J - 1) /    &
             ( drap(jch,J) - drap(jch,J-1) ) * ( plev(J) - plev(J-1) )
        ! *       find cloud radiance at zero crossing to use to get cloud fraction

        fcint = fc(jc,J - 1) + ( fc(jc,J) - fc(jc,J - 1) ) /   &
             ( plev(J) - plev(J - 1) ) * ( ptop(jch) - plev(J - 1) )

        ! *       find cloud fraction based on exact cloud top

        if ( abs(fcint) > 0.d0 )             &
             fcloud(jch) = ( rcal(jc) - radObs(jc) ) / fcint

        fcloud(jch) = min ( fcloud(jch),  1.5d0 )
        fcloud(jch) = max ( fcloud(jch), -0.5d0 )

        if (fcloud(jch) < 0.0d0 .or. fcloud(jch) > 1.0d0 )  cycle channels
      
      end do iteration
     
    end do channels

    deallocate(fc, drap, a_drap )
      
  end subroutine co2_slicing

  !--------------------------------------------------------------------------
  ! seltop
  !--------------------------------------------------------------------------
  subroutine seltop (etop, vtop, ecf, vcf, ngood, he, ht, cf, cfsub, ptop_mb, p0, cldflag, gncldflag)
    !
    ! :Purpose: Select cloud top by averaging co2-slicing results
    !           judged correct. all missing values are -1.
    !
    implicit none

    ! Arguments:
    real(8), intent(in)  :: he       ! equivalent cloud top heights from a window channel (hPa)
    real(8), intent(in)  :: ht(nco2) ! cloud tops from co2-slicing (hPa)
    real(8), intent(in)  :: cf(nco2) ! effective cloud fraction for co2-slicing
    real(8), intent(in)  :: p0       ! surface pressure in (hPa)
    real(8), intent(in)  :: cfsub    ! visible ("subpixel") cloud fraction
    integer, intent(in)  :: cldflag  ! Cloudy flag (0 Clear, 1 Cloudy, -1 undefined)
    integer, intent(in)  :: gncldflag! Garand Nadon cloudy flag (0 Clear, 1 Cloudy, -1 undefined)
    real(8), intent(out) :: etop     ! consensus cloud top (hPa)
    real(8), intent(out) :: vtop     ! corresponding variance on etop (hPa)
    real(8), intent(out) :: ecf      ! consensus effective cloud fraction
    real(8), intent(out) :: vcf      ! corresponding variance on ecf
    integer, intent(out) :: ngood    ! number of good estimates
    real(8), intent(in)  :: ptop_mb  ! height (mb) from cloud_height subroutine

    ! Locals:
    integer    :: n,jch
    real(8)    :: H(nco2),F(nco2)

    etop = -1.d0
    vtop = -1.d0
    ecf  = -1.d0
    vcf  = -1.d0
    ngood = 0

    !     Profile not assimilated if data from 2 windows channels bad
    !     and/or if data from 2 reference co2 channels bad    
    if ( cldflag == -1 ) return

    n = 0
    H(:) = 0.d0
    F(:) = 0.d0

    do jch = 1, nco2

      !     Check for zero cloud fraction

      if ( CF(jch) > -0.9d0 .and. CF(jch) < 1.D-6 ) then
        n = n + 1
        H(n) = p0
        F(n) = 0.d0
      else


        !       Consider only valid values of cloud fraction above some threshold
     
        !       Important logic: for values above 1.0 of co2-slicing cloud fraction,
        !       set it to 1.0 and force the top equal to the effective height he.
        !       co2-slicing not allowed to give estimates below he, which happens
        !        for cloud fraction cf > 1.0.

        if ( ht(jch) > 0.0d0 ) then
          n = n + 1
          H(n) = ht(jch)
          F(n) = min(CF(jch), 1.0d0)
          F(n) = max(F(n), 0.d0)
          if ( CF(jch) > 1.0d0 ) H(n) = he
        end if
      end if

    end do


    ngood = n

    !     Compute mean and variance

    if ( n >= 1 ) then

      call calcul_median_fast(n,H,F,etop,ecf)
    
      vtop = sqrt ( sum((H(1:n) - etop)**2) / n )
      vcf  = sqrt ( sum((F(1:n) - ecf)**2) / n )         

      if ( n == 1 ) then
        vtop = 50.d0
        vcf  = 0.20d0
      end if
       
    else

      !    If no solution from co2-slicing, and not predetermined clear, 
      !    assume cloudy with top equal to effective height he;
      !    however if he is very close to surface pressure p0, assume clear.

      etop = he
      ecf  = 1.0d0
      if (cfsub >= 0.05d0) then
        ecf = cfsub
        etop = min( min(he,ptop_mb) , p0 - 50.0d0)
      end if
      vtop = 50.d0
      vcf  = 0.30d0
      if ( he > (p0 - 10.d0) ) ecf = 0.d0
      if ( gncldflag == 0 ) then
        ecf = 0.0d0
        etop = p0
      end if
    end if

    if ( ecf < 0.05d0 ) then
      ecf = 0.0d0
      etop = p0
    end if
  
  end subroutine seltop
  
  !--------------------------------------------------------------------------
  ! calcul_median_fast
  !--------------------------------------------------------------------------
  subroutine calcul_median_fast(nEstimates, Hin, Fin, ctp, cfr)
    !
    ! :Purpose: Compute cloud fraction and height median.
    !
    ! 
    implicit none

    ! Arguments:
    integer, intent(in)  :: nEstimates  ! Number of Co2 slicing estimates
    real(8), intent(in)  :: Hin(:)      ! Array of Height estimates
    real(8), intent(in)  :: Fin(:)      ! Array of Cloud fraction estimates
    real(8), intent(out) :: ctp         ! Median of  cloud top pressure estimates
    real(8), intent(out) :: cfr         ! Corresponding cloud fraction

    ! Locals:
    integer   :: index(nEstimates)
    real(4)   :: H(nEstimates)
    integer   :: i

    if (nEstimates == 1) then
      ctp = Hin(nEstimates)
      cfr = Fin(nEstimates)
    else
      H(1:nEstimates) = Hin(1:nEstimates)
      call IPSORT(index,H,nEstimates)
      if (mod(nEstimates,2) == 0) then ! N - pair
        i = index(nEstimates / 2)
      else                             ! N - impair
        i = index(1 + nEstimates / 2)
      end if
      ctp = Hin(i)
      cfr = Fin(i)
    end if

  end subroutine calcul_median_fast

  !--------------------------------------------------------------------------
  ! min_pres_new
  !--------------------------------------------------------------------------
  subroutine min_pres_new( maxheight, minp, pmin, dt1, p0, tau, plev, cldflag, &
                           modelTopIndex )
    !
    ! :Purpose: from total transmittance array, find minimum height 
    !           level of sensitivity for a number of profiles and channels.
    !           this may be used to select for assimilation only the
    !           observations without sensitivity to clouds, that is the
    !           response function significant only above cloud level.
    !           the criterion is that dtau/dplev > 0.01 for a 100 mb layer.
    !
    !
    implicit none

    ! Arguments:
    real(8), intent(out)  :: maxHeight(:)                    ! Height (hPa) of the maximum of the weighting function
    integer, intent(out)  :: minp(size(maxHeight))           ! vertical level corresponding to pmin
    real(8), intent(out)  :: pmin(size(maxHeight))           ! minimum height of sensitivity (hPa)
    real(8), intent(out)  :: dt1(size(maxHeight))            ! value of 'dtau/dlogp' at model top
    real(8), intent(in)   :: p0                              ! surface pressure (hPa)
    real(8), intent(in)   :: plev(:)                         ! pressure levels (hPa)
    real(8), intent(in)   :: tau(size(maxHeight),size(plev)) ! layer to space transmittances (0.-1.)
    integer, intent(in)   :: cldflag                         ! Cloudy flag (0 Clear, 1 Cloudy, -1 undefined)
    integer, intent(in)   :: modelTopIndex                   ! rt model level nearest to model top

    ! Locals:   
    real(8) :: maxwf
    integer :: levelIndex, channelIndex, ipos(1), nlev, nchn
    real(8),allocatable :: wfunc(:), rap(:)

    minp(:) = -1
    pmin(:) = -1.d0
    dt1(:)  = -1.d0

    if ( cldflag == -1 ) return

    nlev = size(plev)
    nchn = size(maxHeight)

    allocate( wfunc(nlev-1), rap(nlev-1) )

    do levelIndex = 1, nlev - 1
      rap(levelIndex) = log( plev(levelIndex + 1) / plev(levelIndex) )
    end do

    channels: do channelIndex = 1, nchn

      !       Profile not assimilated if data from 2 windows channels bad
      !       and/or if data from 2 reference co2 channels bad
    
      do levelIndex = 1, nlev
        if ( tau(channelIndex,levelIndex) < 0.d0) cycle channels
      end do

      minp(channelIndex) = nlev
      pmin(channelIndex) = min(plev(nlev),p0)

      !   Compute entire array of dtau/dlog(P)
          
      do levelIndex = 1, nlev - 1
        wfunc(levelIndex) = (tau(channelIndex,levelIndex) - tau(channelIndex,levelIndex + 1)) / rap(levelIndex) 
      end do
       
      dt1(channelIndex) = wfunc(modelTopIndex)

      !   If channel sees the surface, don't recalculate minp and pmin

      if ( tau(channelIndex,nlev) > 0.01d0 ) cycle channels

      ! Recherche du maximum
      ipos = maxloc( wfunc(:) )
      ! Calcul de la valeur du maximum
      maxwf = wfunc(ipos(1))
      ! maximum entre les 2 niveaux puisque WF calculee pour une couche finie ( discutable ?)
      maxheight(channelIndex)= 0.5d0 * ( plev(ipos(1)) +  plev(ipos(1) + 1)  )

      !      If channel doesn't see the surface, see where dtau/dlog(plev) becomes important
      !      for recomputation of minp and pmin.

      do levelIndex = nlev - 1, ipos(1), -1
        if ( ( wfunc(levelIndex)/ maxwf ) > 0.01d0) then
          minp(channelIndex) = levelIndex + 1
          pmin(channelIndex) = min(plev(levelIndex + 1),p0)
          exit
        end if
      end do
     
    end do channels

    deallocate(wfunc, rap )

  end subroutine min_pres_new

  !--------------------------------------------------------------------------
  ! cloud_height
  !--------------------------------------------------------------------------
  subroutine cloud_height ( ptop, ntop, btObs, cldflag, tt, height, p0, plev, &
                            ichref, lev_start, iopt )
    !
    ! :Purpose: 
    !         Computation of cloud top height (above the ground)
    !         based on matching observed brightness temperature at a 
    !         reference surface channel with background temperature profile.
    !         to use with one reference channel. used here on model levels.
    !
    implicit none

    ! Arguments:
    real(8), intent(out)   :: ptop            ! Chosen equivalent cloud tops  (in hpa|m with iopt = 1|2)
    integer, intent(out)   :: ntop            ! Number of possible ptop solutions
    real(8), intent(in)    :: btObs(:)        ! Observed brightness temperature (deg k)
    integer, intent(in)    :: cldflag         ! Cloudy flag (0 Clear, 1 Cloudy, -1 undefined)
    real(8), intent(in)    :: tt(:)           ! Temperature profiles (deg K)
    real(8), intent(in)    :: height(size(tt))! Height profiles above ground (m)
    real(8), intent(in)    :: p0              ! Surface pressure (hPa)
    real(8), intent(in)    :: plev(size(tt))  ! Pressure levels (hPa)
    integer, intent(in)    :: ichref          ! Chosen reference surface channel
    integer, intent(inout) :: lev_start       ! Level to start iteration (ideally tropopause)
    integer, intent(in)    :: iopt            ! Levels using plev (1) or height (2)

    ! Locals:
    integer :: itop
    integer :: nht, nlev
    real(8), allocatable :: ht(:)

    nlev = size(tt)
    allocate(ht(nlev))

    if ( iopt == 1 ) then
     
      ptop = p0
      ntop = 1      

      if ( cldflag == -1 ) return
      
      call get_top ( ht, nht, btObs(ichref), tt, plev, lev_start, iopt ) 

      itop = 1
      if ( nht >= 2 ) itop = 2
      ptop = min ( ht(itop), p0 )
      ntop = nht

    else if ( iopt == 2 ) then
      
      ptop = 0.d0
      ntop = 1      

      if ( cldflag == -1 ) return

      call get_top ( ht,nht, btObs(ichref),tt,height,lev_start,iopt )

      itop = 1
      if ( nht >= 2 ) itop = 2
      ptop = max ( ht(itop), 0.d0 )
      ntop = nht
       
    end if

    deallocate(ht)

  end subroutine cloud_height

  !--------------------------------------------------------------------------
  ! garand1998nadon 
  !--------------------------------------------------------------------------
  subroutine garand1998nadon ( cldflag, btObs, tg, tt, height, &
                               ptop_eq, ntop_eq, ichref )
    !
    ! :Purpose: Determine if the profiles are clear or cloudy based on
    !           the algorithm of garand & nadon 98 j.clim v11 pp.1976-1996
    !           with channel iref.
    implicit none

    ! Arguments:
    integer, intent(inout) :: cldflag         ! Cloudy flag (0 Clear, 1 Cloudy, -1 undefined)
    real(8), intent(in)    :: btObs(:)        ! Observed brightness temperatures (K)
    real(8), intent(in)    :: tg              ! Guess skin temperatures (K)
    real(8), intent(in)    :: tt(:)           ! Guess temperature profiles (K)
    real(8), intent(in)    :: height(size(tt))! Guess height profile above ground (m)
    real(8), intent(in)    :: ptop_eq         ! Chosen equivalent cloud tops (m)
    integer, intent(in)    :: ntop_eq         ! Number of possible ptop_eq solutions
    integer, intent(in)    :: ichref          ! Chosen reference surface channel

    ! Locals:
    integer    :: ninv
    real(8)    :: lev(2)
      
    lev(1) = 222.d0
    lev(2) = 428.d0

    if ( cldflag == -1 ) return

    if ( btObs(ichref) >= tg - 3.d0 .and. btObs(ichref) <= tg + 3.d0 ) then
      cldflag = 0
      return
    end if

    if ( btObs(ichref) >= tg - 4.d0 .and. btObs(ichref) <= tg - 3.d0 ) then
      if ( ptop_eq > 1100.d0 ) then
        cldflag = 1
        return
      else
        cldflag = 0
        return
      end if
    end if
    
    if ( ptop_eq > 728.d0 ) then
      cldflag = 1
      return
    end if

    if ( tg - btObs(ichref) > 8.d0 ) then 
      if ( ntop_eq >= 3 ) then
        if ( ptop_eq > 73.d0 ) then
          cldflag=1
          return
        else
          cldflag=0
          return
        end if
      else
        call monotonic_inversion (ninv, tg,tt,height,lev(1))
        if ( ninv == 1 ) then
          if ( ptop_eq > 222.d0 ) then
            cldflag = 1
            return
          else
            cldflag = 0 
            return
          end if
        else
          cldflag = 0
          return
        end if
      end if
    end if
    
    if ( tg - btObs(ichref) > 5.d0 ) then
      if ( ntop_eq >= 3 ) then
        if ( ptop_eq > 222.d0 ) then
          cldflag = 1
          return
        else
          cldflag = 0
          return
        end if
      else
        call monotonic_inversion (ninv, tg,tt,height,lev(2))
        if ( ninv == 1) then
          if( ptop_eq > 428.d0 ) then
            cldflag = 1
            return
          else
            cldflag = 0
            return
          end if
        else
          cldflag = 0
        end if
      end if
    else
      cldflag = 0
    end if
    
  end subroutine garand1998nadon

  !--------------------------------------------------------------------------
  ! monotonic_inversion
  !--------------------------------------------------------------------------
  subroutine monotonic_inversion ( ninvr, tg, tt, height, lvl )
    !
    ! :Purpose: Determine if there is a presence (ninvr=1) or not (ninvr=0)
    !           of a temperature inversion going from the surface up to the
    !           height lvl
    !
    ! :Arguments:
    !            :ninvr: PRESENCE (1) OR NOT (0) OF A TEMPERATURE INVERSION
    !                    FROM THE SURFACE TO HEIGHT LVL
    !
    implicit none

    ! Arguments:
    real(8), intent(in)  :: tt(:)           ! Temperature profile (K)
    real(8), intent(in)  :: height(size(tt))! Height profile above ground (m)
    real(8), intent(in)  :: tg              ! Skin temperature (K)
    real(8), intent(in)  :: lvl             ! Height to search for temperature inversion (m)
    integer, intent(out) :: ninvr           ! Number of inversions

    ! Locals:
    integer   :: levelIndex, nlevels

    ninvr = 0
    nlevels = size ( tt )
    if ( tg - tt(nlevels) < 0.d0 ) then
      ninvr = 1
      do levelIndex = nlevels - 1, 1, -1
        if ( height(levelIndex) > lvl ) exit
        if ( tt(levelIndex+1) - tt(levelIndex) > 0.d0 ) then
          ninvr = 0
          exit
        end if
      end do
    end if

  end subroutine monotonic_inversion


  !--------------------------------------------------------------------------
  ! estim_ts
  !--------------------------------------------------------------------------
  subroutine estim_ts( ts, tg, emi, rcal, radobs, sfctau, cldflag, &
                       ichref, myCoefs )
    !
    ! :Purpose: Get an estimated skin temperature by inversion of
    !           radiative transfer equation assuming guess t and q profiles
    !           are perfect. designed for a single channel ichref and nprf
    !           profiles. assumes a real tg (guess) over oceans and a tg 
    !           with hypothesis of unity emissivity over land.
    !      
    ! :Note:  Uses rcal = B(TG)*EMI*SFCtau + ATMOS_PART
    !         ts = B(ts)*EMI*SFCtau + ATMOS_PART
    !         SOLVES FOR ts
    !
    implicit none

    ! Arguments:
    integer, intent(in)           :: ichref            ! Reference surface channel (subset values)
    integer, intent(in)           :: cldflag           ! Cloudy flag (0 Clear, 1 Cloudy, -1 undefined)
    real(8), intent(in)           :: tg                ! Guess skin temperature (K)
    real(8), intent(in)           :: emi(:)            ! Surface emissivities from window channel (0.-1.)
    real(8), intent(in)           :: rcal(size(emi))   ! Computed clear radiances (mW/m2/sr/cm-1)
    real(8), intent(in)           :: radobs(size(emi)) ! Observed radiances (mW/m2/sr/cm-1)
    real(8), intent(in)           :: sfctau(size(emi)) ! Surface to space transmittances (0.-1.)
    real(8), intent(out)          :: ts                ! Retrieved skin temperature (-1. for missing)
    type(rttov_coefs), intent(in) :: myCoefs           ! RTTOV coefficients structure

    ! Locals:
    real(8)    :: rtg,radtg
    real(8)    :: radts,tstore,t_effective
  
    ts = -1.d0

    if ( cldflag /= 0 ) return
    if ( ichref == -1 ) return


    !   Transform guess skin temperature to plank radiances 

    t_effective =  myCoefs % coef % ff_bco(ichref) + myCoefs % coef % ff_bcs(ichref) * TG

    radtg =  myCoefs % coef % planck1(ichref) / &
         ( exp( myCoefs % coef % planck2(ichref) / t_effective ) - 1.0d0 )


    !  Compute TOA planck radiances due to guess skin planck radiances

    rtg = radtg * EMI(ichref) * SFCtau(ichref)


    !  Compute true skin planck radiances due to TOA true planck radiances
    
    radts = ( RADOBS(ichref) + rtg - rcal(ichref) ) / &
         ( EMI(ichref) * SFCtau(ichref) )

    if (radts <= 0.d0) then
      write(*,'(A25,1x,8e14.6)') "Warning, negative radts", RADOBS(ichref), rtg, rcal(ichref), EMI(ichref), SFCtau(ichref), &
           ( RADOBS(ichref) + rtg - rcal(ichref) ), ( EMI(ichref) * SFCtau(ichref) ), radts
      write(*,*) "Skipping tskin retrieval."
      return
    end if

    
    !  Transform true skin planck radiances to true skin temperatures

    tstore = myCoefs % coef % planck2(ichref) / log( 1.0d0 + myCoefs % coef % planck1(ichref) / radts )

    ts = ( tstore - myCoefs % coef % ff_bco(ichref) ) / myCoefs % coef % ff_bcs(ichref)
    

  end subroutine estim_ts

  !--------------------------------------------------------------------------
  ! cloud_top
  !--------------------------------------------------------------------------
  subroutine cloud_top ( ptop_bt, ptop_rd, ntop_bt, ntop_rd, btObs,  &
       tt, height, rcal, p0, radObs, cloudyRadiance, plev, &
       cldflag, lev_start, iopt, ihgt, &
       ilist,rejflag_opt, ichref_opt )
    !
    ! :Purpose: Computation of cloud top height (above the ground)
    !           based on matching observed brightness temperature with 
    !           background temperature profiles and/or computed observed
    !           radiances with background radiance profiles.
    !           to use with more than one channel. used here on rttov levels.
    !
    !
    implicit none

    ! Arguments:
    real(8), intent(out)         :: ptop_bt(:)        ! Chosen equivalent cloud tops based on brightness temperatures (in hpa|m with iopt = 1|2)
    real(8), intent(out)         :: ptop_rd(:)        ! Chosen equivalent cloud tops based on radiances (in hpa|m with iopt = 1|2)
    integer,           intent(out)   :: ntop_bt(:)        ! Number of possible ptop_bt solutions
    integer,           intent(out)   :: ntop_rd(:)        ! Number of possible ptop_rd solutions
    real(8),           intent(in)    :: btObs(:)          ! Observed brightness temperautres (K)
    real(8),           intent(in)    :: tt(:)             ! Temperature profiles (K)
    real(8),           intent(in)    :: height(:)         ! Height profiles above ground (m)
    real(8),           intent(in)    :: rcal(:)           ! Computed clear radiances (mW/m2/sr/cm-1)
    real(8),           intent(in)    :: p0                ! Surface pressure (hPa)
    real(8),           intent(in)    :: radObs(:)           ! Computed observed radiances (mW/m2/sr/cm-1)
    real(8),           intent(in)    :: cloudyRadiance(:,:) ! Computed cloud radiances from each level (hPa)
    real(8),           intent(in)    :: plev(:)           ! Pressure levels (hPa)
    integer,           intent(in)    :: cldflag           ! Cloudy flag (0 Clear, 1 Cloudy, -1 undefined)
    integer,           intent(inout) :: lev_start         ! Level to start iteration (ideally tropopause)
    integer,           intent(in)    :: iopt              ! Levels using plev (1) or height (2)
    integer,           intent(in)    :: ihgt              ! Get *_bt* only (0), *_rd* only (1), both (2)
    integer,           intent(in)    :: ilist(:)          ! List of the channel numbers (subset values)
    integer, optional, intent(in)    :: rejflag_opt(1:,0:)! Flags for rejected observations
    integer, optional, intent(in)    :: ichref_opt        ! Reference surface channel (subset value)

    ! Locals:
    integer             :: jch,jc,itop,nht,i10,i,nlev,nch
    real(8),allocatable :: ht(:)
    logical             :: clear, cloudy

    !    Profile not assimilated if data from 2 windows channels bad

    ptop_bt(:) = -10.d0
    ptop_rd(:) = -10.d0
    ntop_bt(:) = 0.d0
    ntop_rd(:) = 0.d0

    if ( cldflag == -1 ) return

    nlev = size( tt )
    i10=1
    do i=2, nlev
      if (plev(i - 1) <= 100.d0 .and. plev(i) > 100.d0) then
        i10 = i
        exit
      end if
    end do

    !    predetermined clear
    if ( present(rejflag_opt)) then
      clear = ( sum( rejflag_opt(ichref_opt,:) ) == 0 )
    else
      clear = ( cldflag == 0 )
    end if

    if ( clear ) then
      
      if ( iopt == 1 ) then
        ptop_bt(:) = min ( plev(nlev), p0 )
        ptop_rd(:) = min ( plev(nlev), p0 )
      else if ( iopt == 2 ) then
        ptop_bt(:) = 0.d0
        ptop_rd(:) = 0.d0
      end if
     
      ntop_bt(:) = 1
      ntop_rd(:) = 1
     
      lev_start = max ( lev_start , i10 )

      return

    end if

    allocate ( ht(nlev) )
    nch = size( ilist)

    channels: do jch = 1, nch
       
      jc = ilist(jch)
       
      !    missing channel ... yes it can happen
      if (jc == -1) cycle channels
      !     gross check failure
      if ( present(rejflag_opt) )  then
        if ( rejflag_opt(jc,9) == 1 ) cycle channels
      else
        if ( btObs(jc) < 150.d0 .or. btObs(jc) > 350.d0) cycle channels
      end if
      !      no clouds if observed radiance warmer than clear estimate

      if ( radObs(jc) > rcal(jc) ) then

        if ( iopt == 1 ) then
          ptop_bt(jc) = min ( plev(nlev), p0 )
          ptop_rd(jc) = min ( plev(nlev), p0 )
        else if ( iopt == 2 ) then
          ptop_bt(jc) = 0.d0
          ptop_rd(jc) = 0.d0
        end if
      
        ntop_bt(jc) = 1
        ntop_rd(jc) = 1

        cycle channels
           
      end if

      !    cloudy

      if ( present(rejflag_opt)) then
        cloudy = ( rejflag_opt(jc,11) == 1 .and. rejflag_opt(jc,23) == 1 )
      else
        cloudy = ( cldflag == 1 )
      end if
      
      if ( cloudy ) then

        if ( iopt == 1 ) then
          
          if ( ihgt == 0 .or. ihgt == 2 ) then
            call get_top ( ht,nht, btObs(jc), tt, plev, lev_start, iopt ) 
            itop = 1
            if ( nht >= 2 ) itop = 2
            ptop_bt(jc) = min ( ht(itop), p0 )
            ntop_bt(jc) = nht
          end if
          
          if ( ihgt == 1 .or. ihgt == 2 ) then
            call get_top ( ht, nht, radObs(jc), cloudyRadiance(jc,:), plev, lev_start, iopt )
            itop = 1
            if ( nht >= 2 ) itop = 2
            ptop_rd(jc) = min ( ht(itop), p0 )
            ntop_rd(jc) = nht
          end if
          
        else if ( iopt == 2 ) then 
          
          if ( ihgt == 0 .or. ihgt == 2 ) then
            call get_top ( ht,nht, btObs(jc),tt,height,lev_start,iopt) 
            itop = 1
            if ( nht >= 2 ) itop = 2
            ptop_bt(jc) = max ( ht(itop), 0.d0 )
            ntop_bt(jc) = nht
          end if
          
          if ( ihgt == 1 .or. ihgt == 2 ) then
            call get_top ( ht,nht, radObs(jc),cloudyRadiance(jc,:),height,lev_start,iopt)
            itop = 1
            if ( nht >= 2 ) itop = 2
            ptop_rd(jc) = max ( ht(itop), 0.d0 )
            ntop_rd(jc) = nht
          end if
          
        end if
       
      end if
      
    end do channels

    deallocate ( ht )

  end subroutine cloud_top

  !--------------------------------------------------------------------------
  ! get_top
  !--------------------------------------------------------------------------
  subroutine get_top (ht, nht, bt, tt ,pp, lev_start, iopt)
    !
    ! :Purpose: Computation of cloud top height and number of possible heights.
    !
    implicit none

    ! Arguments:
    real(8), intent(out)   :: ht(:)            ! Cloud top height in hpa or meters (iopt = 1 or 2)
    integer, intent(out)   :: nht              ! Number of possible cloud height solutions 
    real(8), intent(in)    :: bt               ! Observed brightness temperatures (deg k) or radiance (mw/m2/sr/cm-1)
    real(8), intent(in)    :: tt(:)            ! Temperature profile (deg k) or computed cloud radiance from each level to top
    real(8), intent(in)    :: pp(size(tt))     ! Pressure (hpa) or heights (m) profile (iopt=1 or 2)
    integer, intent(inout) :: lev_start        ! Level to start iteration (ideally tropopause, if <= 0, search & start at coldest level)
    integer, intent(in)    :: iopt             ! Height units in hpa (1) or in meters (2)

    ! Locals:
    integer             :: i, im(1), i10, nlev
    real(8),allocatable :: logp(:)
    real(8)             :: dt, a, b

    ht(:) = -1.
    
    im = lev_start

    nlev = size( tt )

    if ( lev_start <= 0 ) then

      !    Search index im(1) where tt is minimum
      im = minloc ( tt )

      i10 = -1
      do i = 2, nlev
        if (pp(i-1) <= 100.d0 .and. pp(i) > 100.d0) then
          I10 = I
          exit
        end if
      end do
       
      lev_start = im(1)
     
      if ( im(1) == nlev ) then
        lev_start = max(lev_start,i10)
        nht = 1
        ht(1) = pp(nlev)
        return
      end if
       
    end if

    if (iopt == 1) then
      allocate ( logp(nlev) )
      logp(:) = log(pp(:))
    end if

    nht = 0        
    
    do I = im(1), nlev - 1
      dt = tt(I + 1) - tt(I) + 1.D-12
      if ( bt > tt(I) .and. bt <= tt(I + 1) ) then
        
        nht = nht + 1
        
        if (iopt == 1) then
          a = logp(I) + (logp(I + 1) - logp(I)) / dt * ( bt - tt(I))
          ht(nht) = exp(a)
        end if

        if (iopt == 2) then
          b  = pp(I) + (pp(I+1) - pp(I)) / dt * (bt - tt(I))
          ht(nht) = b
        end if
         
      else if ( bt >= tt(I+1) .and. bt < tt(I) ) then
      
        nht = nht + 1
      
        if (iopt == 1) then
          a  = logp(I + 1)- (logp(I + 1)-logp(I)) / dt * (tt(I + 1) - bt)
          ht(nht) = exp(A)
        end if

        if(iopt == 2) then
          b = pp(I + 1)- (pp(I + 1) - pp(I)) / dt * (tt(I + 1) - bt)
          ht(nht) = b
        end if
       
      end if
    end do
    
    
    if ( nht == 0 .and. bt < tt(im(1)) )  then
      nht  = 1
      ht(1) = pp(im(1))
    else if ( nht == 0 .and. bt > tt(nlev) )  then
      nht   = 1
      ht(1) = pp(nlev)
    end if

    if (iopt==1)  deallocate ( logp )
    
  end subroutine get_top

  !--------------------------------------------------------------------------
  ! get_avhrr_emiss
  !--------------------------------------------------------------------------
  subroutine get_avhrr_emiss( iasi_surfem1, freqiasi, nchaniasi, avhrr_surfem1 )
    ! 
    ! :Purpose: choisi l'emissivite d'un canal IASI proche pour AVHRR
    !           a raffiner pour prendre en  compte la largeur  des canaux AVHRR ??
    !
    implicit none

    ! Arguments:
    real(8), intent(in)  :: iasi_surfem1(nchaniasi)! IASI emissivities
    real(8), intent(in)  :: freqiasi(nchaniasi)    ! IASI wavenumbers (cm-1)
    integer, intent(in)  :: nchaniasi              ! Number of IASI channels
    real(8), intent(out) :: avhrr_surfem1(NIR)     ! AVHRR emissivities

    ! Locals:
    real(8),parameter :: freqavhrr(NIR)= (/0.2687000000D+04 , 0.9272000000D+03 , 0.8377000000D+03/)
    integer           :: indxavhrr(NIR)
    integer           :: i, pos(1)

    do I=1,NIR
      pos = minloc ( abs (freqiasi(:) - freqavhrr(I)) )
      indxavhrr(i) = pos(1)
    end do

    do I=1,NIR
      avhrr_surfem1(i) = iasi_surfem1(indxavhrr(i))
    end do
  
  end subroutine get_avhrr_emiss

  !--------------------------------------------------------------------------
  ! tovs_rttov_avhrr_for_IASI
  !--------------------------------------------------------------------------
  subroutine tovs_rttov_avhrr_for_IASI (headerIndex, surfem1_avhrr, idiasi)
    !
    ! :Purpose: Computation of forward radiance with rttov_direct
    !           (for AVHRR).
    !           appel de RTTOV pour le calcul des radiances AVHRR
    !           (non assimilees mais necessaires au background check IASI)
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: headerIndex      ! Location of IASI observation in TOVS structures and obSpaceData
    real(8), intent(in) :: surfem1_avhrr(3) ! AHVRR surface emissivities
    integer, intent(in) :: idiasi           ! iasi (in fact METOP) number

    ! Locals:
    type (rttov_profile), pointer :: profiles(:)
    type (rttov_chanprof)  :: chanprof(3)
    logical :: calcemis  (3)
    integer ::  list_sensor (3),errorstatus,allocStatus
    integer, save :: idiasi_old=-1
    integer :: ich
    integer :: ichan_avhrr (NIR)
    type (rttov_transmission)  :: transmission
    type (rttov_radiance)      :: radiancedata_d
    type (rttov_emissivity)    :: emissivity(3)
    integer :: nchannels
    integer :: nlevels, iptobs(1)

    if (idiasi_old /= idiasi) then
      list_sensor(1) = 10
      list_sensor(2) = idiasi
      list_sensor(3) = 5
      do ich=1,nir
        ichan_avhrr(ich)=ich
      end do
    
      errorstatus = 0

      if (idiasi_old > 0) then
        call rttov_dealloc_coefs(errorstatus, coefs_avhrr )
        if ( errorstatus /= 0) then
          write(*,*) "Probleme dans rttov_dealloc_coefs !"
          call utl_abort("tovs_rttov_avhrr_for_IASI")
        end if
      end if

      call rttov_read_coefs ( errorstatus, &! out
           coefs_avhrr,                    &! out
           tvs_opts(1),                    &! in
           channels=ichan_avhrr,           &! in
           instrument=list_sensor )         ! in
       
      if ( errorstatus /= 0) then
        write(*,*) "Probleme dans rttov_read_coefs !"
        call utl_abort("tovs_rttov_avhrr_for_IASI")
      end if
     
      idiasi_old = idiasi
   
    end if

    call tvs_getProfile(profiles, 'nl')

    iptobs(1) = headerIndex
    nlevels =  profiles(headerIndex)% nlevels

    nchannels = NIR

    calcemis(:) = .false.
    emissivity(1:3)%emis_in = surfem1_avhrr(1:3)
    ! Build the list of channels/profiles indices

    do  ich = 1, nchannels
      chanprof(ich) % prof = 1
      chanprof(ich) % chan = ich
    end do

    allocStatus = 0
    call rttov_alloc_direct(         &
         allocStatus,                &
         asw=1,                      &
         nprofiles=1,                & ! (not used)
         nchanprof=nchannels,        &
         nlevels=nlevels,            &
         opts=tvs_opts(1),           &
         coefs=coefs_avhrr,          &
         transmission=transmission,  &
         radiance=radiancedata_d,    &
         init=.true.)

    if (allocStatus /= 0) then
      write(*,*) "Memory allocation error"
      call utl_abort('tovs_rttov_avhrr_for_IASI')
    end if

    
    call rttov_direct(            &
         errorstatus,             & ! out
         chanprof,                & ! in
         tvs_opts(1),             & ! in
         profiles(iptobs(:)),     & ! in
         coefs_avhrr,             & ! in
         transmission,            & ! inout
         radiancedata_d,          & ! out
         calcemis=calcemis,       & ! in
         emissivity=emissivity)     ! inout
    
    avhrr_bgck(headerIndex)% radclearcalc(NVIS+1:NVIS+NIR) = radiancedata_d % clear(1:NIR)
    avhrr_bgck(headerIndex)% tbclearcalc(NVIS+1:NVIS+NIR)  = radiancedata_d % bt(1:NIR)
    allocate( avhrr_bgck(headerIndex)% radovcalc(nlevels-1,NVIS+1:NVIS+NIR) )
    avhrr_bgck(headerIndex)% radovcalc(1:nlevels-1,NVIS+1:NVIS+NIR) = radiancedata_d % overcast(1:nlevels-1,1:NIR)
    avhrr_bgck(headerIndex)% emiss(NVIS+1:NVIS+NIR) = emissivity(1:NIR)%emis_out
    avhrr_bgck(headerIndex)% transmsurf(NVIS+1:NVIS+NIR) = transmission% tau_total(1:NIR)


    call rttov_alloc_direct(         &
         allocStatus,                &
         asw=0,                      &
         nprofiles=1,                & ! (not used)
         nchanprof=nchannels,        &
         nlevels=nlevels,            &
         opts=tvs_opts(1),           &
         coefs=coefs_avhrr,          &
         transmission=transmission,  &
         radiance=radiancedata_d )

    if (allocStatus /= 0) then
      write(*,*) "Memory deallocation error"
      call utl_abort('tovs_rttov_avhrr_for_IASI')
    end if

    nullify(profiles)
  
  end subroutine tovs_rttov_avhrr_for_IASI

  !--------------------------------------------------------------------------
  ! cor_albedo
  !--------------------------------------------------------------------------
  subroutine cor_albedo(delta, scos)
    !
    ! :Purpose: ce sous-programme calcule un facteur de correction
    !           pour l'albedo a partir du cosinus de l'angle solaire.
    !
    implicit  none

    ! Arguments:
    real(8), intent(in)  :: scos   ! Cosine of solar zenith angle
    real(8), intent(out) :: delta  ! Correction factor

    ! Locals:
    integer  i1, i2, ierr
    real(8)  x1, x2, g1, g2, a, b
    real(8),parameter ::  s(11)=(/00.00d0, 18.19d0, 31.79d0, 41.41d0, 49.46d0, &
                                  56.63d0, 63.26d0, 69.51d0, 75.52d0, 81.37d0, 87.13d0/)
 
    i1  = 12 - ( scos + 0.05d0) * 10.d0 
    i2  = i1 + 1 
    i1  = min(i1,11)
    i2  = min(i2,11)
    x1  = cos ( s(I1) * MPC_RADIANS_PER_DEGREE_R8 )  
    x2  = cos ( s(I2) * MPC_RADIANS_PER_DEGREE_R8 ) 
    g1  = drcld(i1)
    g2  = drcld(i2)
    if (i1 == i2) then
      delta =g1
    else
      call  lineq ( x1, x2, g1, g2, a, b, ierr )
      delta = a * scos + b
    end if
  
  end subroutine cor_albedo

  !--------------------------------------------------------------------------
  ! drcld
  !--------------------------------------------------------------------------
  real(8) function drcld(iz) 
    !
    ! :Purpose: Generaliser pour toutes les plateformes satellitaires.
    !           Ce sous-programme calcule la normalisation due
    !           a l'angle zenith solaire selon 
    !           MINNIS-HARRISSON (COURBE FIG 7), P1038,JCAM 84.  
    !
    ! :Output: facteur de normalisation
    !
    implicit  none

    ! Arguments:
    integer, intent(in) ::  iz ! Index for Sun angle bin

    ! Locals:
    real(8),parameter ::  drf(11)=(/1.000d0, 1.002d0, 1.042d0, 1.092d0, 1.178d0, 1.286d0, &
                                    1.420d0, 1.546d0, 1.710d0, 1.870d0, 2.050d0/) 

    drcld = drf (iz)
    
  end function drcld


  !--------------------------------------------------------------------------
  ! visocn 
  !--------------------------------------------------------------------------
  subroutine visocn(sz, satz, rz, anisot, zlamb, zcloud, ierr)
    !
    ! :Purpose: This routine provides the corrective factors for the anisotropy
    !           of reflectance over clear ocean.
    !                 
    !
    ! :Notes:  Obtained from dr pat minnis,langley , and based on the work
    !          of minnis and harrisson,jcam 1984,p993.
    !          the routine is a look up table along with interpolation on the 
    !          three angles. 
    !
    implicit  none

    ! Arguments:
    real(8), intent(in)  :: sz     ! sun zenith angle in degrees (0 to 90)
    real(8), intent(in)  :: satz   ! satellite zenith angle (0 to 90)
    real(8), intent(in)  :: rz     ! relative angle in degrees (0-180) with 0 as backscattering and 180 as forward scattering
    real(8), intent(out) :: anisot ! anisotropic corrective factor (khi in minnis-harrisson)
    real(8), intent(out) :: zlamb  ! corrective factor for lambertian reflectance (Ocean surface)
    real(8), intent(out) :: zcloud ! Same as zlamb but for cloud surface
    integer, intent(out) :: ierr   ! Error code (0=ok; -1=problem with interpolation)

    ! Locals:
    integer  i1, i2, j1, j2, k1, k2, l, i, n, m, j, k
    real(8) cc, d1, d2, slope, intercept, x1, x2
    real(8) g1, g2, da(2), dd(2) 
    real(8), parameter :: s(11)=(/0.0d0,18.19d0,31.79d0,41.41d0,49.46d0,56.63d0, &
         63.26d0,69.51d0,75.52d0,81.37d0,87.13d0/)    
    real(8), parameter :: r(13)=(/0.0d0, 15.0d0, 30.0d0, 45.0d0, 60.0d0, 75.0d0, 90.0d0, &
         105.0d0, 120.0d0, 135.0d0, 150.0d0, 165.0d0, 180.0d0/)
    real(8), parameter :: v(10)=(/0.0d0, 10.0d0, 20.0d0, 30.0d0, 40.0d0, 50.0d0, 60.0d0, &
         70.0d0, 80.0d0, 90.0d0/)
    real(8) vnorm(11,10,13)

    data ((vnorm(1,j,k),j=1,10),k=1,13)/  &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0/ 
    
    data ((vnorm(2,j,k),j=1,10),k=1,13)/  &
         1.154d0, .960d0, .896d0, .818d0, .748d0, .825d0, .922d0,1.018d0,1.179d0,1.334d0, &
         1.154d0, .954d0, .838d0, .799d0, .735d0, .786d0, .883d0, .960d0,1.128d0,1.250d0, &
         1.514d0, .973d0, .825d0, .786d0, .722d0, .754d0,0.838d0,0.922d0,1.063d0,1.160d0, &
         1.514d0,0.967d0,0.864d0,0.818d0,0.715d0,0.728d0,0.793d0,0.876d0,1.005d0,1.102d0, &
         1.514d0,0.967d0,0.896d0,0.889d0,0.702d0,0.696d0,0.773d0,0.851d0,0.954d0,1.038d0, &
         1.514d0,1.070d0,0.986d0,0.922d0,0.677d0,0.696d0,0.754d0,0.838d0,0.922d0,1.012d0, &
         1.514d0,1.270d0,0.967d0,0.870d0,0.677d0,0.664d0,0.709d0,0.773d0,0.857d0,0.954d0, &
         1.514d0,1.495d0,1.166d0,0.960d0,0.683d0,0.690d0,0.728d0,0.806d0,0.896d0,0.999d0, &
         1.514d0,1.959d0,1.534d0,1.025d0,0.973d0,0.709d0,0.754d0,0.857d0,0.954d0,1.050d0, &
         1.514d0,2.165d0,2.165d0,1.270d0,1.038d0,0.760d0,0.812d0,0.902d0,1.012d0,1.115d0, &
         1.514d0,2.275d0,2.262d0,1.688d0,1.115d0,0.780d0,0.857d0,0.954d0,1.070d0,1.173d0, &
         1.514d0,2.326d0,2.520d0,2.172d0,1.257d0,0.812d0,0.883d0,1.005d0,1.108d0,1.212d0, &
         1.514d0,2.359d0,2.951d0,2.255d0,1.411d0,0.980d0,0.915d0,1.050d0,1.160d0,1.295d0/ 

    data ((vnorm(3,j,k),j=1,10),k=1,13)/   &
         0.897d0,0.792d0,0.765d0,0.765d0,0.778d0,0.897d0,0.996d0,1.095d0,1.306d0,1.431d0, &
         0.897d0,0.712d0,0.739d0,0.745d0,0.765d0,0.891d0,0.970d0,1.069d0,1.214d0,1.359d0, &
         0.897d0,0.666d0,0.699d0,0.745d0,0.759d0,0.811d0,0.917d0,1.042d0,1.148d0,1.306d0, &
         0.897d0,0.646d0,0.693d0,0.739d0,0.693d0,0.752d0,0.858d0,0.989d0,1.102d0,1.234d0, &
         0.897d0,0.686d0,0.679d0,0.726d0,0.679d0,0.693d0,0.792d0,0.924d0,1.049d0,1.154d0, &
         0.897d0,0.660d0,0.673d0,0.693d0,0.646d0,0.660d0,0.759d0,0.858d0,1.003d0,1.102d0, &
         0.897d0,0.673d0,0.765d0,0.792d0,0.712d0,0.600d0,0.699d0,0.811d0,0.963d0,1.055d0, &
         0.897d0,0.706d0,0.772d0,0.917d0,0.904d0,0.613d0,0.726d0,0.858d0,1.055d0,1.121d0, &
         0.897d0,0.825d0,0.924d0,0.996d0,0.989d0,0.686d0,0.778d0,0.937d0,1.115d0,1.181d0, &
         0.897d0,1.036d0,1.253d0,1.286d0,1.260d0,0.778d0,0.858d0,0.996d0,1.181d0,1.260d0, &
         0.897d0,1.201d0,1.788d0,1.986d0,1.827d0,0.884d0,0.851d0,1.062d0,1.227d0,1.333d0, &
         0.897d0,1.530d0,2.249d0,2.546d0,2.381d0,1.352d0,0.891d0,1.108d0,1.286d0,1.405d0, &
         0.897d0,1.854d0,2.401d0,3.325d0,2.559d0,1.590d0,0.937d0,1.168d0,1.214d0,1.425d0/ 

    data ((vnorm(4,j,k),j=1,10),k=1,13)/  &
         0.752d0,0.800d0,0.745d0,0.717d0,0.759d0,0.891d0,1.149d0,1.309d0,1.469d0,1.650d0, &
         0.752d0,0.773d0,0.717d0,0.703d0,0.752d0,0.835d0,1.065d0,1.246d0,1.406d0,1.552d0, &
         0.752d0,0.731d0,0.689d0,0.703d0,0.745d0,0.814d0,0.988d0,1.176d0,1.323d0,1.476d0, &
         0.752d0,0.689d0,0.675d0,0.654d0,0.696d0,0.752d0,0.940d0,1.100d0,1.246d0,1.378d0, &
         0.752d0,0.675d0,0.661d0,0.633d0,0.668d0,0.717d0,0.877d0,1.030d0,1.176d0,1.309d0, &
         0.752d0,0.647d0,0.640d0,0.620d0,0.613d0,0.682d0,0.814d0,0.947d0,1.107d0,1.232d0, &
         0.752d0,0.633d0,0.620d0,0.613d0,0.606d0,0.640d0,0.773d0,0.898d0,1.044d0,1.162d0, &
         0.752d0,0.626d0,0.626d0,0.626d0,0.620d0,0.654d0,0.821d0,0.947d0,1.128d0,1.225d0, &
         0.752d0,0.633d0,0.633d0,0.633d0,0.647d0,0.675d0,0.877d0,1.009d0,1.183d0,1.274d0, &
         0.752d0,0.682d0,0.717d0,0.961d0,1.023d0,0.968d0,0.940d0,1.142d0,1.274d0,1.413d0, &
         0.752d0,0.856d0,1.037d0,1.434d0,1.594d0,1.441d0,1.044d0,1.225d0,1.323d0,1.545d0, &
         0.752d0,1.044d0,1.295d0,2.207d0,1.610d0,2.311d0,1.385d0,1.274d0,1.441d0,1.636d0, &
         0.752d0,1.079d0,1.524d0,2.541d0,3.564d0,3.014d0,1.942d0,1.462d0,1.552d0,1.726d0/ 

    data ((vnorm(5,j,k),j=1,10),k=1,13)/  &
         0.552d0,0.588d0,0.617d0,0.638d0,0.724d0,0.860d0,1.133d0,1.362d0,1.556d0,1.678d0, &
         0.552d0,0.581d0,0.602d0,0.617d0,0.652d0,0.803d0,1.075d0,1.326d0,1.484d0,1.592d0, &
         0.552d0,0.559d0,0.588d0,0.595d0,0.617d0,0.731d0,1.018d0,1.283d0,1.412d0,1.527d0, &
         0.552d0,0.531d0,0.538d0,0.574d0,0.595d0,0.710d0,0.946d0,1.240d0,1.341d0,1.463d0, &
         0.552d0,0.516d0,0.523d0,0.552d0,0.559d0,0.695d0,0.911d0,1.226d0,1.291d0,1.412d0, &
         0.552d0,0.516d0,0.523d0,0.538d0,0.538d0,0.652d0,0.882d0,1.154d0,1.240d0,1.348d0, &
         0.552d0,0.516d0,0.523d0,0.538d0,0.523d0,0.595d0,0.774d0,1.075d0,1.169d0,1.269d0, &
         0.552d0,0.531d0,0.545d0,0.552d0,0.566d0,0.609d0,0.817d0,1.140d0,1.248d0,1.369d0, &
         0.552d0,0.538d0,0.545d0,0.566d0,0.581d0,0.645d0,0.911d0,1.240d0,1.319d0,1.441d0, &
         0.552d0,0.566d0,0.552d0,0.574d0,0.710d0,0.839d0,0.982d0,1.298d0,1.391d0,2.323d0, &
         0.552d0,0.566d0,0.559d0,0.710d0,1.147d0,1.176d0,1.040d0,1.348d0,1.671d0,2.674d0, &
         0.552d0,0.588d0,1.133d0,1.355d0,2.194d0,2.803d0,2.201d0,2.459d0,2.904d0,3.126d0, &
         0.552d0,0.710d0,1.341d0,1.757d0,3.026d0,3.900d0,4.445d0,4.503d0,4.445d0,4.503d0/ 

    data ((vnorm(6,j,k),j=1,10),k=1,13)/  &
         0.551d0,0.627d0,0.665d0,0.734d0,0.826d0,0.971d0,1.231d0,1.537d0,1.721d0,1.866d0, &
         0.551d0,0.604d0,0.619d0,0.665d0,0.765d0,0.895d0,1.185d0,1.476d0,1.568d0,1.652d0, &
         0.551d0,0.597d0,0.604d0,0.619d0,0.734d0,0.849d0,1.101d0,1.346d0,1.453d0,1.568d0, &
         0.551d0,0.581d0,0.589d0,0.597d0,0.665d0,0.795d0,1.032d0,1.262d0,1.346d0,1.445d0, &
         0.551d0,0.558d0,0.558d0,0.566d0,0.612d0,0.727d0,0.987d0,1.201d0,1.262d0,1.399d0, &
         0.551d0,0.505d0,0.505d0,0.512d0,0.566d0,0.696d0,0.925d0,1.117d0,1.185d0,1.308d0, &
         0.551d0,0.474d0,0.497d0,0.512d0,0.535d0,0.673d0,0.864d0,1.048d0,1.124d0,1.216d0, &
         0.551d0,0.497d0,0.505d0,0.520d0,0.551d0,0.681d0,0.902d0,1.124d0,1.201d0,1.323d0, &
         0.551d0,0.535d0,0.535d0,0.551d0,0.566d0,0.711d0,1.017d0,1.201d0,1.269d0,1.422d0, &
         0.551d0,0.535d0,0.543d0,0.558d0,0.704d0,1.193d0,1.247d0,1.285d0,1.346d0,1.950d0, &
         0.551d0,0.543d0,0.551d0,0.581d0,0.994d0,1.545d0,1.583d0,1.354d0,2.019d0,2.883d0, &
         0.551d0,0.566d0,0.612d0,0.788d0,1.468d0,2.233d0,2.340d0,2.531d0,2.983d0,3.365d0, &
         0.551d0,0.658d0,0.665d0,1.101d0,2.134d0,3.120d0,4.221d0,4.856d0,4.956d0,5.613d0/ 

    data ((vnorm(7,j,k),j=1,10),k=1,13)/  &
         0.545d0,0.606d0,0.683d0,0.744d0,0.798d0,0.990d0,1.228d0,1.704d0,1.850d0,2.049d0, &
         0.545d0,0.576d0,0.583d0,0.714d0,0.783d0,0.952d0,1.144d0,1.573d0,1.758d0,1.888d0, &
         0.545d0,0.560d0,0.568d0,0.629d0,0.744d0,0.875d0,1.105d0,1.504d0,1.642d0,1.788d0, &
         0.545d0,0.553d0,0.560d0,0.599d0,0.629d0,0.791d0,1.028d0,1.420d0,1.527d0,1.696d0, &
         0.545d0,0.545d0,0.553d0,0.599d0,0.606d0,0.714d0,0.990d0,1.335d0,1.451d0,1.581d0, &
         0.545d0,0.530d0,0.537d0,0.568d0,0.583d0,0.683d0,0.890d0,1.243d0,1.351d0,1.489d0, &
         0.545d0,0.491d0,0.499d0,0.507d0,0.576d0,0.622d0,0.791d0,1.182d0,1.282d0,1.389d0, &
         0.545d0,0.507d0,0.514d0,0.507d0,0.576d0,0.675d0,0.890d0,1.197d0,1.328d0,1.451d0, &
         0.545d0,0.522d0,0.537d0,0.522d0,0.591d0,0.760d0,0.944d0,1.259d0,1.389d0,1.527d0, &
         0.545d0,0.537d0,0.545d0,0.553d0,0.614d0,0.906d0,1.028d0,1.389d0,1.504d0,2.533d0, &
         0.545d0,0.553d0,0.553d0,0.576d0,0.637d0,1.036d0,1.550d0,1.658d0,1.934d0,3.277d0, &
         0.545d0,0.560d0,0.568d0,0.606d0,1.174d0,1.781d0,2.563d0,3.170d0,3.791d0,4.966d0, &
         0.545d0,0.591d0,0.614d0,1.259d0,2.065d0,2.824d0,3.761d0,4.498d0,5.902d0,6.148d0/ 

    data ((vnorm(8,j,k),j=1,10),k=1,13)/  &
         0.514d0,0.539d0,0.596d0,0.694d0,0.832d0,1.004d0,1.444d0,1.869d0,2.203d0,2.538d0, &
         0.514d0,0.539d0,0.571d0,0.645d0,0.751d0,0.906d0,1.387d0,1.779d0,2.056d0,2.317d0, &
         0.514d0,0.547d0,0.555d0,0.612d0,0.702d0,0.824d0,1.281d0,1.681d0,1.934d0,2.203d0, &
         0.514d0,0.539d0,0.555d0,0.588d0,0.653d0,0.743d0,1.028d0,1.404d0,1.624d0,2.024d0, &
         0.514d0,0.539d0,0.547d0,0.555d0,0.588d0,0.710d0,0.889d0,1.191d0,1.420d0,1.820d0, &
         0.514d0,0.522d0,0.522d0,0.539d0,0.563d0,0.710d0,0.849d0,1.044d0,1.208d0,1.534d0, &
         0.514d0,0.481d0,0.506d0,0.514d0,0.539d0,0.694d0,0.824d0,1.028d0,1.200d0,1.371d0, &
         0.514d0,0.481d0,0.514d0,0.547d0,0.563d0,0.702d0,0.898d0,1.134d0,1.297d0,1.501d0, &
         0.514d0,0.490d0,0.514d0,0.555d0,0.588d0,0.726d0,0.955d0,1.265d0,1.379d0,1.648d0, &
         0.514d0,0.547d0,0.547d0,0.571d0,0.604d0,0.767d0,1.036d0,1.355d0,1.550d0,3.142d0, &
         0.514d0,0.563d0,0.579d0,0.604d0,0.612d0,0.832d0,1.909d0,2.848d0,3.917d0,4.790d0, &
         0.514d0,0.522d0,0.563d0,0.677d0,0.767d0,1.420d0,2.040d0,3.158d0,4.863d0,6.291d0, &
         0.514d0,0.588d0,0.588d0,0.612d0,0.824d0,2.032d0,3.109d0,4.969d0,6.846d0,7.695d0/ 

    data ((vnorm(9,j,k),j=1,10),k=1,13)/  &
         0.572d0,0.608d0,0.679d0,0.751d0,0.831d0,1.001d0,1.377d0,1.913d0,2.512d0,2.879d0, &
         0.572d0,0.572d0,0.608d0,0.679d0,0.760d0,0.930d0,1.243d0,1.707d0,2.369d0,2.700d0, &
         0.572d0,0.563d0,0.590d0,0.644d0,0.706d0,0.831d0,1.171d0,1.618d0,2.190d0,2.378d0, &
         0.572d0,0.554d0,0.563d0,0.599d0,0.662d0,0.760d0,1.010d0,1.502d0,2.011d0,2.235d0, &
         0.572d0,0.545d0,0.563d0,0.590d0,0.626d0,0.715d0,0.885d0,1.323d0,1.815d0,2.119d0, &
         0.572d0,0.527d0,0.554d0,0.572d0,0.608d0,0.670d0,0.724d0,1.144d0,1.618d0,1.868d0, &
         0.572d0,0.545d0,0.572d0,0.572d0,0.599d0,0.662d0,0.724d0,1.117d0,1.484d0,1.761d0, &
         0.572d0,0.554d0,0.590d0,0.599d0,0.608d0,0.679d0,0.760d0,1.216d0,1.582d0,1.922d0, &
         0.572d0,0.572d0,0.599d0,0.608d0,0.635d0,0.715d0,0.822d0,1.377d0,1.707d0,2.056d0, &
         0.572d0,0.590d0,0.608d0,0.635d0,0.662d0,0.742d0,0.912d0,1.529d0,3.075d0,4.693d0, &
         0.572d0,0.590d0,0.626d0,0.644d0,0.670d0,0.760d0,1.109d0,1.564d0,3.111d0,4.702d0, &
         0.572d0,0.599d0,0.644d0,0.662d0,0.688d0,0.822d0,1.788d0,2.816d0,5.346d0,7.295d0, &
         0.572d0,0.608d0,0.662d0,0.670d0,0.715d0,1.851d0,3.227d0,4.810d0,6.669d0,9.557d0/ 
    
    data ((vnorm(10,j,k),j=1,10),k=1,13)/   &
         0.552d0,0.606d0,0.639d0,0.671d0,0.704d0,0.899d0,1.223d0,2.479d0,3.194d0,3.573d0, &
         0.552d0,0.574d0,0.606d0,0.628d0,0.682d0,0.855d0,1.148d0,2.339d0,2.642d0,3.378d0, &
         0.552d0,0.563d0,0.552d0,0.595d0,0.639d0,0.834d0,1.061d0,2.014d0,2.404d0,2.891d0, &
         0.552d0,0.563d0,0.509d0,0.552d0,0.628d0,0.801d0,0.985d0,1.689d0,2.176d0,2.653d0, &
         0.552d0,0.574d0,0.509d0,0.520d0,0.585d0,0.747d0,0.888d0,1.332d0,1.970d0,2.458d0, &
         0.552d0,0.531d0,0.509d0,0.509d0,0.531d0,0.682d0,0.801d0,1.191d0,1.819d0,2.425d0, &
         0.552d0,0.498d0,0.498d0,0.498d0,0.520d0,0.639d0,0.747d0,1.126d0,1.711d0,2.317d0, &
         0.552d0,0.498d0,0.509d0,0.509d0,0.541d0,0.671d0,0.780d0,1.278d0,1.862d0,2.598d0, &
         0.552d0,0.498d0,0.509d0,0.520d0,0.574d0,0.693d0,0.812d0,1.602d0,2.035d0,2.793d0, &
         0.552d0,0.520d0,0.520d0,0.531d0,0.595d0,0.725d0,0.844d0,1.916d0,2.588d0,3.768d0, &
         0.552d0,0.531d0,0.541d0,0.574d0,0.628d0,0.780d0,1.039d0,2.349d0,3.313d0,5.652d0, &
         0.552d0,0.574d0,0.563d0,0.606d0,0.660d0,0.812d0,1.797d0,3.010d0,5.478d0,7.492d0, &
         0.552d0,0.650d0,0.671d0,0.704d0,0.801d0,1.029d0,2.436d0,3.465d0,7.828d0,10.578d0/

    data ((vnorm(11,j,k),j=1,10),k=1,13)/   &
         0.518d0,0.576d0,0.605d0,0.633d0,0.662d0,0.864d0,1.238d0,2.620d0,3.455d0,3.887d0, &
         0.518d0,0.547d0,0.576d0,0.576d0,0.633d0,0.835d0,1.123d0,2.447d0,2.821d0,3.656d0, &
         0.518d0,0.518d0,0.518d0,0.547d0,0.605d0,0.806d0,1.036d0,2.102d0,2.533d0,3.080d0, &
         0.518d0,0.518d0,0.461d0,0.518d0,0.576d0,0.777d0,0.950d0,1.727d0,2.274d0,2.821d0, &
         0.518d0,0.547d0,0.461d0,0.489d0,0.547d0,0.720d0,0.864d0,1.353d0,2.044d0,2.591d0, &
         0.518d0,0.489d0,0.461d0,0.461d0,0.489d0,0.662d0,0.777d0,1.180d0,1.871d0,2.562d0, &
         0.518d0,0.461d0,0.461d0,0.461d0,0.489d0,0.605d0,0.720d0,1.123d0,1.756d0,2.418d0, &
         0.518d0,0.461d0,0.461d0,0.461d0,0.518d0,0.633d0,0.749d0,1.296d0,1.929d0,2.764d0, &
         0.518d0,0.461d0,0.461d0,0.489d0,0.547d0,0.662d0,0.777d0,1.641d0,2.130d0,2.994d0, &
         0.518d0,0.489d0,0.489d0,0.489d0,0.547d0,0.691d0,0.806d0,1.986d0,2.735d0,4.117d0, &
         0.518d0,0.489d0,0.489d0,0.547d0,0.576d0,0.749d0,1.008d0,2.476d0,3.599d0,6.334d0, &
         0.518d0,0.547d0,0.518d0,0.576d0,0.633d0,0.777d0,1.842d0,3.224d0,6.132d0,8.550d0, &
         0.518d0,0.605d0,0.633d0,0.662d0,0.777d0,1.008d0,2.562d0,3.771d0,8.953d0,12.293d0/
 
    !   compute sun zenith bin
    cc  = cos( sz * MPC_RADIANS_PER_DEGREE_R8)
    i1  = 12.d0 - (cc + 0.05d0) * 10.d0
    i2  = i1 + 1 
    if (i1 >= 11) i1 = 11 
    if (i1 == 11) i2 = i1 

    !  compute sat zenith bin 
    j1  = int(satz / 10.d0) + 1 
    j2  = j1 + 1 
    if (j1 == 10) j2 = j1 

    !  compute relative azimuth bin 
    k1  = RZ / 15.d0 + 1.d0
    k2  = k1 + 1 
    if (k1 == 13) k2 = k1 

    !  interpolate
    ierr = 0 
    do l=i1,i2  
      i = l -i1 + 1
       
    !     between r's for constant s
      do n=k1,k2 

        !        between v's for constant r and s 
        m  = n - k1 + 1
        d1 = vnorm(l,j1,n)
        d2 = vnorm(l,j2,n)
        if (d1 == d2) then
          da(m) = d1
        else
          call lineq(V(j1),V(j2),d1,d2,slope,intercept,ierr) 
          da(m) = slope * satz + intercept
        end if
      end do
      if(k1 == k2) then 
        dd(i)  = da(1) 
      else 
        call lineq(R(k1),R(k2),da(1),da(2),slope,intercept,ierr) 
        dd(i) = slope * RZ + intercept
      end if
    end do

    !  between s's using result of other interpolations 
    if(i1 == i2) then
      zlamb  = drm(i1) 
      zcloud = drcld(i1)
      anisot = dd(1)
    else
      x1 = cos(s(i1) * MPC_RADIANS_PER_DEGREE_R8) 
      x2 = cos(s(i2) * MPC_RADIANS_PER_DEGREE_R8) 
      call lineq(x1,x2,dd(1),dd(2),slope,intercept,ierr) 
      anisot = slope * cc + intercept 
      g1 = drm(i1)
      g2 = drm(i2)
      call lineq(x1,x2,g1,g2,slope,intercept,ierr) 
      zlamb  = slope * cc + intercept
      g1 = drcld(i1)
      g2 = drcld(i2)
      call lineq(x1,x2,g1,g2,slope,intercept,ierr) 
      zcloud = slope * cc + intercept 
    end if
    
    if (anisot < 0.) then 
      ierr = -1
      anisot = 1.d0 
      zlamb  = drm(i1) 
      zcloud = drcld(i1)
    end if
    
  end subroutine visocn

  !--------------------------------------------------------------------------
  ! lineq
  !--------------------------------------------------------------------------
  subroutine lineq(x1, x2, y1, y2, a, b, ierr) 
    !
    ! :Purpose: calculate slope and intercept of a line.
    !
    implicit none

    ! Arguments:
    real(8), intent(in)  :: x1   ! coordinate x of point 1
    real(8), intent(in)  :: x2   ! coordinate x of point 2
    real(8), intent(in)  :: y1   ! coordinate y of point 1
    real(8), intent(in)  :: y2   ! coordinate y of point 2
    real(8), intent(out) :: a    ! slope
    real(8), intent(out) :: b    ! intercept
    integer, intent(out) :: ierr ! error code (0=ok)
     
    ierr = 0
    
    if ( (x2 - x1) == 0.d0) then 
      ierr = -1
      return
    end if

    a = ( y2 - y1) / (x2 - x1) 
    b = y1 - a * x1 
    
  end subroutine lineq


  !--------------------------------------------------------------------------
  ! drm
  !--------------------------------------------------------------------------
  real(8) function drm(iz) 
    !
    ! :Purpose: Normalization for sun zenith angle (lambertian)
    !           for ocean.
    !
    ! :Outputs: normalization factor
    !
    implicit none

    ! Arguments:
    integer, intent(in) ::  iz  ! index

    ! Locals:
    real(8),parameter :: drf(11)=(/1.d0,1.0255d0,1.1197d0,1.2026d0,1.3472d0, &
         1.4926d0,1.8180d0,2.1980d0, 2.8180d0,3.8615d0,4.3555d0/)

    drm = drf(IZ) 
  
  end function drm
      

end module multiIRbgck_mod


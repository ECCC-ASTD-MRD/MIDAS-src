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

MODULE biasCorrectionConv_mod
  ! MODULE biasCorrectionConv_mod (prefix="bcc" category='1. High-level functionality')
  !
  ! :Purpose: Performs bias correction for conventional observations
  !
  use utilities_mod
  use obsSpaceData_mod
  use MathPhysConstants_mod
  use midasMpi_mod
  use bufr_mod
  use codtyp_mod
  use timeCoord_mod

  implicit none
  save
  private
  public               :: bcc_applyAIBcor, bcc_applyGPBcor, bcc_applyUABcor
  public               :: bcc_biasActive

  integer, parameter :: nPhases=3, nLevels=5, nAircraftMax=100000
  integer, parameter :: nStationMaxGP = 10000
  integer, parameter :: phaseLevel   = 3
  integer, parameter :: phaseAscent  = 5
  integer, parameter :: phaseDescent = 6
  integer, parameter :: phaseLevelIndex   = 1
  integer, parameter :: phaseAscentIndex  = 2
  integer, parameter :: phaseDescentIndex = 3

  integer, parameter :: nSondesMax = 18, nMandLevs = 16
  integer, parameter :: nStationMaxUA = 1000
  integer, parameter :: Index500mb = 5
  
  ! Missing values in the input bias correction files for AI, GP and UA families
  real(8), parameter :: aiMissingValue = 99.d0
  real(8), parameter :: gpMissingValue = -999.00d0
  real(8), parameter :: uaMissingValue = -99.0d0
  
  integer               :: nbAircrafts, nbGpStations, nlNbSondes
  real(8), allocatable  :: AIttCorrections(:,:,:)
  real(8), allocatable  :: ztdCorrections(:)
  real(8), allocatable  :: ttCorrections(:,:,:,:),    tdCorrections(:,:,:,:)
  real(8), allocatable  :: ttCorrectionsStn(:,:,:,:), tdCorrectionsStn(:,:,:,:)
  
  character(len=9), allocatable             :: aircraftIds(:), gpsStations(:), uaStations(:)
  character(len=8), dimension(nSondesMax)   :: sondeTypes, nlSondeTypes
  integer, dimension(nSondesMax,20)         :: nlSondeCodes
  
  real(8), dimension(nMandLevs)   :: MandLevs, tolPress
  
  logical :: bcc_aiBiasActive, bcc_gpBiasActive, bcc_uaBiasActive
  logical, allocatable :: BiasCorrPresentStype(:,:,:), BiasCorrPresentStn(:,:,:)
  logical :: initialized = .false.
  
  ! Bias correction files (must be in program working directory)
  character(len=8), parameter  :: aiBcFile = "ai_bcors", gpBcFile = "gp_bcors"
  character(len=14), parameter :: uaBcFileStype = "ua_bcors_stype", uaBcFileStn = "ua_bcors_stn"

  integer, external    :: fnom, fclos, newdate
  
  ! Namelist variables
  logical :: aiBiasActive  ! Control if bias correction is applied to aircraft data
  logical :: gpBiasActive  ! Control if bias correction is applied to ground-based GPS data
  logical :: uaBiasActive  ! Control if bias correction is applied to radiosonde data
  logical :: aiRevOnly     ! Don't apply new correction but simply reverse any old corrections for AI
  logical :: gpRevOnly     ! Don't apply new correction but simply reverse any old corrections for GP
  logical :: uaRevOnly     ! Don't apply new correction but simply reverse any old corrections for UA
  logical :: uaRejUnBcor   ! Set DATA QC flag bit 11 on to exclude uncorrected UA observations from assimilation 
  integer :: uaNbiasCat    ! Number of bias profile categories in UA bcor files, e.g. 1, or 2 for "asc" and "desc" phase categoroes
  integer :: uaNlatBands   ! Number of latitude bands in ua_bcors_stype bcor file (= 5 or 1). Set to 1 if there are no latitude bands in file.
  
  ! Min number of bias profiles (sample size) required for a given station/stype/time-of-day to use 
  ! the biases in file "ua_bcors_stn" as corrections.
  integer :: uaNprofsMin 
  
  ! Structure to hold dictionary containing the BUFR sonde type codes associated with each sonde type
  type sonde_type
     character(len=8)        :: name
     integer, dimension(20)  :: codes
  end type sonde_type
  
  type(sonde_type), allocatable  :: rs_types(:)
  
  namelist /nambiasconv/ aiBiasActive,gpBiasActive,aiRevOnly,gpRevOnly,uaBiasActive,uaRevOnly,uaNprofsMin,uaRejUnBcor,uaNbiasCat,uaNlatBands
  namelist /namsondetypes/ nlNbSondes, nlSondeTypes, nlSondeCodes
  
  ! 16 mandatory pressure levels (mb) on which radiosonde bias profiles are defined
  data MandLevs /1000.d0,925.d0,850.d0,700.d0,500.d0,400.d0,300.d0,250.d0,200.d0,150.d0,100.d0,70.d0,50.d0,30.d0,20.d0,10.d0/
  ! +/- tolerance (mb) for matching a radiosonde observation pressure to one of the mandatory levels
  data tolPress /  10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0,  5.d0,  5.d0,  5.d0,  5.d0, 2.d0, 2.d0, 2.d0, 1.d0, 1.d0/
  
CONTAINS

  !-----------------------------------------------------------------------
  ! bcc_UACorrection
  !-----------------------------------------------------------------------
  
  real(8) function bcc_UACorrection(TimeOfDayX,corr_night,corr_day)
    ! :Purpose: Returns a UA bias correction given
    !      TimeOfDayX   = (float)  0.0 (night) <= TimeOfDayX <= 1.0 (day), depends on solar_elev
    !      corr_night   = (float)  night bias correction
    !      corr_day     = (float)  day bias correction
    
    implicit none
    !Arguments:
    real(8), intent(in)  ::  TimeOfDayX, corr_night, corr_day
    
    bcc_UACorrection = MPC_missingValue_R8
    
    if ( TimeOfDayX == 0.0d0 ) then
      bcc_UACorrection = corr_night
    elseif ( TimeOfDayX == 1.0d0 ) then
      bcc_UACorrection = corr_day
    else
      if ( corr_night /= MPC_missingValue_R8 .and. corr_day /= MPC_missingValue_R8 ) then
        bcc_UACorrection = TimeOfDayX*corr_day + (1.0d0-TimeOfDayX)*corr_night
      end if
    end if
    
  end function bcc_UACorrection
  
  !-----------------------------------------------------------------------
  ! bcc_GetUACorrection
  !-----------------------------------------------------------------------
  
  subroutine bcc_GetUACorrection(var,stn_index,stype_index,stype,jCat,TimeOfDayX,latband,plevel,corr,sourceCorr)
    !
    ! :Purpose: Return a TT or TD bias correction (corr) given
    !      var            = (str) 'TT' or 'TD'
    !      stn_index      = (int) station index (in uaStations)  = -1 if station was not found in bcor file
    !      stype_index    = (int) sonde type index (in rs_types) =  0 if sonde-type code is not associated with any type
    !      stype          = (str) sonde-type from rs_types
    !      jCat           = (int) bias profile category = 1 (ascent/none) or 3 (descent)
    !      TimeOfDayX     = (float)  0.0 (night) <= TimeOfDayX <= 1.0 (day), depends on solar_elev
    !      latband        = (int) 1-5
    !      plevel         = (float) pressure level (mb) of observation
    !
    !  Also returns info regarding the source of the correction returned:
    !      sourceCorr     = (str) 'none', 'stn', or 'stype'
    !
    ! Requires TT and TD correction profiles read from UA bcor file ua_bcors_stn and stored in
    !    ttCorrectionsStn(stn_index,stype_index,j,MandLev)
    !    tdCorrectionsStn(stn_index,stype_index,j,MandLev)
    ! with backup corrections by sonde-type read from UA bcor file ua_bcors_stype and stored in
    !    ttCorrections(stype_index,latband,j,MandLev)
    !    tdCorrections(stype_index,latband,j,MandLev)
    ! where MandLev = 1 (1000 mb) to 16 (10 mb)
    ! If uaNbiasCat = 1 (jCat = 1)
    !       j = 1 (night)
    !           2 (day)
    ! If uaNbiasCat = 2 (jCat = 1 [ascent] or 3 [descent])
    !       j = 1 (night-ascent)
    !           2 (day-ascent)
    !           3 (night-descent)
    !           4 (day-descent)
    !
    ! Interpolation of the correction profiles on mandatory levels is used to get the 
    ! correction at the observation level (plevel).
    ! Persistence is applied for observations outside the range of the mandatory levels.
    !
    
    implicit none
    !Arguments:
    integer, intent(in)           ::  stn_index, stype_index, jCat, latband
    real(8), intent(in)           ::  TimeOfDayX
    real(8), intent(in)           ::  plevel
    character(len=2), intent(in)  ::  var
    character(len=8), intent(in)  ::  stype
    real(8), intent(out)          ::  corr
    character(len=5), intent(out) ::  sourceCorr
    !Locals:
    real(8), dimension(16)        ::  corr_profile_stn_day, corr_profile_stn_night, corr_profile_stype_day, corr_profile_stype_night
    real(8)                       ::  corr_day, corr_night
    real(8)                       ::  w1, w2, dp, pa, pb, ca, cb, ppp, can, cbn, cad, cbd
    integer                       ::  i, ii
    logical                       ::  ProfileExistsStn, ProfileExistsStype, doInterp

    ppp = plevel        
    
    sourceCorr = "none"
    
    ! Bias correction by station and sonde-type
    if ( stn_index == -1 ) then
      ProfileExistsStn = .false.
    else
      if ( stype_index > 0 ) then
        ! Check that correction profile for this station, sonde-type and TimeOfDay is available for use
        ProfileExistsStn = BiasCorrPresentStn(stn_index,stype_index,jCat)
      else ! There are no bias profiles for this station/stype combination
        ProfileExistsStn = .false.
      end if
    end if

    ! BACKUP: Bias correction by sonde-type
    ProfileExistsStype = .false.
    if ( trim(stype) /= 'unknown' .and. trim(stype) /= 'Others' .and. trim(stype) /= 'None' ) then
      ! Check if correction profile for this sonde-type,latband,and TimeOfDay is available for use
      ProfileExistsStype = BiasCorrPresentStype(stype_index,latband,jCat)
    end if
       
    if ( .not.ProfileExistsStn  .and. .not.ProfileExistsStype ) then
      corr = MPC_missingValue_R8
      return
    end if
    
    ! Fill the night and day bias correction profiles
    if ( var == 'TT' ) then
      if ( ProfileExistsStn ) then
        corr_profile_stn_night(:)   = ttCorrectionsStn(stn_index,stype_index,jCat,:)
        corr_profile_stn_day(:)     = ttCorrectionsStn(stn_index,stype_index,jCat+1,:)
      end if
      if ( ProfileExistsStype ) then
        corr_profile_stype_night(:) = ttCorrections(stype_index,latband,jCat,:)
        corr_profile_stype_day(:)   = ttCorrections(stype_index,latband,jCat+1,:)
      end if
    elseif ( var == 'TD' ) then
      if ( ProfileExistsStn ) then
         corr_profile_stn_night(:)   = tdCorrectionsStn(stn_index,stype_index,jCat,:)
         corr_profile_stn_day(:)     = tdCorrectionsStn(stn_index,stype_index,jCat+1,:)
      end if
      if ( ProfileExistsStype ) then
         corr_profile_stype_night(:) = tdCorrections(stype_index,latband,jCat,:)
         corr_profile_stype_day(:)   = tdCorrections(stype_index,latband,jCat+1,:)
      end if
    else 
      call utl_abort('bcc_GetUACorrection ERROR: Unsupported var '//var)
    end if
    
    corr = MPC_missingValue_R8
    doInterp = .true.

    !--------------------------------------------------------------------------------------
    ! Get the correction at the observation level (ppp)
    !--------------------------------------------------------------------------------------
    
    ! Check if ppp is outside range of levels (no interpolation possible)
    if ( ppp >= MandLevs(1) ) then
      ii = 1
    elseif ( ppp <= MandLevs(nMandLevs) ) then 
      ii = nMandLevs
    else
      ii = 0
    end if
    
    ! Check if ppp is close to one of the 16 mandatory levels (no interpolation needed)
    if ( ii == 0 ) then
      do i = 1, nMandLevs
        if ( abs(ppp-MandLevs(i)) <= tolPress(i) ) then
          ii = i
          exit
        end if
      end do
    end if
    
    ! If ppp is close to one of the 16 mandatory levels get the correction for the mandatory level: 
    ! Use correction by station with correction by stype as backup
    if ( ii > 0 ) then
      if ( ProfileExistsStn ) then
        corr_night = corr_profile_stn_night(ii)
        corr_day   = corr_profile_stn_day(ii)
        sourceCorr = "stn"
        corr = bcc_UACorrection(TimeOfDayX,corr_night,corr_day)
        if ( corr == MPC_missingValue_R8 ) then
          sourceCorr = "none"
          if ( ProfileExistsStype ) then
            corr_night = corr_profile_stype_night(ii)
            corr_day   = corr_profile_stype_day(ii)
            sourceCorr = "stype"
            corr = bcc_UACorrection(TimeOfDayX,corr_night,corr_day)
            if ( corr == MPC_missingValue_R8 ) sourceCorr = "none"
          end if
        end if
      else
        if ( ProfileExistsStype ) then
          corr_night = corr_profile_stype_night(ii)
          corr_day   = corr_profile_stype_day(ii)
          sourceCorr = "stype"
          corr = bcc_UACorrection(TimeOfDayX,corr_night,corr_day)
          if ( corr == MPC_missingValue_R8 ) sourceCorr = "none"
        end if
      end if
      doInterp = .false.
    end if
    
    ! or interpolate to get correction for observation level ppp:
    ! Use corrections by station with corrections by stype as backup
    if ( doInterp ) then
      ca = MPC_missingValue_R8
      cb = MPC_missingValue_R8
      
      if ( ProfileExistsStn ) then
        do i = 1,nMandLevs-1
          if ( ppp <= MandLevs(i) .and. ppp > MandLevs(i+1) ) then
            pb = MandLevs(i)
            pa = MandLevs(i+1)
            cbn = corr_profile_stn_night(i)
            can = corr_profile_stn_night(i+1)
            cbd = corr_profile_stn_day(i)
            cad = corr_profile_stn_day(i+1)
            exit
          end if
        end do
        sourceCorr = "stn"
        if ( TimeOfDayX == 0.0d0 ) then
          cb = cbn
          ca = can 
        elseif ( TimeOfDayX == 1.0d0 ) then
          cb = cbd
          ca = cad
        else 
          if ( cbn /= MPC_missingValue_R8 .and. cbd /= MPC_missingValue_R8 ) then
            cb = TimeOfDayX*cbd + (1.0d0-TimeOfDayX)*cbn
          end if
          if ( can /= MPC_missingValue_R8 .and. cad /= MPC_missingValue_R8 ) then
            ca = TimeOfDayX*cad + (1.0d0-TimeOfDayX)*can
          end if
        end if
      end if
      
      if ( ca == MPC_missingValue_R8 .or. cb == MPC_missingValue_R8 ) then
        if ( ProfileExistsStype ) then
          do i = 1,nMandLevs-1
            if ( ppp <= MandLevs(i) .and. ppp > MandLevs(i+1) ) then
              pb = MandLevs(i)
              pa = MandLevs(i+1)
              cbn = corr_profile_stype_night(i)
              can = corr_profile_stype_night(i+1)
              cbd = corr_profile_stype_day(i)
              cad = corr_profile_stype_day(i+1)
              exit
            end if
          end do
          sourceCorr = "stype"
          if ( TimeOfDayX == 0.0d0 ) then
            cb = cbn
            ca = can 
          elseif ( TimeOfDayX == 1.0d0 ) then
            cb = cbd
            ca = cad
          else 
            if ( cbn /= MPC_missingValue_R8 .and. cbd /= MPC_missingValue_R8 ) then
              cb = TimeOfDayX*cbd + (1.0d0-TimeOfDayX)*cbn
            end if
            if ( can /= MPC_missingValue_R8 .and. cad /= MPC_missingValue_R8 ) then
              ca = TimeOfDayX*cad + (1.0d0-TimeOfDayX)*can
            end if
          end if
        end if
      end if
      dp = log10(pb)-log10(pa)
      w1 = (log10(pb) - log10(ppp)) / dp
      w2 = (log10(ppp) - log10(pa)) / dp
      if ( ca /= MPC_missingValue_R8 .and. cb /= MPC_missingValue_R8 ) then
        corr = w1*ca + w2*cb
      elseif ( ca /= MPC_missingValue_R8 .and. max(w1,w2) == w1 ) then
        corr = ca
      elseif ( cb /= MPC_missingValue_R8 .and. max(w1,w2) == w2 ) then
        corr = cb
      else
        corr = MPC_missingValue_R8
        sourceCorr = "none"
      end if
    end if
  
  end subroutine bcc_GetUACorrection
  
  !-----------------------------------------------------------------------
  ! bcc_StationIndex
  !-----------------------------------------------------------------------  
  
  integer function bcc_StationIndex(Station)
    !
    ! :Purpose: Return the Station index (order in array uaStations) corresponding to Station
    !
    ! Returns -1 if Station is not found in uaStations.
    !
    implicit none
    !Arguments:
    character(len=9), intent(in)  :: Station
    !Locals:
    integer    :: i
    
    if ( allocated(uaStations) ) then
      bcc_StationIndex = -1
      do i = 1,nStationMaxUA
        if ( trim(uaStations(i)) == trim(Station) ) then
          bcc_StationIndex = i
          exit
        end if
      end do
    else 
      call utl_abort('bcc_StationIndex: ERROR: array uaStations not allocated!')
    end if
    
  end function bcc_StationIndex
  
  !-----------------------------------------------------------------------
  ! bcc_SondeIndex
  !-----------------------------------------------------------------------
  
  integer function bcc_SondeIndex(SondeType)
    !
    ! :Purpose: Return the Sonde Type index (order in array rs_types) corresponding to the SondeType
    !
    ! Requires array of sonde_type structures (rs_types) to be allocated and filled.
    ! Returns -1 if SondeType is not found in rs_types.
    !
    implicit none
    !Arguments:
    character(len=8), intent(in)  :: SondeType
    !Locals:
    integer    :: i, ntypes
    
    if ( allocated(rs_types) ) then
      ntypes = nlNbSondes
    else 
      call utl_abort('bcc_SondeIndex: ERROR: array rs_types not allocated!')
    end if
    
    bcc_SondeIndex = -1
    do i = 1,ntypes
      if ( trim(rs_types(i)%name) == trim(SondeType) ) then
        bcc_SondeIndex = i
        exit
      end if
    end do

  end function bcc_SondeIndex
  
  !-----------------------------------------------------------------------
  ! bcc_GetSondeType
  !-----------------------------------------------------------------------
  subroutine bcc_GetSondeType(code,stype,stype_index)
    !
    ! :Purpose: Returns the sonde type and index given a BUFR table sonde type code (BUFR ele 002011)
    !           Returns stype='unknown', stype_index=0 if code is not found.
    !
    ! Requires array of sonde_type structures (rs_types) to be allocated and filled with the 
    ! sonde type codes associated with each sonde-type (read from namelist).
    !
    implicit none
    !Arguments:
    integer, intent(in)            :: code
    character(len=8), intent(out)  :: stype
    integer, intent(out)           :: stype_index
    !Locals:
    integer  :: i, ntypes, icode
    
    if ( allocated(rs_types) ) then
      ntypes = nlNbSondes
    else 
      call utl_abort('bcc_GetSondeType: ERROR: array rs_types not allocated!')
    end if
    
    if ( code == 190 .or. code == 192 ) then     ! NCAR dropsonde (BUFR only)
      icode = 13                                 ! RS92 code
    elseif ( code == 191 .or. code == 193 ) then ! NCAR dropsonde (BUFR only)
      icode = 41                                 ! RS41 code
    elseif ( code >=100 .and. code < 200 ) then
      icode = code - 100
    else 
      icode = code
    end if
    
    stype = 'unknown'
    stype_index = 0
    do i = 1,ntypes
      if ( ANY(rs_types(i)%codes == icode) ) then
        stype = rs_types(i)%name 
        stype_index = i
        exit
      end if
    end do
    
  end subroutine bcc_GetSondeType
  
  !-----------------------------------------------------------------------
  ! bcc_GetSolarElevation
  !-----------------------------------------------------------------------
  subroutine bcc_GetSolarElevation(lat,lon,date,time,solar_elev)
    !
    ! :Purpose: Returns the solar elevation angle (degrees) given lat,lon,date(yyyymmdd),time(hhmm)
    !
    !    lat,  lon  = obsSpaceData header column OBS_LAT, OBS_LON (radians)
    !    date, time = obsSpaceData header column OBS_DAT (yyyymmdd), OBS_ETM (hhmm)  -or-
    !       datestamp = tim_getDatestamp()   -- date stamp for central (analysis) time
    !       ier  = newdate(datestamp,date,time,-3)
    !       time = time/10000
    !
    implicit none
    !Arguments:
    integer, intent(in)  :: date, time    ! yyyymmdd, hhmm
    real(8), intent(in)  :: lat, lon      ! radians
    real(8), intent(out) :: solar_elev    ! degrees
    !Locals:
    integer, dimension(13)   :: days = (/0,31,28,31,30,31,30,31,31,30,31,30,31/)
    integer, dimension(7)    :: leap_years = (/2016,2020,2024,2028,2032,2036,2040/)
    integer   :: yy, mmdd, mm, dd, hh, nn, doy
    real(8)   :: timeUTC, timeLCL, sol_dec, hour_angle, csz, sza

    yy   = date/10000
    mmdd = date - yy*10000
    mm   = mmdd/100
    dd   = mmdd - mm*100 
    hh   = time/100
    nn   = time - hh*100
    
    if ( ANY(leap_years == yy) ) days(3) = 29
    doy = SUM(days(1:mm)) + dd
    timeUTC = float(hh) + float(nn)/60.0d0
    timeLCL = timeUTC + (lon*MPC_DEGREES_PER_RADIAN_R8)/15.0d0
    if ( timeLCL < 0.0d0  ) timeLCL = 24.0d0 + timeLCL
    if ( timeLCL > 24.0d0 ) timeLCL = timeLCL - 24.0d0
    sol_dec = 0.4093d0 * sin(((2.0d0*MPC_PI_R8)/365.0d0)*(284.0d0 + float(doy)))
    hour_angle = 15.0d0*(timeLCL-12.0d0) * MPC_RADIANS_PER_DEGREE_R8
    csz = sin(sol_dec)*sin(lat) + cos(sol_dec)*cos(lat)*cos(hour_angle)
    sza = MPC_DEGREES_PER_RADIAN_R8 * acos(csz)
    solar_elev = 90.0d0 - sza
  
  end subroutine bcc_GetSolarElevation
  
  !-----------------------------------------------------------------------
  ! bcc_GetTimeOfDay
  !-----------------------------------------------------------------------
  subroutine bcc_GetTimeOfDay(solar_elev,TimeOfDayX)
    !
    ! :Purpose: Returns the time-of-day x value (0.0(night) <= x <= 1.0(day))
    !
    implicit none
    !Arguments:
    real(8), intent(in)                ::  solar_elev     ! degrees
    real(8), intent(out)               ::  TimeOfDayX
    
    if (solar_elev < -7.5d0) then 
      TimeOfDayX = 0.0d0
    elseif (solar_elev < 22.5d0) then
      TimeOfDayX = (solar_elev+7.5d0)/(22.5d0+7.5d0)
    else
      TimeOfDayX = 1.0d0
    end if
    
  end subroutine bcc_GetTimeOfDay
  
  !----------------------------------------------------------------------------------------
  ! bcc_uaPhase
  !----------------------------------------------------------------------------------------
  integer function bcc_uaPhase(CodeType)
    !
    ! :Purpose: Returns the radiosonde phase (1=ascent, 2=descent) given a header code type
    !
    implicit none
    !Arguments:
    integer, intent(in)   ::  CodeType
    
    if (CodeType == 37) then
      bcc_uaPhase = 2
    else
      bcc_uaPhase = 1
    end if
    
  end function bcc_uaPhase  
  
  !-----------------------------------------------------------------------
  ! bcc_LatBand
  !-----------------------------------------------------------------------
  integer function bcc_LatBand(lat)
    !
    ! :Purpose: Returns latitude band number given latitude (radians)
    !
    implicit none
    !Arguments:
    real(8), intent(in) ::  lat    ! radians
    !Locals:
    real(8)             ::  zlat
    
    if ( uaNlatBands /= 5 ) then
      bcc_LatBand = -1
    else
      zlat = MPC_DEGREES_PER_RADIAN_R8 * lat
      if (zlat < -90.0d0 ) bcc_LatBand = 1  ! Should never be the case!
      if   (zlat >= -90.0d0   .and. zlat < -60.0d0) then
         bcc_LatBand = 1
      elseif (zlat >= -60.0d0 .and. zlat < -20.0d0) then
         bcc_LatBand = 2
      elseif (zlat >= -20.0d0 .and. zlat <  20.0d0) then
         bcc_LatBand = 3
      elseif (zlat >= 20.0d0  .and. zlat <  60.0d0) then
         bcc_LatBand = 4
      else
         bcc_LatBand = 5
      end if
    end if
    
  end function bcc_LatBand
    
  
  !-----------------------------------------------------------------------
  ! bcc_readConfig
  !-----------------------------------------------------------------------
  subroutine bcc_readConfig()
    !
    ! :Purpose: Read NAMBIASCONV namelist section and NAMSONDETYPES section if uaBiasActive=true
    !
    ! FOR AN EXAMPLE OF NAMELIST WITH NAMSONDETYPES SECTION SEE
    !   /home/mac003/temp/nml.default_bgck
    !
    implicit none
    !Locals:
    integer  :: ierr, nulnam, i
  
    ! set default values for namelist variables
    aiBiasActive = .false.  ! bias correct AI data (TT)
    gpBiasActive = .false.  ! bias correct GP data (ZTD)
    uaBiasActive = .false.  ! bias correct UA data (TT,ES)
    aiRevOnly    = .false.  ! AI: don't apply new correction but simply reverse any old corrections
    gpRevOnly    = .false.  ! GP: don't apply new correction but simply reverse any old corrections
    uaRevOnly    = .false.  ! UA: don't apply new correction but simply reverse any old corrections
    uaRejUnBcor  = .false.  ! UA: Set DATA QC flag bit 11 on to exclude uncorrected UA observations from assimilation
    uaNprofsMin  = 100      ! UA: If Nprofs for a given stn/stype < uaNprofsMin, flag the bcors as unusable
    uaNbiasCat   = 1        ! UA: Number of bias profile categories in UA bcor files: 1, or 2 for "asc" and "desc" categories
    uaNlatBands  = 1        ! UA: Number of latitude bands in ua_bcors_stype file (1 or 5): 1 = no bands (global biases).
    ! read in the namelist NAMBIASCONV
    if ( utl_isNamelistPresent('nambiasconv','./flnml') ) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=nambiasconv,iostat=ierr)
      if ( ierr /= 0 )  call utl_abort('bcc_readConfig: Error reading namelist section NAMBIASCONV')
      if ( mpi_myid == 0 ) write(*,nml=nambiasconv)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'bcc_readconfig: NAMBIASCONV section is missing in the namelist. The default value will be taken.'
    end if
        
    bcc_aiBiasActive = aiBiasActive
    bcc_gpBiasActive = gpBiasActive
    bcc_uaBiasActive = uaBiasActive
    
    initialized = .true.
    
    ! read in the namelist NAMSONDETYPES
    if ( uaBiasActive .and. .not.uaRevOnly ) then
      if ( utl_isNamelistPresent('namsondetypes','./flnml') ) then
        nlNbSondes = 0
        nlSondeTypes(:)   = 'empty'
        nlSondeCodes(:,:) = -9
        nulnam = 0
        ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
        read(nulnam,nml=namsondetypes,iostat=ierr)
        if ( ierr /= 0 )  call utl_abort('bcc_readConfig: Error reading namelist section NAMSONDETYPES')
        if ( mpi_myid == 0 ) write(*,nml=namsondetypes)
        ierr = fclos(nulnam)
      else 
        write(*,*)
        call utl_abort('bcc_readconfig ERROR: NAMSONDETYPES section is missing in the namelist!')
      end if
      if ( nlNbSondes > nSondesMax ) call utl_abort('bcc_readconfig ERROR: Number of sonde types in namelist exceeds nSondesMax!')
      if ( nlNbSondes <= 0 )         call utl_abort('bcc_readconfig ERROR: Number of sonde types in namelist <= 0!')
      allocate( rs_types(nlNbSondes) )
      do i = 1,nlNbSondes
         rs_types(i)%name     = nlSondeTypes(i)
         rs_types(i)%codes(:) = nlSondeCodes(i,:)
      end do
    end if
    
  end subroutine bcc_readConfig

  !-----------------------------------------------------------------------
  ! bcc_readAIBiases
  !-----------------------------------------------------------------------
  subroutine bcc_readAIBiases(biasEstimateFile)
    !
    ! :Purpose: Read aircraft (AI) TT bias estimates from bias file and fill bias correction array ttCorrections. 
    !           The first line of the file is the number of aircraft plus one.
    !           The rest of the file gives 15 values of Mean O-A for each aircraft, with each (AC,value) line written in format "a9,1x,f6.2".
    !           The order is the same as what is written by genbiascorr script genbc.aircraft_bcor.py.
    !           The first "aircraft" (AC name = BULKBCORS) values are the bulk biases by layer for All-AC (first 5 values), AIREP/ADS 
    !           (second 5 values) and AMDAR/BUFR (last 5 values).
    !           Missing value = 99.0.
    !
    implicit none
    !Arguments:
    character(len=*), intent(in) :: biasEstimateFile

    !Locals:
    integer :: ierr, nulcoeff
    integer :: stationIndex, phaseIndex, levelIndex
    real(8) :: biasEstimate, correctionValue
    character(len=9) :: stationId

    if (.not.initialized) call bcc_readConfig()
    
    if ( aiRevOnly ) return
    
    nulcoeff = 0
    ierr = fnom(nulcoeff, biasEstimateFile, 'FMT+R/O', 0)
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readAIBiases: unable to open airplanes bias correction file ' // biasEstimateFile )
    end if
    read (nulcoeff, '(i5)', iostat=ierr ) nbAircrafts
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readAIBiases: error 1 while reading airplanes bias correction file ' // biasEstimateFile )
    end if
    
    allocate( AIttCorrections(nAircraftMax,nPhases,nLevels) )
    allocate( aircraftIds(nAircraftMax) )

    AIttCorrections(:,:,:) =  MPC_missingValue_R8
   
    do stationIndex=1,nbAircrafts
      do phaseIndex=1,3
        do levelIndex=1,5
          read (nulcoeff, *, iostat=ierr) stationId, biasEstimate
          if ( ierr /= 0 ) then
            call utl_abort('bcc_readAIBiases: error 2 while reading airplanes bias correction file ' // biasEstimateFile )
          end if
          if ( biasEstimate == aiMissingValue ) then
            correctionValue = MPC_missingValue_R8
          else
            correctionValue = -1.0D0*biasEstimate
          end if
          AIttCorrections(stationIndex,phaseIndex,levelIndex) = correctionValue
          aircraftIds(stationIndex)                         = stationId
          !print*, stationIndex, phaseIndex, levelIndex, aircraftIds(stationIndex), AIttCorrections(stationIndex,phaseIndex,levelIndex)
        end do
      end do
    end do
    ierr = fclos(nulcoeff)

    ! Check for bulk bias corrections at start of file
    if ( aircraftIds(1) /= "BULKBCORS" ) then
      call utl_abort('bcc_readAIBiases: ERROR: Bulk bias corrections are missing in bias correction file!' )
    end if

  end subroutine bcc_readAIBiases

  !-----------------------------------------------------------------------
  ! bcc_applyAIBcor
  !-----------------------------------------------------------------------
  subroutine bcc_applyAIBcor(obsSpaceData)
    !
    ! :Purpose:  to apply aircraft (AI) bias corrections to observations in ObsSpaceData
    !
    implicit none
    !Arguments:
    type(struct_obs)        :: obsSpaceData
    !Locals:
    integer  :: headerIndex, bodyIndex, codtyp
    integer  :: flag, phase, bufrCode
    integer  :: phaseIndex, levelIndex, stationIndex, stationNumber
    integer  :: countTailCorrections,  countBulkCorrections
    integer  :: headerFlag
    real(8)  :: corr, tt, oldCorr, pressure
    character(len=9) :: stnid, stnId1, stnId2

    if (.not.initialized) call bcc_readConfig()
    
    if ( .not. aiBiasActive ) return

    write(*,*) "bcc_applyAIBcor: start"

    if ( .not.aiRevOnly ) call bcc_readAIBiases(aiBcFile)
    
    countTailCorrections = 0
    countBulkCorrections = 0

    call obs_set_current_header_list(obsSpaceData,'AI')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER
      
      headerFlag = obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex )
      
      call obs_set_current_body_list(obsSpaceData, headerIndex)

      BODY: do

        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY 

        bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex )

        if ( bufrCode == BUFR_NETT ) then
          tt      = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
          flag    = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
          oldCorr = obs_bodyElem_r(obsSpaceData, OBS_BCOR, bodyIndex )
          corr = MPC_missingValue_R8
          
          if ( tt /= MPC_missingValue_R8 ) then
          
            if ( btest(flag, 6) .and. oldCorr /= MPC_missingValue_R8 ) then
              if ( btest(headerFlag, 15) ) then
                tt = tt - oldCorr
              else
                tt = tt + oldCorr
              end if
              flag = ibclr(flag, 6)
            end if
            if ( aiRevOnly ) corr = 0.0D0
             
            if ( .not.aiRevOnly ) then
                
              pressure = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex ) * MPC_MBAR_PER_PA_R8

              ! Get level index and current (mar 2020) static bulk corrections applied to AI TT data at derivate stage
              if ( (pressure <= 1100.d0) .and. (pressure > 700.d0) ) then
                corr = 0.0D0
                levelIndex = 5
              else if ( (pressure <= 700.d0)  .and. (pressure > 500.d0) ) then
                corr = -0.1D0
                levelIndex = 4
              else if ( (pressure <= 500.d0)  .and. (pressure > 400.d0) ) then
                corr = -0.2D0
                levelIndex = 3
              else if ( (pressure <= 400.d0)  .and. (pressure > 300.d0) ) then
                corr = -0.3D0
                levelIndex = 2
              else if ( (pressure <= 300.d0)  .and. (pressure > 100.d0) ) then
                corr = -0.5D0
                levelIndex = 1
              else 
                levelIndex = 0
                corr = 0.0D0
              end if
 
              codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)

              ! Default bulk corrections read from bcor file (applied if dynamic corrections are not availble for the aircraft)
              select case(  trim( codtyp_get_name(codtyp) ) )
                case('airep','ads')
                  phaseIndex = phaseAscentIndex
                case('amdar','acars')
                  phaseIndex = phaseDescentIndex
                case default
                  write(*,*) 'bcc_applyAIBcor: codtyp=', codtyp
                  call utl_abort('bcc_applyAIBcor: unknown codtyp') 
              end select

              if ( levelIndex /= 0 ) then
                if ( AIttCorrections(1,phaseIndex,levelIndex) /= MPC_missingValue_R8 ) corr = AIttCorrections(1,phaseIndex,levelIndex)
                countBulkCorrections = countBulkCorrections + 1
              end if

              headerIndex = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex  )
              stnid = trim( obs_elem_c(obsSpaceData,'STID',headerIndex) )

              ! on verifie si la station est dans le dictionnaire du fichier de correction de biais
              !---------------------------------------------------------------------------------
              stationNumber = 0
              stnId2 = trim(stnid)
              do stationIndex = 1, nbAircrafts
                stnId1 = trim(aircraftIds(stationIndex))
                if ( stnId2(2:9) == stnId1(1:8) ) stationNumber = stationIndex
              end do
            
              phase =  obs_headElem_i(obsSpaceData, OBS_PHAS, headerIndex  )

              ! If the aircraft is in the bias correction file, get the new correction
              ! AIttCorrections(stationNumber,phaseIndex,levelIndex) where
              !     stationNumber = index for this AC in bias correction file (0 if not found)
              !     phaseIndex   = index for the 3 phases of flight (level, asc, desc)
              !     levelIndex   = index for the 5 layers (100-300, 300-400,400-500,500-700,700-1100)
              !  and use it instead of bulk value (if it is not missing value).
              if ( stationNumber /= 0 ) then 
                phaseIndex = 0
                if ( phase == phaseLevel   ) phaseIndex = phaseLevelIndex
                if ( phase == phaseAscent  ) phaseIndex = phaseAscentIndex
                if ( phase == phaseDescent ) phaseIndex = phaseDescentIndex
                if ( levelIndex /= 0 .and. phaseIndex /= 0 ) then
                  if ( AIttCorrections(stationNumber,phaseIndex,levelIndex) /= MPC_missingValue_R8 ) then
                    corr = AIttCorrections(stationNumber,phaseIndex,levelIndex)
                    countTailCorrections = countTailCorrections + 1
                    countBulkCorrections = countBulkCorrections - 1
                  end if
                end if
              end if
            
              ! Apply the bias correction (bulk or new) and set the "bias corrected" bit in TT data flag ON
              tt = tt + corr
              flag = ibset(flag, 6)
            
            end if
            
          end if

          call obs_bodySet_r( obsSpaceData, OBS_BCOR, bodyIndex, corr )
          call obs_bodySet_r( obsSpaceData, OBS_VAR , bodyIndex, tt   )
          call obs_bodySet_i( obsSpaceData, OBS_FLG , bodyIndex, flag )        

        end if
        
      end do BODY
      
      headerFlag = ibset(headerFlag, 15)
      call obs_headSet_i( obsSpaceData, OBS_ST1, headerIndex, headerFlag )
      
    end do HEADER
    
    if ( countBulkCorrections + countTailCorrections /= 0 ) then
      write (*, '(a50, i10)' ) "bcc_applyAIBcor: Number of obs with TT bulk correction  = ", countBulkCorrections
      write (*, '(a50, i10)' ) "bcc_applyAIBcor: Number of obs with TT tail correction  = ", countTailCorrections
    else
      write(*,*) "bcc_applyAIBcor: No AI data found"
    end if
    
    if ( allocated(AIttCorrections) ) deallocate(AIttCorrections)
    if ( allocated(aircraftIds)   ) deallocate(aircraftIds)
    
    write(*,*) "bcc_applyAIBcor: end"
    
  end subroutine bcc_applyAIBcor

  !-----------------------------------------------------------------------
  ! bcc_readGPBiases
  !-----------------------------------------------------------------------
  subroutine bcc_readGPBiases(biasEstimateFile)
    !
    ! :Purpose: Read GB-GPS bias estimates (mean ZTD O-A [mm] by station) and fill bias correction array ztdCorrections.
    !           Missing value = -999.00
    !
    implicit none
    !Arguments:
    character(len=*), intent(in) :: biasEstimateFile
    !Locals:
    integer :: ierr, nulcoeff
    integer :: stationIndex
    real(8) :: biasEstimate
    character(len=9) :: stationId

    if (.not.initialized) call bcc_readConfig()
    
    if ( gpRevOnly ) return
    
    nulcoeff = 0
    ierr = fnom(nulcoeff, biasEstimateFile, 'FMT+R/O', 0)
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readGPBiases: unable to open GB-GPS bias correction file ' // biasEstimateFile )
    end if
    read (nulcoeff, '(i5)', iostat=ierr ) nbGpStations
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readGPBiases: error 1 while reading GB-GPS bias correction file ' // biasEstimateFile )
    end if

    allocate( ztdCorrections(nStationMaxGP) )
    allocate( gpsStations(nStationMaxGP)  )
    
    ztdCorrections(:) =  MPC_missingValue_R8
    
    do stationIndex=1,nbGpStations
       read (nulcoeff, *, iostat=ierr) stationId, biasEstimate
       if ( ierr /= 0 ) then
          call utl_abort('bcc_readGPBiases: error 2 while reading GB-GPS bias correction file ' // biasEstimateFile )
       end if
       if ( biasEstimate /= gpMissingValue ) ztdCorrections(stationIndex) = -1.0D0*(biasEstimate/1000.0D0)  ! mm to m
       gpsStations(stationIndex) = stationId
    end do
    ierr = fclos(nulcoeff)
    
  end subroutine bcc_readGPBiases

  !-----------------------------------------------------------------------
  ! bcc_applyGPBcor
  !-----------------------------------------------------------------------
  subroutine bcc_applyGPBcor(obsSpaceData)
    !
    ! :Purpose:  to apply GB-GPS (GP) ZTD bias corrections to ZTD observations in ObsSpaceData
    !
    implicit none
    !Arguments:
    type(struct_obs)  :: obsSpaceData
    !Locals:
    integer  :: headerIndex, bodyIndex
    integer  :: flag, bufrCode
    integer  :: stationIndex, stationNumber
    integer  :: nbCorrected
    real(8)  :: corr, ztd, oldCorr
    character(len=9) :: stnid, stnId1, stnId2

    if (.not.initialized) call bcc_readConfig()
    
    if ( .not. gpBiasActive ) return

    write(*,*) "bcc_applyGPBcor: start"

    if ( .not.gpRevOnly ) call bcc_readGPBiases(gpBcFile)
    
    nbCorrected = 0

    call obs_set_current_header_list(obsSpaceData,'GP')
    
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER
      
      call obs_set_current_body_list(obsSpaceData, headerIndex)

      BODY: do

        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY 

        bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex )

        if ( bufrCode == BUFR_NEZD ) then
          
          ztd     = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
          oldCorr = obs_bodyElem_r(obsSpaceData, OBS_BCOR, bodyIndex )
          flag    = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
          
          corr = MPC_missingValue_R8
          
          if ( ztd /= MPC_missingValue_R8 ) then  
            
            ! Remove any previous bias correction
            if ( btest(flag, 6) .and. oldCorr /= MPC_missingValue_R8 ) then
              ztd = ztd - oldCorr
              flag = ibclr(flag, 6)
            end if

            if ( .not. gpRevOnly ) then
              headerIndex = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex  )
              stnid = trim( obs_elem_c(obsSpaceData,'STID',headerIndex) )

              ! on verifie si la station est dans le dictionnaire du fichier de correction de biais
              ! ---------------------------------------------------------------------------------
              stationNumber = 0
              stnId2 = trim(stnid)
              do stationIndex = 1, nbGpStations
                stnId1 = trim(gpsStations(stationIndex))
                if ( stnId2 == stnId1 ) stationNumber = stationIndex
              end do
               
              if (stationNumber /= 0) then 
                corr = ztdCorrections(stationNumber)
              end if
               
              ! Apply the bias correction and set the "bias corrected" bit in ZTD data flag ON
              if ( corr /= MPC_missingValue_R8 ) then
                ztd = ztd + corr
                nbCorrected = nbCorrected + 1
                flag = ibset(flag, 6)
              else 
                corr = 0.0D0
              end if
            else
              corr = 0.0D0
            end if

          end if

          call obs_bodySet_r( obsSpaceData, OBS_BCOR, bodyIndex, corr )
          call obs_bodySet_r( obsSpaceData, OBS_VAR , bodyIndex, ztd  )
          call obs_bodySet_i( obsSpaceData, OBS_FLG , bodyIndex, flag )
           
        end if
        
      end do BODY
    end do HEADER
    
    if ( nbCorrected /= 0 ) then
      write (*, '(a50, i10)' ) "bcc_applyGPBcor: Number of ZTD observations corrected  = ", nbCorrected
    else 
      write(*,*) "bcc_applyGPBcor: No GP data bias corrections made"
    end if
    
    if ( allocated(ztdCorrections) ) deallocate( ztdCorrections )
    if ( allocated(gpsStations) )    deallocate( gpsStations  )
    
    write(*,*) "bcc_applyGPBcor: end"
    
  end subroutine bcc_applyGPBcor


  !-----------------------------------------------------------------------
  ! bcc_readUABcorStype
  !-----------------------------------------------------------------------
  subroutine bcc_readUABcorStype(biasCorrectionFileName,nGroups)
    !
    ! :Purpose: Read night and day TT, TD biases by SONDE TYPE and latitude band on 16 mandatory levels for UA family. 
    !           The first line is the sonde-type (and latitude band if uaNlatBands=5). Sonde type = 'END' for end of file.
    !           nGroups groups of of 16 lines follow, e.g., nGroups = 2: one group for ascending sonde observations and 
    !                                                                    one for descending sonde observations
    !               -- 16 lines are 16 mandatory levels from 100000 Pa to 1000 Pa
    !               -- Lines contain four values: TT night-bias, TT day-bias, TD night-bias, TD day-bias
    !           Missing value = -99.0.
    !
    implicit none
    !Arguments:
    character(len=*), intent(in) :: biasCorrectionFileName
    integer, intent(in)          :: nGroups

    !Locals:
    integer :: ierr, nulcoeff
    integer :: istype, j, ilatband, levelIndex, jj, jmax
    real(8) :: TTbiasNight, TTbiasDay, TDbiasNight, TDbiasDay
    character(len=8) :: stype

    if ( uaRevOnly ) return
    
    nulcoeff = 0
    ierr = fnom(nulcoeff, biasCorrectionFileName, 'FMT+R/O', 0)
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readUABcorStype: unable to open radiosonde bias correction file ' // biasCorrectionFileName )
    end if

    jmax = nGroups*2
    
    allocate( ttCorrections(nSondesMax,uaNlatBands,jmax,nMandLevs) )
    allocate( tdCorrections(nSondesMax,uaNlatBands,jmax,nMandLevs) )
    allocate( BiasCorrPresentStype(nSondesMax,uaNlatBands,jmax) )

    ttCorrections(:,:,:,:) =  MPC_missingValue_R8
    tdCorrections(:,:,:,:) =  MPC_missingValue_R8
    sondeTypes(:) = 'empty'
    
    BiasCorrPresentStype(:,:,:) = .false.
   
    main_loop: do 
      
      if ( uaNlatBands == 1 ) then
        read (nulcoeff, *, iostat=ierr) stype
        ilatband = 1
      else
        read (nulcoeff, *, iostat=ierr) stype, ilatband
      end if
      if ( ierr /= 0 ) then
        call utl_abort('bcc_readUABcorStype: ERROR reading stype, ilatband in radiosonde bias correction file ' // biasCorrectionFileName )
      end if
      if ( trim(stype) == 'END' ) exit main_loop
      istype = bcc_SondeIndex(stype)
      if ( istype < 0 ) call utl_abort('bcc_readUABcorStype ERROR:  Unknown sonde type '//stype)
      if ( ilatband < 0 .or. ilatband > uaNlatBands ) then
        write(*,*) 'stype, latband = ',  stype, ilatband
        call utl_abort('bcc_readUABcorStype ERROR: Latitude band out of range 0-5!')
      end if
      ! Skip stype/latband 0 (globe) if present
      if ( ilatband == 0 ) then
        do j = 1,nGroups
          do levelIndex = 1,nMandLevs
            read (nulcoeff, *, iostat=ierr) TTbiasNight, TTbiasDay, TDbiasNight, TDbiasDay
          end do
        end do
        cycle main_loop
      end if
      sondeTypes(istype) = stype
      do j = 1,nGroups
        jj = (j-1)*2 + 1
        do levelIndex = 1,nMandLevs
          read (nulcoeff, *, iostat=ierr) TTbiasNight, TTbiasDay, TDbiasNight, TDbiasDay
          if ( ierr /= 0 ) then
            call utl_abort('bcc_readUABcorStype: ERROR reading corrections in radiosonde bias correction file ' // biasCorrectionFileName )
          end if
          if ( TTbiasNight == uaMissingValue ) TTbiasNight = MPC_missingValue_R8
          if ( TTbiasDay == uaMissingValue )   TTbiasDay   = MPC_missingValue_R8
          if ( TDbiasNight == uaMissingValue ) TDbiasNight = MPC_missingValue_R8
          if ( TDbiasDay == uaMissingValue )   TDbiasDay   = MPC_missingValue_R8
          if ( TTbiasNight /= MPC_missingValue_R8 ) ttCorrections(istype,ilatband,jj,levelIndex)   = -1.0d0*TTbiasNight
          if ( TTbiasDay /= MPC_missingValue_R8 )   ttCorrections(istype,ilatband,jj+1,levelIndex) = -1.0d0*TTbiasDay
          if ( TDbiasNight /= MPC_missingValue_R8 ) tdCorrections(istype,ilatband,jj,levelIndex)   = -1.0d0*TDbiasNight
          if ( TDbiasDay /= MPC_missingValue_R8 )   tdCorrections(istype,ilatband,jj+1,levelIndex) = -1.0d0*TDbiasDay
        end do
        if ( ttCorrections(i,ilatband,jj,Index500mb)   /= MPC_missingValue_R8 ) BiasCorrPresentStype(istype,ilatband,jj)   = .true.
        if ( ttCorrections(i,ilatband,jj+1,Index500mb) /= MPC_missingValue_R8 ) BiasCorrPresentStype(istype,ilatband,jj+1) = .true.
      end do
      
    end do main_loop
    
    ierr = fclos(nulcoeff)

  end subroutine bcc_readUABcorStype
  
  !-----------------------------------------------------------------------
  ! bcc_readUABcorStn
  !-----------------------------------------------------------------------
  subroutine bcc_readUABcorStn(biasCorrectionFileName,nProfsMin,nGroups)
    !
    ! :Purpose: Read TT, TD biases by STATION/sonde-type on 16 mandatory levels for UA family. 
    !           The first line is the station and sonde-type followed by number of profiles. 
    !           Station = 'END' for end of file.
    !           nGroups groups of of 16 lines follow, e.g., nGroups = 2: one group for ascending sonde observations and 
    !                                                                    one for descending sonde observations 
    !               -- 16 lines are 16 mandatory levels from 100000 Pa to 1000 Pa
    !               -- Lines contain 4 values: TT night-bias, TT day-bias, TD night-bias, TD day-bias
    !           Missing value (file) = -99.0.
    !           Sonde type in file is "None" if type was unknown/missing for the reports from a station.
    !           BiasCorrPresentStn(iStn,iType,TOD) = .true. for all TOD if total Nprofs >= nProfsMin (e.g. 100)
    !           ttCorrectionsStn(iStn,iType,TOD,level) = MPC_missingValue_R8 if TTcorrectionValue == -99.0
    !
    implicit none
    !Arguments:
    character(len=*), intent(in) :: biasCorrectionFileName
    integer, intent(in)          :: nProfsMin, nGroups

    !Locals:
    integer :: ierr, nulcoeff, iStn, iType, N
    integer :: i, j, levelIndex, nProfs, jj, jmax
    real(8) :: TTbiasNight, TTbiasDay, TDbiasNight, TDbiasDay
    character(len=8) :: stype
    character(len=9) :: stn, stn_prev

    if ( uaRevOnly ) return
    
    nulcoeff = 0
    ierr = fnom(nulcoeff, biasCorrectionFileName, 'FMT+R/O', 0)
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readUABcorStn: unable to open radiosonde bias correction file ' // biasCorrectionFileName )
    end if
    
    jmax = nGroups*2

    allocate( ttCorrectionsStn(nStationMaxUA,nSondesMax,jmax,nMandLevs) )
    allocate( tdCorrectionsStn(nStationMaxUA,nSondesMax,jmax,nMandLevs) )
    allocate( uaStations(nStationMaxUA) )
    allocate( BiasCorrPresentStn(nStationMaxUA,nSondesMax,jmax) )

    ttCorrectionsStn(:,:,:,:) =  MPC_missingValue_R8
    tdCorrectionsStn(:,:,:,:) =  MPC_missingValue_R8
    uaStations(:) = 'empty'
    BiasCorrPresentStn(:,:,:) = .false.
   
    iStn = 0
    stn_prev = "toto"
    main_loop: do 
      
      read (nulcoeff, *, iostat=ierr) stn, stype, nProfs
      if ( ierr /= 0 ) then
        call utl_abort('bcc_readUABcorStn: ERROR reading stn line in radiosonde bias correction file ' // biasCorrectionFileName )
      end if
      if ( trim(stn) == 'END' ) exit main_loop
      iType = bcc_SondeIndex(stype)
      if ( iType < 0 ) call utl_abort('bcc_readUABcorStn ERROR:  Unknown sonde type '//stype//' not in rs_types defined in namelist')
      
      ! If new station, add station to uaStations array.
      ! Assumes entries for the same station and different types are consecutive in bcor file
      if ( trim(stn) /= trim(stn_prev) ) then
        iStn = iStn + 1
        uaStations(iStn) = stn
        stn_prev = stn
      end if      
      
      ! Determine which station/stype/time-of-day have enough profiles to use the bias corrections
      ! and set logical BiasCorrPresentStn accordingly.
      
      do j = 1,jmax
        if ( nProfs >= nProfsMin ) BiasCorrPresentStn(iStn,iType,j) = .true.
      end do

      ! Read the bias correction profiles for this station/stype. 
      do j = 1,nGroups
        jj = (j-1)*2 + 1
        do levelIndex = 1,nMandLevs
          read (nulcoeff, *, iostat=ierr) TTbiasNight, TTbiasDay, TDbiasNight, TDbiasDay
          if ( ierr /= 0 ) then
            call utl_abort('bcc_readUABcorStn: ERROR reading corrections in radiosonde bias correction file ' // biasCorrectionFileName )
          end if
          if ( TTbiasNight == uaMissingValue ) TTbiasNight = MPC_missingValue_R8
          if ( TTbiasDay == uaMissingValue )   TTbiasDay   = MPC_missingValue_R8
          if ( TDbiasNight == uaMissingValue ) TDbiasNight = MPC_missingValue_R8
          if ( TDbiasDay == uaMissingValue )   TDbiasDay   = MPC_missingValue_R8
          if ( TTbiasNight /= MPC_missingValue_R8 ) ttCorrectionsStn(iStn,iType,jj,levelIndex)   = -1.0d0*TTbiasNight
          if ( TTbiasDay /= MPC_missingValue_R8 )   ttCorrectionsStn(iStn,iType,jj+1,levelIndex) = -1.0d0*TTbiasDay
          if ( TDbiasNight /= MPC_missingValue_R8 ) tdCorrectionsStn(iStn,iType,jj,levelIndex)   = -1.0d0*TDbiasNight
          if ( TDbiasDay /= MPC_missingValue_R8 )   tdCorrectionsStn(iStn,iType,jj+1,levelIndex) = -1.0d0*TDbiasDay
        end do
      end do

    end do main_loop
    
    ierr = fclos(nulcoeff)

  end subroutine bcc_readUABcorStn
  
  !-----------------------------------------------------------------------
  ! bcc_applyUABcor
  !-----------------------------------------------------------------------
  subroutine bcc_applyUABcor(obsSpaceData)
    !
    ! :Purpose:  To apply bias corrections to radiosonde TT and ES observationa in obsSpaceData 
    !
    !  NOTE: We ADD the correction (negative of the bias) to the raw observation.
    !        We SUBTRACT the old correction from the corrected value to remove (reverse) the old correction!
    !        Value put in OBS_BCOR column is the corection that was ADDED to the raw observation.
    !
    implicit none
    !Arguments:
    type(struct_obs)        :: obsSpaceData
    !Locals:
    integer  :: headerIndex, bodyIndex, codtyp
    integer  :: flag, bufrCode, stype_code, stype_index, stn_index
    integer  :: date, time, latband, jc, phase
    integer  :: countTTCorrections,  countESCorrections, countUnknownStype, countMissingStype
    integer  :: countUnknownStation, countTTCorrByStation, countTTCorrByStype, countTTObs, countRS41
    integer  :: countCat1, countCat2, countNight, countDay, countDD
    integer  :: TTbodyIndex, ESbodyIndex
    real(8)  :: tt, oldCorr, pressure, lat, lon, tt_orig, es, es_orig, td
    real(8)  :: solar_elev, corr, p1, p2, p3, p4
    real(8)  :: TimeOfDayX
    character(len=9)  :: stnid, stnid_prev
    character(len=8)  :: stype
    character(len=10) :: TimeOfDay
    character(len=5)  :: sourceCorr
    logical  :: newStation, SondeTypeFound, debug, StationFound, RealRS41

    if ( .not. uaBiasActive ) return

    write(*,*) "bcc_applyUABcor: start"
    
    debug = .false.
    
    debug = debug .and. ( mpi_myid == 0 )

    ! Read the ascii files containing TT and TD bis profiles by sonde-type and by station
    if ( .not.uaRevOnly ) then
      call bcc_readUABcorStype(uaBcFileStype, uaNbiasCat)
      call bcc_readUABcorStn(uaBcFileStn, uaNprofsMin, uaNbiasCat)
    end if
    
    countTTCorrections    = 0
    countESCorrections    = 0
    countUnknownStype     = 0
    countMissingStype     = 0
    countUnknownStation   = 0
    countTTCorrByStation  = 0
    countTTCorrByStype    = 0
    countTTObs            = 0
    countRS41             = 0
    countCat1             = 0
    countCat2             = 0
    countNight            = 0
    countDay              = 0
    countDD               = 0

    call obs_set_current_header_list(obsSpaceData,'UA')
    
    stnid_prev = 'toto'
    HEADER: do           ! Loop over UA headers (one header per pressure level)
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER
      
      call obs_set_current_body_list(obsSpaceData, headerIndex)

      stnid  = trim(obs_elem_c(obsSpaceData,'STID',headerIndex))
      if ( trim(stnid) /= trim(stnid_prev) ) then
         newStation = .true.
         stnid_prev = stnid
      else 
         newStation = .false.
      end if
      
      ! Get the information needed to apply the bias corrections from the first header for this station
      if ( newStation .and. .not.uaRevOnly ) then
         
         !! Station index in list of stations (uaStations) from ua_bcors_stn file
         StationFound = .false.
         stn_index = bcc_StationIndex(stnid)
         if (debug) write(*,*) 'stnid, index = ', stnid, stn_index
         if ( stn_index > 0 ) then
           StationFound = .true.
         else 
           countUnknownStation = countUnknownStation + 1
           if (debug) write(*,*) 'Unknown station (not in ua_bcors_stn file) '//stnid
         end if
         
         !! Date and lat,lon
         codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex) ! CODE TYPE
         date   = obs_headElem_i(obsSpaceData, OBS_DAT, headerIndex) ! YYYYMMDD
         time   = obs_headElem_i(obsSpaceData, OBS_ETM, headerIndex) ! HHMM
         lat    = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex) ! radians!
         lon    = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex) ! radians!
         
         ! Sonde phase: 1=asc, 2=desc
         phase  = bcc_uaPhase(codtyp)
         
         if (uaNbiasCat == 2) then
            jc = (phase-1)*2 + 1
         else
            jc = 1
         end if
         
         if ( jc == 1 ) then 
            countCat1 = countCat1 + 1
         else
            countCat2 = countCat2 + 1
         end if
         
         !! Sonde type index in rs_types read from namelist
         SondeTypeFound = .false.
         stype_code  = obs_headElem_i(obsSpaceData, OBS_RTP, headerIndex)  ! sonde type BUFR code (WMO table)
         if (debug) then
            write(*,*) ' '
            write(*,*) stnid, stype_code, date, time, lat*MPC_DEGREES_PER_RADIAN_R8
         end if
         if ( stype_code == 0 ) then
            countMissingStype = countMissingStype + 1
            write (*,'(a60)') "bcc_applyUABcor: Missing sonde type at stn "//stnid
         end if
         call bcc_GetSondeType(stype_code,stype,stype_index)
         if ( stype_index == 0 ) then
            countUnknownStype = countUnknownStype + 1
            write (*,'(a70,i4)') "bcc_applyUABcor: Unknown sonde-type code at stn "//stnid//". Code = ", stype_code
         else 
            SondeTypeFound = .true.
            if (debug) write(*,*) 'stype, index = ', stype, stype_index
         end if
         
         RealRS41 = .false.
         if ( trim(stype) == "RS41" ) then
           if ( stnid(1:1) == "5" ) then
             RealRS41 = .false.
           else
             RealRS41 = .true.
           end if
         end if

         !! Time-of-day x-value
         call bcc_GetSolarElevation(lat,lon,date,time,solar_elev)
         call bcc_GetTimeOfDay(solar_elev,TimeOfDayX)
         if ( TimeOfDayX == 0.0 ) then
            countNight = countNight+1
         elseif ( TimeOfDayX == 1.0 ) then
            countDay = countDay+1
         else
            countDD = countDD+1
         end if
         
         !! Latitude band
         if ( uaNlatBands == 1 ) then
            latband = 1
         else
            latband = bcc_LatBand(lat)
            if ( latband == -1 ) then
               write(*,*) "uaNlatBands = ", uaNlatBands
               call utl_abort("bcc_applyUABcor ERROR: uaNlatBands must equal 5")
            end if
         end if

         if (debug) write(*,*) solar_elev, TimeOfDayX, latband
         
      end if
      
      TTbodyIndex = -1
      ESbodyIndex = -1
      
      BODY: do   ! Find the body indices for the TT and ES observations for this header (pressure level)

        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY 

        bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex )

        ! Assumes that normal (non high-precision) codes are used for TT and ES observations in obsSpaceData
        if ( bufrCode == BUFR_NETT ) then
           TTbodyIndex = bodyIndex
        elseif ( bufrCode == BUFR_NEES ) then
           ESbodyIndex = bodyIndex
        end if
        
      end do BODY
      
      ! Reverse any old corrections and apply new corrections (if uaRevOnly=.false.)
      
      ! TT bias correction
      if ( TTbodyIndex >= 0 ) then
         bodyIndex = TTbodyIndex
         tt_orig  = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
         tt       = tt_orig
         pressure = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex ) * MPC_MBAR_PER_PA_R8
         flag     = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
         oldCorr  = obs_bodyElem_r(obsSpaceData, OBS_BCOR, bodyIndex )
         corr     = MPC_missingValue_R8
         if ( debug .and. pressure == 500.0d0 ) write(*,*) 'TT',TTbodyIndex,pressure,tt
         if ( tt /= MPC_missingValue_R8 ) then
            if ( btest(flag, 6) .and. oldCorr /= MPC_missingValue_R8 ) then
               tt = tt - oldCorr
               flag = ibclr(flag, 6)
            end if
            if ( uaRevOnly ) corr = 0.0d0
            if ( .not.uaRevOnly ) then
              ! Get the TT correction for this station, sonde type, time-of-day and level
              ! stn_index   = -1 if station not found in bcor file
              ! stype_index =  0 if sonde-type code is not associated with any of the types in
              !                namelist (rs_types) including stypes "Others" and "None" [code=0]
              countTTObs = countTTObs + 1
              if ( RealRS41 ) countRS41 = countRS41 + 1
              call bcc_GetUACorrection('TT',stn_index,stype_index,stype,jc,TimeOfDayX,latband,pressure,corr,sourceCorr)
              if ( sourceCorr == 'stn' ) then
                countTTCorrByStation = countTTCorrByStation + 1
              elseif ( sourceCorr == 'stype' ) then
                countTTCorrByStype = countTTCorrByStype + 1
              end if
              if ( debug .and. pressure == 500.0d0 ) write(*,*) 'TT corr, source, jc = ',corr, sourceCorr, jc
              if ( corr /= MPC_missingValue_R8 ) then
                 tt = tt + corr
                 flag = ibset(flag, 6)
                 countTTCorrections = countTTCorrections + 1
              else
                 if ( .not.RealRS41 .and. uaRejUnBcor ) flag = ibset(flag, 11)
              end if
            end if
         end if
         call obs_bodySet_r( obsSpaceData, OBS_BCOR, bodyIndex, corr )
         call obs_bodySet_r( obsSpaceData, OBS_VAR , bodyIndex, tt   )
         call obs_bodySet_i( obsSpaceData, OBS_FLG , bodyIndex, flag ) 
      end if
      
      ! ES bias correction
      if ( ESbodyIndex >= 0 .and. TTbodyIndex >= 0 ) then
         bodyIndex = ESbodyIndex
         es_orig  = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
         es       = es_orig
         pressure = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex ) * MPC_MBAR_PER_PA_R8
         flag     = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
         oldCorr  = obs_bodyElem_r(obsSpaceData, OBS_BCOR, bodyIndex )
         corr     = MPC_missingValue_R8
         if ( debug .and. pressure == 500.0d0 ) write(*,*) 'ES',ESbodyIndex,pressure,es
         if ( es /= MPC_missingValue_R8 ) then
            if ( btest(flag, 6) .and. oldCorr /= MPC_missingValue_R8 ) then
               es = es - oldCorr
               flag = ibclr(flag, 6)
            end if
            if ( uaRevOnly ) corr = 0.0d0
            if ( .not.uaRevOnly ) then
              call bcc_GetUACorrection('TD',stn_index,stype_index,stype,jc,TimeOfDayX,latband,pressure,corr,sourceCorr)
              if ( debug .and. pressure == 500.0d0 ) write(*,*) 'TD corr, source, jc = ',corr, sourceCorr, jc
              if ( corr /= MPC_missingValue_R8 ) then
                 td = (tt_orig - es) + corr
                 es = tt - td
                 if ( es < 0.0 ) es =  0.0d0
                 if ( es > 30.0) es = 30.0d0
                 corr =  es - es_orig
                 if ( debug .and. pressure == 500.0d0 ) write(*,*) 'ES corr =',corr
                 flag = ibset(flag, 6)
                 countESCorrections = countESCorrections + 1
              else
                 if ( .not.RealRS41 .and. uaRejUnBcor ) flag = ibset(flag, 11)
              end if
            end if
         end if
         call obs_bodySet_r( obsSpaceData, OBS_BCOR, bodyIndex, corr )
         call obs_bodySet_r( obsSpaceData, OBS_VAR , bodyIndex, es   )
         call obs_bodySet_i( obsSpaceData, OBS_FLG , bodyIndex, flag ) 
      end if     
              
    end do HEADER
    
    if ( .not.uaRevOnly ) then
       write (*, '(a60, i10)' ) "bcc_applyUABcor: Number of TT obs corrected          = ", countTTCorrections
       write (*, '(a60, i10)' ) "bcc_applyUABcor: Number of ES obs corrected          = ", countESCorrections
       write (*, '(a60, i10)' ) "bcc_applyUABcor: Number obs with unknown station     = ", countUnknownStation
       write (*, '(a60, i10)' ) "bcc_applyUABcor: Number obs with missing sonde type  = ", countMissingStype
       write (*, '(a60, i10)' ) "bcc_applyUABcor: Number obs with unknown sonde code  = ", countUnknownStype
       write (*, '(a40)'      ) "----------------------------------------"
       write (*, '(a60, i10)' ) "bcc_applyUABcor: Total number of TT observations     = ", countTTObs
       if ( countTTObs > 0 ) then
         p1 = (float(countTTCorrByStation)/float(countTTObs))*100.d0
         p2 = (float(countTTCorrByStype)/float(countTTObs))*100.d0
         p3 = 100.d0 - (p1+p2)
         p4 = (float(countRS41)/float(countTTObs))*100.d0
         write (*, '(a60, f7.2)' ) "bcc_applyUABcor: Percentage corrected by STATION     = ", p1
         write (*, '(a60, f7.2)' ) "bcc_applyUABcor: Percentage corrected by SONDE-TYPE  = ", p2
         write (*, '(a60, f7.2)' ) "bcc_applyUABcor: Percentage uncorrected              = ", p3
         write (*, '(a60, f7.2)' ) "bcc_applyUABcor: Percentage RS41                     = ", p4
         write (*, '(a60, i10)' )  "bcc_applyUABcor: Number of night profiles            = ", countNight
         write (*, '(a60, i10)' )  "bcc_applyUABcor: Number of day profiles              = ", countDay
         write (*, '(a60, i10)' )  "bcc_applyUABcor: Number of dawn-dusk profiles        = ", countDD
         if (uaNbiasCat == 2) then
           write (*, '(a60, i10)' )  "bcc_applyUABcor: Number of ascent profiles            = ", countCat1
           write (*, '(a60, i10)' )  "bcc_applyUABcor: Number of descent profiles           = ", countCat2
         end if
       end if
    end if
    
    if ( allocated(ttCorrections) )        deallocate(ttCorrections)
    if ( allocated(tdCorrections) )        deallocate(tdCorrections)
    if ( allocated(ttCorrectionsStn) )     deallocate(ttCorrectionsStn)
    if ( allocated(tdCorrectionsStn) )     deallocate(tdCorrectionsStn)    
    if ( allocated(BiasCorrPresentStype) ) deallocate(BiasCorrPresentStype)
    if ( allocated(BiasCorrPresentStn) )   deallocate(BiasCorrPresentStn)
    if ( allocated(uaStations) )           deallocate(uaStations)
    
    if ( allocated(rs_types) ) deallocate(rs_types)
    
    write(*,*) "bcc_applyUABcor: end"
    
  end subroutine bcc_applyUABcor
  
  
  !-----------------------------------------------------------------------
  ! bcc_biasActive
  !-----------------------------------------------------------------------
  logical function bcc_biasActive(obsFam)
    !
    ! :Purpose: returns True if bias correction is active for the given conventional observation family
    !
    ! Used in subroutine obsf_readFiles to see if bias correction element needs to be added to derialt files.
    ! AI derialt files aleady have bias correction element added so AI family is not included here.
    !
    implicit none
    character(len=*),intent(in) :: obsFam

    if (.not.initialized) call bcc_readConfig()
    
    select case(trim(obsFam))
    case('GP')
      bcc_biasActive = bcc_gpBiasActive
    case('UA')
      bcc_biasActive = bcc_uaBiasActive
    case default
      bcc_biasActive = .false.
    end select

  end function bcc_biasActive


end MODULE biasCorrectionConv_mod


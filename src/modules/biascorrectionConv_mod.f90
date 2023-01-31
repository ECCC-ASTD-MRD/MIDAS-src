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
  
  ! This variable is set to .true. when bcc_readConfig() [the routine that reads &NAMBIASCONV namelist]
  ! is called for the first time to read/initialize the bias correction namelist variables.
  logical            :: initialized = .false.

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
  
  integer               :: nbAircrafts, nbGpStations
  real(8), allocatable  :: AIttCorrections(:,:,:)
  real(8), allocatable  :: ztdCorrections(:)
  real(8), allocatable  :: ttCorrections(:,:,:,:),    tdCorrections(:,:,:,:)
  real(8), allocatable  :: ttCorrectionsStn(:,:,:,:), tdCorrectionsStn(:,:,:,:)
  
  character(len=9), allocatable     :: aircraftIds(:), gpsStations(:), uaStations(:)
  character(len=8) :: sondeTypes(nSondesMax)

  real(8) :: mandLevs(nMandLevs), tolPress(nMandLevs)
  
  logical :: bcc_aiBiasActive, bcc_gpBiasActive, bcc_uaBiasActive
  logical, allocatable :: biasCorrPresentStype(:,:,:), biasCorrPresentStn(:,:,:)
  
  ! Bias correction files (must be in program working directory)
  character(len=8), parameter  :: aiBcFile = "ai_bcors", gpBcFile = "gp_bcors"
  character(len=14), parameter :: uaBcFileStype = "ua_bcors_stype", uaBcFileStn = "ua_bcors_stn"

  integer, external    :: fnom, fclos, newdate
  
  ! Namelist variables
  logical :: aiBiasActive ! Control if bias correction is applied to aircraft data
  logical :: gpBiasActive ! Control if bias correction is applied to ground-based GPS data
  logical :: uaBiasActive ! Control if bias correction is applied to radiosonde data
  logical :: aiRevOnly    ! Don't apply new correction but simply reverse any old corrections for AI
  logical :: gpRevOnly    ! Don't apply new correction but simply reverse any old corrections for GP
  logical :: uaRevOnly    ! Don't apply new correction but simply reverse any old corrections for UA
  logical :: uaRejUnBcor  ! Set DATA QC flag bit 11 on to exclude uncorrected UA observations from assimilation 
  integer :: uaNbiasCat   ! Number of bias profile categories in UA bcor files, e.g. 1, or 2 for "asc" and "desc" phase categories
  integer :: uaNlatBands  ! Number of latitude bands in ua_bcors_stype bcor file (= 5 or 1). Set to 1 if there are no latitude bands in file
  integer :: uaNprofsMin  ! Min number of bias profiles required for a station/stype/time-of-day to use biases 'ua_bcors_stn' as corrections
  character(len=8) :: nlSondeTypes(nSondesMax)    ! List of radiosonde type names
  integer          :: nlSondeCodes(nSondesMax,20) ! List of radiosonde type codes
  integer          :: nlNbSondes                  ! Number of radiosonde types in lists

  ! Structure to hold dictionary containing the BUFR sonde type codes associated with each sonde type
  type sondeType
    character(len=8)  :: name
    integer           :: codes(20)
  end type sondeType
  
  type(sondeType), allocatable  :: rsTypes(:)
  
  namelist /nambiasconv/ aiBiasActive,gpBiasActive,aiRevOnly,gpRevOnly,uaBiasActive,uaRevOnly,uaNprofsMin,uaRejUnBcor,uaNbiasCat,uaNlatBands
  namelist /namsondetypes/ nlNbSondes, nlSondeTypes, nlSondeCodes
  
  ! 16 mandatory pressure levels (mb) on which radiosonde bias profiles are defined
  data mandLevs /1000.d0,925.d0,850.d0,700.d0,500.d0,400.d0,300.d0,250.d0,200.d0,150.d0,100.d0,70.d0,50.d0,30.d0,20.d0,10.d0/
  ! +/- tolerance (mb) for matching a radiosonde observation pressure to one of the mandatory levels
  data tolPress /  10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0,  5.d0,  5.d0,  5.d0,  5.d0, 2.d0, 2.d0, 2.d0, 1.d0, 1.d0/
  
CONTAINS

  !-----------------------------------------------------------------------
  ! bcc_readConfig
  !-----------------------------------------------------------------------
  subroutine bcc_readConfig()
    !
    ! :Purpose: Read NAMBIASCONV namelist section and NAMSONDETYPES section if uaBiasActive=true
    !
    implicit none

    !Locals:
    integer  :: ierr, nulnam, sondeIndex
    
    if ( initialized ) then
      write(*,*) "bcc_readConfig has already been called. Returning..."
      return
    end if
  
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
      if ( mmpi_myid == 0 ) write(*,nml=nambiasconv)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'bcc_readconfig: NAMBIASCONV section is missing in the namelist. The default values will be used.'
    end if
        
    bcc_aiBiasActive = aiBiasActive
    bcc_gpBiasActive = gpBiasActive
    bcc_uaBiasActive = uaBiasActive
    
    ! read in the namelist NAMSONDETYPES
    if ( uaBiasActive .and. .not.uaRevOnly ) then
      if ( utl_isNamelistPresent('namsondetypes','./flnml') ) then
        nlNbSondes = MPC_missingValue_INT
        nlSondeTypes(:)   = 'empty'
        nlSondeCodes(:,:) = MPC_missingValue_INT
        nulnam = 0
        ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
        read(nulnam,nml=namsondetypes,iostat=ierr)
        if ( ierr /= 0 )  call utl_abort('bcc_readConfig: Error reading namelist section NAMSONDETYPES')
        if (nlNbSondes /= MPC_missingValue_INT) then
          call utl_abort('bcc_readConfig: check namsondetypes namelist section, you should remove nlNbSondes')
        end if
        nlNbSondes = 0
        do sondeIndex = 1, nSondesMax
          if (trim(nlSondeTypes(sondeIndex)) == 'empty') exit
          nlNbSondes = nlNbSondes + 1
        end do
        if ( mmpi_myid == 0 ) write(*,nml=namsondetypes)
        ierr = fclos(nulnam)
      else 
        write(*,*)
        call utl_abort('bcc_readconfig: NAMSONDETYPES section is missing in the namelist!')
      end if
      if ( nlNbSondes > nSondesMax ) call utl_abort('bcc_readconfig: Number of sonde types in namelist exceeds nSondesMax!')
      if ( nlNbSondes <= 0 )         call utl_abort('bcc_readconfig: Number of sonde types in namelist <= 0!')
      allocate( rsTypes(nlNbSondes) )
      do sondeIndex = 1, nlNbSondes
        rsTypes(sondeIndex)%name     = nlSondeTypes(sondeIndex)
        rsTypes(sondeIndex)%codes(:) = nlSondeCodes(sondeIndex,:)
      end do
    end if
    
    initialized = .true.
    
  end subroutine bcc_readConfig  
  
  !-----------------------------------------------------------------------
  ! bcc_UACorrection
  !-----------------------------------------------------------------------
  function bcc_UACorrection(timeOfDayX,corrNight,corrDay) result(uaCorrection)
    ! :Purpose: Returns a UA bias correction from day and night corrections
    
    implicit none
    
    real(8) :: uaCorrection
    
    !Arguments:
    real(8), intent(in)  :: timeOfDayX ! 0.0 (night) <= timeOfDayX <= 1.0 (day), depends on solar_elev
    real(8), intent(in)  :: corrNight  ! night bias correction
    real(8), intent(in)  :: corrDay    ! day bias correction
    
    uaCorrection = MPC_missingValue_R8
    
    if ( timeOfDayX == 0.0d0 ) then
      uaCorrection = corrNight
    else if ( timeOfDayX == 1.0d0 ) then
      uaCorrection = corrDay
    else
      if ( corrNight /= MPC_missingValue_R8 .and. corrDay /= MPC_missingValue_R8 ) then
        uaCorrection = timeOfDayX*corrDay + (1.0d0-timeOfDayX)*corrNight
      end if
    end if
    
  end function bcc_UACorrection
  
  !-----------------------------------------------------------------------
  ! bcc_GetUACorrection
  !-----------------------------------------------------------------------
  subroutine bcc_GetUACorrection(varName,stnIndex,sondeTypeIndex,sondeType,biasProfileCategory,timeOfDayX,latband,obsPressure,corr,sourceCorr)
    !
    ! :Purpose: Return a TT or TD bias correction (corr) 
    
    !
    !      varName        = (str) 'TT' or 'TD'
    !      stnIndex       = (int) station index (in uaStations) = -1 if station was not found in bcor file
    !      sondeTypeIndex = (int) sonde type index (in rsTypes) =  0 if sonde-type code is not associated with any type
    !      sondeType      = (str) sonde-type from rsTypes
    !      biasProfileCategory = (int) bias profile category = 1 (ascent/none) or 3 (descent)
    !      timeOfDayX     = (float)  0.0 (night) <= TimeOfDayX <= 1.0 (day), depends on solar_elev
    !      latband        = (int) 1-5
    !      obsPressure    = (float) pressure level (mb) of observation
    !
    !  Also returns info regarding the source of the correction returned:
    !      sourceCorr     = (str) 'none', 'stn', or 'stype'
    !
    ! Requires TT and TD correction profiles read from UA bcor file ua_bcors_stn and stored in
    !    ttCorrectionsStn(stnIndex,sondeTypeIndex,j,MandLev)
    !    tdCorrectionsStn(stnIndex,sondeTypeIndex,j,MandLev)
    ! with backup corrections by sonde-type read from UA bcor file ua_bcors_stype and stored in
    !    ttCorrections(sondeTypeIndex,latband,j,MandLev)
    !    tdCorrections(sondeTypeIndex,latband,j,MandLev)
    ! where MandLev = 1 (1000 mb) to 16 (10 mb)
    ! If uaNbiasCat = 1 (biasProfileCategory = 1)
    !       j = 1 (night)
    !           2 (day)
    ! If uaNbiasCat = 2 (biasProfileCategory = 1 [ascent] or 3 [descent])
    !       j = 1 (night-ascent)
    !           2 (day-ascent)
    !           3 (night-descent)
    !           4 (day-descent)
    !
    ! Interpolation of the correction profiles on mandatory levels is used to get the 
    ! correction at the observation level (obsPressure).
    ! Persistence is applied for observations outside the range of the mandatory levels.
    !
    implicit none

    !Arguments:
    integer, intent(in)           ::  stnIndex
    integer, intent(in)           ::  sondeTypeIndex
    integer, intent(in)           ::  biasProfileCategory
    integer, intent(in)           ::  latband
    real(8), intent(in)           ::  timeOfDayX
    real(8), intent(in)           ::  obsPressure
    character(len=*), intent(in)  ::  varName
    character(len=*), intent(in)  ::  sondeType
    real(8), intent(out)          ::  corr
    character(len=*), intent(out) ::  sourceCorr

    !Locals:
    real(8) ::  corrProfileStnDay(16), corrProfileStnNight(16), corrProfileStypeDay(16), corrProfileStypeNight(16)
    real(8) ::  corrDay, corrNight
    real(8) ::  weight1, weight2, deltaPressure, pressureAbove, pressureBelow, corrAbove, corrBelow
    real(8) ::  corrAboveNight, corrBelowNight, corrAboveDay, corrBelowDay
    integer ::  level, levelIndex
    logical ::  profileExistsStn, profileExistsStype, doInterp

    sourceCorr = "none"
    
    ! Bias correction by station and sonde-type
    if ( stnIndex == -1 ) then
      profileExistsStn = .false.
    else
      if ( sondeTypeIndex > 0 ) then
        ! Check that correction profile for this station, sonde-type and TimeOfDay is available for use
        profileExistsStn = biasCorrPresentStn(stnIndex,sondeTypeIndex,biasProfileCategory)
      else ! There are no bias profiles for this station/stype combination
        profileExistsStn = .false.
      end if
    end if

    ! BACKUP: Bias correction by sonde-type
    profileExistsStype = .false.
    if ( trim(sondeType) /= 'unknown' .and. trim(sondeType) /= 'Others' .and. trim(sondeType) /= 'None' ) then
      ! Check if correction profile for this sonde-type,latband,and TimeOfDay is available for use
      profileExistsStype = biasCorrPresentStype(sondeTypeIndex,latband,biasProfileCategory)
    end if
       
    if ( .not.profileExistsStn  .and. .not.profileExistsStype ) then
      corr = MPC_missingValue_R8
      return
    end if
    
    ! Fill the night and day bias correction profiles
    if ( varName == 'TT' ) then
      if ( profileExistsStn ) then
        corrProfileStnNight(:)   = ttCorrectionsStn(stnIndex,sondeTypeIndex,biasProfileCategory,:)
        corrProfileStnDay(:)     = ttCorrectionsStn(stnIndex,sondeTypeIndex,biasProfileCategory+1,:)
      end if
      if ( profileExistsStype ) then
        corrProfileStypeNight(:) = ttCorrections(sondeTypeIndex,latband,biasProfileCategory,:)
        corrProfileStypeDay(:)   = ttCorrections(sondeTypeIndex,latband,biasProfileCategory+1,:)
      end if
    else if ( varName == 'TD' ) then
      if ( profileExistsStn ) then
         corrProfileStnNight(:)   = tdCorrectionsStn(stnIndex,sondeTypeIndex,biasProfileCategory,:)
         corrProfileStnDay(:)     = tdCorrectionsStn(stnIndex,sondeTypeIndex,biasProfileCategory+1,:)
      end if
      if ( profileExistsStype ) then
         corrProfileStypeNight(:) = tdCorrections(sondeTypeIndex,latband,biasProfileCategory,:)
         corrProfileStypeDay(:)   = tdCorrections(sondeTypeIndex,latband,biasProfileCategory+1,:)
      end if
    else 
      call utl_abort('bcc_GetUACorrection: Unsupported varName '//varName)
    end if
    
    corr = MPC_missingValue_R8
    doInterp = .true.

    !--------------------------------------------------------------------------------------
    ! Get the correction at the observation level (obsPressure)
    !--------------------------------------------------------------------------------------
    
    ! Check if obsPressure is outside range of levels (no interpolation possible)
    if ( obsPressure >= mandLevs(1) ) then
      levelIndex = 1
    else if ( obsPressure <= mandLevs(nMandLevs) ) then 
      levelIndex = nMandLevs
    else
      levelIndex = 0
    end if
    
    ! Check if obsPressure is close to one of the 16 mandatory levels (no interpolation needed)
    if ( levelIndex == 0 ) then
      do level = 1, nMandLevs
        if ( abs(obsPressure-mandLevs(level)) <= tolPress(level) ) then
          levelIndex = level
          exit
        end if
      end do
    end if
    
    ! If obsPressure is close to one of the 16 mandatory levels get the correction for the mandatory level: 
    ! Use correction by station with correction by stype as backup
    if ( levelIndex > 0 ) then
      if ( profileExistsStn ) then
        corrNight = corrProfileStnNight(levelIndex)
        corrDay   = corrProfileStnDay(levelIndex)
        sourceCorr = "stn"
        corr = bcc_UACorrection(timeOfDayX,corrNight,corrDay)
        if ( corr == MPC_missingValue_R8 ) then
          sourceCorr = "none"
          if ( profileExistsStype ) then
            corrNight = corrProfileStypeNight(levelIndex)
            corrDay   = corrProfileStypeDay(levelIndex)
            sourceCorr = "stype"
            corr = bcc_UACorrection(timeOfDayX,corrNight,corrDay)
            if ( corr == MPC_missingValue_R8 ) sourceCorr = "none"
          end if
        end if
      else
        if ( profileExistsStype ) then
          corrNight = corrProfileStypeNight(levelIndex)
          corrDay   = corrProfileStypeDay(levelIndex)
          sourceCorr = "stype"
          corr = bcc_UACorrection(timeOfDayX,corrNight,corrDay)
          if ( corr == MPC_missingValue_R8 ) sourceCorr = "none"
        end if
      end if
      doInterp = .false.
    end if
    
    ! or interpolate to get correction for observation level obsPressure:
    ! Use corrections by station with corrections by stype as backup
    if ( doInterp ) then
      corrAbove = MPC_missingValue_R8
      corrBelow = MPC_missingValue_R8
      
      if ( profileExistsStn ) then
        do level = 1, nMandLevs-1
          if ( obsPressure <= mandLevs(level) .and. obsPressure > mandLevs(level+1) ) then
            pressureBelow = mandLevs(level)
            pressureAbove = mandLevs(level+1)
            corrBelowNight = corrProfileStnNight(level)
            corrAboveNight = corrProfileStnNight(level+1)
            corrBelowDay   = corrProfileStnDay(level)
            corrAboveDay   = corrProfileStnDay(level+1)
            exit
          end if
        end do
        sourceCorr = "stn"
        if ( timeOfDayX == 0.0d0 ) then
          corrBelow = corrBelowNight
          corrAbove = corrAboveNight 
        else if ( timeOfDayX == 1.0d0 ) then
          corrBelow = corrBelowDay
          corrAbove = corrAboveDay
        else 
          if ( corrBelowNight /= MPC_missingValue_R8 .and. corrBelowDay /= MPC_missingValue_R8 ) then
            corrBelow = timeOfDayX*corrBelowDay + (1.0d0-timeOfDayX)*corrBelowNight
          end if
          if ( corrAboveNight /= MPC_missingValue_R8 .and. corrAboveDay /= MPC_missingValue_R8 ) then
            corrAbove = timeOfDayX*corrAboveDay + (1.0d0-timeOfDayX)*corrAboveNight
          end if
        end if
      end if
      
      if ( corrAbove == MPC_missingValue_R8 .or. corrBelow == MPC_missingValue_R8 ) then
        if ( profileExistsStype ) then
          do level = 1, nMandLevs-1
            if ( obsPressure <= mandLevs(level) .and. obsPressure > mandLevs(level+1) ) then
              pressureBelow = mandLevs(level)
              pressureAbove = mandLevs(level+1)
              corrBelowNight = corrProfileStypeNight(level)
              corrAboveNight = corrProfileStypeNight(level+1)
              corrBelowDay   = corrProfileStypeDay(level)
              corrAboveDay   = corrProfileStypeDay(level+1)
              exit
            end if
          end do
          sourceCorr = "stype"
          if ( timeOfDayX == 0.0d0 ) then
            corrBelow = corrBelowNight
            corrAbove = corrAboveNight 
          else if ( timeOfDayX == 1.0d0 ) then
            corrBelow = corrBelowDay
            corrAbove = corrAboveDay
          else 
            if ( corrBelowNight /= MPC_missingValue_R8 .and. corrBelowDay /= MPC_missingValue_R8 ) then
              corrBelow = timeOfDayX*corrBelowDay + (1.0d0-timeOfDayX)*corrBelowNight
            end if
            if ( corrAboveNight /= MPC_missingValue_R8 .and. corrAboveDay /= MPC_missingValue_R8 ) then
              corrAbove = timeOfDayX*corrAboveDay + (1.0d0-timeOfDayX)*corrAboveNight
            end if
          end if
        end if
      end if
      deltaPressure = log10(pressureBelow)-log10(pressureAbove)
      weight1 = (log10(pressureBelow) - log10(obsPressure)) / deltaPressure
      weight2 = (log10(obsPressure) - log10(pressureAbove)) / deltaPressure
      if ( corrAbove /= MPC_missingValue_R8 .and. corrBelow /= MPC_missingValue_R8 ) then
        corr = weight1*corrAbove + weight2*corrBelow
      else if ( corrAbove /= MPC_missingValue_R8 .and. max(weight1,weight2) == weight1 ) then
        corr = corrAbove
      else if ( corrBelow /= MPC_missingValue_R8 .and. max(weight1,weight2) == weight2 ) then
        corr = corrBelow
      else
        corr = MPC_missingValue_R8
        sourceCorr = "none"
      end if
    end if
  
  end subroutine bcc_GetUACorrection
  
  !-----------------------------------------------------------------------
  ! bcc_StationIndex
  !-----------------------------------------------------------------------  
  function bcc_StationIndex(station) result(stationIndexOut)
    !
    ! :Purpose: Return the station index (order in array uaStations) corresponding to station
    !           Returns -1 if station is not found in uaStations.
    !
    implicit none
    
    integer  :: stationIndexOut
     
    !Arguments:
    character(len=*), intent(in)  :: station

    !Locals:
    integer    :: stationIndex
    
    if ( allocated(uaStations) ) then
      stationIndexOut = -1
      do stationIndex = 1, nStationMaxUA
        if ( trim(uaStations(stationIndex)) == trim(station) ) then
          stationIndexOut = stationIndex
          exit
        end if
      end do
    else 
      call utl_abort('bcc_StationIndex: array uaStations not allocated!')
    end if
    
  end function bcc_StationIndex
  
  !-----------------------------------------------------------------------
  ! bcc_SondeIndex
  !-----------------------------------------------------------------------
  function bcc_SondeIndex(sondeType) result(sondeIndex)
    !
    ! :Purpose: Return the Sonde Type index (order in array rsTypes) corresponding to the SondeType
    !           Requires array of sondeType structures (rsTypes) to be allocated and filled.
    !           Returns -1 if SondeType is not found in rsTypes.
    !
    implicit none
    
    integer  :: sondeIndex

    !Arguments:
    character(len=*), intent(in)  :: sondeType

    !Locals:
    integer    :: typeIndex, ntypes
    
    if ( allocated(rsTypes) ) then
      ntypes = nlNbSondes
    else 
      call utl_abort('bcc_SondeIndex: array rsTypes not allocated!')
    end if
    
    SondeIndex = -1
    do typeIndex = 1, ntypes
      if ( trim(rsTypes(typeIndex)%name) == trim(sondeType) ) then
        sondeIndex = typeIndex
        exit
      end if
    end do

  end function bcc_SondeIndex
  
  !-----------------------------------------------------------------------
  ! bcc_GetSondeType
  !-----------------------------------------------------------------------
  subroutine bcc_GetSondeType(sondeTypeCode,sondeType,sondeTypeIndex)
    !
    ! :Purpose: Returns the sonde type and index given a BUFR table sondeTypeCode (BUFR elem 002011)
    !           Returns sondeType='unknown', sondeTypeIndex=0 if sondeTypeCode is not found.
    !
    
    ! Requires array of sonde_type structures (rsTypes) to be allocated and filled with the 
    ! sonde type codes associated with each sonde-type (read from namelist).
    !
    implicit none

    !Arguments:
    integer, intent(in)            :: sondeTypeCode
    character(len=*), intent(out)  :: sondeType
    integer, intent(out)           :: sondeTypeIndex

    !Locals:
    integer  :: typeIndex, ntypes, sondeCode
    
    if ( allocated(rsTypes) ) then
      ntypes = nlNbSondes
    else 
      call utl_abort('bcc_GetSondeType: array rsTypes not allocated!')
    end if
    
    if ( sondeTypeCode == 190 .or. sondeTypeCode == 192 ) then      ! NCAR dropsonde (BUFR only)
      sondeCode = 13                                 ! RS92 code
    else if ( sondeTypeCode == 191 .or. sondeTypeCode == 193 ) then ! NCAR dropsonde (BUFR only)
      sondeCode = 41                                 ! RS41 code
    else if ( sondeTypeCode >=100 .and. sondeTypeCode < 200 ) then
      sondeCode = sondeTypeCode - 100
    else 
      sondeCode = sondeTypeCode
    end if
    
    sondeType = 'unknown'
    sondeTypeIndex = 0
    do typeIndex = 1, ntypes
      if ( ANY(rsTypes(typeIndex)%codes == sondeCode) ) then
        sondeType = rsTypes(typeIndex)%name 
        sondeTypeIndex = typeIndex
        exit
      end if
    end do
    
  end subroutine bcc_GetSondeType
  
  !-----------------------------------------------------------------------
  ! bcc_GetSolarElevation
  !-----------------------------------------------------------------------
  subroutine bcc_GetSolarElevation(lat,lon,date,time,solarElev)
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
    integer, intent(in)  :: date          ! yyyymmdd
    integer, intent(in)  :: time          ! hhmm
    real(8), intent(in)  :: lat           ! radians
    real(8), intent(in)  :: lon           ! radians
    real(8), intent(out) :: solarElev     ! degrees

    !Locals:
    integer :: days(13) = (/0,31,28,31,30,31,30,31,31,30,31,30,31/)
    integer :: leap_years(7) = (/2016,2020,2024,2028,2032,2036,2040/)
    integer :: yy, mmdd, mm, dd, hh, nn, doy
    real(8) :: timeUTC, timeLCL, solarDec, hourAngle, csz, sza

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
    solarDec = 0.4093d0 * sin(((2.0d0*MPC_PI_R8)/365.0d0)*(284.0d0 + float(doy)))
    hourAngle = 15.0d0*(timeLCL-12.0d0) * MPC_RADIANS_PER_DEGREE_R8
    csz = sin(solarDec)*sin(lat) + cos(solarDec)*cos(lat)*cos(hourAngle)
    sza = MPC_DEGREES_PER_RADIAN_R8 * acos(csz)
    solarElev = 90.0d0 - sza
  
  end subroutine bcc_GetSolarElevation
  
  !-----------------------------------------------------------------------
  ! bcc_GetTimeOfDay
  !-----------------------------------------------------------------------
  subroutine bcc_GetTimeOfDay(solarElev,timeOfDayX)
    !
    ! :Purpose: Returns the time-of-day x value (0.0(night) <= x <= 1.0(day))
    !
    implicit none
    !Arguments:
    real(8), intent(in)                ::  solarElev     ! degrees
    real(8), intent(out)               ::  timeOfDayX
    
    if (solarElev < -7.5d0) then 
      timeOfDayX = 0.0d0
    else if (solarElev < 22.5d0) then
      timeOfDayX = (solarElev+7.5d0)/(22.5d0+7.5d0)
    else
      timeOfDayX = 1.0d0
    end if
    
  end subroutine bcc_GetTimeOfDay
  
  !----------------------------------------------------------------------------------------
  ! bcc_uaPhase
  !----------------------------------------------------------------------------------------
  function bcc_uaPhase(codeType) result(uaPhase)
    !
    ! :Purpose: Returns the radiosonde phase (1=ascent, 2=descent) given a header code type
    !
    implicit none
    
    integer  :: uaPhase 
    
    !Arguments:
    integer, intent(in)   ::  codeType
    
    if (codeType == codtyp_get_codtyp('tempdrop')) then
      uaPhase = 2
    else
      uaPhase = 1
    end if
    
  end function bcc_uaPhase  
  
  !-----------------------------------------------------------------------
  ! bcc_LatBand
  !-----------------------------------------------------------------------
  function bcc_LatBand(latInRadians) result(latBand)
    !
    ! :Purpose: Returns latitude band number given latitude (radians)
    !
    implicit none
    
    integer :: latBand
    
    !Arguments:
    real(8), intent(in) ::  latInRadians 

    !Locals:
    real(8)             ::  latInDegrees
    
    if ( uaNlatBands /= 5 ) then
      latBand = -1
    else
      latInDegrees = MPC_DEGREES_PER_RADIAN_R8 * latInRadians
      if (latInDegrees < -90.0d0 ) latBand = 1  ! Should never be the case!
      if   (latInDegrees >= -90.0d0   .and. latInDegrees < -60.0d0) then
         latBand = 1
      else if (latInDegrees >= -60.0d0 .and. latInDegrees < -20.0d0) then
         latBand = 2
      else if (latInDegrees >= -20.0d0 .and. latInDegrees <  20.0d0) then
         latBand = 3
      else if (latInDegrees >= 20.0d0  .and. latInDegrees <  60.0d0) then
         latBand = 4
      else
         latBand = 5
      end if
    end if
    
  end function bcc_LatBand
    
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
      call utl_abort('bcc_readAIBiases: Bulk bias corrections are missing in bias correction file!' )
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
      
      ! Set header flag bit 15 on to indicate that bias correction has been ADDED to raw value.
      ! In older versions of this routine, such as the version used in IC-3/GDPS v8.0, the correction was SUBTRACTED.
      headerFlag = ibset(headerFlag, 15)
      call obs_headSet_i( obsSpaceData, OBS_ST1, headerIndex, headerFlag )
      
    end do HEADER
    
    if ( countBulkCorrections + countTailCorrections /= 0 ) then
      write (*, '(a58, i10)' ) "bcc_applyAIBcor: Number of obs with TT bulk correction  = ", countBulkCorrections
      write (*, '(a58, i10)' ) "bcc_applyAIBcor: Number of obs with TT tail correction  = ", countTailCorrections
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
    integer :: sondeTypeIndex, group, latBand, levelIndex, groupIndex, maxGroups
    real(8) :: ttBiasNight, ttBiasDay, tdBiasNight, tdBiasDay
    character(len=8) :: sondeType

    if ( uaRevOnly ) return
    
    nulcoeff = 0
    ierr = fnom(nulcoeff, biasCorrectionFileName, 'FMT+R/O', 0)
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readUABcorStype: unable to open radiosonde bias correction file ' // biasCorrectionFileName )
    end if

    maxGroups = nGroups*2
    
    allocate( ttCorrections(nSondesMax,uaNlatBands,maxGroups,nMandLevs) )
    allocate( tdCorrections(nSondesMax,uaNlatBands,maxGroups,nMandLevs) )
    allocate( biasCorrPresentStype(nSondesMax,uaNlatBands,maxGroups) )

    ttCorrections(:,:,:,:) =  MPC_missingValue_R8
    tdCorrections(:,:,:,:) =  MPC_missingValue_R8
    sondeTypes(:) = 'empty'
    
    biasCorrPresentStype(:,:,:) = .false.
   
    main_loop: do 
      
      if ( uaNlatBands == 1 ) then
        read (nulcoeff, *, iostat=ierr) sondeType
        latBand = 1
      else
        read (nulcoeff, *, iostat=ierr) sondeType, latBand
      end if
      if ( ierr /= 0 ) then
        call utl_abort('bcc_readUABcorStype: error reading sondeType, latBand in radiosonde bias correction file ' // biasCorrectionFileName )
      end if
      if ( trim(sondeType) == 'END' ) exit main_loop
      sondeTypeIndex = bcc_SondeIndex(sondeType)
      if ( sondeTypeIndex < 0 ) call utl_abort('bcc_readUABcorStype: Unknown sonde type '//sondeType)
      if ( latBand < 0 .or. latBand > uaNlatBands ) then
        write(*,*) 'sondeType, latband = ',  sondeType, latBand
        call utl_abort('bcc_readUABcorStype: Latitude band out of range 0-5!')
      end if
      ! Skip sondeType/latband 0 (globe) if present
      if ( latBand == 0 ) then
        do group = 1, nGroups
          do levelIndex = 1, nMandLevs
            read (nulcoeff, *, iostat=ierr) ttBiasNight, ttBiasDay, tdBiasNight, tdBiasDay
          end do
        end do
        cycle main_loop
      end if
      sondeTypes(sondeTypeIndex) = sondeType
      do group = 1, nGroups
        groupIndex = (group-1)*2 + 1
        do levelIndex = 1, nMandLevs
          read (nulcoeff, *, iostat=ierr) ttBiasNight, ttBiasDay, tdBiasNight, tdBiasDay
          if ( ierr /= 0 ) then
            call utl_abort('bcc_readUABcorStype: error reading corrections in radiosonde bias correction file ' // biasCorrectionFileName )
          end if
          if ( ttBiasNight == uaMissingValue ) ttBiasNight = MPC_missingValue_R8
          if ( ttBiasDay == uaMissingValue )   ttBiasDay   = MPC_missingValue_R8
          if ( tdBiasNight == uaMissingValue ) tdBiasNight = MPC_missingValue_R8
          if ( tdBiasDay == uaMissingValue )   tdBiasDay   = MPC_missingValue_R8
          if ( ttBiasNight /= MPC_missingValue_R8 ) ttCorrections(sondeTypeIndex,latBand,groupIndex,levelIndex)   = -1.0d0*ttBiasNight
          if ( ttBiasDay /= MPC_missingValue_R8 )   ttCorrections(sondeTypeIndex,latBand,groupIndex+1,levelIndex) = -1.0d0*ttBiasDay
          if ( tdBiasNight /= MPC_missingValue_R8 ) tdCorrections(sondeTypeIndex,latBand,groupIndex,levelIndex)   = -1.0d0*tdBiasNight
          if ( tdBiasDay /= MPC_missingValue_R8 )   tdCorrections(sondeTypeIndex,latBand,groupIndex+1,levelIndex) = -1.0d0*tdBiasDay
        end do
        if ( ttCorrections(sondeTypeIndex,latBand,groupIndex,Index500mb)   /= MPC_missingValue_R8 ) biasCorrPresentStype(sondeTypeIndex,latBand,groupIndex)   = .true.
        if ( ttCorrections(sondeTypeIndex,latBand,groupIndex+1,Index500mb) /= MPC_missingValue_R8 ) biasCorrPresentStype(sondeTypeIndex,latBand,groupIndex+1) = .true.
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
    !           biasCorrPresentStn(iStn,iType,TOD) = .true. for all TOD if total Nprofs >= nProfsMin (e.g. 100)
    !           ttCorrectionsStn(iStn,iType,TOD,level) = MPC_missingValue_R8 if TTcorrectionValue == -99.0
    !
    implicit none

    !Arguments:
    character(len=*), intent(in) :: biasCorrectionFileName
    integer, intent(in)          :: nProfsMin
    integer, intent(in)          :: nGroups

    !Locals:
    integer :: ierr, nulcoeff, stationIndex, typeIndex
    integer :: group, levelIndex, nProfs, groupIndex, maxGroups
    real(8) :: ttBiasNight, ttBiasDay, tdBiasNight, tdBiasDay
    character(len=8) :: sondeType
    character(len=9) :: station, stnPrev

    if ( uaRevOnly ) return
    
    nulcoeff = 0
    ierr = fnom(nulcoeff, biasCorrectionFileName, 'FMT+R/O', 0)
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readUABcorStn: unable to open radiosonde bias correction file ' // biasCorrectionFileName )
    end if
    
    maxGroups = nGroups*2

    allocate( ttCorrectionsStn(nStationMaxUA,nSondesMax,maxGroups,nMandLevs) )
    allocate( tdCorrectionsStn(nStationMaxUA,nSondesMax,maxGroups,nMandLevs) )
    allocate( uaStations(nStationMaxUA) )
    allocate( biasCorrPresentStn(nStationMaxUA,nSondesMax,maxGroups) )

    ttCorrectionsStn(:,:,:,:) =  MPC_missingValue_R8
    tdCorrectionsStn(:,:,:,:) =  MPC_missingValue_R8
    uaStations(:) = 'empty'
    biasCorrPresentStn(:,:,:) = .false.
   
    stationIndex = 0
    stnPrev = "toto"
    main_loop: do 
      
      read (nulcoeff, *, iostat=ierr) station, sondeType, nProfs
      if ( ierr /= 0 ) then
        call utl_abort('bcc_readUABcorStn: error reading station line in radiosonde bias correction file ' // biasCorrectionFileName )
      end if
      if ( trim(station) == 'END' ) exit main_loop
      typeIndex = bcc_SondeIndex(sondeType)
      if ( typeIndex < 0 ) call utl_abort('bcc_readUABcorStn: Unknown sonde type '//sondeType//' not in rsTypes defined in namelist')
      
      ! If new station, add station to uaStations array.
      ! Assumes entries for the same station and different types are consecutive in bcor file
      if ( trim(station) /= trim(stnPrev) ) then
        stationIndex = stationIndex + 1
        uaStations(stationIndex) = station
        stnPrev = station
      end if      
      
      ! Determine which station/sondeType/time-of-day have enough profiles to use the bias corrections
      ! and set logical biasCorrPresentStn accordingly.
      
      do group = 1, maxGroups
        if ( nProfs >= nProfsMin ) biasCorrPresentStn(stationIndex,typeIndex,group) = .true.
      end do

      ! Read the bias correction profiles for this station/sondeType. 
      do group = 1, nGroups
        groupIndex = (group-1)*2 + 1
        do levelIndex = 1, nMandLevs
          read (nulcoeff, *, iostat=ierr) ttBiasNight, ttBiasDay, tdBiasNight, tdBiasDay
          if ( ierr /= 0 ) then
            call utl_abort('bcc_readUABcorStn: error reading corrections in radiosonde bias correction file ' // biasCorrectionFileName )
          end if
          if ( ttBiasNight == uaMissingValue ) ttBiasNight = MPC_missingValue_R8
          if ( ttBiasDay == uaMissingValue )   ttBiasDay   = MPC_missingValue_R8
          if ( tdBiasNight == uaMissingValue ) tdBiasNight = MPC_missingValue_R8
          if ( tdBiasDay == uaMissingValue )   tdBiasDay   = MPC_missingValue_R8
          if ( ttBiasNight /= MPC_missingValue_R8 ) ttCorrectionsStn(stationIndex,typeIndex,groupIndex,levelIndex)   = -1.0d0*ttBiasNight
          if ( ttBiasDay /= MPC_missingValue_R8 )   ttCorrectionsStn(stationIndex,typeIndex,groupIndex+1,levelIndex) = -1.0d0*ttBiasDay
          if ( tdBiasNight /= MPC_missingValue_R8 ) tdCorrectionsStn(stationIndex,typeIndex,groupIndex,levelIndex)   = -1.0d0*tdBiasNight
          if ( tdBiasDay /= MPC_missingValue_R8 )   tdCorrectionsStn(stationIndex,typeIndex,groupIndex+1,levelIndex) = -1.0d0*tdBiasDay
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
    !            This public routine is called by external routines.
        
    !  NOTE: We ADD the correction (negative of the bias) to the raw observation.
    !        We first SUBTRACT the old correction (if any) from the observation to remove (reverse) the old correction!
    !        Bias correction is put in obsSpaceData OBS_BCOR column.
    !
    !  NOTE: Bias correction is NOT applied if the sonde type is RS41. Observations from this sonde type are 
    !        assumed to be unbiased.
    !
    !  Routine does nothing if uaBiasActive = .false. (just returns to calling routine)
    !
    implicit none

    !Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData

    !Locals:
    integer  :: headerIndex, bodyIndex, codtyp
    integer  :: flag, bufrCode, sondeTypeCode, sondeTypeIndex, stnIndex
    integer  :: date, time, latBand, groupIndex, phase
    integer  :: countTTCorrections,  countESCorrections, countUnknownStype, countMissingStype
    integer  :: countUnknownStation, countTTCorrByStation, countTTCorrByStype, countTTObs, countRS41
    integer  :: countCat1, countCat2, countNight, countDay, countDawnDusk
    integer  :: ttBodyIndex, esBodyIndex
    real(8)  :: tt, oldCorr, pressure, lat, lon, ttOriginal, es, esOriginal, td
    real(8)  :: solarElev, corr, p1, p2, p3, p4
    real(8)  :: timeOfDayX
    character(len=9)  :: stnid, stnidPrev
    character(len=8)  :: sondeType
    character(len=5)  :: sourceCorr
    logical  :: newStation, debug, stationFound, realRS41

    if ( .not. uaBiasActive ) return

    write(*,*) "bcc_applyUABcor: start"
    
    debug = .false.
    
    debug = debug .and. ( mmpi_myid == 0 )

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
    countDawnDusk         = 0

    call obs_set_current_header_list(obsSpaceData,'UA')
    
    stnidPrev = 'toto'
    HEADER: do           ! Loop over UA headers (one header per pressure level)
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER
      
      call obs_set_current_body_list(obsSpaceData, headerIndex)

      stnid  = trim(obs_elem_c(obsSpaceData,'STID',headerIndex))
      if ( trim(stnid) /= trim(stnidPrev) ) then
        newStation = .true.
        stnidPrev = stnid
      else 
        newStation = .false.
      end if
      
      ! Get the information needed to apply the bias corrections from the first header for this station
      if ( newStation .and. .not.uaRevOnly ) then
         
        ! Station index in list of stations (uaStations) from ua_bcors_stn file
        stationFound = .false.
        stnIndex = bcc_StationIndex(stnid)
        if (debug) write(*,*) 'stnid, index = ', stnid, stnIndex
        if ( stnIndex > 0 ) then
          stationFound = .true.
        else 
          countUnknownStation = countUnknownStation + 1
          if (debug) write(*,*) 'Unknown station (not in ua_bcors_stn file) '//stnid
        end if
         
        ! Date and lat,lon
        codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex) ! CODE TYPE
        date   = obs_headElem_i(obsSpaceData, OBS_DAT, headerIndex) ! YYYYMMDD
        time   = obs_headElem_i(obsSpaceData, OBS_ETM, headerIndex) ! HHMM
        lat    = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex) ! radians!
        lon    = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex) ! radians!
         
        ! Sonde phase: 1=asc, 2=desc
        phase  = bcc_uaPhase(codtyp)
         
        if (uaNbiasCat == 2) then
          groupIndex = (phase-1)*2 + 1
        else
          groupIndex = 1
        end if
         
        if ( groupIndex == 1 ) then 
          countCat1 = countCat1 + 1
        else
          countCat2 = countCat2 + 1
        end if
         
        ! Get the sonde type index in rsTypes read from namelist
        sondeTypeCode  = obs_headElem_i(obsSpaceData, OBS_RTP, headerIndex)  ! sonde type BUFR code (WMO table)
        if (debug) then
          write(*,*) 'stnid, sondeTypeCode, date, time, lat'
          write(*,*) stnid, sondeTypeCode, date, time, lat*MPC_DEGREES_PER_RADIAN_R8
        end if
        if ( sondeTypeCode == -999 ) then
          sondeTypeCode = 0
          countMissingStype = countMissingStype + 1
          write (*,'(a60)') "bcc_applyUABcor: Missing sonde type at stn "//stnid
        end if
        call bcc_GetsondeType(sondeTypeCode,sondeType,sondeTypeIndex)
        if ( sondeTypeIndex == 0 ) then
          countUnknownStype = countUnknownStype + 1
          write (*,'(a70,i4)') "bcc_applyUABcor: Unknown sonde-type code at stn "//stnid//". Code = ", sondeTypeCode
        end if
        if (debug) write(*,*) 'sondeType, sondeTypeIndex = ', sondeType, sondeTypeIndex
         
        ! We assume that a sonde type of "RS41" reported by Chinese UA stations is not correct
        realRS41 = .false.
        if ( trim(sondeType) == "RS41" ) then
          if ( stnid(1:1) == "5" ) then
            realRS41 = .false.
          else
            realRS41 = .true.
          end if
        end if
         
        ! Time-of-day x-value
        call bcc_GetSolarElevation(lat,lon,date,time,solarElev)
        call bcc_GetTimeOfDay(solarElev,timeOfDayX)
        if ( timeOfDayX == 0.0 ) then
          countNight = countNight+1
        else if ( timeOfDayX == 1.0 ) then
          countDay = countDay+1
        else
          countDawnDusk = countDawnDusk+1
        end if
         
        ! Latitude band
        if ( uaNlatBands == 1 ) then
          latBand = 1
        else
          latBand = bcc_LatBand(lat)
          if ( latBand == -1 ) then
            write(*,*) "uaNlatBands = ", uaNlatBands
            call utl_abort("bcc_applyUABcor: uaNlatBands must equal 1 or 5")
          end if
        end if

        if (debug) write(*,*) 'solarElev, timeOfDayX, latBand = ',solarElev, timeOfDayX, latBand
         
      end if
      
      ttBodyIndex = -1
      esBodyIndex = -1
      
      BODY: do   ! Find the body indices for the TT and ES observations for this header (pressure level)

        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY 

        bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex )

        ! Assumes that normal (non high-precision) codes are used for TT and ES observations in obsSpaceData
        if ( bufrCode == BUFR_NETT ) then
          ttBodyIndex = bodyIndex
        else if ( bufrCode == BUFR_NEES ) then
          esBodyIndex = bodyIndex
        end if
        
      end do BODY
      
      ! Reverse any old corrections and apply new corrections (if uaRevOnly=.false.)
      
      ! TT bias correction
      if ( ttBodyIndex >= 0 ) then
        bodyIndex = ttBodyIndex
        ttOriginal  = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
        tt          = ttOriginal
        pressure    = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex ) * MPC_MBAR_PER_PA_R8
        flag        = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
        oldCorr     = obs_bodyElem_r(obsSpaceData, OBS_BCOR, bodyIndex )
        corr        = MPC_missingValue_R8
        if ( debug .and. pressure == 500.0d0 ) write(*,*) 'TTin',ttBodyIndex,pressure,tt
        if ( tt /= MPC_missingValue_R8 ) then
          if ( btest(flag, 6) .and. oldCorr /= MPC_missingValue_R8 ) then
            tt = tt - oldCorr
            flag = ibclr(flag, 6)
          end if
          if ( uaRevOnly ) corr = 0.0d0
          if ( .not.uaRevOnly ) then
            ! Get the TT correction for this station, sonde type, time-of-day and level
            ! stnIndex   = -1 if station not found in bcor file
            ! sondeTypeIndex =  0 if sonde-type code is not associated with any of the types in
            !                namelist (rsTypes) including stypes "Others" and "None" [code=0]
            countTTObs = countTTObs + 1
            if ( realRS41 ) then
              corr = MPC_missingValue_R8
              countRS41 = countRS41 + 1
            else
              call bcc_GetUACorrection('TT',stnIndex,sondeTypeIndex,sondeType,groupIndex,timeOfDayX,latBand,pressure,corr,sourceCorr)
              if ( sourceCorr == 'stn' ) then
                 countTTCorrByStation = countTTCorrByStation + 1
              else if ( sourceCorr == 'stype' ) then
                 countTTCorrByStype = countTTCorrByStype + 1
              end if
              if ( debug .and. pressure == 500.0d0 ) write(*,*) 'corrTT, source, groupIndex = ',corr, sourceCorr, groupIndex
            end if
            if ( corr /= MPC_missingValue_R8 ) then
              tt = tt + corr
              flag = ibset(flag, 6)
              countTTCorrections = countTTCorrections + 1
            else
              if ( .not.realRS41 .and. uaRejUnBcor ) flag = ibset(flag, 11)
            end if
          end if
        end if
        if ( debug .and. pressure == 500.0d0 ) write(*,*) 'TTout, corr, flag = ', tt, corr, flag
        call obs_bodySet_r( obsSpaceData, OBS_BCOR, bodyIndex, corr )
        call obs_bodySet_r( obsSpaceData, OBS_VAR , bodyIndex, tt   )
        call obs_bodySet_i( obsSpaceData, OBS_FLG , bodyIndex, flag ) 
      end if
      
      ! ES bias correction
      if ( esBodyIndex >= 0 .and. ttBodyIndex >= 0 ) then
        bodyIndex   = esBodyIndex
        esOriginal  = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
        es          = esOriginal
        pressure    = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex ) * MPC_MBAR_PER_PA_R8
        flag        = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
        oldCorr     = obs_bodyElem_r(obsSpaceData, OBS_BCOR, bodyIndex )
        corr        = MPC_missingValue_R8
        if ( debug .and. pressure == 500.0d0 ) write(*,*) 'ESin',esBodyIndex,pressure,es
        if ( es /= MPC_missingValue_R8 ) then
          if ( btest(flag, 6) .and. oldCorr /= MPC_missingValue_R8 ) then
            es = es - oldCorr
            flag = ibclr(flag, 6)
          end if
          if ( uaRevOnly ) corr = 0.0d0
          if ( .not.uaRevOnly ) then
            if ( realRS41 ) then
              corr = MPC_missingValue_R8
            else
              call bcc_GetUACorrection('TD',stnIndex,sondeTypeIndex,sondeType,groupIndex,timeOfDayX,latBand,pressure,corr,sourceCorr)
              if ( debug .and. pressure == 500.0d0 ) write(*,*) 'corrTD, source, groupIndex = ',corr, sourceCorr, groupIndex
            end if
            if ( corr /= MPC_missingValue_R8 ) then
              td = (ttOriginal - es) + corr
              es = tt - td
              if ( es < 0.0 ) es =  0.0d0
              if ( es > 30.0) es = 30.0d0
              corr =  es - esOriginal
              if ( debug .and. pressure == 500.0d0 ) write(*,*) 'ES corr =',corr
              flag = ibset(flag, 6)
              countESCorrections = countESCorrections + 1
            else
              if ( .not.realRS41 .and. uaRejUnBcor ) flag = ibset(flag, 11)
            end if
          end if
        end if
        if ( debug .and. pressure == 500.0d0 ) write(*,*) 'ESout, corr, flag = ', es, corr, flag
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
        p4 = (float(countRS41)/float(countTTObs))*100.
        write (*, '(a60, f7.2)' ) "bcc_applyUABcor: Percentage corrected by STATION     = ", p1
        write (*, '(a60, f7.2)' ) "bcc_applyUABcor: Percentage corrected by SONDE-TYPE  = ", p2
        write (*, '(a60, f7.2)' ) "bcc_applyUABcor: Percentage uncorrected              = ", p3
        write (*, '(a60, f7.2)' ) "bcc_applyUABcor: Percentage RS41                     = ", p4
        write (*, '(a60, i10)' )  "bcc_applyUABcor: Number of night profiles            = ", countNight
        write (*, '(a60, i10)' )  "bcc_applyUABcor: Number of day profiles              = ", countDay
        write (*, '(a60, i10)' )  "bcc_applyUABcor: Number of dawn-dusk profiles        = ", countDawnDusk
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
    if ( allocated(biasCorrPresentStype) ) deallocate(biasCorrPresentStype)
    if ( allocated(biasCorrPresentStn) )   deallocate(biasCorrPresentStn)
    if ( allocated(uaStations) )           deallocate(uaStations)
    
    if ( allocated(rsTypes) ) deallocate(rsTypes)
    
    write(*,*) "bcc_applyUABcor: end"
    
  end subroutine bcc_applyUABcor
  
  
  !-----------------------------------------------------------------------
  ! bcc_biasActive
  !-----------------------------------------------------------------------
  function bcc_biasActive(obsFam) result(biasActive)
    !
    ! :Purpose: returns True if bias correction is active for the given conventional observation family
    !
    implicit none
 
    character(len=*),intent(in) :: obsFam
    
    logical :: biasActive 

    if (.not.initialized) call bcc_readConfig()
    
    select case(trim(obsFam))
    case('GP')
      biasActive = bcc_gpBiasActive
    case('UA')
      biasActive = bcc_uaBiasActive
    case('AI')
      biasActive = bcc_aiBiasActive
    case default
      biasActive = .false.
    end select

  end function bcc_biasActive


end MODULE biasCorrectionConv_mod


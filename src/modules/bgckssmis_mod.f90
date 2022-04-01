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
!-------------------------------------- LICENCE end --------------------------------------

module bgckssmis_mod
  ! MODULE bgckssmis_mod (prefix='ssbg' category='1. High-level functionality')
  !
  ! :Purpose: Variables for ssmis background check and quality control.
  !
  use mpi_mod
  use MathPhysConstants_mod
  use utilities_mod
  use obsSpaceData_mod
  use tovs_nl_mod
  use obsErrors_mod

  implicit none
  save
  private

  ! Public functions/subroutines
  public :: ssbg_computeSsmisSurfaceType
  public :: ssbg_bgCheckSSMIS
  real    :: ssbg_clwQcThreshold
  logical :: ssbg_debug

  real,    parameter :: ssbg_realMissing=-99. 
  integer, parameter :: ssbg_intMissing=-1
  ! Other variables:
  real,    parameter :: ssbg_zmisg=9.9e09
  real,    parameter :: ssbg_rmisg=-999.0
  real,    parameter :: ssbg_clwThresh=0.02
  integer, parameter :: ssbg_mxval=30
  integer, parameter :: ssbg_maxObsNum=2500
  real,    parameter :: clw_amsu_rej=0.3
  real,    parameter :: clw_amsu_rej_ch3=0.1
  !  Highest peaking AMSU-A like SSMIS channel for ocean-only and CLW filtering
  !    3 = mid-troposphere   (AMSU/operations -- AMSU chan. 5)
  !    4 = upper-troposphere (scat. index used in AMSU/operations -- AMSU chan. 6)
  !       (AMSU-A scat. index cannot be computed here; need AMSU-A channels 1,2)
  integer, parameter :: ipc=4
  ! Module variable

  character(len=128), parameter :: fileGlmg='fstglmg'  ! glace de mer file
  character(len=128), parameter :: fileGlace='bicefil'  ! binaire 0.1degre ice file
  character(len=128), parameter :: fileWentz='wentz_surf.std'  ! surface wentz file
  character(len=128), parameter :: algOption = 'fwentz'
  ! Other NRL thresholds

  integer, parameter :: ssbg_maxNumSat  = 4
  integer, parameter :: ssbg_maxNumChan = 24
  integer, parameter :: ssbg_maxNumTest = 16

  ! namelist variables
  logical                       :: RESETQC                       ! reset Qc flags option
  logical                       :: debug                         ! debug mode


  namelist /nambgck/debug, RESETQC

contains

  subroutine ssbg_init()
    ! :Purpose: This subroutine reads the namelist section NAMBGCK
    !           for the module.

    implicit none

    ! Locals
    integer           :: ierr
    integer           :: nulnam

    ! External functions
    integer, external :: fclos
    integer, external :: fnom

    ! Default values for namelist variables
    debug = .false.
    RESETQC = .false.

    nulnam = 0
    ierr = fnom(nulnam, './flnml','FTN+SEQ+R/O', 0)
    read(nulnam, nml=nambgck, iostat=ierr)
    if (ierr /= 0) call utl_abort('ssbg_init: Error reading namelist')
    if (mpi_myid == 0) write(*, nml=nambgck)
    ierr = fclos(nulnam)

    ssbg_debug = debug

  end subroutine ssbg_init 

  !--------------------------------------------------------------------
  ! ssmis_tb2ta
  !--------------------------------------------------------------------
  subroutine ssmis_tb2ta(numObsToProcess, grossRej, ztb, zta)
    ! :Purpose: Convert Tbs received from UKMO to Tas, by reversing Ta to Tb
    !           spillover correction applied in B. Bell's pre-processing.

    implicit none

    ! Arguments
    integer, intent(in)  :: numObsToProcess  ! Number of obs points to process
    logical, intent(in)  :: grossRej(:)      ! Gross rejection indicator
    real,    intent(in)  :: ztb(:)           ! Tbs from input BURP file
    real,    intent(out) :: zta(:)           ! Tas after conversion

    ! Locals
    integer :: hiIndex
    integer :: loIndex
    integer :: obsIndex

    real    :: spillCoeffs(ssbg_maxNumChan)  ! Spillover correction coefficients

    !  Define spillover correction coefficients
    !!  Row1           ch1  ch2  ch3  ch4   ch5   ch6
    !!  Row2           ch7  ch8  ch9  ch10  ch11  ch12
    !!  Row3           ch13 ch14 ch15 ch16  ch17  ch18
    !!  Row4           ch19 ch20 ch21 ch22  ch23  ch24

    !  Spillover coeff for all channels (from Steve Swadley/NRL 9 June 2010)

    !  spillCoeffs = (/ 0.9850,  0.9850,  0.9850,  0.9850,  0.9850,  0.9815, &
    !                    0.9815,  0.9949,  0.9934,  0.9934,  0.9934,  0.9680, &
    !                    0.9720,  0.9820,  0.9810,  0.9850,  0.9820,  0.9780, &
    !                    0.9815,  0.9815,  0.9815,  0.9815,  0.9815,  0.9815 /)
    !
    !  Spillover coeff for 7 SSM/I-like channels (ch. 12-18) only
    !
    spillCoeffs = (/ 1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000, &
                      1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  0.9680, &
                      0.9720,  0.9820,  0.9810,  0.9850,  0.9820,  0.9780, &
                      1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000 /)

    !  Apply Tb -> Ta conversion

    loIndex = 1
    do obsIndex = 1, numObsToProcess

      hiIndex = obsIndex*ssbg_maxNumChan
      if ( .not. grossRej(obsIndex) ) then
        zta(loIndex:hiIndex) = spillCoeffs(:) * ztb(loIndex:hiIndex)
      else
        zta(loIndex:hiIndex) = ztb(loIndex:hiIndex)
      end if

      loIndex = hiIndex + 1

    end do

  end subroutine ssmis_tb2ta

  !--------------------------------------------------------------------------
  ! f16tdr_remapping
  !--------------------------------------------------------------------------
  subroutine f16tdr_remapping(satId, SSMIS_Ta, Remapped_SSMI_Ta)
    ! :Purpose: Remap SSMIS imaging channel antenna temperature to SSMI Ta

    !       SSMIS         C_Freq          SSMI
    !      -------   -----------------   -------
    !       Chan12    19.35h              Chan2
    !       Chan13    19.35v              Chan1
    !       Chan14    22.235v             Chan3
    !       Chan15    37.0h               Chan5
    !       Chan16    37.0v               Chan4
    !       Chan17    91.65v -> 85.5v     Chan6
    !       Chan18    91.65h -> 85.5h     Chan7

    implicit none

    ! Arguments
    integer, intent(in)  :: satId                             ! Satellite ID
    real,    intent(in)  :: SSMIS_Ta(ssbg_maxNumChan)         ! SSMIS antenna temperature
    real,    intent(out) :: Remapped_SSMI_Ta(ssbg_maxNumChan) ! Remapped SSMI antenna temperature

    ! Locals
    integer, parameter :: f16_id = 1
    integer, parameter :: f17_id = 2
    integer, parameter :: f18_id = 3

    integer(2)         :: channelIndex

    real(8)            :: CP(ssbg_maxNumChan)
    real               :: tbx(ssbg_maxNumChan)

    ! NESDIS Intercept/Slope for F16 SSMIS 7 IMG channels to SSMI linear remapping
                               ! ch1  ch2  ch3  ch4  ch5  ch6  ch7  ch8  ch9  ch10 ch11 (LAS/ENV)
    real(8), parameter :: AP(ssbg_maxNumChan)=(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                   7.44254,7.80472,6.76383,8.55426,7.34409,6.57813,6.45397, &
                                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

    real(8), parameter :: BP(ssbg_maxNumChan)=(/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,&
                                   0.969424,0.967519,0.959808,0.954316,0.958955,0.980339,0.978795, &
                                   1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)

    ! NESDIS F17 and F18 Tb biases with respect to F16
    real(8), parameter :: CP_F17(ssbg_maxNumChan)=(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                      -0.779,  -1.446,  -1.013,  -0.522,  -0.240,   0.735,   0.521,   &
                                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

    real(8), parameter :: CP_F18(ssbg_maxNumChan)=(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                      -0.773,  -0.688,  -1.031,  -0.632,  -0.411,   0.171,   0.928,   &
                                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)


    ! Initialization
    tbx(1:ssbg_maxNumChan) = SSMIS_Ta(1:ssbg_maxNumChan)

    if ( satId == f17_id ) then
      CP = CP_F17
    else if ( satId == f18_id ) then
      CP = CP_F18
    else
      CP = 0.0
    end if

    do channelIndex=1, ssbg_maxNumChan
      Remapped_SSMI_Ta(channelIndex) = AP(channelIndex) + BP(channelIndex)*(tbx(channelIndex)+CP(channelIndex))
    end do

  end subroutine f16tdr_remapping


  !--------------------------------------------------------------------------
  ! ssmi_ta2tb_fweng
  !--------------------------------------------------------------------------  
  subroutine ssmi_ta2tb_fweng(Ta, Tb)
    ! :Purpose: To convert antenna temperature(Ta) to brightness temperature(Tb).

    !  (1) All channel antenna gain spill-over correction
    !  (2) Imaging channel Cross-polarization correction
    !  (3) Doppler correction (to be developed)

    implicit none

    ! Arugments
    real, intent(in)  :: Ta(24) ! Antenna temperature
    real, intent(out) :: Tb(24) ! Brightness temperature

    ! Locals
    real(8), parameter :: AP(24)=(/0.9850,0.9850,0.9850,0.9850,0.9850,0.9790,0.9815,&
                                   0.9949,0.9934,0.9934,0.9934, &
                                   0.9690,0.9690,0.9740,0.9860,0.9860,0.9880,0.9880,&
                                   0.9815,0.9815,0.9815,0.9815,0.9815,0.9815/)
    real(8), parameter :: BP(24)=(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                   0.0, 0.0, 0.0, 0.0, &
                                   0.00415,0.00473,0.0107,0.02612,0.0217,0.01383,0.01947,&
                                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

    integer(4)         :: channelIndex
    real(4)            :: CP(24), DP(24)

    ! All channel antenna gain correction. Note the cross-polarization effects on antenna gain correction.
    do channelIndex=1, 24
      CP(channelIndex) = 1.0/( AP(channelIndex)*(1.0 - BP(channelIndex)) )
      DP(channelIndex) = CP(channelIndex) * BP(channelIndex)
      Tb(channelIndex) = CP(channelIndex)*Ta(channelIndex)
    end do

    ! SSMI IMG channel cross polarization correction
    Tb(12) = Tb(12) - DP(12)*Ta(13)               ! 19H
    Tb(13) = Tb(13) - DP(13)*Ta(12)               ! 19V
    Tb(14) = Tb(14) - DP(13)*(0.65*Ta(13)+96.6)   ! 22V
    Tb(15) = Tb(15) - DP(15)*Ta(16)               ! 37H
    Tb(16) = Tb(16) - DP(16)*Ta(15)               ! 37V
    Tb(17) = Tb(17) - DP(17)*Ta(18)               ! 85V
    Tb(18) = Tb(18) - DP(18)*Ta(17)               ! 85H

  end subroutine ssmi_ta2tb_fweng

  !--------------------------------------------------------------------------
  ! ssmi_ta2tb_fwengtz
  !--------------------------------------------------------------------------  
  subroutine ssmi_ta2tb_fwentz(Ta, Tb)
    ! :Purpose: Convert antenna temperatures to brightness temperatures.

    implicit none

    ! Arguments
    real, intent(in)  :: Ta(:) ! Antenna temperature
    real, intent(out) :: Tb(:) ! Brightness temperature

    ! Locals
    real :: AV19V,AH19V,A019V,AH19H,AV19H,A019H,AV22V,A022V
    real :: AV37V,AH37V,A037V,AH37H,AV37H,A037H
    real :: AV85V,AH85V,A085V,AH85H,AV85H,A085H
    real :: TA19V,TA19H,TA22V,TA37V,TA37H,TA85V,TA85H
    real :: TB19V,TB19H,TB22V,TB37V,TB37H,TB85V,TB85H

    !  Define APC coefficients

    AV19V =  1.0369830
    AH19V = -0.0039359
    A019V = -0.0892273

    AH19H =  1.0384912
    AV19H = -0.0054442
    A019H = -0.0892271

    AV22V =  1.01993
    A022V =  1.994

    AV37V =  1.0368094
    AH37V = -0.0222607
    A037V = -0.0392815

    AH37H =  1.0421693
    AV37H = -0.0276206
    A037H = -0.0392816

    AV85V =  1.0263188
    AH85V = -0.0143165
    A085V = -0.0324062

    AH85H =  1.0321901
    AV85H = -0.0201877
    A085H = -0.0324065

    !  Extract list of Tas from parent array.

    TA19H = Ta(12)
    TA19V = Ta(13)
    TA22V = Ta(14)
    TA37H = Ta(15)
    TA37V = Ta(16)
    TA85V = Ta(17)
    TA85H = Ta(18)

    !  Apply APC corrections: convert TDR's to SDR's.

    TB19V= (AV19V * TA19V) + (AH19V * TA19H) + (A019V)
    TB19H= (AH19H * TA19H) + (AV19H * TA19V) + (A019H)
    TB22V= (AV22V * TA22V) + A022V
    TB37V= (AV37V * TA37V) + (AH37V * TA37H) + (A037V)
    TB37H= (AH37H * TA37H) + (AV37H * TA37V) + (A037H)
    TB85V= (AV85V * TA85V) + (AH85V * TA85H) + (A085V)
    TB85H= (AH85H * TA85H) + (AV85H * TA85V) + (A085H)

    !  Insert list of new Tbs into parent array.

    Tb(12) = TB19H
    Tb(13) = TB19V
    Tb(14) = TB22V
    Tb(15) = TB37H
    Tb(16) = TB37V
    Tb(17) = TB85V
    Tb(18) = TB85H

  end subroutine ssmi_ta2tb_fwentz

  !--------------------------------------------------------------------------
  ! compute_iwv_101 
  !--------------------------------------------------------------------------  
  subroutine compute_iwv_101(Tb, iwv)
    ! :Purpose: Compute integrated water vapor from SSMI brightness temperatures.

    implicit none

    ! Arguments
    real, intent(in) :: Tb(24) ! Brightness temperature
    real, intent(out) :: iwv   ! Integrated water vapor (kg/m**2)

    ! Locals
    integer :: precipScreen    ! = 1: possible presence of precipitation does
                               !      not allow retrieval of IWV and CLW.
                               ! = 0: retrieval is possible

    real    :: IWV_alishouse   ! estimated total precipitable water, with cubic polynomial correction (=iwv)
    real    :: IWV_alishouse0  ! estimated total precipitable water, without cubic polynomial correction
    real    :: precipThresh
    real    :: tb19v,tb22v,tb37v,tb37h

    real    :: ciwv(0:4)

    !  Initializations.

    ciwv = (/232.89393,-0.148596,-1.829125,-0.36954,0.006193/)

    tb19v = Tb(13)
    tb22v = Tb(14)
    tb37h = Tb(15)
    tb37v = Tb(16)

    !  Compute IWV_alishouse.

    IWV_alishouse0 = ciwv(0) + ciwv(1)*tb19v + ciwv(2)*tb22v +   &
          &          ciwv(3)*tb37v + ciwv(4)*tb22v**2

    !  Apply Cubic Polynomial Correction to the CAL/VAL water vapour
    !  algorithm (G.W. Petty).

    IWV_alishouse = -3.75 + 1.507*IWV_alishouse0 - 0.01933*IWV_alishouse0**2 +  &
          &         (2.191e-04)*IWV_alishouse0**3

    !  Compute Precipitation Screen.

    precipThresh = -11.7939 - 0.02727*tb37v + 0.09920*tb37h
    if ( precipThresh > 0.0 ) then
      precipScreen = 1
    else
      precipScreen = 0
    end if

    !  Apply precipitation screen to IWV_alishouse (units are kg/m**2).
    !  Bounds now applied in cld_filter_fweng.ftn90 routinge.

    if( precipScreen == 1 ) then
       IWV_alishouse = ssbg_rmisg
    end if

    !  Store the IWV value.

    iwv = IWV_alishouse

  end subroutine compute_iwv_101

  !--------------------------------------------------------------------------
  ! determ_tpw
  !--------------------------------------------------------------------------  
  subroutine determ_tpw(Tb, sType, seaIce, TPW)
    ! :Purpose: To calculate total precipitable water (in mm).

    implicit none

    ! Arguments
    real,    intent(in)  :: Tb(24)   ! Brightness temperature
    integer, intent(in)  :: sType    ! Surface type
    real,    intent(in)  :: seaIce   ! Sea ice coverage
    real,    intent(out) :: TPW      ! Total precipitable water (mm)

    ! Locals
    integer, parameter :: ocean=0

    real :: SCT
    real :: Tb19V, Tb19H, Tb22V, Tb37V, Tb37H, Tb85V, Tb85H

    ! Extract 7 IMG channels from Tb
    Tb19v = Tb(13)
    Tb19h = Tb(12)
    Tb22v = Tb(14)
    Tb37v = Tb(16)
    Tb37h = Tb(15)
    Tb85v = Tb(17)
    Tb85h = Tb(18)

    TPW = ssbg_rmisg

    ! No TWP over land and sea ice
    if ( sType == ocean ) then

      SCT = -182.7 + 0.75*Tb19V + 2.543*Tb22V - 0.00543*Tb22V*Tb22V - Tb85V

      if ( abs(seaIce) < 70.0 ) then
        TPW = 232.89393 - 0.148596*Tb19V - 1.829125*Tb22V + 0.006193*Tb22V**2 - 0.36954*Tb37V
        TPW = -3.753 + 1.507*TPW - 0.01933*TPW**2 + 0.0002191*TPW**3

        if ( TPW < 0.0 ) TPW = 0.0
        if ( TPW > 80.0 ) TPW = 80.0
      end if

    end if

  end subroutine determ_tpw

  !--------------------------------------------------------------------------
  !  determ_sea_ice
  !--------------------------------------------------------------------------  
  subroutine determ_sea_ice(ocean, Ta, sType, seaIce, latitude)
    ! :Purpose: To calculate sea ice cover (in %).

    implicit none

    ! Arguments
    integer, intent(in)  :: ocean    ! Ocean surface type index
    real,    intent(in)  :: Ta(24)   ! Antenna temperature
    integer, intent(in)  :: sType    ! Surface type
    real,    intent(out) :: seaIce   ! Sea ice coverage
    real,    intent(in)  :: latitude ! Latitude of observation

    ! Locals
    real :: Ta19v, Ta19h, Ta22v, Ta37v, Ta37h, Ta85v, Ta85h

    ! Extract 7 IMG channels from Ta
    Ta19v = Ta(13)
    Ta19h = Ta(12)
    Ta22v = Ta(14)
    Ta37v = Ta(16)
    Ta37h = Ta(15)
    Ta85v = Ta(17)
    Ta85h = Ta(18)

    ! Calculate Sea Ice Coverage
    if ( sType == ocean ) then   ! Over Ocean
      if (latitude > 44.4 .or. latitude < -52.0 ) then
        seaIce = 91.9 - 2.994*Ta22v + 2.846*Ta19v - 0.386*Ta37v + 0.495*Ta85v &
                  + 1.005*Ta19h - 0.904*Ta37h
        if ( seaIce >= 70.0 ) then
          seaIce = 100.0
        else
          seaIce = 0.0
        end if
      else
        seaIce = 0.0
      end if
    else         ! Over Land, there is no sea ice!
      seaIce = ssbg_rmisg
    end if

  end subroutine determ_sea_ice

  !--------------------------------------------------------------------------
  !  determ_clw
  !--------------------------------------------------------------------------  
  subroutine determ_clw(algOption, Ta, Tb, sType, CLW, IWV, latitude)
    ! :Purpose: To calculate cloud liquid water for a single data point (in kg/m**2).

    !    Normally called when sType = 0 (open water).
    !    Also retrieves sea-ice to see if there is sea-ice at water point.

    implicit none

    ! Arguments
    character(len=6), intent(in)    :: algOption ! Algorithm option (fweng, fwentz or nsun)
    real,             intent(in)    :: Ta(24)    ! Antenna temperature
    real,             intent(in)    :: Tb(24)    ! Brightness temperature
    integer,          intent(inout) :: sType     ! Surface type
    real,             intent(out)   :: CLW       ! Cloud liquid water (in kg/m**2)
    real,             intent(inout) :: IWV       ! Integrated water vapor
    real,             intent(in)    :: latitude  ! Latitude of observation

    ! Locals
    integer,    parameter :: isSeaIce = 1
    integer,    parameter :: ocean  = 0
    integer(4), parameter :: RT = 285
    real,       parameter :: clwLimit = 6.0

    real :: ALG1, ALG2, ALG3
    real :: seaIce
    real :: Ta19v, Ta19h, Ta22v, Ta37v, Ta37h, Ta85v, Ta85h
    real :: Tb19v, Tb22v, Tb37v
    real :: TPW

    CLW = ssbg_rmisg

    ! Extract 7 IMG channels from Ta and Tb
    Ta19v = Ta(13)
    Ta19h = Ta(12)
    Ta22v = Ta(14)
    Ta37v = Ta(16)
    Ta37h = Ta(15)
    Ta85v = Ta(17)
    Ta85h = Ta(18)

    ! Call determ_sea_ice to find the seaIce
    !  -- seaIce = 100.0  when sea ice >= 70%
    !             = 0.0    when sea ice  < 70%
    call determ_sea_ice(ocean, Ta, sType, seaIce, latitude)

    ! Calculate CLW Over Ocean
    if ( (sType == ocean) .and. (seaIce /= 100.0) ) then
      ALG1 = ssbg_rmisg
      ALG2 = ssbg_rmisg
      ALG3 = ssbg_rmisg

      if ( trim(algOption) == 'fweng' ) then
        ! Compute IWV using F. Weng algorithm.
        IWV = 232.89 - 0.1486*Tb19v - 0.3695*Tb37v - (1.8291 - 0.006193*Tb22v)*Tb22v
        if ( IWV < 0.0 ) IWV = 0.0
      end if
    if ( trim(algOption) == 'nsun') then
      call determ_tpw(Tb,sType,seaIce,TPW)
      IWV = TPW
    end if

    if ( (Ta19v < RT) .and. (Ta22v < RT) ) then
      ALG1 = -3.20 * ( ALOG(290.0-Ta19v) - 2.80 - 0.42*ALOG(290.0-Ta22v) )          !TA
      ! ALG1 = -3.20 * ( ALOG(290.0-Tb19v) - 2.84 - 0.40*ALOG(290.0-Tb22v) )      !TB
    end if

    if ( (Ta37v < RT) .and. (Ta22v < RT) ) then
      ALG2 = -1.66 * ( ALOG(290.0-Ta37v) - 2.90 - 0.349*ALOG(290.0-Ta22v) )   !TA
      ! ALG2 = -1.66 * ( ALOG(290.0-Tb37v) - 2.99 - 0.32*ALOG(290.0-Tb22v) )    !TB
    end if

    if ( (Ta85h < RT) .and. (Ta22v < RT) ) then
      ALG3 = -0.44 * ( ALOG(290.0-Ta85h) + 1.60 - 1.354*ALOG(290.0-Ta22v) )     !TA
      ! ALG3 = -0.44 * ( ALOG(290.0-Tb85h) + 1.11 - 1.26*ALOG(290.0-Tb22v) )      !TB
    end if

    if ( ALG1 > 0.70 ) then
      CLW = ALG1
    else if ( ALG2 > 0.28 ) then
      CLW = ALG2
    else if ( IWV < 30.0 ) then
      CLW = ALG3
    else
      CLW = ALG2
    end if

    ! Verify CLW is within acceptable upper limit.
    if ( CLW > clwLimit ) CLW = ssbg_rmisg

    ! Force negative CLW values to zero.
    if ( CLW < 0.0 .and. CLW /= ssbg_rmisg ) CLW = 0.0

    else
      ! Sea Ice (>70%) detected from s/r determ_sea_ice but sType was 0 = waterobs (on call)
      CLW = -500.0
      IWV = ssbg_rmisg
      sType = isSeaIce

    end if

  end subroutine determ_clw

  !--------------------------------------------------------------------------
  !  cld_filter_fweng
  !--------------------------------------------------------------------------  
  subroutine cld_filter_fweng(numObsToProcess, obsTb, algOption, waterObs, grossRej,  &
            &                 cloudObs, iwvReject, precipObs, rclw, riwv, iSatId,   &
            &                 obsLatitude, numSeaIceObs)
    ! :Purpose: Compute the cloud liquid water (CLW) from SSMIS channels using the
    !           regression algorithm of Fuzhong Weng and Ninghai Sun.
    !           Retrieve CLW path from F16 SSMIS TDR data

    implicit none

    ! Arguments
    integer,          intent(in)    :: numObsToProcess  ! Number of obs points to process
    real,             intent(in)    :: obsTb(:)         ! Brightness temperature of observations
    character(len=6), intent(in)    :: algOption        ! Algorithm option (fweng, fwentz or nsun)
    logical,          intent(inout) :: waterObs(:)      ! Open water identifier for each obs
    logical,          intent(in)    :: grossRej(:)      ! Logical array of obs with gross error (obs to reject)
    logical,          intent(inout) :: cloudObs(:)      ! Logical array of obs for which CLW > 0.01 kg/m**2 or with precipitations
    logical,          intent(inout) :: iwvReject(:)     ! Logical array of obs for which IWV > 80 kg/m**2
    logical,          intent(inout) :: precipObs(:)     ! Logical array of obs with precipitations (CLW missing)
    real,             intent(inout) :: rclw(:)          ! Real array of CLW
    real,             intent(inout) :: riwv(:)          ! Real array of integrated water vapor (IWV)
    integer,          intent(in)    :: iSatId           ! Satellite identifier
    real,             intent(in)    :: obsLatitude(:)   ! Observation latitudes
    integer,          intent(inout) :: numSeaIceObs     ! Number of observations with sea ice

    ! Locals
    real,    parameter :: iwvThresh = 80.0    ! Upper bound for IWV in kg/m**2
    integer, parameter :: ocean  = 0
    integer, parameter :: seaIce = 1

    integer            :: hiIndex
    integer            :: loIndex
    integer            :: obsIndex
    integer            :: sType

    real               :: clw
    real               :: iwv
    real               :: latitude

    real               :: F16TDR(ssbg_maxNumChan)
    real               :: remappedTa(ssbg_maxNumChan)
    real               :: Tb(ssbg_maxNumChan)
    real               :: zta(ssbg_mxval*ssbg_maxObsNum)

    !---------------------------------------------------

    ! Convert Tbs received from UKMO to Tas, by reversing Ta to Tb
    ! spillover correction applied in B. Bell's pre-processing.
    ! Missing Tbs are also checked for here.
    ! Apply to all obs of current record.

    call ssmis_tb2ta(numObsToProcess,grossRej,obsTb,zta)

    rclw(:) = ssbg_rmisg
    riwv(:) = ssbg_rmisg

    loIndex = 1
    HEADER: do obsIndex = 1, numObsToProcess

      hiIndex = obsIndex*ssbg_maxNumChan
      F16TDR(:) = zta(loIndex:hiIndex)
      latitude = obsLatitude(obsIndex)
      loIndex = hiIndex + 1

      ! Obtain CLW and IWV for obs points over open ocean, and where
      ! Tb values have passed gross filter check.
      ! Initialize IWV and CLW to missing for all other cases.

      clw = ssbg_rmisg
      iwv = ssbg_rmisg

      if ( waterObs(obsIndex) .and. (.not. grossRej(obsIndex)) ) then

        sType = ocean     ! Surface-type=Ocean

        ! Call SSMIS TDR to SSMI TDR remapping subroutine

        call f16tdr_remapping(iSatId, F16TDR, remappedTa)

        ! Call SSM/I Ta to Tb conversion subroutine

        if ( trim(algOption) == 'fweng' ) then
          call ssmi_ta2tb_fweng(remappedTa, Tb)
          ! IWV computed in determ_clw subroutine below.
        else if ( trim(algOption) == 'fwentz' ) then
          call ssmi_ta2tb_fwentz(remappedTa, Tb)
          ! Compute IWV using Alishouse and Petty,
          ! because it won't be computed in determ_clw.
          ! Missing value for IWV means possible precipitation
          ! and so CLW will not be computed
          call compute_iwv_101(Tb,iwv)
          if ( iwv == ssbg_rmisg ) precipObs(obsIndex) = .true.
        else if ( trim(algOption) == 'nsun' ) then
          call ssmi_ta2tb_fweng(remappedTa, Tb)
          ! IWV computed in determ_clw subroutine below.
        else
          write(*,*) ' cld_filter_fweng: Invalid algorithm option !! '
          call abort()
        end if

        ! Call CLW retrieval algorithm subroutine.
        !  -- also computes and returns (output) IWV if algOption=fweng or nsun
        !  -- sType is also changed to seaIce value if sea-ice is detected

        if ( trim(algOption) /= 'fwentz' ) then
          call determ_clw(algOption, remappedTa, Tb, sType, clw, iwv, latitude)
        else
          if ( .not. precipObs(obsIndex) ) then
            call determ_clw(algOption, remappedTa, Tb, sType, clw, iwv, latitude)
          end if
        end if

        ! Check for newly detected sea-ice (sType changed)

        if ( sType == seaIce ) numSeaIceObs = numSeaIceObs + 1

        ! Reject obs with precipitation or cloud amount more than threshold.
        ! Also, reject obs if CLW is missing (undetermined)
        ! Set waterObs flag to false if deter_clw returns "seaIce" value (-500)
        if ( (clw > ssbg_clwThresh) .or. precipObs(obsIndex) ) cloudObs(obsIndex) = .true.
        if ( clw == ssbg_rmisg )  cloudObs(obsIndex) = .true.
        if ( clw == -500.0 ) then
           waterObs(obsIndex) = .false.
           clw = ssbg_rmisg
        end if
        ! Reject obs with IWV value more than threshold and set IWV to missing.
        ! First, set IWV values below zero to 0.0.
        if ( iwv < 0.0 .and. iwv /= ssbg_rmisg ) iwv = 0.0
        if ( iwv > iwvThresh ) then
          if ( .not. cloudObs(obsIndex) ) iwvReject(obsIndex) = .true.
          iwv = ssbg_rmisg
        end if

        ! Store CLW and IWV
        if (ssbg_debug) then
          write(*,*)'cld_filter_fweng: CLOUD BY DETERM_CLW = ', clw
          write(*,*)'cld_filter_fweng: IWV BY DETERM_CLW = ', iwv
        end if
        rclw(obsIndex) = clw
        riwv(obsIndex) = iwv

      end if

    end do HEADER

    120  format(' obsIndex ',2x,i3,2x,' clw ',f4.2,2x,' CLW threshold = ',f4.2,2x,' cloudObs ',l2)
    130  format(' obsIndex ',2x,i3,2x,' iwv ',f4.2,2x,' IWV threshold = ',f4.2,2x,' iwvReject ',l2)

  end subroutine cld_filter_fweng

  !--------------------------------------------------------------------------
  ! copy1Dimto2DimRealArray
  !--------------------------------------------------------------------------
  subroutine copy1Dimto2DimRealArray(oneDimArray, firstDim, secondDim, twoDimArray)
    ! :Purpose: copy 1D real array into 2D real array given firstDim and secondDim

    implicit none

    ! Arguments
    integer, intent(in)    :: firstDim                         ! First dimension
    integer, intent(in)    :: secondDim                        ! Second dimension
    real,    intent(in)    :: oneDimArray(firstDim*secondDim)  ! 1D real array
    real,    intent(inout) :: twoDimArray(firstDim,secondDim)  ! 2D real array
    
    ! Locals
    integer                :: firstDimIndex
    integer                :: productDimIndex
    integer                :: secondDimIndex

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    do secondDimIndex=1,secondDim
      do firstDimIndex=1,firstDim
        productDimIndex = (secondDimIndex-1)*firstDim + firstDimIndex 
        twoDimArray(firstDimIndex,secondDimIndex) = oneDimArray(productDimIndex)
      end do
    end do

  end subroutine copy1Dimto2DimRealArray

  !--------------------------------------------------------------------------
  ! copy1Dimto2DimIntegerArray
  !--------------------------------------------------------------------------
  subroutine copy1Dimto2DimIntegerArray(oneDimArray, firstDim, secondDim, twoDimArray)
    ! :Purpose: copy 1D integer array into 2D Integer array given firstDim and secondDim

    implicit none

    ! Arguments
    integer, intent(in)    :: firstDim                         ! First dimension
    integer, intent(in)    :: secondDim                        ! Second dimension
    integer, intent(in)    :: oneDimArray(firstDim*secondDim)  ! 1D integer array
    integer, intent(inout) :: twoDimArray(firstDim,secondDim)  ! 2D integer array
    
    ! Locals
    integer :: firstDimIndex
    integer :: productDimIndex 
    integer :: secondDimIndex 

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    do secondDimIndex=1,secondDim
      do firstDimIndex=1,firstDim
        productDimIndex = (secondDimIndex-1)*firstDim + firstDimIndex 
        twoDimArray(firstDimIndex,secondDimIndex) = oneDimArray(productDimIndex)
      end do
    end do

  end subroutine copy1Dimto2DimIntegerArray

  !------------------------------------------------------------------------------------
  ! bennartz
  !------------------------------------------------------------------------------------
  subroutine bennartz (ier, numObsToProcess, tb89, tb150, satZenithAngle, landSeaQualifier, scatL, scatW)
    ! :Purpose: Compute the following parameters using 2 AMSU-B channels:
    !           - scattering index (over land and ocean).
    !           The two channels used are: 89Ghz, 150Ghz.

    !           REFERENCES     Bennartz, R., A. Thoss, A. Dybbroe and D. B. Michelson, 
    !                          1999: Precipitation Analysis from AMSU, Nowcasting SAF, 
    !                          Swedish Meteorologicali and Hydrological Institute, 
    !                          Visiting Scientist Report, November 1999.

    implicit none

    ! Arguments
    integer, intent(in)  :: numObsToProcess      ! Number of obs points to process
    integer, intent(out) :: ier(numObsToProcess) ! Error return code
    real,    intent(in)  :: tb89(:)              ! 89Ghz AMSU-B brightness temperature (K)
    real,    intent(in)  :: tb150(:)             ! 150Ghz AMSU-B brightness temperature (K)
    real,    intent(in)  :: satZenithAngle(:)    ! Satellite zenith angle (deg.)
    integer, intent(in)  :: landSeaQualifier(:)  ! Land/sea indicator (0=land; 1=ocean)
    real,    intent(out) :: scatL(:)             ! Scattering index over land
    real,    intent(out) :: scatW(:)             ! Scattering index over water

    ! Locals
    ! Notes: In the case where an output parameter cannot be calculated, the
    ! value of this parameter is the missing value, i.e. -99.
    real, parameter :: zMisg = -99.

    integer :: obsIndex

    ! ____1) Initialise output parameters:**    -------------------------------*

    do obsIndex = 1, numObsToProcess
      scatL(obsIndex) = zMisg
      scatW(obsIndex) = zMisg
    end do

    !____2) Validate input parameters:**    -----------------------------*
    do obsIndex = 1, numObsToProcess
      if ( tb89(obsIndex)               < 120.  .or.     &
            tb89(obsIndex)              > 350.  .or.     &
            tb150(obsIndex)             < 120.  .or.     &
            tb150(obsIndex)             > 350.  .or.     &
            satZenithAngle(obsIndex)    < -90.  .or.     &
            satZenithAngle(obsIndex)    >  90.  .or.     &
            landSeaQualifier(obsIndex)  <   0   .or.     &
            landSeaQualifier(obsIndex)  >   1        ) then
          ier(obsIndex) = 1
      else
          ier(obsIndex) = 0
      end if 
    end do

    !____3) Compute parameters:**    ----------------------*
    do obsIndex = 1, numObsToProcess
      if ( ier(obsIndex) == 0 ) then
        if (landSeaQualifier(obsIndex) == 1 ) then
          scatW(obsIndex) = (tb89(obsIndex)-tb150(obsIndex)) -    &
                     (-39.2010+0.1104*satZenithAngle(obsIndex))
        else
          scatL(obsIndex) = (tb89(obsIndex)-tb150(obsIndex)) -    &
                     (0.158+0.0163*satZenithAngle(obsIndex))
        end if
      else if ( (ier(obsIndex) /= 0  ) .and. (obsIndex <= 100 ) .and. (ssbg_debug)) then
        write(*,*), 'bennartz: Input Parameters are not all valid: '
        write(*,*), ' obsIndex,tb89(obsIndex),tb150(obsIndex),satZenithAngle(obsIndex),landSeaQualifier(obsIndex) = ',     &
                   obsIndex,tb89(obsIndex),tb150(obsIndex),satZenithAngle(obsIndex),landSeaQualifier(obsIndex)
        write(*,*), ' ier(obsIndex),scatL(obsIndex),scatW(obsIndex)=',     &
                   ier(obsIndex),scatL(obsIndex),scatW(obsIndex)
      end if
    end do

  end subroutine bennartz

  !--------------------------------------------------------------------------
  ! ssbg_readGeophysicFieldsAndInterpolate
  !--------------------------------------------------------------------------
  subroutine ssbg_readGeophysicFieldsAndInterpolate(obsLatitude, obsLongitude, modelInterpTer)
    ! :Purpose: Reads geophysical model variable (GZ) and saves for the first time.
    !           GZ is geopotential height (GZ at surface = surface height in dam).
    !           Then interpolates those variables to observation location.

    implicit none

    ! Arguments
    real,              intent(in)  :: obsLatitude(:)    ! Observation latitudes
    real,              intent(in)  :: obsLongitude(:)   ! Observation longitudes
    real, allocatable, intent(out) :: modelInterpTer(:) ! Filtered and interpolated topography (in m)

    ! Locals
    integer, parameter      :: mxLat = 5
    integer, parameter      :: mxLon = 5
    real,    parameter      :: dLat = 0.4
    real,    parameter      :: dLon = 0.6

    character(len=12)       :: etikxx
    character(len=1)        :: grtyp
    character(len=4)        :: nomvxx
    character(len=2)        :: typxx

    integer                 :: boxPointIndex
    integer                 :: boxPointNum
    integer                 :: dataIndex
    integer                 :: dataNum
    integer                 :: ezQkDef
    integer                 :: ezSetOpt
    integer                 :: gdllsval
    integer,           save :: gdgz                   ! topo interpolation param
    integer                 :: idum1, idum2, idum3, idum4
    integer                 :: idum5, idum6, idum7, idum8
    integer                 :: idum9, idum10, idum11, idum12
    integer                 :: idum13, idum14, idum15
    integer                 :: idum16, idum17, idum18
    integer                 :: ier
    integer                 :: ig1, ig2, ig3, ig4
    integer                 :: irec
    integer                 :: iUnGeo
    integer                 :: latIndex
    integer                 :: lonIndex
    integer                 :: ni, nj, nk
    integer                 :: nLat
    integer                 :: nLon
    integer                 :: nObsLat
    integer                 :: nObsLon

    logical,           save :: ifFirstCall = .true.   ! If .True. we read GL, GZ and MG

    real, allocatable, save :: GZ(:)                  ! Modele Topographie (GZ)
    real, allocatable       :: GZIntBox(:,:)
    real, allocatable       :: obsLatBox (:,:)
    real, allocatable       :: obsLonBox (:,:)

    real                    :: topoFact               ! Facteur x topo pour avoir des unites en metre
    real                    :: xLat
    real                    :: xLon

    ! External functions
    integer, external       :: fclos
    integer, external       :: fnom
    integer, external       :: fstfrm
    integer, external       :: fstinf
    integer, external       :: fstlir
    integer, external       :: fstouv
    integer, external       :: fstprm
    integer, external       :: ip1_all

    ! STEP 1: CHECK if obsLatitude AND obsLongitude ARE SAME DIMENSION
    nObsLat = size(obsLatitude)
    nObsLon = size(obsLongitude)
    if (nObsLat /= nObsLon) then
      call utl_abort ('ssbg_readGeophysicFieldsAndInterpolate: OBSERVATION obsLatitude and obsLongitude should have SAME LENGTH')
    else 
      dataNum = nObsLat
    end if

    ! STEP 2: READ GZ from the FST FILE
    if(ifFirstCall) then
      iUnGeo = 0
      ier = fnom(iUnGeo,'trlm_01','STD+RND+R/O',0)
      ier = fstouv(iUnGeo,'RND')

      ! Using hybrid coordinates
      irec = fstinf(iUnGeo,ni,nj,nk,-1,' ',ip1_all(1.0,5),-1,-1,' ','GZ')
      if (irec < 0) then
        call utl_abort('ssbg_readGeophysicFieldsAndInterpolate: LA TOPOGRAPHIE EST INEXISTANTE')
      end if
      topoFact = 10.0  ! dam --> m

      if (allocated(GZ)) deallocate(GZ)
      allocate ( GZ(ni*nj), stat=ier)
      if ( ier /= 0 ) then
        call utl_abort('ssbg_readGeophysicFieldsAndInterpolate: Allocation of array GZ failed')
      end if
      ier = fstlir(GZ,iUnGeo,ni,nj,nk,-1,' ',ip1_all(1.0,5),-1,-1,' ','GZ')

      GZ(:) = GZ(:)*topoFact

      ier = fstprm ( irec, idum1, idum2, idum3, idum4, &
                     idum5, idum6, idum7, idum8, idum9, idum10,  &
                     idum11, typxx, nomvxx, etikxx, grtyp, ig1, &
                     ig2, ig3, ig4, idum12, idum13, idum14,  &
                     idum15, idum16, idum17, idum18 )
      write (*,*) ' GRILLE GZ : ',grtyp,ni,nj, &
                     ig1,ig2,ig3,ig4
      ier  = ezSetOpt('INTERP_DEGREE','LINEAR')
      gdgz = ezQkDef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iUnGeo)

      ier = fstfrm(iUnGeo)
      ier = fclos(iUnGeo)
      ifFirstCall = .false.
    end if

    ! STEP 3:  Interpolation de la glace et le champ terre/mer du modele aux pts TOVS.
    ! N.B.: on examine ces champs sur une boite centree sur chaque obs.
    boxPointNum = mxLat*mxLon
    if(allocated(obsLatBox)) deallocate(obsLatBox)
    allocate (obsLatBox(boxPointNum, dataNum) , stat=ier)
    if(allocated(obsLonBox)) deallocate(obsLonBox)
    allocate (obsLonBox(boxPointNum, dataNum) , stat=ier)
    if(allocated(GZIntBox)) deallocate(GZIntBox)
    allocate (GZIntBox(boxPointNum, dataNum) , stat=ier)

    nLat = (mxLat-1)/2
    nLon = (mxLon-1)/2
    do dataIndex = 1, dataNum
      boxPointIndex = 0
      do latIndex = -nLat, nLat
        xLat = obsLatitude(dataIndex) +latIndex*dLat
        xLat = max(-90.0,min(90.0,xLat))
        do lonIndex = -nLon, nLon
          boxPointIndex = boxPointIndex + 1
          xLon = obsLongitude(dataIndex) +lonIndex*dLon
          if ( xLon < -180. ) xLon = xLon + 360.
          if ( xLon >  180. ) xLon = xLon - 360.
          if ( xLon <    0. ) xLon = xLon + 360.
          obsLatBox(boxPointIndex,dataIndex) = xLat
          obsLonBox(boxPointIndex,dataIndex) = xLon
        end do
      end do
    end do

    ier = gdllsval(gdgz,GZIntBox,GZ,obsLatBox,obsLonBox,boxPointNum*dataNum)
    if (ier < 0) then
      call utl_abort ('ssbg_readGeophysicFieldsAndInterpolate: ERROR in the interpolation of GZ')
    end if

    if(allocated(modelInterpTer)) deallocate(modelInterpTer)
    allocate (modelInterpTer(dataNum) , stat=ier)

    modelInterpTer(:) = 0.0
    do dataIndex = 1, dataNum
      if (ssbg_debug) then
        write(*,*), 'ssbg_readGeophysicFieldsAndInterpolate: infos'
        write(*,*), '   '
        write(*,*), ' dataIndex = ', dataIndex
        write(*,*), '   '
        write(*,*), ' obsLatitude, obsLongitude = ', obsLatitude(dataIndex), obsLongitude(dataIndex)
        write(*,*), '   '
        write(*,*), ' obsLatBox = '
        write(*,*),  (obsLatBox(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        write(*,*), ' obsLonBox = '
        write(*,*),  (obsLonBox(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        write(*,*), ' GZIntBox = '
        write(*,*),  (GZIntBox(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
      end if
      do boxPointIndex=1, boxPointNum
        modelInterpTer(dataIndex) = max(modelInterpTer(dataIndex),GZIntBox(boxPointIndex,dataIndex))
      end do
      if (ssbg_debug) then
        write(*,*), 'ssbg_readGeophysicFieldsAndInterpolate: modelInterpTer = ', modelInterpTer(dataIndex)
      end if
    end do

  end subroutine ssbg_readGeophysicFieldsAndInterpolate

  !--------------------------------------------------------------------
  !  land_ice_mask_ssmis
  !--------------------------------------------------------------------
  subroutine land_ice_mask_ssmis(numObsToProcess, obsLatitude, obsLongitude, landSeaQualifier, &
       &                         terrainType, waterObs)
    ! :Purpose: Determine for each observation point the ice mask value from
    !           the binary file copied to the local work directory.

    !           This may be a user-specified file or it is copied from
    !              /users/dor/afsi/sio/datafiles/data2/ade.maskice10
    !           The ice mask is updated every day at 00 UTC. The binary file has
    !           a resolution of 0.1 deg. Observations with an ice mask value of 0
    !           (=ice; land or sea=-1) are removed. The ice mask value also
    !           determines the terrain-type qualifier (element 13039) for each obs
    !           pt, which is required when writing output to boxed format.

    !           Finally, this routine performs a land/ice consistency check using
    !           the MG and LG fields used by the model which produces the trial field.
    !           The purpose of this check is to remove obs that reside close to coasts
    !           or ice, and so whose TBs may be contaminated.
    !           For the GEM global model, any of the analysis files provides the MG
    !           field, while for the meso-global model a user-specified file is required
    !           to define MG. In either case, the GEM Global analysis file provides the
    !           LG field, and it is interpolated to the resolution of the model.
    !             - MG from climate file :     /data/gridpt/dbase/anal/glbeta2/DATE_000
    !             - LG from analysis file: ex. /data/gridpt/dbase/anal/glbpres2/DATE_000

    !           In the application of this check, a 5x5 mesh, with spacing defined by rLatKm and
    !           rLonKm, is positioned with its center over an obs pt (2 grid pts on either side
    !           of the obs pt; size of mesh is equal to 4*rLatKm x 4*rLonKm). The values of MG
    !           and LG are evaluated at the grid points of this mesh. The maximum value of each
    !           determines whether the obs pt is too close to ice or land to be retained.
    !           **NOTE: the threshold value for MG has a very strong effect on the distance
    !                   from land that is permitted for an obs to be retained


    !      Maximum FOV             x---x---x---x---x     ^
    !         = 75km x 75km        |   |   |   |   |     |
    !         for Meso-sphere CHs  x---x---x---x---x     |
    !         = 74km x 47km        |   |   |   |   |     |
    !         for 19 GHz           x---x---o---x---x     | = 4*rLatKm
    !                              |   |   |   |   |     | = 4*40 km
    !                           ^  x---x---x---x---x     | = 160 km = 80 km north & south
    !                   rLatKm  |  |   |   |   |   |     |
    !                           v  x---x---x---x---x     v
    !                                          <--->
    !                                         rLonKm
    !
    !                              <--------------->
    !                                 = 4*rLonKm
    !                                 = 4*40 km
    !                                 = 160 km = 80 km east & west


    !               MG value = 1.0  ==>  LAND       MG value = 0.0  ==>  OCEAN
    !               LG value = 1.0  ==>  ICE        LG value = 0.0  ==>  NO ICE

    !--------------------------------------------------------------------
    ! Variable Definitions
    ! --------------------
    ! fileGlmg         - input  -  name of file holding model MG and LG fields (external)
    ! numObsToProcess  - input  -  number of input obs pts in record
    ! obsLatitude      - input  -  array holding lat values for all obs pts in record
    ! obsLongitude     - input  -  array holding lon values for all obs pts in record
    ! landSeaQualifier - in/out -  array holding land/sea qualifier values for all obs
    !                              pts of record (0 = land, 1 = sea)
    ! terrainType      - output -  array holding terrain-type values for all obs pts
    !                              of current record
    ! waterObs         - output -  logical array identifying for each obs in current record
    !                              whether it is over open water, far from coast/ice
    ! iMask            -internal-  value determined by interpolating obs pt to
    !                              binary ice mask field
    !                                for land/sea: iMask = -1
    !                                for ice:      iMask = 0
    ! mxLat            -internal-  number of grid pts in lat. direction for mesh
    ! mxLon            -internal-  number of grid pts in lon. direction for mesh
    ! rLatKm           -internal-  spacing desired between mesh grid points in km
    !                              along lat. direction
    ! rLonKm           -internal-  spacing desired between mesh grid points in km
    !                              along lon. direction
    ! dLat             -internal-  spacing between mesh grid points along lon. direction
    !                              in degrees computed from rLatKm
    ! dLon             -internal-  spacing between mesh grid points along lon. direction
    !                              in degrees computed from rLonKm
    ! rKmPerDeg        -internal- distance in km per degree
    !                               = Earth radius * PI/180.0
    !                               = 6371.01 km * PI/180.0
    !                               = 111.195 km
    ! nLat,nLon        -internal-  used to define the lat/lon of the grid pts of mesh
    ! obsLatBox        -internal-  lat values at all grid pts of mesh for all obs pts
    ! obsLonBox        -internal-  lon values at all grid pts of mesh for all obs pts
    ! latMesh          -internal-  lat values at all grid pts of mesh for 1 obs pt
    ! lonMesh          -internal-  lon values at all grid pts of mesh for 1 obs pt
    ! mgIntOb          -internal-  interpolated MG values at all grid pts of mesh for 1 obs pt
    ! lgIntOb          -internal-  interpolated LG values at all grid pts of mesh for 1 obs pt
    ! mgIntrp          -internal-  max. interpolated MG value on mesh for all obs pts
    ! lgIntrp          -internal-  max. interpolated LG value on mesh for all obs pts
    ! MGthresh         -internal-  maximum allowable land fraction for obs to be kept
    ! LGthresh         -internal-  maximum allowable ice  fraction for obs to be kept
    !--------------------------------------------------------------------

    implicit none

    ! Arguments
    integer,              intent(in)    :: numObsToProcess     ! Number of obs points to process
    real,                 intent(in)    :: obsLatitude(:)      ! Observation latitudes
    real,                 intent(in)    :: obsLongitude(:)     ! Observation longitudes
    integer,              intent(inout) :: landSeaQualifier(:) ! Land/sea indicator (0=land; 1=ocean)
    integer, allocatable, intent(out)   :: terrainType(:)      ! Terrain type qualifier
    logical, allocatable, intent(out)   :: waterObs(:)         ! Open water identifier for each obs

    ! Locals
    integer, parameter      :: mxLat=5
    integer, parameter      :: mxLon=5
    real,    parameter      :: LGthresh=0.01
    real,    parameter      :: MGthresh=0.01
    real,    parameter      :: pi=3.141592654
    real,    parameter      :: rKmPerDeg=111.195
    real,    parameter      :: rLatKm=40.0
    real,    parameter      :: rLonKm=40.0

    character(len=12)       :: etikxx
    character(len=1)        :: grtyp
    character(len=1)        :: grtyplg
    character(len=4)        :: nomvxx
    character(len=2)        :: typxx

    integer                 :: boxPointIndex
    integer                 :: ezQkDef
    integer                 :: ezSetOpt
    integer,           save :: gdId
    integer,           save :: gdIdlg
    integer                 :: gdllsval
    integer                 :: idum1, idum2, idum3, idum4
    integer                 :: idum5, idum6, idum7, idum8
    integer                 :: idum9, idum10, idum11, idum12
    integer                 :: idum13, idum14, idum15
    integer                 :: idum16, idum17, idum18
    integer                 :: ier
    integer                 :: ig1, ig2, ig3, ig4
    integer                 :: ig1lg, ig2lg, ig3lg, ig4lg
    integer                 :: iMask
    integer                 :: irec
    integer                 :: iUnGeo
    integer                 :: latIndex
    integer                 :: lonIndex
    integer                 :: ni, nj, nk
    integer                 :: nilg, njlg
    integer                 :: nLat
    integer                 :: nLon
    integer                 :: obsIndex

    logical,           save :: firstCall=.true.

    real, allocatable, save :: lg(:)
    real, allocatable, save :: mg(:)
    real, allocatable       :: latMesh(:)
    real, allocatable       :: lgIntOb(:)
    real, allocatable       :: lgIntrp(:)
    real, allocatable, save :: lonMesh(:)
    real, allocatable       :: mgIntOb(:)
    real, allocatable       :: mgIntrp(:)
    real, allocatable       :: obsLatBox(:,:)
    real, allocatable       :: obsLonBox(:,:)

    real                    :: dLat
    real                    :: dLon
    real                    :: rLatIndex
    real                    :: rLonIndex
    real                    :: xLat
    real                    :: xLatRad
    real                    :: xLon

    ! External functions
    integer, external       :: fclos
    integer, external       :: fnom
    integer, external       :: fstfrm
    integer, external       :: fstinf
    integer, external       :: fstlir
    integer, external       :: fstouv
    integer, external       :: fstprm

    ! Allocate space for arrays holding values on mesh grid pts.
    call utl_reAllocate(latMesh, mxLat*mxLon)
    call utl_reAllocate(lonMesh, mxLat*mxLon)
    call utl_reAllocate(mgIntOb, mxLat*mxLon)
    call utl_reAllocate(lgIntOb, mxLat*mxLon)
    call utl_reAllocate(obsLatBox, mxLat*mxLon, numObsToProcess)
    call utl_reAllocate(obsLonBox, mxLat*mxLon, numObsToProcess)
    call utl_reAllocate(terrainType, numObsToProcess)
    call utl_reAllocate(waterObs, numObsToProcess)

    if (firstCall) then

      firstCall = .false.

      ! Open FST file.
      iUnGeo = 0
      ier = fnom( iUnGeo,fileGlmg,'STD+RND+R/O',0 )
      ier = fstouv( iUnGeo,'RND' )

      ! Read MG field.
      irec = fstinf(iUnGeo,ni,nj,nk,-1,' ',-1,-1,-1,' ' ,'MG')
      if ( irec <  0 ) then
        call utl_abort('land_ice_mask_ssmis: The MG field is MISSING')
      end if

      call utl_reAllocate(mg, ni*nj)

      ier = fstlir(mg,iUnGeo,ni,nj,nk,-1,' ',-1,-1,-1,' ','MG')

      ier = fstprm(irec,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,    &
                   idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyp,ig1,ig2,  &
                   ig3,ig4,idum12,idum13,idum14,idum15,idum16,idum17,      &
                   idum18)

      ! Read LG field.

      irec = fstinf(iUnGeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'LG')
      if ( irec <  0 ) then
        call utl_abort('land_ice_mask_ssmis: The LG field is MISSING ')
      end if
      call utl_reAllocate(lg, nilg*njlg)
      ier = fstlir(lg,iUnGeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','LG')

      ier = fstprm(irec,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,          &
          &        idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyplg,ig1lg,ig2lg,  &
          &        ig3lg,ig4lg,idum12,idum13,idum14,idum15,idum16,idum17,        &
          &        idum18)

      gdId = ezQkDef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iUnGeo)
      gdIdlg = ezQkDef(nilg,njlg,grtyplg,ig1lg,ig2lg,ig3lg,ig4lg,iUnGeo)

      ier = fstfrm(iUnGeo)
      ier = fclos(iUnGeo)

    end if ! firstCall


    ! For each obs pt, define a grid of artificial pts surrounding it.

    nLat = (mxLat-1)/2
    nLon = (mxLon-1)/2

    dLat = rLatKm / rKmPerDeg
    do obsIndex = 1, numObsToProcess
      boxPointIndex = 0

      do latIndex = -nLat, nLat
        rlatIndex = float(latIndex)
        xLat = obsLatitude(obsIndex) + rLatIndex*dLat
        xLat = max( -90.0, min(90.0,xLat) )
        xLatRad = xLat*pi/180.0

        do lonIndex = -nLon, nLon
          dLon = rLonKm / ( rKmPerDeg*cos(xLatRad) )
          rLonIndex = float(lonIndex)
          boxPointIndex = boxPointIndex + 1
          xLon = obsLongitude(obsIndex) + rLonIndex*dLon
          if ( xLon < -180. ) xLon = xLon + 360.
          if ( xLon >  180. ) xLon = xLon - 360.
          if ( xLon <    0. ) xLon = xLon + 360.
          obsLatBox(boxPointIndex,obsIndex) = xLat
          obsLonBox(boxPointIndex,obsIndex) = xLon
        end do

      end do
    end do


    ! Interpolate values from MG and LG field to grid pts of mesh centred over each obs pt.
    ! Determine for each obs pt, the max interpolated MG and LG value within the box
    ! surrounding it.

    ier  = ezSetOpt('INTERP_DEGREE','LINEAR')
 
    call utl_reAllocate(mgIntrp, numObsToProcess)
    call utl_reAllocate(lgIntrp, numObsToProcess)

    mgIntrp(:) = 0.0
    lgIntrp(:) = 0.0
    do obsIndex = 1, numObsToProcess
 
      latMesh = obsLatBox(:,obsIndex)
      lonMesh = obsLonBox(:,obsIndex)

      ier  = gdllsval(gdId,mgIntOb,mg,latMesh,lonMesh,mxLat*mxLon)
      ier  = gdllsval(gdIdlg,lgIntOb,lg,latMesh,lonMesh,mxLat*mxLon)

      mgIntrp(obsIndex) = maxval(mgIntOb(:))
      lgIntrp(obsIndex) = maxval(lgIntOb(:))

    end do

    !  Initialize all obs as being over land and free of ice or snow.
    !  Determine which obs are over open water.

    waterObs(:) = .false.   ! not over open water
    terrainType(:) = -1             ! no ice or snow
    HEADER: do obsIndex = 1, numObsToProcess

      !    Determine for each obs pt, the value of iMask from the 0.1 deg binary file.
      !      - open land/sea: iMask = -1 (missing)
      !      - ice or snow:   iMask = 0  (from LG field --> binary ice mask)
      !    Define the terrain-type qualifier terrainType for each point based on the ice mask values.
      !       terrainType = -1  not defined (open water or snow free land)
      !                      0  sea-ice           (landSeaQualifier = 1)
      !                      1  snow-covered land (landSeaQualifier = 0)

      if (lgintrp(obsIndex) < LGthresh ) then
        iMask = -1
      else 
        iMask = 0
      end if 

      if ( iMask == 0 ) then  ! if ice or snow
        terrainType(obsIndex) = 1 - landSeaQualifier(obsIndex)
      end if

      !    If iMask is -1 (no ice/snow), and this is consistent with the model ice
      !    LG value (ie. < LGthresh), and the max MG value indicates ocean (ie.
      !    < MGthresh), then this is a WATEROBS point.

      if ( iMask == -1 .and. lgintrp(obsIndex) < LGthresh .and. mgintrp(obsIndex) < MGthresh ) then
      !!!! TO RUN WITHOUT BINARY ICE MASK CHECK (I.E. RELY ON LG ONLY TO REMOVE ICE PTS):
      !if ( lgintrp(obsIndex) < LGthresh .and. mgintrp(obsIndex) < MGthresh ) then
        waterObs(obsIndex) = .true.
      end if

      !   Modify land/sea quailifier if not consistent with waterObs (applies to remote small
      !   islands that should be treated as sea points (for RTTOV)):
      !     -- if waterObs=true, land/sea qualifier should be "sea" value (1)

      if ( waterObs(obsIndex) .and. (landSeaQualifier(obsIndex) == 0) ) landSeaQualifier(obsIndex) = 1

    end do HEADER

  end subroutine land_ice_mask_ssmis

  !--------------------------------------------------------------------------
  ! wentz_sfctype_ssmis  
  !--------------------------------------------------------------------------
  subroutine wentz_sfctype_ssmis(numObsToProcess, obsLatitude, obsLongitude, landSeaQualifier)
    ! :Purpose: Determine for each observation point the wentz surface value
    !           from the FST file wentz_surf.std.

    !           This file has a resolution of 0.25 deg, and it discriminates
    !           between the following surface types:
    !                      a) land  (wentz value = 0)
    !                      b) ice   (wentz value = 4)
    !                      c) sea   (wentz value = 5)
    !                      d) coast (wentz value = 6)
    !           This information is used to define the land/sea qualifier
    !           (element 008012) for each obs pt, which is needed by 3D-Var.
    !           Also, the land/sea qualifier is used later to flag obs pts
    !           over land and coast.
    !
    !--------------------------------------------------------------------
    !  Variable Definitions
    !  --------------------
    ! numObsToProcess  - input  -  number of input obs pts in record
    ! obsLatitude      - input  -  array of size ssbg_maxObsNum holding lat values for obs pts
    !                              in record plus undefined pts
    ! obsLongitude     - input  -  array of size ssbg_maxObsNum holding lon values for obs pts
    !                              in record plus undefined pts
    ! landSeaQualifier - output -  array holding land/sea qualifier values for all obs pts of record
    ! xLat             -internal-  array of size numObsToProcess holding lat values for obs pts in record
    ! xLon             -internal-  array of size numObsToProcess holding lon values for obs pts in record
    ! lm               -internal-  array of size 1440x720 holding gridded wentz surface values
    ! wenTyp           -internal-  array of size numObsToProcess holding wentz surface values interpolated to obs pts
    !--------------------------------------------------------------------

    implicit none

    ! Arguments
    integer, intent(in)               :: numObsToProcess      ! Number of obs points to process
    real,    intent(in)               :: obsLatitude(:)       ! Observation latitudes
    real,    intent(in)               :: obsLongitude(:)      ! Observation longitudes
    integer, intent(out), allocatable :: landSeaQualifier(:)  ! Land/sea indicator (0=land; 1=ocean)

    ! Locals
    character(len=12) :: etikxx
    character(len=1)  :: grtyp
    character(len=4)  :: nomvxx
    character(len=2)  :: typxx

    integer           :: ezQkDef
    integer           :: ezSetOpt
    integer           :: gdId
    integer           :: gdllsval
    integer           :: idum1, idum2, idum3, idum4
    integer           :: idum5, idum6, idum7, idum8
    integer           :: idum9, idum10, idum11, idum12
    integer           :: idum13, idum14, idum15
    integer           :: idum16, idum17, idum18
    integer           :: ier
    integer           :: ig1, ig2, ig3, ig4
    integer           :: irec
    integer           :: iUnIn
    integer           :: ni, nj, nk
    integer           :: obsIndex

    real, allocatable :: xLat(:)
    real, allocatable :: xLon(:)
    real, allocatable :: wenTyp(:)
    real, allocatable :: lm(:)

    ! External functions
    integer, external :: fclos
    integer, external :: fnom
    integer, external :: fstfrm
    integer, external :: fstinf
    integer, external :: fstlir
    integer, external :: fstouv
    integer, external :: fstprm

    ! Open Wentz surface field if first call
    iUnIn = 0
    ier = fnom( iUnIn,fileWentz,'STD+RND+R/O',0 )
    ier = fstouv( iUnIn,'RND' )

    irec = fstinf(iUnIn,ni,nj,nk,-1,' ',0,0,0,' ','LM')
    if ( irec <  0 ) then
      call utl_abort('wentz_sfctype_ssmis: The LM field is MISSING ')
    else
      call utl_reAllocate( lm, ni*nj )
      ier = fstlir(lm,iUnIn,ni,nj,nk,-1,' ',-1,-1,-1,' ','LM')
    end if

    ier = fstprm(irec,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,    &
        &        idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyp,ig1,ig2,  &
        &        ig3,ig4,idum12,idum13,idum14,idum15,idum16,idum17,      &
        &        idum18)

    ! Re-define the grid type.
    !    - L-grid: IG's SHOULD define grid spacing over a REGIONAL
    !              domain. Here, IG's are =0 though.

    !    - A-grid: is for global grid. See RPN documentation on
    !              grid-types for proper IG values.
    grtyp ='A'

    ! Transfer data from obsLatitude,obsLongitude to xLat,xLon.
    ! Ensure proper range for lon values.

    call utl_reAllocate( xLat, numObsToProcess)
    call utl_reAllocate( xLon, numObsToProcess)
    call utl_reAllocate( wenTyp, numObsToProcess)
    call utl_reAllocate( landSeaQualifier, numObsToProcess)

    landSeaQualifier(:) = 0
    xLat(:) = obsLatitude(1:numObsToProcess)
    xLon(:) = obsLongitude(1:numObsToProcess)
    do obsIndex = 1, numObsToProcess
      if ( xLon(obsIndex) < 0. ) xLon(obsIndex) = xLon(obsIndex) + 360.
    end do

    ! Interpolate values from LM field to all obs pts of record.
    ier  = ezSetOpt('INTERP_DEGREE','VOISIN')
    gdId = ezQkDef(ni,nj,grTyp,ig1,ig2,ig3,ig4,iUnIn)
    ier  = gdllsval(gdId,wenTyp,lm,xLat,xLon,numObsToProcess)

    ! Define the land/sea qualifier for each point based on wentz surface values.
    do obsIndex = 1, numObsToProcess
      if ( wenTyp(obsIndex) == 0. .or. wenTyp(obsIndex) == 6. ) then
        ! wentz = land/coast --> land
        landSeaQualifier(obsIndex) = 0
      else if ( wenTyp(obsIndex) == 4. .or. wenTyp(obsIndex) == 5. ) then
        ! wentz = sea/sea-ice --> sea
        landSeaQualifier(obsIndex) = 1
      else
        call utl_abort('wentz_sfctype_ssmis: Unexpected Wentz value ')
      end if
    end do
    ier = fstfrm(iUnIn)
    ier = fclos(iUnIn)

  end subroutine wentz_sfctype_ssmis

  !--------------------------------------------------------------------------
  ! ssbg_computeSsmisSurfaceType
  !--------------------------------------------------------------------------
  subroutine ssbg_computeSsmisSurfaceType(obsSpaceData)
    ! :Purpose: Compute surface type element and update obsSpaceData.

    implicit none

    ! Arguments
    type(struct_obs), intent(inout) :: obsSpaceData           ! ObsSpaceData object

    ! Locals
    real, parameter      :: satAzimuthAngle = 210.34
    real, parameter      :: satZenithAngle = 53.1

    integer, allocatable :: landSeaQualifier(:)
    integer              :: codtyp
    integer              :: headerIndex

    logical              :: ssmisDataPresent

    real                 :: obsLatitude(1)
    real                 :: obsLongitude(1)

    write(*,*) 'ssbg_computeSsmisSurfaceType: Starting'

    ssmisDataPresent=.false.
    call obs_set_current_header_list(obsSpaceData,'TO')

    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER0
      codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if ( tvs_isIdBurpInst(codtyp,'ssmis') ) then
        ssmisDataPresent = .true.
        exit HEADER0
      end if
    end do HEADER0
    
    if ( .not. ssmisDataPresent ) then
      write(*,*) 'WARNING: WILL NOT RUN ssbg_computeSsmisSurfaceType since no SSMIS DATA is found'
      return
    end if

    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER1
      obsLatitude(1)  = obs_headElem_r( obsSpaceData, OBS_LAT, headerIndex )
      obsLongitude(1) = obs_headElem_r( obsSpaceData, OBS_LON, headerIndex )
      obsLongitude(1) = obsLongitude(1)*MPC_DEGREES_PER_RADIAN_R8
      if( obsLongitude(1) > 180. ) obsLongitude(1) = obsLongitude(1) - 360.
      obsLatitude(1)  = obsLatitude(1) *MPC_DEGREES_PER_RADIAN_R8
      call wentz_sfctype_ssmis(1, obsLatitude, obsLongitude, landSeaQualifier)
      call obs_headSet_i(obsSpaceData, OBS_STYP, headerIndex, landSeaQualifier(1))
      call obs_headSet_r(obsSpaceData, OBS_SZA, headerIndex,satZenithAngle)
      call obs_headSet_r(obsSpaceData, OBS_AZA, headerIndex, satAzimuthAngle)
    end do HEADER1

    write(*,*) 'ssbg_computeSsmisSurfaceType: Finished'

  end subroutine ssbg_computeSsmisSurfaceType
  
  !--------------------------------------------------------------------------
  ! ssbg_grossValueCheck
  !--------------------------------------------------------------------------
  subroutine ssbg_grossValueCheck(numObsToProcess, obsTb, obsTbMin, obsTbMax, grossRej)
    ! :Purpose: Check obsTb for values that are missing or outside physical limits.

    implicit none

    ! Arguments
    integer, intent(in)               :: numObsToProcess  ! Number of obs points to process
    real,    intent(in)               :: obsTb(:)         ! Brightness temperature of observations
    real,    intent(in)               :: obsTbMin         ! Min(obsTb) threshold for rejection
    real,    intent(in)               :: obsTbMax         ! Max(obsTb) threshold for rejection
    logical, intent(out), allocatable :: grossRej(:)      ! Logical array of obs with gross error (obs to reject)

    ! Locals
    integer :: hiIndex
    integer :: loIndex
    integer :: obsIndex

    call utl_reAllocate(grossRej, numObsToProcess)
    
    grossRej(1:numObsToProcess) = .true.
    loIndex = 1
    do obsIndex = 1, numObsToProcess

      hiIndex = obsIndex*ssbg_maxNumChan
      if ( all( obsTb(loIndex:hiIndex) > obsTbMin ) .and. all( obsTb(loIndex:hiIndex) < obsTbMax ) ) then
        grossRej(obsIndex) = .false.
      end if
      loIndex = hiIndex + 1

    end do

  end subroutine ssbg_grossValueCheck

  !--------------------------------------------------------------------------
  ! ssbg_satqcSsmis
  !--------------------------------------------------------------------------
  subroutine ssbg_satqcSsmis(obsSpaceData, headerIndex, obsToReject)
    ! :Purpose: This program is applied as a first stage of processing to
    !           SSMIS data after it is received from UK MetOffice and
    !           organized into boxes by a program of Jose Garcia. The
    !           processing applied in this program includes:
    !             --  interpolate Wentz surface land mask to each obs pt
    !                 (nearest neighbour) to define land/sea qualifier (008012)
    !             --  interpolate binary ice mask to each obs pt (nearest
    !                 neighbour) to define terrain-type element (013039) where
    !                 0 = sea ice and 1 = snow-covered land
    !             --  interpolate model MG and LG fields to a grid surrounding each obs
    !                 pt to identify obs that are over open water, far from coast/ice
    !             --  identify those obs for which the UKMO rain marker
    !                 is ON (ie. 020029 = 1) indicating poor quality
    !             --  apply a cloud filter to identify those obs in cloudy regions;
    !                 write CLW and IWV (over ocean) to output BURP file
    !             --  re-write data to output BURP file while modifying flags
    !                 for those obs which are not over open water, or have been
    !                 identified in rain/cloud areas, or are of poor quality
    !             --  define satellite zenith angle element (007024) and add
    !                 this and land/sea qualifier and terrain-type elements
    !                 to the output file

    implicit none

    ! Arguments
    type(struct_obs),     intent(inout) :: obsSpaceData           ! ObsSpaceData object
    integer,              intent(in)    :: headerIndex            ! Current header index
    logical, allocatable, intent(out)   :: obsToReject(:)         ! Observations that will be rejected
    
    ! Locals
    ! arrays to get from obsspacedata
    character(len=9)           :: burpFileSatId

    integer, allocatable       :: landSeaQualifier(:)
    integer, allocatable       :: terrainType(:)
    integer, allocatable       :: ukRainObs(:)

    real,    allocatable       :: obsLatitude(:)
    real,    allocatable       :: obsLongitude(:)
    real,    allocatable       :: obsTb(:)
    real,    allocatable       :: rclw(:)
    real,    allocatable       :: riwv(:)
    real,    allocatable       :: satZenithAngle(:)
    real,    allocatable       :: scatW(:)

    ! temporary arrays
    integer, allocatable       :: ier(:)

    integer                    :: actualNumChannel
    integer                    :: bodyIndex
    integer                    :: bodyIndexBeg
    integer                    :: bodyIndexEnd
    integer                    :: codtyp
    integer                    :: currentChannelNumber
    integer                    :: headerCompt
    integer                    :: instr
    integer                    :: instrId
    integer                    :: obsIndex
    integer                    :: iSat
    integer                    :: iSatId
    integer,              save :: numLandObs
    integer,              save :: numUkBadObs
    integer,              save :: numGrossObs
    integer,              save :: numCloudyObs
    integer,              save :: numBadIWVObs
    integer,              save :: numPrecipObs
    integer,              save :: numTotFilteredObs
    integer,              save :: numLandScatObs
    integer,              save :: numSeaScatObs
    integer,              save :: numDryIndexObs
    integer,              save :: numSeaIceObs
    integer,              save :: numScatPrecipObs
    integer,              save :: numObsF16
    integer,              save :: numObsF17
    integer,              save :: numObsF18
    integer                    :: numObsToProcess
    integer                    :: platf
    integer                    :: platfId
    integer                    :: refPosObs
    integer                    :: sensorIndex

    logical, allocatable       :: cloudObs(:)
    logical, allocatable       :: grossRej(:)
    logical, allocatable       :: iwvReject(:)
    logical, allocatable       :: precipObs(:)
    logical, allocatable       :: rainDetectionUKMethod(:)
    logical, allocatable       :: waterObs(:)

    logical,              save :: ifFirstCall = .true.
    logical                    :: sensorIndexFound
    logical                    :: ssmisDataPresent

    real,    allocatable       :: amsubDrynessIndex(:)
    real,    allocatable       :: scatL(:)
    real,    allocatable       :: ztb91(:)
    real,    allocatable       :: ztb150(:)
    real,    allocatable       :: ztb_amsub3(:)
    real,    allocatable       :: ztb_amsub5(:)

    ! Check if its ssmis data:
    codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
    ssmisDataPresent = tvs_isIdBurpInst(codtyp,'ssmis')

    if ( .not. ssmisDataPresent ) then
      write(*,*) 'WARNING: WILL NOT RUN ssbg_satqcSsmis since no SSMIS DATA is found'
      return
    end if


    ! find tvs_sensor index corresponding to current obs

    platf = obs_headElem_i( obsSpaceData, OBS_SAT, headerIndex )
    instr = obs_headElem_i( obsSpaceData, OBS_INS, headerIndex )

    call tvs_mapSat( platf, platfId, iSat )
    call tvs_mapInstrum( instr, instrId )

    sensorIndexFound = .false.
    HEADER0: do sensorIndex = 1, tvs_nsensors
      if ( platfId == tvs_platforms(sensorIndex)  .and. &
           iSat    == tvs_satellites(sensorIndex) .and. &
           instrId == tvs_instruments(sensorIndex) ) then
          sensorIndexFound = .true.
        exit HEADER0
      end if
    end do HEADER0
    if ( .not. sensorIndexFound ) call utl_abort('ssbg_satqcSsmis: sensor Index not found')

    ! find actual Number of channels
    actualNumChannel = tvs_coefs(sensorIndex)%coef%fmv_ori_nchn

    headerCompt = 1
    numObsToProcess = 1

    ! Allocate intent out arrays 

    call utl_reAllocate(obsToReject, numObsToProcess*actualNumChannel)

    ! Allocate Fortran working arrays 

    call utl_reAllocate(amsubDrynessIndex, numObsToProcess)
    call utl_reAllocate(cloudObs, numObsToProcess)
    call utl_reAllocate(grossRej, numObsToProcess)
    call utl_reAllocate(ier, numObsToProcess)
    call utl_reAllocate(iwvReject, numObsToProcess)
    call utl_reAllocate(precipObs, numObsToProcess)
    call utl_reAllocate(rainDetectionUKMethod, numObsToProcess)
    call utl_reAllocate(scatL, numObsToProcess)
    call utl_reAllocate(waterObs, numObsToProcess)
    call utl_reAllocate(ztb91, numObsToProcess)
    call utl_reAllocate(ztb150, numObsToProcess)
    call utl_reAllocate(ztb_amsub3, numObsToProcess)
    call utl_reAllocate(ztb_amsub5, numObsToProcess)

    ! ELEMENTS FROM OBSSPACEDATA
    ! Allocate Header elements
    call utl_reAllocate(landSeaQualifier, numObsToProcess)
    call utl_reAllocate(obsLatitude, numObsToProcess)
    call utl_reAllocate(obsLongitude, numObsToProcess)
    call utl_reAllocate(rclw, numObsToProcess)
    call utl_reAllocate(riwv, numObsToProcess)
    call utl_reAllocate(satZenithAngle, numObsToProcess)
    call utl_reAllocate(scatW, numObsToProcess)
    call utl_reAllocate(terrainType, numObsToProcess)
    call utl_reAllocate(ukRainObs, numObsToProcess)
    ! Allocate Body elements
    call utl_reAllocate(obsTb, numObsToProcess*actualNumChannel)
    !initialization
    obsTb(:) = ssbg_realMissing
    riwv(:) = ssbg_realMissing
    scatW(:) = ssbg_realMissing
    
    ! Lecture dans obsspacedata

    burpFileSatId                 = obs_elem_c    ( obsSpaceData, 'STID'  , headerIndex )
    rclw(headerCompt)             = obs_headElem_r( obsSpaceData, OBS_CLWO, headerIndex )
    ukRainObs(headerCompt)        = obs_headElem_i( obsSpaceData, OBS_RAIN, headerIndex )
    landSeaQualifier(headerCompt) = obs_headElem_i( obsSpaceData, OBS_STYP, headerIndex )
    terrainType(headerCompt)      = obs_headElem_i( obsSpaceData, OBS_TTYP, headerIndex )
    ! If terrain type is missing, set it to -1 for the QC programs
    if (terrainType(headerCompt) ==  99) terrainType(headerCompt) = -1
    obsLatitude(headerCompt)      = obs_headElem_r( obsSpaceData, OBS_LAT, headerIndex )
    obsLongitude(headerCompt)     = obs_headElem_r( obsSpaceData, OBS_LON, headerIndex )
    satZenithAngle(headerCompt)   = obs_headElem_r( obsSpaceData, OBS_SZA, headerIndex )
    ! Convert lat/lon to degrees
    obsLatitude(headerCompt)  = obsLatitude(headerCompt) *MPC_DEGREES_PER_RADIAN_R8
    obsLongitude(headerCompt) = obsLongitude(headerCompt)*MPC_DEGREES_PER_RADIAN_R8
    if( obsLongitude(headerCompt) > 180. ) obsLongitude(headerCompt) = obsLongitude(headerCompt) - 360.

    ! To read body elements
    bodyIndexBeg = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    bodyIndexEnd = bodyIndexBeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) -1

    BODY: do bodyIndex = bodyIndexbeg, bodyIndexEnd
      currentChannelNumber = nint(obs_bodyElem_r( obsSpaceData, OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      obsTb(currentChannelNumber) = obs_bodyElem_r( obsSpaceData,  OBS_VAR, bodyIndex )
    end do BODY

    ! initialization
    if (ifFirstCall) then
      ifFirstCall = .false.
      numBadIWVObs = 0
      numCloudyObs = 0
      numDryIndexObs = 0
      numGrossObs = 0
      numLandObs = 0
      numLandScatObs = 0
      numObsF16 = 0
      numObsF17 = 0
      numObsF18 = 0
      numPrecipObs = 0
      numScatPrecipObs = 0
      numSeaIceObs = 0
      numSeaScatObs = 0
      numTotFilteredObs = 0
      numUkBadObs = 0
    end if 
    waterObs(:) = .false.

    ! Record the total number of obs pts read for each satellite.
    ! Set the satellite ID number.
    select case (burpFileSatId(2:7) )
    case ( 'DMSP16' )
      numObsF16 = numObsF16 + numObsToProcess
      iSatId = 1
    case ( 'DMSP17' )
      numObsF17 = numObsF17 + numObsToProcess
      iSatId = 2
    case ( 'DMSP18' )
      numObsF18 = numObsF18 + numObsToProcess
      iSatId = 3
    end select

    !--------------------------------------------------------------------
    ! Determine which obs pts are over open water (i.e NOT near coasts or
    ! over/near land/ice.
    !    - using model MG field, plus 0.1 deg ice mask and model LG field
    ! Also, define the terrain-type qualifier values (itt).
    ! Also, ensure land/sea qualifier (ilq) = 1 (sea) for points identified
    !   as open water by waterobs array.
    !--------------------------------------------------------------------

    call land_ice_mask_ssmis(numObsToProcess, obsLatitude, obsLongitude, landSeaQualifier, &
                              terrainType, waterObs)

    !--------------------------------------------------------------------
    ! Determine which obs pts have been flagged for rain (using UKMO method of identifying
    ! bad quality SSMIS data for DMSP16). First, initialize all obs as good.
    !--------------------------------------------------------------------

    rainDetectionUKMethod(:) = .false.
    where ( ukRainObs == 1 ) rainDetectionUKMethod = .true.

    !--------------------------------------------------------------------
    ! Check for values of TB that are missing or outside physical limits.
    ! **NOTE: REJECT ALL CHANNELS IF ONE IS FOUND TO BE BAD.
    !--------------------------------------------------------------------

    grossRej(:)  = .false.
    call ssbg_grossValueCheck(numObsToProcess, obsTb, 50., 400., grossRej)

    !--------------------------------------------------------------------
    ! Apply a CLW regression technique to determine cloudy obs pts.
    ! Apply a IWV regression technique to determine IWV obs out of bounds.
    ! Detect Sea_Ice percent from Tb(Ta) and set waterobs() flag to .FALSE. if
    !  amount > 70%
    ! To begin, assume that all obs are good.
    !--------------------------------------------------------------------

    cloudobs(:)  = .false.
    iwvReject(:) = .false.
    precipobs(:) = .false.

    ! --------------------------------------------------------------------

    call cld_filter_fweng(numObsToProcess, obsTb, algOption, waterObs, grossRej, &
                          cloudObs, iwvReject, precipObs, rclw, riwv, iSatId,    &
                          obsLatitude, numSeaIceObs)
    !
    !  NOTE: rclw, riwv missing value is zmisg=9.9e09
    !    --> if ( (clw > clw_thresh) .or. precipobs(obsIndex) ) cloudobs(obsIndex) = .true.
    !    --> if ( clw == zmisg )  cloudobs(obsIndex) = .true.
    !
    !   Extract Tb for channels 18 (91H GHz)  and 8 (150H GHz) for Bennartz SI
    !   Extract Tb for channels 11 (AMSU-B 3) and 9 (AMSU-B 5) for Dryness Index (DI)

    refPosObs = 1
    do obsIndex = 1, numObsToProcess
      ztb91(obsIndex)      = obsTb(refPosObs+17)
      ztb150(obsIndex)     = obsTb(refPosObs+7)
      ztb_amsub3(obsIndex) = obsTb(refPosObs+10)
      ztb_amsub5(obsIndex) = obsTb(refPosObs+8)
      refPosObs = obsIndex*ssbg_maxNumChan + 1
    end do

    !---------------------------------------------------------------------------------
    ! Compute the Bennartz scattering index (SI) for each point from channels 8 and 18
    ! (detects precipitation affected AMSU-B-like channel radiances)
    ! -- using constant satellite zenith angle = 53.1
    ! -- SSMIS channel 18 (91.655 GHz) is used for AMSU-B like channel 1 (89 GHz)
    ! -- SSMIS channel  8 is used for AMSU-B like channel 2 (both at 150 GHz)
    !---------------------------------------------------------------------------------
    call bennartz(ier, numObsToProcess, ztb91, ztb150, satZenithAngle, landSeaQualifier, & 
                 scatL, scatW)

    !--------------------------------------------------------------------
    ! Compute the AMSU-B Dryness Index for all points
    !--------------------------------------------------------------------

    where ( (ztb_amsub3 /= ssbg_realMissing) .and. (ztb_amsub5 /= ssbg_realMissing) )
      amsubDrynessIndex = ztb_amsub3 - ztb_amsub5
    elsewhere
      amsubDrynessIndex = ssbg_realMissing 
    end where

    !--------------------------------------------------------------------
    ! Review all the checks previously made to determine which obs are to be accepted
    ! for assimilation and which are to be flagged for exclusion. First, initialize
    ! all obs for inclusion.
    !   ukBadObs()  = .true. if UKMetO QC sets "rain" flag: solar-intrusion or other anomaly
    !   grossRej()  = .true. if any channel had a gross error at the point
    !   cloudObs()  = .true. if CLW > 0.01 kg/m**2 or precip (over water)
    !   precipObs() = .true. if precip. detected (CLW=missing)
    !   waterObs()  = .true. if open water point (far from land or sea ice)
    !   iwvReject() = .true. if IWV > 80 kg/m**2
    !--------------------------------------------------------------------

    obsToReject(:) = .false.
    HEADER1: do obsIndex = 1, numObsToProcess
      
      !      Reject all channels if UKMet rain flag or gross Tb error detected   
      if ( rainDetectionUKMethod(obsIndex) .or. grossRej(obsIndex) ) then
        !      "BAD" Observations       
        obsToReject(:) = .true.
        numTotFilteredObs = numTotFilteredObs + 1

      else     ! if ( ukBadObs(obsIndex) .or. grossRej(obsIndex) )
        !      "GOOD" Observations       
        !-----------------------------------------------------------------------------------       
        !      OVER LAND OR SEA-ICE, 
        !         -- reject lower tropospheric channels and imager channels
        !         -- check Bennartz SI and DI for AMSU-B like channels
        !----------------------------------------------------------------------------------- 
        !        -- CLW/PRECIP not determined over land
        !        -- surface emissivity effects lower tropospheric channels     
        if  ( .not. waterObs(obsIndex) ) then
          obsToReject(1:ipc) = .true.    ! AMSU-A 3-5(6)
          obsToReject(12:18) = .true.    ! SSMI-like imager 1-7
          obsToReject(8:9) = .true.      ! AMSU-B 2,5
          !        Check BSI and DI for AMSU-B 3-4
          !         Bennartz Scattering Index
          !         Land point           
          if ( scatL(obsIndex) > 0.0  .and. scatL(obsIndex) /= ssbg_realMissing )   obsToReject(10:11) = .true.
          !         Sea-ice point 
          if ( scatW(obsIndex) > 40.0 .and. scatW(obsIndex) /= ssbg_realMissing )   obsToReject(10:11) = .true.
          !         Missing scattering index
          if ( scatW(obsIndex) == ssbg_realMissing .and. scatL(obsIndex) == ssbg_realMissing ) obsToReject(10:11) = .true.
          if ( any(obsToReject(10:11))) then
            numLandScatObs = numLandScatObs + 1
          end if
          !         Dryness index           
          if ( amsubDrynessIndex(obsIndex) > 0.0 ) then
            obsToReject(11) = .true.
          end if
          if ( amsubDrynessIndex(obsIndex) > -10.0 ) then
            obsToReject(10) = .true.
          end if
          if ( amsubDrynessIndex(obsIndex) > -20.0 ) then
            obsToReject(9) = .true.
          end if
          if ( amsubDrynessIndex(obsIndex) > -10.0 ) numDryIndexObs = numDryIndexObs + 1
          numTotFilteredObs = numTotFilteredObs + 1
        end if
        !-----------------------------------------------------------------------------------         
        !      OVER WATER, 
        !        -- reject tropospheric channels and imager channels if cloudy, 
        !           precip, or excessive IWV
        !        -- check Bennartz SI for AMSU-B like channels
        !-----------------------------------------------------------------------------------    
        if  ( waterObs(obsIndex) ) then
          if (cloudObs(obsIndex) .or. iwvReject(obsIndex))  then
            if ( iwvReject(obsIndex) ) then
              !              SSMI-like imager channels
              obsToReject(12:18) = .true.
              !             AMSU-A like channels 3-5(6)
              obsToReject(1:ipc) = .true.
            else  !  ----- CLOUDY OBSERVATION (OR MISSING CLOUD) ------
              !              SSMI-like imager channels
              obsToReject(12:18) = .true.
              !              AMSU-A like channels 3-5(6)
              if ( rclw(obsIndex) > clw_amsu_rej ) then
                obsToReject(2:ipc) = .true.
              end if
              if ( rclw(obsIndex) > clw_amsu_rej_ch3  ) then
                obsToReject(1) = .true.
              end if
            end if
            numTotFilteredObs = numTotFilteredObs + 1
          end if
          !        Check BSI for AMSU B channels 2-5 for all water obs
          !        Bennartz Scattering Index ( AMSU-B 2-5)
          !        Open water point
          if ( scatW(obsIndex) > 15.0 .and. scatW(obsIndex) /= ssbg_realMissing ) then
            obsToReject(8:11) = .true.
            numSeaScatObs = numSeaScatObs + 1
            if ( precipObs(obsIndex) ) numScatPrecipObs = numScatPrecipObs + 1
          end if

        end if   ! if waterObs

      end if   ! if ( ukBadObs(obsIndex) .or. grossRej(obsIndex) [end else] )

      if ( .not. waterObs(obsIndex) ) numLandObs  = numLandObs  + 1
      if ( rainDetectionUKMethod(obsIndex) )  numUkBadObs = numUkBadObs + 1
      if ( grossRej(obsIndex) )  numGrossObs = numGrossObs + 1
      if ( cloudObs(obsIndex)  .and. waterObs(obsIndex) ) numCloudyObs  = numCloudyObs  + 1
      if ( iwvReject(obsIndex) .and. waterObs(obsIndex) .and.   &
      &         (.not. cloudObs(obsIndex)) ) numBadIWVObs  = numBadIWVObs  + 1
      if ( precipObs(obsIndex) .and. waterObs(obsIndex) ) then
        numPrecipObs = numPrecipObs + 1
      end if
    end do HEADER1
      
    !-------------------------------------------------------------------------------
    ! Update ObsspaceData with the computed values of: ZLQ, ZTT, IWV and scatW
    !-------------------------------------------------------------------------------
    call obs_headSet_i(obsSpaceData, OBS_STYP, headerIndex, landSeaQualifier(1))
    call obs_headSet_i(obsSpaceData, OBS_TTYP, headerIndex, terrainType(1))
    call obs_headSet_r(obsSpaceData, OBS_SCAT, headerIndex, scatW(1))
    call obs_headSet_r(obsSpaceData, OBS_CLWO, headerIndex, rclw(1))
    call obs_headSet_r(obsSpaceData, OBS_IWV,  headerIndex, riwv(1))

    if (headerIndex == obs_numHeader(obsSpaceData)) then
      write(*,*) '*******************************************************************'
      write(*,*) '******************* SATQC PROGRAM STATS****************************'
      write(*,*) '*******************************************************************'
      write(*,*)
      write(*,*) 'Number of cloudy obs =    ', numCloudyObs
      write(*,*) 'Number of precip obs =    ', numPrecipObs
      write(*,*) 'Number of gross rej obs = ',  numGrossObs
      write(*,*) 'Number of land  obs =     ',numLandObs
      write(*,*) 'Number of dry obs =       ', numDryIndexObs
      write(*,*) 'Number of Sea Ice obs =   ',numSeaIceObs
      write(*,*) 'Number of F16 obs =       ', numObsF16
      write(*,*) 'Number of F17 obs =       ', numObsF17
      write(*,*) 'Number of F18 obs =       ', numObsF18
      write(*,*) '*******************************************************************'
    end if

  end subroutine ssbg_satqcSsmis 

  !--------------------------------------------------------------------------
  ! ssbg_updateObsSpaceAfterSatQc
  !--------------------------------------------------------------------------
  subroutine ssbg_updateObsSpaceAfterSatQc(obsSpaceData, headerIndex, obsToReject) 
    !:Purpose:      Update obspacedata variables (obstTB and obs flags) after QC

    implicit none

    !Arguments
    type(struct_obs),     intent(inout) :: obsSpaceData           ! ObsSpaceData object
    integer,              intent(in)    :: headerIndex            ! Current header index
    logical,              intent(in)    :: obsToReject(:)         ! Observations that will be rejected

    ! Locals
    integer, allocatable                :: obsFlags(:)
    integer, allocatable                :: obsGlobalFlag(:)
    integer, allocatable                :: satScanPosition(:)

    integer                             :: actualNumChannel
    integer                             :: bodyIndex
    integer                             :: bodyIndexBeg
    integer                             :: bodyIndexEnd
    integer                             :: channelIndex
    integer                             :: currentChannelNumber
    integer                             :: dataIndex
    integer                             :: headerCompt
    integer                             :: instr
    integer                             :: instrId
    integer                             :: iSat
    integer                             :: numObsToProcess
    integer                             :: obsIndex
    integer                             :: platf
    integer                             :: platfId
    integer                             :: sensorIndex

    logical                             :: sensorIndexFound

    ! find tvs_sensor index corresponding to current obs

    platf = obs_headElem_i( obsSpaceData, OBS_SAT, headerIndex )
    instr = obs_headElem_i( obsSpaceData, OBS_INS, headerIndex )

    call tvs_mapSat( platf, platfId, iSat )
    call tvs_mapInstrum( instr, instrId )

    sensorIndexFound = .false.
    HEADER0: do sensorIndex = 1, tvs_nsensors
      if ( platfId == tvs_platforms(sensorIndex)  .and. &
           iSat    == tvs_satellites(sensorIndex) .and. &
           instrId == tvs_instruments(sensorIndex) ) then
          sensorIndexFound = .true.
        exit HEADER0
      end if
    end do HEADER0
    if ( .not. sensorIndexFound ) call utl_abort('ssbg_updateObsSpaceAfterSatQc: sensor Index not found')

    ! find actual Number of channels
    actualNumChannel = tvs_coefs(sensorIndex)%coef%fmv_ori_nchn

    headerCompt = 1
    numObsToProcess = 1
    ! Allocate Header elements
    call utl_reAllocate(obsGlobalFlag, numObsToProcess)
    call utl_reAllocate(satScanPosition, numObsToProcess)
    
    ! Allocation
    call utl_reAllocate(obsFlags, numObsToProcess*actualNumChannel)

    ! Read elements in obsspace

    obsGlobalFlag(headerCompt)   = obs_headElem_i( obsSpaceData, OBS_ST1, headerIndex )
    satScanPosition(headerCompt) = obs_headElem_i( obsSpaceData, OBS_FOV, headerIndex )

    bodyIndexBeg = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    bodyIndexEnd = bodyIndexBeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) -1

    BODY2: do bodyIndex = bodyIndexbeg, bodyIndexEnd
      currentChannelNumber = nint(obs_bodyElem_r( obsSpaceData, OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      obsFlags(currentChannelNumber) = obs_bodyElem_i( obsSpaceData, OBS_FLG, bodyIndex )
    end do BODY2
    
    !------------------------------------------------------------------ 
    ! 1 - Mark global flags if any obsToReject
    !------------------------------------------------------------------ 
    do obsIndex = 1, numObsToProcess
      if ( any(obsToReject) ) obsGlobalFlag(obsIndex) = ibset(obsGlobalFlag(obsIndex), 6)
    end do 

    !------------------------------------------------------------------ 
    ! 2 - Mark obs flags for each value of obsToReject
    !------------------------------------------------------------------ 
    dataIndex = 0
    do obsIndex = 1, numObsToProcess
      do channelIndex = 1, actualNumChannel
        dataIndex = dataIndex+1 
        if (resetQc) obsFlags(dataIndex) = 0
        if (obsToReject(dataIndex)) obsFlags(dataIndex) = ibset(obsFlags(dataIndex),7)
      end do 
    end do

    !-----------------------------------------------------------------
    !    Subtract 270 from FOV values (element 005043).
    !-----------------------------------------------------------------
    
    do obsIndex = 1, numObsToProcess
      if (satScanPosition(obsIndex) > 270) satScanPosition(obsIndex) = satScanPosition(obsIndex) - 270
    end do
 
    ! write elements in obsspace
    call obs_headSet_i(obsSpaceData, OBS_FOV, headerIndex, satScanPosition(1))
    call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex, obsGlobalFlag(1))

    bodyIndexBeg = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    bodyIndexEnd = bodyIndexBeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) -1

    BODY: do bodyIndex = bodyIndexBeg, bodyIndexEnd
      currentChannelNumber=nint(obs_bodyElem_r( obsSpaceData, OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, obsFlags(currentChannelNumber))
    end do BODY

  end subroutine ssbg_updateObsSpaceAfterSatQc 

  !--------------------------------------------------------------------------
  ! ssbg_inovqcSsmis
  !--------------------------------------------------------------------------
  subroutine ssbg_inovqcSsmis(obsSpaceData, headerIndex, flagsInovQc)
!  :Purpose: Identify those observations in SSMIS data that have O-P
!            values greater than a threshold proportional to known
!            standard deviations (computed when the bias correction
!            coefficients were derived). The flags of these observations
!            are adjusted accordingly (ie bit 9 switched ON).
!            Also,
!              -- flag channels for systematic rejection based on UTIL
!                 value in stats_*_errtot file or because flag bit 6 OFF
!                 (uncorrected data)
!              -- reject sets of AMSU-like channels based on O-P for
!                 a single channel
!              -- reject selected AMSU-like channels over land when
!                 model surface height exceeds a specified limit
!                 (topography check)
!
!------------------------------------------------------------------

    implicit none

    ! Arguments
    type(struct_obs),     intent(inout) :: obsSpaceData      ! ObsSpaceData object
    integer,              intent(in)    :: headerIndex       ! Current header index
    integer, allocatable, intent(out)   :: flagsInovQc(:)    ! Flags for assimilation/rejection of obs

    ! Locals
    character(len=9)     :: burpFileSatId            ! Satellite ID

    integer, allocatable :: obsChannels(:)           ! channel numbers
    integer, allocatable :: obsFlags(:)              ! data flags

    integer              :: actualNumChannel         ! actual Num channel
    integer              :: bodyIndex
    integer              :: bodyIndexBeg
    integer              :: bodyIndexEnd 
    integer              :: codtyp                   ! code type
    integer              :: channelIndex
    integer              :: currentChannelNumber
    integer              :: headerCompt
    integer              :: instr
    integer              :: instrId
    integer              :: iSat
    integer              :: numObsToProcess          ! Number of obs points to process
    integer              :: platf
    integer              :: platfId
    integer              :: sensorIndex              ! find tvs_sensor index corresponding to current obs

    logical              :: sensorIndexFound
    logical              :: ssmisDataPresent

    real   , allocatable :: modelInterpTer(:)        ! topo in standard file interpolated to obs point
    real   , allocatable :: obsLatitude(:)           ! obs. point latitude
    real   , allocatable :: obsLongitude(:)          ! obs. point longitude
    real   , allocatable :: ompTb(:)                 ! OMP values

    !-------------------------------------------------------------------------
    ! 1) INOVQC begins
    !-------------------------------------------------------------------------

    ! Check if its ssmis data:
    ssmisDataPresent = .false.
    codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
    if ( tvs_isIdBurpInst(codtyp,'ssmis') ) then
      ssmisDataPresent = .true.
    end if

    if ( .not. ssmisDataPresent ) then
      write(*,*) 'WARNING: WILL NOT RUN ssbg_inovqcSsmis since no SSMIS DATA were found'
      return
    end if

    ! find tvs_sensor index corresponding to current obs

    platf = obs_headElem_i( obsSpaceData, OBS_SAT, headerIndex )
    instr = obs_headElem_i( obsSpaceData, OBS_INS, headerIndex )

    call tvs_mapSat( platf, platfId, iSat )
    call tvs_mapInstrum( instr, instrId )

    sensorIndexFound = .false.
    HEADER1: do sensorIndex = 1, tvs_nsensors
      if ( platfId == tvs_platforms(sensorIndex)  .and. &
           iSat    == tvs_satellites(sensorIndex) .and. &
           instrId == tvs_instruments(sensorIndex) ) then
        sensorIndexFound = .true.
        exit HEADER1
      end if
    end do HEADER1
    if ( .not. sensorIndexFound ) call utl_abort('ssbg_inovqcSsmis: sensor Index not found')

    !--------------------------------------------------------------------
    ! 2) Allocating arrays
    !--------------------------------------------------------------------

    ! find actual Number of channels
    actualNumChannel = tvs_coefs(sensorIndex)%coef%fmv_ori_nchn

    headerCompt = 1
    numObsToProcess = 1

    ! ELEMENTS FROM OBSSPACEDATA
    ! Allocate header elements
    call utl_reAllocate(obsLatitude, numObsToProcess)
    call utl_reAllocate(obsLongitude, numObsToProcess)
    ! Allocate Body elements
    call utl_reAllocate(obsChannels, numObsToProcess*actualNumChannel)
    call utl_reAllocate(obsFlags, numObsToProcess*actualNumChannel)
    call utl_reAllocate(ompTb, numObsToProcess*actualNumChannel)
    ! Allocate intent out array
    call utl_reAllocate(flagsInovQc, numObsToProcess*actualNumChannel)

    ! Lecture dans obsspacedata
    burpFileSatId             = obs_elem_c    ( obsSpaceData, 'STID' , headerIndex )
    obsLatitude(headerCompt)  = obs_headElem_r( obsSpaceData, OBS_LAT, headerIndex )
    obsLongitude(headerCompt) = obs_headElem_r( obsSpaceData, OBS_LON, headerIndex )

    bodyIndexBeg = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    bodyIndexEnd = bodyIndexBeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) -1

    BODY: do bodyIndex = bodyIndexbeg, bodyIndexEnd
      currentChannelNumber = nint(obs_bodyElem_r( obsSpaceData, OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      ompTb(currentChannelNumber)    = obs_bodyElem_r( obsSpaceData, OBS_OMP, bodyIndex )
      obsFlags(currentChannelNumber) = obs_bodyElem_i( obsSpaceData, OBS_FLG, bodyIndex )
    end do BODY

    do channelIndex=1, actualNumChannel
      obsChannels(channelIndex) = channelIndex+tvs_channelOffset(sensorIndex)
    end do

    ! Convert lat/lon to degrees
    obsLatitude(headerCompt)  = obsLatitude(headerCompt) *MPC_DEGREES_PER_RADIAN_R8
    obsLongitude(headerCompt) = obsLongitude(headerCompt)*MPC_DEGREES_PER_RADIAN_R8
    if( obsLongitude(headerCompt) > 180. ) obsLongitude(headerCompt) = obsLongitude(headerCompt) - 360.

    !--------------------------------------------------------------------
    ! 3) Extract model topography from the input GEO file and interpolate
    !    to all observation points in the box
    !--------------------------------------------------------------------

    call ssbg_readGeophysicFieldsAndInterpolate(obsLatitude, obsLongitude, modelInterpTer)

    !--------------------------------------------------------------------
    ! 4) Perform quality control on the data, checking O-P and topography
    !--------------------------------------------------------------------

    call check_stddev(obsChannels, ompTb, flagsInovQc, actualNumChannel, numObsToProcess, &
        &             sensorIndex, burpFileSatId, obsFlags)

    call check_topo(modelInterpTer, flagsInovQc, actualNumChannel, numObsToProcess)

  end subroutine ssbg_inovqcSsmis

  !--------------------------------------------------------------------------
  ! check_stddev
  !--------------------------------------------------------------------------
  subroutine check_stddev(obsChannels, ompTb, flagsInovQc, actualNumChannel, numObsToProcess, &
       &                  sensorIndex, burpFileSatId, obsFlags)
    ! :Purpose: Perform quality control on the radiances by analysing the
    !           magnitude of the residuals.

    !------------------------------------------------------------------
    ! Variable Definitions:
    ! ---------------------
    !   obsChannels      - input  -  channel numbers (residuals)
    !   ompTb            - input  -  radiance residuals (o-p)
    !   obsFlags         - input  -  radiance data flags (bit 7 on for data that
    !                                will not be assimilated)
    !   flagsInovQc      - output -  quality contol indicator for each channel of each
    !                                observation point
    !                                =0  ok
    !                                =1  not checked because FLAG bit 7 ON,
    !                                =2  reject by UTIL value or FLAG bit 6 OFF (data not bias corrected),
    !                                =3  reject by rogue check,
    !                                =4  reject by both UTIL value and rogue check.
    !                                ==> individual channels may be rejected from each obs pt
    !   actualNumChannel - input  -  number of residual channels
    !   numObsToProcess  - input  -  number of groups of NVAL*NELE
    !   sensorIndex      - input  -  number of satellite (index # --> 1-nsat)
    !   burpFileSatId    - input  -  identificateur du satellite
    !   rogueFac         -internal-  constant which determines the severity of the
    !                                quality control applied to the O-P magnitude
    !                                for each channel
    !   productRogueSTD  -internal-  product of rogueFac and standard deviation for
    !                                each channel
    !   oer_tovutil      -external-  UTIL values read from the total error statistics file
    !                                0 = blacklisted channel
    !                                1 = assimilated channel
    !   oer_toverrst     -external-  standard deviation statistics read from the total
    !                                error statistics file
    !------------------------------------------------------------------

    implicit none

    ! Arguments
    integer,          intent(in)    :: obsChannels(:)   ! Channel numbers
    real   ,          intent(in)    :: ompTb(:)         ! Radiance residuals
    integer,          intent(out)   :: flagsInovQc(:)   ! Flags for assimilation/rejection of obs
    integer,          intent(in)    :: actualNumChannel ! Number of channels
    integer,          intent(in)    :: numObsToProcess  ! Number of obs points to process
    integer,          intent(in)    :: sensorIndex      ! Identification number of satellite
    character(len=9), intent(in)    :: burpFileSatId    ! Satellite identification in BURP file
    integer,          intent(in)    :: obsFlags(:)      ! Radiance data flags

    ! Locals
    real, parameter :: factorCh1 = 2.0 ! factor for channel 1 O-P for rejection of channels 1-4
    real, parameter :: maxOmpCh8 = 5.0 ! max Abs(O-P) for channel 8 for rejection of channels 8-11 (units = K)

    integer         :: chanIndex
    integer         :: chanIndex1
    integer         :: chanIndex4
    integer         :: chanIndex8
    integer         :: chanIndex11
    integer         :: chanIndexLast
    integer         :: channelNumber
    integer         :: obsChanIndex
    integer         :: obsIndex

    real            :: productRogueSTD
    real            :: rogueFac(ssbg_maxNumChan)

    !  Initialization. Assume all observations are to be kept.

    ! rogueFac(1:7) should be consistent with bgck.satqc_amsua.f (for AMSUA ch. 3-10)
    ! rogueFac(8:11) should be consistent with bgck.satqc_amsub.f (for AMSUB ch. 2-5)
    ! rogueFac(12:18) should be consistent with ssmi_inovqc.ftn90 (for SSMI ch. 1-7)
    ! rogueFac(19:24) are for upper-atmospheric channels (strato/meso-sphere)
    rogueFac = (/2.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 2.0, &
         &       4.0, 4.0, 4.0, 2.0, 2.0, 2.0, 2.0, 2.0, &
         &       2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0/)

    flagsInovQc(:) = 0

    !  Loop through all observation points in the current record and check
    !  first whether a channel is to be systematically removed and second
    !  whether the O-P value for each channel falls within the acceptable
    !  range defined by rogueFac and oer_toverrst (sensorIndex identifies the satellite).

    !  *** Only do check for data that could be assimilated ***

    !  Also, apply multi-channel rejections for AMSU-like channels using
    !  O-P values from channels 1 and 8.
    !     If Abs(O-P) chan. 1 (AMSU-A 3) > 2*errtot, reject ch. 1-4 (AMSU-A 3-6)
    !     If Abs(O-P) chan. 8 (AMSU-B 2) > 5K, reject ch. 8-11 (AMSU-B 2-5)

    !--------------------------------------------------------------------------
    !  1. Rogue check O-P for all data flagged for assimilation
    !--------------------------------------------------------------------------

    HEADER0: do obsChanIndex = 1, numObsToProcess*actualNumChannel

      channelNumber = obsChannels(obsChanIndex)
      productRogueSTD = rogueFac(channelNumber) * oer_toverrst(channelNumber,sensorIndex)

      if ( .not. btest(obsFlags(obsChanIndex),7) ) then

        if ( oer_tovutil(channelNumber,sensorIndex) == 0 .or. (.not. btest(obsFlags(obsChanIndex),6)) ) then

          ! systematic rejection of this channel
          flagsInovQc(obsChanIndex) = 2
          if ( abs( ompTb(obsChanIndex) ) >= productRogueSTD ) then
            ! rogue check failure as well
            flagsInovQc(obsChanIndex) = 4
          end if

        else

          ! keep channel, now perform rogue check
          if ( abs( ompTb(obsChanIndex) ) >= productRogueSTD ) then
            flagsInovQc(obsChanIndex) = 3
          end if

        end if

      else

        flagsInovQc(obsChanIndex) = 1

      end if

      !  Keep statistics of obs rejected by rogue check.
      if ( (flagsInovQc(obsChanIndex) > 2) .and. (ssbg_debug) ) then
        write(*,*) 'check_stddev: '
        write(*,*) burpFileSatId(2:9),' Rogue check reject: ',   &
             &  ' Obs = ',obsChanIndex,' Channel = ',channelNumber,         &
             &  ' Check value = ',productRogueSTD,             &
             &  ' O-P = ',ompTb(obsChanIndex)
      end if

    end do HEADER0

    !----------------------------------------------------------------------------------
    ! 2. Apply additional criteria for rejections of multiple AMSU-like channels using
    !    Channels 1 and 8 Abs(O-P) (for data that could be assimilated [bit 7 off])
    !----------------------------------------------------------------------------------

    ! Loop over observation points (actualNumChannel = num. channels [24])
    HEADER1: do obsIndex = 1, numObsToProcess

      chanIndexLast = obsIndex*actualNumChannel            ! index of ch.24 obs
      chanIndex1    = chanIndexLast - (actualNumChannel-1) ! index of ch.1 obs
      chanIndex4    = chanIndex1 + 3                       ! index of ch.4 obs
      chanIndex8    = chanIndex4 + 4                       ! index of ch.8 obs
      chanIndex11   = chanIndex8 + 3                       ! index of ch.11 obs

      ! AMSU-A like ch. 1-4
      if ( sum(flagsInovQc(chanIndex1:chanIndex4)) /= 4 ) then
        productRogueSTD = factorCh1 * oer_toverrst(1,sensorIndex)
        if ( abs(ompTb(chanIndex1)) >= productRogueSTD ) then
          do chanIndex = chanIndex1, chanIndex4
            channelNumber = obsChannels(chanIndex)
            if (flagsInovQc(chanIndex) /= 1) flagsInovQc(chanIndex) = max(flagsInovQc(chanIndex),3)
          end do
        end if
      end if

      ! AMSU-B like ch. 8-11
      if ( sum(flagsInovQc(chanIndex8:chanIndex11)) /= 4 ) then
        productRogueSTD = maxOmpCh8
        if ( abs(ompTb(chanIndex8)) >= productRogueSTD ) then
          do chanIndex = chanIndex8, chanIndex11
            channelNumber = obsChannels(chanIndex)
            if (flagsInovQc(chanIndex) /= 1) flagsInovQc(chanIndex) = max(flagsInovQc(chanIndex),3)
          end do
        end if
      end if

    end do HEADER1

  end subroutine check_stddev

  !--------------------------------------------------------------------------
  ! check_topo
  !--------------------------------------------------------------------------
  subroutine check_topo(modelInterpTer, flagsInovQc, actualNumChannel, numObsToProcess)

    ! :Purpose: Perform rejection of observations for selected channels based
    !           on model surface height (for channels assimilated over land)

    !           -- for a single satellite (burpFileSatId,sensorIndex)
    !           -- for a single box (nt observations)

    !------------------------------------------------------------------
    ! Variable Definitions:
    ! ---------------------
    !   modelInterpTer    - input  -  model surface height (m) at each observation point
    !   flagsInovQc       - in/out -  quality control indicator for each channel of each
    !                                 observation point
    !                                 -- on INPUT
    !                                 =0  ok
    !                                 =1  not checked because FLAG bit 7 ON,
    !                                 =2  reject by UTIL value or bit 6 OFF (data not bias corrected),
    !                                 =3  reject by rogue check,
    !                                 =4  reject by both UTIL value and rogue check.
    !                                 -- on OUTPUT,
    !                                 =0  ok
    !                                 =1  not checked because FLAG bit 7 ON,
    !                                 =2  reject by UTIL value,
    !                                 =3  reject by rogue check,
    !                                 =4  reject by both UTIL value and rogue check.
    !                                 =5  rejection due to topography
    !                                 =6  rejection due to topography, reject by UTIL value
    !                                 =7  rejection due to topography, reject by rogue check
    !                                 =8  rejection due to topography, reject by UTIL value,
    !                                     reject by rogue check
    !                                 ==> individual channels may be rejected from each obs pt
    !   actualNumChannel  - input  -  number of residual channels
    !   numObsToProcess   - input  -  number of groups of NVAL*NELE
    !------------------------------------------------------------------

    implicit none

    !  Arguments
    real,    intent(in)    :: modelInterpTer(:)  ! Model surface height (m) for each obs
    integer, intent(inout) :: flagsInovQc(:)     ! Flags for assimilation/rejection of obs
    integer, intent(in)    :: actualNumChannel   ! Number of channels
    integer, intent(in)    :: numObsToProcess    ! Number of obs points to process

    !  Locals
    integer, parameter :: nChanCheck=4             ! number of channels to check

    integer            :: chanIndex
    integer            :: chanIndex1
    integer            :: checkedChan(nChanCheck)
    integer            :: checkedChanIndex
    integer            :: obsIndex
    integer            :: topoReject

    logical            :: debugSupp

    real               :: heightLimit
    real               :: topoHeight
    real               :: topoLimit(nChanCheck)

    !------------------------------------------------------------------
    ! Define channels to check and height limits (m) for rejection
    !------------------------------------------------------------------
    checkedChan = (/     4,     9,    10,    11 /)
    topoLimit   = (/  250., 1000., 2000., 2500. /)

    if (ssbg_debug) then
       write(*,*) 'check_topo: Maximum input GZ = ', maxval(modelInterpTer)
       write(*,*) 'check_topo: numObsToProcess = ', numObsToProcess
    end if

    !------------------------------------------------------------------
    ! Perform the check for the selected channels
    !------------------------------------------------------------------

    !   Loop over all numObsToProcess observation locations

    topoReject = 0
    HEADER: do obsIndex = 1, numObsToProcess
      debugSupp = .false.
      chanIndex1 = (obsIndex-1)*actualNumChannel + 1  ! index of ch.1 obs
      topoHeight = modelInterpTer(obsIndex)           ! model topography height [m]

      if (ssbg_debug) then
        if (topoHeight == maxval(modelInterpTer) .and. topoHeight > minval(topoLimit)) then
          write(*,*) 'check_topo: ****** Max height point! topoHeight = ', topoHeight
          debugSupp = .true.
        end if
      end if

      SUBHEADER: do checkedChanIndex = 1, nChanCheck        ! loop over channels to check (checkedChan)
        heightLimit = topoLimit(checkedChanIndex)   ! height limit [m] for channel checkedChan(checkedChanIndex)
        chanIndex = chanIndex1 + (checkedChan(checkedChanIndex)-1)  ! channel checkedChan(checkedChanIndex) index
        if ( flagsInovQc(chanIndex) /= 1 ) then
          if ( topoHeight > heightLimit ) then
            flagsInovQc(chanIndex) = max(1,flagsInovQc(chanIndex)) + 4
            if (ssbg_debug .and. debugSupp) write(*,*) 'check_topo: Incrementing topoReject for max topoHeight point for ch.= ', checkedChan(checkedChanIndex)
            topoReject = topoReject + 1
            if (ssbg_debug) then
              if ( topoReject <= nChanCheck ) then
                write(*,*) 'check_topo:'
                write(*,*) ' Channel =          ', checkedChan(checkedChanIndex), &
                     &     ' Height limit (m) = ', heightLimit,                   &
                     &     ' Model height (m) = ', topoHeight
              end if
            end if
          end if
        end if
      end do SUBHEADER

    end do HEADER

    if (ssbg_debug .and. (topoReject > 0) ) then
      write(*,*) 'check_topo: Number of topography rejections and observations for this box = ', topoReject, numObsToProcess*actualNumChannel
    end if

  end subroutine check_topo

  !--------------------------------------------------------------------------
  ! ssbg_updateObsSpaceAfterInovQc
  !--------------------------------------------------------------------------
  subroutine ssbg_updateObsSpaceAfterInovQc(obsSpaceData, headerIndex, flagsInovQc)
    ! :Purpose: Update obspacedata variables (obstTB and obs flags) after QC

    implicit none

    !Arguments
    type(struct_obs), intent(inout) :: obsSpaceData    ! ObsSpaceData object
    integer,          intent(in)    :: headerIndex     ! Current header index
    integer,          intent(in)    :: flagsInovQc(:)  ! Flags for assimilation/rejection of obs

    ! Locals
    integer, allocatable :: obsFlags(:)
    integer, allocatable :: obsGlobalFlag(:)
    integer, allocatable :: satScanPosition(:)

    logical              :: sensorIndexFound

    integer              :: actualNumChannel
    integer              :: bodyIndex
    integer              :: bodyIndexBeg
    integer              :: bodyIndexEnd
    integer              :: channelIndex
    integer              :: currentChannelNumber
    integer              :: dataIndex
    integer              :: headerCompt
    integer              :: instr
    integer              :: instrId
    integer              :: iSat
    integer              :: numObsToProcess
    integer              :: obsIndex
    integer              :: platf
    integer              :: platfId
    integer              :: sensorIndex

   ! find tvs_sensor index corresponding to current obs

    platf = obs_headElem_i( obsSpaceData, OBS_SAT, headerIndex )
    instr = obs_headElem_i( obsSpaceData, OBS_INS, headerIndex )

    call tvs_mapSat( platf, platfId, iSat )
    call tvs_mapInstrum( instr, instrId )

    sensorIndexFound = .false.
    HEADER: do sensorIndex = 1, tvs_nsensors
      if ( platfId == tvs_platforms(sensorIndex)  .and. &
           iSat    == tvs_satellites(sensorIndex) .and. &
           instrId == tvs_instruments(sensorIndex) ) then
        sensorIndexFound = .true.
        exit HEADER
      end if
    end do HEADER
    if ( .not. sensorIndexFound ) call utl_abort('ssbg_updateObsSpaceAfterInovQc: sensor Index not found')

    ! find actual Number of channels
    actualNumChannel = tvs_coefs(sensorIndex)%coef%fmv_ori_nchn

    headerCompt = 1
    numObsToProcess = 1
    ! Allocate Header elements
    call utl_reAllocate(obsGlobalFlag, numObsToProcess)
    call utl_reAllocate(satScanPosition, numObsToProcess)

    ! Allocation
    call utl_reAllocate(obsFlags,numObsToProcess*actualNumChannel)

    ! Read elements in obsspace

    obsGlobalFlag(headerCompt)    = obs_headElem_i( obsSpaceData, OBS_ST1, headerIndex )
    satScanPosition(headerCompt)  = obs_headElem_i( obsSpaceData, OBS_FOV, headerIndex )
    bodyIndexBeg                  = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    bodyIndexEnd = bodyIndexBeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) -1

    BODY2: do bodyIndex =  bodyIndexbeg, bodyIndexEnd
      currentChannelNumber = nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      obsFlags(currentChannelNumber) = obs_bodyElem_i( obsSpaceData, OBS_FLG, bodyIndex )
    end do BODY2
    ! Modify flags

    !------------------------------------------------------------------
    ! 1 - Mark global flags if any flagsInovQc > 0
    !------------------------------------------------------------------
    do obsIndex = 1 , numObsToProcess
      if ( maxval(flagsInovQc) > 0 ) obsGlobalFlag(obsIndex) = ibset(obsGlobalFlag(obsIndex), 6)
    end do

    !------------------------------------------------------------------
    ! 2 - Mark obs flags for each value of flagsInovQc
    !------------------------------------------------------------------

    dataIndex = 0
    do obsIndex = 1 , numObsToProcess
      do channelIndex = 1, actualNumChannel
        dataIndex = dataIndex+1
        if (resetQc) obsFlags(dataIndex) = 0

        select case (flagsInovQc(dataIndex))
        case(1)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),9)
        case(2)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),11)
        case(3)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),9)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),16)
        case(4)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),9)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),11)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),16)
        case(5)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),9)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),18)
        case(6)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),9)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),11)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),18)
        case(7)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),9)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),16)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),18)
        case(8)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),9)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),11)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),16)
          obsFlags(dataIndex) = ibset(obsFlags(dataIndex),18)
        end select

      end do
    end do

    !-----------------------------------------------------------------
    !    Subtract 270 from FOV values (element 005043).
    !-----------------------------------------------------------------
    do obsIndex = 1 , numObsToProcess
      if (satScanPosition(obsIndex) > 270) satScanPosition(obsIndex) = satScanPosition(obsIndex) - 270
    end do

    ! write elements in obsspace
    call obs_headSet_i(obsSpaceData, OBS_FOV, headerIndex, satScanPosition(1))
    call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex, obsGlobalFlag(1))

    bodyIndexBeg = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    bodyIndexEnd = bodyIndexBeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) -1
    BODY: do bodyIndex = bodyIndexBeg, bodyIndexEnd
      currentChannelNumber=nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, obsFlags(currentChannelNumber))
    end do BODY

  end subroutine ssbg_updateObsSpaceAfterInovQc

  !--------------------------------------------------------------------------
  ! ssbg_bgCheckSSMIS
  !--------------------------------------------------------------------------
  subroutine ssbg_bgCheckSSMIS(obsSpaceData)
    ! :Purpose: Do the background check for SSMIS data (satQC and inovQC).

    implicit none

    ! Arguments
    type(struct_obs),     intent(inout) :: obsSpaceData           ! ObsSpaceData object

    ! Locals
    integer, allocatable                :: flagsInovQc(:)

    logical, allocatable                :: obsToReject(:)

    integer                             :: codtyp
    integer                             :: dataIndex
    integer                             :: dataIndex1
    integer                             :: headerIndex
    integer                             :: indexFlags
    integer                             :: inovQcSize
    integer                             :: statsInovQcFlags(10)

    logical                             :: otherDataPresent
    logical                             :: ssmisDataPresent

    real                                :: percentInovQcFlags(9)

    write(*,*) 'ssbg_bgCheckSSMIS: Starting'

    call utl_tmg_start(119,'--BgckSSMIS')
    otherDataPresent = .false.
    ssmisDataPresent = .false.
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER0
      codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if ( tvs_isIdBurpInst(codtyp,'ssmis' ) ) then
        ssmisDataPresent = .true.
      else
        otherDataPresent = .true.
      end if
    end do HEADER0

    if ( ssmisDataPresent .and. otherDataPresent ) then
      call utl_abort ('ssbg_bgCheckSSMIS: Other data than SSMIS also included in obsSpaceData')
    endif

    if ( .not. ssmisDataPresent ) then
      write(*,*) 'WARNING: WILL NOT RUN ssbg_bgCheckSSMIS since no SSMIS'
      return
    end if

    statsInovQcFlags(:) = 0
    percentInovQcFlags(:) = 0.0

    ! read nambgck
    call ssbg_init()
    !Quality Control loop over all observations
    !
    ! loop over all header indices of the specified family with surface obs

    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if ( .not. (tvs_isIdBurpInst(codtyp,'ssmis')) ) then
        write(*,*) 'WARNING in ssbg_bgCheckSSMIS: Observation with codtyp = ', codtyp, ' is not SSMIS'
        cycle HEADER
      end if

      !###############################################################################
      ! STEP 1) call satQC SSMIS program                                             !
      !###############################################################################
      if (ssbg_debug) write(*,*) 'ssbg_bgCheckSSMIS: STEP 1) call satQC SSMIS program'
      call ssbg_satqcSsmis(obsSpaceData, headerIndex, obsToReject)

      !###############################################################################
      ! STEP 2) update Flags after satQC SSMIS program                               !
      !###############################################################################
      if (ssbg_debug) write(*,*) 'ssbg_bgCheckSSMIS: STEP 2) update Flags after satQC SSMIS program'
      call ssbg_updateObsSpaceAfterSatQc(obsSpaceData, headerIndex, obsToReject)

      !###############################################################################
      ! STEP 3) call inovQC SSMIS program                                            !
      !###############################################################################
      if (ssbg_debug) write(*,*) 'ssbg_bgCheckSSMIS: STEP 3) call inovQC SSMIS program'
      call ssbg_inovqcSsmis(obsSpaceData, headerIndex, flagsInovQc)

      !###############################################################################
      ! STEP 4) update Flags after inovQC SSMIS program                              !
      !###############################################################################
      if (ssbg_debug) write(*,*) 'ssbg_bgCheckSSMIS: STEP 4) update Flags after inovQC SSMIS program'
      call ssbg_updateObsSpaceAfterInovQc(obsSpaceData, headerIndex, flagsInovQc)

      !###############################################################################
      ! STEP 5) compute statistics of different inovQc flags types                   !
      !###############################################################################
      inovQcSize = size(flagsInovQc)
      if (maxval(flagsInovQc) > 8) call utl_abort('ssbg_bgCheckSSMIS: Problem with flagsInovQc, value greater than 8.')
      do dataIndex = 1,inovQcSize
        dataIndex1 = flagsInovQc(dataIndex)+1
        ! Counting number of flags with value flagsInovQc(dataIndex)
        statsInovQcFlags(dataIndex1) = statsInovQcFlags(dataIndex1) + 1
        ! Counting total number of flags (observations)
        statsInovQcFlags(10) = statsInovQcFlags(10) + 1
      end do

    end do HEADER

    !###############################################################################
    ! STEP 6) displaying statistics of inovQc flags                                !
    !###############################################################################

    do indexFlags = 1,9
      percentInovQcFlags(indexFlags) = float(statsInovQcFlags(indexFlags))/float(statsInovQcFlags(10)-statsInovQcFlags(2))*100
    end do

    write(*,*)   '------------------- Innovation Rejection Statistics -----------------------'
    write(*,*)   '     Flag description                                   Nm. Obs.  Prct (%) '
    write(*,256) ' Total number of observations :                        ', statsInovQcFlags(10)
    write(*,256) ' (1) Not checked because flag bit 7 ON :               ', statsInovQcFlags(2)
    write(*,256) ' Remaining number of observations :                    ', statsInovQcFlags(10)-statsInovQcFlags(2)
    write(*,257) ' (0) Observations that are OK :                        ', statsInovQcFlags(1), percentInovQcFlags(1)
    write(*,257) ' (2) Rejected by UTIL value or bit 6 OFF :             ', statsInovQcFlags(3), percentInovQcFlags(3)
    write(*,257) ' (3) Rejected by rogue check (O-P) :                   ', statsInovQcFlags(4), percentInovQcFlags(4)
    write(*,257) ' (4) Rejected by both UTIL and rogue check :           ', statsInovQcFlags(5), percentInovQcFlags(5)
    write(*,257) ' (5) Topography rejection :                            ', statsInovQcFlags(6), percentInovQcFlags(6)
    write(*,257) ' (6) Topography rejection and by UTIL value :          ', statsInovQcFlags(7), percentInovQcFlags(7)
    write(*,257) ' (7) Topography rejection and by rogue check :         ', statsInovQcFlags(8), percentInovQcFlags(8)
    write(*,257) ' (8) Topography rejection, by UTIL and rogue check :   ', statsInovQcFlags(9), percentInovQcFlags(9)
    write(*,*)   '---------------------------------------------------------------------------'

256 format(A55,i9)
257 format(A55,i9,f7.2,' %')

    call tmg_stop(119)

    write(*,*) 'ssbg_bgCheckSSMIS: Finished'

  end subroutine ssbg_bgCheckSSMIS

end module bgckssmis_mod


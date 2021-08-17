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
  ! :Purpose: Variables for microwave background check and quality control.
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

  character(len=128), parameter :: glmg_file='fstglmg'  ! glace de mer file
  character(len=128), parameter :: fileGlace='bicefil'  ! binaire 0.1degre ice file
  character(len=128), parameter :: wentz_file='wentz_surf.std'  ! surface wentz file
  character(len=128), parameter :: alg_option = 'fwentz'
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
    !
    !:Purpose: This subroutine reads the namelist section NAMBGCK
    !          for the module.
    implicit none

    ! Locals:
    integer :: nulnam, ierr
    integer, external :: fnom, fclos

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

  !--------------------------------------------------------------------------
  ! ISRCHEQI function
  !--------------------------------------------------------------------------
  function ISRCHEQI(KLIST, KLEN, KENTRY) result(ISRCHEQI_out)
    !OBJET          Rechercher un element dans une liste (valeurs entieres).
    !ARGUMENTS      - indx    - output -  position de l'element recherche:
    !                                   =0, element introuvable,
    !                                   >0, position de l'element trouve,
    !               - klist   - input  -  la liste
    !               - klen    - input  -  longueur de la liste
    !               - kentry  - input  -  l'element recherche

    implicit none

    integer :: ISRCHEQI_out
    integer  KLEN, JI

    integer  KLIST(KLEN)
    integer  KENTRY

    ISRCHEQI_out = 0
    do JI=1,KLEN
       if ( KLIST(JI) .EQ. KENTRY ) then
          ISRCHEQI_out = JI
          return
       end if
    end do

  end function ISRCHEQI

  !--------------------------------------------------------------------
  ! ssmis_tb2ta
  !--------------------------------------------------------------------

  subroutine ssmis_tb2ta(npts,grossrej,ztb,zta)
    !--------------------------------------------------------------------
    !  Language: FORTRAN 90
    !
    !  Object: Convert Tbs received from UKMO to Tas, by reversing Ta to Tb
    !          spillover correction applied in B. Bell's pre-processing.
    !
    !  Call:   same as above
    !
    !  History
    !  Version:      Date:      Comment:
    !  --------      -----      --------
    !    0.1       02/03/07     Original code.          D. Anselmo
    !    0.2       28/03/07     Apply correction only to obs that pass
    !                           gross reject filter     D. Anselmo
    !    0.3       09/06/10     Updated spill_coeff     S. Macpherson
    !                           from S. Swadley (NRL)
    !--------------------------------------------------------------------
    ! Variable Definitions:
    ! ---------------------
    ! npts            - input  -  number of obs pts to process
    ! grossrej        - input  -  gross reject indicator
    ! ztb             - input  -  Tbs from input BURP file
    ! zta             - output -  Tas after conversion
    ! spill_coeffs    -internal-  spillover correction coefficients
    !------------------------------------------------------------------
    implicit none

    !  Arguments
    integer, intent(in) :: npts

    real,    intent(in) ::  ztb(:)
    logical, intent(in) ::  grossrej(:)
    real,    intent(out) :: zta(:)

    !  Locals
    integer :: ii, indx1, indx2
    real    :: spill_coeffs(ssbg_maxNumChan)


    !  Define spillover correction coefficients
    !!  Row1           ch1  ch2  ch3  ch4   ch5   ch6
    !!  Row2           ch7  ch8  ch9  ch10  ch11  ch12
    !!  Row3           ch13 ch14 ch15 ch16  ch17  ch18
    !!  Row4           ch19 ch20 ch21 ch22  ch23  ch24

    !  Spillover coeff for all channels (from Steve Swadley/NRL 9 June 2010)

    !  spill_coeffs = (/ 0.9850,  0.9850,  0.9850,  0.9850,  0.9850,  0.9815, &
    !                    0.9815,  0.9949,  0.9934,  0.9934,  0.9934,  0.9680, &
    !                    0.9720,  0.9820,  0.9810,  0.9850,  0.9820,  0.9780, &
    !                    0.9815,  0.9815,  0.9815,  0.9815,  0.9815,  0.9815 /)
    !
    !  Spillover coeff for 7 SSM/I-like channels (ch. 12-18) only
    !
      spill_coeffs = (/ 1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000, &
                        1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  0.9680, &
                        0.9720,  0.9820,  0.9810,  0.9850,  0.9820,  0.9780, &
                        1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000 /)

    !  Apply Tb -> Ta conversion

      indx1 = 1
      do ii = 1, npts

        indx2 = ii*ssbg_maxNumChan
        if ( .not. grossrej(ii) ) then
          zta(indx1:indx2) = spill_coeffs(:) * ztb(indx1:indx2)
        else
          zta(indx1:indx2) = ztb(indx1:indx2)
        end if

        indx1 = indx2 + 1

      end do

    end subroutine ssmis_tb2ta

    !--------------------------------------------------------------------------
    ! f16tdr_remapping 
    !--------------------------------------------------------------------------  
    subroutine f16tdr_remapping(sat_id, SSMIS_Ta, Remapped_SSMI_Ta)
    !--------------------------------------------------------------------
      !
      !  Purpose:
      !    To remap SSMIS imaging channel antenna temperature to SSMI Ta
      !       SSMIS         C_Freq          SSMI
      !      -------   -----------------   -------
      !       Chan12    19.35h              Chan2
      !       Chan13    19.35v              Chan1
      !       Chan14    22.235v             Chan3
      !       Chan15    37.0h               Chan5
      !       Chan16    37.0v               Chan4
      !       Chan17    91.65v -> 85.5v     Chan6
      !       Chan18    91.65h -> 85.5h     Chan7
      !  
      !  Record of Revisions:
      !        Date        Programmer            Description of Change
      !    ============  ==============  ==========================================
      !     2006/07/17     Ninghai Sun    Create original program (for F16).
      !     2007/02/02     Ninghai Sun    Adjust the coefficients
      !     2007/03/02     D. Anselmo     Modify to suit satqc_ssmis.ftn90
      !     2010/07/15     Banghua Yan    Extend to F17 and F18  (Banghua.Yan@noaa.gov, NOAA/NESDIS)
      !                    Ninghai Sun
      !
      ! Reference: 
      !   (1) Yan, B., and F. Weng, 2008: Intercalibration between Special Sensor Microwave Imager and Sounder 
      ! (SSMIS) and Special Sensor Microwave Imager (SSM/I), IEEE Trans. Geosci. Remote Sens, 46, 984-995.
      !   (2) Yan, B., and F. Weng, 2010: A statisical analysis of SSMIS Brightness temperature errors, biases, 
      !  representativeness and Gaussianity from F16 to F18 for Data Assimilation Application (to be submitted)
      !
      !==============================================================================
      IMPLICIT NONE

      INTEGER,PARAMETER :: f16_id = 1, f17_id = 2, f18_id = 3
      REAL, INTENT(IN) :: SSMIS_Ta(ssbg_maxNumChan)
      REAL, INTENT(OUT) :: Remapped_SSMI_Ta(ssbg_maxNumChan)
      INTEGER, INTENT(IN) :: sat_id
      REAL :: tbx(ssbg_maxNumChan)
      INTEGER(2) :: iChan
      REAL(8) :: CP(ssbg_maxNumChan)

      ! NESDIS Intercept/Slope for F16 SSMIS 7 IMG channels to SSMI linear remapping 
                                 ! ch1  ch2  ch3  ch4  ch5  ch6  ch7  ch8  ch9  ch10 ch11 (LAS/ENV)
      REAL(8), PARAMETER :: AP(ssbg_maxNumChan)=(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                     7.44254,7.80472,6.76383,8.55426,7.34409,6.57813,6.45397, &
                                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

      REAL(8), PARAMETER :: BP(ssbg_maxNumChan)=(/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,&
                                     0.969424,0.967519,0.959808,0.954316,0.958955,0.980339,0.978795, &
                                     1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)

      ! NESDIS F17 and F18 Tb biases with respect to F16
      REAL(8), PARAMETER :: CP_F17(ssbg_maxNumChan)=(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                        -0.779,  -1.446,  -1.013,  -0.522,  -0.240,   0.735,   0.521,   &
                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  
      REAL(8), PARAMETER :: CP_F18(ssbg_maxNumChan)=(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                        -0.773,  -0.688,  -1.031,  -0.632,  -0.411,   0.171,   0.928,   &
                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)


      !Initialization
      tbx(1:ssbg_maxNumChan) = SSMIS_Ta(1:ssbg_maxNumChan)

      if ( sat_id == f17_id ) then
        CP = CP_F17
      else if ( sat_id == f18_id ) then
        CP = CP_F18
      else
        CP = 0.0
      endif


      DO iChan=1, ssbg_maxNumChan
        Remapped_SSMI_Ta(iChan) = AP(iChan) + BP(iChan)*(tbx(iChan)+CP(iChan))
      END DO

  end subroutine f16tdr_remapping


  !--------------------------------------------------------------------------
  ! ssmi_ta2tb_fweng
  !--------------------------------------------------------------------------  
  subroutine ssmi_ta2tb_fweng(Ta, Tb)
    !--------------------------------------------------------------------
    !
    !  Purpose:
    !    To convert antenna temperature(Ta) to brightness temperature(Tb).
    !    (1) All channel antenna gain spill-over correction
    !    (2) Imaging channel Cross-polarization correction
    !    (3) Doppler correction ( To be developed )
    !  
    !  Record of Revisions:
    !        Date        Programmer            Description of Change
    !    ============  ==============  ==========================================
    !     2006/07/17     Ninghai Sun    Create original program.
    !     2007/02/02     Ninghai Sun    Adjust the coefficients
    !
    !==============================================================================
    IMPLICIT NONE

    REAL, INTENT(IN) :: Ta(24)
    REAL, INTENT(OUT) :: Tb(24)

    INTEGER(4) :: iChan
    REAL(4) :: CP(24), DP(24)

    REAL(8), PARAMETER :: AP(24)=(/0.9850,0.9850,0.9850,0.9850,0.9850,0.9790,0.9815,&
                                   0.9949,0.9934,0.9934,0.9934, &
                                   0.9690,0.9690,0.9740,0.9860,0.9860,0.9880,0.9880,&
                                   0.9815,0.9815,0.9815,0.9815,0.9815,0.9815/)

    REAL(8), PARAMETER :: BP(24)=(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                   0.0, 0.0, 0.0, 0.0, &
                                   0.00415,0.00473,0.0107,0.02612,0.0217,0.01383,0.01947,&
                                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

    ! All channel antenna gain correction. Note the cross-polarization effects on antenna gain correction.
    DO iChan=1, 24
      CP(iChan) = 1.0/( AP(iChan)*(1.0 - BP(iChan)) )
      DP(iChan) = CP(iChan) * BP(iChan)
      Tb(iChan) = CP(iChan)*Ta(iChan)
    END DO

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
  !--------------------------------------------------------------------
    !  Adapted from: tatotb.f of process_ssmi.ftn90
    !  Original AUTHOR: G. Deblonde ARMA/AES
    !  Code based on decode2.f from F. Wentz.
    !
    !  OBJECT: Convert antenna temperatures to brightness temperatures.
    !
    !   March 23, 1999: coefficients for APC (antenna pattern
    !                   corrections are the same in decode3.f)
    !  History
    !  Version:      Date:      Comment:
    !  --------      -----      --------
    !    0.1       02/03/07     Original adapted code.          D. Anselmo
    !    0.2       02/03/07     Removed DMSP13 along scan bias correction.
    !                                                           D. Anselmo
    !---------------------------------------------------
    implicit none

    !  Arguments
    real, intent(in)  :: Ta(:)
    real, intent(out) :: Tb(:)

    !  Locals
    real :: TB19V,TB19H,TB22V,TB37V,TB37H,TB85V,TB85H
    real :: TA19V,TA19H,TA22V,TA37V,TA37H,TA85V,TA85H

    !  APC coefficients
    real :: AV19V,AH19V,A019V,AH19H,AV19H,A019H,AV22V,A022V
    real :: AV37V,AH37V,A037V,AH37H,AV37H,A037H
    real :: AV85V,AH85V,A085V,AH85H,AV85H,A085H

    !--------------------------

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


  subroutine compute_iwv_101(Tb,iwv)
    !--------------------------------------------------------------------
    !  Original AUTHOR: G. Deblonde ARMA/AES
    !  From: compute_iwv_101.ftn90 of process_ssmi.ftn90
    !
    !  OBJECT: Compute integrated water vapor from SSMI brightness temperatures.
    !
    !     Original version based on fortran routine: alishouse3.f
    !
    !     OBTAIN INTEGRATED WATER VAPOUR FROM BRIGHTNESS TEMPERATURES over OCEAN
    !     CAL/VAL algorithm p7-22 (Hollinger 1991)
    !     ref also: Alishouse et al.(1990) coefficients (f08 data),p.815
    !     ALso, apply Cubic Polynomial Correction to the CAL/VAL water
    !           vapor algorithm (G.W. Petty--published in proceedings ONLY)
    !
    !     INPUT:
    !      Tb - brightness temperature vector (all 24 SSMIS channels)
    !
    !     OUTPUT:
    !      iwv = IWV_alishouse - estimated total precipitable water (units:kg/m**2)
    !                            has cubic polynomial correction
    !
    !     ALSO COMPUTED:
    !      IWV_alishouse0 = estimated total precipitable water (units:kg/m**2)
    !                       does not have cubic polynomial correction
    !      precip_screen = 1: possible presence of precipitation does
    !                         not allow retrieval of IWV and CLW.
    !                      0: retrieval is possible
    !
    !  History
    !  Version:      Date:      Comment:
    !  --------      -----      --------
    !    0.1       06/03/07     Original adapted code.          D. Anselmo
    !-------------------------------------------------------------
    implicit none

    !  Arguments:
    REAL, intent(in) :: Tb(24)
    real, intent(out) :: iwv

    !  Locals:
    integer :: precip_screen

    real :: xx
    real :: tb19v,tb22v,tb37v,tb37h
    real :: IWV_alishouse0,IWV_alishouse
    real :: ciwv(0:4)
    !-----------------------------------------------------------------
    ! write(6,*) ' COMPUTE_IWV_101: IWV calculation    '
    ! write(6,*) ' Technique used = Modified Alishouse '

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

    xx = -11.7939 - 0.02727*tb37v + 0.09920*tb37h
    if ( xx > 0.0 ) then
      precip_screen = 1
    else
      precip_screen = 0
    endif


    !  Apply precipitation screen to IWV_alishouse (units are kg/m**2).
    !  Bounds now applied in cld_filter_fweng.ftn90 routinge.

    if( precip_screen .eq. 1 ) then
       IWV_alishouse = ssbg_rmisg
    endif


    !  Store the IWV value.

    iwv = IWV_alishouse

  end subroutine compute_iwv_101


  !--------------------------------------------------------------------------
  ! DETERM_TPW 
  !--------------------------------------------------------------------------  

  subroutine DETERM_TPW(Tb,Stype,Sea_Ice,TPW)
    !==============================================================================
    !
    !  Purpose:
    !    To calculate total precipitable water .
    !
    !  Paper:
    !    Alishouse et al., 1990, IEEE/new Screening
    !
    !  Output unit:
    !    mm
    !
    !  Record of Revisions:
    !        Date        Programmer            Description of Change
    !    ============  ==============  ==========================================
    !     2006/07/07     Ninghai Sun    Create original program.
    !     2010/07/15     S. Macpherson  Add parameter OCEAN=0
    !
    !==============================================================================
    implicit none

    ! Declare subroutine argurments
    REAL,  INTENT(IN) :: Tb(24)
    REAL, INTENT(IN) :: Sea_Ice
    INTEGER, INTENT(IN) :: Stype
    REAL, INTENT(OUT) :: TPW

    ! Declare local variables
    REAL :: SCT
    REAL :: Tb19V, Tb19H, Tb22V, Tb37V, Tb37H, Tb85V, Tb85H

    INTEGER, PARAMETER :: OCEAN=0

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
    IF ( Stype == OCEAN ) THEN

      SCT = -182.7 + 0.75*Tb19V + 2.543*Tb22V - 0.00543*Tb22V*Tb22V - Tb85V

      IF ( ABS(Sea_Ice) < 70.0 ) THEN
        TPW = 232.89393 - 0.148596*Tb19V - 1.829125*Tb22V + 0.006193*Tb22V**2 - 0.36954*Tb37V

        TPW = -3.753 + 1.507*TPW - 0.01933*TPW**2 + 0.0002191*TPW**3

        IF ( TPW < 0.0 ) TPW = 0.0
        IF ( TPW > 80.0 ) TPW = 80.0
      END IF

      END IF

  end subroutine DETERM_TPW


  !--------------------------------------------------------------------------
  !  determ_sea_ice
  !--------------------------------------------------------------------------  

  subroutine determ_sea_ice(OCEAN, Ta, Stype, Sea_Ice, Latitude)
    !--------------------------------------------------------------------
    !
    !  Purpose:
    !    To calculate sea ice cover.
    !
    !  Paper:
    !    NOAA Testing algorithm
    !
    !  Output unit:
    !    %
    !
    !  Record of Revisions:
    !        Date        Programmer            Description of Change
    !    ============  ==============  ==========================================
    !     2006/07/07     Ninghai Sun    Create original program.
    !     2007/02/02     Ninghai Sun    Adjust input/output terms.
    !     2007/12/04     Ninghai Sun    Add latitude for sea ice threshold
    !
    !==============================================================================
    IMPLICIT NONE

    REAL, INTENT(IN) :: Ta(24)
    INTEGER, INTENT(IN) :: OCEAN,Stype
    REAL, INTENT(OUT) :: Sea_Ice
    REAL, INTENT(IN)  :: Latitude

    REAL :: Ta19v, Ta19h, Ta22v, Ta37v, Ta37h, Ta85v, Ta85h

    ! Extract 7 IMG channels from Ta
    Ta19v = Ta(13)
    Ta19h = Ta(12)
    Ta22v = Ta(14)
    Ta37v = Ta(16)
    Ta37h = Ta(15)
    Ta85v = Ta(17)
    Ta85h = Ta(18)

    ! Calculate Sea Ice Coverage
    IF ( Stype == OCEAN ) THEN   ! Over Ocean
      IF (Latitude > 44.4 .OR. Latitude < -52.0 ) THEN
        Sea_Ice = 91.9 - 2.994*Ta22v + 2.846*Ta19v - 0.386*Ta37v + 0.495*Ta85v &
                  + 1.005*Ta19h - 0.904*Ta37h
        IF ( Sea_Ice >= 70.0 ) THEN
          Sea_Ice = 100.0
        ELSE
          Sea_Ice = 0.0
        END IF
      ELSE
        Sea_Ice = 0.0
      END IF
    ELSE         ! Over Land, there is no sea ice!
      Sea_Ice = ssbg_rmisg
    END IF

  end subroutine determ_sea_ice

  !--------------------------------------------------------------------------
  !  determ_clw
  !--------------------------------------------------------------------------  


  subroutine determ_clw(alg_option, Ta, Tb, Stype, CLW, WVP, Lat)
    !--------------------------------------------------------------------
    !
    !  Purpose:
    !    To calculate cloud liquid water for a single data point.
    !    Normally called when Stype = 0 (open water).
    !    Also retrieves sea-ice to see if there is sea-ice at water point.
    !
    !  Paper:
    !    Weng, F. et al, 1997, Journal of Climate
    !
    !  Output unit:
    !    kg/m^2
    !
    !  Record of Revisions:
    !        Date        Programmer            Description of Change
    !    ============  ==============  ==========================================
    !     2006/07/07     Ninghai Sun    Create original program.
    !     2007/02/02     Ninghai Sun    Adjust input/output terms.
    !     2007/04/14     S. Macpherson  Return CLW = -500 if ice detected
    !     2010/07/15     S. Macpherson  Add Lat argument for determ_sea_ice.
    !                                   Add option to call determ_tpw.
    !     2010/10/29     S. Macpherson  Let program change Stype (for new sea-ice)
    !==============================================================================
    IMPLICIT NONE

    character(len=6), intent(in) :: alg_option
    REAL, INTENT(IN) :: Ta(24), Tb(24)
    INTEGER, INTENT(INOUT) :: Stype
    REAL, INTENT(OUT)   :: CLW
    REAL, INTENT(INOUT) :: WVP
    REAL, INTENT(IN)    :: Lat

    INTEGER(4), PARAMETER :: RT = 285
    INTEGER, PARAMETER :: OCEAN  = 0
    INTEGER, PARAMETER :: SEAICE = 1
    real,    parameter :: clw_limit = 6.0

    REAL :: Ta19v, Ta19h, Ta22v, Ta37v, Ta37h, Ta85v, Ta85h
    REAL :: Tb19v, Tb22v, Tb37v
    REAL :: Sea_Ice
    REAL :: ALG1, ALG2, ALG3
    REAL :: TPW

    CLW = ssbg_rmisg

    ! Extract 7 IMG channels from Ta and Tb
    Ta19v = Ta(13)
    Ta19h = Ta(12)
    Ta22v = Ta(14)
    Ta37v = Ta(16)
    Ta37h = Ta(15)
    Ta85v = Ta(17)
    Ta85h = Ta(18)

    ! Call DETERM_SEA_ICE to find the Sea_Ice
    !  -- Sea_Ice = 100.0  when sea ice >= 70%
    !             = 0.0    when sea ice  < 70%
    CALL determ_sea_ice(OCEAN, Ta, Stype, Sea_Ice, Lat)

    ! Calculate CLW Over Ocean
    IF ( (Stype == OCEAN) .AND. (Sea_Ice /= 100.0) ) THEN
      ALG1 = ssbg_rmisg
      ALG2 = ssbg_rmisg
      ALG3 = ssbg_rmisg

      if ( trim(alg_option) == 'fweng' ) then
        ! Compute IWV using F. Weng algorithm.
        WVP = 232.89 - 0.1486*Tb19v - 0.3695*Tb37v - (1.8291 - 0.006193*Tb22v)*Tb22v
        if ( WVP < 0.0 ) WVP = 0.0
      end if
    if ( trim(alg_option) == 'nsun') then
      CALL determ_tpw(Tb,Stype,Sea_Ice,TPW)
      WVP = TPW
    endif

    IF ( (Ta19v < RT) .AND. (Ta22v < RT) ) THEN
      ALG1 = -3.20 * ( ALOG(290.0-Ta19v) - 2.80 - 0.42*ALOG(290.0-Ta22v) )          !TA
      ! ALG1 = -3.20 * ( ALOG(290.0-Tb19v) - 2.84 - 0.40*ALOG(290.0-Tb22v) )      !TB
    END IF

    IF ( (Ta37v < RT) .AND. (Ta22v < RT) ) THEN
      ALG2 = -1.66 * ( ALOG(290.0-Ta37v) - 2.90 - 0.349*ALOG(290.0-Ta22v) )   !TA
      ! ALG2 = -1.66 * ( ALOG(290.0-Tb37v) - 2.99 - 0.32*ALOG(290.0-Tb22v) )    !TB
    END IF

    IF ( (Ta85h < RT) .AND. (Ta22v < RT) ) THEN
      ALG3 = -0.44 * ( ALOG(290.0-Ta85h) + 1.60 - 1.354*ALOG(290.0-Ta22v) )     !TA
      ! ALG3 = -0.44 * ( ALOG(290.0-Tb85h) + 1.11 - 1.26*ALOG(290.0-Tb22v) )      !TB
    END IF

    IF ( ALG1 > 0.70 ) THEN
      CLW = ALG1
    ELSE IF ( ALG2 > 0.28 ) THEN
      CLW = ALG2
    ELSE IF ( WVP < 30.0 ) THEN
      CLW = ALG3
    ELSE
      CLW = ALG2
    END IF

    ! Verify CLW is within acceptable upper limit.
    ! (Also in compute_clwp_weng.f of process_ssmi.ftn90.)
    IF ( CLW > clw_limit ) CLW = ssbg_rmisg

    ! Force negative CLW values to zero.
    IF ( CLW < 0.0 .and. CLW /= ssbg_rmisg ) CLW = 0.0

    ELSE
      ! Sea Ice (>70%) detected from s/r determ_sea_ice but Stype was 0 = waterobs (on call)
      !    write(6,*) ' DETERM_CLW_SUB: Sea_Ice detected for waterobs point!!!! '
      !    write(6,*) '                 Point will be changed to land/ice point.'
      CLW = -500.0
      WVP = ssbg_rmisg
      Stype = SEAICE

    END IF

  end subroutine determ_clw

  !--------------------------------------------------------------------------
  !  cld_filter_fweng
  !--------------------------------------------------------------------------  


  subroutine cld_filter_fweng(npts,ztb,alg_option,waterobs,grossrej,  &
            &                 cloudobs,iwvreject,precipobs,rclw,riwv,isatid,zlat, &
            &                 iNewIce)
    !--------------------------------------------------------------------
    !
    !  Author of Adapted code:   D. Anselmo  MSC/ARMA
    !  Language: FORTRAN 90
    !
    !  Object:  Compute the cloud liquid water from SSMIS channels using the
    !           regression algorithm of Fuzhong Weng and Ninghai Sun.
    !
    !  Purpose:
    !    Program to retrieve Cloud Liquid Water Path from F16 SSMIS TDR data
    !
    !  Language:
    !    FORTRAN90/95
    !
    !  Internal Subroutines:
    !    F16TDR_REMAPPING    - SSMIS imaging channels (Channel 12 - 18) to SSMI
    !                            Remapping (N. Sun, B. Yan, 2010)
    !    ssmi_ta2tb_fweng    - Antenna temperature (Ta) to brightness temperature (Tb)
    !                            conversion using F. Weng algorithm
    !    ssmi_ta2tb_fwentz   - Antenna temperature (Ta) to brightness temperature (Tb)
    !                            conversion using F. Wentz algorithm
    !    DETERM_CLW          - kg/m**2 (Weng and Grody, 1994, JGR)
    !    DETERM_SEA_ICE      - %  (0-100, NOAA TESTING ALGORITHM)
    !
    !  Question?  Please contact Fuzhong.Weng@noaa and Ninghai.Sun@noaa.gov
    !  Copyright(C) 2007 Fuzhong Weng and Ninghai Sun @ NOAA/NESDIS/STAR/SPB
    
    !  Modifications:
    !    S. Macpherson  13 April 2007
    !       -- added flag for preciptation-affected observations (precipobs)
    !       -- IWV and CLW are now computed for all water observations that passed
    !          gross-error check, regardless of UKMet rain flag
    !       -- IWV is set to missing if iwv > iwv_thresh
    !       -- waterobs() set to FALSE if S/R determ_sea_ice (determ_clw) detects
    !          Sea_ice > 70% from Ta obs
    !
    !    S. Macpherson  15 July 2010
    !       -- added argument isatid (for call to F16TDR_REMAPPING)
    !       -- added argument zlat (for call to DETERM_SEA_ICE)
    !       -- added option to call N. Sun routine DETERM_TPW (in DETERM_CLW)
    !          when alg_option = 'nsun'
    !       -- ssbg_clwThresh now declared in module var_declare
    !
    !    S. Macpherson  29 October 2010
    !       -- changed reporting of newly detected sea-ice at waterobs points;
    !          new input/output argument iNewIce (counter)
    !
    IMPLICIT NONE

    !  Define Arguments
    integer, intent(in)    :: npts, isatid
    real,    intent(in)    :: ztb(:)
    real,    intent(inout) :: rclw(:),riwv(:)
    real,    intent(in)    :: zlat(:)
    integer, intent(inout) :: iNewIce

    character(len=6), intent(in) :: alg_option

    logical, intent(in)    :: grossrej(:)
    logical, intent(inout) :: cloudobs(:), precipobs(:)
    logical, intent(inout) :: iwvreject(:), waterobs(:)

    !  Define Parameters
    real, parameter :: iwv_thresh = 80.0    ! Upper bound for IWV in kg/m**2
    integer, parameter :: OCEAN  = 0
    integer, parameter :: SEAICE = 1

    !  Define local variables
    integer :: ii,indx1,indx2,Stype
    real :: clw,iwv, Lat

    real :: zta(ssbg_mxval*ssbg_maxObsNum)

    REAL :: F16TDR(ssbg_maxNumChan)
    REAL :: RemappedTa(ssbg_maxNumChan)
    REAL :: Tb(ssbg_maxNumChan)

    !---------------------------------------------------

    ! Convert Tbs received from UKMO to Tas, by reversing Ta to Tb
    ! spillover correction applied in B. Bell's pre-processing.
    ! Missing Tbs are also checked for here.
    ! Apply to all obs of current record.

    call ssmis_tb2ta(npts,grossrej,ztb,zta)

    rclw(:) = ssbg_rmisg
    riwv(:) = ssbg_rmisg

    indx1 = 1
    do ii = 1, npts

      indx2 = ii*ssbg_maxNumChan
      F16TDR(:) = zta(indx1:indx2)
      Lat = zlat(ii)
      indx1 = indx2 + 1

      ! Obtain CLW and IWV for obs points over open ocean, and where
      ! Tb values have passed gross filter check.
      ! Initialize IWV and CLW to missing for all other cases.

      ! *** S. Macpherson  13 April 2007: Removed UKMO rain flag check in order
      !     to get more CLW and IWV observations over water.
      !     Rain-flagged data are all rejected anyway.

      clw = ssbg_rmisg
      iwv = ssbg_rmisg

      if ( waterobs(ii) .and. (.not. grossrej(ii)) ) then

        Stype = OCEAN     ! Surface-type=Ocean

        ! Call SSMIS TDR to SSMI TDR remapping subroutine

        CALL f16tdr_remapping(isatid, F16TDR, RemappedTa)

        ! Call SSM/I Ta to Tb conversion subroutine

        if ( trim(alg_option) == 'fweng' ) then
          call ssmi_ta2tb_fweng(RemappedTa, Tb)
          ! IWV computed in determ_clw subroutine below.
        else if ( trim(alg_option) == 'fwentz' ) then
          call ssmi_ta2tb_fwentz(RemappedTa, Tb)
          ! Compute IWV using Alishouse and Petty (from process_ssmi),
          ! because it won't be computed in determ_clw.
          ! Missing value for IWV means possible precipitation
          ! and so CLW will not be computed
          call compute_iwv_101(Tb,iwv)
          if ( iwv == ssbg_rmisg ) precipobs(ii) = .true.
        else if ( trim(alg_option) == 'nsun' ) then
          call ssmi_ta2tb_fweng(RemappedTa, Tb)
          ! IWV computed in determ_clw subroutine below.
        else
          write(6,*) ' CLD_FILTER_FWENG: Invalid algorithm option !! '
          call abort()
        end if
        !     write(6,*) 'Tb= ', Tb

        ! Call CLW retrieval algorithm subroutine.
        !  -- also computes and returns (output) IWV if alg_option=fweng or nsun
        !  -- Stype is also changed to SEAICE value if sea-ice is detected

        if ( trim(alg_option) /= 'fwentz' ) then
          CALL determ_clw(alg_option, RemappedTa, Tb, Stype, clw, iwv, Lat)
        else
          if ( .not. precipobs(ii) ) then
            CALL determ_clw(alg_option, RemappedTa, Tb, Stype, clw, iwv, Lat)
          endif
        endif

        ! Check for newly detected sea-ice (Stype changed)

        if ( Stype == SEAICE ) iNewIce = iNewIce + 1

        ! Reject obs with precipitation or cloud amount more than threshold.
        ! Also, reject obs if CLW is missing (undetermined)
        ! Set waterobs flag to false if deter_clw returns "Sea_Ice" value (-500)
        if ( (clw > ssbg_clwThresh) .or. precipobs(ii) ) cloudobs(ii) = .true.
        if ( clw == ssbg_rmisg )  cloudobs(ii) = .true.
        if ( clw == -500.0 ) then
           waterobs(ii) = .false.
           clw = ssbg_rmisg
        endif
        ! Reject obs with IWV value more than threshold and set IWV to missing.
        ! First, set IWV values below zero to 0.0.
        if ( iwv < 0.0 .and. iwv /= ssbg_rmisg ) iwv = 0.0
        if ( iwv > iwv_thresh ) then
          if ( .not. cloudobs(ii) ) iwvreject(ii) = .true.
          iwv = ssbg_rmisg
        endif

        ! Store CLW and IWV
        if (ssbg_debug) then
          write(*,*)'CLOUD BY DETERM_CLW = ', clw
          write(*,*)'IWV BY DETERM_CLW = ', iwv
        endif
        rclw(ii) = clw
        riwv(ii) = iwv
        !   write(6,130) ii,iwv,iwv_thresh,iwvreject(ii)

      end if

    enddo


    ! Replace values of ssbg_rmisg by ssbg_zmisg for rclw and riwv such that proper
    ! missing values are written to output BURP file.

    !where ( rclw == ssbg_rmisg ) rclw = ssbg_zmisg
    !where ( riwv == ssbg_rmisg ) riwv = ssbg_zmisg


    120  format(' ii ',2x,i3,2x,' clw ',f4.2,2x,' CLW threshold = ',f4.2,2x,' cloudobs ',l2)
    130  format(' ii ',2x,i3,2x,' iwv ',f4.2,2x,' IWV threshold = ',f4.2,2x,' iwvreject ',l2)

  end subroutine cld_filter_fweng


  !--------------------------------------------------------------------------
  ! extractParamForGrodyRun
  !--------------------------------------------------------------------------  

  !--------------------------------------------------------------------------
  ! copy1Dimto2DimRealArray
  !--------------------------------------------------------------------------
  subroutine copy1Dimto2DimRealArray(oneDimArray, firstDim, secondDim, twoDimArray)
    !:Purpose: copy 1 dim Real Array into 2D real array given dim1 and dim2 
    implicit none
    ! Arguments
    integer, intent(in)                 :: firstDim
    integer, intent(in)                 :: secondDim
    real,    intent(in)                 :: oneDimArray(firstDim*secondDim)
    real,    intent(inout)              :: twoDimArray(firstDim, secondDim)
    
    !locals
    integer                             :: firstDimIndex 
    integer                             :: secondDimIndex 
    integer                             :: productDimIndex 

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
    !:Purpose: copy 1 dim Integer Array into 2D Integer array given dim1 and dim2 
    implicit none
    ! Arguments
    integer, intent(in)                 :: firstDim
    integer, intent(in)                 :: secondDim
    integer, intent(in)                 :: oneDimArray(firstDim*secondDim)
    integer, intent(inout)              :: twoDimArray(firstDim, secondDim)
    
    !locals
    integer                             :: firstDimIndex 
    integer                             :: secondDimIndex 
    integer                             :: productDimIndex 

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    do secondDimIndex=1,secondDim
      do firstDimIndex=1,firstDim
        productDimIndex = (secondDimIndex-1)*firstDim + firstDimIndex 
        twoDimArray(firstDimIndex,secondDimIndex) = oneDimArray(productDimIndex)
      end do
    end do

  end subroutine copy1Dimto2DimIntegerArray

  !--------------------------------------------------------------------------
  ! resetQcC
  !--------------------------------------------------------------------------
  subroutine resetQcCases(RESETQC, KCHKPRF, globMarq)
    !:Purpose:        allumer la bit (6) indiquant que l'observation a un element
    !                 rejete par le controle de qualite de l'AO.
    !                 N.B.: si on est en mode resetqc, on remet le marqueur global a
    !                 sa valeur de defaut, soit 1024,  avant de faire la mise a jour.
    implicit none
    !Arguments:
    logical,              intent(in)     :: RESETQC       !reset the quality control flags before adding the new ones ?
    integer,              intent(in)     :: KCHKPRF(:)    !indicateur global controle de qualite tovs. Code:
    integer,              intent(inout)  :: globMarq(:)   !Marqueurs globaux  
    !Locals
    
    integer                              :: dataNum 
    integer                              :: dataIndex

    dataNum = size(globMarq)
    do dataIndex = 1, dataNum
      if (RESETQC) then
        globMarq(dataIndex) = 1024  
      end if
      if ( KCHKPRF(dataIndex) /= 0  ) then
        globMarq(dataIndex) = OR (globMarq(dataIndex),2**6)
      end if
    end do
    if (ssbg_debug) then
      write(*,*) ' KCHKPRF   = ', (KCHKPRF(dataIndex),dataIndex=1,dataNum)
      write(*,*) ' NEW FLAGS = ', (globMarq(dataIndex),dataIndex=1,dataNum)
    end if

  end  subroutine resetQcCases

  !------------------------------------------------------------------------------------
  ! bennartz
  !------------------------------------------------------------------------------------
  subroutine bennartz (ier, knt, tb89, tb150, pangl, ktermer, scatl, scatw)

    !:Purpose: Compute the following parameters using 2 AMSU-B channels:
    !          - scattering index (over land and ocean).*
    !          The two channels used are: 89Ghz, 150Ghz.
    !          REGERENCES     Bennartz, R., A. Thoss, A. Dybbroe and D. B. Michelson, 
    !                         1999: Precipitation Analysis from AMSU, Nowcasting SAF, 
    !                         Swedish Meteorologicali and Hydrological Institute, 
    !                         Visiting Scientist Report, November 1999.
    implicit none
    ! arguments: 
    integer, intent(out) :: ier(knt)              ! error return code:
    integer, intent(in)  :: knt              ! number of points to process
    real,    intent(in)  :: tb89(:)          ! 89Ghz AMSU-B brightness temperature (K)
    real,    intent(in)  :: tb150(:)         ! 150Ghz AMSU-B brightness temperature (K)
    real,    intent(in)  :: pangl(:)         !  satellite zenith angle (deg.)
    integer, intent(in)  :: ktermer(:)       ! land/sea indicator (0=land;1=ocean)
    real,    intent(out) :: scatl(:)         ! scattering index over land
    real,    intent(out) :: scatw(:)         ! scattering index over water

    !     Notes: In the case where an output parameter cannot be calculated, the
    !     value of this parameter is to to the missing value, i.e. -99.
    ! Locals: 
    real, parameter :: zmisg = -99.
    integer :: i 

    ! ____1) Initialise output parameters:**    -------------------------------*

    do i = 1, knt 
      scatl(i) = zmisg
      scatw(i) = zmisg
    end do

    !____2) Validate input parameters:**    -----------------------------*
    do i = 1, knt
      if ( tb89(i)      < 120.  .or.     &
            tb89(i)     > 350.  .or.     &
            tb150(i)    < 120.  .or.     &
            tb150(i)    > 350.  .or.     & 
            pangl(i)    < -90.  .or.     &
            pangl(i)    >  90.  .or.     & 
            ktermer(i)  <   0   .or.     &
            ktermer(i)  >   1        ) then
          ier(i) = 1        
      else
          ier(i) = 0      
      end if 
    enddo

    !____3) Compute parameters:**    ----------------------*
    do i = 1, knt 
      if ( ier(i) == 0 ) then
        if (ktermer(i) == 1 ) then
          scatw(i) = (tb89(i)-tb150(i)) -      &
                     (-39.2010+0.1104*pangl(i))
        else
          scatl(i) = (tb89(i)-tb150(i)) -     &
                     (0.158+0.0163*pangl(i))
        endif
      else if ( (ier(i) /= 0  ) .and. (i <= 100 ) .and. (ssbg_debug)) then
        print *, ' Input Parameters are not all valid: '
        print *, ' i,tb89(i),tb150(i),pangl(i),ktermer(i) = ',     &
                   i,tb89(i),tb150(i),pangl(i),ktermer(i)
        print *, ' ier(i),scatl(i),scatw(i)=',     &
                   ier(i),scatl(i),scatw(i)
      endif
    end do 

end subroutine bennartz

  !--------------------------------------------------------------------------
  ! ssbg_readGeophysicFieldsAndInterpolate
  !--------------------------------------------------------------------------
  subroutine ssbg_readGeophysicFieldsAndInterpolate(zlat, zlon, MTINTRP)

    implicit none

    !:Purpose: Reads Modele Geophysical variables and save for the first time
    !          TOPOGRAPHIE (ME, MX ou GZ):
    !             ME est la topographie filtree avec unites en metres (filtered ME).
    !             MX est la topographie filtree avec unites en m2/s2  (geopotential topography).
    !             GZ is geopotential height; GZ at surface = surface height (dam)
    !               -- looks for surface GZ (hybrid or eta level = 1.0)
    !               -- uses new librmn_008 function ip1_all(level,kind)
    !         Then Interpolate Those variables to observation location
    !Arguments: 
    real,               intent(in)   :: zlat(:)        ! Observation Lats
    real,               intent(in)   :: zlon(:)        ! Observation Lons
    real, allocatable,  intent(out)  :: MTINTRP(:)     ! topographie filtree (en metres) et interpolees
  
    ! Locals:
    real, allocatable, save  :: MT(:)                  ! Modele Topographie (MT)
    logical,           save  :: ifFirstCall = .True.   ! If .True. we read GL, MT and MG
    integer,           save  :: gdmt                   ! topo interpolation param
    integer                  :: gdllsval
    integer                  :: ier, irec
    integer                  :: ezqkdef, ezsetopt
    integer, external        :: FSTINF,FSTPRM,FCLOS
    integer, external        :: FSTLIR,FSTFRM, FNOM, FSTOUV
    integer                  :: NI, NJ, NK, IG1, IG2, IG3, IG4
    integer                  :: IDUM1,IDUM2,IDUM3,IDUM4
    integer                  :: IDUM5,IDUM6,IDUM7,IDUM8
    integer                  :: IDUM9,IDUM10,IDUM11,IDUM12,IDUM13
    integer                  :: IDUM14,IDUM15,IDUM16,IDUM17,IDUM18
    character(len=12)        :: ETIKXX
    character(len=4)         :: CLNOMVAR
    character(len=4)         :: NOMVXX
    character(len=2)         :: TYPXX 
    character(len=1)         :: GRTYP
    integer                  :: NLAT
    integer                  :: NLON
    integer                  :: IUNGEO
    integer, PARAMETER       :: MXLON = 5
    integer, PARAMETER       :: MXLAT = 5
    integer, PARAMETER       :: MXELM = 20
    real,    PARAMETER       :: DLAT = 0.4
    real,    PARAMETER       :: DLON = 0.6
    real                     :: XLAT
    real                     :: XLON
    real                     :: TOPOFACT               ! Facteur x topo pour avoir des unites en metre
    real, allocatable        :: ZLATBOX (:,:)
    real, allocatable        :: ZLONBOX (:,:)
    real, allocatable        :: MTINTBOX(:,:)
    integer                  :: dataIndex
    integer                  :: boxPointIndex
    integer                  :: latIndex
    integer                  :: lonIndex
    integer                  :: zlatNum
    integer                  :: zlonNum
    integer                  :: dataNum
    integer                  :: boxPointNum

    !  External function
    integer, external :: ip1_all

    ! STEP 1: CHECK if ZLAT AND ZLON ARE SAME DIMENSION
    zlatNum = size(zlat)
    zlonNum = size(zlon)
    if (zlatNum .ne. zlonNum) then
      call utl_abort ('ssbg_readGeophysicFieldsAndInterpolate: OBSERVATION ZLAT and ZLON should have SAME LENGTH')
    else 
      dataNum = zlatNum
    end if

    ! STEP 2: READ MT from the FST FILE
    if(ifFirstCall) then
      IUNGEO = 0
      IER = FNOM(IUNGEO,'trlm_01','STD+RND+R/O',0)
      IER = FSTOUV(IUNGEO,'RND')

      ! Using hybrid coordinates
      irec = fstinf(iungeo,ni,nj,nk,-1,' ',ip1_all(1.0,5),-1,-1,' ','GZ')
      if (irec < 0) then
        call utl_abort('ssbg_readGeophysicFieldsAndInterpolate: LA TOPOGRAPHIE EST INEXISTANTE')
      endif
      clnomvar = 'GZ'
      topofact = 10.0  ! dam --> m

      if (allocated(MT)) deallocate(MT)
      allocate ( MT(ni*nj), STAT=ier)
      if ( ier /= 0 ) then
        call utl_abort('ssbg_readGeophysicFieldsAndInterpolate: Allocation of array mt failed')
      end if
      ier = fstlir(MT,iungeo,ni,nj,nk,-1,' ',ip1_all(1.0,5),-1,-1,' ','GZ')

      MT(:) = MT(:)*TOPOFACT

      IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
                     IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10,  &
                     IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1, &
                     IG2, IG3, IG4, IDUM12, IDUM13, IDUM14,  &
                     IDUM15, IDUM16, IDUM17, IDUM18 )
      write (*,*) ' GRILLE MT : ',grtyp,ni,nj, &
                     ig1,ig2,ig3,ig4
      ier  = ezsetopt('INTERP_DEGREE','LINEAR')
      gdmt = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)

      IER = FSTFRM(IUNGEO)
      IER = FCLOS(IUNGEO)
      ifFirstCall = .False.
    end if

    ! STEP 3:  Interpolation de la glace et le champ terre/mer du modele aux pts TOVS.
    ! N.B.: on examine ces champs sur une boite centree sur chaque obs.
    boxPointNum = MXLAT*MXLON
    if(allocated(ZLATBOX)) deallocate(ZLATBOX)
    allocate (ZLATBOX(boxPointNum, dataNum) , STAT=ier) 
    if(allocated(ZLONBOX)) deallocate(ZLONBOX)
    allocate (ZLONBOX(boxPointNum, dataNum) , STAT=ier) 
    if(allocated(MTINTBOX)) deallocate(MTINTBOX)
    allocate (MTINTBOX(boxPointNum, dataNum) , STAT=ier) 

    NLAT = (MXLAT-1)/2
    NLON = (MXLON-1)/2
    do dataIndex = 1, dataNum
      boxPointIndex = 0
      do latIndex = -NLAT, NLAT
        XLAT = ZLAT(dataIndex) +latIndex*DLAT
        XLAT = MAX(-90.0,MIN(90.0,XLAT))
        do lonIndex = -NLON, NLON
          boxPointIndex = boxPointIndex + 1
          XLON = ZLON(dataIndex) +lonIndex*DLON
          if ( XLON < -180. ) XLON = XLON + 360.
          if ( XLON >  180. ) XLON = XLON - 360.
          if ( XLON <    0. ) XLON = XLON + 360.
          ZLATBOX(boxPointIndex,dataIndex) = XLAT
          ZLONBOX(boxPointIndex,dataIndex) = XLON
         end do
      end do
    end do

    ier = gdllsval(gdmt,mtintbox,MT,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
    if (ier < 0) then
      call utl_abort ('ssbg_readGeophysicFieldsAndInterpolate: ERROR in the interpolation of MT')
    end if

    if(allocated(MTINTRP)) deallocate(MTINTRP)
    allocate (MTINTRP(dataNum) , STAT=ier) 

    MTINTRP(:) = 0.0
    do dataIndex = 1, dataNum
      if (ssbg_debug) then
        print *, ' ------------------  '
        print *, ' dataIndex = ', dataIndex
        print *, '   '
        print *, ' zlat,zlon = ', zlat(dataIndex), zlon(dataIndex)
        print *, '   '
        print *, ' ZLATBOX = '
        print *,  (ZLATBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        print *, ' ZLONBOX = '
        print *,  (ZLONBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        print *, ' MTINTBOX = '
        print *,  (MTINTBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
      end if
      do boxPointIndex=1, boxPointNum
        MTINTRP(dataIndex) = MAX(MTINTRP(dataIndex),MTINTBOX(boxPointIndex,dataIndex))
      end do
      if (ssbg_debug) then
        print *, ' MTINTRP = ', MTINTRP(dataIndex)
      end if
    end do

  end subroutine ssbg_readGeophysicFieldsAndInterpolate

  !--------------------------------------------------------------------
  !  land_ice_mask_ssmis
  !--------------------------------------------------------------------
  subroutine land_ice_mask_ssmis(npts,zlat,zlon,zlq,ztt,waterobs)
  !--------------------------------------------------------------------
  !
  ! Author:   D. Anselmo  MSC/ARMA
  ! Language: FORTRAN 90
  ! Adapted from: land_ice_mask.ftn90 of process_ssmi
  !
  ! Object:   Determine for each observation point the ice mask value from
  !           the binary file copied to the local work directory.
  !           This may be a user-specified file or it is copied from
  !              /users/dor/afsi/sio/datafiles/data2/ade.maskice10
  !           The ice mask is updated every day at 00 UTC. The binary file has
  !           a resolution of 0.1 deg. Observations with an ice mask value of 0
  !           (=ice; land or sea=-1) are removed. The ice mask value also
  !           determines the terrain-type qualifier (element 13039) for each obs
  !           pt, which is required when writing output to boxed format.
  !
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
  !
  !           In the application of this check, a 5x5 mesh, with spacing defined by rlat_km and
  !           rlon_km, is positioned with its center over an obs pt (2 grid pts on either side
  !           of the obs pt; size of mesh is equal to 4*rlat_km x 4*rlon_km). The values of MG
  !           and LG are evaluated at the grid points of this mesh. The maximum value of each
  !           determines whether the obs pt is too close to ice or land to be retained.
  !           **NOTE: the threshold value for MG has a very strong effect on the distance
  !                   from land that is permitted for an obs to be retained
  !
  !
  !      Maximum FOV             x---x---x---x---x     ^
  !         = 75km x 75km        |   |   |   |   |     |
  !         for Meso-sphere CHs  x---x---x---x---x     |
  !         = 74km x 47km        |   |   |   |   |     |
  !         for 19 GHz           x---x---o---x---x     | = 4*rlat_km
  !                              |   |   |   |   |     | = 4*40 km
  !                           ^  x---x---x---x---x     | = 160 km = 80 km north & south
  !                   rlat_km |  |   |   |   |   |     |
  !                           v  x---x---x---x---x     v
  !                                          <--->
  !                                         rlon_km
  !
  !                              <--------------->
  !                                 = 4*rlon_km
  !                                 = 4*40 km
  !                                 = 160 km = 80 km east & west
  !
  !
  !               MG value = 1.0  ==>  LAND       MG value = 0.0  ==>  OCEAN
  !               LG value = 1.0  ==>  ICE        LG value = 0.0  ==>  NO ICE
  !
  !
  ! Version:      Date:      Comment:
  ! --------      -----      --------
  !   0.1       04/01/07     Original adapted code.        D. Anselmo
  !   0.2       30/08/08     Added input variable zlq      S. Macpherson
  !                          and code to ensure zlq=1 (sea)
  !                          for waterobs=true points.
  !   0.3       19/11/10     Add terrain-type=1 (snow)     S. Macpherson
  !
  !--------------------------------------------------------------------
  !  Variable Definitions
  !  --------------------
  ! glmg_file  - input  -  name of file holding model MG and LG fields
  ! npts       - input  -  number of input obs pts in record
  ! zlat       - input  -  array holding lat values for all obs pts in record
  ! zlon       - input  -  array holding lon values for all obs pts in record
  ! zlq        - in/out -  array holding land/sea qualifier values for all obs
  !                        pts of record (0 = land, 1 = sea)
  ! ztt        - output -  array holding terrain-type values for all obs pts
  !                        of current record
  ! waterobs   - output -  logical array identifying for each obs in current record
  !                        whether it is over open water, far from coast/ice
  ! imask      -internal-  value determined by interpolating obs pt to
  !                        binary ice mask field
  !                            for land/sea: imask = -1
  !                            for ice:      imask = 0
  ! mxlat      -internal-  number of grid pts in lat. direction for mesh
  ! mxlon      -internal-  number of grid pts in lon. direction for mesh
  ! rlat_km    -internal-  spacing desired between mesh grid points in km
  !                        along lat. direction
  ! rlon_km    -internal-  spacing desired between mesh grid points in km
  !                        along lon. direction
  ! dlat       -internal-  spacing between mesh grid points along lon. direction
  !                        in degrees computed from rlat_km
  ! dlon       -internal-  spacing between mesh grid points along lon. direction
  !                        in degrees computed from rlon_km
  ! rkm_per_deg -internal- distance in km per degree
  !                           = Earth radius * PI/180.0
  !                           = 6371.01 km * PI/180.0
  !                           = 111.195 km
  ! nlat,nlon  -internal-  used to define the lat/lon of the grid pts of mesh
  ! zlatbox    -internal-  lat values at all grid pts of mesh for all obs pts
  ! zlonbox    -internal-  lon values at all grid pts of mesh for all obs pts
  ! latmesh    -internal-  lat values at all grid pts of mesh for 1 obs pt
  ! lonmesh    -internal-  lon values at all grid pts of mesh for 1 obs pt
  ! mgintob    -internal-  interpolated MG values at all grid pts of mesh for 1 obs pt
  ! lgintob    -internal-  interpolated LG values at all grid pts of mesh for 1 obs pt
  ! mgintrp    -internal-  max. interpolated MG value on mesh for all obs pts
  ! lgintrp    -internal-  max. interpolated LG value on mesh for all obs pts
  ! glace      -internal-  integer function that defines the terrain-type an obs
  !                        pt by interpolating its position to a 0.1 deg operational
  !                        ice analysis
  !                            --> function written in C language (glace.c)
  ! MGthresh   -internal-  maximum allowable land fraction for obs to be kept
  ! LGthresh   -internal-  maximum allowable ice  fraction for obs to be kept
  !--------------------------------------------------------------------
  implicit none

  !  Arguments:
  
    integer, intent(in)    :: npts
    real,    intent(in)    :: zlat(:), zlon(:)
    integer, intent(inout) :: zlq(:)
    integer, allocatable, intent(out) :: ztt(:)
  
    logical, allocatable, intent(out) :: waterobs(:)

  !  Locals:
    integer, parameter :: mxlat=5,mxlon=5
    integer            :: iungeo
    logical, save      :: firstCall=.true.

    integer :: ier,key
    integer :: ni,nj,nk,nilg,njlg
    integer :: ig1,ig2,ig3,ig4,ig1lg,ig2lg,ig3lg,ig4lg
    integer :: idum1,idum2,idum3
    integer :: idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11
    integer :: idum12,idum13,idum14,idum15,idum16,idum17,idum18
  
    integer :: indx,ii,jj,kk
    integer :: nlat,nlon,imask
  
    !
    real, parameter :: pi=3.141592654
    real, parameter :: MGthresh=0.01,LGthresh=0.01
    real, parameter :: rlat_km=40.0,rlon_km=40.0
    real, parameter :: rkm_per_deg=111.195
  
    real :: xlat,xlatrad,xlon,rii,rjj
    real :: dlat,dlon
    !
    character(len=12) :: etikxx
    character(len=4)  :: nomvxx
    character(len=2)  :: typxx
    character(len=1)  :: grtyp,grtyplg
    !
    ! F90 allocatable arrays:
    real, allocatable, save :: mg(:),lg(:)
    real, allocatable       :: latmesh(:),lonmesh(:)
    real, allocatable       :: mgintob(:),lgintob(:)
    real, allocatable       :: zlatbox(:,:),zlonbox(:,:)
    real, allocatable       :: mgintrp(:),lgintrp(:)

    !
    ! RPNLIB interpolating functions:
    integer :: ezsetopt,ezqkdef
    integer :: gdllsval
    integer, save :: gdid
    integer, save :: gdidlg
    ! --------------------------------------------------------------------

    ! Define FORTRAN FST functions:
    integer, external :: fstinf,fstprm,fstlir
    integer, external :: fstouv,fstfrm,fnom,fclos

    ! Allocate space for arrays holding values on mesh grid pts.
    call utl_reAllocate(latmesh, mxlat*mxlon)
    call utl_reAllocate(lonmesh, mxlat*mxlon)
    call utl_reAllocate(mgintob, mxlat*mxlon)
    call utl_reAllocate(lgintob, mxlat*mxlon)
    call utl_reAllocate(zlatbox, mxlat*mxlon, npts)
    call utl_reAllocate(zlonbox, mxlat*mxlon, npts)
    call utl_reAllocate(ztt, npts)
    call utl_reAllocate(waterobs, npts)

    if (firstCall) then

      firstCall = .false.

      ! Open FST file.
      iungeo = 0
      ier = fnom( iungeo,glmg_file,'STD+RND+R/O',0 )
      ier = fstouv( iungeo,'RND' )

      ! Read MG field.
      key = fstinf(iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ' ,'MG')
      if ( key <  0 ) then
        call utl_abort('ssbg_landIceMaskAtms: The MG field is MISSING')
      end if

      call utl_reAllocate(mg, ni*nj)

      ier = fstlir(mg,iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ','MG')

      ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,    &
                   idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyp,ig1,ig2,  &
                   ig3,ig4,idum12,idum13,idum14,idum15,idum16,idum17,      &
                   idum18)

      ! Read LG field.

      key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'LG')
      if ( key <  0 ) then
        call utl_abort(' LAND_ICE_MASK_SSMIS: The LG field is MISSING ')
      end if
      call utl_reAllocate(lg, nilg*njlg)
      ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','LG')

      ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,          &
          &        idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyplg,ig1lg,ig2lg,  &
          &        ig3lg,ig4lg,idum12,idum13,idum14,idum15,idum16,idum17,        &
          &        idum18)

      gdid = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
      gdidlg = ezqkdef(nilg,njlg,grtyplg,ig1lg,ig2lg,ig3lg,ig4lg,iungeo)

      ier = fstfrm(iungeo)
      ier = fclos(iungeo)

    end if ! firstCall


    ! For each obs pt, define a grid of artificial pts surrounding it.

    nlat = ( mxlat - 1 ) / 2
    nlon = ( mxlon - 1 ) / 2

    dlat = rlat_km / rkm_per_deg
    do kk = 1, npts
      indx = 0

      do ii = -nlat, nlat
        rii = float(ii)
        xlat = zlat(kk) + rii*dlat
        xlat = max( -90.0, min(90.0,xlat) )
        xlatrad = xlat*pi/180.0

        do jj = -nlon, nlon
          dlon = rlon_km / ( rkm_per_deg*cos(xlatrad) )
          rjj = float(jj)
          indx = indx + 1
          xlon = zlon(kk) + rjj*dlon
          if ( xlon < -180. ) xlon = xlon + 360.
          if ( xlon >  180. ) xlon = xlon - 360.
          if ( xlon <    0. ) xlon = xlon + 360.
          zlatbox(indx,kk) = xlat
          zlonbox(indx,kk) = xlon
        end do

      end do
    end do


    ! Interpolate values from MG and LG field to grid pts of mesh centred over each obs pt.
    ! Determine for each obs pt, the max interpolated MG and LG value within the box
    ! surrounding it.

    ier  = ezsetopt('INTERP_DEGREE','LINEAR')
 
    call utl_reAllocate(mgintrp, npts)
    call utl_reAllocate(lgintrp, npts)

    mgintrp(:) = 0.0
    lgintrp(:) = 0.0
    do kk = 1, npts
 
      latmesh = zlatbox(:,kk)
      lonmesh = zlonbox(:,kk)

      ier  = gdllsval(gdid,mgintob,mg,latmesh,lonmesh,mxlat*mxlon)
      ier  = gdllsval(gdidlg,lgintob,lg,latmesh,lonmesh,mxlat*mxlon)

      mgintrp(kk) = maxval(mgintob(:))
      lgintrp(kk) = maxval(lgintob(:))

    end do

    !  Initialize all obs as being over land and free of ice or snow.
    !  Determine which obs are over open water.

    waterobs(:) = .false.   ! not over open water
    ztt(:) = -1             ! no ice or snow
    do kk = 1, npts

      !    Determine for each obs pt, the value of imask from the 0.1 deg binary file.
      !      - open land/sea: imask = -1 (missing)
      !      - ice or snow:   imask = 0  (from LG field --> binary ice mask)
      !    Define the terrain-type qualifier ztt for each point based on the ice mask values.
      !       ztt = -1  not defined (open water or snow free land)
      !              0  sea-ice     (zlq = 1)
      !              1  snow-covered land (zlq = 0)

      if (lgintrp(kk) < LGthresh ) then 
        imask = -1
      else 
        imask = 0
      end if 

      if ( imask == 0 ) then  ! if ice or snow
        ztt(kk) = 1 - zlq(kk)
      endif

      !    If imask is -1 (no ice/snow), and this is consistent with the model ice
      !    LG value (ie. < LGthresh), and the max MG value indicates ocean (ie.
      !    < MGthresh), then this is a WATEROBS point.

      if ( imask == -1 .and. lgintrp(kk) < LGthresh .and. mgintrp(kk) < MGthresh ) then
      !!!! TO RUN WITHOUT BINARY ICE MASK CHECK (I.E. RELY ON LG ONLY TO REMOVE ICE PTS):
      !if ( lgintrp(kk) < LGthresh .and. mgintrp(kk) < MGthresh ) then
        waterobs(kk) = .true.
      end if

      !   Modify land/sea quailifier if not consistent with waterobs (applies to remote small
      !   islands that should be treated as sea points (for RTTOV)):
      !     -- if waterobs=true, land/sea qualifier should be "sea" value (1)

      if ( waterobs(kk) .and. (zlq(kk) == 0) ) zlq(kk) = 1

    end do

  end subroutine land_ice_mask_ssmis



  !--------------------------------------------------------------------------
  ! wentz_sfctype_ssmis  
  !--------------------------------------------------------------------------

  subroutine wentz_sfctype_ssmis(npts,zlat,zlon,zlq)
    !  Author:   D. Anselmo  MSC/ARMA
    !  Language: FORTRAN 90
    !  Adapted from: wentz_surface_type.ftn90 of process_ssmi
    !
    !  Object:  Determine for each observation point the wentz surface value
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
    !  Version:      Date:      Comment:
    !  --------      -----      --------
    !    0.1       04/01/07     Original adapted code.        D. Anselmo
    !--------------------------------------------------------------------
    !  Variable Definitions
    !  --------------------
    ! npts      - input  -  number of input obs pts in record
    ! zlat      - input  -  array of size ssbg_maxObsNum holding lat values for obs pts
    !                       in record plus undefined pts
    ! zlon      - input  -  array of size ssbg_maxObsNum holding lon values for obs pts
    !                       in record plus undefined pts
    ! zlq       - output -  array holding land/sea qualifier values for all obs
    !                       pts of record
    ! xlat      -internal-  array of size npts holding lat values for obs pts in record
    ! xlon      -internal-  array of size npts holding lon values for obs pts in record
    ! lm        -internal-  array of size 1440x720 holding gridded wentz surface values
    ! wentyp    -internal-  array of size npts holding wentz surface values interpolated to obs pts
    !--------------------------------------------------------------------
    implicit none

    !  Arguments:
    integer, intent(in)  :: npts
    real,    intent(in)  :: zlat(:)
    real,    intent(in)  :: zlon(:)
    integer, intent(out), allocatable :: zlq(:)

    !  Locals:
    integer, parameter :: iunin=40

    integer :: ier,key
    integer :: ni,nj,nk,kk
    integer :: ig1,ig2,ig3,ig4
    integer :: idum1,idum2,idum3
    integer :: idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11
    integer :: idum12,idum13,idum14,idum15,idum16,idum17,idum18

    character(len=12) :: etikxx
    character(len=4)  :: nomvxx
    character(len=2)  :: typxx
    character(len=1)  :: grtyp

    ! F90 allocatable arrays:
    real, allocatable :: xlat(:)
    real, allocatable :: xlon(:)
    real, allocatable :: wentyp(:)
    real, allocatable :: lm(:)
    ! RPNLIB interpolating functions:
    integer :: ezsetopt,ezqkdef
    integer :: gdllsval,gdid
    !--------------------------------------------------------------------
    ! Define FORTRAN FST functions:
    integer, external :: fnom,fclos
    integer, external :: fstinf,fstprm,fstlir
    integer, external :: fstouv,fstfrm


    ! Open Wentz surface field if first call

    ier = fnom( iunin,wentz_file,'STD+RND+R/O',0 )
    ier = fstouv( iunin,'RND' )

    key = fstinf(iunin,ni,nj,nk,-1,' ',0,0,0,' ','LM')
    if ( key <  0 ) then
      call utl_abort(' WENTZ_SFCTYPE_SSMIS: The LM field is MISSING ')
    else
      call utl_reAllocate( lm, ni*nj )
      ier = fstlir(lm,iunin,ni,nj,nk,-1,' ',-1,-1,-1,' ','LM')
    end if

    ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,    &
        &        idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyp,ig1,ig2,  &
        &        ig3,ig4,idum12,idum13,idum14,idum15,idum16,idum17,      &
        &        idum18)

    ! Re-define the grid type.
    !    - L-grid: IG's SHOULD define grid spacing over a REGIONAL
    !              domain. Here, IG's are =0 though.

    !    - A-grid: is for global grid. See RPN documentation on
    !              grid-types for proper IG values.

    grtyp ='A'


    ! Transfer data from zlat,zlon to xlat,xlon.
    ! Ensure proper range for lon values.

    call utl_reAllocate( xlat, npts)
    call utl_reAllocate( xlon, npts)
    call utl_reAllocate( wentyp, npts)
    call utl_reAllocate( zlq, npts)

    zlq(:) = 0
    xlat(:) = zlat(1:npts)
    xlon(:) = zlon(1:npts)
    do kk = 1, npts
      if ( xlon(kk) < 0. ) xlon(kk) = xlon(kk) + 360.
    end do

    ! Interpolate values from LM field to all obs pts of record.

    ier  = ezsetopt('INTERP_DEGREE','VOISIN')
    gdid = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iunin)
    ier  = gdllsval(gdid,wentyp,lm,xlat,xlon,npts)


    ! Define the land/sea qualifier for each point based on wentz surface values.
    do kk = 1, npts
      if ( wentyp(kk) == 0. .or. wentyp(kk) == 6. ) then
        ! wentz = land/coast --> land
        zlq(kk) = 0
      else if ( wentyp(kk) == 4. .or. wentyp(kk) == 5. ) then
        !     wentz = sea/sea-ice --> sea
        zlq(kk) = 1
      else
        call utl_abort(' WENTZ_SFCTYPE_SSMIS: Unexpected Wentz value ')
      end if
    end do
    ier = fstfrm(iunin)
    ier = fclos(iunin)

  end subroutine wentz_sfctype_ssmis

  !--------------------------------------------------------------------------
  ! ssbg_computeSsmisSurfaceType
  !--------------------------------------------------------------------------
  subroutine ssbg_computeSsmisSurfaceType(obsSpaceData)
    !:Purpose:      compute surface type ele and Update obspacedata
    implicit None

    !Arguments
    type(struct_obs),     intent(inout)     :: obsSpaceData           ! obspaceData Object
    ! Locals
    integer                                 :: headerIndex
    real                                    :: obsLatitude(1)
    real                                    :: obsLongitude(1)
    integer, allocatable                    :: landSeaQualifier(:)
    logical                                 :: ssmisDataPresent
    integer                                 :: codtyp
    real,    parameter                      :: satZenithAngle = 53.1 
    real,    parameter                      :: satAzimuthAngle = 210.34


    ssmisDataPresent=.false.
    call obs_set_current_header_list(obsSpaceData,'TO')

    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER0
      codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if ( tvs_isIdBurpInst(codtyp,'ssmis') ) then
        ssmisDataPresent = tvs_isIdBurpInst(codtyp,'ssmis')
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
      call wentz_sfctype_ssmis(1,obsLatitude, obsLongitude, landSeaQualifier)
      call obs_headSet_i(obsSpaceData, OBS_STYP, headerIndex, landSeaQualifier(1))
      call obs_headSet_r(obsSpaceData, OBS_SZA, headerIndex,satZenithAngle)
      call obs_headSet_r(obsSpaceData, OBS_AZA, headerIndex, satAzimuthAngle)
    end do HEADER1

  end subroutine ssbg_computeSsmisSurfaceType
  
  !--------------------------------------------------------------------------
  ! ssbg_grossValueCheck  
  !--------------------------------------------------------------------------

  subroutine ssbg_grossValueCheck(npts,KNO, ztb, ztbThresholdMin, ztbThresholdMax, grossrej)

    !:Purpose: Check Tbs for values that are missing or outside physical limits.
    !          **NOTE: REJECT ALL CHANNELS OF ONE IS FOUND TO BE BAD.
    implicit none

    ! Arguments
    integer, intent(in)               :: npts             ! number of obs pts to process
    integer, intent(in)               :: KNO              ! number of ichannels
    real,    intent(in)               :: ztb(:)           ! bs from input BURP file
    real,    intent(in)               :: ztbThresholdMin  ! ztb threshold for rejection
    real,    intent(in)               :: ztbThresholdMax  ! ztb threshold for rejection
    logical, intent(out), allocatable :: grossrej(:)      ! ogical array defining which obs are to be rejected

    ! Locals
    integer :: ii, indx1, indx2

    call utl_reAllocate(grossrej, npts)
    
    grossrej(1:npts) = .true.
    indx1 = 1
    do ii = 1, npts

      indx2 = ii*KNO
      if ( all( ztb(indx1:indx2) > ztbThresholdMin ) .and. all( ztb(indx1:indx2) < ztbThresholdMax ) ) then
        grossrej(ii) = .false.
      end if
      indx1 = indx2 + 1

    end do

  end subroutine ssbg_grossValueCheck

  !--------------------------------------------------------------------------
  ! ssbg_satqcSsmis
  !--------------------------------------------------------------------------
  subroutine ssbg_satqcSsmis(obsSpaceData, headerIndex, ssmisNewInfoFlag, obsToreject)

    !:Purpose: This program is applied as a first stage of processing to
    !          SSMIS data after it is received from UK MetOffice and
    !          organized into boxes by a program of Jose Garcia. The
    !          processing applied in this program includes:
    !              o  interpolate Wentz surface land mask to each obs pt
    !                 (nearest neighbour) to define land/sea qualifier (008012)
    !              o  interpolate binary ice mask to each obs pt (nearest
    !                 neighbour) to define terrain-type element (013039) where
    !                 0 = sea ice and 1 = snow-covered land
    !              o  interpolate model MG and LG fields to a grid surrounding each obs
    !                 pt to identify obs that are over open water, far from coast/ice
    !              o  identify those obs for which the UKMO rain marker
    !                 is ON (ie. 020029 = 1) indicating poor quality
    !              o  apply a cloud filter to identify those obs in cloudy regions;
    !                 write CLW and IWV (over ocean) to output BURP file
    !              o  re-write data to output BURP file while modifying flags
    !                 for those obs which are not over open water, or have been
    !                 identified in rain/cloud areas, or are of poor quality
    !              o  define satellite zenith angle element (007024) and add
    !                 this and land/sea qualifier and terrain-type elements
    !                 to the output file 
    implicit none
    ! Arguments
    type(struct_obs),     intent(inout) :: obsSpaceData           ! obspaceData Object
    integer,              intent(in)    :: headerIndex
    logical, allocatable, intent(out)   :: obsToreject(:)
    integer, allocatable, intent(out)   :: ssmisNewInfoFlag(:)
    
    !locals
    ! arrays to get from obsspacedata
    character(len=9)                    :: burpFileSatId
    integer, allocatable                :: landSeaQualifier(:)
    real,    allocatable                :: obsLatitude(:)
    real,    allocatable                :: obsLongitude(:)
    real,    allocatable                :: obsTb(:)
    real,    allocatable                :: satZenithAngle(:)
    integer, allocatable                :: terrainType(:)
    real,    allocatable                :: rclw(:)
    real,    allocatable                :: riwv(:)
    real,    allocatable                :: scatw(:)
    integer, allocatable                :: ukRainObs(:)
    ! temporary arrays
    integer, allocatable                :: ier(:)
    real,    allocatable                :: scatl(:)
    logical, allocatable                :: rainDetectionUKMethod(:)
    logical, allocatable                :: grossRej(:)
    logical, allocatable                :: cloudObs(:)
    logical, allocatable                :: iwvreject(:)
    logical, allocatable                :: precipObs(:)
    logical, allocatable                :: waterObs(:)
    real   , allocatable                :: ztb91(:)
    real   , allocatable                :: ztb150(:)
    real   , allocatable                :: ztb_amsub3(:)
    real   , allocatable                :: ztb_amsub5(:)
    real   , allocatable                :: amsubDrynessIndex(:)
    integer, save                       :: numLandObs 
    integer, save                       :: numUkBadObs
    integer, save                       :: numGrossObs
    integer, save                       :: numCloudyObs
    integer, save                       :: numbadIWVObs
    integer, save                       :: numPrecipObs
    integer, save                       :: numTotFilteredObs
    integer, save                       :: numLandScatObs
    integer, save                       :: numSeaScatObs
    integer, save                       :: numDryIndexObs
    integer, save                       :: numSeaIceObs
    integer, save                       :: numScatPrecipObs
    integer, save                       :: numobsF16
    integer, save                       :: numobsF17
    integer, save                       :: numobsF18
    logical, save                       :: ifFirstCall = .true.
    ! variables related to obspacedata
    integer                              :: bodyIndex
    integer                              :: bodyIndexBeg
    integer                              :: bodyIndexEnd
    integer                              :: headerCompt
    integer                              :: currentChannelNumber
    integer                              :: numObsToProcess
    integer                              :: iplatform
    integer                              :: instrum, actualNumChannel
    integer                              :: isat, iplatf
    integer                              :: instr, sensorIndex
    logical                              :: sensorIndexFound
    ! other var
    integer                              :: indx1
    integer                              :: indx2, codtyp
    integer                              :: ii, kk, isatid
    logical                              :: ssmisDataPresent

    ! Check if its ssmis data:
    codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
    ssmisDataPresent = tvs_isIdBurpInst(codtyp,'ssmis')

    if ( .not. ssmisDataPresent ) then
      write(*,*) 'WARNING: WILL NOT RUN ssbg_computeSsmisSurfaceType since no SSMIS DATA is found'
      return
    end if


    ! find tvs_sensor index corresponding to current obs

    iplatf      = obs_headElem_i( obsSpaceData, OBS_SAT, headerIndex )
    instr       = obs_headElem_i( obsSpaceData, OBS_INS, headerIndex )

    call tvs_mapSat( iplatf, iplatform, isat )
    call tvs_mapInstrum( instr, instrum )

    sensorIndexFound = .false.
    do sensorIndex =1, tvs_nsensors
      if ( iplatform ==  tvs_platforms(sensorIndex)  .and. &
           isat      ==  tvs_satellites(sensorIndex) .and. &
           instrum   == tvs_instruments(sensorIndex)       ) then
          sensorIndexFound = .true.
         exit
      end if
    end do
    if ( .not. sensorIndexFound ) call utl_abort('ssbg_satqcSsmis: sensor Index not found')

    ! find actual Number of channels
    actualNumChannel = tvs_coefs(sensorIndex)%coef%fmv_ori_nchn

    headerCompt = 1
    numObsToProcess = 1

    ! Allocate intent out arrays 

    call utl_reAllocate(obsToreject, numObsToProcess*actualNumChannel)
    call utl_reAllocate(ssmisNewInfoFlag, numObsToProcess)

    ! Allocate Fortran working arrays 

    call utl_reAllocate(rainDetectionUKMethod, numObsToProcess)
    call utl_reAllocate(grossRej, numObsToProcess)
    call utl_reAllocate(cloudObs, numObsToProcess)
    call utl_reAllocate(iwvreject, numObsToProcess)
    call utl_reAllocate(precipObs, numObsToProcess)
    call utl_reAllocate(waterObs, numObsToProcess)
    call utl_reAllocate(scatl, numObsToProcess)
    call utl_reAllocate(ier, numObsToProcess)
    call utl_reAllocate(ztb91, numObsToProcess)
    call utl_reAllocate(ztb150, numObsToProcess)
    call utl_reAllocate(ztb_amsub3, numObsToProcess)
    call utl_reAllocate(ztb_amsub5, numObsToProcess)
    call utl_reAllocate(amsubDrynessIndex, numObsToProcess)

    ! ELEMENTS FROM OBSSPACEDATA
    ! Allocate Header elements
    call utl_reAllocate(landSeaQualifier, numObsToProcess)
    call utl_reAllocate(satZenithAngle, numObsToProcess)
    call utl_reAllocate(terrainType, numObsToProcess)
    call utl_reAllocate(obsLatitude, numObsToProcess)
    call utl_reAllocate(obsLongitude, numObsToProcess)
    call utl_reAllocate(rclw, numObsToProcess)
    call utl_reAllocate(riwv, numObsToProcess)
    call utl_reAllocate(scatw, numObsToProcess)
    call utl_reAllocate(ukRainObs, numObsToProcess)
    ! Allocate Body elements
    call utl_reAllocate(obsTb, numObsToProcess*actualNumChannel)
    !initialization
    obsTb(:) = ssbg_realMissing
    riwv(:) = ssbg_realMissing
    scatw(:) = ssbg_realMissing
    
    ! Lecture dans obsspacedata

    burpFileSatId                      = obs_elem_c    ( obsSpaceData, 'STID' , headerIndex )
    rclw(headerCompt)   = obs_headElem_r( obsSpaceData, OBS_CLWO, headerIndex)
    ukRainObs(headerCompt)   = obs_headElem_i( obsSpaceData, OBS_RAIN, headerIndex)
    landSeaQualifier(headerCompt)   = obs_headElem_i( obsSpaceData, OBS_STYP, headerIndex)
    terrainType(headerCompt)     = obs_headElem_i( obsSpaceData, OBS_TTYP, headerIndex)
    ! If terrain type is missing, set it to -1 for the QC programs
    if (terrainType(headerCompt) ==  99) terrainType(headerCompt) = -1
    obsLatitude (headerCompt)          = obs_headElem_r( obsSpaceData, OBS_LAT, headerIndex )
    obsLongitude(headerCompt)          = obs_headElem_r( obsSpaceData, OBS_LON, headerIndex )
    satZenithAngle(headerCompt)        = obs_headElem_r( obsSpaceData, OBS_SZA, headerIndex )
    ! Convert lat/lon to degrees
    obsLongitude(headerCompt) = obsLongitude(headerCompt)*MPC_DEGREES_PER_RADIAN_R8
    if( obsLongitude(headerCompt) > 180. ) obsLongitude(headerCompt) = obsLongitude(headerCompt) - 360.
    obsLatitude(headerCompt)  = obsLatitude(headerCompt) *MPC_DEGREES_PER_RADIAN_R8

    ! To read body elements
    bodyIndexBeg = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    bodyIndexEnd = bodyIndexBeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) -1

    BODY: do bodyIndex =  bodyIndexbeg, bodyIndexEnd
      currentChannelNumber = nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      obsTb(currentChannelNumber)          = obs_bodyElem_r( obsSpaceData,  OBS_VAR, bodyIndex )
    end do BODY

    ! initialization
    if (ifFirstCall) then
      ifFirstCall = .false.
      numLandObs = 0
      numUkBadObs = 0
      numGrossObs = 0
      numCloudyObs = 0
      numbadIWVObs = 0
      numPrecipObs = 0
      numTotFilteredObs = 0
      numLandScatObs = 0
      numSeaScatObs = 0
      numDryIndexObs = 0
      numSeaIceObs = 0
      numScatPrecipObs = 0
      numObsF16 = 0
      numObsF17 = 0
      numObsF18 = 0
    end if 
    waterObs(:) = .false.
    ssmisNewInfoFlag(:) = 0

    ! Record the total number of obs pts read for each satellite.
    ! Set the satellite ID number.
    select case (burpFileSatId(2:7) )
    case ( 'DMSP16' )
      numObsF16 = numObsF16 + numObsToProcess
      isatid = 1
    case ( 'DMSP17' )
      numObsF17 = numObsF17 + numObsToProcess
      isatid = 2
    case ( 'DMSP18' )
      numObsF18 = numObsF18 + numObsToProcess
      isatid = 3
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

    where ( waterobs ) ssmisNewInfoFlag = IBSET(ssmisNewInfoFlag,0)

    !--------------------------------------------------------------------
    ! Determine which obs pts have been flagged for rain (using UKMO method of identifying
    ! bad quality SSMIS data for DMSP16). First, initialize all obs as good.
    !--------------------------------------------------------------------

    rainDetectionUKMethod(:) = .false.
    where ( ukRainObs == 1 ) rainDetectionUKMethod = .true.
    where ( rainDetectionUKMethod ) ssmisNewInfoFlag = IBSET(ssmisNewInfoFlag,10)

    !--------------------------------------------------------------------
    ! Check for values of TB that are missing or outside physical limits.
    ! **NOTE: REJECT ALL CHANNELS IF ONE IS FOUND TO BE BAD.
    !--------------------------------------------------------------------
    !    write(6,*) 'Applying gross value filter.'
    grossRej(:)  = .false.
    call ssbg_grossValueCheck(numObsToProcess, ssbg_maxNumChan, obsTb,50., 400., grossRej)
    where ( grossRej ) ssmisNewInfoFlag = IBSET(ssmisNewInfoFlag,11)

    !--------------------------------------------------------------------
    ! Apply a CLW regression technique to determine cloudy obs pts.
    ! Apply a IWV regression technique to determine IWV obs out of bounds.
    ! Detect Sea_Ice percent from Tb(Ta) and set waterobs() flag to .FALSE. if
    !  amount > 70%
    ! To begin, assume that all obs are good.
    !--------------------------------------------------------------------

    !    write(6,*) 'Applying cloud filter.'
    cloudobs(:)  = .false.
    iwvreject(:) = .false.
    precipobs(:) = .false.

    ! --------------------------------------------------------------------

    call cld_filter_fweng(numObsToProcess, obsTb, alg_option, waterObs, &
                          grossRej, cloudObs, iwvReject, precipObs, rclw, riwv, isatid,&
                          obsLatitude, numSeaIceObs)
    !
    !  NOTE: rclw, riwv missing value is zmisg=9.9e09
    !    --> if ( (clw > clw_thresh) .or. precipobs(ii) ) cloudobs(ii) = .true.
    !    --> if ( clw == zmisg )  cloudobs(ii) = .true.
    !
    where ( iwvReject ) ssmisNewInfoFlag = IBSET(ssmisNewInfoFlag,5)
    where ( precipObs ) ssmisNewInfoFlag = IBSET(ssmisNewInfoFlag,4)
      
    !   Extract Tb for channels 18 (91H GHz)  and 8 (150H GHz) for Bennartz SI
    !   Extract Tb for channels 11 (AMSU-B 3) and 9 (AMSU-B 5) for Dryness Index (DI)

    indx1 = 1
    do ii = 1, numObsToProcess
      indx2 = ii*ssbg_maxNumChan
      ztb91(ii)      = obsTb(indx1+17)
      ztb150(ii)     = obsTb(indx1+7)
      ztb_amsub3(ii) = obsTb(indx1+10)
      ztb_amsub5(ii) = obsTb(indx1+8)
      indx1 = indx2 + 1
    end do

    !---------------------------------------------------------------------------------
    ! Compute the Bennartz scattering index (SI) for each point from channels 8 and 18
    ! (detects precipitation affected AMSU-B-like channel radiances)
    ! -- using constant satellite zenith angle = 53.1
    ! -- SSMIS channel 18 (91.655 GHz) is used for AMSU-B like channel 1 (89 GHz)
    ! -- SSMIS channel  8 is used for AMSU-B like channel 2 (both at 150 GHz)
    !---------------------------------------------------------------------------------
    call bennartz(ier, numObsToProcess, ztb91, ztb150, satZenithAngle, landSeaQualifier, & 
                 scatl, scatw)

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
    !   ukbadobs()  = .true. if UKMetO QC sets "rain" flag: solar-intrusion or other anomaly
    !   grossrej()  = .true. if any channel had a gross error at the point
    !   cloudobs()  = .true. if CLW > 0.01 kg/m**2 or precip (over water)
    !   precipobs() = .true. if precip. detected (CLW=missing)
    !   waterobs()  = .true. if open water point (far from land or sea ice)
    !   iwvreject() = .true. if IWV > 80 kg/m**2
    !--------------------------------------------------------------------

    obsToReject(:) = .false.
    do kk = 1, numObsToProcess
      
      !      Reject all channels if UKMet rain flag or gross Tb error detected   
      if ( rainDetectionUKMethod(kk) .or. grossRej(kk) ) then
        !      "BAD" Observations       
        obsToReject(:) = .true.
        numTotFilteredObs = numTotFilteredObs + 1
        !      Set INFO flag bits (CLW, BSI, DI) for all points without gross data error      
        if ( .not. grossRej(kk) ) then
          !----------------------------------------------------------------------
          if  ( .not. waterobs(kk) ) then   !  LAND or SEA-ICE
            !        Check BSI and DI for AMSU-B 3-4
            !         Bennartz Scattering Index
            !         Land point           
            if ( scatl(kk) > 0.0  .and. scatl(kk) /= ssbg_realMissing ) ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),6)
            !         Sea-ice point 
            if ( scatw(kk) > 40.0 .and. scatw(kk) /= ssbg_realMissing ) ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),6)
            !         Dryness index           
            if ( amsubDrynessIndex(kk) > 0.0 )   ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),9)
            if ( amsubDrynessIndex(kk) > -10.0 ) ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),8)
            if ( amsubDrynessIndex(kk) > -20.0 ) ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),7)
            !----------------------------------------------------------------------
          else                              !  OPEN WATER
          !----------------------------------------------------------------------
            if (cloudobs(kk)) then
              if ( rclw(kk) /= ssbg_realMissing ) ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),1)
              !            AMSU-A like channels 3-5(6)
              if ( rclw(kk) > clw_amsu_rej .and. rclw(kk) /= ssbg_realMissing ) ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),3)
              if ( rclw(kk) > clw_amsu_rej_ch3 .and. rclw(kk) /= ssbg_realMissing ) ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),2)
            endif
            !         Bennartz Scattering Index ( AMSU-B 2-5)
            !         Open water point
            if ( scatw(kk) > 15.0 .and. scatw(kk) /= ssbg_realMissing ) ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),6)
          endif    ! if not waterobs

        endif    ! if not grossrej

      else     ! if ( ukbadobs(kk) .or. grossrej(kk) )
        !      "GOOD" Observations       
        !-----------------------------------------------------------------------------------       
        !      OVER LAND OR SEA-ICE, 
        !         -- reject lower tropospheric channels and imager channels
        !         -- check Bennartz SI and DI for AMSU-B like channels
        !----------------------------------------------------------------------------------- 
        !        -- CLW/PRECIP not determined over land
        !        -- surface emissivity effects lower tropospheric channels     
        if  ( .not. waterobs(kk) ) then
          obsToReject(1:ipc) = .true.    ! AMSU-A 3-5(6)
          obsToReject(12:18) = .true.    ! SSMI-like imager 1-7
          obsToReject(8:9) = .true.      ! AMSU-B 2,5
          !        Check BSI and DI for AMSU-B 3-4
          !         Bennartz Scattering Index
          !         Land point           
          if ( scatl(kk) > 0.0  .and. scatl(kk) /= ssbg_realMissing )   obsToReject(10:11) = .true.
          !         Sea-ice point 
          if ( scatw(kk) > 40.0 .and. scatw(kk) /= ssbg_realMissing )   obsToReject(10:11) = .true.
          !         Missing scattering index
          if ( scatw(kk) == ssbg_realMissing .and. scatl(kk) == ssbg_realMissing ) obsToReject(10:11) = .true.
          if ( any(obsToReject(10:11))) then
            numLandScatObs = numLandScatObs + 1
            ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),6)
          endif
          !         Dryness index           
          if ( amsubDrynessIndex(kk) > 0.0 ) then
            obsToReject(11) = .true.
            ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),9)
          endif
          if ( amsubDrynessIndex(kk) > -10.0 ) then
            obsToReject(10) = .true.
            ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),8)
          endif
          if ( amsubDrynessIndex(kk) > -20.0 ) then
            obsToReject(9) = .true.
            ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),7)
          endif
          if ( amsubDrynessIndex(kk) > -10.0 ) numDryIndexObs = numDryIndexObs + 1
          numTotFilteredObs = numTotFilteredObs + 1
        endif
        !-----------------------------------------------------------------------------------         
        !      OVER WATER, 
        !        -- reject tropospheric channels and imager channels if cloudy, 
        !           precip, or excessive IWV
        !        -- check Bennartz SI for AMSU-B like channels
        !-----------------------------------------------------------------------------------    
        if  ( waterobs(kk) ) then
          if (cloudobs(kk) .or. iwvreject(kk))  then
            if ( iwvreject(kk) ) then
              !              SSMI-like imager channels
              obsToReject(12:18) = .true.
              !             AMSU-A like channels 3-5(6)
              obsToReject(1:ipc) = .true.
            else  !  ----- CLOUDY OBSERVATION (OR MISSING CLOUD) ------
              !              SSMI-like imager channels
              obsToReject(12:18) = .true.
              if ( rclw(kk) /= ssbg_realMissing ) ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),1)
              !              AMSU-A like channels 3-5(6)
              if ( rclw(kk) > clw_amsu_rej ) then
                obsToReject(2:ipc) = .true.
                if ( rclw(kk) /= ssbg_realMissing ) ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),3)
              endif
              if ( rclw(kk) > clw_amsu_rej_ch3  ) then
                obsToReject(1) = .true.
                if ( rclw(kk) /= ssbg_realMissing ) ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),2)
              endif
            endif
            numTotFilteredObs = numTotFilteredObs + 1
          endif
          !        Check BSI for AMSU B channels 2-5 for all water obs
          !        Bennartz Scattering Index ( AMSU-B 2-5)
          !        Open water point
          if ( scatw(kk) > 15.0 .and. scatw(kk) /= ssbg_realMissing ) then
            obsToReject(8:11) = .true.
            ssmisNewInfoFlag(kk) = IBSET(ssmisNewInfoFlag(kk),6)
            numSeaScatObs = numSeaScatObs + 1
            if ( precipobs(kk) ) numScatPrecipObs = numScatPrecipObs + 1
          endif

        endif   ! if waterobs

      endif   ! if ( ukbadobs(kk) .or. grossrej(kk) [end else] )

      if ( .not. waterobs(kk) ) numLandObs  = numLandObs  + 1
      if ( rainDetectionUKMethod(kk) )  numUkBadObs = numUkBadObs + 1
      if ( grossrej(kk) )  numGrossObs = numGrossObs + 1
      if ( cloudobs(kk)  .and. waterobs(kk) ) numCloudyObs  = numCloudyObs  + 1
      if ( iwvreject(kk) .and. waterobs(kk) .and.   &
      &         (.not. cloudobs(kk)) ) numbadIWVObs  = numbadIWVObs  + 1
      if ( precipobs(kk) .and. waterobs(kk) ) then
        numPrecipObs = numPrecipObs + 1
      endif
    end do
      
    !-------------------------------------------------------------------------------
    ! Update ObsspaceData with the computed values of: ZLQ, ZTT, IWV and scatw
    !-------------------------------------------------------------------------------
    call obs_headSet_i(obsSpaceData, OBS_STYP, headerIndex, landSeaQualifier(1))
    call obs_headSet_i(obsSpaceData, OBS_TTYP, headerIndex, terrainType(1))
    call obs_headSet_r(obsSpaceData, OBS_SCAT, headerIndex, scatw(1))
    call obs_headSet_r(obsSpaceData, OBS_CLWO, headerIndex, rclw(1))
    call obs_headSet_r(obsSpaceData, OBS_IWV,  headerIndex, riwv(1))

    if (headerIndex == obs_numHeader(obsSpaceData)) then
      write(*,*) '*******************************************************************'
      write(*,*) '******************* SATQC PROGRAM STATS****************************'
      write(*,*) '*******************************************************************'
      write(*,*)
      write(*,*) 'Number of cloudy obs = ', numCloudyObs 
      write(*,*) 'Number of precip obs = ', numPrecipObs
      write(*,*) 'Number of gross rej obs = ',  numGrossObs
      write(*,*) 'Number of land  obs = ',numLandObs  
      write(*,*) 'Number of dry obs = ', numDryIndexObs 
      write(*,*) 'Number of Sea Ice obs = ',numSeaIceObs
      write(*,*) 'Number of F16 obs = ', numObsF16 
      write(*,*) 'Number of F17 obs = ', numObsF17
      write(*,*) 'Number of F18 obs = ', numObsF18
      write(*,*) '*******************************************************************'
    end if 
  end subroutine ssbg_satqcSsmis 

  !--------------------------------------------------------------------------
  ! ssbg_updateObsSpaceAfterSatQc
  !--------------------------------------------------------------------------
  subroutine ssbg_updateObsSpaceAfterSatQc(obsSpaceData, headerIndex,obsToreject) 

    !:Purpose:      Update obspacedata variables (obstTB and obs flags) after QC
    implicit None

    !Arguments
    type(struct_obs),     intent(inout)     :: obsSpaceData           ! obspaceData Object
    integer,              intent(in)        :: headerIndex            ! current header index
    logical,              intent(in)        :: obsToReject(:)            ! data flags
    ! Locals
    integer, allocatable                    :: obsFlags(:)
    integer, allocatable                    :: obsGlobalFlag(:)
    integer, allocatable                    :: satScanPosition(:)
    integer                              :: bodyIndex
    integer                              :: bodyIndexBeg
    integer                              :: bodyIndexEnd
    integer                              :: headerCompt
    integer                              :: actualNumChannel
    integer                              :: currentChannelNumber
    integer                              :: numObsToProcess
    integer                              :: iplatform
    integer                              :: instrum
    integer                              :: isat, iplatf
    integer                              :: instr
    logical                              :: sensorIndexFound
    integer                              :: sensorIndex
    integer                              :: ntIndex 
    integer                              :: dataIndex
    integer                              :: channelIndex

   ! find tvs_sensor index corresponding to current obs

    iplatf      = obs_headElem_i( obsSpaceData, OBS_SAT, headerIndex )
    instr       = obs_headElem_i( obsSpaceData, OBS_INS, headerIndex )

    call tvs_mapSat( iplatf, iplatform, isat )
    call tvs_mapInstrum( instr, instrum )

    sensorIndexFound = .false.
    do sensorIndex =1, tvs_nsensors
      if ( iplatform ==  tvs_platforms(sensorIndex)  .and. &
           isat      ==  tvs_satellites(sensorIndex) .and. &
           instrum   == tvs_instruments(sensorIndex)       ) then
          sensorIndexFound = .true.
         exit
      end if
    end do
    if ( .not. sensorIndexFound ) call utl_abort('ssbg_updateObsSpaceAfterSatQc: sensor Index not found')

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

    obsGlobalFlag(headerCompt) = obs_headElem_i( obsSpaceData, OBS_ST1, headerIndex )
    satScanPosition(headerCompt)  = obs_headElem_i( obsSpaceData, OBS_FOV , headerIndex)
    bodyIndexBeg = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    bodyIndexEnd = bodyIndexBeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) -1

    BODY2: do bodyIndex =  bodyIndexbeg, bodyIndexEnd
      currentChannelNumber = nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      obsFlags(currentChannelNumber)       = obs_bodyElem_i( obsSpaceData,  OBS_FLG, bodyIndex )
    end do BODY2
    ! Modify flags
    
    !------------------------------------------------------------------ 
    ! 1 - Mark global flags if any obsToReject
    !------------------------------------------------------------------ 
    do ntIndex = 1 , numObsToProcess
      if ( any(obsToReject) ) obsGlobalFlag(ntIndex) = ibset(obsGlobalFlag(ntIndex), 6)
    end do 

    !------------------------------------------------------------------ 
    ! 2 - Mark obs flags for each value of obsToReject
    !------------------------------------------------------------------ 
    dataIndex = 0
    do ntIndex = 1 , numObsToProcess
      do channelIndex = 1, actualNumChannel
        dataIndex = dataIndex+1 
        if (resetQc) obsFlags(dataIndex) = 0
        if (obsToReject(dataIndex)) obsFlags(dataIndex) = ibset(obsFlags(dataIndex),7)
      end do 
    end do

    !-----------------------------------------------------------------
    !    Subtract 270 from FOV values (element 005043).
    !-----------------------------------------------------------------
    
    do ntIndex = 1 , numObsToProcess
      if (satScanPosition(ntIndex) > 270) satScanPosition(ntIndex) = satScanPosition(ntIndex) - 270
    end do
 
    ! write elements in obsspace
    call obs_headSet_i(obsSpaceData, OBS_FOV, headerIndex, satScanPosition(1))
    call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex, obsGlobalFlag(1))
    bodyIndexBeg = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    bodyIndexEnd = bodyIndexBeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) -1
    BODY: do bodyIndex =  bodyIndexbeg, bodyIndexEnd
      currentChannelNumber=nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
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

    !Arguments
    type(struct_obs),     intent(inout) :: obsSpaceData      ! obspaceData Object
    integer,              intent(in)    :: headerIndex       ! header index
    integer, allocatable, intent(out)   :: flagsInovQc(:)    ! flags for assimilation/rejection of obs

    ! Locals
    character(len=9)     :: burpFileSatId            ! Satellite ID

    integer              :: actualNumChannel         ! actual Num channel
    integer              :: bodyIndex
    integer              :: bodyIndexBeg
    integer              :: bodyIndexEnd 
    integer              :: codtyp                   ! code type
    integer              :: channelIndex
    integer              :: currentChannelNumber
    integer              :: headerCompt
    integer              :: instr
    integer              :: instrum
    integer              :: iPlatf
    integer              :: iPlatform
    integer              :: iSat
    integer              :: numObsToProcess          ! number of obs in current report
    integer              :: sensorIndex              ! find tvs_sensor index corresponding to current obs

    integer, allocatable :: obsChannels(:)           ! channel numbers
    integer, allocatable :: obsFlags(:)              ! data flags

    logical              :: sensorIndexFound
    logical              :: ssmisDataPresent

    real   , allocatable :: modelInterpTer(:)        ! topo in standard file interpolated to obs point
    real   , allocatable :: obsLatitude(:)           ! obs. point latitude
    real   , allocatable :: obsLongitude(:)          ! obs. point longitude
    real   , allocatable :: ompTb(:)                 ! OMP values

    ! Temporary arrays
    integer :: rejcnt(ssbg_maxNumChan,ssbg_maxNumSat)
    integer :: rejcnt2(ssbg_maxNumChan,ssbg_maxNumSat)
    integer :: totobs(ssbg_maxNumChan,ssbg_maxNumSat)
    integer :: totobs2(ssbg_maxNumChan,ssbg_maxNumSat)

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

    iPlatf      = obs_headElem_i( obsSpaceData, OBS_SAT, headerIndex )
    instr       = obs_headElem_i( obsSpaceData, OBS_INS, headerIndex )

    call tvs_mapSat( iPlatf, iPlatform, iSat )
    call tvs_mapInstrum( instr, instrum )

    sensorIndexFound = .false.
    do sensorIndex =1, tvs_nsensors
      if ( iPlatform ==  tvs_platforms(sensorIndex)  .and. &
           iSat      ==  tvs_satellites(sensorIndex) .and. &
           instrum   == tvs_instruments(sensorIndex)       ) then
        sensorIndexFound = .true.
        exit
      end if
    end do
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

    ! Lecture dans obsspacedata (entiers)

    burpFileSatId = obs_elem_c    ( obsSpaceData, 'STID' , headerIndex )
    bodyIndexBeg  = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    bodyIndexEnd  = bodyIndexBeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) -1

    ! Lecture dans obsspacedata (reels)

    obsLatitude (headerCompt)          = obs_headElem_r( obsSpaceData, OBS_LAT, headerIndex )
    obsLongitude(headerCompt)          = obs_headElem_r( obsSpaceData, OBS_LON, headerIndex )
    BODY: do bodyIndex =  bodyIndexbeg, bodyIndexEnd
      currentChannelNumber = nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      ompTb(currentChannelNumber)          = obs_bodyElem_r( obsSpaceData,  OBS_OMP, bodyIndex )
      obsFlags(currentChannelNumber)       = obs_bodyElem_i( obsSpaceData,  OBS_FLG, bodyIndex )
    end do BODY

    do channelIndex=1, actualNumChannel
      obsChannels(channelIndex)    = channelIndex+tvs_channelOffset(sensorIndex)
    end do

    ! Convert lat/lon to degrees
    obsLongitude(headerCompt) = obsLongitude(headerCompt)*MPC_DEGREES_PER_RADIAN_R8
    if( obsLongitude(headerCompt) > 180. ) obsLongitude(headerCompt) = obsLongitude(headerCompt) - 360.
    obsLatitude(headerCompt)  = obsLatitude(headerCompt) *MPC_DEGREES_PER_RADIAN_R8

    !  Initialize the counters which tabulate the rejected observation
    !  statistics for each channel and each satellite.
    rejcnt(:,:)  = 0
    totobs(:,:)  = 0
    rejcnt2(:,:) = 0
    totobs2(:,:) = 0

    !--------------------------------------------------------------------
    ! 3) Extract model topography from the input GEO file and interpolate
    !    to all observation points in the box
    !--------------------------------------------------------------------

    call ssbg_readGeophysicFieldsAndInterpolate(obsLatitude,obsLongitude,modelInterpTer)

    !--------------------------------------------------------------------
    ! 4) Perform quality control on the data, checking O-P and topography
    !--------------------------------------------------------------------

    call check_stddev(obsChannels,ompTb,flagsInovQc,actualNumChannel,numObsToProcess,sensorIndex,  &
        &             burpFileSatId,ssbg_maxNumChan,oer_tovutil,oer_toverrst,rejcnt,totobs,obsFlags)

    call check_topo(modelInterpTer,flagsInovQc,actualNumChannel,numObsToProcess,sensorIndex, &
        &           burpFileSatId,rejcnt2,totobs2)

  end subroutine ssbg_inovqcSsmis

  !--------------------------------------------------------------------------
  ! check_stddev
  !--------------------------------------------------------------------------
  subroutine check_stddev(kcanomp,ptbomp,icheck,knomp,knt,knosat,  &
       &                  stnid,mxchn,iutilst,sdstats,rejcnt,totobs,kflags)
    !--------------------------------------------------------------------
    !
    !  This subroutine has been adapted from the following original work:
    !  TITLE:  subroutine TOVCHECK_AMSUA of bgck.satqc_amsua.f
    !  AUTEUR: Jacques Halle (CMC)
    !
    !  Language: FORTRAN 90
    !
    !  Object: Perform quality control on the radiances by analysing the
    !          magnitude of the residuals.
    !
    !  Call:   call check_stddev(icanomp,zomp,kchkprf,icheck,inomp,nt,inosat
    !                            stnid,mxchn,iutilst,sdstats,rejcnt,totobs,iflags)
    !
    !  History
    !  Version:      Date:      Comment:
    !  --------      -----      --------
    !    0.1       05/12/02     Adapted code.                D. Anselmo
    !    0.2       24/03/06     Added check on UTIL values in statistics
    !                           file to systematically filter channels.
    !                                                       D. Anselmo
    !    0.3       26/09/07     Adapted for SSMIS           S. Macpherson
    !    0.4       11/10/07     Modifications for AMSU-like channels
    !                                                       S. Macpherson
    !    0.5       13/02/09     Change to icheck values:
    !                           icheck=1 data not checked
    !                           icheck=2 for UTIL + uncorrected data
    !                                                       S. Macpherson
    !    0.6       06/01/11     Flag bit 7 replaces 9 in satqc_ssmis files
    !                                                       S. Macpherson
    !    0.7       06/11/12     Adapt for new format stats file
    !                                                       S. Macpherson
    !------------------------------------------------------------------
    ! Variable Definitions:
    ! ---------------------
    !   kcanomp - input  -  channel numbers (1-mxchn) (residuals)
    !   ptbomp  - input  -  residus (o-p) (radiance residuals)
    !                       --> ZOMP in subroutine readdata
    !   kflags  - input  -  radiance data flags (bit 7 on for data that
    !                       will not be assimilated)
    !                       --> IFLAGS in subroutine readdata
    !   icheck  - output -  quality contol indicator for each channel of each
    !                       observation point
    !                          =0  ok
    !                          =1  not checked because FLAG bit 7 ON,
    !                          =2  reject by UTIL value or FLAG bit 6 OFF (data not bias corrected),
    !                          =3  reject by rogue check,
    !                          =4  reject by both UTIL value and rogue check.
    !                       icheck is a 1-D array, with dimensions KNT*KNOMP
    !                         ==> individual channels may be rejected from each obs pt
    !   knomp   - input  -  number of residual channels
    !   knt     - input  -  number of groups of NVAL*NELE
    !                       = nt; extracted from bloc header in XTRBLK (ie # of obs pts!!!)
    !   knosat  - input  -  number of satellite (index # --> 1-nsat)
    !                       determined by stnid: DMSP13=1,DMSP14=2,DMSP15=3,DMSP16=4,etc...
    !   stnid   - input  -  identificateur du satellite
    !   mxchn   - input  -  number of channels for SSM/I or SSMIS
    !   iutilst - input  -  UTIL values read from the total error statistics file
    !                         0 = blacklisted channel
    !                         1 = assimilated channel
    !   sdstats - input  -  standard deviation statistics read from the total
    !                       error statistics file
    !   rejcnt  - in/out -  counter of the rejected obs for each channel
    !                       of each satellite (dimension: knomp*mxsat)
    !   totobs  - in/out -  counter of the total obs checked for each channel
    !                       of each satellite (dimension: knomp*mxsat)
    !  roguefac -internal-  constant which determines the severity of the
    !                       quality control applied to the O-P magnitude
    !                       for each channel
    ! xcheckval -internal-  product of roguefac and standard deviation for
    !                       each channel
    !------------------------------------------------------------------
    implicit none

    !  Arguments:
    integer, intent(in) :: mxchn
    integer, intent(in) :: knomp,knt,knosat

    integer, intent(in)    :: iutilst(:,:)
    integer, intent(in)    :: kcanomp(:)
    integer, intent(in)    :: kflags(:)
    integer, intent(out)   :: icheck(:)
    integer, intent(inout) :: rejcnt(:,:)
    integer, intent(inout) :: totobs(:,:)

    real(8), intent(in) :: sdstats(:,:)
    real   , intent(in) :: ptbomp(:)

    character(len=9), intent(in) :: stnid

    !  Locals:
    integer :: ji,jj,ichn,indx1,indx2,iich1,iich8,iich4,iich11
    integer :: j

    real :: roguefac(mxchn)
    real :: xcheckval

    ! Define max Abs(O-P) for channel 8 for rejection of channels 8-11
    ! -- should be consistent with bgck.satqc_amsub.f (AMSU-B)
    real, parameter :: xompch8 = 5.0   ! units = K

    ! Define factor for channel 1 O-P for rejection of channels 1-4
    ! -- should be consistent with bgck.satqc_amsua.f (AMSU-A)
    real, parameter :: xfacch1 = 2.0   ! xfacch1*errtot
    !  Initialization. Assume all observations are to be kept.

    ! roguefac(1:7) should be consistent with bgck.satqc_amsua.f (for AMSUA ch. 3-10)
    ! roguefac(8:11) should be consistent with bgck.satqc_amsub.f (for AMSUB ch. 2-5)
    ! roguefac(12:18) should be consistent with ssmi_inovqc.ftn90 (for SSMI ch. 1-7)
    ! roguefac(19:24) are for upper-atmospheric channels (strato/meso-sphere)

    roguefac = (/2.0,3.0,4.0,4.0,4.0,4.0,4.0,2.0,4.0,4.0,4.0,2.0,2.0,2.0, &
         &             2.0,2.0,2.0,2.0,3.0,3.0,3.0,3.0,3.0,3.0/)

    icheck(:) = 0

    !  Loop through all observation points in the current record and check
    !  first whether a channel is to be systematically removed and second
    !  whether the O-P value for each channel falls within the acceptable
    !  range defined by roguefac and sdstats (knosat identifies the satellite).

    !  *** Only do check for data that could be assimilated ***

    !  Also, apply multi-channel rejections for AMSU-like channels using
    !  O-P values from channels 1 and 8.
    !     If Abs(O-P) chan. 1 (AMSU-A 3) > 2*errtot, reject ch. 1-4 (AMSU-A 3-6)
    !     If Abs(O-P) chan. 8 (AMSU-B 2) > 5K, reject ch. 8-11 (AMSU-B 2-5)

    !--------------------------------------------------------------------------
    !  1. Rogue check O-P for all data flagged for assimilation
    !--------------------------------------------------------------------------

    do ji = 1, knt*knomp

      ichn = kcanomp(ji)
      xcheckval = roguefac(ichn) * sdstats(ichn,knosat)

      if ( .not. btest(kflags(ji),7) ) then

        totobs(ichn,knosat) = totobs(ichn,knosat) + 1

        if ( iutilst(ichn,knosat) == 0 .or. (.not. btest(kflags(ji),6)) ) then

          ! systematic rejection of this channel
          icheck(ji) = 2
          if ( abs( ptbomp(ji) ) >= xcheckval ) then
            ! rogue check failure as well
            icheck(ji) = 4
          endif

        else

          ! keep channel, now perform rogue check
          if ( abs( ptbomp(ji) ) >= xcheckval ) then
            icheck(ji) = 3
          endif

        end if

      else

        icheck(ji) = 1

      endif


      !  Keep statistics of obs rejected by rogue check.

      if ( icheck(ji) > 2 ) then
        rejcnt(ichn,knosat) = rejcnt(ichn,knosat) + 1
        if ( ssbg_debug ) then
          write(6,*) ' CHECK_STDDEV: '
          write(6,*) stnid(2:9),' Rogue check reject: ',   &
               &  ' Obs = ',ji,' Channel = ',ichn,         &
               &  ' Check value = ',xcheckval,             &
               &  ' O-P = ',ptbomp(ji)
        end if
      endif

    end do

    !----------------------------------------------------------------------------------
    ! 2. Apply additional criteria for rejections of multiple AMSU-like channels using
    !    Channels 1 and 8 Abs(O-P) (for data that could be assimilated [bit 7 off])
    !----------------------------------------------------------------------------------
    indx1 = 0
    indx2 = 0
    ! Loop over observation points (knomp = num. channels [24])
    do jj = 1, knt

      indx1 = jj*knomp     ! index of ch.24 obs
      indx2 = indx1 - (knomp-1)
      iich1  = indx2       ! index of ch.1 obs
      iich8  = indx2 + 7   ! index of ch.8 obs
      iich4  = iich1 + 3   ! index of ch.4 obs
      iich11 = iich8 + 3   ! index of ch.11 obs

      ! AMSU-A like ch. 1-4
      if ( sum(icheck(iich1:iich4)) /= 4 ) then
        xcheckval = xfacch1 * sdstats(1,knosat)
        if ( abs(ptbomp(iich1)) >= xcheckval ) then
          do j = iich1, iich4
            ichn = kcanomp(j)
            if (icheck(j) /= 1) then
              if ( icheck(j) < 3 ) rejcnt(ichn,knosat) = rejcnt(ichn,knosat) + 1
              icheck(j) = max(icheck(j),3)
            endif
          end do
        endif
      endif

      ! AMSU-B like ch. 8-11
      if ( sum(icheck(iich8:iich11)) /= 4 ) then
        xcheckval = xompch8
        if ( abs(ptbomp(iich8)) >= xcheckval ) then
          do j = iich8, iich11
            ichn = kcanomp(j)
            if (icheck(j) /= 1) then
              if ( icheck(j) < 3 ) rejcnt(ichn,knosat) = rejcnt(ichn,knosat) + 1
              icheck(j) = max(icheck(j),3)
            endif
          end do
        endif
      endif

    end do  ! end Loop over observation points

  end subroutine check_stddev

  !--------------------------------------------------------------------------
  ! check_topo
  !--------------------------------------------------------------------------
  subroutine check_topo(mtintrp,icheck,knomp,knt,knosat,  &
       &                stnid,rejcnt,totobs)
    !--------------------------------------------------------------------
    !
    !  This subroutine has been adapted from the following original work:
    !  TITLE:  subroutine TOVCHECK_AMSUA of bgck.satqc_amsua.f
    !  AUTEUR: Jacques Halle (CMC)
    !
    !   --- Based on subroutine check_stddev in this file ---
    !
    !  Language: FORTRAN 90
    !
    !  Object: Perform rejection of observations for selected channels based
    !          on model surface height (for channels assimilated over land)
    !          -- for a single satellite (stnid,knosat)
    !          -- for a single box (nt observations)
    !
    !------------------------------------------------------------------
    ! Variable Definitions:
    ! ---------------------
    !   mtintrp - input  -  model surface height (m) at each observation point
    !   icheck  - in/out -  quality contol indicator for each channel of each
    !                       observation point
    !                         -- on INPUT
    !                          =0  ok
    !                          =1  not checked because FLAG bit 7 ON,
    !                          =2  reject by UTIL value or bit 6 OFF (data not bias corrected),
    !                          =3  reject by rogue check,
    !                          =4  reject by both UTIL value and rogue check.
    !                         -- on OUTPUT,
    !                          =0  ok
    !                          =1  not checked because FLAG bit 7 ON,
    !                          =2  reject by UTIL value,
    !                          =3  reject by rogue check,
    !                          =4  reject by both UTIL value and rogue check.
    !                          =5  rejection due to topography
    !                          =6  rejection due to topography, reject by UTIL value
    !                          =7  rejection due to topography, reject by rogue check
    !                          =8  rejection due to topography, reject by UTIL value, reject by rogue check
    !                       icheck is a 1-D array, with dimensions KNT*KNOMP
    !                         ==> individual channels may be rejected from each obs pt
    !   knomp   - input  -  number of residual channels
    !   knt     - input  -  number of groups of NVAL*NELE
    !                       = nt; extracted from bloc header in XTRBLK (ie # of obs pts!!!)
    !                       *** NOT CURRENTLY USED ***
    !   knosat  - input  -  number of satellite (index # --> 1-nsat)
    !                       determined by stnid: DMSP13=1,DMSP14=2,DMSP15=3,DMSP16=4,etc...
    !   stnid   - input  -  identificateur du satellite
    !   rejcnt  - in/out -  counter of the rejected obs for each channel
    !                       of each satellite (dimension: knomp*mxsat)
    !   totobs  - in/out -  counter of the total obs checked for each channel
    !                       of each satellite (dimension: knomp*mxsat)
    !------------------------------------------------------------------
    implicit none

    !  Arguments:
    integer, intent(in) :: knomp,knt,knosat

    integer, intent(inout) :: icheck(:)         ! dimension (knt*knomp)
    integer, intent(inout) :: rejcnt(:,:)
    integer, intent(inout) :: totobs(:,:)

    real, intent(in) :: mtintrp(:)  ! dimension (knt)

    character(len=9), intent(in)    :: stnid

    !  Locals:

    integer, parameter :: nch2chk=4             ! number of channels to check

    integer  ::  ji,jj,indx1,indx2,ii
    integer  ::  itrejcnt

    real     :: zcheckval, zmt

    real     :: topolimit(nch2chk)
    integer  :: mchan(nch2chk)

    logical :: debug2

    !------------------------------------------------------------------
    ! Define channels to check and height limits (m) for rejection
    !------------------------------------------------------------------
    mchan     = (/    4,     9,    10,    11 /)
    topolimit = (/  250., 1000., 2000., 2500. /)

    if (ssbg_debug) then
       write(6,*) ' CHECK_TOPO: Maximum input MT = ', maxval(mtintrp)
       write(6,*) ' CHECK_TOPO: knt = ', knt
       write(6,*) ' CHECK_TOPO: First rejections in box. StnID = ', stnid(2:9)
    end if

    !------------------------------------------------------------------
    ! Perform the check for the selected channels
    !------------------------------------------------------------------

    !   Loop over all knt observation locations

    indx1 = 0
    indx2 = 0
    itrejcnt = 0
    do jj = 1, knt
      debug2 = .false.
      indx1 = jj*knomp
      indx2 = indx1 - (knomp-1)   ! channel     1 index
      zmt = mtintrp(jj)           ! model topography height [m]
      if (ssbg_debug) then
        if (zmt == maxval(mtintrp) .and. zmt > minval(topolimit)) then
          write(6,*) ' CHECK_TOPO: ****** Max height point! zmt = ', zmt
          debug2 = .true.
        endif
      endif
      do ji = 1, nch2chk        ! loop over channels to check (mchan)
        zcheckval = topolimit(ji)   ! height limit [m] for channel mchan(ji)
        ii = indx2 + (mchan(ji)-1)  ! channel mchan(ji) index
        if ( icheck(ii) /= 1 ) then
          totobs(mchan(ji),knosat) = totobs(mchan(ji),knosat) + 1
          if ( zmt > zcheckval ) then
            icheck(ii) = max(1,icheck(ii)) + 4
            if (ssbg_debug .and. debug2) write(6,*) 'Incrementing itrejcnt for max zmt point for ch.= ', mchan(ji)
            itrejcnt = itrejcnt + 1
            rejcnt(mchan(ji),knosat) = rejcnt(mchan(ji),knosat) + 1
            if (ssbg_debug) then
              if ( itrejcnt <= nch2chk ) then
                write(6,*) '   Channel = ', mchan(ji),            &
                     &     ' Height limit (m) = ', zcheckval,     &
                     &     ' Model height (m) = ', zmt
              endif
            endif
          endif
        endif
      end do

    end do

    if (ssbg_debug .and. (itrejcnt > 0) ) then
      write(6,*) '   Number of topography rejections and observations for this box = ', itrejcnt, knt*knomp
    endif

  end subroutine check_topo


  !--------------------------------------------------------------------------
  ! ssbg_updateObsSpaceAfterInovQc
  !--------------------------------------------------------------------------
  subroutine ssbg_updateObsSpaceAfterInovQc(obsSpaceData, headerIndex, flagsInovQc)

    !:Purpose:      Update obspacedata variables (obstTB and obs flags) after QC
    implicit None

    !Arguments
    type(struct_obs),     intent(inout)  :: obsSpaceData           ! obspaceData Object
    integer,              intent(in)     :: headerIndex            ! current header index
    integer,              intent(in)     :: flagsInovQc(:)         ! data flags

    ! Locals
    integer, allocatable                 :: obsFlags(:)
    integer, allocatable                 :: obsGlobalFlag(:)
    integer, allocatable                 :: satScanPosition(:)

    logical                              :: sensorIndexFound

    integer                              :: actualNumChannel
    integer                              :: bodyIndex
    integer                              :: bodyIndexBeg
    integer                              :: bodyIndexEnd
    integer                              :: channelIndex
    integer                              :: currentChannelNumber
    integer                              :: dataIndex
    integer                              :: headerCompt
    integer                              :: instr
    integer                              :: instrum
    integer                              :: iplatf
    integer                              :: iplatform
    integer                              :: isat
    integer                              :: ntIndex
    integer                              :: numObsToProcess
    integer                              :: sensorIndex


   ! find tvs_sensor index corresponding to current obs

    iplatf      = obs_headElem_i( obsSpaceData, OBS_SAT, headerIndex )
    instr       = obs_headElem_i( obsSpaceData, OBS_INS, headerIndex )

    call tvs_mapSat( iplatf, iplatform, isat )
    call tvs_mapInstrum( instr, instrum )

    sensorIndexFound = .false.
    do sensorIndex =1, tvs_nsensors
      if ( iplatform ==  tvs_platforms(sensorIndex)  .and. &
           isat      ==  tvs_satellites(sensorIndex) .and. &
           instrum   == tvs_instruments(sensorIndex)       ) then
        sensorIndexFound = .true.
        exit
      end if
    end do
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

    obsGlobalFlag(headerCompt) = obs_headElem_i( obsSpaceData, OBS_ST1, headerIndex )
    satScanPosition(headerCompt)  = obs_headElem_i( obsSpaceData, OBS_FOV , headerIndex)
    bodyIndexBeg        = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    bodyIndexEnd = bodyIndexBeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) -1

    BODY2: do bodyIndex =  bodyIndexbeg, bodyIndexEnd
      currentChannelNumber = nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      obsFlags(currentChannelNumber)       = obs_bodyElem_i( obsSpaceData,  OBS_FLG, bodyIndex )
    end do BODY2
    ! Modify flags

    !------------------------------------------------------------------
    ! 1 - Mark global flags if any flagsInovQc > 0
    !------------------------------------------------------------------
    do ntIndex = 1 , numObsToProcess
      if ( maxval(flagsInovQc) > 0 ) obsGlobalFlag(ntIndex) = ibset(obsGlobalFlag(ntIndex), 6)
    end do

    !------------------------------------------------------------------
    ! 2 - Mark obs flags for each value of flagsInovQc
    !------------------------------------------------------------------

    dataIndex = 0
    do ntIndex = 1 , numObsToProcess
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
    do ntIndex = 1 , numObsToProcess
      if (satScanPosition(ntIndex) > 270) satScanPosition(ntIndex) = satScanPosition(ntIndex) - 270
    end do

    ! write elements in obsspace


    call obs_headSet_i(obsSpaceData, OBS_FOV, headerIndex, satScanPosition(1))
    call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex, obsGlobalFlag(1))
    bodyIndexBeg = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    bodyIndexEnd = bodyIndexBeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) -1
    BODY: do bodyIndex =  bodyIndexBeg, bodyIndexEnd
      currentChannelNumber=nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, obsFlags(currentChannelNumber))
    end do BODY

  end subroutine ssbg_updateObsSpaceAfterInovQc


  !--------------------------------------------------------------------------
  ! ssbg_bgCheckSSMIS
  !--------------------------------------------------------------------------
  subroutine ssbg_bgCheckSSMIS( obsSpaceData )

    !:Purpose:        do the quality controle for SSMIS
    implicit None

    !Arguments
    type(struct_obs),     intent(inout)  :: obsSpaceData           ! obspaceData Object

    ! Locals
    logical, allocatable                 :: obsToreject(:)

    integer, allocatable                 :: flagsInovQc(:)
    integer, allocatable                 :: ssmisNewInfoFlag(:)

    integer                              :: codtyp
    integer                              :: dataIndex
    integer                              :: dataIndex1
    integer                              :: headerIndex
    integer                              :: indexFlags
    integer                              :: inovQcSize

    integer                              :: statsInovQcFlags(10)
    real                                 :: percentInovQcFlags(9)

    logical                              :: ssmisDataPresent

    call tmg_start(30,'BGCHECK_SSMIS')
    ssmisDataPresent = .false.
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER0
      codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if ( tvs_isIdBurpInst(codtyp,'ssmis' ) ) then
        ssmisDataPresent = .true.
      end if
    end do HEADER0

    if ( .not. ssmisDataPresent ) then
      write(*,*) 'WARNING: WILL NOT RUN ssbg_bgCheckSSMIS since no SSMIS'
      return
    end if

    statsInovQcFlags(:) = 0
    percentInovQcFlags(:) = 0.0

    write(*,*) ' SSMIS  QC PROGRAM STARTS ....'
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
        write(*,*) 'WARNING: Observation with codtyp = ', codtyp, ' is not SSMIS'
        cycle HEADER
      end if

      !###############################################################################
      ! STEP 1) call satQC SSMIS program                                             !
      !###############################################################################
      if (ssbg_debug) write(*,*) 'STEP 1) call satQC SSMIS program'
      call ssbg_satqcSsmis(obsSpaceData, headerIndex, ssmisNewInfoFlag, obsToreject)

      !###############################################################################
      ! STEP 2) update Flags after satQC SSMIS program                               !
      !###############################################################################
      if (ssbg_debug) write(*,*) 'STEP 2) update Flags after satQC SSMIS program'
      call ssbg_updateObsSpaceAfterSatQc(obsSpaceData, headerIndex, obsToreject)

      !###############################################################################
      ! STEP 3) call inovQC SSMIS program                                            !
      !###############################################################################
      if (ssbg_debug) write(*,*) 'STEP 3) call inovQC SSMIS program'
      call ssbg_inovqcSsmis(obsSpaceData, headerIndex, flagsInovQc)

      !###############################################################################
      ! STEP 4) update Flags after inovQC SSMIS program                              !
      !###############################################################################
      if (ssbg_debug) write(*,*) 'STEP 4) update Flags after inovQC SSMIS program'
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

    call tmg_stop(30)
  end subroutine ssbg_bgCheckSSMIS

end module bgckssmis_mod


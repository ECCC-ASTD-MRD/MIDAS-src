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
  use burp_module
  use MathPhysConstants_mod
  use utilities_mod
  use obsSpaceData_mod
  use obsFiles_mod
  use codePrecision_mod
  use obsFilter_mod
  use tovs_nl_mod
  use gridStateVector_mod
  use timeCoord_mod
  use columnData_mod
  use biasCorrectionSat_mod
  use horizontalCoord_mod
  use analysisGrid_mod
  use obsUtil_mod
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
    if (ierr /= 0) call utl_abort('mwbg_init: Error reading namelist')
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

    real,    intent(in),  dimension(:) :: ztb
    logical, intent(in),  dimension(:) :: grossrej
    real,    intent(out), dimension(:) :: zta

    !  Locals
    integer :: ii, indx1, indx2
    real, dimension(ssbg_maxNumChan) :: spill_coeffs


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
      REAL, DIMENSION(ssbg_maxNumChan), INTENT(IN)  :: SSMIS_Ta
      REAL, DIMENSION(ssbg_maxNumChan), INTENT(OUT) :: Remapped_SSMI_Ta
      INTEGER, INTENT(IN) :: sat_id
      REAL,DIMENSION(ssbg_maxNumChan) :: tbx
      INTEGER(2) :: iChan
      REAL(8), DIMENSION(ssbg_maxNumChan) :: CP

      ! NESDIS Intercept/Slope for F16 SSMIS 7 IMG channels to SSMI linear remapping 
                                 ! ch1  ch2  ch3  ch4  ch5  ch6  ch7  ch8  ch9  ch10 ch11 (LAS/ENV)
      REAL(8), DIMENSION(ssbg_maxNumChan) :: AP=(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                     7.44254,7.80472,6.76383,8.55426,7.34409,6.57813,6.45397, &
                                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

      REAL(8), DIMENSION(ssbg_maxNumChan) :: BP=(/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,&
                                     0.969424,0.967519,0.959808,0.954316,0.958955,0.980339,0.978795, &
                                     1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)

      ! NESDIS F17 and F18 Tb biases with respect to F16
      REAL(8), DIMENSION(ssbg_maxNumChan) :: CP_F17=(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                        -0.779,  -1.446,  -1.013,  -0.522,  -0.240,   0.735,   0.521,   &
                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  
      REAL(8), DIMENSION(ssbg_maxNumChan) :: CP_F18=(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
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

    REAL, DIMENSION(24), INTENT(IN) :: Ta
    REAL, DIMENSION(24), INTENT(OUT) :: Tb

    INTEGER(4) :: iChan
    REAL(4), DIMENSION(24) :: CP, DP

    REAL(8), DIMENSION(24) :: AP=(/0.9850,0.9850,0.9850,0.9850,0.9850,0.9790,0.9815,&
                                   0.9949,0.9934,0.9934,0.9934, &
                                   0.9690,0.9690,0.9740,0.9860,0.9860,0.9880,0.9880,&
                                   0.9815,0.9815,0.9815,0.9815,0.9815,0.9815/)

    REAL(8), DIMENSION(24) :: BP=(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
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
    real, intent(in),  dimension(:) :: Ta
    real, intent(out), dimension(:) :: Tb

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
    REAL, DIMENSION(24), INTENT(IN) :: Tb
    real, intent(out) :: iwv

    !  Locals:
    integer :: ii,precip_screen,ngood

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
    REAL, DIMENSION(24), INTENT(IN) :: Tb
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

    REAL, DIMENSION(24), INTENT(IN) :: Ta
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
    REAL, DIMENSION(24), INTENT(IN) :: Ta, Tb
    INTEGER, INTENT(INOUT) :: Stype
    REAL, INTENT(OUT)   :: CLW
    REAL, INTENT(INOUT) :: WVP
    REAL, INTENT(IN)    :: Lat

    INTEGER(4), PARAMETER :: RT = 285
    INTEGER, PARAMETER :: OCEAN  = 0
    INTEGER, PARAMETER :: SEAICE = 1
    real,    parameter :: clw_limit = 6.0

    REAL :: Ta19v, Ta19h, Ta22v, Ta37v, Ta37h, Ta85v, Ta85h
    REAL :: Tb19v, Tb19h, Tb22v, Tb37v, Tb37h, Tb85v, Tb85h
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


  subroutine cld_filter_fweng(npts,ztb,alg_option,waterobs,ukbadobs,grossrej,  &
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
    integer, intent(in)                  :: npts, isatid
    real,    intent(in), dimension(:)    :: ztb
    real,    intent(inout), dimension(:) :: rclw,riwv
    real,    intent(in),  dimension(:)   :: zlat
    integer, intent(inout)               :: iNewIce

    character(len=6), intent(in) :: alg_option

    logical, intent(in),    dimension(:) :: ukbadobs,grossrej
    logical, intent(inout), dimension(:) :: cloudobs, precipobs
    logical, intent(inout), dimension(:) :: iwvreject,waterobs

    !  Define Parameters
    real, parameter :: iwv_thresh = 80.0    ! Upper bound for IWV in kg/m**2
    integer, parameter :: OCEAN  = 0
    integer, parameter :: SEAICE = 1

    !  Define local variables
    INTEGER :: iChan, iLat, iLon, iError
    integer :: ii,indx1,indx2,Stype
    real :: clw,iwv, Lat

    real, dimension(ssbg_mxval*ssbg_maxObsNum) :: zta

    REAL, DIMENSION(ssbg_maxNumChan) :: F16TDR, RemappedTa, Tb

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
        write(*,*)'CLOUD BY DETERM_CLW = ', clw
        write(*,*)'IWV BY DETERM_CLW = ', iwv
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
    logical                              :: debug

    dataNum = size(globMarq)
    do dataIndex = 1, dataNum
      if (RESETQC) then
        globMarq(dataIndex) = 1024  
      end if
      if ( KCHKPRF(dataIndex) /= 0  ) then
        globMarq(dataIndex) = OR (globMarq(dataIndex),2**6)
      end if
    end do
    if (debug) then
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
      else if ( (ier(i) /= 0  ) .and. (i <= 100 ) ) then 
        print *, ' Input Parameters are not all valid: '
        print *, ' i,tb89(i),tb150(i),pangl(i),ktermer(i) = ',     &
                   i,tb89(i),tb150(i),pangl(i),ktermer(i)
        print *, ' ier(i),scatl(i),scatw(i)=',     &
                   ier(i),scatl(i),scatw(i)
      endif
    end do 

end subroutine bennartz

  !--------------------------------------------------------------------------
  ! mwbg_readGeophysicFieldsAndInterpolate 
  !--------------------------------------------------------------------------
  subroutine mwbg_readGeophysicFieldsAndInterpolate(instName, zlat, zlon, MTINTRP, MGINTRP, GLINTRP)

    implicit none

    !:Purpose: Reads Modele Geophysical variables and save for the first time
    !         TOPOGRAPHIE (MF ou MX):
    !             MF est la topographie filtree avec unites en metres (filtered ME).
    !             MX est la topographie filtree avec unites en m2/s2  (geopotential topography).
    !         Glace de Mer (GL)
    !         Masque Terre-Mer (MG)
    !         Then Interpolate Those variables to observation location
    !Arguments: 
    character(*),       intent(in)   :: instName       ! Instrument Name
    real,               intent(in)   :: zlat(:)        ! Obseravtion Lats
    real,               intent(in)   :: zlon(:)        ! Observation Lons
    real, allocatable,  intent(out)  :: MGINTRP(:)     ! Glace de mer interpolees au pt d'obs.
    real, allocatable,  intent(out)  :: MTINTRP(:)     ! topographie filtree (en metres) et interpolees
    real ,allocatable,  intent(out)  :: GLINTRP(:)     ! Glace de mer interpolees au pt d'obs.
  
    ! Locals:
    real, allocatable, save  :: GL(:)                  ! Modele Glace de Mer (GL)
    real, allocatable, save  :: MG(:)                  ! Modele Masque Terre-Mer (MG)
    real, allocatable, save  :: MT(:)                  ! Modele Topographie (MT)
    real,              save  :: TOPOFACT               ! Facteur x topo pour avoir des unites en metre
    logical,           save  :: ifFirstCall = .True.   ! If .True. we read GL, MT and MG
    integer,           save  ::  gdmt                  ! topo interpolation param
    integer,           save  ::  gdmg                  ! mask terre-mer interpolation param
    integer,           save  ::  gdgl                  ! glace interpolation param
    integer                  ::  gdllsval          
    integer                  :: IUNGEO
    logical                  :: readGlaceMask 
    logical                  :: debug 
    integer                  :: ier, irec 
    integer                  :: ezqkdef, ezsetopt
    integer                  :: FSTINF,FSTPRM,FCLOS
    integer                  :: FSTLIR,FSTFRM, FNOM, FSTOUV
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
    integer, PARAMETER       :: MXLON = 5
    integer, PARAMETER       :: MXLAT = 5
    integer, PARAMETER       :: MXELM = 40
    real,    PARAMETER       :: DLAT = 0.4
    real,    PARAMETER       :: DLON = 0.6
    real                     :: XLAT
    real                     :: XLON
    real, allocatable        :: ZLATBOX (:,:)
    real, allocatable        :: ZLONBOX (:,:)
    real, allocatable        :: MGINTBOX(:,:)
    real, allocatable        :: MTINTBOX(:,:)
    real, allocatable        :: GLINTBOX(:,:)
    integer                  :: dataIndex
    integer                  :: boxPointIndex
    integer                  :: latIndex
    integer                  :: lonIndex
    integer                  :: zlatNum
    integer                  :: zlonNum
    integer                  :: dataNum
    integer                  :: boxPointNum

    ! STEP 0: CHECK if ZLAT AND ZLON ARE SAME DIMENSION
    zlatNum = size(zlat)
    zlonNum = size(zlon)
    if (zlatNum .ne. zlonNum) then
      call utl_abort ('bgckMicrowave_mod: ERREUR: OBSERVATION ZLAT and ZLON should have SAME LENGTH')
    else 
      dataNum = zlatNum
    end if

    ! STEP 1: READ MT, GL and MG from the FST FILE 
    debug = ssbg_debug 
    readGlaceMask = .True.
    if (instName == 'ATMS') readGlaceMask = .False.
    if(ifFirstCall) then
      IUNGEO = 0 
      IER = FNOM(IUNGEO,glmg_file,'STD+RND+R/O',0)

      ! 3) Lecture des champs geophysiques (MF/MX) du modele
      IER = FSTOUV(IUNGEO,'RND')

      ! TOPOGRAPHIE (MF ou MX).
      !     MX est la topographie filtree avec unites en m2/s2  (geopotential topography).

      IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','MX')
      if (IREC .GE. 0) then
        TOPOFACT = 9.80616
        CLNOMVAR = 'MX'
        if(allocated(MT)) deallocate(MT)
        allocate ( MT(NI*NJ), STAT=ier)
        IER = FSTLIR(MT,IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1, &
           ' ',CLNOMVAR)
      else
        call utl_abort ('bgckMicrowave_mod: ERREUR: LA TOPOGRAPHIE (MF or MX) EST INEXISTANTE')
      end if
      
      IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
          IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10,  &
          IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1, &
          IG2, IG3, IG4, IDUM12, IDUM13, IDUM14,  &
          IDUM15, IDUM16, IDUM17, IDUM18 )
       write (*,*) ' GRILLE MT : ',grtyp,ni,nj, &
                ig1,ig2,ig3,ig4
      ier  = ezsetopt('INTERP_DEGREE','LINEAR')  
      ier  = ezsetopt('EXTRAP_DEGREE','ABORT')  
      gdmt = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)

      if (readGlaceMask) then 
        ! MG
        IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','MG')
        if (IREC .LT. 0) then
          call utl_abort ('bgckMicrowave_mod: ERREUR: LE MASQUE TERRE-MER EST INEXISTANT')
        end if

        if(allocated(MG)) deallocate(MG)
        allocate ( MG(NI*NJ), STAT=ier)
        IER = FSTLIR(MG,IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,&
                 ' ','MG')

        IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, &
             IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10, &
             IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1,&
             IG2, IG3, IG4, IDUM12, IDUM13, IDUM14, &
             IDUM15, IDUM16, IDUM17, IDUM18 )
        write (*,*) ' GRILLE MG : ',grtyp,ni,nj, &
                ig1,ig2,ig3,ig4
        ier  = ezsetopt('INTERP_DEGREE','LINEAR')  
        ier  = ezsetopt('EXTRAP_DEGREE','ABORT')  
        gdmg = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
        ! GL
        IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','GL')
        if (IREC .LT. 0) then
          call utl_abort ('bgckMicrowave_mod: ERREUR: LE CHAMP GLACE DE MER EST INEXISTANT')
        end if

        if(allocated(GL)) deallocate(GL)
        allocate ( GL(NI*NJ), STAT=ier)
        IER = FSTLIR(GL,IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1, &
                 ' ','GL')

        IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
             IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10, &
             IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1, &
             IG2, IG3, IG4, IDUM12, IDUM13, IDUM14, &
             IDUM15, IDUM16, IDUM17, IDUM18 )
        write (*,*) ' GRILLE GL : ',grtyp,ni,nj, &
                ig1,ig2,ig3,ig4
        ier  = ezsetopt('INTERP_DEGREE','LINEAR')  
        ier  = ezsetopt('EXTRAP_DEGREE','ABORT')  
        gdgl = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
      else 
        gdgl = -1
        gdmg = -1
      end if 
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
    if(allocated(GLINTBOX)) deallocate(GLINTBOX)
    allocate (GLINTBOX(boxPointNum, dataNum) , STAT=ier) 
    if(allocated(MGINTBOX)) deallocate(MGINTBOX)
    allocate (MGINTBOX(boxPointNum, dataNum) , STAT=ier) 
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
    ier = ezsetopt('INTERP_DEGREE','LINEAR')
    ier = gdllsval(gdmt,mtintbox,mt,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
    if (ier < 0) then
      call utl_abort ('bgckMicrowave_mod: ERROR in the interpolation of MT')
    end if
    if(readGlaceMask) then   
      ier = gdllsval(gdmg,mgintbox,mg,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
      if (ier < 0) then
        call utl_abort ('bgckMicrowave_mod: ERROR in the interpolation of MG')
      end if
      ier = gdllsval(gdgl,glintbox,gl,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
      if (ier < 0) then
        call utl_abort ('bgckMicrowave_mod: ERROR in the interpolation of GL')
      end if
    end if 

    if(allocated(MTINTRP)) deallocate(MTINTRP)
    allocate (MTINTRP(dataNum) , STAT=ier) 
    if(allocated(MGINTRP)) deallocate(MGINTRP)
    allocate (MGINTRP(dataNum) , STAT=ier) 
    if(allocated(GLINTRP)) deallocate(GLINTRP)
    allocate (GLINTRP(dataNum) , STAT=ier) 
    do dataIndex = 1, dataNum
      if (DEBUG) then
        print *, ' ------------------  '
        print *, ' dataIndex = ', dataIndex
        print *, '   '
        print *, ' zlat,zlon = ', zlat(dataIndex), zlon(dataIndex)
        print *, '   '
        print *, ' ZLATBOX = '
        print *,  (ZLATBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        print *, ' ZLONBOX = '
        print *,  (ZLONBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        print *, ' MGINTBOX = '
        print *,  (MGINTBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        print *, ' MTINTBOX = '
        print *,  (MTINTBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        print *, ' GLINTBOX = '
        print *,  (GLINTBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
      end if
      MGINTRP(dataIndex) = 0.0
      MTINTRP(dataIndex) = 0.0
      GLINTRP(dataIndex) = 0.0
      do boxPointIndex=1,MXLAT*MXLON
        MTINTRP(dataIndex) = MAX(MTINTRP(dataIndex),MTINTBOX(boxPointIndex,dataIndex)/TOPOFACT)
        if(readGlaceMask) then      
          MGINTRP(dataIndex) = MAX(MGINTRP(dataIndex),MGINTBOX(boxPointIndex,dataIndex))
          GLINTRP(dataIndex) = MAX(GLINTRP(dataIndex),GLINTBOX(boxPointIndex,dataIndex))
        end if
      end do
      if (DEBUG) then
        print *, ' MGINTRP = ', MGINTRP(dataIndex)
        print *, ' MTINTRP = ', MTINTRP(dataIndex)
        print *, ' GLINTRP = ', GLINTRP(dataIndex)
      end if
    end do
  end subroutine mwbg_readGeophysicFieldsAndInterpolate


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
  
    integer, intent(in)  :: npts
    real,    intent(in),  dimension(:) :: zlat,zlon
    integer, intent(inout),  dimension(:) :: zlq
    integer,allocatable,  intent(out) :: ztt(:)
  
    logical, allocatable, intent(out) :: waterobs(:)

  !  Locals:
    integer, parameter :: mxlat=5,mxlon=5
    integer, parameter :: iungeo=50
    logical, save      :: firstCall=.true.

    integer :: ier,key
    integer :: ni,nj,nk,nilg,njlg
    integer :: ig1,ig2,ig3,ig4,ig1lg,ig2lg,ig3lg,ig4lg
    integer :: idum1,idum2,idum3
    integer :: idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11
    integer :: idum12,idum13,idum14,idum15,idum16,idum17,idum18
  
    integer :: indx,ii,jj,kk
    integer :: blat,blon,imask,glace
    integer :: nlat,nlon
  
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
    integer, external :: fstouv,fstfrm,fstinl,fstvoi

    ! Allocate space for arrays holding values on mesh grid pts.
    call utl_reAllocate(latmesh, mxlat*mxlon)
    call utl_reAllocate(lonmesh, mxlat*mxlon)
    call utl_reAllocate(mgintob, mxlat*mxlon)
    call utl_reAllocate(lgintob, mxlat*mxlon)
    call utl_reAllocate(zlatbox, mxlat*mxlon, npts)
    call utl_reAllocate(zlonbox, mxlat*mxlon, npts)
    call utl_reAllocate(ztt, npts)
    call utl_reAllocate(waterobs, npts)

    ! Open FST file.

    ier = fnom( iungeo,glmg_file,'STD+RND+R/O',0 )
    ier = fstouv( iungeo,'RND' )

    if (firstCall) then
      firstCall = .false.
      ! Read MG field.
      key = fstinf(iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ' ,'MG')
      if ( key <  0 ) then
        call utl_abort('mwbg_landIceMaskAtms: The MG field is MISSING')
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

    ier = fstfrm(iungeo)
    ier = fclos(iungeo)

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
    integer, external :: fstouv,fstfrm,fstinl,fstvoi



    ! Open Wentz surface field if first call

    ier = fnom( iunin,wentz_file,'STD+RND',0 )
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
        ssmisDataPresent = .true.
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
  !  mwbg_readObsFromObsSpace
  !--------------------------------------------------------------------------
  subroutine mwbg_readObsFromObsSpace(instName, headerIndex, &
                                      satIdentifier, satZenithAngle, landQualifierIndice, &
                                      terrainTypeIndice, obsLatitude, obsLongitude, &
                                      satScanPosition, obsQcFlag1, satOrbit, & 
                                      obsGlobalMarker, burpFileSatId, obsTb, btClear, &
                                      obsTbBiasCorr, ompTb, obsQcFlag2, obsChannels, &
                                      obsFlags, sensorIndex, actualNumChannel, obsSpaceData)
    
    !:Purpose:        copy headers and bodies from obsSpaceData object to arrays

    implicit None

    !Arguments
    character(len=9),     intent(in)     :: InstName               ! Instrument Name
    integer,              intent(in)     :: headerIndex            ! current header Index 
    integer, allocatable, intent(out)    :: satIdentifier(:)       ! satellite identifier
    real   , allocatable, intent(out)    :: satZenithAngle(:)      ! satellite zenith angle (btyp=3072,ele=7024) 
    integer, allocatable, intent(out)    :: landQualifierIndice(:) ! land/sea qualifier     (btyp=3072,ele=8012)
    integer, allocatable, intent(out)    :: terrainTypeIndice(:)   ! terrain-type (ice)     (btyp=3072,ele=13039)
    real   , allocatable, intent(out)    :: obsLatitude(:)         ! latitude values (btyp=5120,ele=5002)
    real   , allocatable, intent(out)    :: obsLongitude(:)        ! longitude values (btyp=5120,ele=6002)
    integer, allocatable, intent(out)    :: satScanPosition(:)     ! scan position (fov)    (btyp=3072,ele=5043)
    integer, allocatable, intent(out)    :: obsQcFlag1(:,:)        ! flag values for btyp=3072 block ele 033078, 033079, 033080
    integer, allocatable, intent(out)    :: satOrbit(:)            ! orbit number
    integer, allocatable, intent(out)    :: obsGlobalMarker(:)     ! global Marqueur Data
    character(*),intent(out)             :: burpFileSatId          ! Platform Name
    real   , allocatable, intent(out)    :: obsTb(:)               ! brightness temperature (btyp=9248/9264,ele=12163) 
    real   , allocatable, intent(out)    :: btClear(:)             ! clear brightness temperature (btyp=9248/9264,ele=btClearElementId)
    real   , allocatable, intent(out)    :: obsTbBiasCorr(:)       ! bias correction 
    real   , allocatable, intent(out)    :: ompTb(:)               ! OMP values
    integer, allocatable, intent(out)    :: obsQcFlag2(:)          ! flag values for btyp=9248 block ele 033081      
    integer, allocatable, intent(out)    :: obsChannels(:)         ! channel numbers btyp=9248 block ele 5042 (= 1-22)
    integer, allocatable, intent(out)    :: obsFlags(:)            ! data flags
    integer,              intent(out)    :: sensorIndex            ! find tvs_sensor index corresponding to current obs
    integer,              intent(out)    :: actualNumChannel       ! actual Num channel

    type(struct_obs),     intent(inout)  :: obsSpaceData           ! obspaceData Object

    ! Locals
    integer                              :: bodyIndex
    integer                              :: obsNumCurrentLoc
    integer                              :: bodyIndexbeg
    integer                              :: headerCompt 
    integer                              :: currentChannelNumber  
    integer                              :: channelIndex
    integer                              :: numObsToProcess       
    integer                              :: iplatform
    integer                              :: instrum
    integer                              :: isat, iplatf
    integer                              :: instr
    logical                              :: sensorIndexFound

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
    if ( .not. sensorIndexFound ) call utl_abort('mwbg_readObsFromObsSpace: sensor Index not found') 

    ! find actual Number of channels
    actualNumChannel = tvs_coefs(sensorIndex)%coef%fmv_ori_nchn

    headerCompt = 1 
    numObsToProcess = 1
    ! Allocate Header elements
    call utl_reAllocate(satIdentifier, numObsToProcess)
    call utl_reAllocate(satZenithAngle, numObsToProcess)
    call utl_reAllocate(landQualifierIndice, numObsToProcess)
    call utl_reAllocate(terrainTypeIndice, numObsToProcess)
    call utl_reAllocate(obsLatitude, numObsToProcess)
    call utl_reAllocate(obsLongitude, numObsToProcess)
    call utl_reAllocate(satScanPosition, numObsToProcess)
    call utl_reAllocate(obsGlobalMarker, numObsToProcess)
    call utl_reAllocate(satOrbit, numObsToProcess)
    call utl_reAllocate(obsQcFlag1, numObsToProcess,3)
    ! Allocate Body elements
    call utl_reAllocate(obsTb, numObsToProcess*actualNumChannel)
    call utl_reAllocate(btClear, numObsToProcess*actualNumChannel)
    call utl_reAllocate(ompTb, numObsToProcess*actualNumChannel)
    call utl_reAllocate(obsTbBiasCorr, numObsToProcess*actualNumChannel)
    call utl_reAllocate(obsFlags, numObsToProcess*actualNumChannel)
    call utl_reAllocate(obsChannels, numObsToProcess*actualNumChannel)
    call utl_reAllocate(obsQcFlag2, numObsToProcess*actualNumChannel)
    !initialization
    obsTb(:) = ssbg_realMissing
    btClear(:) = ssbg_realMissing
    ompTb(:) = ssbg_realMissing
    obsTbBiasCorr(:) = ssbg_realMissing

        
    burpFileSatId                      = obs_elem_c    ( obsSpaceData, 'STID' , headerIndex ) 
    satIdentifier(headerCompt)         = obs_headElem_i( obsSpaceData, OBS_SAT, headerIndex ) 
    satZenithAngle(headerCompt)        = obs_headElem_r( obsSpaceData, OBS_SZA, headerIndex ) 
    landQualifierIndice(headerCompt)   = obs_headElem_i( obsSpaceData, OBS_STYP, headerIndex) 
    terrainTypeIndice(headerCompt)     = obs_headElem_i( obsSpaceData, OBS_TTYP, headerIndex) 
    ! If terrain type is missing, set it to -1 for the QC programs
    if (terrainTypeIndice(headerCompt) ==  99) terrainTypeIndice(headerCompt) = -1
    obsLatitude (headerCompt)          = obs_headElem_r( obsSpaceData, OBS_LAT, headerIndex ) 
    obsLongitude(headerCompt)          = obs_headElem_r( obsSpaceData, OBS_LON, headerIndex ) 
    ! Convert lat/lon to degrees
    obsLongitude(headerCompt) = obsLongitude(headerCompt)*MPC_DEGREES_PER_RADIAN_R8
    if( obsLongitude(headerCompt) > 180. ) obsLongitude(headerCompt) = obsLongitude(headerCompt) - 360.
    obsLatitude(headerCompt)  = obsLatitude(headerCompt) *MPC_DEGREES_PER_RADIAN_R8
    satScanPosition(headerCompt)       = obs_headElem_i( obsSpaceData, OBS_FOV , headerIndex) 
    obsGlobalMarker(headerCompt)       = obs_headElem_i( obsSpaceData, OBS_ST1, headerIndex ) 
    satOrbit(headerCompt)              = obs_headElem_i( obsSpaceData, OBS_ORBI, headerIndex) 
    if (instName == 'ATMS') then  
      obsQcFlag1(headerCompt,1)        = obs_headElem_i( obsSpaceData, OBS_AQF1, headerIndex) 
      obsQcFlag1(headerCompt,2)        = obs_headElem_i( obsSpaceData, OBS_AQF2, headerIndex) 
      obsQcFlag1(headerCompt,3)        = obs_headElem_i( obsSpaceData, OBS_AQF3, headerIndex) 
    end if

    bodyIndexbeg        = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    obsNumCurrentLoc    = obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex )

    BODY: do bodyIndex =  bodyIndexbeg, bodyIndexbeg + obsNumCurrentLoc - 1
      currentChannelNumber = nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      obsTb(currentChannelNumber)          = obs_bodyElem_r( obsSpaceData,  OBS_VAR, bodyIndex )
      if ( tvs_mwAllskyAssim ) then
        btClear(currentChannelNumber)      = obs_bodyElem_r( obsSpaceData,  OBS_BTCL, bodyIndex )
      end if
      ompTb(currentChannelNumber)          = obs_bodyElem_r( obsSpaceData,  OBS_OMP, bodyIndex )
      obsTbBiasCorr(currentChannelNumber)  = obs_bodyElem_r( obsSpaceData,  OBS_BCOR,bodyIndex)
      obsFlags(currentChannelNumber)       = obs_bodyElem_i( obsSpaceData,  OBS_FLG, bodyIndex )
      obsQcFlag2(currentChannelNumber)     = obs_bodyElem_i( obsSpaceData,  OBS_QCF2, bodyIndex)
      
    end do BODY
    do channelIndex=1, actualNumChannel
      obsChannels(channelIndex)    = channelIndex+tvs_channelOffset(sensorIndex)
    end do
      

  end subroutine mwbg_readObsFromObsSpace 

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
    integer                              :: obsNumCurrentLoc
    integer                              :: bodyIndexbeg
    integer                              :: headerCompt
    integer                              :: currentChannelNumber
    integer                              :: channelIndex
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
    ssmisDataPresent = .false.
    codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
    if ( tvs_isIdBurpInst(codtyp,'ssmis') ) then
      ssmisDataPresent = .true.
    end if

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
    if ( .not. sensorIndexFound ) call utl_abort('mwbg_readObsFromObsSpace: sensor Index not found')

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
    bodyIndexbeg        = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    obsNumCurrentLoc    = obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex )

    BODY: do bodyIndex =  bodyIndexbeg, bodyIndexbeg + obsNumCurrentLoc - 1
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

    call cld_filter_fweng(numObsToProcess, obsTb, alg_option, waterObs, rainDetectionUKMethod, &
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
    integer                              :: obsNumCurrentLoc
    integer                              :: bodyIndexbeg
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
    if ( .not. sensorIndexFound ) call utl_abort('mwbg_readObsFromObsSpace: sensor Index not found')

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
    bodyIndexbeg        = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    obsNumCurrentLoc    = obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex )

    BODY2: do bodyIndex =  bodyIndexbeg, bodyIndexbeg + obsNumCurrentLoc - 1
      currentChannelNumber = nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      obsFlags(currentChannelNumber)       = obs_bodyElem_i( obsSpaceData,  OBS_FLG, bodyIndex )
    end do BODY2
    ! Modify flags
    
    !------------------------------------------------------------------ 
    ! 1 - Mark global flags if any obsToReject
    !------------------------------------------------------------------ 
    do ntIndex = 1 , numObsToProcess
      !write(*,*) 'OLD GLOBAL FLAG VALUE = ', obsGlobalFlag(ntIndex)
      if ( any(obsToReject) ) obsGlobalFlag(ntIndex) = ibset(obsGlobalFlag(ntIndex), 6)
      !write(*,*) 'NEW GLOBAL FLAG VALUE = ', obsGlobalFlag(ntIndex)
    end do 

    !------------------------------------------------------------------ 
    ! 2 - Mark obs flags for each value of obsToReject
    !------------------------------------------------------------------ 
    dataIndex = 0
    do ntIndex = 1 , numObsToProcess
      do channelIndex = 1, actualNumChannel
        dataIndex = dataIndex+1 
        if (resetQc) obsFlags(dataIndex) = 0
        !write(*,*) 'OLD FLAG VALUE = ', obsFlags(dataIndex)
        if (obsToReject(dataIndex)) obsFlags(dataIndex) = ibset(obsFlags(dataIndex),7)
        !write(*,*) 'NEW FLAG VALUE = ', obsFlags(dataIndex)
      end do 
    end do

    !-----------------------------------------------------------------
    !TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY
    !    Subtract 270 from FOV values (element 005043).
    ! TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY
    !-----------------------------------------------------------------
    do ntIndex = 1 , numObsToProcess
      !write(*,*) 'OLD FOV = ', satScanPosition(ntIndex)
      if (satScanPosition(ntIndex) > 270) satScanPosition(ntIndex) = satScanPosition(ntIndex) - 270
      !write(*,*) 'NEW FOV VALUE = ', satScanPosition(ntIndex)
    end do
 
    ! write elements in obsspace


    call obs_headSet_i(obsSpaceData, OBS_FOV, headerIndex, satScanPosition(1))
    call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex, obsGlobalFlag(1))
    bodyIndexbeg        = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    obsNumCurrentLoc    = obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex )
    BODY: do bodyIndex =  bodyIndexbeg, bodyIndexbeg + obsNumCurrentLoc - 1
      currentChannelNumber=nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, obsFlags(currentChannelNumber))
    end do BODY

  end subroutine ssbg_updateObsSpaceAfterSatQc 
   

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
    integer, allocatable                 :: ssmisNewInfoFlag(:)

    integer                              :: codtyp
    integer                              :: headerIndex
    integer, external                    :: exdb, exfin, fnom, fclos
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
      ! STEP 1) call satQC SSMIS program     !
      !###############################################################################
      write(*,*) 'STEP 1) call satQC SSMIS program'
      call ssbg_satqcSsmis(obsSpaceData, headerIndex , ssmisNewInfoFlag, obsToreject) 

      !###############################################################################
      ! STEP 2) update Flags after satQC SSMIS program     !
      !###############################################################################

      write(*,*)' STEP 2) update Flags after satQC SSMIS program'

      call ssbg_updateObsSpaceAfterSatQc(obsSpaceData, headerIndex,obsToreject) 
    end do HEADER

    call tmg_stop(30)
  end subroutine ssbg_bgCheckSSMIS

end module bgckssmis_mod





































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

program midas_bgckAtms
  !
  ! :Purpose: Main program for background check of microwave instruments. 
  !
  use burp_module
  use bgckmicrowave_mod
  use MathPhysConstants_mod

  !  Object: This program is applied to ATMS derialt data.  The processing applied in this
  !          program includes:

  !              o  Compute model-field (MG,LG) based land-sea qualifier and terrain type for 
  !                 determination of OPEN WATER points
  !              o  OVER OPEN WATER ONLY: Apply cloud/precip filter to identify Tb obs in 
  !                 cloudy/precip regions; add to output BURP file CLW (new element 13209)
  !                 and scattering index SI (new element 13208)
  !              o  OVER LAND/ICE: Compute Dryness Index for filtering AMSU-B like channel data
  !              o  Adds IDENT integer (new element 025174) to output BURP file for each data 
  !                 location; filter results are reflected in IDENT bits 0-11
  !              o  Perform QC checks:
  !                 1) Invalid land/sea qualifier,
  !                 2) Invalid terrain type,
  !                 3) Invalid field of view number,
  !                 4) Satellite zenith angle missing or out of range, (> 75 deg),
  !                 5) land/sea qualifier inconsistent with MG field 
  !                 6) Surface inconsistency (internal lsq,trn not the same as input
  !                    values ilq,itt from BURP file)
  !                        - internal lsq= 0,1 (from MG field (0.0 to 1.0) interpolated to obs points)
  !                        - internal trn=-1,0 (from LG field (0.0 to 1.0) interpolated to obs points
  !                                             and updated with SeaIce retrieved from Tb data)
  !                    Inconsistency affects O-P from surface channels at these locations.
  !                 7) ATMS quality flag check (qual. flag elements 33078,33079,33080,33081)
  !              o  Set data QC flag bit 7 ON for all filtered/rejected data
  !
  !  Second part of bgckAtms does O-P based rejections, topography filtering, channel selection and 
  !  some other checks and sets 3D-Var rejection bits in data QC flags (bits 8,9,11). Data with QC flag
  !  bit 7 ON from this program will have bit 9 set ON.
  !
  !
  !  INPUT:  ATMS derialt file in grid box format.
  !          **NOTE: each grid box is treated individually
  !                            
  !------------------------------------------------------------------
  ! Variable Definitions:
  ! ---------------------
  ! nt           -internal-  number of input obs pts in report
  ! zlat         -internal-  array holding lat values for all obs pts in report
  ! zlon         -internal-  array holding lon values for all obs pts in report
  ! ztb          -internal-  array holding TBs for all obs pts in report
  ! ztbcor       -internal-  array holding bias-corrected TBs for all obs pts in report
  ! scanpos      -internal-  array holding scan positions for all obs pts in report
  ! ilq          -internal-  array holding land/sea qualifier values for all obs of report
  ! itt          -internal-  array holding terrain-type values for all obs pts of report
  ! zenith       -internal-  array holding satellite zenith angles for all obs pts of report
  ! ident        -internal-  array holding 14(12)-bit flag to identify all obs pts in report
  !                          as being over land/ice, cloudy, bad IWV, etc.
  ! waterobs     -internal-  logical array identifying which of the obs in current report are
  !                          over open water away from land and sea-ice
  ! cloudobs     -internal-  logical array identifying which of the obs in current report
  !                          have been identified as being cloud-affected
  ! iwvreject    -internal-  logical array identifying which of the obs in current report
  !                          should be rejected because of bad computed IWV value
  ! precipobs    -internal-  logical array identifying which of the obs in current report
  !                          have been identified as being precipitation-affected
  ! lflagchn     -internal-  logical array identifying which channels for each obs in current report
  !                          are to be excluded from assimilation
  ! rclw         -internal-  real array containing CLW values computed for each obs in report
  ! riwv         -internal-  real array containing ECMWF scattering index computed for each obs in report
  ! iNumSeaIce   -internal-  counter for number of waterobs points converted to sea-ice points (in cld_filter)
  ! iRej         -internal-  counter for number of points rejected in nrl(cld)_filter due bad data
  ! lsq          -internal-  array holding land/sea qualifier values computed from model MG field
  ! trn          -internal-  array holding terrain-type values computed from analysis MG field and SeaIce
  !
  !
  ! Description of input/output files:
  !
  !    burp input  =  DERIALT-type ATMS data BURP file
  !    burp output =  BGCKALT-type ATMS data BURP file containing ATMS Tb [with bits 8, 9 and 11 set for 3dvar reject]
  !                 Data FLAG block QC flag bit 7 set for data rejected for CLW, precip, surface-type, dryness, etc. 
  !                 Also contains 3 new INFO block (btyp=3072) elements that can be plotted with SATPLOT
  !                  - cloud water field CLW (13209) over water,
  !                  - ECMWF scattering index (13208) over water, and
  !                  - ident (info) flag (25174) for QC results, 14 bit (12 bits are used)
  !                 like SSMIS files.
  !
  !    mglg_file =  standard file (FST) containing model MG and LG fields (from GEM analysis/trial)
  !                  - can contain GL (discontinuous ice fraction) instead of LG (continuous ice fraction)
  !                    although discontinuities in GL field may cause interpolation issues.

  !
  !    MT_fst    =  Standard file containing filtered model topo fields MF or MX (for TOPO check)
  !                            NOTE: GEM analysis (_000) and 3h trial (_180m) files contain these fields.
  ! stats_atms_assim  =  ATMS observation error file (new 2013 format)
  !                            NOTE: ERBGCK used for rogue O-P check, UTIL column for channel selection
  !
  !   OPTIONS:
  !    -SPADJUST              = scan bias adjust the Tb data before using the data in algorithms
  !    -IRBC coeff_file_atms  = ATMS bias correction coeff file (for -SPADJUST option)
  !    -MODLSQ                = update values of lsq and trn elements in output BURP file
  !                             (for evaluation purposes only)
  !
  ! Define: (land = land or ice)   water: lsq = 1 and trn = -1  (open water away from coast/ice)
  !                                land:  lsq = 0 or  trn = 0   (over or near land or sea-ice)
  !
  !              CMC option (close to current AMSU-A/B)
  ! Channels never assimilated                                      = 1-4, 16
  ! Channels assimilated over water only                            = 5-6, 17-18
  ! Channels assimilated over land and water (with land topo check) = 7-8, 19-22(+dryness check over land)
  ! Channels assimilated over land and water (no land topo check)   = 9-15
  !
  ! Channels sensitive to clouds/precip                             = 5-9, 17-22 (7-9,19-22  over land and water)
  ! Channels insensitive to clouds/precip  (high peaking)           = 10-15
  !
  !               NRL option
  ! Channels never assimilated                                      = 1-3, 16-17
  ! Channels assimilated over water only                            = 4-6, 18-22  CLW/SI check
  ! Channels assimilated over land and water (with land topo check) = 7-8         no CLW/SI check over land
  ! Channels assimilated over land and water (no land topo check)   = 9-15
  !
  ! Assim channels sensitive to clouds/precip                       = 4-9, 18-22  (7-9 over land and water)
  ! Channels insensitive to clouds/precip  (high peaking)           = 10-15       no filters
  !
  !               ECMWF option
  ! Channels never assimilated                                      = 1-5, 16-17
  ! Channels assimilated over water only                            = 6-8, 18-22  LWP check (ch.6-8,18) [Grody 2001]
  !                                                                        ch.16-17 based SI check (18-22) [Bennartz 2002]
  !                                                                   reject all ch. if ch.3 O-P > 5K
  !                                                   
  ! Channels assimilated over land and water                        = 9-15
  !
  ! Assim channels tested for clouds/precip                         = 6-8, 18-22
  ! Channels insensitive to clouds/precip  (high peaking)           = 9-15
  IMPLICIT NONE


  integer, parameter :: MXSAT = 9
  integer, PARAMETER :: MXVAL = 22
  integer, PARAMETER :: MXNT = 3000
  integer, parameter :: nchanAtms=22
  integer, parameter :: mxscan=96
  real, parameter    :: zmisg=9.9e09
  integer, parameter :: iunbc=30
  integer, parameter :: JPNSAT = 9
  integer, parameter :: JPCH = 50
  integer, parameter :: JPMXREJ = 15
  integer, parameter  :: MXCHN = 42   

  

  ! Upper limit for CLW (kg/m**2) for Tb rejection over water
  real, parameter :: clw_atms_nrl_LTrej=0.175      ! lower trop chans 1-6, 16-20
  real, parameter :: clw_atms_nrl_UTrej=0.2        ! upper trop chans 7-9, 21-22

  ! Other NRL thresholds
  real, parameter :: scatec_atms_nrl_LTrej=9.0     ! lower trop chans 1-6, 16-22
  real, parameter :: scatec_atms_nrl_UTrej=18.0    ! upper trop chans 7-9
  real, parameter :: scatbg_atms_nrl_LTrej=10.0    ! lower trop chans 1-6
  real, parameter :: scatbg_atms_nrl_UTrej=15.0    ! upper trop chans 7-9
  real, parameter :: mean_Tb_183Ghz_min=240.0      ! min. value for Mean(Tb) chans. 18-22 

  ! Highest peaking AMSU-A like ATMS channel for ocean-only assimilation
  !   6 = 700 mb   (AMSU/operations == AMSU chan. 5)
  !   7 = 400 mb   (scat. index used in AMSU/operations -- AMSU chan. 6)
  !   8 = 250 mb   (ECMWF)
  integer, parameter :: ipc=6

  character(len=128)                   :: mglg_file

  ! F90 Allocatable arrays (deferred-shape)

  ! For s/r mwbg_readCoeff

  logical :: sp_adj_tb,  modlsqtt, useUnbiasedObsForClw
  logical :: lutb

  ! External functions

  ! Second part of bgckAtms program (atms_inovqc standalone program)
  !OBJET          Effectuer le controle de qualite des radiances level 1b 
  !               ATMS de NPP.
  !
  !NOTES
  !
  !  Must be run after satbcor and 3D-Var (O-P mode)
  !
  !        5 tests are done:                                                      QC flag bits set
  !          1) check for data rejected by first bgckAtms program (QC flag bit 7 ON)          --> bit 9(,7)
  !          2) topography rejection for low-peaking channels (with MF/MX field), --> bits 9,18
  !          3) check for uncorrected radiance (QC flag bit 6 OFF),               --> bit 11           
  !          4) Innovation (O-P) based QC                                         --> bit 9,16
  !          5) channel blacklisting (from UTIL column in stats_atms_assim file)  --> bit 8
  !
  !       *** Array ITEST(5) in SUBROUTINE mwbg_tovCheckAtms used to select tests. ***

  !  (i)   Basic QC tests plus filtering for surface-sensitive channels, cloud water (CLW), 
  !         scattering index, Dryness Index are done in first bgckAtms program.
  !         QC flag bit 7 is set for the rejected data. Test 1 of this program checks
  !         for such data and sets bit 9 ON (for 3D-Var rejection).
  !
  ! (ii)   This program sets data QC flag bit 9 ON in all "data reject" cases except
  !         for test 12 (blacklisting) where bit 8 is set ON.
  !
  INTEGER MXELM
  PARAMETER ( MXELM  =    30 )

  integer                       :: handle
  integer                       :: i
  integer                       :: j
  integer                       :: error
  integer                       :: nobs_tot
  integer                       :: n_bad_reps
  integer                       :: nb_rpts
  integer                       :: n_reps
  integer                       :: reportIndex
  integer                       :: satNumber 
  integer                       :: ntOut 
  integer                       :: nvalOut
  integer                       :: junk 
  integer                       :: exdb
  integer                       :: exfin
  integer                       :: ier
  integer                       :: ISTAT
  integer                       :: nulnam
  integer                       :: INOSAT
  logical                       :: bad_report
  logical                       :: resumeReport
  logical                       :: ifLastReport
  logical                       :: RESETQC 
  logical                       :: debug
  logical                       :: clwQcThreshold
  logical                       :: allowStateDepSigmaObs  
  integer, parameter            :: MISGINT = -1
  character(len=9)              :: STNID
  character(len=9), allocatable :: satelliteId(:)

  
  integer                       :: nlocs
  type(BURP_FILE)               :: File_in
  type(BURP_FILE)               :: File_out
  type(BURP_RPT)                :: Rpt_in
  type(BURP_RPT)                :: Rpt_out

  integer                       :: nblocs
  integer                       :: nsize
  integer                       :: iun_burpin

  character(len=20)             :: opt_missing
  character(len=9)              :: id

  character(len=90)             :: burpFileNameIn
  character(len=90)             :: burpFileNameOut
  character(len=9)              :: ETIKRESU
  real, allocatable             :: MTINTRP(:)
  real, allocatable             :: GLINTRP(:)
  real, allocatable             :: MGINTRP(:)
 
  real,    allocatable          :: zlat(:)
  real,    allocatable          :: zlon(:)
  integer, allocatable          :: ISAT(:)
  real,    allocatable          :: zenith(:)
  integer, allocatable          :: ilq(:)
  integer, allocatable          :: itt(:)
  real,    allocatable          :: ztb(:)
  real,    allocatable          :: ZOMP(:)
  real,    allocatable          :: biasCorr(:)
  integer, allocatable          :: scanpos(:)
  integer, allocatable          :: qcflag1(:,:)
  integer, allocatable          :: qcflag2(:)
  integer, allocatable          :: ican(:)
  integer, allocatable          :: icanomp(:)
  integer, allocatable          :: IMARQ(:)
  integer, allocatable          :: IORBIT(:)
  integer, allocatable          :: globMarq(:)
  integer, allocatable          :: adresses(:)
  integer                       :: IUTILST(JPCH,JPNSAT)
  real                          :: TOVERRST(JPCH,JPNSAT)
  integer                       :: MREJCOD(JPMXREJ,MXCHN,MXSAT)
  integer                       :: MREJCOD2(JPMXREJ,MXCHN,MXSAT)
  integer, allocatable          :: lsq(:)
  integer, allocatable          :: trn(:)
  integer, allocatable          :: icheck(:,:)
  integer, allocatable          :: ichkprf(:)
  logical, allocatable          :: waterobs(:)
  logical, allocatable          :: grossrej(:)
  logical, allocatable          :: cloudobs(:)
  logical, allocatable          :: iwvreject(:)
  logical, allocatable          :: precipobs(:)
  integer, allocatable          :: ident(:)
  real,    allocatable          :: zdi(:)
  real,    allocatable          :: rclw(:)
  real,    allocatable          :: riwv(:)
  real,    allocatable          :: SeaIce(:)
  real,    allocatable          :: scatec(:)
  real,    allocatable          :: scatbg(:)
  logical, allocatable          :: lflagchn(:,:)
  logical, allocatable          :: lqc(:,:)
  real                          :: glbscanb(nchanAtms,mxscan,mxsat)
  real                          :: dglbscanb(nchanAtms,mxscan,mxsat)
  character(len=9)              :: csatid(mxsat)
  real, allocatable             :: ztbcor(:)
  integer                       :: iNumSeaIce
  integer                       :: iRej
  integer                       :: n_reps_tb2misg

  integer                       :: cldcnt


  integer                       :: ierr
  integer                       :: ilnmx
  integer                       :: status
  integer                       :: nele
  integer                       :: nt
  integer                       :: blat
  integer                       :: blon
  integer                       :: idtyp
  integer                       :: flgcnt
  integer                       :: rejcnt
  integer                       :: landcnt
  integer                       :: iwvcnt
  integer                       :: pcpcnt
  integer                       :: drycnt
  integer                       :: indx1, indx2, ii, numbsat, iich, jj, kk


  namelist /nambgck/ debug, sp_adj_tb, modlsqtt, useUnbiasedObsForClw, RESETQC, ETIKRESU 

  burpFileNameIn = './obsatms'
  burpFileNameOut = './obsatms.out'

  mglg_file  = './fstmglg'

  istat = exdb('midas_bgckAtms','DEBUT','NON')

  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM MIDAS-BGCKATMS: --",/,' //   &
            '14x,"-- BACKGROUND CHECK FOR ATMS OBSERVATIONS --",/, ' //&
            '14x,"-- Revision : ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  debug = .false.
  sp_adj_tb = .false.
  modlsqtt = .false.
  useUnbiasedObsForClw = .false.
  RESETQC = .FALSE.
  ETIKRESU = '>>BGCKALT'

  ! Initializations of counters (for total reports/locations in the file).
  flgcnt = 0
  landcnt = 0
  rejcnt = 0
  cldcnt = 0
  iwvcnt = 0
  pcpcnt = 0
  drycnt = 0
  iNumSeaIce = 0
  iRej = 0

  ! read namelist variables
  ! reading namelist
  nulnam = 0
  ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam, nml=nambgck, iostat=ierr)
  if ( ierr /= 0 ) then
    write(*,*) 'midas_bgckAtms: Error reading namelist'
    call abort()
  end if
  write(*,nml=nambgck)
  ierr = fclos(nulnam)

  mwbg_debug = debug
  mwbg_modlsqtt = modlsqtt
  mwbg_useUnbiasedObsForClw = useUnbiasedObsForClw 

   
  if ( modlsqtt ) then
    write(*,*) 'MODLSQ option is activated!'
    write(*,*) '  Output file will contain recomputed values for land/sea qualifier&
                 and terrain type based on LG/MG.'
  endif

  !#################################################################################
  ! STEP 1 ) Lecture des statistiques d'erreur totale pour les  TOVS 
  !#################################################################################

  call mwbg_readStatTovs('./stats_atms_assim', 'ATMS', satelliteId, IUTILST, TOVERRST)
  satNumber = size(satelliteId)

  ! MAIN LOOP through all the reports in the BURP file
  n_bad_reps = 0
  reportIndex = 0
  nobs_tot = 0 
  n_reps_tb2misg = 0 
  REPORTS: do 

    reportIndex = reportIndex + 1

    !###############################################################################
    ! STEP 2 ) Lecture des observations TOVS dans le fichier burpFileNameIn
    !###############################################################################
    write(*,*) ' ==> mwbg_getData: '
    call mwbg_getData(burpFileNameIn, reportIndex, ISAT, zenith, ilq, itt, &
                      zlat, zlon, ztb, biasCorr, ZOMP, scanpos, nvalOut, &
                      ntOut, qcflag1, qcflag2, ican, icanomp, IMARQ, IORBIT, &
                      globMarq, resumeReport, ifLastReport, 'ATMS', STNID)

    if ( .not. resumeReport ) then
      if ( ALL(ZOMP(:) == MPC_missingValue_R4 )) then 
        n_bad_reps = n_bad_reps + 1  
        write(*,*) 'Bad Report '
        cycle REPORTS
      end if
      ! Increment total number of obs pts read
      nobs_tot = nobs_tot + ntOut

    !###############################################################################
    ! STEP 3 ) Determine which obs pts are over open water (i.e NOT near coasts or
    !          over/near land/ice) using model MG and LG fields from glbhyb2 ANAL    
    !###############################################################################
      write(*,*) ' ==> mwbg_landIceMaskAtms: '
      call mwbg_landIceMaskAtms(mglg_file,ntOut,zlat,zlon,ilq,itt,lsq,trn,waterobs)

    !###############################################################################
    ! STEP 4 ) Check for values of TB that are missing or outside physical limits.    
    !###############################################################################

      write(*,*) ' ==> mwbg_grossValueCheck: '
      call mwbg_grossValueCheck(ntOut,ztb,grossrej)
      
      if ( ANY(grossrej) ) then
        write(*,*) ' mwbg_grossValueCheck has detected bad Tb data. Number of affected &
	            locations = ', COUNT(grossrej)
        write(*,*) '   Box lat, lon = ', blat, blon
      endif

    !###############################################################################
    ! STEP 5 ) Preliminary QC checks --> set lqc(nt,nchanAtms)=.true. 
    !          for data that fail QC     
    !###############################################################################

      write(*,*) ' ==> mwbg_firstQcCheckAtms: '
      call mwbg_firstQcCheckAtms(zenith, ilq, itt, zlat, zlon, ztb, scanpos, stnid, &
                                 nvalOut, ntOut, lqc, grossrej, lsq, trn, qcflag1, &
                                 qcflag2, ican, blat, blon, lutb)

      if ( lutb ) n_reps_tb2misg = n_reps_tb2misg + 1

      bad_report = .false.
      if ( COUNT(lqc) == ntOut*nvalOut ) then
        write(*,*) ' mwbg_firstQcCheckAtms has detected a problem with data in this &
                     report!'
        write(*,*) '   Report box lat, lon = ', blat, blon
        bad_report = .true.
        n_bad_reps = n_bad_reps + 1
        cycle REPORTS
      endif
      !  Exclude problem points from further calculations
      do kk = 1,ntOut
        if ( COUNT(lqc(kk,:)) == nchanAtms ) grossrej(kk) = .true.
      end do

    !###############################################################################
    ! STEP 6 )remove scan-dependency of Tb biases for data in the report (OPTION)
    !         or keep initial ztb if no correction    
    !###############################################################################

      if (sp_adj_tb) then    
        call mwbg_correctTbScanDependency(ntOut, numbsat, csatid, stnid, dglbscanb, &
                                          scanpos, ztb, ztbcor)           
      else  ! no correction
        if(allocated(ztbcor)) deallocate(ztbcor)
        allocate(ztbcor(ntOut*nvalOut), stat = ier)
        ztbcor(:) = ztb(:)
      end if

    !###############################################################################
    ! STEP 7 ) mwbg_nrlFilterAtms returns rclw, scatec, scatbg and also does sea-ice 
    !          detection missing value for  rclw, scatec, scatbg  is -99.0 (e.g. over
    !          land or sea-ice).Sets trn=0 (sea ice) for points where retrieved SeaIce
    !          >=0.55. Does nothing if trn=0 (sea ice) and retrieved SeaIce<0.55.
    !###############################################################################
 
      !  
      write(*,*) ' ==> mwbg_nrlFilterAtms: '
      call mwbg_nrlFilterAtms(ntOut,ztbcor,biasCorr, zenith, zlat, lsq, trn, &
                              waterobs, grossrej, rclw, scatec, scatbg, iNumSeaIce, iRej, &
                              SeaIce)

      !###############################################################################
      ! STEP 8 ) Apply NRL cloud filter, scattering index and sea-ice detection algorithms
      !          to OPEN WATER (waterobs=true) points.
      ! Points with SeaIce>0.55 are set to sea-ice points (waterobs --> false)
      !###############################################################################
    
      write(*,*) ' ==> mwbg_flagDataUsingNrlCriteria: '
      call mwbg_flagDataUsingNrlCriteria(ntOut,ztbcor, biasCorr, rclw, scatec, scatbg, &
                                         SeaIce, grossrej, waterobs, useUnbiasedObsForClw, &
                                         scatec_atms_nrl_LTrej, scatec_atms_nrl_UTrej, &
                                         scatbg_atms_nrl_LTrej, scatbg_atms_nrl_UTrej, &
                                         clw_atms_nrl_LTrej,  clw_atms_nrl_UTrej, &
                                         mean_Tb_183Ghz_min, zmisg, iwvreject, cloudobs, &
                                         precipobs, cldcnt , ident, riwv, zdi)




      !###############################################################################
      ! STEP 7 ) ! Review all the checks previously made to determine which obs are to be 
      !            accepted for assimilation and which are to be flagged for exclusion 
      !            (lflagchn). 
      !            grossrej()  = .true. if any channel had a gross error at the point
      !            cloudobs()  = .true. if CLW > clw_atms_nrl_LTrej (0.175) or precipobs
      !            precipobs() = .true. if precip. detected through NRL scattering indices
      !            waterobs()  = .true. if open water point
      !            iwvreject() = .true. if Mean 183 Ghz [ch. 18-22] Tb < 240K (too dry 
      !            for ch.20-22 over land)
      !###############################################################################

      write(*,*) ' ==> mwbg_reviewAllcriteriaforFinalFlags: '
      call mwbg_reviewAllcriteriaforFinalFlags(ntOut,nvalOut, ipc, lqc, grossrej, waterobs, &
                                               precipobs, rclw, scatec, scatbg, iwvreject, riwv, &
                                               IMARQ, globMarq, clw_atms_nrl_LTrej, clw_atms_nrl_UTrej,&
                                               scatec_atms_nrl_LTrej, scatec_atms_nrl_UTrej, &
                                               scatbg_atms_nrl_LTrej, scatbg_atms_nrl_UTrej, zdi, &
                                               ident, zmisg, lflagchn, drycnt, landcnt, rejcnt, &
                                               iwvcnt, pcpcnt, flgcnt)


      !** start second quality control (atms_inovqc standalone program) **
      ! trouver l'indice du satellite

      INOSAT = 0
      DO I = 1, satNumber
        if ( STNID .EQ. '^'//satelliteId(I) ) THEN
          INOSAT = I
          write(*,*)' SATELLITE = ', STNID
          write(*,*)'    INOSAT = ', INOSAT
        end if
      ENDDO
      if ( INOSAT .EQ. 0 ) THEN
        write(*,*) 'SATELLITE ',TRIM(STNID), &
                   ' NOT FOUND IN STATS FILE!'
        call ABORT()
      end if


      !###############################################################################
      ! STEP 8 ) Interpolation de le champ MF/MX (topogrpahy) aux pts TOVS.
      !    N.B.: on examine ce champ sur une boite centree sur chaque obs.
 
      !###############################################################################
      write(*,*) ' ==> mwbg_readGeophysicFieldsAndInterpolate: '
      call mwbg_readGeophysicFieldsAndInterpolate('ATMS', zlat, zlon, MTINTRP, MGINTRP, &
                                                  GLINTRP)

      !###############################################################################
      ! STEP 9 ) Controle de qualite des TOVS. Data QC flags (IMARQ) are modified here!
      !###############################################################################
      write(*,*) ' ==> mwbg_tovCheckAtms: '
      call mwbg_tovCheckAtms(TOVERRST, IUTILST, ISAT, IORBIT, ican, ICANOMP, ztb, biasCorr, &
                             ZOMP, ICHECK, nvalOut, nvalOut, ntOut, ZMISG, INOSAT, ident, &
                             ICHKPRF, scanpos, MTINTRP, IMARQ, MREJCOD, MREJCOD2, STNID, &
                             RESETQC)

      !###############################################################################
      ! STEP 10 ) Accumuler Les statistiques sur les rejets
      !###############################################################################
      write(*,*) ' ==> mwbg_qcStats: '
      call mwbg_qcStats('ATMS', satNumber, ICHECK, ican, INOSAT, nvalOut, &
                        ntOut, satelliteId, .false., MREJCOD, MREJCOD2)
    end if

    !###############################################################################
    ! STEP 11 ) Update the burpfile out burpFileNameIn
    !###############################################################################

    write(*,*) ' ==> mwbg_updateBurpAtms: '
    call mwbg_updateBurpAtms(burpFileNameIn, reportIndex, ETIKRESU, ztb, lsq, trn, riwv, &
                             rclw, ident, globMarq, lflagchn, RESETQC, ICHKPRF,IMARQ, lutb,&
                             burpFileNameout)

    if ( ifLastReport) exit REPORTS

  end do REPORTS

  write(*,*) ' --------------------------------------------------------------- '
  write(*,*) ' Number of obs pts read from BURP file              = ', nobs_tot
  write(*,*) ' Number of BURP file reports                        = ', reportIndex
  write(*,*) ' Number of bad BURP file reports (all data flagged) = ', n_bad_reps
  write(*,*) ' Number of BURP file reports where Tb set to zmisg  = ', n_reps_tb2misg
  write(*,*) ' --------------------------------------------------------------- '
  write(*,*) ' 1. Number of obs pts found over land/ice           = ', landcnt
  write(*,*) ' 2. Number of problem obs pts (Tb err, QCfail)      = ', rejcnt
  write(*,*) ' 3. Number of cloudy obs  (CLW > clw_min)           = ', cldcnt
  write(*,*) ' 4. Number of scatter/precip obs                    = ', pcpcnt
  write(*,*) ' 5. Number of pts with Mean 183 Ghz Tb < 240K       = ', iwvcnt
  write(*,*) ' 6. Number of pts flagged for AMSU-B Dryness Index  = ', drycnt
  write(*,*) ' --------------------------------------------------------------- '
  write(*,*) ' Total number of filtered obs pts                   = ', flgcnt
  write(*,*) ' ----------------------------------------------------------------'
  write(*,*) ' '
  write(*,*) ' -------------------------------------------------------------------------------'
  write(*,*) ' Number of waterobs points converted to sea ice points         = ', iNumSeaIce
  write(*,*) ' Number of points where CLW/SI missing over water due bad data = ', iRej
  write(*,*) ' -------------------------------------------------------------------------------'
  write(*,*) ' '
  write(*,*) ' '
  write(*,*) '   Meaning of IDENT flag bits: '
  write(*,*) ' '
  write(*,*) '      BIT    Meaning'
  write(*,*) '       0     off=land or sea-ice, on=open water away from coast'
  write(*,*) '       1     Mean 183 Ghz [ch. 18-22] is missing'
  write(*,*) '       2     NRL CLW is missing (over water)'
  write(*,*) '       3     NRL > clw_atms_nrl_LTrej (0.175 kg/m2) (cloudobs)'
  write(*,*) '       4     scatec/scatbg > Lower Troposphere limit 9/10 (precipobs)'
  write(*,*) '       5     Mean 183 Ghz [ch. 18-22] Tb < 240K'
  write(*,*) '       6     CLW > clw_atms_nrl_UTrej (0.200 kg/m2)'
  write(*,*) '       7     Dryness Index rejection (for ch. 22)'
  write(*,*) '       8     scatec/scatbg > Upper Troposphere limit 18/15'
  write(*,*) '       9     Dryness Index rejection (for ch. 21)'
  write(*,*) '      10     Sea ice > 0.55 detected'
  write(*,*) '      11     Gross error in Tb (any chan.) or other QC problem (all channels rejected)'
  write(*,*) ' '
  write(*,*) '   New Element 13209 in BURP file = CLW (kg/m2)'
  write(*,*) '   New Element 13208 in BURP file = ECMWF Scattering Index'
  write(*,*) '   New Element 25174 in BURP file = IDENT flag'
  write(*,*) ' '

  !###############################################################################
  ! STEP 12 ) Print the statistics in listing file 
  !###############################################################################
  call mwbg_qcStats('ATMS', satNumber, ICHECK, ican, INOSAT, nvalOut, &
                        ntOut, satelliteId, .TRUE., MREJCOD, MREJCOD2)

  istat = exfin('midas_bgckAtms','FIN','NON')

end program midas_bgckAtms

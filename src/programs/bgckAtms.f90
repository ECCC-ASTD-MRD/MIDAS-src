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

  integer  :: error,nb_rpts,ref_rpt,compteur,nbele,nlocs
  integer  :: handle,nvale,nte,j,i,btyp,bfam

  type(BURP_FILE)        :: File_in,File_out
  type(BURP_RPT)         :: Rpt_in,Rpt_out

  integer,allocatable    :: adresses(:)
  integer                :: ibrptime, nblocs, nsize, iun_burpin, irun

  character(len=20)      :: opt_missing
  character(len=90)      :: brp_in,brp_out
  character(len=9)       :: id

  integer, parameter :: MXSAT = 9
  integer, PARAMETER :: MXVAL = 22
  integer, PARAMETER :: MXNT = 3000
  integer, parameter :: nchanAtms=22
  integer, parameter :: mxscan=96
  real, parameter    :: zmisg=9.9e09
  integer, parameter :: iunbc=30

  integer :: ierr,ilnmx,status
  integer :: istat,nombre
  integer :: nele,nt,blat,blon, idtyp
  integer :: flgcnt,rejcnt
  integer :: cldcnt,landcnt,iwvcnt,pcpcnt
  integer :: drycnt
  integer :: indx1, indx2, ii, numbsat, iich, jj, kk
  integer :: iNumSeaIce, nobs_tot, n_cld, iRej, n_bad_reps, n_reps, n_reps_tb2misg

  integer, dimension(31)               :: alloc_status

  ! These arrays are for subroutine GETBC, called to get bias corrections
  integer, dimension(mxsat)            ::  numor, numan
  integer, dimension(nchanAtms,mxsat)       ::  listan

  real, dimension(nchanAtms,mxscan,mxsat)   ::  glbscanb, dglbscanb

  character(len=9), dimension(mxsat)   ::  csatid

  real, dimension(nchanAtms,mxscan)    ::  zbcor
  real, dimension(mxval*mxnt)          ::  ztbcor
  
  real, dimension(5)                   ::  ztb183

  integer, dimension(mxnt)             :: scanpos
  real, dimension(mxval*mxnt)          :: ztb, biasCorr
  real, dimension(mxnt)                :: zlat,zlon,zenith
  integer, dimension(mxnt)             :: ilq,itt
  integer, dimension(mxval*mxnt)       :: ican, qcflag2
  integer, dimension(mxnt,3)           :: qcflag1
  integer :: reportIndex

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

  character(len=9)                     :: stnid
  character(len=128)                   :: mglg_file

  ! F90 Allocatable arrays (deferred-shape)
  integer, allocatable, dimension(:)   :: ident,iber,err,lsq,trn
  real,    allocatable, dimension(:)   :: rclw,riwv,ascatw,SeaIce
  real,    allocatable, dimension(:)   :: ztb89,ztb150,scatl,scatw,tb23,tb31,tb50,tb53,tb89
  real,    allocatable, dimension(:)   :: bcor23, bcor31, bcor50, bcor89, bcor150
  real,    allocatable, dimension(:)   :: ztb_amsub3, ztb_amsub5, zdi, scatec, scatbg
  real,    allocatable, dimension(:)   :: bcor_amsub3, bcor_amsub5
  logical, allocatable, dimension(:)   :: waterobs,ukbadobs
  logical, allocatable, dimension(:)   :: cloudobs,iwvreject,precipobs,grossrej
  logical, allocatable, dimension(:,:) :: lflagchn, lqc

  ! For s/r mwbg_readCoeff
  INTEGER                                         :: NSAT, NFOV
  INTEGER, PARAMETER                              :: maxpred=6
  character(len=90)                               :: coef_in
  character(len=9), dimension(mxsat)              :: SATNAMES
  integer, dimension(mxsat)                       :: RCNCHAN
  integer, dimension(mxsat,nchanAtms)                 :: CHANNUMS, NPRED
  real,dimension(mxsat,nchanAtms,maxpred+1)           :: COEFF
  real,dimension(mxsat,nchanAtms,mxscan)              :: FOVBIAS
  character(len=5)                                :: CF_INSTRUM
  character(len=2),dimension(mxsat,nchanAtms,maxpred) :: PTYPES

  logical :: sp_adj_tb, debug, modlsqtt, useUnbiasedObsForClw
  logical :: bad_report, lutb, resume_report, qcflag1_check

  ! External functions
  integer, external :: exdb,exfin

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
  INTEGER MXLAT, MXLON
  PARAMETER ( MXELM  =    30 )
  PARAMETER ( MXLAT  =     5 )
  PARAMETER ( MXLON  =     5 )

  integer ezsetopt, ezqkdef
  integer gdllsval, gdmg

  integer :: FSTOUV
  INTEGER FSTINF,FSTPRM,FSTLIR,FSTFRM
  INTEGER NVAL, nvalOut , ntOut
  INTEGER IER,IREC,IREC2,JUNK
  INTEGER JN, JL
  INTEGER IUNGEO, IUNSTAT, INUMSAT
  INTEGER INOSAT
  INTEGER IDUM,IDUM1,IDUM2,IDUM3,IDUM4,IDUM5,IDUM6,IDUM7
  INTEGER IDUM8,IDUM9,IDUM10,IDUM11,IDUM12,IDUM13
  INTEGER IDUM14,IDUM15,IDUM16,IDUM17,IDUM18
  INTEGER IG1,IG2,IG3,IG4
  INTEGER IG1R,IG2R,IG3R,IG4R
  INTEGER NI,NJ,NK,INDX,NLAT,NLON

  INTEGER TBLVAL    (MXELM*MXVAL*MXNT)
  INTEGER KTBLVALN  (MXELM*MXVAL*MXNT)
  INTEGER LSTELE    (MXELM)
  INTEGER KLISTEN   (MXELM)
  INTEGER ELDALT    (MXELM)
  INTEGER IDATA     (MXVAL*MXNT)
  INTEGER ISCNCNT   (MXNT)
  INTEGER ISAT      (MXNT)
  INTEGER IORBIT    (MXNT)
  INTEGER ICANOMP   (MXVAL*MXNT)
  INTEGER ICHECK    (MXVAL*MXNT)
  INTEGER ICHKPRF   (MXNT)
  INTEGER IMARQ     (MXVAL*MXNT)

  CHARACTER *12  ETIKXX
  CHARACTER *9   ETIKRESU
  CHARACTER *4   NOMVXX,CLNOMVAR
  CHARACTER *2   TYPXX 
  CHARACTER *1   GRTYP

  INTEGER, ALLOCATABLE, DIMENSION(:) :: BUF1
  REAL, ALLOCATABLE, DIMENSION(:) :: MT

  REAL  DONIALT (MXELM*MXVAL*MXNT)
  REAL  PRVALN  (MXELM*MXVAL*MXNT)
  REAL  ZDATA   (MXVAL*MXNT)
  REAL  MTINTRP (MXNT)
  REAL  ZOMP    (MXVAL*MXNT)
  REAL  ZOMPNA  (MXVAL*MXNT)
  REAL  ZLATBOX (MXLAT*MXLON,MXNT)
  REAL  ZLONBOX (MXLAT*MXLON,MXNT)
  REAL  MTINTBOX(MXLAT*MXLON,MXNT)
  REAL  XLAT,XLON

  REAL  DLAT, DLON, TOPOFACT

  LOGICAL RESETQC, SKIPENR

  DATA IUNGEO  / 55 /
  DATA IUNSTAT / 60 /
  DATA DLAT   / 0.4 /
  DATA DLON   / 0.6 /

  integer :: nulnam
  namelist /nambgck/ debug, sp_adj_tb, modlsqtt, useUnbiasedObsForClw, RESETQC, ETIKRESU 

  brp_in = './obsatms'
  brp_out = './obsatms.out'

  mglg_file  = './fstmglg'
  coef_in = './bcor'

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

  ! Optional adjustment of Tb for scan-position dependency before computing CLW.
  ! Open and read new UBCOR system ascii file containing mean O-P for each scan position (FOVBIAS).
  if ( sp_adj_tb ) then
    !call mwbg_readCoeff(SATNAMES,CHANNUMS,FOVBIAS,COEFF,NSAT,RCNCHAN,NFOV,NPRED,CF_INSTRUM,maxpred,iunbc,coef_in,PTYPES)
    !if ( NSAT == 0 ) then
    !  write(*,*) 'ERROR: Cannot apply scan position bias correction -- error reading bcor file ', coef_in
    !  call abort()
    !endif
    !if ( TRIM(CF_INSTRUM) /= 'ATMS' ) then
    !  write(*,*) 'ERROR: Wrong instrument type in bcor file! Instrument = ', CF_INSTRUM
    !  call abort()
    !endif
    !write(*,*) ' Tb will be adjusted (internally) for scan-position bias dependency.'

    write(*,*) 'midas_bgckAtms: empty bcor file is used in the unitTest. Modifying before continuing.'
    call abort()
  endif
   
  if ( modlsqtt ) then
    write(*,*) 'MODLSQ option is activated!'
    write(*,*) '  Output file will contain recomputed values for land/sea qualifier and terrain type based on LG/MG.'
  endif

  IER = FNOM(IUNGEO,'./fstgzmx','STD+RND+R/O',0)

  ! 2) Lecture des statistiques d'erreur totale pour les  TOVS 
  IER = FNOM(IUNSTAT,'./stats_atms_assim','SEQ+FMT',0)
  IF(IER.LT.0)THEN
    write(*,*) '(" bgckAtms: Problem opening ", &
          "ATMS total error statistics file ", stats_atms_assim)'
    CALL ABORT ()
  END IF
  CALL mwbg_readStatTovsAtms(IUNSTAT,INUMSAT,CSATID)
  write(*,*) " SATID's = "
  DO I = 1, INUMSAT
    write(*,*) '  ', CSATID(I)
  ENDDO

  ! 3) Lecture des champs geophysiques (MF/MX) du modele
  IER = FSTOUV(IUNGEO,'RND')

  ! TOPOGRAPHIE (MF ou MX).
  !     MF est la topographie filtree avec unites en metres (filtered ME).
  !     MX est la topographie filtree avec unites en m2/s2  (geopotential topography).
  TOPOFACT = 1.0
  IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','MF')
  CLNOMVAR = 'MF'
  IF (IREC .LT. 0) THEN
    TOPOFACT = 9.80616
    IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','MX')
    CLNOMVAR = 'MX'
  ENDIF
  IF (IREC .LT. 0) THEN
    write(*,*) ' LA TOPOGRAPHIE (MF or MX) EST INEXISTANTE' 
    CALL ABORT()
  ELSE
    ALLOCATE ( MT(NI*NJ), STAT=ier)
    IER = FSTLIR(MT,IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1, &
         ' ',CLNOMVAR)
  ENDIF
      
  IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
      IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10,  &
      IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1, &
      IG2, IG3, IG4, IDUM12, IDUM13, IDUM14,  &
      IDUM15, IDUM16, IDUM17, IDUM18 )
  
  write(*,*) ' GRTYP = ', grtyp 

  ! initialisation
  Call BURP_Init(File_in,  F2=File_out,  IOSTAT=error)
  Call BURP_Init(Rpt_in,   R2=Rpt_out,   IOSTAT=error)

  ! Set BURP "missing value" for reals
  opt_missing = 'MISSING'
  Call BURP_Set_Options(REAL_OPTNAME=opt_missing,REAL_OPTNAME_VALUE=zmisg)

  ! ouverture du fichier burp d'entree et de sortie
  Call BURP_New(File_in,  FILENAME= brp_in,  MODE= FILE_ACC_READ,   IOSTAT= error)
  Call BURP_New(File_out, FILENAME= brp_out, MODE= FILE_ACC_CREATE, IOSTAT= error)

  ! Number of reports and maximum report size from input BURP file
  Call BURP_Get_Property(File_in, NRPTS=nb_rpts, IO_UNIT= iun_burpin)
  if (nb_rpts.le.1) then
    write(*,*) 'The input BURP file ''', trim(brp_in), ''' is empty!'
    stop
  end if

  nsize = MRFMXL(iun_burpin)

  write(*,*)
  write(*,*) 'Number of reports containing observations = ', nb_rpts-1
  write(*,*) 'Size of largest report = ', nsize
  write(*,*)

  ! Add nsize to report size to accomodate modified (larger) data blocks
  nsize = nsize*3

  allocate(adresses(nb_rpts), stat=alloc_status(1))
  
  if (alloc_status(1) /= 0) then
    write(*,*) 'ERROR - allocate(adresses(nb_rpts)). alloc_status =' , alloc_status(1)
    call abort()
  endif
  
  adresses(:) = 0

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

  ! Optional scan position adjustment of Tb: 
  !    -- Get global scan bias array glbscanb(nchanAtms,nscan,nsat) from 
  !       bias correction file
  !        numan(numbsat) = number of channels (out)
  !        listan = list of channel numbers (out)
  !        csatid = list of satellites (dimension(3)) (out)
  !    -- Compute Tb adjustments (dglbscanb)
  if (sp_adj_tb) then
    kk = mxscan/2
    dglbscanb(:,:,:) = 0.0
    csatid(:) = 'XXXXXXXXX'
    numbsat = NSAT
    numan(:) = 0
    do j = 1, nchanAtms
      listan(j,:) = j
    end do

    do ii = 1,NSAT
      csatid(ii) = SATNAMES(ii)
      numan(ii)  = RCNCHAN(ii)
      if ( numan(ii) /= nchanAtms ) then
        write(*,*) 'INFO:',csatid(ii),': Number of channels in coeff file is less than ', nchanAtms
        write(*,*) '      Channels are ', CHANNUMS(ii,1:rcnchan(ii))
      endif

      ! Compute the scan position bias (SPB) for channels with bcor data
      do jj = 1, RCNCHAN(ii)
        j = CHANNUMS(ii,jj)
        dglbscanb(j,:,ii) =  FOVBIAS(ii,jj,:) - FOVBIAS(ii,jj,kk)
      end do
    end do

    if (debug) then
      write(*,*) 'Finished reading bcor coeff file.'
      write(*,*) 'Number of satellites in file = ', numbsat
      write(*,*) 'Number of channels per satellite = ', numan(1:numbsat)
      do kk=1,numbsat
        write(*,*) csatid(kk)
        write(*,*) '  Channels are ', CHANNUMS(kk,1:rcnchan(ii))
        write(*,*) '  Scan biases for each channel are:'
        do j = 1, nchanAtms
          write(*,*) j, dglbscanb(j,:,kk)
        end do
      end do
    end if
       
  end if

  ! LOOP OVER ALL REPORTS OF THE INPUT FILE, APPLY PROCESSING, AND WRITE TO OUTPUT FILE.

  ! Initial scan of file to get number of reports and number of data locations.
  ! Store address of each report in array adresses(nb_rpts) for main REPORTS loop
  ref_rpt = 0
  compteur = 0
  nobs_tot = 0

  do
    ref_rpt = BURP_Find_Report(File_in, REPORT= Rpt_in, SEARCH_FROM= ref_rpt, IOSTAT= error)
    if (error /= burp_noerr) call handle_error()
    if (ref_rpt < 0) Exit
    
    Call BURP_Get_Property(Rpt_in,TEMPS=ibrptime,ELEV=nlocs,STNID=id,RUNN=irun)  
    ! ELEV= the number of locations in the data box (for grouped data) ==> nt in each block
    
    if ( id(1:2) .eq. ">>" ) then
      write(*,*) 'Type de fichier a l_entree = ',id 
      if (id .ne. ">>DERIALT") then
        write(*,*) 'WARNING - le type de fichier devrait etre >>DERIALT'
      endif
    elseif (id(1:1) .eq. "^" ) then
      if ( nlocs > mxnt ) then
        write(*,*) 'ERROR: Number of locations (nlocs) in report ',compteur+1, ' exceeds limit (mxnt)!'
        write(*,*) '       nlocs = ', nlocs
        write(*,*) '       mxnt  = ', mxnt
        call handle_error()
      endif
      nobs_tot = nobs_tot + nlocs
    endif
    compteur = compteur+1
    adresses(compteur) = ref_rpt
    
  end do

  write(*,*) ' Scan 1: Number of reports in input BURP file (compteur) = ', compteur
  write(*,*) '         Number of data locations (nobs_tot)             = ', nobs_tot

  ! if no reports ABORT
  if ( compteur == 0 ) call handle_error()

  ! if no observations STOP
  if ( nobs_tot == 0 ) then
    Call BURP_Free(File_in,F2=File_out)
    Call BURP_Free(Rpt_in,R2=Rpt_out)
    STOP
  end if

  ! MAIN LOOP through all the reports in the file
  nobs_tot = 0
  n_reps = 0
  n_bad_reps = 0
  n_reps_tb2misg = 0
  
  REPORTS: do reportIndex = 1, compteur

    resume_report = .false.

    Call BURP_Get_Report(File_in, REPORT= Rpt_in, REF= adresses(reportIndex), IOSTAT= error) 
    if (error /= burp_noerr) call handle_error()

    Call BURP_Get_Property(Rpt_in,STNID=id,IDTYP=idtyp,ELEV=nlocs,LATI=blat,LONG=blon,NBLK=nblocs,HANDLE=handle)

    if ( id(1:2) .eq. ">>" ) then
      resume_report = .true.

      ! change the header
      Call BURP_Set_Property(Rpt_in,STNID=ETIKRESU)  

      Call BURP_Write_Report(File_out,Rpt_in,IOSTAT=error)
      if (error /= burp_noerr)  call handle_error()
      CYCLE REPORTS
    else
      ! Create new report (Rpt_out) to contain modified blocks from Rpt_in
      Call BURP_New(Rpt_out, Alloc_Space = nsize,  IOSTAT=error)
      if (error /= burp_noerr)  call handle_error()
    
      ! initiliser pour ecriture a File_out
      Call BURP_INIT_Report_Write(File_out,Rpt_out,IOSTAT=error)
      if (error /= burp_noerr)  call handle_error()

      !  copier le header du rapport 
      Call BURP_Copy_Header(TO= Rpt_out, FROM= Rpt_in)

      stnid = id

    endif

    IF ( .not. resume_report ) THEN 

      nt = nlocs
    
      ! Increment total number of obs pts read
      nobs_tot = nobs_tot + nt
      n_reps = n_reps + 1

      ! Allocate arrays to hold data for each location of this report
      alloc_status(:) = 0
      allocate( ident(nt),    stat=alloc_status(1) )
      allocate( waterobs(nt), stat=alloc_status(2) )
      allocate( grossrej(nt), stat=alloc_status(3) )
      allocate( cloudobs(nt), stat=alloc_status(4) )
      allocate( iwvreject(nt),stat=alloc_status(5) )
      allocate( lflagchn(nt,nchanAtms), stat=alloc_status(6) )
      allocate( rclw(nt),     stat=alloc_status(7) )
      allocate( riwv(nt),     stat=alloc_status(8) )
      !allocate( isecs(nt),    stat=alloc_status(9) )
      allocate( precipobs(nt), stat=alloc_status(10) )
      allocate( iber(nt),     stat=alloc_status(11) )
      allocate( ztb89(nt),    stat=alloc_status(12) )
      allocate( ztb150(nt),   stat=alloc_status(13) )
      allocate( bcor150(nt) )
      allocate( scatl(nt),    stat=alloc_status(14) )
      allocate( scatw(nt),    stat=alloc_status(15) )
      allocate( ztb_amsub3(nt), stat=alloc_status(16) )
      allocate( bcor_amsub3(nt) )
      allocate( ztb_amsub5(nt), stat=alloc_status(17) )
      allocate( bcor_amsub5(nt) )
      allocate( zdi(nt),      stat=alloc_status(18) )
      allocate( err(nt),      stat=alloc_status(19) )
      allocate( ascatw(nt),   stat=alloc_status(20) )
      allocate( tb23(nt),     stat=alloc_status(21) )
      allocate( bcor23(nt) )
      allocate( tb31(nt),     stat=alloc_status(22) )
      allocate( bcor31(nt) )
      allocate( tb50(nt),     stat=alloc_status(23) )
      allocate( bcor50(nt) )
      allocate( tb53(nt),     stat=alloc_status(24) )
      allocate( tb89(nt),     stat=alloc_status(25) )
      allocate( bcor89(nt) )
      allocate( lsq(nt),      stat=alloc_status(26) )
      allocate( trn(nt),      stat=alloc_status(27) )
      allocate( scatec(nt),   stat=alloc_status(28) )
      allocate( scatbg(nt),   stat=alloc_status(29) )
      allocate( SeaIce(nt),   stat=alloc_status(30) )
      allocate( lqc(nt,nchanAtms), stat=alloc_status(31) )
     
      if( any(alloc_status /= 0) ) then
        write(*,*) ' midas_bgckAtms: Memory allocation error '
        call abort()
      endif
    
      ident(:) = 0        ! filter information flag; set all bits OFF
     
      ! Information flag (ident) values (new BURP element 025174 in header)

      ! BIT    Meaning
      !  0     off=land or sea-ice, on=open water away from coast
      !  1     Mean 183 Ghz [ch. 18-22] is missing
      !  2     CLW is missing (over water)
      !  3     CLW > clw_atms_nrl_LTrej (0.175 kg/m2) (cloudobs)
      !  4     scatec/scatbg > Lower Troposphere limit 9/10 (precipobs)
      !  5     Mean 183 Ghz [ch. 18-22] Tb < 240K
      !  6     CLW > clw_atms_nrl_UTrej (0.200 kg/m2)
      !  7     Dryness Index rejection (for ch. 22)
      !  8     scatec/scatbg > Upper Troposphere limit 18/15
      !  9     Dryness Index rejection (for ch. 21)
      ! 10     Sea ice > 0.55 detected
      ! 11     Gross error in Tb (any chan.)  (all channels rejected)


      !  Get all the required data from the blocks in the report (Rpt_in)
      call mwbg_getData(reportIndex, Rpt_in, ISAT, zenith, ilq, itt, zlat, zlon, ztb, &
                        biasCorr, ZOMP, scanpos, nvalOut, ntOut, qcflag1, qcflag2, &
                        ican, icanomp, IMARQ, IORBIT, 'ATMS')

      if ( ALL(ZOMP(:) == MPC_missingValue_R4 )) then
        n_bad_reps = n_bad_reps + 1  
      
        Call BURP_Free(Rpt_out,IOSTAT=error)
        if (error /= burp_noerr)  call handle_error()

        cycle REPORTS
      end if

      ! Initialize internal land/sea qualifier and terrain type arrays to values
      ! read from file

      lsq(:) = ilq(1:nt)  ! land/sea qualifier
      trn(:) = itt(1:nt)  ! terrain type (sea-ice)

      ! Set WATEROBS() array and reset lsq, trn 

      ! Determine which obs pts are over open water (i.e NOT near coasts or
      ! over/near land/ice) using model MG and LG fields from glbhyb2 ANAL
      !  MG = land/sea mask field (0.0 (water) to 1.0 (land))
      !  LG = ice fraction field  (0.0 - 1.0)
      !  lsq = 0 (land), 1 (water)
      !  trn = -1 (no ice/snow),  0 (ice)
      !  NOTE: mwbg_landIceMaskAtms redefines lsq and trn based on interpolated 
      !  MG, LG fields so that
      !    lsq = 1 (point over water away from land/coast), 0 (land/coast) otherwise
      !    trn = 0 (point over or near sea-ice),           -1 (ice free) otherwise
      !
      !  waterobs(:)=.true. at points where lsq = 1 and trn = -1
      call mwbg_landIceMaskAtms(mglg_file,nt,zlat,zlon,lsq,trn,waterobs)

      ! Check for values of TB that are missing or outside physical limits.
      ! **NOTE: REJECT ALL CHANNELS IF ONE IS FOUND TO BE BAD.

      grossrej(:)  = .false.
      call mwbg_grossValueCheck(nt,ztb,grossrej)
      
      if ( ANY(grossrej) ) then
        write(*,*) ' mwbg_grossValueCheck has detected bad Tb data. Number of affected locations = ', COUNT(grossrej)
        write(*,*) '   Box lat, lon = ', blat, blon
      endif

      ! Preliminary QC checks --> set lqc(nt,nchanAtms)=.true. for data that fail QC

      lqc(:,:) = .false.  ! Flag for preliminary QC checks
      call mwbg_firstQcCheckAtms(zenith, ilq, itt, zlat, zlon, ztb, scanpos, stnid, &
                                 nvalOut, nt, lqc, grossrej, lsq, trn, qcflag1, qcflag2, &
                                 ican, blat, blon, lutb)

      if ( lutb ) n_reps_tb2misg = n_reps_tb2misg + 1

      ! Output of mwbg_firstQcCheckAtms
      !   lqc=true (entire array) --> problem with channels (number of channels or the channel numbers)
      !                               file could be corrupted [abort]
      !   ztb(nchanAtms) = zmisg  (at locations with bad zenith angle and/or lat,lon)
      !   zenith     = zmisg  (if bad zenith angle)
      !   lutb=true  (if ztb set to zmisg at 1 or more locations)
      ! All channels are flagged for all checks except the "data level" QC flag check.

      !if ( lutb ) then
      !  write(*,*) ' Number of Tb data = zmisg after mwbg_firstQcCheckAtms = ', COUNT(ztb == zmisg)
      !  write(*,*) ' Number of Tb data = 330.04              = ', COUNT(ztb == 330.04)
      !  write(*,*) ' Total number of data                    = ', nval*nt
      !endif

      bad_report = .false.
      if ( COUNT(lqc) == nt*nchanAtms ) then
        write(*,*) ' mwbg_firstQcCheckAtms has detected a problem with data in this report!'
        write(*,*) '   Report box lat, lon = ', blat, blon
        bad_report = .true.
        n_bad_reps = n_bad_reps + 1
      endif
     
      IF (.not. bad_report) THEN

        !  Exclude problem points from further calculations

        do kk = 1,nt
          if ( COUNT(lqc(kk,:)) == nchanAtms ) grossrej(kk) = .true.
        enddo

        where ( grossrej ) ident = IBSET(ident,11)

        ! Apply NRL cloud filter, scattering index and sea-ice detection algorithms to 
        !   OPEN WATER (waterobs=true) points.
        ! Points with SeaIce>0.55 are set to sea-ice points (waterobs --> false)

        ! To begin, assume that all obs are good.

        cloudobs(:)  = .false.
        iwvreject(:) = .false.
        precipobs(:) = .false.

        ! First remove scan-dependency of Tb biases for data in the report (OPTION)
        if (sp_adj_tb) then    
          ii = 0

          ! Find satellite index (satellite "^NPP")
          do kk = 1, numbsat
            if ( TRIM(csatid(kk)) == TRIM(stnid(2:9)) ) ii = kk
          end do
          if ( ii == 0 ) then
            write(*,*) ' Error: Satellite not found in bias correction file!'
            write(*,*) '        Satellite = ', stnid(2:4)
            write(*,*) '        Satellites in BCOR file = ', csatid(:)
            call abort()
          end if

          ! Extract the Tb adjustments for each channel for this satellite
          do kk = 1, nchanAtms
            zbcor(kk,:) = dglbscanb(kk,:,ii)
            if ( debug ) &
              write(*,*) 'Scan bias adjustments for channel ', kk, ' : ', zbcor(kk,:)
          end do

          ! Adjust Tb for for each channel according to scan position
          indx1 = 1
          do kk = 1, nt   !  loop over NT locations in report
            indx2 = kk*nchanAtms
            if ( debug ) then
              write(*,*) 'location, indx1, indx2 = ', kk, indx1, indx2
            end if
            do jj = 1, nchanAtms
              ztbcor(indx1+jj-1) = ztb(indx1+jj-1) - zbcor(jj,scanpos(kk))
              if ( debug ) then
                write(*,*) 'scanpos, ztb index = ', scanpos(kk), indx1+jj-1
                write(*,*) 'channel, ztb, zbcor, ztbcor = ', &
                    jj, ztb(indx1+jj-1), zbcor(jj,scanpos(kk)), ztbcor(indx1+jj-1)
              end if
            end do  
            indx1 = indx2 + 1
          end do
           
        else  ! no correction
          ztbcor(:) = ztb(:)
           
        end if
     
        ! extract required channels:
        !  23 Ghz = AMSU-A 1 = ATMS channel 1 
        !  31 Ghz = AMSU-A 2 = ATMS channel 2
        !  50 Ghz = AMSU-A 3 = ATMS channel 3
        !  53 Ghz = AMSU-A 5 = ATMS channel 6
        !  89 Ghz = AMSU-A15 = ATMS channel 16
        ! 150 Ghz = AMSU-B 2 = ATMS channel 17
        !
        ! Extract Tb for channels 16 (AMSU-B 1) and 17 (AMSU-B 2) for Bennartz SI
        ! Extract Tb for channels 22 (AMSU-B 3) and 18 (AMSU-B 5) for Dryness Index (DI)

        indx1 = 1
        do ii = 1, nt
          indx2 = ii*nchanAtms
          tb23(ii)      = ztbcor(indx1)
          bcor23(ii)    = biasCorr(indx1)
          tb31(ii)      = ztbcor(indx1+1)
          bcor31(ii)    = biasCorr(indx1+1)
          tb50(ii)      = ztbcor(indx1+2)
          bcor50(ii)    = biasCorr(indx1+2)
          tb53(ii)      = ztbcor(indx1+5)
          tb89(ii)      = ztbcor(indx1+15)
          ztb89(ii)     = tb89(ii)
          bcor89(ii)    = biasCorr(indx1+15)
          ztb150(ii)    = ztbcor(indx1+16)
          bcor150(ii)    = biasCorr(indx1+16)
          ztb_amsub3(ii) = ztbcor(indx1+21)
          bcor_amsub3(ii) = biasCorr(indx1+21)
          ztb_amsub5(ii) = ztbcor(indx1+17)
          bcor_amsub5(ii) = biasCorr(indx1+17)
          indx1 = indx2 + 1
        end do

        !  mwbg_nrlFilterAtms returns rclw, scatec, scatbg and also does sea-ice detection
        !  Missing value for  rclw, scatec, scatbg  is -99.0 (e.g. over land or sea-ice).
        !  Sets trn=0 (sea ice) for points where retrieved SeaIce>=0.55.
        !  Does nothing if trn=0 (sea ice) and retrieved SeaIce<0.55.

        call mwbg_nrlFilterAtms(err, nt, tb23, bcor23, tb31, bcor31, tb50, bcor50, &
                                tb89, bcor89, ztb150, bcor150, zenith, zlat, lsq, trn, &
                                waterobs, grossrej, rclw, scatec, scatbg, iNumSeaIce, iRej, &
                                SeaIce)
        
        ! Flag data using NRL criteria

        ! Compute Mean 183 Ghz [ch. 18-22] Tb (riwv)
        riwv = -99.0
        indx1 = 1
        do ii = 1, nt
          indx2 = ii*nchanAtms
          if (.not.grossrej(ii)) then
            if ( useUnbiasedObsForClw ) then
              ztb183(1) = ztbcor(indx1+17)
              ztb183(2) = ztbcor(indx1+18)
              ztb183(3) = ztbcor(indx1+19)
              ztb183(4) = ztbcor(indx1+20)
              ztb183(5) = ztbcor(indx1+21)
            else
              ztb183(1) = ztbcor(indx1+17) - biasCorr(indx1+17)
              ztb183(2) = ztbcor(indx1+18) - biasCorr(indx1+18)
              ztb183(3) = ztbcor(indx1+19) - biasCorr(indx1+19)
              ztb183(4) = ztbcor(indx1+20) - biasCorr(indx1+20)
              ztb183(5) = ztbcor(indx1+21) - biasCorr(indx1+21)
            end if
            riwv(ii)  = sum(ztb183)/5.0
            if ( riwv(ii) < mean_Tb_183Ghz_min ) iwvreject(ii) = .true.
          else
            iwvreject(ii) = .true.
          endif
          indx1 = indx2 + 1
        end do

        !  Set bits in ident flag to identify where various data selection criteria are met
        !     precipobs = .true  where ECMWF or BG scattering index > min_threshold (LT)
        !     cloudobs  = .true. where CLW > min_threshold (LT) or if precipobs = .true

        where ( scatec .gt. scatec_atms_nrl_LTrej .or. scatbg .gt. scatbg_atms_nrl_LTrej ) precipobs = .true.
        n_cld = count(rclw .gt. clw_atms_nrl_LTrej)
        cldcnt  = cldcnt  + n_cld
        where ( (rclw .gt. clw_atms_nrl_LTrej) .or. precipobs ) cloudobs = .true.
        where ( waterobs )  ident = IBSET(ident,0)
        where ( iwvreject ) ident = IBSET(ident,5)
        where ( precipobs ) ident = IBSET(ident,4)
        where ( rclw .gt. clw_atms_nrl_LTrej) ident = IBSET(ident,3)
        where ( rclw .gt. clw_atms_nrl_UTrej) ident = IBSET(ident,6)
        where ( scatec .gt. scatec_atms_nrl_UTrej .or. scatbg .gt. scatbg_atms_nrl_UTrej ) ident = IBSET(ident,8)
        where ( SeaIce .ge. 0.55 ) ident = IBSET(ident,10)
        
        where ( waterobs .and. (rclw == -99.) ) ident = IBSET(ident,2)
        where ( riwv == -99.)                   ident = IBSET(ident,1)

        ! Compute the simple AMSU-B Dryness Index zdi for all points = Tb(ch.3)-Tb(ch.5)
        if ( useUnbiasedObsForClw ) then
          where ( .not.grossrej )
            zdi = ztb_amsub3 - ztb_amsub5
          elsewhere
            zdi = zmisg
          end where
        else
          where ( .not.grossrej )
            zdi = (ztb_amsub3 - bcor_amsub3) - (ztb_amsub5 - bcor_amsub5)
          elsewhere
            zdi = zmisg
          end where
        end if

        ! Review all the checks previously made to determine which obs are to be accepted
        ! for assimilation and which are to be flagged for exclusion (lflagchn). 
        !   grossrej()  = .true. if any channel had a gross error at the point
        !   cloudobs()  = .true. if CLW > clw_atms_nrl_LTrej (0.175) or precipobs
        !   precipobs() = .true. if precip. detected through NRL scattering indices
        !   waterobs()  = .true. if open water point
        !   iwvreject() = .true. if Mean 183 Ghz [ch. 18-22] Tb < 240K (too dry for ch.20-22 over land)

        lflagchn(:,:) = lqc(:,:)  ! initialize with flags set in mwbg_firstQcCheckAtms
        do kk = 1, nt
          ! Reject all channels if gross Tb error detected in any channel or other problems 
          if ( grossrej(kk) ) then
            lflagchn(kk,:) = .true.

          else
            ! OVER LAND OR SEA-ICE,
            !    -- CLW/SI not determined over land
            !    -- surface emissivity effects lower tropospheric and window channels     
            !    -- reject window & lower tropospheric channels 1-6, 16-19
            !    -- reject ch. 20-22 if iwvreject = .true.  [ Mean 183 Ghz [ch. 18-22] Tb < 240K ]
            !    -- check DI for AMSU-B like channels
            if  ( .not. waterobs(kk) ) then
              lflagchn(kk,1:ipc)     = .true.      ! AMSU-A 1-6
              lflagchn(kk,16:19)     = .true.      ! AMSU-B (like 1,2,5)
              if ( iwvreject(kk) ) lflagchn(kk,20:22) = .true.  ! AMSU-B (like 4,3)

              ! Dryness index (for AMSU-B channels 19-22 assimilated over land/sea-ice)
              ! Channel AMSUB-3 (ATMS channel 22) is rejected for a dryness index >    0.
              !                (ATMS channel 21) is rejected for a dryness index >   -5.
              ! Channel AMSUB-4 (ATMS channel 20) is rejected for a dryness index >   -8.
              if ( zdi(kk) > 0.0 ) then
                lflagchn(kk,22) = .true.
                ident(kk) = IBSET(ident(kk),7)
              endif
              if ( zdi(kk) > -5.0 ) then
                lflagchn(kk,21) = .true.
                ident(kk) = IBSET(ident(kk),9)
                drycnt = drycnt + 1
              endif
              if ( zdi(kk) > -8.0 ) then
                lflagchn(kk,20) = .true.
              endif
            endif

            ! OVER WATER,
            !    -- reject ch. 1-6, 16-20 if CLW > clw_atms_nrl_LTrej or CLW = -99.0
            !    -- reject ch. 7-9, 21-22 if CLW > clw_atms_nrl_UTrej or CLW = -99.0
            !    -- reject ch. 1-6, 16-22 if scatec > 9  or scatec = -99.0
            !    -- reject ch. 7-9        if scatec > 18 or scatec = -99.0
            !    -- reject ch. 1-6        if scatbg > 10 or scatbg = -99.0
            !    -- reject ch. 7-9        if scatbg > 15 or scatbg = -99.0
            !    -- reject ch. 16-22      if iwvreject = .true.   [ Mean 183 Ghz [ch. 18-22] Tb < 240K ]
            if  ( waterobs(kk) ) then
              if ( rclw(kk)   >  clw_atms_nrl_LTrej )  then
                lflagchn(kk,1:ipc) = .true.
                lflagchn(kk,16:20) = .true. 
              endif
              if ( rclw(kk)   >  clw_atms_nrl_UTrej )  then
                lflagchn(kk,7:9)   = .true.
                lflagchn(kk,21:22) = .true. 
              endif
              if ( scatec(kk) >  scatec_atms_nrl_LTrej ) then
                lflagchn(kk,1:ipc) = .true.
                lflagchn(kk,16:22) = .true.
              endif
              if ( scatec(kk) >  scatec_atms_nrl_UTrej ) lflagchn(kk,7:9) = .true.
              if ( scatbg(kk) >  scatbg_atms_nrl_LTrej ) lflagchn(kk,1:ipc) = .true.
              if ( scatbg(kk) >  scatbg_atms_nrl_UTrej ) lflagchn(kk,7:9) = .true.
              if ( iwvreject(kk) ) lflagchn(kk,16:22) = .true.
              if ( rclw(kk) == -99. ) then
                ident(kk) = IBSET(ident(kk),2)
                lflagchn(kk,1:9)   = .true.
                lflagchn(kk,16:22) = .true.
              endif
              if ( riwv(kk) == -99. ) then     ! riwv = mean_Tb_183Ghz
                ident(kk) = IBSET(ident(kk),1)
                lflagchn(kk,16:22) = .true.
              endif           
            endif
         
          endif

          if ( .not. waterobs(kk) ) landcnt  = landcnt  + 1
          if ( grossrej(kk) )  rejcnt = rejcnt + 1
          if ( iwvreject(kk))  iwvcnt = iwvcnt + 1
          if ( precipobs(kk) .and. waterobs(kk) ) then
            pcpcnt = pcpcnt + 1
          endif
          
          if ( ANY(lflagchn(kk,:)) ) flgcnt = flgcnt + 1

        end do

        ! RESET riwv array to ECMWF scattering index for output to BURP file
        riwv(:) = scatec(:)

        ! Set missing rclw and riwv to BURP missing value (zmisg)
        where (rclw == -99. ) rclw = zmisg
        where (riwv == -99. ) riwv = zmisg


        ! Modify the blocks in Rpt_in and write to Rpt_out
        ! - Modify flag values so that the obs identified above as being over land/ice,
        !   or in cloudy/precip regions, etc. are not assimilated (FLAG block 15392/15408).
        ! - OPTION: update land-sea qualifier and terrain type in INFO block 3072.
        ! - Update Tb data in DATA block 9248/9264 (if Tb was modified).
        ! - Add new elements to INFO block 3072.
        ! - Modify 24bit global flags in 3D block 5120 (if any data rejected).
        call mwbg_writeBlocks(reportIndex, ztb, lsq, trn, riwv, rclw, ident, &
                              lflagchn, IMARQ, lutb, Rpt_in, Rpt_out)

        !** start second quality control (atms_inovqc standalone program) **
        !
        ! trouver l'indice du satellite
        INOSAT = 0
        DO I = 1,MXSAT
          IF ( STNID .EQ. '^'//CSATID(I) ) THEN
            INOSAT = I
          ENDIF
        ENDDO
        IF ( INOSAT .EQ. 0 ) THEN
          write(*,*)'SATELLITE NON-VALIDE', STNID
          CALL ABORT()
        ENDIF

        ! 5) Interpolation de le champ MF/MX (topogrpahy) aux pts TOVS.
        !    N.B.: on examine ce champ sur une boite centree sur chaque obs.
        NLAT = (MXLAT-1)/2
        NLON = (MXLON-1)/2
        DO JN = 1, NT
          INDX = 0
          DO I = -NLAT, NLAT
            XLAT = ZLAT(JN) +I*DLAT
            XLAT = MAX(-90.0,MIN(90.0,XLAT))
            DO J = -NLON, NLON
              INDX = INDX + 1
              XLON = ZLON(JN) +J*DLON
              IF ( XLON .LT. -180. ) XLON = XLON + 360.
              IF ( XLON .GT.  180. ) XLON = XLON - 360.
              IF ( XLON .lt.    0. ) XLON = XLON + 360.
              ZLATBOX(INDX,JN) = XLAT
              ZLONBOX(INDX,JN) = XLON
            ENDDO
          ENDDO
        ENDDO

        ier  = ezsetopt('INTERP_DEGREE','LINEAR')
        gdmg = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
        ier  = gdllsval (gdmg,mtintbox,mt,ZLATBOX,ZLONBOX,MXLAT*MXLON*NT)

        DO JN = 1, NT
          IF (DEBUG) THEN
            PRINT *, ' ------------------  '
            PRINT *, ' JN = ', JN
            PRINT *, '   '
            PRINT *, ' zlat,zlon = ', zlat(jn), zlon(jn)
            PRINT *, '   '
            PRINT *, ' ZLATBOX = '
            PRINT *,  (ZLATBOX(I,JN),I=1,MXLAT*MXLON)
            PRINT *, ' ZLONBOX = '
            PRINT *,  (ZLONBOX(I,JN),I=1,MXLAT*MXLON)
            PRINT *, ' MTINTBOX = '
            PRINT *,  (MTINTBOX(I,JN),I=1,MXLAT*MXLON)
          ENDIF
          MTINTRP(JN) = 0.0
          DO I=1,MXLAT*MXLON
            MTINTRP(JN) = MAX(MTINTRP(JN),MTINTBOX(I,JN)/TOPOFACT)
          ENDDO
          IF (DEBUG) THEN
            PRINT *, ' MTINTRP = ', MTINTRP(JN)
          ENDIF
        ENDDO

        ! 6) Controle de qualite des TOVS. Data QC flags (IMARQ) are modified here!
        CALL mwbg_tovCheckAtms(ISAT, IORBIT, ican, ICANOMP, ztb, biasCorr, &
                               ZOMP, ICHECK, nvalOut, nvalOut, NT, ZMISG, INOSAT, ident, &
                               ICHKPRF, scanpos, MTINTRP, IMARQ, STNID, RESETQC)

        ! Accumuler Les statistiques sur les rejets
        CALL mwbg_qcStatsAtms(INUMSAT, ICHECK, ican, INOSAT, CSATID, nvalOut, &
                              NT, .FALSE.)

        ! 7) Mise a jour des marqueurs.
        CALL mwbg_updatFlgAtms(ICHKPRF, ICHECK, RESETQC, IMARQ, Rpt_out)
      ENDIF

      alloc_status(:) = 0
      deallocate( ident,    stat=alloc_status(1) )
      deallocate( waterobs, stat=alloc_status(2) )
      deallocate( grossrej, stat=alloc_status(3) )
      deallocate( cloudobs, stat=alloc_status(4) )
      deallocate( iwvreject,stat=alloc_status(5) )
      deallocate( lflagchn, stat=alloc_status(6) )
      deallocate( rclw,     stat=alloc_status(7) )
      deallocate( riwv,     stat=alloc_status(8) )
      !deallocate( isecs,    stat=alloc_status(9) )
      deallocate( precipobs, stat=alloc_status(10) )
      deallocate( iber,     stat=alloc_status(11) )
      deallocate( ztb89,    stat=alloc_status(12) )
      deallocate( ztb150,   stat=alloc_status(13) )
      deallocate( bcor150 )
      deallocate( scatl,    stat=alloc_status(14) )
      deallocate( scatw,    stat=alloc_status(15) )
      deallocate( ztb_amsub3, stat=alloc_status(16) )
      deallocate( bcor_amsub3 )
      deallocate( ztb_amsub5, stat=alloc_status(17) )
      deallocate( bcor_amsub5 )
      deallocate( zdi,      stat=alloc_status(18) )
      deallocate( err,      stat=alloc_status(19) )
      deallocate( ascatw,   stat=alloc_status(20) )
      deallocate( tb23,     stat=alloc_status(21) )
      deallocate( bcor23 )
      deallocate( tb31,     stat=alloc_status(22) )
      deallocate( bcor31 )
      deallocate( tb50,     stat=alloc_status(23) )
      deallocate( bcor50 )
      deallocate( tb53,     stat=alloc_status(24) )
      deallocate( tb89,     stat=alloc_status(25) )
      deallocate( bcor89 )
      deallocate( lsq,      stat=alloc_status(26) )
      deallocate( trn,      stat=alloc_status(27) )
      deallocate( scatec,   stat=alloc_status(28) )
      deallocate( scatbg,   stat=alloc_status(29) )
      deallocate( SeaIce,   stat=alloc_status(30) )
      deallocate( lqc,      stat=alloc_status(31) )

      if( any(alloc_status /= 0) ) then
        write(*,*) ' midas_bgckAtms: Memory deallocation error '
        call abort()
      endif

    ENDIF

    ! Write the modified report to the output file
    IF ( .not. resume_report ) THEN
      if (.not. bad_report ) then
        Call BURP_Write_Report(File_out,Rpt_out,IOSTAT=error)
        if (error /= burp_noerr)  call handle_error()
      endif
      Call BURP_Free(Rpt_out,IOSTAT=error)
      if (error /= burp_noerr)  call handle_error()
    ENDIF

  end do REPORTS

  write(*,*) ' --------------------------------------------------------------- '
  write(*,*) ' Number of obs pts read from BURP file              = ', nobs_tot
  write(*,*) ' Number of BURP file reports                        = ', n_reps
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

  ! 9) Fin
  ! Imprimer les statistiques sur les rejets
  CALL mwbg_qcStatsAtms(INUMSAT, ICHECK, ican, INOSAT, CSATID, nvalOut, &
                        NT, .TRUE.)

  Deallocate(adresses)

  istat = exfin('midas_bgckAtms','FIN','NON')

  ! fermeture des fichiers 
  Call BURP_Free(File_in,F2=File_out,IOSTAT=error)
  Call BURP_Free(Rpt_in,R2=Rpt_out,IOSTAT=error)
  ISTAT = FSTFRM(IUNGEO)
  ISTAT = FCLOS (IUNGEO)
  ISTAT = FCLOS (IUNSTAT)

  STOP

  contains


    subroutine handle_error()
      implicit none

      write(*,*) BURP_STR_ERROR()
      write(*,*) "history"
      Call BURP_STR_ERROR_HISTORY()
      Deallocate(adresses)
      Call BURP_Free(File_in,F2=File_out)
      Call BURP_Free(Rpt_in,R2=Rpt_out)
      call abort()
    end subroutine handle_error

end program midas_bgckAtms

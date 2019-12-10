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

program midas_satQCATMS
  !
  ! :Purpose: Main program for background check of microwave instruments. 
  !
  use burp_module

!  Language: FORTRAN 90
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
!  Program atms_inovqc.f does O-P based rejections, topography filtering, channel selection and 
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
!------------------------------------------------------------------
!
!  To compile on AIX:
!
!    >  . s.ssmuse.dot Xlf13.108 rmnlib-dev devtools cmda
!    >  s.compile -src satqc_atms.ftn90 -o satqc_atms -libpriv burp_module -librmn rmn_014_rc1
!
!  To RUN:
!
! >  satqc_atms  -IXENT burpin  -OXSRT burpout  { -SPADJUST -IRBC coeff_file_atms } -MGLGFILE mglg_file
! >  satqc_atms  -IXENT atms_deri_so  -OXSRT atms_deri_so_qc1  -MGLGFILE fstglmg { -MODLSQ }
!
!    burpin    =  DERIALT-type ATMS data BURP file
!    burpout   =  DERIALT-type ATMS data BURP file with data FLAG block QC flag bit 7 set for data 
!                 rejected for CLW, precip, surface-type, dryness, etc. 
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
!   OPTIONS:
!    -SPADJUST              = scan bias adjust the Tb data before using the data in algorithms
!    -IRBC coeff_file_atms  = ATMS bias correction coeff file (for -SPADJUST option)
!    -MODLSQ                = update values of lsq and trn elements in output BURP file
!                             (for evaluation purposes only)
!

!===========================================================================================================
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
!----------------------------------------------------------------------------------------------------------------------
!               NRL option
! Channels never assimilated                                      = 1-3, 16-17
! Channels assimilated over water only                            = 4-6, 18-22  CLW/SI check
! Channels assimilated over land and water (with land topo check) = 7-8         no CLW/SI check over land
! Channels assimilated over land and water (no land topo check)   = 9-15
!
! Assim channels sensitive to clouds/precip                       = 4-9, 18-22  (7-9 over land and water)
! Channels insensitive to clouds/precip  (high peaking)           = 10-15       no filters
!----------------------------------------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------------------------------------

   IMPLICIT NONE

! Define FORTRAN FST functions:
!   integer, external :: fnom,fclos
   integer, external :: fstinf,fstprm,fstlir
   integer, external :: fstouv,fstfrm,fstinl,fstvoi

! Set array limits (ATMS: 22 chan, 96 FOV):
   integer, parameter :: mxele=24,mxval=30,mxnt=2800
   integer, parameter :: nchan=22
   integer, parameter :: mxsat=3
   integer, parameter :: mxan=22
   integer, parameter :: mxor=20
   integer, parameter :: mxscan=96
   
   integer :: idum1,idum2,idum3

! Other variables:
   real, parameter    :: zmisg=9.9e09

   character(len=*), parameter :: VERSION = '2.12'

!----------------------------------------------------------------------------------------------------
!  For new BURP_MODULE routines
!----------------------------------------------------------------------------------------------------
   type(BURP_FILE)        :: File_in,File_out
   type(BURP_RPT)         :: Rpt_in,Rpt_out
   type(BURP_BLOCK)       :: Block_in,Block_copy

   integer(kind=int_def)  :: error,nb_rpts,ref_rpt,compteur,ref_blk,nbele,nlocs,my_idtyp
   integer(kind=int_def)  :: handle,nvale,nte,j,i,btyp,my_nt,my_nval,my_nele,ind_lat,ind_lon,bfam
   integer(kind=int_def)  :: indice, indice1, indice2, indice3, kk
   
   integer,allocatable    :: adresses(:)
   integer                :: ibrptime, nblocs, nsize, iun_burpin, irun

   character(len=20)      :: opt_missing
   character(len=90)      :: brp_in,brp_out
   character(len=9)       :: id
!----------------------------------------------------------------------------------------------------
!
   integer, parameter :: iunbc=30
!
   integer :: ierr,ilnmx,status
   integer :: istat,nombre
   integer :: nele,nval,nt,blat,blon
   integer :: flgcnt,rejcnt
   integer :: cldcnt,landcnt,iwvcnt,pcpcnt
   integer :: drycnt
   integer :: indx1, indx2, ii, numbsat, iich, jj, k
   integer :: iNumSeaIce, nobs_tot, n_cld, iRej, n_bad_reps, n_reps, n_reps_tb2misg

   integer, dimension(mxnt)             :: scanpos
   integer, dimension(31)               :: alloc_status
   
! ----- These arrays are for subroutine GETBC, called to get bias corrections -------
   integer, dimension(mxsat)            ::  numor, numan
   integer, dimension(mxan,mxsat)       ::  listan

   real, dimension(mxor,mxan,mxsat)     ::  amcoeff
   real, dimension(mxan,mxscan,mxsat)   ::  glbscanb, dglbscanb

   character(len=9), dimension(mxsat)   ::  csatid
! -------------------------------------------------------------------------------

   real, dimension(nchan,mxscan)        ::  zbcor
   real, dimension(mxval*mxnt)          ::  ztbcor
   
   real, dimension(5)                   ::  ztb183

! -------------------------------------------------------------------------------
!  Upper limit for CLW (kg/m**2) for Tb rejection over water
   real, parameter :: clw_atms_nrl_LTrej=0.175      ! lower trop chans 1-6, 16-20
   real, parameter :: clw_atms_nrl_UTrej=0.2        ! upper trop chans 7-9, 21-22
!! Other NRL thresholds
   real, parameter :: scatec_atms_nrl_LTrej=9.0     ! lower trop chans 1-6, 16-22
   real, parameter :: scatec_atms_nrl_UTrej=18.0    ! upper trop chans 7-9
   real, parameter :: scatbg_atms_nrl_LTrej=10.0    ! lower trop chans 1-6
   real, parameter :: scatbg_atms_nrl_UTrej=15.0    ! upper trop chans 7-9
   real, parameter :: mean_Tb_183Ghz_min=240.0      ! min. value for Mean(Tb) chans. 18-22 
! -------------------------------------------------------------------------------   
!
!  Highest peaking AMSU-A like ATMS channel for ocean-only assimilation
!    6 = 700 mb   (AMSU/operations == AMSU chan. 5)
!    7 = 400 mb   (scat. index used in AMSU/operations -- AMSU chan. 6)
!    8 = 250 mb   (ECMWF)
   integer, parameter :: ipc=6
   

   real, dimension(mxval*mxnt)          :: zdata,ztb
   real, dimension(mxnt)                :: zlat,zlon,zenith
   integer, dimension(mxnt)             :: ilq,itt
   integer, dimension(mxval*mxnt)       :: ican, qcflag2
   integer, dimension(mxnt,3)           :: qcflag1

   character(len=9)                     :: stnid
   character(len=128)                   :: mglg_file

!  F90 Allocatable arrays (deferred-shape)

   integer, allocatable, dimension(:)   :: ident,iber,err,lsq,trn
   real,    allocatable, dimension(:)   :: rclw,riwv,ascatw,SeaIce
   real,    allocatable, dimension(:)   :: ztb89,ztb150,scatl,scatw,tb23,tb31,tb50,tb53,tb89
   real,    allocatable, dimension(:)   :: ztb_amsub3, ztb_amsub5, zdi, scatec, scatbg
   logical, allocatable, dimension(:)   :: waterobs,ukbadobs
   logical, allocatable, dimension(:)   :: cloudobs,iwvreject,precipobs,grossrej
   logical, allocatable, dimension(:,:) :: lflagchn, lqc

!  For s/r READ_COEFF
   INTEGER                                         :: NSAT, NFOV
   INTEGER, PARAMETER                              :: maxpred=6
   character(len=90)                               :: coef_in
   character(len=9), dimension(mxsat)              :: SATNAMES
   integer, dimension(mxsat)                       :: RCNCHAN
   integer, dimension(mxsat,nchan)                 :: CHANNUMS, NPRED
   real,dimension(mxsat,nchan,maxpred+1)           :: COEFF
   real,dimension(mxsat,nchan,mxscan)              :: FOVBIAS
   character(len=5)                                :: CF_INSTRUM
   character(len=2),dimension(mxsat,nchan,maxpred) :: PTYPES
!
   logical :: sp_adj_tb, debug, modlsqtt
   logical :: bad_report, lutb, resume_report

! External functions
   integer, external :: exdb,exfin

   integer :: nulnam
   namelist /nambgck/ debug, sp_adj_tb, modlsqtt

! Initialize the CHARACTER arrays.

   !clist = (/'L.       ','I.       ','IXENT.   ','OXSRT.   ','SPADJUST.',  'IRBC.    ', 'MGLGFILE.', 'MODLSQ.  '/)
!

  brp_in = './obsatms'
  brp_out = './obsatms.out'

  mglg_file  = './fstmglg'
  coef_in = './bcor'

  ierr = fnom(6,'out.txt','SEQ',0)
   
!--------------------------------------------------------------------
  istat=exdb('midas_satQCATMS',VERSION,'NON')
!--------------------------------------------------------------------

  debug = .false.
  sp_adj_tb = .false.
  modlsqtt = .false.

  ! read namelist variables
  ! reading namelist
  nulnam = 0
  ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam, nml=nambgck, iostat=ierr)
  if ( ierr /= 0 ) then
    write(*,*) 'midas_satQCATMS: Error reading namelist'
    call abort()
  end if
  write(*,nml=nambgck)
  ierr = fclos(nulnam)

! Optional adjustment of Tb for scan-position dependency before computing CLW.
! Open and read new UBCOR system ascii file containing mean O-P for each scan position (FOVBIAS).
  if ( sp_adj_tb ) then
    !call read_coeff(SATNAMES,CHANNUMS,FOVBIAS,COEFF,NSAT,RCNCHAN,NFOV,NPRED,CF_INSTRUM,maxpred,iunbc,coef_in,PTYPES)
    !if ( NSAT == 0 ) then
    !  write(6,*) 'ERROR: Cannot apply scan position bias correction -- error reading bcor file ', coef_in
    !  call abort()
    !endif
    !if ( TRIM(CF_INSTRUM) /= 'ATMS' ) then
    !  write(6,*) 'ERROR: Wrong instrument type in bcor file! Instrument = ', CF_INSTRUM
    !  call abort()
    !endif
    !write(6,*) ' Tb will be adjusted (internally) for scan-position bias dependency.'

    write(6,*) 'midas_satQCATMS: empty bcor file is used in the unitTest. Modifying before continuing.'
    call abort()
  endif
   
  if ( modlsqtt ) then
     write(6,*) 'MODLSQ option is activated!'
     write(6,*) '  Output file will contain recomputed values for land/sea qualifier and terrain type based on LG/MG.'
  endif

! initialisation
! --------------

  Call BURP_Init(File_in,  F2=File_out,  IOSTAT=error)
  Call BURP_Init(Rpt_in,   R2=Rpt_out,   IOSTAT=error)
  Call BURP_Init(Block_in, B2=Block_copy, IOSTAT=error)

! Set BURP "missing value" for reals

  opt_missing = 'MISSING'
  Call BURP_Set_Options(REAL_OPTNAME=opt_missing,REAL_OPTNAME_VALUE=zmisg)

! ouverture du fichier burp d'entree et de sortie
! -----------------------------------------------

  Call BURP_New(File_in,  FILENAME= brp_in,  MODE= FILE_ACC_READ,   IOSTAT= error)
  Call BURP_New(File_out, FILENAME= brp_out, MODE= FILE_ACC_CREATE, IOSTAT= error)

! Number of reports and maximum report size from input BURP file
! --------------------------------------------------------------

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
  nsize = nsize*2

! allouer l'espace
! ----------------

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

!--------------------------------------------------------------------------
!  Optional scan position adjustment of Tb: 
!--------------------------------------------------------------------------
!    -- Get global scan bias array glbscanb(nchan,nscan,nsat) from 
!       bias correction file
!        numan(numbsat) = number of channels (out)
!        listan = list of channel numbers (out)
!        csatid = list of satellites (dimension(3)) (out)
!    -- Compute Tb adjustments (dglbscanb)
!
   if (sp_adj_tb) then
!
        kk = mxscan/2
        dglbscanb(:,:,:) = 0.0
        csatid(:) = 'XXXXXXXXX'
        numbsat = NSAT
        numan(:) = 0
        do j = 1, mxan
          listan(j,:) = j
        end do
        do ii = 1,NSAT
          csatid(ii) = SATNAMES(ii)
          numan(ii)  = RCNCHAN(ii)
          if ( numan(ii) /= nchan ) then
            write(6,*) 'INFO:',csatid(ii),': Number of channels in coeff file is less than ', nchan
            write(6,*) '      Channels are ', CHANNUMS(ii,1:rcnchan(ii))
          endif
        ! Compute the scan position bias (SPB) for channels with bcor data
          do jj = 1, RCNCHAN(ii)
            j = CHANNUMS(ii,jj)
            dglbscanb(j,:,ii) =  FOVBIAS(ii,jj,:) - FOVBIAS(ii,jj,kk)
          end do
        end do
        if (debug) then
          write(6,*) 'Finished reading bcor coeff file.'
          write(6,*) 'Number of satellites in file = ', numbsat
          write(6,*) 'Number of channels per satellite = ', numan(1:numbsat)
          do kk=1,numbsat
            write(6,*) csatid(kk)
            write(6,*) '  Channels are ', CHANNUMS(kk,1:rcnchan(ii))
            write(6,*) '  Scan biases for each channel are:'
            do j = 1, nchan
              write(6,*) j, dglbscanb(j,:,kk)
            end do
          end do
        end if
        
   end if
!--------------------------------------------------------------------------
!  END optional scan position adjustment of Tb

!-------------------------------------------------------------------------------------
! LOOP OVER ALL REPORTS OF THE INPUT FILE, APPLY PROCESSING, AND WRITE TO OUTPUT FILE.
!-------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Initial scan of file to get number of reports and number of data locations.
! Store address of each report in array adresses(nb_rpts) for main REPORTS loop
!-----------------------------------------------------------------------------

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
      Call BURP_Free(Block_in,B2=Block_copy)
      STOP
  end if

!--------------------------------------------------------------------------
!! MAIN LOOP through all the reports in the file
!--------------------------------------------------------------------------

  nobs_tot = 0
  n_reps = 0
  n_bad_reps = 0
  n_reps_tb2misg = 0
  
!!!! ***DEBUG***  
!     compteur = 2  !!!! ***DEBUG***
!!!! ***DEBUG***

!!==============================================================================================================
  REPORTS: do k = 1, compteur
!  REPORTS: do k = 443, 443
!!==============================================================================================================

!write(*,*) 'Report ',k

    resume_report = .false.

    Call BURP_Get_Report(File_in, REPORT= Rpt_in, REF= adresses(k), IOSTAT= error) 
    if (error /= burp_noerr) call handle_error()

    Call BURP_Get_Property(Rpt_in,STNID=id,IDTYP=my_idtyp,ELEV=nlocs,LATI=blat,LONG=blon,NBLK=nblocs,HANDLE=handle)

    if ( id(1:2) .eq. ">>" ) then
       resume_report = .true.
       Call BURP_Write_Report(File_out,Rpt_in,IOSTAT=error)
       if (error /= burp_noerr)  call handle_error()
       CYCLE REPORTS
    else
! Create new report (Rpt_out) to contain modified blocks from Rpt_in
!write(*,*) 'nsize = ', nsize
!write(*,*) 'Before BURP_New(Rpt_out, Alloc_Space = nsize)'
      Call BURP_New(Rpt_out, Alloc_Space = nsize,  IOSTAT=error)
!write(*,*) 'After BURP_New(Rpt_out, Alloc_Space = nsize)'
      if (error /= burp_noerr)  call handle_error()
    
! initiliser pour ecriture a File_out
      Call BURP_INIT_Report_Write(File_out,Rpt_out,IOSTAT=error)
!write(*,*) 'After BURP_INIT_Report_Write'
      if (error /= burp_noerr)  call handle_error()
!  copier le header du rapport 
      Call BURP_Copy_Header(TO= Rpt_out, FROM= Rpt_in)
!write(*,*) 'After BURP_Copy_Header'

      stnid = id

    endif

  IF ( .not. resume_report ) THEN !----------------------------------------------------------------



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
     allocate( lflagchn(nt,nchan), stat=alloc_status(6) )
     allocate( rclw(nt),     stat=alloc_status(7) )
     allocate( riwv(nt),     stat=alloc_status(8) )
!     allocate( isecs(nt),    stat=alloc_status(9) )
     allocate( precipobs(nt), stat=alloc_status(10) )
     allocate( iber(nt),     stat=alloc_status(11) )
     allocate( ztb89(nt),    stat=alloc_status(12) )
     allocate( ztb150(nt),   stat=alloc_status(13) )
     allocate( scatl(nt),    stat=alloc_status(14) )
     allocate( scatw(nt),    stat=alloc_status(15) )
     allocate( ztb_amsub3(nt), stat=alloc_status(16) )
     allocate( ztb_amsub5(nt), stat=alloc_status(17) )
     allocate( zdi(nt),      stat=alloc_status(18) )
     allocate( err(nt),      stat=alloc_status(19) )
     allocate( ascatw(nt),   stat=alloc_status(20) )
     allocate( tb23(nt),     stat=alloc_status(21) )
     allocate( tb31(nt),     stat=alloc_status(22) )
     allocate( tb50(nt),     stat=alloc_status(23) )
     allocate( tb53(nt),     stat=alloc_status(24) )
     allocate( tb89(nt),     stat=alloc_status(25) )
     allocate( lsq(nt),      stat=alloc_status(26) )
     allocate( trn(nt),      stat=alloc_status(27) )
     allocate( scatec(nt),   stat=alloc_status(28) )
     allocate( scatbg(nt),   stat=alloc_status(29) )
     allocate( SeaIce(nt),   stat=alloc_status(30) )
     allocate( lqc(nt,nchan), stat=alloc_status(31) )
     
     if( any(alloc_status /= 0) ) then
       write(6,*) ' midas_satQCATMS: Memory allocation error '
       call abort()
     endif
    
     ident(:) = 0        ! filter information flag; set all bits OFF
     
! Information flag (ident) values (new BURP element 025174 in header)

!   BIT    Meaning
!    0     off=land or sea-ice, on=open water away from coast
!    1     Mean 183 Ghz [ch. 18-22] is missing
!    2     CLW is missing (over water)
!    3     CLW > clw_atms_nrl_LTrej (0.175 kg/m2) (cloudobs)
!    4     scatec/scatbg > Lower Troposphere limit 9/10 (precipobs)
!    5     Mean 183 Ghz [ch. 18-22] Tb < 240K
!    6     CLW > clw_atms_nrl_UTrej (0.200 kg/m2)
!    7     Dryness Index rejection (for ch. 22)
!    8     scatec/scatbg > Upper Troposphere limit 18/15
!    9     Dryness Index rejection (for ch. 21)
!   10     Sea ice > 0.55 detected
!   11     Gross error in Tb (any chan.)  (all channels rejected)


!write(*,*) 'call GetData'

!---------------------------------------------------------------------------
!  Get all the required data from the blocks in the report (Rpt_in)
!=====================================================================
      call GetData
!=====================================================================


! Initialize internal land/sea qualifier and terrain type arrays to values
! read from file

     lsq(:) = ilq(1:nt)  ! land/sea qualifier
     trn(:) = itt(1:nt)  ! terrain type (sea-ice)

!-----------------------------------------------------------------------------
!        *** Set WATEROBS() array and reset lsq, trn *** 
!-----------------------------------------------------------------------------
! Determine which obs pts are over open water (i.e NOT near coasts or
! over/near land/ice) using model MG and LG fields from glbhyb2 ANAL
!  MG = land/sea mask field (0.0 (water) to 1.0 (land))
!  LG = ice fraction field  (0.0 - 1.0)
!  lsq = 0 (land), 1 (water)
!  trn = -1 (no ice/snow),  0 (ice)
!  NOTE: land_ice_mask_atms redefines lsq and trn based on interpolated 
!  MG, LG fields so that
!    lsq = 1 (point over water away from land/coast), 0 (land/coast) otherwise
!    trn = 0 (point over or near sea-ice),           -1 (ice free) otherwise
!
!  waterobs(:)=.true. at points where lsq = 1 and trn = -1
!
!-----------------------------------------------------------------------------
!write(*,*) 'call land_ice_mask_atms'
     call land_ice_mask_atms(mglg_file,nt,zlat,zlon,lsq,trn,waterobs)

!-----------------------------------------------------------------------
! Check for values of TB that are missing or outside physical limits.
! **NOTE: REJECT ALL CHANNELS IF ONE IS FOUND TO BE BAD.
!-----------------------------------------------------------------------

     grossrej(:)  = .false.
     call gross_value_check(nt,ztb,grossrej)
     
     if ( ANY(grossrej) ) then
       write(6,*) ' GROSS_VALUE_CHECK has detected bad Tb data. Number of affected locations = ', COUNT(grossrej)
       write(6,*) '   Box lat, lon = ', blat, blon
     endif

!---------------------------------------------------------------------------
! Preliminary QC checks --> set lqc(nt,nchan)=.true. for data that fail QC
!---------------------------------------------------------------------------

     lqc(:,:) = .false.  ! Flag for preliminary QC checks
!write(*,*) 'call qc_data'
     call qc_data(zenith,ilq,itt,zlat,zlon,ztb,scanpos,stnid,nval,nt,lqc, &
     &            grossrej,lsq,trn,qcflag1,qcflag2,ican,blat,blon,lutb)

     if ( lutb ) n_reps_tb2misg = n_reps_tb2misg + 1


!      if ( k == 443 ) then
!        write (*,*) k, ' Number of points where ILQ not equal LSQ = ', COUNT(lsq(1:nt) /= ilq(1:nt))
!      endif



!! Output of qc_data
!!   lqc=true (entire array) --> problem with channels (number of channels or the channel numbers)
!!                               file could be corrupted [abort]
!!   ztb(nchan) = zmisg  (at locations with bad zenith angle and/or lat,lon)
!!   zenith     = zmisg  (if bad zenith angle)
!!   lutb=true  (if ztb set to zmisg at 1 or more locations)
!! All channels are flagged for all checks except the "data level" QC flag check.

!     if ( lutb ) then
!       write(6,*) ' Number of Tb data = zmisg after QC_DATA = ', COUNT(ztb == zmisg)
!       write(6,*) ' Number of Tb data = 330.04              = ', COUNT(ztb == 330.04)
!       write(6,*) ' Total number of data                    = ', nval*nt
!     endif

     bad_report = .false.
     if ( COUNT(lqc) == nt*nchan ) then
       write(6,*) ' QC_DATA has detected a problem with data in this report!'
       write(6,*) '   Report box lat, lon = ', blat, blon
       bad_report = .true.
       n_bad_reps = n_bad_reps + 1
     endif
     
   IF (.not. bad_report) THEN   !-------------------------------------------------------------------

!  Exclude problem points from further calculations

     do kk = 1,nt
       if ( COUNT(lqc(kk,:)) == nchan ) grossrej(kk) = .true.
     enddo

     where ( grossrej ) ident = IBSET(ident,11)

!---------------------------------------------------------------------------------------------------
! Apply NRL cloud filter, scattering index and sea-ice detection algorithms to 
!   OPEN WATER (waterobs=true) points.
! Points with SeaIce>0.55 are set to sea-ice points (waterobs --> false)
!---------------------------------------------------------------------------------------------------

! To begin, assume that all obs are good.

     cloudobs(:)  = .false.
     iwvreject(:) = .false.
     precipobs(:) = .false.

! -----------------------------------------------------------------------------
!    First remove scan-dependency of Tb biases for data in the report (OPTION)
! -----------------------------------------------------------------------------
     if (sp_adj_tb) then
   
        ii = 0
!      Find satellite index (satellite "^NPP")
        do kk = 1, numbsat
          if ( TRIM(csatid(kk)) == TRIM(stnid(2:9)) ) ii = kk
        end do
        if ( ii == 0 ) then
          write(*,*) ' Error: Satellite not found in bias correction file!'
          write(*,*) '        Satellite = ', stnid(2:4)
          write(*,*) '        Satellites in BCOR file = ', csatid(:)
          call abort()
        end if
!      Extract the Tb adjustments for each channel for this satellite
        do kk = 1, nchan
          zbcor(kk,:) = dglbscanb(kk,:,ii)
          if ( debug ) &
     &       write(6,*) 'Scan bias adjustments for channel ', kk, ' : ', zbcor(kk,:)
        end do
!      Adjust Tb for for each channel according to scan position
        indx1 = 1
        do kk = 1, nt   !  loop over NT locations in report
          indx2 = kk*nchan
          if ( debug ) then
            write(6,*) 'location, indx1, indx2 = ', kk, indx1, indx2
          end if
          do jj = 1, nchan
            ztbcor(indx1+jj-1) = ztb(indx1+jj-1) - zbcor(jj,scanpos(kk))
            if ( debug ) then
              write(6,*) 'scanpos, ztb index = ', scanpos(kk), indx1+jj-1
              write(6,*) 'channel, ztb, zbcor, ztbcor = ', &
      &           jj, ztb(indx1+jj-1), zbcor(jj,scanpos(kk)), ztbcor(indx1+jj-1)
            end if
          end do  
          indx1 = indx2 + 1
        end do
        
     else  ! no correction
     
        ztbcor(:) = ztb(:)
        
     end if
! --------------------------------------------------------------------
     
!!     extract required channels:
!!        23 Ghz = AMSU-A 1 = ATMS channel 1 
!!        31 Ghz = AMSU-A 2 = ATMS channel 2
!!        50 Ghz = AMSU-A 3 = ATMS channel 3
!!        53 Ghz = AMSU-A 5 = ATMS channel 6
!!        89 Ghz = AMSU-A15 = ATMS channel 16
!!       150 Ghz = AMSU-B 2 = ATMS channel 17
!
!   Extract Tb for channels 16 (AMSU-B 1) and 17 (AMSU-B 2) for Bennartz SI
!   Extract Tb for channels 22 (AMSU-B 3) and 18 (AMSU-B 5) for Dryness Index (DI)

     indx1 = 1
     do ii = 1, nt
       indx2 = ii*nchan
       tb23(ii)      = ztbcor(indx1)
       tb31(ii)      = ztbcor(indx1+1)
       tb50(ii)      = ztbcor(indx1+2)
       tb53(ii)      = ztbcor(indx1+5)
       tb89(ii)      = ztbcor(indx1+15)
       ztb89(ii)     = tb89(ii)
       ztb150(ii)    = ztbcor(indx1+16)
       ztb_amsub3(ii) = ztbcor(indx1+21)
       ztb_amsub5(ii) = ztbcor(indx1+17)
       indx1 = indx2 + 1
     end do

!  nrl_filter returns rclw, scatec, scatbg and also does sea-ice detection
!  Missing value for  rclw, scatec, scatbg  is -99.0 (e.g. over land or sea-ice).
!  Sets trn=0 (sea ice) for points where retrieved SeaIce>=0.55.
!  Does nothing if trn=0 (sea ice) and retrieved SeaIce<0.55.
!

!write(*,*) 'call nrl_filter'
     call nrl_filter(err,nt,tb23,tb31,tb50,tb89,ztb150,zenith,zlat,lsq,trn, &
    &                waterobs,grossrej,rclw,scatec,scatbg,iNumSeaIce,iRej,SeaIce)
        
!---------------------------------------------------------------------------
!  Flag data using NRL criteria
!---------------------------------------------------------------------------
!
!    Compute Mean 183 Ghz [ch. 18-22] Tb (riwv)
!
        riwv = -99.0
        indx1 = 1
        do ii = 1, nt
          indx2 = ii*nchan
          if (.not.grossrej(ii)) then
            ztb183(1) = ztbcor(indx1+17)
            ztb183(2) = ztbcor(indx1+18)
            ztb183(3) = ztbcor(indx1+19)
            ztb183(4) = ztbcor(indx1+20)
            ztb183(5) = ztbcor(indx1+21)
            riwv(ii)  = sum(ztb183)/5.0
            if ( riwv(ii) < mean_Tb_183Ghz_min ) iwvreject(ii) = .true.
          else
            iwvreject(ii) = .true.
          endif
          indx1 = indx2 + 1
        end do
!
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
!        
        where ( waterobs .and. (rclw == -99.) ) ident = IBSET(ident,2)
        where ( riwv == -99.)                   ident = IBSET(ident,1)

!----------------------------------------------------------------------------------
! Compute the simple AMSU-B Dryness Index zdi for all points = Tb(ch.3)-Tb(ch.5)
!----------------------------------------------------------------------------------
    
     where ( .not.grossrej )
       zdi = ztb_amsub3 - ztb_amsub5
     elsewhere
       zdi = zmisg
     end where

!-----------------------------------------------------------------------------------
! Review all the checks previously made to determine which obs are to be accepted
! for assimilation and which are to be flagged for exclusion (lflagchn). 
!   grossrej()  = .true. if any channel had a gross error at the point
!   cloudobs()  = .true. if CLW > clw_atms_nrl_LTrej (0.175) or precipobs
!   precipobs() = .true. if precip. detected through NRL scattering indices
!   waterobs()  = .true. if open water point
!   iwvreject() = .true. if Mean 183 Ghz [ch. 18-22] Tb < 240K (too dry for ch.20-22 over land)
!-----------------------------------------------------------------------------------

     lflagchn(:,:) = lqc(:,:)  ! initialize with flags set in qc_data
     do kk = 1, nt
     
!    Reject all channels if gross Tb error detected in any channel or other problems 
       if ( grossrej(kk) ) then
    
         lflagchn(kk,:) = .true.

       else
!-----------------------------------------------------------------------------------       
!      OVER LAND OR SEA-ICE,
!         -- CLW/SI not determined over land
!         -- surface emissivity effects lower tropospheric and window channels     
!         -- reject window & lower tropospheric channels 1-6, 16-19
!         -- reject ch. 20-22 if iwvreject = .true.  [ Mean 183 Ghz [ch. 18-22] Tb < 240K ]
!         -- check DI for AMSU-B like channels
!----------------------------------------------------------------------------------- 
         if  ( .not. waterobs(kk) ) then
           lflagchn(kk,1:ipc)     = .true.      ! AMSU-A 1-6
           lflagchn(kk,16:19)     = .true.      ! AMSU-B (like 1,2,5)
           if ( iwvreject(kk) ) lflagchn(kk,20:22) = .true.  ! AMSU-B (like 4,3)
!    Dryness index (for AMSU-B channels 19-22 assimilated over land/sea-ice)
!!   Channel AMSUB-3 (ATMS channel 22) is rejected for a dryness index >    0.
!!                   (ATMS channel 21) is rejected for a dryness index >   -5.
!!   Channel AMSUB-4 (ATMS channel 20) is rejected for a dryness index >   -8.
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
         endif  ! if .not. waterobs
!-----------------------------------------------------------------------------------         
!      OVER WATER,
!         -- reject ch. 1-6, 16-20 if CLW > clw_atms_nrl_LTrej or CLW = -99.0
!         -- reject ch. 7-9, 21-22 if CLW > clw_atms_nrl_UTrej or CLW = -99.0
!         -- reject ch. 1-6, 16-22 if scatec > 9  or scatec = -99.0
!         -- reject ch. 7-9        if scatec > 18 or scatec = -99.0
!         -- reject ch. 1-6        if scatbg > 10 or scatbg = -99.0
!         -- reject ch. 7-9        if scatbg > 15 or scatbg = -99.0
!         -- reject ch. 16-22      if iwvreject = .true.   [ Mean 183 Ghz [ch. 18-22] Tb < 240K ]
!-----------------------------------------------------------------------------------    
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
         endif   ! if waterobs 
         
       endif   ! if ( grossrej(kk) [end else] )

       if ( .not. waterobs(kk) ) landcnt  = landcnt  + 1
       if ( grossrej(kk) )  rejcnt = rejcnt + 1
       if ( iwvreject(kk))  iwvcnt = iwvcnt + 1
       if ( precipobs(kk) .and. waterobs(kk) ) then
         pcpcnt = pcpcnt + 1
       endif
       
       if ( ANY(lflagchn(kk,:)) ) flgcnt = flgcnt + 1

     end do  ! kk = 1, nt loop

!!  RESET riwv array to ECMWF scattering index for output to BURP file

     riwv(:) = scatec(:)

!!  Set missing rclw and riwv to BURP missing value (zmisg)

     where (rclw == -99. ) rclw = zmisg
     where (riwv == -99. ) riwv = zmisg

!--------------------------------------------------------------------
!       Modify the blocks in Rpt_in and write to Rpt_out
! - Modify flag values so that the obs identified above as being over land/ice,
!   or in cloudy/precip regions, etc. are not assimilated (FLAG block 15392/15408).
! - OPTION: update land-sea qualifier and terrain type in INFO block 3072.
! - Update Tb data in DATA block 9248/9264 (if Tb was modified).
! - Add new elements to INFO block 3072.
! - Modify 24bit global flags in 3D block 5120 (if any data rejected).
!=====================================================================
!write(*,*) 'call WriteBlocks'
     call WriteBlocks
!=====================================================================

   ENDIF  !------------------------ if .not. bad_report ----------------------------------------



! Deallocate arrays.

     alloc_status(:) = 0
     deallocate( ident,    stat=alloc_status(1) )
     deallocate( waterobs, stat=alloc_status(2) )
     deallocate( grossrej, stat=alloc_status(3) )
     deallocate( cloudobs, stat=alloc_status(4) )
     deallocate( iwvreject,stat=alloc_status(5) )
     deallocate( lflagchn, stat=alloc_status(6) )
     deallocate( rclw,     stat=alloc_status(7) )
     deallocate( riwv,     stat=alloc_status(8) )
!     deallocate( isecs,    stat=alloc_status(9) )
     deallocate( precipobs, stat=alloc_status(10) )
     deallocate( iber,     stat=alloc_status(11) )
     deallocate( ztb89,    stat=alloc_status(12) )
     deallocate( ztb150,   stat=alloc_status(13) )
     deallocate( scatl,    stat=alloc_status(14) )
     deallocate( scatw,    stat=alloc_status(15) )
     deallocate( ztb_amsub3, stat=alloc_status(16) )
     deallocate( ztb_amsub5, stat=alloc_status(17) )
     deallocate( zdi,      stat=alloc_status(18) )
     deallocate( err,      stat=alloc_status(19) )
     deallocate( ascatw,   stat=alloc_status(20) )
     deallocate( tb23,     stat=alloc_status(21) )
     deallocate( tb31,     stat=alloc_status(22) )
     deallocate( tb50,     stat=alloc_status(23) )
     deallocate( tb53,     stat=alloc_status(24) )
     deallocate( tb89,     stat=alloc_status(25) )
     deallocate( lsq,      stat=alloc_status(26) )
     deallocate( trn,      stat=alloc_status(27) )
     deallocate( scatec,   stat=alloc_status(28) )
     deallocate( scatbg,   stat=alloc_status(29) )
     deallocate( SeaIce,   stat=alloc_status(30) )
     deallocate( lqc,      stat=alloc_status(31) )

     if( any(alloc_status /= 0) ) then
       write(6,*) ' midas_satQCATMS: Memory deallocation error '
       call abort()
     endif

  ENDIF  !------------------------ if .not.resume_report ---------------------------------------

!------------------------------------------------------------------------
! Write the modified report to the output file
!------------------------------------------------------------------------
IF ( .not. resume_report ) THEN
!write(*,*) ' Call BURP_Write_Report...'
  if (.not. bad_report ) then
    Call BURP_Write_Report(File_out,Rpt_out,IOSTAT=error)
    if (error /= burp_noerr)  call handle_error()
  endif
  Call BURP_Free(Rpt_out,IOSTAT=error)
  if (error /= burp_noerr)  call handle_error()
ENDIF

! proceed to the next report
!!==============================================================================================================
   end do REPORTS
!!==============================================================================================================

   write(6,*) ' --------------------------------------------------------------- '
   write(6,*) ' Number of obs pts read from BURP file              = ', nobs_tot
   write(6,*) ' Number of BURP file reports                        = ', n_reps
   write(6,*) ' Number of bad BURP file reports (all data flagged) = ', n_bad_reps
   write(6,*) ' Number of BURP file reports where Tb set to zmisg  = ', n_reps_tb2misg
   write(6,*) ' --------------------------------------------------------------- '
   write(6,*) ' 1. Number of obs pts found over land/ice           = ', landcnt
   write(6,*) ' 2. Number of problem obs pts (Tb err, QCfail)      = ', rejcnt
   write(6,*) ' 3. Number of cloudy obs  (CLW > clw_min)           = ', cldcnt
   write(6,*) ' 4. Number of scatter/precip obs                    = ', pcpcnt
   write(6,*) ' 5. Number of pts with Mean 183 Ghz Tb < 240K       = ', iwvcnt
   write(6,*) ' 6. Number of pts flagged for AMSU-B Dryness Index  = ', drycnt
   write(6,*) ' --------------------------------------------------------------- '
   write(6,*) ' Total number of filtered obs pts                   = ', flgcnt
   write(6,*) ' ----------------------------------------------------------------'
   write(6,*) ' '
   write(6,*) ' -------------------------------------------------------------------------------'
   write(6,*) ' Number of waterobs points converted to sea ice points         = ', iNumSeaIce
   write(6,*) ' Number of points where CLW/SI missing over water due bad data = ', iRej
   write(6,*) ' -------------------------------------------------------------------------------'
   write(6,*) ' '
   write(6,*) ' '
   write(6,*) '   Meaning of IDENT flag bits: '
   write(6,*) ' '
   write(6,*) '      BIT    Meaning'
   write(6,*) '       0     off=land or sea-ice, on=open water away from coast'
   write(6,*) '       1     Mean 183 Ghz [ch. 18-22] is missing'
   write(6,*) '       2     NRL CLW is missing (over water)'
   write(6,*) '       3     NRL > clw_atms_nrl_LTrej (0.175 kg/m2) (cloudobs)'
   write(6,*) '       4     scatec/scatbg > Lower Troposphere limit 9/10 (precipobs)'
   write(6,*) '       5     Mean 183 Ghz [ch. 18-22] Tb < 240K'
   write(6,*) '       6     CLW > clw_atms_nrl_UTrej (0.200 kg/m2)'
   write(6,*) '       7     Dryness Index rejection (for ch. 22)'
   write(6,*) '       8     scatec/scatbg > Upper Troposphere limit 18/15'
   write(6,*) '       9     Dryness Index rejection (for ch. 21)'
   write(6,*) '      10     Sea ice > 0.55 detected'
   write(6,*) '      11     Gross error in Tb (any chan.) or other QC problem (all channels rejected)'
   write(6,*) ' '
   write(6,*) '   New Element 13209 in BURP file = CLW (kg/m2)'
   write(6,*) '   New Element 13208 in BURP file = ECMWF Scattering Index'
   write(6,*) '   New Element 25174 in BURP file = IDENT flag'
   write(6,*) ' '


   Deallocate(adresses)

   Call BURP_Free(File_in,F2=File_out,IOSTAT=error)
   Call BURP_Free(Rpt_in,R2=Rpt_out,IOSTAT=error)
   Call BURP_Free(Block_in,B2=Block_copy,IOSTAT=error)

!--------------------------------------------------------------------
   istat=exfin('midas_satQCATMS',VERSION,'NON')
!--------------------------------------------------------------------

contains

subroutine GetData
!--------------------------------------------------------------------------------------
! Object:   This routine extracts the needed data from the blocks in the report Rpt_in:
!             zenith(nt)       = satellite zenith angle (btyp=3072,ele=7024)
!             ilq(nt)          = land/sea qualifier     (btyp=3072,ele=8012)
!             itt(nt)          = terrain-type (ice)     (btyp=3072,ele=13039)
!             zlat(nt)                                  (btyp=5120,ele=5002)
!             zlon(nt)                                  (btyp=5120,ele=6002)
!             ztb(nt*nval),    = brightness temperature (btyp=9248/9264,ele=12163) 
!             scanpos(nt),     = scan position (fov)    (btyp=3072,ele=5043) 
!             nval,            = number of channels     (btyp=9248/9264)
!             nt,              = number of locations    (btyp=5120,etc.)
!             qcflag1(nt,3)    = flag values for btyp=3072 block ele 033078, 033079, 033080
!             qcflag2(nt*nval) = flag values for btyp=9248 block ele 033081
!             ican(nt*nval)    = channel numbers btyp=9248 block ele 5042 (= 1-22)

! NOTE:  k = report number (from MAIN program) **** DO NOT MODIFY ****
!       kk = variable for loops over locations (nt)
!        j = variable for loops over nval (nval = 1 or nchan)


! Internal variables
   integer :: ipos

!------------------------------------------------------------------------------------
!  (1) Get the lat,lon from time/location block    BTYP = 5120  (also get nt)
!------------------------------------------------------------------------------------
   ref_blk = 0
   ref_blk = BURP_Find_Block(Rpt_in, &
                   BLOCK       = Block_in, &
                   SEARCH_FROM = ref_blk, &
                   BFAM        = 0, &
                   BTYP        = 5120, &
                   IOSTAT      = error)
   if (error /= burp_noerr) call handle_error()
   if (ref_blk < 0) then
     write(*,*) 'ERREUR -  Location/time (3D) block (btyp=5120) not found in report number ', k
     call abort()
   endif

   call BURP_Get_Property(Block_in, &
                     NELE = my_nele, &
                     NT   = my_nt, &    ! number of locations in the box (report)
                     NVAL = my_nval, IOSTAT=error)
   if (error /= burp_noerr)  call handle_error()
   
   nt = my_nt   ! set nt for MAIN program

   ind_lat = BURP_Find_Element(Block_in,5002,IOSTAT=error)
   ind_lon = BURP_Find_Element(Block_in,6002,IOSTAT=error)
   
   if ( (ind_lat > 0) .AND. (ind_lon > 0) ) then
     j = 1
     do kk =1, my_nt
       zlat(kk) = BURP_Get_Rval(Block_in,ind_lat,j,kk,error)
       zlon(kk) = BURP_Get_Rval(Block_in,ind_lon,j,kk,error)
     end do
   else
     write(*,*) 'ERREUR - lat, lon elements (5002,6002) missing in 3D block (btyp=5120). Report = ', k
     call abort()
   endif

!------------------------------------------------------------------------------------
! (2) Get info elements from the INFO block   BTYP = 3072
!------------------------------------------------------------------------------------
   ref_blk = 0
   ref_blk = BURP_Find_Block(Rpt_in, &
                   BLOCK       = Block_in, &
                   SEARCH_FROM = ref_blk, &
                   BFAM        = 0, &
                   BTYP        = 3072, &
                   IOSTAT      = error)
   if (error /= burp_noerr) call handle_error()
   if (ref_blk < 0) then
     write(*,*) 'ERREUR - INFO block (btyp=3072) not found in report number ', k
     call abort()
   endif

   call BURP_Get_Property(Block_in, &
                     NELE = my_nele, &
                     NT   = my_nt, &    ! number of locations in the box (report)
                     NVAL = my_nval, IOSTAT=error)
   if (error /= burp_noerr)  call handle_error()

!--------------------------------------------------------------------------------
   indice = BURP_Find_Element(Block_in,7024,IOSTAT=error)
   if ( indice > 0 ) then
     j = 1
     do kk =1, my_nt
       zenith(kk) = BURP_Get_Rval(Block_in,indice,j,kk,error)
      end do
   else
     write(*,*) 'ERREUR - satellite zenith angle missing in INFO block (ele=7024). Report = ', k
     call abort()
   endif
!--------------------------------------------------------------------------------   
   indice = BURP_Find_Element(Block_in,8012,IOSTAT=error)
   if ( indice > 0 ) then
     j = 1
     do kk =1, my_nt
       ilq(kk) = BURP_Get_Tblval(Block_in,indice,j,kk,error)
      end do
   else
     write(*,*) 'ERREUR - land/sea qualifier missing in INFO block (ele=8012). Report = ', k
     call abort()
   endif   
!--------------------------------------------------------------------------------
   indice = BURP_Find_Element(Block_in,13039,IOSTAT=error)
   if ( indice > 0 ) then
     j = 1
     do kk =1, my_nt
       itt(kk) = BURP_Get_Tblval(Block_in,indice,j,kk,error)
      end do
   else
     write(*,*) 'ERREUR - terrain-type missing in INFO block (ele=13039). Report = ', k
     call abort()
   endif   
!--------------------------------------------------------------------------------
   indice = BURP_Find_Element(Block_in,5043,IOSTAT=error)
   if ( indice > 0 ) then
     j = 1
     do kk =1, my_nt
       scanpos(kk) = BURP_Get_Tblval(Block_in,indice,j,kk,error)
      end do
   else
     write(*,*) 'ERREUR - scan position missing in INFO block (ele=5043). Report = ', k
     call abort()
   endif   
!--------------------------------------------------------------------------------
   indice = BURP_Find_Element(Block_in,33078,IOSTAT=error)
   if ( indice > 0 ) then
     j = 1
     do kk =1, my_nt
        qcflag1(kk,1) = BURP_Get_Tblval(Block_in,indice,j,kk,error)
      end do
   else
     write(*,*) 'ERREUR - Geolocation quality code missing in INFO block (ele=33078). Report = ', k
     call abort()
   endif  
   indice = BURP_Find_Element(Block_in,33079,IOSTAT=error)
   if ( indice > 0 ) then
     j = 1
     do kk =1, my_nt
        qcflag1(kk,2) = BURP_Get_Tblval(Block_in,indice,j,kk,error)
      end do
   else
     write(*,*) 'ERREUR - Granule level QC flag missing in INFO block (ele=33079). Report = ', k
     call abort()
   endif  
   indice = BURP_Find_Element(Block_in,33080,IOSTAT=error)
   if ( indice > 0 ) then
     j = 1
     do kk =1, my_nt
        qcflag1(kk,3) = BURP_Get_Tblval(Block_in,indice,j,kk,error)
      end do
   else
     write(*,*) 'ERREUR - Scan level QC flag missing in INFO block (ele=33080). Report = ', k
     call abort()
   endif  

!------------------------------------------------------------------------------------
! (3) Get data from the DATA block     BTYP = 9248 or 9264    (also get nval = nchan)
!------------------------------------------------------------------------------------
   ref_blk = 0
   ref_blk = BURP_Find_Block(Rpt_in, &
                   BLOCK       = Block_in, &
                   SEARCH_FROM = ref_blk, &
                   BFAM        = 0, &
                   BTYP        = 9248, &
                   IOSTAT      = error)
   if (error /= burp_noerr) call handle_error()
   if (ref_blk < 0) then
     ref_blk = 0
     ref_blk = BURP_Find_Block(Rpt_in, &
                   BLOCK       = Block_in, &
                   SEARCH_FROM = ref_blk, &
                   BFAM        = 0, &
                   BTYP        = 9264, &
                   IOSTAT      = error)
     if (ref_blk < 0) then
       write(*,*) 'ERREUR - DATA block (btyp 9248 or 9264) not found in report number ', k
       call abort()
     endif
   endif

   call BURP_Get_Property(Block_in, &
                     NELE = my_nele, &
                     NT   = my_nt, &    ! number of locations in the box (report)
                     NVAL = my_nval, IOSTAT=error)
   if (error /= burp_noerr)  call handle_error()

   nval = my_nval    ! set nval (#channels) for MAIN program

!--------------------------------------------------------------------------------   
   indice1 = BURP_Find_Element(Block_in,12163,IOSTAT=error)
   indice2 = BURP_Find_Element(Block_in,33081,IOSTAT=error)
   indice3 = BURP_Find_Element(Block_in,5042,IOSTAT=error)
   
   if ( indice1 > 0 .and. indice2 > 0 .and. indice3 > 0 ) then
     ipos = 0
     do kk = 1, my_nt
       do j = 1, my_nval
         ipos = ipos + 1
         ztb(ipos)     = BURP_Get_Rval(Block_in,indice1,j,kk,error)
         qcflag2(ipos) = BURP_Get_Tblval(Block_in,indice2,j,kk,error)
         ican(ipos)    = BURP_Get_Tblval(Block_in,indice3,j,kk,error)
       end do
     end do
   else
     write(*,*) 'ERREUR - Elements are missing in DATA block. Report = ', k
     if ( indice1 <= 0 )  write(*,*) '       Tb data             (012163) are missing!'
     if ( indice2 <= 0 )  write(*,*) '       Data level QC flags (033081) are missing!'
     if ( indice3 <= 0 )  write(*,*) '       Channel numbers     (005042) are missing!'
     call abort()
   endif 

   return

end subroutine GetData

!--------------------------------------------------------------------


subroutine WriteBlocks
!--------------------------------------------------------------------------------------
! Object:   This routine modifies the blocks in the input Report (Rpt_in) 
!--------------------------------------------------------------------------------------

   integer(kind=int_def)    :: iidata
   integer                  :: ipos

!------------------------------------------------------------------------------
! (1) Read and modify the blocks in Rpt_in and add them to Rpt_out
!------------------------------------------------------------------------------
    ref_blk = 0

!!==========================================
  BLOCKS: do
!!==========================================  

    ref_blk = BURP_Find_Block(Rpt_in,BLOCK= Block_in,SEARCH_FROM= ref_blk,IOSTAT= error)
    if (error /= burp_noerr) call handle_error()
    
    if (ref_blk < 0) Exit

    Call BURP_Get_Property(Block_in, &
                NELE   = nbele, &
                NVAL   = nvale, &       ! 1 or number of channels (obs per location) if Tb data/flag block
                NT     = nte, &         ! 1 or number of locations in block
                BTYP   = btyp, &
                BFAM   = bfam, &
                IOSTAT = error)
    if (error /= burp_noerr) call handle_error()
    
    if (btyp == 5120) then     !  -------- 3D Block ---------------
    ! Set bit 6 in 24-bit global flags if any data rejected

!write(*,*) 'BLOCKS:btyp == 5120'

    ! Extract the global flags, element 55200
    
!    Call BURP_TO_STDOUT(Block_in, CONVERT =.FALSE.)
      
      indice = BURP_Find_Element(Block_in,55200,IOSTAT=error)
      
      if ( indice > 0 ) then
        j = 1
        do kk =1, nte
           iidata = BURP_Get_Tblval(Block_in,indice,j,kk,error)
           if ( ANY(lflagchn(kk,:)) ) iidata = IBSET(iidata,6)
           Call BURP_Set_Tblval(Block_in,indice,j,kk,iidata)
        end do
      else
        write(*,*) 'ERREUR - Global flag missing in 3D block (ele=55200). Report = ', k
        call abort()
      endif 
      
!      Call BURP_TO_STDOUT(Block_in, CONVERT =.FALSE.)

      Block_copy = Block_in
!      write(*,*) 'Before Call BURP_Write_Block'
      Call BURP_Write_Block(Rpt_out,BLOCK=Block_copy,CONVERT_BLOCK=.FALSE.,IOSTAT=error)
      if (error /= burp_noerr)  call handle_error()
!      write(*,*) 'After Call BURP_Write_Block'
    
    elseif (btyp == 3072) then !  -------- INFO block (mix of integer and real data) ---------------

!write(*,*) 'BLOCKS:btyp == 3072'

    ! Add new elements QC indent flag (ident), CLW (rclw) and ECMWF_SI (riwv)
    !   OPTION: replace land-sea qualifier (lsq) and terrain type (trn) with internal values
    !   NOTE: after setting real values (Rval), call BURP_Convert_Block()!!

      Call BURP_Resize_Block(Block_in, ADD_NELE = 3, IOSTAT = error)
      if (error /= burp_noerr)  call handle_error()
      
      Call BURP_Set_Element(Block_in, NELE_IND = nbele+1, ELEMENT = 25174, IOSTAT = error)
      Call BURP_Set_Element(Block_in, NELE_IND = nbele+2, ELEMENT = 13209, IOSTAT = error)
      Call BURP_Set_Element(Block_in, NELE_IND = nbele+3, ELEMENT = 13208, IOSTAT = error)
      Call BURP_Encode_Block(Block_in)   ! encode the element numbers in the block

! Update the REAL elements first

      j = 1
      do kk =1, nte
        iidata = ident(kk)
        Call BURP_Set_Rval  (Block_in, NELE_IND=nbele+2,NVAL_IND=j,NT_IND=kk, RVAL=rclw(kk),IOSTAT=error)
        Call BURP_Set_Rval  (Block_in, NELE_IND=nbele+3,NVAL_IND=j,NT_IND=kk, RVAL=riwv(kk),IOSTAT=error)
      end do
      
      Call BURP_Convert_Block(Block_in)

! Now update the INTEGER elements

      j = 1
      do kk =1, nte
        iidata = ident(kk)
        Call BURP_Set_Tblval(Block_in, NELE_IND=nbele+1,NVAL_IND=j,NT_IND=kk, TBLVAL=iidata,IOSTAT= error)
      end do

      if (modlsqtt) then
         indice1 = BURP_Find_Element(Block_in,  8012, IOSTAT=error)
         if (error /= burp_noerr)  call handle_error()
         indice2 = BURP_Find_Element(Block_in, 13039, IOSTAT=error)
         if (error /= burp_noerr)  call handle_error()
         if ( indice1 > 0 .and. indice2 > 0 ) then
           j = 1
           do kk =1, nte
             iidata = lsq(kk)
             Call BURP_Set_Tblval(Block_in,indice1,j,kk,iidata,error)
             if (error /= burp_noerr)  call handle_error()
             iidata = trn(kk)
             Call BURP_Set_Tblval(Block_in,indice2,j,kk,iidata,error)
             if (error /= burp_noerr)  call handle_error()
           end do
         else
           write(*,*) 'ERREUR - land/sea qualifier (ele=8012) and/or terrain type (ele=13039) not found in INFO block. Report = ', k
           call abort()
         endif                     
      endif

!      if ( debug ) Call BURP_TO_STDOUT(Block_in, CONVERT =.FALSE.)

      Block_copy = Block_in
      Call BURP_Write_Block(Rpt_out,BLOCK=Block_copy,CONVERT_BLOCK=.FALSE.,IOSTAT=error)
      if (error /= burp_noerr)  call handle_error()      

    
    elseif (btyp == 9248 .or. btyp ==9264) then !  -------- DATA block ---------------
    ! Modify Tb data if any data (ztb) were set to zmisg (lutb=.true.)

!write(*,*) 'BLOCKS:btyp == 9248'

       if (lutb) then
          indice = BURP_Find_Element(Block_in, 12163, IOSTAT=error)
          if (error /= burp_noerr)  call handle_error()
          if ( indice > 0 ) then
            ipos = 0
            do kk =1, nte
              do j = 1, nvale
                ipos = ipos + 1
                Call BURP_Set_Rval(Block_in,NELE_IND=indice,NVAL_IND=j,NT_IND=kk,RVAL=ztb(ipos),IOSTAT=error)  
                if (error /= burp_noerr)  call handle_error()
              enddo
            enddo
          else
            write(*,*) 'ERREUR - Cannot find Tb (ele=12163) in DATA block!. Report = ', k
            call abort()
          endif
!          if ( debug ) Call BURP_TO_STDOUT(Block_in, CONVERT =.FALSE.)
          Block_copy = Block_in
          Call BURP_Write_Block(Rpt_out,BLOCK=Block_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
          if (error /= burp_noerr)  call handle_error() 
       else
         Call BURP_Write_Block(Rpt_out,BLOCK=Block_in,CONVERT_BLOCK=.FALSE.,IOSTAT=error)
         if (error /= burp_noerr)  call handle_error() 
       endif
    
    elseif (btyp == 15392 .or. btyp == 15408) then !  -------- FLAG block ---------------
   ! Modify data flag values (set bit 7) for rejected data    

!write(*,*) 'BLOCKS:btyp == 15392'

      indice = BURP_Find_Element(Block_in, 212163, IOSTAT=error)
      if (error /= burp_noerr)  call handle_error()
      if ( indice > 0 ) then 
        do kk =1, nte
          do j = 1, nvale
            iidata = BURP_Get_Tblval(Block_in,indice,j,kk,error)
            if (lflagchn(kk,j)) iidata = IBSET(iidata,7)
            Call BURP_Set_Tblval(Block_in,indice,j,kk,iidata)
          enddo
        enddo
      else
        write(*,*) 'ERREUR - Data QC flags (ele=212163) not found in FLAG block. Report = ', k
        call abort()      
      endif

      Block_copy = Block_in
      Call BURP_Write_Block(Rpt_out,BLOCK=Block_copy,CONVERT_BLOCK=.FALSE.,IOSTAT=error)
      if (error /= burp_noerr)  call handle_error()

    else     ! ------------ OTHER BLOCK --------------
    
      Call BURP_Write_Block(Rpt_out,BLOCK=Block_in,CONVERT_BLOCK=.FALSE.,IOSTAT=error)
      if (error /= burp_noerr)  call handle_error() 
      
    endif      
  
!!==========================================  
  enddo BLOCKS
!!==========================================  

  return

  end subroutine WriteBlocks

!--------------------------------------------------------------------

subroutine land_ice_mask_atms(mglg_file,npts,zlat,zlon,zlq,ztt,waterobs)
!--------------------------------------------------------------------------------------
!
! Author:   S. Macpherson  MSC/ARMA
! Language: FORTRAN 90
! Adapted from: land_ice_mask_ssmis.ftn90 of satqc_ssmis (D. Anselmo, S. Macpherson)
!
! Object:   This routine sets waterobs array by performing a land/ice proximity check using
!           using analysis MG and LG (or GL) fields used by the model which produces the trial field.
!           The purpose of this check is to remove obs that reside close to coasts or ice,
!           and so whose TBs may be contaminated.
!           The GEM Global (glbhyb2) analysis contains MG and LG fields (on different grids).
!
!           NOTE: The 0.1 deg binary ice field check from land_ice_mask_ssmis.ftn90
!           was removed. The land/sea qualifier (zlq) and terrain type (ztt) are modified
!           to indicate proximity to land and sea-ice but are NOT changed in output BURP file.
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
!   0.1       16/08/12     Original adapted code.      S. Macpherson  
!   0.2       01/03/14     Open mglg_file in R/O mode  S. Macpherson
!
!--------------------------------------------------------------------
!  Variable Definitions
!  --------------------
! mglg_file  - input  -  name of file holding model MG and LG (or GL) fields
! npts       - input  -  number of input obs pts in report
! zlat       - input  -  array holding lat values for all obs pts in report
! zlon       - input  -  array holding lon values for all obs pts in report
! zlq        - in/out -  array holding land/sea qualifier values for all obs
!                        pts of report (0 = land, 1 = sea)
! ztt        - in/out -  array holding terrain-type values for all obs pts
!                        of current report (-1 land/open water, 0 = ice)
! waterobs   - output -  logical array identifying for each obs in current report
!                        whether it is over open water, far from coast/ice
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
! MGthresh   -internal-  maximum allowable land fraction for obs to be kept
! LGthresh   -internal-  maximum allowable ice  fraction for obs to be kept
!--------------------------------------------------------------------
!  use var_declare
  implicit none

!  Arguments:
  character(len=128), intent(in) :: mglg_file

  integer, intent(in)                   :: npts
  real,    intent(in),     dimension(:) :: zlat,zlon
  integer, intent(inout),  dimension(:) :: zlq, ztt

  logical, intent(out), dimension(:) :: waterobs

!  Locals:
  integer, parameter :: mxlat=5,mxlon=5
  integer, parameter :: iungeo=50

  integer :: ier,key,istat
  integer :: ni,nj,nk,nilg,njlg
  integer :: ig1,ig2,ig3,ig4,ig1lg,ig2lg,ig3lg,ig4lg
  integer :: idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11
  integer :: idum12,idum13,idum14,idum15,idum16,idum17,idum18

  integer :: indx,ii,jj,kk
  integer :: nlat,nlon

  integer, dimension(10) :: alloc_status
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
  logical  :: llg

! F90 allocatable arrays:
  real, allocatable, dimension(:)   :: mg,lg
  real, allocatable, dimension(:)   :: latmesh,lonmesh
  real, allocatable, dimension(:)   :: mgintob,lgintob
  real, allocatable, dimension(:,:) :: zlatbox,zlonbox
  real, allocatable, dimension(:)   :: mgintrp,lgintrp
!
! RMNLIB interpolating functions:
  integer :: ezsetopt,ezqkdef
  integer :: gdllsval,gdid,gdidlg
!--------------------------------------------------------------------

! Allocate space for arrays holding values on mesh grid pts.

  alloc_status(:) = 0
  allocate ( latmesh(mxlat*mxlon), stat=alloc_status(1) )
  allocate ( lonmesh(mxlat*mxlon), stat=alloc_status(2) )
  allocate ( mgintob(mxlat*mxlon), stat=alloc_status(3) )
  allocate ( lgintob(mxlat*mxlon), stat=alloc_status(4) )
  allocate ( zlatbox(mxlat*mxlon,npts), stat=alloc_status(5) )
  allocate ( zlonbox(mxlat*mxlon,npts), stat=alloc_status(6) )
  if( any(alloc_status /= 0) ) then
    write(6,*) ' LAND_ICE_MASK_ATMS: Memory allocation error '
    call abort()
  endif

! Open FST file.

  ier = fnom( iungeo,mglg_file,'STD+RND+R/O',0 )
  ier = fstouv( iungeo,'RND' )

! Read MG field.

  key = fstinf(iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ' ,'MG')
  if ( key <  0 ) then
    write(6,*) ' LAND_ICE_MASK_ATMS: The MG field is MISSING '
    call abort()
  end if

  allocate ( mg(ni*nj), stat=alloc_status(7) )
  if( any(alloc_status /= 0) ) then
    write(6,*) ' LAND_ICE_MASK_ATMS: Memory allocation error '
    call abort()
  endif
  ier = fstlir(mg,iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ','MG')

  ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,    &
      &        idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyp,ig1,ig2,  &
      &        ig3,ig4,idum12,idum13,idum14,idum15,idum16,idum17,      &
      &        idum18)


! Read LG field. Use GL field as backup.
! **CAUTION**: Discontinuities in GL field may cause interpolation problems! LG field is preferable.

  llg=.false.
  key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'LG')
  if ( key <  0 ) then
!    write(6,*) ' LAND_ICE_MASK_ATMS: The LG field is MISSING. Will try GL. '
    key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'GL')
    if ( key <  0 ) then
      write(6,*) ' LAND_ICE_MASK_ATMS: No ice (LG or GL) fields found. Aborting! '
      call abort()
    else
!      write(6,*) ' LAND_ICE_MASK_ATMS: The GL field was found and will be used.'
    endif
  else
    llg=.true.
  endif

  allocate ( lg(nilg*njlg), stat=alloc_status(8) )
  if( any(alloc_status /= 0) ) then
    write(6,*) ' LAND_ICE_MASK_ATMS: Memory allocation error '
    call abort()
  endif
  
  if ( llg ) then
    ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','LG')
  else
    ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','GL')
  endif

  ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,          &
      &        idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyplg,ig1lg,ig2lg,  &
      &        ig3lg,ig4lg,idum12,idum13,idum14,idum15,idum16,idum17,        &
      &        idum18)

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

  ier    = ezsetopt('INTERP_DEGREE','LINEAR')
  gdid   = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
  gdidlg = ezqkdef(nilg,njlg,grtyplg,ig1lg,ig2lg,ig3lg,ig4lg,iungeo)

  allocate ( mgintrp(npts), stat=alloc_status(9) )
  allocate ( lgintrp(npts), stat=alloc_status(10) )
  if( any(alloc_status /= 0) ) then
    write(6,*) ' LAND_ICE_MASK_ATMS: Memory allocation error '
    call abort()
  endif

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
  ztt(:) = -1             ! no ice (reset terain type)
  zlq(:) = 0              ! land   (reset land/sea qualifier)

  do kk = 1, npts
    if ( mgintrp(kk) < MGthresh ) zlq(kk) = 1  ! ocean point away from coast
    if ( lgintrp(kk) >= LGthresh .and. zlq(kk) == 1 ) ztt(kk) = 0  ! sea-ice affected point
    if ( lgintrp(kk)  < LGthresh .and. zlq(kk) == 1 ) then
      waterobs(kk) = .true.  ! water point not in close proximity to land or sea-ice
    end if
  end do

! Deallocate arrays and close FST file.

  alloc_status(:) = 0
  deallocate ( mgintrp, stat=alloc_status(1) )
  deallocate ( lgintrp, stat=alloc_status(2) )
  deallocate ( mg,      stat=alloc_status(3) )
  deallocate ( lg,      stat=alloc_status(4) )
  deallocate ( latmesh, stat=alloc_status(5) )
  deallocate ( lonmesh, stat=alloc_status(6) )
  deallocate ( mgintob, stat=alloc_status(7) )
  deallocate ( lgintob, stat=alloc_status(8) )
  deallocate ( zlatbox, stat=alloc_status(9) )
  deallocate ( zlonbox, stat=alloc_status(10) )
  if( any(alloc_status /= 0) ) then
    write(6,*) ' LAND_ICE_MASK_ATMS: Memory deallocation error '
    call abort()
  endif
  ier = fstfrm(iungeo)
  ier = fclos(iungeo)

return
end subroutine land_ice_mask_atms

!---------------------------------------------------------------------------------------

subroutine qc_data(zenith,ilq,itt,zlat,zlon,ztb,scanpos,stnid,nval,nt,lqc, &
     &            grossrej,lsq,trn,qcflag1,qcflag2,ican,blat,blon,lutb)

!
! SUBROUTINE QC_DATA
!
!  This routine performs basic quality control checks on the data. It sets array
!  lqc(nt,nchan) elements to .true. to flag data with failed checks.
!
!  The 7 QC checks are:
!                 1) Invalid land/sea qualifier or terrain type,
!                 2) Invalid field of view number,
!                 3) Satellite zenith angle missing or out of range, (> 75 deg),
!                 4) lat,lon check (lat,lon = O(-90.), 0(-180.))
!                 5) Change in (computed) lsq,trn from (input) ilq,itt (from MG,LG fields)
!                      ilq= 0,1 (from hi-res land/sea mask interpolated to obs point [CMDA])
!                      itt=-1,0 (from hi-res ice analysis  interpolated to obs point [CMDA])
!                      lsq= 0,1 (from max interp MG (0.0 to 1.0) in box surrounding obs point)
!                      trn=-1,0 (from max interp LG (0.0 to 1.0) in box surrounding obs point)
!                 6) ATMS quality flag check (qual. flag elements 33078,33079,33080,33081)

!
!  In most cases, lqc(ii,nchan) is set to .true. for all channels at point ii
!  if the check detects a problem. In addition, Tb (ztb) is set to missing_value 
!  for checks 3 and 4 fails.
!

!  use var_declare      ! nchan=22, zmisg, mxscan
  implicit none

  integer, intent(in), dimension(:)    :: ilq, itt, scanpos
  integer, intent(in), dimension(:)    :: ican, qcflag2
  integer, intent(in), dimension(:,:)  :: qcflag1
  integer, intent(in)                  :: nval, nt
  integer, intent(in), dimension(:)    :: lsq, trn

  logical, intent(in), dimension(:)    :: grossrej     ! dim(nt), true if 1 or more Tb fail gross error check
  logical, intent(out)                 :: lutb         ! true if Tb(ztb) are set to missing_value

  real, intent(in), dimension(:)       :: zlat, zlon
  real, intent(inout), dimension(:)    :: ztb, zenith
  
  integer, intent(in)                  :: blat, blon   ! NT box lat,lon (header)
   
  logical, intent(inout), dimension(:,:) :: lqc        ! dim(nt,nchan), lqc = .false. on input
   
  character(len=9), intent(in)         :: stnid
   
!  Locals
  integer :: ii, jj, indx1, icount
  logical :: fail, fail1, fail2


   write(6,*) '============================================================================'
   write(6,*) 'QC_DATA: Processing data box: Stnid, lat, lon = ', stnid, blat, blon
   write(6,*) ' '

   lutb = .false.

!--------------------------------------------------
!  Global rejection checks
!--------------------------------------------------

! Check if number of channels is correct

   if ( nval /= nchan ) then
      write(6,*) 'WARNING: Number of channels (',nval, ') is not equal to nchan (', nchan,')'
      write(6,*) '         All data flagged as bad and returning to calling routine!'
      lqc(:,:) = .true.  ! flag all data in report as bad
      return
   endif

! Check for errors in channel numbers (should be 1-22 for each location ii)
   indx1 = 1
   fail = .false.
   do ii = 1,nt
     do jj = 1,nchan
       if ( ican(indx1+jj-1) /= jj ) fail = .true.
     enddo
     indx1 = indx1 + nchan
   enddo
   if ( fail ) then
      write(6,*) 'WARNING: Bad channel number(s) detected!'
      write(6,*) '         All data flagged as bad and returning to calling routine!'
      write(6,*) '  ican(nt*nchan) array = ', ican(:)
      lqc(:,:) = .true.  ! flag all data in report as bad
      return
   endif

!---------------------------------------------------------
! 1) invalid land/sea qualifier or terrain type
!---------------------------------------------------------

!  ilq = 0 (land),     1 (sea)
!  itt = 0 (sea-ice), -1 otherwise
!  lsq = 1 (sea, away from land/coast [MG]),      0 otherwise
!  trn = 0 (over or near analyzed sea-ice [LG]), -1 otherwise


   do ii = 1,nt
     fail = .false.
     if ( ilq(ii) < 0  .or. ilq(ii) > 2 ) fail = .true.
     if ( itt(ii) < -1 .or. itt(ii) > 1 ) fail = .true.
     if ( fail ) then
       write(6,*) 'WARNING: Invalid land/sea qualifier or terrain type!'
       write(6,*) '  ilq, itt, (lat, lon) = ', ilq(ii), itt(ii), '(',zlat(ii), zlon(ii),')'
     endif
     if ( ilq(ii) == 0 .and. itt(ii) == 0 ) then
        fail = .true.
        write(6,*) 'WARNING: Sea ice point (itt=0) at land point (ilq=0)!'
        write(6,*) ' lat, lon =  ', zlat(ii), zlon(ii)
     endif
     if ( fail ) lqc(ii,:) = .true.
   enddo

   do ii = 1,nt
     fail = .false.
     if ( lsq(ii) < 0  .or. lsq(ii) > 2 ) fail = .true.
     if ( trn(ii) < -1 .or. trn(ii) > 1 ) fail = .true.
     if ( fail ) then
       write(6,*) 'WARNING: Invalid model-based (MG/LG) land/sea qualifier or terrain type!'
       write(6,*) '  lsq, trn, (lat, lon) = ', lsq(ii), trn(ii), '(',zlat(ii), zlon(ii),')'
     endif
     if ( fail ) lqc(ii,:) = .true.
   enddo
   
!---------------------------------------------------------
!  2) invalid field of view number
!---------------------------------------------------------
   
   do ii = 1,nt
     fail = .false.
     if ( scanpos(ii) < 1  .or. scanpos(ii) > mxscan ) then
        fail = .true.
        write(6,*) 'WARNING: Invalid field of view! scanpos, lat, lon = ', scanpos(ii), zlat(ii), zlon(ii)
     endif
     if ( fail ) lqc(ii,:) = .true.
   enddo

!--------------------------------------------------------------
!  3) satellite zenith angle missing or out of range (> 75 deg)
!--------------------------------------------------------------
!  If bad zenith, then set Tb (and zenith) = missing value
   
   indx1 = 1
   do ii = 1,nt
     fail = .false.
     if ( zenith(ii) > 75.0 .or. zenith(ii) < 0. ) then
       fail = .true.
       write(6,*) 'WARNING: Bad or missing zenith angle! zenith, lat, lon = ', zenith(ii), zlat(ii), zlon(ii)
       zenith(ii) = zmisg
       lutb = .true.
     endif
     do jj = 1,nchan
       if ( fail ) then
         lqc(ii,jj) = .true.
         ztb(indx1+jj-1) = zmisg
       endif
     enddo
     indx1 = indx1 + nchan
   enddo

!---------------------------------------------------------
! 4) Lat,lon check
!---------------------------------------------------------

! Check for undecoded BURP file integer values of lat,lon = 0,0
! (usually associated with missing zenith angle and erroneous Tb=330K)

   icount = 0
   indx1 = 1
   do ii = 1,nt
     fail = .false.
     if ( zlat(ii) == -90.0  .and. zlon(ii) == -180.0 ) then
       fail = .true.
       icount =  icount + 1
       lutb = .true.
     endif
     do jj = 1,nchan
       if ( fail ) then
         lqc(ii,jj) = .true.
         ztb(indx1+jj-1) = zmisg
       endif
     enddo
     indx1 = indx1 + nchan
   enddo
   if ( icount > 0 ) write(6,*) 'WARNING: Bad lat,lon pair(s) detected. Number of locations = ', icount

   icount = 0
   indx1 = 1
   do ii = 1,nt
     fail = .false.
     if ( abs(zlat(ii)) > 90.0  .or. abs(zlon(ii)) > 180.0 ) then
       fail = .true.
       icount =  icount + 1
       lutb = .true.
     endif
     do jj = 1,nchan
       if ( fail ) then
         lqc(ii,jj) = .true.
         ztb(indx1+jj-1) = zmisg
       endif
     enddo
     indx1 = indx1 + nchan
   enddo
   if ( icount > 0 ) write(6,*) 'WARNING: Lat or lon out of range! Number of locations = ', icount


!------------------------------------------------------------------------
!  5) Change in land/sea qualifier or terrain-type based on MG,LG fields
!------------------------------------------------------------------------

   icount = 0
   do ii = 1,nt
     fail = .false.
     if ( (ilq(ii) /= lsq(ii)) .or. (itt(ii) /= trn(ii)) ) fail = .true.
     if ( fail ) then
       icount =  icount + 1
     endif
   enddo
   if ( icount > 0 ) write(6,*) 'INFO: Num. pts with land/sea qualifier or terrain type changed (MG,LG) = ', icount

!----------------------------------------------------------------------------
!  6) ATMS quality flag check (qual. flag elements 33078,33079,33080,33081)
!----------------------------------------------------------------------------

!  33078 Geolocation quality code     qcflag1(ii,1)  code value = 0-15 (0= OK, 15=misg)
!  33079 Granule level quality flags  qcflag1(ii,2)  16 bit flag  (start bit 6(2^5)=32) (misg=2^16-1 = 65535)
!  33080 Scan level quality flags     qcflag1(ii,3)  20 bit flag  (start bit 7(2^6)=64) (misg=2^20-1) 
!  33081 Channel data quality flags   qcflag2        12 bit flag  (start bit 3(2^2)=4)  (misg=2^12-1)
!
!  See http://www.wmo.int/pages/prog/www/WMOCodes/WMO306_vI2/2010edition/BUFRver16/BUFR_16_0_0_TableD.pdf
!

   indx1 = 1
   do ii = 1,nt 
     fail1 = .false.
     fail = .false.
     if ( (qcflag1(ii,1) > 0) .or. (qcflag1(ii,2) >= 32) .or. (qcflag1(ii,3) >= 64) ) then
        write(6,*) 'WARNING: INFO BLOCK QC flag(s) indicate problem with data'
        write(6,*) ' ele33078 = ',qcflag1(ii,1),' ele33079 = ',qcflag1(ii,2),' ele33080 = ', qcflag1(ii,3)
        write(6,*) ' lat, lon = ', zlat(ii), zlon(ii)
        fail1 = .true.
        if ( grossrej(ii) ) write(6,*) ' NOTE: grossrej is also true for this point!'
     endif
     do jj = 1,nchan
       fail2 = .false.
       if ( qcflag2(indx1+jj-1) >= 4 ) then
!         write(6,*) 'WARNING: DATA BLOCK QC flag ele33081 = ', qcflag2(indx1+jj-1)
!         write(6,*) '    Lat, lon, channel = ', zlat(ii), zlon(ii), ican(indx1+jj-1)
         fail2 = .true.
         fail = .true.
!         if ( (.not. fail1) .and. grossrej(ii) ) write(6,*) ' NOTE: grossrej is also true for this point!'
       endif
       if ( fail2 .or. fail1 ) lqc(ii,jj) = .true.
     enddo
     if ( fail ) write(6,*) 'WARNING: DATA BLOCK QC flag ele33081 >= 4 for one or more channels! lat, lon = ', zlat(ii), zlon(ii)
     indx1 = indx1 + nchan
   enddo
   
   write(6,*) 'QC_DATA: Total number of data processed in this box = ', nt*nchan
   write(6,*) '         Total number of data flagged in this box   = ', COUNT(lqc)
   write(6,*) ' '

return
end subroutine qc_data

!---------------------------------------------------------------------------------------

subroutine nrl_filter(ier, ni, tb23, tb31, tb50, tb89, tb165, pangl, plat, ilansea, iglace, &
     &                waterobs, grossrej, clw, si_ecmwf, si_bg, iNumSeaIce, iRej, SeaIce)

!**SUBROUTINE  NRL_FILTER 
! 
!AUTEUR         S. Macpherson,  ARMA, dec 2012
!
!REVISIONS
!    REVISION 001   
!  
!LANGAGE        FORTRAN 5 
! 
!OBJET          Compute the following parameters using 5 ATMS channels:
!                  - sea ice, 
!                  - cloud liquid water (clw), 
!                  - 2 scattering indices (si) (ECMWF, Bennartz-Grody)
!               The five channels used are: 23Ghz, 31Ghz, 50Ghz, 89Ghz, and 165Ghz.
!
!NOTES*
!                o  open water points are converted to sea-ice points if sea ice concentration >= 0.55
!                   and iglace (itt or terrain type) is changed accordingly
!                o  clw are missing when out-of-range parameters/Tb detected or grossrej = .true.
!                o  clw and si only computed over open water away from coasts and sea-ice
!                o  clw and si = -99.0 where value cannot be computed.
!
!REFERENCES     Ben Ruston, NRL Monterey
!                  JCSDA Seminar 12/12/12: Impact of NPP Satellite Assimilation in the U.S. Navy Global Modeling System
!
!
!ARGUMENTS      ier         - output - error return code for each location:
!                                        0, ok,  
!                                        1, input parameter out of range or grossrej=.true. 
!               ni          - input  -  number of points to process (= NT)
!               tb23        - input  -  23Ghz brightness temperature (K) -- ch. 1
!               tb31        - input  -  31Ghz brightness temperature (K) -- ch. 2
!               tb50        - input  -  50Ghz brightness temperature (K) -- ch. 3
!               tb89        - input  -  89Ghz brightness temperature (K) -- ch. 16
!               tb165       - input  -  165Ghz brightness temperature (K) -- ch. 17
!               pangl       - input  -  satellite zenith angle (deg.)
!               plat        - input  -  latitude (deg.)
!               ilansea     - input  -  land/sea indicator (0=land, 1=ocean)
!               iglace      - in/out -  terrain type (0=ice, -1 otherwise)
!               waterobs    - in/out -  .true. if open water point (away from coasts and sea-ice)
!               grossrej    - input  -  .true. if any channel had a gross error from gross_value_check
!               clw         - output -  cloud liquid water (kg/m**2) from tb23 & tb31
!               si_ecmwf    - output -  ECMWF scattering index from tb89 & tb165
!               si_bg       - output -  Bennartz-Grody scattering index from tb89 & tb165
!               iNumSeaIce  - in/out -  running counter for number of open water points
!                                       with sea-ice detected (from algorithm)
!               iRej        - in/out -  running counter for number of locations with bad
!                                       pangl, plat, ilansea, or with grossrej=true
!               SeaIce      - output -  computed sea-ice fraction from tb23 & tb50 
!
!               ice         - internal -  sea ice
!             
!
!! Notes: In the case where an output parameter cannot be calculated, the
!!        value of this parameter is set to the missing value, i.e. -99.
!!
    implicit none
!
    integer    ::  i
!
    integer, intent(in)                   ::  ni
    integer, intent(inout)                ::  iNumSeaIce
    integer, intent(in), dimension(:)     ::  ilansea
    integer, intent(out), dimension(:)    ::  ier
    integer, intent(inout), dimension(:)  ::  iglace
    integer, intent(inout)                ::  iRej
    
    
    logical, intent(in), dimension(:)     ::  grossrej
    logical, intent(inout), dimension(:)  ::  waterobs
!
    real, intent(in),  dimension(:)  ::  tb23, tb31, tb50, tb89, tb165, pangl, plat
    real, intent(out), dimension(:)  ::  clw, si_ecmwf, si_bg, SeaIce 

    real, dimension(ni)              ::  ice
!
    real       ::  aa, deltb, abslat, cosz
    real       ::  t23, t31, t50, t89, t165

    logical    ::  debug

    real, parameter  ::  rmisg = -99.


    debug = .false.
    ier = 0
!
!!____1) Initialise parameters:
!!    -------------------------------
!
      do i = 1, ni
         ice(i)      = rmisg
         clw(i)      = rmisg
         si_ecmwf(i) = rmisg
         si_bg(i)    = rmisg
         SeaIce(i)   = 0.0
      enddo
!
!!____2) Validate input parameters:
!!    -----------------------------
!
      do i = 1, ni

         if ( pangl(i)   .lt.   0.  .or. &
     &        pangl(i)   .gt.  70.  .or. &
     &        plat(i)    .lt. -90.  .or. & 
     &        plat(i)    .gt.  90.  .or. &  
     &        ilansea(i) .lt.   0   .or. & 
     &        ilansea(i) .gt.   1        ) then
            ier(i) = 1
!            write(6,*) 'NRL_FILTER: Bad pangl(zenith), plat, ilansea =', pangl(i), plat(i), ilansea(i)
         endif
    ! Skip computations for points where all data are rejected  (bad Tb ANY channel)       
         if ( grossrej(i) ) ier(i) = 1 

      enddo
!
!!____3) Compute parameters:
!!    ----------------------
!
      do i = 1, ni
!
      if ( ier(i) .eq. 0 ) then
!
         abslat = abs(plat(i))
         cosz   = cosd(pangl(i))
         t23 = tb23(i)
         t31 = tb31(i)
         t50 = tb50(i)
         t89 = tb89(i)
         t165 = tb165(i)
         deltb = t89 - t165
!
!    Check for sea-ice over water points. Set terrain type to 0 if ice>=0.55 detected.
!
         if ( ilansea(i) .eq. 1 ) then  ! water point
!
            if ( abslat .lt. 50. ) then
               ice(i) = 0.0
            else
               ice(i) = 2.85 + 0.020*t23 - 0.028*t50
            endif
            
            SeaIce(i) = ice(i)
            
            if ( ice(i) .ge. 0.55 .and. waterobs(i) ) then
              iNumSeaIce = iNumSeaIce + 1
              waterobs(i) = .false.
              iglace(i) = 0
            endif
            
         endif
!
!    Compute CLW and Scattering Indices (over open water only)
!
       if ( waterobs(i) ) then
          if ( t23 .lt. 284. .and. t31 .lt. 284. ) then
            aa = 8.24 - (2.622 - 1.846*cosz)*cosz
            clw(i) = aa + 0.754*alog(285.0-t23) - 2.265*alog(285.0-t31)
            clw(i) = clw(i)*cosz
            if ( clw(i) .lt. 0.0 ) clw(i) = 0.0
          endif
          si_ecmwf(i) = deltb - (-46.94 + 0.248*pangl(i))
          si_bg(i)    = deltb - (-39.201 + 0.1104*pangl(i))
       endif
!

      else  ! ier(i) .eq. 1 case
      
         iRej = iRej + 1

      endif ! if ( ier(i) .eq. 0 )
!
      if ( debug .and. (i .le. 100) ) then
        write(6,*) ' '
        write(6,*) ' i,tb23(i),tb31(i),tb50(i),tb89(i),tb165(i),pangl(i),plat(i), ilansea(i) = ', &
     &             i,tb23(i),tb31(i),tb50(i),tb89(i),tb165(i),pangl(i),plat(i), ilansea(i)
        write(6,*) ' ier(i),ice(i),clw(i),si_ecmwf(i),si_bg(i),iglace(i),waterobs(i) =',ier(i),ice(i),&
     &             clw(i),si_ecmwf(i),si_bg(i),iglace(i),waterobs(i)
      endif
!
      enddo   ! i loop over ni points
!
      return


end subroutine nrl_filter

!---------------------------------------------------------------------------------------

subroutine gross_value_check(npts,ztb,grossrej)
!--------------------------------------------------------------------
!  Language: FORTRAN 90
!
!  Object: Check Tbs for values that are missing or outside physical limits.
!          **NOTE: REJECT ALL CHANNELS OF ONE IS FOUND TO BE BAD.
!
!  Call:   call gross_value_check(nt,ztb,grossrej)
!
!  History
!  Version:      Date:      Comment:
!  --------      -----      --------
!    0.1       28/03/07     Original code.                  D. Anselmo
!    0.2       07/02/13     Modified to catch NaN and Inf   S. Macpherson
!--------------------------------------------------------------------
! Variable Definitions:
! ---------------------
! npts            - input  -  number of obs pts to process
! ztb             - input  -  Tbs from input BURP file
! grossrej        - output -  logical array defining which obs are to be rejected
!------------------------------------------------------------------
!  use var_declare, only : nchan,zmisg
  implicit none

!  Arguments
  integer, intent(in) :: npts

  real,    intent(in),  dimension(:) :: ztb
  logical, intent(out), dimension(:) :: grossrej

!  Locals
  integer :: ii, indx1, indx2
!-----------------------------------------

  grossrej(1:npts) = .true.
  indx1 = 1
  do ii = 1, npts

    indx2 = ii*nchan
    if ( all( ztb(indx1:indx2) > 50.0 ) .and. all( ztb(indx1:indx2) < 380.0 ) ) then
      grossrej(ii) = .false.
    end if
    indx1 = indx2 + 1

  end do


return
end subroutine gross_value_check

!---------------------------------------------------------------------------------------

subroutine read_coeff(sats,chans,fovbias,coeff,nsat,nchan,nfov,npred,cinstrum,maxpred,iun,coeff_file,ptypes)

! max # of satellites, # of channels (24), # of FOV (60)
!use var_declare,  only : mxsat, mxan, mxscan, fnom, fclos

implicit none

! IN
integer                        :: maxpred, iun
character(len=90)              :: coeff_file

! OUT
character(len=9), dimension(mxsat) :: sats        ! dim(maxsat), satellite names
integer*4, dimension(mxsat,mxan)   :: chans       ! dim(maxsat, maxchan), channel numbers
real, dimension(mxsat,mxan,mxscan) :: fovbias     ! dim(maxsat,maxchan,maxfov), bias as F(fov)
real, dimension(mxsat,mxan,maxpred+1) :: coeff    ! dim(maxsat,maxchan,maxpred+1)
integer                            :: nsat, nfov, nbscan
integer, dimension(mxsat)          :: nchan       ! dim(maxsat), number of channels
integer, dimension(mxsat,mxan)     :: npred       ! dim(maxsat, maxchan), number of predictors
character(len=5)                   :: cinstrum    ! string: instrument (e.g. SSMIS)
character(len=2), dimension(mxsat,mxan,maxpred)  :: ptypes ! dim(maxsat,maxchan,maxpred)

! LOCAL
character(len=8)               :: sat
character(len=120)             :: line
integer*4                      :: chan
integer                        :: ndata, nbfov, nbpred, i, j, k, ier, istat, ii
logical                        :: newsat
real                           :: dummy


! RETURNS:

!   sats(nsat)            = satellite names
!   chans(nsat,nchan(i))  = channel numbers of each channel of each satellite i
!   npred(nsat,nchan(i))  = number of predictors for each channel of each satellite i
!   fovbias(i,j,k)        = bias for satellite i, channel j, FOV k   k=1,nfov
!     if FOV not considered for instrument, nfov = 1 and fovbias is global bias for channel
!   coeff(i,j,1)          = regression constant
!   coeff(i,j,2), ..., coeff(i,j,npred(i,j)) = predictor coefficients

!   nsat, nchan, nfov, cinstrum (output) are determined from file
!   if returned nsat = 0, coeff_file was empty

!   maxpred (input) is max number of predictors


 coeff    = 0.0
 fovbias  = 0.0
 sats     = 'XXXXXXXXX'
 cinstrum = 'XXXXX'
 chans    = 0
 npred    = 0
 nsat     = 0
 nchan    = 0
 nfov     = 0
 ptypes   = 'XX'

 nbscan = mxscan

 ier = FNOM(iun,coeff_file,'FMT',0)

IF (ier == 0) THEN

  WRITE(*,*)
  WRITE(*,*) 'Bias correction coefficient file open = ', coeff_file

  READ(iun,*,IOSTAT=istat)
  IF ( istat < 0 ) THEN
    WRITE(*,*) '  ERROR- File appears empty.'
    RETURN
  END IF
  REWIND(iun)

  ii = 0

! Loop over the satellites/channels in the file

  do
    read(iun,'(A)',IOSTAT=istat) line
    if ( istat < 0 ) EXIT
    if ( line(1:3) == 'SAT' ) then
        newsat = .true.
        read(line,'(T53,A8,1X,A5,1X,I6,1X,I8,1X,I2,1X,I3)',IOSTAT=istat) sat, cinstrum, chan, ndata, nbpred, nbfov
        if ( istat /= 0 ) then
          write(*,*) ' ERROR - reading data from SATELLITE line in coeff file!'
          return
        endif
        do i = 1, mxsat
          if ( trim(sats(i)) == trim(sat) ) then
            newsat = .false.
            ii = i
          endif
        end do
        if ( newsat ) then
          ii = ii + 1
          if ( ii > mxsat ) then
             write(*,*) ' ERROR - max number of satellites exceeded in coeff file!'
             return
          endif
          sats(ii) = sat
          if (ii > 1) nchan(ii-1) = j
          j = 1
        else
          j = j + 1
        endif
        chans(ii, j) = chan
        npred(ii, j) = nbpred
        if ( nbpred > maxpred ) then
           write(*,*) ' ERROR - max number of predictors exceeded in coeff file!'
           return
        endif
        read(iun,'(A)',IOSTAT=istat) line
        if ( line(1:3) /= 'PTY' ) then
           write(*,*) ' ERROR - list of predictors is missing in coeff file!'
           return
        endif
        if ( nbpred > 0 ) then
          read(line,'(T8,6(1X,A2))',IOSTAT=istat) (ptypes(ii,j,k),k=1,nbpred)
          if ( istat /= 0 ) then
            write(*,*) ' ERROR - reading predictor types from PTYPES line in coeff file!'
            return
          endif
        endif
        read(iun,*,IOSTAT=istat) (fovbias(ii,j,k),k=1,nbfov)
        if ( istat /= 0 ) then
          write(*,*) ' ERROR - reading fovbias in coeff file!'
          return
        endif
        if ( nbpred > 0 ) then
          read(iun,*,IOSTAT=istat) (coeff(ii,j,k),k=1,nbpred+1)
        else
          read(iun,*,IOSTAT=istat) dummy
        endif
        if ( istat /= 0 ) then
          write(*,*) ' ERROR - reading coeff in coeff file!'
          return
        endif

    else
        EXIT
    endif

  end do

  if ( ii == 0 ) then
     write(*,*) ' ERROR - No data read from coeff file!'
     return
  endif

  nsat      = ii
  nfov      = nbfov
  nchan(ii) = j
  if ( nbscan /= 0 ) then
    if ( nfov /= mxscan ) then
      write(*,*) ' INFO - Number of FOV in coeff file (nfov) does not equal default value (mxscan).'
      write(*,*) '         nfov = ', nfov
      write(*,*) '       mxscan = ', mxscan
!      call abort()
    endif
  else ! nbscan = 0 case
    if ( nfov /= 1 ) then
      write(*,*) ' INFO - Number of FOV in coeff file (nfov) does not equal default value (1).'
      write(*,*) '         nfov = ', nfov
!      call abort()
    endif
  endif

  write(*,*) ' '
  write(*,*) ' ------------- BIAS CORRECTION COEFFICIENT FILE ------------------ '
  write(*,*) ' '
  write(*,*) ' Number of satellites =     ', nsat
  write(*,*) ' Number of FOV =            ', nfov
  write(*,*) ' Max number of predictors = ', maxval(npred)
  write(*,*) ' '
  do i = 1, nsat
    write(*,*) '  Satellite = ' // sats(i)
    write(*,*) '     Number of channels = ', nchan(i)
    write(*,*) '     predictors, fovbias, coeff for each channel: '
    do j = 1, nchan(i)
      write(*,*) i, chans(i,j)
      if ( npred(i,j) > 0 ) then 
        write(*,'(6(1X,A2))') (ptypes(i,j,k),k=1,npred(i,j))
      else
        write(*,'(A)') 'No predictors'
      endif
      write(*,*) (fovbias(i,j,k),k=1,nfov)
      write(*,*) (coeff(i,j,k),k=1,npred(i,j)+1)
    end do
  end do
  write(*,*) ' '

  ier = FCLOS(iun)

ELSE

  write(*,*) 'READ_COEFF: ERROR - Problem opening the coeff file!'

ENDIF

RETURN

end subroutine read_coeff

!------------------------------------- HANDLE_ERROR -----

subroutine handle_error()
implicit none
    write(*,*) BURP_STR_ERROR()
    write(*,*) "history"
    Call BURP_STR_ERROR_HISTORY()
    Deallocate(adresses)
    Call BURP_Free(File_in,F2=File_out)
    Call BURP_Free(Rpt_in,R2=Rpt_out)
    Call BURP_Free(Block_in,B2=Block_copy)
    call abort()
end subroutine handle_error

end program midas_satQCATMS

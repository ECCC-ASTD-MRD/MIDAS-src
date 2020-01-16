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

program midas_bgckmw
  !
  ! :Purpose: Main program for background check of microwave instruments. 
  !
  use burp_module
  use bgckmicrowave_mod

  implicit none

  integer :: error, nb_rpts, ref_rpt, compteur, nbele, nlocs
  integer :: nobs_tot, n_bad_reps, n_reps, reportIndex

  type(BURP_FILE) :: File_in, File_out
  type(BURP_RPT)  :: Rpt_in, Rpt_out

  integer,allocatable :: adresses(:)
  integer             :: ibrptime, nblocs, nsize, iun_burpin, irun

  logical :: bad_report, resume_report, qcflag1_check

  INTEGER MXELM
  INTEGER MXLAT, MXLON
  PARAMETER ( MXELM  =    40 )
  PARAMETER ( MXLAT  =     5 )
  PARAMETER ( MXLON  =     5 )
  integer, parameter :: MXSAT = 9
  integer, PARAMETER :: MXVAL = 22
  integer, PARAMETER :: MXNT = 3000
  integer, parameter :: nchanAtms=22
  integer, parameter :: mxscan=96
  real, parameter    :: zmisg=9.9e09
  integer, parameter :: MISGINT = -1

  integer ezsetopt, ezsetval, ezqkdef
  integer gdllsval, gdxyfll, gdmg, gdmt, gdgl

  INTEGER FSTOUV
  INTEGER FSTINF,FSTPRM,FSTLIR,FSTFRM
  INTEGER HANDLE,ISTAT,NOMBRE,NVAL,NT
  INTEGER EXDB,EXFIN 
  INTEGER NPOSIT, IER,IREC,IREC2,JUNK
  INTEGER I,ILNMX, NELE, J, JN, JL, blat, blon, idtyp
  INTEGER IUNGEO, IUNSTAT, INUMSAT, nulnam 
  INTEGER nvalOut, ntOut, INOSAT
  INTEGER IDUM,IDUM1,IDUM2,IDUM3,IDUM4,IDUM5,IDUM6,IDUM7
  INTEGER IDUM8,IDUM9,IDUM10,IDUM11,IDUM12,IDUM13
  INTEGER IDUM14,IDUM15,IDUM16,IDUM17,IDUM18
  INTEGER IG1GL,IG2GL,IG3GL,IG4GL
  INTEGER IG1MG,IG2MG,IG3MG,IG4MG
  INTEGER IG1MT,IG2MT,IG3MT,IG4MT
  INTEGER NITIC,NJTAC,NIGL,NJGL,NK,INDX,NLAT,NLON
  INTEGER NIMG,NJMG,NIMT,NJMT

  INTEGER ilq       (MXNT)
  INTEGER itt       (MXNT)
  INTEGER ISCNCNT   (MXNT)
  INTEGER scanpos   (MXNT)
  INTEGER ISAT      (MXNT)
  INTEGER IORBIT    (MXNT)
  INTEGER ican      (MXVAL*MXNT)
  INTEGER ICANOMP   (MXVAL*MXNT)
  integer qcflag2   (mxval*mxnt)
  integer qcflag1   (mxnt,3)
  INTEGER ICHECK    (MXVAL*MXNT)
  INTEGER ICHKPRF   (MXNT)
  INTEGER IMARQ     (MXVAL*MXNT)
  real    clw       (MXNT)
  real    clw_avg   (MXNT)
  real    scatw     (MXNT)

  CHARACTER(len=9)   STNID
  CHARACTER(len=9)   CSATID(MXSAT)

  character(len=90)  brp_in,brp_out
  character(len=20)  opt_missing
  CHARACTER(len=12)  ETIKXX
  CHARACTER(len=9)   ETIKRESU, id
  CHARACTER(len=4)   NOMVXX,CLNOMVAR
  CHARACTER(len=2)   TYPXX 
  CHARACTER(len=1)   GRTYPGL,GRTYPMG,GRTYPMT

  INTEGER, ALLOCATABLE, DIMENSION(:) :: BUF1
  REAL, ALLOCATABLE, DIMENSION(:) :: MG
  REAL, ALLOCATABLE, DIMENSION(:) :: MT
  REAL, ALLOCATABLE, DIMENSION(:) :: GL
  REAL, ALLOCATABLE, DIMENSION(:) :: TIC
  REAL, ALLOCATABLE, DIMENSION(:) :: TAC

  REAL  DONIALT (MXELM*MXVAL*MXNT)
  REAL  PRVALN  (MXELM*MXVAL*MXNT)
  REAL  ZDATA   (MXVAL*MXNT)
  REAL  ZLAT    (MXNT)
  REAL  ZLON    (MXNT)
  REAL  zenith  (MXNT)
  REAL  MGINTRP (MXNT)
  REAL  MTINTRP (MXNT)
  REAL  GLINTRP (MXNT)
  REAL  ztb     (MXVAL*MXNT)
  REAL  biasCorr(MXVAL*MXNT)
  REAL  ZOMP    (MXVAL*MXNT)
  REAL  ZLATBOX (MXLAT*MXLON,MXNT)
  REAL  ZLONBOX (MXLAT*MXLON,MXNT)
  REAL  MGINTBOX(MXLAT*MXLON,MXNT)
  REAL  MTINTBOX(MXLAT*MXLON,MXNT)
  REAL  GLINTBOX(MXLAT*MXLON,MXNT)
  REAL  XLAT,XLON

  REAL  DLAT, DLON, TOPOFACT

  LOGICAL RESETQC, SKIPENR

  DATA IUNGEO  / 50 /
  DATA IUNSTAT / 60 /
  DATA DLAT   / 0.4 /
  DATA DLON   / 0.6 /

  EXTERNAL EXDB,EXFIN

  LOGICAL DEBUG, clwQcThreshold

  namelist /nambgck/ debug, RESETQC, ETIKRESU, clwQcThreshold

  JUNK = EXDB('BGCKMW','DEBUT','NON')

  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM MIDAS-BGCKMW: --",/,' //   &
            '14x,"-- BACKGROUND CHECK FOR MW OBSERVATIONS --",/, ' //&
            '14x,"-- Revision : ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  ! 1) Debut
  IER = FNOM(IUNGEO,'./fstglmg' ,'STD+RND+R/O',0)

  ! default values
  debug = .false.
  RESETQC = .FALSE.
  ETIKRESU = '>>BGCKALT'
  clwQcThreshold = 0.3

  ! reading namelist
  nulnam = 0
  ier = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam, nml=nambgck, iostat=ier)
  if ( ier /= 0 ) then 
    write(*,*) 'midas_bgckmw: Error reading namelist'
    call abort()
  end if 
  write(*,nml=nambgck)
  ier = fclos(nulnam)

  mwbg_debug = debug
  mwbg_clwQcThreshold = clwQcThreshold

  brp_in = './obsto_amsua'
  brp_out = './obsto_amsua.out'

  ! initialisation
  Call BURP_Init(File_in,  F2=File_out,  IOSTAT=error)
  Call BURP_Init(Rpt_in,   R2=Rpt_out,   IOSTAT=error)

  ! Set BURP "missing value" for reals
  opt_missing = 'MISSING'
  Call BURP_Set_Options(REAL_OPTNAME=opt_missing,REAL_OPTNAME_VALUE=zmisg)
  IF ( DEBUG ) WRITE(*,*)' MISSING VALUE =', ZMISG

  ! ouverture du fichier burp d'entree et de sortie
  Call BURP_New(File_in,  FILENAME= brp_in,  MODE= FILE_ACC_READ,   IOSTAT= error)
  Call BURP_New(File_out, FILENAME= brp_out, MODE= FILE_ACC_CREATE, IOSTAT= error)

  ! Number of reports and maximum report size from input BURP file
  Call BURP_Get_Property(File_in, NRPTS=nb_rpts, IO_UNIT= iun_burpin)
  if ( nb_rpts <= 1 ) then
    write(*,*) 'The input BURP file ''', trim(brp_in), ''' is empty!'
    stop
  end if

  nsize = MRFMXL(iun_burpin)

  write(*,*)
  write(*,*) 'Number of reports containing observations = ', nb_rpts-1
  write(*,*) 'Size of largest report = ', nsize
  write(*,*)

  ! Add nsize to report size to accomodate modified (larger) data blocks
  nsize = nsize * 3

  allocate(adresses(nb_rpts), stat=error)
  
  if (error /= 0) then
    write(*,*) 'ERROR - allocate(adresses(nb_rpts)). alloc_status =' , error
    call abort()
  endif
  
  adresses(:) = 0

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

  ! 2) Lecture des statistiques d'erreur totale pour les  TOVS 
  IER = FNOM(IUNSTAT,'./stats_amsua_assim','SEQ+FMT',0)
  IF(IER.LT.0)THEN
    WRITE (*,*) '(" bgckMW: Problem opening ", &
           "TOVS total error statistics file ", stats_amsua_assim)'               
    CALL ABORT ()
  END IF
  CALL mwbg_readStatTovs(IUNSTAT,INUMSAT,CSATID)
  WRITE(*,*) " SATID's = "
  DO I = 1, INUMSAT
    WRITE(*,*) '  ', CSATID(I)
  ENDDO

  ! 3) Lecture des champs geophysiques (GL,MG,MT) du modele
  IER = FSTOUV(IUNGEO,'RND')
  ! MG
  IREC = FSTINF(IUNGEO,NIMG,NJMG,NK,-1,' ',-1,-1,-1,' ','MG')
  IF (IREC .LT. 0) THEN
    WRITE (*,*) ' LE MASQUE TERRE-MER EST INEXISTANT' 
    CALL ABORT()
  ENDIF

  ALLOCATE ( mg(NIMG*NJMG), STAT=ier)
  IER = FSTLIR(MG,IUNGEO,NIMG,NJMG,NK,-1,' ',-1,-1,-1,&
               ' ','MG')

  IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, &
       IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10, &
       IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYPMG, IG1MG,&
       IG2MG, IG3MG, IG4MG, IDUM12, IDUM13, IDUM14, &
       IDUM15, IDUM16, IDUM17, IDUM18 )

  ! TOPOGRAPHIE (ME ou MX).
  ! ME est la topographie avec unites en metres.
  ! MX est la topographie avec unites en m2/s2.
  TOPOFACT = 1.0
  IREC = FSTINF(IUNGEO,NIMT,NJMT,NK,-1,' ',-1,-1,-1,' ','ME')
  CLNOMVAR = 'ME'
  IF (IREC .LT. 0) THEN
    TOPOFACT = 9.80616
    IREC = FSTINF(IUNGEO,NIMT,NJMT,NK,-1,' ',-1,-1,-1,' ','MX')
    CLNOMVAR = 'MX'
  ENDIF
  IF (IREC .LT. 0) THEN
    WRITE (*,*) ' LA TOPOGRAPHIE EST INEXISTANTE' 
    CALL ABORT()
  ENDIF

  ALLOCATE ( MT(NIMT*NJMT), STAT=ier)
  IER = FSTLIR(MT,IUNGEO,NIMT,NJMT,NK,-1,' ',-1,-1,-1,&
               ' ',CLNOMVAR)

  IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
       IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10, & 
       IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYPMT, IG1MT, &
       IG2MT, IG3MT, IG4MT, IDUM12, IDUM13, IDUM14, &
       IDUM15, IDUM16, IDUM17, IDUM18 )

  ! GL
  IREC = FSTINF(IUNGEO,NIGL,NJGL,NK,-1,' ',-1,-1,-1,' ','GL')
  IF (IREC .LT. 0) THEN
    WRITE (*,*) 'LE CHAMP GL EST INEXISTANT' 
    CALL ABORT()
  ENDIF

  ALLOCATE ( gl(NIGL*NJGL), STAT=ier)
  IER = FSTLIR(GL,IUNGEO,NIGL,NJGL,NK,-1,' ',-1,-1,-1, &
               ' ','GL')

  IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
       IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10, &
       IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYPGL, IG1GL, &
       IG2GL, IG3GL, IG4GL, IDUM12, IDUM13, IDUM14, &
       IDUM15, IDUM16, IDUM17, IDUM18 )

  ! On regle les parametres pour l'interpolation avec EZSCINT
  ier = ezsetopt('INTERP_DEGREE','LINEAR')
  WRITE (*,*) ' GRILLE MG : ',grtypmg,nimg,njmg, &
              ig1mg,ig2mg,ig3mg,ig4mg
  gdmg= ezqkdef(nimg,njmg,grtypmg,ig1mg,ig2mg,ig3mg,ig4mg,iungeo)
  WRITE (*,*) ' GRILLE MT : ',grtypmt,nimt,njmt, &
              ig1mt,ig2mt,ig3mt,ig4mt
  gdmt= ezqkdef(nimt,njmt,grtypmt,ig1mt,ig2mt,ig3mt,ig4mt,iungeo)
  WRITE (*,*) ' GRILLE GL : ',grtypgl,nigl,njgl, &
              ig1gl,ig2gl,ig3gl,ig4gl
  gdgl= ezqkdef(nigl,njgl,grtypgl,ig1gl,ig2gl,ig3gl,ig4gl,iungeo)

  ! MAIN LOOP through all the reports in the BURP file
  nobs_tot = 0
  n_reps = 0
  n_bad_reps = 0
  
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

      !  Get all the required data from the blocks in the report (Rpt_in)
      call mwbg_getData(reportIndex, Rpt_in, ISAT, zenith, ilq, itt, zlat, zlon, ztb, &
                        biasCorr, ZOMP, scanpos, nvalOut, ntOut, qcflag1, qcflag2, &
                        ican, icanomp, IMARQ, IORBIT, bad_report, qcflag1_check = .FALSE.)
      if ( bad_report ) then
        n_bad_reps = n_bad_reps + 1  

        Call BURP_Free(Rpt_out,IOSTAT=error)
        if (error /= burp_noerr)  call handle_error()

        cycle REPORTS
      end if

      ! trouver l'indice du satellite
      INOSAT = 0
      DO I = 1,MXSAT
        IF ( STNID .EQ. '^'//CSATID(I) ) THEN
          INOSAT = I
        ENDIF
      ENDDO
      IF ( INOSAT .EQ. 0 ) THEN
        WRITE(*,*) 'SATELLITE ',TRIM(STNID), &
                   ' NOT FOUND IN STATS FILE!'
        CALL ABORT()
      ENDIF

      ! 5) Interpolation de la glace et le champ terre/mer du modele aux pts TOVS.
      ! N.B.: on examine ces champs sur une boite centree sur chaque obs.
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

      ier = gdllsval(gdmg,mgintbox,mg,ZLATBOX,ZLONBOX,MXLAT*MXLON*NT)
      IF (ier .lt. 0) THEN
        WRITE(*,*) 'ERROR in the interpolation of MG'
        CALL ABORT()
      ENDIF
      ier = gdllsval(gdmt,mtintbox,mt,ZLATBOX,ZLONBOX,MXLAT*MXLON*NT)
      IF (ier .lt. 0) THEN
        WRITE(*,*) 'ERROR in the interpolation of MT'
        CALL ABORT()
      ENDIF
      ier = gdllsval(gdgl,glintbox,gl,ZLATBOX,ZLONBOX,MXLAT*MXLON*NT)
      IF (ier .lt. 0) THEN
        WRITE(*,*) 'ERROR in the interpolation of GL'
        CALL ABORT()
      ENDIF

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
          PRINT *, ' MGINTBOX = '
          PRINT *,  (MGINTBOX(I,JN),I=1,MXLAT*MXLON)
          PRINT *, ' MTINTBOX = '
          PRINT *,  (MTINTBOX(I,JN),I=1,MXLAT*MXLON)
          PRINT *, ' GLINTBOX = '
          PRINT *,  (GLINTBOX(I,JN),I=1,MXLAT*MXLON)
        ENDIF
        MGINTRP(JN) = 0.0
        MTINTRP(JN) = 0.0
        GLINTRP(JN) = 0.0
        DO I=1,MXLAT*MXLON
          MGINTRP(JN) = MAX(MGINTRP(JN),MGINTBOX(I,JN))
          MTINTRP(JN) = MAX(MTINTRP(JN),MTINTBOX(I,JN)/TOPOFACT)
          GLINTRP(JN) = MAX(GLINTRP(JN),GLINTBOX(I,JN))
        ENDDO
        IF (DEBUG) THEN
          PRINT *, ' MGINTRP = ', MGINTRP(JN)
          PRINT *, ' MTINTRP = ', MTINTRP(JN)
          PRINT *, ' GLINTRP = ', GLINTRP(JN)
        ENDIF
      ENDDO

      ! 6) Controle de qualite des TOVS. Data QC flags (IMARQ) are modified here!
      CALL mwbg_tovCheckAmsua(ISAT, ilq, IORBIT, ican, ICANOMP, ztb, biasCorr, &
                              ZOMP, ICHECK, nvalOut, ntOut, ZMISG, INOSAT, ICHKPRF, &
                              scanpos, MGINTRP, MTINTRP, GLINTRP, itt, zenith, &
                              IMARQ, clw, clw_avg, scatw, STNID, RESETQC, ZLAT)

      ! Accumuler Les statistiques sur les rejets
      CALL mwbg_qcStatsAmsua(INUMSAT, ICHECK, ican, INOSAT, CSATID, nvalOut, &
                             ntOut, .FALSE.)

      ! 7) Mise a jour de rapport.
      CALL mwbg_updateBurpAmsua(clw_avg, scatw, ICHKPRF, ilq, itt, &
                                GLINTRP, ICHECK, RESETQC, IMARQ, &
                                Rpt_in, Rpt_out)

    end if

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

  write(*,*) ' Number of obs pts read from BURP file              = ', nobs_tot
  write(*,*) ' Number of BURP file reports                        = ', n_reps
  write(*,*) ' Number of bad BURP file reports (all data flagged) = ', n_bad_reps

  ! 10) Fin
  ! Imprimer les statistiques sur les rejets
  CALL mwbg_qcStatsAmsua(INUMSAT, ICHECK, ican, INOSAT, CSATID, nvalOut, &
                         ntOut, .TRUE.)

  deallocate(adresses)

  junk  = exfin('bgckMW','FIN','NON')

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

end program midas_bgckmw

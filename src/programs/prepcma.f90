program prepcma
  !
  ! author Peter Houtekamer, Herschel Mitchell and Gerard Pellerin
  !        December 1999 - February 2002 
  ! 
  ! Revision June 2005. Pass the input variables via a namelist (instead
  !         of using the ccard routine).
  !          March 2007. Added scatterometer data.
  !          January 2008. Added GPS/RO observations.
  ! Revision Xingxiu Deng December 2011
  !        Add to output qcvar (to be used in zhacma) in file BRPFORM
  ! Revision Xingxiu Deng January 2014
  !        Add to reject any observation outside the assimilation time window
  ! Revision Xingxiu Deng November 2014
  !        Add atms data (poscatms)
  ! Revision  Bin He      Dec. 2018  
  !        Using Madis  Lib.          
  !
  !object: read observation files that are in bgckalt format 
  !         (i.e. output from the background check) and
  !        transform them to the ObsSpaceData format. All observations will 
  !        collected in a single ObSSpaceData structure. The structure is 
  !        output in binary format. 
  !
  !input:  
  !        the following observation files may have as filename NONE in
  !        which case they will not be read.
  !        brpai: AIREP family,
  !        brpua: Upper Air family,
  !        brpsc: Scatterometer family,
  !        brpsf: Surface family,
  !        brpsw: Satob/satWind family,
  !        brptoa: TOvs family (amsu-a),
  !        brptob: TOvs family (amsu-b),
  !        brpairs: TOvs family (AIRS),
  !        brpiasi: TOvs family (IASI),
  !        brpcris: TOvs family (CrIS),
  !        brppr: profiler family,
  !        brpro: GPS/RO family
  !        brpatms: ATMS family ,
  !        timestamps: CMC stamp of the starting time of the assimilation window
  !        timestampe: CMC stamp of the ending time of the assimilation window
  !        infl_sigo: if 'T': inflate observation errors as in the envar analysis
  !                   if 'F': use regular observation errors as in the background 
  !                           check of the envar.
  !        thinning:  if 'T': perform data thinning.
  !                      'F': no additional data thinning (beyond what has 
  !                           already been done for the input files).
  !        dyn_sw_oer: if 'T': read satwind observation error from brpsw bgckalt file.
  !                       'F': use static satwind observation errors.
  !
  !output: L: standard output file,
  !        OBSOUT: Ascii file with the CMA-info,
  !        CMAHDR: unformatted ObsSpaceData header data,
  !        CMABDY: unformatted ObsSpaceData body data.
  !        CMADIM: Ascii file with the dimension of the ObsSpaceData structure 
  !            (number of stations, number of observations, number
  !             of ensemble members)
  !        BRPFORM: Ascii file with qcvar value to indicate postalt or bgckalt input.
  !
  ! See routine suprep for the list of elements that can be assimilated
  ! by the EnKF.
  !
  use bufr_mod 
  use ObsSpaceData_mod
  use burpread_mod 
  use time_mod
  use MathPhysConstants_mod
  implicit none

  type (struct_obs) :: cma                ! the CMA being prep'ed

  character(len=12) :: filestat
  character(len=11) :: fileform

  ! In 2011 the value below seemed to be generously high. If more observations 
  ! than this are encountered recompile code with a higher value. 
  integer, parameter :: MXOBSTOTAL=72000000, &
                        MXSTNTOTAL=5000000, &
                        nflags=4, &
                        nfilesmax=18

  ! number of pressure ranges used for the thinning of 
  ! aircraft data.
  integer, parameter :: npres_ai = 5
  integer, parameter :: npres_sw = 2
  integer, parameter :: nai_target = 10, &
                        nsc_target = 10, &
                        nsw_target = 6, &
                        nto_target = 6 
  real*8,  dimension(npres_ai) :: nai_pmax = &
         (/ 25000.0, 40000.0, 60000.0, 80000.0, 110000.0/)
  real*8,  dimension(npres_sw) :: nsw_pmax = &
         (/ 60000.0, 110000.0/)
! For a scalar array, no layer selection will be done
  real*8,  dimension(1) :: nsc_pmax = (/ 0.0 /), &
                          nto_pmax = (/ 0.0 /)
! assume that AMSU-B observations in the tropics have error AMSUB_trop_oer 
  real*4, parameter :: AMSUB_trop_oer = 1.0

  integer :: ier,ifiles,inamelist,myip, &
             nbrpform,ncmahdr,ncmahx,        &
             ncmabdy,ncmadim,nfiles,nmaxlen,nobsout, &
             nulout
  integer, dimension(:), allocatable :: nbegintyp,nendtyp           
  integer, dimension(nflags) :: nlistflg
  integer, parameter :: zero_ens=0
  character (len=2)   :: cfamtyp(nfilesmax)
  character (len=256) :: cfilnam(nfilesmax)
  character (len=7) :: type_resume 
     
  real(kind=8), dimension(1,1) :: hx
  real(kind=8) :: lat_obs
 
  integer :: timestamps,timestampe
  integer :: nspace,num_stn,istn,dateobs,timeobs,timestampobs
  integer :: idata,idataend,idum,jdata,idbrp
  logical :: infl_sigo_bin,memaster,qcvar,thinning_bin,dyn_sw_oer_bin

  character(len=1)   :: infl_sigo,thinning,dyn_sw_oer
  character(len=256) :: brpai,brpatms,brppr,brpro,brpsc,brpsf,brpsw,brptoa, &
    brptob,brpairs,brpiasi,brpcris,brpua,brpform,cmahdr,cmabdy,cmadim,l,obsout

  namelist /NAMPRE/ l,obsout,brpai,brpua,brpsf,brpsw,brpsc,brptoa,  &
    brptob,brpairs,brpiasi,brpcris,brppr,brpro,brpatms,timestamps, &
    timestampe,cmahdr,cmabdy,cmadim,brpform,infl_sigo,thinning, &
    dyn_sw_oer

  ! This process is the one and only process. It is therefore like the 
  ! master process in an mpi/environment. Diagnostic messages will be 
  ! printed.
  memaster=.true.
  myip=0

  ncmahx=-1 
  ! ASCII output file
  nulout=7
  ! ASCII file with the contents of the CMA file
  nobsout=90
  ! cma file (header info)
  ncmahdr=13
  ! cma file (body info)
  ncmabdy=ncmahdr+1
  ! ascii file with dimensions for the cma structure
  ncmadim=ncmabdy+1
  ! namelist
  inamelist=ncmadim+1
  ! ASCII  file with qcvar
  nbrpform=inamelist+1

  open(inamelist,FILE='flnml',DELIM='APOSTROPHE')
  read(inamelist,NAMPRE)
  close(inamelist)

  if (infl_sigo.eq.'T') then
    infl_sigo_bin=.true.
    write(*,*) 'Inflate radiance observation error as in envar.'
  else
    infl_sigo_bin=.false.
    write(*,*) 'Use regular radiance observation error ', &
       'as in the background check of the envar.'
  endif
  if (thinning.eq.'T') then
    thinning_bin=.true.
    write(*,*) 'Impose maximum observation densities ', &
       'for various families.'
  else
    thinning_bin=.false.
    write(*,*) 'Do not perform additional data thinning. '
  endif
  if (dyn_sw_oer.eq.'T') then
    dyn_sw_oer_bin=.true.
    write(*,*) 'Use dynamic satwind observation errors '
  else
    dyn_sw_oer_bin=.false.
    write(*,*) 'Use static satwind observation errors. '
  endif


  filestat='NEW'
  fileform='FORMATTED'
! open output file for debugging information
  call openfile(nulout,l,filestat,fileform)
! open output file for observation information
  call openfile(nobsout,obsout,filestat,fileform)

  call openfile(ncmadim,cmadim,filestat,fileform)

  fileform='UNFORMATTED'
  call openfile(ncmahdr,cmahdr,filestat,fileform)  
  call openfile(ncmabdy,cmabdy,filestat,fileform)

! determine the number of files with non-trivial input
  ifiles=0
  if (brpai.ne.'NONE') then
    ifiles=ifiles+1
    cfilnam(ifiles)=brpai
    cfamtyp(ifiles)='AI'
  endif
  if (brpua.ne.'NONE') then 
    ifiles=ifiles+1
    cfilnam(ifiles)=brpua
    cfamtyp(ifiles)='UA'
  endif
  if (brpsc.ne.'NONE') then 
    ifiles=ifiles+1
    cfilnam(ifiles)=brpsc
    cfamtyp(ifiles)='SC'
  endif
  if (brpsf.ne.'NONE') then 
    ifiles=ifiles+1
    cfilnam(ifiles)=brpsf
    cfamtyp(ifiles)='SF'
  endif
  if (brpsw.ne.'NONE') then 
    ifiles=ifiles+1
    cfilnam(ifiles)=brpsw
    cfamtyp(ifiles)='SW'
  endif
  if (brptoa.ne.'NONE') then 
    ifiles=ifiles+1
    cfilnam(ifiles)=brptoa
    cfamtyp(ifiles)='TO'
  endif
  if (brptob.ne.'NONE') then 
    ifiles=ifiles+1
    cfilnam(ifiles)=brptob
    cfamtyp(ifiles)='TO'
  endif
  if (brpairs.ne.'NONE') then 
    ifiles=ifiles+1
    cfilnam(ifiles)=brpairs
    cfamtyp(ifiles)='TO'
  endif
  if (brpiasi.ne.'NONE') then 
    ifiles=ifiles+1
    cfilnam(ifiles)=brpiasi
    cfamtyp(ifiles)='TO'
  endif
  if (brpcris.ne.'NONE') then 
    ifiles=ifiles+1
    cfilnam(ifiles)=brpcris
    cfamtyp(ifiles)='TO'
  endif
  if (brppr.ne.'NONE') then 
    ifiles=ifiles+1
    cfilnam(ifiles)=brppr
    cfamtyp(ifiles)='PR'
  endif
  if (brpro.ne.'NONE') then 
    ifiles=ifiles+1
    cfilnam(ifiles)=brpro
    cfamtyp(ifiles)='RO'
  endif
  if (brpatms.ne.'NONE') then 
    ifiles=ifiles+1
    cfilnam(ifiles)=brpatms
    cfamtyp(ifiles)='TO'
  endif

  nfiles=ifiles
  if (nfiles.eq.0) then 
    write(nulout,*) 'ERROR: NO OBSERVATION FILES PRESENTED!!'
    call qqexit(1)
    stop
  else 
    write(nulout,*) 'List of input observation files'
    do ifiles=1,nfiles
      write(nulout,*) 'family: ',cfamtyp(ifiles),'file: ',cfilnam(ifiles)
    enddo
    allocate(nbegintyp(nfiles),stat=ier)
    allocate(nendtyp(nfiles),stat=ier)
  endif
      
  !    Discussion of rejection criteria for the EnKF as of February 2002.     
  !    The 3d-var, via namelist entry, will only reject on criterion 4 and 5. 
  !    Routine suobs.ftn (from the operational 3d-var) lists (but does 
  !    not use (since the list is overwritten by a namelist entry)) 
  !    criterion 2,4,5,9,11 and 12.
  !    Criterion 10 (element perhaps in error) also looks suspect. 
  ! 2:  element rejected by a selection process (thinning)
  ! 4:  element rejected by either the background check or the qc-var
  ! 5:  element rejected because it is on a blacklist
  ! 9:  doubtful element (this flag will frequently be activated due to
  !       the use of outdated black lists (70's) for in particular surface 
  !       wind observations)
  ! 10: element perhaps erroneous 
  ! 11: element in error
  ! 12: element that exceeds climatological limits (it would seem prudent
  !       to reject observations that are beyond climatological limits   
  !       because there is no qc-var like mechanism in the EnKF to 
  !       gradually decrease the weight given to observations). From an 
  !       inspection by Pierre Koclas (February 2002) it appeared that 
  !       observations with marker 12 that pass the background check and
  !       the qc-var are only barely beyond the first (tightest) set of 
  !       rejection limits (indicating that they are perhaps correct).
  !       There were no observations beyond the second (widest) set of
  !       rejection criteria.
  ! The criteria are listed in a document entitled Traitement des
  !  observations dans l'analyse 3d-var eta, version TT, cmda, nov 2001,
  !  obtained from Gilles Verner). The quality control procedures do
  !  of course continue to evolve and the use of the markers by the 
  !  EnKF will have to be reinvestigated periodically. 
 
  NLISTFLG(1)=2
  NLISTFLG(2)=4
  NLISTFLG(3)=5
  NLISTFLG(4)=-5
 
  call obs_class_initialize('ENKF')
  call obs_initialize(cma, MXSTNTOTAL, MXOBSTOTAL)
      
  ! convert burp input data to cma format
  call selectb(cma,nfiles,nbegintyp,nendtyp,cfamtyp,cfilnam,nulout)

  call obs_print(cma,nobsout)

  ! write the cma info to unformatted files
  ! note that the body information is written in the 
  ! order that it will be used by sekfeta.f

  call obs_write(cma,hx,zero_ens,ncmahdr,ncmabdy,ncmahx,ncmadim)

  close(ncmahdr)
  close(ncmabdy)
  close(ncmadim)

  filestat='NEW'
  fileform='FORMATTED'
  call openfile(nbrpform,brpform,filestat,fileform)
  write(nbrpform,*) qcvar
  close(nbrpform)

  call obs_finalize(cma)

  stop

contains

SUBROUTINE SELECTB(obsdat,nfiles,nbegintyp,nendtyp,cfamtyp,cfilnam, &
         nulout)
  !
  !      PURPOSE: READ CMC BURP FILES FILL UP CMA FILE
  !
  !    ARGUMENTS:
  !            obsdat   - obsdat-file object
  !            input:  nfiles   - number of input files
  !                    nulout   - unit number for standard output
  !                    cfamtyp  - type of observation family
  !                    cfilnam  - filename
  !            output: nbegintyp - first location for family
  !                    nendtyp   - last location for family
  !
  !       AUTHOR: P. KOCLAS(CMC CMDA)
  !       Revision:
  !     NOTE:
  !     BURP FILES ARE ASSUMED TO BE PRESENT IN CURRENT WORKING DIRECTORY
  !
  use ObsSpaceData_mod
  use obsUtil_mod
  use burpread_mod

  IMPLICIT NONE

  type (struct_obs), intent(inout) :: obsdat
  integer, intent(in) :: nfiles,nulout
  integer, intent(out), dimension(:) :: nbegintyp,nendtyp
  character(len=2),   intent(in), dimension(:) :: cfamtyp
  character(len=256), intent(in), dimension(:) :: cfilnam

  INTEGER :: IER,IBEG, iend, nstn1,nstn2
  INTEGER :: EXDB,EXFIN
  INTEGER :: JO,START
  logical :: obs_full
  REAL*4,   PARAMETER  :: PPMIS=-999.
!
  INTEGER NVALS,J
!
! currently supported families of data 'UA' 'AI' 'SC' 'SF' 'SW' 'TO'
!
  WRITE(*,*)' '
  WRITE(*,*)'================================================='
  WRITE(*,*)'                SELECT READBURP BEGIN                     '
  WRITE(*,*)'================================================='
  WRITE(*,*)' '
  IER=EXDB('SELECT READBURP','DEBUT','NON')


  NVALS=NFILES
  ibeg=obs_numbody(obsdat)
  start=obs_numheader(obsdat) +1
  DO J =1,NVALS
!
    call obs_status(obsdat, obs_full, nstn1, ibeg, nulout)
    IBEG=ibeg + 1
    call brpr_readBurp(obsdat,cfamtyp(J),cfilnam(J),J)
    call obs_status(obsdat, obs_full, nstn2, iend, nulout)

    IF ( IBEG .LE. iend ) THEN
      nbegintyp(J)=IBEG
      nendtyp(J) =iend
    ELSE
      nbegintyp(J)=-999
      nendtyp(J)  =-999
    ENDIF
    IF ( trim(CFAMTYP(J)) .ne. 'TO') THEN
      call obsu_windDirectionToUV(obsdat,nstn1+1,nstn2,PPMIS)
      call obsu_adjustHumGZ(obsdat,nstn1+1,nstn2)
      call obsu_computeVertCoordSurfObs(obsdat,nstn1+1,nstn2)
    ENDIF
    start=obs_numHeader(obsdat) +1
    DO JO=nstn1+1,nstn2
      call obs_headSet_i(obsdat,OBS_OTP,JO,J)
      call obs_setFamily(obsdat,trim(cfamtyp(J)),JO)
    END DO

    DO JO=IBEG,IEND
      call obs_bodySet_r(obsdat,OBS_OMA ,JO,PPMIS)
      call obs_bodySet_r(obsdat,OBS_OMP ,JO,PPMIS)
      call obs_bodySet_r(obsdat,OBS_OMP6,JO,PPMIS)
      call obs_bodySet_r(obsdat,OBS_OMA0,JO,PPMIS)
      call obs_bodySet_r(obsdat,OBS_OER ,JO,PPMIS)
      call obs_bodySet_r(obsdat,OBS_HPHT,JO,PPMIS)
      call obs_bodySet_r(obsdat,OBS_HAHT,JO,PPMIS)
    END DO
  END DO

  WRITE(*,*) '  readburp obs_numheader(obsdat)', obs_numheader(obsdat)
  WRITE(*,*) '  readburp obs_numbody(obsdat)  ', obs_numbody  (obsdat)
!
  WRITE(*,*)' '
  WRITE(*,*)'================================================='
  WRITE(*,*)'                SELECT READBURP     END                   '
  WRITE(*,*)'================================================='
  WRITE(*,*)' '
  IER=EXFIN('SELECT READBURP','FIN','NON')
  RETURN

END SUBROUTINE SELECTB

end program prepcma

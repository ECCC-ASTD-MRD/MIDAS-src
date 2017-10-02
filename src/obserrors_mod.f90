!-------------------------------------- LICENCE BEGIN ------------------------------------
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

!--------------------------------------------------------------------------
!! MODULE obsErrors (prefix="oer")
!!
!! *Purpose*: Subroutines to set up the observation error standard deviations.
!!
!--------------------------------------------------------------------------
module obsErrors_mod
  use MathPhysConstants_mod
  use obsSpaceData_mod
  use tovs_nl_mod
  use chem_obserrors_mod
  use codtyp_mod
  use bufr
  use utilities_mod
  use EarthConstants_mod
  use gps_mod
  use columnData_mod
  implicit none
  private

  ! public procedures
  ! -----------------

  public :: oer_setObsErrors, oer_SETERRGPSGB, oer_SETERRGPSRO, oer_setInterchanCorr

  ! TOVS OBS ERRORS
  ! ---------------

  real(8) :: toverrst(tvs_maxChannelNumber,tvs_maxNumberOfSensors)

  !
  ! CONVENTIONAL OBS ERRORS
  ! -----------------------

  real(8) :: xstd_ua_ai_sw(20,11)
  real(8) :: xstd_sf(9,4)
  real(8) :: xstd_pr(2)
  real(8) :: xstd_sc(1)

  save
  character(len=48) :: obserrorMode

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

   SUBROUTINE oer_setInterchanCorr()
!
!  s/r oer_setInterchanCorr : Setup of interchannel observation errors correlations
!           S.  Heilliette
!            - 04/2017 extracted from tvslin_setupallo
      Use rmatrix_mod

      IMPLICIT NONE
!implicits

      INTEGER ::  ISENS

!-----------------------------------------------------------------------

      if (tvs_nsensors == 0) then
         write(*,*) 'oer_setInterchanCorr: tvs_nsensors is zero, skipping setup'
         return
      end if

! Initialization of the correlation matrices
      call rmat_init(tvs_nsensors,tvs_nobtov)
      if (rmat_lnondiagr) then
         do isens = 1, tvs_nsensors
            if (tvs_isReallyPresent(isens) ) call rmat_readCMatrix(tvs_listSensors(:,isens), isens, tvs_ichan(1:tvs_nchan(isens),isens)  )
         end do
      end if

   END SUBROUTINE oer_setInterchanCorr

  subroutine oer_setObsErrors(lobsSpaceData, obserrorMode_in)
    !
    ! s/r set_obsErrors -SET OBSERVATION ERROR FOR ALL DATA 
    !
    ! Author  : S. Laroche  February 2014
    ! Revision: 
    !          Y. Rochon, M. Sitwell and P. Du, March 2017
    !          - Added calls to chm_read_obs_err_stddev and
    !            chm_dealloc_obs_err_stddev
    !
    ! Purpose: read and set observation errors (from former sucovo subroutine).
    !
    !

    type(struct_obs) :: lobsSpaceData
    character(len=*), intent(in) :: obserrorMode_in

    !
    ! Setup Mode
    !
    obserrorMode = obserrorMode_in

    !
    ! Read in the observation stddev errors for radiance data
    !
    if (obs_famexist(lobsSpaceData,'TO')) then
       call read_obs_erreurs_tovs
    else 
       write(*,*) "oer_setObsErrors: No brightness temperature observations found."
    end if

    !
    ! Read in the observation stddev errors for conventional data
    !
    if (obs_famExist(lobsSpaceData,'UA').or.obs_famExist(lobsSpaceData,'AI').or.obs_famExist(lobsSpaceData,'SW').or. &
       obs_famExist(lobsSpaceData,'SF').or.obs_famExist(lobsSpaceData,'GP').or.obs_famExist(lobsSpaceData,'SC').or. &
       obs_famExist(lobsSpaceData,'PR') ) then
       call read_obs_erreurs_conv
    else
       write(*,*) "oer_setObsErrors: No conventional weather observations found."
    end if

    !
    ! Read in the observation stddev error for constituent data
    !
    if (obs_famexist(lobsSpaceData,'CH')) then
       call chm_read_obs_err_stddev
    else
       write(*,*) "oer_setObsErrors: No CH observations found."
    end if
  
    !
    ! Set obs error information in obsSpaceData object
    !
    call fill_obs_erreurs(lobsSpaceData)

    !
    ! Deallocate temporary storage
    !
    if (obs_famExist(lobsspaceData,'CH')) call chm_dealloc_obs_err_stddev

  end subroutine oer_setObsErrors

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  subroutine read_obs_erreurs_tovs
    !
    ! s/r read_obs_erreurs_tovs
    ! - Read the observation erreur statistics and utilization flag for TOVS processing.
    !
    !
    ! Author  : J. Halle !CMDA/AES  May 08, 1996
    !
    ! Revision 01  : J. Halle !CMDA/AES  Oct 1999
    !               - change file name to stats_tovs
    !
    ! Revision 002  : J. Halle !CMDA/AES  dec 2000
    !                adapt to TOVS level 1b.
    !
    ! Revision 003  : J. Halle !CMDA/SMC  may 2002
    !                adapt to RTTOV-7.
    !
    ! Revision 004  : J. Halle !CMDA/SMC  sept 2006
    !                adapt to RTTOV-8.
    !
    ! Revision 005  : A. Beaulne !CMDA/SMC  fevr 2007
    !                adapt utilization flag for AIRS channels
    !
    ! Revision 006  : S. Heilliette
    !                adapt utilization flag for IASI channels
    !
    ! Revision 007: S. Macpherson ARMA, Feb 2013
    !              - add NPP/ATMS codtyp=192
    !
    ! Revision 008  : S. Macpherson      apr 2013
    !              - adapt for new-format stats_tovs file
    !
    ! Revision 008  : S. Laroche      Mar 2014
    !              - upgrade of former f77 subroutine sutovst
    !              - f90 conversion and cleanup 
    !
    !    -------------------
    !    Purpose: Read the observation error statistics and
    !             utilization flag for TOVS processing. This information
    !             resides on an ASCII file and is read using a free format.


    use hirchannels_mod
    use mpivar_mod
    use tovs_nl_mod
    use rmatrix_mod
    implicit none

    integer,external  :: FNOM, FCLOS
    integer  :: IER, ILUTOV, JI, JJ, JK, JL, JM, I, IPOS1, IPOS2, INUMSAT, ISAT, IPLF
    integer, dimension(tvs_maxNumberOfSensors)         :: IPLATFORM, ISATID, IINSTRUMENT, NUMCHN, NUMCHNIN
    integer, dimension(tvs_maxChannelNumber,tvs_maxNumberOfSensors) :: IUTILST, ICHN, ICHNIN

    real :: ZDUM

    real(8), dimension(tvs_maxChannelNumber,2,tvs_maxNumberOfSensors) :: TOVERRIN

    character (len=132) :: CLDUM,CPLATF,CINSTR

    WRITE(*,'(//,10x,"-read_obs_erreurs_tovs: reading observation error statistics required for TOVS processing")')

    !
    !    1. Initialize
    !       ----------
    !
    DO JL = 1, tvs_maxNumberOfSensors  
       DO JI = 1, tvs_maxChannelNumber
          TOVERRST(JI,JL) = 0.0D0
          TOVERRIN(JI,1,JL) = 0.0D0
          TOVERRIN(JI,2,JL) = 0.0D0
          IUTILST   (JI,JL) = 0
       END DO
    END DO

    DO JL = 1, tvs_maxNumberOfSensors
       IPLATFORM(JL) = 0
       NUMCHN(JL) = 0
       NUMCHNIN(JL) = 0
       DO JI = 1, tvs_maxChannelNumber
          ICHN(JI,JL) = 0
          ICHNIN(JI,JL) = 0
       END DO
    END DO

    if (tvs_nobtov == 0) return

    !
    !     2. Open the file
    !        -------------
    !
    ilutov = 0
    IER =  FNOM(ILUTOV,'stats_tovs','SEQ+FMT',0)
    IF(IER.LT.0)THEN
       WRITE ( *, '(" read_obs_erreurs_tovs: Problem opening ","file stats_tovs ")' )
       CALL utl_abort ('read_obs_erreurs_tovs')
    END IF

    !
    !     3. Read number of satellites
    !        -------------------------
    !
    READ (ILUTOV,*)
    READ (ILUTOV,*) INUMSAT
    READ (ILUTOV,*)

    !
    !     4. Read the satellite identification, the number of channels,
    !        the observation errors and the utilization flags
    !        ----------------------------------------------------------
    !
    WRITE(*,'(5X,"read_obs_erreurs_tovs: Reading stats_tovs file: "//)')

    DO JL = 1, INUMSAT

       READ (ILUTOV,*)
       READ (ILUTOV,'(A)') CLDUM
       WRITE(*,'(A)') CLDUM
       CINSTR=CLDUM
       call split(CINSTR," ",CPLATF)
       Write(*,*) "CINSTR: ",CINSTR
       Write(*,*) "CPLATF: ",CPLATF
       READ (ILUTOV,*)
       READ (ILUTOV,*) ISATID(JL), NUMCHNIN(JL)

       DO JI = 1, 3
          READ (ILUTOV,*)
       END DO

       IPLATFORM(JL) =  tvs_getPlatformId(CPLATF)

       IF ( IPLATFORM(JL) .EQ. -1 ) THEN
          WRITE ( *, '(" read_obs_erreurs_tovs: Unknown platform!"/)' )
          CALL utl_abort ('read_obs_erreurs_tovs')
       END IF

       IINSTRUMENT(JL) = tvs_getInstrumentId(CINSTR)

       IF ( IINSTRUMENT(JL) .EQ. -1 ) THEN
          WRITE ( *, '(" read_obs_erreurs_tovs: Unknown instrument!"/)' )
          CALL utl_abort ('read_obs_erreurs_tovs')
       END IF

       DO JI = 1, NUMCHNIN(JL)
          READ (ILUTOV,*) ICHNIN(JI,JL), TOVERRIN(ICHNIN(JI,JL),1,JL), TOVERRIN(ICHNIN(JI,JL),2,JL), IUTILST(ICHNIN(JI,JL),JL), ZDUM
       END DO
       READ (ILUTOV,*)

    END DO

    !
    !   Select input error to use: if ANAL mode, use ERRANAL (JJ=2);
    !   otherwise use ERRBGCK (JJ=1)
    !
    IF ( trim(obserrorMode) == 'analysis' .or. trim(obserrorMode) == 'FSO' ) THEN
       JJ = 2
    ELSE
       JJ = 1
    END IF

    !
    !   Fill the observation error array TOVERRST
    !
    WRITE(*,'(5X,"read_obs_erreurs_tovs: Fill error array TOVERRST: "//)')
    DO JM= 1, INUMSAT
       DO JL = 1, tvs_nsensors
          IF ( tvs_platforms (JL) .EQ. IPLATFORM(JM) .AND. tvs_satellites(JL) .EQ. ISATID(JM) ) THEN
             IF ( tvs_instruments (JL) .EQ. IINSTRUMENT(JM) ) THEN
                NUMCHN(JL)=NUMCHNIN(JM)
                DO JI = 1, tvs_maxChannelNumber
                   TOVERRST(JI,JL) = TOVERRIN(JI,JJ,JM)
                   ICHN(JI,JL) = ICHNIN(JI,JM)
                END DO
                IF ( (trim(obserrorMode) == 'analysis' .or. trim(obserrorMode) == 'FSO') .and. rmat_lnondiagr) then
                  call rmat_setFullRMatrix ( TOVERRST(:,JL), JL, tvs_channelOffset(JL) )
                end if
             END IF
          END IF
       END DO
    END DO

    !
    !  Check that oberservation error statistics have been defined for
    !  all the satellites specified in the namelist.
    !
    DO JL = 1, tvs_nsensors
       IPLF = utl_findArrayIndex( IPLATFORM  , INUMSAT, tvs_platforms  (JL) )
       ISAT =  utl_findArrayIndex( ISATID     , INUMSAT, tvs_satellites (JL) )
       IF ( IPLF .EQ. 0 .OR. ISAT .EQ. 0 ) THEN
          WRITE ( *, '(" read_obs_erreurs_tovs: Observation errors not ","defined for sensor # ", I3)' ) JL
          CALL utl_abort ('read_obs_erreurs_tovs')
       END IF
       IF ( NUMCHN(JL) .EQ. 0 ) THEN 
          WRITE ( *, '(" read_obs_erreurs_tovs: Problem setting errors ","for sensor # ", I3)' ) JL
          CALL utl_abort ('read_obs_erreurs_tovs')
       END IF
    END DO

    !
    !   Utilization flag for AIRS,IASI and CrIS channels (bgck mode only)
    !
    IF ( trim(obserrorMode) == 'bgckIR' ) THEN
       DO JM= 1, INUMSAT
          IF (   IPLATFORM(JM) .EQ. 9  .AND. &
               IINSTRUMENT(JM) .EQ. 11 ) THEN

             call hir_set_assim_chan("AIRS",IUTILST(ICHNIN(1:NUMCHNIN(JM),JM),JM))
          END IF
          IF (   IPLATFORM(JM) .EQ. 10 .AND. &
               IINSTRUMENT(JM) .EQ. 16 ) THEN
             call hir_set_assim_chan("IASI",IUTILST(ICHNIN(1:NUMCHNIN(JM),JM),JM))
          END IF
          IF (   IPLATFORM(JM) .EQ. 17 .AND. &
               IINSTRUMENT(JM) .EQ. 27 ) THEN
             call hir_set_assim_chan("CRIS",IUTILST(ICHNIN(1:NUMCHNIN(JM),JM),JM))
          END IF
       END DO
    END IF

    !
    !    5. Print out observation errors for each sensor
    !       --------------------------------------------
    !
    IF(MPI_MYID.eq.0) THEN
       WRITE(*,'(//1X,"Radiance observation errors read from file")')
       WRITE(*,'(  1X,"------------------------------------------")')
       DO JL = 1, tvs_nsensors
          WRITE(*,'(/1X,"SENSOR #",I2,". Platform: ",A,"Instrument: ",A)') &
               JL, tvs_satelliteName(JL), tvs_instrumentName(JL)
          WRITE(*,'(1X,"Channel",5X,"  error   ")')
          DO JI = 1, NUMCHN(JL)
             WRITE (*,'(1X,I7,1(5X,F10.5))') ICHN(JI,JL),TOVERRST(ICHN(JI,JL),JL)
          END DO
       END DO
    END IF

    !
    !    6. Close the file
    !       --------------
    !
    IER = FCLOS(ILUTOV)
    IF(IER.NE.0)THEN
       CALL utl_abort ('read_obs_erreurs_tovs')
    END IF

  contains

    subroutine compact(str)
      ! Code from Benthien's module: http://www.gbenthien.net/strings/index.html
      ! Converts multiple spaces and tabs to single spaces; deletes control characters;
      ! removes initial spaces.

      character(len=*):: str
      character(len=1):: ch
      character(len=len_trim(str)):: outstr
      integer isp,k,lenstr,i,ich

      str=adjustl(str)
      lenstr=len_trim(str)
      outstr=' '
      isp=0
      k=0

      do i=1,lenstr
         ch=str(i:i)
         ich=iachar(ch)

         select case(ich)
         case(9,32)    ! space or tab character         
            if(isp==0) then
               k=k+1
               outstr(k:k)=' '
            end if
            isp=1
         case(33:)              ! not a space, quote, or control character
            k=k+1
            outstr(k:k)=ch
            isp=0
         end select

      end do

      str=adjustl(outstr)

    end subroutine compact

    subroutine split(str,delims,before)
      ! Code extracted from Benthien's module: http://www.gbenthien.net/strings/index.html
      ! Routine finds the first instance of a character from 'delims' in the
      ! the string 'str'. The characters before the found delimiter are
      ! output in 'before'. The characters after the found delimiter are
      ! output in 'str'. 

      character(len=*) :: str,delims,before
      character :: ch,cha
      integer lenstr,i,k,ipos,iposa
      str=adjustl(str)
      call compact(str)
      lenstr=len_trim(str)

      if(lenstr == 0) return ! string str is empty
      k=0
      before=' '
      do i=1,lenstr
         ch=str(i:i)

         ipos=index(delims,ch)

         if(ipos == 0) then ! character is not a delimiter
            k=k+1
            before(k:k)=ch
            cycle
         end if
         if(ch /= ' ') then ! character is a delimiter that is not a space
            str=str(i+1:)
            exit
         end if

         cha=str(i+1:i+1)  ! character is a space delimiter
         iposa=index(delims,cha)
         if(iposa > 0) then   ! next character is a delimiter 
            str=str(i+2:)
            exit
         else
            str=str(i+1:)
            exit
         end if
      end do
      if(i >= lenstr) str=''

      str=adjustl(str) ! remove initial spaces

      return

    end subroutine split

  end subroutine read_obs_erreurs_tovs

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  subroutine read_obs_erreurs_conv
    !
    ! s/r read_obs_erreurs_conv -READ OBSERVATION ERROR OF CONVENTIONAL DATA 
    !
    ! Author  : S. Laroche  February 2014
    ! Revision: 
    !          
    ! Purpose: read observation errors (modification of former readcovo subroutine).
    !

    implicit none

    integer :: FNOM, FCLOS
    integer :: IERR, JLEV, JELM, icodtyp, nulstat
    logical :: LnewExists

    character (len=128) :: ligne

    EXTERNAL FNOM,FCLOS
    !
    !     CHECK THE EXISTENCE OF THE NEW FILE WITH STATISTICS
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    INQUIRE(FILE='obserr',EXIST=LnewExists)
    IF (LnewExists )then
       WRITE(*,*) '--------------------------------------------------------'
       WRITE(*,*) 'read_obs_errors_conv: reads observation errors in obserr'
       WRITE(*,*) '--------------------------------------------------------'
    else
       CALL utl_abort('read_obs_errors_conv: NO OBSERVATION STAT FILE FOUND!!')     
    END IF
    !
    !     Read observation errors from file obserr for conventional data
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    NULSTAT=0
    IERR=FNOM(NULSTAT,'obserr','SEQ',0)
    IF ( IERR .EQ. 0 ) THEN
       write(*,*) 'read_obs_errors_conv: File =  ./obserr'
       write(*,*) ' opened as unit file ',nulstat
       open(unit=nulstat, file='obserr', status='OLD')
    ELSE
       CALL utl_abort('read_obs_errors_conv:COULD NOT OPEN FILE obserr!!!')
    END IF

    write(*, '(A)') ' '

    do jlev = 1,3
       read(nulstat, '(A)') ligne
       write(*, '(A)') ligne
    enddo

    do jlev = 1, 19
       read(nulstat, * ) (xstd_ua_ai_sw(jlev,jelm), jelm=1,11)
       write(*, '(f6.0,10f6.1)' )  (xstd_ua_ai_sw(jlev,jelm), jelm=1,11)
    enddo

    do jlev = 1,5
       read(nulstat, '(A)') ligne
       write(*, '(A)') ligne
    enddo

    read(nulstat, * ) xstd_pr(1),xstd_pr(2)
    write(*, '(2f6.1)' )  xstd_pr(1),xstd_pr(2)

    do jlev = 1,5
       read(nulstat, '(A)') ligne
       write(*, '(A)') ligne
    enddo

    read(nulstat, * ) xstd_sc(1)
    write(*, '(f8.3)' )  xstd_sc(1)

    read(nulstat, '(A)') ligne
    write(*, '(A)') ligne

    do icodtyp = 1,9
       do jlev = 1,4
          read(nulstat, '(A)') ligne
          write(*, '(A)') ligne
       enddo
       read(nulstat, * ) (xstd_sf(icodtyp,jelm), jelm=1,4)
       write(*, '(f6.2,2f6.1,f8.3)' )  (xstd_sf(icodtyp,jelm), jelm=1,4)
    enddo

    write(*, '(A)') ' '

    CLOSE(UNIT=NULSTAT)
    IERR=FCLOS(NULSTAT)    

  end subroutine read_obs_erreurs_conv

  !-----------------------------------------------------------------------------------------

  subroutine fill_obs_erreurs(lobsSpaceData)
    !
    ! s/r fill_obs_erreurs -FILL OBSERVATION ERRORS IN lobsSpaceData
    !
    ! Author  : S. Laroche  February 2014
    !
    ! Revision: 
    !           Y. Rochon and M. Sitwell, ARQI/AQRD, Nov 2014 - Feb 2015
    !           - Added consideration of 'CH' family
    !           - Removal of 'OZ' family
    !           - Added availability of CSTNID, obs_err_stddev
    !           M. Sitwell, ARQI/AQRD, Apr 2015 
    !           - Added IDATE and ITIME
    !         
    ! Purpose: read observation errors (modification of former readcovo subroutine).
    !
    !--------------------------------------------------------------------------------------

    implicit none

    type(struct_obs) :: lobsSpaceData

    integer :: JN, JI, INDEX_BODY, INDEX_HEADER, ITYP, IFLG, IASS, IDATA, IDATEND, IDBURP
    integer :: ISAT, ICHN, IPLATF, INSTR, IPLATFORM, INSTRUM
    integer :: ILEV,ISTNID,JINT,NLEV,IDATE,ITIME
    integer :: ielem,icodtyp,header_prev
    real(8) :: ZLAT, ZLON, ZLEV, ZVAL, zwb, zwt, obs_err_stddev
    logical :: IFIRST

    CHARACTER(len=2) :: CFAM
    CHARACTER(len=12) :: CSTNID

    !
    !     ========================================================================== 
    !
    WRITE(*,'(10X,"Fill_obs_errors")')
    WRITE(*,'(10X,"-----------------",/)')
    WRITE(*,'(10X,"***********************************")')
    WRITE(*,'(10X,"Fill_obs_errors:",/)')
    WRITE(*,'(10X,"***********************************")')

    !
    !     SET STANDARD DEVIATION ERRORS FOR EACH DATA FAMILY
    !     ---------------------------------------------------
    !
    DO INDEX_HEADER = 1, obs_numheader(lobsSpaceData)

       IDATA   = obs_headElem_i(lobsSpaceData,OBS_RLN,INDEX_HEADER)
       IDATEND = obs_headElem_i(lobsSpaceData,OBS_NLV,INDEX_HEADER) + IDATA - 1
       CFAM    = obs_getFamily(lobsSpaceData,INDEX_HEADER)
       ZLAT    = obs_headElem_r(lobsSpaceData,OBS_LAT,INDEX_HEADER)
       ZLON    = obs_headElem_r(lobsSpaceData,OBS_LON,INDEX_HEADER)
       IDBURP  = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
       IPLATF  = obs_headElem_i(lobsSpaceData,OBS_SAT,INDEX_HEADER)
       INSTR   = obs_headElem_i(lobsSpaceData,OBS_INS,INDEX_HEADER)
       CSTNID  = obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER)
       IDATE  = obs_headElem_i(lobsSpaceData,OBS_DAT,INDEX_HEADER) 
       ITIME  = obs_headElem_i(lobsSpaceData,OBS_ETM,INDEX_HEADER)
       
       NLEV = IDATEND-IDATA+1
       
       DO INDEX_BODY  = IDATA, IDATEND

          ITYP  = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY)
          IFLG  = obs_bodyElem_i(lobsSpaceData,OBS_FLG,INDEX_BODY)
          IASS  = obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY)
          ZVAL  = obs_bodyElem_r(lobsSpaceData,OBS_VAR,INDEX_BODY)

          IF ( IASS .EQ. 1 ) THEN

             !***********************************************************************
             !                           TOVS DATA
             !***********************************************************************

             IF ( CFAM .EQ. 'TO' ) THEN

                IF ( ITYP .EQ. BUFR_NBT1 .OR. &
                     ITYP .EQ. BUFR_NBT2 .OR. &
                     ITYP .EQ. BUFR_NBT3     )THEN

                   ICHN = NINT(obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY))
                   CALL tvs_mapSat(IPLATF,IPLATFORM,ISAT)
                   CALL tvs_mapInstrum(INSTR,INSTRUM)

                   DO JN = 1, tvs_nsensors
                      IF ( IPLATFORM ==  tvs_platforms(JN)  .AND. &
                           ISAT      ==  tvs_satellites(JN) .AND. &
                           INSTRUM   == tvs_instruments(JN)      ) THEN
                         call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,TOVERRST(ICHN,JN))
                      END IF
                   END DO

                END IF

                !***********************************************************************
                !                      RADIOSONDE DATA
                !***********************************************************************

             ELSE IF ( CFAM .EQ. 'UA' ) THEN

                ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)

                IF ( (ITYP .EQ. BUFR_NEUS) .OR. (ITYP .EQ. BUFR_NEVS) )THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(1,4))
                ELSE IF (ITYP .EQ. BUFR_NETS) THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(1,2))
                ELSE IF (ITYP .EQ. BUFR_NESS) THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(1,3))
                ELSE IF (ITYP .EQ. BUFR_NEPS ) THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(1,1))
                ELSE IF (ITYP .EQ. BUFR_NEPN ) THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(1,1))
                ELSE

                   if ( (ITYP .EQ. BUFR_NEUU) .OR. (ITYP .EQ. BUFR_NEVV) ) then
                      ielem = 4
                   else if (ITYP .EQ. BUFR_NETT) then
                      ielem = 2
                   else if (ITYP .EQ. BUFR_NEES) then
                      ielem = 3
                   else if (ITYP .EQ. BUFR_NEGZ) then
                      ielem = 5
                   endif

                   if ( (ZLEV*MPC_MBAR_PER_PA_R8) >= xstd_ua_ai_sw(1,1) ) then

                      call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_ua_ai_sw(1,ielem))

                   else if ( (ZLEV*MPC_MBAR_PER_PA_R8) <= xstd_ua_ai_sw(19,1) ) then

                      call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_ua_ai_sw(19,ielem))

                   else

                      do jn = 1,18
                         if ( (ZLEV*MPC_MBAR_PER_PA_R8) >= xstd_ua_ai_sw(jn+1,1) ) exit
                      end do

                      zwb = log((ZLEV*MPC_MBAR_PER_PA_R8)/xstd_ua_ai_sw(JN,1)) / log(xstd_ua_ai_sw(JN+1,1)/xstd_ua_ai_sw(JN,1))
                      zwt = 1.0D0 - zwb

                      call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,zwt*xstd_ua_ai_sw(JN,ielem) + zwb*xstd_ua_ai_sw(JN+1,ielem))

                   endif

                END IF

                !***********************************************************************
                !                          AMV, AIREP, AMDAR DATA
                !***********************************************************************

             ELSE IF ( CFAM .EQ. 'AI'.OR. CFAM .EQ. 'SW') THEN

                ZLEV=obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)

                IF ( IDBURP .EQ. 188 ) THEN ! AMV
                   if ( (ITYP .EQ. BUFR_NEUU) .OR. (ITYP .EQ. BUFR_NEVV) ) then
                      ielem = 11
                   endif

                ELSE IF (IDBURP .EQ. 128 ) THEN ! AIREP
                   if ( (ITYP .EQ. BUFR_NEUU) .OR. (ITYP .EQ. BUFR_NEVV) ) then
                      ielem = 7
                   else if ( ITYP .EQ. BUFR_NETT ) then
                      ielem = 6
                   endif

                ELSE IF (IDBURP .EQ. 42 .OR. IDBURP .EQ. 157 .OR. IDBURP .EQ. 177) THEN ! AMDAR
                   if ( (ITYP .EQ. BUFR_NEUU) .OR. (ITYP .EQ. BUFR_NEVV) ) then
                      ielem = 10
                   else if ( ITYP .EQ. BUFR_NETT ) then
                      ielem = 8
                   else if ( ITYP .EQ. BUFR_NEES ) then
                      ielem = 9
                   endif

                END IF

                if ( (ZLEV*MPC_MBAR_PER_PA_R8) >= xstd_ua_ai_sw(1,1) ) then
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_ua_ai_sw(1,ielem))
                else if ( (ZLEV*MPC_MBAR_PER_PA_R8) <= xstd_ua_ai_sw(19,1) ) then
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_ua_ai_sw(19,ielem))
                else

                   do jn = 1,18
                      if ( (ZLEV*MPC_MBAR_PER_PA_R8) >= xstd_ua_ai_sw(jn+1,1) ) exit
                   enddo

                   zwb = log((ZLEV*MPC_MBAR_PER_PA_R8)/xstd_ua_ai_sw(JN,1)) / log(xstd_ua_ai_sw(JN+1,1)/xstd_ua_ai_sw(JN,1))
                   zwt = 1.0D0 - zwb 

                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,zwt*xstd_ua_ai_sw(JN,ielem) + zwb*xstd_ua_ai_sw(JN+1,ielem))

                endif

                !***********************************************************************
                !                         SURFACE DATA
                !***********************************************************************

             ELSE IF ( CFAM .EQ. 'SF' ) THEN
                icodtyp = 1   ! Default values
                IF ( IDBURP .EQ. 12  ) icodtyp = 2   ! SYNOP
                IF ( IDBURP .EQ. 13  ) icodtyp = 3   ! SHIP NON-AUTOMATIQUE
                IF ( IDBURP .EQ. 14  ) icodtyp = 4   ! DRIBU
                IF ( IDBURP .EQ. 18  ) icodtyp = 5   ! DRIFTER
                IF ( IDBURP .EQ. 145 ) icodtyp = 6   ! STATION AUTOMATIQUE
                IF ( IDBURP .EQ. 146 ) icodtyp = 7   ! ASYNOP
                IF ( IDBURP .EQ. 147 ) icodtyp = 8   ! ASHIP
                IF ( (ITYP .EQ. BUFR_NEUU) .OR. (ITYP .EQ. BUFR_NEVV) .OR. &
                     (ITYP .EQ. BUFR_NEGZ) .OR. (ITYP .EQ. BUFR_NETT) .OR. (ITYP .EQ. BUFR_NEES) ) icodtyp = 9  ! Others

                IF ( (ITYP .EQ. BUFR_NEUS) .OR. (ITYP .EQ. BUFR_NEVS) )THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(icodtyp,4))
                ELSE IF (ITYP .EQ. BUFR_NETS) THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(icodtyp,2))
                ELSE IF (ITYP .EQ. BUFR_NESS) THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(icodtyp,3))
                ELSE IF (ITYP .EQ. BUFR_NEPS ) THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(icodtyp,1))
                ELSE IF (ITYP .EQ. BUFR_NEPN ) THEN
                   if(icodtyp == 2  .or. icodtyp == 7) then
                      call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(1,1))
                   else
                      call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(icodtyp,1))
                   endif
                ELSE IF ( (ITYP .EQ. BUFR_NEUU) .OR. (ITYP .EQ. BUFR_NEVV) )THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(icodtyp,4))
                ELSE IF (ITYP .EQ. BUFR_NEGZ) THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(icodtyp,1))
                ELSE IF (ITYP .EQ. BUFR_NETT) THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(icodtyp,2))
                ELSE IF (ITYP .EQ. BUFR_NEES) THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(icodtyp,3))
                END IF

                !***********************************************************************
                !                             GPS RO DATA
                !***********************************************************************

             ELSE IF ( CFAM .EQ. 'RO' ) THEN
                !     
                !                 Process only refractivity data (codtyp 169)
                !
                IF ( obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER) .EQ. 169 ) THEN

                   IF ( ITYP .EQ. BUFR_NEPS ) THEN
                      call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY, 50.D0)
                   END IF
                   IF ( ITYP .EQ. BUFR_NETT) THEN
                      call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY, 10.D0)
                   END IF
                   IF ( ITYP .EQ. BUFR_NERF) THEN
                      call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,1001.D0)
                   END IF
                   IF ( ITYP .EQ. BUFR_NEBD) THEN
                      call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,1001.D0)
                   END IF

                END IF

                !***********************************************************************
                !                          GB-GPS SFC MET DATA
                !***********************************************************************

                !              ERRORS ARE SET TO SYNO SFC OBS ERRORS FROM S/R SUCOVO
                !              AND WEIGHTED BY FACTOR YSFERRWGT FOR 3D-VAR FGAT OR 4D-VAR ASSIM.
                !              OF TIME-SERIES (YSFERRWGT = 1.0 FOR 3D THINNING) 
                !
             ELSE IF ( CFAM .EQ. 'GP' ) THEN

                IF ( ITYP .EQ. BUFR_NEPS ) THEN ! Psfc Error (Pa)
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(2,1))
                END IF
                IF ( ITYP .EQ. BUFR_NETS ) THEN ! Tsfc Error (K)
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(2,2))
                END IF
                IF ( ITYP .EQ. BUFR_NESS ) THEN ! T-Td Error (K)
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sf(2,3))
                END IF
                ! ZTD Error (m) (value is formal error, real error set later in s/r seterrgpsgb)
                ! If error is missing, set to dummy value (1 m).
                IF ( ITYP .EQ. BUFR_NEZD ) THEN
                   IF (obs_bodyElem_r(lobsSpaceData,OBS_OER,INDEX_BODY) .LE. 0.0D0) call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY, 1.0D0)
                END IF

                !***********************************************************************
                !        SCATTEROMETER, WIND PROFILER DATA
                !***********************************************************************

             ELSE IF ( CFAM .EQ. 'SC' ) THEN

                call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_sc(1))

             ELSE IF ( CFAM .EQ. 'PR' ) THEN

                ZLEV= obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY) - obs_headElem_r(lobsSpaceData,OBS_ALT,INDEX_HEADER)

                IF ( ZLEV .GE. 6000. ) THEN
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_pr(2))
                ELSE
                   call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY,xstd_pr(1))
                END IF
              
                !***********************************************************************
                !               CONSTITUENT DATA (OZONE AND OTHER CHEMICALS)
                !***********************************************************************

             ELSE IF ( CFAM .EQ. 'CH' ) THEN

                !        Process only retrieved constituent data
                !        Also, exclude BUFR_SCALE_EXPONENT element as a data value!
               
                if (idburp.eq.codtyp_get_codtyp('CHEMREMOTE').or.idburp.eq.codtyp_get_codtyp('CHEMINSITU')) then

                  ifirst = index_header.ne.header_prev
                  if (ifirst) then
                     header_prev = index_header
                 
                     ! Subtract the number of occurences of code BUFR_SCALE_EXPONENT from number of levels
                     do ji=0,nlev-1
                        if (obs_bodyElem_i(lobsSpaceData,OBS_VNM,IDATA+ji).eq.BUFR_SCALE_EXPONENT) then
                           nlev = nlev-1
                        end if
                     end do
                  end if
              
                  ZLEV   = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)

                  ilev = 0
                  do ji=IDATA,INDEX_BODY
                     if (obs_bodyElem_i(lobsSpaceData,OBS_VNM,ji).ne.BUFR_SCALE_EXPONENT) ilev = ilev+1
                  end do

                  obs_err_stddev = chm_get_obs_err_stddev(CSTNID,NLEV,ITYP,ZLAT,ZLON,IDATE,ITIME,ZVAL,ZLEV,ILEV,IFIRST)
                  call obs_bodySet_r(lobsSpaceData,OBS_OER,INDEX_BODY, obs_err_stddev)
              
                end if

             ELSE

                WRITE(*,*)' UNKNOWN DATA FAMILY:',CFAM

             END IF

             !***********************************************************************
             !              Check for case where error should have been set but was
             !              not. 3dvar will abort in this case.
             !***********************************************************************

             IF (obs_bodyElem_r(lobsSpaceData,OBS_OER,INDEX_BODY) .LE. 0.0D0) THEN

                WRITE(*,*)'  PROBLEM OBSERR VARIANCE FAM= ',CFAM

                WRITE(*,'(1X,"STNID= ",A10,"IDBURP= ",I5," LAT= ",F10.2," LON = ",F10.2)') &
                     obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER),    &
                     IDBURP,                                           &
                     ZLAT*MPC_DEGREES_PER_RADIAN_R8,                   &
                     ZLON*MPC_DEGREES_PER_RADIAN_R8

                WRITE(*,'(1X,"ELEMENT= ",I6," LEVEL= ",F10.2," OBSERR = ",E10.2)')         &
                     ITYP,                                             &
                     obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY), &
                     obs_bodyElem_r(lobsSpaceData,OBS_OER,INDEX_BODY)

                CALL utl_abort('fill_obs_erreurs: PROBLEM OBSERR VARIANCE.')

             END IF

          END IF ! end of IASS .EQ. 1

       END DO ! end of INDEX_BODY loop 

    END DO ! end of INDEX_HEADER loop

    WRITE(*,'(10X,"Fill_obs_errors")')
    WRITE(*,'(10X,"---------------",/)')

  end subroutine fill_obs_erreurs

  SUBROUTINE oer_SETERRGPSRO(lcolumnhr,lobsSpaceData)

    !**s/r SETERRGPSRO - Compute estimated errors for GPSRO observations
    IMPLICIT NONE
    !
    type(struct_columnData) :: lcolumnhr
    type(struct_obs)        :: lobsSpaceData
    !
    INTEGER INDEX_HEADER, IDATYP, INDEX_BODY, iProfile
    REAL*8 zLat, Lat, sLat
    REAL*8 zLon, Lon
    REAL*8 zAzm, Azm
    REAL*8, allocatable :: ZPP(:)
    REAL*8, allocatable :: ZDP(:)
    REAL*8, allocatable :: ZTT(:)
    REAL*8, allocatable :: ZHU(:)
    REAL*8, allocatable :: ZUU(:)
    REAL*8, allocatable :: ZVV(:)
    !
    REAL*8 DH,DDH
    INTEGER JL, IAZM, ISAT, JH, NGPSLEV, NWNDLEV
    REAL*8 zMT, Rad, Geo, zP0
    REAL*8 HNH1, SUM0, SUM1, ZMIN
    !
    LOGICAL  ASSIM, L1, L2, L3

    INTEGER NH, NH1
    TYPE(GPS_PROFILE)           :: PRF
    REAL*8       , allocatable :: H   (:),AZMV(:)
    REAL*8       , allocatable :: ZOBS(:),ZREF(:),ZOFF(:),ZERR(:), ZMHX(:)
    TYPE(GPS_DIFF), allocatable :: RSTV(:)

    WRITE(*,*)'ENTER SETERRGPSRO'
    !
    !     * 1.  Initializations
    !     *     ---------------
    !
    NGPSLEV=col_getNumLev(LCOLUMNHR,'TH')
    NWNDLEV=col_getNumLev(LCOLUMNHR,'MM')
    allocate(ZPP (NGPSLEV))
    allocate(ZDP (NGPSLEV))
    allocate(ZTT (NGPSLEV))
    allocate(ZHU (NGPSLEV))
    allocate(ZUU (NGPSLEV))
    allocate(ZVV (NGPSLEV))
    !
    allocate( H    (GPSRO_MAXPRFSIZE) )
    allocate( AZMV (GPSRO_MAXPRFSIZE) )
    allocate( ZOBS (GPSRO_MAXPRFSIZE) )
    allocate( ZREF (GPSRO_MAXPRFSIZE) )
    allocate( ZOFF (GPSRO_MAXPRFSIZE) )
    allocate( ZERR (GPSRO_MAXPRFSIZE) )
    allocate( RSTV (GPSRO_MAXPRFSIZE) )
    allocate( ZMHX (GPSRO_MAXPRFSIZE) )
    !
    !     Loop over all header indices of the 'RO' family:
    !
    call obs_set_current_header_list(lobsSpaceData,'RO')
    HEADER: do
       INDEX_HEADER = obs_getHeaderIndex(lobsSpaceData)
       if (INDEX_HEADER < 0) exit HEADER
       !
       !     *  Process only refractivity data (codtyp 169)
       !
       IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
       IF ( IDATYP .EQ. 169 ) THEN
          !
          !     *     Scan for requested data values of the profile, and count them
          !
          ASSIM = .FALSE.
          NH = 0
          call obs_set_current_body_list(lobsSpaceData, INDEX_HEADER)
          BODY: do 
             index_body = obs_getBodyIndex(lobsSpaceData)
             if (index_body < 0) exit BODY
             IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 ) THEN
                ASSIM = .TRUE.
                NH = NH + 1
             END IF
          END DO BODY
          !
          !     *     If assimilations are requested, prepare and apply the observation operator
          !
          IF (ASSIM) THEN
             iProfile=gps_iprofile_from_index(INDEX_HEADER)
             !
             !     *        Basic geometric variables of the profile:
             !
             IAZM = obs_headElem_i(lobsSpaceData,OBS_AZA,INDEX_HEADER)
             ISAT = obs_headElem_i(lobsSpaceData,OBS_SAT,INDEX_HEADER)
             Rad  = obs_headElem_r(lobsSpaceData,OBS_TRAD,INDEX_HEADER)
             Geo  = obs_headElem_r(lobsSpaceData,OBS_GEOI,INDEX_HEADER)
             zAzm = 0.01d0*IAZM / MPC_DEGREES_PER_RADIAN_R8
             zMT  = col_getHeight(lcolumnhr,NGPSLEV,INDEX_HEADER,'TH')/RG
             !     
             !     *        Profile at the observation location:
             !
             zLat = obs_headElem_r(lobsSpaceData,OBS_LAT,INDEX_HEADER)
             zLon = obs_headElem_r(lobsSpaceData,OBS_LON,INDEX_HEADER)
             Lat  = zLat * MPC_DEGREES_PER_RADIAN_R8
             Lon  = zLon * MPC_DEGREES_PER_RADIAN_R8
             Azm  = zAzm * MPC_DEGREES_PER_RADIAN_R8
             sLat = sin(zLat)
             zMT  = zMT * RG / gpsgravitysrf(sLat)
             zP0  = col_getElem(lcolumnhr,1,INDEX_HEADER,'P0')
             DO JL = 1, NGPSLEV
                !
                !     *           Profile x
                !
                ZPP(JL) = col_getPressure(LCOLUMNHR,JL,INDEX_HEADER,'TH')
                !     *           True implementation of ZDP (dP/dP0)
                ZDP(JL) = col_getPressureDeriv(LCOLUMNHR,JL,INDEX_HEADER,'TH')
                ZTT(JL) = col_getElem(lcolumnhr,JL,INDEX_HEADER,'TT') - p_TC
                ZHU(JL) = col_getElem(lcolumnhr,JL,INDEX_HEADER,'HU')
                ZUU(JL) = 0.d0
                ZVV(JL) = 0.d0
             END DO

             if((col_getPressure(lcolumnhr,1,index_header,'TH') + 1.0d-4) .lt. &
                  col_getPressure(lcolumnhr,1,index_header,'MM')) then
                ! case with top thermo level above top momentum level (Vcode=5002)
                do jl = 1, nwndlev
                   zuu(jl) = col_getElem(lcolumnhr,jl  ,index_header,'UU') * p_knot
                   zvv(jl) = col_getElem(lcolumnhr,jl  ,index_header,'VV') * p_knot
                enddo
             else
                ! case without top thermo above top momentum level or unstaggered (Vcode=5001/4/5)
                do jl = 1, nwndlev-1
                   zuu(jl) = col_getElem(lcolumnhr,jl+1,index_header,'UU') * p_knot
                   zvv(jl) = col_getElem(lcolumnhr,jl+1,index_header,'VV') * p_knot
                enddo
                zuu(nwndlev) = zuu(nwndlev-1)
                zvv(nwndlev) = zuu(nwndlev-1)
             endif
             zuu(ngpslev) = zuu(nwndlev)
             zvv(ngpslev) = zuu(nwndlev)
             !     
             !     *        GPS profile structure:
             !
             call gps_struct1sw(ngpslev,zLat,zLon,zAzm,zMT,Rad,geo,zP0,zPP,zDP,zTT,zHU,zUU,zVV,prf)
             !
             !     *        Prepare the vector of all the observations:
             !
             NH1 = 0
             !
             !     *        Loop over all body indices for this index_header:
             !     *        (start at the beginning of the list)
             !
             call obs_set_current_body_list(lobsSpaceData, index_header)
             BODY_2: do 
                index_body = obs_getBodyIndex(lobsSpaceData)
                if (index_body < 0) exit BODY_2
                IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 ) THEN
                   NH1      = NH1 + 1
                   H(NH1)   = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
                   AZMV(NH1)= zAzm
                   ZOBS(NH1)= obs_bodyElem_r(lobsSpaceData,OBS_VAR,INDEX_BODY)
                   !     *              Reference value:
                   IF (LEVELGPSRO.EQ.1) THEN
                      ZREF(NH1) = 0.025d0*exp(-(H(NH1)-Rad)/6500.d0)
                   ELSE
                      ZREF(NH1) = 300.d0*exp(-H(NH1)/6500.d0)
                   END IF
                END IF
             END DO BODY_2
             !
             !     *        Apply the observation operator:
             !
             IF (LEVELGPSRO.EQ.1) THEN
                CALL GPS_BNDOPV1(H, AZMV, NH, PRF, RSTV)
             ELSE
                CALL GPS_REFOPV (H,       NH, PRF, RSTV)
             END IF
             !
             !     *        Perform the (H(x)-Y)/R operation:
             !
             DO NH1 = 1, NH
                ZMHX(NH1) = RSTV(NH1)%VAR
                !
                !     *           Normalized offset:
                !
                IF ( NUMGPSSATS .GE. 1 ) THEN
                   ZOFF(NH1) = (ZOBS(NH1) - ZMHX(NH1)) / ZREF(NH1)
                ELSE
                   ZOFF(NH1) = (ZOBS(NH1) - ZMHX(NH1)) / ZMHX(NH1)
                END IF
             END DO
             !
             !     *        The procedure below is well tested to collectively
             !     *        create error profiles from the data profile, and
             !     *        intended to be used for these data.
             !
             DH = 5000.d0
             IF (LEVELGPSRO.EQ.1) THEN
                ZMIN=0.01D0
             ELSE
                ZMIN=0.002D0
             END IF

             IF (LEVELGPSRO.EQ.2) THEN
                DO NH1 = 1, NH
                   SUM0=0.d0
                   SUM1=0.d0
                   DO JH = 1, NH
                      DDH=H(JH)-H(NH1)
                      SUM0=SUM0+EXP(-(DDH/DH)**2)
                      SUM1=SUM1+EXP(-(DDH/DH)**2)*ZOFF(JH)**2
                   END DO
                   ZERR(NH1)=(SUM1/SUM0)**0.5D0
                   IF ( NUMGPSSATS .GE. 1 ) THEN
                      IF (ISAT.EQ.3 .OR. ISAT.EQ.4) ZERR(NH1) = 2*ZERR(NH1)
                   END IF
                   IF ( ZERR(NH1) < ZMIN ) ZERR(NH1) = ZMIN
                END DO
             ELSE
                DO NH1 = 1, NH
                   ZERR(NH1)=0.05d0
                   HNH1=H(NH1)-Rad
                   L1=( HNH1.LE.10000.d0 )
                   L2=( HNH1.GT.10000.d0 .AND. HNH1.LT.30000.d0 )
                   L3=( HNH1.GT.30000.d0 )
                   IF ( L1 ) ZERR(NH1)=0.02d0+0.08d0*(10000.d0-HNH1)/10000.d0
                   IF ( L2 ) ZERR(NH1)=0.02d0
                   IF ( L3 ) ZERR(NH1)=0.02d0+0.13d0*(HNH1-30000.d0)/30000.d0
                   IF (ISAT.EQ.3 .OR. ISAT.EQ.4 .OR. ISAT.EQ.5) ZERR(NH1) = 2*ZERR(NH1)
                   IF ( ZERR(NH1) < ZMIN ) ZERR(NH1) = ZMIN
                END DO
             END IF

             NH1 = 0
             !
             !     *        Loop over all body indices for this index_header:
             !     *        (start at the beginning of the list)
             !
             call obs_set_current_body_list(lobsSpaceData, index_header)
             BODY_4: do 
                index_body = obs_getBodyIndex(lobsSpaceData)
                if (index_body < 0) exit BODY_4
                IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 ) THEN
                   NH1 = NH1 + 1
                   !
                   !     *              Observation error    S
                   !
                   IF ( NUMGPSSATS .GE. 1 ) THEN
                      call obs_bodySet_r(lobsSpaceData,OBS_OER,index_body, ZERR(NH1) * ZREF(NH1))
                   ELSE
                      call obs_bodySet_r(lobsSpaceData,OBS_OER,index_body, ZERR(NH1) * ZMHX(NH1))
                   END IF
                END IF
             END DO BODY_4
          END IF
       END IF
    END DO HEADER

    deallocate( RSTV )
    deallocate( ZERR )
    deallocate( ZOFF )
    deallocate( ZREF )
    deallocate( ZOBS )
    deallocate( AZMV )
    deallocate( H    )
    deallocate( ZMHX )

    deallocate(zVV)
    deallocate(zUU)
    deallocate(zHU)
    deallocate(zTT)
    deallocate(zDP)
    deallocate(zPP)

    WRITE(*,*)'EXIT SETERRGPSRO'
    RETURN
  END SUBROUTINE OER_SETERRGPSRO

  SUBROUTINE oer_SETERRGPSGB(columnhr,lobsSpaceData,ldata,analysisMode)
    !
    !s/r SETERRGPSGB - SET OBSERVATION ERROR COVARIANCES FOR GB-GPS ZTD DATA
    !
    !*    Purpose:
    !             - Set the observation errors [OBS_OER] and Std(O-P) [OBS_HPHT] for GB-GPS ZTD data.
    !               (GPS surface met obs errors are set before in observation_erreurs_mod.ftn90.
    !               The ZTD error is also initialized to the "formal error" or to 1.0 if missing.)
    !             - Returns ldata=.false. if there are no GPS ZTD data to assimilate
    !               and also sets the modgpsztd_mod variable numGPSZTD = 0.
    !
    IMPLICIT NONE
    !
    ! NOTE: YZDERRWGT IN modgpsztd_mod (FROM NML FILE) IS USED FOR ERROR WEIGHTING
    !       OF TIME SERIES (FGAT) GPS ZTD OBSERVATIONS TO ACCOUNT FOR TEMPORAL ERROR
    !       CORRELATIONS.
    !
    type(struct_columnData) :: columnhr
    type(struct_obs)        :: lobsSpaceData
    logical                 :: ldata
    logical                 :: analysisMode

    INTEGER INDEX_BODY, INDEX_HEADER, ITYP, IASS, IZTDJ, NBRPDATE, ICOUNT, ICOUNT2
    integer ielem, nlev_T

    LOGICAL LLCZTDE, LLFER, LLFZTDE, LLZTD, LLRZTDE, ASSIM, ERRSET, DEBUG, LESTP
    LOGICAL LLZWD

    !
    REAL*8  ZTDERR, ZZTD, ZMINZDE, ZPSFC, ZHD, ZWD, ZTDOER, ZLEV, ZVAL, ZZWD
    REAL*8  ZBTSFC, ZBPSFC, ZBZSFC, ZDZ, ZSTDOMP

    !
    !     ZZDERMIN = MIN ZTD OER VALUE (M), ZZFERREJ = MAX FERR VALUE (M) FOR REJECTION
    !     ZZDERMAX = MAX ZTD OER VALUE (M)
    !     ZTDERFAC = MULTIPLICATION FACTOR FOR FORMAL ZTD MEASUREMENT ERROR
    !     ZOPEFAC  = FRACTION OF REGRESSION EQUATION SD(O-P) TO USE AS ZTD OBSERVATION ERROR
    !     ----------------------------------------------------------------------------------
    !
    REAL*8 ZZDERMIN, ZZFERREJ, ZZDERMAX, ZTDERFAC, ZOPEFAC
    DATA ZZDERMIN /0.004D0/
    DATA ZZFERREJ /0.015D0/
    DATA ZZDERMAX /0.030D0/
    DATA ZTDERFAC /3.0D0/
    DATA ZOPEFAC  /1.0D0/
    !
    !     FOR ESTIMATION OF PSFC (IF MISSING)
    !       ZGAMMA = (NEG. OF) TEMPERATURE LAPSE RATE (K/M)
    !
    REAL*8 ZGAMMA
    DATA ZGAMMA /0.0065D0/

    !     ----------------------------------------------------------------------------------
    !     LINEAR REGRESSION EQUATION CONSTANTS AND COEFFS FOR ZTD ERROR AND STD(O-P):

    !     ZRCONST, ZRCOEFF:
    !       - From linear regression of Desroziers error estimates binned by observed ZWD.
    !       - Gives ZTDerror (mm) as function of ZWD (m):
    !            ZTDerror(mm) = ZRCONST + ZRCOEFF*ZWD(m)
    !     ZRCONST2, ZRCOEFF2:
    !       - From linear regression of Std(O-P) binned by observed ZWD.
    !       - Gives Std(O-P) (mm) as function of ZWD (m):
    !            Std(O-P)(mm) = ZRCONST2 + ZRCOEFF2*ZWD(m)
    !     ----------------------------------------------------------------------------------
    REAL*8 ZRCONST, ZRCOEFF, ZRCONST2, ZRCOEFF2
    DATA  ZRCONST  /5.12D0/
    DATA  ZRCOEFF  /26.4D0/
    DATA  ZRCONST2 /6.67D0/
    DATA  ZRCOEFF2 /42.6D0/
    !
    !
    WRITE(*,*) 'ENTER SETERRGPSGB'
    !
    DEBUG = .FALSE.

    LLCZTDE = .FALSE.
    LLRZTDE = .FALSE.
    LLFZTDE = .FALSE.
    IF (YZTDERR .LT. 0.0D0) THEN
       LLFZTDE = .TRUE.
    ELSE IF (YZTDERR .GT. 0.0D0) THEN
       LLCZTDE = .TRUE.
    ELSE
       LLRZTDE = .TRUE.
    END IF

    nlev_T = col_getNumLev(columnhr,'TH')

    ldata = .false.
    ICOUNT  = 0
    ICOUNT2 = 0
    !
    !     Loop over all header indices of the 'GP' family:
    !
    call obs_set_current_header_list(lobsSpaceData,'GP')
    HEADER: DO
       index_header = obs_getHeaderIndex(lobsSpaceData)
       if (index_header < 0) exit HEADER
       NBRPDATE  = obs_headElem_i(lobsSpaceData,OBS_DAT,INDEX_HEADER)
       LLZTD     = .FALSE.
       LLFER     = .FALSE.
       LLZWD     = .FALSE.
       ASSIM     = .FALSE.
       ERRSET    = .FALSE.
       ZZTD      = -1.0D0
       ZZWD      = -1.0D0
       ZPSFC     = -1.0D0
       LESTP     = .FALSE.
       ZSTDOMP   = 15.0D0*0.001D0

       !   Get Psfc (Pa), Tsfc (K) and model surface height (m) from background profile

       ZBPSFC = col_getElem(columnhr,1,index_header,'P0')
       ZBTSFC = col_getElem(columnhr,nlev_T,index_header,'TT')
       ZBZSFC = col_getHeight(columnhr,nlev_T,index_header,'TH')/RG
       !
       !    Loop over all body indices of current report; Set the ZTD error if
       !    constant value specified (LLCZTDE=true). Get GPS height and Psfc obs (if any).
       !    Get ZTD obs, ZTD formal error and ZWD observation.
       !
       call obs_set_current_body_list(lobsSpaceData, INDEX_HEADER)
       BODY: DO 
          index_body = obs_getBodyIndex(lobsSpaceData)
          if (index_body < 0) exit BODY
          ITYP   = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY)
          IASS   = obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY)
          ZVAL   = obs_bodyElem_r(lobsSpaceData,OBS_VAR,INDEX_BODY)
          !         Store Psfc
          IF ( ITYP .EQ. BUFR_NEPS ) THEN
             IF ( ZVAL .GT. 0.0D0 ) ZPSFC = ZVAL
          END IF
          !         Set ZTDOER to constant value (if LLCZTDE); get value of ZTD, 
          !         ZTD formal error (OBS_OER) and antenna height (OBS_PPP).
          IF ( ITYP .EQ. BUFR_NEZD ) THEN
             IF ( LLCZTDE ) THEN
                ZTDOER = YZTDERR
                ERRSET = .TRUE.
             END IF
             ZLEV   = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
             ZTDERR = obs_bodyElem_r(lobsSpaceData,OBS_OER,INDEX_BODY)
             IF ( ZTDERR .NE. 1.0D0 ) LLFER = .TRUE.
             IZTDJ = INDEX_BODY
             IF ( ZVAL .GT. 0.0D0 ) THEN
                ZZTD = ZVAL
                LLZTD = .TRUE.
             END IF
             IF ( IASS .EQ. 1 ) ASSIM = .TRUE.
          END IF
          IF ( ITYP .EQ. BUFR_NEZW ) THEN
             IF ( ZVAL .GT. 0.0D0 ) THEN
                ZZWD = ZVAL
                LLZWD = .TRUE.
             END IF
          END IF
       END DO BODY

       !      Initialize Std(O-P) to 15 mm  for ZTD observation (for bgck mode)
       IF ( LLZTD .AND. .NOT.analysisMode) &
            call obs_bodySet_r(lobsSpaceData,OBS_HPHT,IZTDJ,ZSTDOMP)

       !      Replace formal ZTD error with real error for all ZTD to be assimilated.
       !      Set Std(O-P) as function of ZWD for ZTD observation and store in OBS_HPHT. 

       IF ( ASSIM ) THEN
          IF ( LLZTD ) THEN
             ldata = .true.
             ICOUNT = ICOUNT + 1
             IF ( LLZWD ) THEN
                ZWD = ZZWD
             ELSE
                !               If Psfc obs is missing, estimate the pressure from background state
                IF ( ZPSFC .LT. 0.0D0 ) THEN
                   LESTP = .TRUE.
                   ZDZ = ZLEV - ZBZSFC
                   ZPSFC  = ZBPSFC * &
                        (1.0D0-(ZGAMMA/ZBTSFC)*ZDZ)**(RG/(MPC_RGAS_DRY_AIR_R8*ZGAMMA))
                   ICOUNT2 = ICOUNT2 + 1
                END IF
                !                Compute the hydrostatic delay ZHD (m) from Psfc (Pa)
                ZHD = 2.2766D-05 * ZPSFC
                !               Compute the wet delay (m) from ZTD and ZHD. Avoid negative ZWD.
                IF ( ZHD .GT. ZZTD ) THEN
                   ZWD = 0.0D0
                ELSE
                   ZWD = ZZTD - ZHD
                END IF
             END IF
             !              Std(O-P) for background check. Limit to 30 mm in case ZTD obs is bad (too high).             
             ZSTDOMP = (ZRCONST2 + ZRCOEFF2*ZWD)*0.001D0
             ZSTDOMP = MIN(ZZDERMAX, ZSTDOMP)
             !             Compute ZTD error as a function of ZWD using regression coeff (SD(O-P) vs ZWD).
             !             Take fraction ZOPEFAC of computed error and convert from mm to m.
             !             Ensure error is > ZZDERMIN and < ZZDERMAX
             IF ( .NOT. ERRSET ) THEN 
                ZMINZDE = ZRCONST + ZRCOEFF*ZWD
                ZMINZDE = ZMINZDE * ZOPEFAC * 0.001D0
                IF (LLRZTDE) THEN
                   ZTDOER = MAX(ZZDERMIN, ZMINZDE)
                   ZTDOER = MIN(ZZDERMAX, ZTDOER)
                ELSE
                   IF (LLFER) THEN
                      ZTDOER = MAX(ZZDERMIN, ZTDERR*ZTDERFAC)
                   ELSE
                      ZTDOER = MAX(ZZDERMIN, ZMINZDE)
                      ZTDOER = MIN(ZZDERMAX, ZTDOER)
                   END IF
                END IF
                !  Ensure that error is not less than formal error ZTDERR
                IF (LLFER) THEN
                   IF (DEBUG .AND. ICOUNT .LE. 50) THEN
                      WRITE(*,*) obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER), &
                           ' FORMAL ERR, OBS ERROR (mm) = ', &
                           ZTDERR*1000.D0, ZTDOER*1000.D0
                   END IF
                   ZTDOER = MAX(ZTDOER, ZTDERR)
                END IF
             END IF
             !  *** APPLY TIME-SERIES WEIGHTING FACTOR TO OBSERVATION ERROR (YZDERRWGT=1 FOR 3D THINNING)
             call obs_bodySet_r(lobsSpaceData,OBS_OER,IZTDJ, ZTDOER*YZDERRWGT)
             IF (.NOT.analysisMode) call obs_bodySet_r(lobsSpaceData,OBS_HPHT,IZTDJ, ZSTDOMP)
             IF (DEBUG .AND. (ICOUNT2 .LE. 50) .AND. LESTP) THEN
                WRITE(*,*) 'TAG    SITE    ZTD    ERROR    ELEV    PSFC    ZWD     STDOMP'
                WRITE(*,*) 'ERRDEBUG ', obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER), &
                     ZZTD*1000.D0, ZTDOER*1000.D0, ZLEV, ZPSFC/100.D0, ZWD*1000.D0, ZSTDOMP*1000.D0
             END IF
          ELSE
             CALL utl_abort('SETERRGPSGB: ERROR:NEGATIVE ZTD VALUE!')
          END IF
       END IF

       !
    END DO  HEADER

    !      IF (DEBUG) CALL utl_abort('******DEBUG STOP*******')

    IF (.not.ldata) numGPSZTD = 0

    IF (ldata) WRITE(*,*) ' numGPSZTD = ', ICOUNT

    WRITE(*,*) 'EXIT SETERRGPSGB'

  END SUBROUTINE OER_SETERRGPSGB

!-----------------------------------------------------------------------------------------

end module obsErrors_mod

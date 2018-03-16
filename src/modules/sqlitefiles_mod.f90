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

!--------------------------------------------------------------------------
!! MODULE sqliteFiles (prefix="sqlite")
!!
!! *Purpose*: To store the filenames of the sqlite observation files and call
!!            subroutines in readSqlite to read and update sqlite files.
!!
!--------------------------------------------------------------------------
module sqliteFiles_mod
  
  use mathPhysConstants_mod
  use mpivar_mod
  use sqliteread_mod
  use obsSpaceData_mod
  use fSQLite
  use utilities_mod
  use codePrecision_mod
  use obsUtil_mod

  implicit none
  save
  private

  ! public variables
  public :: sqlite_nfiles, sqlite_cfilnam,sqlite_cfamtyp

  ! public procedures
  public :: sqlite_setupfiles, sqlf_updateFile, sqlf_readFile

! type for SQLIte  file handle
  type(fSQL_DATABASE)            :: db
! type for precompiled SQLite statements
  type(fSQL_STATEMENT)           :: stmt,stmt2
  type(FSQL_STATUS)              :: stat

  integer, parameter :: jpfiles=64
  integer :: sqlite_nfiles
  character(len=2)   :: sqlite_cfamtyp(jpfiles)
  integer, parameter :: jmaxsqlitefilename=1060
  character(len=jmaxsqlitefilename) :: sqlite_cfilnam(jpfiles)
  CHARACTER *9         :: RDB4_SCHEMA,STR_DATE
  CHARACTER *128       :: query
  integer :: datesql ,timesql
  character(len=48)  :: sqliteFileMode
  logical :: sqlite_split_L
  character(len=256) :: sqlite_split_mode

contains

SUBROUTINE sqlite_setupfiles(datestamp, sqliteFileMode_in)
  implicit none
!s/r  setup_sqlitefiles -INITIALZE sqlite FILE NAMES and return datestamp
!
!
  character(len=*), intent(in) :: sqliteFileMode_in
  
  INTEGER IER,INBLKS,nulsqlite,JJ
  CHARACTER(len=20) :: CLVALU(JPFILES)
  CHARACTER(len=2) :: CFAMI(JPFILES)
  CHARACTER(len=4) :: cmyidx, cmyidy
  CHARACTER(len=9) :: cmyid
  CHARACTER(len=jmaxsqlitefilename-(1+20+1+9)) :: sqlite_directory
  CHARACTER(len=jmaxsqlitefilename) :: sqlitein   !! the length should be more than len(sqlite_directory)+1+len(clvalu)+1+len(cmyid)
  integer :: length_sqlite_directory, status, length_sqlite_split
  LOGICAL   isExist_L 

  INTEGER KTIME,KDATE,KDATE_RECV,KTIME_RECV
  INTEGER IHANDL,ILONG,DATESTAMP
  INTEGER ITIME,IFLGS,IDBURP,ILAT,ILON,IDX,IDY
  INTEGER IALT,IDELAY,IDATE,IRS,IRUNN,INBLK,ISUP,IXAUX
  INTEGER INSUP,INXAUX

  INTEGER INRECS

  INTEGER NBRPDATE,NBRPHH,ISTAMPOBS,INEWHH,NEWDATE
  REAL*8 DELHH
  INTEGER       IVALS
  CHARACTER*9   CLSTNID
  
  !
  !- Setup the mode
  !
  sqliteFileMode = sqliteFileMode_in

  ier = 0
  call get_environment_variable('OAVAR_BURP_DIRECTORY',sqlite_directory,length_sqlite_directory,ier,.true.)
  if (ier.eq.-1) then
     write(*,*) 'sqlite_setupfiles: Problem when getting the environment variable OAVAR_BURP_DIRECTORY'
     write(*,*) '                 The length of the variable ''sqlite_directory'' is too short to hold the complete value'
     call utl_abort("sqlite_setupfiles")
  end if
  if (ier.gt.1) then
     write(*,*) 'sqlite_setupfiles: Problem when getting the environment variable OAVAR_BURP_DIRECTORY'
     write(*,*) '                 Rely on defaut which is the current directory'
     sqlite_directory = '.'
  end if
  if (ier.eq.1) then
     write(*,*) 'sqlite_setupfiles: The environment variable OAVAR_BURP_DIRECTORY does not exist'
     write(*,*) '                 Rely on defaut which is the current directory'
     sqlite_directory = '.'
  end if
  if (ier.eq.0) then
     write(*,*) 'sqlite_setupfiles: The environment variable OAVAR_BURP_DIRECTORY does exist'
     write(*,*) '                 Retrieve it'
  end if
 
 !
 !- Determine if the observation files are already split by subdomain
 !
 status = 0
    call get_environment_variable('OAVAR_BURP_SPLIT',sqlite_split_mode,length_sqlite_split,status,.true.)

    if ( status == 1 ) then
      write(*,*) 'sqlite_setupfiles: The environment variable OAVAR_BURP_SPLIT has not been detected!'
      write(*,*) '                 The observation files are NOT split'
      sqlite_split_L = .false.  
      sqlite_split_mode = 'DEFAULT'
      ! At this point the code only supports split files
      write(*,*) 'sqlite_setupfiles: Expecting split sqlite files!'
      call utl_abort('sqlite_setupfiles')
    else
      write(*,*)
      write(*,*) 'sqlite_setupfiles: The environment variable OAVAR_BURP_SPLIT has correctly been detected'
      select case ( trim(sqlite_split_mode) )
      case ('yes')
        write(*,*) 'sqlite_setupfiles: The observation files ARE split. Assuming ROUNDROBIN strategy'
        sqlite_split_L = .true.
        sqlite_split_mode = 'ROUNDROBIN'
      case ('no')
        write(*,*) 'sqlite_setupfiles: The observation files are NOT split'
        sqlite_split_L = .false.
        sqlite_split_mode = 'DEFAULT'
        ! At this point the code only supports split files
        write(*,*) 'sqlite_setupfiles: Expecting split sqlite files!'
        call utl_abort('sqlite_setupfiles')
      case ('ROUNDROBIN')
        write(*,*) 'sqlite_setupfiles: The observation files ARE split using the ROUNDROBIN strategy'
        sqlite_split_L = .true.
      case ('LATLONTILES')
        write(*,*) 'sqlite_setupfiles: The observation files ARE split by LAT-LON tiles'
        write(*,*) 'sqlite_setupfiles: LAT-LON TILES NO LONGER PERMITTED, USE ROUND ROBIN!'
        call utl_abort('sqlite_setupfiles')
      case default
        write(*,*) 'sqlite_setupfiles: Unknown sqlite_split_mode ', trim(sqlite_split_mode)
        call utl_abort('sqlite_setupfiles')
      end select
    end if

  write(cmyidy,'(I4.4)') (mpi_npey-mpi_myidy)
  write(cmyidx,'(I4.4)') (mpi_myidx+1)
  cmyid  = trim(cmyidx)//'_'//trim(cmyidy)
  print *, ' CMID=',cmyid

  CLVALU(1) = 'ua.rdb'  
  CLVALU(2) = 'brpai'  
  CLVALU(3) = 'brpsfc'  
  CLVALU(4) = 'brpssmis'
  CLVALU(5) = 'brpairs'
  CLVALU(6) = 'brpto_amsua'
  CLVALU(7) = 'brpto_amsub'
  CLVALU(8) = 'brpcsr'
  CLVALU(9) = 'brpiasi'
  CLVALU(10) = 'brpatms'
  CLVALU(11) = 'brpcris'
  CLVALU(12) = 'brpsw'  
  CLVALU(13) = 'brpsc'  
  CLVALU(14) = 'brppr'  
  CLVALU(15) = 'brpro' 
  CLVALU(16) = 'brpgp' 
  CLVALU(17) = ' '  

  CFAMI(1)   = 'UA' 
  CFAMI(2)   = 'AI' 
  CFAMI(3)   = 'SF' 
  CFAMI(4)   = 'TO' 
  CFAMI(5)   = 'TO' 
  CFAMI(6)   = 'TO' 
  CFAMI(7)   = 'TO' 
  CFAMI(8)   = 'TO' 
  CFAMI(9)   = 'TO' 
  CFAMI(10)  = 'TO' 
  CFAMI(11)  = 'TO' 
  CFAMI(12)  = 'SW' 
  CFAMI(13)  = 'SC' 
  CFAMI(14)  = 'PR' 
  CFAMI(15)  = 'RO' 
  CFAMI(16)  = 'GP' 
  CFAMI(17)  = ' ' 

  IVALS=8
  KDATE=-9999
  KTIME=-9999
  sqlite_nfiles=0
  DO JJ=1,JPFILES 
    IF(CLVALU(JJ) == '') EXIT
    nulsqlite=0
    sqlitein=trim(sqlite_directory)//'/'//trim(CLVALU(JJ))//'_'//trim(cmyid)
    INQUIRE(FILE=trim(sqlitein),EXIST=isExist_L)
    IF (.NOT. isExist_L )THEN
      sqlitein=trim(sqlite_directory)//'/'//trim(CLVALU(JJ))
      INQUIRE(FILE=trim(sqlitein),EXIST=isExist_L)
    else
      sqlite_nfiles=sqlite_nfiles + 1
      sqlite_cfilnam(sqlite_nfiles)=sqlitein
      sqlite_cfamtyp(sqlite_nfiles)=CFAMI(JJ)
      write(*,*) '  SQLITE FILE NAME=',trim(sqlitein)
      CALL fSQL_open( db, sqlitein ,stat)
       query="select schema from rdb4_schema ;"
       RDB4_SCHEMA= SQL_QUERY_CH(db,trim(query))
       print *, ' RDB4_SCHEMA IS =',RDB4_SCHEMA

       query="select date from resume;"
       STR_DATE=SQL_QUERY_CH(db,trim(query))
       read(STR_DATE,*)  datesql
       kdate=datesql

       query="select time from resume;"
       STR_DATE=SQL_QUERY_CH(db,trim(query))
       read(STR_DATE,*)  timesql
       write(*,*) ' DATE and TIME t in FILE=',datesql,timesql
       ktime=timesql

       CALL fSQL_close( db, stat )
    END IF
  END DO

  WRITE(*,*) ' '
  WRITE(*,*)' NUMBER OF SQLITE FILES IS :',sqlite_nfiles
  WRITE(*,*)'TYPE  NAME '
  WRITE(*,*)'----  ---- '
  DO JJ=1,sqlite_nfiles
    WRITE(*,'(1X,A2,1X,A128)' ) sqlite_cfamtyp(JJ),trim(sqlite_cfilnam(JJ))
  END DO

  ! Make sure all mpi tasks have a valid date (important for split sqlite files)
  call rpn_comm_allreduce(kdate,kdate_recv,1,"MPI_INTEGER","MPI_MAX","GRID",ier)
  call rpn_comm_allreduce(ktime,ktime_recv,1,"MPI_INTEGER","MPI_MAX","GRID",ier)


  kdate = kdate_recv
  ktime = ktime_recv

  ier = newdate(istampobs,kdate,ktime,3)
  delhh = 3.0d0
  call INCDATR (datestamp, istampobs, delhh)
  ier = newdate(datestamp,nbrpdate,inewhh,-3)
  nbrphh=KTIME/100
  if (nbrphh .ge. 21 .or. nbrphh .lt. 3) then
    nbrphh = 0
  elseif(nbrphh .ge. 3 .and. nbrphh .lt. 9) then
    nbrphh = 6
  elseif(nbrphh .ge. 9 .and. nbrphh .lt. 15) then
    nbrphh = 12
  else
    nbrphh = 18
  endif
  ier = newdate(datestamp,nbrpdate,nbrphh*1000000,3)
  WRITE(*, *)' SQLITE FILES VALID DATE (YYYYMMDD) : ',nbrpdate
  WRITE(*, *)' SQLITE FILES VALID TIME       (HH) : ',nbrphh

END SUBROUTINE sqlite_setupfiles

subroutine sqlf_readFile(obsdat,fileName,familyType,fileIndex)
!
!
!***********************************************************************
!
!      PURPOSE: READ CMC SQLITE FILES FILL UP OBSSPACEDATA FILE
!
!
!    ARGUMENTS:
!                 INPUT: OBSDAT (CMA_TABLE OBJECTS)
!                OUTPUT: NONE
!
!
!       AUTHOR: P. KOCLAS(CMC/CMDA ) December 2012
!       Revision: Sergey Skachko (ARMA) March 2018
!     NOTE:
!     -SQLITE FILES ARE ASSUMED TO BE PRESENT IN CURRENT WORKING DIRECTORY
!     -IT IS ASSUMED NFILES_SQL as been initialized via NAMELIST (NAMSQLF)
!       IN sqlite_read module.
!
!***********************************************************************
!
      implicit none
      ! arguments
      type (struct_obs), intent(inout) :: obsdat
      character(len=*) :: fileName
      character(len=*) :: familyType
      integer :: fileINdex

      ! locals
      integer :: ibeg, iend, nstn1, nstn2
      integer :: jo
      logical :: obs_full

! SSN     INTEGER                                  :: J,JO
!      CHARACTER *2                             :: FAMILYTYPE
!      CHARACTER *128                           :: FILENAME
!      INTEGER                                  :: IBEG,IEND,START,NSTN1,NSTN2
!      logical                                  :: obs_full
!      INTEGER                                  :: nfiles, nulout,kk

      write(*,*)' '
      write(*,*)'                sqlf_readFile: Starting          '
      write(*,*)' '
!
!SSN      WRITE(*,*) ' sqlite_cfilnam sqlite_nfiles =',  trim(sqlite_cfilnam(1)), sqlite_nfiles
!
!      start=1
!      nfiles=1
!      IBEG=obs_numbody(obsdat)
!      DO J =1,sqlite_nfiles
!
!   GET SCHEMA OF DATABASE FILE FROM resume table
!  --------------------------------
!         FILENAME  = trim(CFILNAM_SQL(J))
!          FILENAME=trim(sqlite_cfilnam(j))
!          FAMILYTYPE=sqlite_cfamtyp(j)
!          WRITE(*,*)' read database -->', j,trim(FILENAME),trim(FAMILYTYPE)
!======================================================================================
!  subroutine obs_status(obsdat, obs_full, numstns_out, numobs_out, kulout)
!
!          IBEG=obs_numbody(obsdat)
          ibeg = obs_numbody(obsdat) + 1
          Nstn1=obs_numheader(obsdat)
         
!	  SELECT CASE(FAMILYTYPE)
         SELECT CASE(trim(familyType)) 
!        ---------------
	    CASE('UA','AI','SW','PR','SF','SC','RO')
!        ---------------
!SSN              call SQL2OBS_CONV(obsdat,familytype, FILENAME, 6 )
              call SQL2OBS_CONV(obsdat,trim(familyType),fileName,fileIndex)
              Nstn2=obs_numheader(obsdat)
              iend=obs_numbody(obsdat)
              call FDTOUV_OBSDAT(obsdat,nstn1+1,nstn2,MPC_missingValue_R4)
              call ADJUST_HUM_GZ  (obsdat,nstn1+1,nstn2)
              call ADJUST_SFVCOORD(obsdat,nstn1+1,nstn2)             
!        ---------------
!        ---------------
	    CASE('TO')
!        ---------------
!SSN              call SQL2OBS_TOVS(obsdat,familytype, FILENAME, 6)
              call SQL2OBS_TOVS(obsdat,trim(familyType),fileName,fileIndex)
              Nstn2=obs_numheader(obsdat)
              iend=obs_numbody(obsdat)

	 END SELECT

        DO jo=nstn1+1,nstn2
!SSN	    call obs_headSet_i(obsdat,OBS_OTP,JO,J)
!            call obs_headSet_i(obsdat,OBS_IDF,JO,J)
!	    call obs_setFamily(obsdat,trim(FAMILYTYPE),JO)
	    call obs_headSet_i(obsdat,OBS_OTP,JO,fileIndex)
            call obs_headSet_i(obsdat,OBS_IDF,JO,fileIndex)
	    call obs_setFamily(obsdat,trim(familyType),JO)

	 END DO
!!
!!=======================================================================================
!!      WRITE(*,*)' readsql obs_numheader(obsdat)',FAMILYTYPE,j,obs_numheader(obsdat)
!!      WRITE(*,*)' readsql obs_numbody(obsdat)',  FAMILYTYPE,j,obs_numbody  (obsdat)
!!=========================================================================================

      write(*,*) ' readsql IBEG,IEND =',ibeg,iend
      do jo=ibeg,iend
        if ( obs_columnActive_RB(obsdat,OBS_OMA) )  call obs_bodySet_r(obsdat,OBS_OMA ,JO,MPC_missingValue_R8)
        if ( obs_columnActive_RB(obsdat,OBS_OMP) )  call obs_bodySet_r(obsdat,OBS_OMP ,JO,MPC_missingValue_R8)
        if ( obs_columnActive_RB(obsdat,OBS_OER) )  call obs_bodySet_r(obsdat,OBS_OER ,JO,MPC_missingValue_R8)
        if ( obs_columnActive_RB(obsdat,OBS_HPHT) ) call obs_bodySet_r(obsdat,OBS_HPHT,JO,MPC_missingValue_R8)
        if ( obs_columnActive_RB(obsdat,OBS_WORK) ) call obs_bodySet_r(obsdat,OBS_WORK,JO,MPC_missingValue_R8)
      enddo

!SSN      END DO
    
    ! For GP family, initialize OBS_OER to element 15032 (ZTD formal error) 
    ! for all ZTD data (element 15031)
    if ( trim(familyType) == 'GP') then
      write(*,*)' Initializing OBS_OER for GB-GPS ZTD to formal error (ele 15032)'
      call set_err_gbgps(obsdat,Nstn1+1,Nstn2)
    end if

    write(*,*) 'sqlf_readFile: obs_numheader(obsdat)', trim(familyType), obs_numheader(obsdat)
    write(*,*) 'sqlf_readFile: obs_numbody(obsdat)  ', trim(familyType), obs_numbody  (obsdat)

       write(*,*)' '
!      WRITE(*,*)'================================================='
       write(*,*)'      sqlf_readFile     END                   '
!      WRITE(*,*)'================================================='
       write(*,*)' '
!      RETURN

end subroutine  sqlf_readFile
!================================

SUBROUTINE sqlf_updateFile(obsSpaceData,fileName,familyType,fileIndex)
!
!     PURPOSE: READ OBSDAT AND UPDATE CMC SQLITE FILES
!
!     ARGUMENTS:
!                   obsSpaceData   - obsdat-file object
!
!       AUTHOR: P. KOCLAS(CMC CMDA)
!       Revision:
!     NOTE:
!     SQLITE FILES ARE ASSUMED TO BE PRESENT IN CURRENT WORKING DIRECTORY
!

      implicit none
      ! arguments
      type (struct_obs), intent(inout) :: obsSpaceData
      character(len=*) :: fileName
      character(len=*) :: familyType
      integer :: fileIndex

      ! locals
      integer :: headerIndex
      !type (struct_obs), intent(inout) :: obsSpaceData
      type(fSQL_DATABASE)                      :: db
      type(FSQL_STATUS)                        :: stat
      !INTEGER*4                                :: J,JO

      !character(len=256) :: FILENAME
      !integer length_burp_split
      !integer   ::  EXDB,EXFIN

      call tmg_start(97,'POST_UPDATESQL')

      write(*,*) 'sqlf_updateFile: Starting'

!      if (trim(sqliteFileMode) == 'analysis') call vint3dfd(obs_oma,obsSpaceData)
!      call vint3dfd(obs_omp,obsSpaceData)
!      if (trim(sqliteFileMode) == 'analysis' .or. trim(sqliteFileMode) == 'FSO') call setassflg(obsSpaceData)
!      call flaguvtofd_obsdat(obsSpaceData)
!
!  ------NOTE----------
! currently supported families of data 'UA' 'AI' 'SC' 'SF' 'SW' 'TO'
!
!     READ DATA FROM FILES CONTAINED IN ARRAY CLVAL.
!
!      WRITE(*,*)' '
!      WRITE(*,*)'================================================='
!      WRITE(*,*)'                sqlf_updateFile BEGIN                '
!      WRITE(*,*)'================================================='
!      WRITE(*,*)' '


!! SSN ????     ! redistribute obs data to how it was just after reading the files
!     call obs_MpiRedistribute(obsSpaceData,OBS_IPF)

!     ISTAMP=EXDB('UPDATE_SQL','DEBUT','NON')
      call SQL2OBS_NML('NAMSQLUPDATE')
!SSN      DO J =1,sqlite_nfiles
!          FILENAME=trim(sqlite_cfilnam(j))

!          WRITE(*,*)'  open database -->', trim(FILENAME)
!         ===============================================
          CALL fSQL_open( db, fileName ,stat)
          if ( fSQL_error(stat) /= FSQL_OK ) then
	    write(*,*) 'fSQL_open: ', fSQL_errmsg(stat)
	  endif

          !call    updsql(db,obsSpaceData,j)
          !call insertsql(db,obsSpaceData,j)
          call    updsql(db,obsSpaceData,familyType,fileName,fileIndex)
          call insertsql(db,obsSpaceData,familyType,fileName,fileIndex)
          
          CALL fSQL_close( db, stat )
!         ===============================================
          WRITE(*,*)'  closed database -->', trim(FILENAME)
        
!!!SSN      END DO
!
      WRITE(*,*)' '
      WRITE(*,*)'================================================='
      WRITE(*,*)'                sqlf_updateFile    END               '
      WRITE(*,*)'================================================='
      WRITE(*,*)' '

      call tmg_stop(97)

end subroutine sqlf_updateFile

end module sqliteFiles_mod

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
  public :: sqlf_getDateStamp, sqlf_updateFile, sqlf_readFile

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

subroutine sqlf_getDateStamp(datestamp, sqliteFileName)

  ! Author     : Pierre Koclas 
  ! Adaptation from sqlite_setupFiles 
  !            : Sergey Skachko (ARMA), March 2018

  implicit none

  ! arguments
  integer :: dateStamp
  character(len=*), intent(in) :: sqliteFileName

  ! locals 
  logical :: isExist_L 
  integer :: ier,ktime,kdate,ivals,kdate_recv,ktime_recv
  integer :: inrecs, mrfopc
  real*8  :: delhh
  integer :: nbrpdate,nbrphh,istampobs,inewhh,newdate
  
  !
  !- Get the date from the sqlite files
  !
  
  ier = mrfopc('MSGLVL','FATAL')

  ivals=8
  kdate=-9999
  ktime=-9999
    
  inquire(file=trim(sqliteFileName),exist=isExist_L)
  if ( isExist_L ) then
    
    write(*,*)' Open File : ',trim(sqliteFileName)
    call fSQL_open( db, trim(sqliteFileName) ,stat)
    query="select date from resume;"
    STR_DATE=SQL_QUERY_CH(db,trim(query))
    read(STR_DATE,*)  datesql
    kdate=datesql

    query="select time from resume;"
    STR_DATE=SQL_QUERY_CH(db,trim(query))
    read(STR_DATE,*)  timesql
    write(*,*) ' DATE and TIME t in FILE=',datesql,timesql
    ktime=timesql
    call fSQL_close( db, stat )

  endif

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
  write(*,*)' SQLITE FILES VALID DATE (YYYYMMDD) : ', nbrpdate
  write(*,*)' SQLITE FILES VALID TIME       (HH) : ', nbrphh
  write(*,*)' SQLITE FILES DATESTAMP             : ', datestamp

end subroutine sqlf_getDateStamp

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
  character(len=*)                 :: fileName
  character(len=*)                 :: familyType
  integer                          :: fileINdex

  ! locals
  integer :: ibeg, iend, nstn1, nstn2
  integer :: jo
  logical :: obs_full
  character(len=*), parameter :: my_name = 'sqlf_readFile'
  character(len=*), parameter :: my_warning = '****** '// my_name //' WARNING: '
  character(len=*), parameter :: my_error   = '******** '// my_name //' ERROR: '
  real(obs_real)  :: misg
  
  ! debug:******************************
  logical, parameter :: ldebug =.false.
  integer :: bodyIndex,headerIndex,varno,iqiv,igav,ilansea,azimuth,inst,ifov,clf,saz,idsat,roqc,id_obs
  character*9 :: stid_l
  real    :: lat,lon,alt,var,channel,ealoc,geoun,sigmao,omp,oma
  !*************************************

  write(*,*)' '
  write(*,*)'                '//trim(my_name)//': Starting          '
  write(*,*)' '
  misg = real(MPC_missingValue_R8,OBS_REAL)
  write(*,*)my_name//': FileName   : ', trim(FileName), mpi_myid
  write(*,*)my_name//': FamilyType : ', FamilyType, mpi_myid

  ibeg = obs_numbody(obsdat) + 1
  Nstn1=obs_numheader(obsdat)
  call sqlr_readSqlite(obsdat,trim(familyType),trim(fileName),fileIndex)
  Nstn2=obs_numheader(obsdat)
  iend=obs_numbody(obsdat)

  if(ldebug) then
    ! debug ###################################################################
    write(*,*)'SSN: debug!!! iend, numproc ', iend,mpi_myid
    call filt_sethind_util(obsdat)
    do bodyIndex = 1, iend
      headerIndex = obs_bodyElem_i(obsdat, OBS_HIND,   bodyIndex)
      call obs_prnthdr(obsdat,headerIndex,100+mpi_myid)
      call obs_prntbdy(obsdat,headerIndex,100+mpi_myid)
      !if ( obs_headElem_i(obsdat, OBS_ONM, headerIndex) /= 0 ) cycle
      stid_l      = obs_Elem_c(obsdat, 'STID',  headerIndex)
      lat         = obs_headElem_r(obsdat, OBS_LAT,  headerIndex)
      lon         = obs_headElem_r(obsdat, OBS_LON,  headerIndex)
      alt         = obs_headElem_r(obsdat, OBS_ALT,  headerIndex)
      channel     = obs_bodyElem_r(obsdat, OBS_PPP,    bodyIndex)
      var         = obs_bodyElem_r(obsdat, OBS_VAR,    bodyIndex)
      varno       = obs_bodyElem_i(obsdat, OBS_VNM,    bodyIndex)
      iqiv        = obs_headElem_i(obsdat, OBS_SZA,  headerIndex)
      igav        = obs_headElem_i(obsdat, OBS_SUN,  headerIndex)
      ilansea     = obs_headElem_i(obsdat, OBS_OFL,  headerIndex)
      azimuth     = obs_headElem_i(obsdat, OBS_AZA,  headerIndex)
      inst        = obs_headElem_i(obsdat, OBS_INS,  headerIndex)
      ifov        = obs_headElem_i(obsdat, OBS_FOV,  headerIndex)
      clf         = obs_headElem_i(obsdat, OBS_CLF,  headerIndex)
      saz         = obs_headElem_i(obsdat, OBS_SAZ,  headerIndex)
      idsat       = obs_headElem_i(obsdat, OBS_SAT,  headerIndex) 
      roqc        = obs_headElem_i(obsdat, OBS_ROQF, headerIndex) 
      geoun       = obs_headElem_r(obsdat, OBS_GEOI, headerIndex)
      ealoc       = obs_headElem_r(obsdat, OBS_TRAD, headerIndex)
      sigmao      = obs_bodyElem_r(obsdat, OBS_OER,    bodyIndex)
      omp         = obs_bodyElem_r(obsdat, OBS_OMP,    bodyIndex)
      oma         = obs_bodyElem_r(obsdat, OBS_OMA,    bodyIndex)
      id_obs      = obs_headElem_i(obsdat, OBS_IDO,  headerIndex)
      !write(100+mpi_myid,'(a6,5f12.4,10i8)') trim(stid_l),lon,lat,alt,channel,var,varno,iqiv,igav,ilansea,azimuth,inst,clf,saz,idsat,roqc
       !write(*,'(a6,5f12.4,11i8)') trim(stid_l),lon,lat,alt,channel,var,varno,iqiv,igav,ilansea,azimuth,inst,clf,saz,idsat,roqc,headerIndex
      !!!!write(100+mpi_myid,'(a12,8f12.4,3i8)') trim(stid_l),lon,lat,alt,channel,var,sigmao,omp,oma,varno,id_obs,headerIndex        
      !write(*,'(a12,7f12.4,3i8)') trim(stid_l),lon,lat,alt,channel,var,omp,oma,varno,id_obs,headerIndex
      !write(100+mpi_myid,'(a6,5f12.4,11i8,2f12.4)') trim(stid_l),lon,lat,alt,channel,var,varno,iqiv,igav,ilansea,azimuth,inst,ifov,clf,saz,idsat,roqc,ealoc,geoun
      !write(100+mpi_myid,'(a6,7f12.4,3i8)') trim(stid_l),lon,lat,alt,channel,var,ealoc,geoun,azimuth,idsat,roqc
    enddo
    write(*,*)'SSN: DONE', mpi_myid
    ! ### FIN debug ##########################################################
  endif  
            

  if ( trim(familyType) /= 'TO' ) then
    call FDTOUV_OBSDAT(obsdat,nstn1+1,nstn2,MPC_missingValue_R4)
    call ADJUST_HUM_GZ  (obsdat,nstn1+1,nstn2)
    call ADJUST_SFVCOORD(obsdat,nstn1+1,nstn2)
  endif             
  do jo=nstn1+1,nstn2
    call obs_headSet_i(obsdat,OBS_OTP,jo,fileIndex)
    call obs_headSet_i(obsdat,OBS_IDF,jo,fileIndex)
    call obs_setFamily(obsdat,trim(familyType),jo)
  enddo

  write(*,*) my_name//': ibeg,iend =',ibeg,iend
  do jo=ibeg,iend
    if ( obs_columnActive_RB(obsdat,OBS_OMA) )  call obs_bodySet_r(obsdat,OBS_OMA ,jo,MISG)
    if ( obs_columnActive_RB(obsdat,OBS_OMP) )  call obs_bodySet_r(obsdat,OBS_OMP ,jo,MISG)
    if ( obs_columnActive_RB(obsdat,OBS_OER) )  call obs_bodySet_r(obsdat,OBS_OER ,jo,MISG)
    if ( obs_columnActive_RB(obsdat,OBS_HPHT) ) call obs_bodySet_r(obsdat,OBS_HPHT,jo,MISG)
    if ( obs_columnActive_RB(obsdat,OBS_WORK) ) call obs_bodySet_r(obsdat,OBS_WORK,jo,MISG)
  enddo
    
  ! For GP family, initialize OBS_OER to element 15032 (ZTD formal error) 
  ! for all ZTD data (element 15031)
  if ( trim(familyType) == 'GP') then
    write(*,*)' Initializing OBS_OER for GB-GPS ZTD to formal error (ele 15032)'
    call set_err_gbgps(obsdat,Nstn1+1,Nstn2)
  end if

  write(*,*) my_name//': obs_numheader(obsdat)', trim(familyType), obs_numheader(obsdat)
  write(*,*) my_name//': obs_numbody(obsdat)  ', trim(familyType), obs_numbody  (obsdat)

  write(*,*)' '
  write(*,*)'      '//trim(my_name)//'     END                   '
  write(*,*)' '

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
!       Revision: S. Skachko ARMA, April 2018
!     NOTE:
!     SQLITE FILES ARE ASSUMED TO BE PRESENT IN CURRENT WORKING DIRECTORY
!

      implicit none
      ! arguments
      type (struct_obs), intent(inout) :: obsSpaceData
      character(len=*)                 :: fileName
      character(len=*)                 :: familyType
      integer                          :: fileIndex

      ! locals
      integer :: headerIndex
      type(fSQL_DATABASE)                      :: db
      type(FSQL_STATUS)                        :: stat
     
      character(len=*), parameter :: my_name = 'sqlf_updateFile'
      character(len=*), parameter :: my_warning = '****** '// my_name //' WARNING: '
      character(len=*), parameter :: my_error   = '******** '// my_name //' ERROR: '

      call tmg_start(97,'POST_UPDATESQL')

      write(*,*) my_name//' Starting'
      write(*,*) my_name//': FileName   : ',trim(fileName)
      write(*,*) my_name//': FamilyType : ',FamilyType

      call SQL2OBS_NML('NAMSQLUPDATE')
      call SQL2OBS_NML('NAMSQLinsert')

      call fSQL_open( db, fileName,stat)
      if ( fSQL_error(stat) /= FSQL_OK ) then
        write(*,*) 'fSQL_open: ', fSQL_errmsg(stat)
        write(*,*) my_error, fSQL_errmsg(stat)
      endif
  
      if (trim(familyType) /='TO' ) then
        call    updsql(db,obsSpaceData,familyType,fileName,fileIndex)
        call insertsql(db,obsSpaceData,familyType,fileName,fileIndex)
      endif    
      write(*,*)'  closed database -->', trim(FileName)
      call fSQL_close( db, stat )
        
      write(*,*)' '
      write(*,*)'================================================='
      write(*,*)'                '//trim(my_name)//'    END               '
      write(*,*)'================================================='
      write(*,*)' '

      call tmg_stop(97)

end subroutine sqlf_updateFile

end module sqliteFiles_mod

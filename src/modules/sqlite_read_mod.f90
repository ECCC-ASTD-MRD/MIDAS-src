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
!! MODULE sqliteRead (prefix= no standard prefix)
!!
!! *Purpose*: To read and update SQLITE observation files. Data is stored in 
!!            obsSpaceData object.
!!
!--------------------------------------------------------------------------
module sqliteread_mod

use codePrecision_mod
use ObsSpaceData_mod
use mpi_mod
use fSQLite
use MathPhysConstants_mod
use mpivar_mod

implicit none
!------------------------------------------------------------------------------------------
private            :: cvt_burp_instrum,handle_error,GEN_HEADER,GEN_INFO,GEN_DATA
public             :: SQL2OBS_TOVS, SQL2OBS_CONV, SQL_QUERY_CH,INSERTSQL,SQL2OBS_NML
!------------------------------------------------------------------------------------------

INTEGER*4          :: NELE_INS,LISTELE_INS(15)
INTEGER*4          :: N_ITEMS
CHARACTER *3       :: ITEMLIST(15)
LOGICAL            :: LABORTFULL
INTEGER, PARAMETER :: N_MAX_FILES=16

INTEGER*4           :: NFILES_SQL

CHARACTER (len=2)   :: CFAMTYP_SQL(N_MAX_FILES)
CHARACTER (len=256) :: CFILNAM_SQL(N_MAX_FILES)

CONTAINS
   subroutine GEN_HEADER(obsdat,STNID,LAT,LON,ELEV,STATUS,DATE,TIME,NOBS)
      !
      !------------------------------------------------------------------
      !s/r GEN_HEADER - Initialze header values for an observation object
      !
      ! Author  : P. Koclas, CMC/CMDA September  2012
      !
      !------------------------------------------------------------------
      !

   type (struct_obs), intent(inout) :: obsdat

   INTEGER*4      :: STATUS,DATE,TIME
   INTEGER*4      :: NOBS
   REAL(OBS_REAL) :: LAT,LON,ELEV

   CHARACTER (len = *)  :: STNID

!=========================================================
   call obs_headSet_r(obsdat,OBS_LAT,nobs,LAT)
   call obs_headSet_r(obsdat,OBS_LON,nobs,LON)
   call obs_headSet_r(obsdat,OBS_ALT,nobs,ELEV)
   call obs_set_c(obsdat,'STID',nobs,trim(STNID))
   call obs_headSet_i(obsdat,OBS_ST1,nobs,STATUS)
   call obs_headSet_i(obsdat,OBS_DAT,nobs,DATE)
   call obs_headSet_i(obsdat,OBS_ETM,nobs,TIME)
!=========================================================

!  ==============================
   END subroutine GEN_HEADER
!  ==============================

   subroutine GEN_INFO(obsdat,ID_SAT,INSTRUMENT,ZENITH,CLOUD_COVER,AZIMUTH,SOLAR_ZENITH,SOLAR_AZIMUTH,LAND_SEA,NOBS)
      !
      !----------------------------------------------------------------------------------------------
      !s/r GEN_INFO   - Initialze info values for an observation object
      !
      ! Author  : P. Koclas, CMC/CMDA September  2012
      !
      !----------------------------------------------------------------------------------------------
      !

   type (struct_obs), intent(inout) :: obsdat
   INTEGER*4    :: NOBS

   INTEGER*4    :: ID_SAT,INSTRUMENT,ZENITH,CLOUD_COVER,AZIMUTH,SOLAR_ZENITH,LAND_SEA,SOLAR_AZIMUTH

!=========================================================
   call obs_headSet_i(obsdat,OBS_INS,nobs,INSTRUMENT )
   call obs_headSet_i(obsdat,OBS_SZA,nobs,ZENITH )
   call obs_headSet_i(obsdat,OBS_SAT,nobs,ID_SAT)
   call obs_headSet_i(obsdat,OBS_CLF,nobs,CLOUD_COVER)
   call obs_headSet_i(obsdat,OBS_AZA,nobs,AZIMUTH)
   call obs_headSet_i(obsdat,OBS_SUN,nobs,SOLAR_ZENITH)
   call obs_headSet_i(obsdat,OBS_OFL,nobs,LAND_SEA)
   call obs_headSet_i(obsdat,OBS_SAZ,nobs,SOLAR_AZIMUTH)

!=========================================================

!  ==============================
   END subroutine GEN_INFO
!  ==============================

   subroutine GEN_DATA(obsdat,VCOORD,OBSVALUE,VARNO,FLAG,VCORDTYP,NOBS,NDATA)
      !
      !----------------------------------------------------------------------------------------------
      !s/r GEN_INFO   - Initialze data values for an observation object
      !
      ! Author  : P. Koclas, CMC/CMDA September  2012
      !
      !----------------------------------------------------------------------------------------------
      !
   type (struct_obs), intent(inout) :: obsdat

!pikpik  REAL*4       :: VCOORD,OBSVALUE
    REAL*8       :: VCOORD,OBSVALUE

   INTEGER*4    :: VARNO,FLAG,VCORDTYP
   INTEGER*4    :: NOBS,NDATA

!=========================================================
   call obs_bodySet_r(obsdat,OBS_PPP,NDATA,VCOORD)
   call obs_bodySet_r(obsdat,OBS_VAR,NDATA,OBSVALUE)
   call obs_bodySet_i(obsdat,OBS_VNM,NDATA,VARNO)
   call obs_bodySet_i(obsdat,OBS_FLG,NDATA,FLAG)
   call obs_bodySet_i(obsdat,OBS_VCO,NDATA,VCORDTYP)
!=========================================================

!===================================
   END subroutine GEN_DATA
!===================================

   subroutine SQL2OBS_NML(NML_SECTION)
      !
      !---------------------------------------------
      !s/r SQL2OBS_NML  - READ NAMELIST FILE SECTION
      !
      !---------------------------------------------
      !
   INTEGER*4          :: NULNAM
   CHARACTER *256     :: NAMFILE
   CHARACTER(len = *) :: NML_SECTION
   !-------------------------------------------
   NAMELIST /NAMDIMO/     LABORTFULL
   NAMELIST /NAMSQLinsert/NELE_INS,LISTELE_INS
   NAMELIST /NAMSQLF/     NFILES_SQL,CFILNAM_SQL,CFAMTYP_SQL
   NAMELIST /NAMSQLUPDATE/N_ITEMS,ITEMLIST

   NAMFILE=trim("flnml")
   WRITE(*,*) ' READ NML_SECTION =',trim(NML_SECTION)
!===============================================================================================
   SELECT CASE(trim(NML_SECTION))
     CASE( 'insert')
       NULNAM=4014
       open(NULNAM,FILE=NAMFILE)
       READ(NULNAM,NML=NAMSQLinsert)
     CASE( 'NAMDIMO')
       NULNAM=4015
       open(NULNAM,FILE=NAMFILE)
       READ(NULNAM,NML=NAMDIMO)
     CASE( 'NAMSQLF')
       NULNAM=4016
       open(NULNAM,FILE=NAMFILE)
       READ(NULNAM,NML=NAMSQLF)
     CASE( 'NAMSQLUPDATE')
       NULNAM=4017
       open(NULNAM,FILE=NAMFILE)
       READ(NULNAM,NML=NAMSQLUPDATE)
     CASE DEFAULT
   WRITE(*,*) ' WRONG NAMELIST SECTION: ',trim(NML_SECTION), ' EXIT ROUTINE SQL2OBS_NML'
   RETURN
   END SELECT
!===============================================================================================
   close(NULNAM)
   RETURN

!  ==============================
   END SUBROUTINE SQL2OBS_NML
!  ==============================

   subroutine SQL2OBS_TOVS(obsdat,familytype,FILENAME,KULOUT)
      !
      !----------------------------------------------------------------
      !
      !   Purpose : QUERY DATA FROM SQLITE FILES
      !             WHICH CONTAIN TOVS  DATA
      !             THEN INSERT THE DATA INTO  ObsSpaceData structures
      !
      ! Author  : P. Koclas, CMC/CMDA September  2012
      ! Adaptation : S. Skachko, ARMA March 2018
      !
      !    ARGUMENTS:
      !                 INPUT:
      !
      !                       -OBSDAT  : ObsSpaceData Structure 
      !                       -FILENAME: SQLITE FILE NAME
      !                       -KULOUT  : STDOUT   FILE UNIT NUMBER 
      !----------------------------------------------------------------
      !
   implicit none

   type (struct_obs), intent(inout) :: obsdat
   character*2,intent(in)           :: familytype
   CHARACTER *128,intent(in)        :: FILENAME
   INTEGER,INTENT(in)               :: KULOUT

   INTEGER              :: NULNAM
   CHARACTER *9         :: RDB4_SCHEMA

   INTEGER              :: ID_OBS,ID_DATA,CODTYP,DATE,TIME,STATUS,FLAG,VARNO
   CHARACTER *9         :: STNID

   REAL                 :: LAT,LON
   REAL(OBS_REAL)       :: RELEV,XLAT,XLON
   REAL                 :: OBSVALUEmat
   REAL*8               :: OBSVALUE
   REAL                 :: VCOORDmat
   REAL*8               :: VCOORD
   INTEGER              :: VCORDTYP

   INTEGER              :: IZENITH
   REAL                 :: ZENITH,SOLAR_ZENITH,AZIMUTH,CLOUD_COVER,SURF_EMISS,SOLAR_AZIMUTH
   REAL*8               :: SURF_EMISS8
   INTEGER              :: ID_SAT,LAND_SEA,TERRAIN_TYPE,INSTRUMENT,SENSOR

   INTEGER              :: I,J,COUNT,NLV,NOBS,NOBS_START
   LOGICAL              :: FINISHED
   REAL*8               :: ZEMFACT,ELEV
   integer              :: N_CHUNKS,chunksize,nsplit
   character*12         :: chmodulo,chCPU,CHPART,chCPUM1

   CHARACTER*10         :: CHTIME
   CHARACTER*512        :: QUERY,QUERY_DAT,QUERY_HDR
   CHARACTER*256        :: CFG,CSQLCRIT
   CHARACTER*256        :: COLUMNS_HDR,COLUMNS_DAT
   character(len=256)   :: splitidobs

   INTEGER              :: BITSFLAGON,BITSFLAGOFF
   INTEGER              :: LN
   INTEGER              :: numstns_out, numobs_out
   LOGICAL              :: obsdat_full

   REAL, ALLOCATABLE    :: matdata  (:,:)
   INTEGER              :: NROWS,NCOLUMNS,last_id

!======F90SQLITE TYPES========================
! type for SQLIte  file handle
   type(fSQL_DATABASE)            :: db

! type for precompiled SQLite statements
   type(fSQL_STATEMENT)           :: stmt,stmt2

!type for error status
   type(fSQL_STATUS)              :: stat,stat2

   CHARACTER *256  :: LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT
   INTEGER         :: NBITSOFF,NBITSON,BITOFF(15),BITON(15)
   CHARACTER *256  :: NAMFILE
!=============================================================================================================
   NAMELIST /NAMSQLamsua/LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLamsub/LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLairs/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLiasi/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLcris/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLssmi/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLgo/   LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLcsr/  LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLatms/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
!=============================================================================================================
   nsplit=mpi_nprocs

   write(*,*) ' SQL2OBS_TOVS  mpi_myid nsplit =',mpi_myid,nsplit
   WRITE(chCPU,'(a1,i3)' )'=', mpi_myid
   WRITE(CHPART,'(a7,i3)' )trim('637776/'), mpi_nprocs
   WRITE(CHPART,'(a8,i3)' )trim('2049887/'), mpi_nprocs
   WRITE(chCPUM1,'(i3)' ) mpi_nprocs-1
   WRITE(chMODULO, '(i4)' ) nsplit
   !pikpik splitidobs=' and id_obs%'//trim(chMODULO)//trim(chCPU)

   NULNAM=4000
   ZEMFACT=0.01
   CALL fSQL_open( db, FILENAME ,stat)
   if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'fSQL_open: ', fSQL_errmsg(stat)
      call qqexit(2)
   endif
!=======================================================================
   query="select schema from rdb4_schema ;"
   RDB4_SCHEMA= SQL_QUERY_CH(db,trim(query))
   write(*,*) '  SUBROUTINE SQL2OBS_TOVS '
   write(*,*) '  RDB4 SCHEMA IS ---> ', trim(RDB4_SCHEMA)
!
   CHTIME=SQL_QUERY_CH(db,"PRAGMA journal_mode = OFF")
   CHTIME=SQL_QUERY_CH(db,"PRAGMA  synchronous = OFF")
   CHTIME=SQL_QUERY_CH(db,"PRAGMA PAGE_SIZE = 2048")
   CHTIME=SQL_QUERY_CH(db,"PRAGMA cache_size=300000")

!--- DEFAULT VALUES----------------------------------------------
   SQLEXTRA_HDR=""
   SQLEXTRA_DAT="" 
   SQLNULL     =""
   SQLLIMIT    =""
   SQLNULL     =" and  obsvalue is not null and vcoord is not null "
   SQLNULL     =" and vcoord is not null "
   LISTELEMENTS="12163"
   NBITSON     =0
   NBITSOFF    =0
!----------------------------------------------------------------

   COLUMNS_DAT="id_data,id_obs,vcoord,varno,obsvalue,flag "

!
!   READ FILTER PARAMETERS FROM NAMELIST FILE OF CORRESPONDING SCHEMA
!  -------------------------------------------------------------------


   NAMFILE=trim("flnml")
   open(NULNAM,FILE=NAMFILE)
   COLUMNS_HDR="id_obs,ID_SAT,ID_STN,lat,lon,codtyp,status,date,time"
   SELECT CASE(trim(RDB4_SCHEMA))
     !-------------------
     CASE( 'airs')
     !-------------------
       COLUMNS_HDR=trim(COLUMNS_HDR)//",LAND_SEA,INSTRUMENT,ZENITH,SOLAR_ZENITH,AZIMUTH,TERRAIN_TYPE,CLOUD_COVER"
       COLUMNS_DAT="id_data,id_obs,vcoord,varno,obsvalue,flag,SURF_EMISS "

       READ(NULNAM,NML=NAMSQLairs)
     !-------------------
     CASE( 'iasi')
     !-------------------
       COLUMNS_HDR=trim(COLUMNS_HDR)//",LAND_SEA,INSTRUMENT,ZENITH,SOLAR_ZENITH,AZIMUTH,TERRAIN_TYPE,CLOUD_COVER,SOLAR_AZIMUTH"
       COLUMNS_DAT="id_data,id_obs,vcoord,varno,obsvalue,flag,SURF_EMISS "
       READ(NULNAM,NML=NAMSQLiasi)
      !33060,33062,IGQISFLAGQUAL, SOLAR_AZIMUTH)
     !---------------------
    CASE( 'cris')
     !-------------------
       COLUMNS_HDR=trim(COLUMNS_HDR)//",LAND_SEA,INSTRUMENT,ZENITH,SOLAR_ZENITH,AZIMUTH,TERRAIN_TYPE,CLOUD_COVER,SOLAR_AZIMUTH"
       COLUMNS_DAT="id_data,id_obs,vcoord,varno,obsvalue,flag,SURF_EMISS "
       READ(NULNAM,NML=NAMSQLcris)

     !--------------------- 
    CASE( 'amsua')
     !---------------------
       COLUMNS_HDR=trim(COLUMNS_HDR)//",LAND_SEA,INSTRUMENT,ZENITH,SOLAR_ZENITH,AZIMUTH,TERRAIN_TYPE,SENSOR"
       READ(NULNAM,NML=NAMSQLamsua)
     !---------------------
     CASE( 'amsub')
     !---------------------
       COLUMNS_HDR=trim(COLUMNS_HDR)//",LAND_SEA,INSTRUMENT,ZENITH,SOLAR_ZENITH,AZIMUTH,TERRAIN_TYPE,SENSOR"
       READ(NULNAM,NML=NAMSQLamsub)
     CASE( 'atms')
     !---------------------
       COLUMNS_HDR=trim(COLUMNS_HDR)//",LAND_SEA,INSTRUMENT,ZENITH,SOLAR_ZENITH,AZIMUTH,TERRAIN_TYPE,SENSOR"
       READ(NULNAM,NML=NAMSQLatms)

     !-------------
     CASE( 'ssmi')
     !-------------
       COLUMNS_HDR=trim(COLUMNS_HDR)//",LAND_SEA,INSTRUMENT,ZENITH,SOLAR_ZENITH,AZIMUTH"
       READ(NULNAM,NML=NAMSQLssmi)

     !-----------
     CASE( 'go')
     !-----------
       COLUMNS_HDR=trim(COLUMNS_HDR)//",INSTRUMENT,ZENITH,SOLAR_ZENITH,LAND_SEA"
       COLUMNS_DAT="id_obs,id_data,vcoord,varno,obsvalue,flag,SURF_EMISS "
       READ(NULNAM,NML=NAMSQLgo)

     !-----------
     CASE( 'csr')
     !-----------
       COLUMNS_HDR=trim(COLUMNS_HDR)//",LAND_SEA,INSTRUMENT,ZENITH,SOLAR_ZENITH"
       READ(NULNAM,NML=NAMSQLcsr)

     !-----------
     CASE DEFAULT
     !-----------
       write(*,*) ' UNSUPPORTED SCHEMA ---> ',trim(RDB4_SCHEMA), ' exit SQL2OBS_TOVS '
       RETURN
   END SELECT
   close(NULNAM)

   WRITE(*,*)' SQLEXTRA_DAT =', SQLEXTRA_DAT
   BITSflagoff=0
   DO I = 1, nbitsoff
      BITSflagoff = IBSET ( BITSflagoff, 13-BITOFF(I) )
   END DO

   BITSflagon=0
   DO I = 1, nbitson
      BITSflagon = IBSET ( BITSflagon,  13-BITON(I)  )
   END DO

   !
   !-----------------------------------------------------------------------------------------
   !  COMPOSE SQL  QUERIES 
   !-----------------------------------------------------------------------------------------
   !
   WRITE(CFG, '(i6)' ) BITSflagoff
   CSQLcrit=trim("   flag & "//CFG)//" =0 "
   WRITE(CFG, '(i6)' ) BITSflagon
   if  (nbitson .gt. 0) then
     WRITE(CFG, '(i6)' ) BITSflagon
     CSQLcrit=trim(CSQLcrit)//trim("  and flag & "//CFG)
     CSQLcrit=trim(CSQLcrit)//" = "//CFG
   endif

   CSQLcrit=trim(CSQLcrit)//" and varno in ( "//trim(LISTELEMENTS)//" )"
!  CSQLcrit=trim(CSQLcrit)//trim("  and flag & "//CFG)
!  CSQLcrit=trim(CSQLcrit)//" = "//CFG

!  CSQLcrit=trim(CSQLcrit)//" and varno in ( "//trim(LISTELEMENTS)//" )"
   CSQLcrit=trim(CSQLcrit)//trim(SQLNull)
   CSQLcrit=trim(CSQLcrit)//trim(SQLEXTRA_DAT)

   !-----------------------------------------------------------------------------------------
   query_dat=                    "SELECT "//COLUMNS_DAT
   query_dat=trim(query_dat)//trim("   FROM DATA where  ")
   query_dat=       trim(query_dat)//"   "
   query_dat=trim(query_dat)//trim(CSQLcrit)
!  query_dat=trim(query_dat)//" order by id_obs "
   query_dat=trim(query_dat)//trim(SQLLIMIT)
   !-----------------------------------------------------------------------------------------
   query_hdr="SELECT "//COLUMNS_HDR
   query_hdr=trim(query_hdr)//" from HEADER "
   query_hdr=trim(query_hdr)//trim(SQLEXTRA_HDR)
   query_hdr=trim(query_hdr)//" order by id_obs "
   !-----------------------------------------------------------------------------------------
   WRITE(*, *)' TOVS QUERY_DAT --> ',trim(query_dat)
   WRITE(*, *)' TOVS QUERY_HDR --> ',trim(query_hdr)
   WRITE(*, *)' ========================================== '
   WRITE(*, *)' ========================================== '

!=========================
   VCORDTYP=3
!=========================
!pikpik   CHTIME=SQL_QUERY_CH(db,"select time('now')")
!pikpik   write(*,*) ' SQL2OBS_TOVS START OF QUERY TIME IS = ', CHTIME,RDB4_SCHEMA

! EXTRACT THE DATA FOM THE DATABASE 
!------------------------------------------------------
   CALL fSQL_prepare( db, trim(query_hdr) , stmt, stat)
   CALL fSQL_prepare( db, trim(query_dat), stmt2, stat2)
   if ( fSQL_error(stat)  /= FSQL_OK ) CALL handle_error(stat ,'fSQL_prepare hdr: ')
   if ( fSQL_error(stat2) /= FSQL_OK ) CALL handle_error(stat2,'fSQL_prepare dat: ')


   NOBS=0  
   COUNT=0
   NLV=0
   NOBS=obs_numHeader(obsdat)
   COUNT=obs_numBody(obsdat)
   NOBS_START=NOBS
 
   call tmg_start(55,'SQLMATRIXFILL')

   CALL fSQL_get_many (  stmt2, NROWS = NROWS , NCOLS = NCOLUMNS ,MODE = FSQL_REAL )
   write(*,*) '  NROWS NCOLUMNS =', NROWS,NCOLUMNS,RDB4_SCHEMA
   WRITE(*, *)' ========================================== '
   WRITE(*, *)' ========================================== '
   ALLOCATE( matdata (NROWS,NCOLUMNS) )
   CALL fSQL_fill_matrix ( stmt2, matdata)
   CALL fSQL_free_mem     (stmt2)
   CALL fSQL_finalize     (stmt2)
   call tmg_stop(55,'SQLMATRIXFILL')

!=============
! HEADER LOOP
!=============
   last_id=1
   HEADER:   DO
!
!    ======================================
!    exit 3DVAR WHEN LABORTFULL is TRUE
!    ======================================
      call obs_status(obsdat, obsdat_full, numstns_out, numobs_out, Kulout)
      IF ( obsdat_full) then

         IF (LABORTFULL) THEN
             WRITE(*,*)' SQL2OBS_TOVS: CMA FILE FULL  -HEADERS- CALL QQEXIT '
             CALL QQEXIT(2)
         ELSE
             WRITE(*,*)' SQL2OBS_TOVS: CMA FILE FULL   -HEADERS- EXIT ROUTINE '
             RETURN
         ENDIF
      ENDIF
!
! Fetch the next row
!======================================
     CALL fSQL_get_row( stmt, finished )
!
!

     IF (finished) then
     ! exit LOOP WHEN LAST ROW HAS BEEN FETCHED
     !========================================
       IF ( NOBS > 1 .and.  nlv > 0   ) THEN
         call obs_setFamily(obsdat,trim(familytype),NOBS)
         call obs_headSet_i(obsdat,OBS_IDO,nobs,ID_OBS)
         call obs_headSet_i(obsdat,OBS_ITY,nobs,CODTYP)
         call obs_headSet_i(obsdat,OBS_ONM,nobs,NOBS)
         CALL GEN_INFO(obsdat,ID_SAT,INSTRUMENT, INT(IZENITH),int(CLOUD_COVER),int(AZIMUTH),int(SOLAR_ZENITH), &
       &  INT(SOLAR_AZIMUTH),LAND_SEA,NOBS)
         XLAT=LAT*MPC_RADIANS_PER_DEGREE_R8
         XLON=LON*MPC_RADIANS_PER_DEGREE_R8
	 RELEV=ELEV
         CALL GEN_HEADER(obsdat,STNID,XLAT,XLON,RELEV,STATUS,DATE,TIME/100,NOBS)
       ENDIF
      EXIT                                           ! exit LOOP
     !========================================
     ENDIF
!

! The query result is inserted into variables
     !=========================================================================================
     TERRAIN_TYPE=MPC_missingValue_INT
     LAND_SEA=MPC_missingValue_INT
     SENSOR=MPC_missingValue_INT
     CLOUD_COVER=MPC_missingValue_R4
     ELEV=000.
     SOLAR_AZIMUTH=MPC_missingValue_R8
     SOLAR_ZENITH=MPC_missingValue_R8
     ZENITH=MPC_missingValue_R8
     CALL fSQL_get_column( stmt, COL_INDEX = 1, INT_VAR   = ID_OBS                     )
     CALL fSQL_get_column( stmt, COL_INDEX = 2, INT_VAR   = ID_SAT                     )
     CALL fSQL_get_column( stmt, COL_INDEX = 3, CHAR_VAR  = STNID                      )
     CALL fSQL_get_column( stmt, COL_INDEX = 4, REAL_VAR  = LAT                        )
     CALL fSQL_get_column( stmt, COL_INDEX = 5, REAL_VAR  = LON                        )
     CALL fSQL_get_column( stmt, COL_INDEX = 6, INT_VAR   = CODTYP                     )
     CALL fSQL_get_column( stmt, COL_INDEX = 7, INT_VAR   = STATUS                     )
     CALL fSQL_get_column( stmt, COL_INDEX = 8, INT_VAR   = DATE                       )
     CALL fSQL_get_column( stmt, COL_INDEX = 9, INT_VAR   = TIME                       )

     CALL fSQL_get_column( stmt, COL_INDEX = 10,INT_VAR   = LAND_SEA , INT_MISSING=MPC_missingValue_INT   )
     CALL fSQL_get_column( stmt, COL_INDEX = 11, INT_VAR  = INSTRUMENT )
     CALL fSQL_get_column( stmt, COL_INDEX = 12, REAL_VAR = ZENITH     )
     CALL fSQL_get_column( stmt, COL_INDEX = 13, REAL_VAR = SOLAR_ZENITH)
     CALL fSQL_get_column( stmt, COL_INDEX = 14, REAL_VAR = AZIMUTH,REAL_MISSING = MPC_missingValue_R4)
     CALL fSQL_get_column( stmt, COL_INDEX = 15, INT_VAR  = TERRAIN_TYPE,INT_MISSING=MPC_missingValue_INT )

     SELECT CASE(trim(RDB4_SCHEMA))
        CASE('airs','iasi')
            CALL fSQL_get_column( stmt, COL_INDEX = 16, REAL_VAR = CLOUD_COVER, REAL_MISSING = MPC_missingValue_R4  )
            CALL fSQL_get_column( stmt, COL_INDEX = 17, REAL_VAR = SOLAR_AZIMUTH, REAL_MISSING = MPC_missingValue_R4  )
        CASE('amsua','amsub','atms')
            CALL fSQL_get_column( stmt, COL_INDEX = 16, INT_VAR  = SENSOR,INT_MISSING=MPC_missingValue_INT         )
!       CASE('ssmi')
!        CASE('csr','go')

     END SELECT

!==========================================
!==========================================
   IF ( INSTRUMENT == 420 ) ID_SAT = 784
!==========================================


    IZENITH     =NINT ( (90.0 + ZENITH)*100       )
    SOLAR_ZENITH=NINT ( (90.0 + SOLAR_ZENITH)*100 )
    AZIMUTH     =NINT ( (AZIMUTH)*100             )
    CLOUD_COVER= NINT ( (CLOUD_COVER)*1           )
    SOLAR_AZIMUTH=NINT( (SOLAR_AZIMUTH)*100       )

!---Is terrain type sea ice (iterrain=0)?, If so, set imask=2.----
           IF ( TERRAIN_TYPE .EQ.  0       ) THEN
             LAND_SEA = 2
            else
           ENDIF

!==================================================================================
!==================================================================================
     IF ( trim(RDB4_SCHEMA)=='amsua' .or. trim(RDB4_SCHEMA)=='amsub' .or. trim(RDB4_SCHEMA)=='atms') THEN
!     CONVERT SENSOR TO INSTRUMENT
!
       if (SENSOR == MPC_missingValue_INT ) then
          SENSOR=0
          if (INSTRUMENT == -99 ) then
             INSTRUMENT=0
          endif
       else
          INSTRUMENT = cvt_burp_instrum2(sensor)
       endif

     ENDIF
!=================================================================================
!=================================================================================

     IF ( ID_SAT < 206 ) THEN
      ZENITH=9000
     ENDIF
     IF ( LON < 0.) LON= LON + 360.
     XLAT=LAT*MPC_RADIANS_PER_DEGREE_R8
     XLON=LON*MPC_RADIANS_PER_DEGREE_R8
     RELEV=ELEV

!==========================================================
!                ------
!    INSERT INTO header   OF OBS_SPACE_DATA STRUCTURE
!                ------
!==========================================================
     NOBS =NOBS  + 1 
     NLV=0
       

!===============================
    DATA: DO J =  last_id, NROWS
!=================================
!-------------------------------------------------------------
    if  ( int(matdata(j,2)) .gt. ID_OBS ) then
!-------------------------------------------------------------
     call obs_setFamily(obsdat,trim(familytype),NOBS)
     call obs_headSet_i(obsdat,OBS_IDO,nobs,ID_OBS)
     call obs_headSet_i(obsdat,OBS_ITY,nobs,CODTYP)
     call obs_headSet_i(obsdat,OBS_ONM,nobs,NOBS)
     CALL GEN_INFO(obsdat,ID_SAT,INSTRUMENT, INT(IZENITH),int(CLOUD_COVER),int(AZIMUTH),int(SOLAR_ZENITH), &
   &  INT(SOLAR_AZIMUTH),LAND_SEA,NOBS)
     CALL GEN_HEADER(obsdat,STNID,XLAT,XLON,RELEV,STATUS,DATE,TIME/100,NOBS)
!    ====
     EXIT
!    ====
!------------------------------------------------------------------
    else if  ( int(matdata(j,2)) == ID_OBS ) then
!------------------------------------------------------------------
   if ( NOBS == NOBS_START .and. NLV == 0  ) then
       call obs_setFamily(obsdat,trim(familytype),NOBS)
       call obs_headSet_i(obsdat,OBS_INS,nobs,INSTRUMENT )
       call obs_headSet_i(obsdat,OBS_IDO,nobs,ID_OBS)
       call obs_headSet_i(obsdat,OBS_ITY,nobs,CODTYP)
       call obs_headSet_i(obsdat,OBS_ONM,nobs,NOBS)
       CALL GEN_INFO(obsdat,ID_SAT,INSTRUMENT, INT(IZENITH),int(CLOUD_COVER),int(AZIMUTH),int(SOLAR_ZENITH), &
   &    INT(SOLAR_AZIMUTH),LAND_SEA,NOBS)
       CALL GEN_HEADER(obsdat,STNID,XLAT,XLON,RELEV,STATUS,DATE,TIME/100,NOBS)
!       print *, ' pikpik lasid=',ID_SAT,INSTRUMENT,INT(IZENITH),int(CLOUD_COVER),int(AZIMUTH),int(SOLAR_ZENITH),INT(SOLAR_AZIMUTH),LAND_SEA,NOBS
 endif

     last_id=j+1
     ID_DATA=int(matdata(j,1)); FLAG=int(matdata(j,6));  VARNO =int(matdata(j,4))
     VCOORDmat =matdata(j,3);   OBSVALUEmat=matdata(j,5)
     if (NCOLUMNS ==7) SURF_EMISS=matdata(j,7) 

!=========================================================================================
 
     COUNT = COUNT  + 1
     NLV= NLV +1

     call obs_status(obsdat, obsdat_full, numstns_out, numobs_out, Kulout)

! exit 3DVAR WHEN LABORTFULL is TRUE
!======================================
     IF ( obsdat_full) then
       IF (LABORTFULL) THEN
          WRITE(*,*)' SQL2OBS_TOVS: CMA FILE FULL - CALL QQEXIT '
          CALL QQEXIT(2)
       ELSE

!         -----------------------------------
          IF ( NLV .gt. 0) THEN
            call obs_headSet_i(obsdat,OBS_NLV,nobs,NLV)

            IF ( NOBS .gt. 1) THEN
              LN= obs_headElem_i(obsdat,OBS_RLN,NOBS-1) +  obs_headElem_i(obsdat,OBS_NLV,NOBS-1)
              call obs_headSet_i(obsdat,OBS_RLN,nobs,LN)
            ELSE
              call obs_headSet_i(obsdat,OBS_RLN,nobs,1)
            ENDIF

            write(*,*) ' CMA FILE FULL: EXIT  SQL2OBS_TOVS:  '
            RETURN

          ELSE
            write(*,*) ' CMA FILE FULL: EXIT  SQL2OBS_TOVS:  '
            RETURN
          ENDIF
!         -----------------------------------

       ENDIF
     ENDIF
!======================================

!pikpikenkf     FLAG = IBCLR(FLAG,12)
!==========================================================
!                -----
!    INSERT INTO data    OF OBS_SPACE_DATA STRUCTURE
!                -----
!==========================================================
     call obs_bodySet_r(obsdat,OBS_SEM,count,SURF_EMISS*ZEMFACT)
     call obs_bodySet_i(obsdat,OBS_IDD,count,ID_DATA)
     VCOORD=VCOORDmat
     OBSVALUE=OBSVALUEmat
     CALL GEN_DATA(obsdat,VCOORD,OBSVALUE,VARNO,FLAG,VCORDTYP,NOBS,count)
!
!----------------------------------------------
   endif !  id_obs   
!----------------------------------------------
!
!===============================
   END DO DATA  !END OF DATA LOOP
!===============================

   IF ( NOBS .eq. 1  ) THEN
      call obs_headSet_i(obsdat,OBS_RLN,nobs,1)
   ENDIF

   IF ( NLV .GT. 0 )THEN
     call obs_headSet_i(obsdat,OBS_NLV,nobs,NLV)
     IF ( NOBS .gt. 1) THEN
       LN= obs_headElem_i(obsdat,OBS_RLN,NOBS-1) +  obs_headElem_i(obsdat,OBS_NLV,NOBS-1)
       call obs_headSet_i(obsdat,OBS_RLN,nobs,LN)
     ENDIF
   ELSE
     NOBS=NOBS-1
   ENDIF

!=====================================
   END DO HEADER  !END OF HEADER LOOP
!=====================================
   deallocate(matdata)
   CALL fSQL_finalize(stmt)

   write(*,*)    ' NUMSTNS_numheader=', obs_numHeader(obsdat)
   write(*,*)    ' NUMOBS _numbody  =', obs_numBody(obsdat)
   write(*,*)    ' ================== END SQL2OBS_TOVS ================= '

! Close The SQLITE FILE
!==========================================
   CALL fSQL_close( db , stat)

!=============================
   END subroutine SQL2OBS_TOVS
!=============================


   subroutine SQL2OBS_CONV(obsdat,FAMILYTYPE,FILENAME,KULOUT)
      !
      !----------------------------------------------------------------
      !
      !   Purpose : QUERY DATA FROM SQLITE FILES
      !             WHICH CONTAIN DATA ON HEIGHT OR PRESSURE LEVELS
      !             OR SURFACE DATA.
      !               
      !             THEN INSERT THE DATA INTO   ObsSpaceData structures
      !
      ! Author  : P. Koclas, CMC/CMDA September  2012
      ! Adaptation : S. Skachko, ARMA March 2018
      !    ARGUMENTS:
      !                 INPUT:
      !
      !                       -OBSDAT  : ObsSpaceData Structure 
      !                       -FILENAME: sqlite file name  
      !                       -KULOUT  : STDOUT   FILE UNIT NUMBER 
      !----------------------------------------------------------------
      !
   implicit none
!----------------------------------------------------------------
!
   type (struct_obs), intent(inout)  :: obsdat
   CHARACTER *2,INTENT(in)           :: FAMILYTYPE
   CHARACTER *128,INTENT(in)         :: FILENAME
   INTEGER,INTENT(in)                :: KULOUT

   INTEGER            :: NULNAM
   CHARACTER *9       :: RDB4_SCHEMA

   INTEGER            :: ID_OBS,ID_DATA,CODTYP,DATE,TIME,STATUS,FLAG,VARNO
   CHARACTER *9       :: STNID

   REAL               :: LAT,LON
   REAL(OBS_REAL)     :: RELEV,XLAT,XLON
   REAL*8             :: OBSVALUE
   INTEGER            :: VCOORD,VCOORDFACT
   INTEGER            :: VCORDTYP

   INTEGER            :: I,J,COUNT,NLV,NOBS,NOBS_START
   LOGICAL            :: FINISHED
   REAL               :: ELEV,ELEVFACT

   CHARACTER *10      :: CHTIME
   CHARACTER *512     :: QUERY,QUERY_DAT,QUERY_HDR
   CHARACTER *256     :: CFG,CSQLCRIT,CSQLNULL,CSQLLIMIT
   CHARACTER *256     :: COLUMNS_HDR,COLUMNS_DAT
   INTEGER            :: BITSFLAGON,BITSFLAGOFF
   INTEGER            :: LN
   INTEGER            :: numstns_out, numobs_out
   LOGICAL            :: obsdat_full
   REAL , ALLOCATABLE :: matdata  (:,:)
   INTEGER            :: NROWS,NCOLUMNS,last_id
   
   REAL*8 GEOID_UNDULATION,EARTH_LOCAL_RAD_CURV
   INTEGER RO_QC_FLAG
!
!===========================================================
!
! type handle for  SQLIte file
   type(fSQL_DATABASE)                        :: db

!  prepared statement for  SQLite
   type(fSQL_STATEMENT)                       :: stmt,stmt2

!   error status 
   type(fSQL_STATUS)                          :: stat
!===========================================================
!
   CHARACTER *256     :: LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT
   INTEGER*4          :: NBITSOFF,NBITSON,BITOFF(15),BITON(15)
   CHARACTER *256     :: NAMFILE
!===========================================================================================================
   NAMELIST /NAMSQLua/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLai/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLsw/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLro/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLsfc/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLsc/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLpr/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
!===========================================================================================================

   write(KULOUT,*) '  SUBROUTINE SQL2OBS_CONV '
   NULNAM=4000
   CALL fSQL_open( db, FILENAME ,stat)
   if ( fSQL_error(stat) /= FSQL_OK ) then
     write(*,*) 'fSQL_open: ', fSQL_errmsg(stat)
     call qqexit(2)
   endif

   CHTIME=SQL_QUERY_CH(db,"select time('now')")

!
!   GET SCHEMA OF DATABASE FILE FROM resume table 
!  --------------------------------
   query="select schema from rdb4_schema ;"
   RDB4_SCHEMA= SQL_QUERY_CH(db,query)
   write(KULOUT,*) ' SQL2OBS_CONV START OF QUERY TIME IS = ', CHTIME,trim(RDB4_SCHEMA)

!--- DEFAULT VALUES----------------------------------------------
   NBITSON=0
   NBITSOFF=0
   LISTELEMENTS="12001,11001,11002,12192,10194"
   SQLEXTRA_HDR="  "
   SQLEXTRA_DAT=""
   SQLLIMIT    =""
   SQLNULL     =""
   VCOORDFACT=1
!----------------------------------------------------------------

   NAMFILE=trim("flnml")
   open(NULNAM,FILE=NAMFILE)
   SELECT CASE(trim(RDB4_SCHEMA))
       CASE ('ua')
          SQLNULL=" and obsvalue is not null  "
          VCOORDFACT=1
          VCORDTYP=2
	  LISTELEMENTS="12004,12203,10051,10004,11011,11012,11215,11216,12001,11001,11002,11003,11004,12192,10194,13210"
          READ(NULNAM,NML=NAMSQLua)
       CASE ('ai')
          SQLNULL=" and obsvalue is not null and vcoord is not null "
          VCOORDFACT=1
          VCORDTYP=2
          LISTELEMENTS="12001,11001,11002,11003,11004,12192"
          READ(NULNAM,NML=NAMSQLai)
       CASE ('sw')
          SQLNULL=" and obsvalue is not null and vcoord is not null "
          VCOORDFACT=1
          VCORDTYP=2
          LISTELEMENTS="11001,11002,11003,11004"
          READ(NULNAM,NML=NAMSQLsw)
       CASE ('pr')
          SQLNULL=" and obsvalue is not null and vcoord is not null "
          VCOORDFACT=1
          VCORDTYP=1
          LISTELEMENTS="11001,11002,11003,11004"
          READ(NULNAM,NML=NAMSQLpr)
       CASE ('ro')
          SQLNULL=" and obsvalue is not null and vcoord is not null "
          VCOORDFACT=1
          VCORDTYP=1
          LISTELEMENTS="15036"
          READ(NULNAM,NML=NAMSQLro)
       CASE ('sfc')
          SQLNULL=" and obsvalue is not null"
          VCOORDFACT=0
          VCORDTYP=1
          LISTELEMENTS="11011,11012,12004,12203,10051,10004"
          READ(NULNAM,NML=NAMSQLsfc)
       CASE ('scat')
          SQLNULL=" and obsvalue is not null"
          VCOORDFACT=0
          VCORDTYP=1
          LISTELEMENTS="11011,11012,11215,11216"
          READ(NULNAM,NML=NAMSQLsc)
       CASE DEFAULT
!===========================================
         write(KULOUT,*) ' Unsupported  SCHEMA ---> ',trim(RDB4_SCHEMA), ' exit SQL2OBS_CONV '
         return
   END SELECT
   WRITE(*,*) ' NULNAM NAMFILE=',NULNAM,NAMFILE
   close(NULNAM)
   WRITE(*,*) ' SQLEXTRA_HDR =',trim(SQLEXTRA_HDR)
   WRITE(*,*) ' SQLNULL      =',trim(SQLNULL)
   WRITE(*,*) ' SQLEXTRA_DAT =', SQLEXTRA_DAT
   WRITE(*,*) ' SQLLIMIT     =', TRIM(SQLLIMIT)

   BITSflagoff=0
   DO I = 1, nbitsoff
      BITSflagoff = IBSET ( BITSflagoff, 13-BITOFF(I) )
   END DO

   BITSflagon=0
   DO I = 1, nbitson
      BITSflagon = IBSET ( BITSflagon,  13-BITON(I)  )
   END DO

   !
   !-----------------------------------------------------------------------------------------
   !  COMPOSE SQL  QUERIES 
   !-----------------------------------------------------------------------------------------
   !
   WRITE(CFG, '(i6)' ) BITSflagoff
   CSQLcrit=trim("  flag & "//CFG)//" =0 "
   if  (nbitson .gt. 0) then
     WRITE(CFG, '(i6)' ) BITSflagon
     CSQLcrit=trim(CSQLcrit)//trim("  and flag & "//CFG)
     CSQLcrit=trim(CSQLcrit)//" = "//CFG
   endif
   CSQLcrit=trim(CSQLcrit)//" and varno in ( "//trim(LISTELEMENTS)//" )"
   CSQLcrit=trim(CSQLcrit)//trim(SQLNull)
   CSQLcrit=trim(CSQLcrit)//trim(SQLEXTRA_DAT)

   !----------------------------------------------------------------------------------------------------
   COLUMNS_HDR="id_obs,lat,lon,codtyp,elev,date,time,status,id_stn"
   COLUMNS_DAT="id_data,id_obs,vcoord,varno,obsvalue,flag "
   !----------------------------------------------------------------------------------------------------

   !-----------------------------------------------------------------------------------------
   query_dat=                    "SELECT "//COLUMNS_DAT
   query_dat=trim(query_dat)//trim("   FROM DATA where  ")
   query_dat=trim(query_dat)//trim(CSQLcrit)
   query_dat=trim(query_dat)//" order by id_obs "
   query_dat=trim(query_dat)//trim(SQLLIMIT)
   !-----------------------------------------------------------------------------------------
   query_hdr="SELECT "//COLUMNS_HDR
   query_hdr=trim(query_hdr)//" from HEADER "
   query_hdr=trim(query_hdr)//trim(SQLEXTRA_HDR)
   query_hdr=trim(query_hdr)//" order by id_obs "
   !-----------------------------------------------------------------------------------------
   WRITE(*,*)' CONV QUERY_DAT --> ',trim(query_dat)
   WRITE(*,*)' CONV QUERY_HDR --> ',trim(query_hdr)
   !-----------------------------------------------------------------------------------------

   if (  trim(RDB4_SCHEMA)=='pr' .or. trim(RDB4_SCHEMA)=='sfc' ) then
     ELEVFACT=1.
   else
     ELEVFACT=0.
   endif

   CALL fSQL_prepare( db, trim(query_hdr) , stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare: ')

   NOBS=0
   COUNT=0
   NOBS=obs_numHeader(obsdat)
   COUNT=obs_numBody(obsdat)
   NOBS_START=NOBS

   NLV=0
   CALL fSQL_prepare  ( db, trim(query_dat), stmt2, stat)
   CALL fSQL_get_many (  stmt2, NROWS = NROWS , NCOLS = NCOLUMNS ,MODE = FSQL_REAL )
   write(KULOUT,*) '  NROWS NCOLUMNS =', NROWS,NCOLUMNS,RDB4_SCHEMA,trim(query_dat)
   ALLOCATE( matdata (NROWS,NCOLUMNS) )
   CALL fSQL_fill_matrix (stmt2, matdata )
   CALL fSQL_free_mem    (stmt2)
   CALL fSQL_finalize    (stmt2)
!=============
! HEADER LOOP
!=============
   last_id=1
   HEADER:   DO
!
!
! Fetch the next row
!======================================
   CALL fSQL_get_row( stmt, finished )
!
!
!========================================
   IF (finished) then
!========================================
     if ( NOBS > 1 .and.  nlv > 0 ) then
       call obs_setFamily(obsdat,trim(FAMILYTYPE),NOBS)
       call obs_headSet_i(obsdat,OBS_IDO,nobs,ID_OBS)
       call obs_headSet_i(obsdat,OBS_ITY,nobs,CODTYP)
       call obs_headSet_i(obsdat,OBS_ONM,nobs,NOBS)
       write(*,*) ' genheader nobsgt1 end',finished
       XLAT=LAT*MPC_RADIANS_PER_DEGREE_R8
       XLON=LON*MPC_RADIANS_PER_DEGREE_R8
       RELEV=ELEV
       write(*,*) 'finished  STNID,XLAT,XLON,RELEV,STATUS,DATE,TIME/100,NOBS=',STNID,XLAT,XLON,RELEV,STATUS,DATE,TIME/100,NOBS
       call GEN_HEADER(obsdat,STNID,XLAT,XLON,RELEV,STATUS,DATE/100,TIME,NOBS)
     endif
!      exit LOOP WHEN LAST ROW HAS BEEN FETCHED
!      ====
     EXIT         ! exit LOOP
!      ====
   ENDIF

   CALL fSQL_get_column( stmt, COL_INDEX = 1, INT_VAR   = ID_OBS                     )
   CALL fSQL_get_column( stmt, COL_INDEX = 2, REAL_VAR  = LAT                       )
   CALL fSQL_get_column( stmt, COL_INDEX = 3, REAL_VAR  = LON                        )
   CALL fSQL_get_column( stmt, COL_INDEX = 4, INT_VAR   = CODTYP                     )
   CALL fSQL_get_column( stmt, COL_INDEX = 5, REAL_VAR  = ELEV                       )
   RELEV=ELEV
   CALL fSQL_get_column( stmt, COL_INDEX = 6, INT_VAR   = DATE                       )
   CALL fSQL_get_column( stmt, COL_INDEX = 7, INT_VAR   = TIME                       )
   CALL fSQL_get_column( stmt, COL_INDEX = 8, INT_VAR   = STATUS                     )
   CALL fSQL_get_column( stmt, COL_INDEX = 9, CHAR_VAR  = STNID                      )

   if (trim(RDB4_SCHEMA)=='ro') then
     CALL fSQL_get_column( stmt, COL_INDEX = 10, INT_VAR  = RO_QC_FLAG   )
     CALL fSQL_get_column( stmt, COL_INDEX = 11, REAL8_VAR = GEOID_UNDULATION     )
     CALL fSQL_get_column( stmt, COL_INDEX = 12, REAL8_VAR = EARTH_LOCAL_RAD_CURV     )
   endif

   call obs_status(obsdat, obsdat_full, numstns_out, numobs_out, Kulout)
   IF (  obsdat_full ) then
! exit 3DVAR WHEN LABORTFULL is TRUE
!======================================
     IF (LABORTFULL) THEN
       WRITE(KULOUT,*)' HEADER FULL CONV',numstns_out,numobs_out
       WRITE(KULOUT,*)' SQL2OBS_CONV: CMA FILE FULL - CALL QQEXIT '
       CALL QQEXIT(2)
     ELSE
       WRITE(KULOUT,*)'HEADER FULL CONV',numstns_out,numobs_out
       WRITE(KULOUT,*)' SQL2OBS_CONV: CMA FILE FULL - CONTINUE '
       RETURN
     ENDIF
   ENDIF
!======================================

!=========================================================================================

   IF ( LON .LT. 0.) LON= LON + 360.
   XLAT=LAT*MPC_RADIANS_PER_DEGREE_R8
   XLON=LON*MPC_RADIANS_PER_DEGREE_R8
!==========================================================
!                ------
!    INSERT INTO header   OF OBS_SPACE_DATA STRUCTURE
!                ------
!==========================================================
   NLV=0
   write(*,*) ' last_id NOBS ID_OBS CODTYP =', last_id,NOBS,ID_OBS,CODTYP

   NOBS =NOBS  + 1 
!===============================
   DATA: DO J =  last_id, NROWS
!===============================
   if  ( int(matdata(j,2)) > ID_OBS ) then
     call obs_setFamily(obsdat,trim(FAMILYTYPE),NOBS)
     call obs_headSet_i(obsdat,OBS_IDO,nobs,ID_OBS)
     call obs_headSet_i(obsdat,OBS_ITY,nobs,CODTYP)
     call obs_headSet_i(obsdat,OBS_ONM,nobs,NOBS)
     RELEV=ELEV
     write(*,*) 'exit  STNID,XLAT,XLON,RELEV,STATUS,DATE,TIME/100,NOBS=',STNID,XLAT,XLON,RELEV,STATUS,DATE,TIME/100,NOBS,NLV
     call GEN_HEADER(obsdat,STNID,XLAT,XLON,RELEV,STATUS,DATE,TIME/100,NOBS)
     !====
     EXIT  DATA
     !====
   else if  ( int(matdata(j,2)) == ID_OBS ) then
     if ( NOBS==NOBS_START .and. NLV==0  ) then
       call obs_setFamily(obsdat,trim(FAMILYTYPE),NOBS)
       call obs_headSet_i(obsdat,OBS_IDO,nobs,ID_OBS)
       call obs_headSet_i(obsdat,OBS_ITY,nobs,CODTYP)
       call obs_headSet_i(obsdat,OBS_ONM,nobs,NOBS)
       RELEV=ELEV
       write(*,*) 'icitt  STNID,XLAT,XLON,RELEV,STATUS,DATE,TIME/100,NOBS=',STNID,XLAT,XLON,RELEV,STATUS,DATE,TIME/100,NOBS
       call GEN_HEADER(obsdat,STNID,XLAT,XLON,RELEV,STATUS,DATE,TIME/100,NOBS)
     endif
     last_id=j+1
     ID_DATA=int(matdata(j,1)); FLAG=int(matdata(j,6));  VARNO =int(matdata(j,4))
     VCOORD=matdata(j,3); OBSVALUE=matdata(j,5)

     call obs_status(obsdat, obsdat_full, numstns_out, numobs_out, Kulout)

! exit 3DVAR WHEN LABORTFULL is TRUE
!
!    ========================
     IF (  obsdat_full ) THEN
!    ========================
!      --------------------------------
       IF (LABORTFULL) THEN
!      --------------------------------
         WRITE(KULOUT,*)' DATA FULL CONV',count,numobs_out
         WRITE(KULOUT,*)' SQL2OBS_CONV: CMA FILE FULL - CALL QQEXIT '
         call qqexit(2)
!      --------------------------------
       ELSE
!      --------------------------------
         WRITE(KULOUT,*)' SQL2OBS_CONV: CMA FILE FULL - CONTINUE '
         WRITE(*,*)' DATA FULL CONV RETURN',numstns_out,numobs_out,nobs
         IF ( NLV .gt. 0) THEN
!        =============================================
           call obs_headSet_i(obsdat,OBS_NLV,nobs,NLV)
           IF ( NOBS .gt. 1) THEN
             LN= obs_headElem_i(obsdat,OBS_RLN,NOBS-1) +  obs_headElem_i(obsdat,OBS_NLV,NOBS-1)
             call obs_headSet_i(obsdat,OBS_RLN,nobs,LN)
           ELSE
             call obs_headSet_i(obsdat,OBS_RLN,nobs,1)
           ENDIF
           RETURN
!        =============================================
         ELSE
           RETURN
         ENDIF
!        =============================================
!      --------------------------------
       ENDIF
!      --------------------------------

!    =====================
     ENDIF
!    =====================

   COUNT = COUNT  + 1
   NLV= NLV +1
!==========================================================
!                -----
!    INSERT INTO data    OF OBS_SPACE_DATA STRUCTURE
!                -----
!==========================================================

   call obs_bodySet_i(obsdat,OBS_IDD,count,ID_DATA)
   CALL GEN_DATA(obsdat,VCOORD*VCOORDFACT + RELEV*ELEVFACT,OBSVALUE,VARNO,FLAG,VCORDTYP,NOBS,count)

!pikpikenkf   FLAG = IBCLR(FLAG,12)

!     ALLOW EXTRA SPACE FOR U V COMPONENTS
!--------------------------------------------------------------
   if ( varno .eq. 11001 .or. varno .eq. 11011) then
     call obs_status(obsdat, obsdat_full, numstns_out, numobs_out, Kulout)

! exit 3DVAR WHEN LABORTFULL is TRUE
!
!-----------------------
     IF (  obsdat_full ) then
!-----------------------
       IF (LABORTFULL) THEN
         WRITE(KULOUT,*)' SQL2OBS_CONV: WINDS CMA FILE FULL - CALL QQEXIT '
         call QQEXIT(2)
       ELSE
         WRITE(*,*)' WINDS CMAFULL',obs_numBody(obsdat)
!        =============================================
         IF ( NLV .gt. 0) THEN
!        =============================================
           call obs_headSet_i(obsdat,OBS_NLV,nobs,NLV)
           IF ( NOBS .gt. 1) THEN
             LN= obs_headElem_i(obsdat,OBS_RLN,NOBS-1) +  obs_headElem_i(obsdat,OBS_NLV,NOBS-1)
             call obs_headSet_i(obsdat,OBS_RLN,nobs,LN)
           ELSE
             call obs_headSet_i(obsdat,OBS_RLN,nobs,1)
           ENDIF
           WRITE(KULOUT,*)' WINDS SQL2OBS_CONV: CMA FILE FULL - RETURN '
           RETURN
!        =============================================
         ENDIF
!        =============================================
         !=====
         RETURN
         !=====
       ENDIF
!-----------------------
     ENDIF
!-----------------------
            
     if ( varno .eq. 11001 ) then
!      U COMPONENT

       call obs_bodySet_i(obsdat,OBS_IDD,count+1,-1)
       CALL GEN_DATA(obsdat,VCOORD*VCOORDFACT + RELEV*ELEVFACT,MPC_missingValue_R8,11003,0,VCORDTYP,NOBS,count+1)

!      V COMPONENT

       call obs_bodySet_i(obsdat,OBS_IDD,count+2,-1)
       CALL GEN_DATA(obsdat,VCOORD*VCOORDFACT + RELEV*ELEVFACT,MPC_missingValue_R8,11004,0,VCORDTYP,NOBS,count+2)
     else if (  varno .eq. 11011) then

!       Us COMPONENT

       call obs_bodySet_i(obsdat,OBS_IDD,count+1,-1)
       CALL GEN_DATA(obsdat,VCOORD*VCOORDFACT + RELEV*ELEVFACT,MPC_missingValue_R8,11215,0,VCORDTYP,NOBS,count+1)

!      Vs COMPONENT

       call obs_bodySet_i(obsdat,OBS_IDD,count+2,-1)
       CALL GEN_DATA(obsdat,VCOORD*VCOORDFACT + RELEV*ELEVFACT,MPC_missingValue_R8,11216,0,VCORDTYP,NOBS,count+2)

     endif
     count = count + 2
     NLV = NLV + 2
   endif

 endif
!================================
 END DO DATA  !END OF DATA LOOP
!================================

!===================================================================

   if ( NOBS==1 ) then
     call obs_headSet_i(obsdat,OBS_RLN,nobs,1)
   endif

   if ( NLV>0 ) then
     call obs_headSet_i(obsdat,OBS_NLV,nobs,NLV)
     if ( NOBS>1) then
       LN= obs_headElem_i(obsdat,OBS_RLN,NOBS-1) +  obs_headElem_i(obsdat,OBS_NLV,NOBS-1)
       call obs_headSet_i(obsdat,OBS_RLN,nobs,LN)
     endif
   else
     NOBS=NOBS-1
!    write(KULOUT,*)  ' enleve obs NOBS',nobs,obs_numHeader(obsdat)
   endif
   write(*,*) ' finished,last_id nrows nobs =',finished,last_id,nrows,nobs
   if ( last_id > nrows) then
     write(*,*)  ' avant exit  numheader  =', obs_numHeader(obsdat)
     write(*,*)  ' avant exit  numbody    =', obs_numBody(obsdat)
     call obs_setFamily(obsdat,trim(FAMILYTYPE),NOBS)
     call obs_headSet_i(obsdat,OBS_IDO,nobs,ID_OBS)
     call obs_headSet_i(obsdat,OBS_ITY,nobs,CODTYP)
     call obs_headSet_i(obsdat,OBS_ONM,nobs,NOBS)

     if (trim(RDB4_SCHEMA)=='ro') then
       if ( obs_columnActive_IH(obsdat,OBS_ROQF) ) call obs_headSet_i(obsdat,OBS_ROQF,nobs,RO_QC_FLAG)
       if ( obs_columnActive_RH(obsdat,OBS_TRAD) ) call obs_headSet_r(obsdat,OBS_TRAD,nobs,EARTH_LOCAL_RAD_CURV)
       if ( obs_columnActive_RH(obsdat,OBS_GEOI) ) call obs_headSet_r(obsdat,OBS_GEOI,nobs,GEOID_UNDULATION)
     endif


     CALL GEN_HEADER(obsdat,STNID,XLAT,XLON,RELEV,STATUS,DATE,TIME/100,NOBS)
     write(*,*)  ' exit  numheader  =', obs_numHeader(obsdat)
     write(*,*)  ' exit  numbody    =', obs_numBody(obsdat)
     exit
   endif


!=========================
  END DO HEADER ! HEADER
!=========================
  deallocate(matdata)
!  CALL fSQL_free_mem(stmt)

  write(KULOUT,*)  ' FIN numheader  =', obs_numHeader(obsdat)
  write(KULOUT,*)  ' FIN numbody    =', obs_numBody(obsdat)
  write(KULOUT,*)  ' fin header '

  CHTIME=SQL_QUERY_CH(db,"select time('now')")
  write(KULOUT,*) ' SQL2OBS_CONV END OF QUERY TIME IS = ', CHTIME,RDB4_SCHEMA

  CALL fSQL_finalize( stmt )
! Close The SQLITE FILE
!==========================================
  CALL fSQL_close( db , stat)

!==============================
END subroutine SQL2OBS_CONV
!==============================

!  ==============================

!===================================================================

   function SQL_QUERY_CH(db,query)
      !
      !    ---------------------------------------------------
      !   Purpose : RETURN RESULT OF A QUERY  TO AND OPENED SQLITE FILE
      !             THE RESULT FROM THIS QURY MUS BE A CHARATER STRING
      !
      ! Author  : P. Koclas, CMC/CMDA September  2012
      !
      !    ARGUMENTS:
      !                 INPUT:
      !
      !                       -db      : SQLITE FILE HANDLE
      !                       -query   : A QUERY
      !    ---------------------------------------------------
      !

   implicit none
   CHARACTER(len = 256 )                          :: SQL_QUERY_CH
   CHARACTER(len = 256 )                          :: CH_RESULT
   
   logical finished
! type handle for  SQLIte file
!
   type(fSQL_DATABASE)                      :: db
!  prepared statement for  SQLite
!
   type(fSQL_STATEMENT)                     :: stmt
!type error status
!
   type(fSQL_STATUS)                        :: stat
   CHARACTER(len = *)                       :: query

   CH_RESULT=''

! prepare query for execution
! ---------------------------------------------------------------------------
   CALL fSQL_prepare( db, trim(query), stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare: ')

! Execute query
! ---------------------------------------------------------------------------
   finished=.FALSE.
   CALL fSQL_get_row( stmt, finished )

!  Put result of query into variable
! ---------------------------------------------------------------------------
   CALL fSQL_get_column( stmt, COL_INDEX = 1, CHAR_VAR   = CH_RESULT    )
   CALL fSQL_get_row( stmt, finished )
   if ( .NOT. finished) then
      write(*,*)' SQL_QUERY_CH ---> QUERY RETURNS MORE THAN ONE ROW...  '
    endif
   CALL fSQL_finalize( stmt )

!===================================
   SQL_QUERY_CH = trim(CH_RESULT)
!===================================

   RETURN

!  ==============================
   end function   SQL_QUERY_CH
!  ==============================
  integer function cvt_burp_instrum2(sensor)
    !
    ! func CVT_BURP_INSTRUM : Map burp satellite sensor indicator (element #2048) to
    !                         the corresponding burp satellite instrument (element
    !                         #2019).
    !
    ! Author  : J. Halle CMDA/SMC May 2002
    ! Revision:
    !           J.W. Blezius ARMA Feb 2011 - converted from subroutine to a function
    !
    !     Purpose:  Map burp satellite sensor indicator (element #2048) to
    !               burp satellite instrument (element #2019). This is a more
    !               complete common element, allowing for future expansion.
    !
    !               Table of BURP satellite sensor indicator element #002048
    !               --------------------------------------------------------
    !               Satellite sensor     BURP satellite sensor indicator
    !               ----------------     -------------------------------
    !               HIRS                 0
    !               MSU                  1
    !               SSU                  2
    !               AMSUA                3
    !               AMSUB                4
    !               AVHRR                5
    !               SSMI                 6
    !               NSCAT                7
    !               SEAWINDS             8
    !               Reserved             9-14
    !               Missing value        15
    !

    implicit none
    integer :: sensor      ! BURP satellite sensor indicator (element #2048)
    integer :: instrument  ! BURP satellite instrument       (element #2019)

    select case (sensor)
      case (000);   instrument=606  ! HIRS
      case (001);   instrument=623  ! MSU
      case (002);   instrument=627  ! SSU
      case (003);   instrument=570  ! AMSUA
      case (004);   instrument=574  ! AMSUB
      case (005);   instrument=591  ! AVHRR
      case (006);   instrument=905  ! SSMI
      case (007);   instrument=312  ! NSCAT
      case (008);   instrument=313  ! SEAWINDS
      case (015);   instrument=2047 ! Missing value
      case default; instrument=2047 ! Unrecognized value
    end select

    cvt_burp_instrum2 = instrument
    return

  end function cvt_burp_instrum2

   integer function cvt_burp_instrum(sensor)
      !
      ! func CVT_BURP_INSTRUM : Map burp satellite sensor indicator (element #2048) to
      !                         the corresponding burp satellite instrument (element
      !                         #2019).
      !
      ! Author  : J. Halle CMDA/SMC May 2002
      ! Revision:
      !           J.W. Blezius ARMA Feb 2011 - converted from subroutine to a function
      !
      !     Purpose:  Map burp satellite sensor indicator (element #2048) to
      !               burp satellite instrument (element #2019). This is a more
      !               complete common element, allowing for future expansion.
      !
      !               Table of BURP satellite sensor indicator element #002048
      !               --------------------------------------------------------
      !               Satellite sensor     BURP satellite sensor indicator
      !               ----------------     -------------------------------
      !               HIRS                 0
      !               MSU                  1
      !               SSU                  2
      !               AMSUA                3
      !               AMSUB                4
      !               AVHRR                5
      !               SSMI                 6
      !               NSCAT                7
      !               SEAWINDS             8
      !               Reserved             9-14
      !               Missing value        15
      !
   implicit none
   integer :: sensor      ! BURP satellite sensor indicator (element #2048)
   integer :: instrument  ! BURP satellite instrument       (element #2019)

   select case (sensor)
      case (000);   instrument=606  ! HIRS
      case (001);   instrument=623  ! MSU
      case (002);   instrument=627  ! SSU
      case (003);   instrument=570  ! AMSUA
      case (004);   instrument=574  ! AMSUB
      case (005);   instrument=591  ! AVHRR
      case (006);   instrument=905  ! SSMI
      case (007);   instrument=312  ! NSCAT
      case (008);   instrument=313  ! SEAWINDS
      case (015);   instrument=2047 ! Missing value
      case default; instrument=2047 ! Unrecognized value
   end select

!==================================
   cvt_burp_instrum = instrument
!==================================
   return
!  =============================
   end function cvt_burp_instrum
!  =============================

   subroutine handle_error(stat, message)
      !
      ! Exit code on f90sqlite 
      ! error message
      ! Author  : P. Koclas CMDA/CMC SEPTEMBER 2012
      !

   type(FSQL_STATUS)  :: stat
   character(len = *) :: message
   write(*,*) message, fSQL_errmsg(stat)
!===================
   call qqexit(2)
!===================

!  ===========================
   END SUBROUTINE handle_error
!  ===========================
   subroutine updsql(db,obsdat,familyType,fileName,file_numb)
      !
      !   Purpose : UPDATE HEADER AND DATA TABLES OF SQLITE FILES
      !
      ! Author  : P. Koclas, CMC/CMDA CMDA/CMC SEPTEMBER 2012
      !
      !    ARGUMENTS:
      !                 INPUT:
      !
      !                       -db        : SQLITE FILE HANDLE
      !                       -obsdat    : ObsSpaceData_mod OBJECT
      !                       -file_numb : FILE NUMBER ASSOCIATED WITH db
      !                       -ITEMLIST  : LIST OF ITEMS ('OMA',OMP','OER','FLG','FGE','VAR')
      !
      !    Rows of sqlite files are updated via a list of items specified in argument itemlist
      !

   implicit none
!
!
   ! arguments
   type(fSQL_DATABASE)              :: db
   type (struct_obs), intent(inout) :: obsdat
   character(len=*) :: fileName
   character(len=*) :: familyType
   integer :: file_numb

!  prepared statement for  SQLite
   type(fSQL_STATEMENT)                     :: stmt

!type error status
   type(fSQL_STATUS)                        :: stat
   
!  local variables 
!   integer*4                        :: file_numb ! 
   integer*4                        :: RLN,NLV,IDF,ID_DATA,ASS,FLAG ! 
   integer*4                        :: IOBS,ID_OBS,STATUS ! 
   integer*4                        :: J,JO !
   integer*4                        :: last_comma ! 
   integer*4                        :: ibegin,ilast,ibeginob,ilastob,NUPD ! 
   CHARACTER *10                    :: CHTIME
   character(len  =   3)            :: item
   integer*4                        :: upd_list(20)
   character(len  =   9)            :: item2
   character(len  = 128)            :: query
   character(len  = 356)            :: CHitem,CHitem2
   integer                          :: UPDLIST
   logical                          :: BACK
   REAL*4                           :: ROMP,OBSV

! Variables FOR LOOPS 
   integer*4                        :: I,K
!----------------------------------------------------------

   write(*,*) '  ======================================================================= '
   write(*,*) '  ==============BEGIN UPDATE============================================== '
   write(*,*) '  ======================================================================= '
   WRITE(*,*) '  ------- UPDSQL ------ ',itemlist(1)
   WRITE(*,*) '  FAMILYTYPE   =', FAMILYTYPE
   WRITE(*,*) '  FileName     =', FileName
   WRITE(*,*) '  MISSING VALUE  = ', MPC_missingValue_R8

   
   !    CREATE QUERY  
   !--------------------------------
   CHitem='  '
   UPDLIST=N_ITEMS
   if (UPDLIST == 0) then
      write(*,*) ' NO UPDATES TO DO : ---- RETURN'
      return
   endif
   do K=1,UPDLIST
      item=itemlist(k)
   
      SELECT CASE(item)
       CASE('OMA')
          upd_list(K) = OBS_OMA
          item2='oma'

       CASE('OMP')
          upd_list(K) = OBS_OMP
          item2='omp'

       CASE('VAR')
          upd_list(K) = OBS_VAR
          item2='obsvalue'

        CASE('OER')
          upd_list(K) = OBS_OER
          item2='obs_error'

        CASE('FGE')
          upd_list(K) = OBS_HPHT
          item2='fg_error'

        CASE('FLG')
          upd_list(K) = OBS_FLG
          item2='flag'

        CASE DEFAULT

!===============================================
         WRITE(*,*)'invalid item: ', item2
         WRITE(*,*)' EXIT ROUTINE: '
         RETURN
!===============================================
       END SELECT
       CHitem = trim(CHitem)//trim(item2)//trim(' = ? ,')
       WRITE(*, "(A9)") trim(item2)
   end do
   back=.true.
   last_comma=SCAN(CHitem, ',', back)
   CHitem2   =CHitem(1:last_comma-1) 
   CHitem    =CHitem2

   query = ' UPDATE data SET flag = ? , '//trim(CHitem)
   query=trim(query)
   query = trim(query)//' WHERE id_data = ?  ;'
   write(*,*) ' QUERY == --->  ',query


   CALL fSQL_prepare( db, query , stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare : ')
   CHTIME=SQL_QUERY_CH(db,"PRAGMA journal_mode = OFF")
   CHTIME=SQL_QUERY_CH(db,"PRAGMA  synchronous = OFF")
   CHTIME=SQL_QUERY_CH(db,"PRAGMA locking_mode =  exclusive ")
   CHTIME=SQL_QUERY_CH(db,"PRAGMA PAGE_SIZE = 4096")
   CHTIME=SQL_QUERY_CH(db,"PRAGMA cache_size=500000")

   write(*,*) ' NUMBER OF HEADERS in obsdat    =',obs_numHeader(obsdat)
!---------------------------------------------------------------------------
! Begin transaction 
!---------------------------------------------------------------------------
   CALL fSQL_begin(db)

!========================================================================================
!
!--------------------------------
!   HEADER LOOP
!--------------------------------
   NUPD=0
   DO JO=1,obs_numHeader(obsdat)

      IDF=obs_headElem_i(obsdat,OBS_IDF,JO)
      if ( IDF .ne. FILE_NUMB) cycle
      ID_OBS=obs_headElem_i(obsdat,OBS_IDO,JO)
      RLN   =obs_headElem_i(obsdat,OBS_RLN,JO)
      NLV   =obs_headElem_i(obsdat,OBS_NLV,JO)

!--------------------------------
! TOP DATA LOOP
!--------------------------------
      DO J = RLN, NLV + RLN -1

          FLAG   =obs_bodyElem_i(obsdat,OBS_FLG,j)
          ID_DATA=obs_bodyElem_i(obsdat,OBS_IDD,j)
          ASS    =obs_bodyElem_i(obsdat,OBS_ASS,j)

!------------------------------------------------------------------------------------
       if ( ASS == 1 ) then
         call fSQL_bind_param(stmt, PARAM_INDEX = 1,   INT_VAR  = flag  )
         DO K =1,UPDLIST  
           OBSV =obs_bodyElem_r(obsdat,OBS_VAR,j)
           if ( OBSV /= MPC_missingValue_R8 ) then  
             ROMP= obs_bodyElem_r(obsdat,UPD_LIST(K),j)
             CALL fSQL_bind_param(stmt, PARAM_INDEX = K+1, REAL_VAR = ROMP   )
             CALL fSQL_bind_param(stmt, PARAM_INDEX = K+2, INT_VAR  = ID_DATA  )
             CALL fSQL_exec_stmt (stmt)
             NUPD=NUPD + 1
           endif
         END DO
       endif
 
!------------------------------------------------------------------------------------
!
       END DO
    END DO
!========================================================================================
!

!---------------------------------------------------------------------------
   CALL fSQL_finalize( stmt )

   WRITE(*,*) ' NUMBER OF DATA ROWS UPDATES =',NUPD
!---------------------------------------------------------------------------
   

! UPDATES FOR THE STATUS FLAGS IN THE HEADER TABBLE
!
!==============================================================================
!
!
   NUPD=0
   query = ' UPDATE header SET status  = ?   WHERE ID_OBS = ? '
   CALL fSQL_prepare( db, query , stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare : ')

   DO JO = 1,obs_numHeader(obsdat)
      IDF=obs_headElem_i(obsdat,OBS_IDF,JO)
      if ( IDF .ne. FILE_NUMB) cycle
      ID_OBS=obs_headElem_i(obsdat,OBS_IDO,JO)
      STATUS=obs_headElem_i(obsdat,OBS_ST1,JO)
      CALL fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = STATUS )
      CALL fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = ID_OBS )
      CALL fSQL_exec_stmt ( stmt)
      NUPD=NUPD + 1
   END DO
   WRITE(*,*) ' NUMBER OF HEADER ROW UPDATES =',NUPD
!
!
!==============================================================================
!
   CALL fSQL_finalize( stmt )
   CALL fSQL_commit(db)

   write(*,*) '  ===============END UPDATE============================================== '
   write(*,*) '  ======================================================================= '

!  =====================
   END SUBROUTINE updsql
!  =====================

   subroutine insertsql(db,obsdat,familyType,fileName,file_numb)
      !
      !   Purpose : INSERT NEW DATA ROWS INTO SQLITE FILES
      !
      ! Author  : P. Koclas, CMC/CMDA September 2012
      !
      !    ARGUMENTS:
      !                 INPUT:
      !
      !                       -db         : SQLITE FILE HANDLE
      !                       -obsdat     : ObsSpaceData OBJECT instance
      !                       -FILE_NUMB  : FILE NUMBER  associated to physical sqlite file
      !
      !    New data rows are inserted into sqlite files via a selected list of BUFR elements
      !     specified in a namlist
      !    Only the data with the id_DATA= MISSING_VALUE(PPMIS) is inserted:
      !                       (elements calculated by the 3dvar)
      !

   implicit none
   ! arguments
   ! type for SQLIte  file handle
   type(fSQL_DATABASE)    :: db
   type(struct_obs)       :: obsdat
   character(len=*)       :: familyType
   character(len=*)       :: fileName
   integer                :: file_numb

! type for precompiled SQLite statements
   type(fSQL_STATEMENT)                     :: stmt

!type for error status
   type(fSQL_STATUS)                        :: stat

!! SSN   character(len  =   2)                    :: famtyp

!---------------------------------------------------------------
!  local  variables 
   integer*4                        :: ID_OBS,VARNO,FLAG,VCOORD_TYPE ! 
   real*4                           :: OBSV,OMA,OMP,OER,FGE,PPP
   integer*4                        :: NINSERT ! 
   integer*4                        :: nlv,rln ,id_data,ilast,idf
   character(len  = 128)            :: query
   logical                          :: LLOK

   integer*4                        :: IERR
   integer*4                        :: IDATA,J,JO,JELE
   

!!!   FAMTYP=sqlite_cfamtyp(FILE_NUMB)
!!! SSN   FAMTYP="SFC"
!   NELE_INS=2
!   LISTELE_INS(1)=11215
!   LISTELE_INS(2)=11216
   if (trim(familyType)=='SF') then
     NELE_INS=2
     LISTELE_INS(1)=11215
     LISTELE_INS(2)=11216
   endif
 
   WRITE(*,*)  ' ---INSERTSQL ---   '
   WRITE(*,*)' FAMILY ---> ' , familyType, '  NSTNS  ----> ', obs_numHeader(obsdat)

!---- Print   new elements to be inserted into database
   WRITE( *,*)' INSERT INTO SQLITE FILE ELEMENTS :--> ',(LISTELE_INS(j),j=1,NELE_INS)
!-------------------------------------- 

   !!! SSN SELECT CASE(FAMTYP)
   SELECT CASE(trim(familyType))
     CASE('SF','SC')
      query = ' insert into data (id_obs,varno,vcoord,obsvalue,flag,oma,omp,fg_error,obs_error) VALUES(?,?,?,?,?,?,?,?,?); '
     CASE DEFAULT
      query = ' insert into data (id_obs,varno,vcoord,VCOORD_TYPE ,obsvalue,flag,oma,omp,fg_error,obs_error) VALUES(?,?,?,?,?,?,?,?,?,?); '
   END SELECT

   query=trim(query)

   write(*,*) ' QUERY    = ',TRIM(query)
   write(*,*) ' FAMILLE  = ',trim(familyType)

! Begin  transaction 
!------------------------------------------------------------------------------------
   CALL fSQL_begin(db)
   CALL fSQL_prepare( db, query , stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare : ')
   WRITE(*,*) 'BEGIN  TRANSACTION '
!------------------------------------------------------------------------------------

   NINSERT=0
!
!--------------------------------
!   HEADER LOOP
!--------------------------------
      NINSERT=0
      do JO=1,obs_numHeader(obsdat)
  
        IDF    =obs_headElem_i(obsdat,OBS_IDF,JO)
        if ( IDF .ne. FILE_NUMB) cycle
        ID_OBS =obs_headElem_i(obsdat,OBS_IDO,JO)
        RLN=obs_headElem_i(obsdat,OBS_RLN,JO)
        NLV=obs_headElem_i(obsdat,OBS_NLV,JO)
!--------------------------------
!      DATA LOOP
!--------------------------------
        DO J = RLN, NLV + RLN -1

           ID_DATA    =obs_bodyElem_i(obsdat,OBS_IDD,j)
           VARNO      =obs_bodyElem_i(obsdat,OBS_VNM,j)
           FLAG       =obs_bodyElem_i(obsdat,OBS_FLG,j)
           VCOORD_TYPE=obs_bodyElem_i(obsdat,OBS_VCO,j)
           OBSV       =obs_bodyElem_r(obsdat,OBS_VAR,j)
           OMA        =obs_bodyElem_r(obsdat,OBS_OMA,j)
           OMP        =obs_bodyElem_r(obsdat,OBS_OMP,j)
           OER        =obs_bodyElem_r(obsdat,OBS_OER,j)
           FGE        =obs_bodyElem_r(obsdat,OBS_HPHT,j)
           PPP        =obs_bodyElem_r(obsdat,OBS_PPP,j)
           LLOK  =.FALSE.
           DO JELE=1,NELE_INS
             IF (VARNO .eq. LISTELE_INS(JELE) ) LLOK=.TRUE.
           END DO
       IF ( ID_DATA == -1 ) then
!           IF ( ID_DATA   .gt. 0 ) then
           IF ( LLOK ) THEN
! print *, ' JO id_obs id_data VARNO LLOK  =', JO,ID_OBS,ID_DATA,VARNO,LLOK 
!=========================================================================================
             !!!! SSN SELECT CASE(FAMTYP)
             SELECT CASE(trim(familyType))
              CASE('SF','SC')
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = ID_OBS)
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = VARNO )
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 3, REAL_VAR = PPP   ) 
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 4, REAL_VAR = OBSV  ) 
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 5, INT_VAR  = FLAG  )
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 6, REAL_VAR = OMA   ) 
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 7, REAL_VAR = OMP   ) 
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 8, REAL_VAR = FGE   ) 
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 9, REAL_VAR = OER   ) 
               CALL fSQL_exec_stmt ( stmt)
              CASE DEFAULT
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = ID_OBS )
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = VARNO  )
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 3, REAL_VAR = PPP    ) 
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 4, INT_VAR  = VCOORD_TYPE) 
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 5, REAL_VAR = OBSV   ) 
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 6, INT_VAR  = FLAG   )
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 7, REAL_VAR = OMA    ) 
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 8, REAL_VAR = OMP    ) 
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 9, REAL_VAR = FGE    ) 
               CALL fSQL_bind_param( stmt, PARAM_INDEX = 10,REAL_VAR = OER    ) 
               CALL fSQL_exec_stmt ( stmt)
             END SELECT
!=========================================================================================
             NINSERT=NINSERT + 1
           ENDIF    ! LLOK
        endif
!------------------------------------------------------
      end DO
      end DO


     WRITE(*,*) 'END OF TRANSACTION '

! END OF TRANSACTION 
!------------------------------------------------------------------------------------
     CALL fSQL_finalize( stmt )
     CALL fSQL_commit(db)
!------------------------------------------------------------------------------------
!!! SSN     WRITE(*,*)' FAMILY ---> ' ,FAMTYP, '  NUMBER OF INSERTIONS ----> ',NINSERT
     WRITE(*,*)' FAMILY ---> ' ,trim(familyType), '  NUMBER OF INSERTIONS ----> ',NINSERT

!  ========================
   END SUBROUTINE INSERTSQL
!  ========================

end module sqliteread_mod

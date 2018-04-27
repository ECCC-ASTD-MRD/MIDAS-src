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
use obsUtil_mod

implicit none
!------------------------------------------------------------------------------------------
private            :: handle_error,GEN_DATA,SAVE_INFO_HEADER
public             :: SQL_QUERY_CH,INSERTSQL,SQL2OBS_NML,sqlr_readSqlite

integer            :: NELE_INS,LISTELE_INS(15)
integer            :: N_ITEMS
character*3        :: ITEMLIST(15)
INTEGER, PARAMETER :: N_MAX_FILES=16
integer            :: NFILES_SQL
character(len=2)   :: CFAMTYP_SQL(N_MAX_FILES)
character(len=256) :: CFILNAM_SQL(N_MAX_FILES)

contains
  
   subroutine GEN_DATA(obsdat,VCOORD,OBSVALUE,VARNO,FLAG,VCORDTYP,NOBS,NDATA)
      !----------------------------------------------------------------------------------------------
      !s/r GEN_INFO   - Initialze data values for an observation object
      ! Author     : P. Koclas,  CMC/CMDA, September  2012
      ! Adaptation : S. Skachko, ARMA,     April      2018
      !----------------------------------------------------------------------------------------------
   implicit none
   ! arguments
   type (struct_obs), intent(inout) :: obsdat
   real(obs_real),    intent(in)    :: vcoord,obsvalue
   integer, intent(in)              :: varno,flag,vcordtyp,nobs,ndata
   call obs_bodySet_r(obsdat,OBS_PPP,NDATA,VCOORD)
   call obs_bodySet_r(obsdat,OBS_VAR,NDATA,OBSVALUE)
   call obs_bodySet_i(obsdat,OBS_VNM,NDATA,VARNO)
   call obs_bodySet_i(obsdat,OBS_FLG,NDATA,FLAG)
   call obs_bodySet_i(obsdat,OBS_VCO,NDATA,VCORDTYP)
   end subroutine GEN_DATA
!===================================

  subroutine SAVE_INFO_HEADER( obsdat, rdb4_schema, FamilyType, nobs, elev, id_sat, azimuth, geoid_undulation &
                             , earth_local_rad_curv, ro_qc_flag, instrument, zenith, cloud_cover, solar_zenith &
                             , solar_azimuth, land_sea, id_obs, lat, lon, codtyp,date, time, status, id_stn )
  ! Purpose: aux.routine to save GEN_INFO and GEN_HEADER, routines by P. Koclas, CMC/CMDA
  ! Author: Sergey Skachko, ARMA, April 2018
  implicit none
  ! arguments
  type (struct_obs),  intent(inout) :: obsdat
  character(len=*),   intent(in)    :: rdb4_schema,id_stn
  character*2,        intent(in)    :: FamilyType
  integer,intent(in)                :: nobs,id_obs,codtyp,date,time,status, azimuth, land_sea, ro_qc_flag
  integer,            intent(in)    :: id_sat, instrument, zenith, cloud_cover, solar_zenith, solar_azimuth
  real(obs_real),     intent(in)    :: geoid_undulation, earth_local_rad_curv, elev, lat, lon
  character(len=*), parameter :: my_name = 'SAVE_INFO_HEADER'
  character(len=*), parameter :: my_warning = '****** '// my_name //' WARNING: '
  character(len=*), parameter :: my_error   = '******** '// my_name //' ERROR: '

  call obs_setFamily( obsdat, trim(FamilyType), nobs      )
  call obs_headSet_i( obsdat, OBS_IDO, nobs, id_obs       )       
  call obs_headSet_i( obsdat, OBS_ONM, nobs, nobs         )
  call obs_headSet_i( obsdat, OBS_ITY, nobs, codtyp       )
  call obs_headSet_r( obsdat, OBS_LAT, nobs, lat          )
  call obs_headSet_r( obsdat, OBS_LON, nobs, lon          )
  call obs_headSet_i( obsdat, OBS_DAT, nobs, date         )
  call obs_headSet_i( obsdat, OBS_ETM, nobs, time         )
  call obs_headSet_i( obsdat, OBS_ST1, nobs, status       )
  call     obs_set_c( obsdat, 'STID' , nobs, trim(id_stn) )
  if ( trim(FamilyType) == 'TO' ) then
     call obs_headSet_i( obsdat, OBS_SAT , nobs, id_sat                 )
     call obs_headSet_i( obsdat, OBS_OFL , nobs, land_sea               )
     call obs_headSet_i( obsdat, OBS_INS , nobs, instrument             )
     call obs_headSet_i( obsdat, OBS_SZA , nobs, zenith                 )
     call obs_headSet_i( obsdat, OBS_SUN , nobs, solar_zenith           )
     if ( trim(rdb4_schema) /= 'csr' ) then
        call obs_headSet_i( obsdat, OBS_AZA , nobs, azimuth             )
     endif
     if ( trim(rdb4_schema) =='airs'.or.trim(rdb4_schema) =='iasi'.or.trim(rdb4_schema) =='cris' ) then
        call obs_headSet_i( obsdat, OBS_CLF , nobs, cloud_cover         )
        call obs_headSet_i( obsdat, OBS_SAZ , nobs, solar_azimuth       )  
     elseif ( trim(rdb4_schema) =='amsua'.or.trim(rdb4_schema) =='amsub'.or.trim(rdb4_schema) =='atms' ) then   
        call obs_headSet_i( obsdat, OBS_SAZ , nobs, solar_azimuth       )
     endif
  else
     call obs_headSet_r( obsdat, OBS_ALT , nobs, elev                   )
     if ( trim(rdb4_schema) == 'ro' ) then
       call obs_headSet_i( obsdat, OBS_ROQF, nobs, ro_qc_flag           )
       call obs_headSet_r( obsdat, OBS_GEOI, nobs, geoid_undulation     )
       call obs_headSet_r( obsdat, OBS_TRAD, nobs, earth_local_rad_curv )
       call obs_headSet_i( obsdat, OBS_SAT , nobs, id_sat               )
       call obs_headSet_i( obsdat, OBS_AZA , nobs, azimuth              )
     endif
  endif
  end subroutine SAVE_INFO_HEADER

   subroutine SQL2OBS_NML(NML_SECTION)
      !---------------------------------------------
      !s/r SQL2OBS_NML  - READ NAMELIST FILE SECTION
      !---------------------------------------------
   integer            :: NULNAM,j
   character*256      :: NAMFILE
   character(len = *) :: NML_SECTION
   NAMELIST /NAMSQLinsert/NELE_INS,LISTELE_INS
   NAMELIST /NAMSQLF/     NFILES_SQL,CFILNAM_SQL,CFAMTYP_SQL
   NAMELIST /NAMSQLUPDATE/N_ITEMS,ITEMLIST
   NAMFILE=trim("flnml")
   WRITE(*,*) ' READ NML_SECTION =',trim(NML_SECTION)
   SELECT CASE(trim(NML_SECTION))
     CASE( 'NAMSQLinsert')
       NULNAM=4014
       open(NULNAM,FILE=NAMFILE)
       READ(NULNAM,NML=NAMSQLinsert)
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
   close(NULNAM)
   RETURN
   END SUBROUTINE SQL2OBS_NML
!===============================================================================================

   subroutine sqlr_readSqlite(obsdat,FamilyType,FileName,FileNumb)
      ! Purpose: QUERY DATA FROM SQLITE DATA FILES
      !          THEN INSERT THE DATA INTO  ObsSpaceData structures
      ! Author : S. Skachko, ARMA, April 2018 from the routines SQL2OBS_TOVS and SQL2OBS_CONV by P. Koclas, CMC/CMDA
   implicit none

   ! arguments
   type (struct_obs), intent(inout) :: obsdat     ! ObsSpaceData Structure
   character(len=*),  intent(in)    :: FamilyType ! Family Type may be TOVS or CONV
   character(len=*),  intent(in)    :: FileName   ! SQLITE filename
   integer,           intent(in)    :: FileNumb   ! Stdout file unit number
   ! locals
   integer,parameter        :: nulnam=4000
   character*9              :: rdb4_schema,id_stn
   integer                  :: ID_OBS,ID_DATA,CODTYP,DATE,TIME,STATUS,FLAG,VARNO
   real(obs_real)           :: RELEV,XLAT,XLON,VCOORD
   real                     :: LAT,LON,elevfact,elev
   integer                  :: AZIMUTH,VCORDTYP,VCOORDFACT
   real                     :: RZENITH,RSOLAR_ZENITH,RCLOUD_COVER,RSOLAR_AZIMUTH
   integer                  :: ZENITH,  SOLAR_ZENITH, CLOUD_COVER, SOLAR_AZIMUTH,ro_qc_flag
   real(obs_real)           :: GEOID_UNDULATION, EARTH_LOCAL_RAD_CURV,razimuth,OBSVALUE,SURF_EMISS
   integer                  :: ID_SAT,LAND_SEA,TERRAIN_TYPE,INSTRUMENT,sensor
   integer                  :: I,J,COUNT,nlv,nobs,nobs_start,BITSFLAGON,BITSFLAGOFF,LN,numstns_out,numobs_out
   real(obs_real),parameter :: ZEMFACT=0.01
   character*10             :: CHTIME
   character*512            :: QUERY,QUERY_DAT,QUERY_HDR
   character*256            :: CFG,CSQLCRIT,COLUMNS_HDR,COLUMNS_DAT
   logical                  :: finished,verbose
   real, allocatable        :: matdata(:,:)
   type(fSQL_DATABASE)      :: db         ! type for SQLIte  file handle
   type(fSQL_STATEMENT)     :: stmt,stmt2 ! type for precompiled SQLite statements
   type(fSQL_STATUS)        :: stat,stat2 !type for error status
   character*256            :: LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,namfile
   integer                  :: NBITSOFF,NBITSON,BITOFF(15),BITON(15),NROWS,NCOLUMNS,last_id
   character(len=*), parameter :: my_name = 'sqlr_readSqlite'
   character(len=*), parameter :: my_warning = '****** '// my_name //' WARNING: '
   character(len=*), parameter :: my_error   = '******** '// my_name //' ERROR: '

   NAMELIST /NAMSQLamsua/LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLamsub/LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLairs/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLiasi/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLcris/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLssmi/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLgo/   LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLcsr/  LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLatms/ LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLua/   LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLai/   LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLsw/   LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLro/   LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLsfc/  LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLsc/   LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   NAMELIST /NAMSQLpr/   LISTELEMENTS,SQLEXTRA_DAT,SQLEXTRA_HDR,SQLNULL,SQLLIMIT,NBITSOFF,BITOFF,NBITSON,BITON
   
   write(*,*) 'Subroutine '//my_name
   write(*,*)my_name//': FileName   : ', trim(FileName)
   write(*,*)my_name//': FamilyType : ', trim(FamilyType)
   call fSQL_open( db, trim(FileName) ,stat)
   if ( fSQL_error(stat) /= FSQL_OK ) then
     write(*,*) my_error//'fSQL_open: ', fSQL_errmsg(stat)
     call qqexit(2)
   endif
   CHTIME=SQL_QUERY_CH(db,"select time('now')")
   query="select schema from rdb4_schema ;"
   rdb4_schema = SQL_QUERY_CH(db,trim(query))
   write(*,*) my_name//' START OF QUERY TIME IS = ', CHTIME, 'rdb4_schema is ---> ', trim(rdb4_schema)
   !--- DEFAULT VALUES----------------------------------------------
   SQLEXTRA_HDR=""
   SQLEXTRA_DAT=""
   SQLLIMIT    =""
   SQLNULL     =""
   if ( trim(FamilyType) == 'TO' ) then
     LISTELEMENTS="12163"
   else
     LISTELEMENTS="11001,11002"
   endif
   NBITSON     =0
   NBITSOFF    =0
   vcoordfact  =1
   !----------------------------------------------------------------
   NAMFILE=trim("flnml")
   open(nulnam,FILE=NAMFILE)
   COLUMNS_DAT="id_data,id_obs,vcoord,varno,obsvalue,flag "
   if ( trim(familyType) == 'TO' ) then      
     COLUMNS_HDR="id_obs,lat,lon,codtyp,date,time,status,id_stn,id_sat,land_sea,instrument,zenith,solar_zenith"
     VCORDTYP=3
   else
     COLUMNS_HDR="id_obs,lat,lon,codtyp,date,time,status,id_stn,elev"  
   endif 
   namfile=trim("flnml")
   open(nulnam,file=namfile)
   SELECT CASE(trim(rdb4_schema))
     CASE ('ua')
       SQLNULL=" and obsvalue is not null  "
       VCOORDFACT=1
       VCORDTYP=2
       LISTELEMENTS="12004,12203,10051,10004,11011,11012,11215,11216,12001,11001,11002,11003,11004,12192,10194,13210"
       read(nulnam,NML=NAMSQLua)
     CASE ('ai')
       SQLNULL=" and obsvalue is not null and vcoord is not null "
       VCOORDFACT=1
       VCORDTYP=2
       LISTELEMENTS="12001,11001,11002,12192"
       read(nulnam,NML=NAMSQLai)
     CASE ('sw')
       SQLNULL=" and obsvalue is not null and vcoord is not null "
       VCOORDFACT=1
       VCORDTYP=2
       LISTELEMENTS="11001,11002"
       read(nulnam,NML=NAMSQLsw)
     CASE ('pr')
       SQLNULL=" and obsvalue is not null and vcoord is not null "
       VCOORDFACT=1
       VCORDTYP=1
       LISTELEMENTS="11001,11002"
       read(nulnam,NML=NAMSQLpr)
     CASE ('ro')     
       COLUMNS_HDR=trim(COLUMNS_HDR)//",ro_qc_flag,geoid_undulation,earth_local_rad_curv,id_sat,azimuth"
       SQLNULL=" and obsvalue is not null and vcoord is not null "
       VCOORDFACT=1
       VCORDTYP=1
       LISTELEMENTS="15036,15037"
       read(nulnam,NML=NAMSQLro)
     CASE ('sf')
       SQLNULL=" and obsvalue is not null"
       VCOORDFACT=0
       VCORDTYP=1
       LISTELEMENTS="11011,11012,12004,12203,10051,10004,15031,15032,15035"
       read(nulnam,NML=NAMSQLsfc)
     CASE ('scat')
       SQLNULL=" and obsvalue is not null"
       SQLEXTRA_DAT=" order by id_obs,varno "
       VCOORDFACT=0
       VCORDTYP=1
       LISTELEMENTS="11011,11012"
       read(nulnam,NML=NAMSQLsc)
     CASE( 'airs')
       COLUMNS_HDR=trim(COLUMNS_HDR)//",azimuth,terrain_type,cloud_cover,solar_azimuth"
       COLUMNS_DAT=trim(COLUMNS_DAT)//",surf_emiss"
       read(nulnam,NML=NAMSQLairs)
     CASE( 'iasi')
       COLUMNS_HDR=trim(COLUMNS_HDR)//",azimuth,terrain_type,cloud_cover,solar_azimuth"
       COLUMNS_DAT=trim(COLUMNS_DAT)//",surf_emiss"
       read(nulnam,NML=NAMSQLiasi)
     CASE( 'cris')
       COLUMNS_HDR=trim(COLUMNS_HDR)//",azimuth,terrain_type,cloud_cover,solar_azimuth"
       COLUMNS_DAT=trim(COLUMNS_DAT)//",surf_emiss"
       read(nulnam,NML=NAMSQLcris)
     CASE( 'amsua')
       COLUMNS_HDR=trim(COLUMNS_HDR)//",azimuth,terrain_type,sensor,solar_azimuth"
       read(nulnam,NML=NAMSQLamsua)
     CASE( 'amsub')
       COLUMNS_HDR=trim(COLUMNS_HDR)//",azimuth,terrain_type,sensor,solar_azimuth"
       read(nulnam,NML=NAMSQLamsub)
     CASE( 'atms')
       COLUMNS_HDR=trim(COLUMNS_HDR)//",azimuth,terrain_type,sensor,solar_azimuth"
       read(nulnam,NML=NAMSQLatms)
     CASE( 'ssmi')
       COLUMNS_HDR=trim(COLUMNS_HDR)//",azimuth,terrain_type"
       read(nulnam,NML=NAMSQLssmi)
     CASE( 'csr')
       read(nulnam,NML=NAMSQLcsr)
     CASE DEFAULT
       write(*,*) my_error//' Unsupported  SCHEMA ---> ',trim(rdb4_schema), ' exit '
       return
   END SELECT
   close(nulnam)
   BITSflagoff=0
   do i = 1, nbitsoff
     BITSflagoff = IBSET ( BITSflagoff, 13-BITOFF(I) )
   enddo
   BITSflagon=0
   do i = 1, nbitson
     BITSflagon = IBSET ( BITSflagon,  13-BITON(I)  )
   enddo
   ! COMPOSE SQL  QUERIES 
   write(CFG, '(i6)' ) BITSflagoff
   CSQLcrit=trim("  flag & "//CFG)//" =0 "
   write(CFG, '(i6)' ) BITSflagon
   if  (nbitson > 0) then
     write(CFG, '(i6)' ) BITSflagon
     CSQLcrit=trim(CSQLcrit)//trim("  and flag & "//CFG)
     CSQLcrit=trim(CSQLcrit)//" = "//CFG
   endif
   CSQLcrit=trim(CSQLcrit)//" and varno in ( "//trim(LISTELEMENTS)//" )"
   CSQLcrit=trim(CSQLcrit)//trim(SQLNull)
   CSQLcrit=trim(CSQLcrit)//trim(SQLEXTRA_DAT)
   query_dat=                    "select "//COLUMNS_DAT
   query_dat=trim(query_dat)//trim("   from data where  ")
   query_dat=       trim(query_dat)//"   "
   query_dat=trim(query_dat)//trim(CSQLcrit)
   query_dat=trim(query_dat)//trim(SQLLIMIT)
   query_hdr="select "//COLUMNS_HDR
   query_hdr=trim(query_hdr)//" from header "
   query_hdr=trim(query_hdr)//trim(SQLEXTRA_HDR)
   query_hdr=trim(query_hdr)//" order by id_obs "
   write(*,*) my_name//': ',trim(rdb4_schema),' QUERY_DAT --> ',trim(query_dat)
   write(*,*) my_name//': ',trim(rdb4_schema),' QUERY_HDR --> ',trim(query_hdr)
   write(*,*)' ========================================== '
   if (  trim(rdb4_schema)=='pr' .or. trim(rdb4_schema)=='sf' ) then
     elevfact=1.
   else
     elevfact=0.
   endif
   ! EXTRACT THE DATA FOM THE DATABASE 
   !------------------------------------------------------
   call fSQL_prepare( db, trim(query_hdr) , stmt, stat)
   call fSQL_prepare( db, trim(query_dat), stmt2, stat2)
   if ( fSQL_error(stat)  /= FSQL_OK ) call handle_error(stat ,'fSQL_prepare hdr: ')
   if ( fSQL_error(stat2) /= FSQL_OK ) call handle_error(stat2,'fSQL_prepare dat: ')
   nlv   = 0
   nobs  = obs_numHeader(obsdat)
   count = obs_numBody(obsdat)
   nobs_start = nobs
   write(*,*)  my_name//' DEBUT numheader  =', obs_numHeader(obsdat)
   write(*,*)  my_name//' DEBUT numbody    =', obs_numBody(obsdat)
   
   call fSQL_get_many (  stmt2, nrows = nrows , ncols = ncolumns , mode = FSQL_REAL )
   write(*,*) '  NROWS NCOLUMNS =', nrows, ncolumns, rdb4_schema, trim(query_dat)
   write(*,*)' ========================================== '
   allocate(matdata(nrows,ncolumns))
   call fSQL_fill_matrix ( stmt2, matdata )
   call fSQL_free_mem    ( stmt2 )
   call fSQL_finalize    ( stmt2 )
 
   ! HEADER LOOP
   last_id=1
   HEADER:   DO
     call fSQL_get_row( stmt, finished )   ! Fetch the next row
     
     if (finished) then ! exit LOOP WHEN LAST ROW HAS BEEN FETCHED
        if ( nobs > 1 .and. nlv > 0 ) &
        call SAVE_INFO_HEADER( obsdat, rdb4_schema, FamilyType, nobs, relev, id_sat, azimuth, geoid_undulation &
                             , earth_local_rad_curv, ro_qc_flag, instrument, zenith, cloud_cover, solar_zenith &
                             , solar_azimuth, land_sea, id_obs, xlat, xlon, codtyp, date, time/100, status, id_stn )
        exit
     endif
     ! The query result is inserted into variables
     !=========================================================================================
     terrain_type   = MPC_missingValue_INT
     land_sea       = MPC_missingValue_INT
     sensor         = MPC_missingValue_INT
     cloud_cover    = MPC_missingValue_INT; rcloud_cover   = MPC_missingValue_R4
     elev = 0.; relev = 0.
     solar_azimuth  = MPC_missingValue_INT; rsolar_azimuth = MPC_missingValue_R4
     solar_zenith   = MPC_missingValue_INT; rsolar_zenith  = MPC_missingValue_R4
     zenith         = MPC_missingValue_INT; rzenith        = MPC_missingValue_R4    
     instrument     = MPC_missingValue_INT
     azimuth        = MPC_missingValue_INT; razimuth       =  MPC_missingValue_R8  
     call fSQL_get_column( stmt, COL_INDEX = 1, INT_VAR   = id_obs                     )
     call fSQL_get_column( stmt, COL_INDEX = 2, REAL_VAR  = lat                        )
     call fSQL_get_column( stmt, COL_INDEX = 3, REAL_VAR  = lon                        )
     call fSQL_get_column( stmt, COL_INDEX = 4, INT_VAR   = codtyp                     )
     call fSQL_get_column( stmt, COL_INDEX = 5, INT_VAR   = date                       )
     call fSQL_get_column( stmt, COL_INDEX = 6, INT_VAR   = time                       )
     call fSQL_get_column( stmt, COL_INDEX = 7, INT_VAR   = status                     )
     call fSQL_get_column( stmt, COL_INDEX = 8, CHAR_VAR  = id_stn                     )
     if ( trim(familyType) == 'TO' ) then
         call fSQL_get_column( stmt, COL_INDEX = 9,  INT_VAR   = id_sat                  )
         call fSQL_get_column( stmt, COL_INDEX = 10, INT_VAR   = land_sea , INT_MISSING=MPC_missingValue_INT  )
         call fSQL_get_column( stmt, COL_INDEX = 11, INT_VAR   = instrument,INT_MISSING=MPC_missingValue_INT  )
         call fSQL_get_column( stmt, COL_INDEX = 12, REAL_VAR  = rzenith, REAL_MISSING=MPC_missingValue_R4    )
         call fSQL_get_column( stmt, COL_INDEX = 13, REAL_VAR  = rsolar_zenith, REAL_MISSING=MPC_missingValue_R4  )
         if ( trim(rdb4_schema) /='csr' ) then
            call fSQL_get_column( stmt, COL_INDEX = 14, REAL8_VAR  = razimuth     )
            call fSQL_get_column( stmt, COL_INDEX = 15, INT_VAR   = terrain_type ,INT_MISSING=MPC_missingValue_INT )
         endif
         if ( trim(rdb4_schema)=='airs' .or. trim(rdb4_schema)=='iasi' .or. trim(rdb4_schema)=='cris') then
            call fSQL_get_column( stmt, COL_INDEX = 16, REAL_VAR = rcloud_cover, REAL_MISSING=MPC_missingValue_R4   )
            call fSQL_get_column( stmt, COL_INDEX = 17, REAL_VAR = rsolar_azimuth, REAL_MISSING=MPC_missingValue_R4 )
         elseif ( trim(rdb4_schema)=='amsua' .or. trim(rdb4_schema)=='amsub' .or. trim(rdb4_schema)=='atms') then
            call fSQL_get_column( stmt, COL_INDEX = 16, INT_VAR  = sensor, INT_MISSING=MPC_missingValue_INT )
            call fSQL_get_column( stmt, COL_INDEX = 17, REAL_VAR = rsolar_azimuth, REAL_MISSING=MPC_missingValue_R4 )
         endif
         if ( instrument == 420 ) id_sat  = 784
         zenith = nint ( (90. + rzenith ) * 100. )
         solar_zenith = nint ( (90. + rsolar_zenith ) * 100. )
         azimuth = nint( razimuth * 100. )
         cloud_cover   = nint ( rcloud_cover * 1.     )
         solar_azimuth = nint ( rsolar_azimuth * 100. )
         if ( terrain_type ==  0 ) land_sea = 2  !---Is terrain type sea ice (iterrain=0)?, If so, set imask=2.----
         if ( sensor == MPC_missingValue_INT ) then
            sensor = 0
            if (instrument == MPC_missingValue_INT ) instrument = 0
         else
            instrument = cvt_obs_instrum(sensor)
         endif
     else  ! familyType = CONV
         call fSQL_get_column( stmt, COL_INDEX = 9, REAL_VAR  = elev, REAL_MISSING=MPC_missingValue_R4 )
         relev=elev
         if ( trim(rdb4_schema)=='ro' ) then
            call fSQL_get_column( stmt, COL_INDEX = 10, INT_VAR   = ro_qc_flag, INT_MISSING=MPC_missingValue_INT           )
            call fSQL_get_column( stmt, COL_INDEX = 11, REAL8_VAR = geoid_undulation     )
            call fSQL_get_column( stmt, COL_INDEX = 12, REAL8_VAR = earth_local_rad_curv )
            call fSQL_get_column( stmt, COL_INDEX = 13, INT_VAR   = id_sat, INT_MISSING=MPC_missingValue_INT               )
            call fSQL_get_column( stmt, COL_INDEX = 14, REAL8_VAR = razimuth )
            azimuth = nint( razimuth * 100 )
         endif
     endif  ! TOVS or CONV
     if ( lon < 0.) lon = lon + 360.
     xlat=lat*MPC_RADIANS_PER_DEGREE_R8
     xlon=lon*MPC_RADIANS_PER_DEGREE_R8
     ! INSERT INTO header   OF OBS_SPACE_DATA STRUCTURE
     nobs = nobs + 1 
     nlv=0
     DATA: do j = last_id,nrows
     ! ---------------------------------------------------
     if ( int(matdata(j,2)) > id_obs ) then 
        call SAVE_INFO_HEADER( obsdat, rdb4_schema, FamilyType, nobs, relev, id_sat, azimuth, geoid_undulation &
                             , earth_local_rad_curv, ro_qc_flag, instrument, zenith, cloud_cover, solar_zenith &
                             , solar_azimuth, land_sea, id_obs, xlat, xlon, codtyp,date, time/100, status, id_stn )
        exit
     elseif ( int(matdata(j,2)) == id_obs ) then
        if ( nobs == nobs_start .and. nlv == 0  ) &
        call SAVE_INFO_HEADER( obsdat, rdb4_schema, FamilyType, nobs, relev, id_sat, azimuth, geoid_undulation &
                             , earth_local_rad_curv, ro_qc_flag, instrument, zenith, cloud_cover, solar_zenith &
                             , solar_azimuth, land_sea, id_obs, xlat, xlon, codtyp,date, time/100, status, id_stn )
        last_id=j+1
        id_data=int(matdata(j,1)); flag=int(matdata(j,6));  varno =int(matdata(j,4))
        vcoord=matdata(j,3); obsvalue=matdata(j,5)
        count = count + 1
        nlv   = nlv + 1
        ! INSERT data INTO OBS_SPACE_DATA STRUCTURE
        call obs_bodySet_i(obsdat,OBS_IDD,count,id_data)
        if (trim(rdb4_schema)=='airs'.or.trim(rdb4_schema)=='iasi'.or.trim(rdb4_schema)=='cris') then
           surf_emiss=matdata(j,7)
           call obs_bodySet_r(obsdat,OBS_SEM,count,surf_emiss*zemfact)
       endif
       if ( trim(FamilyType) == 'TO' ) then
          call GEN_DATA(obsdat,vcoord,obsvalue,varno,flag,vcordtyp,nobs,count)
       else ! CONV
          call GEN_DATA(obsdat,vcoord*vcoordfact+relev*elevfact,obsvalue,varno,flag,vcordtyp,nobs,count)
          ! ALLOW EXTRA SPACE FOR U V COMPONENTS
          if ( varno == 11001 .or. varno == 11011) then
             if ( varno == 11001 ) then
                ! U COMPONENT
                call obs_bodySet_i(obsdat,OBS_IDD,count+1,-1)
                call GEN_DATA(obsdat,vcoord*vcoordfact + relev*elevfact,MPC_missingValue_R8,11003,0,vcordtyp,nobs,count+1)
                ! V COMPONENT
                call obs_bodySet_i(obsdat,OBS_IDD,count+2,-1)
                call GEN_DATA(obsdat,vcoord*vcoordfact + relev*elevfact,MPC_missingValue_R8,11004,0,vcordtyp,nobs,count+2)
             else if ( varno == 11011) then
                ! Us COMPONENT
                call obs_bodySet_i(obsdat,OBS_IDD,count+1,-1)
                call GEN_DATA(obsdat,vcoord*vcoordfact + relev*elevfact,MPC_missingValue_R8,11215,0,vcordtyp,nobs,count+1)
                ! Vs COMPONENT
                call obs_bodySet_i(obsdat,OBS_IDD,count+2,-1)
                call GEN_DATA(obsdat,vcoord*vcoordfact + relev*elevfact,MPC_missingValue_R8,11216,0,vcordtyp,nobs,count+2)
             endif
             count = count + 2
             nlv = nlv + 2
          endif      ! extra space for winds
       endif         ! TOVS or CONV
     endif           !  id_obs   
   END DO DATA  ! END OF DATA LOOP

   if ( nobs == 1 ) call obs_headSet_i(obsdat,OBS_RLN,nobs,1)
   if ( nlv > 0 ) then
      call obs_headSet_i(obsdat,OBS_NLV,nobs,nlv)
      if ( nobs > 1 ) then
         LN = obs_headElem_i(obsdat,OBS_RLN,nobs-1) +  obs_headElem_i(obsdat,OBS_NLV,nobs-1)
         call obs_headSet_i(obsdat,OBS_RLN,nobs,LN)
      endif   
      if ( last_id > nrows) &
      call SAVE_INFO_HEADER( obsdat, rdb4_schema, FamilyType, nobs, relev, id_sat, azimuth, geoid_undulation &
                        , earth_local_rad_curv, ro_qc_flag, instrument, zenith, cloud_cover, solar_zenith &
                        , solar_azimuth, land_sea, id_obs, xlat, xlon, codtyp,date, time/100, status, id_stn )
   else
      nobs=nobs-1
   endif
  END DO HEADER ! HEADER
  deallocate(matdata)
  write(*,*)  my_name//' FIN numheader  =', obs_numHeader(obsdat)                    
  write(*,*)  my_name//' FIN numbody    =', obs_numBody(obsdat)
  write(*,*)  my_name//' fin header '
  CHTIME=SQL_QUERY_CH(db,"select time('now')")
  write(*,*) my_name//' END OF QUERY TIME IS = ', CHTIME,rdb4_schema
  call fSQL_finalize( stmt )
  call fSQL_close( db , stat) ! Close The SQLITE FILE
  write(*,*) my_name//' end subroutine: ', rdb4_schema
  end subroutine sqlr_readSqlite
!===================================================================

   function SQL_QUERY_CH(db,query)
      !   Purpose : RETURN RESULT OF A QUERY  TO AND OPENED SQLITE FILE
      !             THE RESULT FROM THIS QURY MUS BE A CHARATER STRING
      ! Author  : P. Koclas, CMC/CMDA September  2012
      !    ARGUMENTS:
      !                 INPUT:
      !                       -db      : SQLITE FILE HANDLE
      !                       -query   : A QUERY
      !    ---------------------------------------------------
   implicit none
   ! arguments
   type(fSQL_DATABASE)                      :: db   ! type handle for  SQLIte file
   CHARACTER(len = *)                       :: query
   ! locals
   CHARACTER(len = 256 )                          :: SQL_QUERY_CH
   CHARACTER(len = 256 )                          :: CH_RESULT    
   logical finished
   type(fSQL_STATEMENT)                     :: stmt !  prepared statement for  SQLite
   type(fSQL_STATUS)                        :: stat !type error status
   CH_RESULT=''
   CALL fSQL_prepare( db, trim(query), stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare: ')
   finished=.FALSE.
   CALL fSQL_get_row( stmt, finished )
   !  Put result of query into variable
   CALL fSQL_get_column( stmt, COL_INDEX = 1, CHAR_VAR   = CH_RESULT    )
   CALL fSQL_get_row( stmt, finished )
   if ( .NOT. finished) then
      write(*,*)' SQL_QUERY_CH ---> QUERY RETURNS MORE THAN ONE ROW...  '
    endif
   CALL fSQL_finalize( stmt )
   SQL_QUERY_CH = trim(CH_RESULT)
   RETURN
   end function   SQL_QUERY_CH
!  ==============================

   subroutine handle_error(stat, message)
      !
      ! Exit code on f90sqlite 
      ! error message
      ! Author  : P. Koclas CMDA/CMC SEPTEMBER 2012
      !
   type(FSQL_STATUS)  :: stat
   character(len = *) :: message
   write(*,*) message, fSQL_errmsg(stat)
   call qqexit(2)
   END SUBROUTINE handle_error
!  ================================================================================

   subroutine updsql(db,obsdat,familyType,fileName,file_numb)
      !
      !   Purpose : UPDATE HEADER AND DATA TABLES OF SQLITE FILES
      !
      ! Author    : P. Koclas,  CMC/CMDA CMDA/CMC SEPTEMBER 2012
      ! Adaptation: S. Skachko, ARMA, April 2018
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
   implicit none
   ! arguments
   type(fSQL_DATABASE)              :: db
   type (struct_obs), intent(inout) :: obsdat
   character(len=*)                 :: fileName
   character(len=*)                 :: familyType
   integer                          :: file_numb
   type(fSQL_STATEMENT)             :: stmt ! prepared statement for  SQLite
   type(fSQL_STATUS)                :: stat ! type error status
   !  locals 
   integer                          :: RLN,NLV,IDF,ID_DATA,ASS,FLAG ! 
   integer                          :: IOBS,ID_OBS,STATUS ! 
   integer                          :: J,JO !
   integer                          :: last_comma ! 
   integer                          :: ibegin,ilast,ibeginob,ilastob,NUPD ! 
   character*10                     :: CHTIME
   character(len  =   3)            :: item
   integer                          :: upd_list(20),UPDLIST,I,K
   character(len  =   9)            :: item2
   character(len  = 128)            :: query
   character(len  = 356)            :: CHitem,CHitem2
   logical                          :: BACK
   real                             :: ROMP,OBSV

   write(*,*) '  ======================================================================= '
   write(*,*) '  ==============BEGIN UPDATE============================================== '
   write(*,*) '  ======================================================================= '
   write(*,*) '  ------- UPDSQL ------ ',itemlist(1)
   write(*,*) '  FAMILYTYPE   =', FAMILYTYPE
   write(*,*) '  FileName     =', FileName
   write(*,*) '  MISSING VALUE  = ', MPC_missingValue_R8
   !    CREATE QUERY  
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
         WRITE(*,*)'invalid item: ', item2
         WRITE(*,*)' EXIT ROUTINE: '
         RETURN
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
   write(*,*) ' NUMBER OF HEADERS in obsdat    =',obs_numHeader(obsdat)
   call fSQL_begin(db)
   !   HEADER LOOP
   NUPD=0
   DO JO=1,obs_numHeader(obsdat)

      IDF=obs_headElem_i(obsdat,OBS_IDF,JO)
      if ( IDF .ne. FILE_NUMB) cycle
      ID_OBS=obs_headElem_i(obsdat,OBS_IDO,JO)
      RLN   =obs_headElem_i(obsdat,OBS_RLN,JO)
      NLV   =obs_headElem_i(obsdat,OBS_NLV,JO)
      ! TOP DATA LOOP
      DO J = RLN, NLV + RLN -1

        FLAG   =obs_bodyElem_i(obsdat,OBS_FLG,j)
        ID_DATA=obs_bodyElem_i(obsdat,OBS_IDD,j)
        ASS    =obs_bodyElem_i(obsdat,OBS_ASS,j)
        if ( ASS == 1 ) then
          call fSQL_bind_param(stmt, PARAM_INDEX = 1,   INT_VAR  = flag  )
          DO K =1,UPDLIST  
            OBSV =obs_bodyElem_r(obsdat,OBS_VAR,j)
            if ( OBSV /= MPC_missingValue_R8 ) then  
              ROMP= obs_bodyElem_r(obsdat,UPD_LIST(K),j)
              CALL fSQL_bind_param(stmt, PARAM_INDEX = K+1, REAL_VAR = ROMP   )
            endif
          END DO             
          CALL fSQL_bind_param(stmt, PARAM_INDEX = UPDLIST+2, INT_VAR  = ID_DATA  )
          CALL fSQL_exec_stmt (stmt)
          NUPD=NUPD + 1
        endif
      END DO
   END DO
   CALL fSQL_finalize( stmt )
   WRITE(*,*) ' NUMBER OF DATA ROWS UPDATES =',NUPD
   ! UPDATES FOR THE STATUS FLAGS IN THE HEADER TABBLE
   NUPD=0
   query = ' UPDATE header SET status  = ?   WHERE ID_OBS = ? '
   CALL fSQL_prepare( db, query , stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare : ')
   DO JO = 1,obs_numHeader(obsdat)
      IDF=obs_headElem_i(obsdat,OBS_IDF,JO)
      if ( IDF /= FILE_NUMB) cycle
      ID_OBS=obs_headElem_i(obsdat,OBS_IDO,JO)
      STATUS=obs_headElem_i(obsdat,OBS_ST1,JO)
      CALL fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = STATUS )
      CALL fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = ID_OBS )
      CALL fSQL_exec_stmt ( stmt)
      NUPD=NUPD + 1
   END DO
   WRITE(*,*) ' NUMBER OF HEADER ROW UPDATES =',NUPD
   CALL fSQL_finalize( stmt )
   CALL fSQL_commit(db)
   write(*,*) '  ===============END UPDATE============================================== '
   END SUBROUTINE updsql
!  ==========================================================================================

   subroutine insertsql(db,obsdat,FamilyType,FileName,file_numb)
      !   Purpose : INSERT NEW DATA ROWS INTO SQLITE FILES
      !
      ! Author    : P. Koclas , CMC/CMDA, September 2012
      ! Adaptation: S. Skachko, ARMA    , April     2018
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
   implicit none
   ! arguments
   type(fSQL_DATABASE)    :: db   ! type for SQLIte  file handle
   type(struct_obs)       :: obsdat
   character(len=*)       :: familyType
   character(len=*)       :: fileName
   integer                :: file_numb
   type(fSQL_STATEMENT)   :: stmt ! type for precompiled SQLite statements
   type(fSQL_STATUS)      :: stat !type for error status
   !  locals 
   integer                :: ID_OBS,VARNO,FLAG,VCOORD_TYPE ! 
   real                   :: OBSV,OMA,OMP,OER,FGE,PPP
   integer                :: NINSERT,IDATA,J,JO,JELE, IERR
   integer                :: nlv,rln ,id_data,ilast,idf
   character(len  = 128)  :: query
   logical                :: LLOK

   write(*,*)  ' ---INSERTSQL ---   '
   write(*,*)' FAMILY ---> ', trim(FamilyType), '  NSTNS  ----> ', obs_numHeader(obsdat)
   write(*,*)' FileName -> ', trim(FileName)
   write(*,*)' INSERT INTO SQLITE FILE ELEMENTS :--> ',(LISTELE_INS(j),j=1,NELE_INS)
   SELECT CASE(trim(familyType))
     CASE('SF','SC','GP')
      query = ' insert into data (id_obs,varno,vcoord,obsvalue,flag,oma,omp,fg_error,obs_error) VALUES(?,?,?,?,?,?,?,?,?); '
     CASE DEFAULT
      query = ' insert into data (id_obs,varno,vcoord,VCOORD_TYPE ,obsvalue,flag,oma,omp,fg_error,obs_error) VALUES(?,?,?,?,?,?,?,?,?,?); '
   END SELECT
   query=trim(query)
   write(*,*) ' QUERY    = ',TRIM(query)
   write(*,*) ' FAMILLE  = ',trim(familyType)
   CALL fSQL_begin(db)
   CALL fSQL_prepare( db, query , stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare : ')
   WRITE(*,*) 'BEGIN  TRANSACTION '
   NINSERT=0
   ! HEADER LOOP
   NINSERT=0
   do JO=1,obs_numHeader(obsdat)
     IDF    =obs_headElem_i(obsdat,OBS_IDF,JO)
     if ( IDF /= FILE_NUMB) cycle
     ID_OBS =obs_headElem_i(obsdat,OBS_IDO,JO)
     RLN=obs_headElem_i(obsdat,OBS_RLN,JO)
     NLV=obs_headElem_i(obsdat,OBS_NLV,JO)
     ! DATA LOOP
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
         IF (VARNO == LISTELE_INS(JELE) ) LLOK=.TRUE.
       END DO
       IF ( ID_DATA == -1 ) then
         IF ( LLOK ) THEN
           SELECT CASE(trim(familyType))
           CASE('SF','SC','GP')
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
           NINSERT=NINSERT + 1
         ENDIF    ! LLOK
      endif
    end DO
    end DO
    WRITE(*,*) 'END OF TRANSACTION '
    CALL fSQL_finalize( stmt )
    CALL fSQL_commit(db)
     write(*,*)' FAMILY ---> ' ,trim(FamilyType), '  NUMBER OF INSERTIONS ----> ',NINSERT
    END SUBROUTINE INSERTSQL
!  ========================
end module sqliteread_mod

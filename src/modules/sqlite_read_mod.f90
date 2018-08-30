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
!! MODULE sqliteRead (prefix= "sqlr")
!!
!! *Purpose*: To read and update SQLITE observation files. Data is stored in 
!!            obsSpaceData object.
!!
!--------------------------------------------------------------------------
module sqliteRead_mod

use codePrecision_mod
use obsSpaceData_mod
use mpi_mod
use fSQLite
use mathPhysConstants_mod
use mpiVar_mod
use obsUtil_mod
use utilities_mod
use bufr_mod

implicit none
save

private
public :: sqlr_insertSqlite, sqlr_updateSqlite, sqlr_readSqlite, sqlr_query
public :: sqlr_thinSqlite

contains
  
  subroutine sqlr_initData(obsdat, vertCoord, obsValue, obsVarno, obsFlag, vertCoordType, numberData )
    ! Initialze data values for an observation object
    implicit none
    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    real(obs_real)   , intent(in)    :: vertCoord, obsValue
    integer          , intent(in)    :: obsVarno, obsFlag, vertCoordType, numberData

    call obs_bodySet_r(obsdat, OBS_PPP, numberData, vertCoord     )
    call obs_bodySet_r(obsdat, OBS_VAR, numberData, obsValue      )
    call obs_bodySet_i(obsdat, OBS_VNM, numberData, obsVarno      )
    call obs_bodySet_i(obsdat, OBS_FLG, numberData, obsFlag       )
    call obs_bodySet_i(obsdat, OBS_VCO, numberData, vertCoordType )

  end subroutine sqlr_initData


  subroutine sqlr_initHeader( obsdat, rdbSchema, familyType, headerIndex, elev, obsSat, azimuth, geoidUndulation, &
                              earthLocRadCurv, roQcFlag, instrument, zenith, cloudCover, solarZenith, &
                              solarAzimuth, landSea, obsIdo, obsLat, obsLon, codeType, obsDate, obsTime, obsStatus, idStation )
    implicit none
    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    character(len=*) , intent(in)    :: rdbSchema,idStation
    character(len=2) , intent(in)    :: familyType
    integer          , intent(in)    :: headerIndex, obsIdo, codeType, obsDate, obsTime, obsStatus, azimuth, landSea, roQcFlag
    integer          , intent(in)    :: obsSat, instrument, zenith, cloudCover, solarZenith, solarAzimuth
    real(obs_real)   , intent(in)    :: geoidUndulation, earthLocRadCurv, elev, obsLat, obsLon

    call obs_setFamily( obsdat, trim(familyType), headerIndex       )
    call obs_headSet_i( obsdat, OBS_IDO, headerIndex, obsIdo        )       
    call obs_headSet_i( obsdat, OBS_ONM, headerIndex, headerIndex   )
    call obs_headSet_i( obsdat, OBS_ITY, headerIndex, codeType      )
    call obs_headSet_r( obsdat, OBS_LAT, headerIndex, obsLat        )
    call obs_headSet_r( obsdat, OBS_LON, headerIndex, obsLon        )
    call obs_headSet_i( obsdat, OBS_DAT, headerIndex, obsDate       )
    call obs_headSet_i( obsdat, OBS_ETM, headerIndex, obsTime       )
    call obs_headSet_i( obsdat, OBS_ST1, headerIndex, obsStatus     )
    call     obs_set_c( obsdat, 'STID' , headerIndex, trim(idStation) )

    if ( trim(familyType) == 'TO' ) then

      call obs_headSet_i( obsdat, OBS_SAT , headerIndex, obsSat      )
      call obs_headSet_i( obsdat, OBS_OFL , headerIndex, landSea     )
      call obs_headSet_i( obsdat, OBS_INS , headerIndex, instrument  )
      call obs_headSet_i( obsdat, OBS_SZA , headerIndex, zenith      )
      call obs_headSet_i( obsdat, OBS_SUN , headerIndex, solarZenith )

      if ( trim(rdbSchema) /= 'csr' ) then
        call obs_headSet_i( obsdat, OBS_AZA , headerIndex, azimuth )
      end if

      if ( trim(rdbSchema) == 'airs' .or. trim(rdbSchema) == 'iasi' .or. trim(rdbSchema) == 'cris' ) then
        call obs_headSet_i( obsdat, OBS_CLF , headerIndex, cloudCover   )
        call obs_headSet_i( obsdat, OBS_SAZ , headerIndex, solarAzimuth )  
      else if ( trim(rdbSchema) == 'amsua' .or. trim(rdbSchema) == 'amsub' .or. trim(rdbSchema) == 'atms' ) then   
        call obs_headSet_i( obsdat, OBS_SAZ , headerIndex, solarAzimuth )
      end if

    else

      call obs_headSet_r( obsdat, OBS_ALT ,headerIndex , elev )
      if ( trim(rdbSchema) == 'ro' ) then
        call obs_headSet_i( obsdat, OBS_ROQF, headerIndex, roQcFlag        )
        call obs_headSet_r( obsdat, OBS_GEOI, headerIndex, geoidUndulation )
        call obs_headSet_r( obsdat, OBS_TRAD, headerIndex, earthLocRadCurv )
        call obs_headSet_i( obsdat, OBS_SAT , headerIndex, obsSat          )
        call obs_headSet_i( obsdat, OBS_AZA , headerIndex, azimuth         )
      end if

    end if

  end subroutine sqlr_initHeader


  subroutine sqlr_readSqlite(obsdat, familyType, fileName, FileNumb )
    implicit none
    ! arguments
    type (struct_obs), intent(inout) :: obsdat     ! ObsSpaceData Structure
    character(len=*) , intent(in)    :: familyType ! Family Type may be TOVS or CONV
    character(len=*) , intent(in)    :: fileName   ! SQLITE filename
    integer          , intent(in)    :: FileNumb   ! Stdout file unit number
    ! locals
    character(len=9)         :: rdbSchema,idStation
    integer                  :: obsIdo, obsIdd, codeType, obsDate, obsTime, obsStatus, obsFlag, obsVarno
    real(obs_real)           :: elevReal, xlat, xlon, vertCoord
    real                     :: elev    , obsLat, obsLon, elevFact
    integer                  :: azimuth, vertCoordType, vertCoordFact, fnom, fclos, nulnam, ierr 
    real                     :: zenithReal, solarZenithReal, CloudCoverReal, solarAzimuthReal
    integer                  :: zenith    , solarZenith    , CloudCover    , solarAzimuth    , roQcFlag
    real(obs_real)           :: geoidUndulation, earthLocRadCurv, azimuthReal, obsValue, surfEmiss
    integer                  :: obsSat, landSea, terrainType, instrument, sensor, numberElem
    integer                  :: i, rowIndex, obsNlv, headerIndex, headerIndexStart, bodyIndex, bitsFlagOn, bitsFlagOff, reportLocation
    real(obs_real),parameter :: zemFact = 0.01
    character(len=10)        :: timeCharacter
    character(len=512)       :: query, queryData, queryHeader
    character(len=256)       :: cfgSqlite, csqlcrit, columnsHeader, columnsData
    logical                  :: finished
    real, allocatable        :: matdata(:,:)
    type(fSQL_DATABASE)      :: db         ! type for SQLIte  file handle
    type(fSQL_STATEMENT)     :: stmt,stmt2 ! type for precompiled SQLite statements
    type(fSQL_STATUS)        :: stat,stat2 !type for error status
    character(len=256)       :: listElem, sqlExtraDat, sqlExtraHeader, sqlNull, sqlLimit, namelistFileName
    integer                  :: numberBitsOff, numberBitsOn, bitsOff(15), bitsOn(15), numberRows, numberColumns, lastId
    character(len=*), parameter :: myName = 'sqlr_readSqlite'
    character(len=*), parameter :: myWarning = '****** '// myName //' WARNING: '
    character(len=*), parameter :: myError   = '******** '// myName //' ERROR: '
    namelist /NAMSQLamsua/numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLamsub/numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLairs/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLiasi/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLcris/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLssmi/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLgo/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLcsr/  numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLatms/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLua/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLai/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLsw/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLro/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLsfc/  numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLsc/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLpr/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLal/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    
    write(*,*) 'Subroutine '//myName
    write(*,*) myName//': fileName   : ', trim(fileName)
    write(*,*) myName//': familyType : ', trim(familyType)
    call fSQL_open( db, trim(fileName) ,stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) myError//'fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( myError//': fSQL_open' )
    end if

    timeCharacter = sqlr_query(db,"select time('now')")
    query = "select schema from rdb4_schema ;"
    rdbSchema = sqlr_query(db,trim(query))
    write(*,'(4a)') myName//' START OF QUERY TIME IS = ', timeCharacter, 'rdbSchema is ---> ', trim(rdbSchema)

    sqlExtraHeader = ''
    sqlExtraDat    = ''
    sqlLimit       = ''
    sqlNull        = ''
    listElem       = ''
    numberElem = 0
    numberBitsOff = 0
    bitsOff =  0
    numberBitsOn = 0
    bitsOn = 0

    if ( trim(familyType) == 'TO' ) then
      write(listElem,*) bufr_nbt3
    else
      write(listElem,'(i5.5,",",i5.5)') bufr_nedd, bufr_neff
    end if

    vertCoordfact  = 1
    columnsData = " id_data, id_obs, vcoord, varno, obsvalue, flag "
    if ( trim(familyType) == 'TO' ) then      
      columnsHeader = " id_obs, lat, lon, codtyp, date, time, status, id_stn, id_sat, land_sea, instrument, zenith, solar_zenith "
      vertCoordType = 3
    else
      columnsHeader = " id_obs, lat, lon, codtyp, date, time, status, id_stn, elev"  
    end if
    
    nulnam = 0
    ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0 )
    select case(trim(rdbSchema))
      case ( 'ua' )
        vertCoordFact = 1
        vertCoordType = 2
        read(nulnam, nml = NAMSQLua, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLua )
      case ( 'ai' )
        vertCoordFact = 1
        vertCoordType = 2
        read(nulnam, nml = NAMSQLai, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLai )
      case ( 'sw' )
        vertCoordFact = 1
        vertCoordType = 2
        read(nulnam, nml = NAMSQLsw, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLsw )
      case ( 'pr' )
        vertCoordFact = 1
        vertCoordType = 1
        read(nulnam, nml = NAMSQLpr, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLpr )
      case ( 'al' )
        vertCoordFact = 1
        vertCoordType = 1
        read(nulnam, nml = NAMSQLal, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLal )
      case ( 'ro' )     
        columnsHeader = trim(columnsHeader)//",ro_qc_flag, geoid_undulation, earth_local_rad_curv, id_sat, azimuth"
        vertCoordFact = 1
        vertCoordType = 1
        read(nulnam, nml = NAMSQLro, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLro )
      case ( 'sf' )
        vertCoordFact = 0
        vertCoordType = 1
        read(nulnam, nml = NAMSQLsfc, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLsfc )
      case ( 'scat' )
        vertCoordFact = 0
        vertCoordType = 1
        read(nulnam, nml = NAMSQLsc, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLsc )
      case( 'airs' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, cloud_cover, solar_azimuth "
        columnsData = trim(columnsData)//", surf_emiss "
        read(nulnam, nml = NAMSQLairs, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLairs )
      case( 'iasi' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, cloud_cover, solar_azimuth "
        columnsData = trim(columnsData)//", surf_emiss "
        read(nulnam, nml = NAMSQLiasi, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLiasi )
      case( 'cris' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, cloud_cover, solar_azimuth "
        columnsData = trim(columnsData)//", surf_emiss "
        read(nulnam, nml = NAMSQLcris, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLcris )
      case( 'amsua' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, sensor, solar_azimuth "
        read(nulnam, nml = NAMSQLamsua, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLamsua )
      case( 'amsub' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, sensor, solar_azimuth "
        read(nulnam, nml = NAMSQLamsub, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLamsub )
      case( 'atms')
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, sensor, solar_azimuth "
        read(nulnam, nml = NAMSQLatms, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLatms )
      case( 'ssmi' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type "
        read(nulnam, nml = NAMSQLssmi, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLssmi )
      case( 'csr' )
        read(nulnam, nml = NAMSQLcsr, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLcsr )
      case DEFAULT
        write(*,*) myError//' Unsupported  SCHEMA ---> ',trim(rdbSchema), ' ABORT!!! '
        call utl_abort( myError//': Unsupported  SCHEMA in SQLITE file!' )
    end select
    ierr=fclos( nulnam )

    ! Compose SQL queries
    bitsFlagOff=0
    do i = 1, numberBitsOff
      bitsFlagOff = ibset ( bitsFlagOff, 13 - bitsOff(i) )
    end do
    bitsFlagOn=0
    do i = 1, numberBitsOn
      bitsFlagOn = ibset ( bitsFlagOn, 13 - bitsOn(i) )
    end do

    write(cfgSqlite, '(i6)' ) bitsFlagOff
    csqlcrit=trim("  flag & "//cfgSqlite)//" = 0 "
    write(cfgSqlite, '(i6)' ) bitsFlagOn

    if ( numberBitsOn > 0 ) then
      write(cfgSqlite, '(i6)' ) bitsFlagOn
      csqlcrit = trim(CSQLcrit)//trim(" and flag & "//cfgSqlite)
      csqlcrit = trim(CSQLcrit)//" = "//cfgSqlite
    end if

    csqlcrit = trim(CSQLcrit)//" and varno in ( "//trim(listElem)//" )"//trim(SQLNull)//trim(sqlExtraDat)
    queryData= "select "//columnsData
    queryData = trim(queryData)//trim(" from data where ")//trim(csqlcrit)//trim(sqlLimit)//";"
    queryHeader="select "//trim(columnsHeader)//" from header "//trim(sqlExtraHeader)//" order by id_obs;"
    write(*,'(4a)') myName//': ',trim(rdbSchema),' queryData    --> ', trim(queryData)
    write(*,'(4a)') myName//': ',trim(rdbSchema),' queryHeader --> ', trim(queryHeader)
    write(*,*)' ========================================== '

    if ( trim(rdbSchema)=='pr' .or. trim(rdbSchema)=='sf' ) then
      elevFact=1.
    else
      elevFact=0.
    end if

    call fSQL_prepare( db, trim(queryHeader) , stmt, stat)
    call fSQL_prepare( db, trim(queryData), stmt2, stat2)
    if ( fSQL_error(stat)  /= FSQL_OK ) call sqlr_handleError(stat ,'fSQL_prepare hdr: ')
    if ( fSQL_error(stat2) /= FSQL_OK ) call sqlr_handleError(stat2,'fSQL_prepare dat: ')

    headerIndex  = obs_numHeader(obsdat)
    bodyIndex = obs_numBody(obsdat)
    headerIndexStart = headerIndex
    write(*,*) myName//' DEBUT numheader  =', obs_numHeader(obsdat)
    write(*,*) myName//' DEBUT numbody    =', obs_numBody(obsdat)
    call fSQL_get_many (  stmt2, nrows = numberRows , ncols = numberColumns , mode = FSQL_REAL )
    write(*,*) myName//'  numberRows numberColumns =', numberRows, numberColumns, rdbSchema, trim(queryData)
    write(*,*)' ========================================== '
    allocate( matdata(numberRows, numberColumns) )
    call fSQL_fill_matrix ( stmt2, matdata )
    call fSQL_free_mem    ( stmt2 )
    call fSQL_finalize    ( stmt2 )

    obsNlv = 0
    lastId = 1

    HEADER: do

      call fSQL_get_row( stmt, finished )   ! Fetch the next row
      if (finished) then                    ! exit LOOP WHEN LAST ROW HAS BEEN FETCHED
        if ( headerIndex > 1 .and. obsNlv > 0 ) &
          call sqlr_initHeader( obsdat, rdbSchema, familyType, headerIndex, elevReal, obsSat, azimuth, geoidUndulation, &
                                earthLocRadCurv, roQcFlag, instrument, zenith, cloudCover, solarZenith, &
                                solarAzimuth, landSea, obsIdo, xlat, xlon, codeType, obsDate, obsTime/100, obsStatus, idStation )
        exit HEADER
      end if

      ! The query result is inserted into variables
      terrainType = MPC_missingValue_INT
      landSea     = MPC_missingValue_INT
      sensor      = MPC_missingValue_INT
      cloudCover  = MPC_missingValue_INT; cloudCoverReal = MPC_missingValue_R4
      elev = 0.; elevReal = 0.
      solarAzimuth  = MPC_missingValue_INT; solarAzimuthReal = MPC_missingValue_R4
      solarZenith   = MPC_missingValue_INT; solarZenithReal  = MPC_missingValue_R4
      zenith = MPC_missingValue_INT; zenithReal = MPC_missingValue_R4    
      instrument = MPC_missingValue_INT
      azimuth = MPC_missingValue_INT; azimuthReal = MPC_missingValue_R8  

      call fSQL_get_column( stmt, COL_INDEX = 1, INT_VAR   = obsIdo     )
      call fSQL_get_column( stmt, COL_INDEX = 2, REAL_VAR  = obsLat       )
      call fSQL_get_column( stmt, COL_INDEX = 3, REAL_VAR  = obsLon       )
      call fSQL_get_column( stmt, COL_INDEX = 4, INT_VAR   = codeType  )
      call fSQL_get_column( stmt, COL_INDEX = 5, INT_VAR   = obsDate      )
      call fSQL_get_column( stmt, COL_INDEX = 6, INT_VAR   = obsTime      )
      call fSQL_get_column( stmt, COL_INDEX = 7, INT_VAR   = obsStatus    )
      call fSQL_get_column( stmt, COL_INDEX = 8, CHAR_VAR  = idStation )

      if ( trim(familyType) == 'TO' ) then

        call fSQL_get_column( stmt, COL_INDEX = 9,  INT_VAR   = obsSat )
        call fSQL_get_column( stmt, COL_INDEX = 10, INT_VAR   = landSea ,   INT_MISSING=MPC_missingValue_INT  )
        call fSQL_get_column( stmt, COL_INDEX = 11, INT_VAR   = instrument, INT_MISSING=MPC_missingValue_INT  )
        call fSQL_get_column( stmt, COL_INDEX = 12, REAL_VAR  = zenithReal, REAL_MISSING=MPC_missingValue_R4  )
        call fSQL_get_column( stmt, COL_INDEX = 13, REAL_VAR  = solarZenithReal, REAL_MISSING=MPC_missingValue_R4 )

        if ( trim(rdbSchema) /= 'csr' ) then

          call fSQL_get_column( stmt, COL_INDEX = 14, REAL8_VAR  = azimuthReal )
          call fSQL_get_column( stmt, COL_INDEX = 15, INT_VAR   = terrainType, INT_MISSING=MPC_missingValue_INT )

        end if

        if ( trim(rdbSchema) == 'airs' .or. trim(rdbSchema) == 'iasi' .or. trim(rdbSchema) == 'cris' ) then

          call fSQL_get_column( stmt, COL_INDEX = 16, REAL_VAR = cloudCoverReal, REAL_MISSING=MPC_missingValue_R4   )
          call fSQL_get_column( stmt, COL_INDEX = 17, REAL_VAR = solarAzimuthReal, REAL_MISSING=MPC_missingValue_R4 )

        else if ( trim(rdbSchema) == 'amsua' .or. trim(rdbSchema) == 'amsub' .or. trim(rdbSchema) == 'atms' ) then

          call fSQL_get_column( stmt, COL_INDEX = 16, INT_VAR  = sensor, INT_MISSING=MPC_missingValue_INT )
          call fSQL_get_column( stmt, COL_INDEX = 17, REAL_VAR = solarAzimuthReal, REAL_MISSING=MPC_missingValue_R4 )

        end if

        if ( instrument == 420 ) obsSat  = 784
        zenith = nint ( (90. + zenithReal ) * 100. )
        solarZenith = nint ( (90. + solarZenithReal ) * 100. )
        azimuth = nint( azimuthReal * 100. )
        cloudCover   = nint ( cloudCoverReal * 1.     )
        solarAzimuth = nint ( solarAzimuthReal * 100. )
        if ( terrainType ==  0 ) landSea = 2  !---Is terrain type sea ice (iterrain=0)?, If so, set imask=2.----
        if ( sensor == MPC_missingValue_INT ) then
          sensor = 0
          if (instrument == MPC_missingValue_INT ) instrument = 0
        else
          instrument = obsu_cvt_obs_instrum(sensor)
        end if

      else  ! familyType = CONV

        call fSQL_get_column( stmt, COL_INDEX = 9, REAL_VAR  = elev, REAL_MISSING=MPC_missingValue_R4 )
        elevReal=elev

        if ( trim(rdbSchema)=='ro' ) then

          call fSQL_get_column( stmt, COL_INDEX = 10, INT_VAR   = roQcFlag, INT_MISSING=MPC_missingValue_INT )
          call fSQL_get_column( stmt, COL_INDEX = 11, REAL8_VAR = geoidUndulation )
          call fSQL_get_column( stmt, COL_INDEX = 12, REAL8_VAR = earthLocRadCurv )
          call fSQL_get_column( stmt, COL_INDEX = 13, INT_VAR   = obsSat, INT_MISSING=MPC_missingValue_INT )
          call fSQL_get_column( stmt, COL_INDEX = 14, REAL8_VAR = azimuthReal )
          azimuth = nint( azimuthReal * 100 )

        end if

      end if  ! TOVS or CONV

      if ( obsLon < 0. ) obsLon = obsLon + 360.
      xlat = obsLat * MPC_RADIANS_PER_DEGREE_R8
      xlon = obsLon * MPC_RADIANS_PER_DEGREE_R8
      headerIndex = headerIndex + 1 
      
      obsNlv = 0
      DATA: do rowIndex = lastId,numberRows
        ! ---------------------------------------------------
        if ( int(matdata(rowIndex,2)) > obsIdo ) then 

          call sqlr_initHeader( obsdat, rdbSchema, familyType, headerIndex, elevReal, obsSat, azimuth, geoidUndulation, &
                                earthLocRadCurv, roQcFlag, instrument, zenith, cloudCover, solarZenith, &
                                solarAzimuth, landSea, obsIdo, xlat, xlon, codeType, obsDate, obsTime/100, obsStatus, idStation )
          exit DATA

        else if ( int(matdata(rowIndex,2)) == obsIdo ) then

          if ( headerIndex == headerIndexStart .and. obsNlv == 0  ) &
            call sqlr_initHeader( obsdat, rdbSchema, familyType, headerIndex, elevReal, obsSat, azimuth, geoidUndulation, &
                                  earthLocRadCurv, roQcFlag, instrument, zenith, cloudCover, solarZenith, &
                                  solarAzimuth, landSea, obsIdo, xlat, xlon, codeType, obsDate, obsTime/100, obsStatus, idStation )

          lastId = rowIndex + 1
          obsIdd = int(matdata(rowIndex,1))
          vertCoord = matdata(rowIndex,3)
          obsVarno = int(matdata(rowIndex,4))
          obsValue = matdata(rowIndex,5)
          obsFlag = int(matdata(rowIndex,6))

          bodyIndex = bodyIndex + 1
          obsNlv   = obsNlv + 1
          call obs_bodySet_i(obsdat, OBS_IDD, bodyIndex, obsIdd)

          if (trim(rdbSchema) == 'airs' .or. trim(rdbSchema) == 'iasi' .or. trim(rdbSchema) == 'cris' ) then

            surfEmiss = matdata(rowIndex,7)
            call obs_bodySet_r(obsdat, OBS_SEM, bodyIndex, surfEmiss * zemFact)

          end if

          if ( trim(familyType) == 'TO' ) then

            call sqlr_initData(obsdat, vertCoord, obsValue, obsVarno, obsFlag, vertCoordType, bodyIndex)

          else ! CONV
 
           call sqlr_initData( obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                                obsValue, obsVarno, obsFlag,vertCoordType, bodyIndex)

            if ( obsVarno == bufr_nedd .or. obsVarno == bufr_neds ) then ! ALLOW EXTRA SPACE FOR U V COMPONENTS

              if ( obsVarno == bufr_nedd ) then

                ! U COMPONENT
                call obs_bodySet_i( obsdat, OBS_IDD, bodyIndex + 1, -1)
                call sqlr_initData( obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                                    MPC_missingValue_R8, bufr_neuu, 0, vertCoordType, bodyIndex+1)
                ! V COMPONENT
                call obs_bodySet_i(obsdat, OBS_IDD, bodyIndex + 2, -1)
                call sqlr_initData( obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                                    MPC_missingValue_R8, bufr_nevv, 0, vertCoordType, bodyIndex+2)

              else if ( obsVarno == bufr_neds) then

                ! Us COMPONENT
                call obs_bodySet_i(obsdat, OBS_IDD, bodyIndex + 1, -1)
                call sqlr_initData( obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                                    MPC_missingValue_R8, bufr_neus, 0, vertCoordType, bodyIndex+1)
                ! Vs COMPONENT
                call obs_bodySet_i(obsdat, OBS_IDD, bodyIndex + 2, -1)
                call sqlr_initData( obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                                    MPC_missingValue_R8, bufr_nevs, 0, vertCoordType, bodyIndex+2)
              end if

              bodyIndex = bodyIndex + 2
              obsNlv = obsNlv + 2

            end if      ! extra space for winds

          end if       ! TOVS or CONV

        end if          !  obsIdo   

      end do DATA  ! END OF DATA LOOP

      if ( headerIndex == 1 ) call obs_headSet_i(obsdat, OBS_RLN, headerIndex, 1 )

      if ( obsNlv > 0 ) then

        call obs_headSet_i(obsdat, OBS_NLV, headerIndex, obsNlv )

        if ( headerIndex > 1 ) then
          reportLocation = obs_headElem_i(obsdat, OBS_RLN, headerIndex - 1 ) +  obs_headElem_i(obsdat, OBS_NLV, headerIndex - 1 )
          call obs_headSet_i(obsdat, OBS_RLN, headerIndex, reportLocation )
        end if   

        if ( lastId > numberRows ) &
          call sqlr_initHeader( obsdat, rdbSchema, familyType, headerIndex, elevReal, obsSat, azimuth, geoidUndulation, &
                                earthLocRadCurv, roQcFlag, instrument, zenith, cloudCover, solarZenith, &
                                solarAzimuth, landSea, obsIdo, xlat, xlon, codeType,obsDate, obsTime/100, obsStatus, idStation )

      else

        headerIndex = headerIndex - 1

      end if

    end do HEADER ! HEADER

    deallocate(matdata)
    write(*,*)  myName//' FIN numheader  =', obs_numHeader(obsdat)                    
    write(*,*)  myName//' FIN numbody    =', obs_numBody(obsdat)
    write(*,*)  myName//' fin header '
    timeCharacter = sqlr_query(db,"select time('now')")
    write(*,*) myName//' END OF QUERY TIME IS = ', timeCharacter,rdbSchema
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 
    write(*,*) myName//' end subroutine: ', rdbSchema

  end subroutine sqlr_readSqlite


  function sqlr_query(db,query)
    implicit none
    ! arguments
    type(fSQL_DATABASE)  :: db   ! type handle for  SQLIte file
    character(len = *)   :: query
    ! locals
    character(len = 256) :: sqlr_query, result
    logical finished
    type(fSQL_STATEMENT) :: stmt !  prepared statement for  SQLite
    type(fSQL_STATUS)    :: stat !type error status

    result=''
    call fSQL_prepare( db, trim(query), stmt, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat,'fSQL_prepare: ')
    finished=.false.
    call fSQL_get_row( stmt, finished )

    ! Put result of query into variable
    call fSQL_get_column( stmt, COL_INDEX = 1, CHAR_VAR = result )
    call fSQL_get_row( stmt, finished )
    if ( .not. finished ) write(*,*)' SQL QUERY ---> QUERY RETURNS MORE THAN ONE ROW...  '
    call fSQL_finalize( stmt )
    sqlr_query = trim(result)

  end function sqlr_query


  subroutine sqlr_handleError(stat, message)
    implicit none

    type(FSQL_STATUS)  :: stat
    character(len = *) :: message

    write(*,*) message, fSQL_errmsg(stat)
    call utl_abort( trim(message) )

  end subroutine sqlr_handleError

  
  subroutine sqlr_updateSqlite(db, obsdat, familyType, fileName, fileNumber)
    implicit none
    ! arguments
    type(fSQL_DATABASE),intent(inout):: db
    type (struct_obs),  intent(inout) :: obsdat
    character(len=*),   intent(in)    :: fileName   
    character(len=*),   intent(in)    :: familyType ! Observation Family Type 
    integer         ,   intent(in)    :: fileNumber ! FILE NUMBER ASSOCIATED WITH db

    !  locals
    type(fSQL_STATEMENT)             :: stmt ! prepared statement for  SQLite
    type(fSQL_STATUS)                :: stat ! type error status
    integer                          :: obsRln, obsNlv, obsIdf, obsIdd, obsAss, obsFlag, iobs, obsIdo, obsStatus, last_comma
    integer                          :: ibegin, ilast, ibeginob, ilastob, headerIndex, bodyIndex, numberUpdateItems
    character(len =  10)             :: timeCharacter
    character(len =   3)             :: item, itemUpdateList(15)
    integer                          :: itemId, updateList(20), fnom, fclos, nulnam, ierr
    character(len =   9)             :: item2
    character(len = 128)             :: query
    character(len = 356)             :: itemChar,item2Char
    logical                          :: back
    integer                          :: ivalue
    real                             :: romp, obsValue
    character(len=*), parameter      :: myName = 'sqlr_updateSqlite'
    character(len=*), parameter      :: myWarning = '****** '// myName //' WARNING: '
    character(len=*), parameter      :: myError   = '******** '// myName //' ERROR: '
    namelist/namSQLUpdate/ numberUpdateItems, itemUpdateList

    write(*,*) myName//' Starting ===================  '

    ! set default values of namelist variables
    itemUpdateList(:) = ''
    numberUpdateItems = 0

    ! Read the namelist for directives
    nulnam = 0
    ierr   = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0 )
    read(nulnam,nml = namSQLUpdate, iostat = ierr )
    if ( ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
    if ( mpi_myid == 0 ) write(*, nml = namSQLUpdate )
    ierr = fclos( nulnam )

    write(*,'(3a,i8)') myName//': Family Type   = ', trim(familyType), ' number of items to update: ', numberUpdateItems
    write(*,*)         myName//': File Name     = ', trim(fileName)
    write(*,*)         myName//': Missing Value = ', MPC_missingValue_R8    

    ! CREATE QUERY  
    itemChar='  '
    if ( numberUpdateItems == 0) then
      write(*,*) ' NO UPDATES TO DO : ---- RETURN'
      return
    end if

    do itemId = 1, numberUpdateItems
      item = itemUpdateList(itemId)
      write(*,*)'Updating ', itemId, item
      select case(item)
        case('OMA')
          updateList(itemId) = OBS_OMA
          item2='oma'
        case('OMP')
          updateList(itemId) = OBS_OMP
          item2='omp'
        case('VAR')
          updateList(itemId) = OBS_VAR
          item2='obsvalue'
        case('OER')
          updateList(itemId) = OBS_OER
          item2='obs_error'
        case('FGE')
          updateList(itemId) = OBS_HPHT
          item2='fg_error'
        case('FLG')
          updateList(itemId) = OBS_FLG
          item2='flag'
        case DEFAULT
          write(*,*)'invalid item: ', item2,' EXIT sqlr_updateSQL!!!'
          call utl_abort( myError//': invalid item ' )
      end select
      itemChar = trim(itemChar)//trim(item2)//trim(' = ? ,')
    end do

    back=.true.
    last_comma  = scan(itemChar, ',', back)
    item2Char   = itemChar(1:last_comma-1)
    itemChar    = item2Char
    query = ' update data set flag = ? , '//trim(itemChar); query=trim(query)//' where id_data = ?  ;'
    write(*,*) ' Update query --->  ', query
    call fSQL_prepare( db, query , stmt, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat, 'fSQL_prepare : ')
    call fSQL_begin(db)

    HEADER: do headerIndex = 1, obs_numHeader(obsdat)
 
      obsIdf = obs_headElem_i( obsdat,OBS_IDF, headerIndex )
 
      if ( obsIdf /= fileNumber) cycle HEADER
      obsIdo = obs_headElem_i( obsdat, OBS_IDO, headerIndex )
      obsRln = obs_headElem_i( obsdat, OBS_RLN, headerIndex )
      obsNlv = obs_headElem_i( obsdat, OBS_NLV, headerIndex )

      BODY: do bodyIndex = obsRln, obsNlv + obsRln - 1

        obsFlag = obs_bodyElem_i( obsdat, OBS_FLG, bodyIndex )
        obsIdd  = obs_bodyElem_i( obsdat, OBS_IDD, bodyIndex )
        obsAss  = obs_bodyElem_i( obsdat, OBS_ASS, bodyIndex )

        if ( obsAss == 1 ) then
          call fSQL_bind_param(stmt, PARAM_INDEX = 1,   INT_VAR  = obsFlag  )
          ITEMS: do itemId = 1, numberUpdateItems

            obsValue = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)

            if(updateList(itemId) == OBS_FLG)then
              ivalue = obs_bodyElem_i(obsdat, updateList(itemId), bodyIndex )
              call fSQL_bind_param(stmt, PARAM_INDEX = itemId + 1, &
                                   INT_VAR = ivalue )
              
            else ! OBS_FLG
              if ( obsValue /= MPC_missingValue_R8 ) then  
                romp = obs_bodyElem_r(obsdat, updateList(itemId), bodyIndex )
                call fSQL_bind_param(stmt, PARAM_INDEX = itemId + 1, &
                                     REAL_VAR = romp )
              end if
            end if ! OBS_FLG

          end do ITEMS            
          
          call fSQL_bind_param(stmt, PARAM_INDEX = numberUpdateItems + 2, INT_VAR  = obsIdd )
          call fSQL_exec_stmt (stmt)

        end if ! obsAss == 1

      end do BODY

    end do HEADER
    
    call fSQL_finalize( stmt )
    ! UPDATES FOR THE STATUS FLAGS IN THE HEADER TABBLE
    query = ' update header set status  = ? where id_obs = ? '
    call fSQL_prepare( db, query , stmt, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat,'fSQL_prepare : ')
    
    HEADER2: do headerIndex = 1,obs_numHeader(obsdat)

      obsIdf = obs_headElem_i(obsdat, OBS_IDF, headerIndex )
      if ( obsIdf /= fileNumber ) cycle HEADER2
      obsIdo    = obs_headElem_i(obsdat, OBS_IDO, headerIndex )
      obsStatus = obs_headElem_i(obsdat, OBS_ST1, headerIndex )
      call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = obsStatus )
      call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = obsIdo )
      call fSQL_exec_stmt ( stmt)

    end do HEADER2
    
    call fSQL_finalize( stmt )
    call fSQL_commit(db)
    write(*,*) myName//' End ===================  ', trim(familyType)

  end subroutine sqlr_updateSqlite


  subroutine sqlr_insertSqlite( db, obsdat, familyType, fileName, fileNumber )
    implicit none
    ! arguments
    type(fSQL_DATABASE)    :: db   ! type for SQLIte  file handle
    type(struct_obs)       :: obsdat
    character(len=*)       :: familyType
    character(len=*)       :: fileName
    integer                :: fileNumber
    ! locals
    integer                :: itemInsertList(15), numberInsertItems
    type(fSQL_STATEMENT)   :: stmt ! type for precompiled SQLite statements
    type(fSQL_STATUS)      :: stat !type for error status
    integer                :: obsVarno, obsFlag, vertCoordType, fnom, fclos, nulnam, ierr 
    real                   :: obsValue, OMA, OMP, OER, FGE, PPP
    integer                :: numberInsert, idata, headerIndex, bodyIndex, obsNlv, obsRln, obsIdd, obsIdo, ilast, obsIdf, insertItem
    character(len  = 256)  :: query
    logical                :: llok    
    character(len=*), parameter :: myName = 'sqlr_insertSqlite'
    character(len=*), parameter :: myWarning = '****** '// myName //' WARNING: '
    character(len=*), parameter :: myError   = '******** '// myName //' ERROR: '
    namelist/namSQLInsert/ numberInsertItems, itemInsertList

    write(*,*)  myName//' --- Starting ---   '
    write(*,*)' FAMILY ---> ', trim(familyType), '  headerIndex  ----> ', obs_numHeader(obsdat)
    write(*,*)' fileName -> ', trim(fileName)   
   
    ! set default values of namelist variables
    numberInsertItems = 0
    itemInsertList(:) = ''

    nulnam = 0
    ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0 )
    read(nulnam, nml = namSQLInsert, iostat = ierr )
    if (ierr /= 0 ) call utl_abort( myError//': Error reading namelist' )
    if (mpi_myid == 0) write(*, nml = namSQLInsert )
    ierr=fclos( nulnam )

    write(*,*)' INSERT INTO SQLITE FILE ELEMENTS :--> ',( itemInsertList(insertItem), insertItem = 1, numberInsertItems )

    select case( trim( familyType ) )

      case( 'SF', 'SC', 'GP' )

        query = 'insert into data (id_obs,varno,vcoord,obsvalue,flag,oma,omp,fg_error,obs_error) values(?,?,?,?,?,?,?,?,?);'

      case DEFAULT

        query = 'insert into data (id_obs,varno,vcoord,vcoord_type,obsvalue,flag,oma,omp,fg_error,obs_error) values(?,?,?,?,?,?,?,?,?,?);'

    end select

    query=trim(query)
    write(*,*) ' === Family Type === ',trim(familyType)
    write(*,*) ' Insert query = ', trim(query)

    call fSQL_begin(db)
    call fSQL_prepare( db, query , stmt, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat, 'fSQL_prepare : ')

    numberInsert=0
    HEADER: do headerIndex = 1, obs_numHeader(obsdat)

      obsIdf = obs_headElem_i(obsdat, OBS_IDF, headerIndex )
      if ( obsIdf /= fileNumber ) cycle HEADER
      obsIdo = obs_headElem_i(obsdat, OBS_IDO, headerIndex )
      obsRln = obs_headElem_i(obsdat, OBS_RLN, headerIndex )
      obsNlv = obs_headElem_i(obsdat, OBS_NLV, headerIndex )

      BODY: do bodyIndex = obsRln, obsNlv + obsRln -1

        obsIdd        = obs_bodyElem_i(obsdat, OBS_IDD, bodyIndex )
        obsVarno      = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex )
        obsFlag       = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex )
        vertCoordType = obs_bodyElem_i(obsdat, OBS_VCO, bodyIndex )
        obsValue      = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex )
        OMA           = obs_bodyElem_r(obsdat, OBS_OMA, bodyIndex )
        OMP           = obs_bodyElem_r(obsdat, OBS_OMP, bodyIndex )
        OER           = obs_bodyElem_r(obsdat, OBS_OER, bodyIndex )
        FGE           = obs_bodyElem_r(obsdat, OBS_HPHT, bodyIndex )
        PPP           = obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex )
        llok = .false.

        do insertItem = 1, numberInsertItems
          if ( obsVarno == itemInsertList(insertItem) ) llok=.true.
        end do

        if ( obsIdd == -1 ) then
          if ( llok ) then
            select case(trim(familyType))
              case( 'SF', 'SC', 'GP' )
                call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = obsIdo   )
                call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = obsVarno )
                call fSQL_bind_param( stmt, PARAM_INDEX = 3, REAL_VAR = PPP      ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 4, REAL_VAR = obsValue ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 5, INT_VAR  = obsFlag  )
                call fSQL_bind_param( stmt, PARAM_INDEX = 6, REAL_VAR = OMA      ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 7, REAL_VAR = OMP      ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 8, REAL_VAR = FGE      ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 9, REAL_VAR = OER      ) 
                call fSQL_exec_stmt ( stmt)
              case DEFAULT
                call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = obsIdo        )
                call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = obsVarno      )
                call fSQL_bind_param( stmt, PARAM_INDEX = 3, REAL_VAR = PPP           ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 4, INT_VAR  = vertCoordType ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 5, REAL_VAR = obsValue      ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 6, INT_VAR  = obsFlag       )
                call fSQL_bind_param( stmt, PARAM_INDEX = 7, REAL_VAR = OMA           ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 8, REAL_VAR = OMP           ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 9, REAL_VAR = FGE           ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 10,REAL_VAR = OER           ) 
                call fSQL_exec_stmt ( stmt)
            end select
            numberInsert = numberInsert + 1
          end if    ! llok
        end if

      end do BODY

    end do HEADER

    call fSQL_finalize( stmt )
    call fSQL_commit(db)
    write(*,'(3a,i8)') myName//' FAMILY ---> ' ,trim(familyType), '  NUMBER OF INSERTIONS ----> ', numberInsert

  end subroutine sqlr_insertSqlite


  !--------------------------------------------------------------------------
  !!
  !! *Purpose*: to flagged (bit 11 set) observations in an SQL file
  !!
  !--------------------------------------------------------------------------
  subroutine sqlr_thinSqlite(db, obsdat, familyType, fileName, fileNumber)
    implicit none
    ! arguments
    type(fSQL_DATABASE), intent(inout) :: db   ! SQLite file handle
    type(struct_obs),    intent(inout) :: obsdat
    character(len=*),    intent(in) :: familyType
    character(len=*),    intent(in) :: fileName
    integer,             intent(in) :: fileNumber

    ! locals
    character(len=*), parameter :: myName = 'sqlr_thinSqlite'
    character(len=*), parameter :: myWarning = '****** '// myName //' WARNING: '
    character(len=*), parameter :: myError   = '******** '// myName //' ERROR: '

    character(len = 128) :: query
    type(fSQL_STATEMENT) :: statement ! prepared statement for SQLite
    type(fSQL_STATUS)    :: status

    call fSQL_open(db, fileName, status)
    if (fSQL_error(status) /= FSQL_OK) then
      write(*,*) myError, fSQL_errmsg(status)
    end if

    ! Mark for deletion all records with bit 11 (0x800) set
    query = ' delete from data where flag & 800;'
    call fSQL_prepare(db, query, statement, status)
    if (fSQL_error(status) /= FSQL_OK) &
      call sqlr_handleError(status, 'thinning fSQL_prepare : ')
    call fSQL_begin(db)
    call fSQL_exec_stmt(statement)
    call fSQL_finalize(statement)
    call fSQL_commit(db)

    write(*,*)'  closed database -->', trim(FileName)
    call fSQL_close( db, status )
  end subroutine sqlr_thinSqlite

end module sqliteRead_mod

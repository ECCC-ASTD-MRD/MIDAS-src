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

module sqliteRead_mod
  ! MODULE sqliteRead (prefix='sqlr' category='6. Observation input/output')
  !
  ! :Purpose: To read and update SQLITE observation files. Data is stored in 
  !           obsSpaceData object.
  !

use codePrecision_mod
use obsSpaceData_mod
use mpi_mod
use fSQLite
use mathPhysConstants_mod
use mpiVar_mod
use obsUtil_mod
use utilities_mod
use bufr_mod
use ramDisk_mod
use tovs_nl_mod
use codtyp_mod
use obsVariableTransforms_mod
use obsFilter_mod
use sqliteUtilities_mod

implicit none   
 
save

private

public :: sqlr_insertSqlite, sqlr_updateSqlite, sqlr_readSqlite, sqlr_query
public :: sqlr_cleanSqlite, sqlr_writeAllSqlDiagFiles, sqlr_readSqlite_avhrr, sqlr_addCloudParametersandEmissivity

contains
  
  subroutine sqlr_initData(obsdat, vertCoord, obsValue, obsVarno, obsFlag, vertCoordType, numberData )
    !
    ! :Purpose: Initialize data values for an observation object.
    !
    implicit none

    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    real(pre_obsReal), intent(in)    :: vertCoord
    real(pre_obsReal), intent(in)    :: obsValue
    integer          , intent(in)    :: obsVarno
    integer          , intent(in)    :: obsFlag 
    integer          , intent(in)    :: vertCoordType
    integer          , intent(in)    :: numberData

    call obs_bodySet_r(obsdat, OBS_PPP, numberData, vertCoord     )
    call obs_bodySet_r(obsdat, OBS_VAR, numberData, obsValue      )
    call obs_bodySet_i(obsdat, OBS_VNM, numberData, obsVarno      )
    call obs_bodySet_i(obsdat, OBS_FLG, numberData, obsFlag       )
    call obs_bodySet_i(obsdat, OBS_VCO, numberData, vertCoordType )

  end subroutine sqlr_initData


  subroutine sqlr_initHeader( obsdat, rdbSchema, familyType, headerIndex, elev, obsSat, azimuth, geoidUndulation, &
                              earthLocRadCurv, roQcFlag, instrument, zenith, cloudCover, solarZenith, &
                              solarAzimuth, terrainType, landSea, iasiImagerCollocationFlag, iasiGeneralQualityFlag,    &
                              headPrimaryKey, obsLat, obsLon, codeType, obsDate, obsTime, &
                              obsStatus, idStation, idProf, trackCellNum, modelWindSpeed, &
                              obsrzam,  obsrele, obsrans, obsrane,                        &
                              firstBodyIndexOfThisBatch, obsNlv)
    !
    ! :Purpose: To initialize the header information when SQLite files are read.
    !
    implicit none

    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    character(len=*) , intent(in)    :: rdbSchema
    character(len=*) , intent(in)    :: idStation
    character(len=2) , intent(in)    :: familyType
    integer          , intent(in)    :: headerIndex
    integer(8)       , intent(in)    :: headPrimaryKey
    integer          , intent(in)    :: codeType
    integer          , intent(in)    :: obsDate
    integer          , intent(in)    :: obsTime
    integer          , intent(in)    :: obsStatus
    integer          , intent(in)    :: terrainType
    integer          , intent(in)    :: landSea
    integer          , intent(in)    :: iasiImagerCollocationFlag
    integer          , intent(in)    :: iasiGeneralQualityFlag
    integer          , intent(in)    :: roQcFlag
    integer          , intent(in)    :: obsSat
    integer          , intent(in)    :: instrument
    integer          , intent(in)    :: idProf
    integer          , intent(in)    :: trackCellNum
    integer          , intent(in)    :: firstBodyIndexOfThisBatch
    integer          , intent(in)    :: obsNlv
    real(pre_obsReal), intent(in)    :: geoidUndulation
    real(pre_obsReal), intent(in)    :: earthLocRadCurv
    real(pre_obsReal), intent(in)    :: elev
    real(pre_obsReal), intent(in)    :: obsLat
    real(pre_obsReal), intent(in)    :: obsLon
    real(pre_obsReal), intent(in)    :: solarAzimuth
    real(pre_obsReal), intent(in)    :: cloudCover
    real(pre_obsReal), intent(in)    :: solarZenith
    real(pre_obsReal), intent(in)    :: zenith
    real(pre_obsReal), intent(in)    :: azimuth
    real(pre_obsReal), intent(in)    :: modelWindSpeed
    real(pre_obsReal), intent(in)    :: obsrzam
    real(pre_obsReal), intent(in)    :: obsrele
    real(pre_obsReal), intent(in)    :: obsrans
    real(pre_obsReal), intent(in)    :: obsrane

    call obs_setFamily( obsdat, trim(familyType), headerIndex       )
    call obs_setHeadPrimaryKey( obsdat,  headerIndex, headPrimaryKey)
    call obs_headSet_i( obsdat, OBS_ONM, headerIndex, headerIndex   )
    call obs_headSet_i( obsdat, OBS_ITY, headerIndex, codeType      )
    call obs_headSet_r( obsdat, OBS_LAT, headerIndex, obsLat        )
    call obs_headSet_r( obsdat, OBS_LON, headerIndex, obsLon        )
    call obs_headSet_i( obsdat, OBS_DAT, headerIndex, obsDate       )
    call obs_headSet_i( obsdat, OBS_ETM, headerIndex, obsTime       )
    call obs_headSet_i( obsdat, OBS_ST1, headerIndex, obsStatus     )
    call     obs_set_c( obsdat, 'STID' , headerIndex, trim(idStation) )
    call obs_headSet_i( obsdat, OBS_RLN, headerIndex, firstBodyIndexOfThisBatch )
    call obs_headSet_i( obsdat, OBS_NLV, headerIndex, obsNlv )

    if ( trim(familyType) == 'TO' ) then

      call obs_headSet_i( obsdat, OBS_SAT , headerIndex, obsSat      )
      call obs_headSet_i( obsdat, OBS_TTYP, headerIndex, terrainType     )
      call obs_headSet_i( obsdat, OBS_STYP, headerIndex, landSea     )
      call obs_headSet_i( obsdat, OBS_INS , headerIndex, instrument  )
      call obs_headSet_r( obsdat, OBS_SZA , headerIndex, zenith      )
      call obs_headSet_r( obsdat, OBS_SUN , headerIndex, solarZenith )
      if ( obs_columnActive_IH(obsdat,OBS_GQF) ) call obs_headSet_i(obsdat,OBS_GQF,headerIndex,iasiGeneralQualityFlag)
      if ( obs_columnActive_IH(obsdat,OBS_GQL) ) call obs_headSet_i(obsdat,OBS_GQL,headerIndex,iasiImagerCollocationFlag)

      if ( trim(rdbSchema) /= 'csr' ) then
        call obs_headSet_r( obsdat, OBS_AZA , headerIndex, azimuth )
      end if

      if ( trim(rdbSchema) == 'airs' .or. trim(rdbSchema) == 'iasi' .or. trim(rdbSchema) == 'cris' ) then
        call obs_headSet_r( obsdat, OBS_CLF , headerIndex, cloudCover   )
        call obs_headSet_r( obsdat, OBS_SAZ , headerIndex, solarAzimuth )  
      else if ( trim(rdbSchema) == 'amsua' .or. trim(rdbSchema) == 'amsub' .or. trim(rdbSchema) == 'atms' ) then   
        call obs_headSet_r( obsdat, OBS_SAZ , headerIndex, solarAzimuth )
      end if

    else

      call obs_headSet_r( obsdat, OBS_ALT ,headerIndex , elev )
      
      if ( trim(rdbSchema) == 'sst' ) then
        call obs_headSet_r( obsdat, OBS_SUN , headerIndex, solarZenith )
      end if
      
      if ( trim(rdbSchema) == 'ro' ) then
        call obs_headSet_i( obsdat, OBS_ROQF, headerIndex, roQcFlag        )
        call obs_headSet_r( obsdat, OBS_GEOI, headerIndex, geoidUndulation )
        call obs_headSet_r( obsdat, OBS_TRAD, headerIndex, earthLocRadCurv )
        call obs_headSet_i( obsdat, OBS_SAT , headerIndex, obsSat          )
        call obs_headSet_r( obsdat, OBS_AZA , headerIndex, azimuth         )
      else if ( trim(rdbSchema) == 'al' ) then
        call obs_headSet_i( obsdat, OBS_PRFL, headerIndex, idProf          )
      else if ( trim(rdbSchema) == 'gl_ascat' ) then
         call obs_headSet_i( obsdat, OBS_FOV, headerIndex, trackCellNum   )
         call obs_headSet_r( obsdat, OBS_MWS, headerIndex, modelWindSpeed    )
      end if
      if ( trim(rdbSchema) == 'radvel' ) then
          call obs_headSet_r( obsdat, OBS_RZAM, headerIndex, obsrzam      )
          call obs_headSet_r( obsdat, OBS_RELE, headerIndex, obsrele      )
          call obs_headSet_r( obsdat, OBS_RANS, headerIndex, obsrans      )
          call obs_headSet_r( obsdat, OBS_RANE, headerIndex, obsrane      )
      end if   
    end if

  end subroutine sqlr_initHeader


  function sqlr_doesSQLTableExist(db, tableName) result(doesSQLTableExist)
    ! :Purpose: check if the table 'tableName' exists in the SQLite database 'db'
    !    The database 'db' should already have been opened with a call to 'fSQL_open'
    !
    implicit none
    ! arguments:
    type(fSQL_DATABASE)      :: db
    character(len=*)         :: tableName
    ! output
    logical                  :: doesSQLTableExist
    ! locals
    character(len=128)       :: querySqlite, sqliteOutput

    querySqlite = "SELECT count(*) FROM sqlite_master WHERE type='table' AND name like '"// tableName //"' ;"
    sqliteOutput = sqlr_query( db, trim( querySqlite ) )
    if ( trim(sqliteOutput) == '1' ) then
      doesSQLTableExist = .true.
    else
      doesSQLTableExist = .false.
    end if

  end function sqlr_doesSQLTableExist


  subroutine sqlr_readSqlite_avhrr(obsdat, fileName, headerIndexBegin, headerIndexEnd )
    ! :Purpose: To read SQLite avhrr_cloud parameters .
    ! 
    implicit none
    ! Arguments:
    type (struct_obs), intent(inout) :: obsdat     ! ObsSpaceData Structure
    character(len=*) , intent(in)    :: fileName   ! SQLite filename
    integer          , intent(in)    :: headerIndexBegin, headerIndexEnd
    ! locals
    type(fSQL_DATABASE)      :: db   ! type for SQLIte  file handle
    type(fSQL_STATEMENT)     :: stmt ! type for precompiled SQLite statements
    type(fSQL_STATUS)        :: stat ! type for error status
    integer                  :: obsIdo
    character(len=128)       :: querySqlite
    integer                  :: rowIndex, headerIndex, columnIndex
    integer                  :: numberRows ,  numberColumns
    real, allocatable        :: matdata(:,:)
    REAL(pre_obsReal) :: CFRAC,MOYRAD,STDRAD
    character(len=*), parameter :: myName = 'sqlr_readSqlite_avhrr'
    character(len=*), parameter :: myWarning = '****** '// myName //': WARNING: '

    write(*,*) myName//': fileName : ', trim(fileName)

    call fSQL_open( db, trim(fileName) ,stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myName//': fSQL_open '//fSQL_errmsg(stat) )
    end if

    if ( sqlr_doesSQLTableExist(db,'avhrr') ) then
      write(*,*) myName//': Table avhrr does not exist :  ... return  '
      return
    end if
    write(*,*) myName//': Table avhrr exists: insert contents into obsdat '

    querySqlite = ' select mean_radiance,stddev_radiance,fractionClearPixels from avhrr where id_obs = ? '
    call fSQL_prepare( db, querySqlite , stmt, stat )
    write(*,*) myName//': obs_getNchanAvhr=',obs_getNchanAvhrr()
    do headerIndex = headerIndexBegin, headerIndexEnd
      obsIdo = obs_headPrimaryKey( obsdat, headerIndex )
      call fSQL_bind_param(stmt, PARAM_INDEX = 1, INT_VAR  =  obsIdo )
      call fSQL_exec_stmt (stmt)
      call fSQL_get_many (  stmt, nrows = numberRows , ncols = numberColumns , mode = FSQL_REAL )
      allocate( matdata(numberRows, numberColumns) )
      matdata(:,:) = MPC_missingValue_R4
      call fSQL_fill_matrix ( stmt, matdata )

      rowIndex=1
      do columnIndex=OBS_CF1,OBS_CF7
        if(obs_columnActive_RH(obsdat,columnIndex)) then
          CFRAC = matdata(rowIndex,3)
          call obs_headSet_r(obsdat,columnIndex,headerIndex, CFRAC)
          rowIndex = rowIndex+6
        end if
      end do

      rowIndex=1
      do columnIndex=OBS_M1C1,OBS_M7C6
        if(obs_columnActive_RH(obsdat,columnIndex)) then
          MOYRAD = matdata(rowIndex,1) * 100000.d0
          call obs_headSet_r(obsdat,columnIndex,headerIndex, MOYRAD )
          rowIndex = rowIndex+1
        endif
      end do

      rowIndex=1
      do columnIndex=OBS_S1C1,OBS_S7C6
        if(obs_columnActive_RH(obsdat,columnIndex)) then
          STDRAD = matdata(rowIndex,2) * 100000.d0
          call obs_headSet_r(obsdat,columnIndex,headerIndex,  STDRAD)
          rowIndex = rowIndex+1
        end if
      end do

      deallocate(matdata)
      call fSQL_free_mem    ( stmt )

    end do

    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine sqlr_readSqlite_avhrr


  !--------------------------------------------------------------------------
  ! sqlr_readSqlite
  !--------------------------------------------------------------------------
  subroutine sqlr_readSqlite(obsdat, familyType, fileName )
    !
    ! :Purpose: To read SQLite namelist and files.
    ! 
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsdat     ! ObsSpaceData Structure
    character(len=*) , intent(in)    :: familyType ! Family Type 
    character(len=*) , intent(in)    :: fileName   ! SQLite filename

    ! Locals:
    character(len=9)         :: rdbSchema
    character(len=12)        :: idStation
    integer                  :: codeType, obsDate, obsTime, obsStatus, obsFlag, obsVarno
    integer(8)               :: headPrimaryKey,  bodyPrimaryKey, nextHeadPrimaryKey
    integer(8), allocatable  :: bodyHeadKeys(:), bodyPrimaryKeys(:)
    integer(8), allocatable  :: tempBodyKeys(:,:)
    integer                  :: numBodyKeys, firstBodyIndexOfThisBatch
    logical                  :: createHeaderEntry
    real(pre_obsReal)        :: elevReal, xlat, xlon, vertCoord
    real                     :: elev    , obsLat, obsLon, elevFact
    real                     :: obsrzam , obsrans, obsrane, obsrele          
    integer                  :: vertCoordType, vertCoordFact, fnom, fclos, nulnam, ierr, idProf
    real                     :: zenithReal, solarZenithReal, CloudCoverReal, solarAzimuthReal
    integer                  :: roQcFlag
    real(pre_obsReal)        :: geoidUndulation, earthLocRadCurv, obsValue, surfEmiss, biasCorrection
    real(8)                  :: geoidUndulation_R8, earthLocRadCurv_R8, azimuthReal_R8
    integer                  :: trackCellNum
    real(pre_obsReal)        :: modelWindSpeed, range_m
    real(8)                  :: modelWindSpeed_R8
    integer                  :: iasiImagerCollocationFlag, iasiGeneralQualityFlag
    integer                  :: obsSat, landSea, terrainType, instrument, sensor, numberElem
    integer                  :: bitIndex, rowIndex, obsNlv, headerIndex, bodyIndex
    integer                  :: bitsFlagOn, bitsFlagOff, numBody, numHeader
    real(pre_obsReal), parameter :: zemFact = 0.01
    character(len=512)       :: query, queryData, queryHeader, queryIDs
    character(len=256)       :: selectIDs, selectData
    character(len=256)       :: cfgSqlite, csqlcrit, columnsHeader, columnsData
    logical                  :: finished
    real(8), allocatable     :: matdata(:,:)
    type(fSQL_DATABASE)      :: db         ! type for SQLIte  file handle
    type(fSQL_STATEMENT)     :: stmt,stmt2,stmt3 ! type for precompiled SQLite statements
    type(fSQL_STATUS)        :: stat,stat2 !type for error status
    character(len=256)       :: listElem, sqlExtraDat, sqlExtraHeader, sqlNull, sqlLimit 
    character(len=256)       :: sqlDataOrder
    character(len=256),allocatable :: listElemArray(:)
    integer,allocatable            :: listElemArrayInteger(:)
    integer                  :: numberBitsOff, numberBitsOn, bitsOff(15), bitsOn(15), numberRows, numberColumns
    character(len=*), parameter :: myName = 'sqlr_readSqlite'
    namelist /NAMSQLamsua/numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLamsub/numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLairs/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLiasi/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLcris/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLcrisfsr/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
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
    namelist /NAMSQLgl/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLradar/numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLradvel/numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn

    write(*,*) myName//': fileName   : ', trim(fileName)
    write(*,*) myName//': familyType : ', trim(familyType)
    call fSQL_open( db, trim(fileName) ,stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myName//': fSQL_open '//fSQL_errmsg(stat) )
    end if

    query = "select schema from rdb4_schema ;"
    rdbSchema = sqlr_query(db,trim(query))
    write(*,'(4a)') myName//' rdbSchema is ---> ', trim(rdbSchema)

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
      numberElem = 1
    else
      write(listElem,'(i5.5,",",i5.5)') bufr_nedd, bufr_neff
    end if

    vertCoordfact  = 1
    columnsData = " vcoord, varno, obsvalue, flag "
    if ( trim(familyType) == 'TO' ) then      
      columnsHeader = " id_obs, lat, lon, codtyp, date, time, id_stn, status, id_sat, land_sea, instrument, zenith, solar_zenith "
      vertCoordType = 3
    else if ( trim(familyType) == 'GL' ) then
      columnsHeader = " id_obs, lat, lon, codtyp, date, time, id_stn"
    else if ( trim(rdbSchema) == 'ra' ) then
      columnsHeader = " id_obs, lat, lon, codtyp, date, time, id_stn"
      vertCoordType = 1
    else if ( trim(rdbSchema) == 'radvel') then
      columnsHeader = " id_obs, lat, lon, codtyp, date, time, id_stn, antenna_altitude, center_azimuth, center_elevation, range_start, range_end "
      vertCoordType = 1
    else
      columnsHeader = " id_obs, lat, lon, codtyp, date, time, id_stn, status, elev"  
    end if

    nulnam = 0
    ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0 )
    select case(trim(rdbSchema))
      case ( 'ua' )
        vertCoordFact = 1
        vertCoordType = 2
        read(nulnam, nml = NAMSQLua, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLua )
      case ( 'ai' )
        vertCoordFact = 1
        vertCoordType = 2
        read(nulnam, nml = NAMSQLai, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLai )
      case ( 'sw' )
        vertCoordFact = 1
        vertCoordType = 2
        read(nulnam, nml = NAMSQLsw, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLsw )
      case ( 'pr' )
        vertCoordFact = 1
        vertCoordType = 1
        read(nulnam, nml = NAMSQLpr, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLpr )
      case ( 'al' )  
        columnsHeader = trim(columnsHeader)//", id_prof"
        vertCoordFact = 1
        vertCoordType = 1
        read(nulnam, nml = NAMSQLal, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLal )
      case ( 'ro' )     
        columnsHeader = trim(columnsHeader)//",ro_qc_flag, geoid_undulation, earth_local_rad_curv, id_sat, azimuth"
        vertCoordFact = 1
        vertCoordType = 1
        read(nulnam, nml = NAMSQLro, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLro )
      case ( 'sf' )
        vertCoordFact = 0
        vertCoordType = 1
        read(nulnam, nml = NAMSQLsfc, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLsfc )
      case ( 'sst' )
        columnsHeader = trim(columnsHeader)//", solar_zenith "
        vertCoordFact = 0
        vertCoordType = 1
        read(nulnam, nml = NAMSQLsfc, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLsfc )
      case ( 'scat' )
        vertCoordFact = 0
        vertCoordType = 1
        read(nulnam, nml = NAMSQLsc, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLsc )
      case( 'airs' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, cloud_cover, solar_azimuth "
        write(columnsData,'(a,f3.2,a)') trim(columnsData)//", ifnull(surf_emiss,", tvs_defaultEmissivity, "), bias_corr "
        read(nulnam, nml = NAMSQLairs, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLairs )
      case( 'iasi' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, cloud_cover, solar_azimuth, FANION_QUAL_IASI_SYS_IND, INDIC_NDX_QUAL_GEOM "
        columnsData = trim(columnsData)//", surf_emiss, bias_corr "
        read(nulnam, nml = NAMSQLiasi, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLiasi )
      case( 'cris' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, cloud_cover, solar_azimuth "
        columnsData = trim(columnsData)//", surf_emiss, bias_corr "
        read(nulnam, nml = NAMSQLcris, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLcris )
      case( 'crisfsr' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, cloud_cover, solar_azimuth "
        columnsData = trim(columnsData)//", surf_emiss "
        read(nulnam, nml = NAMSQLcrisfsr, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLcrisfsr )
      case( 'amsua' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, sensor, solar_azimuth "
        columnsData = trim(columnsData)//", bias_corr "
        read(nulnam, nml = NAMSQLamsua, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLamsua )
      case( 'amsub' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, sensor, solar_azimuth "
        columnsData = trim(columnsData)//", bias_corr "
        read(nulnam, nml = NAMSQLamsub, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLamsub )
      case( 'atms')
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, sensor, solar_azimuth "
        columnsData = trim(columnsData)//", bias_corr "
        read(nulnam, nml = NAMSQLatms, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLatms )
      case( 'ssmi' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type "
        columnsData = trim(columnsData)//", bias_corr "
        read(nulnam, nml = NAMSQLssmi, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLssmi )
      case( 'csr' )
        columnsData = trim(columnsData)//", bias_corr "
        read(nulnam, nml = NAMSQLcsr, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLcsr )
      case( 'gl' )
        read(nulnam, nml = NAMSQLgl, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLgl )
      case( 'gl_ascat' )
        columnsHeader = trim(columnsHeader)//", track_cell_no, mod_wind_spd "
        read(nulnam, nml = NAMSQLgl, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLgl )
      case( 'ra' )
        read(nulnam, nml = NAMSQLradar, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLradar ) 
      case( 'radvel' )
        columnsHeader = trim(columnsHeader) 
        ! add  range to data columns to read
        columnsData = trim(columnsData)//", range "
        read(nulnam, nml = NAMSQLradvel, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLradvel )
      case DEFAULT
        call utl_abort( myName//': Unsupported  SCHEMA in SQLITE file!'//trim(rdbSchema) )
    end select
    ierr=fclos( nulnam )

    ! Set the observation variable transforms
    if (numberElem > 0) then
      call utl_splitString(listElem,',', & ! IN
                           listElemArray)  ! OUT
      call utl_stringArrayToIntegerArray(listElemArray,      & ! IN
                                         listElemArrayInteger) ! OUT
      call ovt_setup(listElemArrayInteger) ! IN
      deallocate(listElemArrayInteger)
      deallocate(listElemArray)
    end if

    ! Compose SQL queries
    bitsFlagOff=0
    do bitIndex = 1, numberBitsOff
      bitsFlagOff = ibset ( bitsFlagOff, 13 - bitsOff(bitIndex) )
    end do
    bitsFlagOn=0
    do bitIndex = 1, numberBitsOn
      bitsFlagOn = ibset ( bitsFlagOn, 13 - bitsOn(bitIndex) )
    end do

    write(cfgSqlite, '(i6)' ) bitsFlagOff
    csqlcrit=trim(" (flag & "//cfgSqlite)//" = 0 "

    if ( numberBitsOn > 0 ) then
      write(cfgSqlite, '(i6)' ) bitsFlagOn
      csqlcrit = trim(csqlcrit)//trim(" and flag & "//cfgSqlite)
      csqlcrit = trim(csqlcrit)//" = "//cfgSqlite
    end if

    !ordering of data that will get read in matdata, bodyPrimaryKeys and bodyHeadKeys
    sqlDataOrder = ' order by id_obs, varno, id_data' 

    selectIDs  = "select id_data, id_obs"
    selectData = "select "//trim(columnsData)
    csqlcrit = trim(csqlcrit)//" or flag is null) and varno in ( "//trim(listElem)//" )"//trim(sqlExtraDat)//trim(SQLNull)
    !it is very important that queryIDs and queryData be identical except for the column names being selected
    queryIDs  =  trim(selectIDs)//trim(" from data where ")//trim(csqlcrit)//trim(sqlDataOrder)//trim(sqlLimit)//";"
    queryData = trim(selectData)//trim(" from data where ")//trim(csqlcrit)//trim(sqlDataOrder)//trim(sqlLimit)//";"
    write(*,'(4a)') myName//': ',trim(rdbSchema),' queryIDs     --> ', trim(queryIDs)
    write(*,'(4a)') myName//': ',trim(rdbSchema),' queryData    --> ', trim(queryData)
    !the first queryHeader is printed below in the BODY loop

    if ( trim(rdbSchema)=='pr' .or. trim(rdbSchema)=='sf' .or. trim(rdbSchema)=='sst' ) then
      elevFact=1.
    else
      elevFact=0.
    end if

    call fSQL_prepare( db, trim(queryData), stmt2, stat2 )
    if ( fSQL_error(stat2) /= FSQL_OK ) then
      call sqlr_handleError(stat2,'fSQL_prepare has failed for queryData, &
                                   this could be caused by an "order by"  &
                                   statement in sqlExtraDat ')
    end if

    headerIndex  = obs_numHeader(obsdat)
    bodyIndex = obs_numBody(obsdat)
    numHeader = obs_numHeader(obsdat)
    numBody   = obs_numBody(obsdat)
    write(*,*) myName//': DEBUT numheader  =', numHeader
    write(*,*) myName//': DEBUT numbody    =', numBody
    call fSQL_get_many( stmt2, nrows=numberRows, ncols=numberColumns, mode=FSQL_REAL8 )
    write(*,*) myName//':  numberRows numberColumns =', numberRows, numberColumns
    write(*,*) myName//':  rdbSchema = ', rdbSchema
    write(*,*) myName//': =========================================='
    allocate( matdata(numberRows, numberColumns) )
    matdata = 0.0d0
    call fSQL_fill_matrix( stmt2, matdata )
    call fSQL_free_mem( stmt2 )
    call fSQL_finalize( stmt2 )

    ! read id_data and id_obs columns in the body table
    call fSQL_prepare( db, trim(queryIDs) , stmt3, status=stat )
    ! note: "status" not set when getting integers
    allocate( bodyPrimaryKeys(numberRows) )
    allocate( bodyHeadKeys(numberRows) )
    allocate( tempBodyKeys(numberRows,2) )
    call fSQL_get_many( stmt3, nrows=numberRows, ncols=numberColumns, &
                        mode=FSQL_INT8 )
    call fSQL_fill_matrix( stmt3, tempBodyKeys )
    bodyPrimaryKeys(:) = tempBodyKeys(:,1)
    bodyHeadKeys(:)    = tempBodyKeys(:,2)
    deallocate(tempBodyKeys)
    call fSQL_free_mem( stmt3 )
    call fSQL_finalize( stmt3 )

    ! Here is a summary of what is going on in the rest of this routine:
    ! 
    ! for each row of the body table
    !
    !   if this if the first id_data associated with a given id_obs
    !     -> read associated header
    !
    !   -> write body entry in obsspacedata
    !
    !   if this is the last id_data associated with a given id_obs
    !     -> write header entry in obsspacedata
    !
    obsNlv = 0
    numBodyKeys = size(bodyPrimaryKeys)
    if ( numBodyKeys /= numberRows ) then
      call utl_abort( 'sqlr_readSqlite: number of body keys not equal to number of rows in matdata' )
    end if 
    BODY: do rowIndex = 1, numBodyKeys

      !"id_data" and "id_obs" values for this entry
      bodyPrimaryKey = bodyPrimaryKeys(rowIndex)
      headPrimaryKey = bodyHeadKeys(rowIndex)
      bodyIndex = bodyIndex + 1
      obsNlv    = obsNlv + 1
      

      READHEADER: if ( obsNlv == 1 ) then

        !this is the first body entry associated with a new headPrimaryKey
        firstBodyIndexOfThisBatch = bodyIndex

        !we read the associated header
        queryHeader = "select "//trim(columnsHeader)//" from header "//trim(sqlExtraHeader)//" where id_obs = ? "
        if ( rowIndex == 1 ) then
          write(*,'(4a)') myName//': ',trim(rdbSchema),' first queryHeader    --> ', trim(queryHeader)
        end if
        call fSQL_prepare( db, queryHeader, stmt, stat )
        if ( fSQL_error(stat)  /= FSQL_OK ) then
          write(*,*) 'sqlr_readSqlite: problem reading header entry. Query:', trim(queryHeader)
          call sqlr_handleError( stat, 'fSQL_prepare')
        end if
         
        call fSQL_bind_param(stmt, PARAM_INDEX = 1, INT8_VAR = headPrimaryKey )

        call fSQL_get_row( stmt, finished )
        if ( finished ) then
          write(*,*) 'sqlr_readSqlite: problem reading header entry. Query:', trim(queryHeader)
          call utl_abort( 'Problem with fSQL_get_row()' )
        end if

        ! The query result is inserted into variables
        terrainType = MPC_missingValue_INT
        landSea     = MPC_missingValue_INT
        sensor      = MPC_missingValue_INT
        cloudCoverReal = MPC_missingValue_R4
        elev = 0.; elevReal = 0.
        solarAzimuthReal = MPC_missingValue_R4
        solarZenithReal  = MPC_missingValue_R4
        zenithReal = MPC_missingValue_R4
        instrument = MPC_missingValue_INT
        azimuthReal_R8 = MPC_missingValue_R8
        obsrzam = MPC_missingValue_R4
        obsrele = MPC_missingValue_R4
        obsrans = MPC_missingValue_R4
        obsrane = MPC_missingValue_R4

        call fSQL_get_column( stmt, COL_INDEX = 2, REAL_VAR  = obsLat    )
        call fSQL_get_column( stmt, COL_INDEX = 3, REAL_VAR  = obsLon    )
        call fSQL_get_column( stmt, COL_INDEX = 4, INT_VAR   = codeType  )
        call fSQL_get_column( stmt, COL_INDEX = 5, INT_VAR   = obsDate   )
        call fSQL_get_column( stmt, COL_INDEX = 6, INT_VAR   = obsTime   )
        call fSQL_get_column( stmt, COL_INDEX = 7, CHAR_VAR  = idStation )


        FAMILYEXCEPTIONS: if ( trim(familyType) == 'TO' ) then

          call fSQL_get_column( stmt, COL_INDEX =  8,  INT_VAR = obsStatus )
          call fSQL_get_column( stmt, COL_INDEX =  9,  INT_VAR = obsSat    )
          call fSQL_get_column( stmt, COL_INDEX = 10,  INT_VAR = landSea        , INT_MISSING=MPC_missingValue_INT )
          call fSQL_get_column( stmt, COL_INDEX = 11,  INT_VAR = instrument     , INT_MISSING=MPC_missingValue_INT )
          call fSQL_get_column( stmt, COL_INDEX = 12, REAL_VAR = zenithReal     , REAL_MISSING=MPC_missingValue_R4 )
          call fSQL_get_column( stmt, COL_INDEX = 13, REAL_VAR = solarZenithReal, REAL_MISSING=MPC_missingValue_R4 )

          if ( trim(rdbSchema) /= 'csr' ) then

            call fSQL_get_column( stmt, COL_INDEX = 14, REAL8_VAR  = azimuthReal_R8 )
           
            call fSQL_get_column( stmt, COL_INDEX = 15, INT_VAR   = terrainType, INT_MISSING=MPC_missingValue_INT )

          end if

          if ( trim(rdbSchema) == 'airs' .or. trim(rdbSchema) == 'iasi' .or. trim(rdbSchema) == 'cris' ) then

            call fSQL_get_column( stmt, COL_INDEX = 16, REAL_VAR = cloudCoverReal, REAL_MISSING=MPC_missingValue_R4   )
            call fSQL_get_column( stmt, COL_INDEX = 17, REAL_VAR = solarAzimuthReal, REAL_MISSING=MPC_missingValue_R4 )

          else if ( trim(rdbSchema) == 'amsua' .or. trim(rdbSchema) == 'amsub' .or. trim(rdbSchema) == 'atms' ) then

            call fSQL_get_column( stmt, COL_INDEX = 16, INT_VAR  = sensor, INT_MISSING=MPC_missingValue_INT )
            call fSQL_get_column( stmt, COL_INDEX = 17, REAL_VAR = solarAzimuthReal, REAL_MISSING=MPC_missingValue_R4 )

          end if
          if ( trim(rdbSchema) == 'iasi'  ) then
            call fSQL_get_column( stmt, COL_INDEX = 18, INT_VAR  = iasiGeneralQualityFlag,    INT_MISSING=MPC_missingValue_INT )
            call fSQL_get_column( stmt, COL_INDEX = 19, INT_VAR  = iasiImagerCollocationFlag, INT_MISSING=MPC_missingValue_INT )
          end if

          if ( instrument == 420 ) obsSat  = 784
          if ( codeType == 202 .and. instrument == 620 ) instrument = 2046
          if ( sensor == MPC_missingValue_INT ) then
            sensor = 0
            if (instrument == MPC_missingValue_INT ) instrument = 0
          else
            instrument = obsu_cvt_obs_instrum(sensor)
          end if

        else if ( trim(familyType) == 'GL' ) then

          ! It does not have the obsStatus column.

          if ( idStation(1:6) == 'METOP-') then
            call fSQL_get_column( stmt, COL_INDEX = 8, INT_VAR  = trackCellNum )
            if (trackCellNum > 21) trackCellNum = 43 - trackCellNum
            call fSQL_get_column( stmt, COL_INDEX = 9, REAL8_VAR  = modelWindSpeed_R8 )
            modelWindSpeed = modelWindSpeed_R8
          end if

      	else if ( trim(rdbSchema) == 'sst' ) then
        	! satellite SST observations

        	call fSQL_get_column( stmt, COL_INDEX =  8,  INT_VAR = obsStatus                                         )
        	call fSQL_get_column( stmt, COL_INDEX =  9, REAL_VAR = elev           , REAL_MISSING=MPC_missingValue_R4 )
        	call fSQL_get_column( stmt, COL_INDEX = 10, REAL_VAR = solarZenithReal, REAL_MISSING=MPC_missingValue_R4 )

        else if ( trim(familyType) == 'RA' ) then
          if ( trim(rdbSchema) == 'radvel') then
            call fSQL_get_column( stmt, COL_INDEX = 8, REAL_VAR  = elev, REAL_MISSING=MPC_missingValue_R4 )
            elevReal=elev
            call fSQL_get_column( stmt, COL_INDEX = 9,  REAL_VAR  = obsrzam)
            call fSQL_get_column( stmt, COL_INDEX = 10, REAL_VAR  = obsrele)
            call fSQL_get_column( stmt, COL_INDEX = 11, REAL_VAR  = obsrans)
            call fSQL_get_column( stmt, COL_INDEX = 12, REAL_VAR  = obsrane)
          end if

        else  ! remaining families

          call fSQL_get_column( stmt, COL_INDEX = 8,  INT_VAR  = obsStatus )
          call fSQL_get_column( stmt, COL_INDEX = 9, REAL_VAR  = elev, REAL_MISSING=MPC_missingValue_R4 )
          elevReal=elev

          if ( trim(rdbSchema)=='ro' ) then

            call fSQL_get_column( stmt, COL_INDEX = 10, INT_VAR   = roQcFlag, INT_MISSING=MPC_missingValue_INT )
            call fSQL_get_column( stmt, COL_INDEX = 11, REAL8_VAR = geoidUndulation_R8 )
            geoidUndulation = geoidUndulation_R8
            call fSQL_get_column( stmt, COL_INDEX = 12, REAL8_VAR = earthLocRadCurv_R8 )
            earthLocRadCurv = earthLocRadCurv_R8
            call fSQL_get_column( stmt, COL_INDEX = 13, INT_VAR   = obsSat, INT_MISSING=MPC_missingValue_INT )
            call fSQL_get_column( stmt, COL_INDEX = 14, REAL8_VAR = azimuthReal_R8 )

          else if ( trim(rdbSchema)=='al' ) then
            call fSQL_get_column( stmt, COL_INDEX = 10, INT_VAR   = idProf )
          end if

        end if FAMILYEXCEPTIONS

        !we are done reading this header entry
        call fSQL_finalize( stmt )

        if ( obsLon < 0. ) obsLon = obsLon + 360.
        xlat = obsLat * MPC_RADIANS_PER_DEGREE_R8
        xlon = obsLon * MPC_RADIANS_PER_DEGREE_R8

      end if READHEADER

      ! From this point on, we read this body entry and add it to obsspacedata
              
      vertCoord = matdata(rowIndex,1)
      obsVarno = int(matdata(rowIndex,2))
      obsValue = matdata(rowIndex,3)
      obsFlag = int(matdata(rowIndex,4))

      call obs_setBodyPrimaryKey(obsdat, bodyIndex, bodyPrimaryKey)

      if (trim(rdbSchema) == 'airs' .or. trim(rdbSchema) == 'iasi' .or. trim(rdbSchema) == 'cris' ) then

        surfEmiss = matdata(rowIndex,5)
        call obs_bodySet_r(obsdat, OBS_SEM, bodyIndex, surfEmiss * zemFact)

        biasCorrection = matdata(rowIndex,6)

        if ( obs_columnActive_RB(obsdat,OBS_BCOR) ) then
          call obs_bodySet_r(obsdat, OBS_BCOR, bodyIndex, biasCorrection )
        end if

      end if

      if ( trim(rdbSchema) == 'amsua' .or. trim(rdbSchema) == 'amsub' .or. &
           trim(rdbSchema) == 'atms'  .or.   trim(rdbSchema) =="ssmi" .or. &
           trim(rdbSchema) =="csr" ) then

        biasCorrection = matdata(rowIndex,5)

        if ( obs_columnActive_RB(obsdat,OBS_BCOR) ) &
             call obs_bodySet_r(obsdat, OBS_BCOR, bodyIndex, biasCorrection )

      end if

      if (trim(rdbSchema) == 'radvel') then

        !matdata is initialized with 0.0d0 
        !if vcoord is missing in observation file its value will be set to 0.0d0
        !we change it to missing to make it more obvious that this value has to be calculated later
        if (vertCoord == real(0.0d0, pre_obsReal)) then
          vertCoord = real(MPC_missingValue_R8, pre_obsReal)
        end if

        !write standard body values to obsSpaceData
        call sqlr_initData(obsdat, vertCoord, obsValue, obsVarno, obsFlag, vertCoordType, bodyIndex)

        !add range along radar beam
        range_m = matdata(rowIndex,5)
        call obs_bodySet_r(obsdat, OBS_LOCI, bodyIndex, range_m )

      end if
        
      if ( trim(familyType) == 'TO' ) then

        if ( obsValue /= MPC_missingValue_R8 ) then
          call sqlr_initData(obsdat, vertCoord, obsValue, obsVarno, obsFlag, vertCoordType, bodyIndex)
        else
          call sqlr_initData(obsdat, vertCoord, real(MPC_missingValue_R8,pre_obsReal), obsVarno, obsFlag, vertCoordType, bodyIndex)
        endif

      else 

        call sqlr_initData( obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                            obsValue, obsVarno, obsFlag, vertCoordType, bodyIndex)


        if (.not. filt_bufrCodeAssimilated(obsVarno) .and. &
            .not. ovt_bufrCodeSkipped(obsVarno)) then
          
          ! Add an extra row to the obsSpaceData body table
          ! to contain quantity later calculated by ovt_transformObsValue
          call obs_setBodyPrimaryKey( obsdat, bodyIndex+1, -1)
          call sqlr_initData( obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                              obs_missingValue_R, ovt_getDestinationBufrCode(obsVarno), &
                              0, vertCoordType, bodyIndex + 1 )
          bodyIndex = bodyIndex + 1
          obsNlv = obsNlv + 1
          if (ovt_isWindObs(obsVarno)) then
            ! Add an extra row for the other wind component
            call obs_setBodyPrimaryKey( obsdat, bodyIndex+1, -1)
            call sqlr_initData( obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                                obs_missingValue_R, ovt_getDestinationBufrCode(obsVarno,extra_opt=.true.), &
                                0, vertCoordType, bodyIndex + 1 )
            bodyIndex = bodyIndex + 1
            obsNlv = obsNlv + 1
          end if
        end if
      end if     

      !Logic to determine whether we need to create a header entry
      createHeaderEntry = .false.
      if (rowIndex == numBodyKeys) then
        !This is the last body entry, nothing comes after.
        createHeaderEntry = .true.
      else 
        nextHeadPrimaryKey = bodyHeadKeys(rowIndex+1)
        if ( headPrimaryKey /= nextHeadPrimaryKey ) then
          !we just treated the last body entry associated with this id_obs (headPrimaryKey)
          !make header for this id_obs
          createHeaderEntry = .true.
        end if
      end if
      WRITEHEADER: if (createHeaderEntry) then

        headerIndex = headerIndex + 1 

        call sqlr_initHeader( obsdat, rdbSchema, familyType, headerIndex, elevReal, obsSat, &
             real(azimuthReal_R8,kind=pre_obsReal), geoidUndulation, earthLocRadCurv,       &
             roQcFlag, instrument, real(zenithReal,kind=pre_obsReal),                       &
             real(cloudCoverReal,kind=pre_obsReal), real(solarZenithReal,kind=pre_obsReal), &
             real(solarAzimuthReal,kind=pre_obsReal), terrainType, landSea,                 &
             iasiImagerCollocationFlag, iasiGeneralQualityFlag, headPrimaryKey,             &
             xlat, xlon, codeType, obsDate, obsTime/100, obsStatus, idStation, idProf,      &
             trackCellNum, modelWindSpeed,                                                  &
             real(obsrzam,kind=pre_obsReal), real(obsrele,kind=pre_obsReal),                &
             real(obsrans,kind=pre_obsReal), real(obsrane,kind=pre_obsReal),                &
             firstBodyIndexOfThisBatch, obsNlv)

        !reset level counter for next batch of data entries
        obsNlv = 0

      end if WRITEHEADER

    end do BODY

    deallocate(matdata)
    numHeader = obs_numHeader(obsdat)
    numBody   = obs_numBody(obsdat)
    write(*,*) myName//': FIN numheader  =', numHeader
    write(*,*) myName//': FIN numbody    =', numBody
    write(*,*) myName//': fin header '

    call fSQL_close( db, stat ) 
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'sqlr_readSqlite: problem closing sqlite db', trim(fileName)
      call sqlr_handleError( stat, 'fSQL_close')
    end if

  end subroutine sqlr_readSqlite

  !--------------------------------------------------------------------------
  ! sqlr_query
  !--------------------------------------------------------------------------
  function sqlr_query(db,query)
    !
    ! :Purpose: To create a query to read an SQLite file
    !
    implicit none

    ! arguments
    type(fSQL_DATABASE)  :: db    ! type handle for  SQLIte file
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
    if ( .not. finished ) write(*,*) ' SQL QUERY ---> QUERY RETURNS MORE THAN ONE ROW...  '
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

  !--------------------------------------------------------------------------
  ! sqlr_updateSqlite
  !--------------------------------------------------------------------------
  subroutine sqlr_updateSqlite( db, obsdat, familyType, fileName, fileNumber )
    !
    ! :Purpose: update SQLite files. List of items to update is in the 
    !           namSQLUpdate namelist
      
    implicit none
    
    ! arguments
    type(fSQL_database), intent(inout) :: db         ! SQL database
    type (struct_obs)  , intent(inout) :: obsdat     ! obsSpaceData
    character(len=*)   , intent(in)    :: fileName   ! file name  
    character(len=*)   , intent(in)    :: familyType ! Observation Family Type 
    integer            , intent(in)    :: fileNumber ! FILE NUMBER ASSOCIATED WITH db
    
    ! locals
    type(fSQL_statement)        :: stmt ! prepared statement for  SQLite
    type(fSQL_status)           :: stat ! type error status
    integer                     :: obsRln, obsNlv, obsIdf, obsFlag
    integer                     :: obsStatus, last_question, landSea, terrainType
    integer(8)                  :: headPrimaryKey, bodyPrimaryKey
    integer                     :: itemIndex, headPrimaryKeyIndex, landSeaIndex
    integer                     :: headerIndex, bodyIndex, numberUpdateItems, numberUpdateItemsRadar
    character(len =   3)        :: item, itemUpdateList(15), itemUpdateListRadar(15)
    integer                     :: updateList(20), fnom, fclos, nulnam, ierr
    character(len =  10)        :: columnName
    character(len = 128)        :: query
    character(len = 356)        :: itemChar, columnNameChar
    logical                     :: back
    real                        :: romp, obsValue, scaleFactor
    character(len=*), parameter :: myName = 'sqlr_updateSqlite'
    character(len=*), parameter :: myWarning = '****** '// myName //': WARNING: '
    namelist/namSQLUpdate/ numberUpdateItems,      itemUpdateList,     &
                           numberUpdateItemsRadar, itemUpdateListRadar

    write(*,*) myName//': Starting ===================  '

    ! set default values of namelist variables
    itemUpdateList(:) = ''
    itemUpdateListRadar(:) = ''
    numberUpdateItems = 0
    numberUpdateItemsRadar = 0

    ! Read the namelist for directives
    nulnam = 0
    ierr   = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0 )
    read(nulnam,nml = namSQLUpdate, iostat = ierr )
    if ( ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
    if ( mpi_myid == 0 ) write(*, nml = namSQLUpdate )
    ierr = fclos( nulnam )

    ! Append extra sqlite columns to update to itemUpdateList
    if (trim( familyType ) == 'RA') then
      do itemIndex = 1, numberUpdateItemsRadar
        numberUpdateItems = numberUpdateItems + 1
        itemUpdateList(numberUpdateItems) = itemUpdateListRadar( itemIndex )
      end do
    end if

    write(*,*) myName//': family type   = ', trim(familyType)
    write(*,*) myName//': number of items to update: ', numberUpdateItems
    write(*,*) myName//': file name     = ', trim(fileName)
    write(*,*) myName//': missing value = ', MPC_missingValue_R8    

    ! create query
    itemChar='  '

    do itemIndex = 1, numberUpdateItems
    
      item = itemUpdateList( itemIndex )
      write(*,*) myName//': updating ', itemIndex, trim(item)
      select case(item)
        case('OMA')
          updateList( itemIndex ) = OBS_OMA
          columnName = 'oma'
        case('OMP')
          updateList( itemIndex ) = OBS_OMP
          columnName = 'omp'
        case('VAR')
          updateList( itemIndex ) = OBS_VAR
          columnName = 'obsvalue'
        case('OER')
          updateList( itemIndex ) = OBS_OER
          columnName = 'obs_error'
        case('FGE')
          updateList( itemIndex ) = OBS_HPHT
          columnName = 'fg_error'
        case('EMI')
          updateList( itemIndex ) = OBS_SEM
          columnName = 'surf_emiss'
        case('COR')
          updateList( itemIndex ) = OBS_BCOR
          columnName = 'bias_corr'
        case('ALT')
          updateList( itemIndex ) = OBS_PPP
          columnName = 'vcoord'
        case DEFAULT
          call utl_abort( myName//': invalid item '// columnName //' EXIT sqlr_updateSQL!!!' )
      end select
      
      if ( sqlu_sqlColumnExists( fileName, 'data', columnName ) == .true. ) then
        itemChar = trim(itemChar)//','// trim(columnName) // trim(' = ? ')
      else
        write(*,*) myWarning//': column '// columnName// &
                   ' does not exist in the file '//trim(fileName)
      end if	
    end do

    back=.true.
    last_question  = scan(itemChar, '?', back)
    columnNameChar = itemChar(1:last_question)
    itemChar    = columnNameChar
    
    query = ' update data set flag = ? '//trim( itemChar )
    query = trim(query)//' where id_data = ?  ;'
    write(*,*) myName//': update query --->  ', query
    call fSQL_do_many( db, 'PRAGMA  synchronous = OFF; PRAGMA journal_mode = OFF;' )
    call fSQL_prepare( db, query , stmt, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat, 'fSQL_prepare : ')
    call fSQL_begin(db)

    HEADER: do headerIndex = 1, obs_numHeader(obsdat)
 
      obsIdf = obs_headElem_i( obsdat,OBS_IDF, headerIndex )
 
      if ( obsIdf /= fileNumber ) cycle HEADER
      headPrimaryKey = obs_headPrimaryKey( obsdat, headerIndex )
      obsRln = obs_headElem_i( obsdat, OBS_RLN, headerIndex )
      obsNlv = obs_headElem_i( obsdat, OBS_NLV, headerIndex )
	
      BODY: do bodyIndex = obsRln, obsNlv + obsRln - 1

        obsFlag = obs_bodyElem_i( obsdat, OBS_FLG, bodyIndex )
        bodyPrimaryKey  = obs_bodyPrimaryKey(obsdat, bodyIndex)

        call fSQL_bind_param(stmt, PARAM_INDEX = 1, INT_VAR = obsFlag )
        
				ITEMS: do itemIndex = 1, numberUpdateItems
          
	  		  obsValue = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          if ( obsValue /= obs_missingValue_R ) then  
            romp = obs_bodyElem_r(obsdat, updateList( itemIndex ), bodyIndex )
            if ( romp == obs_missingValue_R ) then
              call fSQL_bind_param(stmt, PARAM_INDEX = itemIndex + 1 ) ! sql null values
            else
              scaleFactor=1.0
              if ( updateList( itemIndex ) == OBS_SEM ) scaleFactor=100.0
              call fSQL_bind_param(stmt, PARAM_INDEX = itemIndex + 1, REAL_VAR = romp*scaleFactor )
            end if
          else
            call fSQL_bind_param(stmt, PARAM_INDEX = itemIndex + 1)  ! sql null values
          end if

        end do ITEMS

        call fSQL_bind_param(stmt, PARAM_INDEX = numberUpdateItems + 2, INT8_VAR  = bodyPrimaryKey )
        call fSQL_exec_stmt (stmt)

      end do BODY

    end do HEADER

    call fSQL_finalize( stmt )

    if ( trim( familyType ) /= 'GL'.and. trim( familyType ) /= 'RA' ) then

      ! UPDATES FOR THE STATUS FLAGS and land_sea (for satellites) IN THE HEADER TABLE
      ! The variables 'headPrimaryKeyIndex' and 'landSeaIndex' are defined here and used below.
      ! Any change in this logic must be coherent with the code below!
      if ( trim( familyType ) == 'TO' ) then
        query = ' update header set status  = ?, land_sea= ? where id_obs = ? '
        landSeaIndex = 2
        headPrimaryKeyIndex = 3
      else
        query = ' update header set status  = ?  where id_obs = ? '
        headPrimaryKeyIndex = 2
      end if
      call fSQL_prepare( db, query , stmt, stat)
      if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat,'fSQL_prepare : ')

      HEADER2: do headerIndex = 1,obs_numHeader(obsdat)
         
	terrainType=MPC_missingValue_INT
        obsIdf = obs_headElem_i(obsdat, OBS_IDF, headerIndex )
        if ( obsIdf /= fileNumber ) cycle HEADER2
        headPrimaryKey = obs_headPrimaryKey(obsdat, headerIndex)
        obsStatus = obs_headElem_i(obsdat, OBS_ST1, headerIndex )
        landsea   = obs_headElem_i(obsdat, OBS_STYP,headerIndex )
        call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = obsStatus )
        ! The variables 'headPrimaryKeyIndex' and 'landSeaIndex' are defined above and
        ! they must be coherent with the query designed above
        call fSQL_bind_param( stmt, PARAM_INDEX = headPrimaryKeyIndex, INT8_VAR = headPrimaryKey )
        if ( trim( familyType ) == 'TO' ) then
          call fSQL_bind_param( stmt, PARAM_INDEX = landSeaIndex, INT_VAR  = landSea )
        end if

        call fSQL_exec_stmt ( stmt)

      end do HEADER2
    
      call fSQL_finalize( stmt )

    end if

    call fSQL_commit(db, stat)
    if ( fSQL_error(stat)  /= FSQL_OK ) then
      call sqlr_handleError( stat, myName//'fSQL_commit')
    end if
    write(*,*) myName//': End ===================  ', trim(familyType)

  end subroutine sqlr_updateSqlite

  !--------------------------------------------------------------------------
  ! sqlr_addCloudParametersandEmissivity
  !--------------------------------------------------------------------------
  subroutine sqlr_addCloudParametersandEmissivity( db, obsdat,fileNumber )
    implicit none
    ! arguments
    type(fSQL_DATABASE), intent(inout) :: db   ! SQLite file handle
    type(struct_obs),    intent(in) :: obsdat
    integer,             intent(in) :: fileNumber

    type(fSQL_STATEMENT)   :: stmt ! type for precompiled SQLite statements
    type(fSQL_STATUS)      :: stat !type for error status

    character(len = 256)   :: query
    integer                :: numberInsert, headerIndex, obsIdo, obsIdf
    integer                :: NCO2
    real                   :: ETOP,VTOP,ECF,VCF,HE,ZTSR,ZTM,ZTGM,ZLQM,ZPS
    character(len=*), parameter :: myName    = 'sqlr_addCloudParametersandEmissivity'
    character(len=*), parameter :: myWarning = '****** '// myName //': WARNING: '

    query = 'create table if not exists cld_params( id_obs integer,ETOP real,VTOP real, &
         ECF real,VCF real,HE real,ZTSR real,NCO2 integer,ZTM real,ZTGM real,ZLQM real,ZPS real);'
    query=trim(query)
    write(*,*) myName//': create query = ', trim(query)

    call fSQL_do( db, trim(query), stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat, 'fSQL_do : ')

    query = 'insert into cld_params(id_obs,ETOP,VTOP,ECF,VCF,HE,ZTSR,NCO2,ZTM,ZTGM,ZLQM,ZPS) &
         values(?,?,?,?,?,?,?,?,?,?,?,?);'
    query=trim(query)

    write(*,*) ' Insert query = ', trim(query)

    call fSQL_prepare( db, query, stmt, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat, 'fSQL_prepare : ')
    call fSQL_begin(db)
    numberInsert=0
    HEADER: do headerIndex = 1, obs_numHeader(obsdat)
      obsIdf = obs_headElem_i(obsdat, OBS_IDF, headerIndex )
      if ( obsIdf /= fileNumber ) cycle HEADER
      obsIdo   = obs_headPrimaryKey(obsdat, headerIndex)
      ETOP = obs_headElem_r(obsdat, OBS_ETOP, headerIndex )
      VTOP = obs_headElem_r(obsdat, OBS_VTOP, headerIndex )
      ECF  = obs_headElem_r(obsdat, OBS_ECF,  headerIndex )
      VCF  = obs_headElem_r(obsdat, OBS_VCF,  headerIndex )
      HE   = obs_headElem_r(obsdat, OBS_HE,   headerIndex )
      ZTSR = obs_headElem_r(obsdat, OBS_ZTSR, headerIndex )
      NCO2 = obs_headElem_i(obsdat, OBS_NCO2, headerIndex )
      ZTM  = obs_headElem_r(obsdat, OBS_ZTM,  headerIndex )
      ZTGM = obs_headElem_r(obsdat, OBS_ZTGM, headerIndex )
      ZLQM = obs_headElem_r(obsdat, OBS_ZLQM, headerIndex )
      ZPS  = obs_headElem_r(obsdat, OBS_ZPS,  headerIndex )

      call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = obsIdo )
      call fSQL_bind_param( stmt, PARAM_INDEX = 2, REAL_VAR = ETOP   )
      call fSQL_bind_param( stmt, PARAM_INDEX = 3, REAL_VAR = VTOP   )
      call fSQL_bind_param( stmt, PARAM_INDEX = 4, REAL_VAR = ECF    )
      call fSQL_bind_param( stmt, PARAM_INDEX = 5, REAL_VAR = VCF    )
      call fSQL_bind_param( stmt, PARAM_INDEX = 6, REAL_VAR = HE     )
      call fSQL_bind_param( stmt, PARAM_INDEX = 7, REAL_VAR = ZTSR   )
      call fSQL_bind_param( stmt, PARAM_INDEX = 8, INT_VAR  = NCO2   )
      call fSQL_bind_param( stmt, PARAM_INDEX = 9, REAL_VAR = ZTM    )
      call fSQL_bind_param( stmt, PARAM_INDEX = 10,REAL_VAR = ZTGM   )
      call fSQL_bind_param( stmt, PARAM_INDEX = 11,REAL_VAR = ZLQM   )
      call fSQL_bind_param( stmt, PARAM_INDEX = 12,REAL_VAR = ZPS    )

      call fSQL_exec_stmt ( stmt )
      numberInsert=numberInsert +1
    end do HEADER

    call fSQL_finalize( stmt )
    call fSQL_commit(db)
    write(*,'(a,i8)') myName//': NUMBER OF INSERTIONS ----> ', numberInsert

  end subroutine sqlr_addCloudParametersandEmissivity

  !--------------------------------------------------------------------------
  ! sqlr_insertSqlite
  !--------------------------------------------------------------------------
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
    integer                :: numberInsert, headerIndex, bodyIndex, numHeader
    integer                :: obsNlv, obsRln, obsIdf, insertItem
    integer(8)             :: bodyPrimaryKey, headPrimaryKey
    character(len = 256)   :: query
    logical                :: llok    
    character(len=*), parameter :: myName = 'sqlr_insertSqlite'
    namelist/namSQLInsert/ numberInsertItems, itemInsertList

    write(*,*)  myName//': --- Starting ---   '
    numHeader = obs_numHeader(obsdat)
    write(*,*) ' FAMILY ---> ', trim(familyType), '  headerIndex  ----> ', numHeader
    write(*,*) ' fileName -> ', trim(fileName)

    ! set default values of namelist variables
    numberInsertItems = 0
    itemInsertList(:) = 0

    nulnam = 0
    ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0 )
    read(nulnam, nml = namSQLInsert, iostat = ierr )
    if (ierr /= 0 ) call utl_abort( myName//': Error reading namelist' )
    if (mpi_myid == 0) write(*, nml = namSQLInsert )
    ierr=fclos( nulnam )

    write(*,*) ' INSERT INTO SQLITE FILE ELEMENTS :--> ',( itemInsertList(insertItem), insertItem = 1, numberInsertItems )

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
    call fSQL_prepare( db, query, stmt, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat, 'fSQL_prepare : ')

    numberInsert=0
    HEADER: do headerIndex = 1, obs_numHeader(obsdat)

      obsIdf = obs_headElem_i(obsdat, OBS_IDF, headerIndex )
      if ( obsIdf /= fileNumber ) cycle HEADER
      headPrimaryKey = obs_headPrimaryKey(obsdat, headerIndex)
      obsRln = obs_headElem_i(obsdat, OBS_RLN, headerIndex )
      obsNlv = obs_headElem_i(obsdat, OBS_NLV, headerIndex )

      BODY: do bodyIndex = obsRln, obsNlv + obsRln -1

        bodyPrimaryKey= obs_bodyPrimaryKey(obsdat, bodyIndex)
        obsVarno      = obs_bodyElem_i(obsdat, OBS_VNM , bodyIndex )
        obsFlag       = obs_bodyElem_i(obsdat, OBS_FLG , bodyIndex )
        vertCoordType = obs_bodyElem_i(obsdat, OBS_VCO , bodyIndex )
        obsValue      = obs_bodyElem_r(obsdat, OBS_VAR , bodyIndex )
        OMA           = obs_bodyElem_r(obsdat, OBS_OMA , bodyIndex )
        OMP           = obs_bodyElem_r(obsdat, OBS_OMP , bodyIndex )
        OER           = obs_bodyElem_r(obsdat, OBS_OER , bodyIndex )
        FGE           = obs_bodyElem_r(obsdat, OBS_HPHT, bodyIndex )
        PPP           = obs_bodyElem_r(obsdat, OBS_PPP , bodyIndex )
        llok = .false.

        do insertItem = 1, numberInsertItems
          if ( obsVarno == itemInsertList(insertItem) ) llok=.true.
        end do

        if ( bodyPrimaryKey == -1 ) then
          if ( llok ) then
            select case(trim(familyType))
              case( 'SF', 'SC', 'GP' )
                call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT8_VAR = headPrimaryKey)
                call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = obsVarno )
                call fSQL_bind_param( stmt, PARAM_INDEX = 3, REAL_VAR = PPP      )
                if ( obsValue == obs_missingValue_R ) then          ! sql null values
                  call fSQL_bind_param( stmt, PARAM_INDEX = 4                      )
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 4, REAL_VAR = obsValue )
                end if
                call fSQL_bind_param( stmt, PARAM_INDEX = 5, INT_VAR  = obsFlag  )
                if ( OMA == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 6                    ) 
                else 
                  call fSQL_bind_param( stmt, PARAM_INDEX = 6, REAL_VAR = OMA    ) 
                end if
                if ( OMP == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 7                    ) 
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 7, REAL_VAR = OMP    ) 
                end if
                if ( FGE == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 8                    ) 
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 8, REAL_VAR = FGE    )
                end if
                if ( OER == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 9                    ) 
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 9, REAL_VAR = OER    ) 
                end if
              case DEFAULT
                call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT8_VAR = headPrimaryKey)
                call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = obsVarno      )
                call fSQL_bind_param( stmt, PARAM_INDEX = 3, REAL_VAR = PPP           ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 4, INT_VAR  = vertCoordType )
                if ( obsValue == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 5                         )
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 5, REAL_VAR = obsValue    )
                end if 
                call fSQL_bind_param( stmt, PARAM_INDEX = 6, INT_VAR  = obsFlag       )
                if ( OMA == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 7                         ) 
                else 
                  call fSQL_bind_param( stmt, PARAM_INDEX = 7, REAL_VAR = OMA         )
                end if
                if ( OMP == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 8                         ) 
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 8, REAL_VAR = OMP         ) 
                end if
                if ( FGE == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 9                         ) 
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 9, REAL_VAR = FGE         ) 
                end if
                if ( OER == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 10                        ) 
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 10, REAL_VAR = OER        )
                end if 
            end select
            call fSQL_exec_stmt ( stmt )

            numberInsert = numberInsert + 1
          end if    ! llok
        end if

      end do BODY

    end do HEADER

    call fSQL_finalize( stmt )
    call fSQL_commit(db)
    write(*,'(3a,i8)') myName//': FAMILY ---> ' ,trim(familyType), '  NUMBER OF INSERTIONS ----> ', numberInsert

  end subroutine sqlr_insertSqlite

  !--------------------------------------------------------------------------
  ! sqlr_cleanSqlite
  !--------------------------------------------------------------------------
  subroutine sqlr_cleanSqlite(db, fileName)
    !
    ! :Purpose: Remove flagged (bit 11 set) observations in an SQLite file
    !
    implicit none

    ! arguments
    type(fSQL_DATABASE), intent(inout) :: db   ! SQLite file handle
    character(len=*),    intent(in) :: fileName

    ! locals
    character(len=*), parameter :: myError = 'sqlr_cleanSqlite: ERROR: '

    character(len = 128) :: query
    type(fSQL_STATEMENT) :: statement ! prepared statement for SQLite
    type(fSQL_STATUS)    :: status

    call fSQL_open(db, fileName, status)
    if (fSQL_error(status) /= FSQL_OK) then
      write(*,*) myError, fSQL_errmsg(status)
    end if
    ! Mark for deletion all records with bit 11 (2048) set
    query = ' delete from data where flag & 2048;'
    call fSQL_prepare(db, query, statement, status)
    if (fSQL_error(status) /= FSQL_OK) &
      call sqlr_handleError(status, 'thinning fSQL_prepare : ')
    call fSQL_begin(db)
    call fSQL_exec_stmt(statement)
    call fSQL_finalize(statement)
    call fSQL_commit(db)
    write(*,*) '  closed database -->', trim(FileName)
    call fSQL_close( db, status )
  end subroutine sqlr_cleanSqlite

  !--------------------------------------------------------------------------
  ! getObsFileName
  !--------------------------------------------------------------------------
  function getObsFileName(obsFamily, sfFileName_opt, codetype_opt) result(fileName)
    !
    ! :Purpose: Return the part of the observation file name associated
    !           with the type of observation it contains.
    !
    implicit none

    ! arguments:
    character(len=*)           :: obsFamily
    character(len=*), optional :: sfFileName_opt ! fileName acronym used for surface obs file
    integer, optional          :: codetype_opt
    character(len=20) :: fileName

    if ( obsFamily == 'TO' ) then
      if (.not. present(codetype_opt)) then
        call utl_abort('getObsFileName: codetype_opt must be specified for TO family')
      end if

      if ( codtyp_get_name( codeType_opt ) == 'radianceclear' ) then
        fileName  = 'csr'
      else if ( codtyp_get_name( codeType_opt ) == 'mhs' .or. codtyp_get_name( codeType_opt ) == 'amsub' ) then
        fileName = 'to_amsub'
      else if ( codtyp_get_name( codeType_opt ) == 'amsua' ) then
        if ( tvs_mwAllskyAssim .and. &
                   tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name( codeType_opt ))) ) then
          fileName = 'to_amsua_allsky'
        else
          fileName = 'to_amsua'
        end if
      else if ( codtyp_get_name( codeType_opt ) == 'ssmi' ) then
        fileName = 'ssmis'
      else if ( codtyp_get_name( codeType_opt ) == 'crisfsr' ) then
        fileName = 'cris'
      else
        fileName = codtyp_get_name( codeType_opt )
      end if
    else
      if (.not. present(sfFileName_opt)) then
        call utl_abort('getObsFileName: sfFileName_opt must be specified')
      end if
      call up2low( obsFamily, fileName )
      if ( fileName == 'ra' ) fileName = 'radar'
      if ( fileName == 'sf' ) then
        ! use either 'sf' or 'sfc' for filename with surface obs
        fileName = sfFileName_opt
      end if
    end if

  end function getObsFileName

  !--------------------------------------------------------------------------
  ! sqlr_writeAllSqlDiagFiles
  !--------------------------------------------------------------------------
  subroutine sqlr_writeAllSqlDiagFiles( obsdat, sfFileName, onlyAssimObs, addFSOdiag )
    !
    ! :Purpose: To prepare the writing of obsSpaceData content into SQLite format files
    !  
    implicit none

    ! arguments:
    type(struct_obs)       :: obsdat       ! obsSpaceData object
    character(len=*)       :: sfFileName   ! fileName acronym used for surface obs file
    logical                :: onlyAssimObs ! only write assimilated obs
    logical                :: addFSOdiag   ! only write FSO

    ! locals:
    integer                :: familyIndex, codeTypeIndex, fileIndex
    character(len=2)       :: obsFamilyList(50)
    integer                :: obsFamilyListSize
    integer                :: tovsAllCodeTypeListSize, tovsAllCodeTypeList(30)
    integer                :: tovsCodeTypeListSize, tovsCodeTypeList(10)
    integer                :: tovsFileNameListSize
    character(len=20)      :: tovsFileNameList(30)
    character(len=20)      :: fileName

    ! ensure all mpi tasks have same list of common obs family names
    call getObsFamilyListMpiGlobal(obsdat, obsFamilyListSize, obsFamilyList)

    ! get list of all possible tovs codetype values and unique list of corresponding filenames
    call tvs_getAllIdBurpTovs(tovsAllCodeTypeListSize, tovsAllCodeTypeList)
    write(*,*) 'tovsAllCodeTypeListSize = ', tovsAllCodeTypeListSize
    write(*,*) 'tovsAllCodeTypeList = ', tovsAllCodeTypeList(1:tovsAllCodeTypeListSize)
    
    tovsFileNameListSize = 0
    tovsFileNameList(:) = 'XXXXX'
    do codeTypeIndex = 1, tovsAllCodeTypeListSize
      fileName = getObsFileName('TO', codeType_opt=tovsAllCodeTypeList(codeTypeIndex))
      if ( all(tovsFileNameList(:) /= fileName) ) then
        tovsFileNameListSize = tovsFileNameListSize + 1
        tovsFileNameList(tovsFileNameListSize) = fileName
      end if
    end do
    write(*,*) 'tovsFileNameListSize = ', tovsFileNameListSize
    write(*,*) 'tovsFileNameList = ', tovsFileNameList(1:tovsFileNameListSize)
    
    do familyIndex = 1, obsFamilyListSize

      write(*,*) 'sqlr_writeAllSqlDiagFiles: Family = ', familyIndex, obsFamilyList(familyIndex)

      if ( obsFamilyList(familyIndex) == 'TO' ) then

        do fileIndex = 1, tovsFileNameListSize
          fileName = tovsFileNameList(fileIndex)
          write(*,*) 'tovs filename = ', fileName

          ! get list of codetypes associated with this filename
          tovsCodeTypeListSize = 0
          tovsCodeTypeList(:) = MPC_missingValue_INT
          do codeTypeIndex = 1, tovsAllCodeTypeListSize
            if (fileName == getObsFileName('TO', codeType_opt=tovsAllCodeTypeList(codeTypeIndex))) then
              tovsCodeTypeListSize = tovsCodeTypeListSize + 1
              tovsCodeTypeList(tovsCodeTypeListSize) = tovsAllCodeTypeList(codeTypeIndex)
            end if
          end do

          write(*,*) 'tovsCodeTypeListSize = ', tovsCodeTypeListSize
          write(*,*) 'tovsCodeTypeList = ', tovsCodeTypeList(1:tovsCodeTypeListSize) 
          call sqlr_writeSqlDiagFile(obsdat, 'TO', onlyAssimObs, addFSOdiag, &
                                     tovsFileNameList(fileIndex), &
                                     tovsCodeTypeList(1:tovsCodeTypeListSize)) 
        end do

      else

        fileName = getObsFileName(obsFamilyList(familyIndex), sfFileName_opt=sfFileName)
        call sqlr_writeSqlDiagFile(obsdat, obsFamilyList(familyIndex),  &
                                   onlyAssimObs, addFSOdiag, fileName) 

      end if   
      
    end do

  end subroutine sqlr_writeAllSqlDiagFiles

  !--------------------------------------------------------------------------
  ! getObsFamilyListMpiGlobal
  !--------------------------------------------------------------------------
  subroutine getObsFamilyListMpiGlobal(obsdat, obsFamilyListSizeCommon,  &
                                       obsFamilyListCommon)
    !
    ! :Purpose: Obtain a common set of obs family names over all mpi tasks
    !
    implicit none
      
    ! arguments:
    type(struct_obs) :: obsdat
    integer          :: obsFamilyListSizeCommon
    character(len=*) :: obsFamilyListCommon(:)

    ! locals:
    integer                       :: headerIndex, familyIndex, charIndex, procIndex, nsize, ierr
    integer                       :: obsFamilyListSizeMpiLocal, obsFamilyListSizeMaxMpiLocal, obsFamilyListSizeMax
    character(len=2), allocatable :: obsFamilyListMpiLocal(:)
    character(len=2), allocatable :: obsFamilyListMpiGlobal(:,:)
    character(len=2)              :: currentObsFamily
    integer, allocatable          :: intObsFamilyListMpiLocal(:,:)
    integer, allocatable          :: intObsFamilyListMpiGlobal(:,:,:)
    integer, allocatable          :: allObsFamilyListSizeMpiLocal(:)

    obsFamilyListSizeMax = size(obsFamilyListCommon)
    write(*,*) 'obsFamilyListSizeMax =', obsFamilyListSizeMax

    ! get family list for this mpi task
    obsFamilyListSizeMpiLocal = 0
    allocate(obsFamilyListMpiLocal(obsFamilyListSizeMax))
    obsFamilyListMpiLocal(:) = 'XX'
    HEADER: do headerIndex = 1, obs_numHeader( obsdat )
      currentObsFamily = obs_getFamily( obsdat, headerIndex ) 
      if ( any( obsFamilyListMpiLocal(:) == currentObsFamily ) ) cycle HEADER
      obsFamilyListSizeMpiLocal = obsFamilyListSizeMpiLocal + 1
      obsFamilyListMpiLocal( obsFamilyListSizeMpiLocal ) = currentObsFamily
      write(*,*) 'add the family: ', currentObsFamily
    end do HEADER
    write(*,*) 'obsFamilyListSizeMpiLocal =', obsFamilyListSizeMpiLocal
    write(*,*) 'obsFamilyListMpiLocal = ', obsFamilyListMpiLocal(1:obsFamilyListSizeMpiLocal)

    allocate(allObsFamilyListSizeMpiLocal(mpi_nprocs))
    call rpn_comm_allgather(obsFamilyListSizeMpiLocal,    1, 'mpi_integer',  &
                            allObsFamilyListSizeMpiLocal, 1, 'mpi_integer', 'GRID', ierr)
    call rpn_comm_allreduce(obsFamilyListSizeMpiLocal, obsFamilyListSizeMaxMpiLocal,1,'mpi_integer','mpi_max','GRID',ierr)

    ! convert local family list from characters to integers
    allocate(intObsFamilyListMpiLocal(len(currentObsFamily),obsFamilyListSizeMaxMpiLocal))
    intObsFamilyListMpiLocal(:,:)=0
    do familyIndex = 1, obsFamilyListSizeMpiLocal
      do charIndex = 1, len(currentObsFamily)
        intObsFamilyListMpiLocal(charIndex,familyIndex) =  &
             iachar(obsFamilyListMpiLocal(familyIndex)(charIndex:charIndex))
      end do
    end do

    ! communicate obs family list to all mpi tasks as integers
    allocate(intObsFamilyListMpiGlobal(len(currentObsFamily),obsFamilyListSizeMaxMpiLocal,mpi_nprocs))
    nsize = size(intObsFamilyListMpiLocal)
    call rpn_comm_allgather(intObsFamilyListMpiLocal,  nsize, 'mpi_integer',  &
                            intObsFamilyListMpiGlobal, nsize, 'mpi_integer', 'GRID', ierr)

    ! convert global family lists from integers to characters
    allocate(obsFamilyListMpiGlobal(obsFamilyListSizeMaxMpiLocal,mpi_nprocs))
    obsFamilyListMpiGlobal(:,:) = 'XX'
    do procIndex = 1, mpi_nprocs
      do familyIndex = 1, allObsFamilyListSizeMpiLocal(procIndex)
        do charIndex=1,len(currentObsFamily)
          obsFamilyListMpiGlobal(familyIndex,procIndex)(charIndex:charIndex) =  &
               achar(intObsFamilyListMpiGlobal(charIndex,familyIndex,procIndex))
        end do
      end do
      write(*,*) 'obsFamilyListMpiGlobal = ', procIndex,  &
           obsFamilyListMpiGlobal(1:allObsFamilyListSizeMpiLocal(procIndex),procIndex)
    end do

    ! construct single common list of families to be used for all mpi tasks
    obsFamilyListCommon(:) = 'YY'
    obsFamilyListSizeCommon = 0
    do procIndex = 1, mpi_nprocs
      FAMILY: do familyIndex = 1, obsFamilyListSizeMaxMpiLocal
        if (obsFamilyListMpiGlobal(familyIndex,procIndex) == 'XX') cycle FAMILY
        if (any( obsFamilyListCommon(:) == obsFamilyListMpiGlobal(familyIndex,procIndex) )) cycle FAMILY
        obsFamilyListSizeCommon = obsFamilyListSizeCommon + 1
        obsFamilyListCommon(obsFamilyListSizeCommon) = obsFamilyListMpiGlobal(familyIndex,procIndex)
      end do FAMILY
    end do
    write(*,*) 'obsFamilyListSizeCommon = ', obsFamilyListSizeCommon
    write(*,*) 'obsFamilyListCommon = ', obsFamilyListCommon(1:obsFamilyListSizeCommon)

    deallocate(allObsFamilyListSizeMpiLocal)
    deallocate(obsFamilyListMpiGlobal)
    deallocate(intObsFamilyListMpiGlobal)
    deallocate(intObsFamilyListMpiLocal)
    deallocate(obsFamilyListMpiLocal)

  end subroutine getObsFamilyListMpiGlobal

  !--------------------------------------------------------------------------
  ! sqlr_writeSqlDiagFile
  !--------------------------------------------------------------------------
  subroutine sqlr_writeSqlDiagFile( obsdat, obsFamily, onlyAssimObs, addFSOdiag, &
                                    instrumentFileName, codeTypeList_opt )
    !
    ! :Purpose: To write the obsSpaceData content into SQLite format files
    !
    implicit none

    ! arguments
    type(struct_obs)  :: obsdat
    character(len=*)  :: obsFamily
    logical           :: onlyAssimObs
    logical           :: addFSOdiag
    character(len=*)  :: instrumentFileName
    integer, optional :: codeTypeList_opt(:)

    ! locals
    type(fSQL_DATABASE)    :: db                   ! type for SQLIte  file handle
    type(fSQL_STATEMENT)   :: stmtData, stmtHeader ! type for precompiled SQLite statements
    type(fSQL_STATUS)      :: stat                 ! type for error status
    integer                :: obsVarno, obsFlag, ASS, vertCoordType, codeType, date, time, idObs, idData 
    real                   :: obsValue, OMA, OMP, OER, FGE, PPP, lon, lat, altitude
    real                   :: ensInnovStdDev, ensObsErrStdDev, zhad, fso
    integer                :: numberInsertions, numHeaders, headerIndex, bodyIndex, obsNlv, obsRln
    character(len = 512)   :: queryData, queryHeader, queryCreate 
    character(len = 12 )   :: idStation
    character(len=*), parameter :: myName = 'sqlr_writeSqlDiagFile'
    character(len=*), parameter :: myError = myName //': ERROR: '
    character(len=*), parameter :: myWarning = myName //': WARNING: '
    character(len=30)      :: fileNameExtention
    character(len=256)     :: fileName, fileNameDir
    character(len=4)       :: cmyidx, cmyidy
    logical                :: writeHeader

    ! determine initial idData,idObs to ensure unique values across mpi tasks
    call getInitialIdObsData(obsDat, obsFamily, idObs, idData, codeTypeList_opt)

    ! return if this mpi task does not have any observations of this family
    if (trim(instrumentFileName) == 'XXXXX') return

    ! check if any obs exist for this file, if not return
    numHeaders = 0
    call obs_set_current_header_list( obsdat, obsFamily )
    HEADERCOUNT: do
      headerIndex = obs_getHeaderIndex( obsdat )
      if ( headerIndex < 0 ) exit HEADERCOUNT
      if ( present( codeTypeList_opt ) ) then
        codeType  = obs_headElem_i( obsdat, OBS_ITY, headerIndex )
        if ( all( codeTypeList_opt(:) /= codeType ) ) cycle HEADERCOUNT
      end if
      numHeaders = numHeaders + 1
    end do HEADERCOUNT
    if (numHeaders == 0) return

    fileNameDir = trim(ram_getRamDiskDir())
    if ( fileNameDir == ' ' ) &
    write(*,*) myWarning//' The program may be slow creating many sqlite files in the same directory.'
    write(*,*) myWarning//' Please, use the ram disk option prior to MIDAS run!'

    if ( obs_mpiLocal( obsdat ) ) then
      write(cmyidy,'(I4.4)') ( mpi_myidy + 1 )
      write(cmyidx,'(I4.4)') ( mpi_myidx + 1 )
      fileNameExtention  = trim(cmyidx) // '_' // trim( cmyidy )
    else
      if ( mpi_myid > 0 ) return
      fileNameExtention = ' '
    end if

    fileName = trim(fileNameDir) // 'obs/dia' // trim(instrumentFileName) // '_' // trim( fileNameExtention )

    write(*,*) myName//': Creating file: ', trim(fileName)
    call fSQL_open( db, fileName, stat )
    if ( fSQL_error( stat ) /= FSQL_OK ) write(*,*) myError//'fSQL_open: ', fSQL_errmsg( stat ),' filename: '//trim(fileName)

    ! Create the tables HEADER and DATA
    if (addFSOdiag) then
       queryCreate = 'create table header (id_obs integer primary key, id_stn varchar(50), lat real, lon real, &
                   &codtyp integer, date integer, time integer, elev real); &
                   &create table data (id_data integer primary key, id_obs integer, varno integer, vcoord real, &
                   &vcoord_type integer, obsvalue real, flag integer, oma real, ompt real, oma0 real, omp real, &
                   &an_error real, fg_error real, obs_error real, sigi real, sigo real, zhad real, fso real);'
    else
       queryCreate = 'create table header (id_obs integer primary key, id_stn varchar(50), lat real, lon real, &
                   &codtyp integer, date integer, time integer, elev real); &
                   &create table data (id_data integer primary key, id_obs integer, varno integer, vcoord real, &
                   &vcoord_type integer, obsvalue real, flag integer, oma real, ompt real, oma0 real, omp real, &
                   &an_error real, fg_error real, obs_error real, sigi real, sigo real, zhad real);'
    end if
    call fSQL_do_many( db, queryCreate, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError( stat, 'fSQL_do_many with query: '//trim(queryCreate) )

    if (addFSOdiag) then
       queryData = 'insert into data (id_data, id_obs, varno, vcoord, vcoord_type, obsvalue, flag, oma, oma0, ompt, fg_error, &
    obs_error, sigi, sigo, zhad, fso) values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);'
    else
       queryData = 'insert into data (id_data, id_obs, varno, vcoord, vcoord_type, obsvalue, flag, oma, oma0, ompt, fg_error, &
    obs_error, sigi, sigo, zhad) values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);'
    end if
    queryHeader = ' insert into header (id_obs, id_stn, lat, lon, date, time, codtyp, elev ) values(?,?,?,?,?,?,?,?); '

    write(*,*) myName//': Insert query Data   = ', trim( queryData )
    write(*,*) myName//': Insert query Header = ', trim( queryHeader )

    call fSQL_begin(db)
    call fSQL_prepare( db, queryData, stmtData, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat, 'fSQL_prepare : ')
    call fSQL_prepare( db, queryHeader, stmtHeader, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat, 'fSQL_prepare : ')

    numberInsertions = 0

    call obs_set_current_header_list( obsdat, obsFamily )
    HEADER: do

      headerIndex = obs_getHeaderIndex( obsdat )
      if ( headerIndex < 0 ) exit HEADER
        
      codeType  = obs_headElem_i( obsdat, OBS_ITY, headerIndex )
      if ( present( codeTypeList_opt ) ) then
        if ( all( codeTypeList_opt(:) /= codeType ) ) cycle HEADER
      end if

      obsRln    = obs_headElem_i( obsdat, OBS_RLN, headerIndex )
      obsNlv    = obs_headElem_i( obsdat, OBS_NLV, headerIndex )
      idStation = obs_elem_c    ( obsdat, 'STID' , headerIndex ) 
      altitude  = obs_headElem_r( obsdat, OBS_ALT, headerIndex )      
      lon       = obs_headElem_r( obsdat, OBS_LON, headerIndex ) * MPC_DEGREES_PER_RADIAN_R8
      lat       = obs_headElem_r( obsdat, OBS_LAT, headerIndex ) * MPC_DEGREES_PER_RADIAN_R8
      if (  lon > 180. ) lon = lon - 360.
      date      = obs_headElem_i( obsdat, OBS_DAT, headerIndex )
      time      = obs_headElem_i( obsdat, OBS_ETM, headerIndex ) * 100.

      ! check if at least one observation will be written, if not skip this header
      if (onlyAssimObs) then
        writeHeader = .false.
        BODYCHECK: do bodyIndex = obsRln, obsNlv + obsRln -1
          OMP = obs_bodyElem_r( obsdat, OBS_OMP , bodyIndex )
          FGE = obs_bodyElem_r( obsdat, OBS_HPHT, bodyIndex )
          ASS = obs_bodyElem_i( obsdat, OBS_ASS , bodyIndex )
          if ( (ASS == obs_assimilated) .and.     &
               (OMP /= obs_missingValue_R) .and.  &
               (FGE /= obs_missingValue_R) ) writeHeader = .true.
        end do BODYCHECK
        if (.not. writeHeader) cycle HEADER
      end if

      idObs = idObs + 1
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 1, INT_VAR  = idObs     )
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 2, CHAR_VAR = idStation )
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 3, REAL_VAR = lat       ) 
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 4, REAL_VAR = lon       ) 
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 5, INT_VAR  = date      ) 
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 6, INT_VAR  = time      ) 
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 7, INT_VAR  = codeType  ) 
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 8, REAL_VAR = altitude  ) 
      call fSQL_exec_stmt ( stmtHeader )

      BODY: do bodyIndex = obsRln, obsNlv + obsRln -1
         
        obsVarno      = obs_bodyElem_i( obsdat, OBS_VNM , bodyIndex )
        obsFlag       = obs_bodyElem_i( obsdat, OBS_FLG , bodyIndex )
        vertCoordType = obs_bodyElem_i( obsdat, OBS_VCO , bodyIndex )
        obsValue      = obs_bodyElem_r( obsdat, OBS_VAR , bodyIndex )
        OMA           = obs_bodyElem_r( obsdat, OBS_OMA , bodyIndex )
        OMP           = obs_bodyElem_r( obsdat, OBS_OMP , bodyIndex )
        OER           = obs_bodyElem_r( obsdat, OBS_OER , bodyIndex )
        FGE           = obs_bodyElem_r( obsdat, OBS_HPHT, bodyIndex )
        PPP           = obs_bodyElem_r( obsdat, OBS_PPP , bodyIndex )
        ASS           = obs_bodyElem_i( obsdat, OBS_ASS , bodyIndex )

        ! skip obs if it was not assimilated
        if (onlyAssimObs) then
          if ( (ASS /= obs_assimilated) .or.     &
               (OMP == obs_missingValue_R) .or.  &
               (FGE == obs_missingValue_R) ) cycle BODY
        end if

        if ( obs_columnActive_RB(obsdat, OBS_SIGI) ) then
          ensInnovStdDev = obs_bodyElem_r(obsdat, OBS_SIGI, bodyIndex)
        else
          ensInnovStdDev = obs_missingValue_R
        end if
        if ( obs_columnActive_RB(obsdat, OBS_SIGO) ) then
          ensObsErrStdDev = obs_bodyElem_r(obsdat, OBS_SIGO, bodyIndex)
        else
          ensObsErrStdDev = obs_missingValue_R
        end if
        if ( obs_columnActive_RB(obsdat, OBS_ZHA) ) then
          zhad = obs_bodyElem_r(obsdat, OBS_ZHA, bodyIndex)
        else
          zhad = obs_missingValue_R
        end if
        if (addFSOdiag) then
          if ( obs_columnActive_RB(obsdat, OBS_FSO) ) then
            fso = obs_bodyElem_r(obsdat, OBS_FSO, bodyIndex)
          else
            fso = obs_missingValue_R
          end if
        end if

        select case( obsFamily )
          case ( 'UA', 'AI', 'SW' )
            if ( vertCoordType == 2 ) vertCoordType = 7004
          case ( 'RO' )
            vertCoordType = 7007
          case ( 'PR' )
            vertCoordType = 7006
            PPP = PPP - altitude
          case ( 'TO' )
            vertCoordType = 5042
            if( codeType == 164 .or. codeType == 181 .or. codeType == 182 ) vertCoordType = 2150
          case ( 'SF', 'SC', 'GP' )
            vertCoordType = MPC_missingValue_INT 
        end select

        ! insert order: id_obs,varno,vcoord,vcoord_type,obsvalue,flag,oma,oma0,ompt,fg_error,obs_error,sigi,sigo
        idData = idData + 1
        call fSQL_bind_param( stmtData, PARAM_INDEX = 1, INT_VAR  = idData        )
        call fSQL_bind_param( stmtData, PARAM_INDEX = 2, INT_VAR  = idObs         )
        call fSQL_bind_param( stmtData, PARAM_INDEX = 3, INT_VAR  = obsVarno      )
        call fSQL_bind_param( stmtData, PARAM_INDEX = 4, REAL_VAR = PPP           )
        if ( vertCoordType == MPC_missingValue_INT ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 5                         ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 5, INT_VAR  = vertCoordType ) 
        end if
        call fSQL_bind_param( stmtData, PARAM_INDEX = 6, REAL_VAR = obsValue      ) 
        call fSQL_bind_param( stmtData, PARAM_INDEX = 7, INT_VAR  = obsFlag       )
        if ( OMA == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 8                         ) 
          call fSQL_bind_param( stmtData, PARAM_INDEX = 9                         ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 8, REAL_VAR = OMA         )
          call fSQL_bind_param( stmtData, PARAM_INDEX = 9, REAL_VAR = OMA         )
        end if
        if ( OMP == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 10                         ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 10, REAL_VAR = OMP         )
        end if
        if ( FGE == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 11                         ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 11, REAL_VAR = FGE         )
        end if
        if ( OER == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 12                        ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 12, REAL_VAR = OER        )
        end if 
        if ( ensInnovStdDev == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 13                        ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 13, REAL_VAR = ensInnovStdDev )
        end if 
        if ( ensObsErrStdDev == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 14                        ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 14, REAL_VAR = ensObsErrStdDev )
        end if 
        if ( zhad == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 15                        ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 15, REAL_VAR = zhad )
        end if 
        if (addFSOdiag) then
          if ( fso == obs_missingValue_R ) then
            call fSQL_bind_param( stmtData, PARAM_INDEX = 16                        )
          else
            call fSQL_bind_param( stmtData, PARAM_INDEX = 16, REAL_VAR = fso )
          end if
        end if

        call fSQL_exec_stmt ( stmtData )

        numberInsertions = numberInsertions + 1

      end do BODY
     
    end do HEADER

    write(*,*) myName// ': Observation Family: ', obsFamily,', number of insertions: ', numberInsertions
    call fSQL_finalize( stmtData )
    call fSQL_commit(db)
    call fSQL_close( db, stat )

  end subroutine sqlr_writeSqlDiagFile

  !--------------------------------------------------------------------------
  ! getInitialIdObsData
  !--------------------------------------------------------------------------
  subroutine getInitialIdObsData(obsDat, obsFamily, idObs, idData, codeTypeList_opt)
    !
    ! :Purpose: Compute initial value for idObs and idData that will ensure
    !           unique values over all mpi tasks
    !
    implicit none

    ! arguments:
    type(struct_obs)  :: obsdat
    character(len=*)  :: obsFamily    
    integer           :: idObs, idData
    integer, optional :: codeTypeList_opt(:)

    ! locals:
    integer                :: headerIndex, numHeader, numBody, codeType, ierr
    integer, allocatable   :: allNumHeader(:), allNumBody(:)

    numHeader = 0
    numBody = 0
    call obs_set_current_header_list( obsdat, obsFamily )
    HEADERCOUNT: do
      headerIndex = obs_getHeaderIndex( obsdat )
      if ( headerIndex < 0 ) exit HEADERCOUNT
      if ( present( codeTypeList_opt ) ) then
        codeType  = obs_headElem_i( obsdat, OBS_ITY, headerIndex )
        if ( all( codeTypeList_opt(:) /= codeType ) ) cycle HEADERCOUNT
      end if
      numHeader = numHeader + 1
      numBody = numBody + obs_headElem_i( obsdat, OBS_NLV, headerIndex )
    end do HEADERCOUNT
    allocate(allNumHeader(mpi_nprocs))
    allocate(allNumBody(mpi_nprocs))
    call rpn_comm_allgather(numHeader,1,'mpi_integer',       &
                            allNumHeader,1,'mpi_integer','GRID',ierr)
    call rpn_comm_allgather(numBody,1,'mpi_integer',       &
                            allNumBody,1,'mpi_integer','GRID',ierr)
    if (mpi_myid > 0) then
      idObs = sum(allNumHeader(1:mpi_myid))
      idData = sum(allNumBody(1:mpi_myid))
    else
      idObs = 0
      idData = 0
    end if
    deallocate(allNumHeader)
    deallocate(allNumBody)
  end subroutine getInitialIdObsData

end module sqliteRead_mod


module sqliteRead_mod
  ! MODULE sqliteRead_mod (prefix='sqlr' category='3. Observation input/output')
  !
  !:Purpose:  To read and update SQLITE observation files. Data is stored in 
  !           obsSpaceData object.
  !
  use codePrecision_mod
  use obsSpaceData_mod
  use midasMpi_mod
  use fSQLite
  use mathPhysConstants_mod
  use obsUtil_mod
  use utilities_mod
  use bufr_mod
  use ramDisk_mod
  use codtyp_mod
  use obsVariableTransforms_mod
  use obsFilter_mod
  use sqliteUtilities_mod
  use radialVelocity_mod

  implicit none

  save

  private

  public :: sqlr_insertSqlite, sqlr_updateSqlite, sqlr_readSqlite
  public :: sqlr_cleanSqlite, sqlr_readSqlite_avhrr, sqlr_addCloudParametersandEmissivity
  public :: sqlr_writePseudoSSTobs, sqlr_writeEmptyPseudoSSTobsFile
  public :: sqlr_getColumnValuesDate

  contains
  
  !--------------------------------------------------------------------------
  ! sqlr_readSqlite_avhrr
  !--------------------------------------------------------------------------
  subroutine sqlr_readSqlite_avhrr(obsdat, fileName, headerIndexBegin, headerIndexEnd)
    !
    ! :Purpose: To read SQLite avhrr_cloud parameters.
    ! 
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsdat     ! ObsSpaceData Structure
    character(len=*) , intent(in)    :: fileName   ! SQLite filename
    integer          , intent(in)    :: headerIndexBegin
    integer          , intent(in)    :: headerIndexEnd

    ! Locals:
    type(fSQL_DATABASE)      :: db   ! type for SQLIte  file handle
    type(fSQL_STATEMENT)     :: stmt ! type for precompiled SQLite statements
    type(fSQL_STATUS)        :: stat ! type for error status
    integer                  :: obsIdo
    character(len=128)       :: querySqlite
    integer                  :: rowIndex, headerIndex, columnIndex
    integer                  :: numberAvhrrRows ,  numberAvhrrColumns
    real, allocatable        :: matdata(:,:)
    REAL(pre_obsReal) :: CFRAC,MOYRAD,STDRAD

    write(*,*) 'sqlr_readSqlite_avhrr: fileName : ', trim(fileName)

    if (.not. sqlu_sqlTableExists(trim(fileName), 'avhrr')) then
      write(*,*) 'sqlr_readSqlite_avhrr: Table avhrr does not exist :  ... return  '
      return
    end if
    write(*,*) 'sqlr_readSqlite_avhrr: Table avhrr exists: insert contents into obsdat '

    call fSQL_open(db, trim(fileName) ,stat)
    if (fSQL_error(stat) /= FSQL_OK) then
      call utl_abort('sqlr_readSqlite_avhrr: fSQL_open '//fSQL_errmsg(stat))
    end if

    querySqlite = ' select mean_radiance,stddev_radiance,fractionClearPixels from avhrr where id_obs = ? '
    call fSQL_prepare(db, querySqlite , stmt, stat)
    write(*,*) 'sqlr_readSqlite_avhrr: obs_getNchanAvhr=',obs_getNchanAvhrr()
    do headerIndex = headerIndexBegin, headerIndexEnd
      obsIdo = obs_headPrimaryKey(obsdat, headerIndex)
      call fSQL_bind_param(stmt, param_index = 1, int_var = obsIdo)
      call fSQL_exec_stmt(stmt)
      call fSQL_get_many(stmt, nrows = numberAvhrrRows , ncols = numberAvhrrColumns , mode = FSQL_REAL)
      allocate(matdata(numberAvhrrRows, numberAvhrrColumns))
      matdata(:,:) = MPC_missingValue_R4
      call fSQL_fill_matrix(stmt, matdata)

      rowIndex = 1
      do columnIndex = OBS_CF1, OBS_CF7
        if(obs_columnActive_RH(obsdat,columnIndex)) then
          CFRAC = matdata(rowIndex,3)
          call obs_headSet_r(obsdat,columnIndex,headerIndex, CFRAC)
          rowIndex = rowIndex + 6
        end if
      end do

      rowIndex = 1
      do columnIndex = OBS_M1C1, OBS_M7C6
        if(obs_columnActive_RH(obsdat,columnIndex)) then
          MOYRAD = matdata(rowIndex,1) * 100000.d0
          call obs_headSet_r(obsdat,columnIndex,headerIndex, MOYRAD)
          rowIndex = rowIndex + 1
        endif
      end do

      rowIndex = 1
      do columnIndex = OBS_S1C1, OBS_S7C6
        if(obs_columnActive_RH(obsdat,columnIndex)) then
          STDRAD = matdata(rowIndex,2) * 100000.d0
          call obs_headSet_r(obsdat,columnIndex,headerIndex,  STDRAD)
          rowIndex = rowIndex + 1
        end if
      end do

      deallocate(matdata)
      call fSQL_free_mem(stmt)

    end do

    call fSQL_finalize(stmt)
    call fSQL_close(db, stat)

  end subroutine sqlr_readSqlite_avhrr

  !--------------------------------------------------------------------------
  ! sqlr_readSqlite
  !--------------------------------------------------------------------------
  subroutine sqlr_readSqlite(obsdat, familyType, fileName)
    !
    ! :Purpose: To read SQLite namelist and files.
    ! 
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsdat     ! ObsSpaceData Structure
    character(len=*) , intent(in)    :: familyType ! Family Type 
    character(len=*) , intent(in)    :: fileName   ! SQLite filename

    ! Locals:
    character(len=12)        :: idStation
    integer                  :: codeType, obsDate, obsTime, obsStatus, obsFlag, obsVarno
    integer(8)               :: headPrimaryKey
    integer(8), allocatable  :: bodyHeadKeys(:), bodyPrimaryKeys(:)
    integer(8), allocatable  :: tempBodyKeys(:,:)
    logical                  :: finishedWithHeader
    real(pre_obsReal)        :: elevReal, xlat, xlon, vertCoord
    real                     :: elev    , obsLat, obsLon, elevFact
    real                     :: beamAzimuth, beamRangeStart, beamRangeEnd, beamElevation          
    real(pre_obsReal)        :: beamAzimuthReal, beamElevationReal
    real(pre_obsReal)        :: beamLat, beamLon, beamHeight, beamDistance, beamRange
    integer                  :: vertCoordType, vertCoordFact, fnom, fclos, nulnam, ierr, idProf
    real                     :: zenithReal, solarZenithReal, CloudCoverReal, solarAzimuthReal
    integer                  :: roQcFlag
    real(pre_obsReal)        :: geoidUndulation, earthLocRadCurv, obsValue, obsError
    real(8)                  :: geoidUndulation_R8, earthLocRadCurv_R8, azimuthReal_R8
    integer                  :: trackCellNum, iceChartID
    real(pre_obsReal)        :: modelWindSpeed
    real(8)                  :: modelWindSpeed_R8
    integer                  :: iasiImagerCollocationFlag, iasiGeneralQualityFlag
    integer                  :: obsSat, landSea, terrainType, instrument, sensor
    integer                  :: rowIndex, obsNlv, headerIndex, bodyIndex
    integer                  :: numBody, numHeader
    real(pre_obsReal), parameter :: zemFact = 0.01
    character(len=1024)      :: queryHeader, queryIDs
    character(len=256)       :: csqlcrit, selectIDs
    character(len=1024)      :: columnsHeader
    logical                  :: finished
    type(fSQL_DATABASE)      :: db         ! type for SQLIte  file handle
    type(fSQL_STATEMENT)     :: stmt,stmt2 ! type for precompiled SQLite statements
    type(fSQL_STATUS)        :: stat       !type for error status
    character(len=256)       :: sqlDataOrder, extraQueryData
    character(len=256), allocatable :: listElemArray(:)
    integer, allocatable            :: listElemArrayInteger(:)
    integer                  :: numberBodyRows, numberBodyColumns, numberIDsRows, numberIDsColumns
    integer                  :: columnIndex
    integer, parameter :: lenSqlName    = 60
    character(len=lenSqlName), allocatable :: headSqlNames(:), bodySqlNames(:)
    character(len=lenSqlName), parameter :: headSqlNamesToRead(32) = (/'ID_STN','LAT','LON','CODTYP', &
         'DATE','TIME','STATUS','ELEV','ANTENNA_ALTITUDE','LAND_SEA','ID_SAT','INSTRUMENT','SENSOR', &
         'ZENITH','SOLAR_ZENITH','AZIMUTH','SOLAR_AZIMUTH','TERRAIN_TYPE','CLOUD_COVER', &
         'FANION_QUAL_IASI_SYS_IND','INDIC_NDX_QUAL_GEOM','RO_QC_FLAG','GEOID_UNDULATION', &
         'EARTH_LOCAL_RAD_CURV','CENTER_AZIMUTH','CENTER_ELEVATION','RANGE_START','RANGE_END', &
         'ID_PROF','CHARTINDEX','TRACK_CELL_NO','MOD_WIND_SPD'/)
    real(8),                   allocatable :: bodyValues(:,:), codtypInFileList(:,:)
    logical :: beamRangeFound
    real(pre_obsReal)           :: missingValue

    ! Namelist variables:
    integer                  :: numberElem        ! MUST NOT BE INCLUDED IN NAMELIST!
    character(len=256)       :: listElem          ! list of bufr element ids to read
    character(len=256)       :: sqlExtraDat       ! can be used e.g. for ' and id_obs in (select id_obs from header where...) '
    character(len=256)       :: sqlExtraDat_sat   ! same as sqlExtraDat, but only for SST satellite obs
    character(len=256)       :: sqlExtraHeader    ! can be used e.g. for ' id_stn in (...) '
    character(len=256)       :: sqlExtraHeader_sat! same as sqlExtraHeader, but only for SST satellite obs
    integer                  :: codtyp_sat(10)    ! list of codtyp values of SST satellite obs (default is 88)
    character(len=256)       :: sqlNull           ! can be used e.g. for ' and obsvalue is not null '

    namelist /NAMSQLtovs/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull
    namelist /NAMSQLua/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull
    namelist /NAMSQLai/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull
    namelist /NAMSQLsw/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull
    namelist /NAMSQLro/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull
    namelist /NAMSQLsfc/  numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull
    namelist /NAMSQLsc/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull
    namelist /NAMSQLpr/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull
    namelist /NAMSQLal/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull
    namelist /NAMSQLgl/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull
    namelist /NAMSQLradar/numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull
    namelist /NAMSQLsst/  numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull, &
                          sqlExtraDat_sat,sqlExtraHeader_sat,codtyp_sat

    missingValue = real(MPC_missingValue_R8,pre_obsReal)

    write(*,*) 'sqlr_readSqlite: fileName   : ', trim(fileName)
    write(*,*) 'sqlr_readSqlite: familyType : ', trim(familyType)

    ! Get the codtyp from the file to use in determining how to build the main queries
    call sqlu_getColumnValuesNum(codtypInFileList, fileName=trim(fileName), &
                                 tableName='header', sqlColumnNames=(/'codtyp'/), &
                                 extraQuery_opt='group by codtyp')
    write(*,*) 'sqlr_readSqlite: codtyp array = ', codtypInFileList(:,:)

    ! Default namelist variable values
    sqlExtraHeader     = ''
    sqlExtraDat        = ''
    sqlExtraHeader_sat = ''
    sqlExtraDat_sat    = ''
    codtyp_sat(:)      = MPC_missingValue_INT
    codtyp_sat(1)      = 88
    sqlNull            = ''
    listElem           = ''
    numberElem         = MPC_missingValue_INT
    
    ! Set the type of vertical coordinate
    vertCoordType  = 1
    select case(trim(familyType))
      case('RA','PR','AL','RO','SF','TM','SC','GL','HY','GP')
        vertCoordType = 1
      case('UA','AI','SW')
        vertCoordType = 2
      case('TO')
        vertCoordType = 3
      case default
        call utl_abort('sqlr_readSqlite: unknown family '//trim(familyType))
    end select

    ! Set multiplying factor for vertical coordinate
    select case(trim(familyType))
      case('SF','TM','SC')
        vertCoordFact = 0
      case default
        vertCoordFact = 1
    end select

    ! Set multiplying factor used when adding elevation to vcoord
    select case(trim(familyType))
      case('PR','SF','TM','GP')
        elevFact=1.0
      case default
        elevFact=0.0
    end select

    ! Read appropriate namelist based on obs family
    nulnam = 0
    ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    select case(trim(familyType))
      case('TO')
        read(nulnam, nml = NAMSQLtovs, iostat = ierr)
        if (ierr /= 0) call utl_abort('sqlr_readSqlite: Error reading namelist: NAMSQLtovs')
        if (mmpi_myid == 0) write(*, nml = NAMSQLtovs)
      case('UA')
        read(nulnam, nml = NAMSQLua, iostat = ierr)
        if (ierr /= 0) call utl_abort('sqlr_readSqlite: Error reading namelist: NAMSQLua')
        if (mmpi_myid == 0) write(*, nml = NAMSQLua)
      case ('AI')
        read(nulnam, nml = NAMSQLai, iostat = ierr)
        if (ierr /= 0) call utl_abort('sqlr_readSqlite: Error reading namelist: NAMSQLai')
        if (mmpi_myid == 0) write(*, nml = NAMSQLai)
      case ('SW')
        read(nulnam, nml = NAMSQLsw, iostat = ierr)
        if (ierr /= 0) call utl_abort('sqlr_readSqlite: Error reading namelist: NAMSQLsw')
        if (mmpi_myid == 0) write(*, nml = NAMSQLsw)
      case ('PR')
        read(nulnam, nml = NAMSQLpr, iostat = ierr)
        if (ierr /= 0) call utl_abort('sqlr_readSqlite: Error reading namelist: NAMSQLpr')
        if (mmpi_myid == 0) write(*, nml =  NAMSQLpr)
      case ('AL')  
        read(nulnam, nml = NAMSQLal, iostat = ierr)
        if (ierr /= 0) call utl_abort('sqlr_readSqlite: Error reading namelist: NAMSQLal')
        if (mmpi_myid == 0) write(*, nml =  NAMSQLal)
      case ('RO')     
        read(nulnam, nml = NAMSQLro, iostat = ierr)
        if (ierr /= 0) call utl_abort('sqlr_readSqlite: Error reading namelist: NAMSQLro')
        if (mmpi_myid == 0) write(*, nml = NAMSQLro)
      case ('SF','GP','HY')
        read(nulnam, nml = NAMSQLsfc, iostat = ierr)
        if (ierr /= 0) call utl_abort('sqlr_readSqlite: Error reading namelist: NAMSQLsfc')
        if (mmpi_myid == 0) write(*, nml = NAMSQLsfc)        
      case ('SC')
        read(nulnam, nml = NAMSQLsc, iostat = ierr)
        if (ierr /= 0) call utl_abort('sqlr_readSqlite: Error reading namelist: NAMSQLsc')
        if (mmpi_myid == 0) write(*, nml = NAMSQLsc)
      case ('TM')
        read(nulnam, nml = NAMSQLsst, iostat = ierr)
        if (ierr /= 0) call utl_abort('sqlr_readSqlite: Error reading namelist: NAMSQLsst')
        if (mmpi_myid == 0) write(*, nml = NAMSQLsst)
        do rowIndex = 1, size(codtypInFileList,1)
          if (any(nint(codtypInFileList(rowIndex,1)) == codtyp_sat(:))) then
            write(*,*) 'sqlr_readSqlite: Found satellite SST observation in file:', &
                       nint(codtypInFileList(rowIndex,1))
            sqlExtraHeader = sqlExtraHeader_sat
            sqlExtraDat = sqlExtraDat_sat
          end if
        end do
      case ('GL')
        read(nulnam, nml = NAMSQLgl, iostat = ierr)
        if (ierr /= 0) call utl_abort('sqlr_readSqlite: Error reading namelist: NAMSQLgl')
        if (mmpi_myid == 0) write(*, nml =  NAMSQLgl)
      case('RA')
        read(nulnam, nml = NAMSQLradar, iostat = ierr)
        if (ierr /= 0) call utl_abort('sqlr_readSqlite: Error reading namelist: NAMSQLradar')
        if (mmpi_myid == 0) write(*, nml =  NAMSQLradar) 
      case default
        call utl_abort('sqlr_readSqlite: No namelist read for this family: '//trim(familyType))
    end select
    ierr=fclos(nulnam)
    
    if (numberElem /= MPC_missingValue_R4) then
      call utl_abort('sqlr_readSqlite: check namelist, numberElem should be removed')
    end if
    numberElem = count(transfer(listElem, 'a', len(listElem)) == ',') + 1

    ! Set the observation variable transforms
    if (numberElem > 0) then
      call utl_splitString(listElem,',', listElemArray)
      call utl_stringArrayToIntegerArray(listElemArray, listElemArrayInteger)
      call ovt_setup(listElemArrayInteger)
      deallocate(listElemArrayInteger)
      deallocate(listElemArray)
    end if

    ! Compose SQL queries

    ! Get list of all header table columns in file, then remove those we don't read
    call sqlu_getSqlColumnNames(headSqlNames, fileName=trim(fileName), &
                                tableName='header', dataType='all')
    do columnIndex = 1, size(headSqlNames)
      if (utl_findloc(headSqlNamesToRead,headSqlNames(columnIndex)) == 0) then
        headSqlNames(columnIndex) = ''
      end if
    end do
    call utl_removeEmptyStrings(headSqlNames)
    do columnIndex = 1, size(headSqlNames)
      write(*,*) 'sqlr_readSqlite: headSqlNames   =', columnIndex, &
                 trim(headSqlNames(columnIndex))
    end do
    call utl_combineString(columnsHeader, ',', headSqlNames)

    ! ordering of data that will get read in matdata, bodyPrimaryKeys and bodyHeadKeys
    sqlDataOrder = ' order by id_obs, varno, id_data' 

    selectIDs  = 'select id_data, id_obs'
    csqlcrit = 'varno in ('//trim(listElem)//')'//trim(sqlExtraDat)//trim(SQLNull)

    ! Determine names of body columns present in sqlite file
    call sqlu_getSqlColumnNames(bodySqlNames, fileName=trim(fileName), &
                                tableName='data', dataType='numeric')
    do columnIndex = 1, size(bodySqlNames)
      write(*,*) 'sqlr_readSqlite: bodySqlNames   =', columnIndex, &
                 trim(bodySqlNames(columnIndex))
    end do

    ! Read all of the columns in the file, then we decide what to do with them later
    extraQueryData = 'where '//trim(csqlcrit)//trim(sqlDataOrder)
    call sqlu_getColumnValuesNum(bodyValues, fileName=trim(fileName), &
                                 tableName='data', sqlColumnNames=bodySqlNames, &
                                 extraQuery_opt=extraQueryData)
    numberBodyRows    = size(bodyValues,1)
    numberBodyColumns = size(bodyValues,2)

    ! It is very important that queryIDs and queryData be identical except for the column names being selected
    queryIDs  = trim(selectIDs) //' from data where '//trim(csqlcrit)//trim(sqlDataOrder)//';'
    write(*,'(4a)') 'sqlr_readSqlite: queryIDs     --> ', trim(queryIDs)

    ! Open the sqlite file
    call fSQL_open(db, trim(fileName) ,stat)
    if (fSQL_error(stat) /= FSQL_OK) then
      call utl_abort('sqlr_readSqlite: fSQL_open '//fSQL_errmsg(stat))
    end if

    ! Read id_data and id_obs columns in the body table ("status" not set when getting integers)
    call fSQL_prepare(db, trim(queryIDs) , stmt2, status=stat)
    call fSQL_get_many(stmt2, nrows=numberIDsRows, ncols=numberIDsColumns, &
                       mode=FSQL_INT8)
    allocate(tempBodyKeys(numberIDsRows,numberIDsColumns))
    call fSQL_fill_matrix(stmt2, tempBodyKeys)
    allocate(bodyPrimaryKeys(numberIDsRows))
    allocate(bodyHeadKeys(numberIDsRows))
    bodyPrimaryKeys(:) = tempBodyKeys(:,1)
    bodyHeadKeys(:)    = tempBodyKeys(:,2)
    deallocate(tempBodyKeys)
    call fSQL_free_mem(stmt2)
    call fSQL_finalize(stmt2)

    if (numberIDsRows /= numberBodyRows) then
      call utl_abort('sqlr_readSqlite: number of body keys not equal to number of rows in bodyValues')
    end if

    headerIndex  = obs_numHeader(obsdat)
    bodyIndex = obs_numBody(obsdat)
    numHeader = obs_numHeader(obsdat)
    numBody   = obs_numBody(obsdat)
    write(*,*) 'sqlr_readSqlite: DEBUT numheader  =', numHeader
    write(*,*) 'sqlr_readSqlite: DEBUT numbody    =', numBody
    write(*,*) 'sqlr_readSqlite:  numberBodyRows numberBodyColumns =', numberBodyRows, numberBodyColumns
    write(*,*) 'sqlr_readSqlite: =========================================='

    ! Here is a summary of what is going on in the rest of this routine:
    ! 
    ! for each row of the body table (which has already been read)
    !
    !   if this is the first id_data associated with a given id_obs
    !     -> read associated header and write to obsSpaceData
    !
    !   -> write body entry in obsSpaceData
    !
    obsNlv = 0
    call fSQL_begin(db)
    BODY: do rowIndex = 1, numberIDsRows

      ! id_obs, bodyIndex and obsNlv values for this entry
      headPrimaryKey = bodyHeadKeys(rowIndex)
      bodyIndex = bodyIndex + 1
      obsNlv    = obsNlv + 1
      call obs_setBodyPrimaryKey(obsdat, bodyIndex, bodyPrimaryKeys(rowIndex))

      if (obs_columnActive_RB(obsdat, OBS_OMA))  call obs_bodySet_r(obsdat, OBS_OMA , bodyIndex, missingValue)
      if (obs_columnActive_RB(obsdat, OBS_OMA0)) call obs_bodySet_r(obsdat, OBS_OMA0, bodyIndex, missingValue)
      if (obs_columnActive_RB(obsdat, OBS_OMP))  call obs_bodySet_r(obsdat, OBS_OMP , bodyIndex, missingValue)
      if (obs_columnActive_RB(obsdat, OBS_OMP6)) call obs_bodySet_r(obsdat, OBS_OMP6, bodyIndex, missingValue)
      if (obs_columnActive_RB(obsdat, OBS_OER))  call obs_bodySet_r(obsdat, OBS_OER , bodyIndex, missingValue)
      if (obs_columnActive_RB(obsdat, OBS_HPHT)) call obs_bodySet_r(obsdat, OBS_HPHT, bodyIndex, missingValue)
      if (obs_columnActive_RB(obsdat, OBS_HAHT)) call obs_bodySet_r(obsdat, OBS_HAHT, bodyIndex, missingValue)
      if (obs_columnActive_RB(obsdat, OBS_WORK)) call obs_bodySet_r(obsdat, OBS_WORK, bodyIndex, missingValue)
      if (obs_columnActive_RB(obsdat, OBS_SIGI)) call obs_bodySet_r(obsdat, OBS_SIGI, bodyIndex, missingValue)
      if (obs_columnActive_RB(obsdat, OBS_SIGO)) call obs_bodySet_r(obsdat, OBS_SIGO, bodyIndex, missingValue)
      if (obs_columnActive_RB(obsdat, OBS_ZHA))  call obs_bodySet_r(obsdat, OBS_ZHA , bodyIndex, missingValue)
      if (obs_columnActive_RB(obsdat, OBS_BCOR)) call obs_bodySet_r(obsdat, OBS_BCOR, bodyIndex, missingValue)

      READHEADER: if (obsNlv == 1) then
        headerIndex = headerIndex + 1

        ! This is the first body entry associated with a new headPrimaryKey
        call obs_headSet_i(obsdat, OBS_RLN, headerIndex, bodyIndex)

        ! we read the associated header
        if (trim(sqlExtraHeader) == '') then
          queryHeader = 'select '//trim(columnsHeader)//' from header where id_obs = ? '
        else
          queryHeader = 'select '//trim(columnsHeader)//' from header where '//trim(sqlExtraHeader)//' and id_obs = ? '
        end if

        if (rowIndex == 1) then
          write(*,'(4a)') 'sqlr_readSqlite: first queryHeader    --> ', trim(queryHeader)
        end if
        call fSQL_prepare(db, queryHeader, stmt, stat)
        if (fSQL_error(stat)  /= FSQL_OK) then
          write(*,*) 'sqlr_readSqlite: problem reading header entry. Query: ', trim(queryHeader)
          call sqlu_handleError(stat, 'fSQL_prepare')
        end if

        call fSQL_bind_param(stmt, param_index = 1, int8_var = headPrimaryKey)

        call fSQL_get_row(stmt, finished)
        if (finished) then
          write(*,*) 'sqlr_readSqlite: problem reading header entry. Query:', trim(queryHeader)
          call utl_abort('Problem with fSQL_get_row()')
        end if

        ! The query result is inserted into variables
        terrainType = MPC_missingValue_INT
        landSea     = MPC_missingValue_INT
        sensor      = MPC_missingValue_INT
        cloudCoverReal = MPC_missingValue_R4
        elev = 0.
        elevReal = 0.
        solarAzimuthReal = MPC_missingValue_R4
        solarZenithReal  = MPC_missingValue_R4
        zenithReal = MPC_missingValue_R4
        instrument = MPC_missingValue_INT
        azimuthReal_R8 = MPC_missingValue_R8
        beamAzimuth    = MPC_missingValue_R4
        beamElevation  = MPC_missingValue_R4
        beamRangeStart = MPC_missingValue_R4
        beamRangeEnd   = MPC_missingValue_R4
        obsSat = MPC_missingValue_INT

        ! Set some default values
        if (trim(familyType) == 'TO') then
          call obs_headSet_i(obsdat, OBS_TTYP, headerIndex, terrainType)
        end if
        call obs_headSet_r(obsdat, OBS_ALT, headerIndex, elevReal)

        call obs_setFamily(obsdat, trim(familyType), headerIndex)
        call obs_setHeadPrimaryKey(obsdat, headerIndex, headPrimaryKey)
        call obs_headSet_i(obsdat, OBS_ONM, headerIndex, headerIndex)
        do columnIndex = 1, size(headSqlNames)
          select case(trim(headSqlNames(columnIndex)))
            case('LAT')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real_var  = obsLat)
              xlat = obsLat * MPC_RADIANS_PER_DEGREE_R8
              call obs_headSet_r(obsdat, OBS_LAT, headerIndex, xLat)
            case('LON')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real_var  = obsLon)
              if (obsLon < 0.) obsLon = obsLon + 360.
              xlon = obsLon * MPC_RADIANS_PER_DEGREE_R8
              call obs_headSet_r(obsdat, OBS_LON, headerIndex, xLon)
            case('CODTYP')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = codeType)
              call obs_headSet_i(obsdat, OBS_ITY, headerIndex, codeType)
            case('DATE')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = obsDate)
              call obs_headSet_i(obsdat, OBS_DAT, headerIndex, obsDate)
            case('TIME')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = obsTime)
              obsTime = obsTime/100
              call obs_headSet_i(obsdat, OBS_ETM, headerIndex, obsTime)
            case('ID_STN')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, char_var  = idStation)
              call obs_set_c(obsdat, 'STID' , headerIndex, trim(idStation))
            case('STATUS')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = obsStatus)
              call obs_headSet_i(obsdat, OBS_ST1, headerIndex, obsStatus)
            case('ELEV','ANTENNA_ALTITUDE')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real_var  = elev, REAL_MISSING=MPC_missingValue_R4)
              elevReal=elev
              if (trim(familyType) /= 'TO') then
                call obs_headSet_r(obsdat, OBS_ALT ,headerIndex, elevReal)
              end if
            case('LAND_SEA')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = landSea, INT_MISSING=MPC_missingValue_INT)
              call obs_headSet_i(obsdat, OBS_STYP, headerIndex, landSea)
            case('ID_SAT')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = obsSat, INT_MISSING=MPC_missingValue_INT)
              call obs_headSet_i(obsdat, OBS_SAT , headerIndex, obsSat)
            case('INSTRUMENT')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = instrument, INT_MISSING=MPC_missingValue_INT)
            case('ZENITH')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real_var  = zenithReal, REAL_MISSING=MPC_missingValue_R4)
              call obs_headSet_r(obsdat, OBS_SZA , headerIndex, real(zenithReal,kind=pre_obsReal))
            case('SOLAR_ZENITH')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real_var  = solarZenithReal, REAL_MISSING=MPC_missingValue_R4)
              call obs_headSet_r(obsdat, OBS_SUN , headerIndex, real(solarZenithReal,kind=pre_obsReal))
            case('AZIMUTH')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real8_var = azimuthReal_R8)
              call obs_headSet_r(obsdat, OBS_AZA , headerIndex, real(azimuthReal_R8,kind=pre_obsReal))
            case('TERRAIN_TYPE')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = terrainType, INT_MISSING=MPC_missingValue_INT)
              call obs_headSet_i(obsdat, OBS_TTYP, headerIndex, terrainType)
            case('CLOUD_COVER')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real_var  = cloudCoverReal, REAL_MISSING=MPC_missingValue_R4)
              call obs_headSet_r(obsdat, OBS_CLF , headerIndex, real(cloudCoverReal,kind=pre_obsReal))
            case('SOLAR_AZIMUTH')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real_var  = solarAzimuthReal, REAL_MISSING=MPC_missingValue_R4)
              call obs_headSet_r(obsdat, OBS_SAZ , headerIndex, real(solarAzimuthReal,kind=pre_obsReal))
            case('SENSOR')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = sensor, INT_MISSING=MPC_missingValue_INT)
            case('FANION_QUAL_IASI_SYS_IND')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = iasiGeneralQualityFlag, INT_MISSING=MPC_missingValue_INT)
              if (obs_columnActive_IH(obsdat,OBS_GQF)) then
                call obs_headSet_i(obsdat,OBS_GQF, headerIndex, iasiGeneralQualityFlag)
              end if
            case('INDIC_NDX_QUAL_GEOM')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = iasiImagerCollocationFlag, INT_MISSING=MPC_missingValue_INT)
              if (obs_columnActive_IH(obsdat,OBS_GQL)) then
                call obs_headSet_i(obsdat, OBS_GQL, headerIndex, iasiImagerCollocationFlag)
              end if
            case('RO_QC_FLAG')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = roQcFlag, INT_MISSING=MPC_missingValue_INT)
              call obs_headSet_i(obsdat, OBS_ROQF, headerIndex, roQcFlag)
            case('GEOID_UNDULATION')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real8_var = geoidUndulation_R8)
              geoidUndulation = geoidUndulation_R8
              call obs_headSet_r(obsdat, OBS_GEOI, headerIndex, geoidUndulation)
            case('EARTH_LOCAL_RAD_CURV')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real8_var = earthLocRadCurv_R8)
              earthLocRadCurv = earthLocRadCurv_R8
              call obs_headSet_r(obsdat, OBS_TRAD, headerIndex, earthLocRadCurv)
            case('CENTER_AZIMUTH')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real_var  = beamAzimuth)
              beamAzimuthReal   = beamAzimuth   * MPC_RADIANS_PER_DEGREE_R8
              call obs_headSet_r(obsdat, OBS_RZAM, headerIndex, beamAzimuthReal)
            case('CENTER_ELEVATION')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real_var  = beamElevation)
              beamElevationReal = beamElevation * MPC_RADIANS_PER_DEGREE_R8
              call obs_headSet_r(obsdat, OBS_RELE, headerIndex, beamElevationReal)
            case('RANGE_START')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real_var  = beamRangeStart)
              call obs_headSet_r(obsdat, OBS_RANS, headerIndex, real(beamRangeStart,kind=pre_obsReal))
            case('RANGE_END')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real_var  = beamRangeEnd)
              call obs_headSet_r(obsdat, OBS_RANE, headerIndex, real(beamRangeEnd,kind=pre_obsReal))
            case('ID_PROF')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = idProf)
              call obs_headSet_i(obsdat, OBS_PRFL, headerIndex, idProf)
            case('CHARTINDEX')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = iceChartID)
              call obs_headSet_i(obsdat, OBS_CHID, headerIndex, iceChartID)
            case('TRACK_CELL_NO')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, int_var   = trackCellNum)
              if (trackCellNum > 21) trackCellNum = 43 - trackCellNum
              call obs_headSet_i(obsdat, OBS_FOV, headerIndex, trackCellNum)
            case('MOD_WIND_SPD')
              call fSQL_get_column(stmt, COL_INDEX = columnIndex, real8_var = modelWindSpeed_R8)
              modelWindSpeed = modelWindSpeed_R8
              call obs_headSet_r(obsdat, OBS_MWS, headerIndex, modelWindSpeed)
          end select
        end do

        ! we are done reading this header entry
        call fSQL_finalize(stmt)

        ! adjust some values that were just read
        if (trim(familyType) == 'TO') then
          if (instrument == 420) obsSat  = 784
          if (codeType == codtyp_get_codtyp('crisfsr') .and. instrument == 620) instrument = 2046
          if (sensor == MPC_missingValue_INT) then
            sensor = 0
            if (instrument == MPC_missingValue_INT) instrument = 0
          else
            instrument = obsu_cvt_obs_instrum(sensor)
          end if
          call obs_headSet_i(obsdat, OBS_SAT , headerIndex, obsSat)
          call obs_headSet_i(obsdat, OBS_INS , headerIndex, instrument)
        end if

      end if READHEADER

      ! Initialize some obsSpaceData body values
      if (obs_columnActive_RB(obsdat,OBS_LATD)) then
        call obs_bodySet_r(obsdat, OBS_LATD, bodyIndex, obs_missingValue_R)
      end if
      if (obs_columnActive_RB(obsdat,OBS_LOND)) then
        call obs_bodySet_r(obsdat, OBS_LOND, bodyIndex, obs_missingValue_R)
      end if
      if (obs_columnActive_RB(obsdat,OBS_SEM)) then
        call obs_bodySet_r(obsdat, OBS_SEM, bodyIndex, 0.0)
      end if

      ! Copy body row into obsspacedata by looping over all columns that were read
      beamRangeFound = .false. ! Needed to overwrite obs_ppp, regardless of column order
      do columnIndex = 1, size(bodySqlNames)
        select case(trim(bodySqlNames(columnIndex)))
          case('VCOORD')
            vertCoord = bodyValues(rowIndex,columnIndex)
            if (trim(familyType) /= 'RA' .and. trim(familyType) /= 'TO') then
              vertCoord = vertCoord * vertCoordFact + elevReal * elevFact
            end if
            call obs_bodySet_r(obsdat, OBS_PPP, bodyIndex, vertCoord)
          case('VCOORD_TYPE')
            ! This is set earlier based on familyType
            call obs_bodySet_i(obsdat, OBS_VCO, bodyIndex, vertCoordType)
          case('VARNO')
            obsVarno = int(bodyValues(rowIndex,columnIndex))
            call obs_bodySet_i(obsdat, OBS_VNM, bodyIndex, obsVarno)
          case('OBSVALUE')
            obsValue = bodyValues(rowIndex,columnIndex)
            if (trim(familyType) == 'TO' .and. obsValue == MPC_missingValue_R8) then
              ! Is this really needed???
              obsValue = real(MPC_missingValue_R8,pre_obsReal)
            end if
            call obs_bodySet_r(obsdat, OBS_VAR, bodyIndex, obsValue)
          case('OBS_ERROR')
            obsError = bodyValues(rowIndex,columnIndex)
            if (obsError > 0.0D0) then
              call obs_bodySet_r(obsdat, OBS_OER, bodyIndex, obsError)
            end if
          case('FLAG')
            obsFlag = int(bodyValues(rowIndex,columnIndex))
            if (obsFlag == mpc_missingValue_INT) obsFlag = 0
            call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, obsFlag)
          case('SURF_EMISS')
            if (obs_columnActive_RB(obsdat,OBS_SEM)) then
              call obs_bodySet_r(obsdat, OBS_SEM, bodyIndex, bodyValues(rowIndex,columnIndex) * zemFact)
            end if
          case('BIAS_CORR')
            if (obs_columnActive_RB(obsdat,OBS_BCOR)) then
              call obs_bodySet_r(obsdat, OBS_BCOR, bodyIndex, bodyValues(rowIndex,columnIndex))
            end if
          case('RANGE')
            beamRangeFound = .true.
            beamRange = bodyValues(rowIndex,columnIndex)
            call obs_bodySet_r(obsdat, OBS_LOCI, bodyIndex, beamRange)
            ! elevation and azimuths are converted to radians and pre_obsReal precision
            beamElevationReal = beamElevation * MPC_RADIANS_PER_DEGREE_R8
            beamAzimuthReal   = beamAzimuth   * MPC_RADIANS_PER_DEGREE_R8

            call rdv_getlatlonHRfromRange(xlat, xlon, beamElevationReal, beamAzimuthReal, & !in
                                          elevReal, beamRange,                            & !in
                                          beamLat, beamLon, beamHeight, beamDistance)       !out
            call obs_bodySet_r(obsdat, OBS_LATD, bodyIndex, beamLat)
            call obs_bodySet_r(obsdat, OBS_LOND, bodyIndex, beamLon)
        end select
      end do

      ! Replace the value of obs_ppp if a radar beam range value was found
      if (beamRangeFound) call obs_bodySet_r(obsdat, OBS_PPP, bodyIndex, beamHeight)

      ! Activate the 'rejected by selection process' bit if observed value is missing
      if (obsValue == MPC_missingValue_R8) then
        obsFlag  = obs_bodyElem_i(obsdat, obs_flg, bodyIndex)
        call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
      end if

      ! In some cases we need to add extra row(s) to the BODY table
      if (.not. filt_bufrCodeAssimilated(obsVarno) .and. &
          .not. ovt_bufrCodeSkipped(obsVarno)) then
          
        ! Add an extra row to the obsSpaceData body table
        ! to contain quantity later calculated by ovt_transformObsValue
        call obs_setBodyPrimaryKey(obsdat, bodyIndex+1, -1)
        call sqlr_addExtraDataRow(obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                                  ovt_getDestinationBufrCode(obsVarno), &
                                  vertCoordType, bodyIndex+1)
        bodyIndex = bodyIndex + 1
        obsNlv = obsNlv + 1
        if (ovt_isWindObs(obsVarno)) then
          ! Add an extra row for the other wind component
          call obs_setBodyPrimaryKey(obsdat, bodyIndex+1, -1)
          call sqlr_addExtraDataRow(obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                                    ovt_getDestinationBufrCode(obsVarno,extra_opt=.true.), &
                                    vertCoordType, bodyIndex+1)
          bodyIndex = bodyIndex + 1
          obsNlv = obsNlv + 1
        end if
      end if

      ! Logic to determine whether we are finished with this header entry
      finishedWithHeader = .false.
      if (rowIndex == numberIDsRows) then
        ! This is the last body entry, nothing comes after.
        finishedWithHeader = .true.
      else
        if (headPrimaryKey /= bodyHeadKeys(rowIndex+1)) then
          ! we just treated the last body entry associated with this id_obs (headPrimaryKey)
          ! make header for this id_obs
          finishedWithHeader = .true.
        end if
      end if
      if (finishedWithHeader) then
        ! Record the number of body rows associated with this header entry
        call obs_headSet_i(obsdat, OBS_NLV, headerIndex, obsNlv)
        obsNlv = 0 ! and reset counter for next batch of body entries
      end if

    end do BODY
    call fSQL_commit(db)

    numHeader = obs_numHeader(obsdat)
    numBody   = obs_numBody(obsdat)
    write(*,*) 'sqlr_readSqlite: FIN numheader  =', numHeader
    write(*,*) 'sqlr_readSqlite: FIN numbody    =', numBody

    call fSQL_close(db, stat) 
    if (fSQL_error(stat) /= FSQL_OK) then
      write(*,*) 'sqlr_readSqlite: problem closing sqlite db', trim(fileName)
      call sqlu_handleError(stat, 'fSQL_close')
    end if

  end subroutine sqlr_readSqlite

  !--------------------------------------------------------------------------
  ! sqlr_addExtraDataRow
  !--------------------------------------------------------------------------
  subroutine sqlr_addExtraDataRow(obsdat, vertCoord, obsVarno, vertCoordType, numberData)
    !
    ! :Purpose: Initialize data values for an extra row in data table.
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsdat
    real(pre_obsReal), intent(in)    :: vertCoord
    integer          , intent(in)    :: obsVarno
    integer          , intent(in)    :: vertCoordType
    integer          , intent(in)    :: numberData

    call obs_bodySet_r(obsdat, OBS_PPP,  numberData, vertCoord)
    call obs_bodySet_r(obsdat, OBS_VAR,  numberData, obs_missingValue_R)
    call obs_bodySet_i(obsdat, OBS_VNM,  numberData, obsVarno)
    call obs_bodySet_i(obsdat, OBS_FLG,  numberData, 0)
    call obs_bodySet_i(obsdat, OBS_VCO,  numberData, vertCoordType)
    call obs_bodySet_r(obsdat, OBS_LATD, numberData, obs_missingValue_R)
    call obs_bodySet_r(obsdat, OBS_LOND, numberData, obs_missingValue_R)
    if (obs_columnActive_RB(obsdat, OBS_OMA))  call obs_bodySet_r(obsdat, OBS_OMA , numberData, obs_missingValue_R)
    if (obs_columnActive_RB(obsdat, OBS_OMA0)) call obs_bodySet_r(obsdat, OBS_OMA0, numberData, obs_missingValue_R)
    if (obs_columnActive_RB(obsdat, OBS_OMP))  call obs_bodySet_r(obsdat, OBS_OMP , numberData, obs_missingValue_R)
    if (obs_columnActive_RB(obsdat, OBS_OMP6)) call obs_bodySet_r(obsdat, OBS_OMP6, numberData, obs_missingValue_R)
    if (obs_columnActive_RB(obsdat, OBS_OER))  call obs_bodySet_r(obsdat, OBS_OER , numberData, obs_missingValue_R)
    if (obs_columnActive_RB(obsdat, OBS_HPHT)) call obs_bodySet_r(obsdat, OBS_HPHT, numberData, obs_missingValue_R)
    if (obs_columnActive_RB(obsdat, OBS_HAHT)) call obs_bodySet_r(obsdat, OBS_HAHT, numberData, obs_missingValue_R)
    if (obs_columnActive_RB(obsdat, OBS_WORK)) call obs_bodySet_r(obsdat, OBS_WORK, numberData, obs_missingValue_R)
    if (obs_columnActive_RB(obsdat, OBS_SIGI)) call obs_bodySet_r(obsdat, OBS_SIGI, numberData, obs_missingValue_R)
    if (obs_columnActive_RB(obsdat, OBS_SIGO)) call obs_bodySet_r(obsdat, OBS_SIGO, numberData, obs_missingValue_R)
    if (obs_columnActive_RB(obsdat, OBS_ZHA))  call obs_bodySet_r(obsdat, OBS_ZHA , numberData, obs_missingValue_R)
    if (obs_columnActive_RB(obsdat, OBS_BCOR)) call obs_bodySet_r(obsdat, OBS_BCOR, numberData, obs_missingValue_R)

  end subroutine sqlr_addExtraDataRow

  !--------------------------------------------------------------------------
  ! sqlr_addColumn
  !--------------------------------------------------------------------------
  subroutine sqlr_addColumn(obsSpaceColIndexSource, columnName, tableName, fileName)
    !
    ! :Purpose: Add columns to sqlite tables that does not previously exists.
    !
    implicit none

    ! Arguments:  
    character(len=*)   , intent(in)    :: fileName    
    character(len=*)   , intent(in)    :: tableName   
    character(len=*)   , intent(in)    :: columnName   
    integer,             intent(in)    :: obsSpaceColIndexSource  

    ! Locals: 
    character(len=3000)         :: query
    character(len=10)           :: sqlDataType
    character(len=*), parameter :: myName = 'sqlr_addColumn'
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myName//': fSQL_open '//fSQL_errmsg(stat) )
    end if

    if ( obs_columnDataType(obsSpaceColIndexSource) == 'real' ) then
      sqlDataType = 'double'
    else
      sqlDataType = 'integer'
    end if

    query = 'alter table ' // trim(tableName) // ' add column ' // &
                trim(columnName) // ' ' // trim(sqlDataType) // ';'

    write(*,*) myName//': query = ', trim(query)
    call fSQL_do_many( db, query, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myName//': Problem with fSQL_do_many '//fSQL_errmsg(stat) )
    end if

    ! close the sqlite file
    call fSQL_close( db, stat ) 

  end subroutine sqlr_addColumn  

  !--------------------------------------------------------------------------
  ! sqlr_updateSqlite
  !--------------------------------------------------------------------------
  subroutine sqlr_updateSqlite(db, obsdat, familyType, fileName, fileNumber)
    !
    ! :Purpose: update SQLite files. List of items to update is in the 
    !           namSQLUpdate namelist
    implicit none
    
    ! Arguments:
    type(fSQL_database), intent(inout) :: db         ! SQL database
    type (struct_obs)  , intent(inout) :: obsdat     ! obsSpaceData
    character(len=*)   , intent(in)    :: fileName   ! file name  
    character(len=*)   , intent(in)    :: familyType ! Observation Family Type 
    integer            , intent(in)    :: fileNumber ! FILE NUMBER ASSOCIATED WITH db
    
    ! Locals:
    type(fSQL_statement)        :: stmt ! prepared statement for  SQLite
    type(fSQL_status)           :: stat ! type error status
    integer                     :: obsRln, obsNlv, obsIdf, obsFlag
    integer                     :: obsStatus, last_question, landSea, terrainType
    integer(8)                  :: headPrimaryKey, bodyPrimaryKey
    integer                     :: itemIndex, headPrimaryKeyIndex, landSeaIndex
    integer                     :: headerIndex, bodyIndex
    character(len =   3)        :: item
    integer                     :: updateList(20), fnom, fclos, nulnam, ierr
    character(len =  10)        :: columnName
    character(len = 128)        :: query
    character(len = 356)        :: itemChar, columnNameChar
    logical                     :: back, nonEmptyColumn, nonEmptyColumn_mpiglobal
    real(4)                     :: romp, obsValue, scaleFactor, columnValue
    integer, parameter :: maxNumberUpdate = 15

    ! Namelist variables:
    integer          :: numberUpdateItems                    ! MUST NOT BE INCLUDED IN NAMELIST!
    integer          :: numberUpdateItemsRadar               ! MUST NOT BE INCLUDED IN NAMELIST!
    character(len=3) :: itemUpdateList(maxNumberUpdate)      ! List of columns to be updated (e.g.'OMA','OMP')
    character(len=3) :: itemUpdateListRadar(maxNumberUpdate) ! List of columns to be updated for Radar data

    namelist/namSQLUpdate/ numberUpdateItems,      itemUpdateList,     &
                           numberUpdateItemsRadar, itemUpdateListRadar

    write(*,*) 'sqlr_updateSqlite: Starting ===================  '

    ! set default values of namelist variables
    itemUpdateList(:) = ''
    itemUpdateListRadar(:) = ''
    numberUpdateItems = MPC_missingValue_INT
    numberUpdateItemsRadar = MPC_missingValue_INT

    ! Read the namelist for directives
    nulnam = 0
    ierr   = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam,nml = namSQLUpdate, iostat = ierr)
    if (ierr /= 0) call utl_abort('sqlr_updateSqlite: Error reading namelist')
    if (mmpi_myid == 0) write(*, nml = namSQLUpdate)
    ierr = fclos(nulnam)
    if (numberUpdateItems /= MPC_missingValue_INT) then
      call utl_abort('sqlr_updateSqlite: check namelist section namSQLUpdate, numberUpdateItems should be removed')
    end if
    if (numberUpdateItemsRadar /= MPC_missingValue_INT) then
      call utl_abort('sqlr_updateSqlite: check namelist section namSQLUpdate, numberUpdateItemsRadar should be removed')
    end if
    numberUpdateItems = 0
    do itemIndex = 1, maxNumberUpdate
      if (trim(itemUpdateList(itemIndex)) == '') exit
      numberUpdateItems = numberUpdateItems + 1
    end do
    numberUpdateItemsRadar = 0
    do itemIndex = 1, maxNumberUpdate
      if (trim(itemUpdateListRadar(itemIndex)) == '') exit
      numberUpdateItemsRadar = numberUpdateItemsRadar + 1
    end do
    
    ! Append extra sqlite columns to update to itemUpdateList
    if (trim(familyType) == 'RA') then
      do itemIndex = 1, numberUpdateItemsRadar
        numberUpdateItems = numberUpdateItems + 1
        itemUpdateList(numberUpdateItems) = itemUpdateListRadar(itemIndex)
      end do
    end if

    write(*,*) 'sqlr_updateSqlite: family type   = ', trim(familyType)
    write(*,*) 'sqlr_updateSqlite: number of items to update: ', numberUpdateItems
    write(*,*) 'sqlr_updateSqlite: file name     = ', trim(fileName)
    write(*,*) 'sqlr_updateSqlite: missing value = ', MPC_missingValue_R8    

    ! create query
    itemChar='  '

    do itemIndex = 1, numberUpdateItems
      
      item = itemUpdateList(itemIndex)
      write(*,*) 'sqlr_updateSqlite: updating ', itemIndex, trim(item)

      select case(item)
      case('OMA')
        updateList(itemIndex) = OBS_OMA
        columnName = 'oma'
      case('OMP')
        updateList(itemIndex) = OBS_OMP
        columnName = 'omp'
      case('VAR')
        updateList(itemIndex) = OBS_VAR
        columnName = 'obsvalue'
      case('OER')
        updateList(itemIndex) = OBS_OER
        columnName = 'obs_error'
      case('FGE')
        updateList(itemIndex) = OBS_HPHT
        columnName = 'fg_error'
      case('EMI')
        updateList(itemIndex) = OBS_SEM
        columnName = 'surf_emiss'
      case('COR')
        updateList(itemIndex) = OBS_BCOR
        columnName = 'bias_corr'
      case('ALT')
        updateList(itemIndex) = OBS_PPP
        columnName = 'vcoord'
      case DEFAULT
        call utl_abort('sqlr_updateSqlite: invalid item '// columnName //' EXIT sqlr_updateSQL!!!')
      end select

      ! Check if column exist. If not, add column when corresponding
      ! obsspacedata variable have non-missing values
      if (sqlu_sqlColumnExists(fileName, 'data', columnName)) then
        itemChar = trim(itemChar)//','// trim(columnName) // ' = ? '
      else
        write(*,*) 'sqlr_updateSqlite: WARNING: column '//columnName// &
                   ' does not exist in the file '//trim(fileName)

        nonEmptyColumn = .false.
        ! Check if the ObsSpaceData variable contains non-missing values
        HEADERCHCK: do headerIndex = 1, obs_numHeader(obsdat)

          obsIdf = obs_headElem_i(obsdat,OBS_IDF, headerIndex)
          if (obsIdf /= fileNumber) cycle HEADERCHCK
          
          obsRln = obs_headElem_i(obsdat, OBS_RLN, headerIndex)
          obsNlv = obs_headElem_i(obsdat, OBS_NLV, headerIndex)
    
          BODYCHCK: do bodyIndex = obsRln, obsNlv + obsRln - 1
            columnValue = obs_bodyElem_r(obsdat, updateList(itemIndex), bodyIndex)

            if (columnValue /= obs_missingValue_R) then
              nonEmptyColumn = .true.         
              exit HEADERCHCK
            end if
          
          end do BODYCHCK
        end do HEADERCHCK

        call rpn_comm_allreduce(nonEmptyColumn,nonEmptyColumn_mpiglobal,1, &
                              "MPI_LOGICAL","MPI_LOR","grid",ierr)

        ! Add column into SQLite file if ObsSpaceData value containes non-missing values
        if (nonEmptyColumn_mpiglobal) then
          write(*,*) 'sqlr_updateSqlite: ' // trim(columnName) // ' column is not empty and will be added'
          call sqlr_addColumn(updateList( itemIndex ), columnName, 'data', fileName)
          itemChar = trim(itemChar)//','// trim(columnName) // ' = ? '
        end if

      end if
    end do

    back=.true.
    last_question  = scan(itemChar, '?', back)
    columnNameChar = itemChar(1:last_question)
    itemChar    = columnNameChar
    
    query = ' update data set flag = ? '//trim(itemChar)
    query = trim(query)//' where id_data = ?  ;'
    write(*,*) 'sqlr_updateSqlite: update query --->  ', query
    call fSQL_do_many(db, 'PRAGMA  synchronous = OFF; PRAGMA journal_mode = OFF;')
    call fSQL_prepare(db, query , stmt, stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'fSQL_prepare : ')
    call fSQL_begin(db)

    HEADER: do headerIndex = 1, obs_numHeader(obsdat)
 
      obsIdf = obs_headElem_i(obsdat,OBS_IDF, headerIndex)
 
      if (obsIdf /= fileNumber) cycle HEADER
      headPrimaryKey = obs_headPrimaryKey(obsdat, headerIndex)
      obsRln = obs_headElem_i(obsdat, OBS_RLN, headerIndex)
      obsNlv = obs_headElem_i(obsdat, OBS_NLV, headerIndex)

      BODY: do bodyIndex = obsRln, obsNlv + obsRln - 1

        obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        bodyPrimaryKey  = obs_bodyPrimaryKey(obsdat, bodyIndex)

        call fSQL_bind_param(stmt, param_index = 1, int_var = obsFlag)
        
        ITEMS: do itemIndex = 1, numberUpdateItems

          obsValue = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          if (obsValue /= obs_missingValue_R) then  
            romp = obs_bodyElem_r(obsdat, updateList(itemIndex), bodyIndex)
            if (romp == obs_missingValue_R) then
              call fSQL_bind_param(stmt, param_index = itemIndex + 1) ! sql null values
            else
              scaleFactor=1.0
              if (updateList(itemIndex) == OBS_SEM) scaleFactor=100.0
              call fSQL_bind_param(stmt, param_index = itemIndex + 1, real_var = romp*scaleFactor)
            end if
          else
            call fSQL_bind_param(stmt, param_index = itemIndex + 1)  ! sql null values
          end if

        end do ITEMS

        call fSQL_bind_param(stmt, param_index = numberUpdateItems + 2, int8_var  = bodyPrimaryKey)
        call fSQL_exec_stmt (stmt)

      end do BODY

    end do HEADER

    call fSQL_finalize(stmt)

    if (trim(familyType) /= 'GL'.and. trim(familyType) /= 'RA') then

      ! UPDATES FOR THE STATUS FLAGS and land_sea (for satellites) IN THE HEADER TABLE
      ! The variables 'headPrimaryKeyIndex' and 'landSeaIndex' are defined here and used below.
      ! Any change in this logic must be coherent with the code below!
      if (trim(familyType) == 'TO') then
        query = ' update header set status  = ?, land_sea= ? where id_obs = ? '
        landSeaIndex = 2
        headPrimaryKeyIndex = 3
      else
        query = ' update header set status  = ?  where id_obs = ? '
        headPrimaryKeyIndex = 2
      end if
      call fSQL_prepare(db, query , stmt, stat)
      if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat,'fSQL_prepare : ')

      HEADER2: do headerIndex = 1,obs_numHeader(obsdat)
         
        terrainType=MPC_missingValue_INT
        obsIdf = obs_headElem_i(obsdat, OBS_IDF, headerIndex)
        if (obsIdf /= fileNumber) cycle HEADER2
        headPrimaryKey = obs_headPrimaryKey(obsdat, headerIndex)
        obsStatus = obs_headElem_i(obsdat, OBS_ST1, headerIndex)
        landsea   = obs_headElem_i(obsdat, OBS_STYP,headerIndex)
        call fSQL_bind_param(stmt, param_index = 1, int_var  = obsStatus)
        ! The variables 'headPrimaryKeyIndex' and 'landSeaIndex' are defined above and
        ! they must be coherent with the query designed above
        call fSQL_bind_param(stmt, param_index = headPrimaryKeyIndex, int8_var = headPrimaryKey)
        if (trim(familyType) == 'TO') then
          call fSQL_bind_param(stmt, param_index = landSeaIndex, int_var  = landSea)
        end if

        call fSQL_exec_stmt (stmt)

      end do HEADER2
    
      call fSQL_finalize(stmt)

    end if

    call fSQL_commit(db, stat)
    if (fSQL_error(stat)  /= FSQL_OK) then
      call sqlu_handleError(stat, 'sqlr_updateSqlite: fSQL_commit')
    end if
    write(*,*) 'sqlr_updateSqlite: End ===================  ', trim(familyType)

  end subroutine sqlr_updateSqlite

  !--------------------------------------------------------------------------
  ! sqlr_addCloudParametersandEmissivity
  !--------------------------------------------------------------------------
  subroutine sqlr_addCloudParametersandEmissivity(db, obsdat, fileNumber)
    !
    !:Purpose: Add a new table if it doesn't already exist `cld_params` with
    !          cloud-related information.
    !
    implicit none

    ! Arguments:
    type(fSQL_DATABASE), intent(inout) :: db   ! SQLite file handle
    type(struct_obs),    intent(in)    :: obsdat
    integer,             intent(in)    :: fileNumber

    ! Locals:
    type(fSQL_STATEMENT)   :: stmt ! type for precompiled SQLite statements
    type(fSQL_STATUS)      :: stat !type for error status
    character(len = 256)   :: query
    integer                :: numberInsert, headerIndex, obsIdo, obsIdf
    integer                :: NCO2
    real                   :: ETOP,VTOP,ECF,VCF,HE,ZTSR,ZTM,ZTGM,ZLQM,ZPS

    ! First create the table if it does not already exist
    query = 'create table if not exists cld_params(id_obs integer,ETOP real,VTOP real, &
             ECF real,VCF real,HE real,ZTSR real,NCO2 integer,ZTM real,ZTGM real,ZLQM real,ZPS real);'
    write(*,*) 'sqlr_addCloudParametersandEmissivity: Create query = ', trim(query)

    call fSQL_do(db, trim(query), stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'fSQL_do : ')

    ! Insert values in the table
    query = 'insert into cld_params(id_obs,ETOP,VTOP,ECF,VCF,HE,ZTSR,NCO2,ZTM,ZTGM,ZLQM,ZPS) &
             values(?,?,?,?,?,?,?,?,?,?,?,?);'
    write(*,*) 'sqlr_addCloudParametersandEmissivity: Insert query = ', trim(query)

    call fSQL_prepare(db, query, stmt, stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'fSQL_prepare : ')
    call fSQL_begin(db)
    numberInsert=0
    HEADER: do headerIndex = 1, obs_numHeader(obsdat)
      obsIdf = obs_headElem_i(obsdat, OBS_IDF, headerIndex)
      if (obsIdf /= fileNumber) cycle HEADER
      obsIdo   = obs_headPrimaryKey(obsdat, headerIndex)
      ETOP = obs_headElem_r(obsdat, OBS_ETOP, headerIndex)
      VTOP = obs_headElem_r(obsdat, OBS_VTOP, headerIndex)
      ECF  = obs_headElem_r(obsdat, OBS_ECF,  headerIndex)
      VCF  = obs_headElem_r(obsdat, OBS_VCF,  headerIndex)
      HE   = obs_headElem_r(obsdat, OBS_HE,   headerIndex)
      ZTSR = obs_headElem_r(obsdat, OBS_ZTSR, headerIndex)
      NCO2 = obs_headElem_i(obsdat, OBS_NCO2, headerIndex)
      ZTM  = obs_headElem_r(obsdat, OBS_ZTM,  headerIndex)
      ZTGM = obs_headElem_r(obsdat, OBS_ZTGM, headerIndex)
      ZLQM = obs_headElem_r(obsdat, OBS_ZLQM, headerIndex)
      ZPS  = obs_headElem_r(obsdat, OBS_ZPS,  headerIndex)

      call fSQL_bind_param(stmt, param_index = 1, int_var  = obsIdo)
      call fSQL_bind_param(stmt, param_index = 2, real_var = ETOP)
      call fSQL_bind_param(stmt, param_index = 3, real_var = VTOP)
      call fSQL_bind_param(stmt, param_index = 4, real_var = ECF)
      call fSQL_bind_param(stmt, param_index = 5, real_var = VCF)
      call fSQL_bind_param(stmt, param_index = 6, real_var = HE)
      call fSQL_bind_param(stmt, param_index = 7, real_var = ZTSR)
      call fSQL_bind_param(stmt, param_index = 8, int_var  = NCO2)
      call fSQL_bind_param(stmt, param_index = 9, real_var = ZTM)
      call fSQL_bind_param(stmt, param_index = 10,real_var = ZTGM)
      call fSQL_bind_param(stmt, param_index = 11,real_var = ZLQM)
      call fSQL_bind_param(stmt, param_index = 12,real_var = ZPS)

      call fSQL_exec_stmt (stmt)
      numberInsert=numberInsert +1
    end do HEADER

    call fSQL_finalize(stmt)
    call fSQL_commit(db)
    write(*,'(a,i8)') 'sqlr_addCloudParametersandEmissivity: NUMBER OF INSERTIONS ----> ', numberInsert

  end subroutine sqlr_addCloudParametersandEmissivity

  !--------------------------------------------------------------------------
  ! sqlr_insertSqlite
  !--------------------------------------------------------------------------
  subroutine sqlr_insertSqlite(db, obsdat, familyType, fileName, fileNumber)
    !
    !:Purpose: Insert rows in the sqlite file data table for bufr element IDs
    !          specified in the namelist block `namSQLInsert`.
    !
    implicit none

    ! Arguments:
    type(fSQL_DATABASE), intent(inout) :: db   ! type for SQLIte  file handle
    type(struct_obs),    intent(in)    :: obsdat
    character(len=*),    intent(in)    :: familyType
    character(len=*),    intent(in)    :: fileName
    integer         ,    intent(in)    :: fileNumber

    ! Locals:
    type(fSQL_STATEMENT)   :: stmt ! type for precompiled SQLite statements
    type(fSQL_STATUS)      :: stat !type for error status
    integer                :: obsVarno, obsFlag, vertCoordType, fnom, fclos, nulnam, ierr 
    real                   :: obsValue, OMA, OMP, OER, FGE, PPP
    integer                :: numberInsert, headerIndex, bodyIndex, numHeader, itemIndex
    integer                :: obsNlv, obsRln, obsIdf, insertItem
    integer(8)             :: bodyPrimaryKey, headPrimaryKey
    character(len = 256)   :: query
    logical                :: llok    
    integer, parameter     :: maxNumberInsertItems = 15

    ! Namelist variables
    integer :: itemInsertList(maxNumberInsertItems) ! List of bufr element ids to insert in sql file data table
    integer :: numberInsertItems                    ! MUST NOT BE INCLUDED IN NAMELIST!

    namelist/namSQLInsert/ numberInsertItems, itemInsertList

    write(*,*)  'sqlr_insertSqlite: --- Starting ---   '
    numHeader = obs_numHeader(obsdat)
    write(*,*) ' FAMILY ---> ', trim(familyType), '  headerIndex  ----> ', numHeader
    write(*,*) ' fileName -> ', trim(fileName)

    ! set default values of namelist variables
    numberInsertItems = MPC_missingValue_INT
    itemInsertList(:) = MPC_missingValue_INT

    nulnam = 0
    ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml = namSQLInsert, iostat = ierr)
    if (ierr /= 0) call utl_abort('sqlr_insertSqlite: Error reading namelist')
    if (mmpi_myid == 0) write(*, nml = namSQLInsert)
    ierr=fclos(nulnam)
    if (numberInsertItems /= MPC_missingValue_INT) then
      call utl_abort('sqlr_insertSqlite: check namSQLInsert namelist section, you need to remove numberInsertItems')
    end if
    numberInsertItems = 0
    do itemIndex = 1, maxNumberInsertItems
      if ( itemInsertList(itemIndex) == MPC_missingValue_INT) exit
      numberInsertItems = numberInsertItems + 1
    end do
    
    write(*,*) ' INSERT INTO SQLITE FILE ELEMENTS :--> ',(itemInsertList(insertItem), insertItem = 1, numberInsertItems)

    select case(trim(familyType))

      case('SF', 'SC', 'GP')

        query = 'insert into data (id_obs,varno,vcoord,obsvalue,flag,oma,omp,fg_error,obs_error) values(?,?,?,?,?,?,?,?,?);'

      case DEFAULT

        query = 'insert into data (id_obs,varno,vcoord,vcoord_type,obsvalue,flag,oma,omp,fg_error,obs_error) values(?,?,?,?,?,?,?,?,?,?);'

    end select

    query=trim(query)
    write(*,*) ' === Family Type === ',trim(familyType)
    write(*,*) ' Insert query = ', trim(query)

    call fSQL_begin(db)
    call fSQL_prepare(db, query, stmt, stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'fSQL_prepare : ')

    numberInsert=0
    HEADER: do headerIndex = 1, obs_numHeader(obsdat)

      obsIdf = obs_headElem_i(obsdat, OBS_IDF, headerIndex)
      if (obsIdf /= fileNumber) cycle HEADER
      headPrimaryKey = obs_headPrimaryKey(obsdat, headerIndex)
      obsRln = obs_headElem_i(obsdat, OBS_RLN, headerIndex)
      obsNlv = obs_headElem_i(obsdat, OBS_NLV, headerIndex)

      BODY: do bodyIndex = obsRln, obsNlv + obsRln -1

        bodyPrimaryKey= obs_bodyPrimaryKey(obsdat, bodyIndex)
        obsVarno      = obs_bodyElem_i(obsdat, OBS_VNM , bodyIndex)
        obsFlag       = obs_bodyElem_i(obsdat, OBS_FLG , bodyIndex)
        vertCoordType = obs_bodyElem_i(obsdat, OBS_VCO , bodyIndex)
        obsValue      = obs_bodyElem_r(obsdat, OBS_VAR , bodyIndex)
        OMA           = obs_bodyElem_r(obsdat, OBS_OMA , bodyIndex)
        OMP           = obs_bodyElem_r(obsdat, OBS_OMP , bodyIndex)
        OER           = obs_bodyElem_r(obsdat, OBS_OER , bodyIndex)
        FGE           = obs_bodyElem_r(obsdat, OBS_HPHT, bodyIndex)
        PPP           = obs_bodyElem_r(obsdat, OBS_PPP , bodyIndex)

        llok = .false.
        do insertItem = 1, numberInsertItems
          if (obsVarno == itemInsertList(insertItem)) llok=.true.
        end do

        if (bodyPrimaryKey == -1 .and. llok) then
          select case(trim(familyType))
            case('SF', 'SC', 'GP')
              call fSQL_bind_param(stmt, param_index = 1, int8_var = headPrimaryKey)
              call fSQL_bind_param(stmt, param_index = 2, int_var  = obsVarno)
              call fSQL_bind_param(stmt, param_index = 3, real_var = PPP)
              if (obsValue == obs_missingValue_R) then          ! sql null values
                call fSQL_bind_param(stmt, param_index = 4)
              else
                call fSQL_bind_param(stmt, param_index = 4, real_var = obsValue)
              end if
              call fSQL_bind_param(stmt, param_index = 5, int_var  = obsFlag)
              if (OMA == obs_missingValue_R) then
                call fSQL_bind_param(stmt, param_index = 6) 
              else 
                call fSQL_bind_param(stmt, param_index = 6, real_var = OMA) 
              end if
              if (OMP == obs_missingValue_R) then
                call fSQL_bind_param(stmt, param_index = 7) 
              else
                call fSQL_bind_param(stmt, param_index = 7, real_var = OMP) 
              end if
              if (FGE == obs_missingValue_R) then
                call fSQL_bind_param(stmt, param_index = 8) 
              else
                call fSQL_bind_param(stmt, param_index = 8, real_var = FGE)
              end if
              if (OER == obs_missingValue_R) then
                call fSQL_bind_param(stmt, param_index = 9) 
              else
                call fSQL_bind_param(stmt, param_index = 9, real_var = OER) 
              end if
            case DEFAULT
              call fSQL_bind_param(stmt, param_index = 1, int8_var = headPrimaryKey)
              call fSQL_bind_param(stmt, param_index = 2, int_var  = obsVarno)
              call fSQL_bind_param(stmt, param_index = 3, real_var = PPP) 
              call fSQL_bind_param(stmt, param_index = 4, int_var  = vertCoordType)
              if (obsValue == obs_missingValue_R) then
                call fSQL_bind_param(stmt, param_index = 5)
              else
                call fSQL_bind_param(stmt, param_index = 5, real_var = obsValue)
              end if
              call fSQL_bind_param(stmt, param_index = 6, int_var  = obsFlag)
              if (OMA == obs_missingValue_R) then
                call fSQL_bind_param(stmt, param_index = 7) 
              else 
                call fSQL_bind_param(stmt, param_index = 7, real_var = OMA)
              end if
              if (OMP == obs_missingValue_R) then
                call fSQL_bind_param(stmt, param_index = 8) 
              else
                call fSQL_bind_param(stmt, param_index = 8, real_var = OMP) 
              end if
              if (FGE == obs_missingValue_R) then
                call fSQL_bind_param(stmt, param_index = 9) 
              else
                call fSQL_bind_param(stmt, param_index = 9, real_var = FGE) 
              end if
              if (OER == obs_missingValue_R) then
                call fSQL_bind_param(stmt, param_index = 10) 
              else
                call fSQL_bind_param(stmt, param_index = 10, real_var = OER)
              end if
          end select
          call fSQL_exec_stmt (stmt)

          numberInsert = numberInsert + 1
        end if

      end do BODY

    end do HEADER

    call fSQL_finalize(stmt)
    call fSQL_commit(db)
    write(*,'(3a,i8)') 'sqlr_insertSqlite: FAMILY ---> ' ,trim(familyType), &
         '  NUMBER OF INSERTIONS ----> ', numberInsert

  end subroutine sqlr_insertSqlite

  !--------------------------------------------------------------------------
  ! sqlr_cleanSqlite
  !--------------------------------------------------------------------------
  subroutine sqlr_cleanSqlite(db, fileName)
    !
    ! :Purpose: Remove flagged (bit 11 set) observations in an SQLite file
    !
    implicit none

    ! Arguments:
    type(fSQL_DATABASE), intent(inout) :: db   ! SQLite file handle
    character(len=*),    intent(in)    :: fileName

    ! Locals:
    character(len = 128) :: query
    type(fSQL_STATEMENT) :: statement ! prepared statement for SQLite
    type(fSQL_STATUS)    :: status

    call fSQL_open(db, fileName, status)
    if (fSQL_error(status) /= FSQL_OK) then
      write(*,*) 'sqlr_cleanSqlite: ERROR: ', fSQL_errmsg(status)
    end if
    ! Mark for deletion all records with bit 11 (2048) set
    query = ' delete from data where flag & 2048;'
    call fSQL_prepare(db, query, statement, status)
    if (fSQL_error(status) /= FSQL_OK) &
      call sqlu_handleError(status, 'thinning fSQL_prepare : ')
    call fSQL_begin(db)
    call fSQL_exec_stmt(statement)
    call fSQL_finalize(statement)
    call fSQL_commit(db)
    write(*,*) 'sqlr_cleanSqlite: closed database -->', trim(FileName)
    call fSQL_close(db, status)

  end subroutine sqlr_cleanSqlite

  !--------------------------------------------------------------------------
  ! sqlr_writePseudoSSTobs
  !--------------------------------------------------------------------------
  subroutine sqlr_writePseudoSSTobs(obsData, obsFamily, instrumentFileName)
    !
    ! :Purpose: To write the obsSpaceData content into SQLite format files
    !
    implicit none

    ! Arguments:
    type(struct_obs) , intent(inout) :: obsData
    character(len=*) , intent(in)    :: obsFamily
    character(len=*) , intent(in)    :: instrumentFileName
            
    ! Locals:
    type(fSQL_DATABASE)    :: db                        ! type for SQLIte  file handle
    type(fSQL_STATEMENT)   :: stmtData, stmtHeader      ! type for precompiled SQLite statements
    type(fSQL_STATUS)      :: stat                      ! type for error status
    integer                :: obsVarno, obsFlag, ASS, codeType, date, time, idObs, idData
    integer, parameter     :: obsStatus = 3072 
    real                   :: obsValue, OMA, OMP, OER, FGE, PPP, lon, lat, altitude
    integer                :: numberInsertions, numHeaders, headerIndex, bodyIndex, obsNlv, obsRln
    character(len = 512)   :: queryData, queryHeader, queryCreate
    character(len = 12)    :: idStation
    character(len=256)     :: fileName, fileNameDir
    character(len=4)       :: cmyidx, cmyidy
        
    write(*,*) 'sqlr_writePseudoSSTobs: starting...'
     
    ! determine initial idData,idObs to ensure unique values across mpi tasks
    call sqlu_getInitialIdObsData(obsData, obsFamily, idObs, idData)

    ! return if this mpi task does not have any observations of this family
    if (trim(instrumentFileName) == 'XXXXX') return

    ! check if any obs exist for this file, if not return
    numHeaders = 0
    call obs_set_current_header_list(obsData, obsFamily)
    HEADERCOUNT: do
      headerIndex = obs_getHeaderIndex(obsData)
      if (headerIndex < 0) exit HEADERCOUNT
      numHeaders = numHeaders + 1
    end do HEADERCOUNT
    if (numHeaders == 0) return

    fileNameDir = trim(ram_getRamDiskDir())
    if (fileNameDir == ' ') then
      write(*,*) 'sqlr_writePseudoSSTobs: WARNING! The program may be slow creating many sqlite files in the same directory.'
      write(*,*) 'sqlr_writePseudoSSTobs: WARNING! Please, use the ram disk option prior to MIDAS run!'
    end if  

    if (obs_mpiLocal(obsData)) then
      write(cmyidy,'(I4.4)') (mmpi_myidy + 1)
      write(cmyidx,'(I4.4)') (mmpi_myidx + 1)
      fileName = trim(fileNameDir)//'obs/'//trim(instrumentFileName)//'_'//trim(cmyidx)//'_'//trim(cmyidy)
    else
      if (mmpi_myid > 0) return
      fileName = trim(fileNameDir)//'obs/'//trim(instrumentFileName)
    end if
    

    write(*,*) 'sqlr_writePseudoSSTobs: Creating file: ', trim(fileName)
    call fSQL_open(db, fileName, stat)
    if (fSQL_error(stat) /= FSQL_OK) write(*,*) 'sqlr_writePseudoSSTobs: fSQL_open: ', fSQL_errmsg(stat),' filename: '//trim(fileName)

    ! Create the tables HEADER and DATA
    queryCreate = 'create table header (id_obs integer primary key, id_stn varchar(50), lat real, lon real, &
                  &codtyp integer, date integer, time integer, elev real, status integer); &
                  &create table data (id_data integer primary key, id_obs integer, varno integer, vcoord real, &
                  &vcoord_type integer, obsvalue real, flag integer, oma real, ompt real, oma0 real, omp real, &
                  &an_error real, fg_error real, obs_error real);'
    
    call fSQL_do_many(db, queryCreate, stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'sqlr_writePseudoSSTobs: fSQL_do_many with query: '//trim(queryCreate))
    
    queryHeader = ' insert into header (id_obs, id_stn, lat, lon, date, time, codtyp, elev, status) values(?,?,?,?,?,?,?,?,?); '
    queryData = 'insert into data (id_data, id_obs, varno, vcoord, vcoord_type, obsvalue, flag, oma, oma0, ompt, fg_error, &
                &obs_error) values(?,?,?,?,?,?,?,?,?,?,?,?);'

    write(*,*) 'sqlr_writePseudoSSTobs: Insert query Data   = ', trim(queryData)
    write(*,*) 'sqlr_writePseudoSSTobs: Insert query Header = ', trim(queryHeader)

    call fSQL_begin(db)
    call fSQL_prepare(db, queryData, stmtData, stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'sqlr_writePseudoSSTobs: fSQL_prepare:')
    call fSQL_prepare(db, queryHeader, stmtHeader, stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'sqlr_writePseudoSSTobs: fSQL_prepare:')
     
    numberInsertions = 0

    call obs_set_current_header_list(obsData, obsFamily)
    HEADER: do

      headerIndex = obs_getHeaderIndex(obsData)
      if (headerIndex < 0) exit HEADER

      codeType  = obs_headElem_i(obsData, OBS_ITY, headerIndex)
      obsRln    = obs_headElem_i(obsData, OBS_RLN, headerIndex)
      obsNlv    = obs_headElem_i(obsData, OBS_NLV, headerIndex)
      idStation = obs_elem_c    (obsData, 'STID' , headerIndex)
      altitude  = obs_headElem_r(obsData, OBS_ALT, headerIndex)
      lon       = obs_headElem_r(obsData, OBS_LON, headerIndex) * MPC_DEGREES_PER_RADIAN_R8
      lat       = obs_headElem_r(obsData, OBS_LAT, headerIndex) * MPC_DEGREES_PER_RADIAN_R8
      if (lon > 180.) lon = lon - 360.
      date      = obs_headElem_i(obsData, OBS_DAT, headerIndex)
      time      = obs_headElem_i(obsData, OBS_ETM, headerIndex) * 100.

      idObs = idObs + 1
      call fSQL_bind_param(stmtHeader, param_index = 1, int_var  = idObs)
      call fSQL_bind_param(stmtHeader, param_index = 2, char_var = idStation)
      call fSQL_bind_param(stmtHeader, param_index = 3, real_var = lat) 
      call fSQL_bind_param(stmtHeader, param_index = 4, real_var = lon) 
      call fSQL_bind_param(stmtHeader, param_index = 5, int_var  = date) 
      call fSQL_bind_param(stmtHeader, param_index = 6, int_var  = time) 
      call fSQL_bind_param(stmtHeader, param_index = 7, int_var  = codeType) 
      call fSQL_bind_param(stmtHeader, param_index = 8, real_var = altitude)
      call fSQL_bind_param(stmtHeader, param_index = 9, int_var = obsStatus)
      call fSQL_exec_stmt (stmtHeader)

      BODY: do bodyIndex = obsRln, obsNlv + obsRln -1
         
        obsVarno      = obs_bodyElem_i(obsData, OBS_VNM , bodyIndex)
        obsFlag       = obs_bodyElem_i(obsData, OBS_FLG , bodyIndex)
        obsValue      = obs_bodyElem_r(obsData, OBS_VAR , bodyIndex)
        OMA           = obs_bodyElem_r(obsData, OBS_OMA , bodyIndex)
        OMP           = obs_bodyElem_r(obsData, OBS_OMP , bodyIndex)
        OER           = obs_bodyElem_r(obsData, OBS_OER , bodyIndex)
        FGE           = obs_bodyElem_r(obsData, OBS_HPHT, bodyIndex)
        PPP           = obs_bodyElem_r(obsData, OBS_PPP , bodyIndex)
        ASS           = obs_bodyElem_i(obsData, OBS_ASS , bodyIndex)

        ! insert order: id_obs,varno,vcoord,vcoord_type,obsvalue,flag,oma,oma0,ompt,fg_error,obs_error
        idData = idData + 1
        call fSQL_bind_param(stmtData, param_index = 1, int_var  = idData)
        call fSQL_bind_param(stmtData, param_index = 2, int_var  = idObs)
        call fSQL_bind_param(stmtData, param_index = 3, int_var  = obsVarno)
        call fSQL_bind_param(stmtData, param_index = 4, real_var = PPP)
        call fSQL_bind_param(stmtData, param_index = 5) 
        call fSQL_bind_param(stmtData, param_index = 6, real_var = obsValue) 
        call fSQL_bind_param(stmtData, param_index = 7, int_var  = obsFlag)
        if (OMA == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 8) 
          call fSQL_bind_param(stmtData, param_index = 9) 
        else
          call fSQL_bind_param(stmtData, param_index = 8, real_var = OMA)
          call fSQL_bind_param(stmtData, param_index = 9, real_var = OMA)
        end if
        if (OMP == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 10) 
        else
          call fSQL_bind_param(stmtData, param_index = 10, real_var = OMP)
        end if
        if (FGE == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 11) 
        else
          call fSQL_bind_param(stmtData, param_index = 11, real_var = FGE)
        end if
        if (OER == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 12) 
        else
          call fSQL_bind_param(stmtData, param_index = 12, real_var = OER)
        end if 
        call fSQL_exec_stmt (stmtData)

        numberInsertions = numberInsertions + 1

      end do BODY
     
    end do HEADER
    
    call fSQL_finalize(stmtData)

    write(*,*) 'sqlr_writePseudoSSTobs: Observation Family: ', obsFamily, &
               ', number of insertions: ', numberInsertions

    call fSQL_commit(db)
    call fSQL_close(db, stat)

  end subroutine sqlr_writePseudoSSTobs
  
  !--------------------------------------------------------------------------
  ! sqlr_writeEmptyPseudoSSTobsFile
  !--------------------------------------------------------------------------
  subroutine sqlr_writeEmptyPseudoSSTobsFile(obsData, obsFamily, instrumentFileName)
    !
    ! :Purpose: to generate an empty SQLite SST pseudo obs file for mpi tasks,
    !           with no sea-ice on them.
    !
    implicit none

    ! Arguments:
    type(struct_obs) , intent(inout) :: obsData   
    character(len=*) , intent(in)    :: obsFamily
    character(len=*) , intent(in)    :: instrumentFileName
            
    ! Locals:
    type(fSQL_DATABASE)    :: db                        ! type for SQLIte  file handle
    type(fSQL_STATEMENT)   :: stmtData, stmtHeader      ! type for precompiled SQLite statements
    type(fSQL_STATUS)      :: stat                      ! type for error status
    character(len = 512)   :: queryCreate
    character(len = 512)   :: queryHeader, queryData  
    integer                :: idObs, idData
    character(len=30)      :: fileNameExtention
    character(len=256)     :: fileName, fileNameDir
    character(len=4)       :: cmyidx, cmyidy
        
    ! determine initial idData,idObs to ensure unique values across mpi tasks
    call sqlu_getInitialIdObsData(obsData, obsFamily, idObs, idData)
    
    fileNameDir = trim(ram_getRamDiskDir())
    if (fileNameDir == ' ') &
    write(*,*) 'sqlr_writeEmptyPseudoSSTobsFile: WARNING! The program may be slow creating many sqlite files in the same directory.'
    write(*,*) 'sqlr_writeEmptyPseudoSSTobsFile: WARNING! Please, use the ram disk option prior to MIDAS run!'

    if (obs_mpiLocal(obsData)) then
      write(cmyidy,'(I4.4)') (mmpi_myidy + 1)
      write(cmyidx,'(I4.4)') (mmpi_myidx + 1)
      fileNameExtention  = trim(cmyidx) // '_' // trim(cmyidy)
    else
      if (mmpi_myid > 0) return
      fileNameExtention = ' '
    end if
    
    fileName = trim(fileNameDir) // 'obs/' // trim(instrumentFileName) // '_' // trim(fileNameExtention)

    write(*,*) 'sqlr_writeEmptyPseudoSSTobsFile: Creating file: ', trim(fileName)
    call fSQL_open(db, fileName, stat)
    if (fSQL_error(stat) /= FSQL_OK) write(*,*) 'sqlr_writeEmptyPseudoSSTobsFile: fSQL_open: ', fSQL_errmsg(stat),' filename: '//trim(fileName)

    ! Create the tables HEADER and DATA
    queryCreate = 'create table header (id_obs integer primary key, id_stn varchar(50), lat real, lon real, &
                  &codtyp integer, date integer, time integer, elev real, status integer); &
                  &create table data (id_data integer primary key, id_obs integer, varno integer, vcoord real, &
                  &vcoord_type integer, obsvalue real, flag integer, oma real, ompt real, oma0 real, omp real, &
                  &an_error real, fg_error real, obs_error real);'
    
    call fSQL_do_many(db, queryCreate, stat)
    if (fSQL_error(stat) /= FSQL_OK) then
      call sqlu_handleError(stat, 'sqlr_writeEmptyPseudoSSTobsFile: fSQL_do_many with query: '//trim(queryCreate))
    end if
    
    queryHeader = ' insert into header (id_obs, id_stn, lat, lon, date, time, codtyp, elev, status) values(?,?,?,?,?,?,?,?,?); '
    queryData = 'insert into data (id_data, id_obs, varno, vcoord, vcoord_type, obsvalue, flag, oma, oma0, ompt, fg_error, &
                &obs_error) values(?,?,?,?,?,?,?,?,?,?,?,?);'

    write(*,*) 'sqlr_writeEmptyPseudoSSTobsFile: Insert query Data   = ', trim(queryData)
    write(*,*) 'sqlr_writeEmptyPseudoSSTobsFile: Insert query Header = ', trim(queryHeader)

    call fSQL_begin(db)
    call fSQL_prepare(db, queryData, stmtData, stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'sqlr_writeEmptyPseudoSSTobsFile: fSQL_prepare:')
    call fSQL_prepare(db, queryHeader, stmtHeader, stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'sqlr_writeEmptyPseudoSSTobsFile: fSQL_prepare:')

    call fSQL_commit(db)
    call fSQL_close(db, stat)

  end subroutine sqlr_writeEmptyPseudoSSTobsFile

  !--------------------------------------------------------------------------
  ! sqlr_getColumnValuesDate
  !--------------------------------------------------------------------------
  subroutine sqlr_getColumnValuesDate(columnDateValues, columnTimeValues, fileName)
    !
    ! :Purpose: Read the date and time column values from sqlite file.
    !
    implicit none

    ! Arguments:
    integer, allocatable, intent(out) :: columnDateValues(:)
    integer, allocatable, intent(out) :: columnTimeValues(:)
    character(len=*),     intent(in)  :: fileName

    ! Locals:
    integer              :: numRows, numColumns, rowIndex
    character(len=20), allocatable :: columnValuesStr(:,:)
    character(len=3000)  :: query
    type(fSQL_STATUS)    :: stat ! sqlite error status
    type(fSQL_DATABASE)  :: db   ! sqlite file handle
    type(fSQL_STATEMENT) :: stmt ! precompiled sqlite statements

    ! open the sqlite file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'sqlr_getColumnValuesDate: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'sqlr_getColumnValuesDate: fSQL_open' )
    end if

    ! Get the date and time

    ! build the sqlite query
    query = 'select date, time from header;'
    write(*,*) 'sqlr_getColumnValuesDate: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query), stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'sqlr_getColumnValuesDate: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('sqlr_getColumnValuesDate: problem with fSQL_get_many')
    end if
    write(*,*) 'sqlr_getColumnValuesDate: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValuesStr(numRows,2) )
    call fSQL_fill_matrix( stmt, columnValuesStr )
    allocate( columnDateValues(numRows) )
    allocate( columnTimeValues(numRows) )
    do rowIndex = 1, numRows
      read(columnValuesStr(rowIndex,1),*) columnDateValues(rowIndex)
      read(columnValuesStr(rowIndex,2),*) columnTimeValues(rowIndex)
      columnTimeValues(rowIndex) = columnTimeValues(rowIndex)/100
    end do

    deallocate(columnValuesStr)

    ! close the sqlite file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine sqlr_getColumnValuesDate

end module sqliteRead_mod

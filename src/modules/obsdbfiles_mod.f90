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

module obsdbFiles_mod
  ! MODULE obsdbFiles (prefix='odbf' category='6. Observation input/output')
  !
  ! :Purpose: To read and update sqlite files that are in the new 'obsDB' format.
  !
  use mpi_mod
  use codePrecision_mod
  use mathPhysConstants_mod
  use bufr_mod
  use codtyp_mod
  use obsSpaceData_mod
  use fSQLite
  use utilities_mod
  use clib_interfaces_mod
  use obsUtil_mod
  use obsVariableTransforms_mod

  implicit none
  save
  private

  ! Public subroutines and functions:
  public :: odbf_setup, odbf_isActive, odbf_readFile, odbf_updateFile

 
  ! Arrays used to match obsDB column names with obsSpaceData column names

  integer, parameter :: lenSqlName    = 60
  integer, parameter :: sqlColIndex   = 1
  integer, parameter :: sqlTabIndex   = 1
  integer, parameter :: obsColIndex   = 2
  integer, parameter :: varNoColIndex = 2

  character(len=lenSqlName) :: headTableName = 'Rapport'
  character(len=lenSqlName) :: bodyTableName = 'Observation'

  ! ...for the header table
  integer, parameter :: numHeadMatch = 14
  character(len=lenSqlName) :: headKeySqlName  = 'ID_RAPPORT'
  character(len=lenSqlName) :: headDateSqlName = 'DATE_VALIDITE'
  character(len=lenSqlName) :: headMatchList(2,numHeadMatch) = (/ &
       'ID_STN',              'STID', &
       'TYPE',                'ITY',  &
       'LAT',                 'LAT ', &
       'LON',                 'LON ', &
       'ID_SATELLITE',        'SAT ', &
       'SATELLITEINSTRUMENT', 'INS ', &
       'ZENITHANGLE',         'SZA ', &
       'SOLARZENITHANGLE',    'SUN ', &
       'AZIMUTH',             'AZA ', &
       'SOLARAZIMUTH',        'SAZ ', &
       'FIELDOFVIEW',         'FOV ',  &
       'GEOLOCATIONQUALITY',  'AQF1', &
       'GRANULELEVELQUALITY', 'AQF2', &
       'SCANLEVELQUALITY',    'AQF3' /)

  ! ...for the body table
  integer, parameter :: numBodyMatch = 3
  character(len=lenSqlName) :: bodyKeySqlName = 'ID_OBSERVATION'
  character(len=lenSqlName) :: bodyMatchList(2,numBodyMatch) = (/ &
       'CHANNEL',               'PPP ', &
       'BRIGHTNESSTEMPERATURE', 'VAR ', &
       'CHANQUALITYFLAG',       'QCFL' /)

  ! Also need a dictionary of 'varno' value for each obsDB observation value column
  integer, parameter :: numVarNo = 1
  character(len=lenSqlName) :: varNoList(2,numVarNo) = (/ &
       'BRIGHTNESSTEMPERATURE', '12163' /)

  ! Table names for adding or updating the file
  integer, parameter :: numUpdateTableNames = 4
  character(len=lenSqlName) :: updateTableNames(2,numUpdateTableNames) = (/ &
       'Obs_minus_background',  'OMP ', &
       'Obs_minus_analysis',    'OMA ', &
       'ObsErrorStdDev',        'OER ', &
       'ObsFlags',              'FLG' /)

  ! Other constants
  logical, parameter :: setObsFlagZero = .true.

  ! NAMELIST variables
  logical :: obsDbActive
  integer :: numElemIdList
  integer :: elemIdList(100)

contains

  !--------------------------------------------------------------------------
  ! odbf_setup
  !--------------------------------------------------------------------------
  subroutine odbf_setup()
    !
    ! :Purpose: Read the namelist for obsDB files
    !
    implicit none

    ! locals
    integer           :: nulnam, ierr
    integer, external :: fnom, fclos
    logical, save     :: nmlAlreadyRead = .false.

    namelist /namobsdb/ obsDbActive, numElemIdList, elemIdList

    if ( nmlAlreadyRead ) return

    nmlAlreadyRead = .true.

    ! default values
    obsDbActive = .false.
    numElemIdList = 0
    elemIdList(:) = 0

    if ( .not. utl_isNamelistPresent('NAMOBSDB','./flnml') ) then
      if ( mpi_myid == 0 ) then
        write(*,*) 'odbf_setup: namObsDB is missing in the namelist.'
        write(*,*) '            The default values will be taken.'
      end if
    else
      ! reading namelist variables
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam, nml=namobsdb, iostat=ierr)
      if ( ierr /= 0 ) call utl_abort('odbf_setup: Error reading namelist')
      ierr = fclos(nulnam)
    end if
    if ( mpi_myid == 0 ) write(*, nml=namObsDb)

    if (obsDbActive .and. numElemIdList==0) then
      call utl_abort('odbf_setup: element list is empty')
    end if

  end subroutine odbf_setup

  !--------------------------------------------------------------------------
  ! odbf_isActive
  !--------------------------------------------------------------------------
  function odbf_isActive() result(isActive)
    !
    ! :Purpose: Tell the caller if the namelist indicates obsDB is active
    !
    implicit none

    ! return value
    logical :: isActive

    call odbf_setup()

    isActive = obsDbActive

  end function odbf_isActive

  !--------------------------------------------------------------------------
  ! odbf_readFile
  !--------------------------------------------------------------------------
  subroutine odbf_readFile(obsdat, fileName, familyType, fileIndex)
    !
    ! :Purpose: Read the contents of an obsDB file and put in obsSpaceData
    !
    implicit none

    ! arguments:
    type (struct_obs), intent(inout) :: obsdat
    character(len=*)                 :: fileName
    character(len=*)                 :: familyType
    integer                          :: fileIndex

    ! locals:
    integer :: bodyIndex, bodyIndexBegin, bodyIndexEnd, headIndexBegin
    integer :: headIndexEnd, headIndex, obsRln
    integer :: numBody, numHead, columnIndex
    integer :: headTableIndex, numRowsHeadTable, bodyTableIndex, numRowsBodyTable
    character(len=lenSqlName), allocatable :: headCharSqlNames(:)
    character(len=lenSqlName), allocatable :: headSqlNames(:),    bodySqlNames(:)
    character(len=50),         allocatable :: headCharValues(:,:)
    integer,                   allocatable :: headDateValues(:), headTimeValues(:)
    real(8),                   allocatable :: headValues(:,:), bodyValues(:,:)
    integer(8),                allocatable :: headPrimaryKey(:)
    integer(8),                allocatable :: bodyPrimaryKey(:), bodyHeadKey(:)

    write(*,*)
    write(*,*) 'odbf_readFile: Starting'
    write(*,*)
    write(*,*) 'odbf_readFile: FileName   : ', trim(FileName)
    write(*,*) 'odbf_readFile: FamilyType : ', FamilyType

    !- 0.0 Some initialization
    call ovt_setup(elemIdList(1:numElemIdList))

    !- 1.0 Determine names of columns present in obsDB file

    call odbf_getSqlColumnNames(headCharSqlNames, fileName=trim(fileName), &
                                tableName=headTableName, dataType='varchar')
    call odbf_getSqlColumnNames(headSqlNames, fileName=trim(fileName), &
                                tableName=headTableName, dataType='numeric' )
    call odbf_getSqlColumnNames(bodySqlNames, fileName=trim(fileName), &
                                tableName=bodyTableName, dataType='numeric')

    ! Print all of the column names to the listing
    do columnIndex = 1, size(headCharSqlNames)
      write(*,*) 'odbf_readFile: headCharSqlNames =', columnIndex, &
                 trim(headCharSqlNames(columnIndex))
    end do
    do columnIndex = 1, size(headSqlNames)
      write(*,*) 'odbf_readFile: headSqlNames =', columnIndex, &
                 trim(headSqlNames(columnIndex))
    end do
    do columnIndex = 1, size(bodySqlNames)
      write(*,*) 'odbf_readFile: bodySqlNames   =', columnIndex, &
                 trim(bodySqlNames(columnIndex))
    end do

    !- 1.1 Read most of the contents of the file into local tables

    call odbf_getPrimaryKeys(headPrimaryKey, bodyPrimaryKey, bodyHeadKey, fileName=trim(fileName))

    call odbf_getColumnValues_date(headDateValues, headTimeValues, fileName=trim(fileName), &
                                   tableName=headTableName, sqlColumnName=headDateSqlName)
    call odbf_getColumnValues_char(headCharValues, fileName=trim(fileName), &
                                   tableName=headTableName, sqlColumnNames=headCharSqlNames)
    call odbf_getColumnValues_num (headValues, fileName=trim(fileName), &
                                   tableName=headTableName, sqlColumnNames=headSqlNames)
    call odbf_getColumnValues_num (bodyValues, fileName=trim(fileName), &
                                   tableName=bodyTableName, sqlColumnNames=bodySqlNames)
    numRowsBodyTable = size(bodyValues,1)
    numRowsHeadTable = size(headValues,1)

    ! For debugging, print first 10 rows of each local table to the listing
    do headTableIndex = 1, min(10, numRowsHeadTable)
      write(*,*) 'odbf_readFile: headKeyValues  = ', headPrimaryKey(headTableIndex)
    end do
    do headTableIndex = 1, min(10, numRowsHeadTable)
      write(*,*) 'odbf_readFile: headCharValues = ', headCharValues(headTableIndex,:)
    end do
    do headTableIndex = 1, min(10, numRowsHeadTable)
      write(*,*) 'odbf_readFile: headValues = ', headValues(headTableIndex,:)
    end do
    do bodyTableIndex = 1, min(10, numRowsBodyTable)
      write(*,*) 'odbf_readFile: bodyKeyValues  = ', bodyPrimaryKey(bodyTableIndex)
    end do
    do bodyTableIndex = 1, min(10, numRowsBodyTable)
      write(*,*) 'odbf_readFile: bodyHeadKeyValues  = ', bodyHeadKey(bodyTableIndex)
    end do
    do bodyTableIndex = 1, min(10, numRowsBodyTable)
      write(*,*) 'odbf_readFile: bodyValues = ', bodyValues(bodyTableIndex,:)
    end do

    !- 1.2 Copy values from local tables into obsSpaceData

    ! Starting point for adding rows to obsSpaceData
    bodyIndexBegin = obs_numBody(obsdat) + 1
    headIndexBegin = obs_numHeader(obsdat) + 1

    ! Set the columns related to surface type
    call odbf_setSurfaceType(obsdat, headIndexBegin, fileName=trim(fileName), &
                             tableName=headTableName)

    ! Header date/time values
    call odbf_copyToObsSpaceHeadDate(obsdat, headDateValues, headTimeValues, &
                                     headIndexBegin)

    ! Header-Character values
    call odbf_copyToObsSpaceHeadChar(obsdat, headCharSqlNames, headCharValues, &
                                     headIndexBegin)

    ! Header-Numeric values
    call odbf_copyToObsSpaceHead(obsdat, headSqlNames, headPrimaryKey, &
                                 headValues, headIndexBegin)

    ! Body-Numeric values
    call odbf_copyToObsSpaceBody(obsdat, bodySqlNames, bodyPrimaryKey, bodyHeadKey, &
                                 bodyValues, bodyIndexBegin, headIndexBegin)

    ! Get indexes of last rows added to obsSpaceData
    bodyIndexEnd = obs_numBody(obsdat)
    headIndexEnd = obs_numHeader(obsdat)

    !- 1.3 Set some other quantities in obsSpaceData

    do headIndex = headIndexBegin, headIndexEnd
      call obs_headSet_i(obsdat, OBS_SEN, headIndex, nint(MPC_missingValue_R8))
      call obs_headSet_i(obsdat, OBS_ONM, headIndex, headIndex)
      call obs_headSet_i(obsdat, OBS_OTP, headIndex, fileIndex)
      call obs_headSet_i(obsdat, OBS_IDF, headIndex, fileIndex)
      call obs_setFamily(obsdat, trim(familyType), headIndex)
      if ( headIndex == 1 ) then
        call obs_headSet_i(obsdat, OBS_RLN, headIndex, 1 )
      else
        obsRln = obs_headElem_i(obsdat, OBS_RLN, headIndex-1) + &
                 obs_headElem_i(obsdat, OBS_NLV, headIndex-1)
        call obs_headSet_i(obsdat, OBS_RLN, headIndex, obsRln)
      end if
    end do

    !- 2.0 Additional changes to values after they are in obsSpaceData

    call odbf_adjustValues(obsdat, headIndexBegin, headIndexEnd)

    if ( trim(familyType) /= 'TO' ) then
      call ovt_transformObsValues      (obsdat, headIndexBegin, headIndexEnd )
      call ovt_adjustHumGZ             (obsdat, headIndexBegin, headIndexEnd )
      call obsu_computeVertCoordSurfObs(obsdat, headIndexBegin, headIndexEnd )
    end if

    ! Set a bunch of obsSpaceData body columns to 'missing'
    do bodyIndex = bodyIndexBegin, bodyIndexEnd
      if ( obs_columnActive_RB(obsdat, OBS_OMA) )  call obs_bodySet_r(obsdat, OBS_OMA , bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_OMA0))  call obs_bodySet_r(obsdat, OBS_OMA0, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_OMP) )  call obs_bodySet_r(obsdat, OBS_OMP , bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_OMP6))  call obs_bodySet_r(obsdat, OBS_OMP6, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_OER) )  call obs_bodySet_r(obsdat, OBS_OER , bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_HPHT) ) call obs_bodySet_r(obsdat, OBS_HPHT, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_HAHT))  call obs_bodySet_r(obsdat, OBS_HAHT, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_WORK) ) call obs_bodySet_r(obsdat, OBS_WORK, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_SIGI))  call obs_bodySet_r(obsdat, OBS_SIGI, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_SIGO))  call obs_bodySet_r(obsdat, OBS_SIGO, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_ZHA ))  call obs_bodySet_r(obsdat, OBS_ZHA , bodyIndex, obs_missingValue_R)
    end do

    ! For GP family, initialize OBS_OER to element 15032 (ZTD formal error) 
    ! for all ZTD observations (element 15031)
    if ( trim(familyType) == 'GP') then
      write(*,*)'odbf_readFile: Initializing OBS_OER for GB-GPS ZTD to formal error (ele 15032)'
      call obsu_setGbgpsError(obsdat, headIndexBegin, headIndexEnd )
    end if

    numHead = obs_numHeader(obsdat)
    numBody = obs_numBody(obsdat)
    write(*,*) 'odbf_readFile: after reading file, obs_numHeader = ', numHead
    write(*,*) 'odbf_readFile: after reading file, obs_numBody   = ', numBody
    write(*,*)
    write(*,*) 'odbf_readFile: finished'
    write(*,*)

  end subroutine odbf_readFile

  !--------------------------------------------------------------------------
  ! odbf_updateFile
  !--------------------------------------------------------------------------
  subroutine odbf_updateFile(obsdat, fileName, familyType, fileIndex)
    !
    ! :Purpose: Update the selected quantities in an obsDB file using
    !           values from obsSpaceData. If table does not already
    !           exist, it is created by copying the observation table.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: fileName   
    character(len=*), intent(in)    :: familyType
    integer,          intent(in)    :: fileIndex

    ! locals:
    type(fSQL_STATUS)    :: stat ! sqlite error status
    type(fSQL_DATABASE)  :: db   ! sqlite file handle
    type(fSQL_STATEMENT) :: stmt ! precompiled sqlite statements
    integer(8)           :: obsIdd
    integer              :: obsIdf
    integer              :: tableIndex
    integer              :: headIndex, bodyIndex, bodyIndexBegin, bodyIndexEnd
    integer              :: obsSpaceColIndexSource, fnom, fclos, nulnam, ierr
    real(8)              :: updateValue_r, obsValue
    character(len=4)     :: obsSpaceColumnName
    character(len=lenSqlName) :: tableName, obsValueSqlName
    character(len=3000)  :: query
    logical, save        :: nmlAlreadyRead = .false.

    ! namelist variables
    integer,          save :: numberUpdateItems  ! number of items to use from the list
    character(len=4), save :: updateItemList(15) ! obsSpace column names used to update the file

    namelist/namObsDbUpdate/ numberUpdateItems, updateItemList

    call tmg_start(97,'obdf_updateFile')

    write(*,*)
    write(*,*) 'odbf_updateFile: Starting'
    write(*,*)
    write(*,*) 'odbf_updateFile: FileName   : ', trim(FileName)
    write(*,*) 'odbf_updateFile: FamilyType : ', FamilyType

    if (.not. nmlAlreadyRead) then
      nmlAlreadyRead = .true.

      ! set default values of namelist variables
      updateItemList(:) = ''
      numberUpdateItems = 0

      ! Read the namelist for directives
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam, nml=namObsDbUpdate, iostat=ierr)
      if ( ierr /= 0 ) call utl_abort('odbf_updateFile: Error reading namelist')
      if ( mpi_myid == 0 ) write(*, nml=namObsDbUpdate)
      ierr = fclos(nulnam)

    end if

    ! create tables by copying observation table
    do tableIndex = 1, numberUpdateItems
      tableName = odbf_sqlTableFromObsSpaceName(updateItemList(tableIndex))
      write(*,*) 'odbf_updateFile: creating new sql table for update = ', trim(tableName)

      ! copy bodyTable as template for new table
      call odbf_copySqlTable(fileName, trim(bodyTableName), trim(tableName))

    end do

    ! determine which column(s) will be updated
    obsValueSqlName = odbf_sqlNameFromObsSpaceName('VAR')
    write(*,*) 'odbf_updateFile: column to be updated: ', trim(obsValueSqlName)

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_updateFile: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_updateFile: fSQL_open' )
    end if

    ! update the contents of the new tables
    TABLE: do tableIndex = 1, numberUpdateItems

      ! get obsSpaceData column index for source of updated sql column
      obsSpaceColumnName = updateItemList(tableIndex)
      ierr = clib_toUpper(obsSpaceColumnName)
      obsSpaceColIndexSource = obs_columnIndexFromName(trim(obsSpaceColumnName))

      tableName = odbf_sqlTableFromObsSpaceName(updateItemList(tableIndex))
      write(*,*) 'odbf_updateFile: updating columns in the table: ', trim(tableName)
      write(*,*) 'odbf_updateFile: with contents of obsSpaceData column: ', &
                 trim(obsSpaceColumnName)

      ! create an index for the new table - necessary to speed up the update
      query = 'create index idx_' // trim(tableName) // ' on ' // &
              trim(tableName) // '(' // trim(bodyKeySqlName) // ');'
      write(*,*) 'odbf_copySqlTable: query = ', trim(query)
      call fSQL_do_many( db, query )

      call fSQL_do_many( db,'PRAGMA  synchronous = OFF; PRAGMA journal_mode = OFF;' )

      ! prepare sql update query
      query = 'update ' // trim(tableName) // ' set'
      query = trim(query) // ' ' // trim(obsValueSqlName) // ' = ?'
      query = trim(query)//' where ' // trim(bodyKeySqlName) // ' = ?  ;'
      write(*,*) 'odbf_updateFile: query ---> ', trim(query)

      call fSQL_prepare( db, query , stmt, stat )
      if ( fSQL_error(stat) /= FSQL_OK ) then
        write(*,*) 'odbf_updateFile: fSQL_prepare: ', fSQL_errmsg(stat)
        call utl_abort( 'odbf_updateFile: fSQL_prepare' )
      end if

      call fSQL_begin(db)
      HEADER: do headIndex = 1, obs_numHeader(obsdat)

        obsIdf = obs_headElem_i( obsdat, OBS_IDF, headIndex )
        if ( obsIdf /= fileIndex ) cycle HEADER

        bodyIndexBegin = obs_headElem_i( obsdat, OBS_RLN, headIndex )
        bodyIndexEnd = bodyIndexBegin + &
                       obs_headElem_i( obsdat, OBS_NLV, headIndex ) - 1

        BODY: do bodyIndex = bodyIndexBegin, bodyIndexEnd

          ! do not try to update if the observed value is missing
          obsValue = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          if ( obsValue == obs_missingValue_R ) cycle BODY

          ! update the value, but set to null if it is missing
          updateValue_r = obs_bodyElem_r(obsdat, obsSpaceColIndexSource, bodyIndex)
          if ( updateValue_r == obs_missingValue_R ) then
            call fSQL_bind_param(stmt, PARAM_INDEX=1)  ! sql null values
          else
            call fSQL_bind_param(stmt, PARAM_INDEX=1, REAL8_VAR=updateValue_r)
          end if

          obsIdd  = obs_bodyPrimaryKey( obsdat, bodyIndex )
          call fSQL_bind_param(stmt, PARAM_INDEX=2, INT8_VAR=obsIdd)

          call fSQL_exec_stmt(stmt)

        end do BODY

      end do HEADER

      call fSQL_finalize( stmt )
      call fSQL_commit( db )

      ! drop the index that we just created for the new table
      query = 'drop index idx_' // trim(tableName) // ';'
      write(*,*) 'odbf_copySqlTable: query = ', trim(query)
      call fSQL_do_many( db, query )

    end do TABLE

    ! close the obsDB file
    call fSQL_close( db, stat ) 

    write(*,*)
    write(*,*) 'odbf_updateFile: finished'
    write(*,*)

    call tmg_stop(97)

  end subroutine odbf_updateFile

  !--------------------------------------------------------------------------
  ! odbf_getSqlColumnNames
  !--------------------------------------------------------------------------
  subroutine odbf_getSqlColumnNames(sqlColumnNames, fileName, tableName, dataType)
    !
    ! :Purpose: Read the column names in the obsDB file for the specified table.
    !
    implicit none

    ! arguments:
    character(len=*), allocatable, intent(out) :: sqlColumnNames(:)
    character(len=*),              intent(in)  :: fileName
    character(len=*),              intent(in)  :: tableName
    character(len=*),              intent(in)  :: dataType

    ! locals:
    integer :: numRows, numColumns, rowIndex, ierr
    character(len=100), allocatable :: charValues(:,:)
    character(len=100)       :: dataTypeCriteria
    character(len=3000)      :: query
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle
    type(fSQL_STATEMENT)     :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getSqlColumnNames: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getSqlColumnNames: fSQL_open' )
    end if

    ! read the column names
    if (trim(dataType) == 'varchar') then
      dataTypeCriteria = 'substr(type,1,7)="varchar"'
    else if (trim(dataType) == 'numeric') then
      dataTypeCriteria = 'type="real" or type="REAL" or type="double" or type="DOUBLE" or type="integer" or type="INTEGER"'
    else
      call utl_abort('odbf_getSqlColumnNames: invalid dataType = ' // trim(dataType))
    end if
    query = 'select name from pragma_table_info("' // trim(tableName) // &
            '") where ' // trim(dataTypeCriteria) // ' ;'
    call fSQL_prepare( db, trim(query) , stmt, stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getSqlColumnNames: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getSqlColumnNames: problem with fSQL_get_many')
    end if
    allocate( charValues(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, charValues )

    ! copy to output array and ensure they are upper case
    allocate( sqlColumnNames(numRows) )
    do rowIndex = 1, numRows
      sqlColumnNames(rowIndex) = charValues(rowIndex,1)
      ierr = clib_toUpper(sqlColumnNames(rowIndex))
    end do
    deallocate( charValues )

    ! clean up and close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getSqlColumnNames

  !--------------------------------------------------------------------------
  ! odbf_getPrimaryKeys
  !--------------------------------------------------------------------------
  subroutine odbf_getPrimaryKeys(headPrimaryKey, bodyPrimaryKey, bodyHeadKey, &
                                 fileName)
    !
    ! :Purpose: Read the values from obsDB file for the head and body table
    !           primary keys.
    !
    implicit none

    ! arguments:
    integer(8), allocatable, intent(out) :: headPrimaryKey(:)
    integer(8), allocatable, intent(out) :: bodyPrimaryKey(:)
    integer(8), allocatable, intent(out) :: bodyHeadKey(:)
    character(len=*),        intent(in)  :: fileName

    ! locals:
    integer :: numRows, numColumns
    integer(8), allocatable   :: tempHeadKey(:,:), tempBodyKey(:,:)
    character(len=3000)       :: query
    type(fSQL_STATUS)         :: stat ! sqlite error status
    type(fSQL_DATABASE)       :: db   ! sqlite file handle
    type(fSQL_STATEMENT)      :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getPrimaryKeys: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getPrimaryKeys: fSQL_open' )
    end if

    ! build the sqlite query for the HEADER primary key
    query = 'select ' // trim(headKeySqlName) // ' from ' // &
         trim(headTableName) // ';'
    write(*,*) 'odbf_getPrimaryKeys: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query), stmt, status=stat )
    ! note: "status" not set when getting integers
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_INT8 )
    write(*,*) 'odbf_getPrimaryKeys: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( headPrimaryKey(numRows) )
    allocate( tempHeadKey(numRows,1) )
    call fSQL_fill_matrix( stmt, tempHeadKey )
    headPrimaryKey(:) = tempHeadKey(:,1)
    deallocate(tempHeadKey)

    ! build the sqlite query for the BODY primary key
    query = 'select ' // trim(bodyKeySqlName) // ' from ' // &
            trim(bodyTableName) // ';'
    write(*,*) 'odbf_getPrimaryKeys: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query) , stmt, status=stat )
    ! note: "status" not set when getting integers
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_INT8 )
    write(*,*) 'odbf_getPrimaryKeys: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( bodyPrimaryKey(numRows) )
    allocate( tempBodyKey(numRows,1) )
    call fSQL_fill_matrix( stmt, tempBodyKey )
    bodyPrimaryKey(:) = tempBodyKey(:,1)
    deallocate(tempBodyKey)

    ! build the sqlite query for the BODY-HEAD key
    query = 'select ' // trim(headKeySqlName) // ' from ' // &
            trim(bodyTableName) // ';'
    write(*,*) 'odbf_getPrimaryKeys: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query) , stmt, status=stat )
    ! note: "status" not set when getting integers
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_INT8 )
    write(*,*) 'odbf_getPrimaryKeys: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( bodyHeadKey(numRows) )
    allocate( tempBodyKey(numRows,1) )
    call fSQL_fill_matrix( stmt, tempBodyKey )
    bodyHeadKey(:) = tempBodyKey(:,1)
    deallocate(tempBodyKey)

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getPrimaryKeys

  !--------------------------------------------------------------------------
  ! odbf_getColumnValues_date
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnValues_date(columnDateValues, columnTimeValues, fileName, &
                                       tableName, sqlColumnName)
    !
    ! :Purpose: Read the column values from obsDB file for the specified table
    !           and column names.
    !
    implicit none

    ! arguments:
    integer, allocatable, intent(out) :: columnDateValues(:)
    integer, allocatable, intent(out) :: columnTimeValues(:)
    character(len=*),     intent(in)  :: sqlColumnName
    character(len=*),     intent(in)  :: fileName
    character(len=*),     intent(in)  :: tableName

    ! locals:
    integer              :: numRows, numColumns, rowIndex
    character(len=20), allocatable :: columnValuesStr(:,:)
    character(len=3000)  :: query
    type(fSQL_STATUS)    :: stat ! sqlite error status
    type(fSQL_DATABASE)  :: db   ! sqlite file handle
    type(fSQL_STATEMENT) :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValues_date: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getColumnValues_date: fSQL_open' )
    end if

    ! Get the date and time

    ! build the sqlite query
    query = "select strftime('%Y%m%d'," // trim(sqlColumnName) // "), strftime('%H%M'," // trim(sqlColumnName) // ")"
    query = trim(query) // ' from ' // trim(tableName) // ';'
    write(*,*) 'odbf_getColumnValues_date: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query), stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValues_date: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getColumnValues_date: problem with fSQL_get_many')
    end if
    write(*,*) 'odbf_getColumnValues_date: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValuesStr(numRows,2) )
    call fSQL_fill_matrix( stmt, columnValuesStr )
    allocate( columnDateValues(numRows) )
    allocate( columnTimeValues(numRows) )
    do rowIndex = 1, numRows
      read(columnValuesStr(rowIndex,1),*) columnDateValues(rowIndex)
      read(columnValuesStr(rowIndex,2),*) columnTimeValues(rowIndex)
    end do

    deallocate(columnValuesStr)

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getColumnValues_date

  !--------------------------------------------------------------------------
  ! odbf_setSurfaceType
  !--------------------------------------------------------------------------
  subroutine odbf_setSurfaceType(obsdat, headIndexBegin, fileName, tableName)
    !
    ! :Purpose: Set the surface type based on lat-lon and some external mask files.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer,          intent(in)    :: headIndexBegin
    character(len=*), intent(in)    :: fileName
    character(len=*), intent(in)    :: tableName

    ! locals:
    integer              :: numRows, numColumns, headTableIndex, headIndex
    integer, allocatable :: columnValues(:,:)
    character(len=3000)  :: query
    type(fSQL_STATUS)    :: stat ! sqlite error status
    type(fSQL_DATABASE)  :: db   ! sqlite file handle
    type(fSQL_STATEMENT) :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_setSurfaceType: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_setSurfaceType: fSQL_open' )
    end if

    ! build the sqlite query
    query = 'select mask_mer(lat,lon) from ' // trim(tableName) // ';'
    write(*,*) 'odbf_setSurfaceType: query ---> ', trim(query)

    ! read the values from the query result
    call fSQL_prepare( db, trim(query), stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_INT, status=stat )
    write(*,*) 'odbf_setSurfaceType: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValues(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnValues )

    ! set the values of STYP and TTYP
    do headTableIndex = 1, numRows
      headIndex = headTableIndex + headIndexBegin - 1
      call obs_headSet_i(obsdat, OBS_STYP, headIndex, columnValues(headTableIndex,1))
      call obs_headSet_i(obsdat, OBS_TTYP, headIndex, -1) ! Not sea ice
    end do

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_setSurfaceType

  !--------------------------------------------------------------------------
  ! odbf_getColumnValues_char
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnValues_char(columnValues, fileName, tableName, &
                                       sqlColumnNames)
    !
    ! :Purpose: Read the column values from obsDB file for the specified table
    !           and column names.
    !
    implicit none

    ! arguments:
    character(len=50), allocatable, intent(out) :: columnValues(:,:)
    character(len=*),               intent(in)  :: sqlColumnNames(:)
    character(len=*),               intent(in)  :: fileName
    character(len=*),               intent(in)  :: tableName

    ! locals:
    integer :: numRows, numColumns, columnIndex
    character(len=3000)      :: query
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle
    type(fSQL_STATEMENT)     :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValues_char: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getColumnValues_char: fSQL_open' )
    end if

    ! build the sqlite query
    query = 'select'
    numColumns = size(sqlColumnNames)
    do columnIndex = 1, numColumns
      query = trim(query) // ' ' // trim(sqlColumnNames(columnIndex))
      if (columnIndex < numColumns) query = trim(query) // ','
    end do
    query = trim(query) // ' from ' // trim(tableName) // ';'
    write(*,*) 'odbf_getColumnValues_char: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query), stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValues_char: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getColumnValues_char: problem with fSQL_get_many')
    end if
    write(*,*) 'odbf_getColumnValues_char: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValues(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnValues )

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getColumnValues_char

  !--------------------------------------------------------------------------
  ! odbf_getColumnValues_num
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnValues_num(columnValues, fileName, tableName, &
                                      sqlColumnNames)
    !
    ! :Purpose: Read the column values from obsDB file for the specified table
    !           and column names.
    !
    implicit none

    ! arguments:
    real(8), allocatable, intent(out) :: columnValues(:,:)
    character(len=*),     intent(in)  :: sqlColumnNames(:)
    character(len=*),     intent(in)  :: fileName
    character(len=*),     intent(in)  :: tableName

    ! locals:
    integer :: numRows, numColumns, columnIndex
    character(len=3000)       :: query
    type(fSQL_STATUS)         :: stat ! sqlite error status
    type(fSQL_DATABASE)       :: db   ! sqlite file handle
    type(fSQL_STATEMENT)      :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValues_num: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getColumnValues_num: fSQL_open' )
    end if

    ! build the sqlite query
    query = 'select'
    numColumns = size(sqlColumnNames)
    do columnIndex = 1, numColumns
      query = trim(query) // ' ' // trim(sqlColumnNames(columnIndex))
      if (columnIndex < numColumns) query = trim(query) // ','
    end do
    query = trim(query) // ' from ' // trim(tableName)
    query = trim(query) // ';'
    write(*,*) 'odbf_getColumnValues_num: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query) , stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, mode=FSQL_REAL8, &
                        real8_missing=MPC_missingValue_R8, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValues_num: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getColumnValues_num: problem with fSQL_get_many')
    end if
    write(*,*) 'odbf_getColumnValues_num: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValues(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnValues )

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getColumnValues_num

  !--------------------------------------------------------------------------
  ! odbf_copyToObsSpaceHeadChar
  !--------------------------------------------------------------------------
  subroutine odbf_copyToObsSpaceHeadChar(obsdat, headCharSqlNames, &
                                         headCharValues, headIndexBegin)
    !
    ! :Purpose: Copy character string values from a local table into
    !           obsSpaceData header rows. Currently, only the STATION ID
    !           and OBS_ITY (i.e. codeType from character sql column
    !           containing obs type name).
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: headCharSqlNames(:)
    character(len=*), intent(in)    :: headCharValues(:,:)
    integer,          intent(in)    :: headIndexBegin

    ! locals:
    character(len=lenSqlName) :: stIdSqlName, codeTypeSqlName
    integer :: columnIndex, headTableIndex, headIndex
    integer :: numRowsHeadTable, codeType

    numRowsHeadTable = size(headCharValues,1)

    ! Set the STATION ID
    stIdSqlName = odbf_sqlNameFromObsSpaceName('STID')
    columnIndex = utl_findloc(headCharSqlNames(:), trim(stIdSqlName))
    if (columnIndex == 0) then
      call utl_abort('odbf_copyToObsSpaceHeadChar: Station ID column not found in sql table')
    end if
    do headTableIndex = 1, numRowsHeadTable
      headIndex = headTableIndex + headIndexBegin - 1
      if (headTableIndex == 1) then
        write(*,*) 'odbf_copyToObsSpaceHeadChar: set header char column   : ', &
                   trim(headCharSqlNames(columnIndex))
      end if
      call obs_set_c(obsdat, 'STID', &
                     headIndex, headCharValues(headTableIndex,columnIndex))
    end do

    ! Set the codeType (obs_ity) from a character column containing the obs type
    codeTypeSqlName = odbf_sqlNameFromObsSpaceName('ITY')
    columnIndex = utl_findloc(headCharSqlNames(:), trim(codeTypeSqlName))
    if (columnIndex == 0) then
      call utl_abort('odbf_copyToObsSpaceHeadChar: Obs type column not found in sql table')
    end if
    do headTableIndex = 1, numRowsHeadTable
      headIndex = headTableIndex + headIndexBegin - 1
      if (headTableIndex == 1) then
        write(*,*) 'odbf_copyToObsSpaceHeadChar: set header char column   : ', &
                   trim(headCharSqlNames(columnIndex))
      end if
      codeType = codtyp_get_codtyp(trim(headCharValues(headTableIndex,columnIndex)))
      if (codeType == -1) then
        write(*,*) 'odbf_copyToObsSpaceHeadChar: obs type =', trim(headCharValues(headTableIndex,columnIndex))
        call utl_abort('odbf_copyToObsSpaceHeadChar: codtyp for this obs type not found') 
      end if
      call obs_headSet_i(obsdat, OBS_ITY, headIndex, codeType)
    end do

  end subroutine odbf_copyToObsSpaceHeadChar

  !--------------------------------------------------------------------------
  ! odbf_copyToObsSpaceHeadDate
  !--------------------------------------------------------------------------
  subroutine odbf_copyToObsSpaceHeadDate(obsdat, headDateValues, headTimeValues, headIndexBegin)
    !
    ! :Purpose: Copy date values from a local table into
    !           obsSpaceData header rows.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer,          intent(in)    :: headDateValues(:)
    integer,          intent(in)    :: headTimeValues(:)
    integer,          intent(in)    :: headIndexBegin

    ! locals:
    integer :: headTableIndex, headIndex
    integer :: numRowsHeadTable

    numRowsHeadTable = size(headDateValues,1)

    do headTableIndex = 1, numRowsHeadTable
      headIndex = headTableIndex + headIndexBegin - 1
      if (headTableIndex == 1) then
        write(*,*) 'odbf_copyToObsSpaceHeadDate: set header date/time column'
      end if
      call obs_headSet_i(obsdat, obs_dat, headIndex, headDateValues(headTableIndex))
      call obs_headSet_i(obsdat, obs_etm, headIndex, headTimeValues(headTableIndex))
    end do

  end subroutine odbf_copyToObsSpaceHeadDate

  !--------------------------------------------------------------------------
  ! odbf_copyToObsSpaceHead
  !--------------------------------------------------------------------------
  subroutine odbf_copyToObsSpaceHead(obsdat, headSqlNames, headPrimaryKey, &
                                     headValues, headIndexBegin)
    !
    ! :Purpose: Copy real and integer values from a local table into
    !           obsSpaceData header rows.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: headSqlNames(:)
    integer(8),       intent(in)    :: headPrimaryKey(:)
    real(8),          intent(in)    :: headValues(:,:)
    integer,          intent(in)    :: headIndexBegin

    ! locals:
    integer :: columnIndex, matchIndex, headTableIndex, headIndex
    integer :: numRowsHeadTable, obsColumnIndex

    numRowsHeadTable = size(headValues,1)

    write(*,*) 'odbf_copyToObsSpaceHead: set header primary key'
    do headTableIndex = 1, numRowsHeadTable
      headIndex = headTableIndex + headIndexBegin - 1
      call obs_setHeadPrimaryKey( obsdat,  headIndex, headPrimaryKey(headTableIndex))
    end do

    column_loop: do columnIndex = 1, size(headSqlNames)
      matchIndex = utl_findloc(headMatchList(sqlColIndex,:), headSqlNames(columnIndex))

      if (matchIndex == 0) then

        write(*,*) 'odbf_copyToObsSpaceHead: unknown column name      : ', &
                   trim(headSqlNames(columnIndex))

      else

        obsColumnIndex = obs_columnIndexFromName(trim(headMatchList(obsColIndex,matchIndex)))

        headTable_loop: do headTableIndex = 1, numRowsHeadTable
          headIndex = headTableIndex + headIndexBegin - 1

          if (obs_columnDataType(obsColumnIndex) == 'real') then
            if ( obs_columnActive_RH(obsdat, obsColumnIndex) ) then
              if (headTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceHead: set header real column   : ', trim(headSqlNames(columnIndex))
              end if
              call obs_headSet_r(obsdat, obsColumnIndex, &
                                 headIndex, real(headValues(headTableIndex,columnIndex),pre_obsReal))
            end if
          else if (obs_columnDataType(obsColumnIndex) == 'integer') then
            if ( obs_columnActive_IH(obsdat, obsColumnIndex) ) then
              if (headTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceHead: set header integer column: ', trim(headSqlNames(columnIndex))
              end if
              call obs_headSet_i(obsdat, obsColumnIndex, &
                                 headIndex, nint(headValues(headTableIndex,columnIndex)))
            end if
          else
            call utl_abort('odbf_copyToObsSpaceHead: unknown data type for obs header column')
          end if

        end do headTable_loop

      end if

    end do column_loop

  end subroutine odbf_copyToObsSpaceHead

  !--------------------------------------------------------------------------
  ! odbf_copyToObsSpaceBody
  !--------------------------------------------------------------------------
  subroutine odbf_copyToObsSpaceBody(obsdat, bodySqlNames, &
                                     bodyPrimaryKey, bodyHeadKey, bodyValues, &
                                     bodyIndexBegin, headIndexBegin)
    !
    ! :Purpose: Copy real and integer values from a local table into
    !           obsSpaceData body rows.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: bodySqlNames(:)
    integer(8),       intent(in)    :: bodyPrimaryKey(:)
    integer(8),       intent(in)    :: bodyHeadKey(:)
    real(8),          intent(in)    :: bodyValues(:,:)
    integer,          intent(in)    :: bodyIndexBegin
    integer,          intent(in)    :: headIndexBegin

    ! locals:
    character(len=lenSqlName) :: obsValueSqlName
    integer :: columnIndex, matchIndex, bodyTableIndex, bodyIndex, headIndex
    integer :: numRowsBodyTable, bodyColumnIndexObsValue, obsNlv, obsVarNo
    integer(8) :: lastHeadKey
    integer, allocatable :: matchIndexVec(:)
    integer, allocatable :: obsColumnIndex(:)

    numRowsBodyTable = size(bodyValues,1)

    ! initialize some arrays to save lots of time in the main loops
    allocate(matchIndexVec(size(bodySqlNames)))
    do columnIndex = 1, size(bodySqlNames)
      matchIndexVec(columnIndex) = utl_findloc(bodyMatchList(sqlColIndex,:), &
                                               bodySqlNames(columnIndex))
    end do
    allocate(obsColumnIndex(numBodyMatch))
    do matchIndex = 1, numBodyMatch
      obsColumnIndex(matchIndex) = obs_columnIndexFromName(bodyMatchList(obsColIndex,matchIndex))
    end do

    ! figure out column index for observation value (OBS_VAR)
    obsValueSqlName = odbf_sqlNameFromObsSpaceName('VAR')
    bodyColumnIndexObsValue = utl_findloc(bodySqlNames(:), obsValueSqlName)
    if (bodyColumnIndexObsValue == 0) then
      write(*,*) 'odbf_copyToObsSpaceBody: obsValueSqlName = ', trim(obsValueSqlName)
      call utl_abort('odbf_copyToObsSpaceBody: column with obs value not present')
    end if
    lastHeadKey = 0
    headIndex = headIndexBegin - 1
    bodyIndex = bodyIndexBegin - 1

    ! determine varNo for the observation value
    obsVarNo = odbf_varNoFromSqlName(obsValueSqlName)
    write(*,*) 'odbf_copyToObsSpaceBody: obsVarNo = ', obsVarNo

    bodyIndex_loop: do bodyTableIndex = 1, numRowsBodyTable
      bodyIndex = bodyIndex + 1

      ! check if obs value is null/missing
      if ( bodyValues(bodyTableIndex,bodyColumnIndexObsValue) == MPC_missingValue_R8 ) then
        bodyIndex = bodyIndex - 1
        cycle bodyIndex_loop
      end if

      ! count number of body rows for each header row (OBS_NLV)
      if ( bodyHeadKey(bodyTableIndex) /= lastHeadKey ) then
        headIndex = headIndex + 1
        call obs_headSet_i(obsdat, OBS_NLV, headIndex, 0)
        lastHeadKey = bodyHeadKey(bodyTableIndex)
      end if
      obsNLV = obs_headElem_i(obsdat, OBS_NLV, headIndex)
      call obs_headSet_i(obsdat, OBS_NLV, headIndex, obsNLV + 1)

      ! check that the primary key for header table matches the value in the body table
      if ( obs_headPrimaryKey( obsdat, headIndex ) /= &
           bodyHeadKey(bodyTableIndex) ) then
        write(*,*) 'odbf_copyToObsSpaceBody: primary key in HEADER table = ', &
                   obs_headPrimaryKey( obsdat, headIndex )
        write(*,*) 'odbf_copyToObsSpaceBody: same key in BODY table      = ', &
                   bodyHeadKey(bodyTableIndex)
        call utl_abort('odbf_copyToObsSpaceBody: Primary key of HEADER table not equal ' // &
                       'to value in BODY table')
      end if

      ! copy body primary key to obsSpaceData
      if (bodyTableIndex == 1) then
        write(*,*) 'odbf_copyToObsSpaceBody: set body primary key'
      end if
      call obs_setBodyPrimaryKey(obsdat, bodyIndex, bodyPrimaryKey(bodyTableIndex))

      ! set the varNo, assuming only 1 observed value
      call obs_bodySet_i(obsdat, OBS_VNM, bodyIndex, obsVarNo)

      ! copy real and integer values into obsSpaceData
      columnIndex_loop: do columnIndex = 1, size(bodySqlNames)
        matchIndex = matchIndexVec(columnIndex)

        if (matchIndex == 0) then
          if (bodyTableIndex == 1) then
            write(*,*) 'odbf_copyToObsSpaceBody: unknown column name    : ', &
                       trim(bodySqlNames(columnIndex))
          end if
        else
          if (obs_columnDataType(obsColumnIndex(matchIndex)) == 'real') then
            ! real values
            if ( obs_columnActive_RB(obsdat, obsColumnIndex(matchIndex)) ) then
              if (bodyTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceBody: set body real column   : ', trim(bodySqlNames(columnIndex))
              end if
              call obs_bodySet_r(obsdat, obsColumnIndex(matchIndex), &
                                 bodyIndex, real(bodyValues(bodyTableIndex,columnIndex),pre_obsReal))
            end if
          else if (obs_columnDataType(obsColumnIndex(matchIndex)) == 'integer') then
            ! integer values
            if ( obs_columnActive_IB(obsdat, obsColumnIndex(matchIndex)) ) then
              if (bodyTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceBody: set body integer column: ', trim(bodySqlNames(columnIndex))
              end if
              call obs_bodySet_i(obsdat, obsColumnIndex(matchIndex), &
                                 bodyIndex, nint(bodyValues(bodyTableIndex,columnIndex)))
            end if
          else
            call utl_abort('odbf_copyToObsSpaceBody: unknown data type for obs body column')
          end if
        end if

      end do columnIndex_loop

    end do bodyIndex_loop

    deallocate(matchIndexVec)
    deallocate(obsColumnIndex)

  end subroutine odbf_copyToObsSpaceBody

  !--------------------------------------------------------------------------
  ! odbf_adjustValues
  !--------------------------------------------------------------------------
  subroutine odbf_adjustValues(obsdat, headIndexBegin, headIndexEnd)
    !
    ! :Purpose: Adjust units and other minor modifications of some
    !           obsSpaceData columns after transfer from sqlite files
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer,          intent(in)    :: headIndexBegin
    integer,          intent(in)    :: headIndexEnd

    ! locals:
    integer :: headIndex, bodyIndexStart, bodyIndexEnd, bodyIndex
    integer :: instrument, obsSat, codeType, sensor
    real(8) :: obsLon, obsLat, surfEmiss
    character(len=2) :: obsFamily

    do headIndex = headIndexBegin, headIndexEnd

      obsFamily = obs_getfamily( obsdat, headIndex )
      bodyIndexStart = obs_headElem_i(obsdat, OBS_RLN, headIndex)
      bodyIndexEnd   = obs_headElem_i(obsdat, OBS_NLV, headIndex) + &
                       bodyIndexStart - 1

      ! Convert lon-lat from degrees to radians

      obsLon = obs_headElem_r( obsdat, OBS_LON, headIndex )
      obsLat = obs_headElem_r( obsdat, OBS_LAT, headIndex )

      if ( obsLon < 0.0D0 ) obsLon = obsLon + 360.0D0
      obsLon = obsLon * MPC_RADIANS_PER_DEGREE_R8
      obsLat = obsLat * MPC_RADIANS_PER_DEGREE_R8

      call obs_headSet_r(obsdat, OBS_LON, headIndex, real(obsLon,pre_obsReal))
      call obs_headSet_r(obsdat, OBS_LAT, headIndex, real(obsLat,pre_obsReal))

      ! Set global and observation flags to zero if specified

      if (setObsFlagZero) then
        call obs_headSet_i(obsdat, OBS_ST1, headIndex, 0)
        do bodyIndex = bodyIndexStart, bodyIndexEnd
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, 0)
        end do
      end if

      ! Various adjustment for radiance observations

      if ( obsFamily == 'TO' ) then

        instrument = obs_headElem_i( obsdat, OBS_INS, headIndex )
        obsSat     = obs_headElem_i( obsdat, OBS_SAT, headIndex )
        codeType   = obs_headElem_i( obsdat, OBS_ITY, headIndex )
        sensor     = obs_headElem_i( obsdat, OBS_SEN, headIndex )

        ! set sensor to missing if not amsua/b, mhs or atms

        if ( codeType /= codtyp_get_codtyp('amsua') .and. &
             codeType /= codtyp_get_codtyp('amsub') .and. &
             codeType /= codtyp_get_codtyp('mhs') .and. &
             codeType /= codtyp_get_codtyp('atms') ) then
          sensor = nint(MPC_missingValue_R8)
        end if

        ! modify OBS_SAT, OBS_INS and OBS_SEN

        if ( instrument == 420 ) obsSat  = 784
        if ( codeType == 202 .and. instrument == 620 ) instrument = 2046
        if ( sensor == nint(MPC_missingValue_R8) ) then
          sensor = 0
          if (instrument == nint(MPC_missingValue_R8) ) instrument = 0
        else
          instrument = obsu_cvt_obs_instrum(sensor)
        end if
        call obs_headSet_i(obsdat, OBS_INS, headIndex, instrument)
        call obs_headSet_i(obsdat, OBS_SAT, headIndex, obsSat)
        call obs_headSet_i(obsdat, OBS_SEN, headIndex, sensor)

        ! change units for surface emissivity

        do bodyIndex = bodyIndexStart, bodyIndexEnd
          surfEmiss = obs_bodyElem_r( obsdat, OBS_SEM, bodyIndex )
          surfEmiss = surfEmiss * 0.01D0
          call obs_bodySet_r(obsdat, OBS_SEM, bodyIndex, real(surfEmiss, pre_obsReal))
        end do
        
      end if
      
    end do

  end subroutine odbf_adjustValues

  !--------------------------------------------------------------------------
  ! odbf_sqlNameFromObsSpaceName
  !--------------------------------------------------------------------------
  function odbf_sqlNameFromObsSpaceName(obsSpaceName) result(sqlName)
    !
    ! :Purpose: Return the corresponding sql file column name for a
    !           given obsSpaceData column name from the matching
    !           tables.
    !
    implicit none

    ! arguments:
    character(len=*), intent(in) :: obsSpaceName
    character(len=lenSqlName)    :: sqlName

    ! locals:
    integer                   :: matchIndex

    ! first try the body matching list
    matchIndex = utl_findloc(bodyMatchList(obsColIndex,:), trim(obsSpaceName))
    if (matchIndex > 0) then
      sqlName = bodyMatchList(sqlColIndex,matchIndex)
      return
    end if

    ! now try the header matching list
    matchIndex = utl_findloc(headMatchList(obsColIndex,:), trim(obsSpaceName))
    if (matchIndex > 0) then
      sqlName = headMatchList(sqlColIndex,matchIndex)
      return
    end if

    ! not found in either list, abort
    write(*,*) 'odbf_sqlNameFromObsSpaceName: requested obsSpace name = ', trim(obsSpaceName)
    call utl_abort('odbf_sqlNameFromObsSpaceName: obsSpace name not found in matching list')
    
  end function odbf_sqlNameFromObsSpaceName

  !--------------------------------------------------------------------------
  ! odbf_sqlNameFromObsSpaceName
  !--------------------------------------------------------------------------
  function odbf_sqlTableFromObsSpaceName(obsSpaceName) result(tableName)
    !
    ! :Purpose: Return the corresponding sql file table name for a
    !           given obsSpaceData column name from the matching
    !           tables.
    !
    implicit none

    ! arguments:
    character(len=*), intent(in) :: obsSpaceName
    character(len=lenSqlName)    :: tableName

    ! locals:
    integer                   :: matchIndex

    ! first try the body matching list
    matchIndex = utl_findloc(updateTableNames(obsColIndex,:), trim(obsSpaceName))
    if (matchIndex > 0) then
      tableName = updateTableNames(sqlTabIndex,matchIndex)
      return
    end if

    ! not found, abort
    write(*,*) 'odbf_sqlTableFromObsSpaceName: requested obsSpace name = ', trim(obsSpaceName)
    call utl_abort('odbf_sqlTableFromObsSpaceName: obsSpace name not found in updateTableNames')
    
  end function odbf_sqlTableFromObsSpaceName

  !--------------------------------------------------------------------------
  ! odbf_varNoFromSqlName
  !--------------------------------------------------------------------------
  function odbf_varNoFromSqlName(sqlName) result(varNo)
    !
    ! :Purpose: Return the bufr element id number from the corresponding
    !           sql file column name of an observed value.
    !
    implicit none

    ! arguments:
    character(len=*), intent(in) :: sqlName
    integer                      :: varNo

    ! locals:
    integer                   :: matchIndex
    character(len=10)         :: varNoStr
    
    matchIndex = utl_findloc(varNoList(sqlColIndex,:), trim(sqlName))
    if (matchIndex > 0) then
      varNoStr = varNoList(varNoColIndex,matchIndex)
      read(varNoStr,*) varNo
      return
    end if

    write(*,*) 'odbf_varNoFromSqlName: requested sqlName = ', trim(sqlName)
    call utl_abort('odbf_varNoFromSqlName: not found in varNo list')
    
  end function odbf_varNoFromSqlName

  !--------------------------------------------------------------------------
  ! odbf_copySqlTable
  !--------------------------------------------------------------------------
  subroutine odbf_copySqlTable(fileName, tableNameSource, tableNameNew)
    !
    ! :Purpose: Copy an existing sql table to a new table with a different name
    !
    implicit none

    ! arguments:
    character(len=*),              intent(in)  :: fileName
    character(len=*),              intent(in)  :: tableNameSource
    character(len=*),              intent(in)  :: tableNameNew

    ! locals:
    integer :: numRows, numColumns, rowIndex
    character(len=100), allocatable :: tableValues(:,:)
    character(len=3000)      :: query
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle
    type(fSQL_STATEMENT)     :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getSqlColumnNames: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getSqlColumnNames: fSQL_open' )
    end if

    ! read the column names
    query = 'select name, type from pragma_table_info("' // trim(tableNameSource) // '");'
    call fSQL_prepare( db, trim(query) , stmt, stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_copySqlTable: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_copySqlTable: problem with fSQL_get_many')
    end if
    allocate( tableValues(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, tableValues )

    ! create the new table
    query = 'create table ' // trim(tableNameNew) // ' (' // new_line('A')
    do rowIndex = 1, numRows
      query = trim(query) // '  ' // trim(tableValues(rowIndex,1)) // ' ' // trim(tableValues(rowIndex,2))
      if (rowIndex < numRows) query = trim(query) // ', '
      query = trim(query) // new_line('A')
    end do
    query = trim(query) // ');'
    write(*,*) 'odbf_copySqlTable: query = ', trim(query)
    call fSQL_do_many( db, query )

    ! copy values from original to new table
    query = 'insert into ' // trim(tableNameNew) // ' select * from ' // trim(tableNameSource) // ';'
    write(*,*) 'odbf_copySqlTable: query = ', trim(query)
    call fSQL_do_many( db, query )

    ! clean up and close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 
    deallocate( tableValues )

  end subroutine odbf_copySqlTable

end module obsdbFiles_mod

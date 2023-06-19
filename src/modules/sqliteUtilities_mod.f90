
module sqliteUtilities_mod
  ! MODULE sqliteUtilities (prefix='sqlu' category='3. Observation input/output')
  !
  ! :Purpose: A place to collect utilities for SQLite files.

use fSQLite
use clibInterfaces_mod
use obsSpaceData_mod
use midasMpi_mod
use utilities_mod
use mathPhysConstants_mod

implicit none
save
private
public :: sqlu_sqlColumnExists, sqlu_sqlTableExists, sqlu_getSqlColumnNames
public :: sqlu_query, sqlu_handleError
public :: sqlu_getColumnValuesNum, sqlu_getColumnValuesDateStr, sqlu_getColumnValuesChar
public :: sqlu_getInitialIdObsData

! Arrays used to match SQLite column names with obsSpaceData column names
integer, parameter :: lenSqlName = 60

contains
   
  !--------------------------------------------------------------------------
  ! sqlu_sqlColumnExists
  !--------------------------------------------------------------------------
  function sqlu_sqlColumnExists(fileName, tableName, columnName) result(columnExists)
    !
    ! :Purpose: Check if a column exists in the sqlite file/table
    !
    implicit none

    ! Arguments:
    character(len=*),              intent(in)  :: fileName
    character(len=*),              intent(in)  :: tableName
    character(len=*),              intent(in)  :: columnName
    ! Result:
    logical                                    :: columnExists

    ! Locals:
    integer                     :: ierr
    character(len=3000)         :: query, sqliteOutput
    character(len=lenSqlName)   :: upperColumnName
    type(fSQL_STATUS)           :: stat ! sqlite error status
    type(fSQL_DATABASE)         :: db   ! sqlite file handle
    logical, parameter          :: debug = .false.
    character(len=*), parameter :: myName = 'sqlu_sqlColumnExists'
    
    ! open the SQLite file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myName//fSQL_errmsg(stat) )
    end if

    upperColumnName = trim(columnName)
    ierr = clib_toUpper(upperColumnName)

    query = "select count(*) from pragma_table_info('" // &
            trim(tableName) // "') where upper(name)='" // trim(upperColumnName) // "';"
    if (debug) write(*,*) myName//': query = ', trim(query)

    sqliteOutput = sqlu_query(db,trim(query))
    if (debug) write(*,*) myName//': output = XXX' // trim(sqliteOutput) // 'XXX'
    columnExists = (trim(sqliteOutput) /= '0')

    ! close the sqlite file
    call fSQL_close( db, stat )

  end function sqlu_sqlColumnExists

  !--------------------------------------------------------------------------
  ! sqlu_sqlTableExists
  !--------------------------------------------------------------------------
  function sqlu_sqlTableExists(fileName, tableName) result(tableExists)
    !
    ! :Purpose: Check if a table exists in the sqlite file
    !
    implicit none

    ! Arguments:
    character(len=*),              intent(in)  :: fileName
    character(len=*),              intent(in)  :: tableName
    ! Result:
    logical                                    :: tableExists

    ! Locals:
    integer                     :: ierr
    logical                     :: finished
    character(len=3000)         :: query, sqliteOutput
    character(len=lenSqlName)   :: upperTableName
    type(fSQL_STATUS)           :: stat ! sqlite error status
    type(fSQL_DATABASE)         :: db   ! sqlite file handle
    type(fSQL_STATEMENT)        :: stmt ! precompiled sqlite statements
    logical, parameter          :: debug = .false.
    character(len=*), parameter :: myName = 'sqlu_sqlTableExists'

    ! open the sqlite file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myName//': fSQL_open '//fSQL_errmsg(stat) )
    end if

    upperTableName = trim(tableName)
    ierr = clib_toUpper(upperTableName)

    query = "select upper(name) as uppername from sqlite_master where " // &
            "type='table' and uppername='" // trim(upperTableName) // "';"
    if (debug) write(*,*) myName//': query = ', trim(query)

    call fSQL_prepare( db, trim(query), stmt, stat)
    finished = .false.
    call fSQL_get_row( stmt, finished )
    call fSQL_get_column( stmt, COL_INDEX = 1, CHAR_VAR = sqliteOutput )
    call fSQL_get_row( stmt, finished )
    call fSQL_finalize( stmt )
    if (debug) write(*,*) myName//': output = XXX' // trim(sqliteOutput) // 'XXX'
    tableExists = (trim(sqliteOutput) == trim(upperTableName))

    ! close the sqlite file
    call fSQL_close( db, stat ) 

  end function sqlu_sqlTableExists

  !--------------------------------------------------------------------------
  ! sqlu_getSqlColumnNames
  !--------------------------------------------------------------------------
  subroutine sqlu_getSqlColumnNames(sqlColumnNames, fileName, tableName, dataType)
    !
    ! :Purpose: Read the column names in the sqlite file for the specified table.
    !
    implicit none

    ! Arguments:
    character(len=*), allocatable, intent(out) :: sqlColumnNames(:)
    character(len=*),              intent(in)  :: fileName
    character(len=*),              intent(in)  :: tableName
    character(len=*),              intent(in)  :: dataType

    ! Locals:
    integer :: numRows, numColumns, rowIndex, ierr
    character(len=100), allocatable :: charValues(:,:)
    character(len=200)       :: dataTypeCriteria
    character(len=3000)      :: query
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle
    type(fSQL_STATEMENT)     :: stmt ! precompiled sqlite statements
    character(len=*), parameter :: myName = 'sqlu_getSqlColumnNames'

    ! open the sqlite file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myName//': fSQL_open '//fSQL_errmsg(stat) )
    end if

    ! read the column names
    if (trim(dataType) == 'varchar') then
      dataTypeCriteria = 'substr(type,1,7)="varchar"'
    else if (trim(dataType) == 'numeric') then
      dataTypeCriteria = 'type="real" or type="REAL" or type="double" or ' // &
                         'type="DOUBLE" or type="integer" or type="INTEGER" or ' // &
                         'type="INT" or type=""'
    else if (trim(dataType) /= 'all') then
      call utl_abort( myName//': invalid dataType = ' // trim(dataType) )
    end if
    if (trim(dataType) == 'all') then
      query = 'select name from pragma_table_info("' // trim(tableName) // '");'
    else
      query = 'select name from pragma_table_info("' // trim(tableName) // &
              '") where ' // trim(dataTypeCriteria) // ' ;'
    end if
    call fSQL_prepare( db, trim(query) , stmt, stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myName//': problem with fSQL_get_many '//fSQL_errmsg(stat))
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

    ! clean up and close the sqlite file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine sqlu_getSqlColumnNames

  !--------------------------------------------------------------------------
  ! sqlu_getColumnValuesNum
  !--------------------------------------------------------------------------
  subroutine sqlu_getColumnValuesNum(columnValues, fileName, tableName, &
                                     sqlColumnNames, extraQuery_opt)
    !
    ! :Purpose: Read the column values from sqlite file for the specified table
    !           and column names.
    !
    implicit none

    ! Arguments:
    real(8), allocatable,       intent(out) :: columnValues(:,:)
    character(len=*),           intent(in)  :: sqlColumnNames(:)
    character(len=*),           intent(in)  :: fileName
    character(len=*),           intent(in)  :: tableName
    character(len=*), optional, intent(in)  :: extraQuery_opt

    ! Locals:
    integer :: numRows, numColumns, columnIndex
    character(len=3000)       :: query, extraQuery
    type(fSQL_STATUS)         :: stat ! sqlite error status
    type(fSQL_DATABASE)       :: db   ! sqlite file handle
    type(fSQL_STATEMENT)      :: stmt ! precompiled sqlite statements

    if (present(extraQuery_opt)) then
      extraQuery = trim(extraQuery_opt)
    else
      extraQuery = ''
    end if

    ! open the sqlite file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'sqlu_getColumnValuesNum: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'sqlu_getColumnValuesNum: fSQL_open' )
    end if

    ! build the sqlite query
    query = 'select'
    numColumns = size(sqlColumnNames)
    do columnIndex = 1, numColumns
      query = trim(query) // ' ' // trim(sqlColumnNames(columnIndex))
      if (columnIndex < numColumns) query = trim(query) // ','
    end do
    query = trim(query) // ' from ' // trim(tableName)
    query = trim(query) // ' ' // trim(extraQuery) // ';'
    write(*,*) 'sqlu_getColumnValuesNum: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query) , stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, mode=FSQL_REAL8, &
                        real8_missing=MPC_missingValue_R8, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'sqlu_getColumnValuesNum: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('sqlu_getColumnValuesNum: problem with fSQL_get_many')
    end if
    write(*,*) 'sqlu_getColumnValuesNum: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValues(numRows, numColumns) )
    columnValues(:,:) = 0.0d0
    call fSQL_fill_matrix( stmt, columnValues )

    ! close the sqlite file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine sqlu_getColumnValuesNum

  !--------------------------------------------------------------------------
  ! sqlu_getColumnValuesChar
  !--------------------------------------------------------------------------
  subroutine sqlu_getColumnValuesChar(columnValues, fileName, tableName, &
                                      sqlColumnNames)
    !
    ! :Purpose: Read the column values from sqlite file for the specified table
    !           and column names.
    !
    implicit none

    ! Arguments:
    character(len=50), allocatable, intent(out) :: columnValues(:,:)
    character(len=*),               intent(in)  :: sqlColumnNames(:)
    character(len=*),               intent(in)  :: fileName
    character(len=*),               intent(in)  :: tableName

    ! Locals:
    integer :: numRows, numColumns, columnIndex
    character(len=3000)      :: query
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle
    type(fSQL_STATEMENT)     :: stmt ! precompiled sqlite statements

    ! open the sqlite file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'sqlu_getColumnValuesChar: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'sqlu_getColumnValuesChar: fSQL_open' )
    end if

    ! build the sqlite query
    query = 'select'
    numColumns = size(sqlColumnNames)
    do columnIndex = 1, numColumns
      query = trim(query) // ' ' // trim(sqlColumnNames(columnIndex))
      if (columnIndex < numColumns) query = trim(query) // ','
    end do
    query = trim(query) // ' from ' // trim(tableName) // ';'
    write(*,*) 'sqlu_getColumnValuesChar: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query), stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'sqlu_getColumnValuesChar: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('sqlu_getColumnValuesChar: problem with fSQL_get_many')
    end if
    write(*,*) 'sqlu_getColumnValuesChar: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValues(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnValues )

    ! close the sqlite file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine sqlu_getColumnValuesChar

  !--------------------------------------------------------------------------
  ! sqlu_getColumnValuesDateStr
  !--------------------------------------------------------------------------
  subroutine sqlu_getColumnValuesDateStr(columnDateValues, columnTimeValues, fileName, &
                                         tableName, sqlColumnName)
    !
    ! :Purpose: Read the column values from sqlite file for the specified table
    !           and column names.
    !
    implicit none

    ! Arguments:
    integer, allocatable, intent(out) :: columnDateValues(:)
    integer, allocatable, intent(out) :: columnTimeValues(:)
    character(len=*),     intent(in)  :: sqlColumnName
    character(len=*),     intent(in)  :: fileName
    character(len=*),     intent(in)  :: tableName

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
      write(*,*) 'sqlu_getColumnValuesDateStr: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'sqlu_getColumnValuesDateStr: fSQL_open' )
    end if

    ! Get the date and time

    ! build the sqlite query
    query = "select strftime('%Y%m%d'," // trim(sqlColumnName) // &
            "), strftime('%H%M'," // trim(sqlColumnName) // ") " // &
            "from " // trim(tableName) // ";"
    write(*,*) 'sqlu_getColumnValuesDateStr: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query), stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'sqlu_getColumnValuesDateStr: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('sqlu_getColumnValuesDateStr: problem with fSQL_get_many')
    end if
    write(*,*) 'sqlu_getColumnValuesDateStr: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValuesStr(numRows,2) )
    call fSQL_fill_matrix( stmt, columnValuesStr )
    allocate( columnDateValues(numRows) )
    allocate( columnTimeValues(numRows) )
    do rowIndex = 1, numRows
      read(columnValuesStr(rowIndex,1),*) columnDateValues(rowIndex)
      read(columnValuesStr(rowIndex,2),*) columnTimeValues(rowIndex)
    end do

    deallocate(columnValuesStr)

    ! close the sqlite file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine sqlu_getColumnValuesDateStr

  !--------------------------------------------------------------------------
  ! sqlu_getInitialIdObsData
  !--------------------------------------------------------------------------
  subroutine sqlu_getInitialIdObsData(obsDat, obsFamily, idObs, idData, codeTypeList_opt)
    !
    ! :Purpose: Compute initial value for idObs and idData that will ensure
    !           unique values over all mpi tasks
    !
    implicit none

    ! Arguments:
    type(struct_obs),  intent(inout) :: obsdat
    character(len=*),  intent(in)    :: obsFamily    
    integer         ,  intent(out)   :: idObs
    integer         ,  intent(out)   :: idData
    integer, optional, intent(in)    :: codeTypeList_opt(:)

    ! Locals:
    integer                :: headerIndex, numHeader, numBody, codeType, ierr
    integer, allocatable   :: allNumHeader(:), allNumBody(:)

    numHeader = 0
    numBody = 0
    call obs_set_current_header_list(obsdat, obsFamily)
    HEADERCOUNT: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADERCOUNT
      if (present(codeTypeList_opt)) then
        codeType  = obs_headElem_i(obsdat, OBS_ITY, headerIndex)
        if (all(codeTypeList_opt(:) /= codeType)) cycle HEADERCOUNT
      end if
      numHeader = numHeader + 1
      numBody = numBody + obs_headElem_i(obsdat, OBS_NLV, headerIndex)
    end do HEADERCOUNT
    allocate(allNumHeader(mmpi_nprocs))
    allocate(allNumBody(mmpi_nprocs))
    call rpn_comm_allgather(numHeader,1,'mpi_integer',       &
                            allNumHeader,1,'mpi_integer','GRID',ierr)
    call rpn_comm_allgather(numBody,1,'mpi_integer',       &
                            allNumBody,1,'mpi_integer','GRID',ierr)
    if (mmpi_myid > 0) then
      idObs = sum(allNumHeader(1:mmpi_myid))
      idData = sum(allNumBody(1:mmpi_myid))
    else
      idObs = 0
      idData = 0
    end if
    deallocate(allNumHeader)
    deallocate(allNumBody)

  end subroutine sqlu_getInitialIdObsData

  !--------------------------------------------------------------------------
  ! sqlu_query
  !--------------------------------------------------------------------------
  function sqlu_query(db,query)
    !
    ! :Purpose: To create a query to read an SQLite file
    !
    implicit none

    ! Arguments:
    type(fSQL_DATABASE), intent(inout)  :: db    ! type handle for  SQLIte file
    character(len=*),    intent(in)     :: query
    ! Result:
    character(len=256)   :: sqlu_query

    ! Locals:
    character(len=256)   :: result
    logical :: finished
    type(fSQL_STATEMENT) :: stmt !  prepared statement for  SQLite
    type(fSQL_STATUS)    :: stat !type error status

    result=''
    call fSQL_prepare(db, trim(query), stmt, stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat,'fSQL_prepare: ')
    finished=.false.
    call fSQL_get_row(stmt, finished)

    ! Put result of query into variable
    call fSQL_get_column(stmt, COL_INDEX = 1, char_var = result)
    call fSQL_get_row(stmt, finished)
    if (.not. finished) write(*,*) ' SQL QUERY ---> QUERY RETURNS MORE THAN ONE ROW...  '
    call fSQL_finalize(stmt)
    sqlu_query = trim(result)

  end function sqlu_query

  !--------------------------------------------------------------------------
  ! sqlu_handleError
  !--------------------------------------------------------------------------
  subroutine sqlu_handleError(stat, message)
    implicit none

    ! Arguments:
    type(FSQL_STATUS), intent(in) :: stat
    character(len=*),  intent(in) :: message

    write(*,*) message, fSQL_errmsg(stat)
    call utl_abort(trim(message))

  end subroutine sqlu_handleError

end module sqliteUtilities_mod

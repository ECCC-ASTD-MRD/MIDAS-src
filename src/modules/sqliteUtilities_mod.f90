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

module sqliteUtilities_mod
  ! MODULE sqliteUtilities (prefix='sqlu' category='6. Observation input/output')
  !
  ! :Purpose: A place to collect utilities for SQLite files.

use fSQLite
use clib_interfaces_mod
use utilities_mod

implicit none
save
private
public :: sqlu_sqlColumnExists, sqlu_sqlTableExists, sqlu_getSqlColumnNames

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

    ! arguments:
    character(len=*),              intent(in)  :: fileName
    character(len=*),              intent(in)  :: tableName
    character(len=*),              intent(in)  :: columnName
    logical                                    :: columnExists

    ! locals:
    integer                     :: ierr
    logical                     :: finished
    character(len=3000)         :: query, sqliteOutput
    character(len=lenSqlName)   :: upperColumnName
    type(fSQL_STATUS)           :: stat ! sqlite error status
    type(fSQL_DATABASE)         :: db   ! sqlite file handle
    type(fSQL_STATEMENT)        :: stmt ! precompiled sqlite statements
    logical, parameter          :: debug = .false.
    character(len=*), parameter :: myName = 'sqlu_sqlColumnExists'
    
    ! open the SQLite file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myName//fSQL_errmsg(stat) )
    end if

    upperColumnName = trim(columnName)
    ierr = clib_toUpper(upperColumnName)

    query = "select upper(name) as uppername from pragma_table_info('" // &
            trim(tableName) // "') where uppername='" // trim(upperColumnName) // "';"
    if (debug) write(*,*) myName//': query = ', trim(query)

    call fSQL_prepare( db, trim(query), stmt, stat)
    finished = .false.
    call fSQL_get_row( stmt, finished )
    call fSQL_get_column( stmt, COL_INDEX = 1, CHAR_VAR = sqliteOutput )
    call fSQL_get_row( stmt, finished )
    call fSQL_finalize( stmt )
    if (debug) write(*,*) myName//': output = XXX' // trim(sqliteOutput) // 'XXX'
    columnExists = (trim(sqliteOutput) == trim(upperColumnName))

    ! close the obsDB file
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

    ! arguments:
    character(len=*),              intent(in)  :: fileName
    character(len=*),              intent(in)  :: tableName
    logical                                    :: tableExists

    ! locals:
    integer                     :: ierr
    logical                     :: finished
    character(len=3000)         :: query, sqliteOutput
    character(len=lenSqlName)   :: upperTableName
    type(fSQL_STATUS)           :: stat ! sqlite error status
    type(fSQL_DATABASE)         :: db   ! sqlite file handle
    type(fSQL_STATEMENT)        :: stmt ! precompiled sqlite statements
    logical, parameter          :: debug = .false.
    character(len=*), parameter :: myName = 'sqlu_sqlTableExists'

    ! open the obsDB file
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

    ! close the obsDB file
    call fSQL_close( db, stat ) 

  end function sqlu_sqlTableExists

  !--------------------------------------------------------------------------
  ! sqlu_getSqlColumnNames
  !--------------------------------------------------------------------------
  subroutine sqlu_getSqlColumnNames(sqlColumnNames, fileName, tableName, dataType)
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
    character(len=*), parameter :: myName = 'sqlu_getSqlColumnNames'

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myName//': fSQL_open '//fSQL_errmsg(stat) )
    end if

    ! read the column names
    if (trim(dataType) == 'varchar') then
      dataTypeCriteria = 'substr(type,1,7)="varchar"'
    else if (trim(dataType) == 'numeric') then
      dataTypeCriteria = 'type="real" or type="REAL" or type="double" or ' // &
                         'type="DOUBLE" or type="integer" or type="INTEGER"'
    else
      call utl_abort( myName//': invalid dataType = ' // trim(dataType) )
    end if
    query = 'select name from pragma_table_info("' // trim(tableName) // &
            '") where ' // trim(dataTypeCriteria) // ' ;'
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

    ! clean up and close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine sqlu_getSqlColumnNames
  
end module sqliteUtilities_mod
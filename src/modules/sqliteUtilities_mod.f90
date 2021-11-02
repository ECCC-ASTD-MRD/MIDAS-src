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
public :: sqlu_sqlColumnExists

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
    character(len=*), parameter :: myError   = '******** '// myName //' ERROR: '
    
    ! open the SQLite file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myError//fSQL_errmsg(stat) )
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
  
end module sqliteUtilities_mod

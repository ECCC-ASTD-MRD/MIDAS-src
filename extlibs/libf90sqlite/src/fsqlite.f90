! fsqlite.f90 --
!    Module for interacting with SQLite3
!    Hamid.Benhocine
!
!    General information:
!    inspired by work by Al Danial and Arjen Marcus
!    (http://danial.org).
!    The purpose is to provide a high-level means
!    for Fortran programmers to use the SQLite
!    library by Richard Hipp (http://www.sqlite.org)
!
!    Detailed documentation can be found in
!    the file fSQL.html.
!    This work is inspired by the alda
!
!    TODO:
!    - Think about finalize and reset: error code?
!    - Support NULL (ok!)
!    - Increase Error checking
!    - fSQL_do_script
!    - fSQL_last_rowid, fSQLcol_names, fSQL_col_types
!
!    Implementation notes:
!    The handles to the database or prepared statements
!    are stored in an array of two integers to take
!    care of 64-bits platforms.
!    With the appropriate compilation options (if needed)
!    the code should be thread-safe, as no data are shared.
!

module fSQL_types

   !! constantes {{{1
   integer, parameter            :: LEN_ERR_MSG        = 256 
   integer, parameter            :: LEN_COL_NAME       = 64 
   integer, parameter            :: LEN_COL_TYPE       = 64 
   integer, parameter            :: LEN_MISSING        = 256
   character(len = 4), parameter :: CHAR_NULL_MISSING  = "NULL"
   integer, parameter            :: FSQL_dp            = kind(1.0d00)
   integer, parameter            :: FSQL_INT_KIND_CONST= SELECTED_INT_KIND(8)
   logical, parameter            :: DBG                = .false.

   integer, parameter            :: FSQL_OK            = 0

   integer, parameter            :: SQLITE_INT         = 1
   integer, parameter            :: SQLITE_INT4        = 1
   integer, parameter            :: FSQL_INT           = 1

   integer, parameter            :: SQLITE_INT8        = 2
   integer, parameter            :: FSQL_INT8          = 2

   integer, parameter            :: SQLITE_REAL        = 3
   integer, parameter            :: SQLITE_REAL4       = 3
   integer, parameter            :: FSQL_REAL          = 3

   integer, parameter            :: SQLITE_REAL8       = 4
   integer, parameter            :: FSQL_REAL8         = 4

   integer, parameter            :: SQLITE_CHAR        = 5
   integer, parameter            :: FSQL_CHAR          = 5

   integer, parameter            :: SQLITE_OK          = 0
   integer, parameter            :: SQLITE_ERROR       = 1
   integer, parameter            :: SQLITE_MISUSE      = 21
   integer, parameter            :: SQLITE_ROW         = 100
   integer, parameter            :: SQLITE_DONE        = 101

   CHARACTER, PARAMETER:: INITIAL='*'
   INTEGER             :: ID_LAST=0
   !! fin constantes 1}}}

! FSQL DATATYPES
!
   TYPE T_INIT
            INTEGER   :: ID
            CHARACTER :: DATA
   END TYPE
   type fSQL_STATEMENT ! {{{1
      integer, dimension(2)   :: stmt_handle
      integer                 :: column_count
      integer                 :: row_count
      integer, dimension(2)   :: ptr_handle
      integer                 :: ptr_type 
      character(len=LEN_COL_NAME), dimension(:), pointer :: column_name =>NULL()
      character(len=LEN_COL_TYPE), dimension(:), pointer :: column_type =>NULL()
      type(T_INIT)            :: INIT
   end type! 1}}}
   type fSQL_DATABASE !{{{1
      integer, dimension(2)   :: db_handle
      integer                 :: error
      character(len=LEN_ERR_MSG)      :: errmsg
   end type ! 1}}}
   type fSQL_STATUS !{{{1
      integer                         :: error
      character(len=LEN_ERR_MSG)      :: errmsg
   end type ! 1}}}
   type SQLITE_COLUMN !{{{1
      character(len=40)       :: name
      character(len=40)       :: type
      character(len=40)       :: function
      integer                 :: type_set
      integer                 :: int_value
      integer*8               :: int64_value
      real*4                  :: float_value
      real(kind=FSQL_dp)      :: double_value
      character(len=80)       :: char_value
   end type ! 1}}}



   TYPE(T_INIT), PARAMETER, PUBLIC:: INIT_0 = T_INIT(0,INITIAL)
   PUBLIC :: T_INIT, IS_INIT, INITIALIZE, TEST_INIT, GET_ID,GET_DATA

   CONTAINS

   LOGICAL FUNCTION IS_INIT(INIT)
            TYPE(T_INIT), INTENT(IN) :: INIT
            IS_INIT = (INIT%DATA==INITIAL).AND.INIT%ID>0.AND.INIT%ID<=ID_LAST
   END FUNCTION

   SUBROUTINE INITIALIZE (INIT)
            TYPE(T_INIT), INTENT(INOUT) :: INIT
            ID_LAST = ID_LAST +1
            INIT = T_INIT(ID_LAST,INITIAL)
   END SUBROUTINE

   SUBROUTINE TEST_INIT (INIT)
            TYPE(T_INIT), INTENT(INOUT) :: INIT
            IF(.NOT.IS_INIT(INIT)) CALL INITIALIZE(INIT)
   END SUBROUTINE

   INTEGER FUNCTION GET_ID(INIT)
            TYPE(T_INIT), INTENT(IN) :: INIT
            GET_ID = INIT%ID
   END FUNCTION

   CHARACTER FUNCTION GET_DATA(INIT)
            TYPE(T_INIT), INTENT(IN) :: INIT
            GET_DATA= INIT%DATA
   END FUNCTION

end module fSQL_types

module fSQLite
   use fSQL_types

   implicit none

   private :: stringtof
   private :: stringtoc
   private :: typename

   !
   ! Convenient interfaces
   !
   interface fSQL_bind_param !{{{1
      module procedure fSQL_bind_param_int
      module procedure fSQL_bind_param_int8
      module procedure fSQL_bind_param_real
      module procedure fSQL_bind_param_real8
      module procedure fSQL_bind_param_char
      module procedure fSQL_bind_param_null
   end interface
   !1}}}
   interface fSQL_get_column !{{{1
      module procedure fSQL_get_column_int
      module procedure fSQL_get_column_int_miss
      module procedure fSQL_get_column_int8
      module procedure fSQL_get_column_int8_miss
      module procedure fSQL_get_column_real
      module procedure fSQL_get_column_real_miss
      module procedure fSQL_get_column_real8
      module procedure fSQL_get_column_real8_miss
      module procedure fSQL_get_column_char
      module procedure fSQL_get_column_char_miss
   end interface
   !1}}}
   interface fSQL_fill_matrix !{{{1
      module procedure fSQL_fill_matrix_int
      module procedure fSQL_fill_matrix_int8
      module procedure fSQL_fill_matrix_real
      module procedure fSQL_fill_matrix_real8
      module procedure fSQL_fill_matrix_char
   end interface
   interface fSQL_write_matrix !{{{1
      module procedure fSQL_write_matrix_int
      module procedure fSQL_write_matrix_int8
      module procedure fSQL_write_matrix_real
      module procedure fSQL_write_matrix_real8
      module procedure fSQL_write_matrix_char
   end interface
   !1}}}

   ! interfaces aux fonctions du module c 
   ! toutes les fonctions du module c
   ! ne sont pas necessairement interfaces
   ! 
   interface
      integer function sqlite3_bind_int_c( handle, colidx, value ) ! {{{1
         integer, dimension(*) :: handle
         integer               :: colidx
         integer               :: value
      end function sqlite3_bind_int_c ! 1}}}
   end interface

   interface
      integer function sqlite3_bind_int64_c( handle, colidx, value ) ! {{{1
         integer, dimension(*) :: handle
         integer               :: colidx
         integer*8             :: value
      end function sqlite3_bind_int64_c ! 1}}}
   end interface

   interface
      integer function sqlite3_bind_float_c( handle, colidx, value ) ! {{{1
         use fSQL_types
         integer, dimension(*) :: handle
         integer               :: colidx
         real*4                :: value
      end function sqlite3_bind_float_c ! 1}}}
   end interface

   interface
      integer function sqlite3_bind_double_c( handle, colidx, value ) ! {{{1
         use fSQL_types
         integer, dimension(*) :: handle
         integer               :: colidx
         real(kind=FSQL_dp)         :: value
      end function sqlite3_bind_double_c ! 1}}}
   end interface

   interface
      integer function sqlite3_bind_text_c( handle, colidx, value ) ! {{{1
         integer, dimension(*) :: handle
         integer               :: colidx
         character(len=*)      :: value
      end function sqlite3_bind_text_c ! 1}}}
   end interface

   interface
      integer function sqlite3_bind_null_c( handle, colidx) ! {{{1
         integer, dimension(*) :: handle
         integer               :: colidx
      end function ! 1}}}
   end interface

   interface
      integer function sqlite3_column_int_c( handle, colidx, value ) ! {{{1
         integer, dimension(*) :: handle
         integer               :: colidx
         integer               :: value
      end function sqlite3_column_int_c ! 1}}}
   end interface


   ! interface pour sqlite3_column_int_c avec gestion de
   ! la cellule NULL
   interface
      integer function sqlite3_column_intm_c( handle, colidx, value,missing ) ! {{{1
         integer, dimension(*) :: handle
         integer               :: colidx
         integer               :: value
         integer               :: missing
      end function sqlite3_column_intm_c ! 1}}}
   end interface

   interface
      integer function sqlite3_column_int64_c( handle, colidx, value ) ! {{{1
         integer, dimension(*) :: handle
         integer               :: colidx
         integer*8             :: value
      end function sqlite3_column_int64_c ! 1}}}
   end interface

   ! interface pour sqlite3_column_int64_c avec gestion de
   ! la cellule NULL
   interface
      integer function sqlite3_column_int64m_c( handle, colidx, value,missing ) ! {{{1
         integer, dimension(*) :: handle
         integer               :: colidx
         integer*8             :: value
         integer*8             :: missing
      end function sqlite3_column_int64m_c ! 1}}}
   end interface 

   interface
      integer function sqlite3_column_float_c( handle, colidx, value ) ! {{{1
         use fSQL_types
         integer, dimension(*) :: handle
         integer               :: colidx
         real*4                :: value
      end function ! 1}}}}
   end interface 

   interface
   ! interface pour sqlite3_column_float_c avec gestion de
   ! la cellule NULL
      integer function sqlite3_column_floatm_c( handle, colidx, value ,missing) ! {{{1
         use fSQL_types
         integer, dimension(*) :: handle
         integer               :: colidx
         real*4                :: value
         real*4                :: missing
      end function ! 1}}}}
   end interface 

   interface
      integer function sqlite3_column_double_c( handle, colidx, value ) ! {{{1
         use fSQL_types
         integer, dimension(*) :: handle
         integer               :: colidx
         real(kind=FSQL_dp)         :: value
      end function sqlite3_column_double_c ! 1}}}
   end interface 

   ! interface pour sqlite3_column_double_c avec gestion de
   ! la cellule NULL
   interface
      integer function sqlite3_column_doublem_c( handle, colidx, value,missing ) ! {{{1
         use fSQL_types
         integer, dimension(*) :: handle
         integer               :: colidx
         real(kind=FSQL_dp)         :: value
         real(kind=FSQL_dp)         :: missing
      end function ! 1}}} 
   end interface 

   interface
      integer function sqlite3_column_text_c( handle, colidx, value ) ! {{{1
         integer, dimension(*) :: handle
         integer               :: colidx
         character(len=*)      :: value
      end function ! 1}}}
   end interface 

   ! interface pour sqlite3_column_text_c avec gestion de
   ! la cellule NULL
   interface
      integer function sqlite3_column_textm_c( handle, colidx, value,miss ) ! {{{1
         integer, dimension(*) :: handle
         integer               :: colidx
         character(len=*)      :: value
         character(len=*)      :: miss
      end function  ! 1}}}
   end interface
   interface
      integer function sqlite3_open_c( fnamec, handle ) ! {{{1
         character(len=*)      :: fnamec
         integer, dimension(*) :: handle
      end function sqlite3_open_c ! 1}}}
   end interface
   interface
      integer function sqlite3_close_c( handle ) ! {{{1
         integer, dimension(*) :: handle
      end function sqlite3_close_c ! 1}}}
   end interface
   interface
      integer function sqlite3_prepare_c( db, command, stmt ) ! {{{1
         integer, dimension(*) :: db
         character(len=*)      :: command
         integer, dimension(*) :: stmt
      end function sqlite3_prepare_c  ! 1}}}
   end interface
   interface
      subroutine sqlite3_column_count_c( handle, count ) ! {{{1
         integer, dimension(*) :: handle
         integer               :: count
      end subroutine sqlite3_column_count_c ! 1}}}
   end interface
   interface
      subroutine sqlite3_column_name_type_c( handle, colidx, name, type ) !{{{1
         integer, dimension(*) :: handle
         integer               :: colidx
         character(len=*)      :: name
         character(len=*)      :: type
      end subroutine sqlite3_column_name_type_c ! 1}}}
   end interface
   interface
      subroutine sqlite3_step_c( stmt, completion ) ! {{{1
         integer, dimension(*) :: stmt
         integer               :: completion
      end subroutine sqlite3_step_c ! 1}}}
   end interface

contains


! character(len=40) function typename( column, primary ) {{{1
!    Construct the type and attributes of a column
!    in a new table
! Arguments:
!    column        Column information
!    primary       Name of the primary key
!
character(len=40) function typename( column, primary )
   type(SQLITE_COLUMN), intent(in) :: column
   character(len=*), intent(in)    :: primary

   if ( column%name .ne. primary ) then
      typename = column%type
   else
!      write( typename, '(2a)' ) trim(column%type), ' primary key'
       typename = trim(column%type) // ' primary key'
   endif

end function ! 1}}}

! subroutine stringtof( string ) {{{1
!    Convert a C string to Fortran
! Arguments:
!    string        String to be converted
!
subroutine stringtof( string )
   character(len=*) :: string

   integer          :: last
   last = index( string, char(0) )
   if ( last .gt. 0 ) then
      string(last:) = ' '
   endif

end subroutine ! 1}}}
! subroutine stringtoc( string ) {{{1
!    Convert a Fortran string to C
! Arguments:
!    string        String to be converted
! Note:
!    It is assumed that the last character
!    is a space. As this is a private
!    routine, this should have been taken
!    care of in the caller.
!
subroutine stringtoc( string )
   character(len=*) :: string

   integer          :: last

   last = 1 + len_trim(string)
   string(last:last) = char(0)

end subroutine ! 1}}}

! subroutine sqlite3_column_props( column, name, type, length ) ! {{{1
!    Convenience routine to set the properties of a column
! Arguments:
!    column        Column structure
!    name          Name of the column
!    type          Type of the column
!    length        Length if a string
! Side effects:
!    Fields in column filled
!
subroutine sqlite3_column_props( column, name, type, length )
   type(SQLITE_COLUMN), intent(inout) :: column
   character(len=*), intent(in)       :: name
   integer, intent(in)                :: type
   integer, intent(in), optional      :: length

   integer                            :: length_
   character(len=40)                  :: type_expr

   length_ = 20
   if ( present(length) ) then
      length_ = length
   endif

   column%name     = name
   column%type_set = type

   select case ( type )
   case (SQLITE_INT4)
      column%type = 'INT'
   case (SQLITE_INT8)
      ! column%type = 'INT64'
      column%type = 'INT'
   case (SQLITE_REAL4)
      column%type = 'FLOAT'
   case (SQLITE_REAL8)
      column%type = 'DOUBLE'
   case (SQLITE_CHAR)
      write( column%type, '(a,i0,a)' ) 'CHAR(', length_, ')'
   case default
      column%type = 'UNKNOWN!'
   end select

end subroutine  ! 1}}}

! fSQL_bind_param_XXXX Definbitions  
! fSQL_bind_param_int   -- 
! fSQL_bind_param_int8 -- 
! fSQL_bind_param_real  -- 
! fSQL_bind_param_real8 -- 
! fSQL_bind_param_char  -- 
! fSQL_bind_param_null  -- 
!
!    Convenience routines to set the value of a column
! Arguments:
!    column        Column structure
!    value         The value to be set
! Side effects:
!    Appropriate value field in column set
!
!-----------------------------------------------------------------!
subroutine fSQL_bind_param_int( stmt, param_index, INT_VAR ) ! {{{1
   type(fSQL_STATEMENT), intent(inout) :: stmt
   integer, intent(in)                   :: param_index
   integer, intent(in)                   :: INT_VAR

   integer                               :: rc
   rc = sqlite3_bind_int_c( stmt%stmt_handle, param_index, INT_VAR )
end subroutine fSQL_bind_param_int ! 1}}}

subroutine fSQL_bind_param_int8( stmt, param_index, INT8_VAR ) ! {{{1
   type(fSQL_STATEMENT), intent(inout) :: stmt
   integer, intent(in)                   :: param_index
   integer*8, intent(in)                 :: INT8_VAR

   integer                               :: rc
   rc = sqlite3_bind_int64_c( stmt%stmt_handle, param_index, INT8_VAR )
end subroutine fSQL_bind_param_int8 ! 1}}}

subroutine fSQL_bind_param_real(stmt, param_index, REAL_VAR ) ! {{{1
   type(fSQL_STATEMENT), intent(inout) :: stmt
   integer, intent(in)                   :: param_index
   real*4, intent(in)                    :: REAL_VAR

   integer                               :: rc
   rc = sqlite3_bind_float_c( stmt%stmt_handle, param_index, REAL_VAR )
end subroutine  ! 1}}}

subroutine fSQL_bind_param_real8(stmt, param_index, REAL8_VAR ) ! {{{1
   type(fSQL_STATEMENT), intent(inout) :: stmt
   integer, intent(in)                   :: param_index
   real(kind=FSQL_dp), intent(in)             :: REAL8_VAR

   integer                               :: rc
   rc = sqlite3_bind_double_c( stmt%stmt_handle, param_index, REAL8_VAR )
end subroutine ! 1}}} 

subroutine fSQL_bind_param_char( stmt, param_index, CHAR_VAR ) ! {{{1
   type(fSQL_STATEMENT), intent(inout) :: stmt
   integer, intent(in)                   :: param_index
   character(len=*), intent(in)          :: CHAR_VAR

   integer                               :: rc
   rc = sqlite3_bind_text_c( stmt%stmt_handle, param_index, trim(CHAR_VAR) )
end subroutine ! 1}}} 

subroutine fSQL_bind_param_null( stmt, param_index ) ! {{{1
   type(fSQL_STATEMENT), intent(inout)   :: stmt
   integer, intent(in)                   :: param_index

   integer                               :: rc
   rc = sqlite3_bind_null_c ( stmt%stmt_handle, param_index)
end subroutine ! 1}}} 

! FSQL get_column_XXXX  Definitions 
! fSQL_get_column_int        -- 
! fSQL_get_column_int_miss   -- 
! fSQL_get_column_int8       -- 
! fSQL_get_column_int8_miss  -- 
! fSQL_get_column_real       -- 
! fSQL_get_column_real_miss  -- 
! fSQL_get_column_real8      -- 
! fSQL_get_column_real8_miss -- 
! fSQL_get_column_char       -- 
! fSQL_get_column_char_miss  -- 
!
!    Convenience routines to get the value of a column
! Arguments:
!    col_index     Column index of the result (zero indexed in C)
!    value         Value on return
! Side effects:
!    Value argument will be set
! Note:
!
!-----------------------------------------------------------------!
subroutine fSQL_get_column_int( stmt, COL_INDEX, INT_VAR ) ! {{{1
   type(fSQL_STATEMENT), intent(in)   :: stmt
   integer, intent(in)                :: COL_INDEX
   integer, intent(out)               :: INT_VAR
   integer                            :: rc

   rc = sqlite3_column_int_c( stmt%stmt_handle, (COL_INDEX - 1), INT_VAR )
   if (DBG) write(*,*) 'column int_var called'
end subroutine  ! 1}}}

subroutine fSQL_get_column_int_miss( stmt, COL_INDEX, INT_VAR ,INT_MISSING) ! {{{1
   type(fSQL_STATEMENT), intent(in)   :: stmt
   integer, intent(in)                :: COL_INDEX
   integer, intent(out)               :: INT_VAR
   integer, intent(in)                :: INT_MISSING
   integer                            :: rc
   integer                            :: missingi
   missingi = int_missing

   rc = sqlite3_column_intm_c( stmt%stmt_handle, (COL_INDEX - 1), INT_VAR ,missingi)
   if (DBG) write(*,*) 'column int_var missing called'
end subroutine  ! 1}}}

subroutine fSQL_get_column_int8( stmt, COL_INDEX, INT8_VAR ) ! {{{1
   type(fSQL_STATEMENT), intent(in) :: stmt
   integer, intent(in)                :: COL_INDEX
   integer*8, intent(out)             :: INT8_VAR
   integer                            :: rc

   rc = sqlite3_column_int64_c( stmt%stmt_handle, (COL_INDEX -1), INT8_VAR )
   if (DBG) write(*,*) 'int8_var called'
end subroutine ! 1}}} 

subroutine fSQL_get_column_int8_miss( stmt, COL_INDEX, INT8_VAR,INT8_MISSING ) ! {{{1
   type(fSQL_STATEMENT), intent(in) :: stmt
   integer, intent(in)                :: COL_INDEX
   integer*8, intent(out)             :: INT8_VAR
   integer*8, intent(in)              :: INT8_MISSING
   integer*8                          :: missing64
   integer                            :: rc
   missing64 = int8_missing

   rc = sqlite3_column_int64m_c( stmt%stmt_handle, (COL_INDEX -1), INT8_VAR,missing64 )
   if (DBG) write(*,*) 'int8_var missing called'
end subroutine ! 1}}} 

subroutine fSQL_get_column_real( stmt, COL_INDEX, REAL_VAR ) ! {{{1
   type(fSQL_STATEMENT), intent(in) :: stmt
   integer, intent(in)                :: COL_INDEX
   real, intent(out)                  :: REAL_VAR
   integer                            :: rc

   rc = sqlite3_column_float_c( stmt%stmt_handle, (COL_INDEX -1), REAL_VAR )
   if (DBG) write(*,*) 'real_var called'
end subroutine ! 1}}} 

subroutine fSQL_get_column_real_miss( stmt, COL_INDEX, REAL_VAR,REAL_MISSING ) ! {{{1
   type(fSQL_STATEMENT), intent(in) :: stmt
   integer, intent(in)                :: COL_INDEX
   real, intent(out)                  :: REAL_VAR
   real, intent(in)                   :: REAL_MISSING
   integer                            :: rc

   rc = sqlite3_column_floatm_c( stmt%stmt_handle, (COL_INDEX -1), REAL_VAR, REAL_MISSING )
   if (DBG) write(*,*) 'real_var missing called'
end subroutine ! 1}}} 

subroutine fSQL_get_column_real8( stmt, COL_INDEX, REAL8_VAR ) !{{{1
   type(fSQL_STATEMENT), intent(in) :: stmt
   integer, intent(in)                :: COL_INDEX
   real(kind=FSQL_dp), intent(out)         :: REAL8_VAR 
   integer                            :: rc

   rc = sqlite3_column_double_c( stmt%stmt_handle, (COL_INDEX -1), REAL8_VAR )
   if (DBG) write(*,*) 'real8_var called'
end subroutine ! 1}}} 

subroutine fSQL_get_column_real8_miss( stmt, COL_INDEX, REAL8_VAR, REAL8_MISSING ) ! {{{1
   type(fSQL_STATEMENT), intent(in) :: stmt
   integer, intent(in)                :: COL_INDEX
   real(kind=FSQL_dp), intent(out)         :: REAL8_VAR 
   real(kind=FSQL_dp), intent(in)          :: REAL8_MISSING
   integer                            :: rc

   rc = sqlite3_column_doublem_c( stmt%stmt_handle, (COL_INDEX -1), REAL8_VAR,REAL8_MISSING )
   if (DBG) write(*,*) 'real8_var  missing  called'
end subroutine  ! 1}}}

subroutine fSQL_get_column_char( stmt, COL_INDEX, CHAR_VAR ) ! {{{1
   type(fSQL_STATEMENT), intent(in) :: stmt
   integer, intent(in)                :: COL_INDEX
   character(len=*), intent(out)      :: CHAR_VAR 
   integer                            :: rc

   rc = sqlite3_column_text_c( stmt%stmt_handle, (COL_INDEX -1), CHAR_VAR )
   call stringtof( CHAR_VAR )
   if (DBG) write(*,*) 'char_var called'
end subroutine  ! 1}}}

subroutine fSQL_get_column_char_miss( stmt, COL_INDEX, CHAR_VAR, CHAR_MISSING ) ! {{{1
   type(fSQL_STATEMENT), intent(in) :: stmt
   integer, intent(in)                :: COL_INDEX
   character(len=*), intent(out)      :: CHAR_VAR 
   character(len=*), intent(in)       :: CHAR_MISSING
   integer                            :: rc
   character(len=LEN_MISSING)         :: missingc

   missingc = char_missing
   call stringtoc( missingc )

   rc = sqlite3_column_textm_c( stmt%stmt_handle, (COL_INDEX -1), CHAR_VAR,missingc )
   call stringtof( CHAR_VAR )
   if (DBG) write(*,*) 'char_var  missing called'
end subroutine ! 1}}}

! FSQL error managment
!
! integer function fSQL_error (status) {{{1
!    Return the last error code
! Arguments:
!    db            Structure for the database
! Returns:
!    Last fSQL error code for this library
!
!-----------------------------------------------------------------!
integer function fSQL_error( status )
   type(fSQL_STATUS) :: status

   fSQL_error = status%error
end function  ! 1}}}

! character(len=LEN_ERR_MSG) function fSQL_errmsg( status ) {{{1
!    Return the last error message
! Arguments:
!    status            Structure for the error status
! Returns:
!    Last fSQL error message for this library
!
!-----------------------------------------------------------------!
character(len=LEN_ERR_MSG) function fSQL_errmsg( status )
   type(fSQL_STATUS) :: status

   fSQL_errmsg = status%errmsg
end function ! 1}}}

! integer function fSQL_last_insert_rowid (db) {{{1
!    Return the last inserted key for AUTOINCREMENT PRIMARY integer key
! Arguments:
!    db            Structure for the database
! Returns:
!    Return the last inserted key for AUTOINCREMENT PRIMARY integer key
!
!-----------------------------------------------------------------!
integer function fSQL_last_insert_rowid( db )
   type(fSQL_DATABASE) :: db
   integer*8           :: rowid

   interface
      integer*8 function  sqlite3_last_insert_rowid_c( handle ) 
         integer, dimension(*) :: handle
      end function 
   end interface
   rowid = sqlite3_last_insert_rowid_c(db%db_handle)
   FSQL_last_insert_rowid = INT(rowid,FSQL_INT_KIND_CONST)
end function  ! 1}}}

! FSQL (open, close ) procedures
!
! subroutine fSQL_open(db, fname, status) {{{1
!    Open a database file
! Arguments:
!    fname         Name of the file
!    db            Structure for the database
!    status        error status
! Side effects:
!    The database file is opened and can be
!    used via the db argument
!
subroutine fSQL_open( db, fname, status )
   character(len=*)            :: fname
   type(fSQL_DATABASE)         :: db
   type(fSQL_STATUS), optional :: status

   character(len=len(fname)+1) :: fnamec
   logical                     :: file_exist

   db%db_handle   = 0
   db%error       = 0
   db%errmsg       = ' '
   
   !
   ! check existence of the file
!   inquire(file=fname(1:),exist=file_exist)
!   if (.not.file_exist) then
!      db%error       = -1 
!      db%errmsg = "fsql_open >>"//fname(1:len_trim(fname))// &
!      & " no such file check the path - will be created if mode create !!!"
!      return
!   endif
   fnamec = fname
   call stringtoc( fnamec )

   db%error = sqlite3_open_c( fnamec, db%db_handle )
   if ( db%error /= SQLITE_OK) then
      call sqlite3_errmsg_c( db%db_handle, db%errmsg )
   end if
   if (present(status) ) then
      status%error      = db%error
      status%errmsg     = db%errmsg
   end if
end subroutine ! 1}}} 

! subroutine fSQL_close( db, status)  {{{1 
!    Close a database file
! Arguments:
!    db            Structure for the database
!    status        error status
! Side effects:
!    The database file is closed and can no
!    longer be accessed
!
!-----------------------------------------------------------------!
subroutine fSQL_close( db, status )
   type(fSQL_DATABASE)         :: db
   type(fSQL_STATUS), optional :: status

   db%error       = 0
   db%errmsg       = ' '

   db%error = sqlite3_close_c( db%db_handle )
   if (db%error /= SQLITE_OK) then
      call sqlite3_errmsg_c( db%db_handle, db%errmsg )
   end if

   if (present(status) ) then
      status%error      = db%error
      status%errmsg     = db%errmsg
   end if
   db%db_handle   = 0

end subroutine fSQL_close ! 1}}}

! FSQL DML exec procedures
!
! subroutine fSQL_do ( db, command, status) !{{{1
!    Run a single SQL command
! Arguments:
!    db            Structure for the database
!    command       Complete SQL command
!    status        error status
! Side effects:
!
!-----------------------------------------------------------------!
subroutine fSQL_do( db, command, status )
   type(fSQL_DATABASE)         :: db
   character(len=*)            :: command
   type(fSQL_STATUS), optional :: status

   interface
      integer function sqlite3_do_c( handle, command, errmsg )
         integer, dimension(*) :: handle
         character(len=*)      :: command
         character(len=*)      :: errmsg
      end function sqlite3_do_c
   end interface

   character(len=len(command)+1) :: commandc
   integer                       :: k

   db%error       = 0
   db%errmsg       = ' '

   commandc = command
   call stringtoc( commandc )

   db%error  = sqlite3_do_c( db%db_handle, commandc, db%errmsg )

   if (present(status) ) then
      status%error      = db%error
      status%errmsg     = db%errmsg
   end if
end subroutine ! 1}}}

! subroutine fSQL_do_many(db, command, status) ! {{{1
!    Run a single many SQL commands
! Arguments:
!    db            Structure for the database
!    command       Complete SQL command
!    status        error status
! Side effects:
!
!-----------------------------------------------------------------!
subroutine fSQL_do_many( db, command, status )
   type(fSQL_DATABASE)         :: db
   character(len=*)            :: command
   type(fSQL_STATUS), optional :: status

   interface
      integer function sqlite3_do_many_c( handle, command, errmsg )
         integer, dimension(*) :: handle
         character(len=*)      :: command
         character(len=*)      :: errmsg
      end function sqlite3_do_many_c
   end interface

   character(len=len(command)+1) :: commandc
   integer                       :: k

   db%error       = 0
   db%errmsg       = ' '

   commandc = command
   call stringtoc( commandc )

   db%error  = sqlite3_do_many_c( db%db_handle, commandc, db%errmsg )
   if (present(status) ) then
      status%error      = db%error
      status%errmsg     = db%errmsg
   end if

end subroutine ! 1}}}

! subroutine fSQL_do_script(db, fname, status) ! {{{1
!    Run a single many SQL commands
! Arguments:
!    db            Structure for the database
!    fname         file containing SQL script
!    status        error status
! Side effects:
!
!-----------------------------------------------------------------!
subroutine fSQL_do_script( db,fname, status )
   type(fSQL_DATABASE)         :: db
   character(len=*)            :: fname
   type(fSQL_STATUS), optional :: status

   interface
      integer function sqlite3_do_script_c( handle, fname, errmsg )
         integer, dimension(*) :: handle
         character(len=*)      :: fname
         character(len=*)      :: errmsg
      end function 
   end interface

   character(len=len(fname)+1) :: commandc
   integer                     :: k

   db%error       = 0
   db%errmsg       = ' '

   commandc = fname 
   call stringtoc( commandc )

   db%error  = sqlite3_do_script_c( db%db_handle, commandc, db%errmsg )
   if (present(status) ) then
      status%error      = db%error
      status%errmsg     = db%errmsg
   end if


end subroutine ! 1}}}


! FSQL_Transactions procedures
!
! subroutine fSQL_begin (db, status) {{{1
!    Start a transaction on the given database
! Arguments:
!    db            Structure for the database
!    status        error status
! Note:
!    Should be accompanied by a call to either
!    fSQL_commit or fSQL_rollback
!
!-----------------------------------------------------------------!
subroutine fSQL_begin( db, status )
   type(fSQL_DATABASE)         :: db
   type(fSQL_STATUS), optional :: status

   if (present(status)) then
      status%error = 0
      status%errmsg = ' '
   end if

   if (present(status)) then
      call fSQL_do( db, "BEGIN TRANSACTION" ,status)
   else
      call fSQL_do( db, "BEGIN TRANSACTION" )
   end if

end subroutine ! 1}}}

! subroutine fSQL_commit( db, status) ! {{{ 1
!    Commits a transaction on the given database
! Arguments:
!    db            Structure for the database
!    status        error status
! Note:
!    Accompanies fSQL_begin
!
!-----------------------------------------------------------------!
subroutine fSQL_commit( db, status )
   type(fSQL_DATABASE)         :: db
   type(fSQL_STATUS), optional :: status

   if (present(status)) then
      status%error = 0
      status%errmsg = ' '
   end if

   if (present(status)) then
      call fSQL_do( db, "COMMIT TRANSACTION", status )
   else
      call fSQL_do( db, "COMMIT TRANSACTION" )
   end if

end subroutine ! 1}}}

! subroutine fSQL_rollback(db, status ) {{{1
!    Rolls back any changes to the database since the last commit
! Arguments:
!    db            Structure for the database
!    status        error status
! Note:
!    Accompanies fSQL_begin
!
!-----------------------------------------------------------------!
subroutine fSQL_rollback( db, status )
   type(fSQL_DATABASE)         :: db
   type(fSQL_STATUS), optional :: status

   if (present(status)) then
      status%error = 0
      status%errmsg = ' '
   end if

   if (present(status)) then
      call fSQL_do( db, "ROLLBACK", status )
   else
      call fSQL_do( db, "ROLLBACK" )
   end if

end subroutine ! 1}}}

! subroutine sqlite3_delete_table {{{1
!    Delete a table
! Arguments:
!    db            Structure for the database
!    tablename     Name of the table to be deleted
! Note:
!    The table can not be recovered, unless this
!    is part of a transaction
!
!-----------------------------------------------------------------!
subroutine sqlite3_delete_table( db, tablename )
   type(fSQL_DATABASE) :: db
   character(len=*)    :: tablename

   character(len=20+len(tablename)) :: command

   write( command, "(2A)" ) "DELETE TABLE ", tablename
   call fSQL_do( db, command )

end subroutine ! 1}}} 

! subroutine fSQL_prepare( db, command, stmt, status ) {{{1
!    Prepare a precompiled statatement
!    
! Arguments:
!    dba           Handle to the database
!    stmt          Handle to the prepared statement
!    command       string conataining the SQL command 
!    status        error status
!
!-----------------------------------------------------------------!
subroutine fSQL_prepare( db, command, stmt, status )
   type(fSQL_DATABASE), intent(inout)  :: db
   character(len=*), intent(in)        :: command
   type(fSQL_STATEMENT), intent(inout) :: stmt
   type(fSQL_STATUS), optional         :: status

   integer                             :: count
   integer                             :: i
   character(len=len(command)+1)       :: commandc

   if (present(status)) then
      status%error = 0
      status%errmsg = ' '
   end if

   db%error       = 0
   db%errmsg       = ' '

   if (.not.is_init(stmt%init)) then
     ! write(0,*) "stmt not initialized"
     ! write(0,*) 'init id = ',GET_ID(stmt%init)
     ! write(0,*) 'init data = ',GET_DATA(stmt%init)
      call initialize(stmt%init)
     ! write(0,*) "stmt apres init "
     ! write(0,*) 'init id = ',GET_ID(stmt%init)
     ! write(0,*) 'init data = ',GET_DATA(stmt%init)
      nullify(stmt%column_name)
      nullify(stmt%column_type)
   else
      !write(0,*) "stmt initialezed"
      if (associated(stmt%column_name)) then
         deallocate( stmt%column_name,stat=db%error)
         if (db%error /= FSQL_OK) then
            db%errmsg =  "fSQL_prepare deallocate memory column_name failed "
            write(*,*) db%errmsg 
            return
         endif
         nullify(stmt%column_name)
      endif
      if (associated(stmt%column_type)) then
         deallocate( stmt%column_type,stat=db%error)
         if (db%error /= FSQL_OK) then
            db%errmsg =  "fSQL_prepare deallocate memory column_type failed "
            write(*,*) db%errmsg 
            return
         endif
         nullify(stmt%column_type)
      endif
   end if
   commandc = command
   call stringtoc( commandc )
   db%error = sqlite3_prepare_c( db%db_handle, commandc, stmt%stmt_handle )

   if ( db%error .eq. SQLITE_OK ) then
      call sqlite3_column_count_c( stmt%stmt_handle, stmt%column_count )
      if (DBG) write(*,*) 'nombre de colonnes = ', stmt%column_count



      if (stmt%column_count > 0) then
         allocate(stmt%column_name(stmt%column_count))
         allocate(stmt%column_type(stmt%column_count))
      end if
      count = stmt%column_count


       do i = 1,count
          call sqlite3_column_name_type_c( stmt%stmt_handle, i-1, &
             stmt%column_name(i), stmt%column_type(i) )
          call stringtof( stmt%column_name(i))
          call stringtof( stmt%column_type(i))

      !    select case (columns(i)%type(1:4) )
      !    case( 'INT ', 'INTE' )
      !       columns(i)%type_set = SQLITE_INT4
      !    case( 'FLOA', 'DOUB' )
      !       columns(i)%type_set = SQLITE_REAL8
      !    case( 'CHAR', 'VARC' )
      !       columns(i)%type_set = SQLITE_CHAR
      !    end select

       enddo
   else
      call sqlite3_errmsg_c( db%db_handle, db%errmsg )
   endif
   if (present(status)) then
      status%error      = db%error
      status%errmsg     = db%errmsg
   end if
end subroutine ! 1}}} 

! character(len=LEN_COL_NAME) function fSQL_column_name( stmt, COL_INDEX ) {{{1
!    Return the name of the column COL_INDEX of the result
! Arguments:
!    stmt            Structure for the prepared statement
! Returns:
!    Return the name of the column COL_INDEX of the result
!
!-----------------------------------------------------------------!
character(len=LEN_COL_NAME) function fSQL_column_name( stmt, col_index )
   type(fSQL_STATEMENT), intent(in) :: stmt
   integer, intent(in)              :: col_index

   if (col_index > 0 .AND. col_index <= stmt%column_count) then
      fSQL_column_name = stmt%column_name(col_index)
   end if
end function ! 1}}}

! character(len=LEN_COL_TYPE) function fSQL_column_type( stmt, COL_INDEX ) {{{1
!    Return the declared type of the column COL_INDEX of the result
! Arguments:
!    stmt            Structure for the prepared statement
! Returns:
!    Return the declared type of the column COL_INDEX of the result
!
!-----------------------------------------------------------------!
character(len=LEN_COL_TYPE) function fSQL_column_type( stmt, col_index )
   type(fSQL_STATEMENT), intent(in) :: stmt
   integer, intent(in)              :: col_index

   if (col_index > 0 .AND. col_index <= stmt%column_count) then
      fSQL_column_type = stmt%column_type(col_index)
   end if
end function ! 1}}}

! FSQL  fetch one row
!
! subroutine fSQL_get_row( stmt, finished ) {{{1
!    Get a single row of Data
! Arguments:
!    stmt         db statement
!    finished     logical
! Note:
!    
!    
!
!-----------------------------------------------------------------!
subroutine fSQL_get_row( stmt, finished )
   type(fSQL_STATEMENT)              :: stmt
   logical                           :: finished
   integer                           :: count


   integer                           :: rc
   finished = .false.

   call fSQL_step( stmt, rc )
   if ( rc /= SQLITE_ROW ) then
      finished = .true.
      call fSQL_reset( stmt )
   endif

end subroutine ! 1}}} 

! FSQL precompiled statement finalize
!
! subroutine fSQL_exec_stmt( stmt, status ) {{{1
!    Run the prepared SQL statement
! Arguments:
!    stmt          Handle to the prepared statement
!    completion    Return code, indicating if the command is complete or
!                  not (SQLITE_DONE, SQLITE_MISUSE or SQLITE_ERROR)
!
!-----------------------------------------------------------------!
subroutine fSQL_exec_stmt( stmt, status )
   type(fSQL_STATEMENT), intent(inout)       :: stmt
   type(fSQL_STATUS), optional, intent(out)  :: status

   integer                                   :: completion

   if (present(status)) then
      status%error = 0
      status%errmsg = ' '
   end if

   call sqlite3_step_c( stmt%stmt_handle, completion )
   call fSQL_reset(stmt)
   if (present(status) ) then
      status%error = completion
      status%errmsg= ' in TODO'
   end if

end subroutine ! 1}}}
! subroutine fSQL_finalize( stmt ) {{{1
!    Finalyse precompiled statement
! Arguments:
!    stmt         db statement
! Note:
!    
!-----------------------------------------------------------------!
subroutine fSQL_finalize( stmt )
   type(fSQL_STATEMENT)                      :: stmt

   call sqlite3_finalize_c( stmt%stmt_handle )
   if (.not.is_init(stmt%init)) then
      call initialize(stmt%init)
      nullify(stmt%column_name)
      nullify(stmt%column_type)
   endif
   if (associated(stmt%column_name)) then
      deallocate( stmt%column_name)
      nullify(stmt%column_name)
   endif
   if (associated(stmt%column_type)) then
      deallocate( stmt%column_type)
      nullify(stmt%column_type)
   endif
   stmt%column_count = 0

end subroutine  ! 1}}}

! FSQL precompiled statement used internally
!
! subroutine fSQL_reset( stmt ) {{{1
!    Reset the prepared SQL statement so that it can
!    be used again
! Arguments:
!    stmt          Handle to the prepared statement
!
!-----------------------------------------------------------------!
subroutine fSQL_reset( stmt )
   type(fSQL_STATEMENT)                      :: stmt

   call sqlite3_reset_c( stmt%stmt_handle )

end subroutine  ! 1}}}
! subroutine fSQL_step( stmt, completion ) {{{1
!    Run the prepared SQL statement
! Arguments:
!    stmt          Handle to the prepared statement
!    completion    Return code, indicating if the command is complete or
!                  not (SQLITE_DONE, SQLITE_MISUSE or SQLITE_ERROR)
!
!-----------------------------------------------------------------!
subroutine fSQL_step( stmt, completion )
   type(fSQL_STATEMENT)                        :: stmt
   integer, intent(out)                        :: completion


   call sqlite3_step_c( stmt%stmt_handle, completion )

end subroutine ! 1}}}


! FSQL  fetch many rows 
! in int, int8, real, real8, text types
! count the rows,. With missing management
!
! fSQL_get_many( stmt, nrow, ncols, mode ... ) {{{1
!    
! Arguments:
!    stmt          Handle to the prepared statement
!    NROWS,        total records found
!    NCOLS,        nunber of colums of the result set
!    MODE,         mode must be one of FSQL_INT, FSQL_INT8, FSQL_REAL,
!                                      FSQL_REAL8, FSQL_CHAR
!    char_missing, optional used with FSQL_CHAR mode 
!    int_missing,  optional used with FSQL_INT  mode 
!    int8_missing, optional used with FSQL_INT8 mode 
!    real_missing, optional used with FSQL_REAL mode 
!    real8_missing,optional used with FSQL_REAL8 mode 
!    status        optional error status
!
!-----------------------------------------------------------------!
subroutine fSQL_get_many( stmt, &
                                NROWS,         & 
                                NCOLS,         & 
                                MODE,          & 
                                CHAR_MISSING,  & 
                                INT_MISSING,   & 
                                INT8_MISSING,  & 
                                REAL_MISSING,  & 
                                REAL8_MISSING, & 
                                STATUS )
   type(fSQL_STATEMENT)                        :: STMT
   integer, intent(out)                        :: NROWS
   integer, intent(out)                        :: NCOLS
   integer, intent(in)                         :: MODE
   character(len = *), optional,intent(in)     :: CHAR_MISSING
   integer, optional,intent(in)                :: INT_MISSING
   integer*8, optional,intent(in)              :: INT8_MISSING
   real   ,   optional,intent(in)              :: REAL_MISSING
   real*8   , optional,intent(in)              :: REAL8_MISSING
   type(fSQL_STATUS), optional                 :: STATUS

   type(fSQL_STATUS)                           :: my_status
   character(len=LEN_MISSING)                  :: missingc
   real                                        :: missingr
   real*8                                      :: missingr8
   integer*8                                   :: missingi8
   integer*4                                   :: missingi

   if (present(status)) then
      status%error = 0
      status%errmsg = ' '
   end if

   my_status%error  = -1
   my_status%errmsg = ' '
   select case ( MODE )
   case (FSQL_INT)
      if (present(int_missing) ) then
         missingi = int_missing
         call sqlite3_get_many_intm_c( stmt%stmt_handle,&
         & stmt%ptr_handle,NROWS,NCOLS, my_status%error,missingi )
      else 
         call sqlite3_get_many_int_c( stmt%stmt_handle,&
         & stmt%ptr_handle,NROWS,NCOLS, my_status%error )
      end if

      stmt%row_count    = NROWS
      stmt%column_count = NCOLS
      stmt%ptr_type     = FSQL_INT

   case (FSQL_INT8)
      if (present(int8_missing) ) then
         missingi8 = int8_missing
         call sqlite3_get_many_int8m_c( stmt%stmt_handle,&
         & stmt%ptr_handle,NROWS,NCOLS, my_status%error,missingi8 )

      else 
         call sqlite3_get_many_int8_c( stmt%stmt_handle,&
         & stmt%ptr_handle,NROWS,NCOLS, my_status%error )
      end if

      stmt%row_count    = NROWS
      stmt%column_count = NCOLS
      stmt%ptr_type     = FSQL_INT8

   case (FSQL_REAL)
      if (present(real_missing) ) then
         missingr = real_missing
         call sqlite3_get_many_realm_c( stmt%stmt_handle,&
         & stmt%ptr_handle,NROWS,NCOLS, my_status%error,missingr )
         write(*,*) 'my_status error = ',my_status%error
      else 
         call sqlite3_get_many_real_c( stmt%stmt_handle,&
         & stmt%ptr_handle,NROWS,NCOLS, my_status%error)

      end if

      stmt%row_count    = NROWS
      stmt%column_count = NCOLS
      stmt%ptr_type     = FSQL_REAL

   case (FSQL_REAL8)
      if (present(real8_missing) ) then
         missingr8 = real8_missing
         call sqlite3_get_many_real8m_c( stmt%stmt_handle,&
         & stmt%ptr_handle,NROWS,NCOLS, my_status%error,missingr8 )
      else 
         call sqlite3_get_many_real8_c( stmt%stmt_handle,&
         & stmt%ptr_handle,NROWS,NCOLS, my_status%error )
      end if


      stmt%row_count    = NROWS
      stmt%column_count = NCOLS
      stmt%ptr_type     = FSQL_REAL8

   case (FSQL_CHAR)
      if (present(char_missing)) then
         missingc = char_missing
      else
         missingc = CHAR_NULL_MISSING
      end if
      call stringtoc( missingc )

      call sqlite3_get_many_char_c( stmt%stmt_handle,&
      & stmt%ptr_handle,NROWS,NCOLS, my_status%error , &
      & missingc)
      stmt%row_count    = NROWS
      stmt%column_count = NCOLS
      stmt%ptr_type     = FSQL_CHAR

   case default
      write(*,*) 'unknow mode'
      status%errmsg  = 'unknow mode'
   end select

   if (present(status) )  then
      status%error  = my_status%error
      status%errmsg = ' '
   end if

   if (DBG) write(*,*) 'NROW = ',NROWS,'NCOLS = ',NCOLS
   call fSQL_reset(stmt)
end subroutine ! 1}}} 

! FSQL put the results in a matrix
! called after FSQL_get_many
!
! subroutine fSQL_fill_matrix_int( stmt, int_matrix, status ) {{{1
!    Put the result set in matrix of type int 
!    
! Arguments:
!    stmt          Handle to the prepared statement
!    matrix        matrix where to put the result set
!    status        error status
!
!-----------------------------------------------------------------!
subroutine fSQL_fill_matrix_int( stmt, int_matrix, status )
   type(fSQL_STATEMENT)                    :: stmt
   integer, dimension (:,:), intent(inout) :: int_matrix
   type(FSQL_STATUS), optional             :: status

   integer, dimension(:), pointer          :: vector
   integer                                 :: err
   type(FSQL_STATUS)                       :: my_status

   if (present(status)) then
      status%error = 0
      status%errmsg = ' '
   end if


   my_status%error  = 0
   my_status%errmsg =' ' 
   if ( stmt%ptr_type /= FSQL_INT) then
      my_status%error  = -1
      my_status%errmsg = 'Data fetched is not of type FSQL_INT '
      write(*,*) 'Data fetched is not of type FSQL_INT '
      if (present(status) ) status = my_status
      return
   end if

   allocate(vector(stmt%row_count*stmt%column_count));
   call fsql_vector_int_c(vector,stmt%ptr_handle,stmt%row_count*stmt%column_count,my_status%error)
   if (my_status%error /= 0) then
      write(*,*) "no data"
      my_status%errmsg = 'no Data '
      deallocate(vector)
      if (present(status) ) status = my_status
      return
   end if
   int_matrix = reshape(vector,shape(int_matrix),order = (/2,1/))
   deallocate(vector)
   call fSQL_reset(stmt)
   if (present(status) ) status = my_status
end subroutine ! 1}}} 

! subroutine fSQL_fill_matrix_int8( stmt, int8_matrix, status ) {{{1
!    Put the result set in matrix of type int 
!    
! Arguments:
!    stmt          Handle to the prepared statement
!    matrix        matrix where to put the result set
!    status        error status
!
!-----------------------------------------------------------------!
subroutine fSQL_fill_matrix_int8( stmt, int8_matrix, status )
   type(fSQL_STATEMENT)                      :: stmt
   integer*8, dimension (:,:), intent(inout) :: int8_matrix
   type(FSQL_STATUS), optional               :: status

   integer*8, dimension(:), pointer          :: vector
   integer                                   :: err
   type(FSQL_STATUS)                         :: my_status

   if (present(status)) then
      status%error = 0
      status%errmsg = ' '
   end if

   my_status%error  = 0
   my_status%errmsg =' ' 
   if ( stmt%ptr_type /= FSQL_INT8) then
      my_status%error  = -1
      my_status%errmsg = 'Data fetched is not of type FSQL_INT8 '
      write(*,*) 'Data fetched is not of type FSQL_INT8 '
      if (present(status) ) status = my_status
      return
   end if

   allocate(vector(stmt%row_count*stmt%column_count));
   call fsql_vector_int8_c(vector,stmt%ptr_handle,stmt%row_count*stmt%column_count,my_status%error)
   if (my_status%error /= 0) then
      write(*,*) "no data"
      my_status%errmsg = 'no Data '
      deallocate(vector)
      if (present(status) ) status = my_status
      return
   end if
   int8_matrix = reshape(vector,shape(int8_matrix),order = (/2,1/))
   deallocate(vector)
   call fSQL_reset(stmt)
   if (present(status) ) status = my_status
end subroutine  ! 1}}}

! subroutine fSQL_fill_matrix_real( stmt, real_matrix, status ) {{{1
!    Put the result set in matrix of type real 
!    
! Arguments:
!    stmt          Handle to the prepared statement
!    matrix        matrix where to put the result set
!    status        error status
!
!-----------------------------------------------------------------!
subroutine fSQL_fill_matrix_real( stmt, real_matrix, status )
   type(fSQL_STATEMENT)                 :: stmt
   real, dimension (:,:), intent(inout) :: real_matrix
   type(FSQL_STATUS), optional          :: status

   real, dimension(:), pointer          :: vector
   integer                              :: err
   type(FSQL_STATUS)                    :: my_status

   if (present(status)) then
      status%error = 0
      status%errmsg = ' '
   end if

   my_status%error  = 0
   my_status%errmsg =' ' 
   if ( stmt%ptr_type /= FSQL_REAL) then
      my_status%error  = -1
      my_status%errmsg = 'Data fetched is not of type FSQL_REAL '
      write(*,*) 'Data fetched is not of type FSQL_REAL '
      if (present(status) ) status = my_status
      return
   end if

   if ((stmt%row_count*stmt%column_count) == 0) then
      my_status%error  = -1
      my_status%errmsg = 'nombre enregistrements ou colonnes est nul '
      if (present(status) ) status = my_status
      return
   end if
   allocate(vector(stmt%row_count*stmt%column_count));
   call fsql_vector_real_c(vector,stmt%ptr_handle,stmt%row_count*stmt%column_count,my_status%error)
   if (my_status%error /= 0) then
      write(*,*) "no data"
      my_status%errmsg = 'no Data '
      deallocate(vector)
      if (present(status) ) status = my_status
      return
   end if
   real_matrix = reshape(vector,shape(real_matrix),order = (/2,1/))
   deallocate(vector)
   call fSQL_reset(stmt)
   if (present(status) ) status = my_status
end subroutine  ! 1}}}

! subroutine fSQL_fill_matrix_real8( stmt, real8_matrix, status ) {{{1
!    Put the result set in matrix of type real 
!    
! Arguments:
!    stmt          Handle to the prepared statement
!    matrix        matrix where to put the result set
!    status        error status
!
!-----------------------------------------------------------------!
subroutine fSQL_fill_matrix_real8( stmt, real8_matrix, status )
   type(fSQL_STATEMENT)                   :: stmt
   real*8, dimension (:,:), intent(inout) :: real8_matrix
   type(FSQL_STATUS), optional            :: status

   real*8, dimension(:), pointer          :: vector
   integer                                :: err
   type(FSQL_STATUS)                      :: my_status

   if (present(status)) then
      status%error = 0
      status%errmsg = ' '
   end if

   my_status%error  = 0
   my_status%errmsg =' ' 
   if ( stmt%ptr_type /= FSQL_REAL8) then
      my_status%error  = -1
      my_status%errmsg = 'Data fetched is not of type FSQL_REAL8 '
      write(*,*) 'Data fetched is not of type FSQL_REAL8 '
      if (present(status) ) status = my_status
      return
   end if

   if ((stmt%row_count*stmt%column_count) == 0) then
      my_status%error  = -1
      my_status%errmsg = 'nombre enregistrements ou colonnes est nul '
      if (present(status) ) status = my_status
      return
   end if
   allocate(vector(stmt%row_count*stmt%column_count));
   call fsql_vector_real8_c(vector,stmt%ptr_handle,stmt%row_count*stmt%column_count,my_status%error)
   if (my_status%error /= 0) then
      write(*,*) "no data"
      my_status%errmsg = 'no Data '
      deallocate(vector)
      if (present(status) ) status = my_status
      return
   end if
   real8_matrix = reshape(vector,shape(real8_matrix),order = (/2,1/))
   deallocate(vector)
   call fSQL_reset(stmt)
   if (present(status) ) status = my_status
end subroutine  ! 1}}}

! subroutine fSQL_fill_matrix_char( stmt, char_matrix, status ) {{{
!    Put the result set in matrix of type character 
!    
! Arguments:
!    stmt          Handle to the prepared statement
!    matrix        matrix where to put the result set
!    status        error status
!
!-----------------------------------------------------------------!
subroutine fSQL_fill_matrix_char( stmt, char_matrix, status )
   type(fSQL_STATEMENT)                               :: stmt
   character(len = *), dimension (:,:), intent(inout) :: char_matrix
   type(FSQL_STATUS), optional                        :: status

   character(len = LEN(char_matrix(1,1))), dimension(:), allocatable :: vector
   integer                                            :: err
   type(FSQL_STATUS)                                  :: my_status

   if (present(status)) then
      status%error = 0
      status%errmsg = ' '
   end if

   my_status%error  = 0
   my_status%errmsg =' ' 
   if ( stmt%ptr_type /= FSQL_CHAR) then
      my_status%error  = -1
      my_status%errmsg = 'Data fetched is not of type FSQL_CHAR '
      write(*,*) 'Data fetched is not of type FSQL_CHAR '
      if (present(status) ) status = my_status
      return
   end if

   if ((stmt%row_count*stmt%column_count) == 0) then
      my_status%error  = -1
      my_status%errmsg = 'nombre enregistrements ou colonnes est nul '
      if (present(status) ) status = my_status
      return
   end if
   allocate(vector(stmt%row_count*stmt%column_count));
   call fsql_vector_char_c(vector,stmt%ptr_handle, &
                        & stmt%row_count*stmt%column_count,err)
   if (my_status%error /= 0) then
      write(*,*) "no data"
      my_status%errmsg = 'no Data '
      deallocate(vector)
      if (present(status) ) status = my_status
      return
   end if
   char_matrix = reshape(vector,shape(char_matrix),order = (/2,1/))
   deallocate(vector)
   call fSQL_reset(stmt)
   if (present(status) ) status = my_status
end subroutine ! 1}}} 

! FSQL heap memory management
! called after FSQL_fill_matrix_XXX
!
! fSQL_free_mem( stmt ) {{{1
!    release the memory of the cached result set
!    allocated via fSQL_get_many procedures.
!    
! Arguments:
!    stmt          Handle to the prepared statement
!
!-----------------------------------------------------------------!
subroutine fSQL_free_mem( stmt )
   type(fSQL_STATEMENT)                        :: stmt

   select case ( stmt%ptr_type )
   case (FSQL_INT)
      call fsql_free_int_c(stmt%ptr_handle)
      stmt%ptr_type = 0 
   case (FSQL_INT8)
      call fsql_free_int8_c(stmt%ptr_handle)
      stmt%ptr_type = 0 
   case (FSQL_REAL)
      call fsql_free_real_c(stmt%ptr_handle)
      stmt%ptr_type = 0 
   case (FSQL_REAL8)
      call fsql_free_real8_c(stmt%ptr_handle)
      stmt%ptr_type = 0 
   case (FSQL_CHAR)
      call fsql_free_char_c(stmt%ptr_handle,stmt%row_count*stmt%column_count)
      stmt%ptr_type = 0 
   case default
      write(*,*) 'unknow ptr_type'
   end select
end subroutine ! 1}}} 

! subroutine fSQL_write_matrix_real(a) {{{1
!    dump real matrix to stdout
!    
! Arguments:
!    a             matrix of type real
!
!-----------------------------------------------------------------!
subroutine fsql_write_matrix_real(a)
  real, dimension(:,:) :: a
  integer i,j
  real tmp
  do i = lbound(a,1), ubound(a,1)
     write(*,*)
     do j = lbound(a,2), ubound(a,2)
        write(*,fmt='(1a2,f20.2)', advance = 'no') ' ',a(i,j)
     end do
  end do
end subroutine ! 1}}} 

! subroutine fSQL_write_matrix_real8(a) {{{1
!    dump real matrix to stdout
!    
! Arguments:
!    a             matrix of type real8
!
!-----------------------------------------------------------------!
subroutine fsql_write_matrix_real8(a)
  real*8, dimension(:,:) :: a
  integer i,j
  do i = lbound(a,1), ubound(a,1)
     write(*,*)
     do j = lbound(a,2), ubound(a,2)
        write(*,fmt='(1a2,f20.2)', advance = 'no') '  ',  a(i,j)
     end do
  end do
end subroutine  ! 1}}}

! subroutine fSQL_write_matrix_int(a) {{{1
!    dump real matrix to stdout
!    
! Arguments:
!    a             matrix of type int
!
!-----------------------------------------------------------------!
subroutine fsql_write_matrix_int(a)
  integer, dimension(:,:) :: a
  integer i,j
  do i = lbound(a,1), ubound(a,1)
     write(*,*)
     do j = lbound(a,2), ubound(a,2)
        write(*,fmt='(1a2,i20)', advance = 'no') '  ',  a(i,j)
     end do
  end do
end subroutine  ! 1}}}

! subroutine fSQL_write_matrix_int8(a) {{{1
!    dump real matrix to stdout
!    
! Arguments:
!    a             matrix of type int8
!
!-----------------------------------------------------------------!
subroutine fsql_write_matrix_int8(a)
  integer*8, dimension(:,:) :: a
  integer i,j
  do i = lbound(a,1), ubound(a,1)
     write(*,*)
     do j = lbound(a,2), ubound(a,2)
        write(*,fmt='(1a2,i20)', advance = 'no') '  ',  a(i,j)
     end do
  end do
end subroutine  ! 1}}}

! subroutine fSQL_write_matrix_char(a) {{{1
!    dump real matrix to stdout
!    
! Arguments:
!    a             matrix of type character
!
!-----------------------------------------------------------------!
subroutine fsql_write_matrix_char(a) 
  character(len = *), dimension(:,:) :: a
  integer i,j
  do i = lbound(a,1), ubound(a,1)
     write(*,*)
     do j = lbound(a,2), ubound(a,2)
        write(*,fmt='(1a2,a20)', advance = 'no') '  ',  a(i,j)
     end do
  end do
end subroutine ! 1}}}

end module


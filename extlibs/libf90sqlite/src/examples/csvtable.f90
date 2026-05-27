! csvtable.f90 --
!    Program to read a simple CSV file and put it into a
!    SQLite database, just to demonstrate how the Fortran
!    interface works.
!
!    To keep it simple:
!    - The first line contains the names of the four columns
!    - All lines after that contain the name of the station
!      the date and the two values.
!
program csvtable
   use sqlite

   implicit none

   type(SQLITE_DATABASE)                      :: db
   type(SQLITE_STATEMENT)                     :: stmt
   type(SQLITE_COLUMN), dimension(:), pointer :: column

   integer                                    :: lun = 10
   integer                                    :: i
   integer                                    :: ierr
   character(len=40), dimension(4)            :: name
   real                                       :: salin
   real                                       :: temp
   character(len=40)                          :: station
   character(len=40)                          :: date
   logical                                    :: finished

   !
   ! Read the CSV file and feed the data into the database
   !
   open( lun, file = 'somedata.csv' )
   read( lun, * ) name

   call sqlite3_open( 'somedata.db', db )

   allocate( column(4) )
   call sqlite3_column_props( column(1), name(1), SQLITE_CHAR, 10 )
   call sqlite3_column_props( column(2), name(2), SQLITE_CHAR, 10 )
   call sqlite3_column_props( column(3), name(3), SQLITE_REAL )
   call sqlite3_column_props( column(4), name(4), SQLITE_REAL )
   call sqlite3_create_table( db, 'measurements', column )

   call sqlite3_begin( db )
   do
      read( lun, *, iostat=ierr ) station, date, salin, temp

      if ( ierr .ne. 0 ) exit

      call sqlite3_set_column( column(1), station )
      call sqlite3_set_column( column(2), date    )
      call sqlite3_set_column( column(3), salin   )
      call sqlite3_set_column( column(4), temp    )
      call sqlite3_insert( db, 'measurements', column )

   enddo
   call sqlite3_commit( db )

   close( lun )

   !
   ! We want a simple report, the mean of salinity and temperature
   ! sorted by the station
   !
   ! Note: we need a more convenient routine here ...
   !
   deallocate( column )
   allocate( column(3) )
   call sqlite3_column_props( column(1), 'station', SQLITE_CHAR )
   call sqlite3_column_props( column(2), 'avg('//trim(name(3))//')', SQLITE_REAL )
   call sqlite3_column_props( column(3), 'avg('//trim(name(4))//')', SQLITE_REAL )
!   call sqlite3_prepare_select( db, 'measurements', column, stmt, &
!      'group by station order by station' )
! station,date,salinity,temperature
call sqlite3_prepare( db, 'select station , avg(salinity),avg(temperature) from measurements group by station order by station;', stmt, column)
!              
!call sqlite3_prepare( db, 'select station , salinity from measurements order by station limit 10;', stmt, column)

   write( *, '(3a20)' ) 'Station', 'Mean salinity', 'Mean temperature'
   do
      call sqlite3_next_row( stmt, column, finished )

      if ( finished ) exit

      call sqlite3_get_column( column(1), station )
      call sqlite3_get_column( column(2), salin   )
      call sqlite3_get_column( column(3), temp    )
!temp =0.0

      write( *, '(a20,2f20.3)' ) station, salin, temp
   enddo

   call sqlite3_close( db )
end program

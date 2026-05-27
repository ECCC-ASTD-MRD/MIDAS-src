! airs.f90 --
!    Program to query airs.db 

program airs
   use fSQLite

   implicit none

   type(fSQL_DATABASE)                      :: db
   type(fSQL_STATEMENT)                     :: stmt
   type(fSQL_STATUS)                        :: stat

   integer                                    :: lun = 10
   integer                                    :: i
   integer                                    :: ierr
   character(len=40), dimension(4)            :: name
   character(len=40)            :: var_text 
   integer                                    :: id
   integer*8                                  :: id8,id8m
   integer*8                                  :: miss8
   integer                                    :: canal 
   real*4                                       :: bt 
   real*4                                       :: o_a
   real*4                                       :: o_p
   real*8                                       :: o_p8
   logical                                    :: finished
   integer, dimension(2,5)                    :: ar2d
   integer, dimension(10)                     :: ar1d
   integer, dimension(2)                      :: add
   integer                                    :: r,c
   integer, dimension(:,:), pointer :: matrix
   integer, dimension(:), pointer   :: vector 
   real, dimension(:,:), pointer :: matrix_r
   real*8, dimension(:,:), pointer :: matrix_r8
   integer*8, dimension(:,:), pointer :: matrix_i8
   real, dimension(:), pointer   :: vector_r 
   real*8, dimension(:), pointer   :: vector_r8 
   integer*8, dimension(:), pointer   :: vector_i8 
   integer :: n,m
   character(len = 10), dimension(5)          :: vec_char
   character(len = 20), dimension(:,:),pointer        :: matrix_s 
   character(len  = 500) :: command

   command = " & 
   &  ATTACH DATABASE '/tmp/afsdhmd/airs.dup' AS dup ;   &
   &  CREATE  TABLE dup.toto AS SELECT * FROM airs_body ;&
   &  CREATE  TABLE dup.nono AS SELECT * FROM airs_hdr ; &
   &  CREATE  VIEW dup.V2 AS  SELECT * FROM toto ;"
   write (*,*) command

!   stop
   call fSQL_open( db,'/users/dor/afsd/hmd/nicolas/airs.db', stat)

   call fSQL_do_many( db, command)

   ! a chacune de ses commandes on peut tester si erreur
   if ( fSQL_error(stat) /= 0 ) then
      write(*,*) 'Error: ', fSQL_errmsg(stat)
      stop
   endif

   ! preparer la requette
!   call fSQL_prepare( db, 'SELECT obs_id AS ID,canal_id AS CANAL, BT, O_A, O_P  &
!   & FROM dup.toto WHERE ID = ? AND CANAL BETWEEN ? AND ? LIMIT 10;', stmt)

   call fSQL_prepare( db, 'SELECT obs_id AS ID,canal_id AS CANAL, BT, O_A, O_P  &
   & FROM dup.toto limit 10  ;', stmt, stat)
   if ( fSQL_error(stat) /= 0 ) then
      write(*,*) 'Error mm: ', fSQL_errmsg(stat)
      stop
   endif
   call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR =   5 )
   call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR = 168 )
   call fSQL_bind_param( stmt, PARAM_INDEX = 3, INT_VAR = 169 )
  
if (.false.) then
   write (*,*) 'debut fetch all'
!   call fSQL_get_many_int8(stmt,n,m, INT8_MISSING = int(-99,8), ERR = ierr)
   call fSQL_get_many_int8(stmt,n,m, ERR = ierr)
   write (*,*) 'fin fetch all'
   allocate(matrix_i8(n,m));
   write (*,*) 'debut get all'
   call fSQL_fill_matrix_int8(stmt,matrix_i8)
   write (*,*) 'fin get all'
   call fsql_write_matrix_int8(matrix_i8)
   call fSQL_free_mem(stmt)
   matrix_i8 = 0
   call fSQL_fill_matrix_int8(stmt,matrix_i8)
   call fsql_write_matrix_int8(matrix_i8)
end if

if (.false.) then
   write (*,*) 'debut fetch all'
   call fSQL_get_many_real8(stmt,n,m, REAL8_MISSING = real(-99.0,8))
!   call fSQL_get_many_real(stmt,n,m)
   write (*,*) 'fin fetch all'
   allocate(matrix_r8(n,m));
   write (*,*) 'debut get all'
   call fSQL_fill_matrix_real8(stmt,matrix_r8)
   write (*,*) 'fin get all'
   call fsql_write_matrix_real8(matrix_r8)
   call fSQL_free_mem(stmt)
   matrix_r8 = 0.0
   call fSQL_fill_matrix_real8(stmt,matrix_r8)
   call fsql_write_matrix_real8(matrix_r8)
end if
if (.false.) then
   write (*,*) 'debut fetch all'
!   call fSQL_get_many_real(stmt,n,m,MISSING = 99.0)
   call fSQL_get_many_real(stmt,n,m)
   write (*,*) 'fin fetch all'
   allocate(matrix_r(n,m));
   write (*,*) 'debut get all'
   call fSQL_fill_matrix_real(stmt,matrix_r)
   write (*,*) 'fin get all'
   call fsql_write_matrix_real(matrix_r)
   call fSQL_fill_matrix_real(stmt,matrix_r)
   call fsql_write_matrix_real(matrix_r)
end if

if (.true.) then
   write (*,*) 'debut fetch all'
   call fSQL_get_many_char(stmt,n,m, CHAR_MISSING = '1' , ERR = ierr)
!   call fSQL_get_many_char(stmt,n,m, ERR = ierr)
   write (*,*) 'fin fetch all'
   allocate(matrix_s(n,m));
   allocate(matrix_r(n,m));
   write (*,*) 'debut get all'
   call fSQL_fill_matrix_char(stmt,matrix_s)
   call mat_char_to_real(matrix_r,matrix_s)
   write (*,*) 'fin get all'
   call fsql_write_matrix(matrix_s)
   call fsql_write_matrix(matrix_r)
!   call fSQL_free_mem(stmt)
end if

if (.false.) then
   call fSQL_get_many_int(stmt,n,m,INT_MISSING = -11)
!   call fSQL_get_many_int(stmt,n,m)
   allocate(matrix(n,m));
   allocate(matrix_s(n,m));
   allocate(matrix_r(n,m));
   call fSQL_fill_matrix_char(stmt,matrix_s)
   call fSQL_fill_matrix_real(stmt,matrix_r)
   call fSQL_fill_matrix_int(stmt,matrix)
   call fsql_write_matrix(matrix)
   matrix_r = real(matrix)
   call fsql_write_matrix(matrix_r)
   call fSQL_free_mem(stmt)
   deallocate(matrix);
   allocate(matrix(n,m));
   matrix = 0;
   call fSQL_fill_matrix_int(stmt,matrix)
   call fsql_write_matrix(matrix)
!   call sqlite3_put_c(stmt%ptr_handle)
end if
miss8 =13
write(*,*) 'miss8 = ',miss8
   stop
   write( *, '(5a15)' ) 'ID', 'CANAL', 'BT','O_A','O_P'
   do 
      ! chercher le record suivant
      call fSQL_get_row( stmt, finished )
      if (finished) then
         exit
      endif
      ! On recupere les resulats dans les variables
      call fSQL_get_column_int8( stmt, COL_INDEX = 1 , INT8_VAR  = id8    ) 
!      call fSQL_get_column_char( stmt, COL_INDEX = 1 , CHAR_VAR  = var_text    ) 
      call fSQL_get_column_int( stmt, COL_INDEX = 2 , INT_VAR   = canal ) 
      call fSQL_get_column_real( stmt, COL_INDEX = 3 , REAL_VAR  = bt    ) 
!      call fSQL_get_column_real( stmt, COL_INDEX = 4 , REAL_VAR  = o_a, REAL_MISSING = -77.0   ) 
!      call fSQL_get_column_int8( stmt, COL_INDEX = 4 , INT8_VAR  = id8m, INT8_MISSING = int(13,8)  ) 
      call fSQL_get_column_char( stmt, COL_INDEX = 4 , CHAR_VAR  = var_text, CHAR_MISSING = "TOTO"  ) 
      call fSQL_get_column_real8( stmt, COL_INDEX = 5 , REAL8_VAR = o_p8 , REAL8_MISSING = real(-33.0,8)  ) 
!      call fSQL_get_column_real( stmt, COL_INDEX = 5 , REAL_VAR = o_p,missing = -99.999   ) 
      ! on print vers le stdout
!      write( *, '(1a10,i15,3f20.3)' ) var_text,canal,bt,o_a,o_p8
!      write( *, '(1i10,i15,f20.3,i15,f20.3)' ) id8,canal,bt,id8m,o_p8
      write( *, '(1i10,i15,f20.3,a15,f20.3)' ) id8,canal,bt,var_text,o_p8
!      write( *, '(1i10,i15,f20.3,f20.3,f20.3)' ) id8,canal,bt,o_a,o_p8
   end do

   call fSQL_finalize(stmt)

!!!   call fSQL_prepare( db, 'UPDATE dup.toto SET bt= ?,o_a = ?  &
!!!   & WHERE obs_id = ? AND CANAL_id == ? ;', stmt)
!!!   if ( fSQL_error(db) /= 0 ) then
!!!      write(*,*) 'Error mm: ', fSQL_errmsg(db)
!!!      stop
!!!   end if
!!!   do id = 1,5
!!!      call fSQL_bind_param( stmt, PARAM_INDEX = 1, REAL_VAR =  id*2.0 )
!!!      call fSQL_bind_param( stmt, PARAM_INDEX = 2, REAL_VAR =  id*3.0 )
!!!      call fSQL_bind_param( stmt, PARAM_INDEX = 3, INT_VAR  =   id )
!!!      call fSQL_bind_param( stmt, PARAM_INDEX = 4, INT_VAR  =   168 )
!!!      call fSQL_exec_stmt ( stmt,i)
!!!      write (*,*) i
!!!   end do
!!!   call fSQL_finalize(stmt)
!!!
!!!   ! voir si on a update
!!!   call fSQL_prepare( db, 'SELECT obs_id AS ID,canal_id AS CANAL, BT, O_A, O_P  &
!!!   & FROM dup.toto WHERE CANAL IN (168,169) ;', stmt)
!!!   if ( fSQL_error(db) /= 0 ) then
!!!      write(*,*) 'Error mm: ', fSQL_errmsg(db)
!!!      stop
!!!   endif
!!!   call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR =   5 )
!!!   call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR = 168 )
!!!   call fSQL_bind_param( stmt, PARAM_INDEX = 3, INT_VAR = 169 )
!!!
!!!   write( *, '(5a15)' ) 'ID', 'CANAL', 'BT','O_A','O_P'
!!!   do 
!!!      ! chercher le record suivant
!!!      call fSQL_get_row( stmt, finished )
!!!      if (finished) then
!!!         exit
!!!      endif
!!!      ! On recupere les resulats dans les variables
!!!      call fSQL_get_column( stmt, COL_INDEX = 1 , INT8_VAR  = id8    ) 
!!!      call fSQL_get_column( stmt, COL_INDEX = 2 , INT_VAR   = canal ) 
!!!      call fSQL_get_column( stmt, COL_INDEX = 3 , REAL4_VAR = bt    ) 
!!!      call fSQL_get_column( stmt, COL_INDEX = 4 , REAL4_VAR = o_a   ) 
!!!      call fSQL_get_column( stmt, COL_INDEX = 5 , REAL8_VAR = o_p8   ) 
!!!      ! on print vers le stdout
!!!      write( *, '(2i15,3f20.3)' ) id8,canal,bt,o_a,o_p8
!!!   end do
!!!
!!!   call fSQL_finalize(stmt)
!!!
!!!   ! test insere
!!!
!!!   call fSQL_prepare( db, 'INSERT INTO dup.toto (obs_id ,canal_id , BT, O_A ) VALUES  &
!!!   & (?,?,?,?) ;', stmt)
!!!   if ( fSQL_error(db) /= 0 ) then
!!!      write(*,*) 'Error mm: ', fSQL_errmsg(db)
!!!      stop
!!!   endif
!!!   call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR =   6 )
!!!   call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR = 168 )
!!!   call fSQL_bind_param( stmt, PARAM_INDEX = 3, REAL_VAR = 12.0 )
!!!   call fSQL_bind_param( stmt, PARAM_INDEX = 4, REAL_VAR = 15.0 )
!!!   call fSQL_exec_stmt ( stmt,i)
!!!   call fSQL_finalize(stmt)
!!!
!!!   ! voir si on a insere 
!!!   call fSQL_prepare( db, 'SELECT obs_id AS ID,canal_id AS CANAL, BT, O_A, O_P  &
!!!   & FROM dup.toto WHERE CANAL IN (168,169) ;', stmt)
!!!   if ( fSQL_error(db) /= 0 ) then
!!!      write(*,*) 'Error mm: ', fSQL_errmsg(db)
!!!      stop
!!!   endif
!!!   call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR =   5 )
!!!   call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR = 168 )
!!!   call fSQL_bind_param( stmt, PARAM_INDEX = 3, INT_VAR = 169 )
!!!
!!!   write( *, '(5a15)' ) 'ID', 'CANAL', 'BT','O_A','O_P'
!!!   do 
!!!      ! chercher le record suivant
!!!      call fSQL_get_row( stmt, finished )
!!!      if (finished) then
!!!         exit
!!!      endif
!!!      ! On recupere les resulats dans les variables
!!!      call fSQL_get_column( stmt, COL_INDEX = 1 , INT8_VAR  = id8    ) 
!!!      call fSQL_get_column( stmt, COL_INDEX = 2 , INT_VAR   = canal ) 
!!!      call fSQL_get_column( stmt, COL_INDEX = 3 , REAL4_VAR = bt    ) 
!!!      call fSQL_get_column( stmt, COL_INDEX = 4 , REAL4_VAR = o_a   ) 
!!!      call fSQL_get_column( stmt, COL_INDEX = 5 , REAL8_VAR = o_p8   ) 
!!!      ! on print vers le stdout
!!!      write( *, '(2i15,3f20.3)' ) id8,canal,bt,o_a,o_p8
!!!   end do
!!!
!!!   call fSQL_finalize(stmt)
   call fSQL_do( db, "DETACH DATABASE dup ;" )
   call fSQL_close( db )
!   i = 0
!   do r = 1,2
!      do c =1,5
!         i = i+1
!         ar1d(i) =0 
!      end do
!   end do
!   call print_array_c(ar1d,10)
!   i = 0
!   do r = 1,2
!      do c =1,5
!         i = i+1
!         write(*,*) ar1d(i)
!      end do
!   end do
!   ar2d = reshape(ar1d,(/2,5/),order = (/2,1/))
!   do r = 1,1
!      do c =1,5
!         write(*,*) ar2d(r,c)
!      end do
!   end do
!   allocate(vector(100));
!   call print_array2_c(vector,add,100)
!   do r = 1,1
!      do c =1,5
!         write(*,*) vector(c)
!      end do
!   end do
!   allocate(matrix(10,5));
!   matrix = reshape(vector,shape(matrix),order = (/2,1/))
!   call write_matrix(matrix)
!   do r = 1,1
!      do c =1,5
!         write(*,*) matrix(r,c)
!      end do
!   end do
contains
   !! This subroutine converts a character string to an integer.
   !! adapted for matrix
   subroutine  MAT_CHAR_TO_INTEGER (K, C)
      !! INCOMING: K = an integer variable to receive the converted character string;
      !!           C = the character variable to be converted to an integer.
      integer,  dimension(:,:),intent (out)         :: K
      character (len=*),dimension(:,:), intent (in) :: C
      character (len=LEN(C)) :: TMP

      integer                        :: L
      integer i,j

      L = len (C)

      do i = lbound(k,1), ubound(k,1)
         do j = lbound(k,2), ubound(k,2)
            if (L < 1) then
               print *, "Character value has null length.  Cannot convert to integer."
               cycle
            end  if
            TMP = TRIM(C(i,j))
            
            read (unit=TMP, fmt=*) K(i,j)
         end do
      end do
   end  subroutine  

   !! This subroutine converts a character string to an real.
   !! adapted for matrix
   subroutine  MAT_CHAR_TO_REAL (K, C)
      !! INCOMING: K = an real variable to receive the converted character string;
      !!           C = the character variable to be converted to an integer.
      real   ,  dimension(:,:),intent (out)         :: K
      character (len=*),dimension(:,:), intent (in) :: C
      character (len=LEN(C)) :: TMP

      integer                        :: L
      integer i,j

      L = len (C)

      do i = lbound(k,1), ubound(k,1)
         do j = lbound(k,2), ubound(k,2)
            if (L < 1) then
               print *, "Character value has null length.  Cannot convert to real."
               cycle
            end  if
            TMP = TRIM(C(i,j))
            read (unit=TMP, fmt=*) K(i,j)
         end do
      end do
   end  subroutine  

end program

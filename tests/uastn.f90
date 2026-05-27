! uastn.f90 --
!    Programme qui lira un fichier csv qui contient
!    les stations de radiosondage et leurs lat/lon
!    On va creer un fichier SQLIte ua.db avec 2 tables
!    ua_all et ua_ca. La premiere table contiendra
!    une cle primaire entiere qui identifiera une sation unique
!    de nom name et de latitide lat et de longitude lon.
!    Dans la deuxieme table, on mettra juste les stations
!    Canadiennes (juste la cle ID correpondate Id
!    de la table ua_all et le nom de la station)
!    Cette exercice fait le demo de l'utilsation de l'autoincremnt
!    et la fonction fSQL_last_insert_rwoid().
!    
!    fichier radiosonde.csv
!    "name","lat","lon"
!    "PSGAL",-0.9,-89.6
!    17609,34.88,33.63
!    61995,-20.3,57.5
!    71716,46.9,-71.5
PROGRAM main 
   use fSQLite

   implicit none

   type(fSQL_DATABASE)             :: db
   type(fSQL_STATEMENT)            :: stmt_one
   type(fSQL_STATEMENT)            :: stmt_two

   !! tupe status de l'erreur
   type(fSQL_STATUS)               :: stat

   integer                         :: lun = 10
   integer                         :: i
   integer                         :: ierr
   character(len=40), dimension(3) :: name
   character(len  = 256)           :: query
   character(len  = 256)           :: query_one
   character(len  = 256)           :: query_two
   real                            :: lat
   real                            :: lon
   integer                         :: row_id
   character(len=40)               :: station

   !! matrice de type character pour sauver le resulat de la requette
   character(len = 40), dimension(:,:), pointer:: matrix_char
   integer*4                                   :: n ! nombre de records trouves
   integer*4                                   :: m ! nombre de colonnes du record
   !
   ! Read the CSV file and feed the data into the database
   !
   open( lun, file = 'radiosonde.csv' )
   read( lun, * ) name

   CALL fSQL_open( db,'./ua.db',stat) 
   ! a chacune de ses commandes on peut tester si erreur
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_open : ') 

   ! creer la structure du fichier
   query = 'CREATE TABLE UA_ALL ( ID INTEGER PRIMARY KEY AUTOINCREMENT, name VARCHAR(40), lat REAL, lon REAL);&
         &  CREATE TABLE UA_CA  ( ID INTEGER PRIMARY KEY , name VARCHAR(40), lat REAL, lon REAL);' 

   CALL fSQL_do_many( db,query, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_do_many : ') 

   !! Ajouter des enregistrements
   query_one  = 'INSERT INTO UA_ALL  (name,lat,lon ) VALUES  (?,?,?) ;' 
   query_two  = 'INSERT INTO UA_CA  (ID,name,lat,lon ) VALUES  (?,?,?,?) ;' 

   CALL fSQL_prepare( db, query_one, stmt_one, stat) 
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare stmt_one: ') 

   CALL fSQL_prepare( db, query_two, stmt_two, stat) 
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare stmt_two: ') 

   ! Commencer la transaction 
   CALL fSQL_begin(db) 

   write( *, '(3a20)' ) 'Station', 'latitude', 'longitude'
   do
      read( lun, *, iostat=ierr ) station, lat, lon

      if ( ierr .ne. 0 ) exit
      write( *, '(1a20,2f10.4)' ) station, lat, lon
      CALL fSQL_bind_param( stmt_one, PARAM_INDEX = 1, CHAR_VAR = station )
      CALL fSQL_bind_param( stmt_one, PARAM_INDEX = 2, REAL_VAR = lat )
      CALL fSQL_bind_param( stmt_one, PARAM_INDEX = 3, REAL_VAR = lon )
      ! Executer le statement
      CALL fSQL_exec_stmt ( stmt_one)

      ! si c'est une station qui commence par 71, alors
      ! retrouver son ID de la derniere insertion de la table
      ! ua_all et insere l'enregistement dans la table ua_ca

      if (station(1:2) .EQ. '71' ) then
         row_id = fSQL_last_insert_rowid(db)
         CALL fSQL_bind_param( stmt_two, PARAM_INDEX = 1, INT_VAR = row_id )
         CALL fSQL_bind_param( stmt_two, PARAM_INDEX = 2, CHAR_VAR= station )
         CALL fSQL_bind_param( stmt_two, PARAM_INDEX = 3, REAL_VAR = lat )
         CALL fSQL_bind_param( stmt_two, PARAM_INDEX = 4, REAL_VAR = lon )
         ! Executer le statement
         CALL fSQL_exec_stmt ( stmt_two)
      end if
   enddo

   ! Finir la transaction
   CALL fSQL_commit(db) 

   ! Finalize les statements
   CALL fSQL_finalize(stmt_one)
   CALL fSQL_finalize(stmt_two)

   close( lun )

   ! interoogeons la tables ua_ca pour voir si nous avons insere qqchose
   ! ainsi la distance de cahque station par rapport AEROPORT de DORVAL
   ! AEROPORT D  45.4680555  -73.741111
   query = "SELECT  ID,  NAME, Kms(lat,lon, 45.4680555 ,-73.741111) || ' Kms' FROM ua_ca ORDER BY NAME DESC " 

   !! preparer la requette pour execution
   CALL fSQL_prepare( db, query, stmt_one, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare stmt select: ') 

   ! on cherche les enregistrements les compter et sauver dans le heap 
   CALL fSQL_get_many(stmt_one,NROWS = n,NCOLS = m, MODE = FSQL_CHAR) 

   write (*,*) 'Enregistrements trouves = ',n, 'de ', m, 'Colonnes'

   !! Allouer la matrice avec les bonnes dimensions
   allocate(matrix_char(n,m)); 

   !! Mettre les resultats dans la matrice
   CALL fSQL_fill_matrix(stmt_one,matrix_char)

   !! Sortie des resultats
   write( *, '(2a10,a20)')  'ID', 'NAME','Distance To Dorval'
   CALL fSQL_write_matrix(matrix_char) 
   write( *, *)

   !! liberer les resources 
   CALL fSQL_free_mem(stmt_one) 

   !! Fermer le fichier de Format SQLite 
   CALL fSQL_close( db , stat) 

CONTAINS
subroutine handle_error(stat, message)
   type(FSQL_STATUS)  :: stat
   character(len = *) :: message

   write(*,*) message, fSQL_errmsg(stat)
   CALL fSQL_finalize(stmt_one);
   CALL fSQL_finalize(stmt_two);
   CALL fSQL_close(db);
   stop
end subroutine
END PROGRAM

! query_m1.f90 --
!    Executer un simple Query et recuperer
!    Les enregistrments dans une matrice real4
program query_m1
   use fSQLite

   implicit none

   ! type handle de du fichier SQLIte 
   type(fSQL_DATABASE)                      :: db, db2

   ! tupe statement precompile  du fichier SQLite
   type(fSQL_STATEMENT)                     :: stmt

   ! tupe status de l'erreur
   type(fSQL_STATUS)                        :: stat
   

   ! Declarer les variables de travail du programme 
   ! variables pour le query
   character(len  = 256)            :: query
   integer*4, parameter             :: INTEGER4_KIND_CONST = SELECTED_INT_KIND(8) 
   integer*4, parameter             :: INTEGER8_KIND_CONST = SELECTED_INT_KIND(10) 

   ! matrice de type real pour sauver le resulat de la requette
   real*4 , dimension(:,:), pointer :: matrix_real
   integer*4                        :: n ! nombre de records trouves
   integer*4                        :: m ! nombre de colonnes du record

   integer*4, dimension(:), pointer :: IDs
   integer*4, dimension(:), pointer :: CNs

   ! Variables pour les boucles
   integer*4                        :: I,K
   

   ! Ouvrir le fichier de format SQLite
   CALL fSQL_open( db,'./airs.dup',stat) 
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_open: ') 

   ! soit la requette
   query = "SELECT obs_id AS ID,canal_id AS CANAL,  BT, O_A, O_P  "    
   query = TRIM(query) // ' FROM airs_body WHERE CANAL IN (1,168)  ;' 

   query = "SELECT canal_id as canal, avg(BT), VAR(BT), STDDEV(BT) "    
   query = TRIM(query) // ' FROM airs_body GROUP BY canal  ;' 
   ! preparer la requette pour execution
   CALL fSQL_prepare( db, query, stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare: ') 

   ! on cherche les enregistrements les compter et sauver dans le heap 
   CALL fSQL_get_many(stmt,NROWS = n,NCOLS = m, MODE = FSQL_REAL, REAL_MISSING = 1111.1111) 

   write (*,*) 'Enregestriments trouves = ',n, 'de ', m, 'Colonnes'

   ! Allouer la matrice avec les bonnes dimensions
   allocate(matrix_real(n,m)); 

   ! Mettre les resultats dans la matrice
   CALL fSQL_fill_matrix(stmt,matrix_real) 

   ! Sortie des resultats
!   write( *, '(5a20)' ) 'ID', 'CANAL', 'BT','O_A','O_P'
   write( *, '(5a20)' ) 'Canal', 'AVG(BT)', 'VAR(BT)','STDDEV(BT)'

   CALL fSQL_write_matrix(matrix_real) 

   ! liberer les resources 
   CALL fSQL_free_mem(stmt) 
   
   ! On peut Recuperer les IDS de la matrice dans un vecteur
   allocate(IDs(n))
   IDs = INT(matrix_real(:,1), KIND = INTEGER4_KIND_CONST) 
   write(*,*)
   write(*,'(i10)') (IDs(I), I = lbound(IDs,1), ubound(IDs,1))

   ! On peut Recuperer les CNs de la matrice dans un vecteur
   allocate(CNs(n))
   CNs = INT(matrix_real(:,2), KIND =  INTEGER4_KIND_CONST) 

   ! Sortie stdout
   write(*,*)
   write(*,'(i10)') (CNs(I), I = lbound(CNs,1), ubound(CNs,1))

   ! Deallocate la matrice et les vecteurs
   deallocate(matrix_real) 
   deallocate(IDs)
   deallocate(CNs)
   
   ! Finaliser et fermer le fichier 
   ! finalise la requette pour que stmt puisse etre reutilise
   ! avec une requette (ie avec fSQL_prepare()) et aussi 
   ! l'invoquer avant de fermer la bd
   CALL fSQL_finalize( stmt )

   ! Fermer le fichier de Format SQLite 
   CALL fSQL_close( db , stat)
   
contains
   subroutine handle_error(stat, message)
   type(FSQL_STATUS)  :: stat
   character(len = *) :: message
      write(*,*) message, fSQL_errmsg(stat)
      stop
   end subroutine
end program

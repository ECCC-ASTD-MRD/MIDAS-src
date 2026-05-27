! query_m2.f90 --
!    Creer un nouveau fichier airs.dup2
!    Executer un script sql via la commande fSQl_do_script
!    Qui va creer 2 tables (voir script.sql)
!    La table nono_hdr (mot nono emprunte de chez Pierre.K pour changer
!    des titi et toto)
!    a une cle autoincrement, donc
!    l'engin sql s'occpera de lui aasigner une valeure entiere
!    lors de l'insertion des enregistrements
!    Ensuite on fait une requette pour avoir
!    Les enregistrments dans une matrice real4
program query_m6
   use fSQLite

   implicit none

   ! type handle de du fichier SQLIte 1
   type(fSQL_DATABASE)                      :: db

   ! tupe statement precompile  du fichier SQLite
   type(fSQL_STATEMENT)                     :: stmt

   ! tupe status de l'erreur
   type(fSQL_STATUS)                        :: stat
   

   ! Declarer les variables de travail du programme 1
   ! variables pour le query
   character(len  = 256)            :: query
   integer*4, parameter             :: INTEGER4_KIND_CONST = SELECTED_INT_KIND(8) 
   integer*4, parameter             :: INTEGER8_KIND_CONST = SELECTED_INT_KIND(10) 

   character(len  = 256)            :: one_cmd
   character(len  = 500)            :: many_cmd

   ! matrice de type real pour sauver le resulat de la requette
   real*4 , dimension(:,:), pointer :: matrix_real
   integer*4                        :: n ! nombre de records trouves
   integer*4                        :: m ! nombre de colonnes du record

   integer*4, dimension(:), pointer :: IDs
   integer*4, dimension(:), pointer :: CNs

   ! Variables pour les boucles
   integer*4                        :: I,K
   
   ! Ouvrir le fichier de format SQLite
   CALL fSQL_open( db,'./airs.dup2',stat) 
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_open airs.dup2: ') 

   CALL fSQL_do_script(db,'script.sql',stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_do_script db: ') 

   ! Executer +eurs commandes SQL a la fois
   
   many_cmd ="&
   &  INSERT INTO  nono_hdr  (airs_id,latitude) VALUES (1,2.0) ;&
   &  INSERT INTO  nono_hdr  (airs_id,latitude) VALUES (5,5.0) ;&
   &  INSERT INTO  nono_hdr  (airs_id,latitude) VALUES (9,9.0) ;"

   CALL fSQL_do_many( db, many_cmd, stat) 
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_do_many db: ') 
   
   ! soit la requette

   query = "SELECT airs_id,latitude, longitude, date, time FROM nono_hdr;  "    

   ! preparer la requette pour execution
   CALL fSQL_prepare( db, query, stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare db: ') 


   ! on cherche les enregistrements les compter et sauver dans le heap 
   CALL fSQL_get_many(stmt,NROWS = n,NCOLS = m, MODE = FSQL_REAL, REAL_MISSING = -9999.99) 

   write (*,*) 'Enregestriments trouves = ',n, 'de ', m, 'Colonnes'

   ! Allouer la matrice avec les bonnes dimensions
   allocate(matrix_real(n,m)); 

   ! Mettre les resultats dans la matrice
   CALL fSQL_fill_matrix(stmt,matrix_real) 

   ! Sortie des resultats

   do I= 1, m 
      write( *, fmt='(1a2,a20)', advance = 'no' ) ' ',fSQL_column_name(stmt,I) 
   end do

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

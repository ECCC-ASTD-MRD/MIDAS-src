! query_m3.f90 --
!    Executer un simple Query et recuperer
!    Les enregistrments dans une matrice ide chaines de caracteres
program query_m3
   use fSQLite

   implicit none

   ! Declarer les types du module fSQLite 
   ! type handle de du fichier SQLIte
   type(fSQL_DATABASE)                      :: db

   ! tupe statement precompile  du fichier SQLite
   type(fSQL_STATEMENT)                     :: stmt

   ! tupe status de l'erreur
   type(fSQL_STATUS)                        :: stat
   !

   ! Declarer les variables de travail du programme 
   ! variables pour le query
   character(len  = 256)                       :: query

   ! matrice de type character pour sauver le resulat de la requette
   character(len = 40), dimension(:,:), pointer:: matrix_char
   integer*4                                   :: n ! nombre de records trouves
   integer*4                                   :: m ! nombre de colonnes du record

   integer*4, dimension(:), pointer            :: IDs
   integer*4, dimension(:), pointer            :: CNs
   real*4, dimension(:), allocatable           :: BTs

   ! Variables pour les boucles
   integer*4                                   :: I,K
   !

   ! Ouvrir le fichier de format SQLite
   CALL fSQL_open( db,'./airs.dup',stat) 
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_open: ') 
   

   ! Selectionner les enregistrements
   query = "SELECT obs_id AS ID,canal_id AS CANAL, BT, O_A, O_P "    
   query = TRIM(query) // ' FROM airs_body WHERE CANAL IN (1,168)  ;' 

   ! preparer la requette pour execution
   CALL fSQL_prepare( db, query, stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare: ') 

   ! on cherche les enregistrements les compter et sauver dans le heap 
   CALL fSQL_get_many(stmt,NROWS = n,NCOLS = m, MODE = FSQL_CHAR, CHAR_MISSING = "WHITEHORSE") 

   write (*,*) 'Enregistrements trouves = ',n, 'de ', m, 'Colonnes'

   ! Allouer la matrice avec les bonnes dimensions
   allocate(matrix_char(n,m));  
   

   ! Mettre les resultats dans la matrice
   CALL fSQL_fill_matrix(stmt,matrix_char) 

   ! Sortie des resultats
   write( *, '(5a10)' ) 'ID', 'CANAL', 'BT','O_A','O_P'
   CALL fSQL_write_matrix(matrix_char)  
   

   ! liberer les resources 
   CALL fSQL_free_mem(stmt) ! {{{  

   ! On peut Recuperer les IDS de la matrice dans un vecteur
   allocate(IDs(n))
   CALL VEC_CHAR_TO_INTEGER(IDs, matrix_char(:,1)) 
   write(*,*)
   write(*,'(i10)') (IDs(I), I = lbound(IDs,1), ubound(IDs,1))
   

   ! On peut Recuperer les CNs de la matrice dans un vecteur
   allocate(CNs(n))
   CALL VEC_CHAR_TO_INTEGER(CNs, matrix_char(:,2)) 

   ! Sortie stdout
   write(*,*)
   write(*,'(i10)') (CNs(I), I = lbound(CNs,1), ubound(CNs,1))
   

   ! On peut Recuperer les BTs de la matrice dans un vecteur
   allocate(BTs(n))
   CALL VEC_CHAR_TO_REAL(BTs, matrix_char(:,3)) 

   ! Sortie stdout
   write(*,*)
   write(*,'(f10.4)') (BTs(I), I = lbound(BTs,1), ubound(BTs,1))
   

   ! Deallocate la matrice et les vecteurs 
   deallocate(matrix_char) 
   deallocate(IDs)
   deallocate(CNs)
   deallocate(BTs)
   

   ! Finaliser et fermer le fichier {{{1
   ! finalise la requette pour que stmt puisse etre reutilise
   ! avec une requette (ie avec fSQL_prepare()) et aussi 
   ! l'invoquer avant de fermer la bd
   CALL fSQL_finalize( stmt )

   ! Fermer le fichier de Format SQLite 
   CALL fSQL_close( db , stat)
   

contains
   ! This subroutine converts a character string to an integer.
   ! adapted for a vector 
   subroutine  VEC_CHAR_TO_INTEGER (K, C) 
      ! INCOMING: K = an integer variable to receive the converted character string;
      !           C = the character variable to be converted to an integer.
      integer,  dimension(:),intent (out)         :: K
      character (len=*),dimension(:), intent (in) :: C
      character (len=LEN(C)) :: TMP

      integer                        :: L
      integer i,j

      L = len (C)

      do i = lbound(k,1), ubound(k,1)
         if (L < 1) then
            print *, "Character value has null length.  Cannot convert to integer."
            cycle
         end  if
         TMP = TRIM(C(i))
         
         read (unit=TMP, fmt=*) K(i)
      end do
   end  subroutine  
   

   ! This subroutine converts a character string to an real.
   ! adapted for a vector
   subroutine  VEC_CHAR_TO_REAL (K, C) 
      ! INCOMING: K = an real variable to receive the converted character string;
      !           C = the character variable to be converted to a real
      real   ,  dimension(:),intent (out)         :: K
      character (len=*),dimension(:), intent (in) :: C
      character (len=LEN(C)) :: TMP

      integer                        :: L
      integer i,j

      L = len (C)

      do i = lbound(k,1), ubound(k,1)
         if (L < 1) then
            print *, "Character value has null length.  Cannot convert to real."
            cycle
         end  if
         TMP = TRIM(C(i))
         read (unit=TMP, fmt=*) K(i)
      end do
   end  subroutine 
   

   subroutine handle_error(stat, message)
   type(FSQL_STATUS)  :: stat
   character(len = *) :: message
      write(*,*) message, fSQL_errmsg(stat)
      stop
   end subroutine
end program

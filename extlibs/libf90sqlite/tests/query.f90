! query_one.f90 --
!    Executer un simple Query et recuperer
!    un enregistrement un a la fois au fur
!    a mesure que les resultats sont trouves.
!    On utilse le fichier le fihier airs.dup
!    genere par le programme open_close.f90
!    Utilise le statement precompile

program query_one
   use fSQLite

   implicit none

   ! Declarer les types du module fSQLite 
   ! 
   ! type handle de du fichier SQLIte
   type(fSQL_DATABASE)                      :: db

   ! tupe statement precompile  du fichier SQLite
   type(fSQL_STATEMENT)                     :: stmt

   ! tupe status de l'erreur
   type(fSQL_STATUS)                        :: stat

   ! Declarer les variables de travail du programme 
   ! 
   ! variables pour une commande ou +eurs commandes
   character(len  = 256)                    :: one_cmd
   character(len  = 500)                    :: many_cmd

   ! variables pour le query
   character(len  = 256)                    :: query
   logical                                  :: finished
   integer*4                                :: ID  
   integer*4                                :: CANAL
   real*4                                   :: BT  
   real*4                                   :: O_A
   real*4                                   :: O_P


   ! Ouvrir le fichier de format SQLite
   CALL fSQL_open( db,'./airs.dup',stat) 
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_open: ') 


   ! Selectionner les enregistrements
   query = 'SELECT obs_id AS ID,canal_id AS CANAL, BT, O_A, O_P  ' ! 
   query = TRIM(query) // ' FROM airs_body WHERE CANAL IN (1,168)  ;'

   ! preparer la requette pour execution
   CALL fSQL_prepare( db, query, stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare stmt select: ') 


   ! on cherche les enregistrements aufur et a mesure qu'on
   ! les trouve

   write( *, '(5a15)' ) 'ID', 'CANAL', 'BT','O_A','O_P'
   do 
      ! chercher le record suivant
      CALL fSQL_get_row( stmt, finished )
      if (finished) then
         exit
      endif

      ! On recupere les resulats dans les variables
      CALL fSQL_get_column( stmt, COL_INDEX = 1, INT_VAR   = ID    ) 
      CALL fSQL_get_column( stmt, COL_INDEX = 2, INT_VAR   = CANAL ) 
      CALL fSQL_get_column( stmt, COL_INDEX = 3, REAL_VAR  = BT    ) 
      CALL fSQL_get_column( stmt, COL_INDEX = 4, REAL_VAR  = O_A, REAL_MISSING = -77.0 ) 
      CALL fSQL_get_column( stmt, COL_INDEX = 5, REAL_VAR  = O_P, REAL_MISSING = -99.0 ) 
      ! on print vers le stdout
      write( *, '(1i10,i15,3f20.3)' ) ID,CANAL,BT,O_A,O_P
   end do

   ! finalise la requette pour que stmt puisse etre reutilise
   ! avec une requette (ie avec fSQL_prepare()) et aussi 
   ! l'invoquer avant de fermer la bd
   CALL fSQL_finalize( stmt )

   ! Fermer le fichier de Format SQLite 
   CALL fSQL_close( db , stat) ! 

contains
   subroutine handle_error(stat, message)
   type(FSQL_STATUS)  :: stat
   character(len = *) :: message
      write(*,*) message, fSQL_errmsg(stat)
      CALL fSQL_finalize( stmt )
      stop
   end subroutine
end program

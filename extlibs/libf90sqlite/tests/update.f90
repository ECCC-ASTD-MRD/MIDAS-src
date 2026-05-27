! update.f90 --
!    Executer un simple Query et recuperer
!    un enregistrement un a la fois au fur
!    a mesure que les resultats sont trouves.
!    mettre a jour quelques enregs et faire un query
!    Executer ce programme apres le programme insere.f90
!    Car c'est ce sont les enregistrement ajoutes par insere.f90
!    qui seront mise a jour pour fin de demo.

program update
   use fSQLite

   implicit none

   ! Declarer les types du module fSQLite 
   ! type handle de du fichier SQLIte
   type(fSQL_DATABASE)                      :: db

   ! tupe statement precompile  du fichier SQLite
   type(fSQL_STATEMENT)                     :: stmt

   ! tupe status de l'erreur
   type(fSQL_STATUS)                        :: stat
   

   ! Declarer les variables de travail du programme 
   ! variables pourune commande ou +eurs commandes
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

   integer*4                                :: i_id
   integer*4                                :: i_cn
   

   ! Ouvrir le fichier de format SQLite
   CALL fSQL_open( db,'./airs.dup',stat) 
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_open : ') 

   ! Selectionner les enregistrements
   query = 'SELECT obs_id AS ID,canal_id AS CANAL, BT, O_A, O_P  '  
   query = TRIM(query) // ' FROM airs_body WHERE CANAL IN (1,168)  ;'

   ! preparer la requette pour execution
   CALL fSQL_prepare( db, query, stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare : ') 

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
   

   ! Commencer la transaction 
   CALL fSQL_begin(db) 

   ! Mettre a jour les enregistrements
   query = 'UPDATE airs_body SET  O_A = ?, O_P = ? WHERE obs_id = ? AND canal_id = ? ;' 

   CALL fSQL_prepare( db, query, stmt, stat) 
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare : ') 

   ! on update les obs 6,7 et 8 pour chacun des
   ! canaux 1 et 168
   do i_id = 6,8
      do i_cn = 1,168,167
         CALL fSQL_bind_param( stmt, PARAM_INDEX = 1, REAL_VAR = 1.777 )
         CALL fSQL_bind_param( stmt, PARAM_INDEX = 2, REAL_VAR = -2.888 )
         CALL fSQL_bind_param( stmt, PARAM_INDEX = 3, INT_VAR  = i_id )
         CALL fSQL_bind_param( stmt, PARAM_INDEX = 4, INT_VAR  = i_cn )
         ! Executer le statement
         CALL fSQL_exec_stmt ( stmt)
      end do
   end do
   ! Finalize statement
   CALL fSQL_finalize(stmt)
   

   ! Finir la transaction
   CALL fSQL_commit(db) 

   ! Selectionner les enregistrements
   query = 'SELECT obs_id AS ID,canal_id AS CANAL, BT, O_A, O_P  '  
   query = TRIM(query) // ' FROM airs_body WHERE CANAL IN (1,168)  ;'

   ! preparer la requette pour execution
   CALL fSQL_prepare( db, query, stmt, stat)
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_prepare : ') 

   ! on cherche les enregistrements au fur et a mesure qu'on
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
   CALL fSQL_close( db , stat)  
   

contains
   subroutine handle_error(stat, message)
   type(FSQL_STATUS)  :: stat
   character(len = *) :: message
      write(*,*) message, fSQL_errmsg(stat)
      CALL fSQL_finalize(stmt);
      stop
   end subroutine
end program

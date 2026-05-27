! query2.f90 --
!    Executer un simple Query et recuperer
!    un enregistrement un a la fois au fur
!    a mesure que les resultats sont trouves.
!    utilise le statement precompile avec des parametres
!    Que l"on remplace au runtime

program query2
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
   ! 


   ! Ouvrir le fichier de format SQLite
   call fSQL_open( db,'./airs.dup',stat)
   if ( fSQL_error(stat) /= FSQL_OK ) call handle_error(stat,'fSQL_open: ') 

   ! soit la requette
   query = 'SELECT obs_id AS ID,canal_id AS CANAL, BT, O_A, O_P  &
   & FROM airs_body WHERE CANAL BETWEEN ? AND ?  ;'

   ! preparer la requette pour execution
   call fSQL_prepare( db, query, stmt, stat) 
   if ( fSQL_error(stat) /= FSQL_OK ) call handle_error(stat,'fSQL_prepare: ') 


   ! test bindings precompiled statements 
   call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR =   168 )
   call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR =   170 )

   ! on cherche les enregistrements aufur et a mesure qu'on
   ! les trouve

   write( *, *) ' Apres un 1er bindings'
   write( *, '(5a15)' ) 'ID', 'CANAL', 'BT','O_A','O_P'

   do 
      ! chercher le record suivant
      call fSQL_get_row( stmt, finished )
      if (finished) then
         exit
      endif

      ! On recupere les resulats dans les variables
      call fSQL_get_column( stmt, COL_INDEX = 1, INT_VAR   = ID    ) 
      call fSQL_get_column( stmt, COL_INDEX = 2, INT_VAR   = CANAL ) 
      call fSQL_get_column( stmt, COL_INDEX = 3, REAL_VAR  = BT    ) 
      call fSQL_get_column( stmt, COL_INDEX = 4, REAL_VAR  = O_A, REAL_MISSING = -77.0 ) 
      call fSQL_get_column( stmt, COL_INDEX = 5, REAL_VAR  = O_P, REAL_MISSING = -99.0 ) 
      ! on print vers le stdout
      write( *, '(1i10,i15,3f20.3)' ) ID,CANAL,BT,O_A,O_P
   end do

   ! Faisons un bind avec d'autres params 
   call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR =   21 )
   call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR =   21 )

   ! on cherche les enregistrements aufur et a mesure qu'on
   ! les trouve

   write( *, *) ' Apres un autre bindings'
   write( *, '(5a15)' ) 'ID', 'CANAL', 'BT','O_A','O_P'

   do 
      ! chercher le record suivant
      call fSQL_get_row( stmt, finished )
      if (finished) then
         exit
      endif

      ! On recupere les resulats dans les variables
      call fSQL_get_column( stmt, COL_INDEX = 1, INT_VAR   = ID    ) 
      call fSQL_get_column( stmt, COL_INDEX = 2, INT_VAR   = CANAL ) 
      call fSQL_get_column( stmt, COL_INDEX = 3, REAL_VAR  = BT    ) 
      call fSQL_get_column( stmt, COL_INDEX = 4, REAL_VAR  = O_A, REAL_MISSING = -77.0 ) 
      call fSQL_get_column( stmt, COL_INDEX = 5, REAL_VAR  = O_P, REAL_MISSING = -99.0 ) 
      ! on print vers le stdout
      write( *, '(1i10,i15,3f20.3)' ) ID,CANAL,BT,O_A,O_P
   end do 

   ! on peut changer la requette et les bindings 
   query = 'SELECT obs_id AS ID,canal_id AS CANAL, BT, O_A, O_P  &
   & FROM airs_body WHERE ID = ? AND CANAL == ?   ;'

   ! preparer la requette pour execution 
   call fSQL_prepare( db, query, stmt, stat) 
   if ( fSQL_error(stat) /= FSQL_OK ) call handle_error(stat,'fSQL_prepare: ') 

   call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR =   5 )
   call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR =   21 )

   ! on cherche les enregistrements aufur et a mesure qu'on
   ! les trouve

   write( *, *) ' Apres changement de requette et les params'
   write( *, '(5a15)' ) 'ID', 'CANAL', 'BT','O_A','O_P'
   do 
      ! chercher le record suivant
      call fSQL_get_row( stmt, finished )
      if (finished) then
         exit
      endif
      ! On recupere les resulats dans les variables
      call fSQL_get_column( stmt, COL_INDEX = 1, INT_VAR   = ID    ) 
      call fSQL_get_column( stmt, COL_INDEX = 2, INT_VAR   = CANAL ) 
      call fSQL_get_column( stmt, COL_INDEX = 3, REAL_VAR  = BT    ) 
      call fSQL_get_column( stmt, COL_INDEX = 4, REAL_VAR  = O_A, REAL_MISSING = -77.0 ) 
      call fSQL_get_column( stmt, COL_INDEX = 5, REAL_VAR  = O_P, REAL_MISSING = -99.0 ) 

      ! on print vers le stdout
      write( *, '(1i10,i15,3f20.3)' ) ID,CANAL,BT,O_A,O_P
   end do 



   ! finalise la requette pour que stmt puisse etre reutilise
   ! avec une requette (ie avec fSQL_prepare()) et aussi 
   ! l'invoquer avant de fermer la bd
   call fSQL_finalize( stmt )

   ! Fermer le fichier de Format SQLite 
   call fSQL_close( db , stat)

contains
   subroutine handle_error(stat, message)
   type(FSQL_STATUS)  :: stat
   character(len = *) :: message
      write(*,*) message, fSQL_errmsg(stat)
      stop
   end subroutine
end program

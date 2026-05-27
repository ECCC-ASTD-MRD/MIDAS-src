! open_close.f90
!    open_close fichier SQLite
!    process une ou +eurs  commandes SQL
!    Ce programme sert a genere le fichier SQlite airs.dup
!    a partir du fichier airs.db
!    Montrer comment ouvrir un fichier sqlite, executer 
!    soit une commande ou plusirurs commandes a la fois 
!    Le fichier airs.dup sera utilise dans la suite des
!    autres programmes de demo

program open_close
   use fSQLite

   implicit none

   !! Declarer les types du module fSQLite 
   !! type handle de du fichier SQLIte
   type(fSQL_DATABASE)                      :: db

   !! type statement precompile  du fichier SQLite
   type(fSQL_STATEMENT)                     :: stmt

   !! tupe status de l'erreur
   type(fSQL_STATUS)                        :: stat

   !! Declarer les variables de travail du programme 
   ! variables pourune commande ou +eurs commandes
   character(len  = 256)                    :: one_cmd
   character(len  = 500)                    :: many_cmd

   !! Ouvrir le fichier de format SQLite
   CALL fSQL_open( db,'./airs.db',stat) ! 
   ! a chacune de ses commandes on peut tester si erreur
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_open : ') 

   !! Executer une seule commande SQL
   one_cmd = "ATTACH DATABASE './airs.dup' AS dup ;" ! 
   CALL fSQL_do( db, one_cmd , stat )
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_do : ') 

   !! Executer +eurs commandes SQL a la fois

   many_cmd ="&
   &  CREATE  TABLE dup.airs_body AS SELECT * FROM airs_body ;&
   &  CREATE  TABLE dup.airs_hdr AS SELECT * FROM airs_hdr ;&
   &  CREATE  VIEW dup.V1 AS  SELECT * FROM airs_body WHERE canal_id IN (1,21) ;" 

   CALL fSQL_do_many( db, many_cmd, stat) 
   if ( fSQL_error(stat) /= FSQL_OK ) CALL handle_error(stat,'fSQL_do_many : ') 


   !! detacher la bd dup
   CALL fSQL_do( db, "DETACH DATABASE dup ;" ) ! 

   !! Fermer le fichier de Format SQLite 
   CALL fSQL_close( db ) ! 

contains
   subroutine handle_error(stat, message)
   type(FSQL_STATUS)  :: stat
   character(len = *) :: message
      write(*,*) message, fSQL_errmsg(stat)
      stop
   end subroutine
end program

#include <stdio.h>
#include <regex.h>
#include <assert.h>
#include <unistd.h> /* to get the function 'access' to know if files are already existing */

/* Include pour sqlite3 */
#include <sqlite3.h>

/* Include pour la librairie App dans 'rpn/libs' */
#include "App.h"
/* Include pour la librairie RPN */
#include "rmn.h"

/* Include pour ma librairie de manipulation des fichiers standard RPN */
#include "fstdlib.h"
/* Include pour les constantes OK et NOT_OK */
#include "ok_or_notok.h"
/* Include pour la structure qui definit toutes les options */
#include "options.h"

/* nom des fonctions SQL qui seront creees pour rechercher les observations
 * dans la base de donnees
 */
#define SQLFUNCTION_NAME       "checkgrid_sql"
#define NUMBER_OF_ARGS_FOR_CHECK_GRID         13
#define SQL_VERTICAL_NAME      "checkvertical_sql"
#define NUMBER_OF_ARGS_FOR_CHECK_VERTICAL     4
#define SQL_VERTICAL_GZ_NAME   "checkvertical_gz_sql"
#define NUMBER_OF_ARGS_FOR_CHECK_VERTICAL_GZ  9

#define SQL_BUFFER_SIZE   19865

static int SPLITOBS_SQLITE_VERBOSE = 0;

typedef struct { // to be used as the argument in a callback
  char* split_on_key;
  char* table_list_with_split_key;
  char* table_list_without_split_key;
} sqlite_get_tables_callback_arg;

/*****************************************************/
/******* Prototype des fonctions localement **********/
/*****************************************************/
static int sqlite_schema_callback(void *schema_void, int count, char **data, char **columns);
static int sqlite_get_tables_callback(void *callback_args, int count, char **data, char **columns);

int sqlite_get_tables(char* obsin, char* split_on_key, char* table_list_with_split_key, char* table_list_without_split_key);

void append_table_list_split_key_requests_nsplit(char* requete_sql, char* attached_db_name, char* table_list, char* split_on_key, int nsplit, int id);
void append_table_list_split_key_requests_using_header(char* requete_sql, char* attached_db_name, char* table_list, char* split_on_key, char* header_table, char* data_table);
void append_table_list_without_split_key_requests(char* requete_sql, char* attached_db_name, char* table_list);

void checkgrid_sql(sqlite3_context *context, int argc, sqlite3_value **argv);

void checkvertical_sql(sqlite3_context *context, int argc, sqlite3_value **argv);

void checkvertical_gz_sql(sqlite3_context *context, int argc, sqlite3_value **argv);

/*****************************************************/
/******* Fin des prototype des fonctions *************/
/*****************************************************/

/* Fonction principale pour le splitting en format SQL
 */
int splitobs_sql(options opt, gridtype grid, gridtype grid_gz, int VERBOSE,
                 float* VALEURS_GZ_MIN, float* VALEURS_GZ_MAX) {
  int status, EXIT_STATUS;
  char sqlreq_tables_without_split_key[MAXSTR];
  char table_list_with_split_key[MAXSTR];
  char table_list_without_split_key[MAXSTR];
  sqlite3  *sqldb;
  char *ErrMsg, requete_sql[SQL_BUFFER_SIZE];

  SPLITOBS_SQLITE_VERBOSE = VERBOSE;

  /**********************************************************
   * Cette partie a pour but de manipuler la base de donnees SQL
   * et d'en creer une nouvelle qui ne contient que les observations a
   * l'interieur du domaine defini par la grille donnee plus haut.
   **********************************************************/

  strcpy(table_list_with_split_key,"");
  strcpy(table_list_without_split_key,"");
  status = sqlite_get_tables(opt.obsin, opt.rdb_split_on_key, table_list_with_split_key, table_list_without_split_key);
  if( status != OK ) {
    App_Log(APP_ERROR,"Fonction main: Erreur %d de la fonction sqlite_get_tables pour le fichier '%s'\n", status, opt.obsin);
    return NOT_OK;
  }

  append_table_list_without_split_key_requests(sqlreq_tables_without_split_key,"dbin",table_list_without_split_key);

  /* Si on n'est pas en mode round-robin, alors les fichiers 'grid*.gridid' ont ete ouverts */
  if ( opt.roundrobin == 0 ) {
    /* On ouvre la base de donnees SQL finale */
    status = access(opt.obsout,F_OK);
    if ( status == 0 ) { /* Le fichier existe deja */
      status = sqlite3_open(opt.obsout,&sqldb);
      if ( status != SQLITE_OK ) {
        App_Log(APP_ERROR,"Fonction main: Incapable d'ouvrir la base de donnees: %s\n", sqlite3_errmsg(sqldb));

        status = sqlite3_close(sqldb);
        if( status != SQLITE_OK )
          App_Log(APP_ERROR,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

        return NOT_OK;
      } /* Fin du 'if ( status != SQLITE_OK )' */
    } /* Fin du 'if ( status == 0 )' */
    else {
      /* Si le fichier n'existe pas encore alors on doit le creer et construire le meme schema */
      sqlite3 *sqldbin;
      char sqlschema[MAXSTR*32];

      /* On lit le fichier d'input */
      status = sqlite3_open(opt.obsin,&sqldbin);
      if ( status != SQLITE_OK ) {
        App_Log(APP_ERROR, "Fonction main: Incapable d'ouvrir la base de donnees: %s\n", sqlite3_errmsg(sqldbin));

        status = sqlite3_close(sqldbin);
        if( status != SQLITE_OK )
          App_Log(APP_ERROR,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

        return NOT_OK;
      } /* Fin du 'if ( status != SQLITE_OK )' */

      /* Execution de la requete SQL sur la base de donnees */
      /* Execution de la requete SQL sur la base de donnees */
      /* L'idee est de reproduire la commande UNIX
            echo .schema | sqlite3 obsin | sqlite3 obsout
         Cette requete provient de la documentation http://www.sqlite.org/faq.html#q7
      */
      strcpy(sqlschema,"");
      status = sqlite3_exec(sqldbin, "select * from sqlite_master", sqlite_schema_callback, sqlschema, &ErrMsg);
      if( status != SQLITE_OK ) {
        App_Log(APP_ERROR, "Fonction main: Erreur %d dans la fonction sqlite3_exec: %s\n", status, ErrMsg);
        sqlite3_free(ErrMsg);

        status = sqlite3_close(sqldbin);
        if( status != SQLITE_OK )
          App_Log(APP_ERROR,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

        return NOT_OK;
          } /* Fin du 'if ( status != SQLITE_OK )' */

      status = sqlite3_close(sqldbin);
      if( status != SQLITE_OK ) {
        App_Log(APP_ERROR,"Fonction main: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, opt.obsin);
        return NOT_OK;
      }

      status = sqlite3_open(opt.obsout,&sqldb);
      if ( status != SQLITE_OK ) {
        App_Log(APP_ERROR, "Fonction main: Incapable d'ouvrir la base de donnees pour le fichier '%s': %s\n", opt.obsout, sqlite3_errmsg(sqldb));

        status = sqlite3_close(sqldbin);
        if( status != SQLITE_OK )
          App_Log(APP_ERROR,"Fonction main: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, opt.obsin);
        return NOT_OK;
      } /* Fin du 'if ( status != SQLITE_OK )' */

      status = sqlite3_exec(sqldb, sqlschema, NULL, NULL, &ErrMsg);
      if( status != SQLITE_OK ){
        App_Log(APP_ERROR, "Fonction main: Erreur %d pour le fichier '%s' dans la fonction sqlite3_exec: %s\n", status, opt.obsout, ErrMsg);
        if (strcmp(ErrMsg,"PRIMARY KEY must be unique")==0) {
          App_Log(APP_ERROR,"Cette erreur est probablement due au fait que le fichier de sortie (%s) \n"
                  "n'a pas ete cree avant l'appel a ce programme avec l'utilitaire 'rdbgen'.  \n", opt.obsout);
        }
        sqlite3_free(ErrMsg);

        return NOT_OK;
      }
    } /* Fin du 'else' relie au 'if ( status == 0 )' */

    /* On cree la fonction checkgrid qui verifie si l'observation est a l'interieur
     * du domaine defini par la grille donnee plus haut.
     */
    status = sqlite3_create_function(sqldb, SQLFUNCTION_NAME, NUMBER_OF_ARGS_FOR_CHECK_GRID, SQLITE_UTF8,
                                     NULL, &checkgrid_sql, NULL, NULL);
    if( status != SQLITE_OK ) {
      App_Log(APP_ERROR,"Fonction main: Incapable de creer la fonction %s\n", SQLFUNCTION_NAME);

      status = sqlite3_close(sqldb);
      if( status != SQLITE_OK )
        App_Log(APP_ERROR,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

      return NOT_OK;
    }

    /* On cree la fonction checkvertical qui verifie si l'observation est a l'interieur
     * du domaine defini par la grille donnee plus haut.
     */
    status = sqlite3_create_function(sqldb, SQL_VERTICAL_NAME, NUMBER_OF_ARGS_FOR_CHECK_VERTICAL, SQLITE_UTF8,
                                       NULL, &checkvertical_sql, NULL, NULL);
    if( status != SQLITE_OK ) {
      App_Log(APP_ERROR, "Fonction main: Incapable de creer la fonction %s\n", SQL_VERTICAL_NAME);

      status = sqlite3_close(sqldb);
      if( status != SQLITE_OK )
        App_Log(APP_ERROR, "Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

      return NOT_OK;
    }

    /* On cree la fonction checkvertical_gz qui verifie si l'observation est a l'interieur
     * du domaine defini par la grille donnee plus haut.
     */
    status = sqlite3_create_function(sqldb, SQL_VERTICAL_GZ_NAME, NUMBER_OF_ARGS_FOR_CHECK_VERTICAL_GZ, SQLITE_UTF8,
                                     NULL, &checkvertical_gz_sql, NULL, NULL);
    if( status != SQLITE_OK ) {
      App_Log(APP_ERROR, "Fonction main: Incapable de creer la fonction %s\n", SQL_VERTICAL_NAME);

      status = sqlite3_close(sqldb);
      if( status != SQLITE_OK )
        App_Log(APP_ERROR, "Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

      return NOT_OK;
    }

    /* On cree la requete SQL a l'aide de l'information sur la grille que nous avons */
    if (strlen(opt.channels)==0 && opt.niveau_min == IP1_VIDE && opt.niveau_max == IP1_VIDE)
      /* Aucun filtrage vertical n'est fait */
      snprintf(requete_sql, sizeof(requete_sql), "attach '%s' as dbin; \n"
               "insert into %s select * from dbin.%s where %s(dbin.%s.lat,dbin.%s.lon,%d,%d,%d,%g,%g,%g,%g,%d,%d,%d,%d)=%d;\n"
                "insert into %s select * from dbin.%s where dbin.%s.%s in (select %s from %s);\n",
                opt.obsin, opt.rdb_header_table, opt.rdb_header_table, SQLFUNCTION_NAME, opt.rdb_header_table, opt.rdb_header_table,
                grid.gridid, grid.ni, grid.nj,
                opt.rect.min_i, opt.rect.max_i, opt.rect.min_j, opt.rect.max_j,
                opt.rect.min_i_equal, opt.rect.max_i_equal, opt.rect.min_j_equal, opt.rect.max_j_equal,
                opt.inout,
                opt.rdb_data_table, opt.rdb_data_table, opt.rdb_data_table, opt.rdb_split_on_key, opt.rdb_split_on_key, opt.rdb_header_table);
      else if (strlen(opt.channels)==0 && strlen(opt.gz)==0)
        /* Le filtrage vertical est fait a l'aide d'une hauteur en pression */
        snprintf(requete_sql, sizeof(requete_sql), "attach '%s' as dbin; \n"
                "insert into %s select * from dbin.%s where %s(dbin.%s.lat,dbin.%s.lon,%d,%d,%d,%g,%g,%g,%g,%d,%d,%d,%d)=%d;\n"
                "insert into %s select * from dbin.%s where dbin.%s.%s in (select %s from %s) and \n"
                "  %s(dbin.%s.%s,dbin.%s.vcoord,%d,%d)=1;\n",
                opt.obsin, opt.rdb_header_table, opt.rdb_header_table, SQLFUNCTION_NAME, opt.rdb_header_table, opt.rdb_header_table,
                grid.gridid, grid.ni, grid.nj,
                opt.rect.min_i, opt.rect.max_i, opt.rect.min_j, opt.rect.max_j,
                opt.rect.min_i_equal, opt.rect.max_i_equal, opt.rect.min_j_equal, opt.rect.max_j_equal,
                opt.inout,
                opt.rdb_data_table, opt.rdb_data_table, opt.rdb_data_table, opt.rdb_split_on_key, opt.rdb_split_on_key, opt.rdb_header_table,
                SQL_VERTICAL_NAME, opt.rdb_data_table, opt.rdb_split_on_key, opt.rdb_data_table, opt.niveau_min, opt.niveau_max);
      else if (strlen(opt.channels)==0)
        /* Le filtrage vertical est fait a l'aide d'une hauteur en metre */
        snprintf(requete_sql, sizeof(requete_sql), "attach '%s' as dbin; \n"
                "insert into %s select * from dbin.%s where %s(dbin.%s.lat,dbin.%s.lon,%d,%d,%d,%g,%g,%g,%g,%d,%d,%d,%d)=%d;\n"
                "insert into %s select %s.* from dbin.%s,%s where dbin.%s.%s = %s.%s and \n"
                "  %s(dbin.%s.%s,%s.lat,%s.lon,dbin.%s.vcoord+%s.elev,%d,%d,%d,%d,%d)=1;\n",
                opt.obsin, opt.rdb_header_table, opt.rdb_header_table, SQLFUNCTION_NAME, opt.rdb_header_table, opt.rdb_header_table,
                grid.gridid, grid.ni, grid.nj,
                opt.rect.min_i, opt.rect.max_i, opt.rect.min_j, opt.rect.max_j,
                opt.rect.min_i_equal, opt.rect.max_i_equal, opt.rect.min_j_equal, opt.rect.max_j_equal,
                opt.inout,
                opt.rdb_data_table, opt.rdb_data_table, opt.rdb_data_table, opt.rdb_header_table, opt.rdb_data_table, opt.rdb_split_on_key, opt.rdb_header_table, opt.rdb_split_on_key,
                SQL_VERTICAL_GZ_NAME, opt.rdb_data_table, opt.rdb_split_on_key, opt.rdb_header_table, opt.rdb_header_table, opt.rdb_data_table, opt.rdb_header_table,
                grid_gz.gridid, grid_gz.ni, grid_gz.nj, opt.niveau_min, opt.niveau_max);
      else if (opt.channels_voulus==1) /* On specifie plutot les canaux voulus  */
        snprintf(requete_sql, sizeof(requete_sql), "attach '%s' as dbin; \n"
                "insert into %s select * from dbin.%s where %s(dbin.%s.lat,dbin.%s.lon,%d,%d,%d,%g,%g,%g,%g,%d,%d,%d,%d)=%d;\n"
                "insert into %s select * from dbin.%s where dbin.%s.%s in (select %s from %s) and \n"
                "  dbin.%s.vcoord in (%s);\n",
                opt.obsin, opt.rdb_header_table, opt.rdb_header_table, SQLFUNCTION_NAME, opt.rdb_header_table, opt.rdb_header_table,
                grid.gridid, grid.ni, grid.nj,
                opt.rect.min_i, opt.rect.max_i, opt.rect.min_j, opt.rect.max_j,
                opt.rect.min_i_equal, opt.rect.max_i_equal, opt.rect.min_j_equal, opt.rect.max_j_equal,
                opt.inout,
                opt.rdb_data_table, opt.rdb_data_table, opt.rdb_data_table, opt.rdb_split_on_key, opt.rdb_split_on_key, opt.rdb_header_table,
                opt.rdb_data_table, opt.channels);
      else if (opt.channels_voulus==0) /* On specifie plutot les canaux exclus  */
        snprintf(requete_sql, sizeof(requete_sql), "attach '%s' as dbin; \n"
                "insert into %s select * from dbin.%s where %s(dbin.%s.lat,dbin.%s.lon,%d,%d,%d,%g,%g,%g,%g,%d,%d,%d,%d)=%d;\n"
                "insert into %s select * from dbin.%s where dbin.%s.%s in (select %s from %s) and \n"
                "  dbin.%s.vcoord not in (%s);\n",
                opt.obsin, opt.rdb_header_table, opt.rdb_header_table, SQLFUNCTION_NAME, opt.rdb_header_table, opt.rdb_header_table,
                grid.gridid, grid.ni, grid.nj,
                opt.rect.min_i, opt.rect.max_i, opt.rect.min_j, opt.rect.max_j,
                opt.rect.min_i_equal, opt.rect.max_i_equal, opt.rect.min_j_equal, opt.rect.max_j_equal,
                opt.inout,
                opt.rdb_data_table, opt.rdb_data_table, opt.rdb_data_table, opt.rdb_split_on_key, opt.rdb_split_on_key, opt.rdb_header_table,
                opt.rdb_data_table, opt.channels);
      else {
        App_Log(APP_ERROR, "Fonction main: Incapable de creer la requete SQL\n");

        status = sqlite3_close(sqldb);
        if( status != SQLITE_OK )
          App_Log(APP_ERROR, "Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

        return NOT_OK;
      }

      // On ajoute la requete SQL pour copier les tables qui ne contiennent par la 'split-on-key'.
      strcat(requete_sql,sqlreq_tables_without_split_key);

      append_table_list_split_key_requests_using_header(requete_sql,"dbin",table_list_with_split_key,opt.rdb_split_on_key,opt.rdb_header_table,opt.rdb_data_table);

      strcat(requete_sql,"\ndetach dbin;");

      printf("\n\nVoici la requete SQL effectuee sur la base de donnees:\n");
      printf("%s\n\n", requete_sql);

      /* Execution de la requete SQL sur la base de donnees finale */
      status = sqlite3_exec(sqldb, requete_sql, NULL, NULL, &ErrMsg);
      if( status != SQLITE_OK ){
        App_Log(APP_ERROR, "Fonction main: Erreur %d dans la fonction sqlite3_exec: %s\n", status, ErrMsg);
        if (strcmp(ErrMsg,"PRIMARY KEY must be unique")==0) {
          App_Log(APP_ERROR,"Cette erreur est probablement due au fait que le fichier de sortie (%s) \n"
                  "n'a pas ete cree avant l'appel a ce programme avec l'utilitaire 'rdbgen'.  \n", opt.obsout);
        }
        sqlite3_free(ErrMsg);
        EXIT_STATUS = OK;
      }

      /* On ferme la base de donnees ouverte plus haut avec sqlite3_open */
      status = sqlite3_close(sqldb);
      if( status != SQLITE_OK ) {
        App_Log(APP_ERROR,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);
        EXIT_STATUS = NOT_OK;
      }
    }
    else { /* Ici, opt.roundrobin == 1 */
      /* Il faut travailler les 'npex*npey' fichiers */
      char npex_str[MAXSTR], npey_str[MAXSTR], format_digits[MAXSTR], rdbout[MAXSTR*4];
      char sqlschema[MAXSTR*32];
      int nsplit = opt.npex*opt.npey;
      int ilonband, jlatband;

      strcpy(sqlschema,"");
      snprintf(format_digits, sizeof(format_digits), "%%.%dd",opt.ndigits);

      for (ilonband=0;ilonband<opt.npex;ilonband++) {
        for (jlatband=0;jlatband<opt.npey;jlatband++) {
          int id = ilonband*opt.npey+jlatband;

          /* Si on est en mode 'cherrypick', alors on ne considere que si la tuile est egale a celle voulue */
          if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0)
            if (ilonband != opt.cherrypick_x-1 || jlatband != opt.cherrypick_y-1)
              continue;

          snprintf(npex_str, sizeof(npex_str), format_digits,ilonband+1);
          snprintf(npey_str, sizeof(npey_str), format_digits,jlatband+1);
          snprintf(rdbout, sizeof(rdbout), "%s_%s_%s",opt.obsout,npex_str,npey_str);

          status = access(rdbout,F_OK);
          if ( status == 0 ) { /* Le fichier existe deja */
            /* On ouvre la base de donnees SQL de sortie */
            status = sqlite3_open(rdbout,&sqldb);
            if ( status != SQLITE_OK ) {
              App_Log(APP_ERROR, "Fonction main: Incapable d'ouvrir la base de donnees: %s\n", sqlite3_errmsg(sqldb));
              return NOT_OK;
            } /* Fin du 'if ( status != SQLITE_OK )' */
          } /* Fin du 'if ( status == 0 )' */
          else {
            if (strlen(sqlschema)==0) {
              /* Si le fichier n'existe pas encore alors on doit le creer et construire le meme schema */
              sqlite3 *sqldbin;

              /* On ouvre le fichier d'input */
              status = sqlite3_open(opt.obsin,&sqldbin);
              if ( status != SQLITE_OK ) {
                App_Log(APP_ERROR, "Fonction main: Incapable d'ouvrir la base de donnees: %s\n", sqlite3_errmsg(sqldbin));

                status = sqlite3_close(sqldbin);
                if( status != SQLITE_OK )
                  App_Log(APP_ERROR,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

                return NOT_OK;
              } /* Fin du 'if ( status != SQLITE_OK )' */

              /* Execution de la requete SQL sur la base de donnees */
              /* L'idee est de reproduire la commande UNIX
                 echo .schema | sqlite3 obsin | sqlite3 obsout
                 Cette requete provient de la documentation http://www.sqlite.org/faq.html#q7
              */
              status = sqlite3_exec(sqldbin, "select * from sqlite_master", sqlite_schema_callback, sqlschema, &ErrMsg);
              if( status != SQLITE_OK ) {
                App_Log(APP_ERROR, "Fonction main: Erreur %d pour le fichier '%s' dans la fonction sqlite3_exec: %s\n", status, opt.obsin, ErrMsg);
                sqlite3_free(ErrMsg);

                status = sqlite3_close(sqldbin);
                if( status != SQLITE_OK )
                  App_Log(APP_ERROR,"Fonction main: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, opt.obsin);

                return NOT_OK;
              } /* Fin du 'if ( status != SQLITE_OK )' */

              printf("Voici le schema du fichier d'input: '%s'\n%s\n", opt.obsin, sqlschema);

              status = sqlite3_close(sqldbin);
              if( status != SQLITE_OK ) {
                App_Log(APP_ERROR,"Fonction main: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, opt.obsin);

                return NOT_OK;
              }
            } /* Fin du 'if (strlen(sqlschema)==0)' */

            status = sqlite3_open(rdbout,&sqldb);
            if ( status != SQLITE_OK ) {
              App_Log(APP_ERROR, "Fonction main: Incapable d'ouvrir la base de donnees: %s\n", sqlite3_errmsg(sqldb));
              return NOT_OK;
            } /* Fin du 'if ( status != SQLITE_OK )' */

            status = sqlite3_exec(sqldb, sqlschema, NULL, NULL, &ErrMsg);
            if( status != SQLITE_OK ) {
              App_Log(APP_ERROR, "Fonction main: Erreur %d dans la fonction sqlite3_exec: %s\n", status, ErrMsg);
              if (strcmp(ErrMsg,"PRIMARY KEY must be unique")==0) {
                App_Log(APP_ERROR,"Cette erreur est probablement due au fait que le fichier de sortie (%s) \n"
                        "n'a pas ete cree avant l'appel a ce programme avec l'utilitaire 'rdbgen'.  \n", opt.obsout);
              }
              sqlite3_free(ErrMsg);
              return NOT_OK;
            }
          } /* Fin du 'else' relie au 'if ( status == 0 )' */

          /* On doit fabriquer la requete sql pour faire le splitting */
          snprintf(requete_sql, sizeof(requete_sql), "drop index if exists idx1;\n"
                              "PRAGMA journal_mode = OFF;\n"
                              "PRAGMA  synchronous = OFF;\n"
                              "attach '%s' as dbin; \n"
                              "%s\n", opt.obsin, sqlreq_tables_without_split_key);
          append_table_list_split_key_requests_nsplit(requete_sql,"dbin",table_list_with_split_key,opt.rdb_split_on_key,nsplit,id);
          if ( strcasecmp(opt.rdb_header_table,RDB_HEADER_DEFAUT) == 0   &&
               strcasecmp(opt.rdb_data_table,RDB_DATA_DEFAUT) == 0       &&
               strcasecmp(opt.rdb_split_on_key,RDB_SPLITONKEY_DEFAUT) == 0 ) {
            strcat(requete_sql,"create index idx1 on data(id_obs,vcoord,varno);\n");
          }
          strcat(requete_sql,"detach dbin;");

          printf("\nVoici la requete SQL effectuee sur la base de donnees pour creer le fichier '%s':\n",rdbout);
          printf("%s\n", requete_sql);

          /* Execution de la requete SQL sur la base de donnees finale */
          status = sqlite3_exec(sqldb, requete_sql, NULL, NULL, &ErrMsg);
          if( status != SQLITE_OK ) {
            App_Log(APP_ERROR, "Fonction main: Erreur %d dans la fonction sqlite3_exec: %s\n", status, ErrMsg);
            if (strcmp(ErrMsg,"PRIMARY KEY must be unique")==0) {
              App_Log(APP_ERROR,"Cette erreur est probablement due au fait que le fichier de sortie (%s) \n"
                      "n'a pas ete cree avant l'appel a ce programme avec l'utilitaire 'rdbgen'.  \n", opt.obsout);
            }
            sqlite3_free(ErrMsg);
            EXIT_STATUS = 1;
          }

          /* On ferme la base de donnees ouverte plus haut avec sqlite3_open */
          status = sqlite3_close(sqldb);
          if( status != SQLITE_OK ) {
            App_Log(APP_ERROR,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);
            EXIT_STATUS = 1;
          }
        } /* Fin du 'for (jlatband=0;jlatband<opt.npey;jlatband++)' */
      } /* Fin du 'for (ilonband=0;ilonband<opt.npex;ilonband++)' */
    } /* Fin du 'else' relie au 'if (opt.roundrobin == 0)' */

  return EXIT_STATUS;
} /* Fin de 'int splitobs_sql()' */


  /***************************************************************************
   * fonction: sqlite_schema_callback
   *
   * Cette fonction sert de 'callback' pour la requete executee par 'sqlite3_exec'.
   *    On y imprime les colonnes 'sql' pour obtenir le schema des tables.
   *
   ***************************************************************************/
static int sqlite_schema_callback(void *schema_void, int count, char **data, char **columns) {
  char* schema = schema_void;
  int idx;

  // printf("Dans la fonction 'sqlite_schema_callback'\n");
  // printf("There are %d column(s)\n", count);

  for (idx = 0; idx < count; idx++)
    if (strcmp(columns[idx],"sql")==0)
      if(strcmp(data[idx],"CREATE TABLE sqlite_stat1(tbl,idx,stat)")!=0) {
        // printf("The data in column \"%s\" is: '%s'\n", columns[idx], data[idx]);
        strcat(schema,data[idx]);
        strcat(schema,";\n");
      }

  // printf("\n");

  return 0;
}


  /***************************************************************************
   * fonction: sqlite_get_tables
   *
   * Trouve les tables qui contiennent une colonne 'split_on_key' et celles
   * qui ne contiennent pas la colonne 'split_on_key'.
   *
   ***************************************************************************/
int sqlite_get_tables(char* obsin, char* split_on_key, char* table_list_with_split_key, char* table_list_without_split_key) {
  int status;
  char *ErrMsg;
  sqlite3 *sqldbin;
  sqlite_get_tables_callback_arg callback_arg;

  /* Cette partie sert a trouver la requete pour copier les tables 'resume' et 'rdb4_schema' */
  /* On ouvre le fichier d'input */
  status = sqlite3_open(obsin,&sqldbin);
  if ( status != SQLITE_OK ) {
    App_Log(APP_ERROR, "Fonction sqlite_get_tables Incapable d'ouvrir le fichier '%s' avec l'erreur '%s'\n", obsin, sqlite3_errmsg(sqldbin));

    status = sqlite3_close(sqldbin);
    if( status != SQLITE_OK )
      App_Log(APP_ERROR,"Fonction sqlite_get_tables: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, obsin);

    return NOT_OK;
  } /* Fin du 'if ( status != SQLITE_OK )' */

  strcpy(table_list_with_split_key,"");
  strcpy(table_list_without_split_key,"");

  callback_arg.split_on_key = split_on_key;
  callback_arg.table_list_with_split_key = table_list_with_split_key;
  callback_arg.table_list_without_split_key = table_list_without_split_key;

  /* Execution de la requete SQL sur la base de donnees */
  status = sqlite3_exec(sqldbin, "select * from sqlite_master;", sqlite_get_tables_callback, &callback_arg, &ErrMsg);
  if( status != SQLITE_OK ) {
    App_Log(APP_ERROR, "Fonction sqlite_get_tables: Erreur %d pour le fichier dans la fonction sqlite3_exec: %s\n", status, ErrMsg);
    sqlite3_free(ErrMsg);

    return NOT_OK;
  } /* Fin du 'if ( status != SQLITE_OK )' */


  status = sqlite3_close(sqldbin);
  if( status != SQLITE_OK ) {
    App_Log(APP_ERROR,"Fonction sqlite_get_tables: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, obsin);

    return NOT_OK;
  }

  return OK;
} /* Fin de la fonction 'splitobs_obs' */


  /***************************************************************************
   * fonction: sqlite_get_tables_callback
   *
   * Cette fonction sert de 'callback' pour la requete executee par 'sqlite3_exec'.
   *    On veut savoir si elle contient les tables 'rdb4_schema' et 'resume'
   *
   ***************************************************************************/
static int sqlite_get_tables_callback(void *void_callback_arg, int count, char **data, char **columns) {
  int idx, isTypeTable, found_split_on_key;
  sqlite_get_tables_callback_arg *callback_arg;
  char table_name[MAXSTR];
  char* table_list;
  regex_t regex;
  char regex_errbuf[MAXSTR];
  int regex_err;

  callback_arg = (sqlite_get_tables_callback_arg*) void_callback_arg;

  isTypeTable = 0;
  for (idx = 0; idx < count; idx++) {
    if (strcasecmp(columns[idx],"type")==0)
      if (strcasecmp(data[idx],"table")==0)
        isTypeTable = 1;
  }

  if (isTypeTable == 0) return 0;

  // printf("BEGIN sqlite_get_tables_callback: split_on_key name: '%s'\n", callback_arg->split_on_key);

  regex_err = regcomp(&regex, callback_arg->split_on_key, REG_ICASE);
  if (regex_err!=0) {
    regerror(regex_err, &regex, regex_errbuf, MAXSTR);
    App_Log(APP_ERROR,"sqlite_get_tables_callback: cannot compile regular expression '%s': error '%s'",
            callback_arg->split_on_key, regex_errbuf);
    return 1;
  }

  strcpy(table_name,"");
  found_split_on_key=0;
  for (idx = 0; idx < count; idx++) {
    // printf("sqlite_get_tables_callback: columns[idx]='%s'\n", columns[idx]);
    // printf("sqlite_get_tables_callback: data[idx]='%s'\n", data[idx]);
    if (strcasecmp(columns[idx],"tbl_name")==0) {
      if(strcasecmp(data[idx],"sqlite_stat1")==0) {
        // printf("sqlite_get_tables_callback: ignore this table: '%s'\n", data[idx]);
        regfree(&regex);
        return 0;
      }
      strcpy(table_name,data[idx]);
    }
    else if(strcasecmp(columns[idx],"sql")==0) {
      // Ici, on regarde si 'id_obs' est contenu dans 'data[idx]' qui est de la forme:
      //       CREATE TABLE DATA (
      //       ID_DATA integer primary key   ,
      //       ID_OBS integer,
      //       BURP_BTYP integer,
      //       VCOORD integer,
      //       VARNO integer,
      //       VCOORD_TYPE integer,
      //       OBSVALUE real,
      //       BIAS_CORR real,
      //       SURF_EMISS real,
      //       CLOUD_EMISS real,
      //       FLAG integer,
      //       OMP real,
      //       OMA real,
      //       OBS_ERROR real,
      //       FG_ERROR real,
      //       TEMP_RAD_LOG10 real,
      //       CHAN_QC_FLAG integer,
      //       FSO real
      //       )

      regex_err = regexec(&regex,data[idx],0,(regmatch_t*) NULL,0);
      if (regex_err == 0) {
        // This means that there was a match
        found_split_on_key=1;
      }
      else if (regex_err == REG_NOMATCH) {
        found_split_on_key=0;
      }
      else {
        regerror(regex_err, &regex, regex_errbuf, MAXSTR);
        regfree(&regex);
        App_Log(APP_ERROR,"sqlite_get_tables_callback: Erreur '%s' avec la fonction regexec sur la ligne '%s'\n", regex_errbuf, data[idx]);
        return 1;
      }
    } // End of 'for (idx = 0; idx < count; idx++)'
  } // End of 'for (idx = 0; idx < count; idx++)'

  regfree(&regex);

  if (strlen(table_name)>0) {
    if (found_split_on_key) {
      table_list = callback_arg->table_list_with_split_key;
    }
    else {
      table_list = callback_arg->table_list_without_split_key;
    }

    if (strlen(table_list)>0) {
        strcat(table_list, " ");
    }
    strcat(table_list, table_name);
    /* printf("sqlite_get_tables_callback: found table: '%s'\n", callback_arg->table_list); */
  }
  else {
    App_Log(APP_ERROR,"sqlite_get_tables_callback: found_split_on_key = %d but table_name is empty", found_split_on_key);
    return 1;
  }

  // printf("END sqlite_get_tables_callback: split_on_key name: '%s'\n\n\n", callback_arg->split_on_key);

  return 0;
} /* End of function 'sqlite_get_tables_callback' */


  /***************************************************************************
   * fonction: append_table_list_split_key_requests_nsplit
   *
   * Cette fonction sert a ajouter des requetes SQL pour inclure les rangees des
   * tables qui ont 'split_on_key' comme colonne en se basant sur le critere
   *                 abs(split_on_key) % nsplit == id
   *
   ***************************************************************************/
void append_table_list_split_key_requests_nsplit(char* requete_sql, char* attached_db_name, char* table_list, char* split_on_key, int nsplit, int id) {
  const char separator_char[2] = " ";
  char sqlreqtmp[MAXSTR], table_list_tmp[MAXSTR];
  char *token;

  // Make a copy of 'table_list' input string because 'strtok' is changing in place that string
  strcpy(table_list_tmp,table_list);
  /* get the first token */
  token = strtok(table_list_tmp, separator_char);
  /* walk through other tokens */
  while( token != (char*) NULL ) {
    snprintf(sqlreqtmp, sizeof(sqlreqtmp), "insert into %s select * from %s.%s where abs(%s) %% %d = %d;\n",
            token,attached_db_name,token,split_on_key,nsplit,id);
    strcat(requete_sql,sqlreqtmp);
    token = strtok((char*) NULL, separator_char);
  }
} /* End of function 'append_table_list_split_key_requests_nsplit' */


  /***************************************************************************
   * fonction: append_table_list_split_key_requests_using_header
   *
   * Cette fonction sert a ajouter des requetes SQL pour inclure les rangees des
   * tables qui ont 'split_on_key' comme colonne et qui ont un 'split_on_key' qui
   * est dans le 'header_table'.
   *
   ***************************************************************************/
void append_table_list_split_key_requests_using_header(char* requete_sql, char* attached_db_name, char* table_list, char* split_on_key, char* header_table, char* data_table) {
  const char separator_char[2] = " ";
  char sqlreqtmp[MAXSTR], table_list_tmp[MAXSTR];
  char *token;

  // Make a copy of 'table_list' input string because 'strtok' is changing in place that string
  strcpy(table_list_tmp,table_list);
  /* get the first token */
  token = strtok(table_list_tmp, separator_char);
  /* walk through other tokens */
  while( token != (char*) NULL ) {
    // If the table is 'header' or 'data', then ignore it since they already have been considered in the request
    if ( strcasecmp(token,header_table) != 0 && strcasecmp(token,data_table) != 0) {
      snprintf(sqlreqtmp, sizeof(sqlreqtmp), "insert into %s select * from %s.%s where %s.%s.%s in (select %s from %s);\n",
              token,attached_db_name,token,attached_db_name,token,split_on_key,split_on_key,header_table);
      strcat(requete_sql,sqlreqtmp);
    }
    token = strtok((char*) NULL, separator_char);
  }
} /* End of function 'append_table_list_split_key_requests_using_header' */


  /***************************************************************************
   * fonction: append_table_list_without_split_key_requests
   *
   * Genere une requete SQL pour copier les tables qui n'ont pas la colonne 'split_on_key'.
   *
   ***************************************************************************/
void append_table_list_without_split_key_requests(char* requete_sql, char* attached_db_name, char* table_list) {
  const char separator_char[2] = " ";
  char sqlreqtmp[MAXSTR], table_list_tmp[MAXSTR];
  char *token;

  strcpy(requete_sql,"");
  // Make a copy of 'table_list' input string because 'strtok' is changing in place that string
  strcpy(table_list_tmp,table_list);
  /* get the first token */
  token = strtok(table_list_tmp, separator_char);
  /* walk through other tokens */
  while( token != (char*) NULL ) {
    snprintf(sqlreqtmp, sizeof(sqlreqtmp), "insert into %s select * from %s.%s;\n",token,attached_db_name,token);
    strcat(requete_sql,sqlreqtmp);
    token = strtok((char*) NULL, separator_char);
  }
} // End of function 'append_table_list_without_split_key_requests'

  /***************************************************************************
   * fonction: checkgrid_sql
   *
   * Cette fonction sert a verifier si le point (lat,lon) est a l'interieur d'une grille
   * donnee par les arguments d'entree:
   *     lat: coordonnee de latitude  du point d'observation
   *     lon: coordonnee de longitude du point d'observation
   *     gridid: identifiant de la grille EZSCINT
   *     ni: dimension horizontale de la grille
   *     nj: dimension verticale   de la grille
   *     min_i: 'i' minimal accepte
   *     max_i: 'i' maximal accepte
   *     min_j: 'j' minimal accepte
   *     max_j: 'j' maximal accepte
   ***************************************************************************/
void checkgrid_sql(sqlite3_context *context, int argc, sqlite3_value **argv) {
  int   gridid, ni, nj, status;
  float lat, lon;
  rectangle rect;
  char errmsg[MAXSTR];

  /* On s'assure que le nombre d'arguments est bien de 5 */
  assert( argc==NUMBER_OF_ARGS_FOR_CHECK_GRID );

  /* Lat-Lon de l'observation */
  lat = sqlite3_value_double(argv[0]);
  lon = sqlite3_value_double(argv[1]);
  /* Definition de la grille */
  gridid = sqlite3_value_int(argv[2]);
  ni     = sqlite3_value_int(argv[3]);
  nj     = sqlite3_value_int(argv[4]);
  /* Definition du rectangle */
  rect.min_i = sqlite3_value_double(argv[5]);
  rect.max_i = sqlite3_value_double(argv[6]);
  rect.min_j = sqlite3_value_double(argv[7]);
  rect.max_j = sqlite3_value_double(argv[8]);

  rect.min_i_equal = sqlite3_value_int(argv[ 9]);
  rect.max_i_equal = sqlite3_value_int(argv[10]);
  rect.min_j_equal = sqlite3_value_int(argv[11]);
  rect.max_j_equal = sqlite3_value_int(argv[12]);

  status = checkgrid(gridid, ni, nj, lat, lon, rect, errmsg);

  if (status<0) {
    sqlite3_result_error(context, errmsg, -1);
    return;
  }

  sqlite3_result_int(context, status);

} /* Fin de la fonction checkgrid_sql */


  /***************************************************************************
   * fonction: checkvertical_sql
   *
   * Cette fonction sert a verifier si le point (lat,lon,vcoord) est sous une certaine hauteur en hPa
   *  a partir d'une base de donnees SQL.  On appelle la fonction check_vertical().
   *     id_obs: identificateur de l'observation dans la table
   *     vcoord: hauteur de l'observation
   *     niveau_min: niveau minimum acceptable (en hPa)
   *     niveau_max: niveau maximum acceptable (en hPa)
   *
   ***************************************************************************/
void checkvertical_sql(sqlite3_context *context, int argc, sqlite3_value **argv) {
  int   id_obs, niveau_min, niveau_max, status;
  float vcoord;

  /* On s'assure que le nombre d'arguments est bien de 4 */
  assert( argc==NUMBER_OF_ARGS_FOR_CHECK_VERTICAL );

  /* identificateur de l'observations */
  id_obs = sqlite3_value_int(argv[0]);
  /* hauteur de l'observation */
  vcoord = sqlite3_value_double(argv[1]);
  /* Definition de la grille */
  niveau_min = sqlite3_value_int(argv[2]);
  niveau_max = sqlite3_value_int(argv[3]);

  if ( sqlite3_value_type(argv[1]) == SQLITE_NULL && niveau_min == IP1_VIDE) {
    if (SPLITOBS_SQLITE_VERBOSE>1) {
      printf("debug: id_obs=%d vcoord=NULL niveau_max=%d niveau_min=%d -> ",id_obs,niveau_max,niveau_min);
      printf("Obs acceptee parce qu'a la surface\n");
    }
    sqlite3_result_int(context, 1);
    return;
  }

  if (SPLITOBS_SQLITE_VERBOSE>1) {
    printf("debug: id_obs=%d vcoord=%f niveau_max=%d niveau_min=%d -> ",id_obs,vcoord,niveau_max,niveau_min);
  }

  status = checkvertical(vcoord,niveau_min,niveau_max);

  if (status<0) {
    char errmsg[MAXSTR];
    snprintf(errmsg, sizeof(errmsg),  "Fonction checkvertical_sql: Erreur avec checkvertical pour "
            "id_obs=%d vcoord=%f niveau_min=%d et niveau_max=%d\n", id_obs, vcoord, niveau_min, niveau_max);
    App_Log(APP_ERROR,"%s",errmsg);
    sqlite3_result_error(context, errmsg, -1);
    return;
  }

  sqlite3_result_int(context, status);
  return;

} /* Fin de la fonction checkvertical_sql */


  /***************************************************************************
   * fonction: checkvertical_gz_sql
   *
   * Cette fonction sert a verifier si le point (lat,lon,vcoord) est sous une certaine hauteur en hPa
   *  a partir d'une base de donnees SQL.  On appelle la fonction check_vertical().
   *     id_obs: identificateur de l'observation dans la table
   *     lat: coordonnee de latitude  du point d'observation
   *     lon: coordonnee de longitude du point d'observation
   *     vcoord: hauteur de l'observation
   *     gridid: identifiant de la grille EZSCINT
   *     ni: dimension horizontale de la grille
   *     nj: dimension verticale   de la grille
   *     niveau_min: niveau minimum acceptable (en hPa)
   *     niveau_max: niveau maximum acceptable (en hPa)
   *     gz: nom du fichier standard qui contient le champ GZ estimer la hauteur de la pression
   *
   ***************************************************************************/
void checkvertical_gz_sql(sqlite3_context *context, int argc, sqlite3_value **argv) {
  int   id_obs, gridid, ni, nj, niveau_min, niveau_max, status;
  float lat, lon, vcoord;

  /* On s'assure que le nombre d'arguments est bien de 9 */
  assert( argc==NUMBER_OF_ARGS_FOR_CHECK_VERTICAL_GZ );

  /* identificateur de l'observations */
  id_obs = sqlite3_value_int(argv[0]);
  /* Lat-Lon de l'observation */
  lat = sqlite3_value_double(argv[1]);
  lon = sqlite3_value_double(argv[2]);
  /* hauteur de l'observation */
  vcoord = sqlite3_value_double(argv[3]);
  /* Definition de la grille */
  gridid = sqlite3_value_int(argv[4]);
  ni     = sqlite3_value_int(argv[5]);
  nj     = sqlite3_value_int(argv[6]);
  /* niveaux seuil */
  niveau_min = sqlite3_value_int(argv[7]);
  niveau_max = sqlite3_value_int(argv[8]);

  if ( sqlite3_value_type(argv[3]) == SQLITE_NULL && niveau_min == IP1_VIDE) {
    if (SPLITOBS_SQLITE_VERBOSE>2) {
      printf("debug: id_obs=%d vcoord=NULL niveau_max=%d niveau_min=%d -> ",id_obs,niveau_max,niveau_min);
      printf("Obs acceptee parce qu'a la surface\n");
    }
    sqlite3_result_int(context, 1);
    return;
  }

  if (SPLITOBS_SQLITE_VERBOSE>2) {
    printf("debug: id_obs=%d lat=%f lon=%f vcoord=%f niveau_max=%d niveau_min=%d -> ",id_obs,lat,lon,vcoord,niveau_max,niveau_min);
  }

  status = checkvertical_gz(lat,lon,vcoord,gridid,ni,nj,niveau_min,niveau_max);

  if (status<0) {
    char errmsg[MAXSTR];
    snprintf(errmsg, sizeof(errmsg),  "Fonction checkvertical_gz_sql: Erreur avec checkvertical_gz pour "
             "id_obs=%d lat=%f lon=%f vcoord=%f niveau_min=%d niveau_max=%d\n",id_obs,lat,lon,vcoord,niveau_min,niveau_max);
    App_Log(APP_ERROR,"%s",errmsg);
    sqlite3_result_error(context, errmsg, -1);
    return;
  }

  sqlite3_result_int(context, status);
  return;

} /* Fin de la fonction checkvertical_gz_sql */

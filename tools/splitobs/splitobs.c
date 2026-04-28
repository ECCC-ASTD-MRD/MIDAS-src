#include <unistd.h> /* to get the function 'access' to know if files are already existing */

/* Include pour les librairies RPN */
#include "rmn.h"
#include "App.h"
#include "rmn/rpnmacros.h"

/* Include pour ma librairie de manipulation des fichiers standard RPN */
#include "fstdlib.h"
/* Include pour les constantes OK et NOT_OK */
#include "ok_or_notok.h"
/* Include pour la structure qui definit toutes les options */
#include "options.h"
/* Include pour les fonctions qui permettent de lire le champ GZ pour le clipping vertical */
#include "gzUtils.h"
/* Include pour la fonction 'splitobs_sql' */
#include "splitobs_sql.h"
/* Include pour la fonction 'splitobs_burp' */
#include "splitobs_burp.h"
/* Include pour la fonction 'splitobs_ascii' */
#include "splitobs_ascii.h"

/* Include qui permet d'obtenir la version a la compilation
 * (ce fichier est genere on-the-fly par le Makefile puis efface)
 */
#include "midas_build_info.h"

/********************************/
/*          main                */
/********************************/
int main(int argc, char** argv) {
  int status, EXIT_STATUS = OK, VERBOSE = 0;
  int filetype;
  gridtype grid, grid_gz;
  options  opt = optionsDEFAUT;
  /* Vecteurs qui contiennent les valeurs du champ GZ au niveau voulu
   * pour estimer la hauteur de la pression donnee
   */
  float* valeurs_gz_min = (float*) NULL;
  float* valeurs_gz_max = (float*) NULL;

  /**************************************************************/
  /* Impression de la boite indiquant le demarrage du programme */
  /**************************************************************/

  App_Init(APP_MAIN,PROJECT_NAME_STRING,VERSION,PROJECT_DESCRIPTION_STRING,BUILD_TIMESTAMP);
  App_LogStream("stdout");
  App_Start();

  /***************************************/
  /* On va lire les options au programme */
  /***************************************/
  status = parseOptions(argc,argv,&opt,&VERBOSE);
  if (status == NOT_OK) {
    App_Log(APP_ERROR,"Fonction main: probleme avec la fonction parseOptions\n");
    App_End(-1);exit(1);
  }

  status = access(opt.obsin,F_OK);
  if ( status != 0 ) { /* Le fichier n'existe pas  */
    App_Log(APP_ERROR,"The file '%s' does not exist and should.\n",opt.obsin);
    App_End(-1);exit(1);
  }

  /* Si on n'est pas en mode round-robin, alors on a besoin du fichier 'opt.fstin'. */
  if ( opt.roundrobin == 0 ) {
    /**********************************************************
     * Dans cette partie, on va lire le champ desire dans le fichier
     * opt.fstin tel qu'identifie par les elements dans la structure
     * (options) opt .
     **********************************************************/
    int file_handle = 1;
    status = getGridFromFile(&opt.fst, opt.fstin, &grid, &file_handle);
    if (status == NOT_OK) {
      App_Log(APP_ERROR,"Fonction main: Erreur dans la fonction getGridFromFile avec le fichier %s\n",opt.fstin);
      App_End(-1);exit(1);
    }

    set_domain_rectangle(&opt, grid);

    status = getMinMaxGZ(&opt, &grid_gz, &valeurs_gz_min, &valeurs_gz_max);
    if (status != OK) {
      /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'optptr->fstin'. */
      if ( opt.roundrobin == 0 ) {
        status = c_gdrls(grid.gridid);
        if (status<0)
          App_Log(APP_ERROR,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      }

      App_End(-1);exit(1);
    } /* Fin de 'if (status != OK)' */
  } /* Fin du 'if ( opt.roundrobin == 0 )' */

  filetype = c_wkoffit(opt.obsin,strlen(opt.obsin));

  if ( filetype == WKF_SQLITE3 ) {  /* Alors on traite une base de donnees SQL */
    status = splitobs_sql(opt, grid, grid_gz, valeurs_gz_min, valeurs_gz_max, VERBOSE);
    if ( status<0 ) {
      App_Log(APP_ERROR,"Fonction main: Erreur %d avec la fonction splitobs_sql\n", status);

      EXIT_STATUS = status;
    }
  } /* Fin du  if ( filetype == WKF_SQLITE3 ) */
  else if ( filetype == WKF_BURP ) {  /* Alors on traite un fichier BURP */
    status = splitobs_burp(opt, grid, grid_gz, valeurs_gz_min, valeurs_gz_max, VERBOSE);
    if ( status<0 ) {
      App_Log(APP_ERROR,"Fonction main: Erreur %d avec la fonction splitobs_burp\n", status);

      EXIT_STATUS = status;
    }
  } /* Fin du  if ( filetype == WKF_BURP ) */
  else if ( filetype == WKF_ASCII ) {  /* Alors on traite un fichier ASCII */
    status = splitobs_ascii(opt, grid, grid_gz, valeurs_gz_min, valeurs_gz_max, VERBOSE);
    if ( status<0 ) {
      App_Log(APP_ERROR,"Fonction main: Erreur %d avec la fonction splitobs_ascii\n", status);

      EXIT_STATUS = status;
    }
  } /* Fin du  if ( filetype == WKF_ASCII ) */

  /* Si on n'est pas en mode round-robin, alors on a besoin du fichier 'opt.fstin'. */
  if ( opt.roundrobin == 0 ) {
    /* On ferme la grille EZSCINT allouee pour definir la grille */
    status = c_gdrls(grid.gridid);
    if (status<0) {
      App_Log(APP_ERROR,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      EXIT_STATUS = 1;
    }

    if (strlen(opt.gz)!=0) {
      status = c_gdrls(grid_gz.gridid);
      if (status<0) {
        App_Log(APP_ERROR, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);
        EXIT_STATUS = NOT_OK;
      }

      if (valeurs_gz_min != (float*) NULL)
        free(valeurs_gz_min);
      if (valeurs_gz_max != (float*) NULL)
        free(valeurs_gz_max);
    }
  } /* Fin du 'if ( opt.roundrobin == 0 )'*/

  App_End(EXIT_STATUS);

  return(EXIT_STATUS);

} /* Fin de 'int main(int argc, char** argv)' */

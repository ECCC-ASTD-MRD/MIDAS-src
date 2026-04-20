#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h> /* POSIX says that <sys/types.h> must be included (by the caller) before <regex.h>.  */
#include <regex.h>
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
/* Include pour les fonctions qui permettent de voir si un point dans a l'interieur de la grille de definition */
#include "checkgrid.h"
/* Include pour la fonction 'splitobs_sql' */
#include "splitobs_sql.h"
/* Include pour la fonction 'splitobs_burp' */
#include "splitobs_burp.h"

/* Include qui permet d'obtenir la version a la compilation
 * (ce fichier est genere on-the-fly par le Makefile puis efface)
 */
#include "midas_build_info.h"

/* Prototypes of functions used in 'main' */
int  getGridFromFile(fstparam* fstptr, char* fichier, gridtype* gridptr, int* file_handle);
int  getGZ(char* fichier, gridtype* gridptr, int niveau, float** valeurs);
void set_domain_rectangle(optionsptr optptr, gridtype grid);
int  getMinMaxGZ(optionsptr optptr, gridtype* gridptr,
                 float** valeurs_gz_min_ptr, float** valeurs_gz_max_ptr);

/* Variable globale utilisee pour identifier le niveau de detail
 * que l'on veut dans le listing
 */
static int VERBOSE = 0;

/********************************/
/*          main                */
/********************************/
int main(int argc, char** argv) {
  int status, EXIT_STATUS = OK;
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
    status = splitobs_sql(opt, grid, grid_gz, VERBOSE, valeurs_gz_min, valeurs_gz_max);
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
    FILE* filein = (FILE*) NULL, *fileout = (FILE*) NULL;
    char ligne[MAXSTR], regex_errbuf[MAXSTR];
    char latstr[MAXSTR], lonstr[MAXSTR], altstr[MAXSTR];
    char errmsg[MAXSTR];
    float lat, lon, alt;
    int regex_err;
    regex_t regex;
    regmatch_t regex_match[5];

    if ( opt.roundrobin == 1 ) {
      App_Log(APP_ERROR,"Fonction main: Le mode 'round-robin' n'a pas encore ete "
              "implementee pour des fichiers d'input de type ASCII!\n");

      App_End(-1);exit(1);
    }

    filein=fopen(opt.obsin,"r");
    if (filein == (FILE*) NULL) {
      App_Log(APP_ERROR,"Fonction main: Le fichier ascii %s n'a pu etre ouvert correctement!\n", opt.obsin);

      status = c_gdrls(grid.gridid);
      if (status<0)
        App_Log(APP_ERROR,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

      App_End(-1);exit(1);
    }

    fileout=fopen(opt.obsout,"w");
    if (fileout == (FILE*) NULL) {
      App_Log(APP_ERROR,"Fonction main: Le fichier ascii %s n'a pu etre ouvert correctement!\n", opt.obsout);

      fclose(filein);

      status = c_gdrls(grid.gridid);
      if (status<0)
        App_Log(APP_ERROR,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

      if (strlen(opt.gz)!=0) {
        status = c_gdrls(grid_gz.gridid);
        if (status<0)
          App_Log(APP_ERROR, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

        if (valeurs_gz_min != (float*) NULL)
          free(valeurs_gz_min);
        if (valeurs_gz_max != (float*) NULL)
          free(valeurs_gz_max);
      }

      App_End(-1);exit(1);
    }

#define REGEX_DEFINITION  "^[[:blank:]]*([-]?[0-9]*\\.?[0-9]*)[[:blank:]]+([-]?[0-9]*\\.?[0-9]*)[[:blank:]]+([-]?[0-9]*\\.?[0-9]*).*$"

    regex_err = regcomp(&regex, REGEX_DEFINITION, REG_EXTENDED);
    if (regex_err!=0) {
      regerror(regex_err, &regex, regex_errbuf, MAXSTR);
      App_Log(APP_ERROR,"Fonction main: Erreur avec la fonction regcomp '%s' avec l'expression reguliere '%s'\n",regex_errbuf, REGEX_DEFINITION);

      fclose(filein);
      fclose(fileout);

      status = c_gdrls(grid.gridid);
      if (status<0)
        App_Log(APP_ERROR,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

      if (strlen(opt.gz)!=0) {
        status = c_gdrls(grid_gz.gridid);
        if (status<0)
          App_Log(APP_ERROR, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

        if (valeurs_gz_min != (float*) NULL)
          free(valeurs_gz_min);
        if (valeurs_gz_max != (float*) NULL)
          free(valeurs_gz_max);
      }

      App_End(-1);exit(1);
    }

    while( fgets(ligne,sizeof(ligne),filein) != (char*) NULL ) {

      regex_err = regexec(&regex,ligne,regex.re_nsub+1,regex_match,0);
      if (regex_err!=0) {
        regerror(regex_err, &regex, regex_errbuf, MAXSTR);
        App_Log(APP_ERROR,"Fonction main: Erreur avec la fonction regexec '%s' sur la ligne '%s'\n",regex_errbuf,ligne);
        EXIT_STATUS = 1;
        break;
      }

      strncpy(latstr,&ligne[regex_match[1].rm_so],regex_match[1].rm_eo-regex_match[1].rm_so);
      strncpy(lonstr,&ligne[regex_match[2].rm_so],regex_match[2].rm_eo-regex_match[2].rm_so);
      strncpy(altstr,&ligne[regex_match[3].rm_so],regex_match[3].rm_eo-regex_match[3].rm_so);

      lat = atof(latstr);
      lon = atof(lonstr);
      alt = atof(altstr);

      status = checkgrid(grid.gridid, grid.ni, grid.nj, lat, lon, opt.rect, errmsg, VERBOSE);
      if (status<0) {
        App_Log(APP_ERROR,"Fonction main: Erreur dans la fonction checkgrid pour la ligne '%s' avec le message '%s'\n", ligne, errmsg);
        EXIT_STATUS = 1;
        break;
      }

      if (opt.inout == status) {
        if (strlen(opt.channels)==0 && opt.niveau_min == IP1_VIDE && opt.niveau_max == IP1_VIDE)
          /* Aucun filtrage vertical n'est fait */
          fputs(ligne,fileout);
        else if (strlen(opt.channels)==0 && strlen(opt.gz)==0 && checkvertical(alt,opt.niveau_min,opt.niveau_max,VERBOSE))
          /* Le filtrage vertical est fait a l'aide d'une hauteur en pression */
          fputs(ligne,fileout);
        else if (strlen(opt.channels)==0 && checkvertical_gz(lat,lon,alt,grid_gz.gridid,grid_gz.ni,grid_gz.nj,opt.niveau_min,opt.niveau_max,valeurs_gz_min, valeurs_gz_max, VERBOSE))
          /* Le filtrage vertical est fait a l'aide d'une hauteur en metre */
          fputs(ligne,fileout);
        else if (opt.channels_voulus==checkcanal(alt,opt.channels,VERBOSE))
          fputs(ligne,fileout);
        else
          fputs("\n",fileout);
      }
      else
        fputs("\n",fileout);

    }

    regfree(&regex);
    fclose(filein);
    fclose(fileout);
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

/***************************************************************************
 * fonction: getGridFromFile
 *
 * Cette fonction sert a aller lire la grille horizontale d'un champ.
 *     fstptr: a pointeur vers une structure 'fstparam' qui contient l'information sur le champ a lire.
 *     fichier: nom du fichier standard dans lequel sera lu le champ GZ
 *     gridptr: pointeur a une structure de grille EZSCINT sur laquelle le champ GZ est definit

 *     file_handle: pointeur vers un 'int'.  Si sa valeur est 0,
 *                  alors on ne fermera pas le fichier et on donnera la valeur du
 *                  "file unit" pour que l'usager ferme lui-meme le fichier
 *
 ***************************************************************************/
int getGridFromFile(fstparam* fstptr, char* fichier, gridtype* gridptr, int* file_handle) {
  int iun, status;

  iun = 0;
  status = open_stdfile(iun, fichier, "RND+R/O");
  if (status == NOT_OK) {
    App_Log(APP_ERROR,"Fonction getGridFromFile: Erreur dans la fonction open_stdfile avec le fichier %s\n",fichier);
    App_End(-1);exit(1);
  }

  /* On va chercher la grille definie par le champ identifie avec la structure opt.fstin */
  status = getgrid(iun, gridptr, fstptr, fichier);
  if (status == NOT_OK) {
    App_Log(APP_ERROR,"Fonction getGridFromFile: Erreur dans la fonction getgrid pour les parametres "
            "(%s,%s,%s,%d,%d,%d,%d) dans le fichier %s\n",fstptr->nomvar,fstptr->typvar,fstptr->etiket,
            fstptr->dateo,fstptr->ip1,fstptr->ip2,fstptr->ip3,fichier);

    /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
    close_stdfile(iun, fichier);

    App_End(-1);exit(1);
  }

  if (*file_handle) {
    /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
    status = close_stdfile(iun, fichier);
    if (status == NOT_OK) {
      App_Log(APP_ERROR,"Fonction getGridFromFile: Erreur %d avec la fonction close_stdfile pour le fichier '%s'\n",status,fichier);

      status = c_gdrls(gridptr->gridid);
      if (status<0)
        App_Log(APP_ERROR,"Fonction getGridFromFile: Erreur dans la fonction c_gdrls pour gridid = %d\n", gridptr->gridid);

      return NOT_OK;
    }
  }
  else {
    *file_handle = iun;
  }

  return OK;
} /* Fin de 'int getGridFromFile(optionsptr optptr, char* fichier, gridtype* gridptr, int* file_handle)' */

/***************************************************************************
 * fonction: set_domain_rectangle(
 *
 * Cette fonction initialiser les valeurs de 'optptr->rect' en fonction du parametre 'optptr->pilot'
 *     optptr: un pointeur vers une structure 'options' qui sera mise a jour
 *     grid: information sur la grille qui est lue
 *
 ***************************************************************************/
void set_domain_rectangle(optionsptr optptr, gridtype grid) {
  /* Si on a defini la region avec l'option -pilot alors on definit le rectangle avec cette option */
  if (optptr->pilot!=PILOT_DEFAUT) {
    optptr->rect.min_i =       1+optptr->pilot;
    optptr->rect.max_i = grid.ni-optptr->pilot;
    optptr->rect.min_j =       1+optptr->pilot;
    optptr->rect.max_j = grid.nj-optptr->pilot;
  }
  else if ( optptr->npex != 1 || optptr->npey != 1) {
    /* Dans le cas d'une grille gaussienne, on se doit de considerer
     * la grille complete et d'ajouter 1 en longitude pour avoir les
     * points qui sont pres de la ligne de changement de date.
     */
    if (grid.grtyp[0]=='G') {
      if (optptr->rect.min_i<0) {
        optptr->rect.min_i=0;
        optptr->rect.min_j_equal=1;
      }
      if (optptr->rect.max_i<0) {
        optptr->rect.max_i=grid.ni+1;
        optptr->rect.max_j_equal=1;
      }
      if (optptr->rect.min_j<0) {
        optptr->rect.min_j=0;
        optptr->rect.min_i_equal=1;
      }
      if (optptr->rect.max_j<0) {
        optptr->rect.max_j=grid.nj+1;
        optptr->rect.max_i_equal=1;
      }
    } /* Fin du 'else if ( optptr->npex != 1 || optptr->npey != 1)' */
    else {
      if (optptr->rect.min_i<0)
        optptr->rect.min_i=1;
      if (optptr->rect.max_i<0) {
        optptr->rect.max_i=grid.ni;
        optptr->rect.max_i_equal=1;
      }
      if (optptr->rect.min_j<0) {
        optptr->rect.min_j=0;
        optptr->rect.min_j_equal=1;
      }
      if (optptr->rect.max_j<0)
        optptr->rect.max_j=grid.nj;
    } /* Fin du 'else' relie au if (grid.grtyp[0]=='G') */
  }
  else {
    if (optptr->rect.min_i<0)
      optptr->rect.min_i=1;
    if (optptr->rect.max_i<0) {
      optptr->rect.max_i=grid.ni;
      optptr->rect.max_i_equal=1;
    }
    if (optptr->rect.min_j<0) {
      optptr->rect.min_j=1;
      optptr->rect.min_j_equal=1;
    }
    if (optptr->rect.max_j<0)
      optptr->rect.max_j=grid.nj;
  } /* Fin du else relie au 'if ( optptr->npex != 1 || optptr->npey != 1)' */
} /* Fin de 'void set_domain_rectangle(optionsptr optptr, gridtype grid)' */

/***************************************************************************
 * fonction: getGZ
 *
 * Cette fonction sert a aller chercher le champ GZ dans le fichier donne
 *     optptr: un pointeur vers une structure 'options' qui sera mise a jour
 *     gridptr: pointeur a une structure de grille EZSCINT sur laquelle le champ GZ est definit
 *     valeurs_gz_min_ptr: pointeur de pointeur a un tableau de float pour stocker les valeurs minimales de GZ
 *     valeurs_gz_max_ptr: pointeur de pointeur a un tableau de float pour stocker les valeurs maximales de GZ
 *
 ***************************************************************************/
int getMinMaxGZ(optionsptr optptr, gridtype* gridptr,
                float** valeurs_gz_min_ptr, float** valeurs_gz_max_ptr) {
  int status;

  /**********************************************************
   * Dans cette partie, on va lire, si necessaire, le champ GZ dans le
   * fichier optptr->gz tel qu'identifie par les elements dans la
   * structure (options) opt .
   **********************************************************/
  if (strlen(optptr->gz)!=0) {
    if (optptr->rect.max_i>gridptr->ni) {
      printf("\nLe max_i donne en entree de %g est plus grand que la dimension de la grille ni=%d alors on met max_i a ni\n", optptr->rect.max_i, gridptr->ni);
      optptr->rect.max_i = gridptr->ni;
    }
    if (optptr->rect.max_j>gridptr->nj) {
      printf("\nLe max_j donne en entree de %g est plus grand que la dimension de la grille nj=%d alors on met max_j a nj\n", optptr->rect.max_j, gridptr->nj);
      optptr->rect.max_j = gridptr->nj;
    }
    if (optptr->rect.min_i<0) {
      printf("\nLe min_i donne en entree de %g est plus petit que 0 alors on met min_i a 0\n", optptr->rect.min_i);
      optptr->rect.min_i = 0;
    }
    if (optptr->rect.max_j<0) {
      printf("\nLe min_j donne en entree de %g est plus petit que 0 alors on met min_j a 0\n", optptr->rect.min_j);
      optptr->rect.min_j = 0;
    }

    if (optptr->niveau_min != IP1_VIDE) {
      status = getGZ(optptr->gz,gridptr,optptr->niveau_min,valeurs_gz_min_ptr);
      if (status != OK) {
        return NOT_OK;
      }
    }
    if (optptr->niveau_max != IP1_VIDE) {
      status = getGZ(optptr->gz,gridptr,optptr->niveau_max,valeurs_gz_max_ptr);
      if (status != OK) {
        if (*valeurs_gz_min_ptr != (float*) NULL)
          free(*valeurs_gz_min_ptr);

        return NOT_OK;
      }
    }
  } /* Fin du if (strlen(optptr->gz)!=0) */
  else { /* Si necessaire, on convertit les hPa en Pa en multipliant par 100 */
    if (optptr->niveau_min != IP1_VIDE) optptr->niveau_min *= 100;
    if (optptr->niveau_max != IP1_VIDE) optptr->niveau_max *= 100;
  }

  return OK;
} /* Fin de 'int getMinMaxGZ(int iun, char* fichier, gridtype* gridptr, int niveau, float** valeurs)' */

/***************************************************************************
 * fonction: getGZ
 *
 * Cette fonction sert a aller chercher le champ GZ dans le fichier donne
 *     iun: unite fortran qui sera utilisee pour lire dans le fichier
 *     fichier: nom du fichier standard dans lequel sera lu le champ GZ
 *     gridptr: pointeur a une structure de grille EZSCINT sur laquelle le champ GZ est definit
 *     valeurs: pointeur de pointeur a un tableau de float pour stocker les valeurs de GZ
 *
 ***************************************************************************/
int getGZ(char* fichier, gridtype* gridptr, int niveau, float** valeurs) {
  int iun, ier, key, status, datev;
  fstparam fst = fstparam_DEFAUT;
  double forecast;

  fst.ip1=niveau;
  strcpy(fst.nomvar,"GZ  ");

  status = getGridFromFile(&fst, fichier, gridptr, &iun);
  if (status == NOT_OK) {
    App_Log(APP_ERROR,"Fonction getGZ: Erreur dans la fonction getGridFromFile avec le fichier %s\n",fichier);
    App_End(-1);exit(1);
  }

  /**************************************************************
   * On lit maintenant les valeurs du champ GZ designe
   **************************************************************/
  /* Allocation de la memoire */
  *valeurs = (float*) malloc(gridptr->ni*gridptr->nj*sizeof(float));
  if ( *valeurs == (float*) NULL) {
    App_Log(APP_ERROR, "Fonction getGZ: Incapable d'allouer un vecteur de float de dimension %dx%d=%d elements\n",
            gridptr->ni,gridptr->nj,gridptr->ni*gridptr->nj);

    status = c_gdrls(gridptr->gridid);
    if (status<0)
      App_Log(APP_ERROR, "Fonction getGZ: Erreur dans la fonction c_gdrls pour gridid = %d\n", gridptr->gridid);

    /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
    close_stdfile(iun,fichier);

    return NOT_OK;
  }

  /* On doit incrementer l'heure de recherche pour la date de validite en fonction de fst.dateo et fst.npas, fst.deet */
  forecast = ((double) fst.npas*fst.deet)/3600; /* prevision en heures */
  incdatr_c(&datev,&fst.dateo,&forecast);

  key = c_fstinf(iun,&fst.ni,&fst.nj,&fst.nk,datev,fst.etiket,
                 fst.ip1,fst.ip2,fst.ip3,fst.typvar,fst.nomvar);
  if (key<0) {
    App_Log(APP_ERROR,"Fonction getGZ: Erreur %d avec le fichier %s pour les parametres (%s,%s,%s,%d,%d,%d,%d) dans la fonction c_fstinf\n",
            key,fichier,fst.nomvar,fst.typvar,fst.etiket,fst.dateo,fst.ip1,fst.ip2,fst.ip3);

    status = c_gdrls(gridptr->gridid);
    if (status<0)
      App_Log(APP_ERROR, "Fonction getGZ: Erreur dans la fonction c_gdrls pour gridid = %d\n", gridptr->gridid);

    /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
    close_stdfile(iun,fichier);

    return NOT_OK;
  }

  ier = c_fstluk(*valeurs,key,&gridptr->ni,&gridptr->nj,&gridptr->nk);
  if (ier<0) {
    App_Log(APP_ERROR, "Fonction getGZ: Erreur %d avec le fichier %s pour les parametres "
            "(%s,%s,%s,%d,%d,%d,%d) dans la fonction fstluk\n",
            ier,fichier,fst.nomvar,fst.typvar,fst.etiket,fst.dateo,fst.ip1,fst.ip2,fst.ip3);

    status = c_gdrls(gridptr->gridid);
    if (status<0)
      App_Log(APP_ERROR, "Fonction getGZ: Erreur dans la fonction c_gdrls pour gridid = %d\n", gridptr->gridid);

    /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
    close_stdfile(iun,fichier);

    free(*valeurs);

    return NOT_OK;
  }

  /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
  status = close_stdfile(iun,fichier);
  if (status == NOT_OK) {
    App_Log(APP_ERROR, "Fonction getGZ: Erreur avec la fonction close_stdfile pour le fichier %s\n",fichier);

    status = c_gdrls(gridptr->gridid);
    if (status<0)
      App_Log(APP_ERROR, "Fonction getGZ: Erreur dans la fonction c_gdrls pour gridid = %d\n", gridptr->gridid);

    free(*valeurs);

    return NOT_OK;
  }

  return OK;
} /* Fin de la fonction getGZ() */

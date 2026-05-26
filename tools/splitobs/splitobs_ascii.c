#include <stdio.h>
#include <sys/types.h> /* POSIX says that <sys/types.h> must be included (by the caller) before <regex.h>.  */
#include <regex.h>

/* Include pour les librairies RPN */
#include "rmn.h"
#include "App.h"

/* Include pour les constantes OK et NOT_OK */
#include "ok_or_notok.h"
/* Include pour ma librairie de manipulation des fichiers standard RPN */
#include "fstdlib.h"
/* Include pour la structure qui definit toutes les options */
#include "options.h"
/* Include pour les fonctions qui permettent de voir si un point dans a l'interieur de la grille de definition */
#include "checkgrid.h"

int splitobs_ascii(options opt, gridtype grid, gridtype grid_gz,
                   float* valeurs_gz_min, float* valeurs_gz_max, int VERBOSE) {
    FILE* filein = (FILE*) NULL, *fileout = (FILE*) NULL;
    char ligne[MAXSTR], regex_errbuf[MAXSTR];
    char latstr[MAXSTR], lonstr[MAXSTR], altstr[MAXSTR];
    char errmsg[MAXSTR];
    float lat, lon, alt;
    int regex_err, status, EXIT_STATUS;
    regex_t regex;
    regmatch_t regex_match[5];

    if ( opt.roundrobin == 1 ) {
      App_Log(APP_ERROR,"Fonction splitobs_ascii: Le mode 'round-robin' n'a pas encore ete "
              "implementee pour des fichiers d'input de type ASCII!\n");

      return NOT_OK;
    }

    filein=fopen(opt.obsin,"r");
    if (filein == (FILE*) NULL) {
      App_Log(APP_ERROR,"Fonction splitobs_ascii: Le fichier ascii %s n'a pu etre ouvert correctement!\n", opt.obsin);

      status = c_gdrls(grid.gridid);
      if (status<0)
        App_Log(APP_ERROR,"Fonction splitobs_ascii: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

      return NOT_OK;
    }

    fileout=fopen(opt.obsout,"w");
    if (fileout == (FILE*) NULL) {
      App_Log(APP_ERROR,"Fonction splitobs_ascii: Le fichier ascii %s n'a pu etre ouvert correctement!\n", opt.obsout);

      fclose(filein);

      return NOT_OK;
    }

#define REGEX_DEFINITION  "^[[:blank:]]*([-]?[0-9]*\\.?[0-9]*)[[:blank:]]+([-]?[0-9]*\\.?[0-9]*)[[:blank:]]+([-]?[0-9]*\\.?[0-9]*).*$"

    regex_err = regcomp(&regex, REGEX_DEFINITION, REG_EXTENDED);
    if (regex_err!=0) {
      regerror(regex_err, &regex, regex_errbuf, MAXSTR);
      App_Log(APP_ERROR,"Fonction splitobs_ascii: Erreur avec la fonction regcomp '%s' avec l'expression reguliere '%s'\n",regex_errbuf, REGEX_DEFINITION);

      fclose(filein);
      fclose(fileout);

      return NOT_OK;
    }

    EXIT_STATUS = OK;
    while( fgets(ligne,sizeof(ligne),filein) != (char*) NULL ) {

      regex_err = regexec(&regex,ligne,regex.re_nsub+1,regex_match,0);
      if (regex_err!=0) {
        regerror(regex_err, &regex, regex_errbuf, MAXSTR);
        App_Log(APP_ERROR,"Fonction splitobs_ascii: Erreur avec la fonction regexec '%s' sur la ligne '%s'\n",regex_errbuf,ligne);
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
        App_Log(APP_ERROR,"Fonction splitobs_ascii: Erreur dans la fonction checkgrid pour la ligne '%s' avec le message '%s'\n", ligne, errmsg);
        EXIT_STATUS = NOT_OK;
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
    } /* Fin de 'while( fgets(ligne,sizeof(ligne),filein) != (char*) NULL ) */

    regfree(&regex);
    fclose(filein);
    fclose(fileout);

    return EXIT_STATUS;
} /* End of function 'splitobs_ascii' */

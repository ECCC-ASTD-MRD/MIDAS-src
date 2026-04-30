#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h> /* to get the function 'access' to know if files are already existing */

#include "ok_or_notok.h"
/* Include pour ma librairie de manipulation des fichiers standard RPN */
#include "fstdlib.h"

/* Include pour les librairies RPN */
#include "rmn.h"

/*****************************************************/
/*******        Variables globales        ************/
/*****************************************************/

/* Ce compteur global est utilise dans la fonction stats_field() pour
 * compter le nombre de champs qui sera affiche.
 */
static int compteur_stats = 0;

/***************************************************************************
 * fonction: open_stdfile
 ***************************************************************************/
int open_stdfile(int* iun, char* filename, char* mode) {
  int ier;

  ier = c_fnom(iun,filename,mode,0);
  if (ier<0) {
    fprintf(stderr,"fonction open_stdfile: Erreur %d avec le fichier %s dans la fonction c_fnom\n",ier,filename);
    return NOT_OK;
  }
  ier = c_fstouv(*iun,"RND");
  if (ier<0) {
    fprintf(stderr,"fonction open_stdfile: Erreur %d avec le fichier %s dans la fonction c_fstouv\n",ier,filename);
    ier = c_fclos(*iun);
    if (ier<0)
      fprintf(stderr,"fonction open_stdfile: Erreur %d avec le fichier %s dans la fonction c_fclos\n",ier,filename);

    return NOT_OK;
  }

  return OK;
} /* Fin de la fonction open_stdfile */

/***************************************************************************
 * fonction: open_burpfile
 * This function is a copy of
 *      https://gitlab.science.gc.ca/wxobs-libs/libburpc/-/blob/master/src/burp_api.c?ref_type=heads
 * but where 'iun' is set to the unit file handle generated bu 'fnom'
 ***************************************************************************/
int open_burpfile(int* iun, char* filename, char* op) {
  int ier, number_of_reports;
  char type[30], mode[30];

  if (strcmp(op,"r") == 0) { /* READ ONLY */
    strcpy( type, "RND+OLD+R/O" );
    strcpy( mode, "READ" );
  }
  else if (strcmp(op,"w") == 0) { /* WRITE ONLY */
    strcpy( type, "RND+R/W" );
    strcpy( mode, "CREATE" );
  }
  else if (strcmp(op,"a") == 0) { /* APPEND */
    strcpy( type, "RND+APPEND" );
    strcpy( mode, "APPEND" );

    /* si fichier n'existe pas ouvrir en ecriture */
    ier = access(filename, F_OK);
    if (ier != 0) {
      strcpy( type, "RND+R/W" );
      strcpy( mode, "CREATE" );
    }
  } /* End of 'else if (strcmp(op,"a") == 0)' */
  else {
    fprintf(stderr,"fonction open_burpfile: Le mode '%s' n'est pas supporte pour ouvrir le fichier '%s'.\n",op,filename);

    return NOT_OK;
  }

  *iun = 0;
  ier = c_fnom(iun, filename, type, 0 );
  if ( ier < 0 ) {
    fprintf(stderr,"Unable to open file as %s : %s with 'c_fnom'\n",
            type,  filename);
    c_fclos(*iun);

    return NOT_OK;
  } /* Fin de 'if ( ier < 0 )' */

  /* The function 'c_mrfopn' returns the number of reports in the file */
  number_of_reports = c_mrfopn(*iun, mode);
  if ( ier < 0 ) {
    fprintf(stderr,"Unable to open file as %s : %s with 'c_mrfopn'\n",
            mode,  filename);

    ier = c_mrfcls(*iun);
    ier = c_fclos(*iun);

    return NOT_OK;
  } /* Fin de 'if ( ier < 0 )' */

  return number_of_reports;
} /* Fin de la fonction open_burpfile */


/***************************************************************************
 * fonction: close_stdfile
 ***************************************************************************/
int close_stdfile(int iun, char* filename) {
  int ier;

  /* On ferme le fichier standard ouvert par open_stdfile */
  ier = c_fstfrm(iun);
  if (ier<0) {
    fprintf(stderr,"fonction close_stdfile: Erreur %d avec le fichier %s dans la fonction c_fstfrm\n",ier,filename);
    return NOT_OK;
  }
  ier = c_fclos(iun);
  if (ier<0) {
    fprintf(stderr,"fonction close_stdfile: Erreur %d avec le fichier %s dans la fonction c_fclos\n",ier,filename);
    return NOT_OK;
  }

  return OK;
} /* Fin de la fonction close_stdfile */


/***************************************************************************
 * fonction: getgrid
 ***************************************************************************/
int getgrid(int iun, gridtype* gridptr, fstparam* fst, char* fstin) {
  int key, ier;

  /**********************************************************
   * Dans cette partie, on va lire le champ desire dans le fichier
   * opt.glb tel qu'identifie par les elements dans la structure
   * (options) opt .
   **********************************************************/

  key = c_fstinf(iun,&fst->ni,&fst->nj,&fst->nk,fst->dateo,fst->etiket,
                 fst->ip1,fst->ip2,fst->ip3,fst->typvar,fst->nomvar);
  if (key<0) {
    fprintf(stderr,"Fonction getgrid: Erreur %d avec le fichier %s pour les parametres (%s,%s,%s,%d,%d,%d,%d) dans la fonction c_fstinf\n",
            key,fstin,fst->nomvar,fst->typvar,fst->etiket,fst->dateo,fst->ip1,fst->ip2,fst->ip3);
    return NOT_OK;
  }
  ier = c_fstprm(key,&fst->dateo,&fst->deet,&fst->npas,&fst->ni,&fst->nj,&fst->nk,
                 &fst->nbits,&fst->datyp,&fst->ip1,&fst->ip2,&fst->ip3,fst->typvar,
                 fst->nomvar,fst->etiket,fst->grtyp,&fst->ig1,&fst->ig2,&fst->ig3,
                 &fst->ig4,&fst->swa,&fst->lng,&fst->dltf,&fst->ubc,&fst->extra1,
                 &fst->extra2,&fst->extra3);
  if (ier<0) {
    fprintf(stderr,"Fonction getgrid: Erreur %d avec le fichier %s pour les parametres (%s,%s,%s,%d,%d,%d,%d) dans la fonction c_fstprm\n",
            ier,fstin,fst->nomvar,fst->typvar,fst->etiket,fst->dateo,fst->ip1,fst->ip2,fst->ip3);
    return NOT_OK;
  }
  /**********************************************************
   * Les informations sur le champ meteo ont ete lues avec succes dans
   * le fichier fstin
   **********************************************************/

  gridptr->ni = fst->ni;
  gridptr->nj = fst->nj;

  gridptr->ig1 = fst->ig1;
  gridptr->ig2 = fst->ig2;
  gridptr->ig3 = fst->ig3;
  gridptr->ig4 = fst->ig4;
  strcpy(gridptr->grtyp,fst->grtyp);

  /**********************************************************
   * En ayant lu le fichier standard, on definit une grille EZSCINT
   * en utilisant les parametes de grille du champ en question
   **********************************************************/
  gridptr->gridid = c_ezqkdef(gridptr->ni,gridptr->nj,gridptr->grtyp,gridptr->ig1,gridptr->ig2,gridptr->ig3,gridptr->ig4,iun);
  if ( gridptr->gridid < 0 ) {
    fprintf(stderr,"Fonction getgrid: Probleme %d avec la fonction c_ezqkdef "
            "(nomvar=%s,etiket=%s,ip1=%d,ip2=%d,ip3=%d,date=%d,typvar=%s)\n",
            gridptr->gridid, fst->nomvar,fst->etiket,fst->ip1,fst->ip2,fst->ip3,fst->dateo,fst->typvar);
    return NOT_OK;
  }

  return OK;
} /* Fin de la fonction getgrid */


/***************************************************************************
 * fonction: stats_field
 ***************************************************************************/
int stats_field(float* z, int dim, fstparam* fstptr, statstype* statsptr) {
  int i;
  float moy = 0, var = 0;

  for (i=0;i<dim;i++) {
    moy += z[i];
    var += z[i]*z[i];
  }

  moy /= dim;
  var = var/dim - moy*moy;

  strcpy(statsptr[compteur_stats].nomvar,fstptr->nomvar);
  strcpy(statsptr[compteur_stats].etiket,fstptr->etiket);

  statsptr[compteur_stats].ip1  = fstptr->ip1;
  statsptr[compteur_stats].ip2  = fstptr->ip2;
  statsptr[compteur_stats].ip3  = fstptr->ip3;
  statsptr[compteur_stats].date = fstptr->dateo;

  statsptr[compteur_stats].moyenne  = moy;
  statsptr[compteur_stats].variance = var;

  compteur_stats++;

  return OK;

} /* Fin de la fonction stats_field */


/***************************************************************************
 * fonction: print_stats_field
 ***************************************************************************/
int print_stats_field(statstype* stats, int dim) {
  int i;

  printf("\n\nImpression des statistiques pour les champs\n");
  printf("nomvar\tip1\tip2\tip3\tdate\t\tetiket\t\tmoyenne\t\tvariance\n");
  for (i=0;i<dim;i++)
    printf("%s\t%d\t%d\t%d\t%d\t%s%E\t%E\n", stats[i].nomvar, stats[i].ip1, stats[i].ip2, stats[i].ip3, stats[i].date, stats[i].etiket, stats[i].moyenne, stats[i].variance);

  return OK;

} /* Fin de la fonction stats_field */


/***************************************************************************
 * fonction: padtime
 ***************************************************************************/
int padtime(char* argv) {
  int dat,tim,mode,datev;
  size_t j; // consistency with strlen()
  char dattimstr[MAXSTR], datstr[MAXSTR], timstr[MAXSTR];

  strcpy(dattimstr,argv);
  /* On padde avec des zeros a droite pour completer l'adresse */
  for(j=0;j<MAXSTR_DATETIME-strlen(argv);j++)
    strcat(dattimstr,"0");

  snprintf(datstr, sizeof(datstr), "%.*s", MAXSTR_DATE, dattimstr);
  snprintf(timstr, sizeof(timstr), "%.*s", MAXSTR_TIME, &dattimstr[MAXSTR_DATE]);

  dat = atoi(datstr);
  tim = atoi(timstr);

  mode = 3;
  newdate_c(&datev,&dat,&tim,&mode);  /* printable to CMCstamp */

  return datev;

} /* Fin de la fonction padtime */

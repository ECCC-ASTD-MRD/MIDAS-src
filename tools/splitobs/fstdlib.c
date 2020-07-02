#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ok_or_notok.h"
/* Include pour ma librairie de manipulation des fichiers standard RPN */
#include "fstdlib.h"

extern void f77name(newdate)();  /* rmnlib.a */
extern void f77name(exfin)();  /* rmnlib.a */
extern int c_fnom();
extern int c_fclos();
extern int c_fstouv();
extern int c_fstfrm();
extern int c_fstinf();
extern int c_fstprm();
extern int c_ezqkdef();

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
int open_stdfile(int iun, char* filename, char* mode) {
  int ier;

  ier = c_fnom(iun,filename,mode,0);
  if (ier<0) {
    fprintf(stderr,"fonction open_stdfile: Erreur %d avec le fichier %s dans la fonction c_fnom\n",ier,filename);
    return NOT_OK;
  }
  ier = c_fstouv(iun,"RND");
  if (ier<0) {
    fprintf(stderr,"fonction open_stdfile: Erreur %d avec le fichier %s dans la fonction c_fstouv\n",ier,filename);
    ier = c_fclos(iun);
    if (ier<0)
      fprintf(stderr,"fonction open_stdfile: Erreur %d avec le fichier %s dans la fonction c_fclos\n",ier,filename);

    return NOT_OK;
  }
  
  return OK;
} /* Fin de la fonction open_stdfile */


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
 * fonction: exit_program
 ***************************************************************************/
void exit_program(int status, char* program_name, char* problem, char* version) {

  F2Cl l_program_name = (F2Cl) strlen(program_name);
  F2Cl l_problem      = (F2Cl) strlen(problem);
  F2Cl l_version      = (F2Cl) strlen(version);

  if ( status == 1 ) {
    F2Cl l_fin = (F2Cl) strlen("FIN");
    f77name(exfin)(program_name,problem,"FIN",l_program_name,l_problem,l_fin);
    exit(1);
  }
  else {
    F2Cl l_ok = (F2Cl) strlen("O.K.");
    f77name(exfin)(program_name,version,"O.K.",l_program_name,l_version,l_ok);
  }
}  /* Fin de la fonction exit_program */

/***************************************************************************
 * fonction: padtime
 ***************************************************************************/
int padtime(char* argv) {
  int j,dat,tim,mode,datev;
  char dattimstr[MAXSTR], datstr[MAXSTR], timstr[MAXSTR];

  strcpy(dattimstr,argv);
  /* On padde avec des zeros a droite pour completer l'adresse */
  for(j=0;j<MAXSTR_DATETIME-strlen(argv);j++)
    strcat(dattimstr,"0");
  
  strncpy(datstr,dattimstr,MAXSTR_DATE);
  strncpy(timstr,&dattimstr[MAXSTR_DATE],MAXSTR_TIME);
  
  dat = atoi(datstr);
  tim = atoi(timstr);
  
  mode = 3;
  f77name(newdate)(&datev,&dat,&tim,&mode);  /* printable to CMCstamp */

  return datev;

} /* Fin de la fonction padtime */

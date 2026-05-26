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
/* Include pour les prototypes des fonctions definities dans ce fichier */
#include "gzUtils.h"

/***************************************************************************
 * fonction: set_domain_rectangle
 *
 * Cette fonction initialise les valeurs de 'optptr->rect' en fonction du parametre 'optptr->pilot'
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
    } /* Fin du 'if (grid.grtyp[0]=='G')' */
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
  } /* Fin du 'else if ( optptr->npex != 1 || optptr->npey != 1)' */
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
  } /* Fin du else relie au 'else if ( optptr->npex != 1 || optptr->npey != 1)' */
} /* Fin de 'void set_domain_rectangle(optionsptr optptr, gridtype grid)' */

/***************************************************************************
 * fonction: getGridFromFile
 *
 * Cette fonction sert a aller lire la grille horizontale d'un champ.
 *     fstptr: a pointeur vers une structure 'fstparam' qui contient l'information sur le champ a lire.
 *     fichier: nom du fichier standard dans lequel sera lu le champ GZ
 *     gridptr: pointeur a une structure de grille EZSCINT sur laquelle le champ GZ est definit
 *     file_handle: pointeur vers un 'int'.  Si sa valeur est 0, alors
 *                  on ne fermera pas le fichier et on donnera la
 *                  valeur du "file unit" pour que l'usager ferme
 *                  lui-meme le fichier
 *
 * Cette fonction retourne:
 *            OK si tout s'est bien passe
 *            NOT_OK s'il y a une erreur
 *
 ***************************************************************************/
int getGridFromFile(fstparam* fstptr, char* fichier, gridtype* gridptr, int* file_handle) {
  int iun, status;

  iun = 0;
  status = open_stdfile(&iun, fichier, "RND+R/O");
  if (status == NOT_OK) {
    App_Log(APP_ERROR,"Fonction getGridFromFile: Erreur dans la fonction open_stdfile avec le fichier %s\n",fichier);

    return NOT_OK;
  }

  /* On va chercher la grille definie par le champ identifie avec la structure opt.fstin */
  status = getgrid(iun, gridptr, fstptr, fichier);
  if (status == NOT_OK) {
    App_Log(APP_ERROR,"Fonction getGridFromFile: Erreur dans la fonction getgrid pour les parametres "
            "(%s,%s,%s,%d,%d,%d,%d) dans le fichier %s\n",fstptr->nomvar,fstptr->typvar,fstptr->etiket,
            fstptr->dateo,fstptr->ip1,fstptr->ip2,fstptr->ip3,fichier);

    /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
    close_stdfile(iun, fichier);

    return NOT_OK;
  }

  if (*file_handle != 0) {
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
} /* Fin de 'int getGridFromFile' */

/***************************************************************************
 * fonction: getMinMaxGZ
 *
 * Cette fonction sert a aller chercher le champ GZ dans le fichier donne
 *     optptr: un pointeur vers une structure 'options' qui sera mise a jour
 *     gridptr: pointeur a une structure de grille EZSCINT sur laquelle le champ GZ est definit
 *     valeurs_gz_min_ptr: pointeur de pointeur a un tableau de float pour stocker les valeurs minimales de GZ
 *     valeurs_gz_max_ptr: pointeur de pointeur a un tableau de float pour stocker les valeurs maximales de GZ
 *
 * Cette fonction retourne:
 *            OK si tout s'est bien passe
 *            NOT_OK s'il y a une erreur
 *
 ***************************************************************************/
int getMinMaxGZ(optionsptr optptr, gridtype* gridptr,
                float** valeurs_gz_min_ptr, float** valeurs_gz_max_ptr) {
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
      int status = getGZ(optptr->gz,gridptr,optptr->niveau_min,valeurs_gz_min_ptr);
      if (status != OK) {
        return NOT_OK;
      }
    }
    if (optptr->niveau_max != IP1_VIDE) {
      int status = getGZ(optptr->gz,gridptr,optptr->niveau_max,valeurs_gz_max_ptr);
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
 * Cette fonction retourne:
 *            OK si tout s'est bien passe
 *            NOT_OK s'il y a une erreur
 *
 ***************************************************************************/
int getGZ(char* fichier, gridtype* gridptr, int niveau, float** valeurs) {
  int iun, ier, key, status, datev;
  fstparam fst = fstparam_DEFAUT;
  double forecast;

  fst.ip1=niveau;
  strcpy(fst.nomvar,"GZ  ");

  iun = 0;
  status = getGridFromFile(&fst, fichier, gridptr, &iun);
  if (status == NOT_OK) {
    App_Log(APP_ERROR,"Fonction getGZ: Erreur dans la fonction getGridFromFile avec le fichier %s\n",fichier);

    return NOT_OK;
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

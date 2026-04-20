#include "checkgrid.h"

/* Include pour les librairies RPN */
#include "rmn.h"
#include "App.h"

/***************************************************************************
 * fonction: checkgrid
 *
 * Cette fonction sert a verifier si le point (lat,lon) est a l'interieur d'une grille
 * donnee par les arguments d'entree:
 *     gridid: identifiant de la grille EZSCINT
 *     ni: dimension horizontale de la grille
 *     nj: dimension verticale   de la grille
 *     lat: coordonnee de latitude  du point d'observation
 *     lon: coordonnee de longitude du point d'observation
 *     rect:  rectangle definissant une sous-region de la grille comme domaine
 *
 * Cette fonction retourne:
 *            1 si le point est a l'interieur de la grille
 *            0 si le point est a l'exterieur de la grille
 *           -1 s'il y a une erreur
 ***************************************************************************/
int checkgrid(int gridid, int ni, int nj, float lat, float lon, rectangle rect, char errmsg[MAXSTR], int VERBOSE) {
  int status, criteria;
  float x, y;

  /* appel a la fonction EZSCINT qui permet d'obtenir la coordonnee dans la grille
   * du point (lat,lon) donne en entree
   */
  if (VERBOSE>3)
    printf("Fonction checkgrid: lat=%f  lon=%f\n", lat, lon);

  if(lon<0) lon+=360.;

  status = c_gdxyfll(gridid, &x, &y, &lat, &lon, 1);
  if (status<0) {
    snprintf(errmsg, MAXSTR,  "Fonction checkgrid: Erreur avec c_gdxyfll qui retourne %d "
             "pour lat = %f, lon = %f, ni = %d, nj = %d, gridid = %d\n",
             status, lat, lon, gridid, ni, nj);
    App_Log(APP_ERROR,"%s",errmsg);
    return -1;
  }

  /* Si le point est effectivement dans le rectangle, on retourne 1 sinon 0 */
  if ( (!rect.min_i_equal && x>rect.min_i) || (rect.min_i_equal && x>=rect.min_i) ) {
    if ( (!rect.max_i_equal && x<rect.max_i) || (rect.max_i_equal && x<=rect.max_i) ) {
      if ( (!rect.min_j_equal && y>rect.min_j) || (rect.min_j_equal && y>=rect.min_j) ) {
        if ( (!rect.max_j_equal && y<rect.max_j) || (rect.max_j_equal && y<=rect.max_j) ) {
          if (VERBOSE>2)
            printf("Fonction checkgrid: Obs acceptee: lat = %f  lon = %f x = %f y = %f "
                   "rect.min_i = %f rect.max_i = %f rect.min_j = %f "
                   "rect.max_j = %f\n",lat,lon,x,y,rect.min_i,rect.max_i,rect.min_j,rect.max_j);
          return 1;
        }
        else
          criteria = 4;
      }
      else
        criteria = 3;
    }
    else
      criteria = 2;
  }
  else
    criteria = 1;

  if (VERBOSE>2) {
    printf("Fonction checkgrid: Obs refusee: lat = %f  lon = %f x = %f y= %f "
           "rect.min_i = %f rect.max_i = %f rect.min_j = %f "
           "rect.max_j = %f ",lat,lon,x,y,rect.min_i,rect.max_i,rect.min_j,rect.max_j);
    if (criteria==1)
      printf("(!rect.min_i_equal && x>rect.min_i) || (rect.min_i_equal && x>=rect.min_i)\n");
    else if (criteria==2) {
      printf("(!rect.max_i_equal && x<rect.max_i) || (rect.max_i_equal && x<=rect.max_i)\n");
      printf("(!rect.max_i_equal && x<rect.max_i) = %d       (rect.max_i_equal && x<=rect.max_i) = %d\n", (!rect.max_i_equal && x<rect.max_i), (rect.max_i_equal && x<=rect.max_i));
      printf("%f <= %f = %d max_i_equal = %d\n", x, rect.max_i, x<=rect.max_i, rect.max_i_equal);
    }
    else if (criteria==3)
      printf("(!rect.min_j_equal && y>rect.min_j) || (rect.min_j_equal && y>=rect.min_j)\n");
    else if (criteria==4)
      printf("(!rect.max_j_equal && y<rect.max_j) || (rect.max_j_equal && y<=rect.max_j)\n");
    else
      App_Log(APP_ERROR,"Fonction checkgrid: Le critere '%d' n'est pas possible.  Il faut 1,2, 3 ou 4.  ", criteria);

  }

  return 0;

} /* Fin de la fonction checkgrid */


  /***************************************************************************
   * fonction: checkvertical
   *
   * Cette fonction sert a verifier si le point (lat,lon,vcoord) est sous une certaine hauteur en hPa.
   *     vcoord: hauteur de l'observation
   *     niveau_min: niveau minimum acceptable (en hPa)
   *     niveau_max: niveau maximum acceptable (en hPa)
   *
   ***************************************************************************/
int checkvertical(float vcoord, int niveau_min, int niveau_max, int VERBOSE) {

  if (niveau_min == IP1_VIDE && niveau_max == IP1_VIDE) {
    /* Si aucun niveau n'a ete donne (==-1) alors on ne filtre pas verticalement
     * donc on retourne vrai (1)
     */
    if (VERBOSE>2) {
      printf("debug: Obs acceptee\n");
    }
    return 1;
  }
  else {
    if (VERBOSE>2) {
      printf("debug: vcoord=%f niveau_max=%d niveau_min=%d -> ",vcoord,niveau_max,niveau_min);
    }

    if (niveau_min != IP1_VIDE && niveau_max != IP1_VIDE) { /* On doit filtrer en haut et en bas */
      if (niveau_max <= vcoord && vcoord <= niveau_min) {
        if (VERBOSE>2) {
          printf("Obs acceptee parce que niveau_max=%d <= vcoord=%f et vcoord=%f <= niveau_min=%d\n", niveau_max, vcoord, vcoord, niveau_min);
        }
        return 1;
      }
      else {
        if (VERBOSE>2) {
          printf("Obs refusee parce que niveau_max=%d > vcoord=%f ou vcoord=%f > niveau_min=%d\n", niveau_max, vcoord, vcoord, niveau_min);
        }
        return 0;
      }
    } /* Fin du if (niveau_min != IP1_VIDE && niveau_max != IP1_VIDE) */
    else if (niveau_min != IP1_VIDE) { /* On doit filtrer par le bas */
      if (niveau_min >= vcoord) {
        if (VERBOSE>2) {
          printf("Obs acceptee parce que niveau_min=%d >= vcoord = %f\n", niveau_min, vcoord);
        }
        return 1;
      }
      else {
        if (VERBOSE>2) {
          printf("Obs refusee parce que niveau_min=%d < vcoord = %f\n", niveau_min, vcoord);
        }
        return 0;
      }
    } /* Fin du if (niveau_min != IP1_VIDE) */
    else if (niveau_max != IP1_VIDE) { /* On doit filtrer par le haut */
      if (vcoord >= niveau_max) {
        if (VERBOSE>2) {
          printf("Obs acceptee parce que vcoord = %f >= niveau_max=%d\n", vcoord, niveau_max);
        }
        return 1;
      }
      else {
        if (VERBOSE>2) {
          printf("Obs refusee parce que vcoord = %f < niveau_max=%d\n", vcoord, niveau_max);
        }
        return 0;
      }
    } /* Fin du if (niveau_max != IP1_VIDE) */
    else {
      App_Log(APP_ERROR, "Fonction checkvertical: Erreur pour niveau_min=%d et niveau_max=%d\n", niveau_min, niveau_max);
      return -1;
    }
  } /* Fin du else du if if (niveau_min == IP1_VIDE && niveau_max == IP1_VIDE ) */
} /* Fin de la fonction check_vertical */


  /***************************************************************************
   * fonction: checkvertical_gz
   *
   * Cette fonction sert a verifier si le point (lat,lon,vcoord) est sous une certaine hauteur en hPa
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
int checkvertical_gz(float lat, float lon, float vcoord, int gridid, int ni, int nj, int niveau_min, int niveau_max,
                     float* VALEURS_GZ_MIN, float* VALEURS_GZ_MAX, int VERBOSE) {
  int status;

  if (niveau_min == IP1_VIDE && niveau_max == IP1_VIDE) {
    /* Si aucun niveau n'a ete donne (==-1) alors on ne filtre pas verticalement
     * donc on retourne vrai (1)
     */
    if (VERBOSE>3) {
      printf("debug: Obs acceptee\n");
    }
    return 1;
  }
  else {
    if (VERBOSE>2) {
      printf("debug: lat=%f lon=%f vcoord=%f niveau_max=%d niveau_min=%d -> ",lat,lon,vcoord,niveau_max,niveau_min);
    }

    if (niveau_min != IP1_VIDE && niveau_max != IP1_VIDE) { /* On doit filtrer en haut et en bas */
      float hauteur_min, hauteur_max;

      /* appel a la fonction EZSCINT qui permet d'obtenir la valeur dans la grille
       * du point (lat,lon) donne en entree
       */
      status = c_gdllsval(gridid, &hauteur_min, VALEURS_GZ_MIN, &lat, &lon, 1);
      if (status<0) {
        char errmsg[MAXSTR];
        snprintf(errmsg, sizeof(errmsg),  "Fonction checkvertical_gz: c_gdllsval retourne %d "
                 "pour lat = %f, lon = %f, ni = %d, nj = %d, gridid = %d\n",
                 status, lat, lon, gridid, ni, nj);
        App_Log(APP_ERROR,"%s",errmsg);
        return -1;
      }

      status = c_gdllsval(gridid, &hauteur_max, VALEURS_GZ_MAX, &lat, &lon, 1);
      if (status<0) {
        char errmsg[MAXSTR];
        snprintf(errmsg, sizeof(errmsg),  "Fonction checkvertical_gz: c_gdllsval retourne %d "
                 "pour lat = %f, lon = %f, ni = %d, nj = %d, gridid = %d\n",
                 status, lat, lon, gridid, ni, nj);
        App_Log(APP_ERROR,"%s",errmsg);
        return -1;
      }

      /* On convertit le decametre du GZ en metres */
      hauteur_max *= 10;
      hauteur_min *= 10;

      /* Si le point est effectivement sous le niveau donne, on retourne 1 sinon 0 */
      if ( hauteur_min <= vcoord && vcoord <= hauteur_max ) {
        if (VERBOSE>2) {
          printf("Obs acceptee parce que hauteur_min=%f <= vcoord=%f et vcoord=%f <= hauteur_max=%f\n", hauteur_min, vcoord, vcoord, hauteur_max);
        }
        return 1;
      }
      else {
        if (VERBOSE>2) {
          printf("Obs refusee parce que hauteur_min=%f > vcoord=%f ou vcoord=%f > hauteur_max=%f\n", hauteur_min, vcoord, vcoord, hauteur_max);
        }
        return 0;
      }
    } /* Fin du if (niveau_min != IP1_VIDE && niveau_max != IP1_VIDE) */
    else if (niveau_min != IP1_VIDE) {
      float hauteur_min;

      /* appel a la fonction EZSCINT qui permet d'obtenir la valeur dans la grille
       * du point (lat,lon) donne en entree
       */
      status = c_gdllsval(gridid, &hauteur_min, VALEURS_GZ_MIN, &lat, &lon, 1);
      if (status<0) {
        char errmsg[MAXSTR];
        snprintf(errmsg, sizeof(errmsg),  "Fonction checkvertical_gz: c_gdllsval retourne %d "
                 "pour lat = %f, lon = %f, ni = %d, nj = %d, gridid = %d\n",
                 status, lat, lon, gridid, ni, nj);
        App_Log(APP_ERROR,"%s",errmsg);
        return -1;
      }

      /* On convertit le decametre du GZ en metres */
      hauteur_min *= 10;

      /* Si le point est effectivement sous le niveau donne, on retourne 1 sinon 0 */
      if ( hauteur_min <= vcoord ) {
        if (VERBOSE>2) {
          printf("Obs acceptee parce que hauteur_min=%f <= vcoord=%f\n", hauteur_min, vcoord);
        }
        return 1;
      }
      else {
        if (VERBOSE>2) {
          printf("Obs refusee parce que hauteur_min=%f > vcoord=%f\n", hauteur_min, vcoord);
        }
        return 0;
      }
    } /* Fin du if (niveau_min != IP1_VIDE) */
    else if (niveau_max != IP1_VIDE) { /* On doit filtrer en haut et en bas */
      float hauteur_max;

      /* appel a la fonction EZSCINT qui permet d'obtenir la valeur dans la grille
       * du point (lat,lon) donne en entree
       */

      status = c_gdllsval(gridid, &hauteur_max, VALEURS_GZ_MAX, &lat, &lon, 1);
      if (status<0) {
        char errmsg[MAXSTR];
        snprintf(errmsg, sizeof(errmsg),  "Fonction checkvertical_gz: c_gdllsval retourne %d "
                 "pour lat = %f, lon = %f, ni = %d, nj = %d, gridid = %d\n",
                 status, lat, lon, gridid, ni, nj);
        App_Log(APP_ERROR,"%s",errmsg);
        return -1;
      }

      /* On convertit le decametre du GZ en metres */
      hauteur_max *= 10;

      /* Si le point est effectivement sous le niveau donne, on retourne 1 sinon 0 */
      if (vcoord <= hauteur_max ) {
        if (VERBOSE>2) {
          printf("Obs acceptee parce que vcoord=%f <= hauteur_max=%f\n", vcoord, hauteur_max);
        }
        return 1;
      }
      else {
        if (VERBOSE>2) {
          printf("Obs refusee parce que vcoord=%f > hauteur_max=%f\n", vcoord, hauteur_max);
        }
        return 0;
      }
    } /* Fin du if (niveau_max != IP1_VIDE) */
    else {
      char errmsg[MAXSTR];
      snprintf(errmsg, sizeof(errmsg),  "Fonction checkvertical_gz: niveau_min=%d et niveau_max=%d\n", niveau_min, niveau_max);
      App_Log(APP_ERROR,"%s",errmsg);
      return -1;
    }
  } /* Fin du else du if (niveau_min == IP1_VIDE && niveau_max == IP1_VIDE ) */
} /* Fin de la fonction checkvertical_gz */


  /***************************************************************************
   * fonction: checkcanal
   *
   * Cette fonction sert a verifier si le canal "canal" est dans une liste séparée par des virgules
   *     canal: numero de canal que l'on veut verifier dans la liste "channels"
   *     channels: liste des canaux voulus telle que donnee avec l'option '-channels' ou '-nochannels'
   *
   ***************************************************************************/
int checkcanal(float canal, char* channels, int VERBOSE) {
  char* result;
  char  canalstr[MAXSTR];

  snprintf(canalstr, sizeof(canalstr), "%d", (int) canal);

  result = strstr(channels,canalstr);
  if (result != (char*) NULL) {
    if (VERBOSE>2) {
      printf("canal accepte parce que canal=%d est dans '%s'\n", (int) canal, channels);
    }
    return 1;
  }
  else {
    if (VERBOSE>2) {
      printf("canal refuse parce que canal=%d n'est pas dans '%s'\n", (int) canal, channels);
    }
    return 0;
  }

} /* Fin de la fonction checkcanal */

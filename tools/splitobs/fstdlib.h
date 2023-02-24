
#ifndef __INCLUDE_FSTDLIB_H__
#define __INCLUDE_FSTDLIB_H__

/* Include pour les librairies RPN */
#include "rmn/rpnmacros.h"

/* define pour les differentes longeurs de chaine de caracteres utilisees dans le programme */
#ifndef MAXSTR
#define MAXSTR            1024
#endif
#define MAXSTR_NOMVAR     5
#define MAXSTR_TYPVAR     3
#define MAXSTR_ETIKET     13
#define MAXSTR_GRTYP      2
#define MAXSTR_DATETIME   16
#define MAXSTR_DATE       8
#define MAXSTR_TIME       8   

/* Valeurs par defaut pour les elements de recherche des champ dans un fichier standard
 * avec la fonction RPN fstinf 
 */
#define NOMVAR_VIDE       "    "
#define TYPVAR_VIDE       "  "
#define ETIKET_VIDE       "            "
#define DATEV_VIDE        -1
#define IP1_VIDE          -1
#define IP2_VIDE          -1
#define IP3_VIDE          -1
#define GRTYP_VIDE        " "

/* valeurs par defaut pour certaines structures utilisees dans le programme */
#define fstparam_DEFAUT    {DATEV_VIDE,0,0,1,1,1,0,0,IP1_VIDE,IP2_VIDE,IP3_VIDE,0,0,0,0,0,0,0,0,0,0,0,TYPVAR_VIDE,NOMVAR_VIDE,ETIKET_VIDE,GRTYP_VIDE}

/* structure contenant les elements identifiant un champ dans un fichier standard */
typedef struct {
  int  dateo, deet, npas, ni, nj, nk, nbits, datyp, ip1, ip2, ip3;
  int  ig1, ig2, ig3, ig4, swa, lng, dltf, ubc, extra1, extra2, extra3;
  char typvar[MAXSTR_TYPVAR], nomvar[MAXSTR_NOMVAR], etiket[MAXSTR_ETIKET], grtyp[MAXSTR_GRTYP];
} fstparam;

/* structure contenant les elements definissant une grille avec EZSCINT */
typedef struct {
  int ni, nj, nk, gridid;
  char grtyp[MAXSTR_GRTYP];
  int ig1,ig2,ig3,ig4;
} gridtype;

/* structure contenant l'information pour l'affiche des statistiques des champs */
typedef struct {
  int date, ip1, ip2, ip3;
  char nomvar[MAXSTR_NOMVAR], etiket[MAXSTR_ETIKET];
  float moyenne, variance;
} statstype;


/***************************************************************************
 * fonction: open_stdfile
 * 
 * Cette fonction ouvre proprement le fichier standard filename selon un mode donne
 *
 * En entree, elle prend 3 arguments:
 *      iun: unite fortran identifiant le fichier lu
 *      filename: nom du fichier standard
 *      mode: mode l'ouverture du fichier (read, write, read et write...)
 *
 ***************************************************************************/
int    open_stdfile(int iun,  char* filename, char* mode);
/***************************************************************************
 * fonction: close_stdfile
 *
 * Cette fonction ferme proprement le fichier standard ouvert dans 
 * la fonction open_stdfile.  
 * 
 * En entree, elle prend 2 arguments:
 *      iun: unite fortran identifiant le fichier lu
 *      filename: nom du fichier standard
 *
 ***************************************************************************/
int    close_stdfile(int iun, char* filename);
/***************************************************************************
 * fonction: getgrid
 *
 * Cette fonction genere une representation EZSCINT d'une grille definie en
 * lisant un champ dans un fichier standard.  
 * 
 * En entree, elle prend 5 arguments:
 *      iun: unite fortran identifiant le fichier lu
 *      gridptr: pointeur a une structure de grille qui contiendra l'information sur la grille
 *      nomvar: NOMVAR du champ que l'on veut utiliser pour definir la grille
 *      fstin: nom du fichier standard utilise pour lire le champ definissant la grille
 * 
 ***************************************************************************/
int    getgrid(int iun, gridtype* gridptr, fstparam* fst, char* fstin);

/***************************************************************************
 * fonction: stats_field
 * 
 * Cette fonction imprime des statistiques sur le champ contenu dans "z"
 *
 * En entree, elle prend 4 arguments:
 *          z: un pointeur a un vecteur de float qui contient les valeurs du champ
 *          dim: dimension du champ
 *          fstptr: un pointeur a une structure fstparam qui identifie les parametres du champ
 *          statsptr: un pointeur a une structure statstype qui stockera l'information sur le champ
 * 
 ***************************************************************************/
int    stats_field(float* z, int dim, fstparam* fstptr, statstype* statsptr);
/***************************************************************************
 * fonction: print_stats_field
 * 
 * Cette fonction imprime les statistiques accumulees dans le vecteur
 *     
 * En entree, elle prend 2 arguments:
 *          stats: un pointeur a un vecteur de structures "statstype" contenant l'information statistiques sur tous les champs
 *          dim: dimension de ce vecteur
 * 
 ***************************************************************************/
int    print_stats_field(statstype* stats, int dim);


/***************************************************************************
 * fonction: exit_program
 *
 * Cette fonction imprime la boite indiquant la fin du programme dependant
 * si l'execution s'est deroulee correctement ou non.  
 *
 * En entree, elle prend 1 argument:
 *      status: entier etant egal a OK ou NOT_OK
 * 
 ***************************************************************************/
void exit_program(int status, char* program_name, char* problem, char* version);


/***************************************************************************
 * fonction: padtime
 *
 * Cette fonction convertit une date donnee dans le format YYYYMMDDHHMMSS
 * dans le standard CMCstamp
 * 
 * L'argument donne en entree est une chaine de caracteres dans le format YYYYMMDDHHMMSS
 * et on complete par des 0 a droite si la chaine est incomplete.  
 ***************************************************************************/
int padtime(char* argv);

#endif /* #ifndef __INCLUDE_FSTDLIB_H__ */

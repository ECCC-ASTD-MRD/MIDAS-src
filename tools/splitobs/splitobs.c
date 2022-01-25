#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/types.h> /* POSIX says that <sys/types.h> must be included (by the caller) before <regex.h>.  */
#include <regex.h>
#include <unistd.h> /* to get the function 'access' to know if files are already existing */

/* Include pour sqlite3 */
#include <sqlite3.h>

/* Include pour les librairies RPN */
#include "rpnmacros.h"
/* Include pour la librairie de manipulation de fichiers BURP*/
#include <burp_api.h>

/* Include pour ma librairie de manipulation des fichiers standard RPN */
#include "fstdlib.h"
/* Include pour les constantes OK et NOT_OK */
#include "ok_or_notok.h"
/* Include qui permet d'obtenir la version a la compilation
 * (ce fichier est genere on-the-fly par le Makefile puis efface)
 */
#include "version.h"

/* variables pour les fonction exdb et exfin */
#define PROGRAM_NAME      "splitobs"
#define PROBLEM           "Problem"

/* define pour les differentes longueurs de chaine de caracteres utilisees dans le programme */
#define MAXSTR            1024

/* nom des fonctions SQL qui seront creees pour rechercher les observations
 * dans la base de donnees
 */
#define SQLFUNCTION_NAME       "checkgrid_sql"
#define NUMBER_OF_ARGS_FOR_CHECK_GRID         13
#define SQL_VERTICAL_NAME      "checkvertical_sql"
#define NUMBER_OF_ARGS_FOR_CHECK_VERTICAL     4
#define SQL_VERTICAL_GZ_NAME   "checkvertical_gz_sql"
#define NUMBER_OF_ARGS_FOR_CHECK_VERTICAL_GZ  9

/* valeurs par defaut pour certaines structures utilisees dans le programme */
#define UNIT_NUMBER        200
#define NOMVAR_DEFAUT      "PN"
#define INOUT_DEFAUT       1
#define PILOT_DEFAUT       0
#define NPEX_DEFAUT        1
#define NPEY_DEFAUT        1
#define NDIGITS_DEFAUT     4
#define RECTANGLE_DEFAUT   {-1,-1,-1,-1,0,0,0,0}
#define RDB_HEADER_DEFAUT     "header"
#define RDB_DATA_DEFAUT       "data"
#define RDB_PRIMARYKEY_DEFAUT "header"
#define optionsDEFAUT      {"", /* fstin */ "", /* obsin*/ "", /* obsout*/ "", /* gz*/"", /* channels */ INOUT_DEFAUT, /* inout */ PILOT_DEFAUT, /* pilot */ RECTANGLE_DEFAUT, /* rect */ fstparam_DEFAUT, /* fst */ IP1_VIDE, /* niveau_min */ IP1_VIDE, /* niveau_max */ 1, /* channels_voulus */ NPEX_DEFAUT, /* npex */ NPEY_DEFAUT, /* npey */ NDIGITS_DEFAUT, /* ndigits */ 0, /* check_ua4d */ 0, /* roundrobin */ -1, /* cherrypick_x */ -1, /* cherrypick_y */ RDB_HEADER_DEFAUT, /* rdb_header_table */ RDB_DATA_DEFAUT, /* rdb_data_table */ RDB_PRIMARYKEY_DEFAUT} /* rdb_primarykey

/* differentes options du programme */
#define FSTIN_OPTION       "-fstin"
#define OBSIN_OPTION       "-obsin"
#define OBSOUT_OPTION      "-obsout"
#define RDBIN_OPTION       "-rdbin"   /* Remplacee par '-obsin'  */
#define RDBOUT_OPTION      "-rdbout"  /* Remplacee par '-obsout' */
#define BURPIN_OPTION      "-burpin"  /* Remplacee par '-obsin'  */
#define BURPOUT_OPTION     "-burpout" /* Remplacee par '-obsout' */
#define ASCII_OPTION       "-ascii"   /* Remplacee par '-obsin'  */
#define OUT_OPTION         "-out"     /* Remplacee par '-obsout' */
#define RDB_HEADER_OPTION      "-header"
#define RDB_DATA_OPTION        "-data"
#define RDB_PRIMARYKEY_OPTION  "-primarykey"
#define GZ_OPTION          "-gz"
#define NIVEAU_MIN_OPTION  "-niveau_min"
#define NIVEAU_MAX_OPTION  "-niveau_max"
#define CHANNELS_OPTION    "-channels"
#define NOCHANNELS_OPTION  "-nochannels"
#define INOUT_OPTION       "-inout"
#define PILOT_OPTION       "-pilot"
#define MIN_I_OPTION       "-min_i"
#define MAX_I_OPTION       "-max_i"
#define MIN_J_OPTION       "-min_j"
#define MAX_J_OPTION       "-max_j"
#define NPEX_OPTION        "-npex"    /* Separation du domaine en 'npex' portions egales selon 'i' */
#define NPEY_OPTION        "-npey"    /* Separation du domaine en 'npey' portions egales selon 'j' */
#define CHERRYPICK_X_OPTION "-x"
#define CHERRYPICK_Y_OPTION "-y"
#define NDIGITS_OPTION     "-ndigits"
#define NOMVAR_OPTION      "-nomvar"
#define TYPVAR_OPTION      "-typvar"
#define ETIKET_OPTION      "-etiket"
#define DATEV_OPTION       "-datev"
#define IP1_OPTION         "-ip1"
#define IP2_OPTION         "-ip2"
#define IP3_OPTION         "-ip3"
#define ROUNDROBIN_OPTION  "-round-robin"
#define VERBOSE_OPTION     "-verbose"
#define CHECK_UA4D_OPTION  "-check_ua4d"
#define HELP_OPTION1       "-h"
#define HELP_OPTION2       "-help"
#define HELP_OPTION3       "--help"

#define MAXSTR_CHANNELS    MAXSTR

/*************************************************************/
/******** Definition des structures pour le programme ********/
/*************************************************************/

typedef struct {
  float min_i, max_i;
  float min_j, max_j;
  int min_i_equal, max_i_equal;
  int min_j_equal, max_j_equal;
} rectangle;

/* structure contenant les informations provenant de l'appel au programme */
typedef struct {
  char fstin[MAXSTR], obsin[MAXSTR], obsout[MAXSTR];
  char gz[MAXSTR], channels[MAXSTR_CHANNELS];
  int  inout, pilot;
  rectangle    rect;
  fstparam fst;
  int  niveau_min, niveau_max, channels_voulus;
  int  npex, npey, ndigits;
  int check_ua4d;
  int roundrobin;
  int cherrypick_x, cherrypick_y;
  char rdb_header_table[MAXSTR], rdb_data_table[MAXSTR], rdb_primarykey[MAXSTR];
} options, *optionsptr;


/*****************************************************/
/******* Prototype des fonctions definies ************/
/*****************************************************/
static int sqlite_schema_callback(void *schema_void, int count, char **data, char **columns);
static int sqlite_check_resume_callback(void *is_resume_and_rdb4_schema_present_in_DB_void, int count, char **data, char **columns);
static int sqlite_check_tables_with_id_obs_callback(void *table_list, int count, char **data, char **columns);
void append_id_obs_table_list_requests(char* requete_sql, char* table_list);

int sqlite_add_resume_request(char* obsin, char* requete_sql, char* attached_db_name);
int sqlite_get_tables_with_id_obs(char* obsin, char* table_list);

int    getGZ(int iun, char* fichier, gridtype* gridptr, int niveau, float** valeurs);

void   checkgrid_sql(sqlite3_context *context, int argc, sqlite3_value **argv);
int    checkgrid(int gridid, int ni, int nj, float lat, float lon, rectangle rect, char errmsg[MAXSTR]);

void   checkvertical_sql(sqlite3_context *context, int argc, sqlite3_value **argv);
int    checkvertical(float vcoord, int niveau_min, int niveau_max);

void   checkvertical_gz_sql(sqlite3_context *context, int argc, sqlite3_value **argv);
int    checkvertical_gz(float lat, float lon, float vcoord, int gridid, int ni, int nj, int niveau_min, int niveau_max);

int    which_btyp(int btyp);
int    btypAssociated(int btyp_obs, int btyp);

int    clipping_vertical(BURP_RPT *rptin, optionsptr optptr, gridtype* grid_gz, BURP_RPT *rptout);

int    checkcanal(float canal, char* channels);
int    find_subdomain(int gridid, int ni, int nj, float lat, float lon, rectangle rect, int npex, int npey,
		      int* ilonband, int* jlatband, char errmsg[MAXSTR]);
int    putblk_nt(BURP_RPT *rpt, BURP_BLK *blk, int* t_in_domain, int nt);
int    putblk_nval(BURP_RPT *rpt, BURP_BLK *blk, int* val_in_domain, int nval);

int    extract_data_in_domains_along_nt(optionsptr optptr, gridtype* gridptr, BURP_RPT *rptin,
					int elem_lat, int elem_lon, int* nts, int** t_in_domain_ptr);
int    extract_data_in_domains_along_nval(optionsptr optptr, gridtype* gridptr, BURP_RPT *rptin,
					  int elem_lat, int elem_lon, BURP_BLK *blk_data,
					  int* nvals_in_domain, int* val_in_domain);
int    check_ua4d(BURP_RPT *rptin);
int    find_blk_data_in_rpt(BURP_RPT *rptin, int elem_lat, int elem_lon,
			    int** bknos_data_ptr, int** btyps_data_ptr, int* nblks);
int    find_blk_data_flag_in_rpt(BURP_RPT *rptin, int elem_lat, int elem_lon, int bkno_data,
				 BURP_BLK **blk_data_ptr, BURP_BLK **blk_flags_ptr, 
				 int* colonne_lat_ptr, int* colonne_lon_ptr);
int    fill_rptout_blk(BURP_RPT *rptin, BURP_RPT ** rptout, int* nts, int* t_in_domain,
                       int n, int cherrypick_x, int cherrypick_y, int npey);

int    parseOptions(int argc, char** argv, optionsptr optptr);
void   aide(void);

/* rmnlib.a */
extern void f77name(incdatr)();
extern void f77name(exdb)();
extern int c_gdrls(int);
extern int c_fstinf();
extern int c_fstluk();
extern int c_gdllsval();
extern int c_gdxyfll();
extern int c_mrfbfl(int);
extern int c_ezgprm();
extern wordint c_wkoffit();

/* Vecteur global qui contient les valeurs du champ GZ au niveau voulu
 * pour estimer la hauteur de la pression donnee
 */
float* VALEURS_GZ_MIN = (float*) NULL;
float* VALEURS_GZ_MAX = (float*) NULL;

/* Variable globale utilisee pour identifier le niveau de detail 
 * que l'on veut dans le listing 
 */
int    VERBOSE = 0;

/********************************/
/*          main                */
/********************************/
int f77name(splitobs)(int argc, char** argv) {
  int iun = UNIT_NUMBER, ier, status, EXIT_STATUS = 0, IS_INPUT_BURP_FILE;
  int filetype;
  char *ErrMsg, requete_sql[MAXSTR];
  gridtype grid, grid_gz;
  options  opt = optionsDEFAUT;
  sqlite3  *sqldb;

  /**************************************************************/
  /* Impression de la boite indiquant le demarrage du programme */
  /**************************************************************/

  F2Cl l_program_name, l_version, l_non;

  l_program_name = (F2Cl) strlen(PROGRAM_NAME);
  l_version      = (F2Cl) strlen(VERSION);
  l_non          = (F2Cl) strlen("NON");

  f77name(exdb)(PROGRAM_NAME,VERSION,"NON",l_program_name,l_version,l_non);

  printf("\n%s version: %s (SHA-1 %s)\n\n", PROGRAM_NAME, VERSION, VERSION_SHA1);

  /***************************************/
  /* on va lire les options au programme */
  /***************************************/
  status = parseOptions(argc,argv,&opt);
  if (status == NOT_OK) {
    fprintf(stderr, "Fonction main: probleme avec la fonction parseOptions\n");

    exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
  }

  /* Selon le code dans RMNLIB pour 'wkoffit.c', le chiffre qui sort de 'c_wkoffit' correspond a ceci: */
#define WKF_INEXISTANT            -3
#define WKF_VIDE                  -2
#define WKF_INCONNU               -1
#define WKF_RANDOM89               1
#define WKF_SEQUENTIEL89           2
#define WKF_SEQUENTIELFORTRAN89    3
#define WKF_CCRN                   4
#define WKF_CCRN_RPN               5
#define WKF_BURP                   6
#define WKF_GRIB                   7
#define WKF_BUFR                   8
#define WKF_BLOK                   9
#define WKF_FORTRAN               10
#define WKF_COMPRESS              11
#define WKF_GIF89                 12
#define WKF_GIF87                 13
#define WKF_IRIS                  14
#define WKF_JPG                   15
#define WKF_KMW                   16
#define WKF_PBM                   17
#define WKF_PCL                   18
#define WKF_PCX                   19
#define WKF_PDSVICAR              20
#define WKF_PM                    21
#define WKF_PPM                   22
#define WKF_PS                    23
#define WKF_KMW_                  24
#define WKF_RRBX                  25
#define WKF_SUNRAS                26
#define WKF_TIFF                  27
#define WKF_UTAHRLE               28
#define WKF_XBM                   29
#define WKF_XWD                   30
#define WKF_ASCII                 31
#define WKF_BMP                   32
#define WKF_RANDOM98              33
#define WKF_SEQUENTIEL98          34
#define WKF_NETCDF                35
/* Cette definition est un ajout */
#define WKF_SQLite                36

  status = access(opt.obsin,F_OK);
  if ( status != 0 ) { /* Le fichier existe deja */
      fprintf(stderr,"The file '%s' does not exist and should.\n",opt.obsin);
      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
  }

  filetype = c_wkoffit(opt.obsin,strlen(opt.obsin));
  /* Si le fichier n'est ni BURP ni ASCII alors on suppose que c'est un fichier SQLite */
  if (filetype != WKF_BURP) {
    /* Un fichier SQlite commence avec les caracteres 'SQLite format 3' */
#define SQLITE_FORMAT_STR_LEN     16 /* Longueur de 'SQLite format 3' +1 (pour le caractere nul a la fin de la string) */
    FILE* file;
    char line[SQLITE_FORMAT_STR_LEN];
    char* fgets_return_value;

    file = (FILE*) fopen(opt.obsin,"r");
    fgets_return_value = fgets(line,SQLITE_FORMAT_STR_LEN,file);
    fclose(file);
    if ( fgets_return_value != (char*) NULL ) {
      if (strcmp(line,"SQLite format 3") == 0)
        filetype = WKF_SQLite;
      else
        filetype = WKF_ASCII;
    }
    else {
      fprintf(stderr,"Cannot determine the file type since there was an error reading the file '%s'\n",opt.obsin);
      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }
  }

  /* Si on n'est pas en mode round-robin, alors on a besoin du fichier 'opt.fstin'. */
  if ( opt.roundrobin == 0 ) {
    /**********************************************************
     * Dans cette partie, on va lire le champ desire dans le fichier
     * opt.fstin tel qu'identifie par les elements dans la structure
     * (options) opt .
     **********************************************************/
    status = open_stdfile(iun, opt.fstin, "RND+R/O");
    if (status == NOT_OK) {
      fprintf(stderr, "Fonction main: Erreur dans la fonction open_stdfile avec le fichier %s\n",opt.fstin);

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

    /* On va chercher la grille definie par le champ identifie avec la structure opt.fstin */
    status = getgrid(iun,&grid,&opt.fst,opt.fstin);
    if (status == NOT_OK) {
      fprintf(stderr, "Fonction main: Erreur dans la fonction getgrid pour les parametres "
	      "(%s,%s,%s,%d,%d,%d,%d) dans le fichier %s\n",opt.fst.nomvar,opt.fst.typvar,opt.fst.etiket,
	      opt.fst.dateo,opt.fst.ip1,opt.fst.ip2,opt.fst.ip3,opt.fstin);

      /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
      close_stdfile(iun,opt.fstin);

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

    /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
    status = close_stdfile(iun,opt.fstin);
    if (status == NOT_OK) {
      fprintf(stderr, "Fonction main: Erreur %d avec la fonction close_stdfile pour le fichier '%s'\n",ier,opt.fstin);

      /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
      if ( opt.roundrobin == 0 ) {
	status = c_gdrls(grid.gridid);
	if (status<0)
	  fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      }

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

    /* Si on a defini la region avec l'option -pilot alors on definit le rectangle avec cette option */
    if (opt.pilot!=PILOT_DEFAUT) {
      opt.rect.min_i =       1+opt.pilot;
      opt.rect.max_i = grid.ni-opt.pilot;
      opt.rect.min_j =       1+opt.pilot;
      opt.rect.max_j = grid.nj-opt.pilot;
    }
    else if ( opt.npex != 1 || opt.npey != 1) {
      /* Dans le cas d'une grille gaussienne, on se doit de considerer
       * la grille complete et d'ajouter 1 en longitude pour avoir les
       * points qui sont pres de la ligne de changement de date.
       */
      if (grid.grtyp[0]=='G') {
	if (opt.rect.min_i<0) {
	  opt.rect.min_i=0;
	  opt.rect.min_j_equal=1;
	}
	if (opt.rect.max_i<0) {
	  opt.rect.max_i=grid.ni+1;
	  opt.rect.max_j_equal=1;
	}
	if (opt.rect.min_j<0) {
	  opt.rect.min_j=0;
	  opt.rect.min_i_equal=1;
	}
	if (opt.rect.max_j<0) {
	  opt.rect.max_j=grid.nj+1;
	  opt.rect.max_i_equal=1;
	}
      }
      else {
	if (opt.rect.min_i<0)
	  opt.rect.min_i=1;
	if (opt.rect.max_i<0) {
	  opt.rect.max_i=grid.ni;
	  opt.rect.max_i_equal=1;
	}
	if (opt.rect.min_j<0) {
	  opt.rect.min_j=0;
	  opt.rect.min_j_equal=1;
	}
	if (opt.rect.max_j<0)
	  opt.rect.max_j=grid.nj;
      } /* Fin du else relie au if (grid.grtyp[0]=='G') */
    }
    else {
      if (opt.rect.min_i<0)
	opt.rect.min_i=1;
      if (opt.rect.max_i<0) {
	opt.rect.max_i=grid.ni;
	opt.rect.max_i_equal=1;
      }
      if (opt.rect.min_j<0) {
	opt.rect.min_j=1;
	opt.rect.min_j_equal=1;
      }
      if (opt.rect.max_j<0)
	opt.rect.max_j=grid.nj;
    } /* Fin du else relie au 'if ( opt.npex != 1 || opt.npey != 1)' */

    /**********************************************************
     * Dans cette partie, on va lire, si necessaire, le champ GZ dans le
     * fichier opt.gz tel qu'identifie par les elements dans la
     * structure (options) opt .
     **********************************************************/
    if (strlen(opt.gz)!=0) {
      if (opt.rect.max_i>grid.ni) {
	printf("\nLe max_i donne en entree de %g est plus grand que la dimension de la grille ni=%d alors on met max_i a ni\n", opt.rect.max_i, grid.ni);
	opt.rect.max_i = grid.ni;
      }
      if (opt.rect.max_j>grid.nj) {
	printf("\nLe max_j donne en entree de %g est plus grand que la dimension de la grille nj=%d alors on met max_j a nj\n", opt.rect.max_j, grid.nj);
	opt.rect.max_j = grid.nj;
      }
      if (opt.rect.min_i<0) {
	printf("\nLe min_i donne en entree de %g est plus petit que 0 alors on met min_i a 0\n", opt.rect.min_i);
	opt.rect.min_i = 0;
      }
      if (opt.rect.max_j<0) {
	printf("\nLe min_j donne en entree de %g est plus petit que 0 alors on met min_j a 0\n", opt.rect.min_j);
	opt.rect.min_j = 0;
      }

      if (opt.niveau_min != IP1_VIDE) {
	status = getGZ(iun,opt.gz,&grid_gz,opt.niveau_min,&VALEURS_GZ_MIN);
	if (status == NOT_OK) {
	  /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
	  if ( opt.roundrobin == 0 ) {
	    status = c_gdrls(grid.gridid);
	    if (status<0)
	      fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
	  }

	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
      }
      if (opt.niveau_max != IP1_VIDE) {
	status = getGZ(iun,opt.gz,&grid_gz,opt.niveau_max,&VALEURS_GZ_MAX);
	if (status == NOT_OK) {
	  /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
	  if ( opt.roundrobin == 0 ) {
	    status = c_gdrls(grid.gridid);
	    if (status<0)
	      fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
	  }

	  if (VALEURS_GZ_MIN != (float*) NULL)
	    free(VALEURS_GZ_MIN);

	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
      }
    } /* Fin du if (strlen(opt.gz)!=0) */
    else { /* Si necessaire, on convertit les hPa en Pa en multipliant par 100 */
      if (opt.niveau_min != IP1_VIDE) opt.niveau_min *= 100;
      if (opt.niveau_max != IP1_VIDE) opt.niveau_max *= 100;
    }
  } /* Fin du 'if ( opt.roundrobin == 0 )' */

  if ( filetype == WKF_SQLite ) {  /* Alors on traite une base de donnees SQL */
    char sqlreq_resume[MAXSTR];
    char table_list[MAXSTR];

    /**********************************************************
     * Cette partie a pour but de manipuler la base de donnees SQL
     * et d'en creer une nouvelle qui ne contient que les observations a
     * l'interieur du domaine defini par la grille donnee plus haut.  
     **********************************************************/

    strcpy(table_list,"");
    status = sqlite_get_tables_with_id_obs(opt.obsin, table_list);
    if( status != OK ) {
      fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite_get_tables_with_id_obs pour le fichier '%s'\n", status, opt.obsin);

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

    status = sqlite_add_resume_request(opt.obsin,sqlreq_resume,"dbin");
    if( status != OK ) {
      fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite_add_resume_request pour le fichier '%s'\n", status, opt.obsin);

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

    /* Si on n'est pas en mode round-robin, alors les fichiers 'grid*.gridid' ont ete ouverts */
    if ( opt.roundrobin == 0 ) {
      /* On ouvre la base de donnees SQL finale */
      status = access(opt.obsout,F_OK);
      if ( status == 0 ) { /* Le fichier existe deja */
        status = sqlite3_open(opt.obsout,&sqldb);
        if ( status != SQLITE_OK ) {
          fprintf(stderr, "Fonction main: Incapable d'ouvrir la base de donnees: %s\n", sqlite3_errmsg(sqldb));

          status = sqlite3_close(sqldb);
          if( status != SQLITE_OK )
            fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

          status = c_gdrls(grid.gridid);
          if (status<0)
            fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

          if (strlen(opt.gz)!=0) {
            status = c_gdrls(grid_gz.gridid);
            if (status<0)
              fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

            if (VALEURS_GZ_MIN != (float*) NULL)
              free(VALEURS_GZ_MIN);
            if (VALEURS_GZ_MAX != (float*) NULL)
              free(VALEURS_GZ_MAX);
          }

          exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
        } /* Fin du 'if ( status != SQLITE_OK )' */
      } /* Fin du 'if ( status == 0 )' */
      else {
        /* Si le fichier n'existe pas encore alors on doit le creer et construire le meme schema */
        sqlite3 *sqldbin;
        char sqlschema[MAXSTR*32];

        /* On lit le fichier d'input */
        status = sqlite3_open(opt.obsin,&sqldbin);
        if ( status != SQLITE_OK ) {
          fprintf(stderr, "Fonction main: Incapable d'ouvrir la base de donnees: %s\n", sqlite3_errmsg(sqldbin));

          status = sqlite3_close(sqldbin);
          if( status != SQLITE_OK )
            fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

          status = c_gdrls(grid.gridid);
          if (status<0)
            fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

          if (strlen(opt.gz)!=0) {
            status = c_gdrls(grid_gz.gridid);
            if (status<0)
              fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

            if (VALEURS_GZ_MIN != (float*) NULL)
              free(VALEURS_GZ_MIN);
            if (VALEURS_GZ_MAX != (float*) NULL)
              free(VALEURS_GZ_MAX);
          }

          exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
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
          fprintf(stderr, "Fonction main: Erreur %d dans la fonction sqlite3_exec: %s\n", status, ErrMsg);
          sqlite3_free(ErrMsg);

          status = sqlite3_close(sqldbin);
          if( status != SQLITE_OK )
            fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

          status = c_gdrls(grid.gridid);
          if (status<0)
            fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

          if (strlen(opt.gz)!=0) {
            status = c_gdrls(grid_gz.gridid);
            if (status<0)
              fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

            if (VALEURS_GZ_MIN != (float*) NULL)
              free(VALEURS_GZ_MIN);
            if (VALEURS_GZ_MAX != (float*) NULL)
              free(VALEURS_GZ_MAX);
          }

          exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
        } /* Fin du 'if ( status != SQLITE_OK )' */

        status = sqlite3_close(sqldbin);
        if( status != SQLITE_OK ) {
          fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, opt.obsin);

          status = c_gdrls(grid.gridid);
          if (status<0)
            fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

          if (strlen(opt.gz)!=0) {
            status = c_gdrls(grid_gz.gridid);
            if (status<0)
              fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

            if (VALEURS_GZ_MIN != (float*) NULL)
              free(VALEURS_GZ_MIN);
            if (VALEURS_GZ_MAX != (float*) NULL)
              free(VALEURS_GZ_MAX);
          }

          exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
        }

        status = sqlite3_open(opt.obsout,&sqldb);
        if ( status != SQLITE_OK ) {
          fprintf(stderr, "Fonction main: Incapable d'ouvrir la base de donnees pour le fichier '%s': %s\n", opt.obsout, sqlite3_errmsg(sqldb));

          status = sqlite3_close(sqldbin);
          if( status != SQLITE_OK )
            fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, opt.obsin);

          status = c_gdrls(grid.gridid);
          if (status<0)
            fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

          if (strlen(opt.gz)!=0) {
            status = c_gdrls(grid_gz.gridid);
            if (status<0)

              fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

            if (VALEURS_GZ_MIN != (float*) NULL)
              free(VALEURS_GZ_MIN);
            if (VALEURS_GZ_MAX != (float*) NULL)
              free(VALEURS_GZ_MAX);
          }

          exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
        } /* Fin du 'if ( status != SQLITE_OK )' */

        status = sqlite3_exec(sqldb, sqlschema, (void*) NULL, (void*) NULL, &ErrMsg);
        if( status != SQLITE_OK ){
          fprintf(stderr, "Fonction main: Erreur %d pour le fichier '%s' dans la fonction sqlite3_exec: %s\n", status, opt.obsout, ErrMsg);
          if (strcmp(ErrMsg,"PRIMARY KEY must be unique")==0) {
            fprintf(stderr,"Cette erreur est probablement due au fait que le fichier de sortie (%s) \n"
                    "n'a pas ete cree avant l'appel a ce programme avec l'utilitaire 'rdbgen'.  \n", opt.obsout);
          }
          sqlite3_free(ErrMsg);

          status = c_gdrls(grid.gridid);
          if (status<0)
            fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

          if (strlen(opt.gz)!=0) {
            status = c_gdrls(grid_gz.gridid);
            if (status<0)

              fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

            if (VALEURS_GZ_MIN != (float*) NULL)
              free(VALEURS_GZ_MIN);
            if (VALEURS_GZ_MAX != (float*) NULL)
              free(VALEURS_GZ_MAX);
          }

          exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
        }
      } /* Fin du 'else' relie au 'if ( status == 0 )' */

      /* On cree la fonction checkgrid qui verifie si l'observation est a l'interieur
       * du domaine defini par la grille donnee plus haut.
       */
      status = sqlite3_create_function(sqldb, SQLFUNCTION_NAME, NUMBER_OF_ARGS_FOR_CHECK_GRID, SQLITE_UTF8,
                                       (void*) NULL, &checkgrid_sql, (void*) NULL, (void*) NULL);
      if( status != SQLITE_OK ) {
        fprintf(stderr,"Fonction main: Incapable de creer la fonction %s\n", SQLFUNCTION_NAME);

        status = sqlite3_close(sqldb);
        if( status != SQLITE_OK )
          fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

        status = c_gdrls(grid.gridid);
        if (status<0)
          fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

        if (strlen(opt.gz)!=0) {
          status = c_gdrls(grid_gz.gridid);
          if (status<0)
            fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

          if (VALEURS_GZ_MIN != (float*) NULL)
            free(VALEURS_GZ_MIN);
          if (VALEURS_GZ_MAX != (float*) NULL)
            free(VALEURS_GZ_MAX);
        }

        exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
      }

      /* On cree la fonction checkvertical qui verifie si l'observation est a l'interieur
       * du domaine defini par la grille donnee plus haut.
       */
      status = sqlite3_create_function(sqldb, SQL_VERTICAL_NAME, NUMBER_OF_ARGS_FOR_CHECK_VERTICAL, SQLITE_UTF8,
                                       (void*) NULL, &checkvertical_sql, (void*) NULL, (void*) NULL);
      if( status != SQLITE_OK ) {
        fprintf(stderr, "Fonction main: Incapable de creer la fonction %s\n", SQL_VERTICAL_NAME);

        status = sqlite3_close(sqldb);
        if( status != SQLITE_OK )
          fprintf(stderr, "Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

        status = c_gdrls(grid.gridid);
        if (status<0)
          fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

        if (strlen(opt.gz)!=0) {
          status = c_gdrls(grid_gz.gridid);
          if (status<0)
            fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

          if (VALEURS_GZ_MIN != (float*) NULL)
            free(VALEURS_GZ_MIN);
          if (VALEURS_GZ_MAX != (float*) NULL)
            free(VALEURS_GZ_MAX);
        }

        exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
      }

      /* On cree la fonction checkvertical_gz qui verifie si l'observation est a l'interieur
       * du domaine defini par la grille donnee plus haut.
       */
      status = sqlite3_create_function(sqldb, SQL_VERTICAL_GZ_NAME, NUMBER_OF_ARGS_FOR_CHECK_VERTICAL_GZ, SQLITE_UTF8,
                                       (void*) NULL, &checkvertical_gz_sql, (void*) NULL, (void*) NULL);
      if( status != SQLITE_OK ) {
        fprintf(stderr, "Fonction main: Incapable de creer la fonction %s\n", SQL_VERTICAL_NAME);

        status = sqlite3_close(sqldb);
        if( status != SQLITE_OK )
          fprintf(stderr, "Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

        status = c_gdrls(grid.gridid);
        if (status<0)
          fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

        if (strlen(opt.gz)!=0) {
          status = c_gdrls(grid_gz.gridid);
          if (status<0)
            fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

          if (VALEURS_GZ_MIN != (float*) NULL)
            free(VALEURS_GZ_MIN);
          if (VALEURS_GZ_MAX != (float*) NULL)
            free(VALEURS_GZ_MAX);
        }

        exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
      }

      /* On cree la requete SQL a l'aide de l'information sur la grille que nous avons */
      if (strlen(opt.channels)==0 && opt.niveau_min == IP1_VIDE && opt.niveau_max == IP1_VIDE)
        /* Aucun filtrage vertical n'est fait */
        sprintf(requete_sql,"attach '%s' as dbin; \n"
                "insert into header select * from dbin.header where %s(dbin.header.lat,dbin.header.lon,%d,%d,%d,%g,%g,%g,%g,%d,%d,%d,%d)=%d;\n"
                "insert into data   select * from dbin.data   where dbin.data.id_obs in (select id_obs from header);",
                opt.obsin, SQLFUNCTION_NAME, grid.gridid, grid.ni, grid.nj,
                opt.rect.min_i, opt.rect.max_i, opt.rect.min_j, opt.rect.max_j,
                opt.rect.min_i_equal, opt.rect.max_i_equal, opt.rect.min_j_equal, opt.rect.max_j_equal,
                opt.inout);
      else if (strlen(opt.channels)==0 && strlen(opt.gz)==0)
        /* Le filtrage vertical est fait a l'aide d'une hauteur en pression */
        sprintf(requete_sql,"attach '%s' as dbin; \n"
                "insert into header select * from dbin.header where %s(dbin.header.lat,dbin.header.lon,%d,%d,%d,%g,%g,%g,%g,%d,%d,%d,%d)=%d;\n"
                "insert into data   select * from dbin.data   where dbin.data.id_obs in (select id_obs from header) and \n"
                "  %s(dbin.data.id_obs,dbin.data.vcoord,%d,%d)=1;",
                opt.obsin, SQLFUNCTION_NAME, grid.gridid, grid.ni, grid.nj,
                opt.rect.min_i, opt.rect.max_i, opt.rect.min_j, opt.rect.max_j,
                opt.rect.min_i_equal, opt.rect.max_i_equal, opt.rect.min_j_equal, opt.rect.max_j_equal,
                opt.inout, SQL_VERTICAL_NAME, opt.niveau_min, opt.niveau_max);
      else if (strlen(opt.channels)==0)
        /* Le filtrage vertical est fait a l'aide d'une hauteur en metre */
        sprintf(requete_sql,"attach '%s' as dbin; \n"
                "insert into header select * from dbin.header where %s(dbin.header.lat,dbin.header.lon,%d,%d,%d,%g,%g,%g,%g,%d,%d,%d,%d)=%d;\n"
                "insert into data   select data.* from dbin.data,header where dbin.data.id_obs = header.id_obs and \n"
                "  %s(dbin.data.id_obs,header.lat,header.lon,dbin.data.vcoord+header.elev,%d,%d,%d,%d,%d)=1;",
                opt.obsin, SQLFUNCTION_NAME, grid.gridid, grid.ni, grid.nj,
                opt.rect.min_i, opt.rect.max_i, opt.rect.min_j, opt.rect.max_j,
                opt.rect.min_i_equal, opt.rect.max_i_equal, opt.rect.min_j_equal, opt.rect.max_j_equal,
                opt.inout, SQL_VERTICAL_GZ_NAME, grid_gz.gridid, grid_gz.ni,
                grid_gz.nj, opt.niveau_min, opt.niveau_max);
      else if (opt.channels_voulus==1) /* On specifie plutot les canaux voulus  */
        sprintf(requete_sql,"attach '%s' as dbin; \n"
                "insert into header select * from dbin.header where %s(dbin.header.lat,dbin.header.lon,%d,%d,%d,%g,%g,%g,%g,%d,%d,%d,%d)=%d;\n"
                "insert into data   select * from dbin.data   where dbin.data.id_obs in (select id_obs from header) and \n"
                "  dbin.data.vcoord in (%s);",
                opt.obsin, SQLFUNCTION_NAME, grid.gridid, grid.ni, grid.nj,
                opt.rect.min_i, opt.rect.max_i, opt.rect.min_j, opt.rect.max_j,
                opt.rect.min_i_equal, opt.rect.max_i_equal, opt.rect.min_j_equal, opt.rect.max_j_equal,
                opt.inout, opt.channels);
      else if (opt.channels_voulus==0) /* On specifie plutot les canaux exclus  */
        sprintf(requete_sql,"attach '%s' as dbin; \n"
                "insert into header select * from dbin.header where %s(dbin.header.lat,dbin.header.lon,%d,%d,%d,%g,%g,%g,%g,%d,%d,%d,%d)=%d;\n"
                "insert into data   select * from dbin.data   where dbin.data.id_obs in (select id_obs from header) and \n"
                "  dbin.data.vcoord not in (%s);",
                opt.obsin, SQLFUNCTION_NAME, grid.gridid, grid.ni, grid.nj,
                opt.rect.min_i, opt.rect.max_i, opt.rect.min_j, opt.rect.max_j,
                opt.rect.min_i_equal, opt.rect.max_i_equal, opt.rect.min_j_equal, opt.rect.max_j_equal,
                opt.inout, opt.channels);
      else {
        fprintf(stderr, "Fonction main: Incapable de creer la requete SQL\n");

        status = sqlite3_close(sqldb);
        if( status != SQLITE_OK )
          fprintf(stderr, "Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

        status = c_gdrls(grid.gridid);
        if (status<0)
          fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

        if (strlen(opt.gz)!=0) {
          status = c_gdrls(grid_gz.gridid);
          if (status<0)
            fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

          if (VALEURS_GZ_MIN != (float*) NULL)
            free(VALEURS_GZ_MIN);
          if (VALEURS_GZ_MAX != (float*) NULL)
            free(VALEURS_GZ_MAX);
        }

        exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
      }

      // On ajoute la requete SQL pour copier les tables 'rdb4_schema' et 'resume'.
      strcat(requete_sql,sqlreq_resume);

      append_id_obs_table_list_requests(requete_sql,table_list);

      strcat(requete_sql,"\ndetach dbin;");

      printf("\n\nVoici la requete SQL effectuee sur la base de donnees:\n");
      printf("%s\n\n", requete_sql);

      /* Execution de la requete SQL sur la base de donnees finale */
      status = sqlite3_exec(sqldb, requete_sql, (void*) NULL, (void*) NULL, &ErrMsg);
      if( status != SQLITE_OK ){
        fprintf(stderr, "Fonction main: Erreur %d dans la fonction sqlite3_exec: %s\n", status, ErrMsg);
        if (strcmp(ErrMsg,"PRIMARY KEY must be unique")==0) {
          fprintf(stderr,"Cette erreur est probablement due au fait que le fichier de sortie (%s) \n"
                  "n'a pas ete cree avant l'appel a ce programme avec l'utilitaire 'rdbgen'.  \n", opt.obsout);
        }
        sqlite3_free(ErrMsg);
        EXIT_STATUS = 1;
      }

      /* On ferme la base de donnees ouverte plus haut avec sqlite3_open */
      status = sqlite3_close(sqldb);
      if( status != SQLITE_OK ) {
        fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);
        EXIT_STATUS = 1;
      }

    }
    else { /* Ici, opt.roundrobin == 1 */
      /* Il faut travailler les 'npex*npey' fichiers */
      char npex_str[MAXSTR], npey_str[MAXSTR], format_digits[MAXSTR], rdbout[MAXSTR];
      char sqlschema[MAXSTR*32];
      int nsplit = opt.npex*opt.npey;
      int ilonband, jlatband;

      strcpy(sqlschema,"");
      sprintf(format_digits,"%%.%dd",opt.ndigits);

      for (ilonband=0;ilonband<opt.npex;ilonband++) {
	for (jlatband=0;jlatband<opt.npey;jlatband++) {
	  int id = ilonband*opt.npey+jlatband;

          /* Si on est en mode 'cherrypick', alors on ne considere que si la tuile est egale a celle voulue */
          if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0)
            if (ilonband != opt.cherrypick_x-1 || jlatband != opt.cherrypick_y-1)
              continue;

	  sprintf(npex_str,format_digits,ilonband+1);
	  sprintf(npey_str,format_digits,jlatband+1);
	  sprintf(rdbout,"%s_%s_%s",opt.obsout,npex_str,npey_str);

          status = access(rdbout,F_OK);
          if ( status == 0 ) { /* Le fichier existe deja */
            /* On ouvre la base de donnees SQL de sortie */
            status = sqlite3_open(rdbout,&sqldb);
            if ( status != SQLITE_OK ) {
              fprintf(stderr, "Fonction main: Incapable d'ouvrir la base de donnees: %s\n", sqlite3_errmsg(sqldb));
              exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
            } /* Fin du 'if ( status != SQLITE_OK )' */
          } /* Fin du 'if ( status == 0 )' */
          else {
            if (strlen(sqlschema)==0) {
              /* Si le fichier n'existe pas encore alors on doit le creer et construire le meme schema */
              sqlite3 *sqldbin;

              /* On ouvre le fichier d'input */
              status = sqlite3_open(opt.obsin,&sqldbin);
              if ( status != SQLITE_OK ) {
                fprintf(stderr, "Fonction main: Incapable d'ouvrir la base de donnees: %s\n", sqlite3_errmsg(sqldbin));

                status = sqlite3_close(sqldbin);
                if( status != SQLITE_OK )
                  fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);

                exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
              } /* Fin du 'if ( status != SQLITE_OK )' */

              /* Execution de la requete SQL sur la base de donnees */
              /* L'idee est de reproduire la commande UNIX
                 echo .schema | sqlite3 obsin | sqlite3 obsout
                 Cette requete provient de la documentation http://www.sqlite.org/faq.html#q7
              */
              status = sqlite3_exec(sqldbin, "select * from sqlite_master", sqlite_schema_callback, sqlschema, &ErrMsg);
              if( status != SQLITE_OK ) {
                fprintf(stderr, "Fonction main: Erreur %d pour le fichier '%s' dans la fonction sqlite3_exec: %s\n", status, opt.obsin, ErrMsg);
                sqlite3_free(ErrMsg);
              
                status = sqlite3_close(sqldbin);
                if( status != SQLITE_OK )
                  fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, opt.obsin);

                exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
              } /* Fin du 'if ( status != SQLITE_OK )' */

              printf("Voici le schema du fichier d'input: '%s'\n%s\n", opt.obsin, sqlschema);

              status = sqlite3_close(sqldbin);
              if( status != SQLITE_OK ) {
                fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, opt.obsin);

                exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
              }
            } /* Fin du 'if (strlen(sqlschema)==0)' */

            status = sqlite3_open(rdbout,&sqldb);
            if ( status != SQLITE_OK ) {
              fprintf(stderr, "Fonction main: Incapable d'ouvrir la base de donnees: %s\n", sqlite3_errmsg(sqldb));
              exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
            } /* Fin du 'if ( status != SQLITE_OK )' */

            status = sqlite3_exec(sqldb, sqlschema, (void*) NULL, (void*) NULL, &ErrMsg);
            if( status != SQLITE_OK ){
              fprintf(stderr, "Fonction main: Erreur %d dans la fonction sqlite3_exec: %s\n", status, ErrMsg);
              if (strcmp(ErrMsg,"PRIMARY KEY must be unique")==0) {
                fprintf(stderr,"Cette erreur est probablement due au fait que le fichier de sortie (%s) \n"
                        "n'a pas ete cree avant l'appel a ce programme avec l'utilitaire 'rdbgen'.  \n", opt.obsout);
              }
              sqlite3_free(ErrMsg);
              EXIT_STATUS = 1;
            }
          } /* Fin du 'else' relie au 'if ( status == 0 )' */

          /* On doit fabriquer la requete sql pour faire le splitting */
          sprintf(requete_sql,"drop index if exists idx1;\n"
"PRAGMA journal_mode = OFF;\n"
"PRAGMA  synchronous = OFF;\n"
"attach '%s' as dbin; \n"
"insert into header select * from dbin.header where id_obs %% %d = %d;\n"
"insert into data   select * from dbin.data   where id_obs %% %d = %d;%s\n",
/* "insert into header select * from dbin.header where min (  id_obs/(${maxid}/%d),%d)  = %d;\n" */
/* "insert into data   select * from dbin.data   where min (  id_obs/(${maxid}/%d),%d)  = %d;\n" */
/* "CREATE TABLE rdb4_schema( schema  varchar(9) );\n" */
/* "insert into rdb4_schema values('${TYPE}');\n" */
/* "create table resume(date integer , time integer , run varchar(9)) ;\n" */
/* "insert into resume values(\"$DATE\",\"$HEURE\",\"$RUN\") ;\n" */
                  opt.obsin,nsplit,id,nsplit,id,sqlreq_resume);
          append_id_obs_table_list_requests(requete_sql,table_list);
          strcat(requete_sql,"create index idx1 on data(id_obs,vcoord,varno);");
          strcat(requete_sql,"detach dbin;");

          printf("\nVoici la requete SQL effectuee sur la base de donnees pour creer le fichier '%s':\n",rdbout);
          printf("%s\n", requete_sql);

          /* Execution de la requete SQL sur la base de donnees finale */
          status = sqlite3_exec(sqldb, requete_sql, (void*) NULL, (void*) NULL, &ErrMsg);
          if( status != SQLITE_OK ){
            fprintf(stderr, "Fonction main: Erreur %d dans la fonction sqlite3_exec: %s\n", status, ErrMsg);
            if (strcmp(ErrMsg,"PRIMARY KEY must be unique")==0) {
              fprintf(stderr,"Cette erreur est probablement due au fait que le fichier de sortie (%s) \n"
                      "n'a pas ete cree avant l'appel a ce programme avec l'utilitaire 'rdbgen'.  \n", opt.obsout);
            }
            sqlite3_free(ErrMsg);
            EXIT_STATUS = 1;
          }

          /* On ferme la base de donnees ouverte plus haut avec sqlite3_open */
          status = sqlite3_close(sqldb);
          if( status != SQLITE_OK ) {
            fprintf(stderr,"Fonction main: Erreur %d de la fonction sqlite3_close\n", status);
            EXIT_STATUS = 1;
          }
        } /* Fin du 'for (jlatband=0;jlatband<opt.npey;jlatband++)' */
      } /* Fin du 'for (ilonband=0;ilonband<opt.npex;ilonband++)' */
    } /* Fin du 'else' relie au 'if (opt.roundrobin == 0)' */
  } /* Fin du  if ( filetype == WKF_SQLite ) */
  else if ( filetype == WKF_BURP ) {  /* Alors on traite un fichier BURP */
    int i, i_obs_enrgs, i_enrgs, iout = iun+1, nombre_enregistrements, longueur_max_enregistrement, engrs_resume;
    int *adresses = (int*) NULL, *iouts = (int*) NULL;
    int *t_in_domain = (int*) NULL, *nts = (int*) NULL, *num_obs_per_tile = (int*) NULL;
    int vertical_clipping;
    /* Cette variable sert a identifier si on est en presence du format UA multi-niveaux
     * Si -1, alors on n'a pas encore evaluer le cas, si 0 alors ce sont des UA classiques
     * si 1, ce sont des ua4d.
     */
    int is_ua4d = -1;
    char errmsg[MAXSTR];
    float lat, lon;
    int ilonband, jlatband;
    BURP_RPT *rptin, **rptout;

    /* niveau de tolerance erreur burp */
    /* status = brp_SetOptChar ( "MSGLVL",  "INFORMATIF" ); */
    status = brp_SetOptChar ( "MSGLVL",  "FATAL" );
    if ( status<0 ) {
      fprintf(stderr,"Fonction main: Erreur %d avec la fonction brp_SetOptChar\n", status);

      /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
      if ( opt.roundrobin == 0 ) {
	status = c_gdrls(grid.gridid);
	if (status<0)
	  fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      }

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

    /* Ouverture du fichier burp en entree */
    status = brp_open(iun,opt.obsin,"r");
    if ( status<0 ) {
      fprintf(stderr,"Fonction main: Erreur %d avec la fonction brp_open sur le fichier '%s'\n", status, opt.obsin);

      /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
      if ( opt.roundrobin == 0 ) {
	status = c_gdrls(grid.gridid);
	if (status<0)
	  fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      }

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }
    else if (status == 0) {
      /* Si le fichier BURP est vide alors on doit sortir */
      fprintf(stderr,"Fonction main: Il n'y a aucun enregistrement dans le fichier BURP %s\n",opt.obsin);

      status = brp_close(iun);
      if (status<0)
	fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);

      /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
      if ( opt.roundrobin == 0 ) {
	status = c_gdrls(grid.gridid);
	if (status<0)
	  fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      }

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

    /* Si aucun probleme alors 'status' est le nombre d'enregistrements du fichier BURP */
    nombre_enregistrements = status;
    adresses = (int*) malloc(sizeof(int)*nombre_enregistrements);
    if ( adresses == (int*) NULL) {
      fprintf(stderr,"Fonction main: Incapable d'allouer un vecteur de (int) de dimension %d\n", nombre_enregistrements);

      status = brp_close(iun);
      if (status<0)
	fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);

      /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
      if ( opt.roundrobin == 0 ) {
	status = c_gdrls(grid.gridid);
	if (status<0)
	  fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      }

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);

    }

    nts = (int*) malloc(opt.npex*opt.npey*sizeof(int));
    if ( nts == (int*) NULL) {
      fprintf(stderr,"Fonction main: Incapable d'allouer le vecteur nts de (int) de dimension %dx%d\n", opt.npex, opt.npey);

      status = brp_close(iun);
      if (status<0)
	fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);
	
      /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
      if ( opt.roundrobin == 0 ) {
	status = c_gdrls(grid.gridid);
	if (status<0)
	  fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      }

      free(adresses);
	
      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    } /* Fin du if ( nts == (int*) NULL) */

    num_obs_per_tile = (int*) malloc(opt.npex*opt.npey*sizeof(int));
    if ( num_obs_per_tile == (int*) NULL) {
      fprintf(stderr,"Fonction main: Incapable d'allouer le vecteur num_obs_per_tile de (int) de dimension %dx%d\n", opt.npex, opt.npey);

      status = brp_close(iun);
      if (status<0)
	fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);

      /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
      if ( opt.roundrobin == 0 ) {
	status = c_gdrls(grid.gridid);
	if (status<0)
	  fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      }

      free(adresses);
      free(nts);

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    } /* Fin du if ( nts == (int*) NULL) */

    iouts = (int*) malloc(opt.npex*opt.npey*sizeof(int));
    if ( iouts == (int*) NULL) {
      fprintf(stderr,"Fonction main: Incapable d'allouer le vecteur iouts de (int) de dimension %dx%d\n", opt.npex, opt.npey);

      status = brp_close(iun);
      if (status<0)
	fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);
	
      /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
      if ( opt.roundrobin == 0 ) {
	status = c_gdrls(grid.gridid);
	if (status<0)
	  fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      }

      free(adresses);
      free(nts);
      free(num_obs_per_tile);
	
      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    } /* Fin du if ( iouts == (int*) NULL) */

      /* Ouverture du fichier burp en sortie */
    if ( opt.npex == 1 && opt.npey == 1 ) {

      status = access(opt.obsout,F_OK);
      if ( status==0 ) {
	fprintf(stderr,"Fonction main: Le fichier '%s' existe deja mais il pourrait etre efface "
		"alors il vaut mieux que ce fichier n'existe pas a l'appel du programme\n", opt.obsout);

	status = brp_close(iun);
	if (status<0)
	  fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);

	/* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
	if ( opt.roundrobin == 0 ) {
	  status = c_gdrls(grid.gridid);
	  if (status<0)
	    fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
	}

	free(adresses);
	free(nts);
	free(num_obs_per_tile);
	free(iouts);

	exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
      }


      iouts[0] = iout;
      status = brp_open(iout,opt.obsout,"a");
      if ( status<0 ) {
	fprintf(stderr,"Fonction main: Erreur %d avec la fonction brp_open sur le fichier '%s'\n", status, opt.obsout);
	
	status = brp_close(iun);
	if (status<0)
	  fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);

	status = brp_close(iout);
	if (status<0)
	  fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsout);
	
	/* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
	if ( opt.roundrobin == 0 ) {
	  status = c_gdrls(grid.gridid);
	  if (status<0)
	    fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
	}

	free(adresses);
	free(num_obs_per_tile);
	free(nts);
	free(iouts);

	exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
      }
    }
    else { /* Il faut ouvrir npex*npey fichiers */
      char npex_str[MAXSTR], npey_str[MAXSTR], format_digits[MAXSTR], burpout[MAXSTR];

      sprintf(format_digits,"%%.%dd",opt.ndigits);
      
      for (ilonband=0;ilonband<opt.npex;ilonband++) {
	for (jlatband=0;jlatband<opt.npey;jlatband++) {
	  int id=ilonband*opt.npey+jlatband;

          /* Si on est en mode 'cherrypick', alors on ne considere que si la tuile est egale a celle voulue */
          if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0)
            if (ilonband != opt.cherrypick_x-1 || jlatband != opt.cherrypick_y-1)
              continue;

	  iouts[id] = iout++;

	  sprintf(npex_str,format_digits,ilonband+1);
	  sprintf(npey_str,format_digits,jlatband+1);
	  sprintf(burpout,"%s_%s_%s", opt.obsout,npex_str,npey_str);

	  if (VERBOSE>5)
	    printf("Fonction main: appel de 'brp_open' sur le fichier '%s' pour id=%d\n", burpout, id);

	  status = access(burpout,F_OK);
	  if ( status==0 || iouts[id]>999 ) {
	    int ilonbandtmp, jlatbandtmp;

            if ( status==0 )
              fprintf(stderr,"Fonction main: Le fichier '%s' existe deja mais il pourrait etre efface "
                      "alors il vaut mieux que ce fichier n'existe pas a l'appel du programme\n", burpout);
            else if ( iouts[id]>999 )
              fprintf(stderr,"Fonction main: Comme iouts[%d]=%d, la commande 'brp_open(...)' "
                      "n'acceptera pas cette valeur puisqu'elle est plus grande que 999\n", id, iouts[id]);

	    status = brp_close(iun);
	    if (status<0)
	      fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close pour le fichier %s\n", status, opt.obsin);

	    for (ilonbandtmp=0;ilonbandtmp<ilonband;ilonbandtmp++) {
	      for (jlatbandtmp=0;jlatbandtmp<jlatband;jlatbandtmp++) {
		status = brp_close(iouts[ilonbandtmp*opt.npey+jlatbandtmp]);
		if ( status<0 )
		  fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close "
			  "pour le fichier %s_%d_%d\n", status, opt.obsout,ilonbandtmp+1,jlatbandtmp+1);
	      }
	    }

	    /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
	    if ( opt.roundrobin == 0 ) {
	      status = c_gdrls(grid.gridid);
	      if (status<0)
		fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
	    }

	    free(adresses);
	    free(num_obs_per_tile);
	    free(iouts);
	    free(nts);

	    exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	  } /* Fin du if ( status==0 ) pour l'acces au fichier */

	  status = brp_open(iouts[id],burpout,"a");
	  if ( status<0 ) {
	    int ilonbandtmp, jlatbandtmp;

	    fprintf(stderr,"Fonction main: Erreur %d avec la fonction brp_open sur le fichier %s\n", status, burpout);

	    status = brp_close(iun);
	    if (status<0)
	      fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close pour le fichier %s\n", status, opt.obsin);

	    for (ilonbandtmp=0;ilonbandtmp<ilonband;ilonbandtmp++) {
	      for (jlatbandtmp=0;jlatbandtmp<jlatband;jlatbandtmp++) {
		status = brp_close(iouts[ilonbandtmp*opt.npey+jlatbandtmp]);
		if ( status<0 )
		  fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close "
			  "pour le fichier %s_%d_%d\n", status, opt.obsout,ilonbandtmp+1,jlatbandtmp+1);
	      }
	    }

	    /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
	    if ( opt.roundrobin == 0 ) {
	      status = c_gdrls(grid.gridid);
	      if (status<0)
		fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
	    }

	    free(adresses);
	    free(iouts);
	    free(nts);
	    free(num_obs_per_tile);

	    exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	  } /* Fin du if ( status<0 ) */
	} /* Fin du for (jlatband=0;jlatband<opt.npey;jlatband++) */
      } /* Fin du for (ilonband=0;ilonband<opt.npex;ilonband++) */
    } /* Fin du else pour le if ( opt.npex == 1 && opt.npey == 1 ) */

    /* Toujours initialiser les pointeurs */
    rptout = (BURP_RPT**) malloc(opt.npex*opt.npey*sizeof(BURP_RPT*));
    if ( rptout == (BURP_RPT**) NULL ) {
      fprintf(stderr,"Fonction main: Incapable d'allouer le vecteur rptout de (BURP_RPT*) de dimension %dx%d\n", opt.npex, opt.npey);

      status = brp_close(iun);
      if (status<0)
	fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close pour le fichier %s\n", status, opt.obsin);

      /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
      if ( opt.roundrobin == 0 ) {
	status = c_gdrls(grid.gridid);
	if (status<0)
	  fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      }

      for (ilonband=0;ilonband<opt.npex;ilonband++)
	for (jlatband=0;jlatband<opt.npey;jlatband++) {
	  status = brp_close(iouts[ilonband*opt.npey+jlatband]);
	  if ( status<0 )
	    fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close "
		    "pour le fichier %s_%d_%d\n", status, opt.obsout,ilonband+1,jlatband+1);
	}

      free(adresses);
      free(iouts);
      free(nts);
      free(num_obs_per_tile);

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

    rptin  = brp_newrpt();
    if ( rptin == (BURP_RPT*) NULL ) {
      fprintf(stderr,"Fonction main: Incapable d'allouer rptin de (BURP_RPT*)\n");
      status = brp_close(iun);
      if (status<0)
	fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close pour le fichier %s\n", status, opt.obsin);

      /* Si on n'est pas en mode round-robin, alors on a eu besoin du fichier 'opt.fstin'. */
      if ( opt.roundrobin == 0 ) {
	status = c_gdrls(grid.gridid);
	if (status<0)
	  fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      }

      for (ilonband=0;ilonband<opt.npex;ilonband++)
	for (jlatband=0;jlatband<opt.npey;jlatband++) {
	  status = brp_close(iouts[ilonband*opt.npey+jlatband]);
	  if ( status<0 )
	    fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close "
		    "pour le fichier %s_%d_%d\n", status, opt.obsout,ilonband+1,jlatband+1);
	}

      free(adresses);
      free(iouts);
      free(nts);
      free(num_obs_per_tile);
      free(rptout);

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

    i_enrgs=0;
    RPT_SetHANDLE(rptin, 0);
    while ( brp_findrpt(iun, rptin) >= 0 )
      adresses[i_enrgs++] = RPT_HANDLE(rptin);

    longueur_max_enregistrement = c_mrfbfl(iun);
    if (VERBOSE>0)
	printf("Fonction main: les rapports auront la longueur %d\n", longueur_max_enregistrement);
    /* Il faut absolument ajouter a la longueur maximal d'un
     * enregistrement sinon l'appel a brp_freerpt(rptout) donnera un
     * signal '*** glibc detected *** free(): invalid pointer: 0x08314590 ***'
     */
    for (i=0;i<opt.npex*opt.npey;i++) {
      rptout[i] = brp_newrpt();
      brp_allocrpt(rptout[i], longueur_max_enregistrement);
    }

    /* On regarde si on doit fait le clipping vertical */
    if ( strlen(opt.channels)!=0 || opt.niveau_min != IP1_VIDE || opt.niveau_max != IP1_VIDE)
      vertical_clipping = 1;
    else
      vertical_clipping = 0;

    for (i=0;i<opt.npex*opt.npey;i++) {
      nts[i]=0;
      num_obs_per_tile[i]=0;
    }

    /* Cette variable contient le nombre d'enregistrements d'observation qui ont ete traites jusqu'a maintenant dans la boucle */
    i_obs_enrgs = 0;
    /* Ensuite, on passe chaque enregistrement un a un */
    for (i_enrgs=0;i_enrgs<nombre_enregistrements;i_enrgs++) {
      engrs_resume = 0;

      status = brp_getrpt(iun,adresses[i_enrgs],rptin);
      if (status<0) {
	fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_getrpt pour "
		"le fichier d'entree %s a l'adresse %d (%d rapport)\n", status, opt.obsin, adresses[i_enrgs], i_enrgs);
	continue;
      }

      if (VERBOSE>0)
	printf("stnids = %s  enregistrement: %d\n", RPT_STNID(rptin), i_enrgs);

      /* entete de rapport */
      if (VERBOSE>1) {
	printf ( "\n\n" ) ;
	printf ( "hhmm   =%8d " , RPT_TEMPS(rptin)) ;
	printf ( "flgs   =%6d  ", RPT_FLGS(rptin )) ;
	printf ( "codtyp =%6d  ", RPT_IDTYP(rptin)) ;
	printf ( "stnids =%9s\n", RPT_STNID(rptin)) ;
	printf ( "blat   =%8d " , RPT_LATI(rptin)) ;
	printf ( "blon   =%6d  ", RPT_LONG(rptin )) ;
	printf ( "dx     =%6d  ", RPT_DX(rptin)) ;
	printf ( "dy     =%6d  ", RPT_DY(rptin)) ;
	printf ( "stnhgt =%6d\n", RPT_ELEV(rptin)) ;
	printf ( "yymmdd =%8d " , RPT_DATE(rptin)) ;
	printf ( "oars   =%6d  ", RPT_OARS(rptin)) ;
	printf ( "runn   =%6d  ", RPT_RUNN(rptin)) ;
	printf ( "nblk   =%6d  ", RPT_NBLK(rptin)) ;
	printf ( "dlay   =%6d\n", RPT_DRND(rptin)) ;
	printf ( "\n" ) ;
      }

      if (VERBOSE>4) {
	BURP_BLK* blk = (BURP_BLK*) NULL;

	printf("Impression des blocs au debut\n");

	blk = brp_newblk();
	BLK_SetBKNO(blk, 0);

	while ( brp_findblk( blk, rptin ) >= 0 ) {
	  int thisi, thisj, thisk;
	  BURP_BLK* blkout = (BURP_BLK*) NULL;
	  blkout = brp_newblk();
	  status = brp_readblk(BLK_BKNO(blk), blkout, rptin, 0);

	  printf ( "blkno  =%6d  ", BLK_BKNO(blkout)    ) ;
	  printf ( "nele   =%6d  ", BLK_NELE(blkout)    ) ;
	  printf ( "nval   =%6d  ", BLK_NVAL(blkout)    ) ;
	  printf ( "nt     =%6d  ", BLK_NT(blkout)      ) ;
	  printf ( "bit0   =%6d\n", BLK_BIT0(blkout)    ) ;
	  printf ( "bdesc  =%6d  ", BLK_BDESC(blkout)   ) ;
	  printf ( "btyp   =%6d  ", BLK_BTYP(blkout)    ) ;
	  printf ( "nbit   =%6d  ", BLK_NBIT(blkout)    ) ;
	  printf ( "datyp  =%6d  ", BLK_DATYP(blkout)   ) ;
	  printf ( "bfam   =%6d  ", BLK_BFAM(blkout)    ) ;

	  for ( thisk = 0 ; thisk < BLK_NT(blkout) ; thisk++ ) {
	    printf (  "\nlstele =" ) ;
	    for (thisi=0;thisi<BLK_NELE(blkout);thisi++)
	      printf (  "    %6.6d", BLK_DLSTELE(blkout,thisi) ) ;
	    /* sortie des valeurs des elements */
	    for (  thisj = 0 ; thisj < BLK_NVAL(blkout) ; thisj++ ) {
	      printf (  "\ntblval =" ) ;
	      for (  thisi = 0 ; thisi < BLK_NELE(blkout) ; thisi++ )
		printf (  "%10d", BLK_TBLVAL(blkout,thisi,thisj,thisk) ) ;
	    }
	  }
	  printf("\n");
	  brp_freeblk(blkout);
	} /* Fin du while ( brp_findblk( blk, rptin ) >= 0 ) */
	brp_freeblk(blk);
      } /* Fin du if (VERBOSE>4) */

      /* On verifie si on est en presence d'un ua4d.
       * On ne fait ceci que si on n'est pas en presence d'un enregistrement resume.
       */
      if ((is_ua4d==-1 || opt.check_ua4d) && strncmp(">>",RPT_STNID(rptin),2)!=0) {
	is_ua4d = check_ua4d(rptin);
	if (is_ua4d<0) {
	  fprintf(stderr,"Fonction main: la fonction 'check_ua4d' retourne %d\n"
		  "le fichier d'entree %s a l'adresse %d (%d rapport)\n",
		  is_ua4d, opt.obsin, adresses[i_enrgs], i_enrgs);
	  continue;
	}
      }

      if (strncmp(">>",RPT_STNID(rptin),2)==0) { /* C'est un enregistrement resume */
	engrs_resume = 1;
	if (VERBOSE>0)
	  printf("Fonction main: C'est un enregistrement resume '%s' (i_enrgs=%d)\n", 
		 RPT_STNID(rptin),i_enrgs);
      }

      if (engrs_resume==0 && vertical_clipping==1) {
	BURP_RPT *rptin_clip_vert;

	rptin_clip_vert = brp_newrpt();
	brp_allocrpt(rptin_clip_vert, longueur_max_enregistrement);
	
	status = clipping_vertical(rptin,&opt,&grid_gz,rptin_clip_vert);
	if (status<0) {
	  fprintf(stderr,"Fonction main: Erreur %d dans la fonction clipping_vertical\n", status);
	  EXIT_STATUS=-1;
	  continue;
	}
	/* Apres l'etape du clipping vertical, le 'rptin_clip_vert' devient le 'rptin' */
	brp_freerpt(rptin);
	rptin = rptin_clip_vert;
      } /* Fin du if (engrs_resume==0 && vertical_clipping==1) */

      for (i=0;i<opt.npex*opt.npey;i++) {
	/* On se doit de remettre a 0 tous les elements de 'nts'
	 * puisqu'ils vont nous indiquer combien d'observations
	 * contient chaque domaine.
	 */
	nts[i]=0;
	/* On 'claire' tous les rapports avant de travailler dedans */
	/* Il n'est pas necessaire de faire ceci et on passer 75% du
	   temps sur cette ligne alors on l'enleve. */
	/* brp_clrrpt(rptout[i]); */
      }

      /* On commence par copier le header dans tous les rapports */
      for(i=0;i<opt.npex*opt.npey;i++) {
        /* Si on est en mode 'cherrypick', alors on ne considere que si la tuile est egale a celle voulue */
        if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0) {
          int cherrypick_id = (opt.cherrypick_x-1)*opt.npey+opt.cherrypick_y-1;
          if (i != cherrypick_id)
            continue;
        }

	brp_copyrpthdr(rptout[i],rptin);

	if (VERBOSE>4)
	  printf("Fonction main: appel de 'brp_putrpthdr' pour i=%d, iouts=%d et i_enrgs=%d\n",i,iouts[i],i_enrgs);

	status = brp_putrpthdr(iouts[i],rptout[i]);
	if (status<0) {
	  fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_putrpthdr pour "
		  "le fichier de sortie %s_%d_%d a l'adresse %d (rapport %d) (iun=%d)\n",
		  status, opt.obsout, ilonband+1, jlatband+1, adresses[i_enrgs], i_enrgs,
		  iouts[ilonband*opt.npey+jlatband]);
	  EXIT_STATUS = 1;
	  break;
	}
      }

      if (engrs_resume==1) {
	/* Si c'est un enregistrement resume, */
	/*  on veut alors copier cet enregistrement dans le fichier output */
	if (VERBOSE>0)
	  printf("Fonction main: On copie cet enregistrement resume '%s' (i_enrgs=%d)\n", 
		 RPT_STNID(rptin),i_enrgs);
	status = fill_rptout_blk(rptin,rptout,(int*) NULL,(int*) NULL,opt.npex*opt.npey,opt.cherrypick_x,opt.cherrypick_y,opt.npey);
	if (status<0) {
	  fprintf(stderr,"Fonction main: Erreur %d dans la fonction fill_rptout_blk pour "
		  "l'adresse %d (%d rapport)\n", status, adresses[i_enrgs], i_enrgs);
	  EXIT_STATUS = 1;
	  break;
	}
      } /* Fin du if (strncmp(">>",RPT_STNID(rptin),2)==0) */
      else if (opt.roundrobin==1) {
	/* On fait le splitting en mode round-robin comme l'utilitaire 'reflex' le fait. */
	int bin_id = i_obs_enrgs%(opt.npex*opt.npey);

	if (strncmp("^",RPT_STNID(rptin),1)==0)
	  /* Si c'est un enregistrement resume alors on a le nombre d'observations dans 'elev' */
	  nts[bin_id] = RPT_ELEV(rptin);
	else
	  /* Sinon, chaque enregistrement contient une seule observation. */
	  nts[bin_id] = 1;

	if (VERBOSE>0)
	  printf("Fonction main: L'enregistrement %d sera place dans le fichier %d et contient %d observations.\n",
		 i_enrgs,bin_id,nts[bin_id]);

	num_obs_per_tile[bin_id] += nts[bin_id];
	status = fill_rptout_blk(rptin,rptout,nts,(int*) NULL,opt.npex*opt.npey,opt.cherrypick_x,opt.cherrypick_y,opt.npey);
	if (status<0) {
	  fprintf(stderr,"Fonction main: Erreur %d dans la fonction fill_rptout_blk pour "
		  "l'adresse %d (%d rapport)\n", status, adresses[i_enrgs], i_enrgs);
	  EXIT_STATUS = 1;
	  break;
	}
      } /* Fin du 'else if (opt.roundrobin==1)' */
      /* Cette section traite les enregistrements regroupes qu'ils soient des 'ua4d' ou bien de
       * donnees satellitaires.
       */
      else if (strncmp("^",RPT_STNID(rptin),1)==0) {
	int obs_in_domain=0;

	/* Si c'est un enregistrement regroupe de donnes satellitaires */
	/* 5120 est le btyp du bloc info pour ces donnees */
	/* 5002 est l'element qui donne la latitude de l'observation */
	/* 6002 est l'element qui donne la longitude de l'observation */
	status = extract_data_in_domains_along_nt(&opt,&grid,rptin,5002,6002,nts,&t_in_domain);
	if (status<0) {
	  fprintf(stderr,"Fonction main: Erreur %d dans la fonction extract_data_in_domains_along_nt pour "
		  "l'adresse %d (%d rapport)\n", status, adresses[i_enrgs], i_enrgs);
	  EXIT_STATUS = 1;
	  break;
	}
	
	/* Si aucune observation n'est dans le domaine alors on passe au
	 * prochain enregistrement 
	 */
	for (i=0;i<opt.npex*opt.npey;i++)
	  if (nts[i]!=0) {
	    num_obs_per_tile[i] += nts[i];
	    obs_in_domain=1;
	    if (VERBOSE>1)
	      printf("Fonction main: %d observations sont dans la bande %d pour cet enregistrement\n", nts[i], i);
	  }
	if (obs_in_domain==0) {
	  free(t_in_domain);
	  t_in_domain = (int*) NULL;
	  if (VERBOSE>1)
	    printf("Fonction main: Aucune observation n'a ete selectionnee pour cet enregistrement\n");
	  continue;
	}

	/* Dans le cas des observations regroupees, l'elevation de la
	 * station dans le header contient le nombre d'observations
	 * dans l'enregistrement. Manifestement, ce nombre change
	 * lorsqu'on clippe ou bien on splitte alors il faut mettre a
	 * jour cette information.
	 */
	for (i=0;i<opt.npex*opt.npey;i++) {
	  if (nts[i]==0) continue;

          /* Si on est en mode 'cherrypick', alors on ne considere que si la tuile est egale a celle voulue */
          if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0) {
            int cherrypick_id = (opt.cherrypick_x-1)*opt.npey+opt.cherrypick_y-1;
            if (i != cherrypick_id)
              continue;
          }

	  if (RPT_ELEV(rptout[i]) != nts[i]) {
	    if (VERBOSE>5) {
	      printf("Fonction main: On met a jour le header (RPT_ELEV(rptout[%d])) de %d a %d\n", i, RPT_ELEV(rptout[i]),nts[i]);
	    }
	    RPT_SetELEV(rptout[i],nts[i]);
	    status = brp_updrpthdr(iouts[i],rptout[i]);
	    if (status<0) {
	      fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_updrpthdr pour "
		      "le fichier de sortie %s_%d_%d a l'adresse %d (%d rapport)\n",
		      status, opt.obsout, i/opt.npey+1, i%opt.npey+1, adresses[i_enrgs], i_enrgs);
	      EXIT_STATUS = 1;
	      break;
	    }
	  } /* Fin du if (RPT_ELEV(rptout[i]) != nts[i]) */
	  else if (VERBOSE>5) {
	    printf("Fonction main: On n'a pas pas besoin de mettre a jour le header (RPT_ELEV(rptout[%d])=%d)\n", i, RPT_ELEV(rptout[i]));
	  }
	} /* Fin du for (i=0;i<opt.npex*opt.npey;i++) */
	
	if (EXIT_STATUS)
	  break;

	status = fill_rptout_blk(rptin,rptout,nts,t_in_domain,opt.npex*opt.npey,opt.cherrypick_x,opt.cherrypick_y,opt.npey);
	if (status<0) {
	  fprintf(stderr,"Fonction main: Erreur %d dans la fonction fill_rptout_blk pour "
		  "l'adresse %d (%d rapport)\n", status, adresses[i_enrgs], i_enrgs);
	  EXIT_STATUS = 1;
	  break;
	}
      } /* Fin du else if (strncmp("^",RPT_STNID(rptin),1)==0) */
      else if (is_ua4d==1) { /* Si c'est un 'ua4d' */
	int obs_in_domain=0, i_blk, i_btyp_data, i_btyp, btyp_data, btypnum;
	int *bknos_data, *btyps_data, *nts_for_this_blk, *nvals_in_domain, *val_in_domain;

	bknos_data = (int*) NULL;
	btyps_data = (int*) NULL;
	status = find_blk_data_in_rpt(rptin,5001,6001,&bknos_data,&btyps_data,&btypnum);
	if (status<0) {
	  fprintf(stderr,"Fonction main: Erreur dans la fonction find_blk_data_in_rpt\n");
	  if (bknos_data != (int*) NULL)  free(bknos_data);
	  continue;
	}

	val_in_domain   = (int*) NULL;
	nvals_in_domain = (int*) NULL;
	nvals_in_domain = (int*) malloc(opt.npex*opt.npey*sizeof(int));
	if (nvals_in_domain == (int*) NULL) {
	  fprintf(stderr,"Fonction main: Incapable d'allouer un vecteur de int de dimension %d pour le cas 'ua4d'\n", opt.npex*opt.npey);
	  if (bknos_data != (int*) NULL)  free(bknos_data);
	  continue;
	}
	for (i=0;i<opt.npex*opt.npey;i++)
	  nvals_in_domain[i]=0;

	/* On passe au travers tous les blocs de donnees et on va extraire les donnees */
	for (i_btyp=0;i_btyp<btypnum;i_btyp++) {
	  BURP_BLK* blkdata = (BURP_BLK*) NULL;
	  BURP_BLK* blksearch = (BURP_BLK*) NULL;

	  blkdata = brp_newblk();
	  status = brp_readblk(bknos_data[i_btyp],blkdata,rptin,0);
	  if (status<0) {
	    fprintf(stderr,"Fonction main: Erreur dans la fonction brp_readblk pour bknos_data[%d]=%d\n", i_btyp, bknos_data[i_btyp]);
	    brp_freeblk(blkdata);
	    continue;
	  }
	  btyp_data = BLK_BTYP(blkdata);
	  if (btyp_data != btyps_data[i_btyp])
	    fprintf(stderr,"Fonction main: Probleme potentiel puisque btyp_data=%d est different de btyps_data[%d]=%d",
		    btyp_data, i_btyp, btyps_data[i_btyp]);

	  val_in_domain = (int*) NULL;
	  val_in_domain = (int*) malloc(opt.npex*opt.npey*BLK_NVAL(blkdata)*sizeof(int));
	  if (val_in_domain == (int*) NULL) {
	    fprintf(stderr,"Fonction main: Incapable d'allouer un vecteur de int de dimension %d pour le cas 'ua4d'\n", opt.npex*opt.npey*BLK_NVAL(blkdata));
	    if (bknos_data != (int*) NULL) free(bknos_data);
	    brp_freeblk(blkdata);
	    continue;
	  }
	  for (i=0;i<opt.npex*opt.npey*BLK_NVAL(blkdata);i++)
	    val_in_domain[i]=0;

	  status = extract_data_in_domains_along_nval(&opt,&grid,rptin,5001,6001,blkdata,nvals_in_domain,val_in_domain);
	  if (status<0) {
	    fprintf(stderr,"Fonction main: Erreur %d dans la fonction extract_data_in_domains_along_nval pour "
		    "l'adresse %d (%d rapport) et le bloc %d\n", status, adresses[i_enrgs], i_enrgs, i_btyp);
	    EXIT_STATUS = 1;
	    brp_freeblk(blkdata);
	    free(val_in_domain);
	    break;
	  }

	  /* Si aucune observation n'est dans le domaine alors on passe au
	   * prochain enregistrement
	   */
	  for (i=0;i<opt.npex*opt.npey;i++)
	    if (nvals_in_domain[i]!=0) {
	      num_obs_per_tile[i] += nvals_in_domain[i];
	      obs_in_domain=1;
	      if (VERBOSE>1)
		printf("Fonction main: %d observations sont dans la bande %d pour cet enregistrement\n", nvals_in_domain[i], i);
	    }
	  if (obs_in_domain==0) {
	    free(val_in_domain);
	    val_in_domain = (int*) NULL;
	    if (VERBOSE>1)
	      printf("Fonction main: Aucune observation n'a ete selectionnee pour cet "
		     "enregistrement (rapport %d, bkno=%d, btyp=%d)\n", i_enrgs, bknos_data[i_btyp], btyp_data);
	    continue;
	  }

	  /* A partir d'ici, on a au moins une observation dans le domaine */
	  blksearch = brp_newblk();
	  BLK_SetBKNO(blksearch, 0);
	  while ( brp_findblk( blksearch, rptin ) >= 0 ) {
	    BURP_BLK* blkout = (BURP_BLK*) NULL;
	    blkout = brp_newblk();
	    status = brp_readblk(BLK_BKNO(blksearch), blkout, rptin, 0);
	    if (status<0) {
	      fprintf(stderr,"Fonction main: Erreur dans la fonction brp_readblk pour bkno=%d\n", BLK_BKNO(blksearch));
	      brp_freeblk(blkout);
	      continue;
	    }

	    if ( btyp_data == BLK_BTYP(blkout) || btypAssociated(btyp_data,BLK_BTYP(blkout)) == 1 )
	      if (BLK_NVAL(blkdata) == BLK_NVAL(blkout) ) {
		for (i=0;i<opt.npex*opt.npey;i++) {
                  /* Si on est en mode 'cherrypick', alors on ne
                     considere que si la tuile est egale a celle
                     voulue */
                  if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0) {
                    int cherrypick_id = (opt.cherrypick_x-1)*opt.npey+opt.cherrypick_y-1;
                    if (i != cherrypick_id)
                      continue;
                  }
		  if (nvals_in_domain[i]!=0) {
		    if (VERBOSE>1)
		      printf("Fonction main: Appel de putblk_nval btyp=%d i=%d nvals=%d BLK_NVAL(blkout)*i=%d\n",
			     BLK_BTYP(blkout), i, nvals_in_domain[i],BLK_NVAL(blkout)*i);

		    status = putblk_nval(rptout[i],blkout,&val_in_domain[BLK_NVAL(blkout)*i], nvals_in_domain[i]);
		    if (status<0) {
		      fprintf(stderr,"Fonction main: Erreur %d dans la fonction putblk_nval pour "
			      "l'adresse %d (%d rapport) et le bloc %d (btyp=%d) pour la tuile %d (bloc data)\n",
			      status, adresses[i_enrgs], i_enrgs, i_btyp, BLK_BTYP(blkout), i);
		      EXIT_STATUS = 1;
		      break;
		    }
		  }
                } /* for (i=0;i<opt.npex*opt.npey;i++) */
	      } /* Fin du 'if (btyp_data == BLK_BTYP(blkout) || btypAssociated(btyp_data,BLK_BTYP(blkout)) == 1)' */
	      else {
		fprintf(stderr,"Fonction main: Erreur dans le traitement de l'enregistrement a "
			"l'adresse %d (%d rapport) et le bloc %d (btyp=%d, nval=%d).  Le nval est different "
			"du bloc d'observations associe (bkno=%d btyp=%d nval=%d)\n",
			adresses[i_enrgs], i_enrgs, i_btyp, BLK_BTYP(blkout), BLK_NVAL(blkout),
			BLK_BKNO(blkdata), BLK_BTYP(blkdata), BLK_NVAL(blkdata));
	      }

	    brp_freeblk(blkout);
	  } /* Fin du while ( brp_findblk( blksearch, rptin ) >= 0 ) */

	  free(val_in_domain);
	  val_in_domain = (int*) NULL;

	  for (i=0;i<opt.npex*opt.npey;i++) {
	    nts[i]+=nvals_in_domain[i];
	    nvals_in_domain[i]=0;
	  }

	  brp_freeblk(blkdata);
	  brp_freeblk(blksearch);

	  if (EXIT_STATUS)
	    break;
	} /* Fin du 'for (i_btyp=0;i_btyp<btypnum;i_btyp++)' */

	/* A partir d'ici, tous les blocs de donnees et ceux qui leur sont associes ont ete traites.
	 * Il faut maintenant trouver ceux qui ne l'ont pas ete pour les ecrire dans le fichier
	 */
	{
	  BURP_BLK* blkdata = (BURP_BLK*) NULL;
	  BURP_BLK* blktmp  = (BURP_BLK*) NULL;

	  blkdata = brp_newblk();

	  if (VERBOSE>1)
	    printf("Fonction main: on cherche les blocs qui n'ont pas ete traites\n");

	  /* on trouve les blocs qui ne sont pas associes a aucun bloc */
	  blktmp = brp_newblk();
	  while ( brp_findblk( blktmp, rptin ) >= 0 ) {
	    int associated = 0;
	    BURP_BLK* blkout = (BURP_BLK*) NULL;

	    blkout = brp_newblk();
	    status = brp_readblk(BLK_BKNO(blktmp), blkout, rptin, 0);
	    if (status<0) {
	      fprintf(stderr,"Fonction main: Erreur dans la fonction brp_readblk pour bknos_data[%d]=%d\n", i_btyp, BLK_BKNO(blktmp));
	      brp_freeblk(blkout);
	      continue;
	    }

	    for (i_btyp=0;i_btyp<btypnum;i_btyp++) {
	      status = brp_readblk(bknos_data[i_btyp],blkdata,rptin,0);
	      if (status<0) {
		fprintf(stderr,"Fonction main: Erreur dans la fonction brp_readblk pour bknos_data[%d]=%d\n", i_btyp, bknos_data[i_btyp]);
		brp_freeblk(blkdata);
		continue;
	      }

	      if (BLK_NVAL(blkdata) == BLK_NVAL(blkout) && (btyps_data[i_btyp] == BLK_BTYP(blkout) || btypAssociated(btyps_data[i_btyp],BLK_BTYP(blkout)) == 1)) {
		if (VERBOSE>1)
		  printf("Fonction main: le btyp %d est associe a %d\n", BLK_BTYP(blkout), btyps_data[i_btyp]);
		associated=1;
		break;
	      }
	    }

	    if (!associated) {
	      if (VERBOSE>1)
		printf("Fonction main: le btyp %d n'a associe a aucun autre bloc d'observations\n", BLK_BTYP(blkout));

	      /* Lorsqu'on en trouve qui n'a pas ete associe alors on le copie dans tous les blocs */
	      for (i=0;i<opt.npex*opt.npey;i++) {
                if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0) {
                  int cherrypick_id = (opt.cherrypick_x-1)*opt.npey+opt.cherrypick_y-1;
                  if (i != cherrypick_id)
                    continue;
                }

		status = putblk_nval(rptout[i],blkout,(int*) NULL, 0);
                if (status<0) {
                  fprintf(stderr,"Fonction main: Erreur %d dans la fonction putblk_nval pour "
                          "l'adresse %d (%d rapport) et le bloc %d (btyp=%d) pour la tuile %d (bloc data)\n",
                          status, adresses[i_enrgs], i_enrgs, i_btyp, BLK_BTYP(blktmp), i);
                  EXIT_STATUS = 1;
                  break;
                }
              } /* Fin du for (i=0;i<opt.npex*opt.npey;i++) */
	    } /* Fin du 'if(!associated)' */

	    brp_freeblk(blkout);
	  } /* Fin du 'while ( brp_findblk( blktmp, rptin ) >= 0 )' */

	  brp_freeblk(blkdata);
	} /* Fin du bloc '{ on trouve les blocs qui ne sont pas associes a aucun bloc' */

	free(bknos_data);
	free(btyps_data);
	free(nvals_in_domain);
      } /* Fin du if (is_ua4d==1) */
      else { /* Si ce n'est pas un 'ua4d', ni un enregistrement regroupe ni un enregistrement resume */
	int id;
	/* Dans le cas d'un rapport non-regroupe alors on peut se fier aux valeurs de
	 * latitude et de longitude dans l'entete du rapport.  
	 */
	lat = RPT_LATI(rptin)/100. - 90.;
	lon = RPT_LONG(rptin)/100.;

	if ( opt.npex == 1 && opt.npey == 1 ) {
	  ilonband=1;
	  jlatband=1;
	  status = checkgrid(grid.gridid, grid.ni, grid.nj, lat, lon, opt.rect, errmsg);
	  if (status<0) {
	    fprintf(stderr,"Fonction main: Erreur dans la fonction checkgrid pour le lat=%f "
		    "et lon=%f avec le message '%s'\n", lat, lon, errmsg);
	    EXIT_STATUS = 1;
	    break;
	  }

	  /* Ceci signifie que si opt.inout == 1 alors status == 0 et donc
	   * le point est hors de la grille ce qui n'est pas voulu
	   *
	   * ou bien qui si opt.inout == 0 alors status == 1 et donc le
	   * point est a l'interieur de la grille ce qui n'est pas voulu.
	   */	
	  if (opt.inout != status)
	    continue;
	  id=0;
	  ilonband=1;
	  jlatband=1;
	}
	else {
	  ilonband=-1;
	  jlatband=-1;
	  status = find_subdomain(grid.gridid, grid.ni, grid.nj, lat, lon, opt.rect,
				  opt.npex, opt.npey, &ilonband, &jlatband, errmsg);
	  if (status<0) {
	    fprintf(stderr,"Fonction main: Erreur dans la fonction find_subdomain "
		    "pour le lat=%f et lon=%f avec le message '%s'\n", lat, lon, errmsg);
	    EXIT_STATUS = 1;
	    break;
	  }
	  if (ilonband<1 || ilonband>opt.npex || jlatband<1 || jlatband>opt.npey) {
	    if (VERBOSE>1)
	      printf("Fonction main: cette observation n'est pas dans le domaine npex=%d, npey=%d "
		     "ilonband=%d jlatband=%d lat=%f lon=%f\n", opt.npex, opt.npey, ilonband, jlatband, lat, lon);
	    continue;
	  }

          if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0)
            if (ilonband != opt.cherrypick_x || jlatband != opt.cherrypick_y)
              continue;

	  /* Il faut indiquer dans quel bande cette observation se trouve */
	  id = (ilonband-1)*opt.npey+(jlatband-1);
	  if (VERBOSE>4)
	    printf("Fonction main: cette observation est acceptee: ilonband=%d jlonband=%d nts[%d]=%d\n",
		   ilonband, jlatband, (ilonband-1)*opt.npey+(jlatband-1), nts[(ilonband-1)*opt.npey+(jlatband-1)]);
	}

	nts[id] = 1;
	num_obs_per_tile[id]++; /* On ajoute 1 au compteur total des observations */

	status = fill_rptout_blk(rptin,rptout,nts,(int*) NULL,opt.npex*opt.npey,opt.cherrypick_x,opt.cherrypick_y,opt.npey);
	if (status<0) {
	  fprintf(stderr,"Fonction main: Erreur %d dans la fonction fill_rptout_blk pour "
		  "l'adresse %d (%d rapport)\n", status, adresses[i_enrgs], i_enrgs);
	  EXIT_STATUS = 1;
	  break;
	}
      } /* Fin du else ... */

      if (t_in_domain != (int*) NULL) {
	free(t_in_domain);
	t_in_domain = (int*) NULL;
      }

      if (EXIT_STATUS!=0) {
	fprintf(stderr,"Fonction main: il y a eu une erreur precedemment (dans else if (vertical_clipping==1))\n");
	break;
      }

      /* Si c'est un enregistrement alors on active le compteur pour
	 le nombre d'enregistrements d'observations traites jusqu'a
	 maintenant */
      if (engrs_resume==0) i_obs_enrgs++;

      for (i=0;i<opt.npex*opt.npey;i++) {
	/* Si c'est un enregistrement resume (engrs_resume==1) alors on l'ecrit dans tous les fichiers */
	/* Si ce n'est pas un enregistrement resume (engrs_resume==0) alors
	 * on n'ecrit que dans les fichiers qui contiennent des observations */
	if (VERBOSE>5)
	  printf("Fonction main (juste avant 'brp_writerpt'): engrs_resume=%d nts[%d]=%d\n",
		 engrs_resume, i, nts[i]);

	if (engrs_resume==0 && nts[i]==0)
	  continue;

        /* Si on est en mode 'cherrypick', alors on ne considere que si la tuile est egale a celle voulue */
        if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0) {
          int cherrypick_id = (opt.cherrypick_x-1)*opt.npey+opt.cherrypick_y-1;
          if (i != cherrypick_id)
            continue;
        }

	if (VERBOSE>5) {
	  BURP_BLK* blk = (BURP_BLK*) NULL;

	  printf ( "Fonction main (juste avant 'brp_writerpt'): Entete du rapport rptout[%d]\n",i) ;
	  printf ( "hhmm   =%8d " , RPT_TEMPS(rptout[i])) ;
	  printf ( "flgs   =%6d  ", RPT_FLGS(rptout[i])) ;
	  printf ( "codtyp =%6d  ", RPT_IDTYP(rptout[i])) ;
	  printf ( "stnids =%9s\n", RPT_STNID(rptout[i])) ;
	  printf ( "blat   =%8d " , RPT_LATI(rptout[i])) ;
	  printf ( "blon   =%6d  ", RPT_LONG(rptout[i])) ;
	  printf ( "dx     =%6d  ", RPT_DX(rptout[i])) ;
	  printf ( "dy     =%6d  ", RPT_DY(rptout[i])) ;
	  printf ( "stnhgt =%6d\n", RPT_ELEV(rptout[i])) ;
	  printf ( "yymmdd =%8d " , RPT_DATE(rptout[i])) ;
	  printf ( "oars   =%6d  ", RPT_OARS(rptout[i])) ;
	  printf ( "runn   =%6d  ", RPT_RUNN(rptout[i])) ;
	  printf ( "nblk   =%6d  ", RPT_NBLK(rptout[i])) ;
	  printf ( "dlay   =%6d\n", RPT_DRND(rptout[i])) ;

	  blk = brp_newblk();
	  BLK_SetBKNO(blk, 0);
	  while ( brp_findblk( blk, rptout[i] ) >= 0 ) {
	    int thisi, thisj, thisk;
	    BURP_BLK* blkout = (BURP_BLK*) NULL;

	    blkout = brp_newblk();
	    status = brp_readblk(BLK_BKNO(blk), blkout, rptout[i], 0);

	    printf ( "blkno  =%6d  ", BLK_BKNO(blkout)    ) ;
	    printf ( "nele   =%6d  ", BLK_NELE(blkout)    ) ;
	    printf ( "nval   =%6d  ", BLK_NVAL(blkout)    ) ;
	    printf ( "nt     =%6d  ", BLK_NT(blkout)      ) ;
	    printf ( "bit0   =%6d\n", BLK_BIT0(blkout)    ) ;
	    printf ( "bdesc  =%6d  ", BLK_BDESC(blkout)   ) ;
	    printf ( "btyp   =%6d  ", BLK_BTYP(blkout)    ) ;
	    printf ( "nbit   =%6d  ", BLK_NBIT(blkout)    ) ;
	    printf ( "datyp  =%6d  ", BLK_DATYP(blkout)   ) ;
	    printf ( "bfam   =%6d  ", BLK_BFAM(blkout)    ) ;

	    for ( thisk = 0 ; thisk < BLK_NT(blkout) ; thisk++ ) {
	      printf (  "\nlstele =" ) ;
	      for (thisi=0;thisi<BLK_NELE(blkout);thisi++)
		printf (  "    %6.6d", BLK_DLSTELE(blkout,thisi) ) ;
	      /* sortie des valeurs des elements */
	      for (  thisj = 0 ; thisj < BLK_NVAL(blkout) ; thisj++ ) {
		printf (  "\ntblval =" ) ;
		for (  thisi = 0 ; thisi < BLK_NELE(blkout) ; thisi++ )
		  printf (  "%10d", BLK_TBLVAL(blkout,thisi,thisj,thisk) ) ;
	      }
	    }
	    printf("\n");
	    brp_freeblk(blkout);
	  } /* Fin du while ( brp_findblk( blk, rptout[i] ) >= 0 ) */
	  brp_freeblk(blk);
	} /* Fin du if (VERBOSE>5) */

	if (VERBOSE>4)
	  printf("Fonction main: appel de 'brp_writerpt' pour i=%d, iouts=%d et i_enrgs=%d\n",i,iouts[i],i_enrgs);

	status = brp_writerpt(iouts[i],rptout[i],END_BURP_FILE);
	if (status<0) {
	  fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_writerpt pour "
		  "le fichier de sortie %s_%d_%d a l'adresse %d (%d rapport)\n",
		  status, opt.obsout, i/opt.npey+1, i%opt.npey+1, adresses[i_enrgs], i_enrgs);
	  EXIT_STATUS = 1;
	  break;
	}
      } /* Fin du for (i=0;i<opt.npex*opt.npey;i++) */
    } /* Fin du for (i_enrgs=0;i_enrgs<nombre_enregistrements;i_enrgs++) */

    /* fermeture de fichier burp d'entree */
    status = brp_close(iun);
    if (status<0) {
      fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close pour le fichier %s\n", status, opt.obsin);
      EXIT_STATUS = 1;
    }

    /* fermeture de fichier burp de sortie */
    if ( opt.npex == 1 && opt.npey == 1 ) {
      if (VERBOSE>2)
	  printf("\nClosing BURP file %s iout = %d", opt.obsout, iout);

      printf("\nIl y a %d observations qui ont ete selectionnees\n", num_obs_per_tile[0]);

      status = brp_close(iout);
      if (status<0) {
	fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close pour le fichier %s\n", status, opt.obsout);
	EXIT_STATUS = 1;
      }
      else if (num_obs_per_tile[0]==0) {
	printf("Aucune observation n'a ete garde alors on efface le fichier %s\n", opt.obsout);
	status = remove(opt.obsout);
	if (status!=0) {
	  fprintf(stderr,"Il est impossible d'effacer le fichier %s\n", opt.obsout);
	  EXIT_STATUS = 1;
	}
	else if (VERBOSE>2)
	  printf("On a efface le fichier BURP %s\n", opt.obsout);
      }
      else { /* On imprime le nombre de headers presents dans le domaine */
	FILE* file;
	char burpout_num_headers[MAXSTR];

	sprintf(burpout_num_headers,"%s.num_headers", opt.obsout);
	
	status = access(burpout_num_headers,F_OK);
	if (status==0)
	  fprintf(stderr,"Attention le fichier '%s' sera efface\n", burpout_num_headers);

	file = (FILE*) fopen(burpout_num_headers,"w");
	fprintf(file,"%d\n", num_obs_per_tile[0]);
        fclose(file);
      }
    }
    else {
      char npex_str[MAXSTR], npey_str[MAXSTR], format_digits[MAXSTR], burpout[MAXSTR];
      int max_num_headers=0;

      sprintf(format_digits,"%%.%dd",opt.ndigits);

      for (ilonband=0;ilonband<opt.npex;ilonband++)
	for (jlatband=0;jlatband<opt.npey;jlatband++) {
	  int id=ilonband*opt.npey+jlatband;

          if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0)
            if (ilonband != opt.cherrypick_x-1 || jlatband != opt.cherrypick_y-1)
              continue;

	  sprintf(npex_str,format_digits,ilonband+1);
	  sprintf(npey_str,format_digits,jlatband+1);
	  sprintf(burpout,"%s_%s_%s", opt.obsout,npex_str,npey_str);

	  if (VERBOSE>2)
	    printf("\nClosing correctly BURP file %s iouts[%d] = %d", burpout, id, iouts[id]);

	  status = brp_close(iouts[id]);
	  if ( status<0 )
	    fprintf(stderr,"Fonction main: Erreur %d dans la fonction brp_close "
		    "pour le fichier %s\n", status, burpout);

	  printf("\nIl y a %d observations qui ont ete selectionnees et mise dans le fichier %s\n",
		 num_obs_per_tile[id], burpout);

	  if (num_obs_per_tile[id]==0) {
	    printf("Aucune observation n'a ete garde alors on efface le fichier %s\n", burpout);
	    status = remove(burpout);
	    if (status!=0) {
	      fprintf(stderr,"Il est impossible d'effacer le fichier %s\n", burpout);
	      EXIT_STATUS = 1;
	    }
	    else if (VERBOSE>2)
	      printf("On a efface le fichier BURP %s\n", burpout);
	  }
	  else { /* On imprime le nombre de headers presents dans la tuile id */
	    FILE* file;
	    char burpout_num_headers[MAXSTR];

	    sprintf(burpout_num_headers,"%s.num_headers", burpout);

	    status = access(burpout_num_headers,F_OK);
	    if (status==0)
	      fprintf(stderr,"Attention le fichier '%s' sera efface\n", burpout_num_headers);

	    file = (FILE*) fopen(burpout_num_headers,"w");
	    fprintf(file,"%d\n", num_obs_per_tile[id]);
            fclose(file);
	  }
	  if (num_obs_per_tile[id]>max_num_headers) max_num_headers=num_obs_per_tile[id];
	}

      if (max_num_headers>0) { /* On imprime le nombre maximal de headers */
	FILE* file;
	char burpout_max_num_headers[MAXSTR];

	sprintf(burpout_max_num_headers,"%s.max_num_headers", opt.obsout);

	status = access(burpout_max_num_headers,F_OK);
	if (status==0)
	  fprintf(stderr,"Attention le fichier '%s' sera efface\n", burpout_max_num_headers);

	file = (FILE*) fopen(burpout_max_num_headers,"w");
	fprintf(file,"%d\n", max_num_headers);
        fclose(file);
      }
      else
	printf("Il n'y a aucune observation qui a ete acceptee\n");
    } /* Fin du else associe au if ( opt.npex == 1 && opt.npey == 1 ) */

    brp_freerpt(rptin);
    if (rptout != (BURP_RPT**) NULL) {
      for (i=0;i<opt.npex*opt.npey;i++) {
        /* Si on est en mode 'cherrypick', alors on ne considere que si la tuile est egale a celle voulue */
        if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0) {
          int cherrypick_id = (opt.cherrypick_x-1)*opt.npey+opt.cherrypick_y-1;
          if (i != cherrypick_id)
            continue;
        }
	brp_freerpt(rptout[i]);
      }
      free(rptout);
    }
    free(nts);
    free(num_obs_per_tile);
    free(adresses);
    free(iouts);
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
      fprintf(stderr,"Fonction main: Le mode 'round-robin' n'a pas encore ete "
              "implementee pour des fichiers d'input de type ASCII!\n");

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

    filein=fopen(opt.obsin,"r");
    if (filein == (FILE*) NULL) {
      fprintf(stderr,"Fonction main: Le fichier ascii %s n'a pu etre ouvert correctement!\n", opt.obsin);

      status = c_gdrls(grid.gridid);
      if (status<0)
	fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

    fileout=fopen(opt.obsout,"w");
    if (fileout == (FILE*) NULL) {
      fprintf(stderr,"Fonction main: Le fichier ascii %s n'a pu etre ouvert correctement!\n", opt.obsout);

      fclose(filein);

      status = c_gdrls(grid.gridid);
      if (status<0)
	fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

      if (strlen(opt.gz)!=0) {
	status = c_gdrls(grid_gz.gridid);
	if (status<0)
	  fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

	if (VALEURS_GZ_MIN != (float*) NULL)
	  free(VALEURS_GZ_MIN);
	if (VALEURS_GZ_MAX != (float*) NULL)
	  free(VALEURS_GZ_MAX);
      }

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

#define REGEX_DEFINITION  "^[[:blank:]]*([-]?[0-9]*\\.?[0-9]*)[[:blank:]]+([-]?[0-9]*\\.?[0-9]*)[[:blank:]]+([-]?[0-9]*\\.?[0-9]*).*$"

    regex_err = regcomp(&regex, REGEX_DEFINITION, REG_EXTENDED);
    if (regex_err!=0) {
      regerror(regex_err, &regex, regex_errbuf, MAXSTR);
      fprintf(stderr,"Fonction main: Erreur avec la fonction regcomp '%s' avec l'expression reguliere '%s'\n",regex_errbuf, REGEX_DEFINITION);

      fclose(filein);
      fclose(fileout);

      status = c_gdrls(grid.gridid);
      if (status<0)
	fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);

      if (strlen(opt.gz)!=0) {
	status = c_gdrls(grid_gz.gridid);
	if (status<0)
	  fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);

	if (VALEURS_GZ_MIN != (float*) NULL)
	  free(VALEURS_GZ_MIN);
	if (VALEURS_GZ_MAX != (float*) NULL)
	  free(VALEURS_GZ_MAX);
      }

      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }

    while( fgets(ligne,sizeof(ligne),filein) != (char*) NULL ) {

      regex_err = regexec(&regex,ligne,regex.re_nsub+1,regex_match,0);
      if (regex_err!=0) {
	regerror(regex_err, &regex, regex_errbuf, MAXSTR);
	fprintf(stderr,"Fonction main: Erreur avec la fonction regexec '%s' sur la ligne '%s'\n",regex_errbuf,ligne);
	EXIT_STATUS = 1;
	break;
      }

      strncpy(latstr,&ligne[regex_match[1].rm_so],regex_match[1].rm_eo-regex_match[1].rm_so);
      strncpy(lonstr,&ligne[regex_match[2].rm_so],regex_match[2].rm_eo-regex_match[2].rm_so);
      strncpy(altstr,&ligne[regex_match[3].rm_so],regex_match[3].rm_eo-regex_match[3].rm_so);

      lat = atof(latstr);
      lon = atof(lonstr);
      alt = atof(altstr);

      status = checkgrid(grid.gridid, grid.ni, grid.nj, lat, lon, opt.rect, errmsg);
      if (status<0) {
	fprintf(stderr,"Fonction main: Erreur dans la fonction checkgrid pour la ligne '%s' avec le message '%s'\n", ligne, errmsg);
	EXIT_STATUS = 1;
	break;
      }

      if (opt.inout == status) {
	if (strlen(opt.channels)==0 && opt.niveau_min == IP1_VIDE && opt.niveau_max == IP1_VIDE) 
	  /* Aucun filtrage vertical n'est fait */
	  fputs(ligne,fileout);
	else if (strlen(opt.channels)==0 && strlen(opt.gz)==0 && checkvertical(alt,opt.niveau_min,opt.niveau_max))
	  /* Le filtrage vertical est fait a l'aide d'une hauteur en pression */
	  fputs(ligne,fileout);
	else if (strlen(opt.channels)==0 && checkvertical_gz(lat,lon,alt,grid_gz.gridid,grid_gz.ni,grid_gz.nj,opt.niveau_min,opt.niveau_max))
	  /* Le filtrage vertical est fait a l'aide d'une hauteur en metre */
	  fputs(ligne,fileout);
	else if (opt.channels_voulus==checkcanal(alt,opt.channels))
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
      fprintf(stderr,"Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid.gridid);
      EXIT_STATUS = 1;
    }

    if (strlen(opt.gz)!=0) {
      status = c_gdrls(grid_gz.gridid);
      if (status<0) {
	fprintf(stderr, "Fonction main: Erreur dans la fonction c_gdrls pour gridid = %d\n", grid_gz.gridid);
	EXIT_STATUS = NOT_OK;
      }

      if (VALEURS_GZ_MIN != (float*) NULL)
	free(VALEURS_GZ_MIN);
      if (VALEURS_GZ_MAX != (float*) NULL)
	free(VALEURS_GZ_MAX);
    }
  } /* Fin du 'if ( opt.roundrobin == 0 )'*/

  /* On a termine! */
  if (EXIT_STATUS != OK)
    exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);

  exit_program(OK,PROGRAM_NAME,"OK",VERSION);

  return OK;

} /* Fin du main */


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
   * fonction: sqlite_add_resume_request
   *
   * Genere une requete SQL pour copier les tables resume.
   *
   ***************************************************************************/
int sqlite_add_resume_request(char* obsin, char* requete_sql, char* attached_db_name) {
  int status, is_resume_and_rdb4_schema_present_in_DB;
  char *ErrMsg, sqlreqtmp[MAXSTR];
  sqlite3 *sqldbin;

  /* Cette partie sert a trouver la requete pour copier les tables 'resume' et 'rdb4_schema' */
  /* On ouvre le fichier d'input */
  status = sqlite3_open(obsin,&sqldbin);
  if ( status != SQLITE_OK ) {
    fprintf(stderr, "Fonction sqlite_add_resume_request: Incapable d'ouvrir le fichier '%s' avec l'erreur '%s'\n", obsin, sqlite3_errmsg(sqldbin));

    status = sqlite3_close(sqldbin);
    if( status != SQLITE_OK )
      fprintf(stderr,"Fonction sqlite_add_resume_request: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, obsin);

    return NOT_OK;
  } /* Fin du 'if ( status != SQLITE_OK )' */

  /* Execution de la requete SQL sur la base de donnees */
  /* L'idee est de detecter la presence des tables de resume pour les copier.  */
  is_resume_and_rdb4_schema_present_in_DB = 0;
  status = sqlite3_exec(sqldbin, "select * from sqlite_master", sqlite_check_resume_callback, &is_resume_and_rdb4_schema_present_in_DB, &ErrMsg);
  if( status != SQLITE_OK ) {
    fprintf(stderr, "Fonction sqlite_add_resume_request: Erreur %d pour le fichier dans la fonction sqlite3_exec: %s\n", status, ErrMsg);
    sqlite3_free(ErrMsg);

    return NOT_OK;
  } /* Fin du 'if ( status != SQLITE_OK )' */

  if (is_resume_and_rdb4_schema_present_in_DB != 0 && is_resume_and_rdb4_schema_present_in_DB != 10 &&
      is_resume_and_rdb4_schema_present_in_DB != 1 && is_resume_and_rdb4_schema_present_in_DB != 11) {
    fprintf(stderr, "Fonction sqlite_add_resume_request: On attend 0, 1, 10 ou 11 comme valeur pour 'is_resume_and_rdb4_schema_present_in_DB' et non %d\n", is_resume_and_rdb4_schema_present_in_DB);

    status = sqlite3_close(sqldbin);
    if( status != SQLITE_OK )
      fprintf(stderr,"Fonction sqlite_add_resume_request: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, obsin);

    return NOT_OK;
  }

  strcpy(requete_sql,"");
  /* Si cette variable est egale a 1 ou 11, alors on a trouve une table 'resume' */
  if ( is_resume_and_rdb4_schema_present_in_DB == 1 || is_resume_and_rdb4_schema_present_in_DB == 11 ) {
    /* La table resume existe alors on fait la commande */
    sprintf(sqlreqtmp,"\ninsert into rdb4_schema select * from %s.rdb4_schema;", attached_db_name);
    strcat(requete_sql,sqlreqtmp);
  }

  /* Si cette variable est egale a 10 ou 11, alors on a trouve une table 'rdb4_schema' */
  if ( is_resume_and_rdb4_schema_present_in_DB == 10 || is_resume_and_rdb4_schema_present_in_DB == 11 ) {
    /* La table resume existe alors on fait la commande */
    sprintf(sqlreqtmp,"\ninsert into resume select * from %s.resume;", attached_db_name);
    strcat(requete_sql,sqlreqtmp);
  }

  status = sqlite3_close(sqldbin);
  if( status != SQLITE_OK ) {
    fprintf(stderr,"Fonction sqlite_add_resume_request: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, obsin);

    return NOT_OK;
  }

  return OK;
}


  /***************************************************************************
   * fonction: sqlite_check_resume_callback
   *
   * Cette fonction sert de 'callback' pour la requete executee par 'sqlite3_exec'.
   *    On veut savoir si elle contient les tables 'rdb4_schema' et 'resume'
   *
   ***************************************************************************/
static int sqlite_check_resume_callback(void *is_resume_and_rdb4_schema_present_in_DB_void, int count, char **data, char **columns) {
  int idx, is_resume_and_rdb4_schema_present_in_DB, *is_resume_and_rdb4_schema_present_in_DB_ptr;

  is_resume_and_rdb4_schema_present_in_DB_ptr = (int*) is_resume_and_rdb4_schema_present_in_DB_void;

  // printf("Dans la fonction 'sqlite_check_resume_callback'\n");
  // printf("There are %d column(s)\n", count);

  is_resume_and_rdb4_schema_present_in_DB = 0;
  for (idx = 0; idx < count; idx++) {
    // printf("The data in column \"%s\" is: '%s'\n", columns[idx], data[idx]);
    if (strcmp(columns[idx],"tbl_name")==0)
      if(strcasecmp(data[idx],"rdb4_schema")==0) {
        is_resume_and_rdb4_schema_present_in_DB++;
        // printf("Found 'rdb4schema'\n");
      }
      else if(strcasecmp(data[idx],"resume")==0) {
        is_resume_and_rdb4_schema_present_in_DB += 10;
        // printf("Found 'resume'\n");
      }
  }
  // printf("sqlite_check_resume_callback: is_resume_and_rdb4_schema_present_in_DB = %d\n", is_resume_and_rdb4_schema_present_in_DB);
  *is_resume_and_rdb4_schema_present_in_DB_ptr += is_resume_and_rdb4_schema_present_in_DB;
  // printf("sqlite_check_resume_callback: is_resume_and_rdb4_schema_present_in_DB_ptr = %d\n", *is_resume_and_rdb4_schema_present_in_DB_ptr);

  return 0;
}


  /***************************************************************************
   * fonction: sqlite_get_tables_with_id_obs
   *
   * Trouve les tables qui contiennent une colonne 'id_obs'
   *
   ***************************************************************************/
int sqlite_get_tables_with_id_obs(char* obsin, char* table_list) {
  int status;
  char *ErrMsg;
  sqlite3 *sqldbin;

  /* Cette partie sert a trouver la requete pour copier les tables 'resume' et 'rdb4_schema' */
  /* On ouvre le fichier d'input */
  status = sqlite3_open(obsin,&sqldbin);
  if ( status != SQLITE_OK ) {
    fprintf(stderr, "Fonction sqlite_get_tables_with_id_obs: Incapable d'ouvrir le fichier '%s' avec l'erreur '%s'\n", obsin, sqlite3_errmsg(sqldbin));

    status = sqlite3_close(sqldbin);
    if( status != SQLITE_OK )
      fprintf(stderr,"Fonction sqlite_get_tables_with_id_obs: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, obsin);

    return NOT_OK;
  } /* Fin du 'if ( status != SQLITE_OK )' */

  strcpy(table_list,"");

  /* Execution de la requete SQL sur la base de donnees */
  /* Since the table 'header' and 'data' are already processed in the main */
  /*   request, we must not process them again when adding tables which */
  /*   contains a column 'id_obs'. */
  status = sqlite3_exec(sqldbin, "select * from sqlite_master where lower(name) not in ('header','data');", sqlite_check_tables_with_id_obs_callback, table_list, &ErrMsg);
  if( status != SQLITE_OK ) {
    fprintf(stderr, "Fonction sqlite_get_tables_with_id_obs: Erreur %d pour le fichier dans la fonction sqlite3_exec: %s\n", status, ErrMsg);
    sqlite3_free(ErrMsg);

    return NOT_OK;
  } /* Fin du 'if ( status != SQLITE_OK )' */


  status = sqlite3_close(sqldbin);
  if( status != SQLITE_OK ) {
    fprintf(stderr,"Fonction sqlite_get_tables_with_id_obs: Erreur %d de la fonction sqlite3_close pour le fichier '%s'\n", status, obsin);

    return NOT_OK;
  }

  return OK;
}


  /***************************************************************************
   * fonction: sqlite_check_tables_with_id_obs_callback
   *
   * Cette fonction sert de 'callback' pour la requete executee par 'sqlite3_exec'.
   *    On veut savoir si elle contient les tables 'rdb4_schema' et 'resume'
   *
   ***************************************************************************/
static int sqlite_check_tables_with_id_obs_callback(void *table_list, int count, char **data, char **columns) {
  int idx, isTypeTable, foundID_OBS;
  char table_name[MAXSTR];

  isTypeTable = 0;
  for (idx = 0; idx < count; idx++) {
    if (strcasecmp(columns[idx],"type")==0)
      if (strcasecmp(data[idx],"table")==0)
        isTypeTable = 1;
  }

  if (isTypeTable == 0) return 0;

  strcpy(table_name,"");
  foundID_OBS=0;
  for (idx = 0; idx < count; idx++) {
    if (strcasecmp(columns[idx],"tbl_name")==0)
      strcpy(table_name,data[idx]);
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

      regex_t regex;
      char regex_errbuf[MAXSTR];
      int regex_err;

      regex_err = regcomp(&regex, "id_obs", REG_ICASE);
      if (regex_err!=0) {
	regerror(regex_err, &regex, regex_errbuf, MAXSTR);
        fprintf(stderr,"sqlite_check_tables_with_id_obs_callback: cannot compile regular expression '%s': error '%s'", "id_obs", regex_errbuf);
        return 1;
      }
      regex_err = regexec(&regex,data[idx],0,(regmatch_t*) NULL,0);
      if (regex_err == 0) {
        // This means that there was a match
        foundID_OBS=1;
      }
      else if (regex_err == REG_NOMATCH) {
        foundID_OBS=0;
      }
      else {
	regerror(regex_err, &regex, regex_errbuf, MAXSTR);
        regfree(&regex);
	fprintf(stderr,"sqlite_check_tables_with_id_obs_callback: Erreur '%s' avec la fonction regexec sur la ligne '%s'\n", regex_errbuf, data[idx]);
        return 1;
      }

      regfree(&regex);
    }
  }

  if (foundID_OBS) {
    if (strlen(table_name)>0) {
      if (strlen((char*) table_list)>0) {
        strcat((char*) table_list, " ");
      }
      strcat((char*) table_list, table_name);
    }
    else {
      fprintf(stderr,"sqlite_check_tables_with_id_obs_callback: foundID_OBS = %d but table_name is empty", foundID_OBS);
      return 1;
    }
  }

  return 0;
} /* End of function 'sqlite_check_tables_with_id_obs_callback' */


  /***************************************************************************
   * fonction: append_id_obs_table_list_requests
   *
   * Cette fonction sert a ajouter des requetes SQL pour inclure les
   * tables supplementaires autres que 'header' et 'data' mais qui ont
   * 'id_obs' comme colonne.
   *
   ***************************************************************************/
void append_id_obs_table_list_requests(char* requete_sql, char* table_list) {
  const char separator_char[2] = " ";
  char sqlreqtmp[MAXSTR], table_list_tmp[MAXSTR];
  char *token;

  // Make a copy of 'table_list' input string because 'strtok' is changing in place that string
  strcpy(table_list_tmp,table_list);
  /* get the first token */
  token = strtok(table_list_tmp, separator_char);
  /* walk through other tokens */
  while( token != (char*) NULL ) {
    sprintf(sqlreqtmp,"insert into %s select * from dbin.%s where dbin.%s.id_obs in (select id_obs from header);\n",token,token,token);
    strcat(requete_sql,sqlreqtmp);
    token = strtok((char*) NULL, separator_char);
  }
  // On a termine d'ajouter les requetes pour les autres tables
}


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
int getGZ(int iun, char* fichier, gridtype* gridptr, int niveau, float** valeurs) {
  int ier, key, status, datev;
  fstparam fst = fstparam_DEFAUT;
  double forecast;
  
  fst.ip1=niveau;
  strcpy(fst.nomvar,"GZ  ");

  status = open_stdfile(iun, fichier, "RND+R/O");
  if (status == NOT_OK) {
    fprintf(stderr, "Fonction getGZ: Erreur dans la fonction open_stdfile avec le fichier %s\n",fichier);

    return NOT_OK;
  }
  
  status = getgrid(iun,gridptr,&fst,fichier);
  if (status == NOT_OK) {
    fprintf(stderr, "Fonction getGZ: Erreur dans la fonction getgrid pour les parametres "
	    "(%s,%s,%s,%d,%d,%d,%d) dans le fichier %s\n",
	    fst.nomvar,fst.typvar,fst.etiket,fst.dateo,fst.ip1,fst.ip2,fst.ip3,fichier);

    /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
    close_stdfile(iun,fichier);
    
    return NOT_OK;
  }

  /**************************************************************
   * On lit maintenant les valeurs du champ GZ designe
   **************************************************************/
  /* Allocation de la memoire */
  *valeurs = (float*) malloc(gridptr->ni*gridptr->nj*sizeof(float));
  if ( *valeurs == (float*) NULL) {
    fprintf(stderr, "Fonction getGZ: Incapable d'allouer un vecteur de float de dimension %dx%d=%d elements\n",
	    gridptr->ni,gridptr->nj,gridptr->ni*gridptr->nj);

    status = c_gdrls(gridptr->gridid);
    if (status<0)
      fprintf(stderr, "Fonction getGZ: Erreur dans la fonction c_gdrls pour gridid = %d\n", gridptr->gridid);

    /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
    close_stdfile(iun,fichier);

    return NOT_OK;
  }

  /* On doit incrementer l'heure de recherche pour la date de validite en fonction de fst.dateo et fst.npas, fst.deet */
  forecast = ((double) fst.npas*fst.deet)/3600; /* prevision en heures */
  f77name(incdatr)(&datev,&fst.dateo,&forecast);

  key = c_fstinf(iun,&fst.ni,&fst.nj,&fst.nk,datev,fst.etiket,
		 fst.ip1,fst.ip2,fst.ip3,fst.typvar,fst.nomvar);
  if (key<0) {
    fprintf(stderr,"Fonction getGZ: Erreur %d avec le fichier %s pour les parametres (%s,%s,%s,%d,%d,%d,%d) dans la fonction c_fstinf\n",
	    key,fichier,fst.nomvar,fst.typvar,fst.etiket,fst.dateo,fst.ip1,fst.ip2,fst.ip3);

    status = c_gdrls(gridptr->gridid);
    if (status<0)
      fprintf(stderr, "Fonction getGZ: Erreur dans la fonction c_gdrls pour gridid = %d\n", gridptr->gridid);

    /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
    close_stdfile(iun,fichier);

    return NOT_OK;
  }

  ier = c_fstluk(*valeurs,key,&gridptr->ni,&gridptr->nj,&gridptr->nk);
  if (ier<0) {
    fprintf(stderr, "Fonction getGZ: Erreur %d avec le fichier %s pour les parametres "
	    "(%s,%s,%s,%d,%d,%d,%d) dans la fonction fstluk\n",
	    ier,fichier,fst.nomvar,fst.typvar,fst.etiket,fst.dateo,fst.ip1,fst.ip2,fst.ip3);
    
    status = c_gdrls(gridptr->gridid);
    if (status<0)
      fprintf(stderr, "Fonction getGZ: Erreur dans la fonction c_gdrls pour gridid = %d\n", gridptr->gridid);

    /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
    close_stdfile(iun,fichier);

    free(*valeurs);
    
    return NOT_OK;    
  }

  /* On ferme le fichier standard ouvert pour lire le champ definissant la grille */
  status = close_stdfile(iun,fichier);
  if (status == NOT_OK) {
    fprintf(stderr, "Fonction getGZ: Erreur avec la fonction close_stdfile pour le fichier %s\n",fichier);

    status = c_gdrls(gridptr->gridid);
    if (status<0)
      fprintf(stderr, "Fonction getGZ: Erreur dans la fonction c_gdrls pour gridid = %d\n", gridptr->gridid);

    free(*valeurs);
    
    return NOT_OK;
  }
    
  return OK;
} /* Fin de la fonction getGZ() */


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
int checkgrid(int gridid, int ni, int nj, float lat, float lon, rectangle rect, char errmsg[MAXSTR]) {
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
    sprintf(errmsg, "Fonction checkgrid: Erreur avec c_gdxyfll qui retourne %d "
	    "pour lat = %f, lon = %f, ni = %d, nj = %d, gridid = %d\n", 
	    status, lat, lon, gridid, ni, nj);
    fprintf(stderr,"%s",errmsg);
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
      fprintf(stderr,"Fonction checkgrid: Le critere '%d' n'est pas possible.  Il faut 1,2, 3 ou 4.  ", criteria);

  }

  return 0;

} /* Fin de la fonction checkgrid */


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
    if (VERBOSE>1) {
      printf("debug: id_obs=%d vcoord=NULL niveau_max=%d niveau_min=%d -> ",id_obs,niveau_max,niveau_min);
      printf("Obs acceptee parce qu'a la surface\n");
    }
    sqlite3_result_int(context, 1);
    return;
  }

  if (VERBOSE>1) {
    printf("debug: id_obs=%d vcoord=%f niveau_max=%d niveau_min=%d -> ",id_obs,vcoord,niveau_max,niveau_min);
  }
  
  status = checkvertical(vcoord,niveau_min,niveau_max);

  if (status<0) {
    char errmsg[MAXSTR];
    sprintf(errmsg, "Fonction checkvertical_sql: Erreur avec checkvertical pour "
	    "id_obs=%d vcoord=%f niveau_min=%d et niveau_max=%d\n", id_obs, vcoord, niveau_min, niveau_max);
    fprintf(stderr,"%s",errmsg);
    sqlite3_result_error(context, errmsg, -1);
    return;
  }

  sqlite3_result_int(context, status);
  return;

} /* Fin de la fonction checkvertical_sql */


  /***************************************************************************
   * fonction: checkvertical
   *
   * Cette fonction sert a verifier si le point (lat,lon,vcoord) est sous une certaine hauteur en hPa.
   *     vcoord: hauteur de l'observation
   *     niveau_min: niveau minimum acceptable (en hPa)
   *     niveau_max: niveau maximum acceptable (en hPa)
   *
   ***************************************************************************/
int checkvertical(float vcoord, int niveau_min, int niveau_max) {
  
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
      fprintf(stderr, "Fonction checkvertical: Erreur pour niveau_min=%d et niveau_max=%d\n", niveau_min, niveau_max);
      return -1;
    }
  } /* Fin du else du if if (niveau_min == IP1_VIDE && niveau_max == IP1_VIDE ) */
} /* Fin de la fonction check_vertical */


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
    if (VERBOSE>2) {
      printf("debug: id_obs=%d vcoord=NULL niveau_max=%d niveau_min=%d -> ",id_obs,niveau_max,niveau_min);
      printf("Obs acceptee parce qu'a la surface\n");
    }
    sqlite3_result_int(context, 1);
    return;
  }

  if (VERBOSE>2) {
    printf("debug: id_obs=%d lat=%f lon=%f vcoord=%f niveau_max=%d niveau_min=%d -> ",id_obs,lat,lon,vcoord,niveau_max,niveau_min);
  }
  
  status = checkvertical_gz(lat,lon,vcoord,gridid,ni,nj,niveau_min,niveau_max);

  if (status<0) {
    char errmsg[MAXSTR];
    sprintf(errmsg, "Fonction checkvertical_gz_sql: Erreur avec checkvertical_gz pour "
	    "id_obs=%d lat=%f lon=%f vcoord=%f niveau_min=%d niveau_max=%d\n",id_obs,lat,lon,vcoord,niveau_min,niveau_max);
    fprintf(stderr,"%s",errmsg);
    sqlite3_result_error(context, errmsg, -1);
    return;
  }

  sqlite3_result_int(context, status);
  return;
  
} /* Fin de la fonction checkvertical_gz_sql */


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
int checkvertical_gz(float lat, float lon, float vcoord, int gridid, int ni, int nj, int niveau_min, int niveau_max) {
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
	sprintf(errmsg, "Fonction checkvertical_gz: c_gdllsval retourne %d "
		"pour lat = %f, lon = %f, ni = %d, nj = %d, gridid = %d\n", 
		status, lat, lon, gridid, ni, nj);
	fprintf(stderr,"%s",errmsg);
	return -1;
      }

      status = c_gdllsval(gridid, &hauteur_max, VALEURS_GZ_MAX, &lat, &lon, 1);
      if (status<0) {
	char errmsg[MAXSTR];
	sprintf(errmsg, "Fonction checkvertical_gz: c_gdllsval retourne %d "
		"pour lat = %f, lon = %f, ni = %d, nj = %d, gridid = %d\n", 
		status, lat, lon, gridid, ni, nj);
	fprintf(stderr,"%s",errmsg);
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
	sprintf(errmsg, "Fonction checkvertical_gz: c_gdllsval retourne %d "
		"pour lat = %f, lon = %f, ni = %d, nj = %d, gridid = %d\n", 
		status, lat, lon, gridid, ni, nj);
	fprintf(stderr,"%s",errmsg);
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
	sprintf(errmsg, "Fonction checkvertical_gz: c_gdllsval retourne %d "
		"pour lat = %f, lon = %f, ni = %d, nj = %d, gridid = %d\n", 
		status, lat, lon, gridid, ni, nj);
	fprintf(stderr,"%s",errmsg);
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
      sprintf(errmsg, "Fonction checkvertical_gz: niveau_min=%d et niveau_max=%d\n", niveau_min, niveau_max);
      fprintf(stderr,"%s",errmsg);
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
int checkcanal(float canal, char* channels) {
  char* result;
  char  canalstr[MAXSTR];

  sprintf(canalstr,"%d", (int) canal);

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


  /***************************************************************************
   * fonction: find_subdomain
   *
   * Cette fonction sert a verifier si le point (lat,lon) est a l'interieur d'une grille
   * donnee par les arguments d'entree:
   *     gridid: identifiant de la grille EZSCINT
   *     ni: dimension horizontale de la grille
   *     nj: dimension verticale   de la grille
   *     lat: coordonnee de latitude  du point d'observation
   *     lon: coordonnee de longitude du point d'observation
   *     rect:  rectangle definissant une sous-region de la grille comme domaine
   *     npex: nombre de bandes de separation dans la direction 'x' ou 'i' (longitude)
   *     npey: nombre de bandes de separation dans la direction 'y' ou 'j' (latitude)
   *
   * Cette fonction retourne:
   *            1 si le point est a l'interieur de la grille
   *            0 si le point est a l'exterieur de la grille
   *           -1 s'il y a une erreur
   *
   ***************************************************************************/
int find_subdomain(int gridid, int ni, int nj, float lat, float lon, rectangle rect, int npex, int npey,
		   int* ilonband, int* jlatband, char errmsg[MAXSTR]) {
  int status, criteria;
  float x, y;

  if (VERBOSE>3)
    printf("Fonction find_subdomain: lat=%f  lon=%f npex=%d npey=%d\n", lat, lon, npex, npey);

  if (lon<0) lon+=360;

  status = c_gdxyfll(gridid, &x, &y, &lat, &lon, 1);
  if (status<0) {
    sprintf(errmsg, "Fonction find_subdomain: Erreur avec c_gdxyfll qui retourne %d "
	    "pour lat = %f, lon = %f, gridid = %d, ni = %d, nj = %d\n",
	    status, lat, lon, gridid, ni, nj);
    fprintf(stderr,"%s",errmsg);
    return -1;
  }

  if ( (!rect.min_i_equal && x>rect.min_i) || (rect.min_i_equal && x>=rect.min_i) ) {
    if ( (!rect.max_i_equal && x<rect.max_i) || (rect.max_i_equal && x<=rect.max_i) ) {
      if ( (!rect.min_j_equal && y>rect.min_j) || (rect.min_j_equal && y>=rect.min_j) ) {
	if ( (!rect.max_j_equal && y<rect.max_j) || (rect.max_j_equal && y<=rect.max_j) ) {
	  gridtype input_grid;

	  if (VERBOSE>3)
	    printf("Fonction find_subdomain: Obs dans le domaine: lat=%f  lon=%f x=%f y=%f "
		   "rect.min_i=%f rect.max_i=%f rect.min_j=%f "
		   "rect.max_j=%f\n",lat,lon,x,y,rect.min_i,rect.max_i,rect.min_j,rect.max_j);

	  status = c_ezgprm(gridid, &input_grid.grtyp, &input_grid.ni, &input_grid.nj,
			    &input_grid.ig1, &input_grid.ig2, &input_grid.ig3, &input_grid.ig4);
	  if (status<0) {
	    sprintf(errmsg, "Fonction find_subdomain: Erreur avec c_ezgprm qui retourne %d "
		    "pour gridid = %d, ni = %d, nj = %d\n", status, gridid, ni, nj);
	    fprintf(stderr,"%s",errmsg);
	    return -1;
	  }

	  if (x<1)
	    *ilonband = 1;
	  else if (x>=ni)
	    *ilonband = npex;
	  else {
	    float ilonband_float = (x-1)/(ni/npex)+1;
	    /* reproduce the logic in the OAVAR code routine 'setObsMpiStrategy' */
	    int ilonband_int = (((int) x) - 1)/(ni/npex)+1;

	    *ilonband = ilonband_int;
	  }

	  if (*ilonband<1 || *ilonband>npex)
	    fprintf(stderr,"Fonction find_subdomain: ilonband=%d n'est pas dans l'intervalle permis"
		    " pour lat=%f lon=%f x=%f y=%f npex=%d npey=%d ni=%d nj=%d\n",*ilonband,lat,lon,x,y,npex,npey,ni,nj);

	  if (y<1)
	    *jlatband = 1;
	  else if (y>=nj)
	    *jlatband = npey;
	  else {
	    float jlatband_float;
	    int jlatband_int;

	    /* On doit traiter le cas d'une grille gaussienne
	     * differemment parce que les points prets du pole sont
	     * soit <1 ou >nj
	     */
	    if (input_grid.grtyp[0]=='G') {
	      jlatband_float = y/(nj/npey)+1;
	      /* reproduce the logic in the OAVAR code routine 'setObsMpiStrategy' */
	      jlatband_int = ((int) y)/(nj/npey)+1;
	    }
	    else {
	      jlatband_float = (y-1)/(nj/npey)+1;
	      /* reproduce the logic in the OAVAR code routine 'setObsMpiStrategy' */
	      jlatband_int = (((int) y)-1)/(nj/npey)+1;
	    }

	    *jlatband = jlatband_int;
	  }

	  if (VERBOSE>3) {
	    if (input_grid.grtyp[0]=='G')
	      printf("Fonction find_subdomain: lat=%f lon=%f x=%f y=%f npex=%d npey=%d ilonband=%d jlatband=%d (x-1)/(ni/npex)+1=%f     y/(nj/npey)+1=%f\n",
		     lat,lon,x,y,npex,npey,*ilonband,*jlatband,(x-1)/(ni/npex)+1,y/(nj/npey)+1);
	    else
	      printf("Fonction find_subdomain: lat=%f lon=%f x=%f y=%f npex=%d npey=%d ilonband=%d jlatband=%d (x-1)/(ni/npex)+1=%f (y-1)/(nj/npey)+1=%f\n",
		     lat,lon,x,y,npex,npey,*ilonband,*jlatband,(x-1)/(ni/npex)+1,(y-1)/(nj/npey)+1);
	  }
	  if (*jlatband<1 || *jlatband>npey)
	    fprintf(stderr,"Fonction find_subdomain: jlatband=%d n'est pas dans l'intervalle permis"
		    " pour lat=%f lon=%f x=%f y=%f npex=%d npey=%d ni=%d nj=%d\n",*jlatband,lat,lon,x,y,npex,npey,ni,nj);

	  if (VERBOSE>3)
	    printf("Fonction find_subdomain: ilonband=%d jlatband=%d lat=%f lon=%f x=%f y=%f (90+lat)/(180./%d)=%f "
		   "rect.min_i=%f rect.max_i=%f rect.min_j=%f rect.max_j=%f\n",*ilonband,*jlatband,lat,lon,x,y,npey,
		   (90+lat)/(180./npey),rect.min_i,rect.max_i,rect.min_j,rect.max_j);

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

  if (VERBOSE>3) {
    printf("Fonction find_subdomain: Obs refusee: lat = %f  lon = %f x = %f y= %f "
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
      fprintf(stderr,"Fonction checkgrid: Le critere '%d' n'est pas possible.  Il faut 1,2, 3 ou 4.  ", criteria);

  }

  return 0;
} /* Fin de la fonction find_subdomain */


  /***************************************************************************
   * fonction: which_btyp
   *
   * En entree, cette fonction prend
   *    btyp: btyp d'un bloc
   *
   * En sortie:
   *    0 si c'est un bloc de donnees
   *    1 si c'est un bloc marqueurs
   *    2 si c'est un autre type de bloc (info, OmP, OmA, ...)
   *
   * Cette fonction retourne:
   *           -1 s'il y a une erreur
   *
   * Code python pour imprimer les bits d'un btyp et aider a trouver la formule
   *     btyps_amsub_deri=[5120,3072,9217,15361,9248,15392]
   *     btyps_ua4d_bgckalt=[106,15456,3107,6242,9312,9322,98]
   *
   *     print '\n'.join(["%6d %16s %d" % (btyp,bin(btyp>>11),btyp>>2 & 31==24) for btyp in btyps_amsub])
   ***************************************************************************/
int which_btyp(int btyp) {
  int crit;

  crit=btyp>>2 & 31;
  if (VERBOSE>3)
    printf("Fonction which_btyp: btyp=%d btyp>>2 & 31 = %d\n", btyp, crit);

  if (crit==0 || crit==8 || crit==24) {
    /* alors c'est un bloc de donnees ou un bloc marqueur */
    /* On verifie si c'est un bloc de data */
    crit=btyp>>11 & 3;
    if (VERBOSE>5)
      printf("Fonction which_btyp: btyp>>11 & 3 = %d\n", crit);

    if ( crit == 0 || crit == 2 )
      return 0;
    /* On verifie si c'est un bloc marqueur */
    else if ( crit == 3 )
      return 1;
    else {
      fprintf(stderr,"Fonction which_btyp: Ce n'est ni un bloc marqueur ni un "
	      "bloc de donnees btyp=%d (btyp>>2 & 31 = %d) (btyp>>11 & 3 = %d)\n",
	      btyp, btyp>>2 & 31, btyp>>11 & 3);
      return -1;
    }
  }

  /* Alors c'est un autre type de bloc */
  return 2;
} /* Fin de la fonction which_btyp */


  /***************************************************************************
   * fonction: btypAssociated
   *
   * En entree, cette fonction prend
   *    btyp_obs: btyp d'un bloc d'observations (data)
   *    btyp: un autre btyp pour lequel on va verifier s'il est associe au btyp_obs
   *
   * En sortie:
   *    0 si c'est le btyp n'est pas associe au btyp_obs
   *    1 si c'est btyp est associe au btyp_obs
   *    -1 s'il y a une erreur
   *
   * code Python equivalent:
   *      def is2split(btyp_obs, btyp):
   *         newbtyp_obs = btyp_obs >> 4 & 127;
   *         print 'newbtyp_obs = ', newbtyp_obs,
   *         newbtyp = btyp >> 4 & 127
   *         print ' newbtyp = ', newbtyp,
   *         if newbtyp_obs == newbtyp:
   *             print 'btyp_obs=%d btyp=%d et on splitte' % (btyp_obs,btyp)
   *         else:
   *             print 'btyp_obs=%d btyp=%d et on ne splitte pas' % (btyp_obs,btyp)
   *
   * Pour avoir la fonction 'btyp' qui aide a demeler tous ces aspects, on fait
   *     . ssmuse-sh -d /ssm/net/cmda/base/master
   *
   ***************************************************************************/
int btypAssociated(int btyp_obs, int btyp) {
  int newbtyp_obs, newbtyp, bknat;

  if (btyp_obs == btyp) {
    fprintf(stderr,"Fonction btypAssociated: erreur: btyp_obs = btyp = %d\n", btyp);
    return -1;
  }

  if (VERBOSE>3)
    printf("Fonction btypAssociated: btyp_obs=%d btyp=%d\n", btyp_obs, btyp);

  /* On verifie d'abord si c'est un bloc info */
  bknat = btyp >> 11;
  if (VERBOSE>3)
    printf("Fonction btypAssociated: btyp=%d bknat=%d\n", btyp, bknat);

  if (bknat==1) {
    if (VERBOSE>3)
      printf("Fonction btypAssociated: btyp=%d est un bloc info\n", btyp);

    return 0;
  }

  newbtyp_obs = btyp_obs >> 4 & 127;
  newbtyp = btyp >> 4 & 127;

  if (VERBOSE>3)
    printf("Fonction btypAssociated: btyp_obs>>4&127=%d btyp>>4&127=%d\n", newbtyp_obs, newbtyp);

  if (newbtyp == newbtyp_obs) {
    if (VERBOSE>3)
      printf("Fonction btypAssociated: btyp_obs=%d est associe au btyp=%d\n", btyp_obs, btyp);

    return 1;
  }

  if (VERBOSE>3)
    printf("Fonction btypAssociated: btyp_obs=%d n'est pas associe au btyp=%d\n", btyp_obs, btyp);

  return 0;

} /* Fin de la fonction btypAssociated */


  /***************************************************************************
   * fonction: clipping_vertical
   *
   * En entree, cette fonction prend
   *    rptin: un rapport BURP complet
   *    optptr: un pointeur a une structure 'options' qui permet d'aller chercher les criteres
   *    grid_gz: un pointeur a une structure 'gridtype' qui permet d'avoir l'information sur la grille 'GZ'
   *
   * En sortie:
   *    rptout: un rapport BURP sans les observations hors du domaine vertical
   *
   * Cette fonction retourne:
   *            0 si ce n'est pas un UA multi-niveau (ua4d)
   *           -1 s'il y a une erreur
   *
   ***************************************************************************/
int clipping_vertical(BURP_RPT *rptin, optionsptr optptr, gridtype* grid_gz, BURP_RPT *rptout) {
  int e,v,t,status,trouve_data,trouve_marqueur,EXIT_STATUS = 0;
  int rangee_alt, rangee_lat, rangee_lon;
  BURP_BLK *blk, *blkout, *blk_donnees, *blk_marqueur, *new_blk_donnees, *new_blk_marqueur;

  brp_copyrpthdr(rptout,rptin);

  blk = (BURP_BLK*) NULL;
  blkout = (BURP_BLK*) NULL;

  blk = brp_newblk();

  BLK_SetBKNO(blk, 0);
  while ( brp_findblk( blk, rptin ) >= 0 ) {
    int is_data = 0, is_marqueur = 0, btyp, btypalt, btypsfc;

    blkout = brp_newblk();
    status = brp_readblk(BLK_BKNO(blk), blkout, rptin,0);
    if (status<0) {
      fprintf(stderr,"Fonction clipping_vertical: Erreur %d dans la fonction brp_readblk\n", status);
      EXIT_STATUS = -1;
      brp_freeblk(blkout);
      break;
    }

    btyp = BLK_BTYP(blkout);
    btypalt = btyp>>4 & 31;
    btypsfc = btyp>>1 & 1;

    /* On verifie si c'est un derivate en altitude */
    if ( btypalt == 2 || btypalt == 3 ) {
      /* Si le bloc data est allume alors on le stocke dans le bloc blk_data  */
      /* On verifie si c'est un bloc de data */
      if ( (btyp>>11 & 3) == 0 ) {
	is_data = 1;
	trouve_data = 1;
	blk_donnees = brp_newblk();
	brp_copyblk(blk_donnees, blkout);
      }
      /* Si le bloc marqueur est allume alors on le stocke dans le bloc blk_marqueur  */
      /* On verifie si c'est un bloc marqueur */
      else if ( (btyp>>11 & 3) == 3 ) {
	is_marqueur = 1;
	trouve_marqueur = 1;
	blk_marqueur = brp_newblk();
	brp_copyblk(blk_marqueur, blkout);
      }
    } /* Fin du if ( btypalt == 2 || btypalt ==3 ) */

      /* Fin du if ( btypsfc == 0 ) */

      /* Si ce n'est pas un bloc de donnees ou un bloc marqueur alors, on le copie en entier */
    if (is_data==0 || is_marqueur==0) {
      status = putblk_nt(rptout,blkout,(int*) NULL,0);
      if (status!=0) {
	fprintf(stderr,"Fonction clipping_vertical: Erreur %d dans la fonction putblk_nt (btyp %d)\n", 
		status, btyp);
	brp_freeblk(blkout);
	EXIT_STATUS = 1;
	break;
      }	      
    }
    else if (trouve_data && trouve_marqueur) {
      /* On alloue de nouveaux blocs de la meme grandeur */
      if (BLK_NELE(blk_donnees) != BLK_NELE(blk_marqueur) ||
	  BLK_NVAL(blk_donnees) != BLK_NVAL(blk_marqueur)) {
	fprintf(stderr,"Fonction clipping_vertical: Les blocs data et marqueur "
		"ne sont pas de la meme taille!!!\n");

	brp_freeblk(blk_donnees);
	brp_freeblk(blk_marqueur);

	EXIT_STATUS = 1;
	break;
      }

      new_blk_donnees  = brp_newblk();
      new_blk_marqueur = brp_newblk();

      BLK_SetSTORE_TYPE(new_blk_donnees,BLK_STORE_TYPE(blk_donnees));
      BLK_SetSTORE_TYPE(new_blk_marqueur,BLK_STORE_TYPE(blk_marqueur));

      brp_allocblk(new_blk_donnees,BLK_NELE(blk_donnees),BLK_NVAL(blk_donnees),BLK_NT(blk_donnees));
      brp_allocblk(new_blk_marqueur,BLK_NELE(blk_donnees),BLK_NVAL(blk_donnees),BLK_NT(blk_donnees));

      /* On copie l'info de l'entete du bloc */
      /* bloc data */
      BLK_SetBKNO (new_blk_donnees, BLK_BKNO(blk_donnees) );
      BLK_SetBDESC(new_blk_donnees, BLK_BDESC(blk_donnees));
      BLK_SetBTYP (new_blk_donnees, BLK_BTYP(blk_donnees) );
      BLK_SetNBIT (new_blk_donnees, BLK_NBIT(blk_donnees) );
      BLK_SetDATYP(new_blk_donnees, BLK_DATYP(blk_donnees));
      BLK_SetBFAM (new_blk_donnees, BLK_BFAM(blk_donnees) );
      BLK_SetBKNAT(new_blk_donnees, BLK_BKNAT(blk_donnees));
      BLK_SetBKTYP(new_blk_donnees, BLK_BKTYP(blk_donnees));
      BLK_SetBKSTP(new_blk_donnees, BLK_BKSTP(blk_donnees));
      /* bloc marqueur */
      BLK_SetBKNO (new_blk_marqueur, BLK_BKNO(blk_marqueur) );
      BLK_SetBDESC(new_blk_marqueur, BLK_BDESC(blk_marqueur));
      BLK_SetBTYP (new_blk_marqueur, BLK_BTYP(blk_marqueur) );
      BLK_SetNBIT (new_blk_marqueur, BLK_NBIT(blk_marqueur) );
      BLK_SetDATYP(new_blk_marqueur, BLK_DATYP(blk_marqueur));
      BLK_SetBFAM (new_blk_marqueur, BLK_BFAM(blk_marqueur) );
      BLK_SetBKNAT(new_blk_marqueur, BLK_BKNAT(blk_marqueur));
      BLK_SetBKTYP(new_blk_marqueur, BLK_BKTYP(blk_marqueur));
      BLK_SetBKSTP(new_blk_marqueur, BLK_BKSTP(blk_marqueur));

      /* On trouve l'elements pour l'elevation en pression (7004) ou en metre (7006) */
      for (e=0;e<BLK_NELE(blk_donnees);e++) {
	if (BLK_DLSTELE(blk_donnees,e)==7004 || BLK_DLSTELE(blk_donnees,e)==7006)
	  rangee_alt = e;
	else if (BLK_DLSTELE(blk_donnees,e)==5001 || BLK_DLSTELE(blk_donnees,e)==5002)
	  rangee_lat = e;
	else if (BLK_DLSTELE(blk_donnees,e)==6001 || BLK_DLSTELE(blk_donnees,e)==6002)
	  rangee_lon = e;
	BLK_SetDLSTELE(new_blk_donnees,e,BLK_DLSTELE(blk_donnees,e));
	BLK_SetDLSTELE(new_blk_marqueur,e,BLK_DLSTELE(blk_marqueur,e));
      }

      status = brp_encodeblk(new_blk_donnees);
      if (status<0) {
	fprintf(stderr,"Fonction clipping_vertical: Erreur %d avec la fonction brp_encodeblk pour "
		"le bloc btyp=%d\n", status, BLK_BTYP(new_blk_donnees));
	EXIT_STATUS = 1;
	break;
      }

      status = brp_encodeblk(new_blk_marqueur);
      if (status<0) {
	fprintf(stderr,"Fonction clipping_vertical: Erreur %d avec la fonction brp_encodeblk pour "
		"le bloc btyp=%d\n", status, BLK_BTYP(new_blk_marqueur));
	EXIT_STATUS = 1;
	break;
      }

      /* On construit les nouveaux blocs data et marqueur */
      for (t=0;t<BLK_NT(blk_donnees);t++) {
	int compteur_rangee = 0;
	for (v=0;v<BLK_NVAL(blk_donnees);v++) {
	  int garde_rangee = 0;
	  int alt = BLK_RVAL(blk_donnees,rangee_alt,v,t);

	  if (strlen(optptr->channels)==0 && strlen(optptr->gz)==0 &&
	      checkvertical(alt,optptr->niveau_min,optptr->niveau_max))
	    /* Le filtrage vertical est fait a l'aide d'une hauteur en pression */
	    garde_rangee = 1;
	  else if (strlen(optptr->channels)==0) {
	    float lat = BLK_RVAL(blk_donnees,rangee_lat,v,t);
	    float lon = BLK_RVAL(blk_donnees,rangee_lon,v,t);

	    if(checkvertical_gz(lat,lon,alt,grid_gz->gridid,grid_gz->ni,grid_gz->nj,
				optptr->niveau_min,optptr->niveau_max))
	      /* Le filtrage vertical est fait a l'aide d'une hauteur en metre */
	      garde_rangee = 1;
	  }
	  else if (optptr->channels_voulus==checkcanal(alt,optptr->channels))
	    /* Si on fait une selection a partir d'une liste de canaux */
	    /* Dans ce cas, on verifie les canaux s'ils sont dans la liste ou non selon l'option */
	    garde_rangee = 1;

	  if (garde_rangee) {
	    for (e=0;e<BLK_NELE(blk_donnees);e++) {
	      BLK_SetTBLVAL(new_blk_donnees,e,compteur_rangee,t,BLK_TBLVAL(blk_donnees,e,v,t));
	      BLK_SetTBLVAL(new_blk_marqueur,e,compteur_rangee,t,BLK_TBLVAL(blk_marqueur,e,v,t));
	    }
	    compteur_rangee++;
	  }
	} /* Fin du for (v=0;v<BLK_NVAL(blk_donnees);v++) */
      } /* Fin du for(t=0;t,BLK_NT(blk_donnees);t++) */

      status = putblk_nt(rptout,new_blk_donnees,(int*) NULL,0);
      brp_freeblk(new_blk_donnees);
      if (status!=0) {
	fprintf(stderr,"Fonction clipping_vertical: Erreur %d dans la fonction putblk_nt (btyp %d)\n",
		status, BLK_BTYP(new_blk_donnees));
	brp_freeblk(blkout);
	EXIT_STATUS = 1;
	break;
      }

      status = putblk_nt(rptout,new_blk_marqueur,(int*) NULL,0);
      brp_freeblk(new_blk_marqueur);
      if (status!=0) {
	fprintf(stderr,"Fonction clipping_vertical: Erreur %d dans la fonction putblk_nt (btyp %d)\n",
		status, BLK_BTYP(new_blk_donnees));
	brp_freeblk(blkout);
	EXIT_STATUS = 1;
	break;
      }
    } /* Fin du else if (trouve_data && trouve_marqueur) */
    brp_freeblk(blkout);
  } /* Fin du while ( brp_findblk( blk, rptin ) >= 0 ) */

  return EXIT_STATUS;
} /* Fin de la fonction clipping_vertical */


  /***************************************************************************
   * fonction: check_ua4d
   *
   * En entree, cette fonction prend
   *    rptin: un rapport BURP complet d'entree
   *
   * Cette fonction retourne:
   *            0 si ce n'est pas un UA multi-niveau (ua4d)
   *            1 si c'est un UA multi-niveau (ua4d)
   *            -1 s'il y a une erreur
   *
   ***************************************************************************/
int check_ua4d(BURP_RPT *rptin) {
  int status, n_blk_data;
  int *bknos_data, *btyps_data;

  /* Si le codtyp est different d'un radiosondage, alors ce n'est certainement pas un 'ua4d' */
  /* C'est un 'UA' si le codtyp est 32,33,34,35,36,37,38,132,135,136,137,138,139,140,141,142,159,160,161,162 */
  if (RPT_IDTYP(rptin)<32 ||
      (RPT_IDTYP(rptin)>38 && RPT_IDTYP(rptin)<132) ||
      (RPT_IDTYP(rptin)>132 && RPT_IDTYP(rptin)<135) ||
      (RPT_IDTYP(rptin)>142 && RPT_IDTYP(rptin)<159) ||
      RPT_IDTYP(rptin)>162)
    return 0;

  bknos_data = (int*) NULL;
  btyps_data = (int*) NULL;
  status = find_blk_data_in_rpt(rptin,5001,6001,&bknos_data,&btyps_data,&n_blk_data);
  if (status<0) {
    if (n_blk_data>0) {
      if (VERBOSE>2)
	printf("Fonction check_ua4d: Erreur dans la fonction 'find_blk_data_in_rpt'.  "
	       "Cette derniere retourne %d mais on a trouve %d blocs de donnees\n", status,
	       n_blk_data);
      return -1;
    }
  }
  if(n_blk_data==0) {
    if (VERBOSE>2)
      printf("Fonction check_ua4d: Dans ce rapport, on ne trouve pas la latitude "
	     "et la longitude.  Ce n'est donc pas un UA multi-niveau (ua4d)\n");
    status=0;
  } /* Fin du if (status<0) */
  else {
    if (VERBOSE>2)
      printf("Fonction check_ua4d: Dans le bloc %d (btyp=%d), on a trouve que c'etait un ua4d\n", bknos_data[0],btyps_data[0]);
    status=1;
  }

  if (bknos_data != (int*) NULL) free(bknos_data);
  if (btyps_data != (int*) NULL) free(btyps_data);

  return status;
} /* Fin de la fonction check_ua4d */


  /***************************************************************************
   * fonction: find_blk_data_in_rpt
   *
   * En entree, cette fonction prend
   *    rptin: un rapport BURP complet d'entree
   *    elem_lat: le numero de l'element qui contient la latitude  (5001 ou 5002)
   *    elem_lon: le numero de l'element qui contient la longitude (6001 ou 6002)
   *
   * En sortie, on donne
   *    bknos_data_ptr: un pointeur à un vecteur d'entier contenant les numeros de
   *                  blocs de donnees qui ont les elements 'elem_lat' et 'elem_lon'
   *    btyps_data_ptr: un pointeur à un vecteur d'entier contenant les 'btyp' des
   *                   blocs de donnees associes aux bknos du vecteur precedent.
   *    nblks: le nombre de blocs trouvés (donc la dimension des vecteurs '*bknos_data_ptr' et 'btyps_data_ptr')
   *
   * Cette fonction retourne:
   *            0 si aucune erreur n'a ete detectee
   *            -1 s'il y a une erreur
   *
   * Dans cette implementation, on pourrait eviter de passer deux fois travers les blocs
   * en allouant un premier grand vecteur quitte a en reallouer un nouveau si jamais
   * on devait deborder.
   *
   ***************************************************************************/
int find_blk_data_in_rpt(BURP_RPT *rptin, int elem_lat, int elem_lon,
			 int** bknos_data_ptr, int** btyps_data_ptr, int* nblks) {
  int i, e, btyp, status, nblkstmp, trouve_lat, trouve_lon;
  int trouve_au_moins_un_bloc_avec_lat_lon;
  BURP_BLK *blktmp, *blk;

  if (VERBOSE>5) {
    printf("Fonction find_blk_data_in_rpt: On travaille avec le rapport: \n");
    printf("Fonction find_blk_data_in_rpt: ");
    printf ( "hhmm   =%8d " , RPT_TEMPS(rptin)) ;
    printf ( "flgs   =%6d  ", RPT_FLGS(rptin)) ;
    printf ( "codtyp =%6d  ", RPT_IDTYP(rptin)) ;
    printf ( "stnids =%9s\n", RPT_STNID(rptin)) ;
    printf("Fonction find_blk_data_in_rpt: ");
    printf ( "blat   =%8d " , RPT_LATI(rptin)) ;
    printf ( "blon   =%6d  ", RPT_LONG(rptin)) ;
    printf ( "dx     =%6d  ", RPT_DX(rptin)) ;
    printf ( "dy     =%6d  ", RPT_DY(rptin)) ;
    printf ( "stnhgt =%6d\n", RPT_ELEV(rptin)) ;
    printf("Fonction find_blk_data_in_rpt: ");
    printf ( "yymmdd =%8d " , RPT_DATE(rptin)) ;
    printf ( "oars   =%6d  ", RPT_OARS(rptin)) ;
    printf ( "runn   =%6d  ", RPT_RUNN(rptin)) ;
    printf ( "nblk   =%6d  ", RPT_NBLK(rptin)) ;
    printf ( "dlay   =%6d\n", RPT_DRND(rptin)) ;
  }

  blktmp = brp_newblk();
  blk    = brp_newblk();

  trouve_au_moins_un_bloc_avec_lat_lon=0;

  /* On passe une premiere fois tous les blocs pour trouver le nombre
     de blocs contenant les elements 'elem_lat' et 'elem_lon' */
  BLK_SetBKNO(blktmp, 0);
  nblkstmp = 0;
  while ( brp_findblk( blktmp, rptin ) >= 0 ) {
    status = brp_readblk(BLK_BKNO(blktmp), blk, rptin,0);
    if (status<0) {
      fprintf(stderr,"Fonction find_blk_data_in_rpt: Erreur %d dans la fonction brp_readblk\n", status);
      brp_freeblk(blk);
      brp_freeblk(blktmp);
      return -1;
    }

    if (VERBOSE>5) {
      printf("Fonction find_blk_data_in_rpt: ");
      printf ( "blkno  =%6d  ", BLK_BKNO(blk)    ) ;
      printf ( "nele   =%6d  ", BLK_NELE(blk)    ) ;
      printf ( "nval   =%6d  ", BLK_NVAL(blk)    ) ;
      printf ( "nt     =%6d  ", BLK_NT(blk)      ) ;
      printf ( "bit0   =%6d\n", BLK_BIT0(blk)    ) ;
      printf("Fonction find_blk_data_in_rpt: ");
      printf ( "bdesc  =%6d  ", BLK_BDESC(blk)   ) ;
      printf ( "btyp   =%6d  ", BLK_BTYP(blk)    ) ;
      printf ( "nbit   =%6d  ", BLK_NBIT(blk)    ) ;
      printf ( "datyp  =%6d  ", BLK_DATYP(blk)   ) ;
      printf ( "bfam   =%6d\n", BLK_BFAM(blk)    ) ;
      printf("Fonction find_blk_data_in_rpt: which_btyp(%d)=%d\n", BLK_BTYP(blk), which_btyp(BLK_BTYP(blk)));
    }

    trouve_lat=0;
    trouve_lon=0;
    /* On cherche l'element lat et lon dans ce bloc */
    for (e=0;e<BLK_NELE(blk);e++) {
      /* L'element 5001 indique la valeur de latitude de l'observation */
      if (BLK_DLSTELE(blk,e)==elem_lat)
	trouve_lat=1;
      /* L'element 6001 indique la valeur de longitude de l'observation */
      else if (BLK_DLSTELE(blk,e)==elem_lon)
	trouve_lon=1;

      if (trouve_lat>0 && trouve_lon>0) {
	/* Si on le trouve, on doit confirmer que c'est bien un bloc de donnees */
	if (which_btyp(BLK_BTYP(blk))==0) {
	  trouve_au_moins_un_bloc_avec_lat_lon=1;
	  nblkstmp++;
	  break;
	} /* Fin du if (which_btyp(BLK_BTYP(blk))==0) */
      } /* Fin du if (trouve_lat>0 && trouve_lon>0) */
    } /* Fin du for (e=0;e<BLK_NELE(blk);e++) */
  } /* Fin du while ( brp_findblk( blktmp, rptin ) >= 0 ) */

  if (trouve_au_moins_un_bloc_avec_lat_lon == 0) {
    if (VERBOSE>5)
      printf("Fonction find_blk_data_in_rpt: On ne trouve aucun bloc de donnees\n");
    brp_freeblk(blk);
    brp_freeblk(blktmp);
    *nblks=0;
    return 0;
  }

  if (VERBOSE>5)
    printf("Fonction find_blk_data_in_rpt: on trouve %d blocs de donnees\n", nblkstmp);

  *bknos_data_ptr = (int*) NULL;
  *bknos_data_ptr = (int*) malloc(nblkstmp*sizeof(int));
  if (*bknos_data_ptr == (int*) NULL) {
    fprintf(stderr, "Fonction find_data_flag_in_rpt: Incapable d'allouer le vecteur '*bknos_data_ptr' de dimension %d de (int)\n", nblkstmp);
    brp_freeblk(blk);
    brp_freeblk(blktmp);
    return -1;
  }
  *btyps_data_ptr = (int*) NULL;
  *btyps_data_ptr = (int*) malloc(nblkstmp*sizeof(int));
  if (*btyps_data_ptr == (int*) NULL) {
    fprintf(stderr, "Fonction find_data_flag_in_rpt: Incapable d'allouer le vecteur '*btyps_data_ptr' de dimension %d de (int)\n", nblkstmp);
    free(*bknos_data_ptr);
    brp_freeblk(blk);
    brp_freeblk(blktmp);
    return -1;
  }

  *nblks = nblkstmp;

  nblkstmp = 0;
  BLK_SetBKNO(blktmp, 0);
  while ( brp_findblk( blktmp, rptin ) >= 0 ) {
    status = brp_readblk(BLK_BKNO(blktmp), blk, rptin,0);
    if (status<0) {
      fprintf(stderr,"Fonction find_blk_data_in_rpt: Erreur %d dans la fonction brp_readblk\n", status);
      brp_freeblk(blk);
      brp_freeblk(blktmp);
      free(*bknos_data_ptr);
      free(*btyps_data_ptr);
      return -1;
    }

    if (VERBOSE>5) {
      printf("Fonction find_blk_data_in_rpt: ");
      printf ( "blkno  =%6d  ", BLK_BKNO(blk));
      printf ( "nele   =%6d  ", BLK_NELE(blk));
      printf ( "nval   =%6d  ", BLK_NVAL(blk));
      printf ( "nt     =%6d  ", BLK_NT(blk)  );
      printf ( "bit0   =%6d\n", BLK_BIT0(blk));
      printf("Fonction find_blk_data_in_rpt: ");
      printf ( "bdesc  =%6d  ", BLK_BDESC(blk));
      printf ( "btyp   =%6d  ", BLK_BTYP(blk) );
      printf ( "nbit   =%6d  ", BLK_NBIT(blk) );
      printf ( "datyp  =%6d  ", BLK_DATYP(blk));
      printf ( "bfam   =%6d\n", BLK_BFAM(blk) );
    }

    trouve_lat=0;
    trouve_lon=0;
    /* On cherche l'element lat et lon dans ce bloc */
    for (e=0;e<BLK_NELE(blk);e++) {
      /* L'element 5001 indique la valeur de latitude de l'observation */
      if (BLK_DLSTELE(blk,e)==elem_lat)
	trouve_lat=1;
      /* L'element 6001 indique la valeur de longitude de l'observation */
      else if (BLK_DLSTELE(blk,e)==elem_lon)
	trouve_lon=1;

      if (trouve_lat>0 && trouve_lon>0) {
	/* Si on le trouve, on doit confirmer que c'est bien un bloc de donnees */
	if (which_btyp(BLK_BTYP(blk))==0) {
	  (*bknos_data_ptr)[nblkstmp] = BLK_BKNO(blk);
	  (*btyps_data_ptr)[nblkstmp] = BLK_BTYP(blk);
	  nblkstmp++;
	  if (VERBOSE>5)
	    printf("Fonction find_blk_data_in_rpt: On trouve btyp=%d au blkno=%d\n", BLK_BTYP(blk), BLK_BKNO(blk));
	  break;
	} /* Fin du if (which_btyp(BLK_BTYP(blk))==0) */
      } /* Fin du if (trouve_lat>0 && trouve_lon>0) */
    } /* Fin du for (e=0;e<BLK_NELE(blk);e++) */
  } /* Fin du while ( brp_findblk( blktmp, rptin ) >= 0 ) */

  brp_freeblk(blktmp);
  brp_freeblk(blk);

  if (VERBOSE>5) {
    printf("Fonction find_blk_data_in_rpt: Voici les blocs trouves (bknos, btyp): [");
    i=0;
    printf("%d %d", (*bknos_data_ptr)[i],(*btyps_data_ptr)[i]);
    for (i=1;i<nblkstmp;i++)
      printf(",%d %d", (*bknos_data_ptr)[i],(*btyps_data_ptr)[i]);
    printf("]\n");
  }

  return 0;
} /* Fin de la fonction find_blk_data_in_rpt */


  /***************************************************************************
   * fonction: find_blk_data_flag_in_rpt
   *
   * En entree, cette fonction prend
   *    rptin: un rapport BURP complet d'entree
   *    elem_lat: le numero de l'element qui contient la latitude  (5001 ou 5002)
   *    elem_lon: le numero de l'element qui contient la longitude (6001 ou 6002)
   *    bkno_data: le bkno du bloc de donnees que l'on cherche
   *
   * En sortie, on donne
   *    blk_data_ptr: pointeur qui contient l'adresse de la memoire pour le bloc de donnees (data)
   *    blk_flag_ptr: pointeur qui contient l'adresse de la memoire pour le bloc marqueur (flag)
   *    colonne_lat: qui contient l'indice de l'element pour trouver la latitude
   *    colonne_lon: qui contient l'indice de l'element pour trouver la longitude
   *      si -1, alors on n'a pas trouve la latitude ou la longitude.
   *
   * Cette fonction retourne:
   *            0 si aucune erreur n'a ete detectee
   *            -1 s'il y a une erreur
   *
   ***************************************************************************/
int find_blk_data_flag_in_rpt(BURP_RPT *rptin, int elem_lat, int elem_lon, int bkno_data,
			      BURP_BLK** blk_data_ptr, BURP_BLK** blk_flags_ptr,
			      int* colonne_lat_ptr, int* colonne_lon_ptr) {
  int e, status, btyp, btyp_data, btyp_flags, trouve_flags=0, trouve_data=0;
  BURP_BLK *blktmp, *blk;

  if (VERBOSE>5)
    printf("Fonction find_blk_data_flag_in_rpt: elem_lat=%d elem_lon=%d bkno_data=%d\n", elem_lat, elem_lon, bkno_data);

  *blk_data_ptr = (BURP_BLK*) NULL;
  *blk_flags_ptr = (BURP_BLK*) NULL;
  *colonne_lat_ptr=-1;
  *colonne_lon_ptr=-1;

  blk    = brp_newblk();
  status = brp_readblk(bkno_data, blk, rptin, 0);
  if (status<0) {
    fprintf(stderr,"Fonction find_blk_data_flag_in_rpt: Erreur %d dans la fonction brp_readblk "
	    "pour le bloc %d\n", status, bkno_data);
    brp_freeblk(blk);
    return -1;
  }

  if (VERBOSE>5) {
    printf("Fonction find_blk_data_flag_in_rpt: ");
    printf ( "blkno  =%6d  ", BLK_BKNO(blk));
    printf ( "nele   =%6d  ", BLK_NELE(blk));
    printf ( "nval   =%6d  ", BLK_NVAL(blk));
    printf ( "nt     =%6d  ", BLK_NT(blk)  );
    printf ( "bit0   =%6d\n", BLK_BIT0(blk));
    printf("Fonction find_blk_data_flag_in_rpt: ");
    printf ( "bdesc  =%6d  ", BLK_BDESC(blk));
    printf ( "btyp   =%6d  ", BLK_BTYP(blk) );
    printf ( "nbit   =%6d  ", BLK_NBIT(blk) );
    printf ( "datyp  =%6d  ", BLK_DATYP(blk));
    printf ( "bfam   =%6d\n", BLK_BFAM(blk) );
  }

  btyp_data  = BLK_BTYP(blk);
  /* Le btyp du bloc marqueur est egal au btyp du bloc de donnees plus 6144
   * (les bits 12 et 13 allumés, le bit 1 étant le moins significatif)
   *      2^(12-1)+2^(13-1) = 6144
   * Information de Jose Garcia
   */
  btyp_flags = btyp_data+6144;

  if (VERBOSE>5)
    printf("Fonction find_blk_data_flag_in_rpt: Le bloc de donnes bkno=%d a le btyp=%d.  On cherche le bloc marqueur avec btyp=%d\n", BLK_BKNO(blk), btyp_data, btyp_flags);

  /* On cherche le bloc info qui permettra de connaitre les lat-lon de chaque obs */
  blktmp = brp_newblk();
  BLK_SetBKNO(blktmp, 0);
  while ( brp_findblk( blktmp, rptin ) >= 0 ) {
    status = brp_getblk(BLK_BKNO(blktmp), blk, rptin);
    if (status<0) {
      fprintf(stderr,"Fonction find_blk_data_flag_in_rpt: Erreur %d dans la fonction "
	      "brp_getblk pour bkno=%d\n", status, BLK_BKNO(blktmp));
      brp_freeblk(blk);
      brp_freeblk(blktmp);
      if (*blk_data_ptr != (BURP_BLK*) NULL) brp_freeblk(*blk_data_ptr);
      if (*blk_flags_ptr != (BURP_BLK*) NULL) brp_freeblk(*blk_flags_ptr);
      return -1;
    }

    btyp = BLK_BTYP(blk);

    if (VERBOSE>5)
      printf("Fonction find_blk_data_flag_in_rpt: btyp=%d (btyp>>11 & 3)=%d (btyp>>1 & 1)=%d\n", btyp, btyp>>11 & 3, btyp>>1 & 1);

    if (btyp == btyp_data) {
      trouve_data = 1;

      *blk_data_ptr = brp_newblk();
      brp_copyblk(*blk_data_ptr, blk);

      for (e=0;e<BLK_NELE(*blk_data_ptr);e++) {
	/* L'element 5001 indique la valeur de latitude de l'observation */
	if (BLK_DLSTELE(*blk_data_ptr,e)==elem_lat)
	  *colonne_lat_ptr=e;
	/* L'element 6001 indique la valeur de longitude de l'observation */
	else if (BLK_DLSTELE(*blk_data_ptr,e)==elem_lon)
	  *colonne_lon_ptr=e;
	
	if (*colonne_lat_ptr>=0 && *colonne_lon_ptr>=0)
	  break;
      }
    } /* Fin du 'if (btyp == btyp_data)' */

    if (btyp == btyp_flags) {
      trouve_flags = 1;

      *blk_flags_ptr = brp_newblk();
      brp_copyblk(*blk_flags_ptr, blk);

    } /* Fin du 'if (btyp == btyp_flags)' */

    if (trouve_flags && trouve_data)
      break;
  } /* Fin du while ( brp_findblk( blk, rptin ) >= 0 ) */

  brp_freeblk(blk);
  brp_freeblk(blktmp);

  if (trouve_data == 0) {
    fprintf(stderr,"Fonction find_blk_data_flag_in_rpt: le bloc data n'a pas ete trouve dans cet enregistrement\n");
    if (*blk_data_ptr != (BURP_BLK*) NULL) brp_freeblk(*blk_data_ptr);
    if (*blk_flags_ptr != (BURP_BLK*) NULL) brp_freeblk(*blk_flags_ptr);
    return -1;
  }

  if (trouve_flags == 0) {
    *blk_flags_ptr = (BURP_BLK*) NULL;
    if (VERBOSE>0)
      fprintf(stderr,"Fonction find_blk_data_flag_in_rpt: le bloc marqueur n'a pas ete trouve dans cet enregistrement\n");
  }

  if (*colonne_lat_ptr<0 || *colonne_lon_ptr<0) {
    if (VERBOSE>2)
      fprintf(stdout,"Fonction find_blk_data_flag_in_rpt: On n'a pas trouve "
	      "les lat-lon dans le bloc data de cet enregistrement\n");
    return -1;
  }

  if (VERBOSE>2)
    printf("Fonction find_blk_data_flag_in_rpt: Dans le bloc %d, on trouve la latitude "
	   "a l'element %d et la longitude a l'element %d\n",
	   BLK_BKNO(*blk_data_ptr), *colonne_lat_ptr, *colonne_lon_ptr);

  return 0;
} /* Fin de la fonction find_blk_data_flag_in_rpt */


  /***************************************************************************
   * fonction: fill_rptout_blk
   *
   * En entree, cette fonction prend
   *    rptin: un rapport BURP complet d'entree
   *    rptout: un pointeur a une serie de rapport BURP d'une longueur 'n'
   *    nts: nombre d'observations presentes dans chaque fichier pointe par 'rptout' (longueur 'n')
   *    t_in_domain: vecteur contenant les observations appartenant a chaque enregistrement pointe par 'rptout'
   *    n: longueur des vecteurs 'rptout' et 'nts'
   *    cherrypick_x et cherrypick_y: les tuiles que l'on veut extraire
   *    npey: le nombre de tuiles en y
   *
   * Cette fonction retourne:
   *            0 si la fonction ne rencontre aucune erreur
   *            -1 s'il y a une erreur
   *
   ***************************************************************************/
int fill_rptout_blk(BURP_RPT *rptin, BURP_RPT ** rptout, int* nts, int* t_in_domain, int n, int cherrypick_x,int cherrypick_y, int npey) {
  BURP_BLK *blk, *blkout;
  int i, EXIT_STATUS = 0, status;

  blk = brp_newblk();

  BLK_SetBKNO(blk, 0);
  while ( brp_findblk( blk, rptin ) >= 0 ) {
    blkout = brp_newblk();
    /* On utilise 'readblk' avec docvt = 0 pour ne pas convertir les valeurs puisque
     * cette operation change les valeurs dans certains cas tres particuliers.
     */
    status = brp_readblk(BLK_BKNO(blk), blkout, rptin, 0);
    if (status<0) {
      fprintf(stderr,"Fonction fill_rptout_blk: Erreur %d dans la fonction brp_readblk pour blk_no=%d\n", status, BLK_BKNO(blk));
      brp_freeblk(blkout);
      EXIT_STATUS = -1;
      break;
    }

    for(i=0;i<n;i++) {
      if (cherrypick_x > 0 && cherrypick_y > 0) {
        int cherrypick_id = (cherrypick_x-1)*npey+cherrypick_y-1;
        if (i != cherrypick_id)
          continue;
      }

      if (nts!= (int*) NULL) {
	if (nts[i]==0) {
	  if (VERBOSE>5)
	    printf("Fonction fill_rptout_blk: aucune observation pour id=%d\n", i);
	  continue;
	}
	else
	  if (VERBOSE>5)
	    printf("Fonction fill_rptout_blk: %d observations pour id=%d\n", nts[i], i);
      }

      if (t_in_domain != (int*) NULL)
	status = putblk_nt(rptout[i],blkout,&t_in_domain[BLK_NT(blkout)*i],nts[i]);
      else
	status = putblk_nt(rptout[i],blkout,(int*) NULL, 0);

      if (status!=0) {
	fprintf(stderr,"Fonction fill_rptout_blk: Erreur %d dans la fonction putblk_nt pour btyp=%d "
		"pour le id=%d\n", status, BLK_BTYP(blkout), i);
	EXIT_STATUS = -1;
	brp_freeblk(blkout);
	break;
      }
    } /* Fin du for(i=0;i<opt.npex*opt.npey;i++) */
    brp_freeblk(blkout);
  } /* Fin du 'while ( brp_findblk( blk, rptin ) >= 0 ) */

  brp_freeblk(blk);

  return EXIT_STATUS;
} /* Fin de la fonction fill_rptout_blk */


  /***************************************************************************
   * fonction: extract_data_in_domains_along_nt
   *
   * En entree, cette fonction prend
   *    optptr: un pointeur a une structure 'options' qui permet d'extraire
   *            l'information sur npex et npey et le rectanle du domaine
   *    gridptr: un pointeur a une structure 'gridtyp' qui permet d'utiliser EZSCINT
   *    rptin: un rapport BURP complet d'entree
   *    elem_lat: element donnant la latitude  (5001 ou 5002) (mais surement 5002 puisque ce sont probablement des observations satellitaires)
   *    elem_lon: element donnant la longitude (6001 ou 6002) (mais surement 6002 puisque ce sont probablement des observations satellitaires)
   *
   * En sortie, on donne les elements suivants:
   *    nts: liste de nombre representant le nombre d'observations dans chaque sous-domaine
   *    t_in_domain_ptr: tableau donnant a quel sous-domaine chaque observation appartient
   *
   * Cette fonction retourne:
   *            0 si la fonction ne rencontre aucune erreur
   *            1 s'il y a une erreur
   *
   ***************************************************************************/
int extract_data_in_domains_along_nt(optionsptr optptr, gridtype* gridptr, BURP_RPT *rptin,
				     int elem_lat, int elem_lon, int* nts, int** t_in_domain_ptr) {
  int status, t, colonne_lat, colonne_lon, EXIT_STATUS;
  int ilonband=1, jlatband=1;
  int btypnum, *bknos_data, *bknos_flags;
  float lat, lon;
  char errmsg[MAXSTR];
  BURP_BLK *blk_data, *blk_flags;

  EXIT_STATUS = 0;

  bknos_data = (int*) NULL;
  bknos_flags = (int*) NULL;
  status = find_blk_data_in_rpt(rptin,elem_lat,elem_lon,&bknos_data,&bknos_flags,&btypnum);
  if (status<0) {
    fprintf(stderr,"Fonction extract_data_in_domains_along_nt: Erreur dans la fonction find_blk_data_in_rpt\n");
    return 1;
  }

  if (btypnum!=1) {
    fprintf(stderr,"Fonction extract_data_in_domains_along_nt: Dans la fonction find_blk_data_in_rpt, "
	    "on a trouve %d blocs de donnees qui sont: ", btypnum);
    for (t=0;t<btypnum;t++)
      fprintf(stderr,"%d ", bknos_data[t]);
    fprintf(stderr,"\n");
    fprintf(stderr,"Or, ce programme ne peut traiter des enregistrements contenant plus d'un bloc de donnees regroupees\n");
    free(bknos_data);
    free(bknos_flags);
    return 1;
  }

  colonne_lat = -1;
  colonne_lon = -1;
  blk_data  = (BURP_BLK*) NULL;
  blk_flags = (BURP_BLK*) NULL;
  status = find_blk_data_flag_in_rpt(rptin,elem_lat,elem_lon,bknos_data[0],&blk_data,&blk_flags,&colonne_lat,&colonne_lon);
  if (status<0) {
    fprintf(stderr,"Fonction extract_data_in_domains_along_nt: Erreur dans la fonction find_blk_data_flag_in_rpt\n");
    free(bknos_data);
    free(bknos_flags);
    return 1;
  }

  *t_in_domain_ptr = (int*) malloc(BLK_NT(blk_data)*optptr->npex*optptr->npey*sizeof(int));
  if (*t_in_domain_ptr == (int*) NULL) {
    fprintf(stderr,"Fonction extract_data_in_domains_along_nt: Incapable d'allouer la memoire pour un tableau de "
	    "int de dimension %dx%dx%d=%d\n", BLK_NT(blk_data), optptr->npex, optptr->npey,
	    BLK_NT(blk_data)*optptr->npex*optptr->npey);
    free(bknos_data);
    free(bknos_flags);
    brp_freeblk(blk_data);
    brp_freeblk(blk_flags);
    return 1;
  }

  for (t=0;t<BLK_NT(blk_data);t++) {
    lat=BLK_RVAL(blk_data,colonne_lat,0,t);
    lon=BLK_RVAL(blk_data,colonne_lon,0,t);

    if ( optptr->npex == 1 && optptr->npey == 1 ) {
      /* On verifie si on est dans le domaine */
      status = checkgrid(gridptr->gridid, gridptr->ni, gridptr->nj, lat, lon, optptr->rect, errmsg);
      if (status<0) {
	fprintf(stderr,"Fonction extract_data_in_domains_along_nt: Erreur dans la fonction checkgrid pour le lat=%f "
		"et lon=%f avec le message '%s'\n", lat, lon, errmsg);
	EXIT_STATUS = 1;
	break;
      }

      /* Ceci signifie que si opt.inout == 1 alors status == 0 et donc
       * le point est hors de la grille ce qui n'est pas voulu
       *
       * ou bien qui si opt.inout == 0 alors status == 1 et donc le
       * point est a l'interieur de la grille ce qui n'est pas voulu.
       */
      if (optptr->inout == status) {
	(*t_in_domain_ptr)[nts[0]] = t;
	nts[0]++;
      }
    } /* Fin du if ( optptr->npex == 1 || optptr->npey == 1 ) */
    else {
      ilonband=-1;
      jlatband=-1;
      status = find_subdomain(gridptr->gridid, gridptr->ni, gridptr->nj, lat, lon, optptr->rect,
			      optptr->npex, optptr->npey, &ilonband, &jlatband, errmsg);
      if (status<0) {
	fprintf(stderr,"Fonction extract_data_in_domains_along_nt: Erreur dans la fonction find_subdomain "
		"pour le lat=%f et lon=%f avec le message '%s'\n", lat, lon, errmsg);
	EXIT_STATUS = 1;
	break;
      }
      if (ilonband<1 || ilonband>optptr->npex || jlatband<1 || jlatband>optptr->npey) {
	if (VERBOSE>1)
	  printf("Fonction main: cette observation n'est pas dans le domaine npex=%d, npey=%d "
		 "ilonband=%d jlatband=%d lat=%f lon=%f\n", optptr->npex, optptr->npey, ilonband,
		 jlatband, lat, lon);
      }
      else {
	int id = (ilonband-1)*optptr->npey+(jlatband-1);

	(*t_in_domain_ptr)[nts[id]+id*BLK_NT(blk_data)] = t;
	nts[id]++;

	if (VERBOSE>4)
	  printf("Fonction main: Observation acceptee: ilonband=%d jlatband=%d id=%d t_in_domain[%d]=%d nts[%d]=%d\n",
		 ilonband,jlatband,id,nts[id]+id*BLK_NT(blk_data)-1,t,id,nts[id]);
      }
    } /* Fin du else associe au if ( opt.npex == 1 || opt.npey == 1 ) */
  } /* Fin du for (t=0;t<BLK_NT(blkout);t++)  */

    /* 	printf("nt = %d\n", nt); */
  brp_freeblk(blk_data);
  brp_freeblk(blk_flags);
  free(bknos_data);
  free(bknos_flags);

  if (EXIT_STATUS!=0) {
    free(*t_in_domain_ptr);
    return 1;
  }

  return 0;
} /* Fin de la fonction extract_data_in_domains_along_nt */


  /***************************************************************************
   * fonction: extract_data_in_domains_along_nval
   *
   * En entree, cette fonction prend
   *    optptr: un pointeur a une structure 'options' qui permet d'extraire
   *            l'information sur npex et npey et le rectanle du domaine
   *    gridptr: un pointeur a une structure 'gridtyp' qui permet d'utiliser EZSCINT
   *    rptin: un rapport BURP complet d'entree
   *    elem_lat: element donnant la latitude  (5001 ou 5002)
   *    elem_lon: element donnant la longitude (6001 ou 6002)
   *
   * En sortie, on donne les elements suivants:
   *    nvals_in_domain: liste de nombre representant le nombre d'observations dans chaque sous-domaine
   *    val_in_domain:   tableau donnant a quel sous-domaine chaque observation appartient
   *
   * Cette fonction retourne:
   *            0 si la fonction ne rencontre aucune erreur
   *            1 s'il y a une erreur
   *
   ***************************************************************************/
int extract_data_in_domains_along_nval(optionsptr optptr, gridtype* gridptr, BURP_RPT *rptin,
				       int elem_lat, int elem_lon, BURP_BLK *blk_data,
				       int* nvals_in_domain, int* val_in_domain) {
  int status, i, id, e, v, btyp_data, colonne_lat, colonne_lon, EXIT_STATUS;
  int ilonband=1, jlatband=1;
  float lat, lon;
  char errmsg[MAXSTR];
  BURP_BLK *blk_data_converted = (BURP_BLK*) NULL;

  EXIT_STATUS = 0;

  btyp_data = BLK_BTYP(blk_data);

  blk_data_converted = brp_newblk();
  status = brp_readblk(BLK_BKNO(blk_data), blk_data_converted, rptin, 1);

  if (val_in_domain == (int*) NULL) {
    fprintf(stderr,"Fonction extract_data_in_domains_along_nval: le vecteur 'val_in_domain' doit etre deja alloue\n");
    brp_freeblk(blk_data_converted);
    return 1;
  }

  for (v=0;v<BLK_NVAL(blk_data)*optptr->npex*optptr->npey;v++)
    val_in_domain[v] = 0;

  colonne_lat=-1;
  colonne_lon=-1;
  for (e=0;e<BLK_NELE(blk_data);e++) {
    /* L'element 5001 indique la valeur de latitude de l'observation */
    if (BLK_DLSTELE(blk_data,e)==elem_lat)
      colonne_lat=e;
    /* L'element 6001 indique la valeur de longitude de l'observation */
    else if (BLK_DLSTELE(blk_data,e)==elem_lon)
      colonne_lon=e;
    if (colonne_lat>=0 && colonne_lon>=0)
      break;
  } /* Fin du for (e=0;e<BLK_NELE(blk);e++) */

  if (colonne_lat<0 || colonne_lon<0) {

    if (colonne_lat<0 && colonne_lon<0)
      fprintf(stderr,"Fonction extract_data_in_domains_along_nval: ne trouve pas "
	      "les elements %d et %d dans l'entete du bloc\n", elem_lat, elem_lon);
    else if (colonne_lat<0)
      fprintf(stderr,"Fonction extract_data_in_domains_along_nval: ne trouve pas "
	      "l'element %d dans l'entete du bloc\n", elem_lat);
    else if (colonne_lon<0)
      fprintf(stderr,"Fonction extract_data_in_domains_along_nval: ne trouve pas "
	      "l'element %d dans l'entete du bloc\n", elem_lon);

    brp_freeblk(blk_data_converted);
    return 1;
  }

  for (v=0;v<BLK_NVAL(blk_data);v++) {
    lat=BLK_RVAL(blk_data_converted,colonne_lat,v,0);
    lon=BLK_RVAL(blk_data_converted,colonne_lon,v,0);

    if (VERBOSE>3)
      printf("Fonction extract_data_in_domains_along_nval: appel de 'checkgrid' ou 'find_subdomain' avec lat=%f et lon=%f\n", lat, lon);

    id=-1;

    if ( optptr->npex == 1 && optptr->npey == 1 ) {
      /* On verifie si on est dans le domaine */
      status = checkgrid(gridptr->gridid, gridptr->ni, gridptr->nj, lat, lon, optptr->rect, errmsg);
      if (status<0) {
	fprintf(stderr,"Fonction extract_data_in_domains_along_nval: Erreur dans la fonction checkgrid pour le lat=%f "
		"et lon=%f avec le message '%s'\n", lat, lon, errmsg);
	EXIT_STATUS = 1;
	break;
      }
	
      /* Ceci signifie que si opt.inout == 1 alors status == 0 et donc
       * le point est hors de la grille ce qui n'est pas voulu
       *
       * ou bien qui si opt.inout == 0 alors status == 1 et donc le
       * point est a l'interieur de la grille ce qui n'est pas voulu.
       */
      if (optptr->inout == status)
	id = 0;
    } /* Fin du if ( optptr->npex == 1 || optptr->npey == 1 ) */
    else {
      ilonband=-1;
      jlatband=-1;
      status = find_subdomain(gridptr->gridid, gridptr->ni, gridptr->nj, lat, lon, optptr->rect,
			      optptr->npex, optptr->npey, &ilonband, &jlatband, errmsg);
      if (status<0) {
	fprintf(stderr,"Fonction extract_data_in_domains_along_nval: Erreur dans la fonction find_subdomain "
		"pour le lat=%f et lon=%f avec le message '%s'\n", lat, lon, errmsg);
	EXIT_STATUS = 1;
	break;
      }
      if (ilonband<1 || ilonband>optptr->npex || jlatband<1 || jlatband>optptr->npey) {
	if (VERBOSE>1)
	  printf("Fonction extract_data_in_domains_along_nval: cette observation n'est pas dans le domaine npex=%d, npey=%d "
		 "ilonband=%d jlatband=%d lat=%f lon=%f\n", optptr->npex, optptr->npey, ilonband,
		 jlatband, lat, lon);
      }
      else {
	id = (ilonband-1)*optptr->npey+(jlatband-1);
      }
    } /* Fin du else associe au if ( opt.npex == 1 || opt.npey == 1 ) */

    if (id>=0) {
      val_in_domain[id*BLK_NVAL(blk_data)+nvals_in_domain[id]] = v;
      nvals_in_domain[id]++;

      if (VERBOSE>4)
	printf("Fonction extract_data_in_domains_along_nval: Observation acceptee: "
	       "ilonband=%d jlatband=%d id=%d BLK_NVAL(blk_data)=%d "
	       "vals_in_domain[%d]=%d nvals[%d]=%d\n",
	       ilonband,jlatband,id,BLK_NVAL(blk_data),id*BLK_NVAL(blk_data)+nvals_in_domain[id],
	       v,id,nvals_in_domain[id]);
    }
  } /* Fin du for (v=0;v<BLK_NVAL(blk_data);v++)  */

  if (VERBOSE>1)
    for (id=0;id<optptr->npex*optptr->npey;id++)
      printf("Fonction extract_data_in_domains_along_nval: Il y a %d observations dans la bande %d\n", nvals_in_domain[id], id);

  brp_freeblk(blk_data_converted);

  if (EXIT_STATUS!=0) {
    fprintf(stderr,"Fonction extract_data_in_domains_along_nval: Une erreur dans la "
	    "selection des observations\n");
    return 1;
  }

  return 0;
} /* Fin de la fonction extract_data_in_domains_along_nval */


  /***************************************************************************
   * fonction: putblk_nt
   *
   * Cette fonction sert a inserer un bloc dans un enregistrement pour seulement
   * quelques tranches dans la 3e dimension donnee par le vecteur t_in_domain
   * pour les nt premiers elements
   *
   *      rpt: enregistrement dans lequel le bloc sera insere
   *      blk: bloc initial contenant toutes les observations
   *      t_in_domain: vecteur de (int) contenant les dimensions voulues
   *      nt: dimension du vecteur t_in_domain
   *
   ***************************************************************************/
int putblk_nt(BURP_RPT *rpt, BURP_BLK *blk, int* t_in_domain, int nt) {
  int e,v,t,tt,status = 0;
  BURP_BLK *newblk;

  newblk = brp_newblk();

  if (VERBOSE>2)
    printf("Fonction putblk_nt: btyp=%d nt=%d blk_nele=%d blk_nval=%d blk_nt=%d t_in_domain=%p\n", BLK_BTYP(blk), nt, BLK_NELE(blk),BLK_NVAL(blk),BLK_NT(blk),t_in_domain);
  if (VERBOSE>3 && nt!=0) {
    printf("Fonction putblk_nt: t_in_domain = [%d", t_in_domain[0]);
    for (t=1;t<nt;t++)
      printf(",%d", t_in_domain[t]);
    printf("]\n");
  }

  if (VERBOSE>4) {
    int thisi, thisj, thisk;
    printf("Fonction putblk_nt: Impression du 'blk' BLK_NT(blk)=%d BLK_NELE(blk)=%d BLK_NVAL(blk)=%d",BLK_NT(blk),BLK_NELE(blk),BLK_NVAL(blk));

    for ( thisk = 0 ; thisk < BLK_NT(blk) ; thisk++ ) {
      printf (  "\nlstele =" ) ;
      for (thisi=0;thisi<BLK_NELE(blk);thisi++)
	printf (  "    %6.6d", BLK_DLSTELE(blk,thisi) ) ;
      /* sortie des valeurs des elements */
      for (  thisj = 0 ; thisj < BLK_NVAL(blk) ; thisj++ ) {
	printf (  "\ntblval =" ) ;
	for (  thisi = 0 ; thisi < BLK_NELE(blk) ; thisi++ )
	  printf (  "%10d", BLK_TBLVAL(blk,thisi,thisj,thisk) ) ;
      }
    }
    printf("\nFonction putblk_nt: Impression du 'blk' terminee\n");
  }
  
  if (nt==0) 
    brp_allocblk(newblk,BLK_NELE(blk),BLK_NVAL(blk),BLK_NT(blk));
  else
    brp_allocblk(newblk,BLK_NELE(blk),BLK_NVAL(blk),nt);
  
  BLK_SetBKNO (newblk, BLK_BKNO(blk) );
  BLK_SetBDESC(newblk, BLK_BDESC(blk));
  BLK_SetBTYP (newblk, BLK_BTYP(blk) );
  BLK_SetNBIT (newblk, BLK_NBIT(blk) );
  BLK_SetDATYP(newblk, BLK_DATYP(blk));
  BLK_SetBFAM (newblk, BLK_BFAM(blk) );
  BLK_SetBKNAT(newblk, BLK_BKNAT(blk));
  BLK_SetBKTYP(newblk, BLK_BKTYP(blk));
  BLK_SetBKSTP(newblk, BLK_BKSTP(blk));

  for (e=0;e<BLK_NELE(blk);e++)
    BLK_SetDLSTELE(newblk,e,BLK_DLSTELE(blk,e));

  if (VERBOSE>4)
    printf("Fonction putblk_nt: newblk->datyp = %d (avant encode)\n", BLK_DATYP(newblk));

  status = brp_encodeblk(newblk);
  if (status<0) {
    fprintf(stderr,"Fonction putblk_nt: Erreur %d avec la fonction brp_encodeblk pour le bloc %d\n", status, BLK_BKNO(blk));
    brp_freeblk(newblk);
    return -1;
  }
  if (VERBOSE>4)
    printf("Fonction putblk_nt: newblk->datyp = %d (apres encode)\n", BLK_DATYP(newblk));

  for (tt=0;tt<BLK_NT(newblk);tt++) {
    if (nt==0)
      t=tt;
    else {
      if ( t_in_domain == (int*) NULL )
	t=tt;
      else
	t=t_in_domain[tt];
    }

    for (v=0;v<BLK_NVAL(blk);v++)
      for (e=0;e<BLK_NELE(blk);e++)
	BLK_SetTBLVAL(newblk,e,v,tt,BLK_TBLVAL(blk,e,v,t));
  }

  if (VERBOSE>3)
    fprintf(stdout,"Fonction putblk_nt: Impression du 'newblk' BLK_NT(newblk)=%d BLK_NELE(newblk)=%d BLK_NVAL(newblk)=%d",
	    BLK_NT(newblk),BLK_NELE(newblk),BLK_NVAL(newblk));
  if (VERBOSE>4) {
    int thisi, thisj, thisk;

    for ( thisk = 0 ; thisk < BLK_NT(newblk) ; thisk++ ) {
      printf (  "\nlstele =" ) ;
      for (thisi=0;thisi<BLK_NELE(newblk);thisi++)
	printf (  "    %6.6d", BLK_DLSTELE(newblk,thisi) ) ;
      /* sortie des valeurs des elements */
      for (  thisj = 0 ; thisj < BLK_NVAL(newblk) ; thisj++ ) {
	printf (  "\ntblval =" ) ;
	for (  thisi = 0 ; thisi < BLK_NELE(newblk) ; thisi++ )
	  printf (  "%10d", BLK_TBLVAL(newblk,thisi,thisj,thisk) ) ;
      }
    }
    printf("\nFonction putblk_nt: Impression du 'newblk' terminee");
  }
  else if (VERBOSE>3)
    fprintf(stdout,"\nFonction putblk_nt: Impression du 'newblk' terminee\n");

  status = brp_putblk(rpt,newblk);
  if (VERBOSE>3)
    fprintf(stdout,"Fonction putblk_nt: 'brp_putblk' terminee\n");

  if (status<0) {
    fprintf(stderr,"Fonction putblk_nt: Erreur %d dans la fonction brp_putblk (btyp %d, datyp=%d)\n",
	    status, BLK_BTYP(newblk), BLK_DATYP(newblk));
    brp_freeblk(newblk);
    return -1;
  }

  brp_freeblk(newblk);
  if (VERBOSE>2)
    printf("Fonction putblk_nt: btyp=%d nt=%d blk_nele=%d blk_nval=%d blk_nt=%d t_in_domain=%p return=0\n",
	   BLK_BTYP(blk), nt, BLK_NELE(blk),BLK_NVAL(blk),BLK_NT(blk),t_in_domain);
  
  return 0;
} /* Fin de la fonction putblk_nt */


/***************************************************************************
 * fonction: putblk_nval
 *
 * Cette fonction sert a inserer un bloc dans un enregistrement pour seulement
 * quelques tranches dans la dimension 'v' donnee par le vecteur vals_in_domain
 * pour les nval premiers elements
 *
 *      rpt: enregistrement dans lequel le bloc sera insere
 *      blk: bloc initial contenant toutes les observations
 *      vals_in_domain: vecteur de (int) contenant les dimensions voulues
 *      nval: dimension du vecteur vals_in_domain
 *
 ***************************************************************************/
int putblk_nval(BURP_RPT *rpt, BURP_BLK *blk, int* vals_in_domain, int nval) {
  int e,v,vv,t,status = 0;
  BURP_BLK *newblk;

  newblk = brp_newblk();

  if (VERBOSE>2)
    printf("Fonction putblk_nval: btyp=%d nval=%d blk_nele=%d blk_nval=%d blk_nt=%d vals_in_domain=%p\n", BLK_BTYP(blk), nval, BLK_NELE(blk),BLK_NVAL(blk),BLK_NT(blk),vals_in_domain);
  if (VERBOSE>3 && nval!=0) {
    printf("Fonction putblk_nval: vals_in_domain = [%d", vals_in_domain[0]);
    for (v=1;v<nval;v++)
      printf(",%d", vals_in_domain[v]);
    printf("]\n");
  }

  if (VERBOSE>4) {
    int thisi, thisj, thisk;
    printf("Fonction putblk_nval: Impression du 'blk' BLK_NT(blk)=%d BLK_NELE(blk)=%d BLK_NVAL(blk)=%d",BLK_NT(blk),BLK_NELE(blk),BLK_NVAL(blk));

    for ( thisk = 0 ; thisk < BLK_NT(blk) ; thisk++ ) {
      printf (  "\nlstele =" ) ;
      for (thisi=0;thisi<BLK_NELE(blk);thisi++)
	printf (  "    %6.6d", BLK_DLSTELE(blk,thisi) ) ;
      /* sortie des valeurs des elements */
      for (  thisj = 0 ; thisj < BLK_NVAL(blk) ; thisj++ ) {
	printf (  "\ntblval =" ) ;
	for (  thisi = 0 ; thisi < BLK_NELE(blk) ; thisi++ )
	  printf (  "%10d", BLK_TBLVAL(blk,thisi,thisj,thisk) ) ;
      }
    }
    printf("\nFonction putblk_nval: Impression du 'blk' terminee\n");
  }
  
  if (nval==0) 
    brp_allocblk(newblk,BLK_NELE(blk),BLK_NVAL(blk),BLK_NT(blk));
  else
    brp_allocblk(newblk,BLK_NELE(blk),nval,BLK_NT(blk));
  
  BLK_SetBKNO (newblk, BLK_BKNO(blk) );
  BLK_SetBDESC(newblk, BLK_BDESC(blk));
  BLK_SetBTYP (newblk, BLK_BTYP(blk) );
  BLK_SetNBIT (newblk, BLK_NBIT(blk) );
  BLK_SetDATYP(newblk, BLK_DATYP(blk));
  BLK_SetBFAM (newblk, BLK_BFAM(blk) );
  BLK_SetBKNAT(newblk, BLK_BKNAT(blk));
  BLK_SetBKTYP(newblk, BLK_BKTYP(blk));
  BLK_SetBKSTP(newblk, BLK_BKSTP(blk));

  for (e=0;e<BLK_NELE(blk);e++)
    BLK_SetDLSTELE(newblk,e,BLK_DLSTELE(blk,e));

  if (VERBOSE>4)
    printf("Fonction putblk_nval: newblk->datyp = %d (avant encode)\n", BLK_DATYP(newblk));

  status = brp_encodeblk(newblk);
  if (status<0) {
    fprintf(stderr,"Fonction putblk_nval: Erreur %d avec la fonction brp_encodeblk pour le bloc %d\n", status, BLK_BKNO(blk));
    brp_freeblk(newblk);
    return -1;
  }
  if (VERBOSE>4)
    printf("Fonction putblk_nval: newblk->datyp = %d (apres encode)\n", BLK_DATYP(newblk));

  for (vv=0;vv<BLK_NVAL(newblk);vv++) {
    if (nval==0)
      v=vv;
    else
      v=vals_in_domain[vv];

    for (t=0;t<BLK_NT(blk);t++)
      for (e=0;e<BLK_NELE(blk);e++)
	BLK_SetTBLVAL(newblk,e,vv,t,BLK_TBLVAL(blk,e,v,t));
  }

  if (VERBOSE>3)
    fprintf(stdout,"Fonction putblk_nval: Impression du 'newblk' BLK_NT(newblk)=%d BLK_NELE(newblk)=%d BLK_NVAL(newblk)=%d",BLK_NT(newblk),BLK_NELE(newblk),BLK_NVAL(newblk));
  if (VERBOSE>4) {
    int thisi, thisj, thisk;

    for ( thisk = 0 ; thisk < BLK_NT(newblk) ; thisk++ ) {
      printf (  "\nlstele =" ) ;
      for (thisi=0;thisi<BLK_NELE(newblk);thisi++)
	printf (  "    %6.6d", BLK_DLSTELE(newblk,thisi) ) ;
      /* sortie des valeurs des elements */
      for (  thisj = 0 ; thisj < BLK_NVAL(newblk) ; thisj++ ) {
	printf (  "\ntblval =" ) ;
	for (  thisi = 0 ; thisi < BLK_NELE(newblk) ; thisi++ )
	  printf (  "%10d", BLK_TBLVAL(newblk,thisi,thisj,thisk) ) ;
      }
    }
    printf("\nFonction putblk_nval: Impression du 'newblk' terminee\n");
  }
  else if (VERBOSE>3)
    fprintf(stdout,"Fonction putblk_nval: Impression du 'newblk' terminee\n");

  status = brp_putblk(rpt,newblk);
  if (VERBOSE>3)
    fprintf(stdout,"Fonction putblk_nval: 'brp_putblk' terminee\n");
  if (status<0) {
    fprintf(stderr,"Fonction putblk_nt: Erreur %d dans la fonction blk_putblk (btyp %d)\n", status, BLK_BTYP(blk));
    return -1;
  }
  brp_freeblk(newblk);

  if (VERBOSE>2)
    printf("Fonction putblk_nval: btyp=%d nval=%d blk_nele=%d blk_nval=%d blk_nt=%d vals_in_domain=%p return=0\n",
	   BLK_BTYP(blk), nval, BLK_NELE(blk),BLK_NVAL(blk),BLK_NT(blk),vals_in_domain);

  return 0;
} /* Fin de la fonction putblk_nval */


  /***************************************************************************
   * fonction: parseOptions
   *
   * Cette fonction sert a interpreter les arguments donnes lors de l'appel
   * du programme.
   *
   * Elle prend 3 arguments:
   *    argc: le nombre d'argument a l'appel du programme
   *    argv: pointeur contenant les arguments
   *    optptr: pointeur a une structure existante du type (options)
   ***************************************************************************/
int parseOptions(int argc, char** argv, optionsptr optptr) {
  int i;
  
  if (argc==1) { /* Alors, on veut la documentation */
    aide();
    exit(0);
  }

  strcpy(optptr->fst.nomvar,NOMVAR_DEFAUT);
  
  for (i=1;i<argc;i++) {
    if (argv[i][0]=='-') {
      if (!strcmp(argv[i],HELP_OPTION1) ||
	  !strcmp(argv[i],HELP_OPTION2) ||
	  !strcmp(argv[i],HELP_OPTION3)) { /* Alors, on veut la documentation */
	aide();
	exit(0);
      }
      else if (!strcmp(argv[i],VERBOSE_OPTION)) { /* On indique le niveau de print */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr, "Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  return NOT_OK;
	}
	VERBOSE = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],FSTIN_OPTION)) { /* Alors, on donne le nom du fichier standard input */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->fstin,argv[++i]);
      }
      else if (!strcmp(argv[i],OBSIN_OPTION)) { /* Alors, on donne le nom du fichier d'observations input */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->obsin,argv[++i]);
      }
      else if (!strcmp(argv[i],OBSOUT_OPTION)) { /* Alors, on donne le nom du fichier d'observations output */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->obsout,argv[++i]);
      }
      else if (!strcmp(argv[i],RDBIN_OPTION)) { /* Alors, on donne le nom du fichier RDB input */
        fprintf(stderr,"Fonction parseOptions: L'option %s a ete remplacee par %s\n", argv[i], OBSIN_OPTION);
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->obsin,argv[++i]);
      }
      else if (!strcmp(argv[i],RDBOUT_OPTION)) { /* Alors, on donne le nom du fichier RDB output */
        fprintf(stderr,"Fonction parseOptions: L'option %s a ete remplacee par %s\n", argv[i], OBSOUT_OPTION);
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->obsout,argv[++i]);
      }
      else if (!strcmp(argv[i],BURPIN_OPTION)) { /* Alors, on donne le nom du fichier BURP input */
        fprintf(stderr,"Fonction parseOptions: L'option %s a ete remplacee par %s\n", argv[i], OBSIN_OPTION);
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->obsin,argv[++i]);
      }
      else if (!strcmp(argv[i],BURPOUT_OPTION)) { /* Alors, on donne le nom du fichier BURP output */
        fprintf(stderr,"Fonction parseOptions: L'option %s a ete remplacee par %s\n", argv[i], OBSOUT_OPTION);
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->obsout,argv[++i]);
      }
      else if (!strcmp(argv[i],ASCII_OPTION)) { /* Alors, on donne le nom du fichier ASCII qui contient des couplets lat-lon de points */
        fprintf(stderr,"Fonction parseOptions: L'option %s a ete remplacee par %s\n", argv[i], OBSIN_OPTION);
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->obsin,argv[++i]);
      }
      else if (!strcmp(argv[i],OUT_OPTION)) { /* Alors, on donne le nom du fichier ASCII de sortie */
        fprintf(stderr,"Fonction parseOptions: L'option %s a ete remplacee par %s\n", argv[i], OBSOUT_OPTION);
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->obsout,argv[++i]);
      }
      else if (!strcmp(argv[i],INOUT_OPTION)) { /* On decide si on prend les observations a l'interieur (valeur = 1)
						 * ou a l'exterieur (valeur = 0) du domaine.   
						 */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->inout = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],PILOT_OPTION)) { /* On decide si on utilise une zone de pilotage pour rapetisser le domaine
						 */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->pilot = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],MIN_I_OPTION)) { /* On donne le 'i' minimal pour definir le domaine
						 */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	if (argv[i+1][0]=='=') {
	  optptr->rect.min_i = atof(&(argv[++i][1]));
	  optptr->rect.min_i_equal = 1;
	}
	else
	  optptr->rect.min_i = atof(argv[++i]);
      }
      else if (!strcmp(argv[i],MAX_I_OPTION)) { /* On donne le 'i' maximal pour definir le domaine
						 */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	if (argv[i+1][0]=='=') {
	  optptr->rect.max_i = atof(&(argv[++i][1]));
	  optptr->rect.max_i_equal = 1;
	}
	else
	  optptr->rect.max_i = atof(argv[++i]);
      }
      else if (!strcmp(argv[i],MIN_J_OPTION)) { /* On donne le 'j' minimal pour definir le domaine
						 */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	if (argv[i+1][0]=='=') {
	  optptr->rect.min_j = atof(&(argv[++i][1]));
	  optptr->rect.min_j_equal = 1;
	}
	else
	  optptr->rect.min_j = atof(argv[++i]);
      }
      else if (!strcmp(argv[i],MAX_J_OPTION)) { /* On donne le 'j' maximal pour definir le domaine
						 */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	if (argv[i+1][0]=='=') {
	  optptr->rect.max_j = atof(&(argv[++i][1]));
	  optptr->rect.max_j_equal = 1;
	}
	else
	  optptr->rect.max_j = atof(argv[++i]);
      }
      else if (!strcmp(argv[i],NPEX_OPTION)) { /* On donne le nombre de bandes dont on veut separer en 'i'
						*/
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->npex = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],NPEY_OPTION)) { /* On donne le nombre de bandes dont on veut separer en 'j'
						*/
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->npey = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],CHERRYPICK_X_OPTION)) { /* On donne la bande voulue en 'i' */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->cherrypick_x = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],CHERRYPICK_Y_OPTION)) { /* On donne la bande voulue en 'j' */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->cherrypick_y = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],NDIGITS_OPTION)) { /* On donne le nombre de catacteres que l'on veut
						   * avoir les extensions des fichiers separes en sous-domaines.
						   */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->ndigits = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],CHECK_UA4D_OPTION)) { /* Par defaut, on assume que l'on ne melange pas les enregistrements
						      * ua4d et ua classiques alors on ne fait la verification a chaque
						      * enregistrement.  Si cette option est activee, alors on fera cette
						      * verification a chaque enregistrement.
						      */
	if (i+1<argc && argv[i+1][0]!='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s ne demande aucun argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->check_ua4d = 1;
      }
      else if (!strcmp(argv[i],ROUNDROBIN_OPTION)) { /* Cette option activera le mode round-robin pour faire un splitting
						      * en fichiers sans egard a la position geographique.
						      */
	if (i+1<argc && argv[i+1][0]!='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s ne demande aucun argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->roundrobin = 1;
      }
      else if (!strcmp(argv[i],RDB_HEADER_OPTION)) { /* Cette option enregistrera le nom de la table qui joue le role du 'header' */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->rdb_header_table,argv[++i]);
      }
      else if (!strcmp(argv[i],RDB_DATA_OPTION)) { /* Cette option enregistrera le nom de la table qui joue le role du 'data' */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->rdb_data_table,argv[++i]);
      }
      else if (!strcmp(argv[i],RDB_PRIMARYKEY_OPTION)) { /* Cette option enregistrera le nom de la cle primaire qui lie 'header' et 'data' */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->rdb_primarykey,argv[++i]);
      }
      else if (!strcmp(argv[i],GZ_OPTION)) { /* Fichier standard dans lequel on lira le GZ au niveau voulu */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr, "Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  return NOT_OK;
	}
	strcpy(optptr->gz,argv[++i]);
      }
      else if (!strcmp(argv[i],NIVEAU_MIN_OPTION)) { /* On donne le niveau en hPa voulu pour accepter 
						      * les observations au dessus
						      */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr, "Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  return NOT_OK;
	}
	optptr->niveau_min = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],NIVEAU_MAX_OPTION)) { /* On donne le niveau en hPa voulu pour accepter 
						      * les observations en dessous
						      */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr, "Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  return NOT_OK;
	}
	optptr->niveau_max = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],CHANNELS_OPTION) || !strcmp(argv[i],NOCHANNELS_OPTION))  { 
	/* On specifie les canaux voulus ou que l'on ne veut pas */
	int indice_option = i;

	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr, "Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  return NOT_OK;
	}
	/* Si on donne les canaux avec l'option -nochannels, alors on eviter ces canaux plutot que de les prendre */
	if (!strcmp(argv[i],NOCHANNELS_OPTION)) optptr->channels_voulus = 0;
	i++;
	/* On cherche la prochaine option et on stocke tous les canaux donne dans optptr->channels */
	strcpy(optptr->channels,argv[i]);
	i++;
	while(i<argc && argv[i][0]!='-') {
	  if (strlen(optptr->channels) + strlen(argv[i])>=MAXSTR_CHANNELS) {
	    fprintf(stderr,"Fonction parseOptions: L'option %s ne peut prendre qu'un argument "
		    "d'un maximum de %d caracteres (incluant les espaces)\n", argv[indice_option], MAXSTR_CHANNELS);
	    return NOT_OK;
	  }
	  strcat(optptr->channels,argv[i]);
	  i++;
	}
	i--;
      }
      else if (!strcmp(argv[i],NOMVAR_OPTION)) { /* On a besoin du nomvar du champ */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->fst.nomvar,argv[++i]);
      }
      else if (!strcmp(argv[i],TYPVAR_OPTION)) { /* On a besoin du typvar du champ */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->fst.typvar,argv[++i]);
      }
      else if (!strcmp(argv[i],ETIKET_OPTION)) { /* On a besoin de l'etiquette du champ */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	strcpy(optptr->fst.etiket,argv[++i]);
      }
      else if (!strcmp(argv[i],DATEV_OPTION)) { /* On a besoin de la date du champ */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->fst.dateo = padtime(argv[++i]);
      }
      else if (!strcmp(argv[i],IP1_OPTION)) { /* On a besoin du ip1 du champ */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->fst.ip1 = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],IP2_OPTION)) { /* On a besoin du ip2 du champ */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->fst.ip2 = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],IP3_OPTION)) { /* On a besoin du ip3 du champ */
	if (i+1>=argc || argv[i+1][0]=='-') {
	  fprintf(stderr,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
	  exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
	}
	optptr->fst.ip3 = atoi(argv[++i]);
      }
      else {
	fprintf(stderr,"Fonction parseOptions: Je ne reconnais pas cette option: %s\n", argv[i]);
	exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
      }
    } /* Fin du if (argv[i][0] ...) */
    else {
      fprintf(stderr,"Fonction parseOptions: Erreur dans les arguments d'entree: (%d)\n\n", i);
      exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
    }
  } /* Fin du for (i=1; ...) */

  /* On checke les options recueillies */
  if ( strlen(optptr->fstin)==0 && optptr->roundrobin == 0 ) {
    fprintf(stderr,"Fonction parseOptions: On doit absolument specifier un fichier standard "
	    "en entree avec l'option %s\n", FSTIN_OPTION);
    exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
  }

  printf("On imprime les differentes options choisies pour cette execution\n");

  if ( strlen(optptr->obsin) !=0 && strlen(optptr->obsout) !=0 ) {
    printf("Fichier d'observations en entree:  %s\n", optptr->obsin);
    printf("Fichier d'observations en sortie:  %s\n", optptr->obsout);
  }
  else {
    if ( strlen(optptr->obsin) == 0 ) {
      fprintf(stderr,"Fonction parseOptions: On doit absolument specifier soit un fichier de base de "
              "donnees SQL, un fichier BURP ou ASCII en entree avec l'option %s\n", OBSIN_OPTION);
    }
    if ( strlen(optptr->obsout) == 0 ) {
      fprintf(stderr,"Fonction parseOptions: On doit absolument specifier soit un fichier de base de "
              "donnees SQL, un fichier BURP ou ASCII en sortie avec l'option %s\n", OBSOUT_OPTION);
    }
    exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
  }

  printf("\n");

  if ( optptr->roundrobin == 0 ) {
    if (optptr->inout)
      printf("On filtre pour ne garder que les observations dans le domaine (%s = %d)\n\n", INOUT_OPTION, optptr->inout);
    else
      printf("On filtre pour ne garder que les observations hors du domaine (%s = %d)\n\n", INOUT_OPTION, optptr->inout);
  }

  if (optptr->npex!=1 || optptr->npey!=1) {
    char npex_str[MAXSTR], npey_str[MAXSTR], format_digits[MAXSTR];

    printf("Le domaine sera separee en %d par %d parties egales grace aux options (%s et %s)\n", optptr->npex, optptr->npey, NPEX_OPTION, NPEY_OPTION);
    sprintf(format_digits,"%%.%dd",optptr->ndigits);
    sprintf(npex_str,format_digits,optptr->npex);
    sprintf(npey_str,format_digits,optptr->npey);
    printf("Les extensions auront %d caracteres par exemple %s_%s_%s\n\n", optptr->ndigits, optptr->obsout, npex_str, npey_str);
  }

  if (optptr->cherrypick_x>0 && optptr->cherrypick_y>0) {
    char x_str[MAXSTR], y_str[MAXSTR], format_digits[MAXSTR];
    sprintf(format_digits,"%%.%dd",optptr->ndigits);
    sprintf(x_str,format_digits,optptr->cherrypick_x);
    sprintf(y_str,format_digits,optptr->cherrypick_y);
    printf("On ne va extraire que les observations pour la tuile '%s_%s'\n\n",x_str,y_str);
  }
  else if (!(optptr->cherrypick_x<0 && optptr->cherrypick_y<0)) {
    fprintf(stderr,"Fonction parseOptions: On doit absolument specifier les deux arguments %s et %s en meme temps.\n", CHERRYPICK_X_OPTION, CHERRYPICK_Y_OPTION);
    exit_program(NOT_OK,PROGRAM_NAME,PROBLEM,VERSION);
  }

  if ( optptr->roundrobin == 1 ) {
    printf("On separera le fichier en mode round-robin.\n");
  }
  else {
    if (optptr->check_ua4d)
      printf("Une verification de chaque enregistrement sera faite pour savoir si c'est un fichier UA multi-niveau (UA4D)\n");

    if (optptr->niveau_min != IP1_VIDE || optptr->niveau_max != IP1_VIDE) {
      printf("Un filtrage vertical sera fait ");
      if (strlen(optptr->gz)!=0)
	printf("en lisant le GZ dans le fichier %s\n", optptr->gz);

      printf("On garde les observations ");

      if (optptr->niveau_min != IP1_VIDE)
	printf("\t\t\tplus hautes que %d hPa\n", optptr->niveau_min);

      if (optptr->niveau_max != IP1_VIDE)
	printf("\t\t\tplus basses que %d hPa\n", optptr->niveau_max);
    }
    else if (strlen(optptr->channels)!=0) {
      if (optptr->niveau_min != IP1_VIDE || optptr->niveau_max != IP1_VIDE) {
	fprintf(stderr, "Fonction parseOptions: Si on donne des canaux avec l'option %s, on ne peut utiliser les options %s et %s\n", CHANNELS_OPTION, NIVEAU_MIN_OPTION, NIVEAU_MAX_OPTION);

	return NOT_OK;
      }
      if (optptr->channels_voulus==1)
	printf("On ne gardera que les canaux %s\n", optptr->channels);
      else
	printf("On ne gardera pas les canaux %s ceux-ci seront exclus\n", optptr->channels);
    }
    else {
      printf("Aucun filtrage vertical ne sera effectue");
      if (strlen(optptr->gz)!=0)
	printf(" meme si on a donne un fichier standard %s avec l'option '%s'", optptr->gz, GZ_OPTION);
      printf(".  \n");
    }
    printf("\n");

    if (optptr->pilot != PILOT_DEFAUT &&
	optptr->rect.min_i<0 && optptr->rect.max_i<0 &&
	optptr->rect.min_j<0 && optptr->rect.max_j<0 ) {
      printf("Zone de pilotage: %d points (%s = %d)\n", optptr->pilot, PILOT_OPTION, optptr->pilot);
    }
    else {
      optptr->pilot=PILOT_DEFAUT;
      printf("On choisira la portion suivante de la grille (ce qui est donne avec l'option %s est ignore):\n", PILOT_OPTION);
      printf("   min_i = %g  et  max_i = %g\n", optptr->rect.min_i, optptr->rect.max_i);
      printf("   min_j = %g  et  max_j = %g\n", optptr->rect.min_j, optptr->rect.max_j);
      printf("Si -1 alors on prendra la grandeur de la grille\n");
    }
    printf("\n");
  
    printf("Fichier standard en entree qui contient la grille qui definit le domaine:  %s\n", optptr->fstin);
    printf("Parametres du fichier standard qui definissent le champ voulu\n");
    printf("     nomvar = '%s'\n", optptr->fst.nomvar);
    printf("     typvar = '%s'\n", optptr->fst.typvar);
    printf("     etiket = '%s'\n", optptr->fst.etiket);
    printf("     datev  =  %d\n", optptr->fst.dateo);
    printf("     ip1    =  %d\n", optptr->fst.ip1);
    printf("     ip2    =  %d\n", optptr->fst.ip2);
    printf("     ip3    =  %d\n", optptr->fst.ip3);
    printf("\n");
  } /* Fin du 'else' relie au 'if ( optptr->roundrobin == 1 )' */

  return 0;
}  /* Fin de la fonction parseOptions */


/***************************************************************************
 * fonction: aide
 *
 * Cette fonction imprime une documentation sommaire de l'utilisation de ce programme.   
 ***************************************************************************/
void aide(void) {

  printf("Ce programme %s (version '%s') permet d'extraire les observations d'une base de donnees\n", PROGRAM_NAME, VERSION);
  printf("qui sont a l'interieur du domaine d'une grille definie par une grille d'un fichier standard\n\n");
  
  printf("Les arguments pour ce programme sont:\n");
  printf("  %s        [fichier standard dans lequel on va chercher le champ voulu]\n\n", FSTIN_OPTION);

  printf("  %s        [fichier d'input dans lequel observations sont contenues]\n", OBSIN_OPTION);
  printf("  %s       [fichier dans lequel seront stockees les observations selectionnees]\n\n", OBSOUT_OPTION);

  printf("     Le fichier d'input peut etre un fichier BURP, RDB (SQLite) ou ASCII (avec un format precis)\n");

  printf("  %s  [On utilise la methode 'round-robin' pour separer les enregistrements d'un fichier BURP ou RDB en parties egales.]\n\n", ROUNDROBIN_OPTION);

  printf("  %s        [On decide si on prend les observations a l'interieur (valeur = 1)\n", INOUT_OPTION);
  printf("                 ou a l'exterieur (valeur = 0) du domaine. ] (par defaut, %d, a l'interieur du domaine)\n\n", INOUT_DEFAUT);

  printf("  %s        [On decide si on utilise une zone tampon]\n", PILOT_OPTION);
  printf("                 (par defaut, %d points a l'interieur du domaine)\n\n", PILOT_DEFAUT);

  printf("On choisit une portion de la grille qui definit le domaine.  Pour la definir, on utilise les options suivantes\n");
  printf("  %s        [indice 'i' minimal pour la grille (peut etre un nombre reel)]\n", MIN_I_OPTION);
  printf("  %s        [indice 'i' maximal pour la grille (peut etre un nombre reel)]\n", MAX_I_OPTION);
  printf("  %s        [indice 'j' minimal pour la grille (peut etre un nombre reel)]\n", MIN_J_OPTION);
  printf("  %s        [indice 'j' maximal pour la grille (peut etre un nombre reel)]\n", MAX_J_OPTION);
  printf("         Note: Ces dernieres options sont mutuellement exclusives avec l'option %s\n\n", PILOT_OPTION);

  printf("On peut separer les fichiers en plusieurs bandes de latitude et de longitude avec les options suivantes\n");
  printf("  %s         [nombre de bandes de longitudes (selon i)]\n", NPEX_OPTION);
  printf("  %s         [nombre de bandes de latitudes  (selon j)]\n\n", NPEY_OPTION);

  printf("On peut extraire seulement un seul fichier d'une operation de splitting plutot que de les extraire tous en même temps.\n");
  printf("  %s            [coordonnee en 'x' (selon i)]\n", CHERRYPICK_X_OPTION);
  printf("  %s            [coordonnee en 'y' (selon j)]\n\n", CHERRYPICK_Y_OPTION);

  printf("  %s   [active la verification a chaque enregistrement si on est en presence d'un UA multi-niveaux (ua4d)]\n\n", CHECK_UA4D_OPTION);

  printf("Attention, les options '%s' et '%s' puis '%s' et '%s' viennent en couple"
	 " et sont mutuellement exclusives\n\n", RDBIN_OPTION, RDBOUT_OPTION, BURPIN_OPTION, BURPOUT_OPTION);
  
  printf("On specifie les tables et la cle primaire qui lie les tables ensemble a l'aide des options suivantes\n");
  printf("  %s     [table qui joue le role du 'header' (par defaut '%s')]\n", RDB_HEADER_OPTION, RDB_HEADER_DEFAUT);
  printf("  %s       [table qui joue le role du 'data' (par defaut '%s')]\n", RDB_DATA_OPTION, RDB_DATA_DEFAUT);
  printf("  %s [cle primaire qui lie les tables 'header' et 'data' (par defaut '%s')]\n\n", RDB_PRIMARYKEY_OPTION, RDB_PRIMARYKEY_DEFAUT);

  printf("Les options suivantes servent a identifier le champ qui definit la grille (fonction fstinf):\n\n");

  printf("  %s    [nomvar du champ qu'on veut lire] (par defaut, PN)\n", NOMVAR_OPTION);
  printf("  %s    [typvar du champ qu'on veut lire] (par defaut, vide)\n", TYPVAR_OPTION);
  printf("  %s    [etiket du champ qu'on veut lire] (par defaut, vide)\n", ETIKET_OPTION);
  printf("  %s     [date de validite du champ qu'on veut lire]\n", DATEV_OPTION);
  printf("                   (format: YYYYMMDDHHMMSS avec padding de 0 par la droite, si la date est incomplete:\n");
  printf("                           par exemple, 2004121912 devient 20041219120000)\n");
  printf("                           (par defaut, -1)\n");
  printf("  %s       [ip1 du champ qu'on veut lire] (par defaut, -1)\n", IP1_OPTION);
  printf("  %s       [ip2 du champ qu'on veut lire] (par defaut, -1)\n", IP2_OPTION);
  printf("  %s       [ip3 du champ qu'on veut lire] (par defaut, -1)\n\n", IP3_OPTION);

  printf("On peut aussi filtrer verticalement en imposant un niveau avec l'option\n");
  printf("  %s  [niveau en hPa pour lequel toute observation au dessus sera enlevee] (defaut, -1 donc aucun filtrage vertical)\n", NIVEAU_MAX_OPTION);
  printf("  %s  [niveau en hPa pour lequel toute observation en dessous sera enlevee] (defaut, -1 donc aucun filtrage vertical)\n", NIVEAU_MIN_OPTION);
  printf("  %s          [fichier standard qui contiendra un champ GZ qui donnera la hauteur avec "
	 "laquelle les observations seront filtrees verticalement] (par defaut, aucun filtrage vertical)\n", GZ_OPTION);
  printf("\n");
  printf("Les observations du type 'ai' et 'ua' ont une coordonnee verticale en pression (hPa) donc la selection est directe.  \n");
  printf("Les observations du type 'pr' et 'ro' ont une coordonnee verticale en metres, on donne alors\n");
  printf("    un critere vertical (hPa) en pression et un champ GZ a lire dans le fichier standard specifie avec l'option %s\n", GZ_OPTION);

  printf("Dans le cas des observations satellitaires, on specifie plutot les canaux voulus avec l'option\n");
  printf("  %s   [liste separee par des virgules, sans espace aucun, des canaux voulus (pour un maximum de %d caracteres en tout)]\n", CHANNELS_OPTION, MAXSTR_CHANNELS);
  printf("ou bien plutot ceux que l'on ne veut pas avec l'option\n");
  printf("  %s [liste separee par des virgules, sans espace aucun, des canaux exclus (pour un maximum de %d caracteres en tout)]\n", NOCHANNELS_OPTION, MAXSTR_CHANNELS);

  printf("\n");

  printf("  %s [degre de verbosite (par defaut, 0)]\n", VERBOSE_OPTION);

  printf("  %s ou %s ou %s [affiche cette aide]\n\n", HELP_OPTION1, HELP_OPTION2, HELP_OPTION3);

  printf("Ervig Lapalme, CMDA\n");

} /* Fin de la fonction aide */

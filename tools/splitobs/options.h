#ifndef __INCLUDE_OPTIONS_H__
#define __INCLUDE_OPTIONS_H__

/* valeurs par defaut pour certaines structures utilisees dans le programme */
#define MAXSTR_CHANNELS    MAXSTR
#define NOMVAR_DEFAUT      "PN"
#define INOUT_DEFAUT       1
#define PILOT_DEFAUT       0
#define NPEX_DEFAUT        1
#define NPEY_DEFAUT        1
#define NDIGITS_DEFAUT     4
#define RECTANGLE_DEFAUT   {-1,-1,-1,-1,0,0,0,0}
#define RDB_HEADER_DEFAUT     "header"
#define RDB_DATA_DEFAUT       "data"
#define RDB_SPLITONKEY_DEFAUT "id_obs"
#define optionsDEFAUT      {""                    /* fstin */, \
                            ""                    /* obsin*/, \
                            ""                    /* obsout*/, \
                            ""                    /* gz*/, \
                            ""                    /* channels */, \
                            INOUT_DEFAUT          /* inout */, \
                            PILOT_DEFAUT          /* pilot */, \
                            RECTANGLE_DEFAUT      /* rect */, \
                            fstparam_DEFAUT       /* fst */, \
                            IP1_VIDE              /* niveau_min */, \
                            IP1_VIDE              /* niveau_max */, \
                            1                     /* channels_voulus */, \
                            NPEX_DEFAUT           /* npex */, \
                            NPEY_DEFAUT           /* npey */, \
                            NDIGITS_DEFAUT        /* ndigits */, \
                            0                     /* check_ua4d */, \
                            0                     /* roundrobin */, \
                            -1                    /* cherrypick_x */, \
                            -1                    /* cherrypick_y */, \
                            0                     /* numheaders_files */, \
                            RDB_HEADER_DEFAUT     /* rdb_header_table */, \
                            RDB_DATA_DEFAUT       /* rdb_data_table */, \
                            RDB_SPLITONKEY_DEFAUT /* rdb_split_on_key */ }

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
#define RDB_SPLITONKEY_OPTION  "-split-on-key"
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
#define NUMHEADERS_FILES_OPTION "-numheaders-files"
#define VERBOSE_OPTION     "-verbose"
#define CHECK_UA4D_OPTION  "-check_ua4d"
#define HELP_OPTION1       "-h"
#define HELP_OPTION2       "-help"
#define HELP_OPTION3       "--help"
#define AIDE_OPTION1       "-aide"
#define AIDE_OPTION2       "--aide"

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
  int niveau_min, niveau_max, channels_voulus;
  int npex, npey, ndigits;
  int check_ua4d;
  int roundrobin;
  int cherrypick_x, cherrypick_y;
  int numheaders_files;
  char rdb_header_table[MAXSTR], rdb_data_table[MAXSTR], rdb_split_on_key[MAXSTR];
} options, *optionsptr;

int parseOptions(int argc, char** argv, optionsptr optptr, int* VERBOSE);

#endif /* #ifndef __INCLUDE_OPTIONS_H__ */

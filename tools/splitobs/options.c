#include <stdio.h>

/* Include pour la librairie App dans 'rpn/libs' */
#include "App.h"

/* Include pour ma librairie de manipulation des fichiers standard RPN */
#include "fstdlib.h"
/* Include pour les constantes OK et NOT_OK */
#include "ok_or_notok.h"

/* Include pour les definitions des fonctions definies dans ce fichier */
#include "options.h"

void aide(int VERBOSE);
void help(int VERBOSE);

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
 *
 * Cette fonction retourne:
 *            OK si tout s'est bien passe
 *            NOT_OK s'il y a une erreur
 *
 ***************************************************************************/
int parseOptions(int argc, char** argv, optionsptr optptr, int* VERBOSE) {
  int i;

  if (argc==1) { /* Alors, on veut la documentation */
    help(*VERBOSE);
    exit(0);
  }

  strcpy(optptr->fst.nomvar,NOMVAR_DEFAUT);

  for (i=1;i<argc;i++) {
    if (argv[i][0]=='-') {
      if (!strcmp(argv[i],HELP_OPTION1) ||
          !strcmp(argv[i],HELP_OPTION2) ||
          !strcmp(argv[i],HELP_OPTION3)) { /* Alors, on veut la documentation */
        help(*VERBOSE);
        exit(0);
      }
      if (!strcmp(argv[i],AIDE_OPTION1) ||
          !strcmp(argv[i],AIDE_OPTION2)) { /* Alors, on veut la documentation */
        aide(*VERBOSE);
        exit(0);
      }
      else if (!strcmp(argv[i],VERBOSE_OPTION)) { /* On indique le niveau de print */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR, "Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          return NOT_OK;
        }
        *VERBOSE = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],FSTIN_OPTION)) { /* Alors, on donne le nom du fichier standard input */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->fstin,argv[++i]);
      }
      else if (!strcmp(argv[i],OBSIN_OPTION)) { /* Alors, on donne le nom du fichier d'observations input */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->obsin,argv[++i]);
      }
      else if (!strcmp(argv[i],OBSOUT_OPTION)) { /* Alors, on donne le nom du fichier d'observations output */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->obsout,argv[++i]);
      }
      else if (!strcmp(argv[i],RDBIN_OPTION)) { /* Alors, on donne le nom du fichier RDB input */
        App_Log(APP_WARNING,"Fonction parseOptions: L'option %s a ete remplacee par %s\n", argv[i], OBSIN_OPTION);
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->obsin,argv[++i]);
      }
      else if (!strcmp(argv[i],RDBOUT_OPTION)) { /* Alors, on donne le nom du fichier RDB output */
        App_Log(APP_WARNING,"Fonction parseOptions: L'option %s a ete remplacee par %s\n", argv[i], OBSOUT_OPTION);
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->obsout,argv[++i]);
      }
      else if (!strcmp(argv[i],BURPIN_OPTION)) { /* Alors, on donne le nom du fichier BURP input */
        App_Log(APP_WARNING,"Fonction parseOptions: L'option %s a ete remplacee par %s\n", argv[i], OBSIN_OPTION);
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->obsin,argv[++i]);
      }
      else if (!strcmp(argv[i],BURPOUT_OPTION)) { /* Alors, on donne le nom du fichier BURP output */
        App_Log(APP_WARNING,"Fonction parseOptions: L'option %s a ete remplacee par %s\n", argv[i], OBSOUT_OPTION);
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->obsout,argv[++i]);
      }
      else if (!strcmp(argv[i],ASCII_OPTION)) { /* Alors, on donne le nom du fichier ASCII qui contient des couplets lat-lon de points */
        App_Log(APP_WARNING,"Fonction parseOptions: L'option %s a ete remplacee par %s\n", argv[i], OBSIN_OPTION);
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->obsin,argv[++i]);
      }
      else if (!strcmp(argv[i],OUT_OPTION)) { /* Alors, on donne le nom du fichier ASCII de sortie */
        App_Log(APP_WARNING,"Fonction parseOptions: L'option %s a ete remplacee par %s\n", argv[i], OBSOUT_OPTION);
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->obsout,argv[++i]);
      }
      else if (!strcmp(argv[i],INOUT_OPTION)) { /* On decide si on prend les observations a l'interieur (valeur = 1)
                                                 * ou a l'exterieur (valeur = 0) du domaine.
                                                 */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->inout = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],PILOT_OPTION)) { /* On decide si on utilise une zone de pilotage pour rapetisser le domaine
                                                 */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->pilot = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],MIN_I_OPTION)) { /* On donne le 'i' minimal pour definir le domaine
                                                 */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
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
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
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
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
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
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
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
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->npex = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],NPEY_OPTION)) { /* On donne le nombre de bandes dont on veut separer en 'j'
                                                */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->npey = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],CHERRYPICK_X_OPTION)) { /* On donne la bande voulue en 'i' */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->cherrypick_x = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],CHERRYPICK_Y_OPTION)) { /* On donne la bande voulue en 'j' */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->cherrypick_y = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],NDIGITS_OPTION)) { /* On donne le nombre de catacteres que l'on veut
                                                   * avoir les extensions des fichiers separes en sous-domaines.
                                                   */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->ndigits = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],CHECK_UA4D_OPTION)) { /* Par defaut, on assume que l'on ne melange pas les enregistrements
                                                      * ua4d et ua classiques alors on ne fait la verification a chaque
                                                      * enregistrement.  Si cette option est activee, alors on fera cette
                                                      * verification a chaque enregistrement.
                                                      */
        if (i+1<argc && argv[i+1][0]!='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s ne demande aucun argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->check_ua4d = 1;
      }
      else if (!strcmp(argv[i],ROUNDROBIN_OPTION)) { /* Cette option activera le mode round-robin pour faire un splitting
                                                      * en fichiers sans egard a la position geographique.
                                                      */
        if (i+1<argc && argv[i+1][0]!='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s ne demande aucun argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->roundrobin = 1;
      }
      else if (!strcmp(argv[i],NUMHEADERS_FILES_OPTION)) { /* Cette option indique que l'on veut en sortie les fichiers
                                                            * '*.num_headers' et '*.max_num_headers'.
                                                            */
        if (i+1<argc && argv[i+1][0]!='-') {
          fprintf(stderr,"Fonction parseOptions: L'option %s ne demande aucun argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->numheaders_files = 1;
      }
      else if (!strcmp(argv[i],RDB_HEADER_OPTION)) { /* Cette option enregistrera le nom de la table qui joue le role du 'header' */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->rdb_header_table,argv[++i]);
      }
      else if (!strcmp(argv[i],RDB_DATA_OPTION)) { /* Cette option enregistrera le nom de la table qui joue le role du 'data' */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->rdb_data_table,argv[++i]);
      }
      else if (!strcmp(argv[i],RDB_SPLITONKEY_OPTION)) { /* Cette option enregistrera le nom de la cle qui sera utilisee pour faire le split des fichiers */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->rdb_split_on_key,argv[++i]);
      }
      else if (!strcmp(argv[i],GZ_OPTION)) { /* Fichier standard dans lequel on lira le GZ au niveau voulu */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR, "Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          return NOT_OK;
        }
        strcpy(optptr->gz,argv[++i]);
      }
      else if (!strcmp(argv[i],NIVEAU_MIN_OPTION)) { /* On donne le niveau en hPa voulu pour accepter
                                                      * les observations au dessus
                                                      */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR, "Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          return NOT_OK;
        }
        optptr->niveau_min = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],NIVEAU_MAX_OPTION)) { /* On donne le niveau en hPa voulu pour accepter
                                                      * les observations en dessous
                                                      */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR, "Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          return NOT_OK;
        }
        optptr->niveau_max = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],CHANNELS_OPTION) || !strcmp(argv[i],NOCHANNELS_OPTION))  {
        /* On specifie les canaux voulus ou que l'on ne veut pas */
        int indice_option = i;

        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR, "Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
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
            App_Log(APP_ERROR,"Fonction parseOptions: L'option %s ne peut prendre qu'un argument "
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
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->fst.nomvar,argv[++i]);
      }
      else if (!strcmp(argv[i],TYPVAR_OPTION)) { /* On a besoin du typvar du champ */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->fst.typvar,argv[++i]);
      }
      else if (!strcmp(argv[i],ETIKET_OPTION)) { /* On a besoin de l'etiquette du champ */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        strcpy(optptr->fst.etiket,argv[++i]);
      }
      else if (!strcmp(argv[i],DATEV_OPTION)) { /* On a besoin de la date du champ */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->fst.dateo = padtime(argv[++i]);
      }
      else if (!strcmp(argv[i],IP1_OPTION)) { /* On a besoin du ip1 du champ */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->fst.ip1 = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],IP2_OPTION)) { /* On a besoin du ip2 du champ */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->fst.ip2 = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i],IP3_OPTION)) { /* On a besoin du ip3 du champ */
        if (i+1>=argc || argv[i+1][0]=='-') {
          App_Log(APP_ERROR,"Fonction parseOptions: L'option %s demande un argument\n", argv[i]);
          App_End(-1);exit(1);
        }
        optptr->fst.ip3 = atoi(argv[++i]);
      }
      else {
        App_Log(APP_ERROR,"Fonction parseOptions: Je ne reconnais pas cette option: %s\n", argv[i]);
        App_End(-1);exit(1);
      }
    } /* Fin du if (argv[i][0] ...) */
    else {
      App_Log(APP_ERROR,"Fonction parseOptions: Erreur dans les arguments d'entree: (%d)\n\n", i);
      App_End(-1);exit(1);
    }
  } /* Fin du for (i=1; ...) */

  /* On checke les options recueillies */
  if ( strlen(optptr->fstin)==0 && optptr->roundrobin == 0 ) {
    App_Log(APP_ERROR,"Fonction parseOptions: On doit absolument specifier un fichier standard "
            "en entree avec l'option %s\n", FSTIN_OPTION);
    App_End(-1);exit(1);
  }

  printf("On imprime les differentes options choisies pour cette execution\n");

  if ( strlen(optptr->obsin) !=0 && strlen(optptr->obsout) !=0 ) {
    printf("Fichier d'observations en entree:  %s\n", optptr->obsin);
    printf("Fichier d'observations en sortie:  %s\n", optptr->obsout);
  }
  else {
    if ( strlen(optptr->obsin) == 0 ) {
      App_Log(APP_ERROR,"Fonction parseOptions: On doit absolument specifier soit un fichier de base de "
              "donnees SQL, un fichier BURP ou ASCII en entree avec l'option %s\n", OBSIN_OPTION);
    }
    if ( strlen(optptr->obsout) == 0 ) {
      App_Log(APP_ERROR,"Fonction parseOptions: On doit absolument specifier soit un fichier de base de "
              "donnees SQL, un fichier BURP ou ASCII en sortie avec l'option %s\n", OBSOUT_OPTION);
    }
    App_End(-1);exit(1);
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

    printf("Le domaine sera separe en %d par %d parties egales grace aux options (%s et %s)\n", optptr->npex, optptr->npey, NPEX_OPTION, NPEY_OPTION);
    snprintf(format_digits, sizeof(format_digits), "%%.%dd",optptr->ndigits);
    snprintf(npex_str, sizeof(npex_str), format_digits,optptr->npex);
    snprintf(npey_str, sizeof(npey_str), format_digits,optptr->npey);
    printf("Les extensions auront %d caracteres par exemple %s_%s_%s\n\n", optptr->ndigits, optptr->obsout, npex_str, npey_str);
  }

  if (optptr->cherrypick_x>0 && optptr->cherrypick_y>0) {
    char x_str[MAXSTR], y_str[MAXSTR], format_digits[MAXSTR];
    snprintf(format_digits, sizeof(format_digits), "%%.%dd",optptr->ndigits);
    snprintf(x_str, sizeof(x_str), format_digits,optptr->cherrypick_x);
    snprintf(y_str, sizeof(y_str), format_digits,optptr->cherrypick_y);
    printf("On ne va extraire que les observations pour la tuile '%s_%s'\n\n",x_str,y_str);
  }
  else if (!(optptr->cherrypick_x<0 && optptr->cherrypick_y<0)) {
    App_Log(APP_ERROR,"Fonction parseOptions: On doit absolument specifier les deux arguments %s et %s en meme temps.\n", CHERRYPICK_X_OPTION, CHERRYPICK_Y_OPTION);
    App_End(-1);exit(1);
  }

  if ( optptr->numheaders_files == 1 ) {
    printf("On creera les fichiers '*.num_headers' et '*.max_num_headers' qui indiquent combien d'entetes contient chaque fichier.\n");
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
        App_Log(APP_ERROR, "Fonction parseOptions: Si on donne des canaux avec l'option %s, on ne peut utiliser les options %s et %s\n", CHANNELS_OPTION, NIVEAU_MIN_OPTION, NIVEAU_MAX_OPTION);

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

  return OK;
}  /* Fin de la fonction parseOptions */


/***************************************************************************
 * fonction: aide
 *
 * Cette fonction imprime une documentation sommaire de l'utilisation de ce programme.
 ***************************************************************************/
void aide(int VERBOSE) {

  printf("Ce programme permet d'extraire les observations d'une base de donnees qui sont a l'interieur\n");
  printf("du domaine d'une grille definie par champ d'un fichier standard RPN.\n\n");

  printf("Il permet aussi de separer les observations en plusieurs fichiers.\n\n");

  printf("Les arguments pour ce programme sont:\n");
  printf("  %s        [fichier d'input dans lequel observations sont contenues]\n", OBSIN_OPTION);
  printf("  %s       [fichier dans lequel seront stockees les observations selectionnees]\n\n", OBSOUT_OPTION);

  printf("     Le fichier d'input peut etre un fichier BURP, RDB (SQLite) ou ASCII (avec un format precis)\n\n");

  printf("  %s  [Utilisation de la methode 'round-robin' pour separer les enregistrements d'un fichier BURP ou RDB en parties egales]\n\n", ROUNDROBIN_OPTION);

  printf("  %s  [Activation de la sortie des fichiers '*.num_headers' et '*.max_num_headers' qui indiquent combien d'entetes contient chaque fichier.]\n\n", NUMHEADERS_FILES_OPTION);

  printf("  %s        [fichier standard RPN d'entree dans lequel on va chercher le champ voulu]\n\n", FSTIN_OPTION);

  printf("  %s        [On decide si on prend les observations a l'interieur (si egal a '1') ou a l'exterieur (valeur = 0) du domaine]\n", INOUT_OPTION);
  printf("                 (par defaut, %d, a l'interieur du domaine)\n\n", INOUT_DEFAUT);

  printf("  %s        [On decide si on utilise une zone tampon (par defaut, %d points a l'interieur du domaine)]\n\n", PILOT_OPTION, PILOT_DEFAUT);

  printf("On choisit une portion de la grille qui definit le domaine.  Pour la definir, on utilise les options suivantes:\n");
  printf("  %s        [indice 'i' minimal pour la grille (peut etre un nombre reel)]\n", MIN_I_OPTION);
  printf("  %s        [indice 'i' maximal pour la grille (peut etre un nombre reel)]\n", MAX_I_OPTION);
  printf("  %s        [indice 'j' minimal pour la grille (peut etre un nombre reel)]\n", MIN_J_OPTION);
  printf("  %s        [indice 'j' maximal pour la grille (peut etre un nombre reel)]\n", MAX_J_OPTION);
  printf("         Note: Ces dernieres options sont mutuellement exclusives avec l'option %s\n\n", PILOT_OPTION);

  printf("On peut separer les fichiers en plusieurs fichiers selon 'x' et 'y' avec les options suivantes\n");
  printf("  %s         [nombre de bandes de longitudes (selon i)]\n", NPEX_OPTION);
  printf("  %s         [nombre de bandes de latitudes  (selon j)]\n\n", NPEY_OPTION);

  printf("On peut extraire seulement un seul fichier d'une operation de splitting plutot que de les extraire tous en meme temps avec ces options:\n");
  printf("  %s            [numero de la composante selon 'x']\n", CHERRYPICK_X_OPTION);
  printf("  %s            [numero de la composante selon 'y']\n\n", CHERRYPICK_Y_OPTION);

  printf("Activation de la verification, a chaque enregistrement, si on est en presence d'un UA multi-niveaux (ua4d):\n");
  printf("  %s\n\n", CHECK_UA4D_OPTION);

  printf("Lorsqu'on utilise des fichiers SQLite, on utilise les options suivantes pour specifier ou aller chercher l'information necessaire:\n\n");

  printf("  %s [cle qui lie les tables entre elles et qui est aussi utilisee en mode round-robin\n", RDB_SPLITONKEY_OPTION);
  printf("                  pour selectionner les observations (par defaut '%s')]\n\n", RDB_SPLITONKEY_DEFAUT);

  printf("  %s       [table qui contient l'information sur la latitude (colonne 'lat'), la longitude (colonne 'lon') et\n", RDB_HEADER_OPTION);
  printf("                  l'elevation (colonne 'elev') de la station (si necessaire) (par defaut '%s')]\n\n", RDB_HEADER_DEFAUT);

  printf("  %s         [table qui contient l'information sur la hauteur de chaque observation ou \n", RDB_DATA_OPTION);
  printf("                   le canal satellitaire (colonne 'vcoord') (par defaut '%s')]\n\n", RDB_DATA_DEFAUT);

  printf("Les options suivantes servent a identifier le champ qui definit le domaine (fonction RMNLIB 'fstinf'):\n");
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

  printf("On peut aussi filtrer verticalement en utilisant ces criteres\n");
  printf("  %s  [niveau en hPa pour lequel toute observation au dessus sera enlevee] (defaut, -1 donc aucun niveau maximal)\n", NIVEAU_MAX_OPTION);
  printf("  %s  [niveau en hPa pour lequel toute observation en dessous sera enlevee] (defaut, -1 donc aucun niveau minimal)\n", NIVEAU_MIN_OPTION);
  printf("  %s          [fichier standard RPN qui contiendra un champ GZ qui donnera la hauteur avec\n", GZ_OPTION);
  printf("                   laquelle les observations seront filtrees verticalement]\n");
  printf("                   (doit etre specifie si une conversion de pression en hauteur necessaire)\n\n");

  printf("Les observations du type 'ai' et 'ua' ont une coordonnee verticale en pression (hPa) donc la selection est directe.  \n");
  printf("Les observations du type 'pr' et 'ro' ont une coordonnee verticale en metres, on donne alors\n");
  printf("un critere vertical (hPa) en pression et un champ GZ a lire dans le fichier standard specifie avec l'option %s\n\n", GZ_OPTION);

  printf("Dans le cas des observations satellitaires, on specifie plutot les canaux voulus avec l'option:\n");
  printf("  %s   [liste separee par des virgules, sans espace aucun, des canaux voulus (pour un maximum de %d caracteres en tout)]\n\n", CHANNELS_OPTION, MAXSTR_CHANNELS);

  printf("ou bien plutot ceux que l'on ne veut pas avec l'option:\n");
  printf("  %s [liste separee par des virgules, sans espace aucun, des canaux exclus (pour un maximum de %d caracteres en tout)]\n", NOCHANNELS_OPTION, MAXSTR_CHANNELS);

  printf("\n");

  printf("Degre de verbosite:\n");
  printf("  %s [par defaut, %d]\n\n", VERBOSE_OPTION, VERBOSE);

  printf("  %s ou %s ou [affiche cette aide]\n", AIDE_OPTION1, AIDE_OPTION2);
  printf("  %s or %s or %s [print an English help message]\n\n", HELP_OPTION1, HELP_OPTION2, HELP_OPTION3);

  printf("Ervig Lapalme, CMDA\n");

} /* Fin de la fonction 'aide' */


/***************************************************************************
 * fonction: 'help'
 *
 * This function output a documentation of the program
 ***************************************************************************/
void help(int VERBOSE) {

  printf("This programs extracts observations from a database which are inside \n");
  printf("the domain of a grid as defined by a field in a RPN standard file.\n\n");

  printf("It can also split the observations into several files.\n\n");

  printf("The arguments of this program are:\n");
  printf("  %s        [input observation file]\n", OBSIN_OPTION);
  printf("  %s       [output observation file]\n\n", OBSOUT_OPTION);

  printf("     The input file can be a BURP file, a RDB (SQLite) or ASCII (with a specific format)\n\n");

  printf("  %s  [Use the round-robin method to split the observations in equal parts]\n\n", ROUNDROBIN_OPTION);

  printf("  %s  [Activate the creation of the files '*.num_headers' et '*.max_num_headers' which count how many headers are contained in each file.]\n\n", NUMHEADERS_FILES_OPTION);

  printf("  %s        [input RPN standard file that contains the field which defines the domain]\n\n", FSTIN_OPTION);

  printf("  %s        [Decide whether observations inside the domain are chosen (if equals to '1') or outside the domain (if equals to '0')\n", INOUT_OPTION);
  printf("                 (default, %d, inside the domain)]\n\n", INOUT_DEFAUT);

  printf("  %s        [Decide to use a piloting buffer region (default, %d points inside the domain)]\n\n", PILOT_OPTION, PILOT_DEFAUT);

  printf("It is also possible to choose an arbitrary portion of the domain.  It can be defined using those options:\n");
  printf("  %s        [minimal index 'i' (along 'x' component) of the grid (can be a real number)]\n", MIN_I_OPTION);
  printf("  %s        [maximal index 'i' (along 'x' component) of the grid (can be a real number)]\n", MAX_I_OPTION);
  printf("  %s        [minimal index 'j' (along 'y' component) of the grid (can be a real number)]\n", MIN_J_OPTION);
  printf("  %s        [maximal index 'j' (along 'y' component) of the grid (can be a real number)]\n", MAX_J_OPTION);
  printf("         Note: Those options are mutually exclusive with the option %s\n\n", PILOT_OPTION);

  printf("We can split observations into several files in 'x' and 'y' directions with the following options:\n");
  printf("  %s         [number of files in the 'x' component]\n", NPEX_OPTION);
  printf("  %s         [number of files in the 'y' component]\n\n", NPEY_OPTION);

  printf("We can extract only a single tile instead of all using those options:\n");
  printf("  %s            [tile number along 'x' component]\n", CHERRYPICK_X_OPTION);
  printf("  %s            [tile number along 'y' component]\n\n", CHERRYPICK_Y_OPTION);

  printf("Activation of the verification for each record if it is a multi-level UA profile (ua4d):\n");
  printf("  %s\n\n", CHECK_UA4D_OPTION);

  printf("When processing SQLite files, we can use those options to specify where to find some information:\n\n");

  printf("  %s [key which binds the tables together which is also used in the round-robin mode\n", RDB_SPLITONKEY_OPTION);
  printf("                  to select the observations (default '%s')]\n\n", RDB_SPLITONKEY_DEFAUT);

  printf("  %s       [table which contains the latitude (column 'lat'), the longitude (column 'lon') and\n", RDB_HEADER_OPTION);
  printf("                  elevation (column 'elev') of the station if necessary (default '%s')]\n\n", RDB_HEADER_DEFAUT);

  printf("  %s         [table which contains the information about the height of each observation\n", RDB_DATA_OPTION);
  printf("                   or the satellite channel (column 'vcoord') (default '%s')]\n\n", RDB_DATA_DEFAUT);

  printf("Those options identify the field which defined the domain (RMNLIB function 'fstinf'):\n");
  printf("  %s    [nomvar of the wanted field] (default, PN)\n", NOMVAR_OPTION);
  printf("  %s    [typvar of the wanted field] (default, empty)\n", TYPVAR_OPTION);
  printf("  %s    [etiket of the wanted field] (default, empty)\n", ETIKET_OPTION);
  printf("  %s     [valid date of the wanted field]\n", DATEV_OPTION);
  printf("                   (format: YYYYMMDDHHMMSS padding with 0 at the right, if the date is incomplete\n");
  printf("                           for example, 2004121912 becomes 20041219120000)\n");
  printf("                           (by default, -1)\n");
  printf("  %s       [ip1 of the wanted field] (default, -1)\n", IP1_OPTION);
  printf("  %s       [ip2 of the wanted field] (default, -1)\n", IP2_OPTION);
  printf("  %s       [ip3 of the wanted field] (default, -1)\n\n", IP3_OPTION);

  printf("We can also filter vertically by using those criteria:\n");
  printf("  %s  [vertical level in hPa for which any observation above is not selected] (default, -1 so no maxinum level)\n", NIVEAU_MAX_OPTION);
  printf("  %s  [vertical level in hPa for which any observation below is not selected] (default, -1 so no mininum level)\n", NIVEAU_MIN_OPTION);
  printf("  %s          [RPN standard file which contains the GZ field to convert the pressure into meters\n", GZ_OPTION);
  printf("                to use the pressure levels as vertical criteria] (must be specified pressure to height conversion is needed)\n\n");

  printf("The observations 'ai' and 'ua' have a pressure vertical coordinate (hPa) so the selection is straightforward.\n");
  printf("The observations 'pr' and 'ro' have a vertical coordinate in meters, so we need to give\n");
  printf("pressure criteria (hPa) and a GZ field to read from the file specify with the option %s\n\n", GZ_OPTION);

  printf("For satellite observations, we specify the wanted channels with the option:\n");
  printf("  %s   [comma separated list (without any space) of the wanted channels (for a maximum of %d total characters)]\n\n", CHANNELS_OPTION, MAXSTR_CHANNELS);

  printf("or the channels that are not wanted with the option:\n");
  printf("  %s [comma separated list (without any space) of the channels that are not wanted (for a maximum of %d total characters)]\n", NOCHANNELS_OPTION, MAXSTR_CHANNELS);

  printf("\n");

  printf("Verbosity level:\n");
  printf("  %s [by default, %d]\n\n", VERBOSE_OPTION, VERBOSE);

  printf("  %s or %s or %s [print this help message]\n", HELP_OPTION1, HELP_OPTION2, HELP_OPTION3);
  printf("  %s ou %s ou [affiche cette aide en francais]\n\n", AIDE_OPTION1, AIDE_OPTION2);

  printf("Ervig Lapalme, CMDA\n");

} /* End of function 'help' */

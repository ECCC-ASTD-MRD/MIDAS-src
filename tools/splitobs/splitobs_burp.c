#include <stdio.h>
#include <unistd.h> /* to get the function 'access' to know if files are already existing */

/* Include pour la librairie App dans 'rpn/libs' */
#include "App.h"
/* Include pour la librairie RPN */
#include "rmn.h"
/* Include pour la librairie de manipulation de fichiers BURP*/
#include <burp_api.h>
extern int c_mrfbfl(int); /* This defnition is missing from 'burp_api.h' */

/* Include pour ma librairie de manipulation des fichiers standard RPN */
#include "fstdlib.h"
/* Include pour les constantes OK et NOT_OK */
#include "ok_or_notok.h"
/* Include pour la structure qui definit toutes les options */
#include "options.h"
/* Include pour les fonctions qui permettent de voir si un point dans a l'interieur de la grille de definition */
#include "checkgrid.h"

/**************************************************************/
/******* Prototype des fonctions definies localement **********/
/**************************************************************/

static int which_btyp(int btyp, int VERBOSE);
static int btypAssociated(int btyp_obs, int btyp, int VERBOSE);

static int clipping_vertical(BURP_RPT *rptin, optionsptr optptr, gridtype* grid_gz, BURP_RPT *rptout,
                             float* valeurs_gz_min, float* valeurs_gz_max, int VERBOSE);

static int find_subdomain(int gridid, int ni, int nj, float lat, float lon, rectangle rect, int npex, int npey,
                          int* ilonband, int* jlatband, char errmsg[MAXSTR], int VERBOSE);
static int putblk_nt(BURP_RPT *rpt, BURP_BLK *blk, int* t_in_domain, int nt, int VERBOSE);
static int putblk_nval(BURP_RPT *rpt, BURP_BLK *blk, int* val_in_domain, int nval, int VERBOSE);

static int extract_data_in_domains_along_nt(optionsptr optptr, gridtype* gridptr, BURP_RPT *rptin,
                                            int elem_lat, int elem_lon, int* nts, int** t_in_domain_ptr,
                                            int VERBOSE);
static int extract_data_in_domains_along_nval(optionsptr optptr, gridtype* gridptr, BURP_RPT *rptin,
                                              int elem_lat, int elem_lon, BURP_BLK *blk_data,
                                              int* nvals_in_domain, int* val_in_domain, int VERBOSE);
static int check_ua4d(BURP_RPT *rptin, int VERBOSE);
static int find_blk_data_in_rpt(BURP_RPT *rptin, int elem_lat, int elem_lon,
                                int** bknos_data_ptr, int** btyps_data_ptr, int* nblks, int VERBOSE);
static int find_blk_data_flag_in_rpt(BURP_RPT *rptin, int elem_lat, int elem_lon, int bkno_data,
                                     BURP_BLK **blk_data_ptr, BURP_BLK **blk_flags_ptr,
                                     int* colonne_lat_ptr, int* colonne_lon_ptr, int VERBOSE);
static int fill_rptout_blk(BURP_RPT *rptin, BURP_RPT ** rptout, int* nts, int* t_in_domain,
                           int n, int cherrypick_x, int cherrypick_y, int npey, int VERBOSE);

/*****************************************************/
/******* Fin des prototype des fonctions *************/
/*****************************************************/

int splitobs_burp(options opt, gridtype grid, gridtype grid_gz,
                 float* valeurs_gz_min, float* valeurs_gz_max, int VERBOSE) {
  int i, iun = 200, iout = 201, status, EXIT_STATUS = 0;
  int i_obs_enrgs, i_enrgs, nombre_enregistrements, longueur_max_enregistrement, engrs_resume;
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
      App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d avec la fonction brp_SetOptChar\n", status);
      return NOT_OK;
    }

    /* Ouverture du fichier burp en entree */
    status = brp_open(iun,opt.obsin,"r");
    if ( status<0 ) {
      App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d avec la fonction brp_open sur le fichier '%s'\n", status, opt.obsin);

      return NOT_OK;
    }
  else if (status == 0) {
      /* Si le fichier BURP est vide alors on doit sortir */
      App_Log(APP_WARNING,"Fonction splitobs_burp: Il n'y a aucun enregistrement dans le fichier BURP %s\n",opt.obsin);

      status = brp_close(iun);
      if (status<0)
        App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);

      return NOT_OK;
    }

    /* Si aucun probleme alors 'status' est le nombre d'enregistrements du fichier BURP */
    nombre_enregistrements = status;
    adresses = (int*) malloc(sizeof(int)*nombre_enregistrements);
    if ( adresses == (int*) NULL) {
      App_Log(APP_ERROR,"Fonction splitobs_burp: Incapable d'allouer un vecteur de (int) de dimension %d\n", nombre_enregistrements);

      status = brp_close(iun);
      if (status<0)
        App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);

      return NOT_OK;
    }

    nts = (int*) malloc(opt.npex*opt.npey*sizeof(int));
    if ( nts == (int*) NULL) {
      App_Log(APP_ERROR,"Fonction splitobs_burp: Incapable d'allouer le vecteur nts de (int) de dimension %dx%d\n", opt.npex, opt.npey);

      status = brp_close(iun);
      if (status<0)
        App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);

      free(adresses);

      return NOT_OK;
    } /* Fin du if ( nts == (int*) NULL) */

    num_obs_per_tile = (int*) malloc(opt.npex*opt.npey*sizeof(int));
    if ( num_obs_per_tile == (int*) NULL) {
      App_Log(APP_ERROR,"Fonction splitobs_burp: Incapable d'allouer le vecteur num_obs_per_tile de (int) de dimension %dx%d\n", opt.npex, opt.npey);

      status = brp_close(iun);
      if (status<0)
        App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);

      free(adresses);
      free(nts);

      return NOT_OK;
    } /* Fin du if ( nts == (int*) NULL) */

    iouts = (int*) malloc(opt.npex*opt.npey*sizeof(int));
    if ( iouts == (int*) NULL) {
      App_Log(APP_ERROR,"Fonction splitobs_burp: Incapable d'allouer le vecteur iouts de (int) de dimension %dx%d\n", opt.npex, opt.npey);

      status = brp_close(iun);
      if (status<0)
        App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);

      free(adresses);
      free(nts);
      free(num_obs_per_tile);

      return NOT_OK;
    } /* Fin du if ( iouts == (int*) NULL) */

    /* Ouverture du fichier burp en sortie */
    if ( opt.npex == 1 && opt.npey == 1 ) {

      status = access(opt.obsout,F_OK);
      if ( status==0 ) {
        App_Log(APP_ERROR,"Fonction splitobs_burp: Le fichier '%s' existe deja mais il pourrait etre efface "
        "alors il vaut mieux que ce fichier n'existe pas a l'appel du programme\n", opt.obsout);

        status = brp_close(iun);
      if (status<0)
        App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);

      free(adresses);
      free(nts);
      free(num_obs_per_tile);
      free(iouts);

        return NOT_OK;
      }

      iouts[0] = iout;
      status = brp_open(iout,opt.obsout,"a");
      if ( status<0 ) {
        App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d avec la fonction brp_open sur le fichier '%s'\n", status, opt.obsout);

        status = brp_close(iun);
        if (status<0)
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsin);

        status = brp_close(iout);
        if (status<0)
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close sur le fichier %s\n", status, opt.obsout);

        free(adresses);
        free(num_obs_per_tile);
        free(nts);
        free(iouts);

        return NOT_OK;
      }
    }
    else { /* Il faut ouvrir npex*npey fichiers */
      char npex_str[MAXSTR], npey_str[MAXSTR], format_digits[MAXSTR], burpout[MAXSTR*4];

      snprintf(format_digits, sizeof(format_digits), "%%.%dd",opt.ndigits);

      for (ilonband=0;ilonband<opt.npex;ilonband++) {
        for (jlatband=0;jlatband<opt.npey;jlatband++) {
          int id=ilonband*opt.npey+jlatband;

          /* Si on est en mode 'cherrypick', alors on ne considere que si la tuile est egale a celle voulue */
          if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0)
            if (ilonband != opt.cherrypick_x-1 || jlatband != opt.cherrypick_y-1)
              continue;

          iouts[id] = iout++;

          snprintf(npex_str, sizeof(npex_str), format_digits,ilonband+1);
          snprintf(npey_str, sizeof(npey_str), format_digits,jlatband+1);
          snprintf(burpout, sizeof(burpout),"%s_%s_%s", opt.obsout,npex_str,npey_str);

          if (VERBOSE>5)
            printf("Fonction splitobs_burp: appel de 'brp_open' sur le fichier '%s' pour id=%d\n", burpout, id);

          status = access(burpout,F_OK);
          if ( status==0 || iouts[id]>999 ) {
            int ilonbandtmp, jlatbandtmp;

            if ( status==0 )
              App_Log(APP_ERROR,"Fonction splitobs_burp: Le fichier '%s' existe deja mais il pourrait etre efface "
                      "alors il vaut mieux que ce fichier n'existe pas a l'appel du programme\n", burpout);
            else if ( iouts[id]>999 )
              App_Log(APP_ERROR,"Fonction splitobs_burp: Comme iouts[%d]=%d, la commande 'brp_open(...)' "
                      "n'acceptera pas cette valeur puisqu'elle est plus grande que 999\n", id, iouts[id]);

            status = brp_close(iun);
            if (status<0)
              App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close pour le fichier %s\n", status, opt.obsin);

            for (ilonbandtmp=0;ilonbandtmp<ilonband;ilonbandtmp++) {
              for (jlatbandtmp=0;jlatbandtmp<jlatband;jlatbandtmp++) {
                status = brp_close(iouts[ilonbandtmp*opt.npey+jlatbandtmp]);
                if ( status<0 )
                  App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close "
                                    "pour le fichier %s_%d_%d\n", status, opt.obsout,ilonbandtmp+1,jlatbandtmp+1);
              }
            }

            free(adresses);
            free(num_obs_per_tile);
            free(iouts);
            free(nts);

            return NOT_OK;
          } /* Fin du if ( status==0 ) pour l'acces au fichier */

          status = brp_open(iouts[id],burpout,"a");
          if ( status<0 ) {
            int ilonbandtmp, jlatbandtmp;

            App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d avec la fonction brp_open sur le fichier %s\n", status, burpout);

            status = brp_close(iun);
            if (status<0)
              App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close pour le fichier %s\n", status, opt.obsin);

            for (ilonbandtmp=0;ilonbandtmp<ilonband;ilonbandtmp++) {
              for (jlatbandtmp=0;jlatbandtmp<jlatband;jlatbandtmp++) {
                status = brp_close(iouts[ilonbandtmp*opt.npey+jlatbandtmp]);
                if ( status<0 )
                  App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close "
                                    "pour le fichier %s_%d_%d\n", status, opt.obsout,ilonbandtmp+1,jlatbandtmp+1);
              }
            }

            free(adresses);
            free(iouts);
            free(nts);
            free(num_obs_per_tile);

            return NOT_OK;
          } /* Fin du if ( status<0 ) */
        } /* Fin du for (jlatband=0;jlatband<opt.npey;jlatband++) */
      } /* Fin du for (ilonband=0;ilonband<opt.npex;ilonband++) */
    } /* Fin du else pour le if ( opt.npex == 1 && opt.npey == 1 ) */

    /* Toujours initialiser les pointeurs */
    rptout = (BURP_RPT**) malloc(opt.npex*opt.npey*sizeof(BURP_RPT*));
    if ( rptout == (BURP_RPT**) NULL ) {
      App_Log(APP_ERROR,"Fonction splitobs_burp: Incapable d'allouer le vecteur rptout de (BURP_RPT*) de dimension %dx%d\n", opt.npex, opt.npey);

      status = brp_close(iun);
      if (status<0)
        App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close pour le fichier %s\n", status, opt.obsin);

      for (ilonband=0;ilonband<opt.npex;ilonband++)
      for (jlatband=0;jlatband<opt.npey;jlatband++) {
        status = brp_close(iouts[ilonband*opt.npey+jlatband]);
        if ( status<0 )
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close "
                            "pour le fichier %s_%d_%d\n", status, opt.obsout,ilonband+1,jlatband+1);
      }

      free(adresses);
      free(iouts);
      free(nts);
      free(num_obs_per_tile);

      return NOT_OK;
    }

    rptin  = brp_newrpt();
    if ( rptin == (BURP_RPT*) NULL ) {
      App_Log(APP_ERROR,"Fonction splitobs_burp: Incapable d'allouer rptin de (BURP_RPT*)\n");
      status = brp_close(iun);
      if (status<0)
        App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close pour le fichier %s\n", status, opt.obsin);

      for (ilonband=0;ilonband<opt.npex;ilonband++)
      for (jlatband=0;jlatband<opt.npey;jlatband++) {
        status = brp_close(iouts[ilonband*opt.npey+jlatband]);
        if ( status<0 )
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close "
                            "pour le fichier %s_%d_%d\n", status, opt.obsout,ilonband+1,jlatband+1);
      }

      free(adresses);
      free(iouts);
      free(nts);
      free(num_obs_per_tile);
      free(rptout);

      return NOT_OK;
    }

    i_enrgs=0;
    RPT_SetHANDLE(rptin, 0);
    while ( brp_findrpt(iun, rptin) >= 0 )
      adresses[i_enrgs++] = RPT_HANDLE(rptin);

    longueur_max_enregistrement = c_mrfbfl(iun);
    if (VERBOSE>0)
      printf("Fonction splitobs_burp: les rapports auront la longueur %d\n", longueur_max_enregistrement);
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

    /* Cette variable contient le nombre d'enregistrements d'observation qui ont ete traites jusqu'a splitobs_burptenant dans la boucle */
    i_obs_enrgs = 0;
    /* Ensuite, on passe chaque enregistrement un a un */
    for (i_enrgs=0;i_enrgs<nombre_enregistrements;i_enrgs++) {
      engrs_resume = 0;

      status = brp_getrpt(iun,adresses[i_enrgs],rptin);
      if (status<0) {
        App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_getrpt pour "
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
        is_ua4d = check_ua4d(rptin, VERBOSE);
        if (is_ua4d<0) {
          App_Log(APP_ERROR,"Fonction splitobs_burp: la fonction 'check_ua4d' retourne %d\n"
                  "le fichier d'entree %s a l'adresse %d (%d rapport)\n",
                  is_ua4d, opt.obsin, adresses[i_enrgs], i_enrgs);
          continue;
        }
      }

      if (strncmp(">>",RPT_STNID(rptin),2)==0) { /* C'est un enregistrement resume */
        engrs_resume = 1;
        if (VERBOSE>0)
          printf("Fonction splitobs_burp: C'est un enregistrement resume '%s' (i_enrgs=%d)\n",
                 RPT_STNID(rptin),i_enrgs);
      }

      if (engrs_resume==0 && vertical_clipping==1) {
        BURP_RPT *rptin_clip_vert;

        rptin_clip_vert = brp_newrpt();
        brp_allocrpt(rptin_clip_vert, longueur_max_enregistrement);

        status = clipping_vertical(rptin,&opt,&grid_gz,rptin_clip_vert,valeurs_gz_min,valeurs_gz_max,VERBOSE);
        if (status<0) {
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction clipping_vertical\n", status);
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
          printf("Fonction splitobs_burp: appel de 'brp_putrpthdr' pour i=%d, iouts=%d et i_enrgs=%d\n",i,iouts[i],i_enrgs);

        status = brp_putrpthdr(iouts[i],rptout[i]);
        if (status<0) {
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_putrpthdr pour "
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
          printf("Fonction splitobs_burp: On copie cet enregistrement resume '%s' (i_enrgs=%d)\n",
                 RPT_STNID(rptin),i_enrgs);
        status = fill_rptout_blk(rptin,rptout,(int*) NULL,(int*) NULL,opt.npex*opt.npey,opt.cherrypick_x,opt.cherrypick_y,opt.npey,VERBOSE);
        if (status<0) {
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction fill_rptout_blk pour "
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
          printf("Fonction splitobs_burp: L'enregistrement %d sera place dans le fichier %d et contient %d observations.\n",
                 i_enrgs,bin_id,nts[bin_id]);

        num_obs_per_tile[bin_id] += nts[bin_id];
        status = fill_rptout_blk(rptin,rptout,nts,(int*) NULL,opt.npex*opt.npey,opt.cherrypick_x,opt.cherrypick_y,opt.npey,VERBOSE);
        if (status<0) {
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction fill_rptout_blk pour "
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
        status = extract_data_in_domains_along_nt(&opt,&grid,rptin,5002,6002,nts,&t_in_domain,VERBOSE);
        if (status<0) {
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction extract_data_in_domains_along_nt pour "
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
              printf("Fonction splitobs_burp: %d observations sont dans la bande %d pour cet enregistrement\n", nts[i], i);
          }
        if (obs_in_domain==0) {
          free(t_in_domain);
          t_in_domain = (int*) NULL;
          if (VERBOSE>1)
            printf("Fonction splitobs_burp: Aucune observation n'a ete selectionnee pour cet enregistrement\n");
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
              printf("Fonction splitobs_burp: On met a jour le header (RPT_ELEV(rptout[%d])) de %d a %d\n", i, RPT_ELEV(rptout[i]),nts[i]);
            }
            RPT_SetELEV(rptout[i],nts[i]);
            status = brp_updrpthdr(iouts[i],rptout[i]);
            if (status<0) {
              App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_updrpthdr pour "
                      "le fichier de sortie %s_%d_%d a l'adresse %d (%d rapport)\n",
                      status, opt.obsout, i/opt.npey+1, i%opt.npey+1, adresses[i_enrgs], i_enrgs);
              EXIT_STATUS = 1;
              break;
            }
          } /* Fin du if (RPT_ELEV(rptout[i]) != nts[i]) */
          else if (VERBOSE>5) {
            printf("Fonction splitobs_burp: On n'a pas pas besoin de mettre a jour le header (RPT_ELEV(rptout[%d])=%d)\n", i, RPT_ELEV(rptout[i]));
          }
        } /* Fin du for (i=0;i<opt.npex*opt.npey;i++) */

        if (EXIT_STATUS)
          break;

        status = fill_rptout_blk(rptin,rptout,nts,t_in_domain,opt.npex*opt.npey,opt.cherrypick_x,opt.cherrypick_y,opt.npey,VERBOSE);
        if (status<0) {
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction fill_rptout_blk pour "
                  "l'adresse %d (%d rapport)\n", status, adresses[i_enrgs], i_enrgs);
          EXIT_STATUS = 1;
          break;
        }
      } /* Fin du else if (strncmp("^",RPT_STNID(rptin),1)==0) */
      else if (is_ua4d==1) { /* Si c'est un 'ua4d' */
        int obs_in_domain=0, i_btyp, btyp_data, btypnum;
        int *bknos_data, *btyps_data, *nvals_in_domain, *val_in_domain;

        bknos_data = (int*) NULL;
        btyps_data = (int*) NULL;
        status = find_blk_data_in_rpt(rptin,5001,6001,&bknos_data,&btyps_data,&btypnum,VERBOSE);
        if (status<0) {
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur dans la fonction find_blk_data_in_rpt\n");
          if (bknos_data != (int*) NULL)  free(bknos_data);
          continue;
        }

        val_in_domain   = (int*) NULL;
        nvals_in_domain = (int*) NULL;
        nvals_in_domain = (int*) malloc(opt.npex*opt.npey*sizeof(int));
        if (nvals_in_domain == (int*) NULL) {
          App_Log(APP_ERROR,"Fonction splitobs_burp: Incapable d'allouer un vecteur de int de dimension %d pour le cas 'ua4d'\n", opt.npex*opt.npey);
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
            App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur dans la fonction brp_readblk pour bknos_data[%d]=%d\n", i_btyp, bknos_data[i_btyp]);
            brp_freeblk(blkdata);
            continue;
          }
          btyp_data = BLK_BTYP(blkdata);
          if (btyp_data != btyps_data[i_btyp])
            App_Log(APP_ERROR,"Fonction splitobs_burp: Probleme potentiel puisque btyp_data=%d est different de btyps_data[%d]=%d",
                    btyp_data, i_btyp, btyps_data[i_btyp]);

          val_in_domain = (int*) NULL;
          val_in_domain = (int*) malloc(opt.npex*opt.npey*BLK_NVAL(blkdata)*sizeof(int));
          if (val_in_domain == (int*) NULL) {
            App_Log(APP_ERROR,"Fonction splitobs_burp: Incapable d'allouer un vecteur de int de dimension %d pour le cas 'ua4d'\n", opt.npex*opt.npey*BLK_NVAL(blkdata));
            if (bknos_data != (int*) NULL) free(bknos_data);
            brp_freeblk(blkdata);
            continue;
          }
          for (i=0;i<opt.npex*opt.npey*BLK_NVAL(blkdata);i++)
            val_in_domain[i]=0;

          status = extract_data_in_domains_along_nval(&opt,&grid,rptin,5001,6001,blkdata,nvals_in_domain,val_in_domain,VERBOSE);
          if (status<0) {
            App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction extract_data_in_domains_along_nval pour "
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
                printf("Fonction splitobs_burp: %d observations sont dans la bande %d pour cet enregistrement\n", nvals_in_domain[i], i);
            }
          if (obs_in_domain==0) {
            free(val_in_domain);
            val_in_domain = (int*) NULL;
            if (VERBOSE>1)
              printf("Fonction splitobs_burp: Aucune observation n'a ete selectionnee pour cet "
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
              App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur dans la fonction brp_readblk pour bkno=%d\n", BLK_BKNO(blksearch));
              brp_freeblk(blkout);
              continue;
            }

            if ( btyp_data == BLK_BTYP(blkout) || btypAssociated(btyp_data,BLK_BTYP(blkout),VERBOSE) == 1 ) {
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
                      printf("Fonction splitobs_burp: Appel de putblk_nval btyp=%d i=%d nvals=%d BLK_NVAL(blkout)*i=%d\n",
                             BLK_BTYP(blkout), i, nvals_in_domain[i],BLK_NVAL(blkout)*i);

                    status = putblk_nval(rptout[i],blkout,&val_in_domain[BLK_NVAL(blkout)*i],nvals_in_domain[i],VERBOSE);
                    if (status<0) {
                      App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction putblk_nval pour "
                              "l'adresse %d (%d rapport) et le bloc %d (btyp=%d) pour la tuile %d (bloc data)\n",
                              status, adresses[i_enrgs], i_enrgs, i_btyp, BLK_BTYP(blkout), i);
                      EXIT_STATUS = 1;
                      break;
                    }
                  }
                } /* for (i=0;i<opt.npex*opt.npey;i++) */
              }
              else {
                App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur dans le traitement de l'enregistrement a "
                        "l'adresse %d (%d rapport) et le bloc %d (btyp=%d, nval=%d).  Le nval est different "
                        "du bloc d'observations associe (bkno=%d btyp=%d nval=%d)\n",
                        adresses[i_enrgs], i_enrgs, i_btyp, BLK_BTYP(blkout), BLK_NVAL(blkout),
                        BLK_BKNO(blkdata), BLK_BTYP(blkdata), BLK_NVAL(blkdata));
              }
            }/* Fin du 'if (btyp_data == BLK_BTYP(blkout) || btypAssociated(btyp_data,BLK_BTYP(blkout),VERBOSE) == 1)' */

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
            printf("Fonction splitobs_burp: on cherche les blocs qui n'ont pas ete traites\n");

          /* on trouve les blocs qui ne sont pas associes a aucun bloc */
          blktmp = brp_newblk();
          while ( brp_findblk( blktmp, rptin ) >= 0 ) {
            int associated = 0;
            BURP_BLK* blkout = (BURP_BLK*) NULL;

            blkout = brp_newblk();
            status = brp_readblk(BLK_BKNO(blktmp), blkout, rptin, 0);
            if (status<0) {
              App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur dans la fonction brp_readblk pour bknos_data[%d]=%d\n", i_btyp, BLK_BKNO(blktmp));
              brp_freeblk(blkout);
              continue;
            }

            for (i_btyp=0;i_btyp<btypnum;i_btyp++) {
              status = brp_readblk(bknos_data[i_btyp],blkdata,rptin,0);
              if (status<0) {
                App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur dans la fonction brp_readblk pour bknos_data[%d]=%d\n", i_btyp, bknos_data[i_btyp]);
                brp_freeblk(blkdata);
                continue;
              }

              if (BLK_NVAL(blkdata) == BLK_NVAL(blkout) && (btyps_data[i_btyp] == BLK_BTYP(blkout) || btypAssociated(btyps_data[i_btyp],BLK_BTYP(blkout),VERBOSE) == 1)) {
                if (VERBOSE>1)
                  printf("Fonction splitobs_burp: le btyp %d est associe a %d\n", BLK_BTYP(blkout), btyps_data[i_btyp]);
                associated=1;
                break;
              }
            }

            if (!associated) {
              if (VERBOSE>1)
                printf("Fonction splitobs_burp: le btyp %d n'a associe a aucun autre bloc d'observations\n", BLK_BTYP(blkout));

              /* Lorsqu'on en trouve qui n'a pas ete associe alors on le copie dans tous les blocs */
              for (i=0;i<opt.npex*opt.npey;i++) {
                if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0) {
                  int cherrypick_id = (opt.cherrypick_x-1)*opt.npey+opt.cherrypick_y-1;
                  if (i != cherrypick_id)
                    continue;
                }

                status = putblk_nval(rptout[i],blkout,(int*) NULL,0,VERBOSE);
                if (status<0) {
                  App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction putblk_nval pour "
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
        free(nvals_in_domain)/* Fin du 'if (btyp_data == BLK_BTYP(blkout) || btypAssociated(btyp_data,BLK_BTYP(blkout),VERBOSE) == 1)' */;
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
          status = checkgrid(grid.gridid, grid.ni, grid.nj, lat, lon, opt.rect, errmsg, VERBOSE);
          if (status<0) {
            App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur dans la fonction checkgrid pour le lat=%f "
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
                                  opt.npex, opt.npey, &ilonband, &jlatband, errmsg, VERBOSE);
          if (status<0) {
            App_Log(APP_ERROR,"Fonction main: Erreur dans la fonction find_subdosplitobs_burp "
                    "pour le lat=%f et lon=%f avec le message '%s'\n", lat, lon, errmsg);
            EXIT_STATUS = 1;
            break;
          }
          if (ilonband<1 || ilonband>opt.npex || jlatband<1 || jlatband>opt.npey) {
            if (VERBOSE>1)
              printf("Fonction splitobs_burp: cette observation n'est pas dans le domaine npex=%d, npey=%d "
                     "ilonband=%d jlatband=%d lat=%f lon=%f\n", opt.npex, opt.npey, ilonband, jlatband, lat, lon);
            continue;
          }

          if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0)
            if (ilonband != opt.cherrypick_x || jlatband != opt.cherrypick_y)
              continue;

          /* Il faut indiquer dans quel bande cette observation se trouve */
          id = (ilonband-1)*opt.npey+(jlatband-1);
          if (VERBOSE>4)
            printf("Fonction splitobs_burp: cette observation est acceptee: ilonband=%d jlonband=%d nts[%d]=%d\n",
                   ilonband, jlatband, (ilonband-1)*opt.npey+(jlatband-1), nts[(ilonband-1)*opt.npey+(jlatband-1)]);
        }

        nts[id] = 1;
        num_obs_per_tile[id]++; /* On ajoute 1 au compteur total des observations */

        status = fill_rptout_blk(rptin,rptout,nts,(int*) NULL,opt.npex*opt.npey,opt.cherrypick_x,opt.cherrypick_y,opt.npey,VERBOSE);
        if (status<0) {
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction fill_rptout_blk pour "
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
        App_Log(APP_ERROR,"Fonction splitobs_burp: il y a eu une erreur precedemment (dans else if (vertical_clipping==1))\n");
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
          printf("Fonction splitobs_burp (juste avant 'brp_writerpt'): engrs_resume=%d nts[%d]=%d\n",
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

          printf ( "Fonction splitobs_burp (juste avant 'brp_writerpt'): Entete du rapport rptout[%d]\n",i) ;
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
          printf("Fonction splitobs_burp: appel de 'brp_writerpt' pour i=%d, iouts=%d et i_enrgs=%d\n",i,iouts[i],i_enrgs);

        status = brp_writerpt(iouts[i],rptout[i],END_BURP_FILE);
        if (status<0) {
          App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_writerpt pour "
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
      App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close pour le fichier %s\n", status, opt.obsin);
      EXIT_STATUS = 1;
    }

    /* fermeture de fichier burp de sortie */
    if ( opt.npex == 1 && opt.npey == 1 ) {
      if (VERBOSE>2)
          printf("\nClosing BURP file %s iout = %d", opt.obsout, iout);

      printf("\nIl y a %d observations qui ont ete selectionnees\n", num_obs_per_tile[0]);

      status = brp_close(iout);
      if (status<0) {
        App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close pour le fichier %s\n", status, opt.obsout);
        EXIT_STATUS = 1;
      }
      else if (num_obs_per_tile[0]==0) {
        printf("Aucune observation n'a ete garde alors on efface le fichier %s\n", opt.obsout);
        status = remove(opt.obsout);
        if (status!=0) {
          App_Log(APP_ERROR,"Il est impossible d'effacer le fichier %s\n", opt.obsout);
          EXIT_STATUS = 1;
        }
        else if (VERBOSE>2)
          printf("On a efface le fichier BURP %s\n", opt.obsout);
      }
      else if (opt.numheaders_files == 1) { /* On imprime le nombre de headers presents dans le domaine */
        FILE* file;
        char burpout_num_headers[MAXSTR*2];

        snprintf(burpout_num_headers, sizeof(burpout_num_headers), "%s.num_headers", opt.obsout);

        status = access(burpout_num_headers,F_OK);
        if (status==0)
          App_Log(APP_ERROR,"Attention le fichier '%s' sera efface\n", burpout_num_headers);

        file = (FILE*) fopen(burpout_num_headers,"w");
        fprintf(file,"%d\n", num_obs_per_tile[0]);
        fclose(file);
      }
    }
    else {
      char npex_str[MAXSTR], npey_str[MAXSTR], format_digits[MAXSTR], burpout[MAXSTR*4];
      int max_num_headers=0;

      snprintf(format_digits, sizeof(format_digits), "%%.%dd",opt.ndigits);

      for (ilonband=0;ilonband<opt.npex;ilonband++)
        for (jlatband=0;jlatband<opt.npey;jlatband++) {
          int id=ilonband*opt.npey+jlatband;

          if (opt.cherrypick_x > 0 && opt.cherrypick_y > 0)
            if (ilonband != opt.cherrypick_x-1 || jlatband != opt.cherrypick_y-1)
              continue;

          snprintf(npex_str, sizeof(npex_str), format_digits,ilonband+1);
          snprintf(npey_str, sizeof(npey_str), format_digits,jlatband+1);
          snprintf(burpout, sizeof(burpout), "%s_%s_%s", opt.obsout,npex_str,npey_str);

          if (VERBOSE>2)
            printf("\nClosing correctly BURP file %s iouts[%d] = %d", burpout, id, iouts[id]);

          status = brp_close(iouts[id]);
          if ( status<0 )
            App_Log(APP_ERROR,"Fonction splitobs_burp: Erreur %d dans la fonction brp_close "
                    "pour le fichier %s\n", status, burpout);

          printf("\nIl y a %d observations qui ont ete selectionnees et mise dans le fichier %s\n",
                 num_obs_per_tile[id], burpout);

          if (num_obs_per_tile[id]==0) {
            printf("Aucune observation n'a ete garde alors on efface le fichier %s\n", burpout);
            status = remove(burpout);
            if (status!=0) {
              App_Log(APP_ERROR,"Il est impossible d'effacer le fichier %s\n", burpout);
              EXIT_STATUS = 1;
            }
            else if (VERBOSE>2)
              printf("On a efface le fichier BURP %s\n", burpout);
          }
          else if (opt.numheaders_files == 1) { /* On imprime le nombre de headers presents dans la tuile id */
            FILE* file;
            char burpout_num_headers[MAXSTR*6];

            snprintf(burpout_num_headers, sizeof(burpout_num_headers), "%s.num_headers", burpout);

            status = access(burpout_num_headers,F_OK);
            if (status==0)
              App_Log(APP_ERROR,"Attention le fichier '%s' sera efface\n", burpout_num_headers);

            file = (FILE*) fopen(burpout_num_headers,"w");
            fprintf(file,"%d\n", num_obs_per_tile[id]);
            fclose(file);
          }
          if (num_obs_per_tile[id]>max_num_headers) max_num_headers=num_obs_per_tile[id];
        }

      if (max_num_headers>0) {
        if (opt.numheaders_files == 1) { /* On imprime le nombre maximal de headers */
          FILE* file;
          char burpout_max_num_headers[MAXSTR*2];

          snprintf(burpout_max_num_headers, sizeof(burpout_max_num_headers), "%s.max_num_headers", opt.obsout);

          status = access(burpout_max_num_headers,F_OK);
          if (status==0)
            fprintf(stderr,"Attention le fichier '%s' sera efface\n", burpout_max_num_headers);

          file = (FILE*) fopen(burpout_max_num_headers,"w");
          fprintf(file,"%d\n", max_num_headers);
          fclose(file);
        }
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

    return EXIT_STATUS;
} /* Fin de la fonction 'int splitobs_sql' */



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
                   int* ilonband, int* jlatband, char errmsg[MAXSTR], int VERBOSE) {
  int status, criteria;
  float x, y;

  if (VERBOSE>3)
    printf("Fonction find_subdomain: lat=%f  lon=%f npex=%d npey=%d\n", lat, lon, npex, npey);

  if (lon<0) lon+=360;

  status = c_gdxyfll(gridid, &x, &y, &lat, &lon, 1);
  if (status<0) {
    snprintf(errmsg, MAXSTR,  "Fonction find_subdomain: Erreur avec c_gdxyfll qui retourne %d "
             "pour lat = %f, lon = %f, gridid = %d, ni = %d, nj = %d\n",
             status, lat, lon, gridid, ni, nj);
    App_Log(APP_ERROR,"%s",errmsg);
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

           status = c_ezgprm(gridid, input_grid.grtyp, &input_grid.ni, &input_grid.nj,
                               &input_grid.ig1, &input_grid.ig2, &input_grid.ig3, &input_grid.ig4);
           if (status<0) {
             snprintf(errmsg, MAXSTR,  "Fonction find_subdomain: Erreur avec c_ezgprm qui retourne %d "
                      "pour gridid = %d, ni = %d, nj = %d\n", status, gridid, ni, nj);
             App_Log(APP_ERROR,"%s",errmsg);
             return -1;
           }

           if (x<1)
             *ilonband = 1;
           else if (x>=ni)
             *ilonband = npex;
           else {
             /* reproduce the logic in the OAVAR code routine 'setObsMpiStrategy' */
             int ilonband_int = (((int) x) - 1)/(ni/npex)+1;

             *ilonband = ilonband_int;
           }

           if (*ilonband<1 || *ilonband>npex)
             App_Log(APP_ERROR,"Fonction find_subdomain: ilonband=%d n'est pas dans l'intervalle permis"
                      " pour lat=%f lon=%f x=%f y=%f npex=%d npey=%d ni=%d nj=%d\n",*ilonband,lat,lon,x,y,npex,npey,ni,nj);

           if (y<1)
             *jlatband = 1;
           else if (y>=nj)
             *jlatband = npey;
           else {
             int jlatband_int;

             /* On doit traiter le cas d'une grille gaussienne
              * differemment parce que les points prets du pole sont
              * soit <1 ou >nj
              */
             if (input_grid.grtyp[0]=='G') {
               /* reproduce the logic in the OAVAR code routine 'setObsMpiStrategy' */
               jlatband_int = ((int) y)/(nj/npey)+1;
             }
             else {
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
             App_Log(APP_ERROR,"Fonction find_subdomain: jlatband=%d n'est pas dans l'intervalle permis"
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
      App_Log(APP_ERROR,"Fonction find_subdomain: Le critere '%d' n'est pas possible.  Il faut 1,2, 3 ou 4.  ", criteria);

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
int which_btyp(int btyp, int VERBOSE) {
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
      App_Log(APP_ERROR,"Fonction which_btyp: Ce n'est ni un bloc marqueur ni un "
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
   *     . ssmuse-sh -d eccc/cmd/cmda/utils/${version}
   *
   ***************************************************************************/
int btypAssociated(int btyp_obs, int btyp, int VERBOSE) {
  int newbtyp_obs, newbtyp, bknat;

  if (btyp_obs == btyp) {
    App_Log(APP_ERROR,"Fonction btypAssociated: erreur: btyp_obs = btyp = %d\n", btyp);
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
int clipping_vertical(BURP_RPT *rptin, optionsptr optptr, gridtype* grid_gz, BURP_RPT *rptout,
                      float* valeurs_gz_min, float* valeurs_gz_max, int VERBOSE) {
  int e,v,t,status,trouve_data,trouve_marqueur,EXIT_STATUS = 0;
  int rangee_alt, rangee_lat, rangee_lon;
  BURP_BLK *blk, *blkout, *blk_donnees, *blk_marqueur, *new_blk_donnees, *new_blk_marqueur;

  brp_copyrpthdr(rptout,rptin);

  blk = (BURP_BLK*) NULL;
  blkout = (BURP_BLK*) NULL;

  blk = brp_newblk();

  BLK_SetBKNO(blk, 0);
  while ( brp_findblk( blk, rptin ) >= 0 ) {
    int is_data = 0, is_marqueur = 0, btyp, btypalt;

    blkout = brp_newblk();
    status = brp_readblk(BLK_BKNO(blk), blkout, rptin,0);
    if (status<0) {
      App_Log(APP_ERROR,"Fonction clipping_vertical: Erreur %d dans la fonction brp_readblk\n", status);
      EXIT_STATUS = -1;
      brp_freeblk(blkout);
      break;
    }

    btyp = BLK_BTYP(blkout);
    btypalt = btyp>>4 & 31;

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
      status = putblk_nt(rptout,blkout,(int*) NULL,0,VERBOSE);
      if (status!=0) {
        App_Log(APP_ERROR,"Fonction clipping_vertical: Erreur %d dans la fonction putblk_nt (btyp %d)\n",
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
        App_Log(APP_ERROR,"Fonction clipping_vertical: Les blocs data et marqueur "
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
        App_Log(APP_ERROR,"Fonction clipping_vertical: Erreur %d avec la fonction brp_encodeblk pour "
                "le bloc btyp=%d\n", status, BLK_BTYP(new_blk_donnees));
        EXIT_STATUS = 1;
        break;
      }

      status = brp_encodeblk(new_blk_marqueur);
      if (status<0) {
        App_Log(APP_ERROR,"Fonction clipping_vertical: Erreur %d avec la fonction brp_encodeblk pour "
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
              checkvertical(alt,optptr->niveau_min,optptr->niveau_max,VERBOSE))
            /* Le filtrage vertical est fait a l'aide d'une hauteur en pression */
            garde_rangee = 1;
          else if (strlen(optptr->channels)==0) {
            float lat = BLK_RVAL(blk_donnees,rangee_lat,v,t);
            float lon = BLK_RVAL(blk_donnees,rangee_lon,v,t);

            if(checkvertical_gz(lat,lon,alt,grid_gz->gridid,grid_gz->ni,grid_gz->nj,
                                optptr->niveau_min,optptr->niveau_max,valeurs_gz_min,valeurs_gz_max,VERBOSE))
              /* Le filtrage vertical est fait a l'aide d'une hauteur en metre */
              garde_rangee = 1;
          }
          else if (optptr->channels_voulus==checkcanal(alt,optptr->channels,VERBOSE))
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

      status = putblk_nt(rptout,new_blk_donnees,(int*) NULL,0,VERBOSE);
      brp_freeblk(new_blk_donnees);
      if (status!=0) {
        App_Log(APP_ERROR,"Fonction clipping_vertical: Erreur %d dans la fonction putblk_nt (btyp %d)\n",
                status, BLK_BTYP(new_blk_donnees));
        brp_freeblk(blkout);
        EXIT_STATUS = 1;
        break;
      }

      status = putblk_nt(rptout,new_blk_marqueur,(int*) NULL,0,VERBOSE);
      brp_freeblk(new_blk_marqueur);
      if (status!=0) {
        App_Log(APP_ERROR,"Fonction clipping_vertical: Erreur %d dans la fonction putblk_nt (btyp %d)\n",
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
int check_ua4d(BURP_RPT *rptin, int VERBOSE) {
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
  status = find_blk_data_in_rpt(rptin,5001,6001,&bknos_data,&btyps_data,&n_blk_data,VERBOSE);
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
                         int** bknos_data_ptr, int** btyps_data_ptr, int* nblks, int VERBOSE) {
  int i, e, status, nblkstmp, trouve_lat, trouve_lon;
  int trouve_au_moins_un_bloc_avec_lat_lon;
  BURP_BLK *blktmp, *blk;

  *nblks = 0; // initialisation par defaut

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
      App_Log(APP_ERROR,"Fonction find_blk_data_in_rpt: Erreur %d dans la fonction brp_readblk\n", status);
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
      printf("Fonction find_blk_data_in_rpt: which_btyp(%d)=%d\n", BLK_BTYP(blk), which_btyp(BLK_BTYP(blk), 0));
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
        if (which_btyp(BLK_BTYP(blk),VERBOSE)==0) {
          trouve_au_moins_un_bloc_avec_lat_lon=1;
          nblkstmp++;
          break;
        } /* Fin du if (which_btyp(BLK_BTYP(blk),VERBOSE)==0) */
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
    App_Log(APP_ERROR, "Fonction find_data_flag_in_rpt: Incapable d'allouer le vecteur '*bknos_data_ptr' de dimension %d de (int)\n", nblkstmp);
    brp_freeblk(blk);
    brp_freeblk(blktmp);
    return -1;
  }
  *btyps_data_ptr = (int*) NULL;
  *btyps_data_ptr = (int*) malloc(nblkstmp*sizeof(int));
  if (*btyps_data_ptr == (int*) NULL) {
    App_Log(APP_ERROR, "Fonction find_data_flag_in_rpt: Incapable d'allouer le vecteur '*btyps_data_ptr' de dimension %d de (int)\n", nblkstmp);
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
      App_Log(APP_ERROR,"Fonction find_blk_data_in_rpt: Erreur %d dans la fonction brp_readblk\n", status);
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
        if (which_btyp(BLK_BTYP(blk),VERBOSE)==0) {
          (*bknos_data_ptr)[nblkstmp] = BLK_BKNO(blk);
          (*btyps_data_ptr)[nblkstmp] = BLK_BTYP(blk);
          nblkstmp++;
          if (VERBOSE>5)
            printf("Fonction find_blk_data_in_rpt: On trouve btyp=%d au blkno=%d\n", BLK_BTYP(blk), BLK_BKNO(blk));
          break;
        } /* Fin du if (which_btyp(BLK_BTYP(blk),VERBOSE)==0) */
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
                              int* colonne_lat_ptr, int* colonne_lon_ptr, int VERBOSE) {
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
    App_Log(APP_ERROR,"Fonction find_blk_data_flag_in_rpt: Erreur %d dans la fonction brp_readblk "
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
      App_Log(APP_ERROR,"Fonction find_blk_data_flag_in_rpt: Erreur %d dans la fonction "
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
    App_Log(APP_ERROR,"Fonction find_blk_data_flag_in_rpt: le bloc data n'a pas ete trouve dans cet enregistrement\n");
    if (*blk_data_ptr != (BURP_BLK*) NULL) brp_freeblk(*blk_data_ptr);
    if (*blk_flags_ptr != (BURP_BLK*) NULL) brp_freeblk(*blk_flags_ptr);
    return -1;
  }

  if (trouve_flags == 0) {
    *blk_flags_ptr = (BURP_BLK*) NULL;
    if (VERBOSE>0)
      App_Log(APP_ERROR,"Fonction find_blk_data_flag_in_rpt: le bloc marqueur n'a pas ete trouve dans cet enregistrement\n");
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
int fill_rptout_blk(BURP_RPT *rptin, BURP_RPT ** rptout, int* nts, int* t_in_domain, int n,
                    int cherrypick_x,int cherrypick_y, int npey, int VERBOSE) {
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
      App_Log(APP_ERROR,"Fonction fill_rptout_blk: Erreur %d dans la fonction brp_readblk pour blk_no=%d\n", status, BLK_BKNO(blk));
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
        status = putblk_nt(rptout[i],blkout,&t_in_domain[BLK_NT(blkout)*i],nts[i],VERBOSE);
      else
        status = putblk_nt(rptout[i],blkout,(int*) NULL,0,VERBOSE);

      if (status!=0) {
        App_Log(APP_ERROR,"Fonction fill_rptout_blk: Erreur %d dans la fonction putblk_nt pour btyp=%d "
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
                                     int elem_lat, int elem_lon, int* nts, int** t_in_domain_ptr,
                                     int VERBOSE) {
  int status, t, colonne_lat, colonne_lon, EXIT_STATUS;
  int ilonband=1, jlatband=1;
  int btypnum, *bknos_data, *bknos_flags;
  float lat, lon;
  char errmsg[MAXSTR];
  BURP_BLK *blk_data, *blk_flags;

  EXIT_STATUS = 0;

  bknos_data = (int*) NULL;
  bknos_flags = (int*) NULL;
  status = find_blk_data_in_rpt(rptin,elem_lat,elem_lon,&bknos_data,&bknos_flags,&btypnum,VERBOSE);
  if (status<0) {
    App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nt: Erreur dans la fonction find_blk_data_in_rpt\n");
    return 1;
  }

  if (btypnum!=1) {
    App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nt: Dans la fonction find_blk_data_in_rpt, "
            "on a trouve %d blocs de donnees qui sont: ", btypnum);
    for (t=0;t<btypnum;t++)
      App_Log(APP_ERROR,"%d ", bknos_data[t]);
    App_Log(APP_ERROR,"\n");
    App_Log(APP_ERROR,"Or, ce programme ne peut traiter des enregistrements contenant plus d'un bloc de donnees regroupees\n");
    free(bknos_data);
    free(bknos_flags);
    return 1;
  }

  colonne_lat = -1;
  colonne_lon = -1;
  blk_data  = (BURP_BLK*) NULL;
  blk_flags = (BURP_BLK*) NULL;
  status = find_blk_data_flag_in_rpt(rptin,elem_lat,elem_lon,bknos_data[0],&blk_data,&blk_flags,&colonne_lat,&colonne_lon,VERBOSE);
  if (status<0) {
    App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nt: Erreur dans la fonction find_blk_data_flag_in_rpt\n");
    free(bknos_data);
    free(bknos_flags);
    return 1;
  }

  *t_in_domain_ptr = (int*) malloc(BLK_NT(blk_data)*optptr->npex*optptr->npey*sizeof(int));
  if (*t_in_domain_ptr == (int*) NULL) {
    App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nt: Incapable d'allouer la memoire pour un tableau de "
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
      status = checkgrid(gridptr->gridid, gridptr->ni, gridptr->nj, lat, lon, optptr->rect, errmsg, VERBOSE);
      if (status<0) {
        App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nt: Erreur dans la fonction checkgrid pour le lat=%f "
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
                              optptr->npex, optptr->npey, &ilonband, &jlatband, errmsg, VERBOSE);
      if (status<0) {
        App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nt: Erreur dans la fonction find_subdomain "
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

    /*         printf("nt = %d\n", nt); */
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
                                       int* nvals_in_domain, int* val_in_domain, int VERBOSE) {
  int status, id, e, v, colonne_lat, colonne_lon, EXIT_STATUS;
  int ilonband=1, jlatband=1;
  float lat, lon;
  char errmsg[MAXSTR];
  BURP_BLK *blk_data_converted = (BURP_BLK*) NULL;

  EXIT_STATUS = 0;

  blk_data_converted = brp_newblk();
  status = brp_readblk(BLK_BKNO(blk_data), blk_data_converted, rptin, 1);

  if (val_in_domain == (int*) NULL) {
    App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nval: le vecteur 'val_in_domain' doit etre deja alloue\n");
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
      App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nval: ne trouve pas "
              "les elements %d et %d dans l'entete du bloc\n", elem_lat, elem_lon);
    else if (colonne_lat<0)
      App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nval: ne trouve pas "
              "l'element %d dans l'entete du bloc\n", elem_lat);
    else if (colonne_lon<0)
      App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nval: ne trouve pas "
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
      status = checkgrid(gridptr->gridid, gridptr->ni, gridptr->nj, lat, lon, optptr->rect, errmsg, VERBOSE);
      if (status<0) {
        App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nval: Erreur dans la fonction checkgrid pour le lat=%f "
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
                              optptr->npex, optptr->npey, &ilonband, &jlatband, errmsg, VERBOSE);
      if (status<0) {
        App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nval: Erreur dans la fonction find_subdomain "
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
    App_Log(APP_ERROR,"Fonction extract_data_in_domains_along_nval: Une erreur dans la "
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
int putblk_nt(BURP_RPT *rpt, BURP_BLK *blk, int* t_in_domain, int nt, int VERBOSE) {
  int e,v,t,tt,status = 0;
  BURP_BLK *newblk;

  newblk = brp_newblk();

  if (VERBOSE>2)
    printf("Fonction putblk_nt: btyp=%d nt=%d blk_nele=%d blk_nval=%d blk_nt=%d t_in_domain=%p\n", BLK_BTYP(blk), nt, BLK_NELE(blk),BLK_NVAL(blk),BLK_NT(blk),(void*)t_in_domain);
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
    App_Log(APP_ERROR,"Fonction putblk_nt: Erreur %d avec la fonction brp_encodeblk pour le bloc %d\n", status, BLK_BKNO(blk));
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
    App_Log(APP_ERROR,"Fonction putblk_nt: Erreur %d dans la fonction brp_putblk (btyp %d, datyp=%d)\n",
            status, BLK_BTYP(newblk), BLK_DATYP(newblk));
    brp_freeblk(newblk);
    return -1;
  }

  brp_freeblk(newblk);
  if (VERBOSE>2)
    printf("Fonction putblk_nt: btyp=%d nt=%d blk_nele=%d blk_nval=%d blk_nt=%d t_in_domain=%p return=0\n",
           BLK_BTYP(blk), nt, BLK_NELE(blk),BLK_NVAL(blk),BLK_NT(blk),(void*)t_in_domain);

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
int putblk_nval(BURP_RPT *rpt, BURP_BLK *blk, int* vals_in_domain, int nval, int VERBOSE) {
  int e,v,vv,t,status = 0;
  BURP_BLK *newblk;

  newblk = brp_newblk();

  if (VERBOSE>2)
    printf("Fonction putblk_nval: btyp=%d nval=%d blk_nele=%d blk_nval=%d blk_nt=%d vals_in_domain=%p\n", BLK_BTYP(blk), nval, BLK_NELE(blk),BLK_NVAL(blk),BLK_NT(blk),(void*)vals_in_domain);
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
    App_Log(APP_ERROR,"Fonction putblk_nval: Erreur %d avec la fonction brp_encodeblk pour le bloc %d\n", status, BLK_BKNO(blk));
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
    App_Log(APP_ERROR,"Fonction putblk_nt: Erreur %d dans la fonction blk_putblk (btyp %d)\n", status, BLK_BTYP(blk));
    return -1;
  }
  brp_freeblk(newblk);

  if (VERBOSE>2)
    printf("Fonction putblk_nval: btyp=%d nval=%d blk_nele=%d blk_nval=%d blk_nt=%d vals_in_domain=%p return=0\n",
           BLK_BTYP(blk), nval, BLK_NELE(blk),BLK_NVAL(blk),BLK_NT(blk),(void*)vals_in_domain);

  return 0;
} /* Fin de la fonction putblk_nval */

#ifndef __INCLUDE_CHECKGRID_H__
#define __INCLUDE_CHECKGRID_H__

/* Include pour ma librairie de manipulation des fichiers standard RPN */
#include "fstdlib.h"
#include "options.h" /* To get access to the definition of 'rectangle' */

int checkgrid(int gridid, int ni, int nj, float lat, float lon, rectangle rect, char* errmsg, int VERBOSE);

int checkvertical(float vcoord, int niveau_min, int niveau_max, int VERBOSE);

int checkvertical_gz(float lat, float lon, float vcoord, int gridid, int ni, int nj, int niveau_min, int niveau_max,
                     float* VALEURS_GZ_MIN, float* VALEURS_GZ_MAX, int VERBOSE);

int checkcanal(float canal, char* channels, int VERBOSE);

#endif /* #ifndef __INCLUDE_CHECKGRID_H__ */

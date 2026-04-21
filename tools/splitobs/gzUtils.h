#ifndef __INCLUDE_GZUTILS_H__
#define __INCLUDE_GZUTILS_H__

/* Include pour la structure qui definit toutes les options */
#include "options.h"

void set_domain_rectangle(optionsptr optptr, gridtype grid);
int  getGridFromFile(fstparam* fstptr, char* fichier, gridtype* gridptr, int* file_handle);
int  getGZ(char* fichier, gridtype* gridptr, int niveau, float** valeurs);
int  getMinMaxGZ(optionsptr optptr, gridtype* gridptr,
                 float** valeurs_gz_min_ptr, float** valeurs_gz_max_ptr);

#endif /* #ifndef __INCLUDE_GZUTILS_H__ */

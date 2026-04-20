#ifndef __INCLUDE_SPLITOBS_SQL_H__
#define __INCLUDE_SPLITOBS_SQL_H__

/* Include pour la structure qui definit toutes les options */
#include "options.h"

int splitobs_sql(options opt, gridtype grid, gridtype grid_gz,
                 float* valeurs_gz_min, float* valeurs_gz_max, int VERBOSE);

#endif /* #ifndef __INCLUDE_SPLITOBS_SQL_H__ */

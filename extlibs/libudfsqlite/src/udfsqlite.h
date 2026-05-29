#ifndef _GUARD_UDFSQLITE_H_
#define	_GUARD_UDFSQLITE_H_
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "sqlite3.h"
#include <ctype.h>
#include <math.h>


# include <signal.h>
# include <pwd.h>
# include <unistd.h>
# include <sys/types.h>

#    ifdef    __cplusplus
#      define __BEGIN_DECLS extern "C" {
#      define __END_DECLS              }
#    else
#      define __BEGIN_DECLS /* -vide- */
#      define __END_DECLS   /* -vide- */
#    endif /* __cplusplus  */


 __BEGIN_DECLS 
 int RegisterFuncs(struct sqlite3 *db);
 __END_DECLS

#endif

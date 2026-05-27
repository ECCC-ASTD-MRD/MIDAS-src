/* csqlite.c --
      C wrappers callable from Fortran for the SQLite library
*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "sqlite3.h"
#include "udfsqlite.h"
#define DBG 0
#define FSQL_ALLOC_START 150000
#define FSQL_ALLOC_INCR  10000


static char *my_strncpy(char *dst, const char *src, size_t n) { /* {{{ 1 */
   char *d = dst;
   const char *s = src;
   if (n != 0) {
      do {

         if ((*d++ = *s++) == 0) {
            *(d -1 ) = ' ';
            while (--n != 0)
               *d++ = ' ';
            break;
         }
      } while (--n != 0);
   }
   return dst;
} /* 1}}} */
static int callback(void *NotUsed, int argc, char **argv, char **azColName){ /* {{{ 1 */
  /*
  int i;
  for(i=0; i<argc; i++){
    printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
  }
  printf("\n");
  */
  return 0;
} /* 1}}} */

int  sqlite3_open_c_(char *fname, sqlite3 **db, int len_fname) { /* {{{1 */
   int rc ;


if (DBG) printf("sqlite3_open [%s]\n", fname);
   rc = sqlite3_open(fname, db);
   RegisterFuncs(*db);
if (DBG) printf("sqlite3_open [%d]\n", rc);
   return rc ;
} /* 1}}} */
int  sqlite3_close_c_( sqlite3 **db) { /* {{{ 1 */
   int rc ;

   rc = sqlite3_close(*db);
   return rc ;
} /* 1}}} */

int  sqlite3_do_c_( sqlite3 **db, char *command, char *errmsg,int  len_command, int len_errmsg) { /* {{{ 1 */

    sqlite3_stmt *Stmt;
    const char *zLeftover;
   int   rc = 0  ;
/*   char *msg ;*/

if (DBG) printf("sql_do [%s]\n", command);
if (DBG) printf("sql_do [%d]\n", len_command);

    rc = sqlite3_prepare(*db, command, -1, &Stmt, &zLeftover);
if (DBG) printf("sql_do prepare rc=%d\n", rc);
    if (rc != SQLITE_OK) {
      strncpy( errmsg, sqlite3_errmsg( *db ) , len_errmsg ) ;
if (DBG) printf("sql_do prepare rc=%s\n", errmsg);
       return rc;
    }

    rc = sqlite3_step(Stmt);
if (DBG) printf("sql_do step rc=%d\n", rc);
    if ((rc != SQLITE_OK) && (rc != SQLITE_DONE)) {
      strncpy( errmsg, sqlite3_errmsg( *db ) , len_errmsg ) ;
      return rc;
    }

    rc = sqlite3_finalize(Stmt);
if (DBG) printf("sql_do finalize rc=%d\n", rc);
    if ((rc != SQLITE_OK) ) {
      strncpy( errmsg, sqlite3_errmsg( *db ) , len_errmsg ) ;
      return rc;
    }

    return rc;
} /* 1}}} */
int  sqlite3_do_many_c_( sqlite3 **db, char *command, char *errmsg, int len_command, int len_errmsg) { /* {{{1 */
   char *zErrMsg = 0;
   int   rc = 0  ;

if (DBG) printf("sql_do_many [%s]\n", command);
if (DBG) printf("sql_do_many [%d]\n", len_command);

    rc = sqlite3_exec(*db, command, 0, 0, &zErrMsg);

    if (rc != SQLITE_OK) {
      strncpy( errmsg, zErrMsg , len_errmsg ) ;
    }
if (DBG) printf("sql_do_many [%s]\n", zErrMsg);
    sqlite3_free(zErrMsg);

if (DBG) printf("sql_do_many rc[%d]\n", rc);
    return rc;
} /* 1 }}} */

int  sqlite3_prepare_c_( sqlite3 **db, char *command, sqlite3_stmt **stmt ,int len_command) { /* {{{1 */
   int   rc   ;
   const char *pstr ;
   rc = -1;

   rc = sqlite3_prepare( *db, command, (-1), stmt, &pstr ) ;
if (DBG) fprintf(stdout,"sqlite3_prepare_c called \n");
if (DBG) fprintf(stdout,"sqlite3_prepare_c error %d \n",rc);

   return rc ;
} /* 1}}} */
void  sqlite3_step_c_( sqlite3_stmt **stmt, int *completion) { /* {{{1 */
   *completion = sqlite3_step( *stmt ) ;

   return ;
} /* 1}}} */
void  sqlite3_reset_c_( sqlite3_stmt **stmt) { /* {{{1 */
   int   rc  ;

   rc = sqlite3_reset( *stmt ) ;

   return ;
} /* 1}}} */
void  sqlite3_finalize_c_( sqlite3_stmt **stmt) { /* {{{ 1 */
   int   rc  ;

   rc = sqlite3_finalize( *stmt ) ;

   return ;
} /* 1}}} */

void  sqlite3_errmsg_c_( sqlite3 **db, char *errmsg, int len_errmsg) { /* {{{1 */
   const char *pstr ;

   pstr = sqlite3_errmsg( *db ) ;
   strncpy( errmsg, pstr, len_errmsg ) ;

   return ;
} /* 1}}} */

void  sqlite3_column_count_c_( sqlite3_stmt **stmt, int *count) { /* {{{1 */
   *count = sqlite3_column_count( *stmt ) ;
   return ;
} /* 1}}} */
void  sqlite3_column_name_type_c_( sqlite3_stmt **stmt,int  *colidx,char *name,char *type,int len_name,int len_type) { /* {{{1 */
   const char *pstr ;

   pstr = sqlite3_column_name(*stmt, *colidx ) ;
/*            printf("%s\n", pstr);*/
   strncpy( name, pstr, len_name ) ;
   name[len_name-1] = '\0' ;
   pstr = sqlite3_column_decltype(*stmt, *colidx ) ;
/*            printf("%s\n", pstr);*/
   strncpy( type, (pstr?pstr: "NULL"), len_type ) ;
   type[len_type-1] = '\0' ;
   return ;
} /* 1}}} */

int  sqlite3_bind_int_c_( sqlite3_stmt **stmt, int *colidx, int  *value) { /* {{{1 */
   int   rc   ;

   rc = sqlite3_bind_int(*stmt, *colidx, *value ) ;
   return rc ;
} /* 1}}}} */
int  sqlite3_bind_int64_c_( sqlite3_stmt **stmt, int  *colidx, long long  *value) { /* {{{1 */
   int   rc   ;

   rc = sqlite3_bind_int64(*stmt, *colidx, *value ) ;
   return rc ;
} /* 1}}} */
int  sqlite3_bind_float_c_( sqlite3_stmt **stmt, int  *colidx, float *value) { /* {{{ 1 */
   int   rc   ;

   rc = sqlite3_bind_double(*stmt, *colidx, (double) *value ) ;
   return rc ;
} /* 1}}} */
int  sqlite3_bind_double_c_( sqlite3_stmt **stmt, int *colidx, double *value) { /* {{{1 */
   int   rc   ;

   rc = sqlite3_bind_double(*stmt, *colidx, *value ) ;
   return rc ;
} /* 1}}} */
int  sqlite3_bind_null_c_( sqlite3_stmt **stmt, int *colidx) { /* {{{ 1 */

   int   rc   ;

   rc = sqlite3_bind_null(*stmt, *colidx ) ;
   return rc ;
} /* 1}}} */
int  sqlite3_bind_text_c_( sqlite3_stmt **stmt, int *colidx, char *text, int  len_text) { /* {{{1 */
   int   rc   ;

   rc = sqlite3_bind_text(*stmt, *colidx, text, len_text,
           SQLITE_TRANSIENT ) ;
   return rc ;
} /* 1}}} */

int  sqlite3_column_int_c_( sqlite3_stmt **stmt, int *colidx, int  *value) { /* {{{1 */
   *value = sqlite3_column_int(*stmt, *colidx ) ;
if (DBG) fprintf(stdout,"sqlite3_column_int_c called \n");
   return 0 ;
} /* 1}}} */
int  sqlite3_column_intm_c_( sqlite3_stmt **stmt, int *colidx, int *value, int *miss) { /* {{{ 1 */
   *value  = (( sqlite3_column_text(*stmt, *colidx) !=NULL) ? 
      sqlite3_column_int(*stmt, *colidx) : *miss) ;
if (DBG) fprintf(stdout,"sqlite3_column_intm_c called \n");
   return 0 ;
} /* 1}}} */
int  sqlite3_column_int64_c_( sqlite3_stmt **stmt, int *colidx, long long *value) { /* {{{ 1 */

   *value = sqlite3_column_int64(*stmt, *colidx ) ;
if (DBG) fprintf(stdout,"sqlite3_column_int64_c called \n");
if (DBG) fprintf(stdout,"sqlite3_column_int64_c %lld \n",*value);
   return 0 ;
} /* 1}}} */
int  sqlite3_column_int64m_c_( sqlite3_stmt **stmt, int  *colidx, long long *value, long long  *miss) {/* {{{1 */

   *value  = (( sqlite3_column_text(*stmt, *colidx) !=NULL) ? 
      sqlite3_column_int64(*stmt, *colidx) : *miss) ;
if (DBG) fprintf(stdout,"sqlite3_column_int64m_c called \n");
if (DBG) fprintf(stdout,"sqlite3_column_int64m_c %lld \n",*value);
   return 0 ;
}/* 1}}} */
int  sqlite3_column_float_c_( sqlite3_stmt **stmt, int *colidx, float *value) { /* {{{1 */

   *value = (float) sqlite3_column_double(*stmt, *colidx ) ;
if (DBG) fprintf(stdout,"sqlite3_column_float_c called \n");
   return 0 ;
} /* 1}}} */
int  sqlite3_column_floatm_c_( sqlite3_stmt **stmt, int *colidx, float *value, float *miss) { /* {{{1 */

   *value  = (( sqlite3_column_text(*stmt, *colidx) !=NULL) ? 
      (float)sqlite3_column_double(*stmt, *colidx) : *miss) ;
if (DBG) fprintf(stdout,"sqlite3_column_floatm_c called \n");
   return 0 ;
}/* 1}}} */
int  sqlite3_column_double_c_( sqlite3_stmt **stmt, int *colidx, double *value) { /* {{{1 */

   *value = sqlite3_column_double(*stmt, *colidx ) ;
if (DBG) fprintf(stdout,"sqlite3_column_double_c called \n");
   return 0 ;
} /* 1}}} */
int  sqlite3_column_doublem_c_( sqlite3_stmt **stmt, int *colidx, double *value, double *miss) { /* {{{1 */

   *value  = (( sqlite3_column_text(*stmt, *colidx) !=NULL) ? 
      sqlite3_column_double(*stmt, *colidx) : *miss) ;
if (DBG) fprintf(stdout,"sqlite3_column_doublem_c called \n");
   return 0 ;
}/* 1}}} */
int  sqlite3_column_text_c_( sqlite3_stmt **stmt, int *colidx, char *text, int len_text) { /* {{{1 */

   const unsigned char *pstr ;

   pstr = sqlite3_column_text(*stmt, *colidx ) ;
   if(pstr != NULL) {
      strncpy( text, (const char * )pstr, len_text ) ;
   } else {
      strncpy( text, "NULL", len_text ) ;
   }
   return 0 ;
}/* 1}}}} */
int  sqlite3_column_textm_c_( sqlite3_stmt **stmt, int *colidx, char *text, char *miss, int len_text, int len_miss) { /* {{{1 */

   const unsigned char *pstr ;

   pstr = sqlite3_column_text(*stmt, *colidx ) ;
   if (pstr != NULL) {
      strncpy( text, (const char * )pstr, len_text ) ;
   } else {
      strncpy( text, miss, len_text ) ;
   }
   return 0 ;
} /* 1}}} */

void  sqlite3_get_many_intm_c_( sqlite3_stmt **stmt, int **a,int * nrows, int * ncols, int * completion, int * miss) { /* {{{1 */

   int rc, i,row,ROWS,COLS;
   int *iptr;
   int start = FSQL_ALLOC_START;
   int Step  = FSQL_ALLOC_INCR; 
   row = 0;
   ROWS = 0;
   COLS = sqlite3_column_count(*stmt);
   *completion = 0;

   iptr = (int *)malloc(start*COLS*sizeof(int));
      /*
      ** Tant qu'il ya des enregistrements
      */
      while ((rc = sqlite3_step(*stmt)) == SQLITE_ROW) {
            if (++ROWS > start) {
               start = start + Step;
               iptr = realloc(iptr, start*COLS*sizeof(int));
            };

	 /*
	 ** aller chercher les valeurs de chaque colonne de l'enregistrement
	 */
         for (i = 0; i < sqlite3_column_count(*stmt); ++i) {

             iptr[++row -1]  = (( sqlite3_column_text(*stmt, i) !=NULL) ? 
                sqlite3_column_int(*stmt, i) : *miss) ;
/*            printf("%s\t", sqlite3_column_text(*stmt, i));*/
         }

/*            printf("\n"); */
      }
      *a = iptr;
      *nrows = ROWS;
      *ncols = COLS;
} /* 1}}} */
void  sqlite3_get_many_int_c_( sqlite3_stmt **stmt, int **a,int * nrows, int * ncols, int * completion) { /* {{{1 */

   int rc, i,row,ROWS,COLS;
   int *iptr;
   int start = FSQL_ALLOC_START;
   int Step  = FSQL_ALLOC_INCR; 
   row = 0;
   ROWS = 0;
   COLS = sqlite3_column_count(*stmt);
   *completion = 0;

/*            printf(" csqlite.c : many_int_c\n"); */
   iptr = (int *)malloc(start*COLS*sizeof(int));
      /*
      ** Tant qu'il ya des enregistrements
      */
      while ((rc = sqlite3_step(*stmt)) == SQLITE_ROW) {
            if (++ROWS > start) {
               start = start + Step;
               iptr = realloc(iptr, start*COLS*sizeof(int));
            };

	 /*
	 ** aller chercher les valeurs de chaque colonne de l'enregistrement
	 */
         for (i = 0; i < sqlite3_column_count(*stmt); ++i) {

             iptr[++row -1]  =  sqlite3_column_int(*stmt, i) ;
/*            printf("%s\t", sqlite3_column_text(*stmt, i));*/
         }

/*            printf("\n"); */
      }
      *a = iptr;
      *nrows = ROWS;
      *ncols = COLS;
} /* 1}}} */

void  sqlite3_get_many_int8m_c_( sqlite3_stmt **stmt, long long **a,int *nrows, int *ncols, int *completion, long long *miss) { /* {{{1 */

   int rc, i,row,ROWS,COLS;
   long long *iptr;
   int start = FSQL_ALLOC_START;
   int Step  = FSQL_ALLOC_INCR; 
   row = 0;
   ROWS = 0;
   COLS = sqlite3_column_count(*stmt);
   *completion = 0;

   iptr = (long long *)malloc(start*COLS*sizeof(long long));
      /*
      ** Tant qu'il ya des enregistrements
      */
      while ((rc = sqlite3_step(*stmt)) == SQLITE_ROW) {
            if (++ROWS > start) {
               start = start + Step;
               iptr = realloc(iptr, start*COLS*sizeof(long long));
            };

	 /*
	 ** aller chercher les valeurs de chaque colonne de l'enregistrement
	 */
         for (i = 0; i < sqlite3_column_count(*stmt); ++i) {

             iptr[++row -1]  = (( sqlite3_column_text(*stmt, i) !=NULL) ? 
                sqlite3_column_int64(*stmt, i) : *miss) ;
/*            printf("%s\t", sqlite3_column_text(*stmt, i));*/
         }

/*            printf("\n"); */
      }
      *a = iptr;
      *nrows = ROWS;
      *ncols = COLS;
/*         for (i = 0; i < ROWS*COLS; i++) {*/
/*            if ((i%5) == 0) printf("\n");*/
/*            printf("%d\t", iptr[i]);*/
/*         }*/
} /* 1}}} */
void  sqlite3_get_many_int8_c_( sqlite3_stmt **stmt, long long **a,int *nrows, int *ncols, int *completion) { /* {{{1 */

   int rc, i,row,ROWS,COLS;
   long long *iptr;
   int start = FSQL_ALLOC_START;
   int Step  = FSQL_ALLOC_INCR; 
   row = 0;
   ROWS = 0;
   COLS = sqlite3_column_count(*stmt);
   *completion = 0;

   iptr = (long long *)malloc(start*COLS*sizeof(long long));
      /*
      ** Tant qu'il ya des enregistrements
      */
      while ((rc = sqlite3_step(*stmt)) == SQLITE_ROW) {
            if (++ROWS > start) {
               start = start + Step;
               iptr = realloc(iptr, start*COLS*sizeof(long long));
            };

	 /*
	 ** aller chercher les valeurs de chaque colonne de l'enregistrement
	 */
         for (i = 0; i < sqlite3_column_count(*stmt); ++i) {

             iptr[++row -1]  = sqlite3_column_int64(*stmt, i) ;
/*            printf("%s\t", sqlite3_column_text(*stmt, i));*/
         }

/*            printf("\n"); */
      }
      *a = iptr;
      *nrows = ROWS;
      *ncols = COLS;
/*         for (i = 0; i < ROWS*COLS; i++) {*/
/*            if ((i%5) == 0) printf("\n");*/
/*            printf("%d\t", iptr[i]);*/
/*         }*/
} /* 1}}} */

void  fsql_vector_int_c_( int * vec, int ** add,int * dim, int *completion) { /* {{{1 */
   int i;
   *completion = 0;

   if (*add == NULL) {
      if (DBG) printf("pointer is null \n");
      *completion = -1;
      return;
   }
   if (DBG) printf("dim = %d \n",*dim);
         for (i = 0; i < *dim; i++) {
/*            if ((i%5) == 0) printf("\n");*/
/*            printf("%d\t", iptr[i]);*/
              vec[i] = (*add)[i];
         }
} /* 1}}} */
void  fsql_vector_int8_c_( long long * vec, long long ** add,int * dim, int *completion) { /* {{{1 */
   int i;
   *completion = 0;

   if (*add == NULL) {
      if (DBG) printf("pointer is null \n");
      *completion = -1;
      return;
   }
   if (DBG) printf("dim = %d \n",*dim);
         for (i = 0; i < *dim; i++) {
/*            if ((i%5) == 0) printf("\n");*/
/*            printf("%d\t", iptr[i]);*/
              vec[i] = (*add)[i];
         }
} /* 1}}} */

void  sqlite3_get_many_realm_c_( sqlite3_stmt **stmt, float **a,int *nrows, int *ncols, int *completion, float *miss) { /* {{{1 */

   int rc, i,row,ROWS,COLS;
   float *iptr;
   int start = FSQL_ALLOC_START;
   int Step  = FSQL_ALLOC_INCR; 
   row = 0;
   ROWS = 0;
   COLS = sqlite3_column_count(*stmt);
   *completion = 0;

   iptr = (float *)malloc(start*COLS*sizeof(float));
   if(!iptr) {
      *completion = -1;
      return;
   }
      /*
      ** Tant qu'il ya des enregistrements
      */
      while ((rc = sqlite3_step(*stmt)) == SQLITE_ROW) {
            if (++ROWS > start) {
               start = start + Step;
               iptr = realloc(iptr, start*COLS*sizeof(float));
               if(!iptr) {
                  *completion = -1;
                  return;
               }
            };

	 /*
	 ** aller chercher les valeurs de chaque colonne de l'enregistrement
	 */
         for (i = 0; i < sqlite3_column_count(*stmt); ++i) {

             iptr[++row -1]  = (( sqlite3_column_text(*stmt, i) !=NULL) ? 
                (float)sqlite3_column_double(*stmt, i) : *miss) ;
/*            printf("%s\t", sqlite3_column_text(*stmt, i));*/
         }

/*            printf("\n"); */
      }
      *a = iptr;
      *nrows = ROWS;
      *ncols = COLS;
/*         for (i = 0; i < ROWS*COLS; i++) {*/
/*            if ((i%5) == 0) printf("\n");*/
/*            printf("%d\t", iptr[i]);*/
/*         }*/
} /* 1}}} */
void  sqlite3_get_many_real_c_( sqlite3_stmt **stmt, float **a,int *nrows, int *ncols, int *completion) { /* {{{1 */

   int rc, i,row,ROWS,COLS;
   float *iptr;
   int start = FSQL_ALLOC_START;
   int Step  = FSQL_ALLOC_INCR; 
   row  = 0;
   ROWS = 0;
   COLS = sqlite3_column_count(*stmt);
   *completion = 0;

   if (DBG) fprintf(stderr," csqlite.c : many_real_c\n"); 
   iptr = (float *)malloc(start*COLS*sizeof(float));
   if(!iptr) {
      *completion = -1;
      return;
   }
      /*
      ** Tant qu'il ya des enregistrements
      */
      while ((rc = sqlite3_step(*stmt)) == SQLITE_ROW) {
            if (++ROWS > start) {
               start = start + Step;
               iptr = realloc(iptr, start*COLS*sizeof(float));
               if(!iptr) {
                  *completion = -1;
                  return;
               }
            };

	 /*
	 ** aller chercher les valeurs de chaque colonne de l'enregistrement
	 */
         for (i = 0; i < sqlite3_column_count(*stmt); ++i) {

             iptr[++row -1]  =  (float)sqlite3_column_double(*stmt, i)  ;
/*            printf("%f\t", iptr[row -1]);*/
         }

/*            printf("\n"); */
      }
      *a = iptr;
      *nrows = ROWS;
      *ncols = COLS;
/*         for (i = 0; i < ROWS*COLS; i++) {*/
/*            if ((i%5) == 0) printf("\n");*/
/*            printf("%d\t", iptr[i]);*/
/*         }*/
} /* 1}}} */

void  sqlite3_get_many_real8m_c_( sqlite3_stmt **stmt, double **a,int *nrows, int *ncols, int *completion, double *miss) { /* {{{1 */

   int rc, i,row,ROWS,COLS;
   double *iptr;
   int start = FSQL_ALLOC_START;
   int Step  = FSQL_ALLOC_INCR; 
   row = 0;
   ROWS = 0;
   COLS = sqlite3_column_count(*stmt);
   *completion = 0;

   iptr = (double *)malloc(start*COLS*sizeof(double));
   if(!iptr) {
      *completion = -1;
      return;
   }
      /*
      ** Tant qu'il ya des enregistrements
      */
      while ((rc = sqlite3_step(*stmt)) == SQLITE_ROW) {
            if (++ROWS > start) {
               start = start + Step;
               iptr = realloc(iptr, start*COLS*sizeof(double));
               if(!iptr) {
                  *completion = -1;
                  return;
               }
            };

	 /*
	 ** aller chercher les valeurs de chaque colonne de l'enregistrement
	 */
         for (i = 0; i < sqlite3_column_count(*stmt); ++i) {

             iptr[++row -1]  = (( sqlite3_column_text(*stmt, i) !=NULL) ? 
                sqlite3_column_double(*stmt, i) : *miss) ;
/*            printf("%s\t", sqlite3_column_text(*stmt, i));*/
         }

/*            printf("\n"); */
      }
      *a = iptr;
      *nrows = ROWS;
      *ncols = COLS;
/*         for (i = 0; i < ROWS*COLS; i++) {*/
/*            if ((i%5) == 0) printf("\n");*/
/*            printf("%d\t", iptr[i]);*/
/*         }*/
} /* 1}}} */
void  sqlite3_get_many_real8_c_( sqlite3_stmt **stmt, double **a,int *nrows, int *ncols, int *completion) { /* {{{1 */

   int rc, i,row,ROWS,COLS;
   double *iptr;
   int start = FSQL_ALLOC_START;
   int Step  = FSQL_ALLOC_INCR; 
   row = 0;
   ROWS = 0;
   COLS = sqlite3_column_count(*stmt);
   *completion = 0;

   iptr = (double *)malloc(start*COLS*sizeof(double));
   if(!iptr) {
      *completion = -1;
      return;
   }
      /*
      ** Tant qu'il ya des enregistrements
      */
      while ((rc = sqlite3_step(*stmt)) == SQLITE_ROW) {
            if (++ROWS > start) {
               start = start + Step;
               iptr = realloc(iptr, start*COLS*sizeof(double));
               if(!iptr) {
                  *completion = -1;
                  return;
               }
            };

	 /*
	 ** aller chercher les valeurs de chaque colonne de l'enregistrement
	 */
         for (i = 0; i < sqlite3_column_count(*stmt); ++i) {

             iptr[++row -1]  =  sqlite3_column_double(*stmt, i)  ;
/*            printf("%s\t", sqlite3_column_text(*stmt, i));*/
         }

/*            printf("\n"); */
      }
      *a = iptr;
      *nrows = ROWS;
      *ncols = COLS;
/*         for (i = 0; i < ROWS*COLS; i++) {*/
/*            if ((i%5) == 0) printf("\n");*/
/*            printf("%d\t", iptr[i]);*/
/*         }*/
} /* 1}}} */

void  fsql_vector_real_c_( float * vec, float ** add,int * dim, int *completion) { /* {{{ 1 */

   int i;
   *completion = 0;

   if (*add == NULL) {
      if (DBG) printf("pointer is null \n");
      *completion = -1;
      return;
   }
   if (DBG) printf("dim = %d \n",*dim);
         for (i = 0; i < *dim; i++) {
/*            if ((i%5) == 0) printf("\n");*/
/*            printf("%d\t", iptr[i]);*/
              vec[i] = (*add)[i];
         }
} /* 1}}} */
void  fsql_vector_real8_c_( double * vec, double ** add,int * dim, int *completion) { /* {{{1 */

   int i;
   if (*add == NULL) {
      if (DBG) printf("pointer is null \n");
      *completion = -1;
      return;
   }
   if (DBG) printf("dim = %d \n",*dim);
         for (i = 0; i < *dim; i++) {
/*            if ((i%5) == 0) printf("\n");*/
/*            printf("%d\t", iptr[i]);*/
              vec[i] = (*add)[i];
         }
} /* 1}}} */

void  sqlite3_get_many_char_c_( sqlite3_stmt **stmt, char ***a,int *nrows, int *ncols, int *completion, char *miss, int len_miss) { /* {{{1 */

   int rc, i,row,ROWS,COLS;
   char **iptr;
   int start = FSQL_ALLOC_START;
   int Step  = FSQL_ALLOC_INCR; 
   char **sp;
   char * z;
   row = 0;
   ROWS = 0;
   COLS = sqlite3_column_count(*stmt);
   *completion = 0;

   if (DBG) printf("missing = %s et l = %d \n", miss,len_miss); 
   iptr =malloc(start*COLS*sizeof(char*));
   if(!iptr) {
      *completion = -1;
      return;
   }
      /*
      ** Tant qu'il ya des enregistrements
      */
   while ((rc = sqlite3_step(*stmt)) == SQLITE_ROW) {
      if (++ROWS > start) {
         start = start + Step;
         iptr = realloc(iptr, start*COLS*sizeof(char*));
         if(!iptr) {
            *completion = -1;
            return;
         }
      };

      /*
      ** aller chercher les valeurs de chaque colonne de l'enregistrement
      */
      for (i = 0; i < sqlite3_column_count(*stmt); ++i) {

         z = (char*)sqlite3_column_text(*stmt, i);
         if (z != NULL) {
            iptr[++row -1]  = malloc(strlen(z) +1);
            strcpy(iptr[row -1], z);
         } else {
            iptr[++row -1]  = malloc(len_miss);
            strcpy(iptr[row -1], miss);
         }
/*            printf("%s\t", sqlite3_column_text(*stmt, i));*/
/*            printf("ointer i %p\n",iptr[row-1]); */
      }

/*            printf("\n"); */
   }
   *a = iptr;
   *nrows = ROWS;
   *ncols = COLS;
/*         for (i = 0; i < ROWS*COLS; i++) {*/
/*            if ((i%5) == 0) printf("\n");*/
/*            printf("%s\t", iptr[i]);*/
/*            printf("%p\t", iptr[i]);*/
/*         }*/
} /* 1}}} */
void  fsql_vector_char_c_( char *  vec, char *** add,int * dim, int *completion, int len) { /* {{{1 */

/*   printf("len vec  = %d\n",len);*/
/*   printf("vec  = %s\n",vec);*/
/*   vec[0] = 'A';*/
   int i;
   *completion = 0;

   char * sp = vec;
   printf("len vec  = %d\n",len);
  if (*add == NULL) {
      printf("pointer is null \n");
      *completion = -1;
      return;
   }
/*   printf("pointer vector_c = %p \n",*add);*/
/*   printf("dim = %d \n",*dim);*/
         for (i = 0; i < *dim; i++) {
/*              vec[i] = (*add)[i];*/
              my_strncpy(sp,(*add)[i],len);
              
/*            sp[len -1] = '\0';*/
/*   printf("s = %s \n",(*add)[i]);*/
/*   printf("sc = %s \n",sp);*/
            sp = sp+len;
         }
} /* 1}}} */

void  fsql_free_int_c_( int ** add ) { /* {{{1 */

   if (*add == NULL) return;
   free(*add);
   *add = NULL;
   if (DBG) printf("memory freed\n");
} /* 1}}} */
void  fsql_free_int8_c_( long long ** add ) { /* {{{1 */

   if (*add == NULL) return;
   free(*add);
   *add = NULL;
   if(DBG) printf("memory freed\n");
}/* 1}}} */
void  fsql_free_real_c_( float ** add ) { /* {{{1 */

   if (*add == NULL) return;
   free(*add);
   *add = NULL;
   if(DBG) printf("memory freed\n");
} /* 1}}} */
void  fsql_free_real8_c_( double ** add ) { /* {{{1 */

   if (*add == NULL) return;
   free(*add);
   *add = NULL;
   if(DBG) printf("memory freed\n");
} /* 1}}} */
void  fsql_free_char_c_( char *** add, int * dim ) { /* {{{1 */

   int i = 0;
   if (*add == NULL) return;
   for(;i != (*dim); i++)
   {
      if (DBG) printf("freeing %p\n",(*add)[i]);
      free ((*add)[i]);
      (*add)[i] = 0;
      if (DBG) printf("mow %p\n",(*add)[i]);
   }
   free(*add);
   *add = NULL;
   if (DBG) printf("memory freed\n");
   if (DBG) printf("dim = %d\n",*dim);
}/* 1}}} */
/*
** Return TRUE if the last non-whitespace character in z[] is a semicolon.
** z[] is N characters long.
*/
static int _ends_with_semicolon(const char *z, int N){
  while( N>0 && isspace((unsigned char)z[N-1]) ){ N--; }
  return N>0 && z[N-1]==';';
}

/*
** Test to see if a line consists entirely of whitespace.
*/
static int _all_whitespace(const char *z){
  for(; *z; z++){
    if( isspace(*(unsigned char*)z) ) continue;
    if( *z=='/' && z[1]=='*' ){
      z += 2;
      while( *z && (*z!='*' || z[1]!='/') ){ z++; }
      if( *z==0 ) return 0;
      z++;
      continue;
    }
    if( *z=='-' && z[1]=='-' ){
      z += 2;
      while( *z && *z!='\n' ){ z++; }
      if( *z==0 ) return 1;
      continue;
    }
    return 0;
  }
  return 1;
}

/*
** Return TRUE if the line typed in is an SQL command terminator other
** than a semi-colon.  The SQL Server style "go" command is understood
** as is the Oracle "/".
*/
static int _is_command_terminator(const char *zLine){
  while( isspace(*(unsigned char*)zLine) ){ zLine++; };
  if( zLine[0]=='/' && _all_whitespace(&zLine[1]) ) return 1;  /* Oracle */
  if( tolower(zLine[0])=='g' && tolower(zLine[1])=='o'
         && _all_whitespace(&zLine[2]) ){
    return 1;  /* SQL Server */
  }
  return 0;
}

/*
** This routine reads a line of text from FILE in, stores
** the text in memory obtained from malloc() and returns a pointer
** to the text.  NULL is returned at end of file, or if malloc()
** fails.
**
** The interface is like "readline" but no command-line editing
** is done.
*/
static char *local_getline(char *zPrompt, FILE *in){
  char *zLine;
  int nLine;
  int n;
  int eol;

  if( zPrompt && *zPrompt ){
    printf("%s",zPrompt);
    fflush(stdout);
  }
  nLine = 100;
  zLine = malloc( nLine );
  if( zLine==0 ) return 0;
  n = 0;
  eol = 0;
  while( !eol ){
    if( n+100>nLine ){
      nLine = nLine*2 + 100;
      zLine = realloc(zLine, nLine);
      if( zLine==0 ) return 0;
    }
    if( fgets(&zLine[n], nLine - n, in)==0 ){
      if( n==0 ){
        free(zLine);
        return 0;
      }
      zLine[n] = 0;
      eol = 1;
      break;
    }
    while( zLine[n] ){ n++; }
    if( n>0 && zLine[n-1]=='\n' ){
      n--;
      zLine[n] = 0;
      eol = 1;
    }
  }
  zLine = realloc( zLine, n+1 );
  return zLine;
}

/*
** Retrieve a single line of input text.  "isatty" is true if text
** is coming from a terminal.  In that case, we issue a prompt and
** attempt to use "readline" for command-line editing.  If "isatty"
** is false, use "local_getline" instead of "readline" and issue no prompt.
**
** zPrior is a string of prior text retrieved.  If not the empty
** string, then issue a continuation prompt.
*/
static char *one_input_line(const char *zPrior, FILE *in){
  if( in!=0 ){
    return local_getline(0, in);
  }
   return local_getline(0, in);
}

/*
** Read input from *in and process it.  If *in==0 then input
** is interactive - the user is typing it it.  Otherwise, input
** is coming from a file or device.  A prompt is issued and history
** is saved only if input is interactive.  An interrupt signal will
** cause this routine to exit immediately, unless input is interactive.
*/
int  sqlite3_do_script_c_( sqlite3 **db, char *fname, char *errmsg, int len_command, int len_errmsg) { /* {{{1 */
  char *zLine;
  char *zSql = 0;
  int nSql = 0;
  char *zErrMsg = 0;
  int rc;
  FILE * in;
  int count = 0;
  in = fopen(fname,"rb");
  if( in ){
      fprintf(stderr,"Executing the fSQl script %s\n",fname);
    }
  while( (zLine = one_input_line(zSql, in))!=0 ){
     count++;
    if( 1 ) printf("%s\n", zLine);
    if( (zSql==0 || zSql[0]==0) && _all_whitespace(zLine) ) continue;
    if( zLine && zLine[0]=='.' && nSql==0 ){
/*      int rc = do_meta_command(zLine, p);*/
/*      free(zLine);*/
/*      if( rc ) break;*/
      continue;
    }
    if( _is_command_terminator(zLine) ){
      strcpy(zLine,";");
    }
    if( zSql==0 ){
      int i;
      for(i=0; zLine[i] && isspace((unsigned char)zLine[i]); i++){}
      if( zLine[i]!=0 ){
        nSql = strlen(zLine);
        zSql = malloc( nSql+1 );
        strcpy(zSql, zLine);
      }
    }else{
      int len = strlen(zLine);
      zSql = realloc( zSql, nSql + len + 2 );
      if( zSql==0 ){
        fprintf(stderr,"%s: out of memory!\n", "Argv0");
        exit(1);
      }
      strcpy(&zSql[nSql++], "\n");
      strcpy(&zSql[nSql], zLine);
      nSql += len;
    }
    free(zLine);
    if( zSql && _ends_with_semicolon(zSql, nSql) && sqlite3_complete(zSql) ){
      rc = sqlite3_exec(*db, zSql, 0, 0, &zErrMsg);
      if( rc || zErrMsg ){
        if( in!=0 ) printf("%s\n",zSql);
        if( zErrMsg!=0 ){
          printf("SQL error: %s\n", zErrMsg);
          fprintf(stderr, "SQL error: statement before or near line %02d %s\n",count, zErrMsg);
         strncpy( errmsg, zErrMsg , len_errmsg ) ;
          sqlite3_free(zErrMsg);
          zErrMsg = 0;
          return rc;
        }else{
          fprintf(stderr,"SQL error: %s\n", sqlite3_errmsg(*db));
         strncpy( errmsg, sqlite3_errmsg(*db) , len_errmsg ) ;
          return rc;
        }
      }
      free(zSql);
      zSql = 0;
      nSql = 0;
    }
  }
  if( zSql ){
    if( !_all_whitespace(zSql) ) printf("Incomplete SQL: %s\n", zSql);
    free(zSql);
  }
   fclose(in);
   return rc;
}
long long int sqlite3_last_insert_rowid_c_(sqlite3* * db) { 

   return sqlite3_last_insert_rowid(*db) ;
} 

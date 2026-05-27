#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include <zlib.h>
# include <errno.h>
# include <signal.h>
# include <pwd.h>
# include <unistd.h>
# include <sys/types.h>
#include <rmn.h>
# include "udfsqlite.h"
#include "sqlite3ext.h"
SQLITE_EXTENSION_INIT1


/*
 * Declaration des fonctions (partiel)
 */
void check_latlon(float lat, float lon);
void concat_result(char ** buf, float result1, float result2);
void udf_get_fstvar(char *filename, char *nomvar, char *typvar_arg,
                    char *etiket_arg, char *ez_opt, int ip1_arg, 
                    int ip2_arg, int ip3_arg, int *gdid, float ** field );
static void udf_gdllsval  ( sqlite3_context *context, int argc, sqlite3_value **argv );
static void udf_gdxysval  ( sqlite3_context *context, int argc, sqlite3_value **argv );
static void udf_gdllvval  ( sqlite3_context *context, int argc, sqlite3_value **argv );
static void udf_gdxyvval  ( sqlite3_context *context, int argc, sqlite3_value **argv );
static void udf_gdllwdval ( sqlite3_context *context, int argc, sqlite3_value **argv );
static void udf_gdxywdval ( sqlite3_context *context, int argc, sqlite3_value **argv );

#define AGG_ALLOC_START 150000
#define AGG_ALLOC_INCR  50000

typedef struct VarCtx VarCtx;
struct VarCtx {
  double sum;    /* Sum of terms */
  double sum2;    /* Sum of terms */
  int cnt;                /* Number of elements summed */
  int size;
  double * values;
};

# define max_fstd_fields 20
static int fstd_field_count = 0;

typedef struct fstd_field 
{ 
   char *filename; 
   char *nomvar;
   char *typvar;
   char *etiket;
   char *ez_opt;
   float *field;
   int gdid; 
   int ip1;
   int ip2;
   int ip3;
} fstd_field ;
static fstd_field fstd_fields[max_fstd_fields] = {0};

 /********************************************************************
  * s/r  regrupBox (...)                                             *
  * objet    :    sous-programme pour calculer un numero de boite    *
  * pour storer les observations regroupees de l'ade.                *
  *                                                                  *
  * arguments:    lat  (in)  - float - latitude  en degres [-90 +90] *
  * lon  (in)  - float - longitude en degres [-180 +180]             *
  * dx   (in)  - int   - largeur de la boite en dixiemes de degre    *
  * dy   (in)  - int   - hauteur de la boite en dixiemes de degre    *
  * glat (out) - int   - position inferieure gauche de               *
  * la boite en 100eme degre                                         *
  * glon (out) - int   - position inferieure gauche de               *
  * la boite en 100eme degre                                         *
  * no   (out) - int   - numero de boite                             *
  ********************************************************************/

static void regrupBox ( float lat, float lon, int dx, int dy, int *glat1, int *glon1, int *no1 ) 
{
/*    conversion des lat et lon en 10eme de degres*/
   int nblat,nblon,nolat,nolon; 
   int no  = -1;
   int glat=  0; 
   int glon=  0;
   int latitude = (int) floor(100. *lat +9000 +0.5);
   int longitude = (int)floor (100. * lon + 0.5);
   latitude /= 10; 
   if (longitude < 0)
      longitude = longitude + 36000;
   longitude /=10;

   if ( dx <= 0 || dy <= 0 ) 
      return;

/*   on verifie que dx et dy sont valides, pour cela, il faut que la*/
/*   division par le nombre de degres correspondant soit egale a 0*/

   if ( ((1800%dy)+(3600%dx)) != 0 ) 
      return;

/*    on calcule le nombre de rangees et de colonnes*/

   nblat = 1800 / dy;
   nblon = 3600 / dx;

/*    on verifie si la station est au pole nord*/

   if ( latitude == 1800 ) {
      nolat = 0;
      nolon = 0;
      no = nblat * nblon + 1;
   } else {
      nolat =  latitude   / dy;
      nolon =  longitude  / dx;
      no    = nolat * nblon + nolon + 1;
   }

/*    calculons le point inferieur gauche de la boite en centiemes de degres*/

   glat = nolat * dy * 10;
   glon = nolon * dx * 10;
   *glon1 = glon;
   *glat1 = glat;
   *no1    = no;
}
static void myRegrupBoxFunc(sqlite3_context *context, int argc, sqlite3_value **argv)
{
   int boxLatitude;
   int boxLongitude;
   int boxId;
   int dx = sqlite3_value_int(argv[2]);
   int dy = sqlite3_value_int(argv[3]);
   float latitude  = (float) sqlite3_value_double(argv[0]);
   float longitude = (float) sqlite3_value_double(argv[1]);
   regrupBox ( latitude,  longitude,  dx,  dy,  &boxLatitude,
                  &boxLongitude,  &boxId ); 
   sqlite3_result_int(context, boxId);
}
static void myRegrupLat100Func(sqlite3_context *context, int argc, sqlite3_value **argv)
{
   int boxLatitude;
   int boxLongitude;
   int boxId;
   int dx = sqlite3_value_int(argv[2]);
   int dy = sqlite3_value_int(argv[3]);
   float latitude  = (float) sqlite3_value_double(argv[0]);
   float longitude = (float) sqlite3_value_double(argv[1]);
   regrupBox ( latitude,  longitude,  dx,  dy,  &boxLatitude,
                  &boxLongitude,  &boxId ); 
   sqlite3_result_int(context, boxLatitude);
}
static void myRegrupLatFunc(sqlite3_context *context, int argc, sqlite3_value **argv)
{
   int boxLatitude;
   double lat;
   int boxLongitude;
   int boxId;
   int dx = sqlite3_value_int(argv[2]);
   int dy = sqlite3_value_int(argv[3]);
   float latitude  = (float) sqlite3_value_double(argv[0]);
   float longitude = (float) sqlite3_value_double(argv[1]);
   regrupBox ( latitude,  longitude,  dx,  dy,  &boxLatitude,
                  &boxLongitude,  &boxId ); 
   lat = (boxLatitude - 9000)/100.;
   sqlite3_result_double(context, lat);

}
void myRegrupLon100Func(sqlite3_context *context, int argc, sqlite3_value **argv)
{
   int boxLatitude;
   int boxLongitude;
   int boxId;
   int dx = sqlite3_value_int(argv[2]);
   int dy = sqlite3_value_int(argv[3]);
   float latitude  = (float) sqlite3_value_double(argv[0]);
   float longitude = (float) sqlite3_value_double(argv[1]);
   regrupBox ( latitude,  longitude,  dx,  dy,  &boxLatitude,
                  &boxLongitude,  &boxId ); 
   sqlite3_result_int(context, boxLongitude);
}
static void myRegrupLonFunc(sqlite3_context *context, int argc, sqlite3_value **argv)
{
   int boxLatitude;
   double lon;
   int boxLongitude;
   int boxId;
   int dx = sqlite3_value_int(argv[2]);
   int dy = sqlite3_value_int(argv[3]);
   float latitude  = (float) sqlite3_value_double(argv[0]);
   float longitude = (float) sqlite3_value_double(argv[1]);
   regrupBox ( latitude,  longitude,  dx,  dy,  &boxLatitude,
                  &boxLongitude,  &boxId ); 
   if (boxLongitude <= 18000)
      lon =  boxLongitude/100.0;
   else
      lon = (boxLongitude - 36000)/100.0;
   sqlite3_result_double(context, lon);
}

static void myTimeStampFunc (sqlite3_context *context, int argc, sqlite3_value **argv)
{
    int yymmdd = sqlite3_value_int(argv[0]);
    int hhmmss = sqlite3_value_int(argv[1]);
/*     calcul des donnees individuelles ( on suppose les donnes deja validees )*/
    int year   = yymmdd / 10000;
    int month  = (yymmdd/100)%100; 
    int day    =  yymmdd%100 ;
    int hhmm   = hhmmss / 100;
    int hour   = hhmm / 100;
    int minute = hhmm%100;
    int valtim;
    putenv("TZ=GMT");
/*    fprintf ( stderr, " %d/%d/%d hhmm = %d hour = %d min = %d fuseau hor = %s\n",*/
/*          year, month, day, hhmm,  hour,minute,tzname[0] ) ;*/

    if ( hhmm >= 2100 )  {
      time_t tim ;
      struct tm *temps ;
      time(&tim);
      temps = localtime(&tim);
      temps->tm_year = year-1900;
      temps->tm_mon  = month - 1 ;
      temps->tm_mday = day;
      temps->tm_hour = hour;
      temps->tm_min  = minute;
      tim = mktime(temps);
      tim = tim + 86400;
      temps  = localtime(&tim);
      year   = temps->tm_year +1900;
      month  = temps->tm_mon +1;
      day    = temps->tm_mday;
/*   fprintf ( stderr, " %2.2d:%2.2d:%2.2d %2.2d/%2.2d/%2.2d \n",*/
/*                     temps->tm_hour, temps->tm_min, temps->tm_sec, */
/*                     temps->tm_mday, temps->tm_mon +1, temps->tm_year+1900 ) ;*/
      


/*      temps.tm_min ;*/
/*      temps.tm_sec ;*/
       hour = 00;
    } else if (hhmm < 300 ) {
       hour = 00;
    } else if ( hhmm >= 300 && hhmm < 900) {
       hour = 06;
    } else if (hhmm >= 900 && hhmm < 1500) {
       hour = 12; 
    } else if (hhmm >= 1500) {
       hour = 18;
    }
    valtim = ( ( year * 100 + month ) * 100 + day ) * 100 + hour;
    sqlite3_result_int(context, valtim);

}

static char * ConcatText(char *zIn, char const *zAppend){
   int len;
/*   int i;*/
   int nAppend = strlen(zAppend);
   int nIn = (zIn?strlen(zIn):0);

   len = nAppend+nIn+1;

   zIn = (char *)realloc(zIn, len);
   if( !zIn )
      return 0;

   memcpy(&zIn[nIn], zAppend, nAppend);
   zIn[len-1] = '\0';

   return zIn;
}

/*
** Implementation of the Left functions  
** All three do the same thing.  They return the first non-NULL
** argument.
*/
static void LeftFunc( sqlite3_context *context, int argc, sqlite3_value **argv){
/*   int i;*/
   int Len;
   int pos = sqlite3_value_int(argv[1]) ;
   char *zTmp = 0;
   if( SQLITE_NULL==sqlite3_value_type(argv[0]) ){
      sqlite3_result_null(context);
      return;
   }
   Len = strlen((const char *)sqlite3_value_text(argv[0]));
   if (pos >= Len) {
      sqlite3_result_text(context,(const char*)sqlite3_value_text(argv[0]) , -1, SQLITE_TRANSIENT);
      return;
   }
   if (pos > 0) 
      zTmp =( char *)malloc(pos+1);
   else 
      sqlite3_result_text(context,(const char *)sqlite3_value_text(argv[0]),-1,SQLITE_TRANSIENT);

   if( zTmp ) {
      strncpy((char*)zTmp, (char*)sqlite3_value_text(argv[0]),pos);
      zTmp[pos] = '\0';
      sqlite3_result_text(context, zTmp, -1, SQLITE_TRANSIENT);
      free(zTmp);
   }
}
/*
** Implementation of the concat functions  
** All three do the same thing.  They return the first non-NULL
** argument.
*/
static void ConcatFunc( sqlite3_context *context, int argc, sqlite3_value **argv){
   int i;
   char *zTmp = 0;
   for(i=0; i<argc; i++){
      if( SQLITE_NULL==sqlite3_value_type(argv[i]) ){
         sqlite3_result_null(context);
         if( zTmp ) free(zTmp);
         return;
      }
      zTmp = ConcatText(zTmp,(char const *)sqlite3_value_text(argv[i]));
   }
   sqlite3_result_text(context, zTmp, -1, SQLITE_TRANSIENT);
   if( zTmp ) free(zTmp);
}
/*
** Implementation of the sign() function
*/
static void msignFunc(sqlite3_context *context, int argc, sqlite3_value **argv){
  assert( argc==1 );
  switch( sqlite3_value_type(argv[0]) ){
    case SQLITE_INTEGER: {
      long long int iVal = sqlite3_value_int64(argv[0]);
      iVal = ( iVal > 0) ? 1 : ( iVal < 0 ) ? -1 : 0;
      sqlite3_result_int64(context, iVal);
      break;
    }
    case SQLITE_NULL: {
      sqlite3_result_null(context);
      break;
    }
    default: {
      double rVal = sqlite3_value_double(argv[0]);
      rVal = ( rVal > 0) ? 1 : ( rVal < 0 ) ? -1 : 0;
      sqlite3_result_double(context, rVal);
      break;
    }
  }
}
/*
** Implementation of the m_abs() function
*/
static void mAbsFunc(sqlite3_context *context, int argc, sqlite3_value **argv){
  assert( argc==1 );
  switch( sqlite3_value_type(argv[0]) ){
    case SQLITE_INTEGER: {
      long long int iVal = sqlite3_value_int64(argv[0]);
      iVal = (long long int) (llabs(iVal));
      sqlite3_result_int64(context, iVal);
      break;
    }
    case SQLITE_NULL: {
      sqlite3_result_null(context);
      break;
    }
    default: {
      double rVal = sqlite3_value_double(argv[0]);
      rVal = fabs(rVal);
      sqlite3_result_double(context, rVal);
      break;
    }
  }
}
/*
** Implementation of the PI() function
*/
static void PIFunc(sqlite3_context *context, int argc, sqlite3_value **argv){
   sqlite3_result_double(context, M_PI);
}
/*
** Implementation of the POW() function
*/
static void PowFunc(sqlite3_context *context, int argc, sqlite3_value **argv){
  assert( argc==2 );
  switch( sqlite3_value_type(argv[0]) ){
    case SQLITE_NULL: {
      sqlite3_result_null(context);
      break;
    }
    default: {
      double rVal = sqlite3_value_double(argv[0]);
      sqlite3_result_double(context, pow(rVal, sqlite3_value_double(argv[1])));
      break;
    }
  }
}
/*
** Implementation of the SIN() function
*/
static void SinFunc(sqlite3_context *context, int argc, sqlite3_value **argv){
  assert( argc==1 );
  switch( sqlite3_value_type(argv[0]) ){
    case SQLITE_NULL: {
      sqlite3_result_null(context);
      break;
    }
    default: {
      double rVal = sqlite3_value_double(argv[0]);
      sqlite3_result_double(context, sin(rVal));
      break;
    }
  }
}
/*
** Implementation of the ASIN() function
*/
static void AsinFunc(sqlite3_context *context, int argc, sqlite3_value **argv){
  assert( argc==1 );
  switch( sqlite3_value_type(argv[0]) ){
    case SQLITE_NULL: {
      sqlite3_result_null(context);
      break;
    }
    default: {
      double rVal = sqlite3_value_double(argv[0]);
      sqlite3_result_double(context, asin(rVal));
      break;
    }
  }
}
/*
** Implementation of the cos() function
*/
static void CosFunc(sqlite3_context *context, int argc, sqlite3_value **argv){
  assert( argc==1 );
  switch( sqlite3_value_type(argv[0]) ){
    case SQLITE_NULL: {
      sqlite3_result_null(context);
      break;
    }
    default: {
      double rVal = sqlite3_value_double(argv[0]);
      sqlite3_result_double(context, cos(rVal));
      break;
    }
  }
}
/*
** Implementation of the Acos() function
*/
static void AcosFunc(sqlite3_context *context, int argc, sqlite3_value **argv){
  assert( argc==1 );
  switch( sqlite3_value_type(argv[0]) ){
    case SQLITE_NULL: {
      sqlite3_result_null(context);
      break;
    }
    default: {
      double rVal = sqlite3_value_double(argv[0]);
      sqlite3_result_double(context, acos(rVal));
      break;
    }
  }
}
/*
** Implementation of the SQRT() function
*/
static void SqrtFunc(sqlite3_context *context, int argc, sqlite3_value **argv){
  assert( argc==1 );
  switch( sqlite3_value_type(argv[0]) ){
/*    case SQLITE_INTEGER: {*/
/*      long long int iVal = sqlite3_value_int64(argv[0]);*/
/*      iVal = (long long int) (sqrt(iVal));*/
/*      sqlite3_result_int64(context, iVal);*/
/*      break;*/
/*    }*/
    case SQLITE_NULL: {
      sqlite3_result_null(context);
      break;
    }
    default: {
      double rVal = sqlite3_value_double(argv[0]);
      rVal = sqrt(rVal);
      sqlite3_result_double(context, rVal);
      break;
    }
  }
}
/*
** Implementation of the FLOOR() function
*/
static void FloorFunc(sqlite3_context *context, int argc, sqlite3_value **argv){
  assert( argc==1 );
  switch( sqlite3_value_type(argv[0]) ){
    case SQLITE_INTEGER: {
      long long int iVal = sqlite3_value_int64(argv[0]);
      iVal = (long long int) (floor(iVal));
      sqlite3_result_int64(context, iVal);
      break;
    }
    case SQLITE_NULL: {
      sqlite3_result_null(context);
      break;
    }
    default: {
      double rVal = sqlite3_value_double(argv[0]);
      long long int iVal =  (long long int) (floor(rVal));
      sqlite3_result_int64(context, iVal);
      break;
    }
  }
}
/*
** Implementation of the Ceil() function
*/
static void CeilFunc(sqlite3_context *context, int argc, sqlite3_value **argv){
  assert( argc==1 );
  switch( sqlite3_value_type(argv[0]) ){
    case SQLITE_INTEGER: {
      long long int iVal = sqlite3_value_int64(argv[0]);
      iVal = (long long int) (ceil(iVal));
      sqlite3_result_int64(context, iVal);
      break;
    }
    case SQLITE_NULL: {
      sqlite3_result_null(context);
      break;
    }
    default: {
      double rVal = sqlite3_value_double(argv[0]);
      long long int iVal =  (long long int) (ceil(rVal));
      sqlite3_result_int64(context, iVal);
      break;
    }
  }
}
static void XSTDStep(sqlite3_context *context, int argc, sqlite3_value **argv){
   VarCtx *p;
   int type;
/*   int i;*/
   double rval;
   assert( argc==1 );
   p = sqlite3_aggregate_context(context, sizeof(*p));
   type = sqlite3_value_type(argv[0]);
   if( p && type!=SQLITE_NULL ){
      if (p->cnt == 0) { 
         p->size = AGG_ALLOC_START;
         p->values = (double *) malloc( p->size * sizeof(double));
      }
   }
   if( p && type!=SQLITE_NULL ){
      p->cnt++;
      if (p->cnt > p->size) {
/*         fprintf(stderr,"cnt = %d   p.size = %d\n", p->cnt, p->size);*/
         p->size += AGG_ALLOC_INCR; 
         p->values = realloc(p->values, (p->size) *sizeof(double));
      }
      rval     = sqlite3_value_double(argv[0]);
      (p->values)[p->cnt - 1] = rval;
      p->sum += rval;
   }
}
static void XSTDFinalize(sqlite3_context *context){
   double ave;
/*   double adev;*/
/*   double sdev;*/
   double var; 
/*   double skew; */
/*   double curt;*/
   VarCtx *p;
   int i;
   p = sqlite3_aggregate_context(context, 0);
   if( !p ) return;
   if( p && p->cnt>1 ){
      ave = p->sum/p->cnt;
      p->sum2 = 0;
      for (i = 0; i != p->cnt; ++i) {
         p->sum2 += pow((p->values)[i] - ave,2);
      }
      sqlite3_result_double(context, sqrt(p->sum2/p->cnt));
   }
   if( p && p->cnt <=1 ){
      var = 0.0;
      sqlite3_result_double(context, var);
   }

   free(p->values);
   p->values = 0;

}
static void XVARFinalize(sqlite3_context *context){
   double ave;
/*   double adev;*/
/*   double sdev;*/
   double var; 
/*   double skew; */
/*   double curt;*/
   VarCtx *p;
   int i;
   p = sqlite3_aggregate_context(context, 0);
   if( !p ) return;
   if( p && p->cnt>1 ){
      ave = p->sum/p->cnt;
      p->sum2 = 0;
      for (i = 0; i != p->cnt; ++i) {
         p->sum2 += pow((p->values)[i] - ave,2);
      }
      sqlite3_result_double(context, p->sum2/p->cnt);
   }
   if( p && p->cnt <=1 ){
      var = 0.0;
      sqlite3_result_double(context, var);
   }

   free(p->values);
   p->values = 0;

}
static void VarianceStep(sqlite3_context *context, int argc, sqlite3_value **argv){
   VarCtx *p;
   int type;
/*   int i;*/
   double rval;
   assert( argc==1 );
   p = sqlite3_aggregate_context(context, sizeof(*p));
   type = sqlite3_value_type(argv[0]);
   if( p && type!=SQLITE_NULL ){
      if (p->cnt == 0) { 
         p->size = AGG_ALLOC_START;
         p->values = (double *) malloc( p->size * sizeof(double));
      }
   }
   if( p && type!=SQLITE_NULL ){
      p->cnt++;
      if (p->cnt > p->size) {
         p->size += AGG_ALLOC_INCR; 
         p->values = realloc(p->values, (p->size) *sizeof(double));
      }
      rval     = sqlite3_value_double(argv[0]);
      (p->values)[p->cnt - 1] = rval;
      p->sum += rval;
   }
}
void moment(double *data, int n, double *ave, double *adev, double *sdev, double *var, double *skew, double *curt) {
/********************************************************************************************
 * Given an array of data[1..n], this routine returns its mean ave, average deviation adev, *
 * standard deviation sdev, variance var, skewness skew, and kurtosis curt.                 *
 ********************************************************************************************/
   int j;
   double ep=0.0,s,p;

   s=0.0;    /*First pass to get the mean. */

/*   for (j=0;j<n;j++) s += data[j];*/

/*   *ave=s/n;*/

   *adev=(*var)=0.0;    /* Second pass to get the first (absolute), sec- */

   for (j=0;j<n;j++) {                 /* ond, third, and fourth moments of the */
      *adev += fabs(s=data[j]-(*ave));  /* deviation from the mean. */
      ep += s;
      *var += (p=s*s);
   }
   *adev /= n;
   *var=(*var-ep*ep/n)/(n-1);    /* Corrected two-pass formula. */
   *sdev=sqrt(*var);             /* Put the pieces together according to the con- */
}
static void VarianceFinalize(sqlite3_context *context){
   double ave;
   double adev;
   double sdev;
   double var; 
   double skew; 
   double curt;
   VarCtx *p;
/*   int i;*/
   p = sqlite3_aggregate_context(context, 0);
   if( !p ) return;
   if( p && p->cnt>1 ){
      ave = p->sum/p->cnt;
      moment(p->values, p->cnt, &ave, &adev, &sdev, &var, &skew, &curt);
      sqlite3_result_double(context, var);
   }
   if( p && p->cnt <=1 ){
      var = 0.0;
      sqlite3_result_double(context, var);
   }

   free(p->values);
   p->values = 0;

}
static void StddevFinalize(sqlite3_context *context){
   double ave;
   double adev;
   double sdev;
   double var; 
   double skew; 
   double curt;
   VarCtx *p;
/*   int i;*/
   p = sqlite3_aggregate_context(context, 0);
   if( !p ) return;
   if( p && p->cnt>1 ){
      ave = p->sum/p->cnt;
      moment(p->values, p->cnt, &ave, &adev, &sdev, &var, &skew, &curt);
      sqlite3_result_double(context, sdev);
   }
   if( p && p->cnt <=1 ){
      sdev = 0.0;
      sqlite3_result_double(context, sdev);
   }

   free(p->values);
   p->values = 0;


}
/*static bool make_datetime(date_time_format_types format, TIME *ltime,*/
/*                           String *str)*/
/* {*/
/*   char *buff;*/
/*   CHARSET_INFO *cs= &my_charset_bin;*/
/*   uint length= 30;*/
/* */
/*   if (str->alloc(length))*/
/*     return 1;*/
/*   buff= (char*) str->ptr();*/
/* */
/*   switch (format) {*/
/*   case TIME_ONLY:*/
/*     length= cs->cset->snprintf(cs, buff, length, "%s%02d:%02d:%02d",*/
/*                                ltime->neg ? "-" : "",*/
/*                                ltime->hour, ltime->minute, ltime->second);*/
/*     break;*/
/*   case TIME_MICROSECOND:*/
/*     length= cs->cset->snprintf(cs, buff, length, "%s%02d:%02d:%02d.%06d",*/
/*                                ltime->neg ? "-" : "",*/
/*                                ltime->hour, ltime->minute, ltime->second,*/
/*                                ltime->second_part);*/
/*     break;*/
/*   case DATE_ONLY:*/
/*     length= cs->cset->snprintf(cs, buff, length, "%04d-%02d-%02d",*/
/*                                ltime->year, ltime->month, ltime->day);*/
/*     break;*/
/*   case DATE_TIME:*/
/*     length= cs->cset->snprintf(cs, buff, length,*/
/*                                "%04d-%02d-%02d %02d:%02d:%02d",*/
/*                                ltime->year, ltime->month, ltime->day,*/
/*                                ltime->hour, ltime->minute, ltime->second);*/
/*     break;*/
/*   case DATE_TIME_MICROSECOND:*/
/*     length= cs->cset->snprintf(cs, buff, length,*/
/*                                "%04d-%02d-%02d %02d:%02d:%02d.%06d",*/
/*                                ltime->year, ltime->month, ltime->day,*/
/*                                ltime->hour, ltime->minute, ltime->second,*/
/*                                ltime->second_part);*/
/*     break;*/
/*   }*/
/* */
/*   str->length(length);*/
/*   str->set_charset(cs);*/
/*   return 0;*/
/* }*/

/*
** Implementation of the datetime function  
*/
static void DateTimeFunc( int yyyymmdd, int hhmmss, int * year, int * month, int * day,
                          int * hour, int * minute, int * second)
{
   int mmdd, mmss;

   *year  = yyyymmdd / 10000;
   mmdd  = yyyymmdd % 10000;
   *month = mmdd / 100;
   *day   = mmdd % 100;
   *hour  = hhmmss / 10000;
   mmss  = hhmmss % 10000;
   *minute= mmss / 100;
   *second= mmss % 100;

}
/*
** implementation of the isodatetime function  
*/
static void isoDateTimeFunc( sqlite3_context *context, int argc, sqlite3_value **argv){
   int length = 30;
   char buff[length];
   int year   ;
   int month  ;
   int day    ;
   int hour   ;
   int minute ;
   int second ;
   int yyyymmdd ;
   int yyyymmddhh;
   int hhmmss   ;
   if (argc ==2) {
   if( (SQLITE_NULL==sqlite3_value_type(argv[0])) || (SQLITE_NULL == sqlite3_value_type(argv[1])) ){
      sqlite3_result_null(context);
      return;
   }
      yyyymmdd   = sqlite3_value_int(argv[0]) ;
      hhmmss     = sqlite3_value_int(argv[1]) ;
   } else if (argc == 1) {
   if( (SQLITE_NULL==sqlite3_value_type(argv[0]))  ){
      sqlite3_result_null(context);
      return;
   }
      yyyymmddhh   = sqlite3_value_int(argv[0]) ;
      yyyymmdd     = yyyymmddhh/100;
      hhmmss       = yyyymmddhh % 100;
      hhmmss       *=10000;
   } else {
      sqlite3_result_null(context);
      return;
   }

   DateTimeFunc(yyyymmdd, hhmmss, &year,&month,&day, &hour, &minute, &second);

   length= snprintf( buff, length,
                              "%04d-%02d-%02d %02d:%02d:%02d",
                              year, month, day,
                              hour, minute, second);
   if( length > 0 ) {
      sqlite3_result_text(context, buff, length, SQLITE_TRANSIENT);
   }
}
/*
** implementation of the isodate only function  
*/
static void isoDateFunc( sqlite3_context *context, int argc, sqlite3_value **argv){
   int length = 30;
   char buff[length];
   int year   ;
   int month  ;
   int day    ;
   int hour   ;
   int minute ;
   int second ;
   int yyyymmdd ;
   int hhmmss   ;
   if (argc !=1) {
      sqlite3_result_null(context);
      return;
   }
   yyyymmdd   = sqlite3_value_int(argv[0]) ;
   hhmmss     = 0;
   if( (SQLITE_NULL==sqlite3_value_type(argv[0]))  ){
      sqlite3_result_null(context);
      return;
   }

   DateTimeFunc(yyyymmdd, hhmmss, &year,&month,&day, &hour, &minute, &second);

   length= snprintf( buff, length, "%04d-%02d-%02d", year, month, day);
   if( length > 0 ) {
      sqlite3_result_text(context, buff, length, SQLITE_TRANSIENT);
   }
}
/*
** implementation of the isotime only function  
*/
static void isoTimeFunc( sqlite3_context *context, int argc, sqlite3_value **argv){
   int length = 30;
   char buff[length];
   int year   ;
   int month  ;
   int day    ;
   int hour   ;
   int minute ;
   int second ;
   int yyyymmdd ;
   int hhmmss   ;
   if (argc !=1) {
      sqlite3_result_null(context);
      return;
   }
   yyyymmdd   = 0;
   hhmmss     = sqlite3_value_int(argv[0]) ;
   if( (SQLITE_NULL==sqlite3_value_type(argv[0]))  ){
      sqlite3_result_null(context);
      return;
   }

   DateTimeFunc(yyyymmdd, hhmmss, &year,&month,&day, &hour, &minute, &second);
   length= snprintf( buff, length, "%02d:%02d:%02d", hour, minute, second);

   if( length > 0 ) {
      sqlite3_result_text(context, buff, length, SQLITE_TRANSIENT);
   }
}

/*
** Implementation of the datetime function  
static void IsoDateTimeFunc( sqlite3_context *context, int argc, sqlite3_value **argv){
   int length = 30;
   char buff[length];
   char *zTmp = 0;
   int year   = 2007;
   int month  = 8;
   int day    = 9;
   int hour   = 12;
   int minute = 3;
   int second = 12;
   int yyyymmdd ;
   int hhmmss   ;
   int mmdd, mmss;
   if (argc !=2) {
      sqlite3_result_null(context);
      return;
   }
   yyyymmdd   = sqlite3_value_int(argv[0]) ;
   hhmmss     = sqlite3_value_int(argv[1]) ;
   if( (SQLITE_NULL==sqlite3_value_type(argv[0])) || (SQLITE_NULL == sqlite3_value_type(argv[1])) ){
      sqlite3_result_null(context);
      return;
   }

   year  = yyyymmdd / 10000;
   mmdd  = yyyymmdd % 10000;
   month = mmdd / 100;
   day   = mmdd % 100;
   hour  = hhmmss / 10000;
   mmss  = hhmmss % 10000;
   minute= mmss / 100;
   second= mmss % 100;


   length= snprintf( buff, length,
                              "%04d-%02d-%02d %02d:%02d:%02d",
                              year, month, day,
                              hour, minute, second);
   if( length > 0 ) {
      sqlite3_result_text(context, buff, length, SQLITE_TRANSIENT);
   }
}
*/

/*
** Implementation of the distance function  
**  ref (acurate formula for all disatnces) the same is implemented in IDL
*/
static double distance_rds( double lat1, double lon1, double lat2, double lon2){
   double s;
/*   double dlat = lat2-lat1;*/
   double dlon = lon2-lon1;
   double num;
   double den;

/*   s = 2*asin(sqrt(sin( (lat1-lat2)/2)*sin((lat1-lat2)/2) + cos(lon1)*cos(lon2)*sin((lon1-lon2)/2)*sin((lon1-lon2)/2)));*/
   
   num = pow(cos(lat2)*sin(dlon),2) + pow(cos(lat1) *sin(lat2) - sin(lat1)* cos(lat2) * cos(dlon),2);
   num = sqrt(num);
   den = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)* cos(dlon);
   s = atan2(num,den);

   return(s);

}
/*
** Implementation of the distance function  
**  ref (haversine formula)
*/
static double distance_rds2( double lat1, double lon1, double lat2, double lon2){
   double s;
/*   double dlat = lat2-lat1;*/
/*   double dlon = lon2-lon1;*/
/*   double num;*/
/*   double den;*/

   s = 2*asin(sqrt(sin( (lat1-lat2)/2)*sin((lat1-lat2)/2) + cos(lon1)*cos(lon2)*sin((lon1-lon2)/2)*sin((lon1-lon2)/2)));
   
   return(s);

}
/*
** Implementation of the distance  kms function  
*/
static void KmsFunc( sqlite3_context *context, int argc, sqlite3_value **argv){
   double lat1,lon1,lat2,lon2;
   double s;
/*   double PI = 3.14159265359;*/
/*   double TWOPI = 6.28318530718;*/
   double DE2RA = 0.01745329252;
/*   double RA2DE = 57.2957795129;*/
/*   double ERAD = 6378.135;*/
   double R_EARTH = 6378.2064 ; /*Earth equatorial radius, meters, Clarke 1866 ellipsoid */
/*   double ERADM = 6378135.0;*/
/*   double AVG_ERAD = 6371.0;*/
/*   double EPS = 0.000000000005;*/
/*   double KM2MI = 0.621371;*/
/*   double FLATTENING =  0;*/
   if (argc !=4) {
      sqlite3_result_null(context);
      return;
   }
   lat1   = sqlite3_value_double(argv[0]) ;
   lon1   = sqlite3_value_double(argv[1]) ;
   lat2   = sqlite3_value_double(argv[2]) ;
   lon2   = sqlite3_value_double(argv[3]) ;
   if( (SQLITE_NULL==sqlite3_value_type(argv[0])) || (SQLITE_NULL == sqlite3_value_type(argv[1])) ){
      sqlite3_result_null(context);
      return;
   }
   lat1 *= DE2RA;
   lon1 *= DE2RA;
   lat2 *= DE2RA;
   lon2 *= DE2RA;


   s = distance_rds(lat1,lon1,lat2,lon2);
   sqlite3_result_double(context, R_EARTH*s) ;

}
/*
** Implementation of the distance  miles function  
*/
static void MilesFunc( sqlite3_context *context, int argc, sqlite3_value **argv){
   double lat1,lon1,lat2,lon2;
   double s;
   double DE2RA = 0.01745329252;
   double R_EARTH = 6378.2064 ; /*Earth equatorial radius, meters, Clarke 1866 ellipsoid */
   double KM2MI = 0.621371;
   if (argc !=4) {
      sqlite3_result_null(context);
      return;
   }
   lat1   = sqlite3_value_double(argv[0]) ;
   lon1   = sqlite3_value_double(argv[1]) ;
   lat2   = sqlite3_value_double(argv[2]) ;
   lon2   = sqlite3_value_double(argv[3]) ;
   if( (SQLITE_NULL==sqlite3_value_type(argv[0])) || (SQLITE_NULL == sqlite3_value_type(argv[1])) ){
      sqlite3_result_null(context);
      return;
   }
   lat1 *= DE2RA;
   lon1 *= DE2RA;
   lat2 *= DE2RA;
   lon2 *= DE2RA;


   s = distance_rds(lat1,lon1,lat2,lon2);
   sqlite3_result_double(context, KM2MI*R_EARTH*s) ;

}
/*
** Implementation of the distance  miles function  
*/
static void Miles2Func( sqlite3_context *context, int argc, sqlite3_value **argv){
   double lat1,lon1,lat2,lon2;
   double s;
   double DE2RA = 0.01745329252;
   double R_EARTH = 6378.2064 ; /*Earth equatorial radius, meters, Clarke 1866 ellipsoid */
   double KM2MI = 0.621371;
   if (argc !=4) {
      sqlite3_result_null(context);
      return;
   }
   lat1   = sqlite3_value_double(argv[0]) ;
   lon1   = sqlite3_value_double(argv[1]) ;
   lat2   = sqlite3_value_double(argv[2]) ;
   lon2   = sqlite3_value_double(argv[3]) ;
   if( (SQLITE_NULL==sqlite3_value_type(argv[0])) || (SQLITE_NULL == sqlite3_value_type(argv[1])) ){
      sqlite3_result_null(context);
      return;
   }
   lat1 *= DE2RA;
   lon1 *= DE2RA;
   lat2 *= DE2RA;
   lon2 *= DE2RA;


   s = distance_rds2(lat1,lon1,lat2,lon2);
   sqlite3_result_double(context, KM2MI*R_EARTH*s) ;

}

static void getIntFromBlob ( sqlite3_context *context, int argc, sqlite3_value **argv ) {
    int index, size;
    
    if ( sqlite3_value_type(argv[0]) != SQLITE_BLOB ) {
        sqlite3_result_error(context, "Column type is not BLOB",-1);
        return;
    }
    
    size = sqlite3_value_bytes(argv[0]);
    index = sqlite3_value_int(argv[1]);
    
    if ( index > size ) {
        sqlite3_result_error(context, "Out of bounds",-1);
        return;
    }
    
    unsigned char *buf;
    buf = (unsigned char*) sqlite3_value_blob(argv[0]);

    sqlite3_result_int(context, (int) buf[index-1]);
}

static void getDoubleFromBlob ( sqlite3_context *context, int argc, sqlite3_value **argv ) {
    int index, size;
    
    if ( sqlite3_value_type(argv[0]) != SQLITE_BLOB ) {
        sqlite3_result_error(context, "Column type is not BLOB",-1);
        return;
    }
    
    size = sqlite3_value_bytes(argv[0]);
    index = sqlite3_value_int(argv[1]);
    
    if ( index > size ) {
        sqlite3_result_error(context, "Out of bounds",-1);
        return;
    }
    
    double *buf;
    buf = (double*) sqlite3_value_blob(argv[0]);

    sqlite3_result_double(context, buf[index-1]);
}

static void getCountChar ( sqlite3_context *context, int argc, sqlite3_value **argv ) {
    int count = 0;
    
    if ( sqlite3_value_type(argv[0]) != SQLITE_TEXT ) {
        sqlite3_result_error(context, "Column type is not a string",-1);
        return;
    }
    
    char *str = (char*) sqlite3_value_text(argv[0]);
    char *character = (char*) sqlite3_value_text(argv[1]);
    
    int i;
    for (i=0; i<strlen(str);i++) {
        if (character[0] == str[i]) 
            count++;
    }
    sqlite3_result_int(context, count);
}

static void mask_mer ( sqlite3_context *context, int argc, sqlite3_value **argv ) {
    static FILE *file_termer = NULL;

    char *termer_filename = getenv("DSQLITE_MASK_TERREMER");
    
    int result = -1;
    long addr;
    int  err, j, j2;
    char chaine[500];
    int  istat;
    int blat, blon;
    
    blat = round(sqlite3_value_double(argv[0]) *100) + 9000;
    blon = round(sqlite3_value_double(argv[1]) *100);
    if ( blon < 0 ) {
        blon += 36000;
    }
    
    if ( file_termer == NULL) {
        file_termer = fopen (termer_filename, "rb");
        if (file_termer == NULL) {
            strcpy(chaine, "peux pas ouvrir le fichier terre-mer");
            printf("%s : %s\n", chaine, termer_filename);
            return;
        }
    }
   
    addr = (blon + ((18000 - blat) * 36000));
    fseek(file_termer, addr/8-(addr/8)%4, SEEK_SET);
    istat = fread(&j, 4, 1, file_termer);
    err = addr % 32;

    j2 = ((j & 255) << 24) | ((j>>8 & 255) << 16) | ((j>>16 & 255) << 8) | (j>>24 & 255);
    result = 1 - ( ((j2)>>(31-err)) & 1);

    sqlite3_result_int(context,result);
}

static void crc32_func(sqlite3_context *context, int argc, sqlite3_value **argv){
    const char* argument;
    argument = (const char*)sqlite3_value_text(argv[0]);
    int crc = crc32 (0L, (unsigned char *)argument, strlen(argument));
    sqlite3_result_int(context, crc);
}

static void mask_glace ( sqlite3_context *context, int argc, sqlite3_value **argv ) {
    
    int gdid = -1;
    float *field = NULL;
    float result = 0;
    
    if (argc != 2)
    {
        fprintf(stderr, "Usage : select mask_glace lat lon");
        exit(1);
    }
    
    float lat           = sqlite3_value_double(argv[0]);
    float lon           = sqlite3_value_double(argv[1]);
    char *filename      = getenv("DSQLITE_MASK_GLACE_FSTD");
    char etiket[12+1]   = "";
    int ip1             = -1;
    int ip2             = -1;
    int ip3             = -1;
    char typvar[2+1]    = "";
    char nomvar[4+1]    = "LG";
    char ez_opt[]       = "NEAREST";    
    
    check_latlon(lat, lon);
    udf_get_fstvar(filename, nomvar, typvar, etiket, ez_opt, ip1, ip2, ip3, &gdid, &field);
    
    //Do scalar interpolation at the latlon points
    c_ezsetopt("INTERP_DEGREE", ez_opt);
    if (c_gdllsval(gdid, &result, field, &lat, &lon, 1) < 0)
    {
        fprintf(stderr,"GDLLSVAL ERROR \n");
        exit(1);
    }
    
    sqlite3_result_double(context,result);
}


static void mask_glace_clim ( sqlite3_context *context, int argc, sqlite3_value **argv ) {

    int result = -1;
    long addr;
    int blat, blon;

    FILE *file = NULL;
    static char *buffer = NULL;
    unsigned long fileLen;
    char *maskglace_filename = getenv("DSQLITE_MASK_GLACE");
    char msg[255];

    double lat = sqlite3_value_double(argv[0]);
    double lon = sqlite3_value_double(argv[1]);

    // initialization
    if ( buffer == NULL) {
        //Open file
        file = fopen(maskglace_filename, "rb");
        if (!file)
        {
            sqlite3_result_error_code(context, SQLITE_IOERR);
            return;
        }
        //Get file length
        fseek(file, 0, SEEK_END);
        fileLen = ftell(file);
        fseek(file, 0, SEEK_SET);

        //Allocate memory
        buffer=(char *)malloc(fileLen + 1);
        if (!buffer)
        {
            fclose(file);
            sqlite3_result_error_code(context, SQLITE_NOMEM);
            return;
        }
        //Read file contents into buffer
        fread(buffer, fileLen, 1, file);
        fclose(file);
    }

    // check latitude and longitude
    if (!(lat >= -90 && lat <= 90) || !(lon >= -180 && lon <= 180) ){
        snprintf(msg, sizeof(msg), "La position (%f,%f) est hors limite. ([-90 , 90], [-180, 180])\n", lat, lon);
        sqlite3_result_error(context, msg, -1);
        return;
    }

    blat = round(lat *100) + 9000;
    blon = round(lon *100);
    if ( blon < 0 ) {
        blon += 36000;
    }

    addr = ( (blon/10) + ( (1800 - (blat/10)) * 3600 ) ) ;
    result = buffer[addr];

    sqlite3_result_int(context,result);
}


void check_latlon(float lat, float lon)
{
    char msg[255];
    
    if (!(lat >= -90 && lat <= 90) || !(lon >= -180 && lon <= 180) ){
        snprintf(msg, sizeof(msg), "La position (%f,%f) est hors limite. ([-90 , 90], [-180, 180])\n", lat, lon);
        fprintf(stderr, "Error: %s \n", msg);
        exit(1);
    }
}

void udf_get_fstvar(char *filename, char *nomvar, char *typvar_arg, 
                    char *etiket_arg, char *ez_opt, int ip1_arg, 
                    int ip2_arg, int ip3_arg, int *gdid, float ** field_ptr ){

    int indice = 0;
    int exist = 0;

    if (fstd_field_count > 0 && fstd_field_count < max_fstd_fields)
    {
        while (indice < fstd_field_count && !(strcmp(fstd_fields[indice].filename, filename) == 0 
                                   && strcmp(fstd_fields[indice].nomvar, nomvar) == 0
                                   && strcmp(fstd_fields[indice].typvar, typvar_arg) == 0 
                                   && strcmp(fstd_fields[indice].etiket, etiket_arg) == 0
                                   && strcmp(fstd_fields[indice].ez_opt, ez_opt) == 0
                                   && fstd_fields[indice].ip1 == ip1_arg 
                                   && fstd_fields[indice].ip2 == ip2_arg
                                   && fstd_fields[indice].ip3 == ip3_arg)){
            indice++;
        }
        exist = indice != fstd_field_count;
    }else if (fstd_field_count > max_fstd_fields - 1)
    {
        fprintf(stderr, "Error: Too many files \n");
        fprintf(stderr, "Contact CMDS \n");
        exit(1);
    }
   
    if (exist)
    {
        *field_ptr = fstd_fields[indice].field;
        *gdid = fstd_fields[indice].gdid;
    }else
    {
        char grtyp[3] = {0};
        int ni = -1, nj = -1, nk = -1, ier = -1, key = 0, iun = fstd_field_count + 10,
            ip1 = ip1_arg, ip2 = ip2_arg, ip3 = ip3_arg;
        
        int datev = -1, dateo = -1, deet = -1, npas= -1, nbits = -1, datyp = -1,
            ig1= -1, ig2 = -1, ig3 = -1, ig4 = -1, swa = -1, lng = -1, dltf = -1,
            ubc = -1, extra1= -1, extra2 = -1, extra3 = -1;
        
        char etiket[12+1];
        char typvar[2+1];
        strcpy(etiket, etiket_arg);
        strcpy(typvar, typvar_arg);
       
        //  associate unit number to the file
        if (c_fnom((int *)iun, filename,"STD+RND+OLD", 0) !=0)
        {
            fprintf(stderr,"FNOM ERROR\n");
            exit(1);
        }

        if (c_fstouv(iun,"RND") < 0)
        {
            fprintf(stderr,"FSTOUV ERROR \n");
            exit(1);
        }
 
        //Locates the next record that matches the search keys (datev,etiket,ip1, ip2,ip3,typvar,nomvar)
        key = c_fstinf( iun, &ni, &nj, &nk,datev,etiket,ip1,ip2,ip3,typvar,nomvar );
        if (key < 0)
        {
            fprintf(stderr,"FSTINF ERROR \n");
            exit(1);
        }
        
        // Get all the description information of the record given the handle
        ier = c_fstprm(key, &dateo, &deet, &npas, &ni, &nj, &nk, &nbits, 
                        &datyp, &ip1, &ip2, &ip3, typvar, nomvar, etiket,
                        grtyp, &ig1, &ig2, &ig3, &ig4, &swa, &lng, &dltf,
                        &ubc, &extra1, &extra2, &extra3);
        if (ier < 0)
        {
            fprintf(stderr,"FSTPRM ERROR \n");
            exit(1);
        }
  
        //Read the record at position given by handle
        *field_ptr = malloc(sizeof(float) * ni * nj * nk); 
        if ( *field_ptr == NULL ) 
        {
            fprintf(stderr,"FIELD ALLOCATION ERROR \n");
            return;
        }

        if (c_fstluk(*field_ptr, key, &ni, &nj, &nk) < 0)
        {
            fprintf(stderr,"FSTLUK ERROR \n");
            exit(1);
        }

        // Define grid
        *gdid = c_ezqkdef(ni, nj, grtyp, ig1, ig2, ig3, ig4, iun);
        if (*gdid < 0)
        {
            fprintf(stderr,"EZQKDEF ERROR \n");
            exit(1);
        }
        // add to files list
        fstd_fields[fstd_field_count].filename = strdup(filename);
        fstd_fields[fstd_field_count].nomvar = strdup(nomvar);
        fstd_fields[fstd_field_count].typvar = strdup(typvar_arg);
        fstd_fields[fstd_field_count].etiket = strdup(etiket_arg);
        fstd_fields[fstd_field_count].ez_opt = strdup(ez_opt);
        fstd_fields[fstd_field_count].field = *field_ptr;
        fstd_fields[fstd_field_count].gdid = *gdid;
        fstd_fields[fstd_field_count].ip1 = ip1_arg;
        fstd_fields[fstd_field_count].ip2 = ip2_arg;
        fstd_fields[fstd_field_count].ip3 = ip3_arg;
        
        fstd_field_count++;
    }
}

static void udf_gdllsval ( sqlite3_context *context, int argc, sqlite3_value **argv ) {
    
    if (argc != 10)
    {
        fprintf(stderr, "Error: wrong number of arguments \n");
        fprintf(stderr, "Usage : select gdllsval lat lon filename etiket ip1 ip2 ip3 typevar nomvar ez_opt");
        exit(1);
    }
    
    float lat           = sqlite3_value_double(argv[0]);
    float lon           = sqlite3_value_double(argv[1]);
    char *filename      = (char *)sqlite3_value_text(argv[2]);
    char *etiket        = (char *)sqlite3_value_text(argv[3]);
    int ip1             = sqlite3_value_int(argv[4]);
    int ip2             = sqlite3_value_int(argv[5]);
    int ip3             = sqlite3_value_int(argv[6]);
    char *typvar        = (char *)sqlite3_value_text(argv[7]);
    char *nomvar        = (char *)sqlite3_value_text(argv[8]);
    char *ez_opt        = (char *)sqlite3_value_text(argv[9]);
    
    int gdid = -1;
    float *field = NULL;
    float result[1] ; 
    
    check_latlon(lat, lon);
    udf_get_fstvar(filename, nomvar, typvar, etiket, ez_opt, ip1, ip2, ip3, &gdid, &field);
             
    //Do scalar interpolation at the latlon points
    c_ezsetopt("INTERP_DEGREE", ez_opt);
    if (c_gdllsval(gdid, result, field, &lat, &lon, 1) < 0)
    {
        fprintf(stderr,"GDLLSVAL ERROR \n");
        exit(1);
    }
    
    sqlite3_result_double(context,result[0]);
}

static void udf_gdxysval ( sqlite3_context *context, int argc, sqlite3_value **argv ) {
    
    if (argc != 10)
    {
        fprintf(stderr, "Error: wrong number of arguments \n");
        fprintf(stderr, "Usage : select gdxysval lat lon filename etiket ip1 ip2 ip3 typevar nomvar ez_opt");
        exit(1);
    }
    
    float x             = sqlite3_value_double(argv[0]);
    float y             = sqlite3_value_double(argv[1]);
    char *filename      = (char *)sqlite3_value_text(argv[2]);
    char *etiket        = (char *)sqlite3_value_text(argv[3]);
    int ip1             = sqlite3_value_int(argv[4]);
    int ip2             = sqlite3_value_int(argv[5]);
    int ip3             = sqlite3_value_int(argv[6]);
    char *typvar        = (char *)sqlite3_value_text(argv[7]);
    char *nomvar        = (char *)sqlite3_value_text(argv[8]);
    char *ez_opt        = (char *)sqlite3_value_text(argv[9]);
    
    int gdid = -1;
    float *field = NULL;
    float result[1] ; 
    
    udf_get_fstvar(filename, nomvar, typvar, etiket, ez_opt, ip1, ip2, ip3, &gdid, &field);          

    //Do scalar interpolation at the latlon points
    c_ezsetopt("INTERP_DEGREE", ez_opt);
    if (c_gdxysval(gdid, result, field, &x, &y, 1) < 0)
    {
        fprintf(stderr,"GDLLSVAL ERROR");
        exit(1);
    }

    sqlite3_result_double(context,result[0]);
}

void concat_result(char ** buf, float result1, float result2)
{
    char buffer[50] = {0};
    
    snprintf( buffer, 50, "(%0.3f, %0.3f)", result1, result2);
    
    *buf = strdup(buffer);
}

static void udf_gdllvval ( sqlite3_context *context, int argc, sqlite3_value **argv ) {
    
    if (argc != 9)
    {
        fprintf(stderr, "Error: wrong number of arguments \n");
        fprintf(stderr, "Usage : select gdllvval lat lon filename etiket ip1 ip2 ip3 typevar nomvar ez_opt");
        exit(1);
    }
    
    float lat           = sqlite3_value_double(argv[0]);
    float lon           = sqlite3_value_double(argv[1]);
    char *filename      = (char *)sqlite3_value_text(argv[2]);
    char *etiket        = (char *)sqlite3_value_text(argv[3]);
    int ip1             = sqlite3_value_int(argv[4]);
    int ip2             = sqlite3_value_int(argv[5]);
    int ip3             = sqlite3_value_int(argv[6]);
    char *typvar        = (char *)sqlite3_value_text(argv[7]);
    char *ez_opt        = (char *)sqlite3_value_text(argv[8]);
    
    int gdid = -1;
    float *fieldUU = NULL;
    float *fieldVV = NULL;
    char *nomvarUU  = strdup("UU");
    char *nomvarVV  = strdup("VV");
    float resultUU[1] ;
    float resultVV[1] ;
    
    check_latlon(lat, lon);
    udf_get_fstvar(filename, nomvarUU, typvar, etiket, ez_opt, ip1, ip2, ip3, &gdid, &fieldUU);
    udf_get_fstvar(filename, nomvarVV, typvar, etiket, ez_opt, ip1, ip2, ip3, &gdid, &fieldVV);

    //Do scalar interpolation at the latlon points
    c_ezsetopt("INTERP_DEGREE", ez_opt);
    if (c_gdllvval(gdid, resultUU, resultVV, fieldUU, fieldVV, &lat, &lon, 2) < 0)
    {
        fprintf(stderr,"GDLLSVAL ERROR \n");
        exit(1);
    }
    
    char *buf = NULL;
    concat_result(&buf, resultUU[0], resultVV[0]);

    sqlite3_result_text(context, buf, -1, SQLITE_TRANSIENT);
}


static void udf_gdxyvval ( sqlite3_context *context, int argc, sqlite3_value **argv ) {
    
    if (argc != 9)
    {
        fprintf(stderr, "Error: wrong number of arguments \n");
        fprintf(stderr, "Usage : select gdxyvval lat lon filename etiket ip1 ip2 ip3 typevar ez_opt \n");
        exit(1);
    }
    
    float x             = sqlite3_value_double(argv[0]);
    float y             = sqlite3_value_double(argv[1]);
    char *filename      = (char *)sqlite3_value_text(argv[2]);
    char *etiket        = (char *)sqlite3_value_text(argv[3]);
    int ip1             = sqlite3_value_int(argv[4]);
    int ip2             = sqlite3_value_int(argv[5]);
    int ip3             = sqlite3_value_int(argv[6]);
    char *typvar        = (char *)sqlite3_value_text(argv[7]);
    char *ez_opt        = (char *)sqlite3_value_text(argv[8]);
    
    int gdid = -1;
    float *fieldUU = NULL;
    float *fieldVV = NULL;
    char *nomvarUU  = strdup("UU");
    char *nomvarVV  = strdup("VV");
    float resultUU[1] ;
    float resultVV[1] ;
    
    udf_get_fstvar(filename, nomvarUU, typvar, etiket, ez_opt, ip1, ip2, ip3, &gdid, &fieldUU);
    udf_get_fstvar(filename, nomvarVV, typvar, etiket, ez_opt, ip1, ip2, ip3, &gdid, &fieldVV);

    //Do scalar interpolation at the latlon points
    c_ezsetopt("INTERP_DEGREE", ez_opt);
    if (c_gdxyvval(gdid, resultUU, resultVV, fieldUU, fieldVV, &x, &y, 2) < 0)
    {
        fprintf(stderr,"GDLLSVAL ERROR \n");
        exit(1);
    }
 
    char *buf = NULL;
    concat_result(&buf, resultUU[0], resultVV[0]);
 
    sqlite3_result_text(context, buf, -1, SQLITE_TRANSIENT);
}

static void udf_gdllwdval ( sqlite3_context *context, int argc, sqlite3_value **argv ) {
    
    if (argc != 9)
    {
        fprintf(stderr, "Error: wrong number of arguments \n");
        fprintf(stderr, "Usage : select gdxyvval lat lon filename etiket ip1 ip2 ip3 typevar nomvar ez_opt");
        exit(1);
    }
    
    float lat           = sqlite3_value_double(argv[0]);
    float lon           = sqlite3_value_double(argv[1]);
    char *filename      = (char *)sqlite3_value_text(argv[2]);
    char *etiket        = (char *)sqlite3_value_text(argv[3]);
    int ip1             = sqlite3_value_int(argv[4]);
    int ip2             = sqlite3_value_int(argv[5]);
    int ip3             = sqlite3_value_int(argv[6]);
    char *typvar        = (char *)sqlite3_value_text(argv[7]);
    char *ez_opt        = (char *)sqlite3_value_text(argv[8]);
    
    int gdid = -1;
    float *fieldUU = NULL;
    float *fieldVV = NULL;
    char *nomvarUV  = strdup("UU");
    char *nomvarWD  = strdup("VV");
    float resultUV[1] ;
    float resultWD[1] ;
    
    check_latlon(lat, lon);
    udf_get_fstvar(filename, nomvarUV, typvar, etiket, ez_opt, ip1, ip2, ip3, &gdid, &fieldUU);
    udf_get_fstvar(filename, nomvarWD, typvar, etiket, ez_opt, ip1, ip2, ip3, &gdid, &fieldVV);
    
    //Do scalar interpolation at the latlon points
    c_ezsetopt("INTERP_DEGREE", ez_opt);
    if (c_gdllwdval(gdid, resultUV, resultWD, fieldUU, fieldVV, &lat, &lon, 2) < 0)
    {
        fprintf(stderr,"GDLLSVAL ERROR \n");
        exit(1);
    }
    
    char *buf = NULL;
    concat_result(&buf, resultUV[0], resultWD[0]);
    
    sqlite3_result_text(context, buf, -1, SQLITE_TRANSIENT);
}

static void udf_gdxywdval ( sqlite3_context *context, int argc, sqlite3_value **argv ) {

    if (argc != 9)
    {
        fprintf(stderr, "Error: wrong number of arguments \n");
        fprintf(stderr, "Usage : select gdxyvval lat lon filename etiket ip1 ip2 ip3 typevar nomvar ez_opt");
        exit(1);
    }
    
    float x             = sqlite3_value_double(argv[0]);
    float y             = sqlite3_value_double(argv[1]);
    char *filename      = (char *)sqlite3_value_text(argv[2]);
    char *etiket        = (char *)sqlite3_value_text(argv[3]);
    int ip1             = sqlite3_value_int(argv[4]);
    int ip2             = sqlite3_value_int(argv[5]);
    int ip3             = sqlite3_value_int(argv[6]);
    char *typvar        = (char *)sqlite3_value_text(argv[7]);
    char *ez_opt        = (char *)sqlite3_value_text(argv[8]);
    
    int gdid = -1;
    float *fieldUU = NULL;
    float *fieldVV = NULL;
    char *nomvarUV  = strdup("UU");
    char *nomvarWD  = strdup("VV");
    float resultUV[1] ;
    float resultWD[1] ;
    
    udf_get_fstvar(filename, nomvarUV, typvar, etiket, ez_opt, ip1, ip2, ip3, &gdid, &fieldUU);
    udf_get_fstvar(filename, nomvarWD, typvar, etiket, ez_opt, ip1, ip2, ip3, &gdid, &fieldVV);
    
    //Do scalar interpolation at the latlon points
     c_ezsetopt("INTERP_DEGREE", ez_opt);
    if (c_gdxywdval(gdid, resultUV, resultWD, fieldUU, fieldVV, &x, &y, 2) < 0)
    {
        fprintf(stderr,"GDLLSVAL ERROR");
        exit(1);
    }
    
    char *buf = NULL;
    concat_result(&buf, resultUV[0], resultWD[0]);

    sqlite3_result_text(context, buf, -1, SQLITE_TRANSIENT);
}

int RegisterFuncs(struct sqlite3 *db){
    int rc;

    if ((rc = sqlite3_create_function(db, "crc32", 1, SQLITE_UTF8, NULL, &crc32_func, NULL, NULL) ))
       fprintf(stderr,"Problem with crc32\n");

    if ((rc = sqlite3_create_function(db, "mask_mer", 2, SQLITE_UTF8, NULL, &mask_mer, NULL, NULL) ))
       fprintf(stderr,"Problem with numBox\n");
    
    if ((rc = sqlite3_create_function(db, "mask_glace", 2, SQLITE_UTF8, NULL, &mask_glace, NULL, NULL) ))
       fprintf(stderr,"Problem with maskGlace\n");

    if ((rc = sqlite3_create_function(db, "mask_glace_clim", 2, SQLITE_UTF8, NULL, &mask_glace_clim, NULL, NULL) ))
       fprintf(stderr,"Problem with mask_glace_clim\n");
    
    if ((rc = sqlite3_create_function(db, "gdllsval", 10, SQLITE_UTF8, NULL, &udf_gdllsval, NULL, NULL) ))
       fprintf(stderr,"Problem with udf_gdllsval\n");

    if ((rc = sqlite3_create_function(db, "gdxysval", 10, SQLITE_UTF8, NULL, &udf_gdxysval, NULL, NULL) ))
       fprintf(stderr,"Problem with udf_gdxysval\n");
    
    if ((rc = sqlite3_create_function(db, "gdllvval", 9, SQLITE_UTF8, NULL, &udf_gdllvval, NULL, NULL) ))
       fprintf(stderr,"Problem with udf_gdllvval\n");
    
    if ((rc = sqlite3_create_function(db, "gdxyvval", 9, SQLITE_UTF8, NULL, &udf_gdxyvval, NULL, NULL) ))
       fprintf(stderr,"Problem with udf_gdxyvval\n");
    
    if ((rc = sqlite3_create_function(db, "gdllwdval", 9, SQLITE_UTF8, NULL, &udf_gdllwdval, NULL, NULL) ))
       fprintf(stderr,"Problem with udf_gdllwdval\n");
    
    if ((rc = sqlite3_create_function(db, "gdxywdval", 9, SQLITE_UTF8, NULL, &udf_gdxywdval, NULL, NULL) ))
       fprintf(stderr,"Problem with udf_gdxywdval\n");
    
    if ((rc = sqlite3_create_function(db, "numBox", 4, SQLITE_UTF8, NULL, &myRegrupBoxFunc, NULL, NULL) ))
       fprintf(stderr,"Problem with numBox\n");

    if ((rc = sqlite3_create_function(db, "Lat100Box", 4, SQLITE_UTF8, NULL, &myRegrupLat100Func, NULL, NULL) ))
       fprintf(stderr,"Problem with Lat100Box\n");

    if ((rc = sqlite3_create_function(db, "Lon100Box", 4, SQLITE_UTF8, NULL, &myRegrupLon100Func, NULL, NULL) ))
       fprintf(stderr,"Problem with Lon100Box\n");

    if ((rc = sqlite3_create_function(db, "LatBox", 4, SQLITE_UTF8, NULL, &myRegrupLatFunc, NULL, NULL) ))
       fprintf(stderr,"Problem with LatBox\n");

    if ((rc = sqlite3_create_function(db, "LonBox", 4, SQLITE_UTF8, NULL, &myRegrupLonFunc, NULL, NULL) ))
       fprintf(stderr,"Problem with LonBox\n");
    if ((rc = sqlite3_create_function(db, "StampBox",2,SQLITE_UTF8, NULL, &myTimeStampFunc, NULL, NULL) ))
       fprintf(stderr,"Problem with StampBox\n");

    
    if ((rc = sqlite3_create_function(db, "BlobInt", 2, SQLITE_UTF8, NULL, &getIntFromBlob, NULL, NULL) ))
         fprintf(stderr,"Problem with GET_INT\n");
    if ((rc = sqlite3_create_function(db, "BlobDouble", 2, SQLITE_UTF8, NULL, &getDoubleFromBlob, NULL, NULL) ))
         fprintf(stderr,"Problem with GET_DOUBLE\n");
    
    
    if ((rc = sqlite3_create_function(db, "Kms", -1, SQLITE_UTF8, NULL,
             &KmsFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with KmsFunc\n");

    if ((rc = sqlite3_create_function(db, "Miles", -1, SQLITE_UTF8, NULL,
             &MilesFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with MilesFunc\n");

    if ((rc = sqlite3_create_function(db, "Miles2", -1, SQLITE_UTF8, NULL,
             &Miles2Func, NULL, NULL) != 0))
         fprintf(stderr,"Problem with MilesFunc\n");

    if ((rc = sqlite3_create_function(db, "isodatetime", -1, SQLITE_UTF8, NULL,
             &isoDateTimeFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with isoDateTimeFunc\n");

    if ((rc = sqlite3_create_function(db, "isodate", -1, SQLITE_UTF8, NULL,
             &isoDateFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with isoDateFunc\n");

    if ((rc = sqlite3_create_function(db, "isotime", -1, SQLITE_UTF8, NULL,
             &isoTimeFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with isoDateFunc\n");


    if ((rc = sqlite3_create_function(db, "strleft", 2, SQLITE_UTF8, NULL,
             &LeftFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with strLeft\n");

    if ((rc = sqlite3_create_function(db, "concat", -1, SQLITE_UTF8, NULL,
             &ConcatFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with concat\n");

    if ((rc = sqlite3_create_function(db, "sign", 1, SQLITE_UTF8, NULL,
             &msignFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with sign\n");

    if ((rc = sqlite3_create_function(db, "m_abs", 1, SQLITE_UTF8, NULL,
             &mAbsFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with msign\n");

    if ((rc = sqlite3_create_function(db, "sqrt", 1, SQLITE_UTF8, NULL,
             &SqrtFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with sqrt\n");

    if ((rc = sqlite3_create_function(db, "PI", 0, SQLITE_UTF8, NULL,
             &PIFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with PI\n");

    if ((rc = sqlite3_create_function(db, "pow", 2, SQLITE_UTF8, NULL,
             &PowFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with pow\n");

    if ((rc = sqlite3_create_function(db, "sin", 1, SQLITE_UTF8, NULL,
             &SinFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with sin\n");

    if ((rc = sqlite3_create_function(db, "asin", 1, SQLITE_UTF8, NULL,
             &AsinFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with asin\n");

    if ((rc = sqlite3_create_function(db, "cos", 1, SQLITE_UTF8, NULL,
             &CosFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with cos\n");

    if ((rc = sqlite3_create_function(db, "acos", 1, SQLITE_UTF8, NULL,
             &AcosFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with acos\n");

    if ((rc = sqlite3_create_function(db, "floor", 1, SQLITE_UTF8, NULL,
             &FloorFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with floor\n");

    if ((rc = sqlite3_create_function(db, "ceil", 1, SQLITE_UTF8, NULL,
             &CeilFunc, NULL, NULL) != 0))
         fprintf(stderr,"Problem with ceil\n");

    if ((rc = sqlite3_create_function(db, "variance", 1, SQLITE_UTF8, NULL,
             NULL, &VarianceStep, &VarianceFinalize) != 0))
         fprintf(stderr,"Problem with Variance\n");

    if ((rc = sqlite3_create_function(db, "var", 1, SQLITE_UTF8, NULL,
             NULL, &VarianceStep, &VarianceFinalize) != 0))
         fprintf(stderr,"Problem with Variance\n");

    if ((rc = sqlite3_create_function(db, "STDDEV", 1, SQLITE_UTF8, NULL,
             NULL, &VarianceStep, &StddevFinalize) != 0))
         fprintf(stderr,"Problem with STDDEV\n");

    if ((rc = sqlite3_create_function(db, "XSTD", 1, SQLITE_UTF8, NULL,
             NULL, &XSTDStep, &XSTDFinalize) != 0))
         fprintf(stderr,"Problem with XSTD\n");
    if ((rc = sqlite3_create_function(db, "XVAR", 1, SQLITE_UTF8, NULL,
             NULL, &XSTDStep, &XVARFinalize) != 0))
         fprintf(stderr,"Problem with XVAR\n");
   
    return rc;
}

// make this library available as a Run-Time Loadable Extension
// from http://www.sqlite.org/loadext.html
// change by Chris Malek for Pierre Koclas
#ifdef _WIN32
__declspec(dllexport)
#endif
int sqlite3_extension_init( /* <== Change this name, maybe */
  sqlite3 *db, 
  char **pzErrMsg, 
  const sqlite3_api_routines *pApi
){
  int rc = SQLITE_OK;
  SQLITE_EXTENSION_INIT2(pApi);
  /* insert code to initialize your extension here */
  rc = RegisterFuncs(db);
  return rc;
}

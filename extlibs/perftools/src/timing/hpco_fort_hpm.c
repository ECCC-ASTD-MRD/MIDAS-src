//#include <rpnmacros.h>
#ifndef F2Cl
#define F2Cl int
#endif
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>

/* 
   code for the tmg package, originally written in fortran
   this package is inspired by the hpm timing package that
   is used on IBM Power systems under AIX  (need libhpm_r.a)

   to activate the tmg timing package
   export TMG_ON=YES       (report in seconds)
   export TMG_ON=TSC       (use the time stamp counter on X86 and X86_64 systems - experimental - unstable)

   call tmg_init(pe_id,pname)   tmg_init_c(int pe_id, char *pname)   init package
   call tmg_terminate(pe_id)    tmg_terminate_c(int pe_id)           terminate and print
   call tmg_start(id,sname)     tmg_start_c(int id, char *sname)     start section id
   call tmg_tstart(id,sname)    tmg_tstart_c(int id, char *sname)    start threaded section id
   call tmg_stop(id)            tmg_stop_c(int id)                   end section id
   call tmg_tstop(id)           tmg_tstop_c(int id)                  end threaded section id
   INTEGER pe_id,id
   CHARACTER *(*) pname,sname

   NOTES:
        threaded sections are only supported for AIX  / Power for the time being
        the tmg_xxxx C routines call the equivalent FORTRAN routines
        environment variable HPM_NUM_INST_PTS may be used to set the highest id available (default is 200)
        (this is compatible with native AIX hpm library)

        a direct gateway to the IBM hpm package is available  (see IBM doc for routine behaviour)
        call fort_hpminit(pe_id,pname)           c_hpminit(int pe_id, char *pname)
        call fort_hpmterminate(pe_id)            c_hpmterminate(int  pe_id)
        call fort_hpnstart(id,sname)             c_hpmstart(int id, char *sname)
        call fort_hpntstart(id,sname)            c_hpmtstart(int id, char *sname)
        call fort_hpmstop(id)                    c_hpmstop(int id)
        call fort_hpmtstop(id)                   c_hpmtstop(int id)

   EXTRAS:
   ticks=time_base()                 integer *8 / long long (return time in ticks)
   time=rtools_wtime()               real *8    / double    (return time in seconds, uses timebase or gettimeofday())
   nsecs=time_elapsed()              integer *8 / long long (return time in nanoseconds using gettimeofday())
   nsecs=time_base_to_nsec(ticks)    integer *8 / long long (return nanoseconds from ticks)
   call use_timeofday()              use gettimeofday for timings
   call use_timebase()               use hardware time base if available
*/

static int MAX_ENTRIES=200;

static int tmg_ignore=1;  /* tmg package off by default */
static int not_initialized=1;

static double *tmg_ot, *tmg_ct, c1, r1;
static int *calls;
static char **msg;

static int my_pe=0;

static int use_tsc=0;

void fort_hpminit_(int *taskid, char *name, F2Cl lname);
void fort_hpmterminate_(int *taskid);
void fort_hpmstop_(int *taskid);
void fort_hpmtstop_(int *taskid);
void fort_hpmstart_(int *taskid, char *name, F2Cl lname);
void fort_hpmtstart_(int *taskid, char *name, F2Cl lname);

double rtools_wtime_();
static double ticks_to_sec  = 0.0;
static double ticks_to_nsec = 0.0;
static double ov_cpu_freq   = 0.0;  /* inverse of cpu clock frequency */

void use_timeofday()   { use_tsc=0; }
void use_timeofday_()  { use_tsc=0; }
void use_timeofday__() { use_tsc=0; }

void use_timebase()   { use_tsc=1; }
void use_timebase_()  { use_tsc=1; }
void use_timebase__() { use_tsc=1; }

#if ! defined(_AIX)
int mclock()
{
  clock_t returnvalue=clock();
  return(returnvalue);
}
int mclock_() { return mclock() ;}
int mclock__() { return mclock() ;}
#endif

unsigned long long rtools_clock()
{
  long long t1=clock();
  return(t1);
}
unsigned long long rtools_clock_() { return rtools_clock() ; }
unsigned long long rtools_clock__() { return rtools_clock() ; }

double rtools_cp_time()
{
  double t1;
#if defined(_AIX)
  t1 = clock();
  return(t1 / CLOCKS_PER_SEC);
#else
  struct rusage usage;
  if( getrusage(RUSAGE_SELF,&usage) ) return(0);
  t1  = usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1.0e-6;
  t1 += usage.ru_stime.tv_sec + usage.ru_stime.tv_usec * 1.0e-6;
  return(t1);
#endif
}
double rtools_cp_time_() { return rtools_cp_time() ; }
double rtools_cp_time__() { return rtools_cp_time() ; }

static float get_cpu_freq()  /* return frequency in GHz for Linux, 1.9 otherwise */
{
  char buffer[1024];
  char *pbuf=buffer;
  float freq=1900.0;  /* default at 1900 MHz */
  FILE *cpuinfo;
#ifdef linux
  cpuinfo=fopen("/proc/cpuinfo","r");
  if(cpuinfo != NULL) {
    while(1) {
      fgets(buffer,sizeof(buffer),cpuinfo);
      if(buffer[4]=='M' && buffer[5]=='H' && buffer[6]=='z') {
        while(*pbuf != ':') pbuf++; pbuf++;
        sscanf(pbuf,"%f",&freq);
        break;
      }
    }
  }
#endif
  return .001 * freq;
}

static void init_ticks_to_sec()
{
#if defined(_AIX)
  timebasestruct_t t_s;
  struct timeval tp1, tp2;
  unsigned long long time1, time2;
  double dtime1, dtime2;

  read_real_time(&t_s, TIMEBASE_SZ); /* get the timebasestruct_t initialized properly */
  __fence();
  time1=__mftb();
  __fence();
  gettimeofday(&tp1, NULL);
  gettimeofday(&tp2, NULL);
  __fence();
  time2=__mftb();
  __fence();
  time2=time2 - time1;
  t_s.tb_high= time2 >> 32;
  t_s.tb_low= time2 & 0xFFFFFFFF;
//  printf("high=%d , low=%d \n",t_s.tb_high,t_s.tb_low);
  time_base_to_time(&t_s, TIMEBASE_SZ);
//  printf("high=%d , low=%d \n",t_s.tb_high,t_s.tb_low);
  time1=t_s.tb_high*1000000000 + t_s.tb_low ;
  dtime1 = time1;
  dtime2 = time2;
//  printf(" %d ticks = %d nanoseconds \n",time2,time1);
  dtime2 = dtime1 / dtime2;
//  printf("clock ratio is %f \n",dtime2);
  ticks_to_sec = dtime2 * 1e-9;
#else
  ticks_to_sec = ov_cpu_freq ;
#endif
  return;
}

static void initialize_timers()
{
  char *tmg_is_on;

  if(not_initialized == 0)return;

  tmg_is_on=getenv("TMG_ON");
  if(tmg_is_on != NULL) {
    tmg_ignore = strcmp(tmg_is_on,"YES");
    if(0 == strcmp(tmg_is_on,"TSC")) { /* meaningful only on X86 or X86_64 systems */
      use_tsc=1 ;
      tmg_ignore=0 ;
      fprintf(stderr,"Using hardware time stamp counter \n");
    }
  }
  ov_cpu_freq = 1e-9 / get_cpu_freq();
  init_ticks_to_sec();
  ticks_to_nsec = ticks_to_sec * 1e9;
  not_initialized=0;
printf("Timers are now initialized\n");
}

long long int time_elapsed_()  /* elapsed time in nanoseconds from get time of day */
{
  static long long last=0;  /* kept in nanoseconds, even if gettimeofday has at best microsecond resolution */
  long long int elapsed;
  struct timeval tp;
  gettimeofday(&tp, NULL);
  elapsed=tp.tv_sec;
  elapsed=elapsed*1000000000;
  elapsed=elapsed+tp.tv_usec*1000;
  if (last == 0 ) { last = elapsed ;}  /* initialize last */
  if (elapsed<=last) elapsed=++last;   /* never return same value twice */
  return elapsed;  /* time in nanoseconds */
}
long long int time_elapsed()  { return time_elapsed_() ; }
long long int time_elapsed__() { return time_elapsed_() ; }
long long int c_time_elapsed()  { return time_elapsed_() ; }

unsigned long long time_base_()  /* return raw clock ticks from appropriate registers (if possible) */
{
#if defined(_AIX)
  return __mftb();   /* power series mftb instruction */
#else
  #if defined(i386) || defined(__x86_64)
    unsigned long lo, hi;
    unsigned long long result;  /* use the infamous rdtsc instruction to get the time stamp counter */
    __asm__ __volatile__( "rdtsc" : "=a" (lo), "=d" (hi) );
    result = hi;
    result = lo | (result<<32);
    if (use_tsc) return result ;
    else return time_elapsed() ;
  #else
    return time_elapsed();  /* time in nanoseconds */

  #endif
#endif
}
unsigned long long time_base()  { return time_base_() ; } 
unsigned long long time_base__() { return time_base_() ; } 
unsigned long long c_time_base()  { return time_base_() ; } 

unsigned long long time_base_to_nsec(unsigned long long *ticks)  /* ticks to nanoseconds conversion */
{
  unsigned long long temp=*ticks;
#if defined(_AIX)
  timebasestruct_t t_start;
  int sec, n_sec;
  read_real_time(&t_start, TIMEBASE_SZ);  /* get the timebasestruct_t initialized properly */
  t_start.tb_high= temp >> 32;
  t_start.tb_low= temp & 0xFFFFFFFF;
  time_base_to_time(&t_start, TIMEBASE_SZ);
  sec   = t_start.tb_high;
  n_sec = t_start.tb_low;
  temp = (long long int) sec*1000000000 + (long long int)n_sec;
#else
  #if defined(i386) || defined(__x86_64)
    double total_ticks;
    struct timeval tp;  /* let's do our best conversion effort here */
    if(not_initialized)initialize_timers();
    if(use_tsc) {       /* if time stamp counter is not used, no conversion is needed */
      total_ticks = temp;
      total_ticks=total_ticks*ticks_to_nsec;
      temp = total_ticks;
    }
  #endif
#endif
  return temp;  /* time in nanoseconds */
}
unsigned long long time_base_to_nsec_(unsigned long long *ticks) { return time_base_to_nsec(ticks); }
unsigned long long time_base_to_nsec__(unsigned long long *ticks) { return time_base_to_nsec(ticks); }

/*
  special timing functions needed on non POWERPC/AIX machines
  copied from rtools_wtime and declared static to avoid conflicts
*/
#if defined(_AIX)
double rtools_wtime_()
{
/*
  long long t1 = time_base();
*/
  long long t1;
  __fence();
  t1 = __mftb();
  __fence();
  if(not_initialized)initialize_timers();
  return(t1*ticks_to_sec);
}
#else

#include<time.h>
#include <sys/times.h>
#include <sys/time.h>
#include <sys/resource.h>

//static long ticks_to_sec=0;
static clock_t tim0=0;

double rtools_wtime_()
{
  unsigned long long t1;
  double t2;

  if(not_initialized)initialize_timers();
  if(use_tsc) {
    t1=time_base_();
    t2=t1*ov_cpu_freq;
  }else{
    t1=time_elapsed_();
    t2=t1*1e-9;
  }
  return(t2);
}
#endif
double rtools_wtime__() { return rtools_wtime_() ; }
double rtools_wtime() { return rtools_wtime_() ; }

/*
  the tmg package itself and its various routines
*/
static void dump_tmg_entries()
{
  int i;
  double average;

  printf("TIMING TABLE ENTRIES\n");
  for ( i=0 ; i<MAX_ENTRIES ; i++ ) {
    if(msg[i] != NULL) {
      average = calls[i];
      if(calls[i] != 0) average = tmg_ct[i]/average;
      else average = tmg_ct[i];
      printf("TMG: PE=%4d ,LABEL=%s ,COUNT=%d ,TIME=%20.6f seconds, AVG=%20.6f seconds\n",my_pe,msg[i],calls[i],tmg_ct[i],average);
    }
  }
}
static char *cstring_from_fstring(char *fstring, int lstring)
{
  char *temp = (char *)malloc(lstring+1);
  int i = lstring-1;
  strncpy(temp,fstring,lstring);
  temp[lstring]='\0';
  while( (i > 0) && (temp[i] == ' ') ) { temp[i]='\0' ; i-- ; };
  return(temp);
}

/* the following code was initially a quick and dirty fix to get rid
   once and for all of unsatisfied externals with
   a trailing underscore when using the IBM hpm package in -qextname mode.
   a new set of FORTRAN entry points has been created that call the equivalent
   C entry points. this also gets rid of the FORTRAN macro file f_hpm.h
   fort_hpminit replaces f_hpminit
   fort_hpmterminate replaces f_hpmterminate
   fort_hpmstart replaces f_hpmstart
   fort_hpmtstart replaces f_hpmtstart
   fort_hpmstop replaces f_hpmstop
   fort_hpmtstop replaces f_hpmtstop
   the calling sequences are unchanged

   expanded to provide these entry points on all platforms

   expanded to provide an equivalent set of entry points for c
   c_hpminit c_hpmterminate c_hpmstart c_hpmtstart c_hpmstop c_hpmtstop

   expanded to provide a tmg_xxx_c set of entry points for c
   tmg_init_c tmg_terminate_c tmg_start_c tmg_tstart_c tmg_stop_c tmg_tstop_c
*/

void tmg_init_(int *id, char *string, F2Cl string_len)
{
  if(not_initialized)initialize_timers();
  if(tmg_ignore) return;
  fort_hpminit_(id, string, string_len);
}
void tmg_init__(int *id, char *string, F2Cl string_len) { tmg_init_(id, string, string_len) ; }
void tmg_init(int *id, char *string, F2Cl string_len) { tmg_init_(id, string, string_len) ; }

void tmg_terminate_(int *id)
{
  if(tmg_ignore) return;
  fort_hpmterminate_(id);
}
void tmg_terminate__(int *id) { tmg_terminate_(id) ; }
void tmg_terminate(int *id) { tmg_terminate_(id) ; }

void tmg_start_(int *id, char *string, F2Cl string_len)
{
  if(tmg_ignore) return;
  fort_hpmstart_(id, string, string_len);
}
void tmg_start__(int *id, char *string, F2Cl string_len) { tmg_start_(id, string, string_len) ; }
void tmg_start(int *id, char *string, F2Cl string_len) { tmg_start_(id, string, string_len) ; }

void tmg_tstart_(int *id, char *string, F2Cl string_len)
{
  if(tmg_ignore) return;
  fort_hpmtstart_(id, string, string_len);
}
void tmg_tstart__(int *id, char *string, F2Cl string_len) { tmg_tstart_(id, string, string_len) ; }
void tmg_tstart(int *id, char *string, F2Cl string_len) { tmg_tstart_(id, string, string_len) ; }

void tmg_stop_(int *id)
{
  if(tmg_ignore) return;
  fort_hpmstop_(id);
}
void tmg_stop__(int *id) { tmg_stop_(id) ; }
void tmg_stop(int *id) { tmg_stop_(id) ; }

void tmg_tstop_(int *id)
{
  if(tmg_ignore) return;
  fort_hpmtstop_(id);
}
void tmg_tstop__(int *id) { tmg_tstop_(id) ; }
void tmg_tstop(int *id) { tmg_tstop_(id) ; }

/* tmg_xxx_c series of entry points , they are just gateways to the fortran tmg_xxx entry points*/

void tmg_init_c(int taskid, char *localname)
{
  int ftaskid=taskid;
  F2Cl lname=strlen(localname);
  tmg_init_(&ftaskid,localname,lname);
}
void tmg_terminate_c(int taskid)
{
  int ftaskid=taskid;
  tmg_terminate_(&ftaskid);
}
void tmg_stop_c(int taskid)
{
  int ftaskid=taskid;
  tmg_stop_(&ftaskid);
}
void tmg_tstop_c(int taskid)
{
  int ftaskid=taskid;
  tmg_tstop_(&ftaskid);
}
void tmg_start_c(int taskid, char *localname)
{
  int ftaskid=taskid;
  F2Cl lname=strlen(localname);
  tmg_start_(&ftaskid,localname,lname);
}
void tmg_tstart_c(int taskid, char *localname)
{
  int ftaskid=taskid;
  F2Cl lname=strlen(localname);
  tmg_tstart_(&ftaskid,localname,lname);
}

 void DumhpmInit(int taskid, char *localname)
{
  int i;
  char *num_inst_pts=getenv("HPM_NUM_INST_PTS");
  if(not_initialized)initialize_timers();
  if(num_inst_pts != NULL) {
    sscanf(num_inst_pts,"%d",&MAX_ENTRIES);
    printf("setting MAX_ENTRIES to %d\n",MAX_ENTRIES);
  }
  MAX_ENTRIES++;
  calls  = (int    *)malloc(MAX_ENTRIES*sizeof(int  *));   /* number of timing calls for this region */
  msg    = (char  **)malloc(MAX_ENTRIES*sizeof(char *));   /* name of region */
  tmg_ot = (double *)malloc(MAX_ENTRIES*sizeof(double));   /* start time, acts as an already started flag */
  tmg_ct = (double *)malloc(MAX_ENTRIES*sizeof(double));   /* accumulated time */
  r1 = 1000000.0;
  my_pe = taskid ;
  for ( i=0 ; i< MAX_ENTRIES ; i++ ) { tmg_ot[i] = -1.0 ; tmg_ct[i] = 0.0 ; msg[i] = NULL ; calls[i]=0; }
}

 void DumhpmTerminate(int taskid)
{
  dump_tmg_entries();
}

 void DumhpmStart(int taskid, char *localname)
{
  if( taskid >= MAX_ENTRIES ) return;  /* table is too small */
  if( msg[taskid] == NULL ) {          /* if not already done, take a copy of the name string */
    char *temp = (char *)malloc(strlen(localname)+1) ;
    if(temp != NULL){
      strncpy(temp,localname,strlen(localname)+1);
      msg[taskid] = temp;
    }else{
      fprintf(stderr,"ERROR: cannot allocate string memory\n");
      return;
    }
  }else{
    if(0 != strncmp(msg[taskid],localname,strlen(localname))) {
      fprintf(stderr,"ERROR: id=%d mismatch, expected '%s' , got '%s'\n",taskid,msg[taskid],localname);
      return;
    }
  }
  if( tmg_ot[taskid] >= 0 ) {
    fprintf(stderr,"ERROR: id=%d timing region '%s' already started\n",taskid,msg[taskid]);
    return;    /* error: start already done */
  }
  tmg_ot[taskid] = rtools_wtime_();
  calls[taskid]++;
}
 void DumhpmTstart(int taskid, char *localname) { DumhpmStart(taskid,localname); }

 void DumhpmStop(int taskid)
{
  if( taskid >= MAX_ENTRIES ) return;  /* table is too small */
  if( tmg_ot[taskid] < 0 ) {
    fprintf(stderr,"ERROR: id=%d timing region not started\n",taskid);
    return;     /* error: no start was done */
  }
  tmg_ct[taskid] = tmg_ct[taskid] + (rtools_wtime_() - tmg_ot[taskid]);
  tmg_ot[taskid] = -1.0;
}
 void DumhpmTstop(int taskid)                   { DumhpmStop(taskid); }

#if defined(_AIX)
#if defined(_ARCH_PWR7)
#include <libhpc.h>
#else
#include <libhpm.h>
#endif
#else
static void hpmInit(int taskid, char *localname)
{
  DumhpmInit(taskid, localname);
}

static void hpmTerminate(int taskid)
{
  DumhpmTerminate(taskid);
}

static void hpmStart(int taskid, char *localname)
{
  DumhpmStart(taskid, localname);
}
static void hpmTstart(int taskid, char *localname) { DumhpmStart(taskid,localname); }

static void hpmStop(int taskid)
{
  DumhpmStop(taskid);
}
static void hpmTstop(int taskid)                   { DumhpmStop(taskid); }
#endif

/*
  c_hpmxxx entry points and fort_hpmxxx entry points for IBM and other platforms
  the IBM version of the code calls the hpm package directly
*/
void c_hpminit(int taskid, char *localname)
{
  hpmInit(taskid,localname);
}
void c_hpmterminate(int taskid)
{
  hpmTerminate(taskid);
}
void c_hpmstart(int taskid, char *localname)
{
  hpmStart(taskid,localname);
}
void c_hpmtstart(int taskid, char *localname)
{
  hpmTstart(taskid,localname);
}
void c_hpmstop(int taskid)
{
  hpmStop(taskid);
}
void c_hpmtstop(int taskid)
{
  hpmTstop(taskid);
}

void fort_hpminit_(int *taskid, char *name, F2Cl lname)
{
  char *localname;
  int len = lname;

  localname=cstring_from_fstring(name,len);
  hpmInit(*taskid,localname);
  free(localname);
}
void fort_hpminit__(int *taskid, char *name, F2Cl lname) { fort_hpminit_(taskid, name, lname) ; }
void fort_hpminit(int *taskid, char *name, F2Cl lname) { fort_hpminit_(taskid, name, lname) ; }

void fort_hpmterminate_(int *taskid)
{
  hpmTerminate(*taskid);
}
void fort_hpmterminate__(int *taskid) { fort_hpmterminate_(taskid) ; }
void fort_hpmterminate(int *taskid) { fort_hpmterminate_(taskid) ; }

void  fort_hpmstop_(int *taskid)
{
  hpmStop(*taskid);
}
void  fort_hpmstop__(int *taskid) { fort_hpmstop_(taskid) ; }
void  fort_hpmstop(int *taskid) { fort_hpmstop_(taskid) ; }

void fort_hpmtstop_(int *taskid)
{
  hpmTstop(*taskid);
}
void  fort_hpmtstop__(int *taskid) { fort_hpmtstop_(taskid) ; }
void  fort_hpmtstop(int *taskid) { fort_hpmtstop_(taskid) ; }

void fort_hpmstart_(int *taskid, char *name, F2Cl lname)
{
  char *localname;
  int len = lname;

  localname=cstring_from_fstring(name,len);
  hpmStart(*taskid,localname);
  free(localname);
}
void fort_hpmstart__(int *taskid, char *name, F2Cl lname) { fort_hpmstart_(taskid, name, lname) ; }
void fort_hpmstart(int *taskid, char *name, F2Cl lname) { fort_hpmstart_(taskid, name, lname) ; }

void fort_hpmtstart_(int *taskid, char *name, F2Cl lname)
{
  char *localname;
  int len = lname;

  localname=cstring_from_fstring(name,len);
  hpmTstart(*taskid,localname);
  free(localname);
}
void fort_hpmtstart__(int *taskid, char *name, F2Cl lname) { fort_hpmtstart_(taskid, name, lname) ; }
void fort_hpmtstart(int *taskid, char *name, F2Cl lname) { fort_hpmtstart_(taskid, name, lname) ; }

/* a small idiot test, uses the C entry points because they call the FORTRAN ones */
/*
  AIX compilation
  mpcc_r -I$ARMNLIB/include -I$ARMNLIB/include/AIX -DTEST -DMPI -o zarza hpco_fort_hpm.c -lhpm_r
  Linux compilation
  s.cc -mpi -DTEST -DMPI  hpco_fort_hpm.c -o zarza
  NOTE:
  -DMPI is optional
*/
#if defined(TEST)
#ifdef MPI
#include <mpi.h>
#endif
main(int argc, char **argv)
{
  int junk=1;
  double point1;
  double point2;
  int i;
  double temp , sum;
  double tclk1, tclk2;
//  init_ticks_to_sec();
  point1=rtools_wtime();
  printf("Point1 = %g \n",point1);
#ifdef MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&junk);
  fprintf(stderr,"my rank is %d\n",junk);
#endif
  printf("Before tmg_init_c\n");
  tmg_init_c(junk,"toto");
  printf("Before tmg_start_c 3\n");
  tmg_start_c(3,"titi");
  printf("Before tmg_stop_c 3\n");
  tmg_stop_c(3);
//  tmg_stop_c(3);
  printf("Before tmg_start_c 2\n");
  tmg_start_c(2,"tata");
//  tmg_start_c(2,"tata");
  sleep(1);
  printf("Before tmg_stop_c 2\n");
  tmg_stop_c(2);
  tmg_start_c(3,"titi");
  tmg_stop_c(3);
  tmg_start_c(3,"titi");
  tmg_stop_c(3);
  tmg_start_c(4,"tata");
  tclk1=rtools_cp_time();
  sum = 0.0;
  printf("sum = %lg\n",sum);
  for ( i=0 ; i<100000000 ; i++ ) { temp = i*1.00 ; sum = sum + temp*temp ;}
  printf("sum = %lg\n",sum);
  tmg_stop_c(4);
  tclk2=rtools_cp_time() - tclk1;
  printf("CPU = %g\n",tclk2);
  printf("Before tmg_terminate_c\n");
  tmg_terminate_c(junk);
  point2=rtools_wtime_();
  printf("point2 = %g , point1 = %g \n",point2,point1);
  point2 = point2 - point1;
  printf("Total time = %f seconds , ov_cpu_freq= %g \n",point2,ov_cpu_freq);
#ifdef MPI
  MPI_Finalize();
#endif
}
#endif

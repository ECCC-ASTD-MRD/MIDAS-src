#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>

#define TABLE_SIZE 8200
typedef struct {
  long long value;    /* time tag */
  int pid;            /* process pid */
  char tag[12];       /* 12 character tag */
} table_entry;

static table_entry table[TABLE_SIZE] ;
static int table_index=-1;
static int log_use_elapsed=0;
static int binary_output=1;
static int my_pid=0;
static long long history_offset_nano=0;
static long long history_offset_tsc=0;

void log_use_gettimeofday()
{
  log_use_elapsed=1;
}
void log_use_gettimeofday_()
{
  log_use_elapsed=1;
}
void log_use_gettimeofday__()
{
  log_use_elapsed=1;
}

void log_use_tsc()
{
  log_use_elapsed=0;
}
void log_use_tsc_()
{
  log_use_elapsed=0;
}
void log_use_tsc__()
{
  log_use_elapsed=0;
}

unsigned long long time_elapsed_nanosec()  /* elapsed time in nanoseconds from get time of day */
{
  static long long last=0;  /* kept in nanoseconds, even if gettimeofday has at best microsec resolution */
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
unsigned long long time_elapsed_nanosec_() { return(time_elapsed_nanosec()) ; }
unsigned long long time_elapsed_nanosec__() { return(time_elapsed_nanosec()) ; }

unsigned long long get_history_time()
{
   register char *temp;
   register unsigned long long time;

#if defined(_AIX)
  if(log_use_elapsed)
    time=time_elapsed_nanosec() - history_offset_nano;
  else
    time=__mftb() - history_offset_tsc;
#else
  #if defined(i386) || defined(__x86_64)
    if(log_use_elapsed){
      time=time_elapsed_nanosec() - history_offset_nano;  /* time in nanoseconds */
    }else{
      unsigned long lo, hi;
      __asm__ __volatile__( "rdtsc" : "=a" (lo), "=d" (hi) );
      time = hi;
      time = lo | (time<<32);
      time = time - history_offset_tsc;
    }
  #else
    time=time_elapsed_nanosec() - history_offset_nano;  /* time in nanoseconds */
  #endif

#endif
  return time;
}
unsigned long long get_history_time_()  { return get_history_time() ; }
unsigned long long get_history_time__() { return get_history_time() ; }

void set_history_offset(unsigned long long *offset )
{
  if(log_use_elapsed) history_offset_nano = *offset;
  else                history_offset_tsc  = *offset;
}
void set_history_offset_(unsigned long long *offset)  { set_history_offset(offset); }
void set_history_offset__(unsigned long long *offset) { set_history_offset(offset); }

unsigned long long log_history_trace_(char *tag)
{
   register char *temp;
   unsigned long long time;
   
   if (my_pid==0) {   /* initialization */
      my_pid=getpid();
      temp=getenv("HISTORY_TRACE_CLOCK");
      if(temp != NULL){
         if(*temp == 'S') log_use_tsc();   /* export HISTORY_TRACE_CLOCK=S */
         else             log_use_gettimeofday();
      }
   }
   
   if (table_index>=TABLE_SIZE) return(-1);  // OOPS, table overflow
   table_index++;
   temp=table[table_index].tag;
   *temp++=*tag++; *temp++=*tag++; *temp++=*tag++; *temp++=*tag++;
   *temp++=*tag++; *temp++=*tag++; *temp++=*tag++; *temp++=*tag++;
   *temp++=*tag++; *temp++=*tag++; *temp++=*tag++; *temp++=*tag++;
   time=get_history_time();
   table[table_index].value=time;
   table[table_index].pid=my_pid;
   return(time);
}
unsigned long long log_history_trace(char *tag)
{
  return log_history_trace_(tag);
}
unsigned long long log_history_trace__(char *tag)
{
  return log_history_trace_(tag);
}

static int fd_out=-1234;
static FILE *file_out=NULL;

void write_history_trace_()
{
  int nentries=1+table_index;
  int nbytes;
  char *mpi_child;
  int mpi_rank=-1;
  char msgfile_name[4096];
  int i;

  if(fd_out == -1234){
    char *fname=getenv("HISTORY_TRACE_FILE");
    if(fname == NULL) { fd_out = -2 ; return ; }
    if( *fname == '+' ) { fname++ ; binary_output=0; };

    mpi_child = getenv("MP_CHILD");  // AIX
    if(mpi_child != NULL){
      mpi_rank=atoi(mpi_child);
    }
    mpi_child = getenv("PMI_RANK");  // mpich2
    if(mpi_child != NULL){
      mpi_rank=atoi(mpi_child);
    }
    mpi_child = getenv("OMPI_COMM_WORLD_RANK");  // open mpi
    if(mpi_child != NULL){
      mpi_rank=atoi(mpi_child);
    }
    if(mpi_rank < 0)
      snprintf(msgfile_name,sizeof(msgfile_name)-1,"%s",fname);
    else
      snprintf(msgfile_name,sizeof(msgfile_name)-1,"%s.%4.4d",fname,mpi_rank);

    if(binary_output) {
      fd_out=open(msgfile_name,O_CREAT|O_APPEND|O_WRONLY,0666);
    }else{
      file_out=fopen(msgfile_name,"a");
    }
  }
  if(fd_out < 0 && file_out == NULL) return;
  if(binary_output) {
    nbytes=write(fd_out,&table[0],sizeof(table_entry)*nentries);
  }else{
    for(i=1;i<nentries;i++) { char *c=table[i].tag ; fprintf(file_out,"%c%c%c%c%c%c%c%c: %20.20Ld\n",
                              c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],table[i].value); }
  }
  table_index=-1;
}
void write_history_trace()
{
  write_history_trace_();
}
void write_history_trace__()
{
  write_history_trace_();
}

#ifdef TEST
main()
{
  int i;
  float delta,denom;
  long long ldelta;

  printf("Beginning of test \n");
  log_use_tsc();
  for (i=0 ; i<=4001 ; i+=1) {  log_history_trace_("tototati"); }
  log_use_gettimeofday__();
  for (i=i ; i<=8191 ; i++) {  log_history_trace_("tagadayo"); }
  printf("--------( 1  )------\n");
  printf("--------(tsc)-------\n");
  for (i=0 ; i<7800 ; i+=200) {
     if(i>=4000 & i<4400)continue;
     if(i==4400) printf("--------(time)------\n");
     if(i==2000) printf("--------( 200)------\n");
     ldelta=table[i+200].value-table[i].value;
     denom=200.0;
     if(i<2000) {ldelta=table[i+1].value-table[i].value; denom=1.0;}
     delta=ldelta;
     delta /= denom;
     { char *c=table[i].tag ; printf("%c%c%c%c%c%c%c%c: %4.1f\n",
      c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],delta); }
  }
  printf("\n");
  write_history_trace();
}
#endif

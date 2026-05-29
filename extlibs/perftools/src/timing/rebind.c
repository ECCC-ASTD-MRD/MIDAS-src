#if defined(linux)
#define _GNU_SOURCE
#include <sched.h>
#include <sys/syscall.h>

//extern int Process_hostTaskCount2();
//extern int Process_hostTaskRank2();
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include "weak_openmp.h"

static int initialized=0;
static int hostTaskCount=1;    // tasks per host
static int hostTaskRank=0;     // rank of this task on host
static int taskRank=0;         // rank of this task in WORLD
static int nLCPU=-1;           // number of logical CPUs on host
static int shuffle_factor=1;   // shuffle factor, should be a power of 2 (1, 2 , maybe 4)
static int *list_of_cpus;      // translation index list for cpu number
static char *TARGET_CPU_LIST=NULL;       // contents of environment var TARGET_CPU_LIST
static int mapped_cpus=0;      // number of entries (if any) in environment var TARGET_CPU_LIST

/*
 * convert list of integers in string envar into an array of integers
 * max array size is dimliste. return 0 if sttring pointer is NULL
 * return number of integers stuffed into array
 */
static int strtoilist( char *envar, int *list, int dimliste)
{
  int i=0;
  if(envar == NULL) return 0;
  while(*envar != '\0' && i < dimliste)  // stop if array full or at end of string
  {
    if(*envar==' ') envar++;  // skip blanks
    else
    {
      list[i]=atoi(envar);
//      printf("list[%d]=%d\n",i,list[i]);
      i++;
      while(*envar != ' ' && *envar != '\0') envar++;  // skip blanks
    }
  }
  return i;
}

void set_tasks_per_node(int *n)
{
  hostTaskCount=*n;
}
void set_tasks_per_node_(int *n)  { set_tasks_per_node(n) ; }
void set_tasks_per_node__(int *n) { set_tasks_per_node(n) ; }

void set_rank_on_node(int *n)
{
  hostTaskRank=*n;
}
void set_rank_on_node_(int *n)  { set_rank_on_node(n); }
void set_rank_on_node__(int *n) { set_rank_on_node(n); }

static int compare (const void *e1, const void *e2)
{
  int v1 = *((int *) (e1));
  int v2 = *((int *) (e2));
  return ((v1 < v2 ) ? -1 : (v1 > v2 ) ? 1 : 0);
}

static int TargetCpu(int taskThreadCount, int taskThreadRank)
{
  int targetCPU=-1;
  int thread_stride;
  int task_stride;
  int task_offset;
  int excesscpufactor;

  if( nLCPU == -1) return -1;

  thread_stride = nLCPU / (hostTaskCount * taskThreadCount);
  task_offset = thread_stride;
  thread_stride *= shuffle_factor;
  task_stride = nLCPU / hostTaskCount;
  task_stride *= shuffle_factor;

  targetCPU = taskThreadRank*thread_stride + (hostTaskRank/shuffle_factor)*task_stride;
  if(thread_stride==0) {   /* OOPS, more threads that we have CPUs */
    excesscpufactor=(hostTaskCount * taskThreadCount + nLCPU - 1)/nLCPU;
    targetCPU = taskThreadRank/excesscpufactor + (hostTaskRank/shuffle_factor)*task_stride;
  }
  targetCPU = targetCPU  + (hostTaskRank%shuffle_factor) * task_offset;

  return targetCPU;
}

static void InitializeGeneral(void)  /* all architectures */
{
  char *shuffle_bind=NULL;
  int i;

/*
*   Get the shuffle binding flag
*/
  shuffle_bind=getenv("SHUFFLE_BIND");
  if( shuffle_bind != NULL ) {
    shuffle_factor=(atoi(shuffle_bind));
    if(shuffle_factor <= 1) shuffle_factor = 1;
    if(hostTaskCount%shuffle_factor != 0) {
      fprintf(stderr,"WARNING: cannot shuffle %d tasks by %d\n",hostTaskCount,shuffle_factor);
      shuffle_factor = 1;
    }
    fprintf(stderr,"INFO: shuffle binding order=%d\n",shuffle_factor);
  }
  nLCPU = sysconf(_SC_NPROCESSORS_ONLN);  /* nLCPU is static global */
  list_of_cpus = (int *)calloc( nLCPU , sizeof(int) );
  mapped_cpus = nLCPU;
  for(i=0;i<nLCPU;i++)list_of_cpus[i]=i;   // linear one to one mapping unless told otherwise
  fprintf(stderr,"There are %d active CPUs in node\n",nLCPU);
  TARGET_CPU_LIST = getenv("TARGET_CPU_LIST");
  if(TARGET_CPU_LIST != NULL) {
    mapped_cpus = strtoilist (TARGET_CPU_LIST,list_of_cpus,nLCPU);
    for(i=mapped_cpus;i<nLCPU;i++) list_of_cpus[i] = list_of_cpus[i%mapped_cpus] ; // fill to nLCPU entries
    }
}

#if defined(_AIX)

#include <omp.h>
#include <string.h>
#include <sys/processor.h>
#include <sys/systemcfg.h>
#include <sys/thread.h>
#include <sys/types.h>

static int InitializeSpecific(void)  /* AIX specific */
{
  char *taskRankString;
  char *hostTaskListString;
  char *originalHostTaskListString;
  char *unconv;
  int *taskRankList;
  char *token;
  int i;
/*
 *  Establish the PE rank of the present MPI task
 */

  taskRankString = getenv ("MP_CHILD");
  unconv = 0;
  if(taskRankString!=NULL)taskRank = strtol (taskRankString, &unconv, 0);  /* taskRank is static global */

/*
 *  Establish the number of PE tasks for this job on this host (hostTaskCount)
 */

  originalHostTaskListString = getenv ("MP_COMMON_TASKS");
  hostTaskListString = (char *) malloc ((strlen (originalHostTaskListString) + 1) * sizeof (char));
  strcpy (hostTaskListString, originalHostTaskListString);
  token = strtok (hostTaskListString, ":");
  hostTaskCount = 1 + (int) strtol (token, &unconv, 0);                      /* hostTaskCount is static global */

/*
 *  Establish the relative rank of this PE task on this node (hostTaskRank)
 */

  taskRankList = (int *) malloc ((size_t) (hostTaskCount * sizeof (int)));
  taskRankList[0] = taskRank;
  for (i = 1; i < hostTaskCount; ++i) {
    token = strtok (NULL, ":");
    taskRankList[i] = (int) strtol (token, &unconv, 0);
  }
  free (hostTaskListString);

  qsort (taskRankList, hostTaskCount, (size_t) (sizeof (int)), compare);

  for (i = 0; i < hostTaskCount; ++i)
    if (taskRankList[i] == taskRank) {
      hostTaskRank = i;                                                       /* hostTaskRank is static global */
      break;
    }
  free(taskRankList);
  if(hostTaskRank==0) printf("There %d MPI task(s) running in node\n",hostTaskCount);
}

/*
 * unbind and rebind the current thread
 */

int thread_rebind (void)
{
  int taskThreadCount;
  int taskThreadRank;
  int error=0;
  cpu_t targetCPU;
  tid_t tid;

  if(initialized == 0) /* do this only once */
  {
    InitializeSpecific();
    InitializeGeneral();
    initialized = 1;
  } /* end of do this only once */

/*
 *  Establish the number of threads for this task
 */
    taskThreadCount = OMP_GET_num_threads ();
/*
 *  Unbind at the thread level; if not bound, this will
 *  return with an error, which we will ignore
 */

    tid = thread_self ();
    bindprocessor (BINDTHREAD, tid, PROCESSOR_CLASS_ANY);

/*
 *  Establish the target CPU for this thread
 */

    taskThreadRank = OMP_GET_thread_num ();
    targetCPU = TargetCpu(taskThreadCount, taskThreadRank);
    if(hostTaskCount * taskThreadCount > nLCPU)
      printf("We have more tasks (%d) on node than we have CPUs (%d)\n",hostTaskCount * taskThreadCount,nLCPU);

/*
 *  Bind the thread to the logical processor
 */
//     targetCPU = list_of_cpus[targetCPU];   // remap from list of cpus
    if (bindprocessor (BINDTHREAD, tid, targetCPU) == -1) {
      error = errno;
    }
    printf("thread %2.2d bound to CPU %2.2d\n",OMP_GET_thread_num(),targetCPU);
    if (error != 0) {
      fprintf (stderr,
             "ERROR: thread_rebind: bindprocessor: %s\n",
             strerror (error));
      return -1;
    }

  return targetCPU;
}

int process_rebind(void)  /* unbind and rebind process */
{
  int targetCPU;
  pid_t pid;

/*
 *  Unbind at the process level: if not bound, this will
 *  return with an error, which we will ignore
 */
  pid = getpid ();
  bindprocessor (BINDPROCESS, pid, PROCESSOR_CLASS_ANY);  /* unbind master thread */

  thread_rebind();  /* rebind master thread */
  return 0;
}

void c_process_threads_rebind(void)  /* rebind all C OpenMP threads (AIX only for now) */
{
  int targetCPU;

#pragma omp parallel private (targetCPU)
  {
    targetCPU=thread_rebind();
    printf("thread %2.2d bound to CPU %2.2d\n",OMP_GET_thread_num(),targetCPU);
  }
}

#endif

#if defined(linux)
/* Linux systems */

static cpu_set_t Linux_CPU_ANY;

static void InitializeSpecific(void)  /* Linux specific */
{
  char *taskRankString;
  char *unconv;
  int i;
/*
 *  Establish the PE rank of the present task
 */
  taskRankString = getenv ("PMI_RANK");
  unconv = 0;
  if(taskRankString!=NULL)taskRank = strtol (taskRankString, &unconv, 0);  /* taskRank is static global */
  hostTaskCount = 1;
  hostTaskRank=0;
  CPU_ZERO((&Linux_CPU_ANY)) ;
  for(i=0 ; i<nLCPU ; i++) CPU_SET(i,&(Linux_CPU_ANY));
}

static pid_t mygettid(void)
{
  return syscall(__NR_gettid);
}

int thread_rebind(void) /* unbind and rebind thread */
{
  int taskThreadCount;
  int taskThreadRank;
  int targetCPU;
  pid_t tid;
  cpu_set_t mycpus;

  if(initialized == 0) /* do this only once */
  {
    InitializeSpecific();
    InitializeGeneral();
    initialized = 1;
  } /* end of do this only once */
/*
 *  Establish the number of threads for this task
 */
  taskThreadCount = OMP_GET_num_threads ();
  tid = mygettid();
//printf("tid=%d, pid=%d, nthreads=%d\n",tid,getpid(),taskThreadCount);

  sched_setaffinity(tid,sizeof(cpu_set_t),&Linux_CPU_ANY);  /* unbind thread */
/*
 *  Establish the target CPU for this thread
 */

  taskThreadRank = OMP_GET_thread_num ();
  targetCPU = TargetCpu(taskThreadCount, taskThreadRank);
  targetCPU = list_of_cpus[targetCPU];  // remap from list of cpus
  CPU_ZERO((&mycpus)) ;
  CPU_SET(targetCPU,&(mycpus));
  sched_setaffinity(tid,sizeof(cpu_set_t),&mycpus);   /* bin thread */
  printf("thread %2.2d now bound to CPU %2.2d\n",taskThreadRank,targetCPU);
  return(targetCPU);
}

int process_rebind(void) /* unbind and rebind process */
{
  pid_t pid;
  int thread_rebind();

  pid = getpid();
  sched_setaffinity(pid,sizeof(cpu_set_t),&Linux_CPU_ANY);  /* unbind process */

  thread_rebind();
  return 0;
}

void c_process_threads_rebind(void)  /* will rebind all C OpenMP threads when C OpenMP matches FORTRAN OpemMP */
{
  int targetCPU;

  {
    targetCPU=thread_rebind();
    printf("thread %2.2d bound to CPU %2.2d\n",OMP_GET_thread_num(),targetCPU);
  }
}

#endif

/* Fortran extra entries (AIX and Linux) */
#pragma weak thread_rebind_=thread_rebind
#pragma weak thread_rebind__=thread_rebind

#pragma weak process_rebind_=process_rebind
#pragma weak process_rebind__=process_rebind

#ifdef TEST
main()
{
process_rebind();
//printf("%d %d\n",Process_hostTaskCount2(),Process_hostTaskRank2());

}
#endif

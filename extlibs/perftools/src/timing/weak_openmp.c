#include <stdlib.h>

#if defined(linux)

#pragma weak omp_get_num_threads
int omp_get_num_threads();

int OMP_GET_num_threads() { if(omp_get_num_threads != NULL) return omp_get_num_threads(); return 1; }

#pragma weak omp_get_thread_num
int omp_get_thread_num();

int OMP_GET_thread_num()  { if(omp_get_thread_num != NULL) return omp_get_thread_num(); return 0; }

#endif
#if defined(_AIX)

int OMP_GET_num_threads() { return omp_get_num_threads(); }
int OMP_GET_thread_num()  { return omp_get_thread_num(); }

#endif

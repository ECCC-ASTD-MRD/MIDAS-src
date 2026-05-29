#if defined(_AIX)
#pragma weak _hpmInit_
void _hpmInit_(int a, char* b, int c) { printf("WARNING: hpm library deactivated, use -lhpm_r when linking\n"); DumhpmInit(a,b) ; }
#pragma weak _hpm_start_
void _hpm_start_(int a, int b, char *c, char *e, int f, int g) { DumhpmStart(a,e) ; }
#pragma weak _hpm_stop_
void _hpm_stop_(int a, int b) { DumhpmStop(a); }
#pragma weak _hpm_tstart_
void _hpm_tstart_(int a, int b, char *c, char *e, int f, int g) { DumhpmTstart(a,e) ;; }
#pragma weak _hpm_tstop_
void _hpm_tstop_(int a, int b) { DumhpmTstop(a) ; }
#pragma weak hpmTerminate
void hpmTerminate(int a) { DumhpmTerminate(a) ; }
#endif

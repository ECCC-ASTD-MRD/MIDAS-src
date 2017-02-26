#include <unistd.h> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define BUFFER_SIZE    32768
#define REBM_DONE_STRING "VAR3D_STATUS=REBM_DONE\n"
#define BYTES_TO_READ  23      /* length of string REBM_DONE_STRING */

// On compile avec la commande
//     gcc oavar_monitor_rebm.c -o oavar_monitor_rebm -Wall

int main(int argc, char **argv) {
  int fd;
  char buffer[BUFFER_SIZE];
  char rebm_done[BUFFER_SIZE];
  pid_t pp=getppid();

  if(argc-1 != 2) {
    fprintf(stderr,"argument count must be 2\n");
    exit(1);
  }

  if(fork()) exit(0);

  while(1) {

    if(kill(pp,0)) exit(0);

    if( (fd=open(argv[1],O_RDONLY )) >= 0 ) {
      ssize_t bytes_read;

      strcpy(rebm_done,"");
      while(strcmp(rebm_done,REBM_DONE_STRING)!=0) {
	bytes_read = read(fd, rebm_done, (size_t) BYTES_TO_READ);
	if (bytes_read!=BYTES_TO_READ) {
	  fprintf(stderr,"Could only read %d bytes out of %d bytes in file '%s'\n",(int) bytes_read,BYTES_TO_READ,argv[1]);
	  exit(1);
	}
	sleep(2);
	lseek(fd,(off_t) 0,SEEK_SET);
      }
      close(fd);

      printf("file=%s,cmd=%s\n",argv[1],argv[2]);
      snprintf(buffer,sizeof(buffer)-1,"%s %s",argv[2],argv[1]);
      printf("Executing:%s\n",buffer);
      system(buffer);

      exit(0);
    }
    sleep(2);
  }
}

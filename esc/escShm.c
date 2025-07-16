#include "esc_local.h"
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>


#if defined(__GNU_LIBRARY__) && !defined(_SEM_SEMUN_UNDEFINED)
/* union semun is defined by including <sys/sem.h> */
#else
#if defined(_SEM_SEMUN_UNDEFINED)
/* according to X/OPEN we have to define it ourselves */
union semun {
  int val;                    /* value for SETVAL */
  struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
  unsigned short int *array;  /* array for GETALL, SETALL */
  struct seminfo *__buf;      /* buffer for IPC_INFO */
};
#endif 
#endif

/* Number of elements in shared memory buffer */
#define NUM_ELEM 10

#define SHMEMLENGTH 65536
#define S0 0
#define S1 1

static long int BegOfRecord=0L;
int SemID;
int ShmID;
static struct sembuf bufV={S0, 1,IPC_NOWAIT};
static struct sembuf bufP={S1, -1,SEM_UNDO};
void *lShm;

int ESUpDateCleanFile()
{
  FILE *lF;

  lF	=fopen("escCLEAN","r+");
  if(lF == NULL){
    lF	=fopen("escCLEAN","w+");
    BegOfRecord=0L;
  }
  else{
    fseek(lF,(long int)(-12),SEEK_END);
    BegOfRecord=ftell(lF);
  }
  fprintf(lF,"ipcrm sem %d\n",SemID);
  fprintf(lF,"ipcrm shm %d\n",ShmID);
  fprintf(lF,"kill %d\n",getpid());
  fprintf(lF,"rm escCLEAN\n");
  fclose(lF);
  chmod("escCLEAN",00777);
  printf("escCLEAN - if exists, then use for cleaning system resources\n");
  return(0);
}

int ESRestoreCleanFile()
{
  int k;

  k	=open("escCLEAN",O_RDWR);
  if(k == -1){
    return(0);
  }
  if(BegOfRecord == 0L){
    close(k);
    system("rm escCLEAN");
    return(0);
  }
  lseek(k,BegOfRecord,SEEK_SET);
  write(k,"rm escCLEAN\n",12);
  ftruncate(k,BegOfRecord+12);
  close(k);
  return(0);
}

#ifndef ASTRA
int main(int argc, char **argv)
{
  int Fl=0;
  int Pid;
  int na;
  char ln[32];
  int i,k;
  char elem;
  union semun Mysemun;
  
  strcpy(ln,argv[1]);
  sscanf(argv[2],"%d",&Pid);
  sscanf(argv[3],"%d",&na);

  /* Init Shared memory */
  ShmID	=shmget((key_t)(Pid+1), 0, 0666|IPC_CREAT);

  /* Init semaphores */
  /* Two semaphores (Empty and Full) are created in one set */
  SemID	=semget((key_t)Pid, 0,0666|IPC_CREAT);
  /* attach shared memory to process */
  lShm=shmat(ShmID,NULL,0);
  
#ifdef H
  printf("Esc: SemID=%d %d %lx\n",SemID,ShmID,lShm);
#endif

  semop(SemID, &bufP, 1);   
  esc00_(argv[1],&na);
  semop(SemID, &bufV, 1);   		   

  ESUpDateCleanFile();
  while(1){
    semop(SemID, &bufP, 1);   
    Fl	=*(int*)lShm;
    switch(Fl){
    case 1:
      EscIn();
      break;
    case 2:
      EscEq();
      semop(SemID, &bufV, 1);   		   
      semop(SemID, &bufP, 1);   
      break;
    case 3:
      E2Astra();
      semop(SemID, &bufV, 1);   		   
      semop(SemID, &bufP, 1);   
      break;
    case 4:
      semop(SemID, &bufV, 1);   		   
      esc1_();
      ESRestoreCleanFile();
      exit(0);
      break;
    case 5:
      ESSendProfile();
      semop(SemID, &bufV, 1);   		   
      semop(SemID, &bufP, 1);   
      break;
    case 6:
      ESSendFourierR();
      semop(SemID, &bufV, 1);   		   
      semop(SemID, &bufP, 1);   
      break;
    case 7:
      ESSendFourierZ();
      semop(SemID, &bufV, 1);   		   
      semop(SemID, &bufP, 1);   
      break;
    default:
      break;
    }
    semop(SemID, &bufV, 1);   		   
  }

#ifdef H
  i=0;
  while(i < 26){
    /* Wait(Full) */
    semop(SemID, &bufP, 1);   
    for(k=0; k < 10 && i < 26; k++){
      /* Get and consume element, Example element are chars 'a'-'z' */
      elem	=*((char *)lPtr+(i%NUM_ELEM));
      printf("                             Consumed Elem '%c'\n", elem);
      i++;
    }
    /* Signal(Empty) */
    semop(SemID, &bufV, 1);   		   
  }
#endif
  exit(0);
}
#endif


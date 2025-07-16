#ifdef SHMEM

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>
#include <string.h>

#if defined(__GNU_LIBRARY__) && !defined(_SEM_SEMUN_UNDEFINED)
/* union semun is defined by including <sys/sem.h> */
#else
/* according to X/OPEN we have to define it ourselves */
union semun {
  int val;                    /* value for SETVAL */
  struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
  unsigned short int *array;  /* array for GETALL, SETALL */
  struct seminfo *__buf;      /* buffer for IPC_INFO */
};
#endif

#define SHMEMLENGTH 65536
#define Sem0 0
#define Sem1 1

static int ESNa=20,ESNa1,ESNp,ESNp1,ESFp,ESFp1;

static int Pid;
static int SemID;
static int ShmID;
static union semun Mysemun;
static struct sembuf bufP	={Sem0,-1,SEM_UNDO};
static struct sembuf bufV	={Sem1, 1,IPC_NOWAIT};
static void *lShm;
static int kD,kI;
static double c2gp = 6.2831853071796;

void esc00_(char *Machine, int *na)
{
  int i,k;
  char elem;
  char *lc,ln[64];

  kD	=sizeof(double);
  kI	=sizeof(int);

  /* Initialize */
  Pid	=getpid();

  /* Init semaphores */
  /* Two semaphores (Empty and Full) are created in one set */
  SemID	=semget((key_t)Pid, 2, 0666|IPC_CREAT);

  /* Init Empty to number of elements in shared memory */
  Mysemun.val	=1;
  semctl(SemID, Sem0, SETVAL, Mysemun);

  /* Init Full to zero, no elements are produced yet */
  Mysemun.val	=0;
  semctl(SemID, Sem1, SETVAL, Mysemun);  

  /* Init Shared memory */
  ShmID	=shmget((key_t)(Pid+1),SHMEMLENGTH, 0666 | IPC_CREAT);

  /* attach shared memory to process */
  lShm=shmat(ShmID,NULL,0);

  lc	=ln;
  strcpy(lc,"escShm ");
  lc	+=7;
  if(Machine != NULL){
    while(isspace(*Machine)){
      Machine++;
    }
    while(*Machine != '\0' && !isspace(*Machine)){
      *lc++	=*Machine++;
    }
  }
  else{
    strcpy(lc,"Machine");
    lc	+=7;
  }
  sprintf(lc," %d %d&",Pid,*na);
  semop(SemID, &bufP, 1);   
  semop(SemID, &bufV, 1);   		   

  system(ln);
  semop(SemID, &bufP, 1);   
  ESNa	=*(int*)lShm;
  ESNa1	=ESNa+1;
  ESNp	=*(int*)(lShm+kI);
  ESNp1	=ESNp+1;
  ESFp	=*(int*)(lShm+2*kI);
  ESFp1	=ESFp+1;
  semop(SemID, &bufV, 1);   		   
#ifdef H
  i	=0;
  while(i < 26){
    semop(SemID, &bufP, 1);   
    /* Wait(Empty) */
    for(k=0; k < 10 && i < 26; k++){
      /* Produce element, Example element are chars 'a'-'z' */
      elem='a'+i;
      printf("Produced elem '%c'\n", elem);
      
      /* Put element into shared memory buffer */
      *((char *)lShm+(i%10))=elem;
      i++;
    }    
    /* Signal(Full) */
    semop(SemID, &bufV, 1);   		   
  }
#endif
  return;
}

void esc1_()
{
  int k;
  struct shmid_ds Myshmid_ds;
  semop(SemID, &bufP, 1);   
  k	=4;
  memcpy(lShm,(void*)&k,kI);
  semop(SemID, &bufV, 1);   		   

  semop(SemID, &bufP, 1);   
  /* Remove semaphores */
  semctl(SemID,0, IPC_RMID, Mysemun);
  /* Remove shared memory */
  shmctl(ShmID, IPC_RMID, &Myshmid_ds);
  return;
}

void escin_(double	*pProf	/* p[] - pressure related profile */
	    ,double	*jProf	/* j[] - current related profile */
	    ,double	*aProf	/* a[] - normalized sqrt(Phi) square root 
				   of toroidal flux.
				   If a=NULL (dropped in FORTRAN) - uniform 
				   grid from 0 till 1,
				   otherwise it is considered as an array of 
				   grid points */
	    ,int	*nProf	/* number <= 129 of profile points including \
				   magnetic axis and plasma edge */
	    ,int	*hH	/* First digit h specifies the meaning
				   of pProf:
				   0 - j_p;
				   1 - P =dp/dpsi;
				   2 -p [MPa].
				   Second digit H specifies the meaning
				   of jProf:
				   0 - j_s [MA/m^2];
				   1- j_|| [MA/m^2]=(j dot B)/(R_0 B grad phi),
				                   R_0 = magnetic axis;
				   2- j_||R_0 [MA/m] = (j dot B)/(B grad phi);
				   3 - T=FF';
				   6 - q;
				   7 - 1/q;
				   For Example:
				   26 - p[] and q[] profiles are supplied;
				   21 - p[] and j||[] profiles are supplied;
				   0 -  jp[] and js[] profiles are supplied.

				   Possible combinations are limited to:
				   0,1,2,6,7
				   10,11,12,16,17
				   21,22,26,27
				 */
	    ,double	*Rpv	/* R[m]- plasma boundary */
	    ,double	*Zpv	/* Z[m]- plasma boundary */
	    ,int	*Npv	/* number <= 257 of the plasma-vacuum 
				   points.
				   If *Npv >12, the first and last points
				   coincide */
	    ,double	*RBtor	/* RBtor [m Tesla] outside the 
				   plasma */
	    ,double	*Rext	/* Reference major radius [m] */
	    )
{
  int k;
  int kP,kB;
  void *lc;

  kP	=(*nProf)*kD;
  kB	=(*Npv)*kD;

  semop(SemID, &bufP, 1);   
  lc	=lShm;
  k	=1;
  memcpy(lc,(void*)&k,kI);
  lc	+=kI;
  *(char*)lc++	=aProf == NULL ? 'n' : 'y';
  memcpy(lc,(void*)nProf,kI);
  lc	+=kI;
  if(aProf != NULL){
    memcpy(lc,(void*)aProf,kP);
    lc	+=kP;
  }
  memcpy(lc,(void*)pProf,kP);
  lc	+=kP;
  memcpy(lc,(void*)jProf,kP);
  lc	+=kP;
  memcpy(lc,(void*)hH,kI);
  lc	+=kI;
  memcpy(lc,(void*)Npv,kI);
  lc	+=kI;
  memcpy(lc,(void*)Rpv,kB);
  lc	+=kB;
  memcpy(lc,(void*)Zpv,kB);
  lc	+=kB;
  memcpy(lc,(void*)RBtor,kD);
  lc	+=kD;
  memcpy(lc,(void*)Rext,kD);
  lc	+=kD;
  semop(SemID, &bufV, 1);   		   
  return;
}

void esc_(
	  int *Ffail 		/* Flag of failure:
				   0 - normal operation;
				   1 - some problems;
				   2 - problems near the boundary;
				   4 - problems near the axis;
				   8 - problems in the middle;
				   ....;
				   */
	  ,int *Fjob 		/* Flag for job assignment:
				   0 - nothing special was requested;
				   1 - save the geometry;
				   2 - take the geometry;
				   4 - write ESI data files;
				   ...;
				   */ 
	  ,int *Mpol		/* Working number of Fourier harmonics
				   in $\Psi$ */
	  ,double *sTime	/* time */
	  )
{
  int k;
  void *lc;
  semop(SemID, &bufP, 1);   
  lc	=lShm;
  k	=2;
  memcpy(lc,(void*)&k,kI);
  lc	+=kI;
  memcpy(lc,(void*)Ffail,kI);
  lc	+=kI;
  memcpy(lc,(void*)Fjob,kI);
  lc	+=kI;
  memcpy(lc,(void*)Mpol,kI);
  lc	+=kI;
  memcpy(lc,(void*)sTime,kD);
  lc	+=kI;
  semop(SemID, &bufV, 1);   		   

  semop(SemID, &bufP, 1);   
  *Ffail=*(int*)(lShm+kI);
  k	=0;
  memcpy(lShm,(void*)&k,kI);
  semop(SemID, &bufV, 1);   		   
  return;
}

/* ESC to ASTRA output routine */
void e2astra_(double	*Rmaxis
	      ,double	*Zmaxis
	      ,double	*gFtor
	      ,double	*Fpol
	      ,double	*gM
	      ,double	*G22
	      ,double	*G33
	      ,double	*D2
	      ,double	*D3
	      ,double	*Sside
	      ,double	*Sside1
	      )
{
  int k,kA;
  void *lc;

  kA	=kD*ESNa1;
  semop(SemID, &bufP, 1);   
  lc	=lShm;
  k	=3;
  memcpy(lc,(void*)&k,kI);
  lc	+=kI;
  semop(SemID, &bufV, 1);   		   

  semop(SemID, &bufP, 1);   
  *Rmaxis	=*(double*)lc;
  lc		+=kD;
  *Zmaxis	=*(double*)lc;
  lc		+=kD;
  *gFtor	=*(double*)lc;
  lc		+=kD;
  memcpy((void*)Fpol,lc,kA);
  lc	+=kA;
  memcpy((void*)gM,lc,kA);
  lc	+=kA;
  memcpy((void*)G22,lc,kA);
  lc	+=kA;
  memcpy((void*)G33,lc,kA);
  lc	+=kA;
  memcpy((void*)D2,lc,kA);
  lc	+=kA;
  memcpy((void*)D3,lc,kA);
  lc	+=kA;
  memcpy((void*)Sside,lc,kA);
  lc	+=kA;
  memcpy((void*)Sside1,lc,kA);
  lc	+=kA;
  semop(SemID, &bufV, 1);   		   
  return;
}

/*-------------------------------Get functions -----------------------*/

static char *NmOfProf[6]={
  "bgY",
  "bgF",
  "aF",
  "gm",
  "q",
  "bsp"
};

#ifdef H
strcpy(PlPrP[0].Nm,	"j_PfShl");
strcpy(PlPrP[1].Nm,	"Ps     ");
strcpy(PlPrP[2].Nm,	"p      ");
strcpy(PlPrP[3].Nm,	"p'/a   ");
strcpy(PlPrP[4].Nm,	"ga_CHT ");
strcpy(PlPrP[5].Nm,	"n_e    ");
strcpy(PlPrP[6].Nm,	"n_i    ");
strcpy(PlPrP[7].Nm,	"T_e    ");
strcpy(PlPrP[8].Nm,	"T_i    ");
strcpy(PlPrP[9].Nm,	"p_e    ");
strcpy(PlPrP[10].Nm,	"p_i    ");
strcpy(PlPrP[11].Nm,	"gb_Sh  ");
strcpy(PlPrP[12].Nm,	"gb     ");


strcpy(PlCrP[0].Nm,	"j_GSh  ");
strcpy(PlCrP[1].Nm,	"j_||   ");
strcpy(PlCrP[2].Nm,	"J      ");
strcpy(PlCrP[3].Nm,	"FF'    ");
strcpy(PlCrP[4].Nm,	"<j.B>  ");
strcpy(PlCrP[5].Nm,	"Fs     ");
strcpy(PlCrP[6].Nm,	"q      ");
strcpy(PlCrP[7].Nm,	"1/q    ");
strcpy(PlCrP[8].Nm,	"gY     ");
strcpy(PlCrP[9].Nm,	"gF    ");
strcpy(PlCrP[10].Nm,	"F      ");
strcpy(PlCrP[11].Nm,	"FF     ");
strcpy(PlCrP[12].Nm,	"ir+il  ");
strcpy(PlCrP[13].Nm,	"ir-il  ");
strcpy(PlCrP[14].Nm,	"j_r    ");
strcpy(PlCrP[15].Nm,	"j_l    ");
strcpy(PlCrP[16].Nm,	"dgY/a  ");
strcpy(PlCrP[17].Nm,	"j_||R0 ");
#endif

/**
 * @file 
 * @brief List of Get methods. 
 * The following routines give access equilibrium profile data
 * and geometry after execution (Esc).
 */ 

 
/**
 * Get the radial number of grid points used internally by ESC. 
 * @param Na1 the grid size
 * @author L.E. Zakharov, A. Pletzer
 * @see escGetNp1
 * 
 */
void escGetNa1(int *Na1)
{
  *Na1 = ESNa1;
}
void escgetna1_(int *Na1)
{
  *Na1 = ESNa1;
}
void escgetna1(int *Na1)
{
  *Na1 = ESNa1;
}
void ESCGETNA1(int *Na1)
{
  *Na1 = ESNa1;
}
 
/**
 * Get the number of poloidal rays + 1 used internally by ESC. 
 * @param Np1 the number of poloidal sections + 1
 * @author L.E. Zakharov, A. Pletzer
 * @see escGetNa1
 */
void escGetNp1(int *Np1)
{
  *Np1 = ESNp1;
}
void escgetnp1_(int *Np1)
{
  *Np1 = ESNp1;
}
void escgetnp1(int *Np1)
{
  *Np1 = ESNp1;
}
void ESCGETNP1(int *Np1)
{
  *Np1 = ESNp1;
}
/**
 * Get the poloidal flux/(2*pi) [Wb/rad].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer
 * @see escGetPsi
 */

int ESgetProfile(int k,double *P, int n)
{
  int i;
  void *lc;
  double *ld;

  semop(SemID, &bufP, 1);   
  lc	=lShm;
  i	=5;
  memcpy(lc,(void*)&i,kI);
  lc	+=kI;
  ld	=(double*)lc;
  memcpy(lc,(void*)&k,kI);
  lc	+=kI;
  memcpy(lc,(void*)&n,kI);
  lc	+=kI;
  semop(SemID, &bufV, 1);   		   
  
  semop(SemID, &bufP, 1);   
  ld	=(double*)lc;
  for(i=0; i < n; i++){
    *P++	=*ld++;
  }
  semop(SemID, &bufV, 1);   		   
  return(0);
}

void escGetPsibar(double *prof, int *ns){
  int i;
  ESgetProfile(0,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	=-prof[i];
  }
}
void escgetpsibar_(double *prof, int *ns){
  int i;
  ESgetProfile(0,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	=-prof[i];
  }
}
void escgetpsibar(double *prof, int *ns){
  int i;
  ESgetProfile(0,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	=-prof[i];
  }
}
void ESCGETPSIBAR(double *prof, int *ns){
  int i;
  ESgetProfile(0,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	=-prof[i];
  }
}
 
/**
 * Get the poloidal flux [Wb].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer
 * @see escGetPsibar
 */
void escGetPsi(double *prof, int *ns){
  int i;
  ESgetProfile(0,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	*=-c2gp;
  }
}
void escgetpsi_(double *prof, int *ns){
  int i;
  ESgetProfile(0,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	*=c2gp;
  }
}
void escgetpsi(double *prof, int *ns){
  int i;
  ESgetProfile(0,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	*=c2gp;
  }
}
void ESCGETPSI(double *prof, int *ns){
  int i;
  ESgetProfile(0,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	*=c2gp;
  }
}

/**
 * Get the toroidal flux/(2*pi) [Wb/rad].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer
 * @see escGetPhi
 */
void escGetPhibar(double *prof, int *ns){
  ESgetProfile(1,prof,*ns);
}
void escgetphibar_(double *prof, int *ns){
  ESgetProfile(1,prof,*ns);
}
void escgetphibar(double *prof, int *ns){
  ESgetProfile(1,prof,*ns);
}
void ESCGETPHIBAR(double *prof, int *ns){
  ESgetProfile(1,prof,*ns);
}
 
/**
 * Get the toroidal flux [Wb].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetPhibar
 */
void escGetPhi(double *prof, int *ns){
  int i;
  ESgetProfile(1,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	*=c2gp;
  }
}
void escgetphi_(double *prof, int *ns){
  int i;
  ESgetProfile(1,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	*=c2gp;
  }
}
void escgetphi(double *prof, int *ns){
  int i;
  ESgetProfile(1,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	*=c2gp;
  }
}
void ESCGETPHI(double *prof, int *ns){
  int i;
  ESgetProfile(1,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	*=c2gp;
  }
}
/**
 * Get the covariant toroidal magnetic field function [Tm].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetIota escGetQ escGetP
 */
void escGetG(double *prof, int *ns){
  ESgetProfile(2,prof,*ns);
}
void escgetg_(double *prof, int *ns){
  ESgetProfile(2,prof,*ns);
}
void escgetg(double *prof, int *ns){
  ESgetProfile(2,prof,*ns);
}
void ESCGETG(double *prof, int *ns){
  ESgetProfile(2,prof,*ns);
}

/**
 * Get the 1/q profile [-].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetQ escGetG escGetP
 */ 
void escGetIota(double *prof, int *ns){
  ESgetProfile(3,prof,*ns);
}
void escgetiota_(double *prof, int *ns){
  ESgetProfile(3,prof,*ns);
}
void escgetiota(double *prof, int *ns){
  ESgetProfile(3,prof,*ns);
}
void ESCGETIOTA(double *prof, int *ns){
  ESgetProfile(3,prof,*ns);
}

/**
 * Get the safety factor (q) profile [-].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetIota escGetG escGetP
 */ 

void escGetQ(double *prof, int *ns){
  int i;
  ESgetProfile(3,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	=1./prof[i];
  }
}
void escgetq_(double *prof, int *ns){
  int i;
  ESgetProfile(3,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	=1./prof[i];
  }
}
void escgetq(double *prof, int *ns){
  int i;
  ESgetProfile(3,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	=1./prof[i];
  }
}
void ESCGETQ(double *prof, int *ns){
  int i;
  ESgetProfile(3,prof,*ns);
  for(i=0; i < *ns; i++){
    prof[i]	=1./prof[i];
  }
}

/**
 * Get the pressure in [mu0 Pa].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetIota escGetG escGetQ
 */
void escGetP(double *prof, int *ns){
  ESgetProfile(5,prof,*ns);
}
void escgetp_(double *prof, int *ns){
  ESgetProfile(5,prof,*ns);
}
void escgetp(double *prof, int *ns){
  ESgetProfile(5,prof,*ns);
}
void ESCGETP(double *prof, int *ns){
  ESgetProfile(5,prof,*ns);
}
 
/** 
 * Get the Fourier coefficients of the R coordinate [m].
 * @param rcos the returned cosine coefficients of sizes [ns][nf]
 * @param rsin the returned sine coefficients
 * @param nf the number of Fourier coefficients (input)
 * @param ns the radial size (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetFourierZ
 */

void escGetFourierR(double *rcos, double *rsin, int *nf, int *ns)
{
  int I,ic,inc;
  int i,j;
  void *lc;
  double *ld;

  semop(SemID, &bufP, 1);   
  lc	=lShm;
  i	=6;
  memcpy(lc,(void*)&i,kI);
  lc	+=kI;
  ld	=(double*)lc;
  memcpy(lc,(void*)nf,kI);
  lc	+=kI;
  memcpy(lc,(void*)ns,kI);
  lc	+=kI;
  semop(SemID, &bufV, 1);   		   
  
  semop(SemID, &bufP, 1);   
  inc	=2*(*nf)*kD;
  I	=SHMEMLENGTH-inc;
  ic	=0;
  for(i=0; i < *ns; i++){
    for(j=0; j < *ns; j++){
      *rcos++	=*ld++;
      *rsin++	=*ld++;
    }
    ic	+=inc;
    if(ic > I){
      ld	=(double*)(lShm+kI);
      ic	=0;
      semop(SemID, &bufV, 1);   		   
      semop(SemID, &bufP, 1);   
    }
  }
  semop(SemID, &bufV, 1);   		   
  return;
}
void escgetfourierr_(double *rcos, double *rsin, int *nf, int *ns){
  escGetFourierR(rcos, rsin, nf, ns);
}
void escgetfourierr(double *rcos, double *rsin, int *nf, int *ns){
  escGetFourierR(rcos, rsin, nf, ns);
}
void ESCGETFOURIERR(double *rcos, double *rsin, int *nf, int *ns){
  escGetFourierR(rcos, rsin, nf, ns);
}

/** 
 * Get the R coordinates [m]. The poloidal index varies faster.
 * @param r the returned array of size [ns][nt1].
 * @param nt1 the number of poloidal sections + 1 (input)
 * @param ns the radial size (input)
 * @param clockwise poloidal angle orientation (1 for clockwise, 
 *        0 for counterclockwise)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetZ
 */
void escGetR(double *r, int *nt1, int *ns, int *clockwise){
  double *rcos, *rsin;
  int NT1, NS, i, j, k, ni, nf;
  double sum, theta;
  NT1 = *nt1;
  NS  = *ns;
  nf = ESFp1;

  rcos =(double*) malloc( nf*NS*sizeof(double));
  rsin =(double*) malloc( nf*NS*sizeof(double));

  escGetFourierR(rcos, rsin, &nf, ns);

  if(*clockwise!=1){
    for(i=0; i<NS; i++){
      for(j=0; j<NT1; j++){
	sum = 0.0;
	theta = c2gp*(double)j/(double)(NT1-1);
	for(k=0; k<nf; k++){
	  ni = i*nf + k;
	  sum += rcos[ni]*cos(k*theta) + rsin[ni]*sin(k*theta);
	}
	r[i*NT1 + j] = sum;
      }
    }
  }else{
    for(i=0; i<NS; i++){
      for(j=0; j<NT1; j++){
	sum = 0.0;
	theta = c2gp*(double)(NT1-1-j)/(double)(NT1-1);
	for(k=0; k<nf; k++){
	  ni = i*nf + k;
	  sum += rcos[ni]*cos(k*theta) + rsin[ni]*sin(k*theta);
	}
	r[i*NT1 + j] = sum;
      }
    }
  }
    

  free(rcos);
  free(rsin);
}
void escgetr_(double *r, int *nt1, int *ns, int *clockwise)
{
  escGetR(r, nt1, ns, clockwise);
}
void escgetr(double *r, int *nt1, int *ns, int *clockwise)
{
  escGetR(r, nt1, ns, clockwise);
}
void ESCGETR(double *r, int *nt1, int *ns, int *clockwise)
{
  escGetR(r, nt1, ns, clockwise);
}
 
/** 
 * Get the Fourier coefficients of the Z coordinate [m].
 * @param zcos the returned cosine coefficients of sizes [ns][nf].
 * @param zsin the returned sine coefficients
 * @param nf the number of Fourier coefficients (input)
 * @param ns the radial size (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetFourierR
 */
void escGetFourierZ(double *rcos, double *rsin, int *nf, int *ns)
{
  int i,j;
  void *lc;
  double *ld;

  semop(SemID, &bufP, 1);   
  lc	=lShm;
  i	=7;
  memcpy(lc,(void*)&i,kI);
  lc	+=kI;
  ld	=(double*)lc;
  memcpy(lc,(void*)nf,kI);
  lc	+=kI;
  memcpy(lc,(void*)ns,kI);
  lc	+=kI;
  semop(SemID, &bufV, 1);   		   
  
  semop(SemID, &bufP, 1);   
  for(i=0; i < *ns; i++){
    for(j=0; j < 2; j++){
      *rcos++	=*ld++;
      *rsin++	=*ld++;
    }
    for(j=2; j < *nf; j++){
      *rcos++	=0.;
      *rsin++	=0.;
    }
  }
  semop(SemID, &bufV, 1);   		   
  return;
}

void escgetfourierz_(double *zcos, double *zsin, int *nf, int *ns){
  escGetFourierZ(zcos, zsin, nf, ns);
}
void escgetfourierz(double *zcos, double *zsin, int *nf, int *ns){
  escGetFourierZ(zcos, zsin, nf, ns);
}
void ESCGETFOURIERZ(double *zcos, double *zsin, int *nf, int *ns){
  escGetFourierZ(zcos, zsin, nf, ns);
}

/** 
 * Get the Z coordinate [m]. The poloidal index varies faster.
 * @param r the returned array of size [ns][nt1].
 * @param nt1 the number of poloidal sections + 1 (input)
 * @param ns the radial size (input)
 * @param clockwise poloidal angle orientation (1 for clockwise, 
 *        0 for counterclockwise)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetR
 */
void escGetZ(double *z, int *nt1, int *ns, int *clockwise){
  double *zcos, *zsin;
  int NT1, NS, i, j, k, ni, nf;
  double sum, theta;
  NT1 = *nt1;
  NS  = *ns;
  nf = ESFp1;

  zcos =(double*) malloc( nf*NS*sizeof(double));
  zsin =(double*) malloc( nf*NS*sizeof(double));

  escGetFourierZ(zcos, zsin, &nf, ns);

  if(*clockwise!=1){
    for(i=0; i<NS; i++){
      for(j=0; j<NT1; j++){
	sum = 0.0;
	theta = c2gp*(double)j/(double)(NT1-1);
	for(k=0; k<nf; k++){
	  ni = i*nf + k;
	  sum += zcos[ni]*cos(k*theta) + zsin[ni]*sin(k*theta);
	}
	z[i*NT1 + j] = sum;
      }
    }
  }else{
    for(i=0; i<NS; i++){
      for(j=0; j<NT1; j++){
	sum = 0.0;
	theta = c2gp*(double)(NT1-1-j)/(double)(NT1-1);
	for(k=0; k<nf; k++){
	  ni = i*nf + k;
	  sum += zcos[ni]*cos(k*theta) + zsin[ni]*sin(k*theta);
	}
	z[i*NT1 + j] = sum;
      }
    }
  }    
  free(zcos);
  free(zsin);
}
void escgetz_(double *z, int *nt1, int *ns, int * clockwise){
  escGetZ(z, nt1, ns, clockwise);
}
void escgetz(double *z, int *nt1, int *ns, int * clockwise){
  escGetZ(z, nt1, ns, clockwise);
}
void ESCGETZ(double *z, int *nt1, int *ns, int * clockwise){
  escGetZ(z, nt1, ns, clockwise);
}
#endif


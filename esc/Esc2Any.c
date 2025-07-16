/**
 * @file 
 * @brief ESC API. 
 *
 * This file contains a list of procedures to (in this order):
 *      - initialize ESC
 *      - set the input geometry and profiles
 *      - execute ESC
 *      - de-initialize
 *
 * To enable linking with Fortran on UNIX workstations, each API routine 
 * is in addition shadowed by equivalent routines bearing the same name 
 * but in lower case with and without underescore appended. On Crays, 
 * entry points in upper case are also included.
 *
 * @author L. E. Zakharov, A. Pletzer
 */

#include "esc_local.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifdef XWIN
#include <time.h>
#include <sys/time.h>

extern time_t time(time_t *t);
extern struct tm *localtime(const time_t *timep);
static struct tm *ltm;
static time_t stTime;
#endif

#ifdef SHMEM
#include <unistd.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>

#define SHMEMLENGTH 65536
#define S0 0
#define S1 1
extern int SemID;
extern int ShmID;
extern void *lShm;
static int kI,kD;
static struct sembuf bufV={S0, 1,IPC_NOWAIT};
static struct sembuf bufP={S1, -1,SEM_UNDO};
#endif

#define NDAT 21    /* Number ( <= 257) of data points 
		      taken from Transport code */

#ifndef mzl_ESC
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern int ESNa,ESNa1,ESNp,ESNp1,ESFp,ESFp1,TFCNt,TFCNp,ESLt,ESNt
,ESMp,ESMp1,ES2Mp,ES2Mp1,ESnMp,ESiAx;

extern int ESEqSolvFl,ESEqSolvIt,ESEqSolvRadC,ESEqSolvInPr,ESEqSolvInCr;
extern int ESFail;
extern double ESxDtPr[],ESxDtCr[],*ESaR0,*ESaZ0,*ESsa,*ESpa;
extern int ESnDtPr,ESnDtCr;
extern int ECEqIterFl; 

extern double ESBt,ESRBt,ESRext;
extern double *ESaF,*ESaF1a,*ESaF2a
,*ESgF,*ESgF1a,*ESgF2a
,*ESgm,*ESgm1a,*ESgm2a
,*ESFF,*ESFF1a,*ESFF2a;

extern double *ESgt,*ESsb,*ESsb1a,*ESsb2a
,*rcT,*rcT1a,*rcT2a,*rsT,*rsT1a,*rsT2a;
extern double *EScs1,*ESsn1;

extern double ESxDtPr[],ESxDtCr[];

extern double *R0T,*Z0T,*rcT,*rsT,*bT;

extern double *EZvr,*EZvz;
extern double *ESg22c,*ESg22c1,*ESg22c2,
*ESDPc,*ESDPc1,*ESDPc2,
*ESLc,*ESLc1,*ESLc2,
*ESVc,*ESVc1,*ESVc2;
extern char ESMachineNm[],ESRunID[];
extern double ESTime; 
extern int FlagEFIT;
extern int ESNMSE;
extern double ESrMSE[],ESgiMSE[];

static FILE *fp=NULL;
extern double EcRjet[],EcZjet[];
extern int EcNjet;

#endif

#ifndef mzl_ESC

extern double *ESRC,*ESRC1a,*ESRC2a
,*ESqgF,*ESqgF1a,*ESqgF2a
,*ESqgY,*ESqgY1a,*ESqgY2a
,*ESsh,*ESsh1a,*ESsh2a
,*ESqV,*ESqV1a,*ESqV2a

,*ESjp,*ESjp1a,*ESjp2a
,*ESsp,*ESsp1a,*ESsp2a
,*ESPs,*ESPs1a,*ESPs2a
,*ESgb,*ESgb1a,*ESgb2a
,*ESgB,*ESgB1a,*ESgB2a

,*ESjs,*ESjs1a,*ESjs2a
,*ESjb,*ESjb1a,*ESjb2a
,*ESaJ,*ESaJ1a,*ESaJ2a
,*ESaT,*ESaT1a,*ESaT2a
,*ESjB,*ESjB1a,*ESjB2a

,*ESsq,*ESsq1a,*ESsq2a
,*ESgm,*ESgm1a,*ESgm2a
,*ESgY,*ESgY1a,*ESgY2a
,*ESdgY,*ESdgY1a,*ESdgY2a
,*ESgF,*ESgF1a,*ESgF2a

,*ESaF,*ESaF1a,*ESaF2a
,*ESFF,*ESFF1a,*ESFF2a;

extern double *ESsr,*ESsra,*ESsrt,*ESsrat;
extern double *ESsz,*ESsza,*ESszt,*ESszat;
extern double *ESaB,*ESaBa,*ESaBt,*ESaBat;
extern double *ESgH,*ESgHa,*ESgHt,*ESgHat;

extern int ESXNa,ESXNa1;
extern double *ESXsx,*ESXsr,*ESXsj,*ESXsp,*ESXsq,*ESXsi;
extern int FlagPoints;

#ifndef amzl_ESCTR
static double *rdPV,*zdPV;
static int np=0,nj=0,npv=0,mpv=0;

static double cssp,cssj,cssr;
static double *ydPr,*ydCr,*xdPr,*xdCr;
static int *kdPr,*kdCr;
#endif

static int FlStart=1;

static double *lProf[8],*lProf2[8];
static char *NmOfProf[6]={
  "bgY",
  "bgF",
  "aF",
  "gm",
  "q",
  "bsp"
};

#ifndef amzl_ESC
#ifndef ablb_esc
/* ESC initiation, calling and deinitiation routines.
   Interior of routines is subject to changes by
   L.E. Zakharov without notice */

static int Flesc0=0;

int ESCInitStr()
{
  Flesc0=1;
  ESiAx	=0;
  ESNp	=64;
  ESMp	=12;
  ESMp	=6;
  ESFp	=2*ESMp;
  TFCNt	=0;
  TFCNp	=1;
  ESLt	=0;
  ESNt	=0;

  InitConstants();
  ESInitPlGeomInt();
  ESReInit1DPlV();
  ESReInit2DTPlV();
  ESReInitPlProf(21);
  EcInitArrays();

  ESReInitMetrTnsrCS();
  ECReInitFgY();
  
  ESGetPlVDAddr(&rdPV,&zdPV);
  ESInit1DPlV();
  ESInit2DTPlV();

  cssp	=EZcgm0;
  cssj	=EZcgm0;
  cssr	=1.;
#ifdef SHMEM
  kI	=sizeof(int);
  kD	=sizeof(double);
#endif
  lProf[0]	=ESgY;
  lProf[1]	=ESgF;
  lProf[2]	=ESaF;
  lProf[3]	=ESgm;
  lProf[4]	=ESsq;
  lProf[5]	=ESsp;

  lProf2[0]	=ESgY2a;
  lProf2[1]	=ESgF2a;
  lProf2[2]	=ESaF2a;
  lProf2[3]	=ESgm2a;
  lProf2[4]	=ESsq2a;
  lProf2[5]	=ESsp2a;
  return(0);
}

/** 
 * Another way of initializing Esc. Call this routine instead of 
 * esc0_ if you want to control the internal number of radial zones 
 * in Esc (which is set to 20 when calling esc0_).
 *
 * @param MachineRunID a string of the form '<Machine>.<Runid>' with '.'
 * as separator (input). For instance MachineRunID = "NSTX.123456Z12".
 *
 * @param na no. of radial intervals <=128 (input)
 *
 * @author L. E. Zakharov
 * @see esc0_, esc1_
 *
 */

void esc00_(char *MachineRunID, int *na)
{
  //#include <fenv.h>
#ifdef H
  feenableexcept(FE_INVALID|FE_OVERFLOW|FE_UNDERFLOW|FE_DIVBYZERO);
#endif

  if(Flesc0) return;
  if(MachineRunID != NULL){
    char *p0,*p1;
    int L;

    char *p2;
    p2	=MachineRunID;
    while(isspace(*MachineRunID)){
      MachineRunID++;
    }
    p0	=ESMachineNm;
    p1	=MachineRunID;
    while(*p1 != '.' && *p1 != '\0' && !isspace(*p1)){
      *p0++	=*p1++;
    }
    if(p0-ESMachineNm){
      *p0	='\0';
    }
    while(*p1 != '.' && *p1 != '\0'){
      p1++;
    }
    if(*p1 == '.'){
      p1++;
      while(isspace(*p1)){
	p1++;
      }
      if(*p1 != '\0'){
	CbStrCpy(ESRunID,p1);
      }
    }
#ifdef H
    printf("|%s| |%s|.|%s|\n",MachineRunID,ESMachineNm,ESRunID);
#endif
  }
  if(*na > 128){
    printf("ESC::WARNING too many radial zones %d > 128\n",*na);
    ESNa=128;
  }
  else{
    ESNa	=*na;
  }
  ESCInitStr();
#ifdef SHMEM
  memcpy(lShm,(void*)&ESNa,kI);
  memcpy(lShm+kI,(void*)&ESNp,kI);
  memcpy(lShm+2*kI,(void*)&ESFp,kI);
#endif
  return;
}

#ifdef SHMEM
int ESSendProfile()
{
  int i,k,n;
  double dx;
  double *ld,*lp,*lp2;
  int *li;

  li	=(int*)(lShm+kI);
  ld	=(double*)li;
  k	=*li++;
  n	=*li++;

  lp	=lProf[k];
  lp2	=lProf2[k];

  dx 	=ESsa[ESNa]/(n-1);
  for(i=0; i < n; i++){
    ESSetSplA(dx*i);
    splRA(ld,NULL,lp,lp2);
    ld++;
  }
  return(0);
}

int ESSendFourierR()
{
  int I,ic,inc;
  int i,j,k,n;
  int nfmax,mi;
  double dx;
  double *ld;
  int *li;

  li	=(int*)(lShm+kI);
  ld	=(double*)li;

  k	=*li++;
  n	=*li++;

  nfmax =k > ESFp1 ? ESFp1 : k;

  dx =ESsa[ESNa]/(n-1);

  inc	=2*k*kD;
  I	=SHMEMLENGTH-inc;
  ic	=0;
  for(i=0; i < n; i++){
    ESSetSplA(dx*i);
    mi  =0;
    splRA(ld,NULL,rcT,rcT2a);
    *ld++	+=ESaR0[0];
    *ld++	 =0.0;
    for(j=1; j < nfmax; j++){
      mi	+=ESNa1;
      splRA(ld,NULL,rcT+mi, rcT2a+mi);
      *ld++ *=2.0;
      splRA(ld,NULL,rsT+mi, rsT2a+mi);
      *ld++ *=2.0;
    }
    for(j=nfmax; j < k; j++){
      *ld++	=0.0;
      *ld++	=0.0;
    }
    ic	+=inc;
    if(ic > I){
      ld	=(double*)(lShm+kI);
      ic	=0;
      semop(SemID, &bufV, 1);   		   
      semop(SemID, &bufP, 1);   
    }
  }
}

int ESSendFourierZ()
{
  int i,j,k,n;
  double dx;
  double *ld;
  int *li;

  li	=(int*)(lShm+kI);
  ld	=(double*)li;

  k	=*li++;
  n	=*li++;
  dx =ESsa[ESNa]/(n-1);
  for(i=0; i < n; i++){
    ESSetSplA(dx*i);
    splRA(ld,NULL,rsT,rsT2a);
    *ld++	+=ESaZ0[0];
    splRA(ld,NULL,ESsb,ESsb2a);
    *ld++;
  }
}
#endif

void esc00(char *MachineRunID, int *na)
{ esc00_(MachineRunID, na); }

/** 
 * Initialization routine. This procedure should be called prior 
 * to any other ESC call, or after complete de-initialization.
 * 
 * @param MachineRunID a string of the form '<Machine>.<Runid>' with '.' 
 * as separator (input). For instance MachineRunID = "NSTX.123456Z12".
 *
 * @author L. E. Zakharov
 * @see esc1_
 */

void esc0_(char *MachineRunID)	/* ESC initiation routine */
{
  int na;
  na	=20;
  esc00_(MachineRunID,&na);
  return;
}

void esc0x_(char *MachineNm)	/* ESC initiation routine */
{
  if(Flesc0){
    return;
  }
  Flesc0=1;
  ESiAx	=1;
  ESNa	=19;
  ESNp	=64;
  ESMp	=12;
  ESFp	=2*ESMp;
  TFCNt	=0;
  TFCNp	=1;
  ESLt	=0;
  ESNt	=0;

  InitConstants();
  ESInitPlGeomInt();
  ESReInit1DPlV();
  ESReInit2DTPlV();
  ESReInitPlProf(NDAT);
  EcInitArrays();

  ESReInitMetrTnsrCS();
  ECReInitFgY();
  
  ESGetPlVDAddr(&rdPV,&zdPV);
  ESInit1DPlV();
  ESInit2DTPlV();
}

/**
 * Execution of ESC. This procedure can be called once ESC has been 
 * initialized (esc0_ or esc00_) and the geometry and profiles have 
 * been specified (using escin_, Esc_p_q or Esc_p_jparallel). 
 * Alternatively, instead of, or in addition to specifying the profiles, 
 * a previously saved geometry can be loaded by setting @a Fjob to 2 or 3.
 * 
 * @param Fjob job assignment (input)
 *        -  0 nothing special was requested
 *        -  1 save the geometry in file <MachineRunID>.es file
 *        -  2 load the geometry & profiles from file <MachineRunID>.es
 *        -  3 load from and then save the geometry & profiles in file 
 *             <MachineRunID>.es
 * @param Mpol number of poloidal Fourier modes (<=12) (input). 
 *             Typically @a Mpol should be 3 or 4.
 *
 * @param sTime used when Fjob != 0 to tag the equilibrium data in dump file 
 * <MachineRunID>.es
 *
 * @param Ffail: completion flag (output).
 *       - 0 sucess
 *       - 1 some problems
 *       - 2 problems near the boundary
 *       - 4 problems near the axis
 *       - 8 problems in the middle
 *       - 3, 11, 15 are usually mild problems associated with localized 
 *         inaccuracies.
 *
 * @see Esc_p_q, Esc_p_jparallel, escin_
 */

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
				   4 - Ballooning test;
				   8 - write ESI data files;
				   ...;
				   */ 
	  ,int *Mpol		/* Working number of Fourier harmonics
				   in $\Psi$ */
	  ,double *sTime	/* time */
	  )
{
  int i;
  if(2*(*Mpol) < ESFp){
    ESMp	=*Mpol;
    ESMp1	=ESMp+1;
    ES2Mp	=2*ESMp;
    ES2Mp1	=ES2Mp+1;
    ESnMp	=ESNa1*ES2Mp1;
  }
  else{
    printf("ESC::WARNING too many Fourier modes %d >= %d \n", *Mpol,ESFp/2);
  }
  ESTime	=*sTime; 
#ifndef mzl_ESC
#ifndef blb_Equil
#ifndef mzl_ESC
  if(ESiAx){
    ESReal2RefDataX(0);
    ESPlVdata2PlVX();
    ESGapCoreSepX();
    ESCoreX();
  }
  else{
    ESReal2RefData(0);
    ESReal2RefGeom(0);
    if(FlagPoints){
      ESPlVPoints2PlV(EcRjet,EcZjet,EcNjet);
    }
    else{
      ESPlVdata2PlV();
    }
    ESRef2RealGeom(0);
    if((*Fjob&0x02)){
      ECRestoreESC(ESTime);
    }
  }
  ESRealGeom(0);

#ifndef stg_RadC
#ifdef XWIN
  ESInpRadC();
#endif
#endif
#ifndef stg_PlVac
#ifdef XWIN
  ESInputPlVac();
#endif
#endif
#ifndef stg_EqSolv
#ifdef XWIN
  ESInput2DEqSolv();
#endif
#endif
  ESReInitMetrTnsrCS();
  ECRcs3D2Rcs2D(0);
  ES3DMetrTnsrCS();
  if((*Fjob&0x02)){
    *Ffail	=0;
    if((ESFail&0x0000000F))
      *Ffail	+=1;
    if((ESFail&0x000000F0) || (ESFail&0x0F000000))
      *Ffail	+=8;
    if((ESFail&0x00F00000) || (ESFail&0x0000F000))
      *Ffail	+=4;
    if((ESFail&0x000F0000) || (ESFail&0x00000F00))
      *Ffail	+=2;
    return;
  }
  ECReInitFgY();
  ESInpPlPrCr();
  EC1DEqSolv();
  ESPlProf2Data();
#ifndef stg_PlPr
#ifdef XWIN
  ESInpPlPr();
#endif
#endif
#ifndef stg_PlCr
#ifdef XWIN
  ESInpPlCr();
#endif
#endif
  ECRescue2DGeom(0);
  do{
    if(EC2DEqSolv(0) == 0){
      ECgYcs2dRZcsBound(0);
      ECVirtualPlV2vdrPlV();
      ECgYrc2dRZcs();
      ECvdrz2vrz();
      ECCheckJacob(0);
    }
    else{
      ECRestore2DGeom();
    }
    ECRcs2D2Rcs3D(0);
    if(ESiAx != 0){
      EcMetrTnsrX();
    }
    ESRealGeom(0);
    if(ES3DMetrTnsrCS()){
      ESFail	=0x0FFFFFFF;
      break;
    }
    EC1DEqSolv();
  }while(ECIfEqFound());
  *Ffail	=0;
  if((ESFail&0x0000000F))
    *Ffail	+=1;
  if((ESFail&0x000000F0) || (ESFail&0x0F000000))
    *Ffail	+=8;
  if((ESFail&0x00F00000) || (ESFail&0x0000F000))
    *Ffail	+=4;
  if((ESFail&0x000F0000) || (ESFail&0x00000F00))
    *Ffail	+=2;
  if(*Ffail){
    ECRestore2DGeom();
    ECRcs2D2Rcs3D(0);
    if(ESiAx != 0){
      EcMetrTnsrX();
    }
    ESRealGeom(0);
    ES3DMetrTnsrCS();
    EC1DEqSolv();
    ECEqIterFl	=0;
  }
  else{
    ECEqIterFl	=1;
  }
  ECRefGeom(0);
  if((*Fjob&0x01)){
    ECSaveESC();
  }
  if((*Fjob&0x08)){
    ESBasicFunctions();
#ifdef H
    ESWriteESI0();
#endif
    ESWriteESI(NULL,ESsr,ESsz,ESaB,ESgH,ESNa,ESNp);
  }
#ifndef stg_InitPFC
#endif
#ifndef stg_EqOut
#endif

#ifndef stg_BalSt
  if((*Fjob&4)){
    ESBalStab();
  }
#endif
#endif
#endif

#endif
}

#ifdef SHMEM
void EscEq()
{
  int i;
  int Ffail,Fjob,Mpol;
  double sTime;
  void *lc;

  int k;
  lc	=lShm+kI;
  Ffail	=*(int*)lc;
  lc	+=kI;
  Fjob	=*(int*)lc;
  lc	+=kI;
  Mpol	=*(int*)lc;
  lc	+=kI;
  sTime	=*(double*)lc;

  if(2*(Mpol) < ESFp){
    ESMp	=Mpol;
    ESMp1	=ESMp+1;
    ES2Mp	=2*ESMp;
    ES2Mp1	=ES2Mp+1;
    ESnMp	=ESNa1*ES2Mp1;
  }
  else{
    printf("ESC::WARNING too many Fourier modes %d >= %d \n",Mpol,ESFp/2);
  }
  ESTime	=sTime; 
  if(ESiAx){
    ESReal2RefDataX(0);
    ESPlVdata2PlVX();
    ESGapCoreSepX();
    ESCoreX();
  }
  else{
    ESReal2RefData(0);
    ESReal2RefGeom(0);
    if(FlagPoints){
      ESPlVPoints2PlV(EcRjet,EcZjet,EcNjet);
    }
    else{
      ESPlVdata2PlV();
    }
    ESRef2RealGeom(0);
    if((Fjob&0x02)){
      ECRestoreESC(ESTime);
    }
  }
  ESRealGeom(0);
  ESReInitMetrTnsrCS();
  ECRcs3D2Rcs2D(0);
  ES3DMetrTnsrCS();
  if((Fjob&0x02)){
    Ffail	=0;
    if((ESFail&0x0000000F))
      Ffail	+=1;
    if((ESFail&0x000000F0) || (ESFail&0x0F000000))
      Ffail	+=8;
    if((ESFail&0x00F00000) || (ESFail&0x0000F000))
      Ffail	+=4;
    if((ESFail&0x000F0000) || (ESFail&0x00000F00))
      Ffail	+=2;
    return;
  }
  ECReInitFgY();
  ESInpPlPrCr();
  EC1DEqSolv();
  ESPlProf2Data();
  ECRescue2DGeom(0);
  do{
    if(EC2DEqSolv(0) == 0){
      ECgYcs2dRZcsBound(0);
      ECVirtualPlV2vdrPlV();
      ECgYrc2dRZcs();
      ECvdrz2vrz();
      ECCheckJacob(0);
    }
    else{
      ECRestore2DGeom();
    }
    ECRcs2D2Rcs3D(0);
    if(ESiAx != 0){
      EcMetrTnsrX();
    }
    ESRealGeom(0);
    if(ES3DMetrTnsrCS()){
      ESFail	=0x0FFFFFFF;
      break;
    }
    EC1DEqSolv();
  }while(ECIfEqFound());
  Ffail	=0;
  if((ESFail&0x0000000F))
    Ffail	+=1;
  if((ESFail&0x000000F0) || (ESFail&0x0F000000))
    Ffail	+=8;
  if((ESFail&0x00F00000) || (ESFail&0x0000F000))
    Ffail	+=4;
  if((ESFail&0x000F0000) || (ESFail&0x00000F00))
    Ffail	+=2;
  if(Ffail){
    ECRestore2DGeom();
    ECRcs2D2Rcs3D(0);
    if(ESiAx != 0){
      EcMetrTnsrX();
    }
    ESRealGeom(0);
    ES3DMetrTnsrCS();
    EC1DEqSolv();
    ECEqIterFl	=0;
  }
  else{
    ECEqIterFl	=1;
  }
  ECRefGeom(0);
  if((Fjob&0x01)){
    ECSaveESC();
  }
  if((Fjob&0x08)){
    ESBasicFunctions();
#ifdef H
    ESWriteESI0();
#endif
    ESWriteESI(NULL,ESsr,ESsz,ESaB,ESgH,ESNa,ESNp);
  }
  if((Fjob&4)){
    ESBalStab();
  }
  memcpy(lShm+kI,(void*)&Ffail,kI);
  return;
}
#endif

/**
 * De-initialization. This will deallocate and clean up internal variables.
 * Call this procedure once all desired data have been extracted. You must
 * invoke esc1_ to prevent memory leaks.
 *
 * @author L. E. Zakharov
 * @see esc0_
 */

void esc1_()
{
  ESDeInitEsc();
  ESDeInitEc();
  DeInitSplines();
  DeInitLUdcmp();
  Flesc0=0;
  FlStart=1;
}

void esc0(MachineNm)
     char *MachineNm;
{
  esc0_(MachineNm);
}

void esc0x(MachineNm)
     char *MachineNm;
{
  esc0x_(MachineNm);
}

void ESC0(MachineNm)
     char *MachineNm;
{
  esc0_(MachineNm);
}

void esc(int *Ffail,int *Fjob,int *Mpol,double *sTime)
{
  esc_(Ffail,Fjob,Mpol,sTime);
}

void ESC(int *Ffail,int *Fjob,int *Mpol,double *sTime)
{
  esc_(Ffail,Fjob,Mpol,sTime);
}

void esc1()
{
  esc1_();
}

void ESC1()
{
  esc1_();
}

void ESexit(int i)
{
  esc1_();
  exit(i);
}

#endif
#endif

#ifndef amzl_ESCTR
/* ESC Input interface routines */
#ifndef ablb_PFC
#ifdef PFC

extern double ESaZx;
void tsc2e()
{
  static int Fl=1;

  if(Fl){
    ESEqSolvRadC	=2;
    EcInitPlVacData();
    PFLim2DPlV();
    ESInit2DPlGeom();
    ECRcs3D2Rcs2D(0);
    ECvr2r(0);
    ECvz2z(0);
    Fl	=0;
  }
}
#endif

int JET2escEnd()
{
  EcNjet=0;
  return(0);
}

int JET2esc()
{
  static int Fl=1;
  int k,iT,iB,iL,iR;

  if(Fl){
    fp	=fopen("jet990225.dat","r");
    if(fp == NULL){
      printf("jet990225.dat - file is unavailable%c\n",'\a');
      return(0);
    }
    Fl	=0;
  }
  EcNjet	=0;
  do{
    k	=fscanf(fp,"%lg%lg%lg",EcRjet+EcNjet,EcRjet+EcNjet,EcZjet+EcNjet);
    EcNjet++;
  }while(k == 3);
  fclose(fp);

  iT	=0;
  iB	=0;
  iL	=0;
  iR	=0;
  zdPV[0]	=EcZjet[0];
  rdPV[0]	=EcRjet[0];
  zdPV[1]	=EcZjet[0];
  rdPV[1]	=EcRjet[0];
  zdPV[2]	=EcZjet[0];
  rdPV[2]	=EcRjet[0];
  zdPV[3]	=EcZjet[0];
  rdPV[3]	=EcRjet[0];

  for(k=0; k < EcNjet-1; k++){
    if(zdPV[0] < EcZjet[k]){
      rdPV[0]	=EcRjet[k];
      zdPV[0]	=EcZjet[k];
      iT	=k;
    }
    if(zdPV[1] > EcZjet[k]){
      rdPV[1]	=EcRjet[k];
      zdPV[1]	=EcZjet[k];
      iB	=k;
    }
    if(rdPV[2] > EcRjet[k]){
      rdPV[2]	=EcRjet[k];
      zdPV[2]	=EcZjet[k];
      iL	=k;
    }
    if(rdPV[3] < EcRjet[k]){
      rdPV[3]	=EcRjet[k];
      zdPV[3]	=EcZjet[k];
      iR	=k;
    }
  }
  k	=iT/2;
  zdPV[4]	=EcZjet[k];
  rdPV[4]	=EcRjet[k];
  k	=(iT+iR)/2;
  zdPV[5]	=EcZjet[k];
  rdPV[5]	=EcRjet[k];
  k	=(iR+iB)/2;
  zdPV[6]	=EcZjet[k];
  rdPV[6]	=EcRjet[k];
  k	=(iB+EcNjet-1)/2;
  zdPV[7]	=EcZjet[k];
  rdPV[7]	=EcRjet[k];

  if(ESaR0[0] < EZcr2*(rdPV[3]+rdPV[2])-0.3*EZcr2*(rdPV[3]-rdPV[2]) ||
     ESaR0[0] > EZcr2*(rdPV[3]+rdPV[2])+0.3*EZcr2*(rdPV[3]-rdPV[2])){
    ESaR0[0]	=EZcr2*(rdPV[3]+rdPV[2])+0.1*EZcr2*(rdPV[3]-rdPV[2]);
  }
  if(ESaZ0[0] < EZcr2*(zdPV[0]+zdPV[1])-0.3*EZcr2*(zdPV[0]-zdPV[1]) ||
     ESaZ0[0] > EZcr2*(zdPV[0]+zdPV[1])+0.3*EZcr2*(zdPV[0]-zdPV[1])){
    ESaZ0[0]	=EZcr2*(zdPV[0]+zdPV[1]);
  }

  return(0);
}

/* EFIT-ESC interface. 
   Subject to changes by  F.Paoletti & L.E. Zakharov with no notice. */

static double Ps[257],p[257],q[257];

static int iSplDD,iSplDD1,iSplDDold;
static double crxDD,cHxDD,SplDDx,SplDDxx,SplDDX,SplDDXX;
static double sPs[NDAT],d2sPs[NDAT],sT[NDAT],d2sT[NDAT];


int EFIT2esc()
{
  static int Fl=1;
  int i,I,j,k,K;
  double Rc[32],Zs[32],s,t,dx;
  double *pr,*pz,R0;
  char ch;
  FILE *fp;

  fp=fopen("efit2esc","r");
  if(fp == NULL){
    printf("efit2esc - file is unavailable%c\n",'\a');
    return(-1);
  }
  for(k=0; k < 6; k++){
    while((ch=getc(fp)) != '\n'){
      ;
    }
  }
  while((ch=getc(fp)) != '='){
    putchar(ch);
  }
  fscanf(fp,"%d",&I);
  while((ch=getc(fp)) != '='){
    ;
  }
  fscanf(fp,"%d",&K);
  while((ch=getc(fp)) != '='){
    ;
  }
  while((ch=getc(fp)) != '='){
    ;
  }
  fscanf(fp,"%lg",&R0);
  
  k	=0;
  while(k < K){
    while((ch=getc(fp)) != '='){
      ;
    }
    for(j=0; j < 5 && k < K; j++){
      fscanf(fp,"%lg",Rc+k);
      k++;
    }
  }
  k	=0;
  while(k < K){
    while((ch=getc(fp)) != '='){
      ;
    }
    for(j=0; j < 5 && k < K; j++){
      fscanf(fp,"%lg",Zs+k);
      k++;
    }
  }
  while((ch=getc(fp)) != '='){
    ;
  }
  fscanf(fp,"%lg",&ESBt);
  while((ch=getc(fp)) != '='){
    ;
  }
  fscanf(fp,"%lg",&ESRext);
  ESRBt	=ESBt*ESRext;
  
  while((ch=getc(fp)) != '\n'){
    ;
  }
  while((ch=getc(fp)) != '\n'){
    ;
  }
  fscanf(fp,"%d%lg%lg%lg%lg%lg",&j,&s,q,&s,p,Ps);
  i	=0;
  I	-=2;
  while(i < I){
    fscanf(fp,"%d%lg%lg%lg%lg%lg",&j,&s,q+i,&s,p+i,Ps+i);
    i++;
  }
  fclose(fp);
  if(Fl){
    npv	=12;
    mpv	=npv/2+1;
    ESSetPlVInt(npv,mpv);
  }
#ifdef H
  cssp	=EZcgm0*1e-6;
#endif
  cssp	=-EZcgm0*1e-6*EZc2gp;
  cssj	=1.;
  cssr	=1.;
    
  ESEqSolvFl	=0x31;
  ESEqSolvIt	=20;

#ifdef H
  ESEqSolvRadC	=2;
#endif
#ifdef H
  ESEqSolvInPr	=2;
#endif

  ESEqSolvInPr	=1;
  ESEqSolvInCr	=6;
  np		=NDAT;
  nj		=NDAT;
  ESnDtPr	=np-1;
  ESnDtCr	=nj-1;
  ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
  dx	=1./(np-1);
  for(i=0; i < NDAT; i++){
    ESqgY[i]	=dx*i;
    xdPr[i]	=dx*i;
    xdCr[i]	=dx*i;
    ESxDtPr[i]	=xdPr[i];
    ESxDtCr[i]	=xdCr[i];
  }
  for(i=0; i < ESNa1; i++){
    ESaF[i]	=ESRBt;
    ESaF1a[i]	=0;
    ESaF2a[i]	=0;
  }
  j	=ESNp1*ESNa;
  pr	=EZvr+j;
  pz	=EZvz+j;
  zdPV[0]	=0.;
  zdPV[1]	=0.;
  rdPV[2]	=R0;
  rdPV[3]	=R0;
  for(j=0; j < ESNp1; j++){
    t	=ESgt[j];
    *pr	=R0;
    *pz	=0.;
    for(k=0; k < K; k++){
      s		=(k+1)*t;
      *pr	+=Rc[k]*cos(s);
      *pz	+=Zs[k]*sin(s);
    }
    if(zdPV[0] < *pz){
      rdPV[0]	=*pr;
      zdPV[0]	=*pz;
    }
    if(zdPV[1] > *pz){
      rdPV[1]	=*pr;
      zdPV[1]	=*pz;
    }
    if(rdPV[2] > *pr){
      rdPV[2]	=*pr;
      zdPV[2]	=*pz;
    }
    if(rdPV[3] < *pr){
      rdPV[3]	=*pr;
      zdPV[3]	=*pz;
    }
    pr++;
    pz++;
  }
  j	=ESNp1*ESNa;
  pr	=EZvr+j;
  pz	=EZvz+j;
  ESaR0[0]	=EZcr2*(0.8*rdPV[2]+1.2*rdPV[3]);
  ESaZ0[0]	=0.;
  j	=5;
  rdPV[4]	=pr[j];
  zdPV[4]	=pz[j];
  rdPV[5]	=rdPV[4];
  zdPV[5]	=-zdPV[4];
  j	=ESNp/2-5;
  rdPV[6]	=pr[j];
  zdPV[6]	=pz[j];
  rdPV[7]	=rdPV[6];
  zdPV[7]	=-zdPV[6];
  j	=12;
  rdPV[8]	=pr[j];
  zdPV[8]	=pz[j];
  rdPV[9]	=rdPV[8];
  zdPV[9]	=-zdPV[8];
  j	=ESNp/2-10;
  rdPV[10]	=pr[j];
  zdPV[10]	=pz[j];
  rdPV[11]	=rdPV[10];
  zdPV[11]	=-zdPV[10];

  for(j=0; j < 12; j++){
    printf("??? rdPV[%2d]=%10.3e %10.3e\n",j,rdPV[j],zdPV[j]);
  }

  EFIT2ESPlPrCr();
  if(Fl){
    ESInit2DPlGeom();
    EcInitPlVacData();
    ECRcs3D2Rcs2D(0);
    ECvr2r(0);
    ECvz2z(0);
    ESSwitchInpRadC();
    ESInpPlPrCr();
    ESFirstSetSplDPr();
    ESFirstSetSplDCr();
    Fl	=0;
  }
#ifdef H
#endif
  return(0);
}

int ESFirstSetSplDD()
{
  iSplDD	=0;
  iSplDD1	=1;
  crxDD		=ESxDtPr[1];
  cHxDD		=EZcr6*crxDD*crxDD;
  crxDD		=1./crxDD;
  iSplDDold	=0;
  return(0);
}

int ESSetSplDD(double A)
{
  if(A > ESxDtPr[ESnDtPr])
    A	=ESxDtPr[ESnDtPr];
  if(A < 0.)
    A	=0.;
  while(iSplDD < ESnDtPr && ESxDtPr[iSplDD+1] < A){
    iSplDD++;
  }
  while(iSplDD > 0 && ESxDtPr[iSplDD] > A){
    iSplDD--;
  }
  if(iSplDDold	!= iSplDD){
    iSplDDold	=iSplDD;
    iSplDD1	=iSplDD+1;
    crxDD	=ESxDtPr[iSplDD1]-ESxDtPr[iSplDD];
    cHxDD	=EZcr6*crxDD*crxDD;
    crxDD	=1./crxDD;
  }
  SplDDx	=(A-ESxDtPr[iSplDD])*crxDD;
  SplDDxx	=SplDDx*SplDDx;

  SplDDX	=1.-SplDDx;
  SplDDXX	=SplDDX*SplDDX;
  return(0);
}

int splRDD(double*f, double*df,int k)
{
  double*d0,*d2;
  if(k){
    d0	=sT;
    d2	=d2sT;
  }
  else{
    d0	=sPs;
    d2	=d2sPs;
  }
  *f	=SplDDX*d0[iSplDD]+SplDDx*d0[iSplDD1]
    +(SplDDX*(SplDDXX-1.)*d2[iSplDD]+SplDDx*(SplDDxx-1.)*d2[iSplDD1])*cHxDD;
  if(df != NULL){
    *df	=(d0[iSplDD1]-d0[iSplDD]+
	  ((3.*SplDDxx-1.)*d2[iSplDD1]-(3.*SplDDXX-1.)*d2[iSplDD])*cHxDD)
      *crxDD;
  }
  return(0);
}

int EFIT2ESPlPrCr()
{
  int i;

  ydCr[0]	=sT[0];
  ydPr[0]	=sPs[0];
  for(i=0; i < nj; i++){
    ESSetSplDD(ESqgY[i]);
    splRDD(ydCr+i,NULL,1);
    splRDD(ydPr+i,NULL,0);
  }
  return(0);
}

/* DCON-ESC interface. 
   Subject to changes by L.E. Zakharov with no notice. */

int DCON2ESPlPrCr()
{
  int i,ii,iii,j,k;
  double s,t;

  s	=Ps[0]/(Ps[1]-Ps[0]);
  ydCr[0]	=1./q[0]-(1./q[1]-1./q[0])*s;
  ydPr[0]	=p[0]-(p[1]-p[0])*s;

  Ps[0]	=0.;
  q[0]	=1./ydCr[0];
  p[0]	=ydPr[0];

#ifdef H
  iii	=ESNa1/(nj-1);
  ii	=iii;
  j	=0;
  for(i=1; i < nj; i++){
    s	=ESxDtPr[i];
    ii	+=iii;
    t	=s*s;
    while(Ps[j] < t && j < 128){
      j++;
    }
    k	=j-1;
    s	=(t-Ps[j])/(Ps[k]-Ps[j]);
    ydCr[i]	=1./q[j]+(1./q[k]-1./q[j])*s;
    ydPr[i]	=p[j]+(p[k]-p[j])*s;
  }
#endif
  iii	=ESNa1/(nj-1);
  ii	=iii;
  j	=0;
  for(i=0; i < nj; i++){
    ydCr[i]	=1./q[i];
    ydPr[i]	=p[i];

    xdPr[i]	=sqrt(Ps[i]);
    xdCr[i]	=xdPr[i];
    ESxDtPr[i]	=xdPr[i];
    ESxDtCr[i]	=xdCr[i];
  }
  return(0);
}

int DCON2esc()
{
  static int Fl=1;
  int i,j,na,na1;
  double s;
  double *pr,*pz,a,b,gb;
  char ch,*pc,*str;
  FILE *fp;

  fp=fopen("../dcon_3.10/dcon/binread.EZout","r");
  if(fp == NULL){
    printf("binread.EZout - file is unavailable%c\n",'\a');
    return(-1);
  }
  str	="mpsi";
  pc	=str;
  while(feof(fp) == 0 && *pc != '\0'){
    ch=getc(fp);
    if(ch != *pc){
      pc	=str;
    }
    else{
      pc++;
    }
  }    
  if(*pc != '\0'){
    printf("binread.EZout - cannot find mpsi%c\n",'\a');
    return(-1);
  }
  while(feof(fp) == 0 && ch != '\n'){
    ch	=getc(fp);
  }
  if(ch != '\n' || fscanf(fp,"%d",&na) < 1){
    printf("binread.EZout - cannot get mpsi%c\n",'\a');
    return(-1);
  }
  na1	=na+1;
  printf("na=%d\n",na);

  str	="amean";
  pc	=str;
  while(feof(fp) == 0 && *pc != '\0'){
    ch=getc(fp);
    if(ch != *pc){
      pc	=str;
    }
    else{
      pc++;
    }
  }    
  if(*pc != '\0'){
    printf("binread.EZout - cannot find mpsi%c\n",'\a');
    return(-1);
  }
  while(feof(fp) == 0 && ch != '\n'){
    ch	=getc(fp);
  }
  if(ch != '\n' || fscanf(fp,"%lg%lg%lg%lg",&a,&ESRext,&s,&b) < 4){
    printf("binread.EZout - cannot get a,R0,b%c\n",'\a');
    return(-1);
  }
  b	*=a;
  printf("Rext=%10.3e a=%10.3e b=%10.3e\n",ESRext,a,b);

  str	="betat";
  pc	=str;
  while(feof(fp) == 0 && *pc != '\0'){
    ch=getc(fp);
    if(ch != *pc){
      pc	=str;
    }
    else{
      pc++;
    }
  }    
  if(*pc != '\0'){
    printf("binread.EZout - cannot find mpsi%c\n",'\a');
    return(-1);
  }
  while(feof(fp) == 0 && ch != '\n'){
    ch	=getc(fp);
  }
  if(ch != '\n' || fscanf(fp,"%lg%lg%lg%lg%lg%lg",&s,&s,&s,&gb,&s,&ESBt) < 6){
    printf("binread.EZout - cannot get gb,Bt%c\n",'\a');
    return(-1);
  }
  printf("Bt=%10.3e gb=%10.3e\n",ESBt,gb);
  ESRBt	=ESBt*ESRext;

  str	="ipsi";
  pc	=str;
  while(feof(fp) == 0 && *pc != '\0'){
    ch=getc(fp);
    if(ch != *pc){
      pc	=str;
    }
    else{
      pc++;
    }
  }    
  if(*pc != '\0'){
    printf("binread.EZout - cannot find mpsi%c\n",'\a');
    return(-1);
  }
  while(feof(fp) == 0 && ch != '\n'){
    ch	=getc(fp);
  }
  i	=0;
  while(i < na1){
    if(fscanf(fp,"%d%lg%lg%lg%lg",&j,Ps+i,&s,p+i,q+i) < 5){
      for(j=0; j < i; j++){
	printf("psi[%3d]=%10.3e %12.5e %12.5e %12.5e\n",j,Ps[j],sqrt(Ps[j]),p[j],q[j]);
      }
      printf("binread.EZout - cannot get profile data%c\n",'\a');
      return(-1);
    }
    i++;
  }
  fclose(fp);
#ifdef H
  for(i=0; i < na1; i++){
    printf("psi[%3d]=%10.3e %12.5e %12.5e %12.5e\n",i,Ps[i],sqrt(Ps[i]),p[i],q[i]);
  }
#endif
  ESEqSolvFl	=0x31;
  ESEqSolvIt	=20;
  ESEqSolvRadC	=3;
  ESEqSolvInPr	=2;
  ESEqSolvInCr	=7;
  np		=129;
  nj		=np;
  ESnDtPr	=np-1;
  ESnDtCr	=nj-1;
  ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
  s	=ESsa[ESNa]/ESnDtPr;
#ifdef H
  i	=0;
  xdPr[i]	=s*i;
  xdCr[i]	=s*i;
  ESxDtPr[i]	=xdPr[i];
  ESxDtCr[i]	=xdCr[i];
#endif
  for(i=0; i < np; i++){
    kdPr[i]	=1;
    kdCr[i]	=1;
    xdPr[i]	=s*i;
    xdCr[i]	=s*i;
    ESxDtPr[i]	=xdPr[i];
    ESxDtCr[i]	=xdCr[i];
  }
#ifdef H
  xdPr[i]	=1.;
  xdCr[i]	=1.;
  ESxDtPr[i]	=xdPr[i];
  ESxDtCr[i]	=xdCr[i];
#endif


  for(i=0; i < ESNa1; i++){
    ESaF[i]	=ESRBt;
    ESaF1a[i]	=0;
    ESaF2a[i]	=0;
  }
  j	=ESNp1*ESNa;
  pr	=EZvr+j;
  pz	=EZvz+j;

  rdPV[0]	=ESRext;
  zdPV[0]	=b;
  rdPV[1]	=ESRext;
  zdPV[1]	=-b;
  rdPV[2]	=ESRext-a;
  zdPV[2]	=0.;
  rdPV[3]	=ESRext+a;
  zdPV[3]	=0.;

  for(j=0; j < ESNp1; j++){
    *pr++	=ESRext+a*EScs1[j];
    *pz++	=b*ESsn1[j];
  }
  DCON2ESPlPrCr();
  if(Fl){
    ESaR0[0]	=ESRext;
    ESaZ0[0]	=0.;
    npv	=4;
    mpv	=npv/2+1;
    ESSetPlVInt(npv,mpv);
    ESInit2DPlGeom();
    EcInitPlVacData();
    ECRcs3D2Rcs2D(0);
    ECvr2r(0);
    ECvz2z(0);
    ESSwitchInpRadC();
    ESInpPlPrCr();
    ESFirstSetSplDPr();
    ESFirstSetSplDCr();
    Fl	=0;
  }
  return(0);
}

#endif

#ifndef ablb_ASTRA
static double rc[33],rca[33],rs[33],rsa[33],as[65],as2p[65],as1[65],as12p[65];
static double ch2hp,crHp;
#ifdef ASTRA
/* ASTRA-ESC interface. 
   Subject to changes by  G.V.Pereverzev & L.E. Zakharov with no notice. */

/* ASTRA to ESC input initiation routine */
void astra2e_(double	*Pressure
	      ,double	*Jparallel
	      ,double	*Rpv
	      ,double	*Zpv
	      ,int	*Npv
	      ,double	*RBtor
	      ,double	*Rext
	      ,double	*gFtor
	      )
{
  int i;
 
  cssp	=EZcgm0;
  cssj	=EZcgm0;
  cssr	=1.;
  if(*Npv > 12){
    FlagPoints	=1;
    EcNjet	=*Npv;
    for(i=0; i < EcNjet; i++){
      EcRjet[i]	=Rpv[i];
      EcZjet[i]	=Zpv[i];
    }
  }
  else{
    FlagPoints	=0;
  }
  if(FlStart){
    double dx;

    ch2hp	=EZcr12*ESgt[1]*ESgt[1];
    crHp	=1./ESNp;
    npv		=*Npv < 13 ? *Npv : 12;
    mpv		=FlagPoints ? 8 : npv/2+1;
    ESSetPlVInt(npv,mpv);
    
    ESEqSolvFl	=0x31;
    ESEqSolvIt	=20;
    ESEqSolvRadC	=2;
    ESEqSolvInPr	=2;
    ESEqSolvInCr	=1;
    np		=NDAT;
    nj		=NDAT;
    ESnDtPr	=np-1;
    ESnDtCr	=nj-1;
    ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
    dx	=1./ESnDtPr;
    for(i=0; i < NDAT; i++){
      xdPr[i]	=dx*i;
      xdCr[i]	=dx*i;
      ESxDtPr[i]	=xdPr[i];
      ESxDtCr[i]	=xdCr[i];
    }
    for(i=0; i < ESNa1; i++){
      ESaF[i]	=ESRBt;
      ESaF1a[i]	=0;
      ESaF2a[i]	=0;
    }
  }
  ESRext	=*Rext;
  ESRBt		=*RBtor;
  for(i=0; i < npv; i++){
    rdPV[i]	=Rpv[i]*cssr;
    zdPV[i]	=Zpv[i]*cssr;
  }
  if(FlStart){
    if(FlagPoints){
      double zmx,zmn,rmx,rmn;
      zmx	=EcZjet[0];
      zmn	=zmx;
      rmx	=EcRjet[0];
      rmn	=rmx;
      for(i=1;  i < EcNjet; i++){
	if(zmx < EcZjet[i]){
	  zmx	=EcZjet[i];
	}
	if(zmn > EcZjet[i]){
	  zmn	=EcZjet[i];
	}
	if(rmx < EcRjet[i]){
	  rmx	=EcRjet[i];
	}
	if(rmn > EcRjet[i]){
	  rmn	=EcRjet[i];
	}
      }
      ESaR0[0]	=EZcr2*(0.9*rmn+1.1*rmx);
      ESaZ0[0]	=EZcr2*(zmx+zmn);
      rdPV[0]	=EZcr2*(rmn+rmx);
      zdPV[0]	=zmx;
      rdPV[1]	=EZcr2*(rmn+rmx);
      zdPV[1]	=zmn;
      rdPV[2]	=rmn;
      zdPV[2]	=EZcr2*(zmx+zmn);
      rdPV[3]	=rmx;
      zdPV[3]	=EZcr2*(zmx+zmn);
    }
    else{
      ESaR0[0]	=EZcr2*(0.9*Rpv[2]+1.1*Rpv[3]);
      ESaZ0[0]	=EZcr2*(Zpv[0]+Zpv[1]);
      if(ESiAx != 0){
	ESaZ0[0]	+=0.2*EZcr2*(Zpv[0]-Zpv[1]);
      }
    }
    ESInit2DPlGeom();
    EcInitPlVacData();
    ECRcs3D2Rcs2D(0);
    ECvr2r(0);
    ECvz2z(0);
    ESSwitchInpRadC();
    ESFirstSetSplDPr();
    ESFirstSetSplDCr();
    FlStart	=0;
  }
  for(i=0; i < np; i++){
    ydPr[i]	=Pressure[i]*cssp;
  }
  for(i=0; i < nj; i++){
    ydCr[i]	=Jparallel[i]*cssj;
  }
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
  int i;
  double rR0,R0;
  int j,ji,k,ki,kj;
  double b,ba;
  double ra,rgt,za,zgt;

  R0		=ESaR0[0];
  rR0		=1./R0;
  *Rmaxis	=R0;
  *Zmaxis	=ESaZ0[0];
  *gFtor	=ESgF[ESNa]*EZc2gp;
  for(i=0; i < ESNa1; i++){
    Fpol[i]	=ESaF[i];
    gM[i]	=ESgm[i];
    G22[i]	=ESsa[i]*ESg22c[i];
    G33[i]	=ESsa[i]*ESLc[i]*rR0;
    D2[i]	=ESsa[i]*ESDPc[i];
    D3[i]	=ESsa[i]*R0*(ESLc[i]+ESVc[i]);
  }
  Sside[0]	=0.;
  Sside1[0]	=0.;
  for(i=1; i < ESNa1; i++){
    b		=ESsb[i];
    ba		=ESsb1a[i];
    ki		=i;
    rca[0]	=rcT1a[i];
    rsa[0]	=rsT1a[i];
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rca[k]	=2.*rcT1a[ki];
      rsa[k]	=2.*rsT1a[ki];
      rc[k]	=2.*k*rcT[ki];
      rs[k]	=2.*k*rsT[ki];
    }
    ji	=ESNp1*i;
    for(j=0; j < ESNp; j++){
      ra	=rca[0];
      rgt	=0.;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	ra	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
	rgt	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
      }
      zgt	=b*EScs1[j];
      za	=rsa[0]+ba*ESsn1[j];
      as[j]	=rgt*rgt+zgt*zgt;
      as1[j]	=sqrt(rgt*rgt+zgt*zgt)*ESsr[ji];
      as[j]	*=ESsr[ji]/(ra*zgt-rgt*za);
      ji++;
    }
    as[j]	=as[0];
    as1[j]	=as1[0];
    splP(as,as2p);
    splP(as1,as12p);
    b	=0.;
    ba	=0.;
    for(j=0; j < ESNp; j++){
      b		+=as[j]-ch2hp*as2p[j];
      ba	+=as1[j]-ch2hp*as12p[j];
    }
    Sside[i]	=b*crHp;
    Sside1[i]	=ba*crHp;
  }
}

void EZastra2e(double	*Pressure
	      ,double	*Jparallel
	      ,double	*Rpv
	      ,double	*Zpv
	      ,int	*Npv
	      ,double	*RBtor
	      ,double	*Rext
	      ,double	*gFtor
	      )
{
  astra2e_(Pressure,Jparallel,Rpv,Zpv,Npv,RBtor,Rext,gFtor);
}

void ASTRA2E(double	*Pressure
	      ,double	*Jparallel
	      ,double	*Rpv
	      ,double	*Zpv
	      ,int	*Npv
	      ,double	*RBtor
	      ,double	*Rext
	      ,double	*gFtor
	      )
{
  astra2e_(Pressure,Jparallel,Rpv,Zpv,Npv,RBtor,Rext,gFtor);
}

void e2astra(double	*Rmaxis
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
  e2astra_(Rmaxis,Zmaxis,gFtor,Fpol,gM,G22,G33,D2,D3,Sside,Sside1);
}

void E2ASTRA(double	*Rmaxis
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
  e2astra_(Rmaxis,Zmaxis,gFtor,Fpol,gM,G22,G33,D2,D3,Sside,Sside1);
}

#ifdef SHMEM
void E2Astra()
{
  int i,kA;
  double rR0,R0;
  int j,ji,k,ki,kj;
  double b,ba;
  double ra,rgt,za,zgt;
  double *ld;
  double Sside[ESNa1],Sside1[ESNa1];

  kA	=kD*ESNa1;
  ld	=(double*)(lShm+kI);
  *ld++	=ESaR0[0];
  *ld++	=ESaZ0[0];
  *ld++	=ESgF[ESNa]*EZc2gp;
  memcpy((void*)ld,(void*)ESaF,kA);
  ld	+=ESNa1;
  memcpy((void*)ld,(void*)ESgm,kA);
  ld	+=ESNa1;

  R0		=ESaR0[0];
  rR0		=1./R0;
  /* G22 */
  for(i=0; i < ESNa1; i++){
    *ld++	=ESsa[i]*ESg22c[i];
  }
  /* G33 */
  for(i=0; i < ESNa1; i++){
    *ld++	=ESsa[i]*ESLc[i]*rR0;
  }
  /* D2 */
  for(i=0; i < ESNa1; i++){
    *ld++	=ESsa[i]*ESDPc[i];
  }
  /* D3 */
  for(i=0; i < ESNa1; i++){
    *ld++	=ESsa[i]*R0*(ESLc[i]+ESVc[i]);
  }

  Sside[0]	=0.;
  Sside1[0]	=0.;
  for(i=1; i < ESNa1; i++){
    b		=ESsb[i];
    ba		=ESsb1a[i];
    ki		=i;
    rca[0]	=rcT1a[i];
    rsa[0]	=rsT1a[i];
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rca[k]	=2.*rcT1a[ki];
      rsa[k]	=2.*rsT1a[ki];
      rc[k]	=2.*k*rcT[ki];
      rs[k]	=2.*k*rsT[ki];
    }
    ji	=ESNp1*i;
    for(j=0; j < ESNp; j++){
      ra	=rca[0];
      rgt	=0.;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	ra	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
	rgt	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
      }
      zgt	=b*EScs1[j];
      za	=rsa[0]+ba*ESsn1[j];
      as[j]	=rgt*rgt+zgt*zgt;
      as1[j]	=sqrt(rgt*rgt+zgt*zgt)*ESsr[ji];
      as[j]	*=ESsr[ji]/(ra*zgt-rgt*za);
      ji++;
    }
    as[j]	=as[0];
    as1[j]	=as1[0];
    splP(as,as2p);
    splP(as1,as12p);
    b	=0.;
    ba	=0.;
    for(j=0; j < ESNp; j++){
      b		+=as[j]-ch2hp*as2p[j];
      ba	+=as1[j]-ch2hp*as12p[j];
    }
    Sside[i]	=b*crHp;
    Sside1[i]	=ba*crHp;
  }
  memcpy((void*)ld,(void*)Sside,kA);
  ld	+=ESNa1;
  memcpy((void*)ld,(void*)Sside1,kA);
  ld	+=ESNa1;
  return;
}
#endif

#endif

#ifdef H
/* Our old interface */
sesc_(double*agrid,double*Pressure,double*Jparallel,int*Ngrids,
      double*Rpv,double*Zpv,int*Npv,double*Rmaxis,double*gFtor,
      double*Fpol,double*gM,double*G22,double*G33,double*D2,
      double*D3,double*Sside,double*Sside1,double*TTime)
#endif

#endif



#ifdef ablb_TRANSP

static float *Rc,*Rs,*Zc,*Zs;
static float *pTRANSP,*gmTRANSP;
static double *ydPr,*ydCr,*xdCr,cxx1,cxx2;

void esc2tr0_(Mp,p,gi,rc,rs,zc,zs)	/* ESC initiation routine in TRANSP */
     int *Mp; 	/* highest # of m in poloidal Fourier harmonic in $\Psi$ */
     float *p;	/* address of array with data on pressure profile
		   from TRANSP in Pa */
     float *gi;	/* address of array with data on iota profile
		   from TRANSP */
     float *rc;	/* address of array with data on Rc Fourier coefficients */
     float *rs;	/* address of array with data on Rs Fourier coefficients */
     float *zc;	/* address of array with data on Zc Fourier coefficients */
     float *zs;	/* address of array with data on Zs Fourier coefficients */
{
  /* Customize the following two numbers */
  ESNa	=20;
  /* number of radial intervals in the plasma core */
  ESNp	=64;
  /* number of poloidal intervals in the plasma core
     (for calculations of Fourier coefficients) */
  
  /* ======= Do not touch following (Tue Oct 27 09:53:28 EST 1998) ======= */

  pTRANSP	=p+1;
  gmTRANSP	=gi+1;
  Rc		=rc;
  Rs		=rs;
  Zc		=zc;
  Zs		=zs;
  
  ESMp	=*Mp; /* highest # of Fourier harmonic in $\Psi$; */
  ESFp	=2*ESMp; /* maximum # of poloidal harmonic in representation of r */
  TFCNt	=0;  /* Number of toroidal coils in one EZmain period of the machine */
  TFCNp	=1;
  /* number poloidal intervals in straight piece-wise 
     representation of TF Coils */
  ESLt	=0;
  /* L >= 0;
     Basic toroidal wave number of the system:;
     0 - axisymmetric configurations (tokamaks, RFP,RFT etc.);;
     1 - 1-period stellarators, $x,y \propto \cos(\phi)$;;
     2 - 2-period stellarators, $x,y \propto \cos(2\phi)$;;
     ... etc. */
  
  ESnLm	=0;
  ESNt	=0;  /* ESNt >= 0, Number of toroidal intervals */

  InitConstants();
  ESInitPlGeomInt();
  ESReInit1DPlV();
  ESReInit2DTPlV();

  EcInitArrays();
  ESReInitPlProf();
  ESInit1DPlV();
  ESInit2DTPlV();

  ESInit2DPlGeom();
  EcInitPlVacData();
  ESReInitMetrTnsrCS();
  ECReInitFgY();

  ESEqSolvFl	=0x31;
  ESEqSolvIt	=20;

  ESPlPrAddr2TRANSP(&ydPr,&ydCr,&xdCr,&kdPr,&kdCr);

  {
    double s;
    cxx1	=0.025*0.025;
    cxx2	=0.075*0.075;
    s		=1./(cxx2-cxx1);
    cxx2	*=s;
    cxx1	*=s;
  }
}

void esc2tr1_(rPlV,zPlV,nPlV,Rext,Btext)/* ESC call routine in TRANSP */
     double *rPlV;	/* r[cm] of Plasma-Vacuum in the following order:
			   top, bottom left, right, ... points */
     double *zPlV;	/* z[cm] of Plasma-Vacuum in the following order:
			   top, bottom left, right,  ... points */
     int *nPlV; 	/* number of Plasma-Vacuum points, 3 < nPlV < 13 */
     double *Rext;	/* Rext[cm] - external reference major radius */
     double *Btext;	/* Bext[T] - external field at Rext */
{
  static Fl=0;
  int i;
  double dx;

  dx	=1./(NDAT-1);
  ECGetPlVData(rPlV,zPlV,*nPlV,*Btext);
  ESRext	=*Rext*0.01;
  ESRBt		=ESRext*(*Btext);

 
  do{
    if(Fl == 0){
      R0T[0]	=(0.5*(rPlV[2]+rPlV[3])+0.1*(rPlV[3]-rPlV[2]))*0.01;
      Z0T[0]	=0.5*(zPlV[0]+zPlV[1])*0.01;
      for(i=0; i < NDAT; i++){
	ESxDtPr[i]	=0.05*i;
	ESxDtCr[i]	=0.05*i;
      }
    }
    if(Fl < 4){
      Fl++;
    }
    if(Fl == 3){
      int nz;
      nz	=NDAT-1;
      for(i=0; i < nz; i++){
	xdCr[i]	=dx*(i+0.5);
      }
      xdCr[0]	=0.;
      xdCr[nz]	=1.;
      for(i=0; i < NDAT; i++){
	ESxDtCr[i]	=xdCr[i];
      }
      for(i=0; i < nz; i++){
	ydCr[i]	=gmTRANSP[i];
      }
      ydCr[nz]	=ydCr[nz-1]+0.5*(ydCr[nz-1]-ydCr[nz-2]);
      ydCr[0]	=ydCr[0]*cxx2-ydCr[1]*cxx1;
      ESEqSolvInCr=7;
    }
    if(Fl == 4){
      ESEqSolvInPr=2;
      for(i=0; i < NDAT; i++){
	ydPr[i]	=pTRANSP[i]*1e-6*EZcgm0;
      }
      for(i=0; i < nz; i++){
	ydCr[i]	=gmTRANSP[i];
      }
      ydCr[nz]	=ydCr[nz-1]+0.5*(ydCr[nz-1]-ydCr[nz-2]);
      ydCr[0]	=ydCr[0]*cxx2-ydCr[1]*cxx1;
    }

    ESInpRadC();
    ESInputPlVacTr();
    ECRcs3D2Rcs2D(0);
    ES3DMetrTnsrCS();
    ESInpPlPrTr();
    ESInpPlCrTr();
    do{
      EC2DEqSolv(0);
      ECgYcs2dRZcsBound(0);
      ECVirtualPlV2vdrPlV();
      ECgYrc2dRZcs();
      ECvdrz2vrz();
      ECCheckJacob(0);
      ECRcs2D2Rcs3D(0);
      ESRealGeom(0);
      ES3DMetrTnsrCS();
      switch(ESEqSolvInCr){
      case 0:
	EC1DEqSolv(0);
	break;
      case 6:
      case 7:
	switch(ESEqSolvInPr){
	case 2:
	  EC1DEqpSolv(0);
	  break;
	default :
	  EC1DEqqSolv(0);
	  break;
	}
	break;
      default:
	EC1DEqSolv(0);
	break;
      }
    }while(ECIfEqFound());
  }while(Fl < 4);

  {
    int k,ki;
    float *pRc,*pRs,*pZc,*pZs;

    pRc	=Rc;
    pRs	=Rs;
    pZc	=Zc;
    pZs	=Zs;

    k	=0;
    ki	=0;
    pRc++;
    pRs++;
    pZc++;
    pZs++;
    for(i=0; i < ESNa1; i++){
      *pRc++	=(R0T[0]+rcT[ki])*100.;
      *pRs++	=0.;
      *pZc++	=(Z0T[0]+rsT[ki])*100.;
      *pZs++	=0.;
      ki++;
    }
    k++;
    pRc++;
    pRs++;
    pZc++;
    pZs++;
    for(i=0; i < ESNa1; i++){
      *pRc++	=200.*rcT[ki];
      *pRs++	=200.*rsT[ki];
      *pZs++	=bT[i]*100.;
      *pZc++	=0.;
      ki++;
    }
    while(k < ESMp){
      k++;
      pRc++;
      pRs++;
      pZc++;
      pZs++;
      for(i=0; i < ESNa1; i++){
	*pRc++	=200.*rcT[ki];
	*pRs++	=200.*rsT[ki];
	*pZs++	=0.;
	*pZc++	=0.;
	ki++;
      }
    }
  }
}

void esc2tr2_()	/* ESC deinitiation routine in TRANSP */
{
  ESDeInitEsc();
  ECDeInitFgY();
  EcDeInitArrays();
  DeInitSplines();
}

#endif
#endif

/* TRANSP to ESC intput routine */
void tr2esc_(double	*Pressure	/* p[Pa] */
	      ,double	*Jparallel	/* {<B\cdot j>\over<B\cdot\na\phi>}
					   [cm A/cm^2].
					   Right hand side of the magnetic 
					   flux evolution equation	
					   */
	      ,double	*Rpv		/* R [cm] of *Npv plasma-vacuum 
					   points. 
					   If *Npv < 13, first 4 points 
					   should be taken in the order:
					   top, bottom, left, right
					   */
	      ,double	*Zpv		/* Z [cm] of *Npv plasma-vacuum 
					   points.
					   If *Npv < 13, first 4 points
					   should be taken in the order:
					   top, bottom, left, right
					   */
	      ,int	*Npv		/* number of the plasma-vacuum 
					   points <= 257
					   If *Npv >12, first and last points
					   coincide */
	      ,double	*RBtor		/* RBtor [cm Tesla] outside the 
					   plasma */
	      ,double	*Rext		/* Reference major radius [cm] */
	      ,double	*gFtor		/* Total toroidal flux [Vsec] 
					   through the plasma cross-section. 
					   
					   !!! Is not used in the moment !!!
					   */
	      )
{
  int i;


  cssp  =1.; /* EZcgm0; */
  cssj  =1.; /* EZcgm0; */
  cssr  =1.;
  if(*Npv > 12){
    FlagPoints  =1;
    EcNjet      =*Npv;
    for(i=0; i < EcNjet; i++){
      EcRjet[i] =Rpv[i];
      EcZjet[i] =Zpv[i];
    }
  }
  else{
    FlagPoints  =0;
  }
  if(FlStart){
    double dx;
    ch2hp       =EZcr12*ESgt[1]*ESgt[1];
    crHp        =1./ESNp;
    npv         =*Npv < 13 ? *Npv : 12;
    mpv         =FlagPoints ? 8 : npv/2+1;
    ESSetPlVInt(npv,mpv);

    ESEqSolvFl  =0x31;
    ESEqSolvIt  =20;
    ESEqSolvRadC        =2;
    ESEqSolvInPr        =2;
    ESEqSolvInCr        =1;
    np          =NDAT;
    nj          =NDAT;
    ESnDtPr     =np-1;
    ESnDtCr     =nj-1;
    ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
    dx	=1./ESnDtPr;
    for(i=0; i < NDAT; i++){
      xdPr[i]   =dx*i;
      xdCr[i]   =dx*i;
      ESxDtPr[i]        =xdPr[i];
      ESxDtCr[i]        =xdCr[i];
    }
    for(i=0; i < ESNa1; i++){
      ESaF[i]   =ESRBt;
      ESaF1a[i] =0;
      ESaF2a[i] =0;
    }
  }
  ESRext        =*Rext;
  ESRBt         =*RBtor;
  for(i=0; i < npv; i++){
    rdPV[i]     =Rpv[i]*cssr;
    zdPV[i]     =Zpv[i]*cssr;
  }
  if(FlStart){
    if(FlagPoints){
      double zmx,zmn,rmx,rmn;
      zmx       =EcZjet[0];
      zmn       =zmx;
      rmx       =EcRjet[0];
      rmn       =rmx;
      for(i=1;  i < EcNjet; i++){
        if(zmx < EcZjet[i]){
          zmx   =EcZjet[i];
        }
        if(zmn > EcZjet[i]){
          zmn   =EcZjet[i];
        }
        if(rmx < EcRjet[i]){
          rmx   =EcRjet[i];
        }
        if(rmn > EcRjet[i]){
          rmn   =EcRjet[i];
        }
      }
      ESaR0[0]  =EZcr2*(0.9*rmn+1.1*rmx);
      ESaZ0[0]  =EZcr2*(zmx+zmn);
      rdPV[0]   =EZcr2*(rmn+rmx);
      zdPV[0]   =zmx;
      rdPV[1]   =EZcr2*(rmn+rmx);
      zdPV[1]   =zmn;
      rdPV[2]   =rmn;
      zdPV[2]   =EZcr2*(zmx+zmn);
      rdPV[3]   =rmx;
      zdPV[3]   =EZcr2*(zmx+zmn);
    }
    else{
      ESaR0[0]  =EZcr2*(0.9*Rpv[2]+1.1*Rpv[3]);
      ESaZ0[0]  =EZcr2*(Zpv[0]+Zpv[1]);
      if(ESiAx != 0){
        ESaZ0[0]        +=0.2*EZcr2*(Zpv[0]-Zpv[1]);
      }
    }
    ESInit2DPlGeom();
    EcInitPlVacData();
    ECRcs3D2Rcs2D(0);
    ECvr2r(0);
    ECvz2z(0);
    ESSwitchInpRadC();
    ESFirstSetSplDPr();
    ESFirstSetSplDCr();
    FlStart     =0;
  }
  for(i=0; i < np; i++){
    ydPr[i]     =Pressure[i]*cssp;
  }
  for(i=0; i < nj; i++){
    ydCr[i]     =Jparallel[i]*cssj;
  }
}


/* ESC to TRANSP output routine */
void esc2tr_(double *Rmaxis	/* R [cm] of the magnetic axis */
	     ,double *Zmaxis	/* Z [cm] of the magnetic axis */
	     ,double *gFtor	/* total toroidal flux [Vsec]
				   inside the plasma */
	     ,double *gYpol	/* array of poloidal flux profile [Vsec] */
	     ,double *RBtor	/* array of RBtor profile [cm Tesla] */
	     ,double *gM	/* array of 1/q profile */
	     ,int *nRadGrid	/* number of grid points for plasma profiles,
				   starting from the magnetic axis
				   up to the plasma boundary */
	     ,int *MStorage	/* Second dimensions of the Rcs,Zcs arrays */
	     ,int *NaStorage	/* Third dimensions of the Rcs,Zcs arrays */
	     ,double *Rcs	/* Fourier coefficients
				   for major radius R:
				   R=\sum_{m=0}^{*nMoments}\left(
				   r^c_m\cos m\theta
				   +r^s_m\cos m\theta
				   \right)
				   FORTRAN 3D array for Fourier coefficients:
				   Rcs(i,m,1)=r^c_m(a_i)
				   Rcs(i,m,2)=r^s_m(a_i)
				   */

	     ,double *Zcs	/* Fourier coefficients
				   for Z:
				   Z=\sum_{m=0}^{*nMoments}\left(
				   z^c_m\cos m\theta
				   +z^s_m\cos m\theta
				   \right)
				   FORTRAN 3D array for Fourier coefficients:
				   Zcs(i,m,1)=r^c_m(a_i)
				   Zcs(i,m,2)=r^s_m(a_i)
				   */
	     ,int *nMoments	/* number of moments in Fourier
				   representation of R and Z:
				   R=\sum_{m=0}^{*nMoments}\left(
				   r^c_m\cos m\theta
				   +r^s_m\cos m\theta
				   \right)
				   Z=\sum_{m=0}^{*nMoments}\left(
				   z^c_m(i)\cos m\theta
				   +z^s_m(i)\cos m\theta
				   \right)
				   */
	     )
{
  int i,na1,m,mi,M,N,ni;
  double R0,x,dx,s;
  double *rc,*rs,*zc,*zs;

  M=*nMoments;
  N =*NaStorage;
  rc =Rcs;
  rs =Rcs+*MStorage*N;
  zc =Zcs;
  zs =Zcs+*MStorage*N;

  R0  =ESaR0[0];
  *Rmaxis =R0*100.;
  *Zmaxis =ESaZ0[0]*100.;
  *gFtor =ESgF[ESNa]*EZc2gp;

  na1 =*nRadGrid;
  dx =ESsa[ESNa]/(na1-1);
  for(i=0; i < na1; i++){
    x  =dx*i;
    ESSetSplA(x);
    splRA(RBtor+i,NULL,ESaF,ESaF2a);
    RBtor[i] *=100.;
    splRA(gYpol+i,NULL,ESgY,ESgY2a);
    gYpol[i] *=EZc2gp;
    splRA(gM+i,NULL,ESgm,ESgm2a);
    mi =0;
    ni =i;
    splRA(&s,NULL,rcT,rcT2a);
    rc[i] =(ESaR0[0]+s)*100.;
    rs[i] =0.;
    splRA(&s,NULL,rsT,rsT2a);
    zc[i] =(ESaZ0[0]+s)*100.;
    zs[i] =0.;
    for(m=1; m <= M; m++){
      mi +=ESNa1;
      ni +=N;
      splRA(&s,NULL,rcT+mi,rcT2a+mi);
      rc[ni] =s*200.;
      splRA(&s,NULL,rsT+mi,rsT2a+mi);
      rs[ni] =s*200.;
    }
    splRA(&s,NULL,ESsb,ESsb2a);
    ni =i+N;
    zc[ni] =0.;
    zs[ni] =s*100.;
    for(m=2; m <= M; m++){
      ni +=N;
      zc[ni] =0.;
      zs[ni] =0.;
    }
  }
}

/**
 * General input profile selection routine. 
 * Various combinations of pressure/current profiles can be chosen
 * as equilibrium input. *Note* that some profiles are in rationalized units while
 * others are in MKS units.
 *
 * @param pProf pressure-like profile, see @a hH switch, an array of size @a nProf (input). 
 * @param jProf current-like profile, see @a hH switch, an array of size @a nProf (input).
 * @param aProf normalized sqrt(Phi) square root of toroidal flux radial coordinate (input).
 * @param hH two-digit pressure and current profile specification.  
 *           First digit (i.e., hH/10) can be 0 (j_p = R0 dp/dpsi), 1 (P=dp/dpsi) or 2 (p in MPa).
 *           Second digit (i.e., hH%10) can be 0 (js=F/R0 + R0 dp/dpsi in MA/m^2), 1 (j//=<J.B>/(R0 B.gradPhi> in MA/m^2), 3 (FF'), 6 (q), 
 *           7 (1/q). 
 *           For example set hH=26 to supply  p and q, hH=21 to supply p and j//, hH=0 to supply jp and js. 
 *           Possible choices of hH are 0,1,6,7, 10,11,16,17, 21,26,27.
 *
 * Note: 
 *
 * 1) Unless otherwise stated pressure is in rationalized units (pressure in mu0 Pa)
 *
 * 2) Currents are in MA: 1 MA = 1 [T/m] / (0.4*pi)
 *
 * 3) R0 is the magnetic axis.
 *
 * @param Rpv array[Npv] for the plasma-vacuum boundary radius in [m] (input).
 * @param Zpv array[Npv] for the plasma-vacuum boundary elevation in [m] (input).
 * @param Npv number of the plasma-vacuum input boundary points.
 * @param RBtor covariant toroidal magnetic field in vacuum [T m] (input).
 * @param Rext a reference major radius in [m] (input). Specifies the 
 * reference vacuum magnetic field as Bext= RBext/Rext, which may enter in 
 * definitions of non-invariant variables, like beta. 
 *
 * 
 * @author L. E. Zakharov
 * @see Esc_p_q, Esc_p_jparallel
 */


void escin_(double	*pProf	/* p[] - pressure related profile */
	    ,double	*jProf	/* j[] - current related profile */
	    ,double	*aProf	/* a[] - normalized sqrt(Phi) square root of toroidal flux.
				   If
				   a=NULL (dropped in FORTRAN) - uniform grid from 0 till 1,
				   otherwise it is considered as an array of grid points */
	    ,int	*nProf	/* number <= 257 of profile points including \
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
  int i,k;
  static double *A=NULL;
 
  cssr	=1.;
  if(*Npv > 12){
    FlagPoints	=1;
    EcNjet	=*Npv;
    for(i=0; i < EcNjet; i++){
      EcRjet[i]	=Rpv[i];
      EcZjet[i]	=Zpv[i];
    }
  }
  else{
    FlagPoints	=0;
  }

  i	=(*hH)/10;
  if(ESEqSolvInPr != i){
    np	=0;
    ESEqSolvInPr	=i;
    switch(i){
    case 0:	/* convert jp [MA/m^2] to gm_0 jp */
    case 2:	/* convert p [MPa] to gm_0 p */
      cssp	=EZcgm0;
      break;
    default:
      cssp	=1.;
      break;
    }
  }
  k	=(*hH)%10;
  i	=k != 6 ? k : 7;
  if(ESEqSolvInCr != i){
    np	=0;
    ESEqSolvInCr	=i;
    switch(i){
    case 0:	/* convert js [MA/m^2] to gm_0 js */
    case 1:	/* convert j|| [MA/m^2] to gm_0 j|| */
      cssj	=EZcgm0;
      break;
    default:
      cssj	=1.;
      break;
    }
  }
  if(FlStart || np == 0){
    ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
  }

  ESRext	=*Rext;
  ESRBt		=*RBtor;
  if(FlStart){
    ch2hp	=EZcr12*ESgt[1]*ESgt[1];
    crHp	=1./ESNp;
    npv		=*Npv < 13 ? *Npv : 12;
    mpv		=FlagPoints ? 8 : npv/2+1;
    ESSetPlVInt(npv,mpv);
    ESEqSolvFl	=0x31;
    ESEqSolvIt	=20;
    ESEqSolvRadC=2;
    for(i=0; i < ESNa1; i++){
      ESaF[i]	=ESRBt;
      ESaF1a[i]	=0;
      ESaF2a[i]	=0;
    }
  }
  for(i=0; i < npv; i++){
    rdPV[i]	=Rpv[i]*cssr;
    zdPV[i]	=Zpv[i]*cssr;
  }
  if(FlStart){
    if(FlagPoints){
      double zmx,zmn,rmx,rmn;
      zmx	=EcZjet[0];
      zmn	=zmx;
      rmx	=EcRjet[0];
      rmn	=rmx;
      for(i=1;  i < EcNjet; i++){
	if(zmx < EcZjet[i]){
	  zmx	=EcZjet[i];
	}
	if(zmn > EcZjet[i]){
	  zmn	=EcZjet[i];
	}
	if(rmx < EcRjet[i]){
	  rmx	=EcRjet[i];
	}
	if(rmn > EcRjet[i]){
	  rmn	=EcRjet[i];
	}
      }
      ESaR0[0]	=EZcr2*(0.9*rmn+1.1*rmx);
      ESaZ0[0]	=EZcr2*(zmx+zmn);
      rdPV[0]	=EZcr2*(rmn+rmx);
      zdPV[0]	=zmx;
      rdPV[1]	=EZcr2*(rmn+rmx);
      zdPV[1]	=zmn;
      rdPV[2]	=rmn;
      zdPV[2]	=EZcr2*(zmx+zmn);
      rdPV[3]	=rmx;
      zdPV[3]	=EZcr2*(zmx+zmn);
    }
    else{
      ESaR0[0]	=EZcr2*(0.9*Rpv[2]+1.1*Rpv[3]);
      ESaZ0[0]	=EZcr2*(Zpv[0]+Zpv[1]);
      if(ESiAx != 0){
	ESaZ0[0]	+=0.2*EZcr2*(Zpv[0]-Zpv[1]);
      }
    }
    ESInit2DPlGeom();
    EcInitPlVacData();
    ECRcs3D2Rcs2D(0);
    ECvr2r(0);
    ECvz2z(0);
    ESSwitchInpRadC();
    FlStart	=0;
  }

  if(np != *nProf){
    double dx;
    np	=*nProf;
    nj	=np;
    ESnDtPr	=np-1;
    ESnDtCr	=nj-1;
    dx		=1./ESnDtPr;
    for(i=0; i < np; i++){
      xdPr[i]	=dx*i;
      xdCr[i]	=dx*i;
      ESxDtPr[i]=xdPr[i];
      ESxDtCr[i]=xdCr[i];
      kdPr[i]	=1;
      kdCr[i]	=1;
    }
  }
  if(aProf != NULL){
    for(i=0; i < np; i++){
      xdPr[i]	=aProf[i];
      xdCr[i]	=aProf[i];
      ESxDtPr[i]=aProf[i];
      ESxDtCr[i]=aProf[i];
    }
    ESFirstSetSplDPr();
    ESFirstSetSplDCr();
  }
  else{
    if(A != NULL){
      double dx;
      dx	=1./ESnDtPr;
      for(i=0; i < np; i++){
	xdPr[i]	=dx*i;
	xdCr[i]	=dx*i;
	ESxDtPr[i]=xdPr[i];
	ESxDtCr[i]=xdCr[i];
      }
      ESFirstSetSplDPr();
      ESFirstSetSplDCr();
    }
  }
  A	=aProf;
  for(i=0; i < np; i++){
    ydPr[i]	=pProf[i]*cssp;
    ydCr[i]	=jProf[i]*cssj;
  }
  if(k == 6){
    for(i=0; i < np; i++){
      ydCr[i]	=1./jProf[i];
    }
  }
  return;
}
void escin(double	*pProf 
	    ,double	*jProf
	    ,double	*aProf 
	    ,int	*nProf 
	    ,int	*hH	
	    ,double	*Rpv
	    ,double	*Zpv
	    ,int	*Npv
	    ,double	*RBtor
	    ,double	*Rext
	    )
{escin_(pProf,jProf,aProf,nProf, hH,Rpv,Zpv,Npv,RBtor,Rext);}

#ifdef SHMEM
void EscIn()
{
  int i,k;
  static double *A=NULL;

  double *pProf,*jProf,*aProf;
  int nProf,hH,Npv,kP,kB;
  double *Rpv,*Zpv;
  double RBtor,Rext;
  char ch;
  void *lv;

  lv	=lShm+kI;
  ch	=*(char*)lv;
  lv++;
  nProf	=*(int*)lv;

  kP	=(nProf)*kD;
  lv	+=kI;
  if(ch == 'y'){
    aProf	=(double*)lv;
    lv	+=kP;
  }
  pProf	=(double*)lv;
  lv	+=kP;
  jProf	=(double*)lv;
  lv	+=kP;
  hH	=*(int*)lv;
  lv	+=kI;
  Npv	=*(int*)lv;
  lv	+=kI;
  kB	=(Npv)*kD;
  Rpv	=(double*)lv;
  lv	+=kB;
  Zpv	=(double*)lv;
  lv	+=kB;
  RBtor	=*(double*)lv;
  lv	+=kD;
  Rext	=*(double*)lv;
  cssr	=1.;
  if(Npv > 12){
    FlagPoints	=1;
    EcNjet	=Npv;
    for(i=0; i < EcNjet; i++){
      EcRjet[i]	=Rpv[i];
      EcZjet[i]	=Zpv[i];
    }
  }
  else{
    FlagPoints	=0;
  }
  
  i	=(hH)/10;
  if(ESEqSolvInPr != i){
    np	=0;
    ESEqSolvInPr	=i;
    switch(i){
    case 0:	/* convert jp [MA/m^2] to gm_0 jp */
    case 2:	/* convert p [MPa] to gm_0 p */
      cssp	=EZcgm0;
      break;
    default:
      cssp	=1.;
      break;
    }
  }
  k	=(hH)%10;
  i	=k != 6 ? k : 7;
  if(ESEqSolvInCr != i){
    np	=0;
    ESEqSolvInCr	=i;
    switch(i){
    case 0:	/* convert js [MA/m^2] to gm_0 js */
    case 1:	/* convert j|| [MA/m^2] to gm_0 j|| */
      cssj	=EZcgm0;
      break;
    default:
      cssj	=1.;
      break;
    }
  }
  if(FlStart || np == 0){
    ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
  }

  ESRext	=Rext;
  ESRBt		=RBtor;
  if(FlStart){
    ch2hp	=EZcr12*ESgt[1]*ESgt[1];
    crHp	=1./ESNp;
    npv		=Npv < 13 ? Npv : 12;
    mpv		=FlagPoints ? 8 : npv/2+1;
    ESSetPlVInt(npv,mpv);
    ESEqSolvFl	=0x31;
    ESEqSolvIt	=20;
    ESEqSolvRadC=2;
    for(i=0; i < ESNa1; i++){
      ESaF[i]	=ESRBt;
      ESaF1a[i]	=0;
      ESaF2a[i]	=0;
    }
  }
  for(i=0; i < npv; i++){
    rdPV[i]	=Rpv[i]*cssr;
    zdPV[i]	=Zpv[i]*cssr;
  }
  if(FlStart){
    if(FlagPoints){
      double zmx,zmn,rmx,rmn;
      zmx	=EcZjet[0];
      zmn	=zmx;
      rmx	=EcRjet[0];
      rmn	=rmx;
      for(i=1;  i < EcNjet; i++){
	if(zmx < EcZjet[i]){
	  zmx	=EcZjet[i];
	}
	if(zmn > EcZjet[i]){
	  zmn	=EcZjet[i];
	}
	if(rmx < EcRjet[i]){
	  rmx	=EcRjet[i];
	}
	if(rmn > EcRjet[i]){
	  rmn	=EcRjet[i];
	}
      }
      ESaR0[0]	=EZcr2*(0.9*rmn+1.1*rmx);
      ESaZ0[0]	=EZcr2*(zmx+zmn);
      rdPV[0]	=EZcr2*(rmn+rmx);
      zdPV[0]	=zmx;
      rdPV[1]	=EZcr2*(rmn+rmx);
      zdPV[1]	=zmn;
      rdPV[2]	=rmn;
      zdPV[2]	=EZcr2*(zmx+zmn);
      rdPV[3]	=rmx;
      zdPV[3]	=EZcr2*(zmx+zmn);
    }
    else{
      ESaR0[0]	=EZcr2*(0.9*Rpv[2]+1.1*Rpv[3]);
      ESaZ0[0]	=EZcr2*(Zpv[0]+Zpv[1]);
      if(ESiAx != 0){
	ESaZ0[0]	+=0.2*EZcr2*(Zpv[0]-Zpv[1]);
      }
    }
    ESInit2DPlGeom();
    EcInitPlVacData();
    ECRcs3D2Rcs2D(0);
    ECvr2r(0);
    ECvz2z(0);
    ESSwitchInpRadC();
    FlStart	=0;
  }

  if(np != nProf){
    double dx;
    np	=nProf;
    nj	=np;
    ESnDtPr	=np-1;
    ESnDtCr	=nj-1;
    dx		=1./ESnDtPr;
    for(i=0; i < np; i++){
      xdPr[i]	=dx*i;
      xdCr[i]	=dx*i;
      ESxDtPr[i]=xdPr[i];
      ESxDtCr[i]=xdCr[i];
      kdPr[i]	=1;
      kdCr[i]	=1;
    }
  }
  if(aProf != NULL){
    for(i=0; i < np; i++){
      xdPr[i]	=aProf[i];
      xdCr[i]	=aProf[i];
      ESxDtPr[i]=aProf[i];
      ESxDtCr[i]=aProf[i];
    }
    ESFirstSetSplDPr();
    ESFirstSetSplDCr();
  }
  else{
    if(A != NULL){
      double dx;
      dx	=1./ESnDtPr;
      for(i=0; i < np; i++){
	xdPr[i]	=dx*i;
	xdCr[i]	=dx*i;
	ESxDtPr[i]=xdPr[i];
	ESxDtCr[i]=xdCr[i];
      }
      ESFirstSetSplDPr();
      ESFirstSetSplDCr();
    }
  }
  A	=aProf;
  for(i=0; i < np; i++){
    ydPr[i]	=pProf[i]*cssp;
    ydCr[i]	=jProf[i]*cssj;
  }
  if(k == 6){
    for(i=0; i < np; i++){
      ydCr[i]	=1./jProf[i];
    }
  }
  return;
}
#endif
 
/**
 * Profile and geometry specification.  This procedure must be called
 * immediately prior to execution (Esc).
 *           - Profiles: both @a Pressure and @a Jparallel profiles are arrays
 *             of length 21 and are defined on the uniform sqrt(toroidal flux) mesh.
 *           - Geometry: @a Npv denotes the number of input points @a Rpv and @a Zpv.
 *                      - 4 <= Npv <= 12: the points must be ordered according to
 *                        left, right, top, bottom, then North-east, North-west,
 *                        South-West, South-East, North-east, etc.
 *                      - Npv > 12: boundary points must be adjacent (can go
 *                        clockwise or anticlockwise) with the last point
 *                        identical to the first Rpv[Npv-1] = Rpv[0], and
 *                        likewise for Zpv.
  *
 * Note: consider using escin_ for more flexibility in the input profiles.
 *
* @param Pressure array[21] for the pressure in [mu0 Pa] (input).
 * @param Jparallel array[21]  for <j.B>/(R0 <B.grad phi> in [T/m] (input); R0=magnetic axis.
 * @param Rpv array[Npv] for the plasma-vacuum boundary radius in [m] (input).
 * @param Zpv array[Npv] for the plasma-vacuum boundary elevation in [m] (input).
 * @param RBtor covariant toroidal magnetic field in vacuum [T m] (input).
 * @param Rext major radius estimate in [m] (input). A suitable Rext can
 *        reduce the Fourier spectrum.
 * @author L. E. Zakharov, A. Pletzer
 * @see Esc_p_q, escin_
 */
 
void Esc_p_jparallel(
	      double	*Pressure
	      ,double	*Jparallel
	      ,double	*Rpv
	      ,double	*Zpv
	      ,int	*Npv
	      ,double	*RBtor
	      ,double	*Rext
	      )
{
  int i;
 
 
  cssp	=1.0; /* EZcgm0; */
  cssj	=1.0; /* EZcgm0; */
  cssr	=1.;
  if(*Npv > 12){
    FlagPoints	=1;
    EcNjet	=*Npv;
    for(i=0; i < EcNjet; i++){
      EcRjet[i]	=Rpv[i];
      EcZjet[i]	=Zpv[i];
    }
  }
  else{
    FlagPoints	=0;
  }
 
    ESEqSolvRadC	=2;
    ESEqSolvInPr	=2;
    ESEqSolvInCr	=1;
    ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
 
  if(FlStart){
    ch2hp	=EZcr12*ESgt[1]*ESgt[1];
    crHp	=1./ESNp;
    npv		=*Npv < 13 ? *Npv : 12;
    mpv		=FlagPoints ? 8 : npv/2+1;
    ESSetPlVInt(npv,mpv);
 
    ESEqSolvFl	=0x31;
    ESEqSolvIt	=20;
    /*
    ESEqSolvRadC	=2;
    ESEqSolvInPr	=2;
    ESEqSolvInCr	=1;
    */
    np		=21;
    nj		=21;
    ESnDtPr	=np-1;
    ESnDtCr	=nj-1;
    /*
    ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
    */
    for(i=0; i < 21; i++){
      xdPr[i]	=0.05*i;
      xdCr[i]	=0.05*i;
      ESxDtPr[i]	=xdPr[i];
      ESxDtCr[i]	=xdCr[i];
    }
    for(i=0; i < ESNa1; i++){
      ESaF[i]	=ESRBt;
      ESaF1a[i]	=0;
      ESaF2a[i]	=0;
    }
  }
  ESRext	=*Rext;
  ESRBt		=*RBtor;
  for(i=0; i < npv; i++){
    rdPV[i]	=Rpv[i]*cssr;
    zdPV[i]	=Zpv[i]*cssr;
  }
  if(FlStart){
    if(FlagPoints){
      double zmx,zmn,rmx,rmn;
      zmx	=EcZjet[0];
      zmn	=zmx;
      rmx	=EcRjet[0];
      rmn	=rmx;
      for(i=1;  i < EcNjet; i++){
	if(zmx < EcZjet[i]){
	  zmx	=EcZjet[i];
	}
	if(zmn > EcZjet[i]){
	  zmn	=EcZjet[i];
	}
	if(rmx < EcRjet[i]){
	  rmx	=EcRjet[i];
	}
	if(rmn > EcRjet[i]){
	  rmn	=EcRjet[i];
	}
      }
      ESaR0[0]	=EZcr2*(0.9*rmn+1.1*rmx);
      ESaZ0[0]	=EZcr2*(zmx+zmn);
      rdPV[0]	=EZcr2*(rmn+rmx);
      zdPV[0]	=zmx;
      rdPV[1]	=EZcr2*(rmn+rmx);
      zdPV[1]	=zmn;
      rdPV[2]	=rmn;
      zdPV[2]	=EZcr2*(zmx+zmn);
      rdPV[3]	=rmx;
      zdPV[3]	=EZcr2*(zmx+zmn);
    }
    else{
      ESaR0[0]	=EZcr2*(0.9*Rpv[2]+1.1*Rpv[3]);
      ESaZ0[0]	=EZcr2*(Zpv[0]+Zpv[1]);
      if(ESiAx != 0){
	ESaZ0[0]	+=0.2*EZcr2*(Zpv[0]-Zpv[1]);
      }
    }
    ESInit2DPlGeom();
    EcInitPlVacData();
    ECRcs3D2Rcs2D(0);
    ECvr2r(0);
    ECvz2z(0);
    ESSwitchInpRadC();
    ESFirstSetSplDPr();
    ESFirstSetSplDCr();
    FlStart	=0;
  }
  for(i=0; i < np; i++){
    ydPr[i]	=Pressure[i]*cssp;
  }
  for(i=0; i < nj; i++){
    ydCr[i]	=Jparallel[i]*cssj;
  }
}
void esc_p_jparallel_(
	      double	*Pressure
	      ,double	*Jparallel
	      ,double	*Rpv
	      ,double	*Zpv
	      ,int	*Npv
	      ,double	*RBtor
	      ,double	*Rext
	      )
{
Esc_p_jparallel(
	      Pressure
	      , Jparallel
	      , Rpv
	      , Zpv
	      , Npv
	      , RBtor
	      , Rext
	      );
}
void esc_p_jparallel(
	      double	*Pressure
	      ,double	*Jparallel
	      ,double	*Rpv
	      ,double	*Zpv
	      ,int	*Npv
	      ,double	*RBtor
	      ,double	*Rext
	      )
{Esc_p_jparallel(
	      Pressure
	      , Jparallel
	      , Rpv
	      , Zpv
	      , Npv
	      , RBtor
	      , Rext
	      );
}
void ESC_P_JPARALLEL(
	      double	*Pressure
	      ,double	*Jparallel
	      ,double	*Rpv
	      ,double	*Zpv
	      ,int	*Npv
	      ,double	*RBtor
	      ,double	*Rext
	      )
{Esc_p_jparallel(
	      Pressure
	      , Jparallel
	      , Rpv
	      , Zpv
	      , Npv
	      , RBtor
	      , Rext
	      );
}
 
/**
 * Profile and geometry specification. This procedure must be called
 * immediately prior to execution (Esc).
 *           - Profiles: both @a Pressure and @a qprofile profiles are arrays
 *             of length 21 and are defined on the uniform sqrt(toroidal flux) mesh.
 *           - Geometry: @a Npv denotes the number of input points @a Rpv and @a Zpv.
 *                      - 4 <= Npv <= 12: the points must be ordered according to
 *                        left, right, top, bottom, then North-east, North-west,
 *                        South-West, South-East, North-east, etc.
 *                      - Npv > 12: boundary points must be adjacent (can go
 *                        clockwise or anticlockwise) with the last point
 *                        identical to the first Rpv[Npv-1] = Rpv[0], and
 *                        likewise for Zpv.
 *
 * Note: consider using escin_ for more flexibility in the input profiles.
 *
 * @param Pressure array[21] for the pressure in [mu0 Pa] (input).
 * @param qprofile array[21]  for the q profile (input).
 * @param Rpv array[Npv] for the plasma-vacuum boundary radius in [m] (input).
 * @param Zpv array[Npv] for the plasma-vacuum boundary elevation in [m] (input).
 * @param RBtor covariant toroidal magnetic field in vacuum [T m] (input).
 * @param Rext major radius estimate in [m] (input). A suitable Rext can
 *        reduce the Fourier spectrum.
 * @author L. E. Zakharov, A. Pletzer
 * @see Esc Esc_p_jparallel, esc_in
 */
 
void Esc_p_q(double	*Pressure
	      ,double	*qprofile
	      ,double	*Rpv
	      ,double	*Zpv
	      ,int	*Npv
	      ,double	*RBtor
	      ,double	*Rext
	      )
{
  int i;
 
 
  cssp	=1.0; /* EZcgm0; */
  cssj	=1.0; /* EZcgm0; */
  cssr	=1.;
  if(*Npv > 12){
    FlagPoints	=1;
    EcNjet	=*Npv;
    for(i=0; i < EcNjet; i++){
      EcRjet[i]	=Rpv[i];
      EcZjet[i]	=Zpv[i];
    }
  }
  else{
    FlagPoints	=0;
  }
 
    ESEqSolvRadC	=2;
    ESEqSolvInPr	=2;
    ESEqSolvInCr	=7;
    ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
 
  if(FlStart){
    ch2hp	=EZcr12*ESgt[1]*ESgt[1];
    crHp	=1./ESNp;
    npv		=*Npv < 13 ? *Npv : 12;
    mpv		=FlagPoints ? 8 : npv/2+1;
    ESSetPlVInt(npv,mpv);
 
    ESEqSolvFl	=0x31;
    ESEqSolvIt	=20;
    /*
    ESEqSolvRadC	=2;
	  ESEqSolvInPr	=2;
	  ESEqSolvInCr	=7;
	  */
    np		=21;
    nj		=21;
    ESnDtPr	=np-1;
    ESnDtCr	=nj-1;
    /*
    ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
	  */
    for(i=0; i < 21; i++){
      xdPr[i]	=0.05*i;
      xdCr[i]	=0.05*i;
      ESxDtPr[i]	=xdPr[i];
      ESxDtCr[i]	=xdCr[i];
    }
    for(i=0; i < ESNa1; i++){
      ESaF[i]	=ESRBt;
      ESaF1a[i]	=0;
      ESaF2a[i]	=0;
    }
  }
  ESRext	=*Rext;
  ESRBt		=*RBtor;
  for(i=0; i < npv; i++){
    rdPV[i]	=Rpv[i]*cssr;
    zdPV[i]	=Zpv[i]*cssr;
  }
  if(FlStart){
    if(FlagPoints){
      double zmx,zmn,rmx,rmn;
      zmx	=EcZjet[0];
      zmn	=zmx;
      rmx	=EcRjet[0];
      rmn	=rmx;
      for(i=1;  i < EcNjet; i++){
	if(zmx < EcZjet[i]){
	  zmx	=EcZjet[i];
	}
	if(zmn > EcZjet[i]){
	  zmn	=EcZjet[i];
	}
	if(rmx < EcRjet[i]){
	  rmx	=EcRjet[i];
	}
	if(rmn > EcRjet[i]){
	  rmn	=EcRjet[i];
	}
      }
      ESaR0[0]	=EZcr2*(0.9*rmn+1.1*rmx);
      ESaZ0[0]	=EZcr2*(zmx+zmn);
      rdPV[0]	=EZcr2*(rmn+rmx);
      zdPV[0]	=zmx;
      rdPV[1]	=EZcr2*(rmn+rmx);
      zdPV[1]	=zmn;
      rdPV[2]	=rmn;
      zdPV[2]	=EZcr2*(zmx+zmn);
      rdPV[3]	=rmx;
      zdPV[3]	=EZcr2*(zmx+zmn);
    }
    else{
      ESaR0[0]	=EZcr2*(0.9*Rpv[2]+1.1*Rpv[3]);
      ESaZ0[0]	=EZcr2*(Zpv[0]+Zpv[1]);
      if(ESiAx != 0){
	ESaZ0[0]	+=0.2*EZcr2*(Zpv[0]-Zpv[1]);
      }
    }
    ESInit2DPlGeom();
    EcInitPlVacData();
    ECRcs3D2Rcs2D(0);
    ECvr2r(0);
    ECvz2z(0);
    ESSwitchInpRadC();
    ESFirstSetSplDPr();
    ESFirstSetSplDCr();
    FlStart	=0;
  }
  for(i=0; i < np; i++){
    ydPr[i]	=Pressure[i]*cssp;
  }
  for(i=0; i < nj; i++){
    ydCr[i]	=1.0/qprofile[i];
  }
}
void esc_p_q_(double	*Pressure
	      ,double	*qprofile
	      ,double	*Rpv
	      ,double	*Zpv
	      ,int	*Npv
	      ,double	*RBtor
	      ,double	*Rext
	      )
{
Esc_p_q(
	     Pressure
	     , qprofile
	     , Rpv
	     , Zpv
	     , Npv
	     , RBtor
	     , Rext
	      );
}
void esc_p_q(double	*Pressure
	      ,double	*qprofile
	      ,double	*Rpv
	      ,double	*Zpv
	      ,int	*Npv
	      ,double	*RBtor
	      ,double	*Rext
	      )
{
Esc_p_q(
	     Pressure
	     , qprofile
	     , Rpv
	     , Zpv
	     , Npv
	     , RBtor
	     , Rext
	      );
}
void ESC_P_Q(double	*Pressure
	      ,double	*qprofile
	      ,double	*Rpv
	      ,double	*Zpv
	      ,int	*Npv
	      ,double	*RBtor
	      ,double	*Rext
	      )
{
Esc_p_q(
	     Pressure
	     , qprofile
	     , Rpv
	     , Zpv
	     , Npv
	     , RBtor
	     , Rext
	      );
}
 
/**
 * Profile recovery routine. This procedure can be invoked to recover
 * input profiles used in an equilibrium restored from a dump file
 * (<Machine.Runid>.es).
 *
 * @param Pprofile pressure-like profile (output). The exact type
 * depends on internal settings
 * @param Jprofile current-like profile (output).
 *
 * @author L. E. Zakharov
 */
 
int ESgetDataProfile(double *Pprofile, double *Jprofile)
{
  int i;
 
  for(i=0; i < 21; i++){
    Pprofile[i]	=ydPr[i]/cssp;
    Jprofile[i]	=ydCr[i]/cssj;
  }
  return(0);
}
void esgetdataprofile_(double *Pprofile, double *Jprofile)
{
  ESgetDataProfile(Pprofile,Jprofile);
}
void esgetdataprofile(double *Pprofile, double *Jprofile)
{
  ESgetDataProfile(Pprofile,Jprofile);
}
void ESGETDATAPROFILE(double *Pprofile, double *Jprofile)
{
  ESgetDataProfile(Pprofile,Jprofile);
}

int esc2orbit_()
{
  extern double *ESdgY,*ESdgY2a;
  static int Fic=0;
  int Na,Np;
  int Na1,Np1;
  FILE *lf;
  char FNm[24];
  int i,j,ji,k,kj,ki;
  double a,da,t,dt,gy,sq,sp,aF,qg,L0,rB;
  double dgy,s,cs,sn;
#ifdef XWIN
  double rt,zt,ra,za,B[Np+1],aJ[Np+1],L[Np+1],A[Na+1];
  double r[Np+1],z[Np+1];
  double rc[ESFp1],rs[ESFp1],rca[ESFp1],rsa[ESFp1],rct[ESFp1],rst[ESFp1],b,ba;
#else
  double rt,zt,ra,za,B[129],aJ[129],L[129],A[129];
  double r[128],z[128];
  double rc[33],rs[33],rca[33],rsa[33],rct[33],rst[33],b,ba;

  Na  =20; 
  Np  =63; 
  if(Na > 128){
    Na	=128;
  }
  if(Np > 128){
    Np	=128;
  }
#endif

  sprintf(FNm,"ES2Orbit.%2.2d",Fic);
  lf	=fopen(FNm,"w");
  fputs("      subroutine esc2orbit()\n",lf);
  fputs("      implicit none\n",lf);
  fputs("      include 'orbcom'\n",lf);
  fputs("      integer i,j\n",lf);
  fputs("      open(unit=9,file='ES2Orbit')\n",lf);
  fputs("      do i=1,24\n",lf);
  fputs("        read(9,*)\n",lf);
  fputs("      enddo\n",lf);
  fputs("      read(9,*)lsp,lst\n",lf);
  fputs("      read(9,*)\n",lf);
  fputs("      read(9,*)x1(1,1),z1(1,1),bax,i1(1,1)\n",lf);
  fputs("      read(9,*)\n",lf);
  fputs("      do i=1,lsp\n",lf);
  fputs("        read(9,*)psival(i),qd1(i),pd1(i),gd1(i)\n",lf);
  fputs("      enddo\n",lf);
  fputs("      pw=psival(lsp)\n",lf);
  fputs("      read(9,*)\n",lf);
  fputs("      do i=1,lsp\n",lf);
  fputs("        do j=1,lst\n",lf);
  fputs("          read(9,*)x1(j,i),z1(j,i),b1(j,i),i1(j,i)\n",lf);
  fputs("        enddo\n",lf);
  fputs("      enddo\n",lf);
  fputs("      close(9)\n",lf);
  fputs("      end\n",lf);

  Na1	=Na+1;
  Np1	=Np+1;
  fprintf(lf
	  ,"%4d %4d - Number of radial points (with both ends) and poloidal intervals\n"
	  ,Na1,Np);
  fprintf(lf,"%14s %14s %14s %14s %3s %3s\n"
	  ,"r [m]","z [m]","B0 [T]","I(a,t)/2pi MA","ia","it");
  r[0]	=ESaR0[0];
  z[0]	=ESaZ0[0];
  rB	=r[0]/ESaF[0];
  fprintf(lf,"%14.11f %14.11f %14.11f %14.7e %3d %3d - mag. axis\n"
	  ,r[0],z[0],ESaF[0]/r[0],0.,0,0);

  fprintf(lf,"%14s %14s %14s %14s\n","Psi/2pi/B0 m^2","q","2p/B^2_0"
	  ,"g=rB/B0 [m]");
  da	=ESsa[ESNa]/Na;
  dt	=EZc2gp/Np;
  qg	=ESaR0[0]/ESaF[0];
  qg	*=2.*qg;
  for(i=0; i < Na1; i++){
    gy	=sqrt(da*i);
    if(i){
      if(i < Na){
	do{
	  ESSetSplA(a);
	  splRA(&s,&dgy,ESqgY,ESqgY2a);
	  a		+=(gy-s)/dgy;
	}while(fabs(gy-s) > 1e-12);
      }
      else{
	a	=1.;
      }
    }
    else{
      a	=0.;
    }
    A[i]	=a;

    ESSetSplA(a);
    gy	=da*i*fabs(ESgY[ESNa]);
    splRA(&aF,NULL,ESaF,ESaF2a);
    if(ESEqSolvInPr != 2){
      splRA(&sp,NULL,ESsp,ESsp2a);
    }
    else{
      ESSetSplDPr(a);
      splRDPr(&sp,&s,2);
    }
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&sq,NULL,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&sq,&s,7);
    }
    sq	=1./sq;
    fprintf(lf,"%14.7e %14.7e %14.7e %14.11f\n",gy*rB,sq,sp*qg,aF*rB);
  }

  fprintf(lf,"%14s %14s %14s %14s %3s %3s\n"
	  ,"r [m]","z [m]","B/B0","I(a,t)/2pi/B0 m","ia","it");
  for(j=0; j < Np; j++){
    fprintf(lf,"%14.11f %14.11f %14.11f %14.7e %3d %3d\n"
	    ,r[0],z[0],1.,0.,0,j);
  }
  ji	=Np1;
  for(i=1; i < Na1; i++){
    a	=A[i];
    ESSetSplA(a);
    splRA(&dgy,NULL,ESdgY,ESdgY2a);
    dgy	*=a;
    splRA(&aF,NULL,ESaF,ESaF2a);
    splRA(&b,&ba,ESsb,ESsb2a);
    splRA(rc,rca,rcT,rcT2a);
    rc[0]	+=ESaR0[0];
    splRA(rs,rsa,rsT,rsT2a);
    rs[0]	+=ESaZ0[0];
    ki		=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,rca+k,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      rca[k]	*=2.;
      rst[k]	=-rc[k]*k;
      splRA(rs+k,rsa+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      rct[k]	=rs[k]*k;
      rsa[k]	*=2.;
    }
    L0	=0.;
    for(j=0; j < Np; j++){
      t		=-dt*j;
      cs	=cos(t);
      sn	=sin(t);
      r[j]	=rc[0]+rc[1]*cs+rs[1]*sn;
      ra	=rca[0]+rca[1]*cs+rsa[1]*sn;
      rt	=rct[1]*cs+rst[1]*sn;
      z[j]	=rs[0]+b*sn;
      za	=rsa[0]+ba*sn;
      zt	=b*cs;
      for(k=2; k < ESFp1; k++){
	s	=k*t;
	cs	=cos(s);
	sn	=sin(s);
	r[j]	+=rc[k]*cs+rs[k]*sn;
	ra	+=rca[k]*cs+rsa[k]*sn;
	rt	+=rct[k]*cs+rst[k]*sn;
      }
      qg	=ra*zt-rt*za;
      L[j]	=qg/r[j];
      L0	+=L[j];
      qg	*=r[j];
      s		=aF/r[j];
      ra	=rt*rt+zt*zt;
      za	=-dgy/qg;
      aJ[j]	=ra*za;
      B[j]	=sqrt(ra*za*za+s*s);
      ji++;
    }
    L0	/=Np;
    s	=-aF*aF/dgy;
    for(j=0; j < Np; j++){
      aJ[j]	+=s*(L[j]-L0);
      fprintf(lf,"%14.11f %14.11f %14.11f %14.7e %3d %3d\n"
	      ,r[j],z[j],B[j]*rB,aJ[j]*rB,i,j);
    }
    ji++;
    B[j]	=B[0];
    aJ[j]	=aJ[0];
  }
  fclose(lf);
  printf("%s has been written\n",FNm);
  Fic++;
  return(0);
}

int esc2orbitB()	/* File for Orbit in Boozer coordinates */
{
  extern double *ESdgY,*ESdgY2a,*ESLs,*ESLs2,*ESg22s,*ESg22s2;
  static int Fic=0;
  int Na,Np,Na1,Np1;
  FILE *lf;
  char FNm[24];

  int i,j,ji,k,kj,ki;
  double a,da,t,dt,gy,sq,sp,aJ,aF,qg,rB;
  double dgy,s,cs,sn,s1,gh,qg0;

#ifdef XWIN
  double r,z,rt,zt,ra,za,B[Np+1],A[Na+1];
  double Lc[ES2Mp1],Ls[ES2Mp1]
    ,rc[ESFp1],rs[ESFp1],rca[ESFp1],rsa[ESFp1],rct[ESFp1],rst[ESFp1],b,ba;
  double Gh[ESNp1],Gh2[ESNp1];
#else
  double r,z,rt,zt,ra,za,B[129],A[129];
  double Lc[65],Ls[65],rc[33],rs[33],rca[33],rsa[33],rct[33],rst[33],b,ba;
  double Gh[129],Gh2[129];
#endif
  Na	=20;
  Np	=63;

  sprintf(FNm,"ES2OrbitB.%2.2d",Fic);
  lf	=fopen(FNm,"w");
  fputs("      subroutine esc2orbit()\n",lf);
  fputs("      implicit none\n",lf);
  fputs("      include 'orbcom'\n",lf);
  fputs("      integer i,j\n",lf);
  fputs("      open(unit=9,file='ES2Orbit')\n",lf);
  fputs("      do i=1,24\n",lf);
  fputs("        read(9,*)\n",lf);
  fputs("      enddo\n",lf);
  fputs("      read(9,*)lsp,lst\n",lf);
  fputs("      read(9,*)\n",lf);
  fputs("      read(9,*)x1(1,1),z1(1,1),bax,i1(1,1)\n",lf);
  fputs("      read(9,*)\n",lf);
  fputs("      do i=1,lsp\n",lf);
  fputs("        read(9,*)psival(i),qd1(i),pd1(i),gd1(i)\n",lf);
  fputs("      enddo\n",lf);
  fputs("      pw=psival(lsp)\n",lf);
  fputs("      read(9,*)\n",lf);
  fputs("      do i=1,lsp\n",lf);
  fputs("        do j=1,lst\n",lf);
  fputs("          read(9,*)x1(j,i),z1(j,i),b1(j,i),i1(j,i)\n",lf);
  fputs("        enddo\n",lf);
  fputs("      enddo\n",lf);
  fputs("      close(9)\n",lf);
  fputs("      end\n",lf);

  Na1	=Na+1;
  Np1	=Np+1;
  fprintf(lf
	  ,"%4d %4d - Number of radial points (with both ends) and poloidal intervals\n"
	  ,Na1,Np);
  fprintf(lf,"%14s %14s %14s %14s %3s %3s\n"
	  ,"r [m]","z [m]","B0 [T]","I(a,t)/2pi MA","ia","it");
  r	=ESaR0[0];
  z	=ESaZ0[0];
  rB	=r/ESaF[0];
  fprintf(lf,"%14.11f %14.11f %14.11f %14.7e %3d %3d - mag. axis\n"
	  ,r,z,ESaF[0]/r,0.,0,0);

  fprintf(lf,"%14s %14s %14s %14s\n","Psi/2pi/B0 m^2","q","2p/B^2_0"
	  ,"g=rB/B0 [m]");
  da	=ESsa[ESNa]/Na;
  dt	=EZc2gp/Np;
  qg	=ESaR0[0]/ESaF[0];
  qg	*=2.*qg;
  for(i=0; i < Na1; i++){
    gy	=sqrt(da*i);
    if(i){
      if(i < Na){
	do{
	  ESSetSplA(a);
	  splRA(&s,&dgy,ESqgY,ESqgY2a);
	  a		+=(gy-s)/dgy;
	}while(fabs(gy-s) > 1e-12);
      }
      else{
	a	=1.;
      }
    }
    else{
      a	=0.;
    }
    A[i]	=a;
    ESSetSplA(a);
    gy	=da*i*fabs(ESgY[ESNa]);
    splRA(&aF,NULL,ESaF,ESaF2a);
    if(ESEqSolvInPr != 2){
      splRA(&sp,NULL,ESsp,ESsp2a);
    }
    else{
      ESSetSplDPr(a);
      splRDPr(&sp,&s,2);
    }
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&sq,NULL,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&sq,&s,7);
    }
    sq	=1./sq;
    fprintf(lf,"%14.7e %14.7e %14.7e %14.11f\n",gy*rB,sq,sp*qg,aF*rB);
  }

  fprintf(lf,"%14s %14s %14s %14s %3s %3s\n"
	  ,"r [m]","z [m]","B/B0","I(a,t)/2pi/B0 m","ia","it");
  for(j=0; j < Np; j++){
    fprintf(lf,"%14.11f %14.11f %14.11f %14.7e %3d %3d\n"
	    ,r,z,1.,0.,0,j);
  }
  ji	=Np1;
  for(i=1; i < Na1; i++){
    a	=A[i];
    ESSetSplA(a);
    splRA(&dgy,NULL,ESdgY,ESdgY2a);
    dgy	*=a;
    splRA(&aJ,NULL,ESaJ,ESaJ2a);
    splRA(&aF,NULL,ESaF,ESaF2a);
    s1	=dgy/aF;
    s1	*=s1*ESaR0[0];
    splRA(Lc,NULL,ESLc,ESLc2);
    splRA(&s,NULL,ESg22c,ESg22c2);
    Lc[0]	=2./(Lc[0]+s*s1);
    ki	=0;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      splRA(Lc+k,NULL,ESLc+ki,ESLc2+ki);
      splRA(Ls+k,NULL,ESLs+ki,ESLs2+ki);
      splRA(&s,NULL,ESg22c+ki,ESg22c2+ki);
      Lc[k]	=(Lc[k]+s*s1)*Lc[0];
      splRA(&s,NULL,ESg22s+ki,ESg22s2+ki);
      Ls[k]	=(Ls[k]+s*s1)*Lc[0];
    }
    for(j=0; j < ESNp; j++){
      s		=0.;
      kj	=0;
      for(k=1; k < ES2Mp1; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	s	+=(Lc[k]*ESsn1[kj]+Ls[k]*(1.-EScs1[kj]))/k;
      }
      Gh[j]	=s;
    }
    Gh[j]	=Gh[0];
    splP(Gh,Gh2);

    splRA(&b,&ba,ESsb,ESsb2a);
    splRA(rc,rca,rcT,rcT2a);
    rc[0]	+=ESaR0[0];
    splRA(rs,rsa,rsT,rsT2a);
    rs[0]	+=ESaZ0[0];
    ki	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,rca+k,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      rca[k]	*=2.;
      rst[k]	=-rc[k]*k;
      splRA(rs+k,rsa+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      rct[k]	=rs[k]*k;
      rsa[k]	*=2.;
    }
    t	=EZc2gp;
    qg0	=0.;
    for(j=0; j < Np; j++){
      gh	=EZc2gp-dt*j;
      do{
	ESSetSplP(t);
	splRP(&s,&sp,Gh,Gh2);
	s      +=t-gh;
	t	-=s/(1.+sp);
      }while(fabs(s) > 1e-8);
      s		=t;
      cs	=cos(s);
      sn	=sin(s);

      r		=rc[0]+rc[1]*cs+rs[1]*sn;
      ra	=rca[0]+rca[1]*cs+rsa[1]*sn;
      rt	=rct[1]*cs+rst[1]*sn;
      z		=rs[0]+b*sn;
      za	=rsa[0]+ba*sn;
      zt	=b*cs;
      for(k=2; k < ESFp1; k++){
	s	=k*t;
	cs	=cos(s);
	sn	=sin(s);
	r	+=rc[k]*cs+rs[k]*sn;
	ra	+=rca[k]*cs+rsa[k]*sn;
	rt	+=rct[k]*cs+rst[k]*sn;
      }
      qg	=r*(ra*zt-rt*za);
      s		=aF/r;
      B[j]	=(rt*rt+zt*zt)*dgy*dgy/(qg*qg)+s*s;
      qg0	+=1./B[j];
      fprintf(lf,"%14.11f %14.11f %14.11f %14.7e %3d %3d\n"
	      ,r,z,sqrt(B[j])*rB,aJ*rB,i,j);
      ji++;
    }
    ji++;
  }
  fclose(lf);
#ifdef XWIN
  printf("%s has been written\n",FNm);
#endif
  Fic++;
  return(0);
}

void ESC2ORBIT(){esc2orbit_();}
void esc2orbit(){esc2orbit_();}

void ESC2ORBITB(){esc2orbitB();}
void esc2orbitb_(){esc2orbitB();}
void esc2orbitb(){esc2orbitB();}


/*********  Autonomous Routines for Real Space Equilibrium ********/
static int Na,Na1,Np,Np1,MemFl=0;
static double *sa,*gq,*aF,*daF,*dgF,*d2gF,*dgY,*d2gY,*aT,*daT,*aP,*daP;
static double *sr,*sra,*srq,*sraq;
static double *sz,*sza,*szq,*szaq;
static double *aB,*aBa,*aBq,*aBaq;
static double *gH,*gHa,*gHq,*gHaq;

static int i0=0,i1=1,j00,j10,j01,j11; 
static double A,ha,rha,hq,rhq,cgq0,crgq0;/* period and inverse period in gq */
static double xF,XF,xD,XD,dX,dxD,dXD,yF,YF,yD,YD,dY,dyD,dYD;

int ESFreeRZspaceEq()
{
  free(gH);
  free(aB);
  free(sz);
  free(sr);
  free(gq);
  MemFl=0;
  return(0);
}

int ESWriteESI0()
{
  static int Fic=0;
  int i,j,ji;
  double *f=NULL,*fa,*fq,*faq,rr0;
  FILE *lf,*lfa;
  char FNm[16];

#ifdef XWIN
  double Te,Ti,p,dp,t;
  extern char *CbUserMessage;
  extern char ESMessage[];
  CbUserMessage	=ESMessage;
#endif
  f	=(double*)malloc(4*ESNp1*sizeof(double));
  if(f == NULL){
#ifdef XWIN
    sprintf(ESMessage,"No memory for f");
#endif
    return(1);
  }
  fa	=f+ESNp1;
  fq	=fa+ESNp1;
  faq	=fq+ESNp1;
  sprintf(FNm,"esiBin.00");
  lf	=fopen(FNm,"w");
  if(lf == NULL){
#ifdef XWIN
    sprintf(ESMessage,"%s cannot be open",FNm);
#endif
    free(f);
    return(1);
  }
  sprintf(FNm,"esiA.%2.2d",Fic);
  sprintf(FNm,"esiA.00");
  lfa	=fopen(FNm,"w");
  if(lfa == NULL){
#ifdef XWIN
    sprintf(ESMessage,"%s cannot be open",FNm);
#endif
    fclose(lf);
    free(f);
    return(1);
  }

#ifdef XWIN
  time(&stTime);
  ltm	=localtime(&stTime);
  sprintf(ESMessage,"Date: %2.2d/%2.2d/%2.2d at %2.2d:%2.2d"
	  ,ltm->tm_mon+1,ltm->tm_mday,ltm->tm_year%100,ltm->tm_hour
	  ,ltm->tm_min); 
  if(ESIWriteDensity(ESNa,Fic)){
    fclose(lfa);
    return(1);
  }
#endif

  fprintf(lfa,"!!! Do not edit this file\n");
  fprintf(lfa,"%3d x %3d - numbers of poloidal x radial data points\n"
	  ,ESNp1,ESNa1);
  fwrite(&ESNp1,sizeof(int),1,lf);
  fwrite(&ESNa1,sizeof(int),1,lf);
  fputs("gq:\n",lfa);
  for(j=0; j < ESNp1; j++){
    fprintf(lfa,"%24.16e",ESgt[j]);
    if(j%4 == 3){
      fputc('\n',lfa);
    }
  }
  if(j%4 != 0){
    fputc('\n',lfa);
  }
  fwrite(ESgt,sizeof(double),ESNp1,lf);
  fprintf(lfa,"%24s%24s%24s\n","a","F","Fa");
  fwrite(ESsa,sizeof(double),ESNa1,lf);
  fwrite(ESaF,sizeof(double),ESNa1,lf);
  fwrite(ESaF1a,sizeof(double),ESNa1,lf);
  rr0	=1./ESaR0[0];
  for(i=0; i < ESNa1; i++){
    f[i]	=ESLc[i]*ESaF[i]*rr0;
    fq[i]	=(ESLc1[i]*ESaF[i]+ESLc[i]*ESaF1a[i])*rr0;
    fprintf(lfa,"%24.16e%24.16e%24.16e\n",ESsa[i],ESaF[i],ESaF1a[i]);
  }
  fprintf(lfa,"%24s%24s%24s%24s\n","gF'/a","gF''","gY'/a","(gY'/a)'");
  fwrite(f,sizeof(double),ESNa1,lf);
  fwrite(fq,sizeof(double),ESNa1,lf);
  fwrite(ESdgY,sizeof(double),ESNa1,lf);
  fwrite(ESdgY1a,sizeof(double),ESNa1,lf);
  for(i=0; i < ESNa1; i++){
    fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	    ,f[i],fq[i],ESdgY[i],ESdgY1a[i]);
  }
  fprintf(lfa,"%24s%24s%24s%24s\n","T","Ta","P","Pa");
  fwrite(ESaT,sizeof(double),ESNa1,lf);
  splA1a(ESaT,fa,ESaT2a);
  fwrite(fa,sizeof(double),ESNa1,lf);
  fwrite(ESPs,sizeof(double),ESNa1,lf);
  fwrite(ESPs1a,sizeof(double),ESNa1,lf);
  for(i=0; i < ESNa1; i++){
    fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	    ,ESaT[i],ESaT1a[i],ESPs[i],ESPs1a[i]);
  }

  fprintf(lfa,"%24s%24s%24s%24s\n","r","r'_a","r'_gq","r''_{a,gq}");
  for(i=0; i < ESNa1; i++){
    ji	=ESNp1*i+ESNp1;
    for(j=0; j < ESNp1; j++){
      ji--;
      fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	      ,ESsr[ji],ESsra[ji],-ESsrt[ji],-ESsrat[ji]);
    }
  }

  fprintf(lfa,"%24s%24s%24s%24s\n","z","z'_a","z'_gq","z''_{a,gq}");
  for(i=0; i < ESNa1; i++){
    ji	=ESNp1*i+ESNp1;
    for(j=0; j < ESNp1; j++){
      ji--;
      fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	      ,ESsz[ji],ESsza[ji],-ESszt[ji],-ESszat[ji]);
    }
  }

  fprintf(lfa,"%24s%24s%24s%24s\n","B","B'_a","B'_gq","B''_{a,gq}");
  for(i=0; i < ESNa1; i++){
    ji	=ESNp1*i+ESNp1;
    for(j=0; j < ESNp1; j++){
      ji--;
      fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	      ,ESaB[ji],ESaBa[ji],-ESaBt[ji],-ESaBat[ji]);
    }
  }

  fprintf(lfa,"%24s%24s%24s%24s\n"
	  ,"gh'_gq","gh''_{a,gq}","gh''_{gq,gq}","gh'''_{a,gq,gq}");
  for(i=0; i < ESNa1; i++){
    ji	=ESNp1*i+ESNp1;
    for(j=0; j < ESNp1; j++){
      ji--;
      fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	      ,ESgH[ji],ESgHa[ji],-ESgHt[ji],-ESgHat[ji]);
    }
  }
  for(i=0; i < ESNa1; i++){
    ji	=ESNp1*i;
    j	=ESNp1;
    while(j > 0){
      j--;
      f[j]	=ESsr[ji];
      fa[j]	=ESsra[ji];
      fq[j]	=-ESsrt[ji];
      faq[j]	=-ESsrat[ji];
      ji++;
    }
    fwrite(f,sizeof(double),4*ESNp1,lf);
    ji	=ESNp1*i;
    j	=ESNp1;
    while(j > 0){
      j--;
      f[j]	=ESsz[ji];
      fa[j]	=ESsza[ji];
      fq[j]	=-ESszt[ji];
      faq[j]	=-ESszat[ji];
      ji++;
    }
    fwrite(f,sizeof(double),4*ESNp1,lf);
    ji	=ESNp1*i;
    j	=ESNp1;
    while(j > 0){
      j--;
      f[j]	=ESaB[ji];
      fa[j]	=ESaBa[ji];
      fq[j]	=-ESaBt[ji];
      faq[j]	=-ESaBat[ji];
      ji++;
    }
    fwrite(f,sizeof(double),4*ESNp1,lf);
    ji	=ESNp1*i;
    j	=ESNp1;
    while(j > 0){
      j--;
      f[j]	=ESgH[ji];
      fa[j]	=ESgHa[ji];
      fq[j]	=-ESgHt[ji];
      faq[j]	=-ESgHat[ji];
      ji++;
    }
    fwrite(f,sizeof(double),4*ESNp1,lf);
  }
  free(f);

#ifdef XWIN
  Te	=15.;
  Ti	=15.;
  t	=1./(1.60022e-2*EZcgm0*(Te+Ti));

  {
    double F,f0,f1;
    F	=0.;
    f1	=(ESLc[0]+ESVc[0])*ESsp[0];
    for(i=1; i < ESNa1; i++){
      f0=f1;
      f1=(ESVc[i]+ESLc[i])*ESsp[i];
      F	+=EZcr4*(f0+f1)*(ESpa[i]-ESpa[i-1]);
    }
    F	*=EZc2gp*EZc2gp*ESaR0[0]*t;
    fprintf(lfa,"Te=Ti=15keV N=%8.2e a%24s%24s\n",F*1e+20
	    ,"ne [10^20/m^3]","dne/da [10^20/m^3]");
  }

  ji	=0;
  printf("%2s %8s %8s %12s %12s %12s %12s %12s\n","i"
	 ,"a","r [m]","dr/da","B [T]","dB/da","n[10^20/m^3]","dn/da");
  for(i=0; i < ESNa1; i++){
    fprintf(lfa,"%24.16e%24.16e%24.16e\n",ESsa[i],t*ESsp[i],t*ESsp1a[i]);
    p	=i ? 1./ESsra[ji] : 0.;
    printf("%2d %8.6f %8.6f %12.5e %12.5e %12.5e %12.5e %12.5e\n",i,ESsa[i]
	   ,ESsr[ji],ESsra[ji],ESaB[ji],ESaBa[ji]*p,t*ESsp[i],t*ESsp1a[i]*p);
    ji	+=ESNp1;
  }
  sprintf(ESMessage,"%s has been written",FNm);
#endif
  fclose(lfa);
  fclose(lf);
  Fic++;
#ifdef XWIN
  if(esiread_(FNm)){
    sprintf(ESMessage,"%s Bad file",FNm);
  }
#endif
  return(0);
}

int ESWriteESI(char *CoordNm,double *sr,double*sz,double*aB,double*gH
	       ,int Na,int Np)
{
  static int Fic[3]={0,0,0},ic=0;
  int Na1,Np1;
  int i,j,ji;
  double *srt,*sra,*srat;
  double *szt,*sza,*szat;
  double *aBt,*aBa,*aBat;
  double *gHt,*gHa,*gHat;
  double rr0,a,da,t,dt;
  double F,Fa,dgY,dgYa,L,La,V,Va;

  double Te,Ti,p,dp;

  FILE *lfa;
  char FNm[16];
  extern char ESMessage[];
#ifdef XWIN
  extern char *CbUserMessage;
  CbUserMessage	=ESMessage;
#endif
  Na1	=Na+1;
  Np1	=Np+1;
  da	=(ESsa[ESNa]-ESsa[0])/Na;
  dt	=EZc2gp/Np;

  i	=Np1*Na1;
  sra	=sr	+i;
  srt	=sra	+i;
  srat	=srt	+i;
  sza	=sz	+i;
  szt	=sza	+i;
  szat	=szt	+i;
  aBa	=aB	+i;
  aBt	=aBa	+i;
  aBat	=aBt	+i;
  if(gH != NULL){
    gHa	=gH	+i;
    gHt	=gHa	+i;
    gHat=gHt	+i;
  }
  *ESMessage	=CoordNm != NULL ? *CoordNm : 'A';
  switch(*ESMessage){
  case 'P':
    ic	=1;
    break;
  case 'B':
    ic	=2;
    break;
  default:
    ic	=0;
    break;
  }

  sprintf(FNm,"esi%c.%2.2d",*ESMessage,Fic[ic]);
  sprintf(FNm,"esi%c.00",*ESMessage);
  *ESMessage	='\0';

  lfa	=fopen(FNm,"w");
  if(lfa == NULL){
    sprintf(ESMessage,"%s cannot be open",FNm);
    return(1);
  }

#ifdef XWIN
  time(&stTime);
  ltm	=localtime(&stTime);
  sprintf(ESMessage,"Date: %2.2d/%2.2d/%2.2d at %2.2d:%2.2d"
	  ,ltm->tm_mon+1,ltm->tm_mday,ltm->tm_year%100,ltm->tm_hour
	  ,ltm->tm_min); 
#endif
  fprintf(lfa,"!!! Do not edit this file. %s\n",ESMessage);

  if(CoordNm != NULL){
    fprintf(lfa,"%3d x %3d %s - number of poloidal x radial data points.\n"
	    ,Np1,Na1,CoordNm);
  }
  else{
    fprintf(lfa,"%3d x %3d ESC - numbers of poloidal x radial data points.\n"
	    ,Np1,Na1);
  }

  fputs("gq:\n",lfa);
  for(j=0; j < Np1; j++){
    fprintf(lfa,"%24.16e",dt*j);
    if(j%4 == 3){
      fputc('\n',lfa);
    }
  }
  if(j%4 != 0){
    fputc('\n',lfa);
  }
  fprintf(lfa,"%24s%24s%24s\n","a","F","Fa");
  rr0	=1./ESaR0[0];
  for(i=0; i < Na1; i++){
    a	=da*i;
    ESSetSplA(a);
    splRA(&F,&Fa,ESFF,ESFF2a);
    F	=sqrt(F);
    Fa	/=2.*F;
    fprintf(lfa,"%24.16e%24.16e%24.16e\n",da*i,F,Fa);
  }
  fprintf(lfa,"%24s%24s%24s%24s\n","gF'/a","gF''","gY'/a","(gY'/a)'");
  for(i=0; i < Na1; i++){
    a	=da*i;
    ESSetSplA(a);
    splRA(&F,&Fa,ESFF,ESFF2a);
    splRA(&dgY,&dgYa,ESdgY,ESdgY2a);
    F	=sqrt(F);
    Fa	/=2.*F;
    splRA(&L,&La,ESLc,ESLc2);
    fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	    ,L*F*rr0,(La*F+L*Fa)*rr0,dgY,dgYa);
  }
  fprintf(lfa,"%24s%24s%24s%24s\n","T","Ta","P","Pa");
  for(i=0; i < Na1; i++){
    a	=da*i;
    ESSetSplA(a);
    splRA(&F,&Fa,ESaT,ESaT2a);
    splRA(&dgY,&dgYa,ESPs,ESPs2a);
    fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n",F,Fa,dgY,dgYa);
  }

  fprintf(lfa,"%24s%24s%24s%24s\n","r","r'_a","r'_gq","r''_{a,gq}");
  for(i=0; i < Na1; i++){
    ji	=Np1*i+Np1;
    for(j=0; j < Np1; j++){
      ji--;
      fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	      ,sr[ji],sra[ji],-srt[ji],-srat[ji]);
    }
  }

  fprintf(lfa,"%24s%24s%24s%24s\n","z","z'_a","z'_gq","z''_{a,gq}");
  for(i=0; i < Na1; i++){
    ji	=Np1*i+Np1;
    for(j=0; j < Np1; j++){
      ji--;
      fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	      ,sz[ji],sza[ji],-szt[ji],-szat[ji]);
    }
  }

  fprintf(lfa,"%24s%24s%24s%24s\n","B","B'_a","B'_gq","B''_{a,gq}");
  for(i=0; i < Na1; i++){
    ji	=Np1*i+Np1;
    for(j=0; j < Np1; j++){
      ji--;
      fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	      ,aB[ji],aBa[ji],-aBt[ji],-aBat[ji]);
    }
  }

  fprintf(lfa,"%24s%24s%24s%24s\n"
	  ,"gh'_gq","gh''_{a,gq}","gh''_{gq,gq}","gh'''_{a,gq,gq}");
  if(gH != NULL){
    for(i=0; i < Na1; i++){
      ji	=Np1*i+Np1;
      for(j=0; j < Np1; j++){
	ji--;
	fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
		,gH[ji],gHa[ji],-gHt[ji],-gHat[ji]);
      }
    }
  }
  else{
    for(i=0; i < Na1; i++){
      for(j=0; j < Np1; j++){
	fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n",0.,0.,0.,0.);
      }
    }
  }

  Te	=15.;
  Ti	=15.;
  t	=1./(1.60022e-2*EZcgm0*(Te+Ti));

  {
    double f0,f1,pa;
    F	=0.;
    f1	=(ESLc[0]+ESVc[0])*ESsp[0];
    pa	=0.;
    for(i=1; i < Na1; i++){
      a	=da*i;
      ESSetSplA(a);
      splRA(&L,&La,ESLc,ESLc2);
      splRA(&V,&Va,ESVc,ESVc2);
      splRA(&p,&dp,ESsp,ESsp2a);
      f0=f1;
      V	+=L;
      f1=V*p;
      F	+=EZcr4*(f0+f1)*(a*a-pa);
      pa	=a*a;
    }
    F	*=EZc2gp*EZc2gp*ESaR0[0]*t;
  }

  fprintf(lfa,"Te=Ti=15keV N=%8.2e a%24s%24s\n",F*1e+20
	  ,"ne [10^20/m^3]","dne/da [10^20/m^3]");
  for(i=0; i < Na1; i++){
    a	=da*i;
    ESSetSplA(a);
    splRA(&p,&dp,ESsp,ESsp2a);
    fprintf(lfa,"%24.16e%24.16e%24.16e\n",a,t*p,t*dp);
  }
  fclose(lfa);
  sprintf(ESMessage,"%s has been written",FNm);
#ifdef XWIN
  if(esiread_(FNm)){
    sprintf(ESMessage,"%s Bad file",FNm);
  }
#endif
  Fic[ic]++;
  return(0);
}
#endif

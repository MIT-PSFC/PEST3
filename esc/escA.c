#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "esc.h"

#ifndef mzl_ESmain
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern unsigned int ESmemlev;
extern int EZspl_p(double *g,double *d2g,double *x,int n);
#endif

#ifndef mzl_1Dcore
extern int ESNa,ESNa1;
extern double *ESsa,*ESpa;
#endif

#ifndef mzl_2Dcore
extern int ESNp,ESNp1,ESnAP;
extern int *ESnF,*ESmF,*ESkF,ESNf,ESNf1;
extern int ESFp,ESFp1,ESnAF;
static double rc[64],rs[64],rc1a[64],rs1a[64],rc1t[64],rs1t[64];

extern double *ESgh,*ESgt,*EScs1,*ESsn1,*EZcs2,*EZsn2;
extern double *ECr,*ECz;
#endif

#ifndef mzl_3Dcore
extern int ESLt,ESNt,ESNt1,ESnAT,ESnAPT,ESLt0,ESNLt,ESNLt1;
extern double *rcT,*rcT1a,*rcT2a,*rcT1t,*rsT,*rsT1a,*rsT2a,*rsT1t;
extern double *drcT,*drcT1a,*drcT2a,*drsT,*drsT1a,*drsT2a;

extern double *ESaR0,*ESaZ0,*ESsr,*ESsz,*ESsb,*ESsb1a,*ESsb2a;

extern double *gf1T,*csLT,*snLT,*cs1T,*sn1T;
extern double *EStgF,*EStgF1a,*EStgF3a,ESgFPlV,ESgFPlV1a;
extern double *EStgY,*EStgY1a,*EStgY1a2t,*EStgY3a,*ESgi,*ESgi1a,*ESgi2a;
extern double *ESgF0,*ESgF01a,*ESgY0,*ESgY01a,*ESgY02a;
extern double ESgimn,ESgimx;
extern double *ESaY,*ESaY1a,*ESaY2a;
#endif
extern double *ESFaux,*ESd2Faux;

#ifndef mzl_ESPlV
static double *R0T2t,*Z0T2t,*bT2t; 
static double *gX0,*gX01a,*dR0T,*dZ0T,*dbT,*dbT1a,*dbT2a;

static double *DP,*DP2a,*DP2t,*DP2p,*DP22at,*DP22ap,*DP22tp,*DP222atp;
static double *gX,*gX2a,*gX2t,*gX2p,*gX22at,*gX22ap,*gX22tp,*gX222atp;
static double *BV,*BV2t,*BA,*BA2p,*BF,*BF2p,*BP,*BP2p;
static double *GH,*AA22at,*AA22ap,*AA22tp,*AA222atp;
static double *AF,*AF2t,*AF2p,*AF22at,*AF22tp,*AF222atp;
static double *AP,*AP2a,*AP2t,*AP2p,*AP22ap,*AP22tp;

#ifdef S
static double cha,crha,chha,ch2ha,ch2h0,ch2h1,ch4h0,ch4h1;
#endif

static double chp,crhp,chhp,ch2hp,crHp;
static double cht,crht,chht,ch2ht,crHt;

static int NnMN=0,nMN,nM,nN;
static double *gXc,*gXc2a,*gXs,*gXs2a;
static double *gYRc,*gYRc2a,*gYRs,*gYRs2a;
static int *MgYR,*NgYR,NngYR=0,ngYR=0;

static double *glT=NULL;

static int NaNpNt=0;
static double *ginf=NULL,*g2inf=NULL;

int ESReInit3APTSpl()
{
  if(ReInitArray((void**)&DP,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)DP,"ESReInit3APTSpl() - no memory for DP");
    ESexit(0);
  }
  if(ReInitArray((void**)&DP2a,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)DP2a,"ESReInit3APTSpl() - no memory for DP2a");
    ESexit(0);
  }
  if(ReInitArray((void**)&DP2t,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)DP2t,"ESReInit3APTSpl() - no memory for DP2t");
    ESexit(0);
  }
  if(ReInitArray((void**)&DP2p,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)DP2p,"ESReInit3APTSpl() - no memory for DP2p");
    ESexit(0);
  }
  if(ReInitArray((void**)&DP22at,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)DP22at,"ESReInit3APTSpl() - no memory for DP22at");
    ESexit(0);
  }
  if(ReInitArray((void**)&DP22ap,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)DP22ap,"ESReInit3APTSpl() - no memory for DP22ap");
    ESexit(0);
  }
  if(ReInitArray((void**)&DP22tp,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)DP22tp,"ESReInit3APTSpl() - no memory for DP22tp");
    ESexit(0);
  }
  if(ReInitArray((void**)&DP222atp,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)DP222atp,"ESReInit3APTSpl()-no memory for DP222atp");
    ESexit(0);
  }
  
  if(ReInitArray((void**)&gX,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)gX,"ESReInit3APTSpl() - no memory for gX");
    ESexit(0);
  }
  if(ReInitArray((void**)&gX2a,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)gX2a,"ESReInit3APTSpl() - no memory for gX2a");
    ESexit(0);
  }
  if(ReInitArray((void**)&gX2t,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)gX2t,"ESReInit3APTSpl() - no memory for gX2t");
    ESexit(0);
  }
  if(ReInitArray((void**)&gX2p,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)gX2p,"ESReInit3APTSpl() - no memory for gX2p");
    ESexit(0);
  }
  if(ReInitArray((void**)&gX22at,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)gX22at,"ESReInit3APTSpl() - no memory for gX22at");
    ESexit(0);
  }
  if(ReInitArray((void**)&gX22ap,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)gX22ap,"ESReInit3APTSpl() - no memory for gX22ap");
    ESexit(0);
  }
  if(ReInitArray((void**)&gX22tp,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)gX22tp,"ESReInit3APTSpl() - no memory for gX22tp");
    ESexit(0);
  }
  if(ReInitArray((void**)&gX222atp,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)gX222atp,"ESReInit3APTSpl()-no memory for gX222atp");
    ESexit(0);
  }

  if(ReInitArray((void**)&BV,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)BV,"ESReInit3APTSpl() - no memory for BV");
    ESexit(0);
  }
  if(ReInitArray((void**)&GH,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)GH,"ESReInit3APTSpl() - no memory for GH");
    ESexit(0);
  }
  if(ReInitArray((void**)&BV2t,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)BV2t,"ESReInit3APTSpl() - no memory for BV2t");
    ESexit(0);
  }
  if(ReInitArray((void**)&BF,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)BF,"ESReInit3APTSpl() - no memory for BF");
    ESexit(0);
  }
  if(ReInitArray((void**)&BF2p,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)BF2p,"ESReInit3APTSpl() - no memory for BF2p");
    ESexit(0);
  }
  if(ReInitArray((void**)&AA22at,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AA22at,"ESReInit3APTSpl() - no memory for AA22at");
    ESexit(0);
  }
  if(ReInitArray((void**)&AA22ap,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AA22ap,"ESReInit3APTSpl() - no memory for AA22ap");
    ESexit(0);
  }
  if(ReInitArray((void**)&AA22tp,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AA22tp,"ESReInit3APTSpl() - no memory for AA22tp");
    ESexit(0);
  }
  if(ReInitArray((void**)&AA222atp,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AA222atp,"ESReInit3APTSpl()-no memory for AA222atp");
    ESexit(0);
  }
  
  if(ReInitArray((void**)&AF,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AF,"ESReInit3APTSpl() - no memory for AF");
    ESexit(0);
  }
  if(ReInitArray((void**)&BP,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)BP,"ESReInit3APTSpl() - no memory for BP");
    ESexit(0);
  }

  if(ReInitArray((void**)&AF2t,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AF2t,"ESReInit3APTSpl() - no memory for AF2t");
    ESexit(0);
  }
  if(ReInitArray((void**)&AF2p,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AF2p,"ESReInit3APTSpl() - no memory for AF2p");
    ESexit(0);
  }
  if(ReInitArray((void**)&AF22at,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AF22at,"ESReInit3APTSpl() - no memory for AF22at");
    ESexit(0);
  }
  if(ReInitArray((void**)&BP2p,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)BP2p,"ESReInit3APTSpl() - no memory for BP2p");
    ESexit(0);
  }
  if(ReInitArray((void**)&AF22tp,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AF22tp,"ESReInit3APTSpl() - no memory for AF22tp");
    ESexit(0);
  }
  if(ReInitArray((void**)&AF222atp,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AF222atp,"ESReInit3APTSpl()-no memory for AF222atp");
    ESexit(0);
  }
  
  if(ReInitArray((void**)&AP,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AP,"ESReInit3APTSpl() - no memory for AP");
    ESexit(0);
  }
  if(ReInitArray((void**)&AP2a,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AP2a,"ESReInit3APTSpl() - no memory for AP2a");
    ESexit(0);
  }
  if(ReInitArray((void**)&AP2t,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AP2t,"ESReInit3APTSpl() - no memory for AP2t");
    ESexit(0);
  }
  if(ReInitArray((void**)&AP2p,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AP2p,"ESReInit3APTSpl() - no memory for AP2p");
    ESexit(0);
  }
  if(ReInitArray((void**)&BA,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)BA,"ESReInit3APTSpl() - no memory for BA");
    ESexit(0);
  }
  if(ReInitArray((void**)&AP22ap,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AP22ap,"ESReInit3APTSpl() - no memory for AP22ap");
    ESexit(0);
  }
  if(ReInitArray((void**)&AP22tp,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)AP22tp,"ESReInit3APTSpl() - no memory for AP22tp");
    ESexit(0);
  }
  if(ReInitArray((void**)&BA2p,0,ESnAPT,sizeof(double)) < 0){
    FailureAlarm((char*)BA2p,"ESReInit3APTSpl()-no memory for BA2p");
    ESexit(0);
  }

  if(ReInitArray((void**)&dR0T,0,ESNt1,sizeof(double)) < 0){
    FailureAlarm((char*)dR0T,"ESReInit3APTSpl()-no memory for dR0T");
    ESexit(0);
  }
  if(ReInitArray((void**)&dZ0T,0,ESNt1,sizeof(double)) < 0){
    FailureAlarm((char*)dZ0T,"ESReInit3APTSpl()-no memory for dZ0T");
    ESexit(0);
  }
  if(ReInitArray((void**)&dbT,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)dbT,"ESReInit3APTSpl()-no memory for dbT");
    ESexit(0);
  }
  if(ReInitArray((void**)&dbT1a,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)dbT1a,"ESReInit3APTSpl()-no memory for dbT1a");
    ESexit(0);
  }
  if(ReInitArray((void**)&dbT2a,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)dbT2a,"ESReInit3APTSpl()-no memory for dbT2a");
    ESexit(0);
  }

  if(ReInitArray((void**)&bT2t,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)bT2t,"ESReInit3APTSpl()-no memory for bT2t");
    ESexit(0);
  }
  if(ReInitArray((void**)&EStgF,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)EStgF,"ESReInit3APTSpl() - no memory for EStgF");
    ESexit(0);
  }
  if(ReInitArray((void**)&EStgF1a,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)EStgF1a,"ESReInit3APTSpl() - no memory for EStgF1a");
    ESexit(0);
  }
  if(ReInitArray((void**)&EStgF3a,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)EStgF3a,"ESReInit3APTSpl() - no memory for EStgF3a");
    ESexit(0);
  }
  if(ReInitArray((void**)&EStgY,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)EStgY,"ESReInit3APTSpl() - no memory for EStgY");
    ESexit(0);
  }
  if(ReInitArray((void**)&EStgY1a,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)EStgY1a,"ESReInit3APTSpl() - no memory for EStgY1a");
    ESexit(0);
  }
  if(ReInitArray((void**)&EStgY1a2t,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)EStgY1a2t,"ESReInit3APTSpl() - no memory for EStgY1a2t");
    ESexit(0);
  }
  if(ReInitArray((void**)&EStgY3a,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)EStgY3a,"ESReInit3APTSpl() - no memory for EStgY3a");
    ESexit(0);
  }
  if(ReInitArray((void**)&gX0,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)gX0,"ESReInit3APTSpl() - no memory for gX0");
    ESexit(0);
  }
  if(ReInitArray((void**)&gX01a,0,ESnAT,sizeof(double)) < 0){
    FailureAlarm((char*)gX01a,"ESReInit3APTSpl() - no memory for gX01a");
    ESexit(0);
  }
  if(ReInitArray((void**)&ESgi,0,ESNa1,sizeof(double)) < 0){
    FailureAlarm((char*)ESgi,"ESReInit3APTSpl() - no memory for ESgi");
    ESexit(0);
  }
  if(ReInitArray((void**)&ESgi1a,0,ESNa1,sizeof(double)) < 0){
    FailureAlarm((char*)ESgi1a,"ESReInit3APTSpl() - no memory for ESgi1a");
    ESexit(0);
  }
  if(ReInitArray((void**)&ESgi2a,0,ESNa1,sizeof(double)) < 0){
    FailureAlarm((char*)ESgi2a,"ESReInit3APTSpl() - no memory for ESgi2a");
    ESexit(0);
  }
  if(ReInitArray((void**)&R0T2t,0,ESNt1,sizeof(double)) < 0){
    FailureAlarm((char*)R0T2t,"ESReInit3APTSpl()-no memory for R0T2t");
    ESexit(0);
  }
  if(ReInitArray((void**)&Z0T2t,0,ESNt1,sizeof(double)) < 0){
    FailureAlarm((char*)Z0T2t,"ESReInit3APTSpl()-no memory for Z0T2t");
    ESexit(0);
  }
  if(ReInitArray((void**)&ESgF0,0,ESNa1,sizeof(double)) < 0){
    FailureAlarm((char*)ESgF0,"ESReInit3APTSpl() - no memory for ESgF0");
    ESexit(0);
  }
  if(ReInitArray((void**)&ESgF01a,0,ESNa1,sizeof(double)) < 0){
    FailureAlarm((char*)ESgF01a,"ESReInit3APTSpl() - no memory for ESgF01a");
    ESexit(0);
  }
  if(ReInitArray((void**)&ESgY0,0,ESNa1,sizeof(double)) < 0){
    FailureAlarm((char*)ESgY0,"ESReInit3APTSpl() - no memory for ESgY0");
    ESexit(0);
  }
  if(ReInitArray((void**)&ESgY01a,0,ESNa1,sizeof(double)) < 0){
    FailureAlarm((char*)ESgY01a,"ESReInit3APTSpl() - no memory for ESgY01a");
    ESexit(0);
  }
  if(ReInitArray((void**)&ESgY02a,0,ESNa1,sizeof(double)) < 0){
    FailureAlarm((char*)ESgY02a,"ESReInit3APTSpl() - no memory for ESgY02a");
    ESexit(0);
  }
  if(ReInitArray((void**)&ESaY,0,ESNa1,sizeof(double)) < 0){
    FailureAlarm((char*)ESaY,"ESReInit3APTSpl() - no memory for ESaY");
    ESexit(0);
  }
  if(ReInitArray((void**)&ESaY1a,0,ESNa1,sizeof(double)) < 0){
    FailureAlarm((char*)ESaY1a,"ESReInit3APTSpl() - no memory for ESaY1a");
    ESexit(0);
  }
  if(ReInitArray((void**)&ESaY2a,0,ESNa1,sizeof(double)) < 0){
    FailureAlarm((char*)ESaY2a,"ESReInit3APTSpl() - no memory for ESaY2a");
    ESexit(0);
  }
  NaNpNt	=ESNt1 > ESNa1 ? ESNt1 : ESNa1;
  if(NaNpNt < ESNp1)
    NaNpNt	=ESNp1;
  if(ginf != NULL){
    free(ginf);
  }
  ginf	=(double*)malloc(6*NaNpNt*sizeof(double));
  if(ginf == NULL){
    printf("ESReInit3APTSpl() - no memory for ginf\n");
    ESexit(0);
  }
  g2inf	=ginf+NaNpNt;
  if(ESFaux != NULL){
    free(ESFaux);
  }
  ESFaux	=(double*)malloc(10*NaNpNt*sizeof(double));
  if(ESFaux == NULL){
    printf("ESReInit3APTSpl() - no memory for ESFaux\n");
    ESexit(0);
  }
  ESd2Faux	=ESFaux+NaNpNt;
  
  ESmemlev |=0x00100000;
  return(0);
}

int ESDeInit3APTSpl()
{
  if(ESFaux != NULL){
    free(ESFaux);
    ESFaux	=NULL;
  }
  if(ginf != NULL){
    free(ginf);
    ginf	=NULL;
  }
  free(ESaY2a);
  free(ESaY1a);
  free(ESaY);
  free(ESgY02a);
  free(ESgY01a);
  free(ESgY0);
  free(ESgF01a);
  free(ESgF0);
  free(Z0T2t);
  free(R0T2t);
  free(ESgi2a);
  free(ESgi1a);
  free(ESgi);
  free(gX01a);
  free(gX0);
  free(EStgY3a);
  free(EStgY1a2t);
  free(EStgY1a);
  free(EStgY);
  free(EStgF3a);
  free(EStgF1a);
  free(EStgF);
  free(bT2t);

  free(dbT2a);
  free(dbT1a);
  free(dbT);
  free(dZ0T);
  free(dR0T);
  free(BA2p);
  free(AP22tp);
  free(BA);
  free(AP2p);
  free(AP2t);
  free(AP2a);
  free(AP);
  
  free(AF222atp);
  free(AF22tp);
  free(AF22at);
  free(AF2p);
  free(AF2t);
  free(BP);
  free(AF);
  
  free(AA222atp);
  free(AA22tp);
  free(AA22ap);
  free(AA22at);
  free(BF2p);
  free(BF);
  free(BV2t);
  free(GH);
  free(BV);
  
  free(gX222atp);
  free(gX22tp);
  free(gX22ap);
  free(gX22at);
  free(gX2p);
  free(gX2t);
  free(gX2a);
  free(gX);
  
  free(DP222atp);
  free(DP22tp);
  free(DP22ap);
  free(DP22at);
  free(DP2p);
  free(DP2t);
  free(DP2a);
  free(DP);
  return(0);
}

int ESReInitgXcs()
{
  nMN	=ESNa1*nM*(2*nN+1);
  if(NnMN < nMN){
    if(ReInitArray((void**)&gXc,NnMN,nMN,sizeof(double)) < 0){
      FailureAlarm((char*)gXc,"ESReInitgXcs() - no memory for gXc");
      ESexit(0);
    }
    if(ReInitArray((void**)&gXc2a,NnMN,nMN,sizeof(double)) < 0){
      FailureAlarm((char*)gXc2a,"ESReInitgXcs() - no memory for gXc2a");
      ESexit(0);
    }
    if(ReInitArray((void**)&gXs,NnMN,nMN,sizeof(double)) < 0){
      FailureAlarm((char*)gXs,"ESReInitgXcs() - no memory for gXs");
      ESexit(0);
    }
    if(ReInitArray((void**)&gXs2a,NnMN,nMN,sizeof(double)) < 0){
      FailureAlarm((char*)gXs2a,"ESReInitgXcs() - no memory for gXs2a");
      ESexit(0);
    }
    NnMN=nMN;
  }
  if(NngYR < ngYR){
    int i,ii;
    i	=ESNa1*NngYR;
    ii	=ESNa1*ngYR;

    if(ReInitArray((void**)&gYRc,i,ii,sizeof(double)) < 0){
      FailureAlarm((char*)gYRc,"ESReInitgYRcs() - no memory for gYRc");
      ESexit(0);
    }
    if(ReInitArray((void**)&gYRc2a,i,ii,sizeof(double)) < 0){
      FailureAlarm((char*)gYRc2a,"ESReInitgYRcs() - no memory for gYRc2a");
      ESexit(0);
    }
    if(ReInitArray((void**)&gYRs,i,ii,sizeof(double)) < 0){
      FailureAlarm((char*)gYRs,"ESReInitgYRcs() - no memory for gYRs");
      ESexit(0);
    }
    if(ReInitArray((void**)&gYRs2a,i,ii,sizeof(double)) < 0){
      FailureAlarm((char*)gYRs2a,"ESReInitgYRcs() - no memory for gYRs2a");
      ESexit(0);
    }
    if(ReInitArray((void**)&NgYR,NngYR,ngYR,sizeof(int)) < 0){
      FailureAlarm((char*)NgYR,"ESReInitgYRcs() - no memory for NgYR");
      ESexit(0);
    }
    if(ReInitArray((void**)&MgYR,NngYR,ngYR,sizeof(int)) < 0){
      FailureAlarm((char*)MgYR,"ESReInitgYRcs() - no memory for MgYR");
      ESexit(0);
    }
    NngYR	=ngYR;
  }
  ESmemlev |=0x00400000;
  return(0);
}

int ESDeInitgXcs()
{
  if(NnMN){
    free(gXs2a);
    free(gXs);
    free(gXc2a);
    free(gXc);
    NnMN	=0;
  }
  if(NngYR){
    free(MgYR);
    free(NgYR);
    free(gYRs2a);
    free(gYRs);
    free(gYRc2a);
    free(gYRc);
    NngYR	=0;
  }
  return(0);
}

int splT(double *g,double *d2g)
{
  static int Nv=0;
  static double *v,*d2f,H,rH;
  int i,i1;
  double a1,a2;
  if(Nv < ESNt1){
    if(Nv){
      free(v);
    }
    Nv	=ESNt1;
    v	=(double*)malloc(2*Nv*sizeof(double));
    d2f	=v+Nv;

    d2f[0]	=1.;
    v[0]	=0.;
    i1	=0;
    for(i=1; i < ESNt; i++){
      v[i]	=-1./(4.+v[i1]);
      d2f[i]	=d2f[i1]*v[i];
      i1++;
    }
    i	=ESNt;
    d2f[i]=d2f[0];
    i1	=i-1;
    while(i > 0){
      d2f[i1]	+=v[i1]*d2f[i];
      i--;
      i1--;
    }
    rH	=1./(d2f[1]+4.*d2f[0]+d2f[ESNt-1]);
    H	=6.*crht*crht;
  }

  d2g[0]=(g[1]-2.*g[0]+g[ESNt-1])*crht*crht;
  v[0]	=0.;
  a2	=(g[1]-g[0])*H;
  i1	=0;
  for(i=1; i < ESNt; i++){
    a1	=a2;
    a2	=(g[i+1]-g[i])*H;
    v[i]	=-1./(4.+v[i1]);
    d2g[i]	=(d2g[i1]-a2+a1)*v[i];
    i1++;
  }
  i	=ESNt;
  d2g[i]=d2g[0];
  i1	=i-1;
  while(i > 0){
    d2g[i1]	+=v[i1]*d2g[i];
    i--;
    i1--;
  }
  i1	=ESNt-1;
  a1	=(2.*d2g[0]-d2g[1]-d2g[i1])*rH;
  for(i=0; i < ESNt1; i++){
    d2g[i]	+=a1*d2f[i]; 
  }
  return(0);
}

int f2splT(double *g,double *d2g,double *f,
	   double EZga0,double ga_1,double EZga2,double ga_3)
{
  static int Isize=0;
  static double Cga0=0.,Cga1=0.,Cga2=0.,Cga3=0.,EZga1,EZga3;
  static double w[6],rw[6],RHS[9],H,rH;
  static double t11,t13,t22,t23,t31,t32;
  static double *W,*V,*U0[3],*U1[3],*U2[3],*pV,*pW;
  int i,ii;

  if(Isize < ESNt1){
    if(Isize)
      free(glT);
    Isize=ESNt1;
    glT	=(double*)malloc((24*ESNt1+ESNt1)*sizeof(double));
    W		=glT	+ESNt1;
    V		=W	+6*ESNt1;
    U0[0]	=V	+9*ESNt1;
    U0[1]	=U0[0]	+ESNt1;
    U0[2]	=U0[1]	+ESNt1;
    U1[0]	=U0[2]	+ESNt1;
    U1[1]	=U1[0]	+ESNt1;
    U1[2]	=U1[1]	+ESNt1;
    U2[0]	=U1[2]	+ESNt1;
    U2[1]	=U2[0]	+ESNt1;
    U2[2]	=U2[1]	+ESNt1;
    H		=gf1T[1]*gf1T[1];
    rH		=1./H;
    H		*=EZcr6;


    ESmemlev |=0x00800000;
  }
  if(Cga0 != EZga0 || Cga1 != ga_1 || Cga2 != EZga2 || Cga3 !=ga_3){
    Cga0	=EZga0;
    Cga1	=ga_1;
    Cga2	=EZga2;
    Cga3	=ga_3;
    EZga1		=ga_1*rH;
    EZga3		=ga_3*rH;
    w[0]	=1.+EZga0+2.*EZga1;
    w[1]	=0.;
    w[2]	=2.;
    w[3]	=EZga2+2.*EZga3;
    w[4]	=4.*H;
    w[5]	=0.;
    t11		=EZga1;
    t13		=1.;
    t22		=EZga3;
    t23		=-H;
    t31		=1.;
    t32		=-H;
    rw[0]	=w[0];
    rw[1]	=w[1];
    rw[2]	=w[2];
    rw[3]	=w[3];
    rw[4]	=w[4];
    rw[5]	=w[5];

    U0[0][0]	=1.;
    U0[1][0]	=0.;
    U0[2][0]	=0.;
    U1[0][0]	=0.;
    U1[1][0]	=1.;
    U1[2][0]	=0.;
    U2[0][0]	=0.;
    U2[1][0]	=0.;
    U2[2][0]	=1.;

    pV		=V;
    pW		=W;
    ii	=0;
    for(i=1; i < ESNt; i++){
      EZinv3x3(rw,pW);
      pV[0]	=pW[0]*t11+pW[2]*t31;
      pV[1]	=pW[1]*t22+pW[2]*t32;
      pV[2]	=pW[0]*t13+pW[1]*t23;
      pV[3]	=pW[1]*t11+pW[4]*t31;
      pV[4]	=pW[3]*t22+pW[4]*t32;
      pV[5]	=pW[1]*t13+pW[3]*t23;
      pV[6]	=pW[2]*t11+pW[5]*t31;
      pV[7]	=pW[4]*t22+pW[5]*t32;
      pV[8]	=pW[2]*t13+pW[4]*t23;
      rw[0]	=w[0]-t11*pV[0]-t13*pV[6];
      rw[1]	=    -t11*pV[1]-t13*pV[7];
      rw[2]	=w[2]-t11*pV[2]-t13*pV[8];
      rw[3]	=w[3]-t22*pV[4]-t23*pV[7];
      rw[4]	=w[4]-t22*pV[5]-t23*pV[8];
      rw[5]	=    -t31*pV[2]-t32*pV[5];

      U0[0][i]	=pV[0]*U0[0][ii]+pV[1]*U0[1][ii]+pV[2]*U0[2][ii];
      U0[1][i]	=pV[3]*U0[0][ii]+pV[4]*U0[1][ii]+pV[5]*U0[2][ii];
      U0[2][i]	=pV[6]*U0[0][ii]+pV[7]*U0[1][ii]+pV[8]*U0[2][ii];

      U1[0][i]	=pV[0]*U1[0][ii]+pV[1]*U1[1][ii]+pV[2]*U1[2][ii];
      U1[1][i]	=pV[3]*U1[0][ii]+pV[4]*U1[1][ii]+pV[5]*U1[2][ii];
      U1[2][i]	=pV[6]*U1[0][ii]+pV[7]*U1[1][ii]+pV[8]*U1[2][ii];
      
      U2[0][i]	=pV[0]*U2[0][ii]+pV[1]*U2[1][ii]+pV[2]*U2[2][ii];
      U2[1][i]	=pV[3]*U2[0][ii]+pV[4]*U2[1][ii]+pV[5]*U2[2][ii];
      U2[2][i]	=pV[6]*U2[0][ii]+pV[7]*U2[1][ii]+pV[8]*U2[2][ii];
      pW	+=6;
      pV	+=9;
      ii++;
    }
    U0[0][i]	=1.;
    U0[1][i]	=0.;
    U0[2][i]	=0.;
    U1[0][i]	=0.;
    U1[1][i]	=1.;
    U1[2][i]	=0.;
    U2[0][i]	=0.;
    U2[1][i]	=0.;
    U2[2][i]	=1.;
    while(i > 1){
      pV -=9;
      U0[0][ii]	+=pV[0]*U0[0][i]+pV[1]*U0[1][i]+pV[2]*U0[2][i];
      U0[1][ii]	+=pV[3]*U0[0][i]+pV[4]*U0[1][i]+pV[5]*U0[2][i];
      U0[2][ii]	+=pV[6]*U0[0][i]+pV[7]*U0[1][i]+pV[8]*U0[2][i];

      U1[0][ii]	+=pV[0]*U1[0][i]+pV[1]*U1[1][i]+pV[2]*U1[2][i];
      U1[1][ii]	+=pV[3]*U1[0][i]+pV[4]*U1[1][i]+pV[5]*U1[2][i];
      U1[2][ii]	+=pV[6]*U1[0][i]+pV[7]*U1[1][i]+pV[8]*U1[2][i];

      U2[0][ii]	+=pV[0]*U2[0][i]+pV[1]*U2[1][i]+pV[2]*U2[2][i];
      U2[1][ii]	+=pV[3]*U2[0][i]+pV[4]*U2[1][i]+pV[5]*U2[2][i];
      U2[2][ii]	+=pV[6]*U2[0][i]+pV[7]*U2[1][i]+pV[8]*U2[2][i];
      i--;
      ii--;
    }
    i	=1;
    ii	=ESNt-1;
    RHS[0]	=V[0]*(U0[0][ii]+U0[0][i])+V[1]*(U0[1][ii]+U0[1][i])
      +V[2]*(U0[2][ii]+U0[2][i])-1.;
    RHS[1]	=V[3]*(U0[0][ii]+U0[0][i])+V[4]*(U0[1][ii]+U0[1][i])
      +V[5]*(U0[2][ii]+U0[2][i]);
    RHS[2]	=V[6]*(U0[0][ii]+U0[0][i])+V[7]*(U0[1][ii]+U0[1][i])
      +V[8]*(U0[2][ii]+U0[2][i]);

    RHS[3]	=V[0]*(U1[0][ii]+U1[0][i])+V[1]*(U1[1][ii]+U1[1][i])
      +V[2]*(U1[2][ii]+U1[2][i]);
    RHS[4]	=V[3]*(U1[0][ii]+U1[0][i])+V[4]*(U1[1][ii]+U1[1][i])
      +V[5]*(U1[2][ii]+U1[2][i])-1.;
    RHS[5]	=V[6]*(U1[0][ii]+U1[0][i])+V[7]*(U1[1][ii]+U1[1][i])
      +V[8]*(U1[2][ii]+U1[2][i]);
      
    RHS[6]	=V[0]*(U2[0][ii]+U2[0][i])+V[1]*(U2[1][ii]+U2[1][i])
      +V[2]*(U2[2][ii]+U2[2][i]);
    RHS[7]	=V[3]*(U2[0][ii]+U2[0][i])+V[4]*(U2[1][ii]+U2[1][i])
      +V[5]*(U2[2][ii]+U2[2][i]);
    RHS[8]	=V[6]*(U2[0][ii]+U2[0][i])+V[7]*(U2[1][ii]+U2[1][i])
      +V[8]*(U2[2][ii]+U2[2][i])-1.;

    {
      int indx[3];
      double ac[9],s;

      ac[0]	=RHS[0];
      ac[1]	=RHS[3];
      ac[2]	=RHS[6];
      ac[3]	=RHS[1];
      ac[4]	=RHS[4];
      ac[5]	=RHS[7];
      ac[6]	=RHS[2];
      ac[7]	=RHS[5];
      ac[8]	=RHS[8];
      LUdcmp(ac,3,indx,&s);
      RHS[0]	=1.;
      RHS[1]	=0.;
      RHS[2]	=0.;
      LUbksb(ac,3,indx,RHS);
      RHS[3]	=0.;
      RHS[4]	=1.;
      RHS[5]	=0.;
      LUbksb(ac,3,indx,RHS+3);
      RHS[6]	=0.;
      RHS[7]	=0.;
      RHS[8]	=1.;
      LUbksb(ac,3,indx,RHS+6);
      for(i=0; i < ESNt1; i++){
	ac[0]	=U0[0][i];
	ac[1]	=U0[1][i];
	ac[2]	=U0[2][i];
	ac[3]	=U1[0][i];
	ac[4]	=U1[1][i];
	ac[5]	=U1[2][i];
	ac[6]	=U2[0][i];
	ac[7]	=U2[1][i];
	ac[8]	=U2[2][i];
	U0[0][i]	=RHS[0]*ac[0]+RHS[1]*ac[3]+RHS[2]*ac[6];
	U0[1][i]	=RHS[0]*ac[1]+RHS[1]*ac[4]+RHS[2]*ac[7];
	U0[2][i]	=RHS[0]*ac[2]+RHS[1]*ac[5]+RHS[2]*ac[8];

	U1[0][i]	=RHS[3]*ac[0]+RHS[4]*ac[3]+RHS[5]*ac[6];
	U1[1][i]	=RHS[3]*ac[1]+RHS[4]*ac[4]+RHS[5]*ac[7];
	U1[2][i]	=RHS[3]*ac[2]+RHS[4]*ac[5]+RHS[5]*ac[8];

	U2[0][i]	=RHS[6]*ac[0]+RHS[7]*ac[3]+RHS[8]*ac[6];
	U2[1][i]	=RHS[6]*ac[1]+RHS[7]*ac[4]+RHS[8]*ac[7];
	U2[2][i]	=RHS[6]*ac[2]+RHS[7]*ac[5]+RHS[8]*ac[8];
      }
    }
  }

  pW	=W;
  pV	=V;
  glT[0]=pW[2]*f[0];
  d2g[0]=(f[1]-2.*f[0]+f[ESNt-1])*rH;
  RHS[1]	=W[1]*f[0]-d2g[0];
  RHS[2]	=W[2]*f[0]-glT[0];
  g[0]	=f[0];
  RHS[0]=W[0]*f[0]-g[0];

  ii	=0;
  for(i=1; i < ESNt; i++){
    d2g[i]	=pW[1]*f[i]+pV[3]*g[ii]+pV[4]*d2g[ii]+pV[5]*glT[ii];
    glT[i]	=pW[2]*f[i]+pV[6]*g[ii]+pV[7]*d2g[ii]+pV[8]*glT[ii];
    g[i]	=pW[0]*f[i]+pV[0]*g[ii]+pV[1]*d2g[ii]+pV[2]*glT[ii];
    pW	+=6;
    pV	+=9;
    ii++;
  }
  g[i]	=g[0];
  d2g[i]=d2g[0];
  glT[i]	=glT[0];
  while(i > 1){
    pV -=9;
    g[ii]	+=pV[0]*g[i]+pV[1]*d2g[i]+pV[2]*glT[i];
    d2g[ii]	+=pV[3]*g[i]+pV[4]*d2g[i]+pV[5]*glT[i];
    glT[ii]	+=pV[6]*g[i]+pV[7]*d2g[i]+pV[8]*glT[i];
    i--;
    ii--;
  }

  i	=1;
  ii	=ESNt-1;
  RHS[0]+=V[0]*(g[ii]+g[i])+V[1]*(d2g[ii]+d2g[i])+V[2]*(glT[ii]+glT[i]);
  RHS[1]+=V[3]*(g[ii]+g[i])+V[4]*(d2g[ii]+d2g[i])+V[5]*(glT[ii]+glT[i]);
  RHS[2]+=V[6]*(g[ii]+g[i])+V[7]*(d2g[ii]+d2g[i])+V[8]*(glT[ii]+glT[i]);
  for(i=0; i < ESNt1; i++){
    g[i]	-=RHS[0]*U0[0][i]+RHS[1]*U1[0][i]+RHS[2]*U2[0][i];
    d2g[i]	-=RHS[0]*U0[1][i]+RHS[1]*U1[1][i]+RHS[2]*U2[1][i];
  }
  return(0);
}

static int FlsplAPT=0;
int ESInitSplAPT()
{
  if(FlsplAPT){
    return(0);
  }
  FlsplAPT	=1;

  ESInitSplA();
  ESInitSplP();

  if(ESNt){
    cht	=EZcr6*gf1T[1];
    crht	=1./gf1T[1];
    chht	=EZcr2*gf1T[1];
    ch2ht	=chht*cht;
    crHt	=1./ESNt;
  }
  chp	=EZcr6*ESgt[1];
  crhp	=1./ESgt[1];
  chhp	=EZcr2*ESgt[1];
  ch2hp	=chhp*chp;
  crHp	=1./ESNp;
  return(0);
}

int DeInitf2splT()
{
  if(glT != NULL){
    free(glT);
    glT	=NULL;
  }
  FlsplAPT	=0;
  return(0);
}

int ESGetgFPlV()
{
  int i,j,k,ki,kj;

  double b,r0,EZz0;
  double R,Z,Ar,Art,Bt,At,Az,Azt,gF00;
  double drp,dzp;

  b	=ESsb[ESNa];
  ki	=ESNa;
  r0	=ESaR0[0]+rcT[ki];
  EZz0	=ESaZ0[0]+rsT[ki];
  for(k=1; k < ESFp1; k++){ 
    ki		+=ESNa1;
    rc[k]	=2.*rcT[ki];
    rs[k]	=2.*rsT[ki];
  }
  gF00	=0.;
  for(j=0; j < ESNp; j++){
    R	=r0;
    drp	=0.;
    kj	=0;
    for(k=1; k < ESFp1; k++){
      kj	+=j;
      if(kj >= ESNp)
	kj	-=ESNp;
      R		+=rc[k]*EScs1[kj]+rs[k]*ESsn1[kj];
      drp	+=k*(-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj]);
    }
    Z	=EZz0+b*ESsn1[j];
    ESSpl2ArBtAtAz(&Ar,&Art,&Bt,&At,&Az,&Azt,R,Z,0);
    dzp	=b*EScs1[j];
    gF00	+=Ar*drp+Az*dzp;
  }
  ESgFPlV	=-gF00*crHp;
  ESgFPlV1a	=2.*ESgFPlV;
  return(0);
}

int ES3DCoord2Spl()
{
  int k,ki;
  int i,in,n,n1;

  if(ESNt == 0){
    return(0);
  }

  splT(ESaR0,R0T2t);
  splT(ESaZ0,Z0T2t);

  in	=0;
  for(n=0; n < ESNt1; n++){
    ginf[n]	=ESsb1a[in];
    in	+=ESNa1;
  }
  splT(ginf,g2inf);
  in	=0;
  for(n=0; n < ESNt1; n++){
    bT2t[in]	=g2inf[n];
    in	+=ESNa1;
  }
  for(i=1; i < ESNa1; i++){
    in	=i;
    for(n=0; n < ESNt1; n++){
      ginf[n]	=ESsb[in];
      in	+=ESNa1;
    }
    splT(ginf,g2inf);
    in	=i;
    for(n=0; n < ESNt1; n++){
      bT2t[in]	=g2inf[n];
      in	+=ESNa1;
    }
  }

  for(k=0; k < ESFp1; k++){
    for(i=0; i < ESNa1; i++){
      ki	=ESNa1*k+i;
      for(n=0; n < ESNt1; n++){
	ginf[n]	=rcT[ki];
	ki	+=ESnAF;
      }
      ginf[ESNt]	=ginf[0];
      splT(ginf,g2inf);
      ki	=ESNa1*k+i;
      n1	=1;
      for(n=0; n < ESNt; n++){
	rcT1t[ki]=(ginf[n1]-ginf[n])*crht-(g2inf[n1]+2.*g2inf[n])*cht;
	ginf[n]	=rsT[ki];
	n1++;
	ki	+=ESnAF;
      }
      ginf[ESNt]	=ginf[0];
      rcT1t[ki]	=rcT1t[ki-ESnAF*ESNt];
      splT(ginf,g2inf);
      ki	=ESNa1*k+i;
      n1	=1;
      for(n=0; n < ESNt; n++){
	rsT1t[ki]=(ginf[n1]-ginf[n])*crht-(g2inf[n1]+2.*g2inf[n])*cht;
	n1++;
	ki	+=ESnAF;
      }
      rsT1t[ki]	=rsT1t[ki-ESnAF*ESNt];
    }
  }
  return(0);
}

int ESAextBext2Spl()
{
  int k,ki,kj;
  int i,i1,j,n,n1,in,in1,ji,jn,ji1,jni;
  double *Bf,*Bf2p;
  double *y0,dg0,dg1;
  double R,r0,Z,EZz0,Ar0,Ar,Az0,Az,s,ss;
  double EZdra,EZdza,drp,dzp,EZdrt,EZdzt;
  double x0a,x2ca,x2sa,y0a,b,ba,A,AA;
  double R0t,Z0t,bt;
  double Br,Br0,Brr,Brz,Bt,Bt0,Btr,Btz,Bz,Bz0,Bzr,Bzz;

  Bf	=g2inf+NaNpNt;
  Bf2p	=Bf+ESNp1;
 
  for(n=0; n < ESNLt; n++){
    jni	=ESnAP*n;
    in	=ESNa1*n;
    ji	=jni;
    i		=0;
    if(n){
      n1	=n-1;
      in1	=in-ESNa1;
      bt	=(ESsb1a[in]-ESsb1a[in1])*crht+(2.*bT2t[in]+bT2t[in1])*cht;
      R0t	=(ESaR0[n]-ESaR0[n1])*crht+(2.*R0T2t[n]+R0T2t[n1])*cht;
      Z0t	=(ESaZ0[n]-ESaZ0[n1])*crht+(2.*Z0T2t[n]+Z0T2t[n1])*cht;
    }
    else{
      n1	=n+1;
      in1	=in+ESNa1;
      bt	=(ESsb1a[in1]-ESsb1a[in])*crht-(bT2t[in1]+2.*bT2t[in])*cht;
      R0t	=(ESaR0[1]-ESaR0[0])*crht-(R0T2t[1]+2.*R0T2t[0])*cht;
      Z0t	=(ESaZ0[1]-ESaZ0[0])*crht-(Z0T2t[1]+2.*Z0T2t[0])*cht;
    }
    ba		=ESsb1a[in];
    R		=ESaR0[n];
    Z		=ESaZ0[n];
    ki		=ESnAF*n;
    x0a		=rcT2a[ki];
    y0a		=rsT2a[ki];
    k		=1;
    ki		+=ESNa1;
    rc[k]	=2.*rcT1a[ki];
    rs[k]	=2.*rsT1a[ki];
    k		=2;
    ki		+=ESNa1;
    rc[k]	=rcT2a[ki];
    rs[k]	=rsT2a[ki];
    s		=rc[1]*ba;
    ESSpl2BrBtBz(&Ar0,&Az0,&Br,&Bt,&Bz,R,Z,n);
    ESSpl2BrBtBzBrrBtrBzr(&Br0,&Bt0,&Bz0,&Brr,&Btr,&Bzr,&Brz,&Btz,&Bzz,R,Z,n);
    x2ca	=Bt0*(x0a+2.*rc[2]*ba-rs[1]*y0a);
    x2sa	=Bt0*(2.*rs[2]*ba+rc[1]*y0a);
    Br		=Br0*R-Bt0*R0t;
    Bz	    	=Bz0*R-Bt0*Z0t;
    for(j=0; j < ESNp; j++){
      drp	=-rc[1]*ESsn1[j]+rs[1]*EScs1[j];
      EZdra	=rc[1]*EScs1[j]+rs[1]*ESsn1[j];
      dzp	=ba*EScs1[j];
      EZdza	=ba*ESsn1[j];
      BV[ji]	=x2ca*ESsn1[j]-x2sa*EScs1[j]-(Btr*drp+Btz*dzp)*s;
      BA[ji]	=Br*dzp-Bz*drp;
      BF[ji]	=x2ca*EScs1[j]+x2sa*ESsn1[j]+(Btr*EZdra+Btz*EZdza)*s;
      BP[ji]	=Bz*EZdra-Br*EZdza;
      GH[ji]	=0.;
      AF[ji]	=-BP[ji];
      AP[ji]	=EZcr2*s*Bt0;
      ji++;
    }
    ji1		=ji-ESNp;
    GH[ji]	=GH[ji1];
    AP[ji]	=AP[ji1];
    AF[ji]	=AF[ji1];
    BV[ji]	=BV[ji1];
    BA[ji]	=BA[ji1];
    BF[ji]	=BF[ji1];
    BP[ji]	=BP[ji1];
    ji++;
    s		*=Bt0;
    gX0[in]	=-EZcr2*s;
    EStgF[in]	=EZcr2*s;
    EStgF1a[in]	=s;
    EStgY1a[in]	=Bz*x0a-Br*y0a+Bt0*rs1t[1]*ba
      +EZcr2*((Bz0+Bzr*R-Btr*Z0t)*(rc[1]*rc[1]+rs[1]*rs[1])
	    +((Btr*R0t+(Bzz-Brr)*R-Btz*Z0t)*ba-(Bt0*bt+Br0*ba))*rs[1]
	    +(Btz*R0t-Brz*R)*ba*ba
	    );
    in++;
    for(i=1; i < ESNa1; i++){
      ji	=jni+ESNp1*i;
      b		=ESsb[in];
      ba	=ESsb1a[in];
      ki	=ESnAF*n+i;
      rc[0]	=rcT[ki];
      rs[0]	=rsT[ki];
      rc1a[0]	=rcT1a[ki];
      rs1a[0]	=rsT1a[ki];
      rc1t[0]	=rcT1t[ki];
      rs1t[0]	=rsT1t[ki];
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	rc[k]	=2.*rcT[ki];
	rs[k]	=2.*rsT[ki];
	rc1a[k]	=2.*rcT1a[ki];
	rs1a[k]	=2.*rsT1a[ki];
	rc1t[k]	=2.*rcT1t[ki];
	rs1t[k]	=2.*rsT1t[ki];
      }
      if(n){
	in1	=in-ESNa1;
	bt	=(ESsb[in]-ESsb[in1])*crht+(2.*bT2t[in]+bT2t[in1])*cht;
      }
      else{
	in1	=in+ESNa1;
	bt	=(ESsb[in1]-ESsb[in])*crht-(bT2t[in1]+2.*bT2t[in])*cht;
      }
      r0	=ESaR0[n]+rc[0];
      EZz0	=ESaZ0[n]+rs[0];
      for(j=0; j < ESNp; j++){
	R	=r0;
	drp	=0.;
	EZdrt	=R0t+rc1t[0];
	EZdra	=rc1a[0];
	kj	=0;
	for(k=1; k < ESFp1; k++){
	  kj	+=j;
	  if(kj >= ESNp)
	    kj	-=ESNp;
	  R	+=rc[k]*EScs1[kj]+rs[k]*ESsn1[kj];
	  drp	+=k*(-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj]);
	  EZdrt	+=rc1t[k]*EScs1[kj]+rs1t[k]*ESsn1[kj];
	  EZdra	+=rc1a[k]*EScs1[kj]+rs1a[k]*ESsn1[kj];
	}
	Z	=EZz0+b*ESsn1[j];
	ESSpl2BrBtBz(&Ar,&Az,&Br,&Bt,&Bz,R,Z,n);
	EZdzt	=Z0t+rs1t[0]+bt*ESsn1[j];
	dzp	=b*EScs1[j];
	EZdza	=rs1a[0]+ba*ESsn1[j];
	AP[ji]	=(Ar-Ar0)*drp+(Az-Az0)*dzp;
	BA[ji]	=((Br*R-Bt*EZdrt)*dzp+(Bt*EZdzt-Bz*R)*drp)/ESsa[i];
	BP[ji]	=(Bz*R-Bt*EZdzt)*EZdra+(Bt*EZdrt-Br*R)*EZdza;
	BF[ji]	=Bt*(EZdra*dzp-drp*EZdza)/ESpa[i];
	ji++;
      }
      ji1	=ji-ESNp;
      BA[ji]	=BA[ji1];
      BF[ji]	=BF[ji1];
      BP[ji]	=BP[ji1];
      AP[ji]	=AP[ji1];
      ji	=jni+ESNp1*i;
      splP(BA+ji,BA2p+ji);
      splP(BF+ji,BF2p+ji);
      splP(BP+ji,BP2p+ji);
      splP(AP+ji,AP2p+ji);
      ji1	=ji;	
      AF[ji1]	=0.;
      EZdra	=0.;
      BV[ji1]	=0.;
      x0a	=0.;
      y0a	=0.;
      for(j=0; j < ESNp; j++){
	ji1++;
	AF[ji1]		=AF[ji]+chhp*(BA[ji]+BA[ji1]
				      -ch2hp*(BA2p[ji]+BA2p[ji1]));
	x0a		+=BA[ji]-ch2hp*BA2p[ji];
	BV[ji1]		=BV[ji]+chhp*(BF[ji]+BF[ji1]
				      -ch2hp*(BF2p[ji]+BF2p[ji1]));
	y0a		+=BP[ji]-ch2hp*BP2p[ji];
	EZdra		+=AP[ji]-ch2hp*AP2p[ji];
	ji++;
      }
      gX0[in]	=AF[ji1]*EZcr2gp/ESsa[i];
      x0a	*=crHp;
      EStgF1a[in]=BV[ji1]*EZcr2gp*ESsa[i];
      y0a 	*=crHp;
      EStgY1a[in]=y0a/ESsa[i];
      EStgF[in]	=-EZdra*crHp/ESpa[i];
      s		=BV[ji1]/ESgt[ESNp];
      ba	=AF[ji1]/ESgt[ESNp];
      EZdra	=ESsa[i]/EStgF1a[in];
      ji	=jni+ESNp1*i;
      for(j=0; j < ESNp1; j++){
	AF[ji]	-=ESgt[j]*ba;
	Bf[j]	=AF[ji];
	BV[ji]	-=ESgt[j]*s;
	BA[ji]	-=x0a;
	BF[ji]	-=s;
	BP[ji]	-=y0a;
	GH[ji]	=BV[ji]*EZdra;
	ginf[j]	=BV[ji];
	ji++;
      }
      splP(Bf,Bf2p);
      splP(ginf,g2inf);
      ba=0.;
      s	=0.;
      for(j=0; j < ESNp; j++){
	ba	+=Bf[j]-ch2hp*Bf2p[j];
	s	+=ginf[j]-ch2hp*g2inf[j];
      }
      ba	*=crHp;
      EZdra	=-ESpa[i]*gX0[in]/EStgF1a[in];
      s		*=crHp;
      ji	=jni+ESNp1*i;
      for(j=0; j < ESNp1; j++){
	BV[ji]	-=s;
	AF[ji]	+=BV[ji]*EZdra-ba;
	ji++;
      }
      in++;
    }
    n1	=n+ESNLt;
    while(n1 < ESNt1){
      ji	=ESnAP*n;
      ji1	=ESnAP*n1;
      in	=ESNa1*n;
      in1	=ESNa1*n1;
      for(i=0; i < ESNa1; i++){
	for(j=0; j < ESNp1; j++){
	  BV[ji1]=BV[ji];
	  AF[ji1]=AF[ji];
	  GH[ji1]=GH[ji];
	  BA[ji1]=BA[ji];
	  BP[ji1]=BP[ji];
	  BF[ji1]=BF[ji];
	  AP[ji1]=AP[ji];
	  ji++;
	  ji1++;
	}
	gX0[in1]	=gX0[in];
	EStgF[in1]	=EStgF[in];
	EStgF1a[in1]	=EStgF1a[in];
	EStgY1a[in1]	=EStgY1a[in];
	in++;
	in1++;
      }
      n1	+=ESNLt;
    }
  }
  for(i=0; i < ESNa1; i++){
    in	=i;
    for(n=0; n < ESNt1; n++){
      ginf[n]	=EStgF1a[in];
      in	+=ESNa1;
    }
    splT(ginf,g2inf);
    bt	=0.;
    for(n=0; n < ESNt; n++){
      bt	+=ginf[n]-ch2ht*g2inf[n];
    }
    bt	*=crHt;
    in	=i;
    for(n=0; n < ESNt1; n++){
      ginf[n]	=EStgY1a[in];
      in	+=ESNa1;
    }
    splT(ginf,g2inf);
    n1	=0;
    g2inf[n1]	=ginf[n1]-ch2ht*g2inf[n1];
    ginf[n1]	=0.;
    for(n=0; n < ESNt; n++){
      n1++;
      g2inf[n1]	=ginf[n1]-ch2ht*g2inf[n1];
      ginf[n1]	=ginf[n]+chht*(g2inf[n]+g2inf[n1]);
    }
    ba		=ginf[n1]*EZcr2gp;
    ESgF01a[i]	=bt;
    ESgY01a[i]	=ba;
    ESgi[i]	=ba/bt;
    if(i == 0){
      ESgimn	=ESgi[i];
      ESgimx	=ESgi[i];
    }
    else{
      if(ESgimn > ESgi[i])
	ESgimn	=ESgi[i];
      if(ESgimx < ESgi[i])
	ESgimx	=ESgi[i];
    }
    in	=i;
    for(n=0; n < ESNLt; n++){
      EStgY1a[in]	-=ba;
      gX0[in]	=(ESgFPlV-EStgF[in])/bt;
      ji	=n*ESnAP+i*ESNp1;
      s	=-(ginf[n]-gf1T[n]*ba)/EStgF1a[in];
      for(j=0; j < ESNp1; j++){
	GH[ji]	+=s;
	ji++;
      }
      in	+=ESNa1;
    }
  }

  for(n=0; n < ESNLt; n++){
    n1	=n+ESNLt;
    while(n1 < ESNt1){
      ji	=ESnAP*n;
      ji1	=ESnAP*n1;
      in	=ESNa1*n;
      in1	=ESNa1*n1;
      for(i=0; i < ESNa1; i++){
	for(j=0; j < ESNp1; j++){
	  GH[ji1]=GH[ji];
	  ji++;
	  ji1++;
	}
	gX0[in1]	=gX0[in];
	EStgY1a[in1]	=EStgY1a[in];
	in++;
	in1++;
      }
      n1	+=ESNLt;
    }
  }
  s	=2.*ESgi[ESNa]-ESgi[ESNa-1];
  if(ESgimn > s)
    ESgimn	=s;
  if(ESgimx < s)
    ESgimx	=s;
  for(i=0; i < ESNa1; i++){
    ginf[i]	=ESgY01a[i];
  }  
  s	=0.;
  splA(ginf,g2inf,&s,NULL);
  ESsplAE(EStgY,ginf,g2inf);
  for(i=1; i < ESNa1; i++){
    ji	=i*ESNp1;
    for(j=0; j < ESNp1; j++){
      in	=ji+j;
      jn	=in;
      for(n=0; n < ESNt1; n++){
	ginf[n]	=BV[jn];
	jn	+=ESnAP;
      }
      splT(ginf,g2inf);
      jn	=in;
      BP[jn]	+=ESpa[i]*((ginf[1]-ginf[0])*crht-(g2inf[1]+2.*g2inf[0])*cht);
      jn	+=ESnAP;
      n1	=0;
      for(n=1; n < ESNt; n++){
	BP[jn]	+=ESpa[i]*((ginf[n]-ginf[n1])*crht
			   +(2.*g2inf[n]+g2inf[n1])*cht);
	n1++;
	jn	+=ESnAP;
      }
      BP[jn]	=BP[in];
    }
  }  
  return(0);
}


#ifndef stg_NrmDspl
int ESNormDispl(int iA,int iT,int M,int N)
{
  int k,ki,kj,K;
  int i,j,n,jj,ji,in,Nerr,kR,iR;
  int m,m1,nw,nw1,MN,Np05,Nt05;
  double *gXc0,*gXcc,*gXcs,*gXs0,*gXsc,*gXss,gXcP,gXsP,cs,sn;
  double *gq,*gt1,*d2gt;
  double s,ss,gx,gc,gs,gi,err,gXmx,gYmx,A;
  static double EZga0=1e-5,EZga1=1e-9,EZga2=0.,EZga3=1e-10;
  double r0,EZz0,r1c,r1s,x2c,x2s,x0a,r1ca,r1sa,x2ca,x2sa,y0a,ba,EZdra,drp,EZdza;

  nM	=M;
  nN	=N;
  MN	=nM*(2*nN+1);

  kR	=0;
  m1	=1;
  for(m=0; m < nM; m++){
    nw1	=ESLt0;
    for(nw=0; nw < nN; nw++){
      s	=nw1/((double)m1);
      if(s >= ESgimn && s <= ESgimx){
	kR++;
      }
      s	=-s;
      if(s >= ESgimn && s <= ESgimx){
	kR++;
      }
      nw1	+=ESLt0;
    }
    m1++;
  }
  ngYR	=2*kR;

  ESReInitgXcs();
  Np05	=ESNp/2;
  Nt05	=ESNt/2;
  kR	=0;
  m1	=1;
  for(m=0; m < nM; m++){
    nw1	=ESLt0;
    for(nw=0; nw < nN; nw++){
      s	=nw1/((double)m1);
      if(s >= ESgimn && s <= ESgimx){
	MgYR[kR]	=m1;
	NgYR[kR]	=nw1;
	kR++;
      }
      s	=-s;
      if(s >= ESgimn && s <= ESgimx){
	MgYR[kR]	=m1;
	NgYR[kR]	=-nw1;
	kR++;
      }
      nw1 +=ESLt0;
    }
    nw1	=ESLt0;
    for(nw=0; nw < nN; nw++){
      s	=nw1/((double)m1);
      if(s >= ESgimn && s <= ESgimx){
	MgYR[kR]	=m1;
	NgYR[kR]	=nw1;
	kR++;
      }
      s	=-s;
      if(s >= ESgimn && s <= ESgimx){
	MgYR[kR]	=m1;
	NgYR[kR]	=-nw1;
	kR++;
      }
      nw1 +=ESLt0;
    }
    m1++;
  }

  gq	=g2inf+NaNpNt;
  gt1	=gq+NaNpNt;
  d2gt	=gt1+NaNpNt;

#ifdef H
  {
    double gh0,*pgY,ghmx,ghmn,gYmx,gYmn,x[4],y[2],dy,dx;
    ghmx=0.;
    ghmn=0.;
    gYmx=0.;
    gYmn=0.;
    in	=iA;
    s		=1./ESgF01a[iA];
    for(n=0; n < ESNt1; n++){
      ji	=ESnAP*n+ESNp1*iA;
      gh0	=GH[ji];
      if(ghmx < gh0)
	ghmx=gh0;
      if(ghmn > gh0)
	ghmn=gh0;
      pgY	=AF+ji;
#ifdef H
      f2splP(pgY,g2inf,pgY,0.,1e-6,0.,1e-7);
#endif
      for(j=0; j < ESNp1; j++){
	ss	=pgY[j]*s;
	if(gYmx < ss)
	  gYmx=ss;
	if(gYmn > ss)
	  gYmn=ss;
	ji++;
      }
      in	+=ESNa1;
    }
    ghmx	+=EZc2gp;
    x[0]	=ghmn*EZcr2gp;
    x[1]	=ghmx*EZcr2gp;
    dy		=0.1*(gYmx-gYmn);

    gi	=dy;
    printf("???? gi=%10.3e\n",gi);

    y[0]	=gYmn;
    y[1]	=dy*ESNt+gYmx;
   
    Scale2D(1,2,x,y,2,6,2);
    in	=iA;
    for(n=0; n < ESNt1; n++){
      ji	=ESnAP*n+ESNp1*iA;
      pgY	=AF+ji;
      EZz0	=dy*n;
      for(j=0; j < ESNp1; j++){
	ESgh[j]	=(ESgt[j]+GH[ji])*EZcr2gp;
	ginf[j]	=EZz0+pgY[j]*s;
	ji++;
      }
      Plot2d(1,ESgh,ginf,ESNp1,6,0,0,0);
      in	+=ESNa1;
    }

    gh0		=ESgi[iA]/EZz0;
    y[0]	=0.;
    y[1]	=0.;
    Plot2d(1,x,y,2,6,0,14,0);
    y[0]	=EZz0;
    y[1]	=EZz0;
    Plot2d(1,x,y,2,6,0,14,0);

    i	=x[0]/0.2;
    dx	=0.2*i;
    y[0]	=gYmn;
    x[2]	=0.2*i+gh0*y[0];
    y[1]	=EZz0+gYmx;
    x[3]	=0.2*i+gh0*y[1];

    while((x[0] <= x[2] && x[2] <= x[1]) || (x[0] <= x[3] && x[3] <= x[1])){
      if(x[2] < x[0]){
	x[2]	=x[0];
	y[0]	=(x[2]-0.2*i)/gh0;
      } 
      if(x[2] > x[1]){
	x[2]	=x[1];
	y[0]	=(x[2]-0.2*i)/gh0;
      } 
      if(x[3] < x[0]){
	x[3]	=x[0];
	y[1]	=(x[3]-0.2*i)/gh0;
      } 
      if(x[3] > x[1]){
	x[3]	=x[1];
	y[1]	=(x[3]-0.2*i)/gh0;
      } 
      Plot2d(1,x+2,y,2,6,0,0,0);
      i++;
      y[0]	=gYmn;
      x[2]	=0.2*i+gh0*y[0];
      y[1]	=EZz0+gYmx;
      x[3]	=0.2*i+gh0*y[1];
    }
  }
#endif

  in	=0;
  for(n=0; n < ESNLt; n++){
    for(i=0; i < ESNa1; i++){
      ji	=ESnAP*n+ESNp1*i;
      s		=1./ESgF01a[i];
      ss	=ESsa[i]*s;
      err	=0.;
      splP(AF+ji,AF2p+ji);
      for(j=0; j < ESNp1; j++){
	ginf[j]	=-GH[ji];
	gq[j]	=ESgt[j]+GH[ji];
	ji++;
      }
      EZspl_p(ginf,g2inf,gq,ESNp);
      k	=0;
      for(j=0; j < ESNp; j++){
	gs	=ESgt[j];
	while(gs > gq[ESNp])
	  gs -=EZc2gp;
	while(gs < gq[0])
	  gs +=EZc2gp;
	while(gq[k+1] < gs){
	  k++;
	}
	while(gq[k] > gs){
	  k--;
	}
	splr1(gt1+j,NULL,&gs,ginf,g2inf,gq,&ESNp,&k);
	gt1[j]	+=ESgt[j];
	if(0 && i == iA && n == ESNt-1){
	  printf("??? gt1[%2d]=%10.3e\n",j,gt1[j]);
	}
      }
      gt1[j]	=ESgt[j]+ginf[j];
      if(0 && i == iA && n == ESNt-1){
	printf("??? gt1[%2d]=%10.3e\n",j,gt1[j]);
      }
      ji	=ESnAP*n+ESNp1*i;
      k	=0;
      gc	=0.;
      s		=1./ESgF01a[i];
      for(j=0; j < ESNp; j++){
	gs	=gt1[j];
	while(gs > EZc2gp)
	  gs -=EZc2gp;
	while(gs < 0.)
	  gs +=EZc2gp;
	while(ESgt[k+1] < gs){
	  k++;
	}
	while(ESgt[k] > gs){
	  k--;
	}
	if(0 && i == iA && n == ESNt-1){
	  printf("??? gt[%2d]=%10.3e %10.3e %10.3e \n",k,ESgt[k],gs,ESgt[k+1]);
	}
	splr1(ginf+j,NULL,&gs,AF+ji,AF2p+ji,ESgt,&ESNp,&k);
	ginf[j]	*=s;
	gc	+=ginf[j];
	d2gt[j]	=ginf[j];
      }
      ginf[ESNp]=ginf[0];
      d2gt[ESNp]=ginf[ESNp];
#ifdef H
      f2splP(ginf,g2inf,ginf,0.,1e-6,0.,1e-7);
      splP(ginf,g2inf);
#endif
      gc	*=crHp;
      EZz0	=gc;
      k	=0;	
      ji	=ESnAP*n+ESNp1*i;
      for(j=0; j < ESNp1; j++){
	ginf[j]	-=gc;

	/****/
	AF[ji+j]=ginf[j];
	if(err < fabs(ginf[j])){
	  err	=fabs(ginf[j]);
	}
      }
      m1	=1;
      ji	=ESnAP*n+ESNp1*i;
      for(m=0; m < nM; m++){
	gc	=0.;
	gs	=0.;
	k	=0;	
	for(j=0; j < ESNp; j++){
	  gc	+=ginf[j]*EScs1[k];
	  gs	+=ginf[j]*ESsn1[k];
	  k	+=m1;
	  if(k >= ESNp){
	    k -=ESNp;
	  }
	}
	gX[ji]	=m1 < Np05 ? gc*crHp : EZcr2*gc*crHp;
	ji++;
	gX[ji]	=m1 < Np05 ? gs*crHp : EZcr2*gs*crHp;
	ji++;
	if(0 && i == iA && n == ESNt-1){
	  printf("??? i=%2d m=%2d n=%3d gc=%10.3e gs=%10.3e\n",i,m1,n
		 ,gc*crHp,gs*crHp);
	}
	m1++;
      }

      gs	=0.;
      for(j=0; j < ESNp; j++){
	ji	=ESnAP*n+ESNp1*i;
	gc	=0.;
	m1	=1;
	k	=0;
	for(m=0; m < nM; m++){
	  k	+=j;
	  if(k >= ESNp){
	    k -=ESNp;
	  }
	  gc	+=gX[ji]*EScs1[k];
	  ji++;
	  gc	+=gX[ji]*ESsn1[k];
	  ji++;
	  m1++;
	}
	gt1[j]	=2.*gc;
	gc	=fabs(AF[ESnAP*n+ESNp1*i+j]-2.*gc);
	if(gs < gc){
	  gs	=gc;
	}
      }
      gt1[j]	=gt1[0];

      if(i == iA){
	printf("??? Amx[%3d][%2d]=%10.3e err=%10.3e\n",n,i,err,gs/err);
      }
#ifdef XWIN
      if(i == iA && n == (iT%ESNLt)){
	printf("??? %3d gi=%10.3e\n",n,gi);
	s	=gi*n;
	for(j=0; j < ESNp1; j++){
	  ESgh[j]	=ESgt[j]*EZcr2gp;
	  d2gt[j]	=d2gt[j]*5.+s;
	  gt1[j]	=(gt1[j]+EZz0)*5.+s;
	  ginf[j]		=(ginf[j]+EZz0)*5.+s;
	}
	Plot2d(1,ESgh,d2gt,ESNp1,6,0,4,0);
	Plot2d(1,ESgh,ginf,ESNp1,6,0,5,0);
	Plot2d(1,ESgh,gt1,ESNp1,6,0,14,0);
      }
#endif
      in++;
    }
    nw	=n+ESNLt;
    while(nw < ESNt1){
      for(i=0; i < ESNa1; i++){
	ji	=ESnAP*n+ESNp1*i;
	m1	=ESnAP*nw+ESNp1*i;
	for(m=0; m < nM; m++){
	  gX[m1]	=gX[ji];
	  ji++;
	  m1++;
	  gX[m1]	=gX[ji];
	  ji++;
	  m1++;
	}
      }
      nw	+=ESNLt;
    }
  }
  gXmx	=0.;
  gYmx	=0.;
  ss	=EZcr2*crHt;
  for(i=0; i < ESNa1; i++){
    kR	=0;
    iR	=i;
    if(kR < ngYR){
      gYRc[iR]	=0.;
      gYRs[iR]	=0.;
    }
    gi	=ESgi[i];
    m1	=1;
    for(m=0; m < nM; m++){
      k		=MN*i+(2*nN+1)*m;
      gXc0	=gXc+k;
      gXcc	=gXc0+1;
      gXcs	=gXcc+nN;
      gXs0	=gXs+k;
      gXsc	=gXs0+1;
      gXss	=gXsc+nN;
      ji	=ESNp1*i+2*m;
      gc	=0.;
      for(n=0; n < ESNt; n++){
	ginf[n]	=gX[ji];
	gc	+=ginf[n];
	ji	+=ESnAP;
	if(0 && i == iA && m1 == 1){
	  printf("*?? i=%2d m=%2d n=%3d gc=%10.3e\n",i,m1,n,ginf[n]);
	}
      }
      gXc0[0]	=gc*crHt/gi;
      nw1	=ESLt0;
      for(nw=0; nw < nN; nw++){
	gc	=0.;
	gs	=0.;
	k	=0;
	for(n=0; n < ESNt; n++){
	  gc	+=ginf[n]*cs1T[k];
	  gs	+=ginf[n]*sn1T[k];
	  k	+=nw1;
	  if(k >=ESNt)
	    k -=ESNt;
	}
	if(nw1 == Nt05){
	  gc	*=ss;
	  gs	*=ss;
	}
	else{
	  gc	*=crHt;
	  gs	*=crHt;
	}
	s	=fabs(gc);
	if(gYmx < s)
	  gYmx	=s;
	s	=fabs(gs);
	if(gYmx < s)
	  gYmx	=s;
	s	=nw1/((double)m1);
	if(s < ESgimn || s > ESgimx){
	  gXcc[nw]	=gc/(gi-s);
	  gXss[nw]	=gc/(gi-s);
	  gXsc[nw]	=-gs/(gi-s);
	  gXcs[nw]	=gs/(gi-s);
	}
	else{
	  gXcc[nw]	=0.;
	  gXss[nw]	=0.;
	  gXsc[nw]	=0.;
	  gXcs[nw]	=0.;
	  gYRc[iR]	=gc;
	  gYRs[iR]	=-gs;
	  kR++;
	  iR	+=ESNa1;
	  if(kR < ngYR){
	    gYRc[iR]	=0.;
	    gYRs[iR]	=0.;
	  }
	}
	s	=-s;
	if(s < ESgimn || s > ESgimx){
	  gXcc[nw]	+=gc/(gi-s);
	  gXss[nw]	-=gc/(gi-s);
	  gXsc[nw]	+=gs/(gi-s);
	  gXcs[nw]	+=gs/(gi-s);
	}
	else{
	  gYRc[iR]	=gc;
	  gYRs[iR]	=gs;
	  kR++;
	  iR	+=ESNa1;
	  if(kR < ngYR){
	    gYRc[iR]	=0.;
	    gYRs[iR]	=0.;
	  }
	}
	nw1 +=ESLt0;
      }

      ji	=ESNp1*i+2*m+1;
      gc	=0.;
      for(n=0; n < ESNt; n++){
	ginf[n]	=gX[ji];
	gc	+=ginf[n];
	ji	+=ESnAP;
      }

      gXs0[0]	=gc*crHt/gi;
      nw1	=ESLt0;
      for(nw=0; nw < nN; nw++){
	gc	=0.;
	gs	=0.;
	k	=0;
	for(n=0; n < ESNt; n++){
	  gc	+=ginf[n]*cs1T[k];
	  gs	+=ginf[n]*sn1T[k];
	  k	+=nw1;
	  if(k >=ESNt)
	    k -=ESNt;
	}
	if(nw1 == Nt05){
	  gc	*=ss;
	  gs	*=ss;
	}
	else{
	  gc	*=crHt;
	  gs	*=crHt;
	}
	s	=fabs(gc);
	if(gYmx < s)
	  gYmx	=s;
	s	=fabs(gs);
	if(gYmx < s)
	  gYmx	=s;
	s	=nw1/((double)m1);
	if(s < ESgimn || s > ESgimx){
	  gXcc[nw]	+=gs/(gi-s);
	  gXss[nw]	+=gs/(gi-s);
	  gXcs[nw]	-=gc/(gi-s);
	  gXsc[nw]	+=gc/(gi-s);
	}
	else{
	  gYRs[iR]	=gc;
	  gYRc[iR]	=gs;
	  kR++;
	  iR	+=ESNa1;
	  if(kR < ngYR){
	    gYRc[iR]	=0.;
	    gYRs[iR]	=0.;
	  }
	}
	s	=-s;
	if(s < ESgimn || s > ESgimx){
	  gXcc[nw]	-=gs/(gi-s);
	  gXss[nw]	+=gs/(gi-s);
	  gXcs[nw]	+=gc/(gi-s);
	  gXsc[nw]	+=gc/(gi-s);
	}
	else{
	  gYRs[iR]	=gc;
	  gYRc[iR]	=-gs;
	  kR++;
	  iR	+=ESNa1;
	  if(kR < ngYR){
	    gYRc[iR]	=0.;
	    gYRs[iR]	=0.;
	  }
	}
	nw1 +=ESLt0;
      }
      m1++;
    }
  }
  err	=0.;
  Nerr	=0;
  {
    int jtp,jbt;
    double gXt1p,gXb1p,db,dy0,dx0,dr1c,dr1s,dx2c,dx2s;

    jtp	=ESNp/4;
    jbt	=ESNp-jtp;
    for(n=0; n < ESNLt; n++){
      K	=ESnAF*n;
      for(i=0; i < ESNa1; i++){
	ji	=ESnAP*n+ESNp1*i;
	for(j=0; j < ESNp1; j++){
	  gX[ji]	=0.;
	  ji++;
	}
	gXt1p	=0.;	
	gXb1p	=0.;	
	m1	=1;
	for(m=0; m < nM; m++){
	  k	=MN*i+(2*nN+1)*m;
	  gXc0	=gXc+k;
	  gXcc	=gXc0+1;
	  gXcs	=gXcc+nN;
	  gXs0	=gXs+k;
	  gXsc	=gXs0+1;
	  gXss	=gXsc+nN;
	  gXcP	=gXc0[0];
	  gXsP	=gXs0[0];
	  nw1	=ESLt0*n;
	  if(nw1 >= ESNt){
	    nw1 -=ESNt;
	  }
	  k	=nw1;
	  for(nw=0; nw < nN; nw++){
	    gXcP	+=gXcc[nw]*cs1T[k]+gXcs[nw]*sn1T[k];
	    gXsP	+=gXsc[nw]*cs1T[k]+gXss[nw]*sn1T[k];
	    k	+=nw1;
	    if(k >= ESNt){
	      k -=ESNt;
	    }
	  }
	  gXcP	*=2.;
	  gXsP	*=2.;
	  ji	=ESnAP*n+ESNp1*i;
	  for(j=0; j < ESNp1; j++){
	    s	=m1*(ESgt[j]+GH[ji]);
	    gc	=cos(s);
	    gs	=sin(s);
	    gX[ji]	+=gXcP*gc+gXsP*gs;
	    if(j == jtp){
	      gXt1p	+=m1*(gXsP*gc-gXcP*gs);
	    }
	    if(j == jbt){
	      gXb1p	+=m1*(gXsP*gc-gXcP*gs);
	    }
	    ji++;
	  }
	  if(i == 0){
	    break;
	  }
	  m1++;
	}
	ji	=ESnAP*n+ESNp1*i;
	in	=ESNa1*n+i;
	ba	=ESsb1a[in];
	if(i == 0){
	  ki	=K;
	  rc1a[0]	=0.;
	  rs1a[0]	=0.;
	  ki	+=ESNa1;
	  rc[1]	=2.*rcT1a[ki];
	  rs[1]	=2.*rsT1a[ki];
	  dy0	=EZcr2*ba*(gX[ji+jtp]-gX[ji+jbt]);
	  db	=EZcr2*ba*(gX[ji+jtp]+gX[ji+jbt]);
	  EZz0	=dy0;
	  dZ0T[n]	=EZz0;
	  drsT[K]	=0.;
	  dbT[in]	=0.;
	  gc		=cos(GH[ji]);
	  gs		=sin(GH[ji]);
	  r0		=(gXcP*gc+gXsP*gs)*rc[1]+(gXsP*gc-gXcP*gs)*rs[1];
	  dR0T[n]	=r0;
	  drcT[K]	=0.;
	  ki		=K;
	  for(k=1; k < ESFp1; k++){
	    ki		+=ESNa1;
	    drcT[ki]	=0.;
	    drsT[ki]	=0.;
	  }
	  for(j=0; j < ESNp1; j++){
	    if(n == iT){
	      ECr[j]	=ESaR0[n]+dR0T[n];
	      ECz[j]	=ESaZ0[n]+dZ0T[n];
	    }
	    if(gXmx < fabs(gX[ji])){
	      gXmx	=fabs(gX[ji]);
	    }
	    gX[ji]=0.;
	    ji++;
	  }
	}
	else{
	  A	=2./ESsa[i];
	  ki	=K+i;
	  rc1a[0]	=rcT1a[ki];
	  rs1a[0]	=rsT1a[ki];
	  ki	+=ESNa1;
	  rc[1]	=2.*rcT[ki]/ESsa[i];
	  rs[1]	=2.*rsT[ki]/ESsa[i];

	  s	=ESsa[i]/ESgF01a[i];
	  gXt1p	*=1.+BF[ji+jtp]*s;
	  gXb1p	*=1.+BF[ji+jbt]*s;
	  gc	=ESsa[i]/ESsb[in];
	  ki	=K+i+ESNa1;
	  rc1a[1]=2.*rcT1a[ki];
	  rs1a[1]=2.*rsT1a[ki];
	  for(k=2; k < ESFp1; k++){
	    ki	+=ESNa1;
	    rc[k]	=k*A*rcT[ki];
	    rs[k]	=k*A*rsT[ki];
	    rc1a[k]	=2.*rcT1a[ki];
	    rs1a[k]	=2.*rsT1a[ki];
	  }
	  dy0		=EZcr2*((rs1a[0]+ba)*gX[ji+jtp]+(rs1a[0]-ba)*gX[ji+jbt]);
	  drsT[K+i]	=dy0-EZz0;
	  db		=EZcr2*((rs1a[0]+ba)*gX[ji+jtp]-(rs1a[0]-ba)*gX[ji+jbt]);
	  dbT[in]	=ESsa[i]*db;
	  jj	=ESNp1*i;
	  for(j=0; j < ESNp; j++){
	    EZdza	=rs1a[0]+ba*ESsn1[j];
	    EZdra	=rc1a[0]+rc1a[1]*EScs1[j]+rs1a[1]*ESsn1[j];
	    drp	=-rc[1]*ESsn1[j];
	    kj	=j;
	    for(k=2; k < ESFp1; k++){
	      kj	+=j;
	      if(kj >= ESNp)
		kj	-=ESNp;
	      EZdra	+=rc1a[k]*EScs1[kj]+rs1a[k]*ESsn1[kj];
	      drp	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
	    }
	    drp		*=gc; /* !!!! */
	    gs		=dy0+db*ESsn1[j]-EZdza*gX[ji];
	    if(gXmx < fabs(gX[ji]))
	      gXmx	=fabs(gX[ji]);
	    if(n == iT){
	      ECr[jj]	=ESsr[ji]+EZdra*gX[ji];
	      ECz[jj]	=ESsz[ji]+EZdza*gX[ji];
	      jj++;
	    }
	    gX[ji]	=EZdra*gX[ji]+rs[1]*gs*gc-r0;
	    if(j != jtp){
	      if(j != jbt){
		gX[ji]	+=drp*gs/EScs1[j];
	      }
	      else{
		gX[ji]	-=drp*EZdza*gXb1p;
	      }
	    }
	    else{
	      gX[ji]	+=drp*EZdza*gXt1p;
	    }
	    ji++;
	  }
	  gX[ji]	=gX[ji-ESNp];
	  if(n == iT){
	    ECr[jj]	=ECr[jj-ESNp];
	    ECz[jj]	=ECz[jj-ESNp];
	  }
	  ji++;
	  ESP2F(rc,rs,gX+ESnAP*n+ESNp1*i,ESFp);
	  ki	=K+i;
	  drcT[ki]	=rc[0];
	  for(k=1; k < ESFp1; k++){
	    ki	+=ESNa1;
	    drcT[ki]	=rc[k];
	    drsT[ki]	=rs[k];
	  }
	}
      }
      in	=ESNa1*n;
      gx	=fabs(1.-EStgF[in]/ESgFPlV); 
      if(err < gx){
	err	=gx;
	Nerr	=n;
      }
      for(i=1; i < ESNa1; i++){
	in++;
	gx	=fabs(1.-EStgF[in]/ESgFPlV); 
	if(err < gx){
	  err	=gx;
	  Nerr	=n;
	}
	gx		=ESsa[i]*gX0[in];
	ki		=K+i;
	drcT[ki]	+=rcT1a[ki]*gx;
	drsT[ki]	+=rsT1a[ki]*gx;
	dbT[in]		+=ESsb1a[in]*gx;
	for(k=1; k < ESFp1; k++){
	  ki	+=ESNa1;
	  drcT[ki]	+=rcT1a[ki]*gx;
	  drsT[ki]	+=rsT1a[ki]*gx;
	}
      }
#ifdef H
      if(0 && n == 0){
	double x[2],y[2];
	ki	=K;
	y[0]	=drcT[ki];
	y[1]	=y[0];
	x[0]	=0.;
	x[1]	=1.;
	for(k=0; k < ESFp1; k++){
	  printf("k=%2d\n",k);
	  for(i=0; i < ESNa1; i++){
	    if(y[0] > drcT[ki])
	      y[0]=drcT[ki];
	    if(y[1] < drcT[ki])
	      y[1]=drcT[ki];
	    if(y[0] > drsT[ki])
	      y[0]=drsT[ki];
	    if(y[1] < drsT[ki])
	      y[1]=drsT[ki];
	    printf("dr[%2d]=%10.3e %10.3e*\n",i,drcT[ki],drsT[ki]);
	    ki++;
	  }
	}
	Scale2D(1,2,x,y,2,6,2);
	ki	=K;
	for(k=0; k < ESFp1; k++){
	  Plot2d(1,a,drcT+ki,ESNa1,6,0,0,0);
	  Plot2d(1,a,drsT+ki,ESNa1,6,0,0,0);
	  ki	+=ESNa1;
	}
	CbFlush();
      }
#endif
      ki	=K;
      in	=ESNa1*n;
      A	=0.;
      EZf2spl(dbT+in,dbT2a+in,NULL,NULL,dbT+in,NULL,ESsa,&ESNa
	    ,&EZga0,&EZga1,&EZga2,&EZga3);
      EZf2spl(drcT+ki,drcT2a+ki,&A,NULL,drcT+ki,NULL,ESsa,&ESNa
	    ,&EZga0,&EZga1,&EZga2,&EZga3);
      EZf2spl(drsT+ki,drsT2a+ki,&A,NULL,drsT+ki,NULL,ESsa,&ESNa
	    ,&EZga0,&EZga1,&EZga2,&EZga3);
      splA1a(dbT+in,dbT1a+in,dbT2a+in);
      splA1a(drcT+ki,drcT1a+ki,drcT2a+ki);
      splA1a(drsT+ki,drsT1a+ki,drsT2a+ki);
      ginf[0]	=0.;
      ki	+=ESNa1;
      in	+=ESNa1;
      EZf2spl(drcT+ki,drcT2a+ki,NULL,NULL,drcT+ki,NULL,ESsa,&ESNa,
	    &EZga0,&EZga1,&EZga2,&EZga3);
      EZf2spl(drsT+ki,drsT2a+ki,NULL,NULL,drsT+ki,NULL,ESsa,&ESNa,
	    &EZga0,&EZga1,&EZga2,&EZga3);
      splA1a(drcT+ki,drcT1a+ki,drcT2a+ki);
      splA1a(drsT+ki,drsT1a+ki,drsT2a+ki);
      ki	+=ESNa1;
      for(k=2; k < ESFp1; k++){
	EZf2spl(drcT+ki,drcT2a+ki,&A,NULL,drcT+ki,NULL,ESsa,&ESNa,
	      &EZga0,&EZga1,&EZga2,&EZga3);
	EZf2spl(drsT+ki,drsT2a+ki,&A,NULL,drsT+ki,NULL,ESsa,&ESNa,
	      &EZga0,&EZga1,&EZga2,&EZga3);
	splA1a(drcT+ki,drcT1a+ki,drcT2a+ki);
	splA1a(drsT+ki,drsT1a+ki,drsT2a+ki);
	ki	+=ESNa1;
      }

#ifdef H
      if(n == 0){
	double x[50],y[50];
	ki	=K;
	y[0]	=drcT[ki];
	y[1]	=y[0];
	x[0]	=0.;
	x[1]	=1.;
	for(k=0; k < ESFp1; k++){
	  for(i=0; i < ESNa1; i++){
	    if(y[0] > drcT[ki])
	      y[0]=drcT[ki];
	    if(y[1] < drcT[ki])
	      y[1]=drcT[ki];
	    if(y[0] > drsT[ki])
	      y[0]=drsT[ki];
	    if(y[1] < drsT[ki])
	      y[1]=drsT[ki];
	    ki++;
	  }
	}
	Scale2D(1,2,x,y,2,6,2);
	ki	=K;
	for(i=0; i < ESNa1; i++){
	  x[i]	=drcT[ki];
	  y[i]	=drsT[ki];
	  ki++;
	}
	Plot2d(1,a,x,ESNa1,6,0,4,0);
	Plot2d(1,a,y,ESNa1,6,0,4,0);
	for(k=1; k < ESFp1; k++){
	  for(i=0; i < ESNa1; i++){
	    x[i]	=drcT[ki];
	    y[i]	=drsT[ki];
	    ki++;
	  }
	  Plot2d(1,a,x,ESNa1,6,0,4,0);
	  Plot2d(1,a,y,ESNa1,6,0,4,0);
	}
	CbFlush();
      }
#endif   
      nw		=n+ESNLt;
      while(nw < ESNt1){
	dR0T[nw]	=dR0T[n];
	dZ0T[nw]	=dZ0T[n];
	ji		=ESnAP*n;
	nw1		=ESnAP*nw;
	in		=ESNa1*n;
	m1		=ESNa1*nw;
	for(i=0; i < ESNa1; i++){
	  for(j=0; j < ESNp1; j++){
	    gX[nw1]	=gX[ji];
	    ji++;
	    nw1++;
	  }
	  dbT[m1]	=dbT[in];
	  dbT1a[m1]	=dbT1a[in];
	  dbT2a[m1]	=dbT2a[in];
	  in++;
	  m1++;
	}
	ki	=K;
	m1	=ESnAF*nw;
	for(k=0; k < ESFp1; k++){
	  for(i=0; i < ESNa1; i++){
	    drcT[m1]	=drcT[ki];
	    drcT1a[m1]	=drcT1a[ki];
	    drcT2a[m1]	=drcT2a[ki];
	    drsT[m1]	=drsT[ki];
	    drsT1a[m1]	=drsT1a[ki];
	    drsT2a[m1]	=drsT2a[ki];
	    ki++;
	    m1++;
	  }
	}
	nw	+=ESNLt;
      }
    }
  }
#ifdef H
  {
    double gh0,*pgY,ghmx,ghmn,gYmx,gYmn,x[4],y[2],dy,dx;
    ghmx=0.;
    ghmn=0.;
    gYmx=0.;
    gYmn=0.;
    in	=iA;
    for(n=0; n < ESNt1; n++){
      ji	=ESnAP*n+ESNp1*iA;
      gh0	=GH[ji];
      if(ghmx < gh0)
	ghmx=gh0;
      if(ghmn > gh0)
	ghmn=gh0;
      pgY	=gX+ji;
      for(j=0; j < ESNp1; j++){
	ss	=pgY[j];
	if(gYmx < ss)
	  gYmx=ss;
	if(gYmn > ss)
	  gYmn=ss;
	ji++;
      }
      in	+=ESNa1;
    }
    ghmx	+=EZc2gp;
    x[0]	=ghmn*EZcr2gp;
    x[1]	=ghmx*EZcr2gp;
    dy		=0.1*(gYmx-gYmn);
    y[0]	=gYmn;
    y[1]	=dy*ESNt+gYmx;
   
    Scale2D(2,2,x,y,2,6,2);
    in	=iA;
    for(n=0; n < ESNt1; n++){
      ji	=ESnAP*n+ESNp1*iA;
      pgY	=gX+ji;
      EZz0	=dy*n;
      for(j=0; j < ESNp1; j++){
	ESgh[j]	=(ESgt[j]+GH[ji])*EZcr2gp;
	ginf[j]	=EZz0+pgY[j];
	ji++;
      }
      Plot2d(2,ESgh,ginf,ESNp1,6,0,0,0);
      in	+=ESNa1;
    }
    
    gh0		=ESgi[iA]/EZz0;
    y[0]	=0.;
    y[1]	=0.;
    Plot2d(1,x,y,2,6,0,14,0);
    y[0]	=EZz0;
    y[1]	=EZz0;
    Plot2d(1,x,y,2,6,0,14,0);

    i	=x[0]/0.2;
    dx	=0.2*i;
    y[0]	=gYmn;
    x[2]	=0.2*i+gh0*y[0];
    y[1]	=EZz0+gYmx;
    x[3]	=0.2*i+gh0*y[1];

    while((x[0] <= x[2] && x[2] <= x[1]) || (x[0] <= x[3] && x[3] <= x[1])){
      if(x[2] < x[0]){
	x[2]	=x[0];
	y[0]	=(x[2]-0.2*i)/gh0;
      } 
      if(x[2] > x[1]){
	x[2]	=x[1];
	y[0]	=(x[2]-0.2*i)/gh0;
      } 
      if(x[3] < x[0]){
	x[3]	=x[0];
	y[1]	=(x[3]-0.2*i)/gh0;
      } 
      if(x[3] > x[1]){
	x[3]	=x[1];
	y[1]	=(x[3]-0.2*i)/gh0;
      } 
      Plot2d(2,x+2,y,2,6,0,0,0);
      i++;
      y[0]	=gYmn;
      x[2]	=0.2*i+gh0*y[0];
      y[1]	=EZz0+gYmx;
      x[3]	=0.2*i+gh0*y[1];
    }
  }
#endif
  if(gXmx > ESsa[1]){
    s		=EZcr2*ESsa[1]/gXmx;
    in		=0;
    ji		=0.;
    for(n=0; n < ESNt1; n++){
      dR0T[n]	*=s;
      dZ0T[n]	*=s;
      for(i=0; i < ESNa1; i++){
	dbT[in]		*=s;
	in++;
      }
      ki	=ESnAF*n;
      for(k=0; k < ESFp1; k++){
	for(i=0; i < ESNa1; i++){
	  drcT[ki]	*=s;
	  drcT1a[ki]	*=s;
	  drcT2a[ki]	*=s;
	  drsT[ki]	*=s;
	  drsT1a[ki]	*=s;
	  drsT2a[ki]	*=s;
	  ki++;
	}
      }
    }
  }
  ji	=0;
  in	=0;
  for(n=0; n < ESNLt; n++){
    K		=ESnAF*n;
    ESaR0[n]	+=dR0T[n];
    ESaZ0[n]	+=dZ0T[n];
    in		=ESNa1*n;
    for(i=0; i < ESNa1; i++){
      ESsb[in]	+=dbT[in];
      ESsb1a[in]	+=dbT1a[in];
      ESsb2a[in]	+=dbT2a[in];
      in++;
    }
    ki		=K;
    for(k=0; k < ESFp1; k++){
      for(i=0; i < ESNa1; i++){
	rcT[ki]		+=drcT[ki];
	rcT1a[ki]	+=drcT1a[ki];
	rcT2a[ki]	+=drcT2a[ki];
	rsT[ki]		+=drsT[ki];
	rsT1a[ki]	+=drsT1a[ki];
	rsT2a[ki]	+=drsT2a[ki];
	ki++;
      }
    }
    ESRealGeom(n);

    nw	=n+ESNLt;
    while(nw < ESNt1){
      ESaR0[nw]	=ESaR0[n];
      ESaZ0[nw]	=ESaZ0[n];
      in	=ESNa1*n;
      m1	=ESNa1*nw;
      for(i=0; i < ESNa1; i++){
	ESsb[m1]	=ESsb[in];
	ESsb1a[m1]	=ESsb1a[in];
	ESsb2a[m1]	=ESsb2a[in];
	in++;
	m1++;
      }
      ki	=K;
      m1	=ESnAF*nw;
      for(k=0; k < ESFp1; k++){
	for(i=0; i < ESNa1; i++){
	  rcT[m1]	=rcT[ki];
	  rcT1a[m1]	=rcT1a[ki];
	  rcT2a[m1]	=rcT2a[ki];
	  rsT[m1]	=rsT[ki];
	  rsT1a[m1]	=rsT1a[ki];
	  rsT2a[m1]	=rsT2a[ki];
	  ki++;
	  m1++;
	}
      }
      nw	+=ESNLt;
    }
  }

  in	=ESNa1*Nerr;
  for(i=0; i < ESNa1; i++){
    printf("gF[%2d]=%13.6e gFa=%13.6e gY=%10.3e gYa=%10.3e %10.3e\n",i
	   ,EStgF[in+i],EStgF1a[in+i],EStgY[i],ESgY01a[i],ESgi[i]);
  }
  printf("??? err[%3d]=%10.3e\n",Nerr,err);
  printf("??? Resonant harmonics gYmx=%10.3e gXmx=%10.3e\n",gYmx,gXmx);
  for(kR=0; kR < ngYR; kR++){
    s	=NgYR[kR]/((double)MgYR[kR]);
    printf("??? gi_s[%2d][%2d]=%10.3e\n",MgYR[kR],NgYR[kR],s);
#ifdef J
    iR	=ngYR*kR;
    for(i=0; i < ESNa1; i++){
      printf("??? gYc[%2d]=%10.3e %10.3e\n",i,gYRc[iR],gYRs[iR]);
      iR++;
    }
#endif
  }
  return(0);
}

int ESgP2gF(double*gc,double*gs,double*g,int M)
{
  int m,m1,im,j,k,Np05;
  double s,c;
  
  Np05	=ESNp/2;
  c	=0.;
  for(j=0; j < ESNp; j++){
    c	+=g[j];
  }
  im	=0;
  gc[im]=c*crHp;
  gs[im]=0.;
  im	+=ESNa1;
  m1	=1;
  for(m=0; m < M; m++){
    c	=0.;
    s	=0.;
    k	=0;	
    for(j=0; j < ESNp; j++){
      c	+=g[j]*EScs1[k];
      s	+=g[j]*ESsn1[k];
      k	+=m1;
      if(k >= ESNp){
	k -=ESNp;
      }
    }
    gc[im]	=c*crHp;
    gs[im]	=s*crHp;
    if(m1 == Np05){
      gc[im]	*=EZcr2;
      gs[im]	*=EZcr2;
    }
    im		+=ESNa1;
    m1++;
  }
  return(0);
}

int ESP2F(double*gc,double*gs,double*g,int M)
{
  int m,M1,j,k,Np05;
  double s,c;
  
  M1	=M+1;
  Np05	=ESNp/2;
  c	=0.;
  for(j=0; j < ESNp; j++){
    c	+=g[j];
  }
  m	=0;
  gc[m]	=c*crHp;
  gs[m]	=0.;
  for(m=1; m < M1; m++){
    c	=0.;
    s	=0.;
    k	=0;	
    for(j=0; j < ESNp; j++){
      c	+=g[j]*EScs1[k];
      s	+=g[j]*ESsn1[k];
      k	+=m;
      if(k >= ESNp){
	k -=ESNp;
      }
    }
    gc[m]	=c*crHp;
    gs[m]	=s*crHp;
    if(m == Np05){
      gc[m]	*=EZcr2;
      gs[m]	*=EZcr2;
    }
  }
  return(0);
}

#endif
#endif



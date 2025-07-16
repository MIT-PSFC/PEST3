#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "esc.h"

extern clock_t clock(void);

extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern int ESmemlev;
extern int ESEqSolvInPr,ESEqSolvInCr;

extern int ESNa,ESNa1;
extern double *ESsa;
extern int ESNp,ESNp1;
extern int ES2Mp,ES2Mp1,ESMp,ESMp1;
extern int ESFp1;

extern double *ESsq,*ESsq1a,*ESsq2a
,*ESgm,*ESgm1a,*ESgm2a
,*ESPs,*ESPs1a,*ESPs2a
,*ESaT,*ESaT1a,*ESaT2a
,*ESFF,*ESFF1a,*ESFF2a
,*ESaF,*ESaF1a,*ESaF2a;

extern double *ESDPc,*ESDPc1,*ESDPc2,*ESDPs,*ESDPs1,*ESDPs2;
extern double *ESLc,*ESLs,*ESVc,*ESVs;
extern double *ESLc1,*ESLs1,*ESVc1,*ESVs1;
extern double *ESLc2,*ESLs2,*ESVc2,*ESVs2;
extern double *ESg22c,*ESg22s,*ESg22c1,*ESg22s1,*ESg22c2,*ESg22s2;

extern double *ESg12c,*ESg12s,*ESg12c2,*ESg12s2;
extern double *ESg11c,*ESg11s,*ESg11c2,*ESg11s2;
extern double *EScs1,*ESsn1,*EZcs2,*EZsn2;
extern double *EZdra,*EZdrgt,*EZdza,*EZdzgt;
extern double *ESaR0,*ESsb,*ESsb1a,*ESsb2a
,*rcT,*rsT,*rcT1a,*rsT1a,*rsT2a,*rcT2a;

extern double *ESGaac,*ESGaas,*ESGaac2,*ESGaas2;
extern double *ESGaqc,*ESGaqs,*ESGaqc2,*ESGaqs2;
extern double *ESGqqc,*ESGqqs,*ESGqqc2,*ESGqqs2;
extern double *ESGazc,*ESGazs,*ESGazc2,*ESGazs2;
extern double *ESGqzc,*ESGqzs,*ESGqzc2,*ESGqzs2;
extern double *ESGzzc,*ESGzzs,*ESGzzc2,*ESGzzs2;
extern double *ESsVc,*ESsVs,*ESsVc2,*ESsVs2;
extern int ESStN,ESStM1,ESStM2;
extern int ESStNRec,ESStMRec,ESStKRec;

double ESsAVc[500],ESsAVs[500],ESsAVc2[500],ESsAVs2[500];

static clock_t ESst0;
static double rR0,R0;

static int dM1=6,dM2=6,M1,M2,Mmx,N=1,Mp1,MMp1,M4Mp1;
static int nRes=0,Mres=0,iRes,kRes,*mRes;
static double qmn,qmx,PM1,PM2;
static double *aRes,*qRes,*gmRes,*dgmRes,aResChck;

static double *gyr,*gyi,*aYr,*aYi,gjr,*gji,*Yjr,*Yji;
static double *dgyr,*dgyi,*daYr,*daYi,*dgji,*dYjr,*dYji;

static double *gaar,*gaai,*gaqr,*gaqi,*gqqr,*gqqi
,*gazr,*gazi,*gqzr,*gqzi,*gzzr,*gzzi
,*svr,*svi,*dsvr,*dsvi;
static double *ar,*ai,*Ar,*Ai,*Br,*Bi,*Cr,*Ci,*Kr,*Ki,*Fr,*Fi,*Gr,*Gi
,*Tr,*X1r,*X1i,*X2r,*X2i;
static int *ind;
static int FlStChol=0;

static int kCheck=0;
static double ESsgnX,ESsgnY,gymr,gymi,aYmr,aYmi,ESdet;
static int FlagRes;

static int neq=0;
static double hcur,xout,reltol=1.0e-6,abstol=1e-12;
static double *rwork=NULL,*yh,*vatol;
static int itol=2,itask,istate=1,iopt=0,mf=10,lrw,liw=20;
static int *iwork=NULL;

static int NA2Mp1=0;

static int Mgy,MaY;

static FILE *fpy,*fpt;
static int NRec,MRec,KRec;

int ESReInitMHDMetrTnsr()
{
  int n;
  n	=ESNa1*ES2Mp1;
  if(NA2Mp1 < n){
    if(ReInitArray((void**)&ESGaac,28*NA2Mp1,28*n,sizeof(double)) < 0){
      FailureAlarm((char*)ESGaac,"ESReInitMHDMetrTnsr - no memory for ESGaac");
      return(-1);
    }
    NA2Mp1	=n;
    ESGaac2	=ESGaac	+NA2Mp1;
    ESGaas	=ESGaac2+NA2Mp1;
    ESGaas2	=ESGaas	+NA2Mp1;
    ESGaqc	=ESGaas2+NA2Mp1;
    ESGaqc2	=ESGaqc	+NA2Mp1;
    ESGaqs	=ESGaqc2+NA2Mp1;
    ESGaqs2	=ESGaqs	+NA2Mp1;
    ESGqqc	=ESGaqs2+NA2Mp1;
    ESGqqc2	=ESGqqc	+NA2Mp1;
    ESGqqs	=ESGqqc2+NA2Mp1;
    ESGqqs2	=ESGqqs	+NA2Mp1;
    ESGazc	=ESGqqs2+NA2Mp1;
    ESGazc2	=ESGazc	+NA2Mp1;
    ESGazs	=ESGazc2+NA2Mp1;
    ESGazs2	=ESGazs	+NA2Mp1;
    ESGqzc	=ESGazs2+NA2Mp1;
    ESGqzc2	=ESGqzc	+NA2Mp1;
    ESGqzs	=ESGqzc2+NA2Mp1;
    ESGqzs2	=ESGqzs	+NA2Mp1;
    ESGzzc	=ESGqzs2+NA2Mp1;
    ESGzzc2	=ESGzzc	+NA2Mp1;
    ESGzzs	=ESGzzc2+NA2Mp1;
    ESGzzs2	=ESGzzs	+NA2Mp1;
    ESsVc	=ESGzzs2+NA2Mp1;
    ESsVc2	=ESsVc	+NA2Mp1;
    ESsVs	=ESsVc2	+NA2Mp1;
    ESsVs2	=ESsVs	+NA2Mp1;

    ESmemlev	|=0x04000000;
  }
  return(0);
}

int ESDeInitMHDMetrTnsr()
{
  if(NA2Mp1){
    free(ESGaac);
    NA2Mp1=0;
  }
  return(0);
}

int ESMHDMetrTnsr()
{
  int i,j,k,ki,kj;
  double r0,ra;
  double *g11,*g12,*g22,*g13,*g23,*g33;
  double *dg11,*dg12,*dg22,*dg13,*dg23,*dg33;
  double F,Fa,H,Ha,Haa,q,qa,qaa;
  double *Lc,*Ls,*dLc,*dLs,*d2Lc,*d2Ls,rL0;

  Lc	=ESsVc;
  Ls	=ESsVs;
  dLc	=ESsVc2;
  dLs	=ESsVs2;
  d2Lc	=dLc+ES2Mp1;
  d2Ls	=dLs+ES2Mp1;

  g11	=EZdra	+ESNp1;
  g12	=g11	+ESNp1;
  g22	=g12	+ESNp1;
  dg11	=g22	+ESNp1;
  dg12	=dg11	+ESNp1;
  dg22	=dg12	+ESNp1;
  
  g13	=EZdza	+ESNp1;
  g23	=g13	+ESNp1;
  g33	=g23	+ESNp1;
  dg13	=g33	+ESNp1;
  dg23	=dg13	+ESNp1;
  dg33	=dg23	+ESNp1;

  q	=ESsq[0];
  /* metric coefficients and their radial derivatives at a=0 */
  rL0	=1./ESLc[0];
  H	=q*rL0;
  ki	=ESNa1;
  Lc[1]	=H*ESLc1[ki];
  Ls[1]	=H*ESLs1[ki];
  ki	+=ESNa1;
  Lc[2]	=H*EZcr2*ESLc2[ki];
  Ls[2]	=H*EZcr2*ESLs2[ki];
  F	=R0*rL0;
  H	=-F/q;
  {
    int k1,k2,k3;
    double t1c,t1s,a1c,a1s,ta2c,ta2s,aa2c,aa2s;
    t1c		=Lc[1];
    t1s		=Ls[1];
    a1c		=-Ls[1];
    a1s		=Lc[1];
    ta2c	=Lc[2];
    ta2s	=Ls[2];
    aa2c	=-Ls[2];
    aa2s	=Lc[2];
#ifdef H
    g33	=F;
    g13	=a1*F;
    g23	=t1*F;
    g11	=a1*a1*F;
    g12	=a1*t1*F;
    g22	=t1*t1*F;
    dg33	=t1*H;
    dg13	=aa2*F+a1*t1*H;
    dg23	=ta2*F+t1*t1*H;
    dg11	=2.*aa2*a1*F+a1*a1*t1*H;
    dg12	=(aa2*t1+ta2*a1)*F+a1*t1*t1*H;
    dg22	=2.*ta2*t1*F+t1*t1*t1*H;
#endif
    k1	=ESNa1;
    k2	=k1+ESNa1;
    k3	=k2+ESNa1;

    ESGzzc[0]	=F;
    ESGzzs[0]	=0.;
    ESGazc[k1]	=a1c*F;
    ESGazs[k1]	=a1s*F;
    ESGqzc[k1]	=t1c*F;
    ESGqzs[k1]	=t1s*F;

    ESGaac[0]	=2.*(a1c*a1c+a1s*a1s)*F;
    ESGaas[0]	=0.;
    ESGaac[k2]	=(a1c*a1c-a1s*a1s)*F;
    ESGaas[k2]	=2.*a1c*a1s*F;

    ESGaqc[0]	=2.*(a1c*t1c+a1s*t1s)*F;
    ESGaqs[0]	=0.;
    ESGaqc[k2]	=(a1c*t1c-a1s*t1s)*F;
    ESGaqs[k2]	=(a1c*t1s+a1s*t1c)*F;

    ESGqqc[0]	=2.*(t1c*t1c+t1s*t1s)*F;
    ESGqqs[0]	=0.;
    ESGqqc[k2]	=(t1c*t1c-t1s*t1s)*F;
    ESGqqs[k2]	=2.*t1c*t1s*F;

    ESGzzc2[k1]	=t1c*H;
    ESGzzs2[k1]	=t1s*H;

    ESGazc2[0]	=2.*(a1c*t1c+a1s*t1s)*H;
    ESGazs2[0]	=0.;
    ESGazc2[k2]	=aa2c*F+(a1c*t1c-a1s*t1s)*H;
    ESGazs2[k2]	=aa2s*F+(a1c*t1s+a1s*t1c)*H;

    ESGqzc2[0]	=2.*(t1c*t1c+t1s*t1s)*H;
    ESGqzs2[0]	=0.;
    ESGqzc2[k2]	=ta2c*F+(t1c*t1c-t1s*t1s)*H;
    ESGqzs2[k2]	=ta2s*F+2.*t1c*t1s*H;

    ESGaac2[k1]	=2.*(aa2c*a1c+aa2s*a1s)*F+
      ((3.*a1c*a1c+a1s*a1s)*t1c+2.*a1c*a1s*t1s)*H;
    ESGaas2[k1]	=2.*(aa2s*a1c-aa2c*a1s)*F+
      ((a1c*a1c+3.*a1s*a1s)*t1s+2.*a1c*a1s*t1c)*H;
    ESGaac2[k3]	=2.*(aa2c*a1c-aa2s*a1s)*F+
      ((a1c*a1c-a1s*a1s)*t1c-2.*a1c*a1s*t1s)*H;
    ESGaas2[k3]	=2.*(aa2s*a1c+aa2c*a1s)*F+
      ((a1c*a1c-a1s*a1s)*t1s+2.*a1c*a1s*t1c)*H;

    ESGaqc2[k1]	=(aa2c*t1c+aa2s*t1s+ta2c*a1c+ta2s*a1s)*F+
      ((3.*t1c*t1c+t1s*t1s)*a1c+2.*t1c*t1s*a1s)*H;
    ESGaqs2[k1]	=(aa2s*t1c-aa2c*t1s+ta2s*a1c-ta2c*a1s)*F+
      ((t1c*t1c+3.*t1s*t1s)*a1s+2.*t1c*t1s*a1c)*H;
    ESGaqc2[k3]	=(aa2c*t1c-aa2s*t1s+ta2c*a1c-ta2s*a1s)+
      ((t1c*t1c-t1s*t1s)*a1c-2.*t1c*t1s*a1s)*H;
    ESGaqs2[k3]	=(aa2s*t1c+aa2c*t1s+ta2s*a1c+ta2c*a1s)+
      ((t1c*t1c-t1s*t1s)*a1s+2.*t1c*t1s*a1c)*H;

    ESGqqc2[k1]	=2.*(ta2c*t1c+ta2s*t1s)*F+3.*(t1c*t1c+t1s*t1s)*t1c*H;
    ESGqqs2[k1]	=2.*(ta2s*t1c-ta2c*t1s)*F+3.*(t1c*t1c+t1s*t1s)*t1s*H;
    ESGqqc2[k3]	=2.*(ta2c*t1c-ta2s*t1s)*F+t1c*(t1c*t1c-3.*t1s*t1s)*H;
    ESGqqs2[k3]	=2.*(ta2s*t1c+ta2c*t1s)*F+t1s*(3.*t1c*t1c-t1s*t1s)*H;
  }
  ki	=0;
  for(k=0; k < ES2Mp1; k++){
    if(k != 0){
      ESGzzc[ki]	=0.;
      ESGzzs[ki]	=0.;
      if(k != 2){
	ESGaac[ki]	=0.;
	ESGaas[ki]	=0.;
	ESGaqc[ki]	=0.;
	ESGaqs[ki]	=0.;
	ESGqqc[ki]	=0.;
	ESGqqs[ki]	=0.;
	ESGazc2[ki]	=0.;
	ESGazs2[ki]	=0.;
	ESGqzc2[ki]	=0.;
	ESGqzs2[ki]	=0.;
      }
    }
    if(k != 1){
      ESGazc[ki]	=0.;
      ESGazs[ki]	=0.;
      ESGqzc[ki]	=0.;
      ESGqzs[ki]	=0.;
      ESGzzc2[ki]	=0.;
      ESGzzs2[ki]	=0.;
      if(k != 3){
	ESGaac2[ki]	=0.;
	ESGaas2[ki]	=0.;
	ESGaqc2[ki]	=0.;
	ESGaqs2[ki]	=0.;
	ESGqqc2[ki]	=0.;
	ESGqqs2[ki]	=0.;
      }
    }
    ki	+=ESNa1;
  }
  for(i=1; i < ESNa; i++){
    q		=ESsq[i];
    qa		=ESsq1a[i];
    Lc[0]	=ESLc[i];
    rL0		=1./Lc[0];
    ra		=q*rL0/ESsa[i];
    dLc[0]	=-ESLc1[i]*rL0;
    ki		=i;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      Lc[k]	=2.*ESLc[ki];
      Ls[k]	=2.*ESLs[ki];
      dLc[k]	=2.*ESLc1[ki];
      dLs[k]	=2.*ESLs1[ki];
    }
    for(j=0; j < ESNp; j++){
      F		=0.;
      H		=0.;
      Ha	=0.;
      kj	=0;
      for(k=1; k < ES2Mp1; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	H	+=(Lc[k]*ESsn1[kj]-Ls[k]*EScs1[kj])/k;
	Ha	+=(dLc[k]*ESsn1[kj]-dLs[k]*EScs1[kj])/k;
	F	+=Lc[k]*EScs1[kj]+Ls[k]*ESsn1[kj];
      }
      H		*=rL0;
      Ha	*=rL0;
      g33[j]	=R0/(Lc[0]+F);
      H		=qa*H+q*(Ha+H*dLc[0]);
      F		*=ra;
      g13[j]	=H*g33[j];
      g23[j]	=F*g33[j];
      g11[j]	=H*g13[j];
      g12[j]	=H*g23[j];
      g22[j]	=F*g23[j];
    }
    g11[j]	=g11[0];
    g12[j]	=g12[0];
    g22[j]	=g22[0];
    g13[j]	=g13[0];
    g23[j]	=g23[0];
    g33[j]	=g33[0];
    ESgP2gF(ESGaac+i,ESGaas+i,g11,ES2Mp);
    ESgP2gF(ESGaqc+i,ESGaqs+i,g12,ES2Mp);
    ESgP2gF(ESGqqc+i,ESGqqs+i,g22,ES2Mp);
    ESgP2gF(ESGazc+i,ESGazs+i,g13,ES2Mp);
    ESgP2gF(ESGqzc+i,ESGqzs+i,g23,ES2Mp);
    ESgP2gF(ESGzzc+i,ESGzzs+i,g33,ES2Mp);
  }
  /* Metric coefficients and their derivatives at a=1 */
  i		=ESNa;
  q		=ESsq[i];
  qa		=ESsq1a[i];
  qaa		=ESsq2a[i];
  Lc[0]		=ESLc[i];
  rL0		=1./Lc[0];
  dLc[0]	=-ESLc1[i]*rL0;
  d2Lc[0]	=2.*dLc[0]*dLc[0]-ESLc2[i]*rL0;
  ki		=i;
  for(k=1; k < ES2Mp1; k++){
    ki		+=ESNa1;
    Lc[k]	=2.*ESLc[ki];
    Ls[k]	=2.*ESLs[ki];
    dLc[k]	=2.*ESLc1[ki];
    dLs[k]	=2.*ESLs1[ki];
    d2Lc[k]	=2.*ESLc2[ki];
    d2Ls[k]	=2.*ESLs2[ki];
  }
  for(j=0; j < ESNp; j++){
    H	=0.;
    Ha	=0.;
    Haa	=0.;
    F	=0.;
    Fa	=0.;
    kj	=0;
    for(k=1; k < ES2Mp1; k++){
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      H		+=(Lc[k]*ESsn1[kj]-Ls[k]*EScs1[kj])/k;
      Ha	+=(dLc[k]*ESsn1[kj]-dLs[k]*EScs1[kj])/k;
      Haa	+=(d2Lc[k]*ESsn1[kj]-d2Ls[k]*EScs1[kj])/k;
      F		+=Lc[k]*EScs1[kj]+Ls[k]*ESsn1[kj];
      Fa	+=dLc[k]*EScs1[kj]+dLs[k]*ESsn1[kj];
    }	
    ra	=1./(Lc[0]+F);
    g33[j]	=R0*ra;
    dg33[j]	=-g33[j]*(ESLc1[i]+Fa)*ra;
    H	*=rL0;
    Ha	*=rL0;
    Haa	*=rL0;
    F	*=rL0; 
    Fa	*=rL0; 
    Fa	=qa*F+q*(Fa+F*(dLc[0]-1.));
    F	*=q;
    ra	=Ha+H*dLc[0];
    Haa	=q*(Haa+2.*Ha*dLc[0]+H*d2Lc[0])+2.*qa*ra+qaa*H;
    Ha	=q*ra+qa*H;
    g13[j]	=Ha*g33[j];
    g23[j]	=F*g33[j];
    g11[j]	=Ha*g13[j];
    g12[j]	=Ha*g23[j];
    g22[j]	=F*g23[j];
    
    dg13[j]	=Haa*g33[j]+Ha*dg33[j];
    dg23[j]	=Fa*g33[j]+F*dg33[j];
    dg11[j]	=Haa*g13[j]+Ha*dg13[j];
    dg12[j]	=Haa*g23[j]+Ha*dg23[j];
    dg22[j]	=Fa*g23[j]+F*dg23[j];
  }
  g11[j]	=g11[0];
  g12[j]	=g12[0];
  g22[j]	=g22[0];
  g13[j]	=g13[0];
  g23[j]	=g23[0];
  g33[j]	=g33[0];
  dg11[j]	=dg11[0];
  dg12[j]	=dg12[0];
  dg22[j]	=dg22[0];
  dg13[j]	=dg13[0];
  dg23[j]	=dg23[0];
  dg33[j]	=dg33[0];
  ESgP2gF(ESGaac+i,ESGaas+i,g11,ES2Mp);
  ESgP2gF(ESGaqc+i,ESGaqs+i,g12,ES2Mp);
  ESgP2gF(ESGqqc+i,ESGqqs+i,g22,ES2Mp);
  ESgP2gF(ESGazc+i,ESGazs+i,g13,ES2Mp);
  ESgP2gF(ESGqzc+i,ESGqzs+i,g23,ES2Mp);
  ESgP2gF(ESGzzc+i,ESGzzs+i,g33,ES2Mp);
  ESgP2gF(ESGaac2+i,ESGaas2+i,dg11,ES2Mp);
  ESgP2gF(ESGaqc2+i,ESGaqs2+i,dg12,ES2Mp);
  ESgP2gF(ESGqqc2+i,ESGqqs2+i,dg22,ES2Mp);
  ESgP2gF(ESGazc2+i,ESGazs2+i,dg13,ES2Mp);
  ESgP2gF(ESGqzc2+i,ESGqzs2+i,dg23,ES2Mp);
  ESgP2gF(ESGzzc2+i,ESGzzs2+i,dg33,ES2Mp);

  /* splining Fourier coefficients */
  i	=0;
  j	=ESNa;
  H	=ESGaac2[i];
  F	=ESGaac2[j];
  splA(ESGaac+i,ESGaac2+i,&H,&F);
  H	=ESGaqc2[i];
  F	=ESGaqc2[j];
  splA(ESGaqc+i,ESGaqc2+i,&H,&F);
  H	=ESGqqc2[i];
  F	=ESGqqc2[j];
  splA(ESGqqc+i,ESGqqc2+i,&H,&F);
  H	=ESGazc2[i];
  F	=ESGazc2[j];
  splA(ESGazc+i,ESGazc2+i,&H,&F);
  H	=ESGqzc2[i];
  F	=ESGqzc2[j];
  splA(ESGqzc+i,ESGqzc2+i,&H,&F);
  H	=ESGzzc2[i];
  F	=ESGzzc2[j];
  splA(ESGzzc+i,ESGzzc2+i,&H,&F);
  for(k=1; k < ES2Mp1; k++){
    i 	+=ESNa1;
    j	+=ESNa1;
    H	=ESGaac2[i];
    F	=ESGaac2[j];
    splA(ESGaac+i,ESGaac2+i,&H,&F);
    H	=ESGaas2[i];
    F	=ESGaas2[j];
    splA(ESGaas+i,ESGaas2+i,&H,&F);
    H	=ESGaqc2[i];
    F	=ESGaqc2[j];
    splA(ESGaqc+i,ESGaqc2+i,&H,&F);
    H	=ESGaqs2[i];
    F	=ESGaqs2[j];
    splA(ESGaqs+i,ESGaqs2+i,&H,&F);
    H	=ESGqqc2[i];
    F	=ESGqqc2[j];
    splA(ESGqqc+i,ESGqqc2+i,&H,&F);
    H	=ESGqqs2[i];
    F	=ESGqqs2[j];
    splA(ESGqqs+i,ESGqqs2+i,&H,&F);
    H	=ESGazc2[i];
    F	=ESGazc2[j];
    splA(ESGazc+i,ESGazc2+i,&H,&F);
    H	=ESGazs2[i];
    F	=ESGazs2[j];
    splA(ESGazs+i,ESGazs2+i,&H,&F);
    H	=ESGqzc2[i];
    F	=ESGqzc2[j];
    splA(ESGqzc+i,ESGqzc2+i,&H,&F);
    H	=ESGqzs2[i];
    F	=ESGqzs2[j];
    splA(ESGqzs+i,ESGqzs2+i,&H,&F);
    H	=ESGzzc2[i];
    F	=ESGzzc2[j];
    splA(ESGzzc+i,ESGzzc2+i,&H,&F);
    H	=ESGzzs2[i];
    F	=ESGzzs2[j];
    splA(ESGzzs+i,ESGzzs2+i,&H,&F);
  }
  for(i=0; i < ESNa1; i++){
    ESGaac[i]	+=ESg11c[i];
    ESGaqc[i]	+=ESg12c[i];
    ESGqqc[i]	+=ESg22c[i];
    ESGaac2[i]	+=ESg11c2[i];
    ESGaqc2[i]	+=ESg12c2[i];
    ESGqqc2[i]	+=ESg22c2[i];
  }
  ki	=ESNa1;
  for(k=1; k < ES2Mp1; k++){
    for(i=0; i < ESNa1; i++){
      ESGaac[ki]	+=ESg11c[ki];
      ESGaas[ki]	+=ESg11s[ki];
      ESGaqc[ki]	+=ESg12c[ki];
      ESGaqs[ki]	+=ESg12s[ki];
      ESGqqc[ki]	+=ESg22c[ki];
      ESGqqs[ki]	+=ESg22s[ki];
      ESGaac2[ki]	+=ESg11c2[ki];
      ESGaas2[ki]	+=ESg11s2[ki];
      ESGaqc2[ki]	+=ESg12c2[ki];
      ESGaqs2[ki]	+=ESg12s2[ki];
      ESGqqc2[ki]	+=ESg22c2[ki];
      ESGqqs2[ki]	+=ESg22s2[ki];
      ki++;
    }
  }
  /* Calculation of $\msv$ */
  r0	=R0*R0;
  i	=0;
  rL0	=r0/(ESLc[i]*ESaF[i]);
  ESsVc[i]	=ESLc[i]*rL0;
  ESsVc2[i]	=0.;
  ESsVs[i]	=0.;
  ESsVs2[i]	=0.;
  ki	=ESNa1;
  ESsVc[ki]	=0.;
  ESsVs[ki]	=0.;
  ESsVc2[ki]	=(ESVc1[ki]+ESLc1[ki])*rL0;
  ESsVs2[ki]	=(ESVs1[ki]+ESLs1[ki])*rL0;
  for(k=2; k < ES2Mp1; k++){
    ki		+=ESNa1;
    ESsVc[ki]	=0.;
    ESsVs[ki]	=0.;
    ESsVc2[ki]	=0.;
    ESsVs2[ki]	=0.;
  }
  for(i=1; i < ESNa1; i++){
    rL0		=r0/(ESLc[i]*ESaF[i]);
    ESsVc[i]	=(ESVc[i]+ESLc[i])*rL0;
    ESsVs[i]	=0.;
    ki	=i;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      ESsVc[ki]	=(ESVc[ki]+ESLc[ki])*rL0;
      ESsVs[ki]	=(ESVs[ki]+ESLs[ki])*rL0;
    }
  }
#ifdef H
  ki	=ESNa;
  H	=-ESLc1[ki]/ESLc[ki]-ESaF1a[ki]/ESaF[ki];
  F	=(ESVc1[ki]+ESLc1[ki])*rL0+ESsVc[ki]*H;
  for(i=0; i < ESNa1; i++){
    ESsAVc[i]	=ESsVc[i];
  }
  splXAA0(ESsAVc,ESsAVc2,&F);
  j	=0;
  for(k=1; k < ES2Mp1; k++){
    j	+=ESNa1;
    ki	+=ESNa1;
    Ha	=(ESVc1[ki]+ESLc1[ki])*rL0+ESsVc[ki]*H;
    F	=(ESVs1[ki]+ESLs1[ki])*rL0+ESsVs[ki]*H;
    if(k%2){
      splXAA1(ESsAVc+j,ESsAVc2+j,ESsVc+j,ESsVc2[j],&Ha);
      splXAA1(ESsAVs+j,ESsAVs2+j,ESsVs+j,ESsVs2[j],&F);
    }
    else{
      for(i=0; i < ESNa1; i++){
	ESsAVc[j+i]	=ESsVc[j+i];
	ESsAVs[j+i]	=ESsVs[j+i];
      }
      splXAA0(ESsAVc+j,ESsAVc2+j,&Ha);
      splXAA0(ESsAVs+j,ESsAVs2+j,&F);
    }
  }
#endif
  i	=ESNa;
  H	=-ESLc1[i]/ESLc[i]-ESaF1a[i]/ESaF[i];
  F	=(ESVc1[i]+ESLc1[i])*rL0+ESsVc[i]*H;
  splA(ESsVc,ESsVc2,ESsa,&F);
  j	=0;
  ki	=i;
  for(k=1; k < ES2Mp1; k++){
    j	+=ESNa1;
    ki	+=ESNa1;
    Ha	=ESsVc2[j];
    F	=(ESVc1[ki]+ESLc1[ki])*rL0+ESsVc[ki]*H;
    splA(ESsVc+j,ESsVc2+j,&Ha,&F);
    Ha	=ESsVs2[j];
    F	=(ESVs1[ki]+ESLs1[ki])*rL0+ESsVs[ki]*H;
    splA(ESsVs+j,ESsVs2+j,&Ha,&F);
  }

  return(0);
}

int ESInitStLSODE(int n)
{
  lrw	=liw+16*n;
  iwork	=(int*)malloc(liw*sizeof(int));
  if(iwork == NULL){
    printf("Failure in memory allocation for iwork in ESReInitStLSODE()\n");
    printf("Too many harmonics\n");
    return(-1);
  }
  rwork	=(double*)malloc((lrw+n)*sizeof(double));
  if(rwork == NULL){
    printf("Failure in memory allocation for rwork in ESReInitStLSODE()\n");
    printf("Too many harmonics\n");
    free(iwork);
    return(-1);
  }
  yh	=rwork+liw;
  vatol	=rwork+lrw;
  return(0);
}

int ESDeInitStLSODE()
{
  free(rwork);
  rwork	=NULL;
  free(iwork);
  iwork	=NULL;
  return(0);
}

int ESStResonances()
{
  int i,i1,m,M,it;
  double h,H,a,gm,q,q0,q1,a0,gm0,dgm0,d2gm0,a1,gm1,dgm1,d2gm1
    ,as,gms,dgms,d2gms,ax,gmx,dgmx,d2gmx;

  a1	=ESsa[0];
  gm1	=ESgm[0];
  dgm1	=0.;
  d2gm1	=ESgm2a[0];
  q	=N/gm1;
  m	=q;
  nRes	=0;
  if((double)m == q){
    nRes++;
  }
  for(i=0; i < ESNa; i++){
    i1	=i+1;
    a0	=a1;
    gm0	=gm1;
    dgm0=dgm1;
    d2gm0=d2gm1;
    a1	=ESsa[i1];
    gm1	=ESgm[i1];
    dgm1=ESgm1a[i1];
    d2gm1=ESgm2a[i1];
    if((d2gm0 <= 0 && d2gm1 <= 0) || (d2gm0 >= 0 && d2gm1 >= 0)){
      /* No \gm'' = 0 */
      as	=a1;
      gms	=gm1;
      dgms	=dgm1;
    }
    else{
      /* There is \gm'' = 0 */
      as	=a0+(a1-a0)*d2gm0/(d2gm0-d2gm1);
      ESSetSplA(as);
      splRA(&gms,&dgms,ESgm,ESgm2a);
    }
    while(a0 < a1){
      if((dgm0 <= 0 && dgms <= 0) || (dgm0 >= 0 && dgms >= 0)){
	/* No \gm' = 0 */
	ax	=as;
	gmx	=gms;
      }
      else{
	/* There is \gm' = 0 */
	ax	=EZcr2*(as+a0);
	it	=0;
	do{
	  ESSetSplA(ax);
	  splRA2(&gmx,&dgmx,&d2gmx,ESgm,ESgm2a);
	  h	=dgmx/d2gmx;
	  ax	-=h;
	  it++;
	}while(fabs(h) > 1e-6*(as-a0) && it < 20);
	if(it >= 20){
	  printf("Cannot find gm'(a)=0 in ESStResonances()\n");
	  printf("a=%11.3e %11.3e %11.3e %11.3e\n",ax,gmx,dgmx,d2gmx);
	  return(-1);
	}
	ESSetSplA(ax);
	splRA(&gmx,&dgmx,ESgm,ESgm2a);
	q0	=1./gmx;
	if(qmn > q0){
	  qmn	=q0;
	}
	if(qmx < q0){
	  qmx	=q0;
	}
      }
      while(a0 < as){
	q0	=N/gm0;
	q1	=N/gmx;
	m	=q0;
	if(q0 < q1){
	  m++;
	  while(m <= q1){
	    nRes++;
	    m++;
	  }
	}
	else{
	  while(q1 <= m){
	    nRes++;
	    m--;
	  }
	}
	a0	=ax;
	gm0	=gmx;
	ax	=as;
	gmx	=gms;
      }
      a0	=as;
      gm0	=gms;
      dgm0	=dgms;
      as	=a1;
      gms	=gm1;
      dgms	=dgm1;
    }
  }

  aRes		=(double*)malloc(4*nRes*sizeof(double));
  if(aRes == NULL){
    printf("Failure in memory allocation for aRes in ESStResonances()\n");
    return(-1);
  }
  qRes	=aRes	+nRes;
  gmRes	=qRes	+nRes;
  dgmRes=gmRes	+nRes;

  mRes		=(int*)malloc(nRes*sizeof(int));
  if(mRes == NULL){
    printf("Failure in memory allocation for mRes in ESStResonances()\n");
    return(-1);
  }

  a1	=ESsa[0];
  gm1	=ESgm[0];
  dgm1	=0.;
  d2gm1	=ESgm2a[0];
  q	=N/gm1;
  m	=q;
  nRes	=0;
  if((double)m == q){
    aRes[nRes]	=a1;
    qRes[nRes]	=q;
    mRes[nRes]	=m;
    gmRes[nRes]	=gm1;
    dgmRes[nRes]=dgm1;
    nRes++;
  }
  for(i=0; i < ESNa; i++){
    i1	=i+1;
    a0	=a1;
    gm0	=gm1;
    dgm0=dgm1;
    d2gm0=d2gm1;
    a1	=ESsa[i1];
    gm1	=ESgm[i1];
    dgm1=ESgm1a[i1];
    d2gm1=ESgm2a[i1];
    if((d2gm0 <= 0 && d2gm1 <= 0) || (d2gm0 >= 0 && d2gm1 >= 0)){
      /* No \gm'' = 0 */
      as	=a1;
      gms	=gm1;
      dgms	=dgm1;
    }
    else{
      /* There is \gm'' = 0 */
      as	=a0+(a1-a0)*d2gm0/(d2gm0-d2gm1);
      ESSetSplA(as);
      splRA(&gms,&dgms,ESgm,ESgm2a);
    }
    while(a0 < a1){
      if((dgm0 <= 0 && dgms <= 0) || (dgm0 >= 0 && dgms >= 0)){
	/* No \gm' = 0 */
	ax	=as;
	gmx	=gms;
      }
      else{
	/* There is \gm' = 0 */
	ax	=EZcr2*(as+a0);
	H	=1e-6*(as-a0);
	do{
	  ESSetSplA(ax);
	  splRA2(&gmx,&dgmx,&d2gmx,ESgm,ESgm2a);
	  h	=dgmx/d2gmx;
	  ax	-=h;
	}while(fabs(h) > H);
	ESSetSplA(ax);
	splRA(&gmx,&dgmx,ESgm,ESgm2a);
      }
      while(a0 < as){
	q0	=N/gm0;
	q1	=N/gmx;
	m	=q0;
	a	=EZcr2*(ax+a0);
	H	=1e-6*(ax-a0);
	if(q0 < q1){
	  m++;
	  while(m < q1){
	    gm	=((double)N)/m;
	    it	=0;
	    do{
	      ESSetSplA(a);
	      splRA(&gm0,&dgm0,ESgm,ESgm2a);
	      h		=(gm-gm0)/dgm0;
	      a		+=h;
	      if(a < a0){
		a	=a0;
	      }
	      if(a > ax){
		a	=ax;
	      }
	      it++;
	    }while(fabs(h) > H && it < 20);
	    if(it >= 20){
	      printf("Cannot find gm=%10.3e in ESStResonances()\n",gm);
	      printf("a=%11.3e %11.3e %11.3e\n",a,gm0,dgm0);
	      return(-1);
	    }
	    ESSetSplA(a);
	    splRA(gmRes+nRes,dgmRes+nRes,ESgm,ESgm2a);
	    aRes[nRes]	=a;
	    qRes[nRes]	=1./gmRes[nRes];
	    mRes[nRes]	=m;
	    nRes++;
	    m++;
	  }
	}
	else{
	  while(q1 < m){
	    gm	=((double)N)/m;
	    it	=0;
	    do{
	      ESSetSplA(a);
	      splRA(&gm0,&dgm0,ESgm,ESgm2a);
	      h		=(gm-gm0)/dgm0;
	      a		+=h;
	      if(a < a0){
		a	=a0;
	      }
	      if(a > ax){
		a	=ax;
	      }
	      it++;
	    }while(fabs(h) > H && it < 20);
	    if(it >= 20){
	      printf("?? Cannot find gm=%10.3e in ESStResonances()\n",gm);
	      printf("a=%11.3e %11.3e %11.3e\n",a,gm0,dgm0);
	      return(-1);
	    }
	    ESSetSplA(a);
	    splRA(gmRes+nRes,dgmRes+nRes,ESgm,ESgm2a);
	    aRes[nRes]	=a;
	    qRes[nRes]	=1./gmRes[nRes];
	    mRes[nRes]	=m;
	    nRes++;
	    m--;
	  }
	}
	if((double)m == q1){
	  ESSetSplA(ax);
	  splRA(gmRes+nRes,dgmRes+nRes,ESgm,ESgm2a);
	  aRes[nRes]	=ax;
	  qRes[nRes]	=1./gmx;
	  mRes[nRes]	=m;
	  nRes++;
	}
	a0	=ax;
	gm0	=gmx;
	ax	=as;
	gmx	=gms;
      }
      a0	=as;
      gm0	=gms;
      dgm0	=dgms;
      as	=a1;
      gms	=gm1;
      dgms	=dgm1;
    }
  }
  return(0);
}

int ESReInitMHDStabLowN()
{
  int i;

  MMp1	=Mp1*Mp1;
  M4Mp1	=4*MMp1;
  neq	=M4Mp1;

  if(ESInitStLSODE(neq) < 0){
    return(-1);
  }

  ar	=(double*)malloc((10*MMp1+5*Mp1)*sizeof(double));
  if(ar == NULL){
    printf("Failure in memory allocation for ar in ESReInitMHDStabLowN()\n");
    ESDeInitStLSODE();
    return(-1);
  }
  ai	=ar+MMp1;
  Ar	=ai+MMp1;
  Ai	=Ar+MMp1;
  Br	=Ai+MMp1;
  Bi	=Br+MMp1;
  Cr	=Bi+MMp1;
  Ci	=Cr+MMp1;
  Fr	=Ci+MMp1;
  Fi	=Fr+MMp1;
  Tr	=Fi+MMp1;
  X1r	=Tr+Mp1;
  X1i	=X1r+Mp1;
  X2r	=X1i+Mp1;
  X2i	=X2r+Mp1;

  Kr	=Ar;
  Ki	=Ai;
  Gr	=Br;
  Gi	=Bi;

  gyr		=(double*)malloc(8*MMp1*sizeof(double));
  if(gyr == NULL){
    printf("Failure in memory allocation for gyr in ESReInitMHDStabLowN()\n");
    free(ar);
    ESDeInitStLSODE();
    return(-1);
  }
  gyi	=gyr	+MMp1;
  aYr	=gyi	+MMp1;
  aYi	=aYr	+MMp1;
  dgyr	=aYi	+MMp1;
  dgyi	=dgyr	+MMp1;
  daYr	=dgyi	+MMp1;
  daYi	=daYr	+MMp1;

  ind		=(int*)malloc(2*Mp1*sizeof(int));
  if(ind == NULL){
    printf("Failure in memory allocation for ind in ESReInitMHDStabLowN()\n");
    free(gyr);
    free(ar);
    ESDeInitStLSODE();
    return(-1);
  }

  gaar		=(double*)malloc(16*ES2Mp1*sizeof(double));
  if(gaar == NULL){
    printf("Failure in memory allocation for gaar in ESReInitMHDStabLowN()\n");
    free(ind);
    free(gyr);
    free(ar);
    ESDeInitStLSODE();
    return(-1);
  }
  gaai	=gaar	+ES2Mp1;
  gaqr	=gaai	+ES2Mp1;
  gaqi	=gaqr	+ES2Mp1;
  gqqr	=gaqi	+ES2Mp1;
  gqqi	=gqqr	+ES2Mp1;
  gazr	=gqqi	+ES2Mp1;
  gazi	=gazr	+ES2Mp1;
  gqzr	=gazi	+ES2Mp1;
  gqzi	=gqzr	+ES2Mp1;
  gzzr	=gqzi	+ES2Mp1;
  gzzi	=gzzr	+ES2Mp1;
  svr	=gzzi	+ES2Mp1;
  svi	=svr	+ES2Mp1;
  dsvr	=svi	+ES2Mp1;
  dsvi	=dsvr	+ES2Mp1;

  ESmemlev |=0x08000000;
  return(0);
}

int ESDeInitMHDStabLowN()
{
  free(gaar);
  free(ind);
  free(gyr);
  free(ar);
  ESDeInitStLSODE();

  return(0);
}

int ESInitStFundSolution(int Opt,double x)
{
  double am,bK0,Er,Ei,A;
  double br,bi,sr,si,tr,ti,gdr,gdi;
  double *pr,*pi,*Pr,*Pi,*Yr,*Yi,*yr,*yi,*EZvr,*vi;
  int i,j,jj,k,kk,m,mm;

  bK0	=ESg22c[0];
  k	=2*ESNa1;
  Er	=ESg22c[k];
  Ei	=-ESg22s[k];

  gdi	=-ESsq[0]*x/ESLc[0];
  gdr	=ESLc1[ESNa1]*gdi;
  gdi	*=ESLs1[ESNa1];

  if(M1 <= 0){
    /* m=0 harmonic */
    j	=0;
    k	=(j-M1)*Mp1-M1;
    pr	=gyr+k;
    pi	=gyi+k;
    yr	=aYr+k;
    yi	=aYi+k;
    for(k=M1; k <= M2; k++){
      pr[k]	=0.;
      pi[k]	=0.;
      yr[k]	=0.;
      yi[k]	=0.;
    }
    pr[0]	=x*x;
    yr[0]	=2.*R0;
    /* initialization of vatol and rwork */
    A	=reltol*pr[0];
    kk	=Mp1*(-M1);
    Pr	=vatol+kk;
    Pi	=Pr+MMp1;
    Yr	=Pi+MMp1;
    Yi	=Yr+MMp1;
    for(k=0; k < Mp1; k++){
      *Pr++	=A;
      *Pi++	=A;
      *Yr++	=A;
      *Yi++	=A;
    }
    for(i=0; i < 13; i++){
      Pr	=yh+neq*i+kk;
      Pi	=Pr+MMp1;
      Yr	=Pi+MMp1;
      Yi	=Yr+MMp1;
      for(k=0; k < Mp1; k++){
	*Pr++	=0.;
	*Pi++	=0.;
	*Yr++	=0.;
	*Yi++	=0.;
      }
      if(i == 0){
	yh[kk-M1]	=pr[0];
	yh[kk-M1+2*MMp1]=yr[0];
      }
      if(i == 1){
	yh[kk-M1]	=2.*x;
	yh[kk-M1+2*MMp1]=yr[0];
      }
      if(i == 2){
	yh[kk-M1]	=1.;
      }
    }
    /* M1 <= m=jj and  m=j < |M1| harmonics */
    m		=-M1;
    am		=1.;
    while(j < m){
      j++;
      am	*=x;
      if(am < 1e-100){
	am	=1e-100;
      }
      jj	=-j;
      k		=(j-M1)*Mp1-M1;
      pr	=gyr+k;
      pi	=gyi+k;
      EZvr	=aYr+k;
      vi	=aYi+k;
      k		=(jj-M1)*Mp1-M1;
      Pr	=gyr+k;
      Pi	=gyi+k;
      for(k=M2; k > j; k--){
	pr[k]	=0.;
	pi[k]	=0.;
	Pr[k]	=0.;
	Pi[k]	=0.;
      }
      /* Decomposition loop */
      pr[k]	=am;
      pi[k]	=0.;
      EZvr[k]	=0.;
      vi[k]	=0.;
      k--;
      pr[k]	=0.;
      pi[k]	=0.;
      Pr[-k]	=0.;
      Pi[-k]	=0.;
      k--;
      while(k > jj){
	kk	=k+2;
	A	=(j-k)*(j+k)*bK0;
	bi	=(j+k)*(j+kk);
	br	=bi*Er;
	bi	*=-Ei;
	si	=(j-k)*(j-k+2);
	sr	=si*Er;
	si	*=Ei;
	tr	=A+br*EZvr[kk]-bi*vi[kk];
	ti	=br*vi[kk]+bi*EZvr[kk];
	A	=1./(tr*tr+ti*ti);
	tr	*=-A;
	ti	*=A;
	EZvr[k]	=sr*tr-si*ti;
	vi[k]	=sr*ti+si*tr;
	sr	=br*pr[kk]-bi*pi[kk];
	si	=br*pi[kk]+bi*pr[kk];
	pr[k]	=sr*tr-si*ti;
	pi[k]	=sr*ti+si*tr;
	k--;
	pr[k]	=0.;
	pi[k]	=0.;
	k--;
      }
      pr[k]	=0.;
      pi[k]	=0.;
      Pr[-k]	=0.;
      Pi[-k]	=0.;
      /* Backsubstitution loop */
      kk	=k;
      k++;
      Pr[-k]	=0.;
      Pi[-k]	=0.;
      k++;
      while(k < j){
	pr[k]	+=pr[kk]*EZvr[k]-pi[kk]*vi[k];
	pi[k]	+=pr[kk]*vi[k]+pi[kk]*EZvr[k];
	Pr[-k]	=pr[k];
	Pi[-k]	=-pi[k];
	kk	=k;
	k++;
	Pr[-k]	=0.;
	Pi[-k]	=0.;
	k++;
      }
      Pr[-k]	=pr[k];
      Pi[-k]	=-pi[k];
      for(k=M1; k < jj; k++){
	pr[k]	=0.;
	pi[k]	=0.;
	Pr[k]	=0.;
	Pi[k]	=0.;
      }
      /* Y^j_k=jK_0\gy^j_k+(j+k+2)K^*_2\gy^j_{k+2}+(j-k+2)K_2\gy^j_{k-2} */
      for(k=M1; k <= M2; k++){
	EZvr[k]	=0.;
	vi[k]	=0.;
      }
      for(k=jj+2; k <= j; k +=2){
	EZvr[k]	=j*bK0*pr[k];
	vi[k]	=j*bK0*pi[k];
      }
      k		=jj;
      for(kk=jj+2; kk <= j; kk +=2){
	EZvr[k]	+=(j+kk)*(Er*pr[kk]+Ei*pi[kk]);
	vi[k]	+=(j+kk)*(Er*pi[kk]-Ei*pr[kk]);
	k	+=2;
      }
      k		=j;
      for(kk=j-2; kk >= jj; kk -=2){
	EZvr[k]	+=(j-kk)*(Er*pr[kk]-Ei*pi[kk]);
	vi[k]	+=(j-kk)*(Er*pi[kk]+Ei*pr[kk]);
	k	-=2;
      }
      pr[0]	=gdr*pr[1]-gdi*pi[1];
      pi[0]	=gdr*pi[1]+gdi*pr[1];
      if(M1 < 0){
	pr[0]	+=gdr*pr[-1]+gdi*pi[-1];
	pi[0]	+=gdr*pi[-1]-gdi*pr[-1];
      }
      EZvr[0]	=(j+1)*R0*pr[0]/(x*x);
      vi[0]	=(j+1)*R0*pi[0]/(x*x);
#ifdef gj
#endif
      k		=(jj-M1)*Mp1-M1;
      EZvr	=aYr+k;
      vi	=aYi+k;
      for(k=M1; k <= M2; k++){
	EZvr[k]	=0.;
	vi[k]	=0.;
      }
      for(k=jj; k <= j; k++){
	EZvr[k]	=j*bK0*Pr[k];
	vi[k]	=j*bK0*Pi[k];
      }
      k		=jj;
      for(kk=jj+2; kk <= j; kk +=2){
	EZvr[k]	+=(j+kk)*(Er*Pr[kk]+Ei*Pi[kk]);
	vi[k]	+=(j+kk)*(Er*Pi[kk]-Ei*Pr[kk]);
	k	+=2;
      }
      k		=j;
      for(kk=j-2; kk >= jj; kk -=2){
	EZvr[k]	+=(j-kk)*(Er*Pr[kk]-Ei*Pi[kk]);
	vi[k]	+=(j-kk)*(Er*Pi[kk]+Ei*Pr[kk]);
	k	-=2;
      }
      Pr[0]	=gdr*Pr[1]-gdi*Pi[1];
      Pi[0]	=gdr*Pi[1]+gdi*Pr[1];
      if(M1 < 0){
	Pr[0]	+=gdr*Pr[-1]+gdi*Pi[-1];
	Pi[0]	+=gdr*Pi[-1]-gdi*Pr[-1];
      }
      EZvr[0]	=(j+1)*R0*Pr[0]/(x*x);
      vi[0]	=(j+1)*R0*Pi[0]/(x*x);
#ifdef gj
#endif
      for(mm=0; mm < 2; mm++){
	/* initialization of vatol and rwork */
	kk	=mm ? Mp1*(j-M1) : Mp1*(-j-M1);
	pr	=gyr+kk;
	A	=0.;
	for(k=0; k < Mp1; k++){
	  if(A < fabs(*pr)){
	    A	=fabs(*pr);
	  }
	  pr++;
	}
	A	*=3;
	if(A > reltol){
	  A	=reltol;
	}
	Pr	=vatol+kk;
	Pi	=Pr+MMp1;
	Yr	=Pi+MMp1;
	Yi	=Yr+MMp1;
	for(k=0; k < Mp1; k++){
	  *Pr++	=A;
	  *Pi++	=A;
	  *Yr++	=A;
	  *Yi++	=A;
	}
	A	=1.;
	jj	=j;
	for(i=0; i < 13; i++){
	  pr	=gyr+kk;
	  pi	=gyi+kk;
	  yr	=aYr+kk;
	  yi	=aYi+kk;
	  Pr	=yh+neq*i+kk;
	  Pi	=Pr+MMp1;
	  Yr	=Pi+MMp1;
	  Yi	=Yr+MMp1;
	  for(k=0; k < Mp1; k++){
	    *Pr++	=A*(*pr++);
	    *Pi++	=A*(*pi++);
	    *Yr++	=A*(*yr++);
	    *Yi++	=A*(*yi++);
	  }
	  A		*=(double)jj/(i+1);
	  if(jj){
	    jj--;
	  }
	}
      }
    }
  }
  else{
    /* m=M1 harmonic */
    am	=1.;
    j	=0;
    while(j < M1){
      j++;
      am	*=x;
    }
    if(am < 1e-100){
      am	=1e-100;
    }
    pr	=gyr;
    pi	=gyi;
    yr	=aYr;
    yi	=aYi;
    pr[0]	=am;
    pi[0]	=0.;
    /* Y^j_k=jK_0\gy^j_k+(j+k+2)K^*_2\gy^j_{k+2}+(j-k+2)K_2\gy^j_{k-2} */
    yr[0]	=j*bK0*pr[0];
    yi[0]	=j*bK0*pi[0];
    yr[1]	=0.;
    yi[1]	=0.;
    yr[2]	=j*Er*pr[0];
    yi[2]	=j*Ei*pr[0];
    for(k=3; k < Mp1; k++){
      yr[k]	=0.;
      yi[k]	=0.;
    }	
    A	=3.*am;
    if(A > reltol){
      A	=reltol;
    }
    /* initialization of atol and rwork */
    Pr	=vatol;
    Pi	=Pr+MMp1;
    Yr	=Pi+MMp1;
    Yi	=Yr+MMp1;
    *Pr++	=A;
    *Pi++	=A;
    *Yr++	=A;
    *Yi++	=A;
    for(k=1; k < Mp1; k++){
      pr[k]	=0.;
      pi[k]	=0.;
      *Pr++	=A;
      *Pi++	=A;
      *Yr++	=A;
      *Yi++	=A;
    }	
    A	=1.;
    jj	=j;
    for(i=0; i < 13; i++){
      pr	=gyr;
      pi	=gyi;
      yr	=aYr;
      yi	=aYi;
      Pr	=yh+neq*i;
      Pi	=Pr+MMp1;
      Yr	=Pi+MMp1;
      Yi	=Yr+MMp1;
      for(k=0; k < Mp1; k++){
	*Pr++	=A*(*pr++);
	*Pi++	=A*(*pi++);
	*Yr++	=A*(*yr++);
	*Yi++	=A*(*yi++);
      }
      A		*=(double)jj/(i+1);
      if(jj){
	jj--;
      }
    }
  }
  /* j == |M1| */
  while(j < M2){
    am	*=x;
    if(am < 1e-100){
      am	=1e-100;
    }
    j++;
    k	=(j-M1)*Mp1-M1;
    pr	=gyr+k;
    pi	=gyi+k;
    EZvr	=aYr+k;
    vi	=aYi+k;
    for(k=M2; k > j; k--){
      pr[k]	=0.;
      pi[k]	=0.;
    }
    /* Decomposition loop */
    pr[k]	=am;
    pi[k]	=0.;
    EZvr[k]	=0.;
    vi[k]	=0.;
    k--;
    while(k >= M1){
      pr[k]	=0.;
      pi[k]	=0.;
      k		-=2;
    }
    k		=j-2;
    while(k >= M1){
      kk	=k+2;
      A		=(j-k)*(j+k)*bK0;
      bi	=(j+k)*(j+kk);
      br	=bi*Er;
      bi	*=-Ei;
      si	=(j-k)*(j-k+2);
      sr	=si*Er;
      si	*=Ei;
      tr	=A+br*EZvr[kk]-bi*vi[kk];
      ti	=br*vi[kk]+bi*EZvr[kk];
      A		=1./(tr*tr+ti*ti);
      tr	*=-A;
      ti	*=A;
      EZvr[k]	=sr*tr-si*ti;
      vi[k]	=sr*ti+si*tr;
      sr	=br*pr[kk]-bi*pi[kk];
      si	=br*pi[kk]+bi*pr[kk];
      pr[k]	=sr*tr-si*ti;
      pi[k]	=sr*ti+si*tr;
      k		-=2;
    }
    /* Backsubstitution loop */
    k	+=4;
    kk	=k-2;
    while(k < j){
      pr[k]	+=pr[kk]*EZvr[k]-pi[kk]*vi[k];
      pi[k]	+=pr[kk]*vi[k]+pi[kk]*EZvr[k];
      kk	=k;
      k	+=2;
    }
    /* Y^j_k=jK_0\gy^j_k+(j+k+2)K^*_2\gy^j_{k+2}+(j-k+2)K_2\gy^j_{k-2} */
    for(k=M1; k <= M2; k++){
      EZvr[k]	=0.;
      vi[k]	=0.;
    }
    for(k=j; k >= M1; k	-=2){
      EZvr[k]	=j*bK0*pr[k];
      vi[k]	=j*bK0*pi[k];
    }
    k		=M1;
    for(kk=M1+2; kk <= j; kk++){
      EZvr[k]	+=(j+kk)*(Er*pr[kk]+Ei*pi[kk]);
      vi[k]	+=(j+kk)*(Er*pi[kk]-Ei*pr[kk]);
      k++;
    }
    k		=j;
    for(kk=j-2; kk >= M1; kk--){
      EZvr[k]	+=(j-kk)*(Er*pr[kk]-Ei*pi[kk]);
      vi[k]	+=(j-kk)*(Er*pi[kk]+Ei*pr[kk]);
      k--;
    }
    if(M1 <= 0){
      pr[0]	=gdr*pr[1]-gdi*pi[1];
      pi[0]	=gdr*pi[1]+gdi*pr[1];
      if(M1 < 0){
	pr[0]	+=gdr*pr[-1]+gdi*pi[-1];
	pi[0]	+=gdr*pi[-1]-gdi*pr[-1];
      }
    }
    EZvr[0]	=(j+1)*R0*pr[0]/(x*x);
    vi[0]	=(j+1)*R0*pi[0]/(x*x);

    /* initialization of vatol and rwork */
    kk	=Mp1*(j-M1);
    pr	=gyr+kk;
    A	=0.;
    for(k=0; k < Mp1; k++){
      if(A < fabs(*pr)){
	A	=fabs(*pr);
      }
      pr++;
    }
    A	*=3.;
    if(A > reltol){
      A	=reltol;
    }
    Pr	=vatol+kk;
    Pi	=Pr+MMp1;
    Yr	=Pi+MMp1;
    Yi	=Yr+MMp1;
    for(k=0; k < Mp1; k++){
      *Pr++	=A;
      *Pi++	=A;
      *Yr++	=A;
      *Yi++	=A;
    }
    A	=1.;
    jj	=j;
    for(i=0; i < 13; i++){
      pr	=gyr+kk;
      pi	=gyi+kk;
      yr	=aYr+kk;
      yi	=aYi+kk;
      Pr	=yh+neq*i+kk;
      Pi	=Pr+MMp1;
      Yr	=Pi+MMp1;
      Yi	=Yr+MMp1;
      for(k=0; k < Mp1; k++){
	*Pr++	=A*(*pr++);
	*Pi++	=A*(*pi++);
	*Yr++	=A*(*yr++);
	*Yi++	=A*(*yi++);
      }
      A		*=(double)jj/(i+1);
      if(jj){
	jj--;
      }
    }
  }

  rwork[8-1]	=0.;	/*  */
  rwork[9-1]	=0.;	/*  */
  rwork[10-1]	=0.;	/*  */
  rwork[11-1]	=x;	/* the step size in t last used (successfully) */
  rwork[12-1]	=x;	/* the step size to be attempted on the next step */
  rwork[13-1]	=x;	/* the current value of the independent variable. */
  rwork[14-1]	=2.;	/* a tolerance scale factor, greater than 1.0. */
  rwork[15-1]	=0.;	/*  */
  rwork[16-1]	=0.;	/*  */
  rwork[17-1]	=0.;	/*  */
  rwork[18-1]	=0.;	/*  */
  rwork[19-1]	=0.;	/*  */
  rwork[20-1]	=0.;	/*  */
  iwork[1-1]	=1;	/*# ml, not used */
  iwork[2-1]	=1;	/*# mu, not used */
  iwork[3-1]	=1;	/*#  */
  iwork[4-1]	=1;	/*#  */
  iwork[5-1]	=12;	/* the maximum order to be allowed. */
  iwork[6-1]	=500;	/* maximum number of (internally defined) steps. */
  iwork[7-1]	=10;	/* maximum number of messages printed. */
  iwork[8-1]	=1;	/*#  */
  iwork[9-1]	=1;	/*#  */
  iwork[10-1]	=1;	/*#  */
  iwork[11-1]	=1;	/* the number of steps taken for the problem so far. */
  iwork[12-1]	=0;	/* the number of f evaluations so far. */
  iwork[13-1]	=0;	/* the number of jacobian evaluations. */
  iwork[14-1]	=12;	/* the method order last used (successfully). */
  iwork[15-1]	=12;	/* the order to be attempted on the next step. */
  iwork[16-1]	=1;	/*# the index of the component of largest magnitude */
  iwork[17-1]	=lrw;	/*# the length of rwork actually required. */
  iwork[18-1]	=liw;	/*# the length of iwork actually required. */
  iwork[19-1]	=0;	/*#  */
  iwork[20-1]	=0;	/*#  */

  return(0);
}

#ifdef H
void st2DVrMHD(int *n,double *x,double *Y,double *dY)
{
  double gm,dgm,P,dP,F,dF;
  double sr,si,s;
  double *pr,*pi,*yr,*yi,*dpr,*dpi,*dyr,*dyi;
  double *xr,*xi,*EZvr,*vi,*ur,*ui;

  int i,j,k,kk,m,mm;

  /* Getting metric coefficients{*/
  ESSetSplA(*x);
  sr	=1./(*x);
  splRA(gaar,NULL,ESg11c,ESg11c2);
  gaar[0]	*=sr;
  gaai[0]	=0.;
  splRA(gaqr,NULL,ESg12c,ESg12c2);
  gaqr[0]	=-gaqr[0];
  gaqi[0]	=0.;

  splRA(gqqr,NULL,ESg22c,ESg22c2);
  gqqr[0]	*=*x;
  gqqi[0]	=0.;
  splRA(gzzr,NULL,ESGzzc,ESGzzc2);
  gzzr[0]	*=sr;
  gqqi[0]	=0.;
  splRA(svr,dsvr,ESsVc,ESsVc2);
  svi[0]	=0.;
  dsvi[0]	=0.;
  m	=0;
  for(k=1; k < ES2Mp1; k++){
    m	+=ESNa1;
    splRA(gaar+k,NULL,ESg11c+m,ESg11c2+m);
    splRA(gaai+k,NULL,ESg11s+m,ESg11s2+m);
    gaar[k]	*=sr;
    gaai[k]	*=sr;
    splRA(gaqr+k,NULL,ESg12c+m,ESg12c2+m);
    splRA(gaqi+k,NULL,ESg12s+m,ESg12s2+m);
    gaqr[k]	=-gaqr[k];
    gaqi[k]	=-gaqi[k];
    splRA(gqqr+k,NULL,ESg22c+m,ESg22c2+m);
    splRA(gqqi+k,NULL,ESg22s+m,ESg22s2+m);
    gqqr[k]	*=*x;
    gqqi[k]	*=*x;
    splRA(gzzr+k,NULL,ESGzzc+m,ESGzzc2+m);
    splRA(gzzi+k,NULL,ESGzzs+m,ESGzzs2+m);
    gzzr[k]	*=sr;
    gzzi[k]	*=sr;

    splRA(svr+k,NULL,ESVc+m,ESVc2+m);
    splRA(svi+k,NULL,ESsVs+m,ESsVs2+m);
  }
  /*}*/

  /* Getting 1D profiles{*/
  ESSetSplDPr(*x);
  ESSetSplDCr(*x);
  if(ESEqSolvInCr == 6 || ESEqSolvInCr == 7){
    splRDCr(&gm,&dgm,ESEqSolvInCr);
  }
  else{
    splRA(&gm,&dgm,ESgm,ESgm2a);
  }
  splRA(&F,&dF,ESaT,ESaT2a);
  if(ESEqSolvInPr == 0){
    splRDPr(&P,&dP,0);
    P	*=rR0;
    dP	*=rR0;
  }
  else{
    if(ESEqSolvInPr == 1){
      splRDPr(&P,&dP,1);
    }
    else{
      splRA(&P,&dP,ESPs,ESPs2a);
    }
  }
  /*}*/
  /* Elimination of $\vgl$ {*/
#ifdef H
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1+k;
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1; kk++){
      if(i < ES2Mp1){
	ar[j]	=m*gzzr[i]*mm+N*N*gqqr[i];
	ai[j]	=m*gzzi[i]*mm+N*N*gqqi[i];
      }
      else{
	ar[j]	=0.;
	ai[j]	=0.;
      }
      kkkk	=kk;
      for(kkk=0; kkkk <ES2Mp1; kkk++){
	ar[j]	+=m*N*(gzzr[kkkk]*aZqr[kkk]+gzzi[kkkk]*aZqi[kkk])
	  +mm*N*(aZqr[kkkk]*gzzr[kkk]-aZqi[kkkk]*gzzi[kkk]);
	ai[j]	+=m*N*(gzzr[kkkk]*aZqi[kkk]-gzzi[kkkk]*aZqr[kkk])
	  +mm*N*(aZqr[kkkk]*gzzi[kkk]+aZqi[kkkk]*gzzr[kkk]);
	kkkk++;
      }
      kkkk	=-kk+1;
      for(kkk=1; kkkk < 1 && kkk < ES2Mp1; kkk++){
	ar[j]	+=m*N*(gzzr[kkkk]*aZqr[kkk]-gzzi[kkkk]*aZqi[kkk])
	  +mm*N*(aZqr[kkkk]*gzzr[kkk]+aZqi[kkkk]*gzzi[kkk]);
	ai[j]	+=m*N*(-gzzr[kkkk]*aZqi[kkk]-gzzi[kkkk]*aZqr[kkk])
	  +mm*N*(-aZqr[kkkk]*gzzi[kkk]+aZqi[kkkk]*gzzr[kkk]);
	kkkk++;
      }
      while(kkk < ES2Mp1){
	ar[j]	+=m*N*(gzzr[kkkk]*aZqr[kkk]+gzzi[kkkk]*aZqi[kkk])
	  +mm*N*(aZqr[kkkk]*gzzr[kkk]-aZqi[kkkk]*gzzi[kkk]);
	ai[j]	+=m*N*(-gzzr[kkkk]*aZqi[kkk]+gzzi[kkkk]*aZqr[kkk])
	  +mm*N*(-aZqr[kkkk]*gzzi[kkk]-aZqi[kkkk]*gzzr[kkk]);
	kkk++;
	kkkk++;
      }
      i++;
      j++;
      mm++;
    }
    m++;
  }
  if(CholDeComp(ar,ai,Mp1) < 0){
#endif

  /*}*/

  /* Cholesky decomposition of $\maK${*/
  for(m=0; m < Mp1; m++){
    k	=m*Mp1+m;
    i	=0;
    for(mm=m; mm < Mp1; mm++){
      if(i < ES2Mp1){
	ar[k]	=gqqr[i];
	ai[k]	=gqqi[i];
      }
      else{
	ar[k]	=0.;
	ai[k]	=0.;
      }
      i++;
      k++;
    }
  }
  if(CholDeComp(ar,ai,Mp1) < 0){
#ifdef XWIN
    printf("??? x=%10.3e\n",*x);
    for(m=0; m < Mp1; m++){
      printf("??? Kr[%2d][%2d]=%10.3e\n",m,m,ar[m*Mp1+m]);
    }
#endif
    FlStChol	=1;
    return;
  }
  /*}*/

  pr	=Y;
  pi	=pr+MMp1;
  yr	=pi+MMp1;
  yi	=yr+MMp1;
  dpr	=dY;
  dpi	=dpr+MMp1;
  dyr	=dpi+MMp1;
  dyi	=dyr+MMp1;
  for(j=0; j < Mp1; j++){
    /* Calculation of $dgy${*/
    xr	=dpr;
    xi	=dpi;
    for(k=0; k < Mp1; k++){
      *xr	=yr[k];
      *xi	=yi[k];
      i	=0;
      mm	=M1+k;
      EZvr	=pr+k;
      vi	=pi+k;
      for(kk=k; kk < Mp1 && i < ES2Mp1; kk++){
	*xi	+=mm*(gaqr[i]*(*EZvr)-gaqi[i]*(*vi));
	*xr	-=mm*(gaqr[i]*(*vi)+gaqi[i]*(*EZvr));
	i++;
	mm++;
	EZvr++;
	vi++;
      }
      i		=1;
      kk	=k;
      mm	=M1+k;
      EZvr	=pr+k;
      vi	=pi+k;
      while(kk > 0 && i < ES2Mp1){
	kk--;
	EZvr--;
	vi--;
	mm--;
	*xi	+=mm*(gaqr[i]*(*EZvr)+gaqi[i]*(*vi));
	*xr	-=mm*(gaqr[i]*(*vi)-gaqi[i]*(*EZvr));
	i++;
      }
      xr++;
      xi++;
    }
    CholBksbComp(ar,ai,Mp1,dpr,dpi);
    /*}*/
    /* Calculation of $daY${*/
    xr	=dyr;
    xi	=dyi;
    for(k=0; k < Mp1; k++){
      *xr	=0.;
      *xi	=0.;
      m	=M1+k;
      i	=0;
      mm=m;
      EZvr	=pr+k;
      vi	=pi+k;
      ur	=dpr+k;
      ui	=dpi+k;
      for(kk=k; kk < Mp1 && i < ES2Mp1; kk++){
	*xr	+=mm*(gaar[i]*(*EZvr)-gaai[i]*(*vi));
	*xi	+=mm*(gaar[i]*(*vi)+gaai[i]*(*EZvr));
	*xi	+=gaqr[i]*(*ur)-gaqi[i]*(*ui);
	*xr	-=gaqr[i]*(*ui)+gaqi[i]*(*ur);
	EZvr++;
	vi++;
	ur++;
	ui++;
	i++;
	mm++;
      }
      i		=1;
      kk	=k;
      mm	=m;
      EZvr	=pr+k;
      vi	=pi+k;
      ur	=dpr+k;
      ui	=dpi+k;
      while(kk > 0 && i < ES2Mp1){
	kk--;
	EZvr--;
	vi--;
	ur--;
	ui--;
	mm--;
	*xr	+=mm*(gaar[i]*(*EZvr)+gaai[i]*(*vi));
	*xi	+=mm*(gaar[i]*(*vi)-gaai[i]*(*EZvr));
	*xi	+=gaqr[i]*(*ur)+gaqi[i]*(*ui);
	*xr	-=gaqr[i]*(*ui)-gaqi[i]*(*ur);
	i++;
      }
      *xr	*=m;
      *xi	*=m;
      xr++;
      xi++;
    }
    /*}*/
    pr	+=Mp1;
    pi	+=Mp1;
    yr	+=Mp1;
    yi	+=Mp1;
    dpr	+=Mp1;
    dpi	+=Mp1;
    dyr	+=Mp1;
    dyi	+=Mp1;
  }
}
#endif

/* **************************** */
void st2DVrMHD(int *n,double *x,double *Y,double *dY)
{
  double gm,dgm,P,dP,F,dF;
  double sr,si,s;
  double *pr,*pi,*yr,*yi,*dpr,*dpi,*dyr,*dyi;
  double *xr,*xi,*EZvr,*vi,*ur,*ui;

  int i,j,k,kk,m,mm;

  /* Getting metric coefficients{*/
  ESSetSplA(*x);
  sr	=1./(*x);
  splRA(gaar,NULL,ESGaac,ESGaac2);
  gaar[0]	*=sr;
  gaai[0]	=0.;

  splRA(gaqr,NULL,ESGaqc,ESGaqc2);
  gaqr[0]	=-gaqr[0];
  gaqi[0]	=0.;

  splRA(gqqr,NULL,ESGqqc,ESGqqc2);
  gqqr[0]	*=*x;
  gqqi[0]	=0.;

  splRA(gazr,NULL,ESGazc,ESGazc2);
  gazr[0]	*=-sr;
  gazi[0]	=0.;
  splRA(gqzr,NULL,ESGqzc,ESGqzc2);
  gqzi[0]	=0.;
  splRA(gzzr,NULL,ESGzzc,ESGzzc2);
  gzzr[0]	*=sr;
  gzzi[0]	=0.;

  splRA(svr,dsvr,ESsVc,ESsVc2);
  svi[0]	=0.;
  dsvi[0]	=0.;
  m	=0;
  for(k=1; k < ES2Mp1; k++){
    m	+=ESNa1;
    splRA(gaar+k,NULL,ESGaac+m,ESGaac2+m);
    splRA(gaai+k,NULL,ESGaas+m,ESGaas2+m);
    gaar[k]	*=sr;
    gaai[k]	*=sr;
    splRA(gaqr+k,NULL,ESGaqc+m,ESGaqc2+m);
    splRA(gaqi+k,NULL,ESGaqs+m,ESGaqs2+m);
    gaqr[k]	=-gaqr[k];
    gaqi[k]	=-gaqi[k];
    splRA(gqqr+k,NULL,ESGqqc+m,ESGqqc2+m);
    splRA(gqqi+k,NULL,ESGqqs+m,ESGqqs2+m);
    gqqr[k]	*=*x;
    gqqi[k]	*=*x;

    splRA(gazr+k,NULL,ESGazc+m,ESGazc2+m);
    splRA(gazi+k,NULL,ESGazs+m,ESGazs2+m);
    gazr[k]	*=-sr;
    gazi[k]	*=-sr;
    splRA(gqzr+k,NULL,ESGqzc+m,ESGqzc2+m);
    splRA(gqzi+k,NULL,ESGqzs+m,ESGqzs2+m);
    splRA(gzzr+k,NULL,ESGzzc+m,ESGzzc2+m);
    splRA(gzzi+k,NULL,ESGzzs+m,ESGzzs2+m);
    gzzr[k]	*=sr;
    gzzi[k]	*=sr;

    splRA(svr+k,dsvr+k,ESsVc+m,ESsVc2+m);
    splRA(svi+k,dsvi+k,ESsVs+m,ESsVs2+m);
  }
  /*}*/

  /* Getting 1D profiles{*/
  splRA(&P,&dP,ESFF,ESFF2a);
  dP	*=-EZcr2/P;
  splRA(&F,&dF,ESaT,ESaT2a);
  P	=1./sqrt(P);
  F	*=P;
  dF	=(dF+F*dP)*P;
  if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
    splRA(&gm,&dgm,ESgm,ESgm2a);
  }
  else{
    ESSetSplDCr(*x);
    splRDCr(&gm,&dgm,ESEqSolvInCr);
  }
  m	=M1;
/**************/
  for(k=0; k < Mp1; k++){
    Tr[k]	=m == Mres ? 1./(gm-gmRes[kRes]) : m/(m*gm-N);
    m++;
  }
  if(ESEqSolvInPr > 2){
    splRA(&P,&dP,ESPs,ESPs2a);
  }
  else{
    ESSetSplDPr(*x);
    splRDPr(&P,&dP,ESEqSolvInPr);
    if(ESEqSolvInPr == 0){
      P	*=rR0;
      dP*=rR0;
    }
  }
/**************/
#ifdef H
  P	=0.;
  dP	=0.;
  F	=0.;
  dF	=0.;
#endif

  /*}*/

  /* Elimination of $\vgl$ {*/
  /* Calculation of $\maA=\msa\msa^\dagger$ {*/
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1+k;
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1; kk++){
      if(i < ES2Mp1){
	Ar[j]	=m*gzzr[i]*mm+N*((m+mm)*gqzr[i]+N*gqqr[i]);
	Ai[j]	=m*gzzi[i]*mm+N*((m+mm)*gqzi[i]+N*gqqi[i]);
      }
      else{
	Ar[j]	=0.;
	Ai[j]	=0.;
      }
      i++;
      j++;
      mm++;
    }
    m++;
  }
  if(CholDeComp(Ar,Ai,Mp1) < 0){
#ifdef XWIN
    printf("??? Matrix A inversion failed at x=%10.3e\n",*x);
    for(m=0; m < Mp1; m++){
      printf("??? Kr[%2d][%2d]=%10.3e\n",m,m,ar[m*Mp1+m]);
    }
#endif
    FlStChol	=1;
    return;
  }
  /*}*/
  /* Calculation of $\msa^{-1}\maB$ {*/
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1;
    pr	=Br+j;
    pi	=Bi+j;
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1; kk++){
      if(i < ES2Mp1){
	if(m != 0){
	  pr[kk]=-N*gqqr[i]-mm*gqzr[i];
	  pi[kk]= N*gqqi[i]+mm*gqzi[i];
	}
	else{
	  pr[kk]= N*gqzr[i]+mm*gzzr[i];
	  pi[kk]=-N*gqzi[i]-mm*gzzi[i];
	}
      }
      else{
	pr[kk]	=0.;
	pi[kk]	=0.;
      }
      i++;
      mm++;
    }
    i	=0;
    kk	=k;
    mm	=kk+M1;
    while(kk > 0){
      kk--;
      mm--;
      i++;
      if(i < ES2Mp1){
	if(m != 0){
	  pr[kk]	=-N*gqqr[i]-mm*gqzr[i];
	  pi[kk]	=-N*gqqi[i]-mm*gqzi[i];
	}
	else{
	  pr[kk]	=N*gqzr[i]+mm*gzzr[i];
	  pi[kk]	=N*gqzi[i]+mm*gzzi[i];
	}
      }
      else{
	pr[kk]	=0.;
	pi[kk]	=0.;
      }
    }
    kk		=0;
    for(mm=0; mm < Mp1; mm++){
      sr	=pr[mm];
      si	=pi[mm];
      j		=mm;
      i		=kk+mm;
      while(j > 0){
	j--;
	i--;
	sr	-=Ar[i]*pr[j]-Ai[i]*pi[j];
	si	-=Ar[i]*pi[j]+Ai[i]*pr[j];
      }
      s		=Ar[kk+mm];
      pr[mm]	=s*sr;
      pi[mm]	=s*si;
      kk	+=Mp1;
    }
    m++;
  }
  /*}*/
  /* Calculation of $\maF=\maD-\maB^\dagger\msa^{-1\dagger}\msa^{-1}\maB$ {*/
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1;
    yr	=Fr+j;
    yi	=Fi+j;
    ur	=Br+j;
    ui	=Bi+j;
    EZvr	=Br+j;
    vi	=Bi+j;
    for(kk=k; kk < Mp1; kk++){
      sr	=0.;
      si	=0.;
      for(mm=0; mm < Mp1; mm++){
	sr	+=ur[mm]*(*EZvr)+ui[mm]*(*vi);
	si	+=ur[mm]*(*vi)-ui[mm]*(*EZvr);
	EZvr++;
	vi++;
      }	
      yr[kk]	=-sr;
      yi[kk]	=-si;
    }
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1 && i < ES2Mp1; kk++){
      if(m != 0){
	if(mm != 0){
	  yr[kk]	+=gqqr[i];
	  yi[kk]	+=gqqi[i];
	}
	else{
	  yr[kk]	-=gqzr[i];
	  yi[kk]	-=gqzi[i];
	}
      }
      else{
	if(mm != 0){
	  yr[kk]	-=gqzr[i];
	  yi[kk]	-=gqzi[i];
	}
	else{
	  yr[kk]	+=gzzr[i];
	  yi[kk]	+=gzzi[i];
	}
      }
      i++;
      mm++;
    }
    m++;
  }
  if(CholDeComp(Fr,Fi,Mp1) < 0){
#ifdef XWIN
    printf("??? Matrix F inversion failed at x=%10.3e\n",*x);
    for(m=0; m < Mp1; m++){
      printf("??? Fr[%2d][%2d]=%10.3e\n",m,m,Fr[m*Mp1+m]);
    }
#endif
    FlStChol	=1;
    return;
  }
  /*}*/
  /* Calculation of $\msa^{-1}\maC${*/
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1;
    yr	=Cr+j;
    yi	=Ci+j;
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1; kk++){
      if(i < ES2Mp1){
	if(m != 0){
	  yr[kk]=m*(N*gaqi[i]+mm*gazi[i])+N*P*Tr[k]*svr[i];
	  yi[kk]=m*(N*gaqr[i]+mm*gazr[i])-N*P*Tr[k]*svi[i];
	  if(mm == m){
	    yr[kk]	-=F*mm;
	  }
	}
	else{
	  yr[kk]=(N*gaqi[i]+mm*gazi[i]-P*svr[i])*N;
	  yi[kk]=(N*gaqr[i]+mm*gazr[i]+P*svi[i])*N;
	  if(mm == 0){
	    yr[kk]	-=N*F;
	  }
	}
      }
      else{
	yr[kk]	=0.;
	yi[kk]	=0.;
      }
      i++;
      mm++;
    }
    i	=0;
    kk	=k;
    mm	=kk+M1;
    while(kk > 0){
      kk--;
      mm--;
      i++;
      if(i < ES2Mp1){
	if(m != 0){
	  yr[kk]=m*(-N*gaqi[i]-mm*gazi[i])+N*P*Tr[k]*svr[i];
	  yi[kk]=m*( N*gaqr[i]+mm*gazr[i])+N*P*Tr[k]*svi[i];
	}
	else{
	  yr[kk]=(-N*gaqi[i]-mm*gazi[i]-P*svr[i])*N;
	  yi[kk]=( N*gaqr[i]+mm*gazr[i]-P*svi[i])*N;
	}
      }
      else{
	yr[kk]	=0.;
	yi[kk]	=0.;
      }
    }
    kk		=0;
    for(mm=0; mm < Mp1; mm++){
      sr	=yr[mm];
      si	=yi[mm];
      j		=mm;
      i		=kk+mm;
      while(j > 0){
	j--;
	i--;
	sr	-=Ar[i]*yr[j]-Ai[i]*yi[j];
	si	-=Ar[i]*yi[j]+Ai[i]*yr[j];
      }
      s		=Ar[kk+mm];
      yr[mm]	=s*sr;
      yi[mm]	=s*si;
      kk	+=Mp1;
    }
    m++;
  }
  /*}*/
  /*}*/

  /* Calculation of $\maK=\maE-\maB^\dagger(\msa^\dagger)^{-1}\msa^{-1}\maC${*/
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1;
    yr	=Kr+j;
    yi	=Ki+j;
    ur	=Br+j;
    ui	=Bi+j;
    EZvr	=Cr;
    vi	=Ci;
    for(kk=0; kk < Mp1; kk++){
      sr	=0.;
      si	=0.;
      for(mm=0; mm < Mp1; mm++){
	sr	+=ur[mm]*(*EZvr)+ui[mm]*(*vi);
	si	+=ur[mm]*(*vi)-ui[mm]*(*EZvr);
	EZvr++;
	vi++;
      }
      yr[kk]	=-sr;
      yi[kk]	=-si;
    }
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1 && i < ES2Mp1; kk++){
      if(m != 0){
	if(mm != 0){
	  yr[kk]	+=gaqi[i]*mm;
	  yi[kk]	-=gaqr[i]*mm;
	  if(i){
	    yr[kk]	-=P*svr[i]*Tr[kk];
	    yi[kk]	-=P*svi[i]*Tr[kk];
	  }
	}
	else{
	  yr[kk]	+=N*gaqi[i];
	  yi[kk]	-=N*gaqr[i];
	  if(i){
	    yr[kk]	+=P*svr[i];
	    yi[kk]	+=P*svi[i];
	  }
	}
      }
      else{
	if(mm != 0){
	  yr[kk]	-=gazi[i]*mm;
	  yi[kk]	+=gazr[i]*mm;
	}
	else{
	  yr[kk]	-=N*gazi[i];
	  yi[kk]	+=N*gazr[i];
	}
      }
      i++;
      mm++;
    }
    i	=1;
    kk	=k;
    mm	=kk+M1;
    while(kk > 0 && i < ES2Mp1){
      kk--;
      mm--;
      if(m != 0){
	if(mm != 0){
	  yr[kk]	-=gaqi[i]*mm;
	  yi[kk]	-=gaqr[i]*mm;
	  if(i){
	    yr[kk]	-=P*svr[i]*Tr[kk];
	    yi[kk]	+=P*svi[i]*Tr[kk];
	  }
	}
	else{
	  yr[kk]	-=N*gaqi[i];
	  yi[kk]	-=N*gaqr[i];
	  if(i){
	    yr[kk]	+=P*svr[i];
	    yi[kk]	-=P*svi[i];
	  }
	}
      }
      else{
	if(mm != 0){
	  yr[kk]	+=gazi[i]*mm;
	  yi[kk]	+=gazr[i]*mm;
	}
	else{
	  yr[kk]	+=N*gazi[i];
	  yi[kk]	+=N*gazr[i];
	}
      }
      i++;
    }
    m++;
  }
  /*}*/

  /* Calculation of $\maG=\maH-\maC^\dagger\msa^{-1\dagger}\msa^{-1}\maC$ {*/
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1;
    yr	=Gr+j;
    yi	=Gi+j;
    ur	=Cr+j;
    ui	=Ci+j;
    EZvr	=Cr+j;
    vi	=Ci+j;
    for(kk=k; kk < Mp1; kk++){
      sr	=0.;
      si	=0.;
      for(mm=0; mm < Mp1; mm++){
	sr	+=ur[mm]*(*EZvr)+ui[mm]*(*vi);
	si	+=ur[mm]*(*vi)-ui[mm]*(*EZvr);
	EZvr++;
	vi++;
      }
      yr[kk]	=-sr;
      yi[kk]	=-si;
    }
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1 && i < ES2Mp1; kk++){
      if(m != 0){
	if(mm != 0){
	  yr[kk]	+=m*gaar[i]*mm;
	  yi[kk]	+=m*gaai[i]*mm;
	  if(i){
	    yr[kk]	+=P*Tr[k]*(dgm*svr[i]-gm*dsvr[i])*Tr[kk];
	    yi[kk]	+=P*Tr[k]*(dgm*svi[i]-gm*dsvi[i])*Tr[kk];
	  }
	  if(i == 0){
	    yr[kk]	+=(dP*svr[0]+dF-N*P*dsvr[0]*Tr[k]/m)*Tr[k];
	  }
	}
	else{
	  yr[kk]	+=N*m*gaar[i];
	  yi[kk]	+=N*m*gaai[i];
	  if(i){
	    yr[kk]	-=P*Tr[k]*(dgm*svr[i]-gm*dsvr[i]);
	    yi[kk]	-=P*Tr[k]*(dgm*svi[i]-gm*dsvi[i]);
	  }
	}
      }
      else{
	if(mm != 0){
	  yr[kk]	+=N*gaar[i]*mm;
	  yi[kk]	+=N*gaai[i]*mm;
	  if(i){
	    yr[kk]	-=P*(dgm*svr[i]-gm*dsvr[i])*Tr[kk];
	    yi[kk]	-=P*(dgm*svi[i]-gm*dsvi[i])*Tr[kk];
	  }
	}
	else{
	  yr[kk]	+=N*N*gaar[i]+P*(dgm*svr[0]-gm*dsvr[0])-dF*gm;
	}
      }
      if(i){
	j	=Mp1*kk+k;
	Gr[j]	= yr[kk];
	Gi[j]	=-yi[kk];
      }
      i++;
      mm++;
    }
    m++;
  }
  /*}*/

  pr	=Y;
  pi	=pr+MMp1;
  yr	=pi+MMp1;
  yi	=yr+MMp1;
  dpr	=dY;
  dpi	=dpr+MMp1;
  dyr	=dpi+MMp1;
  dyi	=dyr+MMp1;
  for(j=0; j < Mp1; j++){
    /* Calculation of $\vgy'=\maF^{-1}(\vaY-\maK\vgy)$ {*/
    xr	=dpr;
    xi	=dpi;
    for(k=0; k < Mp1; k++){
      sr	=yr[k];
      si	=yi[k];
      ur	=Kr+k*Mp1;
      ui	=Ki+k*Mp1;
      EZvr	=pr;
      vi	=pi;
      for(kk=0; kk < Mp1; kk++){
	sr	-=(*ur)*(*EZvr)-(*ui)*(*vi);
	si	-=(*ur)*(*vi)+(*ui)*(*EZvr);
	ur++;
	ui++;
	EZvr++;
	vi++;
      }
      *xr++	=sr;
      *xi++	=si;
    }
    CholBksbComp(Fr,Fi,Mp1,dpr,dpi);
    /*}*/

    /* Calculation of $\vaY'=\maK^\dagger\vgy'+\maG\vgy {*/
    for(k=0; k < Mp1; k++){
      sr	=0.;
      si	=0.;
      xr	=dpr;
      xi	=dpi;
      ur	=Gr+k*Mp1;
      ui	=Gi+k*Mp1;
      EZvr	=pr;
      vi	=pi;
      mm	=k;
      for(kk=0; kk < Mp1; kk++){
	sr	+=(*ur)*(*EZvr)-(*ui)*(*vi)+Kr[mm]*(*xr)+Ki[mm]*(*xi);
	si	+=(*ur)*(*vi)+(*ui)*(*EZvr)+Kr[mm]*(*xi)-Ki[mm]*(*xr);
	ur++;
	ui++;
	EZvr++;
	vi++;
	xr++;
	xi++;
	mm	+=Mp1;
      }
      *dyr++	=sr;
      *dyi++	=si;
    }
    /*}*/
    pr	+=Mp1;
    pi	+=Mp1;
    yr	+=Mp1;
    yi	+=Mp1;
    dpr	+=Mp1;
    dpi	+=Mp1;
  }
}

int ESResAsymptotics()
{
  double gm,dgm,P,x;
  double sr,si,s,a11,a12,a21,a22;
  double *pr,*pi,*yr,*yi;
  double *xr,*xi,*EZvr,*vi,*ur,*ui;

  int i,j,k,kk,m,mm;

  /* Getting metric coefficients{*/
  x	=aRes[kRes];

  ESSetSplA(x);
  sr	=1./(x);
  splRA(gaar,NULL,ESGaac,ESGaac2);
  gaar[0]	*=sr;
  gaai[0]	=0.;

  splRA(gaqr,NULL,ESGaqc,ESGaqc2);
  gaqr[0]	=-gaqr[0];
  gaqi[0]	=0.;

  splRA(gqqr,NULL,ESGqqc,ESGqqc2);
  gqqr[0]	*=x;
  gqqi[0]	=0.;

  splRA(gazr,NULL,ESGazc,ESGazc2);
  gazr[0]	*=-sr;
  gazi[0]	=0.;
  splRA(gqzr,NULL,ESGqzc,ESGqzc2);
  gqzi[0]	=0.;
  splRA(gzzr,NULL,ESGzzc,ESGzzc2);
  gzzr[0]	*=sr;
  gzzi[0]	=0.;

  splRA(svr,dsvr,ESsVc,ESsVc2);
  svi[0]	=0.;
  dsvi[0]	=0.;
  m	=0;
  for(k=1; k < ES2Mp1; k++){
    m	+=ESNa1;
    splRA(gaar+k,NULL,ESGaac+m,ESGaac2+m);
    splRA(gaai+k,NULL,ESGaas+m,ESGaas2+m);
    gaar[k]	*=sr;
    gaai[k]	*=sr;
    splRA(gaqr+k,NULL,ESGaqc+m,ESGaqc2+m);
    splRA(gaqi+k,NULL,ESGaqs+m,ESGaqs2+m);
    gaqr[k]	=-gaqr[k];
    gaqi[k]	=-gaqi[k];
    splRA(gqqr+k,NULL,ESGqqc+m,ESGqqc2+m);
    splRA(gqqi+k,NULL,ESGqqs+m,ESGqqs2+m);
    gqqr[k]	*=x;
    gqqi[k]	*=x;

    splRA(gazr+k,NULL,ESGazc+m,ESGazc2+m);
    splRA(gazi+k,NULL,ESGazs+m,ESGazs2+m);
    gazr[k]	*=-sr;
    gazi[k]	*=-sr;
    splRA(gqzr+k,NULL,ESGqzc+m,ESGqzc2+m);
    splRA(gqzi+k,NULL,ESGqzs+m,ESGqzs2+m);
    splRA(gzzr+k,NULL,ESGzzc+m,ESGzzc2+m);
    splRA(gzzi+k,NULL,ESGzzs+m,ESGzzs2+m);
    gzzr[k]	*=sr;
    gzzi[k]	*=sr;

    splRA(svr+k,dsvr+k,ESsVc+m,ESsVc2+m);
    splRA(svi+k,dsvi+k,ESsVs+m,ESsVs2+m);
  }
  /*}*/

  /* Getting 1D profiles{*/
  gm	=gmRes[kRes];
  dgm	=dgmRes[kRes];
  for(k=0; k < Mp1; k++){
    Tr[k]	=0.;
  }
  Tr[Mres-M1]	=1.;

  if(ESEqSolvInPr > 2){
    splRA(&P,NULL,ESPs,ESPs2a);
  }
  else{
    ESSetSplDPr(x);
    splRDPr(&P,NULL,ESEqSolvInPr);
    if(ESEqSolvInPr == 0){
      P	*=rR0;
    }
  }

  /*}*/

  /* Calculation of $\maA=\msa\msa^\dagger$ {*/
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1+k;
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1; kk++){
      if(i < ES2Mp1){
	Ar[j]	=m*gzzr[i]*mm+N*((m+mm)*gqzr[i]+N*gqqr[i]);
	Ai[j]	=m*gzzi[i]*mm+N*((m+mm)*gqzi[i]+N*gqqi[i]);
      }
      else{
	Ar[j]	=0.;
	Ai[j]	=0.;
      }
      i++;
      j++;
      mm++;
    }
    m++;
  }
  if(CholDeComp(Ar,Ai,Mp1) < 0){
#ifdef XWIN
    printf("??? Matrix A inversion failed at x=%10.3e\n",x);
    for(m=0; m < Mp1; m++){
      printf("??? Kr[%2d][%2d]=%10.3e\n",m,m,ar[m*Mp1+m]);
    }
#endif
    FlStChol	=1;
    return(0);
  }
  /*}*/
  /* Calculation of $\msa^{-1}\maB$ {*/
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1;
    pr	=Br+j;
    pi	=Bi+j;
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1; kk++){
      if(i < ES2Mp1){
	if(m != 0){
	  pr[kk]=-N*gqqr[i]-mm*gqzr[i];
	  pi[kk]= N*gqqi[i]+mm*gqzi[i];
	}
	else{
	  pr[kk]= N*gqzr[i]+mm*gzzr[i];
	  pi[kk]=-N*gqzi[i]-mm*gzzi[i];
	}
      }
      else{
	pr[kk]	=0.;
	pi[kk]	=0.;
      }
      i++;
      mm++;
    }
    i	=0;
    kk	=k;
    mm	=kk+M1;
    while(kk > 0){
      kk--;
      mm--;
      i++;
      if(i < ES2Mp1){
	if(m != 0){
	  pr[kk]	=-N*gqqr[i]-mm*gqzr[i];
	  pi[kk]	=-N*gqqi[i]-mm*gqzi[i];
	}
	else{
	  pr[kk]	=N*gqzr[i]+mm*gzzr[i];
	  pi[kk]	=N*gqzi[i]+mm*gzzi[i];
	}
      }
      else{
	pr[kk]	=0.;
	pi[kk]	=0.;
      }
    }
    kk		=0;
    for(mm=0; mm < Mp1; mm++){
      sr	=pr[mm];
      si	=pi[mm];
      j		=mm;
      i		=kk+mm;
      while(j > 0){
	j--;
	i--;
	sr	-=Ar[i]*pr[j]-Ai[i]*pi[j];
	si	-=Ar[i]*pi[j]+Ai[i]*pr[j];
      }
      s		=Ar[kk+mm];
      pr[mm]	=s*sr;
      pi[mm]	=s*si;
      kk	+=Mp1;
    }
    m++;
  }
  /*}*/
  /* Calculation of $\maF=\maD-\maB^\dagger\msa^{-1\dagger}\msa^{-1}\maB$ {*/
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1;
    yr	=Fr+j;
    yi	=Fi+j;
    ur	=Br+j;
    ui	=Bi+j;
    EZvr	=Br+j;
    vi	=Bi+j;
    for(kk=k; kk < Mp1; kk++){
      sr	=0.;
      si	=0.;
      for(mm=0; mm < Mp1; mm++){
	sr	+=ur[mm]*(*EZvr)+ui[mm]*(*vi);
	si	+=ur[mm]*(*vi)-ui[mm]*(*EZvr);
	EZvr++;
	vi++;
      }	
      yr[kk]	=-sr;
      yi[kk]	=-si;
    }
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1 && i < ES2Mp1; kk++){
      if(m != 0){
	if(mm != 0){
	  yr[kk]	+=gqqr[i];
	  yi[kk]	+=gqqi[i];
	}
	else{
	  yr[kk]	-=gqzr[i];
	  yi[kk]	-=gqzi[i];
	}
      }
      else{
	if(mm != 0){
	  yr[kk]	-=gqzr[i];
	  yi[kk]	-=gqzi[i];
	}
	else{
	  yr[kk]	+=gzzr[i];
	  yi[kk]	+=gzzi[i];
	}
      }
      i++;
      mm++;
    }
    m++;
  }
  if(CholDeComp(Fr,Fi,Mp1) < 0){
#ifdef XWIN
    printf("??? Matrix F inversion failed at x=%10.3e\n",x);
    for(m=0; m < Mp1; m++){
      printf("??? Fr[%2d][%2d]=%10.3e\n",m,m,Fr[m*Mp1+m]);
    }
#endif
    FlStChol	=1;
    return(0);
  }
  /*}*/
  /* Calculation of $\msa^{-1}\maC${*/
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1;
    yr	=Cr+j;
    yi	=Ci+j;
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1; kk++){
      if(i < ES2Mp1){
	if(m != 0){
	  yr[kk]	= N*P*Tr[k]*svr[i];
	  yi[kk]	=-N*P*Tr[k]*svi[i];
	}
	else{
	  yr[kk]	=0.;
	  yi[kk]	=0.;
	}
      }
      else{
	yr[kk]	=0.;
	yi[kk]	=0.;
      }
      i++;
      mm++;
    }
    i	=0;
    kk	=k;
    mm	=kk+M1;
    while(kk > 0){
      kk--;
      mm--;
      i++;
      if(i < ES2Mp1){
	if(m != 0){
	  yr[kk]	=N*P*Tr[k]*svr[i];
	  yi[kk]	=N*P*Tr[k]*svi[i];
	}
	else{
	  yr[kk]	=0.;
	  yi[kk]	=0.;
	}
      }
      else{
	yr[kk]	=0.;
	yi[kk]	=0.;
      }
    }
    kk		=0;
    for(mm=0; mm < Mp1; mm++){
      sr	=yr[mm];
      si	=yi[mm];
      j		=mm;
      i		=kk+mm;
      while(j > 0){
	j--;
	i--;
	sr	-=Ar[i]*yr[j]-Ai[i]*yi[j];
	si	-=Ar[i]*yi[j]+Ai[i]*yr[j];
      }
      s		=Ar[kk+mm];
      yr[mm]	=s*sr;
      yi[mm]	=s*si;
      kk	+=Mp1;
    }
    m++;
  }
  /*}*/
  /* Calculation of $\maK=\maE-\maB^\dagger(\msa^\dagger)^{-1}\msa^{-1}\maC${*/
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1;
    yr	=Kr+j;
    yi	=Ki+j;
    ur	=Br+j;
    ui	=Bi+j;
    EZvr	=Cr;
    vi	=Ci;
    for(kk=0; kk < Mp1; kk++){
      sr	=0.;
      si	=0.;
      for(mm=0; mm < Mp1; mm++){
	sr	+=ur[mm]*(*EZvr)+ui[mm]*(*vi);
	si	+=ur[mm]*(*vi)-ui[mm]*(*EZvr);
	EZvr++;
	vi++;
      }
      yr[kk]	=-sr;
      yi[kk]	=-si;
    }
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1 && i < ES2Mp1; kk++){
      if(m != 0){
	if(mm != 0){
	  if(i){
	    yr[kk]	-=P*svr[i]*Tr[kk];
	    yi[kk]	-=P*svi[i]*Tr[kk];
	  }
	}
      }
      i++;
      mm++;
    }

    i	=1;
    kk	=k;
    mm	=kk+M1;
    while(kk > 0 && i < ES2Mp1){
      kk--;
      mm--;
      if(m != 0){
	if(mm != 0){
	  yr[kk]	-=P*svr[i]*Tr[kk];
	  yi[kk]	+=P*svi[i]*Tr[kk];
	}
      }
      i++;
    }
    m++;
  }
  /*}*/
  /* Calculation of $\maG=\maH-\maC^\dagger\msa^{-1\dagger}\msa^{-1}\maC$ {*/
  m	=M1;
  for(k=0; k < Mp1; k++){
    j	=k*Mp1;
    yr	=Gr+j;
    yi	=Gi+j;
    ur	=Cr+j;
    ui	=Ci+j;
    EZvr	=Cr+j;
    vi	=Ci+j;
    for(kk=k; kk < Mp1; kk++){
      sr	=0.;
      si	=0.;
      for(mm=0; mm < Mp1; mm++){
	sr	+=ur[mm]*(*EZvr)+ui[mm]*(*vi);
	si	+=ur[mm]*(*vi)-ui[mm]*(*EZvr);
	EZvr++;
	vi++;
      }
      yr[kk]	=-sr;
      yi[kk]	=-si;
    }
    i	=0;
    mm	=k+M1;
    for(kk=k; kk < Mp1 && i < ES2Mp1; kk++){
      if(m != 0){
	if(mm != 0){
	  if(i){
	    yr[kk]	+=P*Tr[k]*(dgm*svr[i]-gm*dsvr[i])*Tr[kk];
	    yi[kk]	+=P*Tr[k]*(dgm*svi[i]-gm*dsvi[i])*Tr[kk];
	  }
	  if(i == 0){
	    yr[kk]	-=N*P*dsvr[0]*Tr[k]*Tr[k]/m;
	  }
	}
      }
      if(i){
	j	=Mp1*kk+k;
	Gr[j]	= yr[kk];
	Gi[j]	=-yi[kk];
      }
      i++;
      mm++;
    }
    m++;
  }
  /*}*/
  /* Calculation of $\vaX^1=\maF^{-1}\vaY$ {*/
  xr	=X1r;
  xi	=X1i;
  for(k=0; k < Mp1; k++){
    *xr++	=Tr[k];
    *xi++	=0.;
  }
  CholBksbComp(Fr,Fi,Mp1,X1r,X1i);
  /*}*/
  /* Calculation of $a21=\vaX^1_m ,a11=\maK^\dagger\vaX^1{*/
  k	=Mres-M1;
  sr	=0.;
  si	=0.;
  a21	=X1r[k];
  xr	=X1r;
  xi	=X1i;
  mm	=k;
  for(kk=0; kk < Mp1; kk++){
    sr	+=Kr[mm]*(*xr)+Ki[mm]*(*xi);
    si	+=Kr[mm]*(*xi)-Ki[mm]*(*xr);
    xr++;
    xi++;
    mm	+=Mp1;
  }
  a11	=sr;
  /*}*/

  /* Calculation of $\vaX^2=-\maF^{-1}\maK\vaI$ {*/
  xr	=X2r;
  xi	=X2i;
  for(k=0; k < Mp1; k++){
    sr	=0.;
    si	=0.;
    ur	=Kr+k*Mp1;
    ui	=Ki+k*Mp1;
    EZvr	=Tr;
    for(kk=0; kk < Mp1; kk++){
      sr	+=(*ur)*(*EZvr);
      si	+=(*ui)*(*EZvr);
      ur++;
      ui++;
      EZvr++;
    }
    *xr++	=-sr;
    *xi++	=-si;
  }
  CholBksbComp(Fr,Fi,Mp1,X2r,X2i);
  /*}*/
  /* Calculation of $a22=\vaX^2_m ,a12=\maK^\dagger\vaX^2{*/
  k	=Mres-M1;
  sr	=0.;
  si	=0.;
  xr	=X2r;
  xi	=X2i;
  a22	=X2r[k];
  mm	=k;
  for(kk=0; kk < Mp1; kk++){
    sr	+=Kr[mm]*(*xr)+Ki[mm]*(*xi);
    si	+=Kr[mm]*(*xi)-Ki[mm]*(*xr);
    xr++;
    xi++;
    mm	+=Mp1;
  }
  a12	=sr+Gr[k*Mp1+k];
  printf("m=%2d x=%10.3e a11=%10.3e %10.3e %10.3e %10.3e\n"
	 ,Mres,x,a11,a12,a21,a22);
  /*}*/
  a22	-=EZcr2*dgm;
  a22	=a12*a21+a22*a22;
  if(a22 <= 0.){
    printf("m=%2d x=%10.3e %10.3e - Unstable resonance\n",mRes[kRes],x,a22);
    return(1);
  }
  else{
    a22	=sqrt(a22);
    a12	/=fabs(dgm);
    PM1	=EZcr2+a12;
    PM2	=EZcr2-a12;
  }
  return(0);
}

int ESStInitSmallResSolution(int n,double *xcur)
{
  int k,km,i;
  double dgm,gyR,YRr,YRi,x,s;
  double *pr,*pi,*Pr,*Pi,*yr,*yi,*Yr,*Yi;

  x	=aRes[kRes]-(*xcur);
  *xcur	=aRes[kRes]+x;

  km	=Mres-M1;
  k	=Mp1*km;
  pr	=gyr+k;
  pi	=gyi+k;
  yr	=aYr+k;
  yi	=aYi+k;
  for(k=0; k < Mp1; k++){
    yr[k]	=0.;
    yi[k]	=0.;
  }
  dgm	=dgmRes[kRes]*PM1;
  gyR	=pow(x,PM1);
  pr[km]	=gyR;
  pi[km]	=0.;
  YRr		=(dgm-X2r[km])*gyR/X1r[km];
  YRi		=-X2i[km]*gyR/X1r[km];
  yr[km]	=YRr/(dgm*x);
  yi[km]	=YRi/(dgm*x);

  dgm	=1./dgm;
  for(k=0; k < Mp1; k++){
    if(k != km){
      pr[k]	=(X1r[k]*YRr-X1i[k]*YRi+X2r[k]*gyR)*dgm;
      pi[k]	=(X1r[k]*YRi+X1i[k]*YRr+X2i[k]*gyR)*dgm;
    }
  }
  
  x	=1.;
  s	=1.;
  dgm	=PM1;
  for(i=0; i <= n; i++){
    Pr	=yh+neq*i+Mp1*km;
    Pi	=pr+MMp1;
    Yr	=pi+MMp1;
    Yi	=yr+MMp1;
    for(k=0; k < Mp1; k++){
      *Pr++	=pr[k]*x;
      *Pi++	=pi[k]*x;
      *Yr++	=yr[k]*s;
      *Yi++	=yi[k]*s;
    }
    x	*=dgm/(i+1);
    dgm	-=1.;
    s	*=dgm/(i+1);
  }

  return(0);
}

#ifdef XWIN
int ESStDispl(int Opt,double x)
{
  static double a[3],q[3];
  if(Opt == 0){
    int i;
    double h;
    h	=(qmx-qmn)*0.1;
    i	=N*qmn-dM2;
    q[1]=(double)i/N;
    i	=N*qmx+dM2;
    q[0]=(double)i/N;
    a[1]=ESsa[0];
    a[0]=ESsa[ESNa];

    q[1]=0.;
    q[0]=qmx+1.;

    Scale2D(3,2,a,q,2,6,2);
    Plot2d(3,ESsa,ESsq,ESNa1,6,0,1,0);
    for(i=0; i < nRes; i++){
      a[0]	=aRes[i];
      a[1]	=aRes[i];
      a[2]	=a[1];
      q[0]	=qRes[i]-h;
      q[1]	=qRes[i]+h;
      q[2]	=q[1];
      Plot2d(3,a,q,2,6,0,4,0);
    }
    q[1]	=ESsq[0];
    a[1]	=ESsa[0];
  }
  else{
    ESSetSplA(x);
    Plot2d(3,a+1,q+1,2,6,0,15,0);
    q[0]	=q[1];
    a[0]	=a[1];
    a[1]	=x;
    a[2]	=x;
    splRA(q+1,NULL,ESsq,ESsq2a);
    q[2]	=q[1]-0.2;
    Plot2d(3,a,q,3,6,0,14,0);
  }
  CbFlush();
  return(0);
}
#endif

int ESStSetVatol(double x)
{
  int j,k;
  double A;
  double *Vr,*Vi,*EZvr,*vi; 
  
  A	=reltol;
  if(M1 < 0){
    double *Vr1,*Vi1,*vr1,*vi1; 
    j	=-M1*Mp1;
    EZvr	=vatol+j;
    vi	=EZvr+MMp1;
    Vr	=vi+MMp1;
    Vi	=Vr+MMp1;
    vr1	=EZvr-1;
    vi1	=vi-1;
    Vr1	=Vr-1;
    Vi1	=Vi-1;
    for(j=0; j < Mp1; j++){
      *EZvr++	=A;
      *vi++	=A;
      *Vr++	=A;
      *Vi++	=A;
    }
    k	=0;
    while(k < -M1){
      A	*=x;
      k++;
      if(A < 1e-100){
	A	=1e-100;
      }
      for(j=0; j < Mp1; j++){
	*EZvr++	=A;
	*vi++	=A;
	*Vr++	=A;
	*Vi++	=A;
	*vr1--	=A;
	*vi1--	=A;
	*Vr1--	=A;
	*Vi1--	=A;
      }
    }
    A	*=x;
  }
  else{
    EZvr	=vatol;
    vi	=EZvr+MMp1;
    Vr	=vi+MMp1;
    Vi	=Vr+MMp1;
    k	=0;
    while(k < M1){
      A	*=x;
      k++;
    }
    k--;
  }
  while(k < M2){
    if(A < 1e-100){
      A	=1e-100;
    }
    for(j=0; j < Mp1; j++){
      *EZvr++	=A;
      *vi++	=A;
      *Vr++	=A;
      *Vi++	=A;
    }
    A	*=x;
    k++;
  }
  return(0);
}

int ESMHDStCheckLowN(int nord,int mc,double x)
{
  int jc,i,j,j1,jj;
  double *Pr,*Pi,*pr,*pi,*Qr,*Qi,*qr,*qi;
  double *Yr,*Yi,*yr,*yi,*Zr,*Zi,*zr,*zi;
  double s,sr,si,gar,gai;
  static double dr0,di0,dr1,di1;

  KRec	=0;
  jc	=mc-M1;
  /* Low end harmonics */
  Pr	=gyr;
  Pi	=gyi;
  Yr	=aYr;
  Yi	=aYi;
  dr1	=1.;
  di1	=0.;
  for(j=0; j <= jc; j++){
    sr	=Pr[j];
    si	=-Pi[j];
    gar	=dr1*sr+di1*si;
    di1	=di1*sr-dr1*si;
    dr1	=gar;
    s	=dr1*dr1+di1*di1;
    if(s > 1e+30){
      dr1	*=1e-15;
      di1	*=1e-15;
    }	
    if(s < 1e-30){
      dr1	*=1e+15;
      di1	*=1e+15;
    }	
    s	=1./(sr*sr+si*si);
    sr	*=s;
    si	*=s;
    Ar[KRec]	=1.;
    Ai[KRec]	=0.;
    KRec++;
    Qr	=Pr+Mp1;
    Qi	=Pi+Mp1;
    Zr	=Yr+Mp1;
    Zi	=Yi+Mp1;
    for(j1=j+1; j1 < Mp1; j1++){
      gar	=Qi[j]*si-Qr[j]*sr;
      gai	=-Qi[j]*sr-Qr[j]*si;
      Ar[KRec]	=gar;
      Ai[KRec]	=gai;
      KRec++;
      pr	=Pr;
      pi	=Pi;
      yr	=Yr;
      yi	=Yi;
      for(jj=0; jj < Mp1; jj++){
	*Qr++	+=gar*(*pr)-gai*(*pi);
	*Qi++	+=gar*(*pi)+gai*(*pr);
	pr++;
	pi++;
	*Zr++	+=gar*(*yr)-gai*(*yi);
	*Zi++	+=gar*(*yi)+gai*(*yr);
	yr++;
	yi++;
      }
      for(i=0; i <= nord; i++){
	pr	=yh+neq*i+Mp1*j;
	pi	=pr+MMp1;
	yr	=pi+MMp1;
	yi	=yr+MMp1;
	qr	=yh+neq*i+Mp1*j1;
	qi	=qr+MMp1;
	zr	=qi+MMp1;
	zi	=zr+MMp1;
	for(jj=0; jj < Mp1; jj++){
	  *qr++	+=gar*(*pr)-gai*(*pi);
	  *qi++	+=gar*(*pi)+gai*(*pr);
	  pr++;
	  pi++;
	  *zr++	+=gar*(*yr)-gai*(*yi);
	  *zi++	+=gar*(*yi)+gai*(*yr);
	  yr++;
	  yi++;
	}
      }
    }
    Pr	+=Mp1;
    Pi	+=Mp1;
    Yr	+=Mp1;
    Yi	+=Mp1;
  }
  /* High end harmonics */
  j	=MMp1-Mp1;
  Pr	=gyr+j;
  Pi	=gyi+j;
  Yr	=aYr+j;
  Yi	=aYi+j;
  for(j=Mp1-1; j >= jc; j--){
    sr	=Pr[j];
    si	=-Pi[j];
    if(j > jc){
      gar	=dr1*sr+di1*si;
      di1	=di1*sr-dr1*si;
      dr1	=gar;
      s	=dr1*dr1+di1*di1;
      if(s > 1e+30){
	dr1	*=1e-15;
	di1	*=1e-15;
      }	
      if(s < 1e-30){
	dr1	*=1e+15;
	di1	*=1e+15;
      }	
    }
    s	=1./(sr*sr+si*si);
    sr	*=s;
    si	*=s;
    Ar[KRec]	=1.;
    Ai[KRec]	=0.;
    KRec++;
    j1	=j;
    while(j1 > 0){
      j1--;
      jj	=Mp1*j1;
      Qr	=gyr+jj;
      Qi	=gyi+jj;
      Zr	=aYr+jj;
      Zi	=aYi+jj;
      gar	=Qi[j]*si-Qr[j]*sr;
      gai	=-Qi[j]*sr-Qr[j]*si;
      Ar[KRec]	=gar;
      Ai[KRec]	=gai;
      KRec++;
      pr	=Pr;
      pi	=Pi;
      yr	=Yr;
      yi	=Yi;
      for(jj=0; jj < Mp1; jj++){
	*Qr++	+=gar*(*pr)-gai*(*pi);
	*Qi++	+=gar*(*pi)+gai*(*pr);
	pr++;
	pi++;
	*Zr++	+=gar*(*yr)-gai*(*yi);
	*Zi++	+=gar*(*yi)+gai*(*yr);
	yr++;
	yi++;
      }
      for(i=0; i <= nord; i++){
	pr	=yh+neq*i+Mp1*j;
	pi	=pr+MMp1;
	yr	=pi+MMp1;
	yi	=yr+MMp1;
	qr	=yh+neq*i+Mp1*j1;
	qi	=qr+MMp1;
	zr	=qi+MMp1;
	zi	=zr+MMp1;
	for(jj=0; jj < Mp1; jj++){
	  *qr++	+=gar*(*pr)-gai*(*pi);
	  *qi++	+=gar*(*pi)+gai*(*pr);
	  pr++;
	  pi++;
	  *zr++	+=gar*(*yr)-gai*(*yi);
	  *zi++	+=gar*(*yi)+gai*(*yr);
	  yr++;
	  yi++;
	}
      }
    }
    Pr	-=Mp1;
    Pi	-=Mp1;
    Yr	-=Mp1;
    Yi	-=Mp1;
  }
#ifdef H
  sr	=Pr[j];
  si	=-Pi[j];
  gar	=dr1*sr+di1*si;
  di1	=di1*sr-dr1*si;
  dr1	=gar;
  s	=dr1*dr1+di1*di1;
  if(s > 1e+30){
    dr1	*=1e-15;
    di1	*=1e-15;
  }	
  if(s < 1e-30){
    dr1	*=1e+15;
    di1	*=1e+15;
  }	
#endif
  j	=Mp1*jc+jc;
  sr	=gyr[j];
  si	=gyi[j];
  if(kCheck == 0){
    dr0		=dr1;
    di0		=di1;
    aYmr	=aYr[j];
    aYmi	=aYi[j];
    gymr	=sr;
    gymi	=si;
  }
  ESsgnY	=aYr[j]*aYmr+aYi[j]*aYmi;
  ESsgnX	=gymr*gyr[j]+gymi*gyi[j];
  aYmr	=aYr[j];
  aYmi	=aYi[j];
  gymr	=sr;
  gymi	=si;

  ESdet	=dr1*dr0+di1*di0 > 0. ? 1. : -1.;
  EZout("siidddd","SgnXYs",mc,kCheck,x,ESsgnX,ESsgnY,ESdet);

  aYmr	=aYr[j];
  aYmi	=aYi[j];
  gymr	=sr;
  gymi	=si;
  kCheck++;
  return(0);
}

int ESMHDStShiftUpLowN()
{
  int i,j,m,k,kk;
  double x,h,H,A;
  double *Pr,*Pi,*pr,*pi;
  double *Yr,*Yi,*yr,*yi;

  x	=rwork[12];
  /* Getting metric coefficients{*/
  ESSetSplA(x);
  splRA(gqqr,NULL,ESg22c,ESg22c2);
  gqqi[0]	=0.;
  m	=0;
  for(k=1; k < ES2Mp1; k++){
    m	+=ESNa1;
    splRA(gqqr+k,NULL,ESg22c+m,ESg22c2+m);
    splRA(gqqi+k,NULL,ESg22s+m,ESg22s2+m);
  }

  M1++;
  M2++;
  Pr	=gyr;
  Pi	=gyi;
  Yr	=aYr;
  Yi	=aYi;
  pr	=Pr+Mp1;
  pi	=Pi+Mp1;
  yr	=Yr+Mp1;
  yi	=Yi+Mp1;
  for(j=1; j < Mp1; j++){
    pr++;
    pi++;
    yr++;
    yi++;
    for(k=1; k < Mp1; k++){
      *Pr++	=*pr++;
      *Pi++	=*pi++;
      *Yr++	=*yr++;
      *Yi++	=*yi++;
    }
    *Pr++	=0.;
    *Pi++	=0.;
    *Yr++	=0.;
    *Yi++	=0.;
  }

  H	=x;
  j	=1;
  while(2*j <= M2){
    H	*=H;
    j	*=2;
  }
  while(j < M2){
    j++;
    H	*=x;
  }
  if(H < 1e-100){
    H	=1e-100;
  }
  {
    double bK0,Er,Ei;
    double br,bi,sr,si,tr,ti;
    
    bK0	=ESg22c[0];
    k	=2*ESNa1;
    Er	=ESg22c[k];
    Ei	=-ESg22s[k];

    j	=M2;
    Pr	-=M1;
    Pi	-=M1;
    Yr	-=M1;
    Yi	-=M1;
    k	=M2;
    /* Decomposition loop */
    Pr[k]	=H;
    Pi[k]	=0.;
    Yr[k]	=0.;
    Yi[k]	=0.;
    k--;
    while(k >= M1){
      Pr[k]	=0.;
      Pi[k]	=0.;
      k		-=2;
    }
    k		=M2-2;
    while(k >= M1){
      kk	=k+2;
      A		=(j-k)*(j+k)*bK0;
      bi	=(j+k)*(j+kk);
      br	=bi*Er;
      bi	*=-Ei;
      si	=(j-k)*(j-k+2);
      sr	=si*Er;
      si	*=Ei;
      tr	=A+br*Yr[kk]-bi*Yi[kk];
      ti	=br*Yi[kk]+bi*Yr[kk];
      A		=1./(tr*tr+ti*ti);
      tr	*=-A;
      ti	*=A;
      Yr[k]	=sr*tr-si*ti;
      Yi[k]	=sr*ti+si*tr;
      sr	=br*Pr[kk]-bi*Pi[kk];
      si	=br*Pi[kk]+bi*Pr[kk];
      Pr[k]	=sr*tr-si*ti;
      Pi[k]	=sr*ti+si*tr;
      k		-=2;
    }
    /* Backsubstitution loop */
    k	+=4;
    kk	=k-2;
    while(k < M2){
      Pr[k]	+=Pr[kk]*Yr[k]-Pi[kk]*Yi[k];
      Pi[k]	+=Pr[kk]*Yi[k]+Pi[kk]*Yr[k];
      kk	=k;
      k	+=2;
    }
    /* Calculation of $aY${*/
    Yr	+=M1;
    Yi	+=M1;
    Pr	+=M1;
    Pi	+=M1;
    for(k=0; k < Mp1; k++){
      Yr[k]	=0.;
      Yi[k]	=0.;
    }
    Yr[k-1]	=M2*H;
#ifdef H
    for(k=0; k < Mp1; k++){
      kk	=k;
      br	=M2*Pr[k];
      bi	=M2*Pi[k];
      for(i=0; i < ES2Mp1 && kk < Mp1; i++){
	Yr[kk]	+=gqqr[i]*br-gqqi[i]*bi;
	Yi[kk]	+=gqqr[i]*bi+gqqi[i]*br;
	kk++;
      }
      kk	=k;
      for(i=1; i < ES2Mp1 && kk > 0; i++){
	kk--;
	Yr[kk]	+=gqqr[i]*br+gqqi[i]*bi;
	Yi[kk]	+=gqqr[i]*bi-gqqi[i]*br;
      }
    }
#endif
    /*}*/
  }
  for(k=0; k < Mp1; k++){
    Pr[k]	=0.;
    Pi[k]	=0.;
  }

  pr	=vatol;
  pi	=pr+MMp1;
  yr	=pi+MMp1;
  yi	=yr+MMp1;
  for(j=1; j < Mp1; j++){
    A	=pr[Mp1];
    for(k=0; k < Mp1; k++){
      *pr++	=A;
      *pi++	=A;
      *yr++	=A;
      *yi++	=A;
    }
  }
  H	*=reltol;
  for(j=0; j < Mp1; j++){
    *pr++	=H;
    *pi++	=H;
    *yr++	=H;
    *yi++	=H;
  }

  st2DVrMHD(&neq,&x,gyr,dgyr);
  H	=1.;
  h	=rwork[11]/x;
  kk	=iwork[14]+1;
  m	=M2;
  for(i=0; i < kk; i++){
    Pr	=yh+neq*i;
    Pi	=Pr+MMp1;
    Yr	=Pi+MMp1;
    Yi	=Yr+MMp1;
    pr	=Pr+Mp1;
    pi	=Pi+Mp1;
    yr	=Yr+Mp1;
    yi	=Yi+Mp1;
    for(j=1; j < Mp1; j++){
      pr++;
      pi++;
      yr++;
      yi++;
      for(k=1; k < Mp1; k++){
	*Pr++	=*pr++;
	*Pi++	=*pi++;
	*Yr++	=*yr++;
	*Yi++	=*yi++;
      }
      *Pr++	=0.;
      *Pi++	=0.;
      *Yr++	=0.;
      *Yi++	=0.;
    }
    if(i == 0){
      pr	=gyr+MMp1-Mp1;
      pi	=pr+MMp1;
      A		=1.;
    }
    else{
      pr	=dgyr+MMp1-Mp1;
      pi	=pr+MMp1;
      if(i == 1){
	A	=rwork[11];
      }
    }
    yr	=pi+MMp1;
    yi	=yr+MMp1;
    for(j=0; j < Mp1; j++){
      *Pr++	=A*(*pr++);
      *Pi++	=A*(*pi++);
      *Yr++	=H*(*yr++);
      *Yi++	=H*(*yi++);
    }
    H	*=h*m/(i+1);
    if(i){
      A	*=h*m/(i+1);
    }
    if(m){
      m--;
    }
  }
  return(0);
}

int ESMHDStShiftDnLowN()
{
  int i,j,m,k,kk;
  double br,bi,x,h,H,A;
  double *Pr,*Pi,*pr,*pi;
  double *Yr,*Yi,*yr,*yi;

  x	=rwork[12];
  /* Getting metric coefficients{*/
  ESSetSplA(x);
  splRA(gqqr,NULL,ESg22c,ESg22c2);
  gqqi[0]	=0.;
  m	=0;
  for(k=1; k < ES2Mp1; k++){
    m	+=ESNa1;
    splRA(gqqr+k,NULL,ESg22c+m,ESg22c2+m);
    splRA(gqqi+k,NULL,ESg22s+m,ESg22s2+m);
  }

  M1--;
  M2--;
  j	=MMp1-1;
  Pr	=gyr+j;
  Pi	=gyi+j;
  Yr	=aYr+j;
  Yi	=aYi+j;
  pr	=Pr-Mp1;
  pi	=Pi-Mp1;
  yr	=Yr-Mp1;
  yi	=Yi-Mp1;
  for(j=1; j < Mp1; j++){
    pr--;
    pi--;
    yr--;
    yi--;
    for(k=1; k < Mp1; k++){
      *Pr--	=*pr--;
      *Pi--	=*pi--;
      *Yr--	=*yr--;
      *Yi--	=*yi--;
    }
    *Pr--	=0.;
    *Pi--	=0.;
    *Yr--	=0.;
    *Yi--	=0.;
  }
  Pr	=gyr;
  Pi	=gyi;
  Yr	=aYr;
  Yi	=aYi;
  H	=1.;
  j	=0;
  k	=abs(M1);
  if(k){
    H	=x;
    j	=1;
    while(2*j <= k){
      H	*=H;
      j	*=2;
    }
    while(j < k){
      j++;
      H	*=x;
    }
  }
  if(H < 1e-100){
    H	=1e-100;
  }
  if(M1 < 0){
    double bK0,Er,Ei;
    double sr,si,tr,ti;
    
    bK0	=ESg22c[0];
    k	=2*ESNa1;
    Er	=ESg22c[k];
    Ei	=-ESg22s[k];

    j	=-M1;
    Pr	-=M1;
    Pi	-=M1;
    Yr	-=M1;
    Yi	-=M1;
    for(k=M2; k > j; k--){
      Pr[k]	=0.;
      Pi[k]	=0.;
    }
    /* Decomposition loop */
    Pr[k]	=0.;
    Pi[k]	=0.;
    Yr[k]	=0.;
    Yi[k]	=0.;
    k--;
    while(k >= M1){
      Pr[k]	=0.;
      Pi[k]	=0.;
      k		-=2;
    }
    k		=-M1-2;
    while(k >= M1){
      kk	=k+2;
      A		=(j-k)*(j+k)*bK0;
      bi	=(j+k)*(j+kk);
      br	=bi*Er;
      bi	*=-Ei;
      si	=(j-k)*(j-k+2);
      sr	=si*Er;
      si	*=Ei;
      tr	=A+br*Yr[kk]-bi*Yi[kk];
      ti	=br*Yi[kk]+bi*Yr[kk];
      A		=1./(tr*tr+ti*ti);
      tr	*=-A;
      ti	*=A;
      Yr[k]	=sr*tr-si*ti;
      Yi[k]	=sr*ti+si*tr;
      Pr[k]	=0.;
      Pi[k]	=0.;
      k		-=2;
    }
    /* Backsubstitution loop */
    k	+=4;
    kk	=k-2;
    Pr[kk]	=H;
    while(k < M2){
      Pr[k]	+=Pr[kk]*Yr[k]-Pi[kk]*Yi[k];
      Pi[k]	+=Pr[kk]*Yi[k]+Pi[kk]*Yr[k];
      kk	=k;
      k	+=2;
    }
    /*}*/
  }
  else{
    Pr[0]	=H;
    Pi[0]	=0.;
    for(k=1; k < Mp1; k++){
      Pr[k]	=0.;
      Pi[k]	=0.;
    }
  }

  /* Calculation of $aY${*/
  j	=abs(M1);
  Yr	=aYr;
  Yi	=aYi;
  Pr	=gyr;
  Pi	=gyi;
  for(k=0; k < Mp1; k++){
    Yr[k]	=0.;
    Yi[k]	=0.;
  }
  Yr[k-1]	=j == 0 ? 1. : j*H;
#ifdef H
  for(k=0; k < Mp1; k++){
    kk	=k;
    br	=j*Pr[k];
    bi	=j*Pi[k];
    for(i=0; i < ES2Mp1 && kk < Mp1; i++){
      Yr[kk]	+=gqqr[i]*br-gqqi[i]*bi;
      Yi[kk]	+=gqqr[i]*bi+gqqi[i]*br;
      kk++;
    }
    kk	=k;
    for(i=1; i < ES2Mp1 && kk > 0; i++){
      kk--;
      Yr[kk]	+=gqqr[i]*br+gqqi[i]*bi;
      Yi[kk]	+=gqqr[i]*bi-gqqi[i]*br;
    }
  }
#endif
  for(k=0; k < Mp1; k++){
    Pr[k]	=0.;
    Pi[k]	=0.;
  }

  pr	=vatol+MMp1-1;
  pi	=pr+MMp1;
  yr	=pi+MMp1;
  yi	=yr+MMp1;
  for(j=1; j < Mp1; j++){
    A	=pr[-Mp1];
    for(k=0; k < Mp1; k++){
      *pr--	=A;
      *pi--	=A;
      *yr--	=A;
      *yi--	=A;
    }
  }
  H	*=reltol;
  for(j=0; j < Mp1; j++){
    *pr--	=H;
    *pi--	=H;
    *yr--	=H;
    *yi--	=H;
  }

  st2DVrMHD(&neq,&x,gyr,dgyr);
  H	=1.;
  h	=rwork[11]/x;
  kk	=iwork[14]+1;
  m	=M2;
  for(i=0; i < kk; i++){
    Pr	=yh+neq*i+MMp1-1;
    Pi	=Pr+MMp1;
    Yr	=Pi+MMp1;
    Yi	=Yr+MMp1;
    pr	=Pr-Mp1;
    pi	=Pi-Mp1;
    yr	=Yr-Mp1;
    yi	=Yi-Mp1;
    for(j=1; j < Mp1; j++){
      pr--;
      pi--;
      yr--;
      yi--;
      for(k=1; k < Mp1; k++){
	*Pr--	=*pr--;
	*Pi--	=*pi--;
	*Yr--	=*yr--;
	*Yi--	=*yi--;
      }
      *Pr--	=0.;
      *Pi--	=0.;
      *Yr--	=0.;
      *Yi--	=0.;
    }
    if(i == 0){
      pr	=gyr+Mp1-1;
      pi	=pr+MMp1;
      A		=1.;
    }
    else{
      pr	=dgyr+Mp1-1;
      pi	=pr+MMp1;
      if(i == 1){
	A	=rwork[11];
      }
    }
    yr	=pi+MMp1;
    yi	=yr+MMp1;
    for(j=0; j < Mp1; j++){
      *Pr--	=A*(*pr--);
      *Pi--	=A*(*pi--);
      *Yr--	=H*(*yr--);
      *Yi--	=H*(*yi--);
    }
    H	*=h*m/(i+1);
    if(i){
      A	*=h*m/(i+1);
    }
    if(m){
      m--;
    }
  }
  return(0);
}

int ESStSetXres(double *xres, double x)
{
  while(iRes < nRes && aRes[iRes] < x){
    iRes++;
  }
  if(iRes == nRes){
    *xres	=ESsa[ESNa];
    return(0);
  }
  *xres	=aRes[iRes];
  return(1);
}

int ESStSetXcenter(int *mcen, double *xcen, int i)
{
  if(i+1 < nRes){
    *mcen	=mRes[i];
    *xcen	=EZcr2*(aRes[i]+aRes[i+1]);
    Mres	=mRes[i];
    kRes	=i;
    ESResAsymptotics();
    return(0);
  }
  else{
    if(i+1 == nRes){
      *mcen	=mRes[i];
      Mres	=mRes[i];
      kRes	=i;
      ESResAsymptotics();
      if(N*ESsq[ESNa] < mRes[i]+0.5){
	*xcen	=ESsa[ESNa]+ESsa[1];
      }
      else{
	*xcen	=aRes[i]+EZcr2*(ESsa[ESNa]-aRes[i])/(N*ESsq[ESNa]-mRes[i]);
      }
      return(0);
    }
    else{
      if(i == nRes){
	*mcen	=N*ESsq[ESNa]+1;
	*xcen	=ESsa[ESNa]+ESsa[1];
	Mres	=0;
	kRes	=0;
	return(1);
      }
    }
  }
  return(0);
}

int ESStSetCheckPoint(double *xchk,double x,double dx)
{
  static double h=0.;

  *xchk=dx*((int)(x/dx)+1);
  if(*xchk > rwork[0] && FlagRes == 0){
    *xchk	=rwork[0];
  }
  if(*xchk >= aResChck && Mres){
    *xchk	=aResChck;
    FlagRes	=1;
  }

  return(0);
}

#ifdef XWIN
int ES2DStSolv(int Opt,double a1)
{
  double x,xc,xChk,sChk=0.02;
  int nR,mc;

  nR	=aRes[nRes-1] == ESsa[ESNa] ? nRes : nRes+1;

  x	=0.01;
  xChk	=x+sChk;
  iRes	=0;
  ESStSetXcenter(&mc,&xc,iRes);
  if(iRes < nRes){
    if(x > 0.1*aRes[0]){
      x		=0.1*aRes[0];
      xChk	=0.6*aRes[0];
    }
    if(xChk-x > sChk){
      xChk	=x+sChk;
    }
    aResChck	=aRes[0]-(aRes[0]-xChk)*0.001;
  }
  else{
    aResChck	=a1;
  }
  ESStSetXres(rwork+0,x);

  FlStChol	=0;
  
  mf		=10;
  itol		=2;	/* vector tolerance */
  itask		=4;
  itask		=5;
  iopt		=1;	/* optional input ON */
  rwork[6]	=0.;	/* Minimum step */

  rwork[2-1]	=0.;
  rwork[3-1]	=0.;
  rwork[4-1]	=0.;
  rwork[5-1]	=sChk;	/* the step size to be attempted on the first step. */
  rwork[6-1]	=0.1;	/* the maximum absolute step size allowed. */
  rwork[7-1]	=0.;	/* the minimum absolute step size allowed. */

  iwork[4]	=12;	/* Maximum order */
  iwork[5]	=500;
  iwork[6]	=10;
  reltol	=1.0e-6;
  abstol	=1e-12;
  istate	=1;
  
  if(ESInitStFundSolution(0,x) < 0){
    return(-1);
  }
  st2DVrMHD(&neq,&x,gyr,dgyr);

#ifdef XWIN
  ESStDispl(0,x);
  {
    double a[2],q[2];
    ESSetSplA(xChk);
    a[0]	=xChk;
    a[1]	=xChk;
    splRA(q,NULL,ESsq,ESsq2a);
    q[1]	=qmn-1.;
    printf("xChk=%10.3e %10.3e %10.3e\n",xChk,q[0],q[1]);
    Plot2d(3,a,q,2,6,0,12,0);
    ESSetSplA(xc);
    a[0]	=xc;
    a[1]	=xc;
    splRA(q,NULL,ESsq,ESsq2a);
    q[1]	=qmn-1.;
    Plot2d(3,a,q,2,6,0,4,0);
    CbFlush();
  }
#endif

  kCheck=0;
  NRec	=0;
  MRec	=0;
  if(Opt > 2){
    ESStKRec=0;
    fpy	=fopen("LowN.y","w+");
    fpt	=fopen("LowN.t","w+");
  }
  while(iRes < nR && x < a1){
    if(iRes < nRes){
      xout	=aRes[iRes];
      if(iRes){
	if(iRes+1 < nRes){
	  aResChck	=aRes[iRes+1]-aRes[iRes];
	}
	else{
	  aResChck	=fabs(aRes[iRes]-aRes[iRes-1]);
	}
	aResChck	=aRes[iRes]-0.001*aResChck;
      }
    }
    else{
      xout	=a1;
      aResChck	=a1;
    }
    rwork[0]	=xout;	/* tcrit = critical value of t which the solver
				   is not to overshoot. */
    if(rwork[0] > a1){
      xout 	=a1;
      rwork[0]	=a1;
    }
    if(iRes){
      ESStSetCheckPoint(&xChk,x,sChk);
    }
    FlagRes	=0;
    if(Opt > 2){
      fwrite(&x,sizeof(double),1,fpy);
      fwrite(gyr,sizeof(double),neq,fpy);
      NRec++;
    }
    while(x < rwork[0]){
      lsode_(st2DVrMHD,&neq,gyr,&x,&xout,&itol,&reltol,vatol,&itask,
	     &istate,&iopt,rwork,&lrw,iwork,&liw,NULL,&mf);
      EZout("sddd","x=",aRes[iRes],aResChck,x);
      if(istate < 0){
	if(Opt > 2){
	  fclose(fpy);
	  fclose(fpt);
	  fpy	=fopen("LowN.n","w");
	  if(fpy == NULL){
	    printf("file LowN.n cannot be opened by ES2DStSolv()\n");
	  }
	  else{
	    ESStNRec	=NRec;
	    ESStMRec	=MRec;
	    fwrite(&ESStN,sizeof(int),1,fpy);
	    fwrite(&ESStM1,sizeof(int),1,fpy);
	    fwrite(&ESStM2,sizeof(int),1,fpy);
	    fwrite(&ESStNRec,sizeof(int),1,fpy);
	    fwrite(&ESStMRec,sizeof(int),1,fpy);
	    fwrite(&ESStKRec,sizeof(int),1,fpy);
	    fclose(fpy);
	    printf("%6s%6s%6s%6s%6s%6s\n"
		   ,"NRec","MRec","KRec","N","M1","M2");
	    printf("%6d%6d%6d%6d%6d%6d\n"
		   ,ESStNRec,ESStMRec,ESStKRec,ESStN,ESStM1,ESStM2);
	  }
	}
	return(-1);
      }
      if(Opt > 2){
	fwrite(&x,sizeof(double),1,fpy);
	fwrite(gyr,sizeof(double),neq,fpy);
	NRec++;
      }

      {
	int j,k,kk,jj;
	double Wr,Wi,wr,wi,EZvr,vi,W[6];
	double *yr1,*yi1,*Yr1,*Yi1,*yr2,*yi2,*Yr2,*Yi2;
	
	jj	=0;
	for(k=0; k < 1; k++){
	  for(kk=k; kk < 1; kk++){
	    j	=(Mgy+k-M1)*Mp1;
	    yr1	=gyr+j;
	    yi1	=gyi+j;
	    Yr1	=aYr+j;
	    Yi1	=aYi+j;
	    j	=(MaY+kk-M1)*Mp1;
	    yr2	=gyr+j;
	    yi2	=gyi+j;
	    Yr2	=aYr+j;
	    Yi2	=aYi+j;
	    Wr	=0.;
	    Wi	=0.;
	    wr	=0.;
	    wi	=0.;
	    EZvr	=0.;
	    vi	=0.;
	    for(j=0; j < Mp1; j++){
	      Wr	+=(*yr1)*(*Yr2)+(*yi1)*(*Yi2)
		-(*Yr1)*(*yr2)-(*Yi1)*(*yi2);
	      Wi	+=(*yr1)*(*Yi2)-(*yi1)*(*Yr2)
		-(*Yr1)*(*yi2)+(*Yi1)*(*yr2);
	      wr	+=(*yr1)*(*Yr2)+(*yi1)*(*Yi2);
	      wi	+=(*yr1)*(*Yi2)-(*yi1)*(*Yr2);
	      EZvr	+=-(*Yr1)*(*yr2)-(*Yi1)*(*yi2);
	      vi	+=-(*Yi1)*(*yr2)+(*Yr1)*(*yi2);
	      yr1++;
	      yi1++;
	      Yr1++;
	      Yi1++;
	      yr2++;
	      yi2++;
	      Yr2++;
	      Yi2++;
	    }
#ifdef H
	    if(Mgy){
	      EZout("siiddddddd","W",Mgy,MaY,Wr,Wi,wr,wi
		  ,sqrt((Wr*Wr+Wi*Wi)/(wr*wr+wi*wi))
		  ,sqrt((Wr*Wr+Wi*Wi)/(EZvr*EZvr+vi*vi)),x);
	    }
	    else{
	      j	=(Mgy+k-M1)*Mp1;
	      
	      EZout("siidddddd","W",Mgy,MaY,Wr,Wi,wr,wi
		  ,gyr[j-M1],x);
	    }
#endif
	    W[jj]	=Wr;
	    jj++;
	  }
	}
      }
#ifdef XWIN
      ESStDispl(1,x);
      if(x >= xChk){
	double a[2],q[2];
	printf("=========mc=%d============\n",mc);
	ESMHDStCheckLowN(iwork[14],mc,x);
	ESStSetVatol(x);
	if(Opt > 2){
	  if(ESStKRec < KRec){
	    ESStKRec	=KRec;
	  }
	  fwrite(&KRec,sizeof(int),1,fpt);
	  fwrite(&MRec,sizeof(int),1,fpt);
	  fwrite(&NRec,sizeof(int),1,fpt);
	  fwrite(&mc,sizeof(int),1,fpt);
	  fwrite(&FlagRes,sizeof(int),1,fpt);
	  fwrite(Ar,sizeof(double),KRec,fpt);
	  fwrite(Ai,sizeof(double),KRec,fpt);
	  MRec++;
	}
	if(FlagRes){
	  ESStInitSmallResSolution(iwork[14],&x);
	  istate	=1;
	}
	else{
	  ESStSetCheckPoint(&xChk,x,sChk);
	}
	if(xChk <= ESsa[ESNa]){
	  ESSetSplA(xChk);
	  a[0]	=xChk;
	  a[1]	=xChk;
	  splRA(q,NULL,ESsq,ESsq2a);
	  q[1]	=qmn-1.;
	  Plot2d(3,a,q,2,6,0,12,0);
	  a[0]	=x;
	  a[1]	=x;
	  q[0]	=qmx+EZcr2;
	  q[1]	=q[0]+EZcr2*ESdet;
	  Plot2d(3,a,q,2,6,0,14,0);
	  CbFlush();
	}
      }
      if(x > xc){
	int mc0;
	double a[2],q[2];
	printf("-------mc=%d-----------\n",mc);
	ESMHDStCheckLowN(iwork[14],mc,x);
	ESStSetVatol(x);
	if(Opt > 2){
	  if(ESStKRec < KRec){
	    ESStKRec	=KRec;
	  }
	  fwrite(&KRec,sizeof(int),1,fpt);
	  fwrite(&MRec,sizeof(int),1,fpt);
	  fwrite(&NRec,sizeof(int),1,fpt);
	  fwrite(&mc,sizeof(int),1,fpt);
	  fwrite(&FlagRes,sizeof(int),1,fpt);
	  fwrite(Ar,sizeof(double),KRec,fpt);
	  fwrite(Ai,sizeof(double),KRec,fpt);
	  MRec++;
	}
	mc0	=mc;
	ESStSetXcenter(&mc,&xc,iRes);
	if(Opt == 1){
	  if(mc > mc0){
	    Mgy++;
	    MaY++;
	    ESMHDStShiftUpLowN();
	  }
	  if(mc < mc0){
	    Mgy--;
	    MaY--;
	    ESMHDStShiftDnLowN();
	  }
	}
	ESSetSplA(xc);
	a[0]	=xc;
	a[1]	=xc;
	splRA(q,NULL,ESsq,ESsq2a);
	q[1]	=qmn-1.;
	Plot2d(3,a,q,2,6,0,4,0);
	a[0]	=x;
	a[1]	=x;
	q[0]	=qmx+EZcr2;
	q[1]	=q[0]+EZcr2*ESdet;
	Plot2d(3,a,q,2,6,0,4,0);
	CbFlush();
      }
#endif
    }
    iRes++;
  }

  if(Opt > 2){
    fclose(fpy);
    fclose(fpt);
    fpy	=fopen("LowN.n","w");
    if(fpy == NULL){
      printf("file LowN.n cannot be opened by ES2DStSolv()\n");
    }
    else{
      ESStNRec	=NRec;
      ESStMRec	=MRec;
      fwrite(&ESStN,sizeof(int),1,fpy);
      fwrite(&ESStM1,sizeof(int),1,fpy);
      fwrite(&ESStM2,sizeof(int),1,fpy);
      fwrite(&ESStNRec,sizeof(int),1,fpy);
      fwrite(&ESStMRec,sizeof(int),1,fpy);
      fwrite(&ESStKRec,sizeof(int),1,fpy);
      fclose(fpy);
      printf("%6s%6s%6s%6s%6s%6s\n"
	     ,"NRec","MRec","KRec","N","M1","M2");
      printf("%6d%6d%6d%6d%6d%6d\n"
	     ,ESStNRec,ESStMRec,ESStKRec,ESStN,ESStM1,ESStM2);
    }
  }
  return(0);
}

int ESMHDStabLowN()
{
  static int Fh=0,Fc=0,Fd=0,k0=0,k1=11,i0=0,i1=21;
  static int mV=0,m1=0,m2=1,FlV=0;
  static double sa[2]={0.,1.};
  
  double x[2],y[2],Y[2];
  int kW;
  double Wx[500],W1[33][500],*Wy;
  int i,Fl,Opt;

  R0	=ESaR0[0];
  rR0	=1./R0;
  aRes	=ESsa;
  qmn	=ESsq[0];
  qmx	=ESsq[0];
  for(i=1; i < ESNa1; i++){
    if(qmn > ESsq[i]){
      qmn	=ESsq[i];
    }
    if(qmx < ESsq[i]){
      qmx	=ESsq[i];
    }
  }

  ESStNRec	=0;
  ESStMRec	=0;
  ESStKRec	=0;
  Opt	=0;
  Mgy	=N*ESsq[0];
  MaY	=N*ESsq[0];
#ifndef Tbl19a
  if(N < 0){
    N	=-N;
  }
  if(N == 0){
    N	=1;
  }
  if(dM1 < 0){
    dM1	=abs(dM1);
  }
  if(dM2 < 0){
    dM2	=abs(dM2);
  }

  ESst0=clock();
  Fl	=0;
  if((ESmemlev&0x04000000)){
    ESDeInitMHDMetrTnsr();
    ESmemlev -=0x04000000;
  }
  if(aRes != ESsa){ 
    free(mRes);
    free(aRes);
  }
  if(ESStResonances() < 0){
    Fl	|=0x01;
  }
  M1	=N*qmn;
  M1	-=dM1;
  Mmx	=N*qmx+dM2;
  switch(Opt%10){
  case 1:
  case 3:
    M2	=N*qmn+dM2;
    break;
  default:
    M2	=Mmx;
    break;
  }
  Mp1	=M2-M1+1;
  ESStN		=N;
  ESStM1	=M1;
  ESStM2	=M2;

  if(Mgy < M1 || Mgy > M2){
    MaY	-=Mgy;
    Mgy	=(M1+M2)/2;
    MaY	+=Mgy;
  }
  if(MaY < M1 || MaY > M2){
    MaY	=(M1+M2)/2;
  }

  if(Fl == 0 && Opt){
    if(ESReInitMHDStabLowN() < 0){
      Fl	|=0x02;
    }
  }
  if(ESReInitMHDMetrTnsr()){
    Fl	|=0x04;
  }
  if(Fl == 0){
    ESMHDMetrTnsr();
  }
  if(Fl == 0){
    if(Opt && Opt < 5){
      if(ES2DStSolv(Opt,sa[1]) < 0){
	Fl	|=0x08;
      }
    }
    Opt	=0;
  }
  if(Fl && Fl < 0x08){
    printf("No capacity to perform stability test Fl=0x%2.2x%c\n",Fl,7);
    if((Fl&0x01)){
      nRes	=0;
      aRes	=ESsa;
      qRes	=ESsq;
    }
  }
  if((ESmemlev&0x08000000)){
    ESDeInitMHDStabLowN();
    ESmemlev -=0x08000000;
  }
  printf("t=%10.3e sec\n",((double)(clock()-ESst0))/CLOCKS_PER_SEC);
#ifdef TSC
  if(k0 < 0){
    k0	=0;
  }
  if(k0 > ES2Mp){
    k0	=ES2Mp;
  }
  if(k1 < k0+1){
    k1	=k0+1;
  }
  if(k1 > ES2Mp1){
    k1	=ES2Mp1;
  }
  if(i0 < 0){
    i0	=0;
  }
  if(i0 > ESNa-1){
    i0	=ESNa-1;
  }
  if(i1 < i0+2){
    i1	=i0+2;
  }
  if(i1 > ESNa1){
    i1	=ESNa1;
  }

  if(mV < M1){
    mV	=M1;
  }
  if(mV > M2){
    mV	=M2;
  }
  if(m1 < M1){
    m1	=M1;
  }
  if(m1 > M2-1){
    m1	=M2-1;
  }
  if(m2 < m1+1){
    m2	=m1+1;
  }
  if(m2 > M2){
    m2	=M2;
  }
  if(FlV < 0){
    FlV	=0;
  }
  if(FlV > 3){
    FlV	=3;
  }

  {
    int k,ki;
    double *p,*p2,s,t,s1,ds;

    if(Fc){
      if(k0 < 1){
	k0	=1;
      }
      switch(Fh){
      case 0:
	p	=ESGaas;
	p2	=ESGaas2;
	break;
      case 1:
	p	=ESGaqs;
	p2	=ESGaqs2;
	break;
      case 2:
	p	=ESGqqs;
	p2	=ESGqqs2;
	break;
      case 3:
	p	=ESGazs;
	p2	=ESGazs2;
	break;
      case 4:
	p	=ESGqzs;
	p2	=ESGqzs2;
	break;
      case 5:
	p	=ESGzzs;
	p2	=ESGzzs2;
	break;
      case 6:
	p	=ESsVs;
	p2	=ESsVs2;
	break;
      case 7:
	p	=ESVs;
	p2	=ESVs2;
	break;
      case 8:
	p	=ESLs;
	p2	=ESLs2;
	break;
      case 9:
	p	=ESDPs;
	p2	=ESDPs2;
	break;
      case 10:
	p	=ESsVs;
	p2	=ESsVs2;
	break;
      default :
	p	=ESsVs;
	p2	=ESsVs2;
	break;
      }
    }
    else{
      switch(Fh){
      case 0:
	p	=ESGaac;
	p2	=ESGaac2;
	break;
      case 1:
	p	=ESGaqc;
	p2	=ESGaqc2;
	break;
      case 2:
	p	=ESGqqc;
	p2	=ESGqqc2;
	break;
      case 3:
	p	=ESGazc;
	p2	=ESGazc2;
	break;
      case 4:
	p	=ESGqzc;
	p2	=ESGqzc2;
	break;
      case 5:
	p	=ESGzzc;
	p2	=ESGzzc2;
	break;
      case 6:
	p	=ESsVc;
	p2	=ESsVc2;
	break;
      case 7:
	p	=ESVc;
	p2	=ESVc2;
	break;
      case 8:
	p	=ESLc;
	p2	=ESLc2;
	break;
      case 9:
	p	=ESDPc;
	p2	=ESDPc2;
	break;
      case 10:
	p	=ESsVc;
	p2	=ESsVc2;
	break;
      default :
	p	=ESsVc;
	p2	=ESsVc2;
	break;
      }
    }
    Y[0]	=0.;
    Y[1]	=Y[0];
    y[0]	=0.;
    y[1]	=y[0];
    x[0]	=ESsa[i0];
    x[1]	=ESsa[i1-1];
    s		=(x[1]-x[0])/160.;
    t		=x[0];
    kW		=0;
    while(t <= x[1]){
      if(Fh < 10){
	ESSetSplA(t);
      }
      else{
	ESSetSplAA(t);
      }
      Wx[kW]	=t;
      for(k=k0; k < k1; k++){
	ki	=ESNa1*k;
	if(Fd == 1){
	  if(Fh < 10){
	    splRA(&s1,&ds,p+ki,p2+ki);
	  }
	  else{
	    if(k%2){
	      splRAA1(&s1,&ds,p+ki,p2+ki);
	    }
	    else{
	      splRAA0(&s1,&ds,p+ki,p2+ki);
	    }
	  }
	}
	else{
	  if(Fh < 10){
	    splRA(&ds,&s1,p+ki,p2+ki);
	  }
	  else{
	    if(k%2){
	      splRAA1(&ds,&s1,p+ki,p2+ki);
	    }
	    else{
	      splRAA0(&ds,&s1,p+ki,p2+ki);
	    }
	  }
	}
	if(y[0] > ds)
	  y[0]=ds;
	if(y[1] < ds)
	  y[1]=ds;
	if(k == k0){
	  if(Y[0] > ds)
	    Y[0]=ds;
	  if(Y[1] < ds)
	    Y[1]=ds;
	}
	W1[k][kW]	=ds;
      }
      t	+=s;
      kW++;
    }
    if(Fd == 0){
      for(i=i0; i < i1; i++){
	for(k=k0; k < k1; k++){
	  ki	=ESNa1*k;
	  if(Fh < 10 || k%2 == 0){
	    ds	=p[ki+i];
	  }
	  else{
	    ds	=ESsa[i]*p[ki+i];
	  }
	  if(y[0] > ds)
	    y[0]=ds;
	  if(y[1] < ds)
	    y[1]=ds;
	  if(k == k0){
	    if(Y[0] > ds)
	      Y[0]=ds;
	    if(Y[1] < ds)
	      Y[1]=ds;
	  }
	}
      }
    }
    Scale2D(5,2,x,y,2,6,2);
    for(k=k0; k < k1; k++){
      if(Fd == 0 && (Fh < 10 || k%2 == 0)){
	Plot2d(5,ESsa+i0,p+ESNa1*k+i0,i1-i0,6,0,0,0);
      }
      if(k%2)
	Plot2d(5,Wx,W1[k],kW,6,0,14,0);
      else
	Plot2d(5,Wx,W1[k],kW,6,0,4,0);
      if(k == k0)
	Plot2d(5,Wx,W1[k],kW,6,0,14,0);
    }
  }
  Wy	=W1[k0];
#endif
#endif/*Tbl19a*/
  if((ESmemlev&0x08000000)){
    ESDeInitMHDStabLowN();
    ESmemlev -=0x08000000;
  }
  if((ESmemlev&0x04000000)){
    ESDeInitMHDMetrTnsr();
    ESmemlev -=0x04000000;
  }
  if(aRes != ESsa){ 
    free(mRes);
    free(aRes);
  }
  ESStNRec	=NRec;
  ESStMRec	=MRec;
  return(0);
}
#endif /*XWIN*/

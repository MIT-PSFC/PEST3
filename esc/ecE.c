#include "esc_local.h"
#define aaaDEBUG
#define RKS
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

extern clock_t clock(void);
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern int ESNa,ESNa1,ESNp,ESNp1,ESnAP;
extern int ESMp,ESMp1,ES2Mp,ES2Mp1,ESnMp;
extern int ESFp,ESFp1,ESnAF;
extern double ESa0;
extern double *ESsa,*ESpa,*rT1a,*rT1gt,*zT1a,*zT1gt;
extern double *ECb,*ECb1a,*ECb2a;
extern double *EScs1,*ESsn1;

extern double R0,Z0;
extern double *ESaR0,*ESaZ0,*ESsb,*ESsb1a,*ESsb2a;

extern int ESEqSolvFl,ESEqSolvRadC,ESEqSolvInPr,ESEqSolvInCr;
extern double ESEqSolvTol,ESEqTol;
extern int ESFail,ESEqSolvIt,ESEqIt;
/*
   |     nnnn  DDDD dddd  AAAA aaaa  eeee EZpppp
   |____|____||____|____||____|____||____|____|
   p - number of successive failures;
   e - lack of accuracy in the middle;
   a - lack of accuracy near the axis;
   A - lack of accuracy near the edge;
   d - overshooting (negative D) near the axis;
   D - overshooting (negative D) near the edge;
   n - overshootings (negative D) in the middle;
   */

extern double *ESg22c,*ESg22s,*ESg22c1,*ESg22s1,*ESg22c2,*ESg22s2;
extern double *ESr22,*ESr222,*ESr22c,*ESr22s,*ESr22c2,*ESr22s2;
extern double *ESg12c,*ESg12s,*ESg11c,*ESg11s;
extern double *ESg12c2,*ESg12s2,*ESg11c2,*ESg11s2;
extern double *ESLc,*ESLs,*ESVc,*ESVs,*ESDPc,*ESDPs;
extern double *ESLc1,*ESLs1,*ESVc1,*ESVs1,*ESDPc1,*ESDPs1;
extern double *ESLc2,*ESLs2,*ESVc2,*ESVs2,*ESDPc2,*ESDPs2;

extern double *rcT,*rcT1a,*rcT2a,*rsT,*rsT1a,*rsT2a;
extern double *dR0T,*dZ0T,*dbT,*dbT1a,*dbT2a;
extern double *drcT,*drcT1a,*drcT2a,*drsT,*drsT1a,*drsT2a;

extern double ESRBt,ESRext;
extern double ESgbext,ESgbN;
extern double *ESgY0,*ESgY02a,*ESaY,ESgY2a1a;
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
,*ESgF,*ESgF1a,*ESgF2a

,*ESaF,*ESaF1a,*ESaF2a
,*ESFF,*ESFF1a,*ESFF2a
,*ESdgY,*ESdgY1a,*ESdgY2a
,*ESjR,*ESjR1a,*ESjR2a;

extern double *EZrcs,*EZrsn,*EZd1rcs,*EZd1rsn,*EZd2rcs,*EZd2rsn,*EZz0,*EZd1z0,*EZd2z0;
extern double *rcE,*rsE,*rcE1,*rsE1,*rcE2,*rsE2,*bE,*bE1,*bE2,R0E,Z0E;

extern double eps_rel,eps_abs;
extern double Epsrvb,Eps_noise;

extern double*Fcs,*Fsn,*EZxinf,*EZyinf,*d1yinf,*d2yinf,*EZxgt,*EZygt;

extern int M0Mm;
extern int ESiAx;
extern int ECEqIterFl;

static double*gpr,*gpi,*Yr,*Yi,*dgpr,*dgpi,*dYr,*dYi,*Xr,*Xi,*Ur,*Ui,*dUr,*dUi;
static double*gpr0,*gpi0,*Yr0,*Yi0,*dgpr0,*dgpi0,*dYr0,*dYi0;
static double*gpr1,*gpi1,*Yr1,*Yi1;

static int Mr;
static double *drc,*drs,*dz0,*d1dz0;

extern double *EZgper,*EZgpei,*EZdgper,*EZdgpei,*aYer,*aYei,*daYer,*daYei;

static double *LHSr,*LHSi,*RHSr,*RHSi;

extern double *rCc,*rCs,*rSc,*rSs,*drCc,*drCs,*drSc,*drSs;
extern double *bCc,*bSc,*dbCc,*dbSc;

static double *daYr,*daYi;

static double *wMem;
static double *EZd1drc,*EZd1drs;
static double*wr22c,*wr22s,*wg11c,*wg11s,*wg12c,*wg12s,*wg22c,*wg22s;
static double *wuc,*wus,*wduc,*wdus,*wdKc,*wdKs,*wvc,*wvs;
static double *gdjsr,*gdjsi,*gdjsr1,*gdjsi1,gdx,gdx1;

static double rR0,*ac,*Kr,*Ki;
static int *ind;

static int Isize0=0,Isize1=0,Isize2=0,Isize3=0;
static int NIter=0,kFail=0,Fail=0;
static double EcErr=0.;

static int EcCholFl=0;

#ifdef RKS
static void (*peq2D)(int *neq,double *x,double *gyr,double *dgyr);
#endif

#ifndef mzl_2DEq
int ECReInitFgY()
{
  int k,m,M0L0;
  k	=2*ES2Mp+1;
  m	=16*k+4*ESMp1;
  if(Isize0 < m){
    if(Isize0){
      free(ind);
      free(ac);
      free(wMem);
    }
    Isize0	=m;
    wMem	=(double*)malloc(Isize0*sizeof(double));
    if(wMem == NULL){
      printf("Failure in memory allocation for wMem in mem4vgp\n");
      exit(0);
    }
    ac	=(double*)malloc(2*ES2Mp1*ES2Mp1*sizeof(double));
    Kr	=ac;
    Ki	=Kr+ES2Mp1*ES2Mp1;
    ind	=(int*)malloc(ES2Mp1*sizeof(int));
  }
  wr22c	=wMem+ES2Mp;
  wr22s	=wr22c	+k;
  wg11c	=wr22s	+k;
  wg11s	=wg11c	+k;
  wg12c	=wg11s	+k;
  wg12s	=wg12c	+k;
  wg22c	=wg12s	+k;
  wg22s	=wg22c	+k;
  wuc	=wg22s	+k;
  wus	=wuc	+k;
  wduc	=wus	+k;
  wdus	=wduc	+k;
  wdKc	=wdus	+k;
  wdKs	=wdKc	+k;
  wvc	=wdKs	+k;
  wvs	=wvc	+k;
  gdjsr	=wvs	+ES2Mp1;
  gdjsi	=gdjsr	+ESMp1;
  gdjsr1	=gdjsi	+ESMp1;
  gdjsi1	=gdjsr1	+ESMp1;
  for(m=0; m < ESMp1; m++){
    gdjsr[m]	=0.;
    gdjsi[m]	=0.;
    gdjsr1[m]	=0.;
    gdjsi1[m]	=0.;
  }
  M0Mm	=ESMp1*ES2Mp1;
  m	=8*M0Mm+4+12*ES2Mp1;
  if(Isize1 < m){
    if(Isize1) free(gpr);
    Isize1	=m;
    gpr		=(double*)malloc(Isize1*sizeof(double));
    if(gpr == NULL){
      printf("Failure in memory allocation for gpr in ECReInitFgY\n");
      exit(0);
    }
  }
  gpi	=gpr	+M0Mm;
  Yr	=gpi	+M0Mm;
  Yi	=Yr	+M0Mm;
  gpr0	=Yi	+M0Mm	+2;
  gpi0	=gpr0	+ES2Mp1;
  Yr0	=gpi0	+ES2Mp1;
  Yi0	=Yr0	+ES2Mp1;
  dgpr	=Yi0	+ES2Mp1;
  dgpi	=dgpr	+M0Mm;
  dYr	=dgpi	+M0Mm;
  dYi	=dYr	+M0Mm;
  dgpr0	=dYi	+M0Mm	+2;
  dgpi0	=dgpr0	+ES2Mp1;
  dYr0	=dgpi0	+ES2Mp1;
  dYi0	=dYr0	+ES2Mp1;
  gpr1	=dYi0	+ES2Mp1;
  gpi1	=gpr1	+ES2Mp1;
  Yr1	=gpi1	+ES2Mp1;
  Yi1	=Yr1	+ES2Mp1;
  m	=M0Mm;
  if(Isize2 < 6*m){
    if(Isize2){
      free(Xr);
    }
    Isize2= 6*m;
    Xr= (double*)malloc(Isize2*sizeof(double));
    if(Xr==NULL){
      printf("Failure in memory allocation for Xr in ECReInitFgY\n");
      exit(0);
    }
    Xi	=Xr+m;
    Ur	=Xi+m;
    Ui	=Ur+m;
    dUr	=Ui+m;
    dUi	=dUr+m;
  }

  m	=ESNa1*M0Mm;
  Mr	=ESMp1;
  M0L0	=ESMp1*Mr;
  if(Isize3 < 8*m+8*M0L0+4*ESMp1+(2+4*Mr)*ESNa1){
    if(Isize3){
      free(EZgper);
    }
    Isize3	=4*m+8*m+8*M0L0+4*ESMp1+(2+4*Mr)*ESNa1;
    EZgper	=(double*)malloc(Isize3*sizeof(double));
    if(EZgper == NULL){
      printf("Failure in memory allocation for EZgper in ECReInitFgY\n");
      exit(0);
    }
  }
  EZgpei	=EZgper	+m;
  EZdgper	=EZgpei	+m;
  EZdgpei	=EZdgper	+m;
  aYer	=EZdgpei	+m;
  aYei	=aYer	+m;
  daYer	=aYei	+m;
  daYei	=daYer	+m;

  LHSr	=daYei	+m;
  LHSi	=LHSr	+m;
  RHSr	=LHSi	+m;
  RHSi	=RHSr	+m;

  drc	=RHSi	+m;
  drs	=drc	+Mr*ESNa1;
  EZd1drc	=drs	+Mr*ESNa1;
  EZd1drs	=EZd1drc	+Mr*ESNa1;
  rCc	=EZd1drs	+Mr*ESNa1;
  rCs	=rCc	+M0L0;
  rSc	=rCs	+M0L0;
  rSs	=rSc	+M0L0;
  drCc	=rSs	+M0L0;
  drCs	=drCc	+M0L0;
  drSc	=drCs	+M0L0;
  drSs	=drSc	+M0L0;

  bCc	=drSs	+M0L0;
  bSc	=bCc	+ESMp1;
  dbCc	=bSc	+ESMp1;
  dbSc	=dbCc	+ESMp1;
  dz0	=dbSc	+ESMp1;
  d1dz0	=dz0	+ESNa1;
 return(0);
}

int ECDeInitFgY()
{
  if(Isize3){
    free(EZgper);
    Isize3=0;
  }
  if(Isize2){
    free(Xr);
    Isize2=0;
  }
  if(Isize1){
    free(gpr);
    Isize1=0;
  }
  if(Isize0){
    free(ind);
    free(ac);
    free(wr22c-ES2Mp);
    Isize0=0;
  }
  return(0);
}
#ifndef stg_2DInit
int Start2D0(double x)
{
  return(0);
}

int Start2DAAA(double x)
{
  int i,j,k,m;
  double s,sr,si,f,d,am,sj;
  double *pr,*pi,*dpr,*dpi,*xr,*xi,*yr,*yi,*lr,*li;
  double r1r,r1i,z1i;
  double x1r,x1i,X1r,X1i;
  double t1r,t1i,T1r,T1i;
  double tr,ti;

  /*Getting metric tensor{*/
  ESSetSplA(x);
  splRA(wg12c,NULL,ESg12c,ESg12c2);
  wg12s[0]	=0.;
  splRA(wg22c,NULL,ESg22c,ESg22c2);
  wg22c[0]	*=x;
  wg22s[0]	=0.;
  k	=0;
  for(i=1; i < ES2Mp1; i++){
    k	+=ESNa1;
    splRA(wg12c+i,NULL,ESg12c+k,ESg12c2+k);
    wg12c[-i]	=wg12c[i];
    splRA(wg22c+i,NULL,ESg22c+k,ESg22c2+k);
    wg22c[i]	*=x;
    wg22c[-i]	=wg22c[i];
    splRA(wg12s-i,NULL,ESg12s+k,ESg12s2+k);
    wg12s[i]	=-wg12s[-i];
    splRA(wg22s-i,NULL,ESg22s+k,ESg22s2+k);
    wg22s[-i]	*=x;
    wg22s[i]	=-wg22s[-i];
  }
  sj	=ESjs[0];
  z1i	=-EZcr2*ESsb1a[0];
  k	=ESNa1;
  r1r	= rcT1a[k];
  r1i	=-rsT1a[k];
  wuc[1]	=0.;
  wus[1]	=0.;
  wuc[-1]	=0.;
  wus[-1]	=0.;
  s		=EZcr2*x;
  wuc[0]	=s*rcT2a[0];
  wus[0]	=0.;
  for(i=2; i < ES2Mp1; i++){
    k	+=ESNa1;
    wuc[i]	= s*rcT2a[k];
    wus[i]	=-s*rsT2a[k];
    wuc[-i]	= wuc[i];
    wus[-i]	=-wus[i];
  }
  x1r	= r1r-z1i;
  x1i	= r1i;
  X1r	= r1r+z1i;
  X1i	=-r1i;
  t1r	= r1r+z1i;
  t1i	= r1i;
  T1r	= r1r-z1i;
  T1i	=-r1i;
  /*}*/
  for(k=0; k < M0Mm; k++){
    gpr[k]	=0.;
    gpi[k]	=0.;
    dgpr[k]	=0.;
    dgpi[k]	=0.;
    Yr[k]	=0.;
    Yi[k]	=0.;
  }
  pr	=gpr	+ESMp;
  pi	=gpi	+ESMp;
  dpr	=dgpr	+ESMp;
  dpi	=dgpi	+ESMp;
  yr	=Yr	+ESMp;
  yi	=Yi	+ESMp;
  lr	=dYr	+ESMp;
  li	=dYi	+ESMp;

  /*m=0 polynomial solution{*/
  d	=ESaR0[0];
  f	=-EZcr4*x*sj;
  s	=f*(X1r*t1r-X1i*t1i+x1r*T1r-x1i*T1i);
  sr	=f*(x1r*t1r-x1i*t1i);
  si	=f*(x1r*t1i+x1i*t1r);
  
  pr[0]		=x*d*s;
  dpr[0]	=2.*d*s;
  f		=EZcr4*x;
  s		*=f;
  pr[1]		=x*s*r1r;
  dpr[1]	=3.*s*r1r;
  pi[1]		=x*s*r1i;
  dpi[1]	=3.*s*r1i;
  
  pr[2]		=x*d*sr;
  pi[2]		=x*d*si;
  dpr[2]	=2.*d*sr;
  dpi[2]	=2.*d*si;
  sr		*=f;
  si		*=f;
  if(3 < ESMp1){
    s	=sr*r1r-si*r1i;
    pr[3]	+=x*s;
    dpr[3]	+=3.*s;
    s	=si*r1r+sr*r1i;
    pi[3]	+=x*s;
    dpi[3]	+=3.*s;
  }
  s	=sr*r1r+si*r1i;
  pr[1]	+=x*s;
  dpr[1]+=3.*s;
  s	=si*r1r-sr*r1i;
  pi[1]	+=x*s;
  dpi[1]+=3.*s;

  f	=-0.5*x*x*d*sj;
  d	*=-1.5*x*sj;
  for(i=0; i < ESMp1; i++){
    j	=i-1;
    s	=r1r*wuc[j]-r1i*wus[j];
    pr[i]	+=f*s;
    dpr[i]	+=d*s;
    s	=sr*wus[j]+si*wuc[j];
    pi[i]	+=f*s;
    dpi[i]	+=d*s;
    j	=i+1;
    s	=r1r*wuc[j]+r1i*wus[j];
    pr[i]	+=f*s;
    dpr[i]	+=d*s;
    s	=sr*wus[j]-si*wuc[j];
    pi[i]	+=f*s;
    dpi[i]	+=d*s;
  }
  wus[0]=EZcr2*x*rsT2a[0];
  s	=wus[0]*z1i;
  pi[1]	+=f*s;
  dpi[1]+=d*s;
  for(i=1; i < ESMp1; i++){
    pr[-i]	= pr[i];
    pi[-i]	=-pi[i];
    dpr[-i]	= dpr[i];
    dpi[-i]	=-dpi[i];
  }
  /*}*/

  s	=EZcr4/ESaR0[0];
  t1r	*=s;
  t1i	*=s;
  T1r	*=s;
  T1i	*=s;
  yr[0]	= 1.;
  yi[0]	= 0.;
  am	=1.;
  for(m=1; m < ESMp1; m++){/*$m\neq1$ solutions*/
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    am	*=x;
    for(k=-m+1; k < m; k++){
      sr	=yr[k-ES2Mp1];
      si	=yi[k-ES2Mp1];
      /*Second correction O(a);*/
      f	=m*x;
      d	=m*(m+1);
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	s	=sr*wuc[j]-si*wus[j];
	pr[i]	+=f*s;
	dpr[i]	+=d*s;
	s	=sr*wus[j]+si*wuc[j];
	pi[i]	+=f*s;
	dpi[i]	+=d*s;
      }
      /* Main order contribution;*/
      i	=k+1;
      s		=x1r*sr-x1i*si;
      yr[i]	=x*s;
      pr[i]	+=yr[i];
      dpr[i]	+=m*s;
      s		=x1r*si+x1i*sr;
      yi[i]	=x*s;
      pi[i]	+=yi[i];
      dpi[i]	+=m*s;
      i	=k-1;
      s		=X1r*sr-X1i*si;
      yr[i]	=x*s;
      pr[i]	+=yr[i];
      dpr[i]	+=m*s;
      s		=X1r*si+X1i*sr;
      yi[i]	=x*s;
      pi[i]	+=yi[i];
      dpi[i]	+=m*s;
    }
    f	=1./m;
    d	=(m+1)*f;
    f	*=x;
    /*First correction $O(a)$*/
    for(k=-m; k <= m; k++){
      sr	=yr[k];
      si	=yi[k];
      i		=k+1;
      if(i < ESMp1){
	s	=t1r*sr-t1i*si;
	pr[i]	+=f*s;
	dpr[i]	+=d*s;
	s	=t1r*si+t1i*sr;
	pi[i]	+=f*s;
	dpi[i]	+=d*s;
      }
      i	=k-1;
      if(i > -ESMp1){
	s	=T1r*sr-T1i*si;
	pr[i]	+=f*s;
	dpr[i]	+=d*s;
	s	=T1r*si+T1i*sr;
	pi[i]	+=f*s;
	dpi[i]	+=d*s;
      }	
    }

    /*Normalization;$*/
    sr	=pr[m];
    si	=pi[m];
    s	=1./(sr*sr+si*si);
    sr	*=s;
    s	*=si;
    si	=pi[-m]*sr+pr[-m]*s;
    sr	=pr[-m]*sr-pi[-m]*s;
    for(i=-ESMp; i < ESMp1; i++){
      lr[i]	= pr[-i];
      li[i]	=-pi[-i];
    }
    for(i=-ESMp; i < ESMp1; i++){
      pr[i]	-=sr*lr[i]-si*li[i];
      pi[i]	-=sr*li[i]+si*lr[i];
    }
    for(i=-ESMp; i < ESMp1; i++){
      lr[i]	= dpr[-i];
      li[i]	=-dpi[-i];
    }
    for(i=-ESMp; i < ESMp1; i++){
      dpr[i]	-=sr*lr[i]-si*li[i];
      dpi[i]	-=sr*li[i]+si*lr[i];
    }
    if(m == 1){
      sr	=dpr[m];
      si	=dpi[m];
      s	=1./(sr*sr+si*si);
    }
    else{
      sr	=pr[m];
      si	=pi[m];
      s	=am/(sr*sr+si*si);
    }
    sr	*=s;
    si	*=-s;
    for(i=-ESMp; i < ESMp1; i++){
      s	=pr[i];
      pr[i]	=sr*s-si*pi[i];
      pi[i]	=si*s+sr*pi[i];
      s	=dpr[i];
      dpr[i]	=sr*s-si*dpi[i];
      dpi[i]	=si*s+sr*dpi[i];
    }
  }

  /*Normalization of $m=0$;*/
  pr	=gpr	+ESMp;
  pi	=gpi	+ESMp;
  dpr	=dgpr	+ESMp;
  dpi	=dgpi	+ESMp;
  i	=ESMp+2*ES2Mp1;
  yr	=gpr	+i;
  yi	=gpi	+i;
  lr	=dgpr	+i;
  li	=dgpi	+i;
  sr	=yr[2];
  si	=yi[2];
  s	=1./(sr*sr+si*si);
  sr	*=s;	/*$\R{1}{\gy^{*2}_2}$*/
  s	*=si;
  tr	=pr[-2]*sr-pi[-2]*s;	/*$\R{\gy^0_{-2}}{\gy^{*2}_2}$*/
  ti	=pi[-2]*sr+pr[-2]*s;
  si	=pi[2]*sr-pr[2]*s;
  sr	=pr[2]*sr+pi[2]*s;	/*$\R{\gy^0_2}{\gy^2_2}$*/
  for(i=-ESMp; i < ESMp1; i++){
    pr[i]	-=sr*yr[i]-si*yi[i]+tr*yr[-i]+ti*yi[-i];
    pi[i]	-=si*yr[i]+sr*yi[i]+ti*yr[-i]-tr*yi[-i];
    dpr[i]	-=sr*lr[i]-si*li[i]+tr*lr[-i]+ti*li[-i];
    dpi[i]	-=si*lr[i]+sr*li[i]+ti*lr[-i]-tr*li[-i];
  }
  s	=x*x;
  s	*=0.125*s*ESjp[0]*ESVc2[0]/wg22c[0];
  dpr[0]	-=s;
  pr[0]		-=EZcr4*s;

  yr	=Yr	+ESMp;
  yi	=Yi	+ESMp;
  for(m=0; m < ESMp1; m++){/*All $Y_m$*/
    for(i=-ESMp; i < ESMp; i++){
      yr[i]	=0.;
      yi[i]	=0.;
    }
    for(k=-ESMp; k < ESMp1; k++){
      sr	=dpr[k];
      si	=dpi[k];
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	yr[i]	+=sr*wg22c[j]-si*wg22s[j];
	yi[i]	+=si*wg22c[j]+sr*wg22s[j];
      }
      sr	= k*pi[k];
      si	=-k*pr[k];
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	yr[i]	+=sr*wg12c[j]-si*wg12s[j];
	yi[i]	+=si*wg12c[j]+sr*wg12s[j];
      }
    }
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
  }
  /*Averaged solution;*/
  ESaY[0]	=-EZcr2*sj*ESLc[0];
  ESdgY[0]	=ESaY[0]*x/wg22c[0];
  ESgY02a[0]	=ESdgY[0];
  s		=x*x;
  i		=M0Mm;
  Yi[i]		=s*(ESaY[0]-0.125*s*ESjp[0]*ESVc2[0]); /*$Y_0$*/
  i++;
  Yi[i]		=EZcr2*s*ESdgY[0]; /*$bgY_0$*/
  Yi[i]		=-EZcr4*s*x*(sj*ESLc[0]+0.125*s*ESjp[0]*ESVc2[0])/wg22c[0];
  gpr[ESMp]	=Yi[i];
  s	=2./s;
  k	=ESMp;
  return(0);
}

int Start2Dgm(double x)
{
  int i,j,k,m,mm;
  int M,M1;
  double Bf[1024];
  double *bxr,*bxi,x1r,x1i,tyr,by0;
  double *Jr,*Ji;
  double *bKr,*bKi;
  double *bMr,*bMi;
  double *bNr,*bNi;
  double *bDr,*bDi,tD;
  double s,sr,si,tr,ti;
  double K0,Kr,Ki;
  double Cr,Ci;
  double *ur,*ui,*EZvr,*vi;
  double *pr,*pi,*dpr,*dpi;
  double *yr,*yi,*Rr,*Ri;
  double am,a1r,a1i,a2r,a2i;
  
  bxr	=Bf	+ESFp;
  k	=2*ESFp+1;
  bxi	=bxr	+k;
  bKr	=bxi	+k;
  k	=2*ES2Mp+1;
  bKi	=bKr	+k;
  bMr	=bKi	+k;
  bMi	=bMr	+k;
  bNr	=bMi	+k;
  bNi	=bNr	+k;
  bDr	=bNi	+k+2;
  bDi	=bDr	+k+4;
  Jr	=bDi	+k+4;
  Ji	=Jr	+ESMp1;
  EZvr	=Ji	+ES2Mp1;
  vi	=EZvr	+ES2Mp1;
  ur	=vi	+ES2Mp1;
  ui	=ur	+ES2Mp1;
  Rr	=ui	+ES2Mp1;
  Ri	=Rr	+ES2Mp1;

  /*Getting geometry near the origin{*/
  s	=EZcr2*x;
  tyr	=EZcr2*ESsb1a[0];
  bxr[0]=s*rcT2a[0];
  by0	=s*rsT2a[0];
  bxi[0]=0.;
  bxr[1]=0.;
  bxi[1]=0.;
  bxr[-1]=0.;
  bxi[-1]=0.;
  k	=ESNa1;
  x1r	= rcT1a[k];
  x1i	=-rsT1a[k];
  M	=2;
  for(i=2; i < ESFp1; i++){
    k		+=ESNa1;
    bxr[i]	= s*rcT2a[k];
    bxi[i]	=-s*rsT2a[k];
    bxr[i]	=0.;
    bxi[i]	=0.;
    bxr[-i]	= bxr[i];
    bxi[-i]	=-bxi[i];
    if(bxr[i] != 0. || bxi[i] != 0.) M=i;
  }

  bxr[0]=0.;
  by0	=0.;


  M1	=M+1;
  tD	=4.*x1r*tyr;
  s	=1./tD;
  sr	=x1r*x1r;
  si	=x1i*x1i;
  tr	=tyr*tyr;
  K0	=2.*(tr+sr+si)*s;
  Kr	=(tr+si-sr)*s;
  Ki	=2.*x1r*x1i*s;
  /*}*/

  /*Getting corrections to metric tensor{*/
  j	=ES2Mp1+2;
  for(i=0; i < j; i++){
    k	=i-1;
    s	=(2.-k)*tyr;
    bDr[i]	=s*bxr[k];
    bDi[i]	=s*bxi[k];
    k	=i+1;
    if(k < M1){
      s	=(2.+k)*tyr;
      bDr[i]	+=s*bxr[k];
      bDi[i]	+=s*bxi[k];
    }
  }
  bDi[0]	=0.;
  bDr[1]	+=2.*x1i*by0;
  bDi[1]	-=2.*x1r*by0;
  bDr[-1]	= bDr[1];
  bDi[-1]	=-bDi[1];
  bDr[-2]	= bDr[2];
  bDi[-2]	=-bDi[2];
  for(i=0; i < ES2Mp1; i++){
    bKr[i]	=bDr[i]*K0;
    bKi[i]	=bDi[i]*K0;
    bNr[i]	=bDr[i]*K0;
    bNi[i]	=bDi[i]*K0;
    k	=i-2;
    sr	=bDr[k];
    si	=bDi[k];
    tr	=sr*Kr-si*Ki;
    ti	=si*Kr+sr*Ki;
    bKr[i]	-=tr;
    bKi[i]	-=ti;
    bNr[i]	+=tr;
    bNi[i]	+=ti;
    bMr[i]	=-ti;
    bMi[i]	= tr;
    k	=i+2;
    sr	=bDr[k];
    si	=bDi[k];
    tr	=sr*Kr+si*Ki;
    ti	=si*Kr-sr*Ki;
    bKr[i]	-=tr;
    bKi[i]	-=ti;
    bNr[i]	+=tr;
    bNi[i]	+=ti;
    bMr[i]	+=ti;
    bMi[i]	-=tr;
    k	=i-1;
    s	=-2.*k;
    sr	=s*bxr[k];
    si	=s*bxi[k];
    bKr[i]	+=sr*x1r-si*x1i;
    bKi[i]	+=si*x1r+sr*x1i;
    s	=k+2.;
    sr	=-s*bxi[k];
    si	= s*bxr[k];
    bMr[i]	+=sr*x1r-si*x1i;
    bMi[i]	+=si*x1r+sr*x1i;
    sr	=4.*bxr[k];
    si	=4.*bxi[k];
    bNr[i]	+=sr*x1r-si*x1i;
    bNi[i]	+=si*x1r+sr*x1i;
    k	=i+1;
    if(k < M1){
      bNr[i]	+=sr*x1r+si*x1i;
      bNi[i]	+=si*x1r-sr*x1i;
      s	=2.*k;
      sr	=s*bxr[k];
      si	=s*bxi[k];
      bKr[i]	+=sr*x1r+si*x1i;
      bKi[i]	+=si*x1r-sr*x1i;
      s	=k-2.;
      sr	=-s*bxi[k];
      si	= s*bxr[k];
      bMr[i]	+=sr*x1r+si*x1i;
      bMi[i]	+=si*x1r-sr*x1i;
    }
  }
  bMr[1]	+=2.*tyr*by0;
  bNi[1]	-=4.*tyr*by0;
  s	=1./tD;

  for(k=0; k < ES2Mp1; k++){
    bKr[k]	*=s;
    bMr[k]	*=s;
    bNr[k]	*=s;
    if(k){
      bKi[k]	*=s;
      bMi[k]	*=s;
      bNi[k]	*=s;
      i	=-k;
      bDr[i]	= bDr[k];
      bDi[i]	=-bDi[k];
      bKr[i]	= bKr[k];
      bKi[i]	=-bKi[k];
      bMr[i]	= bMr[k];
      bMi[i]	=-bMi[k];
      bNr[i]	= bNr[k];
      bNi[i]	=-bNi[k];
    }
  }
  bKi[0]	=0.;
  bMi[0]	=0.;
  bNi[0]	=0.;
  /*}*/


  for(k=0; k < M0Mm; k++){
    gpr[k]	=0.;
    gpi[k]	=0.;
    dgpr[k]	=0.;
    dgpi[k]	=0.;
    Yr[k]	=0.;
    Yi[k]	=0.;
    dYr[k]	=0.;
    dYi[k]	=0.;
  }
  /*Main order solution{*/
  pr	=gpr+ESMp;
  pi	=gpi+ESMp;
  dpr	=dgpr+ESMp;
  dpi	=dgpi+ESMp;
  yr	=Yr+ESMp;
  yi	=Yi+ESMp;
  Jr[0]	=2.*ESgm[0]*ESaF[0]*K0/(ESaR0[0]*ESaR0[0]);
  Ji[0]	=0.;
  yr[0]	=-EZcr4*ESaR0[0]*tD*Jr[0]/K0;
  pr[0]	=yr[0];
  dpr[0]=2.*yr[0];
  yr	+=ES2Mp1;
  yi	+=ES2Mp1;
  pr	+=ES2Mp1;
  pi	+=ES2Mp1;
  dpr	+=ES2Mp1;
  dpi	+=ES2Mp1;
  yr[1]	=1.;
  pr[1]	=yr[1];
  dpr[1]=yr[1];
  Jr[1]	=0.;
  Ji[1]	=0.;
  yr	+=ES2Mp1;
  yi	+=ES2Mp1;
  pr	+=ES2Mp1;
  pi	+=ES2Mp1;
  dpr	+=ES2Mp1;
  dpi	+=ES2Mp1;
  yr[2]	=1.;
  pr[2]	=yr[2];
  dpr[2]=2.*yr[2];
  s	=-8./(ESaR0[0]*tD);
  Jr[2]	= s*Kr;
  Ji[2]	=-s*Ki;
  for(m=3; m < ESMp1; m++){
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    i	=-m;
    EZvr[i]	=0.;
    vi[i]	=0.;
    for(k=i+2; k < m; k+=2,i+=2){
      if(k != 2){
	s=(m-k)*(m-k+2);
	Cr=s*Kr;
	Ci=s*Ki;
	s	=k ? (m*m-k*k)*K0 : ESaR0[0]*tD;
	sr=s+Cr*EZvr[i]-Ci*vi[i];
	si=Ci*EZvr[i]+Cr*vi[i];
	s	=1./(sr*sr+si*si);
	sr*= s;
	si*=-s;
	tr=Cr*yr[i]-Ci*yi[i];
	ti=Ci*yr[i]+Cr*yi[i];
	yr[k]	=-tr*sr+ti*si;
	yi[k]	=-ti*sr-tr*si;
      }
      else{
	s	=(m*m-k*k)*K0;
	sr=1./s;
	si=0.;
	yr[k]	=0.;
	yi[k]	=0.;
      }
      if(k+2){
	s	=(m+k)*(m+k+2);
	tr= s*Kr;
	ti=-s*Ki;
	EZvr[k]	=-tr*sr+ti*si;
	vi[k]	=-ti*sr-tr*si;
      }
      else{
	EZvr[k]	=0.;
	vi[k]	=0.;
      }
    }
    i	=m;
    yr[i]	=1.;
    yi[i]	=0.;
    j	=2-m;
    while(i > j){
      k	=i-2;
      yr[k]	+=EZvr[k]*yr[i]-vi[k]*yi[i];
      yi[k]	+=vi[k]*yr[i]+EZvr[k]*yi[i];
      i	-=2;
    }
    Jr[m]	=yr[0];
    Ji[m]	=yi[0];
    yr[0]	=0.;
    yi[0]	=0.;
    for(i=-m+2; i <= m; i++){
      pr[i]	=yr[i];
      pi[i]	=yi[i];
      dpr[i]	=m*yr[i];
      dpi[i]	=m*yi[i];
    }
  }  
  /*}*/
  /*Corrections to the EZmain order{*/
  pr	=gpr+ESMp;
  pi	=gpi+ESMp;
  dpr	=dgpr+ESMp;
  dpi	=dgpi+ESMp;
  yr	=Yr+ESMp;
  yi	=Yi+ESMp;
  am	=1.;
  for(m=0; m < ESMp1; m++){/*All fundamental vectors*/
    am	*=x;
    if(m == 1) am	=1.;
    /*RHS;*/
    for(i=-ESMp; i < ESMp1; i++){
      Rr[i]	=-ESaR0[0]*(bDr[i]*Jr[m]-bDi[i]*Ji[m]);
      Ri[i]	=-ESaR0[0]*(bDi[i]*Jr[m]+bDr[i]*Ji[m]);
      for(j=2-m; j <= m; j+=2){
	k	=i-j;
	s	=j*(m+1)+m*i;
	sr	=-m*(m+1)*bKr[k]-s*bMi[k]+i*j*bNr[k];
	si	=-m*(m+1)*bKi[k]+s*bMr[k]+i*j*bNi[k];
	Rr[i]	+=sr*pr[j]-si*pi[j];
	Ri[i]	+=si*pr[j]+sr*pi[j];
      }
      j	=i-1;
      if(0 && 2-m <= j && j <= m){
	sr=x*(m-j)*tyr/ESaR0[0];
	Rr[i]	+=sr*pr[j];
	Ri[i]	+=sr*pi[j];
      }
      j	=i+1;
      if(0 && 2-m <= j && j <= m){
	sr=x*(m+j)*tyr/ESaR0[0];
	Rr[i]	+=sr*pr[j];
	Ri[i]	+=sr*pi[j];
      }
    }
    mm	=m ? m+1 : 3;
    /*Elimination of resonant harmonic in the RHS;*/
    if(mm < ESMp1){ 
      if(m){
	ur	=yr+ES2Mp1;
	ui	=yi+ES2Mp1;
      }
      else{
	ur	=yr+3*ES2Mp1;
	ui	=yi+3*ES2Mp1;
      }
      k	=mm-2;
      sr	=ur[k];
      si	=ui[k];
      k		+=2;
      tr	=mm*K0*ur[k]+Kr*sr-Ki*si;
      ti	=mm*K0*ui[k]+Ki*sr+Kr*si;
      k	=2-mm;
      sr	=Kr*ur[k]+Ki*ui[k];
      si	=Ki*ur[k]-Kr*ui[k];
      s		=EZcr2/(tr*tr+ti*ti-sr*sr-si*si);
      sr	*=-s;
      si	*=-s;
      tr	*=s;
      ti	*=s;
      Cr	=Rr[-mm];
      Ci	=Ri[-mm];
      a1r	=Rr[mm]*tr+Ri[mm]*ti+Cr*sr+Ci*si;
      a1i	=Ri[mm]*tr-Rr[mm]*ti+Ci*sr-Cr*si;
      a2r	=Rr[mm]*tr-Ri[mm]*ti+Cr*sr-Ci*si;
      a2i	=Ri[mm]*tr+Rr[mm]*ti+Ci*sr+Cr*si;
      s	=log(x/ESsa[1]);
      for(i=-mm; i <= mm; i+=2){
	sr	=a1r*ur[i]-a1i*ui[i]+a2r*ur[-i]+a2i*ui[-i];
	si	=a1i*ur[i]+a1r*ui[i]+a2i*ur[-i]-a2r*ui[-i];
	dpr[i]	-=(mm*s+1.)*sr;
	dpi[i]	-=(mm*s+1.)*si;
	pr[i]	-=s*sr;
	pi[i]	-=s*si;
	tr	=2.*mm*K0;
	Rr[i]	+=tr*sr;
	Ri[i]	+=tr*sr;
	k	=i+2;
	if(k < ESMp1){
	  tr	=2.*(mm+1-k)*Kr;
	  ti	=2.*(mm+1-k)*Ki;
	  Rr[k]	+=tr*sr-ti*si;
	  Ri[k]	+=ti*sr+tr*si;
	}
	k	=i-2;
	if(k > -ESMp1){
	  tr	=2.*(mm+1+k)*Kr;
	  ti	=2.*(mm+1+k)*Ki;
	  Rr[k]	+=tr*sr+ti*si;
	  Ri[k]	+=ti*sr-tr*si;
	}
      }
    }
    /*Making corrections with modified RHS;*/
    ur	=vi	+ES2Mp1;
    ui	=ur	+ES2Mp1;
    i	=-ESMp;
    for(k=-ESMp; k < ESMp1; k++){
      if(k != -mm && k != mm){
	i	=k-2;
	sr	=k ? (mm*mm-k*k)*K0 : ESaR0[0]*tD;
	si	=0.;
	if(k != 2 && k+ESMp > 1){
	  s=(mm-k)*(mm-k+2);
	  Cr=s*Kr;
	  Ci=s*Ki;
	  sr	+=Cr*EZvr[i]-Ci*vi[i];
	  si	+=Ci*EZvr[i]+Cr*vi[i];
	  s	=1./(sr*sr+si*si);
	  sr*= s;
	  si*=-s;
	  tr=Rr[k]+Cr*ur[i]-Ci*ui[i];
	  ti=Ri[k]+Ci*ur[i]+Cr*ui[i];
	  ur[k]	=-tr*sr+ti*si;
	  ui[k]	=-ti*sr-tr*si;
	}
	else{
	  sr=1./sr;
	  ur[k]	=-Rr[k]*sr;
	  ui[k]	=-Ri[k]*sr;
	}
	if(k+2){
	  s	=(mm+k)*(mm+k+2);
	  tr= s*Kr;
	  ti=-s*Ki;
	  EZvr[k]	=-tr*sr+ti*si;
	  vi[k]	=-ti*sr-tr*si;
	}
	else{
	  EZvr[k]	=0.;
	  vi[k]	=0.;
	}
      }
      else{
	ur[k]	=0.;
	ui[k]	=0.;
	EZvr[k]	=0.;
	vi[k]	=0.;
      }
    }
    k	=ESMp;
    EZvr[k]	=0.;
    vi[k]	=0.;
    k--;
    EZvr[k]	=0.;
    vi[k]	=0.;
    k	=ESMp1;
    while(k > -ESMp){
      k--;
      i	=k+2;
      ur[k]	+=EZvr[k]*ur[i]-vi[k]*ui[i];
      ui[k]	+=vi[k]*ur[i]+EZvr[k]*ui[i];
      if(k){
	pr[k]	+=ur[k];
	pi[k]	+=ui[k];
	dpr[k]	+=mm*ur[k];
	dpi[k]	+=mm*ui[k];
      }
      pr[k]	*=x*am;
      pi[k]	*=x*am;
      dpr[k]	*=am;
      dpi[k]	*=am;
    }
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
  }
  /*}*/

  /*Calculation of $\vaY${*/
  bKr[-2]	+=Kr;
  bKi[-2]	-=Ki;
  bKr[0]	+=K0;
  bKr[2]	+=Kr;
  bKi[2]	+=Ki;
  bMr[-2]	+=Ki;
  bMi[-2]	+=Kr;
  bMr[2]	+=Ki;
  bMi[2]	-=Kr;
  bNr[-2]	-=Kr;
  bNi[-2]	+=Ki;
  bNr[0]	+=K0;
  bNr[2]	-=Kr;
  bNi[2]	-=Ki;
  pr	=gpr	+ESMp;
  pi	=gpi	+ESMp;
  dpr	=dgpr	+ESMp;
  dpi	=dgpi	+ESMp;
  yr	=Yr	+ESMp;
  yi	=Yi	+ESMp;
  for(m=0; m < ESMp1; m++){/*All $Y^m$*/
    for(i=-ESMp; i < ESMp1; i++){
      yr[i]	=0.;
      yi[i]	=0.;
    }
    for(k=-ESMp; k < ESMp1; k++){
      sr	=x*dpr[k];
      si	=x*dpi[k];
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	yr[i]	+=sr*bKr[j]-si*bKi[j];
	yi[i]	+=si*bKr[j]+sr*bKi[j];
      }
      sr	= k*pi[k];
      si	=-k*pr[k];
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	yr[i]	+=sr*bMr[j]-si*bMi[j];
	yi[i]	+=si*bMr[j]+sr*bMi[j];
      }
    }
    yr[0]	=0.;
    yi[0]	=0.;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
  }
  /*}*/

  /*Averaged solution;*/
  ESdgY[0]	=-ESgm[0]*ESaF[0]*ESLc[0]/ESaR0[0];
  ESaY[0]	=K0*ESdgY[0]/ESaR0[0];
  ESgY02a[0]	=ESdgY[0];
  s		=x*x;
  i		=M0Mm;
  Yi[i]		=s*(ESaY[0]-0.125*s*ESjp[0]*ESVc2[0]);
  i++;
  Yi[i]		=EZcr2*x*x*ESdgY[0];
  return(0);
}

int Start2DgY(double x)
{
  int i,j,k,m,mm;
  int M,M1;
  double Bf[512];
  double x1r,x1i,tyr,by0;
  double *Jr,*Ji;
  double *bKr,*bKi;
  double *bMr,*bMi;
  double tD;
  double s,sr,si,tr,ti;
  double K0,Kr,Ki;
  double Cr,Ci;
  double *EZvr,*vi;
  double *pr,*pi,*dpr,*dpi;
  double *yr,*yi,*Rr,*Ri;
  double d2gY;
  double am;
  
  bKr	=Bf	+ESFp;
  k	=2*ES2Mp+1;
  bKi	=bKr	+k;
  bMr	=bKi	+k;
  bMi	=bMr	+k;
  Jr	=bMi	+k;
  Ji	=Jr	+ESMp1;
  EZvr	=Ji	+ES2Mp1;
  vi	=EZvr	+ES2Mp1;

  /*Getting geometry near the origin{*/
  s	=EZcr2*x;
  tyr	=EZcr2*ESsb1a[0];
  x1r	= rcT1a[k];
  x1i	=-rsT1a[k];

  tD	=4.*x1r*tyr;
  s	=1./tD;
  sr	=x1r*x1r;
  si	=x1i*x1i;
  tr	=tyr*tyr;
  K0	=2.*(tr+sr+si)*s;
  Kr	=(tr+si-sr)*s;
  Ki	=2.*x1r*x1i*s;
  /*}*/

  for(k=0; k < M0Mm; k++){
    gpr[k]	=0.;
    gpi[k]	=0.;
    dgpr[k]	=0.;
    dgpi[k]	=0.;
    Yr[k]	=0.;
    Yi[k]	=0.;
    dYr[k]	=0.;
    dYi[k]	=0.;
  }
  /*Main order solution{*/
  pr	=gpr+ESMp;
  pi	=gpi+ESMp;
  dpr	=dgpr+ESMp;
  dpi	=dgpi+ESMp;
  yr	=Yr+ESMp;
  yi	=Yi+ESMp;

  ESSetSplDCr(0.);
  splRDCr2(Jr,NULL,&d2gY,ESEqSolvInCr);

  Jr[0]	=-2.*K0*d2gY/(ESaR0[0]*ESLc[0]);
  Ji[0]	=0.;
  yr[0]	=EZcr2*d2gY;
  pr[0]	=EZcr2*d2gY;
  dpr[0]=d2gY;

  yr	+=ES2Mp1;
  yi	+=ES2Mp1;
  pr	+=ES2Mp1;
  pi	+=ES2Mp1;
  dpr	+=ES2Mp1;
  dpi	+=ES2Mp1;
  yr[1]	=1.;
  pr[1]	=yr[1];
  dpr[1]=yr[1];
  Jr[1]	=0.;
  Ji[1]	=0.;
  yr	+=ES2Mp1;
  yi	+=ES2Mp1;
  pr	+=ES2Mp1;
  pi	+=ES2Mp1;
  dpr	+=ES2Mp1;
  dpi	+=ES2Mp1;
  yr[2]	=1.;
  pr[2]	=yr[2];
  dpr[2]=2.*yr[2];
  s	=-8./(ESaR0[0]*tD);
  Jr[2]	= s*Kr;
  Ji[2]	=-s*Ki;
  for(m=3; m < ESMp1; m++){
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    i	=-m;
    EZvr[i]	=0.;
    vi[i]	=0.;
    for(k=i+2; k < m; k+=2,i+=2){
      if(k != 2){
	s=(m-k)*(m-k+2);
	Cr=s*Kr;
	Ci=s*Ki;
	s	=k ? (m*m-k*k)*K0 : ESaR0[0]*tD;
	sr=s+Cr*EZvr[i]-Ci*vi[i];
	si=Ci*EZvr[i]+Cr*vi[i];
	s	=1./(sr*sr+si*si);
	sr*= s;
	si*=-s;
	tr=Cr*yr[i]-Ci*yi[i];
	ti=Ci*yr[i]+Cr*yi[i];
	yr[k]	=-tr*sr+ti*si;
	yi[k]	=-ti*sr-tr*si;
      }
      else{
	s	=(m*m-k*k)*K0;
	sr=1./s;
	si=0.;
	yr[k]	=0.;
	yi[k]	=0.;
      }
      if(k+2){
	s	=(m+k)*(m+k+2);
	tr= s*Kr;
	ti=-s*Ki;
	EZvr[k]	=-tr*sr+ti*si;
	vi[k]	=-ti*sr-tr*si;
      }
      else{
	EZvr[k]	=0.;
	vi[k]	=0.;
      }
    }
    i	=m;
    yr[i]	=1.;
    yi[i]	=0.;
    j	=2-m;
    while(i > j){
      k	=i-2;
      yr[k]	+=EZvr[k]*yr[i]-vi[k]*yi[i];
      yi[k]	+=vi[k]*yr[i]+EZvr[k]*yi[i];
      i	-=2;
    }
    Jr[m]	=yr[0];
    Ji[m]	=yi[0];
    yr[0]	=0.;
    yi[0]	=0.;
    for(i=-m+2; i <= m; i++){
      pr[i]	=yr[i];
      pi[i]	=yi[i];
      dpr[i]	=m*yr[i];
      dpi[i]	=m*yi[i];
    }

  }  
  /*}*/

  /*Corrections to the EZmain order{*/
  pr	=gpr+ESMp;
  pi	=gpi+ESMp;
  dpr	=dgpr+ESMp;
  dpi	=dgpi+ESMp;
  yr	=Yr+ESMp;
  yi	=Yi+ESMp;
  am	=1.;
  for(m=0; m < ESMp1; m++){/*All fundamental vectors*/
    am	*=x;
    if(m == 1) am	=1.;
    for(k=-ESMp; k < ESMp1; k++){
      pr[k]	*=x*am;
      pi[k]	*=x*am;
      dpr[k]	*=am;
      dpi[k]	*=am;
    }
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
  }
  /*}*/
  /*Calculation of $\vaY${*/
  bKr[-2]	=Kr;
  bKi[-2]	=-Ki;
  bKr[0]	=K0;
  bKi[0]	=0.;
  bKr[2]	=Kr;
  bKi[2]	=Ki;
  bMr[-2]	=Ki;
  bMi[-2]	=Kr;
  bMr[2]	=Ki;
  bMi[2]	=-Kr;

  pr	=gpr	+ESMp;
  pi	=gpi	+ESMp;
  dpr	=dgpr	+ESMp;
  dpi	=dgpi	+ESMp;
  yr	=Yr	+ESMp;
  yi	=Yi	+ESMp;
  for(m=0; m < ESMp1; m++){/*All $Y^m$*/
    for(i=-ESMp; i < ESMp1; i++){
      yr[i]	=0.;
      yi[i]	=0.;
    }
    for(k=-ESMp; k < ESMp1; k++){
      sr	=x*dpr[k];
      si	=x*dpi[k];
      for(j=-2; j < 3; j+=2){
	i	=k+j;
	if(-ESMp1 < i  && i < ESMp1){
	  yr[i]	+=sr*bKr[j]-si*bKi[j];
	  yi[i]	+=si*bKr[j]+sr*bKi[j];
	}
      }
      sr	= k*pi[k];
      si	=-k*pr[k];
      for(j=-2; j < 3; j+=4){
	i	=k+j;
	if(-ESMp1 < i  && i < ESMp1){
	  yr[i]	+=sr*bMr[j]-si*bMi[j];
	  yi[i]	+=si*bMr[j]+sr*bMi[j];
	}
      }
    }
    yr[0]	=0.;
    yi[0]	=0.;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
  }
  /*}*/

  /*Averaged solution;*/

  ESdgY[0]	=d2gY;
  ESaY[0]	=K0*ESdgY[0]/ESaR0[0];
  ESgY02a[0]	=ESdgY[0];
  s		=x*x;
  i		=M0Mm;
  Yi[i]		=s*(ESaY[0]-0.125*s*ESjp[0]*ESVc2[0]);
  i++;
  Yi[i]		=EZcr2*x*x*ESdgY[0];
  return(0);
}

int Start2Dq(double x)
{
  int i,j,k,m;
  double s,sr,si,f,d,am,sj;
  double *pr,*pi,*dpr,*dpi,*xr,*xi,*yr,*yi,*lr,*li;
  double r1r,r1i,z1i;
  double x1r,x1i,X1r,X1i;
  double t1r,t1i,T1r,T1i;
  double tr,ti;

  double Er,Ei,K0,Br,Bi;
  double *ur,*ui,*EZvr,*vi;

  K0	=ESg22c[0];
  k	=2*ESNa1;
  Er	=ESg22c[k];
  Ei	=-ESg22s[k];
  sj	=2.*ESgm[0]*K0*ESaF[0]/ESaR0[0];

  /*Same as in Start2D{*/
  /*Getting metric tensor{*/
  ESSetSplA(x);
  splRA(wg12c,NULL,ESg12c,ESg12c2);
  wg12s[0]	=0.;
  splRA(wg22c,NULL,ESg22c,ESg22c2);
  wg22c[0]	*=x;
  wg22s[0]	=0.;
  k	=0;
  for(i=1; i < ES2Mp1; i++){
    k	+=ESNa1;
    splRA(wg12c+i,NULL,ESg12c+k,ESg12c2+k);
    wg12c[-i]	=wg12c[i];
    splRA(wg22c+i,NULL,ESg22c+k,ESg22c2+k);
    wg22c[i]	*=x;
    wg22c[-i]	=wg22c[i];
    splRA(wg12s-i,NULL,ESg12s+k,ESg12s2+k);
    wg12s[i]	=-wg12s[-i];
    splRA(wg22s-i,NULL,ESg22s+k,ESg22s2+k);
    wg22s[-i]	*=x;
    wg22s[i]	=-wg22s[-i];
  }
  z1i	=-EZcr2*ESsb1a[0];
  k	=ESNa1;
  r1r	= rcT1a[k];
  r1i	=-rsT1a[k];
  wuc[1]	=0.;
  wus[1]	=0.;
  wuc[-1]	=0.;
  wus[-1]	=0.;
  s		=EZcr2*x;
  wuc[0]	=s*rcT2a[0];
  wus[0]	=0.;
  for(i=2; i < ES2Mp1; i++){
    k	+=ESNa1;
    wuc[i]	= s*rcT2a[k];
    wus[i]	=-s*rsT2a[k];
    wuc[-i]	= wuc[i];
    wus[-i]	=-wus[i];
  }
  x1r	= r1r-z1i;
  x1i	= r1i;
  X1r	= r1r+z1i;
  X1i	=-r1i;
  t1r	= r1r+z1i;
  t1i	= r1i;
  T1r	= r1r-z1i;
  T1i	=-r1i;
  /*}*/  for(k=0; k < M0Mm; k++){
    gpr[k]	=0.;
    gpi[k]	=0.;
    dgpr[k]	=0.;
    dgpi[k]	=0.;
    Yr[k]	=0.;
    Yi[k]	=0.;
  }
  pr	=gpr	+ESMp;
  pi	=gpi	+ESMp;
  dpr	=dgpr	+ESMp;
  dpi	=dgpi	+ESMp;
  yr	=Yr	+ESMp;
  yi	=Yi	+ESMp;
  lr	=dYr	+ESMp;
  li	=dYi	+ESMp;

  /*m=0 polynomial solution{*/
  d	=ESaR0[0];
  f	=-EZcr4*x*sj;
  s	=f*(X1r*t1r-X1i*t1i+x1r*T1r-x1i*T1i);
  sr	=f*(x1r*t1r-x1i*t1i);
  si	=f*(x1r*t1i+x1i*t1r);
  
  pr[0]		=x*d*s;
  dpr[0]	=2.*d*s;
  f		=EZcr4*x;
  s		*=f;
  pr[1]		=x*s*r1r;
  dpr[1]	=3.*s*r1r;
  pi[1]		=x*s*r1i;
  dpi[1]	=3.*s*r1i;
  
  pr[2]		=x*d*sr;
  pi[2]		=x*d*si;
  dpr[2]	=2.*d*sr;
  dpi[2]	=2.*d*si;
  sr		*=f;
  si		*=f;
  if(3 < ESMp1){
    s	=sr*r1r-si*r1i;
    pr[3]	+=x*s;
    dpr[3]	+=3.*s;
    s	=si*r1r+sr*r1i;
    pi[3]	+=x*s;
    dpi[3]	+=3.*s;
  }
  s	=sr*r1r+si*r1i;
  pr[1]	+=x*s;
  dpr[1]+=3.*s;
  s	=si*r1r-sr*r1i;
  pi[1]	+=x*s;
  dpi[1]+=3.*s;

  f	=-0.5*x*x*d*sj;
  d	*=-1.5*x*sj;
  for(i=1; i < ESMp1; i++){
    j	=i-1;
    s	=r1r*wuc[j]-r1i*wus[j];
    pr[i]	+=f*s;
    dpr[i]	+=d*s;
    s	=sr*wus[j]+si*wuc[j];
    pi[i]	+=f*s;
    dpi[i]	+=d*s;
    j	=i+1;
    s	=r1r*wuc[j]+r1i*wus[j];
    pr[i]	+=f*s;
    dpr[i]	+=d*s;
    s	=sr*wus[j]-si*wuc[j];
    pi[i]	+=f*s;
    dpi[i]	+=d*s;
  }
  wus[0]=EZcr2*x*rsT2a[0];
  s	=wus[0]*z1i;
  pi[1]	+=f*s;
  dpi[1]+=d*s;
  for(i=1; i < ESMp1; i++){
    pr[-i]	= pr[i];
    pi[-i]	=-pi[i];
    dpr[-i]	= dpr[i];
    dpi[-i]	=-dpi[i];
  }
  /*}*/

  s	=EZcr4/ESaR0[0];
  t1r	*=s;
  t1i	*=s;
  T1r	*=s;
  T1i	*=s;
  yr[0]	= 1.;
  yi[0]	= 0.;
  am	=1.;
  for(m=1; m < ESMp1; m++){/*$m\neq1$ solutions*/
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    am	*=x;
    for(k=-m+1; k < m; k++){
      sr	=yr[k-ES2Mp1];
      si	=yi[k-ES2Mp1];
      /*Second correction O(a);*/
      f	=m*x;
      d	=m*(m+1);
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	s	=sr*wuc[j]-si*wus[j];
	pr[i]	+=f*s;
	dpr[i]	+=d*s;
	s	=sr*wus[j]+si*wuc[j];
	pi[i]	+=f*s;
	dpi[i]	+=d*s;
      }
      /* Main order contribution;*/
      i	=k+1;
      s		=x1r*sr-x1i*si;
      yr[i]	=x*s;
      pr[i]	+=yr[i];
      dpr[i]	+=m*s;
      s		=x1r*si+x1i*sr;
      yi[i]	=x*s;
      pi[i]	+=yi[i];
      dpi[i]	+=m*s;
      i	=k-1;
      s		=X1r*sr-X1i*si;
      yr[i]	=x*s;
      pr[i]	+=yr[i];
      dpr[i]	+=m*s;
      s		=X1r*si+X1i*sr;
      yi[i]	=x*s;
      pi[i]	+=yi[i];
      dpi[i]	+=m*s;
    }
    f	=1./m;
    d	=(m+1)*f;
    f	*=x;
    /*First correction $O(a)$*/
    for(k=-m; k <= m; k++){
      sr	=yr[k];
      si	=yi[k];
      i		=k+1;
      if(i < ESMp1){
	s	=t1r*sr-t1i*si;
	pr[i]	+=f*s;
	dpr[i]	+=d*s;
	s	=t1r*si+t1i*sr;
	pi[i]	+=f*s;
	dpi[i]	+=d*s;
      }
      i	=k-1;
      if(i > -ESMp1){
	s	=T1r*sr-T1i*si;
	pr[i]	+=f*s;
	dpr[i]	+=d*s;
	s	=T1r*si+T1i*sr;
	pi[i]	+=f*s;
	dpi[i]	+=d*s;
      }	
    }

    /*Normalization;$*/
    sr	=pr[m];
    si	=pi[m];
    s	=1./(sr*sr+si*si);
    sr	*=s;
    s	*=si;
    si	=pi[-m]*sr+pr[-m]*s;
    sr	=pr[-m]*sr-pi[-m]*s;
    for(i=-ESMp; i < ESMp1; i++){
      lr[i]	= pr[-i];
      li[i]	=-pi[-i];
    }
    for(i=-ESMp; i < ESMp1; i++){
      pr[i]	-=sr*lr[i]-si*li[i];
      pi[i]	-=sr*li[i]+si*lr[i];
    }
    for(i=-ESMp; i < ESMp1; i++){
      lr[i]	= dpr[-i];
      li[i]	=-dpi[-i];
    }
    for(i=-ESMp; i < ESMp1; i++){
      dpr[i]	-=sr*lr[i]-si*li[i];
      dpi[i]	-=sr*li[i]+si*lr[i];
    }
    if(m == 1){
      sr	=dpr[m];
      si	=dpi[m];
      s	=1./(sr*sr+si*si);
    }
    else{
      sr	=pr[m];
      si	=pi[m];
      s	=am/(sr*sr+si*si);
    }
    sr	*=s;
    si	*=-s;
    for(i=-ESMp; i < ESMp1; i++){
      s	=pr[i];
      pr[i]	=sr*s-si*pi[i];
      pi[i]	=si*s+sr*pi[i];
      s	=dpr[i];
      dpr[i]	=sr*s-si*dpi[i];
      dpi[i]	=si*s+sr*dpi[i];
    }
  }

  /*Normalization of $m=0$;*/
  pr	=gpr	+ESMp;
  pi	=gpi	+ESMp;
  dpr	=dgpr	+ESMp;
  dpi	=dgpi	+ESMp;
  i	=ESMp+2*ES2Mp1;
  yr	=gpr	+i;
  yi	=gpi	+i;
  lr	=dgpr	+i;
  li	=dgpi	+i;
  sr	=yr[2]; 
  si	=yi[2];
  s	=1./(sr*sr+si*si);
  sr	*=s;	/*$\R{1}{\gy^{*2}_2}$*/
  s	*=si;
  tr	=pr[-2]*sr-pi[-2]*s;	/*$\R{\gy^0_{-2}}{\gy^{*2}_2}$*/
  ti	=pi[-2]*sr+pr[-2]*s;
  si	=pi[2]*sr-pr[2]*s;
  sr	=pr[2]*sr+pi[2]*s;	/*$\R{\gy^0_2}{\gy^2_2}$*/
  for(i=-ESMp; i < ESMp1; i++){
    pr[i]	-=sr*yr[i]-si*yi[i]+tr*yr[-i]+ti*yi[-i];
    pi[i]	-=si*yr[i]+sr*yi[i]+ti*yr[-i]-tr*yi[-i];
    dpr[i]	-=sr*lr[i]-si*li[i]+tr*lr[-i]+ti*li[-i];
    dpi[i]	-=si*lr[i]+sr*li[i]+ti*lr[-i]-tr*li[-i];
  }

  /*Elimination of $\gy^2_0$;*/
  s	=x*x;
  s	*=0.125*s*ESjp[0]*ESVc2[0]/wg22c[0];
  dpr[0]-=s;
  pr[0]	-=EZcr4*s;
  /*}*/

  pr	=gpr	+ESMp;
  pi	=gpi	+ESMp;
  dpr	=dgpr	+ESMp;
  dpi	=dgpi	+ESMp;
  yr	=Yr	+ESMp;
  yi	=Yi	+ESMp;
  lr	=dYr	+ESMp;
  li	=dYi	+ESMp;
  for(m=1; m < ESMp1; m++){/*$m\neq1$ solutions elimination of $gy^m_0$*/
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    am	*=x;
    for(k=-m+1; k < m; k++){
      sr	=yr[k-ES2Mp1];
      si	=yi[k-ES2Mp1];
      /* Main order contribution;*/
      i	=k+1;
      s		=x1r*sr-x1i*si;
      yr[i]	=x*s;
      pr[i]	+=yr[i];
      dpr[i]	+=m*s;
      s		=x1r*si+x1i*sr;
      yi[i]	=x*s;
      pi[i]	+=yi[i];
      dpi[i]	+=m*s;
      i	=k-1;
      s		=X1r*sr-X1i*si;
      yr[i]	=x*s;
      pr[i]	+=yr[i];
      dpr[i]	+=m*s;
      s		=X1r*si+X1i*sr;
      yi[i]	=x*s;
      pi[i]	+=yi[i];
      dpi[i]	+=m*s;
    }
    if(0 && m%2 == 0){
      k	=m;
      EZvr[k]	=0.;
      vi[k]	=0.;
      ur[k]	=0.;
      ui[k]	=0.;
      k	-=2;
      while(k > -m){
	s	=(double)(m+k)*(m+k+2.);
	Bi	=-s*Ei;
	Br	=s*Er;
	sr	=Br*EZvr[k+2]-Bi*vi[k+2]+(m*m-k*k)*K0;
	si	=Br*vi[k+2]+Bi*EZvr[k+2];
	s	=1./(sr*sr+si*si);
	sr	*= s;
	si	*=-s;
	s	=(double)(m-k)*(m-k+2.);
	EZvr[k]	=-(sr*Er-si*Ei)*s;
	vi[k]	=-(sr*Ei+si*Er)*s;
	s	=-sr*Br+si*Bi;
	ui[k]	=-sr*Bi-si*Br;
	if(k == 0){
	  s	-=pr[0];
	  ui[k]	-=pi[0];
	}
	ur[k]	=sr*s+si*ui[k];
	ui[k]	=si*s+sr*ui[k];
	k	-=2;
      }
      ur[k]	=0.;
      ui[k]	=0.;
      s	=m/x;
      while(k < m){
	ur[k]	+=EZvr[k]*ur[k-2]-vi[k]*ui[k-2];
	ui[k]	+=EZvr[k]*ui[k-2]+vi[k]*ur[k-2];
	pr[k]	+=ur[k];
	pi[k]	+=ui[k];
	dpr[k]	+=s*ur[k];
	dpi[k]	+=s*ui[k];
	k	+=2;
      }
    }
    for(k=-m+1; k < m; k++){
      sr	=yr[k-ES2Mp1];
      si	=yi[k-ES2Mp1];
      /*Second correction O(a);*/
      f	=m*x;
      d	=m*(m+1);
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	s	=sr*wuc[j]-si*wus[j];
	pr[i]	+=f*s;
	dpr[i]	+=d*s;
	s	=sr*wus[j]+si*wuc[j];
	pi[i]	+=f*s;
	dpi[i]	+=d*s;
      }
    }
    f	=1./m;
    d	=(m+1)*f;
    f	*=x;
    /*First correction $O(a)$*/
    for(k=-m; k <= m; k++){
      sr	=yr[k];
      si	=yi[k];
      i		=k+1;
      if(i < ESMp1){
	s	=t1r*sr-t1i*si;
	pr[i]	+=f*s;
	dpr[i]	+=d*s;
	s	=t1r*si+t1i*sr;
	pi[i]	+=f*s;
	dpi[i]	+=d*s;
      }
      i	=k-1;
      if(i > -ESMp1){
	s	=T1r*sr-T1i*si;
	pr[i]	+=f*s;
	dpr[i]	+=d*s;
	s	=T1r*si+T1i*sr;
	pi[i]	+=f*s;
	dpi[i]	+=d*s;
      }	
    }
    if(0 && m%2){
      k	=m+1;
      EZvr[k]	=0.;
      vi[k]	=0.;
      ur[k]	=0.;
      ui[k]	=0.;
      k	-=2;
      while(k > -m){
	s	=(double)(m+1+k)*(m+1+k+2.);
	Bi	=-s*Ei;
	Br	=s*Er;
	sr	=Br*EZvr[k+2]-Bi*vi[k+2]+((m+1)*(m+1)-k*k)*K0;
	si	=Br*vi[k+2]+Bi*EZvr[k+2];
	s	=1./(sr*sr+si*si);
	sr	*= s;
	si	*=-s;
	s	=(double)(m+1-k)*(m+1-k+2.);
	EZvr[k]	=-(sr*Er-si*Ei)*s;
	vi[k]	=-(sr*Ei+si*Er)*s;
	s	=-sr*Br+si*Bi;
	ui[k]	=-sr*Bi-si*Br;
	if(k == 0){
	  s	-=pr[0];
	  ui[k]	-=pi[0];
	}
	ur[k]	=sr*s+si*ui[k];
	ui[k]	=si*s+sr*ui[k];
	k	-=2;
      }
      ur[k]	=0.;
      ui[k]	=0.;
      s	=(m+1)/x;
      while(k < m){
	ur[k]	+=EZvr[k]*ur[k-2]-vi[k]*ui[k-2];
	ui[k]	+=EZvr[k]*ui[k-2]+vi[k]*ur[k-2];
	pr[k]	+=ur[k];
	pi[k]	+=ui[k];
	dpr[k]	+=s*ur[k];
	dpi[k]	+=s*ui[k];
	k	+=2;
      }
    }
    else{
      pr[0]	=0.;
      pi[0]	=0.;
      dpr[0]	=0.;
      dpi[0]	=0.;
    }
    /*Normalization;$*/
    sr	=pr[m];
    si	=pi[m];
    s	=1./(sr*sr+si*si);
    sr	*=s;
    s	*=si;
    si	=pi[-m]*sr+pr[-m]*s;
    sr	=pr[-m]*sr-pi[-m]*s;
    for(i=-ESMp; i < ESMp1; i++){
      lr[i]	= pr[-i];
      li[i]	=-pi[-i];
    }
    for(i=-ESMp; i < ESMp1; i++){
      pr[i]	-=sr*lr[i]-si*li[i];
      pi[i]	-=sr*li[i]+si*lr[i];
    }
    for(i=-ESMp; i < ESMp1; i++){
      lr[i]	= dpr[-i];
      li[i]	=-dpi[-i];
    }
    for(i=-ESMp; i < ESMp1; i++){
      dpr[i]	-=sr*lr[i]-si*li[i];
      dpi[i]	-=sr*li[i]+si*lr[i];
    }
    if(m == 1){
      sr	=dpr[m];
      si	=dpi[m];
      s	=1./(sr*sr+si*si);
    }
    else{
      sr	=pr[m];
      si	=pi[m];
      s	=am/(sr*sr+si*si);
    }
    sr	*=s;
    si	*=-s;
    for(i=-ESMp; i < ESMp1; i++){
      s	=pr[i];
      pr[i]	=sr*s-si*pi[i];
      pi[i]	=si*s+sr*pi[i];
      s	=dpr[i];
      dpr[i]	=sr*s-si*dpi[i];
      dpi[i]	=si*s+sr*dpi[i];
    }
  }

  /*Calculation of $\vaY${*/
  pr	=gpr	+ESMp;
  pi	=gpi	+ESMp;
  dpr	=dgpr	+ESMp;
  dpi	=dgpi	+ESMp;
  yr	=Yr	+ESMp;
  yi	=Yi	+ESMp;
  for(m=0; m < ESMp1; m++){/*All $Y^m$*/
    for(i=-ESMp; i < ESMp1; i++){
      yr[i]	=0.;
      yi[i]	=0.;
    }
    for(k=-ESMp; k < ESMp1; k++){
      sr	=dpr[k];
      si	=dpi[k];
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	yr[i]	+=sr*wg22c[j]-si*wg22s[j];
	yi[i]	+=si*wg22c[j]+sr*wg22s[j];
      }
      sr	= k*pi[k];
      si	=-k*pr[k];
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	yr[i]	+=sr*wg12c[j]-si*wg12s[j];
	yi[i]	+=si*wg12c[j]+sr*wg12s[j];
      }
    }
    yr[0]	=0.;
    yi[0]	=0.;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
  }
  /*}*/

  /*Averaged solution;*/
  ESdgY[0]	=-ESgm[0]*ESaF[0]*ESLc[0]/ESaR0[0];
  ESaY[0]	=K0*ESdgY[0];
  ESgY02a[0]	=ESdgY[0];
  s		=x*x;
  i		=M0Mm;
  Yi[i]		=s*(ESaY[0]-0.125*s*ESjp[0]*ESVc2[0]);
  i++;
  Yi[i]		=EZcr2*x*x*ESdgY[0];
  return(0);
}

int ECStartEq2DsjHole()
{
  int m,k;
  double *pr,*pi;
  double *yr,*yi;

  for(k=0; k < ES2Mp1; k++){
    gpr[k]	=0.;
    gpi[k]	=0.;
    Yr[k]	=0.;
    Yi[k]	=0.;
  }    

  pr	=gpr+ESMp;
  pi	=gpi+ESMp;
  yr	=Yr+ESMp;
  yi	=Yi+ESMp;

  for(m=1; m < ESMp1; m++){
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    for(k=-ESMp; k < ESMp1; k++){
      pr[k]	=0.;
      pi[k]	=0.;
      yr[k]	=0.;
      yi[k]	=0.;
    }    
    yr[m]	=1.;
  }

  k		=M0Mm;
  Yi[k]		=0.;
  ESaY[0]	=0.;
  ESdgY[0]	=0.;
  ESgY02a[0]	=0.;
  k++;
  Yi[k]		=0.;
  return(0);
}

int ECStartEq2Dsq(double x)
{
  int m,mm,mp,k,kk,km,kp,K;
  double am,bK0,Er,Ei,A,tr,ti,sr,si,Tr,Ti;
  double *pr,*pi,*dpr,*dpi,*ur,*ui,*EZvr,*vi;
  double *xr,*xi,*yr,*yi,*dyr,*dyi;

  gpi	=gpr+M0Mm;
  Yr	=gpi+M0Mm;
  Yi	=Yr+M0Mm;
  dgpi	=dgpr+M0Mm;
  dYr	=dgpi+M0Mm;
  dYi	=dYr+M0Mm;

  K	=0;
  bK0	=ESg22c[0];
  k	=2*ESNa1;
  Er	=ESg22c[k];
  Ei	=-ESg22s[k];
  
  A	=1./ESsa[1];
  tr	=EZcr6*ESsa[1];
  kk	=0;
  mp	=1;
  sr	=1./ESLc[0];
  for(k=0; k < ES2Mp1; k++){
    if(k%2){
      wduc[k]	=x*ESLc1[kk]*sr;
      wdus[k]	=x*ESLs1[kk]*sr;
      wg11c[k]	=x*((ESg11c[mp]-ESg11c[kk])*A-(ESg11c2[mp]+2.*ESg11c2[kk])*tr);
      wg11s[k]	=x*((ESg11s[mp]-ESg11s[kk])*A-(ESg11s2[mp]+2.*ESg11s2[kk])*tr);
      wg12c[k]	=x*((ESg12c[mp]-ESg12c[kk])*A-(ESg12c2[mp]+2.*ESg12c2[kk])*tr);
      wg12s[k]	=x*((ESg12s[mp]-ESg12s[kk])*A-(ESg12s2[mp]+2.*ESg12s2[kk])*tr);
      wg22c[k]	=x*ESg22c1[kk];
      wg22s[k]	=x*ESg22s1[kk];
    }
    else{
      wduc[k]	=0.;
      wdus[k]	=0.;
      wg11c[k]	=0.;
      wg11s[k]	=0.;
      wg12c[k]	=0.;
      wg12s[k]	=0.;
      wg22c[k]	=0.;
      wg22s[k]	=0.;
    }
    kk	+=ESNa1;
    mp	+=ESNa1;
  }

  k	=2*M0Mm+ESMp;
  ur	=EZgper+k;
  ui	=EZgpei+k;
  EZvr	=ur+M0Mm;
  vi	=ui+M0Mm;
  xr	=EZvr+M0Mm;
  xi	=vi+M0Mm;

  for(m=1; m < ESMp1; m++){
    xr	+=ES2Mp1;
    xi	+=ES2Mp1;
    EZvr	+=ES2Mp1;
    vi	+=ES2Mp1;
    ur	+=ES2Mp1;
    ui	+=ES2Mp1;

    /* Decomposition coefficients{*/
    mm 	=-m;
    k	=ESMp;
    kk	=ESMp1;
    ur[k]	=0.;
    ui[k]	=0.;
    xi[k]	=0.;
    if(k == m){
      EZvr[k]	=0.;
      vi[k]	=0.;
      xr[k]	=0.;
    }
    else{
      A		=1./((m-k)*(m+k)*bK0);
      si	=(m-k)*(m-k+2);
      xr[k]	=A;
      A		*=(double)((m-k)*(m-k+2));
      EZvr[k]	=-Er*A;
      vi[k]	=-Ei*A;
    }
    k--;
    ur[k]	=0.;
    ui[k]	=0.;
    xi[k]	=0.;
    if(k == m){
      EZvr[k]	=0.;
      vi[k]	=0.;
      xr[k]	=0.;
    }
    else{
      A		=1./((m-k)*(m+k)*bK0);
      si	=(m-k)*(m-k+2);
      xr[k]	=A;
      A		*=(double)((m-k)*(m-k+2));
      EZvr[k]	=-Er*A;
      vi[k]	=-Ei*A;
    }
    while(k > -ESMp){
      kk--;
      k--;
      if(k == m || k == 0 || k == mm){
	EZvr[k]	=0.;
	vi[k]	=0.;
	ur[k]	=0.;
	ui[k]	=0.;
	xr[k]	=0.;
	xi[k]	=0.;
      }
      else{
	A	=(m-k)*(m+k)*bK0;
	ti	=(m+k)*(m+kk);
	tr	=ti*Er;
	ti	*=-Ei;
	si	=(m-k)*(m-k+2);
	sr	=si*Er;
	si	*=Ei;
	Tr	=A+tr*EZvr[kk]-ti*vi[kk];
	Ti	=tr*vi[kk]+ti*EZvr[kk];
	A	=1./(Tr*Tr+Ti*Ti);
	Tr	*=-A;
	Ti	*=A;
	xr[k]	=-Tr;
	xi[k]	=-Ti;
	EZvr[k]	=sr*Tr-si*Ti;
	vi[k]	=sr*Ti+si*Tr;
	ur[k]	=tr*Tr-ti*Ti;
	ui[k]	=tr*Ti+ti*Tr;
      }
    }    
    /*}*/
  }

  k	=2*M0Mm+ESMp+ES2Mp1;
  ur	=EZgper+k;
  ui	=EZgpei+k;
  EZvr	=ur+M0Mm;
  vi	=ui+M0Mm;
  xr	=EZvr+M0Mm;
  xi	=vi+M0Mm;

  pr	=gpr+ESMp;
  pi	=gpi+ESMp;
  dpr	=dgpr+ESMp;
  dpi	=dgpi+ESMp;
  yr	=Yr+ESMp;
  yi	=Yi+ESMp;
  dyr	=dYr+ESMp;
  dyi	=dYi+ESMp;
  am	=x;
  for(m=1; m < ESMp1; m++){
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    mm 	=-m;
    /* Homogeneous solution;*/
    for(k=-ESMp; k < ESMp1; k++){
      pr[k]	=0.;
      pi[k]	=0.;
      dpr[k]	=0.;
      dpi[k]	=0.;
    }    
    k		=m;
    pr[k]	=am;
    pi[k]	=0.;
    dpr[k]	=m*am;
    dpi[k]	=0.;
    kk		=k;
    k		-=2;
    if(m%2 == 0)
      mm	=0;
    while(k > mm){
      pr[k]	=ur[k]*pr[kk]-ui[k]*pi[kk];
      pi[k]	=ur[k]*pi[kk]+ui[k]*pr[kk];
      kk	-=2;
      k		-=2;
    }    
    A		=m;
    while(kk < m){
      pr[kk]	+=EZvr[kk]*pr[k]-vi[kk]*pi[k];
      pi[kk]	+=vi[kk]*pr[k]+EZvr[kk]*pi[k];
      dpr[kk]	=A*pr[kk];
      dpi[kk]	=A*pi[kk];
      k		+=2;
      kk	+=2;
    }
    yr[0]	=bK0*dpr[0];
    yi[0]	=bK0*dpi[0];
    mp		=0;
    mm		=0;
    while(mp < m){
      mp++;
      mm--;
      yr[mp]	=bK0*dpr[mp];
      yi[mp]	=bK0*dpi[mp];
      yr[mm]	=bK0*dpr[mm];
      yi[mm]	=bK0*dpi[mm];
    }
    while(mp < ESMp){
      mp++;
      mm--;
      yr[mp]	=0.;
      yi[mp]	=0.;
      yr[mm]	=0.;
      yi[mm]	=0.;
    }
    tr		=dpr[0];
    ti		=dpi[0];
    yr[2]	+=Er*tr-Ei*ti;
    yi[2]	+=Er*ti+Ei*tr;
    yr[-2]	+=Er*tr+Ei*ti;
    yi[-2]	+=Er*ti-Ei*tr;
    k		=0;
    mp		=2;
    while(k < m){
      k++;
      mp++;
      tr	=dpr[k];
      ti	=dpi[k];
      Tr	=k*pr[k];
      Ti	=k*pi[k];
      if(mp < ESMp1){
	sr	=tr-Tr;
	si	=ti-Ti;
	yr[mp]	+=Er*sr-Ei*si;
	yi[mp]	+=Er*si+Ei*sr;
      }
      mm	=k-2;
      sr	=tr+Tr;
      si	=ti+Ti;
      yr[mm]	+=Er*sr+Ei*si;
      yi[mm]	+=Er*si-Ei*sr;
      mm	=-k;
      tr	=dpr[mm];
      ti	=dpi[mm];
      Tr	=mm*pr[mm];
      Ti	=mm*pi[mm];
      sr	=tr-Tr;
      si	=ti-Ti;
      mm	+=2;
      yr[mm]	+=Er*sr-Ei*si;
      yi[mm]	+=Er*si+Ei*sr;
      if(mp < ESMp1){
	mm	-=4;
	sr	=tr+Tr;
	si	=ti+Ti;
	yr[mm]	+=Er*sr+Ei*si;
	yi[mm]	+=Er*si-Ei*sr;
      }
    }
    /*Right Hand Side{*/
    /* Equation for $Y$;*/
    if(m < ESMp){
      sr	=-wg22c[0];
      dyr[0]	=0.;
      dyi[0]	=0.;
      mm	=0;
      for(mp=1; mp < ESMp1; mp++){
	mm--;
	dyr[mp]	=dpr[mp]*sr;
	dyi[mp]	=dpi[mp]*sr;
	dyr[mm]	=dpr[mm]*sr;
	dyi[mm]	=dpi[mm]*sr;
      }
      mm	=0;
      for(mp=1; mp < ESMp1; mp++){
	mm--;
	tr	=dpr[mp];
	ti	=dpi[mp];
	Tr	=dpr[mm];
	Ti	=dpi[mm];
	k	=1;
	kp	=mp+k;
	km	=mm-k;
	while(kp < ESMp1){
	  sr	=-wg22c[k];
	  si	=wg22s[k];
	  k++;
	  dyr[kp]	+=sr*tr-si*ti;
	  dyi[kp]	+=sr*ti+si*tr;
	  kp++;
	  dyr[km]	+=sr*Tr+si*Ti;
	  dyi[km]	+=sr*Ti-si*Tr;
	  km--;
	}
	k	=1;
	kp	=k+mm;
	km	=mp-k;
	while(kp < ESMp1){
	  sr	=-wg22c[k];
	  si	=wg22s[k];
	  k++;
	  dyr[kp]	+=sr*Tr-si*Ti;
	  dyi[kp]	+=sr*Ti+si*Tr;
	  kp++;
	  dyr[km]	+=sr*tr+si*ti;
	  dyi[km]	+=sr*ti-si*tr;
	  km--;
	}
      }

      mm	=0;
      sr	=wg12c[0];
      for(mp=1; mp < ESMp1; mp++){
	mm--;
	si	=sr*mp;
	dyr[mp]	-=si*pi[mp];
	dyi[mp]	+=si*pr[mp];
	dyr[mm]	+=si*pi[mm];
	dyi[mm]	-=si*pr[mm];
      }
      mm	=0;
      for(mp=1; mp < ESMp1; mp++){
	mm--;
	tr	=mp*pr[mp];
	ti	=mp*pi[mp];
	Tr	=mp*pr[mm];
	Ti	=mp*pi[mm];
	k	=1;
	kp	=mp+k;
	km	=mm-k;
	while(kp < ESMp1){
	  sr	=wg12s[k];
	  si	=wg12c[k];
	  k++;
	  dyr[kp]	+=sr*tr-si*ti;
	  dyi[kp]	+=sr*ti+si*tr;
	  kp++;
	  dyr[km]	+=sr*Tr+si*Ti;
	  dyi[km]	+=sr*Ti-si*Tr;
	  km--;
	}
	k	=1;
	kp	=k+mm;
	km	=mp-k;
	while(kp < ESMp1){
	  sr	=wg12s[k];
	  si	=wg12c[k];
	  k++;
	  dyr[kp]	-=sr*Tr-si*Ti;
	  dyi[kp]	-=sr*Ti+si*Tr;
	  kp++;
	  dyr[km]	-=sr*tr+si*ti;
	  dyi[km]	-=sr*ti-si*tr;
	  km--;
	}
      }
      /* Equation for $Y'$;*/
      A		=m+1;
      tr	=m*yr[0];
      ti	=m*yi[0];
      sr	=wg11c[0];
      yr[0]	-=dyr[0];
      yi[0]	-=dyi[0];
      dyr[0]	*=A;
      dyi[0]	*=A;
      mm	=0;
      for(mp=1; mp < ESMp1; mp++){
	mm--;
	yr[mp]	-=dyr[mp];
	yi[mp]	-=dyi[mp];
	yr[mm]	-=dyr[mm];
	yi[mm]	-=dyi[mm];

	k	=mp*mp;
	dyr[mp]	=A*dyr[mp]+wduc[mp]*tr+wdus[mp]*ti+k*pr[mp]*sr;
	dyi[mp]	=A*dyi[mp]+wduc[mp]*ti-wdus[mp]*tr+k*pi[mp]*sr;
	dyr[mm]	=A*dyr[mm]+wduc[mp]*tr-wdus[mp]*ti+k*pr[mm]*sr;
	dyi[mm]	=A*dyi[mm]+wduc[mp]*ti+wdus[mp]*tr+k*pi[mm]*sr;
      }
      mm	=0;
      for(mp=1; mp < ESMp1; mp++){
	mm--;
	tr	=mp*pr[mp];
	ti	=mp*pi[mp];
	Tr	=mm*pr[mm];
	Ti	=mm*pi[mm];
	k	=1;
	kp	=mp+k;
	km	=mm-k;
	while(kp < ESMp1){
	  sr	=wg11c[k];
	  si	=-wg11s[k];
	  k++;
	  dyr[kp]	+=kp*(sr*tr-si*ti);
	  dyi[kp]	+=kp*(sr*ti+si*tr);
	  kp++;
	  dyr[km]	+=km*(sr*Tr+si*Ti);
	  dyi[km]	+=km*(sr*Ti-si*Tr);
	  km--;
	}
	k	=1;
	kp	=k+mm;
	km	=mp-k;
	while(kp < ESMp1){
	  sr	=wg11c[k];
	  si	=-wg11s[k];
	  k++;
	  dyr[kp]	+=kp*(sr*Tr-si*Ti);
	  dyi[kp]	+=kp*(sr*Ti+si*Tr);
	  kp++;
	  dyr[km]	+=km*(sr*tr+si*ti);
	  dyi[km]	+=km*(sr*ti-si*tr);
	  km--;
	}
      }
      si	=wg12c[0];
      mm	=0;
      for(mp=1; mp < ESMp1; mp++){
	mm--;
	dyr[mp]	-=mp*dpi[mp]*si;
	dyi[mp]	+=mp*dpr[mp]*si;
	dyr[mm]	+=mp*dpi[mm]*si;
	dyi[mm]	-=mp*dpr[mm]*si;
      }
      mm	=0;
      for(mp=1; mp < ESMp1; mp++){
	mm--;
	tr	=dpr[mp];
	ti	=dpi[mp];
	Tr	=dpr[mm];
	Ti	=dpi[mm];
	k	=1;
	kp	=mp+k;
	km	=-kp;
	while(kp < ESMp1){
	  sr	=kp*wg12s[k];
	  si	=kp*wg12c[k];
	  k++;
	  dyr[kp]	+=sr*tr-si*ti;
	  dyi[kp]	+=sr*ti+si*tr;
	  kp++;
	  dyr[km]	+=sr*Tr+si*Ti;
	  dyi[km]	+=sr*Ti-si*Tr;
	  km--;
	}
	k	=1;
	kp	=k+mm;
	km	=-kp;
	while(kp < ESMp1){
	  sr	=kp*wg12s[k];
	  si	=kp*wg12c[k];
	  k++;
	  dyr[kp]	+=sr*Tr-si*Ti;
	  dyi[kp]	+=sr*Ti+si*Tr;
	  kp++;
	  dyr[km]	+=sr*tr+si*ti;
	  dyi[km]	+=sr*ti-si*tr;
	  km--;
	}
      }
    }
    /*}*/
    xr	+=ES2Mp1;
    xi	+=ES2Mp1;
    EZvr	+=ES2Mp1;
    vi	+=ES2Mp1;
    ur	+=ES2Mp1;
    ui	+=ES2Mp1;

    if(m < ESMp){
      /* Correction to the RHS because of frozenness; */
      Tr		=0.; 
      Ti		=0.;
      mm		=0;
      for(mp=1; mp < ESMp1; mp++){
	mm--;
	Tr	+=wduc[mp]*(pr[mm]+pr[mp])+wdus[mp]*(pi[mm]-pi[mp]);
	Ti	+=wduc[mp]*(pi[mm]+pi[mp])+wdus[mp]*(pr[mp]-pr[mm]);
      }
      dyr[0]	=-Tr;
      dyi[0]	=-Ti;
      xr[0]	=1.;

      /* Correction to the solution; */
      kk		=ESMp1;
      k		=ESMp;
      sr		=dyr[k];
      si		=dyi[k];
      dyr[k]	=xr[k]*sr-xi[k]*si;
      dyi[k]	=xi[k]*sr+xr[k]*si;
      k--;
      sr		=dyr[k];
      si		=dyi[k];
      dyr[k]	=xr[k]*sr-xi[k]*si;
      dyi[k]	=xi[k]*sr+xr[k]*si;
      while(k > -ESMp){
	kk--;
	k--;
	tr	=dyr[kk];
	ti	=dyi[kk];
	sr	=dyr[k];
	si	=dyi[k];
	dyr[k]	=xr[k]*sr-xi[k]*si+ur[k]*tr-ui[k]*ti;
	dyi[k]	=xi[k]*sr+xr[k]*si+ui[k]*tr+ur[k]*ti;
      }
      A		=m+1;
      pr[k]	+=dyr[k];
      pi[k]	+=dyi[k];
      dpr[k]	+=A*dyr[k];
      dpi[k]	+=A*dyi[k];
      kk		=k+1;
      pr[kk]	+=dyr[kk];
      pi[kk]	+=dyi[kk];
      dpr[kk]	+=A*dyr[kk];
      dpi[kk]	+=A*dyi[kk];
      kk++;
      while(kk < ESMp1){
	dyr[kk]	+=EZvr[kk]*dyr[k]-vi[kk]*dyi[k];
	dyi[kk]	+=vi[kk]*dyr[k]+EZvr[kk]*dyi[k];
	pr[kk]	+=dyr[kk];
	pi[kk]	+=dyi[kk];
	dpr[kk]	+=A*dyr[kk];
	dpi[kk]	+=A*dyi[kk];
	k++;
	kk++;
      }

      yr[0]	+=A*bK0*dyr[0];
      yi[0]	+=A*bK0*dyi[0];
      mp		=1;
      mm		=0;
      while(mp < ESMp1){
	mm--;
	yr[mp]	+=A*bK0*dyr[mp];
	yi[mp]	+=A*bK0*dyi[mp];
	yr[mm]	+=A*bK0*dyr[mm];
	yi[mm]	+=A*bK0*dyi[mm];
	mp++;
      }
      tr		=A*dyr[0];
      ti		=A*dyi[0];
      yr[2]	+=Er*tr-Ei*ti;
      yi[2]	+=Er*ti+Ei*tr;
      yr[-2]	+=Er*tr+Ei*ti;
      yi[-2]	+=Er*ti-Ei*tr;
      kp		=3;
      km		=-2;
      mp		=1;
      mm		=0;
      while(kp < ESMp1){
	mm--;
	km--;
	tr	=(A-mp)*dyr[mp];
	ti	=(A-mp)*dyi[mp];
	Tr	=(A-mp)*dyr[mm];
	Ti	=(A-mp)*dyi[mm];
	yr[kp]	+=Er*tr-Ei*ti;
	yi[kp]	+=Er*ti+Ei*tr;
	yr[km]	+=Er*Tr+Ei*Ti;
	yi[km]	+=Er*Ti-Ei*Tr;
	kp++;
	mp++;
      }
      kp		=-1;
      km		=2;
      mp		=1;
      mm		=0;
      while(mp < ESMp1){
	mm--;
	km--;
	tr	=(A+mp)*dyr[mp];
	ti	=(A+mp)*dyi[mp];
	Tr	=(A+mp)*dyr[mm];
	Ti	=(A+mp)*dyi[mm];
	yr[kp]	+=Er*tr+Ei*ti;
	yi[kp]	+=Er*ti-Ei*tr;
	yr[km]	+=Er*Tr-Ei*Ti;
	yi[km]	+=Er*Ti+Ei*Tr;
	kp++;
	mp++;
      }
    }
    yr[0]	=0.;
    yi[0]	=0.;
    am		*=x;
  }
  for(k=0; k < ES2Mp1; k++){
    gpr[k]	=0.;
    gpi[k]	=0.;
    Yr[k]	=0.;
    Yi[k]	=0.;
  }    
  k		=M0Mm;
  Yi[k]		=-bK0*ESLc[K]*ESaF[0]*ESgm[0]/R0;
  ESaY[0]	=Yi[k];
  ESdgY[0]	=Yi[k]/bK0;
  ESgY02a[0]	=Yi[k]/bK0;
  Yi[k]		*=x*x;
  dgpr[ESMp]	=x*ESgY02a[0];
  k++;
  Yi[k]		=EZcr2*x*dgpr[ESMp];
  gpr[ESMp]	=Yi[k];
  Yr[ESMp]	=x*bK0*dgpr[ESMp];
  tr		=x*Er*dgpr[ESMp];
  ti		=x*Ei*dgpr[ESMp];
  Yr[ESMp+2]	=tr;
  Yi[ESMp+2]	=ti;
  Yr[ESMp-2]	=tr;
  Yi[ESMp-2]	=-ti;
  Yr[ESMp]	=0.;
  return(0);
}

int Start2Daa(double x)
{
  int i,j,k,m;
  double s,sr,si,f,d,am,sj;
  double *pr,*pi,*dpr,*dpi,*xr,*xi,*yr,*yi,*lr,*li;
  double r1r,r1i,tyr;
  double x1r,x1i,X1r,X1i;
  double t1r,t1i,T1r,T1i;
  double tr,ti;

  /*Getting metric tensor{*/
  ESSetSplA(x);
  splRA(wg12c,NULL,ESg12c,ESg12c2);
  wg12s[0]	=0.;
  splRA(wg22c,NULL,ESg22c,ESg22c2);
  wg22c[0]	*=x;
  wg22s[0]	=0.;

  wg22c[0]	=x*ESg22c[0];
  wg12c[0]	=ESg12c[0];
  k	=0;
  for(i=1; i < ES2Mp1; i++){
    k	+=ESNa1;
    splRA(wg12c+i,NULL,ESg12c+k,ESg12c2+k);
    splRA(wg22c+i,NULL,ESg22c+k,ESg22c2+k);
    wg22c[i]	*=x;
    
    wg12c[i]	=ESg12c[k];
    wg22c[i]	=x*ESg22c[k];

    wg12c[-i]	=wg12c[i];
    wg22c[-i]	=wg22c[i];
    splRA(wg12s-i,NULL,ESg12s+k,ESg12s2+k);
    splRA(wg22s-i,NULL,ESg22s+k,ESg22s2+k);
    wg22s[-i]	*=x;

    wg12s[-i]	=ESg12s[k];
    wg22s[-i]	=x*ESg22s[k];

    wg12s[i]	=-wg12s[-i];
    wg22s[i]	=-wg22s[-i];
  }

  sj	=ESjs[0];
  tyr	=EZcr2*ESsb1a[0];
  k	=ESNa1;
  r1r	= rcT1a[k];
  r1i	=-rsT1a[k];
  wuc[1]	=0.;
  wus[1]	=0.;
  wuc[-1]	=0.;
  wus[-1]	=0.;
  s		=EZcr2*x;
  wuc[0]	=s*rcT2a[0];
  wus[0]	=0.;
  for(i=2; i < ES2Mp1; i++){
    k	+=ESNa1;
    wuc[i]	= s*rcT2a[k];
    wus[i]	=-s*rsT2a[k];
    wuc[-i]	= wuc[i];
    wus[-i]	=-wus[i];
  }
  x1r	= r1r+tyr;	/*\tsx+i\tsy*/
  x1i	= r1i;
  X1r	= r1r-tyr;
  X1i	=-r1i;
  t1r	= r1r-tyr;
  t1i	= r1i;
  T1r	= r1r+tyr;
  T1i	=-r1i;
  /*}*/
  for(k=0; k < M0Mm; k++){
    gpr[k]	=0.;
    gpi[k]	=0.;
    dgpr[k]	=0.;
    dgpi[k]	=0.;
    Yr[k]	=0.;
    Yi[k]	=0.;
  }
  pr	=gpr	+ESMp;
  pi	=gpi	+ESMp;
  dpr	=dgpr	+ESMp;
  dpi	=dgpi	+ESMp;
  yr	=Yr	+ESMp;
  yi	=Yi	+ESMp;
  lr	=dYr	+ESMp;
  li	=dYi	+ESMp;
  
  /*m=0 polynomial solution{*/
  d	=ESaR0[0];
  f	=-EZcr4*x*sj;
  s	=f*(X1r*t1r-X1i*t1i+x1r*T1r-x1i*T1i);
  sr	=f*(x1r*t1r-x1i*t1i);
  si	=f*(x1r*t1i+x1i*t1r);
  
  pr[0]		=x*d*s;
  dpr[0]	=2.*d*s;
  f		=EZcr4*x;

  f	=0.;
  s		*=f;
  pr[1]		=x*s*r1r;
  dpr[1]	=3.*s*r1r;
  pi[1]		=x*s*r1i;
  dpi[1]	=3.*s*r1i;
  
  pr[2]		=x*d*sr;
  pi[2]		=x*d*si;
  dpr[2]	=2.*d*sr;
  dpi[2]	=2.*d*si;
  sr		*=f;
  si		*=f;
  if(3 < ESMp1){
    s	=sr*r1r-si*r1i;
    pr[3]	+=x*s;
    dpr[3]	+=3.*s;
    s	=si*r1r+sr*r1i;
    pi[3]	+=x*s;
    dpi[3]	+=3.*s;
  }
  s	=sr*r1r+si*r1i;
  pr[1]	+=x*s;
  dpr[1]+=3.*s;
  s	=si*r1r-sr*r1i;
  pi[1]	+=x*s;
  dpi[1]+=3.*s;

  f	=-0.5*x*x*d*sj;
  d	*=-1.5*x*sj;

  f	=0.;
  d	=0.;

  for(i=0; i < ESMp1; i++){
    j	=i-1;
    s	=r1r*wuc[j]-r1i*wus[j];
    pr[i]	+=f*s;
    dpr[i]	+=d*s;
    s	=sr*wus[j]+si*wuc[j];
    pi[i]	+=f*s;
    dpi[i]	+=d*s;
    j	=i+1;
    s	=r1r*wuc[j]+r1i*wus[j];
    pr[i]	+=f*s;
    dpr[i]	+=d*s;
    s	=sr*wus[j]-si*wuc[j];
    pi[i]	+=f*s;
    dpi[i]	+=d*s;
  }
  wus[0]=EZcr2*x*rsT2a[0];
  s	=-wus[0]*tyr;
  pi[1]	+=f*s;
  dpi[1]+=d*s;
  for(i=1; i < ESMp1; i++){
    pr[-i]	= pr[i];
    pi[-i]	=-pi[i];
    dpr[-i]	= dpr[i];
    dpi[-i]	=-dpi[i];
  }
  /*}*/

  s	=EZcr4/ESaR0[0];
  t1r	*=s;
  t1i	*=s;
  T1r	*=s;
  T1i	*=s;
  yr[0]	= 1.;
  yi[0]	= 0.;
  am	=1.;
  for(m=1; m < ESMp1; m++){/*$m\neq1$ solutions*/
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    am	*=x;
    for(k=-m+1; k < m; k+=2){
      sr	=yr[k-ES2Mp1];
      si	=yi[k-ES2Mp1];

      /*Second correction O(a);*/
#ifdef H
      f	=m*x;
      d	=m*(m+1);
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	s	=sr*wuc[j]-si*wus[j];
	pr[i]	+=f*s;
	dpr[i]	+=d*s;
	s	=si*wuc[j]+sr*wus[j];
	pi[i]	+=f*s;
	dpi[i]	+=d*s;
      }

#endif      
      /* Main order contribution;*/
      i	=k+1;
      s		=sr*x1r-si*x1i;
      dpr[i]	+=m*s;
      s		*=x;
      yr[i]	+=s;
      pr[i]	+=s;
      s		=si*x1r+sr*x1i;
      dpi[i]	+=m*s;
      s		*=x;
      yi[i]	+=s;
      pi[i]	+=s;
      i	=k-1;
      s		=sr*X1r-si*X1i;
      dpr[i]	+=m*s;
      s		*=x;
      yr[i]	+=s;
      pr[i]	+=s;
      s		=si*X1r+sr*X1i;
      dpi[i]	+=m*s;
      s		*=x;
      yi[i]	+=s;
      pi[i]	+=s;
    }

#ifdef H
    f	=1./m;
    d	=(m+1)*f;
    f	*=x;
    /*First correction $O(a)$*/
    for(k=-m; k <= m; k++){
      sr	=yr[k];
      si	=yi[k];
      i		=k+1;
      if(i < ESMp1){
	s	=t1r*sr-t1i*si;
	pr[i]	+=f*s;
	dpr[i]	+=d*s;
	s	=t1r*si+t1i*sr;
	pi[i]	+=f*s;
	dpi[i]	+=d*s;
      }
      i	=k-1;
      if(i > -ESMp1){
	s	=T1r*sr-T1i*si;
	pr[i]	+=f*s;
	dpr[i]	+=d*s;
	s	=T1r*si+T1i*sr;
	pi[i]	+=f*s;
	dpi[i]	+=d*s;
      }	
    }
#endif      

    /*Normalization;*/
    for(i=-ESMp; i < ESMp1; i++){
      lr[i]	= pr[-i];
      li[i]	=-pi[-i];
    }
    sr	=lr[-m];
    si	=li[-m];
    s	=1./(sr*sr+si*si);
    sr	*=s;
    s	*=-si;
    si	=pi[-m]*sr+pr[-m]*s;
    sr	=pr[-m]*sr-pi[-m]*s;
    for(i=-ESMp; i < ESMp1; i++){
      pr[i]	-=sr*lr[i]-si*li[i];
      pi[i]	-=si*lr[i]+sr*li[i];
    }
    for(i=-ESMp; i < ESMp1; i++){
      lr[i]	= dpr[-i];
      li[i]	=-dpi[-i];
    }
    for(i=-ESMp; i < ESMp1; i++){
      dpr[i]	-=sr*lr[i]-si*li[i];
      dpi[i]	-=si*lr[i]+sr*li[i];
    }
    if(m == 1){
      sr	=dpr[m];
      si	=dpi[m];
      s		=1./(sr*sr+si*si);
    }
    else{
      sr	=pr[m];
      si	=pi[m];
      s		=am/(sr*sr+si*si);
    }
    sr	*=s;
    si	*=-s;
    for(i=-ESMp; i < ESMp1; i++){
      s	=pr[i];
      pr[i]	=sr*s-si*pi[i];
      pi[i]	=si*s+sr*pi[i];
      s	=dpr[i];
      dpr[i]	=sr*s-si*dpi[i];
      dpi[i]	=si*s+sr*dpi[i];
    }
  }
  /*Normalization of $m=0$;*/
  pr	=gpr	+ESMp;
  pi	=gpi	+ESMp;
  dpr	=dgpr	+ESMp;
  dpi	=dgpi	+ESMp;
  i	=ESMp+2*ES2Mp1;
  yr	=gpr	+i;
  yi	=gpi	+i;
  lr	=dgpr	+i;
  li	=dgpi	+i;
  sr	=yr[2];
  si	=yi[2];
  s	=1./(sr*sr+si*si);
  sr	*=s;	/*$\R{1}{\gy^{*2}_2}$*/
  s	*=si;
  tr	=pr[-2]*sr-pi[-2]*s;	/*$\R{\gy^0_{-2}}{\gy^{*2}_2}$*/
  ti	=pi[-2]*sr+pr[-2]*s;
  si	=pi[2]*sr-pr[2]*s;
  sr	=pr[2]*sr+pi[2]*s;	/*$\R{\gy^0_2}{\gy^2_2}$*/
  for(i=-ESMp; i < ESMp1; i++){
    pr[i]	-=sr*yr[i]-si*yi[i]+tr*yr[-i]+ti*yi[-i];
    pi[i]	-=si*yr[i]+sr*yi[i]+ti*yr[-i]-tr*yi[-i];
    dpr[i]	-=sr*lr[i]-si*li[i]+tr*lr[-i]+ti*li[-i];
    dpi[i]	-=si*lr[i]+sr*li[i]+ti*lr[-i]-tr*li[-i];
  }
  s	=x*x;
  s	*=0.125*s*ESjp[0]*ESVc2[0]/wg22c[0];
  dpr[0]	-=s;
  pr[0]		-=EZcr4*s;

  yr	=Yr	+ESMp;
  yi	=Yi	+ESMp;
  for(m=0; m < ESMp1; m++){/*All $Y_m$*/
    for(i=-ESMp; i < ESMp; i++){
      yr[i]	=0.;
      yi[i]	=0.;
    }
    for(k=-ESMp; k < ESMp1; k++){
      sr	=dpr[k];
      si	=dpi[k];
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	yr[i]	+=sr*wg22c[j]-si*wg22s[j];
	yi[i]	+=si*wg22c[j]+sr*wg22s[j];
      }
      sr	= k*pi[k];
      si	=-k*pr[k];
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	yr[i]	+=sr*wg12c[j]-si*wg12s[j];
	yi[i]	+=si*wg12c[j]+sr*wg12s[j];
      }
    }
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
  }
  /*Averaged solution;*/
  ESaY[0]	=-EZcr2*sj*ESLc[0];
  ESdgY[0]	=ESaY[0]*x/wg22c[0];
  ESgY02a[0]	=ESdgY[0];
  s		=x*x;
  i		=M0Mm;
  Yi[i]		=s*(ESaY[0]-0.125*s*ESjp[0]*ESVc2[0]); /*$Y_0$*/
  i++;
  Yi[i]		=EZcr2*s*ESdgY[0]; /*$bgY_0$*/
  Yi[i]		=-EZcr4*s*x*(sj*ESLc[0]+0.125*s*ESjp[0]*ESVc2[0])/wg22c[0];
  gpr[ESMp]	=Yi[i];
  s	=2./s;
  k	=ESMp;
  return(0);
}

int Start2Dbb(double x)
{
  extern double *ESg12c1,*ESg12s1,*ESDPc1,*ESDPs1;

  int i,j,k,m,mm;
  int M,M1;
  double Bf[1024];
  double *bxr,*bxi,x1r,x1i,tyr,by0;
  double *Jr,*Ji;
  double *bKr,*bKi;
  double *bMr,*bMi;
  double *bNr,*bNi;
  double *bDr,*bDi,tD;
  double s,sr,si,tr,ti;
  double K0,Kr,Ki;
  double Cr,Ci;
  double *ur,*ui,*EZvr,*vi;
  double *pr,*pi,*dpr,*dpi;
  double *yr,*yi,*Rr,*Ri;
  double am,a1r,a1i,a2r,a2i;
  
  bxr	=Bf	+ESFp;
  k	=2*ESFp+1;
  bxi	=bxr	+k;
  bKr	=bxi	+k;
  k	=2*ES2Mp+1;
  bKi	=bKr	+k;
  bMr	=bKi	+k;
  bMi	=bMr	+k;
  bNr	=bMi	+k;
  bNi	=bNr	+k;
  bDr	=bNi	+k+2;
  bDi	=bDr	+k+4;
  Jr	=bDi	+k+4;
  Ji	=Jr	+ESMp1;
  EZvr	=Ji	+ES2Mp1;
  vi	=EZvr	+ES2Mp1;
  ur	=vi	+ES2Mp1;
  ui	=ur	+ES2Mp1;
  Rr	=ui	+ES2Mp1;
  Ri	=Rr	+ES2Mp1;

  /*Getting geometry near the origin{*/
  s	=EZcr2*x;
  tyr	=EZcr2*ESsb1a[0];
  bxr[0]=s*rcT2a[0];
  by0	=s*rsT2a[0];
  bxi[0]=0.;
  bxr[1]=0.;
  bxi[1]=0.;
  bxr[-1]=0.;
  bxi[-1]=0.;
  k	=ESNa1;
  x1r	= rcT1a[k];
  x1i	=-rsT1a[k];
  M	=2;
  for(i=2; i < ESFp1; i++){
    k		+=ESNa1;
    bxr[i]	= s*rcT2a[k];
    bxi[i]	=-s*rsT2a[k];

  bxr[i]	=0.;
  bxi[i]	=0.;

    bxr[-i]	= bxr[i];
    bxi[-i]	=-bxi[i];
    if(bxr[i] != 0. || bxi[i] != 0.) M=i;
  }

  bxr[0]=0.;
  by0	=0.;
#ifdef H
#endif

  M1	=M+1;
  tD	=4.*x1r*tyr;
  s	=1./tD;
  sr	=x1r*x1r;
  si	=x1i*x1i;
  tr	=tyr*tyr;
  K0	=2.*(tr+sr+si)*s;
  Kr	=(tr+si-sr)*s;
  Ki	=-2.*x1r*x1i*s;
  /*}*/

  /*Getting corrections to metric tensor{*/
  j	=ES2Mp1+2;
  if(j > ESFp1) j =ESFp1;
  for(i=0; i < j; i++){
    k	=i-1;
    s	=(2.-k)*tyr;
    bDr[i]	=s*bxr[k];
    bDi[i]	=s*bxi[k];
    k	=i+1;
    if(k < M1){
      s	=(2.+k)*tyr;
      bDr[i]	+=s*bxr[k];
      bDi[i]	+=s*bxi[k];
    }
  }

  bDi[0]	=0.;
  bDr[1]	+=2.*x1i*by0;
  bDi[1]	-=2.*x1r*by0;
  bDr[-1]	= bDr[1];
  bDi[-1]	=-bDi[1];
  bDr[-2]	= bDr[2];
  bDi[-2]	=-bDi[2];
  for(i=0; i < ES2Mp1; i++){
    bKr[i]	=-bDr[i]*K0;
    bKi[i]	=-bDi[i]*K0;
    bNr[i]	=bKr[i];
    bNi[i]	=bKi[i];
    k	=i-2;
    sr	=bDr[k];
    si	=bDi[k];
    tr	=sr*Kr-si*Ki;
    ti	=si*Kr+sr*Ki;
    bKr[i]	-=tr;
    bKi[i]	-=ti;
    bNr[i]	+=tr;
    bNi[i]	+=ti;
    bMr[i]	=-ti;
    bMi[i]	= tr;
    k	=i+2;
    sr	=bDr[k];
    si	=bDi[k];
    tr	=sr*Kr+si*Ki;
    ti	=si*Kr-sr*Ki;
    bKr[i]	-=tr;
    bKi[i]	-=ti;
    bNr[i]	+=tr;
    bNi[i]	+=ti;
    bMr[i]	-=ti;
    bMi[i]	+=tr;
    k	=i-1;
    s	=-2.*k;
    sr	=s*bxr[k];
    si	=s*bxi[k];
    bKr[i]	+=sr*x1r-si*x1i;
    bKi[i]	+=si*x1r+sr*x1i;
    s	=k+2.;
    sr	=-s*bxi[k];
    si	= s*bxr[k];
    bMr[i]	+=sr*x1r-si*x1i;
    bMi[i]	+=si*x1r+sr*x1i;
    sr	=4.*bxr[k];
    si	=4.*bxi[k];
    bNr[i]	+=sr*x1r-si*x1i;
    bNi[i]	+=si*x1r+sr*x1i;
    k	=i+1;
    if(k < M1){
      bNr[i]	+=sr*x1r+si*x1i;
      bNi[i]	+=si*x1r-sr*x1i;
      s	=2.*k;
      sr	=s*bxr[k];
      si	=s*bxi[k];
      bKr[i]	+=sr*x1r+si*x1i;
      bKi[i]	+=si*x1r-sr*x1i;
      s	=k-2.;
      sr	=-s*bxi[k];
      si	= s*bxr[k];
      bMr[i]	+=sr*x1r+si*x1i;
      bMi[i]	+=si*x1r-sr*x1i;
    }
  }
  bMr[1]	+=2.*tyr*by0;
  bNi[1]	-=4.*tyr*by0;
  s	=1./tD;
  for(k=0; k < ES2Mp1; k++){
    bKr[k]	*=s;
    bMr[k]	*=s;
    bNr[k]	*=s;
    if(k){
      bKi[k]	*=s;
      bMi[k]	*=s;
      bNi[k]	*=s;
      i	=-k;
      bDr[i]	= bDr[k];
      bDi[i]	=-bDi[k];
      bKr[i]	= bKr[k];
      bKi[i]	=-bKi[k];
      bMr[i]	= bMr[k];
      bMi[i]	=-bMi[k];
      bNr[i]	= bNr[k];
      bNi[i]	=-bNi[k];
    }
  }
  bKi[0]	=0.;
  bMi[0]	=0.;
  bNi[0]	=0.;
  /*}*/

  for(k=0; k < M0Mm; k++){
    gpr[k]	=0.;
    gpi[k]	=0.;
    dgpr[k]	=0.;
    dgpi[k]	=0.;
    Yr[k]	=0.;
    Yi[k]	=0.;
    dYr[k]	=0.;
    dYi[k]	=0.;
  }
  /*Main order solution{*/
  pr	=gpr+ESMp;
  pi	=gpi+ESMp;
  dpr	=dgpr+ESMp;
  dpi	=dgpi+ESMp;
  yr	=Yr+ESMp;
  yi	=Yi+ESMp;
  yr[0]	=-EZcr4*ESaR0[0]*tD*ESjs[0]/K0;
  pr[0]	=yr[0];
  dpr[0]=2.*yr[0];
  yr	+=ES2Mp1;
  yi	+=ES2Mp1;
  pr	+=ES2Mp1;
  pi	+=ES2Mp1;
  dpr	+=ES2Mp1;
  dpi	+=ES2Mp1;
  yr[1]	=1.;
  pr[1]	=yr[1];
  dpr[1]=yr[1];

  am	=x*x;
  for(m=2; m < ESMp1; m++){
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    i	=-m;
    EZvr[i]	=0.;
    vi[i]	=0.;
    for(k=i+2; k < m; k+=2,i+=2){
      s=(m-k)*(m-k+2);
      Cr=s*Kr;
      Ci=s*Ki;
      sr=(m*m-k*k)*K0+Cr*EZvr[i]-Ci*vi[i];
      si=Ci*EZvr[i]+Cr*vi[i];
      s=1./(sr*sr+si*si);
      sr*= s;
      si*=-s;
      tr=Cr*yr[i]-Ci*yi[i];
      ti=Ci*yr[i]+Cr*yi[i];
      yr[k]	=-tr*sr+ti*si;
      yi[k]	=-ti*sr-tr*si;
      s	=(m+k)*(m+k+2);
      tr= s*Kr;
      ti=-s*Ki;
      EZvr[k]	=-tr*sr+ti*si;
      vi[k]	=-ti*sr-tr*si;
    }
    i	=m;
    yr[i]=1.;
    yi[i]=0.;
    j	=2-m;
    while(i > j){
      k	=i-2;
      yr[k]	+=EZvr[k]*yr[i]-vi[k]*yi[i];
      yi[k]	+=vi[k]*yr[i]+EZvr[k]*yi[i];
      i	-=2;
    }
    for(i=-m+2; i <= m; i+=2){
      pr[i]	=yr[i];
      pi[i]	=yi[i];
      dpr[i]	=m*yr[i];
      dpi[i]	=m*yi[i];
    }
    am	*=x;
  }  
  /*}*/

  /*Corrections to the EZmain order{*/
  pr	=gpr+ESMp;
  pi	=gpi+ESMp;
  dpr	=dgpr+ESMp;
  dpi	=dgpi+ESMp;
  yr	=Yr+ESMp;
  yi	=Yi+ESMp;
  am	=1.;
  for(m=0; m < ESMp1; m++){/*All fundamental vectors*/
    am	*=x;
    if(m == 1) am	=1.;
    /*RHS;*/
    for(i=-ESMp; i < ESMp1; i++){
      Rr[i]	=0.;
      Ri[i]	=0.;
      for(j=2-m; j <= m; j+=2){
	k	=i-j;
	s	=j*(m+1)+m*i;
	sr	=-m*(m+1)*bKr[k]-s*bMi[k]+i*j*bNr[k];
	si	=-m*(m+1)*bKi[k]+s*bMr[k]+i*j*bNi[k];
	Rr[i]	+=sr*pr[j]-si*pi[j];
	Ri[i]	+=si*pr[j]+sr*pi[j];
      }
      j	=i-1;
      if(0 && 2-m <= j && j <= m){
	sr=x*(m-j)*tyr/ESaR0[0];
	Rr[i]	+=sr*pr[j];
	Ri[i]	+=sr*pi[j];
      }
      j	=i+1;
      if(0 && 2-m <= j && j <= m){
	sr=x*(m+j)*tyr/ESaR0[0];
	Rr[i]	+=sr*pr[j];
	Ri[i]	+=sr*pi[j];
      }
    }
    mm	=m ? m+1 : 3;
    /*Elimination of resonant harmonic in the RHS;*/
    if(mm < ESMp1){ 
      if(m){
	ur	=yr+ES2Mp1;
	ui	=yi+ES2Mp1;
      }
      else{
	ur	=yr+3*ES2Mp1;
	ui	=yi+3*ES2Mp1;
      }
      k	=mm-2;
      sr	=ur[k];
      si	=ui[k];
      k		+=2;
      tr	=mm*K0*ur[k]+Kr*sr-Ki*si;
      ti	=mm*K0*ui[k]+Ki*sr+Kr*si;
      k	=2-mm;
      sr	=Kr*ur[k]+Ki*ui[k];
      si	=Ki*ur[k]-Kr*ui[k];
      s		=EZcr2/(tr*tr+ti*ti-sr*sr-si*si);
      sr	*=-s;
      si	*=-s;
      tr	*=s;
      ti	*=s;
      Cr	=Rr[-mm];
      Ci	=Ri[-mm];
      a1r	=Rr[mm]*tr+Ri[mm]*ti+Cr*sr+Ci*si;
      a1i	=Ri[mm]*tr-Rr[mm]*ti+Ci*sr-Cr*si;
      a2r	=Rr[mm]*tr-Ri[mm]*ti+Cr*sr-Ci*si;
      a2i	=Ri[mm]*tr+Rr[mm]*ti+Ci*sr+Cr*si;
      s	=log(x/ESsa[1]);
      for(i=-mm; i <= mm; i+=2){
	sr	=a1r*ur[i]-a1i*ui[i]+a2r*ur[-i]+a2i*ui[-i];
	si	=a1i*ur[i]+a1r*ui[i]+a2i*ur[-i]-a2r*ui[-i];
	dpr[i]	-=(mm*s+1.)*sr;
	dpi[i]	-=(mm*s+1.)*si;
	pr[i]	-=s*sr;
	pi[i]	-=s*si;
	tr	=2.*mm*K0;
	Rr[i]	+=tr*sr;
	Ri[i]	+=tr*sr;
	k	=i+2;
	if(k < ESMp1){
	  tr	=2.*(mm+1-k)*Kr;
	  ti	=2.*(mm+1-k)*Ki;
	  Rr[k]	+=tr*sr-ti*si;
	  Ri[k]	+=ti*sr+tr*si;
	}
	k	=i-2;
	if(k > -ESMp1){
	  tr	=2.*(mm+1+k)*Kr;
	  ti	=2.*(mm+1+k)*Ki;
	  Rr[k]	+=tr*sr+ti*si;
	  Ri[k]	+=ti*sr-tr*si;
	}
      }
    }

    /*Making corrections with modified RHS;*/
    ur	=vi	+ES2Mp1;
    ui	=ur	+ES2Mp1;
    i	=-ESMp;
    for(k=-ESMp; k < ESMp1; k++){
      if(k != -mm && k != mm){
	i	=k-2;
	sr	=(mm*mm-k*k)*K0;
	si	=0.;
	if(k+ESMp > 1){
	  s=(mm-k)*(mm-k+2);
	  Cr=s*Kr;
	  Ci=s*Ki;
	  sr	+=Cr*EZvr[i]-Ci*vi[i];
	  si	+=Ci*EZvr[i]+Cr*vi[i];
	  s	=1./(sr*sr+si*si);
	  sr*= s;
	  si*=-s;
	  tr=Rr[k]+Cr*ur[i]-Ci*ui[i];
	  ti=Ri[k]+Ci*ur[i]+Cr*ui[i];
	  ur[k]	=-tr*sr+ti*si;
	  ui[k]	=-ti*sr-tr*si;
	}
	else{
	  sr	=1./sr;
	  ur[k]	=-Rr[k]*sr;
	  ui[k]	=-Ri[k]*sr;
	}
	s	=(mm+k)*(mm+k+2);
	tr= s*Kr;
	ti=-s*Ki;
	EZvr[k]	=-tr*sr+ti*si;
	vi[k]	=-ti*sr-tr*si;
      }
      else{
	ur[k]	=0.;
	ui[k]	=0.;
	EZvr[k]	=0.;
	vi[k]	=0.;
      }
    }
    k	=ESMp;
    EZvr[k]	=0.;
    vi[k]	=0.;
    k--;
    EZvr[k]	=0.;
    vi[k]	=0.;
    k	=ESMp1;
    while(k > -ESMp){
      k--;
      i	=k+2;
      ur[k]	+=EZvr[k]*ur[i]-vi[k]*ui[i];
      ui[k]	+=vi[k]*ur[i]+EZvr[k]*ui[i];
      pr[k]	+=ur[k];
      pi[k]	+=ui[k];
      dpr[k]	+=mm*ur[k];
      dpi[k]	+=mm*ui[k];
      pr[k]	*=x*am;
      pi[k]	*=x*am;
      dpr[k]	*=am;
      dpi[k]	*=am;
    }
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
  }
  /*}*/

  /*Calculation of $\vaY${*/

#ifdef HH
  /*Getting metric tensor{*/
  ESSetSplA(x);
  splRA(wg12c,NULL,ESg12c,ESg12c2);
  wg12s[0]	=0.;
  splRA(wg22c,NULL,ESg22c,ESg22c2);
  wg22s[0]	=0.;

  wuc[0]	=ESDPc[0]+x*ESDPc1[0];
  wus[0]	=0.;
  wg22c[0]	=ESg22c[0]+x*ESg22c1[0];
  wg12c[0]	=ESg12c[0]+x*ESg12c1[0];

  k	=0;
  for(i=1; i < ES2Mp1; i++){
    k	+=ESNa1;
    splRA(wg12c+i,NULL,ESg12c+k,ESg12c2+k);
    splRA(wg22c+i,NULL,ESg22c+k,ESg22c2+k);
    wg12c[-i]	=wg12c[i];
    wg22c[-i]	=wg22c[i];
    splRA(wg12s-i,NULL,ESg12s+k,ESg12s2+k);
    splRA(wg22s-i,NULL,ESg22s+k,ESg22s2+k);
    wg12s[i]	=-wg12s[-i];
    wg22s[i]	=-wg22s[-i];

    wuc[i]	=ESDPc[k]+x*ESDPc1[k];
    wus[i]	=-ESDPs[k]-x*ESDPs1[k];
    wg22c[i]	=ESg22c[k]+x*ESg22c1[k];
    wg12c[i]	=ESg12c[k]+x*ESg12c1[k];
    wg22s[i]	=-ESg22s[k]-x*ESg22s1[k];
    wg12s[i]	=-ESg12s[k]-x*ESg12s1[k];
  }

  {
    double Rc[65],Rs[65],Rac[65],Ras[65],b,cs,sn;
    double r[128],ra[128],rt[128],xx[128],xt[128]
      ,za[128],zt[128],D[128],G[128],g[128];

    /* Getting metric coefficients{*/

    b	=ESsb1a[0];
    Rc[0]	=EZcr2*x*rcT2a[0];
    Rs[0]	=EZcr2*x*rsT2a[0];
    Rac[0]	=x*rcT2a[0];
    Ras[0]	=x*rsT2a[0];

    m	=0;
    M1	=2;
    for(k=1; k < ESFp1; k++){
      m	+=ESNa1;
      Rc[k]	=k > 1 ? x*rcT2a[m] : 2.*rcT1a[m];
      Rs[k]	=k > 1 ? x*rsT2a[m] : 2.*rsT1a[m];
      Rac[k]	=k > 1 ? 2.*x*rcT2a[m] : 2.*rcT1a[m];
      Ras[k]	=k > 1 ? 2.*x*rsT2a[m] : 2.*rsT1a[m];
      if(Rc[k] != 0. || Rs[k] != 0.) M1=k+1;
    }
    for(j=0; j < ESNp; j++){
      xx[j]	=Rc[0];
      xt[j]	=0.;
      m	=0;
      for(k=1; k < ESFp1; k++){
	m	+=j;
	if(m >= ESNp) m -=ESNp;
	cs	=EScs1[m];
	sn	=ESsn1[m];
	xx[j]	+=Rc[k]*cs+Rs[k]*sn;
	xt[j]	+=k*(Rs[k]*cs-Rc[k]*sn);
      }
      cs	=EScs1[j];
      sn	=ESsn1[j];

      za[j]	=b*sn;
      zt[j]	=b*cs;

      r[j]	=-x*(Rc[1]*cs+Rs[1]*sn)/ESaR0[0];

#ifdef H
      r		=(1.+r[j])/ESaR0[0];
      ra	=ra[j]+2.*xx[j];
      rt	=rt[j]+xt[j];
      D		=tD+D[j];
      G		=rt[j]*rt[j]+zt[j]*zt[j]+G[j];
      g		=ra[j]*rt[j]+za[j]*zt[j]+g[j];
      d		=(1.+D[j])/tD;
#endif

      ra[j]	=Rc[1]*cs+Rs[1]*sn;
      rt[j]	=Rs[1]*cs-Rc[1]*sn;

      D[j]	=-((2.*xx[j]*cs-xt[j]*sn)*b-rt[j]*Ras[0])/tD;

      G[j]	=(2.*rt[j]*xt[j]+
		  (rt[j]*rt[j]+zt[j]*zt[j])*(r[j]+D[j]))/(tD*ESaR0[0]);
      g[j]	=(ra[j]*xt[j]+2.*xx[j]*rt[j]+Ras[0]*b*cs+
		  (ra[j]*rt[j]+za[j]*zt[j])*(r[j]+D[j]))/(tD*ESaR0[0]);
    }
    ESP2F(wg22c,wg22s,G,ES2Mp);
    ESP2F(wg12c,wg12s,g,ES2Mp);
    for(k=1; k < ES2Mp1; k++){
      wg12s[k]	=-wg12s[k];
      wg22s[k]	=-wg22s[k];
      EZout("sidddddd","M",k,bMr[k],wg12c[k],bMi[k],wg12s[k]
	  ,bMr[k]-wg12c[k],bMi[k]-wg12s[k]);
    }
  }
#endif

  bKr[0]	+=K0;
  bKr[2]	+=Kr;
  bKi[2]	+=Ki;
  bMr[2]	+=Ki;
  bMi[2]	-=Kr;

  sr	=x*x1r/ESaR0[0];
  si	=x*x1i/ESaR0[0];
  sr	=0.;
  si	=0.;
#ifdef H
#endif
  bKr[3]	-=Kr*sr-Ki*si;
  bKi[3]	-=Ki*sr+Kr*si;
  bKr[1]	-=K0*sr+Kr*sr+Ki*si;
  bKi[1]	-=K0*si+Ki*sr-Kr*si;
  bMr[3]	-=Ki*sr+Kr*si;
  bMi[3]	-=Ki*si-Kr*sr;
  bMr[1]	-=Ki*sr-Kr*si;
  bMi[1]	+=Ki*si+Kr*sr;
  s	=1./ESaR0[0];
  bDr[0]	=tD;
  for(k=0; k < ES2Mp1; k++){
    bKr[k]	*=s;
    bKi[k]	*=s;
    bMr[k]	*=s;
    bMi[k]	*=s;
#ifdef H
    EZout("sidddddd","M",k,bMr[k],wg12c[k],bMi[k],wg12s[k]
	,bMr[k]-wg12c[k],bMi[k]-wg12s[k]);
    EZout("sidddddd","K",k,bKr[k],wg22c[k],bKi[k],wg22s[k]
	,bKr[k]-wg22c[k],bKi[k]-wg22s[k]);
    EZout("sidddddd","D",k,bDr[k],wuc[k],bDi[k],wus[k]
	,bDr[k]-wuc[k],bDi[k]-wus[k]);

    bKr[k]	=wg22c[k];
    bKi[k]	=wg22s[k];
    bMr[k]	=wg12c[k];
    bMi[k]	=wg12s[k];
#endif
    i	=-k;
    bKr[i]	= bKr[k];
    bKi[i]	=-bKi[k];
    bMr[i]	= bMr[k];
    bMi[i]	=-bMi[k];
  }

  pr	=gpr	+ESMp;
  pi	=gpi	+ESMp;
  dpr	=dgpr	+ESMp;
  dpi	=dgpi	+ESMp;
  yr	=Yr	+ESMp;
  yi	=Yi	+ESMp;
  for(m=0; m < ESMp1; m++){/*All $Y^m$*/
    for(i=-ESMp; i < ESMp1; i++){
      yr[i]	=0.;
      yi[i]	=0.;
    }
    for(k=-ESMp; k < ESMp1; k++){
      sr	=x*dpr[k];
      si	=x*dpi[k];
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	yr[i]	+=sr*bKr[j]-si*bKi[j];
	yi[i]	+=si*bKr[j]+sr*bKi[j];
      }
      sr	= k*pi[k];
      si	=-k*pr[k];
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	yr[i]	+=sr*bMr[j]-si*bMi[j];
	yi[i]	+=si*bMr[j]+sr*bMi[j];
      }
    }
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
  }
  /*}*/

  /*Averaged solution;*/
  ESdgY[0]	=-EZcr2*ESaR0[0]*tD*ESjs[0]/K0;
  ESaY[0]	=K0*ESdgY[0];
  ESgY02a[0]	=ESdgY[0];
  s		=x*x;
  i		=M0Mm;
  Yi[i]		=s*(ESaY[0]-0.125*s*ESjp[0]*ESVc2[0]);
  i++;
  Yi[i]		=EZcr2*x*x*ESdgY[0];
  return(0);
}

int Start2D(double x)
{
  extern double *ESg12c1,*ESg12s1,*ESDPc1,*ESDPs1;

  int i,j,k,m,mm;
  double Bf[1024];
  double x1r,x1i,tyr;
  double *bKr,*bKi;
  double *bMr,*bMi;
  double tD;
  double s,sr,si,tr,ti;
  double K0,Kr,Ki;
  double Cr,Ci;
  double *ur,*ui,*EZvr,*vi;
  double *pr,*pi,*dpr,*dpi;
  double *yr,*yi;
  double am,a1r,a1i,a2r,a2i;
  
  k	=2*ES2Mp+1;
  bKr	=Bf	+ES2Mp;
  bKi	=bKr	+k;
  bMr	=bKi	+k;
  bMi	=bMr	+k;
  EZvr	=bMi	+k;
  vi	=EZvr	+ES2Mp1;
  ur	=vi	+ES2Mp1;
  ui	=ur	+ES2Mp1;

  /*Getting geometry near the origin{*/
  tyr	=EZcr2*ESsb1a[0];
  k	=ESNa1;
  x1r	= rcT1a[k];
  x1i	=-rsT1a[k];
  tD	=4.*x1r*tyr;
  s	=1./tD;
  sr	=x1r*x1r;
  si	=x1i*x1i;
  tr	=tyr*tyr;
  K0	=2.*(tr+sr+si)*s;
  Kr	=(tr+si-sr)*s;
  Ki	=-2.*x1r*x1i*s;
  /*}*/
  for(k=0; k < M0Mm; k++){
    gpr[k]	=0.;
    gpi[k]	=0.;
    dgpr[k]	=0.;
    dgpi[k]	=0.;
    Yr[k]	=0.;
    Yi[k]	=0.;
    dYr[k]	=0.;
    dYi[k]	=0.;
  }
  /*Main order solution{*/
  pr	=gpr+ESMp;
  pi	=gpi+ESMp;
  dpr	=dgpr+ESMp;
  dpi	=dgpi+ESMp;
  yr	=Yr+ESMp;
  yi	=Yi+ESMp;
  yr[0]	=-EZcr4*ESaR0[0]*tD*ESjs[0]/K0;
  pr[0]	=yr[0];
  dpr[0]=2.*yr[0];
  yr	+=ES2Mp1;
  yi	+=ES2Mp1;
  pr	+=ES2Mp1;
  pi	+=ES2Mp1;
  dpr	+=ES2Mp1;
  dpi	+=ES2Mp1;
  yr[1]	=1.;
  pr[1]	=yr[1];
  dpr[1]=yr[1];

  am	=x*x;
  for(m=2; m < ESMp1; m++){
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    i	=-m;
    EZvr[i]	=0.;
    vi[i]	=0.;
    for(k=i+2; k < m; k+=2,i+=2){
      s=(m-k)*(m-k+2);
      Cr=s*Kr;
      Ci=s*Ki;
      sr=(m*m-k*k)*K0+Cr*EZvr[i]-Ci*vi[i];
      si=Ci*EZvr[i]+Cr*vi[i];
      s=1./(sr*sr+si*si);
      sr*= s;
      si*=-s;
      tr=Cr*yr[i]-Ci*yi[i];
      ti=Ci*yr[i]+Cr*yi[i];
      yr[k]	=-tr*sr+ti*si;
      yi[k]	=-ti*sr-tr*si;
      s	=(m+k)*(m+k+2);
      tr= s*Kr;
      ti=-s*Ki;
      EZvr[k]	=-tr*sr+ti*si;
      vi[k]	=-ti*sr-tr*si;
    }
    i	=m;
    yr[i]=1.;
    yi[i]=0.;
    j	=2-m;
    while(i > j){
      k	=i-2;
      yr[k]	+=EZvr[k]*yr[i]-vi[k]*yi[i];
      yi[k]	+=vi[k]*yr[i]+EZvr[k]*yi[i];
      i	-=2;
    }
    for(i=-m+2; i <= m; i+=2){
      pr[i]	=yr[i];
      pi[i]	=yi[i];
      dpr[i]	=m*yr[i];
      dpi[i]	=m*yi[i];
    }
    am	*=x;
  }  
  /*}*/

  /*Corrections to the EZmain order{*/
  pr	=gpr+ESMp;
  pi	=gpi+ESMp;
  dpr	=dgpr+ESMp;
  dpi	=dgpi+ESMp;
  am	=1.;
  for(m=0; m < ESMp1; m++){/*All fundamental vectors*/
    am	*=x;
    if(m == 1) am	=1.;
    /*Making corrections with modified RHS;*/
    for(k=-ESMp; k < ESMp1; k++){
      pr[k]	*=x*am;
      pi[k]	*=x*am;
      dpr[k]	*=am;
      dpi[k]	*=am;
    }
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
  }
  /*}*/

  /*Calculation of $\vaY${*/
  for(k=0; k < ES2Mp1; k++){
    bKr[k]	=0.;
    bKi[k]	=0.;
    bMr[k]	=0.;
    bMi[k]	=0.;
  }
  s	=1./ESaR0[0];
  bKr[0]	= K0*s;
  bKr[2]	= Kr*s;
  bKi[2]	= Ki*s;
  bMr[2]	= Ki*s;
  bMi[2]	=-Kr*s;
  for(k=0; k < ES2Mp1; k++){
    i	=-k;
    bKr[i]	= bKr[k];
    bKi[i]	=-bKi[k];
    bMr[i]	= bMr[k];
    bMi[i]	=-bMi[k];
  }

  pr	=gpr	+ESMp;
  pi	=gpi	+ESMp;
  dpr	=dgpr	+ESMp;
  dpi	=dgpi	+ESMp;
  yr	=Yr	+ESMp;
  yi	=Yi	+ESMp;
  for(m=0; m < ESMp1; m++){/*All $Y^m$*/
    for(i=-ESMp; i < ESMp1; i++){
      yr[i]	=0.;
      yi[i]	=0.;
    }
    for(k=-ESMp; k < ESMp1; k++){
      sr	=x*dpr[k];
      si	=x*dpi[k];
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	yr[i]	+=sr*bKr[j]-si*bKi[j];
	yi[i]	+=si*bKr[j]+sr*bKi[j];
      }
      sr	= k*pi[k];
      si	=-k*pr[k];
      for(i=-ESMp; i < ESMp1; i++){
	j	=i-k;
	yr[i]	+=sr*bMr[j]-si*bMi[j];
	yi[i]	+=si*bMr[j]+sr*bMi[j];
      }
    }
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
  }
  /*}*/

  /*Averaged solution;*/
  ESdgY[0]	=-EZcr2*ESaR0[0]*tD*ESjs[0]/K0;
  ESaY[0]	=K0*ESdgY[0];
  ESgY02a[0]	=ESdgY[0];
  s		=x*x;
  i		=M0Mm;
  Yi[i]		=s*(ESaY[0]-0.125*s*ESjp[0]*ESVc2[0]);
  i++;
  Yi[i]		=EZcr2*x*x*ESdgY[0];
  return(0);
}

#endif/*stg_2DInit*/

#ifndef stg_2DAdvR

static double *mUr,*mUi;
static double wj_s,dj_s,wj_p,dj_p;

void eq2DVsj(int *neq,double *x,double *gyr,double *dgyr)
{
  double *pr,*pi,*yr,*yi,*dpr,*dpi,*dyr,*dyi;
  double *xr,*xi,*pU;
  int i,j,k,m;

  double sr,si,tr,ti,Tr,Ti,dK0;

  /* Getting metric coefficients{*/
  int M1;
  double Rc[65],Rs[65],Rac[65],Ras[65],b,ba,baa,cs,sn;
  double r[128],ra[128],rt[128],za[128],zt[128],D[128],G[128];

  /* Getting metric coefficients{*/
  ESSetSplA(*x);
  splRA(&b,&ba,ESsb,ESsb2a);
  splRA(Rc,Rac,rcT,rcT2a);
  splRA(Rs,Ras,rsT,rsT2a);
  m	=0;
  M1	=2;
  for(k=1; k < ESFp1; k++){
    m	+=ESNa1;
    splRA(Rc+k,Rac+k,rcT+m,rcT2a+m);
    splRA(Rs+k,Ras+k,rsT+m,rsT2a+m);
    Rc[k]	*=2.;
    Rs[k]	*=2.;
    Rac[k]	*=2.;
    Ras[k]	*=2.;
    if(Rac[k] != 0. || Ras[k] != 0.) M1=k+1;
  }
  for(j=0; j < ESNp; j++){
    r[j]	=ESaR0[0]+Rc[0];
    ra[j]	=Rac[0];
    rt[j]	=0.;
    za[j]	=Ras[0]+ba*ESsn1[j];
    zt[j]	=b*EScs1[j];
    m	=0;
    for(k=1; k < ESFp1; k++){
      m	+=j;
      if(m >= ESNp) m -=ESNp;
      cs	=EScs1[m];
      sn	=ESsn1[m];
      r[j]	+=Rc[k]*cs+Rs[k]*sn;
      ra[j]	+=Rac[k]*cs+Ras[k]*sn;
      rt[j]	+=k*(Rs[k]*cs-Rc[k]*sn);
    }
    D[j]	=ra[j]*zt[j]-rt[j]*za[j];
    G[j]	=(rt[j]*rt[j]+zt[j]*zt[j])/(r[j]*D[j]);
  }
  ESP2F(wg22c,wg22s,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=(ra[j]*rt[j]+za[j]*zt[j])/(r[j]*D[j]);
  }
  ESP2F(wg12c,wg12s,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=(ra[j]*ra[j]+za[j]*za[j])/(r[j]*D[j]);
  }
  ESP2F(wg11c,wg11s,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=D[j]*ESaR0[0]/r[j];
  }
  ESP2F(wuc,wus,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=D[j]*(r[j]/ESaR0[0]-ESaR0[0]/r[j]);
  }
  ESP2F(wvc,wvs,G,ES2Mp);
  for(k=1; k < ES2Mp1; k++){
    i	=-k;
    wg11c[i]	=wg11c[k];
    wg12c[i]	=wg12c[k];
    wg22c[i]	=wg22c[k];
    wg11s[i]	=wg11s[k];
    wg12s[i]	=wg12s[k];
    wg22s[i]	=wg22s[k];
    
    wuc[i]	=wuc[k];
    wus[i]	=wus[k];
    wvc[i]	=wvc[k];
    wvs[i]	=wvs[k];

    wg11s[k]	=-wg11s[i];
    wg12s[k]	=-wg12s[i];
    wg22s[k]	=-wg22s[i];
    wvs[k]	=-wvs[i];
    wus[k]	=-wus[i];
  }
  /*}*/

#ifdef H
  ESSetSplA(*x);
  sr	=1./(*x);
  m	=0;
  splRA(wg11c,NULL,ESg11c+m,ESg11c2+m);
  wg11c[0]	*=sr;
  wg11s[0]	=0.;
  splRA(wg12c,NULL,ESg12c+m,ESg12c2+m);
  wg12s[0]	=0.;
  splRA(wg22c,&dK0,ESg22c+m,ESg22c2+m);
  dK0	=wg22c[0]+(*x)*dK0;
  wg22c[0]	*=*x;
  wg22s[0]	=0.;
  splRA(wuc,NULL,ESLc+m,ESLc2+m);
  wuc[0]	*=*x;
  wus[0]	=0.;
  splRA(wvc,NULL,ESVc+m,ESVc2+m);
  wvc[0]	*=*x;
  wvs[0]	=0.;
  for(k=1; k < ES2Mp1; k++){
    m	+=ESNa1;
    splRA(wg11c+k,NULL,ESg11c+m,ESg11c2+m);
    wg11c[k]	*=sr;
    wg11c[-k]	=wg11c[k];
    splRA(wg12c+k,NULL,ESg12c+m,ESg12c2+m);
    wg12c[-k]	=wg12c[k];
    splRA(wg22c+k,NULL,ESg22c+m,ESg22c2+m);
    wg22c[k]	*=*x;
    wg22c[-k]	=wg22c[k];
    splRA(wuc+k,NULL,ESLc+m,ESLc2+m);
    wuc[k]	*=*x;
    wuc[-k]	=wuc[k];
    splRA(wvc+k,NULL,ESVc+m,ESVc2+m);
    wvc[k]	*=*x;
    wvc[-k]	=wvc[k];
    splRA(wg11s-k,NULL,ESg11s+m,ESg11s2+m);
    wg11s[-k]	*=sr;
    wg11s[k]	=-wg11s[-k];
    splRA(wg12s-k,NULL,ESg12s+m,ESg12s2+m);
    wg12s[k]	=-wg12s[-k];
    splRA(wg22s-k,NULL,ESg22s+m,ESg22s2+m);
    wg22s[-k]	*=*x;
    wg22s[k]	=-wg22s[-k];
    splRA(wus-k,NULL,ESLs+m,ESLs2+m);
    wus[-k]	*=*x;
    wus[k]	=-wus[-k];
    splRA(wvs-k,NULL,ESVs+m,ESVs2+m);
    wvs[-k]	*=*x;
    wvs[k]	=wvs[-k];
  }
#endif

  /*}*/
  /* Getting 1D profiles{*/
  k	=*neq-2;
  Ti	=wg22c[0]/gyr[k];
  ESSetSplDPr(*x);
  ESSetSplDCr(*x);

  switch(ESEqSolvInCr){
  case 0:
    splRDCr(&wj_s,&dj_s,ESEqSolvInCr);
    break;
  case 1:
  case 2:
    splRA(&wj_s,&dj_s,ESjs,ESjs2a);
    break;
  case 3:
    splRDCr(&wj_s,&dj_s,ESEqSolvInCr);
    wj_s	*=rR0;
    dj_s	*=rR0;
    break;
  }

  dj_s	*=Ti;
  switch(ESEqSolvInPr){
  case 0:
    splRDPr(&wj_p,&dj_p,0);
    dj_p	*=Ti;
    break;
  case 1:
    splRDPr(&wj_p,&dj_p,1);
    wj_p	*=R0;
    dj_p	*=R0*Ti;
    break;
  case 2:
    splRDPr2(&si,&wj_p,&dj_p,2);
    Tr	=R0*Ti;
    wj_p	*=Tr;
    dj_p	=(dj_p*Tr+wj_p*(dK0+(wj_s*wuc[0]+wj_p*wvc[0])*Ti)/wg22c[0])*Ti;
    break;
  }
  if(ESEqSolvInCr == 3){
    wj_s	+=wj_p;
    dj_s	+=dj_p;
  }

#ifdef H
  if(ECEqIterFl == 0){
    dj_s	*=-2;
    dj_p	*=-2;
  }
#endif
  dj_s	=0.;
  dj_p	=0.;

  k		=*neq-2;
  dgyr[k+1]	=gyr[k]/wg22c[0];
  dgyr[k]	=-(wj_s*wuc[0]+wj_p*wvc[0]);

  mUr	=Ur+ES2Mp;
  mUi	=Ui+ES2Mp;
  mUr[0]	=-dj_s*wuc[0]-dj_p*wvc[0];
  mUi[0]	=0.;
  for(k=1; k < ES2Mp1; k++){
    mUr[k]	=-dj_s*wuc[k]-dj_p*wvc[k];
    mUi[k]	=-dj_s*wus[k]-dj_p*wvs[k];
    mUr[-k]	= mUr[k];
    mUi[-k]	=-mUi[k];
  }
  /*}*/
  /* Cholesky decomposition of $\maK${*/
  for(i=0; i < ES2Mp1; i++){
    k	=ES2Mp1*i+i;
    for(j=i; j < ES2Mp1; j++){
      Kr[k]	=wg22c[i-j];
      Ki[k]	=wg22s[i-j];
      k++;
    }
  }
  if((k=CholDeComp(Kr,Ki,ES2Mp1))){
#ifdef XWIN
    extern char *CbUserWarning;
    extern char ESMessage[];
    sprintf(ESMessage,"CholDeComp() failed\nat eq2DVsj()\n a=%5.3f i=%d"
	    ,*x,k-1);
    CbUserWarning	=ESMessage;
#endif
#ifdef DEBUG
    printf("??? x=%10.3e\n",*x);
    for(m=0; m < ES2Mp1; m++){
      printf("??? g22[%2d]=%10.3e %10.3e\n",m,wg22c[m],wg22s[m]);
    }
#endif
    EcCholFl=1;
    return;
  }
  /* Calculation of normalized solution with $bY_0=1$;*/
  pr	=gyr+ESMp;
  pi	=pr+M0Mm;
  yr	=pi+M0Mm;
  yi	=yr+M0Mm;
  dpr	=dgyr+ESMp;
  dpi	=dpr+M0Mm;
  dyr	=dpi+M0Mm;
  dyi	=dyr+M0Mm;
  for(m=0; m < ESMp1; m++){
    /* Calculation of $dgy${*/
    for(i=-ESMp; i < ESMp1; i++){
      dpr[i]	=yr[i];
      dpi[i]	=yi[i];
      for(j=-ESMp; j < ESMp1; j++){
	sr	=j*pr[j];
	si	=j*pi[j];
	k	=i-j;
	tr	=-wg12s[k];
	ti	= wg12c[k];
	dpr[i]	+=sr*tr-si*ti;
	dpi[i]	+=si*tr+sr*ti;
      }
      dyr[i]	=0.;
      dyi[i]	=0.;
    }
    CholBksbComp(Kr,Ki,ES2Mp1,dpr-ESMp,dpi-ESMp);
    /*}*/
    /* Calculation of $daY${*/
    for(j=-ESMp; j < ESMp1; j++){
      sr	=dpr[j];
      si	=dpi[j];
      for(i=-ESMp; i < ESMp1; i++){
	k	=i-j;
	tr	=-i*wg12s[k];
	ti	= i*wg12c[k];
	dyr[i]	+=sr*tr-si*ti;
	dyi[i]	+=si*tr+sr*ti;
      }
      if(j){
	sr	=pr[j];
	si	=pi[j];
	for(i=-ESMp; i < ESMp1; i++){
	  k	=i-j;
	  tr	=mUr[k]+i*wg11c[k]*j;
	  ti	=mUi[k]+i*wg11s[k]*j;
	  dyr[i]	+=sr*tr-si*ti;
	  dyi[i]	+=si*tr+sr*ti;
	}
      }
    }
    /*}*/
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;

    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    dyr	+=ES2Mp1;
    dyi	+=ES2Mp1;
  }
  dyr	-=M0Mm;
  dyi	-=M0Mm;
  for(i=-ESMp; i < ESMp1; i++){
    dyr[i]	-=wj_s*wuc[i]+wj_p*wvc[i];
    dyi[i]	-=wj_s*wus[i]+wj_p*wvs[i];
  }
}

void eq2DVjb(int *neq,double *x,double *gyr,double *dgyr)
{
  extern double *ESsr,*ESsra,*ESsz,*ESsza;

  double *pr,*pi,*yr,*yi,*dpr,*dpi,*dyr,*dyi;
  double *mUr,*mUi;
  int i,j,k,m,M1;
  double sr,si,tr,ti,Ti,wj_s,dj_s,jb,djb,bJ,dJ,wj_p,dj_p,dK0,dv0;

  double Rc[65],Rs[65],Rac[65],Ras[65],Raac[65],Raas[65],b,ba,baa,cs,sn;
  double r[128],ra[128],raa[128],rt[128],rat[128]
    ,za[128],zaa,zt[128],zat[128],D[128],Da[128],G[128];

  /* Getting metric coefficients{*/
  ESSetSplA(*x);

  splRA2(&b,&ba,&baa,ESsb,ESsb2a);
  splRA2(Rc,Rac,Raac,rcT,rcT2a);
  splRA2(Rs,Ras,Raas,rsT,rsT2a);
  m	=0;
  M1	=2;
  for(k=1; k < ESFp1; k++){
    m	+=ESNa1;
    splRA2(Rc+k,Rac+k,Raac+k,rcT+m,rcT2a+m);
    splRA2(Rs+k,Ras+k,Raas+k,rsT+m,rsT2a+m);
    Rc[k]	*=2.;
    Rs[k]	*=2.;
    Rac[k]	*=2.;
    Ras[k]	*=2.;
    Raac[k]	*=2.;
    Raas[k]	*=2.;
    if(Rc[k] != 0. || Rs[k] != 0.) M1=k+1;
  }
  for(j=0; j < ESNp; j++){
    r[j]	=ESaR0[0]+Rc[0];
    ra[j]	=Rac[0];
    raa[j]	=Raac[0];
    rt[j]	=0.;
    rat[j]	=0.;
    za[j]	=Ras[0]+ba*ESsn1[j];
    zaa		=Raas[0]+baa*ESsn1[j];
    zt[j]	=b*EScs1[j];
    zat[j]	=ba*EScs1[j];
    m	=0;
    for(k=1; k < M1; k++){
      m	+=j;
      if(m >= ESNp) m -=ESNp;
      cs	=EScs1[m];
      sn	=ESsn1[m];
      r[j]	+=Rc[k]*cs+Rs[k]*sn;
      ra[j]	+=Rac[k]*cs+Ras[k]*sn;
      raa[j]	+=Raac[k]*cs+Raas[k]*sn;
      rt[j]	+=k*(Rs[k]*cs-Rc[k]*sn);
      rat[j]	+=k*(Ras[k]*cs-Rac[k]*sn);
    }
    D[j]	=ra[j]*zt[j]-rt[j]*za[j];
    G[j]	=(rt[j]*rt[j]+zt[j]*zt[j])/(r[j]*D[j]);
    Da[j]	=raa[j]*zt[j]+ra[j]*zat[j]-rat[j]*za[j]-rt[j]*zaa;
  }
  ESP2F(wg22c,wg22s,G,ES2Mp);
  dK0	=0.;
  for(j=0; j < ESNp; j++){
    sr	=1./(r[j]*D[j]);
    dK0	+=(2.*rt[j]*rat[j]+2.*zt[j]*zat[j]
	   -G[j]*(Da[j]*r[j]+D[j]*ra[j]))*sr;
    G[j]	=(ra[j]*rt[j]+za[j]*zt[j])*sr;
  }
  dK0	/=ESNp;
  for(j=0; j < ESNp; j++){
    G[j]	=(ra[j]*rt[j]+za[j]*zt[j])/(r[j]*D[j]);
  }
  ESP2F(wg12c,wg12s,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=(ra[j]*ra[j]+za[j]*za[j])/(r[j]*D[j]);
  }
  ESP2F(wg11c,wg11s,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=D[j]*ESaR0[0]/r[j];
  }
  ESP2F(wuc,wus,G,ES2Mp);
  wduc[0]	=0.;
  for(j=0; j < ESNp; j++){
    wduc[0]	+=(Da[j]*ESaR0[0]-G[j]*ra[j])/r[j];
    G[j]	=D[j]*(r[j]/ESaR0[0]-ESaR0[0]/r[j]);
  }
  wduc[0]	/=ESNp;
  ESP2F(wvc,wvs,G,ES2Mp);
  dv0	=0.;
  for(j=0; j < ESNp; j++){
    dv0	+=Da[j]*(r[j]/ESaR0[0]-ESaR0[0]/r[j])
      +D[j]*(r[j]/ESaR0[0]+ESaR0[0]/r[j])*ra[j]/r[j];
  }
  dv0	/=ESNp;
  if(0 && fabs(*x-ESsa[36]) < 1e-6){
    m	=36;
    for(k=0; k < ES2Mp1; k++){
      EZout("siddddd","==?? Vc",k,ESVc[m]*(*x),wvc[k]
	  ,ESVs[m]*(*x),wvs[k],*x);
      m+=ESNa1;
    }
  }
  for(k=1; k < ES2Mp1; k++){
    wg11c[-k]	=wg11c[k];
    wg12c[-k]	=wg12c[k];
    wg22c[-k]	=wg22c[k];

    wg11s[-k]	=wg11s[k];
    wg12s[-k]	=wg12s[k];
    wg22s[-k]	=wg22s[k];
    
    wuc[-k]	=wuc[k];
    wvc[-k]	=wvc[k];
    wus[-k]	=wus[k];
    wvs[-k]	=wvs[k];

    wg11s[k]	=-wg11s[k];
    wg12s[k]	=-wg12s[k];
    wg22s[k]	=-wg22s[k];
    wus[k]	=-wus[k];
    wvs[k]	=-wvs[k];
  }
  sr	=1./(*x);
#ifdef H
  m	=0;
  splRA(wg11c,NULL,ESg11c+m,ESg11c2+m);
  wg11c[0]	*=sr;
  wg11s[0]	=0.;
  splRA(wg12c,NULL,ESg12c+m,ESg12c2+m);
  wg12s[0]	=0.;
  splRA(wg22c,&dK0,ESg22c+m,ESg22c2+m);
  dK0		=wg22c[0]+(*x)*dK0;
  wg22c[0]	*=*x;
  wg22s[0]	=0.;
  splRA(wuc,wduc,ESLc+m,ESLc2+m);
  wduc[0]	=wuc[0]+(*x)*wduc[0];
  wuc[0]	*=*x;
  wdus[0]	=0.;
  wus[0]	=0.;
  splRA(wvc,&dv0,ESVc+m,ESVc2+m);
  dv0	=wvc[0]+(*x)*dv0;
  wvc[0]	*=*x;
  wvs[0]	=0.;
  for(k=1; k < ES2Mp1; k++){
    m	+=ESNa1;
    splRA(wg11c+k,NULL,ESg11c+m,ESg11c2+m);
    wg11c[k]	*=sr;
    wg11c[-k]	=wg11c[k];
    splRA(wg12c+k,NULL,ESg12c+m,ESg12c2+m);
    wg12c[-k]	=wg12c[k];
    splRA(wg22c+k,NULL,ESg22c+m,ESg22c2+m);
    wg22c[k]	*=*x;
    wg22c[-k]	=wg22c[k];
    splRA(wuc+k,NULL,ESLc+m,ESLc2+m);
    wduc[k]	=wuc[k]+(*x)*wduc[k];
    wduc[-k]	=wduc[k];
    wuc[k]	*=*x;
    wuc[-k]	=wuc[k];
    splRA(wvc+k,wduc+k,ESVc+m,ESVc2+m);
    wvc[k]	*=*x;
    wvc[-k]	=wvc[k];
    splRA(wg11s-k,NULL,ESg11s+m,ESg11s2+m);
    wg11s[-k]	*=sr;
    wg11s[k]	=-wg11s[-k];
    splRA(wg12s-k,NULL,ESg12s+m,ESg12s2+m);
    wg12s[k]	=-wg12s[-k];
    splRA(wg22s-k,NULL,ESg22s+m,ESg22s2+m);
    wg22s[-k]	*=*x;
    wg22s[k]	=-wg22s[-k];
    splRA(wus-k,wdus-k,ESLs+m,ESLs2+m);
    wdus[-k]	=wus[-k]+(*x)*wdus[-k];
    wus[-k]	*=*x;
    wdus[k]	=-wdus[-k];
    wus[k]	=-wus[-k];
    splRA(wvs-k,NULL,ESVs+m,ESVs2+m);
    wvs[-k]	*=*x;
    wvs[k]	=wvs[-k];
  }
  /*}*/
#endif

  /* Getting 1D profiles{*/
  k	=*neq-2;
  bJ	=-gyr[k];
  Ti	=wg22c[0]/gyr[k];
  dgyr[k+1]	=gyr[k]/wg22c[0];
  ESSetSplDPr(*x);
  ESSetSplDCr(*x);
  splRDCr(&jb,&djb,ESEqSolvInCr);
  splRA(&sr,&si,ESaF,ESaF2a);
  if(ESEqSolvInCr == 4){
    ti	=1./(wuc[0]*sr);
    tr	=R0*(wuc[0]+wvc[0])*ti;
    djb	=tr*djb+(R0*(wduc[0]+dv0)-tr*(wuc[0]*si+wduc[0]*sr))*ti*jb;
    jb	*=tr;
  }
  if(ESEqSolvInCr == 2){
    djb	*=rR0;
    jb	*=rR0;
  }
  bJ		/=sr;
  dJ		=wuc[0]*jb/sr;
  dgyr[k]	=-(wuc[0]*jb+bJ*si);

#ifdef H
  dgyr[k]	=-(1.+dF/F)*wuc[0]*jb;
#endif

  switch(ESEqSolvInPr){
  case 0:
    splRDPr(&wj_p,&dj_p,0);
    break;
  case 1:
    splRDPr(&wj_p,&dj_p,1);
    wj_p	*=R0;
    dj_p	*=R0;
    break;
  case 2:
    splRDPr2(&ti,&wj_p,&dj_p,2);
    wj_p	*=R0*Ti;
    dj_p	=R0*Ti*dj_p+wj_p*(dK0-dgyr[k]*Ti)/wg22c[0];
    break;
  }
  ti	=1/(wuc[0]*wg22c[0]);
  si	=R0*ti*bJ;
  sr	=si*bJ;
  si	=2.*dJ*si-sr*(wduc[0]*wg22c[0]+wuc[0]*dK0)*ti;
  tr	=1./(1.+sr);
  dj_s	=1./wuc[0];
  ti	=sr-wvc[0]*dj_s;
  dj_s	=si+(wvc[0]*wduc[0]*dj_s-dv0)*dj_s;
  wj_s	=(jb+ti*wj_p)*tr;
  dj_s	=(djb+ti*dj_p+dj_s*wj_p-wj_s*si)*tr;
  dj_s	*=Ti;
  dj_p	*=Ti;

/*?????????????????*/
  if(0 && ECEqIterFl == 0){
    dj_s	=0.;
    dj_p	=0.;
  }
/*?????????????????*/
  mUr	=Ur+ES2Mp;
  mUi	=Ui+ES2Mp;
  mUr[0]	=-dj_s*wuc[0]-dj_p*wvc[0];
  mUi[0]	=0.;
  for(k=1; k < ES2Mp1; k++){
    mUr[k]	=-dj_s*wuc[k]-dj_p*wvc[k];
    mUi[k]	=-dj_s*wus[k]-dj_p*wvs[k];
    mUr[-k]	= mUr[k];
    mUi[-k]	=-mUi[k];
  }
  /*}*/
  /* Cholesky decomposition of $\maK${*/
  for(i=0; i < ES2Mp1; i++){
    k	=ES2Mp1*i+i;
    for(j=i; j < ES2Mp1; j++){
      Kr[k]	=wg22c[i-j];
      Ki[k]	=wg22s[i-j];
      k++;
    }
  }
  if((k=CholDeComp(Kr,Ki,ES2Mp1))){
#ifdef XWIN
    extern char *CbUserWarning;
    extern char ESMessage[];
    sprintf(ESMessage,"CholDeComp() failed\nat eq2DVjb()\n a=%5.3f i=%d"
	    ,*x,k-1);
    CbUserWarning	=ESMessage;
#endif
#ifdef DEBUG
    printf("??? x=%10.3e\n",*x);
    for(j=0; j < ESNp1; j++){
      sr	=wg22c[0];
      k	=0;
      for(m=0; m < ES2Mp1; m++){
	k	+=j;
	if(k >= ESNp) k	-=ESNp;
	sr	+=wg22c[m]*EScs1[k]+wg22s[m]*ESsn1[k];
      }
      printf("?--- g22[%2d]=%10.3e\n",j,sr);
    }
    printf("??? x=%10.3e\n",*x);
    for(m=0; m < ES2Mp1; m++){
      printf("?+?? g22[%2d]=%10.3e %10.3e\n",m,wg22c[m],wg22s[m]);
    }
#endif
    EcCholFl=1;
    return;
  }
  /*}*/

  pr	=gyr+ESMp;
  pi	=pr+M0Mm;
  yr	=pi+M0Mm;
  yi	=yr+M0Mm;
  dpr	=dgyr+ESMp;
  dpi	=dpr+M0Mm;
  dyr	=dpi+M0Mm;
  dyi	=dyr+M0Mm;
  for(m=0; m < ESMp1; m++){
    /* Calculation of $dgy${*/
    for(i=-ESMp; i < ESMp1; i++){
      dpr[i]	=yr[i];
      dpi[i]	=yi[i];
      for(j=-ESMp; j < ESMp1; j++){
	sr	=-j*pi[j];
	si	= j*pr[j];
	k	=i-j;
	tr	=wg12c[k];
	ti	=wg12s[k];
	dpr[i]	+=sr*tr-si*ti;
	dpi[i]	+=si*tr+sr*ti;
      }
      dyr[i]	=0.;
      dyi[i]	=0.;
    }
    CholBksbComp(Kr,Ki,ES2Mp1,dpr-ESMp,dpi-ESMp);
    /*}*/
    /* Calculation of $daY${*/
    for(j=-ESMp; j < ESMp1; j++){
      sr	=dpr[j];
      si	=dpi[j];
      for(i=-ESMp; i < ESMp1; i++){
	k	=i-j;
	tr	=-i*wg12s[k];
	ti	= i*wg12c[k];
	dyr[i]	+=sr*tr-si*ti;
	dyi[i]	+=si*tr+sr*ti;
      }
      if(j){
	sr	=pr[j];
	si	=pi[j];
	for(i=-ESMp; i < ESMp1; i++){
	  k	=i-j;
	  tr	=mUr[k]+i*wg11c[k]*j;
	  ti	=mUi[k]+i*wg11s[k]*j;
	  dyr[i]	+=sr*tr-si*ti;
	  dyi[i]	+=si*tr+sr*ti;
	}
      }
    }
    /*}*/
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    dyr	+=ES2Mp1;
    dyi	+=ES2Mp1;
  }
  dyr	-=M0Mm;
  dyi	-=M0Mm;
  for(i=-ESMp; i < ESMp1; i++){
    dyr[i]	-=wj_s*wuc[i]+wj_p*wvc[i];
    dyi[i]	-=wj_s*wus[i]+wj_p*wvs[i];
  }
  return;
}

void eq2DVgm(int *neq,double *x,double *gyr,double *dgyr)
{
  double *pr,*pi,*yr,*yi,*dpr,*dpi,*dyr,*dyi;
  double *xr,*xi,*pU;

  int mp,mm,m,k,kp,km,K,n;
  double sr,si,tr,ti,Tr,Ti,wj_s,dj_s,wj_p,dj_p,c0bLxdgy;
  double gm,dgm,bL,dL,bF,dF,bV,bK,dK;


  int M1,i,j;
  double Rc[65],Rs[65],Rac[65],Ras[65],Raac[65],Raas[65],b,ba,baa,cs,sn;
  double r[128],ra[128],raa,rt[128],rat[128]
    ,za[128],zaa,zt[128],zat[128],D[128],Da[128],G[128];


  /* Getting metric coefficients{*/

  ESSetSplA(*x);
  splRA2(&b,&ba,&baa,ESsb,ESsb2a);
  splRA2(Rc,Rac,Raac,rcT,rcT2a);
  splRA2(Rs,Ras,Raas,rsT,rsT2a);
  m	=0;
  M1	=2;
  for(k=1; k < ESFp1; k++){
    m	+=ESNa1;
    splRA2(Rc+k,Rac+k,Raac+k,rcT+m,rcT2a+m);
    splRA2(Rs+k,Ras+k,Raas+k,rsT+m,rsT2a+m);
    Rc[k]	*=2.;
    Rs[k]	*=2.;
    Rac[k]	*=2.;
    Ras[k]	*=2.;
    Raac[k]	*=2.;
    Raas[k]	*=2.;
    if(Rac[k] != 0. || Ras[k] != 0.) M1=k+1;
  }
  for(j=0; j < ESNp; j++){
    r[j]	=ESaR0[0]+Rc[0];
    ra[j]	=Rac[0];
    raa		=Raac[0];
    rt[j]	=0.;
    rat[j]	=0.;
    cs		=EScs1[j];
    sn		=ESsn1[j];
    za[j]	=Ras[0]+ba*sn;
    zaa		=Raas[0]+baa*sn;
    zt[j]	=b*cs;
    zat[j]	=ba*cs;
    m	=0;
    for(k=1; k < M1; k++){
      m		+=j;
      if(m >= ESNp) m -=ESNp;
      cs	=EScs1[m];
      sn	=ESsn1[m];
      r[j]	+=Rc[k]*cs+Rs[k]*sn;
      ra[j]	+=Rac[k]*cs+Ras[k]*sn;
      raa	+=Raac[k]*cs+Raas[k]*sn;
      rt[j]	+=k*(Rs[k]*cs-Rc[k]*sn);
      rat[j]	+=k*(Ras[k]*cs-Rac[k]*sn);
    }
    D[j]	=ra[j]*zt[j]-rt[j]*za[j];
    G[j]	=(rt[j]*rt[j]+zt[j]*zt[j])/(r[j]*D[j]);
    Da[j]	=raa*zt[j]+ra[j]*zat[j]-rat[j]*za[j]-rt[j]*zaa;
  }
  ESP2F(wg22c,wg22s,G,ES2Mp);
  dK	=0.;
  for(j=0; j < ESNp; j++){
    sr	=1./(r[j]*D[j]);
    dK	+=(2.*rt[j]*rat[j]+2.*zt[j]*zat[j]
	   -G[j]*(Da[j]*r[j]+D[j]*ra[j]))*sr;
    G[j]	=(ra[j]*rt[j]+za[j]*zt[j])*sr;
  }
  dK	/=ESNp;
  ESP2F(wg12c,wg12s,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=(ra[j]*ra[j]+za[j]*za[j])/(r[j]*D[j]);
  }
  ESP2F(wg11c,wg11s,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=D[j]*ESaR0[0]/r[j];
  }
  ESP2F(wuc,wus,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=(Da[j]*ESaR0[0]-G[j]*ra[j])/r[j];
  }
  ESP2F(wduc,wdus,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=D[j]*(r[j]/ESaR0[0]-ESaR0[0]/r[j]);
  }
  ESP2F(wvc,wvs,G,ES2Mp);
  for(k=1; k < ES2Mp1; k++){
    i	=-k;
    wg11c[i]	=wg11c[k];
    wg12c[i]	=wg12c[k];
    wg22c[i]	=wg22c[k];
    wg11s[i]	=wg11s[k];
    wg12s[i]	=wg12s[k];
    wg22s[i]	=wg22s[k];
    
    wuc[i]	=wuc[k];
    wus[i]	=wus[k];
    wduc[i]	=wduc[k];
    wdus[i]	=wdus[k];
    wvc[i]	=wvc[k];
    wvs[i]	=wvs[k];

#ifdef H
    wg11s[k]	=-wg11s[i];
    wg12s[k]	=-wg12s[i];
    wg22s[k]	=-wg22s[i];
    wvs[k]	=-wvs[i];
    wus[k]	=-wus[i];
    wdus[k]	=-wdus[i];
#endif
  }
  bK		=wg22c[0];
  bL		=wuc[0];
  dL		=wduc[0];
  bV		=wvc[0];

#ifdef H
  n	=0;
  K	=ESnMp*n;
  ESSetSplA(*x);
  sr	=1./(*x);
  m	=K;
  splRA(wg11c,NULL,ESg11c+m,ESg11c2+m);
  wg11c[0]	*=sr;
  splRA(wg12c,NULL,ESg12c+m,ESg12c2+m);
  splRA(wg22c,&dK,ESg22c+m,ESg22c2+m);
  dK		=wg22c[0]+(*x)*dK;
  wg22c[0]	*=*x;
  bK		=wg22c[0];
  splRA(wuc,wduc,ESLc+m,ESLc2+m);
  wduc[0]	=wuc[0]+(*x)*wduc[0];
  wuc[0]	*=*x;
  bL		=wuc[0];
  dL		=wduc[0];
  splRA(wvc,NULL,ESVc+m,ESVc2+m);
  wvc[0]	*=*x;
  bV		=wvc[0];
  for(k=1; k < ES2Mp1; k++){
    m	+=ESNa1;
    splRA(wg11c+k,NULL,ESg11c+m,ESg11c2+m);
    wg11c[k]	*=sr;
    splRA(wg12c+k,NULL,ESg12c+m,ESg12c2+m);
    splRA(wg22c+k,NULL,ESg22c+m,ESg22c2+m);
    wg22c[k]	*=*x;
    splRA(wuc+k,wduc+k,ESLc+m,ESLc2+m);
    wduc[k]	=wuc[k]+(*x)*wduc[k];
    wuc[k]	*=*x;
    splRA(wvc+k,NULL,ESVc+m,ESVc2+m);
    wvc[k]	*=*x;
    splRA(wg11s+k,NULL,ESg11s+m,ESg11s2+m);
    wg11s[k]	*=sr;
    splRA(wg12s+k,NULL,ESg12s+m,ESg12s2+m);
    splRA(wg22s+k,NULL,ESg22s+m,ESg22s2+m);
    wg22s[k]	*=*x;
    splRA(wus+k,wdus+k,ESLs+m,ESLs2+m);
    wdus[k]	=wus[k]+(*x)*wdus[k];
    wus[k]	*=*x;
    splRA(wvs+k,NULL,ESVs+m,ESVs2+m);
    wvs[k]	*=*x;
  }
#endif
  /*}*/

  /* Getting 1D profiles{*/
  ESSetSplDPr(*x);
  ESSetSplDCr(*x);
  k	=*neq-1;
  splRA(&bF,&dF,ESaF,ESaF2a);
  splRDCr(&gm,&dgm,7);
  dgyr[k]	=-bL*gm*bF*rR0;
  Ti		=gm != 0. ? 1./dgyr[k] : 0.;
  switch(ESEqSolvInPr){
  case 0:
    splRDPr(&wj_p,&dj_p,0);
    break;
  case 1:
    splRDPr(&wj_p,&dj_p,1);
    wj_p	*=R0;
    dj_p	*=R0;
    break;
  case 2:
    splRDPr2(&si,&wj_p,&dj_p,2);
    break;
  }
  if(ESEqSolvInPr == 2){
    Tr		=R0*Ti;
    wj_p	*=Tr;
    wj_s	=-dF/(gm*wuc[0])+wj_p;
    dj_p	=dj_p*Tr+wj_p*(dK+(wj_s*wuc[0]+wj_p*wvc[0])*Ti)/wg22c[0];
  }
  else{
    sr	=bK*bL*gm*gm*rR0;
    si	=1./(1.+sr);
    tr	=((dK*bL+bK*dL)*gm+bK*bL*dgm)*rR0*si;
    ti	=wj_p*(bL+bV)*si;
    wj_s=wj_p-(ti-tr*bF)/bL;
  }
  splRA(&tr,&dj_s,ESjs,ESjs2a);
  dj_p	*=Ti;
  dj_s	*=Ti;
  k--;
  dgyr[k]	=-(wj_s*wuc[0]+wj_p*wvc[0]);

  sr		=1./wuc[0];
  si		=-wduc[0];
  dUr[0]	=-dj_s*wuc[0]-dj_p*wvc[0];
  for(k=1; k < ES2Mp1; k++){
    dUr[k]	=-dj_s*wuc[k]-dj_p*wvc[k];
    dUi[k]	=dj_s*wus[k]+dj_p*wvs[k];
    wuc[k]	*=sr;
    wus[k]	*=sr;
    wduc[k]	=(wduc[k]+wuc[k]*si)*sr;
    wdus[k]	=(wdus[k]+wus[k]*si)*sr;
  }
  /*}*/
  
  /* Cholesky decomposition of $\BbK${*/
  kp		=0;
  k		=0;
  ac[k]		=wg22c[0];
  for(mm=1; mm < ESMp1; mm++){
    k++;
    ac[k]	=2.*wg22c[mm];
    k++;
    ac[k]	=-2.*wg22s[mm];
  }

  for(m=1; m < ESMp1; m++){
    kp++;
    k	=ES2Mp1*kp+kp;
    mp	=m+m;
    km	=0;
    mm	=ESMp1-m;
    while(km < mm){
      ac[k]	=2.*(wg22c[km]+wg22c[mp]);
      k++;
      ac[k]	=-2.*(wg22s[km]+wg22s[mp]);
      k++;
      km++;
      mp++;
    }
    kp++;
    k	=ES2Mp1*kp+kp;
    mm	=ESMp-m;
    mp	=m+m;
    km	=0;
    ac[k]	=2.*(wg22c[km]-wg22c[mp]);
    while(km < mm){
      k++;
      km++;
      mp++;
      ac[k]	=2.*(wg22s[km]-wg22s[mp]);
      k++;
      ac[k]	=2.*(wg22c[km]-wg22c[mp]);
    }
  }
  if(CholDc(ac,ES2Mp1) == -1){
#ifdef XWIN
    printf("??? x=%10.3e\n",*x);
    for(m=0; m < ES2Mp1; m++){
      printf("??? g22[%2d]=%10.3e %10.3e\n",m,wg22c[m],wg22s[m]);
    }
#endif
    EcCholFl=1;
    return;
  }
  /* Calculation of normalized solution with $bY_0=1$;*/
  pU	=Ur+ES2Mp1;
  pU[0]	=1.;
  k	=1;
  mp	=0;
  while(k < ES2Mp1){
    mp++;
    pU[k]	=2.*wuc[mp];
    k++;
    pU[k]	=-2.*wus[mp];
    k++;
  }
  CholBksb(ac,ES2Mp1,pU);
  k	=1;
  kp	=2;
  mp	=0;
  mm	=0;
  c0bLxdgy	=0.;
  while(k < ES2Mp1){
    mp++;
    mm--;
    c0bLxdgy	+=wuc[mp]*pU[k]-wus[mp]*pU[kp];
    k	+=2;
    kp	+=2;
  }
  c0bLxdgy	=1./(2.*c0bLxdgy+pU[0]);
  /*}*/

  xr	=Xr+ESMp;
  xi	=Xi+ESMp;
  pr	=gyr+ESMp;
  pi	=pr+M0Mm;
  yr	=pi+M0Mm;
  yi	=yr+M0Mm;
  dpr	=dgyr+ESMp;
  dpi	=dpr+M0Mm;
  dyr	=dpi+M0Mm;
  dyi	=dyr+M0Mm;
  for(m=0; m < ESMp1; m++){
    /* Calculation of $dgy${*/
    /* Contribution from $tbY$;*/
    xr[0]	=0.;
    xi[0]	=0.;
    sr	=wg12c[0];
    mm	=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      si	=sr*mp;
      xr[mp]	=yr[mp]-si*pi[mp];
      xi[mp]	=yi[mp]+si*pr[mp];
      xr[mm]	=yr[mm]+si*pi[mm];
      xi[mm]	=yi[mm]-si*pr[mm];
    }
    mm	=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      tr	=mp*pr[mp];
      ti	=mp*pi[mp];
      Tr	=mp*pr[mm];
      Ti	=mp*pi[mm];
      k		=1;
      kp	=mp+k;
      km	=mm-k;
      while(kp < ESMp1){
	sr	=wg12s[k];
	si	=wg12c[k];
	k++;
	xr[kp]	+=sr*tr-si*ti;
	xi[kp]	+=sr*ti+si*tr;
	kp++;
	xr[km]	+=sr*Tr+si*Ti;
	xi[km]	+=sr*Ti-si*Tr;
	km--;
      }
      k		=1;
      kp	=k+mm;
      km	=mp-k;
      while(kp < ESMp1){
	sr	=wg12s[k];
	si	=wg12c[k];
	k++;
	xr[kp]	-=sr*Tr-si*Ti;
	xi[kp]	-=sr*Ti+si*Tr;
	kp++;
	xr[km]	-=sr*tr+si*ti;
	xi[km]	-=sr*ti-si*tr;
	km--;
      }
    }
    Ur[0]	=xr[0];
    Ui[0]	=xi[0];
    k	=1;
    mp	=0;
    mm	=0;
    while(k < ES2Mp1){
      mp++;
      mm--;
      Ur[k]	=xr[mp]+xr[mm];
      Ui[k]	=xi[mp]+xi[mm];
      k++;
      Ur[k]	=xi[mp]-xi[mm];
      Ui[k]	=xr[mm]-xr[mp];
      k++;
    }
    CholBksb(ac,ES2Mp1,Ur);
    CholBksb(ac,ES2Mp1,Ui);
    /* Contribution from $bY_0$;*/
    k	=1;
    kp	=2;
    mp	=0;
    mm	=0;
    Tr	=Ur[0];
    Ti	=Ui[0];
    while(k < ES2Mp1){
      mp++;
      mm--;
      Tr	+=wduc[mp]*(pr[mm]+pr[mp])+wdus[mp]*(pi[mm]-pi[mp])
	+2.*(wuc[mp]*Ur[k]-wus[mp]*Ur[kp]);
      Ti	+=wduc[mp]*(pi[mm]+pi[mp])+wdus[mp]*(pr[mp]-pr[mm])
	+2.*(wuc[mp]*Ui[k]-wus[mp]*Ui[kp]);
      k		+=2;
      kp	+=2;
    }
    yr[0]	=-Tr*c0bLxdgy;
    yi[0]	=-Ti*c0bLxdgy;
    if(m == 0){
      yr[0]	+=dgyr[*neq-1]*c0bLxdgy;
    }
    k	=0;
    while(k < ES2Mp1){
      Ur[k]	+=yr[0]*pU[k];
      Ui[k]	+=yi[0]*pU[k];
      k++;
    }
    dpr[0]	=Ur[0];
    dpi[0]	=Ui[0];
    k	=1;
    kp	=2;
    mp	=0;
    mm	=0;
    while(k < ES2Mp1){
      mp++;
      mm--;
      dpr[mp]	=Ur[k]-Ui[kp];
      dpi[mp]	=Ur[kp]+Ui[k];
      dpr[mm]	=Ur[k]+Ui[kp];
      dpi[mm]	=Ui[k]-Ur[kp];
      k		+=2;
      kp	+=2;
    }
    /*}*/
    /* Calculation of $dbY${*/
    tr	=pr[0];
    ti	=pi[0];
    sr	=dUr[0];
    dyr[0]	=tr*sr;
    dyi[0]	=ti*sr;
    mm	=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      dyr[mp]	=pr[mp]*sr;
      dyi[mp]	=pi[mp]*sr;
      dyr[mm]	=pr[mm]*sr;
      dyi[mm]	=pi[mm]*sr;
    }

    tr	=pr[0];
    ti	=pi[0];
    km	=0;
    for(kp=1; kp < ESMp1; kp++){
      km--;
      sr	=dUr[kp];
      si	=dUi[kp];
      dyr[kp]	+=sr*tr-si*ti;
      dyi[kp]	+=sr*ti+si*tr;
      dyr[km]	+=sr*tr+si*ti;
      dyi[km]	+=sr*ti-si*tr;
    }
    mm	=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      tr	=pr[mp];
      ti	=pi[mp];
      Tr	=pr[mm];
      Ti	=pi[mm];
      k		=1;
      kp	=mp+k;
      km	=mm-k;
      while(kp < ESMp1){
	sr	=dUr[k];
	si	=dUi[k];
	k++;
	dyr[kp]	+=sr*tr-si*ti;
	dyi[kp]	+=sr*ti+si*tr;
	kp++;
	dyr[km]	+=sr*Tr+si*Ti;
	dyi[km]	+=sr*Ti-si*Tr;
	km--;
      }
      k		=1;
      kp	=k+mm;
      km	=mp-k;
      while(kp < ESMp1){
	sr	=dUr[k];
	si	=dUi[k];
	k++;
	dyr[kp]	+=sr*Tr-si*Ti;
	dyi[kp]	+=sr*Ti+si*Tr;
	kp++;
	dyr[km]	+=sr*tr+si*ti;
	dyi[km]	+=sr*ti-si*tr;
	km--;
      }
    }

    sr	=wg11c[0];
    mm	=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      k		=mp*mp;
      dyr[mp]	+=k*pr[mp]*sr;
      dyi[mp]	+=k*pi[mp]*sr;
      dyr[mm]	+=k*pr[mm]*sr;
      dyi[mm]	+=k*pi[mm]*sr;
    }
    mm	=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      tr	=mp*pr[mp];
      ti	=mp*pi[mp];
      Tr	=mm*pr[mm];
      Ti	=mm*pi[mm];
      k		=1;
      kp	=mp+k;
      km	=mm-k;
      while(kp < ESMp1){
	sr	=wg11c[k];
	si	=-wg11s[k];
	k++;
	dyr[kp]	+=kp*(sr*tr-si*ti);
	dyi[kp]	+=kp*(sr*ti+si*tr);
	kp++;
	dyr[km]	+=km*(sr*Tr+si*Ti);
	dyi[km]	+=km*(sr*Ti-si*Tr);
	km--;
      }
      k		=1;
      kp	=k+mm;
      km	=mp-k;
      while(kp < ESMp1){
	sr	=wg11c[k];
	si	=-wg11s[k];
	k++;
	dyr[kp]	+=kp*(sr*Tr-si*Ti);
	dyi[kp]	+=kp*(sr*Ti+si*Tr);
	kp++;
	dyr[km]	+=km*(sr*tr+si*ti);
	dyi[km]	+=km*(sr*ti-si*tr);
	km--;
      }
    }

    si	=wg12c[0];
    mm	=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      dyr[mp]	-=mp*dpi[mp]*si;
      dyi[mp]	+=mp*dpr[mp]*si;
      dyr[mm]	+=mp*dpi[mm]*si;
      dyi[mm]	-=mp*dpr[mm]*si;
    }
    tr	=dpr[0];
    ti	=dpi[0];
    km	=0;
    for(kp=1; kp < ESMp1; kp++){
      km--;
      sr	=kp*wg12s[kp];
      si	=kp*wg12c[kp];
      dyr[kp]	+=sr*tr-si*ti;
      dyi[kp]	+=sr*ti+si*tr;
      dyr[km]	+=sr*tr+si*ti;
      dyi[km]	+=sr*ti-si*tr;
    }
    mm	=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      tr	=dpr[mp];
      ti	=dpi[mp];
      Tr	=dpr[mm];
      Ti	=dpi[mm];
      k		=1;
      kp	=mp+k;
      km	=-kp;
      while(kp < ESMp1){
	sr	=kp*wg12s[k];
	si	=kp*wg12c[k];
	k++;
	dyr[kp]	+=sr*tr-si*ti;
	dyi[kp]	+=sr*ti+si*tr;
	kp++;
	dyr[km]	+=sr*Tr+si*Ti;
	dyi[km]	+=sr*Ti-si*Tr;
	km--;
      }
      k		=1;
      kp	=k+mm;
      km	=-kp;
      while(kp < ESMp1){
	sr	=kp*wg12s[k];
	si	=kp*wg12c[k];
	k++;
	dyr[kp]	+=sr*Tr-si*Ti;
	dyi[kp]	+=sr*Ti+si*Tr;
	kp++;
	dyr[km]	+=sr*tr+si*ti;
	dyi[km]	+=sr*ti-si*tr;
	km--;
      }
    }
    tr	=-yr[0];
    ti	=-yi[0];
    Tr	=-dyr[0];
    Ti	=-dyi[0];
    mm	=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      sr	=wduc[mp]*tr+wdus[mp]*ti+wuc[mp]*Tr+wus[mp]*Ti;
      si	=wduc[mp]*ti-wdus[mp]*tr+wuc[mp]*Ti-wus[mp]*Tr;
      dyr[mp]	+=sr;
      dyi[mp]	+=si;
      sr	=wduc[mp]*tr-wdus[mp]*ti+wuc[mp]*Tr-wus[mp]*Ti;
      si	=wduc[mp]*ti+wdus[mp]*tr+wuc[mp]*Ti+wus[mp]*Tr;
      dyr[mm]	+=sr;
      dyi[mm]	+=si;
    }
    yr[0]	=0.;
    yi[0]	=0.;
    dyr[0]	=0.;
    dyi[0]	=0.;
    /*}*/
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    dyr	+=ES2Mp1;
    dyi	+=ES2Mp1;
  }
  dyr	-=M0Mm;
  dyi	-=M0Mm;

  k		=*neq-1;
  tr		=dj_p*gyr[k]-wj_p;
  mm	=0;
  for(mp=1; mp < ESMp1; mp++){
    mm--;
    sr		=(wvc[mp]-wuc[mp]*wvc[0])*tr;
    si		=-(wvs[mp]-wus[mp]*wvc[0])*tr;
    dyr[mp]	+=sr;
    dyi[mp]	+=si;
    dyr[mm]	+=sr;
    dyi[mm]	-=si;
  }
}

void eq2DVgm00(int *neq,double *x,double *gyr,double *dgyr)
{
  double *pr,*pi,*yr,*yi,*dpr,*dpi,*dyr,*dyi;
  double *xr,*xi,*pU;
  double *mUr,*mUi;
  int i,j,k,m;
  double sr,si,tr,ti,Tr,Ti,wj_s,dj_s,wj_p,dj_p,rX0;
  double Y0r,Y0i;
  double gm,dgm,bL,dL,bF,dF,bV,bK,dK;

  int M1;
  double Rc[65],Rs[65],Rac[65],Ras[65],Raac[65],Raas[65],b,ba,baa,cs,sn;
  double r[128],ra[128],raa,rt[128],rat[128]
    ,za[128],zaa,zt[128],zat[128],D[128],Da[128],G[128];

  /* Getting metric coefficients{*/
  ESSetSplA(*x);
  splRA2(&b,&ba,&baa,ESsb,ESsb2a);
  splRA2(Rc,Rac,Raac,rcT,rcT2a);
  splRA2(Rs,Ras,Raas,rsT,rsT2a);

  m	=0;
  M1	=2;
  for(k=1; k < ESFp1; k++){
    m	+=ESNa1;
    splRA2(Rc+k,Rac+k,Raac+k,rcT+m,rcT2a+m);
    splRA2(Rs+k,Ras+k,Raas+k,rsT+m,rsT2a+m);
    Rc[k]	*=2.;
    Rs[k]	*=2.;
    Rac[k]	*=2.;
    Ras[k]	*=2.;
    Raac[k]	*=2.;
    Raas[k]	*=2.;
    if(Rac[k] != 0. || Ras[k] != 0.) M1=k+1;
  }
  for(j=0; j < ESNp; j++){
    r[j]	=ESaR0[0]+Rc[0];
    ra[j]	=Rac[0];
    raa		=Raac[0];
    rt[j]	=0.;
    rat[j]	=0.;
    cs		=EScs1[j];
    sn		=ESsn1[j];
    za[j]	=Ras[0]+ba*sn;
    zaa		=Raas[0]+baa*sn;
    zt[j]	=b*cs;
    zat[j]	=ba*cs;
    m	=0;
    for(k=1; k < M1; k++){
      m		+=j;
      if(m >= ESNp) m -=ESNp;
      cs	=EScs1[m];
      sn	=ESsn1[m];
      r[j]	+=Rc[k]*cs+Rs[k]*sn;
      ra[j]	+=Rac[k]*cs+Ras[k]*sn;
      raa	+=Raac[k]*cs+Raas[k]*sn;
      rt[j]	+=k*(Rs[k]*cs-Rc[k]*sn);
      rat[j]	+=k*(Ras[k]*cs-Rac[k]*sn);
    }
    D[j]	=ra[j]*zt[j]-rt[j]*za[j];
    G[j]	=(rt[j]*rt[j]+zt[j]*zt[j])/(r[j]*D[j]);
    Da[j]	=raa*zt[j]+ra[j]*zat[j]-rat[j]*za[j]-rt[j]*zaa;
  }
  ESP2F(wg22c,wg22s,G,ES2Mp);
  dK	=0.;
  for(j=0; j < ESNp; j++){
    sr	=1./(r[j]*D[j]);
    dK	+=(2.*rt[j]*rat[j]+2.*zt[j]*zat[j]
	   -G[j]*(Da[j]*r[j]+D[j]*ra[j]))*sr;
    G[j]	=(ra[j]*rt[j]+za[j]*zt[j])*sr;
  }
  dK	/=ESNp;
  ESP2F(wg12c,wg12s,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=(ra[j]*ra[j]+za[j]*za[j])/(r[j]*D[j]);
  }
  ESP2F(wg11c,wg11s,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=D[j]*ESaR0[0]/r[j];
  }
  ESP2F(wuc,wus,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=(Da[j]*ESaR0[0]-G[j]*ra[j])/r[j];
  }
  ESP2F(wduc,wdus,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=D[j]*(r[j]/ESaR0[0]-ESaR0[0]/r[j]);
  }
  ESP2F(wvc,wvs,G,ES2Mp);
  for(k=1; k < ES2Mp1; k++){
    i	=-k;
    wg11c[i]	=wg11c[k];
    wg12c[i]	=wg12c[k];
    wg22c[i]	=wg22c[k];
    wg11s[i]	=wg11s[k];
    wg12s[i]	=wg12s[k];
    wg22s[i]	=wg22s[k];
    
    wuc[i]	=wuc[k];
    wus[i]	=wus[k];
    wduc[i]	=wduc[k];
    wdus[i]	=wdus[k];
    wvc[i]	=wvc[k];
    wvs[i]	=wvs[k];

    wg11s[k]	=-wg11s[i];
    wg12s[k]	=-wg12s[i];
    wg22s[k]	=-wg22s[i];
    wvs[k]	=-wvs[i];
    wus[k]	=-wus[i];
    wdus[k]	=-wdus[i];
  }
  bK		=wg22c[0];
  bL		=wuc[0];
  dL		=wduc[0];
  bV		=wvc[0];
  /*}*/

  /* Getting 1D profiles{*/
  ESSetSplDPr(*x);
  ESSetSplDCr(*x);
  k	=*neq-1;
  splRA(&bF,&dF,ESaF,ESaF2a);
  splRDCr(&gm,&dgm,7);
  dgyr[k]	=-bL*gm*bF*rR0;
  Ti		=gm != 0. ? 1./dgyr[k] : 0.;
  switch(ESEqSolvInPr){
  case 0:
    splRDPr(&wj_p,&dj_p,0);
    break;
  case 1:
    splRDPr(&wj_p,&dj_p,1);
    wj_p	*=R0;
    dj_p	*=R0;
    break;
  case 2:
    splRDPr2(&si,&wj_p,&dj_p,2);

    break;
  }
  if(ESEqSolvInPr == 2){
    Tr		=R0*Ti;
    wj_p	*=Tr;
    wj_s	=-dF/(gm*wuc[0])+wj_p;
    dj_p	=dj_p*Tr+wj_p*(dK+(wj_s*wuc[0]+wj_p*wvc[0])*Ti)/wg22c[0];
  }
  else{
    sr	=bK*bL*gm*gm*rR0;
    si	=1./(1.+sr);
    tr	=((dK*bL+bK*dL)*gm+bK*bL*dgm)*rR0*si;
    ti	=wj_p*(bL+bV)*si;
    wj_s=wj_p-(ti-tr*bF)/bL;
  }
  splRA(&tr,&dj_s,ESjs,ESjs2a);
  dj_p	*=Ti;
  dj_s	*=Ti;
  k--;
  dgyr[k]	=-(wj_s*wuc[0]+wj_p*wvc[0]);

  mUr	=Ur+ES2Mp;
  mUi	=Ui+ES2Mp;
  mUr[0]	=-dj_s*wuc[0]-dj_p*wvc[0];
  mUi[0]	=0.;
  sr		=1./wuc[0];
  si		=-wduc[0];
  wuc[0]	=1.;
  wus[0]	=0.;
  wduc[0]	=0.;
  wdus[0]	=0.;
  for(k=1; k < ES2Mp1; k++){
    mUr[k]	=-dj_s*wuc[k]-dj_p*wvc[k];
    mUi[k]	=-dj_s*wus[k]-dj_p*wvs[k];
    wuc[k]	*=sr;
    wus[k]	*=sr;
    wduc[k]	=(wduc[k]+wuc[k]*si)*sr;
    wdus[k]	=(wdus[k]+wus[k]*si)*sr;
    i	=-k;
    mUr[i]	= mUr[k];
    wuc[i]	= wuc[k];
    wduc[i]	= wduc[k];
    mUi[i]	=-mUi[k];
    wus[i]	=-wus[k];
    wdus[i]	=-wdus[k];
  }
  /*}*/
  /* Cholesky decomposition of $\maK${*/
  for(i=0; i < ES2Mp1; i++){
    k	=ES2Mp1*i+i;
    for(j=i; j < ES2Mp1; j++){
      Kr[k]	=wg22c[i-j];
      Ki[k]	=wg22s[i-j];
      k++;
    }
  }
  if((k=CholDeComp(Kr,Ki,ES2Mp1))){
#ifdef XWIN
    extern char *CbUserWarning;
    extern char ESMessage[];
    sprintf(ESMessage,"CholDeComp() failed\nat eq2DVgm()\na=%5.3f i=%d"
	    ,*x,k-1);
    CbUserWarning	=ESMessage;
#endif
#ifdef DEBUG
    printf("??? x=%10.3e\n",*x);
    for(j=0; j < ESNp1; j++){
      sr	=wg22c[0];
      k	=0;
      for(m=0; m < ES2Mp1; m++){
	k	+=j;
	if(k >= ESNp) k	-=ESNp;
	sr	+=wg22c[m]*EScs1[k]+wg22s[m]*ESsn1[k];
      }
      printf("?--- g22[%2d]=%10.3e\n",j,sr);
    }
    printf("??? x=%10.3e\n",*x);
    for(m=0; m < ES2Mp1; m++){
      printf("?+?? g22[%2d]=%10.3e %10.3e\n",m,wg22c[m],wg22s[m]);
    }
#endif
    EcCholFl=1;
    return;
  }
  /*}*/
  
  /* Calculation of normalized solution with $Y_0=1$;*/
  xr	=Xr+ESMp;
  xi	=Xi+ESMp;
  for(k=-ESMp; k < ESMp1; k++){
    xr[k]	=wuc[k];
    xi[k]	=wus[k];
  }
  CholBksbComp(Kr,Ki,ES2Mp1,xr-ESMp,xi-ESMp);
  Tr	= xr[0];
  Ti	=-xi[0];
  sr	=1./(Tr*Tr+Ti*Ti);
  Tr	*=sr;
  Ti	*=sr;
  /*}*/
  
  pr	=gyr+ESMp;
  pi	=pr+M0Mm;
  yr	=pi+M0Mm;
  yi	=yr+M0Mm;
  dpr	=dgyr+ESMp;
  dpi	=dpr+M0Mm;
  dyr	=dpi+M0Mm;
  dyi	=dyr+M0Mm;
  for(m=0; m < ESMp1; m++){
    /* Calculation of $dgy${*/
    for(i=-ESMp; i < ESMp1; i++){
      dpr[i]	=yr[i];
      dpi[i]	=yi[i];
      for(j=-ESMp; j < ESMp1; j++){
	sr	=-j*pi[j];
	si	= j*pr[j];
	k	=i-j;
	tr	=wg12c[k];
	ti	=wg12s[k];
	dpr[i]	+=sr*tr-si*ti;
	dpi[i]	+=si*tr+sr*ti;
      }
      dyr[i]	=0.;
      dyi[i]	=0.;
    }
    CholBksbComp(Kr,Ki,ES2Mp1,dpr-ESMp,dpi-ESMp);
    /* Contribution from $X^0$;*/
    Y0r	=-dpr[0]*Tr+dpi[0]*Ti;
    Y0i	=-dpr[0]*Ti-dpi[0]*Tr;
    if(m == 0){
      k	=*neq-1;
      Y0r	+=dgyr[k]*Tr;
      Y0i	+=dgyr[k]*Ti;
    }
    for(i=-ESMp; i < ESMp1; i++){
      dpr[i]	+=Y0r*xr[i]-Y0i*xi[i];
      dpi[i]	+=Y0i*xr[i]+Y0r*xi[i];
    }
    /*}*/
    /* Calculation of $daY${*/
    for(j=-ESMp; j < ESMp1; j++){
      sr	=dpr[j];
      si	=dpi[j];
      for(i=-ESMp; i < ESMp1; i++){
	if(i){
	  k	=i-j;
	  tr	=-i*wg12s[k];
	  ti	= i*wg12c[k];
	  dyr[i]	+=sr*tr-si*ti;
	  dyi[i]	+=si*tr+sr*ti;
	}
      }
      if(j){
	sr	=pr[j];
	si	=pi[j];
	for(i=-ESMp; i < ESMp1; i++){
	  k	=i-j;
	  tr	=mUr[k]+i*wg11c[k]*j;
	  ti	=mUi[k]+i*wg11s[k]*j;
	  dyr[i]	+=sr*tr-si*ti;
	  dyi[i]	+=si*tr+sr*ti;
	}
      }
    }
    sr	=dyr[0];
    si	=dyi[0];
    for(i=-ESMp; i < ESMp1; i++){
      dyr[i]	-=wuc[i]*sr-wus[i]*si+wduc[i]*Y0r-wdus[i]*Y0i;
      dyi[i]	-=wuc[i]*si+wus[i]*sr+wduc[i]*Y0i+wdus[i]*Y0r;
    }
    yr[0]	=0.;
    yi[0]	=0.;
    dyr[0]	=0.;
    dyi[0]	=0.;
    /*}*/
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    dyr	+=ES2Mp1;
    dyi	+=ES2Mp1;
  }
  dyr	-=M0Mm;
  dyi	-=M0Mm;
  tr		=dj_p*gyr[(*neq-1)]-wj_p;
  for(k=-ESMp; k < ESMp1; k++){
    dyr[k]	+=(wvc[k]-wuc[k]*wvc[0])*tr;
    dyi[k]	+=(wvs[k]-wus[k]*wvc[0])*tr;
  }
}

void eq2DVgY(int *neq,double *x,double *gyr,double *dgyr)
{
  double *pr,*pi,*yr,*yi,*dpr,*dpi,*dyr,*dyi;
  double *xr,*xi,*pU;
  double *mUr,*mUi;
  int i,j,k,m;
  double sr,si,tr,ti,Tr,Ti,wj_s,dj_s,wj_p,dj_p,rX0;
  double Y0r,Y0i;
  double bL,dL,bF,dF,bV,bK,dK;
  double dgY,d2gY;

  int M1;
  double Rc[65],Rs[65],Rac[65],Ras[65],Raac[65],Raas[65],b,ba,baa,cs,sn;
  double r[128],ra[128],raa[128],rt[128],rat[128]
    ,za[128],zaa,zt[128],zat[128],D[128],Da[128],G[128];

  /* Getting metric coefficients{*/
  ESSetSplA(*x);
  splRA2(&b,&ba,&baa,ESsb,ESsb2a);
  splRA2(Rc,Rac,Raac,rcT,rcT2a);
  splRA2(Rs,Ras,Raas,rsT,rsT2a);
  m	=0;
  M1	=2;
  for(k=1; k < ESFp1; k++){
    m	+=ESNa1;
    splRA2(Rc+k,Rac+k,Raac+k,rcT+m,rcT2a+m);
    splRA2(Rs+k,Ras+k,Raas+k,rsT+m,rsT2a+m);
    Rc[k]	*=2.;
    Rs[k]	*=2.;
    Rac[k]	*=2.;
    Ras[k]	*=2.;
    Raac[k]	*=2.;
    Raas[k]	*=2.;
    if(Rc[k] != 0. || Rs[k] != 0.) M1=k+1;
  }
  for(j=0; j < ESNp; j++){
    r[j]	=ESaR0[0]+Rc[0];
    ra[j]	=Rac[0];
    raa[j]	=Raac[0];
    rt[j]	=0.;
    rat[j]	=0.;
    za[j]	=Ras[0]+ba*ESsn1[j];
    zaa		=Raas[0]+baa*ESsn1[j];
    zt[j]	=b*EScs1[j];
    zat[j]	=ba*EScs1[j];
    m	=0;
    for(k=1; k < ESFp1; k++){
      m	+=j;
      if(m >= ESNp) m -=ESNp;
      cs	=EScs1[m];
      sn	=ESsn1[m];
      r[j]	+=Rc[k]*cs+Rs[k]*sn;
      ra[j]	+=Rac[k]*cs+Ras[k]*sn;
      raa[j]	+=Raac[k]*cs+Raas[k]*sn;
      rt[j]	+=k*(Rs[k]*cs-Rc[k]*sn);
      rat[j]	+=k*(Ras[k]*cs-Rac[k]*sn);
    }
    D[j]	=ra[j]*zt[j]-rt[j]*za[j];
    G[j]	=(rt[j]*rt[j]+zt[j]*zt[j])/(r[j]*D[j]);
    Da[j]	=raa[j]*zt[j]+ra[j]*zat[j]-rat[j]*za[j]-rt[j]*zaa;
  }
  ESP2F(wg22c,wg22s,G,ES2Mp);
  dK	=0.;
  for(j=0; j < ESNp; j++){
    sr	=1./(r[j]*D[j]);
    dK	+=(2.*rt[j]*rat[j]+2.*zt[j]*zat[j]
	   -G[j]*(Da[j]*r[j]+D[j]*ra[j]))*sr;
    G[j]	=(ra[j]*rt[j]+za[j]*zt[j])*sr;
  }
  dK	/=ESNp;
  ESP2F(wg12c,wg12s,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=(ra[j]*ra[j]+za[j]*za[j])/(r[j]*D[j]);
  }
  ESP2F(wg11c,wg11s,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=D[j]*ESaR0[0]/r[j];
  }
  ESP2F(wuc,wus,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=(Da[j]*ESaR0[0]-G[j]*ra[j])/r[j];
  }
  ESP2F(wduc,wdus,G,ES2Mp);
  for(j=0; j < ESNp; j++){
    G[j]	=D[j]*(r[j]/ESaR0[0]-ESaR0[0]/r[j]);
  }
  ESP2F(wvc,wvs,G,ES2Mp);
  for(k=1; k < ES2Mp1; k++){
    i	=-k;
    wg11c[i]	=wg11c[k];
    wg12c[i]	=wg12c[k];
    wg22c[i]	=wg22c[k];
    wg11s[i]	=wg11s[k];
    wg12s[i]	=wg12s[k];
    wg22s[i]	=wg22s[k];
    
    wuc[i]	=wuc[k];
    wus[i]	=wus[k];
    wduc[i]	=wduc[k];
    wdus[i]	=wdus[k];
    wvc[i]	=wvc[k];
    wvs[i]	=wvs[k];

    wg11s[k]	=-wg11s[i];
    wg12s[k]	=-wg12s[i];
    wg22s[k]	=-wg22s[i];
    wvs[k]	=-wvs[i];
    wus[k]	=-wus[i];
    wdus[k]	=-wdus[i];
  }
  bK		=wg22c[0];
  bL		=wuc[0];
  dL		=wduc[0];
  bV		=wvc[0];
  /*}*/

  /* Getting 1D profiles{*/
  ESSetSplDPr(*x);
  ESSetSplDCr(*x);
  k	=*neq-1;
  splRDCr2(&sr,&dgY,&d2gY,7);
  dgyr[k]	=dgY;
  Ti		=dgY != 0. ? 1./dgY : 0.;
  switch(ESEqSolvInPr){
  case 0:
    splRDPr(&wj_p,&dj_p,0);
    break;
  case 1:
    splRDPr(&wj_p,&dj_p,1);
    wj_p	*=R0;
    dj_p	*=R0;
    break;
  case 2:
    splRDPr2(&si,&wj_p,&dj_p,2);
    wj_p	*=R0/dgY;
    dj_p	=(dj_p*R0-wj_p*d2gY)/dgY;
    break;
  }
  k--;
  dgyr[k]	=dK*dgY+bK*d2gY;
  wj_s	=-(dgyr[k]+wj_p*bV)/bL;
  splRA(&tr,&dj_s,ESjs,ESjs2a);
  dj_p	*=Ti;
  dj_s	*=Ti;

  mUr	=Ur+ES2Mp;
  mUi	=Ui+ES2Mp;
  mUr[0]	=-dj_s*wuc[0]-dj_p*wvc[0];
  mUi[0]	=0.;
  sr		=1./wuc[0];
  si		=-wduc[0];
  wuc[0]	=1.;
  wus[0]	=0.;
  wduc[0]	=0.;
  wdus[0]	=0.;
  for(k=1; k < ES2Mp1; k++){
    mUr[k]	=-dj_s*wuc[k]-dj_p*wvc[k];
    mUi[k]	=-dj_s*wus[k]-dj_p*wvs[k];
    wuc[k]	*=sr;
    wus[k]	*=sr;
    wduc[k]	=(wduc[k]+wuc[k]*si)*sr;
    wdus[k]	=(wdus[k]+wus[k]*si)*sr;
    i	=-k;
    mUr[i]	= mUr[k];
    wuc[i]	= wuc[k];
    wduc[i]	= wduc[k];
    mUi[i]	=-mUi[k];
    wus[i]	=-wus[k];
    wdus[i]	=-wdus[k];
  }
  /*}*/
  /* Cholesky decomposition of $\maK${*/
  for(i=0; i < ES2Mp1; i++){
    k	=ES2Mp1*i+i;
    for(j=i; j < ES2Mp1; j++){
      Kr[k]	=wg22c[i-j];
      Ki[k]	=wg22s[i-j];
      k++;
    }
  }
  if((k=CholDeComp(Kr,Ki,ES2Mp1))){
#ifdef XWIN
    extern char *CbUserWarning;
    extern char ESMessage[];
    sprintf(ESMessage,"CholDeComp() failed\nat eq2DVgY()\na=%5.3f i=%d"
	    ,*x,k-1);
    CbUserWarning	=ESMessage;
#endif
#ifdef DEBUG
    printf("??? x=%10.3e\n",*x);
    for(j=0; j < ESNp1; j++){
      sr	=wg22c[0];
      k	=0;
      for(m=0; m < ES2Mp1; m++){
	k	+=j;
	if(k >= ESNp) k	-=ESNp;
	sr	+=wg22c[m]*EScs1[k]+wg22s[m]*ESsn1[k];
      }
      printf("?--- g22[%2d]=%10.3e\n",j,sr);
    }
    printf("??? x=%10.3e\n",*x);
    for(m=0; m < ES2Mp1; m++){
      printf("?+?? g22[%2d]=%10.3e %10.3e\n",m,wg22c[m],wg22s[m]);
    }
#endif
    EcCholFl=1;
    return;
  }
  /*}*/
  
  /* Calculation of normalized solution with $Y_0=1$;*/
  xr	=Xr+ESMp;
  xi	=Xi+ESMp;
  for(k=-ESMp; k < ESMp1; k++){
    xr[k]	=wuc[k];
    xi[k]	=wus[k];
  }
  CholBksbComp(Kr,Ki,ES2Mp1,xr-ESMp,xi-ESMp);
  Tr	= xr[0];
  Ti	=-xi[0];
  sr	=1./(Tr*Tr+Ti*Ti);
  Tr	*=sr;
  Ti	*=sr;
  /*}*/
  
  pr	=gyr+ESMp;
  pi	=pr+M0Mm;
  yr	=pi+M0Mm;
  yi	=yr+M0Mm;
  dpr	=dgyr+ESMp;
  dpi	=dpr+M0Mm;
  dyr	=dpi+M0Mm;
  dyi	=dyr+M0Mm;
  for(m=0; m < ESMp1; m++){
    /* Calculation of $dgy${*/
    for(i=-ESMp; i < ESMp1; i++){
      dpr[i]	=yr[i];
      dpi[i]	=yi[i];
      for(j=-ESMp; j < ESMp1; j++){
	sr	=-j*pi[j];
	si	= j*pr[j];
	k	=i-j;
	tr	=wg12c[k];
	ti	=wg12s[k];
	dpr[i]	+=sr*tr-si*ti;
	dpi[i]	+=si*tr+sr*ti;
      }
      dyr[i]	=0.;
      dyi[i]	=0.;
    }
    CholBksbComp(Kr,Ki,ES2Mp1,dpr-ESMp,dpi-ESMp);
    /* Contribution from $X^0$;*/
    Y0r	=-dpr[0]*Tr+dpi[0]*Ti;
    Y0i	=-dpr[0]*Ti-dpi[0]*Tr;
    if(m == 0){
      k	=*neq-1;
      Y0r	+=dgyr[k]*Tr;
      Y0i	+=dgyr[k]*Ti;
    }
    for(i=-ESMp; i < ESMp1; i++){
      dpr[i]	+=Y0r*xr[i]-Y0i*xi[i];
      dpi[i]	+=Y0i*xr[i]+Y0r*xi[i];
    }
    /*}*/
    /* Calculation of $daY${*/
    for(j=-ESMp; j < ESMp1; j++){
      sr	=dpr[j];
      si	=dpi[j];
      for(i=-ESMp; i < ESMp1; i++){
	if(i){
	  k	=i-j;
	  tr	=-i*wg12s[k];
	  ti	= i*wg12c[k];
	  dyr[i]	+=sr*tr-si*ti;
	  dyi[i]	+=si*tr+sr*ti;
	}
      }
      if(j){
	sr	=pr[j];
	si	=pi[j];
	for(i=-ESMp; i < ESMp1; i++){
	  k	=i-j;
	  tr	=mUr[k]+i*wg11c[k]*j;
	  ti	=mUi[k]+i*wg11s[k]*j;
	  dyr[i]	+=sr*tr-si*ti;
	  dyi[i]	+=si*tr+sr*ti;
	}
      }
    }
    sr	=dyr[0];
    si	=dyi[0];
    for(i=-ESMp; i < ESMp1; i++){
      dyr[i]	-=wuc[i]*sr-wus[i]*si+wduc[i]*Y0r-wdus[i]*Y0i;
      dyi[i]	-=wuc[i]*si+wus[i]*sr+wduc[i]*Y0i+wdus[i]*Y0r;
    }
    yr[0]	=0.;
    yi[0]	=0.;
    dyr[0]	=0.;
    dyi[0]	=0.;
    /*}*/
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
    dpr	+=ES2Mp1;
    dpi	+=ES2Mp1;
    dyr	+=ES2Mp1;
    dyi	+=ES2Mp1;
  }
  dyr	-=M0Mm;
  dyi	-=M0Mm;
  tr		=dj_p*gyr[(*neq-1)]-wj_p;
  for(k=-ESMp; k < ESMp1; k++){
    dyr[k]	+=(wvc[k]-wuc[k]*wvc[0])*tr;
    dyi[k]	+=(wvs[k]-wus[k]*wvc[0])*tr;
  }
}
#endif/*stg_2DAdvR*/

#ifndef stg_2DSolv
int EZresete(gyr,yh,nord,gyer,gyei,dgyer,dgyei,icur)
     double*gyr,*yh,*gyer,*gyei,*dgyer,*dgyei;
     int*nord,icur;
{
  int NpY,NYp,Nyh,kmr,kmi,kYr,kYi;
  int i,m,m1,mm,k,kr,ki,nh;
  double s,sr,si,gar,gai;
  double *yr,*yi,*Yr,*Yi;

  yr	=Ur	+ES2Mp1;
  yi	=Ui	+ES2Mp1;
  Yr	=dUr	+ES2Mp1;
  Yi	=dUi	+ES2Mp1;

  NpY	=2*M0Mm;
  NYp	=2*M0Mm+2;
  Nyh	=((*nord)+1)*(NpY+NYp);

  for(m=ESMp; m > 0; m--){
    kmr		=ES2Mp1*m;
    kmi		=kmr+M0Mm;
    kYr		=kmi+M0Mm;
    kYi		=kYr+M0Mm;
    for(k=0; k < ES2Mp1; k++){
      Ur[k]	=gyr[kmr];
      Ui[k]	=gyr[kmi];
      Xr[k]	=gyr[kYr];
      Xi[k]	=gyr[kYi];
      kmr++;
      kmi++;
      kYr++;
      kYi++;
    }
    k		=ESMp+m;
    sr		=Ur[k];
    si		=Ui[k];
    gar		=1./(sr*sr+si*si);
    sr		*=gar;
    si		*=gar;
    for(m1=0; m1 < m; m1++){
      kmr	=ES2Mp1*m1+ESMp-m;
      kmi	=kmr+M0Mm;
      gar	=gyr[kmi]*si-gyr[kmr]*sr;
      gai	=-gyr[kmi]*sr-gyr[kmr]*si;
      gdjsr1[m1]	+=gar*gdjsr1[m]+gai*gdjsi1[m];
      gdjsi1[m1]	+=gai*gdjsr1[m]-gar*gdjsi1[m];
      kmr	=ES2Mp1*m1;
      kmi	=kmr+M0Mm;
      kYr	=kmi+M0Mm;
      kYi	=kYr+M0Mm;
      k		=ES2Mp1;
      kr	=ES2Mp1*(m+1);
      ki	=kr+M0Mm;
      for(mm=0; mm < ES2Mp1; mm++){
	k--;
	kr--;
	ki--;
	gyr[kmr]	+=gar*Ur[k]+gai*Ui[k];
	gyr[kmi]	+=gai*Ur[k]-gar*Ui[k];
	gyr[kYr]	+=gar*Xr[k]+gai*Xi[k];
	gyr[kYi]	+=gai*Xr[k]-gar*Xi[k];
	for(nh=0; nh < Nyh; nh+=NYp){
	  yh[nh+kmr]	+=gar*yh[nh+kr]+gai*yh[nh+ki];
	  yh[nh+kmi]	+=gai*yh[nh+kr]-gar*yh[nh+ki];
	  nh		+=NpY;
	  yh[nh+kmr]	+=gar*yh[nh+kr]+gai*yh[nh+ki];
	  yh[nh+kmi]	+=gai*yh[nh+kr]-gar*yh[nh+ki];
	}
	kmr++;
	kmi++;
	kYr++;
	kYi++;
      }
      for(i=0; i <= icur; i++){
	k	=M0Mm*i+ES2Mp1*(m+1);
	kmr	=M0Mm*i+ES2Mp1*m1;
	for(mm=0; mm < ES2Mp1; mm++){
	  k--;
	  gyer[kmr]	+=gar*gyer[k]+gai*gyei[k];
	  gyei[kmr]	+=gai*gyer[k]-gar*gyei[k];
	  dgyer[kmr]	+=gar*dgyer[k]+gai*dgyei[k];
	  dgyei[kmr]	+=gai*dgyer[k]-gar*dgyei[k];
	  aYer[kmr]	+=gar*aYer[k]+gai*aYei[k];
	  aYei[kmr]	+=gai*aYer[k]-gar*aYei[k];
	  daYer[kmr]	+=gar*daYer[k]+gai*daYei[k];
	  daYei[kmr]	+=gai*daYer[k]-gar*daYei[k];
	  kmr++;
	}
      }
    }
    for(m1=m+1; m1 < ESMp1; m1++){
      kmr	=ES2Mp1*m1+ESMp-m;
      kmi	=kmr+M0Mm;
      gar	=gyr[kmi]*si-gyr[kmr]*sr;
      gai	=-gyr[kmi]*sr-gyr[kmr]*si;
      gdjsr1[m1]+=gar*gdjsr1[m]+gai*gdjsi1[m];
      gdjsi1[m1]+=gai*gdjsr1[m]-gar*gdjsi1[m];
      kmr	=ES2Mp1*m1;
      kmi	=kmr+M0Mm;
      kYr	=kmi+M0Mm;
      kYi	=kYr+M0Mm;
      k		=ES2Mp1;
      kr	=ES2Mp1*(m+1);
      ki	=kr+M0Mm;
      for(mm=0; mm < ES2Mp1; mm++){
	k--;
	kr--;
	ki--;
	gyr[kmr]	+=gar*Ur[k]+gai*Ui[k];
	gyr[kmi]	+=gai*Ur[k]-gar*Ui[k];
	gyr[kYr]	+=gar*Xr[k]+gai*Xi[k];
	gyr[kYi]	+=gai*Xr[k]-gar*Xi[k];
	for(nh=0; nh < Nyh; nh+=NYp){
	  yh[nh+kmr]	+=gar*yh[nh+kr]+gai*yh[nh+ki];
	  yh[nh+kmi]	+=gai*yh[nh+kr]-gar*yh[nh+ki];
	  nh		+=NpY;
	  yh[nh+kmr]	+=gar*yh[nh+kr]+gai*yh[nh+ki];
	  yh[nh+kmi]	+=gai*yh[nh+kr]-gar*yh[nh+ki];
	}
	kmr++;
	kmi++;
	kYr++;
	kYi++;
      }
      for(i=0; i <= icur; i++){
	k	=M0Mm*i+ES2Mp1*(m+1);
	kmr	=M0Mm*i+ES2Mp1*m1;
	for(mm=0; mm < ES2Mp1; mm++){
	  k--;
	  gyer[kmr]	+=gar*gyer[k]+gai*gyei[k];
	  gyei[kmr]	+=gai*gyer[k]-gar*gyei[k];
	  dgyer[kmr]	+=gar*dgyer[k]+gai*dgyei[k];
	  dgyei[kmr]	+=gai*dgyer[k]-gar*dgyei[k];
	  aYer[kmr]	+=gar*aYer[k]+gai*aYei[k];
	  aYei[kmr]	+=gai*aYer[k]-gar*aYei[k];
	  daYer[kmr]	+=gar*daYer[k]+gai*daYei[k];
	  daYei[kmr]	+=gai*daYer[k]-gar*daYei[k];
	  kmr++;
	}
      }
    }
    m1	=m;
    kmr	=ES2Mp1*m1+ESMp-m;
    kmi	=kmr+M0Mm;
    gar	=gyr[kmi]*si-gyr[kmr]*sr;
    gai	=-gyr[kmi]*sr-gyr[kmr]*si;
    s		=gdjsr1[m];
    gdjsr1[m1]	+=gar*s+gai*gdjsi1[m];
    gdjsi1[m1]	+=gai*s-gar*gdjsi1[m];
    kmr	=ES2Mp1*m1;
    kmi	=kmr+M0Mm;
    kYr	=kmi+M0Mm;
    kYi	=kYr+M0Mm;
    k	=ES2Mp1;
    for(mm=0; mm < ES2Mp1; mm++){
      k--;
      gyr[kmr]	+=gar*Ur[k]+gai*Ui[k];
      gyr[kmi]	+=gai*Ur[k]-gar*Ui[k];
      gyr[kYr]	+=gar*Xr[k]+gai*Xi[k];
      gyr[kYi]	+=gai*Xr[k]-gar*Xi[k];
      kmr++;
      kmi++;
      kYr++;
      kYi++;
    }
    for(nh=0; nh < Nyh; nh+=NpY+NYp){
      kmr	=nh+ES2Mp1*m;
      kmi	=kmr+M0Mm;
      kYr	=kmi+M0Mm;
      kYi	=kYr+M0Mm;
      for(mm=0; mm < ES2Mp1; mm++){
	yr[mm]	=yh[kmr];
	yi[mm]	=yh[kmi];
	Yr[mm]	=yh[kYr];
	Yi[mm]	=yh[kYi];
	kmr++;
	kmi++;
	kYr++;
	kYi++;
      }
      kmr	=nh+ES2Mp1*m1;
      kmi	=kmr+M0Mm;
      kYr	=kmi+M0Mm;
      kYi	=kYr+M0Mm;
      k		=ES2Mp1;
      for(mm=0; mm < ES2Mp1; mm++){
	k--;
	yh[kmr]	+=gar*yr[k]+gai*yi[k];
	yh[kmi]	+=gai*yr[k]-gar*yi[k];
	yh[kYr]	+=gar*Yr[k]+gai*Yi[k];
	yh[kYi]	+=gai*Yr[k]-gar*Yi[k];
	kmr++;
	kmi++;
	kYr++;
	kYi++;
      }
    }
    for(i=0; i <= icur; i++){
      k		=M0Mm*i+ES2Mp1*m;
      for(mm=0; mm < ES2Mp1; mm++){
	yr[mm]	=gyer[k];
	yi[mm]	=gyei[k];
	Yr[mm]	=dgyer[k];
	Yi[mm]	=dgyei[k];
	k++;
      }
      kmr	=M0Mm*i+ES2Mp1*m1;
      k		=ES2Mp1;
      for(mm=0; mm < ES2Mp1; mm++){
	k--;
	gyer[kmr]	+=gar*yr[k]+gai*yi[k];
	gyei[kmr]	+=gai*yr[k]-gar*yi[k];
	dgyer[kmr]	+=gar*Yr[k]+gai*Yi[k];
	dgyei[kmr]	+=gai*Yr[k]-gar*Yi[k];
	kmr++;
      }
      k		=M0Mm*i+ES2Mp1*m;
      for(mm=0; mm < ES2Mp1; mm++){
	yr[mm]	=aYer[k];
	yi[mm]	=aYei[k];
	Yr[mm]	=daYer[k];
	Yi[mm]	=daYei[k];
	k++;
      }
      kmr	=M0Mm*i+ES2Mp1*m1;
      k		=ES2Mp1;
      for(mm=0; mm < ES2Mp1; mm++){
	k--;
	aYer[kmr]	+=gar*yr[k]+gai*yi[k];
	aYei[kmr]	+=gai*yr[k]-gar*yi[k];
	daYer[kmr]	+=gar*Yr[k]+gai*Yi[k];
	daYei[kmr]	+=gai*Yr[k]-gar*Yi[k];
	kmr++;
      }
    }
    kmr		=ES2Mp1*m;
    kmi		=kmr+M0Mm;
    kYr		=kmi+M0Mm;
    kYi		=kYr+M0Mm;
    for(k=0; k < ES2Mp1; k++){
      Ur[k]	=gyr[kmr];
      Ui[k]	=gyr[kmi];
      Xr[k]	=gyr[kYr];
      Xi[k]	=gyr[kYi];
      kmr++;
      kmi++;
      kYr++;
      kYi++;
    }
    k		=ESMp+m;
    sr		=Ur[k];
    si		=-Ui[k];
    gar		=1./(sr*sr+si*si);
    sr		*=gar;
    si		*=gar;
    for(m1=0; m1 < m; m1++){
      kmr	=ES2Mp1*m1+ESMp+m;
      kmi	=kmr+M0Mm;
      gar	=gyr[kmi]*si-gyr[kmr]*sr;
      gai	=-gyr[kmi]*sr-gyr[kmr]*si;
      gdjsr1[m1]+=gar*gdjsr1[m]-gai*gdjsi1[m];
      gdjsi1[m1]+=gai*gdjsr1[m]+gar*gdjsi1[m];
      kmr	=ES2Mp1*m1;
      kmi	=kmr+M0Mm;
      kYr	=kmi+M0Mm;
      kYi	=kYr+M0Mm;
      for(mm=0; mm < ES2Mp1; mm++){
	gyr[kmr]	+=gar*Ur[mm]-gai*Ui[mm];
	gyr[kmi]	+=gai*Ur[mm]+gar*Ui[mm];
	gyr[kYr]	+=gar*Xr[mm]-gai*Xi[mm];
	gyr[kYi]	+=gai*Xr[mm]+gar*Xi[mm];
	kmr++;
	kmi++;
	kYr++;
	kYi++;
      }
      for(nh=0; nh < Nyh; nh +=NpY+NYp){
	kmr	=nh+ES2Mp1*m;
	kmi	=kmr+M0Mm;
	kYr	=kmi+M0Mm;
	kYi	=kYr+M0Mm;
	for(mm=0; mm < ES2Mp1; mm++){
	  yr[mm]	=yh[kmr];
	  yi[mm]	=yh[kmi];
	  Yr[mm]	=yh[kYr];
	  Yi[mm]	=yh[kYi];
	  kmr++;
	  kmi++;
	  kYr++;
	  kYi++;
	}
	kmr	=nh+ES2Mp1*m1;
	kmi	=kmr+M0Mm;
	kYr	=kmi+M0Mm;
	kYi	=kYr+M0Mm;
	for(mm=0; mm < ES2Mp1; mm++){
	  yh[kmr]	+=gar*yr[mm]-gai*yi[mm];
	  yh[kmi]	+=gai*yr[mm]+gar*yi[mm];
	  yh[kYr]	+=gar*Yr[mm]-gai*Yi[mm];
	  yh[kYi]	+=gai*Yr[mm]+gar*Yi[mm];
	  kmr++;
	  kmi++;
	  kYr++;
	  kYi++;
	}
      }
      for(i=0; i <=icur; i++){
	kmr	=M0Mm*i+ES2Mp1*m1;
	kYr	=M0Mm*i+ES2Mp1*m;
	for(mm=0; mm < ES2Mp1; mm++){
	  gyer[kmr]	+=gar*gyer[kYr]-gai*gyei[kYr];
	  gyei[kmr]	+=gai*gyer[kYr]+gar*gyei[kYr];
	  dgyer[kmr]	+=gar*dgyer[kYr]-gai*dgyei[kYr];
	  dgyei[kmr]	+=gai*dgyer[kYr]+gar*dgyei[kYr];
	  aYer[kmr]	+=gar*aYer[kYr]-gai*aYei[kYr];
	  aYei[kmr]	+=gai*aYer[kYr]+gar*aYei[kYr];
	  daYer[kmr]	+=gar*daYer[kYr]-gai*daYei[kYr];
	  daYei[kmr]	+=gai*daYer[kYr]+gar*daYei[kYr];
	  kYr++;
	  kmr++;
	}
      }
    }
  }
  return(0);
}

int reseteRKOld(int icur)
{
  int i,j,m,m1,mm,k;
  double s,sr,si,gar,gai,tr,ti,Tr,Ti;
  double *dXr,*dXi;

  if(icur > ESNa) icur=ESNa;

  m	=ESNa1*M0Mm;
  dXr	=Xr	+ES2Mp1;
  dXi	=Xi	+ES2Mp1;
  for(m=ESMp; m > 0; m--){
    j	=ES2Mp1*m;
    for(k=0; k < ES2Mp1; k++){
      Ur[k]	=gpr[j];
      Ui[k]	=gpi[j];
      dUr[k]	=dgpr[j];
      dUi[k]	=dgpi[j];
      Xr[k]	=Yr[j];
      Xi[k]	=Yi[j];
      dXr[k]	=dYr[j];
      dXi[k]	=dYi[j];
      j++;
    }
    k		=ESMp+m;
    sr		=Ur[k];
    si		=Ui[k];
    gar		=1./(sr*sr+si*si);
    sr		*=gar;
    si		*=gar;
    for(m1=0; m1 < ESMp1; m1++){
      j	=ES2Mp1*m1+ESMp-m;
      gar	= gpi[j]*si-gpr[j]*sr;
      gai	=-gpi[j]*sr-gpr[j]*si;
      j		=ES2Mp1*m1;
      k		=ES2Mp1;
      while(k > 0){
	k--;
	gpr[j]	+=gar*Ur[k]+gai*Ui[k];
	gpi[j]	+=gai*Ur[k]-gar*Ui[k];
	dgpr[j]	+=gar*dUr[k]+gai*dUi[k];
	dgpi[j]	+=gai*dUr[k]-gar*dUi[k];
	Yr[j]	+=gar*Xr[k]+gai*Xi[k];
	Yi[j]	+=gai*Xr[k]-gar*Xi[k];
	dYr[j]	+=gar*dXr[k]+gai*dXi[k];
	dYi[j]	+=gai*dXr[k]-gar*dXi[k];
	j++;
      }
      if(m1 != m){
	for(i=0; i <= icur; i++){
	  k	=M0Mm*i+ES2Mp1*(m+1);
	  j	=M0Mm*i+ES2Mp1*m1;
	  for(mm=0; mm < ES2Mp1; mm++){
	    k--;
	    EZgper[j]	+=gar*EZgper[k]+gai*EZgpei[k];
	    EZgpei[j]	+=gai*EZgper[k]-gar*EZgpei[k];
	    EZdgper[j]	+=gar*EZdgper[k]+gai*EZdgpei[k];
	    EZdgpei[j]	+=gai*EZdgper[k]-gar*EZdgpei[k];
	    aYer[j]	+=gar*aYer[k]+gai*aYei[k];
	    aYei[j]	+=gai*aYer[k]-gar*aYei[k];
	    daYer[j]	+=gar*daYer[k]+gai*daYei[k];
	    daYei[j]	+=gai*daYer[k]-gar*daYei[k];
	    j++;
	  }
	}
      }
      else{
	for(i=0; i <= icur; i++){
	  j	=M0Mm*i+ES2Mp1*m+ESMp;
	  tr	=EZgper[j];
	  ti	=EZgpei[j];
	  EZgper[j]	+=gar*tr+gai*ti;
	  EZgpei[j]	+=gai*tr-gar*ti;
	  tr	=EZdgper[j];
	  ti	=EZdgpei[j];
	  EZdgper[j]	+=gar*tr+gai*ti;
	  EZdgpei[j]	+=gai*tr-gar*ti;
	  tr	=aYer[j];
	  ti	=aYei[j];
	  aYer[j]	+=gar*tr+gai*ti;
	  aYei[j]	+=gai*tr-gar*ti;
	  tr	=daYer[j];
	  ti	=daYei[j];
	  daYer[j]	+=gar*tr+gai*ti;
	  daYei[j]	+=gai*tr-gar*ti;
	  k	=j;
	  for(mm=1; mm < ESMp1; mm++){
	    j++;
	    k--;
	    tr	=EZgper[j];	
	    ti	=EZgpei[j];	
	    Tr	=EZgper[k];	
	    Ti	=EZgpei[k];	
	    EZgper[k]	+=gar*tr+gai*ti;
	    EZgpei[k]	+=gai*tr-gar*ti;
	    EZgper[j]	+=gar*Tr+gai*Ti;
	    EZgpei[j]	+=gai*Tr-gar*Ti;
	    tr	=EZdgper[j];	
	    ti	=EZdgpei[j];	
	    Tr	=EZdgper[k];	
	    Ti	=EZdgpei[k];	
	    EZdgper[k]	+=gar*tr+gai*ti;
	    EZdgpei[k]	+=gai*tr-gar*ti;
	    EZdgper[j]	+=gar*Tr+gai*Ti;
	    EZdgpei[j]	+=gai*Tr-gar*Ti;

	    tr	=aYer[j];	
	    ti	=aYei[j];	
	    Tr	=aYer[k];	
	    Ti	=aYei[k];	
	    aYer[k]	+=gar*tr+gai*ti;
	    aYei[k]	+=gai*tr-gar*ti;
	    aYer[j]	+=gar*Tr+gai*Ti;
	    aYei[j]	+=gai*Tr-gar*Ti;
	    tr	=daYer[j];	
	    ti	=daYei[j];	
	    Tr	=daYer[k];	
	    Ti	=daYei[k];	
	    daYer[k]	+=gar*tr+gai*ti;
	    daYei[k]	+=gai*tr-gar*ti;
	    daYer[j]	+=gar*Tr+gai*Ti;
	    daYei[j]	+=gai*Tr-gar*Ti;
	  }
	}
      }
    }
    j		=ES2Mp1*m;
    for(k=0; k < ES2Mp1; k++){
      Ur[k]	=gpr[j];
      Ui[k]	=gpi[j];
      dUr[k]	=dgpr[j];
      dUi[k]	=dgpi[j];
      Xr[k]	=Yr[j];
      Xi[k]	=Yi[j];
      dXr[k]	=dYr[j];
      dXi[k]	=dYi[j];
      j++;
    }
    k		=ESMp+m;
    sr		=Ur[k];
    si		=-Ui[k];
    gar		=1./(sr*sr+si*si);
    sr		*=gar;
    si		*=gar;
    for(m1=0; m1 < m; m1++){
      j	=ES2Mp1*m1+ESMp+m;
      gar	= gpi[j]*si-gpr[j]*sr;
      gai	=-gpi[j]*sr-gpr[j]*si;
      j	=ES2Mp1*m1;
      for(k=0; k < ES2Mp1; k++){
	gpr[j]	+=gar*Ur[k]-gai*Ui[k];
	gpi[j]	+=gai*Ur[k]+gar*Ui[k];
	dgpr[j]	+=gar*dUr[k]-gai*dUi[k];
	dgpi[j]	+=gai*dUr[k]+gar*dUi[k];
	Yr[j]	+=gar*Xr[k]-gai*Xi[k];
	Yi[j]	+=gai*Xr[k]+gar*Xi[k];
	dYr[j]	+=gar*dXr[k]-gai*dXi[k];
	dYi[j]	+=gai*dXr[k]+gar*dXi[k];
	j++;
      }
      for(i=0; i <=icur; i++){
	j	=M0Mm*i+ES2Mp1*m1;
	k	=M0Mm*i+ES2Mp1*m;
	for(mm=0; mm < ES2Mp1; mm++){
	  EZgper[j]	+=gar*EZgper[k]-gai*EZgpei[k];
	  EZgpei[j]	+=gai*EZgper[k]+gar*EZgpei[k];
	  EZdgper[j]	+=gar*EZdgper[k]-gai*EZdgpei[k];
	  EZdgpei[j]	+=gai*EZdgper[k]+gar*EZdgpei[k];
	  aYer[j]	+=gar*aYer[k]-gai*aYei[k];
	  aYei[j]	+=gai*aYer[k]+gar*aYei[k];
	  daYer[j]	+=gar*daYer[k]-gai*daYei[k];
	  daYei[j]	+=gai*daYer[k]+gar*daYei[k];
	  j++;
	  k++;
	}
      }
    }
  }
  return(0);
}

int reseteRK(int icur)
{
  int i,j,m,m1,mm,k;
  double s,sr,si,gar,gai,tr,ti,Tr,Ti;
  double *dXr,*dXi;

  dXr	=Xr	+ES2Mp1;
  dXi	=Xi	+ES2Mp1;
  for(m=ESMp; m > 0; m--){
    j	=ES2Mp1*m+ESMp;
    k	=j+m;
    sr	=gpr[k];
    si	=gpi[k];
    gar	=1./(sr*sr+si*si);
    sr	*=gar;
    si	*=gar;
    k	=j-m;
    gar	= gpi[k]*si-gpr[k]*sr;
    gai	=-gpi[k]*sr-gpr[k]*si;
    mm	=j+ESMp;
    j	-=ESMp;
    for(k=0; k < ES2Mp1; k++,j++,mm--){
      Ur[k]	=gpr[j]+gar*gpr[mm]+gai*gpi[mm];
      Ui[k]	=gpi[j]+gai*gpr[mm]-gar*gpi[mm];
      dUr[k]	=dgpr[j]+gar*dgpr[mm]+gai*dgpi[mm];
      dUi[k]	=dgpi[j]+gai*dgpr[mm]-gar*dgpi[mm];
      Xr[k]	=Yr[j]+gar*Yr[mm]+gai*Yi[mm];
      Xi[k]	=Yi[j]+gai*Yr[mm]-gar*Yi[mm];
      dXr[k]	=dYr[j]+gar*dYr[mm]+gai*dYi[mm];
      dXi[k]	=dYi[j]+gai*dYr[mm]-gar*dYi[mm];
    }
    j	=ES2Mp1*m;
    for(k=0; k < ES2Mp1; k++,j++){
      dgpr[j]	=dUr[k];
      dgpi[j]	=dUi[k];
       gpr[j]	= Ur[k];
       gpi[j]	= Ui[k];
        Yr[j]	= Xr[k];
        Yi[j]	= Xi[k];
       dYr[j]	=dXr[k];
       dYi[j]	=dXi[k];
    }
    for(i=0; i <= icur; i++){
      j		=M0Mm*i+ES2Mp1*m+ESMp;
      tr	=EZgper[j];
      ti	=EZgpei[j];
      EZgper[j]	+=gar*tr+gai*ti;
      EZgpei[j]	+=gai*tr-gar*ti;
      tr	=EZdgper[j];
      ti	=EZdgpei[j];
      EZdgper[j]	+=gar*tr+gai*ti;
      EZdgpei[j]	+=gai*tr-gar*ti;
      tr	=aYer[j];
      ti	=aYei[j];
      aYer[j]	+=gar*tr+gai*ti;
      aYei[j]	+=gai*tr-gar*ti;
      tr	=daYer[j];
      ti	=daYei[j];
      daYer[j]	+=gar*tr+gai*ti;
      daYei[j]	+=gai*tr-gar*ti;
      k	=j;
      for(mm=1; mm < ESMp1; mm++){
	j++;
	k--;
	tr	=EZgper[j];	
	ti	=EZgpei[j];	
	Tr	=EZgper[k];	
	Ti	=EZgpei[k];	
	EZgper[k]	+=gar*tr+gai*ti;
	EZgpei[k]	+=gai*tr-gar*ti;
	EZgper[j]	+=gar*Tr+gai*Ti;
	EZgpei[j]	+=gai*Tr-gar*Ti;
	tr	=EZdgper[j];	
	ti	=EZdgpei[j];	
	Tr	=EZdgper[k];	
	Ti	=EZdgpei[k];	
	EZdgper[k]	+=gar*tr+gai*ti;
	EZdgpei[k]	+=gai*tr-gar*ti;
	EZdgper[j]	+=gar*Tr+gai*Ti;
	EZdgpei[j]	+=gai*Tr-gar*Ti;
	tr	=aYer[j];	
	ti	=aYei[j];	
	Tr	=aYer[k];	
	Ti	=aYei[k];	
	aYer[k]	+=gar*tr+gai*ti;
	aYei[k]	+=gai*tr-gar*ti;
	aYer[j]	+=gar*Tr+gai*Ti;
	aYei[j]	+=gai*Tr-gar*Ti;
	tr	=daYer[j];	
	ti	=daYei[j];	
	Tr	=daYer[k];	
	Ti	=daYei[k];	
	daYer[k]+=gar*tr+gai*ti;
	daYei[k]+=gai*tr-gar*ti;
	daYer[j]+=gar*Tr+gai*Ti;
	daYei[j]+=gai*Tr-gar*Ti;
      }
    }
    k		=ESMp+m;
    sr		=Ur[k];
    si		=Ui[k];
    gar		=1./(sr*sr+si*si);
    sr		*=gar;
    si		*=gar;
    for(m1=0; m1 < m; m1++){
      j	=ES2Mp1*m1+ESMp-m;
      gar	= gpi[j]*si-gpr[j]*sr;
      gai	=-gpi[j]*sr-gpr[j]*si;
      j		=ES2Mp1*m1;
      k		=ES2Mp1;
      while(k > 0){
	k--;
	gpr[j]	+=gar*Ur[k]+gai*Ui[k];
	gpi[j]	+=gai*Ur[k]-gar*Ui[k];
	dgpr[j]	+=gar*dUr[k]+gai*dUi[k];
	dgpi[j]	+=gai*dUr[k]-gar*dUi[k];
	Yr[j]	+=gar*Xr[k]+gai*Xi[k];
	Yi[j]	+=gai*Xr[k]-gar*Xi[k];
	dYr[j]	+=gar*dXr[k]+gai*dXi[k];
	dYi[j]	+=gai*dXr[k]-gar*dXi[k];
	j++;
      }
      for(i=0; i <= icur; i++){
	k	=M0Mm*i+ES2Mp1*(m+1);
	j	=M0Mm*i+ES2Mp1*m1;
	for(mm=0; mm < ES2Mp1; mm++){
	  k--;
	  EZgper[j]	+=gar*EZgper[k]+gai*EZgpei[k];
	  EZgpei[j]	+=gai*EZgper[k]-gar*EZgpei[k];
	  EZdgper[j]	+=gar*EZdgper[k]+gai*EZdgpei[k];
	  EZdgpei[j]	+=gai*EZdgper[k]-gar*EZdgpei[k];
	  aYer[j]	+=gar*aYer[k]+gai*aYei[k];
	  aYei[j]	+=gai*aYer[k]-gar*aYei[k];
	  daYer[j]	+=gar*daYer[k]+gai*daYei[k];
	  daYei[j]	+=gai*daYer[k]-gar*daYei[k];
	  j++;
	}
      }
      j		=ES2Mp1*m1+ESMp+m;
      gar	=-gpi[j]*si-gpr[j]*sr;
      gai	=-gpi[j]*sr+gpr[j]*si;
      j		=ES2Mp1*m1;
      for(k=0; k < ES2Mp1; k++,j++){
	gpr[j]	+=gar*Ur[k]-gai*Ui[k];
	gpi[j]	+=gai*Ur[k]+gar*Ui[k];
	dgpr[j]	+=gar*dUr[k]-gai*dUi[k];
	dgpi[j]	+=gai*dUr[k]+gar*dUi[k];
	Yr[j]	+=gar*Xr[k]-gai*Xi[k];
	Yi[j]	+=gai*Xr[k]+gar*Xi[k];
	dYr[j]	+=gar*dXr[k]-gai*dXi[k];
	dYi[j]	+=gai*dXr[k]+gar*dXi[k];
      }
      for(i=0; i <= icur; i++){
	k	=M0Mm*i+ES2Mp1*m;
	j	=M0Mm*i+ES2Mp1*m1;
	for(mm=0; mm < ES2Mp1; mm++,k++,j++){
	  EZgper[j]	+=gar*EZgper[k]-gai*EZgpei[k];
	  EZgpei[j]	+=gai*EZgper[k]+gar*EZgpei[k];
	  EZdgper[j]	+=gar*EZdgper[k]-gai*EZdgpei[k];
	  EZdgpei[j]	+=gai*EZdgper[k]+gar*EZdgpei[k];
	  aYer[j]	+=gar*aYer[k]-gai*aYei[k];
	  aYei[j]	+=gai*aYer[k]+gar*aYei[k];
	  daYer[j]	+=gar*daYer[k]-gai*daYei[k];
	  daYei[j]	+=gai*daYer[k]+gar*daYei[k];
	}
      }
    }
  }
  return(0);
}

int eqresetOrig(double *gyer,double *gyei,double *dgyer,double *dgyei)
{
  int k,ka;
  int i,m,m1,mm,km,km1,kmm,kmm1;
  double sr,si,gar,gai;
  double *pr,*pi,*dpr,*dpi,*EZvr,*vi,*dvr,*dvi;

  ka		=M0Mm*ESNa;
  for(m=ESMp; m > 0; m--){
    k		=ka+ES2Mp1*m+ESMp+m;
    pr		=gyer+k;
    pi		=gyei+k;
    sr		=*pr;
    si		=*pi;
    gar		=1./(sr*sr+si*si);
    sr		*=gar;
    si		*=gar;
    k		=ka+ESMp-m;
    pr		=gyer+k;
    pi		=gyei+k;
    for(m1=0; m1 < ESMp1; m1++){
      Xr[m1]	= (*pi)*si-(*pr)*sr;
      Xi[m1]	=-(*pi)*sr-(*pr)*si;
      pr	+=ES2Mp1;
      pi	+=ES2Mp1;
    }
    pr		=gyer;
    pi		=gyei;
    dpr		=dgyer;
    dpi		=dgyei;
    k		=ES2Mp1*m;
    EZvr		=gyer+k;
    vi		=gyei+k;
    dvr		=dgyer+k;
    dvi		=dgyei+k;
    for(i=0; i < ESNa1; i++){
      for(k=0; k < ES2Mp1; k++){
	Ur[k]	=*EZvr++;
	Ui[k]	=*vi++;
	dUr[k]	=*dvr++;
	dUi[k]	=*dvi++;
      }
      for(m1=0; m1 < ESMp1; m1++){
	gar	=Xr[m1];
	gai	=Xi[m1];
	mm	=ES2Mp1;
	while(mm > 0){
	  mm--;
	  *pr++		+=gar*Ur[mm]+gai*Ui[mm];
	  *pi++		+=gai*Ur[mm]-gar*Ui[mm];
	  *dpr++	+=gar*dUr[mm]+gai*dUi[mm];
	  *dpi++	+=gai*dUr[mm]-gar*dUi[mm];
	}
      }
      k		=M0Mm-ES2Mp1;
      EZvr	+=k;
      vi	+=k;
      dvr	+=k;
      dvi	+=k;
    }

    k		=ka+ES2Mp1*m+ESMp+m;
    pr		=gyer+k;
    pi		=gyei+k;
    sr		=*pr;
    si		=-(*pi);
    gar		=1./(sr*sr+si*si);
    sr		*=gar;
    si		*=gar;
    k		=ka+ESMp+m;
    pr		=gyer+k;
    pi		=gyei+k;
    for(m1=0; m1 < m; m1++){
      Xr[m1]	= (*pi)*si-(*pr)*sr;
      Xi[m1]	=-(*pi)*sr-(*pr)*si;
      pr	+=ES2Mp1;
      pi	+=ES2Mp1;
    }
    k		=ES2Mp1*m;
    EZvr		=gyer+k;
    vi		=gyei+k;
    dvr		=dgyer+k;
    dvi		=dgyei+k;
    for(i=0; i < ESNa1; i++){
      k		=M0Mm*i;
      pr	=gyer+k;
      pi	=gyei+k;
      dpr	=dgyer+k;
      dpi	=dgyei+k;
      for(k=0; k < ES2Mp1; k++){
	Ur[k]	=*EZvr++;
	Ui[k]	=*vi++;
	dUr[k]	=*dvr++;
	dUi[k]	=*dvi++;
      }
      for(m1=0; m1 < m; m1++){
	gar	=Xr[m1];
	gai	=Xi[m1];
	for(mm=0; mm < ES2Mp1; mm++){
	  *pr++		+=gar*Ur[mm]-gai*Ui[mm];
	  *pi++		+=gai*Ur[mm]+gar*Ui[mm];
	  *dpr++	+=gar*dUr[mm]-gai*dUi[mm];
	  *dpi++	+=gai*dUr[mm]+gar*dUi[mm];
	}
      }
      k		=M0Mm-ES2Mp1;
      EZvr	+=k;
      vi	+=k;
      dvr	+=k;
      dvi	+=k;
    }
  }
  return(0);
}

int eqreset(double *gyer,double *gyei,double *dgyer,double *dgyei)
{
  int k,ka;
  int i,m,m1,mm;
  double sr,si,gar,gai;
  double *pr,*pi,*dpr,*dpi,*EZvr,*vi,*dvr,*dvi;
  double *yr,*yi,*dyr,*dyi,*zr,*zi,*dzr,*dzi;

  ka		=M0Mm*ESNa;
  for(m=ESMp; m > 0; m--){
    k		=ka+ES2Mp1*m+ESMp+m;
    pr		=gyer+k;
    pi		=gyei+k;
    sr		=*pr;
    si		=*pi;
    gar		=1./(sr*sr+si*si);
    sr		*=gar;
    si		*=gar;
    k		=ka+ESMp-m;
    pr		=gyer+k;
    pi		=gyei+k;
    for(m1=0; m1 < ESMp1; m1++){
      Xr[m1]	= (*pi)*si-(*pr)*sr;
      Xi[m1]	=-(*pi)*sr-(*pr)*si;
      pr	+=ES2Mp1;
      pi	+=ES2Mp1;
    }
    pr		=gyer;
    pi		=gyei;
    dpr		=dgyer;
    dpi		=dgyei;
    yr          =aYer;
    yi          =aYei;
    dyr         =daYer;
    dyi         =daYei;
    k		=ES2Mp1*m;
    EZvr		=gyer+k;
    vi		=gyei+k;
    dvr		=dgyer+k;
    dvi		=dgyei+k;
    zr          =aYer+k;
    zi          =aYei+k;
    dzr         =daYer+k;
    dzi         =daYei+k;
    for(i=0; i < ESNa1; i++){
      for(k=0; k < ES2Mp1; k++){
	Ur[k]	=*EZvr++;
	Ui[k]	=*vi++;
	dUr[k]	=*dvr++;
	dUi[k]	=*dvi++;
      }
      for(m1=0; m1 < ESMp1; m1++){
	gar	=Xr[m1];
	gai	=Xi[m1];
	mm	=ES2Mp1;
	while(mm > 0){
	  mm--;
	  *pr++		+=gar*Ur[mm]+gai*Ui[mm];
	  *pi++		+=gai*Ur[mm]-gar*Ui[mm];
	  *dpr++	+=gar*dUr[mm]+gai*dUi[mm];
	  *dpi++	+=gai*dUr[mm]-gar*dUi[mm];
	}
      }
      for(k=0; k < ES2Mp1; k++){
	Ur[k]	=*zr++;
	Ui[k]	=*zi++;
	dUr[k]	=*dzr++;
	dUi[k]	=*dzi++;
      }
      for(m1=0; m1 < ESMp1; m1++){
	gar	=Xr[m1];
	gai	=Xi[m1];
	mm	=ES2Mp1;
	while(mm > 0){
	  mm--;
	  *yr++		+=gar*Ur[mm]+gai*Ui[mm];
	  *yi++		+=gai*Ur[mm]-gar*Ui[mm];
	  *dyr++	+=gar*dUr[mm]+gai*dUi[mm];
	  *dyi++	+=gai*dUr[mm]-gar*dUi[mm];
	}
      }
      k		=M0Mm-ES2Mp1;
      EZvr	+=k;
      vi	+=k;
      dvr	+=k;
      dvi	+=k;
      zr	+=k;
      zi	+=k;
      dzr	+=k;
      dzi	+=k;
    }
    k		=ka+ES2Mp1*m+ESMp+m;
    pr		=gyer+k;
    pi		=gyei+k;
    sr		=*pr;
    si		=-(*pi);
    gar		=1./(sr*sr+si*si);
    sr		*=gar;
    si		*=gar;
    k		=ka+ESMp+m;
    pr		=gyer+k;
    pi		=gyei+k;
    for(m1=0; m1 < m; m1++){
      Xr[m1]	= (*pi)*si-(*pr)*sr;
      Xi[m1]	=-(*pi)*sr-(*pr)*si;
      pr	+=ES2Mp1;
      pi	+=ES2Mp1;
    }
    k		=ES2Mp1*m;
    EZvr		=gyer+k;
    vi		=gyei+k;
    dvr		=dgyer+k;
    dvi		=dgyei+k;
    zr		=aYer+k;
    zi		=aYei+k;
    dzr		=daYer+k;
    dzi		=daYei+k;
    for(i=0; i < ESNa1; i++){
      k		=M0Mm*i;
      pr	=gyer+k;
      pi	=gyei+k;
      dpr	=dgyer+k;
      dpi	=dgyei+k;
      yr	=aYer+k;
      yi	=aYei+k;
      dyr	=daYer+k;
      dyi	=daYei+k;
      for(k=0; k < ES2Mp1; k++){
	Ur[k]	=*EZvr++;
	Ui[k]	=*vi++;
	dUr[k]	=*dvr++;
	dUi[k]	=*dvi++;
      }
      for(m1=0; m1 < m; m1++){
	gar	=Xr[m1];
	gai	=Xi[m1];
	for(mm=0; mm < ES2Mp1; mm++){
	  *pr++		+=gar*Ur[mm]-gai*Ui[mm];
	  *pi++		+=gai*Ur[mm]+gar*Ui[mm];
	  *dpr++	+=gar*dUr[mm]-gai*dUi[mm];
	  *dpi++	+=gai*dUr[mm]+gar*dUi[mm];
	}
      }
      for(k=0; k < ES2Mp1; k++){
	Ur[k]	=*zr++;
	Ui[k]	=*zi++;
	dUr[k]	=*dzr++;
	dUi[k]	=*dzi++;
      }
      for(m1=0; m1 < m; m1++){
	gar	=Xr[m1];
	gai	=Xi[m1];
	for(mm=0; mm < ES2Mp1; mm++){
	  *yr++		+=gar*Ur[mm]-gai*Ui[mm];
	  *yi++		+=gai*Ur[mm]+gar*Ui[mm];
	  *dyr++	+=gar*dUr[mm]-gai*dUi[mm];
	  *dyi++	+=gai*dUr[mm]+gar*dUi[mm];
	}
      }
      k		=M0Mm-ES2Mp1;
      EZvr	+=k;
      vi	+=k;
      dvr	+=k;
      dvi	+=k;
      zr	+=k;
      zi	+=k;
      dzr	+=k;
      dzi	+=k;
    }
  }
  pr	=gyer	+ESMp;
  pi	=gyei	+ESMp;
  dpr	=dgyer	+ESMp;
  dpi	=dgyei	+ESMp;
  yr	=aYer	+ESMp;
  yi	=aYei	+ESMp;
  dyr	=daYer	+ESMp;
  dyi	=daYei	+ESMp;
  for(i=0; i < ESNa1; i++){
    pi[0]	=0.;
    dpi[0]	=0.;
    yi[0]	=0.;
    dyi[0]	=0.;
    for(k=1; k < ESMp1; k++){
      mm	=-k;
      pr[k]	=EZcr2*(pr[k]+pr[mm]);
      pi[k]	=EZcr2*(pi[k]-pi[mm]);
      dpr[k]	=EZcr2*(dpr[k]+dpr[mm]);
      dpi[k]	=EZcr2*(dpi[k]-dpi[mm]);
      yr[k]	=EZcr2*(yr[k]+yr[mm]);
      yi[k]	=EZcr2*(yi[k]-yi[mm]);
      dyr[k]	=EZcr2*(dyr[k]+dyr[mm]);
      dyi[k]	=EZcr2*(dyi[k]-dyi[mm]);
      pr[mm]	= pr[k];
      pi[mm]	=-pi[k];
      dpr[mm]	= dpr[k];
      dpi[mm]	=-dpi[k];
      yr[mm]	= yr[k];
      yi[mm]	=-yi[k];
      dyr[mm]	= dyr[k];
      dyi[mm]	=-dyi[k];
    }
    pr	+=M0Mm;
    pi	+=M0Mm;
    dpr	+=M0Mm;
    dpi	+=M0Mm;
    yr	+=M0Mm;
    yi	+=M0Mm;
    dyr	+=M0Mm;
    dyi	+=M0Mm;
  }
  k	=ka;
  return(0);
}

static int NWRK=0,nWRK;
static double *RKy0,*RKy1,*RKy2,*RKy3,*RKy4;

int ECReInitRKSolv(int n)
{
  nWRK	=5*n;
  if(NWRK < nWRK){
    if(NWRK){
      free(RKy0);
    }
    RKy0	=(double*)malloc(nWRK*sizeof(double));
    RKy1	=RKy0	+n;
    RKy2	=RKy1	+n;
    RKy3	=RKy2	+n;
    RKy4	=RKy3	+n;
    NWRK	=nWRK;
  }
  return(0);
}

int ECDeInitRKSolv()
{
  if(NWRK){
    free(RKy0);
    NWRK	=0;
  }
  return(0);
}

#ifdef RKS
int RKSolvC(void (*AdvODE)(int *,double *, double *,double*),
	   double *Y, double *dY, double *X, double h, int neq)
#else
int RKSolvC(void *AdvODE, double *Y, double *dY, double *X, double h, int neq)
#endif
{
  double x,s;
  double w1,w2,w3,a21,a31,a32;
  int i;
  s	=h/(*X);
  w1	=EZcr6+(EZcr2+EZcr3*s)*s;
  w2	=EZcr3*(2.-4.*s);
  w3	=EZcr6+EZcr3*s;
  a21	=EZcr3/(s+sqrt(s*s+2.*EZcr3*w2));
  a31	=(EZcr2-1.5*s-2.*a21*w2)/w3;
  a32	=(s+a21*w2)/w3;
  w3	*=1.+s;
  s	=1.+EZcr2*s;
  w2	*=s;
  a32	*=s;
  x	=*X;

  for(i=0; i < neq; i++){
    RKy1[i]	=h*dY[i];
    RKy0[i]	=Y[i]+a21*RKy1[i];
  }
  x	=(*X)+EZcr2*h;
  AdvODE(&neq,&x,RKy0,RKy2);
  for(i=0; i < neq; i++){
    RKy2[i]	*=h;
    RKy0[i]	=Y[i]+a31*RKy1[i]+a32*RKy2[i];
  }
  (*X)		+=h;
  AdvODE(&neq,X,RKy0,RKy3);
  for(i=0; i < neq; i++){
    RKy3[i]	*=h;
    Y[i]	+=w1*RKy1[i]+w2*RKy2[i]+w3*RKy3[i];
  }
  AdvODE(&neq,X,Y,dY);
  return(0);
}

#ifdef RKS
int RKSolv(void (*AdvODE)(int *,double *, double *,double*),
	   double *Y, double *dY, double *X, double h, int neq)
#else
int RKSolv(void *AdvODE, double *Y, double *dY, double *X, double h, int neq)
#endif
{
  double x;
  int i;

  for(i=0; i < neq; i++){
    RKy1[i]	=h*dY[i];
    RKy0[i]	=Y[i]+EZcr2*RKy1[i];
  }
  x	=(*X)+EZcr2*h;
  AdvODE(&neq,&x,RKy0,RKy2);
  for(i=0; i < neq; i++){
    RKy2[i]	*=h;
    RKy0[i]	=Y[i]+EZcr2*RKy2[i];
  }
  AdvODE(&neq,&x,RKy0,RKy3);
  for(i=0; i < neq; i++){
    RKy3[i]	*=h;
    RKy0[i]	=Y[i]+RKy3[i];
  }
  (*X)		+=h;
  AdvODE(&neq,X,RKy0,RKy4);
  for(i=0; i < neq; i++){
    RKy4[i]	*=h;
    Y[i]	+=EZcr6*(RKy1[i]+RKy4[i])+EZcr3*(RKy2[i]+RKy3[i]);
  }
  AdvODE(&neq,X,Y,dY);
  return(0);
}

int  ECgY0Corr(int neq)
{
  double *pr,*pi,*yr,*yi,*pr0,*pi0,*yr0,*yi0;

  int mp,mm,m;
  double sr,si,tr,ti;

  pr0	=gpr0+ESMp;
  pi0	=gpi0+ESMp;
  yr0	=Yr0+ESMp;
  yi0	=Yi0+ESMp;

  sr	=pr0[0];
  si	=-pi0[0];
  tr	=1./(sr*sr+si*si);
  sr	*=tr;
  si	*=tr;
  pr0[0]	=1.;
  pi0[0]	=0.;
  tr	=yr0[0];
  ti	=yi0[0];
  yr0[0]	=tr*sr-ti*si;
  yi0[0]	=ti*sr+tr*si;
  mm	=0;
  for(mp=1; mp < ESMp1; mp++){
    mm--;
    tr	=pr0[mp];
    ti	=pi0[mp];
    pr0[mp]	=tr*sr-ti*si;
    pi0[mp]	=ti*sr+tr*si;
    tr	=yr0[mp];
    ti	=yi0[mp];
    yr0[mp]	=tr*sr-ti*si;
    yi0[mp]	=ti*sr+tr*si;
    tr	=pr0[mm];
    ti	=pi0[mm];
    pr0[mm]	=tr*sr-ti*si;
    pi0[mm]	=ti*sr+tr*si;
    tr	=yr0[mm];
    ti	=yi0[mm];
    yr0[mm]	=tr*sr-ti*si;
    yi0[mm]	=ti*sr+tr*si;
  }
  pr	=gpr+ESMp;
  pi	=pr+M0Mm;
  yr	=pi+M0Mm;
  yi	=yr+M0Mm;
  sr	=gpr[neq-1]-pr[0];
  si	=-pi[0];
  pr[0]	=gpr[neq-1];
  pi[0]	=0.;
  mm	=0;
  for(mp=1; mp < ESMp1; mp++){
    mm--;
    tr	=pr0[mp];
    ti	=pi0[mp];
    pr[mp]	+=tr*sr-ti*si;
    pi[mp]	+=ti*sr+tr*si;
    tr	=pr0[mm];
    ti	=pi0[mm];
    pr[mm]	=tr*sr-ti*si;
    pi[mm]	=ti*sr+tr*si;
    tr	=yr0[mp];
    ti	=yi0[mp];
    yr[mp]	=tr*sr-ti*si;
    yi[mp]	=ti*sr+tr*si;
    tr	=yr0[mm];
    ti	=yi0[mm];
    yr[mm]	=tr*sr-ti*si;
    yi[mm]	=ti*sr+tr*si;
  }    
  pr	+=ES2Mp1;
  pi	+=ES2Mp1;
  yr	+=ES2Mp1;
  yi	+=ES2Mp1;
  for(m=1; m < ESMp1; m++){
    sr	=-pr[0];
    si	=-pi[0];
    pr[0]	=0.;
    pi[0]	=0.;
    mm	=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      tr	=pr0[mp];
      ti	=pi0[mp];
      pr[mp]	+=tr*sr-ti*si;
      pi[mp]	+=ti*sr+tr*si;
      tr	=pr0[mm];
      ti	=pi0[mm];
      pr[mm]	=tr*sr-ti*si;
      pi[mm]	=ti*sr+tr*si;
      tr	=yr0[mp];
      ti	=yi0[mp];
      yr[mp]	=tr*sr-ti*si;
      yi[mp]	=ti*sr+tr*si;
      tr	=yr0[mm];
      ti	=yi0[mm];
      yr[mm]	=tr*sr-ti*si;
      yi[mm]	=ti*sr+tr*si;
    }    
    pr	+=ES2Mp1;
    pi	+=ES2Mp1;
    yr	+=ES2Mp1;
    yi	+=ES2Mp1;
  }
  return(0);
}

static int NRKCyl=0;
int RK2DStep(double *Y, double *dY, double *X, double h, int neq)
{
  double x,s;
  double w1,w2,w3,a21,a31,a32;
  int i;
  int k;

  NRKCyl++;
  s	=h/(*X);
  s	=0.;
  w1	=EZcr6+(EZcr2+EZcr3*s)*s;
  w2	=EZcr3*(2.-4.*s);
  w3	=EZcr6+EZcr3*s;
  a21	=EZcr3/(s+sqrt(s*s+2.*EZcr3*w2));
  a31	=(EZcr2-1.5*s-2.*a21*w2)/w3;
  a32	=(s+a21*w2)/w3;
  w3	*=1.+s;
  s	=1.+EZcr2*s;
  w2	*=s;
  a32	*=s;

  x	=*X;
  for(i=0; i < neq; i++){
    RKy1[i]	=h*dY[i];
    RKy0[i]	=Y[i]+a21*RKy1[i];
  }
  x	=(*X)+EZcr2*h;
  peq2D(&neq,&x,RKy0,RKy2);
  for(i=0; i < neq; i++){
    RKy2[i]	*=h;
    RKy0[i]	=Y[i]+a31*RKy1[i]+a32*RKy2[i];
  }
  (*X)		+=h;
  peq2D(&neq,X,RKy0,RKy3);
  for(i=0; i < neq; i++){
    RKy3[i]	*=h;
    Y[i]	+=w1*RKy1[i]+w2*RKy2[i]+w3*RKy3[i];
  }
  peq2D(&neq,X,Y,dY);
  return(0);
}

int EU2DStep(double *Y, double *dY, double *X, double h, int neq)
{
  double x,x0;
  int i;

  x0	=*X;
  for(i=0; i < neq; i++){
    RKy1[i]	=h*dY[i];
    RKy0[i]	=Y[i]+EZcr2*RKy1[i];
  }
  x	=(*X)+EZcr2*h;
  peq2D(&neq,&x,RKy0,RKy2);
  for(i=0; i < neq; i++){
    RKy2[i]	*=h;
    RKy0[i]	=Y[i]+2.*RKy2[i]-RKy1[i];
  }
  (*X)		+=h;
  peq2D(&neq,X,RKy0,RKy3);
  for(i=0; i < neq; i++){
    RKy3[i]	*=h;
    Y[i]	+=EZcr6*(RKy1[i]+4.*RKy2[i]+RKy3[i]);
  }
  peq2D(&neq,X,Y,dY);
  return(0);
}

int RK2DStep000(double *Y, double *dY, double *X, double h, int neq)
{
  double x,x0;
  int i;

  x0	=*X;
  for(i=0; i < neq; i++){
    RKy1[i]	=h*dY[i];
    RKy0[i]	=Y[i]+EZcr2*RKy1[i];
  }
  x	=(*X)+EZcr2*h;
  peq2D(&neq,&x,RKy0,RKy2);
  if(x0 < 0.06*ESsa[1]){
    int m,k,kk;
    double *p,*dp;
    
    p	=RKy0;
    dp=RKy2;
    for(m=1; m < 2; m++){
      kk	=ES2Mp1*m;
      printf("??? m=%2d x=%10.3e\n",m,x);
      for(k=0; k < ES2Mp1; k++){
	printf("??? gy[%3d]=%10.3e %10.3e %10.3e %10.3e\n",k-ESMp
	       ,gpr[kk]/(*X),gpr[kk+M0Mm]/(*X),dgpr[kk],dgpr[kk+M0Mm]);
	printf("*??  p[%3d]=%10.3e %10.3e %10.3e %10.3e\n",k-ESMp
	       ,p[kk],p[kk+M0Mm]/x,dp[kk]/x,dp[kk+M0Mm]);
	  kk++;
      }
    }
  }
  for(i=0; i < neq; i++){
    RKy2[i]	*=h;
    RKy0[i]	=Y[i]+EZcr2*RKy2[i];
  }
  peq2D(&neq,&x,RKy0,RKy3);
  for(i=0; i < neq; i++){
    RKy3[i]	*=h;
    RKy0[i]	=Y[i]+RKy3[i];
  }
  (*X)		+=h;
  peq2D(&neq,X,RKy0,RKy4);
  for(i=0; i < neq; i++){
    RKy4[i]	*=h;
    Y[i]	+=EZcr6*(RKy1[i]+RKy4[i])+EZcr3*(RKy2[i]+RKy3[i]);
  }
  peq2D(&neq,X,Y,dY);
  if(x0 < 0.06*ESsa[1]){
    int m,k,kk;
    double *p,*dp;

    p	=RKy0;
    dp=RKy2;
    for(m=1; m < 2; m++){
      kk	=ES2Mp1*m;
      printf("??? m=%2d x=%10.3e\n",m,x);
      for(k=0; k < ES2Mp1; k++){
	printf("??? gy[%3d]=%10.3e %10.3e %10.3e %10.3e\n",k-ESMp
	       ,gpr[kk]/(*X),gpr[kk+M0Mm]/(*X),dgpr[kk],dgpr[kk+M0Mm]);
	kk++;
      }
    }
  }
  return(0);
}

int RK2DStepE(double *Y, double *dY, double *X, double h, int neq)
{
  double x;
  int i;

  for(i=0; i < neq; i++){
    RKy1[i]	=h*dY[i];
    RKy0[i]	=Y[i]+EZcr2*RKy1[i];
  }
  for(i=0; i < ES2Mp1; i++){
    gpr0[i]	=EZcr2*h*dgpr0[i];
    gpi0[i]	=EZcr2*h*dgpi0[i];
    Yr0[i]	=EZcr2*h*dYr0[i];
    Yi0[i]	=EZcr2*h*dYi0[i];

    gpr1[i]	=EZcr6*h*dgpr0[i];
    gpi1[i]	=EZcr6*h*dgpi0[i];
    Yr1[i]	=EZcr6*h*dYr0[i];
    Yi1[i]	=EZcr6*h*dYi0[i];
  }
  x	=(*X)+EZcr2*h;
  peq2D(&neq,&x,RKy0,RKy2);
  for(i=0; i < neq; i++){
    RKy2[i]	*=h;
    RKy0[i]	=Y[i]+EZcr2*RKy2[i];
  }
  for(i=0; i < ES2Mp1; i++){
    gpr0[i]	=EZcr2*h*dgpr0[i];
    gpi0[i]	=EZcr2*h*dgpi0[i];
    Yr0[i]	=EZcr2*h*dYr0[i];
    Yi0[i]	=EZcr2*h*dYi0[i];

    gpr1[i]	+=EZcr3*h*dgpr0[i];
    gpi1[i]	+=EZcr3*h*dgpi0[i];
    Yr1[i]	+=EZcr3*h*dYr0[i];
    Yi1[i]	+=EZcr3*h*dYi0[i];
  }
  peq2D(&neq,&x,RKy0,RKy3);
  for(i=0; i < neq; i++){
    RKy3[i]	*=h;
    RKy0[i]	=Y[i]+RKy3[i];
  }
  /* Extra function */
  for(i=0; i < ES2Mp1; i++){
    gpr0[i]	=h*dgpr0[i];
    gpi0[i]	=h*dgpi0[i];
    Yr0[i]	=h*dYr0[i];
    Yi0[i]	=h*dYi0[i];

    gpr1[i]	+=EZcr3*gpr0[i];
    gpi1[i]	+=EZcr3*gpi0[i];
    Yr1[i]	+=EZcr3*Yr0[i];
    Yi1[i]	+=EZcr3*Yi0[i];
  }
  /* End of Extra function */
  (*X)		+=h;
  peq2D(&neq,X,RKy0,RKy4);
  for(i=0; i < neq; i++){
    RKy4[i]	*=h;
    Y[i]	+=EZcr6*(RKy1[i]+RKy4[i])+EZcr3*(RKy2[i]+RKy3[i]);
  }
  /* Extra function */
  printf("??? x=%10.3e\n",*X);
  for(i=0; i < ES2Mp1; i++){
    gpr0[i]	=gpr1[i]+EZcr6*h*dgpr0[i];
    gpi0[i]	=gpi1[i]+EZcr6*h*dgpi0[i];
    Yr0[i]	=Yr1[i]+EZcr6*h*dYr0[i];
    Yi0[i]	=Yi1[i]+EZcr6*h*dYi0[i];
    printf("??? gy0[%3d]=%10.3e %10.3e %10.3e %10.3e\n",i-ESMp
	   ,gpr0[i],gpi0[i],Yr0[i],Yi0[i]);
  }
  /* Calculate the corrected solution */
  if((ESEqSolvFl&0x0F) == 1){
    ECgY0Corr(neq);
  }
  for(i=0; i < ES2Mp1; i++){
    gpr0[i]	=0.;
    gpi0[i]	=0.;
    Yr0[i]	=0.;
    Yi0[i]	=0.;
  }
  gdx1		=*X;
  /* End of Extra function */
  peq2D(&neq,X,Y,dY);
  return(0);
}

extern double ESreltol,ESabstol;
static double hcur,xout,reltol=1.0e-6,abstol=1e-12;
static double*rwork=NULL,*yh,*acor,*vtol,*rvtol,eps0,eps1;
static int itol=4,itask,istate=1,iopt=0,mf=10,lrw,liw=20;
static int *iwork=NULL,Niwork=0,Nrwork=0,Nlsode=0;

int ECReInitLSODE(int n)
{
  lrw	=liw+16*n;
  if(Nlsode < n){
    if(iwork != NULL){
      free(iwork);
    }
    Niwork	=liw;
    iwork	=(int*)malloc(Niwork*sizeof(int));
    if(iwork == NULL){
      printf("Failure in memory allocation for iwork in ECReInitLSODE()\n");
      exit(0);
    }
    if(rwork != NULL){
      free(rwork);
    }
    Nrwork	=lrw+2*n;
    rwork	=(double*)malloc(Nrwork*sizeof(double));
    if(rwork == NULL){
      printf("Failure in memory allocation for rwork in ECReInitLSODE()\n");
    }
    Nlsode	=n;
  }
  return(0);
}

int ECDeInitLSODE()
{
  if(Niwork){
    free(iwork);
    iwork	=NULL;
    Niwork=0;
  }
  if(Nrwork){
    free(rwork);
    rwork	=NULL;
    Nrwork=0;
  }
  Nlsode	=0;
  return(0);
}

int EZq2js(double *js,double *djs,double jp,double djp,double t)
{
  double bL,dbL,d2bL,bV,dbV,d2bV,rF,bF,dF,d2F,bK,dbK,d2bK,gm,dgm,d2gm;
  ESSetSplA(t);
  ESSetSplDCr(t);
  
  splRA2(&bL,&dbL,&d2bL,ESLc,ESLc2);
  d2bL	=2.*dbL+t*d2bL;
  dbL	=bL+t*dbL;
  bL	*=t;
  splRA2(&bV,&dbV,&d2bV,ESVc,ESVc2);
  d2bV	=2.*dbV+t*d2bV;
  dbV	=bV+t*dbV;
  bV	*=t;
  splRA2(&bK,&dbK,&d2bK,ESg22c,ESg22c2);
  d2bK	=2.*dbK+t*d2bK;
  dbK	=bK+t*dbK;
  bK	*=t;
  splRDCr2(&gm,&dgm,&d2gm,7);
  splRA2(&bF,&dF,&d2F,ESFF,ESFF2a);
  bF	=sqrt(bF);
  rF	=1./bF;
  dF	*=EZcr2*rF;
  d2F	=(EZcr2*d2F-dF*dF)*rF;

  d2bK	=d2bK*bL+2.*dbK*dbL+bK*d2bL;
  dbK	=dbK*bL+bK*dbL;
  bK	*=bL;

  d2F	=(d2F*gm+2.*dF*dgm+bF*d2gm)*rR0;
  dF	=(dF*gm+bF*dgm)*rR0;
  bF	*=gm*rR0;

  d2bK	=d2bK*bF+2.*dbK*dF+bK*d2F;
  dbK	=dbK*bF+bK*dF;
  bK	*=bF;

  bL	=1./bL;
  dbL	*=-bL*bL;
  
  dbK	-=bV*jp;
  *js	=dbK*bL;
  *djs	=(d2bK-dbV*jp-bV*djp)*bL+dbK*dbL;

  return(0);
}

static clock_t Ect0;
int ECIfEqFound()
{
  NIter++;

  if(EcErr < ESEqSolvTol){
    ESEqTol	=EcErr;
    ESEqIt	=NIter;
    ESFail	=0;

#ifdef DEBUG
    printf("NIter=%d Err=%10.3e\n",NIter,EcErr);
    printf("<gb>_ext=%10.3e gb_N=%10.3e\n",ESgbext,ESgbN);
#endif
    NIter	=0;
    return(0);
  }

  if(EcCholFl > 0){
    ESEqTol	=EcErr;
    ESEqIt	=NIter;
    ESFail	=0x0FFFFFFF;
#ifdef DEBUG
    printf("ATTENTION !!! Too big jump in parameters%c\n",'\a');
    printf("NIter=%d Err=%10.3e\n",NIter,EcErr);
    printf("!! Nit=%d Err=%10.3e %2.2x %2.2x %2.2x %2.2x %2.2x %2.2x %2.2x\n"
	   ,ESEqIt,ESEqTol
	   ,(ESFail&0x0F000000)>>24
	   ,(ESFail&0x00F00000)>>20
	   ,(ESFail&0x000F0000)>>16
	   ,(ESFail&0x0000F000)>>12
	   ,(ESFail&0x00000F00)>>8
	   ,(ESFail&0x000000F0)>>4
	   , ESFail&0x0000000F);
#endif
    NIter	=0;
    kFail	=0;
    EcCholFl	=0;
    return(0);
  }

  if(kFail > 1){
    ESEqTol	=EcErr;
    ESEqIt	=NIter;
#ifdef DEBUG
    printf("ATTENTION !!! Too big jump in parameters%c\n",'\a');
    printf("NIter=%d Err=%10.3e\n",NIter,EcErr);
#endif
    NIter	=0;
    kFail	=0;
    if((ESFail&0x0F) < 0x0F){
      ESFail++;
      if((Fail&0x000000F0)){
	ESFail	+=0x00000010;
      }
      if((Fail&0x00000F00)){
	ESFail	+=0x00000100;
      }
      if((Fail&0x0000F000)){
	ESFail	+=0x00001000;
      }
      if((Fail&0x000F0000)){
	ESFail	+=0x00010000;
      }
      if((Fail&0x00F00000)){
	ESFail	+=0x00100000;
      }
      if((Fail&0x0F000000)){
	ESFail	+=0x01000000;
      }
    }
    if((ESFail&0x0000000F) < 2){
      ESFail	&=0x0F;
      ECSmooth2DCoord(0);
      return(1);
    }
    return(0);
  }
  
  if(ESEqSolvIt == 0){
    return(0);
  }

  if(NIter >= ESEqSolvIt){
    ESEqTol	=EcErr;
    ESEqIt	=NIter;
    if((ESFail&0x0F) < 0x0F){
      ESFail++;
      if((Fail&0x000000F0)){
	ESFail	+=0x00000010;
      }
      if((Fail&0x00000F00)){
	ESFail	+=0x00000100;
      }
      if((Fail&0x0000F000)){
	ESFail	+=0x00001000;
      }
      if((Fail&0x000F0000)){
	ESFail	+=0x00010000;
      }
      if((Fail&0x00F00000)){
	ESFail	+=0x00100000;
      }
      if((Fail&0x0F000000)){
	ESFail	+=0x01000000;
      }
    }
    NIter	=0;
#ifdef DEBUG
    printf("Nit=%d Err=%10.3e %2.2x %2.2x %2.2x %2.2x %2.2x %2.2x %2.2x"
	   " t=%10.3e sec\n"
	   ,ESEqIt,ESEqTol
	   ,(ESFail&0x0F000000)>>24
	   ,(ESFail&0x00F00000)>>20
	   ,(ESFail&0x000F0000)>>16
	   ,(ESFail&0x0000F000)>>12
	   ,(ESFail&0x00000F00)>>8
	   ,(ESFail&0x000000F0)>>4
	   , ESFail&0x0000000F
	   ,(double)(clock()-Ect0)/CLOCKS_PER_SEC);
#endif
    
    if((ESFail&0x0000000F) < 2){
      ECSmooth2DCoord(0);
      ESFail	&=0x0F;
      return(1);
    }
    return(0);
  }
  return(1);
}

int EC2DEqSolv(int n)
{
  double x,h,sr,si;
  int neq,kst,Kst=1,nst;
  int i,ii,k,kk,level,ki,k1,k2,k3,m,mm,km0;
  
  Fail	=0;
  if(NIter == 0){
    Ect0=clock();
  }
  rR0	=1./R0;
  x	=ESa0 == 0. ? ESsa[1]*.25 : ESa0;
  gdx1	=x;
  neq	=4*M0Mm+2;
  EcCholFl	=0; /* Temporary solution */
  switch((ESEqSolvFl&0xF0)){
  case 0x00:
  case 0x10:
  case 0x30:
    switch(ESEqSolvInCr){
    case 0:
    case 3:
      peq2D	=eq2DVsj;
      if(ESa0 == 0.) Start2D(x);
      else ECStartEq2DsjHole();
      break;
    case 1:
    case 2:
    case 4:
      peq2D	=eq2DVjb;
      if(ESa0 == 0.){
	static int ic=0;
	ic++;
#ifdef H
	int i,m,j,k;
	double s,sr,si,tr,ti,am;
	double K0,Kr,Ki;
	double Cr,Ci,Br,Bi,A;
	K0	= ESg22c[0];
	Kr	= ESg22c[2*ESNa1];
	Ki	=-ESg22s[2*ESNa1];
	EZout("sddddd","K0=",K0,Kr,Ki,Kr/K0,Ki/K0);

	Start2Daa(x);
	i	=0;
	am	=x*x;
	for(m=0; m < ESMp1; m++){
	  if(m == 1) am	=x;
	  putchar('\n');
	  for(k=-ESMp; k < ESMp1; k++){
	    s	=(m-k)*(m-k+2);
	    Cr	=s*Kr;
	    Ci	=s*Ki;
	    A	=(m*m-k*k)*K0;
	    s	=(m+k)*(m+k+2);
	    Br	= s*Kr;
	    Bi	=-s*Ki;
	    j	=i-2;
	    if(k-2 < -ESMp){
	      sr	=0.;
	      si	=0.;
	    }
	    else{
	      sr	=gpr[j];
	      si	=gpi[j];
	    }
	    j	=i+2;
	    if(k+2 > ESMp){
	      tr	=0.;
	      ti	=0.;
	    }
	    else{
	      tr	=gpr[j];
	      ti	=gpi[j];
	    }
	    EZout("siidddddd","E",m,k,gpr[i]/am,gpi[i]/am,Yr[i]/am,Yi[i]/am
		,(Cr*sr-Ci*si+A*gpr[i]+Br*tr-Bi*ti)/am
		,(Ci*sr+Cr*si+A*gpi[i]+Bi*tr+Br*ti)/am);

	    i++;
	  }
	  am	*=x;
	}
#endif
	Start2D(x);
#ifdef H
	i	=0;
	am	=x*x;
	for(m=0; m < ESMp1; m++){
	  if(m == 1) am	=x;
	  putchar('\n');
	  for(k=-ESMp; k < ESMp1; k++){
	    s	=(m-k)*(m-k+2);
	    Cr	=s*Kr;
	    Ci	=s*Ki;
	    A	=(m*m-k*k)*K0;
	    s	=(m+k)*(m+k+2);
	    Br	= s*Kr;
	    Bi	=-s*Ki;
	    j	=i-2;
	    if(k-2 < -ESMp){
	      sr	=0.;
	      si	=0.;
	    }
	    else{
	      sr	=gpr[j];
	      si	=gpi[j];
	    }
	    j	=i+2;
	    if(k+2 > ESMp){
	      tr	=0.;
	      ti	=0.;
	    }
	    else{
	      tr	=gpr[j];
	      ti	=gpi[j];
	    }
	    EZout("siidddddd","e",m,k,gpr[i]/am,gpi[i]/am,Yr[i]/am,Yi[i]/am
		,(Cr*sr-Ci*si+A*gpr[i]+Br*tr-Bi*ti)/am
		,(Ci*sr+Cr*si+A*gpi[i]+Bi*tr+Br*ti)/am);
	    i++;
	  }
	  am	*=x;
	}
#endif
      }
      else{
	ECStartEq2DsjHole();
      }
      break;
    case 6:
    case 7:
      peq2D	=eq2DVgm;
      if(ESa0 == 0.){
	int i,j,k;
	k	=0;
	ECStartEq2Dsq(x);
#ifdef H
	Start2Dgm(x);
	for(i=0; i < ESMp1; i++){
	  putchar('\n');
	  for(j=-ESMp; j < ESMp1; j++){
	    EZout("siidddd","+?? gr",i,j,gpr[k],gpi[k],dgpr[k],dgpi[k]);
	    k++;
	  }
	}
	k	=0;
	for(i=0; i < ESMp1; i++){
	  putchar('\n');
	  for(j=-ESMp; j < ESMp1; j++){
	    EZout("siidddd","??? gr",i,j,gpr[k],gpi[k],dgpr[k],dgpi[k]);
	    k++;
	  }
	}
	Start2Dq(x);
#endif
      }
      else{
	ECStartEq2DsjHole();
      }
      break;
    case 8:
      peq2D	=eq2DVgY;
      Start2DgY(x);
      break;
    }
    break;
  default:
    break;
  }
  switch((ESEqSolvFl&0x0F)){
  case 0:
    ECReInitLSODE(neq);
    yh		=rwork+liw;
    acor	=rwork+lrw-n;
    vtol	=rwork+lrw;
    rvtol	=vtol+n;
    itask	=5;
    itask	=4;
    rwork[0]	=ESsa[ESNa];
    istate	=1;
    iopt	=1;
    rwork[4]	=.02;
    rwork[5]	=ESsa[ESNa];
    rwork[6]	=0.;
    iwork[4]	=12;
    iwork[5]	=500;
    iwork[6]	=10;
    break;
  case 1:
    ECReInitRKSolv(neq);
    break;
  }
  if(ESa0 == 0.){
    m	=ESNa1*M0Mm;
    peq2D(&neq,&x,gpr,dgpr);
    for(k=0; k < M0Mm; k++){
      EZgper[k]	=0.;
      EZgpei[k]	=0.;
      EZdgper[k]	=0.;
      EZdgpei[k]	=0.;
      aYer[k]	=0.;
      aYei[k]	=0.;
      daYer[k]	=0.;
      daYei[k]	=0.;
    }
    k	=ESMp;
    EZdgper[k]	=dgpr[k]/x;
    EZdgpei[k]	=dgpi[k]/x;
    k	+=ES2Mp1;
    EZdgper[k]	=dgpr[k]/x;
    EZdgpei[k]	=dgpi[k]/x;
    k	+=ES2Mp1;
    EZdgper[k]	=dgpr[k]/x;
    EZdgpei[k]	=dgpi[k]/x;
    k	=ES2Mp1+1+ESMp;
    EZdgper[k]	=1.;
    EZdgpei[k]	=0.;
  }
  else{
    peq2D(&neq,ESsa,gpr,dgpr);
    for(k=0; k < M0Mm; k++){
      EZgper[k]	=0.;
      EZgpei[k]	=0.;
      EZdgper[k]	=dgpr[k];
      EZdgpei[k]	=dgpi[k];
    }
  }
  ESgY0[0]	=0.;
  switch((ESEqSolvFl&0x0F)){
#ifndef NO_LSODE
  case 0:
    reltol	=ESreltol;
    abstol	=ESabstol;
    i		=1;
    vtol[neq-1]	=0.;
    vtol[neq-2]	=0.;
    rvtol[neq-1]=eps_rel;
    rvtol[neq-2]=eps_rel;
    itol	=1;
    while(i < ESNa1){
      xout=ESsa[i];
      for(m=0;m<ESMp1;m++){
	km0= ESMp+m;
	k= ES2Mp1*m+km0;
	if(gpr[km0]!=0.)
	  eps0= eps_rel*fabs(gpr[ESMp]/gpr[km0]);
	else
	  eps0= 1.;
	if(eps0>1.)
	  eps0= 1.;
	if(Yr[km0]!=0.)
	  eps1= eps_rel*fabs(Yr[ESMp]/Yr[km0]);
	else
	  eps1= 1.;
	if(eps1>1.)
	  eps1= 1.;
	if(i==1){
	  eps0= 1.;
	  eps1= 1.;
	}
	eps0	=eps0*fabs(gpr[k]);
	eps1	=eps1*fabs(Yr[k]);
	k	=ES2Mp1*m;
	k1	=k+M0Mm;
	k2	=k1+M0Mm;
	k3	=k2+M0Mm;
	for(mm=0; mm < ES2Mp1; mm++){
	  kk	=abs(km0-mm);
	  if(kk > ESMp)
	    kk=ESMp;
	  vtol[k]	=eps0;
	  vtol[k1]	=eps0;
	  vtol[k2]	=eps1;
	  vtol[k3]	=eps1;
	  rvtol[k]	=eps_rel;
	  rvtol[k1]	=eps_rel;
	  rvtol[k2]	=eps_rel;
	  rvtol[k3]	=eps_rel;
	  k++;
	  k1++;
	  k2++;
	  k3++;
	}
      }
#ifdef h
      lsode_(peq2D,&neq,gpr,&x,&xout,&itol,rvtol,vtol,&itask,
	     &istate,&iopt,rwork,&lrw,iwork,&liw,NULL,&mf);
#else
      lsode_(peq2D,&neq,gpr,&x,&xout,&itol,&reltol,&abstol,&itask,
	     &istate,&iopt,rwork,&lrw,iwork,&liw,NULL,&mf);
#endif
      peq2D(&neq,&xout,gpr,dgpr);
      gdx1	=ESsa[i];
      ESaY[i]	=gpr[neq-2]/ESpa[i];
      ESgY0[i]	=gpr[neq-1];
      ESgY02a[i]=dgpr[neq-1]/ESsa[i];
      ESdgY[i]	=dgpr[neq-1]/ESsa[i];
      level	=0;
      if(fabs(gpr[ESMp+ESMp]/Yi[M0Mm+1]) > Eps_noise){
	level	=1;
      }
      if(fabs(gpr[ES2Mp1+ESMp+2]/gpr[ES2Mp1+ESMp+1]) > 10.*Eps_noise){
	level	=1;
      }
      k1	=M0Mm*i;
      for(k=0; k < M0Mm; k++){
	EZgper[k1]	=gpr[k];
	EZgpei[k1]	=gpi[k];
	EZdgper[k1]	=dgpr[k];
	EZdgpei[k1]	=dgpi[k];
	aYer[k1]	=Yr[k];
	aYei[k1]	=Yi[k];
	daYer[k1]	=dYr[k];
	daYei[k1]	=dYi[k];
	k1++;
      }
      if((0 && i == 3) || level && x > 0.2){
#ifdef DEBUG
	printf("Redistribution i=%2d\n",i);
#endif
	EZresete(gpr,yh,iwork+14,EZgper,EZgpei,EZdgper,EZdgpei,i);
      }
      if(istate < 0){
	printf("lsode istate=%d\n",istate);
	exit(0);
      }
      i++;
    }
    break;
#endif
  case 1:
    peq2D(&neq,&x,gpr,dgpr);
    Kst	=4;
    if(ESa0 == 0.){
      i		=1;
      kst	=ESMp/3;
      h		=(ESsa[1]-x)/kst;
      if(h > EZcr2*x){
	nst	=2.*h/x;
	if(nst < 2) nst	=2;
	kst	*=nst;
	h	=(ESsa[1]-x)/kst;
      }
      kst	=4;
      h	=(ESsa[1]-x)/kst;
      while(kst > 0){
	kst--;
	RK2DStep(gpr,dgpr,&x,h,neq);
#ifdef H
	RKSolv(peq2D,gpr,dgpr,&x,h,neq);
#endif
      }
      ESaY[i]	=gpr[neq-2]/ESpa[i];
      ESgY0[i]	=gpr[neq-1];
      ESgY02a[i]=dgpr[neq-1]/x;
      ESdgY[i]	=dgpr[neq-1]/x;
      kk	=M0Mm*i;
      for(k=0; k < M0Mm; k++,kk++){
	EZgper[kk]	=gpr[k];
	EZgpei[kk]	=gpi[k];
	EZdgper[kk]	=dgpr[k];
	EZdgpei[kk]	=dgpi[k];
	aYer[kk]	=Yr[k];
	aYei[kk]	=Yi[k];
	daYer[kk]	=dYr[k];
	daYei[kk]	=dYi[k];
      }
    }
    else{
      i	=0;
    }

    while(i < ESNa){
      i++;
      nst	=Kst;
      nst	=(ESsa[i]-x)*ESMp > ESsa[i]*Kst ? 
	(int)((ESsa[i]-x)*ESMp/ESsa[i]) : Kst;
      nst	=6-i;
      if(nst < 2){
	nst	=2;
      }
      if(0 && i == ESNa){
	nst	=3;
      }

      nst	=2;
      if(i > ESNa-2) nst=4;
      h	=(ESsa[i]-x)/nst;
      while(nst > 0){
	nst--;
	RK2DStep(gpr,dgpr,&x,h,neq);
#ifdef H
	RKSolv(peq2D,gpr,dgpr,&x,h,neq);
#endif
	if(EcCholFl){
	  return(1);
	}
      }

      ESaY[i]	=gpr[neq-2]/ESpa[i];
      ESgY0[i]	=gpr[neq-1];
      ESgY02a[i]=dgpr[neq-1]/x;
      ESdgY[i]	=dgpr[neq-1]/x;
    
      level	=0;
      if(fabs(gpr[ESMp+ESMp]/Yi[M0Mm+1]) > Eps_noise){
	level	=1;
      }
      if(fabs(gpr[ES2Mp1+2*ESMp]/Yi[M0Mm+1]) > Eps_noise){
	level	=2;
      }
      kk	=M0Mm*i;
      for(k=0; k < M0Mm; k++,kk++){
	EZgper[kk]	=gpr[k];
	EZdgper[kk]	=dgpr[k];
	EZgpei[kk]	=gpi[k];
	EZdgpei[kk]	=dgpi[k];
	aYer[kk]	=Yr[k];
	aYei[kk]	=Yi[k];
	daYer[kk]	=dYr[k];
	daYei[kk]	=dYi[k];
      }
      if((0 && i == 3) || level && x > 0.1){
#ifdef DEBUG
	printf("Redistribution i=%2d %2d\n",i,level);
#endif
	reseteRK(i);
      }
    }
    reseteRK(ESNa);
    NRKCyl=0;
    break;
  }
#ifdef H
  eqreset(EZgper,EZgpei,EZdgper,EZdgpei);
#endif
  k	=ESMp;
  ESgY02a[0]	=EZdgper[k];
  ESdgY[0]	=EZdgper[k];
  for(i=1; i < ESNa1; i++){
    k	=ESMp+M0Mm*i;
    ESgY0[i]	=EZgper[k];
    ESgY02a[i]	=EZdgper[k]/ESsa[i];
    ESdgY[i]	=EZdgper[k]/ESsa[i];
    ESaY[i]	=aYer[k]/ESpa[i];
  }

#ifdef H
  ESgY2a1a	=dYr[ESMp]/Yr[ESMp]-ESg22c1[ESNa]/ESg22c[ESNa]-1.;
#endif
  k	=ESMp+M0Mm*ESNa;
  ESgY2a1a	=EZdgper[k]/EZgper[k]-ESg22c1[ESNa]/ESg22c[ESNa]-1.;
#ifdef H
  k		=neq-2;
  ESgY2a1a	=dgpr[k]/gpr[k]-ESg22c1[ESNa]/ESg22c[ESNa]-1.;
  ESdgY1a[0]	=ESsa[0] == 0.? 
    0. : (ESLc[0]*ESjs[0]+ESVc[0]*ESjp[0])/(ESa0*ESg22c[0]);
  ESdgY1a[ESNa]	=(ESgY2a1a-1.)*ESdgY[ESNa];
#endif
  splAA(ESdgY,ESdgY1a,ESdgY2a,ESdgY1a,ESdgY1a+ESNa);
  return(0);
}
#endif/*stg_2DSolv*/

#ifndef stg_2DMatching
int ECgYcs2dRZcsBound()
{
  double *gyr,*gyi,*gyar,*gyai,*gxr,*gxi,*gsr,*gsi,*ur,*ui;
  double *pCc,*pCs,*pSc,*pSs;
  double *dCc,*dCs,*dSc,*dSs;

  double Rtr[ESFp1],Rti[ESFp1],Rar[ESFp1],Rai[ESFp1];
  double Ratr[ESFp1],Rati[ESFp1],Raar[ESFp1],Raai[ESFp1];
  double rgY01a,za,zaa,b,ba,baa,rb,barb;

  double sr,si,cr,ci,rr,ri;
  int mv,m,k,kk,mp,mm;

  kk		=ESNa;
  Rar [0]	=rcT1a[kk];
  Raar[0]	=rcT2a[kk];
  za		=2.*rsT1a[kk];
  zaa		=2.*rsT2a[kk];
  for(k=1; k < ESFp1; k++){
    kk		+=ESNa1;
    Rtr [k]	=k*rsT  [kk];
    Rti [k]	=k*rcT  [kk];
    Ratr[k]	=k*rsT1a[kk];
    Rati[k]	=k*rcT1a[kk];
    Rar [k]	=  rcT1a[kk];
    Rai [k]	= -rsT1a[kk];
    Raar[k]	=  rcT2a[kk];
    Raai[k]	= -rsT2a[kk];
  }

  ur	=Ur+ESMp1;
  ui	=Ui+ESMp1;
  gxr	=Xr+ESMp;
  gxi	=Xi+ESMp;
  gsr	=Yr+ESMp;
  gsi	=Yi+ESMp;

  k	=M0Mm*ESNa+ESMp;
  gyr	=EZgper+k;
  gyi	=EZgpei+k;
  gyar	=EZdgper+k;
  gyai	=EZdgpei+k;
  pCc	=rCc;
  pCs	=rCs;
  pSc	=rSc;
  pSs	=rSs;
  dCc	=drCc;
  dCs	=drCs;
  dSc	=drSc;
  dSs	=drSs;

  rgY01a=1./ESgY02a[ESNa];

  k	=ESNa;
  b	=ESsb[k];
  ba	=ESsb1a[k];
  baa	=ESsb2a[k];
  rb	=1./b;
  barb	=ba*rb;
  for(mv=0; mv < ESMp1; mv++){
    gxr[0]	=0.;
    gxi[0]	=0.;
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      gxr[mp]	=-gyr[mp]*rgY01a;
      gxi[mp]	=-gyi[mp]*rgY01a;
      mm--;
      gxr[mm]	=-gyr[mm]*rgY01a;
      gxi[mm]	=-gyi[mm]*rgY01a;
    }
    /* Determination of $\gs$ */
    for(m=-ESMp; m < ESMp1; m++){
      ur[m]	=-za*gxr[m];
      ui[m]	=-za*gxi[m];
    }
    ur[ESMp1]	=0.;
    ui[ESMp1]	=0.;
    ur[-ESMp1]	=0.;
    ui[-ESMp1]	=0.;
    k	=0;
    mp	=0;
    mm	=0;
    for(m=-ESMp; m < ESMp1; m++){
      mp	=m+1;
      cr	=ba*gxi[m];
      ci	=ba*gxr[m];
      ur[mp]	-=cr;
      ui[mp]	+=ci;
      mm	=m-1;
      ur[mm]	+=cr;
      ui[mm]	-=ci;
    }
    mp	=ESMp;
    mm	=-ESMp;
    m	=mp+1;
    k	=mm-1;
    while(mp > 0){
      gsr[mp]	=rb*ur[m];
      gsi[mp]	=rb*ui[m];
      mp--;
      ur[mp]	-=ur[m];
      ui[mp]	-=ui[m];
      m--;
      gsr[mm]	=rb*ur[k];
      gsi[mm]	=rb*ui[k];
      mm++;
      ur[mm]	-=ur[k];
      ui[mm]	-=ui[k];
      k++;
    }
    pCs[0]	=za*gxr[0]-EZcr2*ur[0];
    pSs[0]	=za*gxi[0]-EZcr2*ui[0];
    bCc[mv]	=ba*gxr[0]+EZcr2*(ui[1]-ui[-1]);
    bSc[mv]	=ba*gxi[0]+EZcr2*(ur[-1]-ur[1]);
    gsr[0]	=EZcr2*rb*(ur[1]+ur[-1]);
    gsi[0]	=EZcr2*rb*(ui[1]+ui[-1]);
    /* Contribution of $r'_a\gx$ */
    cr		=gxr[0];
    sr		=gxi[0];
    pCc[0]	=Rar[0]*cr;
    pSc[0]	=Rar[0]*sr;
    for(k=1; k < Mr; k++){
      rr	=Rar[k];
      ri	=Rai[k];
      pCc[k]	=rr*cr;
      pCs[k]	=-ri*cr;
      pSc[k]	=rr*sr;
      pSs[k]	=-ri*sr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gxr[mp]+gxr[mm]);
      ci	=EZcr2*(gxi[mp]-gxi[mm]);
      sr	=EZcr2*(gxi[mp]+gxi[mm]);
      si	=EZcr2*(gxr[mm]-gxr[mp]);
      rr	=Rar[0];
      m		=mp;
      if(m < Mr){
	pCc[m]	+=rr*cr;
	pCs[m]	-=rr*ci;
	pSc[m]	+=rr*sr;
	pSs[m]	-=rr*si;
	m++;
	k	=1;
	while(m < Mr){
	  rr		=Rar[k];
	  ri		=Rai[k];
	  pCc[m]	+=rr*cr-ri*ci;
	  pCs[m]	-=ri*cr+rr*ci;
	  pSc[m]	+=rr*sr-ri*si;
	  pSs[m]	-=ri*sr+rr*si;
	  m++;
	  k++;
	}
      }
      m		=mp;
      k		=0;
      while(m > 0){
 	m--;
	k++;
	rr	=Rar[k];
	ri	=-Rai[k];
	pCc[m]	+=rr*cr-ri*ci;
	pCs[m]	-=ri*cr+rr*ci;
	pSc[m]	+=rr*sr-ri*si;
	pSs[m]	-=ri*sr+rr*si;
      }
      k		=mp;
      m		=0;
      kk	=ESFp1-k;
      if(kk > Mr){
	kk	=Mr;
      }
      while(kk > 0){
	rr	=Rar[k];
	ri	=Rai[k];
	pCc[m]	+=rr*cr+ri*ci;
	pCs[m]	-=ri*cr-rr*ci;
	pSc[m]	+=rr*sr+ri*si;
	pSs[m]	-=ri*sr-rr*si;
 	m++;
	k++;
	kk--;
      }
    }
    /* Contribution of $r'_\gt\gs$ */
    cr		=gsr[0];
    sr		=gsi[0];
    for(k=1; k < Mr; k++){
      rr	=Rtr[k];
      ri	=Rti[k];
      pCc[k]	+=rr*cr;
      pCs[k]	-=ri*cr;
      pSc[k]	+=rr*sr;
      pSs[k]	-=ri*sr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gsr[mp]+gsr[mm]);
      ci	=EZcr2*(gsi[mp]-gsi[mm]);
      sr	=EZcr2*(gsi[mp]+gsi[mm]);
      si	=EZcr2*(gsr[mm]-gsr[mp]);
      m		=mp+1;
      k		=1;
      while(m < Mr){
	rr	=Rtr[k];
	ri	=Rti[k];
	pCc[m]	+=rr*cr-ri*ci;
	pCs[m]	-=ri*cr+rr*ci;
	pSc[m]	+=rr*sr-ri*si;
	pSs[m]	-=ri*sr+rr*si;
	m++;
	k++;
      }
      m		=mp;
      k		=0;
      while(m > 0){
 	m--;
	k++;
	rr	=Rtr[k];
	ri	=-Rti[k];
	pCc[m]	+=rr*cr-ri*ci;
	pCs[m]	-=ri*cr+rr*ci;
	pSc[m]	+=rr*sr-ri*si;
	pSs[m]	-=ri*sr+rr*si;
      }
      k		=mp;
      m		=0;
      kk	=ESFp1-k;
      if(kk > Mr){
	kk	=Mr;
      }
      while(kk > 0){
	rr	=Rtr[k];
	ri	=Rti[k];
	pCc[m]	+=rr*cr+ri*ci;
	pCs[m]	-=ri*cr-rr*ci;
	pSc[m]	+=rr*sr+ri*si;
	pSs[m]	-=ri*sr-rr*si;
 	m++;
	k++;
	kk--;
      }
    }
    /* Derivatives: Contribution of $r''_{aa}\gx$ */
    cr		=gxr[0];
    sr		=gxi[0];
    dCc[0]	=Raar[0]*cr;
    dCs[0]	=0.;
    dSc[0]	=Raar[0]*sr;
    dSs[0]	=0.;
    for(k=1; k < Mr; k++){
      rr	=Raar[k];
      ri	=Raai[k];
      dCc[k]	=rr*cr;
      dCs[k]	=-ri*cr;
      dSc[k]	=rr*sr;
      dSs[k]	=-ri*sr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gxr[mp]+gxr[mm]);
      ci	=EZcr2*(gxi[mp]-gxi[mm]);
      sr	=EZcr2*(gxi[mp]+gxi[mm]);
      si	=EZcr2*(gxr[mm]-gxr[mp]);
      rr	=Raar[0];
      m		=mp;
      if(m < Mr){
	dCc[m]	+=rr*cr;
	dCs[m]	-=rr*ci;
	dSc[m]	+=rr*sr;
	dSs[m]	-=rr*si;
	m++;
	k	=1;
	while(m < Mr){
	  rr	=Raar[k];
	  ri	=Raai[k];
	  dCc[m]	+=rr*cr-ri*ci;
	  dCs[m]	-=ri*cr+rr*ci;
	  dSc[m]	+=rr*sr-ri*si;
	  dSs[m]	-=ri*sr+rr*si;
	  m++;
	  k++;
	}
      }
      m		=mp;
      k		=0;
      while(m > 0){
 	m--;
	k++;
	rr	=Raar[k];
	ri	=-Raai[k];
	dCc[m]	+=rr*cr-ri*ci;
	dCs[m]	-=ri*cr+rr*ci;
	dSc[m]	+=rr*sr-ri*si;
	dSs[m]	-=ri*sr+rr*si;
      }
      k		=mp;
      m		=0;
      kk	=ESFp1-k;
      if(kk > Mr){
	kk	=Mr;
      }
      while(kk > 0){
	rr	=Raar[k];
	ri	=Raai[k];
	dCc[m]	+=rr*cr+ri*ci;
	dCs[m]	-=ri*cr-rr*ci;
	dSc[m]	+=rr*sr+ri*si;
	dSs[m]	-=ri*sr-rr*si;
 	m++;
	k++;
	kk--;
      }
    }
    /* Contribution of $r''_{a\gt}\gs$ */
    cr		=gsr[0];
    sr		=gsi[0];
    for(k=1; k < Mr; k++){
      rr	=Ratr[k];
      ri	=Rati[k];
      dCc[k]	+=rr*cr;
      dCs[k]	-=ri*cr;
      dSc[k]	+=rr*sr;
      dSs[k]	-=ri*sr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gsr[mp]+gsr[mm]);
      ci	=EZcr2*(gsi[mp]-gsi[mm]);
      sr	=EZcr2*(gsi[mp]+gsi[mm]);
      si	=EZcr2*(gsr[mm]-gsr[mp]);
      m		=mp+1;
      k		=1;
      while(m < Mr){
	rr	=Ratr[k];
	ri	=Rati[k];
	dCc[m]	+=rr*cr-ri*ci;
	dCs[m]	-=ri*cr+rr*ci;
	dSc[m]	+=rr*sr-ri*si;
	dSs[m]	-=ri*sr+rr*si;
	m++;
	k++;
      }
      m		=mp;
      k		=0;
      while(m > 0){
 	m--;
	k++;
	rr	=Ratr[k];
	ri	=-Rati[k];
	dCc[m]	+=rr*cr-ri*ci;
	dCs[m]	-=ri*cr+rr*ci;
	dSc[m]	+=rr*sr-ri*si;
	dSs[m]	-=ri*sr+rr*si;
      }
      k		=mp;
      m		=0;
      kk	=ESFp1-k;
      if(kk > Mr){
	kk	=Mr;
      }
      while(kk > 0){
	rr	=Ratr[k];
	ri	=Rati[k];
	dCc[m]	+=rr*cr+ri*ci;
	dCs[m]	-=ri*cr-rr*ci;
	dSc[m]	+=rr*sr+ri*si;
	dSs[m]	-=ri*sr-rr*si;
 	m++;
	k++;
	kk--;
      }
    }
    /* Contribution of derivatives $\gx'_a$ $\gs'_a$ */
    for(mp=-ESMp; mp < ESMp1; mp++){
      gsr[mp]	*=-barb;
      gsi[mp]	*=-barb;
      ur[mp]	=-zaa*gxr[mp];
      ui[mp]	=-zaa*gxi[mp];
    }
    ur[-ESMp1]	=0.;
    ui[-ESMp1]	=0.;
    ur[ESMp1]	=0.;
    ui[ESMp1]	=0.;
    for(m=-ESMp; m < ESMp1; m++){
      mp	=m+1;
      mm	=m-1;
      ur[mp]	-=baa*gxi[m];
      ui[mp]	+=baa*gxr[m];
      ur[mm]	+=baa*gxi[m];
      ui[mm]	-=baa*gxr[m];
    }
    mp	=ESMp;
    mm	=-ESMp;
    m	=mp+1;
    k	=mm-1;
    while(mp > 0){
      gsr[mp]	+=rb*ur[m];
      gsi[mp]	+=rb*ui[m];
      mp--;
      ur[mp]	-=ur[m];
      ui[mp]	-=ui[m];
      m--;
      gsr[mm]	+=rb*ur[k];
      gsi[mm]	+=rb*ui[k];
      mm++;
      ur[mm]	-=ur[k];
      ui[mm]	-=ui[k];
      k++;
    }
    dCs[0]	=-EZcr2*ur[0];
    dSs[0]	=-EZcr2*ui[0];
    dbCc[mv]	=EZcr2*(ui[1]-ui[-1]);
    dbSc[mv]	=EZcr2*(ur[-1]-ur[1]);
    gsr[0]	+=EZcr2*rb*(ur[1]+ur[-1]);
    gsi[0]	+=EZcr2*rb*(ui[1]+ui[-1]);
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      gxr[mp]	=-ESgY2a1a*gxr[mp]-gyar[mp]*rgY01a;
      gxi[mp]	=-ESgY2a1a*gxi[mp]-gyai[mp]*rgY01a;
      mm--;
      gxr[mm]	=-ESgY2a1a*gxr[mm]-gyar[mm]*rgY01a;
      gxi[mm]	=-ESgY2a1a*gxi[mm]-gyai[mm]*rgY01a;
    }
    /* Determination of $\gs'_a$ */
    for(m=-ESMp-1; m < ESMp1+1; m++){
      ur[m]	=0.;
      ui[m]	=0.;
    }
    for(m=-ESMp; m < ESMp1; m++){
      ur[m]	-=za*gxr[m];
      ui[m]	-=za*gxi[m];
      mp	=m+1;
      mm	=m-1;
      ur[mp]	-=ba*gxi[m];
      ui[mp]	+=ba*gxr[m];
      ur[mm]	+=ba*gxi[m];
      ui[mm]	-=ba*gxr[m];
    }
    mp	=ESMp;
    mm	=-ESMp;
    m	=mp+1;
    k	=mm-1;
    while(mp > 0){
      gsr[mp]	+=rb*ur[m];
      gsi[mp]	+=rb*ui[m];
      mp--;
      ur[mp]	-=ur[m];
      ui[mp]	-=ui[m];
      m--;
      gsr[mm]	+=rb*ur[k];
      gsi[mm]	+=rb*ui[k];
      mm++;
      ur[mm]	-=ur[k];
      ui[mm]	-=ui[k];
      k++;
    }
#ifdef H
    gxr[0]	=(EZcr2*(ui[-1]-ui[1])-dbCc[mv])/ba;
    gxi[0]	=(EZcr2*(ur[1]-ur[-1])-dbSc[mv])/ba;
    gxr[0]	=0.;
    gxi[0]	=0.;
#endif

    dCs[0]	+=za*gxr[0]-EZcr2*ur[0];
    dSs[0]	+=za*gxi[0]-EZcr2*ui[0];
    dbCc[mv]	+=ba*gxr[0]+EZcr2*(ui[1]-ui[-1]);
    dbSc[mv]	+=ba*gxi[0]+EZcr2*(ur[-1]-ur[1]);
    gsr[0]	+=EZcr2*rb*(ur[1]+ur[-1]);
    gsi[0]	+=EZcr2*rb*(ui[1]+ui[-1]);

    /* Contribution of $r'_a\gx'_a$ */
    cr		=gxr[0];
    sr		=gxi[0];
    dCc[0]	+=Rar[0]*cr;
    dSc[0]	+=Rar[0]*sr;
    for(k=1; k < Mr; k++){
      rr	=Rar[k];
      ri	=Rai[k];
      dCc[k]	+=rr*cr;
      dCs[k]	-=ri*cr;
      dSc[k]	+=rr*sr;
      dSs[k]	-=ri*sr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gxr[mp]+gxr[mm]);
      ci	=EZcr2*(gxi[mp]-gxi[mm]);
      sr	=EZcr2*(gxi[mp]+gxi[mm]);
      si	=EZcr2*(gxr[mm]-gxr[mp]);
      rr	=Rar[0];
      m		=mp;
      if(m < Mr){
	dCc[m]	+=rr*cr;
	dCs[m]	-=rr*ci;
	dSc[m]	+=rr*sr;
	dSs[m]	-=rr*si;
	m++;
	k	=1;
	while(m < Mr){
	  rr		=Rar[k];
	  ri		=Rai[k];
	  dCc[m]	+=rr*cr-ri*ci;
	  dCs[m]	-=ri*cr+rr*ci;
	  dSc[m]	+=rr*sr-ri*si;
	  dSs[m]	-=ri*sr+rr*si;
	  m++;
	  k++;
	}
      }
      m		=mp;
      k		=0;
      while(m > 0){
 	m--;
	k++;
	rr	=Rar[k];
	ri	=-Rai[k];
	dCc[m]	+=rr*cr-ri*ci;
	dCs[m]	-=ri*cr+rr*ci;
	dSc[m]	+=rr*sr-ri*si;
	dSs[m]	-=ri*sr+rr*si;
      }
      k		=mp;
      m		=0;
      kk	=ESFp1-k;
      if(kk > Mr){
	kk	=Mr;
      }
      while(kk > 0){
	rr	=Rar[k];
	ri	=Rai[k];
	dCc[m]	+=rr*cr+ri*ci;
	dCs[m]	-=ri*cr-rr*ci;
	dSc[m]	+=rr*sr+ri*si;
	dSs[m]	-=ri*sr-rr*si;
 	m++;
	k++;
	kk--;
      }
    }
    /* Contribution of $r'_\gt\gs'_a$ */
    cr		=gsr[0];
    sr		=gsi[0];
    for(k=1; k < Mr; k++){
      rr	=Rtr[k];
      ri	=Rti[k];
      dCc[k]	+=rr*cr;
      dCs[k]	-=ri*cr;
      dSc[k]	+=rr*sr;
      dSs[k]	-=ri*sr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gsr[mp]+gsr[mm]);
      ci	=EZcr2*(gsi[mp]-gsi[mm]);
      sr	=EZcr2*(gsi[mp]+gsi[mm]);
      si	=EZcr2*(gsr[mm]-gsr[mp]);
      m		=mp+1;
      k		=1;
      while(m < Mr){
	rr	=Rtr[k];
	ri	=Rti[k];
	dCc[m]	+=rr*cr-ri*ci;
	dCs[m]	-=ri*cr+rr*ci;
	dSc[m]	+=rr*sr-ri*si;
	dSs[m]	-=ri*sr+rr*si;
	m++;
	k++;
      }
      m		=mp;
      k		=0;
      while(m > 0){
 	m--;
	k++;
	rr	=Rtr[k];
	ri	=-Rti[k];
	dCc[m]	+=rr*cr-ri*ci;
	dCs[m]	-=ri*cr+rr*ci;
	dSc[m]	+=rr*sr-ri*si;
	dSs[m]	-=ri*sr+rr*si;
      }
      k		=mp;
      m		=0;
      kk	=ESFp1-k;
      if(kk > Mr){
	kk	=Mr;
      }
      while(kk > 0){
	rr	=Rtr[k];
	ri	=Rti[k];
	dCc[m]	+=rr*cr+ri*ci;
	dCs[m]	-=ri*cr-rr*ci;
	dSc[m]	+=rr*sr+ri*si;
	dSs[m]	-=ri*sr-rr*si;
 	m++;
	k++;
	kk--;
      }
    }
    gyr	+=ES2Mp1;
    gyi	+=ES2Mp1;
    gyar+=ES2Mp1;
    gyai+=ES2Mp1;
    pCc	+=Mr;
    pCs	+=Mr;
    pSc	+=Mr;
    pSs	+=Mr;
    dCc	+=Mr;
    dCs	+=Mr;
    dSc	+=Mr;
    dSs	+=Mr;
  }
  return(0);
}

int ECgYrc2dRZcs()
{
  double *gyr,*gyi,*gxr,*gxi,*gsr,*gsi,*ur,*ui;

  double Rtr[ESFp1],Rti[ESFp1],Rar[ESFp1],Rai[ESFp1];
  double rgY01a,za,b,ba,rb;
  double cr,ci,rr,ri;
  int i,m,k,kk,mp,mm;

  rgY01a	=1./ESgY02a[0];
  dz0[0]	=2.*ECb1a[0]*EZdgpei[ESMp1]*rgY01a;
  drc[0]	=-4.*(EZd1rcs[ESNa1]*EZdgper[ESMp1]-EZd1rsn[ESNa1]*EZdgpei[ESMp1])
    *rgY01a;
  drs[0]	=0.;
  dbT[0]	=0.;
  for(m=1; m < ESMp1; m++){
    drc[ESNa1*m]	=0.;
    drs[ESNa1*m]	=0.;
  }
  ur	=Ur+ESMp1;
  ui	=Ui+ESMp1;
  gxr	=Xr+ESMp;
  gxi	=Xi+ESMp;
  gsr	=Yr+ESMp;
  gsi	=Yi+ESMp;

  k	=ESMp;
  gyr	=EZgper+k;
  gyi	=EZgpei+k;
  for(i=1; i < ESNa1; i++){/*All magnetic surfaces*/
    gyr		+=M0Mm;
    gyi		+=M0Mm;
    rgY01a	=1./(ESsa[i]*ESgY02a[i]);
    b		=ESsb  [i];
    ba		=ESsb1a[i];
    za		=2.*rsT1a[i];
    rb		=1./b;
    kk		=i;
    Rar[0]	=   rcT1a[kk];
    for(k=1; k < ESFp1; k++){
      kk	+=ESNa1;
      Rtr[k]	=k*rsT  [kk];
      Rti[k]	=k*rcT  [kk];
      Rar[k]	=  rcT1a[kk];
      Rai[k]	= -rsT1a[kk];
    }

    /* Calculation of $\gx$ */
    gxr[0]	=0.;
    gxi[0]	=0.;
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      gxr[mp]	=-gyr[mp]*rgY01a;
      gxi[mp]	=-gyi[mp]*rgY01a;
      mm--;
      gxr[mm]	=-gyr[mm]*rgY01a;
      gxi[mm]	=-gyi[mm]*rgY01a;
    }
    /* Calculation of $\gs$ */
    for(m=-ESMp; m < ESMp1; m++){
      ur[m]	=-za*gxr[m];
      ui[m]	=-za*gxi[m];
    }
    ur[-ESMp1]	=0.;
    ui[-ESMp1]	=0.;
    ur[ESMp1]	=0.;
    ui[ESMp1]	=0.;

    for(m=-ESMp; m < ESMp1; m++){
      mp	=m+1;
      mm	=m-1;
      ur[mp]	-=ba*gxi[m];
      ui[mp]	+=ba*gxr[m];
      ur[mm]	+=ba*gxi[m];
      ui[mm]	-=ba*gxr[m];
    }
    mp	=ESMp;
    mm	=-ESMp;
    m	=mp+1;
    k	=mm-1;
    while(mp > 0){
      gsr[mp]	=rb*ur[m];
      gsi[mp]	=rb*ui[m];
      mp--;
      ur[mp]	-=ur[m];
      ui[mp]	-=ui[m];
      m--;
      gsr[mm]	=rb*ur[k];
      gsi[mm]	=rb*ui[k];
      mm++;
      ur[mm]	-=ur[k];
      ui[mm]	-=ui[k];
      k++;
    }
#ifdef H
    gxr[0]	=EZcr2*(ui[-1]-ui[1])/ba;
    gxi[0]	=0.;
#endif
    dz0[i]	=za*gxr[0]-EZcr2*ur[0];
    dbT[i]	=ba*gxr[0]+EZcr2*(ui[1]-ui[-1]);
    gsr[0]	=EZcr2*rb*(ur[1]+ur[-1]);
    gsi[0]	=EZcr2*rb*(ui[1]+ui[-1]);

    /* Contribution of $r'_a\gx$ to $\gd r_m$ with non-negative $m$*/
    kk		=i;
    cr		=gxr[0];
    drc[kk]	=Rar[0]*cr;
    drs[kk]	=0.;
    for(k=1; k < Mr; k++){
      kk	+=ESNa1;
      rr	=Rar[k];
      ri	=Rai[k];
      drc[kk]	=rr*cr;
      drs[kk]	=-ri*cr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gxr[mp]+gxr[mm]); 
      ci	=EZcr2*(gxi[mp]-gxi[mm]);
      rr	=Rar[0];
      m		=mp;
      if(m < Mr){
	kk	=m*ESNa1+i;
	drc[kk]	+=rr*cr;
	drs[kk]	-=rr*ci;
	m++;
	k	=1;
	while(m < Mr){
	  kk	+=ESNa1;
	  rr	=Rar[k];
	  ri	=Rai[k];
	  drc[kk]	+=rr*cr-ri*ci;
	  drs[kk]	-=ri*cr+rr*ci;
	  m++;
	  k++;
	}
      }
      m		=mp;
      kk	=m*ESNa1+i;
      k		=0;
      while(m > 0){
 	m--;
	kk	-=ESNa1;
	k++;
	rr	=Rar[k];
	ri	=-Rai[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
      }
      k		=mp;
      kk	=i;
      m		=ESFp1-mp;
      if(m > Mr){
	m	=Mr;
      }
      while(m > 0){
	rr	=Rar[k];
	ri	=Rai[k];
	drc[kk]	+=rr*cr+ri*ci;
	drs[kk]	-=ri*cr-rr*ci;
	kk	+=ESNa1;
	k++;
	m--;
      }
    }
    /* Contribution of $r'_\gt\gs$ */
    cr		=gsr[0];
    kk		=i;
    for(k=1; k < Mr; k++){
      kk	+=ESNa1;
      drc[kk]	+=Rtr[k]*cr;
      drs[kk]	-=Rti[k]*cr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gsr[mp]+gsr[mm]);
      ci	=EZcr2*(gsi[mp]-gsi[mm]);
      m		=mp+1;
      kk	=m*ESNa1+i;
      k		=1;
      while(m < Mr){
	rr	=Rtr[k];
	ri	=Rti[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
	m++;
	kk	+=ESNa1;
	k++;
      }
      m		=mp;
      kk	=m*ESNa1+i;
      k		=0;
      while(m > 0){
 	m--;
	kk	-=ESNa1;
	k++;
	rr	=Rtr[k];
	ri	=-Rti[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
      }
      k		=mp;
      kk	=i;
      m	=ESFp1-k;
      if(m > Mr){
	m	=Mr;
      }
      while(m > 0){
	rr	=Rtr[k];
	ri	=Rti[k];
	drc[kk]	+=rr*cr+ri*ci;
	drs[kk]	-=ri*cr-rr*ci;
	kk	+=ESNa1;
	k++;
	m--;
      }
    }
  }
  return(0);
}

int ECgYrc2dRZcsHole1()
{
  double *gyr,*gyi,*gxr,*gxi,*gsr,*gsi,*ur,*ui;

  double Rtr[ESFp1],Rti[ESFp1],Rar[ESFp1],Rai[ESFp1];
  double rgY01a,za,b,ba,rb;
  double cr,ci,rr,ri;
  int i,m,k,kk,mp,mm,k1;

  ur	=Ur+ESMp1;
  ui	=Ui+ESMp1;
  gxr	=Xr+ESMp;
  gxi	=Xi+ESMp;
  gsr	=Yr+ESMp;
  gsi	=Yi+ESMp;

  k	=ESMp;
  gyr	=EZdgper+k;
  gyi	=EZdgpei+k;

  i	=0;
  kk		=0;
  Rar[0]	=   rcT1a[0];
  za		=2.*rsT1a[0];
  for(k=1; k < ESFp1; k++){
    kk	+=ESNa1;
    Rtr[k]	=k*rsT  [kk];
    Rti[k]	=k*rcT  [kk];
    Rar[k]	=  rcT1a[kk];
    Rai[k]	= -rsT1a[kk];
  }
  rgY01a	=-ESg22c[0]/(ESLc[0]*ESjs[0]+ESVc[0]*ESjp[0]);
  b		=ESsb  [0];
  ba		=ESsb1a[0];
  rb		=1./b;
  gxr[0]	=0.;
  gxi[0]	=0.;
  mm		=0;
  for(mp=1; mp < ESMp1; mp++){
    gxr[mp]	=-gyr[mp]*rgY01a;
    gxi[mp]	=-gyi[mp]*rgY01a;
    mm--;
    gxr[mm]	=-gyr[mm]*rgY01a;
    gxi[mm]	=-gyi[mm]*rgY01a;
  }
  /* Determination of $\gs$ */
  ur[0]	=0.;
  ui[0]	=0.;
  mp		=1;
  mm		=-1;
  while(mp < ESMp1){
    ur[mp]	=-za*gxr[mp];
    ui[mp]	=-za*gxi[mp];
    ur[mm]	=-za*gxr[mm];
    ui[mm]	=-za*gxi[mm];
    mp++;
    mm--;
  }
  ur[mp]	=0.;
  ui[mp]	=0.;
  ur[mm]	=0.;
  ui[mm]	=0.;
  k	=0;
  mp	=0;
  mm	=0;
  for(m=0; m < ESMp1; m++){
    mp++;
    mm--;
    ur[mp]	-=ba*gxi[m];
    ui[mp]	+=ba*gxr[m];
    ur[mm]	+=ba*gxi[k];
    ui[mm]	-=ba*gxr[k];
    k--;
  }
  k	=0;
  mp	=0;
  mm	=0;
  for(m=1; m < ESMp1; m++){
    k--;
    ur[mp]	+=ba*gxi[m];
    ui[mp]	-=ba*gxr[m];
    ur[mm]	-=ba*gxi[k];
    ui[mm]	+=ba*gxr[k];
    mp++;
    mm--;
  }
  mp	=ESMp;
  mm	=-ESMp;
  m	=mp+1;
  k	=mm-1;
  while(mp > 0){
    gsr[mp]	=rb*ur[m];
    gsi[mp]	=rb*ui[m];
    mp--;
    ur[mp]	-=ur[m];
    ui[mp]	-=ui[m];
    m--;
    gsr[mm]	=rb*ur[k];
    gsi[mm]	=rb*ui[k];
    mm++;
    ur[mm]	-=ur[k];
    ui[mm]	-=ui[k];
    k++;
  }
  gxr[0]	=-gyr[0]*rgY01a;
  gxi[0]	=-gyi[0]*rgY01a;
  gxr[0]	=EZcr2*(ui[-1]-ui[1])/ba;
  gxi[0]	=0.;
  dz0[0]	=za*gxr[0]-EZcr2*ur[0];
  dbT[0]	=ba*gxr[0]+EZcr2*(ui[1]-ui[-1]);
  gsr[0]	=EZcr2*rb*(ur[1]+ur[-1]);
  gsi[0]	=EZcr2*rb*(ui[1]+ui[-1]);
  /* Contribution of $r'_a\gx$ */
  kk		=i;
  cr		=gxr[0];
  drc[kk]	=Rar[0]*cr;
  drs[kk]	=0.;
  for(k=1; k < Mr; k++){
    kk	+=ESNa1;
    rr	=Rar[k];
    ri	=Rai[k];
    drc[kk]	=rr*cr;
    drs[kk]	=-ri*cr;
  }
  mm		=0;
  for(mp=1; mp < ESMp1; mp++){
    mm--;
    cr	=EZcr2*(gxr[mp]+gxr[mm]); 
    ci	=EZcr2*(gxi[mp]-gxi[mm]);
    rr	=Rar[0];
    m		=mp;
    if(m < Mr){
      kk	=m*ESNa1+i;
      drc[kk]	+=rr*cr;
      drs[kk]	-=rr*ci;
      m++;
      k	=1;
      while(m < Mr){
	kk	+=ESNa1;
	rr	=Rar[k];
	ri	=Rai[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
	m++;
	k++;
      }
    }
    m		=mp;
    kk	=m*ESNa1+i;
    k		=0;
    while(m > 0){
      m--;
      kk	-=ESNa1;
      k++;
      rr	=Rar[k];
      ri	=-Rai[k];
      drc[kk]	+=rr*cr-ri*ci;
      drs[kk]	-=ri*cr+rr*ci;
    }
    k		=mp;
    m		=0;
    kk	=i;
    k1	=ESFp1-k;
    if(k1 > Mr){
      k1	=Mr;
    }
    while(k1 > 0){
      rr	=Rar[k];
      ri	=Rai[k];
      drc[kk]	+=rr*cr+ri*ci;
      drs[kk]	-=ri*cr-rr*ci;
      m++;
      kk	+=ESNa1;
      k++;
      k1--;
    }
  }
  /* Contribution of $r'_\gt\gs$ */
  cr		=gsr[0];
  kk		=i;
  for(k=1; k < Mr; k++){
    kk	+=ESNa1;
    drc[kk]	+=Rtr[k]*cr;
    drs[kk]	-=Rti[k]*cr;
  }
  mm		=0;
  for(mp=1; mp < ESMp1; mp++){
    mm--;
    cr	=EZcr2*(gsr[mp]+gsr[mm]);
    ci	=EZcr2*(gsi[mp]-gsi[mm]);
    m		=mp+1;
    kk	=m*ESNa1+i;
    k		=1;
    while(m < Mr){
      rr	=Rtr[k];
      ri	=Rti[k];
      drc[kk]	+=rr*cr-ri*ci;
      drs[kk]	-=ri*cr+rr*ci;
      m++;
      kk	+=ESNa1;
      k++;
    }
    m		=mp;
    kk	=m*ESNa1+i;
    k		=0;
    while(m > 0){
      m--;
      kk	-=ESNa1;
      k++;
      rr	=Rtr[k];
      ri	=-Rti[k];
      drc[kk]	+=rr*cr-ri*ci;
      drs[kk]	-=ri*cr+rr*ci;
    }
    k		=mp;
    m		=0;
    kk	=i;
    k1	=ESFp1-k;
    if(k1 > Mr){
      k1	=Mr;
    }
    while(k1 > 0){
      rr	=Rtr[k];
      ri	=Rti[k];
      drc[kk]	+=rr*cr+ri*ci;
      drs[kk]	-=ri*cr-rr*ci;
      m++;
      kk	+=ESNa1;
      k++;
      k1--;
    }
  }

  k	=ESMp;
  gyr	=EZgper+k;
  gyi	=EZgpei+k;
  for(i=1; i < ESNa1; i++){
    gyr		+=M0Mm;
    gyi		+=M0Mm;
    kk		=i;
    Rar[0]	=   rcT1a[kk];
    za		=2.*rsT1a[kk];
    for(k=1; k < ESFp1; k++){
      kk	+=ESNa1;
      Rtr[k]	=k*rsT  [kk];
      Rti[k]	=k*rcT  [kk];
      Rar[k]	=  rcT1a[kk];
      Rai[k]	= -rsT1a[kk];
    }
    rgY01a	=1./(ESsa[i]*ESgY02a[i]);
    b		=ESsb  [i];
    ba		=ESsb1a[i];
    rb		=1./b;
    gxr[0]	=0.;
    gxi[0]	=0.;
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      gxr[mp]	=-gyr[mp]*rgY01a;
      gxi[mp]	=-gyi[mp]*rgY01a;
      mm--;
      gxr[mm]	=-gyr[mm]*rgY01a;
      gxi[mm]	=-gyi[mm]*rgY01a;
    }
    /* Determination of $\gs$ */
    ur[0]	=-za*gxr[0];
    ui[0]	=-za*gxi[0];
    mp		=1;
    mm		=-1;
    while(mp < ESMp1){
      ur[mp]	=-za*gxr[mp];
      ui[mp]	=-za*gxi[mp];
      ur[mm]	=-za*gxr[mm];
      ui[mm]	=-za*gxi[mm];
      mp++;
      mm--;
    }
    ur[mp]	=0.;
    ui[mp]	=0.;
    ur[mm]	=0.;
    ui[mm]	=0.;
    k	=0;
    mp	=0;
    mm	=0;
    for(m=0; m < ESMp1; m++){
      mp++;
      mm--;
      ur[mp]	-=ba*gxi[m];
      ui[mp]	+=ba*gxr[m];
      ur[mm]	+=ba*gxi[k];
      ui[mm]	-=ba*gxr[k];
      k--;
    }
    k	=0;
    mp	=0;
    mm	=0;
    for(m=1; m < ESMp1; m++){
      k--;
      ur[mp]	+=ba*gxi[m];
      ui[mp]	-=ba*gxr[m];
      ur[mm]	-=ba*gxi[k];
      ui[mm]	+=ba*gxr[k];
      mp++;
      mm--;
    }
    mp	=ESMp;
    mm	=-ESMp;
    m	=mp+1;
    k	=mm-1;
    while(mp > 0){
      gsr[mp]	=rb*ur[m];
      gsi[mp]	=rb*ui[m];
      mp--;
      ur[mp]	-=ur[m];
      ui[mp]	-=ui[m];
      m--;
      gsr[mm]	=rb*ur[k];
      gsi[mm]	=rb*ui[k];
      mm++;
      ur[mm]	-=ur[k];
      ui[mm]	-=ui[k];
      k++;
    }
    gxr[0]	=-gyr[0]*rgY01a;
    gxi[0]	=-gyi[0]*rgY01a;
    gxr[0]	=EZcr2*(ui[-1]-ui[1])/ba;
    gxi[0]	=0.;
    dz0[i]	=za*gxr[0]-EZcr2*ur[0];
    dbT[i]	=ba*gxr[0]+EZcr2*(ui[1]-ui[-1]);
    gsr[0]	=EZcr2*rb*(ur[1]+ur[-1]);
    gsi[0]	=EZcr2*rb*(ui[1]+ui[-1]);
    /* Contribution of $r'_a\gx$ */
    kk		=i;
    cr		=gxr[0];
    drc[kk]	=Rar[0]*cr;
    drs[kk]	=0.;
    for(k=1; k < Mr; k++){
      kk	+=ESNa1;
      rr	=Rar[k];
      ri	=Rai[k];
      drc[kk]	=rr*cr;
      drs[kk]	=-ri*cr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gxr[mp]+gxr[mm]); 
      ci	=EZcr2*(gxi[mp]-gxi[mm]);
      rr	=Rar[0];
      m		=mp;
      if(m < Mr){
	kk	=m*ESNa1+i;
	drc[kk]	+=rr*cr;
	drs[kk]	-=rr*ci;
	m++;
	k	=1;
	while(m < Mr){
	  kk	+=ESNa1;
	  rr	=Rar[k];
	  ri	=Rai[k];
	  drc[kk]	+=rr*cr-ri*ci;
	  drs[kk]	-=ri*cr+rr*ci;
	  m++;
	  k++;
	}
      }
      m		=mp;
      kk	=m*ESNa1+i;
      k		=0;
      while(m > 0){
 	m--;
	kk	-=ESNa1;
	k++;
	rr	=Rar[k];
	ri	=-Rai[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
      }
      k		=mp;
      m		=0;
      kk	=i;
      k1	=ESFp1-k;
      if(k1 > Mr){
	k1	=Mr;
      }
      while(k1 > 0){
	rr	=Rar[k];
	ri	=Rai[k];
	drc[kk]	+=rr*cr+ri*ci;
	drs[kk]	-=ri*cr-rr*ci;
 	m++;
	kk	+=ESNa1;
	k++;
	k1--;
      }
    }
    /* Contribution of $r'_\gt\gs$ */
    cr		=gsr[0];
    kk		=i;
    for(k=1; k < Mr; k++){
      kk	+=ESNa1;
      drc[kk]	+=Rtr[k]*cr;
      drs[kk]	-=Rti[k]*cr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gsr[mp]+gsr[mm]);
      ci	=EZcr2*(gsi[mp]-gsi[mm]);
      m		=mp+1;
      kk	=m*ESNa1+i;
      k		=1;
      while(m < Mr){
	rr	=Rtr[k];
	ri	=Rti[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
	m++;
	kk	+=ESNa1;
	k++;
      }
      m		=mp;
      kk	=m*ESNa1+i;
      k		=0;
      while(m > 0){
 	m--;
	kk	-=ESNa1;
	k++;
	rr	=Rtr[k];
	ri	=-Rti[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
      }
      k		=mp;
      m		=0;
      kk	=i;
      k1	=ESFp1-k;
      if(k1 > Mr){
	k1	=Mr;
      }
      while(k1 > 0){
	rr	=Rtr[k];
	ri	=Rti[k];
	drc[kk]	+=rr*cr+ri*ci;
	drs[kk]	-=ri*cr-rr*ci;
 	m++;
	kk	+=ESNa1;
	k++;
	k1--;
      }
    }
  }
  return(0);
}

int ECgYrc2dRZcsHole(int n)
{
  double *gyr,*gyi,*gxr,*gxi,*gsr,*gsi,*ur,*ui;

  double Rtr[ESFp1],Rti[ESFp1],Rar[ESFp1],Rai[ESFp1];
  double rgY01a,za,b,ba,rb;
  double cr,ci,rr,ri;
  int i,in,m,k,kk,mp,mm,k1,I,K,M;

  K	=ESnAF*n;
  M	=ESnMp*n;
  I	=ESNa1*n;
  in	=I;

  ur	=Ur+ESMp1;
  ui	=Ui+ESMp1;
  gxr	=Xr+ESMp;
  gxi	=Xi+ESMp;
  gsr	=Yr+ESMp;
  gsi	=Yi+ESMp;
  
  k	=ESMp;
  gyr	=EZdgper+k;
  gyi	=EZdgpei+k;
  for(i=0; i < 1; i++){
    kk		=K+i;
    Rar[0]	=   rcT1a[kk];
    za		=2.*rsT1a[kk];
    for(k=1; k < ESFp1; k++){
      kk	+=ESNa1;
      Rtr[k]	=k*rsT  [kk];
      Rti[k]	=k*rcT  [kk];
      Rar[k]	=  rcT1a[kk];
      Rai[k]	= -rsT1a[kk];
    }
    rgY01a	=-ESg22c[0]/(ESLc[0]*ESjs[0]+ESVc[0]*ESjp[0]);
    b		=ESsb  [in];
    ba		=ESsb1a[in];
    rb		=1./b;
    gxr[0]	=0.;
    gxi[0]	=0.;
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      gxr[mp]	=-gyr[mp]*rgY01a;
      gxi[mp]	=-gyi[mp]*rgY01a;
      mm--;
      gxr[mm]	=-gyr[mm]*rgY01a;
      gxi[mm]	=-gyi[mm]*rgY01a;
    }
    /* Determination of $\gs$ */
    ur[0]	=-za*gxr[0];
    ui[0]	=-za*gxi[0];
    mp		=1;
    mm		=-1;
    while(mp < ESMp1){
      ur[mp]	=-za*gxr[mp];
      ui[mp]	=-za*gxi[mp];
      ur[mm]	=-za*gxr[mm];
      ui[mm]	=-za*gxi[mm];
      mp++;
      mm--;
    }
    ur[mp]	=0.;
    ui[mp]	=0.;
    ur[mm]	=0.;
    ui[mm]	=0.;
    k	=0;
    mp	=0;
    mm	=0;
    for(m=0; m < ESMp1; m++){
      mp++;
      mm--;
      ur[mp]	-=ba*gxi[m];
      ui[mp]	+=ba*gxr[m];
      ur[mm]	+=ba*gxi[k];
      ui[mm]	-=ba*gxr[k];
      k--;
    }
    k	=0;
    mp	=0;
    mm	=0;
    for(m=1; m < ESMp1; m++){
      k--;
      ur[mp]	+=ba*gxi[m];
      ui[mp]	-=ba*gxr[m];
      ur[mm]	-=ba*gxi[k];
      ui[mm]	+=ba*gxr[k];
      mp++;
      mm--;
    }
    mp	=ESMp;
    mm	=-ESMp;
    m	=mp+1;
    k	=mm-1;
    while(mp > 0){
      gsr[mp]	=rb*ur[m];
      gsi[mp]	=rb*ui[m];
      mp--;
      ur[mp]	-=ur[m];
      ui[mp]	-=ui[m];
      m--;
      gsr[mm]	=rb*ur[k];
      gsi[mm]	=rb*ui[k];
      mm++;
      ur[mm]	-=ur[k];
      ui[mm]	-=ui[k];
      k++;
    }
    gxr[0]	=-gyr[0]*rgY01a;
    gxi[0]	=-gyi[0]*rgY01a;
    gxr[0]	=EZcr2*(ui[-1]-ui[1])/ba;
    gxi[0]	=0.;
    dz0[i]	=za*gxr[0]-EZcr2*ur[0];
    dbT[in]	=ba*gxr[0]+EZcr2*(ui[1]-ui[-1]);
    gsr[0]	=EZcr2*rb*(ur[1]+ur[-1]);
    gsi[0]	=EZcr2*rb*(ui[1]+ui[-1]);
    /* Contribution of $r'_a\gx$ */
    kk		=i;
    cr		=gxr[0];
    drc[kk]	=Rar[0]*cr;
    drs[kk]	=0.;
    for(k=1; k < Mr; k++){
      kk	+=ESNa1;
      rr	=Rar[k];
      ri	=Rai[k];
      drc[kk]	=rr*cr;
      drs[kk]	=-ri*cr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gxr[mp]+gxr[mm]); 
      ci	=EZcr2*(gxi[mp]-gxi[mm]);
      rr	=Rar[0];
      m		=mp;
      if(m < Mr){
	kk	=m*ESNa1+i;
	drc[kk]	+=rr*cr;
	drs[kk]	-=rr*ci;
	m++;
	k	=1;
	while(m < Mr){
	  kk	+=ESNa1;
	  rr	=Rar[k];
	  ri	=Rai[k];
	  drc[kk]	+=rr*cr-ri*ci;
	  drs[kk]	-=ri*cr+rr*ci;
	  m++;
	  k++;
	}
      }
      m		=mp;
      kk	=m*ESNa1+i;
      k		=0;
      while(m > 0){
 	m--;
	kk	-=ESNa1;
	k++;
	rr	=Rar[k];
	ri	=-Rai[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
      }
      k		=mp;
      m		=0;
      kk	=i;
      k1	=ESFp1-k;
      if(k1 > Mr){
	k1	=Mr;
      }
      while(k1 > 0){
	rr	=Rar[k];
	ri	=Rai[k];
	drc[kk]	+=rr*cr+ri*ci;
	drs[kk]	-=ri*cr-rr*ci;
 	m++;
	kk	+=ESNa1;
	k++;
	k1--;
      }
    }
    /* Contribution of $r'_\gt\gs$ */
    cr		=gsr[0];
    kk		=i;
    for(k=1; k < Mr; k++){
      kk	+=ESNa1;
      drc[kk]	+=Rtr[k]*cr;
      drs[kk]	-=Rti[k]*cr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gsr[mp]+gsr[mm]);
      ci	=EZcr2*(gsi[mp]-gsi[mm]);
      m		=mp+1;
      kk	=m*ESNa1+i;
      k		=1;
      while(m < Mr){
	rr	=Rtr[k];
	ri	=Rti[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
	m++;
	kk	+=ESNa1;
	k++;
      }
      m		=mp;
      kk	=m*ESNa1+i;
      k		=0;
      while(m > 0){
 	m--;
	kk	-=ESNa1;
	k++;
	rr	=Rtr[k];
	ri	=-Rti[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
      }
      k		=mp;
      m		=0;
      kk	=i;
      k1	=ESFp1-k;
      if(k1 > Mr){
	k1	=Mr;
      }
      while(k1 > 0){
	rr	=Rtr[k];
	ri	=Rti[k];
	drc[kk]	+=rr*cr+ri*ci;
	drs[kk]	-=ri*cr-rr*ci;
 	m++;
	kk	+=ESNa1;
	k++;
	k1--;
      }
    }
  }

  k	=ESMp;
  gyr	=EZgper+k;
  gyi	=EZgpei+k;
  for(i=1; i < ESNa1; i++){
    in++;
    gyr		+=M0Mm;
    gyi		+=M0Mm;
    kk		=K+i;
    Rar[0]	=   rcT1a[kk];
    za		=2.*rsT1a[kk];
    for(k=1; k < ESFp1; k++){
      kk	+=ESNa1;
      Rtr[k]	=k*rsT  [kk];
      Rti[k]	=k*rcT  [kk];
      Rar[k]	=  rcT1a[kk];
      Rai[k]	= -rsT1a[kk];
    }
    rgY01a	=1./(ESsa[i]*ESgY02a[i]);
    b		=ESsb  [in];
    ba		=ESsb1a[in];
    rb		=1./b;
    gxr[0]	=0.;
    gxi[0]	=0.;
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      gxr[mp]	=-gyr[mp]*rgY01a;
      gxi[mp]	=-gyi[mp]*rgY01a;
      mm--;
      gxr[mm]	=-gyr[mm]*rgY01a;
      gxi[mm]	=-gyi[mm]*rgY01a;
    }
    /* Determination of $\gs$ */
    ur[0]	=-za*gxr[0];
    ui[0]	=-za*gxi[0];
    mp		=1;
    mm		=-1;
    while(mp < ESMp1){
      ur[mp]	=-za*gxr[mp];
      ui[mp]	=-za*gxi[mp];
      ur[mm]	=-za*gxr[mm];
      ui[mm]	=-za*gxi[mm];
      mp++;
      mm--;
    }
    ur[mp]	=0.;
    ui[mp]	=0.;
    ur[mm]	=0.;
    ui[mm]	=0.;
    k	=0;
    mp	=0;
    mm	=0;
    for(m=0; m < ESMp1; m++){
      mp++;
      mm--;
      ur[mp]	-=ba*gxi[m];
      ui[mp]	+=ba*gxr[m];
      ur[mm]	+=ba*gxi[k];
      ui[mm]	-=ba*gxr[k];
      k--;
    }
    k	=0;
    mp	=0;
    mm	=0;
    for(m=1; m < ESMp1; m++){
      k--;
      ur[mp]	+=ba*gxi[m];
      ui[mp]	-=ba*gxr[m];
      ur[mm]	-=ba*gxi[k];
      ui[mm]	+=ba*gxr[k];
      mp++;
      mm--;
    }
    mp	=ESMp;
    mm	=-ESMp;
    m	=mp+1;
    k	=mm-1;
    while(mp > 0){
      gsr[mp]	=rb*ur[m];
      gsi[mp]	=rb*ui[m];
      mp--;
      ur[mp]	-=ur[m];
      ui[mp]	-=ui[m];
      m--;
      gsr[mm]	=rb*ur[k];
      gsi[mm]	=rb*ui[k];
      mm++;
      ur[mm]	-=ur[k];
      ui[mm]	-=ui[k];
      k++;
    }
    gxr[0]	=-gyr[0]*rgY01a;
    gxi[0]	=-gyi[0]*rgY01a;
    gxr[0]	=EZcr2*(ui[-1]-ui[1])/ba;
    gxi[0]	=0.;
    dz0[i]	=za*gxr[0]-EZcr2*ur[0];
    dbT[in]	=ba*gxr[0]+EZcr2*(ui[1]-ui[-1]);
    gsr[0]	=EZcr2*rb*(ur[1]+ur[-1]);
    gsi[0]	=EZcr2*rb*(ui[1]+ui[-1]);
    /* Contribution of $r'_a\gx$ */
    kk		=i;
    cr		=gxr[0];
    drc[kk]	=Rar[0]*cr;
    drs[kk]	=0.;
    for(k=1; k < Mr; k++){
      kk	+=ESNa1;
      rr	=Rar[k];
      ri	=Rai[k];
      drc[kk]	=rr*cr;
      drs[kk]	=-ri*cr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gxr[mp]+gxr[mm]); 
      ci	=EZcr2*(gxi[mp]-gxi[mm]);
      rr	=Rar[0];
      m		=mp;
      if(m < Mr){
	kk	=m*ESNa1+i;
	drc[kk]	+=rr*cr;
	drs[kk]	-=rr*ci;
	m++;
	k	=1;
	while(m < Mr){
	  kk	+=ESNa1;
	  rr	=Rar[k];
	  ri	=Rai[k];
	  drc[kk]	+=rr*cr-ri*ci;
	  drs[kk]	-=ri*cr+rr*ci;
	  m++;
	  k++;
	}
      }
      m		=mp;
      kk	=m*ESNa1+i;
      k		=0;
      while(m > 0){
 	m--;
	kk	-=ESNa1;
	k++;
	rr	=Rar[k];
	ri	=-Rai[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
      }
      k		=mp;
      m		=0;
      kk	=i;
      k1	=ESFp1-k;
      if(k1 > Mr){
	k1	=Mr;
      }
      while(k1 > 0){
	rr	=Rar[k];
	ri	=Rai[k];
	drc[kk]	+=rr*cr+ri*ci;
	drs[kk]	-=ri*cr-rr*ci;
 	m++;
	kk	+=ESNa1;
	k++;
	k1--;
      }
    }
    /* Contribution of $r'_\gt\gs$ */
    cr		=gsr[0];
    kk		=i;
    for(k=1; k < Mr; k++){
      kk	+=ESNa1;
      drc[kk]	+=Rtr[k]*cr;
      drs[kk]	-=Rti[k]*cr;
    }
    mm		=0;
    for(mp=1; mp < ESMp1; mp++){
      mm--;
      cr	=EZcr2*(gsr[mp]+gsr[mm]);
      ci	=EZcr2*(gsi[mp]-gsi[mm]);
      m		=mp+1;
      kk	=m*ESNa1+i;
      k		=1;
      while(m < Mr){
	rr	=Rtr[k];
	ri	=Rti[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
	m++;
	kk	+=ESNa1;
	k++;
      }
      m		=mp;
      kk	=m*ESNa1+i;
      k		=0;
      while(m > 0){
 	m--;
	kk	-=ESNa1;
	k++;
	rr	=Rtr[k];
	ri	=-Rti[k];
	drc[kk]	+=rr*cr-ri*ci;
	drs[kk]	-=ri*cr+rr*ci;
      }
      k		=mp;
      m		=0;
      kk	=i;
      k1	=ESFp1-k;
      if(k1 > Mr){
	k1	=Mr;
      }
      while(k1 > 0){
	rr	=Rtr[k];
	ri	=Rti[k];
	drc[kk]	+=rr*cr+ri*ci;
	drs[kk]	-=ri*cr-rr*ci;
 	m++;
	kk	+=ESNa1;
	k++;
	k1--;
      }
    }
  }
  return(0);
}

#ifdef XWIN
int ECMonitor2D()
{
  double *ycs,*ysn;
  static int i,j,k,ki,ki0,ix,iy;
  static double s;
  extern double *EZd1drc,*EZd1drs,*EZrcs,*A2m,*EZxinf,*EZyinf;
  extern int ESMp,ESNa,ESNa1;
  char NmY[32],NmP[32];

  EZxinf[0]	=ESsa[0];
  EZxinf[1]	=ESsa[ESNa];
  EZyinf[0]	=-10.;
  EZyinf[1]	=0.;
  sprintf(NmP,"%d<= k <= %d",0,ESMp);
  SetPlotName("a","log(|dr/a|)",NmP);
  Scale2D(4,2,EZxinf,EZyinf,2,5,2);
  EZyinf[0]	=-3.;
  EZyinf[1]	=-3.;
  Plot2D(4,EZxinf,EZyinf,2,5,0,14,0);
  ycs	=EZd1drc;
  ysn	=EZd1drs;
  for(i=1; i < ESNa1; i++){
    A2m[i]= 1./EZrcs[ESNa1+i];
    A2m[i]= 1./ESsb[i];
  }
  for(k=0; k < ESMp; k++){
    j=0;
    for(i=1,ki=i+ESNa1*k; i < ESNa1; i++,ki++){
      s	=fabs(A2m[i]*ysn[ki]);
      if(s > 1.e-10){
	EZyinf[j]=log10(s);
	EZxinf[j]	=ESsa[i];
	j++;
      }
    }
    Plot2D(4,EZxinf,EZyinf,j,5,0,0,0);
    j=0;
    for(i=1,ki=i+ESNa1*k;i<ESNa1;i++,ki++){
      s= fabs(A2m[i]*ycs[ki]);
      if(s > 1.e-10){
	EZyinf[j]= log10(s);
	EZxinf[j]	=ESsa[i];
	j++;
      }
    }
    Plot2D(4,EZxinf,EZyinf,j,5,0,4,0);
  }
  return(0);
}
#endif

int ECvdrz2vrz()
{
  int i,m,k,ki;
  static double ss,s,*gx,dgx0,dgx1;
  double EZga0=0,EZga1=1e-10,EZga2=0.,EZga3=5e-11;
  EcErr	=0.;

  switch(ESEqSolvRadC){
  case 1:
    gx	=ESVs;
    dgx0=1./ESqV1a[0]-1.;
    dgx1=1./ESqV1a[ESNa]-1.;
    break;
  case 2:
    gx	=ESLs;
    dgx0=1./ESqgF1a[0]-1.;
    dgx1=1./ESqgF1a[ESNa]-1.;
    break;
  case 3:
    gx	=ESg22s;
    dgx0=1./ESqgY1a[0]-1.;
    dgx1=1./ESqgY1a[ESNa]-1.;
    break;
  default :
    gx	=ESDPs;
    dgx0=1./ESsh1a[0]-1.;
    dgx1=1./ESsh1a[ESNa]-1.;
    break;
  }
#ifdef H
  for(i=0; i < ESNa1; i++){
    gx[i]	-=dbT[i]/ESsb1a[i];
  }
  dgx1	-=(dbCc[0]-ESsb2a[ESNa]*dbT[ESNa]/ESsb1a[ESNa])/ESsb1a[ESNa];
#endif
  dR0T[0]	=drc[0];
  dZ0T[0]	=dz0[0];
  R0		=ESaR0[0]+dR0T[0];
  Z0		=ESaZ0[0]+dZ0T[0];
  for(i=0; i < ESNa1; i++){
    drcT[i]	=drc[i]-drc[0]+rcT1a[i]*gx[i];
    drsT[i]	=dz0[i]-dz0[0]+rsT1a[i]*gx[i];
    dbT[i]	=dbT[i]+ESsb1a[i]*gx[i];
  }
  drCc[0]	+=rcT1a[ESNa]*dgx1+rcT2a[ESNa]*gx[ESNa];
  drCs[0]	+=rsT1a[ESNa]*dgx1+rsT2a[ESNa]*gx[ESNa];
  dbCc[0]	+=ESsb1a[ESNa]*dgx1+ESsb2a[ESNa]*gx[ESNa];

#define Smooth00
#ifdef Smooth
  f2splA2(dbT,dbT2a,ESsa,dbCc,dbT);
  splA1a(dbT,dbT1a,dbT2a);

  f2splA(drcT,drcT2a,ESsa,drCc,drcT);
  f2splA(drsT,drsT2a,ESsa,drCs,drsT);
  splA1a(drcT,drcT1a,drcT2a);
  splA1a(drsT,drsT1a,drsT2a);
#else  

#ifdef H
  splAA2(dbT,dbT1a,dbT2a,ESsa,NULL);
  splAA(drcT,drcT1a,drcT2a,ESsa,NULL);
  splAA(drsT,drsT1a,drsT2a,ESsa,NULL);
#endif

  splAA2(dbT,dbT1a,dbT2a,ESsa,dbCc);
  splAA(drcT,drcT1a,drcT2a,ESsa,drCc);
  splAA(drsT,drsT1a,drsT2a,ESsa,drCs);

#endif
  for(i=0; i < ESNa1; i++){
    ECb[i]	=ESsb[i]+dbT[i];
    ECb1a[i]	=ESsb1a[i]+dbT1a[i];
    ECb2a[i]	=ESsb2a[i]+dbT2a[i];
    s		=fabs(dbT1a[i]/ESsb1a[i]);
    if(i){
      ss	=fabs(dbT[i]/ESsb[i]);
      if(s < ss) s	=ss;
    }
    EZrcs[i]	=rcT[i]+drcT[i];
    EZd1rcs[i]	=rcT1a[i]+drcT1a[i];
    EZd2rcs[i]	=rcT2a[i]+drcT2a[i];

#ifdef HMostOfError
    ss	=fabs(drcT1a[i]/rcT1a[i+ESNa1]);
    if(s < ss) s	=ss;
#endif
    if(i){
      ss	=fabs(drcT[i]/rcT[i+ESNa1]);
      if(s < ss) s	=ss;
    }
    EZd1drc[i]	=drcT1a[i];
    EZz0[i]	=rsT[i]+drsT[i];
    EZd1z0[i]	=rsT1a[i]+drsT1a[i];
    EZd2z0[i]	=rsT2a[i]+drsT2a[i];

#ifdef HNextMostOfError
    ss	=fabs(drsT1a[i]/ESsb1a[i]);
    if(s < ss) s	=ss;
#endif
    if(i){
      ss		=fabs(drsT[i]/ESsb[i]);
      if(s < ss) s	=ss;
    }
    drs[i]	=dz0[i];
    EZd1drs[i]	=drsT1a[i];
    if(EcErr < s){
      EcErr	=s;
      if(s >= ESEqSolvTol){
	if(i < 3){
	  Fail	|=0x00000100;
	}
	else{
	  if(i > ESNa-2){
	    Fail|=0x00001000;
	  }
	  else{
	    Fail|=0x00000010;
	  }
	}
      }
    }
  }
  ki	=ESNa1;
  for(m=1; m < Mr; m++){
    k	=ki-1;
    for(i=0; i < ESNa1; i++){
      k++;
      drcT[k]	=drc[k]+rcT1a[k]*gx[i];
      drsT[k]	=drs[k]+rsT1a[k]*gx[i];
    }
    drCc[m]	+=rcT1a[k]*dgx1+rcT2a[k]*gx[ESNa];
    drCs[m]	+=rsT1a[k]*dgx1+rsT2a[k]*gx[ESNa];
    if(m == 1){
#ifdef Smooth
      f2splA2(drcT+ki,drcT2a+ki,NULL,drCc+m,drcT+ki);
      f2splA2(drsT+ki,drsT2a+ki,NULL,drCs+m,drsT+ki);
      splA1a(drcT+ki,drcT1a+ki,drcT2a+ki);
      splA1a(drsT+ki,drsT1a+ki,drsT2a+ki);
#else
#ifdef H
      splAA2(drcT+ki,drcT1a+ki,drcT2a+ki,ESsa,NULL);
      splAA2(drsT+ki,drsT1a+ki,drsT2a+ki,ESsa,NULL);
#endif

      splAA2(drcT+ki,drcT1a+ki,drcT2a+ki,ESsa,drCc+m);
      splAA2(drsT+ki,drsT1a+ki,drsT2a+ki,ESsa,drCs+m);

#endif
    }
    else{
#ifdef H
      EZf2spl(drcT+ki,drcT2a+ki,ESsa,drCc+m,drcT+ki
	    ,NULL,ESsa,&ESNa,&EZga0,&EZga1,&EZga2,&EZga3);
      EZf2spl(drsT+ki,drsT2a+ki,ESsa,drCs+m,drsT+ki
	    ,NULL,ESsa,&ESNa,&EZga0,&EZga1,&EZga2,&EZga3);
#endif
#ifdef Smooth
      f2splA(drcT+ki,drcT2a+ki,ESsa,drCc+m,drcT+ki);
      f2splA(drsT+ki,drsT2a+ki,ESsa,drCs+m,drsT+ki);
      splA1a(drcT+ki,drcT1a+ki,drcT2a+ki);
      splA1a(drsT+ki,drsT1a+ki,drsT2a+ki);
#else
#ifdef H
      splAA(drcT+ki,drcT1a+ki,drcT2a+ki,ESsa,NULL);
      splAA(drsT+ki,drsT1a+ki,drsT2a+ki,ESsa,NULL);
#endif

      splAA(drcT+ki,drcT1a+ki,drcT2a+ki,ESsa,drCc+m);
      splAA(drsT+ki,drsT1a+ki,drsT2a+ki,ESsa,drCs+m);
#endif
    }
    for(i=0; i < ESNa1; i++){
      EZrcs[ki]	=rcT[ki]+drcT[ki];
      EZd1rcs[ki]	=rcT1a[ki]+drcT1a[ki];
      EZd2rcs[ki]	=rcT2a[ki]+drcT2a[ki];
      s		=0.;
#ifdef HNextMostOfError
      s		=fabs(drcT1a[ki]/rcT1a[ESNa1+i]);
#endif
      if(i){
	ss	=fabs(drcT[ki]/rcT[ESNa1+i]);
	if(s < ss) s	=ss;
      }
      EZd1drc[ki]	=drcT1a[ki];
      EZrsn[ki]	=rsT[ki]+drsT[ki];
      EZd1rsn[ki]	=rsT1a[ki]+drsT1a[ki];
      EZd2rsn[ki]	=rsT2a[ki]+drsT2a[ki];
#ifdef HNextMostOfError
      ss	=fabs(drsT1a[ki]/rcT1a[ESNa1+i]);
      if(s < ss) s	=ss;
#endif
      if(i){
	ss	=fabs(drsT[ki]/rcT[ESNa1+i]);
	if(s < ss) s	=ss;
      }
      EZd1drs[ki]	=drsT1a[ki];
      if(EcErr < s){
	EcErr	=s;
	if(s >= ESEqSolvTol){
	  if(i < 3){
	    Fail	|=0x00000100;
	  }
	  else{
	    if(i > ESNa-2){
	      Fail|=0x00001000;
	    }
	    else{
	      Fail|=0x00000010;
	    }
	  }
	}
      }
      ki++;
    }
  }
  return(0);
}

int ECvdrz2vrzHole1()
{
  int i,i1,m,k,km,ki;
  static double ss,s,*gx,dgx0,dgx1;
  
  EcErr	=0.;
  switch(ESEqSolvRadC){
  case 1:
    gx	=ESVs;
    dgx0=1./ESqV1a[0]-1.;
    dgx1=1./ESqV1a[ESNa]-1.;
    break;
  case 2:
    gx	=ESLs;
    dgx0=1./ESqgF1a[0]-1.;
    dgx1=1./ESqgF1a[ESNa]-1.;
    break;
  case 3:
    gx	=ESg22s;
    dgx0=1./ESqgY1a[0]-1.;
    dgx1=1./ESqgY1a[ESNa]-1.;
    break;
  default :
    gx	=ESDPs;
    dgx0=1./ESsh1a[0]-1.;
    dgx1=1./ESsh1a[ESNa]-1.;
    break;
  }

  dR0T[0]	=drc[0];
  dZ0T[0]	=dz0[0];
  R0		=ESaR0[0]+dR0T[0];
  Z0		=ESaZ0[0]+dZ0T[0];
  ki		=0;
  k		=0;
  for(i=0; i < ESNa1; i++){
    drcT[k]	=drc[i]-drc[0]+rcT1a[k]*gx[i];
    drsT[k]	=dz0[i]-dz0[0]+rsT1a[k]*gx[i];
    dbT[i]	=dbT[i]+ESsb1a[i]*gx[i];
    k++;
  }
  k--;
  drCc[0]	+=rcT1a[k]*dgx1+rcT2a[k]*gx[ESNa];
  drCs[0]	+=rsT1a[k]*dgx1+rsT2a[k]*gx[ESNa];
  i--;
  dbT1a[i]	=ESsb1a[ESNa]*dgx1+ESsb2a[ESNa]*gx[ESNa];
  splAA(dbT,dbT1a,dbT2a,NULL,dbT1a+ESNa);
#define SSmooth
#ifdef SSmooth
  f2splA(drcT+ki,drcT2a+ki,NULL,drCc,drcT+ki);
  f2splA(drsT+ki,drsT2a+ki,NULL,drCs,drsT+ki);
  splA1a(drcT+ki,drcT1a+ki,drcT2a+ki);
  splA1a(drsT+ki,drsT1a+ki,drsT2a+ki);
#else  
  splAA(drcT+ki,drcT1a+ki,drcT2a+ki,NULL,drCc);
  splAA(drsT+ki,drsT1a+ki,drsT2a+ki,NULL,drCs);
#endif
  for(i=0; i < ESNa1; i++){
    ECb[i]	=ESsb[i]+dbT[i];
    ECb1a[i]	=ESsb1a[i]+dbT1a[i];
    ECb2a[i]	=ESsb2a[i]+dbT2a[i];
    s		=fabs(dbT1a[i]/ESsb1a[i]);
    if(i){
      ss		=fabs(dbT[i]/ESsb[i]);
      if(s < ss){
	s	=ss;
      }
    }
    EZrcs[i]	=rcT[ki]+drcT[ki];
    EZd1rcs[i]	=rcT1a[ki]+drcT1a[ki];
    EZd2rcs[i]	=rcT2a[ki]+drcT2a[ki];
#ifdef HMostOfError
    ss		=fabs(drcT1a[ki]/rcT1a[ki+ESNa1]);
    if(s < ss){
      s	=ss;
    }
#endif
    if(i){
      ss		=fabs(drcT[ki]/rcT[ki+ESNa1]);
      if(s < ss){
	s	=ss;
      }
    }
    EZd1drc[i]	=drcT1a[ki];

    EZz0[i]	=rsT[ki]+drsT[ki];
    EZd1z0[i]	=rsT1a[ki]+drsT1a[ki];
    EZd2z0[i]	=rsT2a[ki]+drsT2a[ki];
#ifdef HNextMostOfError
    ss		=fabs(drsT1a[ki]/ESsb1a[i]);
    if(s < ss){
      s	=ss;
    }
#endif
    if(i){
      ss		=fabs(drsT[ki]/ESsb[i]);
      if(s < ss){
	s	=ss;
      }
    }
    drs[i]	=dz0[i];
    EZd1drs[i]	=drsT1a[ki];
    if(EcErr < s){
      EcErr	=s;
      if(s >= ESEqSolvTol){
	if(i < 3){
	  Fail	|=0x00000100;
	}
	else{
	  if(i > ESNa-2){
	    Fail|=0x00001000;
	  }
	  else{
	    Fail|=0x00000010;
	  }
	}
      }
    }
    ki++;
  }
  km	=ESNa1;
  for(m=1; m < Mr; m++){
    k	=ki-1;
    for(i=0; i < ESNa1; i++){
      k++;
      drcT[k]	=drc[km+i]+rcT1a[k]*gx[i];
      drsT[k]	=drs[km+i]+rsT1a[k]*gx[i];
    }
    drCc[m]	+=rcT1a[k]*dgx1+rcT2a[k]*gx[ESNa];
    drCs[m]	+=rsT1a[k]*dgx1+rsT2a[k]*gx[ESNa];
#ifdef SSmooth
    f2splA(drcT+ki,drcT2a+ki,NULL,drCc+m,drcT+ki);
    f2splA(drsT+ki,drsT2a+ki,NULL,drCs+m,drsT+ki);
    splA1a(drcT+ki,drcT1a+ki,drcT2a+ki);
    splA1a(drsT+ki,drsT1a+ki,drsT2a+ki);
#else
    splAA(drcT+ki,drcT1a+ki,drcT2a+ki,NULL,drCc+m);
    splAA(drsT+ki,drsT1a+ki,drsT2a+ki,NULL,drCs+m);
#endif
    for(i=0; i < ESNa1; i++){
      EZrcs[km]	=rcT[ki]+drcT[ki];
      EZd1rcs[km]	=rcT1a[ki]+drcT1a[ki];
      EZd2rcs[km]	=rcT2a[ki]+drcT2a[ki];
      s		=0.;
#ifdef HNextMostOfError
      s		=fabs(drcT1a[ki]/rcT1a[ESNa1+i]);
#endif
      if(i){
	ss	=fabs(drcT[ki]/rcT[ESNa1+i]);
	if(s < ss){
	  s	=ss;
	}
      }
      EZd1drc[km]	=drcT1a[ki];
      EZrsn[km]	=rsT[ki]+drsT[ki];
      EZd1rsn[km]	=rsT1a[ki]+drsT1a[ki];
      EZd2rsn[km]	=rsT2a[ki]+drsT2a[ki];
#ifdef HNextMostOfError
      ss	=fabs(drsT1a[ki]/rcT1a[ESNa1+i]);
      if(s < ss){
	s	=ss;
      }
#endif
      if(i){
	ss	=fabs(drsT[ki]/rcT[ESNa1+i]);
	if(s < ss){
	  s	=ss;
	}
      }
      EZd1drs[km]	=drsT1a[ki];
      if(EcErr < s){
	EcErr	=s;
	if(s >= ESEqSolvTol){
	  if(i < 3){
	    Fail	|=0x00000100;
	  }
	  else{
	    if(i > ESNa-2){
	      Fail|=0x00001000;
	    }
	    else{
	      Fail|=0x00000010;
	    }
	  }
	}
      }
      ki++;
      km++;
    }
  }
#ifdef DEBUG
  {
    int n,k0;
    double x[2],y[2],*Wx,**W1,Y[400],Y1[400],gx0[400],gx1[400];
    double *p,*p1,*p2,s,t,s1,sr,si,*pr,*pi;

    ki	=ESNa*8+1;
    Wx	=(double*)malloc(ki*sizeof(double));
    if(Wx == NULL){
      printf("Failure in memory allocation for Wx in EC2DEqSolv()\n");
      exit(0);
    }
    W1	=(double**)malloc(Mr*sizeof(double*));
    if(W1 == NULL){
      printf("Failure in memory allocation for W1 in EC2DEqSolv()\n");
      exit(0);
    }
    for(k=0; k < Mr; k++){
      W1[k]	=(double*)malloc(ki*sizeof(double));
      if(W1[k] == NULL){
	printf("Failure in memory allocation for W1[%2d] in EC2DEqSolv()\n",k);
	exit(0);
      }
    }

    p		=drsT;
    p1		=drcT1a;
    p2		=drsT2a;
    k0		=0;
    y[0]	=p[0];
    y[1]	=y[0];
    x[0]	=0.;
    x[1]	=1.;
    ki		=0;
    s		=0.125*ESsa[1];
    for(k=k0; k < Mr; k++){
      ki	=ESNa1*k;
      n	=0;
      t		=0.;
      while(t <= 1.){
	ESSetSplA(t);
	Wx[n]	=t;
	splRA(W1[k]+n,&s1,p+ki,p2+ki);
	t	+=s;
	if(y[0] > W1[k][n])
	  y[0]	=W1[k][n];
	if(y[1] < W1[k][n])
	  y[1]	=W1[k][n];
	n++;
      }
    }
    Scale2D(3,2,x,y,2,6,2);
    for(k=k0; k < Mr; k++){
      ki	=ESNa1*k;
      Plot2d(3,ESsa,p+ki,ESNa1,6,0,4,0);
      if(k){
	Plot2d(3,Wx,W1[k],n,6,0,0,0);
      }
      else{
	Plot2d(3,Wx,W1[k],n,6,0,14,0);
      }
    }

    y[0]	=ESgY02a[0];
    y[1]	=y[0];
    for(k=0; k < ESNa1; k++){
      Y[k]	=ESgY0[k];
      if(y[0] > Y[k])
	y[0]	=Y[k];
      if(y[1] < Y[k])
	y[1]	=Y[k];
      Y1[k]	=EZgper[ESMp+k*M0Mm];
      if(y[0] > Y1[k])
	y[0]	=Y1[k];
      if(y[1] < Y1[k])
	y[1]	=Y1[k];
    }
    CbFlush();
    for(k=0; k < Mr; k++){
      free(W1[k]);
    }
    free(W1);
    free(Wx);
  }
  {
    double *ycs,*ysn;
    static int i,j,k,ki,ki0,ix,iy;
    static double s;
    extern double *EZd1drc,*EZd1drs,*EZrcs,*A2m,*EZxinf,*EZyinf;
    extern int ESMp,ESNa,ESNa1;

    EZxinf[0]	=ESsa[0];
    EZxinf[1]	=ESsa[ESNa];
    EZyinf[0]	=-10.;
    EZyinf[1]	=0.;
    Scale2D(4,2,EZxinf,EZyinf,2,5,2);
    EZyinf[0]	=-3.;
    EZyinf[1]	=-3.;
    Plot2D(4,EZxinf,EZyinf,2,5,0,14,0);
    ycs	=EZd1drc;
    ysn	=EZd1drs;
    for(i=1; i < ESNa1; i++){
      A2m[i]= 1./EZrcs[ESNa1+i];
    }
    for(k=0; k < ESMp; k++){
      j=0;
      for(i=1,ki= i+ESNa1*k;i<ESNa1;i++,ki++){
	s= fabs(A2m[i]*ysn[ki]);
	if(s > 1.e-10){
	  EZyinf[j]=log10(s);
	  EZxinf[j]	=ESsa[i];
	  j++;
	}
      }
      Plot2D(4,EZxinf,EZyinf,j,5,0,0,0);
      j=0;
      for(i=1,ki=i+ESNa1*k;i<ESNa1;i++,ki++){
	s= fabs(A2m[i]*ycs[ki]);
	if(s > 1.e-10){
	  EZyinf[j]= log10(s);
	  EZxinf[j]	=ESsa[i];
	  j++;
	}
      }
      Plot2D(4,EZxinf,EZyinf,j,5,0,4,0);
    }
    CbFlush();
  }
#endif
  return(0);
}

int ECvdrz2vrzHole(int n)
{
  int i,i1,m,k,km,in,ki,K,N;
  static double ss,s,*gx,dgx0,dgx1;
  
  EcErr	=0.;
  N	=ESNa1*n;
  K	=ESnAF*n;
  switch(ESEqSolvRadC){
  case 1:
    gx	=ESVs;
    dgx0=1./ESqV1a[0]-1.;
    dgx1=1./ESqV1a[ESNa]-1.;
    break;
  case 2:
    gx	=ESLs;
    dgx0=1./ESqgF1a[0]-1.;
    dgx1=1./ESqgF1a[ESNa]-1.;
    break;
  case 3:
    gx	=ESg22s;
    dgx0=1./ESqgY1a[0]-1.;
    dgx1=1./ESqgY1a[ESNa]-1.;
    break;
  default :
    gx	=ESDPs;
    dgx0=1./ESsh1a[0]-1.;
    dgx1=1./ESsh1a[ESNa]-1.;
    break;
  }

  dR0T[n]	=drc[0];
  dZ0T[n]	=dz0[0];
  R0		=ESaR0[n]+dR0T[n];
  Z0		=ESaZ0[n]+dZ0T[n];
  ki		=K;
  k		=ki;
  in		=N;
  for(i=0; i < ESNa1; i++){
    drcT[k]	=drc[i]-drc[0]+rcT1a[k]*gx[i];
    drsT[k]	=dz0[i]-dz0[0]+rsT1a[k]*gx[i];
    dbT[in]	=dbT[in]+ESsb1a[in]*gx[i];
    k++;
    in++;
  }
  k--;
  drCc[0]	+=rcT1a[k]*dgx1+rcT2a[k]*gx[ESNa];
  drCs[0]	+=rsT1a[k]*dgx1+rsT2a[k]*gx[ESNa];
  in--;
  dbT1a[in]	=ESsb1a[ESNa]*dgx1+ESsb2a[ESNa]*gx[ESNa];
  in		=N;
  splAA(dbT,dbT1a,dbT2a,NULL,dbT1a+ESNa);
#define SSmooth
#ifdef SSmooth
  f2splA(drcT+ki,drcT2a+ki,NULL,drCc,drcT+ki);
  f2splA(drsT+ki,drsT2a+ki,NULL,drCs,drsT+ki);
  splA1a(drcT+ki,drcT1a+ki,drcT2a+ki);
  splA1a(drsT+ki,drsT1a+ki,drsT2a+ki);
#else  
  splAA(drcT+ki,drcT1a+ki,drcT2a+ki,NULL,drCc);
  splAA(drsT+ki,drsT1a+ki,drsT2a+ki,NULL,drCs);
#endif
  for(i=0; i < ESNa1; i++){
    ECb[i]	=ESsb[in]+dbT[in];
    ECb1a[i]	=ESsb1a[in]+dbT1a[in];
    ECb2a[i]	=ESsb2a[in]+dbT2a[in];
    s		=fabs(dbT1a[in]/ESsb1a[in]);
    if(i){
      ss		=fabs(dbT[in]/ESsb[in]);
      if(s < ss){
	s	=ss;
      }
    }
    EZrcs[i]	=rcT[ki]+drcT[ki];
    EZd1rcs[i]	=rcT1a[ki]+drcT1a[ki];
    EZd2rcs[i]	=rcT2a[ki]+drcT2a[ki];
#ifdef HMostOfError
    ss		=fabs(drcT1a[ki]/rcT1a[ki+ESNa1]);
    if(s < ss){
      s	=ss;
    }
#endif
    if(i){
      ss		=fabs(drcT[ki]/rcT[ki+ESNa1]);
      if(s < ss){
	s	=ss;
      }
    }
    EZd1drc[i]	=drcT1a[ki];

#ifdef H
/**????????****/
    drsT[ki]	=0.;
    drsT1a[ki]	=0.;
    drsT2a[ki]	=0.;
    dz0[i]	=0.;
/**????????****/
#endif
    EZz0[i]	=rsT[ki]+drsT[ki];
    EZd1z0[i]	=rsT1a[ki]+drsT1a[ki];
    EZd2z0[i]	=rsT2a[ki]+drsT2a[ki];
#ifdef HNextMostOfError
    ss		=fabs(drsT1a[ki]/ESsb1a[in]);
    if(s < ss){
      s	=ss;
    }
#endif
    if(i){
      ss		=fabs(drsT[ki]/ESsb[in]);
      if(s < ss){
	s	=ss;
      }
    }
    drs[i]	=dz0[i];
    EZd1drs[i]	=drsT1a[ki];
    if(EcErr < s){
      EcErr	=s;
      if(s >= ESEqSolvTol){
	if(i < 3){
	  Fail	|=0x00000100;
	}
	else{
	  if(i > ESNa-2){
	    Fail|=0x00001000;
	  }
	  else{
	    Fail|=0x00000010;
	  }
	}
      }
    }
    in++;
    ki++;
  }
  km	=ESNa1;
  for(m=1; m < Mr; m++){
    k	=ki-1;
    for(i=0; i < ESNa1; i++){
      k++;
      drcT[k]	=drc[km+i]+rcT1a[k]*gx[i];
      drsT[k]	=drs[km+i]+rsT1a[k]*gx[i];
    }
    drCc[m]	+=rcT1a[k]*dgx1+rcT2a[k]*gx[ESNa];
    drCs[m]	+=rsT1a[k]*dgx1+rsT2a[k]*gx[ESNa];
    if(m == 1){
#ifdef SSmooth
      f2splA(drcT+ki,drcT2a+ki,NULL,drCc+m,drcT+ki);
      f2splA(drsT+ki,drsT2a+ki,NULL,drCs+m,drsT+ki);
      splA1a(drcT+ki,drcT1a+ki,drcT2a+ki);
      splA1a(drsT+ki,drsT1a+ki,drsT2a+ki);
#else
      splAA(drcT+ki,drcT1a+ki,drcT2a+ki,NULL,drCc+m);
      splAA(drsT+ki,drsT1a+ki,drsT2a+ki,NULL,drCs+m);
#endif
    }
    else{
#ifdef SSmooth
      f2splA(drcT+ki,drcT2a+ki,NULL,drCc+m,drcT+ki);
      f2splA(drsT+ki,drsT2a+ki,NULL,drCs+m,drsT+ki);
      splA1a(drcT+ki,drcT1a+ki,drcT2a+ki);
      splA1a(drsT+ki,drsT1a+ki,drsT2a+ki);
#else
      splAA(drcT+ki,drcT1a+ki,drcT2a+ki,NULL,drCc+m);
      splAA(drsT+ki,drsT1a+ki,drsT2a+ki,NULL,drCs+m);
#endif
    }
    for(i=0; i < ESNa1; i++){
      EZrcs[km]	=rcT[ki]+drcT[ki];
      EZd1rcs[km]	=rcT1a[ki]+drcT1a[ki];
      EZd2rcs[km]	=rcT2a[ki]+drcT2a[ki];
      s		=0.;
#ifdef HNextMostOfError
      s		=fabs(drcT1a[ki]/rcT1a[K+ESNa1+i]);
#endif
      if(i){
	ss	=fabs(drcT[ki]/rcT[K+ESNa1+i]);
	if(s < ss){
	  s	=ss;
	}
      }
      EZd1drc[km]	=drcT1a[ki];
      EZrsn[km]	=rsT[ki]+drsT[ki];
      EZd1rsn[km]	=rsT1a[ki]+drsT1a[ki];
      EZd2rsn[km]	=rsT2a[ki]+drsT2a[ki];
#ifdef HNextMostOfError
      ss	=fabs(drsT1a[ki]/rcT1a[K+ESNa1+i]);
      if(s < ss){
	s	=ss;
      }
#endif
      if(i){
	ss	=fabs(drsT[ki]/rcT[K+ESNa1+i]);
	if(s < ss){
	  s	=ss;
	}
      }
      EZd1drs[km]	=drsT1a[ki];
      if(EcErr < s){
	EcErr	=s;
	if(s >= ESEqSolvTol){
	  if(i < 3){
	    Fail	|=0x00000100;
	  }
	  else{
	    if(i > ESNa-2){
	      Fail|=0x00001000;
	    }
	    else{
	      Fail|=0x00000010;
	    }
	  }
	}
      }
      ki++;
      km++;
    }
  }
#ifdef DEBUG
  {
    int n,k0;
    double x[2],y[2],*Wx,**W1,Y[400],Y1[400],gx0[400],gx1[400];
    double *p,*p1,*p2,s,t,s1,sr,si,*pr,*pi;

    ki	=ESNa*8+1;
    Wx	=(double*)malloc(ki*sizeof(double));
    if(Wx == NULL){
      printf("Failure in memory allocation for Wx in EC2DEqSolv()\n");
      exit(0);
    }
    W1	=(double**)malloc(Mr*sizeof(double*));
    if(W1 == NULL){
      printf("Failure in memory allocation for W1 in EC2DEqSolv()\n");
      exit(0);
    }
    for(k=0; k < Mr; k++){
      W1[k]	=(double*)malloc(ki*sizeof(double));
      if(W1[k] == NULL){
	printf("Failure in memory allocation for W1[%2d] in EC2DEqSolv()\n",k);
	exit(0);
      }
    }

    p		=drsT;
    p1		=drcT1a;
    p2		=drsT2a;
    k0		=0;
    y[0]	=p[0];
    y[1]	=y[0];
    x[0]	=0.;
    x[1]	=1.;
    ki		=0;
    s		=0.125*ESsa[1];
    for(k=k0; k < Mr; k++){
      ki	=ESNa1*k;
      n	=0;
      t		=0.;
      while(t <= 1.){
	ESSetSplA(t);
	Wx[n]	=t;
	splRA(W1[k]+n,&s1,p+ki,p2+ki);
	t	+=s;
	if(y[0] > W1[k][n]) y[0]=W1[k][n];
	if(y[1] < W1[k][n]) y[1]=W1[k][n];
	n++;
      }
    }
    Scale2D(3,2,x,y,2,6,2);
    for(k=k0; k < Mr; k++){
      ki	=ESNa1*k;
      Plot2d(3,ESsa,p+ki,ESNa1,6,0,4,0);
      if(k){
	Plot2d(3,Wx,W1[k],n,6,0,0,0);
      }
      else{
	Plot2d(3,Wx,W1[k],n,6,0,14,0);
      }
    }
    y[0]	=ESgY02a[0];
    y[1]	=y[0];
    for(k=0; k < ESNa1; k++){
      Y[k]	=ESgY0[k];
      if(y[0] > Y[k])
	y[0]	=Y[k];
      if(y[1] < Y[k])
	y[1]	=Y[k];
      Y1[k]	=EZgper[ESMp+k*M0Mm];
      if(y[0] > Y1[k])
	y[0]	=Y1[k];
      if(y[1] < Y1[k])
	y[1]	=Y1[k];
    }
    CbFlush();
    for(k=0; k < Mr; k++){
      free(W1[k]);
    }
    free(W1);
    free(Wx);
  }
  {
    double *ycs,*ysn;
    static int i,j,k,ki,ki0,ix,iy;
    static double s;
    extern double *EZd1drc,*EZd1drs,*EZrcs,*A2m,*EZxinf,*EZyinf;
    extern int ESMp,ESNa,ESNa1;

    EZxinf[0]	=ESsa[0];
    EZxinf[1]	=ESsa[ESNa];
    EZyinf[0]	=-10.;
    EZyinf[1]	=0.;
    Scale2D(4,2,EZxinf,EZyinf,2,5,2);
    EZyinf[0]	=-3.;
    EZyinf[1]	=-3.;
    Plot2D(4,EZxinf,EZyinf,2,5,0,14,0);
    ycs	=EZd1drc;
    ysn	=EZd1drs;
    for(i=1; i < ESNa1; i++){
      A2m[i]= 1./EZrcs[ESNa1+i];
    }
    for(k=0; k < ESMp; k++){
      j=0;
      for(i=1,ki= i+ESNa1*k;i<ESNa1;i++,ki++){
	s= fabs(A2m[i]*ysn[ki]);
	if(s > 1.e-10){
	  EZyinf[j]=log10(s);
	  EZxinf[j]	=ESsa[i];
	  j++;
	}
      }
      Plot2D(4,EZxinf,EZyinf,j,5,0,0,0);
      j=0;
      for(i=1,ki=i+ESNa1*k;i<ESNa1;i++,ki++){
	s= fabs(A2m[i]*ycs[ki]);
	if(s > 1.e-10){
	  EZyinf[j]= log10(s);
	  EZxinf[j]	=ESsa[i];
	  j++;
	}
      }
      Plot2D(4,EZxinf,EZyinf,j,5,0,4,0);
    }
    CbFlush();
  }
#endif
  return(0);
}

int ECSmooth2DCoord(int n)
{
  int i,i1,m,k,km,in,ki,K;

  in	=ESNa1*n;
  K	=ESnAF*n;
  ki	=K;
  km	=ki+ESNa;

  if(ESa0 != 0.){
    f2splA(ESsb,ESsb2a,ESsb1a,ESsb1a+ESNa,ESsb);
    f2splA(rcT+ki,rcT2a+ki,rcT1a+ki,rcT1a+km,rcT+ki);
    f2splA(rsT+ki,rsT2a+ki,rsT1a+ki,rsT1a+km,rsT+ki);
  }
  else{
    f2splA2(ESsb,ESsb2a,ESsa,ESsb1a+ESNa,ESsb);
    f2splA(rcT+ki,rcT2a+ki,ESsa,rcT1a+km,rcT+ki);
    f2splA(rsT+ki,rsT2a+ki,ESsa,rsT1a+km,rsT+ki);
  }
  splA1a(ESsb,ESsb1a,ESsb2a);
  splA1a(rcT+ki,rcT1a+ki,rcT2a+ki);
  splA1a(rsT+ki,rsT1a+ki,rsT2a+ki);

  for(m=1; m < ESFp1; m++){
    ki	+=ESNa1;
    km	+=ESNa1;

    if(ESa0 != 0.){
      f2splA(rcT+ki,rcT2a+ki,rcT1a+ki,rcT1a+km,rcT+ki);
      f2splA(rsT+ki,rsT2a+ki,rsT1a+ki,rsT1a+km,rsT+ki);
    }
    else{
      if(m == 1){
	f2splA2(rcT+ki,rcT2a+ki,ESsa,rcT1a+km,rcT+ki);
	f2splA2(rsT+ki,rsT2a+ki,ESsa,rsT1a+km,rsT+ki);
      }
      else{
	f2splA(rcT+ki,rcT2a+ki,ESsa,rcT1a+km,rcT+ki);
	f2splA(rsT+ki,rsT2a+ki,ESsa,rsT1a+km,rsT+ki);
      }
    }
    splA1a(rcT+ki,rcT1a+ki,rcT2a+ki);
    splA1a(rsT+ki,rsT1a+ki,rsT2a+ki);
  }
#ifdef DEBUG
  {
    int n;
    double x[2],y[2],*Wx,**W1;
    double *p,*p1,*p2,s,t;

    ki	=ESNa*8+1;
    Wx	=(double*)malloc(ki*sizeof(double));
    if(Wx == NULL){
      printf("Failure in memory allocation for Wx in EC2DEqSolv()\n");
      exit(0);
    }
    W1	=(double**)malloc(ESFp1*sizeof(double*));
    if(W1 == NULL){
      printf("Failure in memory allocation for W1 in EC2DEqSolv()\n");
      exit(0);
    }
    for(k=0; k < ESFp1; k++){
      W1[k]	=(double*)malloc(ki*sizeof(double));
      if(W1[k] == NULL){
	printf("Failure in memory allocation for W1[%2d] in EC2DEqSolv()\n",k);
	exit(0);
      }
    }

    p		=rcT;
    p1		=rcT1a;
    p2		=rcT2a;
    y[0]	=p[0];
    y[1]	=y[0];
    x[0]	=0.;
    x[1]	=1.;
    ki		=0;
    s		=0.125*ESsa[1];
    for(k=0; k < ESFp1; k++){
      ki	=ESNa1*k;
      n	=0;
      t		=0.;
      while(t <= 1.){
	ESSetSplA(t);
	Wx[n]	=t;
	splRA(W1[k]+n,NULL,p+ki,p2+ki);
	if(y[0] > W1[k][n])
	  y[0]	=W1[k][n];
	if(y[1] < W1[k][n])
	  y[1]	=W1[k][n];
	t	+=s;
	n++;
      }
    }
    for(k=0; k < ESFp1; k++){
      ki	=ESNa1*k;
      Plot2d(3,ESsa,p+ki,ESNa1,6,0,4,0);
      Plot2d(3,Wx,W1[k],n,6,0,0,0);
    }
    CbFlush();
    for(k=0; k < ESFp1; k++){
      free(W1[k]);
    }
    free(W1);
    free(Wx);
  }
#endif
  return(0);
}

static double Sc=1.;

int ECvdrzCorrection(int n)
{
  int i,m,k,km,ki;
  double s;

  dR0T[n]	*=Sc;
  dZ0T[n]	*=Sc;
  R0		=ESaR0[n]+dR0T[n];
  Z0		=ESaZ0[n]+dZ0T[n];
  ki		=0;
  for(i=0; i < ESNa1; i++){
    dbT[i]	*=Sc;
    dbT1a[i]	*=Sc;
    dbT2a[i]	*=Sc;
    drcT[ki]	*=Sc;
    drcT1a[ki]	*=Sc;
    drcT2a[ki]	*=Sc;
    drsT[ki]	*=Sc;
    drsT1a[ki]	*=Sc;
    drsT2a[ki]	*=Sc;
    ECb[i]	=ESsb[i]	+dbT[i];
    ECb1a[i]	=ESsb1a[i]	+dbT1a[i];
    ECb2a[i]	=ESsb2a[i]	+dbT2a[i];
    EZrcs[i]	=rcT[ki]	+drcT[ki];
    EZd1rcs[i]	=rcT1a[ki]	+drcT1a[ki];
    EZd2rcs[i]	=rcT2a[ki]	+drcT2a[ki];
    EZz0[i]	=rsT[ki]	+drsT[ki];
    EZd1z0[i]	=rsT1a[ki]	+drsT1a[ki];
    EZd2z0[i]	=rsT2a[ki]	+drsT2a[ki];
    ki++;
  }
  km		=ESNa1;
  for(m=1; m < Mr; m++){
    for(i=0; i < ESNa1; i++){
      drcT  [ki]*=Sc;
      drcT1a[ki]*=Sc;
      drcT2a[ki]*=Sc;
      drsT  [ki]*=Sc;
      drsT1a[ki]*=Sc;
      drsT2a[ki]*=Sc;
      EZrcs  [km]	=rcT  [ki]	+drcT  [ki];
      EZd1rcs[km]	=rcT1a[ki]	+drcT1a[ki];
      EZd2rcs[km]	=rcT2a[ki]	+drcT2a[ki];
      EZrsn  [km]	=rsT  [ki]	+drsT  [ki];
      EZd1rsn[km]	=rsT1a[ki]	+drsT1a[ki];
      EZd2rsn[km]	=rsT2a[ki]	+drsT2a[ki];
      ki++;
      km++;
    }
  }
  if(ESiAx != 0){
    ECvdrzCorrectionX();
  }
  return(0);
}

int ECCheckJacob()
{
  int i,j,ji;
  int ki,k,kj;
  double rc[ESFp1],rs[ESFp1],rca[ESFp1],rsa[ESFp1];
  double D,D1,D2,bb,ba,ra,rt,za,zt; 

  Sc	=1.;
  D	=rT1a[0]*zT1gt[0]-zT1a[0]*rT1gt[0];
  D1	=2.*(drcT1a[ESNa1]*ESsb1a[0]+rcT1a[ESNa1]*dbT1a[0]);
  D2	=2.*(drcT1a[ESNa1]*dbT1a[0]);

  if(0.55*D+D1+D2 < 0.){
    za	=sqrt(D1*D1-2.*D*D2);
    if(D1 > 0.){
      ra	=D/(za-D1);
    }
    else{
      ra	=D/(za-D1);
      if(za+D1 < 0.){
	za	=-D/(za+D1);
	if(ra > za){
	  ra	=za;
	}
      }
    }
    if(Sc > ra){
      Sc	=ra;
      Fail	|=0x00010000;
    }
  }

  ji	=ESNp1;
  for(i=1; i < ESNa1; i++){
    bb		=dbT[i];
    ba		=dbT1a[i];
    ki		=i;
    rca[0]	=drcT1a[ki];
    rsa[0]	=drsT1a[ki];
    for(k=1; k < Mr; k++){
      ki	+=ESNa1;
      rc[k]	=2.*k*drcT[ki];
      rs[k]	=2.*k*drsT[ki];
      rca[k]	=2.*drcT1a[ki];
      rsa[k]	=2.*drsT1a[ki];
    }
    for(j=0; j < ESNp; j++){
      rt	=0.;
      ra	=rca[0];
      kj	=0;
      for(k=1; k < Mr; k++){
	kj	+=j;
	if(kj >= ESNp) kj	-=ESNp;
	rt	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
	ra	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
      }
      zt	=bb*EScs1[j];
      za	=rsa[0]+ba*ESsn1[j];
      D		=rT1a[ji]*zT1gt[ji]-rT1gt[ji]*zT1a[ji];
      D1	=rT1a[ji]*zt+ra*zT1gt[ji]-rT1gt[ji]*za-rt*zT1a[ji];
      D2	=ra*zt-rt*za;
      if(D+D1+D2 < 0.){
	za	=sqrt(D1*D1-2.*D*D2);
	ra	=D1 > 0. ? -D/(za+D1) : D/(za-D1);
	if(Sc > ra){
	  Sc	=ra;
	  if(i < 3){
	    Fail	|=0x00010000;
	  }
	  else{
	    if(i > ESNa-2){
	      Fail	|=0x00100000;
	    }
	    else{
	      Fail	|=0x01000000;
	    }
	  }
	}
#ifdef DEBUG
	printf("D[%2d][%2d]=%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n"
	       ,i,j,D,D+D1+D2,D+Sc*(D1+Sc*D2),D+ra*(D1+ra*D2),ra,Sc);
#endif
      }
      ji++;
    }
    ji++;
  }
  
  if(Sc < 1.){
#ifdef DEBUG
    printf("Sc=%10.3e kfail=%d\a\n",Sc,kFail);
#endif
    if(kFail){
      ECRestore2DGeom();
    }
    else{
      ECvdrzCorrection(0);
    }
    if(Sc < 0.1){
      kFail++;
      return(1);
    }
  }
  kFail	=0;
  ESFail&=0x0000FFFF;
  return(0);
}

int ECCheckJacobHole1()
{
  int i,j,ji;
  int ki,k,kj;
  double rc[ESFp1],rs[ESFp1],rca[ESFp1],rsa[ESFp1];
  double *db,*db1a,*drc,*drca,*drs,*drsa;
  double D,D1,D2,bb,ba,ra,rgt,za,zgt;

  ki	=0;
  ji	=0;

  db	=dbT;
  db1a	=dbT1a;

  drc	=drcT;
  drca	=drcT1a;
  drs	=drsT;
  drsa	=drsT1a;

  Sc	=1.;
  for(i=0; i < ESNa1; i++){
    bb		=db[i];
    ba		=db1a[i];
    ki		=i;
    rca[0]	=drca[ki];
    rsa[0]	=drsa[ki];
    for(k=1; k < Mr; k++){
      ki	+=ESNa1;
      rc[k]	=2.*k*drc[ki];
      rs[k]	=2.*k*drs[ki];
      rca[k]	=2.*drca[ki];
      rsa[k]	=2.*drsa[ki];
    }
    for(j=0; j < ESNp; j++){
      rgt	=0.;
      ra	=rca[0];
      kj	=0;
      for(k=1; k < Mr; k++){
	kj	+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	rgt	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
	ra	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
      }
      zgt	=bb*EScs1[j];
      za	=rsa[0]+ba*ESsn1[j];
      D		=rT1a[ji]*zT1gt[ji]-rT1gt[ji]*zT1a[ji];
      D1	=rT1a[ji]*zgt+ra*zT1gt[ji]-rT1gt[ji]*za-rgt*zT1a[ji];
      D2	=ra*zgt-rgt*za;
      if(0.55*D+D1+D2 < 0.){
	za	=sqrt(D1*D1-2.*D*D2);
	if(D1 > 0.){
	  ra	=D/(za-D1);
	}
	else{
	  ra	=D/(za-D1);
	  if(za+D1 < 0.){
	    za	=-D/(za+D1);
	    if(za > 0. && ra > za){
	      ra	=za;
	    }
	  }
	}
	if(Sc > ra){
	  Sc	=ra;
	  if(i < 3){
	    Fail	|=0x00010000;
	  }
	  else{
	    if(i > ESNa-2){
	      Fail	|=0x00100000;
	    }
	    else{
	      Fail	|=0x01000000;
	    }
	  }
	  
#ifdef DEBUG
	  printf("D[%2d][%2d]=%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n"
		 ,i,j,D,D+D1+D2,D+Sc*(D1+Sc*D2),D+ra*(D1+ra*D2),ra,Sc);
#endif
	}
      }
      ji++;
    }
    ji++;
  }
  
  if(Sc < 1.){
#ifdef DEBUG
    printf("Sc=%10.3e kfail=%d %c\n",Sc,kFail,7);
#endif
    if(kFail){
      ECRestore2DGeom();
    }
    else{
      ECvdrzCorrection(0);
    }
    if(Sc < 0.1){
      kFail++;
      return(1);
    }
  }
  kFail	=0;
  ESFail&=0x0000FFFF;
  return(0);
}

int ECCheckJacobHole(int n)
{
  int i,in,j,jni;
  int ki,k,kj;
  double rc[ESFp1],rs[ESFp1],rca[ESFp1],rsa[ESFp1];
  double *db,*db1a,*drc,*drca,*drs,*drsa;
  double D,D1,D2,bb,ba,ra,rgt,za,zgt;

  Sc	=1.;

  in	=ESNa1*n;
  ki	=ESnAF*n;
  jni	=ESnAP*n;

  db	=dbT+in;
  db1a	=dbT1a+in;
  drc	=drcT+ki;
  drca	=drcT1a+ki;
  drs	=drsT+ki;
  drsa	=drsT1a+ki;

  for(i=0; i < ESNa1; i++){
    bb		=db[i];
    ba		=db1a[i];
    ki		=i;
    rca[0]	=drca[ki];
    rsa[0]	=drsa[ki];
    for(k=1; k < Mr; k++){
      ki	+=ESNa1;
      rc[k]	=2.*k*drc[ki];
      rs[k]	=2.*k*drs[ki];
      rca[k]	=2.*drca[ki];
      rsa[k]	=2.*drsa[ki];
    }
    for(j=0; j < ESNp; j++){
      rgt	=0.;
      ra	=rca[0];
      kj	=0;
      for(k=1; k < Mr; k++){
	kj	+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	rgt	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
	ra	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
      }
      zgt	=bb*EScs1[j];
      za	=rsa[0]+ba*ESsn1[j];
      D		=rT1a[jni]*zT1gt[jni]-rT1gt[jni]*zT1a[jni];
      D1	=rT1a[jni]*zgt+ra*zT1gt[jni]-rT1gt[jni]*za-rgt*zT1a[jni];
      D2	=ra*zgt-rgt*za;
      if(0.55*D+D1+D2 < 0.){
	za	=sqrt(D1*D1-2.*D*D2);
	if(D1 > 0.){
	  ra	=D/(za-D1);
	}
	else{
	  ra	=D/(za-D1);
	  if(za+D1 < 0.){
	    za	=-D/(za+D1);
	    if(za > 0. && ra > za){
	      ra	=za;
	    }
	  }
	}	
	if(Sc > ra){
	  Sc	=ra;
	  if(i < 3){
	    Fail	|=0x00010000;
	  }
	  else{
	    if(i > ESNa-2){
	      Fail	|=0x00100000;
	    }
	    else{
	      Fail	|=0x01000000;
	    }
	  }
  
#ifdef DEBUG
	  printf("D[%3d][%2d][%2d]=%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n"
		 ,n,i,j,D,D+D1+D2,D+Sc*(D1+Sc*D2),D+ra*(D1+ra*D2),ra,Sc);
#endif
	}
      }
      jni++;
    }
    jni++;
  }
  
  if(Sc < 1.){
#ifdef DEBUG
    printf("Sc[n=%3d]=%10.3e kfail=%d %c\n",n,Sc,kFail,7);
#endif
    if(kFail){
      ECRestore2DGeom();
    }
    else{
      ECvdrzCorrection(n);
    }
    if(Sc < 0.1){
      kFail++;
      return(1);
    }
  }
  kFail	=0;
  ESFail&=0x0000FFFF;
  return(0);
}

int reseteInit(double*gyr,double*dgyr)
{
  int kmr,kmi,kYr,kYi;
  int i,m,m1,mm,k;
  double s,sr,si,gar,gai;
  double *dXr,*dXi,*yr,*yi,*dyr,*dyi;

  dXr	=Xr	+ES2Mp1;
  dXi	=Xi	+ES2Mp1;
  yr	=Ur	+ES2Mp1;
  yi	=Ui	+ES2Mp1;
  dyr	=dUr	+ES2Mp1;
  dyi	=dUi	+ES2Mp1;
  for(m=ESMp; m > 0; m--){
    kmr		=ES2Mp1*m;
    kmi		=kmr+M0Mm;
    kYr		=kmi+M0Mm;
    kYi		=kYr+M0Mm;
    for(k=0; k < ES2Mp1; k++){
      Ur[k]	=gyr[kmr];
      Ui[k]	=gyr[kmi];
      dUr[k]	=dgyr[kmr];
      dUi[k]	=dgyr[kmi];
      Xr[k]	=gyr[kYr];
      Xi[k]	=gyr[kYi];
      dXr[k]	=dgyr[kYr];
      dXi[k]	=dgyr[kYi];
      kmr++;
      kmi++;
      kYr++;
      kYi++;
    }
    k		=ESMp+m;
    sr		=dUr[k];
    si		=dUi[k];
    gar		=1./(sr*sr+si*si);
    sr		*=gar;
    si		*=gar;
    for(m1=0; m1 < ESMp1; m1++){
      kmr	=ES2Mp1*m1+ESMp-m;
      kmi	=kmr+M0Mm;
      gar	=dgyr[kmi]*si-dgyr[kmr]*sr;
      gai	=-dgyr[kmi]*sr-dgyr[kmr]*si;
      kmr	=ES2Mp1*m1;
      kmi	=kmr+M0Mm;
      kYr	=kmi+M0Mm;
      kYi	=kYr+M0Mm;
      k		=ES2Mp1;
      for(mm=0; mm < ES2Mp1; mm++){
	k--;
	gyr[kmr]	+=gar*Ur[k]+gai*Ui[k];
	gyr[kmi]	+=gai*Ur[k]-gar*Ui[k];
	dgyr[kmr]	+=gar*dUr[k]+gai*dUi[k];
	dgyr[kmi]	+=gai*dUr[k]-gar*dUi[k];
	gyr[kYr]	+=gar*Xr[k]+gai*Xi[k];
	gyr[kYi]	+=gai*Xr[k]-gar*Xi[k];
	dgyr[kYr]	+=gar*dXr[k]+gai*dXi[k];
	dgyr[kYi]	+=gai*dXr[k]-gar*dXi[k];
	kmr++;
	kmi++;
	kYr++;
	kYi++;
      }
    }

    kmr		=ES2Mp1*m;
    kmi		=kmr+M0Mm;
    kYr		=kmi+M0Mm;
    kYi		=kYr+M0Mm;
    for(k=0; k < ES2Mp1; k++){
      Ur[k]	=gyr[kmr];
      Ui[k]	=gyr[kmi];
      dUr[k]	=dgyr[kmr];
      dUi[k]	=dgyr[kmi];
      Xr[k]	=gyr[kYr];
      Xi[k]	=gyr[kYi];
      dXr[k]	=dgyr[kYr];
      dXi[k]	=dgyr[kYi];
      kmr++;
      kmi++;
      kYr++;
      kYi++;
    }
    k		=ESMp+m;
    sr		=dUr[k];
    si		=-dUi[k];
    gar		=1./(sr*sr+si*si);
    sr		*=gar;
    si		*=gar;
    for(m1=0; m1 < m; m1++){
      kmr	=ES2Mp1*m1+ESMp+m;
      kmi	=kmr+M0Mm;
      gar	=dgyr[kmi]*si-dgyr[kmr]*sr;
      gai	=-dgyr[kmi]*sr-dgyr[kmr]*si;
      kmr	=ES2Mp1*m1;
      kmi	=kmr+M0Mm;
      kYr	=kmi+M0Mm;
      kYi	=kYr+M0Mm;
      for(mm= 0;mm<ES2Mp1;mm++){
	gyr[kmr]	+=gar*Ur[mm]-gai*Ui[mm];
	gyr[kmi]	+=gai*Ur[mm]+gar*Ui[mm];
	dgyr[kmr]	+=gar*dUr[mm]-gai*dUi[mm];
	dgyr[kmi]	+=gai*dUr[mm]+gar*dUi[mm];
	gyr[kYr]	+=gar*Xr[mm]-gai*Xi[mm];
	gyr[kYi]	+=gai*Xr[mm]+gar*Xi[mm];
	dgyr[kYr]	+=gar*dXr[mm]-gai*dXi[mm];
	dgyr[kYi]	+=gai*dXr[mm]+gar*dXi[mm];
	kmr++;
	kmi++;
	kYr++;
	kYi++;
      }
    }
  }
  return(0);
}
#endif/*stg_2DMatching*/

#ifndef stg_2DInspect

#ifdef DEBUG
static double gxc[81*32],gxs[81*32],gsc[81*32],gss[81*32],gdz[81],gdb[81];
int ECgY2RZ()
{
  extern double *ESsr,*ESsz;
  double *gyr,*gyi;
  double *dgyr,*dgyi;

  double Rac[ESFp1],Ras[ESFp1];
  double rgY01a,za,b,ba,ra,gx;
  int i,j,k,kk;

  double r[ESNp1],z[ESNp1],r1[ESNp1],z1[ESNp1],*lr,*lz;

  rgY01a=1./ESgY02a[0];
  dz0[0]=2.*ECb1a[0]*EZdgpei[ESMp1]*rgY01a;
  drc[0]=-4.*(EZd1rcs[ESNa1]*EZdgper[ESMp1]-EZd1rsn[ESNa1]*EZdgpei[ESMp1])*rgY01a;

  b	=ESsb[ESNa];
  z[0]	=ESaZ0[0]-1.1*b;
  z[1]	=ESaZ0[0]+1.1*b;
  i	=ESNp1*ESNa;
  r[0]	=ESsr[i+ESNp/2];
  r[1]	=ESsr[i];
  b	=0.5*(r[1]-r[0]);
  r[0]	-=0.1*b;
  r[1]	+=0.1*b;

  Scale2D(1,2,r,z,2,6,0);

  r[0]	=ESaR0[0];
  z[0]	=ESaZ0[0];
  Plot2d(1,r,z,1,6,0,0,0);
  r[0]	+=drc[0];
  z[0]	+=dz0[0];
  Plot2d(1,r,z,1,6,0,4,0);

  k	=ESMp;
  gyr	=EZgper+k;
  gyi	=EZgpei+k;
  lr	=ESsr;
  lz	=ESsz;
  for(i=1; i < ESNa1; i++){
    gyr	+=M0Mm;
    gyi	+=M0Mm;
    lr	+=ESNp1;
    lz	+=ESNp1;
    kk		=i;
    Rac[0]	=rcT1a[kk];
    rgY01a	=1./(ESsa[i]*ESgY02a[i]);
    za		=rsT1a[kk];
    ba		=ESsb1a[i];
    gxc[0]	=0.;
    gxs[0]	=0.;
    for(k=1; k < ESFp1; k++){
      kk	+=ESNa1;
      Rac[k]	=2.*rcT1a[kk];
      Ras[k]	=2.*rsT1a[kk];
    }
    for(k=1; k < ESMp1; k++){
      gxc[k]	=-2.*gyr[k]*rgY01a;
      gxs[k]	=2.*gyi[k]*rgY01a;
    }
    for(j=0; j < ESNp1; j++){
      ra	=Rac[0];
      gx	=gxc[0];
      kk	=0;
      r1[j]	=lr[j]+drc[i];
      for(k=1; k < ESFp1; k++){
	kk	+=j;
	if(kk >= ESNp){
	  kk	-=ESNp;
	}
	ra	+=Rac[k]*EScs1[kk]+Ras[k]*ESsn1[kk];
	if(k < ESMp1){
	  gx	+=gxc[k]*EScs1[kk]+gxs[k]*ESsn1[kk];
	}
	if(k < ESMp){
	  r1[j]	+=drc[i+k*ESNa1]*EScs1[kk]+drs[i+k*ESNa1]*ESsn1[kk];
	}
      }
      r[j]	=lr[j]+ra*gx;
      z[j]	=lz[j]+(za+ba*ESsn1[j])*gx;
      z1[j]	=lz[j]+dz0[i]+dbT[i]*ESsn1[j];
    }
    Plot2d(1,lr,lz,ESNp1,6,0,12,0);
    Plot2d(1,r,z,ESNp1,6,0,4,0);
    Plot2d(1,r1,z1,ESNp1,6,0,14,0);
  }
  return(0);
}

static int Dmv=0,Di0=0,Di1=21,Dm0=-6,Dm1=6,DFl=0;
static double Dx[2]={0.,1.},Dy[2]={0.,0.};

int PlotdR(int k0,int k1,int i0, int i1)
{
  int i,k,n,ki,Fl;
  double dx;
  double x[2],y[2],Wx[ESNa*8+1],W[ESMp1][ESNa*8+1];
  double *p,*p1,*p2,t,f,d1,d2;
  char NmY[32],NmP[32];

  Fl	=(DFl/1)%3;
  dx	=0.125*ESsa[1];
  
  ki	=ESNa1*k0;
  p	=drsT	+ki;
  p1	=Fl ? (Fl == 1 ? drsT1a+ki : drsT2a+ki) : drsT+ki;
  p2	=drsT2a	+ki;
  
  y[0]	=p[0];
  y[1]	=y[0];
  x[0]	=ESsa[i0];
  x[1]	=ESsa[i1-1];
  n	=(i1-i0-1)*8+1;
  t	=ESsa[i0];
  for(i=0; i < n; i++){
    ESSetSplA(t);
    Wx[i]=t;
    ki	=0;
    for(k=k0; k < k1; k++){
      splRA2(&f,&d1,&d2,p+ki,p2+ki);
      if(Fl == 1) f=d1;
      if(Fl == 2) f=d2;
      W[k][i]	=f;
      if(y[0] > f) y[0]=f;
      if(y[1] < f) y[1]=f;
      ki	+=ESNa1;
    }
    t	+=dx;
  }
  if(Fl == 0 && k0 == 0){
    for(i=i0; i < i1; i++){
      if(y[0] > EZz0[i]) y[0]=EZz0[i];
      if(y[1] < EZz0[i]) y[1]=EZz0[i];
    }
  }
  sprintf(NmP,"%d<= k <= %d",k0,k1-1);
  SetPlotName("a","drs",NmP);
  Scale2D(3,2,x,y,2,6,2);
  ki	=i0;
  for(k=k0; k < k1; k++){
    Plot2d(3,ESsa+i0,p1+ki,i1-i0,6,0,4,0);
    if(k > k0){
      Plot2d(3,Wx,W[k],n,6,0,0,0);
    }
    else{
      Plot2d(3,Wx,W[k],n,6,0,14,0);
    }
    ki	+=ESNa1;
  }
  if(Fl == 0 && k0 == 0){
    Plot2d(3,ESsa+i0,EZz0+i0,i1-i0,6,0,9,0);
  }

  ki	=ESNa1*k0;
  p	=drcT	+ki;
  p1	=Fl ? (Fl == 1 ? drcT1a+ki : drcT2a+ki) : drcT+ki;
  p2	=drcT2a	+ki;
  k0	=0;
  y[0]	=p[0];
  y[1]	=y[0];
  n	=(i1-i0-1)*8+1;
  t	=ESsa[i0];
  for(i=0; i < n; i++){
    ESSetSplA(t);
    Wx[i]=t;
    ki	=0;
    for(k=k0; k < k1; k++){
      splRA2(&f,&d1,&d2,p+ki,p2+ki);
      if(Fl == 1) f=d1;
      if(Fl == 2) f=d2;
      W[k][i]	=f;
      if(y[0] > f) y[0]=f;
      if(y[1] < f) y[1]=f;
      ki	+=ESNa1;
    }
    t	+=dx;
  }
  sprintf(NmP,"%d<= k <= %d",k0,k1-1);
  SetPlotName("a","drc",NmP);
  Scale2D(2,2,x,y,2,6,2);
  ki	=i0;
  for(k=k0; k < k1; k++){
    Plot2d(2,ESsa+i0,p1+ki,i1-i0,6,0,4,0);
    i	=k == k0 ? 14 : 0;
    Plot2d(2,Wx,W[k],n,6,0,i,0);
    ki	+=ESNa1;
  }

  y[0]	=0.;
  y[1]	=0.;
  x[0]	=ESsa[i0];
  x[1]	=ESsa[i1-1];
  p	=daYer+M0Mm*i0+ESMp;
  for(i=i0; i < i1; i++){
    ki	=ESNa1*k0+i;
    for(k=k0; k < k1; k++){
      t	=-ESsa[i]*(ESLc[ki]*ESjs[i]+ESVc[ki]*ESjp[i]);
      W[k][i]	=t-p[k];
      if(k)	t -=p[k];
      if(y[0] > t) y[0]	=t;
      if(y[1] < t) y[1]	=t;
      ki	+=ESNa1;
    }
    p	+=M0Mm;
  }
  SetPlotName("a","RHS-LHS","Fourier GSh");
  Scale2d(5,x,y,2,6,2);
  for(k=k0; k < k1; k++){
    i	=k == k0 ? 14 : 0;
    Plot2d(5,ESsa+i0,W[k],i1-i0,6,0,i,0);
  }
  p	=daYei+M0Mm*i0+ESMp;
  for(i=i0; i < i1; i++){
    ki	=ESNa1*k0+i;
    for(k=k0; k < k1; k++){
      t	=ESsa[i]*(ESLs[ki]*ESjs[i]+ESVs[ki]*ESjp[i]);
      W[k][i]	=t-p[k];
      if(k)	t -=p[k];
      if(y[0] > t) y[0]	=t;
      if(y[1] < t) y[1]	=t;
      ki	+=ESNa1;
    }
    p	+=M0Mm;
  }
  for(k=k0; k < k1; k++){
    Plot2d(5,ESsa+i0,W[k],i1-i0,6,0,5,0);
  }
  return(0);
}

int ScalegY()
{
  int Fl,Fd,Fr;
  int i,k,ki;
  double *ld;

  if(Dm0 < -ESMp) Dm0=-ESMp;
  if(Dm1 >  ESMp) Dm1= ESMp;
  if(Dm0 > Dm1) Dm0=Dm1;

  Fr	=DFl%2;
  Fd	=(DFl/2)%3;
  Fl	=DFl/6;

  switch(Fl){
  case 0:
    if(Fd == 0){
      ld	=Fr == 0 ? EZgper : EZgpei;
    }
    else{
      ld	=Fr == 0 ? EZdgper : EZdgpei;
    }
    break;
  default :
  case 1:
    if(Fd == 0){
      ld	=Fr == 0 ? aYer : aYei;
    }
    else{
      ld	=Fr == 0 ? daYer : daYei;
    }
    break;
  }
  ld	+=ES2Mp1*Dmv+ESMp;

  ki	=M0Mm*Di0+Dm0;
  Dy[0]	=ld[ki];
  Dy[1]	=Dy[0];
  Dx[0]	=ESsa[Di0];
  Dx[1]	=ESsa[Di1-1];
  for(i=Di0; i < Di1; i++){
    ki	=M0Mm*i+Dm0;
    for(k=Dm0; k <= Dm1; k++){
      if(k ||  Fr || Dmv){
	if(Dy[0] > ld[ki]) Dy[0]=ld[ki];
	if(Dy[1] < ld[ki]) Dy[1]=ld[ki];
      }
      ki++;
    }
  }

  if(Fr == 0 && Dmv == 0 && Dm0 == 0 && Dm1 == 0){
    for(i=Di0; i < Di1; i++){
      ki	=M0Mm*i;
      if(Dy[0] > ld[ki]) Dy[0]=ld[ki];
      if(Dy[1] < ld[ki]) Dy[1]=ld[ki];
    }
  }
  return(0);
}

void PlotgY()
{
  int Fl,Fd,Fr;
  int k,i,ki;
  double W[ESNa1];
  double *ld;

  Fr	=DFl%2;
  Fd	=(DFl/2)%3;
  Fl	=DFl/6;

  switch(Fl){
  case 0:
    if(Fd == 0){
      ld	=Fr == 0 ? EZgper : EZgpei;
    }
    else{
      ld	=Fr == 0 ? EZdgper : EZdgpei;
    }
    break;
  default :
  case 1:
    if(Fd == 0){
      ld	=Fr == 0 ? aYer : aYei;
    }
    else{
      ld	=Fr == 0 ? daYer : daYei;
    }
    break;
  }
  ld	+=ES2Mp1*Dmv+ESMp;
  ki	=M0Mm*Di0+Dm0;
  Dy[0]	=ld[ki];
  Dy[1]	=Dy[0];
  Dx[0]	=ESsa[Di0];
  Dx[1]	=ESsa[Di1-1];
  for(k=Dm0; k <= Dm1; k++){
    ki	=k+M0Mm*Di0;
    if(k ||  Fr || Dmv){
      for(i=Di0; i < Di1; i++){
	W[i]	=ld[ki];
	ki	+=M0Mm;
      }
      ZPlotPolyLine(ESsa+Di0,W+Di0,Di1-Di0);
    }
  }
  if(Fr == 0 && Dmv == 0){
    ZColor(14);
    ki	=M0Mm*Di0;
    for(i=Di0; i < Di1; i++){
      W[i]	=Dm0 || Dm1 ? ld[ki]*1e-3 : ld[ki];
      ki	+=M0Mm;
    }
    ZPlotPolyLine(ESsa+Di0,W+Di0,Di1-Di0);
    k	=Di1-1;
    for(i=Di0; i < Di1; i++){
      W[i]	-=W[k];
    }
    ZColor(4);
    ZPlotPolyLine(ESsa+Di0,W+Di0,Di1-Di0);
  }
  return;
}

int InsFundSol()
{
  double *xc;
  double *yc,*ys;
  static int mv=0,i0=0,i1=21,k0=0,k1=6,Fl=0,Fd=0,Fr=0;
  double a[ESNa1],Wr[ES2Mp1][ESNa1],Wi[ES2Mp1][ESNa1];

  ECgY2RZ();
#ifndef Tbl_InsFSol
  DFl	=2*(3*Fl+Fd)+Fr;
  ScalegY();
  
  if(k1 >= ESMp) k1 =ESMp-1;
  if(k1 < 0) k1 =0;
  if(k0 < 0) k0 =0;
  if(k0 > k1) k0=k1;
  PlotdR(k0,k1+1,Di0,Di1);

  CbUserPlot	=PlotgY;
#endif /*Tbl_InsFSol*/
  return(0);
}
#else/*DEBUG*/
int InsFundSol()
{
  return(0);
}
#endif

int ESCheckGShFourier00(int i0,int i1,int m0,int m1,int Fl)
{
  extern char ShotNm[];
  int i,k,ki;
  double Kc,Kca,Ks,Ksa;
  double Nct,Nst;
  double Lc,Ls;
  double Vc,Vs;
#ifdef DEBUG
  double LHSc[ESNa1*ESMp1],RHSc[ESNa1*ESMp1];
  double LHSs[ESNa1*ESMp1],RHSs[ESNa1*ESMp1];
#else
  double LHSc[129*17],RHSc[129*17];
  double LHSs[129*17],RHSs[129*17];
#endif
  double K,Ka,Nt,L,V;
  double *lL,*lR,A;
 
  int M1,j,m;
  double Rc[65],Rs[65],Rac[65],Ras[65],Raac[65],Raas[65],b,ba,baa,cs,sn;
  double r[128],ra[128],raa,rt[128],rat[128]
    ,za[128],zaa,zt[128],zat[128],D[128],Da[128],G[128];

  ki	=0;
  for(k=0; k < ESMp1; k++){
    LHSc[ki]	=0.;
    LHSs[ki]	=0.;
    RHSc[ki]	=0.;
    RHSs[ki]	=0.;
    ki	+=ESNa1;
  }
  for(i=1; i < ESNa1; i++){
    A	=ESsa[i];
    ki	=i;

    /* Getting metric coefficients{*/
    b	=ESsb[i];
    ba	=ESsb1a[i];
    baa	=ESsb2a[i];
    Rc[0]	=rcT[ki];
    Rac[0]	=rcT1a[ki];
    Raac[0]	=rcT2a[ki];
    Rs[0]	=rsT[ki];
    Ras[0]	=rsT1a[ki];
    Raas[0]	=rsT2a[ki];
    m	=i;
    M1	=2;
    for(k=1; k < ESFp1; k++){
      m	+=ESNa1;
      Rc[k]	=2.*rcT[m];
      Rac[k]	=2.*rcT1a[m];
      Raac[k]	=2.*rcT2a[m];
      Rs[k]	=2.*rsT[m];
      Ras[k]	=2.*rsT1a[m];
      Raas[k]	=2.*rsT2a[m];
      if(Rc[k] != 0. || Rs[k] != 0.) M1=k+1;
    }
    for(j=0; j < ESNp; j++){
      r[j]	=ESaR0[0]+Rc[0];
      ra[j]	=Rac[0];
      raa	=Raac[0];
      rt[j]	=0.;
      rat[j]	=0.;
      za[j]	=Ras[0]+ba*ESsn1[j];
      zaa	=Raas[0]+baa*ESsn1[j];
      zt[j]	=b*EScs1[j];
      zat[j]	=ba*EScs1[j];
      m		=0;
      for(k=1; k < ESFp1; k++){
	m	+=j;
	if(m >= ESNp) m -=ESNp;
	cs	=EScs1[m];
	sn	=ESsn1[m];
	r[j]	+=Rc[k]*cs+Rs[k]*sn;
	ra[j]	+=Rac[k]*cs+Ras[k]*sn;
	raa	+=Raac[k]*cs+Raas[k]*sn;
	rt[j]	+=k*(Rs[k]*cs-Rc[k]*sn);
	rat[j]	+=k*(Ras[k]*cs-Rac[k]*sn);
      }
      D[j]	=ra[j]*zt[j]-rt[j]*za[j];
      G[j]	=(rt[j]*rt[j]+zt[j]*zt[j])/(r[j]*D[j]);
      Da[j]	=raa*zt[j]+ra[j]*zat[j]-rat[j]*za[j]-rt[j]*zaa;
    }
    ESP2F(wg22c,wg22s,G,ES2Mp);
    for(j=0; j < ESNp; j++){
      G[j]	=(2.*rt[j]*rat[j]+2.*zt[j]*zat[j]
		  -G[j]*(Da[j]*r[j]+D[j]*ra[j]))/(r[j]*D[j]);
    }
    ESP2F(wg11c,wg11s,G,ES2Mp);
    for(j=0; j < ESNp; j++){
      G[j]	=(ra[j]*rt[j]+za[j]*zt[j])/(r[j]*D[j]);
    }
    ESP2F(wg12c,wg12s,G,ES2Mp);
    for(j=0; j < ESNp; j++){
      G[j]	=D[j]*ESaR0[0]/r[j];
    }
    ESP2F(wuc,wus,G,ES2Mp);
    for(j=0; j < ESNp; j++){
      G[j]	=D[j]*(r[j]/ESaR0[0]-ESaR0[0]/r[j]);
    }
    ESP2F(wvc,wvs,G,ES2Mp);
    /*}*/

    Kc	=wg22c[0];
    Ks	=wg22s[0];
    Kca	=wg11c[0];
    Ksa	=wg11s[0];
    
    Nct	=0.;
    Nst	=0.;
    Lc	=wuc[0];
    Ls	=wvs[0];
    Vc	=wvc[0];
    Vs	=wvs[0];

    LHSc[ki]	=(Kc+A*Kca-A*Nct)*ESdgY[i]+A*Kc*ESdgY1a[i];
    m	=ESMp+M0Mm*i;
    LHSc[ki]	=daYer[m];
#ifdef H
#endif
    LHSs[ki]	=0.;
    RHSc[ki]	=-Lc*ESjs[i]-Vc*ESjp[i];
    RHSs[ki]	=0.;
    LHSc[ki]	-=RHSc[ki];
    LHSs[ki]	-=RHSs[ki];
    for(k=1; k < ESMp1; k++){
      ki	+=ESNa1;
      Kc	=wg22c[k];
      Ks	=wg22s[k];
      Kca	=wg11c[k];
      Ksa	=wg11s[k];
      Nct	= k*wg12s[k];
      Nst	=-k*wg12c[k];
      Lc	=wuc[k];
      Ls	=wus[k];
      Vc	=wvc[k];
      Vs	=wvs[k];
      LHSc[ki]	=(Kc+A*(Kca-Nct))*ESdgY[i]+A*Kc*ESdgY1a[i];
      LHSs[ki]	=(Ks+A*(Ksa-Nst))*ESdgY[i]+A*Ks*ESdgY1a[i];

      m	=ESMp+M0Mm*i+k;
      LHSc[ki]	=daYer[m];
      LHSs[ki]	=-daYei[m];
#ifdef H
#endif
      RHSc[ki]	=-Lc*ESjs[i]-Vc*ESjp[i];
      RHSs[ki]	=-Ls*ESjs[i]-Vs*ESjp[i];
      LHSc[ki]	-=RHSc[ki];
      LHSs[ki]	-=RHSs[ki];
    }
  }
#ifdef DEBUG
  {
    extern int ESEqSolvFl;
    double x[2],y[ESNa1],y0[ESNa1];
    double r;
    char ln[32];

    if(i0 < 0){
      i0	=0;
    }
    if(i1 > ESNa1){
      i1	=ESNa1;
    }
    if(m0 < 0){
      m0	=0;
    }
    if(m1 > ESMp1){
      m1	=ESMp1;
    }

    y[1]	=-20.;
    x[0]	=ESsa[i0];
    x[1]	=ESsa[i1-1];
    ki		=ESNa1*ESMp1;
    for(i=0; i < ki; i++){
      r	=log10(fabs(RHSc[i])+1e-20);
      if(y[1] < r) y[1]	=r;
      r	=log10(fabs(RHSs[i])+1e-20);
      if(y[1] < r) y[1]	=r;
      r	=log10(fabs(LHSc[i])+1e-20);
      if(y[1] < r) y[1]	=r;
      r	=log10(fabs(LHSs[i])+1e-20);
      if(y[1] < r) y[1]	=r;
    }
    y[0]	=y[1]-5.;
    if(Fl){
      static int inx=1,iny=1;
      static int kPS=1;
      static char PSFileNm[16];
      if(Fl > 1){
	if(Fl != 4){
	  inx	=Fl-1;
	  iny	=Fl-1;
	}
	else{
	  inx	=4;
	  iny	=4;
	}
      }
      if(kPS && *ShotNm != '\0'){
	sprintf(PSFileNm,"%s.ps",ShotNm);
	PSSetFileName(PSFileNm);
	kPS	=0;
      }
      PSFileOpen(inx,iny);
    }
    sprintf(ln,"i=%d<%d m=%d<%d (%s) %s",i0,i1,m0,m1
	    ,(ESEqSolvFl&0x0F) ? "RungK" : "Lsode",ShotNm);
    SetPlotName("sqrt(gF)","GSh-Fourier log10|LHS-RHS|",ln);
    Scale2d(6,x,y,2,6,2);
    if(Fl) PSNewFrame();

    ki	=ESNa1*m0;
    for(k=m0; k < m1; k++){
      lL	=LHSc+ki;
      lR	=RHSc+ki;
      for(i=i0; i < i1; i++){
	y[i]	=log10(fabs(lL[i])+1e-20);
	y0[i]	=log10(fabs(lR[i])+1e-20);
      }
      if(k < m0+2){
	Plot2d(6,ESsa+i0,y0+i0,i1-i0,6,0,0,0);
      }
      else{
	Plot2d(6,ESsa+i0,y0+i0,i1-i0,6,0,12,0);
      }
      Plot2d(6,ESsa+i0,y+i0,i1-i0,6,0,4,0);
      if(Fl){
	if(k < m0+2){
	  PSPlot2d(ESsa+i0,y0+i0,i1-i0,6,0,0,0);
	}
	else{
	  PSPlot2d(ESsa+i0,y0+i0,i1-i0,6,0,12,0);
	}
	PSPlot2d(ESsa+i0,y+i0,i1-i0,6,0,4,0);
      }

      lL	=LHSs+ki;
      lR	=RHSs+ki;
      for(i=i0; i < i1; i++){
	y[i]	=log10(fabs(lL[i])+1e-20);
	y0[i]	=log10(fabs(lR[i])+1e-20);
      }
      Plot2d(6,ESsa+i0,y0+i0,i1-i0,6,0,12,0);
      Plot2d(6,ESsa+i0, y+i0,i1-i0,6,0,14,0);
      if(Fl){
	PSPlot2d(ESsa+i0,y0+i0,i1-i0,6,0,12,0);
	PSPlot2d(ESsa+i0, y+i0,i1-i0,6,0,14,0);
      }
      ki	+=ESNa1;
    }
    CbFlush();
    if(Fl){
      PSFileClose();
      ESErrorInGSh();
    }
  }
#endif
  return(0);
}

#endif/*stg_2DInspect*/
#endif/*mzl_2DEq*/

#ifdef XWIN
int ESCheckGShFourier(int i0,int i1,int m0,int m1,int Fl)
{
  extern char ShotNm[];

  int L;
  int i,j,k,m,n;

  /* Getting metric coefficients{*/
  int M1;
  int Mp,Mp1;
  int mp,mp1;
  int ki;
  double LHSc[ESNa1*ES2Mp1],RHSc[ESNa1*ES2Mp1];
  double LHSs[ESNa1*ES2Mp1],RHSs[ESNa1*ES2Mp1];
  double *lL,*lR,A,rA;


  L	=1;
  switch(L){
  case 0:
    Mp	=ESMp;
    mp	=ESMp;
    break;
  case 1:
    Mp	=ESMp;
    mp	=0;
    break;
  case 2:
    Mp	=ES2Mp;
    mp	=0;
    break;
  default :
    Mp	=ESMp;
    mp	=ESMp;
    break;
  }
  m1	=Mp1;
  Mp1	=Mp+1;
  mp1	=mp+1;

  ki	=0;
  for(k=0; k < Mp1; k++){
    LHSc[ki]	=0.;
    LHSs[ki]	=0.;
    RHSc[ki]	=0.;
    RHSs[ki]	=0.;
    ki	+=ESNa1;
  }

  {
    double *lLr,*lLi,*lRr,*lRi,*mUr,*mUi;
    double *pr,*pi,*yr,*yi,*dpr,*dpi,*dyr,*dyi;
    double sr,si,tr,ti,wj_s,wj_p,dj_s,dj_p;

    mUr	=Ur+ES2Mp;
    mUi	=Ui+ES2Mp;

    pr	=gpr+ESMp;
    pi	=pr+M0Mm;
    yr	=pi+M0Mm;
    yi	=yr+M0Mm;
    dpr	=dgpr+ESMp;
    dpi	=dpr+M0Mm;
    dyr	=dpi+M0Mm;
    dyi	=dyr+M0Mm;

    for(i=1; i < ESNa1; i++){
      A	=ESsa[i];
      rA	=1./A;

      dj_s	=ESjs1a[i];
      dj_p	=ESjp1a[i];
      wj_s	=ESjs[i];
      wj_p	=ESjp[i];
      dj_s	=0.;
      dj_p	=0.;

      n	=i;
      wuc[0]	=ESsa[i]*ESLc[n];
      wus[0]	=0.;
      wvc[0]	=ESsa[i]*ESVc[n];
      wvs[0]	=0.;
      wg22c[0]=ESsa[i]*ESg22c[n];
      wg22s[0]=0.;
      wg12c[0]=ESg12c[n];
      wg12s[0]=0.;
      wg11c[0]=ESg11c[n]/ESsa[i];
      wg11s[0]=0.;

      for(m=1; m < ES2Mp1; m++){
	n	+=ESNa1;
	k	=-m;
	wuc[k]	=ESsa[i]*ESLc[n];
	wus[k]	=ESsa[i]*ESLs[n];
	wuc[m]	= wuc[k];
	wus[m]	=-wus[k];
	wvc[k]	=ESsa[i]*ESVc[n];
	wvs[k]	=ESsa[i]*ESVs[n];
	wvc[m]	= wvc[k];
	wvs[m]	=-wvs[k];
	wg22c[k]	=ESsa[i]*ESg22c[n];
	wg22s[k]	=ESsa[i]*ESg22s[n];
	wg22c[m]	= wg22c[k];
	wg22s[m]	=-wg22s[k];
	wg12c[k]	=ESg12c[n];
	wg12s[k]	=ESg12s[n];
	wg12c[m]	= wg12c[k];
	wg12s[m]	=-wg12s[k];
	wg11c[k]	=ESg11c[n]/ESsa[i];
	wg11s[k]	=ESg11s[n]/ESsa[i];
	wg11c[m]	= wg11c[k];
	wg11s[m]	=-wg11s[k];
      }

      mUr[0]	=-dj_s*wuc[0]-dj_p*wvc[0];
      mUi[0]	=0.;
      for(k=1; k < ES2Mp1; k++){
	mUr[k]	=-dj_s*wuc[k]-dj_p*wvc[k];
	mUi[k]	=-dj_s*wus[k]-dj_p*wvs[k];
	mUr[-k]	= mUr[k];
	mUi[-k]	=-mUi[k];
      }

      k	=M0Mm*i+ESMp;
      pr	=EZgper+k;
      pi	=EZgpei+k;
      dpr	=EZdgper+k;
      dpi	=EZdgpei+k;
      yr	=aYer+k;
      yi	=aYei+k;
      dyr	=daYer+k;
      dyi	=daYei+k;

      lLr	=LHSr	+Mp;
      lLi	=LHSi	+Mp;
      lRr	=RHSr	+Mp;
      lRi	=RHSi	+Mp;

      for(n=0; n < 1; n++){
	/* Calculation of LHS{*/
	for(m=-Mp; m < Mp1; m++){
	  lLr[m]	=0.;
	  lLi[m]	=0.;
	  lRr[m]	=0.;
	  lRi[m]	=0.;
	  if(n == 0){
	    lRr[m]	=-wj_s*wuc[m]-wj_p*wvc[m];
	    lRi[m]	=-wj_s*wus[m]-wj_p*wvs[m];
	  }
	}
	for(m=-ESMp; m < ESMp1; m++){
	  lLr[m]	=dyr[m];
	  lLi[m]	=dyi[m];
	}
	for(j=-mp; j < mp1; j++){
	  sr	=dpr[j];
	  si	=dpi[j];
	  for(m=-Mp; m < Mp1; m++){
	    if(m){
	      k	=m-j;
	      tr	= m*wg12s[k];
	      ti	=-m*wg12c[k];
	      lLr[m]+=sr*tr-si*ti;
	      lLi[m]+=si*tr+sr*ti;
	    }
	  }
	  if(j){
	    sr	=pr[j];
	    si	=pi[j];
	    for(m=-ESMp; m < ESMp1; m++){
	      if(m){
		k	=m-j;
		tr	=mUr[k]+m*wg11c[k]*j;
		ti	=mUi[k]+m*wg11s[k]*j;
		lLr[m]	-=sr*tr-si*ti;
		lLi[m]	-=si*tr+sr*ti;
	      }
	    }
	  }
	}
	/*}*/

#ifdef H
	for(m=-ESMp; m < ESMp1; m++){
	  EZout("sidddddd","RHS",m,lRi[m],lRr[m],lLi[m],lLr[m]
	      ,lLi[m]-lRi[m],lLr[m]-lRr[m]);
	}
#endif
	ki	=i;
	for(m=0; m < Mp1; m++){
	  RHSc[ki]	=lRr[m];
	  RHSs[ki]	=lRi[m];
	  LHSc[ki]	=lLr[m]-lRr[m];
	  LHSs[ki]	=lLi[m]-lRi[m];
	  ki	+=ESNa1;
	}
	pr	+=ES2Mp1;
	pi	+=ES2Mp1;
	yr	+=ES2Mp1;
	yi	+=ES2Mp1;

	dpr	+=ES2Mp1;
	dpi	+=ES2Mp1;
	dyr	+=ES2Mp1;
	dyi	+=ES2Mp1;
      }
    }
  }
  {
    extern int ESEqSolvFl;
    double x[2],y[ESNa1],y0[ESNa1];
    double r;
    char ln[32];

    if(i0 < 0){
      i0	=0;
    }
    if(i1 > ESNa1){
      i1	=ESNa1;
    }
    if(m0 < 0){
      m0	=0;
    }
    if(m1 > Mp1){
      m1	=Mp1;
    }

    y[1]	=-20.;
    x[0]	=ESsa[i0];
    x[1]	=ESsa[i1-1];
    ki		=ESNa1*Mp1;
    for(i=0; i < ki; i++){
      r	=log10(fabs(RHSc[i])+1e-20);
      if(y[1] < r) y[1]	=r;
      r	=log10(fabs(RHSs[i])+1e-20);
      if(y[1] < r) y[1]	=r;
      r	=log10(fabs(LHSc[i])+1e-20);
      if(y[1] < r) y[1]	=r;
      r	=log10(fabs(LHSs[i])+1e-20);
      if(y[1] < r) y[1]	=r;
    }
    y[0]	=y[1]-10.;
    if(Fl){
      static int inx=1,iny=1;
      static int kPS=1;
      static char PSFileNm[16];
      if(Fl > 1){
	if(Fl != 4){
	  inx	=Fl-1;
	  iny	=Fl-1;
	}
	else{
	  inx	=4;
	  iny	=4;
	}
      }
      if(kPS && *ShotNm != '\0'){
	sprintf(PSFileNm,"%s.ps",ShotNm);
	PSSetFileName(PSFileNm);
	kPS	=0;
      }
      PSFileOpen(inx,iny);
    }
    sprintf(ln,"i=%d<%d m=%d<%d (%s) %s",i0,i1,m0,m1
	    ,(ESEqSolvFl&0x0F) ? "RungK" : "Lsode",ShotNm);
    SetPlotName("sqrt(gF)","GSh-Fourier log10|LHS-RHS|",ln);
    Scale2d(6,x,y,2,6,2);
    if(Fl) PSNewFrame();

    ki	=ESNa1*m0;
    for(k=m0; k < m1; k++){
      lL	=LHSc+ki;
      lR	=RHSc+ki;
      for(i=i0; i < i1; i++){
	y[i]	=log10(fabs(lL[i])+1e-20);
	y0[i]	=log10(fabs(lR[i])+1e-20);
      }
      if(k < m0+2){
	Plot2d(6,ESsa+i0,y0+i0,i1-i0,6,0,0,0);
      }
      else{
	Plot2d(6,ESsa+i0,y0+i0,i1-i0,6,0,12,0);
      }
      Plot2d(6,ESsa+i0,y+i0,i1-i0,6,0,4,0);
      if(Fl){
	if(k < m0+2){
	  PSPlot2d(ESsa+i0,y0+i0,i1-i0,6,0,0,0);
	}
	else{
	  PSPlot2d(ESsa+i0,y0+i0,i1-i0,6,0,12,0);
	}
	PSPlot2d(ESsa+i0,y+i0,i1-i0,6,0,4,0);
      }

      lL	=LHSs+ki;
      lR	=RHSs+ki;
      for(i=i0; i < i1; i++){
	y[i]	=log10(fabs(lL[i])+1e-20);
	y0[i]	=log10(fabs(lR[i])+1e-20);
      }
      Plot2d(6,ESsa+i0,y0+i0,i1-i0,6,0,12,0);
      Plot2d(6,ESsa+i0, y+i0,i1-i0,6,0,14,0);
      if(Fl){
	PSPlot2d(ESsa+i0,y0+i0,i1-i0,6,0,12,0);
	PSPlot2d(ESsa+i0, y+i0,i1-i0,6,0,14,0);
      }
      ki	+=ESNa1;
    }
    CbFlush();
    if(Fl){
      PSFileClose();
      ESErrorInGSh();
    }
  }


  {
    double *lLr,*lLi,*lRr,*lRi,*mUr,*mUi;
    double *pr,*pi,*yr,*yi,*dpr,*dpi,*dyr,*dyi;
    double sr,si,tr,ti,wj_s,wj_p,dj_s,dj_p;

    mUr	=Ur+ES2Mp;
    mUi	=Ui+ES2Mp;

    pr	=gpr+ESMp;
    pi	=pr+M0Mm;
    yr	=pi+M0Mm;
    yi	=yr+M0Mm;
    dpr	=dgpr+ESMp;
    dpi	=dpr+M0Mm;
    dyr	=dpi+M0Mm;
    dyi	=dyr+M0Mm;

    for(i=1; i < ESNa1; i++){
      A	=ESsa[i];
      rA	=1./A;

      n	=i;
      wg22c[0]=ESsa[i]*ESg22c[n];
      wg22s[0]=0.;
      wg12c[0]=ESg12c[n];
      wg12s[0]=0.;
      for(m=1; m < ES2Mp1; m++){
	n	+=ESNa1;
	k	=-m;
	wg22c[k]	=ESsa[i]*ESg22c[n];
	wg22s[k]	=ESsa[i]*ESg22s[n];
	wg22c[m]	= wg22c[k];
	wg22s[m]	=-wg22s[k];
	wg12c[k]	=ESg12c[n];
	wg12s[k]	=ESg12s[n];
	wg12c[m]	= wg12c[k];
	wg12s[m]	=-wg12s[k];
      }

      k	=M0Mm*i+ESMp;
      pr	=EZgper+k;
      pi	=EZgpei+k;
      dpr	=EZdgper+k;
      dpi	=EZdgpei+k;
      yr	=aYer+k;
      yi	=aYei+k;
      dyr	=daYer+k;
      dyi	=daYei+k;

      lLr	=LHSr	+ESMp;
      lLi	=LHSi	+ESMp;
      lRr	=RHSr	+ESMp;
      lRi	=RHSi	+ESMp;

      for(n=0; n < 1; n++){
	/* Calculation of $dgy${*/
	for(m=-Mp; m < Mp1; m++){
	  lLr[m]	=0.;
	  lLi[m]	=0.;
	  lRr[m]	=0.;
	  lRi[m]	=0.;
	}
	for(m=-ESMp; m < ESMp1; m++){
	  lRr[m]	=yr[m];
	  lRi[m]	=yi[m];
	}
	for(m=-Mp; m < Mp1; m++){
	  for(j=-mp; j < mp1; j++){
	    sr	=j*pr[j];
	    si	=j*pi[j];
	    k		=m-j;
	    tr	=-wg12s[k];
	    ti	= wg12c[k];
	    lRr[m]	+=sr*tr-si*ti;
	    lRi[m]	+=si*tr+sr*ti;
	    sr	=dpr[j];
	    si	=dpi[j];
	    tr	=wg22c[k];
	    ti	=wg22s[k];
	    lLr[m]	+=sr*tr-si*ti;
	    lLi[m]	+=si*tr+sr*ti;
	  }
	}

#ifdef H
	for(m=-ESMp; m < ESMp1; m++){
	  EZout("sidddddd","RHS",m,lRi[m],lRr[m],lLi[m],lLr[m]
	      ,lLi[m]-lRi[m],lLr[m]-lRr[m]);
	}
#endif
	ki	=i;
	for(m=0; m < Mp1; m++){
	  RHSc[ki]	=lRr[m];
	  RHSs[ki]	=lRi[m];
	  LHSc[ki]	=lLr[m]-lRr[m];
	  LHSs[ki]	=lLi[m]-lRi[m];
	  ki	+=ESNa1;
	}

	pr	+=ES2Mp1;
	pi	+=ES2Mp1;
	yr	+=ES2Mp1;
	yi	+=ES2Mp1;

	dpr	+=ES2Mp1;
	dpi	+=ES2Mp1;
	dyr	+=ES2Mp1;
	dyi	+=ES2Mp1;

	lLr	+=ES2Mp1;
	lLi	+=ES2Mp1;
	lRr	+=ES2Mp1;
	lRi	+=ES2Mp1;
      }
    }
  }
  {
    extern int ESEqSolvFl;
    double x[2],y[ESNa1],y0[ESNa1];
    double r;
    char ln[32];

    if(i0 < 0){
      i0	=0;
    }
    if(i1 > ESNa1){
      i1	=ESNa1;
    }
    if(m0 < 0){
      m0	=0;
    }
    if(m1 > Mp1){
      m1	=Mp1;
    }

    y[1]	=-20.;
    x[0]	=ESsa[i0];
    x[1]	=ESsa[i1-1];
    ki		=ESNa1*Mp1;
    for(i=0; i < ki; i++){
      r	=log10(fabs(RHSc[i])+1e-20);
      if(y[1] < r) y[1]	=r;
      r	=log10(fabs(RHSs[i])+1e-20);
      if(y[1] < r) y[1]	=r;
      r	=log10(fabs(LHSc[i])+1e-20);
      if(y[1] < r) y[1]	=r;
      r	=log10(fabs(LHSs[i])+1e-20);
      if(y[1] < r) y[1]	=r;
    }
    y[0]	=y[1]-10.;
    if(Fl){
      static int inx=1,iny=1;
      static int kPS=1;
      static char PSFileNm[16];
      if(Fl > 1){
	if(Fl != 4){
	  inx	=Fl-1;
	  iny	=Fl-1;
	}
	else{
	  inx	=4;
	  iny	=4;
	}
      }
      if(kPS && *ShotNm != '\0'){
	sprintf(PSFileNm,"%s.ps",ShotNm);
	PSSetFileName(PSFileNm);
	kPS	=0;
      }
      PSFileOpen(inx,iny);
    }
    sprintf(ln,"i=%d<%d m=%d<%d (%s) %s",i0,i1,m0,m1
	    ,(ESEqSolvFl&0x0F) ? "RungK" : "Lsode",ShotNm);
    SetPlotName("sqrt(gF)","GSh-Fourier log10|LHS1-RHS1|",ln);
    Scale2d(5,x,y,2,6,2);
    if(Fl) PSNewFrame();

    ki	=ESNa1*m0;
    for(k=m0; k < m1; k++){
      lL	=LHSc+ki;
      lR	=RHSc+ki;
      for(i=i0; i < i1; i++){
	y[i]	=log10(fabs(lL[i])+1e-20);
	y0[i]	=log10(fabs(lR[i])+1e-20);
      }
      if(k < m0+2){
	Plot2d(5,ESsa+i0,y0+i0,i1-i0,6,0,0,0);
      }
      else{
	Plot2d(5,ESsa+i0,y0+i0,i1-i0,6,0,12,0);
      }
      Plot2d(5,ESsa+i0,y+i0,i1-i0,6,0,4,0);
      if(Fl){
	if(k < m0+2){
	  PSPlot2d(ESsa+i0,y0+i0,i1-i0,6,0,0,0);
	}
	else{
	  PSPlot2d(ESsa+i0,y0+i0,i1-i0,6,0,12,0);
	}
	PSPlot2d(ESsa+i0,y+i0,i1-i0,6,0,4,0);
      }

      lL	=LHSs+ki;
      lR	=RHSs+ki;
      for(i=i0; i < i1; i++){
	y[i]	=log10(fabs(lL[i])+1e-20);
	y0[i]	=log10(fabs(lR[i])+1e-20);
      }
      Plot2d(5,ESsa+i0,y0+i0,i1-i0,6,0,12,0);
      Plot2d(5,ESsa+i0, y+i0,i1-i0,6,0,14,0);
      if(Fl){
	PSPlot2d(ESsa+i0,y0+i0,i1-i0,6,0,12,0);
	PSPlot2d(ESsa+i0, y+i0,i1-i0,6,0,14,0);
      }
      ki	+=ESNa1;
    }
    CbFlush();
    if(Fl){
      PSFileClose();
      ESErrorInGSh();
    }
  }
  return(0);
}

#endif

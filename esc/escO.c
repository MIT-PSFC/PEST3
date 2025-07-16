#include "esc_local.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

extern clock_t clock(void);
extern time_t time(time_t *t);
extern struct tm *localtime(const time_t *timep);
static struct tm *ltm;
static time_t stTime;

extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern int ESShMem;
extern char ESMessage[];
extern char ShotNm[];
extern int ESNa,ESNa1;
extern int ESNp,ESNp1,ESFp,ESFp1,ESnAP,ESnAF,ESMp,ESMp1,ES2Mp,ES2Mp1;
extern int ESEqSolvFl,ESEqSolvInPr,ESEqSolvInCr;
extern double *ESsa,*ESpa,*ESgt;
extern double ESaRx,ESaZx,*ESaR0,*ESaZ0,*ESsb,*ESsb1a,*ESsb2a;
extern double *rcT,*rcT1a,*rcT2a,*rsT,*rsT1a,*rsT2a;
extern double ESBt,ESRBt,ESRext,ESgbext;
extern double *ESLc,*ESLc1,*ESLc2,*ESLs,*ESLs1,*ESLs2
,*ESVc,*ESVc2,*ESVs,*ESVs2,*ESDPc,*ESDPc1,*ESDPc2;
extern double *ESg22c,*ESg22c1,*ESg22c2,*ESg22s,*ESg22s1,*ESg22s2;
extern double *EScs1,*ESsn1,*EZcs2,*EZsn2,*EZdra,*EZdza,*EZdrgt,*EZdzgt;

extern double *ESsr,*ESsra,*ESsrt,*ESsrat;
extern double *ESsz,*ESsza,*ESszt,*ESszat;
extern double *ESaB,*ESaBa,*ESaBt,*ESaBat;
extern double *ESgH,*ESgHa,*ESgHt,*ESgHat;

extern double *ESgY0;
extern double ESTime; 
extern char ESRunID[];

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

static int FlVm6=1;
static double *rdPV,*zdPV;
int VmNa=20,VmMpol=10;
static int VmNa1=21,FlVm=0,VmNp=64,VmNp1=65;
static double R[65*201],Z[65*201]
,R1[65*201],Z1[65*201],gL1[65*201]
,R2[65*201],Z2[65*201]
,R3[65*201],Z3[65*201]
,R4[4*128*256],Z4[4*128*256],B4[4*128*256],H4[4*128*256],G4[4*128*256];
static clock_t Vmt0;

static int nAjso=64,nA1jso=65,nPjso=64,nP1jso=65,nAPjso=65*65;
static double Rjso[65*65],Zjso[65*65],psibar0;

#ifdef H
static int nAjso=4*64,nA1jso=4*64+1,nPjso=4*64,nP1jso=4*64+1,nAPjso=257*257;
static double Rjso[257*257],Zjso[257*257],psibar0;
#endif

static int FlStart=0;
static int nDjso=21;
static double pprime_s[21],ajb[21];
static char Message[256];

double EZxahg[4000],EZgpahg[4000],EZgp2a[4000];/*data for A.H.Glasser*/
int EZnahg; 			/*data for A.H.Glasser*/

static int *iAMap=NULL;
static int nrMap=20,nzMap=40;
static double rMap[4]={0.,0.,0.,0.},zMap[4]={0.,0.,0.,0.};

#ifdef XWIN
int Dcon2Esc(dfname)
     char dfname[80];
{
  int i,i0,i1;
  double EZxgt[65],EZygt[65],s;
  double x[2],y[2];
  static int mm[2],n1,nt1;
  double *f;

  openftn_(dfname);
  i	=2;
  frdi_(mm,&i);
  n1	=mm[0];
  nt1	=mm[1];
  printf("nr=%4d ntau=%4d\n",n1,nt1);

  printf("Reading Flux coordinates R,Z, [m]\n");
  x[0]	=rdPV[2];
  x[1]	=rdPV[3];
  y[0]	=zdPV[1];
  y[1]	=zdPV[0];
  Scale2D(4,2,x,y,2,6,2);
  for(i=0;i < n1;i++){
    printf("i=%4d x=%12.4e\n",i+1,EZxahg[i]);
    frdd_(EZxgt,&nt1);
    frdd_(EZygt,&nt1);
    Plot2d(4,EZxgt,EZygt,nt1,6,0,0,0);
  }
  CbFlush();
  printf("??? ->");
  getchar();

  f	=(double*)malloc(2*n1*sizeof(double));
  printf("Reading \\Psi, [Vsec]\n");
  frdd_(f,&n1);

  i0=0;
  i1=n1;

#ifdef H
  s	=1./f[n1-1];
#endif
  s	=1./(n1-1);
  for(i=0;i<n1;i++){
    printf("Psi[%3d]=%12.4e\n",i,f[i]);
#ifdef H
    EZxahg[i]=f[i]*s;
#endif
    EZxahg[i]=s*i;
  }
  x[0]	=EZxahg[i0];
  x[1]	=EZxahg[i1-1];
  y[0]	=f[i0];
  y[1]	=f[i1-1];
  Scale2D(4,2,x,y,2,6,2);
  Plot2d(4,EZxahg+i0,f+i0,i1-i0,6,0,0,0);
  CbFlush();
  printf("??? ->");
  getchar();

  printf("Reading q \n");
  frdd_(f,&n1);
  for(i= 0;i<n1;i++)
    printf("q[%3d]=%12.4e\n",i,f[i]);

  x[0]	=EZxahg[i0];
  x[1]	=EZxahg[i1-1];
  y[0]	=f[i0];
  y[1]	=f[i1-1];
  Scale2D(4,2,x,y,2,6,2);
  Plot2d(4,EZxahg+i0,f+i0,i1-i0,6,0,0,0);
  CbFlush();
  printf("??? ->");
  getchar();

  printf("Writing RB_{tor} [mTesla]\n");
  frdd_(f,&n1);
  for(i= 0;i<n1;i++)
    printf("rB[%3d]=%12.4e\n",i,f[i]);

  x[0]	=EZxahg[i0];
  x[1]	=EZxahg[i1-1];
  y[0]	=f[i0];
  y[1]	=f[i1-1];
  Scale2D(4,2,x,y,2,6,2);
  Plot2d(4,EZxahg+i0,f+i0,i1-i0,6,0,0,0);
  CbFlush();
  printf("??? ->");
  getchar();

  printf("Writing {d(RB_{tor})\\over d\\Psi} [mTesla/Vsec]\n");
  frdd_(f,&n1);
  for(i= 0;i<n1;i++)
    printf("FF'(Psi)[%3d]=%12.4e\n",i,f[i]);

  x[0]	=EZxahg[i0];
  x[1]	=EZxahg[i1-1];
  y[0]	=f[i0];
  y[1]	=f[i1-1];
  Scale2D(4,2,x,y,2,6,2);
  Plot2d(4,EZxahg+i0,f+i0,i1-i0,6,0,0,0);
  CbFlush();
  printf("??? ->");
  getchar();

  printf("Writing {dp\\over d\\Psi} [MPa/Vsec]\n");
  frdd_(f,&n1);
  for(i= 0;i<n1;i++)
    printf("p'(Psi)[%3d]=%12.4e\n",i,f[i]);

  x[0]	=EZxahg[i0];
  x[1]	=EZxahg[i1-1];
  y[0]	=f[i0];
  y[1]	=f[i1-1];
  Scale2D(4,2,x,y,2,6,2);
  Plot2d(4,EZxahg+i0,f+i0,i1-i0,6,0,0,0);
  CbFlush();
  printf("??? ->");
  getchar();

  closeftn_();
  free(f);
  printf("%s Interface file for DCON\n",dfname);
  return(0);
}
#endif

int ES4Dcon1()
{
  FILE *fp;
  int i,ii,j,k,ki;
  double a,b,rc[65],rs[65],w[10];
  int na,na1;
  double dx;

  fp=fopen("dc.d","w");

  na	=1024;
  na1	=na+1;
  dx	=ESsa[ESNa]/na;
  spl(EZgpahg,EZgp2a,ESsa,ESdgY+ESNa,EZxahg,&EZnahg);
  ii	=0;
  for(i=0; i < na1; i++){
    a		=dx*i;
    while(EZxahg[ii+1] < a){
      ii++;
    }
    splr1(w,NULL,&a,EZgpahg,EZgp2a,EZxahg,&EZnahg,&ii); /*$\Psi$*/
    ESSetSplA(a);
    splRA(w+1,NULL,ESgm,ESgm2a);
    w[1]	=1./w[1];	/*$q$*/
    splRA(w+2,NULL,ESaF,ESaF2a);/*$RB_{tor}$*/
    splRA(w+3,NULL,ESaT,ESaT2a);
    w[3]	/=w[2];		/*${d(RB_{tor}\over d\Psi}$*/
    splRA(w+4,NULL,ESPs,ESPs2a);		/*${dp\over d\Psi}$*/
    splRA(w+5,NULL,ESsb,ESsb2a);/*$b$*/
    splRA(rc,NULL,rcT,rcT2a);
    splRA(rs,NULL,rsT,rsT2a);
    ki		=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,NULL,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      splRA(rs+k,NULL,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
    }
    fwrite(w,(size_t)sizeof(double),(size_t)6,fp);
    fwrite(rc,(size_t)sizeof(double),(size_t)33,fp);
    fwrite(rs,(size_t)sizeof(double),(size_t)33,fp);
  }
  fclose(fp);
  return(0);
}

#ifdef XWIN
int ES4Dcon(double ttime)
{
  int i,ii,j,k,ki,kj;
  double a,b,rc[65],rs[65],EZxgt[257],EZygt[257],csD[257],snD[257];
  int mm[2],n1,na,na1,np,np1;
  double *f,*xa,dx,t,dt,s,ss,ds;
  char str[12];

  sprintf(str,"ES4Dcon%4.2f",ttime);
  sprintf(str,"ES4DconD");

  na	=128;
  na	=256;
  np	=64;
  na1	=na+1;
  np1	=np+1;
  openftn_(str);
  data4ahg_(EZxahg,EZgpahg,&EZnahg);
  if(EZnahg < 4000){
    printf("EZnahg=%d\n",EZnahg);
  }
  else{
    printf("EZnahg=%d%c\n",EZnahg,7);
    return(0);
  }

  n1	=EZnahg+1;
  i	=EZnahg;
  j	=i-1;
  while(i > 0){
    EZxahg[i]	=EZxahg[j];
    EZgpahg[i]	=EZgpahg[j];
#ifdef H
    printf("i=%4d x=%12.4e %12.4e\n",i,EZxahg[i],EZgpahg[i]);
#endif
    i--;
    j--;
  }
  EZxahg[i]	=0.;
  EZgpahg[i]	=0.;

#ifdef H
  printf("i=%4d x=%12.4e %12.4e\n",i,EZxahg[i],EZgpahg[i]);
#endif

  f	=(double*)malloc(2*na1*sizeof(double));
  xa	=f+na1;
  dx	=ESsa[ESNa]/na;

  spl(EZgpahg,EZgp2a,ESsa,ESdgY+ESNa,EZxahg,&EZnahg);
  j	=2;
#ifdef H
  printf("nr=%4d ntau=%4d\n",ESNa1,ESNp1);
#endif
  mm[0]	=na1;
  mm[1]	=np1;
  fwri_(mm,&j);		/* {\tt write(1)ma,mtau}*/
  
  printf("Writing Flux coordinates R,Z, [m]\n");
  for(j=0; j < np1; j++){
    EZxgt[j]	=ESaR0[0];
    EZygt[j]	=ESaZ0[0];
  }
  fwrd_(EZxgt,&np1);/*$r$, {\tt write(1)(r(itau,ia)),itau=1,mtau)}*/
  fwrd_(EZygt,&np1);/*$r$, {\tt write(1)(z(itau,ia)),itau=1,mtau)}*/
  xa[0]	=0.;
  dt	=EZc2gp/np;
  for(j=0; j < np1; j++){
    t=dt*j;
    csD[j]	=cos(t);
    snD[j]	=sin(t);
  }
  a	=0.;
  for(i=1; i < na1; i++){
    ss	=sqrt(dx*i);
    do{
      ESSetSplA(a);
      splRA(&s,&ds,ESqgY,ESqgY2a);
      a		+=(ss-s)/ds;
    }while(fabs(ss-s) > 1e-8);
    xa[i]	=a;
    ESSetSplA(a);
    splRA(&b,NULL,ESsb,ESsb2a);
    splRA(rc,NULL,rcT,rcT2a);
    splRA(rs,NULL,rsT,rsT2a);
    ki		=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,NULL,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      splRA(rs+k,NULL,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
    }
    for(j=0; j < np1; j++){
      t	=dt*j;
      EZxgt[j]	=ESaR0[0]+rc[0];
      EZygt[j]	=ESaZ0[0]+rs[0]+b*snD[j];
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= np)
	  kj	-=np;
	EZxgt[j]	+=rc[k]*csD[kj]+rs[k]*snD[kj];
      }
    }
    fwrd_(EZxgt,&np1);/*$r$, {\tt write(1)(r(itau,ia)),itau=1,mtau)}*/
    fwrd_(EZygt,&np1);/*$r$, {\tt write(1)(z(itau,ia)),itau=1,mtau)}*/
  }

  printf("Writing \\Psi, [Vsec]\n");
  ii	=0;
  for(i=0; i < na1; i++){
    while(ii < EZnahg && EZxahg[ii+1] < xa[i]){
      ii++;
    }
    splr1(f+i,NULL,xa+i,EZgpahg,EZgp2a,EZxahg,&EZnahg,&ii);
    f[i]	*=EZc2gp;
#ifdef H
    printf("x[%4d]=%12.4e f=%12.4e\n",i,xa[i],f[i]);
#endif
  }
  fwrd_(f,&na1);		/*$\Psi,\ [Vsec]$*/

  printf("Writing q \n");
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(&b,NULL,ESgm,ESgm2a);
    f[i]	=1./b;
  }
  fwrd_(f,&na1);		/*$q$*/

  printf("Writing RB_{tor} [mTesla]\n");
  f[0]	=sqrt(ESaF[0]);
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(f+i,NULL,ESaF,ESaF2a);
  }
  fwrd_(f,&na1);		/*$RB_{tor}$*/

  printf("Writing {d(RB_{tor})\\over d\\Psi} [mTesla/Vsec]\n");
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(&b,NULL,ESaT,ESaT2a);
    f[i]=b*EZcr2gp/f[i];
  }
  fwrd_(f,&na1);		/*${d(RB_{tor}\over d\Psi}$*/

  printf("Writing {dp\\over d\\Psi} [MPa/Vsec]\n");
  b	=5.*EZcr2gp*EZcr2gp;
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(f+i,NULL,ESPs,ESPs2a);
    f[i]	*=b;
  }
  fwrd_(f,&na1);		/*${dp\over d\Psi}$*/
  closeftn_();
  free(f);
  printf("%s Interface file for DCON\n",str);
#ifdef H
  Dcon2Esc(str);
#endif
  return(0);
}

int ES4DconFloat0(double ttime)
{
  int i,ii,j,k,ki,kj;
  double a,b,rc[65],rs[65],EZxgt[257],EZygt[257],csD[257],snD[257];
  int mm[2],n1,na,na1,np,np1;
  double dx,t,dt,s,ss,ds,*xa,*f;
  float xf[257],yf[257],*ff;
  char str[12];

  double ESgY02a[ESNa1];

  sprintf(str,"ES4Dcon%4.2f",ttime);

  na	=128;
  na	=256;
  np	=64;
  na1	=na+1;
  np1	=np+1;
  openftn_(str);
  data4ahg_(EZxahg,EZgpahg,&EZnahg);
  if(EZnahg < 4000){
    printf("EZnahg=%d\n",EZnahg);
  }
  else{
    printf("EZnahg=%d%c\n",EZnahg,7);
    return(0);
  }

  n1	=EZnahg+1;
  i	=EZnahg;
  j	=i-1;
  while(i > 0){
    EZxahg[i]	=EZxahg[j];
    EZgpahg[i]	=EZgpahg[j];
#ifdef H
    printf("i=%4d x=%12.4e %12.4e\n",i,EZxahg[i],EZgpahg[i]);
#endif
    i--;
    j--;
  }
  EZxahg[i]	=0.;
  EZgpahg[i]	=0.;

#ifdef H
  printf("i=%4d x=%12.4e %12.4e\n",i,EZxahg[i],EZgpahg[i]);
#endif

  f	=(double*)malloc(2*na1*sizeof(double));
  xa	=f+na1;
  dx	=ESsa[ESNa]/na;
  ff	=(float*)malloc(na1*sizeof(float));

  spl(EZgpahg,EZgp2a,ESsa,ESdgY+ESNa,EZxahg,&EZnahg);
  spl(ESgY0,ESgY02a,ESsa,ESdgY+ESNa,ESsa,&ESNa);

  j	=2;
#ifdef H
  printf("nr=%4d ntau=%4d\n",ESNa1,ESNp1);
#endif
  mm[0]	=na1;
  mm[1]	=np1;
  fwri_(mm,&j);		/* {\tt write(1)ma,mtau}*/
  
  printf("Writing Flux coordinates R,Z, [m]\n");
  for(j=0; j < np1; j++){
    EZxgt[j]	=ESaR0[0];
    EZygt[j]	=ESaZ0[0];
    xf[j]	=EZxgt[j];
    yf[j]	=EZygt[j];
  }
  fwrf_(xf,&np1);/*$r$, {\tt write(1)(r(itau,ia)),itau=1,mtau)}*/
  fwrf_(yf,&np1);/*$r$, {\tt write(1)(z(itau,ia)),itau=1,mtau)}*/
  xa[0]	=0.;
  dt	=EZc2gp/np;
  for(j=0; j < np1; j++){
    t=dt*j;
    csD[j]	=cos(t);
    snD[j]	=sin(t);
  }
  a	=0.;
  for(i=1; i < na1; i++){
    ss	=sqrt(dx*i);
    do{
      ESSetSplA(a);
      splRA(&s,&ds,ESqgY,ESqgY2a);
      a		+=(ss-s)/ds;
    }while(fabs(ss-s) > 1e-8);
    xa[i]	=a;
    ESSetSplA(a);
    splRA(&b,NULL,ESsb,ESsb2a);
    splRA(rc,NULL,rcT,rcT2a);
    splRA(rs,NULL,rsT,rsT2a);
    ki		=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,NULL,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      splRA(rs+k,NULL,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
    }
    for(j=0; j < np1; j++){
      t	=dt*j;
      EZxgt[j]	=ESaR0[0]+rc[0];
      EZygt[j]	=ESaZ0[0]+rs[0]+b*snD[j];
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= np)
	  kj	-=np;
	EZxgt[j]	+=rc[k]*csD[kj]+rs[k]*snD[kj];
      }
      xf[j]	=EZxgt[j];
      yf[j]	=EZygt[j];
    }
    fwrf_(xf,&np1);/*$r$, {\tt write(1)(r(itau,ia)),itau=1,mtau)}*/
    fwrf_(yf,&np1);/*$r$, {\tt write(1)(z(itau,ia)),itau=1,mtau)}*/
  }

  printf("Writing \\Psi, [Vsec]\n");
  ii	=0;
  for(i=0; i < na1; i++){
    while(ii < EZnahg && EZxahg[ii+1] < xa[i]){
      ii++;
    }
    splr1(f+i,NULL,xa+i,EZgpahg,EZgp2a,EZxahg,&EZnahg,&ii);
    f[i]	*=EZc2gp;
    ff[i]	=f[i];
    ESSetSplA(xa[i]);
    splRA(f+i,NULL,ESgY0,ESgY02a);
    f[i]	*=EZc2gp;
  }
  fwrf_(ff,&na1);		/*$\Psi,\ [Vsec]$*/

  printf("Writing q \n");
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(&b,NULL,ESgm,ESgm2a);
    f[i]	=1./b;
    ff[i]	=f[i];
  }
  fwrf_(ff,&na1);		/*$q$*/

  printf("Writing RB_{tor} [mTesla]\n");
  f[0]	=sqrt(ESaF[0]);
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(f+i,NULL,ESaF,ESaF2a);
    ff[i]	=f[i];
  }
  fwrf_(ff,&na1);		/*$RB_{tor}$*/

  printf("Writing {d(RB_{tor})\\over d\\Psi} [mTesla/Vsec]\n");
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(&b,NULL,ESaT,ESaT2a);
    f[i]=b*EZcr2gp/f[i];
    ff[i]	=f[i];
  }
  fwrf_(ff,&na1);		/*${d(RB_{tor}\over d\Psi}$*/

  printf("Writing {dp\\over d\\Psi} [MPa/Vsec]\n");
  b	=5.*EZcr2gp*EZcr2gp;
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(f+i,NULL,ESPs,ESPs2a);
    f[i]	*=b;
    ff[i]	=f[i];
  }
  fwrf_(ff,&na1);		/*${dp\over d\Psi}$*/
  closeftn_();
  free(ff);
  free(f);
  printf("%s Interface file for DCON\n",str);
#ifdef H
  Dcon2Esc(str);
#endif
  return(0);
}
#endif

int ES4DconFloat(double ttime)
{
  int i,ii,j,k,ki,kj;
  int na=256,na1=na+1,np=64,np1=np+1;
  int mm[2];
  char str[16];
#ifdef XWIN
  double a,b,rc[65],rs[65],EZxgt[np1],EZygt[np1],csD[np1],snD[np1];
  double dx,t,dt,s,ss,ds,xa[na1],f[na1];
  float xf[np1],yf[np1],ff[na1];
  double ESgY02a[ESNa1];
#else
  double a,b,rc[65],rs[65],EZxgt[129],EZygt[129],csD[129],snD[129];
  double dx,t,dt,s,ss,ds,xa[257],f[257];
  float xf[129],yf[129],ff[257];
  double ESgY02a[129];
#endif

  sprintf(str,"ES4Dcon%4.2f",ttime);
  openftn_(str);

  dx	=ESsa[ESNa]/na;
  spl(ESgY0,ESgY02a,ESsa,ESdgY+ESNa,ESsa,&ESNa);

  j	=2;
  mm[0]	=na1;
  mm[1]	=np1;
  fwri_(mm,&j);		/* {\tt write(1)ma,mtau}*/
  
  printf("Writing Flux coordinates R,Z, [m]\n");
  for(j=0; j < np1; j++){
    EZxgt[j]	=ESaR0[0];
    EZygt[j]	=ESaZ0[0];
    xf[j]	=EZxgt[j];
    yf[j]	=EZygt[j];
  }
  fwrf_(xf,&np1);/*$r$, {\tt write(1)(r(itau,ia)),itau=1,mtau)}*/
  fwrf_(yf,&np1);/*$r$, {\tt write(1)(z(itau,ia)),itau=1,mtau)}*/
  xa[0]	=0.;
  dt	=EZc2gp/np;
  for(j=0; j < np1; j++){
    t	=dt*j;
    csD[j]	=cos(t);
    snD[j]	=sin(t);
  }
  a	=sqrt(dx);
  for(i=1; i < na1; i++){
#ifdef H
    ss	=dx*i*ESgY0[ESNa];
    do{
      ESSetSplA(a);
      splRA(&s,&ds,ESgY0,ESgY02a);
      ds	=(ss-s)/ds;
      a		+=ds;
    }while(fabs(ds) > 1e-8);
#endif
    a	=dx*i;

    xa[i]	=a;
    ESSetSplA(a);
    splRA(&b,NULL,ESsb,ESsb2a);
    splRA(rc,NULL,rcT,rcT2a);
    splRA(rs,NULL,rsT,rsT2a);
    ki		=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,NULL,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      splRA(rs+k,NULL,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
    }
    for(j=0; j < np1; j++){
      t	=dt*j;
      EZxgt[j]	=ESaR0[0]+rc[0];
      EZygt[j]	=ESaZ0[0]+rs[0]+b*snD[j];
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= np)
	  kj	-=np;
	EZxgt[j]	+=rc[k]*csD[kj]+rs[k]*snD[kj];
      }
      xf[j]	=EZxgt[j];
      yf[j]	=EZygt[j];
    }
    fwrf_(xf,&np1);/*$r$, {\tt write(1)(r(itau,ia)),itau=1,mtau)}*/
    fwrf_(yf,&np1);/*$r$, {\tt write(1)(z(itau,ia)),itau=1,mtau)}*/
  }

  printf("Writing \\Psi, [Vsec]\n");
  ff[0]	=0.;
  for(i=1; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(&s,NULL,ESgY0,ESgY02a);
    ff[i]	=s*EZc2gp;
  }
  fwrf_(ff,&na1);		/*$\Psi,\ [Vsec]$*/

  printf("Writing q \n");
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(&s,NULL,ESgm,ESgm2a);
    ff[i]	=1./s;
  }
  fwrf_(ff,&na1);		/*$q$*/

  printf("Writing RB_{tor} [mTesla]\n");
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(&s,NULL,ESaF,ESaF2a);
    f[i]	=s;
    ff[i]	=s;
  }
  fwrf_(ff,&na1);		/*$RB_{tor}$*/

  printf("Writing {d(RB_{tor})\\over d\\Psi} [mTesla/Vsec]\n");
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(&b,NULL,ESaT,ESaT2a);
    f[i]	=b*EZcr2gp/f[i];
    ff[i]	=f[i];
  }
  fwrf_(ff,&na1);		/*${d(RB_{tor}\over d\Psi}$*/

  printf("Writing {dp\\over d\\Psi} [MPa/Vsec]\n");
  b	=5.*EZcr2gp*EZcr2gp;
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(f+i,NULL,ESPs,ESPs2a);
    f[i]	*=b;
    ff[i]	=f[i];
  }
  fwrf_(ff,&na1);		/*${dp\over d\Psi}$*/
  closeftn_();
  printf("%s Interface file for DCON\n",str);
#ifdef H
  Dcon2Esc(str);
#endif
  return(0);
}

#ifdef XWIN
int ES4DconASCII(double ttime)
{
  int i,ii,j,k,ki,kj;
  double a,b,rc[65],rs[65],EZxgt[257],EZygt[257],csD[257],snD[257];
  int mm[2],n1,na,na1,np,np1;
  double *f,*xa,dx,t,dt,s,ss,ds;
  char str[16];
  FILE *fp;

  sprintf(str,"ES4Dcon%4.2fA",ttime);
  sprintf(str,"ES4DconA");

  na	=128;
  na	=256;
  np	=64;
  na1	=na+1;
  np1	=np+1;
  fp=fopen(str,"w");
  data4ahg_(EZxahg,EZgpahg,&EZnahg);
  if(EZnahg < 4000){
    printf("EZnahg=%d\n",EZnahg);
  }
  else{
    printf("EZnahg=%d%c\n",EZnahg,7);
    return(0);
  }

  n1	=EZnahg+1;
  i	=EZnahg;
  j	=i-1;
  while(i > 0){
    EZxahg[i]	=EZxahg[j];
    EZgpahg[i]	=EZgpahg[j];
    i--;
    j--;
  }
  EZxahg[i]	=0.;
  EZgpahg[i]	=0.;

  f	=(double*)malloc(2*na1*sizeof(double));
  xa	=f+na1;
  dx	=ESsa[ESNa]/na;

  spl(EZgpahg,EZgp2a,ESsa,ESdgY+ESNa,EZxahg,&EZnahg);
  j	=2;
  fprintf(fp,"%4d %4d\n",na1,np1);
  
  printf("Writing Flux coordinates R,Z, [m]\n");
  for(j=0; j < np1; j++){
    EZxgt[j]	=ESaR0[0];
    EZygt[j]	=ESaZ0[0];
    fprintf(fp,"%22.15e %22.15e\n",EZxgt[j],EZygt[j]);
  }
  xa[0]	=0.;
  dt	=EZc2gp/np;
  for(j=0; j < np1; j++){
    t=dt*j;
    csD[j]	=cos(t);
    snD[j]	=sin(t);
  }
  a	=0.;
  for(i=1; i < na1; i++){
    ss	=sqrt(dx*i);
    do{
      ESSetSplA(a);
      splRA(&s,&ds,ESqgY,ESqgY2a);
      a		+=(ss-s)/ds;
    }while(fabs(ss-s) > 1e-8);
    xa[i]	=a;
    ESSetSplA(a);
    splRA(&b,NULL,ESsb,ESsb2a);
    splRA(rc,NULL,rcT,rcT2a);
    splRA(rs,NULL,rsT,rsT2a);
    ki		=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,NULL,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      splRA(rs+k,NULL,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
    }
    for(j=0; j < np1; j++){
      t	=dt*j;
      EZxgt[j]	=ESaR0[0]+rc[0];
      EZygt[j]	=ESaZ0[0]+rs[0]+b*snD[j];
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= np)
	  kj	-=np;
	EZxgt[j]	+=rc[k]*csD[kj]+rs[k]*snD[kj];
      }
      fprintf(fp,"%22.15e %22.15e\n",EZxgt[j],EZygt[j]);
    }
  }
  printf("Writing \\Psi, [Vsec]\n");
  ii	=0;
  for(i=0; i < na1; i++){
    while(ii < EZnahg && EZxahg[ii+1] < xa[i]){
      ii++;
    }
    splr1(f+i,NULL,xa+i,EZgpahg,EZgp2a,EZxahg,&EZnahg,&ii);
    f[i]	*=EZc2gp;
    fprintf(fp,"%22.15e\n",f[i]);
  }

  printf("Writing q \n");
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(&b,NULL,ESgm,ESgm2a);
    f[i]	=1./b;
    fprintf(fp,"%22.15e\n",f[i]);
  }

  printf("Writing RB_{tor} [mTesla]\n");
  f[0]	=sqrt(ESaF[0]);
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(f+i,NULL,ESaF,ESaF2a);
    fprintf(fp,"%22.15e\n",f[i]);
  }

  printf("Writing {d(RB_{tor})\\over d\\Psi} [mTesla/Vsec]\n");
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(&b,NULL,ESaT,ESaT2a);
    f[i]=b*EZcr2gp/f[i];
    fprintf(fp,"%22.15e\n",f[i]);
  }

  printf("Writing {dp\\over d\\Psi} [MPa/Vsec]\n");
  b	=5.*EZcr2gp*EZcr2gp;
  for(i=0; i < na1; i++){
    ESSetSplA(xa[i]);
    splRA(f+i,NULL,ESPs,ESPs2a);
    f[i]	*=b;
    fprintf(fp,"%22.15e\n",f[i]);
  }
  fclose(fp);
  free(f);
  printf("%s Interface file for DCON\n",str);
  return(0);
}
#endif

int ES4DconF()
{
  int i,ii,j,k,ki,kj;
  double a,b,rc[65],rs[65],csD[257],snD[257];
  float EZxgt[257],EZygt[257],f[257];
  int mm[2],n1,na,na1,np,np1;
  double F,xa[257],dx,t,dt;
  char str[24];
  int ik,ka[7]={64,128,64,128,128,256,256},kp[7]={64,64,128,128,256,128,256};

#ifdef H
  data4ahg_(EZxahg,EZgpahg,&EZnahg);
  if(EZnahg < 4000){
    printf("EZnahg=%d\n",EZnahg);
  }
  else{
    printf("EZnahg=%d%c\n",EZnahg,7);
    return(0);
  }
  
  n1	=EZnahg+1;
  i	=EZnahg;
  j	=i-1;
  while(i > 0){
    EZxahg[i]	=EZxahg[j];
    EZgpahg[i]	=EZgpahg[j];
    i--;
    j--;
  }
  EZxahg[i]	=0.;
  EZgpahg[i]	=0.;
  spl(EZgpahg,EZgp2a,ESsa,ESdgY+ESNa,EZxahg,&EZnahg);
#endif
  for(ik=0; ik < 7; ik++){
    na	=ka[ik];
    np	=kp[ik];
    na1	=na+1;
    np1	=np+1;
    sprintf(str,"ES4Dcon%dx%d",na,np);
    openftn_(str);

    dx	=ESsa[ESNa]/na;
    j	=2;
    mm[0]	=na1;
    mm[1]	=np1;
    fwri_(mm,&j);		/* {\tt write(1)ma,mtau}*/
  
    printf("Writing Flux coordinates R,Z, [m]\n");
    for(j=0; j < np1; j++){
      EZxgt[j]	=ESaR0[0];
      EZygt[j]	=ESaZ0[0];
    }
    fwrf_(EZxgt,&np1);/*$r$, {\tt write(1)(r(itau,ia)),itau=1,mtau)}*/
    fwrf_(EZygt,&np1);/*$r$, {\tt write(1)(z(itau,ia)),itau=1,mtau)}*/
    xa[0]	=0.;
    dt	=EZc2gp/np;
    for(j=0; j < np1; j++){
      t=dt*j;
      csD[j]	=cos(t);
      snD[j]	=sin(t);
    }
    for(i=1; i < na1; i++){
      xa[i]	=dx*i;
      a		=xa[i];
      ESSetSplA(a);
      splRA(&b,NULL,ESsb,ESsb2a);
      splRA(rc,NULL,rcT,rcT2a);
      splRA(rs,NULL,rsT,rsT2a);
      ki		=0;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	splRA(rc+k,NULL,rcT+ki,rcT2a+ki);
	rc[k]	*=2.;
	splRA(rs+k,NULL,rsT+ki,rsT2a+ki);
	rs[k]	*=2.;
      }
      for(j=0; j < np1; j++){
	t	=dt*j;
	EZxgt[j]	=ESaR0[0]+rc[0];
#ifdef H
	EZygt[j]	=ESaZ0[0]+rs[0]+b*ESsn1[j];
#endif
	EZygt[j]	=ESaZ0[0]+rs[0]+b*sin(t);
	kj	=0;
	for(k=1; k < ESFp1; k++){
	  kj	+=j;
	  if(kj >= np)
	    kj	-=np;
#ifdef H
	  EZxgt[j]	+=rc[k]*EScs1[kj]+rs[k]*ESsn1[kj];
#endif
	  EZxgt[j]	+=rc[k]*cos(k*t)+rs[k]*sin(k*t);
	}
      }
      fwrf_(EZxgt,&np1);	/*$r$, {\tt write(1)(r(itau,ia)),itau=1,mtau)}*/
      fwrf_(EZygt,&np1);	/*$r$, {\tt write(1)(z(itau,ia)),itau=1,mtau)}*/
    }
    
    printf("Writing \\Psi, [Vsec]\n");
    ii	=0;
    for(i=0; i < na1; i++){
      while(EZxahg[ii+1] < xa[i]){
	ii++;
      }
      splr1(&F,NULL,xa+i,EZgpahg,EZgp2a,EZxahg,&EZnahg,&ii);
      f[i]	=F*EZc2gp;
#ifdef H
      printf("x[%4d]=%12.4e f=%12.4e\n",i,xa[i],f[i]);
#endif
    }
    fwrf_(f,&na1);		/*$\Psi,\ [Vsec]$*/
    
    printf("Writing q \n");
    for(i=0; i < na1; i++){
      ESSetSplA(xa[i]);
      splRA(&b,NULL,ESgm,ESgm2a);
      f[i]	=1./b;
    }
    fwrf_(f,&na1);		/*$q$*/
    
    printf("Writing RB_{tor} [mTesla]\n");
    f[0]	=sqrt(ESaF[0]);
    for(i=0; i < na1; i++){
      ESSetSplA(xa[i]);
      splRA(&F,NULL,ESaF,ESaF2a);
      f[i]	=F;
    }
    fwrf_(f,&na1);		/*$RB_{tor}$*/
    
    printf("Writing {d(RB_{tor})\\over d\\Psi} [mTesla/Vsec]\n");
    for(i=0; i < na1; i++){
      ESSetSplA(xa[i]);
      splRA(&b,NULL,ESaT,ESaT2a);
      f[i]=b*EZcr2gp/f[i];
    }
    fwrf_(f,&na1);		/*${d(RB_{tor}\over d\Psi}$*/
    
    printf("Writing {dp\\over d\\Psi} [MPa/Vsec]\n");
    b	=5.*EZcr2gp*EZcr2gp;
    for(i=0; i < na1; i++){
      ESSetSplA(xa[i]);
      splRA(&F,NULL,ESPs,ESPs2a);
      f[i]	=F*b;
    }
    fwrf_(f,&na1);		/*${dp\over d\Psi}$*/
    closeftn_();
  }
  printf("%s Interface file for DCON\n",str);
  return(0);
}

int ES4Mars(double ttime)
{
  static double t,gp,df,q,F,Fs,Ps,dt;
  int i,j,n,n1;
  char str[12];
  FILE *fpnt;

  sprintf(str,"ES4Mars%4.2f",ttime);
  fpnt	=fopen(str,"w");
  printf("ES4Mars%4.2f file for MARS",ttime);

  n	=64;
  n1	=n+1;
  fprintf(fpnt,"%5.5s %7.4f %7.5f %12.4e\n",ESRunID,ttime
	  ,ESBt,ESaJ[ESNa]*5e+6);

  fprintf(fpnt,"%3d\n",n1);
  dt	=ESsa[ESNa]/n;
  fprintf(fpnt,"PSI       \n");
  for(i=0; i < n1; i++){
    t	=dt*i;
    ESSetSplA(t);
    splRA(&gp,NULL,ESgY,ESgY2a);
    fprintf(fpnt,"%16.9e\n",gp);
  }
  fprintf(fpnt,"F         \n");
  for(i=0; i < n1; i++){
    t	=dt*i;
    ESSetSplA(t);
    splRA(&F,NULL,ESFF,ESFF2a);
    F	=sqrt(F);
    fprintf(fpnt,"%16.9e\n",F);
  }
  fprintf(fpnt,"FdF/dPSI      \n");
  for(i=0; i < n1; i++){
    t	=dt*i;
    ESSetSplA(t);
    splRA(&Fs,NULL,ESaT,ESaT2a);
#ifdef H
    splRA(&F,NULL,ESFF,ESFF2a);
    Fs	/=sqrt(F);
#endif
    fprintf(fpnt,"%16.9e\n",Fs);
  }
  fprintf(fpnt,"DP/DPSI         \n");
  for(i=0; i < n1; i++){
    t	=dt*i;
    ESSetSplA(t);
    splRA(&Ps,NULL,ESPs,ESPs2a);
    Ps	*=EZcrgm0*1e+6;
    fprintf(fpnt,"%16.9e\n",Ps);
  }

  fprintf(fpnt,"%3d\n",ESNp);
  fprintf(fpnt,"R_BOUNDARY\n");
  i	=ESNp1*ESNa;
  for(j=0; j < ESNp; j++){
    fprintf(fpnt,"%16.9e\n",ESsr[i+j]);
  }
  fprintf(fpnt,"Z_BOUNDARY\n");
  i	=ESNp1*ESNa;
  for(j=0; j < ESNp; j++){
    fprintf(fpnt,"%16.9e\n",ESsz[i+j]);
  }
  fprintf(fpnt,"q         \n");
  for(i=0; i < n1; i++){
    t	=dt*i;
    if(ESEqSolvInCr == 6 || ESEqSolvInCr == 7){
      ESSetSplDCr(t);
      splRDCr(&q,NULL,7);
    }
    else{
      ESSetSplA(t);
      splRA(&q,NULL,ESsq,ESsq2a);
    }
    fprintf(fpnt,"%16.9e\n",q);
  }
  fprintf(fpnt,"beta_ext   \n");
  fprintf(fpnt,"%16.9e\n",ESgbext);
  fclose(fpnt);
  return(0);
}

#ifdef VMEC
int ESInitjso()
{
  int i;
  int ndoub,mn,isym;

  if(FlStart){
    ndoub	=-1;
    mn		=8;
    mn		=6;
  }
  else{
    ndoub	=1;
    mn		=7;
    mn		=5;
  }

  isym	=0;
  jsoinit_(&mn,&ndoub,&isym);

  i	=1; /* 0 => mesh in sqrt(psi)
	       1 => mesh in psi
	       */
  jsosetipsi_(&i);
  
  i	=11; /* max number of inner iterations before
		metric update
		*/
  jsosetnimax_(&i);
      
  i	=1000; /* max number of outer iterations */
  jsosetnumit_(&i);

  i	=0;
  jsosetiverbose_(&i);
  i	=0;
  jsosetieqdsk_(&i);

  return(0);
}

int ESEsc2jso(double *R11,double *Z11)
{
  int i,kmax=12;
  double facimp=-1e-6;
  
  double *pr,*pz;
  double rd[13],zd[13];

  i	=ESNp1*ESNa;
  pr	=ESsr+i;
  pz	=ESsz+i;

  rd[0]	=pr[0];
  zd[0]	=pz[0];
  rd[1]	=pr[6];
  zd[1]	=pz[6];
  rd[2]	=pr[12];
  zd[2]	=pz[12];
  rd[3]	=pr[16];
  zd[3]	=pz[16];
  rd[4]	=pr[20];
  zd[4]	=pz[20];
  rd[5]	=pr[26];
  zd[5]	=pz[26];
  rd[6]	=pr[32];
  zd[6]	=pz[32];
  rd[7]	=pr[38];
  zd[7]	=pz[38];
  rd[8]	=pr[44];
  zd[8]	=pz[44];
  rd[9]	=pr[48];
  zd[9]	=pz[48];
  rd[10]=pr[52];
  zd[10]=pz[52];
  rd[11]=pr[58];
  zd[11]=pz[58];
  rd[12]=pr[64];
  zd[12]=pz[64];

  ESe2jsoPr(pprime_s,ajb);
  ESInitjso();
  jsoexec_(&facimp,&ESRBt,&ESNp,&nDjso,pprime_s,ajb,pr,pz,&psibar0,
	   Rjso,Zjso,&i);
  jsofree_();
  if(i != 0){
    printf("JSOLVER TEST FAILED\n");
    FlStart=0;
  }
  else{
    int i1,j,ji,ji1;
    int ma,mp;
    
    nAjso	=64;
    nA1jso	=nAjso+1;
    nPjso	=64;
    nP1jso	=nPjso+1;
#ifdef H
#endif
    ma	=nA1jso/32;
    mp	=nP1jso/64;

    printf("??? ma=%d mp=%d\n",ma,mp);

    FlStart	=1;
    ji1	=0;
    for(i=0; i < nA1jso; i +=ma){
      ji	=nP1jso*i;
      for(j=0; j < nP1jso; j +=mp){
	R11[ji1]	=Rjso[ji];
	Z11[ji1]	=Zjso[ji];
	ji		+=mp;
	ji1++;
      }
    }
  }
  return(0);
}
#endif

double piota_(double *x)
{
  int i;
  double a,gm;
  
  a	=sqrt(*x);
  if(ESEqSolvInCr == 6 || ESEqSolvInCr == 7){
    ESSetSplDCr(a);
    splRDCr(&gm,NULL,7);
  }
  else{
    ESSetSplA(a);
    splRA(&gm,NULL,ESgm,ESgm2a);
  }
  return(gm);
}

double pmass_(double *x)
{
  int i;
  double a,p;
  
  a	=sqrt(*x);
  if(ESEqSolvInPr == 2){
    ESSetSplDPr(a);
    splRDPr(&p,NULL,2);
  }
  else{
    ESSetSplA(a);
    splRA(&p,NULL,ESsp,ESsp2a);
  }
  return(p);
}

double piotaY_(double *x)
{
  int i;
  double a,gm,ss,ds,s;
  
  ss	=sqrt(*x);
  a	=ss;
  do{
    ESSetSplA(a);
    splRA(&s,&ds,ESqgY,ESqgY2a);
    a		+=(ss-s)/ds;
  }while(fabs(ss-s) > 1e-6);
  if(ESEqSolvInCr == 6 || ESEqSolvInCr == 7){
    ESSetSplDCr(a);
    splRDCr(&gm,NULL,7);
  }
  else{
    ESSetSplA(a);
    splRA(&gm,NULL,ESgm,ESgm2a);
  }
  return(gm);
}

double pmassY_(double *x)
{
  int i;
  double a,p,s,ds,ss;
  
  ss	=sqrt(*x);
  a	=ss;
  do{
    ESSetSplA(a);
    splRA(&s,&ds,ESqgY,ESqgY2a);
    a		+=(ss-s)/ds;
  }while(fabs(ss-s) > 1e-6);
  if(ESEqSolvInPr == 2){
    ESSetSplDPr(a);
    splRDPr(&p,NULL,2);
  }
  else{
    ESSetSplA(a);
    splRA(&p,NULL,ESsp,ESsp2a);
  }
  return(p);
}

int ESe2vm8G(double *rv, double *zv)
{
  int i,j,ji,k,ki,kj;
  double da,ds,a,b,rc[65],rs[65];

  ji	=0;
  for(j=0; j < ESNp1; j++){
    rv[ji]	=ESaR0[0];
    zv[ji]	=ESaZ0[0];
    ji++;
  }

  da	=1./VmNa;
  for(i=1; i < VmNa1; i++){
    a	=sqrt(da*i);
    ESSetSplA(a);
    splRA(&b,&ds,ESsb,ESsb2a);
    splRA(rc,&ds,rcT,rcT2a);
    splRA(rs,&ds,rsT,rsT2a);
    ki		=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,&ds,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      splRA(rs+k,&ds,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
    }
    for(j=0; j < ESNp1; j++){
      rv[ji]	=ESaR0[0]+rc[0];
      zv[ji]	=ESaZ0[0]+rs[0]+b*ESsn1[j];;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	rv[ji]	+=rc[k]*EScs1[kj]+rs[k]*ESsn1[kj];
      }
      ji++;
    }

    splRA(&b,NULL,ESaF,ESaF2a);
    printf("F[%3d]=%12.5e\n",i,b);
  }
  return(0);
}

int ESe2vm8gL(double *gd,double *gd1,double *r,double *z,double *gL,double a)
{
  int i,j,k,ki,kj,jmx,jmn;
  double Lc[33],Ls[33],L0,gh,gt,b,r0,EZz0,rmx,zmx,rmn,zmn;

  ESSetSplA(a);
  splRA(Lc,NULL,ESLc,ESLc2);
  L0	=2./Lc[0];
  ki	=0;
  for(k=1; k < ES2Mp1; k++){
    ki	+=ESNa1;
    splRA(Lc+k,NULL,ESLc+ki,ESLc2+ki);
    Lc[k]	*=L0;
    splRA(Ls+k,NULL,ESLs+ki,ESLs2+ki);
    Ls[k]	*=L0;
  }

  zmx=0;
  zmn=0;
  for(j=0; j < ESNp1; j++){
    if(zmx < z[j]){
      zmx	=z[j];
      rmx	=r[j];
      jmx	=j;
    }
    if(zmn > z[j]){
      zmn	=z[j];
      rmn	=r[j];
      jmn	=j;
    }
  }
  L0	=EZcr2*(z[jmx+1]-z[jmx-1]);
  EZz0	=z[jmx-1]+z[jmx+1]-2.*z[jmx];
  b	=-L0/EZz0;
  zmx	-=EZcr2*L0*L0/EZz0;
  L0	=EZcr2*(r[jmx+1]-r[jmx-1]);
  rmx	+=EZcr2*(r[jmx+1]-r[jmx-1]+(r[jmx-1]+r[jmx+1]-2.*r[jmx])*b)*b;

  L0	=EZcr2*(z[jmn+1]-z[jmn-1]);
  EZz0	=z[jmn-1]+z[jmn+1]-2.*z[jmn];
  b	=-L0/EZz0;
  zmn	-=EZcr2*L0*L0/EZz0;
  rmn	+=EZcr2*(r[jmn+1]-r[jmn-1]+(r[jmn-1]+r[jmn+1]-2.*r[jmn])*b)*b;

  r0	=ESaR0[0];
  EZz0	=EZcr2*(zmx+zmn);
  b	=EZcr2*(zmx-zmn);
  for(j=0; j < ESNp1; j++){
    gt	=(z[j]-EZz0)/b;
    if(gt > 1.){
      gt	=1.;
    }
    if(gt < -1.)
      gt	=-1.;
    if(gt > 0.){
      L0	=(rmx-r0)*(z[j]-EZz0)-(zmx-EZz0)*(r[j]-r0);
    }
    else{
      L0	=(r[j]-r0)*(zmn-EZz0)-(z[j]-EZz0)*(rmn-r0);
    }
    if(j < 2){
      gt=asin(gt);
    }
    else{
      if(j > ESNp-2){
	gt=EZc2gp+asin(gt);
      }
      else{
	if(L0 > 0.){
	  gt=EZcgp-asin(gt);
	}
	else{
	  gt=asin(gt);
	  if(gt < 0.){
	    gt	+=EZc2gp;
	  }
	}
      }
    }
    gh	=0.;
    for(k=1; k < ES2Mp1; k++){
      gh	+=(Lc[k]*sin(k*gt)+Ls[k]*(1.-cos(k*gt)))/k;
    }
#ifdef H
    gd[j]	=gt+gh-ESgt[j];
#endif
    gd[j]	=gh;
    gd1[j]	=ESgt[j]+gL[j]-gt;
  }
  return(0);
}

int  ESe2jsoPr(double *dpdgy,double *jbb)
{
  int i;
  double s,ss,ds,a;

  a	=0.;
  for(i=0; i < 21; i++){
    ss	=sqrt(0.05*i);
    do{
      ESSetSplA(a);
      splRA(&s,&ds,ESqgY,ESqgY2a);
      a		+=(ss-s)/ds;
    }while(fabs(ss-s) > 1e-6);
    ESSetSplA(a);
    splRA(&s,&ds,ESqgY,ESqgY2a);
    if(ESEqSolvInCr == 1){
      ESSetSplDCr(a);
      splRDCr(jbb+i,&ds,ESEqSolvInCr);
    }
    else{
      splRA(jbb+i,&ds,ESjb,ESjb2a);
    }
    jbb[i]	*=ESaR0[0];
    ESSetSplDPr(a);
    switch(ESEqSolvInPr){
    case 0:
      splRDPr(dpdgy+i,&ds,0);
      dpdgy[i]	/=ESaR0[0];
      break;
    case 1:
      splRDPr(dpdgy+i,&ds,1);
      break;
    default:
      splRA(dpdgy+i,&ds,ESPs,ESPs2a);
      break;
    }
    dpdgy[i]	*=ESgY[ESNa];
  }

  return(0);
}

double esjsopp_(double *gp)
{
  static int i=0;
  double s,ss,ds,a;
  
  if(*gp >= 1.){
    a	=1.;
  }
  else{
    if(*gp == 0.){
      a	=0.;
    }
    else{
      ss	=sqrt(*gp);
      while(ESqgY[i+1] < ss){
	i++;
      }
      while(ESqgY[i] > ss){
	i--;
      }
      a	=ESsa[i];
      do{
	ESSetSplA(a);
	splRA(&s,&ds,ESqgY,ESqgY2a);
	a		+=(ss-s)/ds;
      }while(fabs(ss-s) > 1e-6);
    }
  }
  switch(ESEqSolvInPr){
  case 0:
    ESSetSplDPr(a);
    splRDPr(&s,&ds,0);
    s	/=ESaR0[0];
    break;
  case 1:
    ESSetSplDPr(a);
    splRDPr(&s,&ds,1);
    break;
  default:
    ESSetSplA(a);
    splRA(&s,NULL,ESPs,ESPs2a);
    break;
  }
  return(s);
}

double esjsote3h_(double *gp)
{
  static int i=0;
  double s,ss,ds,a;
  
  if(*gp >= 1.){
    a	=1.;
  }
  else{
    if(*gp == 0.){
      a	=0.;
    }
    else{
      ss	=sqrt(*gp);
      while(ESqgY[i+1] < ss){
	i++;
      }
      while(ESqgY[i] > ss){
	i--;
      }
      a	=ESsa[i];
      do{
	ESSetSplA(a);
	splRA(&s,&ds,ESqgY,ESqgY2a);
	a		+=(ss-s)/ds;
      }while(fabs(ss-s) > 1e-6);
    }
  }
  if(ESEqSolvInCr == 1){
    ESSetSplDCr(a);
    splRDCr(&ss,&ds,ESEqSolvInCr);
  }
  else{
    ESSetSplA(a);
    splRA(&ss,&ds,ESjb,ESjb2a);
  }
  return(ss*ESaR0[0]);
}

int ESjso2eG(double *rv, double *zv)
{
  int i,j,ji,k,ki,kj;
  double s,ss,ds,a,b,rc[65],rs[65];

  ji	=0;
  for(j=0; j < ESNp1; j++){
    rv[ji]	=ESaR0[0];
    zv[ji]	=ESaZ0[0];
    ji++;
  }

  a	=ESsa[1];
  for(i=1; i < 33; i++){
    ss	=sqrt((double)i/32.);
    do{
      ESSetSplA(a);
      splRA(&s,&ds,ESqgY,ESqgY2a);
      a		+=(ss-s)/ds;
    }while(fabs(ss-s) > 1e-6);
    ESSetSplA(a);
    splRA(&b,&ds,ESsb,ESsb2a);
    splRA(rc,&ds,rcT,rcT2a);
    splRA(rs,&ds,rsT,rsT2a);
    ki		=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,&ds,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      splRA(rs+k,&ds,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
    }
    for(j=0; j < ESNp1; j++){
      rv[ji]	=ESaR0[0]+rc[0];
      zv[ji]	=ESaZ0[0]+rs[0]+b*ESsn1[j];
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	rv[ji]	+=rc[k]*EScs1[kj]+rs[k]*ESsn1[kj];
      }
      ji++;
    }
  }
  return(0);
}
int ESDconSaveSt(int ntor)
{
  FILE *fp;
  static double R=0.,a0=0.,b0=0.,d0=0.;
  double a,b,d;
  char ch;

#ifdef H
  printf("->(Y/N)");
  ch	='\n';
  while(isspace(ch)){
    ch=getchar();
  }
  a	=EZcr2*(rdPV[3]-rdPV[2]);
  b	=EZcr2*(zdPV[0]-zdPV[1]);
  d	=(EZcr2*(rdPV[3]+rdPV[2])-rdPV[0])/a;
  fp	=fopen("DconSt.d","a+");
  if(R != ESRext || a != a0 || b != b0 || d != d0){
    R	=ESRext;
    a0	=a;
    b0	=b;
    d0	=d;
    fprintf(fp,"%11s%11s%11s%11s\n","---------- R","a","b","d ----------");
    fprintf(fp,"%12.4e%11.4e%11.4e%11.4e\n",ESRext,a0,b0,d0);
    fprintf(fp,"%8s%8s%8s%8s%8s%8s%8s%7s%8s%8s%4s\n","B","gb%","q0","qa"
	    ,"jp(0)","j||(0)","j||(a)","n tor","gb_ext","j_a/j_0","Y/N");
  }
  fprintf(fp,"%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%7d%8.4f%8.4f   %c\n"
	  ,ESBt,ESgbext*100.,ESsq[0],ESsq[ESNa],ESjp[0],ESjb[0],ESjb[ESNa]
	  ,ntor,ESgbext,ESjb[ESNa]/ESjb[0],ch);
  printf("%c%2d%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n"
	 ,ch,ntor,ESBt,ESgbext*100.,ESsq[0],ESsq[ESNa],ESjp[0],ESjb[0]
	 ,ESjb[ESNa],ESgbext,ESjb[ESNa]/ESjb[0]);
  fclose(fp);
#endif
  return(0);
}


int ESDconVsPst(int k)
{
  FILE *fp;
  static double R=0.,a0=0.,b0=0.,d0=0.;
  double a,b,d;
  char ch;

  Vmt0=clock()-Vmt0;
  printf("->");
  ch	='\n';
  while(isspace(ch)){
    ch=getchar();
  }
  a	=EZcr2*(rdPV[3]-rdPV[2]);
  b	=EZcr2*(zdPV[0]-zdPV[1]);
  d	=(EZcr2*(rdPV[3]+rdPV[2])-rdPV[0])/a;
  fp=fopen("DconPst.d","a+");
  if(R != ESRext || a != a0 || b != b0 || d != d0){
    R	=ESRext;
    a0	=a;
    b0	=b;
    d0	=d;
    fprintf(fp,"%11s%11s%11s%11s\n","---------- R","a","b","d ----------");
    fprintf(fp,"%12.4e%11.4e%11.4e%11.4e\n",ESRext,a0,b0,d0);
    fprintf(fp,"%8s%8s%8s%8s%8s%8s%8s%7s%7s%7s\n","B","gb%","q0","qa","jp(0)"
	    ,"jb(0)","jb(a)","DCON","PEST4","PEST12");
  }
  switch(k){
  case 0:
    fprintf(fp,"%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%6.1f%c\n"
	    ,ESBt,ESgbext*100.,ESsq[0],ESsq[ESNa],ESjp[0],ESjb[0],ESjb[ESNa]
	    ,(double)Vmt0/CLOCKS_PER_SEC,ch);
    printf("DCON ");
    break;
  case 1:
    fprintf(fp,"%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%7s%6.1f%c\n"
	    ,ESBt,ESgbext*100.,ESsq[0],ESsq[ESNa],ESjp[0],ESjb[0],ESjb[ESNa]
	    ,"",(double)Vmt0/CLOCKS_PER_SEC,ch);
    printf("Pst4 ");
    break;
  case 2:
    fprintf(fp,"%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%7s%7s%6.1f%c\n"
	    ,ESBt,ESgbext*100.,ESsq[0],ESsq[ESNa],ESjp[0],ESjb[0],ESjb[ESNa]
	    ,"","",(double)Vmt0/CLOCKS_PER_SEC,ch);
    printf("Pst12 ");
    break;
  }
  fclose(fp);

  printf("%7.1f%c%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n"
	 ,(double)Vmt0/CLOCKS_PER_SEC,ch
	 ,ESBt,ESgbext*100.,ESsq[0],ESsq[ESNa],ESjp[0],ESjb[0],ESjb[ESNa]);
  return(0);
}

int ESEditDconInFile(double ttime, int ntor)
{
  int fp;
  char str[16],ch,*p;
#ifdef H
  fp=open("../dcon_3.10/dcon/dcon.in",O_RDWR);
#endif 
  fp=open("vac.in",O_RDWR);
  if(fp == -1){
    system("cp /w/Dcn/vac.in .");
  }
  else{
    close(fp);
  }
  fp=open("drawdcon.in",O_RDWR);
  if(fp == -1){
    system("cp /w/Dcn/drawdcon.in .");
  }
  else{
    close(fp);
  }
  fp=open("drawcrit.in",O_RDWR);
  if(fp == -1){
    system("cp /w/Dcn/drawcrit.in .");
  }
  else{
    close(fp);
  }

  fp=open("equil.in",O_RDWR);
  if(fp == -1){
    system("ln -s /w/Dcn/equil.in equil.in");
    fp=open("equil.in",O_RDWR);
  }
  sprintf(str,"ES4Dcon%4.2f",ttime);
  p	=str;
  ch	='\0';
  while(p-str != 7 && ch != EOF){
    read(fp,&ch,1);
    if(ch == *p){
      p++;
    }
    else{
      p	=str;
    }
  }
  if(ch == EOF){
    close(fp);
    printf("No Image 'ES4Dcon' in equil.in file\n");
    return(-1);
  }
  write(fp,p,4);
  close(fp);

  fp=open("dcon.in",O_RDWR);
  if(fp == -1){
    system("cp /w/Dcn/dcon.in .");
    fp=open("dcon.in",O_RDWR);
  }
  strcpy(str,"nn=");
  p	=str;
  ch	='\0';
  while(p-str != 3 && ch != EOF){
    read(fp,&ch,1);
    if(ch == *p){
      p++;
    }
    else{
      p	=str;
    }
  }
  if(ch == EOF){
    printf("No Image 'nn=' in dcon.in file\n");
    close(fp);
    return(-1);
  }
  switch(ntor){
  case 1:
    ch='1';
    break;
  case 2:
    ch='2';
    break;
  case 3:
    ch='3';
    break;
  case 4:
    ch='4';
    break;
  default:
    ch='1';
    break;
  }
  write(fp,&ch,1);
  close(fp);
  return(0);
}

int ESEditPstInFile(int ntor)
{
  int fp;
  char *str="n=",ch,*p;

#ifdef H
  fp=open("../Pst/fort.26",O_RDWR);
#endif

  fp=open("Pst/fort.20",O_RDWR);
  if(fp == -1){
    system("cp /scratchu/zakh/Pst/fort.20 Pst");
  }
  else{
    close(fp);
  }

  fp=open("Pst/fort.26",O_RDWR);
  if(fp == -1){
    system("cp /scratchu/zakh/Pst/fort.26 Pst");
    fp=open("Pst/fort.26",O_RDWR);
  }

  p	=str;
  ch	='\0';
  while(p-str != 2 && ch != EOF){
    read(fp,&ch,1);
    if(ch == *p){
      p++;
    }
    else{
      p	=str;
    }
  }
  if(ch == EOF){
    printf("No Image n= in fort.26 file\n");
    close(fp);
    return(-1);
  }
  switch(ntor){
  case 1:
    ch='1';
    break;
  case 2:
    ch='2';
    break;
  case 3:
    ch='3';
    break;
  default:
    ch='1';
    break;
  }
  write(fp,&ch,1);
  close(fp);
  return(0);
}

#ifdef XWIN
void ESOutPlot()
{
  int i,j,ji;
  double x[ESNa1],y[ESNa1];
  
  ZColor(14);
  for(j=0; j < ESNp; j++){
    ji	=j;
    for(i=0; i < ESNa1; i++){
      x[i]	=ESsr[ji];
      y[i]	=ESsz[ji];
      ji	+=ESNp1;
    }
    ZPlotPolyLine(x,y,ESNa1);
  }
  return;
}

void ESOutPlot0()
{
  int i,j,ji;
  double x[ESNa1],y[ESNa1];
  
  ZColor(14);
  for(j=0; j < ESNp; j++){
    ji	=j;
    for(i=0; i < ESNa1; i++){
      x[i]	=ESsr[ji];
      y[i]	=ESsz[ji];
      ji	+=ESNp1;
    }
    ZPlotPolyLine(x,y,ESNa1);
  }
  ZColor(4);
  ji	=0;
  for(i=1; i < ESNa1; i++){
    ji	+=ESNp1;
    ZPlotPolyLine(ESsr+ji,ESsz+ji,ESNp1);
  }
  return;
}

void ESOutPlot1()
{
  int i,j,ji;
  double x[VmNa1],y[VmNa1];

  ZColor(14);
  for(j=0; j < VmNp; j++){
    ji	=j;
    for(i=0; i < VmNa1; i++){
      x[i]	=R[ji];
      y[i]	=Z[ji];
      ji	+=VmNp1;
    }
    ZPlotPolyLine(x,y,VmNa1);
  }

  ji	=0;
  for(i=1; i < VmNa1; i++){
    ji	+=VmNp1;
    ZPlotPolyLine(R+ji,Z+ji,VmNp1);
  }

  return;
}

void ESOutPlot2()
{
  int i,j,k,n,n1,errFl;
  double x[8],y[8];
  double r,z,r0,EZz0,r1,z1,dr,dz;
  
  ESSetMapBox();
  dr	=(rMap[1]-rMap[0])/nrMap;
  dz	=(zMap[1]-zMap[0])/nzMap;

  ZColor(12);
  z1	=zMap[0];
  k	=0;
  for(i=0; i < nzMap; i++){
    EZz0	=z1;
    z1	=EZz0+dz;
    r1	=rMap[0];
    for(j=0; j < nrMap; j++){
      r0	=r1;
      r1	=r0+dr;
      if(iAMap[k]){
	x[0]=r0;
	y[0]=EZz0;
	x[1]=r1;
	y[1]=EZz0;
	x[2]=r1;
	y[2]=z1;
	x[3]=r0;
	y[3]=z1;
	x[4]=x[0];
	y[4]=y[0];
	ZPlotPolyLine(x,y,5);
      }
      k++;
    }
  }

  z1	=zMap[0];
  k	=0;
  errFl=0;
  ZColor(14);
  for(i=0; i < nzMap; i++){
    EZz0	=z1;
    z1	=EZz0+dz;
    r1	=rMap[0];
    for(j=0; j < nrMap; j++){
      r0	=r1;
      r1	=r0+dr;
      if(iAMap[k]){
	x[0]=r0;
	y[0]=EZz0;
	x[1]=r1;
	y[1]=EZz0;
	x[2]=r1;
	y[2]=z1;
	x[3]=r0;
	y[3]=z1;
	n	=0;
	switch(iAMap[k]){
	case 1:
#ifdef H
	  errFl	|=ES1DMapZ2gt(x+3,x+4,x+1,EZz0);
	  x[2]=r0;
	  errFl	|=ES1DMapR2gt(x+4,x+3,y+4,y+2,r0);
	  x[3]	=x[1];
	  y[3]	=y[2];
	  ZPlotPolyLine(x+1,y+1,3);
	  n	=3;
#endif
	  errFl	|=ES1DMapZ2gt(x+5,x+4,x+1,EZz0);
	  errFl	|=ES1DMapR2gt(x+4,x+5,y+4,y+3,r0);
	  x[2]	=x[1];
	  y[2]	=y[3];
	  n	=4;

	  break;
	case 2:
#ifdef H
	  errFl	|=ES1DMapZ2gt(x+4,x,x+3,EZz0);
	  errFl	|=ES1DMapR2gt(x+4,y+4,y+3,y+2,r1);
	  x[3]	=x[0];
	  y[3]	=y[2];
	  x[4]	=x[0];
	  y[4]	=y[0];
	  ZPlotPolyLine(x+2,y+2,3);
	  n	=3;
#endif
	  errFl	|=ES1DMapZ2gt(x+4,x,x+5,EZz0);
	  errFl	|=ES1DMapR2gt(x+4,y+4,y+5,y+2,r1);
	  x[3]	=x[0];
	  y[3]	=y[2];
	  n	=4;
	  break;
	case 300:
	  errFl	|=ES1DMapR2gt(x+4,x+4,y+4,y+2,r1);
	  errFl	|=ES1DMapR2gt(x+4,x+4,y+4,y+3,r0);
	  x[4]	=x[0];
	  x[5]	=x[1];
	  y[4]	=y[3] < y[2] ? y[3] : y[2];
	  y[5]	=y[4];
	  ZPlotPolyLine(x+4,y+4,2);

	  y[4]	=y[3] > y[2] ? y[3] : y[2];
	  y[5]	=y[4];
	  ZPlotPolyLine(x+4,y+4,2);
	  n	=4;
	  break;
	case 4:
#ifdef H
	  errFl	|=ES1DMapZ2gt(x+4,x,x+3,z1);
	  y[0]	=z1;
	  errFl	|=ES1DMapR2gt(x+4,y+4,y+1,y+3,r1);
	  x[3]	=x[1];
	  y[3]	=y[1];
	  x[4]	=x[0];
	  y[4]	=y[1];
	  x[5]	=x[0];
	  y[5]	=y[0];
	  ZPlotPolyLine(x+3,y+3,3);
	  n	=3;
#endif
	  errFl	|=ES1DMapZ2gt(x+4,x+3,x+5,z1);
	  errFl	|=ES1DMapR2gt(x+4,y+4,y+1,y+5,r1);
	  x[0]	=x[3];
	  y[0]	=y[1];
	  n	=4;
	  break;
	case 500:
	  errFl	|=ES1DMapZ2gt(x+4,y+4,x+1,EZz0);
	  errFl	|=ES1DMapR2gt(x+6,x+7,y+6,y+2,r1);
	  errFl	|=ES1DMapZ2gt(x+4,x+3,y+4,z1);
	  x[4]	=r0;
	  errFl	|=ES1DMapR2gt(x+6,x+7,y+6,y+4,r0);
	  x[5]	=r0;
	  y[5]	=EZz0;
	  n	=6;
	  break;
	case 600:
	  errFl	|=ES1DMapZ2gt(x+4,x,y+4,EZz0);
	  errFl	|=ES1DMapZ2gt(x+4,x+3,y+4,z1);
	  y[4]	=y[0];
	  y[5]	=y[2];
	  x[4]	=x[3] < x[0] ? x[0] : x[3];
	  x[5]	=x[4];
	  ZPlotPolyLine(x+4,y+4,2);

	  x[4]	=x[3] > x[0] ? x[0] : x[3];
	  x[5]	=x[4];
	  ZPlotPolyLine(x+4,y+4,2);
	  n	=4;
	  break;
	case 7:
#ifdef H
	  errFl	|=ES1DMapZ2gt(x+6,x+3,y+6,z1);
	  x[4]	=r0;
	  errFl	|=ES1DMapR2gt(x+6,y+6,y+6,y+4,r0);
	  n	=5;
#endif
	  n	=4;
	  break;
	case 8:
#ifdef H
	  x[0]	=x[3];
	  y[0]	=y[3];
	  errFl	|=ES1DMapZ2gt(y+4,x+4,x+2,z1);
	  x[1]	=r0;
	  errFl	|=ES1DMapR2gt(x+4,y+4,y+1,y+3,r0);
	  x[3]	=x[2];
	  y[3]	=y[1];
	  x[4]	=x[0];
	  y[4]	=y[1];
	  ZPlotPolyLine(x+2,y+2,3);
	  n	=3;
#endif
	  errFl	|=ES1DMapZ2gt(y+4,x+4,x+2,z1);
	  errFl	|=ES1DMapR2gt(x+4,y+4,y,y+5,r0);
	  x[1]	=x[2];
	  y[1]	=y[0];
	  n	=4;
	  break;
	case 900:
	  errFl	|=ES1DMapZ2gt(x+4,y+4,x+1,EZz0);
	  errFl	|=ES1DMapZ2gt(x+4,y+4,x+2,z1);
	  y[4]	=y[0];
	  y[5]	=y[2];
	  x[4]	=x[1] < x[2] ? x[1] : x[2];
	  x[5]	=x[4];
	  ZPlotPolyLine(x+4,y+4,2);

	  x[4]	=x[1] > x[2] ? x[1] : x[2];
	  x[5]	=x[4];
	  ZPlotPolyLine(x+4,y+4,2);
	  n	=4;
	  break;
	case 10:
	  errFl	|=ES1DMapZ2gt(x+4,x,y+4,EZz0);
	  errFl	|=ES1DMapR2gt(x+6,x+7,y+6,y+2,r1);
	  errFl	|=ES1DMapZ2gt(x+4,y+4,x+3,z1);
	  x[4]	=r0;
	  y[4]	=z1;
	  x[5]	=r0;
	  errFl	|=ES1DMapR2gt(x+6,x+7,y+5,y+6,r0);
	  n	=6;
	  break;
	case 11:
#ifdef H
	  errFl	|=ES1DMapR2gt(x+4,y+4,y+6,y+2,r1);
	  x[4]	=r0;
	  y[4]	=z1;
	  errFl	|=ES1DMapZ2gt(x+6,y+6,x+3,z1);
	  n	=5;
#endif
	  n	=4;
	  break;
	case 1200:
	  errFl	|=ES1DMapR2gt(x+6,x+7,y,y+6,r0);
	  errFl	|=ES1DMapR2gt(x+6,x+7,y+1,y+6,r1);
	  x[4]	=x[2];
	  x[5]	=x[3];
	  y[4]	=y[1] < y[0] ? y[0] : y[1];
	  y[5]	=y[4];
	  ZPlotPolyLine(x+4,y+4,2);

	  y[4]	=y[1] > y[0] ? y[0] : y[1];
	  y[5]	=y[4];
	  ZPlotPolyLine(x+4,y+4,2);
	  n	=4;
	  break;
	case 13:
#ifdef H
	  errFl	|=ES1DMapZ2gt(x+6,y+6,x+1,EZz0);
	  errFl	|=ES1DMapR2gt(x+6,y+6,y+2,y+6,r1);
	  x[3]	=r1;
	  y[3]	=z1;
	  x[4]	=r0;
	  y[4]	=z1;
	  n	=5;
#endif
	  n	=4;
	  break;
	case 14:
#ifdef H
	  errFl	|=ES1DMapZ2gt(x+6,x,y+6,EZz0);
	  x[4]	=r0;
	  errFl	|=ES1DMapR2gt(x+6,y+6,y+4,y+6,r0);
	  n	=5;
#endif
	  n	=4;
	  break;
	case 15:
	  n	=4;
	  break;
	}
	n1	=n+1;
	x[n]	=x[0];
	y[n]	=y[0];
	if(n)ZPlotPolyLine(x,y,n1);
      }
      k++;
    }
  }
  if(errFl) CbStr2UserMessage("Some Cells are not found\n");
  return;
}
#endif

int ES2Orbit(int Na,int Np)
{
  static int Fic=0;
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

  double Bmd[ESNa+ESNa1],Xmd[ESNa+ESNa1],*lBmd,*lXmd;

  lBmd	=Bmd+ESNa;
  lXmd	=Xmd+ESNa;
#else
  double rt,zt,ra,za,B[129],aJ[129],L[129],A[129];
  double r[128],z[128];
  double rc[33],rs[33],rca[33],rsa[33],rct[33],rst[33],b,ba;

  if(Na > 128){
    Na	=128;
  }
  if(Np > 128){
    Np	=128;
  }
#endif
  time(&stTime);
  ltm	=localtime(&stTime);

  sprintf(FNm,"%2.2d%2.2d%2.2dES2Orbit.%2.2d"
	  ,ltm->tm_mon+1,ltm->tm_mday,ltm->tm_year%100,Fic);
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

#ifdef H
    A[i]	=da*i;
#endif

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
    R[j]	=r[0];
    Z[j]	=z[0];
    fprintf(lf,"%14.11f %14.11f %14.11f %14.7e %3d %3d\n"
	    ,r[0],z[0],1.,0.,0,j);
  }
#ifdef XWIN
  lBmd[0]	=ESaF[0]/r[0];
  lXmd[0]	=r[0];
#endif

  R[j]	=r[0];
  Z[j]	=z[0];
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
      R[ji]	=r[j];
      Z[ji]	=z[j];
      ji++;
    }
    L0	/=Np;
    s	=-aF*aF/dgy;
    for(j=0; j < Np; j++){
      aJ[j]	+=s*(L[j]-L0);
      fprintf(lf,"%14.11f %14.11f %14.11f %14.7e %3d %3d\n"
	      ,r[j],z[j],B[j]*rB,aJ[j]*rB,i,j);
    }
    R[ji]	=r[0];
    Z[ji]	=z[0];
    ji++;
    B[j]	=B[0];
    aJ[j]	=aJ[0];

    s	=0.;
    for(j=0; j < Np; j++){
      s	+=aJ[j]*rB;
    }
#ifdef XWIN
    s	/=Np;
    EZout("sid","I",i,s);

    j	=Np/2;
    lBmd[i]	=B[0];
    lBmd[-i]	=B[j];
    lXmd[i]	=r[0];
    lXmd[-i]	=r[j];
#endif
  }
  fclose(lf);
#ifdef XWIN
  printf("%s has been written\n",FNm);
  Fic++;
  j	=ESNa+ESNa1;
  Scale2D(4,2,Xmd,Bmd,j,6,2);
  Plot2d(4,Xmd,Bmd,j,6,0,14,0);
  CbFlush();
#endif
  return(0);
}

int ES2OrbitS(int Na,int Np)
{
  static int Fic=0;
  int Na1,Np1;
  FILE *lf;
  char FNm[24];

  int i,j,ji,k,kj,ki;

  double a,da,t,dt,gy,sq,sp,aJ,aF,qg,rB;
  double dgy,s,cs,sn;
#ifdef XWIN
  double r,z,rt,zt,ra,za,B,A[Na+1];
  double rc[ESFp1],rs[ESFp1],rca[ESFp1],rsa[ESFp1],rct[ESFp1],rst[ESFp1],b,ba;
#else
  double r,z,rt,zt,ra,za,B,A[129];
  double rc[33],rs[33],rca[33],rsa[33],rct[33],rst[33],b,ba;
#endif

  time(&stTime);
  ltm	=localtime(&stTime);
 
  sprintf(FNm,"%2.2d%2.2d%2.2dES2Orbit.%2.2d"
	  ,ltm->tm_mon+1,ltm->tm_mday,ltm->tm_year%100,Fic);
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
    R[j]	=r;
    Z[j]	=z;
    fprintf(lf,"%14.11f %14.11f %14.11f %14.7e %3d %3d\n"
	    ,r,z,1.,0.,0,j);
  }
  R[j]	=r;
  Z[j]	=z;
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
    for(j=0; j < Np; j++){
      t		=-dt*j;
      cs	=cos(t);
      sn	=sin(t);
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
      ra	=(rt*rt+zt*zt);
      za	=-dgy/qg;
      B		=sqrt(ra*za*za+s*s);
      fprintf(lf,"%14.11f %14.11f %14.11f %14.7e %3d %3d\n"
	      ,r,z,B*rB,ra*za*rB,i,j);

      R[ji]	=r;
      Z[ji]	=z;
      ji++;
    }
    R[ji]	=R[ji-Np];
    Z[ji]	=Z[ji-Np];
    ji++;
  }
  fclose(lf);
  printf("%s has been written\n",FNm);
  Fic++;
  return(0);
}

/* Generator of Boozer coordinates */
int ES2Boozer(int Na,int Np)
{
  static int Fic=0;
  int Na1,Np1;
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

  time(&stTime);
  ltm	=localtime(&stTime);
 
  sprintf(FNm,"%2.2d%2.2d%2.2dES2Orbit.%2.2d"
	  ,ltm->tm_mon+1,ltm->tm_mday,ltm->tm_year%100,Fic);
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

/****/
    A[i]	=da*i;


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
    R[j]	=r;
    Z[j]	=z;
    fprintf(lf,"%14.11f %14.11f %14.11f %14.7e %3d %3d\n"
	    ,r,z,1.,0.,0,j);
  }
  R[j]	=r;
  Z[j]	=z;
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
	s	+=(Lc[k]*ESsn1[kj]-Ls[k]*EScs1[kj])/k;
      }
      Gh[j]	=s;
    }
    Gh[j]	=Gh[0];
    if(i == ESNa-1){
      EZout("siddd","gq Booser",j,Gh[j],ESgt[j],(Gh[j]-ESgt[j])*EZcr2gp*ESsq[i]);
    }
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
      rst[k]	=rc[k]*k;
      splRA(rs+k,rsa+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      rct[k]	=-rs[k]*k;
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
      zt	=-b*cs;
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
      R[ji]	=r;
      Z[ji]	=z;
      ji++;
    }
    R[ji]	=R[ji-Np];
    Z[ji]	=Z[ji-Np];
    ji++;
  }
  fclose(lf);
#ifdef XWIN
  printf("%s has been written\n",FNm);
#endif
  Fic++;
  return(0);
}

/* Generator of PEST coordinates */
int ES2PEST(int Na,int Np)
{
  int Na1,Np1;
  char FNm[24];

  int i,j,ji,k,kj,ki,iK;
  double a,da,t,dt;
  double s,s1,s2,cs,sn,gh;

#ifdef XWIN
  double r,z,rt,zt,ra,za,B[Np+1];
  double rc[ESFp1],rs[ESFp1],rca[ESFp1],rsa[ESFp1],rct[ESFp1],rst[ESFp1],b,ba;
  double Ght[ESNp1],Ghtc[ESFp1],Ghts[ESFp1];
#else
  double r,z,rt,zt,ra,za,B[129];
  double rc[33],rs[33],rca[33],rsa[33],rct[33],rst[33],b,ba;
  double Ght[129],Ghtc[33],Ghts[33];
#endif
  Na1	=Na+1;
  Np1	=Np+1;
  da	=ESsa[ESNa]/Na;
  dt	=EZc2gp/Np;

  r	=ESaR0[0];
  z	=ESaZ0[0];
  for(ji=0; ji < Np1; ji++){
    R[ji]	=r;
    Z[ji]	=z;
  }
  for(i=1; i < Na1; i++){
    a	=da*i;
    ESSetSplA(a);
    splRA(&b,&ba,ESsb,ESsb2a);
    splRA(rc,rca,rcT,rcT2a);
    rc[0]	+=ESaR0[0];
    splRA(rs,rsa,rsT,rsT2a);
    rs[0]	+=ESaZ0[0];
    ki	=0;
    iK	=1;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,rca+k,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      rst[k]	=-rc[k]*k;
      rca[k]	*=2.;
      splRA(rs+k,rsa+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      rct[k]	=rs[k]*k;
      rsa[k]	*=2.;
      if(rct[k] != 0. || rca[k] != 0. || rst[k] != 0. || rsa[k] != 0.){
	iK	=k+1;
      }
    }
    s	=0.;
    for(j=0; j < ESNp; j++){
      cs	=EScs1[j];
      sn	=ESsn1[j];
      r		=rc[0];
      ra	=rca[0];
      rt	=0.;
      z		=rs[0]+b*sn;
      za	=rsa[0]+ba*sn;
      zt	=b*cs;
      kj	=0;
      for(k=1; k < iK; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	cs	=EScs1[kj];
	sn	=ESsn1[kj];
	r	+=rc[k]*cs+rs[k]*sn;
	ra	+=rca[k]*cs+rsa[k]*sn;
	rt	+=rct[k]*cs+rst[k]*sn;
      }
      Ght[j]	=(ra*zt-rt*za)/r;
      s		+=Ght[j];
    }
    s	=ESNp/s;
    for(j=0; j < ESNp; j++){
      Ght[j]	=Ght[j]*s-1.;
    }
    Ght[j]	=Ght[0];
    ESP2F(Ghtc,Ghts,Ght,ESFp);
    for(k=1; k < ESFp1; k++){
      Ghtc[k]	*=2.;
      Ghts[k]	*=2.;
    }
    t	=0.;
    for(j=0; j < Np; j++){
      gh	=dt*j;
      kj	=0;
      do{
	s	=t-gh;
	s1	=1.;
	for(k=1; k < ESFp1; k++){
	  cs	=cos(k*t);
	  sn	=sin(k*t);
	  s1	+=Ghtc[k]*cs+Ghts[k]*sn;
	  s	+=(Ghtc[k]*sn+Ghts[k]*(1.-cs))/k;
	}
	t	-=s/s1;
	kj++;
      }while(kj < 10 && fabs(s) > 1e-8);
      if(kj == 10){
	printf("i=%2d j=%d gh=%10.3e - No convergence in ES2PEST()%c\n"
	       ,i,j,gh,7);
      }

      s		=t;
      cs	=cos(s);
      sn	=sin(s);
      r		=rc[0]+rc[1]*cs+rs[1]*sn;
      ra	=rca[0]+rca[1]*cs+rsa[1]*sn;
      rt	=rct[1]*cs+rst[1]*sn;
      z		=rs[0]+b*sn;
      za	=rsa[0]+ba*sn;
      zt	=-b*cs;
      for(k=2; k < ESFp1; k++){
	s	=k*t;
	cs	=cos(s);
	sn	=sin(s);
	r	+=rc[k]*cs+rs[k]*sn;
      }
      R[ji]	=r;
      Z[ji]	=z;
      ji++;
    }
    R[ji]	=R[ji-Np];
    Z[ji]	=Z[ji-Np];
    ji++;
  }
  return(0);
}

/* Generator of Boozer coordinates */
int ES2canonicalB(int Na,int Np,int Fl)
{
  static int Fic=0;
  int Na1,Np1;
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

  time(&stTime);
  ltm	=localtime(&stTime);
 
  sprintf(FNm,"%2.2d%2.2d%2.2dES2Orbit.%2.2d"
	  ,ltm->tm_mon+1,ltm->tm_mday,ltm->tm_year%100,Fic);
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

/****/
    A[i]	=da*i;


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
    R[j]	=r;
    Z[j]	=z;
    fprintf(lf,"%14.11f %14.11f %14.11f %14.7e %3d %3d\n"
	    ,r,z,1.,0.,0,j);
  }
  R[j]	=r;
  Z[j]	=z;
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
      rst[k]	=rc[k]*k;
      splRA(rs+k,rsa+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      rct[k]	=-rs[k]*k;
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
      zt	=-b*cs;
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
      R[ji]	=r;
      Z[ji]	=z;
      ji++;
    }
    R[ji]	=R[ji-Np];
    Z[ji]	=Z[ji-Np];
    ji++;
  }
  fclose(lf);
#ifdef XWIN
  printf("%s has been written\n",FNm);
#endif
  Fic++;
  return(0);
}

/* Generator of Canonical  coordinates: It is impossible */
int ES2Canonical(int Na,int Np,int Fl)
{
  int Na1,Np1;
  int i,j,k,kj,ki;
  double a,da,t,dt,aF;
  double dgy,s,cs,sn,x[2],y[2];
  double r,z,rt,zt,ra,za;

#ifdef XWIN
  double gq[Np+1],bgz[Np+1],dbgz[Np+1];
  double rc[ESFp1],rs[ESFp1],rca[ESFp1],rsa[ESFp1],b,ba;

  Na1	=Na+1;
  Np1	=Np+1;
  da	=ESsa[ESNa]/Na;
  dt	=EZc2gp/Np;

  y[0]	=0.;
  y[1]	=0.;
  for(j=0; j < Np1; j++){
    gq[j]	=-.5+dt*j*EZcr2gp;
    dbgz[j]	=0.;
    bgz[j]	=0.;
  }
  for(i=1; i < Na1; i++){
    a	=da*i;
    ESSetSplA(a);
    splRA(&dgy,NULL,ESdgY,ESdgY2a);
    splRA(&aF,NULL,ESaF,ESaF2a);
    dgy	*=EZcr2*a*da*EZcr2gp/aF;
    ki	=0;
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
      splRA(rs+k,rsa+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      rsa[k]	*=2.;
    }
    for(j=0; j < Np; j++){
      t		=EZc2gp-dt*j;
      cs	=cos(t);
      sn	=sin(t);
      r		=rc[0]+rc[1]*cs+rs[1]*sn;
      ra	=rca[0]+rca[1]*cs+rsa[1]*sn;
      rt	=rc[1]*sn-rs[1]*cs;
      z		=rs[0]+b*sn;
      za	=rsa[0]+ba*sn;
      zt	=-b*cs;
      for(k=2; k < ESFp1; k++){
	s	=k*t;
	cs	=cos(s);
	sn	=sin(s);
	r	+=rc[k]*cs+rs[k]*sn;
	ra	+=rca[k]*cs+rsa[k]*sn;
	rt	+=k*(rc[k]*sn-rs[k]*cs);
      }
      bgz[j]	+=dbgz[j];
      dbgz[j]	=-(ra*rt+za*zt)/(r*(rt*za-ra*zt))*dgy;
      bgz[j]	+=dbgz[j];
      s	=bgz[j];
      if(y[0] > s){
	y[0]	=s;
      }
      if(y[1] < s){
	y[1]	=s;
      }
    }
    bgz[j]	=bgz[0];
  }
  x[0]	=-0.5;
  x[1]	=0.5;

  y[0]	=-0.05;
  y[1]	=0.05;

  if(Fl){
    static int kPS=1;
    static char PSFileNm[16];
    if(kPS && *ShotNm != '\0'){
      sprintf(PSFileNm,"%sCan.ps",ShotNm);
      PSSetFileName(PSFileNm);
      kPS	=0;
    }
    SetPlotName("(zeta-phi)/2pi","theta/2pi","Canonical zeta-phi");
    PSFileOpen(1,1);
  }
  Scale2d(2,x,y,2,6,2);
  if(Fl){
    PSNewFrame();
  }

  for(j=0; j < Np1; j++){
    dbgz[j]	=0.;
    bgz[j]	=0.;
  }
  for(i=1; i < Na1; i++){
    a	=da*i;
    ESSetSplA(a);
    splRA(&dgy,NULL,ESdgY,ESdgY2a);
    splRA(&aF,NULL,ESaF,ESaF2a);
    dgy	*=EZcr2*a*EZcr2gp*da/aF;
    ki	=0;
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
      splRA(rs+k,rsa+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      rsa[k]	*=2.;
    }
    for(j=0; j < Np; j++){
      t		=EZc2gp-(dt*j-EZcgp);
      cs	=cos(t);
      sn	=sin(t);
      r		=rc[0]+rc[1]*cs+rs[1]*sn;
      ra	=rca[0]+rca[1]*cs+rsa[1]*sn;
      rt	=rc[1]*sn-rs[1]*cs;
      z		=rs[0]+b*sn;
      za	=rsa[0]+ba*sn;
      zt	=-b*cs;
      for(k=2; k < ESFp1; k++){
	s	=k*t;
	cs	=cos(s);
	sn	=sin(s);
	r	+=rc[k]*cs+rs[k]*sn;
	ra	+=rca[k]*cs+rsa[k]*sn;
	rt	+=k*(rc[k]*sn-rs[k]*cs);
      }
      bgz[j]	=bgz[j]+dbgz[j];
      dbgz[j]	=-(ra*rt+za*zt)/(r*(rt*za-ra*zt))*dgy;
      bgz[j]	=bgz[j]+dbgz[j];
    }
    bgz[j]	=bgz[0];
    Plot2d(2,gq,bgz,Np1,6,0,0,0);
    if(Fl){
      PSPlot2d(gq,bgz,Np1,6,0,0,0);
    }
  }
  if(Fl){
    PSFileClose();
  }
#endif
  return(0);
}

/* Generator of Orthogonal  coordinates: It is impossible */
int ES2Orthogonal(int Na,int Np)
{
  static int Fic=0;
  int Na1,Np1;
  FILE *lf;
  char FNm[24];

  return(0);
}

#ifdef H
/* Mapper to Boozer coordinates */
int ES4Pst3(int nA,double ttime)
{
  int nA1=NAP+1,nP=NPP,nP1=NPP+1,nP5=NPP+5;
  FILE *Fp,*fp;
  int i,j,k,ki,kj,ji;
  double dgy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double A[NAP+1];
  double *sp,*Ps,sq[NAP+1],dsqgy[NAP+1],*g,*dggy,*fb,*dfbgy;
  double *r,*z,*rt,*zt,*ra,*za;
  double *g22,*p2r,*g11,*g12,*dg12t,*dp2ra,*dg22t,*dp2rt,*qg,*dqga,*gd,*dgda;

  double *rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,EZdza,EZdzt,cs,sn,b,db;
  double *gQ,*gQ2,gq,gqa,gqt,gqat,gqtt,t;
  double *d2rttc,*d2rtts,*d2ratc,*d2rats,d2rtt,d2ztt,d2rat,d2zat;
  double *Lc,*Ls,*dLc,*dLs,D0,dD0,D,dD;
  double dg22,dg12,gh,ght,gha,ghat,ghtt,L0,G0;
  double *rDc,*rDc2,*rDs,*rDs2,*Jc,*dJc,*Js,*dJs;

  double *w1[6],*w2[6],*w3[6];

  nA1	=nA+1;
  i	=nP*nA1;
  for(k=0; k < 6; k++){
    w1[k]	=NULL;
    w2[k]	=NULL;
    w3[k]	=NULL;
    w1[k]	=(double*)malloc(i*sizeof(double));
    w2[k]	=(double*)malloc(i*sizeof(double));
    w3[k]	=(double*)malloc(i*sizeof(double));
    if(w1[k] == NULL || w2[k] == NULL || w3[k] == NULL){
      printf("No memory for w[%d] in ES4Pst1()%c\n",k,7);
      j=k;
      while(j >= 0){
	free(w3[j]);
	free(w2[j]);
	free(w1[j]);
	j--;
      }
      return(0);
    }
  }

  r	=w1[0];
  z	=w1[1];
  rt	=w1[2];
  zt	=w1[3];
  ra	=w1[4];
  za	=w1[5];

  g22	=w2[0];
  p2r	=w2[1];
  g11	=w2[2];
  g12	=w2[3];
  dg12t	=w2[4];
  dp2ra	=w2[5];
  dg22t	=w3[0];
  dp2rt	=w3[1];
  qg	=w3[2];
  dqga	=w3[3];
  gd	=w3[4];
  dgda	=w3[5];
  
  gQ	=(double*)malloc((4*ESNp1+18*ESFp1)*sizeof(double));
  gQ2	=gQ	+ESNp1;
  rc	=gQ2	+ESNp1;
  rs	=rc	+ESFp1;
  drc	=rs	+ESFp1;
  drs	=drc	+ESFp1;
  rct	=drs	+ESFp1;
  rst	=rct	+ESFp1;
  d2ratc	=rst	+ESFp1;
  d2rats	=d2ratc	+ESFp1;
  d2rttc	=d2rats	+ESFp1;
  d2rtts	=d2rttc	+ESFp1;
  Lc		=d2rtts	+ESFp1;
  Ls		=Lc	+ESFp1;
  dLc		=Ls	+ESFp1;
  dLs		=dLc	+ESFp1;
  Jc		=dLs	+ESFp1;
  Js		=Jc	+ESFp1;
  dJc		=Js	+ESFp1;
  dJs		=dJc	+ESFp1;

  if(gQ == NULL){
    printf("No memory for gQ in ES4Pst1()\n");
    return(0);
  }
  sp	=r;
  Ps	=z;
  g	=ra;
  dggy	=za;
  fb	=g22;
  dfbgy	=p2r;

  i	=ESNa1*ESFp1;
  rDc	=(double*)malloc(4*i*sizeof(double));
  rDc2	=rDc	+i;
  rDs	=rDc2	+i;
  rDs2	=rDs	+i;
  if(rDc == NULL){
    printf("No memory for gQ in ES4Pst1()\n");
    return(0);
  }

  i	=ESNa1*ES2Mp1;
  for(ki=0; ki < i; ki++){
    rDc[ki]	=ESLc[ki]	+ESVc[ki];
    rDc2[ki]	=ESLc2[ki]	+ESVc2[ki];
    rDs[ki]	=ESLs[ki]	+ESVs[ki];
    rDs2[ki]	=ESLs2[ki]	+ESVs2[ki];
  }

  rg	=1./ESRBt;
  rBt	=ESRext*rg;
  p2rBt	=rBt*rBt;
  rf	=-ESgY[ESNa]*rBt;
  dgt	=EZc2gp/nP;
  dgy	=1./nA;
  
#ifdef H
  sprintf(Ttl.Title,"ES4PstGm%4.2f",ttime); /* fort.31 */
  sprintf(Ttl.Title,"ES4PstGm");
  Fp	=fopen(Ttl.Title,"w");

  sprintf(Ttl.Title,"ES4PstPr%4.2f",ttime); /* fort.32 */
  sprintf(Ttl.Title,"ES4Pst");
  fp	=fopen(Ttl.Title,"w");
#endif
  fp=fopen("Pst","r");
  if(fp == NULL){
    system("mkdir Pst");
  }
  else{
    fclose(fp);
  }
  sprintf(Ttl.Title,"Pst/fort.31"); /* fort.31 */
  Fp	=fopen(Ttl.Title,"w");

  sprintf(Ttl.Title,"Pst/fort.32"); /* fort.32 */
  fp	=fopen(Ttl.Title,"w");

  printf("Writing plasma profiles\n");
  strcpy(Ttl.Title,"ESC has written this garbage");
  i	=strlen(Ttl.Title);
  while(i < 159){
    Ttl.Title[i]	=' ';
    i++;
  }
  Ttl.Title[i]	='\0';
  g[0]	=0.;
  fwrite((void*)g,(size_t)8,(size_t)1,Fp);
  fwrite((void*)g,(size_t)8,(size_t)1,fp);
  i	=0;
  while(i < 8){
    Ttl.date[i]	='\0';
    i++;
  }
  Ttl.nxx[0]	=nP1; /* number of poloidal intervals (+1 in symmetric case)*/
  Ttl.nxx[1]	=nA1;			/* number of radial intervals +1*/ 
  Ttl.nxx[2]	=401;			/* Number of magnetic surfaces */
  Ttl.nxx[3]	=128;			/* Number of poloidal intervals */ 
  Ttl.nxx[4]	=0;			/* Junk */ 
  Ttl.axx[0]	=rdPV[3]-rdPV[2];
  Ttl.axx[1]	=zdPV[0]-zdPV[1];
  Ttl.axx[2]	=zdPV[0];
  Ttl.axx[3]	=ESaR0[0];
  Ttl.axx[4]	=ESRext;
  Ttl.axx[5]	=ESsp[0]*p2rBt;
  Ttl.axx[6]	=ESRBt*rBt/ESRext;
  Ttl.axx[7]	=EZc2gp*ESgY[ESNa]*rBt;
  Ttl.axx[8]	=0.;
  Ttl.axx[9]	=EZc2gp*ESgY[ESNa]*rBt;
  Ttl.axx[10]	=0.;
  Ttl.axx[11]	=0.;
  Ttl.axx[12]	=1.;	/* constant C in rD=Cr^m\na\bgY^n */

  Ttl.nxy[0]	=0;	/* exponent n in rD=Cr^m\na\bgY^n */
  Ttl.nxy[1]	=2;	/* exponent m in rD=Cr^m\na\bgY^n,
			   2 - straight field lines */
  Ttl.nxy[2]	=0;	/* Junk */ 
  Ttl.nxy[3]	=0;	/* Junk */ 
  Ttl.nxy[4]	=0;	/* Junk */
  Ttl.nxy[5]	=0;	/* Junk */
  Ttl.nxy[6]	=0;	/* Junk */ 
  Ttl.nxy[7]	=0;	/* Junk */
  Ttl.nxy[8]	=NSIZE;	/* Junk */
  Ttl.nxy[9]	=NSIZE;	/* Junk */

  Ttl.axy[0]	=dgt;
  Ttl.axy[1]	=dgy;
  Ttl.axy[2]	=EZcgp;
  Ttl.axy[3]	=ESRext;
  Ttl.axy[4]	=ESRBt/ESRext;
  for(i=5; i < 10; i++){
    Ttl.axy[i]	=0.;
  }
  fwrite((void*)&Ttl,(size_t)8,(size_t)49,Fp);
  fwrite((void*)&Ttl,(size_t)8,(size_t)49,fp);

  a	=0.;
  i	=0;
  A[i]	=a;
  sp[i]	=ESsp[0]*p2rBt;
  Ps[i]	=-ESPs[0]*rBt;
  sq[i]	=1./ESgm[0];
  if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
    dsqgy[i]	=ESgm2a[0]*sq[i]*sq[i]/(ESdgY[0]*rBt);
  }
  else{
    ESSetSplDCr(a);
    splRDCr2(&s,NULL,dsqgy,7);
    dsqgy[i]	*=sq[i]*sq[i]/(ESdgY[0]*rBt);
  }
  g[i]	=ESaF[0]*rg;
  s	=ESaT[0]*rBt;
  dggy[i]=-s/(ESRext*ESRext*g[i]);
  fb[i]	=-ESaR0[0]*ESdgY[0]*rBt/ESLc[0];
  dfbgy[i]	=ESaR0[0]*(ESdgY2a[0]-ESdgY[0]*ESLc2[0]/ESLc[0])/
    (ESdgY[0]*ESLc[0]);
  for(i=1; i < nA1; i++){
    ss	=sqrt(dgy*i);
    do{
      ESSetSplA(a);
      splRA(&s,&dgya,ESqgY,ESqgY2a);
      a		+=(ss-s)/dgya;
    }while(fabs(ss-s) > 1e-8);
    A[i]	=a;
    ESSetSplA(a);
    if(ESEqSolvInPr != 2){
      splRA(sp+i,NULL,ESsp,ESsp2a);
    }
    else{
      ESSetSplDPr(a);
      splRDPr(sp+i,&ds,2);
    }
    sp[i]	*=p2rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,0);
      Ps[i]	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,1);
      break;
    default:
      splRA(Ps+i,NULL,ESPs,ESPs2a);
      break;
    }
    Ps[i]	*=-rBt;
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    splRA(&s,&ss,ESLc,ESLc2);
    fb[i]	=-ESaR0[0]*dgya/s;
    dgya	*=a;
    dfbgy[i]	=(ESaR0[0]*ds+fb[i]*ss)/(dgya*s);
    fb[i]	*=rBt;
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq[i]	=1./s;
    dsqgy[i]	=ds/(dgya*rBt*s*s);
    splRA(g+i,&ds,ESaF,ESaF2a);
    g[i]	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    dggy[i]	=-s/(ESRext*ESRext*g[i]);
  }
  fwrite((void*)sp,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)Ps,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)sq,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dsqgy,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)g,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dggy,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)fb,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dfbgy,(size_t)8,(size_t)nA1,fp);

  a	=0.;
  ss	=sqrt(dgy*0.01);
  do{
    ESSetSplA(a);
    splRA(&s,&dgya,ESqgY,ESqgY2a);
    a		+=(ss-s)/dgya;
  }while(fabs(ss-s) > 1e-8);
  A[0]	=a;
  
  printf("Writing Hamada metric tensor\n");
  
  ji	=0;
  for(i=0; i < nA1; i++){
    a	=A[i];
    ESSetSplA(a);
    splRA(&dgya,&dD0,ESdgY,ESdgY2a);
    rf	=-1./(dgya*a*rBt);
    splRA(Lc,dLc,ESLc,ESLc2);
    L0		=2./Lc[0];
    s		=rf*L0;
    dLc[0]	*=rf/Lc[0];
    ki	=0;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      splRA(Lc+k,dLc+k,ESLc+ki,ESLc2+ki);
      Lc[k]	*=L0;
      dLc[k]	=dLc[k]*s-Lc[k]*dLc[0];
      splRA(Ls+k,dLs+k,ESLs+ki,ESLs2+ki);
      Ls[k]	*=L0;
      dLs[k]	=dLs[k]*s-Ls[k]*dLc[0];
    }

    splRA(Jc,dJc,rDc,rDc2);
    s	=-ESaR0[0]/(dgya*rBt);
    D0	=Jc[0]*s;
    dD0	=(dJc[0]-Jc[0]*dD0/dgya)*s*rf;

    L0		=2./Jc[0];
    s		=rf*L0;
    dJc[0]	*=rf/Jc[0];
    ki	=0;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      splRA(Jc+k,dJc+k,rDc+ki,rDc2+ki);
      Jc[k]	*=L0;
      dJc[k]	=dJc[k]*s-Jc[k]*dJc[0];
      splRA(Js+k,dJs+k,rDs+ki,rDs2+ki);
      Js[k]	*=L0;
      dJs[k]	=dJs[k]*s-Js[k]*dJc[0];
    }

    splRA(&b,&db,ESsb,ESsb2a);
    db	*=rf;
    splRA(rs,drs,rsT,rsT2a);
    drs[0]	*=rf;
    splRA(rc,drc,rcT,rcT2a);
    rc[0]	+=ESaR0[0];
    drc[0]	*=rf;
    ki	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,drc+k,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      drc[k]	*=2.*rf;
      rct[k]	=-k*rc[k];
      d2ratc[k]	=-k*drc[k];
      d2rttc[k]	=k*rct[k];
      splRA(rs+k,drs+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      drs[k]	*=2.*rf;
      rst[k]	=k*rs[k];
      d2rats[k]	=k*drs[k];
      d2rtts[k]	=-k*rst[k];
    }

    for(j=0; j < ESNp; j++){
      s		=0.;
      kj	=0;
      for(k=1; k < ES2Mp1; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	s	+=(Jc[k]*ESsn1[kj]+Js[k]*(1.-EScs1[kj]))/k;
      }
      gQ[j]	=s;
    }
    gQ[j]	=gQ[0];
    splP(gQ,gQ2);
    t	=EZc2gp;
    for(j=0; j < nP; j++){
      gq	=EZc2gp-dgt*j;
      do{
	ESSetSplP(t);
	splRP(&s,&ds,gQ,gQ2);
	s	+=t-gq;
	t	-=s/(1.+ds);
      }while(fabs(s) > 1e-8);
      s		=t;
      cs	=cos(s);
      sn	=sin(s);
      ss	=rc[0]+rc[1]*cs+rs[1]*sn;
      EZdra	=drc[0]+drc[1]*cs+drs[1]*sn;
      EZdrt	=rct[1]*sn+rst[1]*cs;

      gqa	=dJc[1]*sn+dJs[1]*(1.-cs);
      gqt	=1.+Jc[1]*cs+Js[1]*sn;
      gqat	=dJc[1]*cs+dJs[1]*sn;
      gqtt	=-Jc[1]*sn+Js[1]*cs;

      gh	=(Lc[1]-Jc[1])*sn+(Ls[1]-Js[1])*(1.-cs);
      gha	=dLc[1]*sn+dLs[1]*(1.-cs);
      ght	=1.+Lc[1]*cs+Ls[1]*sn;
      ghtt	=-Lc[1]*sn+Ls[1]*cs;
      ghat	=dLc[1]*cs+dLs[1]*sn;

      z[ji]	=rs[0]+b*sn;
      EZdza	=drs[0]+db*sn;
      EZdzt	=b*cs;
      d2rat	=d2ratc[1]*sn+d2rats[1]*cs;
      d2rtt	=d2rttc[1]*cs+d2rtts[1]*sn;
      d2zat	=db*cs;
      d2ztt	=-b*sn;
      for(k=2; k < ESFp1; k++){
	s	=k*t;
	cs	=cos(s);
	sn	=sin(s);
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
	d2rat	+=d2ratc[k]*sn+d2rats[k]*cs;
	d2rtt	+=d2rttc[k]*cs+d2rtts[k]*sn;
	if(k <  ES2Mp1){
	  gqa	+=(dJc[k]*sn+dJs[k]*(1.-cs))/k;
	  gqt	+=Jc[k]*cs+Js[k]*sn;
	  gqat	+=dJc[k]*cs+dJs[k]*sn;
	  gqtt	+=k*(-Jc[k]*sn+Js[k]*cs);

	  gh	+=((Lc[k]-Jc[k])*sn+(Ls[k]-Js[k])*(1.-cs))/k;
	  gha	+=(dLc[k]*sn+dLs[k]*(1.-cs))/k;
	  ght	+=Lc[k]*cs+Ls[k]*sn;
	  ghtt	+=k*(-Lc[k]*sn+Ls[k]*cs);
	  ghat	+=dLc[k]*cs+dLs[k]*sn;
	}
      }
      gd[ji]	=-gh;
      r[ji]	=ss;
      D		=EZdra*EZdzt-EZdrt*EZdza;
      gh	=1./D;
      dD	=(d2rat*EZdzt+EZdra*d2ztt-d2rtt*EZdza-EZdrt*d2zat)*gh;
#ifdef H
      gh	*=gh;
#endif
      gh	*=ss/(D0*gqt);

      g22[ji]	=(EZdrt*EZdrt+EZdzt*EZdzt)*gh;
      dg12	=(d2rat*EZdrt+EZdra*d2rtt+d2zat*EZdzt+EZdza*d2ztt)*gh;
      dg22	=2.*(d2rtt*EZdrt+d2ztt*EZdzt)*gh;

      gh	*=(EZdra*EZdrt+EZdza*EZdzt);

      g12[ji]	=gh*gqt-g22[ji]*gqa;

      s		=1./gqt;
      dg12t[ji]	=(2.*g12[ji]*dD-dg12*gqt-gh*gqtt+dg22*gqa+g22[ji]*gqat)*s;
      dg22t[ji]	=(2.*g22[ji]*dD-dg22)*s;

      rt[ji]	=-EZdrt*s;
      zt[ji]	=-EZdzt*s;

      p2r[ji]	=ss*ss;
      dp2rt[ji]	=2.*r[ji]*rt[ji];
      gqa	*=s;
      EZdra	-=EZdrt*gqa;
      EZdza	-=EZdzt*gqa;
      ra[ji]	=EZdra;
      za[ji]	=EZdza;

      dgda[ji]	=-sq[i]*(gha-ght*gqa)+dsqgy[i]*gd[ji];

      D		*=s;
#ifdef H
      g11[ji]	=(EZdra*EZdra+EZdza*EZdza)/(D*D);
#endif
      g11[ji]	=(EZdra*EZdra+EZdza*EZdza)*ss/(D*D0);
      dp2ra[ji]	=2.*ss*EZdra;

      qg[ji]	=D0;
      dqga[ji]	=dD0;
      ji++;
    }
  }
  free(gQ);
  {
    double *p,*P,d[3]={0.,0.,0.};

    ji	=nP*nA;
    /* $r_b$ into Esc2PstGm*/
    P	=r+ji;
    fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
    fwrite((void*)P,(size_t)8,(size_t)2,Fp);
    /* $z_b$ into Esc2PstGm*/
    P	=z+ji;
    fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
    fwrite((void*)P,(size_t)8,(size_t)2,Fp);
    /* $\bgY_{PEST}$ into Esc2PstGm*/
    s	=ESgY[ESNa]*EZc2gp*rBt;
    ds	=-s*dgy;
    for(i=0; i < nA1; i++){
      A[i]	=s+ds*i;
    }
    fwrite((void*)A,(size_t)8,(size_t)nA1,Fp);
    /* $r$ and $|\na\bgY|^2$*/
    /* $z$ and $r^2$*/
    /* $r'_{\gq_{PEST}}$ and $|\na\gq_{PEST}|^2$*/
    /* $z'_{\gq_{PEST}}$ and $(\na\bgY_{PEST}\cdot\na\gq_{PEST})$ */
    /* $r'_{\bgY_{PEST}}$ and $\R{\pa(\na\bgY\cdot\na\gq)}{\pa\gq}$*/
    /* $z'_{\bgY_{PEST}}$ and $\R{\pa r^2}{\pa a}$ */
    for(k=0; k < 6; k++){
      p	=w1[k];
      P	=w2[k];
      for(i=0; i < nA1; i++){
	fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
	fwrite((void*)P,(size_t)8,(size_t)2,Fp);
	fwrite((void*)d,(size_t)8,(size_t)3,Fp);
	P	+=nP;
	fwrite((void*)p,(size_t)8,(size_t)nP,fp);
	fwrite((void*)p,(size_t)8,(size_t)2,fp);
	fwrite((void*)d,(size_t)8,(size_t)3,fp);
	p	+=nP;
      }
      if(k == 0){
	for(j=0; j < nP5; j++){
	  r[j]	=0.;
	}
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,Fp);
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }

    /* $\R{\pa|\na\bgY|^2}{\pa\gq}$ */
    /* $\R{\pa r^2}{\pa\gq}$ */
    /* $\q{g_{PEST}}$ */
    /* $\R{\pa\q{g_{PEST}}}{\pa a}$ */
    /* $\gd_{PEST}$ */
    /* $\R{\pa\gd_{PEST}}{\pa a}$ */
    for(k=0; k < 6; k++){
      P	=w3[k];
      for(i=0; i < nA1; i++){
	fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
	fwrite((void*)P,(size_t)8,(size_t)2,Fp);
	fwrite((void*)d,(size_t)8,(size_t)3,Fp);
	P	+=nP;
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,Fp);
    }
    for(k=0; k < 6; k++){
      free(w3[k]);
      free(w2[k]);
      free(w1[k]);
    }
  }
  printf("ES4PstGm%4.2f (fort.31) has been written\n",ttime);
  fclose(Fp);
  printf("ES4PstPr%4.2f (fort.32) has been written\n",ttime);
  fclose(fp);
  return(0);
}
#endif

#ifdef XWIN
int ESCheckBasicFunctionsA(int i0,int i1,int j0,int j1,int na1,int np1
			   ,int F	/* choice of function */
			   ,int Fl	/* PostScript Flag */
			   )
{
  int i,j,ji,n;
  double d,s,x[na1],y[na1];
  char ln[64];
  double *ld,*r,*ra,*rq,*za,*zq,da,dt;
  
  n	=np1*na1;
  switch(F/10){
  case 0:
    ld	=R4;
    strcpy(ln," r");
    break;
  case 1:
    ld	=Z4;
    strcpy(ln," z");
    break;
  case 2:
    ld	=B4;
    strcpy(ln," B");
    break;
  case 3:
    ld	=H4;
    strcpy(ln,"gh");
    break;
  case 4:
    ld	=G4;
    strcpy(ln,"D/r");
    r	=R4;
    ra	=r+n;
    rq	=ra+n;
    za	=Z4+n;
    zq	=za+n;
    for(i=0; i < na1; i++){
      s	=0.;
      ji=np1*i;
      for(j=0; j < np1; j++){
	d	=zq[ji]*ra[ji]-za[ji]*rq[ji];
	G4[ji]	=d/r[ji];

	G4[ji]	=B4[ji]*B4[ji]+2.*ESsp[i];
	G4[ji]	=B4[ji]*B4[ji]+2.*ESsp[i];

	if(j){
	  s	+=G4[ji];
	}
	G4[ji+n]=d*r[ji];
	G4[ji+2*n]=d/r[ji];
	G4[ji+3*n]=d*r[ji];
	ji++;
      }
#ifdef H
      s		/=(np1-1);
      ji	=np1*i;

      for(j=0; j < np1; j++){
	G4[ji]	-=s;
	ji++;
      }
#endif
    }
    break;
  }
  switch(F%10){
  case 0:
    break;
  case 1:
    ld	+=n;
    strcpy(ln+2,"'_a");
    break;
  case 2:
    ld	+=2*n;
    strcpy(ln+2,"'_gt");
    break;
  case 3:
    ld	+=3*n;
    strcpy(ln+2,"''_{a,gt}");
    break;
  }

  if(i0 < 0){
    i0	=0;
  }
  if(i1 > na1){
    i1	=na1;
  }
  if(j0 < 0){
    j0	=0;
  }
  if(j1 > np1){
    j1	=np1;
  }

  da	=(ESsa[ESNa]-ESsa[0])/(na1-1);
  dt	=EZc2gp/(np1-1);
  x[0]	=da*i0;
  x[1]	=da*(i1-1);
  j	=np1*i0+j0;
  y[0]	=ld[j];
  y[1]	=y[0];
  for(i=i0; i < i1; i++){
    ji	=np1*i+j0;
    for(j=j0; j < j1; j++){
      s	=ld[ji];
      if(y[0] > s){
	y[0]	=s;
      }
      if(y[1] < s){
	y[1]	=s;
      }
      ji++;
    }
  }	
  
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
  sprintf(ln+strlen(ln)," i0=%d i1=%d j0=%d j1=%d",i0,i1,j0,j1);
  SetPlotName("a",ln,"Basic Functions");
  Scale2d(4,x,y,2,6,2);
  if(Fl){
    PSNewFrame();
  }
  for(j=j0; j < j1; j++){
    for(i=i0; i < i1; i++){
      ji	=np1*i+j;
      x[i]	=da*i;
      y[i]	=ld[ji];
    }
    Plot2d(4,x+i0,y+i0,i1-i0,6,0,14,0);
    if(Fl){
      PSPlot2d(x+i0,y+i0,i1-i0,6,0,14,0);
    }
  }
  if(Fl){
    PSFileClose();
  }
  return(0);
}

int ESCheckBasicFunctionsP(int i0,int i1,int j0,int j1,int na1,int np1
			   ,int F	/* choice of function */
			   ,int Fl	/* PostScript Flag */
			   )
{
  int i,j,ji,n;
  double d,s,x[np1],y[np1];
  char ln[64];
  double *ld,*r,*ra,*rq,*za,*zq,da,dt;
  
  n	=np1*na1;
  switch(F/10){
  case 0:
    ld	=R4;
    strcpy(ln," r");
    break;
  case 1:
    ld	=Z4;
    strcpy(ln," z");
    break;
  case 2:
    ld	=B4;
    strcpy(ln," B");
    break;
  case 3:
    ld	=H4;
    strcpy(ln,"gh");
    break;
  case 4:
    ld	=G4;
    strcpy(ln,"D/r");
    r	=R4;
    ra	=r+n;
    rq	=ra+n;
    za	=Z4+n;
    zq	=za+n;
    for(i=0; i < na1; i++){
      s	=0.;
      ji	=np1*i;
      for(j=0; j < np1; j++){
	d	=zq[ji]*ra[ji]-za[ji]*rq[ji];
	G4[ji]	=d/r[ji];
	if(j){
	  s	+=G4[ji];
	}

	G4[ji]	=B4[ji]*B4[ji]+2.*ESsp[i];

	G4[ji+n]=d*r[ji];
	G4[ji+2*n]=d/r[ji];
	G4[ji+3*n]=d*r[ji];
	ji++;
      }
#ifdef H
      s		/=(np1-1);
      ji	=np1*i;
      for(j=0; j < np1; j++){
	G4[ji]	-=s;
	ji++;
      }
#endif
    }
    break;
  }
  switch(F%10){
  case 0:
    break;
  case 1:
    ld	+=n;
    strcpy(ln+strlen(ln),"'_a");
    break;
  case 2:
    ld	+=2*n;
    strcpy(ln+strlen(ln),"'_gt");
    break;
  case 3:
    ld	+=3*n;
    strcpy(ln+strlen(ln),"''_{a,gt}");
    break;
  }

  if(i0 < 0){
    i0	=0;
  }
  if(i1 > na1){
    i1	=na1;
  }
  if(j0 < 0){
    j0	=0;
  }
  if(j1 > np1){
    j1	=np1;
  }

  da	=(ESsa[ESNa]-ESsa[0])/(na1-1);
  dt	=EZc2gp/(np1-1);
  x[0]	=dt*j0;
  x[1]	=dt*(j1-1);
  j	=np1*i0+j0;
  y[0]	=ld[j];
  y[1]	=y[0];
  for(i=i0; i < i1; i++){
    ji	=np1*i+j0;
    for(j=j0; j < j1; j++){
      s	=ld[ji];
      if(y[0] > s){
	y[0]	=s;
      }
      if(y[1] < s){
	y[1]	=s;
      }
      ji++;
    }
  }	
  
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
  sprintf(ln+strlen(ln)," i0=%d i1=%d j0=%d j1=%d",i0,i1,j0,j1);
  SetPlotName("gt",ln,"Basic Functions");
  Scale2d(3,x,y,2,6,2);
  if(Fl){
    PSNewFrame();
  }
  
  for(i=i0; i < i1; i++){
    ji	=np1*i;
    for(j=j0; j < j1; j++){
      x[j]	=dt*j;
      y[j]	=ld[ji];
      ji++;
    }
    Plot2d(3,x+j0,y+j0,j1-j0,6,0,14,0);
    if(Fl){
      PSPlot2d(x+j0,y+j0,j1-j0,6,0,14,0);
    }
  }
  if(Fl){
    PSFileClose();
  }
  return(0);
}
#ifdef SHMEM

static double *SMsr,*SMsra,*SMsrt,*SMsrat,*SMsz,*SMsza,*SMszt,*SMszat
,*SMaB,*SMaBa,*SMaBt,*SMaBat,*SMgH,*SMgHa,*SMgHt,*SMgHat;

int ESBasicFunctions2ShMem()
{
  int i,j,ji,k,ki,kj,iK;
  double r,ra,rt,rat,raa,rtt,raat,ratt,raatt;
  double z,za,zt,zat,zaa,ztt,zaat,zatt,zaatt;
  double rc[ESFp1],rs[ESFp1];
  double rct[ESFp1],rst[ESFp1],rca[ESFp1],rsa[ESFp1],rcaa[ESFp1],rsaa[ESFp1];
  double b,ba,baa,dgY,d2gY,pF,dpF,gH,gHa;
  double D,Da,Dt,Dat,K,Ka,Kt,Kat,G,Ga,Gt,Gat;
  double cs,sn;

  pF	=ESaF[0];
  G	=1./ESaR0[0];
  b	=ESsb[0];
  ba	=ESsb1a[0];
  baa	=ESsb2a[0];
  rcaa[0]	=rcT2a[0];
  rsaa[0]	=rsT2a[0];
  rca[1]	=2.*rcT1a[ESNa1];
  rsa[1]	=2.*rsT1a[ESNa1];
  ki	=0;
  iK	=2;
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    rcaa[k]	=2.*rcT2a[ki];
    rsaa[k]	=2.*rsT2a[ki];
    if(rcaa[k] != 0. || rsaa[k] != 0.){
      iK	=k+1;
    }
  }
  gH	=0.;
  gHa	=0.;
  for(j=0; j < ESNp; j++){
    cs		=EScs1[j];
    sn		=ESsn1[j];
    ra		=rca[1]*cs+rsa[1]*sn;
    rat		=rsa[1]*cs-rca[1]*sn;
    ratt	=-ra;
    raa		=rcaa[0]+rcaa[1]*cs+rsaa[1]*sn;
    raat	=rsaa[1]*cs-rcaa[1]*sn;
    raatt	=-rcaa[1]*cs-rsaa[1]*sn;
    za		=ba*sn;
    zat		=ba*cs;
    zatt	=-za;
    zaa		=rsaa[0]+baa*sn;
    zaat	=baa*cs;
    zaatt	=-baa*sn;
    kj	=j < ESNp ? j : j-ESNp;
    for(k=2; k < iK; k++){
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      cs	=EScs1[kj];
      sn	=ESsn1[kj];
      raa	+=rcaa[k]*cs+rsaa[k]*sn;
      raat	+=k*(rsaa[k]*cs-rcaa[k]*sn);
      raatt	-=k*k*(rcaa[k]*cs+rsaa[k]*sn);
    }

    SMsr[j]	=ESaR0[0];
    SMsra[j]	=ra;
    SMsrt[j]	=0.;
    SMsrat[j]	=rat;

    SMsz[j]	=ESaZ0[0];
    SMsza[j]	=za;
    SMszt[j]	=0.;
    SMszat[j]	=zat;

    Ga	=-ra*G*G;
    Gt	=0.;
    Gat	=-rat*G*G;
    SMaB[j]	=pF*G;
    SMaBa[j]	=pF*Ga;
    SMaBt[j]	=0.;
    SMaBat[j]	=pF*Gat;

    SMgH[j]	=0.;
    SMgHt[j]	=0.;

    D	=ra*zat-za*rat;

    SMgHa[j]	=(raa*zat-zaa*rat+EZcr2*(ra*zaat-za*raat))/D-ra*G;
    SMgHat[j]	=(raa*zatt-zaa*ratt+EZcr2*(ra*zaatt+zat*raat-rat*zaat-za*raatt))
      /D-rat*G;
    if(j){
      gHa	+=SMgHa[j];
      gH	+=SMgHat[j];
    }
  }
  gHa	/=ESNp;
  gH	/=ESNp;
  for(j=0; j < ESNp; j++){
    SMgHa[j]	-=gHa;
    SMgHat[j]	-=gH;
  }
  SMsr[j]	=SMsr[0];
  SMsra[j]	=SMsra[0];
  SMsrt[j]	=SMsrt[0];
  SMsrat[j]	=SMsrat[0];
  
  SMsz[j]	=SMsz[0];
  SMsza[j]	=SMsza[0];
  SMszt[j]	=SMszt[0];
  SMszat[j]	=SMszat[0];

  SMaB[j]	=SMaB[0];
  SMaBa[j]	=SMaBa[0];
  SMaBt[j]	=SMaBt[0];
  SMaBat[j]	=SMaBat[0];

  SMgH[j]	=SMgH[0];
  SMgHa[j]	=SMgHa[0];
  SMgHt[j]	=SMgHt[0];
  SMgHat[j]	=SMgHat[0];

  ji	=ESNp1;
  for(i=1; i < ESNa1; i++){
    dgY		=ESsa[i]*ESdgY[i];
    d2gY	=ESdgY[i]+ESsa[i]*ESdgY1a[i];
    pF		=ESFF[i];
    dpF		=ESFF1a[i];
    rc[0]	=ESaR0[0]+rcT[i];
    rs[0]	=ESaZ0[0]+rsT[i];
    rca[0]	=rcT1a[i];
    rsa[0]	=rsT1a[i];
    rcaa[0]	=rcT2a[i];
    rsaa[0]	=rsT2a[i];
    b		=ESsb[i];
    ba		=ESsb1a[i];
    baa		=ESsb2a[i];
    ki	=i;
    iK	=1;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=2.*rcT[ki];
      rct[k]	=2.*k*rsT[ki];
      rca[k]	=2.*rcT1a[ki];
      rcaa[k]	=2.*rcT2a[ki];
      rs[k]	=2.*rsT[ki];
      rst[k]	=-2.*k*rcT[ki];
      rsa[k]	=2.*rsT1a[ki];
      rsaa[k]	=2.*rsT2a[ki];
      if(rct[k] != 0. || rca[k] != 0. || rcaa[k] != 0. ||
	 rst[k] != 0. || rsa[k] != 0. || rsaa[k] != 0.){
	iK	=k+1;
      }
    }
    gH	=0.;
    gHa	=0.;

    for(j=0; j < ESNp; j++){
      r		=rc[0];
      ra	=rca[0];
      rt	=0.;
      rat	=0.;
      raa	=rcaa[0];
      rtt	=0.;
      raat	=0.;
      ratt	=0.;
      cs	=EScs1[j];
      sn	=ESsn1[j];
      za	=rsa[0]+ba*sn;
      zt	=b*cs;
      zat	=ba*cs;
      zaa	=rsaa[0]+baa*sn;
      ztt	=-b*sn;
      zaat	=baa*cs;
      zatt	=-ba*sn;
      kj	=0;
      for(k=1; k < iK; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	cs	=EScs1[kj];
	sn	=ESsn1[kj];
	r	+=rc[k]*cs+rs[k]*sn;
	ra	+=rca[k]*cs+rsa[k]*sn;
	rt	+=rct[k]*cs+rst[k]*sn;
	rat	+=k*(rsa[k]*cs-rca[k]*sn);
	raa	+=rcaa[k]*cs+rsaa[k]*sn;
	rtt	+=k*(rst[k]*cs-rct[k]*sn);
	raat	+=k*(rsaa[k]*cs-rcaa[k]*sn);
	ratt	-=k*k*(rca[k]*cs+rsa[k]*sn);
      }
      SMsz[ji]	=rs[0]+b*ESsn1[j];
      SMsza[ji]	=za;
      SMszt[ji]	=zt;
      SMszat[ji]=zat;

      SMsr[ji]	=r;
      SMsra[ji]	=ra;
      SMsrt[ji]	=rt;
      SMsrat[ji]=rat;

      D		=ra*zt-rt*za;
      Da	=raa*zt+ra*zat-rat*za-rt*zaa;
      Dt	=rat*zt+ra*ztt-rtt*za-rt*zat;
      Dat	=raat*zt+raa*ztt+ra*zatt-ratt*za-rtt*zaa-rt*zaat;

      K		=1./r;
      Ka	=-ra*K*K;
      Kt	=-rt*K*K;
      Kat	=(2.*rt*ra*K-rat)*K*K;
      
      SMgH[ji]	=D*K;
      SMgHa[ji]	=Da*K+D*Ka;
      SMgHt[ji]	=Dt*K+D*Kt;
      SMgHat[ji]=Dat*K+Dt*Ka+Da*Kt+D*Kat;
      gH	+=SMgH[ji];
      gHa	+=SMgHa[ji];

      SMaB[ji]	=pF*K*K;
      SMaBa[ji]	=(dpF*K+2.*pF*Ka)*K;
      SMaBt[ji]	=2.*pF*Kt*K;
      SMaBat[ji]=2.*dpF*K*Kt+2.*pF*(Kat*K+Kt*Ka);

#ifdef ForFRC
      SMaB[ji]	=0.;
      SMaBa[ji]	=0.;
      SMaBt[ji]	=0.;
      SMaBat[ji]=0.;
#endif 

      D		=1./D;
      Da	*=D;
      Dt	*=D;
      Dat	=(2.*Dt*Da-Dat*D)*D;
      Da	*=-D;
      Dt	*=-D;

      G		=K*D;
      Ga	=Ka*D+K*Da;
      Gt	=Kt*D+K*Dt;
      Gat	=Kat*D+Kt*Da+Ka*Dt+K*Dat;

      D		=rt*rt+zt*zt;
      Da	=2.*(rat*rt+zat*zt);
      Dt	=2.*(rtt*rt+ztt*zt);
      Dat	=2.*(ratt*rt+rtt*rat+zatt*zt+ztt*zat);

      K		=D*G*G;
      Ka	=(Da*G+2.*D*Ga)*G;
      Kt	=(Dt*G+2.*D*Gt)*G;
      Kat	=(Dat*G+2.*Dt*Ga+2.*Da*Gt+2.*D*Gat)*G+2.*D*Gt*Ga;
      SMaB[ji]	+=K*dgY*dgY;
      SMaBa[ji]	+=(Ka*dgY+2.*K*d2gY)*dgY;
      SMaBt[ji]	+=Kt*dgY*dgY;
      SMaBat[ji]=(Kat*dgY+2.*Kt*d2gY)*dgY;

      D	=sqrt(SMaB[ji]);
      SMaB[ji]	=D;
      D		=EZcr2/D;
      SMaBa[ji]	*=D;
      SMaBt[ji]	*=D;
      SMaBat[ji]=(SMaBat[ji]-2.*SMaBa[ji]*SMaBt[ji])*D;
      ji++;
    }
    G	=1./ESNp;
    gH	*=G;
    gHa	*=G;
    j	=ji-ESNp;
    SMsr[ji]	=SMsr[j];
    SMsra[ji]	=SMsra[j];
    SMsrt[ji]	=SMsrt[j];
    SMsrat[ji]	=SMsrat[j];
    SMsz[ji]	=SMsz[j];
    SMsza[ji]	=SMsza[j];
    SMszt[ji]	=SMszt[j];
    SMszat[ji]	=SMszat[j];
    SMaB[ji]	=SMaB[j];
    SMaBa[ji]	=SMaBa[j];
    SMaBt[ji]	=SMaBt[j];
    SMaBat[ji]	=SMaBat[j];
    SMgH[ji]	=SMgH[j];
    SMgHa[ji]	=SMgHa[j];
    SMgHt[ji]	=SMgHt[j];
    SMgHat[ji]	=SMgHat[j];
    ji++;
    G	=1./gH;
    while(j < ji){
      SMgHa[j]	=(SMgHa[j]-gHa*SMgH[j]*G)*G;
      SMgH[j]	=(SMgH[j]-gH)*G;
      j++;
    }
  }
  return(0);
}

int NewESWriteEqOutShM()
{
  static int n=0;
  extern int *ESShMemI;
  extern double *ESShMemE,*ESShMemP,*ESShMemT;
  double rr0;
  int i,kA;
  char *lc;
  double *ld;
  
  kA	=ESNa1*sizeof(double);

  ESShMemI[7]	=ESMp;
  ESShMemI[9]	=0;
  ESShMemE[0]	=ESaR0[0];
  ESShMemE[1]	=ESaZ0[0];
  lc	=(char*)(ESShMemE+2);
  memcpy(lc,(void*)ESsa,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESsb,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESsb2a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESgY,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESgY2a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESgF,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESgF2a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESaF,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESaF2a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESgm,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESgm2a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESsq,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESsq2a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESsp,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESsp2a,kA);
  lc	+=kA;
  i	=ESNp1*kA;
  memcpy(lc,(void*)ESsr,i);
  lc	+=i;
  memcpy(lc,(void*)ESsz,i);
  lc	+=i;
  i	=kA*ESFp1;
  memcpy(lc,(void*)rcT,i);
  lc	+=i;
  memcpy(lc,(void*)rsT,i);
  lc	+=i;
  memcpy(lc,(void*)rcT2a,i);
  lc	+=i;
  memcpy(lc,(void*)rsT2a,i);
  lc	+=i;

  i	=ESNp1*sizeof(double);
  memcpy(lc,(void*)ESgt,i);
  lc	+=i;
  memcpy(lc,(void*)ESsa,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESaF,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESaF1a,kA);
  lc	+=kA;
  rr0	=1./ESaR0[0];
  ld	=(double*)lc;
  for(i=0; i < ESNa1; i++){
    *ld++	=ESLc[i]*ESaF[i]*rr0;
  }
  for(i=0; i < ESNa1; i++){
    *ld++	=(ESLc1[i]*ESaF[i]+ESLc[i]*ESaF1a[i])*rr0;
  }
  lc	=(void*)ld;
  memcpy(lc,(void*)ESdgY,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESdgY1a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESaT,kA);
  lc	+=kA;
  ld	=(double*)lc;
  splA1a(ESaT,ld,ESaT2a);
  lc	+=kA;
  memcpy(lc,(void*)ESPs,kA);
  lc	+=kA;
  ld	=(double*)lc;
  splA1a(ESPs,ld,ESPs2a);
  lc	+=kA;
  i	=ESNp1*ESNa1;
  SMsr	=(double*)lc;
  SMsra	=SMsr	+i;
  SMsrt	=SMsra	+i;
  SMsrat=SMsrt	+i;
  SMsz	=SMsrat	+i;
  SMsza	=SMsz	+i;
  SMszt	=SMsza	+i;
  SMszat=SMszt	+i;
  SMaB	=SMszat	+i;
  SMaBa	=SMaB	+i;
  SMaBt	=SMaBa	+i;
  SMaBat=SMaBt	+i;
  SMgH	=SMaBat	+i+2*i;
  SMgHa	=SMgH	+i;
  SMgHt	=SMgHa	+i;
  SMgHat=SMgHt	+i;
  ESBasicFunctions2ShMem();
  return(0);
}

int ESWriteEqOutShM()
{
  static int n=0;
  extern int *ESShMemI;
  extern double *ESShMemE,*ESShMemP,*ESShMemT;
  int i,kA;
  double r,rR0,R0;
  int j,k,ki,kj;
  double b,ba,ch2hp,crHp,cs,sn;
  double rc[ESFp1],rca[ESFp1],rs[ESFp1],rsa[ESFp1]
    ,as[ESNp1],as1[ESNp1],as2p[ESNp1],as12p[ESNp1];
  double ra,rgt,za,zgt;
  char *lc;
  double *ld;
  double *a,*F,*gm,*g22,*g33,*D2,*D3,*S,*S1;
  int nP,M1;

  a	=ESShMemP;
  nP	=ESShMemI[6];

  ch2hp	=EZcr12*ESgt[1]*ESgt[1];
  crHp	=1./ESNp;
  
  kA	=ESNa1*sizeof(double);
  R0	=ESaR0[0];
  rR0	=1./R0;

  ESShMemI[7]	=ESMp;
  ESShMemI[9]	=0;
  ESShMemE[0]	=ESaR0[0];
  ESShMemE[1]	=ESaZ0[0];
  lc	=(void*)(ESShMemE+2);
  memcpy(lc,(void*)ESsa,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESsb,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESsb2a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESgY,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESgY2a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESgF,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESgF2a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESaF,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESaF2a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESgm,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESgm2a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESsq,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESsq2a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESsp,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESsp2a,kA);
  lc	+=kA;
  i	=ESNp1*kA;
  memcpy(lc,(void*)ESsr,i);
  lc	+=i;
  memcpy(lc,(void*)ESsz,i);
  lc	+=i;
  i	=kA*ESFp1;
  memcpy(lc,(void*)rcT,i);
  lc	+=i;
  memcpy(lc,(void*)rsT,i);
  lc	+=i;
  memcpy(lc,(void*)rcT2a,i);
  lc	+=i;
  memcpy(lc,(void*)rsT2a,i);
  lc	+=i;

/*****/
  i	=ESNp1*sizeof(double);
  memcpy(lc,(void*)ESgt,i);
  lc	+=i;
  memcpy(lc,(void*)ESsa,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESaF,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESaF1a,kA);
  lc	+=kA;
  ld	=(double*)lc;
  for(i=0; i < ESNa1; i++){
    *ld++	=ESLc[i]*ESaF[i]*rR0;
  }
  for(i=0; i < ESNa1; i++){
    *ld++	=(ESLc1[i]*ESaF[i]+ESLc[i]*ESaF1a[i])*rR0;
  }
  lc	=(void*)ld;
  memcpy(lc,(void*)ESdgY,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESdgY1a,kA);
  lc	+=kA;
  memcpy(lc,(void*)ESaT,kA);
  lc	+=kA;
  ld	=(double*)lc;
  splA1a(ESaT,ld,ESaT2a);
  lc	+=kA;
  memcpy(lc,(void*)ESPs,kA);
  lc	+=kA;
  ld	=(double*)lc;
  splA1a(ESPs,ld,ESPs2a);
  lc	+=kA;
  i	=ESNp1*ESNa1;
  SMsr	=(double*)lc;
  SMsra	=SMsr	+i;
  SMsrt	=SMsra	+i;
  SMsrat=SMsrt	+i;
  SMsz	=SMsrat	+i;
  SMsza	=SMsz	+i;
  SMszt	=SMsza	+i;
  SMszat=SMszt	+i;
  SMaB	=SMszat	+i;
  SMaBa	=SMaB	+i;
  SMaBt	=SMaBa	+i;
  SMaBat=SMaBt	+i;
  SMgH	=SMaBat	+i+2*i;
  SMgHa	=SMgH	+i;
  SMgHt	=SMgHa	+i;
  SMgHat=SMgHt	+i;

  ld	=SMgHat	+i;
#ifdef ASTRA
  *ld++	=ESgF[ESNa]*EZc2gp;
  F	=ld;
  gm	=F+nP;
  g22	=gm+nP;
  g33	=g22+nP;
  D2	=g33+nP;
  D3	=D2+nP;
  S	=D3+nP;
  S1	=S+nP;
  for(i=0; i < nP; i++){
    ESSetSplA(a[i]);
    splRA(F,NULL,ESFF,ESFF2a);
    *F++	=sqrt(*F);
    splRA(gm,NULL,ESgm,ESgm2a);
    *gm++;
    splRA(g22,NULL,ESg22c,ESg22c2);
    *g22++	*=a[i];
    splRA(g33,NULL,ESLc,ESLc2);
    splRA(D3,NULL,ESVc,ESVc2);
    *D3++	=a[i]*ESaR0[0]*(*D3+*g33);
    *g33++	*=a[i]/ESaR0[0];
    splRA(D2,NULL,ESDPc,ESDPc2);
    *D2++	*=a[i];
    if(a[i] != 0.){
      splRA(&b,&ba,ESsb,ESsb2a);
      splRA(rc,rca,rcT,rcT2a);
      splRA(rs,rsa,rsT,rsT2a);
      M1	=2;
      ki	=0;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	splRA(rc+k,rca+k,rcT+ki,rcT2a+ki);
	splRA(rs+k,rsa+k,rsT+ki,rsT2a+ki);
	rca[k]	*=2.;
	rsa[k]	*=2.;
	rc[k]	*=2;
	rs[k]	*=2;
	if(rc[k] != 0. || rs[k] != 0.) M1=k+1;
      }
      for(j=0; j < ESNp; j++){
	r	=ESaR0[0]+rc[0];
	ra	=rca[0];
	rgt	=0.;
	kj	=0;
	for(k=1; k < M1; k++){
	  kj	+=j;
	  if(kj >= ESNp) kj -=ESNp;
	  cs	=EScs1[kj];
	  sn	=ESsn1[kj];
	  r	+=rc[k]*cs+rs[k]*sn;
	  ra	+=rca[k]*cs+rsa[k]*sn;
	  rgt	+=k*(-rc[k]*sn+rs[k]*cs);
	}
	zgt	=b*EScs1[j];
	za	=rsa[0]+ba*ESsn1[j];
	as[j]	=rgt*rgt+zgt*zgt;
	as1[j]	=sqrt(rgt*rgt+zgt*zgt)*r;
	as[j]	*=r/(ra*zgt-rgt*za);
      }
      splP(as,as2p);
      splP(as1,as12p);
      b	=0.;
      ba=0.;
      for(j=0; j < ESNp; j++){
#ifdef H
	b	+=as[j]-ch2hp*as2p[j];
	ba	+=as1[j]-ch2hp*as12p[j];
#endif
	b	+=as[j];
	ba	+=as1[j];
      }
      *S	=b*crHp;
      *S1	=ba*crHp;
      S++;
      S1++;
    }
    else{
      *S++	=0.;
      *S1++	=0.;
    }
  }
#endif
  ESBasicFunctions2ShMem();
  return(0);
}
#endif/*SHMEM*/

double gy00=1.,gy01=0.,gy02=0.,gy03=0.,gy10=0.,gy11=0.,gy12=0.,gy13=0.,
gy20=0.,gy21=0.,gy22=0.,gy23=0.,gy30=0.,gy31=0.,gy32=0.,gy33=0.;

int ESEquilOut()
{
  static int FlB=0,FlIap=0,FlIa=0,FlIi=0,FlIp=0,FlIj=0;
  static int Fla=0,Flp=0;
  static int kMap=0,nMap=0;

  int i,n0,n1,ntor,na1,np1;
  static int nA=400,i0=0,i1=21,j0=0,j1=65,m0=0,m1=0;
  int Fl=0,FlMap,Flw;
  double Bt,Rext,r0,EZz0,rx,zx,Rt,Zt;
  double *pr,*pz,*pvr,*pvz,rD[12],zD[12];
  double gY,gYr,gYz,gy,gyr,gyz;
  char ln[256];
  static double ttime=0.;
  static char FNm[12]="ES4Dcon00.0";

  ESGetPlVDAddr(&rdPV,&zdPV);
  Flw	=0;
  FlMap	=0;
  ntor	=1;

  i1	=ESNa1;
  m1	=ESMp1;
  Rt	=ESaR0[0];
  Zt	=ESaZ0[0]+0.5*ESsb[ESNa];
#ifndef Tbl_EqOut
  if(VmNa > 128){
    VmNa	=128;
  }
  if(VmNa < 0){
    VmNa	=10;
  }
  VmNa1	=VmNa+1;

  if(VmNp > 128){
    VmNp	=128;
  }
  if(VmNp < 12){
    VmNp	=12;
  }
  VmNp1	=VmNp+1;

  if(VmMpol > 20){
    VmMpol	=20;
  }
  if(VmMpol < 3){
    VmMpol	=3;
  }
  pr	=ESsr;
  pz	=ESsz;
  n0	=ESnAP;
  n1	=ESnAP;
  pvr	=ESsr;
  pvz	=ESsz;
  r0	=ESaR0[0];
  EZz0	=ESaZ0[0];
  rx	=ESaRx;
  zx	=ESaZx;
  for(i=0; i < 12; i++){
    rD[i]	=rdPV[i];
    zD[i]	=zdPV[i];
  }
  Bt	=ESRBt/ESRext;
  Rext	=ESRext;

  CbUserPlot=NULL;

  if(kMap){
    if(kMap == 1) ESGetMapBox(rMap,zMap);
    if(kMap < 3){
      i	=nrMap*nzMap;
      if(nMap < i){
	if(iAMap != NULL) free(iAMap);
	nMap	=i;
	iAMap=(int*)malloc(nMap*sizeof(int));
      }
      ESSetInsideInd(iAMap,rMap,zMap,nrMap,nzMap);
      ESgY2Spl(nrMap,nzMap);
    }
    kMap	=3;
    CbUserPlot=ESOutPlot2;
    Fl	=0;
  }
  switch(Fl){
  case 0:
    if(kMap && ESSetSpl4gY(Rt,Zt) == 0){
      ESSpl2gY(&gY,&gYr,&gYz);
      sprintf(ln,"gY=%13.6e %12.5e %12.5e\n",gY,gYr,gYz);
      ESGetgY(&gY,&gYr,&gYz,Rt,Zt);
      sprintf(ln+43,"gY=%13.6e %12.5e %12.5e",gY,gYr,gYz);
      CbUserMessage=ln;
      if(Flw) ESWriteESIRZ(0);
    }
    break;
  case 1:
    for(i=0; i< n0; i++){
      R3[i]	=ESsr[i];
      Z3[i]	=ESsz[i];
    }
    pvr	=R3;
    pvz	=Z3;
    break;
  case 2:
    pvr	=R3;
    pvz	=Z3;
    break;
#ifdef VMEC
  case 3:
    n1	=33*65;
    pr	=R;
    pz	=Z;
    n0	=n1;
    pvr	=R2;
    pvz	=Z2;
    ESjso2eG(R,Z);
    printf("Jsolver started\n");
    Vmt0=clock();
    ESEsc2jso(R2,Z2);
    printf("Jsotime=%10.3e sec\n",(double)(clock()-Vmt0)/CLOCKS_PER_SEC);
    break;
  case 4:
    n1	=33*65;
    pr	=R;
    pz	=Z;
    n0	=n1;
    pvr	=R2;
    pvz	=Z2;
    break;
  case 5:
    FlStart=0;
    break;
  case 6:
    if(FlVm6){
      ESInitTrGrid();
      FlVm6=0;
    }
    i	=ESnAP-ESNp1;
    VMbound(pr+i,pz+i,ESNp,R1+i,Z1+i);
    Vmt0=clock();
    Vm2EC(R1,Z1);
    printf("Vmtime=%10.3e sec\n",(double)(clock()-Vmt0)/CLOCKS_PER_SEC);
    n1	=ESnAP;
    pvr	=R1;
    pvz	=Z1;
#ifdef H
    {
      FILE *fp;
      char *NmFp="esc2vmP",*NmFb="esc2vmB",*NmFg="esc2vmG";
      int n,n1,j,m,m1,k,ki,kj,ji;
      double a,d,gF,gY,p,gm;
      double t,dt,r,z,rc[33],rs[33],cs,sn,b;

      fp	=fopen(NmFp,"w");
      if(fp == NULL){
	printf("%s cannot be written%c\n",NmFp,'\a');
	break;
      }
      n		=40;
      while(n < 100){
	fprintf(fp,"--------- Plasma profiles on %2d %s\n",n,
		"Radial intervals ---------");
	fprintf(fp,"%12s %14s %14s\n","Rref [m]"
		,"Btor ref [T]","RBtor [m x T]");
	fprintf(fp,"%12.5e %14.5e %14.5e\n",ESRext,ESRBt/ESRext,ESRBt);
	fprintf(fp,"%2s %11s %12s %12s %12s %12s\n"," #","\\Phi/\\Phi_0"
		,"\\Phi [Vsec]","\\Psi [Vsec]","p [MPa]","\\iota");
	n1	=n+1;
	d	=ESsa[ESNa]/n;
	for(i=0; i < n1; i++){
	  a	=sqrt(i*d);
	  ESSetSplA(a);
	  splRA(&gF,NULL,ESgF,ESgF2a);
	  splRA(&gY,NULL,ESgY,ESgY2a);
	  if(ESEqSolvInPr == 2){
	    ESSetSplDPr(a);
	    splRDPr(&p,NULL,2);
	  }
	  else{
	    splRA(&p,NULL,ESsp,ESsp2a);
	  }
	  if(ESEqSolvInCr == 6 || ESEqSolvInCr == 7){
	    ESSetSplDCr(a);
	    splRDCr(&gm,NULL,7);
	  }
	  else{
	    splRA(&gm,NULL,ESgm,ESgm2a);
	  }
	  fprintf(fp,"%2d %11.2f %12.5e %12.5e %12.5e %12.5e\n",i,i*d
		  ,gF*EZc2gp,gY*EZc2gp,p*EZcrgm0,gm);
	}
	n	*=2;
	fprintf(fp,"\n");
      }
      fclose(fp);
      printf("%s has been written\n",NmFp);

      fp	=fopen(NmFb,"w");
      if(fp == NULL){
	printf("%s cannot be written%c\n",NmFb,'\a');
	break;
      }
      m		=64;
      while(m < 200){
	fprintf(fp,"--------- Plasma boundary on %2d %s\n",m,
		"Poloidal intervals ---------");
	fprintf(fp,"%3s %12s %12s\n"," #","Rpl-vac","Zpl-vac");
	m1	=m+1;
	dt	=EZc2gp/m;
	ki	=ESNa;
	for(k=0; k < ESFp1; k++){
	  rc[k]	=rcT[ki];
	  rs[k]	=rsT[ki];
	  ki	+=ESNa1;
	}
	for(j=0; j < m1; j++){
	  t	=j*dt;
	  r	=0.;
	  for(k=2; k < ESFp1; k++){
	    cs	=cos(k*t);
	    sn	=sin(k*t);
	    r	+=rc[k]*cs+rs[k]*sn;
	  }
	  cs	=cos(t);
	  sn	=sin(t);
	  r	=ESaR0[0]+rc[0]+2.*(r+rc[1]*cs+rs[1]*sn);
	  z	=ESaZ0[0]+rs[0]+ESsb[ESNa]*sn;
	  fprintf(fp,"%3d %12.5e %12.5e\n",j,r,z);
	}
	m	*=2;
	fprintf(fp,"\n");
      }
      fclose(fp);
      printf("%s has been written\n",NmFb);

      fp	=fopen(NmFg,"w");
      if(fp == NULL){
	printf("%s cannot be written%c\n",NmFg,'\a');
	break;
      }
      n		=20;
      n1	=n+1;
      d	=ESsa[ESNa]/n;
      m		=64;
      m1	=m+1;
      dt	=EZc2gp/m;
      i		=0;
      ji	=0;
      fprintf(fp,"------ Plasma geometry on %2d Radial x %3d Poloidal%s\n"
	      ,n,m," intervals ------");
      while(i < n1){
	a	=sqrt(i*d);
	ESSetSplA(a);
	ki	=0;
	for(k=0; k < ESFp1; k++){
	  splRA(rc+k,NULL,rcT+ki,rcT2a+ki);
	  splRA(rs+k,NULL,rsT+ki,rsT2a+ki);
	  ki	+=ESNa1;
	}
	for(j=0; j < m1; j++){
	  t	=j*dt;
	  r	=0.;
	  for(k=2; k < ESFp1; k++){
	    cs	=cos(k*t);
	    sn	=sin(k*t);
	    r	+=rc[k]*cs+rs[k]*sn;
	  }
	  cs	=cos(t);
	  sn	=sin(t);
	  r	=ESaR0[0]+rc[0]+2.*(r+rc[1]*cs+rs[1]*sn);
	  splRA(&b,NULL,ESsb,ESsb);
	  z	=ESaZ0[0]+rs[0]+b*sn;
	  R[ji]	=r;
	  Z[ji]	=z;
	  ji++;
	}
	i++;
	if(i%3 == 0){
	  fprintf(fp,"%2s     r[%2d] %9sz     r[%2d] %9sz     r[%2d] %9sz\n"
		  ," #",i-3,"",i-2,"",i-1,"");
	  ji	-=m1;
	  kj	=ji-m1;
	  ki	=kj-m1;
	  for(j=0; j < m1; j++){
	    fprintf(fp,"%2d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",j
		    ,R[ki],Z[ki],R[kj],Z[kj],R[ji],Z[ji]);
	    ji++;
	    ki++;
	    kj++;
	  }
	}
      }
      fclose(fp);
      printf("%s has been written\n",NmFg);
    }
#endif
    break;
  case 7:
    if(FlVm6){
      ESInitTrGrid();
      FlVm6=0;
    }
    i	=ESnAP-ESNp1;
    n1	=ESnAP;
    pvr	=R1;
    pvz	=Z1;
    break;
  case 8:
    i	=ESnAP-ESNp1;
    printf("Vmec8 started\n");
    Vmt0=clock();

    vm8B(pr+i,pz+i,ESNp,R1,Z1,gL1);

    printf("Vmtime=%10.3e sec\n",(double)(clock()-Vmt0)/CLOCKS_PER_SEC);
    ESe2vm8G(R,Z);
    FlVm=1;
    n0	=65*VmNa1;
    n1	=65*VmNa1;
    pr	=R;
    pz	=Z;
    pvr	=R1;
    pvz	=Z1;
    {
      extern int FlagPS;
      double gh0,ghmx,ghmn,x[2],y[2],dy,EZz0,gh1[65],ghV[65],da;
      int j,ji,ii;
      ghmx=0.;
      ghmn=0.;
      da	=1./VmNa;
      ii	=VmNa1/40;
      for(i=ii; i < VmNa1; i +=ii){
	ji	=ESNp1*i;
	for(j=0; j < ESNp1; j++){
	  gh0	=gL1[ji];
	  if(ghmx < gh0)
	    ghmx=gh0;
	  if(ghmn > gh0)
	    ghmn=gh0;
	  ji++;
	}
      }
      x[0]	=ESgt[0];
      x[1]	=ESgt[ESNp];
      dy	=0.1*(ghmx-ghmn)/ii;
      y[0]	=ghmn;
      y[1]	=dy*VmNa+ghmx;
      Scale2D(2,2,x,y,2,6,2);
      PSFileOpen();
      PSScale2D(x,y,2,6,2);
      PSPlotFrame(6,2);
      for(i=ii; i < VmNa1; i +=ii){
	EZz0	=dy*i;
	ji	=ESNp1*i;
	ESe2vm8gL(gh1,ghV,R1+ji,Z1+ji,gL1+ji,sqrt(da*(i-0.5)));
	for(j=0; j < ESNp1; j++){
	  gh1[j]	+=EZz0;
	  ghV[j]	+=EZz0;
	}
	Plot2d(2,ESgt,gh1,ESNp1,6,0,0,0);
	PSPlot2D(ESgt,gh1,ESNp1,6,2,0,0);

	Plot2d(2,ESgt,ghV,ESNp1,6,0,4,0);
	PSPlot2D(ESgt,ghV,ESNp1,6,2,4,0);
#ifdef H
	for(j=0; j < ESNp1; j++){
	  gh1[j]=EZz0+gL1[ji];
	  ji++;
	}
	Plot2d(2,ESgt,gh1,ESNp1,6,0,4,0);
	PSPlot2D(ESgt,gh1,ESNp1,6,2,4,0);
#endif
      }
      PSPlot2D(ESgt,gh1,ESNp1,6,2,4,-1);
      CbFlush();
    }
    break;
  case 9:
    i	=ESnAP-ESNp1;
    ESe2vm8G(R,Z);
    n0	=65*VmNa1;
    n1	=65*VmNa1;
    pr	=R;
    pz	=Z;
    pvr	=R1;
    pvz	=Z1;
    break;
#endif
  case 10:
    ttime	+=0.01*0.;
#ifdef H
    ES4DconFloat(ttime);
    if((ESEqSolvFl&0x0F) == 0){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Dcon%4.2f",ttime);
      ES4Dcon(ttime);

      FlMap	=5;
      if(ESEditDconInFile(ttime,ntor) == 0){
	Vmt0=clock();
	system("/scratchu/zakh/dcon_3.10/dcon/dcon;");
	ESDconVsPst(0);
      }
    }
    else{
      printf("Switch Eq Solver to Lsode\n in section EqSolv of ESC%c\n",7);
    }
#endif

    ES4DconFloat(ttime);
    if(ESEditDconInFile(ttime,ntor) == 0){
      Vmt0=clock();
      system("/w/Dcn/dcon;");
      system("/w/Dcn/xdraw crit;");
      system("/w/Dcn/xdraw dcon;");
#ifdef H
      ESDconSaveSt(ntor);
#endif
    }
    break;
  case 11:
    ttime	+=0.01*0.;
    ES4DconFloat(ttime);
    if(ESEditDconInFile(ttime,ntor) == 0){
      Vmt0=clock();
      system("/w/Dcn/Dcon&");
#ifdef H
      system("/w/Dcn/dcon;");
      system("/w/Dcn/xdraw crit;");
      system("/w/Dcn/xdraw dcon;");
      ESDconSaveSt(ntor);
#endif
    }
    break;
  case 12:
    if((ESEqSolvFl&0x0F) == 0){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Dcon%4.2f",ttime);
      ES4DconFloat0(ttime);
    }
    else{
      printf("Switch Eq Solver to Lsode\n in section EqSolv of ESC%c\n",7);
    }
    break;
  case 13:
    if((ESEqSolvFl&0x0F) == 0){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Dcon%4.2f",ttime);
      ES4DconASCII(ttime);
    }
    else{
      printf("Switch Eq Solver to Lsode\n in section EqSolv of ESC%c\n",7);
    }
    break;
  case 14:
    ttime	+=0.01;
    sprintf(FNm,"ES4Mars%4.2f",ttime);
    ES4Mars(ttime);
    break;
  case 15:
#ifdef H
    if(VmNp%2 == 0){
      VmNp	-=1;
      VmNp1	-=1;
    }
#endif
    ES2Orbit(VmNa,VmNp);
    CbUserPlot=ESOutPlot1;
#ifdef H
    ES2eqdskBoozer(ttime);
#endif
    break;
  case 16:
#ifdef H
    if(VmNp%2 == 0){
      VmNp	-=1;
      VmNp1	-=1;
    }
#endif
#ifdef H
    ES2Boozer(VmNa,VmNp);
    CbUserPlot=ESOutPlot1;
    ES2eqdskBoozer(ttime);
#endif
    break;
  case 17:			/* Real Space Equilibrium output */
    switch(Flp){
    case 0:
      ESBasicFunctions();
      FlB	=1;
      if(Flw == 1){
	ESWriteESI0();
      }
      CbUserPlot=ESOutPlot0;
      break;
    case 1:
      ESBasicFunctions();
      FlB	=1;
      if(Flw == 1){
	ESWriteESI("ESC",ESsr,ESsz,ESaB,ESgH,ESNa,ESNp);
      }
      CbUserPlot=ESOutPlot0;
      break;
    case 2:
      ES2PEST(VmNa,VmNp);
      if(Flw == 1){
	ESBasicFunctionsPEST(VmNa,VmNp);
	FlB	=2;
      }
      CbUserPlot=ESOutPlot1;
      break;
    case 3:
      break;
    }
    break;
  case 18:			/* Canonical coordinates */
    ES2Canonical(VmNa,VmNp,0);
    break;
  case 19:			/* Output of canonical coordinates */
    ES2Canonical(VmNa,VmNp,1);
    break;
#ifdef VMEC
  case 20:
    printf("VMEC2000 is not yet imcorporated. Sorry.\n");
    break;
  case 21:
    printf("VMEC2000 is not yet imcorporated. Sorry.\n");
    break;
  case 30:
    if(FlMap != 1){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      FlMap	=1;
    }
    sprintf(FNm,"ES4Pest%4.2f",ttime);
    ES4Pst(nA,ttime);
    printf("Using PEST2-400 version from /scratchu/zakh/Pst\n");
    Vmt0=clock();
    if(ESEditPstInFile(ntor) == 0){
      system("cd Pst;/scratchu/zakh/Pst/go24zb;cd ..");
#ifdef H
      ESDconVsPst(1);
#endif
    }
    break;
  case 31:
    if(FlMap != 1){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      FlMap	=1;
    }
    sprintf(FNm,"ES4Pest%4.2f",ttime);
    ES4Pst(1200,ttime);
    printf("Using PEST2-1200 version from /scratchu/zakh/Pst\n");
    Vmt0=clock();
    if(ESEditPstInFile(ntor) == 0){
      system("cd Pst;/scratchu/zakh/Pst/pst12;cd ..");
#ifdef H
      ESDconVsPst(2);
#endif
    }
    break;
  case 32:
    if(FlMap != 2){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      FlMap	=2;
    }
    sprintf(FNm,"ES4Pest%4.2f",ttime);
    ES4Pst0(nA,ttime);
    printf("Using PEST2-400 version from /scratchu/zakh/Pst\n");
    system("cd Pst;/scratchu/zakh/Pst/go24zb;cd ..");
#ifdef H
    ESDconVsPst(2);
#endif
    break;
  case 33:
    if(FlMap != 1){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      FlMap	=1;
    }
    sprintf(FNm,"ES4Pest%4.2f",ttime);
    ES4Pst0(1200,ttime);
    printf("Using PEST2-1200 version from /scratchu/zakh/Pst\n");
    Vmt0=clock();
    if(ESEditPstInFile(ntor) == 0){
      system("cd Pst;/scratchu/zakh/Pst/pst12;cd ..");
#ifdef H
      ESDconVsPst(2);
#endif
    }
    break;
  case 34:
    if(FlMap != 3){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      FlMap	=3;
    }
    sprintf(FNm,"ES4Pest%4.2f",ttime);
    ES4Pst1(nA,ttime);
    printf("Using PEST2-400 version from /scratchu/zakh/Pst\n");
    system("cd Pst;/scratchu/zakh/Pst/go24zb;cd ..");
    break;
  case 35:
    if(FlMap != 3){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      FlMap	=3;
    }
    sprintf(FNm,"ES4Pest%4.2f",ttime);
    ES4Pst1(1200,ttime);
    printf("Using PEST2-400 version from /scratchu/zakh/Pst\n");
    system("cd Pst;/scratchu/zakh/Pst/pst12;cd ..");
    break;
  case 36:
    if(FlMap != 4){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      FlMap	=4;
    }
    sprintf(FNm,"ES4Pest%4.2f",ttime);
    ES4Pst2(nA,ttime);
    printf("Using PEST2-400 version from /scratchu/zakh/Pst\n");
    system("cd Pst;/scratchu/zakh/Pst/go24zb;cd ..");
    break;
  case 37:
    if(FlMap != 4){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      FlMap	=3;
    }
    sprintf(FNm,"ES4Pest%4.2f",ttime);
    ES4Pst2(1200,ttime);
    printf("Using PEST2-1200 version from /scratchu/zakh/Pst\n");
    system("cd Pst;/scratchu/zakh/Pst/pst12;cd ..");
    break;
#endif
  case 40:
    ECSaveESC();
    break;
  case 41:
    ES2HardDrive(ttime);
    break;
#ifdef HH
  case 5:
    ES2eqdskS(ttime);
    ttime	+=0.01;
    break;
  case 6:
    ES2eqdsk(ttime);
    ttime	+=0.01;
    break;
  case 7:
    if((ESEqSolvFl&0x0F) == 0){
      ES4Dcon(ttime);
      ES4Dcon1();
      ES4DconF();
    }
    else{
      printf("Switch Eq Solver to Lsode%c\n",7);
    }
    break;
  case 8:
    if((ESEqSolvFl&0x0F) == 0){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Dcon%4.2f",ttime);
      ES4Dcon(ttime);
      FlMap	=5;
      if(ESEditDconInFile(FNm,ntor) == 0){
	Vmt0=clock();
	system("/scratchu/zakh/dcon_3.10/dcon/dcon;");
#ifdef H
	system("cd /scratchu/zakh/dcon_3.10/dcon;"
	       "dcon;cd /scratchu/zakh/Tsc");
	ESDconVsPst(0);
#endif
      }
    }
    else{
      printf("Switch Eq Solver to Lsode%c\n",7);
    }
    break;
  case 9:
    break;

    if((ESEqSolvFl&0x0F) == 0){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Dcon%4.2f",ttime);
      ES4Dcon(ttime);
      FlMap	=5;
      if(ESEditDconInFile(FNm,ntor) == 0){
	Vmt0=clock();
	system("cd /scratchu/zakh/dcon_3.10/dcon;"
	       "dconNew;cd /scratchu/zakh/Tsc");
	ESDconVsPst(0);
      }
    }
    else{
      printf("Switch Eq Solver to Lsode%c\n",7);
    }
    break;
  case 10:
    if(FlMap != 1){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst(nA,ttime);
      FlMap	=1;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;gobal; cd /scratchu/zakh/Tsc");
    break;
  case 11:
    if(FlMap != 1){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst(nA,ttime);
      FlMap	=1;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    Vmt0=clock();
    if(ESEditPstInFile(ntor) == 0){
      system("cd /scratchu/zakh/Pst;go24zb; cd /scratchu/zakh/Tsc");
      ESDconVsPst(1);
    }
    break;
  case 12:
    if(FlMap != 1){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst(nA,ttime);
      FlMap	=1;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;gobal");
    if(ESEditPstInFile(ntor) == 0){
      system("go24zb");
    }
    system("cd /scratchu/zakh/Tsc");
    break;
  case 13:
    if(FlMap != 1){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst(1200,ttime);
      FlMap	=1;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    Vmt0=clock();
    if(ESEditPstInFile(ntor) == 0){
      system("cd /scratchu/zakh/Pst;pst12; cd /scratchu/zakh/Tsc");
      ESDconVsPst(2);
    }
    break;
  case 20:
    if(FlMap != 2){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst0(nA,ttime);
      FlMap	=2;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;gobal; cd /scratchu/zakh/Tsc");
    break;
  case 21:
    if(FlMap != 2){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst0(nA,ttime);
      FlMap	=2;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;go24zb; cd /scratchu/zakh/Tsc");
    break;
  case 22:
    if(FlMap != 2){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst0(nA,ttime);
      FlMap	=2;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;gobal;go24zb; cd /scratchu/zakh/Tsc");
    break;
  case 30:
    if(FlMap != 3){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst1(nA,ttime);
      FlMap	=3;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;gobal; cd /scratchu/zakh/Tsc");
    break;
  case 31:
    if(FlMap != 3){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst1(nA,ttime);
      FlMap	=3;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;go24zb; cd /scratchu/zakh/Tsc");
    break;
  case 32:
    if(FlMap != 3){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst1(nA,ttime);
      FlMap	=3;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;gobal;go24zb; cd /scratchu/zakh/Tsc");
    break;
  case 40:
    if(FlMap != 4){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst2(nA,ttime);
      FlMap	=4;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;gobal; cd /scratchu/zakh/Tsc");
    break;
  case 41:
    if(FlMap != 4){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst2(nA,ttime);
      FlMap	=4;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;go24zb; cd /scratchu/zakh/Tsc");
    break;
  case 42:
    if(FlMap != 4){
      if(FlMap == 0){
	ttime	+=0.01;
      }
      sprintf(FNm,"ES4Pest%4.2f",ttime);
      ES4Pst2(nA,ttime);
      FlMap	=4;
    }
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;gobal;go24zb; cd /scratchu/zakh/Tsc");
    break;
  case 50:
    ES2eqdskS(ttime);
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;god24zb263;gobal; cd /scratchu/zakh/Tsc");
    break;
  case 51:
    ES2eqdskS(ttime);
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;god24zb263;go24zb; cd /scratchu/zakh/Tsc");
    break;
  case 52:
    ES2eqdskS(ttime);
    printf("Using PEST2 version from /scratchu/zakh/Pst\n");
    system("cd /scratchu/zakh/Pst;god24zb263;gobal;go24zb;"
	   "cd /scratchu/zakh/Tsc");
    break;
#endif
  default:
    break;
  }
  if(FlB == 0 & FlIap){
    sprintf(Message,"First, use Fl=17\nto generate Basic Functions");
    CbUserWarning	=Message;
    FlIap	=0;
  }
  switch(FlIap){
  case 0:
    ESCheckGShFourier(i0,i1,m0,m1,Flw ? Flw-1 : 0);
    break;
  case 1:
    switch(Flp){
    case 0:
      na1	=ESNa1;
      np1	=ESNp1;
      esiread_("esiA.00");
      break;
    case 1:
      na1	=ESNa1;
      np1	=ESNp1;
      esiread_("esiE.00");
      break;
    case 2:
      na1	=VmNa1;
      np1	=VmNp1;
      esiread_("esiP.00");
      break;
    case 3:
      na1	=VmNa1;
      np1	=VmNp1;
      esiread_("esiB.00");
      break;
    }
    ESICopy(R4,Z4,B4,H4,na1,np1);
    esifree_();
    if(FlIa){
      ESCheckBasicFunctionsA(i0,i1,j0,j1,na1,np1
			     ,(FlIa-1)*10+FlIi,Flw ? Flw-1 : 0);
    }
    if(FlIp){
      ESCheckBasicFunctionsP(i0,i1,j0,j1,na1,np1
			     ,(FlIp-1)*10+FlIj,Flw ? Flw-1 : 0);
    }
    break;
  }
  Flw	=0;
  Fl	=0;
#endif /*Tbl_EqOut*/
  if(iAMap != NULL){
    free(iAMap);
    nMap	=0;
    iAMap	=NULL;
  }
  return(0);
}
#endif

#define NSIZE 263
#define NA 128
#define NP 128

#define NAP 1200
#define NPP 128
typedef struct{
  char		Title[160];
  char		date[8];
  long int	nxx[5];
  double	axx[13];
  long int	nxy[10];
  double	axy[10];
}PPPLEQ;

static PPPLEQ Ttl;

int ES4Pst(int nA,double ttime)
{
  int nA1=NAP+1,nP=NPP,nP1=NPP+1,nP5=NPP+5;
  FILE *Fp,*fp;
  int i,j,k,ki,kj,ji;
  double dgy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double A[NAP+1];
  double *sp,*Ps,sq[NAP+1],dsqgy[NAP+1],*g,*dggy,*fb,*dfbgy;
  double *r,*z,*rt,*zt,*ra,*za;
  double *g22,*p2r,*g11,*g12,*dg12t,*dp2ra,*dg22t,*dp2rt,*qg,*dqga,*gd,*dgda;

  double *rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,EZdza,EZdzt,cs,sn,b,db,d2b;
  double *Gh,*Gh2,gh,t;
  double *d2rttc,*d2rtts,*d2ratc,*d2rats,*d2raac,*d2raas;
  double d2rtt,d2ztt,d2rat,d2zat,d2raa,d2zaa;
  double *Lc,*Ls,*dLc,*dLs,D0,dD0,D,dD,dg22,dg12,dgqt,gqt,L0;

  double *w1[6],*w2[6],*w3[6];


  double dGh[65],dGh2[65],chhp,ch2hp;


  nA1=nA+1;

  chhp	=EZcr2*ESgt[1];
  ch2hp	=chhp*EZcr6*ESgt[1];

  i	=nP*nA1;
  for(k=0; k < 6; k++){
    w1[k]	=NULL;
    w2[k]	=NULL;
    w3[k]	=NULL;
    w1[k]	=(double*)malloc(i*sizeof(double));
    w2[k]	=(double*)malloc(i*sizeof(double));
    w3[k]	=(double*)malloc(i*sizeof(double));
    if(w1[k] == NULL || w2[k] == NULL || w3[k] == NULL){
      printf("No memory for w[%d] in ES4Pst()%c\n",k,7);
      j=k;
      while(j >= 0){
	free(w3[j]);
	free(w2[j]);
	free(w1[j]);
	j--;
      }
      return(0);
    }
  }

  r	=w1[0];
  z	=w1[1];
  rt	=w1[2];
  zt	=w1[3];
  ra	=w1[4];
  za	=w1[5];

  g22	=w2[0];
  p2r	=w2[1];
  g11	=w2[2];
  g12	=w2[3];
  dg12t	=w2[4];
  dp2ra	=w2[5];
  dg22t	=w3[0];
  dp2rt	=w3[1];
  qg	=w3[2];
  dqga	=w3[3];
  gd	=w3[4];
  dgda	=w3[5];
  
  Gh	=(double*)malloc((4*ESNp1+12*ESFp1+4*ES2Mp1)*sizeof(double));
  Gh2	=Gh	+ESNp1;
  rc	=Gh2	+ESNp1;
  rs	=rc	+ESFp1;
  drc	=rs	+ESFp1;
  drs	=drc	+ESFp1;
  rct	=drs	+ESFp1;
  rst	=rct	+ESFp1;
  d2raac	=rst	+ESFp1;
  d2raas	=d2raac	+ESFp1;
  d2ratc	=d2raas	+ESFp1;
  d2rats	=d2ratc	+ESFp1;
  d2rttc	=d2rats	+ESFp1;
  d2rtts	=d2rttc	+ESFp1;
  Lc	=d2rtts	+ESFp1;
  Ls	=Lc	+ES2Mp1;
  dLc	=Ls	+ES2Mp1;
  dLs	=dLc	+ES2Mp1;

  if(Gh == NULL){
    printf("No memory for Gh in ES4Pst()\n");
    fclose(fp);
    return(0);
  }
  sp	=r;
  Ps	=z;
#ifdef H
  sq	=rt;
  dsqgy	=zt;
#endif
  g	=ra;
  dggy	=za;
  fb	=g22;
  dfbgy	=p2r;

  rg	=1./ESRBt;
  rBt	=ESRext*rg;
  p2rBt	=rBt*rBt;
  rf	=-ESgY[ESNa]*rBt;
  dgt	=EZc2gp/nP;
  dgy	=1./nA;
  
#ifdef H
  sprintf(Ttl.Title,"ES4PstGm%4.2f",ttime); /* fort.31 */
  sprintf(Ttl.Title,"ES4PstGm");
  Fp	=fopen(Ttl.Title,"w");

  sprintf(Ttl.Title,"ES4PstPr%4.2f",ttime); /* fort.32 */
  sprintf(Ttl.Title,"ES4Pst");
  fp	=fopen(Ttl.Title,"w");
#endif
  fp=fopen("Pst","r");
  if(fp == NULL){
    system("mkdir Pst");
  }
  else{
    fclose(fp);
  }
  sprintf(Ttl.Title,"Pst/fort.31"); /* fort.31 */
  Fp	=fopen(Ttl.Title,"w");

  sprintf(Ttl.Title,"Pst/fort.32"); /* fort.32 */
  fp	=fopen(Ttl.Title,"w");

  printf("Writing plasma profiles\n");
  strcpy(Ttl.Title,"ESC has written this garbage");
  i	=strlen(Ttl.Title);
  while(i < 159){
    Ttl.Title[i]	=' ';
    i++;
  }
  Ttl.Title[i]	='\0';
  g[0]	=0.;
  fwrite((void*)g,(size_t)8,(size_t)1,Fp);
  fwrite((void*)g,(size_t)8,(size_t)1,fp);
  i	=0;
  while(i < 8){
    Ttl.date[i]	='\0';
    i++;
  }
  Ttl.nxx[0]	=ESNp;
  Ttl.nxx[1]	=ESNa1;			/* number of radial intervals +1*/ 
  Ttl.nxx[2]	=nA1;			/* Number of magnetic surfaces */
  Ttl.nxx[3]	=nP;			/* Number of poloidal intervals */ 
  Ttl.nxx[4]	=0;			/* Junk */ 
  Ttl.axx[0]	=rdPV[3]-rdPV[2];
  Ttl.axx[1]	=zdPV[0]-zdPV[1];
  Ttl.axx[2]	=zdPV[0];
  Ttl.axx[3]	=ESaR0[0];
  Ttl.axx[4]	=ESRext;
  Ttl.axx[5]	=ESsp[0]*p2rBt;
  Ttl.axx[6]	=ESRBt*rBt/ESRext;
  Ttl.axx[7]	=EZc2gp*ESgY[ESNa]*rBt;
  Ttl.axx[8]	=0.;
  Ttl.axx[9]	=EZc2gp*ESgY[ESNa]*rBt;
  Ttl.axx[10]	=0.;
  Ttl.axx[11]	=0.;
  Ttl.axx[12]	=1.;	/* constant C in rD=Cr^m\na\bgY^n */

  Ttl.nxy[0]	=1;	/* exponent n in rD=Cr^m\na\bgY^n */
  Ttl.nxy[1]	=1;	/* exponent m in rD=Cr^m\na\bgY^n,
			   2 - straight field lines */
  Ttl.nxy[2]	=0;	/* Junk */ 
  Ttl.nxy[3]	=0;	/* Junk */ 
  Ttl.nxy[4]	=0;	/* Junk */
  Ttl.nxy[5]	=0;	/* Junk */
  Ttl.nxy[6]	=0;	/* Junk */ 
  Ttl.nxy[7]	=0;	/* Junk */
  Ttl.nxy[8]	=nP1;	/* Junk */
  Ttl.nxy[9]	=nA1;	/* Junk */

  Ttl.axy[0]	=dgt;
  Ttl.axy[1]	=dgy;
  Ttl.axy[2]	=EZcgp;
  Ttl.axy[3]	=ESRext;
  Ttl.axy[4]	=ESRBt/ESRext;
  for(i=5; i < 10; i++){
    Ttl.axy[i]	=0.;
  }
  fwrite((void*)&Ttl,(size_t)8,(size_t)49,Fp);
  fwrite((void*)&Ttl,(size_t)8,(size_t)49,fp);

  a	=0.;
  i	=0;
  A[i]	=a;
  sp[i]	=ESsp[0]*p2rBt;
  Ps[i]	=-ESPs[0]*rBt;
  sq[i]	=1./ESgm[0];
  if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
    dsqgy[i]	=ESgm2a[0]*sq[i]*sq[i]/(ESdgY[0]*rBt);
  }
  else{
    ESSetSplDCr(a);
    splRDCr2(&s,NULL,dsqgy,7);
    dsqgy[i]	*=sq[i]*sq[i]/(ESdgY[0]*rBt);
  }
  g[i]	=ESaF[0]*rg;
  s	=ESaT[0]*rBt;
  dggy[i]=-s/(ESRext*ESRext*g[i]);
  fb[i]	=-ESaR0[0]*ESdgY[0]*rBt/ESLc[0];
  dfbgy[i]	=ESaR0[0]*(ESdgY2a[0]-ESdgY[0]*ESLc2[0]/ESLc[0])/
    (ESdgY[0]*ESLc[0]);
  for(i=1; i < nA1; i++){
    ss	=sqrt(dgy*i);
    do{
      ESSetSplA(a);
      splRA(&s,&dgya,ESqgY,ESqgY2a);
      a		+=(ss-s)/dgya;
    }while(fabs(ss-s) > 1e-8);
    A[i]	=a;
    ESSetSplA(a);
    if(ESEqSolvInPr != 2){
      splRA(sp+i,NULL,ESsp,ESsp2a);
    }
    else{
      ESSetSplDPr(a);
      splRDPr(sp+i,&ds,2);
    }
    sp[i]	*=p2rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,0);
      Ps[i]	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,1);
      break;
    default:
      splRA(Ps+i,NULL,ESPs,ESPs2a);
      break;
    }
    Ps[i]	*=-rBt;
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    splRA(&s,&ss,ESLc,ESLc2);
    fb[i]	=-ESaR0[0]*dgya/s;
    dgya	*=a;
    dfbgy[i]	=(ESaR0[0]*ds+fb[i]*ss)/(dgya*s);
    fb[i]	*=rBt;
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq[i]	=1./s;
    dsqgy[i]	=ds/(dgya*rBt*s*s);
    splRA(g+i,&ds,ESaF,ESaF2a);
    g[i]	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    dggy[i]	=-s/(ESRext*ESRext*g[i]);
  }
  fwrite((void*)sp,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)Ps,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)sq,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dsqgy,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)g,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dggy,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)fb,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dfbgy,(size_t)8,(size_t)nA1,fp);

  a	=0.;
  ss	=sqrt(dgy*0.01);
  do{
    ESSetSplA(a);
    splRA(&s,&dgya,ESqgY,ESqgY2a);
    a		+=(ss-s)/dgya;
  }while(fabs(ss-s) > 1e-8);
  A[0]	=a;
  
  printf("Writing ESC metric tensor\n");
  
  ji	=0;
  for(i=0; i < nA1; i++){
    a	=A[i];
    ESSetSplA(a);
    splRA(&dgya,&dD0,ESdgY,ESdgY2a);
    rf	=-1./(dgya*a*rBt);
    splRA(Lc,dLc,ESLc,ESLc2);
    s		=-1./(ESaR0[0]*dgya*rBt);
    D0		=Lc[0]*s;
    dD0		=(dLc[0]-Lc[0]*dD0/dgya)*s*rf;
    L0		=2./Lc[0];
    s		=rf*L0;
    dLc[0]	*=rf/Lc[0];
    ki	=0;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      splRA(Lc+k,dLc+k,ESLc+ki,ESLc2+ki);
      Lc[k]	*=L0;
      dLc[k]	=dLc[k]*s-Lc[k]*dLc[0];
      splRA(Ls+k,dLs+k,ESLs+ki,ESLs2+ki);
      Ls[k]	*=L0;
      dLs[k]	=dLs[k]*s-Ls[k]*dLc[0];
    }
    splRA(&b,&db,ESsb,ESsb2a);
    db	*=rf;
    splRA(rs,drs,rsT,rsT2a);
    drs[0]	*=rf;
    splRA(rc,drc,rcT,rcT2a);
    rc[0]	+=ESaR0[0];
    drc[0]	*=rf;
    ki	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,drc+k,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      drc[k]	*=2.*rf;
      rct[k]	=-k*rc[k];
      d2ratc[k]	=-k*drc[k];
      d2rttc[k]	=k*rct[k];
      splRA(rs+k,drs+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      drs[k]	*=2.*rf;
      rst[k]	=k*rs[k];
      d2rats[k]	=k*drs[k];
      d2rtts[k]	=-k*rst[k];
    }
#ifdef H
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
#endif
    for(j=0; j < ESNp; j++){
      ss	=rc[0];
      EZdra	=drc[0];
      EZdrt	=0.;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	cs	=EScs1[kj];
	sn	=ESsn1[kj];
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
      }
      dGh[j]	=(EZdra*b*EScs1[j]-(drs[0]+db*ESsn1[j])*EZdrt)/ss;
    }
    dGh[j]	=dGh[0];
    dGh2[j]	=dGh2[0];
    splP(dGh,dGh2);
    k		=0;
    Gh[0]	=0.;
    for(j=0; j < ESNp; j++){
      k++;
      Gh[k]	=Gh[j]+chhp*(dGh[j]+dGh[k])-ch2hp*(dGh2[j]+dGh2[k]);
    }

    L0	=EZc2gp/Gh[k];
    for(j=0; j < ESNp1; j++){
      dGh[j]	*=L0;
      dGh2[j]	*=L0;
      Gh[j]	=Gh[j]*L0-ESgt[j];
    }
    L0	=1./D0;
    splP(Gh,Gh2);
    for(j=0; j < nP; j++){
      t		=EZc2gp-dgt*j;
      s		=t;
      cs	=cos(s);
      sn	=sin(s);
      ss	=rc[0]+rc[1]*cs+rs[1]*sn;
      EZdra	=drc[0]+drc[1]*cs+drs[1]*sn;
      EZdrt	=rct[1]*sn+rst[1]*cs;
      ds	=dLc[1]*sn+dLs[1]*(1.-cs);
      gh	=Lc[1]*sn+Ls[1]*(1.-cs);
      gqt	=1.+Lc[1]*cs+Ls[1]*sn;
      dgqt	=dLc[1]*cs+dLs[1]*sn;
      z[ji]	=rs[0]+b*sn;
      EZdza	=drs[0]+db*sn;
      EZdzt	=b*cs;
      d2rat	=d2ratc[1]*sn+d2rats[1]*cs;
      d2rtt	=d2rttc[1]*cs+d2rtts[1]*sn;
      d2zat	=db*cs;
      d2ztt	=-b*sn;
      for(k=2; k < ESFp1; k++){
	s	=k*t;
	cs	=cos(s);
	sn	=sin(s);
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
	d2rat	+=d2ratc[k]*sn+d2rats[k]*cs;
	d2rtt	+=d2rttc[k]*cs+d2rtts[k]*sn;
	if(k <  ES2Mp1){
	  ds	+=(dLc[k]*sn+dLs[k]*(1.-cs))/k;
	  gh	+=(Lc[k]*sn+Ls[k]*(1.-cs))/k;
	  gqt	+=Lc[k]*cs+Ls[k]*sn;
	  dgqt	+=dLc[k]*cs+dLs[k]*sn;
	}
      }
      r[ji]	=ss;
      p2r[ji]	=ss*ss;
      ra[ji]	=EZdra;
      za[ji]	=EZdza;
      rt[ji]	=-EZdrt;
      zt[ji]	=-EZdzt;
      dp2rt[ji]	=-2.*ss*EZdrt;
      dp2ra[ji]	=2.*ss*EZdra;

      D		=EZdra*EZdzt-EZdrt*EZdza;
      g22[ji]	=(EZdrt*EZdrt+EZdzt*EZdzt)/(D*D);
      g11[ji]	=(EZdra*EZdra+EZdza*EZdza)/(D*D);
      g12[ji]	=(EZdra*EZdrt+EZdza*EZdzt)/(D*D);
      
      dD	=d2rat*EZdzt+EZdra*d2ztt-d2rtt*EZdza-EZdrt*d2zat;
      dg12	=d2rat*EZdrt+EZdra*d2rtt+d2zat*EZdzt+EZdza*d2ztt;
      dg22	=2.*(d2rtt*EZdrt+d2ztt*EZdzt);

      dg12t[ji]	=(2.*g12[ji]*dD-dg12/D)/D;
      dg22t[ji]	=(2.*g22[ji]*dD-dg22/D)/D;
      
      qg[ji]	=ss*D;
#ifdef H
      dqga[ji]	=EZdra*D+ss*(d2raa*d2zat+EZdra*d2zat-d2rat*EZdza-EZdrt*d2zaa)
	-qg[ji]*d2gya;
      qg[ji]	=p2r[ji]*D0*gqt;
#endif

      dqga[ji]	=(dp2ra[ji]*D0+p2r[ji]*dD0)*gqt+p2r[ji]*D0*dgqt;

#ifdef H
      ESSetSplP(t);
      splRP(&s,NULL,Gh,Gh2);
      gd[ji]	=-s;
      dgda[ji]	=-ds;
#endif
      gd[ji]	=-gh;
      dgda[ji]	=-sq[i]*ds-dsqgy[i]*gh;
      ji++;
    }
  }
  free(Gh);
  {
    double *p,*P,d[3]={0.,0.,0.};

    ji	=nP*nA;
    /* $r_b$ into Esc2PstGm*/
    P	=r+ji;
    fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
    fwrite((void*)P,(size_t)8,(size_t)2,Fp);
    /* $z_b$ into Esc2PstGm*/
    P	=z+ji;
    fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
    fwrite((void*)P,(size_t)8,(size_t)2,Fp);
    /* $\bgY_{PEST}$ into Esc2PstGm*/
    s	=ESgY[ESNa]*EZc2gp*rBt;
    ds	=-s*dgy;
    for(i=0; i < nA1; i++){
      A[i]	=s+ds*i;
    }
    fwrite((void*)A,(size_t)8,(size_t)nA1,Fp);
    /* $r$ and $|\na\bgY|^2$*/
    /* $z$ and $r^2$*/
    /* $r'_{\gq_{PEST}}$ and $|\na\gq_{PEST}|^2$*/
    /* $z'_{\gq_{PEST}}$ and $(\na\bgY_{PEST}\cdot\na\gq_{PEST})$ */
    /* $r'_{\bgY_{PEST}}$ and $\R{\pa(\na\bgY\cdot\na\gq)}{\pa\gq}$*/
    /* $z'_{\bgY_{PEST}}$ and $\R{\pa r^2}{\pa a}$ */
    for(k=0; k < 6; k++){
      p	=w1[k];
      P	=w2[k];
      for(i=0; i < nA1; i++){
	fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
	fwrite((void*)P,(size_t)8,(size_t)2,Fp);
	fwrite((void*)d,(size_t)8,(size_t)3,Fp);
	P	+=nP;
	fwrite((void*)p,(size_t)8,(size_t)nP,fp);
	fwrite((void*)p,(size_t)8,(size_t)2,fp);
	fwrite((void*)d,(size_t)8,(size_t)3,fp);
	p	+=nP;
      }
      if(k == 0){
	for(j=0; j < nP5; j++){
	  r[j]	=0.;
	}
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,Fp);
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }

    /* $\R{\pa|\na\bgY|^2}{\pa\gq}$ */
    /* $\R{\pa r^2}{\pa\gq}$ */
    /* $\q{g_{PEST}}$ */
    /* $\R{\pa\q{g_{PEST}}}{\pa a}$ */
    /* $\gd_{PEST}$ */
    /* $\R{\pa\gd_{PEST}}{\pa a}$ */
    for(k=0; k < 6; k++){
      P	=w3[k];
      for(i=0; i < nA1; i++){
	fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
	fwrite((void*)P,(size_t)8,(size_t)2,Fp);
	fwrite((void*)d,(size_t)8,(size_t)3,Fp);
	P	+=nP;
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,Fp);
    }
    for(k=0; k < 6; k++){
      free(w3[k]);
      free(w2[k]);
      free(w1[k]);
    }
  }
  printf("ES4PstGm%4.2f (fort.31) has been written\n",ttime);
  fclose(Fp);
  printf("ES4PstPr%4.2f (fort.32) has been written\n",ttime);
  fclose(fp);
  return(0);
}

/* Mapper with $\gz=\gf$ */
int ES4Pst0(int nA,double ttime)
{
  int nA1=NAP+1,nP=NPP,nP1=NPP+1,nP5=NPP+5;
  FILE *Fp,*fp;
  int i,j,k,ki,kj,ji;
  double dgy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double A[NAP+1];
  double *sp,*Ps,*sq,*dsqgy,*g,*dggy,*fb,*dfbgy;
  double *r,*z,*rt,*zt,*ra,*za;
  double *g22,*p2r,*g11,*g12,*dg12t,*dp2ra,*dg22t,*dp2rt,*qg,*dqga,*gd,*dgda;

  double *rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,EZdza,EZdzt,cs,sn,b,db,dz,fac;
  double *Gh,*Gh2,gh,t;
  double *d2rttc,*d2rtts,*d2ratc,*d2rats,d2rtt,d2ztt,d2rat,d2zat;
  double *Lc,*Ls,*dLc,*dLs,D0,dD0,D,dD,dg22,dg12,dgqt,L0;

  double *w1[6],*w2[6],*w3[6];

  double dGh[65],dGh2[65],chhp,ch2hp;

  nA1=nA+1;
  chhp	=EZcr2*ESgt[1];
  ch2hp	=chhp*EZcr6*ESgt[1];

  i	=nP*nA1;
  for(k=0; k < 6; k++){
    w1[k]	=NULL;
    w2[k]	=NULL;
    w3[k]	=NULL;
    w1[k]	=(double*)malloc(i*sizeof(double));
    w2[k]	=(double*)malloc(i*sizeof(double));
    w3[k]	=(double*)malloc(i*sizeof(double));
    if(w1[k] == NULL || w2[k] == NULL || w3[k] == NULL){
      printf("No memory for w[%d] in ES4Pst0()%c\n",k,7);
      j=k;
      while(j >= 0){
	free(w3[j]);
	free(w2[j]);
	free(w1[j]);
	j--;
      }
      return(0);
    }
  }

  r	=w1[0];
  z	=w1[1];
  rt	=w1[2];
  zt	=w1[3];
  ra	=w1[4];
  za	=w1[5];

  g22	=w2[0];
  p2r	=w2[1];
  g11	=w2[2];
  g12	=w2[3];
  dg12t	=w2[4];
  dp2ra	=w2[5];
  dg22t	=w3[0];
  dp2rt	=w3[1];
  qg	=w3[2];
  dqga	=w3[3];
  gd	=w3[4];
  dgda	=w3[5];
  
  Gh	=(double*)malloc((4*ESNp1+10*ESFp1+4*ES2Mp1)*sizeof(double));
  Gh2	=Gh	+ESNp1;
  rc	=Gh2	+ESNp1;
  rs	=rc	+ESFp1;
  drc	=rs	+ESFp1;
  drs	=drc	+ESFp1;
  rct	=drs	+ESFp1;
  rst	=rct	+ESFp1;
  d2ratc	=rst	+ESFp1;
  d2rats	=d2ratc	+ESFp1;
  d2rttc	=d2rats	+ESFp1;
  d2rtts	=d2rttc	+ESFp1;
  Lc	=d2rtts	+ESFp1;
  Ls	=Lc	+ES2Mp1;
  dLc	=Ls	+ES2Mp1;
  dLs	=dLc	+ES2Mp1;

  if(Gh == NULL){
    printf("No memory for Gh in ES4Pst0()\n");
    fclose(fp);
    return(0);
  }
  sp	=r;
  Ps	=z;
  sq	=rt;
  dsqgy	=zt;
  g	=ra;
  dggy	=za;
  fb	=g22;
  dfbgy	=p2r;

  rg	=1./ESRBt;
  rBt	=ESRext*rg;
  p2rBt	=rBt*rBt;
  rf	=-ESgY[ESNa]*rBt;
  dgt	=EZc2gp/nP;
  dgy	=1./nA;
  
#ifdef H
  sprintf(Ttl.Title,"ES4PstGm%4.2f",ttime); /* fort.31 */
  sprintf(Ttl.Title,"ES4PstGm");
  Fp	=fopen(Ttl.Title,"w");

  sprintf(Ttl.Title,"ES4PstPr%4.2f",ttime); /* fort.32 */
  sprintf(Ttl.Title,"ES4Pst");
  fp	=fopen(Ttl.Title,"w");
#endif
  fp=fopen("Pst","r");
  if(fp == NULL){
    system("mkdir Pst");
  }
  else{
    fclose(fp);
  }
  sprintf(Ttl.Title,"Pst/fort.31"); /* fort.31 */
  Fp	=fopen(Ttl.Title,"w");

  sprintf(Ttl.Title,"Pst/fort.32"); /* fort.32 */
  fp	=fopen(Ttl.Title,"w");

  printf("Writing plasma profiles\n");
  strcpy(Ttl.Title,"ESC has written this garbage");
  i	=strlen(Ttl.Title);
  while(i < 159){
    Ttl.Title[i]	=' ';
    i++;
  }
  Ttl.Title[i]	='\0';
  g[0]	=0.;
  fwrite((void*)g,(size_t)8,(size_t)1,Fp);
  fwrite((void*)g,(size_t)8,(size_t)1,fp);
  i	=0;
  while(i < 8){
    Ttl.date[i]	='\0';
    i++;
  }
  Ttl.nxx[0]	=nP1; /* number of poloidal intervals (+1 in symmetric case)*/
  Ttl.nxx[1]	=nA1;			/* number of radial intervals +1*/ 
  Ttl.nxx[2]	=401;			/* Number of magnetic surfaces */
  Ttl.nxx[3]	=128;			/* Number of poloidal intervals */ 
  Ttl.nxx[4]	=0;			/* Junk */ 
  Ttl.axx[0]	=rdPV[3]-rdPV[2];
  Ttl.axx[1]	=zdPV[0]-zdPV[1];
  Ttl.axx[2]	=zdPV[0];
  Ttl.axx[3]	=ESaR0[0];
  Ttl.axx[4]	=ESRext;
  Ttl.axx[5]	=ESsp[0]*p2rBt;
  Ttl.axx[6]	=ESRBt*rBt/ESRext;
  Ttl.axx[7]	=EZc2gp*ESgY[ESNa]*rBt;
  Ttl.axx[8]	=0.;
  Ttl.axx[9]	=EZc2gp*ESgY[ESNa]*rBt;
  Ttl.axx[10]	=0.;
  Ttl.axx[11]	=0.;
  Ttl.axx[12]	=1.;	/* constant C in rD=Cr^m\na\bgY^n */

  Ttl.nxy[0]	=0;	/* exponent n in rD=Cr^m\na\bgY^n */
  Ttl.nxy[1]	=2;	/* exponent m in rD=Cr^m\na\bgY^n,
			   2 - straight field lines */
  Ttl.nxy[2]	=0;	/* Junk */ 
  Ttl.nxy[3]	=0;	/* Junk */ 
  Ttl.nxy[4]	=0;	/* Junk */
  Ttl.nxy[5]	=0;	/* Junk */
  Ttl.nxy[6]	=0;	/* Junk */ 
  Ttl.nxy[7]	=0;	/* Junk */
  Ttl.nxy[8]	=NSIZE;	/* Junk */
  Ttl.nxy[9]	=NSIZE;	/* Junk */

  Ttl.axy[0]	=dgt;
  Ttl.axy[1]	=dgy;
  Ttl.axy[2]	=EZcgp;
  Ttl.axy[3]	=ESRext;
  Ttl.axy[4]	=ESRBt/ESRext;
  for(i=5; i < 10; i++){
    Ttl.axy[i]	=0.;
  }
  fwrite((void*)&Ttl,(size_t)8,(size_t)49,Fp);
  fwrite((void*)&Ttl,(size_t)8,(size_t)49,fp);

  a	=0.;
  i	=0;
  A[i]	=a;
  sp[i]	=ESsp[0]*p2rBt;
  Ps[i]	=-ESPs[0]*rBt;
  sq[i]	=1./ESgm[0];
  if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
    dsqgy[i]	=ESgm2a[0]*sq[i]*sq[i]/(ESdgY[0]*rBt);
  }
  else{
    ESSetSplDCr(a);
    splRDCr2(&s,NULL,dsqgy,7);
    dsqgy[i]	*=sq[i]*sq[i]/(ESdgY[0]*rBt);
  }

  g[i]	=ESaF[0]*rg;
  s	=ESaT[0]*rBt;
  dggy[i]=-s/(ESRext*ESRext*g[i]);
  fb[i]	=-ESaR0[0]*ESdgY[0]*rBt/ESLc[0];
  dfbgy[i]	=ESaR0[0]*(ESdgY2a[0]-ESdgY[0]*ESLc2[0]/ESLc[0])/
    (ESdgY[0]*ESLc[0]);
  for(i=1; i < nA1; i++){
    ss	=sqrt(dgy*i);
    do{
      ESSetSplA(a);
      splRA(&s,&dgya,ESqgY,ESqgY2a);
      a		+=(ss-s)/dgya;
    }while(fabs(ss-s) > 1e-8);
    A[i]	=a;
    ESSetSplA(a);
    if(ESEqSolvInPr != 2){
      splRA(sp+i,NULL,ESsp,ESsp2a);
    }
    else{
      ESSetSplDPr(a);
      splRDPr(sp+i,&ds,2);
    }
    sp[i]	*=p2rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,0);
      Ps[i]	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,1);
      break;
    default:
      splRA(Ps+i,NULL,ESPs,ESPs2a);
      break;
    }
    Ps[i]	*=-rBt;
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    splRA(&s,&ss,ESLc,ESLc2);
    fb[i]	=-ESaR0[0]*dgya/s;
    dgya	*=a;
    dfbgy[i]	=(ESaR0[0]*ds+fb[i]*ss)/(dgya*s);
    fb[i]	*=rBt;
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq[i]	=1./s;
    dsqgy[i]	=ds/(dgya*rBt*s*s);
    splRA(g+i,&ds,ESaF,ESaF2a);
    g[i]	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    dggy[i]	=-s/(ESRext*ESRext*g[i]);
  }
  fwrite((void*)sp,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)Ps,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)sq,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dsqgy,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)g,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dggy,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)fb,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dfbgy,(size_t)8,(size_t)nA1,fp);

  a	=0.;
  ss	=sqrt(dgy*0.01);
  do{
    ESSetSplA(a);
    splRA(&s,&dgya,ESqgY,ESqgY2a);
    a		+=(ss-s)/dgya;
  }while(fabs(ss-s) > 1e-8);
  A[0]	=a;
  
  printf("Writing PEST metric tensor\n");
  
  ji	=0;
  for(i=0; i < nA1; i++){
    a	=A[i];
    ESSetSplA(a);
    splRA(&dgya,&dD0,ESdgY,ESdgY2a);
    rf	=-1./(dgya*a*rBt);
    splRA(Lc,dLc,ESLc,ESLc2);
    s	=1./(ESaR0[0]*dgya*rBt);
    D0	=-Lc[0]*s;
    dD0	=-(dLc[0]*s+D0*dD0/dgya)*rf;
    L0		=2./Lc[0];
    s		=rf*L0;
    dLc[0]	*=rf/Lc[0];
    ki	=0;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      splRA(Lc+k,dLc+k,ESLc+ki,ESLc2+ki);
      Lc[k]	*=L0;
      dLc[k]	=dLc[k]*s-Lc[k]*dLc[0];
      splRA(Ls+k,dLs+k,ESLs+ki,ESLs2+ki);
      Ls[k]	*=L0;
      dLs[k]	=dLs[k]*s-Ls[k]*dLc[0];
    }
    splRA(&b,&db,ESsb,ESsb2a);
    db	*=rf;
    splRA(rs,drs,rsT,rsT2a);
    drs[0]	*=rf;
    splRA(rc,drc,rcT,rcT2a);
    rc[0]	+=ESaR0[0];
    drc[0]	*=rf;
    ki	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,drc+k,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      drc[k]	*=2.*rf;
      rct[k]	=-k*rc[k];
      d2ratc[k]	=-k*drc[k];
      d2rttc[k]	=k*rct[k];
      splRA(rs+k,drs+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      drs[k]	*=2.*rf;
      rst[k]	=k*rs[k];
      d2rats[k]	=k*drs[k];
      d2rtts[k]	=-k*rst[k];
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
#ifdef H
    for(j=0; j < ESNp; j++){
      ss	=rc[0];
      EZdra	=drc[0];
      EZdrt	=0.;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	cs	=EScs1[kj];
	sn	=ESsn1[kj];
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
      }
      dGh[j]	=(EZdra*b*EScs1[j]-(drs[0]+db*ESsn1[j])*EZdrt)/ss;
    }
    dGh[j]	=dGh[0];
    dGh2[j]	=dGh2[0];
    splP(dGh,dGh2);
    k		=0;
    Gh[0]	=0.;
    for(j=0; j < ESNp; j++){
      k++;
      Gh[k]	=Gh[j]+chhp*(dGh[j]+dGh[k])-ch2hp*(dGh2[j]+dGh2[k]);
    }

    L0	=EZc2gp/Gh[k];
    for(j=0; j < ESNp1; j++){
      dGh[j]	*=L0;
      dGh2[j]	*=L0;
      Gh[j]	=Gh[j]*L0-ESgt[j];
    }
#endif
    L0	=1./D0;
    splP(Gh,Gh2);
    t	=EZc2gp;
    for(j=0; j < nP; j++){
      gh	=EZc2gp-dgt*j;
      do{
	ESSetSplP(t);
	splRP(&s,&ds,Gh,Gh2);
	s	+=t-gh;
	t	-=s/(1.+ds);
      }while(fabs(s) > 1e-8);
      s		=t;
      cs	=cos(s);
      sn	=sin(s);
      ss	=rc[0]+rc[1]*cs+rs[1]*sn;
      EZdra	=drc[0]+drc[1]*cs+drs[1]*sn;
      EZdrt	=rct[1]*sn+rst[1]*cs;
      ds	=dLc[1]*sn+dLs[1]*(1.-cs);
      dgqt	=dLc[1]*cs+dLs[1]*sn;
      z[ji]	=rs[0]+b*sn;
      EZdza	=drs[0]+db*sn;
      EZdzt	=b*cs;
      d2rat	=d2ratc[1]*sn+d2rats[1]*cs;
      d2rtt	=d2rttc[1]*cs+d2rtts[1]*sn;
      d2zat	=db*cs;
      d2ztt	=-b*sn;
      for(k=2; k < ESFp1; k++){
	s	=k*t;
	cs	=cos(s);
	sn	=sin(s);
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
	d2rat	+=d2ratc[k]*sn+d2rats[k]*cs;
	d2rtt	+=d2rttc[k]*cs+d2rtts[k]*sn;
	if(k <  ES2Mp1){
	  ds	+=(dLc[k]*sn+dLs[k]*(1.-cs))/k;
	  dgqt+=dLc[k]*cs+dLs[k]*sn;
	}
      }
      r[ji]	=ss;
      D		=EZdra*EZdzt-EZdrt*EZdza;
      g22[ji]	=(EZdrt*EZdrt+EZdzt*EZdzt)/(D*D);
      dD	=d2rat*EZdzt+EZdra*d2ztt-d2rtt*EZdza-EZdrt*d2zat;
      dg12	=d2rat*EZdrt+EZdra*d2rtt+d2zat*EZdzt+EZdza*d2ztt;
      dg22	=2.*(d2rtt*EZdrt+d2ztt*EZdzt);
      s		=ss/(D*L0);
      dg12t[ji]	=(dg22*ds/(D*D)+g22[ji]*(dgqt-2.*ds*dD/D))*s
	-(dg12-(EZdra*EZdrt+EZdza*EZdzt)*(dD*ss+D*EZdrt)/(D*ss))/(D*D);

      dg22t[ji]	=(2.*g22[ji]*dD-dg22/D)*s/D;


      dg12t[ji]	=-(dg12-(EZdra*EZdrt+EZdza*EZdzt)*(dD*ss+D*EZdrt)/(D*ss))/(D*D);


      dg22t[ji]	=(2.*g22[ji]*dD-dg22/D)*s/D;


      rt[ji]	=-EZdrt*s;
      zt[ji]	=-EZdzt*s;

      p2r[ji]	=r[ji]*r[ji];
      dp2rt[ji]	=2.*r[ji]*rt[ji];
      ds	*=s;
      EZdra	-=EZdrt*ds;
      EZdza	-=EZdzt*ds;
      ra[ji]	=EZdra;
      za[ji]	=EZdza;
      g11[ji]	=(EZdra*EZdra+EZdza*EZdza)*L0*L0/(ss*ss);
      g12[ji]	=(EZdra*EZdrt+EZdza*EZdzt)*L0/(D*ss);
      dp2ra[ji]	=2.*ss*EZdra;
      qg[ji]	=ss*ss*D0;
      dqga[ji]	=dp2ra[ji]*D0+p2r[ji]*dD0;
      gd[ji]	=0.;
      dgda[ji]	=0.;
      ji++;
    }
  }
  free(Gh);
  {
    double *p,*P,d[3]={0.,0.,0.};

    ji	=nP*nA;
    /* $r_b$ into Esc2PstGm*/
    P	=r+ji;
    fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
    fwrite((void*)P,(size_t)8,(size_t)2,Fp);
    /* $z_b$ into Esc2PstGm*/
    P	=z+ji;
    fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
    fwrite((void*)P,(size_t)8,(size_t)2,Fp);
    /* $\bgY_{PEST}$ into Esc2PstGm*/
    s	=ESgY[ESNa]*EZc2gp*rBt;
    ds	=-s*dgy;
    for(i=0; i < nA1; i++){
      A[i]	=s+ds*i;
    }
    fwrite((void*)A,(size_t)8,(size_t)nA1,Fp);
    /* $r$ and $|\na\bgY|^2$*/
    /* $z$ and $r^2$*/
    /* $r'_{\gq_{PEST}}$ and $|\na\gq_{PEST}|^2$*/
    /* $z'_{\gq_{PEST}}$ and $(\na\bgY_{PEST}\cdot\na\gq_{PEST})$ */
    /* $r'_{\bgY_{PEST}}$ and $\R{\pa(\na\bgY\cdot\na\gq)}{\pa\gq}$*/
    /* $z'_{\bgY_{PEST}}$ and $\R{\pa r^2}{\pa a}$ */
    for(k=0; k < 6; k++){
      p	=w1[k];
      P	=w2[k];
      for(i=0; i < nA1; i++){
	fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
	fwrite((void*)P,(size_t)8,(size_t)2,Fp);
	fwrite((void*)d,(size_t)8,(size_t)3,Fp);
	P	+=nP;
	fwrite((void*)p,(size_t)8,(size_t)nP,fp);
	fwrite((void*)p,(size_t)8,(size_t)2,fp);
	fwrite((void*)d,(size_t)8,(size_t)3,fp);
	p	+=nP;
      }
      if(k == 0){
	for(j=0; j < nP5; j++){
	  r[j]	=0.;
	}
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,Fp);
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }

    /* $\R{\pa|\na\bgY|^2}{\pa\gq}$ */
    /* $\R{\pa r^2}{\pa\gq}$ */
    /* $\q{g_{PEST}}$ */
    /* $\R{\pa\q{g_{PEST}}}{\pa a}$ */
    /* $\gd_{PEST}$ */
    /* $\R{\pa\gd_{PEST}}{\pa a}$ */
    for(k=0; k < 6; k++){
      P	=w3[k];
      for(i=0; i < nA1; i++){
	fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
	fwrite((void*)P,(size_t)8,(size_t)2,Fp);
	fwrite((void*)d,(size_t)8,(size_t)3,Fp);
	P	+=nP;
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,Fp);
    }
    for(k=0; k < 6; k++){
      free(w3[k]);
      free(w2[k]);
      free(w1[k]);
    }
  }
  printf("ES4PstGm%4.2f (fort.31) has been written\n",ttime);
  fclose(Fp);
  printf("ES4PstPr%4.2f (fort.32) has been written\n",ttime);
  fclose(fp);
  return(0);
}

/* Mapper with $g_{\gq\gq}=const$ */
int ES4Pst1(int nA,double ttime)
{
  int nA1=NAP+1,nP=NPP,nP1=NPP+1,nP5=NPP+5;
  FILE *Fp,*fp;
  int i,j,k,ki,kj,ji;
  double dgy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double A[NAP+1];
  double *sp,*Ps,sq[NAP+1],dsqgy[NAP+1],*g,*dggy,*fb,*dfbgy;
  double *r,*z,*rt,*zt,*ra,*za;
  double *g22,*p2r,*g11,*g12,*dg12t,*dp2ra,*dg22t,*dp2rt,*qg,*dqga,*gd,*dgda;

  double *rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,EZdza,EZdzt,cs,sn,b,db;
  double *gQ,*gQ2,gq,gqa,gqt,gqat,gqtt,t;
  double *d2rttc,*d2rtts,*d2ratc,*d2rats,d2rtt,d2ztt,d2rat,d2zat;
  double *Lc,*Ls,*dLc,*dLs,D0,dD0,D,dD;
  double dg22,dg12,gh,ght,gha,ghat,ghtt,L0,G0;
  double *g22c,*g22c2,*g22s,*g22s2,*G22c,*dG22c,*G22s,*dG22s;

  double *w1[6],*w2[6],*w3[6];

  nA1	=nA+1;
  i	=nP*nA1;
  for(k=0; k < 6; k++){
    w1[k]	=NULL;
    w2[k]	=NULL;
    w3[k]	=NULL;
    w1[k]	=(double*)malloc(i*sizeof(double));
    w2[k]	=(double*)malloc(i*sizeof(double));
    w3[k]	=(double*)malloc(i*sizeof(double));
    if(w1[k] == NULL || w2[k] == NULL || w3[k] == NULL){
      printf("No memory for w[%d] in ES4Pst1()%c\n",k,7);
      j=k;
      while(j >= 0){
	free(w3[j]);
	free(w2[j]);
	free(w1[j]);
	j--;
      }
      return(0);
    }
  }

  r	=w1[0];
  z	=w1[1];
  rt	=w1[2];
  zt	=w1[3];
  ra	=w1[4];
  za	=w1[5];

  g22	=w2[0];
  p2r	=w2[1];
  g11	=w2[2];
  g12	=w2[3];
  dg12t	=w2[4];
  dp2ra	=w2[5];
  dg22t	=w3[0];
  dp2rt	=w3[1];
  qg	=w3[2];
  dqga	=w3[3];
  gd	=w3[4];
  dgda	=w3[5];
  
  gQ	=(double*)malloc((4*ESNp1+18*ESFp1)*sizeof(double));
  gQ2	=gQ	+ESNp1;
  rc	=gQ2	+ESNp1;
  rs	=rc	+ESFp1;
  drc	=rs	+ESFp1;
  drs	=drc	+ESFp1;
  rct	=drs	+ESFp1;
  rst	=rct	+ESFp1;
  d2ratc	=rst	+ESFp1;
  d2rats	=d2ratc	+ESFp1;
  d2rttc	=d2rats	+ESFp1;
  d2rtts	=d2rttc	+ESFp1;
  Lc		=d2rtts	+ESFp1;
  Ls		=Lc	+ESFp1;
  dLc		=Ls	+ESFp1;
  dLs		=dLc	+ESFp1;
  G22c		=dLs	+ESFp1;
  G22s		=G22c	+ESFp1;
  dG22c		=G22s	+ESFp1;
  dG22s		=dG22c	+ESFp1;

  if(gQ == NULL){
    printf("No memory for gQ in ES4Pst1()\n");
    return(0);
  }
  sp	=r;
  Ps	=z;
  g	=ra;
  dggy	=za;
  fb	=g22;
  dfbgy	=p2r;

  i	=ESNa1*ESFp1;
  g22c	=(double*)malloc(4*i*sizeof(double));
  g22c2	=g22c	+i;
  g22s	=g22c2	+i;
  g22s2	=g22s	+i;
  if(g22c == NULL){
    printf("No memory for gQ in ES4Pst1()\n");
    return(0);
  }

  i		=0;
  b		=ESsb1a[i];
  ki		=ESNa1;
  rc[1]		=2.*rcT1a[ki];
  rs[1]		=2.*rsT1a[ki];
  for(j=0; j < ESNp; j++){
    rt[j]	=-rc[1]*ESsn1[j]+rs[1]*EScs1[j];
    zt[j]	=b*EScs1[j];
    g22[j]	=sqrt(rt[j]*rt[j]+zt[j]*zt[j]);
  }
  g22[j]	=g22[0];
  ESgP2gF(g22c+i,g22s+i,g22,ESFp);

  ki		=2*ESNa1;
  rc[2]		=2.*rcT2a[ki];
  rs[2]		=2.*rsT2a[ki];
  for(j=0; j < ESNp; j++){
    EZdra		=-rc[2]*EZsn2[j]+rs[2]*EZcs2[j];
    g22[j]	=EZdra*rt[j]/g22[j];
  }
  g22[j]	=g22[0];
  ESgP2gF(g22c2+i,g22s2+i,g22,ESFp);
  for(i=1; i < ESNa1; i++){
    db	=1./ESsa[i];
    b	=ESsb[i]*db;
    ki	=i;
    db	*=2.;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=k*rcT[ki]*db;
      rs[k]	=k*rsT[ki]*db;
    }
    for(j=0; j < ESNp; j++){
      rt[j]	=0.;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj		+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	rt[j]	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
      }
      zt[j]	=b*EScs1[j];
      g22[j]	=sqrt(rt[j]*rt[j]+zt[j]*zt[j]);
    }
    g22[j]	=g22[0];
    ESgP2gF(g22c+i,g22s+i,g22,ESFp);
  }
  i		=ESNa;
  ki		=i;
  for(k=1; k < ESFp1; k++){
    ki		+=ESNa1;
    b		=2.*k;
    rc[k]	=b*rcT1a[ki];
    rs[k]	=b*rsT1a[ki];
  }
  b		=ESsb1a[i];
  for(j=0; j < ESNp; j++){
    EZdra	=0.;
    kj	=0;
    for(k=1; k < ESFp1; k++){
      kj	+=j;
      if(kj >= ESNp)
	kj	-=ESNp;
      EZdra	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
    }
    EZdza		=b*EScs1[j];
    g22[j]	=(rt[j]*EZdra+zt[j]*EZdza)/g22[j];
  }
  g22[j]	=g22[0];
  ESgP2gF(g22c2+i,g22s2+i,g22,ESFp);
  dg12	=g22c2[0];
  dg22	=g22c2[ESNa];
  splA(g22c,g22c2,&dg12,&dg22);
  ki	=0;
  for(k=1; k < ESFp1; k++){
    ki 		+=ESNa1;
    i		=ki+ESNa;
    dg12	=g22c2[ki];
    dg22	=g22c2[i];
    splA(g22c+ki,g22c2+ki,&dg12,&dg22);
    dg12	=g22s2[ki];
    dg22	=g22s2[i];
    splA(g22s+ki,g22s2+ki,&dg12,&dg22);
  }

  rg	=1./ESRBt;
  rBt	=ESRext*rg;
  p2rBt	=rBt*rBt;
  rf	=-ESgY[ESNa]*rBt;
  dgt	=EZc2gp/nP;
  dgy	=1./nA;
  
#ifdef H
  sprintf(Ttl.Title,"ES4PstGm%4.2f",ttime); /* fort.31 */
  sprintf(Ttl.Title,"ES4PstGm");
  Fp	=fopen(Ttl.Title,"w");

  sprintf(Ttl.Title,"ES4PstPr%4.2f",ttime); /* fort.32 */
  sprintf(Ttl.Title,"ES4Pst");
  fp	=fopen(Ttl.Title,"w");
#endif
  fp=fopen("Pst","r");
  if(fp == NULL){
    system("mkdir Pst");
  }
  else{
    fclose(fp);
  }
  sprintf(Ttl.Title,"Pst/fort.31"); /* fort.31 */
  Fp	=fopen(Ttl.Title,"w");

  sprintf(Ttl.Title,"Pst/fort.32"); /* fort.32 */
  fp	=fopen(Ttl.Title,"w");

  printf("Writing plasma profiles\n");
  strcpy(Ttl.Title,"ESC has written this garbage");
  i	=strlen(Ttl.Title);
  while(i < 159){
    Ttl.Title[i]	=' ';
    i++;
  }
  Ttl.Title[i]	='\0';
  g[0]	=0.;
  fwrite((void*)g,(size_t)8,(size_t)1,Fp);
  fwrite((void*)g,(size_t)8,(size_t)1,fp);
  i	=0;
  while(i < 8){
    Ttl.date[i]	='\0';
    i++;
  }
  Ttl.nxx[0]	=nP1; /* number of poloidal intervals (+1 in symmetric case)*/
  Ttl.nxx[1]	=nA1;			/* number of radial intervals +1*/ 
  Ttl.nxx[2]	=401;			/* Number of magnetic surfaces */
  Ttl.nxx[3]	=128;			/* Number of poloidal intervals */ 
  Ttl.nxx[4]	=0;			/* Junk */ 
  Ttl.axx[0]	=rdPV[3]-rdPV[2];
  Ttl.axx[1]	=zdPV[0]-zdPV[1];
  Ttl.axx[2]	=zdPV[0];
  Ttl.axx[3]	=ESaR0[0];
  Ttl.axx[4]	=ESRext;
  Ttl.axx[5]	=ESsp[0]*p2rBt;
  Ttl.axx[6]	=ESRBt*rBt/ESRext;
  Ttl.axx[7]	=EZc2gp*ESgY[ESNa]*rBt;
  Ttl.axx[8]	=0.;
  Ttl.axx[9]	=EZc2gp*ESgY[ESNa]*rBt;
  Ttl.axx[10]	=0.;
  Ttl.axx[11]	=0.;
  Ttl.axx[12]	=1.;	/* constant C in rD=Cr^m\na\bgY^n */

  Ttl.nxy[0]	=0;	/* exponent n in rD=Cr^m\na\bgY^n */
  Ttl.nxy[1]	=2;	/* exponent m in rD=Cr^m\na\bgY^n,
			   2 - straight field lines */
  Ttl.nxy[2]	=0;	/* Junk */ 
  Ttl.nxy[3]	=0;	/* Junk */ 
  Ttl.nxy[4]	=0;	/* Junk */
  Ttl.nxy[5]	=0;	/* Junk */
  Ttl.nxy[6]	=0;	/* Junk */ 
  Ttl.nxy[7]	=0;	/* Junk */
  Ttl.nxy[8]	=NSIZE;	/* Junk */
  Ttl.nxy[9]	=NSIZE;	/* Junk */

  Ttl.axy[0]	=dgt;
  Ttl.axy[1]	=dgy;
  Ttl.axy[2]	=EZcgp;
  Ttl.axy[3]	=ESRext;
  Ttl.axy[4]	=ESRBt/ESRext;
  for(i=5; i < 10; i++){
    Ttl.axy[i]	=0.;
  }
  fwrite((void*)&Ttl,(size_t)8,(size_t)49,Fp);
  fwrite((void*)&Ttl,(size_t)8,(size_t)49,fp);

  a	=0.;
  i	=0;
  A[i]	=a;
  sp[i]	=ESsp[0]*p2rBt;
  Ps[i]	=-ESPs[0]*rBt;
  sq[i]	=1./ESgm[0];
  if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
    dsqgy[i]	=ESgm2a[0]*sq[i]*sq[i]/(ESdgY[0]*rBt);
  }
  else{
    ESSetSplDCr(a);
    splRDCr2(&s,NULL,dsqgy,7);
    dsqgy[i]	*=sq[i]*sq[i]/(ESdgY[0]*rBt);
  }
  g[i]	=ESaF[0]*rg;
  s	=ESaT[0]*rBt;
  dggy[i]=-s/(ESRext*ESRext*g[i]);
  fb[i]	=-ESaR0[0]*ESdgY[0]*rBt/ESLc[0];
  dfbgy[i]	=ESaR0[0]*(ESdgY2a[0]-ESdgY[0]*ESLc2[0]/ESLc[0])/
    (ESdgY[0]*ESLc[0]);
  for(i=1; i < nA1; i++){
    ss	=sqrt(dgy*i);
    do{
      ESSetSplA(a);
      splRA(&s,&dgya,ESqgY,ESqgY2a);
      a		+=(ss-s)/dgya;
    }while(fabs(ss-s) > 1e-8);
    A[i]	=a;
    ESSetSplA(a);
    if(ESEqSolvInPr != 2){
      splRA(sp+i,NULL,ESsp,ESsp2a);
    }
    else{
      ESSetSplDPr(a);
      splRDPr(sp+i,&ds,2);
    }
    sp[i]	*=p2rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,0);
      Ps[i]	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,1);
      break;
    default:
      splRA(Ps+i,NULL,ESPs,ESPs2a);
      break;
    }
    Ps[i]	*=-rBt;
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    splRA(&s,&ss,ESLc,ESLc2);
    fb[i]	=-ESaR0[0]*dgya/s;
    dgya	*=a;
    dfbgy[i]	=(ESaR0[0]*ds+fb[i]*ss)/(dgya*s);
    fb[i]	*=rBt;
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq[i]	=1./s;
    dsqgy[i]	=ds/(dgya*rBt*s*s);
    splRA(g+i,&ds,ESaF,ESaF2a);
    g[i]	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    dggy[i]	=-s/(ESRext*ESRext*g[i]);
  }
  fwrite((void*)sp,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)Ps,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)sq,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dsqgy,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)g,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dggy,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)fb,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dfbgy,(size_t)8,(size_t)nA1,fp);

  a	=0.;
  ss	=sqrt(dgy*0.01);
  do{
    ESSetSplA(a);
    splRA(&s,&dgya,ESqgY,ESqgY2a);
    a		+=(ss-s)/dgya;
  }while(fabs(ss-s) > 1e-8);
  A[0]	=a;
  
  printf("Writing equidistant metric tensor\n");
  
  ji	=0;
  for(i=0; i < nA1; i++){
    a	=A[i];
    ESSetSplA(a);
    splRA(&dgya,&dD0,ESdgY,ESdgY2a);
    rf	=-1./(dgya*a*rBt);

    splRA(Lc,dLc,ESLc,ESLc2);
    s	=-1./(ESaR0[0]*dgya*rBt);
    D0	=Lc[0]*s;
    dD0	=(dLc[0]-Lc[0]*dD0/dgya)*s*rf;

    L0		=2./Lc[0];
    s		=rf*L0;
    dLc[0]	*=rf/Lc[0];
    ki	=0;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      splRA(Lc+k,dLc+k,ESLc+ki,ESLc2+ki);
      Lc[k]	*=L0;
      dLc[k]	=dLc[k]*s-Lc[k]*dLc[0];
      splRA(Ls+k,dLs+k,ESLs+ki,ESLs2+ki);
      Ls[k]	*=L0;
      dLs[k]	=dLs[k]*s-Ls[k]*dLc[0];
    }

    splRA(G22c,dG22c,g22c,g22c2);
    L0		=2./G22c[0];
    s		=rf*L0;
    dG22c[0]	*=rf/G22c[0];
    ki	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(G22c+k,dG22c+k,g22c+ki,g22c2+ki);
      G22c[k]	*=L0;
      dG22c[k]	=dG22c[k]*s-G22c[k]*dG22c[0];
      splRA(G22s+k,dG22s+k,g22s+ki,g22s2+ki);
      G22s[k]	*=L0;
      dG22s[k]	=dG22s[k]*s-G22s[k]*dG22c[0];
    }

    splRA(&b,&db,ESsb,ESsb2a);
    db	*=rf;
    splRA(rs,drs,rsT,rsT2a);
    drs[0]	*=rf;
    splRA(rc,drc,rcT,rcT2a);
    rc[0]	+=ESaR0[0];
    drc[0]	*=rf;
    ki	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,drc+k,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      drc[k]	*=2.*rf;
      rct[k]	=-k*rc[k];
      d2ratc[k]	=-k*drc[k];
      d2rttc[k]	=k*rct[k];
      splRA(rs+k,drs+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      drs[k]	*=2.*rf;
      rst[k]	=k*rs[k];
      d2rats[k]	=k*drs[k];
      d2rtts[k]	=-k*rst[k];
    }

    for(j=0; j < ESNp; j++){
      s		=0.;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	s	+=(G22c[k]*ESsn1[kj]+G22s[k]*(1.-EScs1[kj]))/k;
      }
      gQ[j]	=s;
    }
    gQ[j]	=gQ[0];
    splP(gQ,gQ2);
    t	=EZc2gp;
    for(j=0; j < nP; j++){
      gq	=EZc2gp-dgt*j;
      do{
	ESSetSplP(t);
	splRP(&s,&ds,gQ,gQ2);
	s	+=t-gq;
	t	-=s/(1.+ds);
      }while(fabs(s) > 1e-8);
      s		=t;
      cs	=cos(s);
      sn	=sin(s);
      ss	=rc[0]+rc[1]*cs+rs[1]*sn;
      EZdra	=drc[0]+drc[1]*cs+drs[1]*sn;
      EZdrt	=rct[1]*sn+rst[1]*cs;

      gqa	=dG22c[1]*sn+dG22s[1]*(1.-cs);
      gqt	=1.+G22c[1]*cs+G22s[1]*sn;
      gqat	=dG22c[1]*cs+dG22s[1]*sn;
      gqtt	=-G22c[1]*sn+G22s[1]*cs;

      gh	=(Lc[1]-G22c[1])*sn+(Ls[1]-G22s[1])*(1.-cs);
      gha	=dLc[1]*sn+dLs[1]*(1.-cs);
      ght	=1.+Lc[1]*cs+Ls[1]*sn;
      ghtt	=-Lc[1]*sn+Ls[1]*cs;
      ghat	=dLc[1]*cs+dLs[1]*sn;

      z[ji]	=rs[0]+b*sn;
      EZdza	=drs[0]+db*sn;
      EZdzt	=b*cs;
      d2rat	=d2ratc[1]*sn+d2rats[1]*cs;
      d2rtt	=d2rttc[1]*cs+d2rtts[1]*sn;
      d2zat	=db*cs;
      d2ztt	=-b*sn;
      for(k=2; k < ESFp1; k++){
	s	=k*t;
	cs	=cos(s);
	sn	=sin(s);
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
	d2rat	+=d2ratc[k]*sn+d2rats[k]*cs;
	d2rtt	+=d2rttc[k]*cs+d2rtts[k]*sn;
	gqa	+=(dG22c[k]*sn+dG22s[k]*(1.-cs))/k;
	gqt	+=G22c[k]*cs+G22s[k]*sn;
	gqat	+=dG22c[k]*cs+dG22s[k]*sn;
	gqtt	+=k*(-G22c[k]*sn+G22s[k]*cs);
	if(k <  ES2Mp1){
	  gh	+=((Lc[k]-G22c[k])*sn+(Ls[k]-G22s[k])*(1.-cs))/k;
	  gha	+=(dLc[k]*sn+dLs[k]*(1.-cs))/k;
	  ght	+=Lc[k]*cs+Ls[k]*sn;
	  ghtt	+=k*(-Lc[k]*sn+Ls[k]*cs);
	  ghat	+=dLc[k]*cs+dLs[k]*sn;
	}
	else{
	  gh	-=(G22c[k]*sn+G22s[k]*(1.-cs))/k;
	}
      }
      gd[ji]	=-gh;
      r[ji]	=ss;
      D		=EZdra*EZdzt-EZdrt*EZdza;
      gh	=1./D;
      dD	=(d2rat*EZdzt+EZdra*d2ztt-d2rtt*EZdza-EZdrt*d2zat)*gh;
#ifdef H
      gh	*=gh;
#endif
      gh	*=1./(ss*D0*ght);

      g22[ji]	=(EZdrt*EZdrt+EZdzt*EZdzt)*gh;
      dg12	=(d2rat*EZdrt+EZdra*d2rtt+d2zat*EZdzt+EZdza*d2ztt)*gh;
      dg22	=2.*(d2rtt*EZdrt+d2ztt*EZdzt)*gh;

      gh	*=(EZdra*EZdrt+EZdza*EZdzt);

      g12[ji]	=gh*gqt-g22[ji]*gqa;

      s		=1./gqt;
      dg12t[ji]	=(2.*g12[ji]*dD-dg12*gqt-gh*gqtt+dg22*gqa+g22[ji]*gqat)*s;
      dg22t[ji]	=(2.*g22[ji]*dD-dg22)*s;

      rt[ji]	=-EZdrt*s;
      zt[ji]	=-EZdzt*s;

      p2r[ji]	=ss*ss;
      dp2rt[ji]	=2.*r[ji]*rt[ji];
      gqa	*=s;
      EZdra	-=EZdrt*gqa;
      EZdza	-=EZdzt*gqa;
      ra[ji]	=EZdra;
      za[ji]	=EZdza;

      dgda[ji]	=-sq[i]*(gha-ght*gqa)+dsqgy[i]*gd[ji];

      D		*=s;
#ifdef H
      g11[ji]	=(EZdra*EZdra+EZdza*EZdza)/(D*D);
#endif

      dp2ra[ji]	=2.*ss*EZdra;

      qg[ji]	=p2r[ji]*D0*ght*s;
      g11[ji]	=(EZdra*EZdra+EZdza*EZdza)*ss/(D*qg[ji]);

      dqga[ji]	=(dp2ra[ji]*D0+p2r[ji]*dD0)*ght*s
	+p2r[ji]*D0*(ghat*gqt-ght*gqat)*s*s
	  -p2r[ji]*D0*(ghtt*gqt-ght*gqtt)*s*s*gqa;
      ji++;
    }
  }
  free(gQ);
  {
    double *p,*P,d[3]={0.,0.,0.};

    ji	=nP*nA;
    /* $r_b$ into Esc2PstGm*/
    P	=r+ji;
    fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
    fwrite((void*)P,(size_t)8,(size_t)2,Fp);
    /* $z_b$ into Esc2PstGm*/
    P	=z+ji;
    fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
    fwrite((void*)P,(size_t)8,(size_t)2,Fp);
    /* $\bgY_{PEST}$ into Esc2PstGm*/
    s	=ESgY[ESNa]*EZc2gp*rBt;
    ds	=-s*dgy;
    for(i=0; i < nA1; i++){
      A[i]	=s+ds*i;
    }
    fwrite((void*)A,(size_t)8,(size_t)nA1,Fp);
    /* $r$ and $|\na\bgY|^2$*/
    /* $z$ and $r^2$*/
    /* $r'_{\gq_{PEST}}$ and $|\na\gq_{PEST}|^2$*/
    /* $z'_{\gq_{PEST}}$ and $(\na\bgY_{PEST}\cdot\na\gq_{PEST})$ */
    /* $r'_{\bgY_{PEST}}$ and $\R{\pa(\na\bgY\cdot\na\gq)}{\pa\gq}$*/
    /* $z'_{\bgY_{PEST}}$ and $\R{\pa r^2}{\pa a}$ */
    for(k=0; k < 6; k++){
      p	=w1[k];
      P	=w2[k];
      for(i=0; i < nA1; i++){
	fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
	fwrite((void*)P,(size_t)8,(size_t)2,Fp);
	fwrite((void*)d,(size_t)8,(size_t)3,Fp);
	P	+=nP;
	fwrite((void*)p,(size_t)8,(size_t)nP,fp);
	fwrite((void*)p,(size_t)8,(size_t)2,fp);
	fwrite((void*)d,(size_t)8,(size_t)3,fp);
	p	+=nP;
      }
      if(k == 0){
	for(j=0; j < nP5; j++){
	  r[j]	=0.;
	}
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,Fp);
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }

    /* $\R{\pa|\na\bgY|^2}{\pa\gq}$ */
    /* $\R{\pa r^2}{\pa\gq}$ */
    /* $\q{g_{PEST}}$ */
    /* $\R{\pa\q{g_{PEST}}}{\pa a}$ */
    /* $\gd_{PEST}$ */
    /* $\R{\pa\gd_{PEST}}{\pa a}$ */
    for(k=0; k < 6; k++){
      P	=w3[k];
      for(i=0; i < nA1; i++){
	fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
	fwrite((void*)P,(size_t)8,(size_t)2,Fp);
	fwrite((void*)d,(size_t)8,(size_t)3,Fp);
	P	+=nP;
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,Fp);
    }
    for(k=0; k < 6; k++){
      free(w3[k]);
      free(w2[k]);
      free(w1[k]);
    }
  }
  printf("ES4PstGm%4.2f (fort.31) has been written\n",ttime);
  fclose(Fp);
  printf("ES4PstPr%4.2f (fort.32) has been written\n",ttime);
  fclose(fp);
  return(0);
}

/* Mapper to Hamada coordinates */
int ES4Pst2(int nA,double ttime)
{
  int nA1=NAP+1,nP=NPP,nP1=NPP+1,nP5=NPP+5;
  FILE *Fp,*fp;
  int i,j,k,ki,kj,ji;
  double dgy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double A[NAP+1];
  double *sp,*Ps,sq[NAP+1],dsqgy[NAP+1],*g,*dggy,*fb,*dfbgy;
  double *r,*z,*rt,*zt,*ra,*za;
  double *g22,*p2r,*g11,*g12,*dg12t,*dp2ra,*dg22t,*dp2rt,*qg,*dqga,*gd,*dgda;

  double *rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,EZdza,EZdzt,cs,sn,b,db;
  double *gQ,*gQ2,gq,gqa,gqt,gqat,gqtt,t;
  double *d2rttc,*d2rtts,*d2ratc,*d2rats,d2rtt,d2ztt,d2rat,d2zat;
  double *Lc,*Ls,*dLc,*dLs,D0,dD0,D,dD;
  double dg22,dg12,gh,ght,gha,ghat,ghtt,L0,G0;
  double *rDc,*rDc2,*rDs,*rDs2,*Jc,*dJc,*Js,*dJs;

  double *w1[6],*w2[6],*w3[6];

  nA1	=nA+1;
  i	=nP*nA1;
  for(k=0; k < 6; k++){
    w1[k]	=NULL;
    w2[k]	=NULL;
    w3[k]	=NULL;
    w1[k]	=(double*)malloc(i*sizeof(double));
    w2[k]	=(double*)malloc(i*sizeof(double));
    w3[k]	=(double*)malloc(i*sizeof(double));
    if(w1[k] == NULL || w2[k] == NULL || w3[k] == NULL){
      printf("No memory for w[%d] in ES4Pst1()%c\n",k,7);
      j=k;
      while(j >= 0){
	free(w3[j]);
	free(w2[j]);
	free(w1[j]);
	j--;
      }
      return(0);
    }
  }

  r	=w1[0];
  z	=w1[1];
  rt	=w1[2];
  zt	=w1[3];
  ra	=w1[4];
  za	=w1[5];

  g22	=w2[0];
  p2r	=w2[1];
  g11	=w2[2];
  g12	=w2[3];
  dg12t	=w2[4];
  dp2ra	=w2[5];
  dg22t	=w3[0];
  dp2rt	=w3[1];
  qg	=w3[2];
  dqga	=w3[3];
  gd	=w3[4];
  dgda	=w3[5];
  
  gQ	=(double*)malloc((4*ESNp1+18*ESFp1)*sizeof(double));
  gQ2	=gQ	+ESNp1;
  rc	=gQ2	+ESNp1;
  rs	=rc	+ESFp1;
  drc	=rs	+ESFp1;
  drs	=drc	+ESFp1;
  rct	=drs	+ESFp1;
  rst	=rct	+ESFp1;
  d2ratc	=rst	+ESFp1;
  d2rats	=d2ratc	+ESFp1;
  d2rttc	=d2rats	+ESFp1;
  d2rtts	=d2rttc	+ESFp1;
  Lc		=d2rtts	+ESFp1;
  Ls		=Lc	+ESFp1;
  dLc		=Ls	+ESFp1;
  dLs		=dLc	+ESFp1;
  Jc		=dLs	+ESFp1;
  Js		=Jc	+ESFp1;
  dJc		=Js	+ESFp1;
  dJs		=dJc	+ESFp1;

  if(gQ == NULL){
    printf("No memory for gQ in ES4Pst1()\n");
    return(0);
  }
  sp	=r;
  Ps	=z;
  g	=ra;
  dggy	=za;
  fb	=g22;
  dfbgy	=p2r;

  i	=ESNa1*ESFp1;
  rDc	=(double*)malloc(4*i*sizeof(double));
  rDc2	=rDc	+i;
  rDs	=rDc2	+i;
  rDs2	=rDs	+i;
  if(rDc == NULL){
    printf("No memory for gQ in ES4Pst1()\n");
    return(0);
  }

  i	=ESNa1*ES2Mp1;
  for(ki=0; ki < i; ki++){
    rDc[ki]	=ESLc[ki]	+ESVc[ki];
    rDc2[ki]	=ESLc2[ki]	+ESVc2[ki];
    rDs[ki]	=ESLs[ki]	+ESVs[ki];
    rDs2[ki]	=ESLs2[ki]	+ESVs2[ki];
  }

  rg	=1./ESRBt;
  rBt	=ESRext*rg;
  p2rBt	=rBt*rBt;
  rf	=-ESgY[ESNa]*rBt;
  dgt	=EZc2gp/nP;
  dgy	=1./nA;
  
#ifdef H
  sprintf(Ttl.Title,"ES4PstGm%4.2f",ttime); /* fort.31 */
  sprintf(Ttl.Title,"ES4PstGm");
  Fp	=fopen(Ttl.Title,"w");

  sprintf(Ttl.Title,"ES4PstPr%4.2f",ttime); /* fort.32 */
  sprintf(Ttl.Title,"ES4Pst");
  fp	=fopen(Ttl.Title,"w");
#endif
  fp=fopen("Pst","r");
  if(fp == NULL){
    system("mkdir Pst");
  }
  else{
    fclose(fp);
  }
  sprintf(Ttl.Title,"Pst/fort.31"); /* fort.31 */
  Fp	=fopen(Ttl.Title,"w");

  sprintf(Ttl.Title,"Pst/fort.32"); /* fort.32 */
  fp	=fopen(Ttl.Title,"w");

  printf("Writing plasma profiles\n");
  strcpy(Ttl.Title,"ESC has written this garbage");
  i	=strlen(Ttl.Title);
  while(i < 159){
    Ttl.Title[i]	=' ';
    i++;
  }
  Ttl.Title[i]	='\0';
  g[0]	=0.;
  fwrite((void*)g,(size_t)8,(size_t)1,Fp);
  fwrite((void*)g,(size_t)8,(size_t)1,fp);
  i	=0;
  while(i < 8){
    Ttl.date[i]	='\0';
    i++;
  }
  Ttl.nxx[0]	=nP1; /* number of poloidal intervals (+1 in symmetric case)*/
  Ttl.nxx[1]	=nA1;			/* number of radial intervals +1*/ 
  Ttl.nxx[2]	=401;			/* Number of magnetic surfaces */
  Ttl.nxx[3]	=128;			/* Number of poloidal intervals */ 
  Ttl.nxx[4]	=0;			/* Junk */ 
  Ttl.axx[0]	=rdPV[3]-rdPV[2];
  Ttl.axx[1]	=zdPV[0]-zdPV[1];
  Ttl.axx[2]	=zdPV[0];
  Ttl.axx[3]	=ESaR0[0];
  Ttl.axx[4]	=ESRext;
  Ttl.axx[5]	=ESsp[0]*p2rBt;
  Ttl.axx[6]	=ESRBt*rBt/ESRext;
  Ttl.axx[7]	=EZc2gp*ESgY[ESNa]*rBt;
  Ttl.axx[8]	=0.;
  Ttl.axx[9]	=EZc2gp*ESgY[ESNa]*rBt;
  Ttl.axx[10]	=0.;
  Ttl.axx[11]	=0.;
  Ttl.axx[12]	=1.;	/* constant C in rD=Cr^m\na\bgY^n */

  Ttl.nxy[0]	=0;	/* exponent n in rD=Cr^m\na\bgY^n */
  Ttl.nxy[1]	=2;	/* exponent m in rD=Cr^m\na\bgY^n,
			   2 - straight field lines */
  Ttl.nxy[2]	=0;	/* Junk */ 
  Ttl.nxy[3]	=0;	/* Junk */ 
  Ttl.nxy[4]	=0;	/* Junk */
  Ttl.nxy[5]	=0;	/* Junk */
  Ttl.nxy[6]	=0;	/* Junk */ 
  Ttl.nxy[7]	=0;	/* Junk */
  Ttl.nxy[8]	=NSIZE;	/* Junk */
  Ttl.nxy[9]	=NSIZE;	/* Junk */

  Ttl.axy[0]	=dgt;
  Ttl.axy[1]	=dgy;
  Ttl.axy[2]	=EZcgp;
  Ttl.axy[3]	=ESRext;
  Ttl.axy[4]	=ESRBt/ESRext;
  for(i=5; i < 10; i++){
    Ttl.axy[i]	=0.;
  }
  fwrite((void*)&Ttl,(size_t)8,(size_t)49,Fp);
  fwrite((void*)&Ttl,(size_t)8,(size_t)49,fp);

  a	=0.;
  i	=0;
  A[i]	=a;
  sp[i]	=ESsp[0]*p2rBt;
  Ps[i]	=-ESPs[0]*rBt;
  sq[i]	=1./ESgm[0];
  if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
    dsqgy[i]	=ESgm2a[0]*sq[i]*sq[i]/(ESdgY[0]*rBt);
  }
  else{
    ESSetSplDCr(a);
    splRDCr2(&s,NULL,dsqgy,7);
    dsqgy[i]	*=sq[i]*sq[i]/(ESdgY[0]*rBt);
  }
  g[i]	=ESaF[0]*rg;
  s	=ESaT[0]*rBt;
  dggy[i]=-s/(ESRext*ESRext*g[i]);
  fb[i]	=-ESaR0[0]*ESdgY[0]*rBt/ESLc[0];
  dfbgy[i]	=ESaR0[0]*(ESdgY2a[0]-ESdgY[0]*ESLc2[0]/ESLc[0])/
    (ESdgY[0]*ESLc[0]);
  for(i=1; i < nA1; i++){
    ss	=sqrt(dgy*i);
    do{
      ESSetSplA(a);
      splRA(&s,&dgya,ESqgY,ESqgY2a);
      a		+=(ss-s)/dgya;
    }while(fabs(ss-s) > 1e-8);
    A[i]	=a;
    ESSetSplA(a);
    if(ESEqSolvInPr != 2){
      splRA(sp+i,NULL,ESsp,ESsp2a);
    }
    else{
      ESSetSplDPr(a);
      splRDPr(sp+i,&ds,2);
    }
    sp[i]	*=p2rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,0);
      Ps[i]	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,1);
      break;
    default:
      splRA(Ps+i,NULL,ESPs,ESPs2a);
      break;
    }
    Ps[i]	*=-rBt;
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    splRA(&s,&ss,ESLc,ESLc2);
    fb[i]	=-ESaR0[0]*dgya/s;
    dgya	*=a;
    dfbgy[i]	=(ESaR0[0]*ds+fb[i]*ss)/(dgya*s);
    fb[i]	*=rBt;
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq[i]	=1./s;
    dsqgy[i]	=ds/(dgya*rBt*s*s);
    splRA(g+i,&ds,ESaF,ESaF2a);
    g[i]	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    dggy[i]	=-s/(ESRext*ESRext*g[i]);
  }
  fwrite((void*)sp,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)Ps,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)sq,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dsqgy,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)g,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dggy,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)fb,(size_t)8,(size_t)nA1,fp);
  fwrite((void*)dfbgy,(size_t)8,(size_t)nA1,fp);

  a	=0.;
  ss	=sqrt(dgy*0.01);
  do{
    ESSetSplA(a);
    splRA(&s,&dgya,ESqgY,ESqgY2a);
    a		+=(ss-s)/dgya;
  }while(fabs(ss-s) > 1e-8);
  A[0]	=a;
  
  printf("Writing Hamada metric tensor\n");
  
  ji	=0;
  for(i=0; i < nA1; i++){
    a	=A[i];
    ESSetSplA(a);
    splRA(&dgya,&dD0,ESdgY,ESdgY2a);
    rf	=-1./(dgya*a*rBt);
    splRA(Lc,dLc,ESLc,ESLc2);
    L0		=2./Lc[0];
    s		=rf*L0;
    dLc[0]	*=rf/Lc[0];
    ki	=0;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      splRA(Lc+k,dLc+k,ESLc+ki,ESLc2+ki);
      Lc[k]	*=L0;
      dLc[k]	=dLc[k]*s-Lc[k]*dLc[0];
      splRA(Ls+k,dLs+k,ESLs+ki,ESLs2+ki);
      Ls[k]	*=L0;
      dLs[k]	=dLs[k]*s-Ls[k]*dLc[0];
    }

    splRA(Jc,dJc,rDc,rDc2);
    s	=-ESaR0[0]/(dgya*rBt);
    D0	=Jc[0]*s;
    dD0	=(dJc[0]-Jc[0]*dD0/dgya)*s*rf;

    L0		=2./Jc[0];
    s		=rf*L0;
    dJc[0]	*=rf/Jc[0];
    ki	=0;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      splRA(Jc+k,dJc+k,rDc+ki,rDc2+ki);
      Jc[k]	*=L0;
      dJc[k]	=dJc[k]*s-Jc[k]*dJc[0];
      splRA(Js+k,dJs+k,rDs+ki,rDs2+ki);
      Js[k]	*=L0;
      dJs[k]	=dJs[k]*s-Js[k]*dJc[0];
    }

    splRA(&b,&db,ESsb,ESsb2a);
    db	*=rf;
    splRA(rs,drs,rsT,rsT2a);
    drs[0]	*=rf;
    splRA(rc,drc,rcT,rcT2a);
    rc[0]	+=ESaR0[0];
    drc[0]	*=rf;
    ki	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA(rc+k,drc+k,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      drc[k]	*=2.*rf;
      rct[k]	=-k*rc[k];
      d2ratc[k]	=-k*drc[k];
      d2rttc[k]	=k*rct[k];
      splRA(rs+k,drs+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      drs[k]	*=2.*rf;
      rst[k]	=k*rs[k];
      d2rats[k]	=k*drs[k];
      d2rtts[k]	=-k*rst[k];
    }

    for(j=0; j < ESNp; j++){
      s		=0.;
      kj	=0;
      for(k=1; k < ES2Mp1; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	s	+=(Jc[k]*ESsn1[kj]+Js[k]*(1.-EScs1[kj]))/k;
      }
      gQ[j]	=s;
    }
    gQ[j]	=gQ[0];
    splP(gQ,gQ2);
    t	=EZc2gp;
    for(j=0; j < nP; j++){
      gq	=EZc2gp-dgt*j;
      do{
	ESSetSplP(t);
	splRP(&s,&ds,gQ,gQ2);
	s	+=t-gq;
	t	-=s/(1.+ds);
      }while(fabs(s) > 1e-8);
      s		=t;
      cs	=cos(s);
      sn	=sin(s);
      ss	=rc[0]+rc[1]*cs+rs[1]*sn;
      EZdra	=drc[0]+drc[1]*cs+drs[1]*sn;
      EZdrt	=rct[1]*sn+rst[1]*cs;

      gqa	=dJc[1]*sn+dJs[1]*(1.-cs);
      gqt	=1.+Jc[1]*cs+Js[1]*sn;
      gqat	=dJc[1]*cs+dJs[1]*sn;
      gqtt	=-Jc[1]*sn+Js[1]*cs;

      gh	=(Lc[1]-Jc[1])*sn+(Ls[1]-Js[1])*(1.-cs);
      gha	=dLc[1]*sn+dLs[1]*(1.-cs);
      ght	=1.+Lc[1]*cs+Ls[1]*sn;
      ghtt	=-Lc[1]*sn+Ls[1]*cs;
      ghat	=dLc[1]*cs+dLs[1]*sn;

      z[ji]	=rs[0]+b*sn;
      EZdza	=drs[0]+db*sn;
      EZdzt	=b*cs;
      d2rat	=d2ratc[1]*sn+d2rats[1]*cs;
      d2rtt	=d2rttc[1]*cs+d2rtts[1]*sn;
      d2zat	=db*cs;
      d2ztt	=-b*sn;
      for(k=2; k < ESFp1; k++){
	s	=k*t;
	cs	=cos(s);
	sn	=sin(s);
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
	d2rat	+=d2ratc[k]*sn+d2rats[k]*cs;
	d2rtt	+=d2rttc[k]*cs+d2rtts[k]*sn;
	if(k <  ES2Mp1){
	  gqa	+=(dJc[k]*sn+dJs[k]*(1.-cs))/k;
	  gqt	+=Jc[k]*cs+Js[k]*sn;
	  gqat	+=dJc[k]*cs+dJs[k]*sn;
	  gqtt	+=k*(-Jc[k]*sn+Js[k]*cs);

	  gh	+=((Lc[k]-Jc[k])*sn+(Ls[k]-Js[k])*(1.-cs))/k;
	  gha	+=(dLc[k]*sn+dLs[k]*(1.-cs))/k;
	  ght	+=Lc[k]*cs+Ls[k]*sn;
	  ghtt	+=k*(-Lc[k]*sn+Ls[k]*cs);
	  ghat	+=dLc[k]*cs+dLs[k]*sn;
	}
      }
      gd[ji]	=-gh;
      r[ji]	=ss;
      D		=EZdra*EZdzt-EZdrt*EZdza;
      gh	=1./D;
      dD	=(d2rat*EZdzt+EZdra*d2ztt-d2rtt*EZdza-EZdrt*d2zat)*gh;
#ifdef H
      gh	*=gh;
#endif
      gh	*=ss/(D0*gqt);

      g22[ji]	=(EZdrt*EZdrt+EZdzt*EZdzt)*gh;
      dg12	=(d2rat*EZdrt+EZdra*d2rtt+d2zat*EZdzt+EZdza*d2ztt)*gh;
      dg22	=2.*(d2rtt*EZdrt+d2ztt*EZdzt)*gh;

      gh	*=(EZdra*EZdrt+EZdza*EZdzt);

      g12[ji]	=gh*gqt-g22[ji]*gqa;

      s		=1./gqt;
      dg12t[ji]	=(2.*g12[ji]*dD-dg12*gqt-gh*gqtt+dg22*gqa+g22[ji]*gqat)*s;
      dg22t[ji]	=(2.*g22[ji]*dD-dg22)*s;

      rt[ji]	=-EZdrt*s;
      zt[ji]	=-EZdzt*s;

      p2r[ji]	=ss*ss;
      dp2rt[ji]	=2.*r[ji]*rt[ji];
      gqa	*=s;
      EZdra	-=EZdrt*gqa;
      EZdza	-=EZdzt*gqa;
      ra[ji]	=EZdra;
      za[ji]	=EZdza;

      dgda[ji]	=-sq[i]*(gha-ght*gqa)+dsqgy[i]*gd[ji];

      D		*=s;
#ifdef H
      g11[ji]	=(EZdra*EZdra+EZdza*EZdza)/(D*D);
#endif
      g11[ji]	=(EZdra*EZdra+EZdza*EZdza)*ss/(D*D0);
      dp2ra[ji]	=2.*ss*EZdra;

      qg[ji]	=D0;
      dqga[ji]	=dD0;
      ji++;
    }
  }
  free(gQ);
  {
    double *p,*P,d[3]={0.,0.,0.};

    ji	=nP*nA;
    /* $r_b$ into Esc2PstGm*/
    P	=r+ji;
    fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
    fwrite((void*)P,(size_t)8,(size_t)2,Fp);
    /* $z_b$ into Esc2PstGm*/
    P	=z+ji;
    fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
    fwrite((void*)P,(size_t)8,(size_t)2,Fp);
    /* $\bgY_{PEST}$ into Esc2PstGm*/
    s	=ESgY[ESNa]*EZc2gp*rBt;
    ds	=-s*dgy;
    for(i=0; i < nA1; i++){
      A[i]	=s+ds*i;
    }
    fwrite((void*)A,(size_t)8,(size_t)nA1,Fp);
    /* $r$ and $|\na\bgY|^2$*/
    /* $z$ and $r^2$*/
    /* $r'_{\gq_{PEST}}$ and $|\na\gq_{PEST}|^2$*/
    /* $z'_{\gq_{PEST}}$ and $(\na\bgY_{PEST}\cdot\na\gq_{PEST})$ */
    /* $r'_{\bgY_{PEST}}$ and $\R{\pa(\na\bgY\cdot\na\gq)}{\pa\gq}$*/
    /* $z'_{\bgY_{PEST}}$ and $\R{\pa r^2}{\pa a}$ */
    for(k=0; k < 6; k++){
      p	=w1[k];
      P	=w2[k];
      for(i=0; i < nA1; i++){
	fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
	fwrite((void*)P,(size_t)8,(size_t)2,Fp);
	fwrite((void*)d,(size_t)8,(size_t)3,Fp);
	P	+=nP;
	fwrite((void*)p,(size_t)8,(size_t)nP,fp);
	fwrite((void*)p,(size_t)8,(size_t)2,fp);
	fwrite((void*)d,(size_t)8,(size_t)3,fp);
	p	+=nP;
      }
      if(k == 0){
	for(j=0; j < nP5; j++){
	  r[j]	=0.;
	}
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,Fp);
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }

    /* $\R{\pa|\na\bgY|^2}{\pa\gq}$ */
    /* $\R{\pa r^2}{\pa\gq}$ */
    /* $\q{g_{PEST}}$ */
    /* $\R{\pa\q{g_{PEST}}}{\pa a}$ */
    /* $\gd_{PEST}$ */
    /* $\R{\pa\gd_{PEST}}{\pa a}$ */
    for(k=0; k < 6; k++){
      P	=w3[k];
      for(i=0; i < nA1; i++){
	fwrite((void*)P,(size_t)8,(size_t)nP,Fp);
	fwrite((void*)P,(size_t)8,(size_t)2,Fp);
	fwrite((void*)d,(size_t)8,(size_t)3,Fp);
	P	+=nP;
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,Fp);
    }
    for(k=0; k < 6; k++){
      free(w3[k]);
      free(w2[k]);
      free(w1[k]);
    }
  }
  printf("ES4PstGm%4.2f (fort.31) has been written\n",ttime);
  fclose(Fp);
  printf("ES4PstPr%4.2f (fort.32) has been written\n",ttime);
  fclose(fp);
  return(0);
}

int ES2eqdsk(double ttime)
{
  int nA=128,nA1=129,nA5=NSIZE,nP=256,nP1=257,nP4=260,nP5=NSIZE;
  FILE *fp;
  int i,j,ia,j0=2;
  long int nxy[10];
  char Title[160];
  double dgy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double A[2*NSIZE];
  double sp[NSIZE],Ps[NSIZE],sq[NSIZE],dsqgy[NSIZE],gy[NSIZE];
  double g[NSIZE],dggy[NSIZE],fb[NSIZE],dfbgy[NSIZE];

  double aT[NSIZE];

  rg	=1./ESRBt;
  rBt	=ESRext*rg;
  rf	=-ESgY[ESNa]*rBt;
  p2rBt	=rBt*rBt;

  sprintf(Title,"Esc2dsk%4.2f",ttime);
  sprintf(Title,"../Pst/fort.17");
  printf("Writing plasma profiles\n");

  fp	=fopen(Title,"w");
  strcpy(Title,"ESC has written this garbage");
  i	=strlen(Title);
  while(i < 159){
    Title[i]	=' ';
    i++;
  }
  Title[i]	='\0';
  g[0]	=0.;
  fwrite((void*)g,(size_t)8,(size_t)1,fp);
  fwrite((void*)Title,(size_t)8,(size_t)20,fp);
  fwrite((void*)g,(size_t)8,(size_t)1,fp);
  nxy[0]  =nP;
  nxy[1]  =nA1; 
  nxy[2]  =0;
  nxy[3]  =0; 
  nxy[4]  =0; 
  fwrite((void*)nxy,(size_t)8,(size_t)5,fp);

  g[0]	=rdPV[3]-rdPV[2];
  g[1]	=zdPV[0]-zdPV[1];
  g[2]	=zdPV[0];
  g[3]	=ESaR0[0];
  g[4]	=ESRext;
  g[5]	=ESsp[0]*p2rBt;
  g[6]	=ESRBt*rBt/ESRext;
  g[7]	=EZc2gp*ESgY[ESNa]*rBt;
  g[8]	=0.;
  g[9]	=g[7];
  g[10]=0.;
  g[11]=0.;
  g[12]=1.;
  fwrite((void*)g,(size_t)8,(size_t)13,fp);

  nxy[0]	=2; 
  nxy[1]	=0; 
  nxy[2]	=0; 
  nxy[3]	=0; 
  nxy[4]	=2;
  nxy[5]	=0;
  nxy[6]	=1; 
  nxy[7]	=0;	/* symmetry index */
  nxy[8]	=nP5;	/* poloidal allocation size */
  nxy[9]	=nA5;	/* radial allocation size */
  fwrite((void*)nxy,(size_t)8,(size_t)10,fp);
  
  dgt	=EZc2gp/nP;
  dgy	=1./nA;
  g[0]	=dgt;
  g[1]	=dgy;
  g[2]	=EZcgp;
  g[3]	=ESRext;
  g[4]	=ESRBt/ESRext;
  printf("dt=%15.7e %15.7e %15.7e %15.7e %15.7e\n",g[0],g[1],g[2],g[3],g[4]);
  for(i=5; i < 10; i++){
    g[i]	=0.;
  }
  fwrite((void*)g,(size_t)8,(size_t)10,fp);

  a	=0.;
  i	=0;
  ia	=0;
  A[ia]	=a;
  ia++;
  gy[i]	=ESgY[ESNa]*rBt;
  Ps[i]	=-ESPs[0]*rBt;

  sq[i]	=1./ESgm[0];
  if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
    dsqgy[i]	=ESgm2a[0]*sq[i]*sq[i]/(ESdgY[0]*rBt);
  }
  else{
    ESSetSplDCr(a);
    splRDCr2(&s,NULL,dsqgy,7);
    dsqgy[i]	*=sq[i]*sq[i]/(ESdgY[0]*rBt);
  }

  g[i]	=ESaF[0]*rg;
  splRA(&s,NULL,ESaT,ESaT2a);
  s	=ESaT[0]*rBt;
  aT[i]	=s;
  dggy[i]=-s/(ESRext*ESRext*g[i]);
  fb[i]	=-ESaR0[0]*ESdgY[0]*rBt/ESLc[0];
  dfbgy[i]	=ESaR0[0]*(ESdgY2a[0]-ESdgY[0]*ESLc2[0]/ESLc[0])/
    (ESdgY[0]*ESLc[0]);

  ss	=sqrt(dgy*(i+0.5));
  do{
    ESSetSplA(a);
    splRA(&s,&dgya,ESqgY,ESqgY2a);
    a		+=(ss-s)/dgya;
  }while(fabs(ss-s) > 1e-8);
  A[ia]	=a;
  ia++;
  ESSetSplA(a);
  splRA(g+i,&ds,ESaF,ESaF2a);
  g[i]	*=rg;
  if(ESEqSolvInPr != 2){
    splRA(sp+i,NULL,ESsp,ESsp2a);
  }
  else{
    ESSetSplDPr(a);
    splRDPr(sp+i,&ds,2);
  }
  sp[i]	*=p2rBt;
  for(i=1; i < nA1; i++){
    ss	=sqrt(dgy*i);
    do{
      ESSetSplA(a);
      splRA(&s,&dgya,ESqgY,ESqgY2a);
      a		+=(ss-s)/dgya;
    }while(fabs(ss-s) > 1e-8);
    A[ia]	=a;
    ia++;
    ESSetSplA(a);
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    splRA(&s,&ss,ESLc,ESLc2);
    fb[i]	=-ESaR0[0]*dgya/s;
    dgya	*=a;
    dfbgy[i]	=(ESaR0[0]*ds+fb[i]*ss)/(dgya*s);
    fb[i]	*=rBt;
    gy[i]	=ESgY[ESNa]*(1.-dgy*i)*rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,0);
      Ps[i]	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,1);
      break;
    default:
      splRA(Ps+i,NULL,ESPs,ESPs2a);
      break;
    }
    Ps[i]	*=-rBt;
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq[i]	=1./s;
    dsqgy[i]	=ds*sq[i]*sq[i]/(dgya*rBt);
    splRA(g+i,&ds,ESaF,ESaF2a);
    g[i]	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    aT[i]	=s;
    dggy[i]	=-s/(ESRext*ESRext*g[i]);
    if(i < nA){
      ss	=sqrt(dgy*(i+0.5));
      do{
	ESSetSplA(a);
	splRA(&s,&dgya,ESqgY,ESqgY2a);
	a		+=(ss-s)/dgya;
      }while(fabs(ss-s) > 1e-8);
      ESSetSplA(a);
      splRA(g+i,&ds,ESaF,ESaF2a);
      g[i]	*=rg;
      if(ESEqSolvInPr != 2){
	splRA(sp+i,NULL,ESsp,ESsp2a);
      }
      else{
	ESSetSplDPr(a);
	splRDPr(sp+i,&ds,2);
      }
      sp[i]	*=p2rBt;
    }
    else{
      sp[i]	=0.;
      g[i]	=(ESaF[ESNa]+ESaF1a[ESNa]*(a-ESsa[ESNa]))*rg;
    }
    A[ia]	=a;
    ia++;
  }
  ds	=dgy*ESgY[ESNa]/ESgY1a[ESNa];
  A[ia]	=ds;
  ia++;
  sp[i]	=0.;
  Ps[i]	=-(ESPs[ESNa]+ESPs1a[ESNa]*ds)*rBt;
  Ps[i]	=0.;
  sq[i]	=1./(ESgm[ESNa]+ESgm1a[ESNa]*ds);
  sq[i]	=0.;
  dsqgy[i]	=0.;
  g[i]	=(ESaF[ESNa]+ESaF1a[ESNa]*ds)*rg;
  g[i]	=0.;
  dggy[i]=0.;
  fb[i]	=0.;
  dfbgy[i]=0.;
  gy[i]	=ESgY[ESNa]*(1.-dgy*i)*rBt;

  i++;
  while(i < nA5){
    sp[i]	=0.;
    Ps[i]	=0.;
    sq[i]	=0.;
    dsqgy[i]	=0.;
    g[i]	=0.;
    dggy[i]	=0.;
    fb[i]	=0.;
    dfbgy[i]	=0.;
    gy[i]	=0.;
    i++;
  }
  fwrite((void*)sp,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)Ps,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)sq,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)dsqgy,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)g,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)dggy,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)fb,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)dfbgy,(size_t)8,(size_t)nA5,fp);
  i	=0;
  g[i]	=rf;
  s	=-0.5*nA*ESgY[ESNa]*rBt;
  dggy[i]=2.*s;
  for(i=1; i < nA1; i++){
    g[i]	=rf;
    dggy[i]	=s/i;
  }
  while(i < nA5){
    g[i]	=0.;
    dggy[i]	=0.;
    i++;
  }
  fwrite((void*)g,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)dggy,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)gy,(size_t)8,(size_t)nA5,fp);

  {
    int k,ki,jj;
    double *r,*rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,cs,sn,b,db,dz,fac;

    r	=g;
    rc	=sp;
    rs	=Ps;
    drc	=fb;
    drs	=dfbgy;
    rct	=sq;
    rst	=dsqgy;

    printf("Writing geometry ->r[  0]");
    /* r */
    for(j=0; j < nP4; j++){
      r[j]	=ESaR0[0];
    }
    for(j=nP4; j < nP5; j++){
      r[j]	=0.;
    }
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    ia	=2;
    for(i=1; i < nA1; i++){
      printf("\b\b\b\b%3d]",i);
      a		=A[ia];
      ia	+=2;
      ESSetSplA(a);
      ki	=0;
      splRA(rc,NULL,rcT,rcT2a);
      rc[0]	+=ESaR0[0];
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	splRA(rc+k,NULL,rcT+ki,rcT2a+ki);
	splRA(rs+k,NULL,rsT+ki,rsT2a+ki);
      }
      for(j=0; j < nP4; j++){
	s	=dgt*(j-j0);
	ss	=0.;
	for(k=1; k < ESFp1; k++){
	  rf	=k*s;
	  ss	+=rc[k]*cos(rf)+rs[k]*sin(rf);
	}
	r[j]	=rc[0]+2.*ss;
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }
    rf	=A[ia];
    rc[0]	=ESaR0[0]+rcT[ESNa]+rcT1a[ESNa]*rf;
    ki	=ESNa;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=rcT[ki]+rcT1a[ki]*rf;
      rs[k]	=rsT[ki]+rsT1a[ki]*rf;
    }
    for(j=0; j < nP4; j++){
      s		=dgt*(j-j0);
      ss	=0.;
      for(k=1; k < ESFp1; k++){
	rf	=k*s;
	ss	+=rc[k]*cos(rf)+rs[k]*sin(rf);
      }
      r[j]	=rc[0]+2.*ss;
    }
    i++;
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    for(j=0; j < nP4; j++){
      r[j]	=0.;
    }
    while(i < nA5){
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
      i++;
    }
    /* z */
    printf("\b\b\b\b\b\bz[%3d]",0);
    for(j=0; j < nP4; j++){
      r[j]	=ESaZ0[0];
    }
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    ia	=2;
    for(i=1; i < nA1; i++){
      printf("\b\b\b\b%3d]",i);
      a		=A[ia];
      ia	+=2;
      ESSetSplA(a);
      splRA(rs,NULL,rsT,rsT2a);
      rs[0]	+=ESaZ0[0];
      splRA(&b,NULL,ESsb,ESsb2a);
      for(j=0; j < nP4; j++){
	r[j]	=rs[0]+b*sin(dgt*(j-j0));
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }
    b	=ESsb[ESNa]+ESsb1a[ESNa]*A[ia];
    rs[0]=ESaZ0[0]+rsT[ESNa]+rsT1a[ESNa]*A[ia];
    for(j=0; j < nP4; j++){
      r[j]	=rs[0]+b*sin(dgt*(j-j0));
    }
    i++;
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    for(j=0; j < nP4; j++){
      r[j]	=0.;
    }
    while(i < nA5){
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
      i++;
    }
    
    fac	=-ESgY[ESNa]*rBt;
    /* qg05 */
    printf("\b\b\b\b\b\bqg05[%3d]",0);
    s	=-ESaR0[0]*ESDPc[0]*fac/(ESdgY[0]*rBt);
    for(j=0; j < nP4; j++){
      r[j]	=s;
    }
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    ia	=2;
    for(i=1; i < nA1; i++){
      printf("\b\b\b\b%3d]",i);
      a		=A[ia];
      ia	+=2;
      ESSetSplA(a);
      splRA(&rf,&dz,rsT,rsT2a);
      splRA(&dgya,NULL,ESdgY,ESdgY2a);
      rf	=-fac/(dgya*a*rBt);
      splRA(&b,&db,ESsb,ESsb2a);
      db	*=rf;
      dz	*=rf;
      splRA(rc,drc,rcT,rcT2a);
      rc[0]	+=ESaR0[0];
      drc[0]	*=rf;
      ki	=0;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	splRA(rc+k,drc+k,rcT+ki,rcT2a+ki);
	drc[k]	*=rf;
	rct[k]	=-k*rc[k];
	splRA(rs+k,drs+k,rsT+ki,rsT2a+ki);
	drs[k]	*=rf;
	rst[k]	=k*rs[k];
      }
      for(j=0; j < nP4; j++){
	s	=dgt*(j-j0-0.5);
	ss	=0.;
	EZdra	=0.;
	EZdrt	=0.;
	for(k=1; k < ESFp1; k++){
	  rf	=k*s;
	  cs	=cos(rf);
	  sn	=sin(rf);
	  ss	+=rc[k]*cs+rs[k]*sn;
	  EZdra	+=drc[k]*cs+drs[k]*sn;
	  EZdrt	+=rct[k]*sn+rst[k]*cs;
	}
	EZdra	=drc[0]+2.*EZdra;
	EZdrt	*=2.;
	r[j]	=(rc[0]+2.*ss)*(EZdra*b*cos(s)-(dz+db*sin(s))*EZdrt);
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }
    s	=A[ia];
    a	=ESsa[ESNa]+s;
    dgya=ESdgY[ESNa]+ESdgY1a[ESNa]*s;
    rf	=-fac/(dgya*a*rBt);
    b	=ESsb[ESNa]+ESsb1a[ESNa]*s;
    db	=(ESsb1a[ESNa]+ESsb2a[ESNa]*s)*rf;
    dz	=(rsT1a[ESNa]+rsT2a[ESNa]*s)*rf;
    rc[0]	=ESaR0[0]+rcT[ESNa]+rcT1a[ESNa]*s;
    drc[0]	=(rcT1a[ESNa]+rcT2a[ESNa]*s)*rf;
    ki	=ESNa;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=rcT[ki]+rcT1a[ki]*s;
      drc[k]	=(rcT1a[ki]+rcT2a[ki]*s)*rf;
      rct[k]	=-k*rc[k];
      rs[k]	=rsT[ki]+rsT1a[ki]*s;
      drs[k]	=(rsT1a[ki]+rsT2a[ki]*s)*rf;
      rst[k]	=k*rs[k];
    }
    for(j=0; j < nP4; j++){
      s		=dgt*(j-j0-0.5);
      ss	=0.;
      EZdra	=0.;
      EZdrt	=0.;
      for(k=1; k < ESFp1; k++){
	rf	=k*s;
	cs	=cos(rf);
	sn	=sin(rf);
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
      }
      EZdra	=drc[0]+2.*EZdra;
      EZdrt	*=2.;
      r[j]	=(rc[0]+2.*ss)*(EZdra*b*cos(s)-(dz+db*sin(s))*EZdrt);
    }
    i++;
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    for(j=0; j < nP4; j++){
      r[j]	=0.;
    }
    while(i < nA5){
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
      i++;
    }
    /* qg */
    printf("\b\b\b\b\b\b\b\b\bqg[%3d]",0);
    s	=-ESaR0[0]*ESDPc[0]*fac/(ESdgY[0]*rBt);
    for(j=0; j < nP4; j++){
      r[j]	=s;
    }
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    ia	=2;
    for(i=1; i < nA1; i++){
      printf("\b\b\b\b%3d]",i);
      a		=A[ia];
      ia	+=2;
      ESSetSplA(a);
      splRA(&rf,&dz,rsT,rsT2a);
      splRA(&dgya,NULL,ESdgY,ESdgY2a);
      rf	=-fac/(dgya*a*rBt);
      splRA(&b,&db,ESsb,ESsb2a);
      db	*=rf;
      dz	*=rf;
      splRA(rc,drc,rcT,rcT2a);
      rc[0]	+=ESaR0[0];
      drc[0]	*=rf;
      ki	=0;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	splRA(rc+k,drc+k,rcT+ki,rcT2a+ki);
	drc[k]	*=rf;
	rct[k]	=-k*rc[k];
	splRA(rs+k,drs+k,rsT+ki,rsT2a+ki);
	drs[k]	*=rf;
	rst[k]	=k*rs[k];
      }
      for(j=0; j < nP4; j++){
	s	=dgt*(j-j0);
	ss	=0.;
	EZdra	=0.;
	EZdrt	=0.;
	for(k=1; k < ESFp1; k++){
	  rf	=k*s;
	  cs	=cos(rf);
	  sn	=sin(rf);
	  ss	+=rc[k]*cs+rs[k]*sn;
	  EZdra	+=drc[k]*cs+drs[k]*sn;
	  EZdrt	+=rct[k]*sn+rst[k]*cs;
	}
	EZdra	=drc[0]+2.*EZdra;
	EZdrt	*=2.;
	r[j]	=(rc[0]+2.*ss)*(EZdra*b*cos(s)-(dz+db*sin(s))*EZdrt);
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }
    s	=A[ia];
    a	=ESsa[ESNa]+s;
    dgya=ESdgY[ESNa]+ESdgY1a[ESNa]*s;
    rf	=-fac/(dgya*a*rBt);
    b	=ESsb[ESNa]+ESsb1a[ESNa]*s;
    db	=(ESsb1a[ESNa]+ESsb2a[ESNa]*s)*rf;
    dz	=(rsT1a[ESNa]+rsT2a[ESNa]*s)*rf;
    rc[0]	=ESaR0[0]+rcT[ESNa]+rcT1a[ESNa]*s;
    drc[0]	=(rcT1a[ESNa]+rcT2a[ESNa]*s)*rf;
    ki	=ESNa;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=rcT[ki]+rcT1a[ki]*s;
      drc[k]	=(rcT1a[ki]+rcT2a[ki]*s)*rf;
      rct[k]	=-k*rc[k];
      rs[k]	=rsT[ki]+rsT1a[ki]*s;
      drs[k]	=(rsT1a[ki]+rsT2a[ki]*s)*rf;
      rst[k]	=k*rs[k];
    }

    for(j=0; j < nP4; j++){
      s		=dgt*(j-j0);
      ss	=0.;
      EZdra	=0.;
      EZdrt	=0.;
      for(k=1; k < ESFp1; k++){
	rf	=k*s;
	cs	=cos(rf);
	sn	=sin(rf);
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
      }
      EZdra	=drc[0]+2.*EZdra;
      EZdrt	*=2.;
      r[j]	=(rc[0]+2.*ss)*(EZdra*b*cos(s)-(dz+db*sin(s))*EZdrt);
    }
    i++;
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    for(j=0; j < nP4; j++){
      r[j]	=0.;
    }
    while(i < nA5){
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
      i++;
    }
  }
  fclose(fp);
  printf("\r");
  printf("Esc2dsk%4.2f has been written\n",ttime);
  return(0);
}

int ES2eqdskS(double ttime)
{
  int nA=128,nA1=129,nA5=NSIZE,nP=128,nP1=129,nP4=132,nP5=NSIZE;
  FILE *fp;
  int i,j,ia,j0=2;
  long int nxy[10];
  char Title[160];
  double dgy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double A[2*NSIZE];
  double sp[NSIZE],Ps[NSIZE],sq[NSIZE],dsqgy[NSIZE],gy[NSIZE];
  double g[NSIZE],dggy[NSIZE],fb[NSIZE],dfbgy[NSIZE];

  double aT[NSIZE];

  rg	=1./ESRBt;
  rBt	=ESRext*rg;
  rf	=-ESgY[ESNa]*rBt;
  p2rBt	=rBt*rBt;

  sprintf(Title,"Esc2dsk%4.2f",ttime);
  sprintf(Title,"../Pst/fort.17");
  printf("Writing plasma profiles\n");

  fp	=fopen(Title,"w");
  strcpy(Title,"ESC has written this garbage");
  i	=strlen(Title);
  while(i < 159){
    Title[i]	=' ';
    i++;
  }
  Title[i]	='\0';
  g[0]	=0.;
  fwrite((void*)g,(size_t)8,(size_t)1,fp);
  fwrite((void*)Title,(size_t)8,(size_t)20,fp);
  fwrite((void*)g,(size_t)8,(size_t)1,fp);
  nxy[0]  =nP1;	/*  +1 for symmetric case */
  nxy[1]  =nA1; 
  nxy[2]  =0;
  nxy[3]  =0; 
  nxy[4]  =0; 
  fwrite((void*)nxy,(size_t)8,(size_t)5,fp);

  g[0]	=rdPV[3]-rdPV[2];
  g[1]	=zdPV[0]-zdPV[1];
  g[2]	=zdPV[0];
  g[3]	=ESaR0[0];
  g[4]	=ESRext;
  g[5]	=ESsp[0]*p2rBt;
  g[6]	=ESRBt*rBt/ESRext;
  g[7]	=EZc2gp*ESgY[ESNa]*rBt;
  g[8]	=0.;
  g[9]	=g[7];
  g[10]=0.;
  g[11]=0.;
  g[12]=1.;
  fwrite((void*)g,(size_t)8,(size_t)13,fp);

  nxy[0]	=2; 
  nxy[1]	=0; 
  nxy[2]	=0; 
  nxy[3]	=0; 
  nxy[4]	=0;
  nxy[5]	=0;
  nxy[6]	=0; 
  nxy[7]	=0;	/* symmetry index */
  nxy[8]	=nP5;	/* poloidal allocation size */
  nxy[9]	=nA5;	/* radial allocation size */
  fwrite((void*)nxy,(size_t)8,(size_t)10,fp);
  
  dgy	=1./nA;
  dgt	=EZcgp/nP;
  g[0]	=dgt;
  g[1]	=dgy;
  g[2]	=EZcgp;
  g[3]	=ESRext;
  g[4]	=ESRBt/ESRext;
  printf("dt=%15.7e %15.7e %15.7e %15.7e %15.7e\n",g[0],g[1],g[2],g[3],g[4]);
  for(i=5; i < 10; i++){
    g[i]	=0.;
  }
  fwrite((void*)g,(size_t)8,(size_t)10,fp);
  
  for(i=0; i < ESNa1; i++){
    printf("??? ESgm[%2d]=%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",i
	   ,ESgm[i],ESgm1a[i],ESgm2a[i],ESdgY[i],ESdgY1a[i],ESdgY2a[i]);
  }

  a	=0.;
  i	=0;
  ia	=0;
  A[ia]	=a;
  ia++;
  gy[i]	=ESgY[ESNa]*rBt;
  Ps[i]	=-ESPs[0]*rBt;

  sq[i]	=1./ESgm[0];
  if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
    dsqgy[i]	=ESgm2a[0]*sq[i]*sq[i]/(ESdgY[0]*rBt);
  }
  else{
    ESSetSplDCr(a);
    splRDCr2(&s,NULL,dsqgy,7);
    dsqgy[i]	*=sq[i]*sq[i]/(ESdgY[0]*rBt);
  }

  g[i]	=ESaF[0]*rg;
  splRA(&s,NULL,ESaT,ESaT2a);
  s	=ESaT[0]*rBt;
  aT[i]	=s;
  dggy[i]=-s/(ESRext*ESRext*g[i]);
  fb[i]	=-ESaR0[0]*ESdgY[0]*rBt/ESLc[0];
  dfbgy[i]	=ESaR0[0]*(ESdgY2a[0]-ESdgY[0]*ESLc2[0]/ESLc[0])/
    (ESdgY[0]*ESLc[0]);

  ss	=sqrt(dgy*(i+0.5));
  do{
    ESSetSplA(a);
    splRA(&s,&dgya,ESqgY,ESqgY2a);
    a		+=(ss-s)/dgya;
  }while(fabs(ss-s) > 1e-8);
  A[ia]	=a;
  ia++;
  ESSetSplA(a);
  splRA(g+i,&ds,ESaF,ESaF2a);
  g[i]	*=rg;
  if(ESEqSolvInPr != 2){
    splRA(sp+i,NULL,ESsp,ESsp2a);
  }
  else{
    ESSetSplDPr(a);
    splRDPr(sp+i,&ds,2);
  }
  sp[i]	*=p2rBt;
  for(i=1; i < nA1; i++){
    ss	=sqrt(dgy*i);
    do{
      ESSetSplA(a);
      splRA(&s,&dgya,ESqgY,ESqgY2a);
      a		+=(ss-s)/dgya;
    }while(fabs(ss-s) > 1e-8);
    A[ia]	=a;
    ia++;
    ESSetSplA(a);
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    splRA(&s,&ss,ESLc,ESLc2);
    fb[i]	=-ESaR0[0]*dgya/s;
    dgya	*=a;

    dfbgy[i]	=(ESaR0[0]*ds+fb[i]*ss)/(dgya*s);
    fb[i]	*=rBt;
    gy[i]	=ESgY[ESNa]*(1.-dgy*i)*rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,0);
      Ps[i]	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,1);
      break;
    default:
      splRA(Ps+i,NULL,ESPs,ESPs2a);
      break;
    }
    Ps[i]	*=-rBt;
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq[i]	=1./s;
    dsqgy[i]	=ds*sq[i]*sq[i]/(dgya*rBt);
    printf("??? dgya[%3d]=%12.4e %12.4e %12.4e\n",i,dgya,ds,dsqgy[i]);
    splRA(g+i,&ds,ESaF,ESaF2a);
    g[i]	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    aT[i]	=s;
    dggy[i]	=-s/(ESRext*ESRext*g[i]);
    if(i < nA){
      ss	=sqrt(dgy*(i+0.5));
      do{
	ESSetSplA(a);
	splRA(&s,&dgya,ESqgY,ESqgY2a);
	a		+=(ss-s)/dgya;
      }while(fabs(ss-s) > 1e-8);
      ESSetSplA(a);
      splRA(g+i,&ds,ESaF,ESaF2a);
      g[i]	*=rg;
      if(ESEqSolvInPr != 2){
	splRA(sp+i,NULL,ESsp,ESsp2a);
      }
      else{
	ESSetSplDPr(a);
	splRDPr(sp+i,&ds,2);
      }
      sp[i]	*=p2rBt;
    }
    else{
      sp[i]	=0.;
      g[i]	=(ESaF[ESNa]+ESaF1a[ESNa]*(a-ESsa[ESNa]))*rg;
    }
    A[ia]	=a;
    ia++;
  }
  ds	=dgy*ESgY[ESNa]/ESgY1a[ESNa];
  A[ia]	=ds;
  ia++;
  sp[i]	=0.;
  Ps[i]	=-(ESPs[ESNa]+ESPs1a[ESNa]*ds)*rBt;
  Ps[i]	=0.;
  sq[i]	=1./(ESgm[ESNa]+ESgm1a[ESNa]*ds);
  sq[i]	=0.;
  dsqgy[i]	=0.;
  g[i]	=(ESaF[ESNa]+ESaF1a[ESNa]*ds)*rg;
  g[i]	=0.;
  dggy[i]=0.;
  fb[i]	=0.;
  dfbgy[i]=0.;
  gy[i]	=ESgY[ESNa]*(1.-dgy*i)*rBt;

  i++;
  while(i < nA5){
    sp[i]	=0.;
    Ps[i]	=0.;
    sq[i]	=0.;
    dsqgy[i]	=0.;
    g[i]	=0.;
    dggy[i]	=0.;
    fb[i]	=0.;
    dfbgy[i]	=0.;
    gy[i]	=0.;
    i++;
  }

  fwrite((void*)sp,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)Ps,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)sq,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)dsqgy,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)g,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)dggy,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)fb,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)dfbgy,(size_t)8,(size_t)nA5,fp);
  i	=0;
  g[i]	=rf;
  s	=-0.5*nA*ESgY[ESNa]*rBt;
  dggy[i]=2.*s;
  for(i=1; i < nA1; i++){
    g[i]	=rf;
    dggy[i]	=s/i;
  }
  while(i < nA5){
    g[i]	=0.;
    dggy[i]	=0.;
    i++;
  }
  fwrite((void*)g,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)dggy,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)gy,(size_t)8,(size_t)nA5,fp);

  {
    int k,ki;
    double *r,*rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,cs,sn,b,db,dz,fac;

    r	=g;
    rc	=sp;
    rs	=Ps;
    drc	=fb;
    drs	=dfbgy;
    rct	=sq;
    rst	=dsqgy;

    printf("Writing metric tensor\n");
    /* r */
    for(j=0; j < nP4; j++){
      r[j]	=ESaR0[0];
    }
    for(j=nP4; j < nP5; j++){
      r[j]	=0.;
    }
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    ia	=2;
    for(i=1; i < nA1; i++){
      a		=A[ia];
      ia	+=2;
      ESSetSplA(a);
      ki	=0;
      splRA(rc,NULL,rcT,rcT2a);
      rc[0]	+=ESaR0[0];
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	splRA(rc+k,NULL,rcT+ki,rcT2a+ki);
	splRA(rs+k,NULL,rsT+ki,rsT2a+ki);
      }
      for(j=0; j < nP4; j++){
	s	=dgt*(j-j0);
	ss	=0.;
	for(k=1; k < ESFp1; k++){
	  rf	=k*s;
	  ss	+=rc[k]*cos(rf)+rs[k]*sin(rf);
	}
	r[j]	=rc[0]+2.*ss;
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }
    rf	=A[ia];
    rc[0]	=ESaR0[0]+rcT[ESNa]+rcT1a[ESNa]*rf;
    ki	=ESNa;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=rcT[ki]+rcT1a[ki]*rf;
      rs[k]	=rsT[ki]+rsT1a[ki]*rf;
    }
    for(j=0; j < nP4; j++){
      s		=dgt*(j-j0);
      ss	=0.;
      for(k=1; k < ESFp1; k++){
	rf	=k*s;
	ss	+=rc[k]*cos(rf)+rs[k]*sin(rf);
      }
      r[j]	=rc[0]+2.*ss;
    }
    i++;
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    for(j=0; j < nP4; j++){
      r[j]	=0.;
    }
    while(i < nA5){
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
      i++;
    }
    /* z */
    for(j=0; j < nP4; j++){
      r[j]	=ESaZ0[0];
    }
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    ia	=2;
    for(i=1; i < nA1; i++){
      a		=A[ia];
      ia	+=2;
      ESSetSplA(a);
      splRA(rs,NULL,rsT,rsT2a);
      rs[0]	+=ESaZ0[0];
      splRA(&b,NULL,ESsb,ESsb2a);
      for(j=0; j < nP4; j++){
	r[j]	=rs[0]+b*sin(dgt*(j-j0));
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }
    b	=ESsb[ESNa]+ESsb1a[ESNa]*A[ia];
    rs[0]=ESaZ0[0]+rsT[ESNa]+rsT1a[ESNa]*A[ia];
    for(j=0; j < nP4; j++){
      r[j]	=rs[0]+b*sin(dgt*(j-j0));
    }
    i++;
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    for(j=0; j < nP4; j++){
      r[j]	=0.;
    }
    while(i < nA5){
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
      i++;
    }
    
    fac	=-ESgY[ESNa]*rBt;
    /* qg05 */
    s	=-ESaR0[0]*ESDPc[0]*fac/(ESdgY[0]*rBt);
    for(j=0; j < nP4; j++){
      r[j]	=s;
    }
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    ia	=2;
    for(i=1; i < nA1; i++){
      a		=A[ia];
      ia	+=2;
      ESSetSplA(a);
      splRA(&rf,&dz,rsT,rsT2a);
      splRA(&dgya,NULL,ESdgY,ESdgY2a);
      rf	=-fac/(dgya*a*rBt);
      splRA(&b,&db,ESsb,ESsb2a);
      db	*=rf;
      dz	*=rf;
      splRA(rc,drc,rcT,rcT2a);
      rc[0]	+=ESaR0[0];
      drc[0]	*=rf;
      ki	=0;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	splRA(rc+k,drc+k,rcT+ki,rcT2a+ki);
	drc[k]	*=rf;
	rct[k]	=-k*rc[k];
	splRA(rs+k,drs+k,rsT+ki,rsT2a+ki);
	drs[k]	*=rf;
	rst[k]	=k*rs[k];
      }
      for(j=0; j < nP4; j++){
	s	=dgt*(j-j0-0.5);
	ss	=0.;
	EZdra	=0.;
	EZdrt	=0.;
	for(k=1; k < ESFp1; k++){
	  rf	=k*s;
	  cs	=cos(rf);
	  sn	=sin(rf);
	  ss	+=rc[k]*cs+rs[k]*sn;
	  EZdra	+=drc[k]*cs+drs[k]*sn;
	  EZdrt	+=rct[k]*sn+rst[k]*cs;
	}
	EZdra	=drc[0]+2.*EZdra;
	EZdrt	*=2.;
	r[j]	=(rc[0]+2.*ss)*(EZdra*b*cos(s)-(dz+db*sin(s))*EZdrt);
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }
    s	=A[ia];
    a	=ESsa[ESNa]+s;
    dgya=ESdgY[ESNa]+ESdgY1a[ESNa]*s;
    rf	=-fac/(dgya*a*rBt);
    b	=ESsb[ESNa]+ESsb1a[ESNa]*s;
    db	=(ESsb1a[ESNa]+ESsb2a[ESNa]*s)*rf;
    dz	=(rsT1a[ESNa]+rsT2a[ESNa]*s)*rf;
    rc[0]	=ESaR0[0]+rcT[ESNa]+rcT1a[ESNa]*s;
    drc[0]	=(rcT1a[ESNa]+rcT2a[ESNa]*s)*rf;
    ki	=ESNa;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=rcT[ki]+rcT1a[ki]*s;
      drc[k]	=(rcT1a[ki]+rcT2a[ki]*s)*rf;
      rct[k]	=-k*rc[k];
      rs[k]	=rsT[ki]+rsT1a[ki]*s;
      drs[k]	=(rsT1a[ki]+rsT2a[ki]*s)*rf;
      rst[k]	=k*rs[k];
    }
    for(j=0; j < nP4; j++){
      s		=dgt*(j-j0-0.5);
      ss	=0.;
      EZdra	=0.;
      EZdrt	=0.;
      for(k=1; k < ESFp1; k++){
	rf	=k*s;
	cs	=cos(rf);
	sn	=sin(rf);
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
      }
      EZdra	=drc[0]+2.*EZdra;
      EZdrt	*=2.;
      r[j]	=(rc[0]+2.*ss)*(EZdra*b*cos(s)-(dz+db*sin(s))*EZdrt);
    }
    i++;
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    for(j=0; j < nP4; j++){
      r[j]	=0.;
    }
    while(i < nA5){
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
      i++;
    }
    /* qg */
    s	=-ESaR0[0]*ESDPc[0]*fac/(ESdgY[0]*rBt);
    for(j=0; j < nP4; j++){
      r[j]	=s;
    }
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    ia	=2;
    for(i=1; i < nA1; i++){
      a		=A[ia];
      ia	+=2;
      ESSetSplA(a);
      splRA(&rf,&dz,rsT,rsT2a);
      splRA(&dgya,NULL,ESdgY,ESdgY2a);
      rf	=-fac/(dgya*a*rBt);
      splRA(&b,&db,ESsb,ESsb2a);
      db	*=rf;
      dz	*=rf;
      splRA(rc,drc,rcT,rcT2a);
      rc[0]	+=ESaR0[0];
      drc[0]	*=rf;
      ki	=0;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	splRA(rc+k,drc+k,rcT+ki,rcT2a+ki);
	drc[k]	*=rf;
	rct[k]	=-k*rc[k];
	splRA(rs+k,drs+k,rsT+ki,rsT2a+ki);
	drs[k]	*=rf;
	rst[k]	=k*rs[k];
      }
      for(j=0; j < nP4; j++){
	s	=dgt*(j-j0);
	ss	=0.;
	EZdra	=0.;
	EZdrt	=0.;
	for(k=1; k < ESFp1; k++){
	  rf	=k*s;
	  cs	=cos(rf);
	  sn	=sin(rf);
	  ss	+=rc[k]*cs+rs[k]*sn;
	  EZdra	+=drc[k]*cs+drs[k]*sn;
	  EZdrt	+=rct[k]*sn+rst[k]*cs;
	}
	EZdra	=drc[0]+2.*EZdra;
	EZdrt	*=2.;
	r[j]	=(rc[0]+2.*ss)*(EZdra*b*cos(s)-(dz+db*sin(s))*EZdrt);
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }
    s	=A[ia];
    a	=ESsa[ESNa]+s;
    dgya=ESdgY[ESNa]+ESdgY1a[ESNa]*s;
    rf	=-fac/(dgya*a*rBt);
    b	=ESsb[ESNa]+ESsb1a[ESNa]*s;
    db	=(ESsb1a[ESNa]+ESsb2a[ESNa]*s)*rf;
    dz	=(rsT1a[ESNa]+rsT2a[ESNa]*s)*rf;
    rc[0]	=ESaR0[0]+rcT[ESNa]+rcT1a[ESNa]*s;
    drc[0]	=(rcT1a[ESNa]+rcT2a[ESNa]*s)*rf;
    ki	=ESNa;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=rcT[ki]+rcT1a[ki]*s;
      drc[k]	=(rcT1a[ki]+rcT2a[ki]*s)*rf;
      rct[k]	=-k*rc[k];
      rs[k]	=rsT[ki]+rsT1a[ki]*s;
      drs[k]	=(rsT1a[ki]+rsT2a[ki]*s)*rf;
      rst[k]	=k*rs[k];
    }

    for(j=0; j < nP4; j++){
      s		=dgt*(j-j0);
      ss	=0.;
      EZdra	=0.;
      EZdrt	=0.;
      for(k=1; k < ESFp1; k++){
	rf	=k*s;
	cs	=cos(rf);
	sn	=sin(rf);
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
      }
      EZdra	=drc[0]+2.*EZdra;
      EZdrt	*=2.;
      r[j]	=(rc[0]+2.*ss)*(EZdra*b*cos(s)-(dz+db*sin(s))*EZdrt);
    }
    i++;
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    for(j=0; j < nP4; j++){
      r[j]	=0.;
    }
    while(i < nA5){
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
      i++;
    }
  }
  printf("Esc2dsk%4.2f has been written\n",ttime);
  fclose(fp);
  return(0);
}

int ES2eqdskBoozer(double ttime)
{
  static int Fl=0;
  int nA=128,nA1=129,nA5=NSIZE,nP=256,nP1=257,nP4=260,nP5=NSIZE;
  FILE *fp;
  int i,j,ia,j0=2;
  long int nxy[10];
  char Title[160];
  double dgy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double A[2*NSIZE];
  double sp[NSIZE],Ps[NSIZE],sq[NSIZE],dsqgy[NSIZE],gy[NSIZE];
  double g[NSIZE],dggy[NSIZE],fb[NSIZE],dfbgy[NSIZE];

  double aT[NSIZE];

  rg	=1./ESRBt;
  rBt	=ESRext*rg;
  rf	=-ESgY[ESNa]*rBt;
  p2rBt	=rBt*rBt;

  sprintf(Title,"../Pst/fort.17");
  sprintf(Title,"Esc2dsk%4.2f",ttime);
  sprintf(Title,"ES2eqdsk.%2.2d",Fl);
  printf("Writing plasma profiles\n");

  fp	=fopen(Title,"w");
  strcpy(Title,"ESC has written this garbage");
  i	=strlen(Title);
  while(i < 159){
    Title[i]	=' ';
    i++;
  }
  Title[i]	='\0';
  g[0]	=0.;
  fwrite((void*)g,(size_t)8,(size_t)1,fp);
  fwrite((void*)Title,(size_t)8,(size_t)20,fp);
  fwrite((void*)g,(size_t)8,(size_t)1,fp);
  nxy[0]  =nP;
  nxy[1]  =nA1; 
  nxy[2]  =0;
  nxy[3]  =0; 
  nxy[4]  =0; 
  fwrite((void*)nxy,(size_t)8,(size_t)5,fp);

  g[0]	=rdPV[3]-rdPV[2];
  g[1]	=zdPV[0]-zdPV[1];
  g[2]	=zdPV[0];
  g[3]	=ESaR0[0];
  g[4]	=ESRext;
  g[5]	=ESsp[0]*p2rBt;
  g[6]	=ESRBt*rBt/ESRext;
  g[7]	=EZc2gp*ESgY[ESNa]*rBt;
  g[8]	=0.;
  g[9]	=g[7];
  g[10]=0.;
  g[11]=0.;
  g[12]=1.;
  fwrite((void*)g,(size_t)8,(size_t)13,fp);

  nxy[0]	=2; 
  nxy[1]	=0; 
  nxy[2]	=0; 
  nxy[3]	=0; 
  nxy[4]	=2;
  nxy[5]	=0;
  nxy[6]	=1; 
  nxy[7]	=0;	/* symmetry index */
  nxy[8]	=nP5;	/* poloidal allocation size */
  nxy[9]	=nA5;	/* radial allocation size */
  fwrite((void*)nxy,(size_t)8,(size_t)10,fp);
  
  dgt	=EZc2gp/nP;
  dgy	=1./nA;
  g[0]	=dgt;
  g[1]	=dgy;
  g[2]	=EZcgp;
  g[3]	=ESRext;
  g[4]	=ESRBt/ESRext;
  for(i=5; i < 10; i++){
    g[i]	=0.;
  }
  fwrite((void*)g,(size_t)8,(size_t)10,fp);

  a	=0.;
  i	=0;
  ia	=0;
  A[ia]	=a;
  ia++;
  gy[i]	=ESgY[ESNa]*rBt;
  Ps[i]	=-ESPs[0]*rBt;

  sq[i]	=1./ESgm[0];
  if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
    dsqgy[i]	=ESgm2a[0]*sq[i]*sq[i]/(ESdgY[0]*rBt);
  }
  else{
    ESSetSplDCr(a);
    splRDCr2(&s,NULL,dsqgy,7);
    dsqgy[i]	*=sq[i]*sq[i]/(ESdgY[0]*rBt);
  }

  g[i]	=ESaF[0]*rg;
  splRA(&s,NULL,ESaT,ESaT2a);
  s	=ESaT[0]*rBt;
  aT[i]	=s;
  dggy[i]=-s/(ESRext*ESRext*g[i]);
  fb[i]	=-ESaR0[0]*ESdgY[0]*rBt/ESLc[0];
  dfbgy[i]	=ESaR0[0]*(ESdgY2a[0]-ESdgY[0]*ESLc2[0]/ESLc[0])/
    (ESdgY[0]*ESLc[0]);

  ss	=sqrt(dgy*(i+0.5));
  do{
    ESSetSplA(a);
    splRA(&s,&dgya,ESqgY,ESqgY2a);
    a		+=(ss-s)/dgya;
  }while(fabs(ss-s) > 1e-8);
  A[ia]	=a;
  ia++;
  ESSetSplA(a);
  splRA(g+i,&ds,ESaF,ESaF2a);
  g[i]	*=rg;
  if(ESEqSolvInPr != 2){
    splRA(sp+i,NULL,ESsp,ESsp2a);
  }
  else{
    ESSetSplDPr(a);
    splRDPr(sp+i,&ds,2);
  }
  sp[i]	*=p2rBt;
  for(i=1; i < nA1; i++){
    ss	=sqrt(dgy*i);
    do{
      ESSetSplA(a);
      splRA(&s,&dgya,ESqgY,ESqgY2a);
      a		+=(ss-s)/dgya;
    }while(fabs(ss-s) > 1e-8);
    A[ia]	=a;
    ia++;
    ESSetSplA(a);
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    splRA(&s,&ss,ESLc,ESLc2);
    fb[i]	=-ESaR0[0]*dgya/s;
    dgya	*=a;
    dfbgy[i]	=(ESaR0[0]*ds+fb[i]*ss)/(dgya*s);
    fb[i]	*=rBt;
    gy[i]	=ESgY[ESNa]*(1.-dgy*i)*rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,0);
      Ps[i]	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Ps+i,&ds,1);
      break;
    default:
      splRA(Ps+i,NULL,ESPs,ESPs2a);
      break;
    }
    Ps[i]	*=-rBt;
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq[i]	=1./s;
    dsqgy[i]	=ds*sq[i]*sq[i]/(dgya*rBt);
    splRA(g+i,&ds,ESaF,ESaF2a);
    g[i]	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    aT[i]	=s;
    dggy[i]	=-s/(ESRext*ESRext*g[i]);
    if(i < nA){
      ss	=sqrt(dgy*(i+0.5));
      do{
	ESSetSplA(a);
	splRA(&s,&dgya,ESqgY,ESqgY2a);
	a		+=(ss-s)/dgya;
      }while(fabs(ss-s) > 1e-8);
      ESSetSplA(a);
      splRA(g+i,&ds,ESaF,ESaF2a);
      g[i]	*=rg;
      if(ESEqSolvInPr != 2){
	splRA(sp+i,NULL,ESsp,ESsp2a);
      }
      else{
	ESSetSplDPr(a);
	splRDPr(sp+i,&ds,2);
      }
      sp[i]	*=p2rBt;
    }
    else{
      sp[i]	=0.;
      g[i]	=(ESaF[ESNa]+ESaF1a[ESNa]*(a-ESsa[ESNa]))*rg;
    }
    A[ia]	=a;
    ia++;
  }
  ds	=dgy*ESgY[ESNa]/ESgY1a[ESNa];
  A[ia]	=ds;
  ia++;
  sp[i]	=0.;
  Ps[i]	=-(ESPs[ESNa]+ESPs1a[ESNa]*ds)*rBt;
  Ps[i]	=0.;
  sq[i]	=1./(ESgm[ESNa]+ESgm1a[ESNa]*ds);
  sq[i]	=0.;
  dsqgy[i]	=0.;
  g[i]	=(ESaF[ESNa]+ESaF1a[ESNa]*ds)*rg;
  g[i]	=0.;
  dggy[i]=0.;
  fb[i]	=0.;
  dfbgy[i]=0.;
  gy[i]	=ESgY[ESNa]*(1.-dgy*i)*rBt;

  i++;
  while(i < nA5){
    sp[i]	=0.;
    Ps[i]	=0.;
    sq[i]	=0.;
    dsqgy[i]	=0.;
    g[i]	=0.;
    dggy[i]	=0.;
    fb[i]	=0.;
    dfbgy[i]	=0.;
    gy[i]	=0.;
    i++;
  }
  fwrite((void*)sp,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)Ps,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)sq,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)dsqgy,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)g,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)dggy,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)fb,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)dfbgy,(size_t)8,(size_t)nA5,fp);
  i	=0;
  g[i]	=rf;
  s	=-0.5*nA*ESgY[ESNa]*rBt;
  dggy[i]=2.*s;
  for(i=1; i < nA1; i++){
    g[i]	=rf;
    dggy[i]	=s/i;
  }
  while(i < nA5){
    g[i]	=0.;
    dggy[i]	=0.;
    i++;
  }
  fwrite((void*)g,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)dggy,(size_t)8,(size_t)nA5,fp);
  fwrite((void*)gy,(size_t)8,(size_t)nA5,fp);

  {
    int k,ki,jj,kj;
    double *r,*rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,cs,sn,b,db,dz,fac;
    double aF,gh,t,qg0,qg1;
#ifdef XWIN
    double Lc[ES2Mp1],Ls[ES2Mp1];
    double Gh[ESNp1],Gh2[ESNp1];
#else
    double Lc[33],Ls[33];
    double Gh[129],Gh2[129];
#endif

    r	=g;
    rc	=sp;
    rs	=Ps;
    drc	=fb;
    drs	=dfbgy;
    rct	=sq;
    rst	=dsqgy;

    printf("Writing geometry ->r[  0]");
    /* r */
    for(j=0; j < nP4; j++){
      r[j]	=ESaR0[0];
    }
    for(j=nP4; j < nP5; j++){
      r[j]	=0.;
    }
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    ia	=2;
    for(i=1; i < nA1; i++){
      printf("\b\b\b\b%3d]",i);
      fflush(stdout);
      a		=A[ia];
      ia	+=2;
      ESSetSplA(a);

      splRA(&aF,NULL,ESaF,ESaF2a);
      splRA(&dgya,NULL,ESdgY,ESdgY2a);
      dgya	*=a;
      ss	=dgya/aF;
      ss	*=ss*ESaR0[0];
      splRA(Lc,NULL,ESLc,ESLc2);
      splRA(&s,NULL,ESg22c,ESg22c2);
      Lc[0]	=2./(Lc[0]+s*ss);
      ki	=0;
      for(k=1; k < ES2Mp1; k++){
	ki	+=ESNa1;
	splRA(Lc+k,NULL,ESLc+ki,ESLc2+ki);
	splRA(Ls+k,NULL,ESLs+ki,ESLs2+ki);
	splRA(&s,NULL,ESg22c+ki,ESg22c2+ki);
	Lc[k]	=(Lc[k]+s*ss)*Lc[0];
	splRA(&s,NULL,ESg22s+ki,ESg22s2+ki);
	Ls[k]	=(Ls[k]+s*ss)*Lc[0];
      }
      for(j=0; j < ESNp; j++){
	s	=0.;
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

      ki	=0;
      splRA(rc,NULL,rcT,rcT2a);
      rc[0]	+=ESaR0[0];
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	splRA(rc+k,NULL,rcT+ki,rcT2a+ki);
	splRA(rs+k,NULL,rsT+ki,rsT2a+ki);
      }
      t	=-dgt*j0;
      for(j=0; j < nP4; j++){
	gh	=dgt*(j-j0);
	do{
	  ESSetSplP(t);
	  splRP(&s,&ss,Gh,Gh2);
	  s      +=t-gh;
	  t	-=s/(1.+ss);
	}while(fabs(s) > 1e-8);
	s	=t;
	ss	=0.;
	for(k=1; k < ESFp1; k++){
	  rf	=k*s;
	  ss	+=rc[k]*cos(rf)+rs[k]*sin(rf);
	}
	r[j]	=rc[0]+2.*ss;
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }
    rf	=A[ia];
    dgya=(ESdgY[ESNa]+ESdgY1a[ESNa]*rf)*(ESsa[ESNa]+rf);
    aF	=ESaF[ESNa]+ESaF1a[ESNa]*rf;
    ss	=dgya/aF;
    ss	*=ss*ESaR0[0];
    Lc[0]	=ESLc[ESNa]+ESLc1[ESNa]*rf;
    s		=ESg22c[ESNa]+ESg22c1[ESNa]*rf;
    Lc[0]	=2./(Lc[0]+s*ss);
    ki	=ESNa;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      Lc[k]	=ESLc[ki]+ESLc1[ki]*rf;
      Ls[k]	=ESLs[ki]+ESLs1[ki]*rf;
      s		=ESg22c[ki]+ESg22c1[ki]*rf;
      Lc[k]	=(Lc[k]+s*ss)*Lc[0];
      s		=ESg22s[ki]+ESg22s1[ki]*rf;
      Ls[k]	=(Ls[k]+s*ss)*Lc[0];
    }
    for(j=0; j < ESNp; j++){
      s	=0.;
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

    rc[0]	=ESaR0[0]+rcT[ESNa]+rcT1a[ESNa]*rf;
    ki	=ESNa;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=rcT[ki]+rcT1a[ki]*rf;
      rs[k]	=rsT[ki]+rsT1a[ki]*rf;
    }
    t	=-dgt*j0;
    for(j=0; j < nP4; j++){
      gh	=dgt*(j-j0);
      do{
	ESSetSplP(t);
	splRP(&s,&ss,Gh,Gh2);
	s      +=t-gh;
	t	-=s/(1.+ss);
      }while(fabs(s) > 1e-8);
      s	=t;
      ss	=0.;
      for(k=1; k < ESFp1; k++){
	rf	=k*s;
	ss	+=rc[k]*cos(rf)+rs[k]*sin(rf);
      }
      r[j]	=rc[0]+2.*ss;
    }
    i++;
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    for(j=0; j < nP4; j++){
      r[j]	=0.;
    }
    while(i < nA5){
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
      i++;
    }

    /* z */
    printf("\b\b\b\b\b\bz[%3d]",0);
    for(j=0; j < nP4; j++){
      r[j]	=ESaZ0[0];
    }
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    ia	=2;
    for(i=1; i < nA1; i++){
      printf("\b\b\b\b%3d]",i);
      fflush(stdout);
      a		=A[ia];
      ia	+=2;
      ESSetSplA(a);
      splRA(&dgya,NULL,ESdgY,ESdgY2a);
      dgya	*=a;
      splRA(&aF,NULL,ESaF,ESaF2a);
      ss	=dgya/aF;
      ss	*=ss*ESaR0[0];
      splRA(Lc,NULL,ESLc,ESLc2);
      splRA(&s,NULL,ESg22c,ESg22c2);
      Lc[0]	=2./(Lc[0]+s*ss);
      ki	=0;
      for(k=1; k < ES2Mp1; k++){
	ki	+=ESNa1;
	splRA(Lc+k,NULL,ESLc+ki,ESLc2+ki);
	splRA(Ls+k,NULL,ESLs+ki,ESLs2+ki);
	splRA(&s,NULL,ESg22c+ki,ESg22c2+ki);
	Lc[k]	=(Lc[k]+s*ss)*Lc[0];
	splRA(&s,NULL,ESg22s+ki,ESg22s2+ki);
	Ls[k]	=(Ls[k]+s*ss)*Lc[0];
      }
      for(j=0; j < ESNp; j++){
	s	=0.;
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
      splRA(rs,NULL,rsT,rsT2a);
      rs[0]	+=ESaZ0[0];
      splRA(&b,NULL,ESsb,ESsb2a);
      t	=-dgt*j0;
      for(j=0; j < nP4; j++){
	gh	=dgt*(j-j0);
	do{
	  ESSetSplP(t);
	  splRP(&s,&ss,Gh,Gh2);
	  s      +=t-gh;
	  t	-=s/(1.+ss);
	}while(fabs(s) > 1e-8);
	r[j]	=rs[0]+b*sin(t);
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }

    rf	=A[ia];
    dgya=(ESdgY[ESNa]+ESdgY1a[ESNa]*rf)*(ESsa[ESNa]+rf);
    aF	=ESaF[ESNa]+ESaF1a[ESNa]*rf;
    ss	=dgya/aF;
    ss	*=ss*ESaR0[0];
    Lc[0]	=ESLc[ESNa]+ESLc1[ESNa]*rf;
    s		=ESg22c[ESNa]+ESg22c1[ESNa]*rf;
    Lc[0]	=2./(Lc[0]+s*ss);
    ki	=ESNa;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      Lc[k]	=ESLc[ki]+ESLc1[ki]*rf;
      Ls[k]	=ESLs[ki]+ESLs1[ki]*rf;
      s		=ESg22c[ki]+ESg22c1[ki]*rf;
      Lc[k]	=(Lc[k]+s*ss)*Lc[0];
      s		=ESg22s[ki]+ESg22s1[ki]*rf;
      Ls[k]	=(Ls[k]+s*ss)*Lc[0];
    }
    for(j=0; j < ESNp; j++){
      s	=0.;
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

    b	=ESsb[ESNa]+ESsb1a[ESNa]*rf;
    rs[0]=ESaZ0[0]+rsT[ESNa]+rsT1a[ESNa]*rf;
    t	=-dgt*j0;
    for(j=0; j < nP4; j++){
      gh	=dgt*(j-j0);
      do{
	ESSetSplP(t);
	splRP(&s,&ss,Gh,Gh2);
	s      +=t-gh;
	t	-=s/(1.+ss);
      }while(fabs(s) > 1e-8);
      r[j]	=rs[0]+b*sin(t);
    }
    i++;
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    for(j=0; j < nP4; j++){
      r[j]	=0.;
    }
    while(i < nA5){
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
      i++;
    }
    
    /* qg05 */
    printf("\b\b\b\b\b\bqg05[%3d]",0);
    s	=ESgY[ESNa]*ESaR0[0]*ESDPc[0]/ESdgY[0];
    for(j=0; j < nP4; j++){
      r[j]	=s;
    }
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    ia	=1;
    for(i=1; i < nA1; i++){
      printf("\b\b\b\b%3d]",i);
      fflush(stdout);
      a		=A[ia];
      ia	+=2;
      ESSetSplA(a);
      splRA(&dgya,NULL,ESdgY,ESdgY2a);
      dgya	*=a;
      splRA(&aF,NULL,ESaF,ESaF2a);
      ss	=dgya/aF;
      ss	*=ss*ESaR0[0];
      splRA(Lc,NULL,ESLc,ESLc2);
      splRA(&s,NULL,ESg22c,ESg22c2);
      Lc[0]	=2./(Lc[0]+s*ss);
      ki	=0;
      for(k=1; k < ES2Mp1; k++){
	ki	+=ESNa1;
	splRA(Lc+k,NULL,ESLc+ki,ESLc2+ki);
	splRA(Ls+k,NULL,ESLs+ki,ESLs2+ki);
	splRA(&s,NULL,ESg22c+ki,ESg22c2+ki);
	Lc[k]	=(Lc[k]+s*ss)*Lc[0];
	splRA(&s,NULL,ESg22s+ki,ESg22s2+ki);
	Ls[k]	=(Ls[k]+s*ss)*Lc[0];
      }
      splRA(&rf,&dz,rsT,rsT2a);
      splRA(&b,&db,ESsb,ESsb2a);
      splRA(rc,drc,rcT,rcT2a);
      rc[0]	+=ESaR0[0];
      ki	=0;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	splRA(rc+k,drc+k,rcT+ki,rcT2a+ki);
	rct[k]	=-k*rc[k];
	splRA(rs+k,drs+k,rsT+ki,rsT2a+ki);
	rst[k]	=k*rs[k];
      }
      splRA(&fac,NULL,ESDPc,ESDPc2);
      fac	*=ESgY[ESNa]/dgya;
      qg0	=0.;
      for(j=0; j < ESNp; j++){
	s	=0.;
	ss	=0.;
	EZdra	=0.;
	EZdrt	=0.;
	kj	=0;
	for(k=1; k < ES2Mp1; k++){
	  kj	+=j;
	  if(kj >= ESNp){
	    kj	-=ESNp;
	  }
	  s	+=(Lc[k]*ESsn1[kj]+Ls[k]*(1.-EScs1[kj]))/k;
	  ss	+=rc[k]*EScs1[kj]+rs[k]*ESsn1[kj];
	  EZdra	+=drc[k]*EScs1[kj]+drs[k]*ESsn1[kj];
	  EZdrt	+=rct[k]*ESsn1[kj]+rst[k]*EScs1[kj];
	}
	Gh[j]	=s;
	EZdra	=drc[0]+2.*EZdra;
	EZdrt	*=2.;
	ss	=2.*ss+rc[0];
	s	=b*EScs1[j];
	EZdra	=dgya/(ss*(EZdra*s-(dz+db*ESsn1[j])*EZdrt));
	ss	=aF/ss;
	qg0	+=1./((EZdrt*EZdrt+s*s)*EZdra*EZdra+ss*ss);
      }
      fac	*=ESNp/qg0;
      Gh[j]	=Gh[0];
      splP(Gh,Gh2);

      t	=-dgt*j0;
      for(j=0; j < nP4; j++){
	gh	=dgt*(j-j0);
	do{
	  ESSetSplP(t);
	  splRP(&s,&ss,Gh,Gh2);
	  s      +=t-gh;
	  t	-=s/(1.+ss);
	}while(fabs(s) > 1e-8);
	ss	=0.;
	EZdra	=0.;
	EZdrt	=0.;
	for(k=1; k < ESFp1; k++){
	  rf	=k*t;
	  cs	=cos(rf);
	  sn	=sin(rf);
	  ss	+=rc[k]*cs+rs[k]*sn;
	  EZdra	+=drc[k]*cs+drs[k]*sn;
	  EZdrt	+=rct[k]*sn+rst[k]*cs;
	}
	EZdra	=drc[0]+2.*EZdra;
	EZdrt	*=2.;

	ss	=2.*ss+rc[0];
	s	=b*cos(t);
	EZdra	=dgya/(ss*(EZdra*s-(dz+db*sin(t))*EZdrt));
	ss	=aF/ss;
	r[j]	=fac/((EZdrt*EZdrt+s*s)*EZdra*EZdra+ss*ss);
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }
    rf	=A[ia];
    dgya=(ESdgY[ESNa]+ESdgY1a[ESNa]*rf)*(ESsa[ESNa]+rf);
    aF	=ESaF[ESNa]+ESaF1a[ESNa]*rf;
    ss	=dgya/aF;
    ss	*=ss*ESaR0[0];
    Lc[0]	=ESLc[ESNa]+ESLc1[ESNa]*rf;
    s		=ESg22c[ESNa]+ESg22c1[ESNa]*rf;
    Lc[0]	=2./(Lc[0]+s*ss);
    ki	=ESNa;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      Lc[k]	=ESLc[ki]+ESLc1[ki]*rf;
      Ls[k]	=ESLs[ki]+ESLs1[ki]*rf;
      s		=ESg22c[ki]+ESg22c1[ki]*rf;
      Lc[k]	=(Lc[k]+s*ss)*Lc[0];
      s		=ESg22s[ki]+ESg22s1[ki]*rf;
      Ls[k]	=(Ls[k]+s*ss)*Lc[0];
    }
    s	=A[ia];
    a	=ESsa[ESNa]+s;
    dgya=(ESdgY[ESNa]+ESdgY1a[ESNa]*s)*a;
    fac	=(ESDPc[ESNa]+ESDPc1[ESNa]*s)*ESgY[ESNa]/dgya;
    b	=ESsb[ESNa]+ESsb1a[ESNa]*s;
    db	=ESsb1a[ESNa]+ESsb2a[ESNa]*s;
    dz	=rsT1a[ESNa]+rsT2a[ESNa]*s;
    rc[0]	=ESaR0[0]+rcT[ESNa]+rcT1a[ESNa]*s;
    drc[0]	=rcT1a[ESNa]+rcT2a[ESNa]*s;
    ki	=ESNa;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=rcT[ki]+rcT1a[ki]*s;
      drc[k]	=rcT1a[ki]+rcT2a[ki]*s;
      rct[k]	=-k*rc[k];
      rs[k]	=rsT[ki]+rsT1a[ki]*s;
      drs[k]	=rsT1a[ki]+rsT2a[ki]*s;
      rst[k]	=k*rs[k];
    }
    qg0	=0.;
    for(j=0; j < ESNp; j++){
      s	=0.;
      ss	=0.;
      EZdra	=0.;
      EZdrt	=0.;
      kj	=0;
      for(k=1; k < ES2Mp1; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	s	+=(Lc[k]*ESsn1[kj]+Ls[k]*(1.-EScs1[kj]))/k;
	ss	+=rc[k]*EScs1[kj]+rs[k]*ESsn1[kj];
	EZdra	+=drc[k]*EScs1[kj]+drs[k]*ESsn1[kj];
	EZdrt	+=rct[k]*ESsn1[kj]+rst[k]*EScs1[kj];
      }
      Gh[j]	=s;
      EZdra	=drc[0]+2.*EZdra;
      EZdrt	*=2.;
      ss	=2.*ss+rc[0];
      s		=b*EScs1[j];
      EZdra	=dgya/(ss*(EZdra*s-(dz+db*ESsn1[j])*EZdrt));
      ss	=aF/ss;
      qg0	+=1./((EZdrt*EZdrt+s*s)*EZdra*EZdra+ss*ss);
    }
    fac		*=ESNp/qg0;
    Gh[j]	=Gh[0];
    splP(Gh,Gh2);

    t	=-dgt*j0;
    for(j=0; j < nP4; j++){
      gh	=dgt*(j-j0);
      do{
	ESSetSplP(t);
	splRP(&s,&ss,Gh,Gh2);
	s      +=t-gh;
	t	-=s/(1.+ss);
      }while(fabs(s) > 1e-8);
      ss	=0.;
      EZdra	=0.;
      EZdrt	=0.;
      for(k=1; k < ESFp1; k++){
	rf	=k*t;
	cs	=cos(rf);
	sn	=sin(rf);
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
      }
      EZdra	=drc[0]+2.*EZdra;
      EZdrt	*=2.;
      ss	=2.*ss+rc[0];
      s		=b*cos(t);
      EZdra	=dgya/(ss*(EZdra*s-(dz+db*sin(t))*EZdrt));
      ss	=aF/ss;
      r[j]	=fac/((EZdrt*EZdrt+s*s)*EZdra*EZdra+ss*ss);
    }
    i++;
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    for(j=0; j < nP4; j++){
      r[j]	=0.;
    }
    while(i < nA5){
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
      i++;
    }

    /* qg */
    printf("\b\b\b\b\b\b\b\b\bqg[%3d]",0);
    s	=ESgY[ESNa]*ESaR0[0]*ESDPc[0]/ESdgY[0];
    for(j=0; j < nP4; j++){
      r[j]	=s;
    }
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    ia	=2;
    for(i=1; i < nA1; i++){
      printf("\b\b\b\b%3d]",i);
      fflush(stdout);
      a		=A[ia];
      ia	+=2;
      ESSetSplA(a);
      splRA(&dgya,NULL,ESdgY,ESdgY2a);
      dgya	*=a;
      splRA(&aF,NULL,ESaF,ESaF2a);
      ss	=dgya/aF;
      ss	*=ss*ESaR0[0];
      splRA(Lc,NULL,ESLc,ESLc2);
      splRA(&s,NULL,ESg22c,ESg22c2);
      Lc[0]	=2./(Lc[0]+s*ss);
      ki	=0;
      for(k=1; k < ES2Mp1; k++){
	ki	+=ESNa1;
	splRA(Lc+k,NULL,ESLc+ki,ESLc2+ki);
	splRA(Ls+k,NULL,ESLs+ki,ESLs2+ki);
	splRA(&s,NULL,ESg22c+ki,ESg22c2+ki);
	Lc[k]	=(Lc[k]+s*ss)*Lc[0];
	splRA(&s,NULL,ESg22s+ki,ESg22s2+ki);
	Ls[k]	=(Ls[k]+s*ss)*Lc[0];
      }
      splRA(&rf,&dz,rsT,rsT2a);
      splRA(&b,&db,ESsb,ESsb2a);
      splRA(rc,drc,rcT,rcT2a);
      rc[0]	+=ESaR0[0];
      ki	=0;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	splRA(rc+k,drc+k,rcT+ki,rcT2a+ki);
	rct[k]	=-k*rc[k];
	splRA(rs+k,drs+k,rsT+ki,rsT2a+ki);
	rst[k]	=k*rs[k];
      }
      splRA(&fac,NULL,ESDPc,ESDPc2);
      fac	*=ESgY[ESNa]/dgya;
      qg0	=0.;
      for(j=0; j < ESNp; j++){
	s	=0.;
	ss	=0.;
	EZdra	=0.;
	EZdrt	=0.;
	kj	=0;
	for(k=1; k < ES2Mp1; k++){
	  kj	+=j;
	  if(kj >= ESNp){
	    kj	-=ESNp;
	  }
	  s	+=(Lc[k]*ESsn1[kj]+Ls[k]*(1.-EScs1[kj]))/k;
	  ss	+=rc[k]*EScs1[kj]+rs[k]*ESsn1[kj];
	  EZdra	+=drc[k]*EScs1[kj]+drs[k]*ESsn1[kj];
	  EZdrt	+=rct[k]*ESsn1[kj]+rst[k]*EScs1[kj];
	}
	Gh[j]	=s;
	EZdra	=drc[0]+2.*EZdra;
	EZdrt	*=2.;
	ss	=2.*ss+rc[0];
	s	=b*EScs1[j];
	EZdra	=dgya/(ss*(EZdra*s-(dz+db*ESsn1[j])*EZdrt));
	ss	=aF/ss;
	qg0	+=1./((EZdrt*EZdrt+s*s)*EZdra*EZdra+ss*ss);
      }
      fac	*=ESNp/qg0;
      Gh[j]	=Gh[0];
      splP(Gh,Gh2);

      t	=-dgt*j0;
      for(j=0; j < nP4; j++){
	gh	=dgt*(j-j0);
	do{
	  ESSetSplP(t);
	  splRP(&s,&ss,Gh,Gh2);
	  s      +=t-gh;
	  t	-=s/(1.+ss);
	}while(fabs(s) > 1e-8);
	ss	=0.;
	EZdra	=0.;
	EZdrt	=0.;
	for(k=1; k < ESFp1; k++){
	  rf	=k*t;
	  cs	=cos(rf);
	  sn	=sin(rf);
	  ss	+=rc[k]*cs+rs[k]*sn;
	  EZdra	+=drc[k]*cs+drs[k]*sn;
	  EZdrt	+=rct[k]*sn+rst[k]*cs;
	}
	EZdra	=drc[0]+2.*EZdra;
	EZdrt	*=2.;

	ss	=2.*ss+rc[0];
	s	=b*cos(t);
	EZdra	=dgya/(ss*(EZdra*s-(dz+db*sin(t))*EZdrt));
	ss	=aF/ss;
	r[j]	=fac/((EZdrt*EZdrt+s*s)*EZdra*EZdra+ss*ss);
      }
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    }

    rf	=A[ia+1];
    dgya=(ESdgY[ESNa]+ESdgY1a[ESNa]*rf)*(ESsa[ESNa]+rf);
    aF	=ESaF[ESNa]+ESaF1a[ESNa]*rf;
    ss	=dgya/aF;
    ss	*=ss*ESaR0[0];
    Lc[0]	=ESLc[ESNa]+ESLc1[ESNa]*rf;
    s		=ESg22c[ESNa]+ESg22c1[ESNa]*rf;
    Lc[0]	=2./(Lc[0]+s*ss);
    ki	=ESNa;
    for(k=1; k < ES2Mp1; k++){
      ki	+=ESNa1;
      Lc[k]	=ESLc[ki]+ESLc1[ki]*rf;
      Ls[k]	=ESLs[ki]+ESLs1[ki]*rf;
      s		=ESg22c[ki]+ESg22c1[ki]*rf;
      Lc[k]	=(Lc[k]+s*ss)*Lc[0];
      s		=ESg22s[ki]+ESg22s1[ki]*rf;
      Ls[k]	=(Ls[k]+s*ss)*Lc[0];
    }
    s	=A[ia];
    a	=ESsa[ESNa]+s;
    dgya=(ESdgY[ESNa]+ESdgY1a[ESNa]*s)*a;
    fac	=(ESDPc[ESNa]+ESDPc1[ESNa]*s)*ESgY[ESNa]/dgya;
    b	=ESsb[ESNa]+ESsb1a[ESNa]*s;
    db	=ESsb1a[ESNa]+ESsb2a[ESNa]*s;
    dz	=rsT1a[ESNa]+rsT2a[ESNa]*s;
    rc[0]	=ESaR0[0]+rcT[ESNa]+rcT1a[ESNa]*s;
    drc[0]	=rcT1a[ESNa]+rcT2a[ESNa]*s;
    ki	=ESNa;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=rcT[ki]+rcT1a[ki]*s;
      drc[k]	=rcT1a[ki]+rcT2a[ki]*s;
      rct[k]	=-k*rc[k];
      rs[k]	=rsT[ki]+rsT1a[ki]*s;
      drs[k]	=rsT1a[ki]+rsT2a[ki]*s;
      rst[k]	=k*rs[k];
    }
    qg0	=0.;
    for(j=0; j < ESNp; j++){
      s	=0.;
      ss	=0.;
      EZdra	=0.;
      EZdrt	=0.;
      kj	=0;
      for(k=1; k < ES2Mp1; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	s	+=(Lc[k]*ESsn1[kj]+Ls[k]*(1.-EScs1[kj]))/k;
	ss	+=rc[k]*EScs1[kj]+rs[k]*ESsn1[kj];
	EZdra	+=drc[k]*EScs1[kj]+drs[k]*ESsn1[kj];
	EZdrt	+=rct[k]*ESsn1[kj]+rst[k]*EScs1[kj];
      }
      Gh[j]	=s;
      EZdra	=drc[0]+2.*EZdra;
      EZdrt	*=2.;
      ss	=2.*ss+rc[0];
      s		=b*EScs1[j];
      EZdra	=dgya/(ss*(EZdra*s-(dz+db*ESsn1[j])*EZdrt));
      ss	=aF/ss;
      qg0	+=1./((EZdrt*EZdrt+s*s)*EZdra*EZdra+ss*ss);
    }
    fac		*=ESNp/qg0;
    Gh[j]	=Gh[0];
    splP(Gh,Gh2);

    t	=-dgt*j0;
    for(j=0; j < nP4; j++){
      gh	=dgt*(j-j0);
      do{
	ESSetSplP(t);
	splRP(&s,&ss,Gh,Gh2);
	s      +=t-gh;
	t	-=s/(1.+ss);
      }while(fabs(s) > 1e-8);
      ss	=0.;
      EZdra	=0.;
      EZdrt	=0.;
      for(k=1; k < ESFp1; k++){
	rf	=k*t;
	cs	=cos(rf);
	sn	=sin(rf);
	ss	+=rc[k]*cs+rs[k]*sn;
	EZdra	+=drc[k]*cs+drs[k]*sn;
	EZdrt	+=rct[k]*sn+rst[k]*cs;
      }
      EZdra	=drc[0]+2.*EZdra;
      EZdrt	*=2.;
      ss	=2.*ss+rc[0];
      s		=b*cos(t);
      EZdra	=dgya/(ss*(EZdra*s-(dz+db*sin(t))*EZdrt));
      ss	=aF/ss;
      r[j]	=fac/((EZdrt*EZdrt+s*s)*EZdra*EZdra+ss*ss);
    }
    i++;
    fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
    for(j=0; j < nP4; j++){
      r[j]	=0.;
    }
    while(i < nA5){
      fwrite((void*)r,(size_t)8,(size_t)nP5,fp);
      i++;
    }
  }
  fclose(fp);
  Fl++;
  printf("\r");
  printf("Esc2dsk%4.2f has been written\n",ttime);
  return(0);
}

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "esc.h"

#define DEBUG00

#ifndef mzl_ESmain
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
#endif

extern char ShotNm[];

#ifndef mzl_1Dcore
extern int ESNa,ESNa1;
extern double ESa0,*ESsa,*ESpa,*ESgt;
#endif

#ifndef mzl_2Dcore
extern int ESNp,ESNp1,ESnAP;
extern int *ESnF,*ESmF,*ESkF,ESNf,ESNf1;
extern int ESMp,ESMp1,ES2Mp,ES2Mp1,ESnMp;
extern int ESFp,ESFp1,ESnAF;

extern double *ESemK,*ESemL,*ESemV,*ESemW,*ESemU1,*ESemU2;
extern double *ESemK2,*ESemL2,*ESemV2,*ESemW2,*ESemU12,*ESemU22;

extern double *ESr22,*ESr222,*ESr22c,*ESr22s,*ESr22c2,*ESr22s2;
extern double *ESg22c,*ESg22s,*ESg22c1,*ESg22s1,*ESg22c2,*ESg22s2;
extern double *ESg12c,*ESg12s,*ESg11c,*ESg11s;
extern double *ESg12c1,*ESg12s1,*ESg11c1,*ESg11s1;
extern double *ESg12c2,*ESg12s2,*ESg11c2,*ESg11s2;
extern double *ESLc,*ESLs,*ESVc,*ESVs,*ESDPc,*ESDPs;
extern double *ESLc1,*ESLs1,*ESVc1,*ESVs1,*ESDPc1,*ESDPs1;
extern double *ESLc2,*ESLs2,*ESVc2,*ESVs2,*ESDPc2,*ESDPs2;
extern double *ESjp,*ESjp1a,*ESjp2a,*ESjs,*ESjs1a,*ESjs2a,*ESsp,*ESsp1a
,*ESdgY,*ESdgY1a,*ESdgY2a;

extern double *EScs1,*ESsn1,*EZcs2,*EZsn2;
extern double *EZdra,*EZdrgt,*EZdza,*EZdzgt;
extern double ESgbext;
#endif

#ifndef mzl_3Dcore
extern int ESLt,ESNt,ESNt1,ESnAT,ESnAPT,ESLt0,ESNLt,ESNLt1;
extern double *ESsr,*ESsz,*rT1a,*zT1a,*rT1gt,*zT1gt;
extern double *ESaR0,*ESaZ0,*ESsb,*ESsb1a,*ESsb2a
,*rcT,*rsT,*rcT1a,*rsT1a,*rsT2a,*rcT2a;
#endif
extern double *ESFaux,*ESd2Faux;

#ifndef mzl_Esc2RZ
static double Rbox[4],Zbox[4],Rmx=0.,Rmn=0.,Zmx=0.,Zmn=0.;
static double gtrmx=0.,gtrmn=0.,Zrmx=0.,Zrmn=0.;
#endif

#ifdef DEBUG
int ESInspectMetricTensor(double *xc, double *yc, double *ys ,int m
			  ,int i0,int i1,int j0,int j1
			  ,int k0,int k1,int Fl,int Fd,int n)
{
  int i,j,k,K,N;
  int ki,jj,kk;
  double *p,*p1,*p2,*q,*q1,*q2,s,t,s1,ds,d2s;
  double fc[ESFp1],fs[ESFp1],fac[ESFp1],fas[ESFp1],ftc[ESFp1],fts[ESFp1];
  double r,ra,rt,za,zt,cs,sn;
  double g[ESNp1],gt[ESNp1],f[ESNa1],x[2],y[4];
  double W[ESNp1][n*(i1-i0)+1];
  double F[ESNp1][i1-i0+1];
  char *NmP,NmY[32];

  j1++;
  if(j1-j0 > ESNp1) j1=j0+ESNp1;
  N	=n*(i1-i0)+1;
  x[0]	=ESsa[i0];
  x[1]	=ESsa[i1];
  i1++;
  switch(Fl){
  case 0:
    NmP	="r";
    p	=rcT;
    p1	=rcT1a;
    p2	=rcT2a;
    q	=rsT;
    q1	=rsT1a;
    q2	=rsT2a;
    K	=ESFp1;
    for(i=i0; i < i1; i++){
      fc[0]	=rcT[i];
      ki	=i;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	fc[k]	=2.*rcT[ki];
	fs[k]	=2.*rsT[ki];
      }
      for(jj=j0; jj < j1; jj++){
	j	=jj;
	while(j < 0) j +=ESNp;
	ki	=0;
	r	=fc[0];
	for(k=1; k < K; k++){
	  ki	+=j;
	  if(ki >= ESNp) ki -=ESNp;
	  r	+=fc[k]*EScs1[ki]+fs[k]*ESsn1[ki];
	}
	F[jj-j0][i-i0]	=r;
      }
    }
    break;
  case 1:
    NmP	="z";
    p	=rsT;
    p1	=rsT1a;
    p2	=rsT2a;
    q	=ESsb;
    q1	=ESsb1a;
    q2	=ESsb2a;
    K	=1;
    for(i=i0; i < i1; i++){
      for(jj=j0; jj < j1; jj++){
	j	=jj;
	while(j < 0) j +=ESNp;
	F[jj-j0][i-i0]	=rsT[i]+ESsb[i]*ESsn1[j];;
      }
    }
    break;
  case 2:
    NmP	="g11";
    p	=ESg11c;
    p1	=ESg11c1;
    p2	=ESg11c2;
    q	=ESg11s;
    q1	=ESg11s1;
    q2	=ESg11s2;
    K	=ES2Mp1;
    for(i=i0; i < i1; i++){
      fc[0]	=ESaR0[0]+rcT[i];
      fac[0]	=rcT1a[i];
      ftc[0]	=0;
      ki	=i;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	fc[k]	=2.*rcT[ki];
	fs[k]	=2.*rsT[ki];
	fac[k]	=2.*rcT1a[ki];
	fas[k]	=2.*rsT1a[ki];
	if(i){
	  ftc[k]	=k*fs[k]/ESsa[i];
	  fts[k]	=-k*fc[k]/ESsa[i];
	}
	else{
	  ftc[k]	=k*fas[k];
	  fts[k]	=-k*fac[k];
	}
      }
      for(jj=j0; jj < j1; jj++){
	j	=jj;
	while(j < 0) j +=ESNp;
	r	=fc[0];
	ra	=fac[0];
	rt	=0.;
	ki	=0;
	for(k=1; k < K; k++){
	  ki	+=j;
	  if(ki >= ESNp) ki -=ESNp;
	  cs	=EScs1[ki];
	  sn	=ESsn1[ki];
	  r	+=fc[k]*cs+fs[k]*sn;
	  ra	+=fac[k]*cs+fas[k]*sn;
	  rt	+=ftc[k]*cs+fts[k]*sn;
	}
	za	=rsT1a[i]+ESsb1a[i]*ESsn1[j];
	zt	=i ? ESsb[i]*EScs1[j]/ESsa[i] : ESsb1a[i]*EScs1[j];
	F[jj-j0][i-i0]	=(ra*ra+za*za)/(r*(ra*zt-rt*za));
      }
    }
    break;
  case 3:
    NmP	="g12";
    p	=ESg12c;
    p1	=ESg12c1;
    p2	=ESg12c2;
    q	=ESg12s;
    q1	=ESg12s1;
    q2	=ESg12s2;
    K	=ES2Mp1;
    for(i=i0; i < i1; i++){
      fc[0]	=ESaR0[0]+rcT[i];
      fac[0]	=rcT1a[i];
      ftc[0]	=0;
      ki	=i;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	fc[k]	=2.*rcT[ki];
	fs[k]	=2.*rsT[ki];
	fac[k]	=2.*rcT1a[ki];
	fas[k]	=2.*rsT1a[ki];
	if(i){
	  ftc[k]	=k*fs[k]/ESsa[i];
	  fts[k]	=-k*fc[k]/ESsa[i];
	}
	else{
	  ftc[k]	=k*fas[k];
	  fts[k]	=-k*fac[k];
	}
      }
      for(jj=j0; jj < j1; jj++){
	j	=jj;
	while(j < 0) j +=ESNp;
	r	=fc[0];
	ra	=fac[0];
	rt	=0.;
	ki	=0;
	for(k=1; k < K; k++){
	  ki	+=j;
	  if(ki >= ESNp) ki -=ESNp;
	  cs	=EScs1[ki];
	  sn	=ESsn1[ki];
	  r	+=fc[k]*cs+fs[k]*sn;
	  ra	+=fac[k]*cs+fas[k]*sn;
	  rt	+=ftc[k]*cs+fts[k]*sn;
	}
	za	=rsT1a[i]+ESsb1a[i]*ESsn1[j];
	zt	=i ? ESsb[i]*EScs1[j]/ESsa[i] : ESsb1a[i]*EScs1[j];
	F[jj-j0][i-i0]	=(ra*rt+za*zt)/(r*(ra*zt-rt*za));
      }
    }
    break;
  case 4:
    NmP	="g22";
    p	=ESg22c;
    p1	=ESg22c1;
    p2	=ESg22c2;
    q	=ESg22s;
    q1	=ESg22s1;
    q2	=ESg22s2;
    K	=ES2Mp1;
    for(i=i0; i < i1; i++){
      fc[0]	=ESaR0[0]+rcT[i];
      fac[0]	=rcT1a[i];
      ftc[0]	=0.;
      ki	=i;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	fc[k]	=2.*rcT[ki];
	fs[k]	=2.*rsT[ki];
	fac[k]	=2.*rcT1a[ki];
	fas[k]	=2.*rsT1a[ki];
	if(i){
	  ftc[k]	=k*fs[k]/ESsa[i];
	  fts[k]	=-k*fc[k]/ESsa[i];
	}
	else{
	  ftc[k]	=k*fas[k];
	  fts[k]	=-k*fac[k];
	}
      }
      for(jj=j0; jj < j1; jj++){
	j	=jj;
	while(j < 0) j +=ESNp;
	r	=fc[0];
	ra	=fac[0];
	rt	=0.;
	ki	=0;
	for(k=1; k < K; k++){
	  ki	+=j;
	  if(ki >= ESNp) ki -=ESNp;
	  cs	=EScs1[ki];
	  sn	=ESsn1[ki];
	  r	+=fc[k]*cs+fs[k]*sn;
	  ra	+=fac[k]*cs+fas[k]*sn;
	  rt	+=ftc[k]*cs+fts[k]*sn;
	}
	za	=rsT1a[i]+ESsb1a[i]*ESsn1[j];
	zt	=i ? ESsb[i]*EScs1[j]/ESsa[i] : ESsb1a[i]*EScs1[j];
	F[jj-j0][i-i0]	=(rt*rt+zt*zt)/(r*(ra*zt-rt*za));
      }
    }
    break;
  case 5:
    NmP	="V";
    p	=ESVc;
    p1	=ESVc1;
    p2	=ESVc2;
    q	=ESVs;
    q1	=ESVs1;
    q2	=ESVs2;
    K	=ES2Mp1;
    for(i=i0; i < i1; i++){
      fc[0]	=ESaR0[0]+rcT[i];
      fac[0]	=rcT1a[i];
      ftc[0]	=0;
      ki	=i;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	fc[k]	=2.*rcT[ki];
	fs[k]	=2.*rsT[ki];
	fac[k]	=2.*rcT1a[ki];
	fas[k]	=2.*rsT1a[ki];
	if(i){
	  ftc[k]	=k*fs[k]/ESsa[i];
	  fts[k]	=-k*fc[k]/ESsa[i];
	}
	else{
	  ftc[k]	=k*fas[k];
	  fts[k]	=-k*fac[k];
	}
      }
      for(jj=j0; jj < j1; jj++){
	j	=jj;
	while(j < 0) j +=ESNp;
	r	=fc[0];
	ra	=fac[0];
	rt	=0.;
	ki	=0;
	for(k=1; k < K; k++){
	  ki	+=j;
	  if(ki >= ESNp) ki -=ESNp;
	  cs	=EScs1[ki];
	  sn	=ESsn1[ki];
	  r	+=fc[k]*cs+fs[k]*sn;
	  ra	+=fac[k]*cs+fas[k]*sn;
	  rt	+=ftc[k]*cs+fts[k]*sn;
	}
	za	=rsT1a[i]+ESsb1a[i]*ESsn1[j];
	zt	=i ? ESsb[i]*EScs1[j]/ESsa[i] : ESsb1a[i]*EScs1[j];
	F[jj-j0][i-i0]	=(ra*zt-rt*za)*(r/ESaR0[0]-ESaR0[0]/r);
      }
    }
    break;
  case 6:
    NmP	="L";
    p	=ESLc;
    p1	=ESLc1;
    p2	=ESLc2;
    q	=ESLs;
    q1	=ESLs1;
    q2	=ESLs2;
    K	=ES2Mp1;
    for(i=i0; i < i1; i++){
      fc[0]	=ESaR0[0]+rcT[i];
      fac[0]	=rcT1a[i];
      ftc[0]	=0;
      ki	=i;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	fc[k]	=2.*rcT[ki];
	fs[k]	=2.*rsT[ki];
	fac[k]	=2.*rcT1a[ki];
	fas[k]	=2.*rsT1a[ki];
	if(i){
	  ftc[k]	=k*fs[k]/ESsa[i];
	  fts[k]	=-k*fc[k]/ESsa[i];
	}
	else{
	  ftc[k]	=k*fas[k];
	  fts[k]	=-k*fac[k];
	}
      }
      for(jj=j0; jj < j1; jj++){
	j	=jj;
	while(j < 0) j +=ESNp;
	r	=fc[0];
	ra	=fac[0];
	rt	=0.;
	ki	=0;
	for(k=1; k < K; k++){
	  ki	+=j;
	  if(ki >= ESNp) ki -=ESNp;
	  cs	=EScs1[ki];
	  sn	=ESsn1[ki];
	  r	+=fc[k]*cs+fs[k]*sn;
	  ra	+=fac[k]*cs+fas[k]*sn;
	  rt	+=ftc[k]*cs+fts[k]*sn;
	}
	za	=rsT1a[i]+ESsb1a[i]*ESsn1[j];
	zt	=i ? ESsb[i]*EScs1[j]/ESsa[i] : ESsb1a[i]*EScs1[j];
	F[jj-j0][i-i0]	=(ra*zt-rt*za)*ESaR0[0]/r;
      }
    }
    break;
  case 7:
    NmP	="D";
    p	=ESDPc;
    p1	=ESDPc1;
    p2	=ESDPc2;
    q	=ESDPs;
    q1	=ESDPs1;
    q2	=ESDPs2;
    K	=ES2Mp1;
    for(i=i0; i < i1; i++){
      fc[0]	=ESaR0[0]+rcT[i];
      fac[0]	=rcT1a[i];
      ftc[0]	=0;
      ki	=i;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	fc[k]	=2.*rcT[ki];
	fs[k]	=2.*rsT[ki];
	fac[k]	=2.*rcT1a[ki];
	fas[k]	=2.*rsT1a[ki];
	if(i){
	  ftc[k]	=k*fs[k]/ESsa[i];
	  fts[k]	=-k*fc[k]/ESsa[i];
	}
	else{
	  ftc[k]	=k*fas[k];
	  fts[k]	=-k*fac[k];
	}
      }
      for(jj=j0; jj < j1; jj++){
	j	=jj;
	while(j < 0) j +=ESNp;
	r	=fc[0];
	ra	=fac[0];
	rt	=0.;
	ki	=0;
	for(k=1; k < K; k++){
	  ki	+=j;
	  if(ki >= ESNp) ki -=ESNp;
	  cs	=EScs1[ki];
	  sn	=ESsn1[ki];
	  r	+=fc[k]*cs+fs[k]*sn;
	  ra	+=fac[k]*cs+fas[k]*sn;
	  rt	+=ftc[k]*cs+fts[k]*sn;
	}
	za	=rsT1a[i]+ESsb1a[i]*ESsn1[j];
	zt	=i ? ESsb[i]*EScs1[j]/ESsa[i] : ESsb1a[i]*EScs1[j];
	F[jj-j0][i-i0]	=ra*zt-rt*za;
      }
    }
    break;
  default :
    NmP	="r";
    p	=rcT;
    p1	=rcT1a;
    p2	=rcT2a;
    q	=rsT;
    q1	=rsT1a;
    q2	=rsT2a;
    K	=ESFp1;
    for(i=i0; i < i1; i++){
      fc[0]	=rcT[i];
      ki	=i;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	fc[k]	=2.*rcT[ki];
	fs[k]	=2.*rsT[ki];
      }
      for(jj=j0; jj < j1; jj++){
	j	=jj;
	while(j < 0) j +=ESNp;
	ki	=0;
	ds	=fc[0];
	for(k=1; k < K; k++){
	  ki	+=j;
	  if(ki >= ESNp) ki -=ESNp;
	  ds	+=fc[k]*EScs1[ki]+fs[k]*ESsn1[ki];
	}
	F[jj-j0][i-i0]	=ds;
      }
    }
    break;
  }
    
  y[0]	=0.;
  y[1]	=0.;
  y[2]	=0.;
  y[3]	=0.;
  s	=ESsa[1]/n;
  t	=x[0];
  kk	=0;
  while(kk < N){
    ESSetSplA(t);
    xc[kk]	=t;
    for(k=k0; k < k1; k++){
      ki	=ESNa1*k;
      splRA2(&s1,&ds,&d2s,p+ki,p2+ki);
      switch(Fd){
      case 0:
	ds	=s1;
	break;
      case 2:
	ds	=d2s;
	break;
      }
      if(y[0] > ds) y[0]=ds;
      if(y[1] < ds) y[1]=ds;
      W[k][kk]	=ds;
      if(k == m && yc != NULL) yc[kk]	=ds;
      if(k){
	splRA2(&s1,&ds,&d2s,q+ki,q2+ki);
	switch(Fd){
	case 0:
	  ds	=s1;
	  break;
	case 2:
	  ds	=d2s;
	  break;
	}
	if(y[2] > ds) y[2]=ds;
	if(y[3] < ds) y[3]=ds;
	W[k+33][kk]	=ds;
	if(k == m && ys != NULL) ys[kk]	=ds;
      }
    }
    t	+=s;
    kk++;
  }

  for(i=i0; i < i1; i++){
    for(k=k0; k < k1; k++){
      ki	=ESNa1*k;
      switch(Fd){
      case 0:
	ds	=p[ki+i];
	break;
      case 1:
	ds	=p1[ki+i];
	break;
      case 2:
	ds	=p2[ki+i];
	break;
      }
      if(y[0] > ds) y[0]=ds;
      if(y[1] < ds) y[1]=ds;
      if(k){
	switch(Fd){
	case 0:
	  ds	=q[ki+i];
	  break;
	case 1:
	  ds	=q1[ki+i];
	  break;
	case 2:
	  ds	=q2[ki+i];
	  break;
	}
	ds	=q[ki+i];
	if(y[2] > ds) y[2]=ds;
	if(y[3] < ds) y[3]=ds;
      }
    }
  }
  sprintf(NmY,"%s cos [%d <= k < %d]",NmP,k0,k1);
  SetPlotName("a",NmY,NmP);
  Scale2D(2,2,x,y,2,6,2);
  for(k=k0; k < k1; k++){
    ki	=ESNa1*k+i0;
    switch(Fd){
    case 0:
      Plot2d(2,ESsa+i0,p+ki,i1-i0,6,0,0,0);
      break;
    case 1:
      Plot2d(2,ESsa+i0,p1+ki,i1-i0,6,0,0,0);
      break;
    case 2:
      Plot2d(2,ESsa+i0,p2+ki,i1-i0,6,0,0,0);
      break;
    }
    if(k%2)	Plot2d(2,xc,W[k],kk,6,0,14,0);
    else	Plot2d(2,xc,W[k],kk,6,0,4,0);
  }
  sprintf(NmY,"%s sin [%d <= k < %d]",NmP,k0 ? k0 : 1,k1);
  SetPlotName("a",NmY,NmP);
  Scale2D(3,2,x,y+2,2,6,2);
  for(k=k0 ? k0 : 1; k < k1; k++){
    switch(Fd){
    case 0:
      Plot2d(3,ESsa+i0,q+ESNa1*k+i0,i1-i0,6,0,0,0);
      break;
    case 1:
      Plot2d(3,ESsa+i0,q1+ESNa1*k+i0,i1-i0,6,0,0,0);
      break;
    case 2:
      Plot2d(3,ESsa+i0,q2+ESNa1*k+i0,i1-i0,6,0,0,0);
      break;
    }
    if(k%2)Plot2d(3,xc,W[k+33],kk,6,0,14,0);
    else Plot2d(3,xc,W[k+33],kk,6,0,4,0);
  }

  y[0]	=0.;
  y[1]	=0.;
  kk	=0;
  t	=x[0];
  while(kk < N){
    ESSetSplA(t);
    splRA2(fc,&ds,&d2s,p,p2);
    if(Fl == 1) splRA2(fs,&ds,&d2s,q,q2);
    ki	=0;
    for(k=1; k < K; k++){
      ki	+=ESNa1;
      splRA2(fc+k,&ds,&d2s,p+ki,p2+ki);
      splRA2(fs+k,&ds,&d2s,q+ki,q2+ki);
      fc[k]	*=2.;
      fs[k]	*=2.;
    }
    for(jj=j0; jj < j1; jj++){
      j	=jj;
      while(j < 0) j +=ESNp;
      while(j > ESNp) j -=ESNp;
      ki	=0;
      ds	=fc[0];
      for(k=1; k < K; k++){
	ki	+=j;
	if(ki >= ESNp) ki -=ESNp;
	ds	+=fc[k]*EScs1[ki]+fs[k]*ESsn1[ki];
      }
      if(Fl == 1) ds	+=fs[0]*ESsn1[j];
      W[jj-j0][kk]	=ds;
      if(y[0] > ds) y[0]=ds;
      if(y[1] < ds) y[1]=ds;
    }
    t	+=s;
    kk++;
  } 

  sprintf(NmY,"%s [%d <= j < %d]",NmP,j0,j1);
  SetPlotName("a",NmY,NmP);
  Scale2D(5,2,x,y,2,6,2);
  for(jj=j0; jj < j1; jj++){
    gt[jj-j0]	=(double)jj/ESNp;
    Plot2d(5,ESsa+i0,F[jj-j0],i1-i0,6,0,0,0);
    Plot2d(5,xc,W[jj-j0],kk,6,0,4,0);
  }
  x[0]	=gt[0];
  x[1]	=gt[j1-j0-1];
  sprintf(NmY,"%s [%d <= i < %d]",NmP,i0,i1);
  SetPlotName("gt/gp",NmY,NmP);
  Scale2D(6,2,x,y,2,6,2);
  for(i=i0; i < i1; i++){
    j	=0;
    for(jj=j0; jj < j1; jj++){
      g[j]	=F[j][i-i0];
      j++;
    }
    Plot2d(6,gt,g,j,6,0,0,0);
  }
  for(k=0; k < kk; k++){
    j	=0;
    for(jj=j0; jj < j1; jj++){
      g[j]	=W[j][k];
      j++;
    }
    Plot2d(6,gt,g,j,6,0,4,0);
  }
  return(0);
}
#endif

int ES3DMetrTnsrCSHole()
{
  extern double ESa0;
  int i,j,ji,n1,im,im0,m,m1,mm;
  int k,ki,kj;
  double *rc,*rca,*rs,*rsa,*df1;

  double r0,rr0,b;
  double s,ba,rA;
  double raa,ragt,zaa,zagt,dg0,dg1;
  double *D,*L,*V,*G22,*G12,*G11;

  D	=EZdra+ESNp1;
  L	=D+ESNp1;
  rc	=L+ESNp1;

  V	=EZdza+ESNp1;
  G22	=V+ESNp1;
  rca	=G22+ESNp1;

  G12	=EZdrgt+ESNp1;
  rs	=G12+ESNp1;

  G11	=EZdzgt+ESNp1;
  rsa	=G11+ESNp1;

  r0	=ESaR0[0];
  rr0	=1./r0;
  im0	=0;

  i	=0;
  ji	=0;
  rA	=1./ESsa[i];
  b	=ESsb[i];
  ba	=ESsb1a[i];
  ki	=i;
  rca[0]	=rcT1a[ki];
  rsa[0]	=rsT1a[ki];
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    rca[k]	=2.*rcT1a[ki];
    rsa[k]	=2.*rsT1a[ki];
    rc[k]	=2.*k*rcT[ki];
    rs[k]	=2.*k*rsT[ki];
  }
  for(j=0; j < ESNp; j++){
    EZdra[j]	=rca[0];
    EZdrgt[j]	=0.;
    kj	=0;
    for(k=1; k < ESFp1; k++){
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      EZdra[j]	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
      EZdrgt[j]	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
    }
    EZdzgt[j]	=b*EScs1[j];
    EZdza[j]	=rsa[0]+ba*ESsn1[j];
    D[j]	=EZdra[j]*EZdzgt[j]-EZdrgt[j]*EZdza[j];
    if(D[j] <= 0.){
      printf("D[i=%2d][j=%2d]=%10.3e <= 0%c\n",i,j,D[j],7);
      return(1);
    }
    L[j]	=D[j]*r0/ESsr[ji]*rA;
    V[j]	=D[j]*ESsr[ji]*rr0*rA-L[j];
    dg0		=1./(D[j]*ESsr[ji]);
    D[j]	*=rA;
    G11[j]	=ESsa[i]*(EZdra[j]*EZdra[j]+EZdza[j]*EZdza[j])*dg0;
    G12[j]	=(EZdra[j]*EZdrgt[j]+EZdza[j]*EZdzgt[j])*dg0;
    G22[j]	=(EZdrgt[j]*EZdrgt[j]+EZdzgt[j]*EZdzgt[j])*dg0*rA;
    rT1a[ji]	=EZdra[j];
    rT1gt[ji]	=EZdrgt[j];
    zT1a[ji]	=EZdza[j];
    zT1gt[ji]	=EZdzgt[j];
    ji++;
  }
  D[j]	=D[0];
  L[j]	=L[0];
  V[j]	=V[0];
  G11[j]	=G11[0];
  G12[j]	=G12[0];
  G22[j]	=G22[0];
  rT1a[ji]	=EZdra[0];
  rT1gt[ji]	=EZdrgt[0];
  zT1a[ji]	=EZdza[0];
  zT1gt[ji]	=EZdzgt[0];
  ji++;
  im		=0;
  ESgP2gF(ESDPc,ESDPs,D,ES2Mp);
  ESgP2gF(ESLc,ESLs,L,ES2Mp);
  ESgP2gF(ESVc,ESVs,V,ES2Mp);
  ESgP2gF(ESg22c,ESg22s,G22,ES2Mp);
  ESgP2gF(ESg12c,ESg12s,G12,ES2Mp);
  ESgP2gF(ESg11c,ESg11s,G11,ES2Mp);
  
  ji		-=ESNp1;
  ki		=i;
  rca[0]	=rcT2a[ki];
  rsa[0]	=rsT2a[ki];
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    rca[k]	=2.*rcT2a[ki];
    rsa[k]	=2.*rsT2a[ki];
    b		=2.*k;
    rc[k]	=b*rcT1a[ki];
    rs[k]	=b*rsT1a[ki];
  }
  b		=ESsb1a[i];
  ba		=ESsb2a[i];
  for(j=0; j < ESNp; j++){
    raa		=rca[0];
    ragt	=0.;
    kj		=0;
    for(k=1; k < ESFp1; k++){
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      raa	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
      ragt	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
    }
    zaa		=rsa[0]+ba*ESsn1[j];
    zagt	=b*EScs1[j];
    dg0		=-EZdra[j]/ESsr[ji];
    dg1		=raa*EZdzgt[j]+EZdra[j]*zagt-ragt*EZdza[j]-EZdrgt[j]*zaa;
    L[j]	=(dg0-rA)*L[j]+dg1*r0/ESsr[ji];
    V[j]	=(D[j]*EZdra[j]+dg1*ESsr[ji]-D[j]*ESsr[ji]*rA)*rA*rr0-L[j];
    dg0		-=dg1/D[j];

    G12[j]	*=dg0;
    G22[j]	*=(dg0-rA);
    dg0		+=rA;
    G11[j]	*=dg0;
    dg0		=1./(D[j]*ESsr[ji]);
    G11[j]	+=2.*(EZdra[j]*raa+EZdza[j]*zaa)*dg0*ESsa[i];
    G12[j]	+=(raa*EZdrgt[j]+EZdra[j]*ragt+zaa*EZdzgt[j]+EZdza[j]*zagt)*dg0;
    G22[j]	+=2.*(EZdrgt[j]*ragt+EZdzgt[j]*zagt)*dg0*rA;
    D[j]	=(dg1-D[j]*rA)*rA;
    ji++;
  }
  D[j]	=D[0];
  L[j]	=L[0];
  V[j]	=V[0];
  G11[j]	=G11[0];
  G12[j]	=G12[0];
  G22[j]	=G22[0];
  ji++;
  im	=0;
  ESgP2gF(ESDPc2,ESDPs2,D,ES2Mp);
  ESgP2gF(ESLc2,ESLs2,L,ES2Mp);
  ESgP2gF(ESVc2,ESVs2,V,ES2Mp);
  ESgP2gF(ESg22c2,ESg22s2,G22,ES2Mp);
  ESgP2gF(ESg12c2,ESg12s2,G12,ES2Mp);
  ESgP2gF(ESg11c2,ESg11s2,G11,ES2Mp);

  for(i=1; i < ESNa1; i++){
    rA	=1./ESsa[i];
    b	=ESsb[i];
    ba	=ESsb1a[i];
    ki	=i;
    rca[0]	=rcT1a[ki];
    rsa[0]	=rsT1a[ki];
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rca[k]	=2.*rcT1a[ki];
      rsa[k]	=2.*rsT1a[ki];
      rc[k]	=2.*k*rcT[ki];
      rs[k]	=2.*k*rsT[ki];
    }
    for(j=0; j < ESNp; j++){
      EZdra[j]	=rca[0];
      EZdrgt[j]	=0.;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj		+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	EZdra[j]	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
	EZdrgt[j]	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
      }
      EZdzgt[j]	=b*EScs1[j];
      EZdza[j]	=rsa[0]+ba*ESsn1[j];
      D[j]	=EZdra[j]*EZdzgt[j]-EZdrgt[j]*EZdza[j];
      if(D[j] <= 0.){
	printf("D[i=%2d][j=%2d]=%10.3e <= 0%c\n",i,j,D[j],7);
	return(1);
      }
      L[j]	=D[j]*r0/ESsr[ji]*rA;
      V[j]	=D[j]*ESsr[ji]*rr0*rA-L[j];
      dg0	=1./(D[j]*ESsr[ji]);
      D[j]	*=rA;
      G11[j]	=ESsa[i]*(EZdra[j]*EZdra[j]+EZdza[j]*EZdza[j])*dg0;
      G12[j]	=(EZdra[j]*EZdrgt[j]+EZdza[j]*EZdzgt[j])*dg0;
      G22[j]	=(EZdrgt[j]*EZdrgt[j]+EZdzgt[j]*EZdzgt[j])*dg0*rA;
      rT1a[ji]	=EZdra[j];
      rT1gt[ji]	=EZdrgt[j];
      zT1a[ji]	=EZdza[j];
      zT1gt[ji]	=EZdzgt[j];
      ji++;
    }
    D[j]	=D[0];
    L[j]	=L[0];
    V[j]	=V[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    rT1a[ji]	=EZdra[0];
    rT1gt[ji]	=EZdrgt[0];
    zT1a[ji]	=EZdza[0];
    zT1gt[ji]	=EZdzgt[0];
    ji++;
    ESgP2gF(ESDPc+i,ESDPs+i,D,ES2Mp);
    ESgP2gF(ESLc+i,ESLs+i,L,ES2Mp);
    ESgP2gF(ESVc+i,ESVs+i,V,ES2Mp);
    ESgP2gF(ESg22c+i,ESg22s+i,G22,ES2Mp);
    ESgP2gF(ESg12c+i,ESg12s+i,G12,ES2Mp);
    ESgP2gF(ESg11c+i,ESg11s+i,G11,ES2Mp);
  }
  i		=ESNa;
  ji		-=ESNp1;
  ki		=i;
  rca[0]	=rcT2a[ki];
  rsa[0]	=rsT2a[ki];
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    rca[k]	=2.*rcT2a[ki];
    rsa[k]	=2.*rsT2a[ki];
    b		=2.*k;
    rc[k]	=b*rcT1a[ki];
    rs[k]	=b*rsT1a[ki];
  }
  b		=ESsb1a[i];
  ba		=ESsb2a[i];
  for(j=0; j < ESNp; j++){
    raa	=rca[0];
    ragt=0.;
    kj	=0;
    for(k=1; k < ESFp1; k++){
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      raa	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
      ragt	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
    }
    zaa	=rsa[0]+ba*ESsn1[j];
    zagt=b*EScs1[j];
    dg0	=-EZdra[j]/ESsr[ji];
    dg1	=raa*EZdzgt[j]+EZdra[j]*zagt-ragt*EZdza[j]-EZdrgt[j]*zaa;
    L[j]=(dg0-1.)*L[j]+dg1*r0/ESsr[ji];
    V[j]=(D[j]*EZdra[j]+dg1*ESsr[ji]-D[j]*ESsr[ji])*rr0-L[j];
    dg0	-=dg1/D[j];
    G12[j]	*=dg0;
    G22[j]	*=(dg0-1.);
    dg0	+=1.;
    G11[j]	*=dg0;
    dg0	=1./(D[j]*ESsr[ji]);
    G11[j]	+=2.*(EZdra[j]*raa+EZdza[j]*zaa)*dg0*ESsa[i];
    G12[j]	+=(raa*EZdrgt[j]+EZdra[j]*ragt+zaa*EZdzgt[j]+EZdza[j]*zagt)*dg0;
    G22[j]	+=2.*(EZdrgt[j]*ragt+EZdzgt[j]*zagt)*dg0;
    D[j]	=dg1-D[j];
    ji++;
  }
  D[j]	=D[0];
  L[j]	=L[0];
  V[j]	=V[0];
  G11[j]=G11[0];
  G12[j]=G12[0];
  G22[j]=G22[0];
  ji++;
  ESgP2gF(ESDPc2+i,ESDPs2+i,D,ES2Mp);
  ESgP2gF(ESLc2+i,ESLs2+i,L,ES2Mp);
  ESgP2gF(ESVc2+i,ESVs2+i,V,ES2Mp);
  ESgP2gF(ESg22c2+i,ESg22s2+i,G22,ES2Mp);
  ESgP2gF(ESg12c2+i,ESg12s2+i,G12,ES2Mp);
  ESgP2gF(ESg11c2+i,ESg11s2+i,G11,ES2Mp);
  i	=0;
  j	=ESNa;
  dg0	=ESDPc2[i];
  dg1	=ESDPc2[j];
  splAA(ESDPc+i,ESDPc1+i,ESDPc2+i,&dg0,&dg1);
  dg0	=ESLc2[i];
  dg1	=ESLc2[j];
  splAA(ESLc+i,ESLc1+i,ESLc2+i,&dg0,&dg1);
  dg0	=ESVc2[i];
  dg1	=ESVc2[j];
  splAA(ESVc+i,ESVc1+i,ESVc2+i,&dg0,&dg1);
  dg0	=ESg22c2[i];
  dg1	=ESg22c2[j];
  splAA(ESg22c+i,ESg22c1+i,ESg22c2+i,&dg0,&dg1);
  dg0	=ESg12c2[i];
  dg1	=ESg12c2[j];
  splAA(ESg12c+i,ESg12c1+i,ESg12c2+i,&dg0,&dg1);
  dg0	=ESg11c2[i];
  dg1	=ESg11c2[j];
  splAA(ESg11c+i,ESg11c1+i,ESg11c2+i,&dg0,&dg1);
  for(m=1; m < ES2Mp1; m++){
    i 	+=ESNa1;
    j	=i+ESNa;
    dg0	=ESDPc2[i];
    dg1	=ESDPc2[j];
    splAA(ESDPc+i,ESDPc1+i,ESDPc2+i,&dg0,&dg1);
    dg0	=ESDPs2[i];
    dg1	=ESDPs2[j];
    splAA(ESDPs+i,ESDPs1+i,ESDPs2+i,&dg0,&dg1);
    dg0	=ESLc2[i];
    dg1	=ESLc2[j];
    splAA(ESLc+i,ESLc1+i,ESLc2+i,&dg0,&dg1);
    dg0	=ESLs2[i];
    dg1	=ESLs2[j];
    splAA(ESLs+i,ESLs1+i,ESLs2+i,&dg0,&dg1);
    dg0	=ESVc2[i];
    dg1	=ESVc2[j];
    splAA(ESVc+i,ESVc1+i,ESVc2+i,&dg0,&dg1);
    dg0	=ESVs2[i];
    dg1	=ESVs2[j];
    splAA(ESVs+i,ESVs1+i,ESVs2+i,&dg0,&dg1);
    dg0	=ESg22c2[i];
    dg1	=ESg22c2[j];
    splAA(ESg22c+i,ESg22c1+i,ESg22c2+i,&dg0,&dg1);
    dg0	=ESg22s2[i];
    dg1	=ESg22s2[j];
    splAA(ESg22s+i,ESg22s1+i,ESg22s2+i,&dg0,&dg1);
    dg0	=ESg12c2[i];
    dg1	=ESg12c2[j];
    splAA(ESg12c+i,ESg12c1+i,ESg12c2+i,&dg0,&dg1);
    dg0	=ESg12s2[i];
    dg1	=ESg12s2[j];
    splAA(ESg12s+i,ESg12s1+i,ESg12s2+i,&dg0,&dg1);
    dg0	=ESg11c2[i];
    dg1	=ESg11c2[j];
    splAA(ESg11c+i,ESg11c1+i,ESg11c2+i,&dg0,&dg1);
    dg0	=ESg11s2[i];
    dg1	=ESg11s2[j];
    splAA(ESg11s+i,ESg11s1+i,ESg11s2+i,&dg0,&dg1);
  }
  return(0);
}

int ES3DMetrTnsrCSHoleV()
{
  int i,in,j,ji,n,n1,jni,im,im0,m,m1,mm;
  int K,k,ki,kj;
  double *rc,*rca,*rs,*rsa,*df1;

  double r0,rr0,b;
  double s,ba,rA;
  double raa,ragt,zaa,zagt,dg0,dg1;
  double *D,*L,*V,*G22,*G12,*G11;

  D	=EZdra+ESNp1;
  L	=D+ESNp1;
  rc	=L+ESNp1;

  V	=EZdza+ESNp1;
  G22	=V+ESNp1;
  rca	=G22+ESNp1;

  G12	=EZdrgt+ESNp1;
  rs	=G12+ESNp1;

  G11	=EZdzgt+ESNp1;
  rsa	=G11+ESNp1;

  for(n=0; n < ESNLt; n++){
    K		=ESnAF*n;
    r0		=ESaR0[n];
    rr0		=1./r0;
    jni		=ESnAP*n;
    in		=ESNa1*n;
    im0		=ESnMp*n;
    jni		=ESnAP*n;
    i		=0;
    b		=ESsb1a[in];
    ki		=K+ESNa1;
    rc[1]	=2.*rcT1a[ki];
    rs[1]	=2.*rsT1a[ki];
    ba		=b*rc[1];
    for(j=0; j < ESNp; j++){
      EZdrgt[j]	=-rc[1]*ESsn1[j]+rs[1]*EScs1[j];
      EZdzgt[j]	=b*EScs1[j];
      EZdra[j]	=rc[1]*EScs1[j]+rs[1]*ESsn1[j];
      EZdza[j]	=b*ESsn1[j];
      L[j]	=ba;
      dg0	=rr0/ba;
      G11[j]	=(EZdra[j]*EZdra[j]+EZdza[j]*EZdza[j])*dg0;
      G12[j]	=(EZdra[j]*EZdrgt[j]+EZdza[j]*EZdzgt[j])*dg0;
      G22[j]	=(EZdrgt[j]*EZdrgt[j]+EZdzgt[j]*EZdzgt[j])*dg0;
      rT1a[jni]	=EZdra[j];
      rT1gt[jni]=EZdrgt[j];
      zT1a[jni]	=EZdza[j];
      zT1gt[jni]=EZdzgt[j];
      jni++;
    }
    L[j]	=L[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    rT1a[jni]	=EZdra[0];
    rT1gt[jni]	=EZdrgt[0];
    zT1a[jni]	=EZdza[0];
    zT1gt[jni]	=EZdzgt[0];
    jni++;
    im		=im0;
    ESgP2gF(ESg22c+im,ESg22s+im,G22,ES2Mp);
    m1		=im;
    ESg11c[im]	=ESg22c[im];
    ESg11s[im]	=0.;
    ESg12c[im]	=0.;
    ESg12s[im]	=0.;
    ESDPc[im]	=ba;
    ESDPs[im]	=0.;
    ESLc[im]	=ba;
    ESLs[im]	=0.;
    ESVc2[im]	=0.;
    ESVs2[im]	=0.;
    m1	+=ESNa1;
    ESg11c[m1]	=0.;
    ESg11s[m1]	=0.;
    ESg12c[m1]	=0.;
    ESg12s[m1]	=0.;
    ESg22c[m1]	=0.;
    ESg22s[m1]	=0.;
    ESDPc[m1]	=0.;
    ESDPs[m1]	=0.;
    ESLc[m1]	=0.;
    ESLs[m1]	=0.;
    ESVc2[m1]	=ba*rc[1]*rr0;
    ESVs2[m1]	=ba*rs[1]*rr0;
    m1	+=ESNa1;
    ESg12c[m1]	=-ESg22s[m1];
    ESg12s[m1]	=ESg22c[m1];
    ESg11c[m1]	=-ESg22c[m1];
    ESg11s[m1]	=-ESg22s[m1];
    ESDPc[m1]	=0.;
    ESDPs[m1]	=0.;
    ESLc[m1]	=0.;
    ESLs[m1]	=0.;
    ESVc2[m1]	=0.;
    ESVs2[m1]	=0.;
    for(m=3; m < ES2Mp1; m++){
      m1	+=ESNa1;
      ESg11c[m1]	=0.;
      ESg11s[m1]	=0.;
      ESg12c[m1]	=0.;
      ESg12s[m1]	=0.;
      ESDPc[m1]		=0.;
      ESDPs[m1]		=0.;
      ESLc[m1]		=0.;
      ESLs[m1]		=0.;
      ESVc2[m1]		=0.;
      ESVs2[m1]		=0.;
    }
    ki		=K;
    rca[0]	=rcT2a[ki];
    rsa[0]	=rsT2a[ki];
    ki		+=ESNa1;
    k		=1;
    rca[1]	=2.*rcT1a[ki];
    rsa[1]	=2.*rsT1a[ki];
    ki		+=ESNa1;
    k		=2;
    rca[2]	=2.*rcT2a[ki];
    rsa[2]	=2.*rsT2a[ki];
    for(j=0; j < ESNp; j++){
      raa	=rca[0]+rca[2]*EZcs2[j]+rsa[2]*EZsn2[j];
      zaa	=rsa[0];
      ragt	=-rca[2]*EZsn2[j]+rsa[2]*EZcs2[j];
      dg1	=(b*(rca[0]+rca[2])-rsa[1]*rsa[0])*EScs1[j]
	+(b*rsa[2]+rca[1]*rsa[0])*ESsn1[j];
      dg0	=-EZdra[j]*rr0;
      D[j]	=dg1;
      L[j]	=L[j]*dg0+dg1;
      dg0	-=dg1/ba;
      G11[j]	*=dg0;
      G12[j]	*=dg0;
      G22[j]	*=dg0;
      dg0	=rr0/ba;
      G11[j]	+=2.*(EZdra[j]*raa+EZdza[j]*zaa)*dg0;
      G12[j]	+=(raa*EZdrgt[j]+EZdra[j]*ragt+zaa*EZdzgt[j])*dg0;
      G22[j]	+=2.*ragt*EZdrgt[j]*dg0;
    }
    D[j]	=D[0];
    L[j]	=L[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    im		=im0+i;
    ESgP2gF(ESDPc2+im,ESDPs2+im,D,ES2Mp);
    ESgP2gF(ESLc2+im,ESLs2+im,L,ES2Mp);
    m1		=im;
    for(m=0; m < ES2Mp1; m++){
      ESVc[m1]	=0.;
      ESVs[m1]	=0.;
      m1	+=ESNa1;
    }
    ESgP2gF(ESg11c2+im,ESg11s2+im,G11,ES2Mp);
    ESgP2gF(ESg12c2+im,ESg12s2+im,G12,ES2Mp);
    ESgP2gF(ESg22c2+im,ESg22s2+im,G22,ES2Mp);
#ifdef H
    ki		=K;
    rca[0]	=rcT2a[ki];
    rsa[0]	=rsT2a[ki];
    ki	+=ESNa1;
    rc[1]	=2.*rcT1a[ki];
    rs[1]	=2.*rsT1a[ki];
    rca[1]	=2.*rcT2a[ki];
    rsa[1]	=2.*rsT2a[ki];
    for(k=2; k < ESFp1; k++){
      ki	+=ESNa1;
      rca[k]	=2.*rcT2a[ki];
      rsa[k]	=2.*rsT2a[ki];
    }
    for(j=0; j < ESNp; j++){
      EZdrgt[j]	=-rc[1]*ESsn1[j]+rs[1]*EScs1[j];
      ragt	=-rca[1]*ESsn1[j]+rsa[1]*EScs1[j];
      EZdra[j]	=rc[1]*EScs1[j]+rs[1]*ESsn1[j];
      raa	=rca[0];
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	raa	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
	ragt	=k*(-rca[k]*ESsn1[kj]+rsa[k]*EScs1[kj]);
      }
      ragt	*=EZcr2;
      EZdzgt[j]	=ESsb1a[0]*EScs1[j];
      zagt	=EZcr2*ESsb2a[0]*EScs1[j];
      EZdza[j]	=ESsb1a[0]*ESsn1[j];
      zaa	=rsa[0]+ESsb2a[0]*ESsn1[j];
      dg1	=raa*EZdzgt[j]+EZdra[j]*zagt-ragt*EZdza[j]-EZdrgt[j]*zaa;
      D[j]	=dg1;
      dg0	=-EZdra[j]*rr0;
      L[j]	=L[j]*dg0+dg1*rr0;
      dg0	-=dg1/ba;
      G11[j]	*=dg0;
      G12[j]	*=dg0;
      G22[j]	*=dg0;
      dg0	=rr0/ba;
      G11[j]	+=2.*(EZdra[j]*raa+EZdza[j]*zaa)*dg0;
      G12[j]	+=(raa*EZdrgt[j]+EZdra[j]*ragt+zaa*EZdzgt[j])*dg0;
      G22[j]	+=2.*(EZdrgt[j]*ragt+EZdzgt[j]*zagt)*dg0;
    }
    D[j]	=D[0];
    L[j]	=L[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    im		=im0+i;
    ESgP2gF(ESDPc2+im,ESDPs2+im,D,ES2Mp);
    ESgP2gF(ESLc2+im,ESLs2+im,L,ES2Mp);
    m1		=im;
    for(m=0; m < ES2Mp1; m++){
      ESVc[m1]	=0.;
      ESVs[m1]	=0.;
      m1	+=ESNa1;
    }
    ESgP2gF(ESg11c2+im,ESg11s2+im,G11,ES2Mp);
    ESgP2gF(ESg12c2+im,ESg12s2+im,G12,ES2Mp);
    ESgP2gF(ESg22c2+im,ESg22s2+im,G22,ES2Mp);
#endif
    in++;
    for(i=1; i < ESNa1; i++){
      rA	=1./ESsa[i];
      b		=ESsb[in];
      ba	=ESsb1a[in];
      ki	=K+i;
      rca[0]	=rcT1a[ki];
      rsa[0]	=rsT1a[ki];
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	rca[k]	=2.*rcT1a[ki];
	rsa[k]	=2.*rsT1a[ki];
	rc[k]	=2.*k*rcT[ki];
	rs[k]	=2.*k*rsT[ki];
      }
      for(j=0; j < ESNp; j++){
	EZdra[j]	=rca[0];
	EZdrgt[j]	=0.;
	kj	=0;
	for(k=1; k < ESFp1; k++){
	  kj		+=j;
	  if(kj >= ESNp)
	    kj	-=ESNp;
	  EZdra[j]	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
	  EZdrgt[j]	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
	}
	EZdzgt[j]	=b*EScs1[j];
	EZdza[j]	=rsa[0]+ba*ESsn1[j];
	D[j]	=EZdra[j]*EZdzgt[j]-EZdrgt[j]*EZdza[j];
	if(D[j] <= 0.){
#ifdef H
	  printf("D[n=%3d][i=%2d][j=%2d]=%10.3e <= 0%c\n",n,i,j,D[j],7);
#endif
	  return(1);
	}
	L[j]	=D[j]*r0/ESsr[jni]*rA;
	V[j]	=D[j]*ESsr[jni]*rr0*rA-L[j];
	dg0	=1./(D[j]*ESsr[jni]);
	D[j]	*=rA;
	G11[j]	=ESsa[i]*(EZdra[j]*EZdra[j]+EZdza[j]*EZdza[j])*dg0;
	G12[j]	=(EZdra[j]*EZdrgt[j]+EZdza[j]*EZdzgt[j])*dg0;
	G22[j]	=(EZdrgt[j]*EZdrgt[j]+EZdzgt[j]*EZdzgt[j])*dg0*rA;
	rT1a[jni]	=EZdra[j];
	rT1gt[jni]	=EZdrgt[j];
	zT1a[jni]	=EZdza[j];
	zT1gt[jni]	=EZdzgt[j];
	jni++;
      }
      D[j]	=D[0];
      L[j]	=L[0];
      V[j]	=V[0];
      G11[j]	=G11[0];
      G12[j]	=G12[0];
      G22[j]	=G22[0];
      rT1a[jni]	=EZdra[0];
      rT1gt[jni]=EZdrgt[0];
      zT1a[jni]	=EZdza[0];
      zT1gt[jni]=EZdzgt[0];
      jni++;
      im	=im0+i;
      ESgP2gF(ESDPc+im,ESDPs+im,D,ES2Mp);
      ESgP2gF(ESLc+im,ESLs+im,L,ES2Mp);
      ESgP2gF(ESVc+im,ESVs+im,V,ES2Mp);
      ESgP2gF(ESg22c+im,ESg22s+im,G22,ES2Mp);
      ESgP2gF(ESg12c+im,ESg12s+im,G12,ES2Mp);
      ESgP2gF(ESg11c+im,ESg11s+im,G11,ES2Mp);
      in++;
    }
    i		=ESNa;
    jni		-=ESNp1;
    in--;
    ki		=K+i;
    rca[0]	=rcT2a[ki];
    rsa[0]	=rsT2a[ki];
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rca[k]	=2.*rcT2a[ki];
      rsa[k]	=2.*rsT2a[ki];
      b		=2.*k;
      rc[k]	=b*rcT1a[ki];
      rs[k]	=b*rsT1a[ki];
    }
    b		=ESsb1a[in];
    ba		=ESsb2a[in];
    for(j=0; j < ESNp; j++){
      raa	=rca[0];
      ragt	=0.;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	raa	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
	ragt	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
      }
      zaa	=rsa[0]+ba*ESsn1[j];
      zagt	=b*EScs1[j];
      dg0	=-EZdra[j]/ESsr[jni];
      dg1	=raa*EZdzgt[j]+EZdra[j]*zagt-ragt*EZdza[j]-EZdrgt[j]*zaa;
      L[j]	=(dg0-1.)*L[j]+dg1*r0/ESsr[jni];
      V[j]	=(D[j]*EZdra[j]+dg1*ESsr[jni]-D[j]*ESsr[jni])*rr0-L[j];
      dg0	-=dg1/D[j];
      G12[j]	*=dg0;
      G22[j]	*=(dg0-1.);
      dg0	+=1.;
      G11[j]	*=dg0;
      dg0	=1./(D[j]*ESsr[jni]);
      G11[j]	+=2.*(EZdra[j]*raa+EZdza[j]*zaa)*dg0*ESsa[i];
      G12[j]	+=(raa*EZdrgt[j]+EZdra[j]*ragt+zaa*EZdzgt[j]+EZdza[j]*zagt)*dg0;
      G22[j]	+=2.*(EZdrgt[j]*ragt+EZdzgt[j]*zagt)*dg0;
      D[j]	=dg1-D[j];
      jni++;
    }
    D[j]	=D[0];
    L[j]	=L[0];
    V[j]	=V[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    im		=im0+i;
    ESgP2gF(ESDPc2+im,ESDPs2+im,D,ES2Mp);
    ESgP2gF(ESLc2+im,ESLs2+im,L,ES2Mp);
    ESgP2gF(ESVc2+im,ESVs2+im,V,ES2Mp);
    ESgP2gF(ESg22c2+im,ESg22s2+im,G22,ES2Mp);
    ESgP2gF(ESg12c2+im,ESg12s2+im,G12,ES2Mp);
    ESgP2gF(ESg11c2+im,ESg11s2+im,G11,ES2Mp);
    im	=im0;
    in	=im+ESNa;
    dg0	=ESDPc2[im];
    dg1	=ESDPc2[in];
    splAA(ESDPc+im,ESDPc1+im,ESDPc2+im,&dg0,&dg1);
    dg0	=ESLc2[im];
    dg1	=ESLc2[in];
    splAA(ESLc+im,ESLc1+im,ESLc2+im,&dg0,&dg1);
    dg0	=ESVc2[im];
    dg1	=ESVc2[in];
    splAA(ESVc+im,ESVc1+im,ESVc2+im,&dg0,&dg1);
    dg0	=ESg22c2[im];
    dg1	=ESg22c2[in];
    splAA(ESg22c+im,ESg22c1+im,ESg22c2+im,&dg0,&dg1);
    dg0	=ESg12c2[im];
    dg1	=ESg12c2[in];
    splA(ESg12c+im,ESg12c2+im,&dg0,&dg1);
    dg0	=ESg11c2[im];
    dg1	=ESg11c2[in];
    splA(ESg11c+im,ESg11c2+im,&dg0,&dg1);
    for(m=1; m < 2; m++){
      im 	+=ESNa1;
      in	=im+ESNa;
      dg0	=ESDPc2[im];
      dg1	=ESDPc2[in];
      splAA2(ESDPc+im,ESDPc1+im,ESDPc2+im,ESsa,&dg1);
      dg0	=ESDPs2[im];
      dg1	=ESDPs2[in];
      splAA2(ESDPs+im,ESDPs1+im,ESDPs2+im,ESsa,&dg1);
      dg0	=ESLc2[im];
      dg1	=ESLc2[in];
      splAA2(ESLc+im,ESLc1+im,ESLc2+im,ESsa,&dg1);
      dg0	=ESLs2[im];
      dg1	=ESLs2[in];
      splAA2(ESLs+im,ESLs1+im,ESLs2+im,ESsa,&dg1);
      dg0	=ESVc2[im];
      dg1	=ESVc2[in];
      splAA2(ESVc+im,ESVc1+im,ESVc2+im,ESsa,&dg1);
      dg0	=ESVs2[im];
      dg1	=ESVs2[in];
      splAA2(ESVs+im,ESVs1+im,ESVs2+im,ESsa,&dg1);
      dg0	=ESg22c2[im];
      dg1	=ESg22c2[in];
      splAA2(ESg22c+im,ESg22c1+im,ESg22c2+im,ESsa,&dg1);
      dg0	=ESg22s2[im];
      dg1	=ESg22s2[in];
      splAA2(ESg22s+im,ESg22s1+im,ESg22s2+im,ESsa,&dg1);
      dg0	=ESg12c2[im];
      dg1	=ESg12c2[in];
      splA(ESg12c+im,ESg12c2+im,&dg0,&dg1);
      dg0	=ESg12s2[im];
      dg1	=ESg12s2[in];
      splA(ESg12s+im,ESg12s2+im,&dg0,&dg1);
      dg0	=ESg11c2[im];
      dg1	=ESg11c2[in];
      splA(ESg11c+im,ESg11c2+im,&dg0,&dg1);
      dg0	=ESg11s2[im];
      dg1	=ESg11s2[in];
      splA(ESg11s+im,ESg11s2+im,&dg0,&dg1);
    }
    for(m=2; m < ES2Mp1; m++){
      im 	+=ESNa1;
      in	=im+ESNa;
      dg1	=ESDPc2[in];
      splAA(ESDPc+im,ESDPc1+im,ESDPc2+im,ESsa,&dg1);
      dg1	=ESDPs2[in];
      splAA(ESDPs+im,ESDPs1+im,ESDPs2+im,ESsa,&dg1);
      dg1	=ESLc2[in];
      splAA(ESLc+im,ESLc1+im,ESLc2+im,ESsa,&dg1);
      dg1	=ESLs2[in];
      splAA(ESLs+im,ESLs1+im,ESLs2+im,ESsa,&dg1);
      dg1	=ESVc2[in];
      splAA(ESVc+im,ESVc1+im,ESVc2+im,ESsa,&dg1);
      dg1	=ESVs2[in];
      splAA(ESVs+im,ESVs1+im,ESVs2+im,ESsa,&dg1);
      if(m == 3){
	dg1	=ESg12c2[in];
	splAA2(ESg12c+im,ESg22c1+im,ESg12c2+im,ESsa,&dg1);
	dg1	=ESg12s2[in];
	splAA2(ESg12s+im,ESg22c1+im,ESg12s2+im,ESsa,&dg1);
	dg1	=ESg11c2[in];
	splAA2(ESg11c+im,ESg22c1+im,ESg11c2+im,ESsa,&dg1);
	dg1	=ESg11s2[in];
	splAA2(ESg11s+im,ESg22c1+im,ESg11s2+im,ESsa,&dg1);
	dg1	=ESg22c2[in];
	splAA2(ESg22c+im,ESg22c1+im,ESg22c2+im,ESsa,&dg1);
	dg1	=ESg22s2[in];
	splAA2(ESg22s+im,ESg22s1+im,ESg22s2+im,ESsa,&dg1);
      }
      else{
	dg1	=ESg22c2[in];
	splAA(ESg22c+im,ESg22c1+im,ESg22c2+im,ESsa,&dg1);
	dg1	=ESg22s2[in];
	splAA(ESg22s+im,ESg22s1+im,ESg22s2+im,ESsa,&dg1);
	dg1	=ESg12c2[in];
	splA(ESg12c+im,ESg12c2+im,ESsa,&dg1);
	dg1	=ESg12s2[in];
	splA(ESg12s+im,ESg12s2+im,ESsa,&dg1);
	dg1	=ESg11c2[in];
	splA(ESg11c+im,ESg11c2+im,ESsa,&dg1);
	dg1	=ESg11s2[in];
	splA(ESg11s+im,ESg11s2+im,ESsa,&dg1);
      }
    }
    n1	=n+ESNLt;
    while(n1 < ESNt1){
      im	=im0+ESnMp*n1;
      ji	=im0;
      for(m=0; m < ES2Mp1; m++){
	for(i=0; i < ESNa1; i++){
	  ESDPc[im]	=ESDPc[ji];
	  ESDPs[im]	=ESDPs[ji];
	  ESDPc1[im]	=ESDPc1[ji];
	  ESDPs1[im]	=ESDPs1[ji];
	  ESDPc2[im]	=ESDPc2[ji];
	  ESDPs2[im]	=ESDPs2[ji];

	  ESLc[im]	=ESLc[ji];
	  ESLs[im]	=ESLs[ji];
	  ESLc2[im]	=ESLc2[ji];
	  ESLs2[im]	=ESLs2[ji];

	  ESVc[im]	=ESVc[ji];
	  ESVs[im]	=ESVs[ji];
	  ESVc2[im]	=ESVc2[ji];
	  ESVs2[im]	=ESVs2[ji];

	  ESg22c[im]	=ESg22c[ji];
	  ESg22s[im]	=ESg22s[ji];
	  ESg22c1[im]	=ESg22c1[ji];
	  ESg22s1[im]	=ESg22s1[ji];
	  ESg22c2[im]	=ESg22c2[ji];
	  ESg22s2[im]	=ESg22s2[ji];

	  ESg12c[im]	=ESg12c[ji];
	  ESg12s[im]	=ESg12s[ji];
	  ESg12c2[im]	=ESg12c2[ji];
	  ESg12s2[im]	=ESg12s2[ji];

	  ESg11c[im]	=ESg11c[ji];
	  ESg11s[im]	=ESg11s[ji];
	  ESg11c2[im]	=ESg11c2[ji];
	  ESg11s2[im]	=ESg11s2[ji];
	  im++;
	  ji++;
	}
      }
      n1	+=ESNLt;
    }
  }
  return(0);
}

int ES3DMetrTnsrCS()
{
  int i,in,j,ji,n,n1,jni,im,im0,m,m1,mm;
  int K,k,ki,kj;
  double *rc,*rca,*rs,*rsa,*df1;

  double r0,rr0,b;
  double s,ba,rA;
  double raa,ragt,zaa,zagt,dg0,dg1;
  double *D,*L,*V,*G22,*G12,*G11;

  D	=EZdra+ESNp1;
  L	=D+ESNp1;
  rc	=L+ESNp1;

  V	=EZdza+ESNp1;
  G22	=V+ESNp1;
  rca	=G22+ESNp1;

  G12	=EZdrgt+ESNp1;
  rs	=G12+ESNp1;

  G11	=EZdzgt+ESNp1;
  rsa	=G11+ESNp1;

  for(n=0; n < ESNLt; n++){
    K		=ESnAF*n;
    r0		=ESaR0[n];
    rr0		=1./r0;
    jni		=ESnAP*n;
    in		=ESNa1*n;
    im0		=ESnMp*n;
    jni		=ESnAP*n;
    i		=0;
    b		=ESsb1a[in];
    ki		=K+ESNa1;
    rc[1]	=2.*rcT1a[ki];
    rs[1]	=2.*rsT1a[ki];
    ba		=b*rc[1];
    for(j=0; j < ESNp; j++){
      EZdrgt[j]	=-rc[1]*ESsn1[j]+rs[1]*EScs1[j];
      EZdzgt[j]	=b*EScs1[j];
      EZdra[j]	=rc[1]*EScs1[j]+rs[1]*ESsn1[j];
      EZdza[j]	=b*ESsn1[j];
      L[j]	=ba;
      dg0	=rr0/ba;
      G11[j]	=(EZdra[j]*EZdra[j]+EZdza[j]*EZdza[j])*dg0;
      G12[j]	=(EZdra[j]*EZdrgt[j]+EZdza[j]*EZdzgt[j])*dg0;
      G22[j]	=(EZdrgt[j]*EZdrgt[j]+EZdzgt[j]*EZdzgt[j])*dg0;
      rT1a[jni]	=EZdra[j];
      rT1gt[jni]=EZdrgt[j];
      zT1a[jni]	=EZdza[j];
      zT1gt[jni]=EZdzgt[j];
      jni++;
    }
    L[j]	=L[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    rT1a[jni]	=EZdra[0];
    rT1gt[jni]	=EZdrgt[0];
    zT1a[jni]	=EZdza[0];
    zT1gt[jni]	=EZdzgt[0];
    jni++;
    im		=im0;
    ESgP2gF(ESg22c+im,ESg22s+im,G22,ES2Mp);
    m1		=im;
    ESg11c[im]	=ESg22c[im];
    ESg11s[im]	=0.;
    ESg12c[im]	=0.;
    ESg12s[im]	=0.;
    ESDPc[im]	=ba;
    ESDPs[im]	=0.;
    ESLc[im]	=ba;
    ESLs[im]	=0.;
    ESVc2[im]	=0.;
    ESVs2[im]	=0.;
    m1	+=ESNa1;
    ESg11c[m1]	=0.;
    ESg11s[m1]	=0.;
    ESg12c[m1]	=0.;
    ESg12s[m1]	=0.;
    ESg22c[m1]	=0.;
    ESg22s[m1]	=0.;
    ESDPc[m1]	=0.;
    ESDPs[m1]	=0.;
    ESLc[m1]	=0.;
    ESLs[m1]	=0.;
    ESVc2[m1]	=ba*rc[1]*rr0;
    ESVs2[m1]	=ba*rs[1]*rr0;
    m1	+=ESNa1;
    ESg12c[m1]	=-ESg22s[m1];
    ESg12s[m1]	=ESg22c[m1];
    ESg11c[m1]	=-ESg22c[m1];
    ESg11s[m1]	=-ESg22s[m1];
    ESDPc[m1]	=0.;
    ESDPs[m1]	=0.;
    ESLc[m1]	=0.;
    ESLs[m1]	=0.;
    ESVc2[m1]	=0.;
    ESVs2[m1]	=0.;
    for(m=3; m < ES2Mp1; m++){
      m1	+=ESNa1;
      ESg11c[m1]	=0.;
      ESg11s[m1]	=0.;
      ESg12c[m1]	=0.;
      ESg12s[m1]	=0.;
      ESDPc[m1]		=0.;
      ESDPs[m1]		=0.;
      ESLc[m1]		=0.;
      ESLs[m1]		=0.;
      ESVc2[m1]		=0.;
      ESVs2[m1]		=0.;
    }
    ki		=K;
    rca[0]	=rcT2a[ki];
    rsa[0]	=rsT2a[ki];
    ki		+=ESNa1;
    k		=1;
    rca[1]	=2.*rcT1a[ki];
    rsa[1]	=2.*rsT1a[ki];
    ki		+=ESNa1;
    k		=2;
    rca[2]	=2.*rcT2a[ki];
    rsa[2]	=2.*rsT2a[ki];
    for(j=0; j < ESNp; j++){
      raa	=rca[0]+rca[2]*EZcs2[j]+rsa[2]*EZsn2[j];
      zaa	=rsa[0];
      ragt	=-rca[2]*EZsn2[j]+rsa[2]*EZcs2[j];
      dg1	=(b*(rca[0]+rca[2])-rsa[1]*rsa[0])*EScs1[j]
	+(b*rsa[2]+rca[1]*rsa[0])*ESsn1[j];
      dg0	=-EZdra[j]*rr0;
      D[j]	=dg1;
      L[j]	=L[j]*dg0+dg1;
      dg0	-=dg1/ba;
      G11[j]	*=dg0;
      G12[j]	*=dg0;
      G22[j]	*=dg0;
      dg0	=rr0/ba;
      G11[j]	+=2.*(EZdra[j]*raa+EZdza[j]*zaa)*dg0;
      G12[j]	+=(raa*EZdrgt[j]+EZdra[j]*ragt+zaa*EZdzgt[j])*dg0;
      G22[j]	+=2.*ragt*EZdrgt[j]*dg0;
    }
    D[j]	=D[0];
    L[j]	=L[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    im		=im0+i;
    ESgP2gF(ESDPc2+im,ESDPs2+im,D,ES2Mp);
    ESgP2gF(ESLc2+im,ESLs2+im,L,ES2Mp);
    m1		=im;
    for(m=0; m < ES2Mp1; m++){
      ESVc[m1]	=0.;
      ESVs[m1]	=0.;
      m1	+=ESNa1;
    }
    ESgP2gF(ESg11c2+im,ESg11s2+im,G11,ES2Mp);
    ESgP2gF(ESg12c2+im,ESg12s2+im,G12,ES2Mp);
    ESgP2gF(ESg22c2+im,ESg22s2+im,G22,ES2Mp);
#ifdef H
    ki		=K;
    rca[0]	=rcT2a[ki];
    rsa[0]	=rsT2a[ki];
    ki	+=ESNa1;
    rc[1]	=2.*rcT1a[ki];
    rs[1]	=2.*rsT1a[ki];
    rca[1]	=2.*rcT2a[ki];
    rsa[1]	=2.*rsT2a[ki];
    for(k=2; k < ESFp1; k++){
      ki	+=ESNa1;
      rca[k]	=2.*rcT2a[ki];
      rsa[k]	=2.*rsT2a[ki];
    }
    for(j=0; j < ESNp; j++){
      EZdrgt[j]	=-rc[1]*ESsn1[j]+rs[1]*EScs1[j];
      ragt	=-rca[1]*ESsn1[j]+rsa[1]*EScs1[j];
      EZdra[j]	=rc[1]*EScs1[j]+rs[1]*ESsn1[j];
      raa	=rca[0];
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	raa	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
	ragt	=k*(-rca[k]*ESsn1[kj]+rsa[k]*EScs1[kj]);
      }
      ragt	*=EZcr2;
      EZdzgt[j]	=ESsb1a[0]*EScs1[j];
      zagt	=EZcr2*ESsb2a[0]*EScs1[j];
      EZdza[j]	=ESsb1a[0]*ESsn1[j];
      zaa	=rsa[0]+ESsb2a[0]*ESsn1[j];
      dg1	=raa*EZdzgt[j]+EZdra[j]*zagt-ragt*EZdza[j]-EZdrgt[j]*zaa;
      D[j]	=dg1;
      dg0	=-EZdra[j]*rr0;
      L[j]	=L[j]*dg0+dg1*rr0;
      dg0	-=dg1/ba;
      G11[j]	*=dg0;
      G12[j]	*=dg0;
      G22[j]	*=dg0;
      dg0	=rr0/ba;
      G11[j]	+=2.*(EZdra[j]*raa+EZdza[j]*zaa)*dg0;
      G12[j]	+=(raa*EZdrgt[j]+EZdra[j]*ragt+zaa*EZdzgt[j])*dg0;
      G22[j]	+=2.*(EZdrgt[j]*ragt+EZdzgt[j]*zagt)*dg0;
    }
    D[j]	=D[0];
    L[j]	=L[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    im		=im0+i;
    ESgP2gF(ESDPc2+im,ESDPs2+im,D,ES2Mp);
    ESgP2gF(ESLc2+im,ESLs2+im,L,ES2Mp);
    m1		=im;
    for(m=0; m < ES2Mp1; m++){
      ESVc[m1]	=0.;
      ESVs[m1]	=0.;
      m1	+=ESNa1;
    }
    ESgP2gF(ESg11c2+im,ESg11s2+im,G11,ES2Mp);
    ESgP2gF(ESg12c2+im,ESg12s2+im,G12,ES2Mp);
    ESgP2gF(ESg22c2+im,ESg22s2+im,G22,ES2Mp);
#endif
    in++;
    for(i=1; i < ESNa1; i++){
      rA	=1./ESsa[i];
      b		=ESsb[in];
      ba	=ESsb1a[in];
      ki	=K+i;
      rca[0]	=rcT1a[ki];
      rsa[0]	=rsT1a[ki];
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	rca[k]	=2.*rcT1a[ki];
	rsa[k]	=2.*rsT1a[ki];
	rc[k]	=2.*k*rcT[ki];
	rs[k]	=2.*k*rsT[ki];
      }
      for(j=0; j < ESNp; j++){
	EZdra[j]	=rca[0];
	EZdrgt[j]	=0.;
	kj	=0;
	for(k=1; k < ESFp1; k++){
	  kj		+=j;
	  if(kj >= ESNp)
	    kj	-=ESNp;
	  EZdra[j]	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
	  EZdrgt[j]	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
	}
	EZdzgt[j]	=b*EScs1[j];
	EZdza[j]	=rsa[0]+ba*ESsn1[j];
	D[j]	=EZdra[j]*EZdzgt[j]-EZdrgt[j]*EZdza[j];
	if(D[j] <= 0.){
#ifdef H
	  printf("D[n=%3d][i=%2d][j=%2d]=%10.3e <= 0%c\n",n,i,j,D[j],7);
#endif
	  return(1);
	}
	L[j]	=D[j]*r0/ESsr[jni]*rA;
	V[j]	=D[j]*ESsr[jni]*rr0*rA-L[j];
	dg0	=1./(D[j]*ESsr[jni]);
	D[j]	*=rA;
	G11[j]	=ESsa[i]*(EZdra[j]*EZdra[j]+EZdza[j]*EZdza[j])*dg0;
	G12[j]	=(EZdra[j]*EZdrgt[j]+EZdza[j]*EZdzgt[j])*dg0;
	G22[j]	=(EZdrgt[j]*EZdrgt[j]+EZdzgt[j]*EZdzgt[j])*dg0*rA;
	rT1a[jni]	=EZdra[j];
	rT1gt[jni]	=EZdrgt[j];
	zT1a[jni]	=EZdza[j];
	zT1gt[jni]	=EZdzgt[j];
	jni++;
      }
      D[j]	=D[0];
      L[j]	=L[0];
      V[j]	=V[0];
      G11[j]	=G11[0];
      G12[j]	=G12[0];
      G22[j]	=G22[0];
      rT1a[jni]	=EZdra[0];
      rT1gt[jni]=EZdrgt[0];
      zT1a[jni]	=EZdza[0];
      zT1gt[jni]=EZdzgt[0];
      jni++;
      im	=im0+i;
      ESgP2gF(ESDPc+im,ESDPs+im,D,ES2Mp);
      ESgP2gF(ESLc+im,ESLs+im,L,ES2Mp);
      ESgP2gF(ESVc+im,ESVs+im,V,ES2Mp);
      ESgP2gF(ESg22c+im,ESg22s+im,G22,ES2Mp);
      ESgP2gF(ESg12c+im,ESg12s+im,G12,ES2Mp);
      ESgP2gF(ESg11c+im,ESg11s+im,G11,ES2Mp);
      in++;
    }
    i		=ESNa;
    jni		-=ESNp1;
    in--;
    ki		=K+i;
    rca[0]	=rcT2a[ki];
    rsa[0]	=rsT2a[ki];
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rca[k]	=2.*rcT2a[ki];
      rsa[k]	=2.*rsT2a[ki];
      b		=2.*k;
      rc[k]	=b*rcT1a[ki];
      rs[k]	=b*rsT1a[ki];
    }
    b		=ESsb1a[in];
    ba		=ESsb2a[in];
    for(j=0; j < ESNp; j++){
      raa	=rca[0];
      ragt	=0.;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	raa	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
	ragt	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
      }
      zaa	=rsa[0]+ba*ESsn1[j];
      zagt	=b*EScs1[j];
      dg0	=-EZdra[j]/ESsr[jni];
      dg1	=raa*EZdzgt[j]+EZdra[j]*zagt-ragt*EZdza[j]-EZdrgt[j]*zaa;
      L[j]	=(dg0-1.)*L[j]+dg1*r0/ESsr[jni];
      V[j]	=(D[j]*EZdra[j]+dg1*ESsr[jni]-D[j]*ESsr[jni])*rr0-L[j];
      dg0	-=dg1/D[j];
      G12[j]	*=dg0;
      G22[j]	*=(dg0-1.);
      dg0	+=1.;
      G11[j]	*=dg0;
      dg0	=1./(D[j]*ESsr[jni]);
      G11[j]	+=2.*(EZdra[j]*raa+EZdza[j]*zaa)*dg0*ESsa[i];
      G12[j]	+=(raa*EZdrgt[j]+EZdra[j]*ragt+zaa*EZdzgt[j]+EZdza[j]*zagt)*dg0;
      G22[j]	+=2.*(EZdrgt[j]*ragt+EZdzgt[j]*zagt)*dg0;
      D[j]	=dg1-D[j];
      jni++;
    }
    D[j]	=D[0];
    L[j]	=L[0];
    V[j]	=V[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    im		=im0+i;
    ESgP2gF(ESDPc2+im,ESDPs2+im,D,ES2Mp);
    ESgP2gF(ESLc2+im,ESLs2+im,L,ES2Mp);
    ESgP2gF(ESVc2+im,ESVs2+im,V,ES2Mp);
    ESgP2gF(ESg22c2+im,ESg22s2+im,G22,ES2Mp);
    ESgP2gF(ESg12c2+im,ESg12s2+im,G12,ES2Mp);
    ESgP2gF(ESg11c2+im,ESg11s2+im,G11,ES2Mp);
    im	=im0;
    in	=im+ESNa;
    dg0	=ESDPc2[im];
    dg1	=ESDPc2[in];
    splAA(ESDPc+im,ESDPc1+im,ESDPc2+im,&dg0,&dg1);
    dg0	=ESLc2[im];
    dg1	=ESLc2[in];
    splAA(ESLc+im,ESLc1+im,ESLc2+im,&dg0,&dg1);
    dg0	=ESVc2[im];
    dg1	=ESVc2[in];
    splAA(ESVc+im,ESVc1+im,ESVc2+im,&dg0,&dg1);
    dg0	=ESg22c2[im];
    dg1	=ESg22c2[in];
    splAA(ESg22c+im,ESg22c1+im,ESg22c2+im,&dg0,&dg1);
    dg0	=ESg12c2[im];
    dg1	=ESg12c2[in];
    splAA(ESg12c+im,ESg12c1+im,ESg12c2+im,&dg0,&dg1);
    dg0	=ESg11c2[im];
    dg1	=ESg11c2[in];
    splAA(ESg11c+im,ESg11c1+im,ESg11c2+im,&dg0,&dg1);
    for(m=1; m < 2; m++){
      im 	+=ESNa1;
      in	=im+ESNa;
      dg0	=ESDPc2[im];
      dg1	=ESDPc2[in];
      splAA2(ESDPc+im,ESDPc1+im,ESDPc2+im,ESsa,&dg1);
      dg0	=ESDPs2[im];
      dg1	=ESDPs2[in];
      splAA2(ESDPs+im,ESDPs1+im,ESDPs2+im,ESsa,&dg1);
      dg0	=ESLc2[im];
      dg1	=ESLc2[in];
      splAA2(ESLc+im,ESLc1+im,ESLc2+im,ESsa,&dg1);
      dg0	=ESLs2[im];
      dg1	=ESLs2[in];
      splAA2(ESLs+im,ESLs1+im,ESLs2+im,ESsa,&dg1);
      dg0	=ESVc2[im];
      dg1	=ESVc2[in];
      splAA2(ESVc+im,ESVc1+im,ESVc2+im,ESsa,&dg1);
      dg0	=ESVs2[im];
      dg1	=ESVs2[in];
      splAA2(ESVs+im,ESVs1+im,ESVs2+im,ESsa,&dg1);
      dg0	=ESg22c2[im];
      dg1	=ESg22c2[in];
      splAA2(ESg22c+im,ESg22c1+im,ESg22c2+im,ESsa,&dg1);
      dg0	=ESg22s2[im];
      dg1	=ESg22s2[in];
      splAA2(ESg22s+im,ESg22s1+im,ESg22s2+im,ESsa,&dg1);
      dg0	=ESg12c2[im];
      dg1	=ESg12c2[in];
      splAA(ESg12c+im,ESg12c1+im,ESg12c2+im,&dg0,&dg1);
      dg0	=ESg12s2[im];
      dg1	=ESg12s2[in];
      splAA(ESg12s+im,ESg12s1+im,ESg12s2+im,&dg0,&dg1);
      dg0	=ESg11c2[im];
      dg1	=ESg11c2[in];
      splAA(ESg11c+im,ESg11c1+im,ESg11c2+im,&dg0,&dg1);
      dg0	=ESg11s2[im];
      dg1	=ESg11s2[in];
      splAA(ESg11s+im,ESg11s1+im,ESg11s2+im,&dg0,&dg1);
    }
    for(m=2; m < ES2Mp1; m++){
      im 	+=ESNa1;
      in	=im+ESNa;
      dg1	=ESDPc2[in];
      splAA(ESDPc+im,ESDPc1+im,ESDPc2+im,ESsa,&dg1);
      dg1	=ESDPs2[in];
      splAA(ESDPs+im,ESDPs1+im,ESDPs2+im,ESsa,&dg1);
      dg1	=ESLc2[in];
      splAA(ESLc+im,ESLc1+im,ESLc2+im,ESsa,&dg1);
      dg1	=ESLs2[in];
      splAA(ESLs+im,ESLs1+im,ESLs2+im,ESsa,&dg1);
      dg1	=ESVc2[in];
      splAA(ESVc+im,ESVc1+im,ESVc2+im,ESsa,&dg1);
      dg1	=ESVs2[in];
      splAA(ESVs+im,ESVs1+im,ESVs2+im,ESsa,&dg1);
      if(m == 3){
	dg1	=ESg12c2[in];
	splAA2(ESg12c+im,ESg12c1+im,ESg12c2+im,ESsa,&dg1);
	dg1	=ESg12s2[in];
	splAA2(ESg12s+im,ESg12s1+im,ESg12s2+im,ESsa,&dg1);
	dg1	=ESg11c2[in];
	splAA2(ESg11c+im,ESg11c1+im,ESg11c2+im,ESsa,&dg1);
	dg1	=ESg11s2[in];
	splAA2(ESg11s+im,ESg11s1+im,ESg11s2+im,ESsa,&dg1);
	dg1	=ESg22c2[in];
	splAA2(ESg22c+im,ESg22c1+im,ESg22c2+im,ESsa,&dg1);
	dg1	=ESg22s2[in];
	splAA2(ESg22s+im,ESg22s1+im,ESg22s2+im,ESsa,&dg1);
      }
      else{
	dg1	=ESg22c2[in];
	splAA(ESg22c+im,ESg22c1+im,ESg22c2+im,ESsa,&dg1);
	dg1	=ESg22s2[in];
	splAA(ESg22s+im,ESg22s1+im,ESg22s2+im,ESsa,&dg1);
	dg1	=ESg12c2[in];
	splAA(ESg12c+im,ESg12c1+im,ESg12c2+im,ESsa,&dg1);
	dg1	=ESg12s2[in];
	splAA(ESg12s+im,ESg12s1+im,ESg12s2+im,ESsa,&dg1);
	dg1	=ESg11c2[in];
	splAA(ESg11c+im,ESg11c1+im,ESg11c2+im,ESsa,&dg1);
	dg1	=ESg11s2[in];
	splAA(ESg11s+im,ESg11s1+im,ESg11s2+im,ESsa,&dg1);
      }
    }
    n1	=n+ESNLt;
    while(n1 < ESNt1){
      im	=im0+ESnMp*n1;
      ji	=im0;
      for(m=0; m < ES2Mp1; m++){
	for(i=0; i < ESNa1; i++){
	  ESDPc[im]	=ESDPc[ji];
	  ESDPs[im]	=ESDPs[ji];
	  ESDPc1[im]	=ESDPc1[ji];
	  ESDPs1[im]	=ESDPs1[ji];
	  ESDPc2[im]	=ESDPc2[ji];
	  ESDPs2[im]	=ESDPs2[ji];

	  ESLc[im]	=ESLc[ji];
	  ESLs[im]	=ESLs[ji];
	  ESLc2[im]	=ESLc2[ji];
	  ESLs2[im]	=ESLs2[ji];

	  ESVc[im]	=ESVc[ji];
	  ESVs[im]	=ESVs[ji];
	  ESVc2[im]	=ESVc2[ji];
	  ESVs2[im]	=ESVs2[ji];

	  ESg22c[im]	=ESg22c[ji];
	  ESg22s[im]	=ESg22s[ji];
	  ESg22c1[im]	=ESg22c1[ji];
	  ESg22s1[im]	=ESg22s1[ji];
	  ESg22c2[im]	=ESg22c2[ji];
	  ESg22s2[im]	=ESg22s2[ji];

	  ESg12c[im]	=ESg12c[ji];
	  ESg12s[im]	=ESg12s[ji];
	  ESg12c1[im]	=ESg12c1[ji];
	  ESg12c2[im]	=ESg12c2[ji];
	  ESg12s1[im]	=ESg12s1[ji];
	  ESg12s2[im]	=ESg12s2[ji];

	  ESg11c[im]	=ESg11c[ji];
	  ESg11s[im]	=ESg11s[ji];
	  ESg11c1[im]	=ESg11c1[ji];
	  ESg11c2[im]	=ESg11c2[ji];
	  ESg11s1[im]	=ESg11s1[ji];
	  ESg11s2[im]	=ESg11s2[ji];
	  im++;
	  ji++;
	}
      }
      n1	+=ESNLt;
    }
  }
  return(0);
}

int ES3DMetrTnsrCS0()
{
  int i,in,j,ji,n,n1,jni,im,im0,m,m1,mm;
  int K,k,ki,kj;
  double *rc,*rca,*rs,*rsa,*df1;

  double r0,rr0,b;
  double s,ba,rA;
  double raa,ragt,zaa,zagt,dg0,dg1;
  double *D,*L,*V,*G22,*G12,*G11;
  
  D	=EZdra+ESNp1;
  L	=D+ESNp1;
  rc	=L+ESNp1;

  V	=EZdza+ESNp1;
  G22	=V+ESNp1;
  rca	=G22+ESNp1;

  G12	=EZdrgt+ESNp1;
  rs	=G12+ESNp1;

  G11	=EZdzgt+ESNp1;
  rsa	=G11+ESNp1;

  for(n=0; n < ESNLt; n++){
    K		=ESnAF*n;
    r0		=ESaR0[n];
    rr0		=1./r0;
    jni		=ESnAP*n;
    in		=ESNa1*n;
    im0		=ESnMp*n;
    jni		=ESnAP*n;
    i		=0;
    b		=ESsb1a[in];
    ki		=K+ESNa1;
    rc[1]	=2.*rcT1a[ki];
    rs[1]	=2.*rsT1a[ki];
    ba		=b*rc[1];
    for(j=0; j < ESNp; j++){
      EZdrgt[j]	=-rc[1]*ESsn1[j]+rs[1]*EScs1[j];
      EZdzgt[j]	=b*EScs1[j];
      EZdra[j]	=rc[1]*EScs1[j]+rs[1]*ESsn1[j];
      EZdza[j]	=b*ESsn1[j];
      L[j]	=ba;
      dg0	=rr0/ba;
      G11[j]	=(EZdra[j]*EZdra[j]+EZdza[j]*EZdza[j])*dg0;
      G12[j]	=(EZdra[j]*EZdrgt[j]+EZdza[j]*EZdzgt[j])*dg0;
      G22[j]	=(EZdrgt[j]*EZdrgt[j]+EZdzgt[j]*EZdzgt[j])*dg0;
      rT1a[jni]	=EZdra[j];
      rT1gt[jni]=EZdrgt[j];
      zT1a[jni]	=EZdza[j];
      zT1gt[jni]=EZdzgt[j];
      jni++;
    }
    L[j]	=L[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    rT1a[jni]	=EZdra[0];
    rT1gt[jni]	=EZdrgt[0];
    zT1a[jni]	=EZdza[0];
    zT1gt[jni]	=EZdzgt[0];
    jni++;
    im		=im0;
    ESgP2gF(ESg22c+im,ESg22s+im,G22,ES2Mp);
    m1		=im;
    ESg11c[im]	=ESg22c[im];
    ESg11s[im]	=0.;
    ESg12c[im]	=0.;
    ESg12s[im]	=0.;
    ESDPc[im]	=ba;
    ESDPs[im]	=0.;
    ESLc[im]	=ba;
    ESLs[im]	=0.;
    ESVc2[im]	=0.;
    ESVs2[im]	=0.;
    m1	+=ESNa1;
    ESg11c[m1]	=0.;
    ESg11s[m1]	=0.;
    ESg12c[m1]	=0.;
    ESg12s[m1]	=0.;
    ESg22c[m1]	=0.;
    ESg22s[m1]	=0.;
    ESDPc[m1]	=0.;
    ESDPs[m1]	=0.;
    ESLc[m1]	=0.;
    ESLs[m1]	=0.;
    ESVc2[m1]	=ba*rc[1]*rr0;
    ESVs2[m1]	=ba*rs[1]*rr0;
    m1	+=ESNa1;
    ESg12c[m1]	=-ESg22s[m1];
    ESg12s[m1]	=ESg22c[m1];
    ESg11c[m1]	=-ESg22c[m1];
    ESg11s[m1]	=-ESg22s[m1];
    ESDPc[m1]	=0.;
    ESDPs[m1]	=0.;
    ESLc[m1]	=0.;
    ESLs[m1]	=0.;
    ESVc2[m1]	=0.;
    ESVs2[m1]	=0.;
    for(m=3; m < ES2Mp1; m++){
      m1	+=ESNa1;
      ESg11c[m1]	=0.;
      ESg11s[m1]	=0.;
      ESg12c[m1]	=0.;
      ESg12s[m1]	=0.;
      ESDPc[m1]		=0.;
      ESDPs[m1]		=0.;
      ESLc[m1]		=0.;
      ESLs[m1]		=0.;
      ESVc2[m1]		=0.;
      ESVs2[m1]		=0.;
    }
    ki		=K;
    rca[0]	=rcT2a[ki];
    rsa[0]	=rsT2a[ki];
    ki		+=ESNa1;
    k		=1;
    rca[1]	=2.*rcT1a[ki];
    rsa[1]	=2.*rsT1a[ki];
    ki		+=ESNa1;
    k		=2;
    rca[2]	=2.*rcT2a[ki];
    rsa[2]	=2.*rsT2a[ki];
    for(j=0; j < ESNp; j++){
      raa	=rca[0]+rca[2]*EZcs2[j]+rsa[2]*EZsn2[j];
      zaa	=rsa[0];
      ragt	=-rca[2]*EZsn2[j]+rsa[2]*EZcs2[j];
      dg1	=(b*(rca[0]+rca[2])-rsa[1]*rsa[0])*EScs1[j]
	+(b*rsa[2]+rca[1]*rsa[0])*ESsn1[j];
      dg0	=-EZdra[j]*rr0;
      D[j]	=dg1;
      L[j]	=L[j]*dg0+dg1;
      dg0	-=dg1/ba;
      G11[j]	*=dg0;
      G12[j]	*=dg0;
      G22[j]	*=dg0;
      dg0	=rr0/ba;
      G11[j]	+=2.*(EZdra[j]*raa+EZdza[j]*zaa)*dg0;
      G12[j]	+=(raa*EZdrgt[j]+EZdra[j]*ragt+zaa*EZdzgt[j])*dg0;
      G22[j]	+=2.*ragt*EZdrgt[j]*dg0;
    }
    D[j]	=D[0];
    L[j]	=L[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    im		=im0+i;
    ESgP2gF(ESDPc2+im,ESDPs2+im,D,ES2Mp);
    ESgP2gF(ESLc2+im,ESLs2+im,L,ES2Mp);
    m1		=im;
    for(m=0; m < ES2Mp1; m++){
      ESVc[m1]	=0.;
      ESVs[m1]	=0.;
      m1	+=ESNa1;
    }
    ESgP2gF(ESg11c2+im,ESg11s2+im,G11,ES2Mp);
    ESgP2gF(ESg12c2+im,ESg12s2+im,G12,ES2Mp);
    ESgP2gF(ESg22c2+im,ESg22s2+im,G22,ES2Mp);
#ifdef H
    ki		=K;
    rca[0]	=rcT2a[ki];
    rsa[0]	=rsT2a[ki];
    ki	+=ESNa1;
    rc[1]	=2.*rcT1a[ki];
    rs[1]	=2.*rsT1a[ki];
    rca[1]	=2.*rcT2a[ki];
    rsa[1]	=2.*rsT2a[ki];
    for(k=2; k < ESFp1; k++){
      ki	+=ESNa1;
      rca[k]	=2.*rcT2a[ki];
      rsa[k]	=2.*rsT2a[ki];
    }
    for(j=0; j < ESNp; j++){
      EZdrgt[j]	=-rc[1]*ESsn1[j]+rs[1]*EScs1[j];
      ragt	=-rca[1]*ESsn1[j]+rsa[1]*EScs1[j];
      EZdra[j]	=rc[1]*EScs1[j]+rs[1]*ESsn1[j];
      raa	=rca[0];
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	raa	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
	ragt	=k*(-rca[k]*ESsn1[kj]+rsa[k]*EScs1[kj]);
      }
      ragt	*=EZcr2;
      EZdzgt[j]	=ESsb1a[0]*EScs1[j];
      zagt	=EZcr2*ESsb2a[0]*EScs1[j];
      EZdza[j]	=ESsb1a[0]*ESsn1[j];
      zaa	=rsa[0]+ESsb2a[0]*ESsn1[j];
      dg1	=raa*EZdzgt[j]+EZdra[j]*zagt-ragt*EZdza[j]-EZdrgt[j]*zaa;
      D[j]	=dg1;
      dg0	=-EZdra[j]*rr0;
      L[j]	=L[j]*dg0+dg1*rr0;
      dg0	-=dg1/ba;
      G11[j]	*=dg0;
      G12[j]	*=dg0;
      G22[j]	*=dg0;
      dg0	=rr0/ba;
      G11[j]	+=2.*(EZdra[j]*raa+EZdza[j]*zaa)*dg0;
      G12[j]	+=(raa*EZdrgt[j]+EZdra[j]*ragt+zaa*EZdzgt[j])*dg0;
      G22[j]	+=2.*(EZdrgt[j]*ragt+EZdzgt[j]*zagt)*dg0;
    }
    D[j]	=D[0];
    L[j]	=L[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    im		=im0+i;
    ESgP2gF(ESDPc2+im,ESDPs2+im,D,ES2Mp);
    ESgP2gF(ESLc2+im,ESLs2+im,L,ES2Mp);
    m1		=im;
    for(m=0; m < ES2Mp1; m++){
      ESVc[m1]	=0.;
      ESVs[m1]	=0.;
      m1	+=ESNa1;
    }
    ESgP2gF(ESg11c2+im,ESg11s2+im,G11,ES2Mp);
    ESgP2gF(ESg12c2+im,ESg12s2+im,G12,ES2Mp);
    ESgP2gF(ESg22c2+im,ESg22s2+im,G22,ES2Mp);
#endif
    in++;
    for(i=1; i < ESNa1; i++){
      rA	=1./ESsa[i];
      b		=ESsb[in];
      ba	=ESsb1a[in];
      ki	=K+i;
      rca[0]	=rcT1a[ki];
      rsa[0]	=rsT1a[ki];
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	rca[k]	=2.*rcT1a[ki];
	rsa[k]	=2.*rsT1a[ki];
	rc[k]	=2.*k*rcT[ki];
	rs[k]	=2.*k*rsT[ki];
      }
      for(j=0; j < ESNp; j++){
	EZdra[j]	=rca[0];
	EZdrgt[j]	=0.;
	kj	=0;
	for(k=1; k < ESFp1; k++){
	  kj		+=j;
	  if(kj >= ESNp) kj	-=ESNp;
	  EZdra[j]	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
	  EZdrgt[j]	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
	}
	EZdzgt[j]	=b*EScs1[j];
	EZdza[j]	=rsa[0]+ba*ESsn1[j];
	D[j]	=EZdra[j]*EZdzgt[j]-EZdrgt[j]*EZdza[j];
	if(D[j] <= 0.){
#ifdef H
	  printf("D[n=%3d][i=%2d][j=%2d]=%10.3e <= 0%c\n",n,i,j,D[j],7);
#endif
	  return(1);
	}
	L[j]	=D[j]*r0/ESsr[jni]*rA;
	V[j]	=D[j]*ESsr[jni]*rr0*rA-L[j];
	dg0	=1./(D[j]*ESsr[jni]);
	D[j]	*=rA;
	G11[j]	=ESsa[i]*(EZdra[j]*EZdra[j]+EZdza[j]*EZdza[j])*dg0;
	G12[j]	=(EZdra[j]*EZdrgt[j]+EZdza[j]*EZdzgt[j])*dg0;
	G22[j]	=(EZdrgt[j]*EZdrgt[j]+EZdzgt[j]*EZdzgt[j])*dg0*rA;
	rT1a[jni]	=EZdra[j];
	rT1gt[jni]	=EZdrgt[j];
	zT1a[jni]	=EZdza[j];
	zT1gt[jni]	=EZdzgt[j];
	jni++;
      }
      D[j]	=D[0];
      L[j]	=L[0];
      V[j]	=V[0];
      G11[j]	=G11[0];
      G12[j]	=G12[0];
      G22[j]	=G22[0];
      rT1a[jni]	=EZdra[0];
      rT1gt[jni]=EZdrgt[0];
      zT1a[jni]	=EZdza[0];
      zT1gt[jni]=EZdzgt[0];
      jni++;
      im	=im0+i;
      ESgP2gF(ESDPc+im,ESDPs+im,D,ES2Mp);
      ESgP2gF(ESLc+im,ESLs+im,L,ES2Mp);
      ESgP2gF(ESVc+im,ESVs+im,V,ES2Mp);
      ESgP2gF(ESg22c+im,ESg22s+im,G22,ES2Mp);
      ESgP2gF(ESg12c+im,ESg12s+im,G12,ES2Mp);
      ESgP2gF(ESg11c+im,ESg11s+im,G11,ES2Mp);
      in++;
    }
    i		=ESNa;
    jni		-=ESNp1;
    in--;
    ki		=K+i;
    rca[0]	=rcT2a[ki];
    rsa[0]	=rsT2a[ki];
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rca[k]	=2.*rcT2a[ki];
      rsa[k]	=2.*rsT2a[ki];
      b		=2.*k;
      rc[k]	=b*rcT1a[ki];
      rs[k]	=b*rsT1a[ki];
    }
    b		=ESsb1a[in];
    ba		=ESsb2a[in];
    for(j=0; j < ESNp; j++){
      raa	=rca[0];
      ragt	=0.;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	raa	+=rca[k]*EScs1[kj]+rsa[k]*ESsn1[kj];
	ragt	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
      }
      zaa	=rsa[0]+ba*ESsn1[j];
      zagt	=b*EScs1[j];
      dg0	=-EZdra[j]/ESsr[jni];
      dg1	=raa*EZdzgt[j]+EZdra[j]*zagt-ragt*EZdza[j]-EZdrgt[j]*zaa;
      L[j]	=(dg0-1.)*L[j]+dg1*r0/ESsr[jni];
      V[j]	=(D[j]*EZdra[j]+dg1*ESsr[jni]-D[j]*ESsr[jni])*rr0-L[j];
      dg0	-=dg1/D[j];
      G12[j]	*=dg0;
      G22[j]	*=(dg0-1.);
      dg0	+=1.;
      G11[j]	*=dg0;
      dg0	=1./(D[j]*ESsr[jni]);
      G11[j]	+=2.*(EZdra[j]*raa+EZdza[j]*zaa)*dg0*ESsa[i];
      G12[j]	+=(raa*EZdrgt[j]+EZdra[j]*ragt+zaa*EZdzgt[j]+EZdza[j]*zagt)*dg0;
      G22[j]	+=2.*(EZdrgt[j]*ragt+EZdzgt[j]*zagt)*dg0;
      D[j]	=dg1-D[j];
      jni++;
    }
    D[j]	=D[0];
    L[j]	=L[0];
    V[j]	=V[0];
    G11[j]	=G11[0];
    G12[j]	=G12[0];
    G22[j]	=G22[0];
    im		=im0+i;
    ESgP2gF(ESDPc2+im,ESDPs2+im,D,ES2Mp);
    ESgP2gF(ESLc2+im,ESLs2+im,L,ES2Mp);
    ESgP2gF(ESVc2+im,ESVs2+im,V,ES2Mp);
    ESgP2gF(ESg22c2+im,ESg22s2+im,G22,ES2Mp);
    ESgP2gF(ESg12c2+im,ESg12s2+im,G12,ES2Mp);
    ESgP2gF(ESg11c2+im,ESg11s2+im,G11,ES2Mp);
    im	=im0;
    in	=im+ESNa;
    dg0	=ESDPc2[im];
    dg1	=ESDPc2[in];
    splAA(ESDPc+im,ESDPc1+im,ESDPc2+im,&dg0,&dg1);
    dg0	=ESLc2[im];
    dg1	=ESLc2[in];
    splAA(ESLc+im,ESLc1+im,ESLc2+im,&dg0,&dg1);
    dg0	=ESVc2[im];
    dg1	=ESVc2[in];
    splAA(ESVc+im,ESVc1+im,ESVc2+im,&dg0,&dg1);
    dg0	=ESg22c2[im];
    dg1	=ESg22c2[in];
    splAA(ESg22c+im,ESg22c1+im,ESg22c2+im,&dg0,&dg1);
    dg0	=ESg12c2[im];
    dg1	=ESg12c2[in];
    splAA(ESg12c+im,ESg12c1+im,ESg12c2+im,&dg0,&dg1);
    dg0	=ESg11c2[im];
    dg1	=ESg11c2[in];
    splAA(ESg11c+im,ESg11c1+im,ESg11c2+im,&dg0,&dg1);
    for(m=1; m < 2; m++){
      im 	+=ESNa1;
      in	=im+ESNa;
      dg0	=ESDPc2[im];
      dg1	=ESDPc2[in];
      splAA2(ESDPc+im,ESDPc1+im,ESDPc2+im,ESsa,&dg1);
      dg0	=ESDPs2[im];
      dg1	=ESDPs2[in];
      splAA2(ESDPs+im,ESDPs1+im,ESDPs2+im,ESsa,&dg1);
      dg0	=ESLc2[im];
      dg1	=ESLc2[in];
      splAA2(ESLc+im,ESLc1+im,ESLc2+im,ESsa,&dg1);
      dg0	=ESLs2[im];
      dg1	=ESLs2[in];
      splAA2(ESLs+im,ESLs1+im,ESLs2+im,ESsa,&dg1);
      dg0	=ESVc2[im];
      dg1	=ESVc2[in];
      splAA2(ESVc+im,ESVc1+im,ESVc2+im,ESsa,&dg1);
      dg0	=ESVs2[im];
      dg1	=ESVs2[in];
      splAA2(ESVs+im,ESVs1+im,ESVs2+im,ESsa,&dg1);
      dg0	=ESg22c2[im];
      dg1	=ESg22c2[in];
      splAA2(ESg22c+im,ESg22c1+im,ESg22c2+im,ESsa,&dg1);
      dg0	=ESg22s2[im];
      dg1	=ESg22s2[in];
      splAA2(ESg22s+im,ESg22s1+im,ESg22s2+im,ESsa,&dg1);
      dg0	=ESg12c2[im];
      dg1	=ESg12c2[in];
      splAA(ESg12c+im,ESg12c1+im,ESg12c2+im,&dg0,&dg1);
      dg0	=ESg12s2[im];
      dg1	=ESg12s2[in];
      splAA(ESg12s+im,ESg12s1+im,ESg12s2+im,&dg0,&dg1);
      dg0	=ESg11c2[im];
      dg1	=ESg11c2[in];
      splAA(ESg11c+im,ESg11c1+im,ESg11c2+im,&dg0,&dg1);
      dg0	=ESg11s2[im];
      dg1	=ESg11s2[in];
      splAA(ESg11s+im,ESg11s1+im,ESg11s2+im,&dg0,&dg1);
    }
    for(m=2; m < ES2Mp1; m++){
      im 	+=ESNa1;
      in	=im+ESNa;
      dg1	=ESDPc2[in];
      splAA(ESDPc+im,ESDPc1+im,ESDPc2+im,ESsa,&dg1);
      dg1	=ESDPs2[in];
      splAA(ESDPs+im,ESDPs1+im,ESDPs2+im,ESsa,&dg1);
      dg1	=ESLc2[in];
      splAA(ESLc+im,ESLc1+im,ESLc2+im,ESsa,&dg1);
      dg1	=ESLs2[in];
      splAA(ESLs+im,ESLs1+im,ESLs2+im,ESsa,&dg1);
      dg1	=ESVc2[in];
      splAA(ESVc+im,ESVc1+im,ESVc2+im,ESsa,&dg1);
      dg1	=ESVs2[in];
      splAA(ESVs+im,ESVs1+im,ESVs2+im,ESsa,&dg1);
      if(m == 3){
	dg1	=ESg12c2[in];
	splAA2(ESg12c+im,ESg22c1+im,ESg12c2+im,ESsa,&dg1);
	dg1	=ESg12s2[in];
	splAA2(ESg12s+im,ESg22c1+im,ESg12s2+im,ESsa,&dg1);
	dg1	=ESg11c2[in];
	splAA2(ESg11c+im,ESg22c1+im,ESg11c2+im,ESsa,&dg1);
	dg1	=ESg11s2[in];
	splAA2(ESg11s+im,ESg22c1+im,ESg11s2+im,ESsa,&dg1);
	dg1	=ESg22c2[in];
	splAA2(ESg22c+im,ESg22c1+im,ESg22c2+im,ESsa,&dg1);
	dg1	=ESg22s2[in];
	splAA2(ESg22s+im,ESg22s1+im,ESg22s2+im,ESsa,&dg1);
      }
      else{
	dg1	=ESg22c2[in];
	splAA(ESg22c+im,ESg22c1+im,ESg22c2+im,ESsa,&dg1);
	dg1	=ESg22s2[in];
	splAA(ESg22s+im,ESg22s1+im,ESg22s2+im,ESsa,&dg1);
	dg1	=ESg12c2[in];
	splA(ESg12c+im,ESg12c2+im,ESsa,&dg1);
	dg1	=ESg12s2[in];
	splA(ESg12s+im,ESg12s2+im,ESsa,&dg1);
	dg1	=ESg11c2[in];
	splA(ESg11c+im,ESg11c2+im,ESsa,&dg1);
	dg1	=ESg11s2[in];
	splA(ESg11s+im,ESg11s2+im,ESsa,&dg1);
      }
    }
    n1	=n+ESNLt;
    while(n1 < ESNt1){
      im	=im0+ESnMp*n1;
      ji	=im0;
      for(m=0; m < ES2Mp1; m++){
	for(i=0; i < ESNa1; i++){
	  ESDPc[im]	=ESDPc[ji];
	  ESDPs[im]	=ESDPs[ji];
	  ESDPc1[im]	=ESDPc1[ji];
	  ESDPs1[im]	=ESDPs1[ji];
	  ESDPc2[im]	=ESDPc2[ji];
	  ESDPs2[im]	=ESDPs2[ji];

	  ESLc[im]	=ESLc[ji];
	  ESLs[im]	=ESLs[ji];
	  ESLc2[im]	=ESLc2[ji];
	  ESLs2[im]	=ESLs2[ji];

	  ESVc[im]	=ESVc[ji];
	  ESVs[im]	=ESVs[ji];
	  ESVc2[im]	=ESVc2[ji];
	  ESVs2[im]	=ESVs2[ji];

	  ESg22c[im]	=ESg22c[ji];
	  ESg22s[im]	=ESg22s[ji];
	  ESg22c1[im]	=ESg22c1[ji];
	  ESg22s1[im]	=ESg22s1[ji];
	  ESg22c2[im]	=ESg22c2[ji];
	  ESg22s2[im]	=ESg22s2[ji];

	  ESg12c[im]	=ESg12c[ji];
	  ESg12s[im]	=ESg12s[ji];
	  ESg12c1[im]	=ESg12c1[ji];
	  ESg12c2[im]	=ESg12c2[ji];
	  ESg12s1[im]	=ESg12s1[ji];
	  ESg12s2[im]	=ESg12s2[ji];

	  ESg11c[im]	=ESg11c[ji];
	  ESg11s[im]	=ESg11s[ji];
	  ESg11c1[im]	=ESg11c1[ji];
	  ESg11c2[im]	=ESg11c2[ji];
	  ESg11s1[im]	=ESg11s1[ji];
	  ESg11s2[im]	=ESg11s2[ji];
	  im++;
	  ji++;
	}
      }
      n1	+=ESNLt;
    }
  }
  return(0);
}

int ES3DInvG22()
{
  int i,j,ji,n,n1,im,im1,im0,m,m1,mm;
  
  double *p22c,*p22s,*pr22,*pr222;
  double *ac,*f,s;
  int *indx;

  ESReInitPlVeq(ES2Mp1,&ac,&f,&indx);
  
#ifdef H
  im0	=ESNa1*ES2Mp1;
  for(n=0; n < ESNLt; n++){
    im		=im0*n;
    p22c	=ESg22c+im;
    p22s	=ESg22s+im;
    im		*=ES2Mp1;
    pr22	=ESr22+im;
    pr222	=ESr222+im;
    i		=0;
    j		=0;
    im		=i;
    ji		=0;
    ac[ji]	=p22c[im];
    ji++;
    im		+=ESNa1;
    for(mm=0; mm < ESMp; mm++){
      ac[ji]	=2.*p22c[im];
      ji++;
      ac[ji]	=-2.*p22s[im];
      ji++;
      im	+=ESNa1;
    }
    j++;
    for(m=0; m < ESMp; m++){
      ji	=ES2Mp1*j+j;
      m1	=2*m+2;
      im	=i+ESNa1*m1;
      im1	=i;
      for(mm=m; mm < ESMp; mm++){
	if(mm == m){
	  ac[ji]	=2.*(p22c[im]+ac[0]);
	  ji++;
	  ac[ji]	=-2.*p22s[im];
	}
	else{
	  ac[ji]	=2.*(p22c[im]+p22c[im1]);
	  ji++;
	  ac[ji]	=-2.*(p22s[im]+p22s[im1]);
	}
	ji++;
	im1	+=ESNa1;
	im	+=ESNa1;
      }
      j++;
      ji	=ES2Mp1*j+j;
      im	=i+ESNa1*m1;
      im1	=i;
      mm	=m;
      ac[ji]	=2.*(ac[0]-p22c[im]);
      ji++;
      im1	+=ESNa1;
      im	+=ESNa1;
      mm++;
      while(mm < ESMp){
	ac[ji]	=2.*(p22s[im1]-p22s[im]);
	ji++;
	ac[ji]	=2.*(p22c[im1]-p22c[im]);
	ji++;
	im1	+=ESNa1;
	im	+=ESNa1;
	mm++;
      }
      j++;
    }
    if(CholDc(ac,ES2Mp1) == -1){
#ifdef H
      printf("Nonpositive G22 Matrix at n=%3d i=%2d\n",n,i);
      printf("??? ------ ac[%2d] ------\n",i);
      for(m=0; m < ES2Mp1; m++){
	ji	=ES2Mp1*m+m;
	for(mm=0; mm < m; mm++){
	  printf("%10.3s","");
	}
	for(mm=m; mm < ES2Mp1; mm++){
	  printf("%10.3e",ac[ji]);
	  ji++;
	}
	printf("\n");
      }
#endif
    }
    ji	=i;
    for(m=0; m < ES2Mp1; m++){
      for(mm=0; mm <= m ; mm++){
	pr22[ji]	=ac[m*ES2Mp1+mm];
	ji	+=ESNa1;
      }
    }
    p22c	=ESg22c+im0*n;
    p22s	=ESg22s+im0*n;
    for(i=1; i < ESNa1; i++){
      j		=0;
      im	=i;
      ji	=0;
      ac[ji]	=p22c[im];
      ji++;
      im	+=ESNa1;
      for(mm=0; mm < ESMp; mm++){
	ac[ji]	=2.*p22c[im];
	ji++;
	ac[ji]	=-2.*p22s[im];
	ji++;
	im	+=ESNa1;
      }
      j++;
      for(m=0; m < ESMp; m++){
	ji	=ES2Mp1*j+j;
	m1	=2*m+2;
	im	=i+ESNa1*m1;
	im1	=i;
	for(mm=m; mm < ESMp; mm++){
	  if(mm == m){
	    ac[ji]	=2.*(p22c[im]+ac[0]);
	    ji++;
	    ac[ji]	=-2.*p22s[im];
	  }
	  else{
	    ac[ji]	=2.*(p22c[im]+p22c[im1]);
	    ji++;
	    ac[ji]	=-2.*(p22s[im]+p22s[im1]);
	  }
	  ji++;
	  im1	+=ESNa1;
	  im	+=ESNa1;
	}
	j++;
	ji	=ES2Mp1*j+j;
	im	=i+ESNa1*m1;
	im1	=i;
	mm	=m;
	ac[ji]	=2.*(ac[0]-p22c[im]);
	ji++;
	im1	+=ESNa1;
	im	+=ESNa1;
	mm++;
	while(mm < ESMp){
	  ac[ji]	=2.*(p22s[im1]-p22s[im]);
	  ji++;
	  ac[ji]	=2.*(p22c[im1]-p22c[im]);
	  ji++;
	  im1	+=ESNa1;
	  im	+=ESNa1;
	  mm++;
	}
	j++;
      }
      if(CholDc(ac,ES2Mp1) == -1){
#ifdef H
	printf("Nonpositive G22 Matrix at n=%3d i=%2d\n",n,i);
	printf("??? ------ ac[%2d] ------\n",i);
	for(m=0; m < ES2Mp1; m++){
	  ji	=ES2Mp1*m+m;
	  for(mm=0; mm < m; mm++){
	    printf("%10.3s","");
	  }
	  for(mm=m; mm < ES2Mp1; mm++){
	    printf("%10.3e",ac[ji]);
	    ji++;
	  }
	  printf("\n");
	}
#endif
      }
      ji	=i;
      for(m=0; m < ES2Mp1; m++){
	for(mm=0; mm <= m ; mm++){
	  pr22[ji]	=ac[m*ES2Mp1+mm];
	  ji	+=ESNa1;
	}
      }
    }
    ji	=0;
    for(m=0; m < ES2Mp1; m++){
#ifdef J
      printf("??? m=%2d\n",m);
#endif
      for(mm=0; mm <= m ; mm++){
#ifdef J
	printf("???    mm=%2d\n",mm);
	for(i=0; i< ESNa1; i++){
	  printf("??? pr22[%2d]=%10.3e\n",i,*(pr22+ji+i));
	}
#endif
	splA(pr22+ji,pr222+ji,NULL,NULL);
	ji	+=ESNa1;
      }
    }
    n1	=n+ESNLt;
    while(n1 < ESNt1){
      im1	=im0*ES2Mp1*n1;
      ji	=0;
      for(m=0; m < ES2Mp1; m++){
	for(mm=0; mm <= m; mm++){
	  for(i=0; i < ESNa1; i++){
	    ESr22[im1]	=pr22[ji];
	    ESr222[im1]	=pr222[ji];
	    im1++;
	    ji++;
	  }
	}
      }
      n1	+=ESNLt;
    }
  }
#ifdef h
  {
    int kW,k,ki;
    double x[2],y[2],Wx[300],W1[300][200];
    double *p,*p1,*p2,s,t,s1;

    p		=ESr22;
    p1		=rcT1a;
    p2		=ESr222;
    y[0]	=p[0];
    y[1]	=y[0];
    x[0]	=0.;
    x[1]	=1.;
    ki		=0;
    s		=0.125*ESsa[1];
    
    ki	=0;
    k	=0;
    for(m=0; m < ES2Mp1; m++){
      for(mm=0; mm <= m; mm++){
	kW	=0;
	t	=0.;
	while(t <= x[1]+s*0.5){
	  ESSetSplA(t);
	  Wx[kW]	=t;
	  splRA(&s1,NULL,p+ki,p2+ki);
	  if(y[0] > s1)
	    y[0]	=s1;
	  if(y[1] < s1)
	    y[1]	=s1;
	  W1[k][kW]	=s1;
	  t		+=s;
	  kW++;
	}
	k++;
	ki	+=ESNa1;
      }
    }
    Scale2D(4,2,x,y,2,6,2);
    ki	=0;
    k	=0;
    for(m=0; m < ES2Mp1; m++){
      for(mm=0; mm <= m; mm++){
	Plot2d(4,ESsa,p+ki,ESNa1,6,0,4,0);
	Plot2d(4,Wx,W1[k],kW,6,0,0,0);
	k++;
	ki	+=ESNa1;
      }
    }
    CbFlush();
  }
#endif

#endif
  return(0);
}    

int ESInitgYbY(double *gyr,double *gyi,double *Yr,double *Yi,double A, int n)
{
  int i,j,m,M,M0,mm,in,k,k1,k2,kk;
  int K,ki,kj;
  double b,ba,S,sr,si,r1c,r1s;

  double *fr,*fi,*dfr,*dfi,*fmr,*fmi,*dfmr,*dfmi;
  double Fc[128],Fs[128],dFc[128],dFs[128],g22c[128],g22s[128],
  g12c[128],g12s[128];

  double gnr,gni,gnn,gnmr,gnmi,gnnm;
  double s,zr,zi,zmr,zmi;
  double f1r,f1i,f2r,f2i;

  if(ES2Mp1 > 128){
    printf("2*Mp+1=%d > 128\n",ES2Mp1);
    ESexit(0);
  }

  K	=ESnAF*n;
  ESSetSplA(A);
  fr	=ESd2Faux+ESNa1;
  fi	=fr	+ESNp1;
  dfr	=fi	+ESNp1;
  dfi	=dfr	+ESNp1;
  fmr	=dfi	+ESNp1;
  fmi	=fmr	+ESNp1;
  dfmr	=fmi	+ESNp1;
  dfmi	=dfmr	+ESNp1;

  in	=ESNa1*n;
  ki	=K+ESNa1;
  zr	=rcT[ki]+EZcr2*ESsb1a[in];
  zi	=rsT[ki];

  s	=1./(zr*zr+zi*zi);
  zr	*=s;
  zi	*=s;

  b	=ESsb1a[in];
  sr	=2.*rcT[ki];
  s	=sr+b;
  si	=2.*rsT[ki];
  s	=1./(s*s+si*si);
  gnr	=(sr*sr-si*si-b*b)*s;
  gni	=2.*sr*si*s;
  gnn	=gnr*gnr+gni*gni;
  gnmr	=gnr;
  gnmi	=gni;
  gnnm	=gnn;
  zmr	=zr;
  zmi	=zi;
  splRA(&b,&ba,ESsb+in,ESsb2a+in);
  ba	*=A;
  ki	=K;
  for(k=0; k < ESFp1; k++){
    splRA(Fc+k,dFc+k,rcT+ki,rcT2a+ki);
    splRA(Fs+k,dFs+k,rsT+ki,rsT2a+ki);
    ki	+=ESNa1;
  }
  dFc[0]	*=A;
  dFs[0]	*=A;
  S		=2.*A;
  for(k=1; k < ESFp1; k++){
    Fc[k]	*=2.;
    Fs[k]	*=2.;
    dFc[k]	*=S;
    dFs[k]	*=S;
  }
  for(j=0; j < ESNp; j++){
    fr[j]	=Fc[0];
    dfr[j]	=dFc[0];
    kj	=j;
    for(k=1; k < ESFp; k++){
      fr[j]	+=Fc[k]*EScs1[kj]+Fs[k]*ESsn1[kj];
      dfr[j]	+=dFc[k]*EScs1[kj]+dFs[k]*ESsn1[kj];
      kj		+=j;
      if(kj >= ESNp)
	kj	-=ESNp;
    }
    fi[j]	=Fs[0]+b*ESsn1[j];
    dfi[j]	=dFs[0]+ba*ESsn1[j];
    fmr[j]	=fr[j];
    fmi[j]	=fi[j];
    dfmr[j]	=dfr[j];
    dfmi[j]	=dfi[j];
  }
  fr[ESNp]	=fr[0];
  fi[ESNp]	=fi[0];
  dfr[ESNp]	=dfr[0];
  dfi[ESNp]	=dfi[0];
  fmr[ESNp]	=fmr[0];
  fmi[ESNp]	=fmi[0];
  dfmr[ESNp]	=dfmr[0];
  dfmi[ESNp]	=dfmi[0];
  M	=0;
  k	=0;
  for(kk	=0; kk < ES2Mp1; kk++){
    gyr[k]	=0.;
    gyi[k]	=0.;
    Yr[k]	=0.;
    Yi[k]	=0.;
    k++;
  }
  M++;
  ESP2F(Fc,Fs,fmr,ESMp);
  ESP2F(dFc,dFs,dfmr,ESMp);
  k	=ES2Mp1*M+ESMp;
  k2	=k;
  k1	=k;
  mm	=0;
  gyr[k]	=Fc[mm];
  gyi[k]	=Fs[mm];
  Yr[k]		=dFc[mm];
  Yi[k]		=dFs[mm];
  while(mm < ESMp){
    k2++;
    k1--;
    mm++;
    gyr[k1]	=Fc[mm];
    gyi[k1]	=Fs[mm];
    Yr[k1]	=dFc[mm];
    Yi[k1]	=dFs[mm];
    gyr[k2]	=Fc[mm];
    gyi[k2]	=-Fs[mm];
    Yr[k2]	=dFc[mm];
    Yi[k2]	=-dFs[mm];
  }
  ESP2F(Fc,Fs,fmi,ESMp);
  ESP2F(dFc,dFs,dfmi,ESMp);
  k	=ES2Mp1*M+ESMp;
  k2	=k;
  k1	=k;
  mm	=0;
  gyr[k]	-=Fs[mm];
  gyi[k]	+=Fc[mm];
  Yr[k]		-=dFs[mm];
  Yi[k]		+=dFc[mm];
  while(mm < ESMp){
    k2++;
    k1--;
    mm++;
    gyr[k1]	-=Fs[mm];
    gyi[k1]	+=Fc[mm];
    Yr[k1]	-=dFs[mm];
    Yi[k1]	+=dFc[mm];
    gyr[k2]	+=Fs[mm];
    gyi[k2]	+=Fc[mm];
    Yr[k2]	+=dFs[mm];
    Yi[k2]	+=dFc[mm];
  }
  s		=1./(1.-gnnm);
  f1r		=zmr*s;
  f1i		=zmi*s;
  f2r		=(gnmr*zmr+gnmi*zmi)*s;
  f2i		=(gnmi*zmr-gnmr*zmi)*s;
  sr		=gnmr;
  si		=gnmi;
  gnmr		=sr*gnr-si*gni;
  gnmi		=si*gnr+sr*gni;
  gnnm		*=gnn;
  sr		=zmr;
  si		=zmi;
  zmr		=sr*zr-si*zi;
  zmi		=si*zr+sr*zi;
  k		=ES2Mp1*M+ESMp;
  sr		=gyr[k];
  si		=gyi[k];
  gyr[k]	=(f1r-f2r)*sr-(f1i+f2i)*si;
  gyi[k]	=(f1i-f2i)*sr+(f1r+f2r)*si;
  sr		=Yr[k];
  si		=Yi[k];
  Yr[k]		=(f1r-f2r)*sr-(f1i+f2i)*si;
  Yi[k]		=(f1i-f2i)*sr+(f1r+f2r)*si;
  for(mm=1; mm < ESMp1; mm++){
    k2		=k+mm;
    k1		=k-mm;
    sr		=gyr[k2];
    si		=gyi[k2];
    r1c		=gyr[k1];
    r1s		=gyi[k1];
    gyr[k2]	=f1r*sr-f1i*si-f2r*r1c-f2i*r1s;
    gyi[k2]	=f1i*sr+f1r*si-f2i*r1c+f2r*r1s;
    gyr[k1]	=f1r*r1c-f1i*r1s-f2r*sr-f2i*si;
    gyi[k1]	=f1i*r1c+f1r*r1s-f2i*sr+f2r*si;
    sr		=Yr[k2];
    si		=Yi[k2];
    r1c		=Yr[k1];
    r1s		=Yi[k1];
    Yr[k2]	=f1r*sr-f1i*si-f2r*r1c-f2i*r1s;
    Yi[k2]	=f1i*sr+f1r*si-f2i*r1c+f2r*r1s;
    Yr[k1]	=f1r*r1c-f1i*r1s-f2r*sr-f2i*si;
    Yi[k1]	=f1i*r1c+f1r*r1s-f2i*sr+f2r*si;
  }	
  for(M=2; M < ESMp1; M++){
    for(j=0; j < ESNp1; j++){
      sr	=fmr[j];
      si	=fmi[j];
      dfmr[j]	=M*(dfr[j]*sr-dfi[j]*si);
      dfmi[j]	=M*(dfr[j]*si+dfi[j]*sr);
      fmr[j]	=fr[j]*sr-fi[j]*si;
      fmi[j]	=fr[j]*si+fi[j]*sr;
    }
    ESP2F(Fc,Fs,fmr,ESMp);
    ESP2F(dFc,dFs,dfmr,ESMp);
    k		=ES2Mp1*M+ESMp;
    k2		=k;
    k1		=k;
    mm		=0;
    gyr[k]	=Fc[mm];
    gyi[k]	=0.;
    Yr[k]	=dFc[mm];
    Yi[k]	=0.;
    while(mm < ESMp){
      k2++;
      k1--;
      mm++;
      gyr[k1]	=Fc[mm];
      gyi[k1]	=Fs[mm];
      Yr[k1]	=dFc[mm];
      Yi[k1]	=dFs[mm];
      gyr[k2]	=Fc[mm];
      gyi[k2]	=-Fs[mm];
      Yr[k2]	=dFc[mm];
      Yi[k2]	=-dFs[mm];
    }
    ESP2F(Fc,Fs,fmi,ESMp);
    ESP2F(dFc,dFs,dfmi,ESMp);
    k		=ES2Mp1*M+ESMp;
    k2		=k;
    k1		=k;
    mm		=0;
    gyi[k]	+=Fc[mm];
    Yi[k]	+=dFc[mm];
    while(mm < ESMp){
      k2++;
      k1--;
      mm++;
      gyr[k1]	-=Fs[mm];
      gyi[k1]	+=Fc[mm];
      Yr[k1]	-=dFs[mm];
      Yi[k1]	+=dFc[mm];
      gyr[k2]	+=Fs[mm];
      gyi[k2]	+=Fc[mm];
      Yr[k2]	+=dFs[mm];
      Yi[k2]	+=dFc[mm];
    }
    s		=1./(1.-gnnm);
    f1r		=zmr*s;
    f1i		=zmi*s;
    f2r		=(gnmr*zmr+gnmi*zmi)*s;
    f2i		=(gnmi*zmr-gnmr*zmi)*s;
    sr		=gnmr;
    si		=gnmi;
    gnmr	=sr*gnr-si*gni;
    gnmi	=si*gnr+sr*gni;
    gnnm	*=gnn;
    sr		=zmr;
    si		=zmi;
    zmr		=sr*zr-si*zi;
    zmi		=si*zr+sr*zi;
    k		=ES2Mp1*M+ESMp;
    sr		=gyr[k];
    si		=gyi[k];
    gyr[k]	=(f1r-f2r)*sr-(f1i+f2i)*si;
    gyi[k]	=(f1i-f2i)*sr+(f1r+f2r)*si;
    sr		=Yr[k];
    si		=Yi[k];
    Yr[k]	=(f1r-f2r)*sr-(f1i+f2i)*si;
    Yi[k]	=(f1i-f2i)*sr+(f1r+f2r)*si;
    for(mm=1; mm < ESMp1; mm++){
      k2	=k+mm;
      k1	=k-mm;
      sr	=gyr[k2];
      si	=gyi[k2];
      r1c	=gyr[k1];
      r1s	=gyi[k1];
      gyr[k2]	=f1r*sr-f1i*si-f2r*r1c-f2i*r1s;
      gyi[k2]	=f1i*sr+f1r*si-f2i*r1c+f2r*r1s;
      gyr[k1]	=f1r*r1c-f1i*r1s-f2r*sr-f2i*si;
      gyi[k1]	=f1i*r1c+f1r*r1s-f2i*sr+f2r*si;
      sr	=Yr[k2];
      si	=Yi[k2];
      r1c	=Yr[k1];
      r1s	=Yi[k1];
      Yr[k2]	=f1r*sr-f1i*si-f2r*r1c-f2i*r1s;
      Yi[k2]	=f1i*sr+f1r*si-f2i*r1c+f2r*r1s;
      Yr[k1]	=f1r*r1c-f1i*r1s-f2r*sr-f2i*si;
      Yi[k1]	=f1i*r1c+f1r*r1s-f2i*sr+f2r*si;
    }	
  }
  k	=ESnMp*n;
  splRA(g22c,NULL,ESg22c+k,ESg22c2+k);
  splRA(g12c,NULL,ESg12c+k,ESg12c2+k);
  S		=1./A;
  g22c[0]	*=S;
  g12c[0]	*=S;
  g22s[0]	=0.;
  g12s[0]	=0.;
  Fc[0]		=0.;
  Fs[0]		=0.;
  k	 	+=ESNa1;
  for(m=1; m < ES2Mp1; m++){
    Fc[m]	=0.;
    Fs[m]	=0.;
    splRA(g22c+m,NULL,ESg22c+k,ESg22c2+k);
    splRA(g22s+m,NULL,ESg22s+k,ESg22s2+k);
    g22c[m]	*=S;
    g22s[m]	*=S;
    splRA(g12c+m,NULL,ESg12c+k,ESg12c2+k);
    splRA(g12s+m,NULL,ESg12s+k,ESg12s2+k);
    g12c[m]	*=S;
    g12s[m]	*=S;
    k 	+=ESNa1;
  }

  for(M=1; M < ESMp1; M++){
    M0		=ES2Mp1*M;
    for(k=0; k < ES2Mp1; k++){
      mm	=M0+k;
      m		=k-ESMp;
      si	=g12c[0]*m;
      sr	=g22c[0];
      Fc[k]	=sr*Yr[mm]+si*gyi[mm];
      Fs[k]	=sr*Yi[mm]-si*gyr[mm];
      for(kk=1, k1=k-kk; k1 >= 0; kk++,k1--){
	mm	=M0+k1;
	m	=k1-ESMp;
	sr	=m*g12c[kk];
	si	=m*g12s[kk];
	Fc[k]	+=g22c[kk]*Yr[mm]+g22s[kk]*Yi[mm]+sr*gyi[mm]-si*gyr[mm];
	Fs[k]	+=g22c[kk]*Yi[mm]-g22s[kk]*Yr[mm]-sr*gyr[mm]-si*gyi[mm];
      }
      for(kk=1, k1=k+kk; k1 < ES2Mp1; kk++, k1++){
	mm	=M0+k1;
	m	=k1-ESMp;
	sr	=m*g12c[kk];
	si	=m*g12s[kk];
	Fc[k]	+=g22c[kk]*Yr[mm]-g22s[kk]*Yi[mm]+sr*gyi[mm]+si*gyr[mm];
	Fs[k]	+=g22c[kk]*Yi[mm]+g22s[kk]*Yr[mm]-sr*gyr[mm]+si*gyi[mm];
      }
    }
    for(k=0; k < ES2Mp1; k++){
      mm	=M0+k;
      Yr[mm]	=Fc[k];
      Yi[mm]	=Fs[k];
    }
  }
  return(0);
}

#ifndef mzl_Esc2RZ
int ES2DMapFlx2Lab(double *r, double *z
		   ,double *EZdra, double *EZdza
		   ,double *drp, double *dzp
		   ,double a, double gt)
{

  double cs,sn;
  double rc,rca,rs,rsa;
  int k,ki;

  cs	=cos(gt);
  sn	=sin(gt);
  ESSetSplA(a);
  splRA(&rc,&rca,ESsb,ESsb2a);
  splRA(&rs,&rsa,rsT,rsT2a);
  *z	=ESaZ0[0]+rs+rc*sn;
  *EZdza	=rsa+rca*sn;
  *dzp	=rc*cs;

  splRA(&rc,&rca,rcT,rcT2a);
  *r	=ESaR0[0]+rc;
  *EZdra	=rca;
  *drp	=0.;
  ki	=0;
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    splRA(&rc,&rca,rcT+ki,rcT2a+ki);
    splRA(&rs,&rsa,rsT+ki,rsT2a+ki);
    cs	=cos(k*gt);
    sn	=sin(k*gt);
    *r	+=2.*(rc*cs+rs*sn);
    *EZdra+=2.*(rca*cs+rsa*sn);
    *drp+=2.*k*(-rc*sn+rs*cs);
  }
  return(0);
}

int ESSetMapBox()
{
  int j,M;
  int i,k,ki;
  double cs,sn,drp,d2rp,gt;

  Zbox[1]	=ESaZ0[0]+rsT[ESNa]+ESsb[ESNa];
  Zbox[0]	=ESaZ0[0]+rsT[ESNa]-ESsb[ESNa];
 
  gt	=EZcgp;
  i	=0;
  M	=ESFp1;
  j	=0;
  do{
    ki	=ESNa;
    drp	=0.;
    d2rp=0.;
    for(k=1; k < M; k++){
      ki	+=ESNa1;
      cs	=cos(k*gt);
      sn	=sin(k*gt);
      drp	+=2.*k*(rsT[ki]*cs-rcT[ki]*sn);
      d2rp	-=2.*k*k*(rcT[ki]*cs+rsT[ki]*sn);
      if(i == 0 && (rcT[ki] != 0. || rsT[ki] != 0.)) j =k+1;
    }
    if(i== 0) M=j;
    drp	/=-d2rp;
    gt	+=drp;
    i++;
  }while(i < 10 && fabs(drp) > 1e-4);

  ki	=ESNa;
  Rbox[0]	=ESaR0[0]+rcT[ki];
  Rbox[1]	=Rbox[0];
  Rbox[2]	=Rbox[0];
  Rbox[3]	=Rbox[0];
  gtrmn		=gt;
  Zbox[2]	=ESaZ0[0]+rsT[ESNa]+ESsb[ESNa]*sin(gt);
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    cs	=cos(k*gt);
    sn	=sin(k*gt);
    Rbox[0]+=2.*(rcT[ki]*cs+rsT[ki]*sn);
  }
  gt	=0.;
  i	=0;
  do{
    ki	=ESNa;
    drp	=0.;
    d2rp=0.;
    for(k=1; k < M; k++){
      ki	+=ESNa1;
      cs	=cos(k*gt);
      sn	=sin(k*gt);
      drp	+=2.*k*(rsT[ki]*cs-rcT[ki]*sn);
      d2rp	-=2.*k*k*(rcT[ki]*cs+rsT[ki]*sn);
    }
    drp	/=-d2rp;
    gt	+=drp;
    i++;
  }while(i < 10 && fabs(drp) > 1e-4);
  ki	=ESNa;
  gtrmx	=gt;
  Zbox[3]	=ESaZ0[0]+rsT[ESNa]+ESsb[ESNa]*sin(gt);
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    cs	=cos(k*gt);
    sn	=sin(k*gt);
    Rbox[1]+=2.*(rcT[ki]*cs+rsT[ki]*sn);
  }
  cs	=0.;
  sn	=1.;
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    Rbox[2]+=2.*(rcT[ki]*cs-rsT[ki]*sn);
    Rbox[3]+=2.*(rcT[ki]*cs+rsT[ki]*sn);
    gt	=cs;
    cs	=-sn;
    sn	=gt;
  }
  Zmn	=Zbox[0];
  Zmx	=Zbox[1];
  Rmn	=Rbox[0];
  Rmx	=Rbox[1];
  return(0);
}

int ESGetMapBox(double *R,double *Z)
{
  int j,n,M;
  int i,k,ki;
  double cs,sn,drp,d2rp,gt;

  j	=ESNa*ESNp1;
  n	=j+ESNp1;
  Z[1]	=ESaZ0[0]+rsT[ESNa]+ESsb[ESNa];
  Z[0]	=ESaZ0[0]+rsT[ESNa]-ESsb[ESNa];
  
  R[0]	=ESsr[j];
  k	=1;
  for(j++; j < n; j++){
    if(R[0] > ESsr[j]){
      R[0]=ESsr[j];
      i	=j;
    }
  }
  gt	=ESgt[i-ESNa*ESNp1];
  i	=0;
  M	=ESFp1;
  j	=0;
  do{
    ki	=ESNa;
    drp	=0.;
    d2rp=0.;
    for(k=1; k < M; k++){
      ki	+=ESNa1;
      cs	=cos(k*gt);
      sn	=sin(k*gt);
      drp	+=2.*k*(rsT[ki]*cs-rcT[ki]*sn);
      d2rp	-=2.*k*k*(rcT[ki]*cs+rsT[ki]*sn);
      if(i == 0 && (rcT[ki] != 0. || rsT[ki] != 0.)) j =k+1;
    }
    if(i== 0) M=j;
    drp	/=-d2rp;
    gt	+=drp;
    i++;
  }while(i < 10 && fabs(drp) > 1e-4);

  ki	=ESNa;
  R[0]	=ESaR0[0]+rcT[ki];
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    cs	=cos(k*gt);
    sn	=sin(k*gt);
    R[0]+=2.*(rcT[ki]*cs+rsT[ki]*sn);
  }
  j	=ESNa*ESNp1;
  R[1]	=ESsr[j];
  for(j++; j < n; j++){
    if(R[1] < ESsr[j]){
      R[1]=ESsr[j];
      i	=j;
    }
  }
  gt	=ESgt[i-ESNa*ESNp1];
  i	=0;
  do{
    ki	=ESNa;
    drp	=0.;
    d2rp=0.;
    for(k=1; k < M; k++){
      ki	+=ESNa1;
      cs	=cos(k*gt);
      sn	=sin(k*gt);
      drp	+=2.*k*(rsT[ki]*cs-rcT[ki]*sn);
      d2rp	-=2.*k*k*(rcT[ki]*cs+rsT[ki]*sn);
    }
    drp	/=-d2rp;
    gt	+=drp;
    i++;
  }while(i < 10 && fabs(drp) > 1e-4);
  ki	=ESNa;
  R[1]	=ESaR0[0]+rcT[ki];
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    cs	=cos(k*gt);
    sn	=sin(k*gt);
    R[1]+=2.*(rcT[ki]*cs+rsT[ki]*sn);
  }
  return(0);
}

int ES1DMapZ2gt(double *gt,double *R1,double *R2, double z)
{
  double t,cs,sn,s;
  int k,ki;

  ki	=ESNa;
  t	=(z-ESaZ0[0]-rsT[ki])/ESsb[ESNa];
  if(t > 1. || t < -1.){
    return(1);
  }
  *gt	=asin(t);

  *R2	=ESaR0[0]+rcT[ki];
  *R1	=*R2;
  s	=-1;
  for(k=1; k < ESFp1; k++){
    t	=*gt*k;
    cs	=cos(t);
    sn	=sin(t);
    ki	+=ESNa1;
    *R2	+=2.*(rcT[ki]*cs+rsT[ki]*sn);
    *R1	+=2.*(rcT[ki]*cs-rsT[ki]*sn)*s;
    s	=-s;
  }
  return(0);
}

int ES1DMapR2gt(double *t1,double *t2,double *Z1,double *Z2, double R)
{
  double cs,sn,r,drp,gt;
  int i,k,ki,j,M;

  if(R < Rmn || R > Rmx){
    return(1);
  }
  if(R == Rmn){
    *Z1	=Zrmn;
    *Z2	=Zrmn;
    return(0);
  }
  if(R == Rmx){
    *Z1	=Zrmx;
    *Z2	=Zrmx;
    return(0);
  }
  if(R > Rbox[3]){
    gt	=acos((R-Rbox[3])/(Rbox[1]-Rbox[3]));
    gt	=EZcr2*(EZcgp+(2.*gt-EZcgp)*(EZcgp-2.*gtrmx)/EZcgp);
  }
  else{
    gt	=acos((R-Rbox[3])/(Rbox[3]-Rbox[0]));
    gt	=EZcr2*(EZcgp-(2.*gt-EZcgp)*(EZcgp-2.*gtrmn)/EZcgp);
  }
  j	=0;
  M	=ESFp1;
  i	=0;
  do{
    ki	=ESNa;
    r	=ESaR0[0]+rcT[ki]-R;
    drp	=0.;
    for(k=1; k < M; k++){
      ki	+=ESNa1;
      cs	=cos(k*gt);
      sn	=sin(k*gt);
      r		+=2.*(rcT[ki]*cs+rsT[ki]*sn);
      drp	+=2.*k*(-rcT[ki]*sn+rsT[ki]*cs);
      if(i == 0 && (rcT[ki] != 0. || rsT[ki] != 0.)) j =k+1;
    }
    if(i== 0) M=j;
    drp	=-r/drp;
    gt	+=drp;
    i++;
  }while(i < 15 && fabs(drp) > 1e-6);
  if(i	< 15){
    *t2	=gt;
    *Z2	=ESaZ0[0]+rsT[ESNa]+ESsb[ESNa]*sin(gt);
  }
  else{
    return(1);
  }

  if(R > Rbox[3]){
    gt	=acos((R-Rbox[2])/(Rbox[1]-Rbox[2]));
    gt	=-EZcr2*(EZcgp+(2.*gt-EZcgp)*(EZcgp-2.*gtrmx)/EZcgp);
  }
  else{
    gt	=acos((R-Rbox[2])/(Rbox[2]-Rbox[0]));
    gt	=-EZcr2*(EZcgp-(2.*gt-EZcgp)*(EZcgp-2.*gtrmn)/EZcgp);
  }
  i	=0;
  do{
    ki	=ESNa;
    r	=ESaR0[0]+rcT[ki]-R;
    drp	=0.;
    for(k=1; k < M; k++){
      ki	+=ESNa1;
      cs	=cos(k*gt);
      sn	=sin(k*gt);
      r		+=2.*(rcT[ki]*cs+rsT[ki]*sn);
      drp	+=2.*k*(-rcT[ki]*sn+rsT[ki]*cs);
    }
    drp	=-r/drp;
    gt	+=drp;
    i++;
  }while(i < 15 && fabs(drp) > 1e-6);
  if(i	< 15){
    *t1	=gt;
    *Z1	=ESaZ0[0]+rsT[ESNa]+ESsb[ESNa]*sin(gt);
    return(0);
  }
  return(1);
}

int ES2DMapLab2Flx(double *A,double *t,double R, double Z)
{
  int Fl;
  double D,Rt;
  double cs,sn;
  double rc,rca,rc2,rs,rsa,rs2;
  int k,ki,i,M;
  double r,z,EZdra,dr2,EZdza,dz2,drp,dzp,a,gt;

  R	-=ESaR0[0];
  Z	-=ESaZ0[0];
  Fl	=0;
  a	=*A;
  gt	=*t;
  
#ifdef H
  ESSetSplA(a);
  r	=0.;
  Rt	=0.;
  ki	=0;
  EZdra	=Z < 0. ? -1. : 1.;
  M=2;
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    splRA(&rc,NULL,rcT+ki,rcT2a+ki);
    splRA(&rs,NULL,rsT+ki,rsT2a+ki);
    cs	=cos(k*gt);
    sn	=sin(k*gt);
    r	+=(rc*cs+rs*sn);
    switch(k%4){
    case 0:
      Rt	+=rc;
      break;
    case 1:
      Rt	+=EZdra*rs;
      break;
    case 2:
      Rt	-=rc;
      break;
    case 3:
      Rt	-=EZdra*rs;
      break;
    }
  }
  splRA(&rc,NULL,rcT,rcT2a);
  splRA(&rs,NULL,rcT,rcT2a);
  Rt	=rc+2.*Rt;
  r	=rc+2.*r;
  if((R-Rt)*(r-Rt) < 0.) gt=EZcgp-gt;
#endif
  if(a == 0.){
    EZdra	=Z/ESsb1a[0];
    ki	=ESNa1;
    drp	=(R-2.*rsT1a[ki]*EZdra)/(2.*rcT1a[ki]);
    a  	=sqrt(drp*drp+EZdra*EZdra);
    gt	=asin(EZdra/a);
    if(drp < 0.) gt=EZcgp-gt;
  }

  i	=0;
  do{
    ESSetSplA(a);
    cs	=cos(gt);
    sn	=sin(gt);
    splRA2(&rc,&rca,&rc2,ESsb,ESsb2a);
    splRA2(&rs,&rsa,&rs2,rsT,rsT2a);
    z	=rs+rc*sn-Z;
    EZdza	=rsa+rca*sn;
    dz2	=rs2+rc2*sn;
    dzp	=rc*cs;
    
    r	=0.;
    EZdra	=0.;
    dr2	=0.;
    drp	=0.;
    ki	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA2(&rc,&rca,&rc2,rcT+ki,rcT2a+ki);
      splRA2(&rs,&rsa,&rs2,rsT+ki,rsT2a+ki);
      cs	=cos(k*gt);
      sn	=sin(k*gt);
      r		+=(rc*cs+rs*sn);
      EZdra	+=(rca*cs+rsa*sn);
      dr2	+=(rc2*cs+rs2*sn);
      drp	+=k*(-rc*sn+rs*cs);
    }
    splRA2(&rc,&rca,&rc2,rcT,rcT2a);
    r	=rc+2.*r-R;
    EZdra	=rca+2.*EZdra;
    dr2	=rc2+2.*dr2;
    drp	*=2.;
    D	=EZdra*dzp-drp*EZdza;
    drp	=2.*(z*drp-r*dzp);
    drp	/=sqrt(D*D+(dr2*dzp-drp*dz2)*drp)+D;
    a	+=drp;
    EZdra	=(r*EZdza-z*EZdra-EZcr2*(EZdra*dz2-dr2*EZdza)*drp*drp)/D;
    gt	+=EZdra;
    if(a < 0.){
      Fl	=2;
      a		-=drp;
    }
    if(a > 1.){
      a	=1.;
#ifdef H
      return(-1);
#endif
    }
    i++;
  }while(i < 15 && fabs(EZdra)+fabs(drp) > 1e-8);
  if(i < 15){
    while(gt > EZc2gp) gt-=EZc2gp;
    while(gt < 0.) gt+=EZc2gp;
    *A	=a;
    *t	=gt;
  }
  else{
    Fl	=1;
  }
  if(Fl)
    printf("Error[%d]: r=%10.3e z=%10.3e not converted to a=%10.3e gt=%10.3e\n"
	   ,Fl,R,Z,a,gt);
  return(Fl);
}

int ESGetLabMetrics(double *dRa,double *dRp,double *dZa,double *dZp
		    ,double *A,double *t,double R, double Z)
{
  int Fl;
  double D,Rt;
  double cs,sn;
  double rc,rca,rc2,rs,rsa,rs2;
  int k,ki,i;
  double r,z,EZdra,dr2,EZdza,dz2,drp,dzp,a,gt,da,dt;
  
  R	-=ESaR0[0];
  Z	-=ESaZ0[0];
  Fl	=0;
  a	=*A;
  gt	=*t;
  
  ESSetSplA(a);
  r	=0.;
  Rt	=0.;
  ki	=0;
  EZdra	=Z < 0. ? -1. : 1.;
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    splRA(&rc,NULL,rcT+ki,rcT2a+ki);
    splRA(&rs,NULL,rsT+ki,rsT2a+ki);
    cs	=cos(k*gt);
    sn	=sin(k*gt);
    r	+=(rc*cs+rs*sn);
    switch(k%4){
    case 0:
      Rt	+=rc;
      break;
    case 1:
      Rt	+=EZdra*rs;
      break;
    case 2:
      Rt	-=rc;
      break;
    case 3:
      Rt	-=EZdra*rs;
      break;
    }
  }
  splRA(&rc,NULL,rcT,rcT2a);
  splRA(&rs,NULL,rcT,rcT2a);
  Rt	=rc+2.*Rt;
  r	=rc+2.*r;
  if((R-Rt)*(r-Rt) < 0.){
    gt	=EZcgp-gt;
  }
  i	=0;
  do{
    ESSetSplA(a);
    cs	=cos(gt);
    sn	=sin(gt);
    splRA2(&rc,&rca,&rc2,ESsb,ESsb2a);
    splRA2(&rs,&rsa,&rs2,rsT,rsT2a);
    z	=rs+rc*sn-Z;
    EZdza	=rsa+rca*sn;
    dz2	=rs2+rc2*sn;
    dzp	=rc*cs;
    
    r	=0.;
    EZdra	=0.;
    dr2	=0.;
    drp	=0.;
    ki	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA2(&rc,&rca,&rc2,rcT+ki,rcT2a+ki);
      splRA2(&rs,&rsa,&rs2,rsT+ki,rsT2a+ki);
      cs	=cos(k*gt);
      sn	=sin(k*gt);
      r		+=(rc*cs+rs*sn);
      EZdra	+=(rca*cs+rsa*sn);
      dr2	+=(rc2*cs+rs2*sn);
      drp	+=k*(-rc*sn+rs*cs);
    }
    splRA2(&rc,&rca,&rc2,rcT,rcT2a);
    r	=rc+2.*r-R;
    EZdra	=rca+2.*EZdra;
    dr2	=rc2+2.*dr2;
    drp	*=2.;
    D	=EZdra*dzp-drp*EZdza;
    da	=2.*(z*drp-r*dzp);
    da	/=sqrt(D*D+(dr2*dzp-da*dz2)*da)+D;
    a	+=da;
    dt	=(r*EZdza-z*EZdra-EZcr2*(EZdra*dz2-dr2*EZdza)*da*da)/D;
    gt	+=dt;
    if(a < 0.){
      Fl	=2;
      a	-=da;
    }
    if(a > 1.){
      return(-1);
    }
    i++;
  }while(i < 10 && fabs(dt)+fabs(da) > 0.001);
  if(i < 10){
    *A	=a;
    *t	=gt;
    *dRa	=EZdra;
    *dRp	=drp;
    *dZa	=EZdza;
    *dZp	=dzp;
    return(Fl);
  }
  return(1);
}

int ESGetMetrics(double *dRa,double *dRp,double *dZa,double *dZp
		    ,double a,double gt)
{
  int k,ki;
  double cs,sn;
  double rc,rca,rs,rsa;
  double EZdra,drp;
  
  ESSetSplA(a);
  cs	=cos(gt);
  sn	=sin(gt);
  splRA(&rc,&rca,ESsb,ESsb2a);
  splRA(&rs,&rsa,rsT,rsT2a);
  *dZa	=rsa+rca*sn;
  *dZp	=rc*cs;
  EZdra	=0.;
  drp	=0.;
  ki	=0;
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    splRA(&rc,&rca,rcT+ki,rcT2a+ki);
    splRA(&rs,&rsa,rsT+ki,rsT2a+ki);
    cs	=cos(k*gt);
    sn	=sin(k*gt);
    EZdra	+=(rca*cs+rsa*sn);
    drp	+=k*(-rc*sn+rs*cs);
  }
  splRA(&rc,&rca,rcT,rcT2a);
  *dRa	=rca+2.*EZdra;
  *dRp	=2.*drp;
  return(0);
}

int ESGetFullMetrics(double *Ra,double *Rp,double *Raa,double *Rap,double *Rpp
		     ,double *Za,double *Zp,double *Zaa,double *Zap,double *Zpp
		     ,double a,double gt)
{
  int k,ki;
  double cs,sn;
  double rc,rca,rcaa,rs,rsa,rsaa;
  double ra,rp,raa,rap,rpp;
  
  ESSetSplA(a);
  cs	=cos(gt);
  sn	=sin(gt);
  splRA2(&rc,&rca,&rcaa,ESsb,ESsb2a);
  splRA2(&rs,&rsa,&rsaa,rsT,rsT2a);
  *Za	=rsa+rca*sn;
  *Zaa	=rsaa+rcaa*sn;
  *Zp	=rc*cs;
  *Zap	=rca*cs;
  *Zpp	=-rc*sn;
  ra	=0.;
  raa	=0.;
  rp	=0.;
  rap	=0.;
  rpp	=0.;
  ki	=0;
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    splRA2(&rc,&rca,&rcaa,rcT+ki,rcT2a+ki);
    splRA2(&rs,&rsa,&rsaa,rsT+ki,rsT2a+ki);
    cs	=cos(k*gt);
    sn	=sin(k*gt);
    ra	+=rca*cs+rsa*sn;
    raa	+=rcaa*cs+rsaa*sn;
    sn	*=k;
    cs	*=k;
    rp	+=-rc*sn+rs*cs;
    rap	+=-rca*sn+rsa*cs;
    rpp	-=k*(rc*cs+rs*sn);
  }
  splRA2(&rc,&rca,&rcaa,rcT,rcT2a);
  *Ra	=rca+2.*ra;
  *Raa	=rcaa+2.*raa;
  *Rp	=2.*rp;
  *Rap	=2.*rap;
  *Rpp	=2.*rpp;
  return(0);
}

extern int EcNjet;
extern double EcRjet[],EcZjet[];
extern double *ESgY,*ESgY2a;

int FluxAtLabRZ(double *gY,double Br[2],double Bz[2],int nr1, int nz1)
     /* gY_{ij}=gY[nr1*i+j] */
{
  int i,j,j1,k;
  double r,z;
  double R1,R2;
  double a,gt;
  double dr,dz;
  
  dr	=(Br[1]-Br[0])/(nr1-1);
  dz	=(Bz[1]-Bz[0])/(nz1-1);

  EcNjet	=0;
  j1	=0;
  for(i=0; i < nz1;  i++){
    j	=j1;
    j1	+=nr1;
    z	=Bz[0]+i*dz;
    if(ES1DMapZ2gt(&gt,&R2,&R1,z) == 0){
      r	=Br[0];
      while(r < R1 && j < j1){
	gY[j]	=EZc2gp*ESgY[ESNa]; /* EZout of plasma */
	j++;
	r	+=dr;
      }
      gt	=EZcgp-gt;
      a		=0.95;
      while(r < R2 && j < j1){
	k	=ES2DMapLab2Flx(&a,&gt,r,z);
	{
	  int ki,k;
	  double cs,sn,rc,rs;
	  ESSetSplA(a);
	  EcRjet[EcNjet]	=0.;
	  ki	=0;
	  cs	=cos(gt);
	  sn	=sin(gt);
	  splRA(&rc,NULL,ESsb,ESsb2a);
	  splRA(&rs,NULL,rsT,rsT2a);
	  EcZjet[EcNjet]	=ESaZ0[0]+rs+rc*sn;
	  for(k=1; k < ESFp1; k++){
	    ki	+=ESNa1;
	    splRA(&rc,NULL,rcT+ki,rcT2a+ki);
	    splRA(&rs,NULL,rsT+ki,rsT2a+ki);
	    cs	=cos(k*gt);
	    sn	=sin(k*gt);
	    EcRjet[EcNjet]	+=rc*cs+rs*sn;
	  }
	  splRA(&rc,NULL,rcT,rcT2a);
	  EcRjet[EcNjet]	=rc+2.*EcRjet[EcNjet]+ESaR0[0];
	  EcNjet++;
	}

	if(k){
	  printf("Mapping failed\n");
	  return(-1);
	}
	ESSetSplA(a);
	splRA(gY+j,NULL,ESgY,ESgY2a);
	gY[j]	=EZc2gp*(gY[j]-ESgY[ESNa]);
	r	+=dr;
	j++;
      }
      while(j < j1){
	gY[j]	=EZc2gp*ESgY[ESNa];
	r	+=dr;
	j++;
      }
    }
    else{
      while(j < j1){
	gY[j]	=EZc2gp*ESgY[ESNa]; /* EZout of plasma */
	j++;
      }
    }
  }
  return(0);
}

/**********************/

int ES1DMapZ2gtP(double *gt,double *R1,double *R2, double z)
{
  double t,cs,sn,s;
  int k,ki;

  ki =ESNa;
  t =(z-ESaZ0[0]-rsT[ki])/ESsb[ESNa];
  if(t > 1. || t < -1.){
    return(1);
  }
  *gt =asin(t);

  *R2 =ESaR0[0]+rcT[ki];
  *R1 =*R2;
  s =-1.;
  for(k=1; k < ESFp1; k++){
    t =*gt*k;
    cs =cos(t);
    sn =sin(t);
    ki +=ESNa1;
    *R2 +=2.*(rcT[ki]*cs+rsT[ki]*sn);
    *R1 +=2.*(rcT[ki]*cs-rsT[ki]*sn)*s;
    s =-s;
  }
  return(0);
}

int ESIsOutsidePlasma(double *lr1,double *lr2,double r,double z)
{
  double a,gt,t,cs,sn,s,R1,R2;
  int k,ki;

  r	-=ESaR0[0];
  z	-=ESaZ0[0];

  ki =ESNa;
  a =(z-rsT[ki])/ESsb[ESNa];
  if(a > 1.){
    *lr1	=ESsr[ESNp1*ESNa+ESNp/2];
    *lr2	=*lr1;
    return(2);
  }
  if(a < -1.){
    *lr1	=ESsr[ESNp1*ESNa+3*ESNp/2];
    *lr2	=*lr1;
    return(2);
  }
  gt =asin(a);
  R1 	=0.;
  R2 	=0.;
  s	=-1.;
  for(k=1; k < ESFp1; k++){
    t =gt*k;
    cs =cos(t);
    sn =sin(t);
    ki +=ESNa1;
    R1 +=(rcT[ki]*cs-rsT[ki]*sn)*s;
    R2 +=rcT[ki]*cs+rsT[ki]*sn;
    s =-s;
  }
  R1	=2.*R1+rcT[ESNa];
  R2	=2.*R2+rcT[ESNa];
  if(r < R1+1e-6 || R2-1e-6 < r){
    *lr1	=R1+ESaR0[0];
    *lr2	=R2+ESaR0[0];
    return(2);
  }
  if(ESa0 != 0.){
    a =z/ESsb[0];
    if(fabs(a) < 1.){
      gt =asin(a);
      R1 	=0.;
      R2 	=0.;
      ki 	=0;
      s		=-1.;
      for(k=1; k < ESFp1; k++){
	t =gt*k;
	cs =cos(t);
	sn =sin(t);
	ki +=ESNa1;
	R1 +=(rcT[ki]*cs-rsT[ki]*sn)*s;
	R2 +=rcT[ki]*cs+rsT[ki]*sn;
	s =-s;
      }
      R1	=2.*R1+rcT[0];
      R2	=2.*R2+rcT[0];
      *lr1	=R1+ESaR0[0];
      *lr2	=R2+ESaR0[0];
      if(R1 < r && r < R2){
	return(1);
      }
    }
    else{
      *lr1	=a > 1. ? ESsr[ESNp/2] : ESsr[3*ESNp/2];
      *lr2	=*lr1;
    }
  }
  return(0);
}

int ES2DMapLab2FlxP(double *A,double *t,double R, double Z)
{
  int Fl;
  double D,Rt;
  double cs,sn;
  double rc,rca,rcaa,rs,rsa,rsaa;
  int k,ki,i;

  double r,z,ra,dr2,EZdza,dz2,drp,dzp,a,gt;

  R -=ESaR0[0];
  Z -=ESaZ0[0];
  Fl =0;
  a =*A;
  gt =*t;

  ESSetSplA(a);
  r =0.;
  Rt =0.;
  ki =0;
  ra =Z < 0. ? -1. : 1.;
  for(k=1; k < ESFp1; k++){
    ki +=ESNa1;
    splRA(&rc,NULL,rcT+ki,rcT2a+ki);
    splRA(&rs,NULL,rsT+ki,rsT2a+ki);
    cs =cos(k*gt);
    sn =sin(k*gt);
    r +=(rc*cs+rs*sn);
    switch(k%4){
    case 0:
      Rt +=rc;
      break;
    case 1:
      Rt +=ra*rs;
      break;
    case 2:
      Rt -=rc;
      break;
    case 3:
      Rt -=ra*rs;
      break;
    }
  }
  splRA(&rc,NULL,rcT,rcT2a);
  splRA(&rs,NULL,rcT,rcT2a);
  Rt =rc+2.*Rt;
  r =rc+2.*r;
  if((R-Rt)*(r-Rt) < 0.) gt =EZcgp-gt;
  i =0;
  do{
    ESSetSplA(a);
    cs =cos(gt);
    sn =sin(gt);
    splRA2(&rc,&rca,&rcaa,ESsb,ESsb2a);
    splRA2(&rs,&rsa,&rsaa,rsT,rsT2a);
    z =rs+rc*sn-Z;
    EZdza =rsa+rca*sn;
    dz2 =rsaa+rcaa*sn;
    dzp =rc*cs;

    r =0.;
    ra =0.;
    dr2 =0.;
    drp =0.;
    ki =0;
    for(k=1; k < ESFp1; k++){
      ki +=ESNa1;
      splRA2(&rc,&rca,&rcaa,rcT+ki,rcT2a+ki);
      splRA2(&rs,&rsa,&rsaa,rsT+ki,rsT2a+ki);
      cs =cos(k*gt);
      sn =sin(k*gt);
      r  +=(rc*cs+rs*sn);
      ra +=(rca*cs+rsa*sn);
      dr2 +=(rcaa*cs+rsaa*sn);
      drp +=k*(-rc*sn+rs*cs);
    }
    splRA2(&rc,&rca,&rcaa,rcT,rcT2a);
    r =rc+2.*r-R;
    ra =rca+2.*ra;
    dr2 =rcaa+2.*dr2;
    drp *=2.;
    D =ra*dzp-drp*EZdza;
    drp =2.*(z*drp-r*dzp);
    drp /=sqrt(D*D+(dr2*dzp-drp*dz2)*drp)+D;
    a +=drp;
    ra =(r*EZdza-z*ra-EZcr2*(ra*dz2-dr2*EZdza)*drp*drp)/D;
    gt +=ra;
    if(a < 0.){
      Fl =2;
      a-=drp;
    }
    if(a > 1.){
      a	=1.;
    }
    i++;
  }while(i < 10 && fabs(ra)+fabs(drp) > 0.001);
  if(i < 10){
    *A =a;
    *t =gt;
    return(Fl);
  }
  return(1);
}

int FluxAtLabRZP(double *gY,double Br[2],double Bz[2],int nr1, int nz1)
     /* gY_{ij}=gY[nr1*i+j] */
{
  int i,j,j1,k;
  double r,z;
  double R1,R2;
  double a,gt;
  double dr,dz;

  printf("r1=%6g,  EZr2=%6g,  z1=%6g,  EZz2=%6g\n",Br[0],Br[1],Bz[0],Bz[1]);
  dr =(Br[1]-Br[0])/(nr1-1);
  dz =(Bz[1]-Bz[0])/(nz1-1);

  j1 =0;
  for(i=0; i < nz1;  i++){
    j =j1;
    j1 +=nr1;
    z =Bz[0]+i*dz;
    if(ES1DMapZ2gtP(&gt,&R1,&R2,z) == 0){
      r =Br[0];
      while(r < R1 && j < j1){
	gY[j] =0.; /* EZout of plasma */
        printf("A  %d  %10.3e  %10.3e    %10.3e\n",j,r,z,gY[j]);
	j++;
	r +=dr;
      }
      gt =EZcgp-gt;
      a =0.95;
      while(r < R2 && j < j1){
	k =ES2DMapLab2FlxP(&a,&gt,r,z);
	if(k){
	  printf("Mapping failed\n");
	  return(-1);
	}
	ESSetSplA(a);
	splRA(gY+j,NULL,ESgY,ESgY2a);
	gY[j] =EZc2gp*(gY[j]-ESgY[ESNa]);
        printf("B  %d  %10.3e  %10.3e    %10.3e\n",j,r,z,gY[j]);
	r +=dr;
	j++;
      }
      while(j < j1){
	gY[j] =0.;
        printf("C  %d  %10.3e  %10.3e    %10.3e\n",j,r,z,gY[j]);
	r +=dr;
	j++;
      }
    }
    else{
      while(j < j1){
	gY[j] =0.; /* EZout of plasma */
        printf("D  %d  %10.3e %10.3e\n",j,z,gY[j]);
	j++;
      }
    }
  }
  return(0);
}

int escpsi_(double *r1, double *EZr2, double *z1, double *EZz2)
{
  int i, j, j1, nr, nz;
  double Psi[625], Br[2], Bz[2];
  double r, z, dr, dz;

  printf("r1=%6g,  EZr2=%6g,  z1=%6g,  EZz2=%6g\n",*r1,*EZr2,*z1,*EZz2);
  nr = 10;  nz = 10;
  Br[0] = *r1;  Br[1] = *EZr2;  Bz[0] = *z1;  Bz[1] = *EZz2;        
  dr =(Br[1]-Br[0])/(nr-1);
  dz =(Bz[1]-Bz[0])/(nz-1);

  i = FluxAtLabRZP(&Psi[0], Br, Bz, nr, nz);
  i = 0;
  for (j=0; j < nr; j++){
    z =Bz[0]+j*dz;
    r = Br[0];
    for (j1=0; j1 < nz; j1++){
      printf("i=%d,  r=%6g,  z=%6g,  Psi=%g \n", i,r,z,Psi[i]);
      r +=dr;
      i++;
    }
  }
  return(0);
}

int ESrz2agt(double *la,double *lgt,double r,double z)
{
  int Fl;
  int k,ki,i;
  double a,gt,t,cs,sn,s,R1,R2;
  double x,y,EZz0,dr,dz,b,db;

  double D,Rt;
  double rc,rca,rcaa,rs,rsa,rsaa;
  double ra,dr2,EZdza,dz2,drp,dzp,da,dt;

  r	-=ESaR0[0];
  z	-=ESaZ0[0];

  ki =ESNa;
  a =(z-rsT[ki])/ESsb[ESNa];
  if(a > 1. || a < -1.){
    return(2);
  }
  gt =asin(a);
  R1 	=0.;
  R2 	=0.;
  s	=-1.;
  for(k=1; k < ESFp1; k++){
    t =gt*k;
    cs =cos(t);
    sn =sin(t);
    ki +=ESNa1;
    R1 +=(rcT[ki]*cs-rsT[ki]*sn)*s;
    R2 +=rcT[ki]*cs+rsT[ki]*sn;
    s =-s;
  }
  R1	=2.*R1+rcT[ESNa];
  R2	=2.*R2+rcT[ESNa];
  if(r < R1 || R2 < r){
    return(2);
  }
  if(ESa0 != 0.){
    x =z/ESsb[0];
    a	=fabs(x);
    if(a < 1.){
      gt 	=asin(x);
      R1 	=0.;
      R2 	=0.;
      ki 	=0;
      s		=-1.;
      for(k=1; k < ESFp1; k++){
	t =gt*k;
	cs =cos(t);
	sn =sin(t);
	ki +=ESNa1;
	R1 +=(rcT[ki]*cs-rsT[ki]*sn)*s;
	R2 +=rcT[ki]*cs+rsT[ki]*sn;
	s =-s;
      }
      R1	=2.*R1+rcT[0];
      R2	=2.*R2+rcT[0];
      if(R1 < r && r < R2){
	return(1);
      }
      x	=2.*r/(ESsr[0]-ESsr[ESNp/2]);
      y	=z/ESsb[0];
      a	=sqrt(y*y+x*x);	
      t 	=asin(y/a);
      if(r < R1){
	t	=EZcgp-t;
      }
      a	*=ESa0;
      if(a > ESsa[ESNa]){
	a	=ESsa[ESNa]*0.95;
      }
      if(a < ESa0){
	a	=ESa0*1.05;
      }
    }
    else{
      s	=2.*rcT1a[0]/ESsb1a[0];
      y	=z/ESsb[3];
      x	=r/(s*ESsb[3]);
      a	=sqrt(y*y+x*x);	
      if(a > 1.){
	k	=ESNp1*ESNa;
	x	=2.*r/(ESsr[k]-ESsr[k+ESNp/2]);
	y	=z/ESsb[ESNa];
	a	=sqrt(y*y+x*x);	
	t	=asin(y/a);
	a	*=ESsa[ESNa];
	if(a > ESsa[ESNa]){
	  a	=ESsa[ESNa]*0.95;
	}
      }
      else{
	t	=asin(y/a);
	a	*=ESsa[3];
      }
    }
  }
  
  gt	=t;
  ESSetSplA(a);
  dr	=0.;
  Rt	=0.;
  ki	=0;
  s	=z < 0. ? -1. : 1.;
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    splRA(&rc,NULL,rcT+ki,rcT2a+ki);
    splRA(&rs,NULL,rsT+ki,rsT2a+ki);
    cs	=cos(k*gt);
    sn	=sin(k*gt);
    dr	+=(rc*cs+rs*sn);
    switch(k%4){
    case 0:
      Rt	+=rc;
      break;
    case 1:
      Rt	+=s*rs;
      break;
    case 2:
      Rt	-=rc;
      break;
    case 3:
      Rt	-=s*rs;
      break;
    }
  }
  splRA(&rc,NULL,rcT,rcT2a);
  splRA(&rs,NULL,rcT,rcT2a);
  Rt	=rc+2.*Rt;
  dr	=rc+2.*dr;

  if((r-Rt)*(dr-Rt) < 0.){
    gt	=EZcgp-gt;
  }
  i	=0;
  do{
    ESSetSplA(a);
    cs	=cos(gt);
    sn	=sin(gt);
    splRA2(&rc,&rca,&rcaa,ESsb,ESsb2a);
    splRA2(&rs,&rsa,&rsaa,rsT,rsT2a);
    dz	=rs+rc*sn-z;
    EZdza	=rsa+rca*sn;
    dz2	=rsaa+rcaa*sn;
    dzp	=rc*cs;
    
    dr	=0.;
    ra	=0.;
    dr2	=0.;
    drp	=0.;
    ki	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA2(&rc,&rca,&rcaa,rcT+ki,rcT2a+ki);
      splRA2(&rs,&rsa,&rsaa,rsT+ki,rsT2a+ki);
      cs	=cos(k*gt);
      sn	=sin(k*gt);
      dr	+=(rc*cs+rs*sn);
      ra	+=(rca*cs+rsa*sn);
      dr2	+=(rcaa*cs+rsaa*sn);
      drp	+=k*(-rc*sn+rs*cs);
    }
    splRA2(&rc,&rca,&rcaa,rcT,rcT2a);
    dr	=rc+2.*dr-r;
    ra	=rca+2.*ra;
    dr2	=rcaa+2.*dr2;
    drp	*=2.;
    D	=ra*dzp-drp*EZdza;
    da	=2.*(dz*drp-dr*dzp);
    da	/=sqrt(D*D+(dr2*dzp-da*dz2)*da)+D;
    a	+=da;
    dt	=(dr*EZdza-dz*ra-EZcr2*(ra*dz2-dr2*EZdza)*da*da)/D;
    gt	+=dt;
    if(a < ESa0){
      if(ESa0 != 0.){
	a	=ESa0;
      }
      else{
	y	=z/ESsb1a[0];
	x	=r/(2.*rcT1a[0]);
	a	=sqrt(y*y+x*x);	
      }
    }
    if(a > 1.){
      a	=1;
    }
    i++;
  }while(i < 10 && fabs(dt)+fabs(da) > 0.001);
  if(i == 10){
    return(3);
  }
  *la	=a;
  *lgt	=gt;
  return(0);

#ifdef H
  a	=ESa0+(fabs(z)-ESsb[0])*(ESsa[ESNa]-ESa0)/(ESsb[ESNa]-ESsb[0]);
  x	=z > 0. ? a : -a;
  do{
    ESSetSplA(a);
    splRA(&EZz0,&dz,rsT,rsT2a);
    splRA(&b,&db,ESsb,ESsb2a);
    if(x < 0.){
      dz=-dz;
      b	=-b;
    }
    db	=(z-EZz0-b)/(dz+db);
    x	+=db;
    a	=fabs(x);
  }while(fabs(db) > 1e-5);
#endif
}

#ifndef stg_Mapping2R

int ESSetInsideInd(int *iA,double *R, double *Z,int nr,int nz)
{
  int i,j,k,k0,k1,n;
  double r,z,r0,EZz0,r1,z1,dr,dz;
  
  n	=nr*nz;
  dr	=(R[1]-R[0])/nr;
  dz	=(Z[1]-Z[0])/nz;
  z1	=Z[0];
  k	=0;
  for(i=0; i < nz; i++){
    EZz0	=z1;
    z1	=Z[0]+dz*(i+1);
    r1	=R[0];
    k0	=ESIsOutsidePlasma(&r,&z,r1,EZz0);
    k1	=ESIsOutsidePlasma(&r,&z,r1,z1);
    for(j=0; j < nr; j++){
      r0	=r1;
      r1	=R[0]+dr*(j+1);
      iA[k]	=k0 ? 0 : 1;
      if(k1 == 0) iA[k]	|=8;
      k0	=ESIsOutsidePlasma(&r,&z,r1,EZz0);
      if(k0 == 0) iA[k]	|=2;
      k1	=ESIsOutsidePlasma(&r,&z,r1,z1);
      if(k1 == 0) iA[k]	|=4;
      k++;
    }
  }
  return(0);
}

#endif
#endif

int ESCheckGSh(int i0,int i1,int m0,int m1,int Fl)
{
  int i,j,ji;
  int k,ki,kj;
#ifdef DEBUG
  double rc[ESFp1],rca[ESFp1],rc2a[ESFp1],rs[ESFp1],rsa[ESFp1],rs2a[ESFp1];
  double Kc[ESFp1],Kca[ESFp1],Ks[ESFp1],Ksa[ESFp1];
  double Nct[ESFp1],Nst[ESFp1];
  double Lc[ESFp1],Ls[ESFp1];
  double Vc[ESFp1],Vs[ESFp1];
  double LHSc[ES2Mp1],LHSs[ES2Mp1],RHSc[ES2Mp1],RHSs[ES2Mp1];
  double LHS[ESNp1*ESNa1],RHS[ESNp1*ESNa1];
  double KK[ESNp1],KKa[ESNp1];
#else
  double rc[33],rca[33],rc2a[33],rs[33],rsa[33],rs2a[33];
  double Kc[33],Kca[33],Ks[33],Ksa[33];
  double Nct[33],Nst[33];
  double Lc[33],Ls[33];
  double Vc[33],Vs[33];
  double LHSc[33],LHSs[33],RHSc[33],RHSs[33];
  double LHS[65*129],RHS[65*129];
  double KK[65],KKa[65];
#endif
  double b,ba,b2a;
  double r,r0,rr0,s,A,rA;
  double raa,rat,rtt,zaa,zat,ztt;
  double K,Ka,N,Nt,D,L,V;
  double cs,sn;
  double *lL,*lR;

  r0	=ESaR0[0];
  rr0	=1./r0;
  i	=0;
  b	=ESsb1a[i];
  ba	=ESsb1a[i];
  ki	=ESNa1;
  rc[1]	=2.*rcT1a[ki];
  rs[1]	=2.*rsT1a[ki];
  rca[1]=2.*rcT1a[ki];
  rsa[1]=2.*rsT1a[ki];
  D	=b*rc[1];
  r	=r0;
  L	=D;
  s	=1./(r*D);
  for(j=0; j < ESNp; j++){
    cs	=EScs1[j];
    sn	=ESsn1[j];
    EZdrgt[j]	=-rc[1]*sn+rs[1]*cs;
    EZdzgt[j]	=b*cs;
    EZdra[j]	=rca[1]*cs+rsa[1]*sn;
    EZdza[j]	=ba*sn;
    rat		=-rca[1]*sn+rsa[1]*cs;
    rtt		=-rc[1]*cs-rs[1]*sn;
    zat		=ba*cs;
    ztt		=-b*sn;
    K	=(EZdrgt[j]*EZdrgt[j]+EZdzgt[j]*EZdzgt[j])*s;
    N	=(EZdra[j]*EZdrgt[j]+EZdza[j]*EZdzgt[j])*s;
    Nt	=(rat*EZdrgt[j]+EZdra[j]*rtt+zat*EZdzgt[j]+EZdza[j]*ztt)*s;
    LHS[j]	=(2.*K-Nt)*ESdgY[0];
    RHS[j]	=-L*ESjs[0];
    LHS[j]	-=RHS[j];
  }
  LHS[j]	=LHS[0];
  RHS[j]	=RHS[0];
  ESP2F(LHSc,LHSs,LHS,ES2Mp);
  ESP2F(RHSc,RHSs,RHS,ES2Mp);
  ji	=ESNp1;
  for(i=1; i < ESNa1; i++){
    A	=ESsa[i];
    rA	=1./A;
    b	=ESsb[i];
    ba	=ESsb1a[i];
    b2a	=ESsb2a[i];
    ki	=i;
    rc[0]	=rcT[ki];
    rs[0]	=rsT[ki];
    rca[0]	=rcT1a[ki];
    rsa[0]	=rsT1a[ki];
    rc2a[0]	=rcT2a[ki];
    rs2a[0]	=rsT2a[ki];

    Kc[0]	=ESg22c[ki];
    Ks[0]	=ESg22s[ki];
    Kca[0]	=ESg22c1[ki];
    Ksa[0]	=ESg22s1[ki];

    Nct[0]	=0.;
    Nst[0]	=0.;
    Lc[0]	=ESLc[ki];
    Ls[0]	=ESLs[ki];
    Vc[0]	=ESVc[ki];
    Vs[0]	=ESVs[ki];

    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=2.*rcT[ki];
      rs[k]	=2.*rsT[ki];
      rca[k]	=2.*rcT1a[ki];
      rsa[k]	=2.*rsT1a[ki];
      rc2a[k]	=2.*rcT2a[ki];
      rs2a[k]	=2.*rsT2a[ki];

      if(k < ESMp1){
	Kc[k]	=2.*ESg22c[ki];
	Ks[k]	=2.*ESg22s[ki];
	Kca[k]	=2.*ESg22c1[ki];
	Ksa[k]	=2.*ESg22s1[ki];
	
	Nct[k]	=2.*ESg22s[ki];
	Nst[k]	=-2.*ESg22c[ki];
	Lc[k]	=2.*ESLc[ki];
	Ls[k]	=2.*ESLs[ki];
	Vc[k]	=2.*ESVc[ki];
	Vs[k]	=2.*ESVs[ki];
      }
    }
    for(j=0; j < ESNp; j++){
      EZdra[j]	=rca[0];
      raa	=rc2a[0];
      rat	=0.;
      rtt	=0.;
      EZdrgt[j]	=0.;
      K		=Kc[0];
      Ka	=Kca[0];
      Nt	=0.;
      L		=Lc[0];
      V		=Vc[0];
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	cs	=EScs1[kj];
	sn	=ESsn1[kj];
	EZdra[j]	+=rca[k]*cs+rsa[k]*sn;
	raa	+=rc2a[k]*cs+rs2a[k]*sn;
	EZdrgt[j]	+=(-rc[k]*sn+rs[k]*cs)*k;
	rat	+=(-rca[k]*sn+rsa[k]*cs)*k;
	rtt	-=(rc[k]*cs+rs[k]*sn)*k*k;
	if(k < ESMp1){
	  K	+=Kc[k]*cs+Ks[k]*sn;
	  Ka	+=Kca[k]*cs+Ksa[k]*sn;
	}
      }
      cs	=EScs1[j];
      sn	=ESsn1[j];
      EZdzgt[j]	=b*cs;
      EZdza[j]	=rsa[0]+ba*sn;
      zaa	=rs2a[0]+b2a*sn;
      zat	=ba*cs;
      ztt	=-b*sn;
      D		=EZdra[j]*EZdzgt[j]-EZdrgt[j]*EZdza[j];
      r		=ESsr[ji];
      L		=D*r0*rA/r;
      V		=D*r*rr0*rA-L;
      s		=1./(r*D);
      N		=(EZdra[j]*EZdrgt[j]+EZdza[j]*EZdzgt[j])*s;
      Nt	=(rat*EZdrgt[j]+EZdra[j]*rtt+zat*EZdzgt[j]+EZdza[j]*ztt)*s
	-N*(EZdrgt[j]/r+(rat*EZdzgt[j]+EZdra[j]*ztt-rtt*EZdza[j]-EZdrgt[j]*zat)/D);
      LHS[ji]	=(2.*K+A*Ka-Nt)*ESdgY[i]+A*K*ESdgY1a[i];

#ifdef H
      if(i == ESNa){
	K	=(EZdrgt[j]*EZdrgt[j]+EZdzgt[j]*EZdzgt[j])*s*rA;
	Ka	=2.*(EZdrgt[j]*rat+EZdzgt[j]*zat)*s
	  -A*K*(EZdra[j]/r+(raa*EZdzgt[j]+EZdra[j]*zat-rat*EZdza[j]-EZdrgt[j]*zaa)/D);
	KK[j]	=V;
	KKa[j]	=(raa*EZdzgt[j]+EZdra[j]*zat-rat*EZdza[j]-EZdrgt[j]*zaa)*(r*rr0-r0/r)
	  +D*EZdra[j]*(rr0+r0/(r*r));
      }

      LHS[ji]	=(K+Ka-Nt)*ESdgY[i]+A*K*ESdgY1a[i];
#endif

      RHS[ji]	=-L*ESjs[i]-V*ESjp[i];
      LHS[ji]	-=RHS[ji];
      ji++;
    }
    j	=ji-ESNp;
    LHS[ji]	=LHS[j];
    RHS[ji]	=RHS[j];

    ji++;
    if(m0 < m1){
      lL	=LHS+j;
      lR	=RHS+j;
      ESP2F(LHSc,LHSs,lL,m1-1);
      ESP2F(RHSc,RHSs,lR,m1-1);
      for(j=0; j < ESNp; j++){
	r	=0.;
	s	=0.;
	kj	=0;
	for(k=m0 ? m0 : 1; k < m1; k++){
	  kj	+=j;
	  if(kj >= ESNp){
	    kj	-=ESNp;
	  }
	  cs	=EScs1[kj];
	  sn	=ESsn1[kj];
	  s	+=LHSc[k]*cs+LHSs[k]*sn;
	  r	+=RHSc[k]*cs+RHSs[k]*sn;
	}
	s	*=2.;	
	r	*=2.;
	if(m0 == 0){
	  s	+=LHSc[0];
	  r	+=RHSc[0];
	}
	lL[j]	=s;
	lR[j]	=r;
      }
      lL[j]	=lL[0];
      lR[j]	=lR[0];
    }
  }
#ifdef DEBUG
  {
    extern int ESEqSolvFl;
    double x[2],y[ESNp1],y0[ESNp1];
    char ln[32];

    j	=ESNp1*i0;
    j	=0;
    r	=log10(fabs(RHS[j])+1e-20);
    y[0]	=r;
    y[1]	=y[0];
    x[0]	=ESgt[0];
    x[1]	=ESgt[ESNp];
    ji	=ESNp1*i1;
    ji	=ESNp1*ESNa1;
    for(; j < ji; j++){
      r	=log10(fabs(RHS[j])+1e-20);
      if(y[0] > r){
	y[0]	=r;
      }
      if(y[1] < r){
	y[1]	=r;
      }
      r	=log10(fabs(LHS[j])+1e-20);
      if(y[0] > r){
	y[0]	=r;
      }
      if(y[1] < r){
	y[1]	=r;
      }
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
    sprintf(ln,"%d<=i<%d m=%d (%s) %s",i0,i1,ESMp
	    ,(ESEqSolvFl&0x0F) ? "RungK" : "Lsode",ShotNm);
    SetPlotName("sqrt(gF)","GSh log10|LHS-RHS|",ln);
    Scale2d(4,x,y,2,6,2);
    if(Fl){
      PSNewFrame();
    }
    
    if(i0 < 0){
      i0	=0;
    }
    if(i1 > ESNa1){
      i1	=ESNa1;
    }
    ji	=ESNp1*i0;
    for(i=i0; i < i1; i++){
      lL	=LHS+ji;
      lR	=RHS+ji;
      for(j=0; j < ESNp1; j++){
	y[j]	=log10(fabs(lL[j])+1e-20);
	y0[j]	=log10(fabs(lR[j])+1e-20);
      }
      Plot2d(4,ESgt,y,ESNp1,6,0,14,0);
      Plot2d(4,ESgt,y0,ESNp1,6,0,4,0);
      if(Fl){
	PSPlot2d(ESgt,y,ESNp1,6,0,14,0);
	PSPlot2d(ESgt,y0,ESNp1,6,0,4,0);
      }
      ji	+=ESNp1;
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

int ESErrorInGSh()
{
  static int Fl=0;
  int i,j,ji;
  int k,ki,kj;
#ifdef DEBUG
  double rc[ESFp1],rca[ESFp1],rc2a[ESFp1],rs[ESFp1],rsa[ESFp1],rs2a[ESFp1];
#else
  double rc[33],rca[33],rc2a[33],rs[33],rsa[33],rs2a[33];
#endif
  double b,ba,b2a;
  double r,r0,rr0,s,A,rA;
  double ra,rt,za,zt,raa,rat,rtt,zaa,zat,ztt;
  double LHS,RHS,LHS0,RHS0;
  double K,Ka,N,Nt,D,L,V;
  double cs,sn;
  double Err,Max;
  int i0,j0;
  int fpa;
  int offset,Off0;
  static char RunIDa[16];
  char ln[64],*lc[24];
  extern double ESTime;
  extern int ESEqSolvRadC,ESEqSolvInPr,ESEqSolvInCr,ESnDtPr,ESnDtCr,ESFail;

  static double t=0.;
  ESTime	=t;

  offset	=0;
  if(Fl == 0){
    extern char ESRunID[];
#ifdef DEBUG
    sprintf(RunIDa,"%s.err",ShotNm);
#else      
    k	=CbStrCpy(RunIDa,ESRunID);
    strcpy(RunIDa+k,".err");
#endif
    fpa	=open(RunIDa,O_RDWR|O_CREAT|O_TRUNC,S_IRWXU|S_IROTH|S_IXOTH);
    lc[0]	="      subroutine readescerr ! an example reading the file\n";
    lc[1]	="      real*8 t,err\n";
    lc[2]	="      real*8 LHS(65,21),RHS(65,21)\n";
    lc[3]	="      integer na1,nt1\n";
    lc[4]	="      integer i,j,k,u\n";
    lc[5]	="      character*64 ch\n";
    lc[6]	="      u=1\n";
    sprintf(ln,"      open(unit=u,file='%s')\n",RunIDa);
    lc[7]	=ln;
    lc[8]	="      do k=1,26\n";
    lc[9]	="        read(u,*)\n";
    lc[10]	="      enddo\n\n";
    lc[11]	="      read(u,'(a6,i3,a7,i3)')ch,na1,ch,nt1\n";
    lc[12]	="      do k=1,1000\n";
    lc[13]	="        read(u,'(a5,e12.5,a45,e9.2)',end=1)ch,t,ch,err\n";
    lc[14]	="        read(u,*)\n";
    lc[15]	="        do i=1,na1\n";
    lc[16]	="          do j=1,nt1\n";
    lc[17]	="          read(u,'(a5,2e12.4,e10.2)')ch,LHS(j,i),RHS(j,i)\n";
    lc[18]	="          enddo\n";
    lc[19]	="        enddo\n";
    lc[20]	="      enddo\n";
    lc[21]	=" 1    continue\n";
    lc[22]	="      close(u)\n";
    lc[23]	="      end\n\n";
    Fl	=1;
  }
  else{
    fpa	=open(RunIDa,O_RDWR);
    lseek(fpa,(off_t)0,SEEK_END);
    Fl	=2;
  }
  if(fpa == -1){
    printf("%s - ASCII Error file cannot be open\n",RunIDa);
    return(1);
  }
  if(Fl == 1){
    for(i=0; i < 24; i++){
      k	=strlen(lc[i]);
      write(fpa,lc[i],k);
      offset	+=k;
    }
    sprintf(ln,"ESNa1=%3d ESNp1=%3d\n",ESNa1,ESNp1);
    write(fpa,ln,20);
    offset	+=20;
  }
  Err	=0.;
  Max	=0.;
  sprintf(ln,"time=%12.5e RadC=%d InPr=%d InCr=%d nPr=%d nCr=%d 0x%8.8x"
	  ,ESTime,ESEqSolvRadC,ESEqSolvInPr,ESEqSolvInCr,ESnDtPr,ESnDtCr
	  ,ESFail);
  write(fpa,ln,57);
  offset	+=57;
  Off0	=offset;
  sprintf(ln," Err=%9.2e\n",Err);
  write(fpa,ln,15);
  offset+=15;
  sprintf(ln,"%2s %2s %11s %11s %9s\n","i","j","LHS","RHS","LHS-RHS");
  write(fpa,ln,40);
  offset	+=40;

  r0	=ESaR0[0];
  rr0	=1./r0;
  i	=0;
  b	=ESsb1a[i];
  ba	=ESsb1a[i];
  ki	=ESNa1;
  rc[1]	=2.*rcT1a[ki];
  rs[1]	=2.*rsT1a[ki];
  rca[1]=2.*rcT1a[ki];
  rsa[1]=2.*rsT1a[ki];
  D	=b*rc[1];
  r	=r0;
  L	=D;
  s	=1./(r*D);

  j	=0;
  cs	=1.;
  rt	=rs[1];
  zt	=b*cs;
  ra	=rca[1];
  za	=0.;
  rat	=rsa[1];
  rtt	=-rc[1];
  zat	=ba;
  ztt	=0.;
  K	=(rt*rt+zt*zt)*s;
  N	=(ra*rt+za*zt)*s;
  Nt	=(rat*rt+ra*rtt+zat*zt+za*ztt)*s;
  LHS	=(2.*K-Nt)*ESdgY[0];
  RHS	=-L*ESjs[0];
  LHS0	=LHS;
  RHS0	=RHS;
  if(Max < fabs(RHS)){
    Max	=fabs(RHS);
  }
  if(Err < fabs(LHS-RHS)){
    i0	=i;
    j0	=j;
    Err	=fabs(LHS-RHS);
  }
  sprintf(ln,"%2d %2d %11.4e %11.4e %9.2e\n",i,j,LHS,RHS,LHS-RHS);
  write(fpa,ln,40);
  offset	+=40;
  for(j=1; j < ESNp; j++){
    cs	=EScs1[j];
    sn	=ESsn1[j];
    rt	=-rc[1]*sn+rs[1]*cs;
    zt	=b*cs;
    ra	=rca[1]*cs+rsa[1]*sn;
    za	=ba*sn;
    rat	=-rca[1]*sn+rsa[1]*cs;
    rtt	=-rc[1]*cs-rs[1]*sn;
    zat	=ba*cs;
    ztt	=-b*sn;
    K	=(rt*rt+zt*zt)*s;
    N	=(ra*rt+za*zt)*s;
    Nt	=(rat*rt+ra*rtt+zat*zt+za*ztt)*s;
    LHS	=(2.*K-Nt)*ESdgY[0];
    RHS	=-L*ESjs[0];
    sprintf(ln,"%2d %2d %11.4e %11.4e %9.2e\n",i,j,LHS,RHS,LHS-RHS);
    write(fpa,ln,40);
    offset	+=40;
    if(Max < fabs(RHS)){
      Max	=fabs(RHS);
    }
    if(Err < fabs(LHS-RHS)){
      i0	=i;
      j0	=j;
      Err	=fabs(LHS-RHS);
    }
  }
  LHS	=LHS0;
  RHS	=RHS0;
  sprintf(ln,"%2d %2d %11.4e %11.4e %9.2e\n",i,j,LHS,RHS,LHS-RHS);
  write(fpa,ln,40);
  offset	+=40;
  ji	=ESNp1;
  for(i=1; i < ESNa1; i++){
    A	=ESsa[i];
    rA	=1./A;
    b	=ESsb[i];
    ba	=ESsb1a[i];
    b2a	=ESsb2a[i];
    ki	=i;
    rc[0]	=rcT[ki];
    rs[0]	=rsT[ki];
    rca[0]	=rcT1a[ki];
    rsa[0]	=rsT1a[ki];
    rc2a[0]	=rcT2a[ki];
    rs2a[0]	=rsT2a[ki];
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=2.*rcT[ki];
      rs[k]	=2.*rsT[ki];
      rca[k]	=2.*rcT1a[ki];
      rsa[k]	=2.*rsT1a[ki];
      rc2a[k]	=2.*rcT2a[ki];
      rs2a[k]	=2.*rsT2a[ki];
    }
    j	=0;
    ra	=rca[0];
    raa	=rc2a[0];
    rat	=0.;
    rtt	=0.;
    rt	=0.;
    kj	=0;
    for(k=1; k < ESFp1; k++){
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      ra	+=rca[k];
      raa	+=rc2a[k];
      rt	+=rs[k]*k;
      rat	+=rsa[k]*k;
      rtt	-=rc[k]*k*k;
    }
    zt	=b;
    za	=rsa[0];
    zaa	=rs2a[0];
    zat	=ba;
    ztt	=0.;
    D	=ra*zt-rt*za;
    r	=ESsr[ji];
    L	=D*r0*rA/r;
    V	=D*r*rr0*rA-L;
    s	=1./(r*D);
    N	=(ra*rt+za*zt)*s;
    Nt	=(rat*rt+ra*rtt+zat*zt+za*ztt)*s
      -N*(rt/r+(rat*zt+ra*ztt-rtt*za-rt*zat)/D);
    K	=(rt*rt+zt*zt)*s*rA;
    Ka	=2.*(rt*rat+zt*zat)*s
      -A*K*(ra/r+(raa*zt+ra*zat-rat*za-rt*zaa)/D);
    LHS	=(K+Ka-Nt)*ESdgY[i]+A*K*ESdgY1a[i];
    RHS	=-L*ESjs[i]-V*ESjp[i];
    LHS0=LHS;
    RHS0=RHS;
    sprintf(ln,"%2d %2d %11.4e %11.4e %9.2e\n",i,j,LHS,RHS,LHS-RHS);
    write(fpa,ln,40);
    offset	+=40;
    if(Max < fabs(RHS)){
      Max	=fabs(RHS);
    }
    if(Err < fabs(LHS-RHS)){
      i0	=i;
      j0	=j;
      Err	=fabs(LHS-RHS);
    }
    ji++;
    for(j=1; j < ESNp; j++){
      ra	=rca[0];
      raa	=rc2a[0];
      rat	=0.;
      rtt	=0.;
      rt	=0.;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	cs	=EScs1[kj];
	sn	=ESsn1[kj];
	ra	+=rca[k]*cs+rsa[k]*sn;
	raa	+=rc2a[k]*cs+rs2a[k]*sn;
	rt	+=(-rc[k]*sn+rs[k]*cs)*k;
	rat	+=(-rca[k]*sn+rsa[k]*cs)*k;
	rtt	-=(rc[k]*cs+rs[k]*sn)*k*k;
      }
      cs	=EScs1[j];
      sn	=ESsn1[j];
      zt	=b*cs;
      za	=rsa[0]+ba*sn;
      zaa	=rs2a[0]+b2a*sn;
      zat	=ba*cs;
      ztt	=-b*sn;
      D		=ra*zt-rt*za;
      r		=ESsr[ji];
      L		=D*r0*rA/r;
      V		=D*r*rr0*rA-L;
      s		=1./(r*D);
      N		=(ra*rt+za*zt)*s;
      Nt	=(rat*rt+ra*rtt+zat*zt+za*ztt)*s
	-N*(rt/r+(rat*zt+ra*ztt-rtt*za-rt*zat)/D);
      K	=(rt*rt+zt*zt)*s*rA;
      Ka	=2.*(rt*rat+zt*zat)*s
	-A*K*(ra/r+(raa*zt+ra  *zat-rat*za-rt*zaa)/D);
      LHS	=(K+Ka-Nt)*ESdgY[i]+A*K*ESdgY1a[i];
      RHS	=-L*ESjs[i]-V*ESjp[i];
      sprintf(ln,"%2d %2d %11.4e %11.4e %9.2e\n",i,j,LHS,RHS,LHS-RHS);
      write(fpa,ln,40);
      offset	+=40;
      if(Max < fabs(RHS)){
	Max	=fabs(RHS);
      }
      if(Err < fabs(LHS-RHS)){
	i0	=i;
	j0	=j;
	Err	=fabs(LHS-RHS);
      }
      ji++;
    }
    LHS	=LHS0;
    RHS	=RHS0;
    sprintf(ln,"%2d %2d %11.4e %11.4e %9.2e\n",i,j,LHS,RHS,LHS-RHS);
    write(fpa,ln,40);
    offset	+=40;
    ji++;
  }
  offset	=Off0-offset;
  lseek(fpa,(off_t)offset,SEEK_END);
  sprintf(ln," Err=%9.2e\n",Err/Max);
  printf(" Err=%9.2e %9.2e %9.2e\n",Err,Max,Err/Max);
  write(fpa,ln,15);
  close(fpa);
  t	+=1.;
  return(0);
}

int ESFullClean()
{
  int i,k,j,ji,ki;
  
  ESaR0[0]	=1.5;
  ESaZ0[0]	=0.;
  ki	=0;
  for(k=0; k < ESFp1; k++){
    for(i=0; i < ESNa1; i++){
      rcT[ki]	=0.;
      rcT1a[ki]	=0.;
      rcT2a[ki]	=0.;
      rsT[ki]	=0.;
      rsT1a[ki]	=0.;
      rsT2a[ki]	=0.;
      ki++;
    }
  }

  ji	=0;
  ki	=ESNa1;
  for(i=0; i < ESNa1; i++){
    ESsb[i]	=ESsa[i];
    ESsb1a[i]	=1.;
    ESsb2a[i]	=0.;
    rcT[ki]	=0.25*ESsa[i];
    rcT1a[ki]	=0.25;
    for(j=0; j < ESNp1; j++){
      ESsr[ji]	=ESaR0[0]+2.*rcT[ki]*EScs1[j];
      ESsz[ji]	=ESaZ0[0]+ESsb[i]*ESsn1[j];
      ji++;
    }
    ki++;
  }
  return(0);
}


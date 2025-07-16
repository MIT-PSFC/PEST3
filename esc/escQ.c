#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef mzl_MLIB
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern unsigned int ESmemlev;
#endif

extern int ESNa,ESNa1,ESNp,ESNp1,ESnAP,ES2Mp1,ESFp1;
extern double *ESsa,*ESpa,*ESaR0,*ESaZ0;
extern double *ESsb,*ESsb1a,*ESsb2a;
extern double *ESgt,*EScs1,*ESsn1;
extern double *ESjp,*ESjp1a,*ESsp,*ESsp1a,*ESsp2a,*ESsq,*ESsq1a,*ESsq2a,*ESPs,*ESPs1a,*ESPs2a;
extern double *ESjs,*ESjb;
extern double *ESFF,*ESFF1a,*ESFF2a,*ESaF,*ESaF1a,*ESaF2a;
extern double *ESgF,*ESdgY,*ESdgY2a,*ESaJ,*ESaJ2a;
extern double *ESPs,*ESPs1a,*ESPs2a,*ESaT,*ESaT1a,*ESaT2a;
extern double *rcT,*rcT1a,*rcT2a,*rsT,*rsT1a,*rsT2a;
extern double *ESLc,*ESLc1,*ESLc2,*ESLs,*ESLs1,*ESVc,*ESVc1,*ESVc2,*ESVs,*ESVs1;
extern double *ESg22c,*ESg22c1,*ESg22s,*ESg22s1;
extern double *ESdgY,*ESdgY1a;
extern double ESBt,ESRBt,ESRext;

extern double *ESsr,*ESsra,*ESsrt,*ESsrat;
extern double *ESsz,*ESsza,*ESszt,*ESszat;
extern double *ESaB,*ESaBa,*ESaBt,*ESaBat;
extern double *ESgH,*ESgHa,*ESgHt,*ESgHat;

static int Spi0,Spi1,Spj0,Spj1;
static double ha,rha,hq,rhq;
static double gmP;

static double cmi=1.6726e-24; /* [g] - proton mass */
static double cc=2.9979e+10; /* [cm/sec] - speed of light */
static double ce=4.8032e-10; /* [CGS] - proton electric charge */
static double cEeV=1.6022e-12; /* energy of 1 eV */
static double cEkeV=1.6022e-9; /* energy of 1 keV */
double cgWc=0.9579e+8; /* [1/sec] - proton cyclotron frequency 
				 in 1 T field */
static double cV1keV;		/* velocity for the energy 1 keV */
static double cgr1keV;		/* ion larmor radius for 1 keV in 1 T */

static int neq=4;
static double hcur,xout,reltol=1.0e-6,abstol=1e-5;
static double *rwork=NULL,*yh,*vatol;
static int itol=2,itask,istate=1,iopt=0,mf=10,lrw,liw=20;
static int *iwork=NULL;
static double yGC[4],dyGC[4];

static int NWRK=0,nWRK;
static double *RKgr0,*RKgr1,*RKgr2,*RKgr3,*RKgr4;
static double *RKa0,*RKa1,*RKa2,*RKa3,*RKa4;
static double *RKgq0,*RKgq1,*RKgq2,*RKgq3,*RKgq4;
static double *RKgf0,*RKgf1,*RKgf2,*RKgf3,*RKgf4;
static double *rPrt0,*rPrt1,*zPrt0,*zPrt1;
 
static double *Xgr,*Xa,*Xgq,*Xgf,*Xgm,*dXgr,*dXa,*dXgq,*dXgf;

int GCInitRKSolv(int n)
{
  nWRK	=(20+4+9)*n;
  if(NWRK != nWRK){
    if(NWRK){
      free(RKgr0);
    }
    RKgr0	=(double*)malloc(nWRK*sizeof(double));
    if(RKgr0 == NULL){
      printf("No memory for RKgr0 in GCInitRKSolv(%d)\n",n);
      return(1);
    }
    RKgr1	=RKgr0	+n;
    RKgr2	=RKgr1	+n;
    RKgr3	=RKgr2	+n;
    RKgr4	=RKgr3	+n;

    RKa0	=RKgr4	+n;
    RKa1	=RKa0	+n;
    RKa2	=RKa1	+n;
    RKa3	=RKa2	+n;
    RKa4	=RKa3	+n;

    RKgq0	=RKa4	+n;
    RKgq1	=RKgq0	+n;
    RKgq2	=RKgq1	+n;
    RKgq3	=RKgq2	+n;
    RKgq4	=RKgq3	+n;

    RKgf0	=RKgq4	+n;
    RKgf1	=RKgf0	+n;
    RKgf2	=RKgf1	+n;
    RKgf3	=RKgf2	+n;
    RKgf4	=RKgf3	+n;

    rPrt0	=RKgf4	+n;  
    rPrt1	=rPrt0	+n;  
    zPrt0	=rPrt1	+n;  
    zPrt1	=zPrt0	+n;  
    Xgr	=zPrt1	+n;
    Xa	=Xgr	+n;
    Xgq	=Xa	+n;
    Xgf	=Xgq	+n;
    Xgm	=Xgf	+n;
    dXgr=Xgm	+n;
    dXa	=dXgr	+n;
    dXgq=dXa	+n;
    dXgf=dXgq	+n;

    NWRK	=nWRK;
  }
  return(0);
}

int GCDeInitRKSolv()
{
  if(NWRK){
    free(RKgr0);
    NWRK	=0;
  }
  return(0);
}

static int NPart=0,nPart=10;
static double XF[1024],XFa[1024],XgFa[1024],XgFaa[1024]
,XgYa[1024],XgYaa[1024],XT[1024],XTa[1024],XP[1024],XPa[1024]
,Xr[1024],Xra[1024],Xrq[1024],Xz[1024],Xza[1024],Xzq[1024]
,XB[1024],XBa[1024],XBq[1024],Xgh[1024],Xgha[1024],Xghq[1024];

#ifdef XWIN
int ESIget2DFunc(double *r,double *ra,double *rq
		 ,double *z,double *za,double *zq
		 ,double *B,double *Ba,double *Bq
		 ,double *gh,double *gha,double *ghq
		 ,double a,double gq)
{
  int i;
  i	=1;
  esiget2dfunctions_(&a,&gq,&i);
  *r	=Xr[0];
  *ra	=Xra[0];
  *rq	=Xrq[0];
  *z	=Xz[0];
  *za	=Xza[0];
  *zq	=Xzq[0];
  *B	=XB[0];
  *Ba	=XBa[0];
  *Bq	=XBq[0];
  *gh	=Xgh[0];
  *gha	=Xgha[0];
  *ghq	=Xghq[0];
  return(0);
}    

int GCRKStep(double *X,double h)
{
  double x,x0;
  int i,j;

  x0	=*X;
  for(i=0; i < nPart; i++){
    RKgr1[i]	=h*dXgr[i];
    RKgr0[i]	=Xgr[i]+EZcr2*RKgr1[i];
    RKa1[i]	=h*dXa[i];
    RKa0[i]	=Xa[i]+EZcr2*RKa1[i];
    RKgq1[i]	=h*dXgq[i];
    RKgq0[i]	=Xgq[i]+EZcr2*RKgq1[i];
    RKgf1[i]	=h*dXgf[i];
    RKgf0[i]	=Xgf[i]+EZcr2*RKgf1[i];
  }
  x	=(*X)+EZcr2*h;
  gcmotion_(RKgr2,RKa2,RKgq2,RKgf2,RKgr0,RKa0,RKgq0,Xgm,&nPart);
  for(i=0; i < nPart; i++){
    RKgr2[i]	*=h;
    RKa2[i]	*=h;
    RKgq2[i]	*=h;
    RKgf2[i]	*=h;
    RKgr0[i]	=Xgr[i]+EZcr2*RKgr2[i];
    RKa0[i]	=Xa[i]+EZcr2*RKa2[i];
    RKgq0[i]	=Xgq[i]+EZcr2*RKgq2[i];
    RKgf0[i]	=Xgf[i]+EZcr2*RKgf2[i];
  }
  gcmotion_(RKgr3,RKa3,RKgq3,RKgf3,RKgr0,RKa0,RKgq0,Xgm,&nPart);
  for(i=0; i < nPart; i++){
    RKgr3[i]	*=h;
    RKa3[i]	*=h;
    RKgq3[i]	*=h;
    RKgf3[i]	*=h;
    RKgr0[i]	=Xgr[i]+RKgr3[i];
    RKa0[i]	=Xa[i]+RKa3[i];
    RKgq0[i]	=Xgq[i]+RKgq3[i];
    RKgf0[i]	=Xgf[i]+RKgf3[i];
  }
  (*X)		+=h;
  gcmotion_(RKgr4,RKa4,RKgq4,RKgf4,RKgr0,RKa0,RKgq0,Xgm,&nPart);
  for(i=0; i < nPart; i++){
    RKgr4[i]	*=h;
    RKa4[i]	*=h;
    RKgq4[i]	*=h;
    RKgf4[i]	*=h;
    Xgr[i]	+=EZcr6*(RKgr1[i]+RKgr4[i])+EZcr3*(RKgr2[i]+RKgr3[i]);
    Xa[i]	+=EZcr6*(RKa1[i]+RKa4[i])+EZcr3*(RKa2[i]+RKa3[i]);
    Xgq[i]	+=EZcr6*(RKgq1[i]+RKgq4[i])+EZcr3*(RKgq2[i]+RKgq3[i]);
    Xgf[i]	+=EZcr6*(RKgf1[i]+RKgf4[i])+EZcr3*(RKgf2[i]+RKgf3[i]);
  }
  gcmotion_(dXgr,dXa,dXgq,dXgf,Xgr,Xa,Xgq,Xgm,&nPart);
  return(0);
}
#endif

int ESBasicFunctions()
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
  for(j=0; j < ESNp1; j++){
    cs	=EScs1[j];
    sn	=ESsn1[j];
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

    ESsra[j]	=ra;
    ESsrt[j]	=0.;
    ESsrat[j]	=rat;
    ESsza[j]	=za;
    ESszt[j]	=0.;
    ESszat[j]	=zat;

    Ga	=-ra*G*G;
    Gt	=0.;
    Gat	=-rat*G*G;
    ESaB[j]	=pF*G;
    ESaBa[j]	=pF*Ga;
    ESaBt[j]	=0.;
    ESaBat[j]	=pF*Gat;

    ESgH[j]	=0.;
    ESgHt[j]	=0.;

    D	=ra*zat-za*rat;
#ifdef H
    ra	=ra+a*raa;
    rt	=a*(rat+a*EZcr2*raat);

    za	=za+a*zaa;
    zt	=a*(zat+a*EZcr2*zaat);

    1/r	=1-a*ra*G;

    D=(ra+a*raa)*(zat+a*EZcr2*zaat)-(za+a*zaa)*(rat+a*EZcr2*raat);
    D/r=1+a*(raa*zat+a*EZcr2*ra*zaat-a*zaa*rat-a*EZcr2*za*raat)/D0-a*ra*G;

    gH	=a*(raa*zat-zaa*rat+EZcr2*(ra*zaat-za*raat))/D0-a*ra*G;
#endif
    ESgHa[j]	=(raa*zat-zaa*rat+EZcr2*(ra*zaat-za*raat))/D-ra*G;
    ESgHat[j]	=(raa*zatt-zaa*ratt+EZcr2*(ra*zaatt+zat*raat-rat*zaat-za*raatt))/D-rat*G;
    if(j){
      gHa	+=ESgHa[j];
      gH	+=ESgHat[j];
    }
  }
  gHa	/=ESNp;
  gH	/=ESNp;
  for(j=0; j < ESNp1; j++){
    ESgHa[j]	-=gHa;
    ESgHat[j]	-=gH;
  }

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
      ESsz[ji]	=rs[0]+b*ESsn1[j];
      ESsza[ji]	=za;
      ESszt[ji]	=zt;
      ESszat[ji]=zat;
      ESsr[ji]	=r;
      ESsra[ji]	=ra;
      ESsrt[ji]	=rt;
      ESsrat[ji]=rat;

      D		=ra*zt-rt*za;
      Da	=raa*zt+ra*zat-rat*za-rt*zaa;
      Dt	=rat*zt+ra*ztt-rtt*za-rt*zat;
      Dat	=raat*zt+raa*ztt+ra*zatt-ratt*za-rtt*zaa-rt*zaat;

      K		=1./r;
      Ka	=-ra*K*K;
      Kt	=-rt*K*K;
      Kat	=(2.*rt*ra*K-rat)*K*K;
      
      ESgH[ji]	=D*K;
      ESgHa[ji]	=Da*K+D*Ka;
      ESgHt[ji]	=Dt*K+D*Kt;
      ESgHat[ji]=Dat*K+Dt*Ka+Da*Kt+D*Kat;
      gH	+=ESgH[ji];
      gHa	+=ESgHa[ji];

      ESaB[ji]	=pF*K*K;
      ESaBa[ji]	=(dpF*K+2.*pF*Ka)*K;
      ESaBt[ji]	=2.*pF*Kt*K;
      ESaBat[ji]=2.*dpF*K*Kt+2.*pF*(Kat*K+Kt*Ka);

#ifdef ForFRC
      ESaB[ji]	=0.;
      ESaBa[ji]	=0.;
      ESaBt[ji]	=0.;
      ESaBat[ji]=0.;
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
      ESaB[ji]	+=K*dgY*dgY;
      ESaBa[ji]	+=(Ka*dgY+2.*K*d2gY)*dgY;
      ESaBt[ji]	+=Kt*dgY*dgY;
      ESaBat[ji]+=(Kat*dgY+2.*Kt*d2gY)*dgY;

      D	=sqrt(ESaB[ji]);
      ESaB[ji]	=D;
      D		=EZcr2/D;
      ESaBa[ji]	*=D;
      ESaBt[ji]	*=D;
      ESaBat[ji]=(ESaBat[ji]-2.*ESaBa[ji]*ESaBt[ji])*D;
      ji++;
    }
    G	=1./ESNp;
    gH	*=G;
    gHa	*=G;
    j	=ji-ESNp;
    ESsra[ji]	=ESsra[j];
    ESsrt[ji]	=ESsrt[j];
    ESsrat[ji]	=ESsrat[j];
    ESsza[ji]	=ESsza[j];
    ESszt[ji]	=ESszt[j];
    ESszat[ji]	=ESszat[j];
    ESaB[ji]	=ESaB[j];
    ESaBa[ji]	=ESaBa[j];
    ESaBt[ji]	=ESaBt[j];
    ESaBat[ji]	=ESaBat[j];
    ESgH[ji]	=ESgH[j];
    ESgHa[ji]	=ESgHa[j];
    ESgHt[ji]	=ESgHt[j];
    ESgHat[ji]	=ESgHat[j];
    ji++;
    G	=1./gH;
    while(j < ji){
      ESgHa[j]	=(ESgHa[j]-gHa*ESgH[j]*G)*G;
      ESgH[j]	=(ESgH[j]-gH)*G;
      j++;
    }
  }
  return(0);
}

int ESBasicFunctionsPEST(int Na, int Np)
{
  int Np1,Na1;
  int i,j,ji,k,ki,kj,iK;
  double r,ra,rt,rat,raa,rtt,rttt,raat,ratt;
  double z,za,zt,zat,zaa,ztt,zttt,zaat,zatt;
  double rc[ESFp1],rs[ESFp1];
  double rct[ESFp1],rst[ESFp1],rca[ESFp1],rsa[ESFp1],rcaa[ESFp1],rsaa[ESFp1];
  double b,ba,baa,dgY,d2gY,pF,dpF;
  double D,Da,Dt,Dat,Dtt,K,Ka,Kt,Kat,Ktt,G,Ga,Gt,Gat,Gtt;
  double a,da,t,dt,s,s1,s2,cs,sn;
  double L,La;
  double Ght[ESNp1],Ghtc[ESFp1],Ghts[ESFp1];
  double Ghat[ESNp1],Ghatc[ESFp1],Ghats[ESFp1];

  double *sr,*sra,*srt,*srat;
  double *sz,*sza,*szt,*szat;
  double *aB,*aBa,*aBt,*aBat,B,Ba,Bt,Bat,Btt;
  double h,ha,ht,hat,htt;

  Np1	=Np+1;
  Na1	=Na+1;
  da	=(ESsa[ESNa]-ESsa[0])/Na;
  dt	=EZc2gp/Np;

  i	=Np1*Na1;
  sr	=(double*)malloc(4*i*sizeof(double));
  if(sr == NULL){
    printf("Failure in malloc for sr in ESBasicFunctionsPEST()\n");
    return(1);
  }
  sra	=sr	+i;
  srt	=sra	+i;
  srat	=srt	+i;
  
  sz	=(double*)malloc(4*i*sizeof(double));
  if(sz == NULL){
    printf("Failure in malloc for sz in ESBasicFunctionsPEST()\n");
    return(1);
  }
  sza	=sz	+i;
  szt	=sza	+i;
  szat	=szt	+i;
  
  aB	=(double*)malloc(4*i*sizeof(double));
  if(aB == NULL){
    printf("Failure in malloc for aB in ESBasicFunctionsPEST()\n");
    return(1);
  }
  aBa	=aB	+i;
  aBt	=aBa	+i;
  aBat	=aBt	+i;

  /* a=0 surface */
  ki	=ESNa1;
  rca[1]=2.*rcT1a[ki];
  rsa[1]=2.*rsT1a[ki];
  ba	=ESsb1a[0];
  pF	=ESaF[0];
  r	=ESaR0[0];
  K	=1./r;
  ji	=0;
  for(j=0; j < Np1; j++){
    t		=dt*j;
    cs		=cos(t);
    sn		=sin(t);

    ra		=rca[1]*cs+rsa[1]*sn;
    rt		=rsa[1]*cs-rca[1]*sn;
    rat		=rsa[1]*cs-rca[1]*sn;

    za		=ba*sn;
    zt		=ba*cs;
    zat		=ba*cs;

    sr[ji]	=r;
    sra[ji]	=ra;
    srt[ji]	=0.;
    srat[ji]	=rat;
    sz[ji]	=ESaZ0[0];
    sza[ji]	=za;
    szt[j]	=0.;
    szat[j]	=zat;
    
    Ka	=-ra*K*K;
    Kt	=0.;
    Kat	=-rat*K*K;
    
    aB[ji]	=pF*K;
    aBa[ji]	=pF*Ka;
    aBt[ji]	=0.;
    aBat[ji]	=pF*Kat;
    ji++;
  }

  for(i=1; i < Na1; i++){
    a	=da*i;
    ESSetSplA(a);
    splRA2(&b,&ba,&baa,ESsb,ESsb2a);
    splRA2(rc,rca,rcaa,rcT,rcT2a);
    rc[0]	+=ESaR0[0];
    splRA2(rs,rsa,rsaa,rsT,rsT2a);
    rs[0]	+=ESaZ0[0];
    ki	=0;
    iK	=1;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      splRA2(rc+k,rca+k,rcaa+k,rcT+ki,rcT2a+ki);
      rc[k]	*=2.;
      rst[k]	=-rc[k]*k;
      rca[k]	*=2.;
      rcaa[k]	*=2.;
      splRA2(rs+k,rsa+k,rsaa+k,rsT+ki,rsT2a+ki);
      rs[k]	*=2.;
      rct[k]	=rs[k]*k;
      rsa[k]	*=2.;
      rsaa[k]	*=2.;
      if(rct[k] != 0. || rca[k] != 0. || rcaa[k] != 0. ||
	 rst[k] != 0. || rsa[k] != 0. || rsaa[k] != 0.){
	iK	=k+1;
      }
    }
    /* PEST gq */
    L	=0.;
    La	=0.;
    for(j=0; j < ESNp; j++){
      cs	=EScs1[j];
      sn	=ESsn1[j];
      r		=rc[0];
      ra	=rca[0];
      raa	=rcaa[0];
      rt	=0.;
      rat	=0.;
      z		=rs[0]+b*sn;
      za	=rsa[0]+ba*sn;
      zat	=ba*cs;
      zaa	=rsaa[0]+baa*sn;
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
	rat	+=k*(rsa[k]*cs-rca[k]*sn);
	raa	+=rcaa[k]*cs+rsaa[k]*sn;
	rt	+=rct[k]*cs+rst[k]*sn;
      }
      Ght[j]	=(ra*zt-rt*za)/r;
      Ghat[j]	=(raa*zt+ra*zat-rat*za-rt*zaa-ra*Ght[j])/r;
      L		+=Ght[j];
      La	+=Ghat[j];
    }
    La	/=L;
    L	/=ESNp;
    s	=1./L;
    for(j=0; j < ESNp; j++){
      Ght[j]	*=s;
      Ghat[j]	=Ghat[j]*s-Ght[j]*La;
    }
    Ght[j]	=Ght[0];
    Ghat[j]	=Ghat[0];
    ESP2F(Ghtc,Ghts,Ght,ESFp1-1);
    ESP2F(Ghatc,Ghats,Ghat,ESFp1-1);
    Ghtc[0]	=0.;
    Ghatc[0]	=0.;
    for(k=1; k < ESFp1; k++){
      Ghtc[k]	*=2.;
      Ghts[k]	*=2.;
      Ghatc[k]	*=2.;
      Ghats[k]	*=2.;
    }
    splRA(&pF,&dpF,ESFF,ESFF2a);
    splRA(&dgY,&d2gY,ESdgY,ESdgY2a);
    d2gY	=dgY+a*d2gY;
    dgY		*=a;
    d2gY	*=2.*dgY;
    dgY		*=dgY;
    t		=0.;
    for(j=0; j < Np; j++){
      h		=dt*j;
      kj	=0;
      do{
	s	=t-h;
	s1	=1.;
	ha	=0.;
	for(k=1; k < ESFp1; k++){
	  cs	=cos(k*t);
	  sn	=sin(k*t);
	  s1	+=Ghtc[k]*cs+Ghts[k]*sn;
	  s	+=(Ghtc[k]*sn+Ghts[k]*(1.-cs))/k;
	  ha	+=(Ghatc[k]*sn+Ghats[k]*(1.-cs))/k;
	}
	t	-=s/s1;
	kj++;
      }while(kj < 10 && fabs(s) > 1e-8);
      s		=t;
      cs	=cos(s);
      sn	=sin(s);
      r		=rc[0];
      ra	=rca[0];
      rt	=0.;
      rat	=0.;
      raa	=rcaa[0];
      rtt	=0.;
      rttt	=0.;
      raat	=0.;
      ratt	=0.;

      sz[ji]	=rs[0]+b*sn;
      za	=rsa[0]+ba*sn;
      zt	=b*cs;
      zat	=ba*cs;
      zaa	=rsaa[0]+baa*sn;
      ztt	=-b*sn;
      zttt	=-b*cs;
      zaat	=baa*cs;
      zatt	=-ba*sn;
      for(k=1; k < iK; k++){
	s	=k*t;
	cs	=cos(s);
	sn	=sin(s);
	D	=k*k;
	s	=rc[k]*cs+rs[k]*sn;
	r	+=s;
	rtt	-=D*s;
	s	=rca[k]*cs+rsa[k]*sn;
	ra	+=s;
	ratt	-=D*s;
	s	=rct[k]*cs+rst[k]*sn;
	rt	+=s;
	rttt	-=D*s;
	rat	+=k*(rsa[k]*cs-rca[k]*sn);
	raa	+=rcaa[k]*cs+rsaa[k]*sn;
	raat	+=k*(rsaa[k]*cs-rcaa[k]*sn);
      }
      D		=ra*zt-rt*za;
      Da	=raa*zt+ra*zat-rat*za-rt*zaa;
      Dt	=rat*zt+ra*ztt-rtt*za-rt*zat;
      Dat	=raat*zt+raa*ztt+ra*zatt-ratt*za-rtt*zaa-rt*zaat;
      Dtt	=ratt*zt+ra*zttt-rttt*za+2.*(rat*ztt-rtt*zat)-rt*zatt;

      s		=1./(r*L);
      ht	=D*s;
      hat	=(Da-D*(ra/r+La))*s;
      htt	=(Dt-D*rt/r)*s;

      ht	=1./ht;
      ha	*=-ht;
      hat	=-(hat+ha*htt)*ht*ht;

      sza[ji]	=za+zt*ha;
      szt[ji]	=zt*ht;
      szat[ji]	=(zat+ztt*ha)*ht+zt*hat;

      sr[ji]	=r;
      sra[ji]	=ra+rt*ha;;
      srt[ji]	=rt*ht;
      srat[ji]	=(rat+rtt*ha)*ht+rt*hat;

      s		=1./r;
      K		=s*s;
      Ka	=-2.*ra*K*s;
      Kt	=-2.*rt*K*s;
      Kat	=(6.*ra*rt*K-2.*rat*s)*K;
      Ktt	=(6.*rt*rt*K-2.*rtt*s)*K;

      B		=pF*K;
      Ba	=dpF*K+pF*Ka;
      Bt	=pF*Kt;
      Btt	=pF*Ktt;
      Bat	=dpF*Kt+pF*Kat;

      s		=1./D;
      G		=s*s;
      Ga	=-2.*Da*G*s;
      Gt	=-2.*Dt*G*s;
      Gat	=(6.*Da*Da*G-2.*Dat*s)*G;
      Gtt	=(6.*Dt*Dt*G-2.*Dtt*s)*G;

      D		=K*G;
      Da	=Ka*G+K*Ga;
      Dt	=Kt*G+K*Gt;
      Dat	=Kat*G+Kt*Ga+Ka*Gt+K*Gat;
      Dtt	=Ktt*G+2.*Kt*Gt+K*Gtt;

      G		=rt*rt+zt*zt;
      Ga	=2.*(rat*rt+zat*zt);
      Gt	=2.*(rtt*rt+ztt*zt);
      Gat	=2.*(ratt*rt+rtt*rat+zatt*zt+ztt*zat);
      Gtt	=2.*(rttt*rt+rtt*rtt+zttt*zt+ztt*ztt);

      K		=D*G;
      Ka	=Da*G+D*Ga;
      Kt	=Dt*G+D*Gt;
      Kat	=Dat*G+Dt*Ga+Da*Gt+D*Gat;
      Ktt	=Dtt*G+Dt*Gt+Dt*Gt+D*Gtt;

      B		+=K*dgY;
      Ba	+=Ka*dgY+K*d2gY;
      Bt	+=Kt*dgY;
      Bat	+=Kat*dgY+Kt*d2gY;
      Btt	+=Ktt*dgY;

      aB[ji]	=B;
      aBa[ji]	=Ba+Bt*ha;;
      aBt[ji]	=Bt*ht;
      aBat[ji]	=(Bat+Btt*ha)*ht+Bt*hat;

      D		=sqrt(B);
      aB[ji]	=D;
      D		=EZcr2/D;
      aBa[ji]	*=D;
      aBt[ji]	*=D;
      aBat[ji]	=(aBat[ji]-2.*aBa[ji]*aBt[ji])*D;
      ji++;
    }
    j	=ji-Np;
    sr[ji]	=sr[j];
    sra[ji]	=sra[j];
    srt[ji]	=srt[j];
    srat[ji]	=srat[j];

    sz[ji]	=sz[j];
    sza[ji]	=sza[j];
    szt[ji]	=szt[j];
    szat[ji]	=szat[j];

    aB[ji]	=aB[j];
    aBa[ji]	=aBa[j];
    aBt[ji]	=aBt[j];
    aBat[ji]	=aBat[j];
    ji++;
  }
  ESWriteESI("PEST",sr,sz,aB,NULL,Na,Np);
  free(aB);
  free(sz);
  free(sr);
  return(0);
}

int ESFirstBiSpline()
{
  Spi0	=0;
  Spi1	=0;
  Spj0	=0;
  Spj1	=0;
  ha	=ESsa[1]-ESsa[0];
  rha	=1./ha;
  hq	=ESgt[1];
  rhq	=1./hq;
  return(0);
}
double YYa,YYx,YYX,YYq,YYy,YYY;
double XYa,XYx,XYX,XYq,XYy,XYY;
static double XSC1,XSC2,xSC1,xSC2;

int ESBiSplRAP(double *r,double *ra,double *rq
	       ,double *z,double *za,double *zq
	       ,double *B,double *Ba,double *Bq
	       ,double *gh,double *gha,double *ghq
	       ,double a,double gq)
{
  int j00,j01,j10,j11;
  double X,XX,x,xx,Y,YY,y,yy,YSC1,YSC2,ySC1,ySC2;
  double f0,f1,fq0,fq1;
  
  if(a < 0. && a > ESsa[ESNa]){
    return(1);
  }
  Spi0	=a*rha;
  if(Spi0 >= ESNa){
    Spi0	=ESNa-1;
  }
  Spi1	=Spi0+1;
  x	=(a-ESsa[Spi0])*rha;
  X	=1.-x;

#ifdef H
  Spj1	=gq*EZcr2gp;
  if(gq > 0.){
    Spj1++;
  }
  gq	=EZc2gp*Spj1-gq;
  Spj0	=gq*rhq;
  if(Spj0 >= ESNp){
    Spj0	-=ESNp;
    gq		-=EZc2gp;
  }
  if(Spj0 < 0){
    Spj0+=ESNp;
    gq	+=EZc2gp;
  }
  Spj1	=Spj0+1;
#endif
  Spj1	=gq*EZcr2gp;
  if(gq < 0.){
    Spj1--;
  }
  gq	-=EZc2gp*Spj1;
  Spj0	=gq*rhq;
  if(Spj0 >= ESNp){
    Spj0	-=ESNp;
    gq	-=EZc2gp;
  }
  if(Spj0 < 0){
    Spj0	+=ESNp;
    gq	+=EZc2gp;
  }
  Spj1	=Spj0+1;
  Y	=(gq-ESgt[Spj0])*rhq;
  y	=1.-Y;
  Spj0	=ESNp-Spj1;
  Spj1	=Spj0+1;

  XYa	=a;
  XYx	=x;
  XYX	=X;

  XYq	=gq;
  XYy	=y;
  XYY	=Y;

  XX	=X*X;
  xx	=x*x;
  XSC1	=XX*(3.-2.*X);
  XSC2	=XX*x*ha;
  xSC1	=xx*(3.-2.*x);
  xSC2	=-xx*X*ha;

  YY	=Y*Y;
  yy	=y*y;
  YSC1	=YY*(3.-2.*Y);
  YSC2	=YY*y*hq;
  ySC1	=yy*(3.-2.*y);
  ySC2	=-yy*Y*hq;

  xx	=3.*x*X;
  X	-=xx;
  x	-=xx;
  xx	*=2.*rha;

  yy	=3.*y*Y;
  Y	-=yy;
  y	-=yy;
  yy	*=2.*rhq;

  j00	=ESNp1*Spi0+Spj0;
  j01	=j00+1;
  j10	=j00+ESNp1;
  j11	=j10+1;

  f0	=XSC1*ESsr [j00]+xSC1*ESsr [j10]+XSC2*ESsra [j00]+xSC2*ESsra [j10];
  f1	=XSC1*ESsr [j01]+xSC1*ESsr [j11]+XSC2*ESsra [j01]+xSC2*ESsra [j11];
  fq0	=XSC1*ESsrt[j00]+xSC1*ESsrt[j10]+XSC2*ESsrat[j00]+xSC2*ESsrat[j10];
  fq1	=XSC1*ESsrt[j01]+xSC1*ESsrt[j11]+XSC2*ESsrat[j01]+xSC2*ESsrat[j11];
  *r	=YSC1*f0+ySC1*f1+YSC2*fq0+ySC2*fq1;
  *rq	=-(yy*(f1-f0)+Y*fq0+y*fq1);

  f0	=XSC1*ESsz [j00]+xSC1*ESsz [j10]+XSC2*ESsza [j00]+xSC2*ESsza [j10];
  f1	=XSC1*ESsz [j01]+xSC1*ESsz [j11]+XSC2*ESsza [j01]+xSC2*ESsza [j11];
  fq0	=XSC1*ESszt[j00]+xSC1*ESszt[j10]+XSC2*ESszat[j00]+xSC2*ESszat[j10];
  fq1	=XSC1*ESszt[j01]+xSC1*ESszt[j11]+XSC2*ESszat[j01]+xSC2*ESszat[j11];
  *z	=YSC1*f0+ySC1*f1+YSC2*fq0+ySC2*fq1;
  *zq	=-(yy*(f1-f0)+Y*fq0+y*fq1);
  
  f0	=XSC1*ESaB [j00]+xSC1*ESaB [j10]+XSC2*ESaBa [j00]+xSC2*ESaBa [j10];
  f1	=XSC1*ESaB [j01]+xSC1*ESaB [j11]+XSC2*ESaBa [j01]+xSC2*ESaBa [j11];
  fq0	=XSC1*ESaBt[j00]+xSC1*ESaBt[j10]+XSC2*ESaBat[j00]+xSC2*ESaBat[j10];
  fq1	=XSC1*ESaBt[j01]+xSC1*ESaBt[j11]+XSC2*ESaBat[j01]+xSC2*ESaBat[j11];
  *B	=YSC1*f0+ySC1*f1+YSC2*fq0+ySC2*fq1;
  *Bq	=-(yy*(f1-f0)+Y*fq0+y*fq1);

  f0	=XSC1*ESgH [j00]+xSC1*ESgH [j10]+XSC2*ESgHa [j00]+xSC2*ESgHa [j10];
  f1	=XSC1*ESgH [j01]+xSC1*ESgH [j11]+XSC2*ESgHa [j01]+xSC2*ESgHa [j11];
  fq0	=XSC1*ESgHt[j00]+xSC1*ESgHt[j10]+XSC2*ESgHat[j00]+xSC2*ESgHat[j10];
  fq1	=XSC1*ESgHt[j01]+xSC1*ESgHt[j11]+XSC2*ESgHat[j01]+xSC2*ESgHat[j11];
  *gh	=YSC1*f0+ySC1*f1+YSC2*fq0+ySC2*fq1;
  *ghq	=-(yy*(f1-f0)+Y*fq0+y*fq1);

  f0	=xx*(ESsr [j10]-ESsr [j00])+X*ESsra [j00]+x*ESsra [j10];
  f1	=xx*(ESsr [j11]-ESsr [j01])+X*ESsra [j01]+x*ESsra [j11];
  fq0	=xx*(ESsrt[j10]-ESsrt[j00])+X*ESsrat[j00]+x*ESsrat[j10];
  fq1	=xx*(ESsrt[j11]-ESsrt[j01])+X*ESsrat[j01]+x*ESsrat[j11];
  *ra	=YSC1*f0+ySC1*f1+YSC2*fq0+ySC2*fq1;

  f0	=xx*(ESsz [j10]-ESsz [j00])+X*ESsza [j00]+x*ESsza [j10];
  f1	=xx*(ESsz [j11]-ESsz [j01])+X*ESsza [j01]+x*ESsza [j11];
  fq0	=xx*(ESszt[j10]-ESszt[j00])+X*ESszat[j00]+x*ESszat[j10];
  fq1	=xx*(ESszt[j11]-ESszt[j01])+X*ESszat[j01]+x*ESszat[j11];
  *za	=YSC1*f0+ySC1*f1+YSC2*fq0+ySC2*fq1;

  f0	=xx*(ESaB [j10]-ESaB [j00])+X*ESaBa [j00]+x*ESaBa [j10];
  f1	=xx*(ESaB [j11]-ESaB [j01])+X*ESaBa [j01]+x*ESaBa [j11];
  fq0	=xx*(ESaBt[j10]-ESaBt[j00])+X*ESaBat[j00]+x*ESaBat[j10];
  fq1	=xx*(ESaBt[j11]-ESaBt[j01])+X*ESaBat[j01]+x*ESaBat[j11];
  *Ba	=YSC1*f0+ySC1*f1+YSC2*fq0+ySC2*fq1;

  f0	=xx*(ESgH [j10]-ESgH [j00])+X*ESgHa [j00]+x*ESgHa [j10];
  f1	=xx*(ESgH [j11]-ESgH [j01])+X*ESgHa [j01]+x*ESgHa [j11];
  fq0	=xx*(ESgHt[j10]-ESgHt[j00])+X*ESgHat[j00]+x*ESgHat[j10];
  fq1	=xx*(ESgHt[j11]-ESgHt[j01])+X*ESgHat[j01]+x*ESgHat[j11];
  *gha	=YSC1*f0+ySC1*f1+YSC2*fq0+ySC2*fq1;
  return(0);
}

int ESBiSplRA(double *F,double *dgY,double *T,double *P)
{
  *F	=XSC1*ESaF[Spi0]+xSC1*ESaF[Spi1]+XSC2*ESaF1a[Spi0]+xSC2*ESaF1a[Spi1];
  *T	=XSC1*ESaT[Spi0]+xSC1*ESaT[Spi1]+XSC2*ESaT1a[Spi0]+xSC2*ESaT1a[Spi1];
  *P	=XSC1*ESPs[Spi0]+xSC1*ESPs[Spi1]+XSC2*ESPs1a[Spi0]+xSC2*ESPs1a[Spi1];
  *dgY=XSC1*ESdgY[Spi0]+xSC1*ESdgY[Spi1]+XSC2*ESdgY1a[Spi0]+xSC2*ESdgY1a[Spi1];
  
  return(0);
}

static int nst=0,ist=0;
static double Yr[1024],Yra[1024],Yrq[1024],Yz[1024],Yza[1024],Yzq[1024]
,YB[1024],YBa[1024],YBq[1024],Ya[1024],Yq[1024];

#ifdef XWIN
void GCLagrangeEq(int *n,double *x,double *Y,double *dY)
{
  double r,ra,rq,z,za,zq,B,Ba,Bq,dgY,F,D;
  double a11,a12,a21,a22;
  double a,q;
  int i;

  a	=Y[1];
  q	=Y[2];
  if(ESBiSplRAP(&r,&ra,&rq,&z,&za,&zq,&B,&Ba,&Bq,&D,&D,&D,Y[1],Y[2])){
    istate	=3;
    return;
  }

  ESSetSplA(Y[1]);
  splRA(&dgY,NULL,ESdgY,ESdgY2a);
  dgY	*=Y[1];
  splRA(&a21,NULL,ESaF,ESaF2a);
  splRA(dY+2,NULL,ESPs,ESPs2a);
  splRA(dY+3,NULL,ESaT,ESaT2a);

  ESBiSplRA(&a21,&dgY,dY+3,dY+2);
  dgY	*=Y[1];
  dY[3]	*=Y[0];

  a22	=za*rq-zq*ra;
  a12	=a22*((a21+dY[3])/r+Y[0]*r*dY[2]);
  D	=-dgY/(r*a22);
  a11	=(rq*rq+zq*zq)*D;
  D	*=ra*rq+za*zq;
  a22	=dgY*(1.+dY[3]/a21);
  r	=1./(a11*a22-a12*a21);
  a11	*=r;
  a12	*=r;
  a21	*=r;
  a22	*=r;
  D	*=r;
#ifdef H
  r	=cgWc*Y[0]*B;
#endif
  r	=Y[0]*B;

  B	*=r;

#ifdef H
  r	=r*Y[0]+cgWc*gmP;
#endif
  r	=r*Y[0]+gmP;
  Ba	*=r;
  Bq	*=-r;

  dY[0]	=a22*Bq;
  dY[1]	=-a21*Bq;
  dY[2]	=a22*B-a21*Ba;
  dY[3]	=-a12*B+a11*Ba+D*Bq;
  
#ifdef H
  i	=1;
  gcmotion_(Yr,Yra,Yrq,Yz,Y,Y+1,Y+2,&gmP,&i);
  if(Yr[0] != dY[0] || Yra[0] != dY[1] || Yrq[0] !=dY[2] || Yz[0] !=dY[3]){
    EZout("sdddddd","dY=",Yr[0],Yr[0]-dY[0],Yra[0],Yra[0]-dY[1]
	,Yrq[0],Yrq[0]-dY[2]);
    EZout("sdd","dY=",Yz[0],Yz[0]-dY[3]);
  }
#endif  
  return;
}

void GCLagrangeEqA(int *n,double *x,double *Y,double *dY)
{
  double r,ra,rq,z,za,zq,B,Ba,Bq,dgY,F,D;
  double a11,a12,a21,a22;
  double a,q;
  int i;

#ifdef H
  a	=Y[1];
  q	=Y[2];
  if(ESBiSplRAP(&r,&ra,&rq,&z,&za,&zq,&B,&Ba,&Bq,&D,&D,&D,Y[1],Y[2])){
    istate	=3;
    return;
  }

  ESSetSplA(Y[1]);
  splRA(&dgY,NULL,ESdgY,ESdgY2a);
  dgY	*=Y[1];
  splRA(&a21,NULL,ESaF,ESaF2a);
  splRA(dY+2,NULL,ESPs,ESPs2a);
  splRA(dY+3,NULL,ESaT,ESaT2a);

  ESBiSplRA(&a21,&dgY,dY+3,dY+2);
  dgY	*=Y[1];
  dY[3]	*=Y[0];

  a22	=za*rq-zq*ra;
  a12	=a22*((a21+dY[3])/r+Y[0]*r*dY[2]);
  D	=-dgY/(r*a22);
  a11	=(rq*rq+zq*zq)*D;
  D	*=ra*rq+za*zq;
  a22	=dgY*(1.+dY[3]/a21);
  r	=1./(a11*a22-a12*a21);
  a11	*=r;
  a12	*=r;
  a21	*=r;
  a22	*=r;
  D	*=r;
#ifdef H
  r	=cgWc*Y[0]*B;
#endif
  r	=Y[0]*B;

  B	*=r;

#ifdef H
  r	=r*Y[0]+cgWc*gmP;
#endif
  r	=r*Y[0]+gmP;
  Ba	*=r;
  Bq	*=-r;

  dY[0]	=a22*Bq;
  dY[1]	=-a21*Bq;
  dY[2]	=a22*B-a21*Ba;
  dY[3]	=-a12*B+a11*Ba+D*Bq;
#endif  
  i	=1;
  gcmotion_(Yr,Yra,Yrq,Yz,Y,Y+1,Y+2,&gmP,&i);
  dY[0]	=Yr[0];
  dY[1]	=Yra[0];
  dY[2]	=Yrq[0];
  dY[3]	=Yz[0];
#ifdef H
  if(Yr[0] != dY[0] || Yra[0] != dY[1] || Yrq[0] !=dY[2] || Yz[0] !=dY[3]){
    EZout("sdddddd","dY=",Yr[0],Yr[0]-dY[0],Yra[0],Yra[0]-dY[1]
	,Yrq[0],Yrq[0]-dY[2]);
    EZout("sdd","dY=",Yz[0],Yz[0]-dY[3]);
  }
#endif  
  return;
}

int ES4Output(double a0,double a1,double gt0, double gt1,int na, int nt)
{

  return(0);
}

int ESInitOrbLSODE(int n)
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

int ESDeInitOrbLSODE()
{
  free(rwork);
  rwork	=NULL;
  free(iwork);
  iwork	=NULL;
  return(0);
}

void GCOrbitPlot2D()
{
  double t,tout,A,R[256],Z[256],gq0,gf0;
  int n;

  ZColor(14);
  t	=0.;
  tout	=rwork[0];
  gf0	=yGC[3];
  gq0	=yGC[2];
  ES2DMapFlx2Lab(R,Z,&A,&A,&A,&A,yGC[1],-yGC[2]);
  n	=1;
  ist	=0;
  while(t < rwork[0]){
    lsode_(GCLagrangeEq,&neq,yGC,&t,&tout,&itol,&reltol,vatol,&itask,
	   &istate,&iopt,rwork,&lrw,iwork,&liw,NULL,&mf);
    ES2DMapFlx2Lab(R+n,Z+n,&A,&A,&A,&A,yGC[1],-yGC[2]);
    n++;
    if(n == 256){
      ZPlotPolyLine(R,Z,n);
      n--;
      R[0]	=R[n];
      Z[0]	=Z[n];
      n	=1;
    }
  }
  nst	=ist;
  if(n > 1){
    ZPlotPolyLine(R,Z,n);
  }
  EZout("sIIIdd","N steps=",iwork[10],iwork[11],istate
      ,(yGC[2]-gq0)*EZcr2gp,(yGC[3]-gf0)*EZcr2gp);
  return;
}

void GCOrbitPlot2DB()
{
  double t,tout,A,R[256],Z[256],gq0,gf0;
  int n;

  ZColor(14);
  t	=0.;
  tout	=rwork[0];
  gf0	=yGC[3];
  gq0	=yGC[2];
  ES2DMapFlx2Lab(R,Z,&A,&A,&A,&A,yGC[1],-yGC[2]);
  n	=1;
  ist	=0;
  while(t < rwork[0]){
    lsode_(GCLagrangeEqA,&neq,yGC,&t,&tout,&itol,&reltol,vatol,&itask,
	   &istate,&iopt,rwork,&lrw,iwork,&liw,NULL,&mf);
    ES2DMapFlx2Lab(R+n,Z+n,&A,&A,&A,&A,yGC[1],-yGC[2]);
    n++;
    if(n == 256){
      ZPlotPolyLine(R,Z,n);
      n--;
      R[0]	=R[n];
      Z[0]	=Z[n];
      n	=1;
    }
  }
  nst	=ist;
  if(n > 1){
    ZPlotPolyLine(R,Z,n);
  }
  EZout("sIIIdd","N steps=",iwork[10],iwork[11],istate
      ,(yGC[2]-gq0)*EZcr2gp,(yGC[3]-gf0)*EZcr2gp);
  return;
}

void GCOrbitPlot2DA()
{
  double t,h,A;
  int k,n,m;

  ZColor(14);
  t	=0.;
  n	=1;

  m	=1;
  for(k=0; k < nPart; k++){
    ES2DMapFlx2Lab(rPrt0+k,zPrt0+k,&A,&A,&A,&A,Xa[k],-Xgq[k]);
  }

  gcmotion_(dXgr,dXa,dXgq,dXgf,Xgr,Xa,Xgq,Xgm,&nPart);
  h	=20.*EZc2gp;
  ist	=0;
  while(t < rwork[0]){
    GCRKStep(&t,h);
    for(k=0; k < nPart; k++){
      if(m){
	ES2DMapFlx2Lab(rPrt1+k,zPrt1+k,&A,&A,&A,&A,Xa[k],-Xgq[k]);
      }
      else{
	ES2DMapFlx2Lab(rPrt0+k,zPrt0+k,&A,&A,&A,&A,Xa[k],-Xgq[k]);
      }
      ZPlotLine(rPrt1[k],zPrt1[k],rPrt0[k],zPrt0[k]);
    }
    m	^=1;
    ist++;
  }
  nst	=ist;
  return;
}
#endif
int ESIWriteDensity(int Na,int Fic)
{
  int i,Na1;
  double a,da,ne,Te,Ti,p,dp,f;
  FILE *lf;
  char FNm[16];
  extern char ESMessage[];

  sprintf(FNm,"esiT.%2.2d",Fic);
  lf	=fopen(FNm,"w");
  if(lf == NULL){
    sprintf(ESMessage,"%s cannot be open",FNm);
    return(1);
  }
  Na1	=Na+1;
  da	=(ESsa[ESNa]-ESsa[0])/Na;
  Te	=15.;
  Ti	=15.;
  f	=1./(1.60022e-2*EZcgm0*(Te+Ti));

  fprintf(lf,"!!! Do not edit this file. %s\n",ESMessage);
  fprintf(lf,"%24s%24s%24s\n","Te=Ti=15keV     a"
	  ,"ne [10^20/m^3]","dne/da [10^20/m^3]");
  for(i=0; i < Na1; i++){
    a	=da*i;
    ESSetSplA(a);
    splRA(&p,&dp,ESsp,ESsp2a);
    ne	=f*p;
    fprintf(lf,"%24.16e%24.16e%24.16e\n",a,f*p,f*dp);
  }
  fclose(lf);
  return(0);
}

#ifdef XWIN
int ESOrbits()
{
  double X[2],Y[2],R1,R2;
  double r,ra,rq,z,za,zq,B,Ba,Bq,gh,gha,ghq,gq;
  double V,Tc,T0;
  static double A,gt=0.,a=0.5,mi=2.,ei=1.;
  static double dt=1.,t=3000.;
  static double E=1.;
  static double p=0.4,gr,grn,grb,Vn,Vb;
  static int Fl=0,Flmi=1;
  static int n=1000;
  int i,j,k,L,ist;
  
  j	=ESNp1*ESNa;
  i	=ESNp/2;
  X[0]	=ESsr[j+i];
  X[1]	=ESsr[j];
  i	=ESNp/4;
  Y[0]	=ESsz[j+3*i];
  Y[1]	=ESsz[j+i];

  ESBasicFunctions();
  ESFirstBiSpline();

  esiread_("esiB.00");

  ESInitOrbLSODE(neq);
  esilink2c_(XF,XFa,XgFa,XgFaa,XgYa,XgYaa,XT,XTa,XP,XPa,Xr,Xra,Xrq
		     ,Xz,Xza,Xzq,XB,XBa,XBq,Xgh,Xgha,Xghq,NULL);
  L	=0;
  iwork[11-1]	=0; /* the number of steps taken for the problem so far */
  iwork[12-1]	=0; /* the number of f evaluations for the problem so far. */
#ifndef Tbl_Orbits
  switch(Flmi){
  case 0:
    mi	=1.;
    ei	=1.;
    break;
  case 1:
    mi	=2.;
    ei	=1.;
    break;
  case 2:
    mi	=3.;
    ei	=1.;
    break;
  case 3:
    mi	=4.;
    ei	=2.;
    break;
  }
  cgWc		=ce*ei*1e+4/(mi*cmi*cc);
  cV1keV	=sqrt(2.*cEkeV/(mi*cmi))*1e-2;
  cgr1keV	=cV1keV/cgWc;
  T0		=EZc2gp/cgWc;

  if(a < ESsa[0]){
    a	=ESsa[0];
  }
  if(a > ESsa[ESNa]){
    a	=ESsa[ESNa];
  }
  while(gt >= 1.){
    gt	-=1.;
  }
  while(gt < 0.){
    gt	+=1.;
  }
  switch(Fl){
  case 0:
    ES2DMapFlx2Lab(&r,&z,&gq,&gq,&gq,&gq,a,EZc2gp*gt);
    break;
  case 1:
    if(ES1DMapZ2gt(&gq,&R2,&R1,z) == 0 && R1 <= r && r <= R2){
      gq=EZc2gp*gt;
      A	=a;
      if(ES2DMapLab2Flx(&A,&gq,r,z) == 0){
	a	=A;
	gt	=EZcr2gp*gq;
      }
    }
    break;
  }
  gq	=EZc2gp*(1.-gt);
  ESBiSplRAP(&Vb,&ra,&rq,&Vn,&za,&zq,&B,&Ba,&Bq,&gh,&gha,&ghq,a,gq);

  V	=cV1keV*sqrt(E);
  Vb	=p*V;
  Vn	=sqrt(1.-p*p)*V;
  gr	=V/(cgWc*B);
  grb	=Vb/(cgWc*B);
  grn	=Vn/(cgWc*B);

  Tc	=1./B;
  n	=0;
  ist	=0;
  if(L == 0){
    CbUserPlot	=NULL;
  }
  else{
    gmP		=EZcr2*Vn*Vn/(B*cgWc*cgWc);
    yGC[0]	=Vb/(B*cgWc);
    yGC[1]	=a;
    yGC[2]	=gq;
    yGC[3]	=0.;

    switch(L){
    case 1:
      CbUserPlot	=GCOrbitPlot2D;
      break;
    case 2:
      if(GCInitRKSolv(nPart) == 0){
	for(k=0; k < nPart; k++){
	  A	=0.05+0.9*k/nPart;
	  A	=1.-k/nPart;
	  Xa[k] =a;
	  Xgq[k]=gq;
	  Xgf[k]=0.;
	  Xgr[k]=A*gr;
	  Xgm[k]=EZcr2*(1.-A*A)*gr*gr*B;
	}
	CbUserPlot	=GCOrbitPlot2DA;
      }
      break;
    case 3:
      CbUserPlot	=GCOrbitPlot2DB;
      break;
    }
    mf		=10;
    itol	=2;	/* vector tolerance, 1 for scalar toleranc */
    itask	=5;
    iopt	=1;	/* optional input ON */
    rwork[1-1]	=t*EZc2gp;
    rwork[2-1]	=0.;
    rwork[3-1]	=0.;
    rwork[4-1]	=0.;
    rwork[5-1]	=10;/* the step size to be attempted on the first step. */
    rwork[6-1]	=0.1*rwork[0];/* the maximum absolute step size allowed. */
    rwork[7-1]	=0.;	/* the minimum absolute step size allowed. */

    iwork[4]	=12;	/* Maximum order */
    iwork[5]	=500;
    iwork[6]	=10;

    reltol	=0.;
    vatol[0]	=sqrt(grb*grb+grn*grn)*abstol;
    vatol[1]	=abstol;
    vatol[2]	=EZc2gp*abstol;
    vatol[3]	=EZc2gp*abstol;

    istate	=1;

    L	=0;
  }
#endif /*Tbl_Orbits*/
  GCDeInitRKSolv();
  ESDeInitOrbLSODE();
  ESFreeRZspaceEq();
  return(0);
}
#endif
int ESFractionOfTrapped(double *frac,int i)
{
  int j,k,ki,kj,iK;
  double pB[ESNp1],D[ESNp1],D0,aph,Bmx;
  double rc[ESFp1],rs[ESFp1];
  double rct[ESFp1],rst[ESFp1],rca[ESFp1],rsa[ESFp1];
  double r,ra,rt,za,zt;
  double b,ba,dgY,pF;
  double cs,sn;
  double gt[ESNp1];

  if(i == 0){
    return(0);
  }

  dgY		=ESsa[i]*ESdgY[i];
  pF		=ESFF[i];
  rc[0]		=ESaR0[0]+rcT[i];
  rs[0]		=ESaZ0[0]+rsT[i];
  rca[0]	=rcT1a[i];
  rsa[0]	=rsT1a[i];
  b		=ESsb[i];
  ba		=ESsb1a[i];
  ki	=i;
  iK	=1;
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    rc[k]	=2.*rcT[ki];
    rct[k]	=2.*k*rsT[ki];
    rca[k]	=2.*rcT1a[ki];
    rs[k]	=2.*rsT[ki];
    rst[k]	=-2.*k*rcT[ki];
    rsa[k]	=2.*rsT1a[ki];
    if(rc[k] != 0. || rca[k] != 0. || rs[k] != 0. || rsa[k] != 0.){
      iK	=k+1;
    }
  }

  aph	=0.;
  Bmx	=0.;
  D0	=0.;
  for(j=0; j < ESNp; j++){
    r	=rc[0];
    ra	=rca[0];
    rt	=0.;
    cs	=EScs1[j];
    sn	=ESsn1[j];
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
      r		+=rc[k]*cs+rs[k]*sn;
      ra	+=rca[k]*cs+rsa[k]*sn;
      rt	+=rct[k]*cs+rst[k]*sn;
    }
    D[j]	=r*(ra*zt-rt*za);
    pB[j]	=pF/(r*r)+(rt*rt+zt*zt)*dgY*dgY/(D[j]*D[j]);
    aph		+=pB[j]*D[j];
    D0		+=D[j];
    if(Bmx < pB[j]){
      Bmx	=pB[j];
    } 
  }
  D[ESNp]	=D[0];
  pB[ESNp]	=pB[0];
  aph	/=Bmx*D0;
  for(j=0; j < ESNp1; j++){
    pB[j]	=sqrt(pB[j]/Bmx);
    gt[j]	=ESgt[j]*EZcr2gp;
  }
  iK	=41;
  ba	=1./(iK-1);
  ra	=0.;
  pF	=0.;
  for(k=1; k < iK; k++){
    b	=ba*k;
    za	=ra;
    ra	=0.;
    for(j=0; j < ESNp; j++){
      ra	+=sqrt(1.-b*pB[j])*D[j];
    }
    ra	=b*D0/ra;
    pF	+=EZcr2*ba*(ra+za);
  }
  *frac	=1.-0.75*aph*pF;
  return(0);
}

int ESFractionOfTrapped1(double *frac,double a)
{
  int j,k,ki,kj,iK;
  double pB[ESNp1],D[ESNp1],D0,aph,Bmx;
  double rc[ESFp1],rs[ESFp1];
  double rct[ESFp1],rst[ESFp1],rca[ESFp1],rsa[ESFp1];
  double r,ra,rt,za,zt;
  double b,ba,dgY,pF;
  double cs,sn;
  double gt[ESNp1];

  ESSetSplA(a);
  splRA(&dgY,NULL,ESdgY,ESdgY2a);
  dgY		*=a;
  splRA(&pF,NULL,ESFF,ESFF2a);

  splRA(rc,rca,rcT,rcT2a);
  rc[0]		+=ESaR0[0];
  splRA(rs,rsa,rsT,rsT2a);
  rs[0]		+=ESaZ0[0];
  splRA(&b,&ba,ESsb,ESsb2a);
  ki	=0;
  iK	=1;
  for(k=1; k < ESFp1; k++){
    ki	+=ESNa1;
    splRA(rc+k,rca+k,rcT+ki,rcT2a+ki);
    rc[k]	*=2.;
    rca[k]	*=2.;
    rst[k]	=-k*rc[k];
    splRA(rs+k,rsa+k,rsT+ki,rsT2a+ki);
    rs[k]	*=2.;
    rsa[k]	*=2.;
    rct[k]	=k*rs[k];
    if(rc[k] != 0. || rca[k] != 0. || rs[k] != 0. || rsa[k] != 0.){
      iK	=k+1;
    }
  }
  aph	=0.;
  Bmx	=0.;
  D0	=0.;
  for(j=0; j < ESNp; j++){
    r	=rc[0];
    ra	=rca[0];
    rt	=0.;
    cs	=EScs1[j];
    sn	=ESsn1[j];
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
      r		+=rc[k]*cs+rs[k]*sn;
      ra	+=rca[k]*cs+rsa[k]*sn;
      rt	+=rct[k]*cs+rst[k]*sn;
    }
    D[j]	=r*(ra*zt-rt*za);
    pB[j]	=pF/(r*r)+(rt*rt+zt*zt)*dgY*dgY/(D[j]*D[j]);
    aph		+=pB[j]*D[j];
    D0		+=D[j];
    if(Bmx < pB[j]){
      Bmx	=pB[j];
    } 
  }
  D[ESNp]	=D[0];
  pB[ESNp]	=pB[0];
  aph	/=Bmx*D0;
  for(j=0; j < ESNp1; j++){
    pB[j]	=sqrt(pB[j]/Bmx);
    gt[j]	=ESgt[j]*EZcr2gp;
  }
  iK	=41;
  ba	=1./(iK-1);
  ra	=0.;
  pF	=0.;
  for(k=1; k < iK; k++){
    b	=ba*k;
    za	=ra;
    ra	=0.;
    for(j=0; j < ESNp; j++){
      ra	+=sqrt(1.-b*pB[j])*D[j];
    }
    ra	=b*D0/ra;
    pF	+=EZcr2*ba*(ra+za);
  }
  *frac	=1.-0.75*aph*pF;
  return(0);
}

double ESjbs[128]={128*0.};

int ESBootstrap1()
{
  int iboot,i,j,ierr,nout;
  double zeff,rm2avx,epsx,rmajx,qx,ftrapx,
  pvale,ppvale,tempxe,tempxpe,pvali,ppvali,tempxi,tempxpi;
  double amu0=1.256637061e-6,echarge=1.6022e-19,emass=9.1095e-31
    ,pmass=1.6726e-27,xeps0=8.854e-12;
  double zave,xf,alfi,denom,anum1,anum2,tepote,tipoti,xmi,xclog,xdense,xdensi,
  xvthe,xvthi,xtaue,xtaui,xnuse,xnusi,xl31,xl32,col1,col2,alfhh;
  double a13,b13,c13,a23,b23,c23,s,ss;
  double jbs;
  int n1;
  double x[2],Y[65],YY0[65],YY[65],fr[65];
  double A[65]={
    0.0000e+00,
    1.2559e-01,
    1.7571e-01,
    2.1396e-01,
    2.4601e-01,
    2.7399e-01,
    2.9908e-01,
    3.2200e-01,
    3.4308e-01,
    3.6263e-01,
    3.8100e-01,
    3.9833e-01,
    4.1469e-01,
    4.3017e-01,
    4.4503e-01,
    4.5917e-01,
    4.7274e-01,
    4.8571e-01,
    4.9826e-01,
    5.1037e-01,
    5.2206e-01,
    5.3335e-01,
    5.4435e-01,
    5.5506e-01,
    5.6544e-01,
    5.7564e-01,
    5.8553e-01,
    5.9526e-01,
    6.0477e-01,
    6.1414e-01,
    6.2336e-01,
    6.3246e-01,
    6.4142e-01,
    6.5031e-01,
    6.5914e-01,
    6.6791e-01,
    6.7661e-01,
    6.8531e-01,
    6.9401e-01,
    7.0270e-01,
    7.1144e-01,
    7.2017e-01,
    7.2899e-01,
    7.3786e-01,
    7.4686e-01,
    7.5595e-01,
    7.6512e-01,
    7.7450e-01,
    7.8401e-01,
    7.9372e-01,
    8.0362e-01,
    8.1381e-01,
    8.2422e-01,
    8.3497e-01,
    8.4611e-01,
    8.5764e-01,
    8.6936e-01,
    8.8214e-01,
    8.9515e-01,
    9.0877e-01,
    9.2376e-01,
    9.3965e-01,
    9.5680e-01,
    9.7623e-01,
    1.0000e+00
  };
  double Ti2[65],Ti[65]={
    2.1510e+01,
    2.1600e+01,
    2.1630e+01,
    2.1640e+01,
    2.1610e+01,
    2.1560e+01,
    2.1490e+01,
    2.1400e+01,
    2.1300e+01,
    2.1170e+01,
    2.1020e+01,
    2.0860e+01,
    2.0690e+01,
    2.0500e+01,
    2.0290e+01,
    2.0070e+01,
    1.9830e+01,
    1.9590e+01,
    1.9330e+01,
    1.9050e+01,
    1.8760e+01,
    1.8460e+01,
    1.8150e+01,
    1.7830e+01,
    1.7490e+01,
    1.7150e+01,
    1.6790e+01,
    1.6420e+01,
    1.6040e+01,
    1.5650e+01,
    1.5260e+01,
    1.4850e+01,
    1.4430e+01,
    1.4000e+01,
    1.3570e+01,
    1.3120e+01,
    1.2670e+01,
    1.2220e+01,
    1.1750e+01,
    1.1280e+01,
    1.0810e+01,
    1.0330e+01,
    9.8450e+00,
    9.3580e+00,
    8.8690e+00,
    8.3790e+00,
    7.8880e+00,
    7.3970e+00,
    6.9070e+00,
    6.4190e+00,
    5.9330e+00,
    5.4500e+00,
    4.9700e+00,
    4.4940e+00,
    4.0230e+00,
    3.5560e+00,
    3.0930e+00,
    2.6350e+00,
    2.1830e+00,
    1.7380e+00,
    1.3050e+00,
    8.8930e-01,
    5.0730e-01,
    1.8700e-01,
    4.6500e-03
  };
  int n;
  double r1,EZr2,z,EZdra,EZdza,drp,dzp,a,gt,Lc,Vc;
  extern double ESgaP[];

  n	=64;
  n1	=65;
  z	=0.;
  EZf2spl(Ti,Ti2,&z,NULL,Ti,NULL,A,&n,ESgaP,ESgaP+1,ESgaP+2,ESgaP+3);
  ierr	=0;
  nout	=6;

  iboot	=2;

  xmi	=2.;
  xclog	=17.0;

  Y[0]	=0.;
  YY[0]	=0.;
  YY0[0]=0.;
  ESjbs[0]	=0.;
  i	=0;
  zeff	=1.;
  zave	=1.;
  for(i=1; i < 65; i++){
    a	=A[i];
    ESSetSplA(a);

    splRA(&Lc,NULL,ESLc,ESLc2);
    splRA(&Vc,NULL,ESVc,ESVc2);
    rm2avx	=Lc/((Lc+Vc)*ESaR0[0]*ESaR0[0]);
    gt	=0.;
    ES2DMapFlx2Lab(&EZr2,&z,&EZdra,&EZdza,&drp,&dzp,a,gt);
    gt	=EZcgp;
    ES2DMapFlx2Lab(&r1,&z,&EZdra,&EZdza,&drp,&dzp,a,gt);
    rmajx	=EZcr2*(EZr2+r1);
    epsx	=EZcr2*(EZr2-r1)/rmajx;
    splRA(&qx,NULL,ESsq,ESsq2a);
    splRA(&z,NULL,ESsp,ESsp2a);


    s	=1./zave;
    pvale	=z/(1.+s);
    pvali	=z*s/(1.+s);
    if(pvale == 0.){
      pvale	=1e-6;
    }
    if(pvali == 0.){
      pvali	=1e-6;
    }
    splRA(&z,NULL,ESPs,ESPs2a);
    ppvale	=-z/(1.+s);
    ppvali	=-z*s/(1.+s);
    tempxe	=Ti[i];
    tempxi	=Ti[i];

    tempxpe	=0.;
    tempxpi	=0.;
    splr0(&tempxe,&tempxpe,&a,Ti,Ti2,A,&n);
    splRA(&z,NULL,ESdgY,ESdgY2a);
    z	*=a;
    tempxpe	/=-z;
    tempxpi	=tempxpe;

    ESFractionOfTrapped1(&ftrapx,a);

    fr[i]	=ftrapx;
    xf	=ftrapx/(1.-ftrapx);
    alfi=1.172/(1.0+0.462*xf);

    denom=(1.414+zeff)*zeff+xf*(0.754+(2.657+2.*zeff)*zeff
				+xf*(0.348+(1.243+zeff)*zeff));
    anum1	=xf*(0.754+(2.21+zeff)*zeff+xf*(0.348+(1.24+zeff)*zeff));
    anum2	=xf*(0.884+2.07*zeff);

    tepote	=tempxpe/tempxe;
    tipoti	=tempxpi/tempxi;

    YY0[i]	=(anum1*(ESPs[i]+pvali*alfi*tipoti)+anum2*pvale*tepote)
      /(denom*rm2avx);
#ifdef H
    printf("%2d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e\n"
	   ,i,sqrt(epsx),anum1/(xf*denom),anum1/(denom*sqrt(epsx))
	   ,1./(ESaR0[0]*rm2avx),1./(ESaR0[0]*ESaR0[0]*rm2avx)
	   ,sqrt(epsx)*ESPs[i]*EZcrgm0*ESaR0[0],YY0[i]*EZcrgm0/ESaR0[0]);
    printf("%2d frac=%10.3e x=%10.3e\n",i,ftrapx,xf);
    YY0[i]	=EZcrgm0*anum1*ESPs[i]/(denom*rm2avx*ESaR0[0]);
#endif
    if(iboot == 1){
      jbs=-(anum1*(ppvale+ppvali-pvali*alfi*tipoti)-anum2*pvale*tepote)
	/(denom*rm2avx);
      return(0);
    }
#ifdef H
    zave	=pvali != 0. ? (pvale/tempxe)/(pvali/tempxi) : 0.;
    zave	=tempxi/tempxe;
#endif

    if(zeff <= 2.0){
      a13	=1.02-0.13*(zeff-1.0);
      b13	=1.07-0.45*(zeff-1.0);
      c13	=1.07-0.38*(zeff-1.0);
      a23	=0.57-0.05*(zeff-1.0);
      b23	=0.61-0.27*(zeff-1.0);
      c23	=0.61-0.23*(zeff-1.0);
    }
    else{
      a13	=0.89-0.05*(zeff-2.0);
      b13	=0.62-0.03*(zeff-2.0);
      c13	=0.69-0.09*(zeff-2.0);
      a23	=0.52-0.02*(zeff-2.0);
      b23	=0.34-0.005*(zeff-2.0);
      b23	=0.34-0.005*(zeff-2.0);
      c23	=0.38-0.05*(zeff-2.0);
    }
    xdense	=pvale/(amu0*tempxe*echarge*1.e3);
    xdensi	=pvali/(amu0*tempxi*echarge*1.e3);
    xvthe	=sqrt(2.0*echarge*1.e3*tempxe/emass);
    xvthi	=sqrt(2.0*echarge*1.e3*tempxi/(xmi*pmass));

    s	=xeps0*emass*xvthe;
    ss	=echarge*echarge;
    xtaue=(12.0*EZcgp*sqrt(EZcgp)*s*s*xvthe)/(4.0*xdense*zeff*ss*ss*xclog);
    s	=xeps0*xmi*pmass*xvthi;
    ss	=zave*echarge;
    ss	*=ss;
    xtaui=1.414*(12.0*EZcgp*sqrt(EZcgp)*s*s*xvthi)/(4.0*xdensi*ss*ss*xclog);
    s	=epsx*sqrt(epsx);
    xnuse=1.414*(rmajx*qx)/(xtaue*xvthe*s);
    xnusi=1.414*(rmajx*qx)/(xtaui*xvthi*s);
    xl31=anum1/denom;
    xl32=anum2/denom;
    col1=1./((1.+a13*sqrt(xnuse)+b13*xnuse)*(1.+c13*xnuse*s));
    col2=1./((1.+a23*sqrt(xnuse)+b23*xnuse)*(1.+c23*xnuse*s));

    s	*=s;
    ss	=s*xnusi*xnusi;
    s	=1.+ss;
    alfhh=((alfi-0.35*sqrt(xnusi))/(1.+0.7*sqrt(xnusi))-2.1*ss)
      /((1.+ss)*(1.+xnuse*xnuse*s));
    ESjbs[i]=-(xl31*col1*(ppvale+ppvali-alfhh*pvali*tipoti)
	       -pvale*(2.5*xl31*col1-(2.5*xl31-xl32)*col2)*tepote)/rm2avx;
    Y[i]	=EZcrgm0*ESjbs[i];
    ESjbs[i]	/=ESaR0[0];
    YY0[i]	*=EZcrgm0/ESaR0[0];
    EZout("sddd","Te",tempxe,tempxpe,Y[i]);
  }
  return(0);
}

int ESBootstrap()
{
  int iboot,i,j,ierr,nout;
  double zeff,rm2avx,epsx,rmajx,qx,ftrapx,
  pvale,ppvale,tempxe,tempxpe,pvali,ppvali,tempxi,tempxpi;
  double amu0=1.256637061e-6,echarge=1.6022e-19,emass=9.1095e-31
    ,pmass=1.6726e-27,xeps0=8.854e-12;
  double zave,xf,alfi,denom,anum1,anum2,tepote,tipoti,xmi,xclog,xdense,xdensi,
  xvthe,xvthi,xtaue,xtaui,xnuse,xnusi,xl31,xl32,col1,col2,alfhh;
  double a13,b13,c13,a23,b23,c23,s,ss;
  double jbs;
  double x[2],Y[ESNa1],YY0[ESNa1],YY[ESNa1],fr[ESNa1];

  ierr	=0;
  nout	=6;

  iboot	=2;

  xmi	=2.;
  xclog	=17.0;

#ifdef XWIN
  x[0]	=0.;
  x[1]	=1.;
  SetPlotName("#","j_Bootstrap [MA]","Bootstrap");
  Scale2d(2,x,x,2,6,2);
  Y[0]	=0.;
  YY[0]	=0.;
  YY0[0]=0.;
#endif
  ESjbs[0]	=0.;
  i	=0;
#ifdef H
  printf("%2s %11s %11s %11s %11s\n","i"
	 ,"a","jb MA/m^2","jBS MA/m^2","jBS0 MA/m^2");
  printf("%2d %11.4e %11.4e %11.4e %11.4e\n",0,0.,ESjb[0]*EZcrgm0,0.,0.);
#endif
  for(i=1; i < ESNa1; i++){
    zeff	=1.2;
    zeff	=1.;
    zeff	=1.6940;
    zave	=1.21;

    zeff	=1.;
    zave	=1.;

    rm2avx	=ESLc[i]/((ESLc[i]+ESVc[i])*ESaR0[0]*ESaR0[0]);
    j		=ESNp1*i;
    rmajx	=EZcr2*(ESsr[j]+ESsr[j+ESNp/2]);
    epsx	=EZcr2*(ESsr[j]-ESsr[j+ESNp/2])/rmajx;
    qx		=ESsq[i];

    s	=1./zave;
    pvale	=ESsp[i]/(1.+s);
    pvali	=ESsp[i]*s/(1.+s);
    if(pvale == 0.){
      pvale	=1e-6;
    }
    if(pvali == 0.){
      pvali	=1e-6;
    }
    ppvale	=-ESPs[i]/(1.+s);
    ppvali	=-ESPs[i]*s/(1.+s);

    tempxe	=15.;
    tempxi	=15.;
    tempxpe	=0.;
    tempxpi	=0.;
    ESFractionOfTrapped(&ftrapx,i);
#ifdef H
    bootstrappp_(YY+i,&iboot, &j, &nout,
		 &zeff,&rm2avx,&epsx,&rmajx,&qx,&ftrapx,
		 &pvale,&ppvale,&tempxe,&tempxpe,
		 &pvali,&ppvali,&tempxi,&tempxpi);
    YY[i]	*=EZcrgm0/ESaR0[0];
#endif

    YY[i]	=EZcrgm0*1.46*sqrt(epsx)*ESaR0[0]*ESPs[i];

    fr[i]	=ftrapx;
    xf	=ftrapx/(1.-ftrapx);
    alfi=1.172/(1.0+0.462*xf);

    denom=(1.414+zeff)*zeff+xf*(0.754+(2.657+2.*zeff)*zeff
				+xf*(0.348+(1.243+zeff)*zeff));
    anum1	=xf*(0.754+(2.21+zeff)*zeff+xf*(0.348+(1.24+zeff)*zeff));
    anum2	=xf*(0.884+2.07*zeff);

    tepote	=tempxpe/tempxe;
    tipoti	=tempxpi/tempxi;

    YY0[i]	=(anum1*(ESPs[i]+pvali*alfi*tipoti)+anum2*pvale*tepote)
      /(denom*rm2avx);
#ifdef H
    printf("%2d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e\n"
	   ,i,sqrt(epsx),anum1/(xf*denom),anum1/(denom*sqrt(epsx))
	   ,1./(ESaR0[0]*rm2avx),1./(ESaR0[0]*ESaR0[0]*rm2avx)
	   ,1.46*sqrt(epsx)*ESPs[i]*EZcrgm0*ESaR0[0],YY0[i]*EZcrgm0/ESaR0[0]);
    printf("%2d frac=%10.3e x=%10.3e\n",i,ftrapx,xf);
    YY0[i]	=EZcrgm0*anum1*ESPs[i]/(denom*rm2avx*ESaR0[0]);
#endif
    if(iboot == 1){
      jbs=-(anum1*(ppvale+ppvali-pvali*alfi*tipoti)-anum2*pvale*tepote)
	/(denom*rm2avx);
      return(0);
    }
#ifdef H
    zave	=pvali != 0. ? (pvale/tempxe)/(pvali/tempxi) : 0.;
    zave	=tempxi/tempxe;
#endif

    if(zeff <= 2.0){
      a13	=1.02-0.13*(zeff-1.0);
      b13	=1.07-0.45*(zeff-1.0);
      c13	=1.07-0.38*(zeff-1.0);
      a23	=0.57-0.05*(zeff-1.0);
      b23	=0.61-0.27*(zeff-1.0);
      c23	=0.61-0.23*(zeff-1.0);
    }
    else{
      a13	=0.89-0.05*(zeff-2.0);
      b13	=0.62-0.03*(zeff-2.0);
      c13	=0.69-0.09*(zeff-2.0);
      a23	=0.52-0.02*(zeff-2.0);
      b23	=0.34-0.005*(zeff-2.0);
      b23	=0.34-0.005*(zeff-2.0);
      c23	=0.38-0.05*(zeff-2.0);
    }
    xdense	=pvale/(amu0*tempxe*echarge*1.e3);
    xdensi	=pvali/(amu0*tempxi*echarge*1.e3);
    xvthe	=sqrt(2.0*echarge*1.e3*tempxe/emass);
    xvthi	=sqrt(2.0*echarge*1.e3*tempxi/(xmi*pmass));
    s	=xeps0*emass*xvthe;
    ss	=echarge*echarge;
    xtaue=(12.0*EZcgp*sqrt(EZcgp)*s*s*xvthe)/(4.0*xdense*zeff*ss*ss*xclog);
    s	=xeps0*xmi*pmass*xvthi;
    ss	=zave*echarge;
    ss	*=ss;
    xtaui=1.414*(12.0*EZcgp*sqrt(EZcgp)*s*s*xvthi)/(4.0*xdensi*ss*ss*xclog);
    s	=epsx*sqrt(epsx);
    xnuse=1.414*(rmajx*qx)/(xtaue*xvthe*s);
    xnusi=1.414*(rmajx*qx)/(xtaui*xvthi*s);
    xl31=anum1/denom;
    xl32=anum2/denom;
    col1=1./((1.+a13*sqrt(xnuse)+b13*xnuse)*(1.+c13*xnuse*s));
    col2=1./((1.+a23*sqrt(xnuse)+b23*xnuse)*(1.+c23*xnuse*s));

    s	*=s;
    ss	=s*xnusi*xnusi;
    s	=1.+ss;
    alfhh=((alfi-0.35*sqrt(xnusi))/(1.+0.7*sqrt(xnusi))-2.1*ss)
      /((1.+ss)*(1.+xnuse*xnuse*s));
    ESjbs[i]=-(xl31*col1*(ppvale+ppvali-alfhh*pvali*tipoti)
	       -pvale*(2.5*xl31*col1-(2.5*xl31-xl32)*col2)*tepote)/rm2avx;
    ESjbs[i]	/=ESaR0[0];
    Y[i]	=EZcrgm0*ESjbs[i];
    YY0[i]	*=EZcrgm0/ESaR0[0];

    s	=EZcrgm0*1.46*sqrt(epsx);
    ss	=s*ESPs[i]/(rm2avx*ESaR0[0]);
    ss	=s*ESPs[i]*ESaR0[0];


    ESjbs[i]	*=0.75;
#ifdef H
    printf("%2d %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n"
	   ,i,ESsa[i],ESjb[i]*EZcrgm0
	   ,Y[i],YY0[i],ss,1.46*sqrt(epsx)*denom/anum1);
#endif
  }
#ifdef XWIN
  Plot2d(2,ESsa,Y,ESNa1,6,0,0,0);
  Plot2d(2,ESsa,YY,ESNa1,6,0,14,0);
#endif
  return(0);
}

#ifdef H
      subroutine bootstrap(ajavbs, iboot, ierr, nout,
     1  zeff,rm2avx,epsx,rmajx,qx,ftrapx,
     2  pvale,ppvale,tempxe,tempxpe,
     3  pvali,ppvali,tempxi,tempxpi,xnuse,xnusi)
c
c     DESCRIPTION OF RETURNED PARAMETERS:
c     ajavbs   surface averaged parallel current due to bootstrap
c              divided by surface averaged B dot grad phi
c              (divide by mu0=4*pi*1.e-7 to get amperes/meter)
c     ierr     =0 for normal calculation
c              =1 for error exit
c
c     DESCRIPTION OF INPUT PARAMETERS:
c     iboot    =1 for collisionless bootstrap calculation
c              =2 to include collisonal corrections
c     nout     logical unit for error output
c
c     zeff     effective charge
c     rm2avx   surface average of 1/R**2
c     epsx     inverse aspect ratio of this flux surface
c     rmajx    major radius of this flux surface (m)
c     qx       safety factor
c     ftrapx   trapped particle fraction
c     pvale    electron pressure times mu0
c     ppvale   derivative of electron pressure wrt psi  (PF/radian)
c     pvali    ion pressure times mu0
c     ppvali   derivative of ion pressure wrt psi       (PF/radian)
c     tempxe   electron temperature (kev)
c     tempxpe  derivative of electron pressure wrt psi  (PF/radian)
c     tempxi   ion temperature      (kev)
c     tempxpi  derivative of ion temperature wrt psi    (PF/radian)
c
c
      data amu0/1.256637061e-6/
      data pi/3.1415926535/
      data echarge/1.6022e-19/
      data emass/9.1095e-31/
      data pmass/1.6726e-27/
      data xeps0/8.854e-12/
c
      ierr = 0
      if(tempxe.le.0 .or. tempxi.le.0) go to 101
      if(ftrapx.le.0 .or. ftrapx.ge.1) go to 102
      if(pvale .le.0 .or. pvali .le.0) go to 103
      if(rm2avx .le. 0) go to 104
      if(zeff .lt. 1.0) go to 105
      if(iabs(iboot).lt.1 .or. iabs(iboot).gt.2) go to 106
      zave = (pvale/tempxe)/(pvali/tempxi)
      xf=ftrapx/(1.0-ftrapx)
      alfi = 1.172/(1.0+0.462*xf)
      denom = 1.414*zeff + zeff**2 + xf*(0.754+2.657*zeff+2.*zeff**2)
     1                          + xf**2*(0.348+1.243*zeff+   zeff**2)
      anum1 = xf*(0.754 + 2.21*zeff + zeff**2
     1      + xf*(0.348 + 1.24*zeff + zeff**2) )
      anum2 = xf*(0.884 + 2.07*zeff)
      pepope=ppvale/pvale
      pipopi=ppvali/pvali
      tepote=tempxpe/tempxe
      tipoti=tempxpi/tempxi
      if(abs(iboot) .eq. 1) then
      ajavbs=-(pvale/rm2avx)*(anum1*(pepope+(pvali/pvale)
     1          *(pipopi - alfi*tipoti) ) - anum2*tepote) / denom
      return
      endif
c
      if(zeff .le. 2.0) then
      a13=1.02-0.13*(zeff-1.0)
      b13=1.07-0.45*(zeff-1.0)
      c13=1.07-0.38*(zeff-1.0)
      a23=0.57-0.05*(zeff-1.0)
      b23=0.61-0.27*(zeff-1.0)
      c23=0.61-0.23*(zeff-1.0)
      else
      a13=0.89-0.05*(zeff-2.0)
      b13=0.62-0.03*(zeff-2.0)
      c13=0.69-0.09*(zeff-2.0)
      a23=0.52-0.02*(zeff-2.0)
      b23=0.34-0.005*(zeff-2.0)
      c23=0.38-0.05*(zeff-2.0)
      endif
      xmi=1.0
      xclog=17.0
      xdense=pvale/(amu0*tempxe*echarge*1.e3)
      xdensi=pvali/(amu0*tempxi*echarge*1.e3)
      xvthe=sqrt(2.0*echarge*1.e3*tempxe/emass)
      xvthi=sqrt(2.0*echarge*1.e3*tempxi/(xmi*pmass))
      xtaue=(12.0*pi**1.5*xeps0**2*(emass)**2*xvthe**3)/
     +(4.0*xdense*zeff*(echarge)**4*xclog)
      xtaui=1.414*(12.0*pi**1.5*xeps0**2*(xmi*pmass)**2*
     +xvthi**3)/(4.0*xdensi*(zave*echarge)**4*xclog)
      xnuse=1.414*(rmajx*qx)/(xtaue*xvthe*epsx**1.5)
      xnusi=1.414*(rmajx*qx)/(xtaui*xvthi*epsx**1.5)
      xl31=anum1/denom
      xl32=anum2/denom
      col1=(1.0/(1.0+a13*sqrt(xnuse)+b13*xnuse))*
     +(1.0/(1.0+c13*xnuse*epsx**1.5))
      col2=(1.0/(1.0+a23*sqrt(xnuse)+b23*xnuse))*
     +(1.0/(1.0+c23*xnuse*epsx**1.5))
      alfhh=((alfi-0.35*sqrt(xnusi))/(1.0+0.7*sqrt(xnusi))
     +-2.1*xnusi**2*epsx**3)*(1.0/(1.0+xnusi**2*epsx**3))*
     +(1.0/(1.0+xnuse**2*epsx**3))
      ajavbs=-(pvale/rm2avx)*(xl31*col1*(pepope+(pvali/pvale)
     +*(pipopi-alfhh*tipoti))-(2.5*xl31*col1-(2.5*xl31-xl32)*col2)
     +*tepote)
c
      return
c
c.....error exit
  101 continue
      write(nout,1101) tempxe,tempxi
 1101 format(" Error in bootstrap, tempxe,tempxi=",1p2e12.4)
      go to 200
  102 continue
      write(nout,1102) ftrapx
 1102 format(" Error in bootstrap, ftrapx=",1pe12.4)
      go to 200
  103 continue
      write(nout,1103) pvale,pvali
 1103 format(" Error in bootstrap, pvale,pvali =",1p2e12.4)
      go to 200
  104 continue
      write(nout,1104) rm2avx
 1104 format(" Error in bootstrap, rm2avx=",1pe12.4)
      go to 200
  105 continue
      write(nout,1105) zeff
 1105 format(" Error in bootstrap, zeff=",1pe12.4)
      go to 200
  106 continue
      write(nout,1106) iboot
 1106 format(" Error in bootstrap, iboot=",i5)
c
  200 ierr = 1
      return
      end
#endif




#ifdef H
       subroutine bootstrap(ajavbs, iboot, ierr, nout,
      1  zeff,rm2avx,epsx,rmajx,qx,ftrapx,
      2  pvale,ppvale,tempxe,tempxpe,
      3  pvali,ppvali,tempxi,tempxpi )

c                                             6/27/94 scj
c
c     DESCRIPTION OF RETURNED PARAMETERS:
c     ajavbs   surface averaged J dot B due to bootstrap
c              divided by surface averaged B dot grad phi
c              (divide by mu0=4*pi*1.e-7 to get amperes/meter)
c     ierr     =0 for normal calculation
c              =1 for error exit
c
c     DESCRIPTION OF INPUT PARAMETERS:
c     iboot    =1 for collisionless bootstrap calculation
c              =2 to include collisonal corrections
c     nout     logical unit for error output
c
c     zeff     effective charge
c     rm2avx   surface average of 1/R**2
c     epsx     inverse aspect ratio of this flux surface
c     rmajx    major radius of this flux surface (m)
c     qx       safety factor
c     ftrapx   trapped particle fraction
c     pvale    electron pressure times mu0
c     ppvale   derivative of electron pressure wrt psi  (PF/radian)
c     pvali    ion pressure times mu0
c     ppvali   derivative of ion pressure wrt psi       (PF/radian)
c     tempxe   electron temperature (kev)
c     tempxpe  derivative of electron pressure wrt psi  (PF/radian)
c     tempxi   ion temperature      (kev)
c     tempxpi  derivative of ion temperature wrt psi    (PF/radian)
c
c
      data amu0/1.256637061e-6/
      data pi/3.1415926535/
      data echarge/1.6022e-19/
      data emass/9.1095e-31/
      data pmass/1.6726e-27/
      data xeps0/8.854e-12/
      data igoof/0/
c
      ierr = 0
      if(tempxe.le.0 .or. tempxi.le.0) go to 101
      if(ftrapx.le.0 .or. ftrapx.ge.1) go to 102
      if(pvale .le.0 .or. pvali .le.0) go to 103
      if(rm2avx .le. 0) go to 104
      if(zeff .lt. 0.9) go to 105
      if(iboot.lt.1 .or. iboot.gt.2) go to 106
      zave = (pvale/tempxe)/(pvali/tempxi)
      xf=ftrapx/(1.0-ftrapx)

      alfi = 1.172/(1.0+0.462*xf)

      denom = 1.414*zeff + zeff**2 + xf*(0.754+2.657*zeff+2.*zeff**2)
     1                          + xf**2*(0.348+1.243*zeff+   zeff**2)
      anum1 = xf*(0.754 + 2.21*zeff + zeff**2
     1      + xf*(0.348 + 1.24*zeff + zeff**2) )
      anum2 = xf*(0.884 + 2.07*zeff)
      pepope=ppvale/pvale
      pipopi=ppvali/pvali
      tepote=tempxpe/tempxe
      tipoti=tempxpi/tempxi
      if(abs(iboot) .eq. 1) then
      ajavbs=-(pvale/rm2avx)*(anum1*(pepope+(pvali/pvale)
     1          *(pipopi - alfi*tipoti) ) - anum2*tepote) / denom
      return
      endif
c
      if(zeff .le. 2.0) then
      a13=1.02-0.13*(zeff-1.0)
      b13=1.07-0.45*(zeff-1.0)
      c13=1.07-0.38*(zeff-1.0)
      a23=0.57-0.05*(zeff-1.0)
      b23=0.61-0.27*(zeff-1.0)
      c23=0.61-0.23*(zeff-1.0)
      else
      a13=0.89-0.05*(zeff-2.0)
      b13=0.62-0.03*(zeff-2.0)
      c13=0.69-0.09*(zeff-2.0)
      a23=0.52-0.02*(zeff-2.0)
      b23=0.34-0.005*(zeff-2.0)
      b23=0.34-0.005*(zeff-2.0)
      c23=0.38-0.05*(zeff-2.0)
      endif
      xmi=1.0
      xclog=17.0
      xdense=pvale/(amu0*tempxe*echarge*1.e3)
      xdensi=pvali/(amu0*tempxi*echarge*1.e3)
      xvthe=sqrt(2.0*echarge*1.e3*tempxe/emass)
      xvthi=sqrt(2.0*echarge*1.e3*tempxi/(xmi*pmass))
      xtaue=(12.0*pi**1.5*xeps0**2*(emass)**2*xvthe**3)/
     +(4.0*xdense*zeff*(echarge)**4*xclog)
      xtaui=1.414*(12.0*pi**1.5*xeps0**2*(xmi*pmass)**2*
     +xvthi**3)/(4.0*xdensi*(zave*echarge)**4*xclog)
      xnuse=1.414*(rmajx*qx)/(xtaue*xvthe*epsx**1.5)
      xnusi=1.414*(rmajx*qx)/(xtaui*xvthi*epsx**1.5)
      xl31=anum1/denom
      xl32=anum2/denom
      col1=(1.0/(1.0+a13*sqrt(xnuse)+b13*xnuse))*
     +(1.0/(1.0+c13*xnuse*epsx**1.5))
      col2=(1.0/(1.0+a23*sqrt(xnuse)+b23*xnuse))*
     +(1.0/(1.0+c23*xnuse*epsx**1.5))
      alfhh=((alfi-0.35*sqrt(xnusi))/(1.0+0.7*sqrt(xnusi))
     +-2.1*xnusi**2*epsx**3)*(1.0/(1.0+xnusi**2*epsx**3))*
     +(1.0/(1.0+xnuse**2*epsx**3))
      ajavbs=-(pvale/rm2avx)*(xl31*col1*(pepope+(pvali/pvale)
     +*(pipopi-alfhh*tipoti))-(2.5*xl31*col1-(2.5*xl31-xl32)*col2)
     +*tepote)
c
      return
c
c.....error exit
  101 continue
      write(nout,1101) tempxe,tempxi
 1101 format(" Error in bootstrap, tempxe,tempxi=",1p2e12.4)
      go to 200
  102 continue
      write(nout,1102) ftrapx
 1102 format(" Error in bootstrap, ftrapx=",1pe12.4)
      go to 200
  103 continue
      write(nout,1103) pvale,pvali
 1103 format(" Error in bootstrap, pvale,pvali =",1p2e12.4)
      go to 200
  104 continue
      write(nout,1104) rm2avx
 1104 format(" Error in bootstrap, rm2avx=",1pe12.4)
      go to 200
  105 continue
      write(nout,1105) zeff
 1105 format(" Error in bootstrap, zeff=",1pe12.4)
      go to 200
  106 continue
      write(nout,1106) iboot
 1106 format(" Error in bootstrap, iboot=",i5)
c
  200 igoof=igoof+1
      if(igoof.gt.10) ierr=1
      return
      end
#endif


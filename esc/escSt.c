#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fpreproc/f77name.h"

extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern int ESNa,ESNa1;
extern int ESNp,ESNp1,ESFp,ESFp1,ESnAP,ESnAF,ES2Mp,ES2Mp1;
extern int ESEqSolvFl,ESEqSolvInPr,ESEqSolvInCr;
extern double *ESsa,*ESpa,*ESgt;
extern double ESaRx,ESaZx,*ESaR0,*ESaZ0,*ESsr,*ESsz,*ESsb,*ESsb1a,*ESsb2a;
extern double *rcT,*rcT1a,*rcT2a,*rsT,*rsT1a,*rsT2a;
extern double ESRBt,ESRext,ESBt,ESgbext;
extern double *ESLc,*ESLc1,*ESLc2,*ESLs,*ESLs2,*ESVc
,*ESVc1,*ESVc2,*ESVs,*ESVs2,*ESDPc;
extern double *EScs1,*ESsn1,*EZcs2,*EZsn2,*EZdra,*EZdza,*EZdrgt,*EZdzgt;

extern double *ESgY0;

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


void blninit_(int*i,double*a,double*b);  
void bln1_();
void F77NAME(bln2)(double*a, double*b,double*c);
void F77NAME(bln3)();

static double *rdPV,*zdPV;

#define NPP 256
static int nP=128,nP1=128+1,nP2=128+2,nP5=128+5;
static int nP0=128;
static double Bala0=0.1,Bala1=1.,Balda=0.05;
static double R[NPP+1],Z[NPP+1];
static char Message[256];

static double gtmax=16.;

static double *Blnp,*Blnpp,*Blnq,*Blnqp,*Blnfg,*Blngp,*Blnf,*Blnfp
,*Blng22,*Blnp2r,*Blng11,*Blng12,*Blndg12t,*Blndp2ra,*Blndg22t,*Blndp2rt
,*Blnqg,*Blndqga,*Blngd,*Blndgda;

static int (*lBalSt)();

int ESPestBal(int *IsUnStable,double *A,int n1)
{
  int i,j,k,ki,kj;
  double gy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double sq,dsqgy;
  double *rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,EZdza,EZdzt,cs,sn,b,db,d2b;
  double *Gh,*Gh2,gh,t;
  double *d2rttc,*d2rtts,*d2ratc,*d2rats,*d2raac,*d2raas;
  double d2rtt,d2ztt,d2rat,d2zat,d2raa,d2zaa;
  double *Lc,*Ls,*dLc,*dLs,D0,dD0,D,dD,dg22,dg12,dgqt,gqt,L0;

  double dGh[65],dGh2[65],chhp,ch2hp;

  chhp	=EZcr2*ESgt[1];
  ch2hp	=chhp*EZcr6*ESgt[1];

  F77NAME(blninit)(&nP0,&gtmax,&ESRext);

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
    return(0);
  }
  rg	=1./ESRBt;
  rBt	=ESRext*rg;
  p2rBt	=rBt*rBt;
  dgt	=EZc2gp/nP;
  for(i	=0; i < n1; i++){
    a	=A[i];
    if(a < 1e-6){
      a	=1e-6;
    }
    ESSetSplA(a);
    if(ESEqSolvInPr != 2){
      splRA(Blnp,NULL,ESsp,ESsp2a);
    }
    else{
      ESSetSplDPr(a);
      splRDPr(Blnp,&ds,2);
    }
    *Blnp	*=p2rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Blnpp,&ds,0);
      *Blnpp	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Blnpp,&ds,1);
      break;
    default:
      splRA(Blnpp,NULL,ESPs,ESPs2a);
      break;
    }
    *Blnpp	*=-rBt;
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    dD0	=ds;
    rf	=-1./(dgya*a*rBt);
    splRA(&s,&ss,ESLc,ESLc2);
    sq		=-ESaR0[0]*dgya/s;
    *Blnf	=sq*rBt;

    *Blnfp	=(ESaR0[0]*ds+sq*ss)/(dgya*a*s);
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq		=1./s;
    dsqgy	=-ds*rf*sq*sq;
    *Blnq	=sq;
    *Blnqp	=dsqgy;
    splRA(Blnfg,&ds,ESaF,ESaF2a);
    *Blnfg	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    *Blngp	=-s/(ESRext*ESRext*(*Blnfg));
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
    rs[0]	+=ESaZ0[0];
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
      Z[j]	=rs[0]+b*sn;
      ss	=rc[0]+rc[1]*cs+rs[1]*sn;
      EZdra	=drc[0]+drc[1]*cs+drs[1]*sn;
      EZdrt	=rct[1]*sn+rst[1]*cs;
      ds	=dLc[1]*sn+dLs[1]*(1.-cs);
      gh	=Lc[1]*sn+Ls[1]*(1.-cs);
      gqt	=1.+Lc[1]*cs+Ls[1]*sn;
      dgqt	=dLc[1]*cs+dLs[1]*sn;
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
      R[j]	=ss;
      Blnp2r[j]		=ss*ss;
      Blndp2rt[j]	=-2.*ss*EZdrt;
      Blndp2ra[j]	=2.*ss*EZdra;

      D		=EZdra*EZdzt-EZdrt*EZdza;
      Blng22[j]	=(EZdrt*EZdrt+EZdzt*EZdzt)/(D*D);
      Blng11[j]	=(EZdra*EZdra+EZdza*EZdza)/(D*D);
      Blng12[j]	=(EZdra*EZdrt+EZdza*EZdzt)/(D*D);
      
      dD	=d2rat*EZdzt+EZdra*d2ztt-d2rtt*EZdza-EZdrt*d2zat;
      dg12	=d2rat*EZdrt+EZdra*d2rtt+d2zat*EZdzt+EZdza*d2ztt;
      dg22	=2.*(d2rtt*EZdrt+d2ztt*EZdzt);

      Blndg12t[j]	=(2.*Blng12[j]*dD-dg12/D)/D;
      Blndg22t[j]	=(2.*Blng22[j]*dD-dg22/D)/D;
      
      Blnqg[j]	=ss*D;
      Blndqga[j]=(Blndp2ra[j]*D0+Blnp2r[j]*dD0)*gqt+Blnp2r[j]*D0*dgqt;
      Blngd[j]	=-gh;
      Blndgda[j]=-sq*ds-dsqgy*gh;
    }
    R[j]	=R[0];
    Z[j]	=Z[0];
    Blng22[j]	=Blng22[0];
    Blnp2r[j]	=Blnp2r[0];
    Blng11[j]	=Blng11[0];
    Blng12[j]	=Blng12[0];
    Blndg12t[j]	=Blndg12t[0];
    Blndp2ra[j]	=Blndp2ra[0];
    Blndg22t[j]	=Blndg22t[0];
    Blndp2rt[j]	=Blndp2rt[0];
    Blnqg[j]	=Blnqg[0];
    Blndqga[j]	=Blndqga[0];
    Blngd[j]	=Blngd[0];
    Blndgda[j]	=Blndgda[0];
    j++;
    Blng22[j]	=Blng22[1];
    Blnp2r[j]	=Blnp2r[1];
    Blng11[j]	=Blng11[1];
    Blng12[j]	=Blng12[1];
    Blndg12t[j]	=Blndg12t[1];
    Blndp2ra[j]	=Blndp2ra[1];
    Blndg22t[j]	=Blndg22t[1];
    Blndp2rt[j]	=Blndp2rt[1];
    Blnqg[j]	=Blnqg[1];
    Blndqga[j]	=Blndqga[1];
    Blngd[j]	=Blngd[1];
    Blndgda[j]	=Blndgda[1];
    {
      int ifail;
      double di,gl;
      F77NAME(bln2z)(&ifail,&a,&di,&gl);
      IsUnStable[i]	=0;
      if(gl < 0.){
	IsUnStable[i]	|=1;
      }
      if(di > 0.){
	IsUnStable[i]	|=2;
      }
      if(ifail == 1){
	IsUnStable[i]	|=4;
      }
      if(ifail == 2){
	IsUnStable[i]	|=8;
      }
    }
  }
  free(Gh);
  F77NAME(bln3)();
  return(0);
}

void pestbal_(
	      int *IsUnStable  	/* IsUnStable[i] - integer array.
				   If IsUnStable[i] is
				   0 - surface is ballooning Stable;
				   Otherwise, first 4 bits are essential:
				   1 - surface is ballooning Untable;
				   2 - surface is Mercier Untable;
				   4 - large numerical error in the code;
				   8 - no convergence in the code.
				   */
	      ,double *a	/* Array of normalized minor radii, 
				   0 <= a[i] <= 1 */ 
	      ,int *n1		/* Number of test points in arrays 
				   IsUnStable[], a[] */
	      )
{
  nP	=256;
  if(nP > NPP){
    nP	=NPP;
  }
  if(nP < ESNp){
    nP	=ESNp;
  }
  nP1	=nP+1;
  nP2	=nP+2;
  nP5	=nP+5;
  nP0	=nP;
  ESPestBal(IsUnStable,a,*n1);
}

void EZpestbal(
	      int *IsUnStable  	/* IsUnStable[i] - integer array.
				   If IsUnStable[i] is
				   0 - surface is ballooning Stable;
				   Otherwise, first 4 bits are essential:
				   1 - surface is ballooning Untable;
				   2 - surface is Mercier Untable;
				   4 - large numerical error in the code;
				   8 - no convergence in the code.
				   */
	      ,double *a	/* Array of normalized minor radii, 
				   0 <= a[i] <= 1 */ 
	      ,int *n1		/* Number of test points in arrays 
				   IsUnStable[], a[] */
	      )
{
  pestbal_(IsUnStable,a,n1);
}

int ESPstBal()
{
  long int jj;
  int i,j,k,ki,kj;
  double gy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double sq,dsqgy;
  double *rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,EZdza,EZdzt,cs,sn,b,db,d2b;
  double *Gh,*Gh2,gh,t;
  double *d2rttc,*d2rtts,*d2ratc,*d2rats,*d2raac,*d2raas;
  double d2rtt,d2ztt,d2rat,d2zat,d2raa,d2zaa;
  double *Lc,*Ls,*dLc,*dLs,D0,dD0,D,dD,dg22,dg12,dgqt,gqt,L0;

  double dGh[65],dGh2[65],chhp,ch2hp;

  chhp	=EZcr2*ESgt[1];
  ch2hp	=chhp*EZcr6*ESgt[1];

  F77NAME(blninit)(&nP0,&gtmax,&ESRext);

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
    return(0);
  }
  rg	=1./ESRBt;
  rBt	=ESRext*rg;
  p2rBt	=rBt*rBt;
  dgt	=EZc2gp/nP;
  printf("Writing ESC metric tensor\n");
  gy	=Bala0;
  a	=sqrt(gy);

  for(a=Bala0; a <= Bala1+1e-6; a+=Balda){
#ifdef H
  while(gy <= Bala1){
    ss	=sqrt(gy);
    gy	+=Balda;
    if(fabs(gy-Bala1) < 0.01*Balda){
      gy	=Bala1;
    }
    do{
      ESSetSplA(a);
      splRA(&s,&dgya,ESqgY,ESqgY2a);
      a		+=(ss-s)/dgya;
    }while(fabs(ss-s) > 1e-8);
#endif
    ESSetSplA(a);
    if(ESEqSolvInPr != 2){
      splRA(Blnp,NULL,ESsp,ESsp2a);
    }
    else{
      ESSetSplDPr(a);
      splRDPr(Blnp,&ds,2);
    }
    *Blnp	*=p2rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Blnpp,&ds,0);
      *Blnpp	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Blnpp,&ds,1);
      break;
    default:
      splRA(Blnpp,NULL,ESPs,ESPs2a);
      break;
    }
    *Blnpp	*=-rBt;
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    dD0	=ds;
    rf	=-1./(dgya*a*rBt);
    splRA(&s,&ss,ESLc,ESLc2);
    sq		=-ESaR0[0]*dgya/s;
    *Blnf	=sq*rBt;

    *Blnfp	=(ESaR0[0]*ds+sq*ss)/(dgya*a*s);
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq		=1./s;
    dsqgy	=-ds*rf*sq*sq;
    *Blnq	=sq;
    *Blnqp	=dsqgy;
    splRA(Blnfg,&ds,ESaF,ESaF2a);
    *Blnfg	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    *Blngp	=-s/(ESRext*ESRext*(*Blnfg));
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
    rs[0]	+=ESaZ0[0];
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
      Z[j]	=rs[0]+b*sn;
      ss	=rc[0]+rc[1]*cs+rs[1]*sn;
      EZdra	=drc[0]+drc[1]*cs+drs[1]*sn;
      EZdrt	=rct[1]*sn+rst[1]*cs;
      ds	=dLc[1]*sn+dLs[1]*(1.-cs);
      gh	=Lc[1]*sn+Ls[1]*(1.-cs);
      gqt	=1.+Lc[1]*cs+Ls[1]*sn;
      dgqt	=dLc[1]*cs+dLs[1]*sn;
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
      R[j]	=ss;
      Blnp2r[j]		=ss*ss;
      Blndp2rt[j]	=-2.*ss*EZdrt;
      Blndp2ra[j]	=2.*ss*EZdra;

      D		=EZdra*EZdzt-EZdrt*EZdza;
      Blng22[j]	=(EZdrt*EZdrt+EZdzt*EZdzt)/(D*D);
      Blng11[j]	=(EZdra*EZdra+EZdza*EZdza)/(D*D);
      Blng12[j]	=(EZdra*EZdrt+EZdza*EZdzt)/(D*D);
      
      dD	=d2rat*EZdzt+EZdra*d2ztt-d2rtt*EZdza-EZdrt*d2zat;
      dg12	=d2rat*EZdrt+EZdra*d2rtt+d2zat*EZdzt+EZdza*d2ztt;
      dg22	=2.*(d2rtt*EZdrt+d2ztt*EZdzt);

      Blndg12t[j]	=(2.*Blng12[j]*dD-dg12/D)/D;
      Blndg22t[j]	=(2.*Blng22[j]*dD-dg22/D)/D;
      
      Blnqg[j]	=ss*D;
#ifdef H
      Blndqga[j]	=EZdra*D+ss*(d2raa*d2zat+EZdra*d2zat-d2rat*EZdza-EZdrt*d2zaa)
	-Blnqg[j]*d2gya;
      Blnqg[j]	=Blnp2r[j]*D0*gqt;
#endif
      Blndqga[j]=(Blndp2ra[j]*D0+Blnp2r[j]*dD0)*gqt+Blnp2r[j]*D0*dgqt;
#ifdef H
      ESSetSplP(t);
      splRP(&s,NULL,Gh,Gh2);
      Blngd[j]	=-s;
      Blndgda[j]=-ds;
#endif
      Blngd[j]	=-gh;
      Blndgda[j]=-sq*ds-dsqgy*gh;
    }
    R[j]	=R[0];
    Z[j]	=Z[0];
    Blng22[j]	=Blng22[0];
    Blnp2r[j]	=Blnp2r[0];
    Blng11[j]	=Blng11[0];
    Blng12[j]	=Blng12[0];
    Blndg12t[j]	=Blndg12t[0];
    Blndp2ra[j]	=Blndp2ra[0];
    Blndg22t[j]	=Blndg22t[0];
    Blndp2rt[j]	=Blndp2rt[0];
    Blnqg[j]	=Blnqg[0];
    Blndqga[j]	=Blndqga[0];
    Blngd[j]	=Blngd[0];
    Blndgda[j]	=Blndgda[0];
    j++;
    Blng22[j]	=Blng22[1];
    Blnp2r[j]	=Blnp2r[1];
    Blng11[j]	=Blng11[1];
    Blng12[j]	=Blng12[1];
    Blndg12t[j]	=Blndg12t[1];
    Blndp2ra[j]	=Blndp2ra[1];
    Blndg22t[j]	=Blndg22t[1];
    Blndp2rt[j]	=Blndp2rt[1];
    Blnqg[j]	=Blnqg[1];
    Blndqga[j]	=Blndqga[1];
    Blngd[j]	=Blngd[1];
    Blndgda[j]	=Blndgda[1];
    jj	=i+1;
    {
      int ic;
      double di,gl;
      F77NAME(bln2)(&a,&di,&gl);
#ifdef XWIN
      ic	=9;
      if(gl < 0.){
	ic	=13;
	printf("a=%10.3e sq=%10.3e\n",a,sq);
      }
      if(di > 0.){
	ic	=14;
      }
      ZColor(ic);
      ZPlotPolyLine(R,Z,nP1);
#endif
    }
  }
  free(Gh);
  F77NAME(bln3)();
  return(0);
}

void F77NAME(blnaddr)(double*p,double*pp,double*q,double*qp,double*fg,double*gp
	      ,double*f,double*fp
	      ,double *fg22,double *fp2r,double *fg11
	      ,double *fg12,double *fdg12t
	      ,double *fdp2ra,double *fdg22t,double *fdp2rt
	      ,double *fqg,double *fdqga
	      ,double *fgd,double *fdgda)
{
  Blnp	=p;
  Blnpp	=pp;
  Blnq	=q;
  Blnqp	=qp;
  Blnfg	=fg;
  Blngp	=gp;
  Blnf	=f;
  Blnfp	=fp;
  Blng22	=fg22;
  Blnp2r	=fp2r;
  Blng11	=fg11;
  Blng12	=fg12;
  Blndg12t	=fdg12t;
  Blndp2ra	=fdp2ra;
  Blndg22t	=fdg22t;
  Blndp2rt	=fdp2rt;
  Blnqg		=fqg;
  Blndqga	=fdqga;
  Blngd		=fgd;
  Blndgda	=fdgda;
  return;
}

/* Balloning with $\gz=\gf$ */
int ESPstBal0()
{
  long int jj;
  int i,j,k,ki,kj;
  double gy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double sq,dsqgy;
  double *rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,EZdza,EZdzt,cs,sn,b,db;
  double *Gh,*Gh2,gh,t;
  double *d2rttc,*d2rtts,*d2ratc,*d2rats,d2rtt,d2ztt,d2rat,d2zat;
  double *Lc,*Ls,*dLc,*dLs,D0,dD0,D,dD,dg22,dg12,dgqt,L0;

  double dGh[65],dGh2[65],chhp,ch2hp;

  F77NAME(blninit)(&nP0,&gtmax,&ESRext);

  chhp	=EZcr2*ESgt[1];
  ch2hp	=chhp*EZcr6*ESgt[1];

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
    return(0);
  }
  rg	=1./ESRBt;
  rBt	=ESRext*rg;
  p2rBt	=rBt*rBt;
  dgt	=EZc2gp/nP;
  printf("Ballooning calculations using PEST coordinates\n");
  gy	=Bala0;
  a	=sqrt(gy);
  while(gy <= Bala1){
    ss	=sqrt(gy);
    gy	+=Balda;
    if(fabs(gy-Bala1) < 0.01*Balda){
      gy	=Bala1;
    }
    do{
      ESSetSplA(a);
      splRA(&s,&dgya,ESqgY,ESqgY2a);
      a		+=(ss-s)/dgya;
    }while(fabs(ss-s) > 1e-8);
    ESSetSplA(a);
    if(ESEqSolvInPr != 2){
      splRA(Blnp,NULL,ESsp,ESsp2a);
    }
    else{
      ESSetSplDPr(a);
      splRDPr(Blnp,&ds,2);
    }
    *Blnp	*=p2rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Blnpp,&ds,0);
      *Blnpp	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Blnpp,&ds,1);
      break;
    default:
      splRA(Blnpp,NULL,ESPs,ESPs2a);
      break;
    }
    *Blnpp	*=-rBt;
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    dD0	=ds;
    rf	=-1./(dgya*a*rBt);
    splRA(&s,&ss,ESLc,ESLc2);
    sq		=-ESaR0[0]*dgya/s;
    *Blnf	=sq*rBt;

    *Blnfp	=(ESaR0[0]*ds+sq*ss)/(dgya*a*s);
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq		=1./s;
    dsqgy	=-ds*rf*sq*sq;
    *Blnq	=sq;
    *Blnqp	=dsqgy;
    splRA(Blnfg,&ds,ESaF,ESaF2a);
    *Blnfg	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    *Blngp	=-s/(ESRext*ESRext*(*Blnfg));

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
    rs[0]	+=ESaZ0[0];
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
      Z[j]	=rs[0]+b*sn;
      ss	=rc[0]+rc[1]*cs+rs[1]*sn;
      EZdra	=drc[0]+drc[1]*cs+drs[1]*sn;
      EZdrt	=rct[1]*sn+rst[1]*cs;
      ds	=dLc[1]*sn+dLs[1]*(1.-cs);
      dgqt	=dLc[1]*cs+dLs[1]*sn;
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
      R[j]	=ss;
      D		=EZdra*EZdzt-EZdrt*EZdza;
      Blng22[j]	=(EZdrt*EZdrt+EZdzt*EZdzt)/(D*D);
      dD	=d2rat*EZdzt+EZdra*d2ztt-d2rtt*EZdza-EZdrt*d2zat;
      dg12	=d2rat*EZdrt+EZdra*d2rtt+d2zat*EZdzt+EZdza*d2ztt;
      dg22	=2.*(d2rtt*EZdrt+d2ztt*EZdzt);
      s		=ss/(D*L0);
      Blndg12t[j]	=(dg22*ds/(D*D)+Blng22[j]*(dgqt-2.*ds*dD/D))*s
	-(dg12-(EZdra*EZdrt+EZdza*EZdzt)*(dD*ss+D*EZdrt)/(D*ss))/(D*D);

      Blndg22t[j]	=(2.*Blng22[j]*dD-dg22/D)*s/D;


      Blndg12t[j]	=-(dg12-(EZdra*EZdrt+EZdza*EZdzt)*(dD*ss+D*EZdrt)/(D*ss))/(D*D);


      Blndg22t[j]	=(2.*Blng22[j]*dD-dg22/D)*s/D;


      Blnp2r[j]	=ss*ss;
      Blndp2rt[j]	=-2.*ss*EZdrt*s;
      ds	*=s;
      EZdra	-=EZdrt*ds;
      EZdza	-=EZdzt*ds;
      Blng11[j]	=(EZdra*EZdra+EZdza*EZdza)*L0*L0/(ss*ss);
      Blng12[j]	=(EZdra*EZdrt+EZdza*EZdzt)*L0/(D*ss);
      Blndp2ra[j]	=2.*ss*EZdra;
      Blnqg[j]	=ss*ss*D0;
      Blndqga[j]	=Blndp2ra[j]*D0+Blnp2r[j]*dD0;
      Blngd[j]	=0.;
      Blndgda[j]	=0.;
    }
    R[j]	=R[0];
    Z[j]	=Z[0];
    Blng22[j]	=Blng22[0];
    Blnp2r[j]	=Blnp2r[0];
    Blng11[j]	=Blng11[0];
    Blng12[j]	=Blng12[0];
    Blndg12t[j]	=Blndg12t[0];
    Blndp2ra[j]	=Blndp2ra[0];
    Blndg22t[j]	=Blndg22t[0];
    Blndp2rt[j]	=Blndp2rt[0];
    Blnqg[j]	=Blnqg[0];
    Blndqga[j]	=Blndqga[0];
    Blngd[j]	=Blngd[0];
    Blndgda[j]	=Blndgda[0];
    j++;
    Blng22[j]	=Blng22[1];
    Blnp2r[j]	=Blnp2r[1];
    Blng11[j]	=Blng11[1];
    Blng12[j]	=Blng12[1];
    Blndg12t[j]	=Blndg12t[1];
    Blndp2ra[j]	=Blndp2ra[1];
    Blndg22t[j]	=Blndg22t[1];
    Blndp2rt[j]	=Blndp2rt[1];
    Blnqg[j]	=Blnqg[1];
    Blndqga[j]	=Blndqga[1];
    Blngd[j]	=Blngd[1];
    Blndgda[j]	=Blndgda[1];
    jj	=i+1;
    {
      int ic;
      double di,gl;
      F77NAME(bln2)(&a,&di,&gl);
#ifdef XWIN
      ic=9;
      if(gl < 0.){
	ic	=13;
	printf("a=%10.3e sq=%10.3e\n",a,sq);
      }
      if(di > 0.){
	ic	=14;
      }
      ZColor(ic);
      ZPlotPolyLine(R,Z,nP1);
#endif
    }
  }
  free(Gh);
  F77NAME(bln3)();
  return(0);
}

/* Ballooning with $g_{\gq\gq}=const$ */
int ESPstBal1()
{
  long int jj;
  int i,j,k,ki,kj;
  double gy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double sq,dsqgy;
  double *rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,EZdza,EZdzt,cs,sn,b,db;
  double *gQ,*gQ2,gq,gqa,gqt,gqat,gqtt,t;
  double *d2rttc,*d2rtts,*d2ratc,*d2rats,d2rtt,d2ztt,d2rat,d2zat;
  double *Lc,*Ls,*dLc,*dLs,D0,dD0,D,dD;
  double dg22,dg12,gh,ght,gha,ghat,ghtt,L0,G0;
  double *g22c,*g22c2,*g22s,*g22s2,*G22c,*dG22c,*G22s,*dG22s;


 F77NAME(blninit)(&nP0,&gtmax,&ESRext);

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
    ss	=-rc[1]*ESsn1[j]+rs[1]*EScs1[j];
    s	=b*EScs1[j];
    Blnp2r[j]	=ss;
    Blng22[j]	=sqrt(ss*ss+s*s);
  }
  Blng22[j]	=Blng22[0];
  ESgP2gF(g22c+i,g22s+i,Blng22,ESFp);

  ki		=2*ESNa1;
  rc[2]		=2.*rcT2a[ki];
  rs[2]		=2.*rsT2a[ki];
  for(j=0; j < ESNp; j++){
    EZdra	=-rc[2]*EZsn2[j]+rs[2]*EZcs2[j];
    Blng22[j]	=EZdra*Blnp2r[j]/Blng22[j];
  }
  Blng22[j]	=Blng22[0];
  ESgP2gF(g22c2+i,g22s2+i,Blng22,ESFp);
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
      Blnp2r[j]	=0.;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj		+=j;
	if(kj >= ESNp)
	  kj	-=ESNp;
	Blnp2r[j]	+=-rc[k]*ESsn1[kj]+rs[k]*EScs1[kj];
      }
      Blngd[j]	=b*EScs1[j];
      Blng22[j]	=sqrt(Blnp2r[j]*Blnp2r[j]+Blngd[j]*Blngd[j]);
    }
    Blng22[j]	=Blng22[0];
    ESgP2gF(g22c+i,g22s+i,Blng22,ESFp);
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
    Blng22[j]	=(Blnp2r[j]*EZdra+Blngd[j]*EZdza)/Blng22[j];
  }
  Blng22[j]	=Blng22[0];
  ESgP2gF(g22c2+i,g22s2+i,Blng22,ESFp);
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
  printf("Ballooning calculations using equidistant poloidal angle\n");
  gy	=Bala0;
  a	=sqrt(gy);
  while(gy <= Bala1){
    ss	=sqrt(gy);
    gy	+=Balda;
    if(fabs(gy-Bala1) < 0.01*Balda){
      gy	=Bala1;
    }
    do{
      ESSetSplA(a);
      splRA(&s,&dgya,ESqgY,ESqgY2a);
      a		+=(ss-s)/dgya;
    }while(fabs(ss-s) > 1e-8);
    ESSetSplA(a);
    if(ESEqSolvInPr != 2){
      splRA(Blnp,NULL,ESsp,ESsp2a);
    }
    else{
      ESSetSplDPr(a);
      splRDPr(Blnp,&ds,2);
    }
    *Blnp	*=p2rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Blnpp,&ds,0);
      *Blnpp	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Blnpp,&ds,1);
      break;
    default:
      splRA(Blnpp,NULL,ESPs,ESPs2a);
      break;
    }
    *Blnpp	*=-rBt;
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    dD0	=ds;
    rf	=-1./(dgya*a*rBt);
    splRA(&s,&ss,ESLc,ESLc2);
    sq		=-ESaR0[0]*dgya/s;
    *Blnf	=sq*rBt;

    *Blnfp	=(ESaR0[0]*ds+sq*ss)/(dgya*a*s);
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq		=1./s;
    dsqgy	=-ds*rf*sq*sq;
    *Blnq	=sq;
    *Blnqp	=dsqgy;
    splRA(Blnfg,&ds,ESaF,ESaF2a);
    *Blnfg	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    *Blngp	=-s/(ESRext*ESRext*(*Blnfg));

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
    rs[0]	+=ESaZ0[0];
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
      Z[j]	=rs[0]+b*sn;
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
      R[j]	=ss;
      Blngd[j]	=-gh;
      D		=EZdra*EZdzt-EZdrt*EZdza;
      gh	=1./D;
      dD	=(d2rat*EZdzt+EZdra*d2ztt-d2rtt*EZdza-EZdrt*d2zat)*gh;
#ifdef H
      gh	*=gh;
#endif
      gh	*=1./(ss*D0*ght);

      Blng22[j]	=(EZdrt*EZdrt+EZdzt*EZdzt)*gh;
      dg12	=(d2rat*EZdrt+EZdra*d2rtt+d2zat*EZdzt+EZdza*d2ztt)*gh;
      dg22	=2.*(d2rtt*EZdrt+d2ztt*EZdzt)*gh;

      gh	*=(EZdra*EZdrt+EZdza*EZdzt);

      Blng12[j]	=gh*gqt-Blng22[j]*gqa;

      s		=1./gqt;
      Blndg12t[j]=(2.*Blng12[j]*dD-dg12*gqt-gh*gqtt+dg22*gqa+Blng22[j]*gqat)*s;
      Blndg22t[j]=(2.*Blng22[j]*dD-dg22)*s;

      Blnp2r[j]	=ss*ss;
      Blndp2rt[j]	=-2.*ss*EZdrt*s;
      gqa	*=s;
      EZdra	-=EZdrt*gqa;
      EZdza	-=EZdzt*gqa;

      Blndgda[j]	=-sq*(gha-ght*gqa)+dsqgy*Blngd[j];

      D		*=s;
#ifdef H
      Blng11[j]	=(EZdra*EZdra+EZdza*EZdza)/(D*D);
#endif

      Blndp2ra[j]	=2.*ss*EZdra;

      Blnqg[j]	=Blnp2r[j]*D0*ght*s;
      Blng11[j]	=(EZdra*EZdra+EZdza*EZdza)*ss/(D*Blnqg[j]);

      Blndqga[j]	=(Blndp2ra[j]*D0+Blnp2r[j]*dD0)*ght*s
	+Blnp2r[j]*D0*(ghat*gqt-ght*gqat)*s*s
	  -Blnp2r[j]*D0*(ghtt*gqt-ght*gqtt)*s*s*gqa;
    }
    R[j]	=R[0];
    Z[j]	=Z[0];
    Blng22[j]	=Blng22[0];
    Blnp2r[j]	=Blnp2r[0];
    Blng11[j]	=Blng11[0];
    Blng12[j]	=Blng12[0];
    Blndg12t[j]	=Blndg12t[0];
    Blndp2ra[j]	=Blndp2ra[0];
    Blndg22t[j]	=Blndg22t[0];
    Blndp2rt[j]	=Blndp2rt[0];
    Blnqg[j]	=Blnqg[0];
    Blndqga[j]	=Blndqga[0];
    Blngd[j]	=Blngd[0];
    Blndgda[j]	=Blndgda[0];
    j++;
    Blng22[j]	=Blng22[1];
    Blnp2r[j]	=Blnp2r[1];
    Blng11[j]	=Blng11[1];
    Blng12[j]	=Blng12[1];
    Blndg12t[j]	=Blndg12t[1];
    Blndp2ra[j]	=Blndp2ra[1];
    Blndg22t[j]	=Blndg22t[1];
    Blndp2rt[j]	=Blndp2rt[1];
    Blnqg[j]	=Blnqg[1];
    Blndqga[j]	=Blndqga[1];
    Blngd[j]	=Blngd[1];
    Blndgda[j]	=Blndgda[1];
    jj	=i+1;
    {
      int ic;
      double di,gl;
      F77NAME(bln2)(&a,&di,&gl);
#ifdef XWIN
      ic=9;
      if(gl < 0.){
	ic	=13;
	printf("a=%10.3e sq=%10.3e\n",a,sq);
      }
      if(di > 0.){
	ic	=14;
      }
      ZColor(ic);
      ZPlotPolyLine(R,Z,nP1);
#endif
    }
  }
  free(gQ);
  F77NAME(bln3)();
  return(0);
}

/* Ballooning with Hamada coordinates */
int ESPstBal2()
{
  long int jj;
  int i,j,k,ki,kj;
  double gy,dgt,dgya,rf,rg,a,ss,s,ds,rBt,p2rBt;
  double sq,dsqgy;
  double *rc,*rs,*drc,*drs,*rct,*rst,EZdra,EZdrt,EZdza,EZdzt,cs,sn,b,db;
  double *gQ,*gQ2,gq,gqa,gqt,gqat,gqtt,t;
  double *d2rttc,*d2rtts,*d2ratc,*d2rats,d2rtt,d2ztt,d2rat,d2zat;
  double *Lc,*Ls,*dLc,*dLs,D0,dD0,D,dD;
  double dg22,dg12,gh,ght,gha,ghat,ghtt,L0,G0;
  double *rDc,*rDc2,*rDs,*rDs2,*Jc,*dJc,*Js,*dJs;

  F77NAME(blninit)(&nP0,&gtmax,&ESRext);

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
  dgt	=EZc2gp/nP;
  printf("Ballooning calculations using Hamada coordinates\n");
  gy	=Bala0;
  a	=sqrt(gy);
  while(gy <= Bala1){
    ss	=sqrt(gy);
    gy	+=Balda;
    if(fabs(gy-Bala1) < 0.01*Balda){
      gy	=Bala1;
    }
    do{
      ESSetSplA(a);
      splRA(&s,&dgya,ESqgY,ESqgY2a);
      a		+=(ss-s)/dgya;
    }while(fabs(ss-s) > 1e-8);
    ESSetSplA(a);
    if(ESEqSolvInPr != 2){
      splRA(Blnp,NULL,ESsp,ESsp2a);
    }
    else{
      ESSetSplDPr(a);
      splRDPr(Blnp,&ds,2);
    }
    *Blnp	*=p2rBt;
    switch(ESEqSolvInPr){
    case 0:
      ESSetSplDPr(a);
      splRDPr(Blnpp,&ds,0);
      *Blnpp	/=ESaR0[0];
      break;
    case 1:
      ESSetSplDPr(a);
      splRDPr(Blnpp,&ds,1);
      break;
    default:
      splRA(Blnpp,NULL,ESPs,ESPs2a);
      break;
    }
    *Blnpp	*=-rBt;
    splRA(&dgya,&ds,ESdgY,ESdgY2a);
    dD0	=ds;
    rf	=-1./(dgya*a*rBt);
    splRA(&s,&ss,ESLc,ESLc2);
    sq		=-ESaR0[0]*dgya/s;
    *Blnf	=sq*rBt;

    *Blnfp	=(ESaR0[0]*ds+sq*ss)/(dgya*a*s);
    if(ESEqSolvInCr != 6 && ESEqSolvInCr != 7){
      splRA(&s,&ds,ESgm,ESgm2a);
    }
    else{
      ESSetSplDCr(a);
      splRDCr(&s,&ds,7);
    }
    sq		=1./s;
    dsqgy	=-ds*rf*sq*sq;
    *Blnq	=sq;
    *Blnqp	=dsqgy;
    splRA(Blnfg,&ds,ESaF,ESaF2a);
    *Blnfg	*=rg;
    splRA(&s,NULL,ESaT,ESaT2a);
    s		*=rBt;
    *Blngp	=-s/(ESRext*ESRext*(*Blnfg));

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
    rs[0]	+=ESaZ0[0];
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
      Z[j]	=rs[0]+b*sn;
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
      R[j]	=ss;
      Blngd[j]	=-gh;
      D		=EZdra*EZdzt-EZdrt*EZdza;
      gh	=1./D;
      dD	=(d2rat*EZdzt+EZdra*d2ztt-d2rtt*EZdza-EZdrt*d2zat)*gh;
#ifdef H
      gh	*=gh;
#endif
      gh	*=ss/(D0*gqt);

      Blng22[j]	=(EZdrt*EZdrt+EZdzt*EZdzt)*gh;
      dg12	=(d2rat*EZdrt+EZdra*d2rtt+d2zat*EZdzt+EZdza*d2ztt)*gh;
      dg22	=2.*(d2rtt*EZdrt+d2ztt*EZdzt)*gh;

      gh	*=(EZdra*EZdrt+EZdza*EZdzt);

      Blng12[j]	=gh*gqt-Blng22[j]*gqa;

      s		=1./gqt;
      Blndg12t[j]=(2.*Blng12[j]*dD-dg12*gqt-gh*gqtt+dg22*gqa+Blng22[j]*gqat)*s;
      Blndg22t[j]=(2.*Blng22[j]*dD-dg22)*s;

      Blnp2r[j]	=ss*ss;
      Blndp2rt[j]	=-2.*ss*EZdrt*s;
      gqa	*=s;
      EZdra	-=EZdrt*gqa;
      EZdza	-=EZdzt*gqa;

      Blndgda[j]=-sq*(gha-ght*gqa)+dsqgy*Blngd[j];

      D		*=s;
#ifdef H
      Blng11[j]	=(EZdra*EZdra+EZdza*EZdza)/(D*D);
#endif
      Blng11[j]	=(EZdra*EZdra+EZdza*EZdza)*ss/(D*D0);
      Blndp2ra[j]	=2.*ss*EZdra;

      Blnqg[j]	=D0;
      Blndqga[j]	=dD0;
    }
    R[j]	=R[0];
    Z[j]	=Z[0];
    Blng22[j]	=Blng22[0];
    Blnp2r[j]	=Blnp2r[0];
    Blng11[j]	=Blng11[0];
    Blng12[j]	=Blng12[0];
    Blndg12t[j]	=Blndg12t[0];
    Blndp2ra[j]	=Blndp2ra[0];
    Blndg22t[j]	=Blndg22t[0];
    Blndp2rt[j]	=Blndp2rt[0];
    Blnqg[j]	=Blnqg[0];
    Blndqga[j]	=Blndqga[0];
    Blngd[j]	=Blngd[0];
    Blndgda[j]	=Blndgda[0];
    j++;
    Blng22[j]	=Blng22[1];
    Blnp2r[j]	=Blnp2r[1];
    Blng11[j]	=Blng11[1];
    Blng12[j]	=Blng12[1];
    Blndg12t[j]	=Blndg12t[1];
    Blndp2ra[j]	=Blndp2ra[1];
    Blndg22t[j]	=Blndg22t[1];
    Blndp2rt[j]	=Blndp2rt[1];
    Blnqg[j]	=Blnqg[1];
    Blndqga[j]	=Blndqga[1];
    Blngd[j]	=Blngd[1];
    Blndgda[j]	=Blndgda[1];
    jj	=i+1;
    {
      int ic;
      double di,gl;
      F77NAME(bln2)(&a,&di,&gl);
#ifdef XWIN
      ic=9;
      if(gl < 0.){
	ic	=13;
	printf("a=%10.3e sq=%10.3e\n",a,sq);
      }
      if(di > 0.){
	ic	=14;
      }
      ZColor(ic);
      ZPlotPolyLine(R,Z,nP1);
#endif
    }
  }
  free(gQ);
  F77NAME(bln3)();
  return(0);
}

#ifdef XWIN
void ESBalStPlot()
{
  int i,j;

#ifdef PFC
  IbPFIPlot();
#endif
  ZColor(lBalSt == NULL ? 0 : 3);
  j	=0;
  for(i=1; i < ESNa1; i++){
    j	+=ESNp1;
    ZPlotPolyLine(ESsr+j,ESsz+j,ESNp1);
  }
  if(lBalSt != NULL){
    lBalSt();
  }
  return;
}
#endif

int gsvDT(double *f,double *df,double T) /* from ctf.c */
{
  int i,ii,j;
  double a,b,aa,bb,h;
  static int n=41,n1=42;
  static double ccr6= 0.16666666666666;
  static double gsvT[42]={
    0.6,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0,
    5.5,
    6.0,
    6.5,
    7.0,
    7.5,
    8.0,
    8.5,
    9.0,
    9.5,
    10.0,
    15.0,
    20.0,
    25.0,
    30.0,
    35.0,
    40.0,
    45.0,
    50.0,
    55.0,
    60.0,
    65.0,
    70.0,
    75.0,
    80.0,
    85.0,
    90.0,
    95.0,
    100.0,
    150.0,
    200.0,
    250.0,
    300.0,
  };
  static double gsvS[42]={
    1.9564957495362972e-22,
    5.4834999999999999e-21,
    5.8917000000000000e-20,
    2.6275999999999999e-19,
    7.6096000000000001e-19,
    1.7128000000000000e-18,
    3.2740000000000001e-18,
    5.5842999999999996e-18,
    8.7602000000000006e-18,
    1.2891000000000000e-17,
    1.8039000000000001e-17,
    2.4237000000000001e-17,
    3.1494000000000002e-17,
    3.9796999999999999e-17,
    4.9110999999999997e-17,
    5.9390999999999995e-17,
    7.0575999999999998e-17,
    8.2598000000000004e-17,
    9.5384999999999996e-17,
    1.0886000000000000e-16,
    2.6541999999999998e-16,
    4.2434000000000001e-16,
    5.5939999999999998e-16,
    6.6532000000000003e-16,
    7.4479000000000002e-16,
    8.0251999999999997e-16,
    8.4311999999999998e-16,
    8.7050000000000004e-16,
    8.8774999999999998e-16,
    8.9729000000000003e-16,
    9.0096999999999993e-16,
    9.0022000000000006e-16,
    8.9612000000000005e-16,
    8.8950999999999993e-16,
    8.8103000000000004e-16,
    8.7117000000000001e-16,
    8.6033000000000002e-16,
    8.4879000000000003e-16,
    7.2769999999999997e-16,
    6.2782999999999999e-16,
    5.5246000000000001e-16,
    4.9541999999999997e-16,
  };
  static double gsvS2[42]={
    5.5247342629482765e-21,
    1.5196036259651044e-19,
    5.7229139449092497e-19,
    1.1687020594401500e-18,
    1.8174683677491811e-18,
    2.4487844695642135e-18,
    3.0120337539954294e-18,
    3.4814805144558559e-18,
    3.8364441881832558e-18,
    4.0903427328133810e-18,
    4.2149848805656846e-18,
    4.2497177449263956e-18,
    4.2021441397312814e-18,
    4.0457056961508908e-18,
    3.8790330756676228e-18,
    3.6221620011809201e-18,
    3.3523189196110037e-18,
    3.0565623203771328e-18,
    2.7814317988819713e-18,
    2.3297104840968106e-18,
    -1.6910624490066200e-19,
    -1.0868855044940964e-18,
    -1.2097517371235347e-18,
    -1.0677075470124458e-18,
    -8.6741807482733008e-19,
    -6.8022015367876496e-19,
    -5.2290131045800828e-19,
    -4.0097460448950663e-19,
    -3.0440027158423876e-19,
    -2.3182430917369521e-19,
    -1.7470249172115917e-19,
    -1.3256572394171652e-19,
    -9.9034612512089549e-20,
    -7.3695826010012773e-20,
    -5.4982083447849569e-20,
    -3.7575840198654130e-20,
    -2.9914555757549526e-20,
    -1.0765936771162730e-20,
    1.3020516472311116e-20,
    9.6118708819235198e-21,
    7.3320000000007093e-21,
    5.0521291180778995e-21
  };
  if(T < gsvT[0]){
    h		=pow(T,-1./3.);
    *f		=0.7*3.68e-12*exp(-19.94*h)*h*h;

    if(df != NULL){
      *df	=(*f)*(19.94*h-2.)/(3.*T);
    }
    return(0);
  }

  if(T > gsvT[n]) T=gsvT[n];

  ii	=n;
  i	=0;
  while(ii-i > 1){
    j	=(ii+i)/2;
    if(T >= gsvT[j]) i=j;
    else ii=j;
  }
  h	=gsvT[ii]-gsvT[i];
  a	=(gsvT[ii]-T)/h;
  aa	=a*a;
  b	=1.-a;
  bb	=b*b;
  *f	=a*gsvS[i]+b*gsvS[ii]
    +(a*(aa-1.)*gsvS2[i]+b*(bb-1.)*gsvS2[ii])*h*h*ccr6;
  if(df != NULL)
    *df	=(gsvS[ii]-gsvS[i])/h
      +((3.*bb-1.)*gsvS2[ii]-(3.*aa-1.)*gsvS2[i])*h*ccr6;
  return(0);
}

int ESNeutronWallFlux(double *Pcp)
{
  int i,j,ji,k,ki,kj;
  double dV,Pm,Px,Pt,fm1,fm2,fx1,fx2,ft1,ft2;
  double R1,Z1,T,SDT;
  double r,z,rp,zp,ra,za;

  ji	=ESNp1*ESNa;
  Z1	=ESsz[ji+ESNp/4];
  R1	=ESsr[ji+ESNp/2];

  EZout("sdsd","R1=",R1," Z1=",Z1);

  Pt	=0.;
  Pm	=0.;
  Px	=0.;

  r	=ESaR0[0];
  dV	=ESsp[0]*ESsp[0]*r*ESLc[0];

  rp	=r*r-R1*R1;
  zp	=Z1/sqrt(rp+Z1*Z1);
  fm2	=2.*asin(R1/r)*zp*dV;

  rp	=ESaR0[0]-R1;
  zp	=Z1/sqrt(rp*rp+Z1*Z1);
  fx2	=2.*asin(R1/r)*zp*dV;
  ft2	=EZc2gp*dV;

  for(i=1; i < ESNa; i++){
    ji	=ESNp1*i;

    fm1	=fm2;
    fx1	=fx2;
    ft1	=ft2;

    fm2	=0.;
    fx2	=0.;
    ft2	=0.;
    for(j=0; j < ESNp; j++){
      r		=ESsr[ji];
      z		=ESsz[ji];
      za	=rsT1a[i]+ESsb1a[i]*ESsn1[j];
      zp	=ESsb[i]*EScs1[j];
      ra	=0.;
      rp	=0.;
      ki	=i;
      kj	=0;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	kj	+=j;
	if(kj >=ESNp){
	  kj	-=ESNp;
	}
	ra	+=rcT1a[ki]*EScs1[kj]+rsT1a[ki]*ESsn1[kj];
	rp	+=k*(-rcT[ki]*ESsn1[kj]+rsT[ki]*EScs1[kj]);
      }
      rp	*=2.;
      ra	=2.*ra+rcT1a[i];
      dV	=r*(ra*zp-rp*za)/ESsa[i];

      rp	=r*r-R1*R1;
      zp	=Z1-z;
      zp	/=sqrt(rp+zp*zp);
      za	=Z1+z;
      za	/=sqrt(rp+za*za);
      fm2	+=asin(R1/r)*(zp+za)*dV;

      rp	=r-R1;
      zp	=Z1-z;
      zp	/=sqrt(rp*rp+zp*zp);
      za	=Z1+z;
      za	/=sqrt(rp*rp+za*za);
      fx2	+=asin(R1/r)*(zp+za)*dV;

      ft2	+=EZc2gp*dV;
      ji++;
    }
    dV	=ESsp[i]*ESsp[i]*EZc2gp/ESNp;
    fx2	*=dV;
    fm2	*=dV;
    ft2	*=dV;
    
    rp	=EZcr4*(ESpa[i]-ESpa[i-1]);
    Pm	+=rp*(fm2+fm1);
    Px	+=rp*(fx2+fx1);
    Pt	+=rp*(ft2+ft1);
  }
  rp	=EZcr4*(ESpa[i]-ESpa[i-1]);
  Pm	+=rp*fm2;
  Px	+=rp*fx2;
  Pt	+=rp*ft2;

  T	=15.;
  gsvDT(&SDT,NULL,T);
  SDT	*=1.408e+15/T/T;
  rp	=EZcgm0*3.2044e-2;
  r	=(1.+14.1/3.5)*SDT/(rp*rp);
  Pm	*=r;
  Px	*=r;
  Pt	*=r;

  printf("Pmn=%10.3e Pmx=%10.3e Pn=%10.3e Ptot=%10.3e\n"
	 ,Pm,Px,Px*14.1/(3.5+14.1),Pt);
  *Pcp	=Px;
  return(0);
}

int ESFusionParam()
{
  int i,j,k,ki,kj;
  double h2sh,h12p2sh;
  double SDT,T;
  double PDT,pDT,dP0,dP1,d2P0,d2P1,dp0,dp1,d2p0,d2p1;
  double V,dV0,dV1,d2V0,d2V1;
  double S,r,rp,zp,tE,stE,Pxr;

  V		=ESsa[1]-ESsa[0];
  h2sh		=EZcr2*V;
  h12p2sh	=EZcr12*V*V;

  V		=0.;
  dV1		=0.;
  d2V1		=ESLc[0]+ESVc[0];
  
  pDT	=0.;
  dp1	=0.;
  d2p1	=ESsp[0]*ESsp[0]*d2V1;

  PDT	=0.;
  dP1	=0.;
  d2P1	=ESsp[0]*ESsp[0]*d2V1;

  for(i=1; i < ESNa1; i++){
    dV0		=dV1;
    d2V0	=d2V1;
    d2V1	=ESLc[i]+ESVc[i];
    dV1		=ESsa[i]*d2V1;
    d2V1	+=ESsa[i]*(ESLc1[i]+ESVc1[i]);
    V		+=h2sh*(dV0+dV1)+h12p2sh*(d2V0-d2V1);

    dp0		=dp1;
    d2p0	=d2p1;
    dp1		=ESsp[i]*ESsp[i]*dV1;
    d2p1	=ESsp[i]*(2.*ESsp1a[i]*dV1+ESsp[i]*d2V1);
    pDT		+=h2sh*(dp0+dp1)+h12p2sh*(d2p0-d2p1);

    dP0		=dP1;
    d2P0	=d2P1;
    dP1		=ESsp[i]*ESsp[i]*dV1;
    d2P1	=ESsp[i]*(2.*ESsp1a[i]*dV1+ESsp[i]*d2V1);
    PDT		+=h2sh*(dP0+dP1)+h12p2sh*(d2P0-d2P1);
  }
  dV0	=EZc2gp*EZc2gp*ESaR0[0];
  V	*=dV0;
  pDT	*=3.*dV0;

  T	=15.;
  dP0	=EZcgm0*3.2044e-2;
  Pxr	=PDT*(1.+14.1/3.5)*0.0169*sqrt(0.1*T)*dV0/(dP0*T*dP0*T);

  gsvDT(&SDT,NULL,T);
  SDT	*=1.408e+15/T/T;
  PDT	*=SDT*dV0/(dP0*dP0);
  stE	=0.75*EZcrgm0*ESgbext*ESBt*ESBt*V/PDT;
  tE	=PDT > Pxr ? 0.75*EZcrgm0*ESgbext*ESBt*ESBt*V/(PDT-Pxr) : 0.;
  PDT	*=(1.+14.1/3.5);

#ifdef H
  SDT	=1.4*3.68e+3*exp(-19.94*pow(T,-1./3.))/pow(T,8./3.);
  PDT	*=(1.+14.1/3.5)*SDT*dV0/(dP0*dP0);
#endif

  S	=0.;
  i	=ESNa;
  for(j=0; j < ESNp; j++){
    r	=0.;
    rp	=0.;
    ki	=i;
    kj	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      kj	+=j;
      if(kj >=ESNp){
	kj	-=ESNp;
      }
      r	+=rcT[ki]*EScs1[kj]+rsT[ki]*ESsn1[kj];
      rp+=k*(-rcT[ki]*ESsn1[kj]+rsT[ki]*EScs1[kj]);
    }
    rp	*=2.;
    zp	=ESsb[i]*EScs1[j];
    r	=2.*r+ESaR0[0]+rcT[i];
    S	+=r*sqrt(rp*rp+zp*zp);
  }
  S	*=EZc2gp*EZc2gp/ESNp;
  r	=0.75*EZcrgm0*ESgbext*ESBt*ESBt*(3.5+14.1)*V/(PDT*3.5);
#ifdef XWIN
  CbUserMessage	=Message;
#endif
  ESNeutronWallFlux(&rp);
  r	=rp*14.1/(3.5+14.1);
  if(PDT != 0.){
   rp	/=PDT;
  }
  sprintf(Message
	  ,"bet=%10.3e %%\nVol=%10.3e\nS  =%10.3e\nPDT=%10.3e\n"
	  "Pxr=%10.3e\nPnw=%10.3e\nEsc=%10.3e %%\nt^*=%10.3e\nt_E=%10.3e\n"
	  ,ESgbext*100.,V,S,PDT,Pxr,r,rp*100.,stE,tE);
  printf("bet=%10.3e %%\nVol=%10.3e\nS  =%10.3e\nPDT=%10.3e\n"
	 "Pxr=%10.3e\nPnw=%10.3e\nEsc=%10.3e %%\nt^*=%10.3e\nt_E=%10.3e\n"
	 ,ESgbext*100.,V,S,PDT,Pxr,r,rp*100.,stE,tE);
  ESBootstrap();
  return(0);
}

int ESBalStab()
{
  int i,j,n0,n1;
  static int Fl=0,FlMap=0;
  double Bt,Rext,r0,EZz0,rx,zx;
  double *pr,*pz,*pvr,*pvz,rD[12],zD[12];
  static int iW=0;

  ESGetPlVDAddr(&rdPV,&zdPV);


  FlMap	=0;
#ifndef Tbl_BalSt
  if(nP > NPP){
    nP	=NPP;
  }
  if(nP < ESNp){
    nP	=ESNp;
  }
  nP1	=nP+1;
  nP2	=nP+2;
  nP5	=nP+5;
  nP0	=nP;
  if(Bala0 < 1e-6){
    Bala0	=1e-6;
  }
  if(Bala1 < 1e-6){
    Bala1	=1e-6;
  }
  if(Bala1 > 1.){
    Bala1	=1.;
  }
  if(Balda < 1e-6){
    Balda	=1e-6;
  }
  if(Bala0 > Bala1){
    Bala0=Bala1;
  }
  if((Bala1-Bala0)/Balda > 1000.){
    printf("(a1-a0)/Stp=%10.3e > 1000. - too many steps%c\n"
	   ,(Bala1-Bala0)/Balda,7);
    Balda=0.01*(Bala1-Bala0);
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

#ifdef XWIN
  CbUserPlot	=ESBalStPlot;
#endif
  ESFusionParam();
  switch(Fl){
  case 0:
    lBalSt	=NULL;
    break;
  case 1:
    lBalSt	=ESPstBal;
#ifdef H
    printf("%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n"
	   ,ESBt,ESgbext*100.,ESsq[0],ESsq[ESNa],ESjp[0],ESjb[0]
	   ,ESjb[ESNa],ESjb[ESNa]/ESjb[0]);
#endif
    break;
  case 2:
    lBalSt	=ESPstBal0;
    break;
  case 3:
    lBalSt	=ESPstBal1;
    break;
  case 4:
    lBalSt	=ESPstBal2;
    break;
  default:
    lBalSt	=NULL;
    break;
  }
#ifdef H
  Fl	=0;
#endif
  if(iW == 2 && lBalSt != NULL){
#ifdef XWIN
    SetPlotLimits(5);
#endif
    lBalSt();
  }
  if(iW < 2){
    iW++;
  }
#endif/*Tbl_BalSt*/
  return(0);
}

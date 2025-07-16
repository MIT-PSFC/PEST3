#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "esc.h"

extern int ESInspectMetricTensor(double *x, double *yc, double *ys ,int m
			  ,int i0,int i1,int j0,int j1
			  ,int k0,int k1,int Fl,int Fd,int n);

#ifndef mzl_ESmain
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern void *ESShMem;
extern unsigned int ESmemlev;
extern double ESreltol,ESabstol;
#endif

extern int ESkNormJ;
#ifndef mzl_1Dcore
extern int ESNa,ESNa1;
extern double *ESsa,*ESpa,*ESp3a,*ECb,*ECb1a,*ECb2a;
extern double ESa0,ESav;
extern int ESFlHole,ESFlVirtualB,ESiIv;
#endif

#ifndef mzl_2Dcore
extern int ESNp,ESNp1,ESnAP,ESnAF;
extern int ESFp,ESFp1;
extern int *ESnF,*ESmF,*ESkF,ESNf,ESNf1;
extern int ESMp,ESMp1,ES2Mp,ES2Mp1,ESnMp;
extern double *ESemK,*ESemL,*ESemV,*ESemW,*ESemU1,*ESemU2;
extern double *ESemK2,*ESemL2,*ESemV2,*ESemW2,*ESemU12,*ESemU22;
extern double ESTime; 

extern double *ESg22c,*ESg22s,*ESg22c1,*ESg22s1,*ESg22c2,*ESg22s2;
extern double *ESr22,*ESr222,*ESr22c,*ESr22s,*ESr22c2,*ESr22s2;
extern double *ESg12c,*ESg12s,*ESg11c,*ESg11s;
extern double *ESg12c2,*ESg12s2,*ESg11c2,*ESg11s2;
extern double *ESLc,*ESLs,*ESVc,*ESVs,*ESDPc,*ESDPs;
extern double *ESLc1,*ESLs1,*ESVc1,*ESVs1,*ESDPc1,*ESDPs1;
extern double *ESLc2,*ESLs2,*ESVc2,*ESVs2,*ESDPc2,*ESDPs2;

extern double *ESgh,*ESgt,*EScs1,*ESsn1,*EZcs2,*EZsn2;
extern double R0,Z0,*ECr,*ECz,*EZvr,*EZvz,*EZz2;
extern double *EZdra,*EZdrgt,*EZdza,*EZdzgt;
#endif

#ifndef mzl_BndLr
#define NEQ 34
extern double ESaRx,ESaZx;
extern int ESiAx,*k2kPlV;
extern int ESnBL,ESnBL1,ESNax,ESNax1;
extern double ESLx[],ESVx[],ESDx[],ESgFPlVx;
extern double *EZdx0cs,*EZx0cs,*EZdx0sn,*EZx0sn,*EZx0sep,*EZx1sep,*EZx2sep
,*EZrvbcs,*EZrvbsn,*EZrvb,*EZzvb;
extern double *Dx0,*Dx2,*Ddx2,*Dz2,*Dx2z,*Ddx2z,*Dx0c,*Dx0s,*Dx2c,*Dx2s,
*Ddx2c,*Ddx2s,*Dz2c,*Dz2s,*Dx2zc,*Dx2zs,*Ddx2zc,*Ddx2zs;
extern double *EZrcs,*EZrsn,*EZd1rcs,*EZd1rsn,*EZd2rcs,*EZd2rsn,*EZz0,*EZd1z0,*EZd2z0;

extern double *EZxinf;

extern int nf_X;
static double Z0_vb;

extern double EZdx_0x,d2x_0x,EZdx_0o,d2x_0o;
extern double EZd_vb,EZgd_vb,b_X,z_X0,z_X2,z_X2n,t_X,singtX;

#endif

#ifndef mzl_3Dcore
extern int ESLt,ESNt,ESNt1,ESnAT,ESnAPT,ESLt0,ESNLt,ESNLt1;
extern double *ESaR0,*ESaZ0,*ESsr,*ESsz,*ESsb,*ESsb1a,*ESsb2a;

extern double *gf1T,*csLT,*snLT,*EScs1T,*ESsn1T;
extern double *rcT,*rsT,*rcT1a,*rsT1a,*rsT2a,*rcT2a;
extern double *dR0T,*dZ0T,*dbT,*dbT1a,*dbT2a;
extern double *drcT,*drcT1a,*drcT2a,*drsT,*drsT1a,*drsT2a;

extern double ESgFPlV,ESgFPlV1a;
extern int FlagEFIT;

static int NDT=0,NADT=0,nDT,*k2kT,NDFT=0,nDFT;
static double *rPlVdT,*zPlVdT,*rPncrDT,*zPncrDT;
static double *rPlVcT,*rPlVsT,*drPlVcT,*drPlVsT;
#endif

static int kW=0;
static double bW,rcW[32],rsW[32];

extern double EcRjet[],EcZjet[];
extern int EcNjet,FlagPoints;
extern int ESNsep;
extern double ESRsep[],ESZsep[];

#ifndef mzl_PL
extern double ESBt,ESRBt,ESRext;
extern double ESIpl;
#endif

#ifndef mzl_ESPlV
static int NPlVd=0,nPlVd=12,kPlVd,*k2k,MPlVf=0,NPlVf=0,nPlVf=0,NPlVcs=0,nPlVcs;
static double *rPlVd,*zPlVd,*rPlV,*zPlV;

static double bPlV,Z0PlV,R0PlV,*rPlVc,*rPlVs,*drPlVc,*drPlVs;
static int NrC=0,nrC;
static double *aCr,*rCc,*rCc2a,*rCs,*rCs2a;
static int Neq=0,neq,*indx;
static double *ac,*f,*csPlV,*snPlV;

static int nPncr=0,nPsrf=1,KPncr=1,iPncr,LPncr;
static double *rPncr,*zPncr,*xPncr,*yPncr;

extern int VMFlag;
extern int ECEqIterFl;
#endif

extern int ESEqSolvFl,ESEqSolvInPr,ESEqSolvInCr;
extern int ESEqSolvIt;

extern double *ESqgF,*ESqgF1a,*ESqgF2a
,*ESqgY,*ESqgY1a,*ESqgY2a
,*ESsh,*ESsh1a,*ESsh2a
,*ESqV,*ESqV1a,*ESqV2a
,*ESsp,*ESsp1a,*ESsp2a
,*ESsq,*ESsq1a,*ESsq2a
,*ESgm,*ESgm1a,*ESgm2a
,*ESgY,*ESgY1a,*ESgY2a
;

int ESGetPlVd(double *rPlVd0,double *zPlVd0)
{
  int k;

  rPlVd0[0]	=rPlVd[3];
  zPlVd0[0]	=zPlVd[3];
  rPlVd0[1]	=rPlVd[2];
  zPlVd0[1]	=zPlVd[2];
  rPlVd0[2]	=rPlVd[0];
  zPlVd0[2]	=zPlVd[0];
  rPlVd0[3]	=rPlVd[1];
  zPlVd0[3]	=zPlVd[1];
  for(k=4; k < nPlVd; k++){
    rPlVd0[k]	=rPlVd[k];
    zPlVd0[k]	=zPlVd[k];
  }
  return(0);
}

int ESInit1DPlV()
{
  int i;
  double a,b;

  R0	=1.5;
  Z0	=0.;
  a	=0.5;
  b	=1;
  ESRext=R0;
  ESaRx	=R0;
  ESaZx	=-b;
  ESiAx	=0;
  ESBt	=1.;
  ESRBt	=ESBt*ESRext;
  kPlVd	=4;
  rPlVd[0]= R0;
  rPlVd[1]= R0;
  rPlVd[2]= R0-a;
  rPlVd[3]= R0+a;
  zPlVd[0]= Z0+b;
  zPlVd[1]= Z0-b;
  zPlVd[2]= Z0;
  zPlVd[3]= Z0;
  k2k[0]=1;
  k2k[1]=1;
  k2k[2]=1;
  k2k[3]=1;
  
  for(i = kPlVd; i < nPlVd; i++){
    k2k[i]	=0;
    rPlVd[i]=rPlVd[kPlVd-1];
    zPlVd[i]=zPlVd[kPlVd-1];
  }
  return(0);
}

#ifdef TSC
extern double *ESRLm,*ESZLm;
extern int ESnLm,ESnLm1;

static FILE* lF;
static double p0FU=0.001
,gY0FU=0.249
,rgY0FU
,q0FU=0.6
,q1FU=2.5
,dq0FU=0.78
,dq1FU=5.0
,gYsFU,dqFU,dqFU1;

double ESWonchulP(double gY)
{
  return(p0FU*exp(-gY*rgY0FU));
}

double ESWonchulQ(double gY)
{
  if(gY-gYsFU != 0.){
    return(q0FU+gY*(dqFU+dqFU1*(gY-1.)/(gY-gYsFU)));
  }
  else{
    return(q0FU);
  }
}

int ESFindPattern(char *Str)
{
  char ch,*lc;

  ch	='\0';
  lc	=Str;
  while(*lc != '\0' && feof(lF) == 0){
    ch	=(char)fgetc(lF);
    if(*lc == ch){
      lc++;
    }
    else{
      lc	=Str;
    }
  }
  if(feof(lF)){
    return(0);
  }
  return(1);
}

int Wonchul2ES()
{
  double a,R,b,d;
  static int ic=0;

  lF	=fopen("WPark.dat","r");
  if(lF == NULL){
    return(0);
  }

  if(ESFindPattern("Btor=") == 0 || fscanf(lF,"%lg\n",&ESBt) < 1){
    printf("Wrong or absent data on Btor= in WPark.dat file\n");
    fclose(lF);
    return(0);
  }

  if(ESFindPattern("Rext=") == 0 || fscanf(lF,"%lg\n",&ESRext) < 1){
    printf("absent data on Rext= in WPark.dat file\n");
    fclose(lF);
    return(0);
  }
  if(ESFindPattern("R_pl=") == 0 || fscanf(lF,"%lg\n",&R) < 1){
    printf("Wrong or absent data on R_pl= in WPark.dat file\n");
    fclose(lF);
    return(0);
  }
  if(ESFindPattern("a_pl=") == 0 || fscanf(lF,"%lg\n",&a) < 1){
    printf("Wrong or absent data on a_pl= in WPark.dat file\n");
    fclose(lF);
    return(0);
  }
  if(ESFindPattern("Elon=") == 0 || fscanf(lF,"%lg\n",&b) < 1){
    printf("Wrong or absent data on Elon= in WPark.dat file\n");
    fclose(lF);
    return(0);
  }
  if(ESFindPattern("Tria=") == 0 || fscanf(lF,"%lg\n",&d) < 1){
    printf("Wrong or absent data on Tria= in WPark.dat file\n");
    fclose(lF);
    return(0);
  }
  EcNjet	=ESNp1;
  {
    int i;
    b	*=a;
    for(i=0; i < EcNjet; i++){
      EcRjet[i]	=R+a*cos(ESgt[i]+d*ESsn1[i]); 
      EcZjet[i]	=b*ESsn1[i];
    }
    ESaR0[0]	=R;
    rPlVdT[0]	=R+a*cos(EZcr2*EZcgp+d);
    zPlVdT[0]	=b;
    rPlVdT[1]	=rPlVdT[0];
    zPlVdT[1]	=-b;
    rPlVdT[2]	=R-a;
    zPlVdT[2]	=0.;
    rPlVdT[3]	=R+a;
    zPlVdT[3]	=0.;
  }

  if(ESFindPattern("p0Fu=") == 0 || fscanf(lF,"%lg\n",&p0FU) < 1){
    printf("Wrong or absent data on p0FU= in WPark.dat file\n");
    fclose(lF);
  }
  if(ESFindPattern("gY_0=") == 0 || fscanf(lF,"%lg\n",&gY0FU) < 1){
    printf("Wrong or absent data on gY0FU= in WPark.dat file\n");
    fclose(lF);
  }
  if(ESFindPattern("q0FU=") == 0 || fscanf(lF,"%lg\n",&q0FU) < 1){
    printf("Wrong or absent data on q0FU= in WPark.dat file\n");
    fclose(lF);
  }
  if(ESFindPattern("q1FU=") == 0 || fscanf(lF,"%lg\n",&q1FU) < 1){
    printf("Wrong or absent data on q1FU= in WPark.dat file\n");
    fclose(lF);
  }
  if(ESFindPattern("q'_0=") == 0 || fscanf(lF,"%lg\n",&dq0FU) < 1){
    printf("Wrong or absent data on dq0FU= in WPark.dat file\n");
    fclose(lF);
  }
  if(ESFindPattern("q'_1=") == 0 || fscanf(lF,"%lg\n",&dq1FU) < 1){
    printf("Wrong or absent data on =dq1FU in WPark.dat file\n");
    fclose(lF);
  }
  fclose(lF);
  rgY0FU	=gY0FU != 0. ? 1./gY0FU : 1.;
  dqFU	=q1FU-q0FU;
  dqFU1	=dq1FU-dqFU;
  gYsFU	=dqFU1/(dq0FU-dqFU+dqFU1);
  dqFU1	*=1.-gYsFU;

  return(1);
}

int ES2Wonchul(char *wa){
  static int ic=0;
  FILE* lF;
  double a,R,Z,b,d,q0,q1,dq0,dq1;

  lF	=fopen("WPark.dat",wa);
  if(lF == NULL){
    return(0);
  }
  
  R	=EZcr2*(rPlVdT[3]+rPlVdT[2]);
  Z	=EZcr2*(zPlVdT[0]+zPlVdT[1]);

  a	=EZcr2*(rPlVdT[3]-rPlVdT[2]);
  b	=EZcr2*(zPlVdT[0]-zPlVdT[1])/a;
  d	=(R-rPlVdT[0])/a;
  d	=fabs(d) < 1. ? -asin(d) : d/fabs(d);
  if(ESsq[0] != 0.){
    q0FU	=ESsq[0];
    q1FU	=ESgm[ESNa] != 0. ? 1./ESgm[ESNa] : q1FU;
    if(ESqgY1a[0] != 0.){
      dq0FU	=-q0FU*q0FU*EZcr2*ESgm1a[0]/ESqgY1a[0];
    }
    if(ESqgY1a[ESNa] != 0.){
      dq1FU	=-q1FU*q1FU*EZcr2*ESgm1a[ESNa]/ESqgY1a[ESNa];
    }
  }
  if(*wa == 'w'){
    fprintf(lF,"Keep the names (with the '=' sign) in the same order.\n");
    fprintf(lF,"Otherwise, format is not essential.\n");
    ic	=0;
  }
  fprintf(lF,"------------------ %3d --------------------\n",ic);
  fprintf(lF,"Btor=%10.3e\n",ESBt);
  fprintf(lF,"Rext=%10.3e\n",ESRext);
  fprintf(lF,"R_pl=%10.3e\n",R);
  fprintf(lF,"a_pl=%10.3e\n",a);
  fprintf(lF,"Elon=%10.3e\n",b);
  fprintf(lF,"Tria=%10.3e\n",d);
  fprintf(lF,"p0Fu=%10.3e\n",ESsp[0]);
  fprintf(lF,"gY_0=%10.3e\n",gY0FU);
  fprintf(lF,"q0FU=%10.3e\n",q0FU);
  fprintf(lF,"q1FU=%10.3e\n",q1FU);
  fprintf(lF,"q'_0=%10.3e\n",dq0FU);
  fprintf(lF,"q'_1=%10.3e\n",dq1FU);

  ic++;
  fclose(lF);
  return(1);
}

int ESInit1DLim2PlV(int ESnLm, double *ESRLm, double *ESZLm)
{
  int i;
  double a,b;
  double Rmn,Rmx,Zmn,Zmx;

  if(ESRLm == NULL){
    return(0);
  }
  Rmn	=ESRLm[0];
  Rmx	=ESRLm[0];
  Zmn	=ESZLm[0];
  Zmx	=ESZLm[0];
  for(i=1; i < ESnLm; i++){
    if(Rmn > ESRLm[i]){
      Rmn	=ESRLm[i];
    }
    if(Rmx < ESRLm[i]){
      Rmx	=ESRLm[i];
    }
    if(Zmn > ESZLm[i]){
      Zmn	=ESZLm[i];
    }
    if(Zmx < ESZLm[i]){
      Zmx	=ESZLm[i];
    }
  }
  R0	=EZcr2*(Rmx+Rmn);
  Z0	=EZcr2*(Zmx+Zmn);
  a	=EZcr2*(Rmx-Rmn)*0.9;
  b	=EZcr2*(Zmx-Zmn)*0.9;

  ESaRx	=R0;
  ESaZx	=-b;
  ESiAx	=0;
  ESRext=R0;
  ESBt	=1.;
  ESRBt	=ESBt*ESRext;
  kPlVd	=4;
  rPlVd[0]= R0;
  rPlVd[1]= R0;
  rPlVd[2]= R0-a;
  rPlVd[3]= R0+a;
  zPlVd[0]= Z0+b;
  zPlVd[1]= Z0-b;
  zPlVd[2]= Z0;
  zPlVd[3]= Z0;
  k2k[0]=1;
  k2k[1]=1;
  k2k[2]=1;
  k2k[3]=1;
  for(i = kPlVd; i < nPlVd; i++){
    k2k[i]	=0;
    rPlVd[i]=rPlVd[kPlVd-1];
    zPlVd[i]=zPlVd[kPlVd-1];
  }
  ESRef2RealData(0);
  return(0);
}
#endif

int ESReInit1DPlV()
{
  if(NPlVd < nPlVd){
    int i;
    if(ReInitArray((void**)&k2k,NPlVd,nPlVd,sizeof(int)) < 0){
      FailureAlarm((char*)k2k,"ESReInit1DPlV() - no memory for k2k");
      ESexit(0);
    }
    for(i=NPlVd; i < nPlVd; i++){
      k2k[i]	=0;
    }
    if(ReInitArray((void**)&rPlVd,NPlVd,nPlVd,sizeof(double)) < 0){
      FailureAlarm((char*)rPlVd,"ESReInit1DPlV() - no memory for rPlVd");
      ESexit(0);
    }
    if(ReInitArray((void**)&zPlVd,NPlVd,nPlVd,sizeof(double)) < 0){
      FailureAlarm((char*)zPlVd,"ESReInit1DPlV() - no memory for zPlVd");
      ESexit(0);
    }
    NPlVd =nPlVd;
    for(i=NPlVd; i < nPlVd; i++){
      k2k[i]	=0;
      rPlVd[i]=rPlVd[NPlVd-1];
      zPlVd[i]=zPlVd[NPlVd-1];
    }
    ESmemlev |=0x00000040;
  }
  return(0);
}

int ESDeInit1DPlV()
{
  if(NPlVd){
    free(zPlVd);
    free(rPlVd);
    free(k2k);
    NPlVd=0;
  }
  return(0);
}

int ESReInit2DAPrC()
{
  if(NrC < nrC){
    if(ReInitArray((void**)&aCr,NrC,nrC,sizeof(double)) < 0){
      FailureAlarm((char*)aCr,"ESReInit2DAPrC() - no memory for aCr");
      ESexit(0);
    }
    if(ReInitArray((void**)&rCc,NrC,nrC,sizeof(double)) < 0){
      FailureAlarm((char*)rCc,"ESReInit2DAPrC() - no memory for rCc");
      ESexit(0);
    }
    if(ReInitArray((void**)&rCc2a,NrC,nrC,sizeof(double)) < 0){
      FailureAlarm((char*)rCc2a,"ESReInit2DAPrC() - no memory for rCc2a");
      ESexit(0);
    }
    if(ReInitArray((void**)&rCs,NrC,nrC,sizeof(double)) < 0){
      FailureAlarm((char*)rCs,"ESReInit2DAPrC() - no memory for rCs");
      ESexit(0);
    }
    if(ReInitArray((void**)&rCs2a,NrC,nrC,sizeof(double)) < 0){
      FailureAlarm((char*)rCs2a,"ESReInit2DAPrC() - no memory for rCs2a");
      ESexit(0);
    }
    NrC =nrC; 
    ESmemlev |=0x00200000;
  }
  return(0);
}

int ESDeInit2DAPrC()
{
  if(NrC){
    free(rCs2a);
    free(rCs);
    free(rCc2a);
    free(rCc);
    free(aCr);
    NrC=0;
  }
  return(0);
}

#ifndef mzl_3Dcore
int ESRef2RealData(int iT)
{
  int j,*k2;
  double *pr,*pz;

  iT	%=ESNLt;
  while(iT < ESNt1){
    pr		=rPlVdT+nPlVd*iT;
    pz		=zPlVdT+nPlVd*iT;
    k2		=k2kT+nPlVd*iT;
    ESaR0[iT]	=R0;
    ESaZ0[iT]	=Z0;
    for(j=0; j < nPlVd; j++){
      pr[j]	=rPlVd[j];
      pz[j]	=zPlVd[j];
      k2[j]	=k2k[j];
    }
    iT	+=ESNLt;
  }
  return(0);
}

int ESReal2RefData(int iT)
{
  int j,*k2;
  double *pr,*pz;

  pr	=rPlVdT+nPlVd*iT;
  pz	=zPlVdT+nPlVd*iT;
  k2	=k2kT+nPlVd*iT;
  R0	=ESaR0[iT];
  Z0	=ESaZ0[iT];
  for(j=0; j < nPlVd; j++){
    rPlVd[j]	=pr[j];
    zPlVd[j]	=pz[j];
    k2k[j]	=k2[j] ? 1 : 0;
  }
  return(0);
}

int ESReal2RefDataX(int iT)
{
  int j,*k2;
  double *pr,*pz;

  pr	=rPlVdT+nPlVd*iT;
  pz	=zPlVdT+nPlVd*iT;
  k2	=k2kT+nPlVd*iT;
  R0	=ESaR0[iT];
  Z0	=ESaZ0[iT];
  for(j=0; j < nPlVd; j++){
    rPlVd[j]	=pr[j];
    zPlVd[j]	=pz[j];
    k2k[j]	=k2[j] ? 1 : 0;
    k2kPlV[j]	=j;
  }

  if(ESaZx > 0.){
    ESaRx	=rPlVd[0];
    ESaZx	=zPlVd[0];
  }
  else{
    k2kPlV[0]	=1;
    k2kPlV[1]	=0;
    ESaRx	=rPlVd[1];
    ESaZx	=zPlVd[1];
  }
  return(0);
}

#ifdef ESC
int ESPncr2Sort(int n,int iPsrf)
{
  int i,ibot,il,ir,ii,ii1,j,j1,nd;
  double u[2],v[2],pr[8],pz[8],x[100],*px,*py;
  double s,ss,sss,dx,dy,dx1,dy1,dx2,dy2,z1,EZz2;
  
  ESGetPncrAddr(&nPncr,&nPsrf,&rPncr,&zPncr,&xPncr,&yPncr);
  if(nPncr == 0) return(0);
  KPncr	=nPncr/nPsrf;
  LPncr	=ESNt/ESLt0;
  i	=n == ESNt ? 0 : n;
  ii	=i+KPncr*iPsrf;
  pr[0]	=rPncr[ii];
  pz[0]	=zPncr[ii];
  pr[1]	=rPncr[ii];
  pz[1]	=zPncr[ii];
  nd	=0;
  while(i < KPncr){
    xPncr[nd]	=rPncr[ii];
    yPncr[nd]	=zPncr[ii];
    if(pz[0] < yPncr[nd]){
      pr[0]	=xPncr[nd];
      pz[0]	=yPncr[nd];
      z1	=xPncr[0];
      xPncr[0]	=xPncr[nd];
      xPncr[nd]	=z1;
      z1	=yPncr[0];
      yPncr[0]	=yPncr[nd];
      yPncr[nd]	=z1;
    }
    if(pz[1] > yPncr[nd]){
      pr[1]	=xPncr[nd];
      pz[1]	=yPncr[nd];
    }
    nd++;
    i +=LPncr;
    ii +=LPncr;
  }
  xPncr[nd]	=xPncr[0];
  yPncr[nd]	=yPncr[0];
  j	=0;
  dy	=pz[0]-pz[1];
  s	=dy*dy;
  ss	=s;
  ii	=0;
  ii1	=0;
  i	=1;
  while(i < nd){
    dx2	=xPncr[i]-xPncr[j];
    dy2	=yPncr[i]-yPncr[j];
    sss	=dx2*dx2+dy2*dy2;
    if(xPncr[i] < xPncr[0]){
      if(s > sss){
	s	=sss;
	ii	=i;
      }
    }
    else{
      if(ss > sss){
	ss	=sss;
	ii1	=i;
      }
    }
    i++;
  }
  j	=1;
  z1		=xPncr[j];
  xPncr[j]	=xPncr[ii];
  xPncr[ii]	=z1;
  z1		=yPncr[j];
  yPncr[j]	=yPncr[ii];
  yPncr[ii]	=z1;
  j1	=nd-1;
  z1		=xPncr[j1];
  xPncr[j1]	=xPncr[ii1];
  xPncr[ii1]	=z1;
  z1		=yPncr[j1];
  yPncr[j1]	=yPncr[ii1];
  yPncr[ii1]	=z1;
  while(j1-j > 1){
    i	=j+1;
    dx2	=xPncr[i]-xPncr[j];
    dy2	=yPncr[i]-yPncr[j];
    s	=dx2*dx2+dy2*dy2;
    dx2	=xPncr[i]-xPncr[j1];
    dy2	=yPncr[i]-yPncr[j1];
    ss	=dx2*dx2+dy2*dy2;
    ii	=i;
    ii1	=i;
    i++;
    while(i < j1){
      dx2	=xPncr[i]-xPncr[j];
      dy2	=yPncr[i]-yPncr[j];
      sss	=dx2*dx2+dy2*dy2;
      if(s > sss){
	s	=sss;
	ii	=i;
      }
      dx2	=xPncr[i]-xPncr[j1];
      dy2	=yPncr[i]-yPncr[j1];
      sss	=dx2*dx2+dy2*dy2;
      if(ss > sss){
	ss	=sss;
	ii1	=i;
      }
      i++;
    }
    if((j1 -j) > 2){
      dx	=xPncr[j]-xPncr[j-1];
      dy	=yPncr[j]-yPncr[j-1];
      dx1	=xPncr[j]-xPncr[ii];
      dy1	=yPncr[j]-yPncr[ii];
      s		=-(dx1*dx+dy1*dy)/(dx*dx+dy*dy);
      dx1	+=s*dx;
      dy1	+=s*dy;
      s		=dx1*dx1+dy1*dy1;
      dx	=xPncr[j1]-xPncr[j1+1];
      dy	=yPncr[j1]-yPncr[j1+1];
      dx1	=xPncr[j1]-xPncr[ii1];
      dy1	=yPncr[j1]-yPncr[ii1];
      ss	=-(dx1*dx+dy1*dy)/(dx*dx+dy*dy);
      dx1	+=ss*dx;
      dy1	+=ss*dy;
      ss	=dx1*dx1+dy1*dy1;
      if(s < ss){
	j++;
	z1	=xPncr[j];
	xPncr[j]	=xPncr[ii];
	xPncr[ii]	=z1;
	z1	=yPncr[j];
	yPncr[j]	=yPncr[ii];
	yPncr[ii]	=z1;
      }
      else{
	j1--;
	z1	=xPncr[j1];
	xPncr[j1]	=xPncr[ii1];
	xPncr[ii1]=z1;
	z1	=yPncr[j1];
	yPncr[j1]	=yPncr[ii1];
	yPncr[ii1]=z1;
      }
    }
    else{
      j++;
      z1	=xPncr[j];
      xPncr[j]	=xPncr[ii];
      xPncr[ii]	=z1;
      z1	=yPncr[j];
      yPncr[j]	=yPncr[ii];
      yPncr[ii]	=z1;
    }
  }
  ii	=nd+1;
  j	=0;
  x[0]	=0.;
  pz[1]	=yPncr[0];
  ibot	=0;
  for(i=1; i < ii; i++){
    dx	=xPncr[i]-xPncr[j];
    dy	=yPncr[i]-yPncr[j];
    x[i]=x[j]+sqrt(dx*dx+dy*dy);
    if(pz[1] > yPncr[i]){
      pz[1]=yPncr[i];
      ibot=i;
    }
    j++;
  }
  z1	=EZcr2*(pz[0]+pz[1]);
  pz[2]	=pz[0]-z1;
  pz[3]	=pz[2];
  il	=0;
  ir	=0;
  s	=1./x[nd];
  for(i=1; i < ii; i++){
    x[i]	*=s;
    ss	=fabs(yPncr[i]-z1);
    if(i < ibot && pz[2] > ss){
      il	=i;
      pz[2]	=ss;
    }
    if(i > ibot && pz[3] > ss){
      ir	=i;
      pz[3]	=ss;
    }
  }
  x[nd]	=1.;
  EZspl_p(xPncr,xPncr+ii,x,nd);
  EZspl_p(yPncr,yPncr+ii,x,nd);
  i	=0;
  j	=1;
  dx	=x[j]-x[i];
  dx1	=EZcr6*dx;
  dy	=(yPncr[j]-yPncr[i])/dx;
  dy1	=dy+(2.*yPncr[j+ii]+yPncr[i+ii])*dx1;
  while(dy1 > 0.){
    i++;
    j++;
    dx	=x[j]-x[i];
    dx1	=EZcr6*dx;
    dy	=(yPncr[j]-yPncr[i])/dx;
    dy1	=dy+(2.*yPncr[j+ii]+yPncr[i+ii])*dx1;
  }
  dy1	=dy-(yPncr[j+ii]+2.*yPncr[i+ii])*dx1;
  if(dy1 < 0.){
    i	=nd;
    j	=i+1;
  }
  while(dy1 < 0.){
    i--;
    j--;
    dx	=x[j]-x[i];
    dx1	=EZcr6*dx;
    dy	=(yPncr[j]-yPncr[i])/dx;
    dy1	=dy-(yPncr[j+ii]+2.*yPncr[i+ii])*dx1;
  }
  ss	=yPncr[j+ii];
  s	=yPncr[i+ii];
  {
    double a,b,c;
    c	=dy+EZcr4*(s-ss)*dx1;
    b	=-1.5*(ss+s)*dx1;
    a	=3.*(ss-s)*dx1;
    s	=sqrt(b*b-a*c);
    dx1	=b > 0. ? c/(b+s) : c/(b-s);
  }
  dx	=x[i]+dx*(EZcr2+dx1);
  SetIspl(x,dx,nd);
  splR(pz,&s,yPncr,yPncr+ii);
  splR(pr,&s,xPncr,xPncr+ii);
  j	=ibot+1;
  dx	=x[j]-x[ibot];
  dx1	=EZcr6*dx;
  dy	=(yPncr[j]-yPncr[ibot])/dx;
  dy1	=dy+(2.*yPncr[j+ii]+yPncr[ibot+ii])*dx1;
  while(dy1 < 0.){
    ibot++;
    j++;
    dx	=x[j]-x[ibot];
    dx1	=EZcr6*dx;
    dy	=(yPncr[j]-yPncr[ibot])/dx;
    dy1	=dy+(2.*yPncr[j+ii]+yPncr[ibot+ii])*dx1;
  }
  dy1	=dy-(yPncr[j+ii]+2.*yPncr[ibot+ii])*dx1;
  while(dy1 > 0.){
    ibot--;
    j--;
    dx	=x[j]-x[ibot];
    dx1	=EZcr6*dx;
    dy	=(yPncr[j]-yPncr[ibot])/dx;
    dy1	=dy-(yPncr[j+ii]+2.*yPncr[ibot+ii])*dx1;
  }
  ss	=yPncr[j+ii];
  s	=yPncr[ibot+ii];
  {
    double a,b,c;
    c	=dy+EZcr4*(s-ss)*dx1;
    b	=-1.5*(ss+s)*dx1;
    a	=3.*(ss-s)*dx1;
    s	=sqrt(b*b-a*c);
    dx1	=b > 0. ? c/(b+s) : c/(b-s);
  }
  dx	=x[ibot]+dx*(EZcr2+dx1);
  SetIspl(x,dx,nd);
  splR(pz+1,&s,yPncr,yPncr+ii);
  splR(pr+1,&s,xPncr,xPncr+ii);
  for(i=2; i < 8; i++){
    switch(i){
    case 2:
      z1	=EZcr2*(pz[0]+pz[1]);
      u[0]	=x[il];
      u[1]	=x[il+1];
      break;
    case 3:
      z1	=EZcr2*(pz[0]+pz[1]);
      u[0]	=x[ir];
      u[1]	=x[ir+1];
      break;
    case 4:
      z1	=EZcr2*(pz[0]+pz[1])+EZcr4*(pz[0]-pz[1]);
      u[0]	=x[il/2];
      u[1]	=x[il/2+1];
      break;
    case 5:
      z1	=EZcr2*(pz[0]+pz[1])-EZcr4*(pz[0]-pz[1]);
      u[0]	=x[(il+ibot)/2];
      u[1]	=x[(il+ibot)/2+1];
      break;
    case 6:
      z1	=EZcr2*(pz[0]+pz[1])-EZcr4*(pz[0]-pz[1]);
      u[0]	=x[(ir+ibot)/2];
      u[1]	=x[(ir+ibot)/2+1];
      break;
    case 7:
      z1	=EZcr2*(pz[0]+pz[1])+EZcr4*(pz[0]-pz[1]);
      u[0]	=x[(ir+nd)/2];
      u[1]	=x[(ir+nd)/2+1];
      break;
    }
    j	=0;
    do{
      EZzero(&j,u,v,&dx,s);
      SetIspl(x,dx,nd);
      splR(pz+i,&ss,yPncr,yPncr+ii);
      s	=pz[i]-z1;
      ss= j == 1 ? 1. : fabs(dx-u[0]);
    }while(j > 0 && j < 15 && ss > 1e-4);
    splR(pr+i,&s,xPncr,xPncr+ii);
  }

  dx	=EZcr4/nd;
  px	=xPncr+2*ii;
  py	=yPncr+2*ii;
  j	=4*nd+1;
  for(i=0; i < j; i++){
    s	=dx*i;
    if(s > 1.){
      s	-=1.;
    }
    SetIspl(x,s,nd);
    splR(px+i,&ss,xPncr,xPncr+ii);
    splR(py+i,&ss,yPncr,yPncr+ii);
  }
  Plot2D(3,px,py,j,6,0,0,0);
  Plot2D(3,pr,pz,8,6,1,14,0);
  return(nd);
}
#endif

int ESPncr2PlVData(int n,int iPsrf)
{
#ifdef H
  int i,ibot,il,ir,ii,ii1,j,j1,nd,*k2;
  double *pr,*pz,*x;
  double u[2],v[2];
  double b,s,ss,sss,dx,dy,dx1,dy1,dx2,dy2,z1,EZz2;
  
  k2	=k2kT+nPlVd*n;
  pr	=rPlVdT+nPlVd*n;
  pz	=zPlVdT+nPlVd*n;

  while(dy1 > 0.){
    ibot--;
    j--;
    dx	=x[j]-x[ibot];
    dx1	=EZcr6*dx;
    dy	=(yPncr[j]-yPncr[ibot])/dx;
    dy1	=dy-(yPncr[j+ii]+2.*yPncr[ibot+ii])*dx1;
  }
  ss	=yPncr[j+ii];
  s	=yPncr[ibot+ii];
  {
    double a,b,c;
    c	=dy+EZcr4*(s-ss)*dx1;
    b	=-1.5*(ss+s)*dx1;
    a	=3.*(ss-s)*dx1;
    s	=sqrt(b*b-a*c);
    dx1	=b > 0. ? c/(b+s) : c/(b-s);
  }
  dx	=x[ibot]+dx*(EZcr2+dx1);
  SetIspl(x,dx,nd);
  splR(pz+1,&s,yPncr,yPncr+ii);
  splR(pr+1,&s,xPncr,xPncr+ii);

  Z0	=EZcr2*(pz[0]+pz[1]);
  b	=EZcr2*(pz[0]-pz[1]);
  for(i=2; i < 12; i++){
    switch(i){
    case 2:
      z1	=Z0;
      j		=il;
      break;
    case 3:
      z1	=Z0;
      j		=ir;
      break;
    case 4:
      z1	=Z0+EZcr2*b;
      j		=il/2;
      break;
    case 5:
      z1	=Z0-EZcr2*b;
      j		=(il+ibot)/2;
      break;
    case 6:
      z1	=Z0-EZcr2*b;
      j		=(ir+ibot)/2;
      break;
    case 7:
      z1	=Z0+EZcr2*b;
      j		=(ir+nd)/2;
      break;
    case 8:
      z1	=Z0+0.86*b;
      z1	=Z0+0.8*b;
      j		=il/4;
      break;
    case 9:
      z1	=Z0-0.86*b;
      z1	=Z0-0.8*b;
      j		=(il+ibot)/2+(ibot-il)/4;
      break;
    case 10:
      z1	=Z0-0.86*b;
      z1	=Z0-0.8*b;
      j		=ibot+(ir-ibot)/4;
      break;
    case 11:
      z1	=Z0+0.86*b;
      z1	=Z0+0.8*b;
      j		=(ir+nd)/2+(nd-ir)/4;
      break;
    }
    u[0]	=x[j];
    u[1]	=x[j+1];
    j	=0;
    do{
      EZzero(&j,u,v,&dx,s);
      SetIspl(x,dx,nd);
      splR(pz+i,&ss,yPncr,yPncr+ii);
      s		=pz[i]-z1;
      ss	= j == 1 ? 1. : fabs(dx-u[0]);
    }while(j > 0 && j < 15 && ss > 1e-4);
    splR(pr+i,&s,xPncr,xPncr+ii);
  }
  ESaR0[n]	=EZcr2*(pr[2]+pr[3]);
  ESaZ0[n]	=Z0;
  for(j=0; j < 12; j++){
    k2[j]	=1;
  }
  for(j=12; j < nPlVd; j++){
    k2[j]	=0;
    pr[j]	=pr[3];
    pz[j]	=pz[0];
  }
  return(nd);
#endif
  return(0);
}

#ifdef ESC
int ESPncr2PlVData0(int n,int iPsrf)
{
  int i,ii,j,ibot,nd,*k2;
  double *pr,*pz;
  double b,s,ss,dx,dy,dx1,dy1,z1,EZz2;
 
  k2	=k2kT+nPlVd*n;
  pr	=rPlVdT+nPlVd*n;
  pz	=zPlVdT+nPlVd*n;
  nd	=0;
  i	=n == ESNt ? 0 : n;
  ii	=i+KPncr*iPsrf;
  while(i < KPncr){
    xPncr[nd]	=rPncr[ii];
    yPncr[nd]	=zPncr[ii];
    nd++;
    i +=LPncr;
    ii +=LPncr;
  }
  
  pr[0]	=xPncr[0];
  pz[0]	=yPncr[0];
  for(i=0; i < nd; i++){
    if(pz[0] < yPncr[i]){
      pr[0]	=xPncr[i];
      pz[0]	=yPncr[i];
      z1	=xPncr[0];
      xPncr[0]	=xPncr[i];
      xPncr[i]	=z1;
      z1	=yPncr[0];
      yPncr[0]	=yPncr[i];
      yPncr[i]	=z1;
    }
  }
  j	=0;
  i	=j+1;
  while(i < nd){
    if(xPncr[i] < xPncr[0]){
      break;
    }
    i++;
  }
  dx	=xPncr[i]-xPncr[0];
  dy	=yPncr[i]-yPncr[0];
  s	=dx*dx+dy*dy;

  ii	=i;
  while(i < nd){
    if(xPncr[i] < xPncr[0]){
      dx1	=xPncr[i]-xPncr[0];
      dy1	=yPncr[i]-yPncr[0];
      ss	=dx1*dx1+dy1*dy1;
      if(s > ss){
	s	=ss;
	ii	=i;
      }
    }
    i++;
  }
  j++;
  pr[1]	=pr[0];
  pz[1]	=pz[0];
  ibot	=0;
  while(j < nd){
    z1	=xPncr[j];
    xPncr[j]	=xPncr[ii];
    xPncr[ii]	=z1;
    z1	=yPncr[j];
    yPncr[j]	=yPncr[ii];
    yPncr[ii]	=z1;
    dx	=xPncr[j]-xPncr[j-1];
    dy	=yPncr[j]-yPncr[j-1];
    i	=j+1;
    dx1	=xPncr[i]-xPncr[j];
    dy1	=yPncr[i]-yPncr[j];
    s	=dx1*dx1+dy1*dy1;
    ii	=i;
    while(i < nd){
      dx1	=xPncr[i]-xPncr[j];
      dy1	=yPncr[i]-yPncr[j];
      ss	=dx1*dx1+dy1*dy1;
      if(dx*dx1+dy*dy1 > 0. && s > ss){
	s	=ss;
	ii	=i;
      }
      i++;
    }
    if(pz[1] > yPncr[j]){
      pr[1]	=xPncr[j];
      pz[1]	=yPncr[j];
      ibot	=j;
    }
    j++;
  }
  pr[2]	=xPncr[1]-pr[0];
  pz[2]	=(yPncr[1]-pz[0])/pr[2];
 
  i	=nd-1;
  pr[3]	=xPncr[i]-pr[0];
  pz[3]	=(yPncr[i]-pz[0])/pr[3];
  z1	=pz[3]*pr[2]-pz[2]*pr[3];
  EZz2	=pz[2]-pz[3];
  pr[0] -=EZcr2*z1/EZz2;
  pz[0]	-=EZcr4*z1*z1/(EZz2*(pr[2]-pr[3]));

  i	=ibot+1;
  pr[2]	=xPncr[i]-pr[1];
  pz[2]	=(yPncr[i]-pz[1])/pr[2];
  i	=ibot-1;
  pr[3]	=xPncr[i]-pr[1];
  pz[3]	=(yPncr[i]-pz[1])/pr[3];
  z1	=pz[3]*pr[2]-pz[2]*pr[3];
  EZz2	=pz[2]-pz[3];
  pr[1] -=EZcr2*z1/EZz2;
  pz[1]	-=EZcr4*z1*z1/(EZz2*(pr[2]-pr[3]));

  b	=EZcr2*(pz[0]-pz[1]);
  Z0	=EZcr2*(pz[0]+pz[1]);

  ii	=0;
  s	=fabs(yPncr[ii]);
  for(i=ii+1; i < ibot; i++){
    ss	=fabs(yPncr[i]);
    if(s > ss){
      s	=ss;
      ii=i;
    }
  }
  pr[2]	=xPncr[ii];
  pz[2]	=yPncr[ii];
  
  j	=ii;
  z1	=Z0+b*EZcr2;
  ii	=0;
  s	=fabs(yPncr[ii]-z1);
  for(i=ii+1; i < j; i++){
    ss	=fabs(yPncr[i]-z1);
    if(s > ss){
      s	=ss;
      ii=i;
    }
  }
  pr[4]	=xPncr[ii];
  pz[4]	=yPncr[ii];

  z1	=Z0-b*EZcr2;
  ii	=j;
  s	=fabs(yPncr[ii]-z1);
  for(i=ii+1; i < ibot; i++){
    ss	=fabs(yPncr[i]-z1);
    if(s > ss){
      s	=ss;
      ii=i;
    }
  }
  pr[5]	=xPncr[ii];
  pz[5]	=yPncr[ii];

  ii	=ibot;
  s	=fabs(yPncr[ii]);
  for(i=ii+1; i < nd; i++){
    ss	=fabs(yPncr[i]);
    if(s > ss){
      s	=ss;
      ii=i;
    }
  }
  pr[3]	=xPncr[ii];
  pz[3]	=yPncr[ii];

  j	=ii;
  z1	=Z0-b*EZcr2;
  ii	=ibot;
  s	=fabs(yPncr[ii]-z1);
  for(i=ii+1; i < j; i++){
    ss	=fabs(yPncr[i]-z1);
    if(s > ss){
      s	=ss;
      ii=i;
    }
  }
  pr[6]	=xPncr[ii];
  pz[6]	=yPncr[ii];

  z1	=Z0+b*EZcr2;
  ii	=j;
  s	=fabs(yPncr[ii]-z1);
  for(i=ii+1; i < nd; i++){
    ss	=fabs(yPncr[i]-z1);
    if(s > ss){
      s	=ss;
      ii=i;
    }
  }
  pr[7]	=xPncr[ii];
  pz[7]	=yPncr[ii];
  ESaR0[n]	=EZcr2*(pr[2]+pr[3]);
  ESaZ0[n]	=Z0;
  for(j=0; j < 8; j++){
    k2[j]	=1;
  }
  for(j=8; j < nPlVd; j++){
    k2[j]	=0;
    pr[j]	=pr[3];
    pz[j]	=pz[0];
  }
  return(0);
}
#endif

int ESInit2DTPlV()
{
  int n;
  
  if(nPncr){
    for(n=0; n < ESNt1; n++){
      ESPncr2PlVData(n,nPsrf-1);
    }
  }
  else{
    for(n=0; n < ESNt1; n++){
      ESRef2RealData(n);
    } 
  }
  return(0);
}

int ESReInit2DTPlV()
{
  nDT	=ESNt1*nPlVd;
  if(NDT < nDT){
    int i;
    if(ReInitArray((void**)&k2kT,NDT,nDT,sizeof(int)) < 0){
      FailureAlarm((char*)k2kT,"ESReInit2DTPlV() - no memory for k2kT");
      ESexit(0);
    }
    for(i=NDT; i < nDT; i++){
      k2kT[i]	=0;
    }
    if(ReInitArray((void**)&rPlVdT,NDT,nDT,sizeof(double)) < 0){
      FailureAlarm((char*)rPlVdT,"ESReInit2DTPlV() - no memory for rPlVdT");
      ESexit(0);
    }
    if(ReInitArray((void**)&zPlVdT,NDT,nDT,sizeof(double)) < 0){
      FailureAlarm((char*)zPlVdT,"ESReInit2DTPlV() - no memory for zPlVdT");
      ESexit(0);
    }
    i	=nDT*ESNa1;
    if(ReInitArray((void**)&rPncrDT,NADT,i,sizeof(double)) < 0){
      FailureAlarm((char*)rPncrDT,"ESReInit2DTPlV() - no memory for rPncrDT");
      ESexit(0);
    }
    if(ReInitArray((void**)&zPncrDT,NADT,i,sizeof(double)) < 0){
      FailureAlarm((char*)zPncrDT,"ESReInit2DTPlV() - no memory for zPncrDT");
      ESexit(0);
    }
    NADT =i;
    NDT =nDT;
    ESmemlev |=0x00000080;
  }
  return(0);
}

int ESDeInit2DTPlV()
{
  if(NDT){
    free(zPncrDT);
    free(rPncrDT);
    free(zPlVdT);
    free(rPlVdT);
    free(k2kT);
    NDT=0;
    NADT=0;
  }
  return(0);
}
#endif

#ifndef mzl_2Dcore
int ESReInit1DFPlV()
{
  if(NPlVf < MPlVf){
    int m;
    if(ReInitArray((void**)&rPlVc,NPlVf,MPlVf,sizeof(double)) < 0){
      FailureAlarm((char*)rPlVc,"ESReInit1DFPlV() - no memory for rPlVc");
      ESexit(0);
    }
    if(ReInitArray((void**)&rPlVs,NPlVf,MPlVf,sizeof(double)) < 0){
      FailureAlarm((char*)rPlVs,"ESReInit1DFPlV() - no memory for rPlVs");
      ESexit(0);
    }
    if(ReInitArray((void**)&drPlVc,NPlVf,MPlVf,sizeof(double)) < 0){
      FailureAlarm((char*)drPlVc,"ESReInit1DFPlV() - no memory for drPlVc");
      ESexit(0);
    }
    if(ReInitArray((void**)&drPlVs,NPlVf,MPlVf,sizeof(double)) < 0){
      FailureAlarm((char*)drPlVs,"ESReInit1DFPlV() - no memory for drPlVs");
      ESexit(0);
    }
    for(m=NPlVf; m < MPlVf; m++){
      rPlVc[m]	=0.;
      rPlVs[m]	=0.;
      drPlVc[m]	=0.;
      drPlVs[m]	=0.;
    }
    NPlVf=MPlVf;
    ESmemlev |=0x00000100;
  }
  return(0);
}

int ESDeInit1DFPlV()
{
  if(NPlVf){
    free(drPlVs);
    free(drPlVc);
    free(rPlVs);
    free(rPlVc);
    NPlVf=0;
  }
  return(0);
}
#endif

#ifndef mzl_3Dcore
int ESReInit2DFTPlV()
{
  if(NDFT < nDFT){
    if(ReInitArray((void**)&rPlVcT,NDFT,nDFT,sizeof(double)) < 0){
      FailureAlarm((char*)rPlVcT,"ESReInit2DFTPlV() - no memory for rPlVcT");
      ESexit(0);
    }
    if(ReInitArray((void**)&rPlVsT,NDFT,nDFT,sizeof(double)) < 0){
      FailureAlarm((char*)rPlVsT,"ESReInit2DFTPlV() - no memory for rPlVsT");
      ESexit(0);
    }
    if(ReInitArray((void**)&drPlVcT,NDFT,nDFT,sizeof(double)) < 0){
      FailureAlarm((char*)drPlVcT,"ESReInit2DFTPlV() - no memory for drPlVcT");
      ESexit(0);
    }
    if(ReInitArray((void**)&drPlVsT,NDFT,nDFT,sizeof(double)) < 0){
      FailureAlarm((char*)drPlVsT,"ESReInit2DFTPlV() - no memory for drPlVsT");
      ESexit(0);
    }
    NDFT =nDFT;
    ESmemlev |=0x00000200;
  }
  return(0);
}

int ESDeInit2DFTPlV()
{
  if(NDFT){
    free(drPlVsT);
    free(drPlVcT);
    free(rPlVcT);
    free(rPlVsT);
    NDFT=0;
  }
  return(0);
}
#endif

#ifndef mzl_2Dcore
int ESInit2DFPPlV()
{
  int j,m,jm;
  double s,t;
  jm	=0;
  for(j =0; j < ESNp1; j++){
    s	=ESgt[j];
    for(m=1; m < MPlVf; m++){
      t	=m*s;
      csPlV[jm]	=2.*cos(t);
      snPlV[jm]	=2.*sin(t);
      jm++;
    }
  }

  nPlVf=MPlVf;
  return(0);
}

int ESReInit2DFPPlV()
{
  if(NPlVcs < nPlVcs){
    if(ReInitArray((void**)&csPlV,NPlVcs,nPlVcs,sizeof(double)) < 0){
      FailureAlarm((char*)csPlV,"ESReInit2DFPPlV() - no memory for csPlV");
      ESexit(0);
    }
    if(ReInitArray((void**)&snPlV,NPlVcs,nPlVcs,sizeof(double)) < 0){
      FailureAlarm((char*)snPlV,"ESReInit2DFPPlV() - no memory for snPlV");
      ESexit(0);
    }
    NPlVcs=nPlVcs;
    ESmemlev |=0x00000400;
  }
  return(0);
}

int ESDeInit2DFPPlV()
{
  if(NPlVcs){
    free(csPlV);
    free(snPlV);
    NPlVcs=0;
  }
  return(0);
}

int ESReInitPlVeq(int n, double **a, double **ff,int **i)
{
  if(Neq < n){
    if(ReInitArray((void**)&indx,Neq,n,sizeof(int)) < 0){
      FailureAlarm((char*)indx,"ESReIniteq() - no memory for indx");
      ESexit(0);
    }
    if(ReInitArray((void**)&ac,Neq*Neq,n*n,sizeof(double)) < 0){
      FailureAlarm((char*)ac,"ESReIniteq() - no memory for ac");
      ESexit(0);
    }
    if(ReInitArray((void**)&f,Neq,n,sizeof(double)) < 0){
      FailureAlarm((char*)f,"ESReIniteq() - no memory for f");
      ESexit(0);
    }
    Neq	=n;
    ESmemlev |=0x00000800;
  }
  *i	=indx;
  *a	=ac;
  *ff	=f;
  return(0);
}

int ESDeInitPlVeq()
{
  if(Neq){
    free(f);
    free(ac);
    free(indx);
    Neq=0;
  }
  return(0);
}
#endif

#ifndef mzl_3Dcore
int ESRealGeom(int iT)
{
  int i,i0,j,ji,in,k,ki,kj;
  double *pr,*pz,b;
  double *rc,*rs,*prc,*prs;

  pr	=ESsr+ESnAP*iT;
  pz	=ESsz+ESnAP*iT;
  rc	=EZdra;
  rs	=rc+ESFp1;
  in	=ESNa1*iT;
  i	=ESnAF*iT;
  prc	=rcT+i;
  prs	=rsT+i;

  ji	=0;
  i0	=0;
  if(ESa0 == 0.){
    for(j=0; j < ESNp1; j++){
      pr[ji]	=ESaR0[iT];
      pz[ji]	=ESaZ0[iT];
      ji++;
    }
    in++;
    i0++;
  }
  for(i=i0; i < ESNa1; i++){
    b	=ESsb[in];
    k	=0;
    ki	=i;
    rc[0]=ESaR0[iT]+prc[ki];
    rs[0]=ESaZ0[iT]+prs[ki];
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      rc[k]	=2.*prc[ki];
      rs[k]	=2.*prs[ki];
    }
    for(j=0; j < ESNp; j++){
      pr[ji]	=rc[0];
      kj	=0;
      for(k=1; k < ESFp1; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	pr[ji]	+=rc[k]*EScs1[kj]+rs[k]*ESsn1[kj];
      }
      pz[ji]	=rs[0]+b*ESsn1[j];
      ji++;
    }
    kj		=ji-ESNp;
    pr[ji]	=pr[kj];
    pz[ji]	=pz[kj];
    ji++;
    in++;
  }
  k	=iT+ESNLt;
  while(k < ESNt1){
    ji	=ESnAP*iT;
    ki	=ESnAP*k;
    for(i=0; i < ESNa1; i++){
      for(j=0; j < ESNp1; j++){
	ESsr[ki]	=ESsr[ji];
	ESsz[ki]	=ESsz[ji];
	ji++;
	ki++;
      }
    }
    k	+=ESNLt;
  }
  return(0);
}

int ESRef2RealGeom(int iT)
{
  int i,in,in1,n,m,k,kk,ki;
#ifdef XWIN
  double dgc,dgs,db,A[ESNa1];
#else
  double dgc,dgs,db,A[129];
#endif

  iT	%=ESNLt;
  drPlVc[0]	+=ESaR0[iT]-R0;
  drPlVs[0]	+=ESaZ0[iT]-Z0;

  ESaR0[iT]	=R0;
  ESaZ0[iT]	=Z0;
  in		=ESNa1*iT;
  db		=bPlV-ESsb[in+ESNa];
  for(i=0; i < ESNa1; i++){

    ESsb[in]	+=ESsa[i]*db;
    ESsb1a[in]	+=db;
#ifdef H
    ESsb[in]	+=ESsa[i]*db*(0.6+0.4*ESpa[i]);
    ESsb1a[in]	+=db*(0.6+1.2*ESpa[i]);
    ESsb2a[in]	+=db*2.4*ESsa[i];
#endif

    in++;
  }
  ki	=0;
  k	=0; 
  for(i=0; i < ESNa1; i++){
    db	=2./(ESpa[ESNa]-ESpa[0]);
    rcT2a[ki]	+=db*drPlVc[k];
    rsT2a[ki]	+=db*drPlVs[k];
    db	*=ESsa[i];
    rcT1a[ki]	+=db*drPlVc[k];
    rsT1a[ki]	+=db*drPlVs[k];
    db	=(ESpa[i]-ESpa[0])/(ESpa[ESNa]-ESpa[0]);
    rcT[ki]	+=db*drPlVc[k];
    rsT[ki]	+=db*drPlVs[k];
    ki++;
  }
  k=1;
  for(i=0; i < ESNa1; i++){
    rcT[ki]	+=ESsa[i]*drPlVc[k];
    rsT[ki]	+=ESsa[i]*drPlVs[k];
    rcT1a[ki]	+=drPlVc[k];
    rsT1a[ki]	+=drPlVs[k];
    ki++;
  }
  k=2;
  for(i=0; i < ESNa1; i++){
    rcT[ki]	+=ESpa[i]*drPlVc[k];
    rsT[ki]	+=ESpa[i]*drPlVs[k];
    rcT1a[ki]	+=2.*ESsa[i]*drPlVc[k];
    rsT1a[ki]	+=2.*ESsa[i]*drPlVs[k];
    rcT2a[ki]	+=2.*drPlVc[k];
    rsT2a[ki]	+=2.*drPlVs[k];
    A[i]	=ESp3a[i];
    ki++;
  }
  k=3;
  if(k < MPlVf){
    for(i=0; i < ESNa1; i++){
      rcT[ki]	+=ESp3a[i]*drPlVc[k];
      rsT[ki]	+=ESp3a[i]*drPlVs[k];
      rcT1a[ki]	+=3.*ESpa[i]*drPlVc[k];
      rsT1a[ki]	+=3.*ESpa[i]*drPlVs[k];
      rcT2a[ki]	+=6.*ESsa[i]*drPlVc[k];
      rsT2a[ki]	+=6.*ESsa[i]*drPlVs[k];
      ki++;
    }
  }  
  for(k=4; k < MPlVf; k++){
    dgc		=A[0]*k;
    dgs		=dgc*drPlVs[k];
    dgc		*=drPlVc[k];
    kk		=ki;
    for(i=0; i < ESNa1; i++){
      A[i]	*=ESsa[i];
      rcT[ki]	+=A[i]*drPlVc[k];
      rsT[ki]	+=A[i]*drPlVs[k];
      ki++;
    }
    db		=rcT1a[kk+ESNa]+k*drPlVc[k];
    splAA(rcT+kk,rcT1a+kk,rcT2a+kk,&dgc,&db);
    db		=rsT1a[kk+ESNa]+k*drPlVs[k];
    splAA(rsT+kk,rsT1a+kk,rsT2a+kk,&dgs,&db);
  }
  n	=iT+ESNLt;
  while(n < ESNt1){
    ESaR0[n]	=R0;
    ESaZ0[n]	=Z0;
    in		=ESNa1*iT;
    in1		=ESNa1*n;
    for(i=0; i < ESNa1; i++){
      ESsb[in1]	=ESsb[in];
      ESsb1a[in1]	=ESsb1a[in];
      ESsb2a[in1]	=ESsb2a[in];
      in++;
      in1++;
    }
    ki	=ESnAF*iT;
    kk	=ESnAF*n;
    for(m=0; m < MPlVf; m++){
      for(i=0; i < ESNa1; i++){
	rcT[kk]=rcT[ki];
	rsT[kk]=rsT[ki];
	rcT1a[kk]=rcT1a[ki];
	rsT1a[kk]=rsT1a[ki];
	rcT2a[kk]=rcT2a[ki];
	rsT2a[kk]=rsT2a[ki];
	ki++;
	kk++;
      }
    }
    n	+=ESNLt;
  }
  return(0);
}

int ESReal2RefGeom(int iT)
{
  int k,ki;
  
  R0		=ESaR0[iT];
  Z0		=ESaZ0[iT];
  ki		=ESNa1*ESFp1*iT+ESNa;
  bPlV		=ESsb[ESNa1*iT+ESNa];
  rPlVc[0]	=R0+rcT[ki];
  rPlVs[0]	=Z0+rsT[ki];
  for(k=1; k < MPlVf; k++){
    ki	+=ESNa1;
    rPlVc[k]	=rcT[ki];
    rPlVs[k]	=rsT[ki];
  }
  return(0);
}

int ESEraseCoord()
{
  int i,k;

  for(i=0; i < ESnAP; i++){
    ECr[i]	=0.;
    ECz[i]	=0.;
  }

  for(i=0; i < ESnAT; i++){ 
    ESsb[i]	=ESsa[i];
    ESsb1a[i]	=1.;
    ESsb2a[i]	=0.;
  }
  for(i=0; i < ESnAPT; i++){
    ESsr[i]	=0.;
    ESsz[i]	=0.;
  }
  for(i=0; i < NPlVf; i++){
    rPlVc[i]	=0.;
    rPlVs[i]	=0.;
    drPlVc[i]	=0.;
    drPlVs[i]	=0.;
  }

  k	=ESnAF*ESNt1;
  for(i=0; i < k; i++){
    rcT[i]	=0.;
    rcT1a[i]	=0.;
    rcT2a[i]	=0.;
    rsT[i]	=0.;
    rsT1a[i]	=0.;
    rsT2a[i]	=0.;
  }
  return(0);
}

int ESInit3DCoord()
{
  int n,n1,i,j,ji,in,ii,m,jm,k,kk;
  double *prc,*prc1a,*prc2a,*prs,*prs1a,*prs2a;

  ESEraseCoord();
  if(nPsrf > 1){
    int im,nP1,nP2,nP3,in,id,jd,nAD;
    double *b,*b2,*A,*AA,EZga0,EZga1,EZga2,EZga3,S,s,rc,rs;
    EZga0	=0.;
    EZga1	=1e-9;
    EZga1	=1e-3;
    EZga2	=0.;
    EZga3	=1e-11;
    EZga3	=1e-3;
    
    nP1	=nPsrf+1;
    nP2	=nP1+nP1;
    nP3	=nP2+nP1;
    A	=aCr;
    AA	=aCr+nP1;
    b	=aCr+nP2;
    b2	=aCr+nP3;
    
    nAD	=nPlVd*ESNa1;
    
    for(n=0; n < ESNLt; n++){
      id	=nAD*n;
      jd	=nPlVd*n;
      for(i=0; i < nPsrf; i++){
	ESPncr2PlVData(n,i);
	for(j=0; j < nPlVd; j++){
	  rPncrDT[id+ESNa1*j]	=rPlVdT[jd+j];
	  zPncrDT[id+ESNa1*j]	=zPlVdT[jd+j];
	}
	id++;
      }
      n1	=n+ESNLt;
      while(n1 < ESNt1){
	id	=nAD*n;
	ii	=nAD*n1;
	for(j=0; j < nPlVd; j++){
	  ji	=id+ESNa1*j;
	  jm	=ii+ESNa1*j;
	  for(i=0; i < nPsrf; i++){
	    rPncrDT[jm]	=rPncrDT[ji];
	    zPncrDT[jm]	=zPncrDT[ji];
	    ji++;
	    jm++;
	  }
	}
	n1	+=ESNLt;
      }
    }
    ESRext	=0.;
    for(n=0; n < ESNLt; n++){
      m		=ESnAF*n;
      prc	=rcT+m;
      prc1a	=rcT1a+m;
      prc2a	=rcT2a+m;
      prs	=rsT+m;
      prs1a	=rsT1a+m;
      prs2a	=rsT2a+m;
      rPlVc[0]	=ESaR0[n];
      rPlVs[0]	=ESaZ0[n];
      for(m=1; m < MPlVf; m++){
	rPlVc[m]	=0.;
	rPlVs[m]	=0.;
      }
      A[0]	=0.;
      AA[0]	=0.;
      b[0]	=0.;
      im	=0;
      for(m=0; m < MPlVf; m++){
	rCc[im]	=0.;
	rCs[im]	=0.;
	im	+=nP1;
      }

      id	=nAD*n;
      jd	=nPlVd*n;
      for(i=1; i < nP1; i++){
	for(j=0; j < nPlVd; j++){
	  rPlVd[j]	=rPncrDT[id+ESNa1*j];
	  zPlVd[j]	=zPncrDT[id+ESNa1*j];
	  if(i == nPsrf){
	    rPlVdT[jd+j]	=rPlVd[j];
	    zPlVdT[jd+j]	=zPlVd[j];
	  }
	}
	id++;
	ESPlVdata2PlV();
	AA[i]	=bPlV*rPlVc[1];
	b[i]	=bPlV;
	im	=i;
	for(m=0; m < MPlVf; m++){
	  rCc[im]	=rPlVc[m];
	  rCs[im]	=rPlVs[m];
	  im	+=nP1;
	}
      }
      S		=1./AA[nPsrf];
      for(i=1; i < nP1; i++){
	AA[i]	*=S;
	A[i]	=sqrt(AA[i]);
      }
      A[nPsrf]		=1.;
      AA[nPsrf]		=1.;
      S	=1./(AA[2]-AA[1]);
      ESaR0[n]	=(AA[2]*rCc[1]-AA[1]*rCc[2])*S;
      ESaZ0[n]	=(AA[2]*rCs[1]-AA[1]*rCs[2])*S;
      for(i=1; i < nP1; i++){
	rCc[i]	-=ESaR0[n];
	rCs[i]	-=ESaZ0[n];
      }
      S	=0.;
      EZf2spl(b,b2,&S,NULL,b,A,A,&nPsrf,&EZga0,&EZga1,&EZga2,&EZga3);
      EZf2spl(rCc+nP1,rCc2a+nP1,&S,NULL,rCc+nP1,A,A,&nPsrf,&EZga0,&EZga1,&EZga2,&EZga3);
      EZf2spl(rCs+nP1,rCs2a+nP1,&S,NULL,rCs+nP1,A,A,&nPsrf,&EZga0,&EZga1,&EZga2,&EZga3);
      EZf2spl(rCc,rCc2a,&S,NULL,rCc,AA,A,&nPsrf,&EZga0,&EZga1,&EZga2,&EZga3);
      EZf2spl(rCs,rCs2a,&S,NULL,rCs,AA,A,&nPsrf,&EZga0,&EZga1,&EZga2,&EZga3);
      for(m=2; m < MPlVf; m++){
	im	=m*nP1;
	EZf2spl(rCc+im,rCc2a+im,&S,NULL,rCc+im,AA,A,&nPsrf,&EZga0,&EZga1,&EZga2,&EZga3);
	EZf2spl(rCs+im,rCs2a+im,&S,NULL,rCs+im,AA,A,&nPsrf,&EZga0,&EZga1,&EZga2,&EZga3);
	for(i=1; i < nP1; i++){
	  AA[i]	*=A[i];
	}
      }
      in	=ESNa1*n;
      i		=0;
      {
	SetIspl(A,ESsa[i],nPsrf);
	ESsb[in]	=0.;
	ESsb1a[in]=b[0];
	k	=i;
	ii	=0;
	m	=0;
	prc[k]	=0.;
	prs[k]	=0.;
	prc1a[k]=0.;
	prs1a[k]=0.;
	ii	+=nP1;
	k	+=ESNa1;
	m	=1;
	prc[k]	=0.;
	prs[k]	=0.;
	prc1a[k]=rCc[ii];
	prs1a[k]=rCs[ii];
	ii	+=nP1;
	k	+=ESNa1;
	for(m=2; m < MPlVf; m++){
	  prc[k]*=0.;
	  prs[k]*=0.;
	  prc1a[k]*=0.;
	  prs1a[k]*=0.;
	  k	+=ESNa1;
	}
      }
      in++;
      for(i=1; i < ESNa1; i++){
	SetIspl(A,ESsa[i],nPsrf);
	splR(ESsb+in,ESsb1a+in,b,b2);
	ESsb1a[in]	=ESsb[in]+ESsa[i]*ESsb1a[in];
	ESsb[in]		*=ESsa[i];
	k	=i;
	ii	=0;
	m	=0;
	splR(prc+k,prc1a+k,rCc+ii,rCc2a+ii);
	splR(prs+k,prs1a+k,rCs+ii,rCs2a+ii);
	prc[k]	*=ESpa[i];
	prs[k]	*=ESpa[i];
	ii	+=nP1;
	k	+=ESNa1;
	S	=ESsa[i];
	for(m=1; m < MPlVf; m++){
	  splR(prc+k,prc1a+k,rCc+ii,rCc2a+ii);
	  splR(prs+k,prs1a+k,rCs+ii,rCs2a+ii);
	  prc[k]*=S;
	  prs[k]*=S;
	  ii	+=nP1;
	  k	+=ESNa1;
	  S	*=ESsa[i];
	}
	in++;
      }
      in	=ESNa1*n;
      splAA2(ESsb+in,ESsb1a+in,ESsb2a+in,ESsa[0],ESsb1a+in+ESNa);
      k		=0;
      m		=0;
      kk	=k+ESNa;
      prc1a[k]	=0.;
      prc1a[kk]	+=2.*prc[kk];
      splAA(prc+k,prc1a+k,prc2a+k,prc1a+k,prc1a+kk);
      prs1a[k]	=0.;
      prs1a[kk]	+=2.*prs[kk];
      splAA(prs+k,prs1a+k,prs2a+k,prs1a+k,prs1a+kk);
      k		+=ESNa1;
      m		=1;
      kk	=k+ESNa;
      prc1a[kk]	+=prc[kk];
      splAA(prc+k,prc1a+k,prc2a+k,prc1a+k,prc1a+kk);
      prs1a[kk]	+=prs[kk];
      splAA(prs+k,prs1a+k,prs2a+k,prs1a+k,prs1a+kk);
      k	+=ESNa1;
      for(m=2; m < MPlVf; m++){
	kk	=k+ESNa;
	prc1a[k]=0.;
	prc1a[kk]+=m*prc[kk];
	splAA(prc+k,prc1a+k,prc2a+k,prc1a+k,prc1a+kk);
	prs1a[k]=0.;
	prs1a[kk]+=m*prs[kk];
	splAA(prs+k,prs1a+k,prs2a+k,prs1a+k,prs1a+kk);
	k	+=ESNa1;
      }
      ESRealGeom(n);
      ESRext	+=ESaR0[n];
    }
    ESRext	/=ESNLt;
    ESBt	=1.;
    ESRBt	=ESBt*ESRext;
  }
  else{
    for(n=0; n < ESNLt; n++){
      ESReal2RefData(n);
      ESReal2RefGeom(n);
      switch(FlagPoints){
      case 0:
	ESPlVdata2PlV();
	break;
      case 1:
	ESPlVPoints2PlV(EcRjet,EcZjet,EcNjet);
	break;
      case 2:
	ESPlVSamePlV();
	break;
      }
      ESRef2RealGeom(n);
      ESRealGeom(n);
    }
  }
  for(n=0; n < ESNLt; n++){
    n1	=n+ESNLt;
    while(n1 < ESNt1){
      ESaR0[n1]	=ESaR0[n];
      ESaZ0[n1]	=ESaZ0[n];
      in	=ESNa1*n;
      ii	=ESNa1*n1;
      kk	=ESnAF*n1;
      for(i=0; i < ESNa1; i++){
	ESsb[ii]	=ESsb[in];
	ESsb1a[ii]=ESsb1a[in];
	ESsb2a[ii]=ESsb2a[in];
	in++;
	ii++;
      }
      k		=0;
      kk	=ESnAF*n1;
      for(m=0; m < MPlVf; m++){
	for(i=0; i < ESNa1; i++){
	  rcT[kk]	=prc[k];
	  rcT1a[kk]	=prc1a[k];
	  rcT2a[kk]	=prc2a[k];
	  rsT[kk]	=prs[k];
	  rsT1a[kk]	=prs1a[k];
	  rsT2a[kk]	=prs2a[k];
	  k++;
	  kk++;
	}
      }
      n1	+=ESNLt;
    }
  }
  return(0);
}
#endif

#ifndef mzl_ESPlV
int ESPlVdata2PlV()
{
  int i,j,k,M,m;
  double t,s;

  M	=2*MPlVf-1;
  bPlV	=EZcr2*(zPlVd[0]-zPlVd[1]);
  Z0PlV	=EZcr2*(zPlVd[0]+zPlVd[1]);
  R0PlV	=EZcr2*(rPlVd[2]+rPlVd[3]);
  for(i=0; i < M; i++){
    j	=neq*i;
    for(k=0; k < M; k++){
      ac[j]=0.;
      j++;
    }
    k	=(i+1)/2;
    k	*=k;
    ac[neq*i+i]=(double)(k*k);
    f[i]= 0;
  }
  for(i=M; i < neq; i++){
    j	=neq*i+M;
    for(k=M; k < neq; k++){
      ac[j]= 0.;
      j++;
    }
  }
  for(k=0; k < kPlVd; k++){
    f[M+k]=rPlVd[k]-R0PlV;
    t	=(zPlVd[k]-Z0PlV)/bPlV;
    if(t > 1.)      t=1.;
    if(t < -1.)
      t=-1.;
    if(t > 0.){
      s=(rPlVd[0]-R0PlV)*(zPlVd[k]-Z0PlV)-
	(zPlVd[0]-Z0PlV)*(rPlVd[k]-R0PlV);
    }
    else{
      s=(rPlVd[k]-R0PlV)*(zPlVd[1]-Z0PlV)-
	(zPlVd[k]-Z0PlV)*(rPlVd[1]-R0PlV);
    }
    if(s > 0.){
      t=EZcgp-asin(t);
    }
    else{
      t=asin(t);
    }
    j=neq*(M+k);
    for(m=0; m < MPlVf; m++){
      if(m){
	ac[j]=2.*cos(m*t);
	ac[neq*(2*m-1)+M+k]=ac[j];
	j++;
	ac[j]=2.*sin(m*t);
	ac[neq*(2*m)+M+k]=ac[j];
	j++;
      }
      else{
	ac[j]=1.;
	ac[M+k]=1.;
	j++;
      }
    }
  }
  LUdcmp(ac,neq,indx,&s);
  LUbksb(ac,neq,indx,f);
  
  R0PlV	+=f[0];
  drPlVc[0]=R0PlV-rPlVc[0];
  drPlVs[0]=Z0PlV-rPlVs[0];
  rPlVc[0]=R0PlV;
  rPlVs[0]=Z0PlV;
  i	=1;
  for(m=1; m < MPlVf; m++){
    drPlVc[m]=f[i]-rPlVc[m];
    rPlVc[m]=f[i];
    i++;
    drPlVs[m]=f[i]-rPlVs[m];
    rPlVs[m]=f[i];
    i++;
  }
  return(0);
}

int ESMakePlVdata()
{
  int Fl;
  int i,j,k,ki,kj;

  MPlVf	=3;
  i	=ESNa+2*ESNa1;
  for(k=3; k < ESFp1; k++){
    i		+=ESNa1;
    if(fabs(rcT[i])+fabs(rsT[i]) > 1e-10) MPlVf=k+1;
  }
  kPlVd	=MPlVf > 6 ? 12: 2*(MPlVf-1);

  for(i=0; i < kPlVd; i++) k2kT[i]=1;
  while(i < 12){
    k2kT[i]	=0;
    i++;
  }

  Fl	=0;
  if(NPlVf < MPlVf){
    ESReInit1DFPlV();
    Fl	=1;
  }
  nDFT	=MPlVf*ESNt1;
  if(NDFT < nDFT){
    ESReInit2DFTPlV();
    Fl	=1;
  }
  if(nPlVf != MPlVf) Fl=1;
  nPlVcs	=ESNp1*MPlVf;
  if(NPlVcs < nPlVcs){
    ESReInit2DFPPlV();
    Fl	=1;
  }
  neq	=2*MPlVf-1+kPlVd;
  if(Neq < neq) ESReInitPlVeq(neq,&ac,&f,&indx);
  if(Fl || nPlVf != MPlVf) ESInit2DFPPlV();
  i	=ESNa;

  rPlVc[0]	=ESaR0[0]+rcT[i];
  rPlVs[0]	=ESaZ0[0]+rsT[i];
  for(k=1; k < MPlVf; k++){
    i		+=ESNa1;
    rPlVc[k]	=rcT[i];
    rPlVs[k]	=rsT[i];
  }
  
  EcNjet	=ESNp1;

  FlagPoints	=2;
  for(j=0; j < ESNp1; j++){
    ki		=ESNa;
    EcRjet[j]	=ESaR0[0]+rcT[ki];
    EcZjet[j]	=ESaZ0[0]+rsT[ki]+ESsb[ki]*ESsn1[j];
    kj	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      kj	+=j;
      if(kj >= ESNp) kj-=ESNp;
      EcRjet[j]	+=2.*(rcT[ki]*EScs1[kj]+rsT[ki]*ESsn1[kj]);
    }
  }
#ifdef H
  for(i=0; i < 12; i++){
    switch(i){
    case 0:
      j	=ESNp/4;
      break;
    case 1:
      j	=3*ESNp/4;
      break;
    case 2:
      j	=ESNp/2;
      break;
    case 3:
      j	=0;
      break;
    case 4:
      j	=ESNp/4-ESNp/12;
      break;
    case 5:
      j	=3*ESNp/4+ESNp/12;
      break;
    case 6:
      j	=ESNp/4+ESNp/12;
      break;
    case 7:
      j	=3*ESNp/4-ESNp/12;
      break;
    case 8:
      j	=ESNp/12;
      break;
    case 9:
      j	=ESNp-ESNp/12;
      break;
    case 10:
      j	=ESNp/2-ESNp/12;
      break;
    case 11:
      j	=ESNp/2+ESNp/12;
      break;
    }
    zPlVdT[i]	=rPlVs[0]+ESsb[ESNa]*ESsn1[j];
    rPlVdT[i]	=rPlVc[0];
    ki	=ESNa;
    kj	=0;
    for(k=1; k < ESFp1; k++){
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      ki	+=ESNa1;
      rPlVdT[i]	+=2.*(rcT[ki]*EScs1[kj]+rsT[ki]*ESsn1[kj]);
    }
  }
#endif
  return(0);
}

int ESPlVSamePlV()
{
  int i,k;
  R0PlV		=ESaR0[0]+rcT[ESNa];
  Z0PlV		=ESaZ0[0]+rsT[ESNa];
  drPlVc[0]	=0.;
  drPlVs[0]	=0.;
  rPlVc[0]	=R0PlV;
  rPlVs[0]	=Z0PlV;
  i	=ESNa;
  for(k=1; k < MPlVf; k++){
    i	+=ESNa1;
    drPlVc[k]	=0.;
    rPlVc[k]	=rcT[i];
    drPlVs[k]	=0.;
    rPlVs[k]	=rsT[i];
  }
  return(0);
}

int ESPlVPoints2PlV(double *Rb,double *Zb,int Nb)
{
  int i,j,k,kk,M,N;
  double t,s,ss,ga;
  int jE[4];
  double z[4],r[4],tE[4];
  double Cs[256],Sn[256],Cs1[256],Sn1[256],T[256];
  double *paC,*paS;

  j	=Nb-1;
  if(Rb[j] == Rb[0] && Zb[j] == Zb[0]) Nb--;
  ga	=1e-10;
  M	=2*MPlVf-1;
  N	=M+4;
  ESReInitPlVeq(N,&ac,&f,&indx);

  jE[0]	=0;
  jE[1]	=0;
  jE[2]	=0;
  jE[3]	=0;
  z[0]	=Zb[0];
  z[1]	=z[0];
  r[2]	=Rb[0];
  r[3]	=r[2];
  for(j=1; j < Nb; j++){
    if(z[0] < Zb[j]){
      z[0]	=Zb[j];
      jE[0]	=j;
    }
    if(z[1] > Zb[j]){
      z[1]	=Zb[j];
      jE[1]	=j;
    }
    if(r[2] > Rb[j]){
      r[2]	=Rb[j];
      jE[2]	=j;
    }
    if(r[3] < Rb[j]){
      r[3]	=Rb[j];
      jE[3]	=j;
    }
  }

  for(kk=0; kk < 2; kk++){
    j	=jE[kk];
    i	=j ? j-1 : Nb-1;
    k	=j+1;
    if(k == Nb) k=0;
    s	=EZcr2*(Zb[k]-Zb[i]);
    t	=Zb[k]+Zb[i]-2.*Zb[j];
    if(t != 0.){
      t	=-s/t;
      z[kk]	+=EZcr2*t*s;
      r[kk]	=Rb[j]+EZcr2*((Rb[k]-Rb[i])+(Rb[k]+Rb[i]-2.*Rb[j])*t)*t;
    }
    else{
      r[kk]	=Rb[j];
      z[kk]	=Zb[j];
    }
    rPlVdT[kk]	=r[kk];
    zPlVdT[kk]	=z[kk];
  }
  for(kk=2; kk < 4; kk++){
    j	=jE[kk];
    i	=j ? j-1 : Nb-1;
    k	=j+1;
    if(k == Nb) k=0;
    s	=EZcr2*(Rb[k]-Rb[i]);
    t	=Rb[k]+Rb[i]-2.*Rb[j];
    if(t != 0.){
      t	=-s/t;
      r[kk]	+=EZcr2*t*s;
      z[kk]	=Zb[j]+EZcr2*((Zb[k]-Zb[i])+(Zb[k]+Zb[i]-2.*Zb[j])*t)*t;
    }
    else{
      r[kk]	=Rb[j];
      z[kk]	=Zb[j];
    }
    rPlVdT[kk]	=r[kk];
    zPlVdT[kk]	=z[kk];
  }
  zPlVd[0]	=z[0];
  zPlVd[1]	=z[1];
  bPlV	=EZcr2*(z[0]-z[1]);
  Z0PlV	=EZcr2*(z[0]+z[1]);
  R0PlV	=EZcr2*(r[3]+r[2]);
  t	=0.05*(r[3]-r[2]);
  for(j=0; j < Nb; j++){
    t	=(Zb[j]-Z0PlV)/bPlV;
    if(t > 1.){
      t=1.;
    }
    if(t < -1.){
      t=-1.;
    }
    if(t > 0.){
      s=(r[0]-R0PlV)*(Zb[j]-Z0PlV)-(z[0]-Z0PlV)*(Rb[j]-R0PlV);
    }
    else{
      s=(Rb[j]-R0PlV)*(z[1]-Z0PlV)-(Zb[j]-Z0PlV)*(r[1]-R0PlV);
    }
    if(s > 0.){
      t=EZcgp-asin(t);
    }
    else{
      t=asin(t);
    }
    T[j]	=t;
  }
  tE[2]	=EZcgp-asin((z[2]-Z0PlV)/bPlV);
  tE[3]	=asin((z[3]-Z0PlV)/bPlV);
  
  ac[0]	=Nb;
  s	=0.;
  for(j=0; j < Nb; j++){
    s		+=Rb[j];
  }
  f[0]	=s;
  for(kk=1; kk < MPlVf; kk++){
    s	=0.;
    ss	=0.;
    for(j=0; j < Nb; j++){
      t		=kk*T[j];
      s		+=cos(t);
      ss	+=sin(t);
    }
    ac[2*kk-1]	=2.*s;
    ac[2*kk]	=2.*ss;
  }
  j	=M;
  ac[j++]	=1.;
  ac[j++]	=1.;
  ac[j++]	=1.;
  ac[j]		=1.;
  for(k=1; k < MPlVf; k++){
    i	=2*k-1;
    paC	=ac+N*i;
    paS	=ac+N*(i+1);
    s	=0.;
    ss	=0.;
    for(j=0; j < Nb; j++){
      t	=k*T[j];
      Cs[j]	=2.*cos(t);
      Sn[j]	=2.*sin(t);
      s		+=Rb[j]*Cs[j];
      ss	+=Rb[j]*Sn[j];
    }
    f[i]	=s;
    f[i+1]	=ss;
    kk	=k;
    t	=ga*k*k*k*k;
    s	=t;
    ss	=t;
    t	=0.;
    for(j=0; j < Nb; j++){
      s		+=Cs[j]*Cs[j];
      t		+=Cs[j]*Sn[j];
      ss	+=Sn[j]*Sn[j];
    }
    paC[i++]	=s;
    paC[i]	=t;
    paS[i++]	=ss;
    for(kk++; kk < MPlVf; kk++){
      s	=0.;
      ss=0.;
      for(j=0; j < Nb; j++){
	t	=cos(kk*T[j]);
	s	+=Cs[j]*t;
	ss	+=Sn[j]*t;
      }
      paC[i]	=2.*s;
      paS[i++]	=2.*ss;
      s	=0.;
      ss=0.;
      for(j=0; j < Nb; j++){
	t	=sin(kk*T[j]);
	s	+=Cs[j]*t;
	ss	+=Sn[j]*t;
      }
      paC[i]	=2.*s;
      paS[i++]	=2.*ss;
    }
    t	=k*EZcr2*EZcgp;
    paC[i]	=2.*cos(t);
    paS[i++]	=2.*sin(t);
    t	=-k*EZcr2*EZcgp;
    paC[i]	=2.*cos(t);
    paS[i++]	=2.*sin(t);
    t	=k*tE[2];
    paC[i]	=2.*cos(t);
    paS[i++]	=2.*sin(t);
    t	=k*tE[3];
    paC[i]	=2.*cos(t);
    paS[i]	=2.*sin(t);
  }
  for(i=M; i < N; i++){
    for(j=i; j < N; j++){
      ac[i*N+j]	=0.;
    }
    f[i]	=r[i-M];
  }
  for(i=0; i < N; i++){
    for(j=0; j < i; j++){
      ac[N*i+j]	=ac[N*j+i];
    }
  }

  LUdcmp(ac,N,indx,&s);
  LUbksb(ac,N,indx,f);
#ifdef H
  if(CholDc(ac,N) == -1){
    for(k=0; k < N; k++){
      printf("CholDc failed ac[%2d]=%10.3e\n",k,ac[k*N+k]);
    }
  }
  CholBksb(ac,N,f);
#endif
  R0PlV		=f[0];
  drPlVc[0]	=R0PlV-rPlVc[0];

  if(ESaR0[0] < rPlVdT[2] || rPlVdT[3] < ESaR0[0]){
    t	=EZcr2*(rPlVdT[2]+rPlVdT[3])-ESaR0[0];
    ESaR0[0]	+=t;
    drPlVc[0]	-=t;
    R0	+=t;
  }

  drPlVs[0]	=Z0PlV-rPlVs[0];
  rPlVc[0]	=R0PlV;
  rPlVs[0]	=Z0PlV;
  i	=1;
  for(k=1; k < MPlVf; k++){
    drPlVc[k]	=f[i]-rPlVc[k];
    rPlVc[k]	=f[i++];
    drPlVs[k]	=f[i]-rPlVs[k];
    rPlVs[k]	=f[i++];
  }
  return(0);
}

int ESPlVdata2PlVX()
{
  static double t,s,dW;
  static double gg_x;
  static int K,n,nn,i,j,k,kk,ki,ndata,neq;
  static double ac[NEQ*NEQ],f[NEQ],f0[NEQ],f1[NEQ],f2[NEQ];
  static int index[NEQ];

  ndata	=kPlVd+2;
  if(2*nf_X < ndata+1){
    nf_X=(ndata+2)/2;
  }
  K	=2*nf_X;
  neq	=K+ndata;
  if(neq > NEQ){
    nf_X=(NEQ-ndata)/2;
    K	=2*nf_X;
    neq	=K+ndata;
  }  

  if(ESaZx > 0.){
    k2kPlV[0]	=0;
    k2kPlV[1]	=1;
    t_X		=EZcr2*EZcgp;
    singtX	=1.;
  }
  else{
    k2kPlV[0]	=1;
    k2kPlV[1]	=0;
    t_X		=-EZcr2*EZcgp;
    singtX	=-1.;
  }	
  z_X0		=EZcr2*(zPlVd[0]+zPlVd[1]);
  b_X		=EZcr2*(zPlVd[0]-zPlVd[1]);

  /*\bf Matrix formation*/
  for(k=0; k < K; k++){
    j	=neq*k;
    for(kk=0; kk < K; kk++){
      ac[j]	=0.;
      j++;
    }
    j	=k/2;
    ac[neq*k+k]=(double)(j*j*j*j);
    f0[k]	=0.;
    f1[k]	=0.;
    f2[k]	=0.;
  }
  for(n=K; n < neq; n++){
    j	=neq*n+K;
    for(nn=K; nn < neq; nn++){
      ac[j]	=0.;
      j++;
    }
  }
  for(n=0; n < kPlVd; n++){
    i		=k2kPlV[n];
    f0[K+n]	=rPlVd[i]-ESaRx;
    f1[K+n]	=0.;
    f2[K+n]	=0.;
    if(n){
      t=(zPlVd[i]-z_X0)/b_X;
      if(t > 1.)
	t=1.;
      if(t < -1.)
	t=-1.;
      s=(rPlVd[i]-rPlVd[0])*(zPlVd[1]-zPlVd[0])-
	(rPlVd[1]-rPlVd[0])*(zPlVd[i]-zPlVd[0]);
      t=s < 0 ? asin(t) : EZcgp-asin(t);
      s=fabs(sin(EZcr2*(t_X-t)));
      j	=neq*(K+n);
      for(k=0; k < nf_X; k++){
	if(k){
	  ac[j]		=2.*cos(k*t)*s;
	  ac[j+1]	=2.*sin(k*t)*s;
	}
	else{
	  ac[j]	=s;
	  ac[j+1]	=sin(t)-singtX;
	}
	ac[neq*(2*k)+K+n]	=ac[j];
	ac[neq*(2*k+1)+K+n]	=ac[j+1];
	j +=2;
      }
    }
    else{
      t	=t_X;
      j=neq*K;
      for(k=0; k < nf_X; k++){
	if(k){
	  ac[j]	=2.*cos(k*t);
	  ac[j+1]	=2.*sin(k*t);
	}
	else{
	  ac[j]	=1.;
	  ac[j+1]	=0.;
	}
	ac[neq*(2*k)+K]	=ac[j];
	ac[neq*(2*k+1)+K]	=ac[j+1];
	j +=2;
      }
    }
  }
  i	=neq-2;
  f0[i]	=0.;
  f1[i]	=-b_X;
  f2[i]	=0.;
  j	=neq*i;
  for(k=0; k < K; k++){
    ac[j]	=k == 1 ? 1. : 0.;
    ac[neq*k+i]=ac[j];
    j++;
  }
  i++;
  f0[i]	=0.;
  f1[i]	=0.;
  f2[i]	=-b_X*singtX;
  t	=t_X;
  j=neq*i;
  for(k=0; k < nf_X; k++){
    if(k){
      ac[j]	=-2.*k*sin(k*t);
      ac[j+1]	=2.*k*cos(k*t);
    }
    else{
      ac[j]	=0.;
      ac[j+1]	=0.;
    }
    ac[neq*(2*k)+i]	=ac[j];
    ac[neq*(2*k+1)+i]	=ac[j+1];
    j +=2;
  }
  LUdcmp(ac,neq,index,&s);
  LUbksb(ac,neq,index,f0);
  LUbksb(ac,neq,index,f1);
  LUbksb(ac,neq,index,f2);
  ac[0]	=0.3;
  ac[1]	=-0.3;
  n=0;
  do{
    EZzero(&n,ac,ac+2,&gg_x,dW);
    s	=sqrt(1.+gg_x*gg_x);
    dW	=0.;
    for(k=0,i=0; k < nf_X; k++){
      f[i]	=f0[i]+f1[i]*gg_x+f2[i]*s;
      dW +=k*k*k*k*f[i]*(f1[i]+f2[i]*gg_x/s);
      i++;
      f[i]	=f0[i]+f1[i]*gg_x+f2[i]*s;
      dW +=k*k*k*k*f[i]*(f1[i]+f2[i]*gg_x/s);
      i++;
    }
  }while(n > 0 && n < 15);
  for(k=0,i=0; k < nf_X; k++){
    EZdx0cs[k]	=f[i]-EZx0cs[k];
    EZx0cs[k]	=f[i];
    i++;
    EZdx0sn[k]	=f[i]-EZx0sn[k];
    EZx0sn[k]	=f[i];
    i++;
  }
  return(0);
}

static double rBL[1500],zBL[1500];
int ESGapCoreSepX()
{
  static int i,j,k,ki,ji,kj,iX,iO;
  static double x_0,dx_0,dgdx_0,d2x_0,d2gdx_0;
  static double cs,sn,d,s,sqs,x_1,x_2,z_1;

#ifdef H
  EZd_vb	=(ESsb[ESNa-ESnBL]-ESsb[ESNa])/b_X; /*$ d_{vb}=(b_{vb}-b_X)/b_X$*/
  EZd_vb	=(ESNa-ESNax)/(double)ESNax;
  EZd_vb	=-0.1;
#endif
  Z0_vb	=z_X0+b_X*(EZd_vb+EZd_vb*EZd_vb*z_X2)*singtX;
  for(i=0; i < ESnBL1; i++){
    ESDx[i]	=0.;
    ESLx[i]	=0.;
    ESVx[i]	=0.;
  }
  i	=ESNp/4;
  j	=3*ESNp/4;
  if(ESaZx > 0.){
    iX	=i;
    iO	=j;
  }
  else{
    iX	=j;
    iO	=i;
  }
  EZdx_0x		=0.;
  d2x_0x	=0.;
  EZdx_0o		=0.;
  d2x_0o	=0.;
  ki		=iX;
  kj		=iO;
  for(k=1; k < nf_X; k++){
    cs		=EScs1[ki];
    sn		=ESsn1[ki];
    ki		+=iX;
    if(ki >= ESNp){
      ki	-=ESNp;
    }
    EZdx_0x	+=k*(-EZdx0cs[k]*sn+EZdx0sn[k]*cs);
    d2x_0x	+=k*k*(EZdx0cs[k]*cs+EZdx0sn[k]*sn);
    cs		=EScs1[kj];
    sn		=ESsn1[kj];
    kj		+=iO;
    if(kj >= ESNp){
      kj	-=ESNp;
    }
    EZdx_0o	+=k*(-EZdx0cs[k]*sn+EZdx0sn[k]*cs);
    d2x_0o	+=k*k*(EZdx0cs[k]*cs+EZdx0sn[k]*sn);
  }
  EZdx_0x		*=2.*singtX;
  d2x_0x	*=-2.;
  EZdx_0o		*=2.*singtX;
  d2x_0o	*=-2.;
  for(j=0; j < ESNp1; j++){
    dx_0	=0.;
    d2x_0	=0.;
    x_0		=0.;
    dgdx_0	=0.;
    d2gdx_0	=0.;
    kj		=0;
    for(k=1; k <nf_X; k++){
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      cs	=EScs1[kj];
      sn	=ESsn1[kj];
      dx_0	+=k*(-EZx0cs[k]*sn+EZx0sn[k]*cs);
      d2x_0	+=k*k*(EZx0cs[k]*cs+EZx0sn[k]*sn);
      x_0	+=EZdx0cs[k]*cs+EZdx0sn[k]*sn;
      dgdx_0	+=k*(-EZdx0cs[k]*sn+EZdx0sn[k]*cs);
      d2gdx_0	+=k*k*(EZdx0cs[k]*cs+EZdx0sn[k]*sn);
    }
    dx_0	*=2.*singtX;
    d2x_0	*=-2.;
    x_0		*=2.;
    dgdx_0	*=2.*singtX;
    d2gdx_0	*=-2.;
    x_0		+=EZdx0cs[0];

    EZx0sep[j]	+=x_0;
    cs		=EScs1[j];
    sn		=ESsn1[j]*singtX;
    EZx1sep[j]	=EZx0sep[j]*(3.-sn)+2.*dx_0*cs;
    if(j == iX){
      x_2	=13.*d2x_0x;
    }
    else{
      if(j == iO){
	x_2	=-x_0*(2.-z_X2)+4.*d2x_0o*z_X2;
      }
      else{
	x_2	=d2gdx_0*(1.+sn)+x_0*(EZcr4*(sn-5.)-z_X2)
	  -(dgdx_0*cs+x_0)/(1.-sn)
	    +(2.*z_X2*(dgdx_0-EZdx_0o)*(1.-sn)
	      -4.*(dgdx_0-EZdx_0x)*(1.+sn))/cs;
      }
    }
    EZx2sep[j]	+=x_2;
    s		=1./(EZcr2*(1.-sn)-2.*EZd_vb);
    sqs		=sqrt(s);
    x_1		=3.*dx_0*(1.-sn*singtX)+(2.*d2x_0-EZx0sep[j])*cs;

    Dz2[j]	=(dx_0/s+EZcr4*(EZd_vb*(x_1+EZx1sep[j]*cs*EZcr4*s)-EZx0sep[j]*cs))*sqs;
    Ddx2[j]	=1.+sn*singtX;
    Dx2[j]	=(2.+EZd_vb*s)*EZd_vb*(1.+EZd_vb)*cs;
    Dx0[j]	=(1.+EZd_vb)*(-EZx0sep[j]+EZcr4*EZx1sep[j]*(1.+EZd_vb*s))*sqs*cs-
      Dz2[j]*Ddx2[j];
    Dz2[j]	*=-2.*EZd_vb;
    Ddx2[j]	*=-EZd_vb*EZd_vb;

    for(i=0; i < ESnBL1; i++){
      d		=EZd_vb*sqrt((ESnBL-i)/(double)ESnBL);
      s		=EZcr2*(1.-sn)-2.*d;
      sqs	=sqrt(s);
      if(i != ESnBL){
	x_1	=EZcr4*(EZx1sep[j]+d*EZx2sep[j])/s;
      }
      else{
	x_1	=0.;
      }
      z_1	=(d+d*d*z_X2)*singtX+(1.+d)*sn*singtX;
      ji	=j+ESNp1*i;
      rBL[ji]	=ESaRx+(EZx0sep[j]+d*x_1)*sqs+EZx0sn[0]*(z_1-singtX);
      zBL[ji]	=z_X0+b_X*z_1;

      if(j != iX){
	x_1	=EZcr2*(dx_0*(4.*z_X2*(1.-sn)-8.*(1.+sn))
		      +(2.*d2x_0*(1.+sn)-2.*EZx2sep[j]-2.*EZx0sep[j]*z_X2
			-(2.*dx_0*cs+2.*EZx0sep[j])/(1.-sn)
			-EZcr2*EZx0sep[j]*(5.-sn)
			)*cs)/sqs;
	s	=ESgt[j]-ESgt[iX];
	x_2	=16.*EZdx_0x/sqrt(s*s-8.*d);
	ESDx[i]	+=x_1+x_2;
	ESLx[i]	+=x_1/rBL[ji]+x_2/ESaRx;
	ESVx[i]	+=x_1*rBL[ji]+x_2*ESaRx;
      }	
    }
  }
  s	=1./ESNp;
  for(i=0; i < ESnBL1; i++){
    ESDx[i]	*=s;
    ESLx[i]	*=s;
    ESVx[i]	*=s;
  }
  x_1	=-16.*EZdx_0x;
  d	=EZd_vb;
  ESgFPlVx=-d*d*b_X*ESRBt*(EZcr4*x_1*(log(-2.*d/(EZcgp*EZcgp))-EZcr2)/ESaRx
			-0.125*(ESLx[0]+ESLx[1]));
  for(i=0; i < ESNp1; i++){
    EZrvb[i]=rBL[i];
    EZzvb[i]=zBL[i];
  }

  EZfour(Dx2c,Dx2s,Dx2);
  EZfour(Dz2c,Dz2s,Dz2);
  EZfour(EZrvbcs,EZrvbsn,EZrvb);
  EZrvbcs[0]	-= R0;
  
  ki	=ESNa;
  for(k=0; k < ESFp1; k++){
    EZrvbcs[k]	-=rcT[ki];
    EZrvbsn[k]	-=rsT[ki];
    ki		+=ESNa1;
  }
  EZrvbsn[0]	=0.;
  return(0);
}

int ESCoreX()
{
  int i,ki,k,j,ji,kj;
#ifdef XWIN
  double zb0,drc0,drs0,s,db,A[ESNa1];
#else
  double zb0,drc0,drs0,s,db,A[129];
#endif
  zb0	=Z0_vb-Z0-rsT[ESNa];
  drc0	=EZrvbcs[0];

  ESaR0[0]	=R0;
  ESaZ0[0]	=Z0;
  db		=b_X+EZd_vb*b_X-ESsb[ESNa];
  for(i=0; i < ESNa1; i++){
    ESsb[i]	+=ESsa[i]*db;
    ESsb1a[i]	+=db;
  }

  for(i=0; i < ESNa1; i++){
    rcT[i]	+=ESpa[i]*drc0;
    rsT[i]	+=ESpa[i]*zb0;
    s		=2.;
    rcT2a[i]	+=s*drc0;
    rsT2a[i]	+=s*zb0;
    s		*=ESsa[i];
    rcT1a[i]	+=s*drc0;
    rsT1a[i]	+=s*zb0;
  }		

  k	=1;
  ki	=ESNa1*k;
  drc0	=EZrvbcs[k];
  drs0	=EZrvbsn[k];
  for(i=0; i < ESNa1; i++){
    rcT[ki]	+=ESsa[i]*drc0;
    rsT[ki]	+=ESsa[i]*drs0;
    rcT1a[ki]	+=drc0;
    rsT1a[ki]	+=drs0;
    ki++;
  }		

  k	=2;
  drc0	=EZrvbcs[k];
  drs0	=EZrvbsn[k];
  for(i=0; i < ESNa1; i++){
    rcT[ki]	+=ESpa[i]*drc0;
    rsT[ki]	+=ESpa[i]*drs0;
    rcT1a[ki]	+=2.*ESsa[i]*drc0;
    rsT1a[ki]	+=2.*ESsa[i]*drs0;
    rcT2a[ki]	+=2.*drc0;
    rsT2a[ki]	+=2.*drs0;
    A[i]	=ESp3a[i];
    ki++;
  }

  k	=3;
  drc0	=EZrvbcs[k];
  drs0	=EZrvbsn[k];
  for(i=0; i < ESNa1; i++){
    rcT[ki]	+=ESp3a[i]*drc0;
    rsT[ki]	+=ESp3a[i]*drs0;
    rcT1a[ki]	+=3.*ESpa[i]*drc0;
    rsT1a[ki]	+=3.*ESpa[i]*drs0;
    rcT2a[ki]	+=6.*ESsa[i]*drc0;
    rsT2a[ki]	+=6.*ESsa[i]*drs0;
    ki++;
  }

  for(k=4; k < ESFp1; k++){
    ki		=ESNa1*k;
    kj		=ki;
    drc0	=EZrvbcs[k];
    drs0	=EZrvbsn[k];
    for(i=0; i < ESNa1; i++){
      A[i]	*=ESsa[i];
      rcT[ki]	+=A[i]*drc0;
      rsT[ki]	+=A[i]*drs0;
      ki++;
    }	
    db		=rcT1a[kj+ESNa]+k*drc0;
    splAA(rcT+kj,rcT1a+kj,rcT2a+kj,ESsa,&db);
    db		=rsT1a[kj+ESNa]+k*drs0;
    splAA(rsT+kj,rsT1a+kj,rsT2a+kj,ESsa,&db);
  }

#ifdef H
  for(j=0; j < ESNp1; j++){
    ESsr[j]	=R0;
    ESsz[j]	=Z0;
  }
  ji	=ESNp1;
  for(i=1; i < ESNa1; i++){
    for(j=0; j < ESNp1; j++){
      rBL[ji]	=R0+EZrcs[i];
      zBL[ji]	=Z0+EZz0[i]+ESsb[i]*ESsn1[j];
      kj	=0;
      ki	=i;
      s		=0.;
      for(k=1; k < ESFp1; k++){
	ki	+=ESNa1;
	kj	+=j;
	if(kj >= ESNp){
	  kj -=ESNp;
	}
	s	+=EZrcs[ki]*EScs1[kj]+EZrsn[ki]*ESsn1[kj];
      }
      rBL[ji]	+=2.*s;
      ji++;
    }
  }
#endif
  return(0);
}

int ESPlVdata2PlVPolar()
{
  int i,j,k,M,m;
  double t,s,dx,dy;

  M	=2*MPlVf-1;
  for(i=0; i < M; i++){
    j	=neq*i;
    for(k=0; k < M; k++){
      ac[j]=0.;
      j++;
    }
    k	=(i+1)/2;
    k	*=k;
    ac[neq*i+i]=(double)(k*k);
    f[i]= 0;
  }
  for(i=M; i < neq; i++){
    j	=neq*i+M;
    for(k=M; k < neq; k++){
      ac[j]=0.;
      j++;
    }
  }
  for(k=0; k < kPlVd; k++){
    dx	=rPlVd[k]-R0;
    dy	=zPlVd[k]-Z0;
    f[M+k]=dx*dx+dy*dy;
    if(dx > 0.){
      t	=atan(dy/dx);
      if(dy < 0.){
	t+=EZc2gp;
      }
    }
    else{
      if(dx < 0.){
	t=atan(dy/dx)+EZcgp;
      }
      else{
	t	=EZcr2*EZcgp;
	if(dy < 0.){
	  t+=EZcgp;
	}
      }
    }
    j	=neq*(M+k);
    for(m=0; m < MPlVf; m++){
      if(m){
	ac[j]=2.*cos(m*t);
	ac[neq*(2*m-1)+M+k]=ac[j];
	j++;
	ac[j]=2.*sin(m*t);
	ac[neq*(2*m)+M+k]=ac[j];
	j++;
      }
      else{
	ac[j]=1.;
	ac[M+k]=1.;
	j++;
      }
    }
  }
  LUdcmp(ac,neq,indx,&s);
  LUbksb(ac,neq,indx,f);
  
  drPlVc[0]	=f[0]-rPlVc[0];
  rPlVc[0]	=f[i];
  i	=1;
  for(m=1; m < MPlVf; m++){
    drPlVc[m]	=f[i]-rPlVc[m];
    rPlVc[m]	=f[i];
    i++;
    drPlVs[m]	=f[i]-rPlVs[m];
    rPlVs[m]	=f[i];
    i++;
  }
  return(0);
}

int ESComprPlVacData(double *prD,double *pzD,int *k2,int n)
{
  int i,j,k;
  double s;

  k	=0;
  while(k < n){
    if(k2[k] == 0){
      i		=k+1;
      while(i < n){
	if(k2[i]){
	  k2[k]	=k2[i];
	  k2[i]	=0;
	  s	=prD[k];
	  prD[k]	=prD[i];
	  prD[i]	=s;
	  s		=pzD[k];
	  pzD[k]	=pzD[i];
	  pzD[i]	=s;
	  break;
	}
	i++;
      }
    }
    if(k2[k]){
      k++;
    }
    else{
      break;
    }
  }
  return(k);
}

int ESCheckIntData()
{
  if(2*MPlVf-1 < kPlVd){
    MPlVf=kPlVd/2+1;
  }
  if(MPlVf < 3){
    MPlVf=3;
  }
  return(0);
}

int ESGetPlVDAddr(double **rdPV, double **zdPV)
{
  *rdPV	=rPlVdT;
  *zdPV	=zPlVdT;
  return(0);
}

int ESSetPlVInt(int npv,int mpv)
{
  int i;

  kPlVd	=npv;
  MPlVf	=mpv;
  for(i=0; i < kPlVd; i++){
    k2kT[i]	=1;
  }
  while(i < nPlVd){
    k2kT[i]	=0;
    i++;
  }
  ESCheckIntData();
  return(0);
}

int ESInit2DPlGeom()
{
  static int i,iT=0;
  int Fl,kP,*k2;
  double *prD,*pzD,*pr,*pz;

  Fl	=0;
  prD	=rPlVdT;
  pzD	=zPlVdT;
  k2	=k2kT;
  pr	=ESsr;
  pz	=ESsz;
  if(ESMp < 3){
    ESMp=3;
  }
  if(ESMp*4 > ESNp){
    ESMp=ESNp/4;
  }
  if(ESMp > ESFp1){
    ESMp=ESFp1;
  }
  ESMp1	=ESMp+1;
  ES2Mp	=2*ESMp;
  ES2Mp1=ES2Mp+1;
  ESnMp	=ESNa1*ES2Mp1;
  kPlVd	=ESComprPlVacData(prD,pzD,k2,nPlVd);
  ESCheckIntData();

  if(NPlVf < MPlVf){
    ESReInit1DFPlV();
    Fl=1;
  }
  nDFT	=MPlVf*ESNt1;
  if(NDFT < nDFT){
    ESReInit2DFTPlV();
    Fl=1;
  }
  if(nPlVf != MPlVf){ 
    Fl=1;
  }
  nPlVcs	=ESNp1*MPlVf;
  if(NPlVcs < nPlVcs){
    ESReInit2DFPPlV();
    Fl=1;
  }
  neq=2*MPlVf-1+kPlVd;
  if(Neq < neq){
    ESReInitPlVeq(neq,&ac,&f,&indx);
  }
  nrC	=MPlVf*(nPsrf+1);
  if(NrC < nrC){
    ESReInit2DAPrC();
  }
  if(Fl){
    ESInit2DFPPlV();
    ESInit3DCoord();
    Fl=0;
  }
  ESReal2RefData(iT);
  ESReal2RefGeom(iT);
  return(0);
}

#ifdef H
int ESSavePlVata(FILE *lf)
{
  int j,ji;

  fprintf(lf,"%9.6f %9.6f %s\n",ESBt,ESRext,
	  "- vacuum toroidal field and its reference radius");

  fprintf(lf,"%9.6f %10.6f - R and Z (approximate) of the magnetic axis\n"
	  ,ESaR0[0],ESaZ0[0]);
  fprintf(lf,"%3d - number of boundary points\n",ESNp1);
  fprintf(lf,"%3s %9s %9s\n","#","Rbound[#]","Zbound[#]");
  ji	=ESNa*ESNp1;
  for(j=0; j < ESNp1; j++){
    fprintf(lf,"%3d %9.6f %9.6f\n",j,ESsr[ji],ESsz[ji]);
    ji++;
  }
  return(0);
}
#endif

int ESUsePlVDataFile()
{
  puts("Key names in the PlVac.ai file:\n");
  puts("(order of data is not important)\n");
  puts("RBtor=...  Rext=... R0=... Z0=...\n");
  puts("Z0=...\n");
  puts("nPlV=...\n");
  puts("rPlV=...\n");
  puts("zPlV=...\n");
  return(0);
}

int ESReadInt(FILE *lf,char *Nm,int *i,int n)
{
#ifdef XWIN
  if(PFIfPatternFound(lf,Nm) == 0){
    printf("PlVac.ai: No \"%s\" pattern%c\n",Nm,'\a');
    return(1);
  }
  else{
    int j;
    for(j=0; j < n; j++){
      if(fscanf(lf,"%d",i+j) < 1){
	printf("PlVac.ai: No data for \"%s\"%c\n",Nm,'\a');
	return(2);
      }
    }
  }
#endif
  return(0);
}

int ESReadDbl(FILE *lf,char *Nm,double *a,int n)
{
#ifdef XWIN
  if(PFIfPatternFound(lf,Nm) == 0){
    printf("PlVac.ai: No \"%s\" pattern%c\n",Nm,'\a');
    return(1);
  }
  else{
    int j;
   
    for(j=0; j < n; j++){
      if(fscanf(lf,"%lg",a+j) < 1){
	printf("PlVac.ai: No data for \"%s\"%c\n",Nm,'\a');
	return(2);
      }
    }
  }
#endif
  return(0);
}

int ESWritePlVacAO(FILE *lf)
{
  int j,ji;

  fprintf(lf,"RBtor=%9.6f Rext=%9.6f %s\n",ESRBt,ESRext,
	  "- vacuum R x Btor toroidal field and its reference radius");

  fprintf(lf,"R0=%9.6f Z0=%10.6f -R and Z (approximate) of the magnetic axis\n"
	  ,ESaR0[0],ESaZ0[0]);
  fprintf(lf,"nPlV=%3d - number of boundary points\n",ESNp1);
  ji	=ESNa*ESNp1;
  fprintf(lf,"rPlV=");
  for(j=0; j < ESNp1; j++){
    fprintf(lf,"%9.6f",ESsr[ji]);
    if(j%8 == 0){
      fputc('\n',lf);
    }
    ji++;
  }
  ji	=ESNa*ESNp1;
  fprintf(lf,"zPlV=");
  for(j=0; j < ESNp1; j++){
    fprintf(lf,"%9.6f",ESsz[ji]);
    if(j%8 == 0){
      fputc('\n',lf);
    }
    ji++;
  }
  return(0);
}

int ESReadPlVacAI(FILE *lf)
{
  int j,k,L;
  char ch;

  L	=0;
  L	|=ESReadDbl(lf,"RBtor=",&ESRBt,1);
  rewind(lf);
  L	|=ESReadDbl(lf,"Rext=",&ESRext,1);
  ESBt	=ESRBt/ESRext;
  rewind(lf);
  L	|=ESReadDbl(lf,"R0=",ESaR0,1);
  rewind(lf);
  L	|=ESReadDbl(lf,"Z0=",ESaZ0,1);
  rewind(lf);
  k	=ESReadInt(lf,"nPlV=",&EcNjet,1);
  L	|=k;
  if(k == 0){
    L	|=ESReadDbl(lf,"rPlV=",EcRjet,EcNjet);
    rewind(lf);
    L	|=ESReadDbl(lf,"zPlV=",EcZjet,EcNjet);
    rewind(lf);
  }
  if(EcNjet > 12) FlagPoints	=1;
  if(L == 1)	 ESUsePlVDataFile();
  return(0);
}

#ifdef SHMEM
int ESReadPlVacShM()
{
  static int FStart=0;
  extern int *ESShMemI;
  extern double *ESShMemT;
  extern double *ESShMemB;

  ESRBt	=ESShMemT[1];
  ESRext=ESShMemT[2];
  ESBt	=ESRBt/ESRext;

  return(0);

  EcNjet=ESShMemI[5];
  if(EcNjet > 12){
    FlagPoints=1;
    memcpy((void *)EcRjet,(void *)(ESShMemB),EcNjet*sizeof(double));
    memcpy((void *)EcZjet,(void *)(ESShMemB+EcNjet),EcNjet*sizeof(double));
    if(FStart == 0){
      int i;
      zPlVdT[0]	=EcZjet[0];
      zPlVdT[1]	=EcZjet[0];
      rPlVdT[2]	=EcRjet[0];
      rPlVdT[3]	=EcRjet[0];
      for(i=1; i < EcNjet; i++){
	if(zPlVdT[0] < EcZjet[i]) zPlVdT[0]=EcZjet[i];
	if(zPlVdT[1] > EcZjet[i]) zPlVdT[1]=EcZjet[i];
	if(rPlVdT[3] < EcZjet[i]) rPlVdT[3]=EcZjet[i];
	if(rPlVdT[2] > EcZjet[i]) rPlVdT[2]=EcZjet[i];
      }
      ESaR0[0]=EZcr2*(rPlVdT[2]+rPlVdT[3]);
      ESaZ0[0]=EZcr2*(zPlVdT[0]+zPlVdT[1]);
      FStart	=1;
    }
  }
  else{
    int j;
    FlagPoints=0;
    memcpy((void *)rPlVdT,(void *)(ESShMemB),EcNjet*sizeof(double));
    memcpy((void *)zPlVdT,(void *)(ESShMemB+EcNjet),EcNjet*sizeof(double));
    for(j=0; j < EcNjet; j++){
      k2k[j]=1;
      k2kT[j]=1;
    }
    for(; j < 12; j++){
      k2k[j]=0;
      k2kT[j]=0;
    }

    if(FStart == 0){
      ESaR0[0]=EZcr2*(rPlVdT[2]+rPlVdT[3]);
      ESaZ0[0]=EZcr2*(zPlVdT[0]+zPlVdT[1]);
      FStart	=1;
    }
  }

  return(0);
}
#endif/*SHMEM*/

int ReSetVirtualBSize()
{
  int i,k,ki;
  double *lc,*ls;
#ifdef XWIN
  double s,ds,g[ESNa1],a[ESNa1],drc[ESFp1],drs[ESFp1];
#else
  double s,ds,g[129],a[129],drc[33],drs[33];
#endif 
  ds	=(ESsa[ESNa]-ESa0)/ESNa;
  for(i=0; i < ESNa1; i++){
    a[i]	=ESa0+ds*i;
  }
  ESSetSplA(ESa0);
  splRA(g,ESsb1a,ESsb,ESsb2a);
  ki	=ESNa1;
  lc	=rcT1a+ESNa1;
  ls	=rsT1a+ESNa1;
  drc[0]	=0.;
  drs[0]	=0.;
  for(k=1; k < ESFp1; k++){
    splRA(lc,drc+k,rcT+ki,rcT2a+ki);
    splRA(ls,drs+k,rsT+ki,rsT2a+ki);
    ki	+=ESNa1;
    lc	+=ESNa1;
    ls	+=ESNa1;
  }
  for(i=1; i < ESNa; i++){
    ESSetSplA(a[i]);
    splRA(g+i,NULL,ESsb,ESsb2a);
    ki	=0;
    lc	=rcT1a+i;
    ls	=rsT1a+i;
    for(k=0; k < ESFp1; k++){
      splRA(lc,NULL,rcT+ki,rcT2a+ki);
      splRA(ls,NULL,rsT+ki,rsT2a+ki);
      ki	+=ESNa1;
      lc	+=ESNa1;
      ls	+=ESNa1;
    }
  }
  for(i=0; i < ESNa; i++){
    ESsb[i]	=g[i];
  }
  ki	=ESNa1;
  for(k=1; k < ESFp1; k++){
    for(i=0; i < ESNa; i++){
      rcT[ki]	=rcT1a[ki];
      rsT[ki]	=rsT1a[ki];
      ki++;
    }
    ki++;
  }
  spl(ESsb,ESsb2a,ESsb1a,ESsb1a+ESNa,a,&ESNa);
  ds	=a[1]-a[0];
  s	=1./ds;
  ds	*=EZcr6;
  for(i=1; i < ESNa; i++){
    ESsb1a[i]	=(ESsb[i+1]-ESsb[i])*s-(ESsb2a[i+1]+2.*ESsb2a[i])*ds;
  }
  ki	=0;
  lc	=rcT1a;
  ls	=rsT1a;
  for(k=0; k < ESFp1; k++){
    spl(rcT+ki,rcT2a+ki,drc+k,lc+ESNa,a,&ESNa);
    spl(rsT+ki,rsT2a+ki,drs+k,ls+ESNa,a,&ESNa);
    for(i=0; i < ESNa; i++){
      *lc++	=(rcT[ki+1]-rcT[ki])*s-(rcT2a[ki+1]+2.*rcT2a[ki])*ds;
      *ls++	=(rsT[ki+1]-rsT[ki])*s-(rsT2a[ki+1]+2.*rsT2a[ki])*ds;
      ki++;
    }
    ki++;
    lc++;
    ls++;
  }
  return(0);
}

int ReSetREScs1()
{
  int i,k,ki;
  double *lc,*ls;
#ifdef XWIN
  double s,ds,g[ESNa1],a[ESNa1],drc[ESFp1],drs[ESFp1];
#else
  double s,ds,g[129],a[129],drc[33],drs[33];
#endif 
 
  ds	=(ESsa[ESNa]-ESa0)/ESNa;
  for(i=0; i < ESNa1; i++){
    a[i]	=ESa0+ds*i;
  }
  ESSetSplA(ESa0);
  splRA(g,ESsb1a,ESsb,ESsb2a);
  ki	=ESNa1;
  lc	=rcT1a+ESNa1;
  ls	=rsT1a+ESNa1;
  drc[0]	=0.;
  drs[0]	=0.;
  for(k=1; k < ESFp1; k++){
    splRA(lc,drc+k,rcT+ki,rcT2a+ki);
    splRA(ls,drs+k,rsT+ki,rsT2a+ki);
    ki	+=ESNa1;
    lc	+=ESNa1;
    ls	+=ESNa1;
  }
  for(i=1; i < ESNa; i++){
    ESSetSplA(a[i]);
    splRA(g+i,NULL,ESsb,ESsb2a);
    ki	=0;
    lc	=rcT1a+i;
    ls	=rsT1a+i;
    for(k=0; k < ESFp1; k++){
      splRA(lc,NULL,rcT+ki,rcT2a+ki);
      splRA(ls,NULL,rsT+ki,rsT2a+ki);
      ki	+=ESNa1;
      lc	+=ESNa1;
      ls	+=ESNa1;
    }
  }
  for(i=0; i < ESNa; i++){
    ESsb[i]	=g[i];
  }
  ki	=ESNa1;
  for(k=1; k < ESFp1; k++){
    for(i=0; i < ESNa; i++){
      rcT[ki]	=rcT1a[ki];
      rsT[ki]	=rsT1a[ki];
      ki++;
    }
    ki++;
  }
  spl(ESsb,ESsb2a,ESsb1a,ESsb1a+ESNa,a,&ESNa);
  ds	=a[1]-a[0];
  s	=1./ds;
  ds	*=EZcr6;
  for(i=1; i < ESNa; i++){
    ESsb1a[i]	=(ESsb[i+1]-ESsb[i])*s-(ESsb2a[i+1]+2.*ESsb2a[i])*ds;
  }
  ki	=0;
  lc	=rcT1a;
  ls	=rsT1a;
  for(k=0; k < ESFp1; k++){
    spl(rcT+ki,rcT2a+ki,drc+k,lc+ESNa,a,&ESNa);
    spl(rsT+ki,rsT2a+ki,drs+k,ls+ESNa,a,&ESNa);
    for(i=0; i < ESNa; i++){
      *lc++	=(rcT[ki+1]-rcT[ki])*s-(rcT2a[ki+1]+2.*rcT2a[ki])*ds;
      *ls++	=(rsT[ki+1]-rsT[ki])*s-(rsT2a[ki+1]+2.*rsT2a[ki])*ds;
      ki++;
    }
    ki++;
    lc++;
    ls++;
  }
  return(0);
}

int ReSetRcs2()
{
  int i,k,ki;
  double *lc,*ls;
  double f1,f2,f3;
  double EZga0=0,EZga2=0,EZga1=1e-8,EZga3=1e-8;
#ifdef XWIN
  double s,ds0,ds,g[ESNa1],a[ESNa1],drc[ESFp1],drs[ESFp1],pa[ESNa1];
#else
  double s,ds0,ds,g[129],a[129],drc[33],drs[33],pa[129];
#endif 

  ds	=(1.-ESa0)/ESNa;
  for(i=0; i < ESNa1; i++){
    a[i]	=ESa0+ds*i;
    pa[i]	=a[i]*a[i];
  }
  s	=1./ESsa[0];
  f3	=ESsb[0]*s-ESsb1a[0];
  f2	=-3.*f3*s-ESsb2a[0];
  f3	=(f3*s+EZcr2*ESsb2a[0])*s;
  f1	=3.*ESsb[0]*s-2.*ESsb1a[0]+EZcr2*ESsb2a[0]*ESsa[0];
  i	=0;
  while(a[i] < ESsa[0]){
    g[i]	=a[i]*(f1+a[i]*(f2+f3*a[i]));
    i++;
  }
  while(i < ESNa){
    ESSetSplA(a[i]);
    splRA(g+i,NULL,ESsb,ESsb2a);
    i++;
  }
  g[i]	=ESsb[ESNa];
  ds	=ESsb1a[ESNa];
  ds0	=f1+(2.*f2+3.*f3*ESa0)*ESa0;
  for(i=0; i < ESNa1; i++){
    ESsb[i]	=g[i];
  }
  EZf2spl(ESsb,ESsb2a,&ds0,&ds,g,NULL,a,&ESNa,&EZga0,&EZga1,&EZga2,&EZga3);

  s	=1./ESsa[0];

  ki	=0;
  lc	=rcT1a;
  ls	=rsT1a;
  drc[0]	=ESa0*s*lc[0];
  drs[0]	=ESa0*s*ls[0];
  *lc	=0.;
  *ls	=0.;
  f2	=(ESpa[ESNa]-ESpa[0])/(pa[ESNa]-pa[0]);
  for(i=1; i < ESNa1; i++){
    f1	=sqrt(ESpa[0]+(pa[i]-pa[0])*f2);
    ESSetSplA(f1);
    splRA(lc+i,NULL,rcT,rcT2a);
    splRA(ls+i,NULL,rsT,rsT2a);
  }

  k	=1;
  ki	+=ESNa1;
  lc	+=ESNa1;
  ls	+=ESNa1;

  f3	=rcT[ki]*s-lc[0];
  f2	=-3.*f3*s-rcT2a[ki];
  f3	=(f3*s+EZcr2*rcT2a[ki])*s;
  f1	=3.*rcT[ki]*s-2.*lc[0]+EZcr2*rcT2a[ki]*ESsa[0];
  i	=0;
  while(a[i] < ESsa[0]){
    lc[i]	=a[i]*(f1+a[i]*(f2+f3*a[i]));
    i++;
  }
  drc[k]	=f1+(2.*f2+3.*f3*ESa0)*ESa0;
  f3	=rsT[ki]*s-ls[0];
  f2	=-3.*f3*s-rsT2a[ki];
  f3	=(f3*s+EZcr2*rsT2a[ki])*s;
  f1	=3.*rsT[ki]*s-2.*ls[0]+EZcr2*rsT2a[ki]*ESsa[0];
  i	=0;
  while(a[i] < ESsa[0]){
    ls[i]	=a[i]*(f1+a[i]*(f2+f3*a[i]));
    i++;
  }
  drs[k]=f1+(2.*f2+3.*f3*ESa0)*ESa0;
  for(k=2; k < ESFp1; k++){
    ki	+=ESNa1;
    lc	+=ESNa1;
    ls	+=ESNa1;
    f3	=(lc[0]-2.*rcT[ki]*s)*s*s;	

    f2	=(3.*rcT[ki]*s-lc[0])*s;
    i	=0;
    while(a[i] < ESsa[0]){
      lc[i]	=pa[i]*(f2+f3*a[i]);
      i++;
    }
    drc[k]	=(2.*f2+3.*f3*ESa0)*ESa0;
    f3	=(ls[0]-2.*rsT[ki]*s)*s*s;	
    f2	=(3.*rsT[ki]*s-ls[0])*s;
    i	=0;
    while(a[i] < ESsa[0]){
      ls[i]	=pa[i]*(f2+f3*a[i]);
      i++;
    }
    drs[k]	=(2.*f2+3.*f3*ESa0)*ESa0;
  }
  ki	=0;
  lc	=rcT1a;
  ls	=rsT1a;
  while(i < ESNa){
    ESSetSplA(a[i]);
    ki	=0;
    lc	=rcT1a+i;
    ls	=rsT1a+i;
    for(k=0; k < ESFp1; k++){
      splRA(lc,NULL,rcT+ki,rcT2a+ki);
      splRA(ls,NULL,rsT+ki,rsT2a+ki);
      ki	+=ESNa1;
      lc	+=ESNa1;
      ls	+=ESNa1;
    }
    i++;
  }
  ki	=0;
  for(k=0; k < ESFp1; k++){
    for(i=0; i < ESNa; i++){
      rcT[ki]	=rcT1a[ki];
      rsT[ki]	=rsT1a[ki];
      ki++;
    }
    ki++;
  }
  ds	=a[1]-a[0];
  s	=1./ds;
  ds	*=EZcr6;
  ki	=0;
  lc	=rcT1a;
  ls	=rsT1a;
  for(k=0; k < ESFp1; k++){
    EZf2spl(rcT+ki,rcT2a+ki,drc+k,lc+ESNa,rcT+ki,NULL,a,&ESNa
	  ,&EZga0,&EZga1,&EZga2,&EZga3);
    EZf2spl(rsT+ki,rsT2a+ki,drs+k,ls+ESNa,rsT+ki,NULL,a
	  ,&ESNa,&EZga0,&EZga1,&EZga2,&EZga3);
    for(i=0; i < ESNa; i++){
      *lc++	=(rcT[ki+1]-rcT[ki])*s-(rcT2a[ki+1]+2.*rcT2a[ki])*ds;
      *ls++	=(rsT[ki+1]-rsT[ki])*s-(rsT2a[ki+1]+2.*rsT2a[ki])*ds;
      ki++;
    }
    ki++;
    *lc++;
    *ls++;
  }
  return(0);
}

int ESHoleEdgeBoundary1()
{
  int i,j,k,ki,kj,m,itr;
  extern double *ESjs,*ESjp,*ESaF,*ESgF;
  double ra,zt,d2gy,L0,K0,err;
#ifdef XWIN
  int indx[ES2Mp1];
  double r[ESNp1],rt[ESNp1],g[ESNp1],f[ESNp1],D[ESNp1],d[ESNp1]
    ,U[ESNp1],V[ESNp1];
  double Vr[ES2Mp1],Vi[ES2Mp1];
  double aa[ES2Mp1*ES2Mp1],h[ES2Mp1],*Ar[ESMp1],*Ai[ESMp1];
#else
  int indx[33];
  double r[129],rt[129],g[129],f[129],D[129],d[129],U[129],V[129];
  double Vr[33],Vi[33];
  double aa[33*33],h[33],*Ar[17],*Ai[17];
#endif 

  j	=0;
  Ar[0]	=aa;
  for(m=1; m < ESMp1; m++){
    j	+=ES2Mp1;
    Ar[m]	=aa+j;
    j	+=ES2Mp1;
    Ai[m]	=aa+j;
  }

  for(j=0; j < ESNp; j++){
    ra		=0.;
    r[j]	=0.;
    rt[j]	=0.;
    ki	=0;
    kj	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      r[j]	+=rcT[ki]*EScs1[kj]+rsT[ki]*ESsn1[kj];
      rt[j]	+=k*(rsT[ki]*EScs1[kj]-rcT[ki]*ESsn1[kj]);
      if(k >= ESMp){
	ra	+=rcT1a[ki]*EScs1[kj]+rsT1a[ki]*ESsn1[kj];
      }
    }
    zt		=ESsb[0]*EScs1[j];
    d[j]	=2.*ra*zt;
    r[j]	=2.*r[j]+ESaR0[0]+rcT[0];
    rt[j]	*=2.;
    g[j]	=(rt[j]*rt[j]+zt*zt)/r[j];
    f[j]	=(ESjs[0]-ESjp[0])*ESaR0[0]/r[j]+r[j]/ESaR0[0]*ESjp[0];
  }
  r[ESNp]	=r[0];
  rt[ESNp]	=rt[0];
  g[ESNp]	=g[0];
  f[ESNp]	=f[0];

  itr	=0;
  do{
    L0	=0.;
    K0	=0.;
    d2gy=0.;
    for(j=0; j < ESNp; j++){
      /* Normalization */
      ra	=0.;
      ki	=0;
      kj	=0;
      for(k=1; k < ESMp; k++){
	ki	+=ESNa1;
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	ra	+=rcT1a[ki]*EScs1[kj]+rsT1a[ki]*ESsn1[kj];
      }
      ra	=2.*ra+rcT1a[0];
      D[j]	=ra*ESsb[0]*EScs1[j]-rt[j]*(ESsb1a[0]*ESsn1[j]+rsT1a[0]);
      ra	=D[j]+d[j];
      L0	+=ra/r[j];
      K0	+=g[j]/ra;
      d2gy	+=f[j]*ra;
    }
    L0		=2.*ESgF[ESNa]*ESa0*ESNp/(ESaF[0]*L0);
    d2gy	*=-L0*L0/K0;

    ESsb1a[0]	*=L0;
    ki	=0;
    for(k=0; k < ESFp1; k++){
      rcT1a[ki]	*=L0;
      rsT1a[ki]	*=L0;
      ki	+=ESNa1;
    }
    for(j=0; j < ESNp; j++){
      D[j]	*=L0;
      d[j]	*=L0;
      ra	=D[j]+d[j];
      K0	=g[j]*d2gy/ra;
      U[j]	=K0/ra-f[j];
      V[j]	=K0+f[j]*ra;
    }

    /* Calculation of D_m */
    ESP2F(Vr,Vi,V,ES2Mp);
    for(m=1; m < ES2Mp1; m++){
      Vi[m]	=-Vi[m];
    }
    i	=0;
    h[i]	=Vr[0];
    i++;
    for(m=1; m < ESMp1; m++){
      h[i]	=Vr[m];
      i++;
      h[i]	=Vi[m];
      i++;
    }
    for(j=0; j < ESNp; j++){
      V[j]	=-rt[j]*U[j];
    }
    ESP2F(Vr,Vi,V,ES2Mp);
    j	=ES2Mp1-1;
    Ar[0][j]	=Vr[0];
    Ar[0][j]	=0.;
    for(m=1; m < ESMp1; m++){
      Ar[m][j]	=Vr[m];
      Ai[m][j]	=-Vi[m];
    }
    for(j=0; j < ESNp; j++){
      V[j]	*=ESsn1[j];
    }
    ESP2F(Vr,Vi,V,ES2Mp);
    j	=ES2Mp1-2;
    Ar[0][j]	=Vr[0];
    for(m=1; m < ESMp1; m++){
      Ar[m][j]	=Vr[m];
      Ai[m][j]	=-Vi[m];
    }
    for(j=0; j < ESNp; j++){
      V[j]	=U[j]*EScs1[j]*ESsb[0];
    }

    ESP2F(Vr,Vi,V,ES2Mp);
    for(m=1; m < ES2Mp1; m++){
      Vi[m]	=-Vi[m];
    }
    j	=0;
    Ar[0][j++]	=Vr[0];
    for(k=1; k < ESMp; k++){
      Ar[0][j++]	=2.*Vr[k];
      Ar[0][j++]	=2.*Vi[k];
    }
    for(m=1; m < ESMp1; m++){
      j	=0;
      Ar[m][j]	=Vr[m];
      Ai[m][j]	=Vi[m];
      j++;
      for(k=1; k < m; k++){
	Ar[m][j]	=Vr[m+k]+Vr[m-k];
	Ai[m][j]	=Vi[m+k]+Vi[m-k];
	j++;
	Ar[m][j]	=Vi[m+k]-Vi[m-k];
	Ai[m][j]	=Vr[m-k]-Vr[m+k];
	j++;
      }
      for(k=m; k < ESMp; k++){
	Ar[m][j]	=Vr[m+k]+Vr[k-m];
	Ai[m][j]	=Vi[m+k]-Vi[k-m];
	j++;
	Ar[m][j]	=Vi[m+k]+Vi[k-m];
	Ai[m][j]	=Vr[k-m]-Vr[m+k];
	j++;
      }
    }
    LUdcmp(aa,ES2Mp1,indx,&L0);
    LUbksb(aa,ES2Mp1,indx,h);
    i	=0;
    ki	=0;
    err	=fabs(h[i]);
    rcT1a[0]	+=h[i++];
    for(k=1; k < ESMp; k++){
      ki	+=ESNa1;
      err	+=fabs(h[i]);
      rcT1a[ki]	+=h[i++];
      err	+=fabs(h[i]);
      rsT1a[ki]	-=h[i++];
    }
    err	+=fabs(h[i]);
    ESsb1a[0]	+=h[i++];
    err	+=fabs(h[i]);
    rsT1a[0]	+=h[i];
    err		/=rcT1a[ESNa1];
    itr++;
  }while(itr < 10 && err > 1e-5);

  h[0]	=ESsb1a[0];
  h[1]	=ESsb1a[ESNa];
  splAA(ESsb,ESsb1a,ESsb2a,h,h+1);
  ki	=0;
  for(m=0; m < ESMp; m++){
    h[0]	=rcT1a[ki];
    h[1]	=rcT1a[ki+ESNa];
    splAA(rcT+ki,rcT1a+ki,rcT2a+ki,h,h+1);
    h[0]	=rsT1a[ki];
    h[1]	=rsT1a[ki+ESNa];
    splAA(rsT+ki,rsT1a+ki,rsT2a+ki,h,h+1);
    ki	+=ESNa1;
  }
  return(0);
}

int ESHoleEdgeBoundary()
{
  int i,j,k,ki,kj,m,itr;
  double ra,zt,L0,K0,D,err,rb0;
  extern double *ESjs,*ESjp,*ESaF,*ESgF;
#ifdef XWIN
  int indx[ES2Mp1];
  double r[ESNp],rt[ESNp],g[ESNp],f[ESNp],d[ESNp],U[ESNp];
  double K[ESNp],dK[ESNp],L[ESNp];
  double Kr[ESMp1],Ki[ESMp1],Lr[ESMp1],Li[ESMp1];
  double Ur[ES2Mp1],Ui[ES2Mp1];
  double aa[ES2Mp1*ES2Mp1],H[ES2Mp1],h[ES2Mp1],*Ar[ESMp1],*Ai[ESMp1];
#else
  int indx[33];
  double r[128],rt[128],g[128],f[128],d[128],U[128];
  double K[128],dK[128],L[128];
  double Kr[17],Ki[17],Lr[17],Li[17];
  double Ur[33],Ui[33];
  double aa[33*33],H[33],h[33],*Ar[17],*Ai[17];
#endif 

#ifdef H
  j	=0;
  Ar[0]	=aa;
  for(m=1; m < ESMp1; m++){
    j	+=ES2Mp1;
    Ar[m]	=aa+j;
    j	+=ES2Mp1;
    Ai[m]	=aa+j;
  }

  rb0	=2./ESsb1a[0];
  L0	=0.;
  for(j=0; j < ESNp; j++){
    D		=0.;
    ra		=0.;
    r[j]	=0.;
    rt[j]	=0.;
    ki	=0;
    kj	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      r[j]	+=rcT[ki]*EScs1[kj]+rsT[ki]*ESsn1[kj];
      rt[j]	+=k*(rsT[ki]*EScs1[kj]-rcT[ki]*ESsn1[kj]);
      err	=rcT1a[ki]*EScs1[kj]+rsT1a[ki]*ESsn1[kj];
      D		+=err;
      if(k >= ESMp){
	ra	+=err;
      }
    }
    rt[j]	*=2.;
    r[j]	=2.*r[j]+ESaR0[0]+rcT[0];
    zt		=ESsb[0]*EScs1[j];
    d[j]	=rb0*ra*zt-rt[j]*ESsn1[j];
    D		=(2.*D+rcT1a[0])*zt-rt[j]*(ESsb1a[0]*ESsn1[j]+rsT1a[0]);
    ra		=r[j]/ESaR0[0];
    r[j]	=1./r[j];
    g[j]	=(rt[j]*rt[j]+zt*zt)*r[j];
    f[j]	=(ESjs[0]-ESjp[0])*ESaR0[0]*r[j]+ESjp[0]*ra;
    L0		+=D*r[j];
  }
  L0		=2.*ESgF[ESNa]*ESa0*ESNp/(ESaF[0]*L0);

  j	=0;
  ki	=0;
  H[j++]	=rcT1a[ki]*L0;
  for(k=1; k < ESMp; k++){
    ki	+=ESNa1;
    H[j++]	=rcT1a[ki]*L0;
    H[j++]	=rsT1a[ki]*L0;
  }
  H[j++]	=ESsb1a[0]*L0;
  H[j]		=rsT1a[0]*L0;
  rb0		=1./ESsb1a[0];


  itr	=0;
  do{
    L0	=0.;
    for(j=0; j < ESNp; j++){
      L0	+=r[j]*d[j];
      U[j]	=r[j]*ESsb[0]*EScs1[j];
    }
    ESP2F(Ur,Ui,U,ESMp);
    j	=0;
    aa[j++]	=Ur[0];
    for(k=1; k < ESMp; k++){
      aa[j++]	=2.*Ur[k];
      aa[j++]	=-2.*Ui[k];
    }
    aa[j++]	=L0/ESNp;
    aa[j++]	=0.;

    for(j=0; j < ESNp; j++){
      ra	=0.;
      ki	=1;
      kj	=0;
      for(k=1; k < ESMp; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	ra	+=H[ki]*EScs1[kj]+H[ki+1]*ESsn1[kj];
	ki	+=2;
      }
      ra	=2.*ra+H[0];
      ki	=ES2Mp1-2;
      D		=ESsb[0]*ra*EScs1[j]+H[ki]*d[j]-rt[j]*H[ki+1];
      L[j]	=f[j]*D;
      K[j]	=g[j]/D;
      dK[j]	=K[j]/D;
    }
    ESP2F(Kr,Ki,K,ESMp);
    ESP2F(Lr,Li,L,ESMp);
    k	=0;
    h[k++]	=0.;
    K0	=Kr[0];
    L0	=Lr[0];
    for(m=1; m < ESMp1; m++){
      Li[m]	=-Li[m];
      Ki[m]	=-Ki[m];
      h[k++]	=L0*Kr[m]-K0*Lr[m];
      h[k++]	=L0*Ki[m]-K0*Li[m];
    }
    K0	=0.;
    for(j=0; j < ESNp; j++){
      K0	-=f[j]*d[j];
      U[j]	=-f[j]*ESsb[0]*EScs1[j];
    }
    K0	/=ESNp;
    ESP2F(Ur,Ui,U,ESMp);
    for(k=1; k < ESMp; k++){
      Ur[k]	*=2.;
      Ui[k]	*=-2.;
    }
    for(m=1; m < ESMp1; m++){
      j	=0;
      Ar[m][j]	=Kr[m]*Ur[0];
      Ai[m][j]	=Ki[m]*Ur[0];
      j++;
      for(k=1; k < ESMp; k++){
	Ar[m][j]	=Kr[m]*Ur[k];
	Ai[m][j]	=Ki[m]*Ur[k];
	j++;
	Ar[m][j]	=Kr[m]*Ui[k];
	Ai[m][j]	=Ki[m]*Ui[k];
	j++;
      }
      Ar[m][j]	=Kr[m]*K0;
      Ai[m][j]	=Ki[m]*K0;
      j++;
      Ar[m][j]	=0.;
      Ai[m][j]	=0.;
    }
    K0	=0.;
    L0	=0.;
    for(j=0; j < ESNp; j++){
      K0	-=dK[j]*d[j];
      L0	+=dK[j]*rt[j];
      U[j]	=-dK[j]*ESsb[0]*EScs1[j];
    }
    K0	/=ESNp;
    L0	/=ESNp;
    ESP2F(Ur,Ui,U,ESMp);
    for(k=1; k < ESMp; k++){
      Ur[k]	*=2.;
      Ui[k]	*=-2.;
    }
    for(m=1; m < ESMp1; m++){
      j	=0;
      Ar[m][j]		+=Lr[m]*Ur[0];
      Ai[m][j]		+=Li[m]*Ur[0];
      j++;
      for(k=1; k < ESMp; k++){
	Ar[m][j]	+=Lr[m]*Ur[k];
	Ai[m][j]	+=Li[m]*Ur[k];
	j++;
	Ar[m][j]	+=Lr[m]*Ui[k];
	Ai[m][j]	+=Li[m]*Ui[k];
	j++;
      }
      Ar[m][j]	+=Lr[m]*K0;
      Ai[m][j]	+=Li[m]*K0;
      j++;
      Ar[m][j]	+=Lr[m]*L0;
      Ai[m][j]	+=Li[m]*L0;
    }
#ifdef H
    for(m=1; m < ESMp1; m++){
      for(j=0; j < ES2Mp1; j++){
	EZout("siidd","Ar",m,j,Ar[m][j],Ai[m][j]);
      }
    }
    return(0);
#endif

    L0	=Lr[0];
    K0	=Kr[0];
    for(j=0; j < ESNp; j++){
      U[j]	=L0*dK[j]+K0*f[j];
      K[j]	=U[j]*d[j];
      L[j]	=-U[j]*rt[j];
      U[j]	*=ESsb[0]*EScs1[j];
    }
    ESP2F(Lr,Li,L,ESMp);
    ESP2F(Kr,Ki,K,ESMp);
    ESP2F(Ur,Ui,U,ES2Mp);
    for(m=1; m < ES2Mp1; m++){
      Ui[m]	=-Ui[m];
    }
    for(m=1; m < ESMp1; m++){
      j	=0;
      Ar[m][j]	+=Ur[m];
      Ai[m][j]	+=Ui[m];
      j++;
      for(k=1; k < m; k++){
	Ar[m][j]	+=Ur[m+k]+Ur[m-k];
	Ai[m][j]	+=Ui[m+k]+Ui[m-k];
	j++;
	Ar[m][j]	+=Ui[m+k]-Ui[m-k];
	Ai[m][j]	+=Ur[m-k]-Ur[m+k];
	j++;
      }
      for(k=m; k < ESMp; k++){
	Ar[m][j]	+=Ur[m+k]+Ur[k-m];
	Ai[m][j]	+=Ui[m+k]-Ui[k-m];
	j++;
	Ar[m][j]	+=Ui[m+k]+Ui[k-m];
	Ai[m][j]	+=Ur[k-m]-Ur[m+k];
	j++;
      }
      Ar[m][j]	+=Kr[m];
      Ai[m][j]	-=Ki[m];
      j++;
      Ar[m][j]	+=Lr[m];
      Ai[m][j]	-=Li[m];
    }

    LUdcmp(aa,ES2Mp1,indx,&L0);
    LUbksb(aa,ES2Mp1,indx,h);

    i	=0;
    err	=fabs(h[i]);
    H[i]	+=h[i];
    i++;
    for(k=1; k < ESMp; k++){
      err	+=fabs(h[i]);
      H[i]	+=h[i];
      i++;
      err	+=fabs(h[i]);
      H[i]	-=h[i];
      i++;
    }
    err	+=fabs(h[i]);
    H[i]	+=h[i];
    i++;
    err	+=fabs(h[i]);
    H[i]	+=h[i];

#ifdef H
    i	=0;
    err	=fabs(h[i]);
    H[i]	+=h[i++];
    for(k=1; k < ESMp; k++){
      err	+=fabs(h[i]);
      H[i]	+=h[i];
      i++;
      err	+=fabs(h[i]);
      H[i]	-=h[i];
      i++;
    }
    err	+=fabs(h[i]);
    H[i]	+=h[i];
    i++;
    err	+=fabs(h[i]);
    H[i]	+=h[i];
#endif
    err		/=H[1];
    itr++;
  }while(itr < 10 && err > 1e-5);
  h[0]	=ESsb1a[0];
  h[1]	=ESsb1a[ESNa];
  splAA(ESsb,ESsb1a,ESsb2a,h,h+1);
  ki	=0;
  for(m=0; m < ESMp; m++){
    h[0]	=rcT1a[ki];
    h[1]	=rcT1a[ki+ESNa];
    splAA(rcT+ki,rcT1a+ki,rcT2a+ki,h,h+1);
    h[0]	=rsT1a[ki];
    h[1]	=rsT1a[ki+ESNa];
    splAA(rsT+ki,rsT1a+ki,rsT2a+ki,h,h+1);
    ki	+=ESNa1;
  }
#endif
  return(0);
}

int ESHoleEdgeBoundary2()
{
  int i,j,k,ki,kj,m,itr;
  double ra,zt,L0,K0,D,err,rb0;
  extern double *ESjs,*ESjp,*ESaF,*ESgF;
#ifdef XWIN
  int indx[ES2Mp1];
  double r[ESNp],rt[ESNp],g[ESNp],f[ESNp],d[ESNp],U[ESNp];
  double K[ESNp],dK[ESNp],L[ESNp];
  double Kc[ESMp1],Ks[ESMp1],Lc[ESMp1],Ls[ESMp1];
  double Uc[ES2Mp1],Us[ES2Mp1];
  double aa[ES2Mp1*ES2Mp1],H[ES2Mp1],h[ES2Mp1],*Ac[ESMp1],*As[ESMp1];
#else
  int indx[33];
  double r[128],rt[128],g[128],f[128],d[128],U[128];
  double K[128],dK[128],L[128];
  double Kc[17],Ks[17],Lc[17],Ls[17];
  double Uc[33],Us[33];
  double aa[33*33],H[33],h[33],*Ac[17],*As[17];
#endif 

  j	=0;
  Ac[0]	=aa;
  for(m=1; m < ESMp1; m++){
    j	+=ES2Mp1;
    Ac[m]	=aa+j;
    j	+=ES2Mp1;
    As[m]	=aa+j;
  }

  rb0	=2./ESsb1a[0];
  L0	=0.;
  for(j=0; j < ESNp; j++){
    D		=0.;
    ra		=0.;
    r[j]	=0.;
    rt[j]	=0.;
    ki	=0;
    kj	=0;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      r[j]	+=rcT[ki]*EScs1[kj]+rsT[ki]*ESsn1[kj];
      rt[j]	+=k*(rsT[ki]*EScs1[kj]-rcT[ki]*ESsn1[kj]);
      err	=rcT1a[ki]*EScs1[kj]+rsT1a[ki]*ESsn1[kj];
      D		+=err;
      if(k >= ESMp){
	ra	+=err;
      }
    }
    rt[j]	*=2.;
    r[j]	=2.*r[j]+ESaR0[0]+rcT[0];
    zt		=ESsb[0]*EScs1[j];
    d[j]	=rb0*ra*zt-rt[j]*ESsn1[j];
    D		=(2.*D+rcT1a[0])*zt-rt[j]*(ESsb1a[0]*ESsn1[j]+rsT1a[0]);
    ra		=r[j]/ESaR0[0];
    r[j]	=1./r[j];
    g[j]	=(rt[j]*rt[j]+zt*zt)*r[j];
    f[j]	=(ESjs[0]-ESjp[0])*ESaR0[0]*r[j]+ESjp[0]*ra;
    L0		+=D*r[j];
  }
  L0		=2.*ESgF[ESNa]*ESa0*ESNp/(ESaF[0]*L0);
  j	=0;
  ki	=0;
  H[j++]	=rcT1a[ki]*L0;
  for(k=1; k < ESMp; k++){
    ki	+=ESNa1;
    H[j++]	=rcT1a[ki]*L0;
    H[j++]	=rsT1a[ki]*L0;
  }
  H[j++]	=ESsb1a[0]*L0;
  H[j]		=rsT1a[0]*L0;

  itr	=0;
  do{
    L0	=0.;
    for(j=0; j < ESNp; j++){
      L0	+=r[j]*d[j];
      U[j]	=r[j]*ESsb[0]*EScs1[j];
    }
    ESP2F(Uc,Us,U,ESMp);
    j	=0;
    aa[j++]	=Uc[0];
    for(k=1; k < ESMp; k++){
      aa[j++]	=2.*Uc[k];
      aa[j++]	=2.*Us[k];
    }
    aa[j++]	=L0/ESNp;
    aa[j++]	=0.;

    for(j=0; j < ESNp; j++){
      ra	=0.;
      ki	=1;
      kj	=0;
      for(k=1; k < ESMp; k++){
	kj	+=j;
	if(kj >= ESNp){
	  kj	-=ESNp;
	}
	ra	+=H[ki]*EScs1[kj]+H[ki+1]*ESsn1[kj];
	ki	+=2;
      }
      ra	=2.*ra+H[0];
      D		=ESsb[0]*ra*EScs1[j]+H[ki]*d[j]-rt[j]*H[ki+1];
      L[j]	=f[j]*D;
      K[j]	=g[j]/D;
      dK[j]	=K[j]/D;
    }
    ESP2F(Kc,Ks,K,ESMp);
    ESP2F(Lc,Ls,L,ESMp);
    k	=0;
    h[k++]	=0.;
    K0	=Kc[0];
    L0	=Lc[0];
    for(m=1; m < ESMp1; m++){
      h[k++]	=L0*Kc[m]-K0*Lc[m];
      h[k++]	=L0*Ks[m]-K0*Ls[m];
    }
    K0	=0.;
    for(j=0; j < ESNp; j++){
      K0	-=f[j]*d[j];
      U[j]	=-f[j]*ESsb[0]*EScs1[j];
    }
    K0	/=ESNp;
    ESP2F(Uc,Us,U,ESMp);
    for(k=1; k < ESMp; k++){
      Uc[k]	*=2.;
      Us[k]	*=2.;
    }
    for(m=1; m < ESMp1; m++){
      j	=0;
      Ac[m][j]	=Kc[m]*Uc[0];
      As[m][j]	=Ks[m]*Uc[0];
      j++;
      for(k=1; k < ESMp; k++){
	Ac[m][j]	=Kc[m]*Uc[k];
	As[m][j]	=Ks[m]*Uc[k];
	j++;
	Ac[m][j]	=Kc[m]*Us[k];
	As[m][j]	=Ks[m]*Us[k];
	j++;
      }
      Ac[m][j]	=Kc[m]*K0;
      As[m][j]	=Ks[m]*K0;
      j++;
      Ac[m][j]	=0.;
      As[m][j]	=0.;
    }

    K0	=0.;
    L0	=0.;
    for(j=0; j < ESNp; j++){
      K0	-=dK[j]*d[j];
      L0	+=dK[j]*rt[j];
      U[j]	=-dK[j]*ESsb[0]*EScs1[j];
    }
    K0	/=ESNp;
    L0	/=ESNp;
    ESP2F(Uc,Us,U,ESMp);
    for(k=1; k < ESMp; k++){
      Uc[k]	*=2.;
      Us[k]	*=2.;
    }
    for(m=1; m < ESMp1; m++){
      j	=0;
      Ac[m][j]	+=Lc[m]*Uc[0];
      As[m][j]	+=Ls[m]*Uc[0];
      j++;
      for(k=1; k < ESMp; k++){
	Ac[m][j]	+=Lc[m]*Uc[k];
	As[m][j]	+=Ls[m]*Uc[k];
	j++;
	Ac[m][j]	+=Lc[m]*Us[k];
	As[m][j]	+=Ls[m]*Us[k];
	j++;
      }
      Ac[m][j]	+=Lc[m]*K0;
      As[m][j]	+=Ls[m]*K0;
      j++;
      Ac[m][j]	+=Lc[m]*L0;
      As[m][j]	+=Ls[m]*L0;
    }
    L0	=Lc[0];
    K0	=Kc[0];
    for(j=0; j < ESNp; j++){
      U[j]	=L0*dK[j]+K0*f[j];
      K[j]	=U[j]*d[j];
      L[j]	=-U[j]*rt[j];
      U[j]	*=ESsb[0]*EScs1[j];
    }
    ESP2F(Lc,Ls,L,ESMp);
    ESP2F(Kc,Ks,K,ESMp);
    ESP2F(Uc,Us,U,ES2Mp);
    for(m=1; m < ESMp1; m++){
      j	=0;
      Ac[m][j]	+=Uc[m];
      As[m][j]	+=Us[m];
      j++;
      for(k=1; k < m; k++){
	Ac[m][j]	+=Uc[m+k]+Uc[m-k];
	As[m][j]	+=Us[m+k]+Us[m-k];
	j++;
	Ac[m][j]	+=Us[m+k]-Us[m-k];
	As[m][j]	+=Uc[m-k]-Uc[m+k];
	j++;
      }
      for(k=m; k < ESMp; k++){
	Ac[m][j]	+=Uc[m+k]+Uc[k-m];
	As[m][j]	+=Us[m+k]-Us[k-m];
	j++;
	Ac[m][j]	+=Us[m+k]+Us[k-m];
	As[m][j]	+=Uc[k-m]-Uc[m+k];
	j++;
      }
      Ac[m][j]	+=Kc[m];
      As[m][j]	+=Ks[m];
      j++;
      Ac[m][j]	+=Lc[m];
      As[m][j]	+=Ls[m];
    }

    LUdcmp(aa,ES2Mp1,indx,&L0);
    LUbksb(aa,ES2Mp1,indx,h);

    i	=0;
    err	=fabs(h[i]);
    EZout("siddd","rc=",i,H[i],h[i],err);
    H[i]	+=h[i];
    i++;
    for(k=1; k < ESMp; k++){
      EZout("siddddd","H=",k,H[i],h[i],H[i+1],h[i+1],err);

      err	+=fabs(h[i]);
      H[i]	+=h[i];
      i++;
      err	+=fabs(h[i]);
      H[i]	+=h[i];
      i++;
    }
    err	+=fabs(h[i]);
    H[i]	+=h[i];
    i++;
    err	+=fabs(h[i]);
    H[i]	+=h[i];
    err		/=H[1];
    itr++;
  }while(itr < 10 && err > 1e-5);

#ifdef H
  i	=0;
  ki	=0;
  h[0]	=H[i++];
  h[1]	=rcT1a[ki+ESNa];
  splAA(rcT+ki,rcT1a+ki,rcT2a+ki,h,h+1);
  for(m=1; m < ESMp; m++){
    ki		+=ESNa1;
    h[0]	=H[i++];
    h[1]	=rcT1a[ki+ESNa];
    splAA(rcT+ki,rcT1a+ki,rcT2a+ki,h,h+1);
    h[0]	=H[i++];
    h[1]	=rsT1a[ki+ESNa];
    splAA(rsT+ki,rsT1a+ki,rsT2a+ki,h,h+1);
  }
  h[0]	=H[i++];
  h[1]	=rsT1a[ESNa];
  splAA(rsT,rsT1a,rsT2a,h,h+1);
  h[0]	=H[i];
  h[1]	=ESsb1a[ESNa];
  splAA(ESsb,ESsb1a,ESsb2a,h,h+1);
#endif

#ifdef H
  ki	=0;
  i	=0;
  k	=0;
  h[k++]=H[i++]-rcT1a[ki];
  h[k++]=H[ES2Mp1-1]-rsT1a[ki];
  for(m=1; m < ESMp; m++){
    ki	+=ESNa1;
    h[k++]	=H[i++]-rcT1a[ki];
    h[k++]	=H[i++]-rsT1a[ki];
  }
  h[k]	=H[i]-ESsb1a[0];
  rb0	=H[i]*rb0-1.;
  ki	=0;
  for(m=0; m < ESFp1; m++){
    r[m]	=rcT1a[ki]*rb0;
    d[m]	=rsT1a[ki]*rb0;
    ki	+=ESNa1;
  }
  for(i=0; i < ESNa1; i++){
    ra	=(ESsa[i]-ESa0)/(ESsa[ESNa]-ESa0);
    K0	=(ESsa[i]-ESa0)*(1.-ra);
    L0	=1.-2.*ra;
    ra	=-2./(ESsa[ESNa]-ESa0);
    ki	=i;
    k	=0;
    for(m=0; m < ESMp; m++){
      rcT[ki]	+=h[k]*K0;
      rcT1a[ki]	+=h[k]*L0;
      rcT2a[ki]	+=h[k++]*ra;
      rsT[ki]	+=h[k]*K0;
      rsT1a[ki]	+=h[k]*L0;
      rsT2a[ki]	+=h[k++]*ra;
      ki	+=ESNa1;
    }
    ESsb[i]	+=h[k]*K0;
    ESsb1a[i]	+=h[k]*L0;
    ESsb2a[i]	+=h[k]*ra;
    for(; m < ESFp1; m++){
      rcT[ki]	+=r[m]*K0;
      rcT1a[ki]	+=r[m]*L0;
      rcT2a[ki]	+=r[m]*ra;
      rsT[ki]	+=d[m]*K0;
      rsT1a[ki]	+=d[m]*L0;
      rsT2a[ki]	+=d[m]*ra;
      ki	+=ESNa1;
    }
  }
#endif
  return(0);
}

int ESD2vr()
{
  int i,j,k,ki,kj,m,kk,itr;
  double ra,zt,d2gy,L0,K0,err;
  extern double *ESjs,*ESjp,*ESaF,*ESgF;

#ifdef XWIN
  int indx[ES2Mp1];
  double r[ESNp1],rt[ESNp1],g[ESNp1],f[ESNp1],D[ESNp1],U[ESNp1],V[ESNp1];
  double Vr[ES2Mp1],Vi[ES2Mp1],Ur[ES2Mp1],Ui[ES2Mp1],rtr[ES2Mp1],rti[ES2Mp1];
  double aa[ES2Mp1*ES2Mp1],h[ES2Mp1],*Ar[ESMp1],*Ai[ESMp1];
#else
  int indx[33];
  double r[129],rt[129],g[129],f[129],D[129],U[129],V[129];
  double Vr[33],Vi[33],Ur[33],Ui[33],rtr[33],rti[33];
  double aa[33*33],h[33],*Ar[17],*Ai[17];
#endif 
  
  j	=0;
  Ar[0]	=aa;
  for(i=1; i < ESMp1; i++){
    j	+=ES2Mp1;
    Ar[i]	=aa+j;
    j	+=ES2Mp1;
    Ai[i]	=aa+j;
  }
  ki	=0;
  rtr[0]	=0.;
  rti[0]	=0.;
  for(k=1; k < ES2Mp1; k++){
    ki	+=ESNa1;
    rtr[k]	=k*rcT[ki];
    rti[k]	=-k*rsT[ki];
  }

  for(j=0; j < ESNp; j++){
    ra	=0.;
    rt[j]	=0.;
    ki	=0;
    kj	=0;
    for(k=1; k < ESFp1; k++){
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      ki	+=ESNa1;
      ra	+=rcT1a[ki]*EScs1[kj]+rsT1a[ki]*ESsn1[kj];
      rt[j]	+=k*(rsT[ki]*EScs1[kj]-rcT[ki]*ESsn1[kj]);
    }
    ra	=2.*ra+rcT1a[0];
    rt[j]	*=2.;
    D[j]	=ra*ESsb[0]*EScs1[j]-rt[j]*(ESsb1a[0]*ESsn1[j]+rsT1a[0]);
  }
  ESP2F(Ur,Ui,D,ESMp);
  for(k=1; k < ESMp1; k++){
    Ui[k]	=-Ui[k];
  }  
  itr	=0;
  do{
    /* Conversion of D into dr */
    L0	=EZcr2*ESsb[0];
    memset((void*)aa,'\0',(size_t)(sizeof(double)*ES2Mp1*ES2Mp1));
    Ar[0][1]	=ESsb[0];
    Ar[0][ES2Mp1-2]=rtr[1];

    Vr[0]	=ESsb[0]*rcT1a[ESNa1]+rtr[1]*ESsb1a[0];
    Vi[0]	=0.;
    for(k=1; k < ESMp1; k++){
      Vr[k]	=0.;
      Vi[k]	=0.;
      kk	=k-1;
      ki	=ESNa1*kk;
      if(kk == 0){
	Ar[k][0]	=L0;
	Vr[k]	+=L0*rcT1a[ki];
      }
      else{
	j	=2*kk-1;
	Ar[k][j]	=L0;
	Vr[k]	+=L0*rcT1a[ki];
	j++;
	Ai[k][j]	=L0;
	Vi[k]	-=L0*rsT1a[ki];
      }
      kk	=k+1;
      ki	=ESNa1*kk;
      if(kk < ESMp){
	j	=2*kk-1;
	Ar[k][j]	=L0;
	Vr[k]	+=L0*rcT1a[ki];
	j++;
	Ai[k][j]	=L0;
	Vi[k]	-=L0*rsT1a[ki];
      }
      else{
	Ur[k]	-=L0*rcT1a[ki];
	Ui[k]	+=L0*rsT1a[ki];
      }
      j		=ES2Mp1-2;
      kk	=k-1;
      Ar[k][j]	=-EZcr2*rtr[kk];
      Ai[k][j]	=-EZcr2*rti[kk];
      Vr[k]	-=EZcr2*rtr[kk]*ESsb1a[0];
      Vi[k]	-=EZcr2*rti[kk]*ESsb1a[0];
      kk	=k+1;
      if(kk < ESMp){
	Ar[k][j]	+=EZcr2*rtr[kk];
	Ai[k][j]	+=EZcr2*rti[kk];
	Vr[k]	+=EZcr2*rtr[kk]*ESsb1a[0];
	Vi[k]	+=EZcr2*rti[kk]*ESsb1a[0];
      }
      else{
	Ur[k]	-=EZcr2*rtr[kk]*ESsb1a[0];
	Ui[k]	-=EZcr2*rti[kk]*ESsb1a[0];
      }
      j++;
      if(k < ESMp){
	Ar[k][j]=rti[k];
	Ai[k][j]=-rtr[k];
	Vr[k]	+=rti[k]*rsT1a[0];
	Vi[k]	-=rtr[k]*rsT1a[0];
      }
      else{
	Ur[k]	-=rti[k]*rsT1a[0];
	Ui[k]	+=rtr[k]*rsT1a[0];
      }
    }
    i	=0;
    h[i++]	=Ur[0];
    for(k=1; k < ESMp1; k++){
      h[i++]	=Ur[k];
      h[i++]	=Ui[k];
    }
    LUdcmp(aa,ES2Mp1,indx,&L0);
    LUbksb(aa,ES2Mp1,indx,h);

    i	=0;
    ki	=0;
    err	=fabs(rcT1a[0]-h[i]);
    i++;
    for(k=1; k < ESMp; k++){
      ki	+=ESNa1;
      err	+=fabs(rcT1a[ki]-h[i]);
      i++;
      err	+=fabs(rsT1a[ki]+h[i]);
      i++;
    }
    err	+=fabs(ESsb1a[0]-h[i]);
    i++;
    err	+=fabs(rsT1a[0]-h[i]);
    err		/=fabs(rcT1a[ESNa1]*ESsb[0]+rcT[ESNa1]*ESsb1a[0]);
    itr++;
  }while(itr < 1 && err > 1e-3);
  return(0);
}

int ReSetHoleSize()
{
  if(ESsa[0] < ESa0){
    ReSetREScs1();
  }
  else{
    ReSetRcs2();
  }
  return(0);
}

int ESInputPlVac()
{
  static int Fbn=0;
  static int FlCntrP=1,iT=0;
  static int kT=0,kT1=0,kT2=0,kTE=0;
  int k2jet[256];

  static double ttime=0.,time1=0.,time2=0.,av=1.;

  int nD;
  int i,j,Fl,Fl0,Fl1,*k2,VMnAP,VMnAP1;
  double *prD,*pzD,*pr,*pz,*pvr,*pvz;
  extern void IbPFIPlot();

  Fl	=0;

#ifdef SHMEM
  if(ESShMem) ESReadPlVacShM();
#endif/*SHMEM*/

#ifdef PFC
  if(FlCntrP){
    ESGetPlVac(rPlVdT,zPlVdT,k2kT);
  }
#endif

  Fbn	=FlagPoints;
#ifndef Tbl_2DGeom
  if(VMFlag == 30){
#ifndef AscIPlV
#undef AscIPlV
#endif/*AscIPlV*/
    Fbn	=FlagPoints;
    VMFlag	=0;
  }
  FlagPoints	=Fbn;
  prD	=rPlVdT;
  pzD	=zPlVdT;
  k2	=k2kT;
  nD	=12;
  kPlVd	=ESComprPlVacData(prD,pzD,k2,nPlVd);
  ESCheckIntData();

  pr	=ESsr;
  pz	=ESsz;
  VMnAP	=ESnAP;
  VMnAP1=ESnAP;
  VMnAP	=0;
  pvr	=ECr;
  pvz	=ECz;

  if(NPlVf < MPlVf){
    ESReInit1DFPlV();
    Fl	=1;
  }
  nDFT	=MPlVf*ESNt1;
  if(NDFT < nDFT){
    ESReInit2DFTPlV();
    Fl	=1;
  }
  if(nPlVf != MPlVf) Fl=1;
  nPlVcs	=ESNp1*MPlVf;
  if(NPlVcs < nPlVcs){
    ESReInit2DFPPlV();
    Fl	=1;
  }
  neq=2*MPlVf-1+kPlVd;
  if(Neq < neq) ESReInitPlVeq(neq,&ac,&f,&indx);
  if(Fl || nPlVf != MPlVf){
    ESInit2DFPPlV();
    Fl	=0;
  }
#ifdef PFC
  CbUserPlot	=IbPFIPlot;
  if(av >= 1.){
    av	=1.;
  }
  if(av < 0.){
    av	=0.;
  }
  if(ESav != av){
    ESa0	=0.;
    ESav	=av;
    ESFlVirtualB=63;
    ESiIv	=0;
    while(ESiIv < ESNa && ESsa[ESiIv+1] < ESav){
      ESiIv++;
    }
  }
  if((ESFlVirtualB&1)){
    ReSetVirtualBSize();
    ESFlVirtualB&=~1;
  }
  if((ESFlVirtualB&2)){
    ESInit1ACore();
    ESFirstSetSplA();
    ESFlVirtualB&=~2;
  }

  if(ESa0 != ESsa[0]){
    ESFlHole	=63;
  }
  if((ESFlHole&1)){
    ReSetHoleSize();
    ESFlHole	&=~1;
  }	
  if((ESFlHole&2)){
    ESInit1ACore();
    ESFirstSetSplA();
    ESFlHole	&=~2;
  }
#endif

#ifdef TSC
  if(VMFlag == 40 || kT1 || kT2 || kT){
    ECSetRestoreESCTime(&time1,&time2,&ttime,kT1,kT2,kT);
    kT1	=0;
    kT2	=0;
  }
  if(kTE) ESTime=ttime;
  if(kT){
    ECRestoreESC(ttime);
    VMFlag	=0;
  }
  switch(VMFlag){
  case 0:
    ECEqIterFl=1; /* Newton iterations */
    break;
  case 1:
    ECEqIterFl=0; /* Simple iterations */
    break;
#ifdef PFC
  case 2:
    ESGetPlVac(rPlVdT,zPlVdT,k2kT);
    break;
  case 3:
    ESSetCntrPoints(rPlVdT,zPlVdT,k2kT);
    FlagPoints	=0;
    break;
#endif
  case 4:
    ESkNormJ=1;
    break;
  case 5:
    ESkNormJ=0;
    break;
  case 8:
    DCON2esc();
    i	=ESnAP-ESNp1;
    VMnAP	=ESNp1;
    pvr	=EZvr+i;
    pvz	=EZvz+i;
    VMFlag	=0;
    break;
  case 9:
    if(ECSoloviev(prD,pzD,kPlVd) == 0){
      VMnAP	=ESnAP;
      pvr	=EZvr;
      pvz	=EZvz;
    }
    break;
  case 13:
    {
      int j;
      FILE *lf;

      kW	=MPlVf;
      j	=ESNa;
      rcW[0]	=ESaR0[0]+rcT[j];
      rsW[0]	=ESaZ0[0]+rsT[j];
      bW	=ESsb[j];
      for(i=1; i < kW; i++){
	j	+=ESNa1;
	rcW[i]	=rcT[j];
	rsW[i]	=rsT[j];
      }	
      lf	=fopen("Boundary","w");
      j	=0;
      fprintf(lf,"%2d - number of PlV control points\n",kPlVd);
      for(i=0; i < nPlVd; i++){
	if(k2[i]){
	  fprintf(lf,"%2d %2d %13.6e %13.6e\n",j,i,prD[i],pzD[i]);
	  j++;
	}
      }
      fprintf(lf,"%2d %13.6e - number of Fourier coefficients and bW\n",kW,bW);
      for(i=0; i < kW; i++){
	fprintf(lf,"%2d %13.6e %13.6e\n",i,rcW[i],rsW[i]);
      }
      printf("Boundary has been recorded\n");
      fclose(lf);
    }
  case 16:
    {
      int nr1=10,nz1=10;
      double gY[21*21];
      double Br[2]={1.8,4.},Bz[2]={-1.87,1.87};
      EcNjet	=0;
      escpsi_(Br,Br+1,Bz,Bz+1);
#ifdef H
      FluxAtLabRZ(gY,Br,Bz,nr1,nz1);
#endif
      VMnAP	=EcNjet;
      pvr	=EcRjet;
      pvz	=EcZjet;
    }
    break;
  case 20:
    EcNjet	=ESNp1;
    {
      double EZz0,r0,a,b,gd;
      r0=EZcr2*(prD[3]+prD[2]);
      a	=EZcr2*(prD[3]-prD[2]);
      EZz0=EZcr2*(pzD[0]+pzD[1]);
      b	=EZcr2*(pzD[0]-pzD[1]);
      gd	=(prD[0]-r0)/a;
      gd	=-asin(gd);
      prD[1]	=prD[0];
      for(i=0; i < EcNjet; i++){
	EcRjet[i]	=r0+a*cos(ESgt[i]+gd*ESsn1[i]); 
	EcZjet[i]	=EZz0+b*ESsn1[i];
      }
    }
    FlagPoints	=1;
    VMnAP	=EcNjet;
    pvr	=EcRjet;
    pvz	=EcZjet;
    break;
  case 40:
    i=ECRestoreESC(ttime);
    VMFlag	=0;
    break;
  case 41:
    ESFromHardDrive(ttime);
    VMFlag	=0;
    break;
  }
#endif
  if(ESiAx){
    ESReal2RefDataX(iT);
    ESPlVdata2PlVX();
    ESGapCoreSepX();
    ESCoreX();
    VMnAP	=ESNp1*ESnBL1;
    pvr	=rBL;
    pvz	=zBL;
  }
  else{
    ESReal2RefData(iT);
    ESReal2RefGeom(iT);
    switch(Fbn){
    case 1:
      ESPlVPoints2PlV(EcRjet,EcZjet,EcNjet);
      VMnAP	=EcNjet;
      pvr	=EcRjet;
      pvz	=EcZjet;
      break;
    case 2:
      ESPlVSamePlV();
      break;
    case 0:
    default:
      ESPlVdata2PlV();
      break;
    }
    if(ESa0 != 0.){
      ESHoleEdgeBoundary();
    }
    ESRef2RealGeom(iT);
  }
  ESRealGeom(iT);
  ESRBt	=ESBt*ESRext;
#endif/*Tbl_2DGeom*/
#ifdef TSCW
  ESiAx	=0;
#endif
#ifdef PFC
  if(FlCntrP){
    ESSetCntrPoints(rPlVdT,zPlVdT,k2kT);
    PFSetPlotPoints();
    FlCntrP	=0;
  }
#endif
#ifndef AscOPlV
#undef AscOPlV
#endif/*AscOPlV*/
  return(0);
}

int ESInput2DEqSolv()
{
  static int Fl=0,Fc=0,Fd=0,k0=0,k1=11,i0=0,i1=20,j0=0,j1=64,m=0,n=4;
  int i,j,Fl0,Fl1;
  double *ga,g[4],x[2]={0.,1.,},Y[2]={0.,1.,},y[2];
  double xc[512],yc[512],ys[512];
  int Nc,Ns;

#ifdef DEBUG
  double *Wx,*W1[64];
  j	=(ESNa1*8+1)*sizeof(double);
  Wx	=(double*)malloc(j);

  for(i=0; i < 33; i++){
    W1[i]	=(double*)malloc(j);
    if(W1[i] == NULL){
      printf("W1[%d] in ESInput2DEqSolv() memory falure\n",i);
      ESexit(0);
    }
  }  
#endif

  ESGetSmoothAddr(&ga);
  for(i=0; i < 4; i++) g[i]=ga[i];

  Fl0	=ESEqSolvFl&0x0F;
  if(Fl0 > 1) Fl0=1;
  Fl1	=(ESEqSolvFl&0xF0)>>4;
#ifndef Tbl_EqSolv
  if(ESEqSolvIt < 0){
    ESEqSolvIt	=1;
  }
  if(Fl1 < 0){
    Fl1	=0;
  }
  if(Fl1 > 0xF0){
    Fl1	=0xF0;
  }
  ESEqSolvFl	=(Fl1<<4)+Fl0;

  if(ESMp < 2) ESMp=2;
  if(ESMp*4 > ESNp) ESMp=ESNp/4;
  if(ESMp > ESFp) ESMp=ESFp;
  ESMp1	=ESMp+1;
  ES2Mp	=2*ESMp;
  ES2Mp1=ES2Mp+1;
  ESnMp	=ESNa1*ES2Mp1;

  if(k0 < 0) k0	=0;
  if(k0 > ES2Mp) k0	=ES2Mp;
  if(k1 < k0+1) k1	=k0+1;
  if(k1 > ES2Mp1) k1	=ES2Mp1;
  if(i0 < 0) i0	=0;
  if(i0 > ESNa-1) i0	=ESNa-1;
  if(i1 < i0+1) i1	=i0+1;
  if(i1 > ESNa) i1	=ESNa;
  if(j1 > ESNp) j1	=ESNp;
  j	=0;
  if(g[2] != ga[2]){
    j	=1;
  }
  for(i=0; i < 4; i++){
    g[i]	=ga[i];
  }
#ifdef DEBUG
  if(j){
    int k,ki;
    double *p,*p1,*p2,s,t;
    ECSmooth2DCoord(0);
    p		=rcT;
    p1	=rcT1a;
    p2	=rcT2a;
    y[0]	=p[0];
    y[1]	=y[0];
    x[0]	=0.;
    x[1]	=1.;
    ki		=0;
    s		=0.125*ESsa[1];
    for(k=0; k < ESFp1; k++){
      ki	=ESNa1*k;
      kW	=0;
      t		=0.;
      while(t <= ESsa[1]){
	ESSetSplA(t);
	Wx[kW]	=t;
	splRA(W1[k]+kW,NULL,p+ki,p2+ki);
	if(y[0] > W1[k][kW]){
	  y[0]	=W1[k][kW];
	}
	if(y[1] < W1[k][kW]){
	  y[1]	=W1[k][kW];
	}
	t	+=s;
	kW++;
      }
    }
    Scale2D(3,2,x,y,2,6,2);
    for(k=0; k < ESFp1; k++){
      ki	=ESNa1*k;
      Plot2d(3,ESsa,p+ki,ESNa1,6,0,14,0);
      Plot2d(3,Wx,W1[k],kW,6,0,14,0);
    }
    CbFlush();
  }
#endif
  ESReInitMetrTnsrCS();
  ECRcs3D2Rcs2D(0);
  if(ESa0 == 0.){
    ES3DMetrTnsrCS();
  }
  else{
    if(ESav == 0.){
      ESHoleEdgeBoundary();
      ES3DMetrTnsrCSHole();
    }
    else{
      ES3DMetrTnsrCSHoleV();
    }
  }
  Nc	=0;
  Ns	=0;
#ifdef DEBUG
  ESInspectMetricTensor(xc,yc,ys,m,i0,i1,j0,j1,k0,k1,Fl,Fd,n);
  Nc	=n*(i1-i0)+1;
  Y[0]	=0.;
  Y[1]	=0.;
  x[0]	=xc[0];
  x[1]	=xc[Nc-1];
  Ns	=m ? Nc : 0;
  for(i=0; i < Nc; i++){
    if(Y[0] > yc[i]) Y[0]	=yc[i];
    if(Y[1] < yc[i]) Y[1]	=yc[i];
    if(m){
      if(Y[0] > ys[i]) Y[0]	=ys[i];
      if(Y[1] < ys[i]) Y[1]	=ys[i];
    }
  }
#endif
#endif/*Tbl_EqSolv*/

#ifdef DEBUG
  for(i=0; i < 33; i++){
    free(W1[i]);
  }  
  free(Wx);
#endif
  ESMp1	=ESMp+1;
  ES2Mp	=2*ESMp;
  ES2Mp1=ES2Mp+1;
  ESnMp	=ESNa1*ES2Mp1;
  return(0);
}

int ESInitPlGeom()
{
  static int i,iT=0,iT1;
  int Fl,kP,*k2;
  double *prD,*pzD,*pr,*pz;

  iT1=iT;
  ESGetPncrAddr(&nPncr,&nPsrf,&rPncr,&zPncr,&xPncr,&yPncr);
  KPncr	=nPncr/nPsrf;
  LPncr	=ESNLt;
  Fl	=0;
#ifndef Tbl_Geom
  if(iT < 0){
    iT=0;
  }
  if(iT > ESNt){
    iT=ESNt;
  }
  if(ESMp < 3){
    ESMp=3;
  }
  if(ESMp*4 > ESNp){
    ESMp=ESNp/4;
  }
  if(ESMp > ESFp1){
    ESMp=ESFp1;
  }
  ESMp1	=ESMp+1;
  ES2Mp	=2*ESMp;
  ES2Mp1=ES2Mp+1;
  ESnMp	=ESNa1*ES2Mp1;
  if(iT1 != iT){
    ESRef2RealData(iT1);
    ESPlVdata2PlV();
    ESRef2RealGeom(iT1);
    ESRealGeom(iT1);
    iT1	=iT;
  }
  prD	=rPlVdT+nPlVd*iT;
  pzD	=zPlVdT+nPlVd*iT;
  k2	=k2kT+nPlVd*iT;
  pr	=ESsr+ESnAP*iT;
  pz	=ESsz+ESnAP*iT;
  kPlVd	=ESComprPlVacData(prD,pzD,k2,nPlVd);
  ESCheckIntData();

  if(NPlVf < MPlVf){
    ESReInit1DFPlV();
    Fl=1;
  }
  nDFT	=MPlVf*ESNt1;
  if(NDFT < nDFT){
    ESReInit2DFTPlV();
    Fl=1;
  }
  if(nPlVf != MPlVf){ 
    Fl=1;
  }
  nPlVcs	=ESNp1*MPlVf;
  if(NPlVcs < nPlVcs){
    ESReInit2DFPPlV();
    Fl=1;
  }
  neq=2*MPlVf-1+kPlVd;
  if(Neq < neq){
    ESReInitPlVeq(neq,&ac,&f,&indx);
  }
  nrC	=MPlVf*(nPsrf+1);
  if(NrC < nrC){
    ESReInit2DAPrC();
  }
  if(Fl){
    ESInit2DFPPlV();
    ESInit3DCoord();
    Fl=0;
  }
  ESReal2RefData(iT);
  ESReal2RefGeom(iT);
  R0	=ESaR0[iT];
  Z0	=ESaZ0[iT];
  ESPlVdata2PlV();
  ESRef2RealGeom(iT);
  ESRealGeom(iT);
  if(nPncr){
    kP	=0;
    for(i=iT+KPncr*(nPsrf-1); i < nPncr; i +=LPncr){
      xPncr[kP]	=rPncr[i];
      yPncr[kP]	=zPncr[i];
      kP++;
    }
  }
#endif/*Tbl_Geom*/
#ifdef H
  rPlV	=r+ESNp1*ESNa;
  zPlV	=z+ESNp1*ESNa;
  for(i=0; i < ESNt1; i++){
    ESReal2RefData(i);
    ESReal2RefGeom(i);
    ESPlVdata2PlV();
    ESRef2RealGeom(i);
    ESRealGeom(i);
  }
#endif
  return(0);
}

int ESA2RefRZ()
{
  static int iT=0,iA=20,nM=8,nN=50;
  int kP;
  double *pr,*pz,Rmn,Rmx,Zmn,Zmx;
  
  ESGetgFPlV();
#ifndef Tbl_RefRTZ
  ESgFPlV1a=2.*ESgFPlV;
  if(iT < 0){
    iT=0;
  }
  if(iT > ESNt){
    iT=ESNt;
  }
  pr	=ESsr+ESnAP*iT;
  pz	=ESsz+ESnAP*iT;
  if(iA < 0){
    iA=0;
  }
  if(iA > ESNa){
    iA=ESNa;
  }
  if(nM < 1){
    nM=1;
  }
  if(nM > ESNp/2){
    nM=ESNp/2;
  }
  if(nN < 0){
    nN=0;
  }
  if(nN > ESNLt/2){
    nN=ESNLt/2;
  }

#ifdef H
  ESSmoothRZ(nM,nN);
#endif

  ES3DCoord2Spl();
  ES3DMetrTnsrCS();
  ES3DInvG22();

  ESAextBext2Spl();
  ESNormDispl(iA,iT,nM,nN);
  if(nPncr){
    int i,is,ii;
    kP	=0;
    Rmn=rPncr[iT];
    Rmx=rPncr[iT];
    Zmn=zPncr[iT];
    Zmx=zPncr[iT];
    for(is=0; is < nPsrf; is++){
      ii	=iT+KPncr*is;
      for(i=iT; i < KPncr; i +=LPncr){
	xPncr[kP]	=rPncr[ii];
	yPncr[kP]	=zPncr[ii];
	if(Rmn > xPncr[kP])
	  Rmn=xPncr[kP];
	if(Rmx < xPncr[kP])
	  Rmx=xPncr[kP];
	if(Zmn > yPncr[kP])
	  Zmn=yPncr[kP];
	if(Zmx < yPncr[kP])
	  Zmx=yPncr[kP];
	kP++;
	ii	+=LPncr;
      }
    }
  }
#endif/*Tbl_RefRTZ*/
  return(0);
}

#endif

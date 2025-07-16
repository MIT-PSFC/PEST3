#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern int ESNa,ESNp,ESNa1,ESNp1,ESnAP;
extern int ESFp,ESFp1,ESnAF;
extern int ESFail;
extern int FlagPoints;

extern double *ESgF0,*ESgF01a
,*ESgY0,*ESgY01a,*ESgY02a
,*ESgi,*ESgi1a,*ESgi2a
,*ESaY,*ESaY1a,*ESaY2a;
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

extern int ESEqSolvRadC,ESEqSolvInPr,ESEqSolvInCr;
extern int ESnDtPr,ESnDtCr;
extern double ESxDtPr[],ESxDtCr[];
extern double ESBt,ESRext,ESRBt;

extern double *ECr,*ECz,*EZdra,*EZdza;
extern double *ESgt;
extern double R0,Z0;
extern double *ESaR0,*ESaZ0,*ESsb,*ESsb1a,*ESsb2a;
extern double *EScs1,*ESsn1;
extern double EcRjet[],EcZjet[];
extern int EcNjet;
extern double *ECb,*ECb1a,*ECb2a;

extern double *rcT,*rcT1a,*rcT2a,*rsT,*rsT1a,*rsT2a;
extern double *drcT,*drcT1a,*drcT2a,*drsT,*drsT1a,*drsT2a;

extern int nf_X,ind_sep,NPlVac,*k2kPlV,ESiAx;
extern double *rPlVd,*zPlVd;
extern double b_X,r_X,z_X0,z_X2,z_X2n,t_X;
extern double *EZx0sep,*EZx1sep,*EZx2sep,*X2vb,*EZx0cs,*EZx0sn,*EZdx0cs,*EZdx0sn,
*X2vbc,*X2vbs;
extern double *EZrvb,*EZzvb,*EZrvbcs,*EZrvbsn,*drvbcs,*drvbsn;

extern double EZd_vb,b_X,z_X0,z_X2;
extern double EZdx_0x,d2x_0x,EZdx_0o,d2x_0o;

extern double d_vbE,b_XE,z_X0E,z_X2E;
extern double dx_0xE,d2x_0xE,dx_0oE,d2x_0oE;
extern double *x0sepE,*x1sepE,*x2sepE;

extern double *Dx0,*Dx2,*Ddx2,*Dz2,*Dx2z,*Ddx2z,*Dx0c,*Dx0s,*Dx2c,*Dx2s,
*Ddx2c,*Ddx2s,*Dz2c,*Dz2s,*Dx2zc,*Dx2zs,*Ddx2zc,*Ddx2zs;

extern double *EZrcs,*EZrsn,*EZd1rcs,*EZd1rsn,*EZd2rcs,*EZd2rsn,*EZz0,*EZd1z0,*EZd2z0;
extern double *rcE,*rsE,*rcE1,*rsE1,*rcE2,*rsE2,*bE,*bE1,*bE2,R0E,Z0E;
extern double*Fcs,*Fsn,*dFcs,*dFsn,*EZxinf,*A2m,*EZyinf,*d1yinf,*d2yinf,*EZxgt,*EZygt;
extern double Eps_tr;

extern int ESRecoveryFlag;
extern char ESRunID[];
extern double ESTime; 

int EcInitArrays()
{
  int i,j,k;
  double t;

  t_X	=2.*EZcgp_4;
  EZx0sep	=(double*)malloc((22*ESFp1+15*ESNp1)*sizeof(double));
  if(EZx0sep == NULL){
    printf("Failure in memory allocation for EZx0sep in EcInitArrays\n");
    exit(0);
  }
  EZx1sep	=EZx0sep	+ESNp1;
  EZx2sep	=EZx1sep	+ESNp1;
  x0sepE=EZx2sep	+ESNp1;
  x1sepE=x0sepE	+ESNp1;
  x2sepE=x1sepE	+ESNp1;
  X2vb	=x2sepE	+ESNp1;
  EZrvb	=X2vb	+ESNp1;
  EZzvb	=EZrvb	+ESNp1;
  Dx0	=EZzvb	+ESNp1;
  Dx2	=Dx0	+ESNp1;
  Ddx2	=Dx2	+ESNp1;
  Dz2	=Ddx2	+ESNp1;
  Dx2z	=Dz2	+ESNp1;
  Ddx2z	=Dx2z	+ESNp1;

  EZx0cs	=Ddx2z	+ESNp1;
  EZx0sn	=EZx0cs	+ESFp1;
  EZdx0cs	=EZx0sn	+ESFp1;
  EZdx0sn	 =EZdx0cs	+ESFp1;
  EZrvbcs	=EZdx0sn	+ESFp1;
  EZrvbsn	=EZrvbcs	+ESFp1;
  drvbcs=EZrvbsn	+ESFp1;
  drvbsn=drvbcs	+ESFp1;
  Dx0c	=drvbsn	+ESFp1;
  Dx0s	=Dx0c	+ESFp1;
  Dx2c	=Dx0s	+ESFp1;
  Dx2s	=Dx2c	+ESFp1;
  Ddx2c	=Dx2s	+ESFp1;
  Ddx2s	=Ddx2c	+ESFp1;
  Dz2c	=Ddx2s	+ESFp1;
  Dz2s	=Dz2c	+ESFp1;
  Dx2zc	=Dz2s	+ESFp1;
  Dx2zs	=Dx2zc	+ESFp1;
  Ddx2zc=Dx2zs	+ESFp1;
  Ddx2zs=Ddx2zc	+ESFp1;
  X2vbc	=Ddx2zs	+ESFp1;
  X2vbs	=X2vbc	+ESFp1;
  for(j=0; j < ESNp1; j++){
    EZx0sep[j]=0.;
    EZx1sep[j]=0.;
    EZx2sep[j]=0.;
  }
  EZrcs	=(double*)malloc((12*ESnAF+6*ESNa1)*sizeof(double));
  if(EZrcs == NULL){
    printf("Failure in memory allocation for EZrcs in EcInitArrays\n");
    exit(0);
  }
  EZrsn	=EZrcs	+ESnAF;
  EZd1rcs	=EZrsn	+ESnAF;
  EZd2rcs	=EZd1rcs	+ESnAF;
  EZd1rsn	=EZd2rcs	+ESnAF;
  EZd2rsn	=EZd1rsn	+ESnAF;
  EZz0	=EZd2rsn	+ESnAF;
  EZd1z0	=EZz0	+ESNa1;
  EZd2z0	=EZd1z0	+ESNa1;
  rcE	=EZd2z0	+ESNa1;
  rsE	=rcE	+ESnAF;
  rcE1	=rsE	+ESnAF;
  rcE2	=rcE1	+ESnAF;
  rsE1	=rcE2	+ESnAF;
  rsE2	=rsE1	+ESnAF;
  bE	=rsE2	+ESnAF;
  bE1	=bE	+ESNa1;
  bE2	=bE1	+ESNa1;
  Fcs	=(double*)malloc((4*ESFp1+18*ESNa1+2*ESNp1)*sizeof(double));
  if(Fcs==NULL){
    printf("Failure in memory allocation for Fcs in EcInitArrays\n");
    exit(0);
  }
  Fsn	=Fcs	+ESFp1;
  dFcs	=Fsn	+ESFp1;
  dFsn	=dFcs	+ESFp1;
  A2m	=dFsn	+ESFp1;
  EZxinf	=A2m	+ESNa1;
  EZyinf	=EZxinf	+4*ESNa1;
  d1yinf=EZyinf	+4*ESNa1;
  d2yinf=d1yinf	+4*ESNa1;
  EZxgt	=d2yinf	+4*ESNa1;
  EZygt	=EZxgt	+ESNp1;

  NPlVac= 18;
  rPlVd= (double*)malloc(2*NPlVac*sizeof(double));
  if(rPlVd==NULL){
    printf("Allocation failure for rPlVd in InitPVac()\n");
    exit(0);
  }
  zPlVd= rPlVd+NPlVac;
  k2kPlV	=(int*)malloc(NPlVac*sizeof(int));
  return(0);
}

int EcDeInitArrays()
{
  free(k2kPlV);
  free(rPlVd);
  free(Fcs);
  free(EZrcs);
  free(EZx0sep);
  return(0);
}

int ESDeInitEc()
{
  ECDeInitLSODE();
  ECDeInitRKSolv();
  ECDeInitFgY();
  ECDeVirtualPlV2vdrPlV();
  EcDeInitArrays();
  return(0);
}

int ECRefGeom(int n)
{
  int i,j,ji,in,k,ki,kj;
  double *pr,*pz,b;
  double *rc,*rs,*prc,*prs;

  pr	=ECr;
  pz	=ECz;
  rc	=EZdra;
  rs	=rc+ESFp1;
  in	=ESNa1*n;
  ki	=ESnAF*n;
  prc	=rcT+ki;
  prs	=rsT+ki;

  ji	=0;
  for(j=0; j < ESNp1; j++){
    pr[ji]	=ESaR0[n];
    pz[ji]	=ESaZ0[n];
    ji++;
  }
  in++;
  for(i=1; i < ESNa1; i++){
    b	=ESsb[in];
    k	=0;
    ki	=i;
    rc[0]=ESaR0[n]+prc[ki];
    rs[0]=ESaZ0[n]+prs[ki];
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
	if(kj >= ESNp)
	  kj	-=ESNp;
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
  return(0);
}

int ECvr2r(int n)
{
  int i,j,k,ji,ki,kj;
  double rc,rs;

  ki	=ESnAF*n;
  ji	=0;
  for(i=0; i < ESNa1;i++){
    rc		=ESaR0[n]+rcT[ki];
    for(j=0; j < ESNp1; j++){
      ECr[ji]	=rc;
      ji++;
    }
    ki++;
  }
  for(k=1; k < ESFp1; k++){
    ji		=ESNp1;
    ki++;
    for(i=1; i < ESNa1; i++){
      rc	=2.*rcT[ki];
      rs	=2.*rsT[ki];
      kj	=0;
      for(j=0; j < ESNp; j++){
	ECr[ji]	+=rc*EScs1[kj]+rs*ESsn1[kj];
	kj	+=k;
	if(kj >= ESNp)
	  kj	-= ESNp;
	ji++;
      }
      ECr[ji]	=ECr[ji-ESNp];
      ji++;
      ki++;
    }
  }
  return(0);
}

int ECvz2z(int n)
{
  int i,j,ji,ki,in;
  double y0,b;

  in	=ESNa1*n;
  ki	=ESnAF*n;
  ji	=0;
  for(i=0; i < ESNa1;i++){
    y0		=ESaZ0[n]+rsT[ki];
    b		=ESsb[in];
    for(j=0; j < ESNp1; j++){
      ECz[ji]	=y0+b*ESsn1[j];
      ji++;
    }
    ki++;
    in++;
  }
  return(0);
}

#ifdef Hspline
int four_s(gc,gs,m,g,d2g,x,n)
     double gc[],gs[],g[],d2g[],x[];
     int*m,*n;
{
  int i,ii,k,kk;
  double h,rh,cs,sn,t;
  double dg1,dg2,d3g,gg,rg;
  double*fcs,*fsn;
  i	=(*m)+1;
  fcs	=(double*)malloc(2*i*sizeof(double));
  if(fcs==NULL){
    printf("Failure in memory allocation in four_s\n");
    exit(0);
  }
  fsn= fcs+i;
  k= *m;
  t= EZc2gp*x[0];
  cs= cos(t);
  sn= sin(t);
  t*= (double)k;
  fcs[k]= cos(t);
  fsn[k]= sin(t);
  gc[k]= 0.;
  gs[k]= 0.;
  for(kk= k,k--;k>=0;k--,kk--){
    gc[k]= 0.;
    gs[k]= 0.;
    fcs[k]= fcs[kk]*cs+fsn[kk]*sn;
    fsn[k]= fsn[kk]*cs-fcs[kk]*sn;
  }
  h= x[1]-x[0];
  dg2= (g[1]-g[0])/h-(d2g[1]+2.*d2g[0])*h*EZcr6;

  for(i= 0,ii= 1;i<(*n);i++,ii++){
    h= x[ii]-x[i];
    rh= 1./h;
    dg1= dg2;
    dg2= (g[ii]-g[i])*rh+(2.*d2g[ii]+d2g[i])*h*EZcr6;
    d3g= (d2g[ii]-d2g[i])*rh;
    for(k= 1;k<=(*m);k++){
      gg= EZc2gp*(double)k;
      rg= 1./gg;
      gc[k]-= (g[i]*fsn[k]+(dg1*fcs[k]-(d2g[i]*fsn[k]
					+d3g*fcs[k]*rg)*rg)*rg)*rg;
      gs[k]-= (-g[i]*fcs[k]+(dg1*fsn[k]+(d2g[i]*fcs[k]
					 -d3g*fsn[k]*rg)*rg)*rg)*rg;
    }
    k= *m;
    t= EZc2gp*x[ii];
    cs= cos(t);
    sn= sin(t);
    t*= (double)k;
    fcs[k]= cos(t);
    fsn[k]= sin(t);
    for(kk= k,k--;k>=0;k--,kk--)
      {
	fcs[k]= fcs[kk]*cs+fsn[kk]*sn;
	fsn[k]= fsn[kk]*cs-fcs[kk]*sn;
      }
    gc[0]+= (g[ii]+(-0.5*dg2+EZcr6*(d2g[ii]-0.25*d3g*h)*h)*h)*h;
    for(k= 1;k<=*m;k++)
      {
	gg= EZc2gp*(double)k;
	rg= 1./gg;
	gc[k]+= (g[ii]*fsn[k]+(dg2*fcs[k]-(d2g[ii]*fsn[k]
					   +d3g*fcs[k]*rg)*rg)*rg)*rg;
	gs[k]+= (-g[ii]*fcs[k]+(dg2*fsn[k]+(d2g[ii]*fcs[k]
					    -d3g*fsn[k]*rg)*rg)*rg)*rg;
      }
  }
  free(fcs);
  return(0);
}
#endif

int EZfour(double *gc,double *gs,double *g)
{
  int j,kj,k;
  double h;

  for(k=0; k < ESFp1; k++){
    gc[k]	=0.;
    gs[k]	=0.;
  }
  for(j=0; j < ESNp; j++){
    gc[0]	+=g[j];
    kj		=0;
    for(k=1; k < ESFp1; k++){
      kj	+=j;
      if(kj >= ESNp)
	kj	-= ESNp;
      gc[k]	+=g[j]*EScs1[kj];
      gs[k]	+=g[j]*ESsn1[kj];
    }
  }
  h	=1./(double)ESNp;
  gc[0]	*=h;
  for(k=1; k < ESFp1; k++){
    gc[k]	*=h;
    gs[k]	*=h;
  }
  if(2*ESFp == ESNp){
    gc[ESFp]		*=EZcr2;
  }
  return(0);
}

int ECRescue2DGeomX()
{
  int j;
  for(j=0; j < ESNp1; j++){
    x0sepE[j]	=EZx0sep[j];
    x1sepE[j]	=EZx1sep[j];
    x2sepE[j]	=EZx2sep[j];
  }

  d_vbE	=EZd_vb;
  b_XE 	=b_X;
  z_X0E	=z_X0;
  z_X2E	=z_X2;
  dx_0xE	=EZdx_0x;
  d2x_0xE	=d2x_0x;
  dx_0oE	=EZdx_0o;
  d2x_0oE	=d2x_0o;

  return(0);
}

int ECRescue2DGeom(int n)
{
  int i,m,k,km,in,ki,K;

  in		=ESNa1*n;
  K		=ESnAF*n;
  R0E		=ESaR0[n];
  Z0E		=ESaZ0[n];
  ki		=K;
  for(i=0; i < ESNa1; i++){
    bE[i]	=ESsb  [in];
    bE1[i]	=ESsb1a[in];	
    bE2[i]	=ESsb2a[in];
    rcE[i]	=rcT  [ki];
    rcE1[i]	=rcT1a[ki];
    rcE2[i]	=rcT2a[ki];
    rsE[i]	=rsT  [ki];
    rsE1[i]	=rsT1a[ki];
    rsE2[i]	=rsT2a[ki];
    in++;
    ki++;
  }
  km		=ESNa1;
  for(m=1; m < ESFp1; m++){
    for(i=0; i < ESNa1; i++){
      rcE [km]	=rcT  [ki];
      rcE1[km]	=rcT1a[ki];
      rcE2[km]	=rcT2a[ki];
      rsE [km]	=rsT  [ki];
      rsE1[km]	=rsT1a[ki];
      rsE2[km]	=rsT2a[ki];
      ki++;
      km++;
    }
  }
  if(ESiAx != 0){
    ECRescue2DGeomX();
  }
  return(0);
}

int ECSaveESC()
{
  FILE *fpa,*fpb;
  static char RunIDb[16],RunIDa[16];
  static int Fl=0,L;
  int long Fpnt;
  int k;

  if(Fl == 0){
#ifdef XWIN
    {
      char *lD,*lN;
      CbGetCurrentCbName(&lD,&lN);
      lN	=ESRunID;
      while(*lD != '\0'){
	*lN	=*lD++;
	if(*lN == '/'){
	  lN	=ESRunID;
	}
	else{
	  lN++;
	}
      }
      *lN	='\0';
    }
#endif
    L	=CbStrCpy(RunIDb,ESRunID);
    L	=CbStrCpy(RunIDa,ESRunID);
    strcpy(RunIDb+L,".es");
    strcpy(RunIDa+L,".esa");
    Fl	=1;
  }
  fpa	=fopen(RunIDa,"r");
  fpb	=fopen(RunIDb,"r");
  Fl=fpa != NULL && fpb != NULL ? 3 : 1;
  if(fpa != NULL) fclose(fpa);
  if(fpb != NULL) fclose(fpb);

  if(Fl == 1){
    fpa	=fopen(RunIDa,"w");
    fpb	=fopen(RunIDb,"w");
    Fl	=2;
  }
  else{
    fpa	=fopen(RunIDa,"a");
    fpb	=fopen(RunIDb,"a");
  }
  if(fpa == NULL){
    printf("%s - ASCII recovery file cannot be open\n",RunIDa);
    return(-1);
  }
  if(fpb == NULL){
    printf("%s - BINARY recovery file cannot be open\n",RunIDb);
    fclose(fpa);
    return(-1);
  }
  if(Fl == 2){
    fprintf(fpa,"ESNa1=%d ESFp1=%d\n",ESNa1,ESFp1);
    Fl	=3;
  }
  Fpnt	=ftell(fpb);
  fprintf(fpa,"time=%12.5e RadC=%d InPr=%d InCr=%d "
	  "nPr=%d nCr=%d Fpnt=0x%8.8x 0x%8.8x\n"
	  ,ESTime,ESEqSolvRadC,ESEqSolvInPr,ESEqSolvInCr
	  ,ESnDtPr,ESnDtCr,Fpnt,ESFail);
  {
    double *ydPr,*ydCr,*xdPr,*xdCr;
    int *kdPr,*kdCr;
    ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
    fwrite(xdPr,sizeof(double),ESnDtPr+1,fpb);
    fwrite(ydPr,sizeof(double),ESnDtPr+1,fpb);
    fwrite(kdPr,sizeof(int),ESnDtPr+1,fpb);
    
    fwrite(xdCr,sizeof(double),ESnDtCr+1,fpb);
    fwrite(ydCr,sizeof(double),ESnDtCr+1,fpb);
#ifdef H
    if(1 || ESEqSolvInCr != 6){
    }
    else{
      double y[ESnDtCr+1];
      int i;
      k	=ESnDtCr+1;
      for(i=0; i < k; i++){
	y[i]	=1./ydCr[i];
      }
      fwrite(y,sizeof(double),k,fpb);
    }
#endif
    fwrite(kdCr,sizeof(int),ESnDtCr+1,fpb);
  }    
  fwrite(ESaF,sizeof(double),ESNa1,fpb);
  fwrite(ESaF1a,sizeof(double),ESNa1,fpb);
  fwrite(ESaF2a,sizeof(double),ESNa1,fpb);

  fwrite(ESgY,sizeof(double),ESNa1,fpb);
  fwrite(ESgY1a,sizeof(double),ESNa1,fpb);
  fwrite(ESgY2a,sizeof(double),ESNa1,fpb);

  fwrite(ESgF0,sizeof(double),ESNa1,fpb);
  fwrite(ESgF01a,sizeof(double),ESNa1,fpb);
  fwrite(ESqgF,sizeof(double),ESNa1,fpb);
  fwrite(ESqgF1a,sizeof(double),ESNa1,fpb);
  fwrite(ESqgF2a,sizeof(double),ESNa1,fpb);
  fwrite(ESgF,sizeof(double),ESNa1,fpb);
  fwrite(ESgF1a,sizeof(double),ESNa1,fpb);
  fwrite(ESgF2a,sizeof(double),ESNa1,fpb);
  fwrite(ESdgY,sizeof(double),ESNa1,fpb);
  fwrite(ESdgY1a,sizeof(double),ESNa1,fpb);
  fwrite(ESdgY2a,sizeof(double),ESNa1,fpb);
  fwrite(ESgY0,sizeof(double),ESNa1,fpb);
  fwrite(ESgY01a,sizeof(double),ESNa1,fpb);
  fwrite(ESgY02a,sizeof(double),ESNa1,fpb);
  fwrite(ESFF,sizeof(double),ESNa1,fpb);
  fwrite(ESFF1a,sizeof(double),ESNa1,fpb);
  fwrite(ESFF2a,sizeof(double),ESNa1,fpb);
  fwrite(ESsq,sizeof(double),ESNa1,fpb);
  fwrite(ESsq1a,sizeof(double),ESNa1,fpb);
  fwrite(ESsq2a,sizeof(double),ESNa1,fpb);
  fwrite(ESgm,sizeof(double),ESNa1,fpb);
  fwrite(ESgm1a,sizeof(double),ESNa1,fpb);
  fwrite(ESgm2a,sizeof(double),ESNa1,fpb);
  fwrite(ESaT,sizeof(double),ESNa1,fpb);
  fwrite(ESaT1a,sizeof(double),ESNa1,fpb);
  fwrite(ESaT2a,sizeof(double),ESNa1,fpb);
  fwrite(ESjB,sizeof(double),ESNa1,fpb);
  fwrite(ESjB1a,sizeof(double),ESNa1,fpb);
  fwrite(ESjB2a,sizeof(double),ESNa1,fpb);
  fwrite(ESaJ,sizeof(double),ESNa1,fpb);
  fwrite(ESaJ1a,sizeof(double),ESNa1,fpb);
  fwrite(ESaJ2a,sizeof(double),ESNa1,fpb);
  fwrite(ESjs,sizeof(double),ESNa1,fpb);
  fwrite(ESjs1a,sizeof(double),ESNa1,fpb);
  fwrite(ESjs2a,sizeof(double),ESNa1,fpb);
  fwrite(ESjp,sizeof(double),ESNa1,fpb);
  fwrite(ESjp1a,sizeof(double),ESNa1,fpb);
  fwrite(ESjp2a,sizeof(double),ESNa1,fpb);
  fwrite(ESPs,sizeof(double),ESNa1,fpb);
  fwrite(ESPs1a,sizeof(double),ESNa1,fpb);
  fwrite(ESPs2a,sizeof(double),ESNa1,fpb);
  fwrite(ESsp,sizeof(double),ESNa1,fpb);
  fwrite(ESsp1a,sizeof(double),ESNa1,fpb);
  fwrite(ESsp2a,sizeof(double),ESNa1,fpb);
  fwrite(ESjb,sizeof(double),ESNa1,fpb);
  fwrite(ESjb1a,sizeof(double),ESNa1,fpb);
  fwrite(ESjb2a,sizeof(double),ESNa1,fpb);
  fwrite(ESgb,sizeof(double),ESNa1,fpb);
  fwrite(ESgb1a,sizeof(double),ESNa1,fpb);
  fwrite(ESgb2a,sizeof(double),ESNa1,fpb);
  fwrite(ESgB,sizeof(double),ESNa1,fpb);
  fwrite(ESgB1a,sizeof(double),ESNa1,fpb);
  fwrite(ESgB2a,sizeof(double),ESNa1,fpb);
  fwrite(ESqV,sizeof(double),ESNa1,fpb);
  fwrite(ESqV1a,sizeof(double),ESNa1,fpb);
  fwrite(ESqV2a,sizeof(double),ESNa1,fpb);
  fwrite(ESqgY,sizeof(double),ESNa1,fpb);
  fwrite(ESqgY1a,sizeof(double),ESNa1,fpb);
  fwrite(ESqgY2a,sizeof(double),ESNa1,fpb);

  fwrite(&ESRBt,sizeof(double),1,fpb);
  fwrite(&ESRext,sizeof(double),1,fpb);
  fwrite(ESaR0,sizeof(double),1,fpb);
  fwrite(ESaZ0,sizeof(double),1,fpb);
  fwrite(ESsb  ,sizeof(double),ESNa1,fpb);
  fwrite(ESsb1a,sizeof(double),ESNa1,fpb);
  fwrite(ESsb2a,sizeof(double),ESNa1,fpb);
  k	=ESNa1*ESFp1;
  fwrite(rcT  ,sizeof(double),k,fpb);
  fwrite(rcT1a,sizeof(double),k,fpb);
  fwrite(rcT2a,sizeof(double),k,fpb);
  fwrite(rsT  ,sizeof(double),k,fpb);
  fwrite(rsT1a,sizeof(double),k,fpb);
  fwrite(rsT2a,sizeof(double),k,fpb);
  fclose(fpa);
  fclose(fpb);
  return(0);
}

int ECSetRestoreESCTime(double *t1,double *t2,double *t,int k1,int k2,int k)
{
  static double tm[256];
  double s,ss;
  FILE *fpa;
  static char RunIDa[16];
  static int n=0,i1=0,i2=0,i0=0;
  int i;
  char ch;

  if(n == 0 || (k1 == 1 && k2 == 1)){
#ifdef XWIN
    {
      char *lD,*lN;
      CbGetCurrentCbName(&lD,&lN);
      lN	=ESRunID;
      while(*lD != '\0'){
	*lN	=*lD++;
	if(*lN == '/'){
	  lN	=ESRunID;
	}
	else{
	  lN++;
	}
      }
      *lN	='\0';
    }
#endif
    i	=CbStrCpy(RunIDa,ESRunID);
    strcpy(RunIDa+i,".esa");
    fpa	=fopen(RunIDa,"r");

    if(fpa == NULL){
      printf("%s - ASCII recovery file cannot be open\n",RunIDa);
      return(-1);
    }
    while(feof(fpa) == 0 && (ch=getc(fpa)) != '\n') ;
    n	=0;
    while(fscanf(fpa,"time=%lg",tm+n) == 1){
      while(feof(fpa) == 0 && (ch=getc(fpa)) != '\n') ;
      n++;
    }
    fclose(fpa);
  }
  EZout("sI","Number of time slices",n);
  if(n == 0)  return(1);

  switch(k1){
  case 1:
    i1	=0;
    break;
  case 2:
    if(i1) i1--;
    break;
  case 3:
    i1++;
    if(i1 >= n)  i1=n-1;
    break;
  default:
    break;
  }
  switch(k2){
  case 1:
    i2	=n-1;
    break;
  case 2:
    if(i2 > i1) i2--;
    break;
  case 3:
    if(i2 < n-1) i2++;
    break;
  default:
    break;
  }

  switch(k){
  case 1:
    ss	=fabs(*t-tm[0]);
    i0	=0;
    for(i=1; i < n; i++){
      s	=fabs(*t-tm[i]);
      if(ss > s){
	ss	=s;
	i0	=i;
      }
    }
    break;
  case 2:
    if(i0 > i1) i0--;
    break;
  case 3:
    if(i0 < i2) i0++;
    break;
  default:
    break;
  }
  if(i2 <  i1) i2=i1;
  if(i0 <  i1) i0=i1;
  if(i0 >  i2) i0=i2;

  *t1	=tm[i1];
  *t2	=tm[i2];
  *t	=tm[i0];
  return(0);
}

void ecssetpointflag_()
{
  extern int FlagPoints;
  FlagPoints=1; /* 0 - number of boundary points <= 12;
		   1 - number of boundary points > 12;
		   2 - boundary is specified by Fourier harmonics
		   */
  return;
}

int ECRestoreESC(double ttime)
{
  FILE *fpa,*fpb;
  static char RunIDb[16],RunIDa[16];
  static int Fl=0,L;
  int long Fpnt;
  int Na1,Fp1,k,i,j,kk;
  double t;

  if(Fl == 0){
#ifdef XWIN
    {
      char *lD,*lN;
      CbGetCurrentCbName(&lD,&lN);
      printf("lN=<%s><%s><%s>\n",lN,lD,ESRunID);
      lN	=ESRunID;
      while(*lD != '\0'){
	*lN	=*lD++;
	if(*lN == '/'){
	  lN	=ESRunID;
	}
	else{
	  lN++;
	}
      }
      *lN	='\0';
    }
#endif
    L	=CbStrCpy(RunIDb,ESRunID);
    L	=CbStrCpy(RunIDa,ESRunID);
    strcpy(RunIDb+L,".es");
    strcpy(RunIDa+L,".esa");
    Fl	=1;
  }
  fpa	=fopen(RunIDa,"r");
  fpb	=fopen(RunIDb,"r");
  if(fpa == NULL){
    printf("%s - ASCII recovery file cannot be open\n",RunIDa);
    return(-1);
  }
  if(fpb == NULL){
    printf("%s - BINARY recovery file cannot be open\n",RunIDb);
    fclose(fpa);
    return(-1);
  }
  fscanf(fpa,"ESNa1=%d ESFp1=%d\n",&Na1,&Fp1);
#ifdef H
  printf("ESNa1=%d ESFp1=%d time=%10.3e\n",Na1,Fp1,ttime);
#endif
  t	=-1e+20;
  Fpnt	=0.;
  k	=8;
  while(fabs(t-ttime) > 1e-6 && k == 8){
    k=fscanf(fpa,"time=%lg RadC=%1d InPr=%2d InCr=%2d "
	     "nPr=%d nCr=%d Fpnt=0x%x 0x%x\n"
	     ,&t,&ESEqSolvRadC,&ESEqSolvInPr,&ESEqSolvInCr
	     ,&ESnDtPr,&ESnDtCr,&Fpnt,&ESFail);
#ifdef H
    printf("k=%d\n",k);
    printf("time=%12.5e RadC=%1d InPr=%2d InCr=%2d "
	   "nPr=%2d nCr=%2d Fpnt=0x%8.8x\n"
	   ,ESTime,ESEqSolvRadC,ESEqSolvInPr,ESEqSolvInCr
	   ,ESnDtPr,ESnDtCr,Fpnt);
#endif
  }
  if(fabs(t-ttime) > 1e-6){
    printf("time=%12.5e is absent in %s - recovery file\n",ttime,RunIDa);
    fclose(fpa);
    fclose(fpb);
    return(-1);
  }
  if(k != 8){
    printf("Fpnt is absent in %s - recovery file\n",RunIDa);
    fclose(fpa);
    fclose(fpb);
    return(-1);
  }
  fseek(fpb,Fpnt,SEEK_SET);

  {
    double *ydPr,*ydCr,*xdPr,*xdCr;
    int *kdPr,*kdCr;

    k	=ESnDtPr+1;
    ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
    fread(xdPr,sizeof(double),k,fpb);
    fread(ydPr,sizeof(double),k,fpb);
    fread(kdPr,sizeof(int),k,fpb); 
    for(i=0; i < k; i++) ESxDtPr[i]=xdPr[i];
    while(i < 257){
      kdPr[i]	=0;
      i++;
    }
    k	=ESnDtCr+1;
    fread(xdCr,sizeof(double),k,fpb);
    fread(ydCr,sizeof(double),k,fpb);
    fread(kdCr,sizeof(int),k,fpb);
    for(i=0; i < k; i++) ESxDtCr[i]=xdCr[i];
    while(i < 257){
      kdCr[i]	=0;
      i++;
    }
    ESInpPr2Spl();
  }    
  fread(ESaF,sizeof(double),ESNa1,fpb);
  fread(ESaF1a,sizeof(double),ESNa1,fpb);
  fread(ESaF2a,sizeof(double),ESNa1,fpb);

  fread(ESgY,sizeof(double),ESNa1,fpb);
  fread(ESgY1a,sizeof(double),ESNa1,fpb);
  fread(ESgY2a,sizeof(double),ESNa1,fpb);

  fread(ESgF0,sizeof(double),ESNa1,fpb);
  fread(ESgF01a,sizeof(double),ESNa1,fpb);
  fread(ESqgF,sizeof(double),ESNa1,fpb);
  fread(ESqgF1a,sizeof(double),ESNa1,fpb);
  fread(ESqgF2a,sizeof(double),ESNa1,fpb);
  fread(ESgF,sizeof(double),ESNa1,fpb);
  fread(ESgF1a,sizeof(double),ESNa1,fpb);
  fread(ESgF2a,sizeof(double),ESNa1,fpb);
  fread(ESdgY,sizeof(double),ESNa1,fpb);
  fread(ESdgY1a,sizeof(double),ESNa1,fpb);
  fread(ESdgY2a,sizeof(double),ESNa1,fpb);
  fread(ESgY0,sizeof(double),ESNa1,fpb);
  fread(ESgY01a,sizeof(double),ESNa1,fpb);
  fread(ESgY02a,sizeof(double),ESNa1,fpb);
  fread(ESFF,sizeof(double),ESNa1,fpb);
  fread(ESFF1a,sizeof(double),ESNa1,fpb);
  fread(ESFF2a,sizeof(double),ESNa1,fpb);
  fread(ESsq,sizeof(double),ESNa1,fpb);
  fread(ESsq1a,sizeof(double),ESNa1,fpb);
  fread(ESsq2a,sizeof(double),ESNa1,fpb);
  fread(ESgm,sizeof(double),ESNa1,fpb);
  fread(ESgm1a,sizeof(double),ESNa1,fpb);
  fread(ESgm2a,sizeof(double),ESNa1,fpb);
  fread(ESaT,sizeof(double),ESNa1,fpb);
  fread(ESaT1a,sizeof(double),ESNa1,fpb);
  fread(ESaT2a,sizeof(double),ESNa1,fpb);
  fread(ESjB,sizeof(double),ESNa1,fpb);
  fread(ESjB1a,sizeof(double),ESNa1,fpb);
  fread(ESjB2a,sizeof(double),ESNa1,fpb);
  fread(ESaJ,sizeof(double),ESNa1,fpb);
  fread(ESaJ1a,sizeof(double),ESNa1,fpb);
  fread(ESaJ2a,sizeof(double),ESNa1,fpb);
  fread(ESjs,sizeof(double),ESNa1,fpb);
  fread(ESjs1a,sizeof(double),ESNa1,fpb);
  fread(ESjs2a,sizeof(double),ESNa1,fpb);
  fread(ESjp,sizeof(double),ESNa1,fpb);
  fread(ESjp1a,sizeof(double),ESNa1,fpb);
  fread(ESjp2a,sizeof(double),ESNa1,fpb);
  fread(ESPs,sizeof(double),ESNa1,fpb);
  fread(ESPs1a,sizeof(double),ESNa1,fpb);
  fread(ESPs2a,sizeof(double),ESNa1,fpb);
  fread(ESsp,sizeof(double),ESNa1,fpb);
  fread(ESsp1a,sizeof(double),ESNa1,fpb);
  fread(ESsp2a,sizeof(double),ESNa1,fpb);
  fread(ESjb,sizeof(double),ESNa1,fpb);
  fread(ESjb1a,sizeof(double),ESNa1,fpb);
  fread(ESjb2a,sizeof(double),ESNa1,fpb);
  fread(ESgb,sizeof(double),ESNa1,fpb);
  fread(ESgb1a,sizeof(double),ESNa1,fpb);
  fread(ESgb2a,sizeof(double),ESNa1,fpb);
  fread(ESgB,sizeof(double),ESNa1,fpb);
  fread(ESgB1a,sizeof(double),ESNa1,fpb);
  fread(ESgB2a,sizeof(double),ESNa1,fpb);
  fread(ESqV,sizeof(double),ESNa1,fpb);
  fread(ESqV1a,sizeof(double),ESNa1,fpb);
  fread(ESqV2a,sizeof(double),ESNa1,fpb);
  fread(ESqgY,sizeof(double),ESNa1,fpb);
  fread(ESqgY1a,sizeof(double),ESNa1,fpb);
  fread(ESqgY2a,sizeof(double),ESNa1,fpb);

  fread(&ESRBt,sizeof(double),1,fpb);
  fread(&ESRext,sizeof(double),1,fpb);
  ESBt	=ESRBt/ESRext;
  fread(ESaR0,sizeof(double),1,fpb);
  fread(ESaZ0,sizeof(double),1,fpb);
  fread(ESsb  ,sizeof(double),ESNa1,fpb);
  fread(ESsb1a,sizeof(double),ESNa1,fpb);
  fread(ESsb2a,sizeof(double),ESNa1,fpb);
  k	=Na1*Fp1;
  fread(rcT  ,sizeof(double),k,fpb);
  fread(rcT1a,sizeof(double),k,fpb);
  fread(rcT2a,sizeof(double),k,fpb);
  fread(rsT  ,sizeof(double),k,fpb);
  fread(rsT1a,sizeof(double),k,fpb);
  fread(rsT2a,sizeof(double),k,fpb);

  fclose(fpa);
  fclose(fpb);

  R0		=ESaR0[0];
  Z0		=ESaZ0[0];
  for(i=0; i < ESNa1; i++){
    ECb  [i]	=ESsb[i];
    ECb1a[i]	=ESsb1a[i];	
    ECb2a[i]	=ESsb2a[i];
    EZrcs  [i]	=rcT[i];
    EZd1rcs[i]	=rcT1a[i];
    EZd2rcs[i]	=rcT2a[i];
    EZz0  [i]	=rsT[i];
    EZd1z0[i]	=rsT1a[i];
    EZd2z0[i]	=rsT2a[i];
  }
  kk		=ESNa1;
  for(k=1; k < ESFp1; k++){
    for(i=0; i < ESNa1; i++){
      EZrcs  [kk]	=rcT [kk];
      EZd1rcs[kk]	=rcT1a[kk];
      EZd2rcs[kk]	=rcT2a[kk];
      EZrsn  [kk]	=rsT [kk];
      EZd1rsn[kk]	=rsT1a[kk];
      EZd2rsn[kk]	=rsT2a[kk];
      kk++;
    }
  }
  ESMakePlVdata();
  return(0);
}

int ES2HardDrive(double ttime)
{
  FILE *fpa,*fpb;
  static char RunIDb[16],RunIDa[16];
  static int Fl=0,L;
  int long Fpnt;
  int k;

  if(Fl == 0){
    L	=CbStrCpy(RunIDb,ESRunID);
    L	=CbStrCpy(RunIDa,ESRunID);
    strcpy(RunIDb+L,".lz");
    strcpy(RunIDa+L,".lza");
    Fl	=1;
  }
  if(Fl == 1){
    fpa	=fopen(RunIDa,"w");
    fpb	=fopen(RunIDb,"w");
    Fl	=2;
  }
  else{
    fpa	=fopen(RunIDa,"a");
    fpb	=fopen(RunIDb,"a");
  }
  if(fpa == NULL){
    printf("%s - ASCII recovery file cannot be open\n",RunIDa);
    return(-1);
  }
  if(fpb == NULL){
    printf("%s - BINARY recovery file cannot be open\n",RunIDb);
    fclose(fpa);
    return(-1);
  }
  if(Fl == 2){
    fprintf(fpa,"ESNa1=%d ESFp1=%d\n",ESNa1,ESFp1);
    Fl	=3;
  }
  Fpnt	=ftell(fpb);
  fprintf(fpa,"time=%12.5e RadC=%1d InPr=%2d InCr=%2d "
	  "nPr=%2d nCr=%2d Fpnt=0x%8.8x\n"
	  ,ttime,ESEqSolvRadC,ESEqSolvInPr,ESEqSolvInCr
	  ,ESnDtPr,ESnDtCr,Fpnt);
  {
    double *ydPr,*ydCr,*xdPr,*xdCr;
    int *kdPr,*kdCr;
    ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
    fwrite(xdPr,sizeof(double),ESnDtPr+1,fpb);
    fwrite(ydPr,sizeof(double),ESnDtPr+1,fpb);
    fwrite(kdPr,sizeof(int),ESnDtPr+1,fpb);
    
    fwrite(xdCr,sizeof(double),ESnDtCr+1,fpb);
    fwrite(ydCr,sizeof(double),ESnDtCr+1,fpb);
    fwrite(kdCr,sizeof(int),ESnDtCr+1,fpb);
  }    
  fwrite(ESaF,sizeof(double),ESNa1,fpb);
  fwrite(ESaF1a,sizeof(double),ESNa1,fpb);
  fwrite(ESaF2a,sizeof(double),ESNa1,fpb);

  fwrite(ESgY,sizeof(double),ESNa1,fpb);
  fwrite(ESgY1a,sizeof(double),ESNa1,fpb);
  fwrite(ESgY2a,sizeof(double),ESNa1,fpb);

  fwrite(ESgF0,sizeof(double),ESNa1,fpb);
  fwrite(ESgF01a,sizeof(double),ESNa1,fpb);
  fwrite(ESqgF,sizeof(double),ESNa1,fpb);
  fwrite(ESqgF1a,sizeof(double),ESNa1,fpb);
  fwrite(ESqgF2a,sizeof(double),ESNa1,fpb);
  fwrite(ESgF,sizeof(double),ESNa1,fpb);
  fwrite(ESgF1a,sizeof(double),ESNa1,fpb);
  fwrite(ESgF2a,sizeof(double),ESNa1,fpb);
  fwrite(ESdgY,sizeof(double),ESNa1,fpb);
  fwrite(ESdgY1a,sizeof(double),ESNa1,fpb);
  fwrite(ESdgY2a,sizeof(double),ESNa1,fpb);
  fwrite(ESgY0,sizeof(double),ESNa1,fpb);
  fwrite(ESgY01a,sizeof(double),ESNa1,fpb);
  fwrite(ESgY02a,sizeof(double),ESNa1,fpb);
  fwrite(ESFF,sizeof(double),ESNa1,fpb);
  fwrite(ESFF1a,sizeof(double),ESNa1,fpb);
  fwrite(ESFF2a,sizeof(double),ESNa1,fpb);
  fwrite(ESsq,sizeof(double),ESNa1,fpb);
  fwrite(ESsq1a,sizeof(double),ESNa1,fpb);
  fwrite(ESsq2a,sizeof(double),ESNa1,fpb);
  fwrite(ESgm,sizeof(double),ESNa1,fpb);
  fwrite(ESgm1a,sizeof(double),ESNa1,fpb);
  fwrite(ESgm2a,sizeof(double),ESNa1,fpb);
  fwrite(ESaT,sizeof(double),ESNa1,fpb);
  fwrite(ESaT1a,sizeof(double),ESNa1,fpb);
  fwrite(ESaT2a,sizeof(double),ESNa1,fpb);
  fwrite(ESjB,sizeof(double),ESNa1,fpb);
  fwrite(ESjB1a,sizeof(double),ESNa1,fpb);
  fwrite(ESjB2a,sizeof(double),ESNa1,fpb);
  fwrite(ESaJ,sizeof(double),ESNa1,fpb);
  fwrite(ESaJ1a,sizeof(double),ESNa1,fpb);
  fwrite(ESaJ2a,sizeof(double),ESNa1,fpb);
  fwrite(ESjs,sizeof(double),ESNa1,fpb);
  fwrite(ESjs1a,sizeof(double),ESNa1,fpb);
  fwrite(ESjs2a,sizeof(double),ESNa1,fpb);
  fwrite(ESjp,sizeof(double),ESNa1,fpb);
  fwrite(ESjp1a,sizeof(double),ESNa1,fpb);
  fwrite(ESjp2a,sizeof(double),ESNa1,fpb);
  fwrite(ESPs,sizeof(double),ESNa1,fpb);
  fwrite(ESPs1a,sizeof(double),ESNa1,fpb);
  fwrite(ESPs2a,sizeof(double),ESNa1,fpb);
  fwrite(ESsp,sizeof(double),ESNa1,fpb);
  fwrite(ESsp1a,sizeof(double),ESNa1,fpb);
  fwrite(ESsp2a,sizeof(double),ESNa1,fpb);
  fwrite(ESjb,sizeof(double),ESNa1,fpb);
  fwrite(ESjb1a,sizeof(double),ESNa1,fpb);
  fwrite(ESjb2a,sizeof(double),ESNa1,fpb);
  fwrite(ESgb,sizeof(double),ESNa1,fpb);
  fwrite(ESgb1a,sizeof(double),ESNa1,fpb);
  fwrite(ESgb2a,sizeof(double),ESNa1,fpb);
  fwrite(ESgB,sizeof(double),ESNa1,fpb);
  fwrite(ESgB1a,sizeof(double),ESNa1,fpb);
  fwrite(ESgB2a,sizeof(double),ESNa1,fpb);
  fwrite(ESqV,sizeof(double),ESNa1,fpb);
  fwrite(ESqV1a,sizeof(double),ESNa1,fpb);
  fwrite(ESqV2a,sizeof(double),ESNa1,fpb);
  fwrite(ESqgY,sizeof(double),ESNa1,fpb);
  fwrite(ESqgY1a,sizeof(double),ESNa1,fpb);
  fwrite(ESqgY2a,sizeof(double),ESNa1,fpb);

  fwrite(&ESRBt,sizeof(double),1,fpb);
  fwrite(&ESRext,sizeof(double),1,fpb);

  fwrite(ESaR0,sizeof(double),1,fpb);
  fwrite(ESaZ0,sizeof(double),1,fpb);
  fwrite(ESsb  ,sizeof(double),ESNa1,fpb);
  fwrite(ESsb1a,sizeof(double),ESNa1,fpb);
  fwrite(ESsb2a,sizeof(double),ESNa1,fpb);
  k	=ESNa1*ESFp1;
  fwrite(rcT  ,sizeof(double),k,fpb);
  fwrite(rcT1a,sizeof(double),k,fpb);
  fwrite(rcT2a,sizeof(double),k,fpb);
  fwrite(rsT  ,sizeof(double),k,fpb);
  fwrite(rsT1a,sizeof(double),k,fpb);
  fwrite(rsT2a,sizeof(double),k,fpb);
  fclose(fpa);
  fclose(fpb);
  return(0);
}

void es2harddrive_(double *ttime)
{
  ES2HardDrive(*ttime);
}

void esfromharddrive_(double *ttime)
{
  ESFromHardDrive(*ttime);
}

int ESFromHardDrive(double ttime)
{
  FILE *fpa,*fpb;
  static char RunIDb[16],RunIDa[16];
  static int Fl=0,L;
  int long Fpnt;
  int Na1,Fp1,k,i,kk,np,nj;
  double t;

  if(Fl == 0){
    L	=CbStrCpy(RunIDb,ESRunID);
    L	=CbStrCpy(RunIDa,ESRunID);
    strcpy(RunIDb+L,".lz");
    strcpy(RunIDa+L,".lza");
    Fl	=1;
  }
  fpa	=fopen(RunIDa,"r");
  fpb	=fopen(RunIDb,"r");
  if(fpa == NULL){
    printf("%s - ASCII recovery file cannot be open\n",RunIDa);
    return(-1);
  }
  if(fpb == NULL){
    printf("%s - BINARY recovery file cannot be open\n",RunIDb);
    fclose(fpa);
    return(-1);
  }
  fscanf(fpa,"ESNa1=%d ESFp1=%d\n",&Na1,&Fp1);
#ifdef H
  printf("ESNa1=%d ESFp1=%d time=%10.3e\n",Na1,Fp1,ttime);
#endif
  t	=-1e+20;
  Fpnt	=0.;
  k	=7;
  while(fabs(t-ttime) > 1e-6 && k == 7){
    k=fscanf(fpa,"time=%lg RadC=%1d InPr=%2d InCr=%2d "
	     "nPr=%2d nCr=%2d Fpnt=0x%x\n"
	     ,&t,&ESEqSolvRadC,&ESEqSolvInPr,&ESEqSolvInCr
	     ,&ESnDtPr,&ESnDtCr,&Fpnt);
    np	=ESnDtPr+1;
    nj	=ESnDtCr+1;
#ifdef H
    printf("k=%d\n",k);
    printf("time=%12.5e RadC=%1d InPr=%2d InCr=%2d "
	   "nPr=%2d nCr=%2d Fpnt=0x%8.8x\n"
	   ,ESTime,ESEqSolvRadC,ESEqSolvInPr,ESEqSolvInCr
	   ,ESnDtPr,ESnDtCr,Fpnt);
#endif
  }
  if(fabs(t-ttime) > 1e-6){
    printf("time=%12.5e is absent in %s - recovery file\n",ttime,RunIDa);
    fclose(fpa);
    fclose(fpb);
    return(-1);
  }
  if(k != 7){
    printf("Fpnt is absent in %s - recovery file\n",RunIDa);
    fclose(fpa);
    fclose(fpb);
    return(-1);
  }
  fseek(fpb,Fpnt,SEEK_SET);

  {
    double *ydPr,*ydCr,*xdPr,*xdCr;
    int *kdPr,*kdCr;

    ESGetPlPrfAddr(&ydPr,&ydCr,&xdPr,&xdCr,&kdPr,&kdCr);
    fread(xdPr,sizeof(double),ESnDtPr+1,fpb);
    fread(ydPr,sizeof(double),ESnDtPr+1,fpb);
    fread(kdPr,sizeof(int),ESnDtPr+1,fpb);
    k	=ESnDtPr+1;
    for(i=0; i < k; i++){
      ESxDtPr[i]	=xdPr[i];
    }
    while(i < 257){
      kdPr[i]	=0;
      i++;
    }
    fread(xdCr,sizeof(double),ESnDtCr+1,fpb);
    fread(ydCr,sizeof(double),ESnDtCr+1,fpb);
    fread(kdCr,sizeof(int),ESnDtCr+1,fpb);
    k	=ESnDtCr+1;
    for(i=0; i < k; i++){
      ESxDtCr[i]	=xdCr[i];
    }
    while(i < 257){
      kdCr[i]	=0;
      i++;
    }
    ESInputPrData2Pr();
    ESInpPlPrCr();
  }    
  fread(ESaF,sizeof(double),ESNa1,fpb);
  fread(ESaF1a,sizeof(double),ESNa1,fpb);
  fread(ESaF2a,sizeof(double),ESNa1,fpb);
  fread(ESgY,sizeof(double),ESNa1,fpb);
  fread(ESgY1a,sizeof(double),ESNa1,fpb);
  fread(ESgY2a,sizeof(double),ESNa1,fpb);

  fread(ESgF0,sizeof(double),ESNa1,fpb);
  fread(ESgF01a,sizeof(double),ESNa1,fpb);
  fread(ESqgF,sizeof(double),ESNa1,fpb);
  fread(ESqgF1a,sizeof(double),ESNa1,fpb);
  fread(ESqgF2a,sizeof(double),ESNa1,fpb);
  fread(ESgF,sizeof(double),ESNa1,fpb);
  fread(ESgF1a,sizeof(double),ESNa1,fpb);
  fread(ESgF2a,sizeof(double),ESNa1,fpb);
  fread(ESdgY,sizeof(double),ESNa1,fpb);
  fread(ESdgY1a,sizeof(double),ESNa1,fpb);
  fread(ESdgY2a,sizeof(double),ESNa1,fpb);
  fread(ESgY0,sizeof(double),ESNa1,fpb);
  fread(ESgY01a,sizeof(double),ESNa1,fpb);
  fread(ESgY02a,sizeof(double),ESNa1,fpb);
  fread(ESFF,sizeof(double),ESNa1,fpb);
  fread(ESFF1a,sizeof(double),ESNa1,fpb);
  fread(ESFF2a,sizeof(double),ESNa1,fpb);
  fread(ESsq,sizeof(double),ESNa1,fpb);
  fread(ESsq1a,sizeof(double),ESNa1,fpb);
  fread(ESsq2a,sizeof(double),ESNa1,fpb);
  fread(ESgm,sizeof(double),ESNa1,fpb);
  fread(ESgm1a,sizeof(double),ESNa1,fpb);
  fread(ESgm2a,sizeof(double),ESNa1,fpb);
  fread(ESaT,sizeof(double),ESNa1,fpb);
  fread(ESaT1a,sizeof(double),ESNa1,fpb);
  fread(ESaT2a,sizeof(double),ESNa1,fpb);
  fread(ESjB,sizeof(double),ESNa1,fpb);
  fread(ESjB1a,sizeof(double),ESNa1,fpb);
  fread(ESjB2a,sizeof(double),ESNa1,fpb);
  fread(ESaJ,sizeof(double),ESNa1,fpb);
  fread(ESaJ1a,sizeof(double),ESNa1,fpb);
  fread(ESaJ2a,sizeof(double),ESNa1,fpb);
  fread(ESjs,sizeof(double),ESNa1,fpb);
  fread(ESjs1a,sizeof(double),ESNa1,fpb);
  fread(ESjs2a,sizeof(double),ESNa1,fpb);
  fread(ESjp,sizeof(double),ESNa1,fpb);
  fread(ESjp1a,sizeof(double),ESNa1,fpb);
  fread(ESjp2a,sizeof(double),ESNa1,fpb);
  fread(ESPs,sizeof(double),ESNa1,fpb);
  fread(ESPs1a,sizeof(double),ESNa1,fpb);
  fread(ESPs2a,sizeof(double),ESNa1,fpb);
  fread(ESsp,sizeof(double),ESNa1,fpb);
  fread(ESsp1a,sizeof(double),ESNa1,fpb);
  fread(ESsp2a,sizeof(double),ESNa1,fpb);
  fread(ESjb,sizeof(double),ESNa1,fpb);
  fread(ESjb1a,sizeof(double),ESNa1,fpb);
  fread(ESjb2a,sizeof(double),ESNa1,fpb);
  fread(ESgb,sizeof(double),ESNa1,fpb);
  fread(ESgb1a,sizeof(double),ESNa1,fpb);
  fread(ESgb2a,sizeof(double),ESNa1,fpb);
  fread(ESgB,sizeof(double),ESNa1,fpb);
  fread(ESgB1a,sizeof(double),ESNa1,fpb);
  fread(ESgB2a,sizeof(double),ESNa1,fpb);
  fread(ESqV,sizeof(double),ESNa1,fpb);
  fread(ESqV1a,sizeof(double),ESNa1,fpb);
  fread(ESqV2a,sizeof(double),ESNa1,fpb);
  fread(ESqgY,sizeof(double),ESNa1,fpb);
  fread(ESqgY1a,sizeof(double),ESNa1,fpb);
  fread(ESqgY2a,sizeof(double),ESNa1,fpb);

  fread(&ESRBt,sizeof(double),1,fpb);
  fread(&ESRext,sizeof(double),1,fpb);
  ESBt	=ESRBt/ESRext;

  fread(ESaR0,sizeof(double),1,fpb);
  fread(ESaZ0,sizeof(double),1,fpb);
  fread(ESsb  ,sizeof(double),ESNa1,fpb);
  fread(ESsb1a,sizeof(double),ESNa1,fpb);
  fread(ESsb2a,sizeof(double),ESNa1,fpb);
  k	=Na1*Fp1;
  fread(rcT  ,sizeof(double),k,fpb);
  fread(rcT1a,sizeof(double),k,fpb);
  fread(rcT2a,sizeof(double),k,fpb);
  fread(rsT  ,sizeof(double),k,fpb);
  fread(rsT1a,sizeof(double),k,fpb);
  fread(rsT2a,sizeof(double),k,fpb);

  fclose(fpa);
  fclose(fpb);

  R0		=ESaR0[0];
  Z0		=ESaZ0[0];
  for(i=0; i < ESNa1; i++){
    ECb  [i]	=ESsb[i];
    ECb1a[i]	=ESsb1a[i];	
    ECb2a[i]	=ESsb2a[i];
    EZrcs  [i]	=rcT[i];
    EZd1rcs[i]	=rcT1a[i];
    EZd2rcs[i]	=rcT2a[i];
    EZz0  [i]	=rsT[i];
    EZd1z0[i]	=rsT1a[i];
    EZd2z0[i]	=rsT2a[i];
  }
  kk		=ESNa1;
  for(k=1; k < ESFp1; k++){
    for(i=0; i < ESNa1; i++){
      EZrcs  [kk]	=rcT [kk];
      EZd1rcs[kk]	=rcT1a[kk];
      EZd2rcs[kk]	=rcT2a[kk];
      EZrsn  [kk]	=rsT [kk];
      EZd1rsn[kk]	=rsT1a[kk];
      EZd2rsn[kk]	=rsT2a[kk];
      kk++;
    }
  }
  ESMakePlVdata();
  return(0);
}

int ECRestore2DGeomX()
{
  int j;
  for(j=0; j < ESNp1; j++){
    EZx0sep[j]	=x0sepE[j];
    EZx1sep[j]	=x1sepE[j];
    EZx2sep[j]	=x2sepE[j];
  }

  b_X 	=b_XE;
  EZd_vb	=(bE[ESNa]-b_X)/b_X;
  z_X0	=z_X0E;
  z_X2	=z_X2E;
  EZdx_0x	=dx_0xE;
  d2x_0x=d2x_0xE;
  EZdx_0o	=dx_0oE;
  d2x_0o=d2x_0oE;

  return(0);
}

int ECRestore2DGeom()
{
  int i,m,k,km,in;

  R0		=R0E;
  Z0		=Z0E;
  for(i=0; i < ESNa1; i++){
    ECb  [i]	=bE[i];
    ECb1a[i]	=bE1[i];	
    ECb2a[i]	=bE2[i];
    EZrcs  [i]	=rcE[i];
    EZd1rcs[i]	=rcE1[i];
    EZd2rcs[i]	=rcE2[i];
    EZz0  [i]	=rsE[i];
    EZd1z0[i]	=rsE1[i];
    EZd2z0[i]	=rsE2[i];
  }
  km		=ESNa1;
  for(m=1; m < ESFp1; m++){
    for(i=0; i < ESNa1; i++){
      EZrcs  [km]	=rcE [km];
      EZd1rcs[km]	=rcE1[km];
      EZd2rcs[km]	=rcE2[km];
      EZrsn  [km]	=rsE [km];
      EZd1rsn[km]	=rsE1[km];
      EZd2rsn[km]	=rsE2[km];
      km++;
    }
  }
  if(ESiAx != 0){
    ECRestore2DGeomX();
  }
#ifdef XWIN
  {
    extern char *CbUserWarning,ESMessage[];
    sprintf(ESMessage,"Solution not found\nGeometry is restored\n");
    CbUserWarning	=ESMessage;
  }
#endif
  return(0);
}

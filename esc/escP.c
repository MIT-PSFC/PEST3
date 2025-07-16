#include "esc_local.h"
#define aaSHMEMMEM

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "esc.h"
#ifdef PFC
#include "pfc.h"
#endif
#ifndef mzl_ESmain
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern int ESShMem;
extern unsigned int ESmemlev;
#endif

#ifndef mzl_1Dcore
extern int ESkNormJ,ESkNormP;
extern int ESNa,ESNa1,ESNp1;
extern double *ESsa,*ESpa,*ESsr,*ESaR0;
extern double ESa0,ESav;
extern int ESFlHole;
/* Bitwise Flag of modications related to change in 
   the hole size:  
   1 - interpolate or extrapolate Fourier harmonics 
   of coordinates;
   2 - change the radial coordinate;
   4 - change the xd-coordinates in the profile data;
   4 - change the F-profile;
   8 - interpolate or extrapolate js,jp-profiles;
   16 - calculate the derivatives of coordinates at the
   hole boundary;
   32 - recalculate the plasma profiles;
   */

extern double *ESg22c,*ESg22c2,*ESLc,*ESLc2,*ESVc,*ESVc2;
extern double ESBt,ESRBt,ESRext;
extern double ESgbext,ESgbN;
extern PlProf PRdC,PlPrP[],PlCrP[];
extern RadCoord RdC[];
extern double ESgaP[];
extern double ESgaC[];
static int NPrP=0,NCrP=0,NRdC=0,NPlPD=0;
static int iPIn=3,iCIn=2,iXIn=0;
static int ioPrP[ESNPRPR],ioCrP[ESNCRPR],*kaRdC;
#endif

#ifndef mzl_PROFL
extern double ESIpl;
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

#ifdef PFC
extern void EFITjpPlot();
extern void EFITjsPlot();

extern int MSEFl,nMSE;
extern DATAMSE MSE[64];
#endif

#ifdef MSE
extern int ESXNa,ESXNa1,*ESXsxK;
extern double *ESXsxD,*ESXsrD,*ESXsjD,*ESXspD,*ESXsqD,*ESXsiD;
extern double *ESXsx,*ESXsr,*ESXsj,*ESXsp,*ESXsq,*ESXsi;
extern int ESNMSE;
extern double ESrMSE[],ESgiMSE[];
#endif

extern int ESnDtPr,ESnDtCr;
extern double ESxDtPr[],ESxDtCr[];

static int iSplDP,iSplDP1,iSplDPold;
static double crxDP,cHxDP,SplDPx,SplDPxx,SplDPX,SplDPXX;

static int iSplDC,iSplDC1,iSplDCold;
static double crxDC,cHxDC,SplDCx,SplDCxx,SplDCX,SplDCXX;
#endif
extern int ESEqSolvRadC,ESEqSolvInPr,ESEqSolvInCr,ESEqSolvInRadC;
extern double *ESg22s,*ESLs,*ESVs,*ESDPs;

static int PrDataRadC=2,CrDataRadC=2;


#ifdef STUPID
double ESWonchulP(double gY);
double ESWonchulQ(double gY);

static int Njpol=0,Nppol=0,IPol[21];
static double Apol[484],Pjk[21],Ppk[21],Xjk[21],Xpk[21];
static double Xja=0.,d1Xja=0.,d2Xja=0.,Xpa=0.,d1Xpa=0.,d2Xpa=0.;

int ESData2Poly(double *x,double *y, int K)
{
  int i,k,ik;
  double pkx,xi;
  char str[3];

  ik	=0;
  for(i=0; i < K; i++){
    pkx	=1.;
    xi	=x[i]*x[i];
    for(k=0; k < K; k++){
      Apol[ik++]=pkx;
      pkx	*=xi;
    }
  }
  if(LUdcmp(Apol,K,IPol,&xi)){
    return(-1);
  }
  LUbksb(Apol,K,IPol,y);
  switch(ESEqSolvInCr){
  case 0:
    sprintf(str,"js");
    break;
  case 1:
    sprintf(str,"jb");
    break;
  case 3:
    sprintf(str," T");
    break;
  case 6:
    sprintf(str," q");
    break;
  case 7:
    sprintf(str,"gm");
    break;
  default:
    sprintf(str," Y");
    break;
  }
  for(i=0; i < K; i++){
    if(i){
      printf("   %+16.9e*gY^%d\n",y[i],i);
    }
    else{
      printf("===========================\n");
      printf("%s=%16.9e\n",str,y[i]);
    }
  }
  return(0);
}
#endif

int ESSolov2ES(double Ps,double T)
{
  int i,nd;

  ESEqSolvInPr	=1;
  nd	=PlPrP[1].nd;
  for(i=0; i < nd; i++){
    PlPrP[1].yd[i]	=Ps;
  }
  ESEqSolvInCr	=3;
  nd	=PlCrP[3].nd;
  for(i=0; i < nd; i++){
    PlCrP[3].yd[i]	=T;
  }
  for(i=0; i < ESNa1; i++){
    ESPs[i]	=Ps;
    ESPs1a[i]	=0.;
    ESPs2a[i]	=0.;
    ESaT[i]	=T;
    ESaT1a[i]	=0.;
    ESaT2a[i]	=0.;
  }
  return(0);
}

int ESGetPlPrfAddr(double **Pr,double **Cr,double **xPr,double **xCr
		   ,int **kPr,int **kCr)
{
  *Pr	=PlPrP[ESEqSolvInPr].yd;
  *xPr	=PlPrP[ESEqSolvInPr].xd;
  *kPr	=PlPrP[ESEqSolvInPr].kk;
  *Cr	=PlCrP[ESEqSolvInCr].yd;
  *xCr	=PlCrP[ESEqSolvInCr].xd;
  *kCr	=PlCrP[ESEqSolvInCr].kk;
  return(0);
}

int ESPlPrAddr2TRANSP(double **p,double **gm,double **xgm)
{
  *p	=PlPrP[2].yd;
  *gm	=PlCrP[7].yd;
  *xgm	=PlCrP[7].xd;
  return(0);
}

int ESPrSplData2Pr(PlProf *Prof)
{
  int i,i1,i2,nd1,inx;
  double dg0,dg1,x,dx,dy;
  double *xd,*d0,*d2,*x0,*y0;

  xd	=Prof->xd;
  d0	=Prof->d0;
  d2	=Prof->d2;
  nd1	=Prof->nd-1;
  x0	=Prof->x0;
  y0	=Prof->y0;
  inx	=0;
  
  i1	=0;
  for(i=0; i < ESNa1; i++){
    if(inx){
      i2	=i == 0 ? 0 : i-1;
#ifdef H
      splr1(&x,&dx,ESsa+i,RadC[inx].x,RadC[inx].d2x,ESsa,&na,&i2);
#endif
    }
    else{
      x		=ESsa[i];
      dx	=1.;
    }
    while(i1 < nd1 && xd[i1+1] < x) i1++;

    splr1(y0+i,&dy,&x,d0,d2,xd,&nd1,&i1);
    if(i == 0) dg0=dy*dx;
    x0[i]	=x;
  }
  dg1	=dy*dx;
  splAA(y0,Prof->y1,Prof->y2,&dg0,&dg1);
  return(0);
}

int ESReInitPlProf(int ND)
{
  int i,j,nd,n;
  double d,s;

  if(NPrP){
    return(0);
  }

  n	=ESNa1+2;
  nd	=ND;

  NPrP	=ESNPRPR;
  NCrP	=ESNCRPR;
  NRdC	=ESNRDCR;
  NPlPD	=ESND;

  d	=(1.-ESa0)/5.;
  for(i=0; i < NRdC; i++){
    RdC[i].x0	=(double*)malloc(6*n*sizeof(double));
    if(RdC[i].x0 == NULL){
      FailureAlarm((char*)RdC[i].x0
		   ,"ESReInitPlProf() - no memory for RdC[i].x0");
      ESexit(0);
    }
    RdC[i].x1	=RdC[i].x0+n;
    RdC[i].x2	=RdC[i].x1+n;
    RdC[i].a0	=RdC[i].x2+n;
    RdC[i].a1	=RdC[i].a0+n;
    RdC[i].a2	=RdC[i].a1+n;
    RdC[i].X	=1.;
    RdC[i].io	=0;
    for(j=0; j < ESNa1; j++){
      RdC[i].x0[j]	=ESsa[j];
      RdC[i].x1[j]	=1.;
      RdC[i].x2[j]	=0.;
      RdC[i].a0[j]	=ESsa[j];
      RdC[i].a1[j]	=1.;
      RdC[i].a2[j]	=0.;
    }
    RdC[i].nd	=6;
    for(j=0; j < NPlPD; j++){
      if(j < 6){
	RdC[i].ad[j]	=ESa0+d*j;
	RdC[i].xd[j]	=ESa0+d*j;
	RdC[i].kk[j]	=1;
      }
      else{
	RdC[i].ad[j]	=0.;
	RdC[i].xd[j]	=0.;
	RdC[i].kk[j]	=0;
      }
    }
  }
  RdC[iXIn].io	=1;
  
  strcpy(RdC[0].Nm,"b");
  ESsh		=RdC[0].x0;
  ESsh1a	=RdC[0].x1;
  ESsh2a	=RdC[0].x2;
  strcpy(RdC[1].Nm,"qaV");
  ESqV		=RdC[1].x0;
  ESqV1a	=RdC[1].x1;
  ESqV2a	=RdC[1].x2;
  strcpy(RdC[2].Nm,"qgF");
  ESqgF		=RdC[2].x0;
  ESqgF1a	=RdC[2].x1;
  ESqgF2a	=RdC[2].x2;
  strcpy(RdC[3].Nm,"qgY");
  ESqgY		=RdC[3].x0;
  ESqgY1a	=RdC[3].x1;
  ESqgY2a	=RdC[3].x2;

  kaRdC		=(int*)malloc(NPlPD*sizeof(int));
  if(kaRdC == NULL){
    FailureAlarm((char*)kaRdC,"ESReInitPlProf() - no memory for kaRdC");
    ESexit(0);
  }
  PRdC.x0	=(double*)malloc(4*n*sizeof(double));
  if(PRdC.x0 == NULL){
    FailureAlarm((char*)PRdC.x0,"ESReInitPlProf() - no memory for PRdC.x0");
    ESexit(0);
  }
  PRdC.y0	=PRdC.x0+n;
  PRdC.y1	=PRdC.y0+n;
  PRdC.y2	=PRdC.y1+n;
  for(i=0; i < ESNa1; i++){
    PRdC.x0[i]	=ESsa[i];
    PRdC.y0[i]	=ESsa[i];
    PRdC.y1[i]	=1.;
    PRdC.y2[i]	=0.;
  }
  PRdC.io	=0;
  PRdC.xin	=2;
  PRdC.nd	=6;
  for(j=0; j < NPlPD; j++){
    if(j < 6){
      PRdC.xd[j]	=ESa0+d*j;
      PRdC.yd[j]	=ESa0+d*j;
      PRdC.kk[j]	=1;
    }
    else{
      PRdC.xd[j]	=0.;
      PRdC.yd[j]	=0.;
      PRdC.kk[j]	=0;
    }
  }
  ESnDtPr=nd-1;
  ESnDtCr=nd-1;
  d	=(1.-ESa0)/ESnDtPr;
  for(j=0; j < NPlPD; j++){
    if(j < nd){
      ESxDtPr[j]	=ESa0+d*j;
      ESxDtCr[j]	=ESa0+d*j;
    }
    else{
      ESxDtPr[j]	=0.;
      ESxDtCr[j]	=0.;
    }
  }
  for(i=0; i < NPrP; i++){
    PlPrP[i].x0	=(double*)malloc(4*n*sizeof(double));
    if(PlPrP[i].x0 == NULL){
      FailureAlarm((char*)PlPrP[i].x0
		   ,"ESReInitPlProf() - no memory for PlPrP[i].x0");
      ESexit(0);
    }
    PlPrP[i].y0	=PlPrP[i].x0+n;
    PlPrP[i].y1	=PlPrP[i].y0+n;
    PlPrP[i].y2	=PlPrP[i].y1+n;
    for(j=0; j < ESNa1; j++){
      PlPrP[i].x0[j]	=ESsa[j];
      PlPrP[i].y0[j]	=0.;
      PlPrP[i].y1[j]	=0.;
      PlPrP[i].y2[j]	=0.;
    }

    PlPrP[i].io	=0;
    PlPrP[i].xin=2;
    PlPrP[i].nd	=nd;
    ioPrP[i]	=0;
    for(j=0; j < NPlPD; j++){
      if(j < nd){
	PlPrP[i].xd[j]	=ESa0+d*j;
	PlPrP[i].yd[j]	=0.;
	PlPrP[i].kk[j]	=1;
      }
      else{
	PlPrP[i].xd[j]	=0.;
	PlPrP[i].yd[j]	=0.;
	PlPrP[i].kk[j]	=0;
      }
    }
  }
  for(i=0; i < NCrP; i++){
    PlCrP[i].x0	=(double*)malloc(4*n*sizeof(double));
    if(PlCrP[i].x0 == NULL){
      FailureAlarm((char*)PlCrP[i].x0
		   ,"ESReInitPlProf() - no memory for PlCrP[i].x0");
      ESexit(0);
    }
    PlCrP[i].y0	=PlCrP[i].x0+n;
    PlCrP[i].y1	=PlCrP[i].y0+n;
    PlCrP[i].y2	=PlCrP[i].y1+n;
    for(j=0; j < ESNa1; j++){
      PlCrP[i].x0[j]	=ESsa[j];
      PlCrP[i].y0[j]	=0.;
      PlCrP[i].y1[j]	=0.;
      PlCrP[i].y2[j]	=0.;
    }
    PlCrP[i].io	=0;
    PlCrP[i].xin=2;
    PlCrP[i].nd	=nd;
    ioCrP[i]	=0;
    for(j=0; j < NPlPD; j++){
      if(j < nd){
	PlCrP[i].xd[j]	=ESa0+d*j;
	PlCrP[i].yd[j]	=0.;
	PlCrP[i].kk[j]	=1;
      }
      else{
	PlCrP[i].xd[j]	=0.;
	PlCrP[i].yd[j]	=0.;
	PlCrP[i].kk[j]	=0;
      }
    }
  }

  strcpy(PRdC.Nm,	"RC     ");
  ESRC	=PRdC.y0;
  ESRC1a=PRdC.y1;
  ESRC2a=PRdC.y2;

  strcpy(PlPrP[0].Nm,	"j_PfShl");
  ESjp	=PlPrP[0].y0;
  ESjp1a=PlPrP[0].y1;
  ESjp2a=PlPrP[0].y2;
  strcpy(PlPrP[1].Nm,	"Ps     ");
  ESPs	=PlPrP[1].y0;
  ESPs1a=PlPrP[1].y1;
  ESPs2a=PlPrP[1].y2;
  strcpy(PlPrP[2].Nm,	"p      ");
  ESsp	=PlPrP[2].y0;
  ESsp1a=PlPrP[2].y1;
  ESsp2a=PlPrP[2].y2;
  strcpy(PlPrP[3].Nm,	"p'/a   ");
  strcpy(PlPrP[4].Nm,	"ga_CHT ");
  strcpy(PlPrP[5].Nm,	"n_e    ");
  strcpy(PlPrP[6].Nm,	"n_i    ");
  strcpy(PlPrP[7].Nm,	"T_e    ");
  strcpy(PlPrP[8].Nm,	"T_i    ");
  strcpy(PlPrP[9].Nm,	"p_e    ");
  strcpy(PlPrP[10].Nm,	"p_i    ");
  strcpy(PlPrP[11].Nm,	"gb_Sh  ");
  ESgb	=PlPrP[11].y0;
  ESgb1a=PlPrP[11].y1;
  ESgb2a=PlPrP[11].y2;
  strcpy(PlPrP[12].Nm,	"gb     ");
  ESgB	=PlPrP[12].y0;
  ESgB1a=PlPrP[12].y1;
  ESgB2a=PlPrP[12].y2;

  strcpy(PlCrP[0].Nm,	"j_GSh  ");
  ESjs	=PlCrP[0].y0;
  ESjs1a=PlCrP[0].y1;
  ESjs2a=PlCrP[0].y2;

  strcpy(PlCrP[1].Nm,	"j_||   ");
  ESjb	=PlCrP[1].y0;
  ESjb1a=PlCrP[1].y1;
  ESjb2a=PlCrP[1].y2;

  strcpy(PlCrP[17].Nm,	"J      ");
  ESaJ	=PlCrP[17].y0;
  ESaJ1a=PlCrP[17].y1;
  ESaJ2a=PlCrP[17].y2;

  strcpy(PlCrP[3].Nm,	"FF'    ");
  ESaT	=PlCrP[3].y0;
  ESaT1a=PlCrP[3].y1;
  ESaT2a=PlCrP[3].y2;

  strcpy(PlCrP[4].Nm,	"<j.B>  ");
  ESjB	=PlCrP[4].y0;
  ESjB1a=PlCrP[4].y1;
  ESjB2a=PlCrP[4].y2;

  strcpy(PlCrP[5].Nm,	"Fs     ");

  strcpy(PlCrP[6].Nm,	"q      ");
  ESsq	=PlCrP[6].y0;
  ESsq1a=PlCrP[6].y1;
  ESsq2a=PlCrP[6].y2;

  strcpy(PlCrP[7].Nm,	"1/q    ");
  ESgm	=PlCrP[7].y0;
  ESgm1a=PlCrP[7].y1;
  ESgm2a=PlCrP[7].y2;

  strcpy(PlCrP[8].Nm,	"gY     ");
  ESgY	=PlCrP[8].y0;
  ESgY1a=PlCrP[8].y1;
  ESgY2a=PlCrP[8].y2;

  strcpy(PlCrP[9].Nm,	"gF    ");
  ESgF	=PlCrP[9].y0;
  ESgF1a=PlCrP[9].y1;
  ESgF2a=PlCrP[9].y2;

  strcpy(PlCrP[10].Nm,	"F      ");
  ESaF	=PlCrP[10].y0;
  ESaF1a=PlCrP[10].y1;
  ESaF2a=PlCrP[10].y2;
  strcpy(PlCrP[11].Nm,	"FF     ");
  ESFF	=PlCrP[11].y0;
  ESFF1a=PlCrP[11].y1;
  ESFF2a=PlCrP[11].y2;

  strcpy(PlCrP[12].Nm,	"ir+il  ");
  strcpy(PlCrP[13].Nm,	"ir-il  ");
  strcpy(PlCrP[14].Nm,	"j_r    ");
  strcpy(PlCrP[15].Nm,	"j_l    ");

  strcpy(PlCrP[16].Nm,	"dgY/a  ");
  ESdgY	=PlCrP[16].y0;
  ESdgY1a=PlCrP[16].y1;
  ESdgY2a=PlCrP[16].y2;

  strcpy(PlCrP[2].Nm,	"j_||R0 ");
  ESjR	=PlCrP[2].y0;
  ESjR1a=PlCrP[2].y1;
  ESjR2a=PlCrP[2].y2;

  ioPrP[0]	=1;
  PlPrP[0].io	=1;
  d	=1./(nd-1);
  for(i=0; i < nd; i++){
    s	=(d*i-ESa0)/(1.-ESa0);
    s	=1.-s*s;
    PlPrP[0].yd[i]	=0.2*s;
    PlPrP[2].yd[i]	=0.001*s*s;
  }
  ESPrData2Spl(PlPrP,PlPrP[0].kk,PlPrP[0].xd,PlPrP[0].yd,ESgaP+1,ESgaP+3);
  ESPrSplData2Pr(PlPrP);
  ESPrData2Spl(PlPrP+2,PlPrP[2].kk,PlPrP[2].xd,PlPrP[2].yd,ESgaP+1,ESgaP+3);
  ESPrSplData2Pr(PlPrP+2);

  ioCrP[0]	=1;
  PlCrP[0].io	=1;
  for(i=0; i < nd; i++){
    s	=(d*i-ESa0)/(1.-ESa0);
    s	=1.-s*s;
    PlCrP[0].yd[i]	=s;
    PlCrP[1].yd[i]	=s;
  }
  ESPrData2Spl(PlCrP,PlCrP[0].kk,PlCrP[0].xd,PlCrP[0].yd,ESgaC+1,ESgaC+3);
  ESPrSplData2Pr(PlCrP);
  ESPrData2Spl(PlCrP+1,PlCrP[1].kk,PlCrP[1].xd,PlCrP[1].yd,ESgaC+1,ESgaC+3);
  ESPrSplData2Pr(PlCrP+1);

#ifdef MSE
  ESXsx	=(double*)malloc(12*NPlPD*sizeof(double));
  if(ESXsx == NULL){
    FailureAlarm((char*)ESXsx,"ESReInitPlProf() - no memory for ESXsx");
    ESexit(0);	
  }
  ESXsr	=ESXsx	+NPlPD;
  ESXsj	=ESXsr	+NPlPD;
  ESXsp	=ESXsj	+NPlPD;
  ESXsq	=ESXsp	+NPlPD;
  ESXsi	=ESXsq	+NPlPD;
  ESXsxD	=ESXsi	+NPlPD;
  ESXsrD	=ESXsxD	+NPlPD;
  ESXsjD	=ESXsrD	+NPlPD;
  ESXspD	=ESXsjD	+NPlPD;
  ESXsqD	=ESXspD	+NPlPD;
  ESXsiD	=ESXsqD	+NPlPD;
  ESXsxK	=(int*)malloc(NPlPD*sizeof(int));

  d	=2./ESXNa;
  for(i=0; i < NPlPD; i++){
    ESXsx[i]	=i < ESXNa1 ?  -1.+d*i : 0.;
    ESXsr[i]	=ESXsx[i];
    ESXsj[i]	=0.;
    ESXsp[i]	=0.;
    ESXsq[i]	=0.;
    ESXsi[i]	=0.;
    ESXsxK[i]	=0;	
    ESXsxD[i]	=0.;
    ESXsrD[i]	=0.;
    ESXsjD[i]	=0.;
    ESXspD[i]	=0.;
    ESXsqD[i]	=0.;
    ESXsiD[i]	=0.;
  }
#endif
  ESmemlev |=0x01000000;
  return(0);
}

int ESDeInitPlProf()
{
  int i;

  if(ESNa1 == 0)
    return(0);
#ifdef MSE
  free(ESXsxK);
  free(ESXsx);
#endif
  for(i=0; i < NCrP; i++){
    free(PlCrP[i].x0);
  }
  for(i=0; i < NPrP; i++){
    free(PlPrP[i].x0);
  }
  free(PRdC.x0);
  free(kaRdC);
  for(i=0; i < NRdC; i++){
    free(RdC[i].x0);
  }
  NPrP	=0;
  NCrP	=0;
  NRdC	=0;
  return(0);
}

#ifdef STUPID
double SplDPdA=1.,SplDPd2A=0.;
double SplDCdA=1.,SplDCd2A=0.;
static double dFDat[257];
#endif


int ESFirstSetSplDPr()
{
  iSplDP	=0;
  iSplDP1	=1;
  crxDP		=ESxDtPr[1]-ESxDtPr[0];
  cHxDP		=EZcr6*crxDP*crxDP;
  crxDP		=1./crxDP;
  iSplDPold	=0;
  return(0);
}

int ESSetSplDPr(double A)
{
#ifdef STUPID
  if(PrDataRadC != ESEqSolvRadC){
    ESSetSplA(A);
    splRA2(&A,&SplDPdA,&SplDPd2A,RdC[PrDataRadC].x0,RdC[PrDataRadC].x2);
  }
#endif
  if(A > ESxDtPr[ESnDtPr]){
    A	=ESxDtPr[ESnDtPr];
  }
  if(A < ESxDtPr[0]){
    A	=ESxDtPr[0];
  }
  while(iSplDP < ESnDtPr && ESxDtPr[iSplDP+1] < A) iSplDP++;
  while(iSplDP > 0 && ESxDtPr[iSplDP] > A) iSplDP--;
  if(iSplDPold	!= iSplDP){
    iSplDPold	=iSplDP;
    iSplDP1	=iSplDP+1;
    crxDP	=ESxDtPr[iSplDP1]-ESxDtPr[iSplDP];
    cHxDP	=EZcr6*crxDP*crxDP;
    crxDP	=1./crxDP;
  }
  SplDPx	=(A-ESxDtPr[iSplDP])*crxDP;
  SplDPxx	=SplDPx*SplDPx;

  SplDPX	=1.-SplDPx;
  SplDPXX	=SplDPX*SplDPX;
  return(0);
}

int splRDPr(double*f, double*df,int k)
{
  double *d0,*d2;
  d0	=PlPrP[k].d0;
  d2	=PlPrP[k].d2;

  *f	=SplDPX*d0[iSplDP]+SplDPx*d0[iSplDP1]
    +(SplDPX*(SplDPXX-1.)*d2[iSplDP]+SplDPx*(SplDPxx-1.)*d2[iSplDP1])*cHxDP;
  if(df != NULL){
    *df	=(d0[iSplDP1]-d0[iSplDP]+
	  ((3.*SplDPxx-1.)*d2[iSplDP1]-(3.*SplDPXX-1.)*d2[iSplDP])*cHxDP)
      *crxDP;
  }

#ifdef STUPID
  if(PrDataRadC != ESEqSolvRadC && df != NULL){
    *df	*=SplDPdA;
  }
#endif
  return(0);
}

int splRDPr2(double*f, double*df, double*d2f, int k)
{
  double dy,*d0,*d2;

  d0	=PlPrP[k].d0;
  d2	=PlPrP[k].d2;
  *f	=SplDPX*d0[iSplDP]+SplDPx*d0[iSplDP1]
    +(SplDPX*(SplDPXX-1.)*d2[iSplDP]+SplDPx*(SplDPxx-1.)*d2[iSplDP1])*cHxDP;
  if(df != NULL){
    *df	=(d0[iSplDP1]-d0[iSplDP]
	  +((3.*SplDPxx-1.)*d2[iSplDP1]-(3.*SplDPXX-1.)*d2[iSplDP])*cHxDP)
      *crxDP;
  }
  if(d2f != NULL){
    *d2f	=SplDPX*d2[iSplDP]+SplDPx*d2[iSplDP1];
  }
#ifdef STUPID
  if(PrDataRadC != ESEqSolvRadC){
    if(d2f != NULL){
      if(df != NULL){
	dy	=*df;
      }
      else{
	dy =(d0[iSplDP1]-d0[iSplDP]
	     +((3.*SplDPxx-1.)*d2[iSplDP1]-(3.*SplDPXX-1.)*d2[iSplDP])*cHxDP)
	  *crxDP;
      }
      *d2f	=(*d2f)*SplDPdA*SplDPdA+dy*SplDPd2A;
    }
    if(df != NULL){
      *df	*=SplDPdA;
    }
  }
#endif
  return(0);
}

int ESFirstSetSplDCr()
{
  iSplDC	=0;
  iSplDC1	=1;
  crxDC		=ESxDtCr[1]-ESxDtCr[0];
  cHxDC		=EZcr6*crxDC*crxDC;
  crxDC		=1./crxDC;
  iSplDCold	=0;
  return(0);
}

int ESSetSplDCr(double A)
{
#ifdef STUPID
  if(CrDataRadC != ESEqSolvRadC){
    ESSetSplA(A);
    splRA2(&A,&SplDCdA,&SplDCd2A,RdC[CrDataRadC].x0,RdC[CrDataRadC].x2);
  }
  if(Njpol){
    ESSetSplA(A);
    splRA2(&Xja,&d1Xja,&d2Xja,ESqgY,ESqgY2a);
    return(0);
  }
#endif
  if(A > ESxDtCr[ESnDtCr]) A=ESxDtCr[ESnDtCr];
  if(A < ESxDtCr[0]) A=ESxDtCr[0];
  while(iSplDC < ESnDtCr && ESxDtCr[iSplDC+1] < A) iSplDC++;
  while(iSplDC > 0 && ESxDtCr[iSplDC] > A) iSplDC--;
  if(iSplDCold	!= iSplDC){
    iSplDCold	=iSplDC;
    iSplDC1	=iSplDC+1;
    crxDC	=ESxDtCr[iSplDC1]-ESxDtCr[iSplDC];
    cHxDC	=EZcr6*crxDC*crxDC;
    crxDC	=1./crxDC;
  }
  SplDCx	=(A-ESxDtCr[iSplDC])*crxDC;
  SplDCxx	=SplDCx*SplDCx;

  SplDCX	=1.-SplDCx;
  SplDCXX	=SplDCX*SplDCX;
  return(0);
}

int splRDCr(double*f,double*df,int k)
{
  double *d0,*d2;

#ifdef STUPID
  if(Njpol){
    int i;
    double p2ix,px;

    *f	=Pjk[0];
    p2ix=1.;
    px	=Xja*Xja;
    for(i=1; i < Njpol; i++){
      p2ix	*=px;
      *f	+=Pjk[i]*p2ix;
    }
    if(df != NULL){
      *df	=0.;
      p2ix=1.;
      for(i=1; i < Njpol; i++){
	*df	+=Pjk[i]*i*p2ix;
      }
      *df	*=2.*Xja*d1Xja;
    }
    if(ESEqSolvInCr == 6){
      *f	=1./(*f);
      if(df != NULL){
	*df	*=-(*f)*(*f);
      }
    }
    return(0);
  }
#endif
  d0	=PlCrP[k].d0;
  d2	=PlCrP[k].d2;
  *f	=SplDCX*d0[iSplDC]+SplDCx*d0[iSplDC1]
    +(SplDCX*(SplDCXX-1.)*d2[iSplDC]+SplDCx*(SplDCxx-1.)*d2[iSplDC1])*cHxDC;
  if(df != NULL){
    *df	=(d0[iSplDC1]-d0[iSplDC]+
	  ((3.*SplDCxx-1.)*d2[iSplDC1]-(3.*SplDCXX-1.)*d2[iSplDC])*cHxDC)
      *crxDC;
  }
#ifdef STUPID
  if(CrDataRadC != ESEqSolvRadC && df != NULL){
    *df	*=SplDCdA;
  }
#endif
  return(0);
}

int splRDCr2(double*f, double*df, double*d2f, int k)
{
  double *d0,*d2,dy,d2y;

#ifdef STUPID
  if(Njpol){
    int i;
    double p2ix,px;

    *f	=Pjk[0];
    dy	=0.;
    p2ix	=1.;
    px	=Xja*Xja;
    for(i=1; i < Njpol; i++){
      dy	+=Pjk[i]*i*p2ix;
      p2ix	*=px;
      *f	+=Pjk[i]*p2ix;
    }
    d2y		=0.;
    p2ix	=1.;
    for(i=2; i < Njpol; i++){
      d2y	+=Pjk[i]*i*(i-1)*p2ix;
      p2ix	*=px;
    }
    p2ix	=2.*Xja*d1Xja;
    d2y		=d2y*p2ix*p2ix+2.*dy*(Xja*d2Xja+d1Xja*d1Xja);
    dy		*=p2ix;
    if(ESEqSolvInCr == 6){
      *f	=1./(*f);
      px	=dy*(*f);
      if(df != NULL){
	*df	=-px*(*f);
      }
      if(d2f != NULL){
	*d2f	=(2.*px*px-d2y*(*f))*(*f);
      }
    }
    return(0);
  }
#endif
  d0	=PlCrP[k].d0;
  d2	=PlCrP[k].d2;
  *f	=SplDCX*d0[iSplDC]+SplDCx*d0[iSplDC1]
    +(SplDCX*(SplDCXX-1.)*d2[iSplDC]+SplDCx*(SplDCxx-1.)*d2[iSplDC1])*cHxDC;
  if(df != NULL){
    *df	=(d0[iSplDC1]-d0[iSplDC]
	  +((3.*SplDCxx-1.)*d2[iSplDC1]-(3.*SplDCXX-1.)*d2[iSplDC])*cHxDC)
      *crxDC;
  }
  if(d2f != NULL){
    *d2f	=SplDCX*d2[iSplDC]+SplDCx*d2[iSplDC1];
  }
#ifdef STUPID
  if(CrDataRadC != ESEqSolvRadC){
    if(d2f != NULL){
      if(df != NULL){
	dy	=*df;
      }
      else{
	dy=(d0[iSplDC1]-d0[iSplDC]
	     +((3.*SplDCxx-1.)*d2[iSplDC1]-(3.*SplDCXX-1.)*d2[iSplDC])*cHxDC)
	  *crxDC;
      }
      *d2f	=(*d2f)*SplDCdA*SplDCdA+dy*SplDCd2A;
    }
    if(df != NULL){
      *df	*=SplDCdA;
    }
  }
#endif
  return(0);
}

int ESPrSplData2PrPr(int k)
{
  int i,j;
  static double dy0=0.;
  double *x0,*y0,*y1,*y2;

  x0	=PlPrP[k].x0;
  y0	=PlPrP[k].y0;
  y1	=PlPrP[k].y1;
  y2	=PlPrP[k].y2;
  for(i=0; i < ESNa1; i++){
    x0[i]	=ESsa[i];
    ESSetSplDPr(x0[i]);
    splRDPr2(y0+i,y1+i,y2+i,k);
  }
  return(0);
}

int ESScalePrData(double s)
{
  int i;
  PlProf *lP;

  lP	=PlPrP+ESEqSolvInPr;
  for(i=0; i < ESnDtPr; i++){
    lP->yd[i]	*=s;
    lP->d0[i]	*=s;
    lP->d2[i]	*=s;
  }
  return(0);
}

int ESScaleCrData(double s)
{
  int i;
  PlProf *lP;

  lP	=PlCrP+ESEqSolvInCr;
  for(i=0; i < ESnDtCr; i++){
    lP->yd[i]	*=s;
    lP->d0[i]	*=s;
    lP->d2[i]	*=s;
  }
  return(0);
}

int ESScalePrCrData(double sP,double sJ)
{
  int i;
  PlProf *lP;

  lP	=PlPrP+ESEqSolvInPr;
  for(i=0; i < ESnDtPr; i++){
    lP->yd[i]	*=sP;
    lP->d0[i]	*=sP;
    lP->d2[i]	*=sP;
  }
  lP	=PlCrP+ESEqSolvInCr;
  for(i=0; i < ESnDtCr; i++){
    lP->yd[i]	*=sJ;
    lP->d0[i]	*=sJ;
    lP->d2[i]	*=sJ;
  }
  return(0);
}

int ESPrSplData2CrPr(int k)
{
  int i;
  static double dy0=0.;
  double dy1,*x0,*y0,*y1,*y2;

  x0	=PlCrP[k].x0;
  y0	=PlCrP[k].y0;
  y1	=PlCrP[k].y1;
  y2	=PlCrP[k].y2;
  for(i=0; i < ESNa1; i++){
    x0[i]	=ESsa[i];
    ESSetSplDCr(x0[i]);
    splRDCr2(y0+i,y1+i,y2+i ,k);
  }
  return(0);
}

int ESData2Order(double *xd,double *yd,int *k2k,int nd)
{
  int i,j,k;
  double s;
  
  for(k=0; k < nd; k++){
    if(k2k[k] == 0){
      i	=k+1;
      while(i < nd && k2k[i] == 0) i++;
      if(i == nd){
	return(k);
      }
      else{
	s	=xd[i];
	xd[i]	=xd[k];
	xd[k]	=s;
	s	=yd[i];
	yd[i]	=yd[k];
	yd[k]	=s;
	k2k[k]	=1;
	k2k[i]	=0;
      }
    }
    i	=k+1;
    while(i < nd){
      if(k2k[i]){
	if(xd[i] < xd[k]){
	  s	=xd[i];
	  xd[i]	=xd[k];
	  xd[k]	=s;
	  s	=yd[i];
	  yd[i]	=yd[k];
	  yd[k]	=s;
	}
	else{
	  if(xd[i] == xd[k]){
	    k2k[i]	=0;
	  }
	}
      }
      i++;
    }
  }
  return(nd);
}

int ESPrData2Spl(PlProf *Prof,int *k2k,double *xd,double *yd,
	     double *EZga1,double *EZga3)
{
  extern double *ESg22c1,*ESLc1;
  int i,j,k;
  double dg0,dg1,s;

  k	=0;
  j	=0;
  s	=-1e+20;
  while(j < NPlPD){
    if(k2k[j]){
      Prof->xd[k]	=xd[j];
      Prof->yd[k]	=yd[j];
      if(s < fabs(yd[j])) s=fabs(yd[j]);
      k++;
    }
    j++;
  }
  Prof->nd	=k;
  k--;

#ifdef H
  if(Prof->xd[0] != ESa0) Prof->xd[0]	=ESa0;
  if(Prof->xd[k] != 1.) Prof->xd[k]	=1.;
#endif
  
  if(fabs(Prof->yd[0]) < 1e-4*s || ESa0 != 0.){
    if(ESa0 == 0.) Prof->yd[0]=0.;
    EZf2spl(Prof->d0,Prof->d2,NULL,NULL,Prof->yd,NULL,Prof->xd,&k,
	  ESgaC,EZga1,ESgaP+2,EZga3);
#ifdef H
    spl(Prof->d0,Prof->d2,NULL,NULL,Prof->xd,&k);
#endif
  }
  else{
    dg0	=0.;
    EZf2spl(Prof->d0,Prof->d2,&dg0,NULL,Prof->yd,NULL,Prof->xd,&k,
	  ESgaC,EZga1,ESgaP+2,EZga3);
#ifdef H
    if(Prof != PlCrP+7){
      EZf2spl(Prof->d0,Prof->d2,&dg0,NULL,Prof->yd,NULL,Prof->xd,&k,
	    ESgaC,EZga1,ESgaP+2,EZga3);
    }
    else{/*gm*/
      dg1	=ESg22c[ESNa]*ESLc[ESNa];
      dg1=(ESjb[ESNa]*ESLc[ESNa]*ESaR0[0]/ESaF[ESNa]
	   -Prof->yd[k]*(2.*dg1+ESg22c1[ESNa]*ESLc[ESNa]
			 +ESg22c[ESNa]*ESLc1[ESNa]))/dg1;
      EZf2spl(Prof->d0,Prof->d2,&dg0,&dg1,Prof->yd,NULL,Prof->xd,&k,
	    ESgaC,EZga1,ESgaP+2,EZga3);
    }
#endif
#ifdef H
    for(i=0; i < k+1; i++){
      Prof->d0[i]	=Prof->yd[i];
    }
    spl(Prof->d0,Prof->d2,&dg0,NULL,Prof->xd,&k);
#endif
  }
  return(0);
}

int ESPlProf2Data()
{
  int i,k,nd,xin;
  double *xd,x,a,da;

  nd	=PlPrP[ESEqSolvInPr].nd;
  xd	=PlPrP[ESEqSolvInPr].xd;
  xin	=PlPrP[ESEqSolvInPr].xin;

  for(i=0; i < nd; i++){
    if(xin == ESEqSolvRadC){
      x	=xd[i];
    }
    else{
    }
    x	=xd[i];
    ESSetSplA(x);
    for(k=0; k < NPrP; k++){
      if(k != ESEqSolvInPr){
	PlPrP[k].xd[i]	=x;
	PlPrP[k].kk[i]	=1;
	if(PlPrP[k].kk[i]){
	  splRA(PlPrP[k].yd+i,NULL,PlPrP[k].y0,PlPrP[k].y2);
	}
      }
    }
  }
  while(i < NPlPD){
    for(k=0; k < NPrP; k++){
      if(k != ESEqSolvInPr) PlPrP[k].kk[i] =0;
    }
    i++;
  }
  
  nd	=PlCrP[ESEqSolvInCr].nd;
  xd	=PlCrP[ESEqSolvInCr].xd;
  for(i=0; i < nd; i++){
    x	=xd[i];
    ESSetSplA(x);
    for(k=0; k < NCrP; k++){
      if(k != ESEqSolvInCr){
	PlCrP[k].xd[i]=x;
	PlCrP[k].kk[i]=1;
	if(PlCrP[k].kk[i]){
	  splRA(PlCrP[k].yd+i,NULL,PlCrP[k].y0,PlCrP[k].y2);
	}
      }
    }
  }
  while(i < NPlPD){
    for(k=0; k < NPrP; k++){
      if(k != ESEqSolvInCr){
	PlCrP[k].kk[i]=0;
      }
    }
    i++;
  }
  return(0);
}

int ESHoleReMapInputProfData()
{
  int i;
  double x0,s;
  double *xd;

  xd	=PlPrP[ESEqSolvInPr].xd;
  x0	=xd[0];
  s	=(ESsa[ESNa]-ESa0)/(ESsa[ESNa]-x0);
  for(i=0; i < NPlPD; i++){
    xd[i]	=ESa0+(xd[i]-x0)*s;
  }
  xd	=PlCrP[ESEqSolvInCr].xd;
  x0	=xd[0];
  s	=(ESsa[ESNa]-ESa0)/(ESsa[ESNa]-x0);
  for(i=0; i < NPlPD; i++){
    xd[i]	=ESa0+(xd[i]-x0)*s;
  }
  return(0);
}

int ESHoleReMapRadCData()
{
  int i,k;
  double x0,s;
  double *xd;

  x0	=RdC[ESEqSolvRadC].xd[0];
  s	=(ESsa[ESNa]-ESa0)/(ESsa[ESNa]-x0);

  for(k=0; k < NRdC; k++){
    xd	=RdC[k].xd;
    for(i=0; i < ESNa1; i++){
      RdC[k].x0[i]	=ESsa[i];
      RdC[k].x1[i]	=1.;
      RdC[k].x2[i]	=0.;
      RdC[k].a0[i]	=ESsa[i];
      RdC[k].a1[i]	=1.;
      RdC[k].a2[i]	=0.;
    }
    for(i=0; i < NPlPD; i++){
      xd[i]	=ESa0+(xd[i]-x0)*s;
      RdC[k].ad[i]	=xd[i];
    }
  }
  return(0);
}

int ESRdCSplData2RdC(RadCoord *pRdC)
{
  int i,i1,i2,nd1;
  double dx0,x,dx;
  double *xd,*d0,*d2,*x0;

  nd1	=PRdC.nd-1;
  xd	=PRdC.xd;
  d0	=PRdC.d0;
  d2	=PRdC.d2;
  x0	=pRdC->x0;
  i1	=0;
  for(i=0; i < ESNa1; i++){
    x		=ESsa[i];
    while(i1 < nd1 && xd[i1+1] < x){
      i1++;
    }
    splr1(x0+i,&dx,&x,d0,d2,xd,&nd1,&i1);
    if(i == 0){
      dx0	=dx;
    }
  }
  splAA(x0,pRdC->x1,pRdC->x2,&dx0,&dx);
  return(0);
}

int ESRdCSpl2Data()
{
  int i,k,nd,*kk;
  double *ad;

  nd	=RdC[ESEqSolvRadC].nd;
  ad	=RdC[ESEqSolvRadC].ad;
  kk	=RdC[ESEqSolvRadC].kk;
  for(i=0; i < NPlPD; i++){
    ESSetSplA(ad[i]);
    for(k=0; k < NRdC; k++){
      if(1 || k != ESEqSolvRadC){
	RdC[k].ad[i]=ad[i];
	if(kk[i] == 1){
	  RdC[k].kk[i]=kk[i];
	}
	if(RdC[k].kk[i]){
	  splRA(RdC[k].xd+i,NULL,RdC[k].x0,RdC[k].x2);
	}
      }
    }
  }
  return(0);
}

int ESSwitchInpRadC()
{
  RadCoord *pRdC;

  pRdC		=RdC+ESEqSolvRadC;
  pRdC->nd	=ESData2Order(pRdC->ad,pRdC->xd,pRdC->kk,NPlPD);
  ESRdCSpl2Data();
  return(0);
}

int ConvertRadC(int xin2,int xin1,double *xd,int nd)
{
  int i;

  if(xin1 == ESEqSolvRadC){
    for(i=0; i < nd; i++){
      ESSetSplA(xd[i]);
      splRA(xd+i,NULL,RdC[xin2].x0,RdC[xin2].x2);
    }
    return(0);
  }
  {
    double x,x1,dx;
    for(i=0; i < nd; i++){
      x	=xd[i];
      do{
	ESSetSplA(x);
	splRA(&x1,&dx,RdC[xin1].x0,RdC[xin1].x2);
	dx	=(xd[i]-x1)/dx;
	x	+=dx;
      }while(fabs(dx) > 1e-6);      
      ESSetSplA(x);
      splRA(xd+i,NULL,RdC[xin2].x0,RdC[xin2].x2);
    }
  }
  return(0);
}

int GenerateRadC(double *xd2[ESNRDCR],double *xd,int xin,int nd)
{
  int i,k;

  if(xin == ESEqSolvRadC){
    for(i=0; i < nd; i++){
      ESSetSplA(xd[i]);
      for(k=0; k < NRdC; k++){
	if(k != xin){
	  splRA(xd2[k]+i,NULL,RdC[k].x0,RdC[k].x2);
	}
	else{
	  xd2[k][i]	=xd[i];
	}
      }
    }
    return(0);
  }
  {
    double x,x1,dx;
    for(i=0; i < nd; i++){
      x	=xd[i];
      do{
	ESSetSplA(x);
	splRA(&x1,&dx,RdC[xin].x0,RdC[xin].x2);
	dx	=(xd[i]-x1)/dx;
	x	+=dx;
      }while(fabs(dx) > 1e-6);      
      ESSetSplA(x);
      for(k=0; k < NRdC; k++){
	if(k != xin){
	  splRA(xd2[k]+i,NULL,RdC[k].x0,RdC[k].x2);
	}
	else{
	  xd2[k][i]	=xd[i];
	}
      }
    }
  }
  return(0);
}

int ESInpRadC()
{
  int i,k,nd;
  RadCoord *pRdC;
  double *ad;

  kaRdC[0]	=ESEqSolvRadC;
  pRdC		=RdC+ESEqSolvRadC;
  ad		=pRdC->ad;
#ifndef Tbl_RadC
  pRdC->nd	=ESData2Order(pRdC->ad,pRdC->xd,pRdC->kk,NPlPD);
#ifdef Later
  ESPrData2Spl(&PRdC,pRdC->kk,pRdC->ad,pRdC->xd,ESgaC+1,ESgaC+3);
  ESRdCSplData2RdC(pRdC);
#endif
  ESRdCSpl2Data();
  if(kaRdC[0] != ESEqSolvRadC){
    ESEqSolvRadC	=kaRdC[0];
    pRdC	=RdC+ESEqSolvRadC;
    ad		=pRdC->ad;
  }
  for(i=1; i < NPlPD; i++){
    kaRdC[i]	=kaRdC[0];
  }
#endif/*Tbl_RadC*/
  return(0);
}

int ESInputPrData2Pr()
{
  PlProf *pP;

  pP	=PlPrP+ESEqSolvInPr;
  ESPrData2Spl(pP,pP->kk,pP->xd,pP->yd,ESgaP+1,ESgaP+3);
  pP	=PlCrP+ESEqSolvInCr;
  ESPrData2Spl(pP,pP->kk,pP->xd,pP->yd,ESgaC+1,ESgaC+3);

  ESFirstSetSplDPr();
  ESPrSplData2PrPr(ESEqSolvInPr);
  ESFirstSetSplDCr();
  ESPrSplData2CrPr(ESEqSolvInCr);
  return(0);
}

int ESSavePrData(FILE *lf)
{
  int i;
  double y[3];
  
  fprintf(lf,"%3d %s - number of radial points and the name of profile\n"
	  ,ESNa1,PlPrP[ESEqSolvInPr].Nm);
  for(i=0; i < ESNa1; i++){
    ESSetSplDPr(ESsa[i]);
    splRDPr2(y,y+1,y+2,ESEqSolvInPr);
#ifdef H
    fprintf(lf,"%3d %6.4f %14.7e %14.7e %14.7e\n",i,ESsa[i],y[0],y[1],y[2]);
#endif
    fprintf(lf,"%3d %6.4f %14.7e\n",i,ESsa[i],y[0]);
  }
  return(0);
}

int ESReadPrData(FILE *lf)
{
  char ch;
  int i,k,n,L;
  PlProf *lP;

  lP	=PlPrP+ESEqSolvInPr;
  L	=1;
  if(fscanf(lf,"%d",&n) < 1){
    printf("PlPr.ai:%d: No data on number of points%c\n",L,'\a');
    return(1);
  }
  while((ch=getc(lf)) != '\n'){
    ;
  }
  L++;

  lP->nd	=n;
  lP->xin	=ESEqSolvInRadC;
  for(i=0; i < n; i++){
    if(fscanf(lf,"%d %lg %lg",&k,lP->xd+i,lP->yd+i) < 3){
      printf("PlPr.ai:%d: No x[%d], Pr[%d]%c\n",L,i,i,'\a');
      return(1);
    }
    L++;
    ESxDtPr[i]	=lP->xd[i];
    lP->kk[i]	=1;
  }
  ESnDtPr	=n-1;
  ESFirstSetSplDPr();
  return(0);
}

#ifdef SHMEM
static int kAstR=1,nAst;
static double AstIpl,AstgF;
static double *AstA,*AstP,*AstJ;
int ESReadPlPrShM()
{
  static int n=0;
  extern int *ESShMemI;
  extern double *ESShMemP;
  int i;
  PlProf *lP,*lJ;
  double s;

  n	=ESShMemI[6];
  AstA	=ESShMemP;
  AstP	=AstA+n;
  AstJ	=AstP+n;
  kAstR	=ESShMemI[4]/1000 ? 1 : 0;
  if(n > ESND) n=ESND;
  nAst	=n;
  s	=1./AstA[n-1];
  for(i=0; i < n; i++){
    ESxDtPr[i]	=AstA[i]*s;
  }
  ESEqSolvInRadC=(ESShMemI[4]%1000)/100;
  ESEqSolvRadC=(ESShMemI[4]%1000)/100;
  CrDataRadC=ESEqSolvInRadC;
  PrDataRadC=ESEqSolvInRadC;

  ESEqSolvInPr	=(ESShMemI[4]%100)/10;
  ESEqSolvInCr	=ESShMemI[4]%10;

  lP	=PlPrP+ESEqSolvInPr;
  lP->xin=PrDataRadC;
  lP->nd=n;

  lJ	=ESEqSolvInCr != 8 ? PlCrP+ESEqSolvInCr : PlCrP+7;
  lJ->xin	=CrDataRadC;
  lJ->nd	=n;
  memcpy((void *)lP->yd,(void *)AstP,n*sizeof(double));
  memcpy((void *)lJ->yd,(void *)AstJ,n*sizeof(double));

  if(ESxDtPr[0] != 0.){
    double pa0;
    pa0	=ESxDtPr[0]*ESxDtPr[0];
    pa0	/=(ESxDtPr[1]*ESxDtPr[1]-pa0);
    lP->yd[0]	-=(lP->yd[1]-lP->yd[0])*pa0;
    lJ->yd[0]	-=(lJ->yd[1]-lJ->yd[0])*pa0;
    ESxDtPr[0]	=0.;
  }
  memcpy((void *)lP->xd,(void *)ESxDtPr,n*sizeof(double));
  memcpy((void *)lJ->xd,(void *)ESxDtPr,n*sizeof(double));
  memcpy((void *)ESxDtCr,(void *)ESxDtPr,n*sizeof(double));
  for(i=0; i < n; i++){
    lP->kk[i]=1;
    lJ->kk[i]=1;
  }
  for(; i < ESND; i++){
    lP->kk[i]=0;
    lJ->kk[i]=0;
  }
  ESnDtPr	=n-1;
  ESnDtCr	=n-1;

  if(ESEqSolvInCr == 8){
    double dg0,dg1,Ipl;
    dg0	=0.;
    ESEqSolvInCr=7;
    i	=ESnDtCr-1;
    dg1	=(lJ->yd[ESnDtCr]-lJ->yd[i])/(lJ->xd[ESnDtCr]-lJ->xd[i]);
    dg1	=AstJ[n]*AstA[n-1];
    ESjb[ESNa]	=AstJ[n+1];
    AstIpl	=EZcrgm0*ESg22c[ESNa]*dg1;
    AstgF	=AstA[ESnDtCr]*AstA[ESnDtCr]*ESBt/2.;
    EZf2spl(lJ->d0,lJ->d2,&dg0,&dg1,lJ->yd,NULL,lJ->xd,&ESnDtCr,
	  ESgaC,ESgaC+1,ESgaC+2,ESgaC+3);
    ESFirstSetSplDCr();

    s	=EZcr2gp/(ESBt*AstA[ESnDtCr]*AstA[ESnDtCr]);
    ESSetSplDCr(ESxDtCr[0]);
    splRDCr2(&dg0,NULL,lJ->yd,ESEqSolvInCr);
    lJ->yd[0]	*=s;
    AstJ[0]	=lJ->yd[0];
    i	=0;
    for(i=1; i < n; i++){
      ESSetSplDCr(ESxDtCr[i]);
      splRDCr(&dg0,lJ->yd+i,ESEqSolvInCr);
      lJ->yd[i]	*=s/ESxDtCr[i];
      AstJ[i]	=lJ->yd[i];
    }
  }
#ifdef H  
  ESPrData2Spl(lP,lP->kk,lP->xd,lP->yd,ESgaP+1,ESgaP+3);
  ESFirstSetSplDPr();
  ESPrSplData2PrPr(ESEqSolvInPr);
  ESPrData2Spl(lJ,lJ->kk,lJ->xd,lJ->yd,ESgaC+1,ESgaC+3);
  ESFirstSetSplDCr();
  ESPrSplData2CrPr(ESEqSolvInCr);
#endif

  return(0);
}
#endif/*SHMEM*/

int ESSaveCrData(FILE *lf)
{
  int i;
  double y[3];
  
  fprintf(lf,"%3d %s\n",ESNa1,PlCrP[ESEqSolvInCr].Nm);
  for(i=0; i < ESNa1; i++){
    ESSetSplDCr(ESsa[i]);
    if(ESEqSolvInCr != 6){
      splRDCr2(y,y+1,y+2,ESEqSolvInCr);
    }
    else{
      splRDCr2(y,y+1,y+2,7);
      y[0]	=1./y[0];
    }
#ifdef H
    fprintf(lf,"%3d %6.4f %14.7e %14.7e %14.7e\n",i,ESsa[i],y[0],y[1],y[2]);
#endif
    fprintf(lf,"%3d %6.4f %14.7e\n",i,ESsa[i],y[0]);
  }
  return(0);
}

int ESReadCrData(FILE *lf)
{
  char ch;
  int i,k,n,L;
  PlProf *lP;
  double s;
  char *lc,ln[128];


  lP	=PlCrP+ESEqSolvInCr;
  L	=1;
  if(fscanf(lf,"%d",&n) < 1){
    printf("PlCr.ai:%d: No data on number of points%c\n",L,'\a');
    return(1);
  }
  while((ch=getc(lf)) != '\n') ;

  L++;
  lP->nd	=n;
  lP->xin	=ESEqSolvInRadC;
  for(i=0; i < n; i++){
    lc	=ln;
    while((*lc=getc(lf)) != '\n') lc++;
    *lc	='\0';
    if(sscanf(ln,"%d %lg %lg",&k,lP->xd+i,lP->yd+i) < 3){
      printf("PlCr.ai:%d: No x[%d], Cr[%d]%c\n",L,i,i,'\a');
      return(1);
    }
    L++;
#ifdef H
    printf("F[%2d]=%12.5e %12.5e %12.5e\n"
	   ,i,lP->xd[i],lP->yd[i]*ESRBt,lP->yd[i]*lP->yd[i]*ESRBt*ESRBt);
    lP->yd[i]	*=EZcr2*lP->yd[i]*ESRBt*ESRBt;
    lP->d0[i]	=lP->yd[i];
#endif
    ESxDtCr[i]	=lP->xd[i];
    lP->kk[i]	=1;
  }
  ESnDtPr	=n-1;

#ifdef H
  spl(lP->d0,lP->d2,NULL,NULL,lP->xd,&ESnDtPr);
  for(i=0; i < n; i++){
    SetIspl(lP->xd,lP->xd[i],ESnDtPr);
    splR(&s,dFDat+i,lP->d0,lP->d2);
  }
  for(i=0; i < n; i++){
    lP->xd[i]	=sqrt(lP->xd[i]);
    ESxDtCr[i]	=lP->xd[i];
  }
#endif
  ESFirstSetSplDCr();
  return(0);
}

#ifdef PFC
int ESInpPlPr4PF()
{
  PlProf *lP;
  int i,k,n;
  
  for(i=0; i < NPrP; i++){
    PlPrP[i].io	=i == ESEqSolvInPr ? 1 : 0;
    ioPrP[i]	=PlPrP[i].io;
  }
  lP	=PlPrP+ESEqSolvInPr;
  lP->nd	=ESData2Order(lP->xd,lP->yd,lP->kk,NPlPD);
  n	=lP->nd;
  PrDataRadC	=lP->xin;
  ESnDtPr	=lP->nd;
  for(i=0; i < ESnDtPr; i++){
    ESxDtPr[i]	=lP->xd[i];
  }
  ESPrData2Spl(lP,lP->kk,lP->xd,lP->yd,ESgaP+1,ESgaP+3);
  ESnDtPr--;
  ESFirstSetSplDPr();
  ESPrSplData2PrPr(ESEqSolvInPr);
  EC1DEqSolv();
  ESPlProf2Data();
  return(0);
}
#endif

int ESInpPlPr()
{
  PlProf *lP;
  static int FlP=0,FlY=0,FlSt=0;
  int i,k,n,*kk,InP,n2,nd,xin;
  double X[256],Y[256],X1[2];
  double *xd,*x;

  static int ii[12]={0,2,1,0,0,0,0,0,0,0,3,4};
  int n1=3*ESNa1+1;
  double *ld,*lld[5],da;
  double a[n1],f0[n1],f1[n1],f2[n1],f3[n1],f4[n1];

  lld[0]	=f0;
  lld[1]	=f1;
  lld[2]	=f2;
  lld[3]	=f3;
  lld[4]	=f4;
  
  n2	=n1+2;
  if(FlSt == 0){
    double s,f0,f1;
    s	=0.;
    for(i=0; i < ESNa1; i++){
      ESaF[i]	=ESRBt;
      ESaF1a[i]	=0.;
      ESaF2a[i]	=0.;
      ESFF[i]	=ESRBt*ESRBt;
      ESFF1a[i]	=0.;
      ESFF2a[i]	=0.;
      f1	=ESLc[i];
      if(i){
	f0	=f1;
	s	+=0.25*(f0+f1)*(ESpa[i]-ESpa[i-1]);
      }
    }
    ESgF[ESNa]=ESRBt*s/ESaR0[0];
    FlSt	=1;
  }

#ifdef SHMEM
  if(ESShMem) ESReadPlPrShM();
#endif/*SHMEM*/

  InP	=ESEqSolvInPr;
  xin	=PrDataRadC;
  if(FlP%10 != 9){
    FlP	=ESEqSolvInPr+(FlP/10)*10;
  }

  if((ESFlHole&4)){
    ESHoleReMapRadCData();
    ESHoleReMapInputProfData();
    ESFlHole	&=~4;
  }
#ifndef Tbl_PlPr
  if(ESa0 != 0.){
    ESHoleEdgeBoundary();
    ES3DMetrTnsrCSHole();
  }
  if((FlP%10) == 9){
#ifdef STUPID
    double x;
    for(i=0; i < 21; i++){
      ESSetSplA(PlPrP[2].xd[i]);
      splRA(&x,NULL,ESqgY,ESqgY2a);
      PlPrP[2].yd[i]	=ESWonchulP(x*x);
      PlPrP[2].kk[i]	=1;
    }
    ESEqSolvInPr	=2;
#endif
  }
  else{
    ESEqSolvInPr	=FlP%10;
  }
  if(FlP/10 == 3){
#ifndef AscIPr
#undef AscIPr
#endif/*AscIPr*/
    FlP	%=10;
  }  
  if(ESEqSolvInPr < 0) ESEqSolvInPr=0;
  if(ESEqSolvInPr > 1) ESEqSolvInPr=2;
  if(ESEqSolvInPr == 2 && (ESEqSolvInCr == 0 || ESEqSolvInCr == 3)){
    printf("given p(a) is not allowed for InCr == %d%c\n",ESEqSolvInCr,'\a');
    ESEqSolvInPr	=InP;
  }
  for(i=0; i < NPrP; i++){
    PlPrP[i].io	=i == ESEqSolvInPr ? 1 : 0;
    ioPrP[i]	=PlPrP[i].io;
  }
  lP	=PlPrP+ESEqSolvInPr;

  nd	=lP->nd;
  lP->nd=ESData2Order(lP->xd,lP->yd,lP->kk,NPlPD);
  nd	=lP->nd;
#ifdef SHMEM
#ifdef ASTRA
  if(ESShMem && kAstR){
    PlProf *lJ;
    double s;
    lJ	=PlCrP+ESEqSolvInCr;
    s	=sqrt(0.5*ESBt/ESgF[ESNa]);
    for(i=0; i < nd; i++){
      lP->xd[i]	=AstA[i]*s;
      lJ->xd[i]	=lP->xd[i];
      ESxDtCr[i]=lP->xd[i];
    }
    lP->xd[0]	=0.;
    lJ->xd[0]	=0.;

    while(lP->kk[i]){
      lP->kk[i]	=0;
      lJ->kk[i]	=0;
      i++;
    }
    lJ->nd	=nAst;
    i	=nAst-1;
    ESjb[ESNa]	=AstJ[nAst+1];
    if(lJ->xd[i] < 1.){
      lJ->xd[nAst]	=1.;
      ESxDtCr[nAst]	=0.;
      ESSetSplA(lJ->xd[i]);
      lJ->yd[nAst]=lJ->yd[i];
      splRA(&s,NULL,ESg22c,ESg22c2);
      lJ->yd[nAst]	*=s;
      splRA(&s,NULL,ESLc,ESLc2);
      lJ->yd[nAst]	*=s;
      splRA(&s,NULL,ESaF,ESaF2a);
      lJ->yd[nAst]	*=s/(ESg22c[ESNa]*ESLc[ESNa]*ESaF[ESNa]);
#ifdef H
      lJ->yd[nAst]=0.2*AstIpl*ESaR0[0]/(ESg22c[ESNa]*ESLc[ESNa]*ESaF[ESNa]);
      lJ->yd[i]	*=s/ESxDtCr[i];
      AstJ[i]	=lJ->yd[i];
#endif
      lJ->kk[nAst]	=1;
      lJ->nd	=nAst+1;
      ESjb[ESNa]=0.;
    }
    ESPrData2Spl(lJ,lJ->kk,lJ->xd,lJ->yd,ESgaC+1,ESgaC+3);
    ESFirstSetSplDCr();
    ESPrSplData2CrPr(ESEqSolvInCr);
  }
#endif
#endif
#ifdef H
  /* Smoothing near the center */
  if(ESEqSolvInPr == 2){
    double gb,gB;
    gb	=ESjs[0]*lP->xd[1];
    gb	=4.*(lP->yd[0]-lP->yd[1])/(gb*gb);
    gB	=0.25*ESaR0[0]/lP->xd[1];
    EZout("sdd","gb=",gb,gB);
    if(gb > gB){
      lP->kk[1] =0;
      gB	*=(lP->yd[0]-lP->yd[1])/gb;
      gb	=lP->yd[0]-lP->yd[1];
      lP->yd[0]+=16./13.*(gB-gb);
      lP->yd[1] =lP->yd[0]-gB;
    }
  }
#endif
  if(lP->xin != xin){
    ConvertRadC(xin,lP->xin,lP->xd,lP->nd);
    lP->xin	=xin;
  }
  n	=lP->nd+2;
  for(i=0; i < n; i++) Y[i]	=0.;
#ifdef PFC
  if(FlY){
    if(ESEqSolvInPr == 0){
      PFSetPProfile(lP->xd,Y+2,n-2);
      Y[0]	=0.;
      Y[1]	=0.;
    }
    if(FlP/10 == 2){
      for(i=2; i < n; i++){
	lP->yd[i-2]	=Y[i];
      }
    }
    if(FlP/10 == 6){
      EFITaP2Jp(lP->xd,lP->yd,lP->nd);
    }
  }
#endif
  PrDataRadC	=lP->xin;
  kk	=lP->kk;
  ESnDtPr	=lP->nd;
  for(i=0; i < ESnDtPr; i++)  ESxDtPr[i]=lP->xd[i];

  X[0]	=0.;
  X[1]	=ESa0;
  for(i=2; i < n; i++){
    X[i]	=lP->xd[i-2];
  }
#ifdef STUPID
  if(PrDataRadC != ESEqSolvRadC){
    double a,da;
    for(i=0; i < ESnDtPr; i++){
      do{
	ESSetSplA(X[i]);
	splRA(&a,&da,RdC[PrDataRadC].x0,RdC[PrDataRadC].x2);
	da	=(lP->xd[i]-a)/da;
	X[i]	+=da;
      }while(fabs(da) > 1e-8);
    }
    printf("Ps[0]=%12.5e %12.5e %12.5e\n",ESPs[0]*0.091135*EZcgm0/ESsp[0]
	   ,ESsp[0]*EZcrgm0*1e+6,ESPs[0]);
  }
#endif
  xd	=X+2;
  x	=lP->x0;
  X1[0]	=0.;
  X1[1]	=ESa0;

  ESPrData2Spl(lP,lP->kk,lP->xd,lP->yd,ESgaP+1,ESgaP+3);
  ESnDtPr--;
  ESFirstSetSplDPr();
  ESPrSplData2PrPr(ESEqSolvInPr);
#ifdef XWIN
  EC1DEqSolv();
  ESPlProf2Data();
#endif

#ifdef MSE
  if(FlP >= 10){
    for(i=0; i < ESXNa1; i++){
      ESXsxD[i]	=ESXsx[i];
      ESXsrD[i]	=ESXsr[i];
      ESXsjD[i]	=ESXsj[i];
      ESXspD[i]	=ESXsp[i];
      ESXsqD[i]	=ESXsq[i];
      ESXsiD[i]	=ESXsi[i];
    }
    FlP	=ESEqSolvInPr;
  }
  Scale2D(4,2,ESXsx,ESXsi,ESXNa1,6,2);
  Plot2D(4,ESXsx,ESXsi,ESXNa1,6,0,0,0);
  Plot2D(4,ESXsxD,ESXsiD,ESXNa1,6,1,0,0);
#endif
  da	=(ESsa[ESNa]-ESsa[0])/(n1-1);
  for(i=0; i < n1; i++){
    a[i]	=ESa0+da*i;
    ESSetSplA(a[i]);
    splRA(f0+i,NULL,ESjp,ESjp2a);
    splRA(f1+i,NULL,ESsp,ESsp2a);
    splRA(f2+i,NULL,ESPs,ESPs2a);
    splRA(f3+i,NULL,ESgB,ESgB2a);
    splRA(f4+i,NULL,ESgb,ESgb2a);
  }
  ld	=lld[ii[ESEqSolvInPr]];
  for(i=0; i < n1; i++){
    a[i]	=ESa0+da*i;
    ESSetSplDPr(a[i]);
    splRDPr(ld+i,NULL,ESEqSolvInPr);
  }
#ifdef PFC
  CbUserPlot	=EFITjpPlot;
#endif
#endif/*Tbl_PlPr*/

#ifndef AscOPr
#undef AscOPr
#endif/*AscOPr*/
  return(0);
}

#define TTTTT
#ifdef TTT
int nnd;
double eSq0=1.45,*pyd,*pyd0,*pyd2;
#endif

#ifdef PFC
static int MSEnSpl=0,MSEna=0;
static double MSEr[64],MSEBz[64],MSEd2Bz[64];
static double MSEa[64],MSEdgY[64],MSEdgY2a[64];
static double MSEdx=0.;

void ESInpMSEPlot()
{
  int i;
  double r[33],Bz[33],dr;
  
  ZColor(0);
  r[0]	=0.;
  r[1]	=0.;
  for(i=0; i < nMSE; i++){
    if((MSE[i].kA&0x03) == 3){
      ZPlotMark(MSE[i].R,MSE[i].Bz+MSE[i].dBz,0x0110);
      ZPlotMark(MSE[i].R,MSE[i].Bz-MSE[i].dBz,0x0110);
      ZPlotMark(MSE[i].R,MSE[i].Bz,0x0004);
    }
    if((MSE[i].kA&0x08)){
      if(r[0] > MSE[i].Bz-MSE[i].dBz){
	r[0]	=MSE[i].Bz-MSE[i].dBz;
      }

      if(r[1] < MSE[i].Bz+MSE[i].dBz){
	r[1]	=MSE[i].Bz+MSE[i].dBz;
      }
    }
  }
  ZColor(4);
  for(i=0; i < nMSE; i++){
    if((MSE[i].kA&0x03) == 2){
      ZPlotMark(MSE[i].R,MSE[i].Bz+MSE[i].dBz,0x0110);
      ZPlotMark(MSE[i].R,MSE[i].Bz-MSE[i].dBz,0x0110);
      ZPlotMark(MSE[i].R,MSE[i].Bz,0x0004);
    }
  }
  ZColor(10);
  for(i=0; i < nMSE; i++){
    if((MSE[i].kA&0x0b) == 0x08){
      ZPlotMark(MSE[i].R,MSE[i].Bz+MSE[i].dBz,0x0110);
      ZPlotMark(MSE[i].R,MSE[i].Bz-MSE[i].dBz,0x0110);
      ZPlotMark(MSE[i].R,MSE[i].Bz,0x0004);
    }
  }
  i	=ESNp1*ESNa;
  ZPlotLine(ESsr[i],r[0],ESsr[i],r[1]);
  ZColor(14);
  for(i=0; i < nMSE; i++){
    if((MSE[i].kA&0x03) == 1){
      ZPlotMark(MSE[i].R,MSE[i].Bz+MSE[i].dBz,0x0110);
      ZPlotMark(MSE[i].R,MSE[i].Bz-MSE[i].dBz,0x0110);
      ZPlotMark(MSE[i].R,MSE[i].Bz,0x0004);
    }
  }
  ZPlotLine(ESaR0[0],r[0],ESaR0[0],r[1]);
  if(MSEnSpl == 0){
    return;
  }
  dr	=(MSEr[MSEnSpl]-MSEr[0])/32.;
  for(i=0; i < 33; i++){
    r[i]	=MSEr[0]+dr*i;
    SetIspl(MSEr,r[i],MSEnSpl);
    splR(Bz+i,NULL,MSEBz,MSEd2Bz);
  }
  ZPlotPolyLine(r,Bz,33);

  return;
}

int MSEPitchAngleApproxPlot(int W)
{
  int i,j,n;
  double r[33],Bz[33],dr;
  double gY,d,Re,x[nMSE],y[nMSE],y1[nMSE],y2[nMSE];
  double dRa,dRp,dZa,dZp;
  double Bt,rqg;
  
  Re	=ESsr[ESNp1*ESNa];

  Bz[0]	=0.;
  Bz[1]	=0.;
  j	=0;
  for(i=0; i < nMSE; i++){
    if((MSE[i].kA&0x08)){
      d		=MSE[i].Bz;
      y[j]	=MSE[i].Bz;
      x[j]	=MSE[i].R;
      if((MSE[i].kA&0x02)){
	if(Bz[0] > d){
	  Bz[0]	=d;
	}
	if(Bz[1] < d){
	  Bz[1]	=d;
	}
      }
      ESGetMetrics(&dRa,&dRp,&dZa,&dZp,MSE[i].a,MSE[i].gt);
      splRA(&gY,&d,ESgY,ESgY2a);
      d		*=dZp/(MSE[i].R*(dRa*dZp-dRp*dZa));
      if(Bz[0] > d){
	Bz[0]	=d;
      }
      if(Bz[1] < d){
	Bz[1]	=d;
      }
      y1[j]	=d;
      SetIspl(MSEa,MSE[i].a,MSEna);
      splR(&d,NULL,MSEdgY,MSEdgY2a);
      y2[j]	=d*dZp/(MSE[i].R*(dRa*dZp-dRp*dZa));
      j++;
    }
  }
  i	=nMSE-1;
  r[0]	=MSE[0].R < MSE[i].R ? MSE[0].R : MSE[i].R;
  r[1]	=MSE[0].R < MSE[i].R ? MSE[i].R : MSE[0].R;
  if(r[1] < Re){
    r[1]	=Re;
  }

  SetPlotName("R","Bz/Bt'","MSE Pitch angle");
  Scale2d(W,r,Bz,2,6,2);
  Bz[2]	=0.;
  Bz[3]	=0.;
  r[2]	=ESaR0[0];
  r[3]	=ESaR0[0];
  r[4]	=Re;
  r[5]	=Re;
  Plot2d(W,r,Bz+2,2,6,0,9,0);
  Plot2d(W,r+2,Bz,2,6,0,9,0);
  r[2]	+=MSEdx;
  r[3]	+=MSEdx;
  Plot2d(W,r+2,Bz,2,6,0,12,0);
  Plot2d(W,r+4,Bz,2,6,0,9,0);
  Plot2d(W,x,y,j,6,0x08,0,0);
  Plot2d(W,x,y1,j,6,0x08,4,0);
  Plot2d(W,x,y2,j,6,0,14,0);

  if(MSEnSpl == 0){
    return(0);
  }
  dr	=(MSEr[MSEnSpl]-MSEr[0])/32.;
  j	=0;
  for(i=0; i < 33; i++){
    r[j]	=MSEr[0]+dr*i;
    d		=r[j]+MSEdx;
    if(MSEr[0] <= d && d <=  MSEr[MSEnSpl]){
      SetIspl(MSEr,d,MSEnSpl);
      splR(Bz+j,NULL,MSEBz,MSEd2Bz);
      j++;
    }
  }
  Plot2d(W,r,Bz,j,6,0,4,0);
  j	=0;
  for(i=0; i < 33; i++){
    r[j]	=MSEr[0]+dr*i;
    d		=r[j];
    SetIspl(MSEr,d,MSEnSpl);
    splR(Bz+j,NULL,MSEBz,MSEd2Bz);
    j++;
  }
  Plot2d(W,r,Bz,j,6,0,12,0);
  return(0);
}

int MSEdgYApproxPlot(int W)
{
  int i,j;
  double a[33],dgY[33],da;
  double d,x[nMSE],y[nMSE],y1[nMSE];
  
  dgY[0]	=0.;
  dgY[1]	=0.3;
  j	=0;
  for(i=0; i < nMSE; i++){
    if((MSE[i].kA&0x08) && fabs(MSE[i].gt) <= EZcr2*EZcgp){
      d		=MSE[i].dgY;
      if(dgY[0] > d){
	dgY[0]	=d;
      }		
      x[j]	=MSE[i].a;
      y[j]	=MSE[i].dgY;
      if(MSE[i].a >= ESa0){
	ESSetSplA(MSE[i].a);
	splRA(&d,y1+j,ESgY,ESgY2a);
      }
      else{
	y1[j]	=0.;
      }
      if(dgY[0] > y1[j]){
	dgY[0]	=y1[j];
      }		
      j++;
    }
  }
  a[0]	=0.;
  a[1]	=ESsa[ESNa];
  SetPlotName("R","gY'","MSE gY'");
  Scale2d(W,a,dgY,2,6,2);
  dgY[2]	=0.;
  dgY[3]	=0.;
  Plot2d(W,a,dgY+2,2,6,0,9,0);
  Plot2d(W,x,y,j,6,0x08,0,0);
  Plot2d(W,x,y1,j,6,0x08,4,0);
  j	=0;
  for(i=0; i < nMSE; i++){
    if((MSE[i].kA&0x08) && fabs(MSE[i].gt) > EZcr2*EZcgp){
      x[j]	=MSE[i].a;
      y[j]	=MSE[i].dgY;
      ESSetSplA(MSE[i].a);
      splRA(&d,y1+j,ESgY,ESgY2a);
      j++;
    }
  }
  Plot2d(W,x,y,j,6,0x08,14,0);

  if(MSEna == 0){
    return(0);
  }
  da	=MSEa[MSEna]/32.;
  for(i=0; i < 33; i++){
    a[i]	=da*i;
    SetIspl(MSEa,a[i],MSEna);
    splR(dgY+i,NULL,MSEdgY,MSEdgY2a);
  }
  Plot2d(W,a,dgY,33,6,0,4,0);
  return(0);
}

int MSEBz2Spline(double EZga1,double EZga3,double dX,double dY)
{
  int i,j,k,n;
  double s,d,EZga0=0.,EZga2=0.;

  if(nMSE == 0){
    return(0);
  }
  n	=0;
  for(i=0; i < nMSE && i <64; i++){
    if((MSE[i].kA&0x02)){
      MSEr[n]	=MSE[i].R+dX;
      MSEBz[n]=MSE[i].Bz+dY;
      n++;
    }
  }
  MSEnSpl	=n-1;
  /* ordering */
  for(i=0; i < n; i++){
    k	=i;
    for(j=i+1; j < n; j++){
      if(MSEr[j] < MSEr[k]){
	k	=j;
      }
    }
    if(k != i){
      s	=MSEr[k];
      MSEr[k]	=MSEr[i];
      MSEr[i]	=s;
      s	=MSEBz[k];
      MSEBz[k]	=MSEBz[i];
      MSEBz[i]	=s;
    }
  }

  EZf2spl(MSEBz,MSEd2Bz,NULL,NULL,MSEBz,NULL,MSEr,&MSEnSpl,&EZga0,&EZga1,&EZga2,&EZga3);

  i=1;
  while(i < n &&  MSEBz[i]*MSEBz[i-1] > 0.){
    i++;
  }
  s	=MSEBz[i]-MSEBz[i-1];
  MSEdx	= s != 0. ? MSEr[i-1]-(MSEr[i]-MSEr[i-1])*MSEBz[i-1]/s : 
    EZcr2*(MSEr[i]-MSEr[i-1]);
  MSEdx	-=ESaR0[0];
  do{
    SetIspl(MSEr,ESaR0[0]+MSEdx,MSEnSpl);
    splR(&s,&d,MSEBz,MSEd2Bz);
    s	/=-d;
    MSEdx	+=s;
  }while(fabs(s) > 1e-8);
  return(n);
}

int MSEdgY2Spline(double EZga1,double EZga3)
{
  int i,j,k,n;
  double dRa,dRp,dZa,dZp;
  double s,d,EZga0=0.,EZga2=0.;

  if(nMSE == 0){
    return(0);
  }
  n	=0;
  MSEa[n]	=0.;
  MSEdgY[n]	=0.;
  n++;
  for(i=0; i < nMSE && i <64; i++){
    if((MSE[i].kA&0x08)){
      MSEa[n]	=MSE[i].a;
      s		=MSE[i].R+MSEdx;
      if(MSEr[0] <= s && s <=  MSEr[MSEnSpl]){
	SetIspl(MSEr,s,MSEnSpl);
	splR(&s,NULL,MSEBz,MSEd2Bz);
	ESGetMetrics(&dRa,&dRp,&dZa,&dZp,MSE[i].a,MSE[i].gt);
	MSEdgY[n]=s*MSE[i].R*(dRa*dZp-dRp*dZa)/dZp;
      }
      else{
	if(MSE[i].a >= ESa0){
	  ESSetSplA(MSE[i].a);
	  splRA(&s,&d,ESgY,ESgY2a);
	  MSEdgY[n]=d;
	}
	else{
	  MSEdgY[n]=0.;
	}
      }
      n++;
    }
  }
  dRa	=MSEa[n-1];
  i	=ESNa*dRa;
  while(i < ESNa1){
    if(ESsa[i] > dRa){
      MSEa[n]	=ESsa[i];
      ESSetSplA(ESsa[i]);
      splRA(&s,&d,ESgY,ESgY2a);
      MSEdgY[n]=d;
      n++;
    }
    i++;
  }

  /* ordering */
  for(i=0; i < n; i++){
    k	=i;
    for(j=i+1; j < n; j++){
      if(MSEa[j] < MSEa[k]){
	k	=j;
      }
    }
    if(k != i){
      s	=MSEa[k];
      MSEa[k]	=MSEa[i];
      MSEa[i]	=s;
      s	=MSEdgY[k];
      MSEdgY[k]	=MSEdgY[i];
      MSEdgY[i]	=s;
    }
  }
  /* unfoolding */
  for(i=0,k=1; k < n; i++,k++){
    while(fabs(MSEa[i]-MSEa[k]) < 0.1*ESsa[1]){
      MSEa[i]	=EZcr2*(MSEa[i]+MSEa[k]);
      MSEdgY[i]	=EZcr2*(MSEdgY[i]+MSEdgY[k]);
      n--;
      for(j=k; j < n; j++){
	MSEa[j]		=MSEa[j+1];
	MSEdgY[j]	=MSEdgY[j+1];
      }
    }
  }
  MSEna	=n-1;
  EZf2spl(MSEdgY,MSEdgY2a,NULL,NULL,MSEdgY,NULL,MSEa,&MSEna,&EZga0,&EZga1,&EZga2,&EZga3);
  return(n);
}

int MSEdgY2js()
{
  int i,j,k,n;
  double a,g,dg,jp,L,V,dgY,d2gY;

  n	=PlCrP[0].nd;
  a	=ESa0;
  SetIspl(MSEa,a,MSEna);
  splR(&dgY,&d2gY,MSEdgY,MSEdgY2a);
  g	=ESg22c[0];
  L	=ESLc[0];
  PlCrP[0].yd[0]	=-2.*g*d2gY/L;
  for(i=1; i < n; i++){
    a	=PlCrP[0].xd[i];
    SetIspl(MSEa,a,MSEna);
    splR(&dgY,&d2gY,MSEdgY,MSEdgY2a);
    ESSetSplA(a);
    splRA(&g,&dg,ESg22c,ESg22c2);
    splRA(&L,NULL,ESLc,ESLc2);
    splRA(&V,NULL,ESVc,ESVc2);
    splRA(&jp,NULL,ESjp,ESjp2a);
    PlCrP[0].yd[i]	=-(g*dgY/a+g*d2gY+dg*dgY+V*jp)/L;
  }
  return(0);
}

int ESInpMSE()
{
  static int Fl=1;
  static double EZga0=0.,EZga1=1e-9,EZga2=0.,EZga3=5e-8,Yshft=0.,Xshft=0.;
  int i,N,n,kk[nMSE];
  double x[nMSE],y[nMSE];

  if(nMSE == 0){
    return(0);
  }
  MSEtangg();
  MSEPitchAngle();
  n	=0;
  for(i=0; i < nMSE && i < 64; i++){
    if((MSE[i].kA&0x08)){
      x[n]	=MSE[i].R;
      y[n]	=MSE[i].Bz;
      kk[n]	=(MSE[i].kA&0x02) ? 1 : 0;
      n++;
    }
  }
#ifndef Tbl_MSE
#ifndef AscIMSE
#undef AscIMSE
#endif/*AscIMSE*/
  if(Fl){
    n	=0;
    for(i=0; i < nMSE; i++){
      MSE[i].kA	&=0x09;
      if((MSE[i].kA&0x01)){
	MSE[i].kA	|=0x02;
      }
      if((MSE[i].kA&0x08)){
	kk[n]	=(MSE[i].kA&0x02) ? 1 : 0;
	n++;
      }
    }
    Fl	=0;
  }
  n	=0;
  for(i=0; i < nMSE && i <64; i++){
    if((MSE[i].kA&0x08)){
      MSE[i].kA	=kk[n] ? MSE[i].kA|0x02 : MSE[i].kA&0x09;
      n++;
    }
  }
  N	=MSEBz2Spline(EZga1,EZga3,Xshft,Yshft);

#ifdef PFC11
  MSEPitchAnglePLot(4);
  MSEDataPlot(6,0);
#endif

  CbUserPlot	=ESInpMSEPlot;
#endif/*Tbl_MSE*/

#ifndef AscOMSE
#undef AscOMSE
#endif/*AscOMSE*/
  return(0);
}

int ESInpPlCr4PF()
{
  PlProf *lP;
  int i,k,n;
  double r0,Ipl,da;

  lP	=PlCrP+ESEqSolvInCr;

  for(i=0; i < NCrP; i++){
    PlCrP[i].io	=i == ESEqSolvInCr ? 1 : 0;
    ioCrP[i]=PlCrP[i].io;
  }
  lP->nd	=ESData2Order(lP->xd,lP->yd,lP->kk,NPlPD);
  n	=lP->nd;
  CrDataRadC	=lP->xin;
  for(i=0; i < n; i++){
    ESxDtCr[i]	=lP->xd[i];
  }
  ESnDtCr	=n-1;
  ESFirstSetSplDCr();
  if(ESEqSolvInCr == 6){
    double dg0,dg1;
    PlCrP[7].nd	=lP->nd;
    for(i=0; i < lP->nd; i++){
      PlCrP[7].kk[i]	=1;
      PlCrP[7].xd[i]	=lP->xd[i];
    }
    for(i=0; i < lP->nd; i++){
      PlCrP[7].yd[i]	=lP->yd[i] > 1e-4 ? 1./lP->yd[i] 
	: 1./(lP->yd[i]+1e-4);
    }
    ESPrData2Spl(PlCrP+7,lP->kk,lP->xd,PlCrP[7].yd,ESgaC+1,ESgaC+3);
    ESPrSplData2CrPr(7);
    for(i=0; i < ESNa1; i++){
      ESsq[i]	=ESgm[i] > 1e-4 ? 1./ESgm[i] 
	: 1./(ESgm1a[i]*0.01*ESsa[ESNa]);
    }
    dg0	=ESa0 != 0. ? 0. : -ESgm1a[0]*ESsq[0]*ESsq[0];
    dg1	=-ESgm1a[ESNa]*ESsq[ESNa]*ESsq[ESNa];
    splAA(ESsq,ESsq1a,ESsq2a,&dg0,&dg1);
  }
  else{
    ESPrData2Spl(lP,lP->kk,lP->xd,lP->yd,ESgaC+1,ESgaC+3);
    ESPrSplData2CrPr(ESEqSolvInCr);
#ifdef H
    if(ESkNormJ && (ESEqSolvInCr == 0 || ESEqSolvInCr == 1)){
      double u[2],v[2],ss,s,aJ,dJ;
      int j;

      EC1DEqSolv();
      ss	=1.;
      aJ	=5.*ESaJ[ESNa];
      u[0]	=1.;
      u[1]	=ESIpl/aJ;
      j		=0;
      do{
	EZzero(&j,u,v,&s,dJ);
	s	/=ss;
	for(i=0; i < n; i++){
	  lP->yd[i]	*=s;
	  lP->d0[i]	*=s;
	  lP->d2[i]	*=s;
	}
	ss	*=s;
	ESPrSplData2CrPr(ESEqSolvInCr);
	EC1DEqSolv();
	dJ	=ESIpl-5.*ESaJ[ESNa];
      }while(j > 0 && j < 10 && fabs(dJ/ESIpl) > 1e-6);
    }
#endif
  }
  if(ESEqSolvInCr == 7){
    double dg0,dg1;
    for(i=0; i < lP->nd; i++){
      PlCrP[6].kk[i]	=lP->kk[i];
      PlCrP[6].xd[i]	=lP->xd[i];
      PlCrP[6].yd[i]	=lP->yd[i] > 1e-4 ? 1./lP->yd[i] 
	: 1./(ESgm1a[i]*0.01*ESsa[ESNa]);
    }
    for(i=0; i < ESNa1; i++){
      ESsq[i]	=ESgm[i] > 1e-4 ? 1./ESgm[i] 
	: 1./(ESgm1a[i]*0.01*ESsa[ESNa]);
    }
    dg0	=ESa0 != 0. ? 0. : -ESgm1a[0]*ESsq[0]*ESsq[0];
    dg1	=-ESgm1a[ESNa]*ESsq[ESNa]*ESsq[ESNa];
    splAA(ESsq,ESsq1a,ESsq2a,&dg0,&dg1);
  }
  if(ESEqSolvInCr == 1){
    r0	=ESaR0[0];
    for(i=0; i < lP->nd; i++){
      PlCrP[2].kk[i]	=lP->kk[i];
      PlCrP[2].xd[i]	=lP->xd[i];
      PlCrP[2].yd[i]	=r0*lP->yd[i];
    }
    PlCrP[2].nd	=lP->nd;
    for(i=0; i < ESNa1; i++){
      ESjR[i]	=ESjb[i]*r0;
      ESjR1a[i]	=ESjb1a[i]*r0;
      ESjR2a[i]	=ESjb2a[i]*r0;
    }
  }
  if(ESEqSolvInCr == 2){
    r0	=1./ESaR0[0];
    for(i=0; i < lP->nd; i++){
      PlCrP[1].kk[i]	=lP->kk[i];
      PlCrP[1].xd[i]	=lP->xd[i];
      PlCrP[1].yd[i]	=r0*lP->yd[i];
    }
    PlCrP[1].nd	=lP->nd;
    for(i=0; i < ESNa1; i++){
      ESjb[i]	=ESjR[i]*r0;
      ESjb1a[i]	=ESjR1a[i]*r0;
      ESjb2a[i]	=ESjR2a[i]*r0;
    }
  }
  EC1DEqSolv();
  ESPlProf2Data();
  Ipl	=ESkNormJ && (ESEqSolvInCr == 0 || ESEqSolvInCr == 1) ? 
    ESIpl : 5.*ESaJ[ESNa];
  return(0);
}
#endif

int ESInpPlCr()
{
  extern char ESMessage[];
  extern double ESjbs[];
  PlProf *lP;
  static int FlC=0,FlY=0,xin;
  int i,k,n,io,InC,Nax,nd;
  int N=257;
  int kk[N];
  double XD[ESNRDCR][N],*Xd[ESNRDCR],X[N+2],Y[N+2],X1[2]={0.,0.};
  extern int FlagPS;
  double *xd,*x,r0,Ipl,da;
  static double gg1=1e-9,gg3=5e-8;
  static double gm1,am1=0.5,gm2=-0.1;

  static int ii[18]={0,1,2,10,3,0,5,6,7,8,9,0,0,0,0,0,11,4};
  int n1=3*ESNa1+1,n2;
  double *ld,*lld[12];
  double *ad[ESNRDCR],aD[ESNRDCR][n1],a[n1],f0[n1],f1[n1],f2[n1],f3[n1]
    ,f4[n1],f5[n1],f6[n1],f7[n1],f8[n1],f9[n1],f10[n1],f11[n1];

  lld[0]	=f0;
  lld[1]	=f1;
  lld[2]	=f2;
  lld[3]	=f3;
  lld[4]	=f4;
  lld[5]	=f5;
  lld[6]	=f6;
  lld[7]	=f7;
  lld[8]	=f8;
  lld[9]	=f9;
  lld[10]	=f10;
  lld[11]	=f11;

  n2	=n1+2;
  X1[1]	=ESa0;
  for(i=0; i < ESNRDCR; i++){
    Xd[i]	=XD[i];
    ad[i]	=aD[i];
  }
  for(i=0; i < N; i++){
    kk[i]	=0;
  }
  Nax	=ESNa1;
  InC	=ESEqSolvInCr;
  xin	=CrDataRadC;
  if(FlC%10 != 9){
    FlC	=ESEqSolvInCr+(FlC/10)*10;
  }
#ifndef Tbl_PlCr
  *ESMessage	='\0';
  if(ESa0 != 0.){
    ESHoleEdgeBoundary();
    ES3DMetrTnsrCSHole();
  }
  if(FlC/10 == 3){
#ifndef AscICr
#undef AscICr
#endif/*AscICr*/
    FlC	%=10;
  }
  if((FlC%10) == 9){
#ifdef STUPID
    double x;
    for(i=0; i < 21; i++){
      ESSetSplA(PlCrP[7].xd[i]);
      splRA(&x,NULL,ESqgY,ESqgY2a);
      PlCrP[7].yd[i]	=1./ESWonchulQ(x*x);
      PlCrP[7].kk[i]	=1;
    }
    ESEqSolvInCr	=7;
#endif
  }
  else{
    ESEqSolvInCr	=FlC%10;
  }
  if(ESEqSolvInCr == 5 || ESEqSolvInCr > 7){
    printf("InC == 5 or > 7 are not allowed%c\n",'\a');
    ESEqSolvInCr	=InC;
  }
  if(ESEqSolvInPr == 2 && (ESEqSolvInCr == 0 || ESEqSolvInCr == 3)){
    printf("InC == %d is not allowed for given p(a) == 2%c\n"
	   ,ESEqSolvInCr,'\a');
    ESEqSolvInCr	=InC;
  }
  if(ESEqSolvInCr > ESNCRPR-1){
    ESEqSolvInCr	=ESNCRPR-1;
  }
  io	=ESEqSolvInCr;
  lP	=PlCrP+ESEqSolvInCr;
  for(i=0; i < NCrP; i++){
    PlCrP[i].io	=i == ESEqSolvInCr ? 1 : 0;
    ioCrP[i]	=PlCrP[i].io;
  }
#ifdef STUPID
  Njpol	=0;
  for(i=0; i < N; i++){
    if(lP->kk[i] == 2 && Njpol < 20){
      Xjk[Njpol]	=lP->xd[i];
      Pjk[Njpol]	=lP->yd[i];
      Njpol++;
    }
    qgy[i]	=lP->xd[i];
  }
  if(Njpol){
    if(ESData2Poly(Xjk,Pjk, Njpol)){
      for(i=0; i < NPlPD; i++){
	if(lP->kk[i] == 2){
	  lP->kk[i]	=1;
	}
      }
      Njpol	=0;
    }
  }
#endif
  lP->nd=ESData2Order(lP->xd,lP->yd,lP->kk,NPlPD);
  nd	=lP->nd;
  if(lP->xin != xin){
    ConvertRadC(xin,lP->xin,lP->xd,lP->nd);
    lP->xin	=xin;
  }
  if(ESkNormJ == 2){
    ESIpl	=Ipl;
    ESkNormJ	=0;
  }
  if(ESkNormJ == 3){
    ESIpl	=Ipl;
    ESkNormJ	=1;
  }
  n	=lP->nd+2;
  for(i=0; i < n; i++) Y[i]=0.;
#ifdef PFC
  if(FlY){
    if(ESEqSolvInCr == 0){
      PFSetJProfile(lP->xd,Y+2,n-2);
      Y[0]	=0.;
      Y[1]	=0.;
    }
    if(FlC/10 == 1){
      for(i=2; i < n; i++) lP->yd[i-2]=Y[i];
      ESkNormJ=0;
    }
    if(FlC/10 == 6){
      EFITaT2Js(lP->xd,lP->yd,lP->nd);
    }
  }
#endif
  CrDataRadC	=lP->xin;
  ESnDtCr	=lP->nd-1;
  for(i=0; i < lP->nd; i++){
    ESxDtCr[i]	=lP->xd[i];
  }
  X[0]	=0.;
  X[1]	=ESa0;
  for(i=2; i < n; i++){
    X[i]	=lP->xd[i-2];
  }
  xd	=lP->xd;
  xd	=X+2;
  ESFirstSetSplDCr();
  if(ESEqSolvInCr == 6){
    double dg0,dg1;
    PlCrP[7].nd	=lP->nd;
    for(i=0; i < lP->nd; i++){
      PlCrP[7].kk[i]	=1;
      PlCrP[7].xd[i]	=lP->xd[i];
    }
    for(i=0; i < lP->nd; i++){
      PlCrP[7].yd[i]	=lP->yd[i] > 1e-4 ? 1./lP->yd[i]
	: 1./(lP->yd[i]+1e-4);
    }
    ESPrData2Spl(PlCrP+7,lP->kk,lP->xd,PlCrP[7].yd,ESgaC+1,ESgaC+3);
    ESPrSplData2CrPr(7);
    for(i=0; i < ESNa1; i++){
      ESsq[i]	=ESgm[i] > 1e-4 ? 1./ESgm[i] 
	: 1./(ESgm1a[i]*0.01*ESsa[ESNa]);
    }
    dg0	=-ESgm1a[0]*ESsq[0]*ESsq[0];
    dg1	=-ESgm1a[ESNa]*ESsq[ESNa]*ESsq[ESNa];
    splAA(ESsq,ESsq1a,ESsq2a,&dg0,&dg1);
  }
  else{
    ESPrData2Spl(lP,lP->kk,lP->xd,lP->yd,ESgaC+1,ESgaC+3);
    ESPrSplData2CrPr(ESEqSolvInCr);
    EC1DEqSolv();
#ifdef H
    if(ESkNormJ && ESEqSolvInCr < 3){
      double u[2],v[2],ss,s,aJ,dJ;
      int j;

      EC1DEqSolv();
      ss	=1.;
      aJ	=5.*ESaJ[ESNa];
      u[0]	=1.;
      u[1]	=ESIpl/aJ;
      j		=0;
      do{
	EZzero(&j,u,v,&s,dJ);
	s	/=ss;
	for(i=0; i < ESnDtCr; i++){
	  lP->yd[i]	*=s;
	  lP->d0[i]	*=s;
	  lP->d2[i]	*=s;
	}
	ss	*=s;
	ESPrSplData2CrPr(ESEqSolvInCr);
	EC1DEqSolv();
	dJ	=ESIpl-5.*ESaJ[ESNa];
      }while(j > 0 && j < 10 && fabs(dJ/ESIpl) > 1e-6);
    }
#endif
  }

  if(ESEqSolvInCr == 7){
    double dg0,dg1;
    for(i=0; i < lP->nd; i++){
      PlCrP[6].kk[i]	=lP->kk[i];
      PlCrP[6].xd[i]	=lP->xd[i];
      PlCrP[6].yd[i]	=lP->yd[i] > 1e-4 ? 1./lP->yd[i]
	: 1./(ESgm1a[i]*0.01*ESsa[ESNa]);
    }
    for(i=0; i < ESNa1; i++){
      ESsq[i]	=ESgm[i] > 1e-4 ? 1./ESgm[i] 
	: 1./(ESgm1a[i]*0.01*ESsa[ESNa]);
    }
    dg0	=-ESgm1a[0]*ESsq[0]*ESsq[0];
    dg1	=-ESgm1a[ESNa]*ESsq[ESNa]*ESsq[ESNa];
    splAA(ESsq,ESsq1a,ESsq2a,&dg0,&dg1);
  }

  InC	=ESEqSolvInCr;
  EC1DEqSolv();
  ESPlProf2Data();
#ifdef XWIN
  x	=lP->x0;
  Ipl	=ESkNormJ && ESEqSolvInCr < 3 ? ESIpl : 5.*ESaJ[ESNa];
#endif

#ifdef PFC
  if(nMSE){
    if(ESa0 == 0.){
      MSEtangg();
#ifdef H
      MSEPitchAngle();
      MSEdgY2Spline(gg1,gg3);
      MSEPitchAngleApproxPlot(5);
      MSEdgYApproxPlot(6);
#endif
    }
    else{
      MSEtanggHole();
    }
    MSEDataPlot(4,0);
    if(FlC/10 == 2){
      MSEdgY2js();
    }
  }
  CbUserPlot	=EFITjsPlot;
#endif
 
  GenerateRadC(Xd,lP->xd,lP->xin,lP->nd);
  da	=(ESsa[ESNa]-ESsa[0])/(n1-1);
  for(i=0; i < n1; i++){
    a[i]	=ESa0+da*i;
    ESSetSplA(a[i]);
    splRA(f0+i,NULL,ESjs,ESjs2a);
    splRA(f1+i,NULL,ESjb,ESjb2a);
    splRA(f2+i,NULL,ESjR,ESjR2a);
    splRA(f3+i,NULL,ESjB,ESjB2a);
    splRA(f4+i,NULL,ESaJ,ESaJ2a);
    splRA(f5+i,NULL,ESsq,ESsq2a);
    splRA(f6+i,NULL,ESgm,ESgm2a);
#ifdef H
    f4[i]	=f5[i] > 1e-4 ? 1./f5[i] : 0.;
#endif
    splRA(f7+i,NULL,ESgY,ESgY2a);
    splRA(f8+i,NULL,ESgF,ESgF2a);
    splRA(f9+i,NULL,ESaF,ESaF2a);
    splRA(f10+i,NULL,ESaT,ESaT2a);
    splRA(f11+i,NULL,ESdgY,ESdgY2a);
  }
  ld	=lld[ii[ESEqSolvInCr]];
  for(i=0; i < n1; i++){
    a[i]	=ESa0+da*i;
    ESSetSplDCr(a[i]);
    splRDCr(ld+i,NULL,ESEqSolvInCr);
  }
  GenerateRadC(ad,a,ESEqSolvRadC,n1);

#ifdef XWIN
#ifdef H
  {
    int j;
    double z,js,jp;
    double x[2],y[2],rx[2*n1],jx[2*n1]; 
    y[0]	=0.;
    y[1]	=0.;
    for(i=0; i < n1; i++){
      ESSetSplA(a[i]);
      splRA(&jp,NULL,ESjp,ESjp2a);
      js	=f0[i];
      j	=n1-1+i;
      ES2DMapFlx2Lab(rx+j,&z,&z,&z,&z,&z,a[i],0.);
      z	=ESaR0[0]/rx[j];
      jx[j]	=js*z+jp*(1./z-z);
      if(y[0] > jx[j]){
	y[0]	=jx[j];
      }
      if(y[1] < jx[j]){
	y[1]	=jx[j];
      }
      j	=n1-1-i;
      ES2DMapFlx2Lab(rx+j,&z,&z,&z,&z,&z,a[i],EZcgp);
      z	=ESaR0[0]/rx[j];
      jx[j]	=js*z+jp*(1./z-z);
      if(y[0] > jx[j]){
	y[0]	=jx[j];
      }
      if(y[1] < jx[j]){
	y[1]	=jx[j];
      }
    }
    x[0]	=rx[0];
    x[1]	=rx[2*(n1-1)];
    Scale2D(5,2,x,y,2,6,2);
    Plot2d(5,rx,jx,2*n1-1,6,0,14,0);
  }
#endif
#endif

  gm1	=ESgm[0];
  k	=0;
  for(i=1; i < ESNa1 && ESsa[i] < am1; i++){
    if(gm1 > ESgm[i]){
      gm1=ESgm[i];
      k	=i;
    }
  }
#ifdef XWIN
  if(gm1 < gm2){
    CbUserWarning	=ESMessage;
    sprintf(ESMessage,"Too small gm at a=%5.3f:\n gm[%d]=%10.3e <%10.3e"
	    ,ESsa[k],k,ESgm[k],gm2);
  }		
#endif

#endif/*Tbl_PlCr*/
#ifndef AscOCr
#undef AscOCr
#endif/*AscOCr*/
  return(0);
}


int ESInpPr2Spl()
{
  int i,k;
  PlProf *lP,*lJ;

  lP	=PlPrP+ESEqSolvInPr;
  ESPrData2Spl(lP,lP->kk,lP->xd,lP->yd,ESgaP+1,ESgaP+3);
  ESFirstSetSplDPr();
  ESPrSplData2PrPr(ESEqSolvInPr);

  lJ	=PlCrP+ESEqSolvInCr;
  k	=lJ->nd;
  ESPrData2Spl(lJ,lJ->kk,lJ->xd,lJ->yd,ESgaC+1,ESgaC+3);
  ESFirstSetSplDCr();
  ESPrSplData2CrPr(ESEqSolvInCr);
  return(0);
}

int ESInpPlPrCr()
{
  PlProf *lP,*lJ;
  int i,k;
  double *x,r0;
  
  if((ESFlHole&4)){
    ESHoleReMapRadCData();
    ESHoleReMapInputProfData();
    ESFlHole	&=~4;
  }
  if(ESa0 != 0.){
    ESHoleEdgeBoundary();
    ES3DMetrTnsrCSHole();
  }
  for(i=0; i < NPrP; i++){
    PlPrP[i].io	=i == ESEqSolvInPr ? 1 : 0;
    ioPrP[i]	=PlPrP[i].io;
  }
  lP	=PlPrP+ESEqSolvInPr;

  lP->nd=ESData2Order(lP->xd,lP->yd,lP->kk,NPlPD);

  lJ	=PlCrP+ESEqSolvInCr;
  for(i=0; i < NCrP; i++){
    PlCrP[i].io	=i == ESEqSolvInCr ? 1 : 0;
    ioCrP[i]=PlCrP[i].io;
  }
  lJ->nd	=ESData2Order(lJ->xd,lJ->yd,lJ->kk,NPlPD);

#ifdef SHMEM
#ifdef ASTRA
  if(ESShMem && kAstR){
    int nd;
    double s;
    nd	=lP->nd;
    s	=sqrt(0.5*ESBt/ESgF[ESNa]);
    for(i=0; i < nd; i++){
      lP->xd[i]	=AstA[i]*s;
    }
    lP->xd[0]	=0.;
    memcpy((void*)ESxDtPr,(void*)lP->xd,nd*sizeof(double));
    memcpy((void*)ESxDtCr,(void*)lP->xd,nd*sizeof(double));

    memcpy((void*)lJ->xd,(void*)ESxDtCr,lJ->nd*sizeof(double));
    i	=nAst;
    while(lJ->kk[i]){
      lJ->kk[i]	=0;
      i++;
    }
    lJ->xd[0]	=0.;
    lJ->nd	=nAst;
    i	=nAst-1;
    ESjb[ESNa]	=AstJ[nAst+1];
    if(lJ->xd[i] < 1.){
      lJ->xd[nAst]=1.;
      ESxDtCr[nAst]	=1.;
      ESnDtCr	=nAst;

      ESSetSplA(lJ->xd[i]);
      lJ->yd[nAst]=lJ->yd[i];
      splRA(&s,NULL,ESg22c,ESg22c2);
      lJ->yd[nAst]	*=s;
      splRA(&s,NULL,ESLc,ESLc2);
      lJ->yd[nAst]	*=s;
      splRA(&s,NULL,ESaF,ESaF2a);
      lJ->yd[nAst]	*=s/(ESg22c[ESNa]*ESLc[ESNa]*ESaF[ESNa]);

#ifdef H
      lJ->yd[nAst]=0.2*AstIpl*ESaR0[0]/(ESg22c[ESNa]*ESLc[ESNa]*ESaF[ESNa]);
#endif
      lJ->kk[nAst]	=1;
      lJ->nd	=nAst+1;
      ESjb[ESNa]	=0.;
    }
  }
#endif
#endif/*SHMEM*/

  ESPrData2Spl(lP,lP->kk,lP->xd,lP->yd,ESgaP+1,ESgaP+3);
  ESFirstSetSplDPr();
  ESPrSplData2PrPr(ESEqSolvInPr);

  ESPrData2Spl(lJ,lJ->kk,lJ->xd,lJ->yd,ESgaC+1,ESgaC+3);
  ESFirstSetSplDCr();
  ESPrSplData2CrPr(ESEqSolvInCr);

  if(ESEqSolvInCr == 6){
    double dg0,dg1;
    PlCrP[7].nd	=lJ->nd;
    for(i=0; i < lJ->nd; i++){
      PlCrP[7].kk[i]	=1;
      PlCrP[7].xd[i]	=lJ->xd[i];
    }
    for(i=0; i < lJ->nd; i++){
      PlCrP[7].yd[i]	=lJ->yd[i] > 1e-4 ? 1./lJ->yd[i]
	: 1./(lJ->yd[i]+1e-4);
    }
    ESPrData2Spl(PlCrP+7,lJ->kk,lJ->xd,PlCrP[7].yd,ESgaC+1,ESgaC+3);
    ESPrSplData2CrPr(7);
    for(i=0; i < ESNa1; i++){
      ESsq[i]	=ESgm[i] > 1e-4 ? 1./ESgm[i] 
	: 1./(ESgm1a[i]*0.01*ESsa[ESNa]);
    }
    dg0	=-ESgm1a[0]*ESsq[0]*ESsq[0];
    dg1	=-ESgm1a[ESNa]*ESsq[ESNa]*ESsq[ESNa];
    splAA(ESsq,ESsq1a,ESsq2a,&dg0,&dg1);
  }

  if(ESEqSolvInCr == 7){
    double dg0,dg1;
    for(i=0; i < lJ->nd; i++){
      PlCrP[6].kk[i]	=lJ->kk[i];
      PlCrP[6].xd[i]	=lJ->xd[i];
      PlCrP[6].yd[i]	=1./lJ->yd[i];
    }
    for(i=0; i < ESNa1; i++){
      ESsq[i]	=ESgm[i] > 1e-4 ? 1./ESgm[i] : 1./(ESgm[i]+1e-4);
    }
    dg0	=-ESgm1a[0]*ESsq[0]*ESsq[0];
    dg1	=-ESgm1a[ESNa]*ESsq[ESNa]*ESsq[ESNa];
    splAA(ESsq,ESsq1a,ESsq2a,&dg0,&dg1);
  }
  EC1DEqSolv();
  ESPlProf2Data();
  return(0);
}

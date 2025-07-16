#include "esc.h"

#ifndef mzl_Extrn
/* thise constants are set in zkhL.c (pletzer)
   double EZcr2=0.5,EZcr3,EZcr4=0.25,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
*/

#ifndef blb_esc
#ifndef mzl_Cntr
int ESShMem=0;
void *ESlShMem=NULL;
int *ESShMemI=NULL;
double *ESShMemT=NULL;
double *ESShMemB=NULL;
double *ESShMemP=NULL;
double *ESShMemE=NULL;

unsigned int ESmemlev=0,ESmemlevP=0;
int ESEqSolvFl=0x31,ESEqSolvRadC=0,ESEqSolvInPr=0,ESEqSolvInCr=0
,ESEqSolvInRadC=2;
int FlagPoints=0;

double ESEqSolvTol=1e-3,ESEqTol;
int ESFail=0,ESEqSolvIt=5,ESEqIt;
int ESRecoveryFlag=0;
int ESkNormJ=0;
int ESkNormP=0;
int VMFlag=0;
int ESNMSE=0;
double ESrMSE[256],ESgiMSE[256];

double ESgaG[4]={0.,1e-9,1e-9,1e-9};
double ESgaP[4]={0.,1e-10,1e-15,1e-11};
double ESgaC[4]={0.,1e-10,1e-15,1e-11};

char ESMachineNm[16]=
{'M','a','c','h','i','n','e','N','a','m','e','\0','\0','\0','\0','\0'};
char ESRunID[16]=
{'0','0','0','0','0','0','0','a','0','0','\0','\0','\0','\0','\0','\0'};
char ESMessage[256];
double ESTime=0.; 
#endif

#ifndef mzl_3Dcore
int ESLt=0,ESNt=0,ESNt1=1,ESnAT,ESnAPT,ESLt0,ESNLt,ESNLt1;
double *ESaR0,*ESaZ0,*ESsb,*ESsb1a,*ESsb2a
,*rT1a,*zT1a,*rT1gt,*zT1gt;

double *ESsr,*ESsra,*ESsrt,*ESsrat;
double *ESsz,*ESsza,*ESszt,*ESszat;
double *ESaB,*ESaBa,*ESaBt,*ESaBat;
double *ESgH,*ESgHa,*ESgHt,*ESgHat;

double *gf1T,*csLT,*snLT,*cs1T,*sn1T;
double *rcT,*rcT1a,*rcT2a,*rcT1t,*rsT,*rsT1a,*rsT2a,*rsT1t;
double *dR0T,*dZ0T,*dbT,*dbT1a,*dbT2a;
double *drcT,*drcT1a,*drcT2a,*drsT,*drsT1a,*drsT2a;
#endif

#ifndef mzl_2Dcore
int ESNp=64,ESNp1=65,ESnAP,ESnAF;
int ESFp=32,ESFp1=33;
int *ESnF,*ESmF,*ESkF,ESNf=4,ESNf1=5;
int ESMp=5,ESMp1=6,ES2Mp=10,ES2Mp1=11,ESnMp=0;
double R0,Z0,*ECr,*ECz,*EZvr,*EZvz;
double *EZdra,*EZdrgt,*EZdza,*EZdzgt;
double *ESgh,*ESgt,*EScs1,*ESsn1,*EZcs2,*EZsn2;
double ESaRx,ESaZx;
int ESNsep=0;
double ESRsep[256],ESZsep[256];
int ESiAx=0,ESnBL=1,ESnBL1=2,ESNax,ESNax1;
double ESLx[3],ESVx[3],ESDx[3],ESgFPlVx;
double EcRjet[256],EcZjet[256];
int EcNjet=0;
#endif

#ifndef mzl_1Dcore
int ESNa=20,ESNa1=21;
double *ESsa,*ESqa,*ESpa,*ESp3a,*ECb,*ECb1a,*ECb2a;
double ESa0=0.,ESav=1.;
int ESFlHole=0,ESFlVirtualB=0,ESiIv=0;

int ESBNa=20,ESBNa1=21;
double *ESBsa,*ESBpa,*ESBsb,*ESBsb1a,*ESBsb2a;
#endif

#ifndef mzl_TFC
int TFCLt=12,TFCNt=1,TFCNp=1,TFCNp1=2;
double *xTFC,*yTFC,*zTFC,*gfTFC,*rTFC,ITFC=-1.;
#endif

#ifndef mzl_2DPlPr
double *EStgF,*EStgF1a,*EStgF3a;
double *EStgY,*EStgY1a,*EStgY3a;
double *EStgY1a2t;
#endif

#ifndef mzl_1DPlPr
PlProf PRdC,PlPrP[ESNPRPR],PlCrP[ESNCRPR];
RadCoord RdC[ESNRDCR];

double *ESgF0,*ESgF01a
,*ESgY0,*ESgY01a,*ESgY02a
,*ESgi,*ESgi1a,*ESgi2a
,*ESaY,*ESaY1a,*ESaY2a;
double *ESRC,*ESRC1a,*ESRC2a
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

#ifdef MSE
int ESXNa,ESXNa1,*ESXsxK;
double *ESXsxD,*ESXsrD,*ESXsjD,*ESXspD,*ESXsqD,*ESXsiD;
double *ESXsx,*ESXsr,*ESXsj,*ESXsp,*ESXsq,*ESXsi;
#endif

int ESnDtPr,ESnDtCr;
double ESxDtPr[ESND],ESxDtCr[ESND];
#endif

#ifndef mzl_PlCnst
double ESgFPlV=0,ESgFPlV1a,ESgY2a1a;
double ESRBt;
double ESRext;
double ESBt=1.;
double ESIpl=1.;
double ESgimn;
double ESgimx;
double ESgbN;
double ESgbext;
#endif

#ifndef mzl_3Dgik
double *ESemK,*ESemL,*ESemV,*ESemW,*ESemU1,*ESemU2;
double *ESemK2,*ESemL2,*ESemV2,*ESemW2,*ESemU12,*ESemU22;

double *ESg22c,*ESg22s,*ESg22c1,*ESg22s1,*ESg22c2,*ESg22s2;
double *ESr22,*ESr222,*ESr22c,*ESr22s,*ESr22c2,*ESr22s2;

double *ESg12c,*ESg12s,*ESg12c1,*ESg12s1,*ESg12c2,*ESg12s2
,*ESg11c,*ESg11s,*ESg11c1,*ESg11s1,*ESg11c2,*ESg11s2;

double *ESLc,*ESLs,*ESLc1,*ESLs1,*ESLc2,*ESLs2
,*ESVc,*ESVs,*ESVc1,*ESVs1,*ESVc2,*ESVs2
,*ESDPc,*ESDPs,*ESDPc1,*ESDPs1,*ESDPc2,*ESDPs2;
#endif

#ifndef mzl_SPL
double *ESFaux=NULL,*ESd2Faux;
#ifndef blb_Machine
#ifndef mzl_Tsc2

#ifndef blb_Year
#ifndef mzl_Tsc2
#ifndef blb_Shot
#ifndef mzl_Tsc2
#ifndef blb_Equil
#ifndef mzl_Tsc2
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#ifndef mzl_stab
#ifndef blb_MHD
#ifndef mzl_stab
double *ESGaac,*ESGaas,*ESGaac2,*ESGaas2;
double *ESGaqc,*ESGaqs,*ESGaqc2,*ESGaqs2;
double *ESGqqc,*ESGqqs,*ESGqqc2,*ESGqqs2;
double *ESGazc,*ESGazs,*ESGazc2,*ESGazs2;
double *ESGqzc,*ESGqzs,*ESGqzc2,*ESGqzs2;
double *ESGzzc,*ESGzzs,*ESGzzc2,*ESGzzs2;
double *ESsVc,*ESsVs,*ESsVc2,*ESsVs2;

int ESStN,ESStM1,ESStM2;
int ESStNRec,ESStMRec,ESStKRec;
#endif
#endif
#endif

#endif
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "esc.h"

#ifndef mzl_ESmain
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern unsigned int ESmemlev;
#endif

#ifndef mzl_1Dcore
extern int ESRecoveryFlag;
extern int ESNa,ESNa1;

#ifdef MSE
extern int ESXNa,ESXNa1;
#endif

extern double ESa0;
extern double *ESsa,*ESqa,*ESpa,*ESp3a,*ECb,*ECb1a,*ECb2a;

static int NA=0;
static int NaD=0,naD=0;
#endif

#ifndef mzl_2Dcore
extern int ESNp,ESNp1,ESnAP,ESnAF;
extern int ESFp,ESFp1;
extern int *ESnF,*ESmF,*ESkF,ESNf,ESNf1;
extern int ESMp,ESMp1,ES2Mp,ES2Mp1,ESnMp;
extern double *ESemK,*ESemL,*ESemV,*ESemW,*ESemU1,*ESemU2;
extern double *ESemK2,*ESemL2,*ESemV2,*ESemW2,*ESemU12,*ESemU22;

extern double *ESg22c,*ESg22s,*ESg22c1,*ESg22s1,*ESg22c2,*ESg22s2;
extern double *ESr22,*ESr222,*ESr22c,*ESr22s,*ESr22c2,*ESr22s2;
extern double *ESg12c,*ESg12s,*ESg11c,*ESg11s;
extern double *ESg12c1,*ESg12s1,*ESg11c1,*ESg11s1;
extern double *ESg12c2,*ESg12s2,*ESg11c2,*ESg11s2;
extern double *ESLc,*ESLs,*ESVc,*ESVs,*ESDPc,*ESDPs;
extern double *ESLc1,*ESLs1,*ESVc1,*ESVs1,*ESDPc1,*ESDPs1;
extern double *ESLc2,*ESLs2,*ESVc2,*ESVs2,*ESDPc2,*ESDPs2;

extern double *ESgh,*ESgt,*EScs1,*ESsn1,*EZcs2,*EZsn2;
extern double R0,Z0,*ECr,*ECz,*EZvr,*EZvz;
extern double *EZdra,*EZdrgt,*EZdza,*EZdzgt;
extern double *ESFaux,*ESd2Faux;

extern int ESnBL,ESnBL1,ESNax,ESNax1;

static int NP=0,NpD=-1,NAP=0,NF=0,nF=49;
#endif
extern int ESiIv;

#ifndef mzl_3Dcore
extern int ESLt,ESNt,ESNt1,ESnAT,ESnAPT,ESLt0,ESNLt,ESNLt1;
extern double *ESaR0,*ESaZ0,*ESsb,*ESsb1a,*ESsb2a
,*rT1a,*zT1a,*rT1gt,*zT1gt;

extern double *ESsr,*ESsra,*ESsrt,*ESsrat;
extern double *ESsz,*ESsza,*ESszt,*ESszat;
extern double *ESaB,*ESaBa,*ESaBt,*ESaBat;
extern double *ESgH,*ESgHa,*ESgHt,*ESgHat;

extern double *gf1T,*csLT,*snLT,*cs1T,*sn1T;
extern double *rcT,*rsT,*rcT1a,*rsT1a,*rsT2a,*rcT2a,*rcT1t,*rsT1t;
extern double *dR0T,*dZ0T,*dbT,*dbT1a,*dbT2a;
extern double *drcT,*drcT1a,*drcT2a,*drsT,*drsT1a,*drsT2a;

static int NT=0,LtD=-1,NtD=-1,NAT=0,NAPT=0,NAFT=0;
static int NAFTMp=0,nAFTMp;
#endif

#ifndef mzl_TFC
extern int TFCNt,TFCNp;
#endif


#ifndef mzl_1Dcore
int ESInit1ACore()
{
  int i;
  double h;
  
  h	=(1.-ESa0)/ESNa;
  for(i=0; i < ESNa1; i++){
    ESsa[i]	=ESa0+h*i;
    ESqa[i]	=sqrt(ESsa[i]);
    ESpa[i]	=ESsa[i]*ESsa[i];
    ESp3a[i]	=ESpa[i]*ESsa[i];
  }
  NaD	=ESNa;
  return(0);
}

int ESReInit1ACore()
{
  int n;

  n	=ESNa1+2;
  if(NA < n){
    if(ReInitArray((void**)&ESsa,2*NA,2*n,sizeof(double)) < 0){
      FailureAlarm((char*)ESsa,"ESReInit1ACore() - no memory for ESsa");
      ESexit(0);
    }
    ESqa	=ESsa+n;
    if(ReInitArray((void**)&ESpa,NA,n,sizeof(double)) < 0){
      FailureAlarm((char*)ESpa,"ESReInit1ACore() - no memory for ESpa");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESp3a,NA,n,sizeof(double)) < 0){
      FailureAlarm((char*)ESp3a,"ESReInit1ACore() - no memory for ESp3a");
      ESexit(0);
    }
    if(ReInitArray((void**)&ECb,NA,n,sizeof(double)) < 0){
      FailureAlarm((char*)ECb,"ESReInit1ACore() - no memory for ECb");
      ESexit(0);
    }
    if(ReInitArray((void**)&ECb1a,NA,n,sizeof(double)) < 0){
      FailureAlarm((char*)ECb1a,"ESReInit1ACore() - no memory for ECb1a");
      ESexit(0);
    }
    if(ReInitArray((void**)&ECb2a,NA,n,sizeof(double)) < 0){
      FailureAlarm((char*)ECb2a,"ESReInit1ACore() - no memory for ECb2a");
      ESexit(0);
    }
    NA =n;
    ESmemlev |=0x00000001;
  }
  return(0);
}

int ESDeInit1ACore()
{
  if(NA){
    free(ECb2a);
    free(ECb1a);
    free(ECb);
    free(ESp3a);
    free(ESpa);
    free(ESsa);
    NA=0;
    NaD	=0;
  }
  return(0);
}
#endif

#ifndef mzl_2Dcore
int ESInit1PCore()
{
  int j,k,kk;
  double h,t;

  h	=EZc2gp/ESNp;
  for(j=0; j < ESNp; j++){
    t		=j*h;
    ESgt[j]	=t;
    ESgh[j]	=t;
    EScs1[j]	=cos(t);
    ESsn1[j]	=sin(t);
    t		*=2.;
    EZcs2[j]	=cos(t);
    EZsn2[j]	=sin(t);
  }
  ESgt[j]	=EZc2gp;
  ESgh[j]	=ESgt[j];
  EScs1[j]	=EScs1[0];
  ESsn1[j]	=ESsn1[0];
  EZcs2[j]	=EZcs2[0];
  EZsn2[j]	=EZsn2[0];
  NpD		=ESNp;
  return(0);
}

int ESReInit1PCore()
{
  if(NP < ESNp1){
    int k,kk;
    if(ReInitArray((void**)&ESgh,NP,ESNp1,sizeof(double)) < 0){
      FailureAlarm((char*)ESgh,"ESReInit1PCore() - no memory for ESgh");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESgt,NP,ESNp1,sizeof(double)) < 0){
      FailureAlarm((char*)ESgt,"ESReInit1PCore() - no memory for ESgt");
      ESexit(0);
    }
    if(ReInitArray((void**)&EScs1,NP,ESNp1,sizeof(double)) < 0){
      FailureAlarm((char*)EScs1,"ESReInit1PCore() - no memory for EScs1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESsn1,NP,ESNp1,sizeof(double)) < 0){
      FailureAlarm((char*)ESsn1,"ESReInit1PCore() - no memory for ESsn1");
      ESexit(0);
    }
    if(ReInitArray((void**)&EZcs2,NP,ESNp1,sizeof(double)) < 0){
      FailureAlarm((char*)EZcs2,"ESReInit1PCore() - no memory for EZcs2");
      ESexit(0);
    }
    if(ReInitArray((void**)&EZsn2,NP,ESNp1,sizeof(double)) < 0){
      FailureAlarm((char*)EZsn2,"ESReInit1PCore() - no memory for EZsn2");
      ESexit(0);
    }
    NP =ESNp1;
    ESmemlev |=0x00000002;
  }
  return(0);
}

int ESDeInit1PCore()
{
  if(NP){
    free(EZsn2);
    free(EZcs2);
    free(ESsn1);
    free(EScs1);
    free(ESgt);
    free(ESgh);
    NP	=0;
    NpD	=0;
  }
  return(0);
}
#endif

#ifndef mzl_2Dcore
int ESReInit2PACore()
{
  if(NAP < ESnAP){
    if(ReInitArray((void**)&ECr,NAP,ESnAP,sizeof(double)) < 0){
      FailureAlarm((char*)ECr,"ESReInit2PACore() - no memory for ECr");
      ESexit(0);
    }
    if(ReInitArray((void**)&EZdra,NAP,ESnAP,sizeof(double)) < 0){
      FailureAlarm((char*)EZdra,"ESReInit1ACore() - no memory for EZdra");
      ESexit(0);
    }
    if(ReInitArray((void**)&EZdrgt,NAP,ESnAP,sizeof(double)) < 0){
      FailureAlarm((char*)EZdrgt,"ESReInit1ACore() - no memory for EZdrgt");
      ESexit(0);
    }
    if(ReInitArray((void**)&ECz,NAP,ESnAP,sizeof(double)) < 0){
      FailureAlarm((char*)ECz,"ESReInit2PACore() - no memory for ECz");
      ESexit(0);
    }
    if(ReInitArray((void**)&EZdza,NAP,ESnAP,sizeof(double)) < 0){
      FailureAlarm((char*)EZdza,"ESReInit1ACore() - no memory for EZdza");
      ESexit(0);
    }
    if(ReInitArray((void**)&EZdzgt,NAP,ESnAP,sizeof(double)) < 0){
      FailureAlarm((char*)EZdzgt,"ESReInit1ACore() - no memory for EZdzgt");
      ESexit(0);
    }
    if(ReInitArray((void**)&EZvr,NAP,ESnAP,sizeof(double)) < 0){
      FailureAlarm((char*)EZvr,"ESReInit2PACore() - no memory for EZvr");
      ESexit(0);
    }
    if(ReInitArray((void**)&EZvz,NAP,ESnAP,sizeof(double)) < 0){
      FailureAlarm((char*)EZvz,"ESReInit2PACore() - no memory for EZvz");
      ESexit(0);
    }
    NAP =ESnAP;
    ESmemlev |=0x00000010;
  }
  return(0);
}

int ESDeInit2PACore()
{
  if(NAP){
    free(EZvz);
    free(EZvr);
    free(EZdzgt);
    free(EZdza);
    free(ECz);
    free(EZdrgt);
    free(EZdra);
    free(ECr);
    NAP=0;
  }
  return(0);
}
#endif


int ESReInit1FCoreInt()
{
  if(NF < nF){
    int i;
    if(ReInitArray((void**)&ESnF,NF,nF,sizeof(int)) < 0){
      FailureAlarm((char*)ESnF,"ESReInit1FCoreInt() - no memory for ESnF");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESmF,NF,nF,sizeof(int)) < 0){
      FailureAlarm((char*)ESmF,"ESReInit1FCoreInt() - no memory for ESmF");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESkF,NF,nF,sizeof(int)) < 0){
      FailureAlarm((char*)ESkF,"ESReInit1FCoreInt() - no memory for ESkF");
      ESexit(0);
    }
    if(NF == 0){
      ESInit1FCoreInt();
    }
    else{
      for(i=NF; i < nF; i++){
	ESkF[i]	=0;
	ESmF[i]	=i;
	ESnF[i]	=10;
      }
    }
    NF =nF;
    ESmemlev |=0x00001000;
  }
  return(0);
}

int ESDeInit1FCoreInt()
{
  if(NF){
    free(ESkF);
    free(ESmF);
    free(ESnF);
    NF=0;
  }
  return(0);
}

int ESCheck1FGeomInt()
{
  int i,j;
  
  ESNf1=0;
  while(ESNf1 < nF){
    if(ESkF[ESNf1]){
      i=0;
      while(i < ESNf1){
	if(ESnF[i] == ESnF[ESNf1] && ESmF[i] == ESmF[ESNf1]){
	  ESkF[ESNf1]	=0;
	  break;
	}
	i++;
      }
    }
    if(ESkF[ESNf1] == 0){
      i=ESNf1+1;
      while(i < nF){
	if(ESkF[i]){
	  ESkF[ESNf1]	=ESkF[i];
	  ESkF[i]	=0;
	  break;
	}
	i++;
      }
    }
    if(ESkF[ESNf1]){
      ESNf1++;
    }
    else{
      break;
    }
  }
  ESNf=ESNf1-1;
  return(0);
}

int ESInit1FCoreInt()
{
  int j;

  for(j=0; j < ESNf1; j++){
    ESnF[j]	=0;
    ESmF[j]	=j;
    ESkF[j]	=1;
  }
  for(; j < nF; j++){
    ESnF[j]	=3;
    ESmF[j]	=j;
    ESkF[j]	=0;
  }
  return(0);
}

int ESInit1FEqWaveNum()
{
  int j,n,jn;
  
  if(ESLt == 0){
    return(0);
  }
  jn	=ESNf1;
  for(j=0; j < ESNf1 && jn < nF; j++){
    if(ESnF[j]	== 0 && ESmF[j]){
      ESmF[jn]	=ESmF[j];
      ESnF[jn]	=ESLt*ESmF[j];
      jn++;
    }
  }
  return(0);
}
#ifndef mzl_3Dcore
int ESInit1TCore()
{
  int n;
  double h;

  if(ESNt){
    h	=EZc2gp/ESNt;
    for(n=0; n < ESNt; n++){
      gf1T[n]	=n*h;
      csLT[n]	=cos(ESLt*gf1T[n]);
      snLT[n]	=sin(ESLt*gf1T[n]);
      cs1T[n]	=cos(gf1T[n]);
      sn1T[n]	=sin(gf1T[n]);
    }
    gf1T[n]	=EZc2gp;
    csLT[n]	=csLT[0];
    snLT[n]	=snLT[0];
    cs1T[n]	=cs1T[0];
    sn1T[n]	=sn1T[0];
  }
  else{
    gf1T[0]	=0.;
    csLT[0]	=1.;
    snLT[0]	=0.;
    cs1T[0]	=1.;
    sn1T[0]	=0.;
  }
  NtD	=ESNt;
  LtD	=ESLt;
  return(0);
}

int ESReInit1TCore()
{
  if(NT < ESNt1){
    if(ReInitArray((void**)&ESaR0,NT,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)ESaR0,"ESReInit1TCore() - no memory for ESaR0");
      ESexit(0);
    }
    if(ReInitArray((void**)&dR0T,NT,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)dR0T,"ESReInit1TCore() - no memory for dR0T");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESaZ0,NT,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)ESaZ0,"ESReInit1TCore() - no memory for ESaZ0");
      ESexit(0);
    }
    if(ReInitArray((void**)&dZ0T,NT,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)dZ0T,"ESReInit1TCore() - no memory for dZ0T");
      ESexit(0);
    }
    if(ReInitArray((void**)&gf1T,NT,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)gf1T,"ESReInit1TCore() - no memory for gf1T");
      ESexit(0);
    }
    if(ReInitArray((void**)&csLT,NT,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)csLT,"ESReInit1TCore() - no memory for csLT");
      ESexit(0);
    }
    if(ReInitArray((void**)&snLT,NT,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)snLT,"ESReInit1TCore() - no memory for snLT");
      ESexit(0);
    }
    if(ReInitArray((void**)&cs1T,NT,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)cs1T,"ESReInit1TCore() - no memory for cs1T");
      ESexit(0);
    }
    if(ReInitArray((void**)&sn1T,NT,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)sn1T,"ESReInit1TCore() - no memory for sn1T");
      ESexit(0);
    }
    NT =ESNt1;
    ESmemlev |=0x00000004;
  }
  return(0);
}

int ESDeInit1TCore()
{
  if(NT){
    free(sn1T);
    free(cs1T);
    free(snLT);
    free(csLT);
    free(gf1T);
    free(dZ0T);
    free(ESaZ0);
    free(dR0T);
    free(ESaR0);
    NT=0;
  }
  return(0);
}

int ESReInit2ATCore()
{
  int n;
  
  n	=ESNt1*(ESNa1+2);
  if(NAT < n){
    if(ReInitArray((void**)&ESsb,NAT,n,sizeof(double)) < 0){
      FailureAlarm((char*)ESsb,"ESReInit2ATCore() - no memory for ESsb");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESsb1a,NAT,n,sizeof(double)) < 0){
      FailureAlarm((char*)ESsb1a,"ESReInit2ATCore() - no memory for ESsb1a");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESsb2a,NAT,n,sizeof(double)) < 0){
      FailureAlarm((char*)ESsb2a,"ESReInit2ATCore() - no memory for ESsb2a");
      ESexit(0);
    }
    if(ReInitArray((void**)&dbT,NAT,n,sizeof(double)) < 0){
      FailureAlarm((char*)dbT,"ESReInit2ATCore() - no memory for dbT");
      ESexit(0);
    }
    if(ReInitArray((void**)&dbT1a,NAT,n,sizeof(double)) < 0){
      FailureAlarm((char*)dbT1a,"ESReInit2ATCore() - no memory for dbT1a");
      ESexit(0);
    }
    if(ReInitArray((void**)&dbT2a,NAT,n,sizeof(double)) < 0){
      FailureAlarm((char*)dbT2a,"ESReInit2ATCore() - no memory for dbT2a");
      ESexit(0);
    }
    NAT =n;
    ESmemlev |=0x00000008;
  }
  return(0);
}

int ESDeInit2ATCore()
{
  if(NAT){
    free(dbT2a);
    free(dbT1a);
    free(dbT);
    free(ESsb2a);
    free(ESsb1a);
    free(ESsb);
    NAT=0;
  }
  return(0);
}
#endif
#ifndef mzl_2Dcore
int ESReInit3PATCore()
{
  int i,j;
  if(NAPT < ESnAPT){
    if(NAPT){
      free(ESsr);
    }
    ESsr	=(double*)malloc(8*ESnAPT*sizeof(double));
    ESsra	=ESsr+ESnAPT;
    ESsrt	=ESsra+ESnAPT;
    ESsrat	=ESsrt+ESnAPT;
    ESaB	=ESsrat+ESnAPT;
    ESaBa	=ESaB+ESnAPT;
    ESaBt	=ESaBa+ESnAPT;
    ESaBat	=ESaBt+ESnAPT;

    ESsz	=(double*)malloc(8*ESnAPT*sizeof(double));
    ESsza	=ESsz+ESnAPT;
    ESszt	=ESsza+ESnAPT;
    ESszat	=ESszt+ESnAPT;
    ESgH	=ESszat+ESnAPT;
    ESgHa	=ESgH+ESnAPT;
    ESgHt	=ESgHa+ESnAPT;
    ESgHat	=ESgHt+ESnAPT;

    if(ReInitArray((void**)&rT1a,NAPT,ESnAPT,sizeof(double)) < 0){
      FailureAlarm((char*)rT1a,"ESReInit3PATCore() - no memory for rT1a");
      ESexit(0);
    }
    if(ReInitArray((void**)&rT1gt,NAPT,ESnAPT,sizeof(double)) < 0){
      FailureAlarm((char*)rT1gt,"ESReInit3PATCore() - no memory for rT1gt");
      ESexit(0);
    }
    if(ReInitArray((void**)&zT1a,NAPT,ESnAPT,sizeof(double)) < 0){
      FailureAlarm((char*)zT1a,"ESReInit3PATCore() - no memory for zT1a");
      ESexit(0);
    }
    if(ReInitArray((void**)&zT1gt,NAPT,ESnAPT,sizeof(double)) < 0){
      FailureAlarm((char*)zT1gt,"ESReInit3PATCore() - no memory for zT1gt");
      ESexit(0);
    }
    NAPT	=ESnAPT;
    ESmemlev |=0x00000020;
  }
  j	=ESFp1*(ESNa1+2)*ESNt1;
  ESnAPT=ESNa1*ESNp1*ESNt1;
  if(NAFT < j){
    if(ReInitArray((void**)&rcT,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)rcT,"ESReInit3PATCore() - no memory for rcT");
      ESexit(0);
    }
    if(ReInitArray((void**)&rcT1a,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)rcT1a,"ESReInit3PATCore() - no memory for rcT1a");
      ESexit(0);
    }
    if(ReInitArray((void**)&rcT2a,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)rcT2a,"ESReInit3PATCore() - no memory for rcT2a");
      ESexit(0);
    }
    if(ReInitArray((void**)&rcT1t,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)rcT1t,"ESReInit3PATCore() - no memory for rcT1t");
      ESexit(0);
    }
    if(ReInitArray((void**)&drcT,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)drcT,"ESReInit3PATCore() - no memory for drcT");
      ESexit(0);
    }
    if(ReInitArray((void**)&drcT1a,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)drcT1a,"ESReInit3PATCore() - no memory for drcT1a");
      ESexit(0);
    }
    if(ReInitArray((void**)&drcT2a,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)drcT2a,"ESReInit3PATCore() - no memory for drcT2a");
      ESexit(0);
    }
    if(ReInitArray((void**)&rsT,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)rsT,"ESReInit3PATCore() - no memory for rsT");
      ESexit(0);
    }
    if(ReInitArray((void**)&rsT1a,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)rsT1a,"ESReInit3PATCore() - no memory for rsT1a");
      ESexit(0);
    }
    if(ReInitArray((void**)&rsT2a,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)rsT2a,"ESReInit3PATCore() - no memory for rsT2a");
      ESexit(0);
    }
    if(ReInitArray((void**)&rsT1t,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)rsT1t,"ESReInit3PATCore() - no memory for rsT1t");
      ESexit(0);
    }
    if(ReInitArray((void**)&drsT,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)drsT,"ESReInit3PATCore() - no memory for drsT");
      ESexit(0);
    }
    if(ReInitArray((void**)&drsT1a,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)drsT1a,"ESReInit3PATCore() - no memory for drsT1a");
      ESexit(0);
    }
    if(ReInitArray((void**)&drsT2a,NAFT,j,sizeof(double)) < 0){
      FailureAlarm((char*)drsT2a,"ESReInit3PATCore() - no memory for drsT2a");
      ESexit(0);
    }
    for(i=NAFT; i < j; i++){
      rcT[i]	=0.;
      rcT1a[i]	=0.;
      rcT2a[i]	=0.;
      drcT[i]	=0.;
      drcT1a[i]	=0.;
      drcT2a[i]	=0.;
      rcT1t[i]	=0.;
      rsT[i]	=0.;
      rsT1a[i]	=0.;
      rsT2a[i]	=0.;
      drsT[i]	=0.;
      drsT1a[i]	=0.;
      drsT2a[i]	=0.;
      rsT1t[i]	=0.;
    }
    NAFT	=j;
    ESmemlev |=0x00000020;
  }
  return(0);
}

int ESDeInit3PATCore()
{
  if(NAFT){
    free(drsT2a);
    free(drsT1a);
    free(drsT);
    free(rsT1t);
    free(rsT2a);
    free(rsT1a);
    free(rsT);
    free(drcT2a);
    free(drcT1a);
    free(drcT);
    free(rcT1t);
    free(rcT2a);
    free(rcT1a);
    free(rcT);
    NAFT	=0;
  }
  if(NAPT){
    free(zT1gt);
    free(zT1a);
    free(rT1gt);
    free(rT1a);
    free(ESsz);
    free(ESsr);
    NAPT=0;
  }
  return(0);
}
#endif

int ESReInitMetrTnsrCS()
{
  nAFTMp=ESNt1*ESNa1*ES2Mp1;
  if(NAFTMp < nAFTMp){
    int i;
    i	=nAFTMp*ES2Mp1;
    if(NAFTMp){
      free(ESr22);
    }
    ESr22	=(double*)malloc(2*i*sizeof(double));
    if(ESr22 == NULL){
      FailureAlarm((char*)ESr22,"ESReInitMetrTnsrCS() - no memory for ESr22");
      ESexit(0);
    }
    ESr222	=ESr22+i;
    if(ReInitArray((void**)&ESr22c,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESr22c,"ESReInitMetrTnsrCS()- no memory for ESr22c");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESr22c2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESr22c2,"ESReInitMetrTnsrCS-no memory for ESr22c2");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESr22s,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESr22s,"ESReInitMetrTnsrCS - no memory for ESr22s");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESr22s2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESr22s2,"ESReInitMetrTnsrCS- no memory for ESr22s2");
      ESexit(0);
    }

    if(ReInitArray((void**)&ESg22c,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg22c,"ESReInitMetrTnsrCS()- no memory for ESg22c");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg22c1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg22c1,"ESReInitMetrTnsrCS-no memory for ESg22c1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg22c2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg22c2,"ESReInitMetrTnsrCS-no memory for ESg22c2");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg22s,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg22s,"ESReInitMetrTnsrCS - no memory for ESg22s");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg22s1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg22s1,"ESReInitMetrTnsrCS- no memory for ESg22s1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg22s2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg22s2,"ESReInitMetrTnsrCS- no memory for ESg22s2");
      ESexit(0);
    }

    if(ReInitArray((void**)&ESg12c,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg12c,"ESReInitMetrTnsrCS - no memory for ESg12c");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg12c1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg12c1,"ESReInitMetrTnsrCS -no memory for ESg12c1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg12c2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg12c2,"ESReInitMetrTnsrCS- no memory for ESg12c2");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg12s,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg12s,"ESReInitMetrTnsrCS - no memory for ESg12s");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg12s1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg12s1,"ESReInitMetrTnsrCS- no memory for ESg12s1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg12s2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg12s2,"ESReInitMetrTnsrCS- no memory for ESg12s2");
      ESexit(0);
    }

    if(ReInitArray((void**)&ESg11c,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg11c,"ESReInitMetrTnsrCS - no memory for ESg11c");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg11c1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg11c1,"ESReInitMetrTnsrCS- no memory for ESg11c1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg11c2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg11c2,"ESReInitMetrTnsrCS- no memory for ESg11c2");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg11s,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg11s,"ESReInitMetrTnsrCS- no memory for ESg11s");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg11s1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg11s1,"ESReInitMetrTnsrCS- no memory for ESg11s1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESg11s2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESg11s2,"ESReInitMetrTnsrCS- no memory for ESg11s2");
      ESexit(0);
    }

    if(ReInitArray((void**)&ESLc,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESLc,"ESReInitMetrTnsrCS() - no memory for ESLc");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESLc1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESLc1,"ESReInitMetrTnsrCS() - no memory for ESLc1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESLc2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESLc2,"ESReInitMetrTnsrCS() - no memory for ESLc2");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESLs,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESLs,"ESReInitMetrTnsrCS() - no memory for ESLs");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESLs1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESLs1,"ESReInitMetrTnsrCS() - no memory for ESLs1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESLs2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESLs2,"ESReInitMetrTnsrCS() - no memory for ESLs2");
      ESexit(0);
    }

    if(ReInitArray((void**)&ESVc,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESVc,"ESReInitMetrTnsrCS() - no memory for ESVc");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESVc1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESVc1,"ESReInitMetrTnsrCS() - no memory for ESVc1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESVc2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESVc2,"ESReInitMetrTnsrCS() - no memory for ESVc2");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESVs,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESVs,"ESReInitMetrTnsrCS() - no memory for ESVs");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESVs1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESVs1,"ESReInitMetrTnsrCS() - no memory for ESVs1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESVs2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESVs2,"ESReInitMetrTnsrCS() - no memory for ESVs2");
      ESexit(0);
    }

    if(ReInitArray((void**)&ESDPc,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESDPc,"ESReInitMetrTnsrCS() - no memory for ESDPc");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESDPc2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESDPc2,"ESReInitMetrTnsrCS()- no memory for ESDPc2");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESDPc1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESDPc1,"ESReInitMetrTnsrCS()- no memory for ESDPc1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESDPs,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESDPs,"ESReInitMetrTnsrCS() - no memory for ESDPs");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESDPs1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESDPs1,"ESReInitMetrTnsrCS()- no memory for ESDPs1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESDPs2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESDPs2,"ESReInitMetrTnsrCS()- no memory for ESDPs2");
      ESexit(0);
    }
#ifdef H
    if(ReInitArray((void**)&ESemK,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESemK,"ESReInitMetrTnsrCS() - no memory for ESemK");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESemL,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESemL,"ESReInitMetrTnsrCS() - no memory for ESemL");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESemV,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESemV,"ESReInitMetrTnsrCS() - no memory for ESemV");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESemW,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESemW,"ESReInitMetrTnsrCS() - no memory for ESemW");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESemU1,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESemU1,"ESReInitMetrTnsrCS() - no memory for ESemU1");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESemU2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESemU2,"ESReInitMetrTnsrCS() - no memory for ESemU2");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESemK2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESemK2,"ESReInitMetrTnsrCS() - no memory for ESemK2");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESemL2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESemL2,"ESReInitMetrTnsrCS() - no memory for ESemL2");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESemV2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESemV2,"ESReInitMetrTnsrCS() - no memory for ESemV2");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESemW2,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESemW2,"ESReInitMetrTnsrCS()- no memory for ESemW2");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESemU12,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESemU12,"ESReInitMetrTnsrCS()- no memory for ESemU12");
      ESexit(0);
    }
    if(ReInitArray((void**)&ESemU22,NAFTMp,nAFTMp,sizeof(double)) < 0){
      FailureAlarm((char*)ESemU22,"ESReInitMetrTnsrCS() - no memory for ESemU22");
      ESexit(0);
    }
#endif
    NAFTMp	=nAFTMp;
    ESmemlev |=0x02000000;
  }
  return(0);
}

int ESDeInitMetrTnsrCS()
{
  if(NAFTMp){
    free(ESDPs2);
    free(ESDPs1);
    free(ESDPs);
    free(ESDPc2);
    free(ESDPc1);
    free(ESDPc);
    free(ESVs2);
    free(ESVs1);
    free(ESVs);
    free(ESVc2);
    free(ESVc1);
    free(ESVc);
    free(ESLs2);
    free(ESLs1);
    free(ESLs);
    free(ESLc2);
    free(ESLc1);
    free(ESLc);
    free(ESg11s2);
    free(ESg11s1);
    free(ESg11s);
    free(ESg11c2);
    free(ESg11c1);
    free(ESg11c);
    free(ESg12s2);
    free(ESg12s1);
    free(ESg12s);
    free(ESg12c2);
    free(ESg12c1);
    free(ESg12c);
    free(ESg22s2);
    free(ESg22s1);
    free(ESg22s);
    free(ESg22c2);
    free(ESg22c1);
    free(ESg22c);

    free(ESr22s2);
    free(ESr22s);
    free(ESr22c2);
    free(ESr22c);
    free(ESr22);
#ifdef H
    free(ESemU22);
    free(ESemU12);
    free(ESemW2);
    free(ESemV2);
    free(ESemL2);
    free(ESemK2);
    free(ESemU2);
    free(ESemU1);
    free(ESemW);
    free(ESemV);
    free(ESemL);
    free(ESemK);
#endif
    NAFTMp=0;
  }
  return(0);
}

int ESCheckIntGeom()
{
  ESLt0	=abs(ESLt);
  if(ESLt == 0){
    ESNt	=0;
    ESNLt	=1;
    ESNLt1	=1;
  }
  else{
    ESNLt	=ESNt/ESLt0;
    ESNLt1	=ESNLt+1;
  }
  if(ESNt < ESLt0){
    ESNt=ESLt0;
  }
  if(ESNa <= 0){
    ESNa=10;
  }
  if(ESNp <= 0){
    ESNp=12;
  }
  return(0);
}

int ESInitPlGeomInt()
{
#ifdef PFC
  extern int ESnLm,nPFB;
  extern void readefn(char *NmOut,char *Nm);
  static int kVa=0;

  int kk[2]={0,0};
  char str[4][32]={"./","./","mhdin.dat","input"};
  char *NmD[2],*NmF[2];
  char Nm1[64],Nm2[64];
  char Dir[64],Nm[64];

  NmD[0]	=str[0];
  NmD[1]	=str[1];
  NmF[0]	=str[2];
  NmF[1]	=str[3];
#endif

#ifndef Tbl_Machine
  ESCheckIntGeom();
  if(NF < nF){
    ESReInit1FCoreInt();
  }
  ESCheck1FGeomInt();
  ESNt1	=ESNt+1;
  if(NT < ESNt1){
    ESReInit1TCore();
  }
  if(NtD != ESNt || LtD != ESLt){
    ESInit1TCore();
    ESInit1FEqWaveNum();
  }
  ESNa1	=ESNa+1;
  ESiIv	=ESNa;

#ifdef MSE
  ESXNa	=2*ESNa;
  ESXNa1=ESXNa+1;
#endif

  ESNax	=ESNa+ESnBL;
  ESNax1=ESNax+1;
  if(NA < ESNa1){
    ESReInit1ACore();
  }
  if(NaD != ESNa){
    ESInit1ACore();
  }
  ESNp1	=ESNp+1;
  if(ESMp > ESNp/4){
    ESMp	=ESNp/4;
  }
  ESFp	=2*ESMp;

  ESFp	=ESNp/2;

  ESMp1	=ESMp+1;
  ES2Mp	=2*ESMp;
  ES2Mp1=ES2Mp+1;
  ESFp1	=ESFp+1;
  ESnAF	=ESFp1*ESNa1;
  if(NP < ESNp1){
    ESReInit1PCore();
  }
  if(NpD != ESNp){
    ESInit1PCore();
  }

#ifdef PFC
  if(kk[0] || kk[1]){
    char *lD,*lN;
    CbGetCurrentCbName(&lD,&lN);
    sprintf(Dir,"%s/",lD);
  }
  if(kk[0]){
    sprintf(Nm1,"%s%s",NmD[0],NmF[0]);
    efitpfc_(&kVa,Dir,Nm1);
    kk[0]	=0;
  }
  if(kk[1]){
    sprintf(Nm2,"%s%s",NmD[1],NmF[1]);
    efitlim_(Dir,Nm2);
    kk[1]	=0;
  }
#endif
#endif/*Tbl_Machine*/
  ESnAT	=ESNa1*ESNt1;
  ESReInit2ATCore();
  ESnAP	=ESNa1*ESNp1;
  if(NAP < ESnAP){
    ESReInit2PACore();
  }
  ESnAPT=ESnAP*ESNt1;
  ESReInit3PATCore();
  ESReInit3APTSpl();
  ESInitSplAPT();
  return(0);
}

int ESDeInitEsc()
{
  if((ESmemlev&0x80000000)){
    ESmemlev -=0x80000000;
  }
  if((ESmemlev&0x40000000)){
    ESmemlev -=0x40000000;
  }
  if((ESmemlev&0x20000000)){
    ESmemlev -=0x20000000;
  }
  if((ESmemlev&0x10000000)){
    ESDeInitMapSpl();
    ESmemlev -=0x10000000;
  }
  if((ESmemlev&0x08000000)){
    ESDeInitMHDStabLowN();
    ESmemlev -=0x08000000;
  }
  if((ESmemlev&0x04000000)){
    ESDeInitMHDMetrTnsr();
    ESmemlev -=0x04000000;
  }
  if((ESmemlev&0x02000000)){
    ESDeInitMetrTnsrCS();
    ESmemlev -=0x02000000;
  }
  if((ESmemlev&0x01000000)){
    ESDeInitPlProf();
    ESmemlev -=0x01000000;
  }
  if((ESmemlev&0x00800000)){
    DeInitf2splT();
    DeInitf2splP();
    DeInitf2splA();
    ESmemlev -=0x00800000;
  }
  if((ESmemlev&0x004000000)){
    ESDeInitgXcs();
    ESmemlev -=0x00400000;
  }
  if((ESmemlev&0x00200000)){
    ESDeInit2DAPrC();
    ESmemlev -=0x00200000;
  }
  if((ESmemlev&0x00100000)){
    ESDeInit3APTSpl();
    ESmemlev -=0x00100000;
  }
  if((ESmemlev&0x00080000)){
    ESDeInitEnvTFC();
    ESmemlev -=0x00080000;
  }
  if((ESmemlev&0x00040000)){
    ESDeInitDmTFC();
    ESmemlev -=0x00040000;
  }
  if((ESmemlev&0x00020000)){
    ESDeInitPncr();
    ESmemlev -=0x00020000;
  }
  if((ESmemlev&0x00010000)){
    ESDeInitLsode();
    ESmemlev -=0x00010000;
  }
  if((ESmemlev&0x00008000)){
    ESDeInit1LTFC();
    ESmemlev -=0x00008000;
  }
  if((ESmemlev&0x00004000)){
    ESDeInit1DTFC();
    ESmemlev -=0x00004000;
  }
  if((ESmemlev&0x00002000)){
    ESDeInit3DTFC();
    ESmemlev -=0x00002000;
  }
  if((ESmemlev&0x00001000)){
    ESDeInit1FCoreInt();
    ESmemlev -=0x00001000;
  }
  if((ESmemlev&0x00000800)){
    ESDeInitPlVeq();
    ESmemlev -=0x00000800;
  }
  if((ESmemlev&0x00000400)){
    ESDeInit2DFPPlV();
    ESmemlev -=0x00000400;
  }
  if((ESmemlev&0x00000200)){
    ESDeInit2DFTPlV();
    ESmemlev -=0x00000200;
  }
  if((ESmemlev&0x00000100)){
    ESDeInit1DFPlV();
    ESmemlev -=0x00000100;
  }
  if((ESmemlev&0x00000080)){
    ESDeInit2DTPlV();
    ESmemlev -=0x00000080;
  }
  if((ESmemlev&0x00000040)){
    ESDeInit1DPlV();
    ESmemlev -=0x00000040;
  }
  if((ESmemlev&0x00000020)){
    ESDeInit3PATCore();
    ESmemlev -=0x00000020;
  }
  if((ESmemlev&0x00000010)){
    ESDeInit2PACore();
    ESmemlev -=0x00000010;
  }
  if((ESmemlev&0x00000008)){
    ESDeInit2ATCore();
    ESmemlev -=0x00000008;
  }
  if((ESmemlev&0x00000004)){
    ESDeInit1TCore();
    ESmemlev -=0x00000004;
  }
  if((ESmemlev&0x00000002)){
    ESDeInit1PCore();
    ESmemlev -=0x00000002;
  }
  if((ESmemlev&0x00000001)){
    ESDeInit1ACore();
    ESmemlev -=0x00000001;
  }
  if(ESFaux != NULL){
    free(ESFaux);
    ESFaux	=NULL;
  }
  return(0);
}

int ESSaveESC(double time)
{
  FILE *fpb,*fpa;
  static int Fl	=0;
 
  fpa=fopen("ESState.ES","a");
  fpb=fopen("ESState.es","a");

  if(fpa == NULL || fpb == NULL){
    if(fpa != NULL) fclose(fpa);
    if(fpb != NULL) fclose(fpb);
    fpa=fopen("ESState.ES","w");
    fpb=fopen("ESState.es","w");
  }

#ifdef H
  if(Fl){
    fpa=fopen("ESState.ES","a");
    fpb=fopen("ESState.es","a");
  }
  else{
    fpa=fopen("ESState.ES","w");
    fpb=fopen("ESState.es","w");
  }
#endif
  fprintf(fpa,"--- time[sec]=%12.5ea ---------------------\n",time);
  fwrite(&ESNa,sizeof(int),1,fpb);
  fwrite(&ESNp,sizeof(int),1,fpb);

  fclose(fpa);
  fclose(fpb);
  return(0);
}

int ESRecoverESC(double time)
{
  FILE *fpb,*fpa;
  fpa=fopen("ESState.ES","r");
  fpb=fopen("ESState.es","r");
  fclose(fpa);
  fclose(fpb);
  return(0);
}

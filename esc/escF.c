#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "esc.h"

#ifndef mzl_ESmain
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern int ESmemlev;
extern double difftime(time_t, time_t);
extern clock_t clock(void);

static time_t sTime1,sTime2;
#endif

#ifndef mzl_2Dcore
extern int ESNp,ESNp1;
extern int *ESnF,*ESmF,*ESkF,ESNf,ESNf1;
#endif

#ifndef mzl_3Dcore
extern int ESLt,ESNt,ESNt1;
extern double *gf1T,*csLT,*snLT,*cs1T,*sn1T;
#endif

#ifndef mzl_TFC
extern int TFCLt,TFCNt,TFCNp,TFCNp1;
extern double *xTFC,*yTFC,*zTFC,*gfTFC,*rTFC,ITFC;

static double *gfL,*csL,*snL;
static int nTP,TFCNt1;
static double Rmn,Rmx,Zmn,Zmx,*gf0TFC;
static double *x1TFC=NULL,*x2TFC,*y1TFC,*y2TFC,*r1TFC,*r2TFC,*gf1TFC,*gf2TFC;

static int *jRmnTFC,*jRmxTFC,*jZmnTFC,*jZmxTFC,NDm=0;

static double *x1Env,*x2Env,*y1Env,*y2Env,*zEnv,*r1Env,*r2Env;
static double drSpl,dtSpl,dzSpl;

static double rCell,r1Cell,tCell,t1Cell,zCell,z1Cell,rhrCell,rhtCell,rhzCell;
int irCell,ktCell,jzCell,
J000=0,J100=0,J010=0,J001=0,J110=0,J101=0,J011=0,J111=0;
static double Pr0,Pr1,Pr2,Pr3,dPr0,dPr1,dPr2,dPr3
,Pt0,Pt1,Pt2,Pt3,dPt0,dPt1,dPt2,dPt3
,Pz0,Pz1,Pz2,Pz3,dPz0,dPz1,dPz2,dPz3;

static int SplNt=0,*SplNz,SplN=0;

static int **Nr1Spl,**Nr2Spl,*Nz1Spl,*Nz2Spl
,NrSpl=0,NzSpl=0,nrSpl=20,nzSpl=20,nrSpl1,nzSpl1;

static int ***k2kSpl;
static double *xSpl,*ySpl,*zSpl,*rSpl;
static double *AR,*ARr,*ARt,*ARz,*ARrt,*ARrz,*ARtz,*ARrtz;
static double *AT,*ATr,*ATt,*ATz,*ATrt,*ATrz,*ATtz,*ATrtz;
static double *AZ,*AZr,*AZt,*AZz,*AZrt,*AZrz,*AZtz,*AZrtz;
static double *BR,*BRr,*BRt,*BRz,*BRrt,*BRrz,*BRtz,*BRrtz;
static double *BT,*BTr,*BTt,*BTz,*BTrt,*BTrz,*BTtz,*BTrtz;
static double *BZ,*BZr,*BZt,*BZz,*BZrt,*BZrz,*BZtz,*BZrtz;

#endif

#ifndef mzl_SPL
int ESInitDmTFC()
{
  int j,n,jn;

  jn	=0;
  for(n=0; n < TFCNt; n++){
    Rmn	=rTFC[jn];
    Rmx	=rTFC[jn];
    Zmn	=zTFC[jn];
    Zmx	=zTFC[jn];
    jRmnTFC[n]	=0;
    jRmxTFC[n]	=0;
    jZmnTFC[n]	=0;
    jZmxTFC[n]	=0;
    for(j=0; j < TFCNp1; j++){
      if(Rmn > rTFC[jn]){
	Rmn=rTFC[jn];
	jRmnTFC[n]	=j;
      }
      if(Rmx < rTFC[jn]){
	Rmx=rTFC[jn];
	jRmxTFC[n]	=j;
      }
      if(Zmn > zTFC[jn]){
	Zmn=zTFC[jn];
	jZmnTFC[n]	=j;
      }
      if(Zmx < zTFC[jn]){
	Zmx=zTFC[jn];
	jZmxTFC[n]	=j;
      }
      jn++;
    }
  }
  return(0);
}

int ESReInitDmTFC()
{
  if(NDm < TFCNt){
    if(ReInitArray((void**)&jRmnTFC,NDm,TFCNt,sizeof(int)) < 0){
      FailureAlarm((char*)jRmnTFC,"ESReInitDmTFC() - no memory for jRmnTFC");
      ESexit(0);
    }
    if(ReInitArray((void**)&jRmxTFC,NDm,TFCNt,sizeof(int)) < 0){
      FailureAlarm((char*)jRmxTFC,"ESReInitDmTFC() - no memory for jRmxTFC");
      ESexit(0);
    }
    if(ReInitArray((void**)&jZmnTFC,NDm,TFCNt,sizeof(int)) < 0){
      FailureAlarm((char*)jZmnTFC,"ESReInitDmTFC() - no memory for jZmnTFC");
      ESexit(0);
    }
    if(ReInitArray((void**)&jZmxTFC,NDm,TFCNt,sizeof(int)) < 0){
      FailureAlarm((char*)jZmxTFC,"ESReInitDmTFC() - no memory for jZmxTFC");
      ESexit(0);
    }
    NDm =TFCNt;
    ESmemlev |=0x00040000;
  }
  return(0);
}

int ESDeInitDmTFC()
{
  if(NDm){
    free(jZmxTFC);
    free(jZmnTFC);
    free(jRmxTFC);
    free(jRmnTFC);
    NDm=0;
  }
  return(0);
}

int ESInitEnvTFC(double z)
{
  int j1,j2,j,jj,n,jn,jz[8];
  double *pz,s,csL,snL;
  double r1,EZr2,gf;

  for(n=0; n < TFCNt; n++){
    jn	=TFCNp1*n;
    jj	=0;
    j		=jn+jZmnTFC[n];
    if(z <= zTFC[j]){
      jj=1;
    }
    j	=jn+jZmxTFC[n];
    if(zTFC[j] <= z){
      jj=1;
    }
    if(jj){
      x1TFC[n]	=xTFC[j];
      x2TFC[n]	=xTFC[j];
      y1TFC[n]	=yTFC[j];
      y2TFC[n]	=yTFC[j];
    }    
    else{
      pz	=zTFC+jn;
      jj	=0;
      for(j=0; j < TFCNp; j++){
	j1=j+1;
	if(pz[j] != pz[j1]){
	  if((z-pz[j] <= 0. && z-pz[j1] > 0.) || 
	     (z-pz[j] >= 0. && z-pz[j1] < 0.)){
	    jz[jj]	=j;
	    jj++;
	    if(jj > 7){
	      printf("Too many intersections in TFC %d%c\n",n,7);
	      jj=7;
	    }
	  }
	}
      }
      j	=jn+jz[0];
      r1	=rTFC[j];
      j1	=jn+jz[1];
      EZr2	=rTFC[j1];
      if(EZr2 < r1){
	j	=jz[0];
	jz[0]	=jz[1];
	jz[1]	=j;
      }
      if(jj != 2){
	j	=jn+jz[0];
	r1	=rTFC[j];
	j	=jn+jz[1];
	EZr2	=rTFC[j];
	for(j1=2; j1 < jj; j1++){
	  j	=jn+jz[j1];
	  if(r1 > rTFC[j]){
	    r1	=rTFC[j];
	    jz[0]=jz[j1];
	  }
	  j	=jn+jz[j1];
	  if(EZr2 < rTFC[j]){
	    EZr2	=rTFC[j];
	    jz[1]=jz[j1];
	  }
	}
      }
      j	=jn+jz[0];
      j1=j+1;
      s	=(z-zTFC[j])/(zTFC[j1]-zTFC[j]);
      x1TFC[n]	=xTFC[j]+(xTFC[j1]-xTFC[j])*s;
      y1TFC[n]	=yTFC[j]+(yTFC[j1]-yTFC[j])*s;
      j	=jn+jz[1];
      j1=j+1;
      s	=(z-zTFC[j])/(zTFC[j1]-zTFC[j]);
      x2TFC[n]	=xTFC[j]+(xTFC[j1]-xTFC[j])*s;
      y2TFC[n]	=yTFC[j]+(yTFC[j1]-yTFC[j])*s;
    }
    r1TFC[n]	=sqrt(x1TFC[n]*x1TFC[n]+y1TFC[n]*y1TFC[n]);
    r2TFC[n]	=sqrt(x2TFC[n]*x2TFC[n]+y2TFC[n]*y2TFC[n]);
    gf1TFC[n]	=asin(y1TFC[n]/r1TFC[n]);
    if(y1TFC[n] >= 0.){
      if(x1TFC[n] < 0.){
	gf1TFC[n] =EZcgp-gf1TFC[n];
      }
    }
    else{
      if(x1TFC[n] >= 0.){
	gf1TFC[n] +=EZc2gp;
      }
      else{
	gf1TFC[n] =EZcgp-gf1TFC[n];
      }
    }
    if(gf1TFC[n] > EZc2gp)
      gf1TFC[n] -=EZc2gp;
    gf2TFC[n]	=asin(y2TFC[n]/r2TFC[n]);
    if(y2TFC[n] >= 0.){
      if(x2TFC[n] < 0.){
	gf2TFC[n] =EZcgp-gf2TFC[n];
      }
    }
    else{
      if(x2TFC[n] >= 0.){
	gf2TFC[n] +=EZc2gp;
      }
      else{
	gf2TFC[n] =EZcgp-gf2TFC[n];
      }
    }
    if(gf2TFC[n] > EZc2gp)
      gf2TFC[n] -=EZc2gp;
    zEnv[n]	=z;
  }
  jn	=TFCNt;
  for(j=1; j < TFCLt; j++){
    s	=EZc2gp*j/TFCLt;
    csL	=cos(s);
    snL	=sin(s);
    for(n=0; n < TFCNt; n++){
      x1TFC[jn]	=x1TFC[n]*csL-y1TFC[n]*snL;
      y1TFC[jn]	=x1TFC[n]*snL+y1TFC[n]*csL;
      r1TFC[jn]	=r1TFC[n];
      gf1TFC[jn]=gf1TFC[n]+s;
      if(gf1TFC[jn] > EZc2gp)
	gf1TFC[jn] -=EZc2gp;
      x2TFC[jn]	=x2TFC[n]*csL-y2TFC[n]*snL;
      y2TFC[jn]	=x2TFC[n]*snL+y2TFC[n]*csL;
      r2TFC[jn]	=r2TFC[n];
      gf2TFC[jn]=gf2TFC[n]+s;
      if(gf2TFC[jn] > EZc2gp)
	gf2TFC[jn] -=EZc2gp;
      zEnv[jn]	=z;
      jn++;
    }
  }
  n	=0;
  x1TFC[jn]	=x1TFC[n];
  y1TFC[jn]	=y1TFC[n];
  r1TFC[jn]	=r1TFC[n];
  gf1TFC[jn]	=gf1TFC[n]+EZc2gp;
  if(gf1TFC[jn] > EZc2gp)
    gf1TFC[jn] -=EZc2gp;
  x2TFC[jn]	=x2TFC[n];
  y2TFC[jn]	=y2TFC[n];
  r2TFC[jn]	=r2TFC[n];
  gf2TFC[jn]	=gf2TFC[n]+EZc2gp;
  if(gf2TFC[jn] > EZc2gp)
    gf2TFC[jn] -=EZc2gp;
  zEnv[jn]	=z;
  return(0);
}

int ESInitSplGrid1()
{
  int j1,j2,j,jj,n,jn,jz[8],iz;
  double *pz,s,csL,snL;
  double r1,EZr2,gf,z;
  
  for(iz=0; iz < nzSpl1; iz++){
    z	=Zmn+dzSpl*iz;
    for(n=0; n < TFCNt; n++){
      jn	=TFCNp1*n;
      jj	=0;
      j		=jn+jZmnTFC[n];
      if(z <= zTFC[j]){
	jj=1;
      }
      j	=jn+jZmxTFC[n];
      if(zTFC[j] <= z){
	jj=1;
      }
      if(jj){
	x1TFC[n]	=xTFC[j];
	x2TFC[n]	=xTFC[j];
	y1TFC[n]	=yTFC[j];
	y2TFC[n]	=yTFC[j];
      }    
      else{
	pz	=zTFC+jn;
	jj	=0;
	for(j=0; j < TFCNp; j++){
	  j1=j+1;
	  if(pz[j] != pz[j1]){
	    if((z-pz[j] <= 0. && z-pz[j1] > 0.) || 
	       (z-pz[j] >= 0. && z-pz[j1] < 0.)){
	      jz[jj]	=j;
	      jj++;
	      if(jj > 7){
		printf("Too many intersections in TFC %d%c\n",n,7);
		jj=7;
	      }
	    }
	  }
	}
	j	=jn+jz[0];
	r1	=rTFC[j];
	j1	=jn+jz[1];
	EZr2	=rTFC[j1];
	if(EZr2 < r1){
	  j	=jz[0];
	  jz[0]	=jz[1];
	  jz[1]	=j;
	}
	if(jj != 2){
	  j	=jn+jz[0];
	  r1	=rTFC[j];
	  j	=jn+jz[1];
	  EZr2	=rTFC[j];
	  for(j1=2; j1 < jj; j1++){
	    j	=jn+jz[j1];
	    if(r1 > rTFC[j]){
	      r1	=rTFC[j];
	      jz[0]=jz[j1];
	    }
	    j	=jn+jz[j1];
	    if(EZr2 < rTFC[j]){
	      EZr2	=rTFC[j];
	      jz[1]=jz[j1];
	    }
	  }
	}
	j	=jn+jz[0];
	j1=j+1;
	s	=(z-zTFC[j])/(zTFC[j1]-zTFC[j]);
	x1TFC[n]	=xTFC[j]+(xTFC[j1]-xTFC[j])*s;
	y1TFC[n]	=yTFC[j]+(yTFC[j1]-yTFC[j])*s;
	j	=jn+jz[1];
	j1=j+1;
	s	=(z-zTFC[j])/(zTFC[j1]-zTFC[j]);
	x2TFC[n]	=xTFC[j]+(xTFC[j1]-xTFC[j])*s;
	y2TFC[n]	=yTFC[j]+(yTFC[j1]-yTFC[j])*s;
      }
      r1TFC[n]	=sqrt(x1TFC[n]*x1TFC[n]+y1TFC[n]*y1TFC[n]);
      r2TFC[n]	=sqrt(x2TFC[n]*x2TFC[n]+y2TFC[n]*y2TFC[n]);
      gf1TFC[n]	=asin(y1TFC[n]/r1TFC[n]);
      if(y1TFC[n] >= 0.){
	if(x1TFC[n] < 0.){
	  gf1TFC[n] =EZcgp-gf1TFC[n];
	}
      }
      else{
	if(x1TFC[n] >= 0.){
	  gf1TFC[n] +=EZc2gp;
	}
	else{
	  gf1TFC[n] =EZcgp-gf1TFC[n];
	}
      }
      if(gf1TFC[n] > EZc2gp)
	gf1TFC[n] -=EZc2gp;
      gf2TFC[n]	=asin(y2TFC[n]/r2TFC[n]);
      if(y2TFC[n] >= 0.){
	if(x2TFC[n] < 0.){
	  gf2TFC[n] =EZcgp-gf2TFC[n];
	}
      }
      else{
	if(x2TFC[n] >= 0.){
	  gf2TFC[n] +=EZc2gp;
	}
	else{
	  gf2TFC[n] =EZcgp-gf2TFC[n];
	}
      }
      if(gf2TFC[n] > EZc2gp)
	gf2TFC[n] -=EZc2gp;
    }
    jn	=TFCNt;
    for(j=1; j < TFCLt; j++){
      s	=EZc2gp*j/TFCLt;
      csL	=cos(s);
      snL	=sin(s);
      for(n=0; n < TFCNt; n++){
	x1TFC[jn]	=x1TFC[n]*csL-y1TFC[n]*snL;
	y1TFC[jn]	=x1TFC[n]*snL+y1TFC[n]*csL;
	r1TFC[jn]	=r1TFC[n];
	gf1TFC[jn]=gf1TFC[n]+s;
	if(gf1TFC[jn] > EZc2gp)
	  gf1TFC[jn] -=EZc2gp;
	x2TFC[jn]	=x2TFC[n]*csL-y2TFC[n]*snL;
	y2TFC[jn]	=x2TFC[n]*snL+y2TFC[n]*csL;
	r2TFC[jn]	=r2TFC[n];
	gf2TFC[jn]=gf2TFC[n]+s;
	if(gf2TFC[jn] > EZc2gp)
	  gf2TFC[jn] -=EZc2gp;
	jn++;
      }
    }
    n	=0;
    x1TFC[jn]	=x1TFC[n];
    y1TFC[jn]	=y1TFC[n];
    r1TFC[jn]	=r1TFC[n];
    gf1TFC[jn]	=gf1TFC[n]+EZc2gp;
    if(gf1TFC[jn] > EZc2gp)
      gf1TFC[jn] -=EZc2gp;
    x2TFC[jn]	=x2TFC[n];
    y2TFC[jn]	=y2TFC[n];
    r2TFC[jn]	=r2TFC[n];
    gf2TFC[jn]	=gf2TFC[n]+EZc2gp;
    if(gf2TFC[jn] > EZc2gp)
      gf2TFC[jn] -=EZc2gp;
    
    for(n=0; n < ESNt1; n++){
      r1Env[n]	=Rmn;	
      r2Env[n]	=Rmx;	
      Nr1Spl[n][iz]=-1;
      Nr2Spl[n][iz]=-1;
      zEnv[n]	=z;
    }
    jn	=0;
    j1	=0;
    j2	=0;
    for(j=0; j < TFCLt; j++){
      for(n=0; n < TFCNt; n++){
	gf	=gf1TFC[jn];
	while(j1 < ESNt && gf > gf1T[j1+1]){
	  j1++;
	}
	while(j1 > 0 && gf < gf1T[j1]){
	  j1--;
	}
	if(r1Env[j1] <= r1TFC[jn]){
	  jj	=(r1TFC[jn]-Rmn)/drSpl+1;
	  r1Env[j1]=Rmn+drSpl*jj;
	  Nr1Spl[j1][iz]=jj;
	}
	gf	=gf2TFC[jn];
	while(j2 < ESNt && gf > gf1T[j2+1]){
	  j2++;
	}
	while(j2 > 0 && gf < gf1T[j2]){
	  j2--;
	}
	if(r2Env[j2] > r2TFC[jn]){
	  jj	=(r2TFC[jn]-Rmn)/drSpl;
	  r2Env[j2]=Rmn+drSpl*jj;
	  Nr2Spl[j2][iz]=jj;
	}
	jn++;
      }
    }
    j1	=0;
    while(Nr1Spl[j1][iz] == -1){
      j1++;
    }
    jj	=j1;
    j2	=j1+1;
    while(Nr1Spl[j2][iz] == -1){
      j2++;
    }
    j	=j1+1;
    while(j != jj){
      if(Nr1Spl[j][iz] == -1){
	if(j2 > j1){
	  s	=(gf1T[j]-gf1T[j1])/(gf1T[j2]-gf1T[j1]);
	}
	else{
	  if(j > j1)
	    s	=(gf1T[j]-gf1T[j1])/(gf1T[j2]+EZc2gp-gf1T[j1]);
	  else
	    s	=(gf1T[j]+EZc2gp-gf1T[j1])/(gf1T[j2]+EZc2gp-gf1T[j1]);
	}
	Nr1Spl[j][iz]=Nr1Spl[j1][iz]+(Nr1Spl[j2][iz]-Nr1Spl[j1][iz])*s;
	r1Env[j]=Rmn+drSpl*Nr1Spl[j][iz];
	j++;
      }
      else{
	j1	=j2;
	j		=j1+1;
	j2	=j;
	if(j2 == ESNt){
	  j2=0;
	}
	while(Nr1Spl[j2][iz] == -1){
	  j2++;
	  if(j2 == ESNt){
	    j2=0;
	  }
	}
      }
      if(j == ESNt){
	j=0;
      }
    }
    Nr1Spl[ESNt][iz]=Nr1Spl[0][iz];
    r1Env[ESNt]=r1Env[0];
    
    j1=0;
    while(Nr2Spl[j1][iz] == -1){
      j1++;
    }
    jj	=j1;
    j2	=j1+1;
    while(Nr2Spl[j2][iz] == -1){
      j2++;
    }
    j	=j1+1;
    while(j != jj){
      if(Nr2Spl[j][iz] == -1){
	if(j2 > j1){
	  s	=(gf1T[j]-gf1T[j1])/(gf1T[j2]-gf1T[j1]);
	}
	else{
	  if(j > j1)
	    s	=(gf1T[j]-gf1T[j1])/(gf1T[j2]+EZc2gp-gf1T[j1]);
	  else
	    s	=(gf1T[j]+EZc2gp-gf1T[j1])/(gf1T[j2]+EZc2gp-gf1T[j1]);
	}
	Nr2Spl[j][iz]=Nr2Spl[j1][iz]+(Nr2Spl[j2][iz]-Nr2Spl[j1][iz])*s;
	r2Env[j]=Rmn+drSpl*Nr2Spl[j][iz];
	j++;
      }
      else{
	j1	=j2;
	j		=j1+1;
	j2	=j;
	if(j2 == ESNt){
	  j2=0;
	}
	while(Nr2Spl[j2][iz] == -1){
	  j2++;
	  if(j2 == ESNt){
	    j2=0;
	  }
	}
      }
      if(j == ESNt){
	j=0;
      }
    }
    Nr2Spl[ESNt][iz]=Nr2Spl[0][iz];
    r2Env[ESNt]=r2Env[0];
    for(n=0; n < ESNt1; n++){
      if(Nr2Spl[n][iz] <= Nr1Spl[n][iz]){
	Nr2Spl[n][iz]	=Nr1Spl[n][iz];
	r2Env[n]	=r1Env[n];	
      }
    }
  }
  {
    int i,i1,i2,*k2;
    jj	=0;
    for(n=0; n < ESNt1; n++){
      iz=0;
      while(Nr2Spl[n][iz] == Nr1Spl[n][iz]){
	iz++;
      }
      Nz1Spl[n]	=iz;
      j1	=iz;
      iz	=nzSpl;
      while(Nr2Spl[n][iz] == Nr1Spl[n][iz]){
	iz--;
      }
      Nz2Spl[n]	=iz;
      j2	=iz+1;
      if(SplNt){
	if(SplNz[n] = 0){
	  for(j=0; j < SplNz[n]; j++){
	    free(k2kSpl[n][j]);
	  }
	}
	free(k2kSpl[n]);
      }
      SplNz[n]	=j2-j1;
      k2kSpl[n]	=(int**)malloc((j2-j1)*sizeof(int*));
      if(k2kSpl[n] == NULL){
	FailureAlarm((char*)k2kSpl[n],
		     "ESInitSplGrid() - no memory for k2kSpl[n]");
	ESexit(0);
      }
      jn	=0;
      for(j=j1; j < j2; j++){
	i1	=Nr1Spl[n][j];
	i2	=Nr2Spl[n][j]+1;
	k2	=(int*)malloc((i2-i1)*sizeof(int));
	k2kSpl[n][jn]	=k2;
	for(i=i1; i < i2; i++){
	  k2[i-i1]	=jj;
	  jj++;
	}
	jn++;
      }
    }
    if(SplN){
      free(BR);
      free(AR);
      free(xSpl);
    }
    xSpl=(double*)malloc(4*jj*sizeof(double));
    if(xSpl == NULL){
      FailureAlarm((char*)xSpl,"ESInitSplGrid() - no memory for xSpl");
      ESexit(0);
    }
    ySpl=xSpl+jj;
    zSpl=ySpl+jj;
    rSpl=zSpl+jj;

    AR	=(double*)malloc(24*jj*sizeof(double));
    if(AR == NULL){
      FailureAlarm((char*)AR,"ESInitSplGrid() - no memory for AR");
      ESexit(0);
    }
    ARr	=AR+jj;
    ARt	=ARr+jj;
    ARz	=ARt+jj;
    ARrt=ARz+jj;
    ARrz=ARrt+jj;
    ARtz=ARrz+jj;
    ARrtz=ARtz+jj;
    AT	=ARrtz+jj;
    ATr	=AT+jj;
    ATt	=ATr+jj;
    ATz	=ATt+jj;
    ATrt=ATz+jj;
    ATrz=ATrt+jj;
    ATtz=ATrz+jj;
    ATrtz=ATtz+jj;
    AZ	=ATrtz+jj;
    AZr	=AZ+jj;
    AZt	=AZr+jj;
    AZz	=AZt+jj;
    AZrt=AZz+jj;
    AZrz=AZrt+jj;
    AZtz=AZrz+jj;
    AZrtz=AZtz+jj;

    BR	=(double*)malloc(24*jj*sizeof(double));
    if(BR == NULL){
      FailureAlarm((char*)BR,"ESInitSplGrid() - no memory for BR");
      ESexit(0);
    }
    BRr	=BR+jj;
    BRt	=BRr+jj;
    BRz	=BRt+jj;
    BRrt=BRz+jj;
    BRrz=BRrt+jj;
    BRtz=BRrz+jj;
    BRrtz=BRtz+jj;

    BT	=BRrtz+jj;
    BTr	=BT+jj;
    BTt	=BTr+jj;
    BTz	=BTt+jj;
    BTrt=BTz+jj;
    BTrz=BTrt+jj;
    BTtz=BTrz+jj;
    BTrtz=BTtz+jj;

    BZ	=BTrtz+jj;
    BZr	=BZ+jj;
    BZt	=BZr+jj;
    BZz	=BZt+jj;
    BZrt=BZz+jj;
    BZrz=BZrt+jj;
    BZtz=BZrz+jj;
    BZrtz=BZtz+jj;
    SplN	=jj;
    n	=SplN*(4+6*8);
    printf("??? Vector=%8d Words=%8d %3d Mb+%4d Kb+%4d b=%9d\n",SplN,n,
	   (n*8)/0xFFFFF,((n*8)%0xFFFFF)/0x4FF,((n*8)%0xFFFFF)%0x4FF,n*8);
    for(n=0; n < ESNt1; n++){
      j1	=Nz1Spl[n];
      j2	=Nz2Spl[n]+1;
      for(j=j1; j < j2; j++){
	i1	=Nr1Spl[n][j];
	i2	=Nr2Spl[n][j]+1;
	k2	=k2kSpl[n][j-j1];
	for(i=i1; i < i2; i++){
	  jj	=k2[i-i1];
	  rSpl[jj]	=Rmn+drSpl*i;
	  zSpl[jj]	=Zmn+dzSpl*j;
	  xSpl[jj]	=rSpl[jj]*cs1T[n];
	  ySpl[jj]	=rSpl[jj]*sn1T[n];
	}
      }
    }
  }
  return(0);
}

int ESInitSplGrid()
{
  int j1,j2,j,jj,n,jn,jz[8],iz;
  double *pz,s,csL,snL;
  double r1,EZr2,gf,z;
  
  for(iz=0; iz < nzSpl1; iz++){
    z	=Zmn+dzSpl*iz;
    for(n=0; n < TFCNt; n++){
      jn	=TFCNp1*n;
      jj	=0;
      j		=jn+jZmnTFC[n];
      if(z <= zTFC[j]){
	jj=1;
      }
      j	=jn+jZmxTFC[n];
      if(zTFC[j] <= z){
	jj=1;
      }
      if(jj){
	x1TFC[n]	=xTFC[j];
	x2TFC[n]	=xTFC[j];
	y1TFC[n]	=yTFC[j];
	y2TFC[n]	=yTFC[j];
      }    
      else{
	pz	=zTFC+jn;
	jj	=0;
	for(j=0; j < TFCNp; j++){
	  j1=j+1;
	  if(pz[j] != pz[j1]){
	    if((z-pz[j] <= 0. && z-pz[j1] > 0.) || 
	       (z-pz[j] >= 0. && z-pz[j1] < 0.)){
	      jz[jj]	=j;
	      jj++;
	      if(jj > 7){
		printf("Too many intersections in TFC %d%c\n",n,7);
		jj=7;
	      }
	    }
	  }
	}
	j	=jn+jz[0];
	r1	=rTFC[j];
	j1	=jn+jz[1];
	EZr2	=rTFC[j1];
	if(EZr2 < r1){
	  j	=jz[0];
	  jz[0]	=jz[1];
	  jz[1]	=j;
	}
	if(jj != 2){
	  j	=jn+jz[0];
	  r1	=rTFC[j];
	  j	=jn+jz[1];
	  EZr2	=rTFC[j];
	  for(j1=2; j1 < jj; j1++){
	    j	=jn+jz[j1];
	    if(r1 > rTFC[j]){
	      r1	=rTFC[j];
	      jz[0]=jz[j1];
	    }
	    j	=jn+jz[j1];
	    if(EZr2 < rTFC[j]){
	      EZr2	=rTFC[j];
	      jz[1]=jz[j1];
	    }
	  }
	}
	j	=jn+jz[0];
	j1=j+1;
	s	=(z-zTFC[j])/(zTFC[j1]-zTFC[j]);
	x1TFC[n]	=xTFC[j]+(xTFC[j1]-xTFC[j])*s;
	y1TFC[n]	=yTFC[j]+(yTFC[j1]-yTFC[j])*s;
	j	=jn+jz[1];
	j1=j+1;
	s	=(z-zTFC[j])/(zTFC[j1]-zTFC[j]);
	x2TFC[n]	=xTFC[j]+(xTFC[j1]-xTFC[j])*s;
	y2TFC[n]	=yTFC[j]+(yTFC[j1]-yTFC[j])*s;
      }
      r1TFC[n]	=sqrt(x1TFC[n]*x1TFC[n]+y1TFC[n]*y1TFC[n]);
      r2TFC[n]	=sqrt(x2TFC[n]*x2TFC[n]+y2TFC[n]*y2TFC[n]);
      gf1TFC[n]	=asin(y1TFC[n]/r1TFC[n]);
      if(y1TFC[n] >= 0.){
	if(x1TFC[n] < 0.){
	  gf1TFC[n] =EZcgp-gf1TFC[n];
	}
      }
      else{
	if(x1TFC[n] >= 0.){
	  gf1TFC[n] +=EZc2gp;
	}
	else{
	  gf1TFC[n] =EZcgp-gf1TFC[n];
	}
      }
      if(gf1TFC[n] > EZc2gp)
	gf1TFC[n] -=EZc2gp;
      gf2TFC[n]	=asin(y2TFC[n]/r2TFC[n]);
      if(y2TFC[n] >= 0.){
	if(x2TFC[n] < 0.){
	  gf2TFC[n] =EZcgp-gf2TFC[n];
	}
      }
      else{
	if(x2TFC[n] >= 0.){
	  gf2TFC[n] +=EZc2gp;
	}
	else{
	  gf2TFC[n] =EZcgp-gf2TFC[n];
	}
      }
      if(gf2TFC[n] > EZc2gp)
	gf2TFC[n] -=EZc2gp;
    }
    jn	=TFCNt;
    for(j=1; j < TFCLt; j++){
      s	=EZc2gp*j/TFCLt;
      csL	=cos(s);
      snL	=sin(s);
      for(n=0; n < TFCNt; n++){
	x1TFC[jn]	=x1TFC[n]*csL-y1TFC[n]*snL;
	y1TFC[jn]	=x1TFC[n]*snL+y1TFC[n]*csL;
	r1TFC[jn]	=r1TFC[n];
	gf1TFC[jn]=gf1TFC[n]+s;
	if(gf1TFC[jn] > EZc2gp)
	  gf1TFC[jn] -=EZc2gp;
	x2TFC[jn]	=x2TFC[n]*csL-y2TFC[n]*snL;
	y2TFC[jn]	=x2TFC[n]*snL+y2TFC[n]*csL;
	r2TFC[jn]	=r2TFC[n];
	gf2TFC[jn]=gf2TFC[n]+s;
	if(gf2TFC[jn] > EZc2gp)
	  gf2TFC[jn] -=EZc2gp;
	jn++;
      }
    }
    n	=0;
    x1TFC[jn]	=x1TFC[n];
    y1TFC[jn]	=y1TFC[n];
    r1TFC[jn]	=r1TFC[n];
    gf1TFC[jn]	=gf1TFC[n]+EZc2gp;
    if(gf1TFC[jn] > EZc2gp){
      gf1TFC[jn] -=EZc2gp;
    }
    x2TFC[jn]	=x2TFC[n];
    y2TFC[jn]	=y2TFC[n];
    r2TFC[jn]	=r2TFC[n];
    gf2TFC[jn]	=gf2TFC[n]+EZc2gp;
    if(gf2TFC[jn] > EZc2gp)
      gf2TFC[jn] -=EZc2gp;
    
    for(n=0; n < ESNt1; n++){
      r1Env[n]	=Rmn;	
      r2Env[n]	=Rmx;	
      Nr1Spl[n][iz]=-1;
      Nr2Spl[n][iz]=-1;
      zEnv[n]	=z;
    }
    jn	=0;
    j1	=0;
    j2	=0;
    for(j=0; j < TFCLt; j++){
      for(n=0; n < TFCNt; n++){
	gf	=gf1TFC[jn];
	while(j1 < ESNt && gf > gf1T[j1+1]-1e-5){
	  j1++;
	}
	while(j1 > 0 && gf < gf1T[j1]+1e-5){
	  j1--;
	}
	if(j1 == 0 && gf < gf1T[j1]+1e-5){
	  j1	=ESNt-1;
	}
	if(r1Env[j1] <= r1TFC[jn]){
	  jj	=(r1TFC[jn]-Rmn)/drSpl+1;
	  r1Env[j1]=Rmn+drSpl*jj;
	  Nr1Spl[j1][iz]=jj;
	}
	gf	=gf2TFC[jn];
	while(j2 < ESNt && gf > gf1T[j2+1]-1e-5){
	  j2++;
	}
	while(j2 > 0 && gf < gf1T[j2]+1e-5){
	  j2--;
	}
	if(j2 == 0 && gf < gf1T[j2]+1e-5){
	  j2	=ESNt-1;
	}
	if(r2Env[j2] > r2TFC[jn]){
	  jj	=(r2TFC[jn]-Rmn)/drSpl;
	  r2Env[j2]=Rmn+drSpl*jj;
	  Nr2Spl[j2][iz]=jj;
	}
	jn++;
      }
    }
    j1	=0;
    while(Nr1Spl[j1][iz] == -1){
      j1++;
    }
    jj	=j1;
    j2	=j1+1;
    while(Nr1Spl[j2][iz] == -1){
      j2++;
    }
    j	=j1+1;
    while(j != jj){
      if(Nr1Spl[j][iz] == -1){
	if(j2 > j1){
	  s	=(gf1T[j]-gf1T[j1])/(gf1T[j2]-gf1T[j1]);
	}
	else{
	  if(j > j1)
	    s	=(gf1T[j]-gf1T[j1])/(gf1T[j2]+EZc2gp-gf1T[j1]);
	  else
	    s	=(gf1T[j]+EZc2gp-gf1T[j1])/(gf1T[j2]+EZc2gp-gf1T[j1]);
	}
	Nr1Spl[j][iz]=0;
	r1Env[j]=Rmn;
	j++;
      }
      else{
	j1	=j2;
	j	=j1+1;
	j2	=j;
	if(j2 == ESNt){
	  j2=0;
	}
	while(Nr1Spl[j2][iz] == -1){
	  j2++;
	  if(j2 == ESNt){
	    j2=0;
	  }
	}
      }
      if(j == ESNt){
	j=0;
      }
    }
    Nr1Spl[ESNt][iz]=Nr1Spl[0][iz];
    r1Env[ESNt]=r1Env[0];
    
    j1=0;
    while(Nr2Spl[j1][iz] == -1){
      j1++;
    }
    jj	=j1;
    j2	=j1+1;
    while(Nr2Spl[j2][iz] == -1){
      j2++;
    }
    j	=j1+1;
    while(j != jj){
      if(Nr2Spl[j][iz] == -1){
	if(j2 > j1){
	  s	=(gf1T[j]-gf1T[j1])/(gf1T[j2]-gf1T[j1]);
	}
	else{
	  if(j > j1)
	    s	=(gf1T[j]-gf1T[j1])/(gf1T[j2]+EZc2gp-gf1T[j1]);
	  else
	    s	=(gf1T[j]+EZc2gp-gf1T[j1])/(gf1T[j2]+EZc2gp-gf1T[j1]);
	}
	Nr2Spl[j][iz]=nrSpl;
	r2Env[j]=Rmn+drSpl*Nr2Spl[j][iz];
	j++;
      }
      else{
	j1	=j2;
	j	=j1+1;
	j2	=j;
	if(j2 == ESNt){
	  j2=0;
	}
	while(Nr2Spl[j2][iz] == -1){
	  j2++;
	  if(j2 == ESNt){
	    j2=0;
	  }
	}
      }
      if(j == ESNt){
	j=0;
      }
    }
    Nr2Spl[ESNt][iz]=Nr2Spl[0][iz];
    r2Env[ESNt]=r2Env[0];
    for(n=0; n < ESNt1; n++){
      if(Nr2Spl[n][iz] <= Nr1Spl[n][iz]){
	Nr2Spl[n][iz]	=Nr1Spl[n][iz];
	r2Env[n]	=r1Env[n];	
      }
    }
  }
  {
    int i,i1,i2,*k2;
    jj	=0;
    for(n=0; n < ESNt1; n++){
      iz	=0;
      while(Nr2Spl[n][iz] == Nr1Spl[n][iz]){
	iz++;
      }
      Nz1Spl[n]	=iz;
      j1	=iz;
      iz	=nzSpl;
      while(Nr2Spl[n][iz] == Nr1Spl[n][iz]){
	iz--;
      }
      Nz2Spl[n]	=iz;
      j2	=iz+1;
      if(SplNt){
	if(SplNz[n] = 0){
	  for(j=0; j < SplNz[n]; j++){
	    free(k2kSpl[n][j]);
	  }
	}
	free(k2kSpl[n]);
      }
      SplNz[n]	=j2-j1;
      k2kSpl[n]	=(int**)malloc((j2-j1)*sizeof(int*));
      if(k2kSpl[n] == NULL){
	FailureAlarm((char*)k2kSpl[n],
		     "ESInitSplGrid() - no memory for k2kSpl[n]");
	ESexit(0);
      }
      jn	=0;
      for(j=j1; j < j2; j++){
	i1	=Nr1Spl[n][j];
	i2	=Nr2Spl[n][j]+1;
	k2	=(int*)malloc((i2-i1)*sizeof(int));
	k2kSpl[n][jn]	=k2;
	for(i=i1; i < i2; i++){
	  k2[i-i1]	=jj;
	  jj++;
	}
	jn++;
      }
    }
    if(SplN){
      free(BR);
      free(AR);
      free(xSpl);
    }
    xSpl=(double*)malloc(4*jj*sizeof(double));
    if(xSpl == NULL){
      FailureAlarm((char*)xSpl,"ESInitSplGrid() - no memory for xSpl");
      ESexit(0);
    }
    ySpl=xSpl+jj;
    zSpl=ySpl+jj;
    rSpl=zSpl+jj;

    AR	=(double*)malloc(24*jj*sizeof(double));
    if(AR == NULL){
      FailureAlarm((char*)AR,"ESInitSplGrid() - no memory for AR");
      ESexit(0);
    }
    ARr	=AR+jj;
    ARt	=ARr+jj;
    ARz	=ARt+jj;
    ARrt=ARz+jj;
    ARrz=ARrt+jj;
    ARtz=ARrz+jj;
    ARrtz=ARtz+jj;
    AT	=ARrtz+jj;
    ATr	=AT+jj;
    ATt	=ATr+jj;
    ATz	=ATt+jj;
    ATrt=ATz+jj;
    ATrz=ATrt+jj;
    ATtz=ATrz+jj;
    ATrtz=ATtz+jj;
    AZ	=ATrtz+jj;
    AZr	=AZ+jj;
    AZt	=AZr+jj;
    AZz	=AZt+jj;
    AZrt=AZz+jj;
    AZrz=AZrt+jj;
    AZtz=AZrz+jj;
    AZrtz=AZtz+jj;

    BR	=(double*)malloc(24*jj*sizeof(double));
    if(BR == NULL){
      FailureAlarm((char*)BR,"ESInitSplGrid() - no memory for BR");
      ESexit(0);
    }
    BRr	=BR+jj;
    BRt	=BRr+jj;
    BRz	=BRt+jj;
    BRrt=BRz+jj;
    BRrz=BRrt+jj;
    BRtz=BRrz+jj;
    BRrtz=BRtz+jj;

    BT	=BRrtz+jj;
    BTr	=BT+jj;
    BTt	=BTr+jj;
    BTz	=BTt+jj;
    BTrt=BTz+jj;
    BTrz=BTrt+jj;
    BTtz=BTrz+jj;
    BTrtz=BTtz+jj;

    BZ	=BTrtz+jj;
    BZr	=BZ+jj;
    BZt	=BZr+jj;
    BZz	=BZt+jj;
    BZrt=BZz+jj;
    BZrz=BZrt+jj;
    BZtz=BZrz+jj;
    BZrtz=BZtz+jj;
    SplN	=jj;
    n	=SplN*(4+6*8);
    printf("??? Vector=%8d Words=%8d %3d Mb+%4d Kb+%4d b=%9d\n",SplN,n,
	   (n*8)/0xFFFFF,((n*8)%0xFFFFF)/0x4FF,((n*8)%0xFFFFF)%0x4FF,n*8);
    for(n=0; n < ESNt1; n++){
      j1	=Nz1Spl[n];
      j2	=Nz2Spl[n]+1;
      for(j=j1; j < j2; j++){
	i1	=Nr1Spl[n][j];
	i2	=Nr2Spl[n][j]+1;
	k2	=k2kSpl[n][j-j1];
	for(i=i1; i < i2; i++){
	  jj	=k2[i-i1];
	  rSpl[jj]	=Rmn+drSpl*i;
	  zSpl[jj]	=Zmn+dzSpl*j;
	  xSpl[jj]	=rSpl[jj]*cs1T[n];
	  ySpl[jj]	=rSpl[jj]*sn1T[n];
	}
      }
    }
  }
  return(0);
}

int ESReInitSplGrid()
{
  if(NzSpl < nzSpl1){
    int n;
    for(n=0; n < ESNt1; n++){
      if(NzSpl){
	free(Nr1Spl[n]);
      }
      Nr1Spl[n]	=(int*)malloc(2*nzSpl1*sizeof(int));
      if(Nr1Spl[n] == NULL){
	FailureAlarm((char*)Nr1Spl[n],
		     "ESReInitSplGrid() - no memory for Nr1Spl[n]");
	ESexit(0);
      }
      Nr2Spl[n]	=Nr1Spl[n]+nzSpl1;
    }
    NzSpl =nzSpl1;
  }
  return(0);
}

int ESDeInitSplGrid()
{
  if(NzSpl){
    int n;
    for(n=0; n < ESNt1; n++){
      free(Nr1Spl[n]);
    }
    NzSpl=0;
  }
  return(0);
}

int ESReInitEnvTFC()
{
  if(SplNt < ESNt1){
    int n;
    if(x1TFC != NULL){
      free(x1TFC);
    }
    n	=TFCNt*TFCLt+1;
    x1TFC=(double*)malloc(8*n*sizeof(double));
    x2TFC=x1TFC	+n;
    y1TFC=x2TFC	+n;
    y2TFC=y1TFC	+n;
    r1TFC=y2TFC	+n;
    r2TFC=r1TFC	+n;
    gf1TFC=r2TFC	+n;
    gf2TFC=gf1TFC	+n;
    if(ReInitArray((void**)&k2kSpl,SplNt,ESNt1,sizeof(int**)) < 0){
      FailureAlarm((char*)k2kSpl,"ESReInitEnvTFC() - no memory for k2kSpl");
      ESexit(0);
    }
    if(ReInitArray((void**)&Nr1Spl,SplNt,ESNt1,sizeof(int*)) < 0){
      FailureAlarm((char*)Nr1Spl,"ESReInitEnvTFC() - no memory for Nr1Spl");
      ESexit(0);
    }
    if(ReInitArray((void**)&Nr2Spl,SplNt,ESNt1,sizeof(int*)) < 0){
      FailureAlarm((char*)Nr2Spl,"ESReInitEnvTFC() - no memory for Nr2Spl");
      ESexit(0);
    }
    if(ReInitArray((void**)&SplNz,SplNt,ESNt1,sizeof(int)) < 0){
      FailureAlarm((char*)SplNz,"ESReInitEnvTFC() - no memory for SplNz");
      ESexit(0);
    }
    if(ReInitArray((void**)&Nz1Spl,SplNt,ESNt1,sizeof(int)) < 0){
      FailureAlarm((char*)Nz1Spl,"ESReInitEnvTFC() - no memory for Nz1Spl");
      ESexit(0);
    }
    if(ReInitArray((void**)&Nz2Spl,SplNt,ESNt1,sizeof(int)) < 0){
      FailureAlarm((char*)Nz2Spl,"ESReInitEnvTFC() - no memory for Nz2Spl");
      ESexit(0);
    }
    if(ReInitArray((void**)&x1Env,SplNt,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)x1Env,"ESReInitEnvTFC() - no memory for x1Env");
      ESexit(0);
    }
    if(ReInitArray((void**)&x2Env,SplNt,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)x2Env,"ESReInitEnvTFC() - no memory for x2Env");
      ESexit(0);
    }
    if(ReInitArray((void**)&y1Env,SplNt,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)y1Env,"ESReInitEnvTFC() - no memory for y1Env");
      ESexit(0);
    }
    if(ReInitArray((void**)&y2Env,SplNt,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)y2Env,"ESReInitEnvTFC() - no memory for y2Env");
      ESexit(0);
    }
    if(ReInitArray((void**)&zEnv,SplNt,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)zEnv,"ESReInitEnvTFC() - no memory for zEnv");
      ESexit(0);
    }
    if(ReInitArray((void**)&r1Env,SplNt,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)r1Env,"ESReInitEnvTFC() - no memory for r1Env");
      ESexit(0);
    }
    if(ReInitArray((void**)&r2Env,SplNt,ESNt1,sizeof(double)) < 0){
      FailureAlarm((char*)r2Env,"ESReInitEnvTFC() - no memory for r2Env");
      ESexit(0);
    }
    for(n=SplNt; n < ESNt1; n++){
      SplNz[n]	=0;
      k2kSpl[n]	=NULL;
    }
    SplNt =ESNt1;
    ESReInitSplGrid();
    ESmemlev |=0x00080000;
  }
  return(0);
}

int ESDeInitEnvTFC()
{
  if(SplNt){
    int n,j;
    for(n=0; n < SplNt; n++){
      if(SplNz[n] = 0){
	for(j=0; j < SplNz[n]; j++){
	  free(k2kSpl[n][j]);
	}
      }
      free(k2kSpl[n]);
    }
    ESDeInitSplGrid();
    free(r2Env);
    free(r1Env);
    free(zEnv);
    free(y2Env);
    free(y1Env);
    free(x2Env);
    free(x1Env);
    free(x1TFC);
    free(Nz2Spl);
    free(Nz1Spl);
    free(SplNz);
    free(Nr2Spl);
    free(Nr1Spl);
    free(k2kSpl);
    SplNt=0;
  }
  if(SplN){
    free(rSpl);
    free(zSpl);
    free(ySpl);
    free(xSpl);
    SplN=0;
  }
  return(0);
}

extern double Ax,Ay,Az,Bx,By,Bz;

int ESInitArtzBrtzSpl()
{
  int N,J,J1,J2,JJ,I,I1,I2;
  double X,Z;
  double aX,aXx,aXz,aXy,aXxy,aXxz,aXyy,aXyz,aXzz,aXxyy,aXxyz,aXxzz,aXyyz,aXyzz
    ,aXxyyz,aXxyzz;
  double aY,aYx,aYy,aYz,aYxx,aYxy,aYxz,aYyz,aYzz,aYxxy,aYxxz,aYxyz,aYxzz,aYyzz
    ,aYxxyz,aYxyzz;
  double aZ,aZx,aZy,aZz,aZxx,aZxy,aZxz,aZyy,aZyz,aZxxy,aZxxz,aZxyy,aZxyz,aZyyz
    ,aZxxyz,aZxyyz;
  double F,Fx,Fy,Fz,Fxx,Fxy,Fxz,Fyy,Fyz,Fzz,Fxxy,Fxxz,Fxyy,Fxyz,Fxzz,Fyyz,Fyzz
    ,Fxxyz,Fxyyz,Fxyzz,rF;

  int j,n,jn,L;
  double csLL,snLL,csT,snT,x0,y0,x1,y1,z1,x2,y2,EZz2;
  double Ix,Iy,Iz,dL,R,rR,Z1,Rx,Ry,Rz;
  double FxFx,FyFy,FzFz,FxFy,FxFz,FyFz;

  time(&sTime2);
  printf("%s - Start of  ESInitArtzBrtzSpl()\n",ctime(&sTime2));
  clock();

  JJ	=0;
  for(N=0; N < ESNt; N++){
    j=ESNt/4;
    if(N%j == 0){
      sTime1=sTime2;
      time(&sTime2);
      printf("%10.3e secs\n",difftime(sTime2,sTime1));
      printf("%10.3e secs of CPU\n",(double)(clock()/CLOCKS_PER_SEC));
    }
#undef Ins_SPL2
    csT	=cs1T[N];
    snT =sn1T[N];
    J1	=Nz1Spl[N];
    J2	=Nz2Spl[N]+1;
    for(J=J1; J < J2; J++){
      Z		=Zmn+dzSpl*J;
      I1	=Nr1Spl[N][J];
      I2	=Nr2Spl[N][J]+1;
      for(I=I1; I < I2; I++){
	X	=Rmn+drSpl*I;
#undef Ins_SPL3
	aX	=0.;
	aXx	=0.;
	aXz	=0.;
	aXy	=0.;
	aXxy	=0.;
	aXxz	=0.;
	aXyy	=0.;
	aXyz	=0.;
	aXzz	=0.;
	aXxyy	=0.;
	aXxyz	=0.;
	aXxzz	=0.;
	aXyyz	=0.;
	aXyzz	=0.;
	aXxyyz	=0.;
	aXxyzz	=0.;
	aY	=0.;
	aYx	=0.;
	aYy	=0.;
	aYz	=0.;
	aYxx	=0.;
	aYxy	=0.;
	aYxz	=0.;
	aYyz	=0.;
	aYzz	=0.;
	aYxxy	=0.;
	aYxxz	=0.;
	aYxyz	=0.;
	aYxzz	=0.;
	aYyzz	=0.;
	aYxxyz	=0.;
	aYxyzz	=0.;
	aZ	=0.;
	aZx	=0.;
	aZy	=0.;
	aZz	=0.;
	aZxx	=0.;
	aZxy	=0.;
	aZxz	=0.;
	aZyy	=0.;
	aZyz	=0.;
	aZxxy	=0.;
	aZxxz	=0.;
	aZxyy	=0.;
	aZxyz	=0.;
	aZyyz	=0.;
	aZxxyz	=0.;
	aZxyyz	=0.;
	for(L=0; L < TFCLt; L++){
	  csLL	=csL[L]*csT+snL[L]*snT;
	  snLL	=snL[L]*csT-csL[L]*snT;
	  jn	=0;
	  for(n=0; n < TFCNt; n++){
	    j	=0;
	    do{
	      x0	=xTFC[jn];
	      y0	=yTFC[jn];
	      x1	=x0*csLL-y0*snLL;
	      y1	=x0*snLL+y0*csLL;	
	      z1	=zTFC[jn];
	      jn++;
	      j++;
	    }while(x1 == X && y1 == 0. && z1 == Z);
	    Rx	=x1-X;
	    Ry	=y1;
	    Rz	=z1-Z;
	    R	=sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
	    rR	=1./R;
	    Rx	*=rR;
	    Ry	*=rR;
	    Rz	*=rR;
	    while(j < TFCNp1){
	      x0	=xTFC[jn];
	      y0	=yTFC[jn];
	      x2	=x0*csLL-y0*snLL;
	      y2	=x0*snLL+y0*csLL;	
	      EZz2	=zTFC[jn];
	      if(x2 == X && y2 == 0. && EZz2 == Z){
		jn++;
		j++;
		continue;
	      }
	      Ix	=x2-x1;
	      Iy	=y2-y1;
	      Iz	=EZz2-z1;
	      dL	=1./sqrt(Ix*Ix+Iy*Iy+Iz*Iz);
	      Ix	*=dL;
	      Iy	*=dL;
	      Iz	*=dL;
	      Z1	=(Ix*Rx+Iy*Ry+Iz*Rz)*R;
	      if(Z1 < 0.){
		Z1	=-Z1;
		Ix	=-Ix;
		Iy	=-Iy;
		Iz	=-Iz;
	      }
	      rF	=1./(Z1+R);
	      Fx	=(Rx+Ix)*rF;
	      Fy	=(Ry+Iy)*rF;
	      Fz	=(Rz+Iz)*rF;

	      FxFx	=Rx*rR;
	      FyFy	=Ry*rR;
	      FzFz	=Rz*rR;

	      Fxx	=Rx*FxFx-rR;
	      Fxy	=Rx*FyFy;
	      Fxz	=Rz*FxFx;
	      Fyy	=Ry*FyFy-rR;
	      Fyz	=Ry*FzFz;
	      Fzz	=Rz*FzFz-rR;

	      F		=rR*rF;
	      dL	=2.*Fxy;
	      Fxxyz	=(dL*Fxz+Fxx*Fyz)*F;
	      Fxyyz	=(dL*Fyz+Fyy*Fxz)*F;
	      Fxyzz	=(2.*Fxz*Fyz+Fxy*Fzz)*F;

	      Fxx	*=rF;
	      Fxy	*=rF;
	      Fxz	*=rF;
	      Fyy	*=rF;
	      Fyz	*=rF;
	      Fzz	*=rF;

	      F		=2.*Fxy;
	      Fxxy	=F*FxFx+Fxx*FyFy;
	      Fxyy	=F*FyFy+Fyy*FxFx;
	      F		=2.*Fxz;
	      Fxxz	=F*FxFx+Fxx*FzFz;
	      Fxzz	=F*FzFz+Fzz*FxFx;
	      Fxyz	=Fxz*FyFy+Fyz*FxFx+Fxy*FzFz;


	      F		=2.*Fyz;
	      Fyyz	=F*FyFy+Fyy*FzFz;
	      Fyzz	=F*FzFz+Fzz*FyFy;

	      F		=2.*Fxyz;
	      Fxxyz	+=F*FxFx+Fxxy*FzFz+Fxxz*FyFy;
	      Fxyyz	+=F*FyFy+Fxyy*FzFz+Fyyz*FxFx;
	      Fxyzz	+=F*FzFz+Fxzz*FyFy+Fyzz*FxFx;
	      
	      FxFx	=Fx*Fx;
	      FxFy	=Fx*Fy;
	      FxFz	=Fx*Fz;
	      FyFy	=Fy*Fy;
	      FyFz	=Fy*Fz;
	      FzFz	=Fz*Fz;

	      F		=6.*FxFy*Fz;
	      Fxxyz	+=F*Fx+Fxxz*Fy+Fxxy*Fz+Fxx*Fyz
		+2.*(2.*(Fxz*FxFy+Fxy*FxFz)+Fyz*FxFx+Fxx*FyFz+Fxyz*Fx+Fxy*Fxz);

	      Fxyyz	+=F*Fy+Fyyz*Fx+Fxyy*Fz+Fyy*Fxz
		+2.*(2.*(Fyz*FxFy+Fxy*FyFz)+Fxz*FyFy+Fyy*FxFz+Fxyz*Fy+Fxy*Fyz);

	      Fxyzz	+=F*Fz+Fxzz*Fy+Fyzz*Fx+Fxy*Fzz
		+2.*(2.*(Fxz*FyFz+Fyz*FxFz)+Fzz*FxFy+Fxy*FzFz+Fxyz*Fz+Fxz*Fyz);

	      Fxx	+=FxFx;
	      Fyy	+=FyFy;
	      Fzz	+=FzFz;

	      Fxxy	+=Fxx*Fy+Fxy*Fx;
  	      Fxxz	+=Fxx*Fz+Fxz*Fx;
	      Fxzz	+=Fzz*Fx+Fxz*Fz;
	      Fxyy	+=Fyy*Fx+Fxy*Fy;
	      Fyyz	+=Fyy*Fz+Fyz*Fy;
	      Fyzz	+=Fzz*Fy+Fyz*Fz;
	      Fxyz	+=Fyz*Fx;

	      Fxy	+=FxFy;
	      Fxz	+=FxFz;
	      Fyz	+=FyFz;

	      Fxxy	+=Fxy*Fx;
  	      Fxxz	+=Fxz*Fx;
	      Fxzz	+=Fxz*Fz;
	      Fxyy	+=Fxy*Fy;
	      Fyyz	+=Fyz*Fy;
	      Fyzz	+=Fyz*Fz;
	      Fxyz	+=Fxy*Fz+Fxz*Fy;

	      aXx   	+=Ix*Fx   ;
	      aXy   	+=Ix*Fy   ;
	      aXz   	+=Ix*Fz   ;
	      aXxy  	+=Ix*Fxy  ;
	      aXxz  	+=Ix*Fxz  ;
	      aXyy  	+=Ix*Fyy  ;
	      aXyz  	+=Ix*Fyz  ;
	      aXzz  	+=Ix*Fzz  ;
	      aXxyy 	+=Ix*Fxyy ;
	      aXxyz 	+=Ix*Fxyz ;
	      aXxzz 	+=Ix*Fxzz ;
	      aXyyz 	+=Ix*Fyyz ;
	      aXyzz 	+=Ix*Fyzz ;
	      aXxyyz	+=Ix*Fxyyz;
	      aXxyzz	+=Ix*Fxyzz;

	      aYx   	+=Iy*Fx   ;
	      aYy   	+=Iy*Fy   ;
	      aYz   	+=Iy*Fz   ;
	      aYxx  	+=Iy*Fxx  ;
	      aYxy  	+=Iy*Fxy  ;
	      aYxz  	+=Iy*Fxz  ;
	      aYyz  	+=Iy*Fyz  ;
	      aYzz  	+=Iy*Fzz  ;
	      aYxxy 	+=Iy*Fxxy ;
	      aYxxz 	+=Iy*Fxxz ;
	      aYxyz 	+=Iy*Fxyz ;
	      aYxzz 	+=Iy*Fxzz ;
	      aYyzz 	+=Iy*Fyzz ;
	      aYxxyz	+=Iy*Fxxyz;
	      aYxyzz	+=Iy*Fxyzz;

	      aZx   	+=Iz*Fx   ;
	      aZy   	+=Iz*Fy   ;
	      aZz   	+=Iz*Fz   ;
	      aZxx  	+=Iz*Fxx  ;
	      aZxy  	+=Iz*Fxy  ;
	      aZxz  	+=Iz*Fxz  ;
	      aZyy  	+=Iz*Fyy  ;
	      aZyz  	+=Iz*Fyz  ;
	      aZxxy 	+=Iz*Fxxy ;
	      aZxxz 	+=Iz*Fxxz ;
	      aZxyy 	+=Iz*Fxyy ;
	      aZxyz 	+=Iz*Fxyz ;
	      aZyyz 	+=Iz*Fyyz ;
	      aZxxyz	+=Iz*Fxxyz;
	      aZxyyz	+=Iz*Fxyyz;

	      Rx	=x2-X;
	      Ry	=y2;
	      Rz	=EZz2-Z;
	      Z1	=Ix*Rx+Iy*Ry+Iz*Rz;
	      R		=sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
	      rR	=1./R;
	      Rx	*=rR;
	      Ry	*=rR;
	      Rz	*=rR;
	      F		=log((Z1+R)*rF);
	      aX	+=Ix*F;
	      aY	+=Iy*F;
	      aZ	+=Iz*F;
	      rF	=1./(Z1+R);

	      Fx	=(Rx+Ix)*rF;
	      Fy	=(Ry+Iy)*rF;
	      Fz	=(Rz+Iz)*rF;

	      FxFx	=Rx*rR;
	      FyFy	=Ry*rR;
	      FzFz	=Rz*rR;

	      Fxx	=Rx*FxFx-rR;
	      Fxy	=Rx*FyFy;
	      Fxz	=Rz*FxFx;
	      Fyy	=Ry*FyFy-rR;
	      Fyz	=Ry*FzFz;
	      Fzz	=Rz*FzFz-rR;

	      F		=rR*rF;
	      dL	=2.*Fxy;
	      Fxxyz	=(dL*Fxz+Fxx*Fyz)*F;
	      Fxyyz	=(dL*Fyz+Fyy*Fxz)*F;
	      Fxyzz	=(2.*Fxz*Fyz+Fxy*Fzz)*F;

	      Fxx	*=rF;
	      Fxy	*=rF;
	      Fxz	*=rF;
	      Fyy	*=rF;
	      Fyz	*=rF;
	      Fzz	*=rF;

	      F		=2.*Fxy;
	      Fxxy	=F*FxFx+Fxx*FyFy;
	      Fxyy	=F*FyFy+Fyy*FxFx;
	      F		=2.*Fxz;
	      Fxxz	=F*FxFx+Fxx*FzFz;
	      Fxzz	=F*FzFz+Fzz*FxFx;
	      Fxyz	=Fxz*FyFy+Fyz*FxFx+Fxy*FzFz;
	      F		=2.*Fyz;
	      Fyyz	=F*FyFy+Fyy*FzFz;
	      Fyzz	=F*FzFz+Fzz*FyFy;

	      F		=2.*Fxyz;
	      Fxxyz	+=F*FxFx+Fxxy*FzFz+Fxxz*FyFy;
	      Fxyyz	+=F*FyFy+Fxyy*FzFz+Fyyz*FxFx;
	      Fxyzz	+=F*FzFz+Fxzz*FyFy+Fyzz*FxFx;
	      
	      FxFx	=Fx*Fx;
	      FxFy	=Fx*Fy;
	      FxFz	=Fx*Fz;
	      FyFy	=Fy*Fy;
	      FyFz	=Fy*Fz;
	      FzFz	=Fz*Fz;

	      F		=6.*FxFy*Fz;
	      Fxxyz	+=F*Fx+Fxxz*Fy+Fxxy*Fz+Fxx*Fyz
		+2.*(2.*(Fxz*FxFy+Fxy*FxFz)+Fyz*FxFx+Fxx*FyFz+Fxyz*Fx+Fxy*Fxz);

	      Fxyyz	+=F*Fy+Fyyz*Fx+Fxyy*Fz+Fyy*Fxz
		+2.*(2.*(Fyz*FxFy+Fxy*FyFz)+Fxz*FyFy+Fyy*FxFz+Fxyz*Fy+Fxy*Fyz);

	      Fxyzz	+=F*Fz+Fxzz*Fy+Fyzz*Fx+Fxy*Fzz
		+2.*(2.*(Fxz*FyFz+Fyz*FxFz)+Fzz*FxFy+Fxy*FzFz+Fxyz*Fz+Fxz*Fyz);

	      Fxx	+=FxFx;
	      Fyy	+=FyFy;
	      Fzz	+=FzFz;

	      Fxxy	+=Fxx*Fy+Fxy*Fx;
  	      Fxxz	+=Fxx*Fz+Fxz*Fx;
	      Fxzz	+=Fzz*Fx+Fxz*Fz;
	      Fxyy	+=Fyy*Fx+Fxy*Fy;
	      Fyyz	+=Fyy*Fz+Fyz*Fy;
	      Fyzz	+=Fzz*Fy+Fyz*Fz;
	      Fxyz	+=Fyz*Fx;

	      Fxy	+=FxFy;
	      Fxz	+=FxFz;
	      Fyz	+=FyFz;

	      Fxxy	+=Fxy*Fx;
  	      Fxxz	+=Fxz*Fx;
	      Fxzz	+=Fxz*Fz;
	      Fxyy	+=Fxy*Fy;
	      Fyyz	+=Fyz*Fy;
	      Fyzz	+=Fyz*Fz;
	      Fxyz	+=Fxy*Fz+Fxz*Fy;

	      aXx   	-=Ix*Fx   ;
	      aXy   	-=Ix*Fy   ;
	      aXz   	-=Ix*Fz   ;
	      aXxy  	-=Ix*Fxy  ;
	      aXxz  	-=Ix*Fxz  ;
	      aXyy  	-=Ix*Fyy  ;
	      aXyz  	-=Ix*Fyz  ;
	      aXzz  	-=Ix*Fzz  ;
	      aXxyy 	-=Ix*Fxyy ;
	      aXxyz 	-=Ix*Fxyz ;
	      aXxzz 	-=Ix*Fxzz ;
	      aXyyz 	-=Ix*Fyyz ;
	      aXyzz 	-=Ix*Fyzz ;
	      aXxyyz	-=Ix*Fxyyz;
	      aXxyzz	-=Ix*Fxyzz;

	      aYx   	-=Iy*Fx   ;
	      aYy   	-=Iy*Fy   ;
	      aYz   	-=Iy*Fz   ;
	      aYxx  	-=Iy*Fxx  ;
	      aYxy  	-=Iy*Fxy  ;
	      aYxz  	-=Iy*Fxz  ;
	      aYyz  	-=Iy*Fyz  ;
	      aYzz  	-=Iy*Fzz  ;
	      aYxxy 	-=Iy*Fxxy ;
	      aYxxz 	-=Iy*Fxxz ;
	      aYxyz 	-=Iy*Fxyz ;

	      aYxzz 	-=Iy*Fxzz ;
	      aYyzz 	-=Iy*Fyzz ;
	      aYxxyz	-=Iy*Fxxyz;
	      aYxyzz	-=Iy*Fxyzz;

	      aZx   	-=Iz*Fx   ;
	      aZy   	-=Iz*Fy   ;
	      aZz   	-=Iz*Fz   ;
	      aZxx  	-=Iz*Fxx  ;
	      aZxy  	-=Iz*Fxy  ;
	      aZxz  	-=Iz*Fxz  ;
	      aZyy  	-=Iz*Fyy  ;
	      aZyz  	-=Iz*Fyz  ;
	      aZxxy 	-=Iz*Fxxy ;
	      aZxxz 	-=Iz*Fxxz ;
	      aZxyy 	-=Iz*Fxyy ;
	      aZxyz 	-=Iz*Fxyz ;
	      aZyyz 	-=Iz*Fyyz ;
	      aZxxyz	-=Iz*Fxxyz;
	      aZxyyz	-=Iz*Fxyyz;
	      z1	=EZz2;
	      x1	=x2;
	      y1	=y2;
	      j++;
	      jn++;
	    }
	  }
	}
	AR[JJ]		=0.1*aX;
	ARr[JJ]		=0.1*aXx;
	ARt[JJ]		=0.1*(X*aXy+aY);
	ARrt[JJ]	=0.1*(X*aXxy+aXy+aYx);
	ARz[JJ]		=0.1*aXz;
	ARrz[JJ]	=0.1*aXxz;
	ARtz[JJ]	=0.1*(X*aXyz+aYz);
	ARrtz[JJ]	=0.1*(X*aXxyz+aXyz+aYxz);

	AT[JJ]		=0.1*X*aY;
	ATr[JJ]		=0.1*(X*aYx+aY);
	ATt[JJ]		=0.1*X*(X*aYy-aX);
	ATrt[JJ]	=0.1*(X*(X*aYxy+aYy+aYy-aXx)-aX);
	ATz[JJ]		=0.1*X*aYz;
	ATrz[JJ]	=0.1*(X*aYxz+aYz);
	ATtz[JJ]	=0.1*X*(X*aYyz-aXz);
	ATrtz[JJ]	=0.1*(X*(X*aYxyz+aYyz+aYyz-aXxz)-aXz);

	AZ[JJ]		=0.1*aZ;
	AZr[JJ]		=0.1*aZx;
	AZt[JJ]		=0.1*X*aZy;
	AZrt[JJ]	=0.1*(X*aZxy+aZy);
	AZz[JJ]		=0.1*aZz;
	AZrz[JJ]	=0.1*aZxz;
	AZtz[JJ]	=0.1*X*aZyz;
	AZrtz[JJ]	=0.1*(X*aZxyz+aZyz);

	BR[JJ]		=0.1*(aZy-aYz);
	BRr[JJ]		=0.1*(aZxy-aYxz);
	BRt[JJ]		=0.1*(X*(aZyy-aYyz)+aXz-aZx);
	BRrt[JJ]	=0.1*(X*(aZxyy-aYxyz)+aXxz-aZxx+aZyy-aYyz);
	BRz[JJ]		=0.1*(aZyz-aYzz);
	BRrz[JJ]	=0.1*(aZxyz-aYxzz);
	BRtz[JJ]	=0.1*(X*(aZyyz-aYyzz)+aXzz-aZxz);
	BRrtz[JJ]	=0.1*(X*(aZxyyz-aYxyzz)+aXxzz-aZxxz+aZyyz-aYyzz);

	BT[JJ]		=0.1*(aXz-aZx);
	BTr[JJ]		=0.1*(aXxz-aZxx);
	BTt[JJ]		=0.1*(X*(aXyz-aZxy)-aZy+aYz);
	BTrt[JJ]	=0.1*(X*(aXxyz-aZxxy)-aZxy+aYxz+aXyz-aZxy);


	BTz[JJ]		=0.1*(aXzz-aZxz);
	BTrz[JJ]	=0.1*(aXxzz-aZxxz);
	BTtz[JJ]	=0.1*(X*(aXyzz-aZxyz)-aZyz+aYzz);
	BTrtz[JJ]	=0.1*(X*(aXxyzz-aZxxyz)-aZxyz+aYxzz+aXyzz-aZxyz);

	BZ[JJ]		=0.1*(aYx-aXy);
	BZr[JJ]		=0.1*(aYxx-aXxy);
	BZt[JJ]		=0.1*X*(aYxy-aXyy);
	BZrt[JJ]	=0.1*(X*(aYxxy-aXxyy)+aYxy-aXyy);

	BZz[JJ]		=0.1*(aYxz-aXyz);
	BZrz[JJ]	=0.1*(aYxxz-aXxyz);
	BZtz[JJ]	=0.1*X*(aYxyz-aXyyz);
	BZrtz[JJ]	=0.1*(X*(aYxxyz-aXxyyz)+aYxyz-aXyyz);
	JJ++;
      }
    }
  }
  J1	=Nz1Spl[0];
  J2	=Nz2Spl[0]+1;
  j	=0;
  for(J=J1; J < J2; J++){
    I1	=Nr1Spl[0][J];
    I2	=Nr2Spl[0][J]+1;
    for(I=I1; I < I2; I++){
      AR[JJ]    =AR[j];
      ARr[JJ]   =ARr[j];
      ARt[JJ]   =ARt[j];
      ARrt[JJ]  =ARrt[j];
      ARz[JJ]   =ARz[j];
      ARrz[JJ]  =ARrz[j];
      ARtz[JJ]  =ARtz[j];
      ARrtz[JJ] =ARrtz[j];
                            
      AT[JJ]    =AT[j];
      ATr[JJ]   =ATr[j];
      ATt[JJ]   =ATt[j];
      ATrt[JJ]  =ATrt[j];
      ATz[JJ]   =ATz[j];
      ATrz[JJ]  =ATrz[j];
      ATtz[JJ]  =ATtz[j];
      ATrtz[JJ] =ATrtz[j];
                            
      AZ[JJ]    =AZ[j];
      AZr[JJ]   =AZr[j];
      AZt[JJ]   =AZt[j];
      AZrt[JJ]  =AZrt[j];
      AZz[JJ]   =AZz[j];
      AZrz[JJ]  =AZrz[j];
      AZtz[JJ]  =AZtz[j];
      AZrtz[JJ] =AZrtz[j];
                            
      BR[JJ]    =BR[j];
      BRr[JJ]   =BRr[j];
      BRt[JJ]   =BRt[j];
      BRrt[JJ]  =BRrt[j];
      BRz[JJ]   =BRz[j];
      BRrz[JJ]  =BRrz[j];
      BRtz[JJ]  =BRtz[j];
      BRrtz[JJ] =BRrtz[j];
                            
      BT[JJ]    =BT[j];
      BTr[JJ]   =BTr[j];
      BTt[JJ]   =BTt[j];
      BTrt[JJ]  =BTrt[j];
                            
                            
      BTz[JJ]   =BTz[j];
      BTrz[JJ]  =BTrz[j];
      BTtz[JJ]  =BTtz[j];
      BTrtz[JJ] =BTrtz[j];
                            
      BZ[JJ]    =BZ[j];
      BZr[JJ]   =BZr[j];
      BZt[JJ]   =BZt[j];
      BZrt[JJ]  =BZrt[j];
                            
      BZz[JJ]   =BZz[j];
      BZrz[JJ]  =BZrz[j];
      BZtz[JJ]  =BZtz[j];
      BZrtz[JJ] =BZrtz[j];
      
      j++;
      JJ++;
    }
  }
  time(&sTime2);
  printf("%s - End  of  ESInitArtzBrtzSpl() JJ=%d\n",ctime(&sTime2),JJ);

  return(0);
}

int ESSetSplRgFZ(double r,double gf,double z)
{
  int iR,kT,jZ;
  int i,j,kT1;
  double x,X,xx,XX;
  
  iR	=irCell;
  kT	=ktCell;
  jZ	=jzCell;
  while(r1Cell <= r){
    rCell	=r1Cell;
    r1Cell	+=drSpl;
    irCell++;
  }
  while(rCell > r){
    r1Cell	=rCell;
    rCell	-=drSpl;
    irCell--;
  }
  while(t1Cell <= gf){
    tCell	=t1Cell;
    t1Cell	+=dtSpl;
    ktCell++; 
   if(ktCell == ESNt){
      ktCell	=0;
    }
  }
  while(tCell > gf){
    t1Cell	=tCell;
    tCell	-=dtSpl;
    if(ktCell == 0){
      ktCell	=ESNt;
    }
    ktCell--;
  }
  while(z1Cell <= z){
    zCell	=z1Cell;
    z1Cell	+=dzSpl;
    jzCell++;
  }
  while(zCell > z){
    z1Cell	=zCell;
    zCell	-=dzSpl;
    jzCell--;
  }
  if(iR != irCell || kT != ktCell || jZ != jzCell){
    iR	=irCell;
    kT	=ktCell;
    jZ	=jzCell;
    
    kT1=kT+1;
    if(kT1 == ESNt)
      kT1=0;
    if(jZ >= Nz2Spl[kT] || jZ >= Nz2Spl[kT1]){
      printf("??? jZ=%3d >= Nz2[%3d]=%3d Nz2[%3d]=%3d - Out of spline grid%c\n"
	     ,jZ,kT,Nz2Spl[kT],kT1,Nz2Spl[kT1],7);
      printf("??? r=%10.3e gf=%10.3e %10.3e z=%10.3e kT kT1\n"
	     ,r,gf,gf*ESNt/EZc2gp,z);
      return(1);
    }
    j	=jZ-Nz1Spl[kT];
    if(j < 0 || iR >= Nr2Spl[kT][jZ]){
      if(j < 0)
	printf("??? jZ=%3d < Nz1[%3d]=%3d - Out of spline grid%c\n"
	       ,jZ,kT,Nz1Spl[kT],7);
      if(iR >= Nr2Spl[kT][jZ])
	printf("??? iR >= %3d < Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT,jZ,Nr2Spl[kT][jZ],7);
      printf("??? r=%10.3e gf=%10.3e %10.3e z=%10.3e kT\n"
	     ,r,gf,gf*ESNt/EZc2gp,z);
      return(1);
    }
    i	=iR-Nr1Spl[kT][jZ];
    if(i < 0){
      printf("??? iR < %3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	     ,iR,kT,jZ,Nr1Spl[kT][jZ],7);
      printf("??? r=%10.3e gf=%10.3e %10.3e z=%10.3e kT jZ\n"
	     ,r,gf,gf*ESNt/EZc2gp,z);
      return(1);
    }
    J000	=k2kSpl[kT][j][i];
    J100	=J000+1;
    i	=iR-Nr1Spl[kT][jZ+1];
    if(i < 0 || iR >= Nr2Spl[kT][jZ+1]){
      if(i < 0)
	printf("??? iR=%3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT,jZ+1,Nr1Spl[kT][jZ+1],7);
      if(iR >= Nr2Spl[kT][jZ+1])
	printf("??? iR=%3d >= Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT,jZ+1,Nr2Spl[kT][jZ+1],7);
      printf("??? r=%10.3e gf=%10.3e %10.3e z=%10.3e kT jZ+1\n"
	     ,r,gf,gf*ESNt/EZc2gp,z);
      return(1);
    }
    J001	=k2kSpl[kT][j+1][i];
    J101	=J001+1;
    j	=jZ-Nz1Spl[kT1];
    if(j < 0 || iR >= Nr2Spl[kT1][jZ]){
      if(j < 0)
	printf("??? jZ=%3d < Nz1[%3d]=%3d - Out of spline grid%c\n"
	       ,jZ,kT1,Nz1Spl[kT1],7);
      if(iR >= Nr2Spl[kT1][jZ])
	printf("??? iR >= %3d < Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ,Nr2Spl[kT1][jZ],7);
      printf("??? r=%10.3e gf=%10.3e %10.3e z=%10.3e kT1 jZ\n"
	     ,r,gf,gf*ESNt/EZc2gp,z);
      return(1);
    }
    i	=iR-Nr1Spl[kT1][jZ];
    if(i < 0){
      if(i < 0)
	printf("??? iR=%3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ,Nr1Spl[kT1][jZ],7);
      printf("??? r=%10.3e gf=%10.3e %10.3e z=%10.3e kT1 jZ\n"
	     ,r,gf,gf*ESNt/EZc2gp,z);
      return(1);
    }
    J010	=k2kSpl[kT1][j][i];
    J110	=J010+1;
    i	=iR-Nr1Spl[kT1][jZ+1];
    if(i < 0 || iR >= Nr2Spl[kT1][jZ+1]){
      if(i < 0)
	printf("??? iR=%3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ+1,Nr1Spl[kT1][jZ+1],7);
      if(iR >= Nr2Spl[kT1][jZ+1])
	printf("??? iR >= %3d < Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ+1,Nr2Spl[kT1][jZ+1],7);
      printf("??? r=%10.3e gf=%10.3e %10.3e z=%10.3e kT1 jZ+1\n"
	     ,r,gf,gf*ESNt/EZc2gp,z);
      return(1);
    }
    J011	=k2kSpl[kT1][j+1][i];
    J111	=J011+1;
  }
  x	=(r-rCell)*rhrCell;
  xx	=x*x;
  X	=1.-x;
  XX	=X*X;
  Pr0=XX*(3.-2.*X);
  Pr1=xx*(3.-2.*x);
  Pr2=XX*x*drSpl;
  Pr3=-xx*X*drSpl;
  
  x	=(gf-tCell)*rhtCell;
  xx	=x*x;
  X	=1.-x;
  XX	=X*X;
  Pt0=XX*(3.-2.*X);
  Pt1=xx*(3.-2.*x);
  Pt2=XX*x*dtSpl;
  Pt3=-xx*X*dtSpl;
  
  x	=(z-zCell)*rhzCell;
  xx	=x*x;
  X	=1.-x;
  XX	=X*X;
  Pz0=XX*(3.-2.*X);
  Pz1=xx*(3.-2.*x);
  Pz2=XX*x*dzSpl;
  Pz3=-xx*X*dzSpl;
  return(0);
}

int ESSetSplDerRTZ(double r,double gf,double z)
{
  int iR,kT,jZ;
  int i,j,kT1;
  double x,X,xx,XX,Xx;
  
  iR	=irCell;
  kT	=ktCell;
  jZ	=jzCell;
  while(r1Cell <= r){
    rCell	=r1Cell;
    r1Cell	+=drSpl;
    irCell++;
  }
  while(rCell > r){
    r1Cell	=rCell;
    rCell	-=drSpl;
    irCell--;
  }
  while(t1Cell <= gf){
    tCell	=t1Cell;
    t1Cell	+=dtSpl;
    if(ktCell == ESNt){
      ktCell	=0;
    }
    ktCell++;
  }
  while(tCell > gf){
    t1Cell	=tCell;
    tCell	-=dtSpl;
    if(ktCell == 0){
      ktCell	=ESNt;
    }
    ktCell--;
  }
  while(z1Cell <= z){
    zCell	=z1Cell;
    z1Cell	+=dzSpl;
    jzCell++;
  }
  while(zCell > z){
    z1Cell	=zCell;
    zCell	-=dzSpl;
    jzCell--;
  }
  if(iR != irCell || kT != ktCell || jZ != jzCell){
    iR	=irCell;
    kT	=ktCell;
    jZ	=jzCell;
    kT1	=kT+1;
    if(kT1 == ESNt)
      kT1=0;
  
    if(jZ >= Nz2Spl[kT] || jZ >= Nz2Spl[kT1]){
      printf("??? jZ=%3d >= Nz2[%3d]=%3d Nz2[%3d]=%3d - Out of spline grid%c\n"
	     ,jZ,kT,Nz2Spl[kT],kT1,Nz2Spl[kT1],7);
      return(1);
    }
    j	=jZ-Nz1Spl[kT];
    if(j < 0 || iR >= Nr2Spl[kT][jZ]){
      if(j < 0)
	printf("??? jZ=%3d < Nz1[%3d]=%3d - Out of spline grid%c\n"
	       ,jZ,kT,Nz1Spl[kT],7);
      if(iR >= Nr2Spl[kT][jZ])
	printf("??? iR >= %3d < Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT,jZ,Nr2Spl[kT][jZ],7);
      return(1);
    }
    i	=iR-Nr1Spl[kT][jZ];
    if(i < 0){
      printf("??? iR < %3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	     ,iR,kT,jZ,Nr1Spl[kT][jZ],7);
      return(1);
    }
    J000	=k2kSpl[kT][j][i];
    J100	=J000+1;
    i	=iR-Nr1Spl[kT][jZ+1];
    if(i < 0 || iR >= Nr2Spl[kT][jZ+1]){
      if(i < 0)
	printf("??? iR=%3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT,jZ+1,Nr1Spl[kT][jZ+1],7);
      if(iR >= Nr2Spl[kT][jZ+1])
	printf("??? iR=%3d >= Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT,jZ+1,Nr2Spl[kT][jZ+1],7);
      return(1);
    }
    J001	=k2kSpl[kT][j+1][i];
    J101	=J001+1;
    j	=jZ-Nz1Spl[kT1];
    if(j < 0 || iR >= Nr2Spl[kT1][jZ]){
      if(j < 0)
	printf("??? jZ=%3d < Nz1[%3d]=%3d - Out of spline grid%c\n"
	       ,jZ,kT1,Nz1Spl[kT1],7);
      if(iR >= Nr2Spl[kT1][jZ])
	printf("??? iR >= %3d < Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ,Nr2Spl[kT1][jZ],7);
      return(1);
    }
    i	=iR-Nr1Spl[kT1][jZ];
    if(i < 0){
      if(i < 0)
	printf("??? iR=%3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ,Nr1Spl[kT1][jZ],7);
      return(1);
    }
    J010	=k2kSpl[kT1][j][i];
    J110	=J010+1;
    i	=iR-Nr1Spl[kT1][jZ+1];
    if(i < 0 || iR >= Nr2Spl[kT1][jZ+1]){
      if(i < 0)
	printf("??? iR=%3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ+1,Nr1Spl[kT1][jZ+1],7);
      if(iR >= Nr2Spl[kT1][jZ+1])
	printf("??? iR >= %3d < Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ+1,Nr2Spl[kT1][jZ+1],7);
      return(1);
    }
    J011	=k2kSpl[kT1][j+1][i];
    J111	=J011+1;
  }
  x	=(r-rCell)*rhrCell;
  xx	=x*x;
  X	=1.-x;
  XX	=X*X;
  Xx	=2.*x*X;
  Pr0	=XX*(3.-2.*X);
  Pr1	=xx*(3.-2.*x);
  Pr2	=XX*x*drSpl;
  Pr3	=-xx*X*drSpl;

  dPr1	=3.*rhrCell*Xx;
  dPr0	=-dPr1;
  dPr2	=(XX-Xx);
  dPr3	=(xx-Xx);
  
  x	=(gf-tCell)*rhtCell;
  xx	=x*x;
  X	=1.-x;
  XX	=X*X;
  Xx	=2.*x*X;
  Pt0	=XX*(3.-2.*X);
  Pt1	=xx*(3.-2.*x);
  Pt2	=XX*x*dtSpl;
  Pt3=	-xx*X*dtSpl;
  dPt1	=3.*rhtCell*Xx;
  dPt0	=-dPt1;
  dPt2	=(XX-Xx);
  dPt3	=(xx-Xx);
  
  x	=(z-zCell)*rhzCell;
  xx	=x*x;
  X	=1.-x;
  XX	=X*X;
  Xx	=2.*x*X;
  Pz0	=XX*(3.-2.*X);
  Pz1	=xx*(3.-2.*x);
  Pz2	=XX*x*dzSpl;
  Pz3	=-xx*X*dzSpl;
  dPz1	=3.*rhzCell*Xx;
  dPz0	=-dPz1;
  dPz2	=(XX-Xx);
  dPz3	=(xx-Xx);
  return(0);
}

extern double ESBgf;
void ESMFLineSplStep1(nq,x,yeq,dy)
     int *nq;
     double *x,*yeq,*dy;
{
  double r,Bt,Art,Arz,Atr,Atz,Azr,Azt;

  ESSetSplDerRTZ(yeq[0],*x,yeq[2]);

  Atr	=(Pt0*(dPr0*(Pz0*AT[J000]
		      +Pz1*AT[J001]
		      +Pz2*ATz[J000]
		      +Pz3*ATz[J001]
		      )
		 +dPr1*(Pz0*AT[J100]
		       +Pz1*AT[J101]
		       +Pz2*ATz[J100]
		       +Pz3*ATz[J101]
		       )
		 +dPr2*(Pz0*ATr[J000]
		       +Pz1*ATr[J001]
		       +Pz2*ATrz[J000]
		       +Pz3*ATrz[J001]
		       )
		 +dPr3*(Pz0*ATr[J100]
		       +Pz1*ATr[J101]
		       +Pz2*ATrz[J100]
		       +Pz3*ATrz[J101]
		       )
		 )
	    +Pt1*(dPr0*(Pz0*AT[J010]
		       +Pz1*AT[J011]
		       +Pz2*ATz[J010]
		       +Pz3*ATz[J011]
		       )
		  +dPr1*(Pz0*AT[J110]
			+Pz1*AT[J111]
			+Pz2*ATz[J110]
			+Pz3*ATz[J111]
			)
		  +dPr2*(Pz0*ATr[J010]
			+Pz1*ATr[J011]
			+Pz2*ATrz[J010]
			+Pz3*ATrz[J011]
			)
		  +dPr3*(Pz0*ATr[J110]
			+Pz1*ATr[J111]
			+Pz2*ATrz[J110]
			+Pz3*ATrz[J111]
			)
		  )
	    +Pt2*(dPr0*(Pz0*ATt[J000]
		       +Pz1*ATt[J001]
		       +Pz2*ATtz[J000]
		       +Pz3*ATtz[J001]
		       )
		  +dPr1*(Pz0*ATt[J100]
			+Pz1*ATt[J101]
			+Pz2*ATtz[J100]
			+Pz3*ATtz[J101]
			)
		  +dPr2*(Pz0*ATrt[J000]
			+Pz1*ATrt[J001]
			+Pz2*ATrtz[J000]
			+Pz3*ATrtz[J001]
			)
		  +dPr3*(Pz0*ATrt[J100]
			+Pz1*ATrt[J101]
			+Pz2*ATrtz[J100]
			+Pz3*ATrtz[J101]
			)
		  )
	    +Pt3*(dPr0*(Pz0*ATt[J010]
		       +Pz1*ATt[J011]
		       +Pz2*ATtz[J010]
		       +Pz3*ATtz[J011]
		       )
		  +dPr1*(Pz0*ATt[J110]
			+Pz1*ATt[J111]
			+Pz2*ATtz[J110]
			+Pz3*ATtz[J111]
			)
		  +dPr2*(Pz0*ATrt[J010]
			+Pz1*ATrt[J011]
			+Pz2*ATrtz[J010]
			+Pz3*ATrtz[J011]
			)
		  +dPr3*(Pz0*ATrt[J110]
			+Pz1*ATrt[J111]
			+Pz2*ATrtz[J110]
			+Pz3*ATrtz[J111]
			)
		  )
	    );
  

  Atz	=(Pt0*(Pr0*(dPz0*AT[J000]
		    +dPz1*AT[J001]
		    +dPz2*ATz[J000]
		    +dPz3*ATz[J001]
		    )
	       +Pr1*(dPz0*AT[J100]
		     +dPz1*AT[J101]
		     +dPz2*ATz[J100]
		     +dPz3*ATz[J101]
		     )
	       +Pr2*(dPz0*ATr[J000]
		     +dPz1*ATr[J001]
		     +dPz2*ATrz[J000]
		     +dPz3*ATrz[J001]
		     )
	       +Pr3*(dPz0*ATr[J100]
		     +dPz1*ATr[J101]
		     +dPz2*ATrz[J100]
		     +dPz3*ATrz[J101]
		     )
	       )
	  +Pt1*(Pr0*(dPz0*AT[J010]
		     +dPz1*AT[J011]
		     +dPz2*ATz[J010]
		     +dPz3*ATz[J011]
		     )
		+Pr1*(dPz0*AT[J110]
		      +dPz1*AT[J111]
		      +dPz2*ATz[J110]
		      +dPz3*ATz[J111]
		      )
		+Pr2*(dPz0*ATr[J010]
		      +dPz1*ATr[J011]
		      +dPz2*ATrz[J010]
		      +dPz3*ATrz[J011]
		      )
		+Pr3*(dPz0*ATr[J110]
		      +dPz1*ATr[J111]
		      +dPz2*ATrz[J110]
		      +dPz3*ATrz[J111]
		      )
		)
	  +Pt2*(Pr0*(dPz0*ATt[J000]
		     +dPz1*ATt[J001]
		     +dPz2*ATtz[J000]
		     +dPz3*ATtz[J001]
		     )
		+Pr1*(dPz0*ATt[J100]
		      +dPz1*ATt[J101]
		      +dPz2*ATtz[J100]
		      +dPz3*ATtz[J101]
		      )
		+Pr2*(dPz0*ATrt[J000]
		      +dPz1*ATrt[J001]
		      +dPz2*ATrtz[J000]
		      +dPz3*ATrtz[J001]
		      )
		+Pr3*(dPz0*ATrt[J100]
		      +dPz1*ATrt[J101]
		      +dPz2*ATrtz[J100]
		      +dPz3*ATrtz[J101]
		      )
		)
	  +Pt3*(Pr0*(dPz0*ATt[J010]
		     +dPz1*ATt[J011]
		     +dPz2*ATtz[J010]
		     +dPz3*ATtz[J011]
		     )
		+Pr1*(dPz0*ATt[J110]
		      +dPz1*ATt[J111]
		      +dPz2*ATtz[J110]
		      +dPz3*ATtz[J111]
		      )
		+Pr2*(dPz0*ATrt[J010]
		      +dPz1*ATrt[J011]
		      +dPz2*ATrtz[J010]
		      +dPz3*ATrtz[J011]
		      )
		+Pr3*(dPz0*ATrt[J110]
		      +dPz1*ATrt[J111]
		      +dPz2*ATrtz[J110]
		      +dPz3*ATrtz[J111]
		      )
		)
	  );
  

  Art	=(dPt0*(Pr0*(Pz0*AR[J000]
		     +Pz1*AR[J001]
		     +Pz2*ARz[J000]
		     +Pz3*ARz[J001]
		     )
		+Pr1*(Pz0*AR[J100]
		      +Pz1*AR[J101]
		      +Pz2*ARz[J100]
		      +Pz3*ARz[J101]
		      )
		+Pr2*(Pz0*ARr[J000]
		      +Pz1*ARr[J001]
		      +Pz2*ARrz[J000]
		      +Pz3*ARrz[J001]
		      )
		+Pr3*(Pz0*ARr[J100]
		      +Pz1*ARr[J101]
		      +Pz2*ARrz[J100]
		      +Pz3*ARrz[J101]
		      )
		)
	  +dPt1*(Pr0*(Pz0*AR[J010]
		      +Pz1*AR[J011]
		      +Pz2*ARz[J010]
		      +Pz3*ARz[J011]
		      )
		 +Pr1*(Pz0*AR[J110]
		       +Pz1*AR[J111]
		       +Pz2*ARz[J110]
		       +Pz3*ARz[J111]
		       )
		 +Pr2*(Pz0*ARr[J010]
		       +Pz1*ARr[J011]
		       +Pz2*ARrz[J010]
		       +Pz3*ARrz[J011]
		       )
		 +Pr3*(Pz0*ARr[J110]
		       +Pz1*ARr[J111]
		       +Pz2*ARrz[J110]
		       +Pz3*ARrz[J111]
		       )
		 )
	  +dPt2*(Pr0*(Pz0*ARt[J000]
		      +Pz1*ARt[J001]
		      +Pz2*ARtz[J000]
		      +Pz3*ARtz[J001]
		      )
		 +Pr1*(Pz0*ARt[J100]
		       +Pz1*ARt[J101]
		       +Pz2*ARtz[J100]
		       +Pz3*ARtz[J101]
		       )
		 +Pr2*(Pz0*ARrt[J000]
		       +Pz1*ARrt[J001]
		       +Pz2*ARrtz[J000]
		       +Pz3*ARrtz[J001]
		       )
		 +Pr3*(Pz0*ARrt[J100]
		       +Pz1*ARrt[J101]
		       +Pz2*ARrtz[J100]
		       +Pz3*ARrtz[J101]
		       )
		 )
	  +dPt3*(Pr0*(Pz0*ARt[J010]
		      +Pz1*ARt[J011]
		      +Pz2*ARtz[J010]
		      +Pz3*ARtz[J011]
		      )
		 +Pr1*(Pz0*ARt[J110]
		       +Pz1*ARt[J111]
		       +Pz2*ARtz[J110]
		       +Pz3*ARtz[J111]
		       )
		 +Pr2*(Pz0*ARrt[J010]
		       +Pz1*ARrt[J011]
		       +Pz2*ARrtz[J010]
		       +Pz3*ARrtz[J011]
		       )
		 +Pr3*(Pz0*ARrt[J110]
		       +Pz1*ARrt[J111]
		       +Pz2*ARrtz[J110]
		       +Pz3*ARrtz[J111]
		       )
		 )
	  );
  

  Arz	=(Pt0*(Pr0*(dPz0*AR[J000]
		    +dPz1*AR[J001]
		    +dPz2*ARz[J000]
		    +dPz3*ARz[J001]
		    )
	       +Pr1*(dPz0*AR[J100]
		     +dPz1*AR[J101]
		     +dPz2*ARz[J100]
		     +dPz3*ARz[J101]
		     )
	       +Pr2*(dPz0*ARr[J000]
		     +dPz1*ARr[J001]
		     +dPz2*ARrz[J000]
		     +dPz3*ARrz[J001]
		     )
	       +Pr3*(dPz0*ARr[J100]
		     +dPz1*ARr[J101]
		     +dPz2*ARrz[J100]
		     +dPz3*ARrz[J101]
		     )
	       )
	  +Pt1*(Pr0*(dPz0*AR[J010]
		     +dPz1*AR[J011]
		     +dPz2*ARz[J010]
		     +dPz3*ARz[J011]
		     )
		+Pr1*(dPz0*AR[J110]
		      +dPz1*AR[J111]
		      +dPz2*ARz[J110]
		      +dPz3*ARz[J111]
		      )
		+Pr2*(dPz0*ARr[J010]
		      +dPz1*ARr[J011]
		      +dPz2*ARrz[J010]
		      +dPz3*ARrz[J011]
		      )
		+Pr3*(dPz0*ARr[J110]
		      +dPz1*ARr[J111]
		      +dPz2*ARrz[J110]
		      +dPz3*ARrz[J111]
		      )
		)
	  +Pt2*(Pr0*(dPz0*ARt[J000]
		     +dPz1*ARt[J001]
		     +dPz2*ARtz[J000]
		     +dPz3*ARtz[J001]
		     )
		+Pr1*(dPz0*ARt[J100]
		      +dPz1*ARt[J101]
		      +dPz2*ARtz[J100]
		      +dPz3*ARtz[J101]
		      )
		+Pr2*(dPz0*ARrt[J000]
		      +dPz1*ARrt[J001]
		      +dPz2*ARrtz[J000]
		      +dPz3*ARrtz[J001]
		      )
		+Pr3*(dPz0*ARrt[J100]
		      +dPz1*ARrt[J101]
		      +dPz2*ARrtz[J100]
		      +dPz3*ARrtz[J101]
		      )
		)
	  +Pt3*(Pr0*(dPz0*ARt[J010]
		     +dPz1*ARt[J011]
		     +dPz2*ARtz[J010]
		     +dPz3*ARtz[J011]
		     )
		+Pr1*(dPz0*ARt[J110]
		      +dPz1*ARt[J111]
		      +dPz2*ARtz[J110]
		      +dPz3*ARtz[J111]
		      )
		+Pr2*(dPz0*ARrt[J010]
		      +dPz1*ARrt[J011]
		      +dPz2*ARrtz[J010]
		      +dPz3*ARrtz[J011]
		      )
		+Pr3*(dPz0*ARrt[J110]
		      +dPz1*ARrt[J111]
		      +dPz2*ARrtz[J110]
		      +dPz3*ARrtz[J111]
		      )
		)
	  );
  
  Azr	=(Pt0*(dPr0*(Pz0*AZ[J000]
		     +Pz1*AZ[J001]
		     +Pz2*AZz[J000]
		     +Pz3*AZz[J001]
		     )
	       +dPr1*(Pz0*AZ[J100]
		      +Pz1*AZ[J101]
		      +Pz2*AZz[J100]
		      +Pz3*AZz[J101]
		      )
	       +dPr2*(Pz0*AZr[J000]
		      +Pz1*AZr[J001]
		      +Pz2*AZrz[J000]
		      +Pz3*AZrz[J001]
		      )
	       +dPr3*(Pz0*AZr[J100]
		      +Pz1*AZr[J101]
		      +Pz2*AZrz[J100]
		      +Pz3*AZrz[J101]
		      )
	       )
	  +Pt1*(dPr0*(Pz0*AZ[J010]
		      +Pz1*AZ[J011]
		      +Pz2*AZz[J010]
		      +Pz3*AZz[J011]
		      )
		+dPr1*(Pz0*AZ[J110]
		       +Pz1*AZ[J111]
		       +Pz2*AZz[J110]
		       +Pz3*AZz[J111]
		       )
		+dPr2*(Pz0*AZr[J010]
		       +Pz1*AZr[J011]
		       +Pz2*AZrz[J010]
		       +Pz3*AZrz[J011]
		       )
		+dPr3*(Pz0*AZr[J110]
		       +Pz1*AZr[J111]
		       +Pz2*AZrz[J110]
		       +Pz3*AZrz[J111]
		       )
		)
	  +Pt2*(dPr0*(Pz0*AZt[J000]
		      +Pz1*AZt[J001]
		      +Pz2*AZtz[J000]
		      +Pz3*AZtz[J001]
		      )
		+dPr1*(Pz0*AZt[J100]
		       +Pz1*AZt[J101]
		       +Pz2*AZtz[J100]
		       +Pz3*AZtz[J101]
		       )
		+dPr2*(Pz0*AZrt[J000]
		       +Pz1*AZrt[J001]
		       +Pz2*AZrtz[J000]
		       +Pz3*AZrtz[J001]
		       )
		+dPr3*(Pz0*AZrt[J100]
		       +Pz1*AZrt[J101]
		       +Pz2*AZrtz[J100]
		       +Pz3*AZrtz[J101]
		       )
		)
	  +Pt3*(dPr0*(Pz0*AZt[J010]
		      +Pz1*AZt[J011]
		      +Pz2*AZtz[J010]
		      +Pz3*AZtz[J011]
		      )
		+dPr1*(Pz0*AZt[J110]
		       +Pz1*AZt[J111]
		       +Pz2*AZtz[J110]
		       +Pz3*AZtz[J111]
		       )
		+dPr2*(Pz0*AZrt[J010]
		       +Pz1*AZrt[J011]
		       +Pz2*AZrtz[J010]
		       +Pz3*AZrtz[J011]
		       )
		+dPr3*(Pz0*AZrt[J110]
		       +Pz1*AZrt[J111]
		       +Pz2*AZrtz[J110]
		       +Pz3*AZrtz[J111]
		       )
		)
	  );

  Azt	=(dPt0*(Pr0*(Pz0*AZ[J000]
		     +Pz1*AZ[J001]
		     +Pz2*AZz[J000]
		     +Pz3*AZz[J001]
		     )
		+Pr1*(Pz0*AZ[J100]
		      +Pz1*AZ[J101]
		      +Pz2*AZz[J100]
		      +Pz3*AZz[J101]
		      )
		+Pr2*(Pz0*AZr[J000]
		      +Pz1*AZr[J001]
		      +Pz2*AZrz[J000]
		      +Pz3*AZrz[J001]
		      )
		+Pr3*(Pz0*AZr[J100]
		      +Pz1*AZr[J101]
		      +Pz2*AZrz[J100]
		      +Pz3*AZrz[J101]
		      )
		)
	  +dPt1*(Pr0*(Pz0*AZ[J010]
		      +Pz1*AZ[J011]
		      +Pz2*AZz[J010]
		      +Pz3*AZz[J011]
		      )
		 +Pr1*(Pz0*AZ[J110]
		       +Pz1*AZ[J111]
		       +Pz2*AZz[J110]
		       +Pz3*AZz[J111]
		       )
		 +Pr2*(Pz0*AZr[J010]
		       +Pz1*AZr[J011]
		       +Pz2*AZrz[J010]
		       +Pz3*AZrz[J011]
		       )
		 +Pr3*(Pz0*AZr[J110]
		       +Pz1*AZr[J111]
		       +Pz2*AZrz[J110]
		       +Pz3*AZrz[J111]
		       )
		 )
	  +dPt2*(Pr0*(Pz0*AZt[J000]
		      +Pz1*AZt[J001]
		      +Pz2*AZtz[J000]
		      +Pz3*AZtz[J001]
		      )
		 +Pr1*(Pz0*AZt[J100]
		       +Pz1*AZt[J101]
		       +Pz2*AZtz[J100]
		       +Pz3*AZtz[J101]
		       )
		 +Pr2*(Pz0*AZrt[J000]
		       +Pz1*AZrt[J001]
		       +Pz2*AZrtz[J000]
		       +Pz3*AZrtz[J001]
		       )
		 +Pr3*(Pz0*AZrt[J100]
		       +Pz1*AZrt[J101]
		       +Pz2*AZrtz[J100]
		       +Pz3*AZrtz[J101]
		       )
		 )
	  +dPt3*(Pr0*(Pz0*AZt[J010]
		      +Pz1*AZt[J011]
		      +Pz2*AZtz[J010]
		      +Pz3*AZtz[J011]
		      )
		 +Pr1*(Pz0*AZt[J110]
		       +Pz1*AZt[J111]
		       +Pz2*AZtz[J110]
		       +Pz3*AZtz[J111]
		       )
		 +Pr2*(Pz0*AZrt[J010]
		       +Pz1*AZrt[J011]
		       +Pz2*AZrtz[J010]
		       +Pz3*AZrtz[J011]
		       )
		 +Pr3*(Pz0*AZrt[J110]
		       +Pz1*AZrt[J111]
		       +Pz2*AZrtz[J110]
		       +Pz3*AZrtz[J111]
		       )
		 )
	  );
  ESBgf	=Arz-Azr;
  r	=1./ESBgf;
  dy[0]	=r*(Azt-Atz);
  dy[1]	=1.;
  dy[2]	=r*(Atr-Art);
}

int ESSplineArtz(double R,double T,double Z,double *Ar,double *At,double *Az
		 ,double *Br,double *Bt,double *Bz)
{
  ESSetSplRgFZ(R,T,Z);

  *At	=(Pt0*(Pr0*(Pz0*AT[J000]
		    +Pz1*AT[J001]
		    +Pz2*ATz[J000]
		    +Pz3*ATz[J001]
		    )
	       +Pr1*(Pz0*AT[J100]
		     +Pz1*AT[J101]
		     +Pz2*ATz[J100]
		     +Pz3*ATz[J101]
		     )
	       +Pr2*(Pz0*ATr[J000]
		     +Pz1*ATr[J001]
		     +Pz2*ATrz[J000]
		     +Pz3*ATrz[J001]
		     )
	       +Pr3*(Pz0*ATr[J100]
		     +Pz1*ATr[J101]
		     +Pz2*ATrz[J100]
		     +Pz3*ATrz[J101]
		     )
	       )
	  +Pt1*(Pr0*(Pz0*AT[J010]
		     +Pz1*AT[J011]
		     +Pz2*ATz[J010]
		     +Pz3*ATz[J011]
		     )
		+Pr1*(Pz0*AT[J110]
		      +Pz1*AT[J111]
		      +Pz2*ATz[J110]
		      +Pz3*ATz[J111]
		      )
		+Pr2*(Pz0*ATr[J010]
		      +Pz1*ATr[J011]
		      +Pz2*ATrz[J010]
		      +Pz3*ATrz[J011]
		      )
		+Pr3*(Pz0*ATr[J110]
		      +Pz1*ATr[J111]
		      +Pz2*ATrz[J110]
		      +Pz3*ATrz[J111]
		      )
		)
	  +Pt2*(Pr0*(Pz0*ATt[J000]
		     +Pz1*ATt[J001]
		     +Pz2*ATtz[J000]
		     +Pz3*ATtz[J001]
		     )
		+Pr1*(Pz0*ATt[J100]
		      +Pz1*ATt[J101]
		      +Pz2*ATtz[J100]
		      +Pz3*ATtz[J101]
		      )
		  +Pr2*(Pz0*ATrt[J000]
			+Pz1*ATrt[J001]
			+Pz2*ATrtz[J000]
			+Pz3*ATrtz[J001]
			)
		  +Pr3*(Pz0*ATrt[J100]
			+Pz1*ATrt[J101]
			+Pz2*ATrtz[J100]
			+Pz3*ATrtz[J101]
			)
		  )
	    +Pt3*(Pr0*(Pz0*ATt[J010]
		       +Pz1*ATt[J011]
		       +Pz2*ATtz[J010]
		       +Pz3*ATtz[J011]
		       )
		  +Pr1*(Pz0*ATt[J110]
			+Pz1*ATt[J111]
			+Pz2*ATtz[J110]
			+Pz3*ATtz[J111]
			)
		  +Pr2*(Pz0*ATrt[J010]
			+Pz1*ATrt[J011]
			+Pz2*ATrtz[J010]
			+Pz3*ATrtz[J011]
			)
		  +Pr3*(Pz0*ATrt[J110]
			+Pz1*ATrt[J111]
			+Pz2*ATrtz[J110]
			+Pz3*ATrtz[J111]
			)
		  )
	    );
  
  *Ar	=(Pt0*(Pr0*(Pz0*AR[J000]
		    +Pz1*AR[J001]
		    +Pz2*ARz[J000]
		    +Pz3*ARz[J001]
		    )
		 +Pr1*(Pz0*AR[J100]
		       +Pz1*AR[J101]
		       +Pz2*ARz[J100]
		       +Pz3*ARz[J101]
		       )
		 +Pr2*(Pz0*ARr[J000]
		       +Pz1*ARr[J001]
		       +Pz2*ARrz[J000]
		       +Pz3*ARrz[J001]
		       )
		 +Pr3*(Pz0*ARr[J100]
		       +Pz1*ARr[J101]
		       +Pz2*ARrz[J100]
		       +Pz3*ARrz[J101]
		       )
		 )
	    +Pt1*(Pr0*(Pz0*AR[J010]
		       +Pz1*AR[J011]
		       +Pz2*ARz[J010]
		       +Pz3*ARz[J011]
		       )
		  +Pr1*(Pz0*AR[J110]
			+Pz1*AR[J111]
			+Pz2*ARz[J110]
			+Pz3*ARz[J111]
			)
		  +Pr2*(Pz0*ARr[J010]
			+Pz1*ARr[J011]
			+Pz2*ARrz[J010]
			+Pz3*ARrz[J011]
			)
		  +Pr3*(Pz0*ARr[J110]
			+Pz1*ARr[J111]
			+Pz2*ARrz[J110]
			+Pz3*ARrz[J111]
			)
		  )
	    +Pt2*(Pr0*(Pz0*ARt[J000]
		       +Pz1*ARt[J001]
		       +Pz2*ARtz[J000]
		       +Pz3*ARtz[J001]
		       )
		  +Pr1*(Pz0*ARt[J100]
			+Pz1*ARt[J101]
			+Pz2*ARtz[J100]
			+Pz3*ARtz[J101]
			)
		  +Pr2*(Pz0*ARrt[J000]
			+Pz1*ARrt[J001]
			+Pz2*ARrtz[J000]
			+Pz3*ARrtz[J001]
			)
		  +Pr3*(Pz0*ARrt[J100]
			+Pz1*ARrt[J101]
			+Pz2*ARrtz[J100]
			+Pz3*ARrtz[J101]
			)
		  )
	    +Pt3*(Pr0*(Pz0*ARt[J010]
		       +Pz1*ARt[J011]
		       +Pz2*ARtz[J010]
		       +Pz3*ARtz[J011]
		       )
		  +Pr1*(Pz0*ARt[J110]
			+Pz1*ARt[J111]
			+Pz2*ARtz[J110]
			+Pz3*ARtz[J111]
			)
		  +Pr2*(Pz0*ARrt[J010]
			+Pz1*ARrt[J011]
			+Pz2*ARrtz[J010]
			+Pz3*ARrtz[J011]
			)
		  +Pr3*(Pz0*ARrt[J110]
			+Pz1*ARrt[J111]
			+Pz2*ARrtz[J110]
			+Pz3*ARrtz[J111]
			)
		  )
	    );
  
  *Az	=(Pt0*(Pr0*(Pz0*AZ[J000]
		      +Pz1*AZ[J001]
		      +Pz2*AZz[J000]
		      +Pz3*AZz[J001]
		      )
		 +Pr1*(Pz0*AZ[J100]
		       +Pz1*AZ[J101]
		       +Pz2*AZz[J100]
		       +Pz3*AZz[J101]
		       )
		 +Pr2*(Pz0*AZr[J000]
		       +Pz1*AZr[J001]
		       +Pz2*AZrz[J000]
		       +Pz3*AZrz[J001]
		       )
		 +Pr3*(Pz0*AZr[J100]
		       +Pz1*AZr[J101]
		       +Pz2*AZrz[J100]
		       +Pz3*AZrz[J101]
		       )
		 )
	    +Pt1*(Pr0*(Pz0*AZ[J010]
		       +Pz1*AZ[J011]
		       +Pz2*AZz[J010]
		       +Pz3*AZz[J011]
		       )
		  +Pr1*(Pz0*AZ[J110]
			+Pz1*AZ[J111]
			+Pz2*AZz[J110]
			+Pz3*AZz[J111]
			)
		  +Pr2*(Pz0*AZr[J010]
			+Pz1*AZr[J011]
			+Pz2*AZrz[J010]
			+Pz3*AZrz[J011]
			)
		  +Pr3*(Pz0*AZr[J110]
			+Pz1*AZr[J111]
			+Pz2*AZrz[J110]
			+Pz3*AZrz[J111]
			)
		  )
	    +Pt2*(Pr0*(Pz0*AZt[J000]
		       +Pz1*AZt[J001]
		       +Pz2*AZtz[J000]
		       +Pz3*AZtz[J001]
		       )
		  +Pr1*(Pz0*AZt[J100]
			+Pz1*AZt[J101]
			+Pz2*AZtz[J100]
			+Pz3*AZtz[J101]
			)
		  +Pr2*(Pz0*AZrt[J000]
			+Pz1*AZrt[J001]
			+Pz2*AZrtz[J000]
			+Pz3*AZrtz[J001]
			)
		  +Pr3*(Pz0*AZrt[J100]
			+Pz1*AZrt[J101]
			+Pz2*AZrtz[J100]
			+Pz3*AZrtz[J101]
			)
		  )
	    +Pt3*(Pr0*(Pz0*AZt[J010]
		       +Pz1*AZt[J011]
		       +Pz2*AZtz[J010]
		       +Pz3*AZtz[J011]
		       )
		  +Pr1*(Pz0*AZt[J110]
			+Pz1*AZt[J111]
			+Pz2*AZtz[J110]
			+Pz3*AZtz[J111]
			)
		  +Pr2*(Pz0*AZrt[J010]
			+Pz1*AZrt[J011]
			+Pz2*AZrtz[J010]
			+Pz3*AZrtz[J011]
			)
		  +Pr3*(Pz0*AZrt[J110]
			+Pz1*AZrt[J111]
			+Pz2*AZrtz[J110]
			+Pz3*AZrtz[J111]
			)
		  )
	    );


  *Bt	=(Pt0*(Pr0*(Pz0*BT[J000]
		    +Pz1*BT[J001]
		    +Pz2*BTz[J000]
		    +Pz3*BTz[J001]
		    )
	       +Pr1*(Pz0*BT[J100]
		     +Pz1*BT[J101]
		     +Pz2*BTz[J100]
		     +Pz3*BTz[J101]
		     )
	       +Pr2*(Pz0*BTr[J000]
		     +Pz1*BTr[J001]
		     +Pz2*BTrz[J000]
		     +Pz3*BTrz[J001]
		     )
	       +Pr3*(Pz0*BTr[J100]
		     +Pz1*BTr[J101]
		     +Pz2*BTrz[J100]
		     +Pz3*BTrz[J101]
		     )
	       )
	  +Pt1*(Pr0*(Pz0*BT[J010]
		     +Pz1*BT[J011]
		     +Pz2*BTz[J010]
		     +Pz3*BTz[J011]
		     )
		+Pr1*(Pz0*BT[J110]
		      +Pz1*BT[J111]
		      +Pz2*BTz[J110]
		      +Pz3*BTz[J111]
		      )
		+Pr2*(Pz0*BTr[J010]
		      +Pz1*BTr[J011]
		      +Pz2*BTrz[J010]
		      +Pz3*BTrz[J011]
		      )
		+Pr3*(Pz0*BTr[J110]
		      +Pz1*BTr[J111]
		      +Pz2*BTrz[J110]
		      +Pz3*BTrz[J111]
		      )
		)
	  +Pt2*(Pr0*(Pz0*BTt[J000]
		     +Pz1*BTt[J001]
		     +Pz2*BTtz[J000]
		     +Pz3*BTtz[J001]
		     )
		+Pr1*(Pz0*BTt[J100]
		      +Pz1*BTt[J101]
		      +Pz2*BTtz[J100]
		      +Pz3*BTtz[J101]
		      )
		  +Pr2*(Pz0*BTrt[J000]
			+Pz1*BTrt[J001]
			+Pz2*BTrtz[J000]
			+Pz3*BTrtz[J001]
			)
		  +Pr3*(Pz0*BTrt[J100]
			+Pz1*BTrt[J101]
			+Pz2*BTrtz[J100]
			+Pz3*BTrtz[J101]
			)
		  )
	    +Pt3*(Pr0*(Pz0*BTt[J010]
		       +Pz1*BTt[J011]
		       +Pz2*BTtz[J010]
		       +Pz3*BTtz[J011]
		       )
		  +Pr1*(Pz0*BTt[J110]
			+Pz1*BTt[J111]
			+Pz2*BTtz[J110]
			+Pz3*BTtz[J111]
			)
		  +Pr2*(Pz0*BTrt[J010]
			+Pz1*BTrt[J011]
			+Pz2*BTrtz[J010]
			+Pz3*BTrtz[J011]
			)
		  +Pr3*(Pz0*BTrt[J110]
			+Pz1*BTrt[J111]
			+Pz2*BTrtz[J110]
			+Pz3*BTrtz[J111]
			)
		  )
	    );
  
  *Br	=(Pt0*(Pr0*(Pz0*BR[J000]
		    +Pz1*BR[J001]
		    +Pz2*BRz[J000]
		    +Pz3*BRz[J001]
		    )
		 +Pr1*(Pz0*BR[J100]
		       +Pz1*BR[J101]
		       +Pz2*BRz[J100]
		       +Pz3*BRz[J101]
		       )
		 +Pr2*(Pz0*BRr[J000]
		       +Pz1*BRr[J001]
		       +Pz2*BRrz[J000]
		       +Pz3*BRrz[J001]
		       )
		 +Pr3*(Pz0*BRr[J100]
		       +Pz1*BRr[J101]
		       +Pz2*BRrz[J100]
		       +Pz3*BRrz[J101]
		       )
		 )
	    +Pt1*(Pr0*(Pz0*BR[J010]
		       +Pz1*BR[J011]
		       +Pz2*BRz[J010]
		       +Pz3*BRz[J011]
		       )
		  +Pr1*(Pz0*BR[J110]
			+Pz1*BR[J111]
			+Pz2*BRz[J110]
			+Pz3*BRz[J111]
			)
		  +Pr2*(Pz0*BRr[J010]
			+Pz1*BRr[J011]
			+Pz2*BRrz[J010]
			+Pz3*BRrz[J011]
			)
		  +Pr3*(Pz0*BRr[J110]
			+Pz1*BRr[J111]
			+Pz2*BRrz[J110]
			+Pz3*BRrz[J111]
			)
		  )
	    +Pt2*(Pr0*(Pz0*BRt[J000]
		       +Pz1*BRt[J001]
		       +Pz2*BRtz[J000]
		       +Pz3*BRtz[J001]
		       )
		  +Pr1*(Pz0*BRt[J100]
			+Pz1*BRt[J101]
			+Pz2*BRtz[J100]
			+Pz3*BRtz[J101]
			)
		  +Pr2*(Pz0*BRrt[J000]
			+Pz1*BRrt[J001]
			+Pz2*BRrtz[J000]
			+Pz3*BRrtz[J001]
			)
		  +Pr3*(Pz0*BRrt[J100]
			+Pz1*BRrt[J101]
			+Pz2*BRrtz[J100]
			+Pz3*BRrtz[J101]
			)
		  )
	    +Pt3*(Pr0*(Pz0*BRt[J010]
		       +Pz1*BRt[J011]
		       +Pz2*BRtz[J010]
		       +Pz3*BRtz[J011]
		       )
		  +Pr1*(Pz0*BRt[J110]
			+Pz1*BRt[J111]
			+Pz2*BRtz[J110]
			+Pz3*BRtz[J111]
			)
		  +Pr2*(Pz0*BRrt[J010]
			+Pz1*BRrt[J011]
			+Pz2*BRrtz[J010]
			+Pz3*BRrtz[J011]
			)
		  +Pr3*(Pz0*BRrt[J110]
			+Pz1*BRrt[J111]
			+Pz2*BRrtz[J110]
			+Pz3*BRrtz[J111]
			)
		  )
	    );
  *Bz	=(Pt0*(Pr0*(Pz0*BZ[J000]
		      +Pz1*BZ[J001]
		      +Pz2*BZz[J000]
		      +Pz3*BZz[J001]
		      )
		 +Pr1*(Pz0*BZ[J100]
		       +Pz1*BZ[J101]
		       +Pz2*BZz[J100]
		       +Pz3*BZz[J101]
		       )
		 +Pr2*(Pz0*BZr[J000]
		       +Pz1*BZr[J001]
		       +Pz2*BZrz[J000]
		       +Pz3*BZrz[J001]
		       )
		 +Pr3*(Pz0*BZr[J100]
		       +Pz1*BZr[J101]
		       +Pz2*BZrz[J100]
		       +Pz3*BZrz[J101]
		       )
		 )
	    +Pt1*(Pr0*(Pz0*BZ[J010]
		       +Pz1*BZ[J011]
		       +Pz2*BZz[J010]
		       +Pz3*BZz[J011]
		       )
		  +Pr1*(Pz0*BZ[J110]
			+Pz1*BZ[J111]
			+Pz2*BZz[J110]
			+Pz3*BZz[J111]
			)
		  +Pr2*(Pz0*BZr[J010]
			+Pz1*BZr[J011]
			+Pz2*BZrz[J010]
			+Pz3*BZrz[J011]
			)
		  +Pr3*(Pz0*BZr[J110]
			+Pz1*BZr[J111]
			+Pz2*BZrz[J110]
			+Pz3*BZrz[J111]
			)
		  )
	    +Pt2*(Pr0*(Pz0*BZt[J000]
		       +Pz1*BZt[J001]
		       +Pz2*BZtz[J000]
		       +Pz3*BZtz[J001]
		       )
		  +Pr1*(Pz0*BZt[J100]
			+Pz1*BZt[J101]
			+Pz2*BZtz[J100]
			+Pz3*BZtz[J101]
			)
		  +Pr2*(Pz0*BZrt[J000]
			+Pz1*BZrt[J001]
			+Pz2*BZrtz[J000]
			+Pz3*BZrtz[J001]
			)
		  +Pr3*(Pz0*BZrt[J100]
			+Pz1*BZrt[J101]
			+Pz2*BZrtz[J100]
			+Pz3*BZrtz[J101]
			)
		  )
	    +Pt3*(Pr0*(Pz0*BZt[J010]
		       +Pz1*BZt[J011]
		       +Pz2*BZtz[J010]
		       +Pz3*BZtz[J011]
		       )
		  +Pr1*(Pz0*BZt[J110]
			+Pz1*BZt[J111]
			+Pz2*BZtz[J110]
			+Pz3*BZtz[J111]
			)
		  +Pr2*(Pz0*BZrt[J010]
			+Pz1*BZrt[J011]
			+Pz2*BZrtz[J010]
			+Pz3*BZrtz[J011]
			)
		  +Pr3*(Pz0*BZrt[J110]
			+Pz1*BZrt[J111]
			+Pz2*BZrtz[J110]
			+Pz3*BZrtz[J111]
			)
		  )
	    );
  return(0);
}

void ESMFLineSplStep(nq,x,yeq,dy)
     int *nq;
     double *x,*yeq,*dy;
{
  double r;

  if(ESSetSplRgFZ(yeq[0],*x,yeq[2])){
    printf("??? Out of spline grid%c\n",7);
  }
  
  ESBgf	=(Pt0*(Pr0*(Pz0*BT[J000]
		    +Pz1*BT[J001]
		    +Pz2*BTz[J000]
		    +Pz3*BTz[J001]
		    )
	       +Pr1*(Pz0*BT[J100]
		     +Pz1*BT[J101]
		     +Pz2*BTz[J100]
		     +Pz3*BTz[J101]
		     )
	       +Pr2*(Pz0*BTr[J000]
		     +Pz1*BTr[J001]
		     +Pz2*BTrz[J000]
		     +Pz3*BTrz[J001]
		     )
	       +Pr3*(Pz0*BTr[J100]
		     +Pz1*BTr[J101]
		     +Pz2*BTrz[J100]
		     +Pz3*BTrz[J101]
		     )
	       )
	  +Pt1*(Pr0*(Pz0*BT[J010]
		     +Pz1*BT[J011]
		     +Pz2*BTz[J010]
		     +Pz3*BTz[J011]
		     )
		+Pr1*(Pz0*BT[J110]
		      +Pz1*BT[J111]
		      +Pz2*BTz[J110]
		      +Pz3*BTz[J111]
		      )
		+Pr2*(Pz0*BTr[J010]
		      +Pz1*BTr[J011]
		      +Pz2*BTrz[J010]
		      +Pz3*BTrz[J011]
		      )
		+Pr3*(Pz0*BTr[J110]
		      +Pz1*BTr[J111]
		      +Pz2*BTrz[J110]
		      +Pz3*BTrz[J111]
		      )
		)
	  +Pt2*(Pr0*(Pz0*BTt[J000]
		     +Pz1*BTt[J001]
		     +Pz2*BTtz[J000]
		     +Pz3*BTtz[J001]
		     )
		+Pr1*(Pz0*BTt[J100]
		      +Pz1*BTt[J101]
		      +Pz2*BTtz[J100]
		      +Pz3*BTtz[J101]
		      )
		+Pr2*(Pz0*BTrt[J000]
		      +Pz1*BTrt[J001]
		      +Pz2*BTrtz[J000]
		      +Pz3*BTrtz[J001]
		      )
		+Pr3*(Pz0*BTrt[J100]
		      +Pz1*BTrt[J101]
		      +Pz2*BTrtz[J100]
		      +Pz3*BTrtz[J101]
		      )
		)
	  +Pt3*(Pr0*(Pz0*BTt[J010]
		     +Pz1*BTt[J011]
		     +Pz2*BTtz[J010]
		     +Pz3*BTtz[J011]
		     )
		+Pr1*(Pz0*BTt[J110]
		      +Pz1*BTt[J111]
		      +Pz2*BTtz[J110]
		      +Pz3*BTtz[J111]
		      )
		+Pr2*(Pz0*BTrt[J010]
		      +Pz1*BTrt[J011]
		      +Pz2*BTrtz[J010]
		      +Pz3*BTrtz[J011]
		      )
		+Pr3*(Pz0*BTrt[J110]
		      +Pz1*BTrt[J111]
		      +Pz2*BTrtz[J110]
		      +Pz3*BTrtz[J111]
		      )
		)
	  );
  r	=yeq[0]/ESBgf;
  dy[0]	=r*(Pt0*(Pr0*(Pz0*BR[J000]
		      +Pz1*BR[J001]
		      +Pz2*BRz[J000]
		      +Pz3*BRz[J001]
		      )
		 +Pr1*(Pz0*BR[J100]
		       +Pz1*BR[J101]
		       +Pz2*BRz[J100]
		       +Pz3*BRz[J101]
		       )
		 +Pr2*(Pz0*BRr[J000]
		       +Pz1*BRr[J001]
		       +Pz2*BRrz[J000]
		       +Pz3*BRrz[J001]
		       )
		 +Pr3*(Pz0*BRr[J100]
		       +Pz1*BRr[J101]
		       +Pz2*BRrz[J100]
		       +Pz3*BRrz[J101]
		       )
		 )
	    +Pt1*(Pr0*(Pz0*BR[J010]
		       +Pz1*BR[J011]
		       +Pz2*BRz[J010]
		       +Pz3*BRz[J011]
		       )
		  +Pr1*(Pz0*BR[J110]
			+Pz1*BR[J111]
			+Pz2*BRz[J110]
			+Pz3*BRz[J111]
			)
		  +Pr2*(Pz0*BRr[J010]
			+Pz1*BRr[J011]
			+Pz2*BRrz[J010]
			+Pz3*BRrz[J011]
			)
		  +Pr3*(Pz0*BRr[J110]
			+Pz1*BRr[J111]
			+Pz2*BRrz[J110]
			+Pz3*BRrz[J111]
			)
		  )
	    +Pt2*(Pr0*(Pz0*BRt[J000]
		       +Pz1*BRt[J001]
		       +Pz2*BRtz[J000]
		       +Pz3*BRtz[J001]
		       )
		  +Pr1*(Pz0*BRt[J100]
			+Pz1*BRt[J101]
			+Pz2*BRtz[J100]
			+Pz3*BRtz[J101]
			)
		  +Pr2*(Pz0*BRrt[J000]
			+Pz1*BRrt[J001]
			+Pz2*BRrtz[J000]
			+Pz3*BRrtz[J001]
			)
		  +Pr3*(Pz0*BRrt[J100]
			+Pz1*BRrt[J101]
			+Pz2*BRrtz[J100]
			+Pz3*BRrtz[J101]
			)
		  )
	    +Pt3*(Pr0*(Pz0*BRt[J010]
		       +Pz1*BRt[J011]
		       +Pz2*BRtz[J010]
		       +Pz3*BRtz[J011]
		       )
		  +Pr1*(Pz0*BRt[J110]
			+Pz1*BRt[J111]
			+Pz2*BRtz[J110]
			+Pz3*BRtz[J111]
			)
		  +Pr2*(Pz0*BRrt[J010]
			+Pz1*BRrt[J011]
			+Pz2*BRrtz[J010]
			+Pz3*BRrtz[J011]
			)
		  +Pr3*(Pz0*BRrt[J110]
			+Pz1*BRrt[J111]
			+Pz2*BRrtz[J110]
			+Pz3*BRrtz[J111]
			)
		  )
	    );
  
  dy[1]	=1.;
  dy[2]	=r*(Pt0*(Pr0*(Pz0*BZ[J000]
		      +Pz1*BZ[J001]
		      +Pz2*BZz[J000]
		      +Pz3*BZz[J001]
		      )
		 +Pr1*(Pz0*BZ[J100]
		       +Pz1*BZ[J101]
		       +Pz2*BZz[J100]
		       +Pz3*BZz[J101]
		       )
		 +Pr2*(Pz0*BZr[J000]
		       +Pz1*BZr[J001]
		       +Pz2*BZrz[J000]
		       +Pz3*BZrz[J001]
		       )
		 +Pr3*(Pz0*BZr[J100]
		       +Pz1*BZr[J101]
		       +Pz2*BZrz[J100]
		       +Pz3*BZrz[J101]
		       )
		 )
	    +Pt1*(Pr0*(Pz0*BZ[J010]
		       +Pz1*BZ[J011]
		       +Pz2*BZz[J010]
		       +Pz3*BZz[J011]
		       )
		  +Pr1*(Pz0*BZ[J110]
			+Pz1*BZ[J111]
			+Pz2*BZz[J110]
			+Pz3*BZz[J111]
			)
		  +Pr2*(Pz0*BZr[J010]
			+Pz1*BZr[J011]
			+Pz2*BZrz[J010]
			+Pz3*BZrz[J011]
			)
		  +Pr3*(Pz0*BZr[J110]
			+Pz1*BZr[J111]
			+Pz2*BZrz[J110]
			+Pz3*BZrz[J111]
			)
		  )
	    +Pt2*(Pr0*(Pz0*BZt[J000]
		       +Pz1*BZt[J001]
		       +Pz2*BZtz[J000]
		       +Pz3*BZtz[J001]
		       )
		  +Pr1*(Pz0*BZt[J100]
			+Pz1*BZt[J101]
			+Pz2*BZtz[J100]
			+Pz3*BZtz[J101]
			)
		  +Pr2*(Pz0*BZrt[J000]
			+Pz1*BZrt[J001]
			+Pz2*BZrtz[J000]
			+Pz3*BZrtz[J001]
			)
		  +Pr3*(Pz0*BZrt[J100]
			+Pz1*BZrt[J101]
			+Pz2*BZrtz[J100]
			+Pz3*BZrtz[J101]
			)
		  )
	    +Pt3*(Pr0*(Pz0*BZt[J010]
		       +Pz1*BZt[J011]
		       +Pz2*BZtz[J010]
		       +Pz3*BZtz[J011]
		       )
		  +Pr1*(Pz0*BZt[J110]
			+Pz1*BZt[J111]
			+Pz2*BZtz[J110]
			+Pz3*BZtz[J111]
			)
		  +Pr2*(Pz0*BZrt[J010]
			+Pz1*BZrt[J011]
			+Pz2*BZrtz[J010]
			+Pz3*BZrtz[J011]
			)
		  +Pr3*(Pz0*BZrt[J110]
			+Pz1*BZrt[J111]
			+Pz2*BZrtz[J110]
			+Pz3*BZrtz[J111]
			)
		  )
	    );
}

int ESSetSplRZ(double r,double z,int n)
{
  int iR,kT,jZ;
  int i,j,kT1;
  double x,X,xx,XX;
  
  iR	=irCell;
  kT	=ktCell;
  jZ	=jzCell;
  while(r1Cell <= r){
    rCell	=r1Cell;
    r1Cell	+=drSpl;
    irCell++;
  }
  while(rCell > r){
    r1Cell	=rCell;
    rCell	-=drSpl;
    irCell--;
  }
  tCell=gf1T[n];
  t1Cell=tCell+dtSpl;
  ktCell= n == ESNt ? 0 : n; 
  while(z1Cell <= z){
    zCell	=z1Cell;
    z1Cell	+=dzSpl;
    jzCell++;
  }
  while(zCell > z){
    z1Cell	=zCell;
    zCell	-=dzSpl;
    jzCell--;
  }
  if(iR != irCell || kT != ktCell || jZ != jzCell){
    iR	=irCell;
    kT	=ktCell;
    jZ	=jzCell;
   
    kT1=kT+1;
    if(kT1 == ESNt)
      kT1=0;
    if(jZ >= Nz2Spl[kT] || jZ >= Nz2Spl[kT1]){
      printf("??? jZ=%3d >= Nz2[%3d]=%3d Nz2[%3d]=%3d - Out of spline grid%c\n"
	     ,jZ,kT,Nz2Spl[kT],kT1,Nz2Spl[kT1],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT kT1\n",r,n,z);
      return(1);
    }
    j	=jZ-Nz1Spl[kT];
    if(j < 0 || iR >= Nr2Spl[kT][jZ]){
      if(j < 0)
	printf("??? jZ=%3d < Nz1[%3d]=%3d - Out of spline grid%c\n"
	       ,jZ,kT,Nz1Spl[kT],7);
      if(iR >= Nr2Spl[kT][jZ])
	printf("??? iR >= %3d < Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT,jZ,Nr2Spl[kT][jZ],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT\n",r,n,z);
      return(1);
    }
    i	=iR-Nr1Spl[kT][jZ];
    if(i < 0){
      printf("??? iR < %3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	     ,iR,kT,jZ,Nr1Spl[kT][jZ],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT jZ\n",r,n,z);
      return(1);
    }
    J000	=k2kSpl[kT][j][i];
    J100	=J000+1;
    i	=iR-Nr1Spl[kT][jZ+1];
    if(i < 0 || iR >= Nr2Spl[kT][jZ+1]){
      if(i < 0)
	printf("??? iR=%3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT,jZ+1,Nr1Spl[kT][jZ+1],7);
      if(iR >= Nr2Spl[kT][jZ+1])
	printf("??? iR=%3d >= Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT,jZ+1,Nr2Spl[kT][jZ+1],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT jZ+1\n",r,n,z);
      return(1);
    }
    J001	=k2kSpl[kT][j+1][i];
    J101	=J001+1;
    j	=jZ-Nz1Spl[kT1];
    if(j < 0 || iR >= Nr2Spl[kT1][jZ]){
      if(j < 0)
	printf("??? jZ=%3d < Nz1[%3d]=%3d - Out of spline grid%c\n"
	       ,jZ,kT1,Nz1Spl[kT1],7);
      if(iR >= Nr2Spl[kT1][jZ])
	printf("??? iR >= %3d < Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ,Nr2Spl[kT1][jZ],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT1 jZ\n",r,n,z);
      return(1);
    }
    i	=iR-Nr1Spl[kT1][jZ];
    if(i < 0){
      if(i < 0)
	printf("??? iR=%3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ,Nr1Spl[kT1][jZ],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT1 jZ\n",r,n,z);
      return(1);
    }
    J010	=k2kSpl[kT1][j][i];
    J110	=J010+1;
    i	=iR-Nr1Spl[kT1][jZ+1];
    if(i < 0 || iR >= Nr2Spl[kT1][jZ+1]){
      if(i < 0)
	printf("??? iR=%3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ+1,Nr1Spl[kT1][jZ+1],7);
      if(iR >= Nr2Spl[kT1][jZ+1])
	printf("??? iR >= %3d < Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ+1,Nr2Spl[kT1][jZ+1],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT1 jZ+1\n",r,n,z);
      return(1);
    }
    J011	=k2kSpl[kT1][j+1][i];
    J111	=J011+1;
  }
  x	=(r-rCell)*rhrCell;
  xx	=x*x;
  X	=1.-x;
  XX	=X*X;
  Pr0=XX*(3.-2.*X);
  Pr1=xx*(3.-2.*x);
  Pr2=XX*x*drSpl;
  Pr3=-xx*X*drSpl;

  x	=(z-zCell)*rhzCell;
  xx	=x*x;
  X	=1.-x;
  XX	=X*X;
  Pz0=XX*(3.-2.*X);
  Pz1=xx*(3.-2.*x);
  Pz2=XX*x*dzSpl;
  Pz3=-xx*X*dzSpl;
  return(0);
}

int ESSetSplDerRZ(double r,double z,int n)
{
  int iR,kT,jZ;
  int i,j,kT1;
  double x,X,xx,XX,Xx;
  
  iR	=irCell;
  kT	=ktCell;
  jZ	=jzCell;
  while(r1Cell <= r){
    rCell	=r1Cell;
    r1Cell	+=drSpl;
    irCell++;
  }
  while(rCell > r){
    r1Cell	=rCell;
    rCell	-=drSpl;
    irCell--;
  }
  tCell=gf1T[n];
  t1Cell=tCell+dtSpl;
  ktCell= n == ESNt ? 0 : n; 
  while(z1Cell <= z){
    zCell	=z1Cell;
    z1Cell	+=dzSpl;
    jzCell++;
  }
  while(zCell > z){
    z1Cell	=zCell;
    zCell	-=dzSpl;
    jzCell--;
  }
  if(iR != irCell || kT != ktCell || jZ != jzCell){
    iR	=irCell;
    kT	=ktCell;
    jZ	=jzCell;
   
    kT1=kT+1;
    if(kT1 == ESNt)
      kT1=0;
 
    if(jZ >= Nz2Spl[kT] || jZ >= Nz2Spl[kT1]){
      printf("??? jZ=%3d >= Nz2[%3d]=%3d Nz2[%3d]=%3d - Out of spline grid%c\n"
	     ,jZ,kT,Nz2Spl[kT],kT1,Nz2Spl[kT1],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT kT1\n",r,n,z);
      return(1);
    }
    j	=jZ-Nz1Spl[kT];
    if(j < 0 || iR >= Nr2Spl[kT][jZ]){
      if(j < 0)
	printf("??? jZ=%3d < Nz1[%3d]=%3d - Out of spline grid%c\n"
	       ,jZ,kT,Nz1Spl[kT],7);
      if(iR >= Nr2Spl[kT][jZ])
	printf("??? iR >= %3d < Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT,jZ,Nr2Spl[kT][jZ],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT\n",r,n,z);
      return(1);
    }
    i	=iR-Nr1Spl[kT][jZ];
    if(i < 0){
      printf("??? iR < %3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	     ,iR,kT,jZ,Nr1Spl[kT][jZ],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT jZ\n",r,n,z);
      return(1);
    }
    J000	=k2kSpl[kT][j][i];
    J100	=J000+1;
    i	=iR-Nr1Spl[kT][jZ+1];
    if(i < 0 || iR >= Nr2Spl[kT][jZ+1]){
      if(i < 0)
	printf("??? iR=%3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT,jZ+1,Nr1Spl[kT][jZ+1],7);
      if(iR >= Nr2Spl[kT][jZ+1])
	printf("??? iR=%3d >= Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT,jZ+1,Nr2Spl[kT][jZ+1],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT jZ+1\n",r,n,z);
      return(1);
    }
    J001	=k2kSpl[kT][j+1][i];
    J101	=J001+1;
    j	=jZ-Nz1Spl[kT1];
    if(j < 0 || iR >= Nr2Spl[kT1][jZ]){
      if(j < 0)
	printf("??? jZ=%3d < Nz1[%3d]=%3d - Out of spline grid%c\n"
	       ,jZ,kT1,Nz1Spl[kT1],7);
      if(iR >= Nr2Spl[kT1][jZ])
	printf("??? iR >= %3d < Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ,Nr2Spl[kT1][jZ],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT1 jZ\n",r,n,z);
      return(1);
    }
    i	=iR-Nr1Spl[kT1][jZ];
    if(i < 0){
      if(i < 0)
	printf("??? iR=%3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ,Nr1Spl[kT1][jZ],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT1 jZ\n",r,n,z);
      return(1);
    }
    J010	=k2kSpl[kT1][j][i];
    J110	=J010+1;
    i	=iR-Nr1Spl[kT1][jZ+1];
    if(i < 0 || iR >= Nr2Spl[kT1][jZ+1]){
      if(i < 0)
	printf("??? iR=%3d < Nr1[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ+1,Nr1Spl[kT1][jZ+1],7);
      if(iR >= Nr2Spl[kT1][jZ+1])
	printf("??? iR >= %3d < Nr2[%3d][%2d]=%3d - Out of spline grid%c\n"
	       ,iR,kT1,jZ+1,Nr2Spl[kT1][jZ+1],7);
      printf("??? r=%10.3e n=%d z=%10.3e kT1 jZ+1\n",r,n,z);
      return(1);
    }
    J011	=k2kSpl[kT1][j+1][i];
    J111	=J011+1;
  }

  x	=(r-rCell)*rhrCell;
  xx	=x*x;
  X	=1.-x;
  XX	=X*X;
  Pr0=XX*(3.-2.*X);
  Pr1=xx*(3.-2.*x);
  Pr2=XX*x*drSpl;
  Pr3=-xx*X*drSpl;
  Xx	=2.*x*X;
  dPr1	=3.*rhrCell*Xx;
  dPr0	=-dPr1;
  dPr2	=XX-Xx;
  dPr3	=xx-Xx;

  x	=(z-zCell)*rhzCell;
  xx	=x*x;
  X	=1.-x;
  XX	=X*X;
  Pz0=XX*(3.-2.*X);
  Pz1=xx*(3.-2.*x);
  Pz2=XX*x*dzSpl;
  Pz3=-xx*X*dzSpl;
  Xx	=2.*x*X;
  dPz1	=3.*rhzCell*Xx;
  dPz0	=-dPz1;
  dPz2	=XX-Xx;
  dPz3	=xx-Xx;
  return(0);
}

int ESSpl2ArBtAtAz(double *Ar, double *Ar1t, double *Bt,
		   double *At, double *Az, double *Az1t,
		   double R, double Z, int n)
{
  double r;

  if(ESSetSplRZ(R,Z,n)){
    printf("??? Out of spline grid%c\n",7);
  }
  *Ar	=(Pr0*(Pz0*AR[J000]
	       +Pz1*AR[J001]
	       +Pz2*ARz[J000]
	       +Pz3*ARz[J001]
	       )
	  +Pr1*(Pz0*AR[J100]
		+Pz1*AR[J101]
		+Pz2*ARz[J100]
		+Pz3*ARz[J101]
		)
	  +Pr2*(Pz0*ARr[J000]
		+Pz1*ARr[J001]
		+Pz2*ARrz[J000]
		+Pz3*ARrz[J001]
		)
	  +Pr3*(Pz0*ARr[J100]
		+Pz1*ARr[J101]
		+Pz2*ARrz[J100]
		+Pz3*ARrz[J101]
		)
	  );
  
  *Ar1t	=(Pr0*(Pz0*ARt[J000]
	       +Pz1*ARt[J001]
	       +Pz2*ARtz[J000]
	       +Pz3*ARtz[J001]
	       )
	  +Pr1*(Pz0*ARt[J100]
		+Pz1*ARt[J101]
		+Pz2*ARtz[J100]
		+Pz3*ARtz[J101]
		)
	  +Pr2*(Pz0*ARrt[J000]
		+Pz1*ARrt[J001]
		+Pz2*ARrtz[J000]
		+Pz3*ARrtz[J001]
		)
	  +Pr3*(Pz0*ARrt[J100]
		+Pz1*ARrt[J101]
		+Pz2*ARrtz[J100]
		+Pz3*ARrtz[J101]
		)
	  );

  *Bt	=(Pr0*(Pz0*BT[J000]
	       +Pz1*BT[J001]
	       +Pz2*BTz[J000]
	       +Pz3*BTz[J001]
	       )
	  +Pr1*(Pz0*BT[J100]
		+Pz1*BT[J101]
		+Pz2*BTz[J100]
		+Pz3*BTz[J101]
		)
	  +Pr2*(Pz0*BTr[J000]
		+Pz1*BTr[J001]
		+Pz2*BTrz[J000]
		+Pz3*BTrz[J001]
		)
	  +Pr3*(Pz0*BTr[J100]
		+Pz1*BTr[J101]
		+Pz2*BTrz[J100]
		+Pz3*BTrz[J101]
		)
	  );

  *At	=(Pr0*(Pz0*AT[J000]
	       +Pz1*AT[J001]
	       +Pz2*ATz[J000]
	       +Pz3*ATz[J001]
	       )
	  +Pr1*(Pz0*AT[J100]
		+Pz1*AT[J101]
		+Pz2*ATz[J100]
		+Pz3*ATz[J101]
		)
	  +Pr2*(Pz0*ATr[J000]
		+Pz1*ATr[J001]
		+Pz2*ATrz[J000]
		+Pz3*ATrz[J001]
		)
	  +Pr3*(Pz0*ATr[J100]
		+Pz1*ATr[J101]
		+Pz2*ATrz[J100]
		+Pz3*ATrz[J101]
		)
	  );

  *Az	=(Pr0*(Pz0*AZ[J000]
	       +Pz1*AZ[J001]
	       +Pz2*AZz[J000]
	       +Pz3*AZz[J001]
	       )
	  +Pr1*(Pz0*AZ[J100]
		+Pz1*AZ[J101]
		+Pz2*AZz[J100]
		+Pz3*AZz[J101]
		)
	  +Pr2*(Pz0*AZr[J000]
		+Pz1*AZr[J001]
		+Pz2*AZrz[J000]
		+Pz3*AZrz[J001]
		)
	  +Pr3*(Pz0*AZr[J100]
		+Pz1*AZr[J101]
		+Pz2*AZrz[J100]
		+Pz3*AZrz[J101]
		)
	  );

  *Az1t	=(Pr0*(Pz0*AZt[J000]
	       +Pz1*AZt[J001]
	       +Pz2*AZtz[J000]
	       +Pz3*AZtz[J001]
	       )
	  +Pr1*(Pz0*AZt[J100]
		+Pz1*AZt[J101]
		+Pz2*AZtz[J100]
		+Pz3*AZtz[J101]
		)
	  +Pr2*(Pz0*AZrt[J000]
		+Pz1*AZrt[J001]
		+Pz2*AZrtz[J000]
		+Pz3*AZrtz[J001]
		)
	  +Pr3*(Pz0*AZrt[J100]
		+Pz1*AZrt[J101]
		+Pz2*AZrtz[J100]
		+Pz3*AZrtz[J101]
		)
	  );

  return(0);
}

int ESSpl2BrBtBz(double *Ar, double *Az,double *Br, double *Bt, double *Bz, double R, double Z, int n)
{
  double r;

  if(ESSetSplRZ(R,Z,n)){
    printf("??? Out of spline grid%c\n",7);
  }
  
  *Ar	=(Pr0*(Pz0*AR[J000]
	       +Pz1*AR[J001]
	       +Pz2*ARz[J000]
	       +Pz3*ARz[J001]
	       )
	  +Pr1*(Pz0*AR[J100]
		+Pz1*AR[J101]
		+Pz2*ARz[J100]
		+Pz3*ARz[J101]
		)
	  +Pr2*(Pz0*ARr[J000]
		+Pz1*ARr[J001]
		+Pz2*ARrz[J000]
		+Pz3*ARrz[J001]
		)
	  +Pr3*(Pz0*ARr[J100]
		+Pz1*ARr[J101]
		+Pz2*ARrz[J100]
		+Pz3*ARrz[J101]
		)
	  );

  *Az	=(Pr0*(Pz0*AZ[J000]
	       +Pz1*AZ[J001]
	       +Pz2*AZz[J000]
	       +Pz3*AZz[J001]
	       )
	  +Pr1*(Pz0*AZ[J100]
		+Pz1*AZ[J101]
		+Pz2*AZz[J100]
		+Pz3*AZz[J101]
		)
	  +Pr2*(Pz0*AZr[J000]
		+Pz1*AZr[J001]
		+Pz2*AZrz[J000]
		+Pz3*AZrz[J001]
		)
	  +Pr3*(Pz0*AZr[J100]
		+Pz1*AZr[J101]
		+Pz2*AZrz[J100]
		+Pz3*AZrz[J101]
		)
	  );
  
  *Br	=(Pr0*(Pz0*BR[J000]
	       +Pz1*BR[J001]
	       +Pz2*BRz[J000]
	       +Pz3*BRz[J001]
	       )
	  +Pr1*(Pz0*BR[J100]
		+Pz1*BR[J101]
		+Pz2*BRz[J100]
		+Pz3*BRz[J101]
		)
	  +Pr2*(Pz0*BRr[J000]
		+Pz1*BRr[J001]
		+Pz2*BRrz[J000]
		+Pz3*BRrz[J001]
		)
	  +Pr3*(Pz0*BRr[J100]
		+Pz1*BRr[J101]
		+Pz2*BRrz[J100]
		+Pz3*BRrz[J101]
		)
	  );

  *Bt	=(Pr0*(Pz0*BT[J000]
	       +Pz1*BT[J001]
	       +Pz2*BTz[J000]
	       +Pz3*BTz[J001]
	       )
	  +Pr1*(Pz0*BT[J100]
		+Pz1*BT[J101]
		+Pz2*BTz[J100]
		+Pz3*BTz[J101]
		)
	  +Pr2*(Pz0*BTr[J000]
		+Pz1*BTr[J001]
		+Pz2*BTrz[J000]
		+Pz3*BTrz[J001]
		)
	  +Pr3*(Pz0*BTr[J100]
		+Pz1*BTr[J101]
		+Pz2*BTrz[J100]
		+Pz3*BTrz[J101]
		)
	  );

  *Bz	=(Pr0*(Pz0*BZ[J000]
	       +Pz1*BZ[J001]
	       +Pz2*BZz[J000]
	       +Pz3*BZz[J001]
	       )
	  +Pr1*(Pz0*BZ[J100]
		+Pz1*BZ[J101]
		+Pz2*BZz[J100]
		+Pz3*BZz[J101]
		)
	  +Pr2*(Pz0*BZr[J000]
		+Pz1*BZr[J001]
		+Pz2*BZrz[J000]
		+Pz3*BZrz[J001]
		)
	  +Pr3*(Pz0*BZr[J100]
		+Pz1*BZr[J101]
		+Pz2*BZrz[J100]
		+Pz3*BZrz[J101]
		)
	  );

  return(0);
}

int ESSpl2BrBtBzBrrBtrBzr(double *Br, double *Bt, double *Bz,
		 double *Brr, double *Btr, double *Bzr,
		 double *Brz, double *Btz, double *Bzz,
		 double R, double Z, int n)
{
  double r,p0,p1,p2,p3;

  if(ESSetSplDerRZ(R,Z,n)){
    printf("??? Out of spline grid%c\n",7);
  }
  
  p0	=Pz0*BR[J000]+Pz1*BR[J001]+Pz2*BRz[J000]+Pz3*BRz[J001];
  p1	=Pz0*BR[J100]+Pz1*BR[J101]+Pz2*BRz[J100]+Pz3*BRz[J101];
  p2	=Pz0*BRr[J000]+Pz1*BRr[J001]+Pz2*BRrz[J000]+Pz3*BRrz[J001];
  p3	=Pz0*BRr[J100]+Pz1*BRr[J101]+Pz2*BRrz[J100]+Pz3*BRrz[J101];
  
  *Br	=Pr0*p0+Pr1*p1+Pr2*p2+Pr3*p3;
  *Brr	=dPr0*p0+dPr1*p1+dPr2*p2+dPr3*p3;
  *Brz	=(Pr0*(dPz0*BR[J000]+dPz1*BR[J001]+dPz2*BRz[J000]+dPz3*BRz[J001])
	  +Pr1*(dPz0*BR[J100]+dPz1*BR[J101]+dPz2*BRz[J100]+dPz3*BRz[J101])
	  +Pr2*(dPz0*BRr[J000]+dPz1*BRr[J001]+dPz2*BRrz[J000]+dPz3*BRrz[J001])
	  +Pr3*(dPz0*BRr[J100]+dPz1*BRr[J101]+dPz2*BRrz[J100]+dPz3*BRrz[J101])
	  );

  p0	=Pz0*BT[J000]+Pz1*BT[J001]+Pz2*BTz[J000]+Pz3*BTz[J001];
  p1	=Pz0*BT[J100]+Pz1*BT[J101]+Pz2*BTz[J100]+Pz3*BTz[J101];
  p2	=Pz0*BTr[J000]+Pz1*BTr[J001]+Pz2*BTrz[J000]+Pz3*BTrz[J001];
  p3	=Pz0*BTr[J100]+Pz1*BTr[J101]+Pz2*BTrz[J100]+Pz3*BTrz[J101];
  
  *Bt	=Pr0*p0+Pr1*p1+Pr2*p2+Pr3*p3;
  *Btr	=dPr0*p0+dPr1*p1+dPr2*p2+dPr3*p3;
  *Btz	=(Pr0*(dPz0*BT[J000]+dPz1*BT[J001]+dPz2*BTz[J000]+dPz3*BTz[J001])
	  +Pr1*(dPz0*BT[J100]+dPz1*BT[J101]+dPz2*BTz[J100]+dPz3*BTz[J101])
	  +Pr2*(dPz0*BTr[J000]+dPz1*BTr[J001]+dPz2*BTrz[J000]+dPz3*BTrz[J001])
	  +Pr3*(dPz0*BTr[J100]+dPz1*BTr[J101]+dPz2*BTrz[J100]+dPz3*BTrz[J101])
	  );

  p0	=Pz0*BZ[J000]+Pz1*BZ[J001]+Pz2*BZz[J000]+Pz3*BZz[J001];
  p1	=Pz0*BZ[J100]+Pz1*BZ[J101]+Pz2*BZz[J100]+Pz3*BZz[J101];
  p2	=Pz0*BZr[J000]+Pz1*BZr[J001]+Pz2*BZrz[J000]+Pz3*BZrz[J001];
  p3	=Pz0*BZr[J100]+Pz1*BZr[J101]+Pz2*BZrz[J100]+Pz3*BZrz[J101];
  
  *Bz	=Pr0*p0+Pr1*p1+Pr2*p2+Pr3*p3;
  *Bzr	=dPr0*p0+dPr1*p1+dPr2*p2+dPr3*p3; 
  *Bzz	=(Pr0*(dPz0*BZ[J000]+dPz1*BZ[J001]+dPz2*BZz[J000]+dPz3*BZz[J001])
	  +Pr1*(dPz0*BZ[J100]+dPz1*BZ[J101]+dPz2*BZz[J100]+dPz3*BZz[J101])
	  +Pr2*(dPz0*BZr[J000]+dPz1*BZr[J001]+dPz2*BZrz[J000]+dPz3*BZrz[J001])
	  +Pr3*(dPz0*BZr[J100]+dPz1*BZr[J101]+dPz2*BZrz[J100]+dPz3*BZrz[J101])
	  );
  return(0);
}

int ESInit3DSpl()
{
  int n,jn,j,NEp1,NEp2,NEp3,nT,nS;
  static int iT=0,iT1,kT=0,iz=0,kT1,ir;
  double s,*px,*py,*pr,*pz,*pr1,*pz1,gF;
  double R[7],Z[7],X[7],Y[7];
  double rC[7],zC[7],xC[7],yC[7];

  if(TFCNt == 0){
    return(-1);
  }
  if(NDm < TFCNt){
    ESReInitDmTFC();
  }
  TFCNt1=TFCNt+1;
  n	=TFCNt*TFCLt+1;
  NEp1	=n+5;
  NEp2	=NEp1+n;
  NEp3	=NEp2+TFCNp1;
  if(SplNt < ESNt1){
    ESReInitEnvTFC();
  }
  ESInitDmTFC();
  ESGetTFCAddr(&gf0TFC,&gfL,&csL,&snL);

  nTP=TFCNt*TFCNp1;

  Rmn=rTFC[jRmnTFC[0]];
  Rmx=rTFC[jRmxTFC[0]];
  Zmn=zTFC[jZmnTFC[0]];
  Zmx=zTFC[jZmxTFC[0]];
  for(n=1; n < TFCNt; n++){
    jn	=n*TFCNp1;
    j	=jn+jRmnTFC[n];
    if(Rmn > rTFC[j]){
      Rmn=rTFC[j];
    }
    j	=jn+jRmxTFC[n];
    if(Rmx < rTFC[j]){
      Rmx=rTFC[j];
    }
    j	=jn+jZmnTFC[n];
    if(Zmn > zTFC[j]){
      Zmn=zTFC[j];
    }
    j	=jn+jZmxTFC[n];
    if(Zmx < zTFC[j]){
      Zmx=zTFC[j];
    }
  }
  R[0]	=EZcr2*(Rmn+Rmx);
  R[1]	=Rmn;
  R[2]	=Rmx;
  Z[0]	=EZcr2*(Zmn+Zmx);
  Z[1]	=Zmn;
  Z[2]	=Zmx;
  X[0]	=R[0];
  Y[0]	=0.;
#ifndef Tbl_SPL
  nT	=ESNt;
  nrSpl1=nrSpl+1;
  nzSpl1=nzSpl+1;
  if(NzSpl < nzSpl1){
    ESReInitSplGrid();
  }
  Rmn	=R[1];
  Rmx	=R[2];
  Zmn	=Z[1];
  Zmx	=Z[2];
  drSpl	=(Rmx-Rmn)/nrSpl;
  dtSpl	=gf1T[1];
  dzSpl	=(Zmx-Zmn)/nzSpl;
  rhrCell=1./drSpl;
  rhtCell=1./dtSpl;
  rhzCell=1./dzSpl;

  iz	=(Z[0]-Zmn)/dzSpl;
  R[0]	=sqrt(X[0]*X[0]+Y[0]*Y[0]);
  ir	=(R[0]-Rmn)/drSpl;
  
  if(iz < 0){
    iz	=0;
  }
  if(iz > nzSpl){
    iz	=nzSpl;
  }
  if(X[0] == 0. && Y[0] == 0.){
    gF	=0.;
  }
  else{
    if(fabs(Y[0]) < fabs(X[0])){
      gF	=atan(Y[0]/X[0]);
      if(X[0] < 0.){
	gF	+=EZcgp;
      }
      else{
	if(Y[0] < 0.){
	  gF	+=EZc2gp;
	}
      }
    }
    else{
      gF	=EZcr2*EZcgp-atan(X[0]/Y[0]);
      if(Y[0] < 0.){
	gF	+=EZcgp;
      }
    }
  }
  
  kT	=gF*ESNt/EZc2gp;
  if(kT < 0)
    kT=0;
  if(kT >= ESNt)
    kT -=ESNt;
  kT1	=kT+1;
  if(kT1 >= ESNt)
    kT1 -=ESNt;
 
  gF	=gf0TFC[0];
  while(gF+gfL[1] <= gf1T[kT]){
    gF	+=gfL[1];
  }
  while(gF > gf1T[kT]){
    gF	-=gfL[1];
  }
  gF	=gf1T[kT]-gF+gf0TFC[0];
  iT1	=iT+1;
  while(iT < TFCNt && gf0TFC[iT1] <= gF){
    iT1++;
    iT++;
  }
  while(iT && gf0TFC[iT] > gF){
    iT1--;
    iT--;
  }
  if(iT1 == TFCNt){
    iT1	=0;
  }
  px	=xTFC+TFCNp1*iT;
  py	=yTFC+TFCNp1*iT;
  pr	=rTFC+TFCNp1*iT;
  pz	=zTFC+TFCNp1*iT;
  pr1	=rTFC+TFCNp1*iT1;
  pz1	=zTFC+TFCNp1*iT1;

  ESInitSplGrid();
  ESInitEnvTFC(Z[0]);
  R[3]	=r1TFC[iT];
  R[4]	=r2TFC[iT];
  R[5]	=r1TFC[iT1];
  R[6]	=r2TFC[iT1];
  Z[3]	=Z[0];
  Z[4]	=Z[0];
  Z[5]	=Z[0];
  Z[6]	=Z[0];
  X[1]	=-Rmx;
  X[2]	=Rmx;
  X[3]	=x1TFC[iT];
  X[4]	=x2TFC[iT];
  X[5]	=x1TFC[iT1];
  X[6]	=x2TFC[iT1];
  Y[1]	=-Rmx;
  Y[2]	=Rmx;
  Y[3]	=y1TFC[iT];
  Y[4]	=y2TFC[iT];
  Y[5]	=y1TFC[iT1];
  Y[6]	=y2TFC[iT1];
  n	=TFCNt*TFCLt+1;
#undef Ins_SPL1
#endif /*Tbl_SPL*/

  nT	=-1;
  kT1	=1;
#ifndef BnrI_SPL
#endif /*BnrI_SPL*/
  if(kT1){
    ESInitArtzBrtzSpl();
#ifndef BnrO_SPL
#endif /*BnrO_SPL*/
  }
#ifdef H
  rCell=Rmn;
  tCell=0.;
  zCell=Zmn;
  r1Cell=rCell+drSpl;
  t1Cell=tCell+dtSpl;
  z1Cell=zCell+dzSpl;
  irCell=0;
  ktCell=0;
  jzCell=0;
#endif
  rCell=Rmn-drSpl;
  tCell=0.;
  zCell=Zmn-dzSpl;
  r1Cell=rCell+drSpl;
  t1Cell=tCell+dtSpl;
  z1Cell=zCell+dzSpl;
  irCell=-1;
  ktCell=0;
  jzCell=-1;
  return(0);
}
#endif

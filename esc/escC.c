#include <math.h>
#include <stdio.h>
#include <string.h>
#include "esc.h"

#ifndef mzl_ESmain
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern int ESmemlev;
#endif

#ifndef mzl_1Dcore
#endif

#ifndef mzl_2Dcore
extern int ESNp,ESNp1;
extern int *ESnF,*ESmF,*ESkF,ESNf,ESNf1;
extern double *EZdra,*EZdrt,*EZdza,*EZdzt;
#endif

#ifndef mzl_3Dcore
extern int ESLt,ESNt,ESNt1;
extern double *gf1T,*csLT,*snLT,*cs1T,*sn1T;
#endif

#ifndef mzl_TFC
extern int TFCLt,TFCNt,TFCNp,TFCNp1;
extern double *xTFC,*yTFC,*zTFC,*gfTFC,*rTFC,ITFC;

double ESBgf;

static int NTP=0,nTP=0,NT=0,nT=0,LT=0;
static double *hTFC,*gf0TFC;
static double *gfL,*csL,*snL;
static double Rmn,Rmx,Zmn,Zmx,Hmn,Hmx;

static int NPncr=0,nPncr=0,nPsrf=1,KPncr=1,iPncr=0,LPncr;
static double *rPncr,*zPncr,*xPncr,*yPncr,gfPncr0,gfPncr1;

static double Bextz=0.,Bgf;

double Ax,Ay,Az,Bx,By,Bz;

#endif

#ifndef mzl_Lsode
static int lrw=68,liw=20;
static double *yeq,*rwork,*yh,*acor,*vtol,*rvtol;
static int iwork[20],istate,Neq=0,neq=3,itol=1,itask=4,mf=10,iopt=0;
static double ErrR=1e-8,ErrA=1e-12;

int ESReInitLsode()
{
  if(Neq < lrw+2*neq){
    lrw		=liw+16*neq;
    if(ReInitArray((void**)&rwork,Neq,lrw+3*neq,sizeof(double)) < 0){
      FailureAlarm((char*)rwork,"ESReInitLsode() - no memory for rwork");
      ESexit(0);
    }
    Neq	=lrw+2*neq;
    yh		=rwork+20;
    acor	=rwork+lrw-neq;
    vtol	=rwork+lrw;
    rvtol	=vtol+neq;
    yeq		=rvtol+neq;
    ESmemlev	|=0x00010000;
  }
  return(0);
}

int ESDeInitLsode()
{
  free(rwork);
  return(0);
}
#endif

#ifndef mzl_TFC
#ifndef stg_DataTFC
int ESGetTFCAddr(double **gf0,double **gfL0,double **csL0,double **snL0)
{
  *gf0=gf0TFC;
  *gfL0=gfL;
  *csL0=csL;
  *snL0=snL;
  return(0);
}

int ESReInit1DTFC()
{
  if(NT < TFCNt){
    if(ReInitArray((void**)&gf0TFC,NT,TFCNt+1,sizeof(double)) < 0){
      FailureAlarm((char*)gf0TFC,"ESReInit1DTFC() - no memory for gf0TFC");
      ESexit(0);
    }
    NT =TFCNt+1;
    ESmemlev |=0x00004000;
  }
  return(0);
}

int ESDeInit1DTFC()
{
  if(NT){
    free(gf0TFC);
    NT=0;
  }
  return(0);
}

int ESInit1LTFC()
{
  int L;
  double t,dt;

  dt	=EZc2gp/abs(ESLt);
  for(L=0; L < abs(ESLt); L++){
    t		=L*dt;
    gfL[L]	=t;
    csL[L]	=cos(t);
    snL[L]	=sin(t);
  }
  gfL[L]	=EZc2gp;
  csL[L]	=1.;
  snL[L]	=0.;
  return(0);
}

int ESInit3DTFC()
{
  int i,n,in;
  double cgp05,x,y,z,gf,r,cs,sn,X,Y;

  cgp05	=EZcr2*EZcgp;
  Rmn=1e+20;
  Rmx=-Rmn;
  Zmn=Rmn;
  Zmx=Rmx;
  Hmn=Rmn;
  Hmx=Rmx;
  in=0;
  for(n=0; n < TFCNt; n++){
    x	=0.;
    y	=0.;
    for(i=0; i < TFCNp; i++){
      x	+=xTFC[in];
      y	+=yTFC[in];
      in++;
    }
    in++;
    if(x == 0. && y == 0.){
      gf=0.;
    }
    else{
      if(fabs(y) < fabs(x)){
	gf	=atan(y/x);
	if(x < 0.){
	  gf	+=EZcgp;
	}
      }
      else{
	gf	=cgp05-atan(x/y);
	if(y < 0.){
	  gf	+=EZcgp;
	}
      }
    }
    gf0TFC[n]	=n && gf < gf0TFC[n-1] ? gf+EZc2gp : gf;
  }
  in	=0;
  for(n=0; n < TFCNt; n++){
    cs	=cos(gf0TFC[n]);
    sn	=sin(gf0TFC[n]);
    for(i=0; i < TFCNp1; i++){
      x	=xTFC[in];
      y	=yTFC[in];
      z	=zTFC[in];
      r	=sqrt(x*x+y*y);
      X	=x*cs+y*sn;
      Y	=-x*sn+y*cs;
      if(fabs(Y) < fabs(X)){
	gf	=atan(Y/X);
	if(X < 0.){
	  if(Y < 0.)
	    gf	-=EZcgp;
	  else
	    gf	+=EZcgp;
	}
      }
      else{
	gf	=cgp05-atan(X/Y);
	if(Y < 0.){
	  gf	-=EZcgp;
	}
      }
      gfTFC[in]	=gf+gf0TFC[n];
      rTFC[in]	=r;
      if(Rmn > r){
	Rmn=r;
      }
      if(Rmx < r){
	Rmx=r;
      }
      if(Zmn > z){
	Zmn=z;
      }
      if(Zmx < z){
	Zmx=z;
      }
      in++;
    }
  }
  for(n=0; n <TFCNt; n++){
    in=TFCNp1*n;
    gf=0;
    for(i=0; i < TFCNp1; i++){
      gf =gfTFC[in]-gf0TFC[n];
      hTFC[in]=rTFC[in]*sin(gf);
      if(Hmn > hTFC[in]){
	Hmn=hTFC[in];
      }
      if(Hmx < hTFC[in]){
	Hmx=hTFC[in];
      }
      in++;
    }
  }
  gf0TFC[n]	=gf0TFC[0]+EZc2gp/abs(ESLt);
  return(0);
}

int ESReInit3DTFC()
{
  if(NTP < nTP){
    if(ReInitArray((void**)&hTFC,NTP,nTP,sizeof(double)) < 0){
      FailureAlarm((char*)hTFC,"ESReInit3DTFC() - no memory for hTFC");
      ESexit(0);
    }
    if(ReInitArray((void**)&gfTFC,NTP,nTP,sizeof(double)) < 0){
      FailureAlarm((char*)gfTFC,"ESReInit3DTFC() - no memory for gfTFC");
      ESexit(0);
    }
    if(ReInitArray((void**)&rTFC,NTP,nTP,sizeof(double)) < 0){
      FailureAlarm((char*)rTFC,"ESReInit3DTFC() - no memory for rTFC");
      ESexit(0);
    }
    if(ReInitArray((void**)&zTFC,NTP,nTP,sizeof(double)) < 0){
      FailureAlarm((char*)zTFC,"ESReInit3DTFC() - no memory for zTFC");
      ESexit(0);
    }
    if(ReInitArray((void**)&yTFC,NTP,nTP,sizeof(double)) < 0){
      FailureAlarm((char*)yTFC,"ESReInit3DTFC() - no memory for yTFC");
      ESexit(0);
    }
    if(ReInitArray((void**)&xTFC,NTP,nTP,sizeof(double)) < 0){
      FailureAlarm((char*)xTFC,"ESReInit3DTFC() - no memory for xTFC");
      ESexit(0);
    }
    NTP =nTP;
    ESmemlev |=0x00002000;
  }
  return(0);
}

int ESDeInit3DTFC()
{
  if(NTP){
    free(xTFC);
    free(yTFC);
    free(zTFC);
    free(rTFC);
    free(gfTFC);
    free(hTFC);
    NTP=0;
  }
  return(0);
}

int ESInitTFC()
{
  int i;
  static int iT=0;
  double s,*pr,*ph,*pz;

  if(ESLt == 0){
    return(0);
  }
#ifndef AscI_TFC
#undef Ins1_TFC
  if(nTP == 0){
    return(-1);
  }
  TFCLt=abs(ESLt);
  TFCNp1=TFCNp+1;;
  TFCNt=nTP/TFCNp1;

  if(NT < TFCNt){
    ESReInit1DTFC();
  }
  if(NTP < nTP){
    ESReInit3DTFC();
  }
#undef Ins2_TFC
#endif /*AscI_TFC*/
  ESInit3DTFC();
#ifndef Tbl_TFC
  if(iT < 0)
    iT=0;
  if(iT >= TFCNt)
    iT=TFCNt-1;
  pr	=rTFC+TFCNp1*iT;
  ph	=hTFC+TFCNp1*iT;
  pz	=zTFC+TFCNp1*iT;
#endif /*Tbl_TFC*/
  if(LT < TFCLt+1){
    ESReInit1LTFC();
  }
  ESInit1LTFC();
  return(0);
}
#endif

#ifndef stg_VacMF
int ESReInit1LTFC()
{
  if(LT < TFCLt+1){
    if(ReInitArray((void**)&gfL,LT,TFCLt+1,sizeof(double)) < 0){
      FailureAlarm((char*)gfL,"ESReInit1LTFC() - no memory for gfL");
      ESexit(0);
    }
    if(ReInitArray((void**)&csL,LT,TFCLt+1,sizeof(double)) < 0){
      FailureAlarm((char*)csL,"ESReInit1LTFC() - no memory for csL");
      ESexit(0);
    }
    if(ReInitArray((void**)&snL,LT,TFCLt+1,sizeof(double)) < 0){
      FailureAlarm((char*)snL,"ESReInit1LTFC() - no memory for snL");
      ESexit(0);
    }
    LT =TFCLt+1;
    ESmemlev |=0x00008000;
  }
  return(0);
}

int ESDeInit1LTFC()
{
  if(LT){
    free(snL);
    free(csL);
    free(gfL);
    LT=0;
  }
  return(0);
}

int ESReInitPncr()
{
  if(NPncr < nPncr){
    if(ReInitArray((void**)&rPncr,NPncr,nPncr,sizeof(double)) < 0){
      FailureAlarm((char*)rPncr,"ESReInitPncr() - no memory for rPncr");
      ESexit(0);
    }
    if(ReInitArray((void**)&zPncr,NPncr,nPncr,sizeof(double)) < 0){
      FailureAlarm((char*)zPncr,"ESReInitPncr() - no memory for zPncr");
      ESexit(0);
    }
    if(ReInitArray((void**)&xPncr,NPncr,nPncr,sizeof(double)) < 0){
      FailureAlarm((char*)xPncr,"ESReInitPncr() - no memory for xPncr");
      ESexit(0);
    }
    if(ReInitArray((void**)&yPncr,NPncr,nPncr,sizeof(double)) < 0){
      FailureAlarm((char*)yPncr,"ESReInitPncr() - no memory for yPncr");
      ESexit(0);
    }
    NPncr	=nPncr;
    ESmemlev	|=0x00020000;
  }
  return(0);
}

int ESDeInitPncr()
{
  free(zPncr);
  free(rPncr);
  return(0);
}

int ESGetPncrAddr(int *nPncr0,int *nPsrf0,double **rPncr0,double **zPncr0,
		  double **xPncr0,double **yPncr0)
{
  *nPncr0=nPncr;
  *nPsrf0=nPsrf;
  *rPncr0=rPncr;
  *zPncr0=zPncr;
  *xPncr0=xPncr;
  *yPncr0=yPncr;
  return(0);
}


#endif
#endif

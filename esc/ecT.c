#include <math.h>
#include <stdio.h>
#include <time.h>

#define AAASAS
extern clock_t clock(void);

extern int ESkNormP,ESkNormJ,ESnDtPr,ESnDtCr;

extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern int ESNa,ESNa1,ESNp,ESNp1,ESnAP;
extern int ESMp,ESMp1,ES2Mp1,ESnMp;
extern int ESFp,ESFp1,ESnAF;
extern double *ESsa,*ESpa;
extern double *ESsr,*ESsz;
extern double ESa0;
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
extern double *ECb,*ECb1a,*ECb2a;
extern double *EScs1,*ESsn1;

extern double R0,Z0;
extern double *ESsa,*ESsb,*ESaR0;

extern int ESEqSolvFl,ESEqSolvRadC,ESEqSolvInPr,ESEqSolvInCr;
extern int ESEqSolvIt;
extern double ESEqSolvTol;

extern double *ESg22c,*ESg22s,*ESg22c1,*ESg22s1,*ESg22c2,*ESg22s2;
extern double *ESr22,*ESr222,*ESr22c,*ESr22s,*ESr22c2,*ESr22s2;
extern double *ESg12c,*ESg12s,*ESg11c,*ESg11s;
extern double *ESg12c2,*ESg12s2,*ESg11c2,*ESg11s2;
extern double *ESLc,*ESLs,*ESVc,*ESVs,*ESDPc,*ESDPs;
extern double *ESLc1,*ESLs1,*ESVc1,*ESVs1,*ESDPc1,*ESDPs1;
extern double *ESLc2,*ESLs2,*ESVc2,*ESVs2,*ESDPc2,*ESDPs2;

extern double *rcT,*rcT1a,*rcT2a,*rsT,*rsT1a,*rsT2a;
extern double *dR0T,*dZ0T,*dbT,*dbT1a,*dbT2a;
extern double *drcT,*drcT1a,*drcT2a,*drsT,*drsT1a,*drsT2a;

extern double ESBt,ESRBt,ESRext,ESgbext,ESgbN;
extern double *ESgY0,*ESgY02a,*ESaY,ESgY2a1a;
extern double ESgFPlV;

#ifdef MSE
extern int ESXNa,ESXNa1;
extern double *ESXsxD,*ESXsrD;
extern double *ESXsx,*ESXsr,*ESXsj,*ESXsp,*ESXsq,*ESXsi;
#endif

extern int ESnBL,ESiAx,ESaJNorm;
extern double ESgFPlVx,EZd_vb,EZgd_vb;
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

extern double EZga0,EZga1,EZga2,EZga3;
extern double *EZrcs,*EZrsn,*EZd1rcs,*EZd1rsn,*EZd2rcs,*EZd2rsn,*EZz0,*EZd1z0,*EZd2z0;
extern double *rcE,*rsE,*rcE1,*rsE1,*rcE2,*rsE2,*bE,*bE1,*bE2,R0E,Z0E;

extern double eps_rel,eps_abs;
extern double Epsrvb,Eps_noise;

extern double*Fcs,*Fsn,*EZxinf,*EZyinf,*d1yinf,*d2yinf,*EZxgt,*EZygt;

extern int M0Mm;

extern double *EZd1drc,*EZd1drs;

extern double *EZgper,*EZgpei,*EZdgper,*EZdgpei;
extern double *rCc,*rCs,*rSc,*rSs,*drCc,*drCs,*drSc,*drSs;
extern double *bCc,*bSc,*dbCc,*dbSc;

static double *gpi,*Yr,*Yi,*dgpi,*dYr,*dYi;

static double gBext,gBN;
static double ECVol,rR0;
static double RK1Dh;

#ifdef AAASAS
static void (*peq1D)(int ,double *, double *,double*);
#endif

int ESRK4D1(double *Y, double *dY, double *X)
{
  double RKy0,RKy1,RKy2,RKy3,RKy4,x;

  RKy1	=RK1Dh*(*dY);
  RKy0	=(*Y)+EZcr2*RKy1;
  x	=(*X)+EZcr2*RK1Dh;
  peq1D(1,&x,&RKy0,&RKy2);
  RKy2	*=RK1Dh;
  RKy0	=(*Y)+EZcr2*RKy2;
  peq1D(1,&x,&RKy0,&RKy3);
  RKy3	*=RK1Dh;
  RKy0	=(*Y)+RKy3;
  (*X)	+=RK1Dh;
  peq1D(1,X,&RKy0,&RKy4);
  RKy4	*=RK1Dh;
  *Y	+=EZcr6*(RKy1+RKy4)+EZcr3*(RKy2+RKy3);
  peq1D(1,X,Y,dY);
  return(0);
}

int ES2Tr6gip(float *r0,float *EZz0,float *ip,float *rbt
	      ,float *xz,float *xc,float *p,float *gi)
{
  int k,i;
  double s,x;

  *r0		=R0*100.;
  *EZz0		=Z0*100.;
  *ip		=ESaJ[ESNa]*5e+6;
  *rbt		=ESRBt*100.;
  if(ESEqSolvInPr == 2){
    k		=1;
    for(i=0; i < ESNa; i++){
      x		=xz[k];
      ESSetSplDPr(x);
      splRDPr(&s,NULL,2);
      p[k]	=s*EZcrgm0*1e+6;
      k++;
    }
    x	=1.;
    ESSetSplDPr(x);
    splRDPr(&s,NULL,2);
    p[k]	=s*EZcrgm0*1e+6;
  }
  else{
    k		=1;
    for(i=0; i < ESNa; i++){
      x		=xz[k];
      ESSetSplA(x);
      splRA(&s,NULL,ESsp,ESsp2a);
      p[k]	=s*EZcrgm0*1e+6;
      k++;
    }
    p[k]	=ESsp[ESNa]*EZcrgm0*1e+6;
  }

  if(ESEqSolvInCr == 6 || ESEqSolvInCr == 7){
    k		=1;
    for(i=0; i < ESNa; i++){
      x		=xc[i];
      ESSetSplDCr(x);
      splRDCr(&s,NULL,7);
      gi[k]	=s;
      k++;
    }
    gi[k]	=ESgm[ESNa];
    gi[0]	=gi[1];
  }
  else{
    k		=1;
    for(i=0; i < ESNa; i++){
      x		=xc[i];
      ESSetSplA(x);
      splRA(&s,NULL,ESgm,ESgm2a);
      gi[k]	=s;
      k++;
    }
    gi[k]	=ESgm[ESNa];
    gi[0]	=gi[1];
  }

  printf("Ipl=%13.6e [MA]\n",ESaJ[ESNa]*5.);
  printf(" i %13s %13s %13s %13s\n","sqrt(Phi)","Phi","iota","p [MPa]");
  for(i=0; i < ESNa1; i++){
    printf("%2d %11.4e %11.4e %13.6e %13.6e\n"
	   ,i,ESsa[i],ESpa[i],ESgm[i],ESsp[i]*EZcrgm0);
  }

  return(0);
}

#ifdef MSE
int ESXProf()
{
  double x,r1t,r1a,z1a,z1t,gtt,D;
  int i,k,ki;

  ESXsx[ESNa]	=0.;
  ESXsr[ESNa]	=ESaR0[0];
  ESXsj[ESNa]	=ESjs[0];
  ESXsp[ESNa]	=ESsp[0];
  ESXsq[ESNa]	=ESsq[0];
  ESXsi[ESNa]	=0.;
  for(i=1; i < ESNa1; i++){
    /*--- tau=0 ---*/  
    ki	=i;
    x	=rcT[ki];
    r1t	=0.;
    r1a	=rcT1a[ki];
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      x		+=2.*rcT[ki]; 
      r1t	+=2.*k*rsT[ki];
      r1a	+=2.*rcT1a[ki]; 
    }
    ki		=i+ESNa;
    ESXsx[ki]	=x;
    ESXsr[ki]	=x+ESaR0[0];
    D		=ESXsr[ki]/ESaR0[0];
    ESXsj[ki]	=ESjs[i]/D+ESjp[i]*(D-1./D);
    ESXsp[ki]	=ESsp[i];
    ESXsq[ki]	=ESsq[i];
    z1a		=rsT1a[i];
    z1t		=ESsb[i];
    gtt		=r1t*r1t+z1t*z1t;
    D		=r1a*z1t-r1t*z1a;
    ESXsi[ki]	=sqrt(gtt)*ESsa[i]*ESLc[i]*ESgm[i]/(ESaR0[0]*D); 
 
    /*--- tau=Pi ---*/ 
    ki	=i;
    x	=rcT[ki];
    r1t	=0;
    r1a	=rcT1a[ki];
    D	=-2.;
    for(k=1; k < ESFp1; k++){
      ki	+=ESNa1;
      x		+=D*rcT[ki]; 
      r1t	+=D*k*rsT[ki];
      r1a	+=D*rcT1a[ki]; 
      D		*=-1.;
    }
    ki		=ESNa-i;
    ESXsx[ki]	=x;
    ESXsr[ki]	=x+ESaR0[0];
    D		=ESXsr[ki]/ESaR0[0];
    ESXsj[ki]	=ESjs[i]/D+ESjp[i]*(D-1./D);
    ESXsp[ki]	=ESsp[i];
    ESXsq[ki]	=ESsq[i];
    z1a	=rsT1a[i];
    z1t	=-ESsb[i];
    gtt	=r1t*r1t+z1t*z1t;
    D	=r1a*z1t-r1t*z1a;
    ESXsi[ki]=-sqrt(gtt)*ESsa[i]*ESLc[i]*ESgm[i]/(ESaR0[0]*D); 
  }
  for(i=0; i < ESXNa1; i++){
    ESXsxD[i]	=ESXsrD[i]-ESaR0[0];
  }
  return(0);
}
#endif

int EC1DRadCSolvHole()
{
  int i,j;
  double rV,rgF,rgY,rB,R2R0;
  double *gxF,*gxV,*gxY,*gxB;

  rR0	=1./R0;
  R2R0	=2.*R0*R0;
  
  gxF	=ESLs;
  gxV	=ESVs;
  gxY	=ESg22s;
  gxB	=ESDPs;

  /*Calculation of Radial coordinates and other parameters;*/
  {
    int i1;
    double h2sh,h12p2sh;
    double dV0,dV1,d2V0,d2V1;
    double ds0,ds1,d2s0,d2s1,S;
    double dP0,dP1,d2P0,d2P1,P;
    double dgF0,dgF1,d2gF0,d2gF1;
    double dgb0,dgb1,d2gb0,d2gb1;
    
    h12p2sh	=ESsa[1]-ESa0;
    h2sh	=EZcr2*h12p2sh;
    h12p2sh	=EZcr12*h12p2sh*h12p2sh;

    dgF0	=0.;
    ECVol	=0.;
    for(j=0; j < ESNp; j++){
      dgF0	+=log(ESsr[j]*rR0)*EScs1[j];
      ECVol	+=(ESsr[j]*ESsr[j]-R0*R0)*EScs1[j];
      S	+=(ESsr[j]-R0)*EScs1[j];
    }
    ESgF[0]	=dgF0*ESsb[0]*ESaF[0]/ESNp;
    ESqV[0]	=EZcr2*ESsb[0]*ECVol*rR0/ESNp;
    S	=ESsb[0]*rcT[ESNa1];

    i	=0;
    d2V1	=ESLc[0]+ESVc[0];
    dV1		=ESsa[0]*d2V1;
    d2V1	+=ESsa[0]*(ESLc1[0]+ESVc1[0]);
    ESqV1a[0]	=dV1;

    d2gF1	=ESaF[0]*ESLc[0]*rR0;
    dgF1	=ESsa[0]*d2gF1;
    d2gF1	+=ESsa[0]*(ESaF1a[0]*ESLc[0]+ESaF[0]*ESLc1[0])*rR0;
    ESgF1a[0]	=dgF1;
    ESqgF1a[0]	=d2gF1;

    ESgB[0]	=R2R0*ESsp[0]/ESFF[0];
    ESgB1a[0]	=0.;

    dgb1	=ESsp[0]*dV1;
    d2gb1	=ESsp[0]*d2V1;
    gBext	=ESsp[0]*ESqV[0];
   
    ds1		=ESsa[0]*ESDPc[0];
    d2s1	=ESDPc[0]+ESsa[0]*ESDPc1[0];
    S		=ESsb[0]*rcT[ESNa1];

    P		=ESjs[0]*ESLc[0]+ESjp[0]*ESVc[0];
    dP1		=0.;
    d2P1	=-ESjp[0]*S*P/ESg22c[0];
    ESgb[0]	=2.*S*ESjp[0]*rR0/(ESpa[0]*ESg22c[0]*P);
    ESgb1a[0]	=0.;
    P		=0.;
    i1	=0;
    for(i=1; i < ESNa1; i++){
      dV0	=dV1;
      d2V0	=d2V1;
      d2V1	=ESLc[i]+ESVc[i];
      dV1	=ESsa[i]*d2V1;
      d2V1	+=ESsa[i]*(ESLc1[i]+ESVc1[i]);
      ESqV1a[i]	=dV1;
      ESqV[i]	=ESqV[i1]+h2sh*(dV0+dV1)+h12p2sh*(d2V0-d2V1);

      dgF0	=dgF1;
      d2gF0	=d2gF1;
      d2gF1	=ESaF[i]*ESLc[i]*rR0;
      dgF1	=ESsa[i]*d2gF1;
      d2gF1	+=ESsa[i]*(ESaF1a[i]*ESLc[i]+ESaF[i]*ESLc1[i])*rR0;
      ESgF1a[i]	=dgF1;
      ESgF[i]	=ESgF[i1]+h2sh*(dgF0+dgF1)+h12p2sh*(d2gF0-d2gF1);

      ESgB[i]	=R2R0*ESsp[i]/ESFF[i];
      ESgB1a[i]	=(R2R0*ESsp1a[i]-ESgB[i]*ESFF1a[i])/ESFF[i];

      dgb0	=dgb1;
      d2gb0	=d2gb1;
      dgb1	=ESsp[i]*dV1;
      d2gb1	=ESsp1a[i]*dV1+ESsp[i]*d2V1;
      gBext	+=h2sh*(dgb0+dgb1)+h12p2sh*(d2gb0-d2gb1);

      ds0	=ds1;
      d2s0	=d2s1;
      ds1	=ESsa[i]*ESDPc[i];
      d2s1	=ESDPc[i]+ESsa[i]*ESDPc1[i];
      S		+=h2sh*(ds0+ds1)+h12p2sh*(d2s0-d2s1);

      dP0	=dP1;
      d2P0	=d2P1;
      dP1	=ESjp[i]*ESgY1a[i]*S;
      d2P1	=ESjp1a[i]*ESgY1a[i]*S+ESjp[i]*(ESgY2a[i]*S+ESgY1a[i]*ds1);
      P	+=h2sh*(dP0+dP1)+h12p2sh*(d2P0-d2P1);
      ESgb[i]	=-4.*P*rR0/(ESaJ[i]*ESaJ[i]);
      i1++;
    }
    ESgb1a[ESNa]	=-4.*rR0*dP1/(ESaJ[ESNa]*ESaJ[ESNa])
      -2.*ESgb[ESNa]*ESaJ1a[ESNa]/ESaJ[ESNa];
  }
  splAA(ESgb,ESgb1a,ESgb2a,ESgb1a,ESgb1a+ESNa);
  splA(ESgF,ESgF2a,ESgF1a,ESgF1a+ESNa);
  splA(ESgB,ESgB2a,ESgB1a,ESgB1a+ESNa);
  ESgFPlV	=ESgF[ESNa];
  if(ESiAx){
    rB	=ESnBL/(double)ESNa;
    rB	=rB*(rB+2.)*ESgFPlV/ESgFPlVx;
    EZgd_vb	=EZd_vb*sqrt(rB)-EZd_vb;
  }
  ECVol		=ESqV[ESNa];
  gBext		*=2.*ESRext*ESRext/(ECVol*ESFF[ESNa]);
  i	=ESNp1*ESNa;
  gBN		=100.*gBext*ESBt*EZcr2*(ESsr[i]-ESsr[i+ESNp/2])/(5.*ESaJ[ESNa]);
#ifdef XWIN__
  printf("<gb>_ext=%10.3e <gb>_N=%10.3e\n",gBext,gBN);
#endif

  i		=0;
  rB		=1./ECb[ESNa];
  ESsh[i]	=ECb[i]*rB;
  ESsh1a[i]	=ECb1a[i]*rB;
  ESsh2a[i]	=ECb2a[i]*rB;
  rV		=1./ECVol;
  ESqV[i]	=sqrt(ESqV[i]*rV);
  ESqV1a[i]	*=EZcr2*rV/ESqV[i];
  rgF		=1./ESgF[ESNa];
  ESqgF[i]	=sqrt(ESgF[i]*rgF);
  ESqgF1a[i]	=EZcr2*ESgF1a[i]*rgF/ESqgF[i];
  rgY		=1./ESgY[ESNa];
  ESqgY[i]	=0.;
  ESqgY1a[i]	=sqrt(EZcr2*ESgY2a[0]*rgY);
  for(i=1; i < ESNa1; i++){
    ESsh[i]	=ECb[i]*rB;
    ESsh1a[i]	=ECb1a[i]*rB;
    ESsh2a[i]	=ECb2a[i]*rB;

    ESqV[i]	=sqrt(ESqV[i]*rV);
    ESqV1a[i]	*=EZcr2*rV/ESqV[i];
    ESqgF[i]	=sqrt(ESgF[i]*rgF);
    ESqgF1a[i]	=EZcr2*ESgF1a[i]*rgF/ESqgF[i];
    ESqgY[i]	=sqrt(ESgY[i]*rgY);
    ESqgY1a[i]	=EZcr2*ESgY1a[i]*rgY/ESqgY[i];
  }
  splA(ESqV,ESqV2a,ESqV1a,ESqV1a+ESNa);
  splA(ESqgF,ESqgF2a,ESqgF1a,ESqgF1a+ESNa);
  splA(ESqgY,ESqgY2a,ESqgY1a,ESqgY1a+ESNa);
  for(i=0; i < ESNa1; i++){
    gxB[i]	=(ESsa[i]-ESsh[i])/ESsh1a[i];
    gxV[i]	=(ESsa[i]-ESqV[i])/ESqV1a[i];
    gxF[i]	=(ESsa[i]-ESqgF[i])/ESqgF1a[i];
    gxY[i]	=(ESsa[i]-ESqgY[i])/ESqgY1a[i];
  }
#ifdef MSE
  ESXProf();
#endif
  return(0);
}

int EC1DRadCSolv()
{
  int i;
  double rV,rgF,rgY,rB,R2R0;
  double *gxF,*gxV,*gxY,*gxB;

  if(ESa0 != 0.){
    EC1DRadCSolvHole();
    return(0);
  }

  rR0	=1./R0;
  R2R0	=2.*R0*R0;

  gxF	=ESLs;
  gxV	=ESVs;
  gxY	=ESg22s;
  gxB	=ESDPs;

  /*Calculation of Radial coordinates and other parameters;*/
  {
    int i1;
    double h2sh,h12p2sh;
    double dV0,dV1,d2V0,d2V1;
    double ds0,ds1,d2s0,d2s1,S;
    double dS0,dS1,d2S0,d2S1,aps;
    double dgF0,dgF1,d2gF0,d2gF1;
    double dgb0,dgb1,d2gb0,d2gb1;
    
    h2sh	=EZcr2*ESsa[1];
    h12p2sh	=EZcr12*ESpa[1];

    ESgF[0]	=0.;
    ESgF1a[0]	=0.;
    dgF1	=0.;
    d2gF1	=ESaF[0]*ESLc[0]*rR0;
    ESqgF1a[0]	=d2gF1;

    ESqV[0]	=0.;
    ESqV1a[0]	=0.;
    dV1		=0.;
    d2V1	=ESLc[0]+ESVc[0];
    ds1		=0.;
    d2s1	=0.;
    S		=0.;

    dS1		=ESDPc[0]*ESjs[0];
    ESgb[0]	=-2.*(ESsp2a[0]*ESDPc[0]+ESsp[0]*ESDPc1[0])/(dS1*dS1);
    ESgb1a[0]	=0.;
    dS1		=0.;
    d2S1	=0.;
    aps		=0.;

    dgb1	=0.;
    d2gb1	=ESsp[0]*d2V1;
    ESgB[0]	=R2R0*ESsp[0]/ESFF[0];
    ESgB1a[0]	=0.;
    gBext	=0.;
    i1	=0;
    for(i=1; i < ESNa1; i++){
      dV0	=dV1;
      d2V0	=d2V1;
      d2V1	=ESLc[i]+ESVc[i];
      dV1	=ESsa[i]*d2V1;
      d2V1	+=ESsa[i]*(ESLc1[i]+ESVc1[i]);
      ESqV1a[i]	=dV1;
      ESqV[i]	=ESqV[i1]+h2sh*(dV0+dV1)+h12p2sh*(d2V0-d2V1);

      ds0	=ds1;
      d2s0	=d2s1;
      d2s1	=ESDPc1[i];
      ds1	=d2s1*ESpa[i];
      d2s1	=2.*d2s1*ESsa[i]+ESDPc2[i]*ESpa[i];
      S		+=h2sh*(ds0+ds1)+h12p2sh*(d2s0-d2s1);

      dS0	=dS1;
      d2S0	=d2S1;
      d2S1	=ESsp1a[i]*ESDPc[i]+ESsp[i]*ESDPc1[i];

      dS1	=d2S1*ESpa[i];
      d2S1	=2.*d2S1*ESsa[i]
	+(ESsp2a[i]*ESDPc[i]+2.*ESsp1a[i]*ESDPc1[i]+ESsp[i]*ESDPc2[i])*ESpa[i];
      aps	+=h2sh*(dS0+dS1)+h12p2sh*(d2S0-d2S1);

      ESgb[i]	=-2.*(aps-ESsp[i]*S)/(ESaJ[i]*ESaJ[i]);
      dgF0	=dgF1;
      d2gF0	=d2gF1;
      dgb0	=dgb1;
      d2gb0	=d2gb1;
      d2gF1	=ESaF[i]*ESLc[i]*rR0;
      dgF1	=ESsa[i]*d2gF1;
      d2gF1	+=ESsa[i]*(ESaF1a[i]*ESLc[i]+ESaF[i]*ESLc1[i])*rR0;
      ESgF1a[i]	=dgF1;
      ESgF[i]	=ESgF[i1]+h2sh*(dgF0+dgF1)+h12p2sh*(d2gF0-d2gF1);
      dgb1	=ESsp[i]*dV1;
      d2gb1	=ESsp1a[i]*dV1+ESsp[i]*d2V1;
      gBext	+=h2sh*(dgb0+dgb1)+h12p2sh*(d2gb0-d2gb1);

      ESgB[i]	=R2R0*ESsp[i]/ESFF[i];
      ESgB1a[i]	=(R2R0*ESsp1a[i]-ESgB[i]*ESFF1a[i])/ESFF[i];
      i1++;
    }
    dS0	=1./ESaJ[ESNa];
    dS1	-=ESsp[ESNa]*ds1;
    ESgb1a[ESNa]=-(2.*dS1*dS0+2.*ESgb[ESNa]*ESaJ1a[ESNa])*dS0;
  }
  splA(ESqV,ESqV2a,ESqV1a,ESqV1a+ESNa);
  splAA(ESgb,ESgb1a,ESgb2a,ESgb1a,ESgb1a+ESNa);
  splA(ESgF,ESgF2a,ESgF1a,ESgF1a+ESNa);
  splA(ESgB,ESgB2a,ESgB1a,ESgB1a+ESNa);
  ESgFPlV	=ESgF[ESNa];
  if(ESiAx != 0){
    rB	=ESnBL/(double)ESNa;
    rB	=rB*(rB+2.)*ESgFPlV/ESgFPlVx;
    EZgd_vb	=EZd_vb*sqrt(rB)-EZd_vb;
  }

  ECVol		=ESqV[ESNa];

  rB		=gBext/ECVol;
  gBext		*=2.*ESRext*ESRext/(ECVol*ESFF[ESNa]);
  i	=ESNp1*ESNa;
  gBN		=100.*gBext*ESBt*EZcr2*(ESsr[i]-ESsr[i+ESNp/2])/(5.*ESaJ[ESNa]);
#ifdef XWIN__
  rV	=4.*EZcgp*EZcgp*ESaR0[0];
  printf("<gb>_ext=%10.3e gb_N=%10.3e V=%10.3e <p>=%10.3e\n"
	 ,gBext,gBN,ECVol*rV,rB);
#endif
  rB		=1./ECb[ESNa];
  rV		=1./ECVol;
  ESqV1a[0]	=sqrt(EZcr2*ESLc[0]*rV);
  ESqV1a[ESNa]	*=EZcr2*rV;
  rgF		=1./ESgF[ESNa];
  ESqgF1a[0]	=sqrt(EZcr2*ESqgF1a[0]*rgF);
  ESqgF1a[ESNa]	=EZcr2*ESgF1a[ESNa]*rgF;
  rgY		=1./ESgY[ESNa];
  ESqgY1a[0]	=EZcr2*sqrt(-ESjs[0]*ESLc[0]*rgY/ESg22c[0]);
  ESqgY1a[ESNa]	=EZcr2*ESgY1a[ESNa]*rgY;
  i		=0;
  ESsh[i]	=0.;
  ESsh1a[i]	=ECb1a[i]*rB;
  ESsh2a[i]	=ECb2a[i]*rB;
  ESqV[i]	=0.;
  ESqgF[i]	=0.;
  ESqgY[i]	=0.;
  for(i=1; i < ESNa1; i++){
    ESsh[i]	=ECb[i]*rB;
    ESsh1a[i]	=ECb1a[i]*rB;
    ESsh2a[i]	=ECb2a[i]*rB;
    ESqV[i]	=sqrt(ESqV[i]*rV);
    ESqV1a[i]	*=EZcr2*rV/ESqV[i];
    ESqgF[i]	=sqrt(ESgF[i]*rgF);
    ESqgF1a[i]	=EZcr2*ESgF1a[i]*rgF/ESqgF[i];
    ESqgY[i]	=sqrt(ESgY[i]*rgY);
    ESqgY1a[i]	=EZcr2*ESgY1a[i]*rgY/ESqgY[i];
  }
  splA(ESqV,ESqV2a,ESqV1a,ESqV1a+ESNa);
  splA(ESqgF,ESqgF2a,ESqgF1a,ESqgF1a+ESNa);
  splA(ESqgY,ESqgY2a,ESqgY1a,ESqgY1a+ESNa);
  i		=0;
  gxB[0]	=0.;
  gxV[0]	=0.;
  gxF[0]	=0.;
  gxY[0]	=0.;
  for(i=1; i < ESNa1; i++){
    gxB[i]	=(ESsa[i]-ESsh[i])/ESsh1a[i];
    gxV[i]	=(ESsa[i]-ESqV[i])/ESqV1a[i];
    gxF[i]	=(ESsa[i]-ESqgF[i])/ESqgF1a[i];
    gxY[i]	=(ESsa[i]-ESqgY[i])/ESqgY1a[i];
  }
  ESgm[1]	=ESgm[0]+(ESgm[2]-ESgm[0])*ESpa[1]/ESpa[2];
  for(i=1; i < ESNa1; i++){
    ESsq[i]	=1./ESgm[i];
  }
  splAA(ESgm,ESgm1a,ESgm2a,ESgm1a,ESgm1a+ESNa);
  splAA(ESsq,ESsq1a,ESsq2a,ESsq1a,ESsq1a+ESNa);
#ifdef MSE
  ESXProf();
#endif
  return(0);
}

int ESEqSolv1DjsjpHole()
{
  int i;

  rR0	=1./R0;
  switch(ESEqSolvInPr){
  case 1:
    for(i=0; i < ESNa1; i++){
      ESjp[i]	=R0*ESPs[i];
      ESjp1a[i]	=R0*ESPs1a[i];
      ESjp2a[i]	=R0*ESPs2a[i];
    }
    break;
  default :
    for(i=0; i < ESNa1; i++){
      ESPs[i]	=rR0*ESjp[i];
      ESPs1a[i]	=rR0*ESjp1a[i];
      ESPs2a[i]	=rR0*ESjp2a[i];
    }
    break;
  }
  switch(ESEqSolvInCr){
  case 3:
    for(i=0; i < ESNa1; i++){
      ESjs[i]	=ESjp[i]+rR0*ESaT[i];
      ESjs1a[i]	=ESjp1a[i]+rR0*ESaT1a[i];
      ESjs2a[i]	=ESjp2a[i]+rR0*ESaT2a[i];
    }
    break;
  default :
    for(i=0; i < ESNa1; i++){
      ESaT[i]	=R0*(ESjs[i]-ESjp[i]);
      ESaT1a[i]	=R0*(ESjs1a[i]-ESjp1a[i]);
      ESaT2a[i]	=R0*(ESjs2a[i]-ESjp2a[i]);
    }
    break;
  }
  {
    int i1;
    double h2sh,h12p2sh;
    double dJ0,dJ1,d2J0,d2J1;
    double dgY0,dgY1,d2gY0,d2gY1;
    double dsp0,dsp1,d2sp0,d2sp1;
    double dFF0,dFF1,d2FF0,d2FF1;

    h12p2sh	=ESsa[1]-ESsa[0];
    h2sh	=EZcr2*h12p2sh;
    h12p2sh	*=EZcr12*h12p2sh;

    i		=0;
    d2J1	=ESjs[i]*ESLc[i]+ESjp[i]*ESVc[i];
    dJ1		=ESsa[i]*d2J1;
    d2J1	+=ESsa[i]*(ESjs1a[i]*ESLc[i]+ESjp1a[i]*ESVc[i]
			   +ESjs[i]*ESLc1[i]+ESjp[i]*ESVc1[i]);
    ESaJ[i]	=0.;
    ESaJ1a[i]	=dJ1;

    ESgY[0]	=0.;
    ESgY1a[0]	=0.;
    dgY1	=0.;
    d2gY1	=-dJ1/(ESsa[0]*ESg22c[0]);

    ESsp[0]	=0.;
    ESsp1a[0]	=0.;
    dsp1	=0.;
    d2sp1	=ESPs[0]*d2gY1;

    ESFF[0]	=0.;
    ESFF1a[0]	=0.;
    dFF1	=0.;
    d2FF1	=2.*ESaT[0]*d2gY1;

    ESgm[i]	=0.;
    ESgm1a[i]	=-R0*d2gY1/(ESsa[i]*ESLc[i]);
    i1	=0;
    for(i=1; i < ESNa1; i++){
      dJ0	=dJ1;
      d2J0	=d2J1;
      d2J1	=ESjs[i]*ESLc[i]+ESjp[i]*ESVc[i];
      dJ1	=ESsa[i]*d2J1;
      d2J1	+=ESsa[i]*(ESjs1a[i]*ESLc[i]+ESjp1a[i]*ESVc[i]
			   +ESjs[i]*ESLc1[i]+ESjp[i]*ESVc1[i]);
      ESaJ1a[i]	=dJ1;
      ESaJ[i]	=ESaJ[i1]+h2sh*(dJ0+dJ1)+h12p2sh*(d2J0-d2J1);

      dgY0	=dgY1;
      d2gY0	=d2gY1;
      d2gY1	=1./(ESsa[i]*ESg22c[i]);
      dgY1	=-ESaJ[i]*d2gY1;
      d2gY1	*=ESaJ[i]*(1./ESsa[i]+ESg22c1[i]/ESg22c[i])-dJ1;
      ESgY1a[i]	=dgY1;
      ESgY[i]	=ESgY[i1]+h2sh*(dgY0+dgY1)+h12p2sh*(d2gY0-d2gY1);

      dsp0	=dsp1;
      d2sp0	=d2sp1;
      dsp1	=ESPs[i]*dgY1;
      d2sp1	=ESPs1a[i]*dgY1+ESPs[i]*d2gY1;
      ESsp1a[i]	=dsp1;
      ESsp[i]	=ESsp[i1]+h2sh*(dsp0+dsp1)+h12p2sh*(d2sp0-d2sp1);

      dFF0	=dFF1;
      d2FF0	=d2FF1;
      d2FF1	=2.*ESaT[i];
      dFF1	=d2FF1*dgY1;
      d2FF1	=d2FF1*d2gY1+2.*ESaT1a[i]*dgY1;
      ESFF1a[i]	=dFF1;
      ESFF[i]	=ESFF[i1]+h2sh*(dFF0+dFF1)+h12p2sh*(d2FF0-d2FF1);

      ESgm[i]	=-R0*dgY1/(ESsa[i]*ESLc[i]);
      ESgm1a[i]	=R0*(dgY1*(1./ESsa[i]+ESLc1[i]/ESLc[i])-d2gY1)
	/(ESsa[i]*ESLc[i]);
      i1++;
    }
    dsp1	=-ESsp[ESNa];
    dFF1	=ESRBt*ESRBt-ESFF[ESNa];
    for(i=0; i < ESNa1; i++){
      ESsp[i]	+=dsp1;
      ESFF[i]	+=dFF1;
      ESaF[i]	=sqrt(ESFF[i]);
      dFF0	=1./ESaF[i];
      ESaF1a[i]	=EZcr2*ESFF1a[i]*dFF0;
      ESgm[i]	*=dFF0;
      ESgm1a[i]	=(ESgm1a[i]-ESgm[i]*ESaF1a[i])*dFF0;
    }
  }
  splA(ESaJ,ESaJ2a,ESaJ1a,ESaJ1a+ESNa);
  splA(ESgY,ESgY2a,ESgY1a,ESgY1a+ESNa);
  splA(ESsp,ESsp2a,ESsp1a,ESsp1a+ESNa);
  splA(ESFF,ESFF2a,ESFF1a,ESFF1a+ESNa);
  splA(ESaF,ESaF2a,ESaF1a,ESaF1a+ESNa);
  splAA(ESgm,ESgm1a,ESgm2a,ESgm1a,ESgm1a+ESNa);

  for(i=0; i < ESNa1; i++){
    ESsq[i]	=ESgm[i] > 1e-4 ? 1./ESgm[i] : 1./(ESgm[i]+1e-4);
    ESsq1a[i]	=-ESgm1a[i]*ESsq[i]*ESsq[i];

    ESjb2a[i]	=ESg22c[i]*ESaT[i]*ESgY1a[i]*ESgY1a[i]/ESFF[i];
    ESjb1a[i]	=(ESVc[i]*ESjp[i]+ESjb2a[i])/ESLc[i];
    ESjb[i]	=ESjs[i]+ESjb1a[i];

    ESjb1a[i]	=ESjs1a[i]+
      (ESVc[i]*ESjp1a[i]+ESVc1[i]*ESjp[i]
       +((ESaT1a[i]*ESg22c[i]*ESgY1a[i]
	  +ESaT[i]*(ESg22c1[i]*ESgY1a[i]+2.*ESg22c[i]*ESgY2a[i]))*ESgY1a[i]
	 -ESjb2a[i]*ESFF1a[i]
	 )/ESFF[i]
       -ESjb1a[i]*ESLc1[i]
       )/ESLc[i];

    ESjB1a[i]	=ESLc[i]*ESaF[i]*rR0/(ESLc[i]+ESVc[i]);
    ESjB[i]	=ESjb[i]*ESjB1a[i];
    ESjB1a[i]	=ESjb1a[i]*ESjB1a[i]
      +ESjb[i]*(ESLc1[i]*ESaF[i]
		+ESLc[i]*(ESaF1a[i]
			  -ESaF[i]*(ESLc1[i]+ESVc1[i])/(ESLc[i]+ESVc[i]))
		)*rR0/(ESLc[i]+ESVc[i]);
  }
  splA(ESsq,ESsq2a,ESsq1a,ESsq1a+ESNa);
  splA(ESjb,ESjb2a,ESjb1a,ESjb1a+ESNa);
  splA(ESjB,ESjB2a,ESjB1a,ESjB1a+ESNa);
  for(i=0; i < ESNa1; i++){
    ESjR[i]	=ESjb[i]*R0;
    ESjR1a[i]	=ESjb1a[i]*R0;
    ESjR2a[i]	=ESjb2a[i]*R0;
  }
  EC1DRadCSolvHole();
  return(0);
}

int ESEqSolv1Djsjp()
{
  int i;

  if(ESa0 != 0.){
    ESEqSolv1DjsjpHole();
    return(0);
  }

  rR0	=1./R0;
  switch(ESEqSolvInPr){
  case 1:
    for(i=0; i < ESNa1; i++){
      ESjp[i]	=R0*ESPs[i];
      ESjp1a[i]	=R0*ESPs1a[i];
      ESjp2a[i]	=R0*ESPs2a[i];
    }
    break;
  default :
    for(i=0; i < ESNa1; i++){
      ESPs[i]	=rR0*ESjp[i];
      ESPs1a[i]	=rR0*ESjp1a[i];
      ESPs2a[i]	=rR0*ESjp2a[i];
    }
    break;
  }
  switch(ESEqSolvInCr){
  case 3:
    for(i=0; i < ESNa1; i++){
      ESjs[i]	=ESjp[i]+rR0*ESaT[i];
      ESjs1a[i]	=ESjp1a[i]+rR0*ESaT1a[i];
      ESjs2a[i]	=ESjp2a[i]+rR0*ESaT2a[i];
    }
    break;
  default :
    for(i=0; i < ESNa1; i++){
      ESaT[i]	=R0*(ESjs[i]-ESjp[i]);
      ESaT1a[i]	=R0*(ESjs1a[i]-ESjp1a[i]);
      ESaT2a[i]	=R0*(ESjs2a[i]-ESjp2a[i]);
    }
    break;
  }
  {
    int i1;
    double h2sh,h12p2sh;
    double dJ0,dJ1,d2J0,d2J1;
    double dgY0,dgY1,d2gY0,d2gY1;
    double dsp0,dsp1,d2sp0,d2sp1;
    double dFF0,dFF1,d2FF0,d2FF1;

    h2sh	=EZcr2*ESsa[1];
    h12p2sh	=EZcr12*ESpa[1];

    ESaJ[0]	=0.;
    ESaJ1a[0]	=0.;
    dJ1		=0.;
    d2J1	=ESjs[0]*ESLc[0];
    ESgY[0]	=0.;
    ESgY1a[0]	=0.;
    dgY1	=0.;
    d2gY1	=-EZcr2*d2J1/ESg22c[0];
    ESsp[0]	=0.;
    ESsp1a[0]	=0.;
    dsp1	=0.;
    d2sp1	=ESPs[0]*d2gY1;
    ESFF[0]	=0.;
    ESFF1a[0]	=0.;
    dFF1	=0.;
    d2FF1	=2.*ESaT[0]*d2gY1;
    ESsq[0]	=2.*rR0*ESg22c[0]/ESjs[0];
    ESsq1a[0]	=0.;
    i1	=0;
    for(i=1; i < ESNa1; i++){
      dJ0	=dJ1;
      d2J0	=d2J1;
      dgY0	=dgY1;
      d2gY0	=d2gY1;
      dsp0	=dsp1;
      d2sp0	=d2sp1;
      dFF0	=dFF1;
      d2FF0	=d2FF1;
      d2J1	=ESjs[i]*ESLc[i]+ESjp[i]*ESVc[i];
      dJ1	=ESsa[i]*d2J1;
      d2J1	+=ESsa[i]*(ESjs1a[i]*ESLc[i]+ESjp1a[i]*ESVc[i]
			+ESjs[i]*ESLc1[i]+ESjp[i]*ESVc1[i]);
      ESaJ1a[i]	=dJ1;
      ESaJ[i]	=ESaJ[i1]+h2sh*(dJ0+dJ1)+h12p2sh*(d2J0-d2J1);

      d2gY1	=ESpa[i]*ESg22c[i]*ESLc[i]/ESaJ[i];
      ESsq[i]	=rR0*d2gY1;
      ESsq1a[i]	=rR0*(2.*ESsa[i]*ESg22c[i]*ESLc[i]-d2gY1*dJ1
		      +ESpa[i]*(ESg22c1[i]*ESLc[i]+ESg22c[i]*ESLc1[i]))/ESaJ[i];
      d2gY1	=1./(ESsa[i]*ESg22c[i]);
      dgY1	=-ESaJ[i]*d2gY1;
      d2gY1	=(ESaJ[i]*(1./ESsa[i]+ESg22c1[i]/ESg22c[i])-dJ1)*d2gY1;
      ESgY1a[i]	=dgY1;
      ESgY[i]	=ESgY[i1]+h2sh*(dgY0+dgY1)+h12p2sh*(d2gY0-d2gY1);
      dsp1	=ESPs[i]*dgY1;
      d2sp1	=ESPs1a[i]*dgY1+ESPs[i]*d2gY1;
      ESsp1a[i]	=dsp1;
      ESsp[i]	=ESsp[i1]+h2sh*(dsp0+dsp1)+h12p2sh*(d2sp0-d2sp1);
      d2FF1	=2.*ESaT[i];
      dFF1	=d2FF1*dgY1;
      d2FF1	=d2FF1*d2gY1+2.*ESaT1a[i]*dgY1;
      ESFF1a[i]	=dFF1;
      ESFF[i]	=ESFF[i1]+h2sh*(dFF0+dFF1)+h12p2sh*(d2FF0-d2FF1);

      ESgm[i]	=-R0*dgY1/(ESsa[i]*ESLc[i]);
      ESgm1a[i]	=R0*(dgY1*(1./ESsa[i]+ESLc1[i]/ESLc[i])-d2gY1)
	/(ESsa[i]*ESLc[i]);

      i1++;
    }
  }
  splA(ESaJ,ESaJ2a,ESaJ1a,ESaJ1a+ESNa);
  splA(ESsp,ESsp2a,ESsp1a,ESsp1a+ESNa);
  splA(ESFF,ESFF2a,ESFF1a,ESFF1a+ESNa);
  splA(ESgY,ESgY2a,ESgY1a,ESgY1a+ESNa);
  i		=0;
  ESsp[0]	=-ESsp[ESNa];
  ESFF[0]	=ESRBt*ESRBt-ESFF[ESNa];
  if(ESFF[0] < 0.){
    ESaF[0]	=0.;
    ESgm[0]	=1./ESsq[0];
    ESsq[0]	=0.;
  }
  else{
    ESaF[0]	=sqrt(ESFF[0]);
    ESsq[0]	*=ESaF[0];
    ESgm[0]	=1./ESsq[0];
  }
  ESaF1a[0]	=0.;
  ESgm1a[0]	=0.;
  ESjb[0]	=ESjs[0];
  ESjb1a[0]	=0.;
  ESjB[0]	=ESjb[0]*ESaF[0]*rR0;
  ESjB1a[0]	=0.;
  for(i=1; i < ESNa1; i++){
    ESsp[i]	+=ESsp[0];
    ESFF[i]	+=ESFF[0];
    if(ESFF[i] < 0.){
      ESaF[i]	=0.;
      ESaF1a[i]	=0.;
      ESgm[i]	=ESgm[i-1];
      ESgm1a[i]	=0.;
      ESsq[i]	=0.;
      ESsq1a[i]	=0.;
    }
    else{
      ESaF[i]	=sqrt(ESFF[i]);
      ESaF1a[i]	=EZcr2*ESFF1a[i]/ESaF[i];
      ESsq1a[i]	=ESsq1a[i]*ESaF[i]+ESsq[i]*ESaF1a[i];
      ESsq[i]	*=ESaF[i];
      ESgm[i]	=1./ESsq[i];
      ESgm1a[i]	=-ESsq1a[i]*ESgm[i]*ESgm[i];
    }

    ESjb2a[i]	=ESg22c[i]*ESaT[i]*ESgY1a[i]*ESgY1a[i]/ESFF[i];
    ESjb1a[i]	=(ESVc[i]*ESjp[i]+ESjb2a[i])/ESLc[i];
    ESjb[i]	=ESjs[i]+ESjb1a[i];
    ESjb1a[i]	=ESjs1a[i]+
      (ESVc[i]*ESjp1a[i]+ESVc1[i]*ESjp[i]
       +((ESaT1a[i]*ESg22c[i]*ESgY1a[i]
	  +ESaT[i]*(ESg22c1[i]*ESgY1a[i]+2.*ESg22c[i]*ESgY2a[i]))*ESgY1a[i]
	 -ESjb2a[i]*ESFF1a[i]
	 )/ESFF[i]
       -ESjb1a[i]*ESLc1[i]
       )/ESLc[i];
    ESjB1a[i]	=ESLc[i]*ESaF[i]*rR0/(ESLc[i]+ESVc[i]);
    ESjB[i]	=ESjb[i]*ESjB1a[i];
    ESjB1a[i]	=ESjb1a[i]*ESjB1a[i]
      +ESjb[i]*(ESLc1[i]*ESaF[i]
		+ESLc[i]*(ESaF1a[i]
			  -ESaF[i]*(ESLc1[i]+ESVc1[i])/(ESLc[i]+ESVc[i]))
		)*rR0/(ESLc[i]+ESVc[i]);
  }
  splA(ESaF,ESaF2a,ESaF1a,ESaF1a+ESNa);
  splA(ESsq,ESsq2a,ESsq1a,ESsq1a+ESNa);
  splA(ESgm,ESgm2a,ESgm1a,ESgm1a+ESNa);
  splA(ESjb,ESjb2a,ESjb1a,ESjb1a+ESNa);
  splA(ESjB,ESjB2a,ESjB1a,ESjB1a+ESNa);
  for(i=0; i < ESNa1; i++){
    ESjR[i]	=ESjb[i]*R0;
    ESjR1a[i]	=ESjb1a[i]*R0;
    ESjR2a[i]	=ESjb2a[i]*R0;
  }
  EC1DRadCSolv();
  return(0);
}

int ESEqSolv1DaFjp()
{
  int i,i1;
  double ainv;
  double h2sh,h12p2sh;
  double dJ0,dJ1,d2J0,d2J1;
  double dgY0,dgY1,d2gY0,d2gY1;
  double dsp0,dsp1,d2sp0,d2sp1;
  double dFF0,dFF1,d2FF0,d2FF1;

  rR0	=1./R0;
  switch(ESEqSolvInPr){
  case 1:
    for(i=0; i < ESNa1; i++){
      ESjp[i]	=R0*ESPs[i];
      ESjp1a[i]	=R0*ESPs1a[i];
      ESjp2a[i]	=R0*ESPs2a[i];
    }
    break;
  default :
    for(i=0; i < ESNa1; i++){
      ESPs[i]	=rR0*ESjp[i];
      ESPs1a[i]	=rR0*ESjp1a[i];
      ESPs2a[i]	=rR0*ESjp2a[i];
    }
    break;
  }
  ainv	=0.;
  for(i=0; i < ESNa1; i++){
    ESjs[i]	=ESjp[i]+rR0*ESaT[i];
    ESjs1a[i]	=ESjp1a[i]+rR0*ESaT1a[i];
    ESjs2a[i]	=ESjp2a[i]+rR0*ESaT2a[i];
    if(i && ESaT[i-1] > 0. && ESaT[i] <= 0.){
      ainv	=ESsa[i]-(ESsa[i]-ESsa[i-1])*ESaT[i]/(ESaT[i]-ESaT[i-1]);
    }
  }

  if(ainv > 0.){
    do{
      ESSetSplA(ainv);
      splRA(&dJ0,&dJ1,ESaT,ESaT2a);
      dJ0	/=dJ1;
      ainv	-=dJ0;
    }while(fabs(dJ0) > 1e-6);
  }    

  h2sh	=EZcr2*ESsa[1];
  h12p2sh	=EZcr12*ESpa[1];
  
  ESaJ[0]	=0.;
  ESaJ1a[0]	=0.;
  dJ1		=0.;
  d2J1	=ESjs[0]*ESLc[0];
  ESgY[0]	=0.;
  ESgY1a[0]	=0.;
  dgY1	=0.;
  d2gY1	=-EZcr2*d2J1/ESg22c[0];
  ESsp[0]	=0.;
  ESsp1a[0]	=0.;
  dsp1	=0.;
  d2sp1	=ESPs[0]*d2gY1;
  ESFF[0]	=0.;
  ESFF1a[0]	=0.;
  dFF1	=0.;
  d2FF1	=2.*ESaT[0]*d2gY1;
  ESsq[0]	=2.*rR0*ESg22c[0]/ESjs[0];
  ESsq1a[0]	=0.;
  i1	=0;
  for(i=1; i < ESNa1; i++){
    dJ0	=dJ1;
    d2J0	=d2J1;
    dgY0	=dgY1;
    d2gY0	=d2gY1;
    dsp0	=dsp1;
    d2sp0	=d2sp1;
    dFF0	=dFF1;
    d2FF0	=d2FF1;
    d2J1	=ESjs[i]*ESLc[i]+ESjp[i]*ESVc[i];
    dJ1	=ESsa[i]*d2J1;
    d2J1	+=ESsa[i]*(ESjs1a[i]*ESLc[i]+ESjp1a[i]*ESVc[i]
			   +ESjs[i]*ESLc1[i]+ESjp[i]*ESVc1[i]);
    ESaJ1a[i]	=dJ1;
    ESaJ[i]	=ESaJ[i1]+h2sh*(dJ0+dJ1)+h12p2sh*(d2J0-d2J1);

    d2gY1	=ESpa[i]*ESg22c[i]*ESLc[i]/ESaJ[i];
    ESsq[i]	=rR0*d2gY1;
    ESsq1a[i]	=rR0*(2.*ESsa[i]*ESg22c[i]*ESLc[i]-d2gY1*dJ1
		      +ESpa[i]*(ESg22c1[i]*ESLc[i]+ESg22c[i]*ESLc1[i]))/ESaJ[i];
    d2gY1	=1./(ESsa[i]*ESg22c[i]);
    dgY1	=-ESaJ[i]*d2gY1;
    d2gY1	=(ESaJ[i]*(1./ESsa[i]+ESg22c1[i]/ESg22c[i])-dJ1)*d2gY1;
    ESgY1a[i]	=dgY1;
    ESgY[i]	=ESgY[i1]+h2sh*(dgY0+dgY1)+h12p2sh*(d2gY0-d2gY1);
    dsp1	=ESPs[i]*dgY1;
    d2sp1	=ESPs1a[i]*dgY1+ESPs[i]*d2gY1;
    ESsp1a[i]	=dsp1;
    ESsp[i]	=ESsp[i1]+h2sh*(dsp0+dsp1)+h12p2sh*(d2sp0-d2sp1);
    d2FF1	=2.*ESaT[i];
    dFF1	=d2FF1*dgY1;
    d2FF1	=d2FF1*d2gY1+2.*ESaT1a[i]*dgY1;
    ESFF1a[i]	=dFF1;
    ESFF[i]	=ESFF[i1]+h2sh*(dFF0+dFF1)+h12p2sh*(d2FF0-d2FF1);

    ESgm[i]	=-R0*dgY1/(ESsa[i]*ESLc[i]);
    ESgm1a[i]	=R0*(dgY1*(1./ESsa[i]+ESLc1[i]/ESLc[i])-d2gY1)
      /(ESsa[i]*ESLc[i]);
    i1++;
  }
  splA(ESaJ,ESaJ2a,ESaJ1a,ESaJ1a+ESNa);
  splA(ESsp,ESsp2a,ESsp1a,ESsp1a+ESNa);
  splA(ESFF,ESFF2a,ESFF1a,ESFF1a+ESNa);
  splA(ESgY,ESgY2a,ESgY1a,ESgY1a+ESNa);
  i		=0;
  ESsp[0]	=-ESsp[ESNa];

  ESSetSplA(ainv);
  splRA(&dJ0,NULL,ESFF,ESFF2a);
  ESFF[0]	-=dJ0;
  ESaF[0]	=sqrt(ESFF[0]);
  ESsq[0]	*=ESaF[0];
  ESgm[0]	=1./ESsq[0];

  ESaF1a[0]	=0.;
  ESgm1a[0]	=0.;
  ESjb[0]	=ESjs[0];
  ESjb1a[0]	=0.;
  ESjB[0]	=ESjb[0]*ESaF[0]*rR0;
  ESjB1a[0]	=0.;
  for(i=1; i < ESNa1; i++){
    ESsp[i]	+=ESsp[0];
    ESFF[i]	-=dJ0;
    ESaF[i]	=sqrt(ESFF[i]);
    ESaF1a[i]	=EZcr2*ESFF1a[i]/ESaF[i];
    ESsq1a[i]	=ESsq1a[i]*ESaF[i]+ESsq[i]*ESaF1a[i];
    ESsq[i]	*=ESaF[i];
    if(ESsa[i] > ainv){
      ESaF[i]	=-ESaF[i];
      ESaF1a[i]	=-ESaF1a[i];
      ESsq1a[i]	=-ESsq1a[i];
      ESsq[i]	=-ESsq[i];
    }
    if(ESsq[i] != 0.){
      ESgm[i]	=1./ESsq[i];
      ESgm1a[i]	=-ESsq1a[i]*ESgm[i]*ESgm[i];
    }
    else{
      ESgm[i]	=0.;
      ESgm1a[i]	=0.;
    }
    ESjb2a[i]	=ESg22c[i]*ESaT[i]*ESgY1a[i]*ESgY1a[i]/ESFF[i];
    ESjb1a[i]	=(ESVc[i]*ESjp[i]+ESjb2a[i])/ESLc[i];
    ESjb[i]	=ESjs[i]+ESjb1a[i];
    ESjb1a[i]	=ESjs1a[i]+
      (ESVc[i]*ESjp1a[i]+ESVc1[i]*ESjp[i]
       +((ESaT1a[i]*ESg22c[i]*ESgY1a[i]
	  +ESaT[i]*(ESg22c1[i]*ESgY1a[i]+2.*ESg22c[i]*ESgY2a[i]))*ESgY1a[i]
	 -ESjb2a[i]*ESFF1a[i]
	 )/ESFF[i]
       -ESjb1a[i]*ESLc1[i]
       )/ESLc[i];
    ESjB1a[i]	=ESLc[i]*ESaF[i]*rR0/(ESLc[i]+ESVc[i]);
    ESjB[i]	=ESjb[i]*ESjB1a[i];
    ESjB1a[i]	=ESjb1a[i]*ESjB1a[i]
      +ESjb[i]*(ESLc1[i]*ESaF[i]
		+ESLc[i]*(ESaF1a[i]
			  -ESaF[i]*(ESLc1[i]+ESVc1[i])/(ESLc[i]+ESVc[i]))
		)*rR0/(ESLc[i]+ESVc[i]);
  }
  ESRBt		=ESaF[ESNa];
  ESBt		=ESRBt/ESRext;
  splA(ESaF,ESaF2a,ESaF1a,ESaF1a+ESNa);
  splA(ESsq,ESsq2a,ESsq1a,ESsq1a+ESNa);
  splA(ESgm,ESgm2a,ESgm1a,ESgm1a+ESNa);
  splA(ESjb,ESjb2a,ESjb1a,ESjb1a+ESNa);
  splA(ESjB,ESjB2a,ESjB1a,ESjB1a+ESNa);
  for(i=0; i < ESNa1; i++){
    ESjR[i]	=ESjb[i]*R0;
    ESjR1a[i]	=ESjb1a[i]*R0;
    ESjR2a[i]	=ESjb2a[i]*R0;
  }
  EC1DRadCSolv();
  return(0);
}

int ESEqSolv1Djssp()
{
  int i;

  rR0	=1./R0;
  {
    int i1;
    double h2sh,h12p2sh;
    double dJ0,dJ1,d2J0,d2J1;
    double dgY0,dgY1,d2gY0,d2gY1;
    double dsp0,dsp1,d2sp0,d2sp1;
    double dFF0,dFF1,d2FF0,d2FF1;

    h2sh	=EZcr2*ESsa[1];
    h12p2sh	=EZcr12*ESpa[1];

    ESaJ[0]	=0.;
    ESaJ1a[0]	=0.;
    dJ1		=0.;
    d2J1	=ESjs[0]*ESLc[0];
    ESjp[0]	=-2.*R0*ESg22c[0]*ESsp2a[0]/d2J1;
    ESjp1a[0]	=0.;
    ESgY[0]	=0.;
    ESgY1a[0]	=0.;
    dgY1	=0.;
    d2gY1	=-EZcr2*d2J1/ESg22c[0];
    ESsp[0]	=0.;
    ESsp1a[0]	=0.;
    dsp1	=0.;
    d2sp1	=ESPs[0]*d2gY1;
    ESFF[0]	=0.;
    ESFF1a[0]	=0.;
    dFF1	=0.;
    d2FF1	=2.*ESaT[0]*d2gY1;
    ESsq[0]	=2.*rR0*ESg22c[0]/ESjs[0];
    ESsq1a[0]	=0.;
    ESgm[0]	=R0*ESjs[0]/ESg22c[0];
    ESgm1a[0]	=0.;
    i1	=0;
    for(i=1; i < ESNa1; i++){
      dJ0	=dJ1;
      d2J0	=d2J1;
      dgY0	=dgY1;
      d2gY0	=d2gY1;
      dsp0	=dsp1;
      d2sp0	=d2sp1;
      dFF0	=dFF1;
      d2FF0	=d2FF1;
      d2J1	=ESjs[i]*ESLc[i]+ESjp[i]*ESVc[i];
      dJ1	=ESsa[i]*d2J1;
      d2J1	+=ESsa[i]*(ESjs1a[i]*ESLc[i]+ESjp1a[i]*ESVc[i]
			+ESjs[i]*ESLc1[i]+ESjp[i]*ESVc1[i]);
      ESaJ1a[i]	=dJ1;
      ESaJ[i]	=ESaJ[i1]+h2sh*(dJ0+dJ1)+h12p2sh*(d2J0-d2J1);

      d2gY1	=ESpa[i]*ESg22c[i]*ESLc[i]/ESaJ[i];
      ESsq[i]	=rR0*d2gY1;
      ESsq1a[i]	=rR0*(2.*ESsa[i]*ESg22c[i]*ESLc[i]-d2gY1*dJ1
		      +ESpa[i]*(ESg22c1[i]*ESLc[i]+ESg22c[i]*ESLc1[i]))/ESaJ[i];
      d2gY1	=1./(ESsa[i]*ESg22c[i]);
      dgY1	=-ESaJ[i]*d2gY1;
      d2gY1	=(ESaJ[i]*(1./ESsa[i]+ESg22c1[i]/ESg22c[i])-dJ1)*d2gY1;
      ESgY1a[i]	=dgY1;
      ESgY[i]	=ESgY[i1]+h2sh*(dgY0+dgY1)+h12p2sh*(d2gY0-d2gY1);

      dsp1	=ESPs[i]*dgY1;
      d2sp1	=ESPs1a[i]*dgY1+ESPs[i]*d2gY1;
      ESsp1a[i]	=dsp1;
      ESsp[i]	=ESsp[i1]+h2sh*(dsp0+dsp1)+h12p2sh*(d2sp0-d2sp1);
      d2FF1	=2.*ESaT[i];
      dFF1	=d2FF1*dgY1;
      d2FF1	=d2FF1*d2gY1+2.*ESaT1a[i]*dgY1;
      ESFF1a[i]	=dFF1;
      ESFF[i]	=ESFF[i1]+h2sh*(dFF0+dFF1)+h12p2sh*(d2FF0-d2FF1);
      ESgm[i]	=-R0*dgY1/(ESsa[i]*ESLc[i]);
      ESgm1a[i]	=R0*(dgY1*(1./ESsa[i]+ESLc1[i]/ESLc[i])-d2gY1)
	/(ESsa[i]*ESLc[i]);
      i1++;
    }
  }
  splA(ESaJ,ESaJ2a,ESaJ1a,ESaJ1a+ESNa);
  splA(ESjp,ESjp2a,ESjp1a,ESjp1a+ESNa);
  splA(ESFF,ESFF2a,ESFF1a,ESFF1a+ESNa);
  splA(ESgY,ESgY2a,ESgY1a,ESgY1a+ESNa);
  for(i=0; i < ESNa1; i++){
    ESaT[i]	=R0*(ESjs[i]-ESjp[i]);
    ESaT1a[i]	=R0*(ESjs1a[i]-ESjp1a[i]);
    ESaT2a[i]	=R0*(ESjs2a[i]-ESjp2a[i]);
  }
  i		=0;
  ESsp[0]	=-ESsp[ESNa];
  ESFF[0]	=ESRBt*ESRBt-ESFF[ESNa];
  if(ESFF[0] < 0.){
    ESaF[0]	=0.;
    ESgm[0]	=1./ESsq[0];
    ESsq[0]	=0.;
  }
  else{
    ESaF[0]	=sqrt(ESFF[0]);
    ESsq[0]	*=ESaF[0];
    ESgm[0]	=1./ESsq[0];
  }
  ESaF1a[0]	=0.;
  ESgm1a[0]	=0.;
  ESjb[0]	=ESjs[0];
  ESjb1a[0]	=0.;
  ESjB[0]	=ESjb[0]*ESaF[0]*rR0;
  ESjB1a[0]	=0.;
  for(i=1; i < ESNa1; i++){
    ESsp[i]	+=ESsp[0];
    ESFF[i]	+=ESFF[0];
    if(ESFF[i] < 0.){
      ESaF[i]	=0.;
      ESaF1a[i]	=0.;
      ESgm[i]	=ESgm[i-1];
      ESgm1a[i]	=0.;
      ESsq[i]	=0.;
      ESsq1a[i]	=0.;
    }
    else{
      ESaF[i]	=sqrt(ESFF[i]);
      ESaF1a[i]	=EZcr2*ESFF1a[i]/ESaF[i];
      ESsq1a[i]	=ESsq1a[i]*ESaF[i]+ESsq[i]*ESaF1a[i];
      ESsq[i]	*=ESaF[i];
      ESgm[i]	=1./ESsq[i];
      ESgm1a[i]	=-ESsq1a[i]*ESgm[i]*ESgm[i];
    }

    ESjb2a[i]	=ESg22c[i]*ESaT[i]*ESgY1a[i]*ESgY1a[i]/ESFF[i];
    ESjb1a[i]	=(ESVc[i]*ESjp[i]+ESjb2a[i])/ESLc[i];
    ESjb[i]	=ESjs[i]+ESjb1a[i];
    ESjb1a[i]	=ESjs1a[i]+
      (ESVc[i]*ESjp1a[i]+ESVc1[i]*ESjp[i]
       +((ESaT1a[i]*ESg22c[i]*ESgY1a[i]
	  +ESaT[i]*(ESg22c1[i]*ESgY1a[i]+2.*ESg22c[i]*ESgY2a[i]))*ESgY1a[i]
	 -ESjb2a[i]*ESFF1a[i]
	 )/ESFF[i]
       -ESjb1a[i]*ESLc1[i]
       )/ESLc[i];
    ESjB1a[i]	=ESLc[i]*ESaF[i]*rR0/(ESLc[i]+ESVc[i]);
    ESjB[i]	=ESjb[i]*ESjB1a[i];
    ESjB1a[i]	=ESjb1a[i]*ESjB1a[i]
      +ESjb[i]*(ESLc1[i]*ESaF[i]
		+ESLc[i]*(ESaF1a[i]
			  -ESaF[i]*(ESLc1[i]+ESVc1[i])/(ESLc[i]+ESVc[i]))
		)*rR0/(ESLc[i]+ESVc[i]);
  }
  splA(ESaF,ESaF2a,ESaF1a,ESaF1a+ESNa);
  splA(ESsq,ESsq2a,ESsq1a,ESsq1a+ESNa);
  splA(ESgm,ESgm2a,ESgm1a,ESgm1a+ESNa);
  splA(ESjb,ESjb2a,ESjb1a,ESjb1a+ESNa);
  splA(ESjB,ESjB2a,ESjB1a,ESjB1a+ESNa);
  for(i=0; i < ESNa1; i++){
    ESjR[i]	=ESjb[i]*R0;
    ESjR1a[i]	=ESjb1a[i]*R0;
    ESjR2a[i]	=ESjb2a[i]*R0;
  }
  EC1DRadCSolv();
  return(0);
}

int ESEq1Djbjp()
{
  double s,ds,rs;
  double h2sh,h12p2sh;
  double bJ,dJ0,dJ1,d2J0,d2J1;
  double dF0,dF1,d2F0,d2F1;
  int i,i1;
  
  h2sh		=EZcr2*ESsa[1];
  h12p2sh	=EZcr12*ESpa[1];

  ESjs[0]	=ESjb[0];
  ESjs1a[0]	=0.;
  ESaJ[0]	=0.;
  ESaJ1a[0]	=0.;
  bJ		=0.;
  dJ1		=0.;
  d2J1		=ESLc[0]*ESjb[0]/ESaF[0];

  ESaF[0]	=ESaF[0];
  ESaF1a[0]	=0.;
  dF1		=0.;
  d2F1		=EZcr2*R0*d2J1*(ESjp[0]-ESjb[0])/ESg22c[0];

  i1		=0;
  for(i=1; i < ESNa1; i++){
    dJ0		=dJ1;
    d2J0	=d2J1;
    dF0		=dF1;
    d2F0	=d2F1;
    d2J1	=ESLc[i]*ESjb[i]/ESaF[i];
    dJ1		=ESsa[i]*d2J1;
    d2J1	+=ESsa[i]*(ESLc[i]*ESjb1a[i]+ESLc1[i]*ESjb[i]
			   -d2J1*ESaF1a[i])/ESaF[i];
    bJ		+=h2sh*(dJ0+dJ1)+h12p2sh*(d2J0-d2J1);

    ESaJ[i]	=bJ*ESaF[i];
    ESaJ1a[i]	=dJ1*ESaF[i]+bJ*ESaF1a[i];
    dJ0		=1./ESLc[i];
    d2J0	=ESg22c1[i]/ESg22c[i];
    ds		=R0*bJ*dJ0/(ESpa[i]*ESg22c[i]);
    s		=bJ*ds;
    ds		=2.*dJ1*ds-s*(d2J0+ESLc1[i]*dJ0+2./ESsa[i]);
    rs		=1./(1.+s);
    d2F1	=(1.+ESVc[i]*dJ0)*ESjp[i]-ESjb[i];
    dF1		=R0*bJ*d2F1*rs/(ESsa[i]*ESg22c[i]);
    
    ESjs1a[i]	=(ESVc1[i]-ESVc[i]*ESLc1[i]*dJ0)*dJ0*ESjp[i];
    d2F1	=-dF1*(1./ESsa[i]+d2J0+ds*rs)
      +R0*(dJ1*d2F1
	   +bJ*((1.+ESVc[i]*dJ0)*ESjp1a[i]-ESjb1a[i]+ESjs1a[i])
	   )*rs/(ESsa[i]*ESg22c[i]);
    ESaF[i]	=ESaF[i1]+h2sh*(dF0+dF1)+h12p2sh*(d2F0-d2F1);
    ESaF1a[i]	=dF1;
    
    ESjs[i]	=(ESjb[i]+(s-ESVc[i]*dJ0)*ESjp[i])*rs;
    ESjs1a[i]	=(ESjb1a[i]+(s-ESVc[i]*dJ0)*ESjp1a[i]-ESjs1a[i]
		  +ds*(ESjp[i]-ESjs[i]))*rs;
    i1++;
  }
  ds	=ESRBt-ESaF[ESNa];
  for(i=0; i < ESNa1; i++){
    ESaF[i]	+=ds;
  }
  
  splA(ESjs,ESjs2a,ESjs1a,ESjs1a+ESNa);
  splA(ESaJ,ESaJ2a,ESaJ1a,ESaJ1a+ESNa);
  for(i=0; i < ESNa1; i++){
    ESjB1a[i]	=ESLc[i]*ESaF[i]*rR0/(ESLc[i]+ESVc[i]);
    ESjB[i]	=ESjb[i]*ESjB1a[i];
    ESjB1a[i]	=ESjb1a[i]*ESjB1a[i]
      +ESjb[i]*(ESLc1[i]*ESaF[i]
		+ESLc[i]*(ESaF1a[i]
			  -ESaF[i]*(ESLc1[i]+ESVc1[i])/(ESLc[i]+ESVc[i]))
		)*rR0/(ESLc[i]+ESVc[i]);

    ESFF[i]	=ESaF[i]*ESaF[i];
    ESFF1a[i]	=2.*ESaF[i]*ESaF1a[i];

    ESaT[i]	=R0*(ESjs[i]-ESjp[i]);
    ESaT1a[i]	=R0*(ESjs1a[i]-ESjp1a[i]);
    ESaT2a[i]	=R0*(ESjs2a[i]-ESjp2a[i]);
  }
  splA(ESjB,ESjB2a,ESjB1a,ESjB1a+ESNa);
  splA(ESFF,ESFF2a,ESFF1a,ESFF1a+ESNa);
  splA(ESaF,ESaF2a,ESaF1a,ESaF1a+ESNa);

  return(0);
}

int ESEq1DjBjp()
{
  double s,ds,rs;
  double h2sh,h12p2sh;
  double bJ,dJ0,dJ1,d2J0,d2J1;
  double dF0,dF1,d2F0,d2F1;
  int i,i1;
  
  h2sh		=EZcr2*ESsa[1];
  h12p2sh	=EZcr12*ESpa[1];
  
  ESjb[0]	=R0*(ESLc[0]+ESVc[0])*ESjB[0]/(ESLc[0]*ESaF[0]);
  ESjb1a[0]	=0.;
  ESjs[0]	=ESjb[0];
  ESjs1a[0]	=0.;
  ESaJ[0]	=0.;
  ESaJ1a[0]	=0.;
  bJ		=0.;
  dJ1		=0.;
  d2J1		=ESLc[0]*ESjb[0]/ESaF[0];

  ESaF[0]	=ESaF[0];
  ESaF1a[0]	=0.;
  dF1		=0.;
  d2F1		=EZcr2*R0*d2J1*(ESjp[0]-ESjb[0])/ESg22c[0];

  i1		=0;
  for(i=1; i < ESNa1; i++){
    dJ0		=dJ1;
    d2J0	=d2J1;
    dF0		=dF1;
    d2F0	=d2F1;

    dJ1		=1./(ESLc[i]*ESaF[i]);
    d2J1	=ESLc[i]+ESVc[i];
    ESjb[i]	=R0*d2J1*ESjB[i]*dJ1;
    ESjb1a[i]	=R0*((ESLc1[i]+ESVc1[i])*ESjB[i]+d2J1*ESjB1a[i])*dJ1
      -ESjb[i]*(ESLc1[i]*ESaF[i]+ESLc[i]*ESaF1a[i])*dJ1;

    dF1		=1./ESaF[i];
    d2J1	=ESLc[i]*ESjb[i]*dF1;
    dJ1		=ESsa[i]*d2J1;
    d2J1	+=ESsa[i]*(ESLc[i]*ESjb1a[i]+ESLc1[i]*ESjb[i]
			   -d2J1*ESaF1a[i])*dF1;
    bJ		+=h2sh*(dJ0+dJ1)+h12p2sh*(d2J0-d2J1);

    ESaJ[i]	=bJ*ESaF[i];
    ESaJ1a[i]	=dJ1*ESaF[i]+bJ*ESaF1a[i];
    
    dJ0		=1./ESLc[i];
    d2J0	=ESg22c1[i]/ESg22c[i];
    ds		=R0*bJ*dJ0/(ESpa[i]*ESg22c[i]);
    s		=bJ*ds;
    ds		=2.*dJ1*ds-s*(d2J0+ESLc1[i]*dJ0+2./ESsa[i]);
    rs		=1./(1.+s);
    d2F1	=(1.+ESVc[i]*dJ0)*ESjp[i]-ESjb[i];
    dF1		=R0*bJ*d2F1*rs/(ESsa[i]*ESg22c[i]);
    
    ESjs1a[i]	=(ESVc1[i]-ESVc[i]*ESLc1[i]*dJ0)*dJ0*ESjp[i];
    d2F1	=-dF1*(1./ESsa[i]+d2J0+ds*rs)
      +R0*(dJ1*d2F1
	   +bJ*((1.+ESVc[i]*dJ0)*ESjp1a[i]-ESjb1a[i]+ESjs1a[i])
	   )*rs/(ESsa[i]*ESg22c[i]);
    ESaF[i]	=ESaF[i1]+h2sh*(dF0+dF1)+h12p2sh*(d2F0-d2F1);
    ESaF1a[i]	=dF1;
    
    ESjs[i]	=(ESjb[i]+(s-ESVc[i]*dJ0)*ESjp[i])*rs;
    ESjs1a[i]	=(ESjb1a[i]+(s-ESVc[i]*dJ0)*ESjp1a[i]-ESjs1a[i]
		  +ds*(ESjp[i]-ESjs[i]))*rs;
    i1++;
  }
  ds	=ESRBt-ESaF[ESNa];
  for(i=0; i < ESNa1; i++){
    ESaF[i]	+=ds;
  }
  
  splA(ESjs,ESjs2a,ESjs1a,ESjs1a+ESNa);
  splA(ESaJ,ESaJ2a,ESaJ1a,ESaJ1a+ESNa);
  for(i=0; i < ESNa1; i++){
    ESFF[i]	=ESaF[i]*ESaF[i];
    ESFF1a[i]	=2.*ESaF[i]*ESaF1a[i];

    ESaT[i]	=R0*(ESjs[i]-ESjp[i]);
    ESaT1a[i]	=R0*(ESjs1a[i]-ESjp1a[i]);
    ESaT2a[i]	=R0*(ESjs2a[i]-ESjp2a[i]);
  }
  splA(ESjb,ESjb2a,ESjb1a,ESjb1a+ESNa);
  splA(ESFF,ESFF2a,ESFF1a,ESFF1a+ESNa);
  splA(ESaF,ESaF2a,ESaF1a,ESaF1a+ESNa);

  return(0);
}

int ESEqSolv1Djbjp()
{
  int i;

  RK1Dh	=ESsa[1];
  rR0	=1./R0;
  switch(ESEqSolvInPr){
  case 1:
    for(i=0; i < ESNa1; i++){
      ESjp[i]	=R0*ESPs[i];
      ESjp1a[i]	=R0*ESPs1a[i];
      ESjp2a[i]	=R0*ESPs2a[i];
    }
    break;
  case 2:
    for(i=0; i < ESNa1; i++){
      ESjp[i]	=R0*ESPs[i];
      ESjp1a[i]	=R0*ESPs1a[i];
      ESjp2a[i]	=R0*ESPs2a[i];
    }
    break;
  default :
    for(i=0; i < ESNa1; i++){
      ESPs[i]	=rR0*ESjp[i];
      ESPs1a[i]	=rR0*ESjp1a[i];
      ESPs2a[i]	=rR0*ESjp2a[i];
    }
    break;
  }
  switch(ESEqSolvInCr){
  case 1:
    ESEq1Djbjp();
    for(i=0; i < ESNa1; i++){
      ESjR[i]	=R0*ESjb[i];
      ESjR1a[i]	=R0*ESjb1a[i];
      ESjR2a[i]	=R0*ESjb2a[i];
    }
    break;
  case 2:
    for(i=0; i < ESNa1; i++){
      ESjb[i]	=rR0*ESjR[i];
      ESjb1a[i]	=rR0*ESjR1a[i];
      ESjb2a[i]	=rR0*ESjR2a[i];
    }
    ESEq1Djbjp();
    break;
  case 4:
    ESEq1DjBjp();
    break;
  }
  
  {
    int i1;
    double h2sh,h12p2sh;
    double dgY0,dgY1,d2gY0,d2gY1;
    double dsp0,dsp1,d2sp0,d2sp1;

    h2sh	=EZcr2*ESsa[1];
    h12p2sh	=EZcr12*ESpa[1];

    ESgY[0]	=0.;
    ESgY1a[0]	=0.;
    dgY1	=0.;
    d2gY1	=-EZcr2*ESjs[0]*ESLc[0]/ESg22c[0];
    ESsp[0]	=0.;
    ESsp1a[0]	=0.;
    dsp1	=0.;
    d2sp1	=ESPs[0]*d2gY1;
    ESsq[0]	=2.*rR0*ESg22c[0]*ESaF[0]/ESjs[0];

    ESsq1a[0]	=0.;
    ESgm[0]	=1./ESsq[0];
    ESgm1a[0]	=0.;
    i1	=0;
    for(i=1; i < ESNa1; i++){
      dgY0	=dgY1;
      d2gY0	=d2gY1;
      dsp0	=dsp1;
      d2sp0	=d2sp1;
      d2gY1	=ESpa[i]*ESg22c[i]*ESLc[i]/ESaJ[i];
      ESsq[i]	=rR0*d2gY1;
      ESsq1a[i]	=rR0*(2.*ESsa[i]*ESg22c[i]*ESLc[i]-d2gY1*ESaJ1a[i]
		  +ESpa[i]*(ESg22c1[i]*ESLc[i]+ESg22c[i]*ESLc1[i]))/ESaJ[i];
      ESsq1a[i]	=ESsq1a[i]*ESaF[i]+ESsq[i]*ESaF1a[i];
      ESsq[i]	*=ESaF[i];
      ESgm[i]	=1./ESsq[i];
      ESgm1a[i]	=-ESsq1a[i]*ESgm[i]*ESgm[i];
      
      d2gY1	=1./(ESsa[i]*ESg22c[i]);
      dgY1	=-ESaJ[i]*d2gY1;
      d2gY1	=(ESaJ[i]*(1./ESsa[i]+ESg22c1[i]/ESg22c[i])-ESaJ1a[i])*d2gY1;
      ESgY1a[i]	=dgY1;
      ESgY[i]	=ESgY[i1]+h2sh*(dgY0+dgY1)+h12p2sh*(d2gY0-d2gY1);
      dsp1	=ESPs[i]*dgY1;
      d2sp1	=ESPs1a[i]*dgY1+ESPs[i]*d2gY1;
      ESsp1a[i]	=dsp1;
      ESsp[i]	=ESsp[i1]+h2sh*(dsp0+dsp1)+h12p2sh*(d2sp0-d2sp1);
      i1++;
    }
  }
  splA(ESsp,ESsp2a,ESsp1a,ESsp1a+ESNa);
  splA(ESgY,ESgY2a,ESgY1a,ESgY1a+ESNa);
  splA(ESsq,ESsq2a,ESsq1a,ESsq1a+ESNa);
  splA(ESgm,ESgm2a,ESgm1a,ESgm1a+ESNa);
  ESsp[0]	=-ESsp[ESNa];
  for(i=1; i < ESNa1; i++){
    ESsp[i]	+=ESsp[0];
  }
  EC1DRadCSolv();
  return(0);
}

int ESEq1Djbsp()
{
  double s,ds,rs;
  double h2sh,h12p2sh;
  double bJ,dJ0,dJ1,d2J0,d2J1;
  double dF0,dF1,d2F0,d2F1;
  int i,i1;
  
  h2sh		=EZcr2*ESsa[1];
  h12p2sh	=EZcr12*ESpa[1];
  ESjp[0]	=-2.*R0*ESg22c[0]*ESsp2a[0]/(ESjb[0]*ESLc[0]);
  ESjp1a[0]	=0.;
  ESjs[0]	=ESjb[0];
  ESjs1a[0]	=0.;
  ESaJ[0]	=0.;
  ESaJ1a[0]	=0.;
  bJ		=0.;
  dJ1		=0.;
  d2J1		=ESLc[0]*ESjb[0]/ESaF[0];

  ESaF1a[0]	=0.;
  dF1		=0.;
  d2F1		=EZcr2*R0*d2J1*(ESjp[0]-ESjb[0])/ESg22c[0];

  i1		=0;
  for(i=1; i < ESNa1; i++){
    dJ0		=dJ1;
    d2J0	=d2J1;
    dF0		=dF1;
    d2F0	=d2F1;
    d2J1	=ESLc[i]*ESjb[i]/ESaF[i];
    dJ1		=ESsa[i]*d2J1;
    d2J1	+=ESsa[i]*(ESLc[i]*ESjb1a[i]+ESLc1[i]*ESjb[i]
			-d2J1*ESaF1a[i])/ESaF[i];
    bJ		+=h2sh*(dJ0+dJ1)+h12p2sh*(d2J0-d2J1);

    ESaJ[i]	=bJ*ESaF[i];
    ESaJ1a[i]	=dJ1*ESaF[i]+bJ*ESaF1a[i];
    d2J0	=1./ESaJ[i];
    dJ0		=-R0*ESsa[i]*ESg22c[i]*d2J0;
    ESjp[i]	=dJ0*ESsp1a[i];
    ESjp1a[i]	=dJ0*ESsp2a[i]
      -(R0*(ESg22c[i]+ESsa[i]*ESg22c1[i])*ESsp1a[i]+ESjp[i]*ESaJ1a[i])*d2J0;
    dJ0		=1./ESLc[i];
    d2J0	=ESg22c1[i]/ESg22c[i];
    ds		=R0*bJ*dJ0/(ESpa[i]*ESg22c[i]);
    s		=bJ*ds;
    ds		=2.*dJ1*ds-s*(d2J0+ESLc1[i]*dJ0+2./ESsa[i]);
    rs		=1./(1.+s);
    d2F1	=(1.+ESVc[i]*dJ0)*ESjp[i]-ESjb[i];
    dF1		=R0*bJ*d2F1*rs/(ESsa[i]*ESg22c[i]);
    
    ESjs1a[i]	=(ESVc1[i]-ESVc[i]*ESLc1[i]*dJ0)*dJ0*ESjp[i];
    d2F1	=-dF1*(1./ESsa[i]+d2J0+ds*rs)
      +R0*(dJ1*d2F1
	   +bJ*((1.+ESVc[i]*dJ0)*ESjp1a[i]-ESjb1a[i]+ESjs1a[i])
	   )*rs/(ESsa[i]*ESg22c[i]);
    ESaF[i]	=ESaF[i1]+h2sh*(dF0+dF1)+h12p2sh*(d2F0-d2F1);
    ESaF1a[i]	=dF1;
    
    ESjs[i]	=(ESjb[i]+(s-ESVc[i]*dJ0)*ESjp[i])*rs;
    ESjs1a[i]	=(ESjb1a[i]+(s-ESVc[i]*dJ0)*ESjp1a[i]-ESjs1a[i]
		  +ds*(ESjp[i]-ESjs[i]))*rs;
    i1++;
  }
  ds	=ESRBt-ESaF[ESNa];
  for(i=0; i < ESNa1; i++){
    ESaF[i]	+=ds;
  }
  splA(ESjs,ESjs2a,ESjs1a,ESjs1a+ESNa);
  splA(ESjp,ESjp2a,ESjp1a,ESjp1a+ESNa);
  splA(ESaJ,ESaJ2a,ESaJ1a,ESaJ1a+ESNa);
  for(i=0; i < ESNa1; i++){
    ESjB1a[i]	=ESLc[i]*ESaF[i]*rR0/(ESLc[i]+ESVc[i]);
    ESjB[i]	=ESjb[i]*ESjB1a[i];
    ESjB1a[i]	=ESjb1a[i]*ESjB1a[i]
      +ESjb[i]*(ESLc1[i]*ESaF[i]
		+ESLc[i]*(ESaF1a[i]
			  -ESaF[i]*(ESLc1[i]+ESVc1[i])/(ESLc[i]+ESVc[i]))
		)*rR0/(ESLc[i]+ESVc[i]);

    ESFF[i]	=ESaF[i]*ESaF[i];
    ESFF1a[i]	=2.*ESaF[i]*ESaF1a[i];

    ESaT[i]	=R0*(ESjs[i]-ESjp[i]);
    ESaT1a[i]	=R0*(ESjs1a[i]-ESjp1a[i]);
    ESaT2a[i]	=R0*(ESjs2a[i]-ESjp2a[i]);
    ESPs[i]	=rR0*ESjp[i];
    ESPs1a[i]	=rR0*ESjp1a[i];
    ESPs2a[i]	=rR0*ESjp2a[i];
  }
  splA(ESjB,ESjB2a,ESjB1a,ESjB1a+ESNa);
  splA(ESFF,ESFF2a,ESFF1a,ESFF1a+ESNa);
  splA(ESaF,ESaF2a,ESaF1a,ESaF1a+ESNa);
  return(0);
}

int ESEq1DjBsp()
{
  double s,ds,rs;
  double h2sh,h12p2sh;
  double bJ,dJ0,dJ1,d2J0,d2J1;
  double dF0,dF1,d2F0,d2F1;
  int i,i1;
  
  h2sh		=EZcr2*ESsa[1];
  h12p2sh	=EZcr12*ESpa[1];
  
  ESjb[0]	=R0*(ESLc[0]+ESVc[0])*ESjB[0]/(ESLc[0]*ESaF[0]);
  ESjb1a[0]	=0.;
  ESjs[0]	=ESjb[0];
  ESjs1a[0]	=0.;
  ESjp[0]	=-2.*R0*ESg22c[0]*ESsp2a[0]/(ESjb[0]*ESLc[0]);
  ESjp1a[0]	=0.;
  ESaJ[0]	=0.;
  ESaJ1a[0]	=0.;
  bJ		=0.;
  dJ1		=0.;
  d2J1		=ESLc[0]*ESjb[0]/ESaF[0];

  ESaF[0]	=ESaF[0];
  ESaF1a[0]	=0.;
  dF1		=0.;
  d2F1		=EZcr2*R0*d2J1*(ESjp[0]-ESjb[0])/ESg22c[0];

  i1		=0;
  for(i=1; i < ESNa1; i++){
    dJ0		=dJ1;
    d2J0	=d2J1;
    dF0		=dF1;
    d2F0	=d2F1;

    dJ1		=1./(ESLc[i]*ESaF[i]);
    d2J1	=ESLc[i]+ESVc[i];
    ESjb[i]	=R0*d2J1*ESjB[i]*dJ1;
    ESjb1a[i]	=R0*((ESLc1[i]+ESVc1[i])*ESjB[i]+d2J1*ESjB1a[i])*dJ1
      -ESjb[i]*(ESLc1[i]*ESaF[i]+ESLc[i]*ESaF1a[i])*dJ1;

    dF1		=1./ESaF[i];
    d2J1	=ESLc[i]*ESjb[i]*dF1;
    dJ1		=ESsa[i]*d2J1;
    d2J1	+=ESsa[i]*(ESLc[i]*ESjb1a[i]+ESLc1[i]*ESjb[i]
			   -d2J1*ESaF1a[i])*dF1;
    bJ		+=h2sh*(dJ0+dJ1)+h12p2sh*(d2J0-d2J1);

    ESaJ[i]	=bJ*ESaF[i];
    ESaJ1a[i]	=dJ1*ESaF[i]+bJ*ESaF1a[i];

    ESaJ[i]	=bJ*ESaF[i];
    ESaJ1a[i]	=dJ1*ESaF[i]+bJ*ESaF1a[i];
    d2J0	=1./ESaJ[i];
    dJ0		=-R0*ESsa[i]*ESg22c[i]*d2J0;
    ESjp[i]	=dJ0*ESsp1a[i];
    ESjp1a[i]	=dJ0*ESsp2a[i]
      -(R0*(ESg22c[i]+ESsa[i]*ESg22c1[i])*ESsp1a[i]+ESjp[i]*ESaJ1a[i])*d2J0;
    
    dJ0		=1./ESLc[i];
    d2J0	=ESg22c1[i]/ESg22c[i];
    ds		=R0*bJ*dJ0/(ESpa[i]*ESg22c[i]);
    s		=bJ*ds;
    ds		=2.*dJ1*ds-s*(d2J0+ESLc1[i]*dJ0+2./ESsa[i]);
    rs		=1./(1.+s);
    d2F1	=(1.+ESVc[i]*dJ0)*ESjp[i]-ESjb[i];
    dF1		=R0*bJ*d2F1*rs/(ESsa[i]*ESg22c[i]);
    
    ESjs1a[i]	=(ESVc1[i]-ESVc[i]*ESLc1[i]*dJ0)*dJ0*ESjp[i];
    d2F1	=-dF1*(1./ESsa[i]+d2J0+ds*rs)
      +R0*(dJ1*d2F1
	   +bJ*((1.+ESVc[i]*dJ0)*ESjp1a[i]-ESjb1a[i]+ESjs1a[i])
	   )*rs/(ESsa[i]*ESg22c[i]);
    ESaF[i]	=ESaF[i1]+h2sh*(dF0+dF1)+h12p2sh*(d2F0-d2F1);
    ESaF1a[i]	=dF1;
    
    ESjs[i]	=(ESjb[i]+(s-ESVc[i]*dJ0)*ESjp[i])*rs;
    ESjs1a[i]	=(ESjb1a[i]+(s-ESVc[i]*dJ0)*ESjp1a[i]-ESjs1a[i]
		  +ds*(ESjp[i]-ESjs[i]))*rs;
    i1++;
  }
  ds	=ESRBt-ESaF[ESNa];
  for(i=0; i < ESNa1; i++){
    ESaF[i]	+=ds;
  }
  
  splA(ESjs,ESjs2a,ESjs1a,ESjs1a+ESNa);
  splA(ESjp,ESjp2a,ESjp1a,ESjp1a+ESNa);
  splA(ESaJ,ESaJ2a,ESaJ1a,ESaJ1a+ESNa);
  for(i=0; i < ESNa1; i++){
    ESFF[i]	=ESaF[i]*ESaF[i];
    ESFF1a[i]	=2.*ESaF[i]*ESaF1a[i];

    ESaT[i]	=R0*(ESjs[i]-ESjp[i]);
    ESaT1a[i]	=R0*(ESjs1a[i]-ESjp1a[i]);
    ESaT2a[i]	=R0*(ESjs2a[i]-ESjp2a[i]);
    ESPs[i]	=rR0*ESjp[i];
    ESPs1a[i]	=rR0*ESjp1a[i];
    ESPs2a[i]	=rR0*ESjp2a[i];
  }
  splA(ESjb,ESjb2a,ESjb1a,ESjb1a+ESNa);
  splA(ESFF,ESFF2a,ESFF1a,ESFF1a+ESNa);
  splA(ESaF,ESaF2a,ESaF1a,ESaF1a+ESNa);
  return(0);
}

int ESEqSolv1Djbsp()
{
  int i;

  RK1Dh	=ESsa[1];
  rR0	=1./R0;

  switch(ESEqSolvInCr){
  case 1:
    ESEq1Djbsp();
     for(i=0; i < ESNa1; i++){
      ESjR[i]	=R0*ESjb[i];
      ESjR1a[i]	=R0*ESjb1a[i];
      ESjR2a[i]	=R0*ESjb2a[i];
    }
   break;
  case 2:
    for(i=0; i < ESNa1; i++){
      ESjb[i]	=rR0*ESjR[i];
      ESjb1a[i]	=rR0*ESjR1a[i];
      ESjb2a[i]	=rR0*ESjR2a[i];
    }
    ESEq1Djbsp();
    break;
  case 4:
    ESEq1DjBsp();
    break;
  }
  {
    int i1;
    double h2sh,h12p2sh;
    double dgY0,dgY1,d2gY0,d2gY1;

    h2sh	=EZcr2*ESsa[1];
    h12p2sh	=EZcr12*ESpa[1];

    ESgY[0]	=0.;
    ESgY1a[0]	=0.;
    dgY1	=0.;
    d2gY1	=-EZcr2*ESjs[0]*ESLc[0]/ESg22c[0];
    ESsq[0]	=2.*rR0*ESg22c[0]*ESaF[0]/ESjs[0];
    ESsq1a[0]	=0.;
    ESgm[0]	=1./ESsq[0];
    ESgm1a[0]	=0.;
    i1	=0;
    for(i=1; i < ESNa1; i++){
      dgY0	=dgY1;
      d2gY0	=d2gY1;
      d2gY1	=ESpa[i]*ESg22c[i]*ESLc[i]/ESaJ[i];
      ESsq[i]	=rR0*d2gY1;
      ESsq1a[i]	=rR0*(2.*ESsa[i]*ESg22c[i]*ESLc[i]-d2gY1*ESaJ1a[i]
		      +ESpa[i]*(ESg22c1[i]*ESLc[i]+ESg22c[i]*ESLc1[i])
		      )/ESaJ[i];
      ESsq1a[i]	=ESsq1a[i]*ESaF[i]+ESsq[i]*ESaF1a[i];
      ESsq[i]	*=ESaF[i];
      ESgm[i]	=1./ESsq[i];
      ESgm1a[i]	=-ESsq1a[i]*ESgm[i]*ESgm[i];
      
      d2gY1	=1./(ESsa[i]*ESg22c[i]);
      dgY1	=-ESaJ[i]*d2gY1;
      d2gY1	=(ESaJ[i]*(1./ESsa[i]+ESg22c1[i]/ESg22c[i])-ESaJ1a[i])*d2gY1;
      ESgY1a[i]	=dgY1;
      ESgY[i]	=ESgY[i1]+h2sh*(dgY0+dgY1)+h12p2sh*(d2gY0-d2gY1);
      i1++;
    }
  }

  splA(ESgY,ESgY2a,ESgY1a,ESgY1a+ESNa);
  splA(ESsq,ESsq2a,ESsq1a,ESsq1a+ESNa);
  splA(ESgm,ESgm2a,ESgm1a,ESgm1a+ESNa);

  for(i=0; i < ESNa1; i++){
    ESjR[i]	=ESjb[i]*R0;
    ESjR1a[i]	=ESjb1a[i]*R0;
    ESjR2a[i]	=ESjb2a[i]*R0;
  }

  EC1DRadCSolv();
  return(0);
}

int ESEq1Dgm2bF()
{
  double RKy0,RKy1,RKy2,RKy3,RKy4,x;
  double gm,dgm,jp,bL,dL,bV,bK,dK;
  double K,R,s,ss;
  int i;

  i	=ESNa;
  x	=ESsa[i];
  ESaF[i]=ESRBt;
  gm	=ESgm[i];
  dgm	=ESgm1a[i];
  bL	=ESLc[i];
  dL	=ESLc1[i];
  bK	=ESg22c[i];
  dK	=ESg22c1[i];
  bV	=ESVc[i];
  jp	=ESjp[i];
  ss	=bK*bL*gm*gm*rR0;
  s	=1./(1.+x*x*ss);
  K	=(2.*ss+x*((dK*bL+bK*dL)*gm+bK*bL*dgm)*gm*rR0)*s;
  R	=gm*jp*(bL+bV)*s;
  ESaF1a[i]	=R-K*ESaF[i];
  ESjs[i]	=jp-ESaF1a[i]/(gm*bL);
  ESaF1a[i]	*=x;
  while(i > 0){
    RKy1	=-RK1Dh*ESaF1a[i];
    RKy0	=ESaF[i]+EZcr2*RKy1;
    x		-=EZcr2*RK1Dh;
    ESSetSplA(x);
    ESSetSplDPr(x);
    ESSetSplDCr(x);
    splRDCr(&gm,&dgm,7);
    splRDPr(&jp,NULL,ESEqSolvInPr);
    if(ESEqSolvInPr == 1){
      jp	*=R0;
    }
    splRA(&bL,&dL,ESLc,ESLc2);
    splRA(&bV,NULL,ESVc,ESVc2);
    splRA(&bK,&dK,ESg22c,ESg22c2);
    ss	=bK*bL*gm*gm*rR0;
    s	=x/(1.+x*x*ss);
    K	=(2.*ss+x*((dK*bL+bK*dL)*gm+bK*bL*dgm)*gm*rR0)*s;
    R	=gm*jp*(bL+bV)*s;
    RKy2=-RK1Dh*(R-K*RKy0);
    RKy0=ESaF[i]+EZcr2*RKy2;
    RKy3=-RK1Dh*(R-K*RKy0);
    RKy0=ESaF[i]+RKy3;
    i--;
    x	=ESsa[i];
    ESSetSplDCr(x);
    splRDCr(&gm,&dgm,7);
    bL	=ESLc[i];
    dL	=ESLc1[i];
    bK	=ESg22c[i];
    dK	=ESg22c1[i];
    bV	=ESVc[i];
    jp	=ESjp[i];

    ss	=bK*bL*gm*gm*rR0;
    s	=1./(1.+x*x*ss);
    K	=(2.*ss+x*((dK*bL+bK*dL)*gm+bK*bL*dgm)*gm*rR0)*s;
    R	=gm*jp*(bL+bV)*s;
    RKy4=-RK1Dh*(R-K*RKy0)*x;
    ESaF[i]	=ESaF[i+1]+EZcr6*(RKy1+RKy4)+EZcr3*(RKy2+RKy3);
    ESaF1a[i]	=R-K*ESaF[i];

    ESjs[i]	=gm > 1e-4 ? jp-ESaF1a[i]/(gm*bL)
      : x*bK*bL*dgm*rR0*ESaF[i]/bL-jp*bV/bL;
    ESaF1a[i]	*=x;
  }
  return(0);
}

int ESEqSolv1DgmjpHole()
{
  double s,ss;
  int i;
  
  rR0	=1./R0;
  RK1Dh	=ESsa[1]-ESsa[0];

  switch(ESEqSolvInPr){
  case 1:
    for(i=0; i < ESNa1; i++){
      ESjp[i]	=R0*ESPs[i];
      ESjp1a[i]	=R0*ESPs1a[i];
      ESjp2a[i]	=R0*ESPs2a[i];
    }
    break;
  default :
    for(i=0; i < ESNa1; i++){
      ESPs[i]	=rR0*ESjp[i];
      ESPs1a[i]	=rR0*ESjp1a[i];
      ESPs2a[i]	=rR0*ESjp2a[i];
    }
    break;
  }
  /* Determination $\bbF$ from 1D EqEq{*/
  ESEq1Dgm2bF();
  splAA(ESaF,ESaF1a,ESaF2a,ESaF1a,ESaF1a+ESNa);
  s	=ESg22c[0]*ESaF[0]*ESLc[0]*rR0;
  ss	=((ESg22c1[0]*ESaF[0]+ESg22c[0]*ESaF1a[0])*ESLc[0]
	  +ESg22c[0]*ESaF[0]*ESLc1[0])*rR0;
  ESjs1a[0]	=(3.*s+2.*ESa0*ss)*ESgm1a[0]+ESa0*s*ESgm2a[0];

  ESjs1a[ESNa]	=ESjp1a[ESNa]-ESaF2a[ESNa]/(ESgm[ESNa]*ESLc[ESNa])
    -(ESgm1a[ESNa]/ESgm[ESNa]+ESLc1[ESNa]/ESLc[ESNa])*ESjs[ESNa];
  splAA(ESjs,ESjs1a,ESjs2a,ESjs1a,ESjs1a+ESNa);
  /*}*/
  
  for(i=0; i < ESNa1; i++){
    ESaT[i]	=R0*(ESjs[i]-ESjp[i]);
    ESaT1a[i]	=R0*(ESjs1a[i]-ESjp1a[i]);
    ESaT2a[i]	=R0*(ESjs2a[i]-ESjp2a[i]);
    ESFF[i]	=ESaF[i]*ESaF[i];
    ESFF1a[i]	=2.*ESaF[i]*ESaF1a[i];
  }
  splAA(ESFF,ESFF1a,ESFF2a,ESFF1a,ESFF1a+ESNa);

  {
    int i1;
    double h2sh,h12p2sh;
    double dgY0,dgY1,d2gY0,d2gY1;
    double dsp0,dsp1,d2sp0,d2sp1;

    h12p2sh	=ESsa[1]-ESsa[0];
    h2sh	=EZcr2*h12p2sh;
    h12p2sh	*=EZcr12*h12p2sh;

    ESgY[0]	=0.;
    ESgY1a[0]	=0.;
    dgY1	=0.;
    d2gY1	=-ESa0*ESLc[0]*ESaF[0]*ESgm1a[0]*rR0;

    ESaJ[0]	=0.;
    ESaJ1a[0]	=-ESa0*ESg22c[0]*d2gY1;

    ESsp[0]	=0.;
    ESsp1a[0]	=0.;
    dsp1	=0.;
    d2sp1	=ESPs[0]*d2gY1;
    i1	=0;
    for(i=1; i < ESNa1; i++){
      dgY0	=dgY1;
      d2gY0	=d2gY1;
      dsp0	=dsp1;
      d2sp0	=d2sp1;

      dsp1	=ESLc[i]*ESaF[i]*ESgm[i]*rR0;
      d2sp1	=((ESLc1[i]*ESaF[i]+ESLc[i]*ESaF1a[i])*ESgm[i]
			+ESLc[i]*ESaF[i]*ESgm1a[i])*rR0;
      d2gY1	=-dsp1;
      ESaJ[i]	=ESpa[i]*ESg22c[i]*dsp1;
      ESaJ1a[i]	=(2.*ESsa[i]*ESg22c[i]+ESpa[i]*ESg22c1[i])*dsp1
	+ESpa[i]*ESg22c[i]*d2sp1;

      dgY1	=ESsa[i]*d2gY1;
      d2gY1	-=ESsa[i]*d2sp1;
      ESgY1a[i]	=dgY1;
      ESgY[i]	=ESgY[i1]+h2sh*(dgY0+dgY1)+h12p2sh*(d2gY0-d2gY1);

      dsp1	=ESPs[i]*dgY1;
      d2sp1	=ESPs1a[i]*dgY1+ESPs[i]*d2gY1;
      ESsp1a[i]	=dsp1;
      ESsp[i]	=ESsp[i1]+h2sh*(dsp0+dsp1)+h12p2sh*(d2sp0-d2sp1);
      i1++;
    }
    dsp1	=-ESsp[ESNa];
    for(i=0; i < ESNa1; i++){
      ESsp[i]	+=dsp1;
    }
  }
  splAA(ESaJ,ESaJ1a,ESaJ2a,ESaJ1a,ESaJ1a+ESNa);
  splAA(ESsp,ESsp1a,ESsp2a,ESsp1a,ESsp1a+ESNa);
  splAA(ESgY,ESgY1a,ESgY2a,ESgY1a,ESgY1a+ESNa);

  ESjb[0]	=ESjs[0];
  ESjb1a[0]	=0.;
  ESjB[0]	=ESjb[0]*ESFF[0]*rR0;
  ESjB1a[0]	=0.;

  for(i=0; i < ESNa1; i++){
    ESjb2a[i]	=ESg22c[i]*ESaT[i]*ESgY1a[i]*ESgY1a[i]/ESFF[i];
    ESjb1a[i]	=(ESVc[i]*ESjp[i]+ESjb2a[i])/ESLc[i];
    ESjb[i]	=ESjs[i]+ESjb1a[i];
    ESjb1a[i]	=ESjs1a[i]+
      (ESVc[i]*ESjp1a[i]+ESVc1[i]*ESjp[i]
       +((ESaT1a[i]*ESg22c[i]*ESgY1a[i]
	  +ESaT[i]*(ESg22c1[i]*ESgY1a[i]+2.*ESg22c[i]*ESgY2a[i]))*ESgY1a[i]
	 -ESjb2a[i]*ESFF1a[i]
	 )/ESFF[i]
       -ESjb1a[i]*ESLc1[i]
       )/ESLc[i];
    ESjB1a[i]	=ESLc[i]*ESFF[i]*rR0/(ESLc[i]+ESVc[i]);
    ESjB[i]	=ESjb[i]*ESjB1a[i];
    ESjB1a[i]	=ESjb1a[i]*ESjB1a[i]
      +ESjb[i]*(ESLc1[i]*ESFF[i]
		+ESLc[i]*(ESFF1a[i]
			  -ESFF[i]*(ESLc1[i]+ESVc1[i])/(ESLc[i]+ESVc[i]))
		)*rR0/(ESLc[i]+ESVc[i]);
  }
  splAA(ESjb,ESjb1a,ESjb2a,ESjb1a,ESjb1a+ESNa);
  splAA(ESjB,ESjB1a,ESjB2a,ESjB1a,ESjB1a+ESNa);

  for(i=0; i < ESNa1; i++){
    ESjR[i]	=ESjb[i]*R0;
    ESjR1a[i]	=ESjb1a[i]*R0;
    ESjR2a[i]	=ESjb2a[i]*R0;
  }

  EC1DRadCSolvHole();
  return(0);
}

int ESEqSolv1Dgmjp()
{
  double s,ss;
  int i;
  
  rR0	=1./R0;
  RK1Dh	=ESsa[1]-ESsa[0];

  if(ESa0 != 0.){
    ESEqSolv1DgmjpHole();
    return(0);
  }

  switch(ESEqSolvInPr){
  case 1:
    for(i=0; i < ESNa1; i++){
      ESjp[i]	=R0*ESPs[i];
      ESjp1a[i]	=R0*ESPs1a[i];
      ESjp2a[i]	=R0*ESPs2a[i];
    }
    break;
  default :
    for(i=0; i < ESNa1; i++){
      ESPs[i]	=rR0*ESjp[i];
      ESPs1a[i]	=rR0*ESjp1a[i];
      ESPs2a[i]	=rR0*ESjp2a[i];
    }
    break;
  }
  /* Determination $\bbF$ from 1D EqEq{*/
  ESEq1Dgm2bF();
  splA(ESaF,ESaF2a,ESaF1a,ESaF1a+ESNa);
  
  ESjs1a[0]	=0.;
  ESjs1a[ESNa]	=ESjp1a[ESNa]-ESaF2a[ESNa]/(ESgm[ESNa]*ESLc[ESNa])
    -(ESgm1a[ESNa]/ESgm[ESNa]+ESLc1[ESNa]/ESLc[ESNa])*ESjs[ESNa];
  splAA(ESjs,ESjs1a,ESjs2a,ESjs1a,ESjs1a+ESNa);
  /*}*/
  
  for(i=0; i < ESNa1; i++){
    ESaT[i]	=R0*(ESjs[i]-ESjp[i]);
    ESaT1a[i]	=R0*(ESjs1a[i]-ESjp1a[i]);
    ESaT2a[i]	=R0*(ESjs2a[i]-ESjp2a[i]);
    ESFF[i]	=ESaF[i]*ESaF[i];
    ESFF1a[i]	=2.*ESaF[i]*ESaF1a[i];
  }
  splA(ESFF,ESFF2a,ESFF1a,ESFF1a+ESNa);

  {
    int i1;
    double h2sh,h12p2sh;
    double dgY0,dgY1,d2gY0,d2gY1;
    double dsp0,dsp1,d2sp0,d2sp1;

    h2sh	=EZcr2*ESsa[1];
    h12p2sh	=EZcr12*ESpa[1];

    ESaJ[0]	=0.;
    ESaJ1a[0]	=0.;

    ESgY[0]	=0.;
    ESgY1a[0]	=0.;
    dgY1	=0.;
    d2gY1	=-ESLc[0]*ESaF[0]*ESgm[0]*rR0;

    ESsp[0]	=0.;
    ESsp1a[0]	=0.;
    dsp1	=0.;
    d2sp1	=ESPs[0]*d2gY1;
    i1	=0;
    for(i=1; i < ESNa1; i++){
      dgY0	=dgY1;
      d2gY0	=d2gY1;
      dsp0	=dsp1;
      d2sp0	=d2sp1;

      dsp1	=ESLc[i]*ESaF[i]*ESgm[i]*rR0;
      d2sp1	=((ESLc1[i]*ESaF[i]+ESLc[i]*ESaF1a[i])*ESgm[i]
			+ESLc[i]*ESaF[i]*ESgm1a[i])*rR0;
      d2gY1	=-dsp1;
      ESaJ[i]	=ESpa[i]*ESg22c[i]*dsp1;
      ESaJ1a[i]	=(2.*ESsa[i]*ESg22c[i]+ESpa[i]*ESg22c1[i])*dsp1
	+ESpa[i]*ESg22c[i]*d2sp1;

      dgY1	=ESsa[i]*d2gY1;
      d2gY1	-=ESsa[i]*d2sp1;
      ESgY1a[i]	=dgY1;
      ESgY[i]	=ESgY[i1]+h2sh*(dgY0+dgY1)+h12p2sh*(d2gY0-d2gY1);

      dsp1	=ESPs[i]*dgY1;
      d2sp1	=ESPs1a[i]*dgY1+ESPs[i]*d2gY1;
      ESsp1a[i]	=dsp1;
      ESsp[i]	=ESsp[i1]+h2sh*(dsp0+dsp1)+h12p2sh*(d2sp0-d2sp1);
      i1++;
    }
  }
  splA(ESaJ,ESaJ2a,ESaJ1a,ESaJ1a+ESNa);
  splA(ESsp,ESsp2a,ESsp1a,ESsp1a+ESNa);
  splA(ESgY,ESgY2a,ESgY1a,ESgY1a+ESNa);
  ESsp[0]	=-ESsp[ESNa];
  ESjb[0]	=ESjs[0];
  ESjb1a[0]	=0.;
  ESjB[0]	=ESjb[0]*ESFF[0]*rR0;
  ESjB1a[0]	=0.;
  for(i=1; i < ESNa1; i++){
    ESsp[i]	+=ESsp[0];

    ESjb2a[i]	=ESg22c[i]*ESaT[i]*ESgY1a[i]*ESgY1a[i]/ESFF[i];
    ESjb1a[i]	=(ESVc[i]*ESjp[i]+ESjb2a[i])/ESLc[i];
    ESjb[i]	=ESjs[i]+ESjb1a[i];
    ESjb1a[i]	=ESjs1a[i]+
      (ESVc[i]*ESjp1a[i]+ESVc1[i]*ESjp[i]
       +((ESaT1a[i]*ESg22c[i]*ESgY1a[i]
	  +ESaT[i]*(ESg22c1[i]*ESgY1a[i]+2.*ESg22c[i]*ESgY2a[i]))*ESgY1a[i]
	 -ESjb2a[i]*ESFF1a[i]
	 )/ESFF[i]
       -ESjb1a[i]*ESLc1[i]
       )/ESLc[i];
    ESjB1a[i]	=ESLc[i]*ESFF[i]*rR0/(ESLc[i]+ESVc[i]);
    ESjB[i]	=ESjb[i]*ESjB1a[i];
    ESjB1a[i]	=ESjb1a[i]*ESjB1a[i]
      +ESjb[i]*(ESLc1[i]*ESFF[i]
		+ESLc[i]*(ESFF1a[i]
			  -ESFF[i]*(ESLc1[i]+ESVc1[i])/(ESLc[i]+ESVc[i]))
		)*rR0/(ESLc[i]+ESVc[i]);
  }
  splA(ESjb,ESjb2a,ESjb1a,ESjb1a+ESNa);
  splA(ESjB,ESjB2a,ESjB1a,ESjB1a+ESNa);
  for(i=0; i < ESNa1; i++){
    ESjR[i]	=ESjb[i]*R0;
    ESjR1a[i]	=ESjb1a[i]*R0;
    ESjR2a[i]	=ESjb2a[i]*R0;
  }
  EC1DRadCSolv();
  return(0);
}

int ESEq1Dgm2FF()
{
  double RKy0,RKy1,RKy2,RKy3,RKy4,x;
  double gm,dgm,dp,bL,dL,bV,bK,dK,p2R0;
  double h2sh,h12p2sh;
  double dgY0,dgY1,d2gY0,d2gY1;
  double K,R,s,ss;
  int i,i1;

  
  RK1Dh	=ESsa[1]-ESsa[0];
  rR0	=1./R0;
  s	=ESsa[1]-ESsa[0];
  h2sh		=-EZcr2*s;
  h12p2sh	=EZcr12*s*s;

  i	=ESNa;
  x	=ESsa[i];
  ESSetSplDPr(x);
  ESSetSplDCr(x);
  splRDCr2(&gm,&dgm,&d2gY1,7);
  splRDPr(&s,&dp,2);
  ESaF[i]=ESRBt;
  ESFF[i]=ESRBt*ESRBt;
  ESgm[i]	=gm;
  ESgm1a[i]	=dgm;

  bL	=ESLc[i];
  dL	=ESLc1[i];
  bK	=ESg22c[i];
  dK	=ESg22c1[i];
  bV	=ESVc[i];
  p2R0	=R0*R0;
  dp	=ESsp1a[i];

  ss	=bK*bL*gm*gm*rR0;
  s	=2./(1.+x*x*ss);
  K	=(2.*ss+x*((dK*bL+bK*dL)*gm+bK*bL*dgm)*gm*rR0)*s;
  R	=-p2R0*dp*(1.+bV/bL)*s/x;
  ESFF1a[i]	=R-K*ESFF[i];
  ESaF1a[i]	=EZcr2*ESFF1a[i]/ESRBt;
 
  s	=dL/bL;
  ss	=dK*gm+bK*dgm+bK*gm*s;
  dgY1		=2.*bK*gm+ss;
  ESjb[i]  	=rR0*ESaF[i]*dgY1;
  ESjb1a[i]  	=rR0*(ESaF1a[i]*dgY1
		      +ESaF[i]*((2.+s)*(dK*gm+bK*dgm)+ss
				+ESg22c2[i]*gm+dK*dgm
				+dK*dgm+bK*d2gY1
				+bK*gm*(ESLc2[i]/bL-s*s)
				));
  ESjb1a[0]  	=0.;
  
  ESjB[i]	=ESjb[i]*bL*ESFF[i]*rR0/(bL+bV);
  ESjB1a[i]	=(rR0*(ESjb1a[i]*bL*ESFF[i]+ESjb[i]*dL*ESFF[i]
		       +ESjb[i]*bL*ESFF1a[i])-ESjB[i]*(dL+ESVc1[i]))/(bL+bV);
  ESjB1a[0]	=0.;

  d2gY1	=-gm*bL*ESaF[i]*rR0;
  dgY1	=d2gY1;
  d2gY1	-=(gm*(bL*ESaF1a[i]+dL*ESaF[i])+dgm*bL*ESaF[i])*rR0;
  ESaJ[i]	=-bK*dgY1;
  ESgY[i]	=0.;
  ESgY1a[i]	=dgY1;
  RKy0		=1./dgY1;
  ESPs[i]	=dp*RKy0;
  ESaT[i]	=EZcr2*ESFF1a[i]*RKy0;
  i1		=ESNa;
  while(i > 0){
    RKy1	=-RK1Dh*ESFF1a[i];
    RKy0	=ESFF[i]+EZcr2*RKy1;
    x		-=EZcr2*RK1Dh;
    ESSetSplA(x);
    ESSetSplDPr(x);
    ESSetSplDCr(x);
    splRDCr(&gm,&dgm,7);
    splRDPr(&s,&dp,2);
    splRA(&bL,&dL,ESLc,ESLc2);
    splRA(&bV,NULL,ESVc,ESVc2);
    splRA(&bK,&dK,ESg22c,ESg22c2);
    ss	=bK*bL*gm*gm*rR0;
    s	=2./(1.+x*x*ss);
    K	=x*(2.*ss+x*((dK*bL+bK*dL)*gm+bK*bL*dgm)*gm*rR0)*s;
    R	=-p2R0*dp*(1.+bV/bL)*s;
    RKy2=-RK1Dh*(R-K*RKy0);
    RKy0=ESFF[i]+EZcr2*RKy2;
    RKy3=-RK1Dh*(R-K*RKy0);
    RKy0=ESFF[i]+RKy3;
    i--;
    x	=ESsa[i];
    ESSetSplDPr(x);
    ESSetSplDCr(x);
    splRDCr(&gm,&dgm,7);
    if(i){
      splRDPr(&s,&dp,2);
      dp	/=x;
    }
    else{
      splRDPr2(&s,NULL,&dp,2);
    }
    ESgm[i]	=gm;
    ESgm1a[i]	=dgm;
    bL	=ESLc[i];
    dL	=ESLc1[i];
    bK	=ESg22c[i];
    dK	=ESg22c1[i];
    bV	=ESVc[i];
    ss	=bK*bL*gm*gm*rR0;
    s	=2./(1.+x*x*ss);
    K	=(2.*ss+x*((dK*bL+bK*dL)*gm+bK*bL*dgm)*gm*rR0)*s;
    R		=-p2R0*dp*(1.+bV/bL)*s;
    RKy4	=-RK1Dh*(R-K*RKy0)*x;
    ESFF[i]	=ESFF[i+1]+EZcr6*(RKy1+RKy4)+EZcr3*(RKy2+RKy3);
    ESaF[i]	=sqrt(ESFF[i]);
    ESFF1a[i]	=R-K*ESFF[i];
    dgY0	=dgY1;
    d2gY0	=d2gY1;
    d2gY1	=-gm*bL*ESaF[i]*rR0;
    RKy0	=1./d2gY1;
    ESPs[i]	=dp*RKy0;
    ESaT[i]	=EZcr2*ESFF1a[i]*RKy0;
    ESFF1a[i]	*=x;
    ESaF1a[i]	=EZcr2*ESFF1a[i]/ESaF[i];

    ESjb[i]  	=rR0*ESaF[i]*(2.*bK*gm+x*(dK*gm+bK*dgm+bK*gm*dL/bL));

    dgY1	=x*d2gY1;
    d2gY1	-=x*(gm*(bL*ESaF1a[i]+dL*ESaF[i])+dgm*bL*ESaF[i])*rR0;
    ESaJ[i]	=-x*bK*dgY1;
    ESgY1a[i]	=dgY1;
    ESgY[i]	=ESgY[i1]+h2sh*(dgY0+dgY1)+h12p2sh*(d2gY0-d2gY1);
    i1--;
  }

  return(0);
}

int ESEqSolv1Dgmsp()
{
  double s,ss,RR;
  int i;
  
  rR0	=1./R0;
  RR	=R0*R0;


  ESEq1Dgm2FF();
  splA(ESgY,ESgY2a,ESgY1a,ESgY1a+ESNa);
  splA(ESjb,ESjb2a,ESjb1a,ESjb1a+ESNa);
  splA(ESjB,ESjB2a,ESjB1a,ESjB1a+ESNa);

#ifdef H
  EZout("sII","ESEqSolvInPr=",ESEqSolvInPr,ESEqSolvInCr);
  for(i=0; i < ESnDtPr; i++){
    EZout("sidd","Pr",i,lP->xd[i],lP->yd[i]);
  }


  ss	=ESg22c[0]*ESLc[0]*ESgm[0]*ESgm[0]*rR0;
  ESFF2a[0]	=4.*ss;
  ESFF1a[0]	=-2.*RR*ESsp2a[0]*(1.+ESVc[0]/ESLc[0]);
  for(i=1; i < ESNa1; i++){
    ss	=ESg22c[i]*ESLc[i]*ESgm[i]*ESgm[i]*rR0;
    s	=2./(1.+ESpa[i]*ss);
    ESFF2a[i]	=(2.*ss+((ESg22c1[i]*ESLc[i]
			  +ESg22c[i]*ESLc1[i])*ESgm[i]
			 +ESg22c[i]*ESLc[i]*ESgm1a[i]
			 )*ESgm[i]*ESsa[i]*rR0)*s;
    ESFF1a[i]	=-RR*ESsp1a[i]*(1.+ESVc[i]/ESLc[i])*s/ESsa[i];
  }
  ESFF[ESNa]	=ESRBt*ESRBt;
  ESsplA1D(ESFF,ESFF2a,ESFF2a,ESFF1a);
  splA1a(ESFF,ESFF1a,ESFF2a);
  if(ESFF[0] < 0.)
    return(-1);
  ESaF[0]	=sqrt(ESFF[0]);
  s		=ESLc[0]*ESaF[0]*ESgm[0]*rR0;
  ESaJ[0]	=0.;
  ESgY[0]	=-s;

  ESPs[0]	=-ESsp2a[0]/s;
  ESaT[0]	=-EZcr2*ESFF2a[0]/s;
  for(i=1; i < ESNa1; i++){
    if(ESFF[i] < 0.)
      return(-1);
    ESaF[i]	=sqrt(ESFF[i]);
    s		=ESLc[i]*ESaF[i]*ESgm[i]*rR0;
    ESaJ[i]	=ESpa[i]*s*ESg22c[i];
    ESgY[i]	=-s;
    s		=-1./(s*ESsa[i]);
    ESPs[i]	=ESsp1a[i]*s;
    ESaT[i]	=EZcr2*ESFF1a[i]*s;
  }
  ESaF1a[0]	=0.;
  ESaF1a[ESNa]	=EZcr2*ESFF1a[ESNa]/ESaF[ESNa];
  splAA(ESaF,ESaF1a,ESaF2a,ESaF1a,ESaF1a+ESNa);
#endif

  splAA(ESFF,ESFF1a,ESFF2a,ESFF1a,ESFF1a+ESNa);
  splAA(ESaF,ESaF1a,ESaF2a,ESaF1a,ESaF1a+ESNa);
  ESaJ1a[0]	=0.;
  ESaJ1a[ESNa]	=(2.+ESLc1[ESNa]/ESLc[ESNa]+
		  ESaF1a[ESNa]/ESaF[ESNa]+ESgm1a[ESNa]/ESgm[ESNa]
		  +ESg22c1[ESNa]/ESg22c[ESNa])*ESaJ[ESNa];
  splAA(ESaJ,ESaJ1a,ESaJ2a,ESaJ1a,ESaJ1a+ESNa);
  ESPs1a[0]	=0.;
  ESPs1a[ESNa]	=(ESsp2a[ESNa]-ESPs[ESNa]*ESgY2a[ESNa])/ESgY1a[ESNa];
  splAA(ESPs,ESPs1a,ESPs2a,ESPs1a,ESPs1a+ESNa);
  ESaT1a[0]	=0.;
  ESaT1a[ESNa]	=(EZcr2*ESFF2a[ESNa]-ESaT[ESNa]*ESgY2a[ESNa])/ESgY1a[ESNa];
  splAA(ESaT,ESaT1a,ESaT2a,ESaT1a,ESaT1a+ESNa);
  s	=ESgY[0];

  for(i=0; i < ESNa1; i++){
    ESgY[i]	-=s;
    ESjp[i]	=R0*ESPs[i];
    ESjp1a[i]	=R0*ESPs1a[i];
    ESjp2a[i]	=R0*ESPs2a[i];
    ESjs[i]	=ESjp[i]+rR0*ESaT[i];
    ESjs1a[i]	=ESjp1a[i]+rR0*ESaT1a[i];
    ESjs2a[i]	=ESjp2a[i]+rR0*ESaT2a[i];
    ESjR[i]	=ESjb[i]*R0;
    ESjR1a[i]	=ESjb1a[i]*R0;
    ESjR2a[i]	=ESjb2a[i]*R0;
  }
  EC1DRadCSolv();
  return(0);
}

extern double *EZvr,*EZvz;
double ESLR1,ESLR2,ESLrR2,ESLrZ2,ESLaP,ESLaT,ESLbgYb,ESLbgY0,ESLR0;

int ECSoloviev(double *rd,double *zd,int nd)
{
  double p2ri,p2re,p2rt,p2zt;
  double b,b0,c,c0,p2z,EZz0,sb0,sb,db;
  double *pr,*pz,gY;
  int i,j,j1,j2,j3;

  rd[1]	=rd[0];
  p2ri	=rd[2]*rd[2];
  p2re	=rd[3]*rd[3];
  p2rt	=rd[0]*rd[0];
  sb0	=EZcr2*(zd[0]-zd[1]);
  p2zt	=sb0*sb0;
  EZz0	=EZcr2*(zd[0]+zd[1]);
  zd[2]	=EZz0;
  zd[3]	=EZz0;
  
  ESLrZ2=p2ri+p2re-2.*p2rt;
  if(ESLrZ2 <= 0.){
    printf("rt=%10.3e parameter in Soloviev equilibrium >= %10.3e%c\n"
	   ,rd[0],sqrt(EZcr2*(p2ri+p2re)),7);
    return(1);
  }
  ESLR1		=(p2ri*p2re-p2rt*p2rt)/ESLrZ2;
  ESLR2		=p2ri+p2re-ESLR1;
  ESLrR2	=1./ESLR2;
  ESLrZ2	*=ESLrR2/p2zt;

  ESLaP		=2.*ESLrZ2;
  ESLaT		=-ESLaP*ESLR1;
  ESLaP		+=8.*ESLrR2;
  ESLR0		=EZcr2*(p2ri+p2re);
  ESLbgY0	=(1.-ESLR0*ESLrR2)*(ESLR0-ESLR1);
  ESLbgYb	=(1.-p2re*ESLrR2)*(p2re-ESLR1);
  ESLR0		=sqrt(ESLR0);

  printf("%5s %13s %10s %10s %10s %10s %10s\n","P=1."
	 ,"T","R0","bgY_int","bsp_0","F^2_0","F_0");
  db	=(ESLbgY0-ESLbgYb)/ESLaP;
  sb	=2.*ESLaT/ESLaP*db+ESRBt*ESRBt;
  if(sb > 0.){
    printf("%5.3f %13.6e %10.4e %10.4e %10.4e %10.4e %10.4e\n",1.
	   ,ESLaT/ESLaP,ESLR0,db,db,sb,sqrt(sb));
  }
  else{
    printf("%5.3f %13.6e %10.3e %10.3e %10.3e %10.3e F^2_0 < 0\n",1.
	   ,ESLaT/ESLaP,ESLR0,db,db,sb);
  }
  ESSolov2ES((double)1.,ESLaT/ESLaP);
  db	=sb0/ESNa;
  b0	=EZcr2*(1.+ESLR1*ESLrR2);
  for(j=0; j < ESNp1; j++){
    EZvz[j]	=EZz0;
    EZvr[j]	=ESLR0;
  }
  pr	=EZvr;
  pz	=EZvz;
  for(i=1; i < ESNa1; i++){
    pr		+=ESNp1;
    pz		+=ESNp1;
    sb		=db*i;
    if(ESEqSolvRadC){
      sb	=ESsb[i];
    }
    p2zt	=sb*sb;
    p2rt	=EZcr2*(p2ri+p2re-p2zt*ESLR2*ESLrZ2);
    gY		=(1.-p2rt*ESLrR2-p2zt*ESLrZ2)*(p2rt-ESLR1);
    c0		=gY+ESLR1;
    j1		=ESNp/2;
    j2		=j1;
    j3		=ESNp;
    for(j=0; j < j1; j++){
      p2z	=sb*ESsn1[j];
      pz[j]	=EZz0+p2z;
      pz[j1]	=EZz0+p2z;
      pz[j2]	=EZz0-p2z;
      pz[j3]	=EZz0-p2z;
      p2z	*=p2z*ESLrZ2;
      b		=b0-EZcr2*p2z;
      c		=c0-p2z*ESLR1;
      p2z	=sqrt(b*b-c*ESLrR2);
      pr[j]	=sqrt((b+p2z)*ESLR2);
      pr[j1]	=sqrt((b-p2z)*ESLR2);
      pr[j2]	=pr[j1];
      pr[j3]	=pr[j];
      j1--;
      j2++;
      j3--;
    }
    p2zt	=sqrt(p2zt);
    p2rt	=sqrt(p2rt);
    pz[j]	=EZz0+p2zt;
    pr[j]	=p2rt;
    pz[j2]	=EZz0-p2zt;
    pr[j2]	=pr[j];
  }
  j	=ESNp/8;
  rd[4]	=pr[j];
  zd[4]	=pz[j];
  j1	=ESNp/2-j;
  rd[5]	=pr[j1];
  zd[5]	=pz[j1];
  j1	=ESNp/2+j;
  rd[6]	=pr[j1];
  zd[6]	=pz[j1];
  j1	=ESNp-j;
  rd[7]	=pr[j1];
  zd[7]	=pz[j1];

  j1	=ESNp/2-6;
  rd[5]	=pr[j1];
  zd[5]	=pz[j1];
  j1	=ESNp/2+6;
  rd[6]	=pr[j1];
  zd[6]	=pz[j1];
  j1	=ESNp/2-12;
  rd[8]	=pr[j1];
  zd[8]	=pz[j1];
  j1	=ESNp/2+12;
  rd[9]	=pr[j1];
  zd[9]	=pz[j1];

  j	=6;
  rd[4]	=pr[j];
  zd[4]	=pz[j];
  j	=ESNp-6;
  rd[7]	=pr[j];
  zd[7]	=pz[j];
  j	=12;
  rd[10]=pr[j];
  zd[10]=pz[j];
  j	=ESNp-12;
  rd[11]=pr[j];
  zd[11]=pz[j];

  return(0);
}

int EC1DEqSolv()
{
  int (*lEqS)();

  switch(ESEqSolvInCr){
  case 0:
    if(ESEqSolvInPr == 2){
      lEqS	=ESEqSolv1Djssp;
      break;
    }
    lEqS	=ESEqSolv1Djsjp;
    break;
  case 1:
  case 2:
  case 4:
    if(ESEqSolvInPr == 2){
      lEqS	=ESEqSolv1Djbsp;
    }
    else{
     lEqS	= ESEqSolv1Djbjp;
    }
    break;
  case 3:
    lEqS	=ESEqSolv1Djsjp;
    break;
  case 6:
  case 7:
    if(ESEqSolvInPr == 2){
      lEqS	=ESEqSolv1Dgmsp;
    }
    else{
      lEqS	=ESEqSolv1Dgmjp;
    }
    break;
  case 10:
    lEqS	=ESEqSolv1DaFjp;
    break;
  default:
    if(ESEqSolvInPr == 2){
      lEqS	=ESEqSolv1Djssp;
      break;
    }
    lEqS	=ESEqSolv1Djsjp;
    break;
  }
  lEqS();
  
  if(ESkNormP == 0){
    if(ESkNormJ && ESEqSolvInCr < 3){
      int j;
      double u[2],v[2],xJ,sJ,aJ,dJ;
      xJ	=1.;
      aJ	=5.*ESaJ[ESNa];
      u[0]	=1.;
      u[1]	=ESIpl/aJ;
      if(fabs(u[1]-1.) < 1e-2){
	u[1]	=1.2;
      }
      j		=0;
      do{
	EZzero(&j,u,v,&sJ,dJ);
	ESScaleCrData(sJ/xJ);
	xJ	=sJ;
	ESPrSplData2CrPr(ESEqSolvInCr);
	lEqS();
	dJ	=ESIpl-5.*ESaJ[ESNa];
      }while(j > 0 && j < 10 && fabs(dJ/ESIpl) > 1e-6);
    }
    ESgbext	=gBext;
    ESgbN	=gBN;
  }
  else{
    int j;
    double u[3],v[3],xP,sP,db,b,gb0,errP;
    double *ld;
    
    switch(ESkNormP){
    case 1:
      gb0	=ESgbext;
      ld	=&gBext;
      errP	=1e-6;
      break;
    case 2:
      gb0	=ESgbN;
      ld	=&gBN;
      errP	=1e-4;
      break;
    }
    xP	=1.;
    u[0]=1.;
    u[1]=gb0/(*ld);
    if(fabs(u[1]-1.) < 1e-2){
      u[1]	=1.2;
    }
    j	=0;
    if(ESkNormJ == 0 || ESEqSolvInCr > 2){
      do{
	EZzero(&j,u,v,&sP,db);
	ESScalePrData(sP/xP);
	xP	=sP;
	ESPrSplData2PrPr(ESEqSolvInPr);
	lEqS();
	db	=gb0-(*ld);
      }while(j > 0 && j < 10 && fabs(db) > errP);
    }
    else{
      double xJ,sJ,aJ,dJ,err;
      xJ	=1.;
      aJ	=5.*ESaJ[ESNa];
      u[2]	=1.;
      v[0]	=1.;
      v[1]	=1.;
      v[2]	=ESIpl/aJ;
      if(fabs(v[2]-1.) < 1e-2){
	v[2]	=1.2;
      }
      err	=1e-4;
      do{
	zero2D(&j,u,v,&sP,&sJ,db,dJ,err);
	ESScalePrCrData(sP/xP,sJ/xJ);
	xP	=sP;
	xJ	=sJ;
	ESPrSplData2PrPr(ESEqSolvInPr);
	ESPrSplData2CrPr(ESEqSolvInCr);
	lEqS();
	db	=gb0-(*ld);
	dJ	=ESIpl-5.*ESaJ[ESNa];
      }while(j > 0 && j < 10 && (fabs(dJ/ESIpl) > 1e-5 || fabs(db) > errP));
    }
    switch(ESkNormP){
    case 1:
      ESgbN	=gBN;
      break;
    case 2:
      ESgbext	=gBext;
      break;
    }
  }
  return(0);
}

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "esc.h"

#include <time.h>
#include <sys/time.h>

extern time_t time(time_t *t);
extern struct tm *localtime(const time_t *timep);
static struct tm *ltm;
static time_t stTime;


#ifndef mzl_ESmain
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern unsigned int ESmemlev;
#endif

#ifndef mzl_1Dcore
extern int ESNa,ESNa1;
extern double *ESsa,*ESqa,*ESpa,ESa0;
extern double *ESdgY,*ESdgY1a,*ESdgY2a;
extern double *ESFF,*ESFF1a,*ESFF2a;
#endif

#ifndef mzl_2Dcore
extern int ESNp,ESNp1,ESnAP;
extern double *ESgt;
extern double *ESsr,*ESsz;
#endif

extern double ESgaG[4];

static double cr3ha,cha,crha,c2rha,cHa,crHa,chha,ch2ha,ch2h1,ch4h0,ch4h1;
static int iSplA=0,iSplA1;
static double SplAx,SplAxx,SplAX,SplAXX,*glA=NULL,*vA,*vXAA0,*vXAA1;

static double chp,crhp,chhp,ch2hp,cHp,crHp,cx6prhp,crd2fP,cprhp,cr6php;
static double *vP=NULL,*d2fP,*glP,*WP,*VP,*U0P[3],*U1P[3],*U2P[3];
static int iSplP=0,iSplP1;
static double SplPx,SplPxx,SplPX,SplPXX;

int ESGetSmoothAddr(double **EZga0)
{
  *EZga0=ESgaG;
  return(0);
}

int ESFirstSetSplA()
{
  double h;

  h	=ESsa[1]-ESsa[0];
  cha	=EZcr6*h;
  cr3ha	=EZcr3*h;
  crha	=1./h;
  c2rha	=2.*crha;
  chha	=EZcr2*h;
  ch2ha	=chha*cha;
  cHa	=EZcr6*h*h;
  crHa	=6.*crha*crha; 
  ch2h1	=EZcr3*h*h;
  ch4h1	=-EZcr3*h*h*h*h/30.;
  ch4h0	=7.*EZcr4*ch4h1;
  ch4h1	*=2.;
  return(0);
}

int ESInitSplA()
{
  ESFirstSetSplA();
  if(glA == NULL){
    glA		=(double*)malloc(12*ESNa1*sizeof(double));
    vA		=glA+ESNa1;
    vXAA0	=vA+9*ESNa1;
    vXAA1	=vXAA0+ESNa1;
  }
  ESFirstSetSplXAA();
  return(0);
}

int DeInitf2splA()
{
  if(glA != NULL){
    free(glA);
    glA		=NULL;
    iSplA	=0;
  }
  return(0);
}

int splA(double*g,double*d2g,double*dg0,double*dg1)
{
  static int i,ii,i1;
  double r0,r1;

  r0= (g[1]-g[0])*crHa;
  r1= (g[2]-g[1])*crHa;

  if(dg0 != NULL){
    vA[0]	=-EZcr2;
    d2g[0]	=EZcr2*r0-3.*(*dg0)*crha;
    vA[1]	=-1./(4.+vA[0]);
    d2g[1]	=(r0-r1+d2g[0])*vA[1];
  }
  else{
    vA[0]	=2.;
    d2g[0]	=0.;
    vA[1]	=0.;
    d2g[1]	=EZcr6*(r1-r0);
  }

  i1	=1;
  for(i=2; i < ESNa; i++){
    r0		=r1;
    r1		=(g[i+1]-g[i])*crHa;
    vA[i]	=-1./(4.+vA[i1]);
    d2g[i]	=(r0-r1+d2g[i1])*vA[i];
    i1++;
  }

  if(dg1 == NULL){
    d2g[i1]	=EZcr6*(r1-r0);
    vA[i1]	=0.;
    ii	=i1-1;
    d2g[i]	=(2.-vA[ii])*d2g[i1]-d2g[ii];
  }
  else{
    d2g[i]	=((*dg1)/cha-r1-d2g[i1])/(vA[i1]+2.);
  }
  while(i > 0){
    d2g[i1]	+=vA[i1]*d2g[i];
    i--;
    i1--;
  }
  if(dg0 == NULL){
    d2g[0]	-=d2g[2];
  }
  return(0);
}

int splA05(double*g,double*d2g,double*dg0,double *g05)
{
  static int i,ii,i1;
  double r0,r1;

  r0= (g[1]-g[0])*crHa;
  r1= (g[2]-g[1])*crHa;

  if(dg0 != NULL){
    vA[0]	=-EZcr2;
    d2g[0]	=EZcr2*r0-3.*(*dg0)*crha;
    vA[1]	=-1./(4.+vA[0]);
    d2g[1]	=(r0-r1+d2g[0])*vA[1];
  }
  else{
    vA[0]	=2.;
    d2g[0]	=0.;
    vA[1]	=0.;
    d2g[1]	=EZcr6*(r1-r0);
  }

  i1	=1;
  for(i=2; i < ESNa; i++){
    r0		=r1;
    r1		=(g[i+1]-g[i])*crHa;
    vA[i]	=-1./(4.+vA[i1]);
    d2g[i]	=(r0-r1+d2g[i1])*vA[i];
    i1++;
  }

  if(g05 == NULL){
    d2g[i1]	=EZcr6*(r1-r0);
    vA[i1]	=0.;
    ii	=i1-1;
    d2g[i]	=(2.-vA[ii])*d2g[i1]-d2g[ii];
  }
  else{
    d2g[i]	=(8.*crha*crha*(g[i]+g[i1]-2.*(*g05))-d2g[i1])/(vA[i1]+1.);
  }
  while(i > 0){
    d2g[i1]	+=vA[i1]*d2g[i];
    i--;
    i1--;
  }
  if(dg0 == NULL){
    d2g[0]	-=d2g[2];
  }
  return(0);
}

int f2splA(double*g,double*d2g,double*dg0,double*dg1,double*f)
{
  static int i,ii,k,Ind0;

  static double s,sc0,sc1,sc11,sd0,sd1,sd11,rh1,rh11,sq;
  static double r1,EZr2,r3;
  static double w[6],rw[6];
  static double t11,t13,t22,t23,t31,t32;
  static double bc00= 0.,bc01= 0.,bd00= 0.,bd01= 0.,be0= 1.,bf0= 0.;
  static double bcn0= 0.,bcn1= 0.,bdn0= 0.,bdn1= 0.,ben= 1.,bfn= 0.;

  i	=ESNa;
  rh1=ESsa[i]-ESsa[i-1];

  if(dg1!=NULL){
    bfn= (*dg1);
    bcn0= 1./rh1;
    bcn1= -bcn0;
    bdn0= EZcr3*rh1;
    bdn1= EZcr6*rh1;
    ben= 0.;
  }
  else{
    bfn= 0.;
    bcn0= 0.;
    bcn1= 0.;
    bdn0= 0.;
    bdn1= 0.;
    ben= 1.;
  }
  bfn-= bcn0*f[i];
  bcn0= 0.;
  rh1= ESsa[1]-ESsa[0];
  if(dg0!=NULL){
    bf0= (*dg0);
    bc01= 1./rh1;
    bc00= -bc01;
    bd00= -EZcr3*rh1;
    bd01= -EZcr6*rh1;
    be0= 0.;
  }
  else{
    bf0= 0.;
    bc01= 0.;
    bc00= 0.;
    bd00= 0.;
    bd01= 0.;
    be0= 1.;
  }
  for(i= 0;i<6;i++){
    w[i]= 0;
  }
  Ind0	=1;

  r1= 0.;
  EZr2= 0.;
  r3= 0.;
  for(i= 0,ii= 1,sc1= 0.,sd1= 0.,rh1= 0.;i < ESNa1;i++,ii++){
    k= 9*i;
    rh11= rh1;
    sc11= sc1;
    sd11= sd1;
    if(i<ESNa){
      rh1= ESsa[ii]-ESsa[i];
      sd1= EZcr6*rh1;
      sc1= 1./rh1;
      rh1= sc1*sc1;
    }
    if(i==0){
      sq= EZcr2*(ESsa[ii]-ESsa[i]);
      w[0]= ESgaG[1]*sc1+sq*(ESgaG[0]+1.);
      w[2]= bc00;
      w[3]= sq*ESgaG[2]+ESgaG[3]*sc1;
      w[4]= bd00;
      w[5]= be0;
      r3= bf0;
      if(Ind0){
	w[2]= 0.;
	r3-= bc00*f[0];
      }
    }
    else if(i<ESNa){
      sq= EZcr2*(ESsa[ii]-ESsa[i-1]);
      w[0]+= ESgaG[1]*(sc1+sc11)+sq*(ESgaG[0]+1.);
      w[2]+= sc1+sc11;
      w[3]+= sq*ESgaG[2]+ESgaG[3]*(sc1+sc11);
      w[4]+= EZcr3*(ESsa[ii]-ESsa[i-1]);
      if(i==1&&Ind0){
	r1+= ESgaG[1]*sc11*f[0];
	r3+= sc11*f[0];
      }
    }
    else{
      sq= EZcr2*(ESsa[i]-ESsa[i-1]);
      w[0]+= sq;
      w[3]+= sq*ESgaG[2]+ESgaG[3]*sc11;
      w[4]+= bdn0;
      w[5]+= ben;
      r3+= bfn;
    }
    r1+= sq*f[i];
    if(i==0){
      t11= -ESgaG[1]*sc1;
      t13= -sc1;
      t23= sd1;
      t31= bc01;
      t32= bd01;
      if(Ind0){
	w[0]= sq;
	t13= 0.;
	t11= 0.;
      }
    }
    else if(ii<ESNa){
      t11= -ESgaG[1]*sc1;
      t13= -sc1;
      t31= -sc1;
      t23= sd1;
      t32= sd1;
    }
    else{
      t11= 0.;
      t13= bcn1;
      t23= bdn1;
      t31= 0.;
      t32= sd1;
      if(i<ESNa){
	r1+= ESgaG[1]*sc1*f[ESNa];
	r3+= sc1*f[ESNa];
      }
    }
    t22= -ESgaG[3]*sc1;

    EZinv3x3(w,rw);

    g[i]= rw[0]*r1+rw[1]*EZr2+rw[2]*r3;
    d2g[i]= rw[1]*r1+rw[3]*EZr2+rw[4]*r3;
    glA[i]= rw[2]*r1+rw[4]*EZr2+rw[5]*r3;
    if(i<ESNa){
      vA[k+0]= -(rw[0]*t11+rw[2]*t31);
      vA[k+3]= -(rw[1]*t11+rw[4]*t31);
      vA[k+6]= -(rw[2]*t11+rw[5]*t31);
      vA[k+1]= -(rw[1]*t22+rw[2]*t32);
      vA[k+4]= -(rw[3]*t22+rw[4]*t32);
      vA[k+7]= -(rw[4]*t22+rw[5]*t32);
      vA[k+2]= -(rw[0]*t13+rw[1]*t23);
      vA[k+5]= -(rw[1]*t13+rw[3]*t23);
      vA[k+8]= -(rw[2]*t13+rw[4]*t23);

      w[0]= t11*vA[k+0]+t31*vA[k+6];
      w[1]= t11*vA[k+1]+t31*vA[k+7];
      w[2]= t11*vA[k+2]+t31*vA[k+8];
      w[3]= t22*vA[k+4]+t32*vA[k+7];
      w[4]= t13*vA[k+1]+t23*vA[k+4];
      w[5]= t13*vA[k+2]+t23*vA[k+5];


      r1= -(t11*g[i]+t31*glA[i]);
      EZr2= -(t22*d2g[i]+t32*glA[i]);
      r3= -(t13*g[i]+t23*d2g[i]);
    }
  }
  for(ii= ESNa,i= ii-1;ii>0;ii--,i--){
    k= 9*i;
    g[i]= g[i]+vA[k+0]*g[ii]+vA[k+1]*d2g[ii]+vA[k+2]*glA[ii];
    d2g[i]= d2g[i]+vA[k+3]*g[ii]+vA[k+4]*d2g[ii]+vA[k+5]*glA[ii];
    glA[i]= glA[i]+vA[k+6]*g[ii]+vA[k+7]*d2g[ii]+vA[k+8]*glA[ii];
  }
  return(0);
}

int f2splA2(double*g,double*d2g,double*d2g0,double*d1g1,double*f)
{
  static int i,ii,k,Ind0;

  static double s,sc0,sc1,sc11,sd0,sd1,sd11,rh1,rh11,sq;
  static double r1,EZr2,r3;
  static double w[6],rw[6];
  static double t11,t13,t22,t23,t31,t32;
  static double bc00= 0.,bc01= 0.,bd00= 0.,bd01= 0.,be0= 1.,bf0= 0.;
  static double bcn0= 0.,bcn1= 0.,bdn0= 0.,bdn1= 0.,ben= 1.,bfn= 0.;

  i	=ESNa;
  rh1=ESsa[i]-ESsa[i-1];

  if(d1g1!=NULL){
    bfn= (*d1g1);
    bcn0= 1./rh1;
    bcn1= -bcn0;
    bdn0= EZcr3*rh1;
    bdn1= EZcr6*rh1;
    ben= 0.;
  }
  else{
    bfn= 0.;
    bcn0= 0.;
    bcn1= 0.;
    bdn0= 0.;
    bdn1= 0.;
    ben= 1.;
  }
  bfn-= bcn0*f[i];
  bcn0= 0.;
  rh1= ESsa[1]-ESsa[0];
  if(d2g0 != NULL){
    bf0		=*d2g0;
    bc01	=0.;
    bc00	=0.;
    bd00	=1.;
    bd01	=0.;
    be0		=0.;
  }
  else{
    bf0		=0.;
    bc01	=0.;
    bc00	=0.;
    bd00	=0.;
    bd01	=0.;
    be0		=1.;
  }
  for(i= 0;i<6;i++){
    w[i]= 0;
  }
  Ind0	=1;

  r1= 0.;
  EZr2= 0.;
  r3= 0.;
  for(i= 0,ii= 1,sc1= 0.,sd1= 0.,rh1= 0.;i < ESNa1;i++,ii++){
    k= 9*i;
    rh11= rh1;
    sc11= sc1;
    sd11= sd1;
    if(i<ESNa){
      rh1= ESsa[ii]-ESsa[i];
      sd1= EZcr6*rh1;
      sc1= 1./rh1;
      rh1= sc1*sc1;
    }
    if(i==0){
      sq= EZcr2*(ESsa[ii]-ESsa[i]);
      w[0]= ESgaG[1]*sc1+sq*(ESgaG[0]+1.);
      w[2]= bc00;
      w[3]= sq*ESgaG[2]+ESgaG[3]*sc1;
      w[4]= bd00;
      w[5]= be0;
      r3= bf0;
      if(Ind0){
	w[2]= 0.;
	r3-= bc00*f[0];
      }
    }
    else if(i<ESNa){
      sq= EZcr2*(ESsa[ii]-ESsa[i-1]);
      w[0]+= ESgaG[1]*(sc1+sc11)+sq*(ESgaG[0]+1.);
      w[2]+= sc1+sc11;
      w[3]+= sq*ESgaG[2]+ESgaG[3]*(sc1+sc11);
      w[4]+= EZcr3*(ESsa[ii]-ESsa[i-1]);
      if(i==1&&Ind0){
	r1+= ESgaG[1]*sc11*f[0];
	r3+= sc11*f[0];
      }
    }
    else{
      sq= EZcr2*(ESsa[i]-ESsa[i-1]);
      w[0]+= sq;
      w[3]+= sq*ESgaG[2]+ESgaG[3]*sc11;
      w[4]+= bdn0;
      w[5]+= ben;
      r3+= bfn;
    }
    r1+= sq*f[i];
    if(i==0){
      t11= -ESgaG[1]*sc1;
      t13= -sc1;
      t23= sd1;
      t31= bc01;
      t32= bd01;
      if(Ind0){
	w[0]= sq;
	t13= 0.;
	t11= 0.;
      }
    }
    else if(ii<ESNa){
      t11= -ESgaG[1]*sc1;
      t13= -sc1;
      t31= -sc1;
      t23= sd1;
      t32= sd1;
    }
    else{
      t11= 0.;
      t13= bcn1;
      t23= bdn1;
      t31= 0.;
      t32= sd1;
      if(i<ESNa){
	r1+= ESgaG[1]*sc1*f[ESNa];
	r3+= sc1*f[ESNa];
      }
    }
    t22= -ESgaG[3]*sc1;

    EZinv3x3(w,rw);

    g[i]= rw[0]*r1+rw[1]*EZr2+rw[2]*r3;
    d2g[i]= rw[1]*r1+rw[3]*EZr2+rw[4]*r3;
    glA[i]= rw[2]*r1+rw[4]*EZr2+rw[5]*r3;
    if(i<ESNa){
      vA[k+0]= -(rw[0]*t11+rw[2]*t31);
      vA[k+3]= -(rw[1]*t11+rw[4]*t31);
      vA[k+6]= -(rw[2]*t11+rw[5]*t31);
      vA[k+1]= -(rw[1]*t22+rw[2]*t32);
      vA[k+4]= -(rw[3]*t22+rw[4]*t32);
      vA[k+7]= -(rw[4]*t22+rw[5]*t32);
      vA[k+2]= -(rw[0]*t13+rw[1]*t23);
      vA[k+5]= -(rw[1]*t13+rw[3]*t23);
      vA[k+8]= -(rw[2]*t13+rw[4]*t23);

      w[0]= t11*vA[k+0]+t31*vA[k+6];
      w[1]= t11*vA[k+1]+t31*vA[k+7];
      w[2]= t11*vA[k+2]+t31*vA[k+8];
      w[3]= t22*vA[k+4]+t32*vA[k+7];
      w[4]= t13*vA[k+1]+t23*vA[k+4];
      w[5]= t13*vA[k+2]+t23*vA[k+5];


      r1= -(t11*g[i]+t31*glA[i]);
      EZr2= -(t22*d2g[i]+t32*glA[i]);
      r3= -(t13*g[i]+t23*d2g[i]);
    }
  }
  for(ii= ESNa,i= ii-1;ii>0;ii--,i--){
    k= 9*i;
    g[i]= g[i]+vA[k+0]*g[ii]+vA[k+1]*d2g[ii]+vA[k+2]*glA[ii];
    d2g[i]= d2g[i]+vA[k+3]*g[ii]+vA[k+4]*d2g[ii]+vA[k+5]*glA[ii];
    glA[i]= glA[i]+vA[k+6]*g[ii]+vA[k+7]*d2g[ii]+vA[k+8]*glA[ii];
  }
  return(0);
}

int splA1a(double*g,double*d1g,double*d2g)
{
  static int i,ii;
  
  ii	=1;
  i	=0;
  d1g[0]=(g[ii]-g[i])*crha-(d2g[ii]+2.*d2g[i])*cha;
  while(ii < ESNa1){
    d1g[ii]	=(g[ii]-g[i])*crha+(2.*d2g[ii]+d2g[i])*cha;
    i++;
    ii++;
  }
  return(0);
}

int splAA(double*g,double*d1g,double*d2g,double*dg0,double*dg1)
{
  static int i,ii,i1;
  double r0,r1;

  r0= (g[1]-g[0])*crHa;
  r1= (g[2]-g[1])*crHa;

  if(dg0 != NULL){
    vA[0]	=-EZcr2;
    d2g[0]	=EZcr2*r0-3.*(*dg0)*crha;
    vA[1]	=-1./(4.+vA[0]);
    d2g[1]	=(r0-r1+d2g[0])*vA[1];
  }
  else{
    vA[0]	=2.;
    d2g[0]	=0.;
    vA[1]	=0.;
    d2g[1]	=EZcr6*(r1-r0);
  }

  i1	=1;
  for(i=2; i < ESNa; i++){
    r0		=r1;
    r1		=(g[i+1]-g[i])*crHa;
    vA[i]	=-1./(4.+vA[i1]);
    d2g[i]	=(r0-r1+d2g[i1])*vA[i];
    i1++;
  }

  if(dg1 == NULL){
    d2g[i1]	=EZcr6*(r1-r0);
    vA[i1]	=0.;
    ii	=	i1-1;
    d2g[i]	=(2.-vA[ii])*d2g[i1]-d2g[ii];
  }
  else{
    d2g[i]	=((*dg1)/cha-r1-d2g[i1])/(vA[i1]+2.);
  }
  while(i > 0){
    d2g[i1]	+=vA[i1]*d2g[i];
    i--;
    i1--;
  }
  i1	=1;
  if(dg0 == NULL){
    d2g[0]	-=d2g[2];
    d1g[0]	=(g[i1]-g[i])*crha-(d2g[i1]+2.*d2g[i])*cha;
  }
  else{
    d1g[0]	=*dg0;
  }
  while(i1 < ESNa){
    i++;
    i1++;
    d1g[i]	=(g[i1]-g[i])*crha-(d2g[i1]+2.*d2g[i])*cha;
  }
  if(dg1 == NULL){
    d1g[i1]	=(g[i1]-g[i])*crha+(2.*d2g[i1]+d2g[i])*cha;
  }
  else{
    d1g[i1]	=*dg1;
  }
  return(0);
}

int splAA2(double*g,double*d1g,double*d2g,double *d2g0,double*dg1)
{
  static int i,ii,i1;
  double r0,r1;

  r1		=(g[1]-g[0])*crHa;
  vA[0]		=0.;
  d2g[0]	=*d2g0;

  i1	=0;
  for(i=1; i < ESNa; i++){
    r0		=r1;
    r1		=(g[i+1]-g[i])*crHa;
    vA[i]	=-1./(4.+vA[i1]);
    d2g[i]	=(r0-r1+d2g[i1])*vA[i];
    i1++;
  }
  if(dg1 == NULL){
    d2g[i1]	=EZcr6*(r1-r0);
    vA[i1]	=0.;
    ii	=i1-1;
    d2g[i]	=(2.-vA[ii])*d2g[i1]-d2g[ii];
  }
  else{
    d2g[i]	=((*dg1)/cha-r1-d2g[i1])/(vA[i1]+2.);
  }
  while(i1 > 0){
    d2g[i1]	+=vA[i1]*d2g[i];
    i--;
    i1--;
  }
  i	=0;
  i1	=1;
  while(i < ESNa){
    d1g[i]	=(g[i1]-g[i])*crha-(d2g[i1]+2.*d2g[i])*cha;
    i++;
    i1++;
  }
  if(dg1 == NULL){
    i1		=i-1;
    d1g[i]	=(g[i]-g[i1])*crha+(2.*d2g[i]+d2g[i1])*cha;
  }
  else{
    d1g[i]	=*dg1;
  }
  return(0);
}

int ESsplA1D(double*g,double*d2g,double*K,double*R)
{
  double RKy0,RKy1,RKy2,RKy3,RKy4;
  double x,r0,r1,h;
  int i,ii,i1;

  h	=-ESsa[1];
  i	=ESNa;
  x	=ESsa[i];
  ii	=i-1;
  i1	=i-2;
  while(ii > 0){
    RKy1	=h*x*(R[i]-K[i]*g[i]);
    RKy0	=g[i]+EZcr2*RKy1;
    x		+=EZcr2*h;
    r0		=0.125*(6.*R[ii]+3.*R[i]-R[i1]);
    r1		=0.125*(6.*K[ii]+3.*K[i]-K[i1]);
    RKy2	=h*x*(r0-r1*RKy0);
    RKy0	=g[i]+EZcr2*RKy2;
    RKy3	=h*x*(r0-r1*RKy0);
    RKy0	=g[i]+RKy3;
    x		=ESsa[ii];
    RKy4	=h*x*(R[ii]-K[ii]*RKy0);
    g[ii]	=g[i]+EZcr6*(RKy1+RKy4)+EZcr3*(RKy2+RKy3);
    i--;
    ii--;
    i1--;
  }
  i1	=2;
  RKy1	=h*x*(R[i]-K[i]*g[i]);
  RKy0	=g[i]+EZcr2*RKy1;
  x	+=EZcr2*h;
  r0	=0.125*(6.*R[i]+3.*R[ii]-R[i1]);
  r1	=0.125*(6.*K[i]+3.*K[ii]-K[i1]);
  RKy2	=h*x*(r0-r1*RKy0);
  RKy0	=g[i]+EZcr2*RKy2;
  RKy3	=h*x*(r0-r1*RKy0);
  RKy0	=g[i]+RKy3;
  x	=ESsa[ii];
  RKy4	=h*x*(R[ii]-K[ii]*RKy0);
  g[ii]	=g[i]+EZcr6*(RKy1+RKy4)+EZcr3*(RKy2+RKy3);

  r0	=0.;
  r1	=R[ESNa]-K[ESNa]*g[ESNa];

  r0= (g[1]-g[0])*crHa;
  r1= (g[2]-g[1])*crHa;

  vA[0]		=-EZcr2;
  d2g[0]	=EZcr2*r0;
  vA[1]		=-1./(4.+vA[0]);
  d2g[1]	=(r0-r1+d2g[0])*vA[1];

  i1	=1;
  for(i=2; i < ESNa; i++){
    r0		=r1;
    r1		=(g[i+1]-g[i])*crHa;
    vA[i]	=-1./(4.+vA[i1]);
    d2g[i]	=(r0-r1+d2g[i1])*vA[i];
    i1++;
  }
  d2g[i]	=((R[ESNa]-K[ESNa]*g[ESNa])/cha-r1-d2g[i1])/(vA[i1]+2.);
  while(i > 0){
    d2g[i1]	+=vA[i1]*d2g[i];
    i--;
    i1--;
  }
  return(0);
}

int ESCsplA1D(double*g,double*d2g,double*K,double*R)
{
  static int i,ii;
  double *v,*vv,D,w0,w1,w2,w3,h;

  /* The method of the routine is bad. Do not use it!!!;*/

  h		=ESsa[1];
  vv		=vA;
  D		=1./(1.-2.*K[0]*cHa);
  vv[0]		=D;
  vv[1]		=-cHa*D;
  vv[2]		=-K[0]*D;
  vv[3]		=K[0]*cHa*D;
  g[0]		=2.*vv[1]*R[0];
  d2g[0]	=R[0]*D;
  i	=0;
  ii	=1;
  while(ii < ESNa){
    v		=vv;
    vv		+=4;
    w2		=h*ESsa[ii]*K[ii];
    w0		=1.+w2-v[0]+cHa*v[2];
    D		=2.*cHa;
    w1		=D-v[1]+cHa*v[3];
    w2		-=1.;
    w3		=-D;
    D		=1./(w0*w3-w1*w2);
    vv[0]	=w1*D;
    vv[1]	=-w1*cHa;
    vv[2]	=-w0*D;
    vv[3]	=w0*cHa;
    d2g[ii]	=h*ESsa[ii]*R[ii]*D;
    D		=d2g[ii]+(g[i]-cHa*d2g[i])*D;
    g[ii]	=w3*D-w1*d2g[ii];
    d2g[ii]	=-w2*D+w0*d2g[ii];
    i++;
    ii++;
  }
  w1		=1./(2.*cHa-vv[1]+cHa*vv[3]);
  w0		=-(1.+h*K[ii]-vv[0]+cHa*vv[2])*w1;
  d2g[ii]	=w1*(h*R[ii]+g[i]-cHa*d2g[i])+w0*g[ii];
  while(ii > 0){
    g[i]	+=vv[0]*g[ii]+vv[1]*d2g[ii];
    d2g[i]	+=vv[2]*g[ii]+vv[3]*d2g[ii];
    ii--;
    i--;
    vv	-=4;
  }
  return(0);
}

int ESFC0HrmA(double *g,double *f,double *df)
{
  int i,ii;

  g[0]	=0.;
  i	=0;
  for(ii=1; ii < ESNa1; ii++){
    g[ii]	=g[i]+chha*(f[i]+f[ii])+EZcr2*cHa*(df[i]-df[ii]);
    i++;
  }
  return(0);
}

int ESsplA1Es(double *g,double *d2g,double *R)
{
  int i,ii;
  double s,ss;

  d2g[0]=R[0];
  g[0]	=0.;
  i	=0;
  ss	=0.;
  for(ii=1; ii < ESNa1; ii++){
    s		=ss;
    ss		=ESsa[ii]*R[ii];
    g[ii]	=g[i]+cr3ha*(2.*s+ss)+cHa*d2g[i];
    d2g[ii]	=c2rha*(ss-s)-d2g[i];
    i++;
  }
  return(0);
}

int ESsplA1E(double *g,double *d2g,double *R)
{
  double r0,r1;
  int i,ii,i1;

  r1	=R[0];
  d2g[0]=R[0];
  d2g[ESNa]	=R[ESNa]/cha;
  g[0]	=0.;
  i	=0;
  i1	=2;
  for(ii=1; ii < ESNa; ii++){
    r0	=r1;
    r1	=R[ii];
    g[ii]	=g[i]+cha*(ESsa[i]*r0+ESsa[ii]*r1
			   +EZcr4*(ESsa[i]+ESsa[ii])*(6.*r1+3.*r0-R[i1]));
    i++;
    i1++;
  }
  g[ii]	=g[i]+cha*(ESsa[i]*r1+ESsa[ii]*R[ii]
		   +EZcr4*(ESsa[i]+ESsa[ii])*(6.*r1+3.*R[ii]-r0));
  r0= (g[1]-g[0])*crHa;
  r1= (g[2]-g[1])*crHa;
  vA[0]	=0.;
  vA[1]	=-1./(4.+vA[0]);
  d2g[1]	=(r0-r1+d2g[0])*vA[1];
  i1	=1;
  for(i=2; i < ESNa; i++){
    r0		=r1;
    r1		=(g[i+1]-g[i])*crHa;
    vA[i]	=-1./(4.+vA[i1]);
    d2g[i]	=(r0-r1+d2g[i1])*vA[i];
    i1++;
  }
  d2g[i]	=(d2g[ESNa]-r1-d2g[i1])/(vA[i1]+2.);
  while(i > 0){
    d2g[i1]	+=vA[i1]*d2g[i];
    i--;
    i1--;
  }
  return(0);
}

int ESsplAE(double *f,double *g,double *d2g)
{
  int i,ii;
  double eg;

  eg	=0.;
  f[0]	=EZcr2*g[0];
  ii	=1;
  for(i=0; i < ESNa; i++){
    eg	+=ESsa[i]*chha*(g[i]+g[ii]-ch2ha*(d2g[i]+d2g[ii]));
    eg	+=cHa*g[i]+ch2h1*g[ii]+ch4h0*d2g[i]+ch4h1*d2g[ii];
    f[ii]=eg/(ESpa[ii]);
    ii++;
  }
  return(0);
}

int ESSetSplA(double A)
{
  int iA;
  int i,j,kT1;
  
  if(A > 1.) A=1.;
  if(A < ESa0) A=ESa0;
  while(iSplA < ESNa && ESsa[iSplA+1] < A) iSplA++;
  while(iSplA > 0 && ESsa[iSplA] > A) iSplA--;
  iSplA1	=iSplA+1;
  SplAx		=(A-ESsa[iSplA])*crha;
  SplAxx	=SplAx*SplAx;
  SplAX		=1.-SplAx;
  SplAXX	=SplAX*SplAX;
  return(0);
}

int ESSetSplEA(double *ef, double*g, double*d2g)
{
  int i,ii;
  double g0,g1,d2g0,d2g1;

  ef[0]	=0.;
  g1	=g[0];
  d2g1	=d2g[0];
  for(i=1,ii=0; i < ESNa1; i++,ii++){
    g0		=g1;
    d2g0	=d2g1;
    g1		=g[i];
    d2g1	=d2g[i];
    ef[i]	=ef[ii]+chha*ESsa[ii]*(g0+g1-ch2ha*(d2g0+d2g1))
      +cHa*((g0+2.*g1)-cHa*(.7*d2g0+.8*d2g1));
  }
  return(0);
}

int splREA(double*ef, double*f, double*df, double*eg, double*g, double*d2g)
{
  double a,g0,g1,d2g0,d2g1;

  g0	=g[iSplA];
  g1	=g[iSplA1];
  d2g0	=d2g[iSplA];
  d2g1	=d2g[iSplA1];

  a	=ESsa[iSplA];
  *ef	=eg[iSplA]+
    +a*chha*(g0*SplAx*(1.+SplAX)+g1*SplAxx
	     -ch2ha*SplAxx*(d2g0*(1.+SplAX)*(1.+SplAX)+d2g1*(2.-SplAxx)))
      +cHa*SplAxx*(g0*(3.-2.*SplAx)+2.*g1*SplAx
		   +cHa*SplAx*(d2g0*(4.5*SplAx-1.2*SplAxx-4.)
			       +d2g1*(1.2*SplAxx-2.)));
  if(f != NULL) *f=SplAX*g0+SplAx*g1
    +cHa*(SplAX*(SplAXX-1.)*d2g0+SplAx*(SplAxx-1.)*d2g1);
  if(df !=NULL) *df=(g1-g0+((3.*SplAxx-1.)*d2g1-(3.*SplAXX-1.)*d2g0)*cHa)*crha;
  return(0);
}

int splRA(double*f, double*df, double*g, double*d2g)
{
  *f	=SplAX*g[iSplA]+SplAx*g[iSplA1]
    +(SplAX*(SplAXX-1.)*d2g[iSplA]+SplAx*(SplAxx-1.)*d2g[iSplA1])*cHa;
  if(df != NULL){
    *df	=(g[iSplA1]-g[iSplA]
	  +((3.*SplAxx-1.)*d2g[iSplA1]-(3.*SplAXX-1.)*d2g[iSplA])*cHa)*crha;
  }
  return(0);
}

int splRA2(double*f, double*df, double*d2f, double*g, double*d2g)
{
  *f	=SplAX*g[iSplA]+SplAx*g[iSplA1]
    +(SplAX*(SplAXX-1.)*d2g[iSplA]+SplAx*(SplAxx-1.)*d2g[iSplA1])*cHa;
  if(df != NULL){
    *df	=(g[iSplA1]-g[iSplA]
	  +((3.*SplAxx-1.)*d2g[iSplA1]-(3.*SplAXX-1.)*d2g[iSplA])*cHa)*crha;
  }
  if(d2f != NULL){
    *d2f	=SplAX*d2g[iSplA]+SplAx*d2g[iSplA1];
  }
  return(0);
}

int ESInitSplP()
{
  int i,i1;
  if(vP == NULL){
    vP	=(double*)malloc((2*ESNp1+24*ESNp1+ESNp1)*sizeof(double));
    d2fP	=vP+ESNp1;
    glP		=d2fP+ESNp1;
    WP		=glP	+ESNp1;
    VP		=WP	+6*ESNp1;
    U0P[0]	=VP	+9*ESNp1;
    U0P[1]	=U0P[0]	+ESNp1;
    U0P[2]	=U0P[1]	+ESNp1;
    U1P[0]	=U0P[2]	+ESNp1;
    U1P[1]	=U1P[0]	+ESNp1;
    U1P[2]	=U1P[1]	+ESNp1;
    U2P[0]	=U1P[2]	+ESNp1;
    U2P[1]	=U2P[0]	+ESNp1;
    U2P[2]	=U2P[1]	+ESNp1;
  }
  chp	=EZcr6*ESgt[1];
  crhp	=1./ESgt[1];
  chhp	=EZcr2*ESgt[1];
  ch2hp	=chhp*chp;
  cHp	=chp*ESgt[1];
  crHp	=1./ESNp;
  cx6prhp	=6.*crhp*crhp;
  cprhp		=crhp*crhp;
  cr6php	=EZcr6*ESgt[1]*ESgt[1];

  d2fP[0]	=1.;
  vP[0]	=0.;
  i1	=0;
  for(i=1; i < ESNp; i++){
    vP[i]	=-1./(4.+vP[i1]);
    d2fP[i]	=d2fP[i1]*vP[i];
    i1++;
  }
  i	=ESNp;
  d2fP[i]=d2fP[0];
  i1	=i-1;
  while(i > 0){
    d2fP[i1]	+=vP[i1]*d2fP[i];
    i--;
    i1--;
  }
  crd2fP	=1./(d2fP[1]+4.*d2fP[0]+d2fP[ESNp-1]);
  
  ESmemlev |=0x00800000;
  return(0);
}

int DeInitf2splP()
{
  if(vP != NULL){
    free(vP);
    vP	=NULL;
  }
  return(0);
}

int splP(double *g,double *d2g)
{
  int i,i1;
  double a1,a2;

  d2g[0]=(g[1]-2.*g[0]+g[ESNp-1])*crhp*crhp;
  vP[0]	=0;
  a2	=(g[1]-g[0])*cx6prhp;
  i1	=0;
  for(i=1; i < ESNp; i++){
    a1		=a2;
    a2		=(g[i+1]-g[i])*cx6prhp;
    vP[i]	=-1./(4.+vP[i1]);
    d2g[i]	=(d2g[i1]+a1-a2)*vP[i];
    i1++;
  }
  i	=ESNp;
  d2g[i]=d2g[0];
  i1	=i-1;
  while(i > 0){
    d2g[i1]	+=vP[i1]*d2g[i];
    i--;
    i1--;
  }
  i1	=ESNp-1;
  a1	=(2.*d2g[0]-d2g[1]-d2g[i1])*crd2fP;
  for(i=0; i < ESNp1; i++){
    d2g[i]	+=a1*d2fP[i]; 
  }
  return(0);
}

int f2splP(double *g,double *d2g,double *f,
	   double EZga0,double ga_1,double EZga2,double ga_3)
{
  static double Cga0=0.,Cga1=0.,Cga2=0.,Cga3=0.,EZga1,EZga3;
  static double w[6],rw[6],RHS[9];
  static double t11,t13,t22,t23,t31,t32;
  static double *pV,*pW;
  int i,ii;

  if(Cga0 != EZga0 || Cga1 != ga_1 || Cga2 != EZga2 || Cga3 !=ga_3){
    int indx[3];
    double ac[9];

    Cga0	=EZga0;
    Cga1	=ga_1;
    Cga2	=EZga2;
    Cga3	=ga_3;
    EZga1		=ga_1*cprhp;
    EZga3		=ga_3*cprhp;
    w[0]	=1.+EZga0+2.*EZga1;
    w[1]	=0.;
    w[2]	=2.;
    w[3]	=EZga2+2.*EZga3;
    w[4]	=4.*cr6php;
    w[5]	=0.;
    t11		=EZga1;
    t13		=1.;
    t22		=EZga3;
    t23		=-cr6php;
    t31		=1.;
    t32		=-cr6php;
    rw[0]	=w[0];
    rw[1]	=w[1];
    rw[2]	=w[2];
    rw[3]	=w[3];
    rw[4]	=w[4];
    rw[5]	=w[5];

    U0P[0][0]	=1.;
    U0P[1][0]	=0.;
    U0P[2][0]	=0.;
    U1P[0][0]	=0.;
    U1P[1][0]	=1.;
    U1P[2][0]	=0.;
    U2P[0][0]	=0.;
    U2P[1][0]	=0.;
    U2P[2][0]	=1.;

    pV		=VP;
    pW		=WP;
    ii	=0;
    for(i=1; i < ESNp; i++){
      EZinv3x3(rw,pW);
      pV[0]	=pW[0]*t11+pW[2]*t31;
      pV[1]	=pW[1]*t22+pW[2]*t32;
      pV[2]	=pW[0]*t13+pW[1]*t23;
      pV[3]	=pW[1]*t11+pW[4]*t31;
      pV[4]	=pW[3]*t22+pW[4]*t32;
      pV[5]	=pW[1]*t13+pW[3]*t23;
      pV[6]	=pW[2]*t11+pW[5]*t31;
      pV[7]	=pW[4]*t22+pW[5]*t32;
      pV[8]	=pW[2]*t13+pW[4]*t23;
      rw[0]	=w[0]-t11*pV[0]-t13*pV[6];
      rw[1]	=    -t11*pV[1]-t13*pV[7];
      rw[2]	=w[2]-t11*pV[2]-t13*pV[8];
      rw[3]	=w[3]-t22*pV[4]-t23*pV[7];
      rw[4]	=w[4]-t22*pV[5]-t23*pV[8];
      rw[5]	=    -t31*pV[2]-t32*pV[5];

      U0P[0][i]	=pV[0]*U0P[0][ii]+pV[1]*U0P[1][ii]+pV[2]*U0P[2][ii];
      U0P[1][i]	=pV[3]*U0P[0][ii]+pV[4]*U0P[1][ii]+pV[5]*U0P[2][ii];
      U0P[2][i]	=pV[6]*U0P[0][ii]+pV[7]*U0P[1][ii]+pV[8]*U0P[2][ii];

      U1P[0][i]	=pV[0]*U1P[0][ii]+pV[1]*U1P[1][ii]+pV[2]*U1P[2][ii];
      U1P[1][i]	=pV[3]*U1P[0][ii]+pV[4]*U1P[1][ii]+pV[5]*U1P[2][ii];
      U1P[2][i]	=pV[6]*U1P[0][ii]+pV[7]*U1P[1][ii]+pV[8]*U1P[2][ii];
      
      U2P[0][i]	=pV[0]*U2P[0][ii]+pV[1]*U2P[1][ii]+pV[2]*U2P[2][ii];
      U2P[1][i]	=pV[3]*U2P[0][ii]+pV[4]*U2P[1][ii]+pV[5]*U2P[2][ii];
      U2P[2][i]	=pV[6]*U2P[0][ii]+pV[7]*U2P[1][ii]+pV[8]*U2P[2][ii];
      pW	+=6;
      pV	+=9;
      ii++;
    }
    U0P[0][i]	=1.;
    U0P[1][i]	=0.;
    U0P[2][i]	=0.;
    U1P[0][i]	=0.;
    U1P[1][i]	=1.;
    U1P[2][i]	=0.;
    U2P[0][i]	=0.;
    U2P[1][i]	=0.;
    U2P[2][i]	=1.;
    while(i > 1){
      pV -=9;
      U0P[0][ii]	+=pV[0]*U0P[0][i]+pV[1]*U0P[1][i]+pV[2]*U0P[2][i];
      U0P[1][ii]	+=pV[3]*U0P[0][i]+pV[4]*U0P[1][i]+pV[5]*U0P[2][i];
      U0P[2][ii]	+=pV[6]*U0P[0][i]+pV[7]*U0P[1][i]+pV[8]*U0P[2][i];

      U1P[0][ii]	+=pV[0]*U1P[0][i]+pV[1]*U1P[1][i]+pV[2]*U1P[2][i];
      U1P[1][ii]	+=pV[3]*U1P[0][i]+pV[4]*U1P[1][i]+pV[5]*U1P[2][i];
      U1P[2][ii]	+=pV[6]*U1P[0][i]+pV[7]*U1P[1][i]+pV[8]*U1P[2][i];

      U2P[0][ii]	+=pV[0]*U2P[0][i]+pV[1]*U2P[1][i]+pV[2]*U2P[2][i];
      U2P[1][ii]	+=pV[3]*U2P[0][i]+pV[4]*U2P[1][i]+pV[5]*U2P[2][i];
      U2P[2][ii]	+=pV[6]*U2P[0][i]+pV[7]*U2P[1][i]+pV[8]*U2P[2][i];
      i--;
      ii--;
    }
    i	=1;
    ii	=ESNp-1;
    RHS[0]	=VP[0]*(U0P[0][ii]+U0P[0][i])+VP[1]*(U0P[1][ii]+U0P[1][i])
      +VP[2]*(U0P[2][ii]+U0P[2][i])-1.;
    RHS[1]	=VP[3]*(U0P[0][ii]+U0P[0][i])+VP[4]*(U0P[1][ii]+U0P[1][i])
      +VP[5]*(U0P[2][ii]+U0P[2][i]);
    RHS[2]	=VP[6]*(U0P[0][ii]+U0P[0][i])+VP[7]*(U0P[1][ii]+U0P[1][i])
      +VP[8]*(U0P[2][ii]+U0P[2][i]);

    RHS[3]	=VP[0]*(U1P[0][ii]+U1P[0][i])+VP[1]*(U1P[1][ii]+U1P[1][i])
      +VP[2]*(U1P[2][ii]+U1P[2][i]);
    RHS[4]	=VP[3]*(U1P[0][ii]+U1P[0][i])+VP[4]*(U1P[1][ii]+U1P[1][i])
      +VP[5]*(U1P[2][ii]+U1P[2][i])-1.;
    RHS[5]	=VP[6]*(U1P[0][ii]+U1P[0][i])+VP[7]*(U1P[1][ii]+U1P[1][i])
      +VP[8]*(U1P[2][ii]+U1P[2][i]);
      
    RHS[6]	=VP[0]*(U2P[0][ii]+U2P[0][i])+VP[1]*(U2P[1][ii]+U2P[1][i])
      +VP[2]*(U2P[2][ii]+U2P[2][i]);
    RHS[7]	=VP[3]*(U2P[0][ii]+U2P[0][i])+VP[4]*(U2P[1][ii]+U2P[1][i])
      +VP[5]*(U2P[2][ii]+U2P[2][i]);
    RHS[8]	=VP[6]*(U2P[0][ii]+U2P[0][i])+VP[7]*(U2P[1][ii]+U2P[1][i])
      +VP[8]*(U2P[2][ii]+U2P[2][i])-1.;

    ac[0]	=RHS[0];
    ac[1]	=RHS[3];
    ac[2]	=RHS[6];
    ac[3]	=RHS[1];
    ac[4]	=RHS[4];
    ac[5]	=RHS[7];
    ac[6]	=RHS[2];
    ac[7]	=RHS[5];
    ac[8]	=RHS[8];
    LUdcmp(ac,3,indx,RHS);
    RHS[0]	=1.;
    RHS[1]	=0.;
    RHS[2]	=0.;
    LUbksb(ac,3,indx,RHS);
    RHS[3]	=0.;
    RHS[4]	=1.;
    RHS[5]	=0.;
    LUbksb(ac,3,indx,RHS+3);
    RHS[6]	=0.;
    RHS[7]	=0.;
    RHS[8]	=1.;
    LUbksb(ac,3,indx,RHS+6);
    for(i=0; i < ESNp1; i++){
      ac[0]	=U0P[0][i];
      ac[1]	=U0P[1][i];
      ac[2]	=U0P[2][i];
      ac[3]	=U1P[0][i];
      ac[4]	=U1P[1][i];
      ac[5]	=U1P[2][i];
      ac[6]	=U2P[0][i];
      ac[7]	=U2P[1][i];
      ac[8]	=U2P[2][i];
      U0P[0][i]	=RHS[0]*ac[0]+RHS[1]*ac[3]+RHS[2]*ac[6];
      U0P[1][i]	=RHS[0]*ac[1]+RHS[1]*ac[4]+RHS[2]*ac[7];
      U0P[2][i]	=RHS[0]*ac[2]+RHS[1]*ac[5]+RHS[2]*ac[8];

      U1P[0][i]	=RHS[3]*ac[0]+RHS[4]*ac[3]+RHS[5]*ac[6];
      U1P[1][i]	=RHS[3]*ac[1]+RHS[4]*ac[4]+RHS[5]*ac[7];
      U1P[2][i]	=RHS[3]*ac[2]+RHS[4]*ac[5]+RHS[5]*ac[8];

      U2P[0][i]	=RHS[6]*ac[0]+RHS[7]*ac[3]+RHS[8]*ac[6];
      U2P[1][i]	=RHS[6]*ac[1]+RHS[7]*ac[4]+RHS[8]*ac[7];
      U2P[2][i]	=RHS[6]*ac[2]+RHS[7]*ac[5]+RHS[8]*ac[8];
    }
  }

  pW	=WP;
  pV	=VP;
  glP[0]	=pW[2]*f[0];
  d2g[0]	=(f[1]-2.*f[0]+f[ESNp-1])*cprhp;
  RHS[1]	=WP[1]*f[0]-d2g[0];
  RHS[2]	=WP[2]*f[0]-glP[0];
  g[0]		=f[0];
  RHS[0]	=WP[0]*f[0]-g[0];
  ii	=0;
  for(i=1; i < ESNp; i++){
    d2g[i]	=pW[1]*f[i]+pV[3]*g[ii]+pV[4]*d2g[ii]+pV[5]*glP[ii];
    glP[i]	=pW[2]*f[i]+pV[6]*g[ii]+pV[7]*d2g[ii]+pV[8]*glP[ii];
    g[i]	=pW[0]*f[i]+pV[0]*g[ii]+pV[1]*d2g[ii]+pV[2]*glP[ii];
    pW	+=6;
    pV	+=9;
    ii++;
  }
  g[i]	=g[0];
  d2g[i]=d2g[0];
  glP[i]	=glP[0];
  while(i > 1){
    pV -=9;
    g[ii]	+=pV[0]*g[i]+pV[1]*d2g[i]+pV[2]*glP[i];
    d2g[ii]	+=pV[3]*g[i]+pV[4]*d2g[i]+pV[5]*glP[i];
    glP[ii]	+=pV[6]*g[i]+pV[7]*d2g[i]+pV[8]*glP[i];
    i--;
    ii--;
  }

  i	=1;
  ii	=ESNp-1;
  RHS[0]+=VP[0]*(g[ii]+g[i])+VP[1]*(d2g[ii]+d2g[i])+VP[2]*(glP[ii]+glP[i]);
  RHS[1]+=VP[3]*(g[ii]+g[i])+VP[4]*(d2g[ii]+d2g[i])+VP[5]*(glP[ii]+glP[i]);
  RHS[2]+=VP[6]*(g[ii]+g[i])+VP[7]*(d2g[ii]+d2g[i])+VP[8]*(glP[ii]+glP[i]);
  for(i=0; i < ESNp1; i++){
    g[i]	-=RHS[0]*U0P[0][i]+RHS[1]*U1P[0][i]+RHS[2]*U2P[0][i];
    d2g[i]	-=RHS[0]*U0P[1][i]+RHS[1]*U1P[1][i]+RHS[2]*U2P[1][i];
  }
  return(0);
}

int ESSetSplP(double T)
{
  int iT;
  int i,j,kT1;
  
  if(T > EZc2gp){
    T	-=EZc2gp;
  }
  if(T < 0.){
    T	+=EZc2gp;
  }
  while(iSplP < ESNp && ESgt[iSplP+1] < T){
    iSplP++;
  }
  while(iSplP > 0 && ESgt[iSplP] > T){
    iSplP--;
  }
  iSplP1	=iSplP+1;
  SplPx		=(T-ESgt[iSplP])*crhp;
  SplPxx	=SplPx*SplPx;
  SplPX		=1.-SplPx;
  SplPXX	=SplPX*SplPX;
  return(0);
}

int splRP(double*f, double*df, double*g, double*d2g)
{
  *f	=SplPX*g[iSplP]+SplPx*g[iSplP1]
    +(SplPX*(SplPXX-1.)*d2g[iSplP]+SplPx*(SplPxx-1.)*d2g[iSplP1])*cHp;
  if(df != NULL){
    *df	=(g[iSplP1]-g[iSplP]
	  +((3.*SplPxx-1.)*d2g[iSplP1]-(3.*SplPXX-1.)*d2g[iSplP])*cHp)*crhp;
  }
  return(0);
}

/* Splines with respect to pa=a*a */
static int iSplAA,iSplAA1,iSplAAold;
static double xAA,xxAA,crxAA,cHxAA,SplAAx,SplAAxx,SplAAX,SplAAXX;

int ESFirstSetSplXAA()
{
  static int i,ii,i1;
  static double t1,t0,w0,w1;

  iSplAA	=0;
  iSplAA1	=1;
  crxAA		=ESpa[1]-ESpa[0];
  cHxAA		=EZcr6*crxAA*crxAA;
  crxAA		=1./crxAA;
  iSplAAold	=0;
  xAA	=0.;
  xxAA	=0.;

  w0	=ESpa[1]*EZcr3;
  t0	=EZcr2*w0;
  w1	=(ESpa[2]-ESpa[1])*EZcr3;
  t1	=EZcr2*w1;
  vXAA0[0]	=2.;
  vXAA0[1]	=(t0-t1)/(w0+w1+t0*vXAA0[0]);
  vXAA1[0]	=2.;
  t0	+=w0;
  vXAA1[1]	=(t0-t1)/(w1+t0+t0*vXAA1[0]);
  for(i=2,ii=3,i1=1; i < ESNa; i++,ii++,i1++){
    w0		=w1;
    w1		=(ESpa[ii]-ESpa[i])*EZcr3;
    t0		=t1;
    t1		=EZcr2*w1;
    vXAA0[i]	=-t1/(w0+w1+t0*vXAA0[i1]);
    vXAA1[i]	=-t1/(w0+w1+t0*vXAA1[i1]);
  }

  return(0);
}

int ESSetSplAA(double A)
{
  if(A > ESsa[ESNa]){
    A	=ESsa[ESNa];
  }
  if(A < ESa0){
    A	=ESa0;
  }
  xAA	=A;
  xxAA	=A*A;
  while(iSplAA < ESNa && ESsa[iSplAA+1] < A){
    iSplAA++;
  }
  while(iSplAA > 0 && ESsa[iSplAA] > A){
    iSplAA--;
  }
  if(iSplAAold	!= iSplAA){
    iSplAAold	=iSplAA;
    iSplAA1	=iSplAA+1;
    crxAA	=ESpa[iSplAA1]-ESpa[iSplAA];
    cHxAA	=EZcr6*crxAA*crxAA;
    crxAA	=1./crxAA;
  }
  SplAAx	=(xxAA-ESpa[iSplAA])*crxAA;
  SplAAxx	=SplAAx*SplAAx;

  SplAAX	=1.-SplAAx;
  SplAAXX	=SplAAX*SplAAX;
  return(0);
}

int splRAA0(double*f, double*df, double*d0, double*d2)
{
  *f	=SplAAX*d0[iSplAA]+SplAAx*d0[iSplAA1]
    +(SplAAX*(SplAAXX-1.)*d2[iSplAA]+SplAAx*(SplAAxx-1.)*d2[iSplAA1])*cHxAA;
  if(df != NULL){
    *df	=(d0[iSplAA1]-d0[iSplAA]+
	  ((3.*SplAAxx-1.)*d2[iSplAA1]-(3.*SplAAXX-1.)*d2[iSplAA])*cHxAA
	  )*crxAA*2.*xAA;
  }
  return(0);
}

int splRAA1(double*f, double*df, double*d0, double*d2)
{
  *f	=SplAAX*d0[iSplAA]+SplAAx*d0[iSplAA1]
    +(SplAAX*(SplAAXX-1.)*d2[iSplAA]+SplAAx*(SplAAxx-1.)*d2[iSplAA1])*cHxAA;
  if(df != NULL){
    *df	=(d0[iSplAA1]-d0[iSplAA]+
	  ((3.*SplAAxx-1.)*d2[iSplAA1]-(3.*SplAAXX-1.)*d2[iSplAA])*cHxAA
	  )*crxAA*2.*xxAA+(*f);
  }
  *f	*=xAA;
  return(0);
}

int splRAA20(double*f, double*df, double*d2f, double*d0, double*d2)
{
  double ds;
  *f	=SplAAX*d0[iSplAA]+SplAAx*d0[iSplAA1]
    +(SplAAX*(SplAAXX-1.)*d2[iSplAA]+SplAAx*(SplAAxx-1.)*d2[iSplAA1])*cHxAA;
  ds	=(d0[iSplAA1]-d0[iSplAA]
	  +((3.*SplAAxx-1.)*d2[iSplAA1]-(3.*SplAAXX-1.)*d2[iSplAA])*cHxAA
	  )*crxAA;
  if(df != NULL){
    (*df)	=ds*2.*xAA;
  }
  if(d2f != NULL){
    *d2f	=4.*xxAA*(SplAAX*d2[iSplAA]+SplAAx*d2[iSplAA1])+2.*ds;
  }
  return(0);
}

int splRAA21(double*f, double*df, double*d2f, double*d0, double*d2)
{
  double ds;
  *f	=SplAAX*d0[iSplAA]+SplAAx*d0[iSplAA1]
    +(SplAAX*(SplAAXX-1.)*d2[iSplAA]+SplAAx*(SplAAxx-1.)*d2[iSplAA1])*cHxAA;
  ds	=(d0[iSplAA1]-d0[iSplAA]
	  +((3.*SplAAxx-1.)*d2[iSplAA1]-(3.*SplAAXX-1.)*d2[iSplAA])*cHxAA
	  )*crxAA;
  if(df != NULL){
    *df	=2.*xxAA*ds+(*f);
  }
  if(d2f != NULL){
    *d2f	=xAA*(4.*xxAA*(SplAAX*d2[iSplAA]+SplAAx*d2[iSplAA1])+6.*ds);
  }
  *f	*=xAA;
  return(0);
}

int splXAA0(double *g, double *d2g, double *df1)
{
  static int i,ii,i1;
  static double h,t1,t0,r0,r1,w0,w1;
  double *pV;

  w0	=ESpa[1]*EZcr3;
  t0	=EZcr2*w0;
  h	=ESpa[2]-ESpa[1];
  w1	=h*EZcr3;
  t1	=EZcr2*w1;
  r0		=(g[1]-g[0])/ESpa[1];
  r1		=(g[2]-g[1])/h;
  pV		=vXAA0;
  d2g[0]	=0.;
  d2g[1]	=(r1-r0-t0*d2g[0])/(w0+w1+t0*pV[0]);

  for(i=2,ii=3,i1=1; i <ESNa ;i++,ii++,i1++){
    h		=ESpa[ii]-ESpa[i];
    w0		=w1;
    w1		=h*EZcr3;
    t0		=t1;
    t1		=EZcr2*w1;
    r0		=r1;
    r1		=(g[ii]-g[i])/h;
    d2g[i]	=(r1-r0-t0*d2g[i1])/(w0+w1+t0*pV[i1]);
  }

  ii	=ESNa;
  i	=ii-1;
  i1	=i-1;
  if(df1 == NULL){
    t1		=pV[i1]-2.;
    d2g[ii]	=-(d2g[i1]+t1*d2g[i])/(t1*pV[i]+1.);
  }
  else{
    d2g[ii]	=(EZcr2*(*df1)/ESsa[ii]-r1-t1*d2g[i])/(t1*pV[i]+w1);
  }
  while(ii > 0){
    d2g[i]	+=pV[i]*d2g[ii];
    i--;
    ii--;
  }
  d2g[0]	-=d2g[2];
  return(0);
}

int splXAA1(double *g, double *d2g, double *f, double df0, double *df1)
{
  static int i,ii,i1;
  static double h,t1,t0,r0,r1,w0,w1;
  double *pV;

  g[0]	=df0;
  w0	=ESpa[1]*EZcr3;
  t0	=EZcr2*w0;
  h	=ESpa[2]-ESpa[1];
  w1	=h*EZcr3;
  t1	=EZcr2*w1;
  g[1]	=f[1]/ESsa[1];
  g[2]	=f[2]/ESsa[2];
  r0	=(g[1]-g[0])/ESpa[1];
  r1	=(g[2]-g[1])/h;
  pV	=vXAA0;
  d2g[0]	=0.;
  d2g[1]	=(r1-r0-t0*d2g[0])/(w0+w1+t0*pV[0]);
  for(i=2,ii=3,i1=1; i < ESNa; i++,ii++,i1++){
    h		=ESpa[ii]-ESpa[i];
    w0		=w1;
    w1		=h*EZcr3;
    t0		=t1;
    t1		=EZcr2*w1;
    r0		=r1;
    g[ii]	=f[ii]/ESsa[ii];
    r1		=(g[ii]-g[i])/h;
    d2g[i]	=(r1-r0-t0*d2g[i1])/(w0+w1+t0*pV[i1]);
  }
  ii	=ESNa;
  i	=ii-1;
  i1	=i-1;
  if(df1 == NULL){
    t1		=pV[i1]-2.;
    d2g[ii]	=-(d2g[i1]+t1*d2g[i])/(t1*pV[i]+1.);
  }
  else{
    d2g[ii]	=(EZcr2*((*df1)-g[ii])/ESpa[ii]-r1-t1*d2g[i])/(t1*pV[i]+w1);
  }
  while(ii > 0){
    d2g[i]	+=pV[i]*d2g[ii];
    i--;
    ii--;
  }
  d2g[0]	-=d2g[2];
  return(0);
}

#ifndef mzl_2Dcore

static int nrMap,nzMap,nrMap1,nzMap1,nMap,nMap1;
static int *iAMap=NULL,*k00Map=NULL,irMap,izMap,k00,k01,k10,k11;
static double *f0Map,*fxMap,*fyMap,*f2Map;
static double *F0Map,*FxMap,*FyMap,*F2Map;
static double rBox[2],zBox[2],drMap,dzMap,rdrMap,rdzMap;
static double *gY;
static double X[4],dX[4],Y[4],dY[4];

int ESDeInitMapSpl()
{
  if(iAMap != NULL){
    free(f0Map);
    free(iAMap);
    iAMap	=NULL;
  }
  return(0);
}

/*************/
extern double *ESaR0;
extern double gy00,gy01,gy02,gy03,gy10,gy11,gy12,gy13,
gy20,gy21,gy22,gy23,gy30,gy31,gy32,gy33;

int ESSetgY()
{
  int i,ii;
  double g0,g1,d2g0,d2g1;

#ifdef H
  for(i=0; i < ESNa1; i++){
    ESdgY[i]	=gy00+ESsa[i]*(gy01+ESsa[i]*(gy02+ESsa[i]*gy03));
    ESdgY2a[i]	=2.*gy02+6.*ESsa[i]*gy03;
  }
#endif

  gY[0]	=0.;
  g1	=ESdgY[0];
  d2g1	=ESdgY2a[0];
  for(i=1,ii=0; i < ESNa1; i++,ii++){
    g0		=g1;
    d2g0	=d2g1;
    g1		=ESdgY[i];
    d2g1	=ESdgY2a[i];
    gY[i]	=gY[ii]+chha*ESsa[ii]*(g0+g1-ch2ha*(d2g0+d2g1))
      +cHa*(g0+2.*g1-cHa*(0.7*d2g0+0.8*d2g1));
  }
  return(0);
}

int ESSetSpl4gY(double r,double z)
{
  int i;
  static double h,rh,w,rw;
  double t,tt,T,TT,tT;

  if(r < rBox[0] || r > rBox[1] || z < zBox[0] || z > zBox[1]) return(1); 
  t	=(r-rBox[0])*rdrMap;
  i	=t;
  if(i == nrMap) i--;
  t	-=i;
  T	=1.-t;
  tt	=t*t;
  TT	=T*T;
  tT	=2.*t*T;

  X[0]	=TT*(3.-2.*T);
  X[1]	=tt*(3.-2.*t);
  X[2]	=TT*t;
  X[3]	=-tt*T;
  dX[1]	=3.*tT;
  dX[0]	=-dX[1];
  dX[2]	=TT-tT;
  dX[3]	=tt-tT;

  irMap	=i;
  t	=(z-zBox[0])*rdzMap;
  i	=t;
  if(i == nzMap) i--;
  t	-=i;
  T	=1.-t;
  tt	=t*t;
  TT	=T*T;
  tT	=2.*t*T;

  Y[0]	=TT*(3.-2.*T);
  Y[1]	=tt*(3.-2.*t);
  Y[2]	=TT*t;
  Y[3]	=-tt*T;
  dY[1]	=3.*tT;
  dY[0]	=-dY[1];
  dY[2]	=TT-tT;
  dY[3]	=tt-tT;
  izMap	=i;

  i	=nrMap*izMap+irMap;
  if(iAMap[i] == 0) return(1);

  k00	=k00Map[i];
  k01	=k00+1;
  k10	=iAMap[i] == 0x0f ? k00+nrMap1 : k00+2;
  k11	=k10+1;
  return(0);
}

int SetSplX(double dr)
{
  double t,tt,T,TT,tT;

  t	=dr*rdrMap;
  T	=1.-t;
  tt	=t*t;
  TT	=T*T;
  tT	=2.*t*T;
  X[0]	=TT*(3.-2.*T);
  X[1]	=tt*(3.-2.*t);
  X[2]	=TT*t;
  X[3]	=-tt*T;
  dX[1]	=3.*tT;
  dX[0]	=-dX[1];
  dX[2]	=TT-tT;
  dX[3]	=tt-tT;
  return(0);
}

int SetSplY(double dz)
{
  double t,tt,T,TT,tT;

  t	=dz*rdzMap;
  T	=1.-t;
  tt	=t*t;
  TT	=T*T;
  tT	=2.*t*T;
  Y[0]	=TT*(3.-2.*T);
  Y[1]	=tt*(3.-2.*t);
  Y[2]	=TT*t;
  Y[3]	=-tt*T;
  dY[1]	=3.*tT;
  dY[0]	=-dY[1];
  dY[2]	=TT-tT;
  dY[3]	=tt-tT;
  return(0);
}

int ESGetgY(double *gy,double *gYr,double *gYz,double r, double z)
{
  int i,j,k;
  double dx,dy,d;
  double a,gt;
  double ra,raa,rp,rap,rpp,za,zaa,zp,zap,zpp,gYa,gYaa;

#ifdef H
  r	=r-ESaR0[0];
  *gy  =gy00+r*(gy01+r*(gy02+r*gy03))
    +z*((gy10+r*(gy11+r*(gy12+r*gy13)))
	+z*((gy20+r*(gy21+r*(gy22+r*gy23)))
	    +z*(gy30+r*(gy31+r*(gy32+r*gy33)))));
  *gYr  =gy01+r*(2.*gy02+3.*r*gy03)
    +z*(gy11+r*(2.*gy12+3.*r*gy13)
	+z*(gy21+r*(2.*gy22+3.*r*gy23)
	    +z*(gy31+r*(2.*gy32+3.*r*gy33))));
  *gYz  =gy10+r*(gy11+r*(gy12+r*gy13))
    +z*(2.*(gy20+r*(gy21+r*(gy22+r*gy23)))
	+3.*z*(gy30+r*(gy31+r*(gy32+r*gy33))));

  EZout("sdddd","Flx gY=",*gy,*gYr,*gYz,
      gy11+r*(2.*gy12+3.*r*gy13)
      +z*(2.*(gy21+r*(2.*gy22+3.*r*gy23))
	  +3.*z*(gy31+r*(2.*gy32+3.*r*gy33)))
      );
#endif

  if(r < rBox[0] || r > rBox[1] || z < zBox[0] || z > zBox[1]){
    *gy	=0.;
    *gYr=0.;
    *gYz=0.;
    return(1); 
  }

  k	=0;
  d	=1e+20;
  for(i=0; i < ESNa1; i++){
    k	=ESNp1*i;
    for(j=0; j < ESNp; j++){
      dx	=r-ESsr[k];
      dy	=z-ESsz[k];
      dx	=dx*dx+dy*dy;
      if(d > dx){
	d	=dx;
	a	=ESsa[i];
	gt	=ESgt[j];
      }
      k++;
    }
  }
  ES2DMapLab2Flx(&a,&gt,r,z);
  ESGetFullMetrics(&ra,&rp,&raa,&rap,&rpp,&za,&zp,&zaa,&zap,&zpp,a,gt);
  splREA(gy,&gYa,&gYaa,gY,ESdgY,ESdgY2a);

  gYaa	=a*gYaa+gYa;
  gYa	*=a;
  d	=ra*zp-za*rp;
  if(d != 0.)	d=1./d;
  *gYr	=zp*gYa*d;
  *gYz	=-rp*gYa*d;

  ra	*=d;
  rp	*=d;
  za	*=d;
  zp	*=d;
  dx	=raa*zp+ra*zap-zaa*rp-za*rap;
  dy	=rap*zp+ra*zpp-zap*rp-za*rpp;
  dx	=((dx*rp-dy*ra)*zp+(ra*zpp-rp*zap)*d)*gYa-rp*zp*gYaa;
  return(0);
}

int ESGetFullFgY(double *F,double *dgy,double a, double gt)
{
  int i,j,k;
  double dx,dy,d;
  double ra,raa,rp,rap,rpp,za,zaa,zp,zap,zpp,gYa,gYaa;

#ifdef H
  r	=r-ESaR0[0];
  *gy  =gy00+r*(gy01+r*(gy02+r*gy03))
    +z*((gy10+r*(gy11+r*(gy12+r*gy13)))
	+z*((gy20+r*(gy21+r*(gy22+r*gy23)))
	    +z*(gy30+r*(gy31+r*(gy32+r*gy33)))));
  *gYr  =gy01+r*(2.*gy02+3.*r*gy03)
    +z*(gy11+r*(2.*gy12+3.*r*gy13)
	+z*(gy21+r*(2.*gy22+3.*r*gy23)
	    +z*(gy31+r*(2.*gy32+3.*r*gy33))));
  *gYz  =gy10+r*(gy11+r*(gy12+r*gy13))
    +z*(2.*(gy20+r*(gy21+r*(gy22+r*gy23)))
	+3.*z*(gy30+r*(gy31+r*(gy32+r*gy33))));

  EZout("sdddd","Flx gY=",*gy,*gYr,*gYz,
      gy11+r*(2.*gy12+3.*r*gy13)
      +z*(2.*(gy21+r*(2.*gy22+3.*r*gy23))
	  +3.*z*(gy31+r*(2.*gy32+3.*r*gy33)))
      );
#endif
  ESGetFullMetrics(&ra,&rp,&raa,&rap,&rpp,&za,&zp,&zaa,&zap,&zpp,a,gt);
  splREA(dgy,&gYa,&gYaa,gY,ESdgY,ESdgY2a);
  gYaa	=a*gYaa+gYa;
  gYa	*=a;
  d	=ra*zp-za*rp;
  if(d != 0.)	d=1./d;
  ra	*=d;
  rp	*=d;
  za	*=d;
  zp	*=d;
  dgy[1]= zp*gYa;
  dgy[2]=-rp*gYa;
  dx	=raa*zp+ra*zap-zaa*rp-za*rap;
  dy	=rap*zp+ra*zpp-zap*rp-za*rpp;
  dx	=((dx*rp-dy*ra)*zp+(ra*zpp-rp*zap)*d);
  dy	=rp*zp;
  dgy[3]=dx*gYa-dy*gYaa;

  splRA2(F,F+1,F+2,ESFF,ESFF2a);
  F[0]	=sqrt(F[0]);
  F[1]	*=0.5/F[0];	
  F[2]	=(0.5*F[2]-F[1]*F[1])/F[0];	

  F[3]	=dx*F[1]-dy*F[2];
  F[2]	=-rp*F[1];
  F[1]	*=zp;
  return(0);
}

int ESSpl2gY(double *gy,double *gYr,double *gYz)
{
  *gy=(Y[0]*(X[0]*f0Map[k00]+X[1]*f0Map[k01]+X[2]*fxMap[k00]+X[3]*fxMap[k01])
       +Y[1]*(X[0]*f0Map[k10]+X[1]*f0Map[k11]+X[2]*fxMap[k10]+X[3]*fxMap[k11])
       +Y[2]*(X[0]*fyMap[k00]+X[1]*fyMap[k01]+X[2]*f2Map[k00]+X[3]*f2Map[k01])
       +Y[3]*(X[0]*fyMap[k10]+X[1]*fyMap[k11]+X[2]*f2Map[k10]+X[3]*f2Map[k11])
       );
  *gYr=rdrMap*
   ( Y[0]*(dX[0]*f0Map[k00]+dX[1]*f0Map[k01]+dX[2]*fxMap[k00]+dX[3]*fxMap[k01])
    +Y[1]*(dX[0]*f0Map[k10]+dX[1]*f0Map[k11]+dX[2]*fxMap[k10]+dX[3]*fxMap[k11])
    +Y[2]*(dX[0]*fyMap[k00]+dX[1]*fyMap[k01]+dX[2]*f2Map[k00]+dX[3]*f2Map[k01])
    +Y[3]*(dX[0]*fyMap[k10]+dX[1]*fyMap[k11]+dX[2]*f2Map[k10]+dX[3]*f2Map[k11])
    );
  *gYz=rdzMap*
    ( dY[0]*(X[0]*f0Map[k00]+X[1]*f0Map[k01]+X[2]*fxMap[k00]+X[3]*fxMap[k01])
     +dY[1]*(X[0]*f0Map[k10]+X[1]*f0Map[k11]+X[2]*fxMap[k10]+X[3]*fxMap[k11])
     +dY[2]*(X[0]*fyMap[k00]+X[1]*fyMap[k01]+X[2]*f2Map[k00]+X[3]*f2Map[k01])
     +dY[3]*(X[0]*fyMap[k10]+X[1]*fyMap[k11]+X[2]*f2Map[k10]+X[3]*f2Map[k11])
     );
#ifdef H
  EZout("sidddd","Spl k00",k00,f0Map[k00],fxMap[k00]*rdrMap,fyMap[k00]*rdzMap
      ,f2Map[k00]*rdrMap*rdzMap);
  EZout("sidddd","Spl k01",k01,f0Map[k01],fxMap[k01]*rdrMap,fyMap[k01]*rdzMap
      ,f2Map[k01]*rdrMap*rdzMap);
  EZout("sidddd","Spl k10",k10,f0Map[k10],fxMap[k10]*rdrMap,fyMap[k10]*rdzMap
      ,f2Map[k10]*rdrMap*rdzMap);
  EZout("sidddd","Spl k11",k11,f0Map[k11],fxMap[k11]*rdrMap,fyMap[k11]*rdzMap
      ,f2Map[k11]*rdrMap*rdzMap);

  EZout("sdddd","Spl gy=",*gy,*gYr,*gYz,
   rdrMap*rdzMap*
     ( dY[0]*(dX[0]*f0Map[k00]+dX[1]*f0Map[k01]+dX[2]*fxMap[k00]+dX[3]*fxMap[k01])
      +dY[1]*(dX[0]*f0Map[k10]+dX[1]*f0Map[k11]+dX[2]*fxMap[k10]+dX[3]*fxMap[k11])
      +dY[2]*(dX[0]*fyMap[k00]+dX[1]*fyMap[k01]+dX[2]*f2Map[k00]+dX[3]*f2Map[k01])
      +dY[3]*(dX[0]*fyMap[k10]+dX[1]*fyMap[k11]+dX[2]*f2Map[k10]+dX[3]*f2Map[k11])
      )
      );
#endif
  return(0);
}

int ESSpl2F(double *F,double *Fr,double *Fz)
{
  *F=(Y[0]*(X[0]*f0Map[k00]+X[1]*f0Map[k01]+X[2]*fxMap[k00]+X[3]*fxMap[k01])
       +Y[1]*(X[0]*f0Map[k10]+X[1]*f0Map[k11]+X[2]*fxMap[k10]+X[3]*fxMap[k11])
       +Y[2]*(X[0]*fyMap[k00]+X[1]*fyMap[k01]+X[2]*f2Map[k00]+X[3]*f2Map[k01])
       +Y[3]*(X[0]*fyMap[k10]+X[1]*fyMap[k11]+X[2]*f2Map[k10]+X[3]*f2Map[k11])
       );
  *Fr=rdrMap*
   ( Y[0]*(dX[0]*f0Map[k00]+dX[1]*f0Map[k01]+dX[2]*fxMap[k00]+dX[3]*fxMap[k01])
    +Y[1]*(dX[0]*f0Map[k10]+dX[1]*f0Map[k11]+dX[2]*fxMap[k10]+dX[3]*fxMap[k11])
    +Y[2]*(dX[0]*fyMap[k00]+dX[1]*fyMap[k01]+dX[2]*f2Map[k00]+dX[3]*f2Map[k01])
    +Y[3]*(dX[0]*fyMap[k10]+dX[1]*fyMap[k11]+dX[2]*f2Map[k10]+dX[3]*f2Map[k11])
    );
  *Fz=rdzMap*
    ( dY[0]*(X[0]*f0Map[k00]+X[1]*f0Map[k01]+X[2]*fxMap[k00]+X[3]*fxMap[k01])
     +dY[1]*(X[0]*f0Map[k10]+X[1]*f0Map[k11]+X[2]*fxMap[k10]+X[3]*fxMap[k11])
     +dY[2]*(X[0]*fyMap[k00]+X[1]*fyMap[k01]+X[2]*f2Map[k00]+X[3]*f2Map[k01])
     +dY[3]*(X[0]*fyMap[k10]+X[1]*fyMap[k11]+X[2]*f2Map[k10]+X[3]*f2Map[k11])
     );
  return(0);
}

int ESgY2Spl(int nr, int nz)
{
  int m;
  int i,ii,j,jj,k,kk,n,n1,errFl;
  double r,z,r0,EZz0,r1,z1,dr,dz;
  double ra,raa,rp,rap,rpp,za,zaa,zp,zap,zpp,gYa,gYaa;
  double a,gt,d,dx,dy,f0,f1,fx,fy;
  double t,tt,T,TT,tT;
  double b1,b2,b3,a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double f00,f01,f10,f11,fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11
    ,f200,f210,f201,f211;
  double F00,F01,F10,F11,Fx00,Fx10,Fx01,Fx11,Fy00,Fy10,Fy01,Fy11
    ,F200,F210,F201,F211;
  int inA[4];
  double A[16],B[4],F[4];
  double err;
  double x;
  
  nrMap	=nr;
  nzMap	=nz;
  nrMap1=nrMap+1;
  nzMap1=nzMap+1;

  ESSetMapBox();
  ESGetMapBox(rBox,zBox);
  nMap	=nrMap*nzMap;
  if(iAMap != NULL){
    free(f0Map);
    free(iAMap);
  }
  iAMap	=(int*)malloc(2*nMap*sizeof(int));
  k00Map=iAMap+nMap;
  ESSetInsideInd(iAMap,rBox,zBox,nrMap,nzMap);
  k	=0;
  for(i=0; i < nMap; i++){
    if(iAMap[i] != 0 && iAMap[i] != 0x0f) k+=4; 
    k00Map[i]	=0;
  }
  k	+=nrMap1*nzMap1;
  nMap1	=k;
  f0Map	=(double*)calloc(8*k+ESNa1,sizeof(double));
  fxMap	=f0Map+k;
  fyMap	=fxMap+k;
  f2Map	=fyMap+k;
  F0Map	=f2Map+k;
  FxMap	=F0Map+k;
  FyMap	=FxMap+k;
  F2Map	=FyMap+k;
  gY	=F2Map+k;
  ESmemlev |=0x10000000;
  ESSetgY();

  dr	=(rBox[1]-rBox[0])/nrMap;
  dz	=(zBox[1]-zBox[0])/nzMap;
  drMap	=dr;
  dzMap	=dz;
  rdrMap	=1./dr;
  rdzMap	=1./dz;

  m	=0;
  k	=0;
  errFl=0;
  err	=0.;
  for(i=0; i < nzMap1; i++){
    z	=zBox[0]+dz*i;
    ES1DMapZ2gtP(&gt,&r0,&r1,z);
    for(j=0; j < nrMap1; j++){
      r	=rBox[0]+dr*j;
      if(r0 <= r && r <= r1){
	d	=1e+20;
	for(ii=0; ii < ESNa1; ii++){
	  kk	=ESNp1*ii;
	  for(jj=0; jj < ESNp; jj++){
	    dx	=r-ESsr[kk];
	    dy	=z-ESsz[kk];
	    dx	=dx*dx+dy*dy;
	    if(d > dx){
	      d	=dx;
	      a	=ESsa[ii];
	      gt=ESgt[jj];
	    }
	    kk++;
	  }
	}
	if(ES2DMapLab2Flx(&a,&gt,r,z)) putchar('\a');
	kk	=0;
	k	=nrMap*i+j;
	if(i < nzMap){
	  if(j < nrMap){
	    if((iAMap[k]&0x01) == 0) kk=1;
	  }
	  else{
	    k--;
	    if((iAMap[k]&0x02) == 0) kk =2;
	  }
	}
	else{
	  k	-=nrMap;
	  if(j < nrMap){
	    if((iAMap[k]&0x08) == 0) kk =4;
	  }
	  else{
	    k--;
	    if((iAMap[k]&0x04) == 0) kk=3;
	  }
	}
	if(kk) EZout("I<I I<I siI",i,nzMap1,j,nrMap1,"iA",kk-1,iAMap[k]);
#ifdef H
	ES2DMapFlx2Lab(&dx,&dy,&ra,&za,&rp,&zp,a,gt);
	dx	-=r;
	dy	-=z;
	dx	=fabs(dx)+fabs(dy);
	if(err < dx) err=dx;
#endif
	ESGetFullFgY(F,B,a,gt);
	f0Map[m]	=B[0];
	fxMap[m]	=B[1]*drMap;
	fyMap[m]	=B[2]*dzMap;
	f2Map[m]	=B[3]*drMap*dzMap;
	F0Map[m]	=F[0];
	FxMap[m]	=F[1]*drMap;
	FyMap[m]	=F[2]*dzMap;
	F2Map[m]	=F[3]*drMap*dzMap;
      }
      m++;
    }
  }
#ifdef H
  EZout("sd","Accuracy in Mapping=",err);
#endif  
  k00	=nrMap1*nzMap1;
  k	=0;
  errFl=0;
  for(i=0; i < nzMap; i++){
    EZz0	=zBox[0]+dz*i;
    z1	=zBox[0]+dz*(i+1);
    for(j=0; j < nrMap; j++){
      r0	=rBox[0]+dr*j;
      r1	=rBox[0]+dr*(j+1);
      m	=nrMap1*i+j;
      k00Map[k]=m;
      if(iAMap[k] != 0 && iAMap[k] != 15){
	a	=1.;
	k01	=k00+1;
	k10	=k00+2;
	k11	=k10+1;
	if((iAMap[k]&0x01)){
	  f0Map[k00]	=f0Map[m];
	  fxMap[k00]	=fxMap[m];
	  fyMap[k00]	=fyMap[m];
	  f2Map[k00]	=f2Map[m];
	  f00	=f0Map[k00];
	  fx00	=fxMap[k00];
	  fy00	=fyMap[k00];
	  f200	=f2Map[k00];

	  F0Map[k00]	=F0Map[m];
	  FxMap[k00]	=FxMap[m];
	  FyMap[k00]	=FyMap[m];
	  F2Map[k00]	=F2Map[m];
	  F00	=F0Map[k00];
	  Fx00	=FxMap[k00];
	  Fy00	=FyMap[k00];
	  F200	=F2Map[k00];
	}
	m++;
	if((iAMap[k]&0x02)){
	  f0Map[k01]=f0Map[m];
	  fxMap[k01]=fxMap[m];
	  fyMap[k01]=fyMap[m];
	  f2Map[k01]=f2Map[m];
	  f01	=f0Map[k01];
	  fx01	=fxMap[k01];
	  fy01	=fyMap[k01];
	  f201	=f2Map[k01];

	  F0Map[k01]=F0Map[m];
	  FxMap[k01]=FxMap[m];
	  FyMap[k01]=FyMap[m];
	  F2Map[k01]=F2Map[m];
	  F01	=F0Map[k01];
	  Fx01	=FxMap[k01];
	  Fy01	=FyMap[k01];
	  F201	=F2Map[k01];
	}
	m	+=nrMap1;
	if((iAMap[k]&0x04)){
	  f0Map[k11]=f0Map[m];
	  fxMap[k11]=fxMap[m];
	  fyMap[k11]=fyMap[m];
	  f2Map[k11]=f2Map[m];
	  f11	=f0Map[k11];
	  fx11	=fxMap[k11];
	  fy11	=fyMap[k11];
	  f211	=f2Map[k11];

	  F0Map[k11]=F0Map[m];
	  FxMap[k11]=FxMap[m];
	  FyMap[k11]=FyMap[m];
	  F2Map[k11]=F2Map[m];
	  F11	=F0Map[k11];
	  Fx11	=FxMap[k11];
	  Fy11	=FyMap[k11];
	  F211	=F2Map[k11];
	}
	m--;
	if((iAMap[k]&0x08)){
	  f0Map[k10]	=f0Map[m];
	  fxMap[k10]	=fxMap[m];
	  fyMap[k10]	=fyMap[m];
	  f2Map[k10]	=f2Map[m];
	  f10	=f0Map[k10];
	  fx10	=fxMap[k10];
	  fy10	=fyMap[k10];
	  f210	=f2Map[k10];

	  F0Map[k10]	=F0Map[m];
	  FxMap[k10]	=FxMap[m];
	  FyMap[k10]	=FyMap[m];
	  F2Map[k10]	=F2Map[m];
	  F10	=F0Map[k10];
	  Fx10	=FxMap[k10];
	  Fy10	=FyMap[k10];
	  F210	=F2Map[k10];
	}
	switch(iAMap[k]){
	case 1:
	  errFl	|=ES1DMapZ2gt(&gt,&b2,&r,EZz0);
	  SetSplX(r-r0);
	  a11	=X[1];
	  a12	=X[3];
	  a21	=dX[1];
	  a22	=dX[3];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  B[0]	-= X[0]*f00 + X[2]*fx00;
	  B[1]	-=dX[0]*f00 +dX[2]*fx00;
	  B[2]	-= X[0]*fy00+ X[2]*f200;
	  B[3]	-=dX[0]*fy00+dX[2]*f200;
	  f01	=a22*B[0]-a12*B[1];
	  fx01	=a11*B[1]-a21*B[0];
	  fy01	=a22*B[2]-a12*B[3];
	  f201	=a11*B[3]-a21*B[2];
	  f0Map[k01]	= f01;
	  fxMap[k01]	=fx01;
	  fyMap[k01]	=fy01;
	  f2Map[k01]	=f201;

	  F[0]	-= X[0]*F00 + X[2]*Fx00;
	  F[1]	-=dX[0]*F00 +dX[2]*Fx00;
	  F[2]	-= X[0]*Fy00+ X[2]*F200;
	  F[3]	-=dX[0]*Fy00+dX[2]*F200;
	  F01	=a22*F[0]-a12*F[1];
	  Fx01	=a11*F[1]-a21*F[0];
	  Fy01	=a22*F[2]-a12*F[3];
	  F201	=a11*F[3]-a21*F[2];
	  F0Map[k01]	= F01;
	  FxMap[k01]	=Fx01;
	  FyMap[k01]	=Fy01;
	  F2Map[k01]	=F201;

	  errFl	|=ES1DMapR2gt(&b1,&gt,&b2,&z,r0);
	  SetSplY(z-EZz0);
	  a11	=Y[1];
	  a12	=Y[3];
	  a21	=dY[1];
	  a22	=dY[3];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-=Y[0]*f00+Y[2]*fy00;
	  B[2]	-=dY[0]*f00+dY[2]*fy00;
	  B[1]	-=Y[0]*fx00+Y[2]*f200;
	  B[3]	-=dY[0]*fx00+dY[2]*f200;
	   f10	=a22*B[0]-a12*B[2];
	  fy10	=a11*B[2]-a21*B[0];
	  fx10	=a22*B[1]-a12*B[3];
	  f210	=a11*B[3]-a21*B[1];
	  f0Map[k10]	= f10;
	  fxMap[k10]	=fx10;
	  fyMap[k10]	=fy10;
	  f2Map[k10]	=f210;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-=Y[0]*F00+Y[2]*Fy00;
	  F[2]	-=dY[0]*F00+dY[2]*Fy00;
	  F[1]	-=Y[0]*Fx00+Y[2]*F200;
	  F[3]	-=dY[0]*Fx00+dY[2]*F200;
	   F10	=a22*F[0]-a12*F[2];
	  Fy10	=a11*F[2]-a21*F[0];
	  Fx10	=a22*F[1]-a12*F[3];
	  F210	=a11*F[3]-a21*F[1];
	  F0Map[k10]	= F10;
	  FxMap[k10]	=Fx10;
	  FyMap[k10]	=Fy10;
	  F2Map[k10]	=F210;
	  r	=0.5*(r+r0);
	  z	=0.5*(z+EZz0);
	  SetSplX(r-r0);
	  SetSplY(z-EZz0);
	  d	=1e+20;
	  for(ii=0; ii < ESNa1; ii++){
	    kk	=ESNp1*ii;
	    for(jj=0; jj < ESNp; jj++){
	      dx	=r-ESsr[kk];
	      dy	=z-ESsz[kk];
	      dx	=dx*dx+dy*dy;
	      if(d > dx){
		d	=dx;
		a	=ESsa[ii];
		gt	=ESgt[jj];
	      }
	      kk++;
	    }
	  }
	  if(ES2DMapLab2Flx(&a,&gt,r,z)){
	    errFl	|=1;
	    putchar('\a');
	  }
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  A[0]	=Y[1]*X[1];
	  A[1]	=Y[1]*X[3];
	  A[2]	=Y[3]*X[1];
	  A[3]	=Y[3]*X[3];
	  B[0]-=(Y[0]*(X[0]*f00+X[1]*f01+X[2]*fx00+X[3]*fx01)
		 +Y[1]*(X[0]*f10+X[2]*fx10)
		 +Y[2]*(X[0]*fy00+X[1]*fy01+X[2]*f200+X[3]*f201)
		 +Y[3]*(X[0]*fy10+X[2]*f210)
		 );
	  A[4]	=Y[1]*dX[1];
	  A[5]	=Y[1]*dX[3];
	  A[6]	=Y[3]*dX[1];
	  A[7]	=Y[3]*dX[3];
	  B[1]-=(Y[0]*(dX[0]*f00+dX[1]*f01+dX[2]*fx00+dX[3]*fx01)
		 +Y[1]*(dX[0]*f10+dX[2]*fx10)
		 +Y[2]*(dX[0]*fy00+dX[1]*fy01+dX[2]*f200+dX[3]*f201)
		 +Y[3]*(dX[0]*fy10+dX[2]*f210)
		 );
	  A[8]	=dY[1]*X[1];
	  A[9]	=dY[1]*X[3];
	  A[10]	=dY[3]*X[1];
	  A[11]	=dY[3]*X[3];
	  B[2]-=(dY[0]*(X[0]*f00+X[1]*f01+X[2]*fx00+X[3]*fx01)
		 +dY[1]*(X[0]*f10+X[2]*fx10)
		 +dY[2]*(X[0]*fy00+X[1]*fy01+X[2]*f200+X[3]*f201)
		 +dY[3]*(X[0]*fy10+X[2]*f210)
		 );
	  A[12]	=dY[1]*dX[1];
	  A[13]	=dY[1]*dX[3];
	  A[14]	=dY[3]*dX[1];
	  A[15]	=dY[3]*dX[3];
	  B[3]-=(dY[0]*(dX[0]*f00+dX[1]*f01+dX[2]*fx00+dX[3]*fx01)
		 +dY[1]*(dX[0]*f10+dX[2]*fx10)
		 +dY[2]*(dX[0]*fy00+dX[1]*fy01+dX[2]*f200+dX[3]*f201)
		 +dY[3]*(dX[0]*fy10+dX[2]*f210)
		 );

	  LUdcmp(A,4,inA,&d);
	  LUbksb(A,4,inA,B);
	  f0Map[k11]	=B[0];
	  fxMap[k11]	=B[1];
	  fyMap[k11]	=B[2];
	  f2Map[k11]	=B[3];
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]-=(Y[0]*(X[0]*F00+X[1]*F01+X[2]*Fx00+X[3]*Fx01)
		 +Y[1]*(X[0]*F10+X[2]*Fx10)
		 +Y[2]*(X[0]*Fy00+X[1]*Fy01+X[2]*F200+X[3]*F201)
		 +Y[3]*(X[0]*Fy10+X[2]*F210)
		 );
	  F[1]-=(Y[0]*(dX[0]*F00+dX[1]*F01+dX[2]*Fx00+dX[3]*Fx01)
		 +Y[1]*(dX[0]*F10+dX[2]*Fx10)
		 +Y[2]*(dX[0]*Fy00+dX[1]*Fy01+dX[2]*F200+dX[3]*F201)
		 +Y[3]*(dX[0]*Fy10+dX[2]*F210)
		 );
	  F[2]-=(dY[0]*(X[0]*F00+X[1]*F01+X[2]*Fx00+X[3]*Fx01)
		 +dY[1]*(X[0]*F10+X[2]*Fx10)
		 +dY[2]*(X[0]*Fy00+X[1]*Fy01+X[2]*F200+X[3]*F201)
		 +dY[3]*(X[0]*Fy10+X[2]*F210)
		 );
	  F[3]-=(dY[0]*(dX[0]*F00+dX[1]*F01+dX[2]*Fx00+dX[3]*Fx01)
		 +dY[1]*(dX[0]*F10+dX[2]*Fx10)
		 +dY[2]*(dX[0]*Fy00+dX[1]*Fy01+dX[2]*F200+dX[3]*F201)
		 +dY[3]*(dX[0]*Fy10+dX[2]*F210)
		 );
	  LUbksb(A,4,inA,F);
	  F0Map[k11]	=F[0];
	  FxMap[k11]	=F[1];
	  FyMap[k11]	=F[2];
	  F2Map[k11]	=F[3];
	  break;
	case 2:
	  errFl	|=ES1DMapZ2gt(&gt,&r,&b2,EZz0);
	  gt	=EZcgp-gt;
	  SetSplX(r-r0);
	  a11	=X[0];
	  a12	=X[2];
	  a21	=dX[0];
	  a22	=dX[2];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-= X[1]*f01 + X[3]*fx01;
	  B[1]	-=dX[1]*f01 +dX[3]*fx01;
	  B[2]	-= X[1]*fy01+ X[3]*f201;
	  B[3]	-=dX[1]*fy01+dX[3]*f201;
	  f00	=a22*B[0]-a12*B[1];
	  fx00	=a11*B[1]-a21*B[0];
	  fy00	=a22*B[2]-a12*B[3];
	  f200	=a11*B[3]-a21*B[2];
	  f0Map[k00]	= f00;
	  fxMap[k00]	=fx00;
	  fyMap[k00]	=fy00;
	  f2Map[k00]	=f200;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-= X[1]*F01 + X[3]*Fx01;
	  F[1]	-=dX[1]*F01 +dX[3]*Fx01;
	  F[2]	-= X[1]*Fy01+ X[3]*F201;
	  F[3]	-=dX[1]*Fy01+dX[3]*F201;
	  F00	=a22*F[0]-a12*F[1];
	  Fx00	=a11*F[1]-a21*F[0];
	  Fy00	=a22*F[2]-a12*F[3];
	  F200	=a11*F[3]-a21*F[2];
	  F0Map[k00]	= F00;
	  FxMap[k00]	=Fx00;
	  FyMap[k00]	=Fy00;
	  F2Map[k00]	=F200;

	  errFl	|=ES1DMapR2gt(&b1,&gt,&b2,&z,r1);
	  SetSplY(z-EZz0);
	  a11	=Y[1];
	  a12	=Y[3];
	  a21	=dY[1];
	  a22	=dY[3];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-=Y[0]*f01+Y[2]*fy01;
	  B[2]	-=dY[0]*f01+dY[2]*fy01;
	  B[1]	-=Y[0]*fx01+Y[2]*f201;
	  B[3]	-=dY[0]*fx01+dY[2]*f201;
	   f11	=a22*B[0]-a12*B[2];
	  fy11	=a11*B[2]-a21*B[0];
	  fx11	=a22*B[1]-a12*B[3];
	  f211	=a11*B[3]-a21*B[1];
	  f0Map[k11]	= f11;
	  fxMap[k11]	=fx11;
	  fyMap[k11]	=fy11;
	  f2Map[k11]	=f211;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-=Y[0]*F01+Y[2]*Fy01;
	  F[2]	-=dY[0]*F01+dY[2]*Fy01;
	  F[1]	-=Y[0]*Fx01+Y[2]*F201;
	  F[3]	-=dY[0]*Fx01+dY[2]*F201;
	   F11	=a22*F[0]-a12*F[2];
	  Fy11	=a11*F[2]-a21*F[0];
	  Fx11	=a22*F[1]-a12*F[3];
	  F211	=a11*F[3]-a21*F[1];
	  F0Map[k11]	= F11;
	  FxMap[k11]	=Fx11;
	  FyMap[k11]	=Fy11;
	  F2Map[k11]	=F211;

	  r	=0.5*(r+r1);
	  z	=0.5*(z+EZz0);
	  SetSplX(r-r0);
	  SetSplY(z-EZz0);
	  d	=1e+20;
	  for(ii=0; ii < ESNa1; ii++){
	    kk	=ESNp1*ii;
	    for(jj=0; jj < ESNp; jj++){
	      dx	=r-ESsr[kk];
	      dy	=z-ESsz[kk];
	      dx	=dx*dx+dy*dy;
	      if(d > dx){
		d	=dx;
		a	=ESsa[ii];
		gt	=ESgt[jj];
	      }
	      kk++;
	    }
	  }
	  if(ES2DMapLab2Flx(&a,&gt,r,z)){
	    errFl	|=1;
	    putchar('\a');
	  }
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  A[0]	=Y[1]*X[0];
	  A[1]	=Y[1]*X[2];
	  A[2]	=Y[3]*X[0];
	  A[3]	=Y[3]*X[2];
	  B[0]-=(Y[0]*(X[0]*f00+X[1]*f01+X[2]*fx00+X[3]*fx01)
		 +Y[1]*(X[1]*f11+X[3]*fx11)
		 +Y[2]*(X[0]*fy00+X[1]*fy01+X[2]*f200+X[3]*f201)
		 +Y[3]*(X[1]*fy11+X[3]*f211)
		 );
	  A[4]	=Y[1]*dX[0];
	  A[5]	=Y[1]*dX[2];
	  A[6]	=Y[3]*dX[0];
	  A[7]	=Y[3]*dX[2];
	  B[1]-=
	    (Y[0]*(dX[0]*f00+dX[1]*f01+dX[2]*fx00+dX[3]*fx01)
	     +Y[1]*(dX[1]*f11+dX[3]*fx11)
	     +Y[2]*(dX[0]*fy00+dX[1]*fy01+dX[2]*f200+dX[3]*f201)
	     +Y[3]*(dX[1]*fy11+dX[3]*f211)
	     );
	  A[8]	=dY[1]*X[0];
	  A[9]	=dY[1]*X[2];
	  A[10]	=dY[3]*X[0];
	  A[11]	=dY[3]*X[2];
	  B[2]-=(dY[0]*(X[0]*f00+X[1]*f01+X[2]*fx00+X[3]*fx01)
		 +dY[1]*(X[1]*f11+X[3]*fx11)
		 +dY[2]*(X[0]*fy00+X[1]*fy01+X[2]*f200+X[3]*f201)
		 +dY[3]*(X[1]*fy11+X[3]*f211)
		 );
	  A[12]	=dY[1]*dX[0];
	  A[13]	=dY[1]*dX[2];
	  A[14]	=dY[3]*dX[0];
	  A[15]	=dY[3]*dX[2];
	  B[3]-=(dY[0]*(dX[0]*f00+dX[1]*f01+dX[2]*fx00+dX[3]*fx01)
		 +dY[1]*(dX[1]*f11+dX[3]*fx11)
		 +dY[2]*(dX[0]*fy00+dX[1]*fy01+dX[2]*f200+dX[3]*f201)
		 +dY[3]*(dX[1]*fy11+dX[3]*f211)
		 );

	  LUdcmp(A,4,inA,&d);
	  LUbksb(A,4,inA,B);
	  f0Map[k10]	=B[0];
	  fxMap[k10]	=B[1];
	  fyMap[k10]	=B[2];
	  f2Map[k10]	=B[3];

	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]-=(Y[0]*(X[0]*F00+X[1]*F01+X[2]*Fx00+X[3]*Fx01)
		 +Y[1]*(X[1]*F11+X[3]*Fx11)
		 +Y[2]*(X[0]*Fy00+X[1]*Fy01+X[2]*F200+X[3]*F201)
		 +Y[3]*(X[1]*Fy11+X[3]*F211)
		 );
	  F[1]-=
	    (Y[0]*(dX[0]*F00+dX[1]*F01+dX[2]*Fx00+dX[3]*Fx01)
	     +Y[1]*(dX[1]*F11+dX[3]*Fx11)
	     +Y[2]*(dX[0]*Fy00+dX[1]*Fy01+dX[2]*F200+dX[3]*F201)
	     +Y[3]*(dX[1]*Fy11+dX[3]*F211)
	     );
	  F[2]-=(dY[0]*(X[0]*F00+X[1]*F01+X[2]*Fx00+X[3]*Fx01)
		 +dY[1]*(X[1]*F11+X[3]*Fx11)
		 +dY[2]*(X[0]*Fy00+X[1]*Fy01+X[2]*F200+X[3]*F201)
		 +dY[3]*(X[1]*Fy11+X[3]*F211)
		 );
	  F[3]-=(dY[0]*(dX[0]*F00+dX[1]*F01+dX[2]*Fx00+dX[3]*Fx01)
		 +dY[1]*(dX[1]*F11+dX[3]*Fx11)
		 +dY[2]*(dX[0]*Fy00+dX[1]*Fy01+dX[2]*F200+dX[3]*F201)
		 +dY[3]*(dX[1]*Fy11+dX[3]*F211)
		 );
	  LUbksb(A,4,inA,F);
	  F0Map[k10]	=F[0];
	  FxMap[k10]	=F[1];
	  FyMap[k10]	=F[2];
	  F2Map[k10]	=F[3];
	  break;
	case 3:
	  errFl	|=ES1DMapR2gt(&b1,&gt,&b2,&z,r0);
	  SetSplY(z-EZz0);
	  a11	=Y[1];
	  a12	=Y[3];
	  a21	=dY[1];
	  a22	=dY[3];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-=Y[0]*f00+Y[2]*fy00;
	  B[2]	-=dY[0]*f00+dY[2]*fy00;
	  B[1]	-=Y[0]*fx00+Y[2]*f200;
	  B[3]	-=dY[0]*fx00+dY[2]*f200;
	   f10	=a22*B[0]-a12*B[2];
	  fy10	=a11*B[2]-a21*B[0];
	  fx10	=a22*B[1]-a12*B[3];
	  f210	=a11*B[3]-a21*B[1];
	  f0Map[k10]	= f10;
	  fxMap[k10]	=fx10;
	  fyMap[k10]	=fy10;
	  f2Map[k10]	=f210;

	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-=Y[0]*F00+Y[2]*Fy00;
	  F[2]	-=dY[0]*F00+dY[2]*Fy00;
	  F[1]	-=Y[0]*Fx00+Y[2]*F200;
	  F[3]	-=dY[0]*Fx00+dY[2]*F200;
	   F10	=a22*F[0]-a12*F[2];
	  Fy10	=a11*F[2]-a21*F[0];
	  Fx10	=a22*F[1]-a12*F[3];
	  F210	=a11*F[3]-a21*F[1];
	  F0Map[k10]	= F10;
	  FxMap[k10]	=Fx10;
	  FyMap[k10]	=Fy10;
	  F2Map[k10]	=F210;

	  errFl	|=ES1DMapR2gt(&b1,&gt,&b2,&z,r1);
	  SetSplY(z-EZz0);
	  a11	=Y[1];
	  a12	=Y[3];
	  a21	=dY[1];
	  a22	=dY[3];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-=Y[0]*f01+Y[2]*fy01;
	  B[2]	-=dY[0]*f01+dY[2]*fy01;
	  B[1]	-=Y[0]*fx01+Y[2]*f201;
	  B[3]	-=dY[0]*fx01+dY[2]*f201;
	   f11	=a22*B[0]-a12*B[2];
	  fy11	=a11*B[2]-a21*B[0];
	  fx11	=a22*B[1]-a12*B[3];
	  f211	=a11*B[3]-a21*B[1];
	  f0Map[k11]	= f11;
	  fxMap[k11]	=fx11;
	  fyMap[k11]	=fy11;
	  f2Map[k11]	=f211;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-=Y[0]*F01+Y[2]*Fy01;
	  F[2]	-=dY[0]*F01+dY[2]*Fy01;
	  F[1]	-=Y[0]*Fx01+Y[2]*F201;
	  F[3]	-=dY[0]*Fx01+dY[2]*F201;
	   F11	=a22*F[0]-a12*F[2];
	  Fy11	=a11*F[2]-a21*F[0];
	  Fx11	=a22*F[1]-a12*F[3];
	  F211	=a11*F[3]-a21*F[1];
	  F0Map[k11]	= F11;
	  FxMap[k11]	=Fx11;
	  FyMap[k11]	=Fy11;
	  F2Map[k11]	=F211;
	  break;
	case 4:
	  errFl	|=ES1DMapZ2gt(&gt,&r,&b2,z1);
	  gt	=EZcgp-gt;
	  SetSplX(r-r0);
	  a11	=X[0];
	  a12	=X[2];
	  a21	=dX[0];
	  a22	=dX[2];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-= X[1]*f11 + X[3]*fx11;
	  B[1]	-=dX[1]*f11 +dX[3]*fx11;
	  B[2]	-= X[1]*fy11+ X[3]*f211;
	  B[3]	-=dX[1]*fy11+dX[3]*f211;
	  f10	=a22*B[0]-a12*B[1];
	  fx10	=a11*B[1]-a21*B[0];
	  fy10	=a22*B[2]-a12*B[3];
	  f210	=a11*B[3]-a21*B[2];
	  f0Map[k10]	= f10;
	  fxMap[k10]	=fx10;
	  fyMap[k10]	=fy10;
	  f2Map[k10]	=f210;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-= X[1]*F11 + X[3]*Fx11;
	  F[1]	-=dX[1]*F11 +dX[3]*Fx11;
	  F[2]	-= X[1]*Fy11+ X[3]*F211;
	  F[3]	-=dX[1]*Fy11+dX[3]*F211;
	  F10	=a22*F[0]-a12*F[1];
	  Fx10	=a11*F[1]-a21*F[0];
	  Fy10	=a22*F[2]-a12*F[3];
	  F210	=a11*F[3]-a21*F[2];
	  F0Map[k10]	= F10;
	  FxMap[k10]	=Fx10;
	  FyMap[k10]	=Fy10;
	  F2Map[k10]	=F210;

	  errFl	|=ES1DMapR2gt(&gt,&b1,&z,&b2,r1);
	  SetSplY(z-EZz0);
	  a11	=Y[0];
	  a12	=Y[2];
	  a21	=dY[0];
	  a22	=dY[2];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-=Y[1]*f11+Y[3]*fy11;
	  B[2]	-=dY[1]*f11+dY[3]*fy11;
	  B[1]	-=Y[1]*fx11+Y[3]*f211;
	  B[3]	-=dY[1]*fx11+dY[3]*f211;
	   f01	=a22*B[0]-a12*B[2];
	  fy01	=a11*B[2]-a21*B[0];
	  fx01	=a22*B[1]-a12*B[3];
	  f201	=a11*B[3]-a21*B[1];
	  f0Map[k01]	= f01;
	  fxMap[k01]	=fx01;
	  fyMap[k01]	=fy01;
	  f2Map[k01]	=f201;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-=Y[1]*F11+Y[3]*Fy11;
	  F[2]	-=dY[1]*F11+dY[3]*Fy11;
	  F[1]	-=Y[1]*Fx11+Y[3]*F211;
	  F[3]	-=dY[1]*Fx11+dY[3]*F211;
	   F01	=a22*F[0]-a12*F[2];
	  Fy01	=a11*F[2]-a21*F[0];
	  Fx01	=a22*F[1]-a12*F[3];
	  F201	=a11*F[3]-a21*F[1];
	  F0Map[k01]	= F01;
	  FxMap[k01]	=Fx01;
	  FyMap[k01]	=Fy01;
	  F2Map[k01]	=F201;

	  r	=0.5*(r+r1);
	  z	=0.5*(z+z1);
	  SetSplX(r-r0);
	  SetSplY(z-EZz0);
	  d	=1e+20;
	  for(ii=0; ii < ESNa1; ii++){
	    kk	=ESNp1*ii;
	    for(jj=0; jj < ESNp; jj++){
	      dx	=r-ESsr[kk];
	      dy	=z-ESsz[kk];
	      dx	=dx*dx+dy*dy;
	      if(d > dx){
		d	=dx;
		a	=ESsa[ii];
		gt	=ESgt[jj];
	      }
	      kk++;
	    }
	  }
	  if(ES2DMapLab2Flx(&a,&gt,r,z)){
	    errFl	|=1;
	    putchar('\a');
	  }
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  A[0]	=Y[0]*X[0];
	  A[1]	=Y[0]*X[2];
	  A[2]	=Y[2]*X[0];
	  A[3]	=Y[2]*X[2];
	  B[0]-=(Y[0]*(X[1]*f01+X[3]*fx01)
		 +Y[1]*(X[0]*f10+X[1]*f11+X[2]*fx10+X[3]*fx11)
		 +Y[2]*(X[1]*fy01+X[3]*f201)
		 +Y[3]*(X[0]*fy10+X[1]*fy11+X[2]*f210+X[3]*f211)
		 );
	  A[4]	=Y[0]*dX[0];
	  A[5]	=Y[0]*dX[2];
	  A[6]	=Y[2]*dX[0];
	  A[7]	=Y[2]*dX[2];
	  B[1]-=
	    (Y[0]*(dX[1]*f01+dX[3]*fx01)
	     +Y[1]*(dX[0]*f10+dX[1]*f11+dX[2]*fx10+dX[3]*fx11)
	     +Y[2]*(dX[1]*fy01+dX[3]*f201)
	     +Y[3]*(dX[0]*fy10+dX[1]*fy11+dX[2]*f210+dX[3]*f211)
	     );
	  A[8]	=dY[0]*X[0];
	  A[9]	=dY[0]*X[2];
	  A[10]	=dY[2]*X[0];
	  A[11]	=dY[2]*X[2];
	  B[2]-=(dY[0]*(X[1]*f01+X[3]*fx01)
		 +dY[1]*(X[0]*f10+X[1]*f11+X[2]*fx10+X[3]*fx11)
		 +dY[2]*(X[1]*fy01+X[3]*f201)
		 +dY[3]*(X[0]*fy10+X[1]*fy11+X[2]*f210+X[3]*f211)
		 );
	  A[12]	=dY[0]*dX[0];
	  A[13]	=dY[0]*dX[2];
	  A[14]	=dY[2]*dX[0];
	  A[15]	=dY[2]*dX[2];
	  B[3]-=(dY[0]*(dX[1]*f01+dX[3]*fx01)
		 +dY[1]*(dX[0]*f10+dX[1]*f11+dX[2]*fx10+dX[3]*fx11)
		 +dY[2]*(dX[1]*fy01+dX[3]*f201)
		 +dY[3]*(dX[0]*fy10+dX[1]*fy11+dX[2]*f210+dX[3]*f211)
		 );

	  LUdcmp(A,4,inA,&d);
	  LUbksb(A,4,inA,B);
	  f0Map[k00]	=B[0];
	  fxMap[k00]	=B[1];
	  fyMap[k00]	=B[2];
	  f2Map[k00]	=B[3];

	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]-=(Y[0]*(X[1]*F01+X[3]*Fx01)
		 +Y[1]*(X[0]*F10+X[1]*F11+X[2]*Fx10+X[3]*Fx11)
		 +Y[2]*(X[1]*Fy01+X[3]*F201)
		 +Y[3]*(X[0]*Fy10+X[1]*Fy11+X[2]*F210+X[3]*F211)
		 );
	  F[1]-=
	    (Y[0]*(dX[1]*F01+dX[3]*Fx01)
	     +Y[1]*(dX[0]*F10+dX[1]*F11+dX[2]*Fx10+dX[3]*Fx11)
	     +Y[2]*(dX[1]*Fy01+dX[3]*F201)
	     +Y[3]*(dX[0]*Fy10+dX[1]*Fy11+dX[2]*F210+dX[3]*F211)
	     );
	  F[2]-=(dY[0]*(X[1]*F01+X[3]*Fx01)
		 +dY[1]*(X[0]*F10+X[1]*F11+X[2]*Fx10+X[3]*Fx11)
		 +dY[2]*(X[1]*Fy01+X[3]*F201)
		 +dY[3]*(X[0]*Fy10+X[1]*Fy11+X[2]*F210+X[3]*F211)
		 );
	  F[3]-=(dY[0]*(dX[1]*F01+dX[3]*Fx01)
		 +dY[1]*(dX[0]*F10+dX[1]*F11+dX[2]*Fx10+dX[3]*Fx11)
		 +dY[2]*(dX[1]*Fy01+dX[3]*F201)
		 +dY[3]*(dX[0]*Fy10+dX[1]*Fy11+dX[2]*F210+dX[3]*F211)
		 );
	  LUbksb(A,4,inA,F);
	  F0Map[k00]	=F[0];
	  FxMap[k00]	=F[1];
	  FyMap[k00]	=F[2];
	  F2Map[k00]	=F[3];
	  break;
	case 5:
	  errFl	|=2;
	  break;
	case 6:
	  errFl	|=ES1DMapZ2gt(&gt,&r,&b2,EZz0);
	  gt	=EZcgp-gt;
	  SetSplX(r-r0);
	  a11	=X[0];
	  a12	=X[2];
	  a21	=dX[0];
	  a22	=dX[2];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-= X[1]*f01 + X[3]*fx01;
	  B[1]	-=dX[1]*f01 +dX[3]*fx01;
	  B[2]	-= X[1]*fy01+ X[3]*f201;
	  B[3]	-=dX[1]*fy01+dX[3]*f201;
	  f00	=a22*B[0]-a12*B[1];
	  fx00	=a11*B[1]-a21*B[0];
	  fy00	=a22*B[2]-a12*B[3];
	  f200	=a11*B[3]-a21*B[2];
	  f0Map[k00]	= f00;
	  fxMap[k00]	=fx00;
	  fyMap[k00]	=fy00;
	  f2Map[k00]	=f200;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-= X[1]*F01 + X[3]*Fx01;
	  F[1]	-=dX[1]*F01 +dX[3]*Fx01;
	  F[2]	-= X[1]*Fy01+ X[3]*F201;
	  F[3]	-=dX[1]*Fy01+dX[3]*F201;
	  F00	=a22*F[0]-a12*F[1];
	  Fx00	=a11*F[1]-a21*F[0];
	  Fy00	=a22*F[2]-a12*F[3];
	  F200	=a11*F[3]-a21*F[2];
	  F0Map[k00]	= F00;
	  FxMap[k00]	=Fx00;
	  FyMap[k00]	=Fy00;
	  F2Map[k00]	=F200;

	  errFl	|=ES1DMapZ2gt(&gt,&r,&b2,z1);
	  gt	=EZcgp-gt;
	  SetSplX(r-r0);
	  a11	=X[0];
	  a12	=X[2];
	  a21	=dX[0];
	  a22	=dX[2];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-= X[1]*f11 + X[3]*fx11;
	  B[1]	-=dX[1]*f11 +dX[3]*fx11;
	  B[2]	-= X[1]*fy11+ X[3]*f211;
	  B[3]	-=dX[1]*fy11+dX[3]*f211;
	  f10	=a22*B[0]-a12*B[1];
	  fx10	=a11*B[1]-a21*B[0];
	  fy10	=a22*B[2]-a12*B[3];
	  f210	=a11*B[3]-a21*B[2];
	  f0Map[k10]	= f10;
	  fxMap[k10]	=fx10;
	  fyMap[k10]	=fy10;
	  f2Map[k10]	=f210;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-= X[1]*F11 + X[3]*Fx11;
	  F[1]	-=dX[1]*F11 +dX[3]*Fx11;
	  F[2]	-= X[1]*Fy11+ X[3]*F211;
	  F[3]	-=dX[1]*Fy11+dX[3]*F211;
	  F10	=a22*F[0]-a12*F[1];
	  Fx10	=a11*F[1]-a21*F[0];
	  Fy10	=a22*F[2]-a12*F[3];
	  F210	=a11*F[3]-a21*F[2];
	  F0Map[k10]	= F10;
	  FxMap[k10]	=Fx10;
	  FyMap[k10]	=Fy10;
	  F2Map[k10]	=F210;
	  break;
	case 7:
	  errFl	|=ES1DMapR2gt(&b1,&gt,&b2,&z,r0);
	  SetSplY(z-EZz0);
	  a11	=Y[1];
	  a12	=Y[3];
	  a21	=dY[1];
	  a22	=dY[3];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-=Y[0]*f00+Y[2]*fy00;
	  B[2]	-=dY[0]*f00+dY[2]*fy00;
	  B[1]	-=Y[0]*fx00+Y[2]*f200;
	  B[3]	-=dY[0]*fx00+dY[2]*f200;
	   f10	=a22*B[0]-a12*B[2];
	  fy10	=a11*B[2]-a21*B[0];
	  fx10	=a22*B[1]-a12*B[3];
	  f210	=a11*B[3]-a21*B[1];
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-=Y[0]*F00+Y[2]*Fy00;
	  F[2]	-=dY[0]*F00+dY[2]*Fy00;
	  F[1]	-=Y[0]*Fx00+Y[2]*F200;
	  F[3]	-=dY[0]*Fx00+dY[2]*F200;
	   F10	=a22*F[0]-a12*F[2];
	  Fy10	=a11*F[2]-a21*F[0];
	  Fx10	=a22*F[1]-a12*F[3];
	  F210	=a11*F[3]-a21*F[1];

	  errFl	|=ES1DMapZ2gt(&gt,&r,&b2,z1);
	  gt	=EZcgp-gt;
	  SetSplX(r-r0);
	  a11	=X[0];
	  a12	=X[2];
	  a21	=dX[0];
	  a22	=dX[2];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-= X[1]*f11 + X[3]*fx11;
	  B[1]	-=dX[1]*f11 +dX[3]*fx11;
	  B[2]	-= X[1]*fy11+ X[3]*f211;
	  B[3]	-=dX[1]*fy11+dX[3]*f211;
	  f0Map[k10]	=0.5*(f10+a22*B[0]-a12*B[1]);
	  fxMap[k10]	=0.5*(fx10+a11*B[1]-a21*B[0]);
	  fyMap[k10]	=0.5*(fy10+a22*B[2]-a12*B[3]);
	  f2Map[k10]	=0.5*(f210+a11*B[3]-a21*B[2]);
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-= X[1]*F11 + X[3]*Fx11;
	  F[1]	-=dX[1]*F11 +dX[3]*Fx11;
	  F[2]	-= X[1]*Fy11+ X[3]*F211;
	  F[3]	-=dX[1]*Fy11+dX[3]*F211;
	  F0Map[k10]	=0.5*(F10+a22*F[0]-a12*F[1]);
	  FxMap[k10]	=0.5*(Fx10+a11*F[1]-a21*F[0]);
	  FyMap[k10]	=0.5*(Fy10+a22*F[2]-a12*F[3]);
	  F2Map[k10]	=0.5*(F210+a11*F[3]-a21*F[2]);
	  break;
	case 8:
	  errFl	|=ES1DMapZ2gt(&gt,&b1,&r,z1);
	  SetSplX(r-r0);
	  a11	=X[1];
	  a12	=X[3];
	  a21	=dX[1];
	  a22	=dX[3];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-= X[0]*f10 + X[2]*fx10;
	  B[1]	-=dX[0]*f10 +dX[2]*fx10;
	  B[2]	-= X[0]*fy10+ X[2]*f210;
	  B[3]	-=dX[0]*fy10+dX[2]*f210;
	  f11	=a22*B[0]-a12*B[1];
	  fx11	=a11*B[1]-a21*B[0];
	  fy11	=a22*B[2]-a12*B[3];
	  f211	=a11*B[3]-a21*B[2];
	  f0Map[k11]	= f11;
	  fxMap[k11]	=fx11;
	  fyMap[k11]	=fy11;
	  f2Map[k11]	=f211;

	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-= X[0]*F10 + X[2]*Fx10;
	  F[1]	-=dX[0]*F10 +dX[2]*Fx10;
	  F[2]	-= X[0]*Fy10+ X[2]*F210;
	  F[3]	-=dX[0]*Fy10+dX[2]*F210;
	  F11	=a22*F[0]-a12*F[1];
	  Fx11	=a11*F[1]-a21*F[0];
	  Fy11	=a22*F[2]-a12*F[3];
	  F211	=a11*F[3]-a21*F[2];
	  F0Map[k11]	= F11;
	  FxMap[k11]	=Fx11;
	  FyMap[k11]	=Fy11;
	  F2Map[k11]	=F211;

	  errFl	|=ES1DMapR2gt(&gt,&b1,&z,&b2,r0);
	  SetSplY(z-EZz0);
	  a11	=Y[0];
	  a12	=Y[2];
	  a21	=dY[0];
	  a22	=dY[2];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-=Y[1]*f10+Y[3]*fy10;
	  B[2]	-=dY[1]*f10+dY[3]*fy10;
	  B[1]	-=Y[1]*fx10+Y[3]*f210;
	  B[3]	-=dY[1]*fx10+dY[3]*f210;
	   f00	=a22*B[0]-a12*B[2];
	  fy00	=a11*B[2]-a21*B[0];
	  fx00	=a22*B[1]-a12*B[3];
	  f200	=a11*B[3]-a21*B[1];
	  f0Map[k00]	= f00;
	  fxMap[k00]	=fx00;
	  fyMap[k00]	=fy00;
	  f2Map[k00]	=f200;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-=Y[1]*F10+Y[3]*Fy10;
	  F[2]	-=dY[1]*F10+dY[3]*Fy10;
	  F[1]	-=Y[1]*Fx10+Y[3]*F210;
	  F[3]	-=dY[1]*Fx10+dY[3]*F210;
	   F00	=a22*F[0]-a12*F[2];
	  Fy00	=a11*F[2]-a21*F[0];
	  Fx00	=a22*F[1]-a12*F[3];
	  F200	=a11*F[3]-a21*F[1];
	  F0Map[k00]	= F00;
	  FxMap[k00]	=Fx00;
	  FyMap[k00]	=Fy00;
	  F2Map[k00]	=F200;

	  r	=0.5*(r+r0);
	  z	=0.5*(z+z1);
	  SetSplX(r-r0);
	  SetSplY(z-EZz0);
	  d	=1e+20;
	  for(ii=0; ii < ESNa1; ii++){
	    kk	=ESNp1*ii;
	    for(jj=0; jj < ESNp; jj++){
	      dx	=r-ESsr[kk];
	      dy	=z-ESsz[kk];
	      dx	=dx*dx+dy*dy;
	      if(d > dx){
		d	=dx;
		a	=ESsa[ii];
		gt	=ESgt[jj];
	      }
	      kk++;
	    }
	  }
	  if(ES2DMapLab2Flx(&a,&gt,r,z)){
	    errFl	|=1;
	    putchar('\a');
	  }
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  A[0]	=Y[0]*X[1];
	  A[1]	=Y[0]*X[3];
	  A[2]	=Y[2]*X[1];
	  A[3]	=Y[2]*X[3];
	  B[0]-=(Y[0]*(X[0]*f00+X[2]*fx00)
		 +Y[1]*(X[0]*f10+X[1]*f11+X[2]*fx10+X[3]*fx11)
		 +Y[2]*(X[0]*fy00+X[2]*f200)
		 +Y[3]*(X[0]*fy10+X[1]*fy11+X[2]*f210+X[3]*f211)
		 );
	  A[4]	=Y[0]*dX[1];
	  A[5]	=Y[0]*dX[3];
	  A[6]	=Y[2]*dX[1];
	  A[7]	=Y[2]*dX[3];
	  B[1]-=
	    (Y[0]*(dX[0]*f00+dX[2]*fx00)
	     +Y[1]*(dX[0]*f10+dX[1]*f11+dX[2]*fx10+dX[3]*fx11)
	     +Y[2]*(dX[0]*fy00+dX[2]*f200)
	     +Y[3]*(dX[0]*fy10+dX[1]*fy11+dX[2]*f210+dX[3]*f211)
	     );
	  A[8]	=dY[0]*X[1];
	  A[9]	=dY[0]*X[3];
	  A[10]	=dY[2]*X[1];
	  A[11]	=dY[2]*X[3];
	  B[2]-=(dY[0]*(X[0]*f00+X[2]*fx00)
		 +dY[1]*(X[0]*f10+X[1]*f11+X[2]*fx10+X[3]*fx11)
		 +dY[2]*(X[0]*fy00+X[2]*f200)
		 +dY[3]*(X[0]*fy10+X[1]*fy11+X[2]*f210+X[3]*f211)
		 );
	  A[12]	=dY[0]*dX[1];
	  A[13]	=dY[0]*dX[3];
	  A[14]	=dY[2]*dX[1];
	  A[15]	=dY[2]*dX[3];
	  B[3]-=(dY[0]*(dX[0]*f00+dX[2]*fx00)
		 +dY[1]*(dX[0]*f10+dX[1]*f11+dX[2]*fx10+dX[3]*fx11)
		 +dY[2]*(dX[0]*fy00+dX[2]*f200)
		 +dY[3]*(dX[0]*fy10+dX[1]*fy11+dX[2]*f210+dX[3]*f211)
		 );

	  LUdcmp(A,4,inA,&d);
	  LUbksb(A,4,inA,B);
	  f0Map[k01]	=B[0];
	  fxMap[k01]	=B[1];
	  fyMap[k01]	=B[2];
	  f2Map[k01]	=B[3];
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]-=(Y[0]*(X[0]*F00+X[2]*Fx00)
		 +Y[1]*(X[0]*F10+X[1]*F11+X[2]*Fx10+X[3]*Fx11)
		 +Y[2]*(X[0]*Fy00+X[2]*F200)
		 +Y[3]*(X[0]*Fy10+X[1]*Fy11+X[2]*F210+X[3]*F211)
		 );
	  F[1]-=
	    (Y[0]*(dX[0]*F00+dX[2]*Fx00)
	     +Y[1]*(dX[0]*F10+dX[1]*F11+dX[2]*Fx10+dX[3]*Fx11)
	     +Y[2]*(dX[0]*Fy00+dX[2]*F200)
	     +Y[3]*(dX[0]*Fy10+dX[1]*Fy11+dX[2]*F210+dX[3]*F211)
	     );
	  F[2]-=(dY[0]*(X[0]*F00+X[2]*Fx00)
		 +dY[1]*(X[0]*F10+X[1]*F11+X[2]*Fx10+X[3]*Fx11)
		 +dY[2]*(X[0]*Fy00+X[2]*F200)
		 +dY[3]*(X[0]*Fy10+X[1]*Fy11+X[2]*F210+X[3]*F211)
		 );
	  F[3]-=(dY[0]*(dX[0]*F00+dX[2]*Fx00)
		 +dY[1]*(dX[0]*F10+dX[1]*F11+dX[2]*Fx10+dX[3]*Fx11)
		 +dY[2]*(dX[0]*Fy00+dX[2]*F200)
		 +dY[3]*(dX[0]*Fy10+dX[1]*Fy11+dX[2]*F210+dX[3]*F211)
		 );

	  LUbksb(A,4,inA,F);
	  F0Map[k01]	=F[0];
	  FxMap[k01]	=F[1];
	  FyMap[k01]	=F[2];
	  F2Map[k01]	=F[3];
	  break;
	case 9:
	  errFl	|=ES1DMapZ2gt(&gt,&b2,&r,EZz0);
	  SetSplX(r-r0);
	  a11	=X[1];
	  a12	=X[3];
	  a21	=dX[1];
	  a22	=dX[3];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-= X[0]*f00 + X[2]*fx00;
	  B[1]	-=dX[0]*f00 +dX[2]*fx00;
	  B[2]	-= X[0]*fy00+ X[2]*f200;
	  B[3]	-=dX[0]*fy00+dX[2]*f200;
	  f01	=a22*B[0]-a12*B[1];
	  fx01	=a11*B[1]-a21*B[0];
	  fy01	=a22*B[2]-a12*B[3];
	  f201	=a11*B[3]-a21*B[2];
	  f0Map[k01]	= f01;
	  fxMap[k01]	=fx01;
	  fyMap[k01]	=fy01;
	  f2Map[k01]	=f201;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-= X[0]*F00 + X[2]*Fx00;
	  F[1]	-=dX[0]*F00 +dX[2]*Fx00;
	  F[2]	-= X[0]*Fy00+ X[2]*F200;
	  F[3]	-=dX[0]*Fy00+dX[2]*F200;
	  F01	=a22*F[0]-a12*F[1];
	  Fx01	=a11*F[1]-a21*F[0];
	  Fy01	=a22*F[2]-a12*F[3];
	  F201	=a11*F[3]-a21*F[2];
	  F0Map[k01]	= F01;
	  FxMap[k01]	=Fx01;
	  FyMap[k01]	=Fy01;
	  F2Map[k01]	=F201;

	  errFl	|=ES1DMapZ2gt(&gt,&b1,&r,z1);
	  SetSplX(r-r0);
	  a11	=X[1];
	  a12	=X[3];
	  a21	=dX[1];
	  a22	=dX[3];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-= X[0]*f10 + X[2]*fx10;
	  B[1]	-=dX[0]*f10 +dX[2]*fx10;
	  B[2]	-= X[0]*fy10+ X[2]*f210;
	  B[3]	-=dX[0]*fy10+dX[2]*f210;
	  f11	=a22*B[0]-a12*B[1];
	  fx11	=a11*B[1]-a21*B[0];
	  fy11	=a22*B[2]-a12*B[3];
	  f211	=a11*B[3]-a21*B[2];
	  f0Map[k11]	= f11;
	  fxMap[k11]	=fx11;
	  fyMap[k11]	=fy11;
	  f2Map[k11]	=f211;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-= X[0]*F10 + X[2]*Fx10;
	  F[1]	-=dX[0]*F10 +dX[2]*Fx10;
	  F[2]	-= X[0]*Fy10+ X[2]*F210;
	  F[3]	-=dX[0]*Fy10+dX[2]*F210;
	  F11	=a22*F[0]-a12*F[1];
	  Fx11	=a11*F[1]-a21*F[0];
	  Fy11	=a22*F[2]-a12*F[3];
	  F211	=a11*F[3]-a21*F[2];
	  F0Map[k11]	= F11;
	  FxMap[k11]	=Fx11;
	  FyMap[k11]	=Fy11;
	  F2Map[k11]	=F211;
	  break;
	case 10:
	  errFl	|=2;
	  break;
	case 11:
	  errFl	|=ES1DMapZ2gt(&gt,&b1,&r,z1);
	  SetSplX(r-r0);
	  a11	=X[1];
	  a12	=X[3];
	  a21	=dX[1];
	  a22	=dX[3];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-= X[0]*f10 + X[2]*fx10;
	  B[1]	-=dX[0]*f10 +dX[2]*fx10;
	  B[2]	-= X[0]*fy10+ X[2]*f210;
	  B[3]	-=dX[0]*fy10+dX[2]*f210;
	  f11	=a22*B[0]-a12*B[1];
	  fx11	=a11*B[1]-a21*B[0];
	  fy11	=a22*B[2]-a12*B[3];
	  f211	=a11*B[3]-a21*B[2];
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-= X[0]*F10 + X[2]*Fx10;
	  F[1]	-=dX[0]*F10 +dX[2]*Fx10;
	  F[2]	-= X[0]*Fy10+ X[2]*F210;
	  F[3]	-=dX[0]*Fy10+dX[2]*F210;
	  F11	=a22*F[0]-a12*F[1];
	  Fx11	=a11*F[1]-a21*F[0];
	  Fy11	=a22*F[2]-a12*F[3];
	  F211	=a11*F[3]-a21*F[2];

	  errFl	|=ES1DMapR2gt(&b1,&gt,&b2,&z,r1);
	  SetSplY(z-EZz0);
	  a11	=Y[1];
	  a12	=Y[3];
	  a21	=dY[1];
	  a22	=dY[3];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-=Y[0]*f01+Y[2]*fy01;
	  B[2]	-=dY[0]*f01+dY[2]*fy01;
	  B[1]	-=Y[0]*fx01+Y[2]*f201;
	  B[3]	-=dY[0]*fx01+dY[2]*f201;
	  f0Map[k11]	=0.5*(f11+a22*B[0]-a12*B[2]);
	  fyMap[k11]	=0.5*(fy11+a11*B[2]-a21*B[0]);
	  fxMap[k11]	=0.5*(fx11+a22*B[1]-a12*B[3]);
	  f2Map[k11]	=0.5*(f211+a11*B[3]-a21*B[1]);
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-=Y[0]*F01+Y[2]*Fy01;
	  F[2]	-=dY[0]*F01+dY[2]*Fy01;
	  F[1]	-=Y[0]*Fx01+Y[2]*F201;
	  F[3]	-=dY[0]*Fx01+dY[2]*F201;
	  F0Map[k11]	=0.5*(F11+a22*F[0]-a12*F[2]);
	  FyMap[k11]	=0.5*(Fy11+a11*F[2]-a21*F[0]);
	  FxMap[k11]	=0.5*(Fx11+a22*F[1]-a12*F[3]);
	  F2Map[k11]	=0.5*(F211+a11*F[3]-a21*F[1]);
	  break;
	case 12:
	  errFl	|=ES1DMapR2gt(&gt,&b1,&z,&b2,r0);
	  SetSplY(z-EZz0);
	  a11	=Y[0];
	  a12	=Y[2];
	  a21	=dY[0];
	  a22	=dY[2];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-=Y[1]*f10+Y[3]*fy10;
	  B[2]	-=dY[1]*f10+dY[3]*fy10;
	  B[1]	-=Y[1]*fx10+Y[3]*f210;
	  B[3]	-=dY[1]*fx10+dY[3]*f210;
	   f00	=a22*B[0]-a12*B[2];
	  fy00	=a11*B[2]-a21*B[0];
	  fx00	=a22*B[1]-a12*B[3];
	  f200	=a11*B[3]-a21*B[1];
	  f0Map[k00]	= f00;
	  fxMap[k00]	=fx00;
	  fyMap[k00]	=fy00;
	  f2Map[k00]	=f200;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-=Y[1]*F10+Y[3]*Fy10;
	  F[2]	-=dY[1]*F10+dY[3]*Fy10;
	  F[1]	-=Y[1]*Fx10+Y[3]*F210;
	  F[3]	-=dY[1]*Fx10+dY[3]*F210;
	   F00	=a22*F[0]-a12*F[2];
	  Fy00	=a11*F[2]-a21*F[0];
	  Fx00	=a22*F[1]-a12*F[3];
	  F200	=a11*F[3]-a21*F[1];
	  F0Map[k00]	= F00;
	  FxMap[k00]	=Fx00;
	  FyMap[k00]	=Fy00;
	  F2Map[k00]	=F200;

	  errFl	|=ES1DMapR2gt(&gt,&b1,&z,&b2,r1);
	  SetSplY(z-EZz0);
	  a11	=Y[0];
	  a12	=Y[2];
	  a21	=dY[0];
	  a22	=dY[2];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-=Y[1]*f11+Y[3]*fy11;
	  B[2]	-=dY[1]*f11+dY[3]*fy11;
	  B[1]	-=Y[1]*fx11+Y[3]*f211;
	  B[3]	-=dY[1]*fx11+dY[3]*f211;
	   f01	=a22*B[0]-a12*B[2];
	  fy01	=a11*B[2]-a21*B[0];
	  fx01	=a22*B[1]-a12*B[3];
	  f201	=a11*B[3]-a21*B[1];
	  f0Map[k01]	= f01;
	  fxMap[k01]	=fx01;
	  fyMap[k01]	=fy01;
	  f2Map[k01]	=f201;
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-=Y[1]*F11+Y[3]*Fy11;
	  F[2]	-=dY[1]*F11+dY[3]*Fy11;
	  F[1]	-=Y[1]*Fx11+Y[3]*F211;
	  F[3]	-=dY[1]*Fx11+dY[3]*F211;
	   F01	=a22*F[0]-a12*F[2];
	  Fy01	=a11*F[2]-a21*F[0];
	  Fx01	=a22*F[1]-a12*F[3];
	  F201	=a11*F[3]-a21*F[1];
	  F0Map[k01]	= F01;
	  FxMap[k01]	=Fx01;
	  FyMap[k01]	=Fy01;
	  F2Map[k01]	=F201;
	  break;
	case 13:
	  errFl	|=ES1DMapZ2gt(&gt,&b2,&r,EZz0);
	  SetSplX(r-r0);
	  a11	=X[1];
	  a12	=X[3];
	  a21	=dX[1];
	  a22	=dX[3];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-= X[0]*f00 + X[2]*fx00;
	  B[1]	-=dX[0]*f00 +dX[2]*fx00;
	  B[2]	-= X[0]*fy00+ X[2]*f200;
	  B[3]	-=dX[0]*fy00+dX[2]*f200;
	  f01	=a22*B[0]-a12*B[1];
	  fx01	=a11*B[1]-a21*B[0];
	  fy01	=a22*B[2]-a12*B[3];
	  f201	=a11*B[3]-a21*B[2];
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-= X[0]*F00 + X[2]*Fx00;
	  F[1]	-=dX[0]*F00 +dX[2]*Fx00;
	  F[2]	-= X[0]*Fy00+ X[2]*F200;
	  F[3]	-=dX[0]*Fy00+dX[2]*F200;
	  F01	=a22*F[0]-a12*F[1];
	  Fx01	=a11*F[1]-a21*F[0];
	  Fy01	=a22*F[2]-a12*F[3];
	  F201	=a11*F[3]-a21*F[2];

	  errFl	|=ES1DMapR2gt(&gt,&b1,&z,&b2,r1);
	  SetSplY(z-EZz0);
	  a11	=Y[0];
	  a12	=Y[2];
	  a21	=dY[0];
	  a22	=dY[2];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-=Y[1]*f11+Y[3]*fy11;
	  B[2]	-=dY[1]*f11+dY[3]*fy11;
	  B[1]	-=Y[1]*fx11+Y[3]*f211;
	  B[3]	-=dY[1]*fx11+dY[3]*f211;
	  f0Map[k01]	=0.5*(f01+a22*B[0]-a12*B[2]);
	  fyMap[k01]	=0.5*(fy01+a11*B[2]-a21*B[0]);
	  fxMap[k01]	=0.5*(fx01+a22*B[1]-a12*B[3]);
	  f2Map[k01]	=0.5*(f201+a11*B[3]-a21*B[1]);
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-=Y[1]*F11+Y[3]*Fy11;
	  F[2]	-=dY[1]*F11+dY[3]*Fy11;
	  F[1]	-=Y[1]*Fx11+Y[3]*F211;
	  F[3]	-=dY[1]*Fx11+dY[3]*F211;
	  F0Map[k01]	=0.5*(F01+a22*F[0]-a12*F[2]);
	  FyMap[k01]	=0.5*(Fy01+a11*F[2]-a21*F[0]);
	  FxMap[k01]	=0.5*(Fx01+a22*F[1]-a12*F[3]);
	  F2Map[k01]	=0.5*(F201+a11*F[3]-a21*F[1]);
	  break;
	case 14:
	  errFl	|=ES1DMapZ2gt(&gt,&r,&b2,EZz0);
	  gt	=EZcgp-gt;
	  SetSplX(r-r0);
	  a11	=X[0];
	  a12	=X[2];
	  a21	=dX[0];
	  a22	=dX[2];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-= X[1]*f01 + X[3]*fx01;
	  B[1]	-=dX[1]*f01 +dX[3]*fx01;
	  B[2]	-= X[1]*fy01+ X[3]*f201;
	  B[3]	-=dX[1]*fy01+dX[3]*f201;
	  f00	=a22*B[0]-a12*B[1];
	  fx00	=a11*B[1]-a21*B[0];
	  fy00	=a22*B[2]-a12*B[3];
	  f200	=a11*B[3]-a21*B[2];
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-= X[1]*F01 + X[3]*Fx01;
	  F[1]	-=dX[1]*F01 +dX[3]*Fx01;
	  F[2]	-= X[1]*Fy01+ X[3]*F201;
	  F[3]	-=dX[1]*Fy01+dX[3]*F201;
	  F00	=a22*F[0]-a12*F[1];
	  Fx00	=a11*F[1]-a21*F[0];
	  Fy00	=a22*F[2]-a12*F[3];
	  F200	=a11*F[3]-a21*F[2];

	  errFl	|=ES1DMapR2gt(&gt,&b1,&z,&b2,r0);
	  SetSplY(z-EZz0);
	  a11	=Y[0];
	  a12	=Y[2];
	  a21	=dY[0];
	  a22	=dY[2];
	  b1	=1./(a11*a22-a12*a21);
	  a11	*=b1;
	  a12	*=b1;
	  a21	*=b1;
	  a22	*=b1;
	  ESGetFullFgY(F,B,a,gt);
	  B[1]	*=dr;
	  B[2]	*=dz;
	  B[3]	*=dr*dz;
	  B[0]	-=Y[1]*f10+Y[3]*fy10;
	  B[2]	-=dY[1]*f10+dY[3]*fy10;
	  B[1]	-=Y[1]*fx10+Y[3]*f210;
	  B[3]	-=dY[1]*fx10+dY[3]*f210;
	  f0Map[k00]	=0.5*(f00+a22*B[0]-a12*B[2]);
	  fyMap[k00]	=0.5*(fy00+a11*B[2]-a21*B[0]);
	  fxMap[k00]	=0.5*(fx00+a22*B[1]-a12*B[3]);
	  f2Map[k00]	=0.5*(f200+a11*B[3]-a21*B[1]);
	  F[1]	*=dr;
	  F[2]	*=dz;
	  F[3]	*=dr*dz;
	  F[0]	-=Y[1]*F10+Y[3]*Fy10;
	  F[2]	-=dY[1]*F10+dY[3]*Fy10;
	  F[1]	-=Y[1]*Fx10+Y[3]*F210;
	  F[3]	-=dY[1]*Fx10+dY[3]*F210;
	  F0Map[k00]	=0.5*(F00+a22*F[0]-a12*F[2]);
	  FyMap[k00]	=0.5*(Fy00+a11*F[2]-a21*F[0]);
	  FxMap[k00]	=0.5*(Fx00+a22*F[1]-a12*F[3]);
	  F2Map[k00]	=0.5*(F200+a11*F[3]-a21*F[1]);
	  break;
	default:
	  break;
	}
	k00Map[k]	=k00;
	k00+=4;
      }
      k++;
    }
  }
#ifdef XWIN
  if(errFl) CbStr2UserMessage("Some Cells are not found\n");
#endif
  return 0;
}

int ESWriteESIRZ(int iRec)
{
  int i;

  FILE *lfa;
  char FNm[16];

  extern char ESMessage[];
#ifdef XWIN
  extern char *CbUserMessage;
  CbUserMessage	=ESMessage;
#endif
  sprintf(FNm,"esiRZ.%2.2d",iRec);
  lfa	=fopen(FNm,"w");
  if(lfa == NULL){
    sprintf(ESMessage,"%s cannot be open",FNm);
    return(1);
  }

  time(&stTime);
  ltm	=localtime(&stTime);
  sprintf(ESMessage,"Date: %2.2d/%2.2d/%2.2d at %2.2d:%2.2d"
	  ,ltm->tm_mon+1,ltm->tm_mday,ltm->tm_year%100,ltm->tm_hour
	  ,ltm->tm_min); 
  fprintf(lfa,"!!! Do not edit this file. %s\n",ESMessage);
  fprintf(lfa,"%3d x %3d %5d - number of radial x vertical"
	  " intervals and data points\n",nrMap,nzMap,nMap1);
  fprintf(lfa,"%24s%24s%24s%24s\n","rBox[0]","zBox[0]","rBox[1]","zBox[1]");
  fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	  ,rBox[0],zBox[0],rBox[1],zBox[1]);

  fprintf(lfa,"%6s 0x%2s\n","k00","iA");
  for(i=0; i < nMap; i++){
    fprintf(lfa,"%6d %2.2x",k00Map[i],iAMap[i]);
    if((i+1)%8 == 0) putc('\n',lfa);
  }

  fprintf(lfa,"%24s%24s%24s%24s\n","gY","gY'_r","gY'_z","gY''_{rz}");
  for(i=0; i < nMap1; i++){
    fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	    ,f0Map[i],fxMap[i],fyMap[i],f2Map[i]);
  }
  fprintf(lfa,"%24s%24s%24s%24s\n","F","F'_r","F'_z","F''_{rz}");
  for(i=0; i < nMap1; i++){
    fprintf(lfa,"%24.16e%24.16e%24.16e%24.16e\n"
	    ,F0Map[i],FxMap[i],FyMap[i],F2Map[i]);
  }
  fclose(lfa);
  sprintf(ESMessage,"%s has been written",FNm);
    
#ifdef aXWIN
  if(esireadrz_("FNm")){
    sprintf(ESMessage,"%s Bad file",FNm);
  }
#endif
  return(0);
}
#endif

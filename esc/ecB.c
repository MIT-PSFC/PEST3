#include <math.h>
#include <stdio.h>
#include <stdlib.h>
extern double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
extern int ESNa,ESNa1,ESNp,ESNp1;
extern int ESMp,ESMp1,ES2Mp1;
extern int ESFp,ESFp1,ESnAF;
extern int NPlVac,*k2kPlV;
extern double*rPlVd,*zPlVd;
extern double *EScs1,*ESsn1;

extern double *ESaR0,*ESaZ0,*ESsa,*ESsb,*ESsb1a,*ESsb2a;
extern double *rcT,*rcT1a,*rcT2a,*rsT,*rsT1a,*rsT2a;

extern double *ECb,*ECb1a,*ECb2a;
extern double *ECr,*ECz;
extern double *ESgt;
extern double R0,Z0,ESRext,ESRBt;
extern int ESiAx;

#ifndef mzl_PL
extern int ESnLm,ESnLm1;
extern double*ESRLm,*ESZLm;
#endif

extern int NPlVac,*k2kPlV;
extern double *rPlVd,*zPlVd;

extern int ESnBL,ESnBL1,nf_X;
extern double ESaRx,ESaZx,ESDx[],ESLx[],ESVx[],ESgFPlVx;
extern double EZd_vb,EZgd_vb,b_X,z_X0,z_X2,z_X2n,t_X,singtX;
extern double *EZx0sep,*EZx1sep,*EZx2sep,*X2vb,*EZx0cs,*EZx0sn,*EZdx0cs,*EZdx0sn
,*X2vbc,*X2vbs;
extern double *EZrvb,*EZzvb,*EZrvbcs,*EZrvbsn,*drvbcs,*drvbsn;
extern double*Dx0,*Dx2,*Ddx2,*Dz2,*Dx0c,*Dx0s,*Dx2c,*Dx2s,*Ddx2c,*Ddx2s
,*Dz2c,*Dz2s;

extern double *EZrcs,*EZrsn,*EZd1rcs,*EZd1rsn,*EZd2rcs,*EZd2rsn,*EZz0,*EZd1z0,*EZd2z0;

extern double*Fcs,*Fsn,*dFcs,*dFsn,*EZxinf,*A2m,*EZyinf,*d1yinf,*d2yinf,*EZxgt,*EZygt;
extern double Eps_tr;

int EcInitPlVacData()
{
  int i,j,k,ki;
  double a,b;

  Eps_tr	=1e-4;
  R0		=ESaR0[0];
  Z0		=ESaZ0[0];
  z_X0		=Z0;
  rPlVd[0]	=R0;
  rPlVd[1]	=R0;
  a		=R0/3.;
  rPlVd[2]	=R0-a;
  rPlVd[3]	=R0+a;
  b		=a;
  zPlVd[0]	=b;
  zPlVd[1]	=-b;
  zPlVd[2]	=Z0;
  zPlVd[3]	=Z0;
  ESaRx		=R0;
  ESaZx		=-b;
  for(i=4 ; i < NPlVac; i++){
    rPlVd[i]	=rPlVd[4];
    zPlVd[i]	=zPlVd[2];
  }
  for(i=0; i < ESNa1; i++){
    EZz0[i]	=0.;
    EZd1z0[i]	=0.;
    EZd2z0[i]	=0.;
  }
  ki	=0;
  for(k=0; k < ESFp1; k++){
    for(i=0; i < ESNa1; i++){
      EZrcs[ki]	=0.;
      EZrsn[ki]	=0.;
      EZd1rcs[ki]	=0.;
      EZd1rsn[ki]	=0.;
      EZd2rcs[ki]	=0.;
      EZd2rsn[ki]	=0.;
      ki++;
    }
  }

  for(j= 0;j < ESNp1;j++){
    EZx0sep[j]= 0.;
    EZx1sep[j]= 0.;
    EZx2sep[j]= 0.;
  }
  for(k= 0; k < ESFp1; k++){
    EZx0cs[k]	=0.;
    EZx0sn[k]	=0.;
    EZrvbcs[k]	=0.;
    EZrvbsn[k]	=0.;
  }
  z_X2= 2.;
  for(i=0; i < 4; i++)
    k2kPlV[i]=1;
  while(i < NPlVac){
    k2kPlV[i]=0;
    rPlVd[i]=rPlVd[3];
    zPlVd[i]=zPlVd[3];
    i++;
  }
  return(0);
}

int ECRcs3D2Rcs2D(int n)
{
  int i,ki,k,km;
  R0	=ESaR0[n];
  Z0	=ESaZ0[n];
  ki	=ESNa1*n;
  for(i=0; i < ESNa1; i++){
    ECb[i]	=ESsb[ki];
    ECb1a[i]	=ESsb1a[ki];
    ECb2a[i]	=ESsb2a[ki];
    ki++;
  }
  ki	=ESnAF*n;
  km	=0;
  for(i=0; i < ESNa1; i++){
    EZrcs[km]	=rcT[ki];
    EZd1rcs[km]	=rcT1a[ki];
    EZd2rcs[km]	=rcT2a[ki];
    EZz0[km]	=rsT[ki];
    EZd1z0[km]	=rsT1a[ki];
    EZd2z0[km]	=rsT2a[ki];
    km++;
    ki++;
  }
  for(k=1; k < ESFp1; k++){
    for(i=0; i < ESNa1; i++){
      EZrcs[km]	=rcT[ki];
      EZd1rcs[km]	=rcT1a[ki];
      EZd2rcs[km]	=rcT2a[ki];
      EZrsn[km]	=rsT[ki];
      EZd1rsn[km]	=rsT1a[ki];
      EZd2rsn[km]	=rsT2a[ki];
      km++;
      ki++;
    }
  }
  return(0);
}

int ECRcs2D2Rcs3D(int n)
{
  int i,ki,k,km,in;

  ESaR0[n]	=R0;
  ESaZ0[n]	=Z0;
  in		=ESNa1*n;
  ki		=ESnAF*n;
  km		=0;
  for(i=0; i < ESNa1; i++){
    ESsb[in]	=ECb[i];
    ESsb1a[in]	=ECb1a[i];
    ESsb2a[in]	=ECb2a[i];
    rcT[ki]	=EZrcs[km];
    rcT1a[ki]	=EZd1rcs[km];
    rcT2a[ki]	=EZd2rcs[km];
    rsT[ki]	=EZz0[i];
    rsT1a[ki]	=EZd1z0[i];
    rsT2a[ki]	=EZd2z0[i];
    km++;
    ki++;
    in++;
  }
  for(k=1; k < ESFp1; k++){
    for(i=0; i < ESNa1; i++){
      rcT[ki]	=EZrcs[km];
      rcT1a[ki]	=EZd1rcs[km];
      rcT2a[ki]	=EZd2rcs[km];
      rsT[ki]	=EZrsn[km];
      rsT1a[ki]	=EZd1rsn[km];
      rsT2a[ki]	=EZd2rsn[km];
      km++;
      ki++;
    }
  }
  return(0);
}

extern double*rCc,*rCs,*rSc,*rSs,*drCc,*drCs,*drSc,*drSs;
extern double*EZgper,*EZgpei,*EZdgper,*EZdgpei;
extern int M0Mm,M0Lm,m;
static double *ac=NULL,*f,*G,*G1,*F;
static int *index;

int ECVirtualPlV2vdrPlV()
{
  static int i,j,k,ki,ji,kj,Isize=0;
  static int m,mm,mi,mj,ii,jj;
  static double cs,sn,d,s,ss,am,x_0,x_1,x_2;
  static int nX,neq;

  if(ESiAx == 0){
    for(i=0; i < ESNa1; i++){
      ki	=M0Mm*i+ESMp;
      mj	=M0Mm*i+ESMp;
      EZgpei[ki]	=0.;
      for(k=1; k <  ESMp1; k++){
	kj	= ki-k;
	jj	= ki+k;
	EZgper[jj]	=EZcr2*(EZgper[jj]+EZgper[kj]);
	EZgpei[jj]	=EZcr2*(EZgpei[jj]-EZgpei[kj]);
	EZdgper[jj]	=EZcr2*(EZdgper[jj]+EZdgper[kj]);
	EZdgpei[jj]	=EZcr2*(EZdgpei[jj]-EZdgpei[kj]);
	EZgper[kj]	=EZgper[jj];
	EZgpei[kj]	=-EZgpei[jj];
	EZdgper[kj]	=EZdgper[jj];
	EZdgpei[kj]	=-EZdgpei[jj];
      }
    }
    return(0);
  }
  nX	=ESFp1;
  neq	=2*(ESMp+nX);
  k	= neq*(neq+1)+3*ESNp1;
  if(Isize < k){
    if(Isize){
      free(ac);
    }
    Isize	= k;
    ac		= (double*)malloc(Isize*sizeof(double)+neq*sizeof(int));
    if(ac == NULL){
      printf("Allocation failure for ac in ECVirtualPlV2vdrPlV\n");
      exit(0);
    }
    f	=ac+neq*neq;
    G	=f+neq;
    G1	=G+ESNp1;
    F	=G1+ESNp1;
    index	=(int*)(F+ESNp1);
  }
  EZd_vb	+=EZgd_vb;
  am	=b_X+EZd_vb*b_X-ESsb[ESNa];
  for(i=0; i < ESNa1; i++){
    ESsb[i]	+=ESsa[i]*am;
    ESsb1a[i]	+=am;
  }

  ki	=ESNa;
  for(k=0; k < ESMp; k++){
    Fcs[k]=rcT1a[ki]+drCc[k];
    Fsn[k]=rsT1a[ki]+drCs[k];
    ki	+=ESNa1;
  }
  am	=ESsb1a[ESNa]/b_X;
  for(j=0; j < ESNp1; j++){
    ss		=EZcr2*(1.-ESsn1[j]*singtX)-2.*EZd_vb;
    s		=sqrt(ss);
    EZx2sep[j]	=s/EZd_vb;
    F[j]	=EZcr4*am*(2.*ss+EZd_vb)/ss;
    x_0		=EZx0sep[j];
    x_1		=EZcr4*EZx1sep[j];
    G[j]	=-x_0*s-EZd_vb*x_1/s;
    G1[j]	=am*((x_0-EZd_vb*x_1/ss)-x_1)/s;
  }
  EZfour(EZrvbcs,EZrvbsn,G1);
  EZrvbcs[0]	-=am*EZx0sn[0]*singtX;
  EZrvbsn[0]	=-ESsb1a[ESNa]*singtX;
  EZrvbsn[1]	-=EZcr2*am*EZx0sn[0];
  for(mi=0; mi < ESMp; mi++){
    i		=mi;
    ii		=i+ESMp;
    f[mi]	=EZrvbcs[i]+Fcs[i];
    f[ii]	=EZrvbsn[i]+Fsn[i];
    i		*=neq;
    ii		*=neq;
    for(mj=0,j=2*nX; mj < ESMp; mj++,j++){
      jj		=(mj+1)*ESMp1+mi;
      ac[i+j]		=-drCc[jj];
      ac[i+j+ESMp]	=-drSc[jj];
      ac[ii+j]		=-drCs[jj];
      ac[ii+j+ESMp]	=-drSs[jj];
    }
  }
  EZfour(Fcs,Fsn,F);
  ac[0]		=Fcs[0];
  ac[ESMp*neq]	=0.;
  for(k=1; k < nX; k++){
    ac[k]		=2.*Fcs[k];
    ac[nX+k]		=2.*Fsn[k];
    ac[ESMp*neq+k]	=0.;
    ac[ESMp*neq+k+nX]	=0.;
  }
  for(mi=1; mi < ESMp; mi++){
    i	=mi;
    ii	=mi+ESMp;
    i	*=neq;
    ii	*=neq;
    for(mj=0; mj < nX; mj++){
      ac[i+mj]		=0.;
      ac[i+mj+nX]	=0.;
      ac[ii+mj]		=0.;
      ac[ii+mj+nX]	=0.;
    }
    for(k=1,mm=mi+k; k <ESFp1 && mm < nX; k++,mm++){
      ac[i+mm]		+=Fcs[k];
      ac[i+mm+nX]	+=Fsn[k];
      ac[ii+mm]		-=Fsn[k];
      ac[ii+mm+nX]	+=Fcs[k];
    }
    for(k=0,mm=mi-k; k < ESFp1; k++,mm--){
      mj	=abs(mm);
      if(mj < nX){
	if(mm > 0){
	  ac[i+mm]	+=Fcs[k];
	  ac[i+mm+nX]	-=Fsn[k];
	  ac[ii+mm]	+=Fsn[k];
	  ac[ii+mm+nX]	+=Fcs[k];
	}
	else{
	  ac[i+mj]	+=Fcs[k];
	  ac[i+mj+nX]	+=Fsn[k];
	  ac[ii+mj]	+=Fsn[k];
	  ac[ii+mj+nX]	-=Fcs[k];
	}
      }
    }
    ac[i+nX]	=0.;
    ac[ii+nX]	=0.;
  }
  ac[nX]		=2.*am*EZx0sn[0]*EZd_vb*singtX;
  ac[ESMp*neq+nX]	=2.*ESsb1a[ESNa]*EZd_vb*singtX;

  EZfour(EZrvbcs,EZrvbsn,G);
  EZrvbcs[0]	+=ESaR0[0]-ESaRx+EZx0sn[0]*(1.-EZd_vb)*singtX;
  EZrvbsn[0]	=ESaZ0[0]-z_X0-EZd_vb*b_X*singtX;
  EZrvbsn[1]	-=EZcr2*ESsb[ESNa]*EZx0sn[0]/b_X;
  for(mi=0,k=0; mi < nX; mi++,k++){
    i		=2*ESMp+mi;
    ii		=i+nX;
    f[i]	=EZrvbcs[k]+rcT[ESNa1*k+ESNa];
    f[ii]	=EZrvbsn[k]+rsT[ESNa1*k+ESNa];
    i		*= neq;
    ii		*= neq;
    for(mj=0,j=2*nX; mj < ESMp; mj++,j++){
      if(k < ESMp1){
	jj		=(mj+1)*ESMp1+k;
	ac[i+j]		=-rCc[jj];
	ac[i+j+ESMp]	=-rSc[jj];
	ac[ii+j]	=-rCs[jj];
	ac[ii+j+ESMp]	=-rSs[jj];
      }
      else{
	ac[i+j]		=0.;
	ac[i+j+ESMp]	=0.;
	ac[ii+j]	=0.;
	ac[ii+j+ESMp]	=0.;
      }
    }
  }
  am	=EZcr4*EZd_vb;
  for(mi=0,k=0; mi < nX; mi++,k++){
    i	=2*ESMp+mi;
    ii	=i+nX;
    i	*=neq;
    ii	*=neq;
    for(mj=0; mj < nX; mj++){
      ac[i+mj]		=0.;
      ac[i+mj+nX]	=0.;
      ac[ii+mj]		=0.;
      ac[ii+mj+nX]	=0.;
    }
    ac[i+k]		=am;
    ac[ii+k+nX]		=am;
  }
  ac[2*ESMp*neq+nX]	=EZx0sn[0]*EZd_vb*EZd_vb*singtX;
  ac[(2*ESMp+nX)*neq+nX]=b_X*EZd_vb*EZd_vb*singtX;

  LUdcmp(ac,neq,index,&am);
  LUbksb(ac,neq,index,f);

  z_X2	=f[nX];
  for(j=0; j < ESNp; j++){
    x_2		=0.;
    kj		=0;
    for(k=1; k < nX; k++){
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      x_2	+=f[k]*EScs1[kj]+f[nX+k]*ESsn1[kj];
    }
    EZx2sep[j]	*=(f[0]+2.*x_2);
  }
  EZx2sep[ESNp]	=EZx2sep[0];
  for(i=0; i <ESNa1; i++){
    mj		=M0Mm*i+ESMp;
    EZgpei[mj]	=0.;
    for(k=1; k < ESMp1; k++){
      kj	=mj-k;
      jj	=mj+k;
      EZgper[jj]	=EZcr2*(EZgper[jj]+EZgper[kj]);
      EZgpei[jj]	=EZcr2*(EZgpei[jj]-EZgpei[kj]);
      EZdgper[jj]	=EZcr2*(EZdgper[jj]+EZdgper[kj]);
      EZdgpei[jj]	=EZcr2*(EZdgpei[jj]-EZdgpei[kj]);
    }
    for(k=0; k < ESMp1; k++){
      mj	=M0Mm*i+ESMp+k;
      kj	=M0Mm*i+ESMp-k;
      for(m=1,mi=2*nX; m < ESMp1; m++,mi++){
	mm	=mi+ESMp;
	j	=mj+ES2Mp1*m;
	jj	=kj+ES2Mp1*m;
	EZgper[mj]+=EZcr2*(f[mi]*(EZgper[j]+EZgper[jj])+f[mm]*(EZgpei[j]+EZgpei[jj]));
	EZgpei[mj]+=EZcr2*(f[mi]*(EZgpei[j]-EZgpei[jj])-f[mm]*(EZgper[j]-EZgper[jj]));
	EZdgper[mj]+=EZcr2*(f[mi]*(EZdgper[j]+EZdgper[jj])+f[mm]*(EZdgpei[j]+EZdgpei[jj]));
	EZdgpei[mj]+=EZcr2*(f[mi]*(EZdgpei[j]-EZdgpei[jj])-f[mm]*(EZdgper[j]-EZdgper[jj]));
      }
      EZgper[kj]	=EZgper[mj];
      EZgpei[kj]	=-EZgpei[mj];
      EZdgper[kj]	=EZdgper[mj];
      EZdgpei[kj]	=-EZdgpei[mj];
    }
  }

  for(k=0; k < ESMp1; k++){
    for(m=1,mi=2*nX; m < ESMp1; m++,mi++){
      mm	=mi+ESMp;
      j		=ESMp1*m+k;
      rCc[k]	+=f[mi]*rCc[j]+f[mm]*rSc[j];
      drCc[k]	+=f[mi]*drCc[j]+f[mm]*drSc[j];
      rCs[k]	+=f[mi]*rCs[j]+f[mm]*rSs[j];
      drCs[k]	+=f[mi]*drCs[j]+f[mm]*drSs[j];
    }
  }
  return(0);
}

int ECDeVirtualPlV2vdrPlV()
{
  if(ac != NULL){
    free(ac);
    ac	=NULL;
  }
  return(0);
}

int EcMetrTnsrX()
{
  static int i,j,k,ki,kj,iX,iO;
  static double r,x_0,dx_0,d2x_0;
  static double cs,sn,d,s,sqs,x_1,x_2,z_1;
  extern double EZdx_0x,d2x_0x,EZdx_0o,d2x_0o;

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
  for(j=0; j < ESNp1; j++){
    dx_0	=0.;
    d2x_0	=0.;
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
    }
    dx_0	*=2.*singtX;
    d2x_0	*=-2.;
    cs		=EScs1[j];
    sn		=ESsn1[j]*singtX;
    s		=1./(EZcr2*(1.-sn)-2.*EZd_vb);
    sqs		=sqrt(s);
    x_1		=3.*dx_0*(1.-sn*singtX)+(2.*d2x_0-EZx0sep[j])*cs;
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
      r		=ESaRx+(EZx0sep[j]+d*x_1)*sqs+EZx0sn[0]*(z_1-singtX);
      if(j != iX){
	x_1	=EZcr2*(dx_0*(4.*z_X2*(1.-sn)-8.*(1.+sn))
		      +(2.*d2x_0*(1.+sn)-2.*EZx2sep[j]-2.*EZx0sep[j]*z_X2
			-(2.*dx_0*cs+2.*EZx0sep[j])/(1.-sn)
			-EZcr2*EZx0sep[j]*(5.-sn)
			)*cs)/sqs;
	s	=ESgt[j]-ESgt[iX];
	x_2	=16.*EZdx_0x/sqrt(s*s-8.*d);
	ESDx[i]	+=x_1+x_2;
	ESLx[i]	+=x_1/r+x_2/ESaRx;
	ESVx[i]	+=x_1*r+x_2*ESaRx;
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
  printf("??? dgaF=%10.3e EsiAx=%d\n",ESgFPlVx,ESiAx);
  return(0);
}

int ECVirtualPlV00()
{
  static int i,j,k,ki,ji,kj;
  static int m,mm,mi,mj,ii,jj;
  static double cs,sn,d,s,ss,am,x_0,x_1,x_2;
  static int nX,neq;

  if(ESiAx == 0){
    return(0);
  }

  nX	=ESFp1;
  neq	=2*nX;

  for(j=0; j < ESNp1; j++){
    cs		=EScs1[j];
    sn		=ESsn1[j];
    ss		=EZcr2*(1.-singtX*sn)-2.*EZd_vb;
    s		=sqrt(ss);
    EZx2sep[j]	=s/EZd_vb;
    x_0		=EZx0sep[j];
    x_1		=EZcr4*EZx1sep[j];
    G[j]	=-x_0*s-EZd_vb*x_1/s;
  }
  EZfour(EZrvbcs,EZrvbsn,G);
  z_X2	=(ESaZ0[0]+rsT[ESNa]-z_X0-EZd_vb*b_X*singtX)/b_X;
  am		=4./EZd_vb;
  EZrvbcs[0]	+=ESaR0[0]+rcT[ESNa]-ESaRx-EZx0sn[0]*(1.-EZd_vb-z_X2)*singtX;
  f[0]		=EZrvbcs[0]*am;
  z_X2		/=EZd_vb*EZd_vb;
  EZrvbsn[1]	+=EZcr2*ESsb[ESNa]*EZx0sn[0]/b_X;
  for(k=1; k < nX; k++){
    f[k]	=(EZrvbcs[k]+rcT[ESNa1*k+ESNa])*am;
    f[k+nX]	=(EZrvbsn[k]+rsT[ESNa1*k+ESNa])*am;
  }

  for(j=0; j < ESNp; j++){
    x_2		=0.;
    kj		=0;
    for(k=1; k < nX; k++){
      kj	+=j;
      if(kj >= ESNp){
	kj	-= ESNp;
      }
      x_2	+=f[k]*EScs1[kj]+f[nX+k]*ESsn1[kj];
    }
    EZx2sep[j]	*=(f[0]+2.*x_2);
  }
  EZx2sep[ESNp]	=EZx2sep[0];
  return(0);
}

ECvdrzCorrectionX()
{
  static int i,j,k,ki,ji,kj;
  static int m,mm,mi,mj,ii,jj;
  static double cs,sn,d,s,ss,am,x_0,x_1,x_2;
  static int nX,neq;

  nX	=ESFp1;
  neq	=2*nX;
  EZd_vb	=(ECb[ESNa]-b_X)/b_X;

  for(j=0; j < ESNp1; j++){
    ss		=EZcr2*(1.-ESsn1[j]*singtX)-2.*EZd_vb;
    s		=sqrt(ss);
    EZx2sep[j]	=s/EZd_vb;
    G[j]	=-EZx0sep[j]*s-EZd_vb*EZcr4*EZx1sep[j]/s;
  }
  EZfour(EZrvbcs,EZrvbsn,G);
  z_X2	=(Z0+EZz0[ESNa]-z_X0-EZd_vb*b_X*singtX)/b_X;
  am		=4./EZd_vb;
  EZrvbcs[0]	+=R0+EZrcs[ESNa]-ESaRx-EZx0sn[0]*(1.-EZd_vb-z_X2)*singtX;
  f[0]		=EZrvbcs[0]*am;
  z_X2		/=EZd_vb*EZd_vb;
  EZrvbsn[1]	-=EZcr2*ECb[ESNa]*EZx0sn[0]/b_X;
  for(k=1; k < nX; k++){
    f[k]	=(EZrvbcs[k]+EZrcs[ESNa1*k+ESNa])*am;
    f[k+nX]	=(EZrvbsn[k]+EZrsn[ESNa1*k+ESNa])*am;
  }
  for(j=0; j < ESNp; j++){
    x_2		=0.;
    kj		=0;
    for(k=1; k < nX; k++){
      kj	+=j;
      if(kj >= ESNp){
	kj	-=ESNp;
      }
      x_2	+=f[k]*EScs1[kj]+f[nX+k]*ESsn1[kj];
    }
    EZx2sep[j]	*=(f[0]+2.*x_2);
  }
  EZx2sep[ESNp]	=EZx2sep[0];
  return(0);
}


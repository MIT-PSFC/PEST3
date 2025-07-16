#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef mzl_SPLINES

static double ccr2= 0.5,ccr4= 0.25,ccr5= 0.2,ccr3= 0.3333333333333,
ccr6= 0.16666666666666,ccr12= 0.083333333333333,ccr7= 0.1428571428571,
ccr30= 0.03333333333333;

static int Nv=0,Ispl=0,Ispl1=1;
static double rhSpl,hSpl,hhSpl,xSpl,xxSpl,XSpl,XXSpl;
static double *v0=NULL;
static double *v1=NULL;
static double *v2=NULL;
static double *v4=NULL;
static double *vp=NULL,*d2f,*u;
static double *gl=NULL,*v3;

int DeInitSplines()
{
  if(vp != NULL){
    free(vp);
    vp	=NULL;
  }
  if(v0 != NULL){
    free(v0);
    v0	=NULL;
  }
  if(v1 != NULL){
    free(v1);
    v1	=NULL;
  }
  if(v2 != NULL){
    free(v2);
    v2	=NULL;
  }
  if(gl != NULL){
    free(gl);
    gl	=NULL;
  }
  if(v4 != NULL){
    free(v4);
    v4	=NULL;
  }
  return(0);
}

int spl(double *g, double *d2g, double *dg0, double *dg1, double *x, int *n)
{
  static int i,ii,i1,Ind0,Isize= 0;
  static double h,t1,t0,r0,r1,w,w0,w1,s;
  if(Isize < (*n)+1 || v0 == NULL){
    if(v0 != NULL){
      free(v0);
    }
    Isize	=(*n)+1;
    v0	=(double*)malloc(Isize*sizeof(double));
  }
  h	=x[1]-x[0];
  r0	=(g[1]-g[0])/h;
  w0	=h*ccr3;
  t0	=ccr2*w0;
  h	=x[2]-x[1];
  r1	=(g[2]-g[1])/h;
  w1	=h*ccr3;
  t1	=ccr2*w1;
  if(dg0 != NULL){
    v0[0]	=-ccr2;
    d2g[0]	=(r0-(*dg0))/w0;
    w		=1./(w0+w1+t0*v0[0]);
    v0[1]	=-t1*w;
    d2g[1]	=(r1-r0-t0*d2g[0])*w;
  }
  else{
    v0[0]	=2.;
    d2g[0]	=0.;
    w		=1./(w0+w1+t0*v0[0]);
    v0[1]	=-(t1-t0)*w;
    d2g[1]	=(r1-r0-t0*d2g[0])*w;
  }
  for(i=2,ii=3,i1=1; i <*n ;i++,ii++,i1++){
    h		=x[ii]-x[i];
    w0		=w1;
    w1		=h*ccr3;
    t0		=t1;
    t1		=ccr2*w1;
    r0		=r1;
    r1		=(g[ii]-g[i])/h;
    w		=1./(w0+w1+t0*v0[i1]);
    v0[i]	=-t1*w;
    d2g[i]	=(r1-r0-t0*d2g[i1])*w;
  }
  ii	=*n;
  i	=ii-1;
  i1	=i-1;
  if(dg1 == NULL){
    t1		=v0[i1]-2.;
    d2g[ii]	=-(d2g[i1]+t1*d2g[i])/(t1*v0[i]+1.);
  }
  else{
    d2g[ii]	=((*dg1)-r1-t1*d2g[i])/(t1*v0[i]+w1);
  }
  while(ii > 0){
    d2g[i]	+=v0[i]*d2g[ii];
    i--;
    ii--;
  }
  if(dg0 == NULL){
    d2g[0]	-=d2g[2];
  }
  return(0);
}

int spl2(double *g, double *d2g, double *d2g0, double *d2g1, double *x, int *n)
{
  static int i,ii,i1,Ind0,Isize= 0;
  static double h,t1,t0,r0,r1,w,w0,w1,s;
  if(Isize < (*n)+1 || v0 == NULL){
    if(v0 != NULL)
      free(v0);
    Isize	=(*n)+1;
    v0	=(double*)malloc(Isize*sizeof(double));
  }
  h	=x[1]-x[0];
  r0	=(g[1]-g[0])/h;
  w0	=h*ccr3;
  t0	=ccr2*w0;
  h	=x[2]-x[1];
  r1	=(g[2]-g[1])/h;
  w1	=h*ccr3;
  t1	=ccr2*w1;
  if(d2g0 != NULL){
    v0[0]	=0.;
    d2g[0]	=*d2g0;
    w		=1./(w0+w1+t0*v0[0]);
    v0[1]	=-t1*w;
    d2g[1]	=(r1-r0-t0*d2g[0])*w;
  }
  else{
    v0[0]	=2.;
    d2g[0]	=0.;
    w		=1./(w0+w1+t0*v0[0]);
    v0[1]	=-(t1-t0)*w;
    d2g[1]	=(r1-r0-t0*d2g[0])*w;
  }
  for(i=2,ii=3,i1=1; i <*n ;i++,ii++,i1++){
    h		=x[ii]-x[i];
    w0		=w1;
    w1		=h*ccr3;
    t0		=t1;
    t1		=ccr2*w1;
    r0		=r1;
    r1		=(g[ii]-g[i])/h;
    w		=1./(w0+w1+t0*v0[i1]);
    v0[i]	=-t1*w;
    d2g[i]	=(r1-r0-t0*d2g[i1])*w;
  }
  ii	=*n;
  i	=ii-1;
  i1	=i-1;
  if(d2g1 == NULL){
    t1		=v0[i1]-2.;
    d2g[ii]	=-(d2g[i1]+t1*d2g[i])/(t1*v0[i]+1.);
  }
  else{
    d2g[ii]	=*d2g1;
  }
  while(ii > 0){
    d2g[i]	+=v0[i]*d2g[ii];
    i--;
    ii--;
  }
  if(d2g0 == NULL){
    d2g[0]	-=d2g[2];
  }

  return(0);
}

int EZspl_p(double *g,double *d2g,double *x,int n)
{
  int i;
  double a1,a2,b,c1,c2,h1,h2,s;
  if(Nv < n+1 || vp == NULL){
    if(vp != NULL){
      free(vp);
    }
    Nv	=n+1;
    vp	=(double*)malloc(3*Nv*sizeof(double));
    d2f	=vp+Nv;
    u	=d2f+Nv;
  }
  d2g[0]=0.;
  d2f[0]=1.;
  vp[0]	=0;
  u[0]	=0;
  h2	=(x[1]-x[0])*ccr3;
  a2	=0.5*h2;
  c2	=(g[1]-g[0])/(x[1]-x[0]);
  for(i=1; i<n; i++){
    h1	=h2;
    h2	=(x[i+1]-x[i])*ccr3;
    a1	=a2;
    a2	=0.5*h2;
    b	=h1+h2;
    c1	=c2;
    c2	=(g[i+1]-g[i])/(x[i+1]-x[i]); 
    s		=1./(b+a1*vp[i-1]);
    vp[i]	=-a2*s;
    d2g[i]	=(c2-c1-a1*d2g[i-1])*s;
    s		=1./(b+a1*u[i-1]);
    u[i]	=-a2*s;
    d2f[i]	=-a1*d2f[i-1]*s; 
  }
  for(i=n,d2g[i]=d2g[0], d2f[i]=d2f[0]; i>0; i--){
    d2g[i-1]	+=vp[i-1]*d2g[i];
    d2f[i-1]	+=u[i-1]*d2f[i];
  }
  h1	=x[1]-x[0];
  a1	=(g[1]-g[0])/h1-(0.5*d2g[1]+d2g[0])*h1*ccr3;
  c1	=-(0.5*d2f[1]+d2f[0])*h1*ccr3;
  i	=n;
  h2	=x[i]-x[i-1];
  a2	=(g[i]-g[i-1])/h2+(d2g[i]+0.5*d2g[i-1])*h2*ccr3;
  c2	=(d2f[i]+0.5*d2f[i-1])*h2*ccr3;
  b	=(a2-a1)/(c1-c2);
  for(i=0; i<=n; i++){
    d2g[i]	+=b*d2f[i]; 
  }
  return(0);
}

int splr0(f,df,t0,g,d2g,x,n)
     double*f,*df,*t0,*g,*d2g,*x;
     int*n;
{
  int i,ii,j;
  double a,b,aa,bb,h,t;
  t	=*t0;
  ii	=(*n);
  i	=0;
  if(t > x[ii]) t=x[ii];
  if(t < x[i])  t=x[i];
  while(ii-i > 1){
    j=(ii+i)/2;
    if(t >= x[j]) i=j;
    else ii=j;
  }
  h=x[ii]-x[i];
  a=(x[ii]-t)/h;
  aa=a*a;
  b=1.-a;
  bb=b*b;
  *f=a*g[i]+b*g[ii]+(a*(aa-1.)*d2g[i]+b*(bb-1.)*d2g[ii])*h*h*ccr6;
  if(df != NULL)
    *df=(g[ii]-g[i])/h+((3.*bb-1.)*d2g[ii]-(3.*aa-1.)*d2g[i])*h*ccr6;
  return(0);
}

int SetIspl(double *x,double s,int n)
{
  if(Ispl1 >= n){
    Ispl	=n-1;
  }
  while(x[Ispl] > s){
    Ispl--;
  }
  while(x[Ispl+1] < s){
    Ispl++;
  }
  Ispl1	=Ispl+1;
  hSpl	=x[Ispl1]-x[Ispl];
  rhSpl	=1/hSpl;
  XSpl	=(x[Ispl1]-s)*rhSpl;
  XXSpl	=XSpl*XSpl;
  xSpl	=1.-XSpl;
  xxSpl	=xSpl*xSpl;
  hhSpl	=hSpl*hSpl*ccr6;
  return(0);
}

int splR(double*f, double*df, double*g, double*d2g)
{
  *f=XSpl*g[Ispl]+xSpl*g[Ispl1]
    +(XSpl*(XXSpl-1.)*d2g[Ispl]+xSpl*(xxSpl-1.)*d2g[Ispl1])*hhSpl;
  if(df!=NULL){
    *df=(g[Ispl1]-g[Ispl]
	 +((3.*xxSpl-1.)*d2g[Ispl1]-(3.*XXSpl-1.)*d2g[Ispl])*hhSpl)*rhSpl;
  }
  return(0);
}

int splr1(double*f,double*df,double*t0,double*g,double*d2g,double*x
	  ,int*n,int*it)
{
  int i,ii;
  double a,b,aa,bb,h,t;
  t	=*t0;
  i	=*it;
  ii	=i+1;
  h	=x[ii]-x[i];
  a	=(x[ii]-t)/h;
  aa	=a*a;
  b	=1.-a;
  bb	=b*b;
  *f	=a*g[i]+b*g[ii]+(a*(aa-1.)*d2g[i]+b*(bb-1.)*d2g[ii])*h*h*ccr6;
  if(df != NULL)
    *df	=(g[ii]-g[i])/h+((3.*bb-1.)*d2g[ii]-(3.*aa-1.)*d2g[i])*h*ccr6;
  return(0);
}

int splr2(double*f,double*df,double*d2f,double*t0,double*g,double*d2g,double*x
	  ,int*n,int*it)
{
  int i,ii;
  double a,b,aa,bb,h,t;
  t	=*t0;
  i	=*it;
  ii	=i+1;
  h	=x[ii]-x[i];
  a	=(x[ii]-t)/h;
  aa	=a*a;
  b	=1.-a;
  bb	=b*b;
  *f	=a*g[i]+b*g[ii]+(a*(aa-1.)*d2g[i]+b*(bb-1.)*d2g[ii])*h*h*ccr6;
  if(df != NULL)
    *df	=(g[ii]-g[i])/h+((3.*bb-1.)*d2g[ii]-(3.*aa-1.)*d2g[i])*h*ccr6;
  if(d2f != NULL)
    *d2f=a*d2g[i]+b*d2g[ii];
  return(0);
}

int spl0(g,d2g,x,n)
     double g[],d2g[],x[];
     int*n;
{
  static int i,Isize= 0;
  static double a1,a2,b,c1,c2,h1,h2,s;
  if(Isize<(*n)+1 || v1 == NULL){
    if(v1 != NULL)
      free(v1);
    Isize= (*n)+1;
    v1=(double*)malloc(Isize*sizeof(double));
  }
  d2g[0]= 0.;
  v1[0]= 0;
  h2= (x[1]-x[0])*ccr3;
  a2= 0.5*h2;
  c2= (g[1]-g[0])/(x[1]-x[0]);
  for(i= 1;i<(*n);i++)
    {
      h1= h2;
      h2= (x[i+1]-x[i])*ccr3;
      a1= a2;
      a2= 0.5*h2;
      b= h1+h2;
      c1= c2;
      c2= (g[i+1]-g[i])/(x[i+1]-x[i]);
      s= 1./(b+a1*v1[i-1]);
      v1[i]= -a2*s;
      d2g[i]= (c2-c1-a1*d2g[i-1])*s;
    }
  for(i= (*n),d2g[i]= 0.;i>0;i--){
    d2g[i-1]+= v1[i-1]*d2g[i];
  }
  return(0);
}

int spl_f0(g,d2g,p,f,x,n,dg0,dg1)
     double*g,*d2g,*p,*f,*x,*dg0,*dg1;
     int*n;
{
  static int i,ii,i1,Ind0,Isize= 0;
  static double h,t1,t0,r0,r1,w,w0,w1,s;
  if(Isize<(*n)+1 || v2 == NULL){
    if(v2 != NULL)	
      free(v2);
    Isize= (*n)+1;
    v2= (double*)malloc(Isize*sizeof(double));
  }
  if(p!=NULL&&p[0]/p[*n]<1e-6)Ind0= 0;
  else Ind0= 1;
  h= x[2]-x[1];
  r1= (g[2]-g[1])/h;
  w1= h*ccr3;
  t1= ccr2*w1;
  h= x[1]-x[0];
  w0= h*ccr3;
  t0= ccr2*w0;
  w= w0+w1;
  if(p==NULL){
    g[1]= f[1];
    g[2]= f[2];
  }
  else{
    g[1]= f[1]/p[1];
    g[2]= f[2]/p[2];
  }
  if(Ind0){
    if(p==NULL)g[0]= f[0];
    else g[0]= f[0]/p[0];
    r0= (g[1]-g[0])/h;
    v2[0]= -ccr2;
    d2g[0]= (r0-(*dg0))/w0;
  }
  else{
    r0= *dg0;
    w+= t0;
    t0+= w0;
    s= 1./(t0-t1);
    v2[0]= -s*(w+w1);
    d2g[0]= s*(r1-r0);
  }
  s= 1./(w+t0*v2[0]);
  v2[1]= -s*t1;
  d2g[1]= s*(r1-r0-t0*d2g[0]);
  for(i= 2,ii= 3,i1= 1;i<(*n);i++,ii++,i1++){
    if(p==NULL)g[ii]= f[ii];
    else g[ii]= f[ii]/p[ii];
    h= x[ii]-x[i];
    w0= w1;
    w1= h*ccr3;
    t0= t1;
    t1= ccr2*w1;
    w= w0+w1;
    r0= r1;
    r1= (g[ii]-g[i])/h;
    if(ii==*n)t0-= t1;
    s= 1./(w+t0*v2[i1]);
    v2[i]= -s*t1;
    d2g[i]= (r1-r0-t0*d2g[i1])*s;
  }
  ii= *n;
  i= ii-1;
  if(dg1==NULL)d2g[ii]= d2g[i]/(ccr2-v2[i]);
  else d2g[ii]= ((*dg1)-r1-t1*d2g[i])/(t1*v2[i]+w1);
  while(ii>0){
    d2g[i]+= v2[i]*d2g[ii];
    i--;
    ii--;
  }
  if(dg1==NULL)d2g[*n]-= d2g[(*n)-2];
  if(Ind0==0){
    h= x[1]-x[0];
    g[0]= -h*(h*(ccr3*d2g[0]+ccr6*d2g[1])+(*dg0))-g[1];
  }
  return(0);
}

int EZinv3x3(w,rw)
     double w[],rw[];
{
  double d;
  rw[0]= w[3]*w[5]-w[4]*w[4];
  rw[1]= w[4]*w[2]-w[1]*w[5];
  rw[2]= w[1]*w[4]-w[3]*w[2];
  d= 1./(w[0]*rw[0]+w[1]*rw[1]+w[2]*rw[2]);
  rw[0]*= d;
  rw[1]*= d;
  rw[2]*= d;
  rw[3]= (w[0]*w[5]-w[2]*w[2])*d;
  rw[4]= (w[1]*w[2]-w[0]*w[4])*d;
  rw[5]= (w[3]*w[0]-w[1]*w[1])*d;
  return(0);
}

int EZf2spl(double*g,double*d2g,double*dg0,double*dg1,double*f,double*p,double*x,
int*nx,double*EZga0,double*ga_1,double*ga_2,double*ga_3)
{
  static int i,ii,k,Ind0,Isize= 0;
  static double EZga1,EZga2,EZga3,g0;
  static double s,sc0,sc1,sc11,sd0,sd1,sd11,rh1,rh11,sq;
  static double r1,EZr2,r3;
  static double w[6],rw[6];
  static double t11,t13,t22,t23,t31,t32;
  static double bc00= 0.,bc01= 0.,bd00= 0.,bd01= 0.,be0= 1.,bf0= 0.;
  static double bcn0= 0.,bcn1= 0.,bdn0= 0.,bdn1= 0.,ben= 1.,bfn= 0.;

  if(Isize<(*nx)+1 || gl == NULL){
    if(gl != NULL)
      free(gl);
    Isize= (*nx)+1;
    gl= (double*)malloc(10*Isize*sizeof(double));
    v3= gl+Isize;
  }
  i= *nx;
  rh1= x[i]-x[i-1];
  if(dg1 != NULL){
    bfn= (*dg1);
    bcn0= 1./rh1;
    bcn1= -bcn0;
    bdn0= ccr3*rh1;
    bdn1= ccr6*rh1;
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
  if(p != NULL){
    bfn-= bcn0*f[i]/p[i];
  }
  else{
    bfn-= bcn0*f[i];
  }
  bcn0= 0.;
  rh1= x[1]-x[0];
  if(dg0 != NULL){
    bf0= (*dg0);
    bc01= 1./rh1;
    bc00= -bc01;
    bd00= -ccr3*rh1;
    bd01= -ccr6*rh1;
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
  s= x[*nx]-x[0];
  s*= s;
  EZga1= (*ga_1)*s;
  EZga2= (*ga_2)*s*s;
  EZga3= (*ga_3)*s*s*s;
  for(i=0; i <6; i++){
    w[i]=0.;
  }
  if(p != NULL && fabs(p[0]/p[*nx]) < 1e-6){
    Ind0= 0.;
  }
  else{
    if(p!=NULL){
      g0= f[0]/p[0];
    }
    else{
      g0= f[0];
    }
    Ind0= 1;
  }
  r1= 0.;
  EZr2= 0.;
  r3= 0.;
  for(i= 0,ii= 1,sc1= 0.,sd1= 0.,rh1= 0.;i<=*nx;i++,ii++){
    k= 9*i;
    rh11= rh1;
    sc11= sc1;
    sd11= sd1;
    if(i<*nx){
      rh1= x[ii]-x[i];
      sd1= ccr6*rh1;
      sc1= 1./rh1;
      rh1= sc1*sc1;
    }
    if(i==0){
      sq= ccr2*(x[ii]-x[i]);
      if(p!=NULL)w[0]= EZga1*sc1+sq*((*EZga0)+p[i]*p[i]);
      else w[0]= EZga1*sc1+sq*((*EZga0)+1.);
      w[2]= bc00;
      w[3]= sq*EZga2+EZga3*sc1;
      w[4]= bd00;
      w[5]= be0;
      r3= bf0;
      if(Ind0){
	w[2]= 0.;
	r3-= bc00*g0;
      }
    }
    else if(i<*nx){
      sq= ccr2*(x[ii]-x[i-1]);
      if(p!=NULL)w[0]+= EZga1*(sc1+sc11)+sq*((*EZga0)+p[i]*p[i]);
      else w[0]+= EZga1*(sc1+sc11)+sq*((*EZga0)+1.);
      w[2]+= sc1+sc11;
      w[3]+= sq*EZga2+EZga3*(sc1+sc11);
      w[4]+= ccr3*(x[ii]-x[i-1]);
      if(i==1&&Ind0){
	r1+= EZga1*sc11*g0;
	r3+= sc11*g0;
      }
    }
    else{
      sq= ccr2*(x[i]-x[i-1]);
      if(p!=NULL)w[0]+= sq*p[i]*p[i];
      else w[0]+= sq;
      w[3]+= sq*EZga2+EZga3*sc11;
      w[4]+= bdn0;
      w[5]+= ben;
      r3+= bfn;
    }
    if(p!=NULL)r1+= sq*p[i]*f[i];
    else r1+= sq*f[i];
    if(i==0){
      t11= -EZga1*sc1;
      t13= -sc1;
      t23= sd1;
      t31= bc01;
      t32= bd01;
      if(Ind0){
	if(p!=NULL)w[0]= sq*p[i]*p[i];
	else w[0]= sq;
	t13= 0.;
	t11= 0.;
      }
    }
    else if(ii<*nx){
      t11= -EZga1*sc1;
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
      if(i<*nx){
	if(p!=NULL){
	  r1+= EZga1*sc1*f[*nx]/p[*nx];
	  r3+= sc1*f[*nx]/p[*nx];
	}
	else{
	  r1+= EZga1*sc1*f[*nx];
	  r3+= sc1*f[*nx];
	}
      }
    }
    t22= -EZga3*sc1;

    EZinv3x3(w,rw);

    g[i]= rw[0]*r1+rw[1]*EZr2+rw[2]*r3;
    d2g[i]= rw[1]*r1+rw[3]*EZr2+rw[4]*r3;
    gl[i]= rw[2]*r1+rw[4]*EZr2+rw[5]*r3;
    if(i<*nx){

      v3[k+0]= -(rw[0]*t11+rw[2]*t31);
      v3[k+3]= -(rw[1]*t11+rw[4]*t31);
      v3[k+6]= -(rw[2]*t11+rw[5]*t31);
      v3[k+1]= -(rw[1]*t22+rw[2]*t32);
      v3[k+4]= -(rw[3]*t22+rw[4]*t32);
      v3[k+7]= -(rw[4]*t22+rw[5]*t32);
      v3[k+2]= -(rw[0]*t13+rw[1]*t23);
      v3[k+5]= -(rw[1]*t13+rw[3]*t23);
      v3[k+8]= -(rw[2]*t13+rw[4]*t23);

      w[0]= t11*v3[k+0]+t31*v3[k+6];
      w[1]= t11*v3[k+1]+t31*v3[k+7];
      w[2]= t11*v3[k+2]+t31*v3[k+8];
      w[3]= t22*v3[k+4]+t32*v3[k+7];
      w[4]= t13*v3[k+1]+t23*v3[k+4];
      w[5]= t13*v3[k+2]+t23*v3[k+5];


      r1= -(t11*g[i]+t31*gl[i]);
      EZr2= -(t22*d2g[i]+t32*gl[i]);
      r3= -(t13*g[i]+t23*d2g[i]);
    }
  }
  for(ii= *nx,i= ii-1;ii>0;ii--,i--){
    k= 9*i;
    g[i]= g[i]+v3[k+0]*g[ii]+v3[k+1]*d2g[ii]+v3[k+2]*gl[ii];
    d2g[i]= d2g[i]+v3[k+3]*g[ii]+v3[k+4]*d2g[ii]+v3[k+5]*gl[ii];
    gl[i]= gl[i]+v3[k+6]*g[ii]+v3[k+7]*d2g[ii]+v3[k+8]*gl[ii];
  }
  return(0);
}

int EZf2spl2(double*g,double*d2g
	   ,double*d2g0,double*d2g1
	   ,double*f,double*p,double*x
	   ,int*nx,double*EZga0,double*ga_1,double*ga_2,double*ga_3)
{
  static int i,ii,k,Ind0,Isize= 0;
  static double EZga1,EZga2,EZga3,g0;
  static double s,sc0,sc1,sc11,sd0,sd1,sd11,rh1,rh11,sq;
  static double r1,EZr2,r3;
  static double w[6],rw[6];
  static double t11,t13,t22,t23,t31,t32;
  static double bc00= 0.,bc01= 0.,bd00= 0.,bd01= 0.,be0= 1.,bf0= 0.;
  static double bcn0= 0.,bcn1= 0.,bdn0= 0.,bdn1= 0.,ben= 1.,bfn= 0.;

  if(Isize<(*nx)+1 || gl == NULL){
    if(gl != NULL)
      free(gl);
    Isize= (*nx)+1;
    gl= (double*)malloc(10*Isize*sizeof(double));
    v3= gl+Isize;
  }
  i= *nx;
  rh1= x[i]-x[i-1];
  if(d2g1 != NULL){
    bfn		=*d2g1;
    bcn0	=0.;
    bcn1	=0.;
    bdn0	=1.;
    bdn1	=0.;
    ben		=0.;
  }
  else{
    bfn		=0.;
    bcn0	=0.;
    bcn1	=0.;
    bdn0	=0.;
    bdn1	=0.;
    ben		=1.;
  }
  if(p != NULL){
    bfn-= bcn0*f[i]/p[i];
  }
  else{
    bfn-= bcn0*f[i];
  }
  bcn0= 0.;
  rh1= x[1]-x[0];
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
  s= x[*nx]-x[0];
  s*= s;
  EZga1= (*ga_1)*s;
  EZga2= (*ga_2)*s*s;
  EZga3= (*ga_3)*s*s*s;
  for(i=0; i <6; i++){
    w[i]=0.;
  }
  if(p != NULL && fabs(p[0]/p[*nx]) < 1e-6){
    Ind0= 0.;
  }
  else{
    if(p!=NULL){
      g0= f[0]/p[0];
    }
    else{
      g0= f[0];
    }
    Ind0= 1;
  }
  r1= 0.;
  EZr2= 0.;
  r3= 0.;
  for(i= 0,ii= 1,sc1= 0.,sd1= 0.,rh1= 0.;i<=*nx;i++,ii++){
    k= 9*i;
    rh11= rh1;
    sc11= sc1;
    sd11= sd1;
    if(i<*nx){
      rh1= x[ii]-x[i];
      sd1= ccr6*rh1;
      sc1= 1./rh1;
      rh1= sc1*sc1;
    }
    if(i==0){
      sq= ccr2*(x[ii]-x[i]);
      if(p!=NULL)w[0]= EZga1*sc1+sq*((*EZga0)+p[i]*p[i]);
      else w[0]= EZga1*sc1+sq*((*EZga0)+1.);
      w[2]= bc00;
      w[3]= sq*EZga2+EZga3*sc1;
      w[4]= bd00;
      w[5]= be0;
      r3= bf0;
      if(Ind0){
	w[2]= 0.;
	r3-= bc00*g0;
      }
    }
    else if(i<*nx){
      sq= ccr2*(x[ii]-x[i-1]);
      if(p!=NULL)w[0]+= EZga1*(sc1+sc11)+sq*((*EZga0)+p[i]*p[i]);
      else w[0]+= EZga1*(sc1+sc11)+sq*((*EZga0)+1.);
      w[2]+= sc1+sc11;
      w[3]+= sq*EZga2+EZga3*(sc1+sc11);
      w[4]+= ccr3*(x[ii]-x[i-1]);
      if(i==1&&Ind0){
	r1+= EZga1*sc11*g0;
	r3+= sc11*g0;
      }
    }
    else{
      sq= ccr2*(x[i]-x[i-1]);
      if(p!=NULL)w[0]+= sq*p[i]*p[i];
      else w[0]+= sq;
      w[3]+= sq*EZga2+EZga3*sc11;
      w[4]+= bdn0;
      w[5]+= ben;
      r3+= bfn;
    }
    if(p!=NULL)r1+= sq*p[i]*f[i];
    else r1+= sq*f[i];
    if(i==0){
      t11= -EZga1*sc1;
      t13= -sc1;
      t23= sd1;
      t31= bc01;
      t32= bd01;
      if(Ind0){
	if(p!=NULL)w[0]= sq*p[i]*p[i];
	else w[0]= sq;
	t13= 0.;
	t11= 0.;
      }
    }
    else if(ii<*nx){
      t11= -EZga1*sc1;
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
      if(i<*nx){
	if(p!=NULL){
	  r1+= EZga1*sc1*f[*nx]/p[*nx];
	  r3+= sc1*f[*nx]/p[*nx];
	}
	else{
	  r1+= EZga1*sc1*f[*nx];
	  r3+= sc1*f[*nx];
	}
      }
    }
    t22= -EZga3*sc1;

    EZinv3x3(w,rw);

    g[i]= rw[0]*r1+rw[1]*EZr2+rw[2]*r3;
    d2g[i]= rw[1]*r1+rw[3]*EZr2+rw[4]*r3;
    gl[i]= rw[2]*r1+rw[4]*EZr2+rw[5]*r3;
    if(i<*nx){

      v3[k+0]= -(rw[0]*t11+rw[2]*t31);
      v3[k+3]= -(rw[1]*t11+rw[4]*t31);
      v3[k+6]= -(rw[2]*t11+rw[5]*t31);
      v3[k+1]= -(rw[1]*t22+rw[2]*t32);
      v3[k+4]= -(rw[3]*t22+rw[4]*t32);
      v3[k+7]= -(rw[4]*t22+rw[5]*t32);
      v3[k+2]= -(rw[0]*t13+rw[1]*t23);
      v3[k+5]= -(rw[1]*t13+rw[3]*t23);
      v3[k+8]= -(rw[2]*t13+rw[4]*t23);

      w[0]= t11*v3[k+0]+t31*v3[k+6];
      w[1]= t11*v3[k+1]+t31*v3[k+7];
      w[2]= t11*v3[k+2]+t31*v3[k+8];
      w[3]= t22*v3[k+4]+t32*v3[k+7];
      w[4]= t13*v3[k+1]+t23*v3[k+4];
      w[5]= t13*v3[k+2]+t23*v3[k+5];


      r1= -(t11*g[i]+t31*gl[i]);
      EZr2= -(t22*d2g[i]+t32*gl[i]);
      r3= -(t13*g[i]+t23*d2g[i]);
    }
  }
  for(ii= *nx,i= ii-1;ii>0;ii--,i--){
    k= 9*i;
    g[i]= g[i]+v3[k+0]*g[ii]+v3[k+1]*d2g[ii]+v3[k+2]*gl[ii];
    d2g[i]= d2g[i]+v3[k+3]*g[ii]+v3[k+4]*d2g[ii]+v3[k+5]*gl[ii];
    gl[i]= gl[i]+v3[k+6]*g[ii]+v3[k+7]*d2g[ii]+v3[k+8]*gl[ii];
  }
  return(0);
}

int spl_d(g,d2g,f,x,nx,dg1)
     double*g,*d2g,*f,*x,*dg1;
     int*nx;
{
  static int i,ii,i1,k,Isize= 0;
  static double ss,sc1,sc11,sd1,sd11,rh1,rh11;
  static double w[4],rw[4];
  static double r1,EZr2,cr360;
  static double t0,t1,t2,t3,s0,s1,s2,s3;
  if(Isize<(*nx)+1 || v4 == NULL){
    if(v4 != NULL)
      free(v4);
    Isize= (*nx)+1;
    v4= (double*)malloc(4*Isize*sizeof(double));
    cr360= 1./360.;
  }
  for(i= 0;i<4;i++)w[i]= 0;
  for(i= 0,ii= 1,i1= -1,sc1= 0.,sd1= 0.,rh1= 0.;i<=*nx;i++,ii++,i1++){
    rh11= rh1;
    sc11= sc1;
    sd11= sd1;
    if(i<*nx){
      rh1= x[ii]-x[i];
      sd1= ccr6*rh1;
      sc1= 1./rh1;
    }
    if(i==0){
      ss= -cr360*rh1*rh1;
      w[0]= 1.;
      w[1]= 0.;
      w[2]= sc1;
      w[3]= 2.*sd1;
      t0= 0.;
      t1= 0.;
      t2= -sc1;
      t3= sd1;
      r1= 2.*f[0];
      EZr2= 0.;
      s0= ccr6*(x[ii]+2.*x[i]);
      s1= ss*(7.*x[ii]+8.*x[i]);
      s2= -sc1;
      s3= sd1;
    }
    else if(i<*nx){
      ss= -cr360*rh11*rh11;
      w[0]+= ccr6*(2.*x[i]+x[i1]);
      w[1]+= ss*(8.*x[i]+7.*x[i1]);
      w[2]+= sc1+sc11;
      w[3]+= ccr3*(x[ii]-x[i1]);
      t0= 0.;
      t1= 0.;
      t2= -sc1;
      t3= sd1;
      r1+= sc11*(x[i]*x[i]*f[i]-x[i1]*x[i1]*f[i1]);
      ss= -cr360*rh1*rh1;
      s0= ccr6*(x[ii]+2.*x[i]);
      s1= ss*(7.*x[ii]+8.*x[i]);
      if(ii==*nx&&dg1==NULL){
	k= 4*i1;
	s2= v4[k+2];
	s3= -2.+v4[k+3];
      }
      else{
	s2= -sc1;
	s3= sd1;
      }
    }
    else{
      ss= -cr360*rh11*rh11;
      w[0]+= ccr6*(2.*x[i]+x[i1]);
      w[1]+= ss*(8.*x[i]+7.*x[i1]);
      w[2]+= 0.;
      w[3]+= 1.;
      r1+= sc11*(x[i]*x[i]*f[i]-x[i1]*x[i1]*f[i1]);
    }
    EZinv2x2(w,rw);
    g[i]= rw[0]*r1+rw[1]*EZr2;
    d2g[i]= rw[2]*r1+rw[3]*EZr2;
    k= 4*i;
    if(i<*nx){
      v4[k+0]= -(rw[0]*t0+rw[1]*t2);
      v4[k+1]= -(rw[0]*t1+rw[1]*t3);
      v4[k+2]= -(rw[2]*t0+rw[3]*t2);
      v4[k+3]= -(rw[2]*t1+rw[3]*t3);
      w[0]= s0*v4[k+0]+s1*v4[k+2];
      w[1]= s0*v4[k+1]+s1*v4[k+3];
      w[2]= s2*v4[k+0]+s3*v4[k+2];
      w[3]= s2*v4[k+1]+s3*v4[k+3];
      r1= -(s0*g[i]+s1*d2g[i]);
      EZr2= -(s2*g[i]+s3*d2g[i]);
      if(ii==*nx){
	if(dg1==NULL)EZr2-= d2g[i1];
	else EZr2-= *dg1;
      }
    }
  }
  for(ii= *nx,i= ii-1;ii>0;ii--,i--){
    k= 4*i;
    g[i]+= v4[k+0]*g[ii]+v4[k+1]*d2g[ii];
    d2g[i]+= v4[k+2]*g[ii]+v4[k+3]*d2g[ii];
  }
  return(0);
}

int EZinv2x2(w,rw)
     double w[],rw[];
{
  double d;
  d= 1./(w[0]*w[3]-w[1]*w[2]);
  rw[0]= w[3]*d;
  rw[1]= -w[1]*d;
  rw[2]= -w[2]*d;
  rw[3]= w[0]*d;
  return(0);
}

int spli(double *f,double *g,double *d2g,double *x,int *nx)
{
  int i,ii;
  double h,hh,eg;

  h= x[1]-x[0];
  eg= 0.;
  f[0]= ccr2*g[0];
  for(i= 0,ii= 1;i<(*nx);i++,ii++){
    h= x[ii]-x[i];
    hh= h*h;
    eg+= x[i]*ccr2*h*(g[i]+g[ii]-ccr12*hh*(d2g[i]+d2g[ii]));
    eg+= hh*((ccr6*g[i]+ccr3*g[ii])-(7.*ccr12*d2g[i]+2.*ccr3*d2g[ii])
	     *hh*ccr30);
    f[ii]= eg/(x[ii]*x[ii]);
  }
  return(0);
}
#endif

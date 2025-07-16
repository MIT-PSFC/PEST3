#include <string.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>

#ifndef mzl_MLIB
double EZcr2,EZcr3,EZcr4,EZcr6,EZcr12,EZcgp,EZc2gp,EZcgp_4,EZcr2gp,EZcgm0,EZcrgm0;
#ifdef AAASAS
int EZout(char *str)
#else
int EZout(char *str,...)
#endif
{
  va_list ap;
  va_start(ap,str);
  while(*str != '\0'){
    switch(*str){
    case 'i':
      printf("[%2d] ",va_arg(ap,int));
      break;
    case 'd':
      printf("%10.3e ",va_arg(ap,double));
      break;
    case 's':
      printf("%s",va_arg(ap,char*));
      break;
    case 'I':
      printf("%d ",va_arg(ap,int));
      break;
    case 'D':
      printf("%12.5e ",va_arg(ap,double));
      break;
    case 'x':
      printf("0x%8.8x ",va_arg(ap,long));
      break;
    case 'l':
      printf("0x%ld ",va_arg(ap,long));
      break;
    case 'c':
      printf("%c ",va_arg(ap,int));
      break;
    }
    str++;
  }
  putchar('\n');
  va_end(ap);
  return(0);
}

int IfPatternFound(FILE *lF,char *Str)
{
  char ch,*lc;

  ch	='\0';
  lc	=Str;
  while(*lc != '\0' && feof(lF) == 0){
    ch	=(char)fgetc(lF);
    if(*lc == ch){
      lc++;
    }
    else{
      lc	=Str;
    }
  }
  if(feof(lF)){
    return(0);
  }
  return(1);
}

int InitConstants()
{
  EZcr2= 0.5;
  EZcr3= 1./3.;
  EZcr4= 0.25;
  EZcr6= 1./6.;
  EZcr12= 1./12.;
  EZcgp_4= atan((double)1.);
  EZcgp= 4.*EZcgp_4;
  EZc2gp= 8.*EZcgp_4;
  EZcr2gp= 1./EZc2gp;
  EZcgm0= 0.4*EZcgp;
  EZcrgm0= 1./EZcgm0;
  return(0);
}

int EZCmsqrt(double*gr,double*gi,double ggr,double ggi)
{
  static double s;
  s= sqrt(ggr*ggr+ggi*ggi);
  *gr= ggr>=0.?sqrt(EZcr2*(s+ggr)):sqrt(EZcr2*ggi*ggi/(s-ggr));
  *gi= (*gr)!=0.?EZcr2*ggi/(*gr):sqrt(fabs(ggr));
  return(0);
}

#define TINY 1.0e-50;
static double *vv=NULL;
static int Nsys=0;

int ReInitLUdcmp(int n)
{
  if(n > Nsys || vv == NULL){
    if(vv != NULL){
      free(vv);
    }
    Nsys =n;
    vv =(double*)malloc(Nsys*sizeof(double));
    if(vv == NULL){
      printf("failure to allocate memory for vv in ReInitLUdcmp\n");
      exit(1);
    }
  }
  return(0);
}

int DeInitLUdcmp()
{
  if(vv != NULL){
    free(vv);
    vv	=NULL;
  }
  return(0);
}

int LUdcmp(a,n,indx,d)
     int n,*indx;
     double*a,*d;
{
  int i,imax,j,k;
  double big,dum,sum,temp;

  if(n > Nsys || vv == NULL){
    ReInitLUdcmp(n);
  }

  *d= 1.0;
  for(i=0; i < n; i++){
    big		=0.0;
    for(j=0; j < n; j++){
      if((temp=fabs(a[i*n+j])) > big){
	big	=temp;
      }
    }
    if(big == 0.0){
      printf("i=%d Singular matrix in routine LUdcmp\n",i);
      return(-1);
    }
    vv[i]	=1.0/big;
  }
  for(j=0; j < n; j++){
    for(i=0; i < j; i++){
      sum	=a[i*n+j];
      for(k=0; k < i; k++){
	sum	-=a[i*n+k]*a[k*n+j];
      }
      a[i*n+j]	=sum;
    }
    big	=0.0;
    for(i=j; i < n; i++){
      sum	=a[i*n+j];
      for(k=0; k < j; k++){
	sum	-=a[i*n+k]*a[k*n+j];
      }
      a[i*n+j]	=sum;
      if((dum=vv[i]*fabs(sum)) >= big){
	big	=dum;
	imax	=i;
      }
    }
    if(j != imax){
      for(k=0; k < n; k++){
	dum		=a[imax*n+k];
	a[imax*n+k]	=a[j*n+k];
	a[j*n+k]	=dum;
      }
      *d	=-(*d);
      vv[imax]	=vv[j];
    }
    indx[j]	=imax;
    if(a[j*n+j] == 0.0){
      a[j*n+j]	=TINY;
    }
    if(j != n-1){
      dum	=1.0/(a[j*n+j]);
      for(i=j+1; i < n; i++){
	a[i*n+j]	*=dum;
      }
    }
  }
  return(0);
}

int LUbksb(a,n,indx,b)
     int n,*indx;
     double*a,b[];
{
  int i,ii,ip,j;
  double sum;
  ii= 0;
  for(i= 0;i<n;i++){
    ip= indx[i];
    sum= b[ip];
    b[ip]= b[i];
    if(ii){
      for(j=ii-1;j<i;j++){
	sum-= a[i*n+j]*b[j];
      }
    }
    else{
      if(sum){
	ii= i+1;
      }
    }
    b[i]= sum;
  }
  for(i= n-1;i>=0;i--){
    sum= b[i];
    for(j= i+1;j<n;j++){
      sum-= a[i*n+j]*b[j];
    }
    b[i]= sum/a[i*n+i];
  }
  return(0);
}

int LUdcmpComp(ar,ai,n,indx)
     int n,*indx;
     double *ar,*ai;
{
  int i,imax,j,k;
  double big,dum,temp;
  double sr,si,tr,ti,Sr,Si;
  

  if(n > Nsys || vv == NULL){
    ReInitLUdcmp(n);
  }
  for(i=0; i < n; i++){
    big		=0.;
    for(j=0; j < n; j++){
      sr	=ar[i*n+j];
      si	=ai[i*n+j];
      temp	=sr*sr+si*si;
      if(temp > big){
	big	=temp;
      }
    }
    if(big == 0.){
      printf("i=%d Singular matrix in routine LUdcmpComp\n",i);
      return(-1);
    }
    vv[i]	=1./sqrt(big);
  }

  for(j=0; j < n; j++){
    for(i=0; i < j; i++){
      Sr	=ar[i*n+j];
      Si	=ai[i*n+j];
      for(k=0; k < i; k++){
	sr	=ar[i*n+k];
	si	=ai[i*n+k];
	tr	=ar[k*n+j];
	ti	=ai[k*n+j];
	Sr	-=sr*tr-si*ti;
	Si	-=sr*ti+si*tr;
      }
      ar[i*n+j]	=Sr;
      ai[i*n+j]	=Si;
    }

    big	=0.;
    for(i=j; i < n; i++){
      Sr	=ar[i*n+j];
      Si	=ai[i*n+j];
      for(k=0; k < j; k++){
	sr	=ar[i*n+k];
	si	=ai[i*n+k];
	tr	=ar[k*n+j];
	ti	=ai[k*n+j];
	Sr	-=sr*tr-si*ti;
	Si	-=sr*ti+si*tr;
      }
      ar[i*n+j]	=Sr;
      ai[i*n+j]	=Si;
      dum	=vv[i]*vv[i]*(Sr*Sr+Si*Si);
      if(dum >= big){
	big	=dum;
	imax	=i;
      }
    }
    if(j != imax){
      for(k=0; k < n; k++){
	sr		=ar[imax*n+k];
	si		=ai[imax*n+k];
	ar[imax*n+k]	=ar[j*n+k];
	ai[imax*n+k]	=ai[j*n+k];
	ar[j*n+k]	=sr;
	ai[j*n+k]	=si;
      }
      vv[imax]	=vv[j];
    }
    indx[j]	=imax;
    if(ar[j*n+j] == 0.0 && ai[j*n+j] == 0.0){
      ar[j*n+j]	=TINY;
    }
    if(j != n-1){
      sr	=ar[j*n+j];
      si	=ai[j*n+j];
      dum	=1.0/(sr*sr+si*si);
      sr	*=dum;
      si	*=-dum;
      for(i=j+1; i < n; i++){
	tr	=ar[i*n+j];
	ti	=ai[i*n+j];
	ar[i*n+j]	=sr*tr-si*ti;
	ai[i*n+j]	=sr*ti+si*tr;
      }
    }
  }
  return(0);
}

int LUbksbComp(ar,ai,n,indx,br,bi)
     int n,*indx;
     double *ar,*ai,*br,*bi;
{
  int i,ii,ip,j;
  double sr,si,t,Sr,Si;

  ii	=0;
  for(i=0; i < n; i++){
    ip	=indx[i];
    Sr	=br[ip];
    Si	=bi[ip];
    br[ip]=br[i];
    bi[ip]=bi[i];
    if(ii){
      for(j=ii-1; j < i; j++){
	sr	=ar[i*n+j];
	si	=ai[i*n+j];
	Sr	-=sr*br[j]-si*bi[j];
	Si	-=sr*bi[j]+si*br[j];
      }
    }
    else{
      if(Sr !=0. || Si != 0.){
	ii	=i+1;
      }
    }
    br[i]	=Sr;
    bi[i]	=Si;
  }
  for(i=n-1; i >= 0; i--){
    Sr	=br[i];
    Si	=bi[i];
    for(j=i+1; j < n; j++){
      sr	=ar[i*n+j];
      si	=ai[i*n+j];
      Sr	-=sr*br[j]-si*bi[j];
      Si	-=sr*bi[j]+si*br[j];
    }
    sr	=ar[i*n+i];
    si	=ai[i*n+i];
    t	=1./(sr*sr+si*si);
    sr	*=t;
    si	*=-t;
    br[i]=Sr*sr-Si*si;
    bi[i]=Sr*si+Si*sr;
  }
  return(0);
}

int luBKSB(double *a,int n,int *indx,double *bb,int K)
{
  static int i,k,ii,in,ip,j;
  static double sum,*b;

  for(k=0; k < K; k++){
    b	=bb+n*k;
    ii	=0;
    for(i=0; i < n; i++){
      ip	=indx[i];
      sum	=b[ip];
      b[ip]	=b[i];
      if(ii){
	j	=ii-1;
	in	=i*n+j;
	while(j < i){
	  sum	-=a[in]*b[j];
	  j++;
	  in++;
	}
      }
      else{
	if(sum){
	  ii	=i+1;
	}
      }
      b[i]	=sum;
    }
    i	=n;
    in	=n*n;
    while(i > 0){
      j		=i;
      i--;
      in	-=n;
      sum	=b[i];
      while(j < n){
	sum	-=a[in+j]*b[j];
	j++;
      }
      b[i]	=sum/a[in+i];
    }
  }
  return(0);
}

int EZMatinvert(a,ra,n)
     double*a,*ra;
     int n;
{
  static int i,ii,ip,j,k;
  static double s,sum,*v;
  static int*index= NULL,Nsys= 0;

  if(n>Nsys){
    if(index!=NULL){
      free(index);
      free(v);
    }
    Nsys= n;
    index= (int*)malloc(Nsys*sizeof(int));
    v= (double*)malloc(Nsys*sizeof(double));
    if(index==NULL||v==NULL){
      printf("failure to allocate memory for index in EZMatinvert()\n");
      exit(1);
    }
  }
  if(LUdcmp(a,n,index,&sum)==0)
    return(0);
  for(k= 0;k<n;k++){
    for(i= 0;i<n;i++)
      v[i]= 0.;
    v[k]= 1.;
    ii= 0;
    for(i= 0;i<n;i++){
      ip= index[i];
      sum= v[ip];
      v[ip]= v[i];
      if(ii)
	for(j= ii-1;j<i;j++)
	  sum-= a[i*n+j]*v[j];
      else
	if(sum)
	  ii= i+1;
      v[i]= sum;
    }
    for(i= n-1;i>=0;i--){
      sum= v[i];
      for(j= i+1;j<n;j++)
	sum-= a[i*n+j]*v[j];
      v[i]= sum/a[i*n+i];
    }
    for(i= 0;i<n;i++)
      ra[i*n+k]= v[i];
  }
  return(1);
}

#ifdef H
int CholDc(double *a, int n)
{
  int i,ii,j,jj,k,kk,kj,m;
  double s;

  m	=0;
  ii	=0;
  for(i=0; i < n; i++){
    jj	=ii;
    for(j=i; j < n; j++){
      s		=a[ii+j];
      k		=i;
      kk	=ii+i;
      kj	=jj+i;
      while(k){
	k--;
	kk--;
	kj--;
	s	-=a[kk]*a[kj];
      }
      if(i == j){
	if(s <= 0.){
#ifdef H
	  printf("s[%2d]=%10.3e <= 0. in CholDc()%c\n",i,s,7);
#endif
	  m	=-1;
	  s	=-s;
	}
	a[ii+i]	=sqrt(s);
      }
      else{
	a[jj+i]=s/a[ii+i];
      }
      jj	+=n;
    }
    ii	+=n;
  }
  return(m);
}

int CholBksb(double *a, int n, double *b)
{
  int i,ii,k,kk;
  double s;
  ii	=0;
  for(i=0; i < n; i++){
    s	=b[i];
    k	=i;
    kk	=ii+i;
    while(k > 0){
      k--;
      kk--;
      s	-=a[kk]*b[k];
    }
    b[i]	=s/a[ii+i];
    ii	+=n;
  }
  while(i > 0){
    i--;
    kk	=ii+i;
    ii	-=n;
    s	=b[i];
    k	=i+1;
    while(k < n){
      s	-=a[kk]*b[k];
      k++;
      kk	+=n;
    }
    b[i]	=s/a[ii+i];
  }
  return(0);
}

int CholBksb(double *K, int n, double *b)
{
  static Fl=1;

  int i,ii,j,jj;
  double s;
  ii	=0;
  for(i=0; i < n; i++){
    s	=b[i];
    j	=i;
    jj	=ii+i;
    while(j > 0){
      j--;
      jj--;
      s	-=K[jj]*b[j];
    }
    b[i]	=s/K[ii+i];
    ii	+=n;
  }

  while(i > 0){
    i--;
    jj	=ii+i;
    ii	-=n;
    s	=b[i];
    j	=i+1;
    while(j < n){
      s	-=K[jj]*b[j];
      j++;
      jj	+=n;
    }
    b[i]	=s/K[ii+i];
  }
  return(0);
}
#endif

int CholDc(double *K, int n)
{     
  int i,j,k,m;
  double s,d,*p,*q;
     
  m	=0;
  p	=K;
  for(i=0; i < n; i++){
    j	=i;
    k	=i;
    s	=p[i];
    while(k){
      k--;
      s	-=p[k]*p[k];
    }
    if(s <= 0.){
#ifdef H
      printf("s[%2d]=%10.3e <= 0. in CholDc()%c\n",i,s,7);
#endif
      m	=-1;
      s	=-s;
      if(s == 0.){
	s	=1e-20;
      }
    }
    p[i]=sqrt(s);
    d	=1./p[i];
    q	=p;
    for(j=i+1; j < n; j++){
      q		+=n;
      k		=i;
      s		=p[j];
      while(k){
	k--;
	s	-=p[k]*q[k];
      }
      q[i]	=s*d;
    }
    p	+=n;
  }
  return(m);
}

int CholBksb(double *K, int n, double *b)
{
  static int Fl=1;

  int i,ii,j,jj;
  double s;
  ii	=0;
  for(i=0; i < n; i++){
    s	=b[i];
    j	=i;
    jj	=ii+i;
    while(j > 0){
      j--;
      jj--;
      s	-=K[jj]*b[j];
    }
    b[i]	=s/K[ii+i];
    ii	+=n;
  }

  while(i > 0){
    i--;
    jj	=ii+i;
    ii	-=n;
    s	=b[i];
    j	=i+1;
    while(j < n){
      s	-=K[jj]*b[j];
      j++;
      jj	+=n;
    }
    b[i]	=s/K[ii+i];
  }
  return(0);
}

int CholeskyMatInv(double *a, int n, double *p)
{
  int i,ii,j,jj,k,kk,kj;
  double s;
  ii	=0;
  for(i=0; i < n; i++){
    jj	=ii;
    for(j=i; j < n; j++){
      s		=a[ii+j];
      k		=i;
      kk	=ii+i;
      kj	=jj+i;
      while(k > 0){
	k--;
	kk--;
	kj--;
	s	-=a[kk]*a[kj];
      }
      if(i == j){
	p[i]	=sqrt(s);
      }
      else{
	a[jj+i]=s/p[i];
      }
      jj	+=n;
    }
    ii	+=n;
  }
  ii	=0;
  for(i=0; i < n; i++){
    a[ii+i]	=1./p[i];
    jj	=ii+n;
    for(j=i+1; j < n; j++){
      s		=0.;
      kk	=ii+i;
      kj	=jj+i;
      for(k=i; k < j; k++){
	s	-=a[kj]*a[kk];
	kk	+=n;
	kj++;
      }
      a[jj+i]	=s/p[i];
      jj	+=n;
    }
    ii	+=n;
  }
  return(0);
}

int CholDeComp(double *Kr, double *Ki, int n)
{     
  int i,j,k,m;
  double sr,si,*pr,*qr,*pi,*qi;
     
  m	=0;
  pr	=Kr;
  pi	=Ki;
  for(i=0; i < n; i++){
    j	=i;
    k	=i;
    sr	=pr[i];
    while(k){
      k--;
      sr-=pr[k]*pr[k]+pi[k]*pi[k];
    }
    if(sr <= 0.){
      printf("sr[%2d]=%10.3e <= 0. in CholDc()%c\n",i,sr,7);
#ifdef H
#endif
      m	=-1;
      sr=-sr;
      if(sr == 0.){
	sr	=1e-20;
      }
    }
    pr[i]=1./sqrt(sr);
    qr	=pr;
    qi	=pi;
    for(j=i+1; j < n; j++){
      qr	+=n;
      qi	+=n;
      k		=i;
      sr	=pr[j];
      si	=pi[j];
      while(k){
	k--;
	sr	-=pr[k]*qr[k]+pi[k]*qi[k];
	si	-=pi[k]*qr[k]-pr[k]*qi[k];
      }
      qr[i]	=sr*pr[i];
      qi[i]	=-si*pr[i];
    }
    pr	+=n;
    pi	+=n;
  }
  return(m);
}

int CholBksbComp(double *Kr, double *Ki, int n, double *br, double *bi)
{
  int i,ii,j,jj;
  double sr,si;

  ii	=0;
  for(i=0; i < n; i++){
    sr	=br[i];
    si	=bi[i];
    j	=i;
    jj	=ii+i;
    while(j > 0){
      j--;
      jj--;
      sr	-=Kr[jj]*br[j]-Ki[jj]*bi[j];
      si	-=Kr[jj]*bi[j]+Ki[jj]*br[j];
    }
    br[i]	=sr*Kr[ii+i];
    bi[i]	=si*Kr[ii+i];
    ii	+=n;
  }

  while(i > 0){
    i--;
    jj	=ii+i;
    ii	-=n;
    sr	=br[i];
    si	=bi[i];
    j	=i+1;
    while(j < n){
      sr	-=Kr[jj]*br[j]+Ki[jj]*bi[j];
      si	-=Kr[jj]*bi[j]-Ki[jj]*br[j];
      j++;
      jj	+=n;
    }
    br[i]	=sr*Kr[ii+i];
    bi[i]	=si*Kr[ii+i];
  }
  return(0);
}

void EZiround(ch,k,ne,i5)
     char*ch;
     int*k,*ne,*i5;
{
  int i;
  i= *k-1;
  if(ch[*k]>=53+*i5){
    ch[i]++;
    *i5= 0;
    if(ch[i]==53)
      *i5= 1;
  }
  while(ch[i]==58){
    if(i>0){
      ch[i]= 48;
      i--;
      ch[i]++;
      *i5= 0;
      if(ch[i]==53)
	*i5= 1;
      (*k)--;
      (*ne)++;
    }
    else{
      ch[i]= 49;(*ne)++;
    }
  }
  while(ch[(*k)-1]==48){
    if((*k)>1){
      (*k)--;
      (*ne)++;
    }
  }
}

int EZnumstr(x,str,l)
     char*str;
     double x;
     int l;
{
  static int n= 12;
  int i,ne,is,ll,k,i5;
  char ch[30];
  
  for(i= 0;i<l;i++)str[i]= ' ';
  str[l]= '\0';
  if(x==0.){
    str[l-1]= '0';
    return(0);
  }
  sprintf(ch,"%+.*e",n,x);
  i= sscanf((ch+n+4),"%d",&ne);
  if(ch[0]=='-')is= 1;else is= 0;
  ll= l-is;
  if(ll<=0){
    str[l-1]= '*';
    return(1);
  }
  k= n+2;
  while(ch[k]=='0')k--;
  ch[0]= ch[1];
  for(i= 1;i<k;ch[i]= ch[i+2],i++);
  k--;
  ch[k]= '\0';
  ne+= (1-k);
  i5= 0;
  if(k>ll){
    ne= ne+k-ll;
    k= ll;
    EZiround(ch,&k,&ne,&i5);
  }
  while(k){
    if(ne>=0){
      if(ne+k<=ll){
	for(i= l-1;ne>0;ne--,str[i--]= '0');
	for(k--;k>=0;str[i--]= ch[k--]);
	if(is)str[i]= '-';
	return(0);
      }
      else{
	i= 2;
	if(ne>9)i++;
	if(ne>99)i++;
	if(k+i>ll){
	  k--;
	  if(k==0){
	    str[l-1]= '*';
	    return(1);
	  }
	  ne++;
	  EZiround(ch,&k,&ne,&i5);
	}
	else{
	  i= l-i;
	  sprintf(str+i,"e%d",ne);
	  for(--k;k>=0;str[--i]= ch[k--]);
	  if(is)str[--i]= '-';
	  return(0);
	}
      }
    }
    else{
      if(ll>-ne){
	if(k-ll){
	  i= k+ne;
	  if(i>=0){
	    i= l;
	    while(ne){
	      k--;
	      i--;
	      str[i]= ch[k];
	      ne++;
	    }
	    i--;
	    str[i]= '.';
	    while(k){
	      k--;
	      i--;
	      str[i]= ch[k];
	    }
	  }
	  else{
	    i= l;
	    while(k){
	      k--;
	      i--;
	      str[i]= ch[k];
	      ne++;
	    }
	    while(ne){
	      i--;
	      str[i]= '0';
	      ne++;
	    }
	    i--;
	    str[i]= '.';
	  }
	  if(is){
	    i--;
	    str[i]= '-';
	  }
	  return(0);
	}
	else{
	  k--;
	  ne++;
	  EZiround(ch,&k,&ne,&i5);
	}
      }
      else{
	i= 3;
	if(ne<-9)i= 4;
	if(ne<-99)i= 5;
	if(ll>=k+i){
	  i= l-i;
	  sprintf(str+i,"e%d",ne);
	  while(k){
	    i--;
	    k--;
	    str[i]= ch[k];
	  }
	  if(is){
	    i--;
	    str[i]= '-';
	  }
	  return(0);
	}
	k--;
	if(k==0){
	  str[l-1]= '*';
	  return(1);
	}
	ne++;
	EZiround(ch,&k,&ne,&i5);
      }
    }

  }
  return(0);
}

int u_choice(double*xmin,int*imin,double*xmax,int*imax,double*unitx)
{
  double s,u;
  int i,j,k;
  s	=(*xmax)-(*xmin);
  u	=fabs(*xmax)+fabs(*xmin);
  if(u == 0. || s < 1e-12*u){
    if(*xmax > 0.){
      *xmin	=0.;
    }
    else{
      if(*xmax == 0.){
	*xmax	=1.;
      }
      else{
	*xmax	=0.;
      }
    }
  }
#ifdef H
  if(*xmax == *xmin){
    if(*xmax > 0.){
      *xmin	=0.;
    }
    else{
      if(*xmax == 0.){
	*xmax	=1.;
      }
      else{
	*xmax	=0.;
      }
    }
  }
#endif
  u	=*xmax-*xmin;
  s	=log10(u);
  i	=s;
  if(s < 0.)
    i--;
  u	=pow(10.,(double)i);
  s	=(*xmax-*xmin)/u;
  if(s > 1.2){
    if(s > 2.5){
      if(s > 5.){
	u	*= 2.;
      }
    }
    else{
      u	*= 0.5;
    }
  }
  else{
    u	*= 0.2;
  }
  *unitx	=u;
  s=(*xmin)/u+1.e-6;
  if((*xmin)>0){
    j= (int)s;
    i= j+1;
  }
  else{
    i=(int)s;
    j= i-1;
  }
  s-= j;
  if(s<0.25){
    *xmin= j*u;
    *imin= j;
  }
  else{
    *imin= i;
    if(s < 0.5)
      *xmin= (j+0.25)*u;
    else if(s<0.75)
      *xmin= (j+0.5)*u;
    else if(s<1)
      *xmin= (j+0.75)*u;
  }
  s=(*xmax)/u-1.e-6;
  if((*xmax) > 0){
    j= (int)s;
    i= j+1;
  }
  else{
    i= (int)s;
    j= i-1;
  }
  s= i-s;
  if(s < 0.25){
    *xmax=i*u;
    *imax=i;
  }
  else{
    *imax= j;
    if(s<0.5)
      *xmax=(j+0.75)*u;
    else if(s < 0.75)
      *xmax=(j+0.5)*u;
    else if(s < 1.)
      *xmax= (j+0.25)*u;
  }
  return(0);
}

#ifdef H
void main0()
{
  char s[80],s1[80];
  double xmn,xmx,u;
  int imn,imx;
  int l= 4;
  while(1){
    printf("xmn, xmx =?\n");
    gets(s);
    if(s[0]=='c')break;
    imn= sscanf(s,"%lf %lf",&xmn,&xmx);
    u_choice(&xmn,&imn,&xmx,&imx,&u);
    printf("xmn=%8e %8e; %8e xmx=%8e\n",xmn,u*imn,xmx,u*imx);
    EZnumstr(xmn,s1,l);
    printf("%s, x=%8g \n",s1,xmn);
  }
}
#endif

int cpy(Targ,Source)
     char*Targ,*Source;
{
  static int i,j;
  j= 0;
  while(Source[j]==' ')
    j++;
  i= 0;
  while(Source[j]!=' '&&Source[j]!='\0'){
    Targ[i]= Source[j];
    j++;
    i++;
  }
  return(i);
}

int strcpycore(p,pp)
     char*p,*pp;
{
  while(isspace(*pp)&&*pp!='\0')
    pp++;
  *p= *pp;
  while(*pp!='\0'&&!isspace(*pp)){
    *p= *pp;
    p++;
    pp++;
  }
  *p= '\0';
  return(0);
}

int strncpyleft(p,pp,Ln)
     char*p,*pp;
     int Ln;
{
  int i;
  i= 0;
  while(*pp!='\0'){
    *p= *pp;
    p++;
    pp++;
    i++;
  }
  while(i<Ln){
    *p= ' ';
    p++;
    i++;
  }
  *p= '\0';
  return(0);
}

int EZzero(int*iter,double *x,double *f,double *x0,double f0)
{
  static double s0,s1,dx;
  (*iter)++;
  if((*iter) == 1){
    (*x0)	=x[0];
    return(0);
  }
  if((*iter) == 2){
    (*x0)	=x[1];
    f[0]	=f0;
    dx		=fabs(x[1]-x[0]);
    return(0);
  }
  x[1]	=(*x0);
  f[1]	=f0;
  (*x0)	=x[0]-(x[1]-x[0])*f[0]/(f[1]-f[0]);
  s0	=f[0] < 0. ? -f[0] : f[0];
  s1	=f[1] < 0. ? -f[1] : f[1];
  if(s0 > s1){
    x[0]	=x[1];
    f[0]	=f[1];
  }
  s0	=fabs((*x0)-x[0])/dx;
  if(s0 < 1e-10){
    (*iter)	=-1;
  }
  s0	=fabs((*x0)-x[0])/dx;
  if(s0 > 100.){
    *iter	=-2;
  }
  return(0);
}

int EZZero(int*iter,double *x,double *f,double *x0,double f0)
{
  static double s0,s1,dx;

  (*iter)++;
  if(*iter == 1){
    *x0		=x[0];
    return(1);
  }
  if(*iter == 2){
    *x0		=x[1];
    f[0]	=f0;
    dx		=fabs(x[1]-x[0]);
    return(2);
  }

  x[1]	=*x0;
  f[1]	=f0;
  *x0	=x[0]-(x[1]-x[0])*f[0]/(f[1]-f[0]);
  s0	=f[0] < 0. ? -f[0] : f[0];
  s1	=f[1] < 0. ? -f[1] : f[1];
  if(s0 > s1){
    x[0]	=x[1];
    f[0]	=f[1];
  }
  s0	=fabs((*x0)-x[0])/dx;

  if(s0 < 1e-10){
    return(0);
  }
  s0	=fabs((*x0)-x[0])/dx;
  if(s0 > 100.){
    return(-1);
  }
  return(*iter);
}

int EZzero1(iter,inst,x,f,vx,vf,racc)
     int*iter,inst;
     double*x,*f,*vx,vf,racc;
{
  static double s,dx,xm,dxf,d2dx,b0,b1,xlim,s0,df,d2f;
  static int fl;

  (*iter)++;
  if((*iter)==1){
    (*vx)= x[0];
    fl= 0;
    return(0);
  }
  if((*iter)==2){
    if(inst){
      (*iter)--;
      fl= 1;
      xlim= *vx;
      x[0]+= x[0]-x[1];
      x[1]= EZcr2*(x[0]+xlim);
      (*vx)= x[0];
      return(0);
    }
    if(vf<0.){
      (*iter)--;
      x[0]+= x[0]-x[1];
      x[1]= *vx;
      (*vx)= x[0];
      return(0);
    }
    f[0]= vf;
    (*vx)= x[1];
    dx= fabs(x[1]-x[0]);
    if(dx<fabs(x[1]))dx= fabs(x[1]);
    if(dx<fabs(x[0]))dx= fabs(x[0]);
    return(0);
  }

  if((*iter)==3){
    if(inst){
      (*iter)--;
      fl= 1;
      xlim= *vx;
      x[1]= EZcr2*(x[0]+x[1]);
      (*vx)= x[1];
      return(0);
    }
    if(vf>0.){
      (*iter)--;
      if(fl){
	x[1]= EZcr2*(x[1]+xlim);
      }
      else
	x[1]= x[1]-(x[0]-x[1]);
      x[0]= *vx;
      f[0]= vf;
      (*vx)= x[1];
      return(0);
    }
    f[1]= vf;
    (*vx)= x[0]-(x[1]-x[0])*f[0]/(f[1]-f[0]);
    return(0);
  }
  b0= (f[0]-vf)/(x[0]-(*vx));
  b1= (f[1]-vf)/(x[1]-(*vx));
  d2f= (b1-b0)/(x[1]-x[0]);
  df= EZcr2*(b0-d2f*(x[0]-(*vx)));
  s= sqrt(df*df-vf*d2f);
  b0= -vf/(df+s);
  b1= -vf/(df-s);
  s= b0+(*vx);
  s= x[1]<=s&&s<=x[0]?b0:b1;
  if(vf>0.){
    x[0]= *vx;
    f[0]= vf;
  }
  else{
    x[1]= *vx;
    f[1]= vf;
  }
  if(fabs(x[0]-x[1])>4.*fabs(s))
    s*= 2.;
  *vx+= s;
  if(fabs(x[1]-x[0])/dx<racc){
#ifdef H
    printf("end=%2d xdx=%11.4e%11.4e%11.4e %11.4e%11.4e%11.4e\n",*iter,
	   x[0],x[1],dx,f[0],f[1],*vx);
#endif
    (*iter)= -1;
  }
  if(fabs((*vx))/dx>100.)
    (*iter)= -2;
  if(*iter==-1)
    return(0);
  return(0);
}

int zero2D(int*iter,double x[],double y[],double *x0,double *y0,double f0
	   ,double g0,double eps) 
{
  static double a11,a12,a21,a22,b1,b2,dxf,dxg,s,s0,s1,dx,dy;
  static double f[3],g[3];
  (*iter)++;
  if((*iter)==1){
    (*x0)= x[0];
    (*y0)= y[0];
    return(0);
  }
  if((*iter)==2){
    (*x0)= x[1];
    (*y0)= y[1];
    f[0]= f0;
    g[0]= g0;
    dx= fabs(x[1]-x[0]);
    dy= fabs(y[1]-y[0]);
    return(0);
  }
  if((*iter)==3){
    (*x0)= x[2];
    (*y0)= y[2];
    f[1]= f0;
    g[1]= g0;
    if(dx<fabs(x[2]-x[0]))dx= fabs(x[2]-x[0]);
    if(dy<fabs(y[2]-y[0]))dy= fabs(y[2]-y[0]);
    return(0);
  }
  
  x[2]= (*x0);
  y[2]= (*y0);
  f[2]= f0;
  g[2]= g0;
  a11= f[1]-f[0];
  a12= g[1]-g[0];
  a21= f[2]-f[0];
  a22= g[2]-g[0];
  s= 1./(a22*a11-a12*a21);
  b1= x[1]-x[0];
  b2= x[2]-x[0];
  dxf= (b1*a22-b2*a12)*s;
  dxg= (b2*a11-b1*a21)*s;
  *x0= x[0]-dxf*f[0]-dxg*g[0];
  b1= y[1]-y[0];
  b2= y[2]-y[0];
  dxf= (b1*a22-b2*a12)*s;
  dxg= (b2*a11-b1*a21)*s;
  *y0= y[0]-dxf*f[0]-dxg*g[0];
  
  b2= (x[2]-(*x0))/dx;
  b2*= b2;
  s= (y[2]-(*y0))/dy;
  s*= s;
  b2+= s;
  b1= (x[1]-(*x0))/dx;
  b1*= b1;
  s= (y[1]-(*y0))/dy;
  s*= s;
  b1+= s;
  if(b2<b1){
    s= x[1];
    x[1]= x[2];
    x[2]= s;
    s= y[1];
    y[1]= y[2];
    y[2]= s;
    s= f[1];
    f[1]= f[2];
    f[2]= s;
    s= g[1];
    g[1]= g[2];
    g[2]= s;
  }
  b2= (x[2]-(*x0))/dx;
  b2*= b2;
  s= (y[2]-(*y0))/dy;
  s*= s;
  b2+= s;
  b1= (x[0]-(*x0))/dx;
  b1*= b1;
  s= (y[0]-(*y0))/dy;
  s*= s;
  b1+= s;
  if(b2<b1){
    s= x[0];
    x[0]= x[2];
    x[2]= s;
    s= y[0];
    y[0]= y[2];
    y[2]= s;
    s= f[0];
    f[0]= f[2];
    f[2]= s;
    s= g[0];
    g[0]= g[2];
    g[2]= s;
  }
  b2= (x[2]-(*x0))/dx;
  b2*= b2;
  s= (y[2]-(*y0))/dy;
  s*= s;
  b2+= s;
  if(b2<eps*eps)(*iter)= -1;
  if(b2>100.)*iter= -2;
  return(0);
}

int guess2D(x,y,f,g,xv,yv,dxf,dxg,dyf,dyg)
     double*x,*y,*f,*g,*xv,*yv,*dxf,*dxg,*dyf,*dyg;
{
  static double a11,a12,a21,a22,b1,b2,s;
  a11= f[1]-f[0];
  a12= g[1]-g[0];
  a21= f[2]-f[0];
  a22= g[2]-g[0];
  s= a22*a11-a12*a21;
  s= 1./s;
  b1= x[1]-x[0];
  b2= x[2]-x[0];
  *dxf= (b1*a22-b2*a12)*s;
  *dxg= (b2*a11-b1*a21)*s;
  b1= y[1]-y[0];
  b2= y[2]-y[0];
  *dyf= (b1*a22-b2*a12)*s;
  *dyg= (b2*a11-b1*a21)*s;
  *xv= x[0]-(*dxf)*f[0]-(*dxg)*g[0];
  *yv= y[0]-(*dyf)*f[0]-(*dyg)*g[0];
  return(0);
}

int zero2D1(iter,inst,x,y,x0,y0,f0,g0,eps)
     int*iter,inst;
     double x[],y[],*x0,*y0,f0,g0,eps;
{
  static double a11,a12,a21,a22,b1,b2,dxf,dxg,dyf,dyg,s,s0,s1,dx,dy;
  static double f[3],g[3],d[3],xc,yc,df,dg,xlim,ylim;
  static int i,i0,imax,fl;
  (*iter)++;
  if((*iter)==1){
    (*x0)= x[0];
    (*y0)= y[0];
    fl= 0;
    return(0);
  }
  if((*iter)==2){
    if(inst){
      (*iter)--;
      xlim= *x0;
      ylim= *y0;
      fl= 1;
      x[0]*= 2.;
      *x0= x[0];
      return(0);
    }
    f[0]= f0;
    g[0]= g0;
    if(fl){
      x[1]= EZcr2*(xlim+x[0]);
      x[2]= x[1];
    }
    (*x0)= x[1];
    (*y0)= y[1];
    dx= fabs(x[1]-x[0]);
    dy= fabs(y[1]-y[0]);
    return(0);
  }
  if((*iter)==3){
    if(inst){
      (*iter)--;
      xlim= *x0;
      ylim= *y0;
      fl= 1;
      x[1]= 0.5*(x[0]+x[1]);
      y[1]= 0.5*(y[0]+y[1]);
      x[2]= x[1];
      *x0= x[1];
      *y0= y[1];
      return(0);
    }
    f[1]= f0;
    g[1]= g0;
    (*x0)= x[2];
    (*y0)= y[2];
    if(dx<fabs(x[2]-x[0]))dx= fabs(x[2]-x[0]);
    if(dy<fabs(y[2]-y[0]))dy= fabs(y[2]-y[0]);
    return(0);
  }
  if((*iter)==4){
    if(inst){
      xlim= *x0;
      ylim= *y0;
      fl= 1;
      x[2]= 0.5*(x[0]+x[2]);
      y[2]= 0.5*(y[0]+y[2]);
      *x0= x[2];
      *y0= y[2];
      (*iter)--;
      return(0);
    }
    x[2]= (*x0);
    y[2]= (*y0);
    f[2]= f0;
    g[2]= g0;
    guess2D(x,y,f,g,x0,y0,&dxf,&dxg,&dyf,&dyg);
    return(0);
  }
  if(inst){
    xlim= *x0;
    ylim= *y0;
    fl= 1;
    *x0= EZcr2*(x[0]+(*x0));
    *y0= EZcr2*(y[0]+(*y0));
    return(0);
  }
  for(i= 0;i<3;i++){
    b1= x[i]-(*x0);
    b2= y[i]-(*y0);
    if(i){
      if(b1*b1+b2*b2<s)
	i0= i;
      if(f[i]*f[i]+g[i]*g[i]>s0)
	imax= i;
    }
    else{
      s= b1*b1+b2*b2;
      s0= f[i]*f[i]+g[i]*g[i];
      i0= i;
      imax= i;
    }
  }
  x[imax]= *x0;
  y[imax]= *y0;
  f[imax]= f0;
  g[imax]= g0;
  guess2D(x,y,f,g,x0,y0,&dxf,&dxg,&dyf,&dyg);
  b2= (fabs(x[0]-(*x0))+fabs(x[1]-(*x0))+fabs(x[2]-(*x0)))/dx;
  b2+= (fabs(y[0]-(*y0))+fabs(y[1]-(*y0))+fabs(y[2]-(*y0)))/dy;
  if(b2<eps)(*iter)= -1;
  if(b2>100.)*iter= -2;
  return(0);
}

int zero2D2(iter,inst,x,y,x0,y0,f0,g0,eps,dx,dy)
     int*iter,inst;
     double x[],y[],*x0,*y0,f0,g0,eps,dx,dy;
{
  static double b1,b2,s,dxf,dxg,dyf,dyg;
  static double f[4],g[4],ff[4],gg[4],xx[4],yy[4],sx[2],sy[2],dist;
  static int ifg,kf[4]= {0,0,2,2},kf1[4]= {1,0,3,2},
  kg[4]= {3,1,1,3},kg1[4]= {3,2,1,0},i,i1,j,iworst,iworst1;

  printf("Set of points Iter=%2d\n",*iter);
  printf("0 %11.4e%11.4e %11.4e%11.4e\n",f[0],g[0],x[0],y[0]);
  printf("1 %11.4e%11.4e %11.4e%11.4e\n",f[1],g[1],x[1],y[1]);
  printf("2 %11.4e%11.4e %11.4e%11.4e\n",f[2],g[2],x[2],y[2]);
  if(inst){
    *x0= EZcr2*((*x0)+x[0]);
    *y0= EZcr2*((*y0)+y[0]);
    (*iter)--;
    return(0);
  }
  b1= 0;
  for(i= 0;i<4;i++){
    b2= f[i]*f[i]+g[i]*g[i];
    if(b2>b1){
      b1= b2;
      iworst= i;
    }
  }
  b1= 0;
  for(i= 0,j= 0;i<4;i++){
    if(i-iworst){
      xx[j]= x[i];
      yy[j]= y[i];
      ff[j]= f[i];
      gg[j]= g[i];
      j++;
      b2= f[i]*f[i]+g[i]*g[i];
      if(b2>b1){
	b1= b2;
	iworst1= i;
      }
    }
  }
  if(*iter==5){
    guess2D(xx,yy,ff,gg,x0,y0,&dxf,&dxg,&dyf,&dyg);
    printf("new %11.4e%11.4e %1d\n",*x0,*y0,i1);
    return(0);
  }
  b1= 0;
  for(i= 0;i<3;i++){
    b2= ff[i]*ff[i]+gg[i]*gg[i];
    if(b2>b1){
      b1= b2;
      j= i;
    }
  }
  xx[j]= *x0;
  yy[j]= *y0;
  ff[j]= f0;
  gg[j]= g0;
  if(f0>0.){
    ifg= g0>=0?0:3;
  }
  else{
    if(f0==0.){
      ifg= g0>=0?0:3;
    }
    else{
      ifg= g0>=0?1:2;
    }
  }
  i= ifg;
  b2= (f[i]-f0)*(f[i]-f0)+(g[i]-g0)*(g[i]-g0);
  x[ifg]= *x0;
  y[ifg]= *y0;
  f[ifg]= f0;
  g[ifg]= g0;
  guess2D(xx,yy,ff,gg,x0,y0,&dxf,&dxg,&dyf,&dyg);
  i= iworst;
  b1= (f[i]-f0)*(f[i]-f0)+(g[i]-g0)*(g[i]-g0);
  if(ifg-iworst&&b1>25.*b2){
    b1= EZcr4*f[iworst];
    b2= EZcr4*g[iworst];
    *x0= (*x0)+dxf*b1+dxg*b2;
    *y0= (*y0)+dyf*b1+dyg*b2;
    printf("attempt %1d\n",iworst);
  }
  printf("new %11.4e%11.4e %1d %1d\n",*x0,*y0,i1,j);
  b1= (x[0]-x[2])/dx;
  b2= (y[0]-y[2])/dy;
  dist= b1*b1+b2*b2;
  printf("ifg=%1d disc1=%11.4e ",ifg,dist);
  if(dist<eps*eps)(*iter)= -1;
  b1= (x[1]-x[3])/dx;
  b2= (y[1]-y[3])/dy;
  dist= b1*b1+b2*b2;
  printf("disc2=%11.4e\n",dist);
  if(dist<eps*eps)(*iter)= -1;
  if(dist>1e+4)(*iter)= -2;
  return(0);
}

int zero2D22(iter,inst,x,y,x0,y0,f0,g0,eps,dx,dy)
     int*iter,inst;
     double x[],y[],*x0,*y0,f0,g0,eps,dx,dy;
{
  static double b1,b2,s,dxf,dxg,dyf,dyg;
  static double f[4],g[4],ff[4],gg[4],xx[4],yy[4],sx[2],sy[2],dist;
  static int ifg,kf[4]= {0,0,2,2},kf1[4]= {1,0,3,2},
  kg[4]= {3,1,1,3},kg1[4]= {3,2,1,0},i,i1,j,iworst,iworst1;
  typedef struct Guess{
    int fl;
    double x,y,f,g;
  }Guess;
  static Guess v[5];
  (*iter)++;
  if(*iter==1){
    *x0= x[0];
    *y0= y[0];
    v[0].fl= 1;
    v[1].fl= 1;
    v[2].fl= 1;
    v[3].fl= 1;
    v[4].fl= 1;
    return(0);
  }

  printf("itr=%2d\n",*iter);
  printf("fgxyIn =%11.4e%11.4e%11.4e%11.4e\n",f0,g0,*x0,*y0);
  if(*iter<6){
    printf("flags=%2d %2d %2d %2d %2d\n",
	   v[0].fl,v[1].fl,v[2].fl,v[3].fl,v[4].fl);
  }

  if(*iter==2){
    if(inst){
      sx[0]= x[0];
      sy[0]= y[0];
      v[4].fl= 0;
      x[0]= 2.*x[0];
      (*iter)--;
      *x0= x[0];
      *y0= y[0];
    }
    else{
      if(f0>=0.)
	ifg= g0>=0?0:3;
      else
	ifg= g0>=0?1:2;
      if(v[ifg].fl){
	v[ifg].fl= 0;
	v[ifg].f= f0;
	v[ifg].g= g0;
	v[ifg].x= *x0;
	v[ifg].y= *y0;
      }
      else{
	if(f0*f0+g0*g0<v[ifg].f*v[ifg].f+v[ifg].g*v[ifg].g){
	  v[ifg].fl= 0;
	  v[ifg].f= f0;
	  v[ifg].g= g0;
	  v[ifg].x= *x0;
	  v[ifg].y= *y0;
	}
      }
      if(ifg==1||ifg==2){
	x[0]= 2.*x[0];
	(*iter)--;
	*x0= x[0];
	*y0= y[0];
      }
      else{
	if(v[1].fl&&v[2].fl){
	  *iter= 2;
	  if(v[4].fl){
	    x[2]= x[2];
	  }
	  else{
	    x[2]= EZcr2*(sx[0]+x[0]);
	  }
	  *x0= x[2];
	  *y0= y[2];
	}
	else{
	  if(v[0].fl||v[3].fl){
	    *iter= 3;
	    x[1]= EZcr2*(x[0]+x[2]);
	    *x0= x[1];
	    *y0= y[1];
	  }
	  else{
	    if(v[1].fl||v[2].fl){
	      *iter= 4;
	      x[3]= EZcr2*(x[0]+x[2]);
	      *x0= x[3];
	      *y0= y[3];
	    }
	    else{
	      *iter= 5;
	      for(i= 0;i<4;i++){
		x[i]= v[i].x;
		y[i]= v[i].y;
		f[i]= v[i].f;
		g[i]= v[i].g;
	      }
	    }
	  }
	}
      }
    }
    printf("ifg=%2d\n",ifg);
    printf("flags=%2d %2d %2d %2d %2d\n",
	   v[0].fl,v[1].fl,v[2].fl,v[3].fl,v[4].fl);
    if(*iter<5)return(0);
  }
  if(*iter==3){
    if(inst){
      sx[0]= x[2];
      sy[0]= y[2];
      v[4].fl= 0;
      x[2]= EZcr2*(x[0]+x[2]);
      (*iter)--;
      *x0= x[2];
      *y0= y[2];
    }
    else{
      if(f0>=0.)
	ifg= g0>=0?0:3;
      else
	ifg= g0>=0?1:2;
      if(v[ifg].fl){
	v[ifg].fl= 0;
	v[ifg].f= f0;
	v[ifg].g= g0;
	v[ifg].x= *x0;
	v[ifg].y= *y0;
      }
      else{
	if(f0*f0+g0*g0<v[ifg].f*v[ifg].f+v[ifg].g*v[ifg].g){
	  v[ifg].fl= 0;
	  v[ifg].f= f0;
	  v[ifg].g= g0;
	  v[ifg].x= *x0;
	  v[ifg].y= *y0;
	}
      }
      if(ifg==0||ifg==3){
	(*iter)--;
	if(v[4].fl){
	  x[2]-= EZcr2*(x[0]-x[2]);;
	}
	else{
	  x[2]= EZcr2*(sx[0]+x[2]);
	  y[2]= EZcr2*(sy[0]+y[2]);
	}
	*x0= x[2];
	*y0= y[2];
      }
      else{
	if(v[0].fl||v[3].fl){
	  *iter= 3;
	  x[1]= EZcr2*(x[0]+x[2]);
	  *x0= x[1];
	  *y0= y[1];
	}
	else{
	  if(v[1].fl||v[2].fl){
	    *iter= 4;
	    x[3]= EZcr2*(x[0]+x[2]);
	    *x0= x[3];
	    *y0= y[3];
	  }
	  else{
	    *iter= 5;
	    for(i= 0;i<4;i++){
	      x[i]= v[i].x;
	      y[i]= v[i].y;
	      f[i]= v[i].f;
	      g[i]= v[i].g;
	    }
	  }
	}
      }
    }
    printf("ifg=%2d\n",ifg);
    printf("flags=%2d %2d %2d %2d %2d\n",
	   v[0].fl,v[1].fl,v[2].fl,v[3].fl,v[4].fl);
    if(*iter<5)return(0);
  }
  if(*iter==4){
    if(inst){
      sx[0]= x[1];
      sy[0]= y[1];
      v[4].fl= 0;
      x[1]= EZcr2*(x[0]+x[1]);
      (*iter)--;
      *x0= x[1];
      *y0= y[1];
    }
    else{
      if(f0>=0.){
	ifg= g0>=0?0:3;
	if(v[ifg].fl==0){
	  if(fabs(v[ifg].g)<1e-7*fabs(g0)){
	    i= ifg==0?3:0;
	    v[i].fl= v[ifg].fl;
	    v[i].f= v[ifg].f;
	    v[i].g= v[ifg].g;
	    v[i].x= v[ifg].x;
	    v[i].y= v[ifg].y;
	    v[ifg].fl= 1;
	  }
	  if(fabs(g0)<1e-6*fabs(v[ifg].g))
	    ifg= ifg==0?3:0;
	  else{
	    g[1]= v[ifg].g;
	    y[1]= v[ifg].y;
	  }
	}
      }
      else
	ifg= g0>=0?1:2;
      if(v[ifg].fl){
	v[ifg].fl= 0;
	v[ifg].f= f0;
	v[ifg].g= g0;
	v[ifg].x= *x0;
	v[ifg].y= *y0;
      }
      else{
	if(f0*f0+g0*g0<v[ifg].f*v[ifg].f+v[ifg].g*v[ifg].g){
	  v[ifg].fl= 0;
	  v[ifg].f= f0;
	  v[ifg].g= g0;
	  v[ifg].x= *x0;
	  v[ifg].y= *y0;
	}
      }
      if(ifg==1||ifg==2){
	(*iter)--;
	x[1]= EZcr2*(x[1]+x[0]);
	y[1]= EZcr2*(y[1]+y[0]);
	*x0= x[1];
	*y0= y[1];
      }
      else{
	if(v[0].fl||v[3].fl){
	  (*iter)--;
	  b1= fabs(g[1])<fabs(v[ifg].g)?-g[1]:-v[ifg].g;
	  y[1]= v[ifg].y+(v[ifg].y-y[1])*(b1-v[ifg].g)/(v[ifg].g-g[1]);
	  *x0= x[1];
	  *y0= y[1];
	}
	else{
	  if(v[1].fl||v[2].fl){
	    *iter= 4;
	    x[3]= EZcr2*(x[0]+x[2]);
	    *x0= x[3];
	    *y0= y[3];
	  }
	  else{
	    *iter= 5;
	    for(i= 0;i<4;i++){
	      x[i]= v[i].x;
	      y[i]= v[i].y;
	      f[i]= v[i].f;
	      g[i]= v[i].g;
	    }
	  }
	}
      }
    }
    printf("ifg=%2d\n",ifg);
    printf("flags=%2d %2d %2d %2d %2d\n",
	   v[0].fl,v[1].fl,v[2].fl,v[3].fl,v[4].fl);
    if(*iter<5)return(0);
  }
  if(*iter==5&&(v[1].fl||v[2].fl)){
    if(inst){
      (*iter)--;
      sx[0]= x[3];
      sy[0]= y[3];
      v[4].fl= 0;
      x[3]= EZcr2*(x[0]+x[3]);
      *x0= x[3];
      *y0= y[3];
    }
    else{
      if(f0>=0.)
	ifg= g0>=0?0:3;
      else{
	ifg= g0>=0?1:2;
	if(v[ifg].fl==0){
	  if(fabs(v[ifg].g)<1e-7*fabs(g0)){
	    i= ifg==1?2:1;
	    v[i].fl= v[ifg].fl;
	    v[i].f= v[ifg].f;
	    v[i].g= v[ifg].g;
	    v[i].x= v[ifg].x;
	    v[i].y= v[ifg].y;
	    v[ifg].fl= 1;
	  }
	  if(fabs(g0)<1e-6*fabs(v[ifg].g))
	    ifg= ifg==1?2:1;
	  else{
	    g[3]= v[ifg].g;
	    y[3]= v[ifg].y;
	  }
	}
      }
      if(v[ifg].fl){
	v[ifg].fl= 0;
	v[ifg].f= f0;
	v[ifg].g= g0;
	v[ifg].x= *x0;
	v[ifg].y= *y0;
      }
      else{
	if(f0*f0+g0*g0<v[ifg].f*v[ifg].f+v[ifg].g*v[ifg].g){
	  v[ifg].fl= 0;
	  v[ifg].f= f0;
	  v[ifg].g= g0;
	  v[ifg].x= *x0;
	  v[ifg].y= *y0;
	}
      }
      if(ifg==0||ifg==3){
	(*iter)--;
	x[3]= EZcr2*(x[3]+x[2]);
	y[3]= EZcr2*(y[3]+y[2]);
	*x0= x[3];
	*y0= y[3];
      }
      else{
	if(v[1].fl||v[2].fl){
	  (*iter)--;
	  b1= fabs(g[3])<fabs(v[ifg].g)?-g[3]:-v[ifg].g;
	  y[3]= v[ifg].y+(v[ifg].y-y[3])*(b1-v[ifg].g)/(v[ifg].g-g[3]);
	  *x0= x[3];
	  *y0= y[3];
	}
	else{
	  *iter= 5;
	  for(i= 0;i<4;i++){
	    x[i]= v[i].x;
	    y[i]= v[i].y;
	    f[i]= v[i].f;
	    g[i]= v[i].g;
	  }
	}
      }
    }
    printf("ifg=%2d\n",ifg);
    printf("flags=%2d %2d %2d %2d %2d\n",
	   v[0].fl,v[1].fl,v[2].fl,v[3].fl,v[4].fl);
    if(*iter<5)return(0);
  }
  printf("Set of points Iter=%2d\n",*iter);
  printf("0 %11.4e%11.4e %11.4e%11.4e\n",f[0],g[0],x[0],y[0]);
  printf("1 %11.4e%11.4e %11.4e%11.4e\n",f[1],g[1],x[1],y[1]);
  printf("2 %11.4e%11.4e %11.4e%11.4e\n",f[2],g[2],x[2],y[2]);
  printf("3 %11.4e%11.4e %11.4e%11.4e\n",f[3],g[3],x[3],y[3]);
  if(inst){
    *x0= EZcr2*((*x0)+x[0]);
    *y0= EZcr2*((*y0)+y[0]);
    (*iter)--;
    return(0);
  }
  b1= 0;
  for(i= 0;i<4;i++){
    b2= f[i]*f[i]+g[i]*g[i];
    if(b2>b1){
      b1= b2;
      iworst= i;
    }
  }
  b1= 0;
  for(i= 0,j= 0;i<4;i++){
    if(i-iworst){
      xx[j]= x[i];
      yy[j]= y[i];
      ff[j]= f[i];
      gg[j]= g[i];
      j++;
      b2= f[i]*f[i]+g[i]*g[i];
      if(b2>b1){
	b1= b2;
	iworst1= i;
      }
    }
  }
  if(*iter==5){
    guess2D(xx,yy,ff,gg,x0,y0,&dxf,&dxg,&dyf,&dyg);
    printf("new %11.4e%11.4e %1d\n",*x0,*y0,i1);
    return(0);
  }
  b1= 0;
  for(i= 0;i<3;i++){
    b2= ff[i]*ff[i]+gg[i]*gg[i];
    if(b2>b1){
      b1= b2;
      j= i;
    }
  }
  xx[j]= *x0;
  yy[j]= *y0;
  ff[j]= f0;
  gg[j]= g0;
  if(f0>0.){
    ifg= g0>=0?0:3;
  }
  else{
    if(f0==0.){
      ifg= g0>=0?0:3;
    }
    else{
      ifg= g0>=0?1:2;
    }
  }
  i= ifg;
  b2= (f[i]-f0)*(f[i]-f0)+(g[i]-g0)*(g[i]-g0);
  x[ifg]= *x0;
  y[ifg]= *y0;
  f[ifg]= f0;
  g[ifg]= g0;
  guess2D(xx,yy,ff,gg,x0,y0,&dxf,&dxg,&dyf,&dyg);
  i= iworst;
  b1= (f[i]-f0)*(f[i]-f0)+(g[i]-g0)*(g[i]-g0);
  if(ifg-iworst&&b1>16.*b2){
    b1= EZcr4*f[iworst];
    b2= EZcr4*g[iworst];
    *x0= (*x0)+dxf*b1+dxg*b2;
    *y0= (*y0)+dyf*b1+dyg*b2;
    printf("attempt %1d\n",iworst);
  }
  printf("new %11.4e%11.4e %1d %1d\n",*x0,*y0,i1,j);
  b1= (x[0]-x[2])/dx;
  b2= (y[0]-y[2])/dy;
  dist= b1*b1+b2*b2;
  printf("ifg=%1d disc1=%11.4e ",ifg,dist);
  if(dist<eps*eps)(*iter)= -1;
  b1= (x[1]-x[3])/dx;
  b2= (y[1]-y[3])/dy;
  dist= b1*b1+b2*b2;
  printf("disc2=%11.4e\n",dist);
  if(dist<eps*eps)(*iter)= -1;
  if(dist>1e+4)(*iter)= -2;
  return(0);
}

int ReadStrData(line,yd,nd)
     char*line;
     float*yd;
     int nd;
{
  int i,k;
  i= 0;
  while(isspace(line[i])){
    if(line[i]=='\n')
      return(0);
    i++;
  }
  while(isdigit(line[i]))
    i++;
  if(line[i]=='\n')
    return(0);
  if(!isspace(line[i])&&line[i]-'.')
    return(0);
  if(line[i]=='.')
    i= 0;
  k= 0;
  while(k<nd){
    if(sscanf(line+i,"%g",yd+k)){
      k++;
      while(isspace(line[i])){
	if(line[i]=='\n')
	  return(k);
	i++;
      }
      while(!isspace(line[i]))
	i++;
      if(line[i]=='\n')
	return(k);
    }
    else
      return(k);
  }
  return(k);
}

int FailureAlarm(char*p,char*RoutineName)
{
  if(p==NULL){
    printf("Failure in routine %s\n",RoutineName);
#ifdef H
    exit(0);
#endif
  }
  return(0);
}

int CbStrCmp(char*str,char*ln)
{
  while(isspace(*str)){
    str++;
  }
  while(isspace(*ln)){
    ln++;
  }
  while(*str!='\0'&&*str==*ln){
    str++;
    ln++;
  }
  if(*str != *ln){
    if(*str == '\0' && isspace(*ln)){
      while(isspace(*ln))
	ln++;
      if(*ln=='\0'){
	return(0);
      }
      else{
	return(1);
      }
    }

    if(*ln=='\0'&&isspace(*str)){
      while(isspace(*str))
	str++;
      if(*str=='\0'){
	return(0);
      }
      else{
	return(1);
      }
    }
    return(1);
  }
  else{
    return(0);
  }
}

int CbStrCpy0(char*str,char*ln)
{
  char *lc;
  while(isspace(*ln)){
    ln++;
  }
  lc	=str;
  while(*ln != '\0'){
    *lc++= *ln++;
  }
  *lc= '\0';
  return((int)(lc-str));
}

int CbStrCpy(char*str,char*ln)
{
  char *lc;
  while(isspace(*ln)){
    ln++;
  }
  lc	=str;
  while(*ln != '\0' && isspace(*ln) == 0){
    *lc++= *ln++;
  }
  *lc	='\0';
  return((int)(lc-str));
}

int CbStrCpy1(char*str,char*ln)
{
  int i,L;
  char *p;

  while(isspace(*ln)){
    ln++;
  }
  L	=0;
  while(*ln != '\0'){
    if(isspace(*ln)){
      if(*(str-1) != ' '){
	*str++	=' ';
	L++;
      }
    }
    else{
      *str++	=*ln;
      L++;
    }
    ln++;
  }
  if(*(str-1) == ' '){
    *str--;
    L--;
  }
  *str= '\0';
  return(L);
}

int CbStrNCpy1(char*str,char*ln,int n)
{
  int i,L;
  char *p;

  while(isspace(*ln)) 
    ln++;
  L	=0;
  while(*ln != '\0' && L < n){
    if(isspace(*ln)){
      if(*(str-1) != ' '){
	*str++	=' ';
	L++;
      }
      ln++;
    }
    else{
      *str++	=*ln++;
      L++;
    }
  }
  if(*(str-1) == ' '){
    *str--;
    L--;
  }
  *str= '\0';
  return(L);
}

#ifdef H
int CbTxtCmp(char*str,char*ln)
{
  while(*str != '\0'){
    while(isspace(*str)){
      str++;
    }
    while(isspace(*ln)){
      ln++;
    }
    while(*str !='\0' && *str == *ln){
      str++;
      ln++;
    }
    if(*str != *ln){
      if(*str == '\0' && isspace(*ln)){
	while(isspace(*ln)){
	  ln++;
	}
	if(*ln=='\0'){
	  return(0);
	}
	else{
	  return(1);
	}
      }

      if(*ln == '\0' && isspace(*str)){
	while(isspace(*str)){
	  str++;
	}
	if(*str=='\0'){
	  return(0);
	}
	else{
	  return(1);
	}
      }
      return(1);
    }
    else{
      return(0);
    }
  }
}
#endif
int EZppp()
{
  printf("???\n");
  return(0);
}

int EZpppp()
{
  printf("????\n");
  return(0);
}

int EZppppp()
{
  printf("?????\n");
  return(0);
}

int EZpppppp()
{
  printf("??????\n");
  return(0);
}

int EZppppppp()
{
  printf("???????\n");
  return(0);
}

#ifndef stg_INIT
static long int Mmem=0;
int ReInitArray(void **p, int N, int n, int L)
{
  if(N < n){
    if(N > 0){
      *p=(void*)realloc((void*)(*p),(size_t)(n*L));
    }
    else{
      *p	=(void*)malloc((size_t)(n*L));
    }
    if(*p == NULL){
      FailureAlarm((char*)(*p),"ReInitAray()- no memory for *p");
      printf("number allocated=%d number requested=%d, unit size=%d\n",N,n,L);
      printf("Total allocated memory=%ld requested=%ld\n",Mmem,Mmem+n*L);
      exit(0);
      return(-1);
    }
    Mmem	+=n*L;
  }
  return(0);
}

int GetMemAmmount(long int *M)
{
  *M	=Mmem;
  return(0);
}

#endif
#endif

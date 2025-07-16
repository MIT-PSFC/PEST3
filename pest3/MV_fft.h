#ifndef __MV_fft__
#define __MV_fft__


#include "MV_vector.h"
#include "MV_matrix.h"

extern "C" int rffti(int n, double wspace[] );
extern "C" int rfftf(int n, double r[], double wspace[] );
extern "C" int rfftb(int n, double r[], double wspace[] );

typedef MV_ColMat<double> Mat;
typedef MV_Vector<double> Vec;
typedef MV_Vector<int> Vec_int;


/*

A. Pletzer February 1998
updated on Sept 8 98

fft + Cos and Sin coefficients

*/


/**************************************************************************

   Fast Fourier Transform

   If x is the sampling input vector of length N with spacing dt then:

                N/2
    x(n) = a0 + sum a(k)*cos(2*pi*k*t(n)/(N*dt))+b(k)*sin(2*pi*k*t(n)/(N*dt))
                k=1
   
    and

    a0    =   X(0)
    a(k+1)=   X(2*k+1)
    b(k+1)=  -X(2*k+2)

    where 

    X(k) is the vector return by the fft call. The sampling can take
    place along the x or y direction according to whether the fft0 or
    fft1 routines are called.

    Note that the backward ifft's are normalized (by 1/N) so that the
    sequential call of ifft0(fft0(a)) produces the matrix a. 
    
    ***********************************************************************/

/* vector FFT */

Vec fft(const Vec &a) {
  unsigned int i;
  unsigned int an = a.size();
  Vec b(an);

  double *wspace = new double[2*an+15];
  double *f = new double[an];

  rffti(an, wspace); // initialization

  for( i=0; i<an; ++i) *(f+i) = a(i);
  rfftf(an, f, wspace);
  for( i=0; i<an; ++i) b(i) = *(f+i);

  return b;
}
      

/* Fourier transform along rows (SECOND index) */

Mat fft1(const Mat &a) {
  double *wspace;
  double *f;
  unsigned int anr = a.size(0);
  unsigned int anc = a.size(1);
  Mat b(anr, anc);

  wspace = new double[2*anc+15];
  f = new double[anc];

  // initialization

  rffti(anc, wspace);

  Vec x(anc), bt(anc);
  Vec_int J = range(0, anc);
  for(unsigned int i=0; i<anr; ++i){
    x = a(i, J);
    bt = fft( x );
    for (unsigned int j=0; j<anc; ++j) b(i,j) = bt(j);
  }
  return b;
}



/* Fourier transform along columns (FIRST index) direction */

Mat fft0(const Mat &a) {
  double *wspace;
  double *f;
  unsigned int anr = a.size(0);
  unsigned int anc = a.size(1);
  Mat b(anr, anc);

  wspace = new double[2*anr+15];
  f = new double[anr];

  // initialization

  rffti(anr, wspace);

  Vec x(anr), bt(anr);
  Vec_int I=range(0, anr);
  for(unsigned int j=0; j<anc; ++j){
    x = a(I, j);
    bt = fft( x );
    for(unsigned int i=0; i<anr; ++i) b(i,j) = bt(i);
  }
  return b;
}
      


/* Inverse vector FFT/ normalized by vector length */

Vec ifft(const Vec &a) {
  unsigned int i;
  unsigned int an = a.size();
  Vec b(an);

  double *wspace = new double[2*an+15];
  double *f = new double[an];

  rffti(an, wspace); // initialization

  for(i=0; i<an; ++i) *(f+i) = a(i);
  rfftb(an, f, wspace);
  for(i=0; i<an; ++i) b(i) = *(f+i)/double(an);

  return b;
}

           

/* Inverse Fourier transform along columns (SECOND index) */

Mat ifft1(const Mat &a) {
  double *wspace;
  double *f;
  unsigned int anr = a.size(0);
  unsigned int anc = a.size(1);
  Mat b(anr, anc);

  wspace = new double[2*anc+15];
  f = new double[anc];

  // initialization

  rffti(anc, wspace);

  Vec x(anc), bt(anc);
  Vec_int J = range(0, anc);
  for(unsigned int i=0; i<anr; ++i){
    x = a(i,J);
    bt = ifft( x );
    for( unsigned int j=0; j<anc; ++j) b(i,j) = bt(j);
  }
  return b;
}



/* Inverse Fourier transform along rows (FIRST index) direction */

Mat ifft0(const Mat &a) {
  double *wspace;
  double *f;
  unsigned int anr = a.size(0);
  unsigned int anc = a.size(1);
  Mat b(anr, anc);

  wspace = new double[2*anr+15];
  f = new double[anr];

  // initialization

  rffti(anr, wspace);

  Vec x(anr), bt(anr);
  Vec_int I = range(0, anr);
  for(unsigned int j=0; j<anc; ++j){
    x = a(I,j);
    bt = ifft( x );
    for(unsigned int i=0; i<anr; ++i) b(i,j) = bt(i);
  }
  return b;
}
      
/* cos coefficients */

Vec fcos(const Vec &a){
  unsigned int an = a.size();
  Vec b = fft(a);
  Vec c(an/2+1);
  c(0) = b(0)/double(an);
  for(unsigned int i=1; i<an/2; ++i) c(i) = 2.0*b(2*i-1)/double(an);
  c(an/2) = b(an-1)/double(an);
  return c;
}
  

/* get cosine and sine coefficients from matrices returned by fft0 and 
   fft1 */

Mat fcos0(const Mat &a) {

  unsigned int anr = a.size(0);
  unsigned int anc = a.size(1);
  Mat b(anr, anc/2+1, 0.0);

  for(unsigned int i=0; i<anr; ++i) {
    b(i,0) = a(i,0);
    for(unsigned int j=1; j<anc/2+1; ++j) {
      b(i,j) = a(i,2*j-1);
    }
  }
  return b;
}

/* sin coefficients */

Vec fsin(const Vec &a){
  unsigned int an = a.size();
  Vec b = fft(a);
  Vec c(an/2+1);
  c(0) = 0.0;
  for(unsigned int i=1; i<an/2; ++i) c(i) = -2.0*b(2*i)/double(an);
  c(an/2) = 0.0;
  return c;
}

Mat fcos1(const Mat &a) {

  unsigned int anr = a.size(0);
  unsigned int anc = a.size(1);
  Mat b(anr/2+1, anc, 0.0);

  for(unsigned int j=0; j<anc; ++j) {
    b(0,j) = a(0,j);
    for(unsigned int i=1; i<anr/2+1; ++i) {
      b(i,j) = a(2*i-1,j);
    }
  }
  return b;
}

Mat fsin0(const Mat &a) {

  unsigned int anr = a.size(0);
  unsigned int anc = a.size(1);
  Mat b(anr, anc/2+1, 0.0);

  for(unsigned int i=0; i<anr; ++i) {
    b(i,0) = 0.;
    for(unsigned int j=1; j<anc/2; ++j) {
      b(i,j) = -a(i,2*j);
    }
    if( 2*(anc/2) != anc) b(i,anc/2) = -2.*a(i,anc-1); // little fix
  }
  return b;
}

Mat fsin1(const Mat &a) {

  unsigned int anr = a.size(0);
  unsigned int anc = a.size(1);
  Mat b(anr/2+1, anc, 0.0);

  for(unsigned int j=0; j<anc; ++j) {
    b(0,j) = 0.;
    for(unsigned int i=1; i<anr/2; ++i) {
      b(i,j) = -a(2*i,j);
    }
    if( 2*(anr/2) != anr) b(anr/2,j) = -2.*a(anr-1,j); // little fix
  }
  return b;
}


#endif /* __MV_fft__ */


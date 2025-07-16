#include <stdlib.h>
#ifndef mzl_PL
#define ESND 257
#define ESNPRPR 13
#define ESNCRPR 18
#define ESNRDCR 4
typedef struct PlProf{
  int kk[ESND];
  double d0[ESND];
  double d1[ESND];
  double d2[ESND];
  double xd[ESND];
  double yd[ESND];
  double *x0;
  double *y0;
  double *y1;
  double *y2;
  int io;
  int nd;
  int xin;
  char Nm[8];
}PlProf;

typedef struct RadCoord{
  int kk[ESND];
  double ad[ESND];
  double xd[ESND];
  double *x0;
  double *x1;
  double *x2;
  double *a0;
  double *a1;
  double *a2;
  double X;
  int io;
  int nd;
  char Nm[8];
}RadCoord;

typedef struct EqC{
  int Na;
  int Np;
  int Fp;
  int Ina;
  int Inp;
  int Inj;
  int Ndp;
  int Ndj;
  double *a;
  double *b;
  double *b1a;
  double *b2a;

  double *cR;
  double *cR1a;
  double *cR2a;


  double *xdp;
  double *ydp;
  double *xdj;
  double *ydj;

}EqC;

#endif

/**
 * @file 
 * @brief List of Get methods. 
 * The following routines give access equilibrium profile data
 * and geometry after execution (Esc).
 */ 

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

extern int ESNp1, ESNa1, ESNa;

const double ESC_TwoPi = 6.2831853071796;
 
/**
 * Get the radial number of grid points used internally by ESC. 
 * @param Na1 the grid size
 * @author L.E. Zakharov, A. Pletzer
 * @see escGetNp1
 * 
 */
void escGetNa1(int *Na1){*Na1 = ESNa1;}
void escgetna1_(int *Na1){escGetNa1(Na1);}
void escgetna1(int *Na1){escGetNa1(Na1);}
void ESCGETNA1(int *Na1){escGetNa1(Na1);}
 
/**
 * Get the number of poloidal rays + 1 used internally by ESC. 
 * @param Np1 the number of poloidal sections + 1
 * @author L.E. Zakharov, A. Pletzer
 * @see escGetNa1
 */
void escGetNp1(int *Np1){*Np1 = ESNp1;};
void escgetnp1_(int *Np1){escGetNp1(Np1);}
void escgetnp1(int *Np1){escGetNp1(Np1);}
void ESCGETNP1(int *Np1){escGetNp1(Np1);}
 
/**
 * Get the poloidal flux/(2*pi) [Wb/rad].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer
 * @see escGetPsi
 */
void escGetPsibar(double *prof, int *ns){
  extern double *ESgY, *ESgY2a, *ESsa;
  int i;
  double dx, x;
  dx =ESsa[ESNa]/(*ns-1);
  for(i=0; i < *ns; i++){
      x		=dx*i;
      ESSetSplA(x);
      splRA(prof+i,NULL,ESgY,ESgY2a);
      prof[i] = -prof[i];
  }
}
void escgetpsibar_(double *prof, int *ns){escGetPsibar(prof, ns);}
void escgetpsibar(double *prof, int *ns){escGetPsibar(prof, ns);}
void ESCGETPSIBAR(double *prof, int *ns){escGetPsibar(prof, ns);}
 
/**
 * Get the poloidal flux [Wb].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer
 * @see escGetPsibar
 */
void escGetPsi(double *prof, int *ns){
  int i;
  escGetPsibar(prof, ns);
  for(i=0; i < *ns; i++) prof[i] *= ESC_TwoPi;
}
void escgetpsi_(double *prof, int *ns){escGetPsi(prof, ns);}
void escgetpsi(double *prof, int *ns){escGetPsi(prof, ns);}
void ESCGETPSI(double *prof, int *ns){escGetPsi(prof, ns);}
 

/**
 * Get the toroidal flux/(2*pi) [Wb/rad].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer
 * @see escGetPhi
 */
void escGetPhibar(double *prof, int *ns){
  extern double *ESgF, *ESgF2a, *ESsa;
  int i;
  double dx, x;
  dx =ESsa[ESNa]/(*ns-1);
  for(i=0; i < *ns; i++){
      x		=dx*i;
      ESSetSplA(x);
      splRA(prof+i,NULL,ESgF,ESgF2a);
  }
}
void escgetphibar_(double *prof, int *ns){escGetPhibar(prof, ns);}
void escgetphibar(double *prof, int *ns){escGetPhibar(prof, ns);}
void ESCGETPHIBAR(double *prof, int *ns){escGetPhibar(prof, ns);}
 
/**
 * Get the toroidal flux [Wb].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetPhibar
 */
void escGetPhi(double *prof, int *ns){
  int i;
  escGetPhibar(prof, ns);
  for(i=0; i < *ns; i++) prof[i] *= ESC_TwoPi;
}
void escgetphi_(double *prof, int *ns){escGetPhi(prof, ns);}
void escgetphi(double *prof, int *ns){escGetPhi(prof, ns);}
void ESCGETPHI(double *prof, int *ns){escGetPhi(prof, ns);}
 
/**
 * Get the covariant toroidal magnetic field function [Tm].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetIota escGetQ escGetP
 */
void escGetG(double *prof, int *ns){
  extern double *ESaF, *ESaF2a, *ESsa;
  int i;
  double dx, x;
  dx =ESsa[ESNa]/(*ns-1);
  for(i=0; i < *ns; i++){
      x		=dx*i;
      ESSetSplA(x);
      splRA(prof+i,NULL,ESaF,ESaF2a);
  }
}
void escgetg_(double *prof, int *ns){escGetG(prof, ns);}
void escgetg(double *prof, int *ns){escGetG(prof, ns);}
void ESCGETG(double *prof, int *ns){escGetG(prof, ns);}

/**
 * Get the 1/q profile [-].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetQ escGetG escGetP
 */ 
void escGetIota(double *prof, int *ns){
  extern double *ESgm, *ESgm2a, *ESsa;
  int i;
  double dx, x;
  dx =ESsa[ESNa]/(*ns-1);
  for(i=0; i<ESNa1; i++){
  }

  for(i=0; i < *ns; i++){
      x		=dx*i;
      ESSetSplA(x);
      splRA(prof+i,NULL,ESgm,ESgm2a);
  }
}
void escgetiota_(double *prof, int *ns){escGetIota(prof, ns);}
void escgetiota(double *prof, int *ns){escGetIota(prof, ns);}
void ESCGETIOTA(double *prof, int *ns){escGetIota(prof, ns);}

/**
 * Get the safety factor (q) profile [-].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetIota escGetG escGetP
 */ 
void escGetQ(double *prof, int *ns){
    int i;
    escGetIota(prof, ns);
    for(i=0; i < *ns; i++) prof[i] =  1.0/prof[i];
}
void escgetq_(double *prof, int *ns){escGetQ(prof, ns);}
void escgetq(double *prof, int *ns){escGetQ(prof, ns);}
void ESCGETQ(double *prof, int *ns){escGetQ(prof, ns);}

/**
 * Get the pressure in [mu0 Pa].
 * @param prof the returned profile
 * @param ns the size of the returned profile (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetIota escGetG escGetQ
 */
void escGetP(double *prof, int *ns){
  extern double *ESsp, *ESsp2a, *ESsa;
  int i;
  double dx, x;
  dx =ESsa[ESNa]/(*ns-1);
  for(i=0; i < *ns; i++){
      x		=dx*i;
      ESSetSplA(x);
      splRA(prof+i,NULL,ESsp,ESsp2a);
  }
}
void escgetp_(double *prof, int *ns){escGetP(prof, ns);}
void escgetp(double *prof, int *ns){escGetP(prof, ns);}
void ESCGETP(double *prof, int *ns){escGetP(prof, ns);}
 
/** 
 * Get the Fourier coefficients of the R coordinate [m].
 * @param rcos the returned cosine coefficients of sizes [ns][nf]
 * @param rsin the returned sine coefficients
 * @param nf the number of Fourier coefficients (input)
 * @param ns the radial size (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetFourierZ
 */
void escGetFourierR(double *rcos, double *rsin, int *nf, int *ns){
  extern double *ESsa;
  extern double *rcT, *rcT2a;
  extern double *rsT, *rsT2a;
  extern double *ESaR0;
  extern int ESFp1;

  int i, j, nfmax, mi, ni;
  double dx, x, s;

  nfmax = (*nf > ESFp1) ? ESFp1 : *nf;

  dx =ESsa[ESNa]/(*ns-1);
  for(i=0; i < *ns; i++){
      x		=dx*i;
      ESSetSplA(x);
      j = 0;
      mi  = j * ESNa1;
      ni  = i* (*nf) + j;
      splRA(&rcos[ni], NULL, rcT+mi, rcT2a+mi);
      rcos[ni] += ESaR0[0];
      rsin[ni] = 0.0;
      
      for(j=1; j < nfmax; j++){
	mi  = j * ESNa1;
	ni  = i* (*nf) + j;
	splRA(&rcos[ni], NULL, rcT+mi, rcT2a+mi);
	splRA(&rsin[ni], NULL, rsT+mi, rsT2a+mi);
	rcos[ni] *= 2.0;
	rsin[ni] *= 2.0;
      }
      for(j=nfmax; j < *nf; j++){
	rcos[ni] = 0.0;
	rsin[ni] = 0.0;
      }
  }
}
void escgetfourierr_(double *rcos, double *rsin, int *nf, int *ns){
  escGetFourierR(rcos, rsin, nf, ns);
}
void escgetfourierr(double *rcos, double *rsin, int *nf, int *ns){
  escGetFourierR(rcos, rsin, nf, ns);
}
void ESCGETFOURIERR(double *rcos, double *rsin, int *nf, int *ns){
  escGetFourierR(rcos, rsin, nf, ns);
}


/** 
 * Get the R coordinates [m]. The poloidal index varies faster.
 * @param r the returned array of size [ns][nt1].
 * @param nt1 the number of poloidal sections + 1 (input)
 * @param ns the radial size (input)
 * @param clockwise poloidal angle orientation (1 for clockwise, 
 *        0 for counterclockwise)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetZ
 */
void escGetR(double *r, int *nt1, int *ns, int *clockwise){
  extern int ESFp1;
  double *rcos, *rsin;
  int NT1, NS, i, j, k, ni, nf;
  double sum, theta;
  NT1 = *nt1;
  NS  = *ns;
  nf = ESFp1;

  rcos =(double*) malloc( nf*NS*sizeof(double));
  rsin =(double*) malloc( nf*NS*sizeof(double));

  escGetFourierR(rcos, rsin, &nf, ns);

  if(*clockwise!=1){
    for(i=0; i<NS; i++){
      for(j=0; j<NT1; j++){
	sum = 0.0;
	theta = ESC_TwoPi*(double)j/(double)(NT1-1);
	for(k=0; k<nf; k++){
	  ni = i*nf + k;
	  sum += rcos[ni]*cos(k*theta) + rsin[ni]*sin(k*theta);
	}
	r[i*NT1 + j] = sum;
      }
    }
  }else{
    for(i=0; i<NS; i++){
      for(j=0; j<NT1; j++){
	sum = 0.0;
	theta = ESC_TwoPi*(double)(NT1-1-j)/(double)(NT1-1);
	for(k=0; k<nf; k++){
	  ni = i*nf + k;
	  sum += rcos[ni]*cos(k*theta) + rsin[ni]*sin(k*theta);
	}
	r[i*NT1 + j] = sum;
      }
    }
  }
    

  free(rcos);
  free(rsin);
}
void escgetr_(double *r, int *nt1, int *ns, int *clockwise)
{
  escGetR(r, nt1, ns, clockwise);
}
void escgetr(double *r, int *nt1, int *ns, int *clockwise)
{
  escGetR(r, nt1, ns, clockwise);
}
void ESCGETR(double *r, int *nt1, int *ns, int *clockwise)
{
  escGetR(r, nt1, ns, clockwise);
}
 
/** 
 * Get the Fourier coefficients of the Z coordinate [m].
 * @param zcos the returned cosine coefficients of sizes [ns][nf].
 * @param zsin the returned sine coefficients
 * @param nf the number of Fourier coefficients (input)
 * @param ns the radial size (input)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetFourierR
 */
void escGetFourierZ(double *zcos, double *zsin, int *nf, int *ns){
  extern double *ESsa;
  extern double *ESsb, *ESsb2a, *rsT, *rsT2a;
  extern double *ESaZ0;
  extern int ESFp1;

  int i, j, mi, ni;
  double dx, x, s;

  dx =ESsa[ESNa]/(*ns-1);
  for(i=0; i < *ns; i++){
      x		=dx*i;
      ESSetSplA(x);
      j = 0;
      ni  = i* (*nf) + j;
      /* 0-th Fourier of Z-cos stored in R-sin */
      splRA(&zcos[ni], NULL, rsT, rsT2a);
      zcos[ni] += ESaZ0[0];
      zsin[ni] = 0.0;
      j = 1;
      mi  = j * ESNa1;
      ni  = i* (*nf) + j;
      zcos[ni] = 0.0;
      splRA(&zsin[ni], NULL, ESsb, ESsb2a);
      
      for(j=2; j < *nf; j++){
	ni  = i* (*nf) + j;
	zcos[ni] = 0.0;
	zsin[ni] = 0.0;
      }
  }
}
void escgetfourierz_(double *zcos, double *zsin, int *nf, int *ns){
  escGetFourierZ(zcos, zsin, nf, ns);
}
void escgetfourierz(double *zcos, double *zsin, int *nf, int *ns){
  escGetFourierZ(zcos, zsin, nf, ns);
}
void ESCGETFOURIERZ(double *zcos, double *zsin, int *nf, int *ns){
  escGetFourierZ(zcos, zsin, nf, ns);
}

/** 
 * Get the Z coordinate [m]. The poloidal index varies faster.
 * @param r the returned array of size [ns][nt1].
 * @param nt1 the number of poloidal sections + 1 (input)
 * @param ns the radial size (input)
 * @param clockwise poloidal angle orientation (1 for clockwise, 
 *        0 for counterclockwise)
 * @author L.E. Zakharov, A. Pletzer 
 * @see escGetR
 */
void escGetZ(double *z, int *nt1, int *ns, int *clockwise){
  extern int ESFp1;
  double *zcos, *zsin;
  int NT1, NS, i, j, k, ni, nf;
  double sum, theta;
  NT1 = *nt1;
  NS  = *ns;
  nf = ESFp1;

  zcos =(double*) malloc( nf*NS*sizeof(double));
  zsin =(double*) malloc( nf*NS*sizeof(double));

  escGetFourierZ(zcos, zsin, &nf, ns);

  if(*clockwise!=1){
    for(i=0; i<NS; i++){
      for(j=0; j<NT1; j++){
	sum = 0.0;
	theta = ESC_TwoPi*(double)j/(double)(NT1-1);
	for(k=0; k<nf; k++){
	  ni = i*nf + k;
	  sum += zcos[ni]*cos(k*theta) + zsin[ni]*sin(k*theta);
	}
	z[i*NT1 + j] = sum;
      }
    }
  }else{
    for(i=0; i<NS; i++){
      for(j=0; j<NT1; j++){
	sum = 0.0;
	theta = ESC_TwoPi*(double)(NT1-1-j)/(double)(NT1-1);
	for(k=0; k<nf; k++){
	  ni = i*nf + k;
	  sum += zcos[ni]*cos(k*theta) + zsin[ni]*sin(k*theta);
	}
	z[i*NT1 + j] = sum;
      }
    }
  }    
  free(zcos);
  free(zsin);
}
void escgetz_(double *z, int *nt1, int *ns, int * clockwise){
  escGetZ(z, nt1, ns, clockwise);
}
void escgetz(double *z, int *nt1, int *ns, int * clockwise){
  escGetZ(z, nt1, ns, clockwise);
}
void ESCGETZ(double *z, int *nt1, int *ns, int * clockwise){
  escGetZ(z, nt1, ns, clockwise);
}

/**
 * RGA: additional support routines
 *  ecsavesc() - write out current time.  This is 0.0 before any equilibrium is solved.
 *  ecrestoreesc(time) - read in an equilibrium, return 0 on success.
 *
 */
int ecsaveesc()  { return ECSaveESC() ; } ;
int ecsaveesc_() { return ECSaveESC() ; } ;
int ECSAVEESC()  { return ECSaveESC() ; } ;

int ecrestoreesc(double *time)  { return ECRestoreESC(*time) ; } ;
int ecrestoreesc_(double *time) { return ECRestoreESC(*time) ; } ;
int ECRESTOREESC(double *time)  { return ECRestoreESC(*time) ; } ;


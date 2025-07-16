#ifndef __mapper__
#define __mapper__

#include <string>
#include <stdlib.h>
#include <math.h>
#include "ezcdf.hh"
#include "MV_vector.h"
#include "MV_matrix.h"

/* Mapper class. A. Pletzer Feb-23-2000
   Read various equilibrium output file formats and produce
   metric data on new radial and poloidal mesh. */

extern "C" {

  void dmpinit_(void);

  void dmpsetinputformat_(int *);
  void dmpsetnths0_(int *);
  void dmpsetnsf0_(int *);
  void dmpsetmjac_(int *);
  void dmpsetnp_(int *);
  void dmpsetigrid_(int *);
  
  void dmprinp1b_(int *);
  void dmpreqdsk_(void);
  void dmpreadinp1_ (void);

  void dmpcheck_positiveness_(double *, int *, int *);

  void dmpgrid_(void);
  void dmparrange_(void);
  void dmpnewtheta_(void);
  
  void dmpwmpout1_(void);
  void dmpwmapdsk_(void);

  void dmpgetrscale_(double *);
  void dmpgetpsimin_(double *);
  void dmpgetpsilim_(double *);

  void dmpgetp_(double *);
  void dmpgetpp_(double *);
  void dmpgetq_(double *);
  void dmpgetqp_(double *);
  void dmpgetg_(double *);
  void dmpgetgp_(double *);
  void dmpgetf_(double *);
  void dmpgetfp_(double *);
  void dmpgetx_(double *);
  void dmpgetz_(double *);
  void dmpgetxdth_(double *);
  void dmpgetzdth_(double *);
  void dmpgetxdps_(double *);
  void dmpgetzdps_(double *);

  void dmpgetxinf_(double *);
  void dmpgetzinf_(double *);
  void dmpgetpsibig_(double *);
  void dmpgetxsq_(double *);
  void dmpgetgrpssq_(double *);
  void dmpgetgrthsq_(double *);
  void dmpgetgrpsth_(double *);
  void dmpgetgrptdth_(double *);
  void dmpgetxsqdps_(double *);
  void dmpgetgrpsdth_(double *);
  void dmpgetxsqdth_(double *);
  void dmpgetxjacob_(double *);
  void dmpgetxjprym_(double *);
  void dmpgetdelta_(double *);
  void dmpgetqdelp_(double *);
  

  void dmpfree_(void);

}

const std::string mapper_inputError   = "--Mapper Input Error     --";
const std::string mapper_execError    = "--Mapper Execution Error --";
const std::string mapper_unknownError = "--Mapper Unknown Error   --";


class mapper {
public:
  int mjac, np, nosurf, mth, nsf, nths, ier;
  int nErrors;

  int inputFormat; /*  -1 Chease INP1
		       0 Jsolver netCDF
		       1 Chease netCDF */

  int outputFormat; /* 0 no output files
		       1 save in MAPDSK.cdf & MPOUT1.cdf */

  int igrid; /* 1 ~ sqrt(psi)
		2 ~ psi
		3 original */

  double rscale, psimin, psilim;
  Vec xinf, zinf, psibig;
  Mat xsq, grpssq, grthsq, grpsth, grptdth, 
    xsqdps, grpsdth, xsqdth, xjacob, xjprym, delta, qdelp; 
  Vec p, pp, q, qp, g, gp, f, fp;
  Mat x, z, xdth, zdth, xdps, zdps;

  // c'tor
  mapper(){
    mjac=2;
    np  =0;
    nosurf=601;
    mth   =128;
    nsf   = nosurf + 1;
    nths  = mth + 5;
    inputFormat = 0;
    outputFormat = 1;
    ier = 0;
    nErrors = 0;
    igrid=3;
  } 
    
  // d'tor
  ~mapper(){}
  
  // run

  int run(void);
};

int mapper::run(void){

  /* go */

  dmpinit_();
  
  // set input

  nsf = nosurf + 1;
  nths = mth + 5;

  if(inputFormat==1){

    // adjust poloidal number of rays if necessary
    // input file is inp1.cdf

    EZcdf inp1((std::string)"inp1.cdf", "r");
    int nxx[3];
    inp1.cdfGetVar((std::string) "nxx", &nxx[0]);
    mth = nxx[0];
    nosurf = nxx[1];
  }

  cout << " Mapper input:\n";
  cout << " inputFormat=" << inputFormat << endl;
  cout << " mth        =" << mth << endl;
  cout << " nosurf     =" << nosurf << endl;
  cout << " mjac       =" << mjac << endl;
  cout << " np         =" << np << endl;
  cout << " igrid      =" << igrid << endl;
  dmpsetinputformat_(&inputFormat);
  dmpsetnths0_(&mth);
  dmpsetnsf0_(&nosurf);
  dmpsetmjac_(&mjac);
  dmpsetnp_(&np);
  dmpsetigrid_(&igrid);

  // allocate

//   if(p.size()==0) p.alloc(nosurf);
//   if(pp.size()==0) pp.alloc(nosurf);
//   if(q.size()==0) q.alloc(nosurf);
//   if(qp.size()==0) qp.alloc(nosurf);
//   if(g.size()==0) g.alloc(nosurf);
//   if(gp.size()==0) gp.alloc(nosurf);
//   if(f.size()==0) f.alloc(nosurf);
//   if(fp.size()==0) fp.alloc(nosurf);
//   if(x.size(0)*x.size(1)==0) x.alloc(nths, nsf);
//   if(z.size(0)*z.size(1)==0) z.alloc(nths, nsf);
//   if(xdth.size(0)*xdth.size(1)==0) xdth.alloc(nths, nsf);
//   if(zdth.size(0)*zdth.size(1)==0) zdth.alloc(nths, nsf);
//   if(xdth.size(0)*xdth.size(1)==0) xdth.alloc(nths, nsf);
//   if(zdth.size(0)*zdth.size(1)==0) zdth.alloc(nths, nsf);

//   if(xinf.size()==0) xinf.alloc(mth);
//   if(zinf.size()==0) zinf.alloc(mth);
//   if(psibig.size()==0) psibig.alloc(nosurf);
//   if(xsq.size(0)*xsq.size(1)==0) xsq.alloc(nths, nsf);
//   if(grpssq.size(0)*grpssq.size(1)==0) grpssq.alloc(nths, nsf);
//   if(grthsq.size(0)*grthsq.size(1)==0) grthsq.alloc(nths, nsf);
//   if(grpsth.size(0)*grpsth.size(1)==0) grpsth.alloc(nths, nsf);
//   if(grptdth.size(0)*grptdth.size(1)==0) grptdth.alloc(nths, nsf);
//   if(xsqdps.size(0)*xsqdps.size(1)==0) xsqdps.alloc(nths, nsf);
//   if(grpsdth.size(0)*grpsdth.size(1)==0) grpsdth.alloc(nths, nsf);
//   if(xsqdth.size(0)*xsqdth.size(1)==0) xsqdth.alloc(nths, nsf);
//   if(xjacob.size(0)*xjacob.size(1)==0) xjacob.alloc(nths, nsf);
//   if(xjprym.size(0)*xjprym.size(1)==0) xjprym.alloc(nths, nsf);
//   if(delta.size(0)*delta.size(1)==0) delta.alloc(nths, nsf);
//   if(qdelp.size(0)*qdelp.size(1)==0) qdelp.alloc(nths, nsf);


  // read input file

  switch(inputFormat){
  case -1:
    // Chease INP1
    cout << " Looking for file INP1\n";
    dmprinp1b_(&ier);
    if (ier !=0) return (++nErrors);
    break;
  case 0:
    // Jsolver eqdsk.cdf
    dmpreqdsk_();
    break;
  case 1:
    dmpreadinp1_();
    break;
  default:
    cerr << mapper_inputError << " Invalid input format selection \n";
    return (++nErrors);
  }

  // regrid

  dmpgrid_();
  
  // re-arrange nodes

  dmparrange_();

  // compute new theta angle

  dmpnewtheta_();

  // access remapped data

  dmpgetrscale_(&rscale);
  dmpgetpsimin_(&psimin);
  dmpgetpsilim_(&psilim);

  // Note: the following routines extract and interpolate the
  // onto the new grid. They CANNOT be called in conjunction 
  // the dmpwpout1 and dmpwmapdsk routines!

//    cout << "9 1\n";
//    dmpgetp_(&p(0));
//    cout << "9 2\n";
//    dmpgetpp_(&pp(0));
//    cout << "9 3\n";
//    dmpgetq_(&q(0));
//    cout << "9 4\n";
//    dmpgetqp_(&qp(0));
//    cout << "9 5\n";
//    dmpgetg_(&g(0));
//    cout << "9 6\n";
//    dmpgetgp_(&gp(0));
//    cout << "9 7\n";
//    dmpgetf_(&f(0));
//    cout << "9 8\n";
//    dmpgetfp_(&fp(0));
//    cout << "9 9\n";
//    cout <<" x.size(0)=" << x.size(0) << " x.size(1)=" << x.size(1) << endl;
//    dmpgetx_(&x(0,0));
//    cout << "x=" << x << endl;
//    cout << "9 10\n";
//    dmpgetz_(&z(0,0));
//    cout << "z=" << z << endl;
//   cout << "9 11\n";
//   dmpgetxdth_(&xdth(0,0));
//   cout << "9 12\n";
//   dmpgetzdth_(&zdth(0,0));
//   cout << "9 13\n";
//   dmpgetxdps_(&xdps(0,0));
//   cout << "9 14\n";
//   dmpgetzdps_(&zdps(0,0));

//   cout << "9 15\n";
//   dmpgetxinf_(&xinf(0));
//   cout << "9 16\n";
//   dmpgetzinf_(&zinf(0));
//   cout << "9 17\n";
//   dmpgetpsibig_(&psibig(0));
//   cout << "9 18\n";
//   dmpgetxsq_(&xsq(0,0));
//   cout << "9 19\n";
//   dmpgetgrpssq_(&grpssq(0,0));
//   cout << "9 20\n";
//   dmpgetgrthsq_(&grthsq(0,0));
//   cout << "9 21\n";
//   dmpgetgrpsth_(&grpsth(0,0));
//   cout << "9 22\n";
//   dmpgetgrptdth_(&grptdth(0,0));
//   cout << "9 23\n";
//   dmpgetxsqdps_(&xsqdps(0,0));
//   cout << "9 24\n";
//   dmpgetgrpsdth_(&grpsdth(0,0));
//   cout << "9 25\n";
//   dmpgetxsqdth_(&xsqdth(0,0));
//   cout << "9 26\n";
//   dmpgetxjacob_(&xjacob(0,0));
//   cout << "9 27\n";
//   dmpgetxjprym_(&xjprym(0,0));
//   cout << "9 28\n";
//   dmpgetdelta_(&delta(0,0));
//   cout << "9 29\n";
//   dmpgetqdelp_(&qdelp(0,0));

  if(outputFormat==1){
    // interpolate and save data into 
    // MPOUT1 and MAPDSK files. The order
    // here matters. 
    cout << " Saving metric data in files MPOUT1.cdf and MAPDSK.cdf \n";
    dmpwmpout1_();
    dmpwmapdsk_();
  }
  
  // that's it, bye...

  dmpfree_();

  return (nErrors);
}
  

#endif // __mapper__

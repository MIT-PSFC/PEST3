#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "pest3.hh"

using std::cout;
using std::endl;

// driver for pest3 

int main(int narg, char *argc[]){

  pest3 p;
  p.inputFormat = 1;
  int i;

  if( p.parseArgs(narg, argc) !=0 ){
    return (p.nErrors);
  }
  if( p.checkArgs() !=0 ){
    return (p.nErrors);
  }

  cout << " PEST3\n";
  cout << " resonant poloidal mode    ms = " << p.ms << endl;
  cout << " toroidal mode             n  = " << p.n  << endl;
  cout << " number of finite elements mm = ";
  for(i=0; i<pest3_kmax; i++) {
    if(p.mm[i]!=0) cout << p.mm[i] << ' ';
  }
  cout << endl;
  cout << " Fourier modes " << -p.lmax1 << ", ..., " << p.lmax1 << endl;
  cout << " normalized wall distance  b  = " << p.b << endl;

  if(p.run() != 0){
    return p.nErrors;
  }

  cout << " ==========================================================\n";
  cout << " equilibrium grid dimension " << p.mth << '*' << p.nosurf << endl;
  cout << " minor/major radius     = " << p.minorRadius << '/' << 
    p.majorRadius << " =" << p.minorRadius/p.majorRadius << endl;
  cout << " R-magnetic             = " << p.rMagnetic << endl; 
  cout << " elongation             = " << p.elongation << endl;
  cout << " triangularity          = " << p.triangularity << endl;
  cout << " separatrix             = " << p.separatrix << endl;
  cout << " q axis/edge            = " << p.q0   << '/' << p.qa   << endl;
  cout << " q min/max              = " << p.qmin << '/' << p.qmax << endl;
  cout << " beta-poloidal          = " << p.betaPoloidal << endl;
  cout << " beta-toroidal          = " << p.betaToroidal << endl;
  cout << " beta-N                 = " << p.TroyonG << endl;
  cout << " B0^2                   = " << p.b0SquareCentre << endl;
  cout << " li                     = " << p.fourPiInductance << endl;
  cout << " Ip                     = " << p.totalToroidalCurrent << endl;
  cout << p.nIdealInstabilities << " ideal instabilitie(s)\n";
  cout << " ==========================================================\n";
  for(i=0; i < p.nosing; i++){
    cout << i << "  Rat. surf. psi/psi_a = " << p.psis[i] << 
      '/' << p.psia << endl;
    cout << "    safety factor       = " << p.qs[i] << endl;
    cout << "    sqrt(-D_I)          = " << p.mu[i] << endl;
    cout << "    Mercier D_R         = " << p.dr[i] << endl;
    cout << "    <rs d psi/d r>      = " << p.rdpdr[i] << " (estimate)\n";
    cout << "    Lambdas [psi_s norm]= " << p.lambdas[i] << endl;
    for(int j = 0; j < p.nosing; j++){
      cout << i << ' ' << j << " psi_s^(2 mu) Delta' = " << 
	p.re_deltap[i][j] << " + i " << 
	p.im_deltap[i][j] << " +/- " << 
	p.err_deltap[i][j] << endl; 
      cout << i << ' ' << j << " psi_s^(2 mu) Gamma' = " << 
	p.re_gamtap[i][j] << " + i " << 
	p.im_gamtap[i][j] << " +/- " << 
	p.err_gamtap[i][j] << endl; 
    }
    cout << " ----------------------------------------------------------\n";
  }

  return 0;
}
	    

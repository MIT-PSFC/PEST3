int NPlVac=18,*k2kPlV;
int nf_X=2;

double *rPlVd,*zPlVd;

double EZd_vb=-0.08,EZgd_vb=1.,b_X,z_X0,z_X2,z_X2n,t_X,singtX;
double *EZx0sep,*EZx1sep,*EZx2sep,*X2vb,*EZx0cs,*EZx0sn,*EZdx0cs,*EZdx0sn,*X2vbc,*X2vbs;
double EZdx_0x,d2x_0x,EZdx_0o,d2x_0o;

double d_vbE,b_XE,z_X0E,z_X2E;
double dx_0xE,d2x_0xE,dx_0oE,d2x_0oE;
double *x0sepE,*x1sepE,*x2sepE;

double *EZrvb,*EZzvb,*EZrvbcs,*EZrvbsn,*drvbcs,*drvbsn;

double *Dx0,*Dx2,*Ddx2,*Dz2,*Dx2z,*Ddx2z,*Dx0c,*Dx0s,*Dx2c,*Dx2s,
*Ddx2c,*Ddx2s,*Dz2c,*Dz2s,*Dx2zc,*Dx2zs,*Ddx2zc,*Ddx2zs;

double *EZrcs,*EZrsn,*EZd1rcs,*EZd1rsn,*EZd2rcs,*EZd2rsn,*EZz0,*EZd1z0,*EZd2z0;
double *rcE,*rsE,*rcE1,*rsE1,*rcE2,*rsE2,*bE,*bE1,*bE2,R0E,Z0E;

double*Fcs,*Fsn,*dFcs,*dFsn,*EZxinf,*A2m,*EZyinf,*d1yinf,*d2yinf,*EZxgt,*EZygt;
double Eps_tr=1e-4,Eps_noise=10.;

double *EZgper,*EZgpei,*EZdgper,*EZdgpei,*aYer,*aYei,*daYer,*daYei;
double *rCc,*rCs,*rSc,*rSs,*drCc,*drCs,*drSc,*drSs;
double *bCc,*bSc,*dbCc,*dbSc;

double EZga0=0,EZga1=1e-10,EZga2=0.,EZga3=5e-11;
int M0Mm;
double eps_rel=1e-6;
int nPIn1=0,nXIn1=0,nPIn2=0,nXIn2=0;

int ECEqIterFl=0;
double ESreltol=1e-6,ESabstol=1e-12;


#ifndef __MV_solve__
#define __MV_solve__

#include "MV_matrix.h"


/*

inverse, system of equations, etc

Solve(a, b)      x = a^(-1) * b using lapack (b can be MV_Vect)
Inverse(a)

*/

extern MV_ColMat<double> Solve(MV_ColMat<double> &a, MV_ColMat<double> &b);
extern MV_Vector<double> Solve(MV_ColMat<double> &a, MV_Vector<double> &b);
extern MV_ColMat<double> Inverse(MV_ColMat<double> &a); 



#endif /* __MV_solve__ */

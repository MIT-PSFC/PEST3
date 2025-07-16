/*
 * mds__short   -- Return mds descriptor type for INTEGER*2
 * mds__long    -- Return mds descriptor type for INTEGER
 * mds__cstring -- Return mds descriptor type for CHARACTER
 * mds__float   -- Return mds descriptor type for REAL*4
 * mds__double  -- Return mds descriptor type for REAL*8
 *
 *
 * see "mdsdescrip.h" for DTYPE_F or DTYPE_FS float logic
 * see "mdsdescrip.h" for DTYPE_G, DTYPE_D, or DTYPE_FT double logic
 *
 * 03/15/1999 ler
 * 06/22/2000 ler For mds__float  use DTYPE_F which is VMS F_FLOATING
 *                    mds__double use DTYPE_D which is VMS D_FLOATING
 */
#ifdef __NOMDSPLUS
#define DTYPE_L 8
#define DTYPE_T 14
#define DTYPE_F 10
#define DTYPE_D 11
#define DTYPE_W 7
#else
#include <mdsdescrip.h>
#endif
 
#include "fpreproc/f77name.h"
 
int F77NAME(mds__long) ()    { return DTYPE_L;      }
int F77NAME(mds__cstring) () { return DTYPE_T;      }
int F77NAME(mds__float) ()   { return DTYPE_F;      }
int F77NAME(mds__double) ()  { return DTYPE_D;      }
int F77NAME(mds__short) ()   { return DTYPE_W;      }

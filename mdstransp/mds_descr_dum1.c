#include <stdio.h>
/* moved from mdsdummy.c to here by Pletzer */

#ifdef __NOMDSPLUS

int descr()
{
  printf(" __mdsdummy.c:  dummy descr call (MDSplus).\n");
  return(-1);
}
int descr_()
{descr();
return(-1);}

#else

int mds_descr_dum1()
{
  return(-1);
}

#endif

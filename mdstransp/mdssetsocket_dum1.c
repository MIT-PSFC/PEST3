#include <stdio.h>

/* moved from mdsdummy.c to here by Pletzer */

#ifdef __NOMDSPLUS
int MdsSetSocket()
{
  printf(" __mdsdummy.c:  dummy MdsSetSocket call.\n");
  printf("   link with a more recent version of MdsLib.\n");
  return(-1);
}
int mdssetsocket()
{
  MdsSetSocket();
  return(-1);
}
int mdssetsocket_()
{
  MdsSetSocket();
  return(-1);
}

#else
int mdssetsocket_dum1()
{
  return(-1);
}
#endif

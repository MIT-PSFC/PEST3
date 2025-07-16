/* csleep.c */
/* sleep callable from fortran */
 
#include <unistd.h>
 
/* function names link differently depending on OS */
 
#include "fpreproc/f77name.h"
#if (WIN32)
#include <windows.h>
#endif
 
void F77NAME(csleep )(iarg)
     int *iarg;   /* integer form -- ascii character to write */
{
  int j;
 
#if (WIN32)
  Sleep ( *iarg );
#else
  j = sleep ( *iarg );
#endif
 
}

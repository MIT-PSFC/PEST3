/* mdsfconnect.c
**
** Fortran callable wrapper for MdsConnect.
**
** Fortran, depending on the operating system, expects names
** like "mdsconnect" or "mdsconnect_". So this wrapper has a
** Fortran-like name and then calls the C-like named function.
**
** **** WARNING ****
**  SOCKET defined here instead of ipdesc.h
**  Mds prototypes defined here instead of MdsLib.h
** *****************
**
** 04/16/1999 ler
*/
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "fpreproc/f77name.h"
 
#define SOCKET int           /* from MdsLib.h */
 
SOCKET MdsConnect(char *host);
SOCKET MdsCacheConnect(char *host);
 
SOCKET F77NAME(mdsfconnect) (char *host)
{
#ifndef __NOMDSPLUS
  SOCKET id;
  int loop, n, j, wait, k, ic;
  char *p;
  time_t    now;
  struct tm *t;
  char date[20];

  ic = -1 ;
  for (k=strlen(host)-1 ; k>=0 && isspace(host[k]) ; k--) ic=k ;
  if (ic>=0) host[ic]='\0' ;  /* trim host name */

  n=1;
  id=0;
  wait=10;
  if ( p = getenv("MDS_TRANSP_RETRY") ) {
    n=atoi(p);
  }

  for ( loop=1; loop <= n && id <=0; loop++) {
    now = time(NULL);
    id = MdsCacheConnect (host);
    if (id <=0) {
      t = localtime(&now);
      /*      sprintf(date,"%d/%d %d:%d:%d \0",t->tm_mon,t->tm_mday,t->tm_hour,t->tm_min,t->tm_sec); */
      sprintf(date,"    %d:%d:%d \0",t->tm_hour,t->tm_min,t->tm_sec); 
      printf("\n%s ... connection to %s failed\n",date,host);
      if (loop < n) {
	printf("\n... retry in %d sec\n",wait); 
	j = sleep ( wait );
      }
    }
  }

  if (ic>=0) host[ic]=' ' ;  /* restore host name */
  return(id);
#else
  /* if no mdsplus */
  printf("\n...dummy mdsfconnect (no mds library)\n");
  return 0;
#endif
}

/*
  Replacement for MdsConnect(char*) and MdsDisconnect() which caches and reuses sockets.
*/

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>

#include "fpreproc/f77name.h"

#define SOCKET int           /* from MdsLib.h */
#define INVALID_SOCKET -1

static SOCKET* mds_cache_connect_socket  = NULL ;       /* stored sockets */
static char**  mds_cache_connect_server  = NULL ;       /* stored servers in order of requested connection */
static size_t  mds_cache_connect_max     = 0 ;          /* allocated size of cache */
static size_t  mds_cache_connect_num     = 0 ;          /* number of stored sockets in the cache */

static size_t  mds_cache_connect_max_default = 100 ;    /* default value of mds_cache_connect_max, 100 means it was not set */
static int     mds_cache_connect_max_verbose = 0 ;      /* nonzero for verbose messages */

#define VERBOSE mds_cache_connect_max_verbose

SOCKET MdsConnect(char *host);
void   MdsDisconnect();
SOCKET MdsSetSocket(SOCKET *socket) ;

/*
  Set the default cache size before the first mdsCacheConnect.  This overrides the
  MDS_CACHE_CONNECT environment variable.  A cache size <0 will use the absolute value
  and turn on verbose messages.
*/
void MdsCacheConnectDefault(int* pmax) {
  int imax ;

  imax = *pmax ;
  mds_cache_connect_max_verbose = imax<0 ;

  if (imax<0)    imax=-imax ;
  if (imax>=100) imax=0 ;
  
  mds_cache_connect_max_default = (size_t)imax ;  /* final value must be <100 to mark it as set */
}


/*
  Print out the socket cache if writing verbose messages.
*/
void MdsCacheConnectShow() {
  size_t k ;

  if (VERBOSE && mds_cache_connect_socket!=NULL && mds_cache_connect_max>0) {
    if (mds_cache_connect_num==0) {
      printf("%%MdsCacheConnectShow: empty cache\n") ;
    }
    else {
      for (k=0 ; k<mds_cache_connect_num ; k++) {  
	printf("%%MdsCacheConnectShow: [%2d]  socket=%-4d  %s\n",k,
	       mds_cache_connect_socket[k],mds_cache_connect_server[k]) ;
      }
    }
  }
}


/*
  Drop in replacement for MdsConnect(char*) which caches and reuses sockets identified
  by case insensitive server name.  The size of the cache is set by the MDS_CACHE_CONNECT environment
  variable and defaults to mds_cache_connect_max_default if not set.  MDS_CACHE_CONNECT can be 0 
  to prevent caching.  Setting MDS_CACHE_CONNECT<0 will use the absolute value and turn on 
  verbose messages.  abs(MDS_CACHE_CONNECT)>=100 is equivalent to 0. 

  A server name starting with ! will always execute MdsConnect() on the server name without
  the leading !.

  A server name of only ! will clear out the cache and leave mdsplus disconnected.
*/
SOCKET MdsCacheConnect(char* aserver) {
  char*  menv ;        /* value of environment variable */
  int    i ;           /* temp */
  size_t sz,k,kk ;     /* temp */
  int    iconnect ;    /* true when need to do MdsConnect */
  char*  server ;      /* actual server name to use */
  char*  terver ;      /* lower case server name */
  size_t lserver ;     /* length of server,terver name */
  char*  kerver ;      /* loop server name */
  size_t lkerver ;     /* length of kerver name */
  char   c ;           /* temp */
  SOCKET sock ;        /* working socket */

  if (mds_cache_connect_socket==NULL) {
    /* initialize the cache storage */

    mds_cache_connect_max = mds_cache_connect_max_default ;  

    if (mds_cache_connect_max>=100) {
      mds_cache_connect_max = 0 ;

      menv = getenv("MDS_CACHE_CONNECT") ;
      if (menv!=NULL) {
	i = atoi(menv) ;
	mds_cache_connect_max_verbose = i<0 ;

	if (VERBOSE) printf("%%MdsCacheConnect: found MDS_CACHE_CONNECT = %d\n",i) ;

	if (i<0)    i=-i ;
	if (i>=100) i=0 ;
	mds_cache_connect_max = (size_t)i ;
      }
    }
    
    if (VERBOSE) printf("%%MdsCacheConnect: cache size = %d\n",mds_cache_connect_max) ;

    sz = (1+mds_cache_connect_max)*sizeof(SOCKET) ;   /* keep >0 */
    mds_cache_connect_socket = malloc(sz) ;

    sz = (1+mds_cache_connect_max)*sizeof(char*) ;
    mds_cache_connect_server = malloc(sz) ;

    mds_cache_connect_num = 0 ;
  }

  iconnect = mds_cache_connect_max==0 ;
  server   = aserver ;

  if (VERBOSE) printf("%%MdsCacheConnect: server = '%s'\n",server) ;

  while((*server)!=0) {
    if (!isblank(*server)) break ;        /* remove leading blanks */
    server++ ;
  }

  lserver = strlen(server) ;              /* find last non blank */
  for (k=lserver; k>0 ; k--) {
    if (!isblank(server[k-1])) break ;
    lserver-- ;
  }

  if (lserver>0 && server[0]=='!') {      /* first character might be ! to force an mdsconnect */
    iconnect = 1 ;
    server++ ;       /* remove ! */
    lserver-- ;

    if (lserver==0 && mds_cache_connect_max>0) {    /* server was '!' so clear the cache */

      if (VERBOSE) printf("%%MdsCacheConnect: clearing the cache and disconnecting the sockets\n") ;

      for (k=0 ; k<mds_cache_connect_num ; k++) {
	if (VERBOSE) printf("%%MdsCacheConnect: disconnecting from '%s', socket=%d\n",
			    mds_cache_connect_server[k],mds_cache_connect_socket[k]) ;
	free(mds_cache_connect_server[k]) ;
	MdsSetSocket(&(mds_cache_connect_socket[k])) ;
	MdsDisconnect() ;
      }
      mds_cache_connect_num = 0 ;

      sock = INVALID_SOCKET ;
      MdsSetSocket(&sock) ;       /* current socket in mdsplus was just disconnected */

      MdsCacheConnectShow() ;
      return INVALID_SOCKET ;
    }
  }

  terver = malloc((lserver+1)*sizeof(char)) ;
  strncpy(terver,server,lserver) ;
  terver[lserver]=0 ;
  
  for (k=0 ; k<lserver ; k++) {           /* lower case server name */
    terver[k] = (char)tolower(terver[k]) ;
  }

  if (VERBOSE) printf("%%MdsCacheConnect: lower case server = '%s'\n",terver) ;

  if (mds_cache_connect_max>0) {                     /* look for server in cache */
    for (k=0 ; k<mds_cache_connect_num ; k++) {

      kerver  = mds_cache_connect_server[k] ;
      lkerver = strlen(kerver) ;
      
      if (lserver==lkerver) {
	if (lserver==0 || strncmp(kerver,terver,lserver)==0) { 

	  /* found server called kerver at index k in the cache */

	  sock = mds_cache_connect_socket[k] ;    

	  if (!iconnect) {
	    /* reuse socket from cache */
	    free(terver) ;                           /* don't need this anymore */
	    terver = NULL ;
	    
	    for (kk=k ; kk>0 ; kk--) {            
	      mds_cache_connect_socket[kk] = mds_cache_connect_socket[kk-1] ;  /* move sockets up */
	      mds_cache_connect_server[kk] = mds_cache_connect_server[kk-1] ;
	    }
	    mds_cache_connect_socket[0] = sock ;                               /* put found socket at head of cache */
	    mds_cache_connect_server[0] = kerver ;
	    
	    if (VERBOSE) printf("%%MdsCacheConnect: found server '%s', reusing socket=%d\n",kerver,sock) ;
	    
	    MdsSetSocket(&sock) ; 
	    
	    MdsCacheConnectShow() ;
	    return sock ;
	  }
	  else {
	    /* will be reconnecting, disconnect and remove this socket from the cache */
	    if (VERBOSE) printf("%%MdsCacheConnect: disconnecting and removing cached server '%s', socket=%d\n",
				kerver,sock) ;

	    free(mds_cache_connect_server[k]) ;

	    MdsSetSocket(&sock) ;
	    MdsDisconnect() ;

	    mds_cache_connect_num-- ;
	    for (kk=k ; kk<mds_cache_connect_num ; kk++) {                    /* shrink the list */
	      mds_cache_connect_socket[kk] = mds_cache_connect_socket[kk+1] ;
	      mds_cache_connect_server[kk] = mds_cache_connect_server[kk+1] ;
	    }	    
	  }
	  break ;  /* end found server, break out of for loop */

	}
      }

    } /* end for */
  }

  /* still not connected */
  if (VERBOSE) printf("%%MdsCacheConnect: calling MdsConnect on server '%s'\n",server) ;

  if (mds_cache_connect_max>0) {
    sock = INVALID_SOCKET ;               /* prevent disconnection of existing cached socket */
    MdsSetSocket(&sock) ;
  }

  sock = MdsConnect(server) ;

  if (sock!=INVALID_SOCKET && mds_cache_connect_max>0) {         /* store new socket */

    for (k=0 ; k<mds_cache_connect_num ; k++) {    /* if new socket matches an existing one, replace */
      if (mds_cache_connect_socket[k]==sock) {
	if (VERBOSE) printf("!MdsCacheConnect: duplicate socket, removing cached server '%s', socket=%d\n",
			    mds_cache_connect_server[k],sock) ;
	free(mds_cache_connect_server[k]) ;        /* keep the most recent name */

	for (kk=k ; kk>0 ; kk--) {                 /* move this socket to the head of the list */
	  mds_cache_connect_socket[kk] = mds_cache_connect_socket[kk-1] ;
	  mds_cache_connect_server[kk] = mds_cache_connect_server[kk-1] ;
	}
	mds_cache_connect_socket[0] = sock ;
	mds_cache_connect_server[0] = terver ;

	terver = NULL ;

	if (VERBOSE) printf("%%MdsCacheConnect: adding server to cache '%s', socket=%d\n",
			    mds_cache_connect_server[0],mds_cache_connect_socket[0]);

	MdsCacheConnectShow() ;
	return sock ;
      }
    }

    if (mds_cache_connect_num==mds_cache_connect_max) {     /* cache is full, remove last one */
      k = mds_cache_connect_num-1 ;

      if (VERBOSE) printf("%%MdsCacheConnect: full cache, removing cached server '%s', disconnecting socket=%d\n",
			  mds_cache_connect_server[k],mds_cache_connect_socket[k]) ;

      free(mds_cache_connect_server[k]) ;
      MdsSetSocket(&(mds_cache_connect_socket[k])) ;        /* set so it can be disconnected */
      MdsDisconnect() ;
      MdsSetSocket(&sock) ;                                 /* restore current connected socket */
      mds_cache_connect_num = k ;
    }

    /* put new socket at beginning of list */
    if (VERBOSE) printf("%%MdsCacheConnect: adding server to cache '%s', socket=%d\n",terver,sock);

    for (k=mds_cache_connect_num ; k>0 ; k--) {
      mds_cache_connect_socket[k] = mds_cache_connect_socket[k-1] ;
      mds_cache_connect_server[k] = mds_cache_connect_server[k-1] ;
    }
    mds_cache_connect_socket[0] = sock ;
    mds_cache_connect_server[0] = terver ;
    mds_cache_connect_num++ ;

    terver = NULL ;
  }
  else if (sock==INVALID_SOCKET) {
    printf("?MdsCacheConnect: failed to connect to '%s'\n",server);
  }
  else if (VERBOSE) {
    printf("%%MdsCacheConnect: connected but server was NOT cached '%s', socket=%d\n",terver,sock);
  }

  if (terver!=NULL) free(terver) ;

  MdsCacheConnectShow() ;
  return sock ;
}


/*
  Drop in replacement for MdsDisconnect().  This disconnects the current socket along with
  removing the socket from the cache.
*/
void MdsCacheDisconnect() {
  SOCKET sock, isock ;   
  size_t k,kk ;

  if (mds_cache_connect_socket!=NULL && mds_cache_connect_max>0 && mds_cache_connect_num>0) {
    isock = INVALID_SOCKET ;
    sock  = MdsSetSocket(&isock) ;

    for (k=0 ; k<mds_cache_connect_num ; k++) {    /* look for the socket in the cache */
      if (mds_cache_connect_socket[k]==sock) {
	if (VERBOSE) printf("%%MdsCacheDisconnect: for current socket, removing cached server '%s', socket=%d\n",
			    mds_cache_connect_server[k],sock) ;
	free(mds_cache_connect_server[k]) ;              

	mds_cache_connect_num-- ;
	for (kk=k ; kk<mds_cache_connect_num ; kk++) {                    /* shrink the list */
	  mds_cache_connect_socket[kk] = mds_cache_connect_socket[kk+1] ;
	  mds_cache_connect_server[kk] = mds_cache_connect_server[kk+1] ;
	}
	break ;
      }
    }

    MdsSetSocket(&sock) ;  /* restore socket prior to disconnect */
  }
  if (VERBOSE) printf("%%MdsCacheDisconnect: disconnecting current socket\n") ;

  MdsDisconnect() ;
  MdsCacheConnectShow() ;
}


/*
  Set the mdsplus socket to INVALID_SOCKET but do not disconnect unless not caching sockets.
*/
void MdsCacheUnconnect() {
  SOCKET isock ;   

  if (mds_cache_connect_socket!=NULL && mds_cache_connect_max>0 && mds_cache_connect_num>0) {
    if (VERBOSE) printf("%%MdsCacheUnconnect: unconnecting current socket\n") ;

    isock = INVALID_SOCKET ;
    MdsSetSocket(&isock) ;

    MdsCacheConnectShow() ;
  }
  else {
    if (VERBOSE) printf("%%MdsCacheUnconnect: disconnecting current socket\n") ;

    MdsDisconnect() ;
  }
}


SOCKET F77NAME(c_mdscacheconnect)(char* server) {
  return MdsCacheConnect(server) ;
}

void F77NAME(mdscachedisconnect)() { 
  MdsCacheDisconnect() ;
}

void F77NAME(mdscacheunconnect)() { 
  MdsCacheUnconnect() ;
}

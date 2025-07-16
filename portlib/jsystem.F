      integer function jsystem(str)
c
c  Unix "system" call
c
c-Mods
c   20Jan2010    jim.conboy@ccfe.ac.uk
c                Log error returns from call
c
      use logmod
c
#include "fpreproc/byte_declare.h"
c
c
#if __NAGWARE
      use f90_unix_proc       ! Linux NagWare
 
      character*(*) str
      integer status
 
      call system(str,status)
      jsystem = status
 
#elif __CRAY
      character*(*) str
      integer l
      character*(255) line
      integer CSYSTEM
 
      l = len(str)
      if(l.gt.255) then
         write(6,1001) str(1:l)
 1001    format(' ?jsystem:  command too long:  '/1x,a)
         jsystem = 2
         return
      end if
      line = str//char(0)
      jsystem = CSYSTEM(line,l)
 
#elif __UNIX
      character*(*) str
      BYTE_DECLARE line(256)
      integer l,str_length
      integer fsystem
 
      l = str_length(str)
      if(l.gt.255) then
         write(6,1001) str(1:l)
 1001    format(' ?jsystem:  command too long:  '/1x,a)
         jsystem = 2
         return
      end if
      !!DMC call infolog('jsystem -I- :'//str(1:l) )
      call cstring(str(1:l),line,'2C')
      jsystem = fsystem(line)
      if( jsystem /= 0 )     then
         call setlogbuffered(.false.)
         call errorlog( 'jsystem -E- :'//trim(str)//': ' , jsystem )
      else
         call infolog( 'jsystem -I- :'//trim(str) )
         continue
      endif
c
#elif __VMS
      character*(*) str     ! VMS stub
      integer str_length
 
      jsystem = lib$spawn(str(1:str_length(str)))
#endif
      return
      end

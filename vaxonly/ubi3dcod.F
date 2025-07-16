      subroutine ubi3dcod(ilbuf,istart,ivalue)
C
C  decode a 3 byte sequence as an integer, hi order bits first
C
      implicit NONE
C
      integer ivalue
 
#include "fpreproc/byte_declare.h"
      BYTE_DECLARE ilbuf(*)
 
      integer istart,ibyte(3)
C
      integer ict,is
C
      ict=0
      do is=istart,istart+2
         ict=ict+1
#if __CRAY || __SX
         if(ichar(ilbuf(is)).lt.0) then
            ibyte(ict)=256+ichar(ilbuf(is))
         else
            ibyte(ict)=ichar(ilbuf(is))
         endif
#else
         if(ilbuf(is).lt.0) then
            ibyte(ict)=256+ilbuf(is)
         else
            ibyte(ict)=ilbuf(is)
         endif
#endif
      enddo
C
      ivalue=256*(256*ibyte(1)+ibyte(2))+ibyte(3)
C
      return
      end

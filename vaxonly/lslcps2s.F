C  dmc-- force byte(4)/integer equivalence through subroutine
C        call interface; NagWare won't let me use EQUIVALENCE
C
#include "fpreproc/byte_declare.h"
C--------------
      subroutine copyint4(intb4,intout)
      implicit NONE
      integer*4 intb4,intout
C
      intout=intb4
      return
      end
C--------------
      subroutine byte4out(byte4,ii,byteo)
C
      implicit NONE
      integer ii
      BYTE_DECLARE byte4(4),byteo
C
      byteo=byte4(ii)
      return
      end

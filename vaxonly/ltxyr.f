C------------------------------------------------------------------
C
C RGA: LTX year from Craig Jacobson email 14Nov2013  
C

C
      SUBROUTINE LTXYR(NSHOT,ZYEAR)
C
      implicit NONE
C
      integer NSHOT
      CHARACTER*(*) ZYEAR
C
      integer ilen
C
C
      zyear=' '
      ilen=len(zyear)
C
      if(ilen.gt.3) then
         IF(NSHOT.GE. 108000) then
            ZYEAR(ilen-3:ilen)='2013'
         else
            ZYEAR(ilen-3:ilen)='2012'
         endif
      else
         IF(NSHOT.GE.108000) then
            ZYEAR(ilen-1:ilen)='13'
         else
            ZYEAR(ilen-1:ilen)='12'
         endif
      endif
C
      RETURN
      END

 

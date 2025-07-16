C------------------------------------------------------------------
C
C  ARIS vsn adapted from TFTRYR.FOR
C  06/08/12 C. Ludescher-Furth

C
      SUBROUTINE ARISYR(NSHOT,ZYEAR)
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
         IF(NSHOT.GE. 10000) then
            ZYEAR(ilen-3:ilen)='2012'
         else
            ZYEAR(ilen-3:ilen)='2011'
         endif
      else
         IF(NSHOT.GE.10000) then
            ZYEAR(ilen-1:ilen)='12'
         else
            ZYEAR(ilen-1:ilen)='11'
         endif
      endif
C
      RETURN
      END

 

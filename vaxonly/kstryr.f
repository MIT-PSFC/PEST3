C------------------------------------------------------------------
C  KSTRYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
C
C  KSTR vsn adapted from TFTRYR.FOR
C  04/15/13 C. Ludescher-Furth
C           shot info from Robert Budny

C
      SUBROUTINE KSTRYR(NSHOT,ZYEAR)
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
         IF(NSHOT.GE.8356) then
            ZYEAR(ilen-3:ilen)='2013'
         else IF(NSHOT.GE.6471) then
            ZYEAR(ilen-3:ilen)='2012'
         else IF(NSHOT.GE.4469) then
            ZYEAR(ilen-3:ilen)='2011'
         else IF(NSHOT.GE.2372) then
            ZYEAR(ilen-3:ilen)='2010'
         else
            ZYEAR(ilen-3:ilen)='2009'
         endif
      else
         IF(NSHOT.GE.8356) then
            ZYEAR(ilen-1:ilen)='13'
         else IF(NSHOT.GE.6471) then
            ZYEAR(ilen-1:ilen)='12'
         ELSE IF(NSHOT.GE.4469) then
            ZYEAR(ilen-1:ilen)='11'
         else IF(NSHOT.GE.6471) then
            ZYEAR(ilen-1:ilen)='10'
         else
            ZYEAR(ilen-1:ilen)='09'
         endif
      endif
C
      RETURN
      END

 

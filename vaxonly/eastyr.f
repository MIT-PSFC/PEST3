C------------------------------------------------------------------
C  EASTYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
C
C  EAST vsn adapted from TFTRYR.FOR
C  10/12/04 C. Ludescher-Furth
C           shot info from Siye Ding

C
      SUBROUTINE EASTYR(NSHOT,ZYEAR)
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
         IF(NSHOT.GE.44327) then
            ZYEAR(ilen-3:ilen)='2013'
         else IF(NSHOT.GE.36813) then
            ZYEAR(ilen-3:ilen)='2012'
         else IF(NSHOT.GE.14705) then
            ZYEAR(ilen-3:ilen)='2010'
         else IF(NSHOT.GE.9989) then
            ZYEAR(ilen-3:ilen)='2009'
         else IF(NSHOT.GE.5672) then
            ZYEAR(ilen-3:ilen)='2008'
         else IF(NSHOT.GE.3080) then
            ZYEAR(ilen-3:ilen)='2007'
         else IF(NSHOT.GE.1138) then
            ZYEAR(ilen-3:ilen)='2006'
         else
            ZYEAR(ilen-3:ilen)='2005'
         endif
      else
         IF(NSHOT.GE.44327) then
            ZYEAR(ilen-1:ilen)='13'
         else IF(NSHOT.GE.36813) then
            ZYEAR(ilen-1:ilen)='12'
         ELSE IF(NSHOT.GE.14705) then
            ZYEAR(ilen-1:ilen)='10'
         else IF(NSHOT.GE.9989) then
            ZYEAR(ilen-1:ilen)='09'
         else IF(NSHOT.GE.5672) then
            ZYEAR(ilen-1:ilen)='08'
         else IF(NSHOT.GE.3080) then
            ZYEAR(ilen-1:ilen)='07'
         else IF(NSHOT.GE.1138) then
            ZYEAR(ilen-1:ilen)='06'
         else
            ZYEAR(ilen-1:ilen)='05'
         endif
      endif
C
      RETURN
      END

 

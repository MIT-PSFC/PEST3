C******************** START FILE NSTXYR.FOR ; GROUP TFTRG1 *************
C------------------------------------------------------------------
C& NSTXYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
C
C  NSTX vsn adapted from TFTRG1:TFTRYR.FOR
C  I don't have all the calender year boundary shot numbers!
C
      SUBROUTINE NSTXYR(NSHOT,ZYEAR)
C
      implicit NONE
C
      integer NSHOT
      CHARACTER*(*) ZYEAR
C
      integer ilen,ilen1,ilen3
C
C  A CAREFUL EXAMINATION OF THE RECORD YIELDS...
C
      zyear=' '
      ilen=len(zyear)
C
      ilen1 = max(1,ilen-1) ! RGA, quick fix for out of range error
      ilen3 = max(1,ilen-3)
      if(ilen.gt.3) then
         ZYEAR(ilen3:ilen)='2000'
         IF(NSHOT.GE.104515) ZYEAR(ilen3:ilen)='2001'
         IF(NSHOT.GE.106641) ZYEAR(ilen3:ilen)='2002'
         IF(NSHOT.GE.109080) ZYEAR(ilen3:ilen)='2003'
         IF(NSHOT.GE.110945) ZYEAR(ilen3:ilen)='2004'
         IF(NSHOT.GE.115000) ZYEAR(ilen3:ilen)='2005'
         IF(NSHOT.GE.119000) ZYEAR(ilen3:ilen)='2006'
         IF(NSHOT.GE.122001) ZYEAR(ilen3:ilen)='2007'
         IF(NSHOT.GE.125858) ZYEAR(ilen3:ilen)='2008'
         IF(NSHOT.GE.130890) ZYEAR(ilen3:ilen)='2009'
         IF(NSHOT.GE.136452) ZYEAR(ilen3:ilen)='2010'
         IF(NSHOT.GE.142633) ZYEAR(ilen3:ilen)='2013'
      else
         ZYEAR(ilen1:ilen)='00'
         IF(NSHOT.GE.104515) ZYEAR(ilen1:ilen)='01'
         IF(NSHOT.GE.106641) ZYEAR(ilen1:ilen)='02'
         IF(NSHOT.GE.109080) ZYEAR(ilen1:ilen)='03'
         IF(NSHOT.GE.110945) ZYEAR(ilen1:ilen)='04'
         IF(NSHOT.GE.115000) ZYEAR(ilen1:ilen)='05'
         if(nshot.ge.119000) ZYEAR(ilen1:ilen)='06'
         IF(NSHOT.GE.122001) ZYEAR(ilen1:ilen)='07'
         IF(NSHOT.GE.125858) ZYEAR(ilen1:ilen)='08'
         IF(NSHOT.GE.130890) ZYEAR(ilen1:ilen)='09'
         IF(NSHOT.GE.136452) ZYEAR(ilen1:ilen)='10'
         IF(NSHOT.GE.142633) ZYEAR(ilen1:ilen)='13'
      endif
C
      RETURN
      END
C******************** END FILE NSTXYR.FOR ; GROUP TFTRG1 ***************
 

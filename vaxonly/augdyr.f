C------------------------------------------------------------------
C  AUGDYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
C
C  AUGD vsn adapted from TFTRYR.FOR
C  10/12/04 C. Ludescher-Furth
C           shot info from G. Tardini (IPP)

C
      SUBROUTINE AUGDYR(NSHOT,ZYEAR)
C
      implicit NONE
C
      integer NSHOT
      CHARACTER*(*) ZYEAR
C
      integer ilen
C
C  A CAREFUL EXAMINATION OF THE RECORD YIELDS...
C
      zyear=' '
      ilen=len(zyear)
C
      if(ilen.gt.3) then
         ZYEAR(ilen-3:ilen)='1999'
         IF(NSHOT.GE.12980) ZYEAR(ilen-3:ilen)='2000'
         IF(NSHOT.GE.14049) ZYEAR(ilen-3:ilen)='2001'
         IF(NSHOT.GE.15011) ZYEAR(ilen-3:ilen)='2002'
         IF(NSHOT.GE.16624) ZYEAR(ilen-3:ilen)='2003'
         IF(NSHOT.GE.18244) ZYEAR(ilen-3:ilen)='2004'
         IF(NSHOT.GE.19630) ZYEAR(ilen-3:ilen)='2005'
         IF(NSHOT.GE.20713) ZYEAR(ilen-3:ilen)='2006'
         IF(NSHOT.GE.21482) ZYEAR(ilen-3:ilen)='2007'
         IF(NSHOT.GE.22586) ZYEAR(ilen-3:ilen)='2008'
         IF(NSHOT.GE.24190) ZYEAR(ilen-3:ilen)='2009'
         IF(NSHOT.GE.25890) ZYEAR(ilen-3:ilen)='2010'
         IF(NSHOT.GE.26085) ZYEAR(ilen-3:ilen)='2011'
         IF(NSHOT.GE.27405) ZYEAR(ilen-3:ilen)='2012'
         IF(NSHOT.GE.29145) ZYEAR(ilen-3:ilen)='2013'
      else
         ZYEAR(ilen-1:ilen)='99'
         IF(NSHOT.GE.12980) ZYEAR(ilen-1:ilen)='00'
         IF(NSHOT.GE.14049) ZYEAR(ilen-1:ilen)='01'
         IF(NSHOT.GE.15011) ZYEAR(ilen-1:ilen)='02'
         IF(NSHOT.GE.16624) ZYEAR(ilen-1:ilen)='03'
         IF(NSHOT.GE.18244) ZYEAR(ilen-1:ilen)='04'
         IF(NSHOT.GE.19630) ZYEAR(ilen-1:ilen)='05'
         IF(NSHOT.GE.20713) ZYEAR(ilen-1:ilen)='06'
         IF(NSHOT.GE.21482) ZYEAR(ilen-1:ilen)='07'
         IF(NSHOT.GE.22586) ZYEAR(ilen-1:ilen)='08'
         IF(NSHOT.GE.24190) ZYEAR(ilen-1:ilen)='09'
         IF(NSHOT.GE.25890) ZYEAR(ilen-1:ilen)='10'
         IF(NSHOT.GE.26085) ZYEAR(ilen-1:ilen)='11'
         IF(NSHOT.GE.27405) ZYEAR(ilen-1:ilen)='12'
         IF(NSHOT.GE.29145) ZYEAR(ilen-1:ilen)='13'
      endif
C
      RETURN
      END

 

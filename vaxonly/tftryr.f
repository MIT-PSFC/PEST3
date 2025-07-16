C******************** START FILE TFTRYR.FOR ; GROUP TFTRG1 *************
C------------------------------------------------------------------
C& TFTRYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
C
      SUBROUTINE TFTRYR(NSHOT,ZYEAR)
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
         ZYEAR(ilen-3:ilen)='1983'
         IF(NSHOT.GE.8839) ZYEAR(ilen-3:ilen)='1984'
         IF(NSHOT.GE.12000) ZYEAR(ilen-3:ilen)='1985'
         IF(NSHOT.GE.17220) ZYEAR(ilen-3:ilen)='1986'
         IF(NSHOT.GE.26699) ZYEAR(ilen-3:ilen)='1987'
         IF(NSHOT.GE.31951) ZYEAR(ilen-3:ilen)='1988'
         IF(NSHOT.GE.39117) ZYEAR(ilen-3:ilen)='1989'
         IF(NSHOT.GE.45121) ZYEAR(ilen-3:ilen)='1990'
         IF(NSHOT.GE.56486) ZYEAR(ilen-3:ilen)='1991'
         IF(NSHOT.GE.60108) ZYEAR(ilen-3:ilen)='1992'
         IF(NSHOT.GE.69346) ZYEAR(ilen-3:ilen)='1993'
         IF(NSHOT.GE.73489) ZYEAR(ilen-3:ilen)='1994'
         IF(NSHOT.GE.81729) ZYEAR(ilen-3:ilen)='1995'
         IF(NSHOT.GE.90000) ZYEAR(ilen-3:ilen)='1996'
         if(nshot.ge.101951) zyear(ilen-3:ilen)='1997'
      else
         ZYEAR(ilen-1:ilen)='83'
         IF(NSHOT.GE.8839) ZYEAR(ilen-1:ilen)='84'
         IF(NSHOT.GE.12000) ZYEAR(ilen-1:ilen)='85'
         IF(NSHOT.GE.17220) ZYEAR(ilen-1:ilen)='86'
         IF(NSHOT.GE.26699) ZYEAR(ilen-1:ilen)='87'
         IF(NSHOT.GE.31951) ZYEAR(ilen-1:ilen)='88'
         IF(NSHOT.GE.39117) ZYEAR(ilen-1:ilen)='89'
         IF(NSHOT.GE.45121) ZYEAR(ilen-1:ilen)='90'
         IF(NSHOT.GE.56486) ZYEAR(ilen-1:ilen)='91'
         IF(NSHOT.GE.60108) ZYEAR(ilen-1:ilen)='92'
         IF(NSHOT.GE.69346) ZYEAR(ilen-1:ilen)='93'
         IF(NSHOT.GE.73489) ZYEAR(ilen-1:ilen)='94'
         IF(NSHOT.GE.81729) ZYEAR(ilen-1:ilen)='95'
         IF(NSHOT.GE.90000) ZYEAR(ilen-1:ilen)='96'
         if(nshot.ge.101951) zyear(ilen-1:ilen)='97'
      endif
C
      RETURN
      END
C******************** END FILE TFTRYR.FOR ; GROUP TFTRG1 ***************

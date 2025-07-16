C******************** START FILE D3DYR.FOR ; GROUP TFTRG1 *************
C------------------------------------------------------------------
C& D3DYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
C
C  D3D vsn adapted from TFTRG1:TFTRYR.FOR
C  I don't have all the calender year boundary shot numbers!
C
      SUBROUTINE D3DYR(NSHOT,ZYEAR)
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
         ZYEAR(ilen3:ilen)='1983'
         IF(NSHOT.GE.40560) ZYEAR(ilen3:ilen)='1984'
         IF(NSHOT.GE.50531) ZYEAR(ilen3:ilen)='1985'
         IF(NSHOT.GE.50690) ZYEAR(ilen3:ilen)='1986'
         IF(NSHOT.GE.53583) ZYEAR(ilen3:ilen)='1987'
         IF(NSHOT.GE.57504) ZYEAR(ilen3:ilen)='1988'
         IF(NSHOT.GE.62639) ZYEAR(ilen3:ilen)='1989'
         IF(NSHOT.GE.68163) ZYEAR(ilen3:ilen)='1990'
         IF(NSHOT.GE.70801) ZYEAR(ilen3:ilen)='1991'
         IF(NSHOT.GE.74182) ZYEAR(ilen3:ilen)='1992'
         IF(NSHOT.GE.76295) ZYEAR(ilen3:ilen)='1993'
         IF(NSHOT.GE.80293) ZYEAR(ilen3:ilen)='1994'
         IF(NSHOT.GE.84468) ZYEAR(ilen3:ilen)='1995'
         IF(NSHOT.GE.88061) ZYEAR(ilen3:ilen)='1996'
         IF(NSHOT.GE.90895) ZYEAR(ilen3:ilen)='1997'
         IF(NSHOT.GE.94458) ZYEAR(ilen3:ilen)='1998'
         IF(NSHOT.GE.97674) ZYEAR(ilen3:ilen)='1999'
         IF(NSHOT.GE.100513) ZYEAR(ilen3:ilen)='2000'
         IF(NSHOT.GE.104793) ZYEAR(ilen3:ilen)='2001'
         IF(NSHOT.GE.108779) ZYEAR(ilen3:ilen)='2002'
         IF(NSHOT.GE.111982) ZYEAR(ilen3:ilen)='2003'
         IF(NSHOT.GE.116638) ZYEAR(ilen3:ilen)='2004'
         IF(NSHOT.GE.121238) ZYEAR(ilen3:ilen)='2005'
         if(nshot.ge.124000) ZYEAR(ilen3:ilen)='2006'
         IF(NSHOT.GE.127146) ZYEAR(ilen3:ilen)='2007'
         IF(NSHOT.GE.131027) ZYEAR(ilen3:ilen)='2008'
         IF(NSHOT.GE.134884) ZYEAR(ilen3:ilen)='2009'
         IF(NSHOT.GE.140817) ZYEAR(ilen3:ilen)='2010'
         IF(NSHOT.GE.143160) ZYEAR(ilen3:ilen)='2011'
         IF(NSHOT.GE.147780) ZYEAR(ilen3:ilen)='2012'
         IF(NSHOT.GE.151236) ZYEAR(ilen3:ilen)='2013'
      else
         ZYEAR(ilen1:ilen)='83'
         IF(NSHOT.GE.40560) ZYEAR(ilen1:ilen)='84'
         IF(NSHOT.GE.50531) ZYEAR(ilen1:ilen)='85'
         IF(NSHOT.GE.50690) ZYEAR(ilen1:ilen)='86'
         IF(NSHOT.GE.53583) ZYEAR(ilen1:ilen)='87'
         IF(NSHOT.GE.57504) ZYEAR(ilen1:ilen)='88'
         IF(NSHOT.GE.62639) ZYEAR(ilen1:ilen)='89'
         IF(NSHOT.GE.68163) ZYEAR(ilen1:ilen)='90'
         IF(NSHOT.GE.70801) ZYEAR(ilen1:ilen)='91'
         IF(NSHOT.GE.74182) ZYEAR(ilen1:ilen)='92'
         IF(NSHOT.GE.76295) ZYEAR(ilen1:ilen)='93'
         IF(NSHOT.GE.80293) ZYEAR(ilen1:ilen)='94'
         IF(NSHOT.GE.84468) ZYEAR(ilen1:ilen)='95'
         IF(NSHOT.GE.88061) ZYEAR(ilen1:ilen)='96'
         IF(NSHOT.GE.90895) ZYEAR(ilen1:ilen)='97'
         IF(NSHOT.GE.94458) ZYEAR(ilen1:ilen)='98'
         IF(NSHOT.GE.97674) ZYEAR(ilen1:ilen)='99'
         IF(NSHOT.GE.100513) ZYEAR(ilen1:ilen)='00'
         IF(NSHOT.GE.104793) ZYEAR(ilen1:ilen)='01'
         IF(NSHOT.GE.108779) ZYEAR(ilen1:ilen)='02'
         IF(NSHOT.GE.111982) ZYEAR(ilen1:ilen)='03'
         IF(NSHOT.GE.116638) ZYEAR(ilen1:ilen)='04'
         IF(NSHOT.GE.121238) ZYEAR(ilen1:ilen)='05'
         if(nshot.ge.124000) ZYEAR(ilen1:ilen)='06'
         IF(NSHOT.GE.127146) ZYEAR(ilen1:ilen)='07'
         IF(NSHOT.GE.131027) ZYEAR(ilen1:ilen)='08'
         IF(NSHOT.GE.134884) ZYEAR(ilen1:ilen)='09'
         IF(NSHOT.GE.140817) ZYEAR(ilen1:ilen)='10'
         IF(NSHOT.GE.143160) ZYEAR(ilen1:ilen)='11'
         IF(NSHOT.GE.147780) ZYEAR(ilen1:ilen)='12'
         IF(NSHOT.GE.151236) ZYEAR(ilen1:ilen)='13'
      endif
C
      RETURN
      END
C******************** END FILE D3DYR.FOR ; GROUP TFTRG1 ***************
 

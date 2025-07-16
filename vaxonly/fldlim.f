C******************** START FILE FLDLIM.FOR ; GROUP URDUTS *************
C------------------------------
C&   FLDLIM OLD UREAD - FIELD LIMITS, UNPACKED CHAR BUFFER
C
C  FIND BEGINNING AND END OF FIRST STRING OF NON-BLANKS IN
C  DBUFF (1 CHAR/WORD)
C
      SUBROUTINE FLDLIM(DBUFF,NID1,NID2,NUM)
C
      implicit none
C
      integer imaxc
      PARAMETER (IMAXC=512)
C
      CHARACTER*1 DBUFF(IMAXC),MBLANK
      integer nid1,nid2,num
C
      integer i
C
      DATA MBLANK/' '/
C
C-------------------------------------------
C
      DO 10 I=1,IMAXC
        IF(DBUFF(I).NE.MBLANK) GO TO 20
 10   CONTINUE
      NUM=0
      NID1=0
      NID2=0
      RETURN
 20   CONTINUE
      NID1=I
      DO 30 I=NID1,IMAXC
        IF(DBUFF(I).EQ.MBLANK) GO TO 40
 30   CONTINUE
      I=IMAXC+1
 40   CONTINUE
      NID2=I-1
      NUM=NID2-NID1+1
      RETURN
      	END
C******************** END FILE FLDLIM.FOR ; GROUP URDUTS ***************

C	****************** START FILE  TRMFLNAM.FOR  ;  GROUP TRREAD_LIB ******
 
 
      SUBROUTINE TRMFLNAM ( DISK, DIR, RUNID, SUFFIX, FILENAME)
C
C  unix-ification -- dmc 6 Jul 1992
c       dmc-- also allow var. length DISK and DIR arguments
C
C       SUBROUTINE TRMFLNAM MAKES A FILE NAME, GIVEN THE DISK, DIRECTORY, THE
C	    RUNID AND THE FILE SUFFIX. "RUNDATA" IS THE ASSUMED LOGICAL
C	    NAME OF THE DISK IF "DISK" = " "
 
C	.LAST CHANGED:
c       aug 1996  dmc  RUNID -> 10 chars
C	04/09/90  TBT  CHANGED RUNID TO C*8 FROM C*6. DEFINED RUNID9 AS C*9 TO
C                      HANDLE AN 8 CHARACTER RUNID.
C		       DEFINED MAXCHARS=64 AS THE MAXIMUM # OF CHARACTERS IN
C                      THE NEW FILENAME.
 
C	ARGUMENTS:
C	    INPUT:  DISK    C*64 DISK     , E.G. RUNDATA (= DEFAULT IF " ")
C		    DIR	    C*64 DIRECTORY, E.G. TRANSP.TFTR.86
C		    RUNID   C*8  RUNID    , E.G. 4321  OR 10000A01
C		    SUFFIX  C*6	 SUFFIX   , E.G. TF.PLN
C
C	    OUTPUT: FILENAME C*64 FILE NAME, E.G.
C		             RUNDATA:[TRANSP.TFTR.86]4321TF.PLN
C
      IMPLICIT NONE
 
      INTEGER    MAXCHARS
      PARAMETER (MAXCHARS=64)
	
      CHARACTER*(*)    DISK
      CHARACTER*(*)	DIR
      CHARACTER*10	RUNID
      CHARACTER*11	RUNID9	    ! ALWAYS HAS RIGHTMOST CHAR BLANK
      CHARACTER*6	SUFFIX
      CHARACTER*64	FILENAME
 
      INTEGER		IPT,      I
      INTEGER		LUOUT
C
      data LUOUT/6/
C	--------------------------------------------------------------
 
      RUNID9 = RUNID              ! MAKE SURE THERE IS ONE BLANK IN VARIABLE
 
      IPT=0
 
      FILENAME=' '
 
#ifdef __VMS
      IF (DISK(1:1) .EQ. ' ')  THEN
! ASSUME DEFAULT
         FILENAME = 'RUNDATA:['
         IPT = 9
      ELSE
          I = INDEX(DISK, ':') - 1
          IF(I.LE.0) I = INDEX(DISK, ' ') - 1
          if(i.le.0) i = len(disk)
          FILENAME = DISK(1:I) // ':['
          IPT = I+2
      END IF
 
      IF(DIR(1:1).EQ.'[') IPT=IPT-1
      IF(DIR(1:1).EQ.'<') IPT=IPT-1
#endif
 
      I = INDEX(DIR, ' ') - 1
      IF (I .LE. 0)  i = len(dir)
      IF (IPT+I+1.GT. MAXCHARS)  GO TO 999
      FILENAME(IPT+1:IPT+I) = DIR(1:I)
      IPT = IPT+I
C
#ifdef __VMS
      IF((FILENAME(IPT:IPT).NE.']').AND.
     >     (FILENAME(IPT:IPT).NE.'>')) THEN
        FILENAME(IPT+1:IPT+1) = ']'
        IPT = IPT+1
      ENDIF
#endif
C
#ifdef __UNIX
      IF(IPT.GT.0) THEN
        IF(FILENAME(IPT:IPT).NE.'/') THEN
          IPT=IPT+1
          FILENAME(IPT:IPT)='/'
        ENDIF
      ENDIF
#endif
C
 
      I = INDEX(RUNID9,' ') - 1
      IF (I .LE. 0) GO TO 300
      IF (IPT+I+6 .GT. MAXCHARS)  GO TO 999
      FILENAME(IPT+1:IPT+I) = RUNID(1:I)
      IPT = IPT+I
 
  300 CONTINUE
      FILENAME(IPT+1:IPT+6) = SUFFIX
 
c  debug write commented out, dmc 22 Jan 1996
cx        WRITE(LUOUT, 9555) IPT+6, filename(1:IPT+6)
cx 9555	Format('TRMFLNAM: filename length =', I6 / ' --->', A, '<---')
 
      RETURN
 
C	-----------------------------------------------------------------
C	.ERROR RETURN - NEW FILENAME TOO LONG.
 
  999   WRITE (LUOUT, 9000, ERR=998) I,DISK,DIR,RUNID,
     1		  SUFFIX, FILENAME
 9000   FORMAT(/ ' ?TRMFLNAM: NEW FILE NAME > 64 ',
     1           'CHARACTERS. LENGTH = ',
     1           I5 / 1X, A64 / 1X, A64 / 1X, A10 / 1X, A6 / 1X, A64//)
  998   CONTINUE
        RETURN
 
C	****************** END FILE  TRMFLNAM.FOR  ;  GROUP TRREAD_LIB ******
 
      END

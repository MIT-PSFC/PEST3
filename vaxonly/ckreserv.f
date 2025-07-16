C-----------------------------------------------------------------------
C  CKRESERV:  CHECK NODE RESERVATION STATUS
C
C   CALL TO SEE IF A GIVEN NODE IS RESERVED AND CANNOT BE USED BY
C   A GIVEN QUEUED TRANSP JOB
C
      SUBROUTINE CKRESERV(ZRESERV,ZNAME1,ZNAME2,ISTAT)
C
C  INPUT
C    ZRESERV -- RESERVATION STRING, BLANK MEANS NO RESERVATIONS
C        ZRESERV contains comma or blank delimitted list of NAMEs or ~NAMEs.
C        *caution* only one delimitter character btw names allowed!
C        *noted* dmc 10 Dec 1997
C    ZNAME1,ZNAME2 -- TOKAMAK AND OWNER NAME ASSOCIATED WITH CURRENT JOB
C
C    if ZRESERV is blank, any job can use node
C    if ZRESERV is non-blank, job's associated names must match one
C       of the NAMEs in ZRESERV, and must not match any of the ~NAMEs
C       in ZRESERV.
C
C  OUTPUT
C    ISTAT=0:  JOB IS FREE TO USE NODE
C    ISTAT=1:  JOB MAY NOT USE NODE -- RESERVED FOR OTHER OWNER OR
C                   TOKAMAK NAMES
C
      implicit NONE
C
      CHARACTER*(*) ZRESERV,ZNAME1,ZNAME2
C
      INTEGER ISTAT
C
C-----------------------------------
      integer istat1,ict1,istat2,ilen,ic,ilenb,ic1,ic2,ic1p
C-----------------------------------
C
      call uupper(zreserv)
      call uupper(zname1)
      call uupper(zname2)
C
      ISTAT=0
C
      IF(ZRESERV.EQ.' ') RETURN
C
C  NON-BLANK RESERVATION STRING; NODE MAY BE RESERVED
C
C  a name (ZNAME1 or ZNAME2) must match one of the names
C  (not starting with ~) in ZRESERV.  If there are such names
C  in ZRESERV, ICT1 gets incremented.  If one of these names
C  matches, ISTAT1 gets cleared.
C
      ISTAT1=1
        ICT1=0
C
C  a name (ZNAME1 or ZNAME2) must *not* match names given in the
C  form ~NAME in the ZRESERV string.  If a match of this type is
C  detected, ISTAT2 gets set
        ISTAT2=0
C
C  find, nonblank length of ZRESERV, store in ILENB
      ILEN=LEN(ZRESERV)
      DO IC=ILEN,1,-1
        IF(ZRESERV(IC:IC).NE.' ') GO TO 10
      ENDDO
C
 10   CONTINUE
      ILENB=IC
C
C-------------------
C  loop through names...
C
      IC1=1
 20   CONTINUE
      DO IC=IC1+1,ILENB
        IF(ZRESERV(IC:IC).EQ.',') GO TO 30
          IF(ZRESERV(IC:IC).EQ.' ') GO TO 30
      ENDDO
      IC=ILENB+1
C
 30   CONTINUE
      IC2=IC-1
C
        if(zreserv(ic1:ic1).ne.'~') then
           ict1=ict1+1
           IF(ZNAME1.EQ.ZRESERV(IC1:IC2)) ISTAT1=0
           IF(ZNAME2.EQ.ZRESERV(IC1:IC2)) ISTAT1=0
        else
           ic1p=ic1+1
           if(ic1p.le.ic2) then
              IF(ZNAME1.EQ.ZRESERV(IC1P:IC2)) ISTAT2=1
              IF(ZNAME2.EQ.ZRESERV(IC1P:IC2)) ISTAT2=1
           endif
        endif
C
      IC1=IC2+2
      IF(IC1.GT.ILENB) go to 100
C
      GO TO 20  ! re-enter loop
C
C---------------------
C  exit; return "LOGICAL AND" of both conditions:
C
 100    continue
C
        if(ict1.eq.0) istat1=0
        istat = istat1 + istat2
C
        return
C
      END

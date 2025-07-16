C******************** START FILE IRNS.FOR ; GROUP RNG ******************
      FUNCTION IRNS(IDUMMY)
C  11/29/78
C  THIS FUNCTION WAS WRITTEN BY HARRY H. TOWNER OF THE PRINCETON PLASMA
C  PHYSICS LAB.  IRNS WILL RETURN A UNIQUE INTEGER SUITABLE FOR SETTING
C  THE RANDOM NUMBER SEED.
c  5/01 RGA -- replace with f90 standard call
      integer iv(8)
      call date_and_time(values=iv)
      secnds = 60*(60*iv(5)+iv(6))+iv(7)+iv(8)/1000.
 
      IRNS=1.E4*SECNDS
      IRNS = 2*IRNS + 1	!MAKE SURE IT'S ODD, RTM 26 FEB 86
      RETURN
      END
C******************** END FILE IRNS.FOR ; GROUP RNG ******************

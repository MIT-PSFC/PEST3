C******************** START FILE IERFCN.FOR ; GROUP FILTR6 ******************
C=================================================================
C  IERFCN   TEST Y(J)-YSM(J)
C   IERFCN=1==> EPSLON VIOLATION
C   IERFCN=0==> NO VIOLATION
C
C   ISIGN=SIGN OF VIOLATION
C   DY=MAGNITUDE OF VIOLATION
C
      INTEGER FUNCTION IERFCN(Y,YSM,N,JP,EPS,NE,DY,ISIGN)
      DIMENSION Y(N),YSM(N),EPS(NE)
C
      IERFCN=0
      ZEPS=EPS(NE)
      IF(JP.LT.NE) ZEPS=EPS(JP)
      DY=Y(JP)-YSM(JP)
      ISIGN=SIGN(1.0,DY)
      DY=ABS(DY)
      IF(DY.LE.ZEPS) RETURN
      IERFCN=1
      RETURN
      	END
C******************** END FILE IERFCN.FOR ; GROUP FILTR6 ******************

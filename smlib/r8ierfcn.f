C******************** START FILE IERFCN.FOR ; GROUP FILTR6 ******************
C=================================================================
C  IERFCN   TEST Y(J)-YSM(J)
C   IERFCN=1==> EPSLON VIOLATION
C   IERFCN=0==> NO VIOLATION
C
C   ISIGN=SIGN OF VIOLATION
C   DY=MAGNITUDE OF VIOLATION
C
      INTEGER FUNCTION r8ierfcn(Y,YSM,N,JP,EPS,NE,DY,ISIGN)
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n,jp,ne,isign
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 ysm,eps,dy,y,zeps
C============
      DIMENSION Y(N),YSM(N),EPS(NE)
C
      r8ierfcn=0
      ZEPS=EPS(NE)
      IF(JP.LT.NE) ZEPS=EPS(JP)
      DY=Y(JP)-YSM(JP)
      ISIGN=SIGN(1.0E0_R8,DY)
      DY=ABS(DY)
      IF(DY.LE.ZEPS) RETURN
      r8ierfcn=1
      RETURN
      	END
C******************** END FILE IERFCN.FOR ; GROUP FILTR6 ******************
! 03Dec1999 fgtok -s r8_precision.sub r8smlib.sub "r8con.csh conversion"

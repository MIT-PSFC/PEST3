C******************** START FILE WXCINT.FOR ; GROUP FILTR6 ******************
C==========================================================
C
C  WXCINT
C
C  COMPUTE C*INTEGRAL OF WEIGHTING FUNCTION FROM X1 TO X2
C
C  DENOMINATOR IN WEIGHTED AVERAGE FORMULA
C
C  W(T)=1-ABS((XCEN-T)/DELTA)   XCEN-DELTA .LE. T .LE. XCEN+DELTA
C
C  TRIANGULAR SHAPED WEIGHTING FUNCTION CENTERED AT XCEN
C    X1 SHOULD BE .LT. X2 AND X1 AND X2 SHOULD
C    BE WITHIN DELTA OF XCEN
C-------------
      REAL*8 FUNCTION r8wxcint(X1,X2,C,XCEN,DELTA)
C
C============
C idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 x2,c,xcen,delta,x1,xa,xb,zdi
C============
      r8wxcint=C
      XA=max(-DELTA,(X1-XCEN))
      XB=min(DELTA,(X2-XCEN))
      IF((XA.GT.DELTA).OR.(XB.LT.-DELTA)) GO TO 1000
      IF(XB.LT.XA) GO TO 1000
C
C  NUMERICS CHECK; DELTA <<< X1
      IF(XCEN.EQ.(XCEN+DELTA)) RETURN
C
      ZDI=1.E0_R8/DELTA
      IF((XA.GT.0.E0_R8).OR.(XB.LT.0.E0_R8)) GO TO 100
C
      r8wxcint=C*ZDI*((XB-XA)-.5E0_R8*(XA*XA+XB*XB)*ZDI)
      RETURN
C
 100  CONTINUE
      IF(XA.GT.0.E0_R8)
     >   r8wxcint=C*ZDI*((XB-XA)-.5E0_R8*(XA+XB)*(XB-XA)*ZDI)
      IF(XB.LT.0.E0_R8)
     >   r8wxcint=C*ZDI*((XB-XA)-.5E0_R8*(XA+XB)*(XA-XB)*ZDI)
      RETURN
C
 1000 CONTINUE
      write(6,9000)
 9000 FORMAT(' SUBROUTINE WXCINT FROM FILFN6 FROM FILTR6.FOR:'/
     >         ' ILLEGAL INTEGRATION LIMITS')
      call bad_exit
      	END
C******************** END FILE WXCINT.FOR ; GROUP FILTR6 ******************

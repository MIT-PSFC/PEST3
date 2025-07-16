C******************** START FILE FILFN6.FOR ; GROUP FILTR6 ******************
C
C
C------------------------------------------------------------
C  FILFN6
C
C  COMPUTE WEIGHTED AVERAGE ABOUT PT X(J0).
C
C  END CONDITIONS:
C     XEND?=0==> FIX END PT. OF REGION BEING SMOOTHED.
C                THE DATA CURVE IS REFLECTED ABOUT THE PT.
C                [X(1),Y(1)] ( OR [X(N),Y(N)])       (E0)
C               SO THAT BY SYMMETRY OF THE WEIGHTING
C                FUNCTION ==> THE END PTS ARE FIXED
C
C      XEND?=1==> DON'T FIX END PT. OF REGION. END PT
C                IS RESET TO WEIGHTED AVERAGE OF NEARBY PTS.
C                THE DATA IS REFLECTED ABOUT THE LINE X=X(1)
C                FOR THE LEFT END, ABOUT THE LINE X=X(N)
C                FOR  THE RIGHT END.                 (E1)
C
C
C  IF XEND DOES NOT EQUAL 0 OR 1 XEND SHOULD BE BETWEEN 0 AND 1;
C   THEN THE 'EXTENSION FORMULA' IS E=(1-XEND)*E0+XEND*E1
C---------------------------------------------
C
      REAL*8 FUNCTION r8filfn6(X,Y,N,J0,DELTA0,XEND1,XEND2,ICOD)
C
C  DMC NOV 1985 - ICOD=1 : NORMAL TRIANGULAR WEIGHTED AVG INTEGRAL
C
C    ICOD=2 : BROKEN TWO PIECE TRIANGLE FOR SAWTOOTH CORRELATION
C             INTEGRAL
C
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n,j0,icod,j,jr
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 y,delta0,xend1,xend2,x,delta,xlim1,xlim2,denom,r8wxcint,
     > zsign,wasum,xa,ya,xb,yb,ya0,ya1,xa1,r8wxfint,yb0,yb1,xb1
C============
      DIMENSION X(N),Y(N),DELTA0(N)
C
C  START OF EXECUTABLE CODE --
C
C
C  CHECK FOR NUMERICALLY ZERO OR NEGATIVE DELTA
      DELTA=min((X(N)-X(1)),DELTA0(J0))
      IF((DELTA.LE.0.0E0_R8).OR.((X(J0)+DELTA).EQ.X(J0))) THEN
        IF(ICOD.EQ.1) THEN
          r8filfn6=Y(J0)
        ELSE
          r8filfn6=0.0E0_R8
        ENDIF
        RETURN
      ENDIF
C
      IF(ICOD.EQ.1) THEN
C  CHECK FOR FIXED ENDPOINTS
        r8filfn6=Y(J0)
        IF((J0.EQ.1).AND.(XEND1.EQ.0.0E0_R8)) RETURN
        IF((J0.EQ.N).AND.(XEND2.EQ.0.0E0_R8)) RETURN
      ENDIF
C
C  EVALUATE DENOMINATOR OF WEIGHTED AVERAGE FORMULA
C
      XLIM1=X(J0)-DELTA
      XLIM2=X(J0)+DELTA
C
      DENOM=r8wxcint(XLIM1,XLIM2,1.E0_R8,X(J0),DELTA)
C
C  EVALUATE NUMERATOR OF WGTED. AVG.
C
C  (A) PTS TO THE LEFT OF X(J0)
C
      ZSIGN=1.0E0_R8
      WASUM=0.E0_R8
C
      J=J0
      XA=X(J0)
      YA=Y(J0)
 10   J=J-1
      XB=XA
      YB=YA
      IF (J) 11,11,12
C
C  REFLECT PTS ABOUT X(1),Y(1)
C
 11   CONTINUE
      JR=min(N,(2-J))
      IF((X(J0)-XB).GE.DELTA.or.(J.LE.(1-N))) GO TO 50
      XA=X(1)-(X(JR)-X(1))
      YA0=Y(1)-(Y(JR)-Y(1))
      YA1=Y(JR)
      YA=(1.E0_R8-XEND1)*YA0+XEND1*YA1
      GO TO 15
C
 12   CONTINUE
      IF((X(J0)-XB).GE.DELTA) GO TO 50
      XA=X(J)
      YA=Y(J)
C
C  COMPUTE WEIGHTED AVERAGE INTEGRAL FROM XA1 TO XB
C
 15   CONTINUE
      XA1=max(XA,(X(J0)-DELTA))
      WASUM=WASUM+ZSIGN*r8wxfint(XA1,XB,XA,YA,XB,YB,X(J0),DELTA)
C
      GO TO 10
C
C (B) COMPUTE CONTRIBUTION OF DATA TO RIGHT OF X(J0)
C
 50   CONTINUE
C
      IF(ICOD.EQ.2) ZSIGN=-1.0E0_R8
C
      J=J0-1
      XB=X(J0)
      YB=Y(J0)
 60   J=J+1
      XA=XB
      YA=YB
      IF(N-J) 61,61,62
C
C  REFLECT ABOUT X(N),Y(N)
C
 61   CONTINUE
      JR=max(1,(N+N-J-1))
      IF((XA-X(J0)).GE.DELTA.or.(J.GE.(2*N-1))) GO TO 100
      XB=X(N)+(X(N)-X(JR))
      YB0=Y(N)+(Y(N)-Y(JR))
      YB1=Y(JR)
      YB=(1.E0_R8-XEND2)*YB0+XEND2*YB1
      GO TO 65
C
 62   CONTINUE
      IF((XA-X(J0)).GE.DELTA) GO TO 100
      XB=X(J+1)
      YB=Y(J+1)
C
C  INTEGRATE FROM XA TO XB1
C
 65   XB1=min(XB,(X(J0)+DELTA))
C
      WASUM=WASUM+ZSIGN*r8wxfint(XA,XB1,XA,YA,XB,YB,X(J0),DELTA)
C
      GO TO 60
C
C
C  RETURN COMPUTED WEIGHTED AVERAGE
C
 100  CONTINUE
C
      r8filfn6=WASUM/DENOM
C
      RETURN
      	END
C******************** END FILE FILFN6.FOR ; GROUP FILTR6 ******************
! 03Dec1999 fgtok -s r8_precision.sub r8smlib.sub "r8con.csh conversion"

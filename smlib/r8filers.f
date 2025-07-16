C******************** START FILE FILERS.FOR ; GROUP FILTR6 ******************
C=============================================================
C  FILERS    SEEK OUT NEXT REGION WHERE SMOOTHED DATA VIOLATES
C            EPSLON ERROR BARS ON ORIGINAL DATA
C
C  JS=LOCATION AT WHICH TO START SEARCH
C
C  IND(6)  INDEX DESCRIPTION OF ERROR  IND(1)=0 ==> NO ERROR
C          FOUND  (OUTPUT)
C          IND(3),IND(4) DELIMIT ACTUAL VIOLATIONS WITHIN A
C          REGION WHERE THE DIVERGENCE (Y-YSM) IS OF THE SAME
C          SIGN   (X,Y (DIM(N)) ARE RAW DATA, YSM IS SMOOTHED DATA)
C          IND(2),IND(5) DELIMIT REGION WHERE (Y-YSM) IS OF THE
C          SAME SIGN
C          IND(1),IND(6) DELIMIT REGION WITHIN DELTA(IND(2))
C            OF X(IND(3)) (ON THE RIGHT); AND DELTA(IND(5))
C            OF X(IND(4)) (ON THE LEFT) ;  THIS IS THE REGION
C            WHERE A CORRECTION WILL BE APPLIED
C
C  EPS(NE) ARRAY CONTAINS ERROR BARS FOR POINTS 1 THRU NE; IF
C          J>NE, EPS(NE) IS USED AS THE ERROR BAR FOR
C          THE JTH DATA POINT.
C
C   EPS2(NE2):  IF A REGION OF SAME SIGN DIVERGENCE (Y-YSM)
C              CONTAINS ND LESS THAN NE2 DATA POINTS, THEN THE
C            ERROR BAR WILL BE MULTIPLIED BY THE FACTOR
C              EPS2(ND) ... THIS ALLOWS A LESS STIFF ERROR
C              BAR CRITERION FOR VIOLATIONS OF SHORT WIDTH
C
C----------------------------
C
      SUBROUTINE r8filers(X,Y,YSM,N,DELTA,EPS,NE,EPS2,NE2,IND,JS)
C
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n,ne2,ind,js,j,jlp,ierrf,isign0,je1,je2,js2,jd2,js1,jd1,
     > j2,jlp2,nse,jst,ne,isign
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 y,ysm,delta,eps,eps2,x,dymax,dy,x0,zemax,zetes
C============
      DIMENSION X(N),Y(N),YSM(N),DELTA(N),EPS(NE),EPS2(1),
     >             IND(6)
C
      integer r8ierfcn
C
C  CLEAR ERROR FLAG
C
      IND(1)=0
C
      IF(JS.GT.N) RETURN
C  FIND FIRST EPSLON VIOLATION
      DO 100 J=JS,N
      DYMAX=0.E0_R8
      JLP=J
      IERRF=r8ierfcn(Y,YSM,N,JLP,EPS,NE,DY,ISIGN0)
      IF(IERRF.EQ.0) GO TO 100
C  POSSIBLE ERROR DETECTED
      DYMAX=max(DYMAX,DY)
      IF(J.LT.N) GO TO 10
      JE1=N
      JE2=N
      JS2=N
      JD2=N
      GO TO 40
C
 10   CONTINUE
      JE1=J
      IF(JE1.GT.1) GO TO 20
      JS1=1
      JD1=1
C  FIND  JE2=IND(4),  JS2=IND(5)
 20   CONTINUE
      JE2=JE1
      DO 30 J2=JE1,N
      JLP2=J2
      IERRF=r8ierfcn(Y,YSM,N,JLP2,EPS,NE,DY,ISIGN)
      DYMAX=max(DYMAX,DY)
      IF(ISIGN.NE.ISIGN0) GO TO 35
      IF(IERRF.EQ.1) JE2=J2
 30   CONTINUE
      J2=N+1
 35   CONTINUE
      JS2=J2-1
C  FIND JD2 (IND(6))
      X0=X(JE2)
      DO 37 J2=JE2,N
      IF((X(J2)-X0).GT.DELTA(J2)) GO TO 38
 37   CONTINUE
      J2=N
 38   CONTINUE
      JD2=J2
C
      IF(JE1.EQ.1) GO TO 60
C  FIND JS1  (IND(2))
 40   CONTINUE
      DO 50 J2=JE1-1,1,-1
      JLP2=J2
      IERRF=r8ierfcn(Y,YSM,N,JLP2,EPS,NE,DY,ISIGN)
      IF(ISIGN.NE.ISIGN0) GO TO 55
 50   CONTINUE
      J2=0
 55   CONTINUE
      JS1=J2+1
C  FIND JD1 (IND(1))
      X0=X(JE1)
      DO 57 J2=JE1,1,-1
      IF((X0-X(J2)).GT.DELTA(J2)) GO TO 58
 57   CONTINUE
      J2=1
 58   CONTINUE
      JD1=J2
C  DOES EPS2 ARRAY COME TO BEAR ?
 60   CONTINUE
      NSE=JS2-JS1+1
      IF(NSE.GT.NE2) GO TO 80
C  MAYBE
      ZEMAX=EPS(NE)*EPS2(NSE)
      IF(JE1.GE.NE) GO TO 75
      JST=min(JE2,NE)
      ZEMAX=0.E0_R8
      DO 70 J2=JE1,JST
      ZETES=EPS(J2)*EPS2(NSE)
      ZEMAX=max(ZEMAX,ZETES)
 70   CONTINUE
 
 75   CONTINUE
      IF(ZEMAX.GE.DYMAX) GO TO 100
C  ERROR DEFINATELY DETECTED
 80   CONTINUE
      IND(1)=JD1
      IND(2)=JS1
      IND(3)=JE1
      IND(4)=JE2
      IND(5)=JS2
      IND(6)=JD2
      GO TO 200
C
C  END OF LOOPT
C
 100  CONTINUE
 200  CONTINUE
      RETURN
      	END
C******************** END FILE FILERS.FOR ; GROUP FILTR6 ******************
! 03Dec1999 fgtok -s r8_precision.sub r8smlib.sub "r8con.csh conversion"

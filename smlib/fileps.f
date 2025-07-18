C******************** START FILE FILEPS.FOR ; GROUP FILTR6 ******************
C------------------------------------------------------------
C   FILEPS   IMPLIMENT ERROR BAR FEATURE OF ALGORITHM
C
C   THE SMOOTHING RADIUS ARRAY DELTA IS LOCALLY
C   ADJUSTED WHEREVER THE SMOOTHED DATA FUNCTION IS
C   OUTSIDE THE ERROR BARS OF THE ORIGINAL DATA, UNTIL
C   ALL ERROR BAR VIOLATIONS ARE REMOVED
C
C   THE DATA INITIALIZED VARIABLE RITER IS THE FACTOR BY
C   WHICH THE LOCAL DELTA IS TO BE REDUCED WHERE VIOLATIONS
C   ARE DETECTED.  RITER SHOULD BE SET BETWEEN 0. AND 1.
C   VALUES CLOSER TO ONE MAY PRODUCE BETTER SMOOTHING EFECT,
C   VALUES CLOSER TO ZERO WILL SPEED THE ALGORITHM.
C
C
      SUBROUTINE FILEPS(X,Y,YSM,N,DELTA,EPS,NE,EPS2,NE2,
     >                    XEND1,XEND2)
C
      DIMENSION X(N),Y(N),YSM(N),DELTA(N),EPS(NE) ,EPS2(1)
C
      DIMENSION IND0(6),IND(6)
C
C  THE "IND" ARRAYS WILL CONTAIN INFORMATION ON
C   ANY ERROR BAR VIOLATIONS.
C   THEY POINT TO LOCATIONS IN X,Y, AND YSM.
C    IND(3)=STARTING POINT OF AN ERROR BAR VIOLATION
C    IND(4)=END PT OF SAME VIOLATION IN CONTIGUOUS REGION
C           OF DIVERGENCE (Y-YSM) OF SAME SIGN
C    IND(2)=START POINT OF REGION OF SAME SIGN DIVERGENCE
C    IND(5)=END POINT OF REGION
C    IND(1)=LOCATION OF REGION DELTA(IND(3)) DISTANT FROM X(IND(3))
C           DEFINING END OF REGION WHERE DATA WILL BE RESMOOTHED
C    IND(6)=SAME AS IND(1), RIGHT ENDPT.
C  IND IS GENERALLY IN INCREASING ORDER
C
C
      DATA RITER/.9/
C  RITER= CONVERGENCE FACTOR IN REDIFINING DELTA ARRAY
C
      DATA XISOL/1.0/
C  XISOL= ISOLATION PARAMETER IN DETERMINING SEPERATION
C         OF NEIGHBOURING REGIONS OF EPSILON VIOLATIONS
C
C
C  START OF EXECUTABLE CODE
C
      JS0=1
C
      DO 6 IJ=1,6
      IND(IJ)=0
      IND0(IJ)=0
 6    CONTINUE
 10   CONTINUE
C  SEEK ERROR
CT	TYPE 9901
 9901 FORMAT(' FILERS--A')
      CALL FILERS(X,Y,YSM,N,DELTA,EPS,NE,EPS2,NE2,IND0,JS0)
CT	TYPE *,IND0
C
C  IND0(0)=0==> NO MORE ERRORS
C
      IF(IND0(1).EQ.0) RETURN
C
      JC1=IND0(3)
C  CHECK RELATIVE LOCATION OF LEFT ENDPOINT OF DATA
      IF(JS0.GT.1) GO TO 20
      D1=AMAX1(DELTA(IND0(1)),DELTA(1))
      IF((X(IND0(1))-X(1)).LT.(XISOL*(D1))) JC1=1
C
 20   CONTINUE
C  SEEK NEXT ERROR
 25   CONTINUE
      JS=IND0(5)+1
CT	TYPE 9902
 9902 FORMAT(' FILERS--B')
      CALL FILERS(X,Y,YSM,N,DELTA,EPS,NE,EPS2,NE2,IND,JS)
CT	TYPE *,IND
      IF(IND(1).NE.0) GO TO 30
C  NO NEXT ERROR
      JC2=IND0(4)
C  CHECK RIGHT ENDPT OF DATA
      D1=AMAX1(DELTA(IND0(6)),DELTA(N))
      IF((X(N)-X(IND0(6))).LT.(XISOL*(D1))) JC2=N
      GO TO 90
C  YES NEXT ERROR:  JOINT TREATMENT WITH LAST ERROR ?
 30   CONTINUE
      JC2=IND0(4)
      D2=AMAX1(DELTA(IND0(6)),DELTA(IND(1)) )
      IF((X(IND(1))-X(IND0(6))).GE.(XISOL*(D2))) GO TO 90
C  YES:  JOINT TREATMENT
      DO 40 I=4,6
      IND0(I)=IND(I)
 40   CONTINUE
C  CHECK FOR ANOTHER ERROR
      GO TO 25
C
C  ERROR TREATMENT
C
 90   CONTINUE
C
CT	TYPE *,IND0
CT	TYPE *,JC1,JC2
CT	PAUSE 'IND0,JC1,JC2'
      JD1=JC1-1
      IF(JC1.EQ.1) GO TO 95
      JD1=IND0(1)
      IF(JD1.GE.(JC1-1) ) GO TO 95
      ZS=X(JC1)-X(JC1-1)
      ZDS=DELTA(JC1)-ZS
      ZDSB=DELTA(JC1)*RITER-ZS
      IF(ZDSB.LE.0.) GO TO 95
      ZP=ZDSB/(RITER*ZDS)
      DO 94 J=JD1+1,JC1-1
      ZS=X(JC1)-X(J)
      ZDS=DELTA(J)-ZS
      ZZ=ZP*ZDS/DELTA(J)
      RITERL=ZS*(1.+ZZ)/DELTA(J)
      RITERM=AMIN1(1.,AMAX1(RITERL,RITER))
CT	TYPE *,J,JC1,RITERL,RITERM
      DELTA(J)=DELTA(J)*RITERM
 94   CONTINUE
CT	PAUSE 'TO THE LEFT'
C
 95   CONTINUE
      DO 96 J=JC1,JC2
      DELTA(J)=DELTA(J)*RITER
 96   CONTINUE
C
      JD2=JC2+1
      IF(JC2.EQ.N) GO TO 100
      JD2=IND0(6)
      IF(JD2.LE.(JC2+1)) GO TO 100
      ZS=X(JC2+1)-X(JC2)
      ZDS=DELTA(JC2)-ZS
      ZDSB=RITER*DELTA(JC2)-ZS
      IF(ZDSB.LE.0.) GO TO 100
      ZP=ZDSB/(RITER*ZDS)
      DO 97 J=JC2+1,JD2-1
      ZS=X(J)-X(JC2)
      ZDS=DELTA(J)-ZS
      ZZ=ZP*ZDS/DELTA(J)
      RITERL=ZS*(1.+ZZ)/DELTA(J)
      RITERM=AMIN1(1.,AMAX1(RITER,RITERL))
CT	TYPE *,J,JC2,RITERL,RITERM
      DELTA(J)=RITERM*DELTA(J)
 97   CONTINUE
CT	PAUSE 'TO THE RIGHT'
C
 100  CONTINUE
C  REEVALUATE SMOOTHED DATA IN REGION OF ERROR BAR VIOLATION(S)
      J1=MIN0(JC1,(JD1+1))
      J2=MAX0(JC2,(JD2-1))
CT	TYPE *,J1,J2,N
CT	PAUSE 'J1,J2,N'
      DO 110 J=J1,J2
C  DO NOT ALLOW DELTA(J) TO GET TOO SMALL -- RESET TO ZERO FOR SPEED
      JLP=J
      IF(JLP.EQ.1) ZDX=X(2)-X(1)
      IF(JLP.EQ.N) ZDX=X(N)-X(N-1)
      IF((JLP.GT.1).AND.(JLP.LT.N)) ZDX=0.5*(X(JLP+1)-X(JLP-1))
      IF(DELTA(J).LT.0.01*ZDX) DELTA(J)=0.0
      YSM(J)=FILFN6(X,Y,N,JLP,DELTA,XEND1,XEND2,1)
 110  CONTINUE
CT	CALL GRFADH(X,DELTA,DELTA,N,-11,0,
CT     >  30HSMOOTHING DELTA                  ,
CT     >  30H   X                             ,
CT    >  10H (FILEPS) ,0,30H                               )
      GO TO 10
C
C
      	END
C******************** END FILE FILEPS.FOR ; GROUP FILTR6 ******************

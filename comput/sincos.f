C******************** START FILE SINCOS.FOR ; GROUP TKODE2 ******************
C.............................................................
 
      SUBROUTINE SINCOS(ZTHETA,JANG,SNTHTK,CSTHTK)
 
C
C	THIS SUBROUTINE CALCULATES:
C		SIN(ZTHETA)
C		COS(ZTHETA)
C		SIN(2*ZTHETA)
C		COS(2*ZTHETA)
C	ETC., UP TO JANG*ZTHETA
C
C
C		ZTHETA IS INPUT ANGLE
C		JANG IS HIGHEST N*ZTHETA TO CALCULATE
C		SNTHTK IS THE ARRAY CONTAINING SIN(N*ZTHETA)
C		CSTHTK IS THE ARRAY CONTAINING COS(N*ZTHETA)
C
      DIMENSION SNTHTK(JANG),CSTHTK(JANG)
C
C  LOCAL MEMORY (DMC 6 JUL 1994)
C
      DATA ZTHETAP/0.0/
      DATA ZSINP/0.0/
      DATA ZCOSP/1.0/
C
      SAVE ZTHETAP,ZSINP,ZCOSP
C
C--------------------------------------------------------------------
C
C  DMC -- USE LOCAL MEMORY FOR SPEED
C
      IF(ZTHETA.NE.ZTHETAP) THEN
C
C  EVALUATE SIN,COS
C
        SNTHTK(1)=SIN(ZTHETA)
        CSTHTK(1)=COS(ZTHETA)
        ZTHETAP=ZTHETA
        ZSINP=SNTHTK(1)
        ZCOSP=CSTHTK(1)
      ELSE
C
C  REUSE PREVIOUS RESULTS
C
        SNTHTK(1)=ZSINP
        CSTHTK(1)=ZCOSP
      ENDIF
C
      DO 100 I=2,JANG
        SNTHTK(I)=SNTHTK(I-1)*CSTHTK(1)+CSTHTK(I-1)*SNTHTK(1)
        CSTHTK(I)=CSTHTK(I-1)*CSTHTK(1)-SNTHTK(I-1)*SNTHTK(1)
100   CONTINUE
C
      RETURN
      END
C
C******************** END FILE SINCOS.FOR ; GROUP TKODE2 ******************

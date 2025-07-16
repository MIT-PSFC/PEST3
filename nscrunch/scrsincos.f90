!******************** START FILE SINCOS.FOR ; GROUP TKODE2 ******************
!.............................................................
 
      SUBROUTINE scrSINCOS(ZTHETA,JANG,SNTHTK,CSTHTK)
      IMPLICIT NONE
       INTEGER 	i ! <sincos.f>
       REAL*8 	csthtk ! <sincos.f>
       REAL*8 	snthtk ! <sincos.f>
       REAL*8 	ztheta ! <sincos.f>
       REAL*8 	zcosp ! <sincos.f>
       REAL*8 	zsinp ! <sincos.f>
       REAL*8 	zthetap ! <sincos.f>
       INTEGER 	jang ! <sincos.f>
 
!
!	THIS SUBROUTINE CALCULATES:
!		SIN(ZTHETA)
!		COS(ZTHETA)
!		SIN(2*ZTHETA)
!		COS(2*ZTHETA)
!	ETC., UP TO JANG*ZTHETA
!
!
!		ZTHETA IS INPUT ANGLE
!		JANG IS HIGHEST N*ZTHETA TO CALCULATE
!		SNTHTK IS THE ARRAY CONTAINING SIN(N*ZTHETA)
!		CSTHTK IS THE ARRAY CONTAINING COS(N*ZTHETA)
!
	DIMENSION SNTHTK(JANG),CSTHTK(JANG)
!
!  LOCAL MEMORY (DMC 6 JUL 1994)
!
	DATA ZTHETAP/0.0D0/
	DATA ZSINP/0.0D0/
	DATA ZCOSP/1.0D0/
!
	SAVE ZTHETAP,ZSINP,ZCOSP
!
!--------------------------------------------------------------------
!
!  DMC -- USE LOCAL MEMORY FOR SPEED
!
	IF(ZTHETA.NE.ZTHETAP) THEN
!
!  EVALUATE SIN,COS
!
	  SNTHTK(1)=SIN(ZTHETA)
	  CSTHTK(1)=COS(ZTHETA)
	  ZTHETAP=ZTHETA
	  ZSINP=SNTHTK(1)
	  ZCOSP=CSTHTK(1)
	ELSE
!
!  REUSE PREVIOUS RESULTS
!
	  SNTHTK(1)=ZSINP
	  CSTHTK(1)=ZCOSP
	ENDIF
!
	DO 100 I=2,JANG
	  SNTHTK(I)=SNTHTK(I-1)*CSTHTK(1)+CSTHTK(I-1)*SNTHTK(1)
	  CSTHTK(I)=CSTHTK(I-1)*CSTHTK(1)-SNTHTK(I-1)*SNTHTK(1)
100	CONTINUE
!
	RETURN
	END
!
!******************** END FILE SINCOS.FOR ; GROUP TKODE2 ******************
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"

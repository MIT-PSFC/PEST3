C******************** START FILE ZKBOLT.FOR ; GROUP KAPAI ******************
C
C
C@@@
C  BOLTONS KI(NU-HAT,EPS) FUNCTION (THESIS PP 74--78)
C
      REAL*8 FUNCTION ZKBOLT(ZNUHAT,ZEPS)
C
C  N-C COEFFICIENTS FROM TABLE VI
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 zeps,znuhat,zanc,zbnc,zcnc,zdnc,zaps,zbps,zcps,zdps
      REAL*8 zse,zeps32,ztanc,ztbnc,ztcnc,ztdnc,znh,zkinc,ztaps
      REAL*8 ztbps,ztcps,ztdps,zkips
!============
      DIMENSION ZANC(3),ZBNC(3),ZCNC(3),ZDNC(3)
C  P-S COEFFICIENTS FROM TABLE VI
      DIMENSION ZAPS(3),ZBPS(3),ZCPS(3),ZDPS(3)
C
      DATA ZANC/2.441_R8,-3.87_R8,2.19_R8/
      DATA ZBNC/.9362_R8,-3.109_R8,4.087_R8/
      DATA ZCNC/.241_R8,3.40_R8,-2.54_R8/
      DATA ZDNC/.2664_R8,-.352_R8,.44_R8/
C
      DATA ZAPS/.364_R8,-2.76_R8,2.21_R8/
      DATA ZBPS/.553_R8,2.41_R8,-3.42_R8/
      DATA ZCPS/1.18_R8,.292_R8,1.07_R8/
      DATA ZDPS/.0188_R8,.180_R8,-.127_R8/
C
      ZSE=SQRT(ZEPS)
      ZEPS32=ZSE*ZSE*ZSE
      ZTANC=.66_R8+((ZANC(3)*ZSE+ZANC(2))*ZSE+ZANC(1))*ZSE
      ZTBNC=((ZBNC(3)*ZEPS+ZBNC(2))*ZEPS+ZBNC(1))/(ZEPS**.75_R8)
      ZTCNC=((ZCNC(3)*ZEPS+ZCNC(2))*ZEPS+ZCNC(1))/(ZEPS32)
      ZTDNC=((ZDNC(3)*ZEPS+ZDNC(2))*ZEPS+ZDNC(1))/(ZEPS32)
C
C  N-C PART
C
      ZNH=SQRT(ZNUHAT)
      ZKINC=ZTANC/
     >  (1._R8+ZTBNC*ZNH+ZTCNC*ZNUHAT+ZTDNC*ZNUHAT*ZNUHAT)
C
      ZTAPS=((ZAPS(3)*ZSE+ZAPS(2))*ZSE+ZAPS(1))*ZEPS32
      ZTBPS=((ZBPS(3)*ZSE+ZBPS(2))*ZSE+ZBPS(1))*ZEPS32
      ZTCPS=((ZCPS(3)*ZSE+ZCPS(2))*ZSE+ZCPS(1))
      ZTDPS=((ZDPS(3)*ZSE+ZDPS(2))*ZSE+ZDPS(1))
C
C  P-S PART
C
      ZKIPS=1.57_R8*ZEPS32 + (ZTAPS+ZTBPS*ZNH)/
     >          (1._R8+ZTCPS*ZNUHAT**1.5_R8+ZTDPS*ZNUHAT**2.5_R8)
C
      ZKBOLT=ZKINC+ZKIPS
C
      RETURN
      	END
C******************** END FILE ZKBOLT.FOR ; GROUP KAPAI ******************
! 20jan2003 fgtok -s r8_precision.sub all.sub "r8con.csh conversion"

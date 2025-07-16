C******************** START FILE SPLCRD.FOR ; GROUP TKSERVE ************
C.............................................................
 
        SUBROUTINE r8splcrd (ZXI99, IJACI,
     >                     LCENTR, MIMOM, NMOM, MTBL, XIBLO,
     >                     RMB0,RMB0S,RMB2,RMB2S,YMB2,YMB2S,
     >			   ZRMB0X, ZRMB2X, ZYMB2X,
     >		           ZRMB0P, ZRMB2P, ZYMB2P)
 
C
c  dmc 24 Dec 1996:  this version does not depend on TRCOM
c  ** updown symmetric SPLCRD **
c
c  call from MOMRY/MOMRZ; caller does arg error checking!
c
c  input args:
c
c   ZXI99 -- "xi" (radial coordinate) **** must be in range ****
c     ==> interpolate moments set to "xi"
c   IJACI -- IJAC flag as input (see below)
c
c   LCENTR -- 1 + index offset in XIBLO, RMB0, RMB2, YMB2, rel. to.
c      RMB0S,RMB2S,YMB2S
c      (for TRANSP compatibility; set to 1 if all arrays are on
c      the same indexing grid)
c
c   MIMOM -- max no. of moments (array dimension)
c   NMOM  -- actual no. of moments in use
c   MTBL --- max no. of radial points including "bloat" (array dimension)
c
c   XIBLO(..) -- "xi" grid points, including bloat region beyond plasma bdy
c
c   RMB0(..) -- 0'th R moment array
c   RMB0S(..) -- 0'th R moment spline coeffs
c   RMB2(..) -- higher R moments array
c   RMB2S(..) -- higher R moment spline coeffs
c   YMB2(..) -- higher Y moments array
c   YMB2S(..) -- higher Y moment spline coeffs
c
c  output args:
c   ZRMB0X -- interpolated 0'th R moment
c   ZRMB2X(..) -- interpolated higher R moments
c   ZYMB2X(..) -- interpolated higher Y moments
c  (if IJAC set):
c   ZRMB0P -- 0'th R moment dR0/dxi
c   ZRMB2P(..) -- higher R moments derivatives w.r.t. xi
c   ZYMB2P(..) -- higher Y moments derivatives w.r.t. xi
c
C	EVALUATE SPLINE COEFFICIENTS DIRECTLY ... GIVEN A RADIAL
C	COORD "ZXI99", RETURN THE MHD FOURIER COEFFICIENTS (AND THEIR
C	DERIVATIVES) DESCRIBING THE PLASMA SURFACE ON THAT COORD.
C
C  MOD DMC APRIL 1990 -- LINEAR INTERPOLATION OPTION
C
C  IF IJAC=0 OR IJAC=1 USE THE SPLINES
C  IF IJAC=2 OR IJAC=3 DO LINEAR INTERPOLATION INSTEAD
C
C  IF IJAC=1 OR IJAC=3 EVALUATE DATA FOR THE JACOBIAN
C
C   input arrays:
 
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER ijaci,lcentr,mimom,nmom,mtbl,ijac,is,isl1,j,isp1
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zrmb0x,zrmb2x,zymb2x,zrmb0p,zrmb2p,zymb2p,zxi99,zdxii
      REAL*8 zdx99,zdx99i
!============
        REAL*8 RMB0(MTBL),RMB0S(MTBL,3),XIBLO(MTBL)
        REAL*8 RMB2(MTBL,MIMOM),RMB2S(MTBL,MIMOM,3)
        REAL*8 YMB2(MTBL,MIMOM),YMB2S(MTBL,MIMOM,3)
 
C   output arrays:
 
        DIMENSION ZRMB2X(MIMOM), ZYMB2X(MIMOM)
        DIMENSION ZRMB2P(MIMOM), ZYMB2P(MIMOM)
 
C
C-------------------------------------
C
        IJAC=abs(IJACI)
C
        ZDXII= 1.0D0/ (XIBLO(LCENTR+1)-XIBLO(LCENTR))
        IS = ZXI99 * ZDXII + 1
C
 20     continue
        ISL1=IS+LCENTR-1
        if(zxi99.le.xiblo(isl1)) then
           is=is-1
           go to 20
        endif
        if(zxi99.gt.xiblo(isl1+1)) then
           is=is+1
           go to 20
        endif
C
        ZDX99 = ZXI99 - XIBLO(ISL1)
C
C  EVALUATE MOMENTS INTERPOLATION
C
        IF(IJAC.LE.1) THEN
 
C	EVALUATE SPLINE
C
C	  R0(X)
C
           ZRMB0X = RMB0(ISL1) + ZDX99*(RMB0S(IS,1) +
     >			    ZDX99*(RMB0S(IS,2) +
     >			    ZDX99*RMB0S(IS,3)))
 
 
           DO 100 J = 1, NMOM
 
C	      RM(X)
            ZRMB2X(J) = RMB2(ISL1,J) + ZDX99*(RMB2S(IS,J,1) +
     >			   ZDX99*(RMB2S(IS,J,2) + ZDX99*RMB2S(IS,J,3)))
C	      YM(X)
            ZYMB2X(J) = YMB2(ISL1,J) + ZDX99*(YMB2S(IS,J,1) +
     >			   ZDX99*(YMB2S(IS,J,2) + ZDX99*YMB2S(IS,J,3)))
 
 100       CONTINUE
 
        ELSE
C
C  EVALUATE BY LINEAR INTERPOLATION
C
           ISP1=ISL1+1
           ZDX99I=ZDX99*ZDXII
C
           ZRMB0X=RMB0(ISL1)+ZDX99I*(RMB0(ISP1)-RMB0(ISL1))
C
           DO 110 J = 1, NMOM
 
            ZRMB2X(J)= RMB2(ISL1,J)+ZDX99I*(RMB2(ISP1,J)-RMB2(ISL1,J))
            ZYMB2X(J)= YMB2(ISL1,J)+ZDX99I*(YMB2(ISP1,J)-YMB2(ISL1,J))
 
 110       CONTINUE
 
        ENDIF
 
C..................................
 
      IF (IJAC.EQ.0) RETURN
      IF (IJAC.EQ.2) RETURN
 
C..................................
 
C
C  EVALUATE DERIVATIVES
C
 
        IF(IJAC.LE.1) THEN
C
C  SPLINES
C
C	  D[R0(X)]/DX
           ZRMB0P = RMB0S(IS,1) + ZDX99*(2.D0*RMB0S(IS,2) +
     >		     3.D0*ZDX99*RMB0S(IS,3))
 
           DO 200 J = 1, NMOM
 
C	  D[RM(X)]/DX
            ZRMB2P(J) = RMB2S(IS,J,1) + ZDX99*(2.D0*RMB2S(IS,J,2) +
     >			   3.D0*ZDX99*RMB2S(IS,J,3))
C	  D[YM(X)]/DX
            ZYMB2P(J) = YMB2S(IS,J,1) + ZDX99*(2.D0*YMB2S(IS,J,2) +
     >			   3.D0*ZDX99*YMB2S(IS,J,3))
 
 200       CONTINUE
 
        ELSE
C
C  LINEAR INTERPOLATION
C
           ZRMB0P=(RMB0(ISP1)-RMB0(ISL1))*ZDXII
C
           DO 210 J = 1, NMOM
 
            ZRMB2P(J) = (RMB2(ISP1,J)-RMB2(ISL1,J))*ZDXII
            ZYMB2P(J) = (YMB2(ISP1,J)-YMB2(ISL1,J))*ZDXII
 
 210       CONTINUE
 
        ENDIF
C
        RETURN
        END
C******************** END FILE SPLCRD.FOR ; GROUP TKSERVE **************
! 11Jan2003 fgtok -s r8_precision.sub misc.sub "r8con.csh conversion"

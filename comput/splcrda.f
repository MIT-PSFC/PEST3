C******************** START FILE SPLCRD.FOR ; GROUP TKSERVE ************
C.............................................................
 
        SUBROUTINE SPLCRDA (ZXI99, IJACI,
     >                     LCENTR, MIMOM, NMOM, MTBL, XIBLO,
     >                     RMCX,RMCXS,YMCX,YMCXS,
     >			   ZRMCX1,ZRMCX2,ZYMCX1,ZYMCX2,
     >		           ZRMCXP1,ZRMCXP2,ZYMCXP1,ZYMCXP2)
 
C
C  B.BALET - JET  MAR 1994 : MODIFIED VERSION OF SPLCRD.FOR TO SUPPORT
C                            UP-DOWN ASYMMETRIC CASE; SHOULD BE CALLED
C                            ONLY IF LEVGEO = 6 or 7 or 8 or 9 (NLSYM = F)
c
c  dmc 24 Dec 1996:  this version does not depend on TRCOM:  all args passed.
c  ** updown asymmetric SPLCRDA **
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
c   RMCX(..) -- R moments array
c   RMCXS(..) -- R moments spline coeffs
c   YMCX(..) -- Y moments array
c   YMCXS(..) -- Y moments spline coeffs
c
c  output args:
c   ZRMCX1(..) -- R cos moments interpolated
c   ZRMCX2(..) -- R sin moments interpolated
c   ZYMCX1(..) -- Y cos moments interpolated
c   ZYMCX2(..) -- Y sin moments interpolated
c  (if IJAC set):  derivatives w.r.t. xi:
c   ZRMCXP1(..) -- R cos moments derivative
c   ZRMCXP2(..) -- R sin moments derivative
c   ZYMCXP1(..) -- Y cos moments derivative
c   ZYMCXP2(..) -- Y sin moments derivative
c
C
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
 
        REAL XIBLO(MTBL)
        REAL RMCX(MTBL,0:MIMOM,2),RMCXS(MTBL,0:MIMOM,2,3)
        REAL YMCX(MTBL,0:MIMOM,2),YMCXS(MTBL,0:MIMOM,2,3)
 
C   output arrays:
 
        DIMENSION ZRMCX1(0:MIMOM), ZRMCX2(0:MIMOM)
        DIMENSION ZYMCX1(0:MIMOM), ZYMCX2(0:MIMOM)
        DIMENSION ZRMCXP1(0:MIMOM), ZRMCXP2(0:MIMOM)
        DIMENSION ZYMCXP1(0:MIMOM), ZYMCXP2(0:MIMOM)
 
C
C-------------------------------------
C
        IJAC=abs(IJACI)
C
        ZDXII= 1.0 / (XIBLO(LCENTR+1)-XIBLO(LCENTR))
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
        IF(IJAC.LE.1) THEN
 
C	EVALUATE SPLINE
C
           DO 100 J = 0, NMOM
 
C	      RM(X) :
C
            ZRMCX1(J) = RMCX(ISL1,J,1) + ZDX99*(RMCXS(IS,J,1,1) +
     >		       ZDX99*(RMCXS(IS,J,1,2) + ZDX99*RMCXS(IS,J,1,3)))
              ZRMCX2(J) = RMCX(ISL1,J,2) + ZDX99*(RMCXS(IS,J,2,1) +
     >                 ZDX99*(RMCXS(IS,J,2,2) + ZDX99*RMCXS(IS,J,2,3)))
C
C	      YM(X) :
C
              ZYMCX1(J) = YMCX(ISL1,J,1) + ZDX99*(YMCXS(IS,J,1,1) +
     >                 ZDX99*(YMCXS(IS,J,1,2) + ZDX99*YMCXS(IS,J,1,3)))
              ZYMCX2(J) = YMCX(ISL1,J,2) + ZDX99*(YMCXS(IS,J,2,1) +
     >                 ZDX99*(YMCXS(IS,J,2,2) + ZDX99*YMCXS(IS,J,2,3)))
 
 100       CONTINUE
 
        ELSE
C
C  EVALUATE BY LINEAR INTERPOLATION
C
           ISP1=ISL1+1
           ZDX99I=ZDX99*ZDXII
C
           DO 110 J = 0, NMOM
 
              ZRMCX1(J)= RMCX(ISL1,J,1) +
     >           ZDX99I*(RMCX(ISP1,J,1)-RMCX(ISL1,J,1))
              ZRMCX2(J)= RMCX(ISL1,J,2) +
     >           ZDX99I*(RMCX(ISP1,J,2)-RMCX(ISL1,J,2))
              ZYMCX1(J)= YMCX(ISL1,J,1) +
     >           ZDX99I*(YMCX(ISP1,J,1)-YMCX(ISL1,J,1))
              ZYMCX2(J)= YMCX(ISL1,J,2) +
     >           ZDX99I*(YMCX(ISP1,J,2)-YMCX(ISL1,J,2))
 
 110       CONTINUE
 
        ENDIF
 
C..................................
 
        IF (IJAC.EQ.0) RETURN
        IF (IJAC.EQ.2) RETURN
 
C..................................
 
        IF(IJAC.LE.1) THEN
C
C  SPLINES
C
           DO 200 J = 0, NMOM
 
C     D[RM(X)]/DX
              ZRMCXP1(J) = RMCXS(IS,J,1,1) + ZDX99*(2.*RMCXS(IS,J,1,2) +
     >           3.*ZDX99*RMCXS(IS,J,1,3))
              ZRMCXP2(J) = RMCXS(IS,J,2,1) + ZDX99*(2.*RMCXS(IS,J,2,2) +
     >           3.*ZDX99*RMCXS(IS,J,2,3))
 
C     D[YM(X)]/DX
              ZYMCXP1(J) = YMCXS(IS,J,1,1) + ZDX99*(2.*YMCXS(IS,J,1,2) +
     >           3.*ZDX99*YMCXS(IS,J,1,3))
              ZYMCXP2(J) = YMCXS(IS,J,2,1) + ZDX99*(2.*YMCXS(IS,J,2,2) +
     >           3.*ZDX99*YMCXS(IS,J,2,3))
 
 200       CONTINUE
 
        ELSE
C
C  LINEAR INTERPOLATION
C
           DO 210 J = 0, NMOM
 
              ZRMCXP1(J) = (RMCX(ISP1,J,1)-RMCX(ISL1,J,1))*ZDXII
              ZRMCXP2(J) = (RMCX(ISP1,J,2)-RMCX(ISL1,J,2))*ZDXII
              ZYMCXP1(J) = (YMCX(ISP1,J,1)-YMCX(ISL1,J,1))*ZDXII
              ZYMCXP2(J) = (YMCX(ISP1,J,2)-YMCX(ISL1,J,2))*ZDXII
 
 210       CONTINUE
 
        ENDIF
C
        RETURN
        END
C******************** END FILE SPLCRD.FOR ; GROUP TKSERVE **************

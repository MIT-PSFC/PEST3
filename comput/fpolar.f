C******************** START FILE FPOLAR.FOR ; GROUP TKBNDRY ******************
C.....................................................................
C  REV DMC -- REMOVED TRANSP COMMON
C   SO I CAN USE THIS IN RPLOT
C
C  REV BB / MAY 94 -- USE ATAN2 INSTEAD OF ATAN TO EVALUATE FPOLAR
C
C   output in range [0,twopi] instead of range [-pi,pi]
C   atan2 and fpolar are the same in the upper half plane; fpolar
C   is twopi greater in the lower half plane.
C
      REAL FUNCTION FPOLAR ( ZR77, ZZ77 )
C
C	CALC. THE POLAR ANGLE DEFINED BY TAN(THETA)=ZZ77/ZR77
C
      DATA TWOPI/6.2831853071795862E+00/
      DATA PI/3.1415926535897931E+00/
C
C--------------------------------------------
C
      if (zz77 .ge. 0.0) then
C
C  upper half plane  --  including z=0 line
C
       IF (ZR77 .NE. 0.) THEN
C
          FPOLAR = ATAN2 ( ZZ77 , ZR77 )
C
       ELSE
C
C  on the vertical axis, at or above the z=0 line
C
        FPOLAR = PI*0.5
C
       END IF
C
      else
C
C  lower half plane  --  excluding z=0 line
C
          FPOLAR = ATAN2 ( ZZ77 , ZR77 ) + TWOPI
C
      endif
C
      RETURN
      END
C******************** END FILE FPOLAR.FOR ; GROUP TKBNDRY ******************

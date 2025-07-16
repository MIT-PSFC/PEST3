C******************** START FILE LOKJAC.FOR ; GROUP TKBLOAT ******************
 
C-- END FILE #BLOAT# ************************************
CFORMFEEDC-- START FILE #LOKJAC# ************************************
       LOGICAL FUNCTION LOKJAC(R0,DR0DXI,RJ,DRJDXI,YJ,DYJDXI,NMOM)
C
C
C       BOB MCCANN, 19-APR-85
C
C
C       --------
C       COMMENTS
C       --------
C
C
C	THIS CHECKS TO SEE IF THE JACOBIAN D(R<Y)/D(XI,XI) CHANGES SIGN
C    AROUND THE GIVEN FLUX SURFACE.  LOKJAC=.T. IF THE SIGN IS UNCHANGED,
C    AND LOKJAC=.F. IF THE SIGN CHANGES
C    (THE TRANSFORMATION WENT SINGULAR).
C
C
C       -------------------
C       GLOBAL DECLARATIONS
C       -------------------
C
C
C***
C
C
C       ---------------
C       INPUT VARIABLES
C       ---------------
C
C
      INTEGER	NMOM		!# OF MOMENTS OF R, Y
C
C    THE FOLLOWING REFER TO A GIVEN FLUX SURFACE:
C
      REAL	R0		!ZEROTH MOMENT
      REAL	DR0DXI		!D(R0)/DXI
      REAL	RJ(NMOM)	!J=1,NMOM MOMENTS OF R
      REAL	DRJDXI(NMOM)	!D(RJ)/DXI
C
      REAL	YJ(NMOM)	!K=1, NMOM MOMENTS OF Y
      REAL	DYJDXI(NMOM)	!D(YJ)/DXI
C
C
C       ----------------
C       OUTPUT VARIABLES
C       ----------------
C
C
C	LOGICAL LOKJAC	!=.T. IF THE TRANSFORMATION IS OK
C			!=.F. IF THE COORDINATES ARE SINGULAR
C
C
C       ------------------
C       LOCAL DECLARATIONS
C       ------------------
C
C
      PARAMETER	(NTH = 127)	!# ANGLES TO USE IN FINDING EXTREMA
C
      REAL	ZC(NTH)		!COS(J*TH) ASSUMING NMOM.LE.NTH
      REAL	ZS(NTH)
C
      REAL	ZR77		!R(TH)
      REAL	ZY77		!Y(TH)
      REAL	ZDRDTH		!DR/DTH
      REAL	ZDYDTH		!DY/DTH
      REAL	ZDRDXI		!DR/DXI
      REAL	ZDYDXI		!DY/DXI
      REAL	ZJ		!JACOBIAN
C
C
C       ---------------
C       DATA STATEMENTS
C       ---------------
C
C
C***
C
C
C       -------------------
C       STATEMENT FUNCTIONS
C       -------------------
C
C
C***
C
C
C=========================================
C
C
C       ----------------------
C       0.1     INITIALIZATION
C       ----------------------
C
C
      ZPI = 4.0*ATAN(1.0)
C
C    LOOP OVER SOME ANGLES FROM 0 TO PI:
C
      DO 69000 I = 1, NTH
C
        ZTH = ZPI*FLOAT(I-1)/FLOAT(NTH-1)
C
        DO 100 J = 1, NMOM
          ZC(J) = COS(J*ZTH)
          ZS(J) = SIN(J*ZTH)
100     CONTINUE
C
C
C	----------------------------------------------
C	1.0	COMPUTE VARIOUS DERIVATIVES OF R AND Y
C	----------------------------------------------
C
C
10000 CONTINUE
C
C    HANDLE ZEROTH MOMENT AND INITIALIZE THE MOMENTS LOOP
C
      ZR77 = R0
      ZDRDXI = DR0DXI
      ZDRDTH = 0.0
C
      ZY77 = 0.0
      ZDYDTH = 0.0
      ZDYDXI = 0.0
C
C    ADD UP THE CONTRIBUTION FROM EACH MOMENT
C
      DO 10100 J = 1, NMOM
C
        ZR77 = ZR77 + RJ(J)*ZC(J)
C
        ZDRDXI = ZDRDXI + DRJDXI(J)*ZC(J)
        ZDRDTH = ZDRDTH - RJ(J)*FLOAT(J)*ZS(J)
C
        ZY77 = ZY77 + YJ(J)*ZS(J)
C
        ZDYDXI = ZDYDXI + DYJDXI(J)*ZS(J)
        ZDYDTH = ZDYDTH + YJ(J)*FLOAT(J)*ZC(J)
C
10100 CONTINUE
C
C
C       -------------------------------------------------
C       2.0     COMPUTE THE JACOBIAN OF THIS FLUX SURFACE
C       -------------------------------------------------
C
C
20000   CONTINUE
C
      ZJ = ZDRDXI*ZDYDTH - ZDRDTH*ZDYDXI	!JACOBIAN
C
C
C	--------------------------------------------------
C	3.0	FIND THE MAXIMA AND MINIMA OF THE JACOBIAN
C	--------------------------------------------------
C
C
30000 CONTINUE
C
      IF(I.EQ.1)	THEN	!INITIALIZE MAX AND MIN
C
        ZJMIN = ZJ
        THJMIN = ZTH
        ZJMAX = ZJ
        THJMAX = ZTH
C
      ELSE			!FIND THE TRUE MIN AND MAX
C
        IF(ZJ.LT.ZJMIN)	THEN
          ZJMIN = ZJ
          THJMIN = ZTH
        ELSEIF(ZJ.GT.ZJMAX)	THEN
          ZJMAX = ZJ
          THJMAX = ZTH
        ENDIF
C
      ENDIF
C
69000 CONTINUE		!END OF ANGLE LOOP
C
C
C       ------------------------
C       7.0     CLEANUP & RETURN
C       ------------------------
C
C
70000   CONTINUE
C
C
      IF(ZJMIN*ZJMAX.GT.0.0)	THEN	!THEY HAVE THE SAME SIGN
          LOKJAC = .TRUE.
      ELSE
          LOKJAC = .FALSE.
      ENDIF
C
C
79000   CONTINUE        !ERROR RETURNS HERE
C
        RETURN
C
C
C       ----------------------
C       8.0     ERROR HANDLING
C       ----------------------
C
C
C***
C
C
        END
C******************** END FILE LOKJAC.FOR ; GROUP TKBLOAT ******************

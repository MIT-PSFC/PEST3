       LOGICAL FUNCTION LOKJACA
     >                  (RJC,RJS,DRCDXI,DRSDXI,YJC,YJS,DYCDXI,DYSDXI,NM)
C
C
C       B. BALET, JAN 1994 : MODIFICATIONS TO THE FUNCTION LOKJAC WRITTEN BY
C                            BOB MCCANN FOR THE UP-DOWN ASYMMETRIC CASE
C
C       --------
C       COMMENTS
C       --------
C
C
C	THIS CHECKS TO SEE IF THE JACOBIAN D(R<Y)/D(XI,XI) CHANGES SIGN
C    AROUND THE GIVEN FLUX SURFACE.  LOKJACA=.T. IF THE SIGN IS UNCHANGED,
C    AND LOKJACA=.F. IF THE SIGN CHANGES
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
      INTEGER	NM  		!# OF MOMENTS OF R, Y
C
C    THE FOLLOWING REFER TO A GIVEN FLUX SURFACE:
C
      REAL	RJC(0:NM) 	!J=0,NM MOMENTS OF R : COS TERMS
        REAL	RJS(0:NM) 	!J=0,NM MOMENTS OF R : SIN TERMS
      REAL	DRCDXI(0:NM)	!D(RJC)/DXI
      REAL	DRSDXI(0:NM)	!D(RJS)/DXI
C
        REAL    YJC(0:NM)         !J=0,NM MOMENTS OF Y : COS TERMS
        REAL    YJS(0:NM)         !J=0,NM MOMENTS OF Y : SIN TERMS
        REAL    DYCDXI(0:NM)      !D(YJC)/DXI
        REAL    DYSDXI(0:NM)      !D(YJS)/DXI
C
C
C       ----------------
C       OUTPUT VARIABLES
C       ----------------
C
C
C	LOGICAL LOKJACA	!=.T. IF THE TRANSFORMATION IS OK
C			!=.F. IF THE COORDINATES ARE SINGULAR
C
C
C       ------------------
C       LOCAL DECLARATIONS
C       ------------------
C
C
      PARAMETER	(NTH = 254)	!# ANGLES TO USE IN FINDING EXTREMA
C
      REAL	ZC(0:NTH)		!COS(J*TH) ASSUMING NM.LE.NTH
      REAL	ZS(0:NTH)         !SIN(J*TH) ASSUMING NM.LE.NTH
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
      ZTWOPI = 8.0*ATAN(1.0)
C
C    LOOP OVER SOME ANGLES FROM 0 TO TWO*PI:
C
      DO 69000 I = 1, NTH
C
        ZTH = ZTWOPI*FLOAT(I-1)/FLOAT(NTH-1)
C
        DO 100 J = 0, NM
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
C    INITIALIZE THE MOMENTS LOOP
C
      ZR77 = 0.0
      ZDRDXI = 0.0
      ZDRDTH = 0.0
C
      ZY77 = 0.0
      ZDYDTH = 0.0
      ZDYDXI = 0.0
C
C    ADD UP THE CONTRIBUTION FROM EACH MOMENT
C
      DO 10100 J = 0, NM
C
        ZR77 = ZR77 + RJC(J)*ZC(J) + RJS(J)*ZS(J)
C
        ZDRDXI = ZDRDXI + DRCDXI(J)*ZC(J) + DRSDXI(J)*ZS(J)
        ZDRDTH =
     >    ZDRDTH - RJC(J)*FLOAT(J)*ZS(J) + RJS(J)*FLOAT(J)*ZC(J)
C
          ZY77 = ZY77 + YJC(J)*ZC(J) + YJS(J)*ZS(J)
C
          ZDYDXI = ZDYDXI + DYCDXI(J)*ZC(J) + DYSDXI(J)*ZS(J)
          ZDYDTH =
     >    ZDYDTH - YJC(J)*FLOAT(J)*ZS(J) + YJS(J)*FLOAT(J)*ZC(J)
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
          LOKJACA = .TRUE.
      ELSE
          LOKJACA = .FALSE.
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

C******************** START FILE ECRAD.FOR ; GROUP ECRAD ******************
C---------------------------------------------------------------------
C  ECRAD
C
C  REAL FUNCTION, DEFINE ECE RADIUS FROM ECE FREQUENCY
C   HARMONIC AND R*TOROIDAL FIELD STRENGTH PASSED THRU COMMON
C
      REAL FUNCTION ECRAD(FREQ)
C
C--------------------------------------------------------------------
C
C  COMMON BLOCK ECRADI-- TO CONTAIN ECE HARMONIC AND R*TOROIDAL FIELD
C
      COMMON/ECRADI/ RBT,JHARM
C
C  DMC JAN 1984-- FUDGE FACTOR 28.0 PUTS ECRAD IN CM
C
C       RBT IS IN TESLA*CM
C       FREQ IS IN GHZ
C       HARMONIC JHARM = 1, 2, 3, ...
C
C  FUDGE FACTOR COMES FROM STANDARD FORMULA
C
C                 (FIELD STRENGTH)*(ELECTRON CHARGE)
C      OMEGA(CE)= ----------------------------------
C                  (ELECTRON MASS)*(SPEED OF LIGHT)
C
C  MULTIPLY THRU BY MAJOR RADIUS/OMEGA(CE), AND NOTE THAT
C  FREQUENCY = OMEGA/(2*PI)
C
C--------------------------------------------------------------------
C
      ECRAD=28.0*RBT*JHARM/FREQ
C
      RETURN
C
      END
C******************** END FILE ECRAD.FOR ; GROUP ECRAD ******************

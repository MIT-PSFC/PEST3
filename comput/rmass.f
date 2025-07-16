      SUBROUTINE RMASS(NREACT,KK,AM1,AM2,AM3,AM4)
C
C...5/10/76
C...This subroutine was written by Harry H. Towner of the Princeton
C...Plasma Physics lab.  This routine will return the masses for
C...a given reaction.  The reaction is assumed to be in the form of:
C		AM1(AM2,AM3)AM4
C...PARAMETERS:
C...NREACT - Specifies the reaction.
C		1)	T(D,N)HE4	2)   HE3(D,P)HE4
C		3)	D(D,P)T		4)   D(D,N)HE3
C		5)	T(T,2N)HE4 	6)   T(HE3,X)Y
C...KK - If KK=1 then the above reaction is reversed.
C...AM1-4 - The mass for the given reaction.  (keV/C**2)
C...Note that indeterminate masses are given a value of -1.e00.
C...All masses have the units of kev/c**2.
C
      DIMENSION A1(6),A2(6),A3(6),A4(6)
      DATA A1/2.80888259E+6,2.808352957E+6,2*1.875587331E+6,
     >2*2.80888259E+6/
      DATA A2/4*1.875587331E+6,2.80888259E+6,2.808352957E+6/
      DATA A3/9.395526709E+5,2*9.382592233E+5,9.395526709E+5,2*-1.E00/
      DATA A4/2*3.727327439E+6,2.80888259E+6,2.808352957E+6,
     >3.727327439E+6,-1.E00/
C
C************************************************************************
C
      IF (KK .EQ. 0) THEN
        AM1=A1(NREACT)
        AM2=A2(NREACT)
      ELSE IF (KK .EQ. 1) THEN
        AM1=A2(NREACT)
        AM2=A1(NREACT)
      ELSE
        WRITE(6,10)
10      FORMAT(T2,I2,' is an illegal value for KK!')
        call bad_exit
      END IF
C
      AM3=A3(NREACT)
      AM4=A4(NREACT)
      RETURN
      END

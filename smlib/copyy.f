C******************** START FILE COPYY.FOR ; GROUP GSMOO1 ******************
C------------------------------------------------------------
C  COPYY
C
C  COPY DATA FROM ARRAY TO ARRAY, FACTOR IN ZMUL
C
      SUBROUTINE COPYY(Y1,ZMUL,Y2,N)
      DIMENSION Y1(N),Y2(N)
C
      DO 10 I=1,N
      Y2(I)=ZMUL*Y1(I)
 10   CONTINUE
      RETURN
      END
C******************** END FILE COPYY.FOR ; GROUP GSMOO1 ******************

C******************** START FILE FILTER.FOR ; GROUP FILTR6 ******************
C------------------------------------------------------
C  FILTER
C
C  USE LOCAL AVERAGING FUNCTION OF RADIUS DELTA
C
C  TO SMOOTH DATA CURVE.
C
C-----------------
      SUBROUTINE FILTER(X,Y0,Y,N,DELTA,XEND1,XEND2)
C
      DIMENSION X(N),Y0(N),Y(N),DELTA(N)
C
C  START OF EXECUTABLE CODE --
C
C
      DO 100 J=1,N
C
      JLP=J
      Y(J)=FILFN6(X,Y0,N,JLP,DELTA,XEND1,XEND2,1)
C
 100  CONTINUE
C
      RETURN
C
 
      	END
C******************** END FILE FILTER.FOR ; GROUP FILTR6 ******************

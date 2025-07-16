C******************** START FILE FILTER.FOR ; GROUP FILTR6 ******************
C------------------------------------------------------
C  FILTER
C
C  USE LOCAL AVERAGING FUNCTION OF RADIUS DELTA
C
C  TO SMOOTH DATA CURVE.
C
C-----------------
      SUBROUTINE r8filter(X,Y0,Y,N,DELTA,XEND1,XEND2)
C
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n,j,jlp
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 y0,y,delta,xend1,xend2,x,r8filfn6
C============
      DIMENSION X(N),Y0(N),Y(N),DELTA(N)
C
C  START OF EXECUTABLE CODE --
C
C
      DO 100 J=1,N
C
      JLP=J
      Y(J)=r8filfn6(X,Y0,N,JLP,DELTA,XEND1,XEND2,1)
C
 100  CONTINUE
C
      RETURN
C
 
      	END
C******************** END FILE FILTER.FOR ; GROUP FILTR6 ******************
! 03Dec1999 fgtok -s r8_precision.sub r8smlib.sub "r8con.csh conversion"

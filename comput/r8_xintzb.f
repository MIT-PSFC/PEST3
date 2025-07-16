C******************** START FILE R8_XINTZB.FOR ; GROUP IXCALC ****************
C-------------------------------------------------------------
C  R8_XINTZB
C
C  INTERPOLATE ZONE-CENTERED VARIABLE TO ZONE BOUNDARY VALUES
C
C  ONE SLIGHT EXTRAPOLATION AT OUTER EDGE
C
      SUBROUTINE R8_XINTZB(FZC,FZB,N)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER n,inm1,i,ip1
!============
      REAL*8 FZC(N),FZB(N)
C
      INM1=N-1
      DO 10 I=1,INM1
         IP1=I+1
         FZB(I)=0.5D0*(FZC(I)+FZC(IP1))
 10   CONTINUE
 
C       Linear EXTRAPOLATION AT EDGE
      FZB(N)=FZC(N)+0.5D0*(FZC(N)-FZC(INM1))
C
      RETURN
      END
C******************** END FILE R8_XINTZB.FOR ; GROUP IXCALC ******************

C******************** START FILE XINTZB.FOR ; GROUP IXCALC ******************
C-------------------------------------------------------------
C  XINTZB
C
C  INTERPOLATE ZONE-CENTERED VARIABLE TO ZONE BOUNDARY VALUES
C
C  ONE SLIGHT EXTRAPOLATION AT OUTER EDGE
C
      SUBROUTINE XINTZB(FZC,FZB,N)
C
      REAL FZC(N),FZB(N)
C
      INM1=N-1
      DO 10 I=1,INM1
         IP1=I+1
         FZB(I)=0.5*(FZC(I)+FZC(IP1))
 10   CONTINUE
 
C       Linear EXTRAPOLATION AT EDGE
      FZB(N)=FZC(N)+0.5*(FZC(N)-FZC(INM1))
C
      RETURN
      END
C******************** END FILE XINTZB.FOR ; GROUP IXCALC ******************

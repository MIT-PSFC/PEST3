C******************** START FILE XINTZ0.FOR ; GROUP IXCALC ******************
C-----------------------------------------------------------------
C  XINTZ0
C
C  INTERPOLATE ZONE BDY FCN TO VALUES AT ZONE CENTER.
C   Assume value is 0.0 at very center FZB(0) (if it existed)
C     tbt 6/1
 
 
      SUBROUTINE XINTZ0(FZB,FZC,N)
 
 
C	Updates:
C	tbt 06/01/95  At GA. Copied from XintzC
 
C	-----------------
      REAL FZC(N),FZB(N)
      Real Zfactor
C	-----------------
 
C       Interpolation AT CENTER - assume = 0
        FZB0 = 0.0
 
      FZC(1)= 0.5*(FZB0+FZB(1))
C
      DO 20 I=2,N
      IM1=I-1
      FZC(I)=0.5*(FZB(IM1)+FZB(I))
 20   CONTINUE
      RETURN
      END
C******************** END FILE XINTZ0.FOR ; GROUP IXCALC ******************

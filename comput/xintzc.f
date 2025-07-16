C******************** START FILE XINTZC.FOR ; GROUP IXCALC ******************
C-----------------------------------------------------------------
C  XINTZC
C
C  INTERPOLATE ZONE BDY FCN TO VALUES AT ZONE CENTER.
C   ONE SLIGHT EXTRAPOLATION AT CENTER-- ASSUME APPROACHING A ZERO
C   RATE OF CHANGE THERE
 
 
      SUBROUTINE XINTZC(FZB,FZC,N)
 
 
C	Updates:
C	tbt 03/08/94 - Changed Zfactor to .125 from .3333333
 
C	-----------------
      REAL FZC(N),FZB(N)
      Real Zfactor
C	-----------------
 
C       EXTRAPOLATION AT CENTER
      Zfactor = .125          ! Was .33333  - for parabolic fit.
 
      FZC(1)=FZB(1)- Zfactor*(FZB(2)-FZB(1))
C
      DO 20 I=2,N
      IM1=I-1
      FZC(I)=0.5*(FZB(IM1)+FZB(I))
 20   CONTINUE
      RETURN
      END
C******************** END FILE XINTZC.FOR ; GROUP IXCALC ******************

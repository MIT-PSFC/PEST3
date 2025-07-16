      FUNCTION SIGE(E)
C...3/8/77
      COMMON /CSIGMA/ A1,A2,A3,A4,A5,FE
C...This function is used to calculate the fusion cross section
C...after a call to sigma has been made.  The functional form of
C...the fusion cross sections is given in the paper by B. H. Duane
C...in BNWL-1685.  The value of the cross section is returned by
C...SIGE and has the units of m**2.
C
C---------------------------
C  MOD DMC 11 SEP 1991 -- SUPPORT BOSCH NEW BOSCH FITS.
C  SEE COMMENTS IN SUBROUTINE SIGMA (SIGMA.FOR) FOR FULL EXPLANATION
C
C NEW COMMON BLOCK:
C
      REAL BG,EFAC,ELOW,EMID,EHI,SLOPE
      REAL ALOW(5),BLOW(4),AHI(5),BHI(4)
      LOGICAL LBOSCH
      COMMON /CSIGMA2/ BG,EFAC,ELOW,EMID,EHI,
     >    ALOW,BLOW,AHI,BHI,SLOPE,LBOSCH
C
C---------------------------
C
C...PARAMETER
C...E - THE ENERGY OF THE INCIDENT PARTICLE.  (KEV)
C...I NOW DEFINE THE FUNCTION USED TO EVALUATE SIGMA.
      FUN(B,A2,A3,A4,A5,E)=B*(A2/(1.E00+(A3*E-A4)**2)+A5)/(E*(1.E00-B
     >))
C
C  NEW 11 SEP 1991
C  FBOSCH is the analytic form of H.-S. Bosch's cross section fits.
       FBOSCH(E,A1,A2,A3,A4,A5,B1,B2,B3,B4)=
     1 (A1+E*(A2+E*(A3+E*(A4+E*A5))))/(1.+E*(B1+E*(B2+E*(B3+E*B4))))
C
C***********************************************************************
C
      IF(LBOSCH) THEN
C----------------------
C  BEGINNING OF NEW BOSCH FIT CODE -- DMC 11 SEP 1991
C
        SIGE=0.0
        ENERGY=EFAC*E  ! CONVERT TO CTR OF MASS FRAME
C
        EINPUT=MIN(ENERGY,EHI)
        IF(EINPUT.GE.ELOW) THEN
          IF(EINPUT.LE.EMID) THEN
            ZS=FBOSCH(EINPUT,ALOW(1),ALOW(2),ALOW(3),ALOW(4),ALOW(5),
     >                  BLOW(1),BLOW(2),BLOW(3),BLOW(4))
          ELSE
            ZS=FBOSCH(EINPUT,AHI(1),AHI(2),AHI(3),AHI(4),AHI(5),
     >                  BHI(1),BHI(2),BHI(3),BHI(4))
          ENDIF
          IF(ENERGY.GT.EHI) ZS=ZS+SLOPE*(ENERGY-EHI)
          SIGE=1.E-31*(ZS*EXP(-BG/SQRT(ENERGY))/ENERGY)
        ENDIF
C
C  END OF NEW BOSCH FIT CODE
C-----------------------
      ELSE
C-----------------------
C  BEGINNING OF OLD FIT CODE
        IF (E .GT. 1.E-3) THEN
          ENERGY=E*FE*1.E+3
C...Note that the energy has been corrected for a possible reverse
C...reaction.
          ADD=A1/SQRT(ENERGY)
C
          IF (ADD .LT. 87.5) THEN
            B=EXP(-ADD)
            SIGE=FUN(B,A2,A3,A4,A5,ENERGY)
C
            IF (SIGE .GT. 1.E-10) THEN
      	SIGE = SIGE*1.E-28
            ELSE
      	SIGE = 0.
            END IF
C
          ELSE
            SIGE = 0.
          END IF
C
        ELSE
          SIGE = 0.
        END IF
C
C  END OF OLD FIT CODE
C-----------------------
      ENDIF
C
      RETURN
      END

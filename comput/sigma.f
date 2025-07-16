c---------------------------------------------------------
c  dmc May 2001:  real*8 interface to the old "sigma" cross-section routine
c
      real*8 function r8_sigma(n,kk,ipunit,er8)
      implicit NONE
      integer n,kk,ipunit
      real*8 er8
c
      real e,zsig,sigma
c
      e=er8
      zsig=sigma(n,kk,ipunit,e)
      r8_sigma=zsig
      return
      end
 
c---------------------------------------------------------
 
      FUNCTION SIGMA(N,KK,IPUNIT,E)
 
 
      COMMON /CSIGMA/ A1,A2,A3,A4,A5,FE
 
C  10/09/94 TBT
C           Add Lithium reactions 7-11
 
C...3/08/77
C...This  function was written by Harry H. Towner of the
C...Princeton Plasma Physics Lab. This function will initialize
C...parameters in common block csigma and then call function
C...SIGE to calculate the desired fusion cross-section (m**2).
C...if the only parameter which varies is e then SIGE should
C...be used after an initializing call to sigma is made.
C
C-----------------------
C UPGRADE - DMC - 11 SEP 1991
C  CODE IS ADDED TO USE THE NEW BOSCH FITTING PROCEDURE FOR THE
C  REACTIONS 1--4, SUPERSEDING THE OLDER DUANE/JARMIE FITS USED
C  HERETOFORE.  HOWEVER, THE OLD DUANE FITTING CODE IS RETAINED
C  AND MAY BE REVIVED BY CHANGING THE WAY THE LOGICAL VARIABLE
C  LBOSCH IS SET.
C
C THE BOSCH FITS AND FORMULAE ARE GIVEN IN THE PAPER
C
C  "Review of Data and Formulas for Fusion Cross-Sections"
C  by Hans-Stephan Bosch   IPP I/252, Garching, September 1990
C
C The Gamov factor "BG", sqrt(KeV), formula is on p. 2.
C A conversion to center-of-mass frame is required, cf. p. 2.
C The fit coefficients are on p. 29.
C
C Note that for D-T and D-He3 the fitting is split into two
C energy ranges; for simplicity the code supports two energy
C ranges for all reactions
C
C the call to SIGMA sets up COMMON block CSIGMA2 for subsequent
C calls to SIGE to actually evaluate the cross-sections, as before.
C
C NEW COMMON BLOCK:
C
      COMMON /CSIGMA2/ BG,EFAC,ELOW,EMID,EHI,
     >    ALOW,BLOW,AHI,BHI,SLOPE,LBOSCH
      LOGICAL LBOSCH
      REAL BG,EFAC,ELOW,EMID,EHI,SLOPE
      REAL ALOW(5),BLOW(4),AHI(5),BHI(4)
C
C NEW DATA:
C
      REAL BGDATA(4)   ! GAMOV FACTORS
      REAL ERDATA(3,4) ! ENERGY FIT RANGES
      REAL ALDATA(5,4),AHDATA(5,4),BLDATA(4,4),BHDATA(4,4) !COEFFS
C
C  DATA STATEMENTS BELOW...
C
C-----------------------
C
      Parameter (NREACTS=11)
C...PARAMETERS:
C...N - Tells what reaction you want
C	  =1 you get the cross section for the reaction  T(D,N)HE4
C	  =2 you get the cross section for the reaction  HE3(D,P)HE4
C	  =3 you get the cross section for the reaction  D(D,P)T
C	  =4 you get the cross section for the reaction  D(D,N)HE3
C	  =5 you get the cross section for the reaction  T(T,2N)HE4
C	  =6 you get the cross section for the reaction  T(HE3,X)Y
C         =7 Li
C         =8 Li    fudge tbt
C         =9 Li
C         =10Li
C         =11Li
C...KK - If KK=1 than you get the reverse reaction ie if N=1 and
C	 KK=1 you get the cross section for the reaction D(T,N)HE4
C...IPUNIT - the logical unit used for printing.
C...E - The energy of the incident particle.  (KEV)
C
      DIMENSION AA1(NREACTS),AA2(NREACTS),AA3(NREACTS)
        Dimension AA4(NREACTS),AA5(NREACTS),FFE(NREACTS)
 
C       Lithium Fudge tbt (last 5 repeated.
      DATA AA1/1453.E00,2823.E00,1457.7E00,1514.E00,1214.E00,
     >3893.E00,
     >2823.E00,1457.7E00,1514.E00,1214.E00,
     >3893.E00/      ! tbt
 
      DATA AA2/502.E+5,259.E+5,372.E+3,482.E+3,448.E+3,1125.E+4,
     >                   259.E+5,372.E+3,482.E+3,448.E+3,1125.E+4/
      DATA AA3/1368.E-8,398.E-8,436.E-9,308.E-9,102.E-8,0.E00,
     >                    398.E-8,436.E-9,308.E-9,102.E-8,0.E00/
      DATA AA4/1.076E00,1.297E00,1.220E00,1.177E00,2.09E00,0.E00,
     >                    1.297E00,1.220E00,1.177E00,2.09E00,0.E00/
 
      DATA AA5/409.E+3,647.E+3,0.E00,0.E00,0.E00,0.E00,
     >                   647.E+3,0.E00,0.E00,0.E00,0.E00/
 
      DATA FFE/1.497466E00,1.497456E00,1.E00,1.E00,1.E00,1.000007E00,
     >                       1.497456E00,1.E00,1.E00,1.E00,1.000007E00/
 
C-----------------------------
C  DATA STATEMENTS FOR FIT COEFFS FOR NEW BOSCH FITS
C  DMC 11 SEP 1991 -- SEE COMMENTS, ABOVE...
C
      DATA BGDATA/34.3827,68.7508,31.3970,31.3970/
      DATA ERDATA/
     >    0.5, 550.0, 4700.0,
     >    0.3, 900.0, 4800.0,
     >    0.5, 5000.0, 5000.0,
     >    0.5, 4900.0, 4900.0/
      DATA ALDATA/
     >    6.9270E+04, 7.4540E+08, 2.0500E+06, 5.2002E+04, 0.0,
     >    5.7501E+06, 2.5226E+03, 4.5566E+01, 0.0, 0.0,
     >    5.5576E+04, 2.1054E+02,-3.2638E-02, 1.4987E-06, 1.8181E-10,
     >    5.3701E+04, 3.3027E+02,-1.2706E-01, 2.9327E-05,-2.5151E-09/
      DATA BLDATA/
     >    6.3800E+01,-9.9500E-01, 6.9810E-05, 1.728E-04,
     >   -3.1995E-03,-8.5530E-06, 5.9014E-08, 0.0,
     >   0.0, 0.0, 0.0, 0.0,
     >   0.0, 0.0, 0.0, 0.0/
      DATA AHDATA/
     >   -1.4714E+06, 0.0, 0.0, 0.0, 0.0,
     >   -8.3993E+05, 0.0, 0.0, 0.0, 0.0,
     >   0.0, 0.0, 0.0, 0.0, 0.0,
     >   0.0, 0.0, 0.0, 0.0, 0.0/
      DATA BHDATA/
     >   -8.4127E-03, 4.7983E-06,-1.0748E-09,8.5184E-14,
     >   -2.6830E-03, 1.1633E-06,-2.1332E-10,1.4250E-14,
     >   0.0, 0.0, 0.0, 0.0,
     >   0.0, 0.0, 0.0, 0.0/
C
C  NEW 11 SEP 1991
C  FBOSCH is the analytic form of H.-S. Bosch's cross section fits.
       FBOSCH(E,A1,A2,A3,A4,A5,B1,B2,B3,B4)=
     1 (A1+E*(A2+E*(A3+E*(A4+E*A5))))/(1.+E*(B1+E*(B2+E*(B3+E*B4))))
C
C***********************************************************************
C
C
      IF (N .GE. 1 .AND. N .LE. NREACTS) THEN
C
C  DMC 11 SEP 1991 -- DECIDE BTW OLD DUANE FITS AND NEW BOSCH FITS
C
        IF(N.LE.4) THEN
C  SET UP NEW BOSCH FITS
          LBOSCH=.TRUE.
          BG=BGDATA(N)
          CALL RMASS(N,KK,AM1,AM2,AM3,AM4)
          EFAC=AM1/(AM1+AM2)
          ELOW=ERDATA(1,N)
          EMID=ERDATA(2,N)
          EHI=ERDATA(3,N)
          DO I=1,5
            ALOW(I)=ALDATA(I,N)
            AHI(I)=AHDATA(I,N)
          END DO
          DO I=1,4
            BLOW(I)=BLDATA(I,N)
            BHI(I)=BHDATA(I,N)
          END DO
C
C  COMPUTE SLOPE FOR HIGH ENERGY EXTRAPOLATION
C  I HAVE TESTED THESE TO 25 MeV USING TREXE:NUKTEST.EXE
C  DMC 12 SEP 1991
C
          ZE2=EHI
          ZE1=0.99*EHI
          IF(EHI.GT.EMID) THEN
            ZS1=FBOSCH(ZE1,AHI(1),AHI(2),AHI(3),AHI(4),AHI(5),
     >                  BHI(1),BHI(2),BHI(3),BHI(4))
            ZS2=FBOSCH(ZE2,AHI(1),AHI(2),AHI(3),AHI(4),AHI(5),
     >                  BHI(1),BHI(2),BHI(3),BHI(4))
          ELSE
            ZS1=FBOSCH(ZE1,ALOW(1),ALOW(2),ALOW(3),ALOW(4),ALOW(5),
     >                  BLOW(1),BLOW(2),BLOW(3),BLOW(4))
            ZS2=FBOSCH(ZE2,ALOW(1),ALOW(2),ALOW(3),ALOW(4),ALOW(5),
     >                  BLOW(1),BLOW(2),BLOW(3),BLOW(4))
          ENDIF
C
          SLOPE=(ZS2-ZS1)/(ZE2-ZE1)
C  MODIFICATION SEEMS BEST FOR D-T (DMC 12 SEP)
          IF(N.EQ.1) SLOPE=0.15*SLOPE
C
        ELSE
          LBOSCH=.FALSE.
C  SET UP OLD FITS
          FE=1.E00
          IF(KK.EQ.1) FE=1.E00/FFE(N)
          A1=AA1(N)
          A2=AA2(N)
          A3=AA3(N)
          A4=AA4(N)
          A5=AA5(N)
        ENDIF
C
        SIGMA=SIGE(E)
        RETURN
C
      ELSE
C  ERROR CHECK
        WRITE(IPUNIT,10) N
10      FORMAT(1X,'This value of n is not allowed ',I3,
     >	          ' returned from function sigma')
        call bad_exit
      ENDIF
C
      END

      SUBROUTINE GEOTRN(MODE, XIN, EIN, DIN, EOUT, DOUT)
C***********************************************************************
C*****GEOTRN TRANSFORMS ELONGATION AND TRIANGULARITY FROM A MOMENTS TO *
C*****A GEOMETRICAL REPRESENTATION AND VICE VERSA.                     *
C*****REFERENCES:                                                      *
C*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *
C*****CALCULATED PARAMETERS:                                           *
C*****EOUT-OUTPUT ELONGATION.                                          *
C*****DOUT-OUTPUT TRIANGULARITY.                                       *
C*****INPUT PARAMETERS:                                                *
C*****MODE-DESIGNATES DIRECTION OF TRANSFORMATION.                     *
C*****    =1 GEOMETRICAL => MOMENTS.                                   *
C*****    =2 MOMENTS => GEOMETRICAL.                                   *
C*****XIN-INPUT REDUCED MINOR RADIUS                                   *
C*****EIN-INPUT ELONGATION.                                            *
C*****DIN-INPUT TRIANGULARITY.                                         *
C***********************************************************************
      IF (XIN.LE.0.0) GO TO 30
      IF (MODE.EQ.2) GO TO 20
C*****GEOMETRICAL ==> MOMENTS REPRESENTATION.
      DOUT = DIN/4.0
      CTC = 0.0
      DO 10 I=1,10
           CTC = 4.0*DOUT/(SQRT(XIN**2+32.0*DOUT**2)+XIN)
           DOUT = XIN*DIN/(4.0-6.0*CTC**2)
   10 CONTINUE
      EOUT = XIN*EIN/(SQRT(1.0-CTC*CTC)*(XIN+2.0*DOUT*CTC))
      RETURN
C*****MOMENTS ==> GEOMETRICAL REPRESENTATION.
   20 CTC = 4.0*DIN/(SQRT(XIN**2+32.0*DIN**2)+XIN)
      STC = SQRT(1.0-CTC*CTC)
      S2TC = 2.0*STC*CTC
      C2TC = 2.0*CTC*CTC - 1.0
      EOUT = EIN*(XIN*STC+DIN*S2TC)/XIN
      DOUT = DIN*(1.0-3.0*C2TC)/XIN
      RETURN
   30 EOUT = EIN
      DOUT = DIN
      RETURN
      END

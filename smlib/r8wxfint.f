C******************** START FILE WXFINT.FOR ; GROUP FILTR6 ******************
C---------------------------------------------------------------
C   WXFINT
 
C
C  COMPUTES INTEGRAL OF W(T)*F(T)*DT
C   WHERE W(T)=1-ABS(XCEN-T)/DELTA IS THE WEIGHTING FUNCTION
C   AND F(T) IS THE PIECEWISE LINEAR INTERPOLATION FUCNTION
C   OF THE DATA SERIES.
C
C   INTEGRATE FROM X1 TO X2
C   (XF1,YF1),(XF2,YF2) ARE THE ENDPOINTS OF A LINEAR
C   SEGMENT OF F(T);  XCEN SHOULD NOT LIE BETWEEN XF1 AND
C   XF2, X1 AND X2 SHOULD EQUAL OR LIE
C   BETWEEN XF1 AND XF2
C
      REAL*8 FUNCTION r8wxfint(X1,X2,XF1,YF1,XF2,YF2,XCEN,DELTA)
C
C============
C idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 x2,xf1,yf1,xf2,yf2,xcen,delta,x1,zs,xa,xb,xfa,zc0,zdi,zc1,
     > zc2
      integer isign
C============
      ZS=(YF2-YF1)/(XF2-XF1)
      XA=max((X1-XCEN),-DELTA)
      XB=min((X2-XCEN),DELTA)
      XFA=XF1-XCEN
C
      ZC0=YF1-XFA*ZS
      ZDI=1.E0_R8/DELTA
      ISIGN=-1
      IF(XA.LT.0.E0_R8) ISIGN=1
      ZC1=.5E0_R8*(ZS+ISIGN*ZC0*ZDI)
      ZC2=.3333333333333333333E0_R8*ISIGN*ZS*ZDI
C
      r8wxfint=ZDI*(XB-XA)*(ZC0+
     >           (XB+XA)*ZC1
     >          +(XB*XB+XA*XB+XA*XA)*ZC2)
C
      RETURN
      	END
C******************** END FILE WXFINT.FOR ; GROUP FILTR6 ******************
! 03Dec1999 fgtok -s r8_precision.sub r8smlib.sub "r8con.csh conversion"

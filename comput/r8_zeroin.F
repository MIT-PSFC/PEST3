#include "fpreproc/fp_const.h"
C******************** START FILE R8_ZEROIN.FOR ; GROUP TRKRLIB ************
CC
      REAL*8 FUNCTION R8_ZEROIN(AX,BX,F,TOL,EPSLON)
C
C============
C idecl:  explicitize implicit REAL declarations:
!
      IMPLICIT NONE
      REAL*8, PARAMETER :: ZERO = 0.0D0 ! _fgtok_constant_
      REAL*8, PARAMETER :: HALF = 0.5D0 ! _fgtok_constant_
      REAL*8, PARAMETER :: ONE = 1.0D0 ! _fgtok_constant_
      REAL*8, PARAMETER :: TWO = 2.0D0 ! _fgtok_constant_
      REAL*8, PARAMETER :: THREE = 3.0D0 ! _fgtok_constant_
!
      REAL*8 epslon,zansr
C============
      REAL*8 AX,BX,F,TOL
      EXTERNAL F
C
C      A ZERO OF THE FUNCTION  F(X)  IS COMPUTED IN THE INTERVAL AX,BX .
C
C  INPUT..
C
C  AX     LEFT ENDPOINT OF INITIAL INTERVAL
C  BX     RIGHT ENDPOINT OF INITIAL INTERVAL
C  F      FUNCTION SUBPROGRAM WHICH EVALUATES F(X) FOR ANY X IN
C         THE INTERVAL  AX,BX
C  TOL    DESIRED LENGTH OF THE INTERVAL OF UNCERTAINTY OF THE
C         FINAL RESULT ( .GE. 0.0)
C
C
C  OUTPUT..
C
C  R8_ZEROIN ABCISSA APPROXIMATING A ZERO OF  F  IN THE INTERVAL AX,BX
C
C
C      IT IS ASSUMED  THAT   F(AX)   AND   F(BX)   HAVE  OPPOSITE  SIGNS
C  WITHOUT  A  CHECK.  R8_ZEROIN  RETURNS A ZERO  X  IN THE GIVEN INTERVAL
C  AX,BX  TO WITHIN A TOLERANCE  4*MACHEPS*ABS(X) + TOL, WHERE MACHEPS
C  IS THE RELATIVE MACHINE PRECISION.
C      THIS FUNCTION SUBPROGRAM IS A SLIGHTLY  MODIFIED  TRANSLATION  OF
C  THE ALGOL 60 PROCEDURE  ZERO  GIVEN IN  RICHARD BRENT, ALGORITHMS FOR
C  MINIMIZATION WITHOUT DERIVATIVES, PRENTICE - HALL, INC. (1973).
C
C  dmc 10 Nov 1998 -- added a test to enforces AX.le.R8_ZEROIN.le.BX
C    or AX.ge.R8_ZEROIN.ge.BX on exit
C
C
      REAL*8  A,B,C,D,E,EPS,FA,FB,FC,TOL1,XM,P,Q,R,S
C
C  COMPUTE EPS, THE RELATIVE MACHINE PRECISION
C
C      EPS = 1.0
C   10 EPS = EPS/2.0
C      TOL1 = 1.0 + EPS
C      IF (TOL1 .GT. 1.0) GO TO 10
       EPS = EPSLON
C
C INITIALIZATION
C
      A = AX
      B = BX
      FA = F(A)
      FB = F(B)
C
C BEGIN STEP
C
   20 C = A
      FC = FA
      D = B - A
      E = D
   30 IF (ABS(FC) .GE. ABS(FB)) GO TO 40
      A = B
      B = C
      C = A
      FA = FB
      FB = FC
      FC = FA
C
C CONVERGENCE TEST
C
   40 TOL1 = TWO*EPS*ABS(B) + HALF*TOL
      XM = HALF*(C - B)
      IF (ABS(XM) .LE. TOL1) GO TO 90
      IF (FB .EQ. ZERO) GO TO 90
C
C IS BISECTION NECESSARY
C
      IF (ABS(E) .LT. TOL1) GO TO 70
      IF (ABS(FA) .LE. ABS(FB)) GO TO 70
C
C IS QUADRATIC INTERPOLATION POSSIBLE
C
      IF (A .NE. C) GO TO 50
C
C LINEAR INTERPOLATION
C
      S = FB/FA
      P = TWO*XM*S
      Q = ONE- S
      GO TO 60
C
C INVERSE QUADRATIC INTERPOLATION
C
   50 Q = FA/FC
      R = FB/FC
      S = FB/FA
      P = S*(TWO*XM*Q*(Q - R) - (B - A)*(R - ONE))
      Q = (Q - ONE)*(R - ONE)*(S - ONE)
C
C ADJUST SIGNS
C
   60 IF (P .GT. ZERO) Q = -Q
      P = ABS(P)
C
C IS INTERPOLATION ACCEPTABLE
C
      IF ((TWO*P) .GE. (THREE*XM*Q - ABS(TOL1*Q))) GO TO 70
      IF (P .GE. ABS(HALF*E*Q)) GO TO 70
      E = D
      D = P/Q
      GO TO 80
C
C BISECTION
C
   70 D = XM
      E = D
C
C COMPLETE STEP
C
   80 A = B
      FA = FB
      IF (ABS(D) .GT. TOL1) B = B + D
      IF (ABS(D) .LE. TOL1) B = B + SIGN(TOL1, XM)
      FB = F(B)
      IF ((FB*(FC/ABS(FC))) .GT. ZERO) GO TO 20
      GO TO 30
C
C DONE
C
   90 continue
      zansr = B                    ! dmc -- guarantee answer is in range
      zansr = max(zansr,min(ax,bx))
      zansr = min(zansr,max(ax,bx))
      R8_ZEROIN = zansr
C
      RETURN
      END
C******************** END FILE R8_ZEROIN.FOR ; GROUP TRKRLIB **************
! 11May2001 fgtok -s r8_precision.sub "r8con.csh conversion"
! 11May2001 fgtok

      SUBROUTINE scrclear
      USE mintrp
      USE scrunch_inc1
      USE scrunch_inc2
 
!...begin set dynamic arrays to zero
      IMPLICIT NONE
       rin = 0.0D0
       zin = 0.0D0
       xvec = 0.0D0
       gvec = 0.0D0
       xdot = 0.0D0
       xstore = 0.0D0
       m1 = 0
       n1 = 0
       rmomb = 0.0D0
       zmomb = 0.0D0
       mm = 0.0D0
       dm1 = 0.0D0
       faccon = 0.0D0
       xmpq = 0.0D0
       result = 0.0D0
!...end set dynamic arrays to zero
       RETURN
       END
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"

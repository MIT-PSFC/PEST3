!**********************************************************************
       SUBROUTINE pstboundy(is)
!**********************************************************************
!      purpose: to implement appropriate boundary conditions.
!
 USE pstcom

 USE l22com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IS


!................
!
!
!      wall on plasma surface  ??
!
      if( .not. wall) go to 10
!
!      first the m-1 ' th block
!
      CALL pstwallck
   10 continue
!
!      now work on the m'th block. add vacuum or fix for wall
!
      CALL pstaddvac(is)
!
      return
      end


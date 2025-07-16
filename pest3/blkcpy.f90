      SUBROUTINE pstblkcpy ( a,b,num )
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER NUM

!
!      routine to do blockcopies.... temporary.....
!
!***      dimension a(1),b(1)
!***
!
 COMPLEX*16, DIMENSION(*) 	 :: a
 COMPLEX*16, DIMENSION(*) 	 :: b

 b(1:num) = a(1:num)

!
      return
      end



      function pstvf (mn)
!     lcm ( vf )
 USE pstcom
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 pstVF
      INTEGER MN
      real*8 Y


!
!      linear piecewise elements...
!
      pstvf =  0.0_r8 
      if ( .not. lcub )  return
!
!      hermite cubic elements...
!
      y = psiv - psinod(mn)
      if ( y  >   0.0_r8  )  y = y / dtent(mn)
      if ( y  <   0.0_r8  )  y = y / dtent(mn-1)
      pstvf = y * ( abs(y)- 1.0_r8  ) ** 2
      return
      end



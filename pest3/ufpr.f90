      function pstufpr ( mn )
!     lcm ( ufpr )
 USE pstcom
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 pstUFPR
      INTEGER MN
      real*8 PSIS
      real*8 Y
      real(r8),parameter ::  epsilon=1.12487e-6
!
      psis = psiv
      if(psiv  <  psivax) psiv = psivax
      y = psiv - psinod(mn)
      if ( lcub )  go to 20
!
!      linear piecewise elements...
!
      if ( y  >   0.0_r8  )  pstufpr = - 1.0_r8  / dtent(mn)
      if ( y  <   0.0_r8  )  pstufpr =   1.0_r8  / dtent(mn-1)
      if ( abs(y)  <  half*bit )  pstufpr =  0.0_r8 
!
      if ( mn  ==  1 )  then 
         pstufpr =  0.5_r8 *(  1.0_r8 - 3.0_r8 *psiv/dtent(1) ) &
  &                          / (sqrt(psiv )+epsilon)
      endif
      psiv = psis
      return
   20 continue
!
!      hermite cubic elements...
!
      if ( y  >=   0.0_r8  )  go to 30
      y = abs(y) / dtent(mn-1)
      pstufpr = - 6.0_r8 *y*( y- 1.0_r8  ) / dtent(mn-1)
      go to 40
   30 y = abs(y) / dtent(mn)
      pstufpr =  6.0_r8 *y*( y- 1.0_r8  ) / dtent(mn)
   40 continue
!
      psiv = psis
      return
      end



      function pstuf (mn)
!------------------------------------------------------------
!     lcm(uf)
 USE pstcom
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 pstUF
      INTEGER MN
      real*8 PSIS
      real*8 Y


!
      psis = psiv
      if(psiv  <  psivax) psiv = psivax
      y = psiv - psinod(mn)
      if ( y  >   0.0_r8  )  y = y / dtent(mn)
      if ( y  <   0.0_r8  )  y = y / dtent(mn-1)
      if ( lcub )  go to 20
!
!      piecewise linear elements...
!
      pstuf =  1.0_r8  - abs(y)
!
      if ( mn  ==  1 )  pstuf = sqrt(psiv)*(  1.0_r8 -psiv/dtent(1) )
      psiv = psis
!
      return
   20 continue
!
!      hermite cubic elements...
!
      y = abs(y)
      pstuf = ( y- 1.0_r8  )**2 * ( 2.0_r8 *y+ 1.0_r8 )
!
      psiv = psis
      return
      end



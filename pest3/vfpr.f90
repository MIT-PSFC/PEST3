      function pstvfpr (mn)
!     lcm ( vfpr)
 USE pstcom
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 pstVFPR
      INTEGER MN
      real*8 Y
      real*8 DTE


!
!      linear piecewise elements...
!
      pstvfpr =  0.0_r8 
      if ( .not. lcub )  return
!
!      hermite cubic elements...
!
      y = psiv - psinod(mn)
      if ( y  >=   0.0_r8  )  go to 20
      y = abs(y) / dtent(mn-1)
      dte = dtent(mn-1)
      go to 30
   20 y = abs(y) / dtent(mn)
      dte = dtent(mn)
   30 pstvfpr = (  3.0_r8 *y- 1.0_r8  ) * ( y- 1.0_r8  ) / dte
!
      return
      end



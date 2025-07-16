      SUBROUTINE pstloadw( kth1, zf, zw )
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER KTH1
      INTEGER KTH2

!
!
 REAL*8, DIMENSION(*) 	 :: zf
 REAL*8, DIMENSION(*) 	 :: zw
!
      kth2 = kth1 + 1
      zw(1:kth2) =  0.0_r8 
!
      zw(1:kth1) = zf(1:kth1)
!
      return
	  end



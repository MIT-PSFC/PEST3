!--------------------------------------------------------------
      subroutine pstsolve( n, ndim, a, b, ip )
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER N
      INTEGER NDIM
      INTEGER NM1
      INTEGER K
      INTEGER KP1
      INTEGER M
      real*8 T
      INTEGER I
      INTEGER KB
      INTEGER KM1

!--------------------------------------------------------------
!     lcm(solve)
 REAL*8, DIMENSION(ndim,ndim) 	 :: a
 REAL*8, DIMENSION(ndim) 	 :: b
 INTEGER, DIMENSION(ndim) 	 :: ip
      if( n  ==  1 ) go to 9
      nm1 = n - 1
      do k = 1, nm1
         kp1 = k + 1
         m = ip(k)
         t = b(m)
         b(m) = b(k)
         b(k) = t
         do   i = kp1, n
            b(i) = b(i) + a(i,k) * t
         end do
      end do
      do kb = 1, nm1
         km1 = n - kb
         k = km1 + 1
         b(k) = b(k) / a(k,k)
         t = -b(k)
         do   i = 1, km1
            b(i) = b(i) + a(i,k) * t
         end do
      end do
    9 b(1) = b(1) / a(1,1)
      return
      end



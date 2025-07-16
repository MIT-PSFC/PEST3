!-------------------------------------------------------------------
      SUBROUTINE pstdecomp( n, ndim, a, ip )
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER N
      INTEGER NDIM
      INTEGER K
      INTEGER KP1
      INTEGER M
      INTEGER I
      real*8 T
      INTEGER J

!-------------------------------------------------------------------
!     lcm(decomp)
!     communications of the acm, volume 15, #4, april, 1972, page 274
!     algorithm 423 by cleve b. moler
!
 REAL*8, DIMENSION(ndim,ndim) 	 :: a
 INTEGER, DIMENSION(ndim) 	 :: ip
      ip(n) = 1
      do 6   k = 1, n
      if( k  ==  n ) go to 5
      kp1 = k + 1
      m = k
      do 1   i = kp1, n
      if( abs(a(i,k))  >   abs(a(m,k))) m = i
    1 continue
      ip(k) = m
      if(m  /=  k) ip(n) = - ip(n)
      t = a(m,k)
      a(m,k) = a(k,k)
      a(k,k) = t
      if( t  ==   0.0_r8  ) go to 5
      do i = kp1, n
         a(i,k) = -a(i,k)/t
      end do
      do 4   j = kp1, n
      t = a(m,j)
      a(m,j) = a(k,j)
      a(k,j) = t
      if( t  ==   0.0_r8  ) go to 4
      do    i = kp1, n
         a(i,j) = a(i,j) + a(i,k) * t
      end do
    4 continue
    5 if( a(k,k)  ==   0.0_r8  ) ip(n) = 0
    6 continue
      return
      end



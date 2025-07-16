      SUBROUTINE pstdecomc( n, ndim, a, ip )
!-------------------------------------------------------------------
!     communications of the acm, volume 15, #4, april, 1972, page 274
!     algorithm 423 by cleve b. moler
!
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER N
      INTEGER NDIM
      INTEGER K
      INTEGER KP1
      INTEGER M
      INTEGER I
      INTEGER J

 COMPLEX*16, DIMENSION(ndim,ndim) 	 :: a
 COMPLEX*16	 	 	 :: t
 INTEGER, DIMENSION(ndim) 	 :: ip
!
 ip(n) = 1
 do   k = 1, n
    if( k  ==  n ) go to 5
    kp1 = k + 1
    m = k
    do   i = kp1, n
       if( abs(a(i,k))  >   abs(a(m,k))) m = i
    end do
    ip(k) = m
    if(m  /=  k) ip(n) = - ip(n)
    t = a(m,k)
    a(m,k) = a(k,k)
    a(k,k) = t
    if( t  ==  ( 0._r8 , 0._r8 )  ) go to 5
    do   i = kp1, n
       a(i,k) = -a(i,k)/t
    end do
    do j = kp1, n
       t = a(m,j)
       a(m,j) = a(k,j)
       a(k,j) = t
       do i = kp1, n
          a(i,j) = a(i,j) + a(i,k) * t
       end do
    end do
5   if( a(k,k)  ==  ( 0._r8 , 0._r8 )  ) ip(n) = 0
 end do
      return
      end



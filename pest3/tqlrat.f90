      subroutine pstTqlrat(n,d,e2,ierr)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)

!
!
!     this SUBROUTINE pstis a translation of the algol procedure tqlrat,
!     algorithm 464, comm. acm 16, 689(1973) by reinsch.
!
!     this SUBROUTINE pstfinds the eigenvalues of a symmetric
!     tridiagonal matrix by the rational ql method.
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e2 contains the squares of the subdiagonal elements of the
!          input matrix in its last n-1 positions.  e2(1) is arbitrary.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        e2 has been destroyed.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august  1983._r8 
!
!     ------------------------------------------------------------------
!
 INTEGER	 	 	 :: n
 REAL*8, DIMENSION(n) 	 :: d
 REAL*8, DIMENSION(n) 	 :: e2
 REAL*8	 	 	 :: b
 REAL*8	 	 	 :: c
 REAL*8	 	 	 :: f
 REAL*8	 	 	 :: g
 REAL*8	 	 	 :: h
 REAL*8	 	 	 :: p
 REAL*8	 	 	 :: r
 REAL*8	 	 	 :: s
 REAL*8	 	 	 :: t
 REAL*8	 	 	 :: epslon
 REAL*8	 	 	 :: pythag
 INTEGER	 	 	 :: i
 INTEGER	 	 	 :: j
 INTEGER	 	 	 :: l
 INTEGER	 	 	 :: m
 INTEGER	 	 	 :: ii
 INTEGER	 	 	 :: l1
 INTEGER	 	 	 :: mml
 INTEGER	 	 	 :: ierr
      ierr = 0
      if (n  ==  1) go to 1001
!
      do i = 2, n
        e2(i-1) = e2(i)
      end do
!
      f =  0.0E0_r8  
      t =  0.0E0_r8  
      e2(n) =  0.0E0_r8  
!
      do 290 l = 1, n
         j = 0
         h = abs(d(l)) + sqrt(e2(l))
         if (t  >  h) go to 105
         t = h
         b = epslon(t)
         c = b * b
!     .......... look for small squared sub-diagonal element ..........
  105    do 110 m = l, n
            if (e2(m)  <=  c) go to 120
!     .......... e2(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
!
  120    if (m  ==  l) go to 210
  130    if (j  ==  30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         s = sqrt(e2(l))
         g = d(l)
         p = (d(l1) - g) / ( 2.0E0_r8   * s)
         r = pythag(p, 1.0e0_r8 ) 
         d(l) = s / (p + sign(r,p))
         h = g - d(l)
!
         do i = l1, n
          d(i) = d(i) - h
         end do
!
         f = f + h
!     .......... rational ql transformation ..........
         g = d(m)
         if (g  ==   0.0E0_r8  ) g = b
         h = g
         s =  0.0E0_r8  
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            p = g * h
            r = p + e2(i)
            e2(i+1) = s * r
            s = e2(i) / r
            d(i+1) = h + s * (h + d(i))
            g = d(i) - e2(i) / g
            if (g  ==   0.0E0_r8  ) g = b
            h = g * p / r
  200    continue
!
         e2(l) = s * g
         d(l) = h
!     .......... guard against underflow in convergence test ..........
         if (h  ==   0.0E0_r8  ) go to 210
         if (abs(e2(l))  <=  abs(c/h)) go to 210
         e2(l) = h * e2(l)
         if (e2(l)  /=   0.0E0_r8  ) go to 130
  210    p = d(l) + f
!     .......... order eigenvalues ..........
         if (l  ==  1) go to 250
!     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p  >=  d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
!
  250    i = 1
  270    d(i) = p
  290 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end



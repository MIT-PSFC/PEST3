      subroutine pstTql2(nm,n,d,e,z,ierr)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
!
!
!     this SUBROUTINE pstis a translation of the algol procedure tql2,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this SUBROUTINE pstfinds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2, ..._r8 , ierr- 1._r8 
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
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
 INTEGER	 	 	 :: nm
 REAL*8, DIMENSION(n) 	 :: d
 REAL*8, DIMENSION(n) 	 :: e
 REAL*8, DIMENSION(nm,n) 	 :: z
 REAL*8	 	 	 :: c
 REAL*8	 	 	 :: c2
 REAL*8	 	 	 :: c3
 REAL*8	 	 	 :: dl1
 REAL*8	 	 	 :: el1
 REAL*8	 	 	 :: f
 REAL*8	 	 	 :: g
 REAL*8	 	 	 :: h
 REAL*8	 	 	 :: p
 REAL*8	 	 	 :: r
 REAL*8	 	 	 :: s
 REAL*8	 	 	 :: s2
 REAL*8	 	 	 :: tst1
 REAL*8	 	 	 :: tst2
 REAL*8	 	 	 :: pythag
 INTEGER	 	 	 :: i
 INTEGER	 	 	 :: j
 INTEGER	 	 	 :: k
 INTEGER	 	 	 :: l
 INTEGER	 	 	 :: m
 INTEGER	 	 	 :: ii
 INTEGER	 	 	 :: l1
 INTEGER	 	 	 :: l2
 INTEGER	 	 	 :: mml
 INTEGER	 	 	 :: ierr
      ierr = 0
      if (n  ==  1) go to 1001
!
      do i = 2, n
        e(i-1) = e(i)
      end do
!
      f =  0.0E0_r8  
      tst1 =  0.0E0_r8  
      e(n) =  0.0E0_r8  
!
      do 240 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (tst1  <  h) tst1 = h
!     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + abs(e(m))
            if (tst2  ==  tst1) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
!
  120    if (m  ==  l) go to 220
  130    if (j  ==  30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / ( 2.0E0_r8   * e(l))
         r = pythag(p, 1.0e0_r8 ) 
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2  >  n) go to 145
!
         do i = l2, n
          d(i) = d(i) - h
         end do
!
  145    f = f + h
!     .......... ql transformation ..........
         p = d(m)
         c =  1.0E0_r8  
         c2 = c
         el1 = e(l1)
         s =  0.0E0_r8  
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
!     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
!
  200    continue
!
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2  >  tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
!     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
!
         do 260 j = ii, n
            if (d(j)  >=  p) go to 260
            k = j
            p = d(j)
  260    continue
!
         if (k  ==  i) go to 300
         d(k) = d(i)
         d(i) = p
!
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
!
  300 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end



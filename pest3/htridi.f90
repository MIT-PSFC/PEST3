      SUBROUTINE psthtridi(nm,n,ar,ai,d,e,e2,tau)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
!
!
!     this SUBROUTINE pstis a translation of a complex anlogue of
!     the algol procedure tred1, num. math. 11, 181-195(1968)
!     by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this SUBROUTINE pstreduces a complex hermitian matrix
!     to a real symmetric tridiagonal matrix using
!     unitary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex hermitian input matrix.
!          only the lower triangle of the matrix need be supplied.
!
!     on output
!
!        ar and ai contain information about the unitary trans-
!          formations used in the reduction in their full lower
!          triangles.  their strict upper triangles and the
!          diagonal of ar are unaltered.
!
!        d contains the diagonal elements of the the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!        tau contains further information about the transformations.
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
 REAL*8, DIMENSION(nm,n) 	 :: ar
 REAL*8, DIMENSION(nm,n) 	 :: ai
 REAL*8, DIMENSION(n) 	 :: d
 REAL*8, DIMENSION(n) 	 :: e
 REAL*8, DIMENSION(n) 	 :: e2
 REAL*8, DIMENSION(2,n) 	 :: tau
 REAL*8	 	 	 :: f
 REAL*8	 	 	 :: g
 REAL*8	 	 	 :: h
 REAL*8	 	 	 :: fi
 REAL*8	 	 	 :: gi
 REAL*8	 	 	 :: hh
 REAL*8	 	 	 :: si
 REAL*8	 	 	 :: scale
 REAL*8	 	 	 :: pythag
 INTEGER	 	 	 :: i
 INTEGER	 	 	 :: j
 INTEGER	 	 	 :: k
 INTEGER	 	 	 :: l
 INTEGER	 	 	 :: ii
 INTEGER	 	 	 :: jp1
      tau(1,n) =  1.0E0_r8  
      tau(2,n) =  0.0E0_r8  
!
      do i = 1, n
         d(i) = ar(i,i)
      end do
!     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h =  0.0E0_r8  
         scale =  0.0E0_r8  
         if (l  <  1) go to 130
!     .......... scale row (algol tol then not needed) ..........
         do k = 1, l
           scale = scale + abs(ar(i,k)) + abs(ai(i,k))
         end do
!
         if (scale  /=   0.0E0_r8  ) go to 140
         tau(1,l) =  1.0E0_r8  
         tau(2,l) =  0.0E0_r8  
  130    e(i) =  0.0E0_r8  
         e2(i) =  0.0E0_r8  
         go to 290
!
  140    do 150 k = 1, l
            ar(i,k) = ar(i,k) / scale
            ai(i,k) = ai(i,k) / scale
            h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
  150    continue
!
         e2(i) = scale * scale * h
         g = sqrt(h)
         e(i) = scale * g
         f = pythag(ar(i,l),ai(i,l))
!     .......... form next diagonal element of matrix t ..........
         if (f  ==   0.0E0_r8  ) go to 160
         tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
         si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
         h = h + f * g
         g =  1.0E0_r8   + g / f
         ar(i,l) = g * ar(i,l)
         ai(i,l) = g * ai(i,l)
         if (l  ==  1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         ar(i,l) = g
  170    f =  0.0E0_r8  
!
         do 240 j = 1, l
            g =  0.0E0_r8  
            gi =  0.0E0_r8  
!     .......... form element of a*u ..........
            do 180 k = 1, j
               g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
               gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
  180       continue
!
            jp1 = j + 1
            if (l  <  jp1) go to 220
!
            do 200 k = jp1, l
               g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
               gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
  200       continue
!     .......... form element of p ..........
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
  240    continue
!
         hh = f / (h + h)
!     .......... form reduced a ..........
         do j = 1, l
            f = ar(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = -ai(i,j)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
!
            do k = 1, j
               ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k) &
                           + fi * tau(2,k) + gi * ai(i,k)
               ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k) &
                           - fi * e(k) - gi * ar(i,k)
            end do
         end do
!
  270    do 280 k = 1, l
            ar(i,k) = scale * ar(i,k)
            ai(i,k) = scale * ai(i,k)
  280    continue
!
         tau(2,l) = -si
  290    hh = d(i)
         d(i) = ar(i,i)
         ar(i,i) = hh
         ai(i,i) = scale * sqrt(h)
  300 continue
!
      return
      end



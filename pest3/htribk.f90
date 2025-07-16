      SUBROUTINE psthtribk(nm,n,ar,ai,tau,m,zr,zi)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)

!
!
!     this SUBROUTINE pstis a translation of a complex anlogue of
!     the algol procedure trbak1, num. math. 11, 181-195(1968)
!     by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this SUBROUTINE pstforms the eigenvectors of a complex hermitian
!     matrix by back transforming those of the corresponding
!     real symmetric tridiagonal matrix determined by  htridi.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        ar and ai contain information about the unitary trans-
!          formations used in the reduction by  htridi  in their
!          full lower triangles except for the diagonal of ar.
!
!        tau contains further information about the transformations.
!
!        m is the number of eigenvectors to be back transformed.
!
!        zr contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
!
!     note that the last component of each returned vector
!     is real and that vector euclidean norms are preserved.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august  1983._r8 
!
!     ------------------------------------------------------------------
!
 INTEGER	 	 	 :: m
 INTEGER	 	 	 :: n
 INTEGER	 	 	 :: nm
 REAL*8, DIMENSION(nm,n) 	 :: ar
 REAL*8, DIMENSION(nm,n) 	 :: ai
 REAL*8, DIMENSION(2,n) 	 :: tau
 REAL*8, DIMENSION(nm,m) 	 :: zr
 REAL*8, DIMENSION(nm,m) 	 :: zi
 REAL*8	 	 	 :: h
 REAL*8	 	 	 :: s
 REAL*8	 	 	 :: si
 INTEGER	 	 	 :: i
 INTEGER	 	 	 :: j
 INTEGER	 	 	 :: k
 INTEGER	 	 	 :: l
      if (m  ==  0) go to 200
!     .......... transform the eigenvectors of the real symmetric
!                tridiagonal matrix to those of the hermitian
!                tridiagonal matrix. ..........
      do k = 1, n
!
         do j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
         end do
      end do
!
      if (n  ==  1) go to 200
!     .......... recover and apply the householder matrices ..........
      do 140 i = 2, n
         l = i - 1
         h = ai(i,i)
         if (h  ==   0.0E0_r8  ) go to 140
!
         do 130 j = 1, m
            s =  0.0E0_r8  
            si =  0.0E0_r8  
!
            do 110 k = 1, l
               s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
               si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
  110       continue
!     .......... double divisions avoid possible underflow ..........
            s = (s / h) / h
            si = (si / h) / h
!
            do 120 k = 1, l
               zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
               zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
  120       continue
!
  130    continue
!
  140 continue
!
  200 return
      end




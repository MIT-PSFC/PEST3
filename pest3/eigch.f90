      SUBROUTINE psteigch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)

      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
!
!
!     this SUBROUTINE pstcalls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a complex hermitian matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a=(ar,ai).
!
!        ar  and  ai  contain the real and imaginary parts,
!        respectively, of the complex hermitian matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        zr  and  zi  contain the real and imaginary parts,
!        respectively, of the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1, fv2, and  fm1  are temporary storage arrays.
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
 REAL*8, DIMENSION(n) 	 :: w
 REAL*8, DIMENSION(nm,n) 	 :: zr
 REAL*8, DIMENSION(nm,n) 	 :: zi
 REAL*8, DIMENSION(n) 	 :: fv1
 REAL*8, DIMENSION(n) 	 :: fv2
 REAL*8, DIMENSION(2,n) 	 :: fm1
 INTEGER	 	 	 :: i
 INTEGER	 	 	 :: j
 INTEGER	 	 	 :: ierr
 INTEGER	 	 	 :: matz
      if (n  <=  nm) go to 10
      ierr = 10 * n
      go to 50
!
   10 CALL psthtridi(nm,n,ar,ai,w,fv1,fv2,fm1)
      if (matz  /=  0) go to 20
!     .......... find eigenvalues only ..........
      CALL psttqlrat(n,w,fv2,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
!
         do 30 j = 1, n
            zr(j,i) =  0.0E0_r8  
   30    continue
!
         zr(i,i) =  1.0E0_r8  
   40 continue
!
      CALL psttql2(nm,n,w,fv1,zr,ierr)
      if (ierr  /=  0) go to 50
      CALL psthtribk(nm,n,ar,ai,fm1,n,zr,zi)
   50 return
      end



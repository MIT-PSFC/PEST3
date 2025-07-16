 REAL*8 function epslon (x)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)

!
!     estimate unit roundoff in quantities of size x.
!
!
!     this program should function properly on all systems
!     satisfying the following two assumptions,
!         1._r8   the base used in representing floating point
!            numbers is not a power of three.
!         2._r8   the quantity  a  in statement 10 is represented to
!            the accuracy used in floating point variables
!            that are stored in memory.
!     the statement number 10 and the go to 10 are intended to
!     force optimizing compilers to generate code satisfying
!     assumption  2._r8 
!     under these assumptions, it should be true that,
!            a  is not exactly equal to four-thirds,
!            b  has a zero for its last bit or digit,
!            c  is not exactly equal to one,
!            eps  measures the separation of  1.0_r8  from
!                 the next larger floating point number.
!     the developers of eispack would appreciate being informed
!     about any systems where these assumptions do not hold.
!
!     this version dated 4/6/ 83._r8 
!
 REAL*8	 	 	 :: x
 REAL*8	 	 	 :: a
 REAL*8	 	 	 :: b
 REAL*8	 	 	 :: c
 REAL*8	 	 	 :: eps
      a =  4.0E0_r8  / 3.0E0_r8  
   10 b = a -  1.0E0_r8  
      c = b + b + b
      eps = abs(c- 1.0E0_r8  )
      if (eps  ==   0.0E0_r8  ) go to 10
      epslon = eps*abs(x)
      return
      end



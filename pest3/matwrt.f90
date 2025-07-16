!
!......................................................................
      SUBROUTINE pstmatwrt ( a, maxj1, maxj2, jmx1, jmx2, label )
!
!        writes out the elements of matrix a.  label must be
!        20 letters long.
!
!     integer label ( 2 )
!
 USE pstcom
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER MAXJ1
      INTEGER MAXJ2
      INTEGER JMX1
      INTEGER JMX2
      INTEGER J1
      INTEGER J2


!
      character*(*) label
!
 REAL*8, DIMENSION(maxj1,maxj2) 	 :: a
 INTEGER	 	 	 :: out
      out = outmod
!
      write ( out, 10 ) label
   10 format ( //, 5x, "matrix elements of  ", a40 )
!
      do 30 j1 = 1, jmx1
      write ( out, 20 ) ( a(j1,j2),j2 = 1, jmx2 )
   20 format ( 10(8e11.4,/) )
   30 continue
!
      return
      end


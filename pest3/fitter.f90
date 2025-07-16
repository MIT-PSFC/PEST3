      SUBROUTINE pstfitter(npts,e0,psif,arrf,coef)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER NPTS
      real*8 E0
      INTEGER I
      real*8 PEXP
      INTEGER J
      real*8 DET
      INTEGER IDET
      real*8 S

!......................................................................
!      fits npts points to a polynomial in sqrt(psi). the leading term
!      of the polynomial is (psi)**(e0/2). the coefficients coef are
!      determined by matching the values of the array arrf at the points
!      psif.      jm. march 1981
!......................................................................
!      load the work matrix a
 REAL*8, DIMENSION(10,10) 	 :: a
 REAL*8, DIMENSION(10) 	 :: b
 REAL*8, DIMENSION(*) 	 :: psif
 REAL*8, DIMENSION(*) 	 :: arrf
 REAL*8, DIMENSION(*) 	 :: coef
 REAL*8, DIMENSION(10,10) 	 :: as
 INTEGER, DIMENSION(10) 	 :: ip
      do 20 i=1, npts
!      set the exponent
      pexp = (e0 + i - 1._r8 ) /  2._r8 
      do 10 j=1, npts
      a(i,j) =  1.0_r8 
      if(pexp  ==   0._r8 ) go to 10
      a(i,j) = (psif(j))**pexp
   10 continue
   20 continue
!
!      now decompose a
!
      CALL pstdecomp(npts,10,a,ip)
!      check for singularity
!
      det =  0.0_r8 
      do i=1, npts
         det = det + log10(abs(a(i,i)))
      end do
      idet = INT(det)
      det  =  10._r8 **(det-idet)
      s    =  1.0_r8 
      do i=1, npts
         s    = s * a(i,i)/abs(a(i,i))
      end do
      s    = s * ip(npts)
      if(s  >    0._r8 ) go to 50
      write(6,8000)s,det,idet
!     CALL psterrmes(6,6hfitter)
      return
   50 continue
!
!      compute coefficients
!
      do 70 i=1, npts
      do 60 j=1, npts
      b(j) =  0.0_r8 
      if(i  ==  j)b(j) =  1.0_r8 
   60 continue
!
      call pstsolve(npts,10,a,b,ip)
!
      do j=1, npts
         as(j,i) = b(j)
      end do
!
   70 continue
!
      do j=1, npts
         coef(j) =  0.0_r8 
         do i=1, npts
            coef(j) = coef(j) + as(i,j)*arrf(i)
         end do
      end do
!
      return
!
!
 8000 format(' *** fitter: determinant = ',f4.1,' l.e. 0 ',   &
     f5.2,i5/)
      end



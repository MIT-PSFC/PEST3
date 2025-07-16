!     ..................................................................
!
!        description of parameters
!           a - original matrix (symmetric), destroyed in computation.
!               resultant eigenvalues are developed in diagonal of
!               matrix a in descending order.
!           r - resultant matrix of eigenvectors (stored columnwise,
!               in same sequence as eigenvalues)
!           n - order of matrices a and r
!           mv- input code
!                   0   compute eigenvalues and eigenvectors
!                   1   compute eigenvalues only (r need not be
!                       dimensioned but must still appear in calling
!                       sequence)
!
!        remarks
!           original matrix a must be real symmetric (storage mode=1)
!           matrix a cannot be in the same location as matrix r
!
!        subroutines and function subprograms required
!           none
!
!        method
!           diagonalization method originated by jacobi and adapted
!           by von neumann for large computers as found in "mathematical
!           methods for digital computers", edited by a. ralston and
!           h.s. wilf, john wiley and sons, new york, 1962, chapter 7
!
!     ..................................................................
!
      SUBROUTINE psteigen(a,r,n,mv)
!     lcm  (eigen)
!
!        ...............................................................
!
!        if a double precision version of this routine is desired, the
!        c in column 1 should be removed from the double precision
!        statement which follows.
!
!
!        the c must also be removed from double precision statements
!        appearing in other routines used in conjunction with this
!        routine.
!
!        the double precision version of this SUBROUTINE pstmust also
!        contain double precision fortran functions.  sqrt in statements
!        40, 68, 75, and 78 must be changed to dsqrt.  abs in statement
!        62 must be changed to abs. the constant in statement 5 should
!        be changed to 1.0d- 12._r8 
!
!        ...............................................................
!
!        generate identity matrix
!
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER N
      INTEGER MV
      real*8 RANGE
      INTEGER IQ
      INTEGER J
      INTEGER I
      INTEGER IJ
      real*8 ANORM
      INTEGER IA
      real*8 ANRMX
      INTEGER IND
      real*8 THR
      INTEGER L
      INTEGER M
      INTEGER MQ
      INTEGER LQ
      INTEGER LM
      INTEGER LL
      INTEGER MM
      real*8 X
      real*8 Y
      real*8 YP
      real*8 SINX
      real*8 SINX2
      real*8 COSX
      real*8 COSX2
      real*8 SINCS
      INTEGER ILQ
      INTEGER IMQ
      INTEGER IM
      INTEGER IL
      INTEGER ILR
      INTEGER IMR
      INTEGER JQ
      INTEGER K

 REAL*8, DIMENSION(*) 	 :: a
 REAL*8, DIMENSION(*) 	 :: r
    5 range= 1.0E-12_r8  
      if(mv-1) 10,25,10
   10 iq=-n
      do 20 j=1,n
      iq=iq+n
      do 20 i=1,n
      ij=iq+i
      r(ij)= 0.0_r8 
      if(i-j) 20,15,20
   15 r(ij)= 1.0_r8 
   20 continue
!
!        compute initial and final norms (anorm and anormx)
!
   25 anorm= 0.0_r8 
      do 35 i=1,n
      do 35 j=i,n
      if(i-j) 30,35,30
   30 ia=i+(j*j-j)/2
      anorm=anorm+a(ia)*a(ia)
   35 continue
      if(anorm) 165,165,40
   40 anorm= 1.414_r8 *sqrt(anorm)
      anrmx=anorm*range/float(n)
!
!        initialize indicators and compute threshold, thr
!
      ind=0
      thr=anorm
   45 thr=thr/float(n)
   50 l=1
   55 m=l+1
!
!        compute sin and cos
!
   60 mq=(m*m-m)/2
      lq=(l*l-l)/2
      lm=l+mq
   62 if( abs(a(lm))-thr) 130,65,65
   65 ind=1
      ll=l+lq
      mm=m+mq
      x= 0.5_r8 *(a(ll)-a(mm))
   68 y=-a(lm)/ sqrt(a(lm)*a(lm)+x*x)
      if(x) 70,75,75
   70 y=-y
   75 continue
      yp =  1.0_r8  - y*y
      if ( abs(yp)  <   1.0E-10_r8   ) yp =  0.0_r8 
      sinx=y/ sqrt( 2.0_r8 *( 1.0_r8 +sqrt(yp+ 0.0_r8 )))
      sinx2=sinx*sinx
   78 cosx= sqrt( 1.0_r8 -sinx2)
      cosx2=cosx*cosx
      sincs =sinx*cosx
!
!        rotate l and m columns
!
      ilq=n*(l-1)
      imq=n*(m-1)
      do 125 i=1,n
      iq=(i*i-i)/2
      if(i-l) 80,115,80
   80 if(i-m) 85,115,90
   85 im=i+mq
      go to 95
   90 im=m+iq
   95 if(i-l) 100,105,105
  100 il=i+lq
      go to 110
  105 il=l+iq
  110 x=a(il)*cosx-a(im)*sinx
      a(im)=a(il)*sinx+a(im)*cosx
      a(il)=x
  115 if(mv-1) 120,125,120
  120 ilr=ilq+i
      imr=imq+i
      x=r(ilr)*cosx-r(imr)*sinx
      r(imr)=r(ilr)*sinx+r(imr)*cosx
      r(ilr)=x
  125 continue
      x= 2.0_r8 *a(lm)*sincs
      y=a(ll)*cosx2+a(mm)*sinx2-x
      x=a(ll)*sinx2+a(mm)*cosx2+x
      a(lm)=(a(ll)-a(mm))*sincs+a(lm)*(cosx2-sinx2)
      a(ll)=y
      a(mm)=x
!
!        tests for completion
!
!        test for m = last column
!
  130 if(m-n) 135,140,135
  135 m=m+1
      go to 60
!
!        test for l = second from last column
!
  140 if(l-(n-1)) 145,150,145
  145 l=l+1
      go to 55
  150 if(ind-1) 160,155,160
  155 ind=0
      go to 50
!
!        compare threshold with final norm
!
  160 if(thr-anrmx) 165,165,45
!
!        sort eigenvalues and eigenvectors
!
  165 iq=-n
      do 185 i=1,n
      iq=iq+n
      ll=i+(i*i-i)/2
      jq=n*(i-2)
      do 185 j=i,n
      jq=jq+n
      mm=j+(j*j-j)/2
      if(a(ll)-a(mm)) 170,185,185
  170 x=a(ll)
      a(ll)=a(mm)
      a(mm)=x
      if(mv-1) 175,185,175
  175 do 180 k=1,n
      ilr=iq+k
      imr=jq+k
      x=r(ilr)
      r(ilr)=r(imr)
  180 r(imr)=x
  185 continue
      return
      end


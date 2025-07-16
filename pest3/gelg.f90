!        SUBROUTINE pstgelg
!
!        purpose
!           to solve a general system of simultaneous linear equations.
!
!        usage
!           CALL pstgelg(r,a,m,n,eps,ier)
!
!        description of parameters
!           r      - the m by n matrix of right hand sides.  (destroyed)
!                    on return r contains the solution of the equations.
!           a      - the m by m coefficient matrix.  (destroyed)
!           m      - the number of equations in the system.
!           n      - the number of right hand side vectors.
!           eps    - an input constant which is used as relative
!                    tolerance for test on loss of significance.
!           ier    - resulting error parameter coded as follows
!                    ier=0  - no error,
!                    ier=-1 - no result because of m less than 1 or
!                             pivot element at any elimination step
!                             equal to 0,
!                    ier=k  - warning due to possible loss of signifi-
!                             cance indicated at elimination step k+1,
!                             where pivot element was less than or
!                             equal to the internal tolerance eps times
!                             absolutely greatest element of matrix a.
!
!        remarks
!           input matrices r and a are assumed to be stored columnwise
!           in m*n resp. m*m successive storage locations. on return
!           solution matrix r is stored columnwise too.
!           the procedure gives results if the number of equations m is
!           greater than 0 and pivot elements at all elimination steps
!           are different from  0._r8  however warning ier=k - if given -
!           indicates possible loss of significance. in case of a well
!           scaled matrix a and appropriate tolerance eps, ier=k may be
!           interpreted that matrix a has the rank k. no warning is
!           given in case m= 1._r8 
!
!        subroutines and function subprograms required
!           none
!
!        method
!           solution is done by means of gauss-elimination with
!           complete pivoting.
!
!     ..................................................................
!
      SUBROUTINE pstgelg(r,a,m,n,eps,ier)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER M
      INTEGER N
      real*8 EPS
      INTEGER IER
      real*8 PIV
      INTEGER MM
      INTEGER NM
      INTEGER L
      real*8 TB
      INTEGER I
      real*8 TOL
      INTEGER LST
      INTEGER K
      real*8 PIVI
      INTEGER J
      INTEGER LL
      INTEGER LEND
      INTEGER II
      INTEGER IST

!     lcm (gelg)
!
!
 REAL*8, DIMENSION(*) 	 :: a
 REAL*8, DIMENSION(*) 	 :: r
!
      if(m)23,23,1
!
!     search for greatest element in matrix a
    1 ier=0
      piv= 0._r8 
      mm=m*m
      nm=n*m
      do 3 l=1,mm
      tb=abs(a(l))
      if(tb-piv)3,3,2
    2 piv=tb
      i=l
    3 continue
      tol=eps*piv
!     a(i) is pivot element. piv contains the absolute value of a(i).
!
!
!     start elimination loop
      lst=1
      do 17 k=1,m
!
!     test on singularity
      if(piv)23,23,4
    4 if(ier)7,5,7
    5 if(piv-tol)6,6,7
    6 ier=k-1
    7 pivi= 1._r8 /a(i)
      j=(i-1)/m
      i=i-j*m-k
      j=j+1-k
!     i+k is row-index, j+k column-index of pivot element
!
!     pivot row reduction and row interchange in right hand side r
      do 8 l=k,nm,m
      ll=l+i
      tb=pivi*r(ll)
      r(ll)=r(l)
    8 r(l)=tb
!
!     is elimination terminated
      if(k-m)9,18,18
!
!     column interchange in matrix a
    9 lend=lst+m-k
      if(j)12,12,10
   10 ii=j*m
      do 11 l=lst,lend
      tb=a(l)
      ll=l+ii
      a(l)=a(ll)
   11 a(ll)=tb
!
!     row interchange and pivot row reduction in matrix a
   12 do 13 l=lst,mm,m
      ll=l+i
      tb=pivi*a(ll)
      a(ll)=a(l)
   13 a(l)=tb
!
!     save column interchange information
      a(lst)=j
!
!     element reduction and next pivot search
      piv= 0._r8 
      lst=lst+1
      j=0
      do 16 ii=lst,lend
      pivi=-a(ii)
      ist=ii+m
      j=j+1
      do 15 l=ist,mm,m
      ll=l-j
      a(l)=a(l)+pivi*a(ll)
      tb=abs(a(l))
      if(tb-piv)15,15,14
   14 piv=tb
      i=l
   15 continue
      do 16 l=k,nm,m
      ll=l+j
   16 r(ll)=r(ll)+pivi*r(l)
   17 lst=lst+m
!     end of elimination loop
!
!
!     back substitution and back interchange
   18 if(m-1)23,22,19
   19 ist=mm+m
      lst=m+1
      do 21 i=2,m
      ii=lst-i
      ist=ist-lst
      l=ist-m
      l=a(l)+ .5_r8 
      do 21 j=ii,nm,m
      tb=r(j)
      ll=j
      do 20 k=ist,mm,m
      ll=ll+1
   20 tb=tb-a(k)*r(ll)
      k=j+l
      r(j)=r(k)
   21 r(k)=tb
   22 return
!
!
!     error return
   23 ier=-1
      return
      end



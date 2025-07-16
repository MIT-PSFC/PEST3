      MODULE mintrp
      IMPLICIT NONE
	integer nu, ncmax
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rin ! (nu,ncmax)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: zin ! (nu,ncmax)
	integer ntheta0, numc                 ! no of points in contours
	REAL*8 raxis,zaxis                   ! the axis location
        REAL*8 facmid, facedg, facdx0
!															
!  recommended value:  moption=3
	integer moption                    ! control flag
!
!  rotated low order polynomial description of theta line
!  (for internal use only)
!
	REAL*8 r0,z0,r1,z1,r2,z2
	REAL*8 a,b,c,d,xbrk,a2,b2,c2,d2,xmax,zcosa,zsina
	REAL*8 zcosa2,zsina2
!
        REAL*8 xspl(5),yspl(4,5),xpkg(5,4) ! for moption=3
!
!  surface whose intercepts are wanted; last nearby theta interval
!  (for internal use only)
!
	integer ksurf,kth
!
!--------------------------------------------------------------------
        END MODULE mintrp
        MODULE scrunch_inc1
!       begin SCRUNCH_INC1
!       make mpol = mpolDIM + 1 , where mpolDIM is defined in SCRUNCH_INC0
      IMPLICIT NONE
      INTEGER mnd
      INTEGER mpnt
      INTEGER nphi2
      INTEGER nphi
      INTEGER mpol
      INTEGER nznt
      INTEGER ntheta
      INTEGER mpol2
      INTEGER mpol3
      INTEGER mpol4
      INTEGER n2
      INTEGER nfp
      INTEGER nthetax
      INTEGER mpolx
      INTEGER nresets
        parameter(nphi=1,nphi2=1+nphi/2)
	logical ioflagc
      INTEGER, DIMENSION(:), ALLOCATABLE :: mm ! (mpol)
      INTEGER nn(nphi)
      INTEGER, DIMENSION(:), ALLOCATABLE :: m1 ! (mnd)
      INTEGER, DIMENSION(:), ALLOCATABLE :: n1 ! (mnd)
      REAL*8, DIMENSION(:), ALLOCATABLE :: dm1 ! (mpol)
      REAL*8 r0n(nphi)
      REAL*8 z0n(nphi)
      REAL*8 raxis3d(nphi)
      REAL*8 zaxis3d(nphi)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rin3d ! (ntheta,nphi)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: zin3d ! (ntheta,nphi)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: angle ! (ntheta,nphi)
      REAL*8, DIMENSION(:), ALLOCATABLE :: faccon ! (mpol)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: xmpq ! (mpol,4)
      REAL*8 gnorm
      REAL*8 specw
      REAL*8 delt
      REAL*8 deltf
      REAL*8 pexp
      REAL*8 qexp
      REAL*8 dnorm
      REAL*8 elongate
      REAL*8 twopi
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: result ! (nphi2,mpol,4)
      REAL*8, DIMENSION(:), ALLOCATABLE :: xvec ! (n2)
      REAL*8, DIMENSION(:), ALLOCATABLE :: gvec ! (n2)
      REAL*8, DIMENSION(:), ALLOCATABLE :: xdot ! (n2)
      REAL*8, DIMENSION(:), ALLOCATABLE :: xstore ! (n2)
!       end SCRUNCH_INC1
 
      integer :: lunmsg = 6
 
      END MODULE scrunch_inc1
      MODULE scrunch_inc2
!
      IMPLICIT NONE
      LOGICAL lasym
      LOGICAL lzakharov
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rmomb ! (mpol,2)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: zmomb ! (mpol,2)
      REAL*8 rcent
      REAL*8 zcent
      END MODULE scrunch_inc2
 

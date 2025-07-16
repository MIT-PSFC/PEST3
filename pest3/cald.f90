         SUBROUTINE pstcald(a,eps,mb1,mb2,nsing,neg, knode)
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 EPS
      INTEGER MB1
      INTEGER MB2
      INTEGER NSING
      INTEGER NEG
      INTEGER KNODE
      INTEGER N
      INTEGER IKD
      real*8 AD
      INTEGER JIB
      INTEGER M
      INTEGER LOPBND
      INTEGER IJ
      INTEGER JJB
      INTEGER ITOP
      INTEGER IELEM
      INTEGER JRECT
      INTEGER MF
      INTEGER ILEFT


!
!
!     lcm (cald)
!
!***         dimension a(1)
!***
!
!     initialize
!
 COMPLEX*16, DIMENSION(*) 	 :: a
 COMPLEX*16	 	 	 :: diag
 COMPLEX*16	 	 	 :: top
 INTEGER, DIMENSION(100) 	 :: ineg
 COMMON /c2f90ineg/ ineg 
         nsing=0
         n=mb2+1
         if (mb1  ==  0) n=mb2
         ikd=0
         ad=abs(a(1))*eps
!
!     scan over the whole length of a
!
         do  30  jib=2,n
         m=mb1+mb2-jib+2
         diag=a(ikd+1)
         if ( real(diag) <   0._r8 ) then
         neg=neg+1
         ineg(neg) = knode
         end if
!
!     test for zero pivot
!
         if (abs(diag) > ad) go to 5
         if (nsing == -1) go to 5
         nsing=-1
         return
!
!     restriction of loop for not exceeding band matrix
!
   5     continue
         lopbnd=m
!
!     diagonal element before gauss elimination
!
         ij=ikd+m+1
         ad=abs(a(ij))*eps
!
!     sets the row of the transposed left hand side matrix lt
!
         do  20  jjb=2,lopbnd
         itop=ikd+jjb
         top=a(itop)
         a(itop)=a(itop)/diag
         ielem=itop
!
!     gauss rectangular rule going downwards
!
         do  10  jrect=2,jjb
         mf=m-jrect+1
         ielem=ielem+mf
         ileft=ikd+jrect
!***         a(ielem)=a(ielem)-top*a(ileft)
         a(ielem)=a(ielem)-top* conjg( a(ileft) )
!***
  10     continue
  20     continue
         ikd=ikd+m
  30     continue
         if (nsing  ==  0) go to 40
         nsing=0
         return
  40     continue
!
!     last diagonal element
!
         if (mb1  /=  0) return
!***         if (a(ikd+1)  <   0._r8 ) neg=neg+1
         if ( real(a(ikd+1))  <   0._r8 ) neg=neg+1
!***         if (abs(a(ikd+1))  >  ad) return
         if (abs(a(ikd+1))  >  ad) return
         nsing=-1
!
         return
         end




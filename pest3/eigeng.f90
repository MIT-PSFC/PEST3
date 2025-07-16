         SUBROUTINE psteigeng(eps)
!     lcm (eigeng)
 USE pstcom

 USE combla

 USE comivi

 USE l34com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 EPS
      INTEGER IS
      real*8 ALA
      real*8 ALB
      INTEGER I
      INTEGER MJ
      INTEGER NADRES
      INTEGER IX
      INTEGER IN
      real*8 SL
      INTEGER M21
      INTEGER J2


!
!     eigenvalue alam=(x,a*x)/(x,b*x)
!
      is = 1
         ala= 0._r8 
         alb= 0._r8 
         do  10  i=1,m12
         va(i)= 0._r8 
         vb(i)= 0._r8 
  10     continue
         mj=ml
      nadres = 1
      CALL pstblkcpy (xt(1),u(1),mf)
         ix=mf+1
         do  30  in=1,ng
         if (in == ng) mj=m12
      CALL pstblkcpy (xt(ix),u(mf+1),ml)
         ix=ix+ml
    a(1:nlong) = wpot(nadres:nadres+nlong-1)
         CALL pstvax(mj)
         CALL pstalva(sl,va,u,mj)
         ala=ala+sl
    b(1:nlong) = akin(nadres:nadres+nlong-1)
      nadres = nadres + nlong
         CALL pstvbx(mj)
         CALL pstalva(sl,vb,u,mj)
         m21=ml+1

u(m21-ml:m12-ml) = u(m21:m12)
!!$         do 40 j1=m21,m12
!!$         i1=j1-ml
!!$   40    u(i1)=u(j1)

         alb=alb+sl
!
!     eigenvalue vector u
!
         if (nit < nitmax) go to 30
         do  20  j2=1,mj
         if (abs(vb(j2)) > eps) go to 15
         vb(j2)= 0._r8 
         go to 20
  15     continue
!***
! not sure what to do with this line ap  17.10_r8 .95
!***
         vb(j2)=va(j2)/vb(j2)
  20     continue
  30     continue
         alam=ala/alb
         return
!
!       trouble
 7000 CALL psterrmes(ndlt,'eigen',is)
         end



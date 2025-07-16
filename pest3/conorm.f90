         SUBROUTINE pstconorm(eps)
!     lcm (conorm)
 USE pstcom

 USE combla

 USE comivi
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 EPS
      INTEGER IX
      INTEGER IU
      INTEGER MJ
      INTEGER J1
      INTEGER MI
      INTEGER MB
      INTEGER IN
      INTEGER J2
      INTEGER M21


!
         nconv=0
!
!     first mf components
!
         ix=mf+1
         iu=ml+1
      CALL pstblkcpy( xt(1),x(1),mf)
      CALL pstblkcpy(ut(1),u(1),ml)
         mj=ml
         do  10  j1=1,mf
         u(j1)=u(j1)/xnorm
      if (abs(abs(u(j1))-abs(x(j1))) <  eps) nconv=nconv+1
  10     continue
      CALL pstblkcpy(u(1),vt(1),mf)
         mi=ml-mf
         if (mi  <=  0) go to 30
!
!     store ml-mf components of u
!

u(1:mi) = u(1+mf:mi+mf)
!!$         do  20  j1=1,mi
!!$         i2=mf+j1
!!$   20    u(j1)=u(i2)

  30     continue
         mb=mi+ml
!
!     scan over all radial blocks
!
         do  60  in=1,ng
         if (in  ==  ng-1) mb=mb+mf
      CALL pstblkcpy( xt(ix),x(1),ml)
      if (in  /=  ng)CALL pstblkcpy(ut(iu),u(mi+1),m12)
         do  40  j2=1,mj
         u(j2)=u(j2)/xnorm
      if (abs(abs(u(j2))-abs(x(j2))) < eps) nconv=nconv+1
  40     continue
      CALL pstblkcpy(u(1),vt(ix),ml)
         ix=ix+ml
         iu=iu+ml
         if ( in  ==  ng )   go to 60
         if (mb-ml  <=  0) go to 60
!
!     store ml-mf components of u
!
         m21=ml+1

u(m21-ml:mb-ml) = u(m21:mb)
!!$         do  50  j1=m21,mb
!!$         i2=j1-ml
!!$         u(i2)=u(j1)
!!$  50     continue

  60     continue
         return
         end



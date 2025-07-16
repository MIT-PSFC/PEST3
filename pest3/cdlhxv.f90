!
         SUBROUTINE pstcdlhxv
!
 USE pstcom

 USE combla

 USE comivi

 USE l34com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IS
      INTEGER NADRES
      INTEGER IU
      INTEGER MJ
      INTEGER IKD
      real*8 UABS
      INTEGER MJ1
      INTEGER J2
      INTEGER I2
      INTEGER J12
      INTEGER I12
      INTEGER IN
      INTEGER I1

!
!     last triangle
!
 COMPLEX*16	 	 	 :: ui2
!
       is = 1
      nadres = (ng-1) * nlong + 1
!!      CALL pstzrd(nds,a(1),nlong,nadres,is, 7000)
 if ( nds == inpot ) then
    a(1:nlong) = wpot(nadres:nadres+nlong-1)
 else if ( nds == inkin ) then
    a(1:nlong) = akin(nadres:nadres+nlong-1)
 else
    a(1:nlong) = ashift(nadres:nadres+nlong-1)
 endif
         iu=ml*(ng-1)+1
      CALL pstblkcpy( ut(iu),u(1),m12)
         mj=m12
         ikd=nlong
         u(m12)=u(m12)/a(ikd)
!***         uabs=u(m12)**2
         uabs= abs( u(m12) )**2
!***
         xnorm=xnorm+uabs
         mj1=mf-1
         if (mf  <=  1) go to 25
         do  20  j2=1,mj1
         ikd=ikd-j2-1
         i2=m12-j2
         ui2=u(i2)/a(ikd)
         do  10  j12=1,j2
         i12=ikd+j12
         iu=i2+j12
         ui2=ui2-a(i12)*u(iu)
  10     continue
         u(i2)=ui2
!***         uabs=ui2*ui2
         uabs= abs(ui2)**2
!***
         xnorm=xnorm+uabs
  20     continue
  25     continue
!
!     ng radial blocks
!
         do  60  in=1,ng
!
!     scan over ml rows
!
         ikd=mf*(ml-1)+(ml*(ml+1))/2
         do  40  j2=1,ml
         i2=ml-j2+1
         ui2=u(i2)/a(ikd)
         i1=j2+mf-1
         do  30  j12=1,i1
         i12=ikd+j12
         iu=i2+j12
         ui2=ui2-a(i12)*u(iu)
  30     continue
         ikd=ikd-mf-j2-1
         u(i2)=ui2
!***         uabs=ui2*ui2
         uabs= abs(ui2)**2
!**
         xnorm=xnorm+uabs
  40     continue
!
!     write new vector
!
         iu=ml*(ng-in)+1
      CALL pstblkcpy(u(1),ut(iu),mj)
         iu=iu-ml
         mj=ml

         if (in == ng) return
!
!     replace mf vector components
!
         do  50  j2=1,mf
         i2=ml+j2
         u(i2)=u(j2)
  50     continue
!
!     read old vector
!
      CALL pstblkcpy( ut(iu),u(1),ml)
!
!     read matrix blocks
!
      nadres = nadres - nlong
!!      CALL pstzrd(nds,a(1),nlong,nadres,is, 7000)
 if ( nds == inpot ) then
    a(1:nlong) = wpot(nadres:nadres+nlong-1)
 else if ( nds == inkin ) then
    a(1:nlong) = akin(nadres:nadres+nlong-1)
 else
    a(1:nlong) = ashift(nadres:nadres+nlong-1)
 endif
 

  60     continue

         return
!
!       trouble
 7000 CALL psterrmes(ndlt,'cdlhxv',is)
         end



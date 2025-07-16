!-----------------------------------------------------------------
         SUBROUTINE pstvbx(mj)
!
! up-down asym. check aplet  17.10_r8 .95
!
!     lcm (vbx)
 USE pstcom

 USE combla

 USE comivi
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER MJ
      INTEGER IKD
      INTEGER J2
      INTEGER I2
      INTEGER J12
      INTEGER I12


!
!     multiply vb=b*x
!
vb(1:mf) = vb(1+ml:mf+ml)
vb(m11:m12) =  0._r8 

!!$         do  10  i=1,mf
!!$         im2=i+ml
!!$         vb(i)=vb(im2)
!!$  10     continue
!!$         do  20  i=m11,m12
!!$  20     vb(i)= 0._r8 
         ikd=1
         do  40  j2=1,mj
         m=m12-j2+1
         vb(j2)=vb(j2)+b(ikd)*u(j2)
         if (j2  ==  m12) go to 40
         i2=j2+1
         do  30  j12=i2,m12
         i12=ikd+j12-i2+1
         vb(j2)=vb(j2)+b(i12)*u(j12)
         vb(j12)=vb(j12)+b(i12)*u(j2)
  30     continue
         ikd=ikd+m
  40     continue
         return
         end



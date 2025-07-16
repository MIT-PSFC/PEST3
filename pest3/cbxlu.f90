         SUBROUTINE pstcbxlu(idp)
 USE pstcom

 USE combla

 USE comivi

 USE l34com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IDP
      INTEGER IS
      INTEGER NADRES
      INTEGER MJ
      INTEGER IX
      INTEGER IN
      INTEGER IKD
      INTEGER J2
      INTEGER I2
      INTEGER J12
      INTEGER I12
      INTEGER IU
      INTEGER I


!
!
!     initialize
!
      is = 1
      nadres = 1
!
    u(1:m12)= 0._r8 
!
!     read x for first mf components
!
      CALL pstblkcpy( xt(1),x(1),mf)
!
!     scan over all radial blocks
!
      m11 = mf + 1
         mj=ml
         ix=mf+1
         do  70  in=1,ng
!
!     read b and x for a radial block
!

      CALL pstblkcpy( xt(ix),x(mf+1),ml)

         ix=ix+ml
!
      if ( idp  /=  0 )  go to 35
!
!!      CALL pstzrd(ndb,b(1),nlong,nadres,is, 7000)
 b(1:nlong) = akin(nadres:nadres+nlong-1)
!
!     scan over ml rows
!
         if (in  ==  ng) mj=m12
         ikd=1
         do  30  j2=1,mj
         m=m12-j2+1
         u(j2)=u(j2)+b(ikd)*x(j2)
         if (j2  ==  m12) go to 30
         i2=j2+1
         do  20  j12=i2,m12
         i12=ikd+j12-i2+1
!***
! not sure whether these lines need to be modified aplet...
!
         u(j2)=u(j2)+b(i12)*x(j12)
         u(j12)=u(j12)+b(i12)*x(j2)
!***
  20     continue
         ikd=ikd+m
  30     continue
      go to 39
!
   35 continue
!
!      for delta' calculation...
!
      if ( in  ==  ng )  mj = m12
      ikd = 1
      do 37 j2 = 1, mj
      m = m12 - j2 + 1
      u(j2) = u(j2) + x(j2)
      ikd = ikd + m
   37 continue
!
   39 continue
!
!     read a
!
!!      CALL pstzrd(nds,a(1),nlong,nadres,is, 7000)
 if( nds == inpot ) then
    a(1:nlong) = wpot(nadres:nadres+ nlong-1)
 else if ( nds == inkin ) then
    a(1:nlong) = akin(nadres:nadres+ nlong-1)
 else
    a(1:nlong) = ashift(nadres:nadres+ nlong-1)
 end if

      nadres = nadres + nlong

!
!     backsubstitution l*x=u
!
         ikd=1
         do  50  j2=1,mj
         if (j2  ==  m12) go to 50
         i2=j2+1
         m=m12-j2+1
         do  40  j12=i2,m12
         i12=ikd+j12-i2+1
!***         u(j12)=u(j12)-a(i12)*u(j2)
         u(j12)=u(j12)- conjg( a(i12) )*u(j2)
!***
  40     continue
         ikd=ikd+m
  50     continue
!
!     store u
!
         iu=ml*(in-1)+1
      CALL pstblkcpy(u(1),ut(iu),mj)

         if (in ==  ng) return
!
!     set mf components of u and x
!
         do  60  i=1,mf
         im2=i+ml
         u(i)=u(im2)
         x(i)=x(im2)
  60     continue

!
!     ml u-components to zero
!
     u(m11:m12)= 0.0_r8 

  70     continue
!
!
!       trouble
 7000 CALL psterrmes(ndlt,'cbxlu',is )
         end



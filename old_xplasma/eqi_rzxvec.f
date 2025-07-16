      subroutine eqi_rzxvec(ivec,zR,zZ,zphi,ibflg,bvec,gbtensr,
     >   nlist,ivecd,fvals,fgrads)
C
C  transform [B] and grad[B] and list of scalar functions f and grad(f)
C  from (R,Z,phi) coordinates to (x,y,z) coordinates
C
C    x = R*cos(phi)
C    y = R*sin(phi)
C    z = Z
C
C  grad in (R,Z,phi) is (d/dR,d/dZ,(1/R)d/dphi)
C  grad in (x,y,z)   is (d/dx,d/dy,d/dz)
C
      implicit NONE
C
      integer ivec                      ! vector dimensioning
      real*8 zR(ivec),zZ(ivec),zphi(ivec) ! input coordinates
C
      integer ibflg                     ! flag =1 to transform B
      real*8 bvec(3,ivec)               ! B vector to transform
      real*8 gbtensr(3,3,ivec)          ! grad[B] tensor to transform
C
      integer nlist                     ! number of scalar functions (0 OK)
      integer ivecd                     ! output vector dimensioning
      real*8 fvals(ivecd,*)             ! function values
      real*8 fgrads(3,ivecd,*)          ! function gradients
C
      real*8, parameter :: CZERO = 0.0d0
C
C----------------------------------------------
C
      real*8 gradphi(3),gradBR(3),gradBphi(3),gradBZ(3)
      real*8 tmp(3),btmp(3),zcosphi,zsinphi
      integer i,j
C
C----------------------------------------------
C  scalar function values are correct and do not need to be transformed...
C
 
      do i=1,ivec
         zcosphi=cos(zphi(i))
         zsinphi=sin(zphi(i))
C
C  transform gradients of functions
C
         do j=1,nlist
            tmp(3)=fgrads(2,i,j)
            tmp(1)=fgrads(1,i,j)*zcosphi -fgrads(3,i,j)*zsinphi
            tmp(2)=fgrads(1,i,j)*zsinphi +fgrads(3,i,j)*zcosphi
            fgrads(1:3,i,j)=tmp
         enddo
C
         if(ibflg.eq.1) then
C
C  transform B field
C
            btmp(3)=bvec(2,i)
            btmp(1)=bvec(1,i)*zcosphi -bvec(3,i)*zsinphi
            btmp(2)=bvec(1,i)*zsinphi +bvec(3,i)*zcosphi
C
C  transform grad[B]
C
C  grad(Bz):
            gradBZ(3)=gbtensr(2,2,i)
            gradBZ(1)=gbtensr(1,2,i)*zcosphi -gbtensr(3,2,i)*zsinphi
            gradBZ(2)=gbtensr(1,2,i)*zsinphi +gbtensr(3,2,i)*zcosphi
C
C  grad(Bx) and grad(By):  first get grad(phi),grad(BR),grad(BZ)
C                          in x-y-z form
C
            if(zR(i).ne.czero) then
               gradphi(1)= -zsinphi/zR(i)
               gradphi(2)=  zcosphi/zR(i)
            else
               gradphi(1)= CZERO
               gradphi(2)= CZERO
            endif
C
            gradphi(3)=  CZERO
C
            gradBR(3)=gbtensr(2,1,i)
            gradBR(1)=gbtensr(1,1,i)*zcosphi -gbtensr(3,1,i)*zsinphi
            gradBR(2)=gbtensr(1,1,i)*zsinphi +gbtensr(3,1,i)*zcosphi
C
            gradBphi(3)=gbtensr(2,3,i)
            gradBphi(1)=gbtensr(1,3,i)*zcosphi -gbtensr(3,3,i)*zsinphi
            gradBphi(2)=gbtensr(1,3,i)*zsinphi +gbtensr(3,3,i)*zcosphi
C
C  grad(Bz):
            gbtensr(1:3,3,i)=gradBZ
C
C  grad(Bx):
            gbtensr(1:3,1,i)=gradBR*zcosphi - bvec(1,i)*zsinphi*gradphi
     >         - gradBphi*zsinphi - bvec(3,i)*zcosphi*gradphi
C
C  grad(By):
            gbtensr(1:3,2,i)=gradBr*zsinphi + bvec(1,i)*zcosphi*gradphi
     >         + gradbphi*zcosphi - bvec(3,i)*zsinphi*gradphi
C
C  and store vector Bx,By,Bz
C
            bvec(1:3,i)=btmp
C
         endif
      enddo
C
      return
      end

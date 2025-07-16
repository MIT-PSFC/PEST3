      SUBROUTINE scrgctnorm(r,y,n,i,vec,inorm)
!
!  return normal vector, normalized to unit length if inorm=1
!
!  inputs:  (R(j),Y(j)), j going from 1 to n, describes a closed contour,
!           with periodicity, i.e. one can think of
!                  (R(n+1),Y(n+1)) = (R(1),Y(1))
!           i     -- the point where the normal vector is desired
!           inorm -- flags whether output vector is to be normalized
!
!  output:  vec(2) gives the numerically evaluated vector normal to the
!           contour at point i.
!
      IMPLICIT NONE
       REAL*8 	zd ! <gctnorm.f>
       INTEGER 	inorm ! <gctnorm.f>
       REAL*8 	zdely ! <gctnorm.f>
       REAL*8 	zdelr ! <gctnorm.f>
       INTEGER 	imins ! <gctnorm.f>
       INTEGER 	i ! <gctnorm.f>
       INTEGER 	iplus ! <gctnorm.f>
       INTEGER 	n ! <gctnorm.f>
      REAL*8 r(n),y(n)
      REAL*8 vec(2)
!
!-----------------------------
!
      iplus=i+1
      if(iplus.gt.n) iplus=iplus-n
      imins=i-1
      if(imins.lt.1) imins=imins+n
!
      zdelr=r(iplus)-r(imins)
      zdely=y(iplus)-y(imins)
!
      vec(1)=zdely
      vec(2)=-zdelr
!
      if(inorm.eq.1) then
         zd=sqrt(vec(1)**2+vec(2)**2)
         vec(1)=vec(1)/zd
         vec(2)=vec(2)/zd
      endif
!
      return
      end
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"

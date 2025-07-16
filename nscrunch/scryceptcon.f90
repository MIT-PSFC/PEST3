!--------------------------
      REAL*8 function scrYceptcon(x,Rtest,Ztest)
!
!  evaluate y(x):  theta line fit in some reference frame
!  also return the (x,y) rotated back to (R,Z) space.
!
      USE mintrp
      IMPLICIT NONE
       REAL*8 	ztest ! <scrYceptcon.f>
       REAL*8 	rtest ! <scrYceptcon.f>
       REAL*8 	zdydx ! <scrYceptcon.f>
       REAL*8 	zdx ! <scrYceptcon.f>
       REAL*8 	zref ! <scrYceptcon.f>
       REAL*8 	rref ! <scrYceptcon.f>
       REAL*8 	zsin ! <scrYceptcon.f>
       REAL*8 	zcos ! <scrYceptcon.f>
       REAL*8 	y ! <scrYceptcon.f>
       REAL*8 	zxx ! <scrYceptcon.f>
       REAL*8 	x ! <scrYceptcon.f>
       REAL*8 	zx ! <scrYceptcon.f>
!
!  for spvec (moption=3)...
!
       integer ict(3),iwarn,ier
       data ict/1,0,0/
!
!---------------------------------------------
!
      zx=max(xspl(1),min(xmax,x))
      if(zx.gt.xbrk) then
         zxx=zx-xbrk
         y=d2+zxx*(c2+zxx*(b2+zxx*a2))
         zcos=zcosa2
         zsin=zsina2
         rref=r1
         zref=z1
         zx=zxx
      else
         if(moption.eq.2) then
            y=d+zx*(c+zx*(b+zx*a))
         else
            call r8spvec(ict,1,zx,1,y,5,xpkg,yspl,iwarn,ier)
            if(ier.ne.0) call errmsg_exit('?scryceptcon: spvec error.')
         endif
         zcos=zcosa
         zsin=zsina
         rref=r0
         zref=z0
      endif
!
!  linear extrapolation
!
      if(x.gt.xmax) then
         zdx=x-xmax
         if(xmax.gt.xbrk) then
            zx=xmax-xbrk
            zdydx=c2+2.0D0*zx*(b2+1.5D0*zx*a2)
         else
            zx=x
            zdydx=c+2.0D0*xmax*(b+1.5D0*xmax*a)
         endif
      else
         zdx=0.0D0
         zdydx=0.0D0
      endif
!
      y=y+zdx*zdydx
      scrYceptcon=y
!
!  transform to (R,Z)
!
      Rtest = rref + zcos*(zx+zdx) -zsin*y
      Ztest = zref + zsin*(zx+zdx) +zcos*y
!
      return
      end
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"

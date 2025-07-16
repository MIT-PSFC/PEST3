      SUBROUTINE scrfsetup(ith,icc,rwork,zwork,inpts,idim1,zxmax,ier)
!
!
      USE mintrp
!
!  set up theta line fit function for theta line #ith
!  data points scatterred in rwork and zwork arrays (from
!  previous scrunch operations by caller)
!
!  location specified by icc(...) and ith -- depending on
!  algorithm selected by moption
!
!  experimental code...
!
      IMPLICIT NONE
       REAL*8 	zxmax ! <fsetup.f>
       REAL*8 	dfdx0 ! <fsetup.f>
       REAL*8 	zsinb ! <fsetup.f>
       REAL*8 	zcosb ! <fsetup.f>
       REAL*8 	dfdx1 ! <fsetup.f>
       REAL*8 	zs0 ! <fsetup.f>
       REAL*8 	zadjust ! <fsetup.f>
       REAL*8 	zdz2 ! <fsetup.f>
       REAL*8 	zdr2 ! <fsetup.f>
       REAL*8 	zdz1 ! <fsetup.f>
       REAL*8 	zdr1 ! <fsetup.f>
       REAL*8 	zs1 ! <fsetup.f>
       REAL*8 	zdy ! <fsetup.f>
       REAL*8 	zdx ! <fsetup.f>
       REAL*8 	zdz ! <fsetup.f>
       REAL*8 	zdr ! <fsetup.f>
       REAL*8 	zxd ! <fsetup.f>
       REAL*8 	zdelz ! <fsetup.f>
       REAL*8 	zdelr ! <fsetup.f>
       REAL*8 	y2 ! <fsetup.f>
       REAL*8 	x2 ! <fsetup.f>
       REAL*8 	znorm ! <fsetup.f>
       INTEGER 	ith ! <fsetup.f>
       INTEGER 	ic2 ! <fsetup.f>
       INTEGER 	ics ! <fsetup.f>
       INTEGER 	ic1 ! <fsetup.f>
       INTEGER 	ic0 ! <fsetup.f>
       INTEGER 	idim1 ! <fsetup.f>
       INTEGER 	inpts ! <fsetup.f>
      REAL*8 rwork(inpts+1,idim1),zwork(inpts+1,idim1),zxax1,zwk(5)
      integer icc(6),i,ilinx,ier
      REAL*8 enorm(2)
      REAL*8, parameter :: ZERO = 0.0D0
!
      integer lunmsg
!
!--------------------------------
!
      ier = 0  ! assume normal status initially...
!
      call scrgetlun(lunmsg)
!
      ic0=icc(1)
      if(moption.eq.2) then
         ic1=icc(2)
         ics=icc(3)
         ic2=icc(4)
      else if(moption.eq.3) then
         ic1=icc(4)
         ics=icc(5)
         ic2=icc(6)
      endif
!
      r0=rwork(ith,ic0)
      z0=zwork(ith,ic0)
!
      r2=rwork(ith,ic2)
      z2=zwork(ith,ic2)
!
!  normal vector
!
      enorm(1)=rwork(ith,ic0)-raxis
      enorm(2)=zwork(ith,ic0)-zaxis
      znorm=sqrt(enorm(1)**2+enorm(2)**2)
      zxax1=znorm
      enorm(1)=enorm(1)/znorm
      enorm(2)=enorm(2)/znorm
!
      if((moption.eq.2).or.(moption.eq.3)) then
!
!  will use a two piece curve, C(1), with "x" defined along two legs
!  of a triangle
!    for moption.eq.3 the piece nearer the axis is a spline.
!
         r1=rwork(ith,ic1)
         z1=zwork(ith,ic1)
!
         zdelr=r1-r0
         zdelz=z1-z0
         xbrk=sqrt(zdelr*zdelr+zdelz*zdelz)
         zcosa=zdelr/xbrk
         zsina=zdelz/xbrk
!
         zdelr=r2-r1
         zdelz=z2-z1
         zxd=sqrt(zdelr*zdelr+zdelz*zdelz)
         zcosa2=zdelr/zxd
         zsina2=zdelz/zxd
!
         xmax=xbrk+zxd
!
!  segment x = xbrk to xmax -- use low order fit
!
         zdr=r2-rwork(ith,ics)
         zdz=z2-zwork(ith,ics)
         zdx= zcosa2*zdr +zsina2*zdz
         if(zdx.le.ZERO) then
           CALL scrabwarn(lunmsg,'?momintrp: bdy slope ^ ref angle too large')
           ier=99
           return
         endif
         zdy=-zsina2*zdr +zcosa2*zdz
!
         zs1=zdy/zdx
!
         zdy=-zdy
!
!  direction in (R,Z) space if (xbrk,xmax) fit is parabalic
!
         zdr1= zcosa2*zdx -zsina2*zdy
         zdz1= zsina2*zdx +zcosa2*zdy
         znorm=sqrt(zdr1**2+zdz1**2)
         zdr1=zdr1/znorm
         zdz1=zdz1/znorm
!
!  direction that parallels the segment (ctr) (R0,Z0) -- (R2,Z2) (bdy)
!
         zdr2=(r2-r0)
         zdz2=(z2-z0)
         znorm=sqrt(zdr2**2+zdz2**2)
         zdr2=zdr2/znorm
         zdz2=zdz2/znorm
!
!  set zadjust=ZERO for parabolic fit if (xbrk,xmax) segment
!  set zadjust=1.0 to use (R0,Z0) -- (R2,Z2) segment for direction
!
         zadjust=0.5D0
         zdr=(1.0D0-zadjust)*zdr1 + zadjust*zdr2
         zdz=(1.0D0-zadjust)*zdz1 + zadjust*zdz2
!
!  back to ref from for (xbrk,xmax) segment
!
         zdx= zcosa2*zdr +zsina2*zdz
         if(zdx.le.ZERO) then
           CALL scrabwarn(lunmsg,'?momintrp: brk slope x ref angle too large')
           ier=99
           return
         endif
         zdy=-zsina2*zdr +zcosa2*zdz
!
         zs0=zdy/zdx
!
!  here is the cubic fit  (zadjust=ZERO => zs0=-zs1 & parabolic fit results)
!
         a2=(zs0+zs1)/(zxd*zxd)
         b2=-(2*zs0+zs1)/zxd
         c2=zs0
         d2=ZERO
!
!  use same (zdr,zdz) in the (0,xbrk) segment ref. frame.
!
         zdx= zcosa*zdr +zsina*zdz
         if(zdx.le.ZERO) then
           CALL scrabwarn(lunmsg,'?momintrp: brk slope ^ ref angle too large')
           ier=99
           return
         endif
         zdy=-zsina*zdr +zcosa*zdz
!
         dfdx1=zdy/zdx
!
!  slope of theta line at ctr -- segment frame
!
         if(moption.eq.2) then
            zcosb= zcosa*enorm(1) +zsina*enorm(2)
            zsinb=-zsina*enorm(1) +zcosa*enorm(2)
            if(zcosb.le.ZERO) then
               CALL scrabwarn(lunmsg,'?momintrp:  normal dot ref is .le. 0.0')
               ier=99
               return
            endif
            dfdx0=zsinb/zcosb
!
!  use cubic fit -- connecting segment reference frame
!
            a=(dfdx1+dfdx0)/(xbrk*xbrk)
            b=-(2*dfdx0+dfdx1)/xbrk
            c=dfdx0
            d=ZERO
            xspl=ZERO
         else
!
!  moption=3 -- use spline for interior region
!
            xspl(2)=ZERO;   yspl(1,2)=ZERO   ! icc(1) reference pt.
            xspl(5)=xbrk;   yspl(1,5)=ZERO   ! icc(4) reference pt.
!
! axis
!
            xspl(1)=  -zxax1*( zcosa*enorm(1) +zsina*enorm(2))
            yspl(1,1)=-zxax1*(-zsina*enorm(1) +zcosa*enorm(2))
!
! intermediate surfaces
!
            do i=2,3
               zdr=rwork(ith,icc(i))-r0
               zdz=zwork(ith,icc(i))-z0
               xspl(i+1)=   zcosa*zdr +zsina*zdz
               yspl(1,i+1)=-zsina*zdr +zcosa*zdz
            enddo
!
            if((xspl(2)-xspl(1)).lt.1.0D-3*(xspl(3)-xspl(2))) then
               ier=88  ! fit will not work - too much change in direction
            else
!
! setup spline (& uneven x axis)
!
               call r8cspline(xspl,5,yspl,0,ZERO,1,dfdx1,zwk,5,ilinx,ier)
               if(ier.ne.0) then
                  call scrabwarn(lunmsg,'?fsetup: cspline failure.')
                  ier=99
                  return
               endif
               call r8genxpkg(5,xspl,xpkg,0,0,0,ZERO,2,ier)
               if(ier.ne.0) then
                  call scrabwarn(lunmsg,'?fsetup: genxpkg failure.')
                  ier=99
                  return
               endif
            endif
         endif  ! (moption.eq.3)
!
      else
         CALL scrabwarn(lunmsg,'?fsetup:  illegal moption value.')
         ier=99
         return
      endif
!
      zxmax=xmax
!
      return
      end

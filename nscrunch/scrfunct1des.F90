#include "fpreproc/library_names.h"
#include "fpreproc/blas_names.h"
      SUBROUTINE scrfunct1des(rmc,rms,zmc,zms,xpts,grc,grs, &
     &  gzc,gzs,gpts,fsq,r10,mmax)
      USE scrunch_inc1
      USE scrunch_inc2
      IMPLICIT NONE
       REAL*8 	fac ! <funct1des.f>
       INTEGER 	i ! <funct1des.f>
       REAL*8 	Sdot ! <funct1des.f>
       REAL*8 	tcon ! <funct1des.f>
       REAL*8 	r10 ! <funct1des.f>
       REAL*8 	fsq ! <funct1des.f>
       REAL*8 	t1 ! <funct1des.f>
       REAL*8 	t2 ! <funct1des.f>
       INTEGER 	mmax ! <funct1des.f>
       INTEGER 	m ! <funct1des.f>
       INTEGER 	l ! <funct1des.f>
       INTEGER 	j ! <funct1des.f>
       REAL*8 	denom ! <funct1des.f>
!**     jpj@jet 28-feb-1996 modif. following declarations
!       REAL*8 yc(mpol),ys(mpol),gcon(ntheta),gtt(ntheta),r1(ntheta),
!    >  z1(ntheta),rt1(ntheta),zt1(ntheta),arg(ntheta),rcon(ntheta),
!    >  zcon(ntheta),cosa(ntheta,mpol),sina(ntheta,mpol),rmc(1),rms(1),
!    >  zmc(1),zms(1),xpts(1),grc(1),grs(1),gzc(1),gzs(1),gpts(1)
!**
        REAL*8 yc(mpol),ys(mpol),gcon(ntheta),gtt(ntheta),r1(ntheta), &
     &  z1(ntheta),rt1(ntheta),zt1(ntheta),arg(ntheta),rcon(ntheta), &
     &  zcon(ntheta),cosa(ntheta,mpol),sina(ntheta,mpol),rmc(*),rms(*), &
     &  zmc(*),zms(*),xpts(*),grc(*),grs(*),gzc(*),gzs(*),gpts(*)
!
!
!	FORCES: dW/dRmn = SUM(i)[ R(i) - RIN3D(i)]*cos(m*u[i]) ....
!		dW/dZmn = SUM(i)[ Z(i) - ZIN3D(i)]*sin(m*u[i]) ....
!		dW/du[i]=    rt(i)*[R(i)-RIN3D(i)] + zt(i)*[Z(i) - ZIN3D(i)]
!	THE NORM ON THE ANGLE FORCE (dW/du[i]) FOLLOWS FROM NEWTON S
!	LAW AND IS APPROXIMATELY GTT = RT**2 + ZT**2
!
        denom = 0.D0
        gnorm = 0.D0
        specw = 0.D0
        do 10 j = 1,nthetax
        r1(j) = -rin3d(j,1)
        z1(j) = -zin3d(j,1)
        rcon(j) = 0.D0
        zcon(j) = 0.D0
        rt1(j)  = 0.D0
 10     zt1(j)  = 0.D0
        do 20 l = 1,n2
 20     grc(l) = 0.D0
!**********************************************************************
!       COMPUTE SPECTRAL WIDTH OF CURVE
!**********************************************************************
        do 30 m = 2,mmax
        t2 = (rmc(m)*rmc(m)+zmc(m)*zmc(m)+ &
     &      rms(m)*rms(m)+zms(m)*zms(m))*xmpq(m,1)
        denom = denom+t2
 30     specw = specw+xmpq(m,2)*t2
        specw = specw/denom
!**********************************************************************
!       COMPUTE CURVE AND CONSTRAINT FORCES
!**********************************************************************
        do 40 m = 1,mmax
        do 40 l = 1,nthetax
        arg(l)    = dm1(m) * xpts(l)
        cosa(l,m) = cos(arg(l))
        sina(l,m) = sin(arg(l))
        gtt(l)  = rmc(m)*cosa(l,m) + rms(m)*sina(l,m)
        gcon(l) = zmc(m)*cosa(l,m) + zms(m)*sina(l,m)
        r1(l)   = r1(l)  + gtt(l)
        z1(l)   = z1(l)  + gcon(l)
        rt1(l)  = rt1(l) + dm1(m)*(rms(m)*cosa(l,m) - rmc(m)*sina(l,m))
        zt1(l)  = zt1(l) + dm1(m)*(zms(m)*cosa(l,m) - zmc(m)*sina(l,m))
        rcon(l) = rcon(l) + xmpq(m,4)*gtt(l)
        zcon(l) = zcon(l) + xmpq(m,4)*gcon(l)
 40     continue
        do 50 l = 1,nthetax
        gtt(l)  = rt1(l)*rt1(l) + zt1(l)*zt1(l)
        gpts(l) = r1(l)*rt1(l)+z1(l)*zt1(l)
        gcon(l) = rcon(l)*rt1(l) + zcon(l)*zt1(l)
 50     continue
	t1 = -1.0D0
	do l = 1,nthetax
	t1 = max(t1,gtt(l))
	enddo
!**********************************************************************
!       COMPUTE MEAN-SQUARE DEVIATION BETWEEN POINTS AND FIT
!**********************************************************************
        fsq = 0.D0
        do 60 l = 1,nthetax
	gpts(l) = gpts(l)/t1
 60     fsq = fsq + .5D0*(r1(l)**2 + z1(l)**2)
        fsq = sqrt(dnorm*fsq)/r10
!**********************************************************************
!       FILTER CONSTRAINT FORCE TO REMOVE ALIASING
!**********************************************************************
        tcon = 1.D0/(t1*sqrt(1.D0+xmpq(mmax,3)))
        do 70 m = 2,mmax-1
        yc(m) = Sdot(nthetax,cosa(1,m),1,gcon,1)*faccon(m)*tcon
 70     ys(m) = Sdot(nthetax,sina(1,m),1,gcon,1)*faccon(m)*tcon
        do 80 l = 1,nthetax
 80     gcon(l) = 0.D0
        do 90 m = 2,mmax-1
        do 90 l = 1,nthetax
 90     gcon(l) = gcon(l) + yc(m)*cosa(l,m) + ys(m)*sina(l,m)
!**********************************************************************
!       ADD CURVE AND CONSTRAINT FORCES
!**********************************************************************
        do 100 m = 1,mmax
        do 100 l = 1,nthetax
        rcon(l) = r1(l)  + gcon(l)*rt1(l)*xmpq(m,3)
        zcon(l) = z1(l)  + gcon(l)*zt1(l)*xmpq(m,3)
        grc(m)  = grc(m) + cosa(l,m)*rcon(l)
	if ((lasym.and. .not.lzakharov).or.m.eq.1) &
     &	    gzc(m)  = gzc(m) + cosa(l,m)*zcon(l)
        if (lasym) grs(m)  = grs(m) + sina(l,m)*rcon(l)
	if (.not.lzakharov .or. m.le.2) gzs(m)  = gzs(m) + &
     &	    sina(l,m)*zcon(l)
  100   continue
!**********************************************************************
!       COMPUTE m=1 CONSTRAINT (ZC(m=1)=RS(m=1)) to the Forces
!**********************************************************************
	if(.not.lzakharov) then
	  gzc(2) = grs(2) + gzc(2)
	  grs(2) = gzc(2)
	  grc(1) = 0.5D0*grc(1)
	end if
 
	do i = 1,mpol4
	grc(i) = dnorm*grc(i)
	enddo
        do 110 j=1,mmax
 110    gnorm = gnorm + grc(j)*grc(j) + gzc(j)*gzc(j) &
     &                + grs(j)*grs(j) + gzs(j)*gzs(j)
        gnorm = gnorm/r10**2
        do 120 j = 1,nthetax
 120    gnorm = gnorm + dnorm*gpts(j)*gpts(j)
        do 130 m = 2,mmax
	fac = 1.D0/(1.D0+ tcon*xmpq(m,3))
        grc(m) = grc(m)*fac
        gzc(m) = gzc(m)*fac
        grs(m) = grs(m)*fac
 130    gzs(m) = gzs(m)*fac
        return
        end
 
 
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"
! 29Oct2006 fgtok

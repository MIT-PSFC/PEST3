#include "fpreproc/f77_dcomplx.h"
      SUBROUTINE scrgetangldes(rval,zval,xpts,rcenter,zcenter)
      USE scrunch_inc1
      IMPLICIT NONE
       REAL*8 	del_ang ! <getangldes.f>
       REAL*8 	denom ! <getangldes.f>
       REAL*8 	dnum ! <getangldes.f>
       REAL*8 	yc ! <getangldes.f>
       REAL*8 	xc ! <getangldes.f>
       INTEGER 	iterate ! <getangldes.f>
       INTEGER 	j ! <getangldes.f>
       INTEGER 	i ! <getangldes.f>
        REAL*8 rval(ntheta,*),zval(ntheta,*),xpts(ntheta,*), &
     &  rcenter(*),zcenter(*),rcos(nphi),rsin(nphi), &
     &  zcos(nphi),zsin(nphi),phi_ang(nphi)
!**********************************************************************
!       Compute angle offset consistent with constraint Z1n = Z1,-n
!       Note: This is done iteratively, since elongation is unknown
!**********************************************************************
        do 10 i = 1,nphi
        do 10 j = 1,nthetax
 10     xpts(j,i) = twopi*(j-1)/AREAL(nthetax)
        do 100 iterate = 1,5
        do 20 i = 1,nphi
        rcos(i) = 0.D0
        rsin(i) = 0.D0
        zcos(i) = 0.D0
        zsin(i) = 0.D0
        do 20 j = 1,nthetax
        xc = rval(j,i) - rcenter(i)
        yc = zval(j,i) - zcenter(i)
        rcos(i) = rcos(i) + cos(xpts(j,i))*xc
        rsin(i) = rsin(i) + sin(xpts(j,i))*xc
        zcos(i) = zcos(i) + cos(xpts(j,i))*yc
        zsin(i) = zsin(i) + sin(xpts(j,i))*yc
 20     continue
!**********************************************************************
!       Compute new angles starting from offset phi_ang(i)
!**********************************************************************
	dnum = 0.D0
	denom= 0.D0
	do i = 1,nphi
	dnum = dnum + zsin(i)
	denom= denom+ rcos(i)
	enddo
        elongate = dnum/denom
        del_ang = 0.D0
        do 30 i = 1,nphi
        phi_ang(i) = atan2(elongate*zcos(i)-rsin(i), &
     &                     elongate*zsin(i)+rcos(i))
        del_ang = max(del_ang,abs(phi_ang(i)))
        do 30 j = 1,nthetax
 30     xpts(j,i) = xpts(j,i) + phi_ang(i)
        if(del_ang.lt.0.02D0)go to 40
 100    continue
 40     if (ioflagc) write(lunmsg,*)' Average elongation = ',elongate
        return
        end
 
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"
! 29Oct2006 fgtok

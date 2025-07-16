#include "fpreproc/library_names.h"
#include "fpreproc/f77_dcomplx.h"
#include "fpreproc/blas_names.h"
      SUBROUTINE scramplituddes(rcenter,zcenter,angin,rmc,rms, &
     &  zmc,zms,xpts,xin,yin,ier)
      USE scrunch_inc1
      USE scrunch_inc2
      IMPLICIT NONE
       INTEGER 	ier ! <amplituddes.f>
       INTEGER 	m ! <amplituddes.f>
       REAL*8 	yi ! <amplituddes.f>
       REAL*8 	xi ! <amplituddes.f>
       REAL*8 	arg ! <amplituddes.f>
       INTEGER 	j ! <amplituddes.f>
       REAL*8 	xmult ! <amplituddes.f>
       REAL*8 	zcenter ! <amplituddes.f>
       REAL*8 	rcenter ! <amplituddes.f>
      REAL*8 rmc(*),rms(*),zmc(*),zms(*), &
     &       angin(*),xpts(*),xin(*),yin(*)
!*****************************************************************
!       This subroutine assigns initial guesses for angles and
!       Fourier mode amplitudes to the appropriate components of
!       the xvec array ... use Fourier Analysis to project out
!       each coefficient, below.
!*****************************************************************
        rmc(1) = rcenter
        zmc(1) = zcenter
        CALL Scopy(ntheta,angin,1,xpts,1)
        xmult = 2.D0/AREAL(nthetax)
        do 10 j = 1,nthetax
        arg = angin(j)
        xi = xmult*(xin(j) - rcenter)
        yi = xmult*(yin(j) - zcenter)
        do 10 m = 2,mpolx
        rmc(m) = rmc(m) + cos((m-1)*arg)*xi
        if (lasym) rms(m) = rms(m) + sin((m-1)*arg)*xi
	if (lasym.and. .not.lzakharov) &
     &	    zmc(m) = zmc(m) + cos((m-1)*arg)*yi
	if (.not.lzakharov.or. m.eq.2) zms(m) = zms(m) + &
     &	    sin((m-1)*arg)*yi
 10   continue
!
!	INITIALIZE RMS = ZMC FOR M=1
!
	rms(2) = 0.5D0*(rms(2) + zmc(2))
	if (.not. lzakharov) zmc(2) = rms(2)
        return
        end
 
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"
! 29Oct2006 fgtok
! 29Oct2006 fgtok

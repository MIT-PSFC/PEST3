#include "fpreproc/f77_dcomplx.h"
      SUBROUTINE scringuessdes (ier)
      USE scrunch_inc1
      IMPLICIT NONE
       INTEGER 	ier ! <inguessdes.f>
       INTEGER 	j ! <inguessdes.f>
       INTEGER 	i ! <inguessdes.f>
!**********************************************************************
!       This subroutine obtains initial guesses for the centroid at each
!       toroidal plane. By default, the polar axis is set equal to this
!       centroid.  It is imperative that the polar axis be enclosed by
!       the surface (otherwise the computed theta angles will not span
!       [0,2pi]). For certain complex cross-sections, this simple estimate
!       of the polar axis may fail, and the user must supply values for
!       raxis3d(i),zaxis3d(i).  In addition, for non-starlike domains, the
!       points along the curve must be monitonically increasing as a
!       function of the arclength from any fixed point on the curve. This
!       ordering is attempted by the subroutine ORDER, but may fail in
!       general if too few theta points are used.
!**********************************************************************
!       COMPUTE CENTROID
!**********************************************************************
        do 10 i = 1,nphi
        r0n(i) = 0.0D0
        z0n(i) = 0.0D0
        do 20 j = 1,nthetax
        r0n(i) = r0n(i)+rin3d(j,i)/AREAL(nthetax)
 20     z0n(i) = z0n(i)+zin3d(j,i)/AREAL(nthetax)
        raxis3d(i) = r0n(i)
 10     zaxis3d(i) = z0n(i)
!**********************************************************************
!       ORDER POINTS ON FLUX SURFACE AT EACH TOROIDAL ANGLES
!**********************************************************************
        if (ioflagc) write(lunmsg,*)'ORDERING SURFACE POINTS'
        DO i = 1,nphi
           CALL scrorderdes(rin3d(1,i),zin3d(1,i), &
     &                   raxis3d(i),zaxis3d(i),ier)
        END DO
!**********************************************************************
!       COMPUTE OPTIMUM ANGLE OFFSETS FOR M=1 MODES
!**********************************************************************
      CALL scrgetangldes(rin3d,zin3d,angle,r0n,z0n)
        return
        end
 
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"
! 29Oct2006 fgtok

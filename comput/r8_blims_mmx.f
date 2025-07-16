C-----------------------------------------------
c  utility subroutine for blims -- find min/max singularity
c  use parabolic fit from 3 data pts
c
      subroutine r8_blims_mmx(zth0,zdth,zbm,zb0,zbp,zbsing,zthsing)
c
c input:
c  zth0 -- theta of ctrmost data pt (zb0)
c  zdth -- delta(theta) from zth0 to zbm (-) and zbp (+) data points
c  zbm,zb0,zbp -- sequence of three B field data pts
c
c  expected:  either zbm.ge.zb0.and.zbp.ge.zb0
c               or   zbm.le.zb0.and.zbp.le.zb0
c
c output:
c  zbsing -- min (or max) B based on parabolic fit
c  zthsing -- theta location (normalized to range (-pi,pi)) of the
c            singularity.
c
c fit:
c
c  th from zth0-zdth to zth0+zdth
c
c  f(th) = zb0 + (zbp-zbm)/(2*zdth) * (th-zth0)
c              + (zbp+zbm-2*zb0)/(2*zdth**2) * (th-zth0)**2
c
c  f'(th) = (zbp-zbm)/(2*zdth) + (zbp+zbm-2*zb0)/(zdth**2) * (th-zth0)
c
c  solve f'(th)=0 to find zthsing; evaluate f(zthsing).
c
c-------------------
c  1st check for degenerate case
c
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      REAL*8 zdth,zbm,zb0,zbp,zbsing,zthsing,zth0,xpi,twopi,zanum
      REAL*8 za,zb,zdelth
!============
      data xpi/3.1415926535897931D+00/
      data twopi/6.2831853071795862D+00/
c
      zanum=(zbp+zbm-2*zb0)
      if(zanum.eq.0.0D0) then
         zbsing=zb0
         zthsing=zth0
         go to 10
      endif
c
      za=zanum/(2*zdth**2)
      zb=(zbp-zbm)/(2*zdth)
c
      zdelth = -zdth*zb/(2*za)
      zthsing= zth0 + zdelth
c
      zbsing = zb0 + zb*zdelth + za*zdelth**2
c
c  standardize zthsing range
c
 10   continue
      if(zthsing.lt.-xpi) then
         zthsing=min(xpi,zthsing+twopi)
         go to 10
      else if(zthsing.gt.xpi) then
         zthsing=max(-xpi,zthsing-twopi)
         go to 10
      endif
c
      return
      end

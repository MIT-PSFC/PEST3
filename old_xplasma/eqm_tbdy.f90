subroutine eqm_tbdy(ilines,rl,zl,thl,icircs,rc,zc,rad,ierr)

  !  set up a TRANSP-style machine boundary specification.  This consists
  !  two lists:  a list of circles (R,Z) & radius, and a list of infinite
  !  lines (R,Z) and inclination angle, thl (degrees).

  !  for each circle, the plasma is excluded from the space inside
  !  the circle (sign=-1.) or is excluded from the space outside
  !  (sign=+1.).  The sign is determined by whether the plasma magnetic
  !  axis point (R0,Z0) is inside or outside the circle.

  !  for each infinite line, a start point and orientation angle are
  !  given.  the angle (thl) is in degrees, with convention:

  !    thl(i)=0.0d0 -- horizontal orientation
  !    thl(i)=45.0d0 -- line goes up & to right (+R,+Z) & down to left (-R,-Z)
  !    thl(i)=90.0d0 -- vertical orientation

  !  **CAUTION TRANSP USERS** in TRANSP limiters theta=45.0d is up and to
  !  the left (-R,+Z). -- this is a different convention ***

  !  each line partitions space into two regions, from one of which the
  !  plasma is excluded.  Using (R0,Z0) to test, a sign value is assigned
  !  which determines which half plane is excluded.

  !  the result:  excluded space = union of all excluded spaces of
  !                 individual circle/line "limiters".

  !               included space = all space - excluded space

  !  the plasma boundary must lie entirely inside the included space.
  !  this routine test for this-- if this test fails the error flag is
  !  set.

  !  NOTE -- this routine also makes a "sanity limiter" with ample space
  !  around the plasma equilibrium boundary (which is presumed known).
  !  Because of this, it would be well to call this routine after each
  !  equilibrium update, if it is used at all.
  !-----------------------

  use xplasma_obj_instance
  use eq_module
 
  implicit NONE

  !  input:

  integer ilines                    ! # of line limiters (can be zero)
  real*8 rl(ilines),zl(ilines)      ! (R,Z) starting point of line(s)
  real*8 thl(ilines)                ! in DEGREES, line orientation

  integer icircs                    ! # of circle limiters (can be zero)
  real*8 rc(icircs),zc(icircs)      ! (R,Z) center of circle(s)
  real*8 rad(icircs)                ! radius of circle(s)

  !  output:

  integer ierr                      ! completion code, 0=OK
      
  !-----------------------
  integer iwarn
  real*8 zrmin,zrmax,zzmin,zzmax,zdelr,zdelz
  !-----------------------

  call eq_glimrz(cone,zrmin,zrmax,zzmin,zzmax,ierr)
  if(ierr.ne.0) return

  if((ilines.le.0).and.(icircs.le.0)) then
     call eq_errmsg(' ?eqm_tbdy:  no limiters specified.')
     ierr=2
     return
  endif

  !  4 line-type "sanity check" limiters:
  !  horizontal, vertical from expanded (zrmin,zzmin), &
  !  horizontal, vertical from expanded (zrmax,zzmax)

  zdelr=chalf*(zrmax-zrmin)
  zdelz=chalf*(zzmax-zzmin)

  zrmin=max(chalf*zrmin,zrmin-zdelr)
  zrmax=zrmax+zdelr

  zzmax=zzmax+min(zdelr,zdelz)
  zzmin=zzmin-min(zdelr,zdelz)

  !-----------------------------------------

  call xplasma_mklim_special(s,iwarn,ierr, &
       Rmin=zRmin, Rmax=zRmax, Zmin=zZmin, Zmax=zZmax, &
       nlines=ilines, Rl=rl, Zl=zl, thl=thl, &
       ncircs=icircs, Rc=rc, Zc=zc, rad=rad)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_tbdy: xplasma_mklim_special error:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eqm_tbdy

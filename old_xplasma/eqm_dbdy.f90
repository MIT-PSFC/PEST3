subroutine eqm_dbdy(zdist,zrmini,zrmaxi,zzmini,zzmaxi,ierr)

  !  f77 xplasma legacy interface.

  !  set up a limiter which is a combination of a box and
  !  a fixed distance set back from the plasma last closed flux surface

  !  The equilibrium must already be known to xplasma.

  !  A sanity box limiter is imposed -- one midplane half-width away from the
  !  R and Z min and max of the plasma equilibrium boundary.  The user
  !  supplied box limiter (zrmini, zrmaxi, zzmini, zzmaxi) is constrained
  !  by this.
  !-----------------------

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  !  input:

  real*8 zdist                      ! setback distance of limiter
  real*8 zrmini,zrmaxi              ! R limits of enclosing box
  real*8 zzmini,zzmaxi              ! Z limits of enclosing box

  !  output:

  integer ierr                      ! completion code, 0=OK

  !-----------------------
  integer :: iwarn   ! ignored...
  real*8 zrmin,zrmax,zzmin,zzmax,zdelr,zdelz
  real*8 zrminu,zrmaxu,zzminu,zzmaxu
  !-----------------------

  call eq_glimrz(cone,zrmin,zrmax,zzmin,zzmax,ierr)
  if(ierr.ne.0) return

  !  4 line-type "sanity check" limiters:
  !  horizontal, vertical from expanded (zrmin,zzmin), &
  !  horizontal, vertical from expanded (zrmax,zzmax)

  zdelr=chalf*(zrmax-zrmin)
  zdelz=chalf*(zzmax-zzmin)

  zrmin=max(chalf*zrmin,zrmin-zdelr)
  zrmax=zrmax+zdelr

  zzmax=zzmax+min(zdelr,zdelz)
  zzmin=zzmin-min(zdelr,zdelz)

  !  The RZ box limiter used is the input limiter constrained by the
  !  sanity limiter.  This prevents the limiter from being drawn "too far"
  !  from the plasma.  It is still the caller's responsibility to make sure
  !  that the limiter does not cut into the plasma itself.

  zrminu=max(zrmin,min(zrmax,zrmini))
  zrmaxu=max(zrmin,min(zrmax,zrmaxi))

  zzminu=max(zzmin,min(zzmax,zzmini))
  zzmaxu=max(zzmin,min(zzmax,zzmaxi))

  call xplasma_mklim_special(s,iwarn,ierr, &
       dist = zdist, &
       Rmin = zRminu, Rmax = zRmaxu, Zmin = zZminu, Zmax = zZmaxu)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_dbdy: xplasma_mklim_special error:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eqm_dbdy

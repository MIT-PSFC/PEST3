subroutine eq_glimb1(rho,phi,b_min,chi_b_min,b_max,chi_b_max,ierr)

  use xplasma_obj_instance
  use eq_module

  !  find the minimum and maximum B on the given surface

  IMPLICIT NONE

  REAL*8 rho                        ! rho of surface (in)
  !  if a negative rho is supplied, reversed theta is implied!

  REAL*8 phi                        ! toroidal angle (in)
  REAL*8 b_min                      ! minimum mod(B) found (out)
  real*8 chi_b_min                  ! chi angle location of min(B)
  REAL*8 b_max                      ! maximum mod(B) found (out)
  real*8 chi_b_max                  ! chi angle location of max(B)
  integer ierr                      ! completion code, 0=OK

  real*8 rho_axis,rho_bdy,zrho

  logical :: axisymm,ccwflag
  real*8, parameter :: ZERO=0.0d0
  real*8, parameter :: ZTOL=1.0d-14
  !---------------------------

  call xplasma_global_info(s,ierr, axisymm=axisymm)
  if(ierr.ne.0) return

  if(.not.axisymm) then
     ierr=107
     call xplasma_error(s,ierr,lunerr)
     call eq_errmsg(' ?eq_glimb1:  non-axisymmetry not supported!')
     return
  endif

  ccwflag = rho.ge.ZERO
  zrho = abs(rho)

  call eq_rholim(rho_axis,rho_bdy)

  if(zrho.gt.rho_bdy+ZTOL) then
     call eq_errmsg(' ?eq_glimb1:  rho argument out of range:')
     write(lunerr,*) ' |rho|=',zrho,' rho_axis=',rho_axis,' rho_bdy=',rho_bdy
     ierr=9999
     return
  else
     zrho=min(rho_bdy,zrho)
  endif

  call xplasma_Bminmax(s,zrho,ierr, phi=phi, &
       ccw_theta=ccwflag, i2pi=2, bmin = B_min, bmax = B_max, &
       thbmin = chi_B_min, thbmax = chi_B_max)

end subroutine eq_glimb1

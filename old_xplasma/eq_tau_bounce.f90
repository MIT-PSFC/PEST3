subroutine eq_tau_bounce(zrho,inum,zx,ztol,zcut, &
     zbmin,zbmax,zvtaub,ierr)
 
  !  evaluate zero banana width approximation to v*tau_bounce
  !
  !  this is the integral over the limits of the orbit of
  !     dtheta*(dlp/dtheta)*(B/Bp)*(1/sqrt(1-B/Brefl))
  !
  !  where Brefl is a trapped orbit's reflection point (vpll=0,v=vperp)
  !  for passing orbits (Brefl > Bmax) the integrand is over the
  !  entire surface contour; for trapped orbits the integrand is
  !  twice the value from theta_lower to theta_upper where these
  !  limits need to be found by a root finder so that
  !     B(theta_lower)=B(theta_upper)=Brefl
  !
  !  the integrand formula follows from the assumptions that v
  !  and mu=vperp**2/B are constants of motion
 
  !  axisymmetry is assumed...
 
  use xplasma_obj_instance
  use eq_module

  implicit NONE
 
  real*8, intent(in) :: zrho   ! surface on which to calculate
  !  note -- must not be magnetic axis; non-singular surface required.
 
  integer, intent(in) :: inum  ! no. of evaluation points
  real*8, intent(in) :: zx(inum) ! evaluation points:
 
  !   zx(j) = (Brefl(j)-Bmin)/(Bmax-Bmin)
  ! where Bmin and Bmax, the min and mod(B) on the zrho surface, will
  ! be found; zx > 0 is required; zx(j+1) > zx(j) required for all j < inum.
 
  real*8, intent(in) :: ztol  ! relative error tolerance for integral
  ! (recommended value 1.0d-4)

  real*8, intent(in) :: zcut  ! minimum allowed value for (1-B/Brefl)
  ! (recommended value 1.0d-8)
 
  !---
 
  real*8, intent(out) :: zbmin,zbmax  ! Bmin and Bmax on zrho surface (T)
  real*8, intent(out) :: zvtaub(inum) ! integral values returned (m)
  !  (divide by ion velocity (m/sec) to get bounce time, seconds).
 
  integer, intent(out) :: ierr  ! completion code, 0=OK
 
  !-------------------------------------------

  call xplasma_vtaub(s,zrho,zx(1:inum),zvtaub(1:inum),ierr, &
       tol=ztol, bcut=zcut, bmin=zbmin, bmax=zbmax)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_tau_bounce: error detected:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_tau_bounce

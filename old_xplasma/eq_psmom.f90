!
! Return the Pfirsch-Schluter moments
!
! Below RHOMIN, the axis value is returned.
! The ZERROR argument returns the relative error of the moments
! sum = |sum(Fm) - <(n.grad(B))**2>/<B**2>| / max(<(n.grad(B))**2>/<B**2>)
!
! MOD DMC Feb 2008 -- used the new xplasma_psmom routine from xplasma2

subroutine eq_psmom(ivec, rho, rhomin, nmom, phi1, idmom, psmom, gamma, ngrdb2, b2, zerror, ier)

  use xplasma_obj_instance
  implicit NONE

  integer, intent(in) :: ivec       ! number of surfaces * sign(gamma)
  real*8,  intent(in) :: rho(abs(ivec))  ! radial coordinate
  real*8,  intent(in) :: rhomin     ! rho values below this are considered on axis
  integer, intent(in) :: nmom       ! number of moments to compute
  real*8,  intent(in) :: phi1       ! toroidal coordinate, same one used for all surfaces
  integer, intent(in) :: idmom      ! dimension of the psmom array

  real*8,  intent(out) :: psmom(idmom,abs(ivec))  ! PS moments
  real*8,  intent(out) :: gamma(abs(ivec))        ! wayne's gamma = 2*pi/integral_0_2pi[Bmod/B.grad(theta)]_dtheta
  real*8,  intent(out) :: ngrdb2(abs(ivec))       ! <(n.grad(B))**2>
  real*8,  intent(out) :: b2(abs(ivec))           ! <B**2>
  real*8,  intent(out) :: zerror(abs(ivec))       ! |sum(Fm) - <(n.grad(B))**2>/<B**2>| / max(<(n.grad(B))**2>/<B**2>)
  integer, intent(out) :: ier                ! nonzero on error

  !----------------------

  integer :: is    ! sign of argument chi, this was ichi_ccw but now ichi_ccw=+1 to reflect the storage in xplasma

  integer :: iveca ! abs(ivec)
  !----------------------
  iveca=abs(ivec)

  is=1
  if(ivec.lt.0) is=-1

  if(nmom.gt.idmom) then
     call eq_errmsg('eq_psmom: nmom exceeds array dimension (idmom).')
     ier=99
     return
  endif

  call xplasma_psmom(s, &
       rho,rhomin,psmom(1:nmom,1:iveca),gamma,ngrdb2,b2,zerror, &
       ier)

  if(ier.ne.0) then
     call eq_errmsg('eq_psmom: error in xplasma_psmom:')
  else
     gamma(1:iveca)  = is*gamma(1:iveca)
  endif

end subroutine eq_psmom

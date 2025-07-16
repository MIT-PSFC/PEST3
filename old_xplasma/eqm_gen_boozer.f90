subroutine eqm_gen_boozer(id_iboozer,id_nuboozer,ierr)

  ! generate profiles:
  !   IBoozer:  T*m*radians
  !      I(theta,rho) = Ib(rho)*thetab(theta,rho)
  !      where Ib(rho) = [Enclosed current, amps, vs. rho]*mu/2pi
  !      thetab(rho) = Boozer poloidal angle theta, 0:2pi
  !
  !   NuBoozer:  radians
  !      nu(theta,rho) -- periodic, 0 @ 1st & last theta point
  !        = zeta - phi, where zeta is the Boozer toroidal component.
  !
  !   These profiles satisfy
  !      mod(Bpol) = 
  !        (dLp/dtheta)**-1 * ( d(Iboozer)/dtheta + g(rho)*d(nu)/dtheta )
  !
  !      where g(rho) = R*B_phi, T*m.
  !
  !   The profiles are computed on the grids of the equilibrium 
  !   {R,Z}(theta,rho)
  !
  !   Note that I(theta,rho) has a branch cut at the end points of the
  !   equilibrium theta grid.  Since other theta grids could have different
  !   endpoints, user code may need to know this...

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  integer, intent(out) :: id_iboozer   ! I(theta,rho) xplasma ID
  integer, intent(out) :: id_nuboozer  ! nu(theta,rho) xplasma ID
  integer, intent(out) :: ierr         ! status code (0=OK)

  !---------------------------

  call xplasma_gen_boozer(s,id_iboozer,id_nuboozer,ierr)

  if(ierr.ne.0) then
     write(lunerr,*) &
          ' ?eqm_gen_boozer: xplasma_gen_boozer call failed:'
     call xplasma_error(s,ierr,lunerr)
  endif
 
end subroutine eqm_gen_boozer

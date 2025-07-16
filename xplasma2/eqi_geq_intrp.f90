subroutine eqi_geq_intrp(psi_geq, f_geq, ngeq, psi, ns, f, ier)
  !
  !  EFIT 1d profile interpolation routine
  !  (Hermite interpolation used, default boundary conditions)
  !
 
  use ezspline_obj
  use ezspline
 
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
 
  ! ------------------ arguments ---------------------------------
 
  integer, intent(in) :: ngeq        ! no. of pts in G-EQDSK 1d profiles
  real(r8), intent(in) :: psi_geq(ngeq)  ! G-EQDSK psi grid (strict ascending)
  real(r8), intent(in) :: f_geq(ngeq)    ! G-EQDSK f(psi) profile
 
  integer, intent(in) :: ns          ! no. of psi-surfaces wanted
  real(r8), intent(in) :: psi(ns)    ! psi values at psi-surfaces
 
  real(r8), intent(out) :: f(ns)     ! profile values at psi-surfaces
 
  integer, intent(out) :: ier        ! completion code:  0 = OK
 
  ! ------------------ local items -------------------------------
 
  type(ezspline1_r8) :: fspl
  integer iok
 
  ! ------------------ executable code ---------------------------
 
  ier=0
 
  call ezspline_init(fspl, ngeq, (/0,0/), iok)
  call ezspline_error(iok)
  if(iok.ne.0) ier=ier+1
 
  ! fspl%isHermite = 1  ! comment out to use spline
  fspl%x1 = psi_geq
 
  call ezspline_setup(fspl, f_geq, iok)
  call ezspline_error(iok)
  if(iok.ne.0) ier=ier+1
 
  call ezspline_interp(fspl, ns, psi, f, iok)
  call ezspline_error(iok)
  if(iok.ne.0) ier=ier+1
 
  call ezspline_free(fspl,iok)
  call ezspline_error(iok)
  if(iok.ne.0) ier=ier+1
 
  if(ier.ne.0) ier=120
 
  return
 
end subroutine eqi_geq_intrp

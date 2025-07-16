subroutine eq_rholim(zrho_axis,zrho_bdy)

  !  return the radial coordinate value at:
  !    the magnetic axis (zrho_axis)
  !    the plasma bdy    (zrho_bdy)

  !----------------------------------------------

  IMPLICIT NONE

  REAL*8, intent(out) :: zrho_axis,zrho_bdy

  !----------------------------------------------

  zrho_axis = 0
  zrho_bdy = 1

end subroutine eq_rholim

!================================================================

subroutine eq_rhoxlim(zrho_bdy,zrho_extend)

  use xplasma_obj_instance
  use eq_module

  !  return the external rho range covered (if any)
  !  if none, return (0,0)
  !----------------------------------------------

  IMPLICIT NONE

  REAL*8, intent(out) :: zrho_bdy,zrho_extend
  

  !----------------------------------------------
  integer :: id,ierr
  !----------------------------------------------

  zrho_bdy=0
  zrho_extend=0

  call xplasma_find_item(s,"__RHOX",id,ierr,nf_noerr=.TRUE.)
  if(id.eq.0) return

  zrho_bdy=1
  call xplasma_grid_info(s,id,ierr, xmax=zrho_extend)

end subroutine eq_rhoxlim

!================================================================

subroutine eq_rholima(zrho_axis,zrho_extend)

  use xplasma_obj_instance
  use eq_module

  !  return the rho range covered, including extrapolated region
  !----------------------------------------------

  IMPLICIT NONE

  real*8, intent(out) :: zrho_axis,zrho_extend

  !----------------------------------------------
  real*8 zrho_bdy
  !----------------------------------------------

  call eq_rhoxlim(zrho_bdy,zrho_extend)
  if((zrho_bdy.eq.czero).and.(zrho_extend.eq.czero)) then
     call eq_rholim(zrho_axis,zrho_extend)
  else
     call eq_rholim(zrho_axis,zrho_bdy)
  endif

end subroutine eq_rholima

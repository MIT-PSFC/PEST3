subroutine eq_glimrz(rho,r_min,r_max,z_min,z_max,ierr)

  use xplasma_obj_instance
  use eq_module

  !  find the minimum and maximum R and Z on the given surface

  IMPLICIT NONE
  REAL*8 rho                        ! rho of surface (in)
  REAL*8 R_min                      ! minimum R found (out)
  REAL*8 R_max                      ! maximum R found (out)
  REAL*8 Z_min                      ! minimum Z found (out)
  REAL*8 Z_max                      ! maximum Z found (out)
  integer ierr                      ! completion code, 0=OK (out)

  !---------------------------

  call xplasma_RZminmax(s,rho,ierr, rmin=r_min,rmax=r_max, &
       zmin=z_min,zmax=z_max)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?error in eq_glimrz.  Input argument rho=',rho,'.'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_glimrz

subroutine eq_glimrz1(rho,phi,r_min,r_max,z_min,z_max,ierr)

  use xplasma_obj_instance
  use eq_module

  !  find the minimum and maximum R and Z on the given surface
  !  at a particular toroidal angle (for axisymmetric cases this is the same
  !  as eq_glimrz).

  IMPLICIT NONE
  REAL*8 rho                        ! rho of surface (in)
  REAL*8 phi                        ! toroidal angle (in)
  REAL*8 R_min                      ! minimum R found (out)
  REAL*8 R_max                      ! maximum R found (out)
  REAL*8 Z_min                      ! minimum Z found (out)
  REAL*8 Z_max                      ! maximum Z found (out)
  integer ierr                      ! completion code, 0=OK (out)

  !---------------------------

  call xplasma_RZminmax(s,rho,ierr, phi=phi, rmin=r_min,rmax=r_max, &
       zmin=z_min,zmax=z_max)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?error in eq_glimrz1.'
     write(lunerr,*) '  Input arguments rho=',rho,' phi=',phi,'.'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_glimrz1

subroutine eqm_cbdy(ipts,rlim,zlim,ierr)

  !  xplasma f77 legacy interface
  !  set up an axisymmetric piecewise linear boundary / limiter contour
  !-----------------------

  use xplasma_obj_instance
  use eq_module
  
  implicit NONE

  !  input:

  integer ipts                      ! # of (R,Z) points defining contour
  real*8 rlim(ipts),zlim(ipts)      ! (R,Z) points

  !  output:

  integer ierr                      ! completion code, 0=OK

  !-----------------------
  integer :: iwarn   ! warning code is ignored...
  !-----------------------

  call xplasma_mklim_contour(s,rlim,zlim,iwarn,ierr)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_cbdy: xplasma_mklim_contour error:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eqm_cbdy

subroutine eq_limcon(imax,inum,rlim,zlim,dtol,ierr)

  !  for any limiter configuration, convert to a closed piecewise
  !  linear representation.  This will contain up to (imax) pieces
  !  but will be reduced by eliminating points that are colinear to
  !  tolerance (dtol); to prevent a reduction & get the maximum number
  !  of points, set dtol.le.0.0d0

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  !  input:

  integer imax                      ! max no. of contour points output

  !  output:

  integer inum                      ! actual no. of contour points output
  real*8 rlim(imax),zlim(imax)      ! the contour points

  !  input:

  real*8 dtol                       ! tolerance for estimating colinearity

  !  output:

  integer ierr                      ! completion code, 0=OK
  !---------------------------------------

  call xplasma_limcon(s,rlim,zlim,inum,ierr, &
       tol=dtol,maptype = 2)

end subroutine eq_limcon

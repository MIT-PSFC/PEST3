subroutine eqx_dist(ivec,xx,yy,zz,zdist,ilpt,rlim,zlim,philim,ierr)

  !  **vectorized**
  !  compute distance from limiter: .gt.0 means outside bdy
  !  .lt.0 means inside, # gives actual distance in (R,Z) units
  !  relative to nearest limiting object.

  !  if switch is set, also return coordinates of nearest contact
  !  point on the limiting object.

  implicit NONE

  ! input:
  integer ivec                      ! vector dimension
  real*8 xx(ivec),yy(ivec),zz(ivec) ! locations

  integer ilpt                      ! =1 to calculate nearest limiter pt

  ! output:
  real*8 zdist(ivec)                ! signed distance from bdy/limiter
  integer ierr                      ! completion code, 0=OK

  ! set if ilpt=1:
  real*8 rlim(ivec),zlim(ivec),philim(ivec) ! nearest limiter contact pt

  !--------------------------
  real*8 zR(ivec),zphi(ivec)
  !--------------------------

  call eq_rcyl(ivec,xx,yy,zR,zphi)
  call eq_rzdist(ivec,zR,zZ,zphi,zdist,ilpt,rlim,zlim,philim,ierr)

end subroutine eqx_dist

!----------------------------------------------------------------------
subroutine eq_rzdist(ivec,zR,zZ,zphi,zdist,ilpt,rlim,zlim,philim,ierr)

  !  **vectorized**
  !  compute distance from limiter: .gt.0 means outside limiter
  !  .lt.0 means inside, # gives actual distance in (R,Z) units
  !  relative to nearest limiting object.

  !  if switch is set, also return coordinates of nearest contact
  !  point on the limiting object.

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  ! input:
  integer ivec                      ! vector dimension
  real*8 zR(ivec),zZ(ivec),zphi(ivec) ! location

  integer ilpt                      ! =1 to calculate nearest limiter pt

  ! output:
  real*8 zdist(ivec)                ! signed distance from bdy/limiter
  integer ierr                      ! completion code, 0=OK

  real*8 rlim(ivec),zlim(ivec),philim(ivec) ! nearest limiter contact pt

  !---------------------------------------------------
  logical :: axisymm
  !---------------------------------------------------

  call xplasma_global_info(s,ierr, axisymm=axisymm)
  if(.not.axisymm) then
     write(lunerr,*) ' ?eq_rzdist: non-axisymmetric geometry not supported!'
     ierr=1
     return
  endif

  call xplasma_lim_distance(s,zR,zZ,zdist,ierr, rlim=rlim,zlim=zlim)

  philim = 0
  if(ierr.ne.0) then
     rlim = 0
     zlim = 0
     write(lunerr,*) ' ?error detected in eq_rzdist:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_rzdist

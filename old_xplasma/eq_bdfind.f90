subroutine eqx_bdfind(ivec,xx,yy,zz,zrho_in,zchi_out,zphi_out,zdist,ierr)

  !  accurately find the point on a boundary surface which is
  !  closest to a target point (x,y,z).

  !  input:

  IMPLICIT NONE

  integer ivec
  REAL*8 xx(abs(ivec)),yy(abs(ivec)),zz(abs(ivec)) ! cartesian target points
  REAL*8 zrho_in                    ! surface to be searched

  !  output:

  REAL*8 zchi_out(abs(ivec))             ! chi of nearest approach
  REAL*8 zphi_out(abs(ivec))             ! phi of nearest approach

  REAL*8 zdist(abs(ivec))                ! distance outside surface;

  integer ierr                      ! completion code, 0=OK

  !--------------------------
  real*8 zR(abs(ivec)),zphi(abs(ivec))
  !--------------------------

  call eq_rcyl(abs(ivec),xx,yy,zR,zphi)
  call eq_bdfind(ivec,zR,zZ,zphi,zrho_in,zchi_out,zphi_out,zdist,ierr)

end subroutine eqx_bdfind

!--------------------------------------------------------------------
subroutine eq_bdfind(ivec,zR,zZ,zphi,zrho_in, zchi_out,zphi_out,zdist,ierr)

  !  **vectorized** dmc Apr 2000

  !  accurately find the point on a boundary surface which is
  !  closest to a target point (zR,zZ,zphi)

  !  to begin the search, evaluate a coarse grid of points and look
  !  for the distance minima; use a root finder to zero in exactly
  !  (to machine precision).

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  !  input:

  integer ivec                      ! vector length
  REAL*8 zR(abs(ivec)),zZ(abs(ivec)),zphi(abs(ivec)) ! target point
  REAL*8 zrho_in                    ! the (1) surface to be searched

  !  output:

  REAL*8 zchi_out(abs(ivec))             ! chi of nearest approach
  REAL*8 zphi_out(abs(ivec))             ! phi of nearest approach

  REAL*8 zdist(abs(ivec))                ! distance outside surface;

  integer ierr                      ! completion code, 0=OK

  !  if zdist.lt.0.0 then point is *inside* the surface;
  !  if zdist.eq.0.0 then the point is *on* the surface.
  !------------------------------------------------------------
  logical :: iccw
  !------------------------------------------------------------

  iccw = ivec.gt.0

  ierr = 9999

  call xplasma_bdfind(s,abs(ivec),zR,zZ,ierr, &
       phi_in=zphi,rho_in=zrho_in,ccw_theta=iccw,i2pi=2, &
       theta_out=zchi_out,phi_out=zphi_out,dist=zdist)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_rcyl: error in xplasma_bdfind:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_bdfind

subroutine eq_xget(ivec,zrho,zchi,zphi,xx,yy,zz,nregion,ierr)

  !  vectorized
  !  mapping:  (rho,chi,phi) -> (R,Z,phi) -> (x,y,z)

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  integer ivec                      ! vector dimensioning (- for reversed chi)
  real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! mag. coordinates
  real*8 xx(abs(ivec)),yy(abs(ivec)),zz(abs(ivec)) ! cartesian coordinates
  
  integer nregion(abs(ivec))             ! output region code
  integer ierr                      ! output error code

  !  on output,
  !      nregion = 1 means:  rho_axis.le.rho.le.rho_bdy  (plasma)
  !      nregion = 2 means:  rho_bdy.lt.rho.le.rho_vac    (scrapeoff.eq.1)
  !      nregion = 3 means:  rho.gt.rho_vac -- out of bounds

  !      ierr=0:  normal exit.  ierr=1 only on a serious error.
  !--------------------------------------------------------------

  REAL*8 zR(abs(ivec))

  !--------------------------------------------------------------

  call eq_rzget(ivec,zrho,zchi,zphi,zR,zZ,nregion,ierr)
  if(ierr.eq.0) call eq_xcyl(abs(ivec),zR,zphi,xx,yy)

end subroutine eq_xget

!--------------------------------------------------------------
subroutine eq_rzget(ivec,zrho,zchi,zphi,zR,zZ,nregion,ierr)

  !  mapping:  (rho,chi,phi) -> (R,Z,[phi])

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  integer ivec                      ! vector dimensioning (- for reversed chi)
  real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! mag. coordinates
  real*8 zR(abs(ivec)),zZ(abs(ivec))          ! Cylindric coordinates (w/phi)

  integer nregion(abs(ivec))             ! output region code
  integer ierr                      ! output error code

  !  on output,
  !      nregion = 1 means:  rho_axis.le.rho.le.rho_bdy  (plasma)
  !      nregion = 2 means:  rho_bdy.lt.rho.le.rho_vac    (scrapeoff.eq.1)
  !      nregion = 3 means:  rho.gt.rho_vac (out of bounds)

  !      ierr=0:  normal exit.  ierr=1 means a serious error.

  !--------------------------------------------------------------
  logical :: iccw
  !--------------------------------------------------------------

  iccw = ivec.gt.0

  call xplasma_ctrans(s,xp_use_coord_vecsize,ierr, &
       rho_in=zrho, theta_in=zchi, phi_in=zphi, ccw_theta=iccw, i2pi=2, &
       r_out=zR, z_out=zZ, nregion=nregion)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqi_rzget: error in xplasma_ctrans:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_rzget

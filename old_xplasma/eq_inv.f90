subroutine eqx_inv(ivec,xx,yy,zz,tol,zrho,zchi,zphi,nregion,ierr)

  !  **vectorized**
  !  mapping:  (x,y,z) -> (R,Z,phi) -> (rho,chi,phi)

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  !  sign of ivec indicates chi reversal.  set ivec = -[# of pts] for the
  !  poloidal angle orientation to be reversed, i.e. tracing clockwise along
  !  a flux surface drawn to the right of the machine centerline...

  integer ivec                      ! vector dimension:  # of pts
  real*8 xx(abs(ivec)),yy(abs(ivec)),zz(abs(ivec)) ! input (x,y,z) target pts
  REAL*8 tol                        ! input relative error tolerance

  REAL*8 zrho(abs(ivec)),zchi(abs(ivec)) ! output (rho,chi) flux coordinates
  real*8 zphi(abs(ivec))                 ! output (phi) toroidal angle
  integer nregion(abs(ivec))             ! output region code

  integer ierr                      ! output error code

  !  on output,
  !     nregion(i) = 1 means:  rho_axis.le.rho.le.rho_bdy (core plasma)
  !     nregion(i) = 2 means:  rho_bdy.lt.rho.le.rho_vac (SOL)
  !     nregion(i) = 3 means:  rho.gt.rho_vac or rho.lt.rho_axis, out of bounds

  !      ierr.ne.0 means a serious error occurred.

  !------------------------
  real*8 zR(abs(ivec))
  !------------------------

  call eq_rcyl(abs(ivec),xx,yy,zR,zphi)
  call eq_inv(ivec,zR,zZ,zphi,tol,zrho,zchi,nregion,ierr)

end subroutine eqx_inv

!-----------------------------------------------------------------
subroutine eq_inv(ivec,zR,zZ,zphi,tol,zrho,zchi,nregion,ierr)

  !  **vectorized**
  !  mapping:  (R,Z,phi) -> (rho,chi,phi)

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  !  sign of ivec indicates chi reversal.  set ivec = -[# of pts] for the
  !  poloidal angle orientation to be reversed, i.e. tracing clockwise along
  !  a flux surface drawn to the right of the machine centerline...

  integer ivec
  REAL*8 zR(abs(ivec)),zZ(abs(ivec)),zphi(abs(ivec)) ! input (R,Z,phi) target
  REAL*8 tol                             ! input relative error tolerance

  !  if tol.gt.1.0d-3 use maptype=2 approximation; otherwise use maptype=1
  !  which is slower but more accurate.

  REAL*8 zrho(abs(ivec)),zchi(abs(ivec)) ! output (rho,chi) flux coords.

  !  on input, if the fast inverse map is disabled, (R,Z,phi) contains the
  !  initial guess.  a good initial guess improves performance, but is
  !  not required for accurate results.

  integer nregion(abs(ivec))             ! output region code
  integer ierr                           ! output error code; 0=OK

  !  on output,
  !      nregion = 1 means:  rho_axis.le.rho.le.rho_bdy
  !      nregion = 2 means:  rho_bdy.lt.rho.le.rho_vac
  !      nregion = 3 or 4 means:  rho.gt.rho_vac or rho.lt.rho_axis

  !   ierr .ne. zero means a serious error in the code
  !--------------------------------------------------------------
  logical :: iccw
  integer :: imaptype
  !--------------------------------------------------------------

  iccw = ivec.gt.0

  if(tol.ge.1.0d-3) then
     imaptype=xp_imap_polar
  else
     imaptype=xp_imap_newton
  endif

  call xplasma_ctrans(s,xp_use_coord_vecsize,ierr, &
       R_in=zR, z_in=zZ, phi_in=zphi, ccw_theta = iccw, i2pi=2, &
       maptype = imaptype, &
       tol=tol, rho_out=zrho, theta_out=zchi, nregion=nregion)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_inv: error in xplasma_ctrans:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_inv

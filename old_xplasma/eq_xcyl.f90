subroutine eq_xcyl(ivec,zR,zphi,xx,yy)

  !  vectorized
  !  cylindric to cartesian map.

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  integer ivec       ! vector dimension; ivec>0: normal Phi; ivec<0 -- reversed
  REAL*8 zR(abs(ivec)),zphi(abs(ivec))        ! input: (R,[Z],phi)
  REAL*8 xx(abs(ivec)),yy(abs(ivec))          ! output:  (x,y,[z])

  !-----------
  integer :: ier
  logical :: iccw
  !-----------

  iccw=ivec.gt.0
  call xplasma_ctrans(s,xp_use_coord_vecsize,ier, &
       r_in=zR,phi_in=zphi, ccw_phi=iccw, &
       x_out=xx,y_out=yy)

  if(ier.ne.0) then
     write(lunerr,*) ' ?eq_xcyl: error in xplasma_ctrans:'
     call xplasma_error(s,ier,lunerr)
  endif

end subroutine eq_xcyl

!----------------------------------------------------------------

subroutine eq_xcylcs(ivec,zR,zphi,xx,yy,zcosphi,zsinphi)

  !  vectorized.
  !  convert from cylindric coordinates (R,Z,phi) to
  !  cartesian coordinates (x,y,z).  also return cos(phi),sin(phi)

  !  z->Z
  !  x = R*cos(phi)
  !  y = R*sin(phi)

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  integer ivec       ! vector dimension; ivec>0: normal Phi; ivec<0 -- reversed
  REAL*8 zR(abs(ivec)),zphi(abs(ivec))         ! input: (R,phi)
  REAL*8 xx(abs(ivec)),yy(abs(ivec))           ! output:  (x,y)
  real*8 zcosphi(abs(ivec)),zsinphi(abs(ivec)) ! output:  cos(phi),sin(phi)

  !-------------------------------
  integer :: ier
  logical :: iccw
  !-----------

  iccw=ivec.gt.0
  call xplasma_ctrans(s,xp_use_coord_vecsize,ier, &
       r_in=zR,phi_in=zphi, ccw_phi=iccw, &
       x_out=xx,y_out=yy, cosphi=zcosphi,sinphi=zsinphi)

  if(ier.ne.0) then
     write(lunerr,*) ' ?eq_xcylcs: error in xplasma_ctrans:'
     call xplasma_error(s,ier,lunerr)
  endif

end subroutine eq_xcylcs

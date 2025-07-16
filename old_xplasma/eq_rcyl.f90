subroutine eq_rcyl(ivec,xx,yy,zR,zphi)

  !  *vectorized*
  !  cartesian to cylindric conversion; see eq_rcylcs, below
  !    ivec > 0: normal Phi; < 0: Phi orientation reversed

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  integer ivec                 ! vector dimension
  REAL*8 xx(abs(ivec)),yy(abs(ivec))    ! input:  (x,y,[z])
  REAL*8 zR(abs(ivec)),zphi(abs(ivec))  ! output: (R,phi,[z])  ! [z] unaffected

  !------------
  integer :: ier
  logical :: iccw
  !------------

  iccw = ivec.gt.0
  call xplasma_ctrans(s,xp_use_coord_vecsize,ier, &
       x_in=xx, y_in=yy, ccw_phi = iccw, i2pi=2, &
       r_out=zR, phi_out=zphi)

  if(ier.ne.0) then
     write(lunerr,*) ' ?eq_rcyl: error in xplasma_ctrans:'
     call xplasma_error(s,ier,lunerr)
  endif

end subroutine eq_rcyl

!--------------------------------------------------------------------
subroutine eq_rcylcs(ivec,xx,yy,zR,zphi,zcosphi,zsinphi)

  !  *vectorized*
  !  convert from cartesian coordinates (x,y,z) to
  !  cylindric coordinates (R,Z,phi)

  !  also return cos(phi), sin(phi)

  !  sqrt(x**2+y**2)->R
  !  y=0, x>0 half line corresponds to phi=0; phi=atan2(y,x)

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  integer ivec                      ! vector dimension
  REAL*8 xx(abs(ivec)),yy(abs(ivec))    ! input:  (x,y,[z])
  REAL*8 zR(abs(ivec)),zphi(abs(ivec))  ! output: (R,phi,[z])  ! [z] unaffected

  real*8 zcosphi(abs(ivec)),zsinphi(abs(ivec)) ! cos(phi),zin(phi)

  !------------
  integer :: ier
  logical :: iccw
  !-------------------------------

  iccw = ivec.gt.0
  call xplasma_ctrans(s,xp_use_coord_vecsize,ier, &
       x_in=xx, y_in=yy, ccw_phi = iccw, i2pi=2, &
       r_out=zR, phi_out=zphi, cosphi=zcosphi, sinphi=zsinphi)

  if(ier.ne.0) then
     write(lunerr,*) ' ?eq_rcylcs: error in xplasma_ctrans:'
     call xplasma_error(s,ier,lunerr)
  endif

end subroutine eq_rcylcs

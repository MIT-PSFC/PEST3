subroutine eqx_dfast(ivec,xx,yy,zz,zdist,ierr)

  !  using fast map, find approximate distance to plasma boundary

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  integer ivec
  REAL*8 xx(ivec),yy(ivec),zz(ivec) ! cartesian target points

  !  output:

  REAL*8 zdist(ivec)                ! distance outside surface;

  integer ierr                      ! completion code, 0=OK

  !  if zdist.lt.0.0 then point is *inside* the surface;
  !  if zdist.eq.0.0 then the point is *on* the surface.
  !--------------------------
  integer :: maptype
  real*8 :: zR(ivec)
  !--------------------------
  
  call xplasma_max_maptype(s,maptype,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' xplasma_max_maptype call error in eqx_dfast:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  call xplasma_ctrans(s,xp_use_coord_vecsize,ierr, &
       x_in=xx, y_in=yy, r_out=zR)
  if(ierr.ne.0) then
     write(lunerr,*) ' xplasma_ctrans call error in eqx_dfast:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  call xplasma_bdfind(s,ivec,zR,zZ,ierr, maptype=maptype, dist=zdist)
  if(ierr.ne.0) then
     write(lunerr,*) ' xplasma_bdfind call error in eqx_dfast:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

end subroutine eqx_dfast

!--------------------------------------------------------------------
subroutine eq_dfast(ivec,zR,zZ,zphi,zdist,ierr)

  !  **vectorized** dmc Apr 2000
  !  using fast map, find approximate distance to plasma boundary

  !  to begin the search, evaluate a coarse grid of points and look
  !  for the distance minima; use a root finder to zero in exactly
  !  (to machine precision).

  !  input:

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  integer ivec                        ! vector length
  REAL*8 zR(ivec),zZ(ivec),zphi(ivec) ! target point

  !  output:

  REAL*8 zdist(ivec)                  ! distance outside surface;

  integer ierr                        ! completion code, 0=OK

  !  if zdist.lt.0.0 then point is *inside* the surface;
  !  if zdist.eq.0.0 then the point is *on* the surface.
  !------------------------------------------------------------
  integer :: maptype
  !--------------------------
  
  call xplasma_max_maptype(s,maptype,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' xplasma_max_maptype call error in eq_dfast:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  call xplasma_bdfind(s,ivec,zR,zZ,ierr, maptype=maptype, dist=zdist)
  if(ierr.ne.0) then
     write(lunerr,*) ' xplasma_bdfind call error in eq_dfast:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

end subroutine eq_dfast

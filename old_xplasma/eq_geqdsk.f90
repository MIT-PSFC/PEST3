subroutine eq_geqdsk(lun_geqdsk,geqdsk_lbl, &
     Rmin,Rmax,Zmin,Zmax,zcur, &
     id_p,id_q,nh,nv,nb,ierr)
  
  !  write GEQDSK file from xplasma module

  !  **correction** dmc 9 Nov 2001:
  !  follow Lang Lao sign convention for G-EQDSK files:
  !  -- psi always increasing
  !  -- direction of current specified in pcur scalar *only*

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

  !  input:

  integer lun_geqdsk                ! logical unit where to write file

  character*48 geqdsk_lbl           ! 48 character label for GEQDSK file
  real*8 Rmin,Rmax                  ! (Rmin,Rmax) of Psi(R,Z)
  real*8 Zmin,Zmax                  ! (Zmin,Zmax) of Psi(R,Z)
  real*8 zcur                       ! plasma current (amps)

  !  [Rmin,Rmax]x[Zmin,Zmax] must contain the core plasma but not exceed the
  !  available (R,Z) grids.

  integer id_p                      ! xplasma id:  Pressure profile
  integer id_q                      ! xplasma id:  q profile

  integer nh                        ! #of GEQDSK horizontal pts
  integer nv                        ! #of GEQDSK vertical pts

  !  note nh also controls the number of pts in the 1d profiles

  integer nb                        ! #of pts in bdy contour

  !  output:

  integer ierr                      ! completion code (0=OK, file written).

  !------------------------------------------------------
  ! (old f77 xplasma eq_geqdsk code moved into core modules; most arguments
  ! made optional).

  call xplasma_wr_geqdsk(s,ierr, &
       lun_geqdsk=lun_geqdsk, label=geqdsk_lbl, &
       Rmin=Rmin, Rmax=Rmax, Zmin=Zmin, Zmax=Zmax, &
       cur=zcur, id_pprof=id_p, id_qprof=id_q, &
       nh=nh, nv=nv, nbdy=nb)

  if(ierr.ne.0) then
     write(lunerr,*) ' %eq_geqdsk: error in xplasma_wr_geqdsk:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_geqdsk

subroutine eq_geqdsk_ipq(lun_geqdsk,geqdsk_lbl,zcur,id_p,id_q,ierr)

  ! another version, with a simplified calling argument list:  I, P, Q...
  !  input:

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

  integer lun_geqdsk                ! logical unit where to write file

  character*(*) geqdsk_lbl           ! 48 character label for GEQDSK file
  real*8 zcur                       ! plasma current (amps)

  integer id_p                      ! xplasma id:  Pressure profile
  integer id_q                      ! xplasma id:  q profile

  !  output:

  integer ierr                      ! completion code (0=OK, file written).

  !------------------------------------------------------
  ! (old f77 xplasma eq_geqdsk code moved into core modules; most arguments
  ! made optional).

  call xplasma_wr_geqdsk(s,ierr, &
       lun_geqdsk=lun_geqdsk, label=geqdsk_lbl, &
       cur=zcur, id_pprof=id_p, id_qprof=id_q)

  if(ierr.ne.0) then
     write(lunerr,*) ' %eq_geqdsk_ipq: error in xplasma_wr_geqdsk:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_geqdsk_ipq


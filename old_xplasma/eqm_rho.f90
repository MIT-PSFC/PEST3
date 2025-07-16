subroutine eqm_rho(zrho,ibdy,ztol,id,ierr)

  use xplasma_obj_instance
  use eq_module

  !  set up rho grid -- must be strict ascending
  !
  !  must be the first discretization of the "rho" coordinate
  !  must span the range [0,1]

  IMPLICIT NONE

  !  input:

  integer, intent(in) :: ibdy       ! #of grid bdy pts -- out to boundary
  REAL*8, intent(in) :: zrho(ibdy)  ! the grid (zone bdys, ibdy-1 zones).
  REAL*8, intent(in) ::  ztol       ! even spacing tolerance

  !  because of spline requirements:
  !    ibdy.ge.4 required

  !  the core plasma is covered by zrho(1:ibdy)

  !  output:

  integer, intent(out) :: id        ! axis id code (returned)
  integer, intent(out) :: ierr      ! =0: OK

  !-----------------

  integer i,inrho,iertmp
  REAL*8 zzrho(ibdy)

  !-------------------------------------------------------------
  !  error checks:

  ierr=0
  id=0

  if(.not.eq_module_init) then
     write(lunerr,*) ' ?eqm_rho: initialization: call eqm_select first!'
     ierr=1
     return
  else 
     call xplasma_coord_info(s,xplasma_rho_coord,ierr, ngrids=inrho)
     if(ierr.ne.0) then
        call xplasma_error(s,ierr,lunerr)
     else if(inrho.gt.0) then
        write(lunerr,*) ' ?eqm_rho: rho already defined!'
        write(lunerr,*) '  re-initialization: call eqm_select first.'
        ierr=1
        return
     endif
  endif
  if(ierr.ne.0) return

  zzrho=zrho
  if(abs(zzrho(1)).le.max(ceps4,ztol)) zzrho(1)=0
  if(abs(zzrho(ibdy)-1).le.max(ceps4,ztol)) zzrho(ibdy)=1
  call eqi_evenx(zzrho,ibdy,ztol)

  call xoi_author_set(ierr)

  call xplasma_create_grid(s,'__RHO',xplasma_rho_coord,zzrho, &
       id,ierr)

  call xoi_author_clear(iertmp)

  return
end subroutine eqm_rho

subroutine eqm_rzmag2(rzmapper,id_R,id_Z,ierr)

  !  legacy xplasma f77 interface...
  !  create 2d functions R(rho,chi), Z(rho,chi)
  !  using passed external routine "rzmapper"
  !  rzmapper is a dummy argument for a subroutine which MUST have the
  !  following interface:

  !      subroutine rzmapper(rho,chi,R,Z,ierr)

  !      real*8 rho,chi      ! INPUT (rho,chi) - radial, poloid coordinates
  !      real*8 R,Z          ! OUTPUT R(rho,chi), Z(rho,chi)
  !      integer IERR        ! OUTPUT completion code, 0=OK

  !--------------------------------------------------------------

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  ! input:
  external rzmapper                 ! passed subroutine

  ! output:
  integer id_R                      ! id for R spline
  integer id_Z                      ! id for Z spline
  integer IERR                      ! completion code, 0-OK
  
  !--------------------------
  !  local

  integer :: id_rho,id_chi,inrho,inchi,iertmp,icount
  real*8, parameter :: zdiff=1.0d0
  integer, parameter :: ninfin=1000000000
  logical :: iforce=.FALSE.
  !--------------------------

  id_R=0
  id_Z=0

  call eqm_rzmag_check(ninfin,ninfin,1,id_rho,id_chi,inrho,inchi,ierr)
  if(ierr.ne.0) return

  call xoi_author_set(iertmp)

  call xplasma_rzmagf(s,rzmapper,id_rho,id_chi, &
       id_R,id_Z,ierr, zdiff=zdiff)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_rzmag2: error in xplasma_rzmagf call.'
     call xplasma_error(s,ierr,lunerr)
  else
     call xplasma_eqcheck(s,iforce,icount,ierr)
     if(icount.gt.0) then
        write(lunerr,*) ' %eqm_rzmag2: equilibrium specification completed.'
     endif
  endif

  call xoi_author_clear(iertmp)

end subroutine eqm_rzmag2

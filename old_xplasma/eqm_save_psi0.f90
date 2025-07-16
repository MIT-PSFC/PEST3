subroutine eqm_save_psi0(psi0,ierr)

  ! store Psi0 -- EFIT Psi offset, mag. axis to machine axis

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  real*8, intent(in) :: psi0        ! Psi0 value, Wb/rad
  integer, intent(out) :: ierr      ! error status code; 0=OK

  !----------------------
  integer :: iertmp
  !----------------------

  call xoi_author_set(ierr)
  if(ierr.ne.0) return

  call xplasma_save_psi0(s,psi0,ierr)

  call xoi_author_clear(iertmp)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_save_psi0: error detected in xplasma: '
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eqm_save_psi0

subroutine eq_get_psi0(psi0,ifind,ierr)

  ! retrieve Psi0 -- EFIT Psi offset, mag. axis to machine axis
  !  ...if available

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  real*8, intent(out) :: psi0       ! Psi0 value, Wb/rad
  integer, intent(out) :: ifind     ! =1 if value is found, =0 otherwise
  integer, intent(out) :: ierr      ! error status code; 0=OK

  ! if no Psi0 data exists, ifind=ierr=0 are returned.

  !----------------------

  call xplasma_get_psi0(s,psi0,ifind,ierr)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_get_psi0: error detected in xplasma: '
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_get_psi0

subroutine eq_rhopsin(npsi,rhovals,ierr)

  !  generate an evenly spaced psi grid and return the corresponding
  !  rho values.  The first point in the psi grid corresponds to the
  !  magnetic axis; the last point, the plasma boundary.

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer npsi                      ! no. of psi pts; rho(psi) pts returned
  real*8 rhovals(npsi)              ! the rho values returned

  integer ierr                      ! completion code (0=OK)

  !--------------------------------------

  call xplasma_rhopsi_gen(s,npsi,ierr, rhovals=rhovals)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_rhopsin: xplasma_rhopsi_gen failed: '
     call xplasma_error(s,ierr,lunerr)
  endif
end subroutine eq_rhopsin

!----------------------------------------------------------------------
subroutine eq_rhopsi(npsi,psivals,rhovals,ztol,ierr)

  ! given a vector of psi values find the corresponding rho values
  ! -- to the indicated accuracy tolerance in rho

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  !  input:

  integer npsi                      ! no. of psivals (vector size)
  real*8 psivals(npsi)              ! the psi's for which rho's are wanted.

  !  output:
  real*8 rhovals(npsi)              ! the rho values returned

  !  input:
  real*8 ztol                       ! accuracy tolerance

  !  output:
  integer ierr                      ! completion code, 0=OK

  !------------------------------------

  call xplasma_rhopsi_find(s,psivals,rhovals,ierr, tol=ztol)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_rhopsin: xplasma_rhopsi failed: '
     call xplasma_error(s,ierr,lunerr)
  endif
end subroutine eq_rhopsi

subroutine eqm_gen_q(qname,id,ierr)
 
  use xplasma_obj_instance
  use eq_module

  implicit NONE
 
  ! use XPLASMA to generate a q profile
  ! q = phi'/(2pi*psi') = (dV/drho)*<1/R**2>*g / (dpsi/drho) / (2pi)**2
 
  !-------------------------------------------------------
  character*(*), intent(in) :: qname   ! name to assign to this q profile
  integer, intent(out) :: id           ! xplasma id of profile created
  integer, intent(out) :: ierr         ! completion code, 0=OK
 
  !  if an error occors, id=0 on output.
  !-------------------------------------------------------
  integer :: iorder
  !-------------------------------------------------------
 
  iorder=1  ! Generate a Hermite profile

  call xplasma_gen_q(s,qname,iorder,id,ierr)

  if(ierr.ne.0) then
     write(lunerr,*) &
          ' ?eqm_gen_q: xplasma_gen_q call failed:'
     call xplasma_error(s,ierr,lunerr)
  endif
 
end subroutine eqm_gen_q

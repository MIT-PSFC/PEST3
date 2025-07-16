subroutine eq_restore(filename,ier)
 
  !  f77 interface -- restore xplasma from file

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  character*(*), intent(in) :: filename
  integer, intent(out) :: ier

  !-----------------------
  call eqm_clear

  call xplasma_read(s,filename,ier)
  if(ier.ne.0) then

     call xplasma_error(s,ier,lunerr)

  endif

end subroutine eq_restore

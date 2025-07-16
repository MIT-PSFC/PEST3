subroutine eq_save(filename,ier)
 
  !  f77 interface -- save xplasma to file

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  character*(*), intent(in) :: filename
  integer, intent(out) :: ier

  call xplasma_write(s,filename,ier)
  if(ier.ne.0) then

     call xplasma_error(s,ier,lunerr)

  endif

end subroutine eq_save

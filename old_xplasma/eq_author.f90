subroutine eq_author_set(zname,ierr)

  ! set the author name

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  character*(*), intent(in) :: zname
  integer, intent(out) :: ierr

  call xplasma_author_set(s,zname,ierr)

end subroutine eq_author_set


subroutine eq_author_clear(zname,ierr)

  ! clear the author name

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  character*(*), intent(in) :: zname
  integer, intent(out) :: ierr

  call xplasma_author_clear(s,zname,ierr)

end subroutine eq_author_clear

subroutine eq_inquire_author(id,zauth,ier)

  ! find the author of an existing item

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  integer, intent(in) :: id
  character*(*), intent(out) :: zauth
  integer, intent(out) :: ier

  !---------------------------------------
  integer :: iertmp
  !---------------------------------------

  call xplasma_get_item_info(s,id,ier, author=zauth)
  if(ier.ne.0) then
     write(lunerr,*) ' ?eq_inquire_author: error code detected.'
     call xplasma_error(s,ier,lunerr)
     return
  endif

end subroutine eq_inquire_author

subroutine eqm_set_units(id,zuns,ier)

  ! add units label to xplasma item

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  integer, intent(in) :: id
  character*(*), intent(in) :: zuns
  integer, intent(out) :: ier

  !---------------------------------------
  integer :: iertmp
  character*32 zauth
  !---------------------------------------

  call xplasma_get_item_info(s,id,ier, author=zauth)
  if(ier.ne.0) then
     write(lunerr,*) ' ?eqm_set_units: error code detected.'
     call xplasma_error(s,ier,lunerr)
     return
  endif

  call xplasma_author_set(s,zauth,iertmp)

  call xplasma_label_item(s,id,ier, units=zuns)

  call xplasma_author_clear(s,zauth,iertmp)

  if(ier.ne.0) then
     write(lunerr,*) ' ?eqm_set_units: error code detected.'
     call xplasma_error(s,ier,lunerr)
  endif

end subroutine eqm_set_units

subroutine eqm_flabel(id,label,units,ier)
  !
  ! set descriptive and units label for item (usually a profile function)

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  integer, intent(in) :: id   ! id of item to be labeled
  character*(*), intent(in) :: label   ! desired label
  character*(*), intent(in) :: units   ! desired units label

  integer, intent(out) :: ier ! completion code, 0=OK

  !-----------------------
  integer :: iertmp
  !-----------------------

  call xoi_author_set(iertmp)

  call xplasma_label_item(s,id,ier, &
       label=label,units=units)

  if(ier.ne.0) call xplasma_error(s,ier,lunerr)

  call xoi_author_clear(iertmp)

end subroutine eqm_flabel

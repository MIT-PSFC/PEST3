subroutine eqm_mark_pmhd(id,ier)

  ! mark specified profile as being the MHD pressure...

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  integer, intent(in) :: id
  integer, intent(out) :: ier

  call xplasma_set_pressure_id(s,id,ier)
  if(ier.ne.0) then

     call xplasma_error(s,ier,lunerr)

  endif

end subroutine eqm_mark_pmhd

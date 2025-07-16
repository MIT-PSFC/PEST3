subroutine eqi_fromg_err(zmsg)

  !  append error or warning message & increment counter

  use xplasma_definitions
  use eqi_rzbox_module
  use eqi_geq_mod  ! contains jwarn counter
  
  implicit NONE

  character*(*), intent(in) :: zmsg

  !---------------------------

  jwarn = jwarn + 1
  if(jwarn.eq.0) then
     call xplasma_errmsg_append(sp,' --> messages from eqi_fromgeqdsk:')
  endif

  call xplasma_errmsg_append(sp,zmsg)

end subroutine eqi_fromg_err

subroutine eqi_fromg_warn(iwarn)

  use eqi_geq_mod  ! contains jwarn counter

  iwarn = jwarn

end subroutine eqi_fromg_warn

subroutine eqm_clear

  ! delete contents of xplasma object "s"; free up memory

  use xplasma_obj_instance
  use eq_module

  integer :: ier

  if(.not.eq_module_init) then

     call xoi_init(ier)
     if(ier.ne.0) write(lunerr,*) &
          ' ?eqm_clear: xplasma_init returned ier = ',ier
     eq_module_init = .TRUE.

  else

     call xplasma_free(s,ier)

  endif

  eq_premsg_text = ' '

end subroutine eqm_clear

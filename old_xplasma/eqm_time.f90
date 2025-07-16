subroutine eqm_time(ztime)

  !  set the time (seconds)

  use xplasma_obj_instance
  use eq_module

  real*8, intent(in) :: ztime

  !-------------------
  integer :: ier
  !-------------------

  call xplasma_time_set(s,ztime,ier)
  if(ier.ne.0) then
     write(lunerr,*) ' %eqm_time: failed to set xplasma time, ier=',ier
     call xplasma_error(s,ier,lunerr)
  endif

end subroutine eqm_time

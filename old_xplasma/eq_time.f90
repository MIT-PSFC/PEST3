subroutine eq_time(ztime)

  !  fetch the time (seconds)

  use xplasma_obj_instance
  use eq_module

  real*8, intent(out) :: ztime

  !-------------------
  integer :: ier
  !-------------------

  call xplasma_time_get(s,ztime,ier)
  if(ier.ne.0) then
     write(lunerr,*) ' %eq_time: failed to get xplasma time, ier=',ier
     call xplasma_error(s,ier,lunerr)
  endif

end subroutine eq_time

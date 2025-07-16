subroutine eqi_geq_getlun(ilun)

  !  get GEQ module I/O unit number

  use eqi_geq_mod
  implicit NONE

  integer :: ilun

  !--------------------------

  call geq_getlun(ilun)

end subroutine eqi_geq_getlun

module eqi_rzbox_module

  use xplasma_obj
  implicit NONE

  private
  public :: sp,nsrch0,nsrch1

  type (xplasma), pointer :: sp

  integer, parameter :: nsrch0=60       ! #of poloidal zones for search
  integer, parameter :: nsrch1=nsrch0+1 ! #of zone bdys (1st & last coincide)

end module eqi_rzbox_module

module eqi_geq_axis_data

  ! module to contain spline variables used by eqi_geq_axis.f90
  ! and additional routines passed to a rootfinder by eqi_geq_axis
  ! DMC Jan 2011 -- effort to improve reliability of root finder

  implicit NONE
  SAVE

  real*8, dimension(:,:), allocatable :: rpkg,zpkg
  real*8, dimension(:,:,:), allocatable :: psirz_spl

  ! for single point evaluation:

  integer :: ii(1),jj(1)
  real*8 :: rparam(1),zparam(1),hr(1),hri(1),hz(1),hzi(1)

  real*8 :: eqi_eta,eqi_eps

end module eqi_geq_axis_data


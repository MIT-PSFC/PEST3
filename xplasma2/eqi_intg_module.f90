module eqi_intg_module

  !  a public module with some XPLASMA pointers
  !  enable communication to integrand routines via integrators, to 
  !  debug plotting routines, etc., where passing a pointer through 
  !  the entire call chain is not practical or desirable...

  use xplasma_obj
  implicit NONE

  private
  public :: sp,sp_debug

  type (xplasma), pointer :: sp,sp_debug

  integer, parameter, public :: eqi_intg_nfi=10

  integer, parameter, public :: eqi_intg_vol=51         ! Volume
  integer, parameter, public :: eqi_intg_area=52        ! Area
  integer, parameter, public :: eqi_intg_dvdrho=53      ! dV/drho
  integer, parameter, public :: eqi_intg_r2i=54         ! <1/R^2>
  integer, parameter, public :: eqi_intg_itor=55        ! enclosed current
  integer, parameter, public :: eqi_intg_lpol=56
  integer, parameter, public :: eqi_intg_surf=57
  integer, parameter, public :: eqi_intg_dvol=58
  integer, parameter, public :: eqi_intg_darea=59
  integer, parameter, public :: eqi_intg_g2r2i=60
  integer, parameter, public :: eqi_intg_ri=61
  integer, parameter, public :: eqi_intg_r3i=62
  integer, parameter, public :: eqi_intg_r2=63
  integer, parameter, public :: eqi_intg_r1=64
  integer, parameter, public :: eqi_intg_bi=65
  integer, parameter, public :: eqi_intg_b2i=66
  integer, parameter, public :: eqi_intg_b2=67
  integer, parameter, public :: eqi_intg_b1=68
  integer, parameter, public :: eqi_intg_bz2=69
  integer, parameter, public :: eqi_intg_g2b2i=70
  integer, parameter, public :: eqi_intg_g1=71
  integer, parameter, public :: eqi_intg_g2=72
  integer, parameter, public :: eqi_intg_g2r3i=73
  integer, parameter, public :: eqi_intg_g2r2=74
  integer, parameter, public :: eqi_intg_giri=75
  integer, parameter, public :: eqi_intg_fhnclass=76
  integer, parameter, public :: eqi_intg_fhnc_tsc=77

  integer, parameter, public :: eqi_intg_booz_qovg=78

  integer, parameter, public :: eqi_intg_fhmmx=79

  integer, parameter, public :: eqi_intg_br2=80

  integer, parameter, public :: eqi_intg_r=1
  integer, parameter, public :: eqi_intg_drdth=2
  integer, parameter, public :: eqi_intg_drdrho=3

  integer, parameter, public :: eqi_intg_z=4
  integer, parameter, public :: eqi_intg_dzdth=5
  integer, parameter, public :: eqi_intg_dzdrho=6

  integer, parameter, public :: eqi_intg_modb=7
  integer, parameter, public :: eqi_intg_br=8
  integer, parameter, public :: eqi_intg_bz=9

  integer, parameter, public :: eqi_intg_2pirdetj=10

  integer, parameter, public :: eqi_intg_bign=1000000
end module eqi_intg_module

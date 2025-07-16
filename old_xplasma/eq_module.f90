module eq_module

  ! private module used in implementation of original f77 style XPLASMA
  ! interface... i.e. private to the xplasma library, not meant for use
  ! by subroutines outside xplasma.

  ! NOTE: this module is meant to be used in the xplasma library along with
  ! the public module "xplasma_obj_instance" which makes a shared xplasma
  ! object available.

  implicit NONE
  save

  logical :: eq_module_init = .FALSE.

  integer :: lunerr = 6

  character*132 :: eq_premsg_text = ' '

  !--------------------------
  !  constants...

  real*8, parameter :: c2pi = 6.2831853071795862D+00
  real*8, parameter :: cpi = 3.1415926535897931D+00
  real*8, parameter :: ceps4 = 1.0d-4
  real*8, parameter :: czero = 0.0d0
  real*8, parameter :: chalf = 0.5d0
  real*8, parameter :: cone = 1.0d0
  real*8, parameter :: ctwo = 2.0d0

end module eq_module

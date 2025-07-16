subroutine eq_volume(ivec,zrho,iwant,zval,ierr)

  !  Volume(rho) interpolation routine -- with error checking
  !  units: m**3

  !  see also subroutine eq_area (very similar)

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  !  input:

  integer ivec                      ! vector dimension
  REAL*8 zrho(ivec)                 ! argument
  integer iwant                     ! 0:  value, 1: df/drho, 2: d2f/drho2

  !  output:

  REAL*8 zval(ivec)                 ! result of interpolation
  integer ierr                      ! =0: OK, =1: zrho out of range

  !--------------------------

  call xplasma_volume(s,zrho,zval,ierr, ideriv=iwant)
  if(ierr.ne.0) then

     write(lunerr,*) ' ?eq_volume: error detected:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_volume

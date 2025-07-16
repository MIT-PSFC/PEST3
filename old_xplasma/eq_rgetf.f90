subroutine eq_rgetf(ivec,zrho,ifcn,iwant,zval,ierr)

  !  f(rho) interpolation routine -- with error checking

  use xplasma_obj_instance
  use eq_module
  IMPLICIT NONE

  !  input:

  integer ivec                      ! vector dimension
  REAL*8 zrho(ivec)                 ! argument

  integer ifcn                      ! function id number
  integer iwant                     ! 0:  value, 1: df/drho, 2: d2f/drho2

  !  output:
  REAL*8 zval(ivec)                 ! result of interpolation
  integer ierr                      ! =0: OK, =1: zrho out of range
  !  ... or some other error...
  !----------------------------

  integer :: ioutside

  !----------------------------

  call xplasma_eval_prof(s,ifcn,zrho,zval,ierr, &
       ideriv1=iwant, n_out_of_bounds=ioutside)

  if(ierr.ne.0) then
     write(lunerr,*) ' %eq_rgetf: error flag set: ',ierr
     call xplasma_error(s,ierr,lunerr)
     if(ioutside.gt.0) write(lunerr,*) &
          '  number of target points out of bounds: ',ioutside
  endif
  
  return
end subroutine eq_rgetf

subroutine eq_rgetff(ivec,zrho,nlist,ifcns,iwant,ivecd,zvals,ierr)

  !  f(rho) interpolation routine -- with error checking

  use xplasma_obj_instance
  use eq_module
  IMPLICIT NONE

  !  input:

  integer ivec                      ! input vector dimension
  REAL*8 zrho(ivec)                 ! argument

  integer nlist                     ! # of fcns to evaluate
  integer ifcns(nlist)              ! function numbers
  integer iwant                     ! 0:  value, 1: df/drho, 2: d2f/drho2

  !  output:
  integer ivecd                     ! output vector dimension, .ge. ivec
  REAL*8 zvals(ivecd,nlist)         ! result of interpolations

  integer ierr                      ! =0: OK, =1: zrho out of range
  !  ... or some other error...
  !----------------------------

  integer :: ioutside

  !----------------------------

  if(ivecd.lt.ivec) then
     ierr=99
     write(lunerr,*) ' ?eq_rgetff: dimension arguments: ivecd < ivec: ', &
          ivecd,ivec
     return
  endif

  call xplasma_eval_prof(s,ifcns,zrho,zvals(1:ivec,1:nlist),ierr, &
       ideriv1=iwant, n_out_of_bounds=ioutside)

  if(ierr.ne.0) then
     write(lunerr,*) ' %eq_rgetf: error flag set: ',ierr
     call xplasma_error(s,ierr,lunerr)
     if(ioutside.gt.0) write(lunerr,*) &
          '  number of target points out of bounds: ',ioutside
  endif
  
  return
end subroutine eq_rgetff

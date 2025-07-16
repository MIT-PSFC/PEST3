subroutine eq_xgetf(ivec,zx,ifcn,iwant,zval,ierr)

  !  f(x) interpolation routine -- with error checking

  use xplasma_obj_instance
  use eq_module
  IMPLICIT NONE

  !  input:

  integer ivec                    ! vector dimension
  REAL*8 zx(abs(ivec))                 ! argument
  integer ifcn                    ! function id number

  !  input:
  integer iwant                   ! 0:  value, 1: df/dx, 2: d2f/dx2

  !  output:
  REAL*8 zval(abs(ivec))               ! result of interpolation
  integer ierr                    ! =0: OK, =1: zx out of range

  !---------------------------
  !  the caller provides the name of the function desired, but it also
  !  provides a number which is writable.  If the caller saves this number
  !  and repeats the call, it saves having to lookup the function name,
  !  a performance enhancement.

  !---------------------------
  integer :: irank,idg,iertmp,ioutside,inum
  logical :: ccwflag
  character*32 fname
  !---------------------------

  call xplasma_prof_info(s,ifcn,ierr, rank=irank, gridId1=idg)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_xgetf: xplasma_prof_info error detected:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  if(irank.ne.1) then
     call xplasma_prof_info(s,ifcn,iertmp, name=fname)
     write(lunerr,*) ' ?eq_xgetf: not a rank1 profile: "',trim(fname),'".'
     ierr=1
     return
  endif

  ccwflag = (ivec.gt.0) 
  inum=abs(ivec)

  call xplasma_eval_prof(s,ifcn,zx(1:inum),zval(1:inum),ierr, &
       ideriv1=iwant,ccwflag1=ccwflag,n_out_of_bounds=ioutside)

  if(ioutside.gt.0) then
     ierr=1
     write(lunerr,*) ' %eq_xgetf: error flag set, ',ioutside, &
          ' points out of range.'
  else if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_xgetf: xplasma_eval_prof error detected:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_xgetf

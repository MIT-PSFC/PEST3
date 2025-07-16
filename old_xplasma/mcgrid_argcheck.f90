subroutine mcgrid_argcheck(subname,ist_th,ist_ph,iextend,ierr)
 
  ! partial argument checker for mcgrid_get/putobj subroutines

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  character*(*), intent(in) :: subname

  integer, intent(in) :: ist_th,ist_ph,iextend
  integer, intent(out) :: ierr

  !---------------------------------------------------------

  ierr = 0

  if((iextend.ne.0).and.(iextend.ne.1)) then
     write(lunerr,*) ' ?'//trim(subname)//':  iextend = ',iextend
     write(lunerr,*) &
          '  legal values are 0 (inside plasma) and 1 (inside and outside).'
     ierr=1
  endif
 
  if((ist_th.ne.0).and.(ist_th.ne.-1)) then
     write(lunerr,*) ' ?'//trim(subname)//':  theta start parameter ist_th = ',ist_th
     write(lunerr,*) &
          '  legal values are 0 (start at 0) and -1 (start at -pi)'
     ierr=1
  endif
 
  if((ist_ph.ne.0).and.(ist_ph.ne.-1)) then
     write(lunerr,*) ' ?'//trim(subname)//':  theta start parameter ist_ph = ',ist_ph
     write(lunerr,*) &
          '  legal values are 0 (start at 0) and -1 (start at -pi)'
     ierr=1
  endif

end subroutine mcgrid_argcheck

subroutine eqm_wr_geqdsk(filename,ier)
 
!
!  write a G-EQDSK file; data read earlier by eqm_fromgeqdsk
!
  use eqi_geq_mod
  implicit NONE
 
  character*(*), intent(in) :: filename
  integer, intent(out) :: ier             ! completion code, 0=OK
 
!----------------
 
  ier = 0
  call geq_writ(geq, filename, ier)
 
  return
 
end subroutine eqm_wr_geqdsk

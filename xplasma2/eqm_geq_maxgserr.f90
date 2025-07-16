subroutine eqm_geq_maxGSerr(zmaxerr)
  use eqi_geq_mod
  implicit NONE
 
  real*8, intent(out) :: zmaxerr
 
  zmaxerr = gs_errmax

end subroutine eqm_geq_maxGSerr

subroutine eqi_jacheck_maxvar_Set(value)

  !  set the maximum allowed variation of det[J] (metric Jacobian) on a
  !  flux surface

  use xplasma_obj_instance
  use eq_module

  real*8, intent(in) :: value

  !-------------------
  integer :: ier
  !-------------------

  call xplasma_ajac_maxvar_set(s,value,ier)
  if(ier.ne.0) then
     write(lunerr,*) ' %eqi_jacheck_maxvar_set: failed to set xplasma max Jacobian variation, ier=',ier
     call xplasma_error(s,ier,lunerr)
  endif

end subroutine eqi_jacheck_maxvar_Set

subroutine eqi_jacheck_maxvar_Get(value)

  !  fetch the maximum allowed variation of det[J] (metric Jacobian) on a
  !  flux surface

  use xplasma_obj_instance
  use eq_module

  real*8, intent(out) :: value

  !-------------------
  integer :: ier
  !-------------------

  call xplasma_global_info(s,ier, ajac_maxVar=value)
  if(ier.ne.0) then
     write(lunerr,*) ' %eqi_jacheck_maxvar_get: failed to retrieve xplasma max Jacobian variation, ier=',ier
     call xplasma_error(s,ier,lunerr)
  endif

end subroutine eqi_jacheck_maxvar_Get

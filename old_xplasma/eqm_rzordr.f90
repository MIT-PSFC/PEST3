subroutine eqm_rzordr_set(iorder)

  !  set fit order for R(rho,chi), Z(rho,chi), and B(rho,chi) component

  !      iorder = 1 ==> Hermite fits
  !      iorder = 2 ==> (default) Spline fits
 
  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: iorder

  !-----------------------
  integer :: iord_use,ier
  !-----------------------

  iord_use = max(1,min(2,iorder))

  call xplasma_rzOrder_set(s,iord_use,ier)

  if(ier.ne.0) then
     write(lunerr,*) ' ?eqm_rzordr_set: unexpected error.'
     call xplasma_error(s,ier,lunerr)
  endif

end subroutine eqm_rzordr_set

subroutine eq_rzordr_get(iorder)

  !  get fit order for R(rho,chi), Z(rho,chi), and B(rho,chi) component

  !      iorder = 1 ==> Hermite fits
  !      iorder = 2 ==> (default) Spline fits
 
  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(out) :: iorder

  !-----------------------
  integer :: ier
  !-----------------------

  call xplasma_global_info(s,ier, rzOrder=iorder)

  if(ier.ne.0) then
     write(lunerr,*) ' ?eqm_rzordr_set: unexpected error.'
     call xplasma_error(s,ier,lunerr)
  endif

end subroutine eq_rzordr_get

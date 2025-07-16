subroutine eq_find_intg(inumth,idi,inrhoi,inthi,ierr)

  !  xplasma2 f77 interface-- internal routine
  !  find integrator

  use xplasma_obj_instance
  use eq_module
  IMPLICIT NONE

  integer, intent(in) :: inumth  ! theta dimension of output array
  !  if =1, assume flux surface integrals; if <>1 assume zonal integrals

  integer, intent(out) :: idi    ! integrator data set ID or 0 if error

  integer, intent(out) :: inrhoi ! #rho surfaces in integrator (0 if error)
  integer, intent(out) :: inthi  ! #theta surfaces in integrator (0 if error)

  integer, intent(out) :: ierr   ! completion code 0=OK

  ! #rho zones = inrhoi-1; #theta zones = inthi-1; for flux surface averages
  !                                                inthi=2, #theta zones = 1.

  !----------------------------------
  !----------------------------------

  idi = 0
  inrhoi = 0
  inthi = 0

  if(inumth.eq.1) then
     call xplasma_find_item(s,"__F77_INTEGRATOR_1D",idi,ierr)
     if((ierr.ne.0).or.(idi.eq.0)) then
        write(lunerr,*) ' ?eq_flxint: 1D integrator dataset not found.'
        write(lunerr,*) '  prior call to "eq_flxint_init" required.'
        ierr=1

     else
        inthi=2
        call xplasma_integ_info(s,idi,ierr, n_rho_int=inrhoi)

     endif

  else

     call xplasma_find_item(s,"__F77_INTEGRATOR_2D",idi,ierr)
     if((ierr.ne.0).or.(idi.eq.0)) then
        write(lunerr,*) ' ?eq_flxint_arr2d: 2D integrator dataset not found.'
        write(lunerr,*) '  prior call to "eq_flxint_chinit" required.'
        return
        
     else
        call xplasma_integ_info(s,idi,ierr, n_rho_int=inrhoi, n_th_int=inthi)

     endif
  endif

end subroutine eq_find_intg

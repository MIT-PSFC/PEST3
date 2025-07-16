subroutine fluxav_nzones_get(nthzones,nzones)

  !  get number of poloidal and radial zones in current f77 xplasma integrator

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(out) :: nthzones,nzones

  !----------------------------------
  integer :: idi,ierr,inum,inumth
  !----------------------------------

  nthzones=0
  nzones=0

  call xplasma_find_item(s,"__F77_INTEGRATOR_2D",idi,ierr,nf_noerr=.TRUE.)

  if(idi.eq.0) then
     call xplasma_find_item(s,"__F77_INTEGRATOR_1D",idi,ierr,nf_noerr=.TRUE.)

     if(idi.eq.0) then
        write(lunerr,*) ' ?fluxav_nzones_get: no F77 integrator found; call eq_flxint first.'

     else
        call xplasma_integ_info(s,idi,ierr,n_rho_int=inum)
        if(ierr.eq.0) then
           nthzones=1
           nzones=inum-1
        endif
     endif

  else

     call xplasma_integ_info(s,idi,ierr, n_rho_int=inum, n_th_int=inumth)
     if(ierr.eq.0) then
        nthzones=inumth-1
        nzones=inum-1
     endif

  endif

  if(ierr.ne.0) then

     write(lunerr,*) ' ?fluxav_nzones_get: unexpected error ierr=',ierr
     call xplasma_error(s,ierr,lunerr)

  endif

end subroutine fluxav_nzones_get

subroutine fluxav_nsurfs_get(nthzones,nsurfs)

  !  get number of poloidal zones and radial surfaces in current f77 
  !  xplasma integrator

  implicit NONE

  integer, intent(out) :: nthzones,nsurfs

  call fluxav_nzones_get(nthzones,nsurfs)
  if(nthzones.gt.0) nsurfs = nsurfs + 1

end subroutine fluxav_nsurfs_get

subroutine xmoments_kmom_set(ikmom)
 
  !  set no. of moments for Fourier representation of 2d equilibrium
  !  xplasma will accept 2 <= ikmom <= 64 (dmc: as of May 2006)
  !  There is no error check, this routine assumes "s" has been initialized.

  use xplasma_obj_instance
  use eq_module
  implicit NONE
 
  integer, intent(in) :: ikmom
 
  !------------------
  integer :: ier
  !------------------

  call xoi_init(ier)
  call xoi_kmom_set(ikmom)
  call xplasma_kmom_set(s,ikmom,ier)  ! error not checked...

end subroutine xmoments_kmom_set
 
!--------------------
 
subroutine xmoments_kmom_get(ikmom)
 
  !  get no. of moments for Fourier representation of 2d equilibrium
  !  (if there is an error, ikmom=0 is returned).

  use xplasma_obj_instance
  use eq_module
  implicit NONE
 
  integer, intent(out) :: ikmom
 
  !------------------
  integer :: ier
  !------------------

  call xoi_init(ier)
  call xoi_kmom_get(ikmom)
  call xplasma_kmom_set(s,ikmom,ier)  ! error not checked...
 
end subroutine xmoments_kmom_get

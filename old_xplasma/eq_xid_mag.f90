subroutine eq_xid_mag(ifnum,id_rank,id_rho,id_chi,id_phi,iorder,ierr)

  !  return information on a function vs. mag. coordinates

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: ifnum      ! function id

  integer, intent(out) :: id_rank   ! dimensionality of object
  integer, intent(out) :: id_rho    ! rho axis id (0 if none)
  integer, intent(out) :: id_chi    ! chi axis id (0 if none)
  integer, intent(out) :: id_phi    ! phi axis id (0 if none)
  integer, intent(out) :: iorder    ! object's fit order

  integer, intent(out) :: ierr      ! completion code, (0=OK)

  !---------------------------------------------
  integer :: iertmp
  !---------------------------------------------

  id_rho=0
  id_chi=0
  id_phi=0

  call xplasma_prof_info(s,ifnum,ierr, &
       rank=id_rank, splineType=iorder)
  if(ierr.ne.0) then
     id_rank=0
     iorder=-99
     write(lunerr,*) ' ?xplasma error detected in eq_xid_mag.f90:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  call xplasma_prof_gridInfo(s,ifnum,xplasma_rho_coord,id_rho,iertmp)
  call xplasma_prof_gridInfo(s,ifnum,xplasma_theta_coord,id_chi,iertmp)
  call xplasma_prof_gridInfo(s,ifnum,xplasma_phi_coord,id_phi,iertmp)

end subroutine eq_xid_mag

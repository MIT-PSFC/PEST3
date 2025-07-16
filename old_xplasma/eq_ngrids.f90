subroutine eq_ngrids(inumrho,inumchi,inumphi,inumR,inumZ)

  !  return the sizes of "common" grids-- this means, the first defined
  !  grid for each coordinate.

  !  zeroes are returned for undefined grids.
  !  for axisymmetric cases inumphi=1 is returned.

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer inumrho                   ! radial flux coordinate grid
  integer inumchi                   ! poloidal angle grid
  integer inumphi                   ! toroidal angle grid

  integer inumR                     ! R grid
  integer inumZ                     ! Z grid

  !----------------------------------
  logical :: axisymm
  integer :: ierr
  !----------------------------------

  inumrho = 0
  inumchi = 0
  inumphi = 0

  inumR = 0
  inumZ = 0

  call xplasma_global_info(s,ierr, axisymm=axisymm); if(ierr.ne.0) return

  if(axisymm) then
     inumphi=1
  else
     call get1('__phi_coord',inumphi)
  endif

  call get1('__rho_coord',inumrho)
  call get1('__theta_coord',inumchi)

  call get1('__R_coord',inumR)
  call get1('__Z_coord',inumZ)

  contains

    subroutine get1(cname,ingr1)

      !  find size of 1st grid discretizing named coordinate...

      character*(*), intent(in) :: cname   ! name of coordinate
      integer, intent(out) :: ingr1        ! size of 1st grid

      integer :: id_coord,id_grid,inum

      ingr1=0

      call xplasma_find_item(s,cname,id_coord,ierr)
      if(ierr.ne.0) then
         write(lunerr,*) ' ?eq_ngrids: unexpected error:'
         call xplasma_error(s,ierr,lunerr)
         return
      endif

      call xplasma_coord_info(s,id_coord,ierr, ngrids=inum)
      if(ierr.ne.0) then
         write(lunerr,*) ' ?eq_ngrids: unexpected error:'
         call xplasma_error(s,ierr,lunerr)
         return
      endif

      if(inum.eq.0) return  ! no grid; ingr1=0 stands.

      call xplasma_coord_gridId(s,id_coord,1,id_grid,ierr)
      if(ierr.ne.0) then
         write(lunerr,*) ' ?eq_ngrids: unexpected error:'
         call xplasma_error(s,ierr,lunerr)
         return
      endif

      call xplasma_grid_size(s,id_grid,ingr1,ierr)
      if(ierr.ne.0) then
         ingr1=0
         write(lunerr,*) ' ?eq_ngrids: unexpected error:'
         call xplasma_error(s,ierr,lunerr)
         return
      endif

    end subroutine get1

end subroutine eq_ngrids

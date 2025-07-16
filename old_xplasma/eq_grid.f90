subroutine eq_grid(id_axis,zdata,nmax,ngot,ierr)

  !  retrieve the axis grid identified by id_axis

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  !  input:

  integer, intent(in) :: id_axis    ! xplasma id of desired grid
  integer, intent(in) :: nmax       ! size of buffer to hold grid data

  !  output:

  real*8 zdata(nmax)                ! buffer for grid data
  integer, intent(out) :: ngot      ! number of points in grid, retrieved
  integer, intent(out) :: ierr      ! completion code, 0=OK

  !  ierr gets set to #pts in grid, if #pts in grid exceeds nmax
  !---------------------------------------------

  call eq_ngrid(id_axis,ngot)
  if(ngot.eq.0) then
     ierr=77
     write(lunerr,*) ' ?eq_grid:  invalid axis id:  ',id_axis
     zdata=czero

  else if(ngot.gt.nmax) then
     write(lunerr,*) ' ?eq_grid:  data buffer too small, have=',nmax, &
          ' need=',ngot
     ierr=ngot
     ngot=0
     zdata=czero

  else
     call xplasma_grid(s,id_axis,zdata(1:ngot),ierr)
     if(ierr.ne.0) then
        write(lunerr,*) ' ?eq_grid: unexpected error in xplasma_grid:'
        call xplasma_error(s,ierr,lunerr)
        ngot=0
        zdata=czero
     else
        if(ngot.lt.nmax) zdata(ngot+1:nmax)=czero
     endif

  endif

end subroutine eq_grid

subroutine eq_grid_zc(id_axis,zdata,nmax,ngot,ierr)

  !  retrieve the axis grid "zone centers" as identified by id_axis
  !  the grid itself is assumed to define boundaries of zones.  There
  !  are 1 fewer zone centers; zone center j = 0.5*(bdy(j)+bdy(j+1))

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id_axis    ! xplasma id of desired grid
  integer, intent(in) :: nmax       ! size of buffer to hold grid data

  !  output:

  real*8 zdata(nmax)                ! buffer for grid zc data
  integer, intent(out) :: ngot      ! number of points in grid, retrieved
  integer, intent(out) :: ierr      ! completion code, 0=OK

  !  ierr gets set to #pts in grid, if #pts in grid exceeds nmax
  !---------------------------------------------

  integer :: i
  real*8, dimension(:), allocatable :: zdata2

  !---------------------------------------------

  call eq_ngrid(id_axis,ngot)
  if(ngot.eq.0) then
     ierr=77
     write(lunerr,*) ' ?eq_grid_zc:  invalid axis id:  ',id_axis
     zdata=czero

  else if(ngot-1.gt.nmax) then
     write(lunerr,*) ' ?eq_grid_zc:  data buffer too small, have=',nmax, &
          ' need=',ngot-1
     ierr=ngot-1
     ngot=0
     zdata=czero

  else
     ierr=0
     allocate(zdata2(ngot))

     call xplasma_grid(s,id_axis,zdata2,ierr)
     if(ierr.ne.0) then
        write(lunerr,*) ' ?eq_grid_zc: unexpected error in xplasma_grid:'
        call xplasma_error(s,ierr,lunerr)
        ngot=0
        zdata=czero
     else
        do i=1,ngot-1
           zdata(i)=chalf*(zdata2(i)+zdata2(i+1))
        enddo
        ngot=ngot-1
     endif
     deallocate(zdata2)
  endif

end subroutine eq_grid_zc

subroutine eq_ngrid(id_axis,isize)

  !  return the grid size, or, 0 if id_axis is invalid.

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id_axis     ! xplasma axis id code
  integer, intent(out) :: isize      ! no. of pts in axis grid

  !-----------------------------------
  integer :: itype,ierr
  !-----------------------------------

  isize=0
  call xplasma_get_item_info(s,id_axis,ierr, itype=itype)
  if(ierr.ne.0) return
  if(itype.ne.xplasma_gridType) return

  call xplasma_grid_size(s,id_axis,isize,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_ngrid: unexpected error in xplasma_grid_size:'
     call xplasma_error(s,ierr,lunerr)
     isize=0
  endif

end subroutine eq_ngrid

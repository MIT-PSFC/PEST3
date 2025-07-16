subroutine eq_ganum(zname,ianum)

  use xplasma_obj_instance
  use eq_module

  !  given the name of a grid (axis) return its number.  if there is no such
  !  grid, silently return ianum=0

  IMPLICIT NONE

  character*(*), intent(in) :: zname   ! (input) name of grid
  integer, intent(out) :: ianum        ! (output) grid id or ZERO

  !----------------------------
  integer :: itype,id,ierr
  !----------------------------

  call xplasma_gridId(s,zname,ianum)

end subroutine eq_ganum
 
!==================================================

subroutine eq_get_aname(ianum,zname)

  use xplasma_obj_instance
  use eq_module

  !  given id, get name of axis grid.  if an invalid id is given,
  !  a blank name is returned.

  IMPLICIT NONE

  integer, intent(in) :: ianum          ! grid id number
  character*(*), intent(out) :: zname   ! name of grid

  !-------------------------------
  integer :: itype,ierr
  !-------------------------------

  zname=' '

  call xplasma_get_item_info(s,ianum,ierr, nf_noerr=.TRUE., &
       name=zname, itype=itype)

  if((ierr.ne.0).or.(itype.ne.xplasma_gridType)) then
     zname=' '
  endif

end subroutine eq_get_aname

subroutine eq_gfnum(zname,ifnum)

  use xplasma_obj_instance
  use eq_module

  !  given the name of a profile return its number.  if there is no such
  !  profile, silently return ifnum=0

  IMPLICIT NONE

  character*(*), intent(in) :: zname   ! (input) name of profile
  integer, intent(out) :: ifnum        ! (output) profile id or ZERO

  !----------------------------
  integer :: itype,id,ierr
  !----------------------------

  call xplasma_profId(s,zname,ifnum)

end subroutine eq_gfnum
 
!==================================================

subroutine eq_get_fname(ifnum,zname)

  use xplasma_obj_instance
  use eq_module

  !  given id, get name of axis profile.  if an invalid id is given,
  !  a blank name is returned.

  IMPLICIT NONE

  integer, intent(in) :: ifnum          ! profile id number
  character*(*), intent(out) :: zname   ! name of profile

  !-------------------------------
  integer :: itype,ierr
  !-------------------------------

  call xplasma_get_item_info(s,ifnum,ierr, nf_noerr=.TRUE., &
       name=zname, itype=itype)

  if((ierr.ne.0).or.(itype.ne.xplasma_profType)) then
     zname=' '
  endif

end subroutine eq_get_fname

subroutine eq_gfitem(zname,ifnum)

  ! return ID of item of indeterminate type
  ! return ID=0 if no such item exists

  use xplasma_obj_instance
  use eq_module


  character*(*), intent(in) :: zname   ! (input) name of profile
  integer, intent(out) :: ifnum        ! (output) item id or ZERO

  !----------------------------
  integer :: ierr
  !----------------------------

  idtmp = 0
  call xplasma_find_item(s,zname,ifnum,ierr,nf_noerr=.TRUE.)

end subroutine eq_gfitem

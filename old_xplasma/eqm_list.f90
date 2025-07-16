subroutine eqm_list_create(zlist_name,zlist_label,inum,elem_names,id,ierr)
  !
  !  create new list object in XPLASMA

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  character*(*), intent(in) :: zlist_name   ! name of list
  character*(*), intent(in) :: zlist_label  ! label for list
  integer, intent(in) :: inum               ! size of list (1 or more)
  character*(*), dimension(inum), intent(in) :: elem_names ! list element names

  integer, intent(out) :: id                ! id of newly created list
  integer, intent(out) :: ierr              ! completion code: 0=OK

  ! id=0 on exit if an error occurs.
  !
  ! list name must be unique.  Duplicate name => error.
  ! inum > 0 required.
  !
  ! internally, all names are set to uppercase; comparisons e.g. for
  ! duplication are all done in uppercase.
  !-------------------------------------------------
  integer :: iertmp
  !-------------------------------------------------

  ierr=0
  id=0

  call xplasma_find_item(s,zlist_name,id,ierr)
  if(ierr.eq.0) then
     !  error -- item already exists
     write(lunerr,*) ' ?eqm_list_create: cannot make list named "', &
          trim(zlist_name),'":'
     write(lunerr,*) '  name is already in use.'
     ierr=1
     return
  endif

  call xoi_author_set(iertmp)

  call eqm_list_replace(zlist_name,zlist_label,inum,elem_names,id,ierr)

  call xoi_author_clear(iertmp)

end subroutine eqm_list_create

subroutine eqm_list_replace(zlist_name,zlist_label,inum,elem_names,id,ierr)
  !
  !  replace list object in XPLASMA or if does not exist
  !  create new list object in XPLASMA

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  character*(*), intent(in) :: zlist_name   ! name of list
  character*(*), intent(in) :: zlist_label  ! label for list
  integer, intent(in) :: inum               ! size of list (1 or more)
  character*(*), dimension(inum), intent(in) :: elem_names ! list element names

  integer, intent(out) :: id                ! id of newly created list
  integer, intent(out) :: ierr              ! completion code: 0=OK

  ! id=0 on exit if an error occurs.
  !
  ! list name must be unique.  Duplicate name => error.
  ! inum > 0 required.
  !
  ! internally, all names are set to uppercase; comparisons e.g. for
  ! duplication are all done in uppercase.
  !-------------------------------------------------
  integer :: iertmp
  !-------------------------------------------------

  ierr=0
  id=0

  call xoi_author_set(iertmp)

  call xplasma_create_list(s,zlist_name,elem_names,id,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_list_create: error creating list named "', &
          trim(zlist_name),'":'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  call xplasma_label_item(s,id,ierr, label=zlist_label)

  call xoi_author_clear(iertmp)

end subroutine eqm_list_replace

subroutine eqm_list_merge(id1,id2,idelete,ierr)
  ! merge lists
  ! list (id2) is appended to list (id1); when done (id2) is deleted
  !   if (idelete) is set.
  ! if an element in (id2) has the same name as an element in the original
  ! (id1), then the values in the original are replaced with the values in
  ! the element from (id2).

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id1  ! id of list to be expanded
  integer, intent(in) :: id2  ! id of list to be added to list (id1)
  logical, intent(in) :: idelete  ! .TRUE. to delete (id2) when done.
  integer, intent(out) :: ierr    ! completion code (0=OK)

  !  ierr < 0 means a warning: there were duplicate names btw (id1) and
  !           (id2); (id2) list elements replaced.

  !----------------
  integer :: iertmp
  !-------------------------------------------------

  call xoi_author_set(iertmp)

  call xplasma_merge_lists(s,id1,id2,ierr, idelete)
  if(ierr.gt.0) then
     write(lunerr,*) ' ?eqm_list_merge: list merge failure, list ids:'
     write(lunerr,*) '  id1=',id1,'; id2=',id2
     call xplasma_error(s,ierr,lunerr)
  endif

  call xoi_author_clear(iertmp)

end subroutine eqm_list_merge

subroutine eqm_list_chval(ident,inum,zchvals,ierr)

  ! set list element character data for list #idlist

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: ident   ! list ID
  integer, intent(in) :: inum    ! number of elements in zchvals
  character*(*), dimension(inum), intent(in) :: zchvals  ! character data in
  integer, intent(out) :: ierr   ! completion code, 0=OK

  !-------------------------------
  integer :: iertmp
  !-------------------------------------------------

  call xoi_author_set(iertmp)

  call xplasma_setlist(s,ident,zchvals,ierr)
  if(ierr.ne.0) call eqm_list_errors(ident,ierr,'character string')

  call xoi_author_clear(iertmp)

end subroutine eqm_list_chval

subroutine eqm_list_ival(ident,inum,ivals,ierr)

  ! set list element integer data for list #idlist

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: ident   ! list ID
  integer, intent(in) :: inum    ! number of elements in ivals
  integer, dimension(inum), intent(in) :: ivals  ! integer data in
  integer, intent(out) :: ierr   ! completion code, 0=OK

  !-------------------------------
  integer :: iertmp
  !-------------------------------------------------

  call xoi_author_set(iertmp)

  call xplasma_setlist(s,ident,ivals,ierr)
  if(ierr.ne.0) call eqm_list_errors(ident,ierr,'integer')

  call xoi_author_clear(iertmp)

end subroutine eqm_list_ival

subroutine eqm_list_r8val(ident,inum,r8vals,ierr)

  ! set list element real*8 data for list #idlist

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: ident   ! list ID
  integer, intent(in) :: inum    ! number of elements in r8vals
  real*8, dimension(inum), intent(in) :: r8vals  ! real*8 data in
  integer, intent(out) :: ierr   ! completion code, 0=OK

  !-------------------------------
  integer :: iertmp
  !-------------------------------------------------

  call xoi_author_set(iertmp)

  call xplasma_setlist(s,ident,r8vals,ierr)
  if(ierr.ne.0) call eqm_list_errors(ident,ierr,'r8 floating point')

  call xoi_author_clear(iertmp)

end subroutine eqm_list_r8val

subroutine eqm_list_r4val(ident,inum,r4vals,ierr)

  ! set list element real*8 data for list #idlist (real data passed)

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: ident   ! list ID
  integer, intent(in) :: inum    ! number of elements in r8vals
  real, dimension(inum), intent(in) :: r4vals  ! floating point data in
  integer, intent(out) :: ierr   ! completion code, 0=OK

  !-------------------------------
  integer :: iertmp
  !-------------------------------------------------

  call xoi_author_set(iertmp)

  call xplasma_setlist(s,ident,r4vals,ierr)
  if(ierr.ne.0) call eqm_list_errors(ident,ierr,'r4 floating point')

  call xoi_author_clear(iertmp)

end subroutine eqm_list_r4val

subroutine eqm_list_errors(id,ierr,zstr)

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  !  error reporting routine for eqm_list.f90

  integer, intent(in) :: id   ! ID of list where error occurred
  integer, intent(in) :: ierr ! INPUT error status code
  character*(*), intent(in) :: zstr   ! type of data intended for insertion

  !-------------------------------
  character*32 znamel
  integer :: iertmp
  !-------------------------------

  if(ierr.ne.0) then
     call xplasma_get_item_info(s,id,iertmp, name=znamel)
     write(lunerr,*) &
          ' ?eqm_list_chval: error inserting ',trim(zstr),' data in list "', &
          trim(znamel),'":'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eqm_list_errors

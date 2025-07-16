subroutine eq_nlist(inum_lists)

  !  Return the number of known lists in the current XPLASMA.

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(out) :: inum_lists

  integer :: ierr

  call xplasma_num_items(s,ierr, num_lists=inum_lists)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_nlist: xplasma_num_lists returned ierr=',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_nlist

subroutine eq_nlist_p(inum_lists)

  !  Return the number of lists containing plottable data

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(out) :: inum_lists

  integer :: ierr

  call xplasma_num_items(s,ierr, num_plists=inum_lists)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_nlist: xplasma_num_lists returned ierr=',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_nlist_p

subroutine eq_lists(inum_lists,list_names,list_ids,ierr)

  !  Return sorted list of list names and associated list id numbers.
  !  This is a "list of lists".

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: inum_lists  ! max no. of lists (array dimension)
  character*(*), intent(out) :: list_names(inum_lists)  ! sorted names returned
  integer, intent(out) :: list_ids(inum_lists)  ! list ids returned
  integer, intent(out) :: ierr       ! completion code, 0=OK

  !--------------------------------------
  integer :: inum_act,inum_tot,i,ii,ilist,itype
  !--------------------------------------

  list_names = ' '
  list_ids = 0

  call eq_nlist(inum_act)
  if(inum_act.gt.inum_lists) then
     ierr=1
     write(lunerr,*) ' ?eq_lists: array dimension = ',inum_lists,' is not'
     write(lunerr,*) '  big enough for actual # of lists: ',inum_act
     return
  else
     ierr = 0
  endif

  if(inum_act.eq.0) return

  call xplasma_contents(s,ierr, id_lists=list_ids)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_lists: xplasma_contents returned ierr=',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

  do i=1,inum_act
     call xplasma_get_item_info(s,list_ids(i),ierr, name=list_names(i))
  enddo

end subroutine eq_lists

subroutine eq_lists_p(inum_lists,list_names,list_ids,ierr)

  !  Return sorted list of list names and associated list id numbers.
  !  ** return only lists containing plottable datasets **
  !  This is a "list of lists".

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: inum_lists  ! max no. of lists (array dimension)
  character*(*), intent(out) :: list_names(inum_lists)  ! sorted names returned
  integer, intent(out) :: list_ids(inum_lists)  ! list ids returned
  integer, intent(out) :: ierr       ! completion code, 0=OK

  !-------------------------------
  integer :: inum_act,inum_tot,i,ii,ilist,itype
  !--------------------------------------

  list_names = ' '
  list_ids = 0

  call eq_nlist_p(inum_act)
  if(inum_act.gt.inum_lists) then
     ierr=1
     write(lunerr,*) ' ?eq_lists_p: array dimension = ',inum_lists,' is not'
     write(lunerr,*) '  big enough for actual # of lists: ',inum_act
     return
  else
     ierr = 0
  endif

  if(inum_act.eq.0) return

  call xplasma_contents(s,ierr, id_plists=list_ids)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_lists: xplasma_contents returned ierr=',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

  do i=1,inum_act
     call xplasma_get_item_info(s,list_ids(i),ierr, name=list_names(i))
  enddo

end subroutine eq_lists_p

subroutine eq_glistnum(list_name,list_id)

  !  return id of a specific named list
  !  return 0 if no such list

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  character*(*), intent(in) :: list_name
  integer, intent(out) :: list_id

  !--------------------

  integer :: itype,ierr

  !--------------------

  call xplasma_listId(s,list_name,list_id)

end subroutine eq_glistnum

subroutine eq_list_size(id,list_size)

  !  return size of individual list specified by "id".
  !  return size=0 if no such list

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id    ! list id
  integer, intent(out) :: list_size

  !--------------------
  integer :: iertmp

  call xplasma_getlist(s,id,list_size,iertmp)
  if(iertmp.ne.0) then
     list_size = 0
  endif

end subroutine eq_list_size

subroutine eq_list_label(id,zlabel)

  !  return label of individual list specified by "id".
  !  return blank label if no such list

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id    ! list id
  character*(*), intent(out) :: zlabel

  !--------------------
  integer :: ilen,itype,ierr
  !--------------------

  zlabel=' '
  call xplasma_get_item_info(s,id,ierr, itype=itype)
  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
  else if(itype.ne.xplasma_listType) then
     write(lunerr,*) ' ?eq_list_label: item# ',id,' not a list.'
  else
     call xplasma_get_item_info(s,id,ierr, label=zlabel)
  endif

end subroutine eq_list_label

subroutine eq_list_maxsize(list_size)

  !  return the maximum size of any list in XPLASMA

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(out) :: list_size

  integer :: i,ierr,inum_tot,isize
  integer,dimension(:), allocatable :: id_lists

  list_size = 0

  !  cycle through all items in alphabetic order, get list item names only...

  call eq_nlist(inum_tot)
  if(inum_tot.eq.0) return

  allocate(id_lists(inum_tot))
  call xplasma_contents(s,ierr, id_lists=id_lists)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_lists: xplasma_contents returned ierr=',ierr
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  do i=1,inum_tot
     call eq_list_size(id_lists(i),isize)
     list_size=max(list_size,isize)
  enddo

  deallocate(id_lists)

end subroutine eq_list_maxsize

subroutine eq_list_enames(id,list_size,enames,ierr)

  !  return names of list elements in enames(1:list_size)
  !  return ierr.ne.0 & enames blank if list_size < actual size of list.

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id    ! list id
  integer, intent(in) :: list_size                 ! size of list passed in
  character*(*), intent(out) :: enames(list_size)  ! element names out
  integer, intent(out) :: ierr ! completion code, 0=OK

  !--------------------
  integer :: isize
  !--------------------

  call eq_list_size(id,isize)
  isize=min(max(1,isize),list_size)

  call xplasma_getlist_names(s,id,enames(1:isize),ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_list_enames -- access to list element names failed.'
     write(lunerr,*) '  list id: ',id,' ierr: ',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_list_enames

subroutine eq_list_chvals(id,list_size,chvals,ierr)

  !  return character data values of list elements in chvals(1:list_size)
  !  return ierr.ne.0 & chvals blank if list_size < actual size of list.

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id    ! list id
  integer, intent(in) :: list_size                 ! size of list passed in
  character*(*), intent(out) :: chvals(list_size)  ! character data values out
  integer, intent(out) :: ierr ! completion code, 0=OK

  !--------------------
  integer :: isize
  !--------------------

  call eq_list_size(id,isize)
  isize=min(max(1,isize),list_size)

  call xplasma_getlist(s,id,chvals(1:isize),ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_list_chvals-- access to list character data failed.'
     write(lunerr,*) '  list id: ',id,' ierr: ',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_list_chvals

subroutine eq_list_r8vals(id,list_size,r8vals,ierr)

  !  return floating pt data values of list elements in r8vals(1:list_size)
  !  return ierr.ne.0 & r8vals blank if list_size < actual size of list.

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id    ! list id
  integer, intent(in) :: list_size                 ! size of list passed in
  real*8, intent(out) :: r8vals(list_size)         ! character data values out
  integer, intent(out) :: ierr ! completion code, 0=OK

  !--------------------
  integer :: isize
  !--------------------

  call eq_list_size(id,isize)
  isize=min(max(1,isize),list_size)

  call xplasma_getlist(s,id,r8vals(1:isize),ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_list_r8vals-- access to list r8 float data failed.'
     write(lunerr,*) '  list id: ',id,' ierr: ',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_list_r8vals

subroutine eq_list_r4vals(id,list_size,r4vals,ierr)

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id    ! list id
  integer, intent(in) :: list_size                 ! size of list passed in
  real*4, intent(out) :: r4vals(list_size)         ! real*4 data values out
  integer, intent(out) :: ierr ! completion code, 0=OK
  
  !--------------------
  integer :: isize
  !--------------------

  call eq_list_size(id,isize)
  isize=min(max(1,isize),list_size)

  call xplasma_getlist(s,id,r4vals(1:isize),ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_list_r4vals -- access to list r4 float data failed.'
     write(lunerr,*) '  list id: ',id,' ierr: ',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_list_r4vals

subroutine eq_list_ivals(id,list_size,ivals,ierr)

  !  return integer data values of list elements in ivals(1:list_size)
  !  return ierr.ne.0 & ivals blank if list_size < actual size of list.

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id    ! list id
  integer, intent(in) :: list_size                 ! size of list passed in
  integer, intent(out) :: ivals(list_size)         ! integer data values out
  integer, intent(out) :: ierr ! completion code, 0=OK

  !--------------------
  integer :: isize
  !--------------------

  call eq_list_size(id,isize)
  isize=min(max(1,isize),list_size)

  call xplasma_getlist(s,id,ivals(1:isize),ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_list_ivals -- access to list integer data failed.'
     write(lunerr,*) '  list id: ',id,' ierr: ',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_list_ivals

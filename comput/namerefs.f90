module namerefs

  implicit NONE
  save

  !  expandable array of names
  character*64, dimension(:), pointer :: clist,clist_tmp

  !  work strings
  character*64 :: wk1,wk2,wk3

  !  expandable ordering array
  integer, dimension(:), pointer :: chord,chord_tmp

  !  expandable reference counter arrays
  integer, dimension(:), pointer :: nrefs,nrefs_tmp

  !  list size (in use, available space) (starts at zero)
  integer :: nsize = 0
  integer :: nmax = 0

  integer :: iupper = 1  ! set to zero to make upper/lower case distinct

contains

  subroutine expand

    integer :: init1

    if(nmax.eq.0) then
       init1=1          ! initial list -- length 64
       nmax=64
       allocate(clist(nmax))
       allocate(chord(nmax))
       allocate(nrefs(nmax))
    else
       init1=nmax       ! expand list by factor of 2
       nmax=4*nmax
       allocate(clist_tmp(nmax))
       allocate(chord_tmp(nmax))
       allocate(nrefs_tmp(nmax))
       clist_tmp(1:init1)=clist
       chord_tmp(1:init1)=chord
       nrefs_tmp(1:init1)=nrefs
       deallocate(clist,chord,nrefs)
       clist => clist_tmp
       chord => chord_tmp
       nrefs => nrefs_tmp
       nullify(clist_tmp)
       nullify(chord_tmp)
       nullify(nrefs_tmp)
       init1=init1+1
    endif

    clist(init1:nmax)=' '
    chord(init1:nmax)=0
    nrefs(init1:nmax)=0

  end subroutine expand

end module namerefs

subroutine eq_gflist(nlist,fnames,inames,ierr)

  !  for a list of function names, generate the corresponding set
  !  of function numbers.  set ierr if any names cannot be translated.

  IMPLICIT NONE

  integer nlist                     ! size of list (in)
  character*(*) fnames(nlist)       ! the function names (in)
  integer inames(nlist)             ! the function numbers (out)
  integer ierr                      ! exit code, 0=OK

  !  ierr is set, if any of the names in fnames are unknown.  If
  !  fnames(j) is unknown, inames(j)=0 is returned.

  !  this routine writes no messages.  The caller should check ierr.

  !------------------------------------------------------------------
  integer i
  !------------------------------------------------------------------

  ierr=0
  do i=1,nlist
     call eq_gfnum(fnames(i),inames(i))
     if(inames(i).eq.0) ierr=1
  enddo

end subroutine eq_gflist

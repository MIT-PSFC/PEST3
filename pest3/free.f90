subroutine pstFree

  ! Free the memory

  USE pstcom
  IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
  integer ier
  !
  call i2mex_free(ier)
  call i2mex_error(ier)

  call psttidyer

  CALL pstfree_mem(1)
  !CALL pstfree_mem(2)


  return
end subroutine pstFree

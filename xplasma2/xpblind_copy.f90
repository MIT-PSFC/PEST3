subroutine xpblind_copy(isize,afrom,ato)

  ! blind contiguous copy from array (afrom) to array (ato)

  implicit NONE
  integer, intent(in) :: isize
  real*8,  intent(in) :: afrom(isize)
  real*8, intent(out) :: ato(isize)

  ato = afrom

end subroutine xpblind_copy

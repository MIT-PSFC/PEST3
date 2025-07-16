Subroutine fort_flush(iu)

  ! flush output on sequential file on FORTRAN unit iu

  implicit none
  Integer, intent(in) :: iu

#if __IBM || __RS6000 || __AIX__ || HAVE_FC_INTRINSIC_FLUSH_
  Call flush_(iu)
#else
  Call flush(iu)
#endif

End subroutine fort_flush

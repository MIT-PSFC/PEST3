!...........................................................
      SUBROUTINE pstmatelt
!...........................................................
!
!      this routine organises the calculation of matrix elements...
!
 USE pstcom
 USE newmet

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)


!
      mbeg = 1
      nsing1 = nosing + 1
!
!      this version uses ke=true to omit singular functions...
!
      if ( .not. ke )  go to 5
      nsing1 = 1
      msing(2) = m + 1
    5 continue
!
      mend = m + 1
!
! regular-regular contribution to potential energy...
!
      CALL pstdelwrr
!
! scan over rational surfaces...
!
      do 101 ms = 1, nsing1
!
      mend = msing(ms + 1)
!
      if( .not. lsing ) go to 101
      if( ms == nsing1 ) go to 11
!
! compute prescribed big solution and associated forcing terms...
!
      CALL pstfrobe3

      CALL pstBigSolution
!
! regular-big solution contribution to energy...
!
      CALL pstdelwrb
!
! big-big solution contribution to energy...
!
      CALL pstdelwbb
!
 11   continue
!
      mbeg = mend + 1
!
 101  continue
!
      return
      end



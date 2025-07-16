!      1.1_r8    defmod
!.......................................................................
      SUBROUTINE pstdefmod
!.......................................................................
 USE pstcom

 USE l22com

 USE mtrik1
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER I
      logical :: lreal


!      1.1.3   default values.
!
!
!
!      set for regular ideal mhd...
!
      nosing = 0
      lsub = .false.
do i = 1, nfe
         msub(i) = 2
enddo
!
      return
      end



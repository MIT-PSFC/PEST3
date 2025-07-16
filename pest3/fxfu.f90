         SUBROUTINE pstfxfu
!
 USE pstcom

 USE combla

 USE comivi
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER I
      INTEGER IX
      INTEGER J


!
!
         do 1 i=1,mf
         x(i)= 1._r8 
    1    continue
      CALL pstblkcpy(x(1),xt(1),mf)
!
!
         ix=mf+1
         do 2 i=1,ng
         do 3 j=1,ml
         x(j)= 1._r8 
    3    continue
      CALL pstblkcpy( x(1),xt(ix),ml)
         ix=ix+ml
    2    continue
!
!
!
         return
         end



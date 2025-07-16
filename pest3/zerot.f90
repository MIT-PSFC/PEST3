!--------------------------------------------------------------
      SUBROUTINE pstzerot( ileft,ir )
!--------------------------------------------------------------
 USE pstcom

 USE r33com

 USE l22com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER ILEFT
      INTEGER IR
      INTEGER I


!
      do i = ileft, ir
         tw(i,:,:) =  0.0_r8 
         ty(i,:,:) =  0.0_r8 
         tz(i,:,:) =  0.0_r8 
         td(i,:,:) =  0.0_r8 
      end do
!
      return
      end



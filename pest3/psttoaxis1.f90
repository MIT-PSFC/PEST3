subroutine pstToAxis1(si, zinout)
!
! Cubic extrpolation to axis
!
! pletzer@pppl.gov Wed May 10 09:12:24 EDT 2000

use pstcom
implicit none
integer, parameter :: r8 = selected_real_kind(12,100)
real(r8), intent(in) :: si(*)
real(r8), intent(inout) :: zinout(*)

include 'CUCCCC.inc'

       zinout(1) = FCCCC0(zinout(2), zinout(3), zinout(4), zinout(5), &
     &                 si(2),  si(3),  si(4),  si(5), &
     & si(1))

end subroutine pstToAxis1

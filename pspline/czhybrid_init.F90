! -*-f90-*-
! $Id: czhybrid_init.F90 1042 2011-11-08 21:47:39Z pletzer $
!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Initialization of ezspline

#include "cpp/czspline_handle_size.h"

#define _S czhybrid_init2_r4
#define _O czspline2_r4
subroutine _S(handle, n1, n2, hspline, BCS1, BCS2, ier)
   use ezspline_obj
   use ezspline
   use czspline_pointer_types
   implicit none
   integer, intent(inout) :: handle(_ARRSZ)
   integer, intent(in) :: n1, n2
   integer, intent(in) :: hspline(2)
   integer, intent(in) :: BCS1(2), BCS2(2)
   integer, intent(out) :: ier
   type(_O) :: self
   allocate(self % ptr)
   call ezhybrid_init(self % ptr, n1, n2, hspline, ier, BCS1, BCS2)
   handle = 0
   handle = transfer(self, handle)
end subroutine _S
#undef _S
#undef _O

#define _S czhybrid_init3_r4
#define _O czspline3_r4
subroutine _S(handle, n1, n2, n3, hspline, BCS1, BCS2, BCS3, ier)
   use ezspline_obj
   use ezspline
   use czspline_pointer_types
   implicit none
   integer, intent(inout) :: handle(_ARRSZ)
   integer, intent(in) :: n1, n2, n3
   integer, intent(in) :: hspline(3)
   integer, intent(in) :: BCS1(2), BCS2(2), BCS3(2)
   integer, intent(out) :: ier
   type(_O) :: self
   allocate(self % ptr)
   call ezhybrid_init(self % ptr, n1, n2, n3, hspline, ier, BCS1, BCS2, BCS3)
   handle = 0
   handle = transfer(self, handle)
end subroutine _S
#undef _S
#undef _O

#define _S czhybrid_init2_r8
#define _O czspline2_r8
subroutine _S(handle, n1, n2, hspline, BCS1, BCS2, ier)
   use ezspline_obj
   use ezspline
   use czspline_pointer_types
   implicit none
   integer, intent(inout) :: handle(_ARRSZ)
   integer, intent(in) :: n1, n2
   integer, intent(in) :: hspline(2)
   integer, intent(in) :: BCS1(2), BCS2(2)
   integer, intent(out) :: ier
   type(_O) :: self
   allocate(self % ptr)
   call ezhybrid_init(self % ptr, n1, n2, hspline, ier, BCS1, BCS2)
   handle = 0
   handle = transfer(self, handle)
end subroutine _S
#undef _S
#undef _O

#define _S czhybrid_init3_r8
#define _O czspline3_r8
subroutine _S(handle, n1, n2, n3, hspline, BCS1, BCS2, BCS3, ier)
   use ezspline_obj
   use ezspline
   use czspline_pointer_types
   implicit none
   integer, intent(inout) :: handle(_ARRSZ)
   integer, intent(in) :: n1, n2, n3
   integer, intent(in) :: hspline(3)
   integer, intent(in) :: BCS1(2), BCS2(2), BCS3(2)
   integer, intent(out) :: ier
   type(_O) :: self
   allocate(self % ptr)
   call ezhybrid_init(self % ptr, n1, n2, n3, hspline, ier, BCS1, BCS2, BCS3)
   handle = 0
   handle = transfer(self, handle)
end subroutine _S
#undef _S
#undef _O

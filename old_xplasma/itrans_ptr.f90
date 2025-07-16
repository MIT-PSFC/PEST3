subroutine itrans_ptr(isize,iarray)

  ! convert the old_xplasma instance pointer "s" to an integer array
  ! this allows the pointer value to be held in routines that do not
  ! use xplasma definition modules.

  ! the caller must provide an integer array of sufficient size.
  ! for codes using Plasma State definitions, the parameter ps_container_size
  ! defines a safe value.

  use xplasma_obj_instance
  implicit NONE

  integer, intent(in) :: isize   ! integer array size
  integer :: iarray(isize)       ! array to receive pointer data

  !---------------------------
  type :: interp_object

     type (xplasma), pointer :: s_xpobj => NULL()

  end type interp_object
  type (interp_object) :: xobj
  !---------------------------

  xobj%s_xpobj => s

  iarray = 0
  iarray = transfer(xobj, iarray)

end subroutine itrans_ptr

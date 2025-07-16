module xplasma_obj_instance_mod

  ! declare an instance of an xplasma object; make this static, shared
  ! global instance publicly available.

  use xplasma_obj

  implicit NONE
  save
  private

  !-----------------------
  ! # of Fourier moments
  !   in Fourier spline equilibrium representation

  integer :: kmom = 16

  !-----------------------
  ! base object pointer

  type (xplasma), public, pointer :: s => NULL()

  ! a couple of spares...

  type (xplasma), public, pointer :: s1 => NULL()
  type (xplasma), public, pointer :: s2 => NULL()

  !----------------------------------
  ! a 100-high stack (private)

  private :: ptr_container

  type :: ptr_container 
     type(xplasma), pointer :: s => NULL()
  end type ptr_container

  integer, parameter, private :: nmax_sstack = 100
  type (ptr_container), dimension(nmax_sstack), private :: sstack
  integer, private :: ncur_sstack = 0

  !-----------------------
  ! private flag if "F77_USER" has been asserted as an author

  logical :: iflag_f77_user = .FALSE.
  
  !-----------------------
  !  the objects themselves...
  !   (these are large objects and so just the pointers will be passed
  !   between subroutine calls).

  type (xplasma), target :: s_obj,s1_obj,s2_obj

  !-----------------------
  
  logical, private :: xplasma_obj_init = .FALSE.

  public :: xoi_init, grab_f77_xplasma, release_f77_xplasma
  public :: xoi_kmom_get, xoi_kmom_set
  public :: xoi_author_set, xoi_author_clear

  contains

    subroutine xoi_init(ier)

      !  check if xplasma instances and pointers are initialized...

      integer, intent(out) :: ier  ! completion code

      integer :: iertmp
      iertmp = 0

      ier = 0
      if(.not.xplasma_obj_init) then

         s => s_obj
         s1 => s1_obj
         s2 => s2_obj

         xplasma_obj_init = .TRUE.

      endif

    end subroutine xoi_init

    subroutine grab_f77_xplasma(swant,ierr)

      ! reset f77 xplasma pointer to point to "swant"; push old pointer
      ! onto stack.  Possible error: stack overflow.

      ! calls to this routine should be paired with "release_f77_xplasma"
      ! calls.

      type (xplasma), pointer :: swant
      integer, intent(out) :: ierr

      !------------------------

      if(.not.xplasma_obj_init) then
         call xoi_init(ierr)
         if(ierr.ne.0) return
      endif

      if(ncur_sstack.eq.nmax_sstack) then
         ierr = 1  ! stack overflow

      else
         ierr = 0
         ncur_sstack = ncur_sstack + 1
         sstack(ncur_sstack)%s => s
         s => swant
      endif

    end subroutine grab_f77_xplasma

    subroutine release_f77_xplasma(ierr)

      ! restore prior xplasma pointer
      ! possible error: stack underflow.  Calls to this routine should
      ! be paired with "grab_f77_xplasma" calls...

      integer, intent(out) :: ierr

      !------------------------

      if(.not.xplasma_obj_init) then
         call xoi_init(ierr)
         if(ierr.ne.0) return
      endif

      if(ncur_sstack.eq.0) then
         ierr = 1  ! stack underflow

      else
         ierr = 0
         s => sstack(ncur_sstack)%s
         ncur_sstack = ncur_sstack - 1
      endif

    end subroutine release_f77_xplasma

    !----------------------------------------------------

    subroutine xoi_kmom_get(ikmom)

      integer, intent(out) :: ikmom

      ikmom = kmom

    end subroutine xoi_kmom_get

    subroutine xoi_kmom_set(ikmom)

      integer, intent(in) :: ikmom

      kmom = max(2,min(64,ikmom))

    end subroutine xoi_kmom_set

    !----------------------------------------------------

    subroutine xoi_author_set(ierr)

      !  see if a write-enabled author is defined; if not, define as F77_USER

      integer, intent(out) :: ierr

      character*32 :: zauth
      logical :: ienable

      call xplasma_author_query(s,zauth,ienable,ierr)
      if(ierr.ne.0) return

      if(.not.ienable) then
         iflag_f77_user = .TRUE.
         call xplasma_author_set(s,'F77_USER',ierr)
      endif

    end subroutine xoi_author_set

    subroutine xoi_author_clear(ierr)

      ! clear F77_USER author, if flag is set

      integer, intent(out) :: ierr

      ierr = 0

      if(iflag_f77_user) then
         iflag_f77_user = .FALSE.
         call xplasma_author_clear(s,'F77_USER',ierr)
      endif

    end subroutine xoi_author_clear
         
end module xplasma_obj_instance_mod

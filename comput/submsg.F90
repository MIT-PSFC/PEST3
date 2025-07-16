
!
! defines a submsg type which can be used for returning a multi-line error message from
! a fortran 90 subroutine.  A default submsg is defined with C accessor functions.
!
! Mods
!    22Dec2009    Jim.Conboy@cffe.ac.uk
!                 cstring(trim(),  => cstr_2C(trim(), 
!                 [ Former invalid under debug compilation ]
!
module submsg_mod
  implicit none

  integer, parameter     :: nsubmsg = 132   ! size of one line of the msg
  character(len=nsubmsg) :: mline           ! single line for temporary use in external code

  !
  ! --- submsg type ---
  !
  type :: submsg
     private 
     integer :: num = -1                                          ! number of logical lines or -1 for unallocated memory
     character*(nsubmsg), pointer, dimension(:) :: msg => null()  ! the message
  end type submsg

  !
  ! --- default submsg ---
  !
  logical,private,save              :: igot=.false.   ! true when default_msg has been constructed
  type(submsg),pointer,private,save :: default_msg    ! default message for passing to C routines

  !
  ! --- interfaces ---
  !
  interface add_msg
     module procedure add_msg_char, add_msg_msg
  end interface

  interface copy_msg
     module procedure copy_msg_char, copy_msg_msg
  end interface

contains
  !
  ! ----- new_msg -------
  ! initializes a new submsg type -- use this right after declaration or allocation.
  ! The optional 'ns' argument is a desired storage size -- the default is unallocated.
  !
  ! Update: not needed anymore because of fortran 95 initialization of types
  !
  subroutine new_msg(m, ns)
    type(submsg), intent(inout)   :: m   ! message
    integer, optional, intent(in) :: ns  ! optional initial size as the number of lines

    if (present(ns)) then
       if (ns>0) call resize_msg(m,ns)
    end if
  end subroutine new_msg

  !
  ! ------- size_msg ------
  ! return the logical size of the message
  !
  function size_msg(m) result(n)
    type(submsg),intent(in) :: m    ! message
    integer :: n                    ! logical size

    n = max(0, m%num)    ! unallocated returns 0
  end function size_msg

  !
  ! ---- resize_msg ----
  ! force the internal storage of a submsg to a fixed size.  The logical
  ! size in submsg%num is not changed unless the physical size drops
  ! below the logical size.  If the ns size argument is nonnegative, use this
  ! as the final size of the storage.  If the ns size argument is negative, ensure
  ! the storage has at least |ns| free spaces.
  !
  subroutine resize_msg(m, ns)
    type(submsg),intent(inout) :: m    ! message
    integer,     intent(in)    :: ns   ! ns>=0 -> new size of the storage
                                       ! ns<0  -> need at least |ns| free spaces

    integer :: n                                     ! required size
    integer :: k                                     ! overlap between old and new storage
    character*(nsubmsg), pointer, dimension(:) :: p  ! new storage

    ! find n
    if (ns>=0) then
       n = ns
    else
       if (m%num<=0) then
          n = abs(ns)
       else
          n = m%num + abs(ns)            ! minimum required size of storage
          if (n<=size(m%msg)) return     ! there is at least |ns| free lines
          n = max(n+3,4*n/3)             ! need to increase the storage with some slop
       end if
    end if
          
    if (m%num<0) then
       !print *, 'allocated new ',n
       allocate(m%msg(max(0,n)))     ! get at least something
       m%num=0                       ! logically empty
    else
       if (n==size(m%msg)) return    ! nothing to do
       allocate(p(n))
       !print *, 'allocated ',n
       k = min(size(p),m%num) 
       p(1:k) = m%msg(1:k)           ! copy old lines over
       deallocate(m%msg)             ! free up old storage
       m%msg => p
       m%num = k                     ! might chop lines off
    end if
  end subroutine resize_msg

  !
  ! ------ trim_msg ------
  ! remove lines by reducing the logical size of the message
  !
  subroutine trim_msg(m, ns)
    type(submsg), intent(inout) :: m    ! message
    integer,      intent(in)    :: ns   ! number of lines to trim

    if (m%num<=0 .or. ns<=0) return
    m%num = max(0,m%num-ns)
  end subroutine trim_msg

  !
  ! ------- clear_msg --------
  ! empty the submsg deallocating the memory if necessary
  !
  subroutine clear_msg(m)
    type(submsg), intent(inout) :: m    ! message

    if (m%num>=0) then
       deallocate(m%msg)   ! only when allocated
       m%num=-1
    end if
  end subroutine clear_msg

  !
  ! -------- add_msg_char --------
  ! append a new string to the message
  !
  subroutine add_msg_char(m, c)
    type(submsg),  intent(inout) :: m    ! message
    character*(*), intent(in)    :: c    ! string to add to message
    integer :: n   ! new size

    call resize_msg(m,-1)
    n = m%num+1
    m%msg(n) = c
    m%num = n
  end subroutine add_msg_char


  !
  ! -------- add_msg_msg --------
  ! append all the strings in one message to the target
  !
  subroutine add_msg_msg(m, mc)
    type(submsg), intent(inout) :: m     ! target message
    type(submsg), intent(in)    :: mc    ! source of strings to add to message

    integer :: i,n     ! temps

    if (mc%num<=0) return    ! empty source
    call resize_msg(m,-mc%num)

    do i=1, mc%num
       n = m%num+1
       m%msg(n) = mc%msg(i)
       m%num = n
    end do
  end subroutine add_msg_msg

  !
  ! -------- copy_msg_char --------
  ! set the message argument to the single string
  !
  subroutine copy_msg_char(m, c)
    type(submsg),  intent(inout) :: m     ! target message
    character*(*), intent(in)    :: c     ! string to add to message

    if (m%num<=0) call resize_msg(m,-1)
    m%msg(1) = c
    m%num = 1
  end subroutine copy_msg_char

  !
  ! -------- copy_msg_msg --------
  ! set the message argument to the source message
  !
  subroutine copy_msg_msg(m, mc)
    type(submsg),  intent(inout) :: m     ! target message
    type(submsg),  intent(in)    :: mc    ! source of strings to copy to message

    integer :: i    ! loop variable

    if (mc%num<=0) then
       if (m%num>0) m%num=0    ! gone
    else
       if (m%num<mc%num) call resize_msg(m,mc%num)
       do i=1,mc%num
          m%msg(i) = mc%msg(i)
       end do
       m%num = mc%num
    end if
  end subroutine copy_msg_msg

  !
  ! ------- line_msg -------
  ! return a single line in the msg
  !
  subroutine line_msg(m, i, line)
    type(submsg), intent(in)    :: m     ! target message
    integer,      intent(in)    :: i     ! desired line
    character*(*),intent(inout) :: line  ! return line here

    integer :: n                     ! number of available lines
    character*(nsubmsg) :: errmsg    ! errormsg
    
    n = max(0,m%num)
    if (i<=0 .or. i>n) then
       write(errmsg,'(a,i4,a,i4,a)') '?line_msg: index ',i,' is out of the fortran range [1:',n,']'
       line = errmsg
    else
       line = trim(m%msg(i))
    end if
  end subroutine line_msg

  !
  ! ------- write_msg -------
  ! write a message to a logical unit number
  !
  subroutine write_msg(m, n)
    type(submsg), intent(in) :: m     ! target message
    integer,      intent(in) :: n     ! logical unit number for write
    integer  :: i   ! loop variable

    if (m%num<=0) return
    do i = 1, m%num
       write(n,'(a)') trim(m%msg(i))
    end do
  end subroutine write_msg

  !
  ! ------- print_msg -------
  ! write a message using the print statment
  !
  subroutine print_msg(m)
    type(submsg), intent(in) :: m     ! target message
    integer  :: i   ! loop variable

    if (m%num<=0) return
    do i = 1, m%num
       print '(a)', trim(m%msg(i))
    end do
  end subroutine print_msg

  !
  ! ---------- get_defmsg() ----------
  ! returns a pointer to the default message
  !
  function get_defmsg() result(m)
    type(submsg), pointer :: m

    if (.not. igot) then
       allocate(default_msg)
       call new_msg(default_msg)
       igot = .true.
    end if
    m => default_msg
  end function get_defmsg

  !
  ! ---------- replace_defmsg() ----------
  ! replace the default message with the values in the argument
  !
  subroutine replace_defmsg(m) 
    type(submsg), intent(in) :: m

    if (.not. igot) then
       allocate(default_msg)
       call new_msg(default_msg)
       igot = .true.
    end if
    call clear_msg(default_msg)
    call copy_msg(default_msg,m)
  end subroutine replace_defmsg

end module submsg_mod


!
! =================================== C ============================
!
! C callable routines for getting to the default message returned by 
! get_defmsg() in submsg_mod.  The way to structure the C interface to
! your fortran code is as follows,
!   subroutine c_mycode(ier)
!     use submsg_mod
!
!     integer :: ier
!     type(submsg), pointer :: errmsg
!
!     errmsg => get_defmsg()
!     call mycode(errmsg, ier)
!   end subroutine c_mycode
!
! And then the C++ caller would
!   void do_mycode() {
!     int ier=0 ;
!     c_mycode(ier) ;
!     if (ier!=0) c_print_defmsg() ;   // or construct an error message and throw an exception
!   } ;
!
! ---- c_clear_defmsg ---
! extern "C" void F77NAME(c_clear_defmsg)() ;
!
! clear out the default message
!
subroutine c_clear_defmsg
  use submsg_mod
  implicit none
  type(submsg),pointer :: d    ! default message

  d => get_defmsg()
  call clear_msg(d)
end subroutine c_clear_defmsg

!
! ---- c_print_defmsg ---
! extern "C" void F77NAME(c_print_defmsg)() ;
!
! print out the default message
!
subroutine c_print_defmsg
  use submsg_mod
  implicit none
  type(submsg),pointer :: d    ! default message

  d => get_defmsg()
  call print_msg(d)
end subroutine c_print_defmsg

!
! ---- c_size_defmsg ----
! extern "C" int F77NAME(c_size_defmsg)() ;
!
! return the logical size of the default message
!
function c_size_defmsg() result(n)
  use submsg_mod
  implicit none
  integer :: n                 ! size of the default
  type(submsg),pointer :: d    ! default message

  d => get_defmsg()
  n = size_msg(d)
end function c_size_defmsg

!
! ---- c_size_line ----
! extern "C" int F77NAME(c_size_line)() ;
!
! return the maximum line size nsubmsg
!
function c_size_line() result(n)
  use submsg_mod
  implicit none
  integer :: n                 ! size of the default
  n = nsubmsg
end function c_size_line

!
! ----- c_resize_defmsg ----
! extern "C" void F77NAME(c_resize_defmsg)(int*) ;
!
! resize the default message, see resize_msg() subroutine
!
subroutine c_resize_defmsg(ns)
  use submsg_mod
  implicit none
  integer :: ns                        ! ns>=0  -> new size of the storage
                                       ! ns<0   -> need at least |ns| free spaces
  type(submsg),pointer :: d            ! default message

  d => get_defmsg()
  call resize_msg(d, ns)
end subroutine c_resize_defmsg

!
! ----- c_trim_defmsg ----
! extern "C" void F77NAME(c_trim_defmsg)(int*) ;
!
! remove lines from the default message
!
subroutine c_trim_defmsg(ns)
  use submsg_mod
  implicit none
  integer :: ns                        ! number of lines to remove
  type(submsg),pointer :: d            ! default message

  d => get_defmsg()
  call trim_msg(d, ns)
end subroutine c_trim_defmsg


!
! ----- c_line_defmsg ----
! extern "C" void F77NAME(c_line_defmsg)(int*, int*, char*) ;
!
! fetch the i'th line in the default message -- uses "C" indexing !!!
!
subroutine c_line_defmsg(i, n, c)
  use submsg_mod
  implicit none
#include "fpreproc/byte_declare.h"
  integer       :: i       ! C index of line which is desired
  integer       :: n       ! size of character buffer not including end 0, must be >0
  BYTE_DECLARE  :: c(*)    ! will contain null terminated string on output

  character*(nsubmsg)  :: p     ! buffer
  type(submsg),pointer :: d     ! default message

  if (n<=0) return

  d => get_defmsg()
  call line_msg(d,i+1,p)
  call cstr_2C(trim(p(1:n)), c, "2C")
end subroutine c_line_defmsg

!
! ---- c_add_defmsg ---
! extern "C" void F77NAME(c_add_defmsg)(char*) ;
!
! add a null terminated character string to the default message
!
subroutine c_add_defmsg(c)
  use submsg_mod
  implicit none
#include "fpreproc/byte_declare.h"
  BYTE_DECLARE :: c(*)    ! will contain null terminated string on output

  character*(nsubmsg)  :: p     ! buffer
  type(submsg),pointer :: d     ! default message

  call cstring(p,c,"2F")

  d => get_defmsg()
  call add_msg(d,p)
end subroutine c_add_defmsg

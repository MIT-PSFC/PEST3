      subroutine sset_env(varname,value,ierr)
C
C  set the value of an environment variable (UNIX)
C
C Mods
C    22Dec2009    Jim.Conboy@cffe.ac.uk
C                 remove trim from cstring arg list ( fails debug compilation )
C
      implicit NONE
C
      character*(*), intent(in) :: varname
      character*(*), intent(in) :: value
      integer, intent(out) :: ierr    ! completion code (0=OK)
C
C-------------------------
#include "fpreproc/byte_declare.h"
C
      integer :: ilenvar,ilenval
C
      BYTE_DECLARE, dimension(:), allocatable :: cvar,cval
      integer :: f77_setenv
C-------------------------
C
      ilenvar = len(trim(varname))
      ilenval = len(trim(value))
C
      allocate(cvar(ilenvar+1),cval(ilenval+1))  ! leave room for trailing null
C
C-      call cstring(trim(varname),cvar,'2C')
C-      call cstring(trim(value),cval,'2C')
C
      call cstring(varname(:ilenvar),cvar,'2C')
      call cstring(value(:ilenval)  ,cval,'2C')
C
      ierr = f77_setenv(cvar,cval)
      ierr = abs(ierr)
C
      return
      end

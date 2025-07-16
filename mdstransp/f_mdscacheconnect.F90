#include "fpreproc/byte_declare.h"

!
! fortran version of MdsCacheConnect(char*)
!
function MdsCacheConnect(server) result(istat)
  character*(*), intent(in) :: server   ! server for cached connection
  
  integer      :: istat                        ! socket
  BYTE_DECLARE :: cserver(len_trim(server)+1)  ! C string

  integer, external :: c_MdsCacheConnect

  call cstring(trim(server), cserver, '2C')

  istat = c_MdsCacheConnect(cserver)
end function MdsCacheConnect


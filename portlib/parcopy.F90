      subroutine parcopy(npcopy_comm,myid,read_id,write_id,inputfile, &
           outputfile,ierr)

      implicit none

      integer :: npcopy_comm  ! MPI communicator
      integer :: myid         ! current process ID
      integer :: read_id      ! ID of process that should READ file
      integer :: write_id     ! ID of process that should WRITE file

      character*(*),intent(in) :: inputfile,outputfile  ! filenames

      integer, intent(out) :: ierr

      !---------------
      integer :: portlib_parcopy
      logical :: idebug = .FALSE.
      !---------------

#include "fpreproc/byte_declare.h"

      BYTE_DECLARE cinput(1+len(inputfile))
      BYTE_DECLARE coutput(1+len(outputfile))

      integer str_length

      if(idebug) then
         if(myid.eq.read_id) then 
            write(0,*) ' I, ',myid,' am reader: '
            write(0,*) '     inputfile: '//trim(inputfile)
            write(0,*) '    outputfile: '//trim(outputfile)
         endif
         if(myid.eq.write_id) then 
            write(0,*) ' I, ',myid,' am writer: '
            write(0,*) '     inputfile: '//trim(inputfile)
            write(0,*) '    outputfile: '//trim(outputfile)
         endif
      endif

      call cstring(inputfile(1:str_length(inputfile)),cinput,'2C')
      call cstring(outputfile(1:str_length(outputfile)),coutput,'2C')

      ierr = portlib_parcopy(npcopy_comm,myid,read_id,write_id,cinput,coutput)

      return 
!     stop
      end
!


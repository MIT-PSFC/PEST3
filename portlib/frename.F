      subroutine frename_OK(oldname,newname,ier)
C
C  dmc 30 June 2011 -- rename a file if it exists & can be opened; otherwise
C  do nothing & return 0 error status 
C
      implicit NONE
C
      character*(*), intent(in) :: oldname   ! old filename
      character*(*), intent(in) :: newname   ! new filename
      integer, intent(out) :: ier
C
C-------------------------------------------
      integer :: io,istat
C-------------------------------------------
      ier=0
C
      call find_io_unit(io)
      open(unit=io,file=oldname,status='old',iostat=istat)
      if(istat.ne.0) then
                                ! return silently
         ier=0
         return
      endif

      !  OK, file is accessible; proceed with rename
      close(unit=io)
      call frename(oldname,newname,ier)

      return
      end

C======================================================
C
      subroutine frename(oldname,newname,ier)
C
C  dmc 18 June 1996 -- rename file
C
C  encapsulate VMS/UNIX OS dependence here.
C
      implicit NONE

      character*(*), intent(in) :: oldname   ! old filename
      character*(*), intent(in) :: newname   ! new filename
      integer, intent(out) :: ier  
C
C  ** rename a file **
C
C  return error code
C    ier=0 -- success
C    ier=1 -- error
C
      integer str_length,ilo,iln,istat
C
#ifdef __UNIX
#include "fpreproc/byte_declare.h"
      BYTE_DECLARE fnbuf1(500),fnbuf2(500)
      integer cdelete,crename,cccopy
#endif
C----------------------------------
C
      ier=0
C
      ilo=str_length(oldname)
      iln=str_length(newname)
      if(min(ilo,iln).le.0) then
         write(6,*) ' ?frename:  filename argument is blank.'
         ier=1
         return
      endif
C
      if(oldname.eq.newname) then
         write(6,*) ' %frename:  oldname & newname match: '//
     >        trim(oldname)
         write(6,*) '  ...no action taken.'
         ier=0
         return
      endif
C
C  UNIX code:
C
#ifdef __UNIX
      call cstring(oldname(1:ilo),fnbuf1,'2C')
      call cstring(newname(1:iln),fnbuf2,'2C')
      ier=cdelete(fnbuf2)
      ier=crename(fnbuf1,fnbuf2)
      if(ier.ne.0) then
         write(6,*) ' %frename: retry as copy & delete...'
         ier=cccopy(fnbuf1,fnbuf2)
         if(ier.eq.0) then
            ier=cdelete(fnbuf1)
            if(ier.ne.0) then
               write(6,*) ' ?frename: delete of oldname failed: '//
     >              trim(oldname)
            endif
         else
            write(6,*) ' ?frename: copy failed: '
            write(6,*) '    oldname: '//trim(oldname)
            write(6,*) '    newname: '//trim(newname)
         endif
      endif
#endif
C
C  VMS code:
C
#ifdef __VMS
      istat=lib$rename_file(oldname(1:ilo),newname(1:iln))
      if(istat.ne.1) then
        ier=1
      else
        ier=0
      endif
#endif
C
      return
      end

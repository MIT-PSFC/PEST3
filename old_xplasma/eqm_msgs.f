      subroutine eqm_msgs(ilun,zfile,ierr)
C
      use eq_module
C
C  set error message i/o unit, and, optionally, open a file
C  *** caller should check error code returned ***
C
      IMPLICIT NONE
C
C  input:
C
      integer ilun                      ! unit on which to write messages
      character*(*) zfile               ! (optional) file for messages
C
C  output:
C
      integer ierr                      ! error code (open iostat, 0=normal)
C
c$cface_input:  zfile
C------------------------------------
C
      ierr=0
      lunerr=ilun
C
      if(zfile.ne.' ') then
         open(unit=lunerr,file=zfile,status='unknown',iostat=ierr)
      endif
C
      return
      end

C  open a file for write using UFOPEN.
C  on a unix system, first remove any existing file that might be
C  in the way.
C
C  careful, misuse of this routine can be dangerous to your health
C  it deletes files!
C
      subroutine xfopen(ilun,zfile,ier,ibin)
C
      implicit NONE
C
      integer ilun               ! fortran logical unit number (input)
      character*(*) zfile        ! filename (input)
      integer ier                ! ufopen completion code (output)
      integer ibin               ! ufopen format argument (input)
C
C----------------
C  UNIX code:
C
#ifdef __UNIX
      open(unit=ilun,file=zfile,status='old',err=10)
C
C  no error:  delete file!
C
      close(unit=ilun,status='DELETE')
 10   continue
#endif
C
C  VMS & UNIX code:
C   hardwire "2" for write access only.
C
      call ufopen(ilun,zfile,2,ier,ibin)
C
      return
      end
 

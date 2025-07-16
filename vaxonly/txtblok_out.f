      subroutine txtblok_out(lun,filename,txtarray,nsize,istat)
c
c  write a file from a text array in memory out to disk
c
      implicit NONE
c
c input:
      integer lun                       ! logical unit no. for i/o
      character*(*) filename            ! filename to write
c
      integer nsize                     ! size of text arrary
      character*(*) txtarray(nsize)     ! text array
c
c output:
c
      integer istat                     ! status code on exit
c
c  istat = 0  -- normal
c
c  istat = 1  -- file open error
c
      integer ier,i,iln,ilnurd
c
c----------------------------------------
c
      istat=0
c
      call ufopen(lun,filename,2,ier,0)
      if(ier.ne.0) then
         istat=1
         return
      endif
c
      do i=1,nsize
         iln=ilnurd(txtarray(i))
         write(lun,'(A)') txtarray(i)(1:iln)
      enddo
c
      close(unit=lun)
      return
      end

      subroutine txtblok_in(lun,filename,txtarray,nmax,ngot,istat,
     >   maxwid)
c
c  read a file into a text array in memory
c
      implicit NONE
C
c input:
      integer lun                       ! input logical unit no.
      character*(*) filename            ! file to open and read
c
      integer nmax                      ! size of text array
c
c output:
      character*(*) txtarray(nmax)      ! array to hold text read
c
      integer ngot                      ! number of text lines read
c
      integer istat                     ! output status code (see below)
c
      integer maxwid                    ! maximum non-blank line width
c
c on exit:
c
c  istat=0  --  normal completion; entire file fits in array.
c   ...ngot=<# of lines read in>
c  istat=1  --  file open error
c   ...ngot=0
c  istat=2  --  file read error
c   ...ngot=<# of lines read in prior to error>
c
c  istat=10 --  error:  txtarray too short to contain entire file
c   ...ngot=nmax
c
c  istat=11 --  warning:  truncation:  text line in file is wider than
c               txtarray array element can store
c   ...ngot=<# of lines read in>
c
      integer ier,iline,ilb,ilnurd
      character*180 ibuff               ! input line buffer
c
c----------------------------
c
      ngot=0
      istat=0
      maxwid=0
c
      call ufopen(lun,filename,1,ier,0)
c
      if(ier.ne.0) then
         istat=1
         go to 1100
      endif
c
      iline=len(txtarray(nmax))
c
c  read loop
c
 10   continue
      read(lun,'(A)',end=80,err=90) ibuff
c
      ilb=ilnurd(ibuff)
      maxwid=max(ilb,maxwid)
      if(ilb.gt.iline) then
         istat=11
      endif
c
      ilb=min(iline,ilb)
c
      if(ngot.eq.nmax) then
         istat=10
         go to 1000
      endif
 
      ngot=ngot+1
      txtarray(ngot)=ibuff(1:ilb)
c
      go to 10
c
c------------------------
c
 80   continue
c
c  end of file (normal)
c
      go to 1000
c
 90   continue
c
c  READ ERROR
c
      istat=2
      go to 1000
c
c---------------------
c
c  exit
c
 1000 continue
      close(unit=lun)
c
 1100 continue
      return
      end

      subroutine mkfrid(runid, start, end)
C
C     Get start and end of runid
C
C     11/04/96 CAL
C
      implicit NONE
C
      character*(*) runid
      integer       start, end
 
      integer       ic, l
 
      l = len(runid)
      ic = 1
      do 10 while (runid(ic:ic) .eq. ' ' .and. ic .lt. l)
        ic = ic+1
 10   continue
      start = ic
C
      if (runid(l:l) .ne. ' ') then
         end = l
      else
         end = l-1
      end if
      return
      end
 
 
 

c-----------------------------------------------------------------------
c mds__leniftrim
c return length of string if trimmed of rightmost blanks
c
c 03/11/1999 ler
c-----------------------------------------------------------------------
      integer function mds__leniftrim(string)
      implicit none
 
      character*(*) string
      integer i, l
 
      l = len(string)
      mds__leniftrim = l
      do i=l,1,-1
         if (string(i:i).ne.' ') then
            goto 100
         else
            mds__leniftrim = mds__leniftrim - 1
         end if
      end do
 100  continue
 
      return
      end

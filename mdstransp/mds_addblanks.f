c-----------------------------------------------------------------------
c mds_addblanks
c replace null with trailing blanks in string
c
c 04/06/1999 ler
c-----------------------------------------------------------------------
      subroutine mds_addblanks(string)
      implicit none
 
      character*(*) string
      integer i, is, l
 
      l = len(string)
      is = index(string, char(0))
      if (is.ne.0) then
        do i=is,l
           string(i:i) = ' '
        end do
      end if
 
      end

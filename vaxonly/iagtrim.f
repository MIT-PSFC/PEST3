      integer function IAGTRIM(str)
c  return length of a character string minus any trailing blanks
c  if string is blank, returned length is 0
      implicit NONE
      character str*(*)
 
      integer jt,jl
 
      jt=0		!trimmed length
      jl=len(str)	!allocated length
      do while (jl.gt.jt)
        if(ichar(str(jl:jl)).gt.ichar(' ')) jt=jl
        jl=jl-1
      enddo
      IAGTRIM=jt
c
      return
      end

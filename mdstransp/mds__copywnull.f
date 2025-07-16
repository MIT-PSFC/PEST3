c-----------------------------------------------------------------------
c mds__copywnull
c copy source string into new string with trailing blanks removed and
c  null appended
c
c Input
c   source    -- source string to copy
c Output
c   newstring -- where to put string with null
c
c 1999-04-26 ler
c-----------------------------------------------------------------------
c Rationale:
c
c This is in support of avoiding: "string(1:l) // char(0)"
c Run-time concatenation is not supported in g77.
c
c Example of use:
c
c     integer mds_open(tree, shot)
c     character*(*) tree
c     integer shot
c     character*4097 newstring
c     call mds_copywnull(newstring, tree)
c     mds_open = MdsOpen(newstring, shot)
c     end
c
c Note: dynamic arrays not used to avoid tripping over potential
c       lack of consistency between fortran compilers
c-----------------------------------------------------------------------
      subroutine mds__copywnull(newstring, source)
      implicit none
 
      character*(*) newstring
      character*(*) source
 
      integer i, l, n
 
      newstring = source
 
      l = len(newstring)
      n = l
      do i=l,1,-1
         if (newstring(i:i).ne.' ') then
            goto 100
         else
            n = n - 1
         end if
      end do
 100  continue
 
      n = min(n,l-1)
      newstring(n+1:n+1) = char(0)
 
      end

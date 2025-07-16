      subroutine r8bsmoox(zarray,istride,inum,znorm)
c
c  extract array elements separated by istride into contiguous work array,
c  call smoother (bsmoo subroutine)
c  replace array elements with smoothed values back from work array
c
      use r8bsmoo_mod
      implicit NONE
c
      real*8 zwork(mj)
c
      integer istride,inum
      real*8 zarray(istride,inum)
      real*8 znorm(mj)
c
      integer i
c
c-----------------------------
c
      do i=1,inum
         zwork(i)=zarray(1,i)
      enddo
c
      call r8bsmoo(zwork,znorm)
c
      do i=1,inum
         zarray(1,i)=zwork(i)
      enddo
c
      return
      end

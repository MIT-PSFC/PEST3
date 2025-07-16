#include "fpreproc/library_names.h"
      subroutine filtr90(x,y1,y,n,epslon,delta,
     >   iddf,ddv,iblf,blv,
     >   xenda,xendb)
c
c  a simplified interface to filtr6 -- fortran-90
c
c  all of the arguments correspond to filtr6.for arguments (see the
c  extensive comments in filtr6.for).
c
c  the following simplifications are present:
c    1.  epslon is a scalar; epslon=0.0 => epslon feature not used;
c        there is no "eps2" feature.
c    2.  delta is a scalar
c    3.  the output quantities sm,ndrop,nblrs are dropped.
c
c...............
c
      implicit none
c
      integer n
      real x(n)                         ! indep. coordinate (input)
c
c  x must be strict ascending
c
      real y1(n)                        ! data to be smoothed (input)
      real y(n)                         ! the data as smoothed (output)
c
c  y and y1 need to be distinct arrays
c
      real epslon                       ! limit on abs(y(i)-y1(i))
c
c  epslon.le.0  =>  epslon = infinity
c
      real delta                        ! smoothing parameter
c
c  on output each y(j) is a triangular weighted average of the piecewise
c  linear interpolant of y1 over the range x(j)-delta to x(j)+delta.
c
      integer iddf                      ! data drop flag (=1 to enable)
      real ddv                          ! drip criterion
c
      integer iblf                      ! baseline flag (=1 to enable)
      real blv                          ! baseline value
c
c  **bdy conditions**
c
      real xenda                        ! LHS boundary condition
      real xendb                        ! RHS boundary condition
c
c  xend[a,b]=0.0 => endpoint is fixed, e.g. y(1)=y1(1) on output
c  xend[a,b]=1.0 => dy/dx --> 0 as x--> the boundary
c  intermediate values => intermediate action at end point.
c
c------------------------------
c
c  local
      real sm
      integer ndrop,nblrs
c
      real epslim,r1mach
      real, pointer :: epsarr(:), delarr(:)
c
      integer stat_alloc
c------------------------------
c
      allocate(epsarr(n),stat=stat_alloc)
      if(stat_alloc.ne.0) then
         call errmsg_exit(' ?filtr90:  epsarr array allocation failed!')
      endif
c
      allocate(delarr(n),stat=stat_alloc)
      if(stat_alloc.ne.0) then
         call errmsg_exit(' ?filtr90:  delarr array allocation failed!')
      endif
c
      epslim=r1mach(2)                  ! a large number
      if(epslon.le.0.0) then
         epsarr=epslim
      else
         epsarr=epslon
      endif
c
      delarr=delta
c
      call filtr6(x,y1,y,n,epsarr(1),n,epslim,0,delarr(1),iddf,ddv,
     >   iblf,blv,xenda,xendb,sm,ndrop,nblrs)
c
      deallocate(epsarr)
      deallocate(delarr)
c
      return
      end

      subroutine r8filtr90(x,y1,y,n,epslon,delta,
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
c
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer n
      REAL*8 x(n)                         ! indep. coordinate (input)
c
c  x must be strict ascending
c
      REAL*8 y1(n)                        ! data to be smoothed (input)
      REAL*8 y(n)                         ! the data as smoothed (output)
c
c  y and y1 need to be distinct arrays
c
      REAL*8 epslon                       ! limit on abs(y(i)-y1(i))
c
c  epslon.le.0  =>  epslon = infinity
c
      REAL*8 delta                        ! smoothing parameter
c
c  on output each y(j) is a triangular weighted average of the piecewise
c  linear interpolant of y1 over the range x(j)-delta to x(j)+delta.
c
      integer iddf                      ! data drop flag (=1 to enable)
      REAL*8 ddv                          ! drop criterion
c
      integer iblf                      ! baseline flag (=1 to enable)
      REAL*8 blv                          ! baseline value
c
c  **bdy conditions**
c
      REAL*8 xenda                        ! LHS boundary condition
      REAL*8 xendb                        ! RHS boundary condition
c
c  xend[a,b]=0.0 => endpoint is fixed, e.g. y(1)=y1(1) on output
c  xend[a,b]=1.0 => dy/dx --> 0 as x--> the boundary
c  intermediate values => intermediate action at end point.
c
c------------------------------
c
c  local
      REAL*8 sm
      integer ndrop,nblrs
c
      REAL*8 epslim
 
      real*8, parameter :: r8large = 1.0d36
 
      integer, parameter :: r8kind = kind(epslim)
      REAL(r8kind), pointer :: epsarr(:), delarr(:)
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
      epslim=r8large
      if(epslon.le.0.0E0_R8) then
         epsarr=epslim
      else
         epsarr=epslon
      endif
c
      delarr=delta
c
      call r8filtr6(x,y1,y,n,epsarr,n,epsarr,0,delarr,iddf,ddv,
     >   iblf,blv,xenda,xendb,sm,ndrop,nblrs)
c
      deallocate(epsarr)
      deallocate(delarr)
c
      return
      end

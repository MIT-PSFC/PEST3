#include "fpreproc/fp_const.h"
      subroutine zridderx (ivec, iok0, x1, x2, eps, eta, subr, x, ifail,
     >   ivecd,zinput,ninput,zoutput,noutput)
c
c  *** see also zriddery.for ***
c
c  **vectorized dmc Apr 2000**
c  real*8 root finder -- extended to allow extra arguments to
c  subr:  zinput(ninput) & zoutput(noutput) [dmc 4 Nov 1999]
c
      IMPLICIT NONE
c
      integer ivec                      ! vector size
      integer ifail                     ! error code
      REAL*8 x1(ivec), x2(ivec), eps, eta, x(ivec)
      logical iok0(ivec)                ! mask -- .true. to skip element
      external subr
c
      integer ivecd
      INTEGER ninput,noutput
      REAL*8 zinput(ivecd,ninput),zoutput(ivecd,noutput)
c
c     returns an approximate zero for the function (cf subr) given an
c     initial interval (x1,x2) which brackets a root.
c     The approximate zero is determined so that x is within the accuracy
c     of eps or |func(x)| < eta.
c     The algorithm is based on Ridders' method
c     see Numerical Recipes 2nd ed. p351
c
c  external subroutine header:
c     subroutine subr(ivec,iok,x,ansr,ivecd,zinput,ninput,zoutput,noutput)
c        integer ivec
c        logical iok(ivec)  ! if iok(i) TRUE then SKIP element i of vector
c        real*8 x(ivec)     ! input vector {x(i)}
c        real*8 ansr(ivec)  ! output vector {f(x(i))}
c        integer ivecd      ! auxilliary information vector size (.ge.ivec)
c        integer ninput,noutput
c        real*8 zinput(ivecd,ninput),zoutput(ivecd,noutput)
c
c ----------------------------------------------------------------------
c
      real*8, dimension(:), allocatable :: xl,xh,fl,fh,xm,fm,xnew,fnew,s
      logical, dimension(:), allocatable :: iok
c
      integer nc
c
      integer    itermax
      parameter (itermax = 100)
      REAL*8       rbig
      parameter (rbig    = 1.024D99)
      integer    stderr
      parameter (stderr  = 6)
      integer    iter,imsg,i,isave
c
      character  errmsg(3)*80
      data       errmsg(1)/'ZRIDDER: root not bracketed'/
      data       errmsg(2)/'ZRIDDER: anomalous behavior of func=0.(?)'/
      data       errmsg(3)/'ZRIDDER: exceeded maximum (100) iterations'/
      save       errmsg
c
c---------------------------------------------------------
      allocate(iok(ivec))
      allocate(xl(ivec),xh(ivec),fl(ivec),fh(ivec))
      allocate(xm(ivec),fm(ivec),xnew(ivec),fnew(ivec),s(ivec))
c
      nc=0
      do i=1,ivec
         iok(i)=iok0(i)
         if(.not.iok(i)) nc=nc+1
         xl(i) = x1(i)
         xh(i) = x2(i)
      enddo
      if(nc.eq.0) then
         ifail=0
         go to 10000                    ! all done
      endif
c
c  eval at endpts of intervals; check for acceptance of endpts...
c
      call subr(ivec,iok,xl,fl,ivecd,zinput,ninput,zoutput,noutput)
      nc=0
      do i=1,ivec
         if(.not.iok(i)) then
            if (ABS (fl(i)) .le. eta) then
               iok(i)=.TRUE.
               x(i)     = xl(i)
            else
               nc=nc+1
            endif
         endif
      enddo
      if(nc.eq.0) then
         ifail=0
         go to 10000                    ! all done
      endif
c
      call subr(ivec,iok,xh,fh,ivecd,zinput,ninput,zoutput,noutput)
      nc=0
      imsg=0
      do i=1,ivec
         if(.not.iok(i)) then
            if (ABS (fh(i)) .le. eta) then
               iok(i)=.TRUE.
               x(i)     = xh(i)
            else
               x(i) = rbig
               nc=nc+1
               if((fl(i).gt.0.D0).and.(fh(i).gt.0.D0).or.
     >            (fl(i).lt.0.D0).and.(fh(i).lt.0.D0))
     >            then
c
                  imsg=1
                  isave=i
c
               endif
            end if
         endif
      enddo
      if(imsg.eq.1) go to 9999
      if(nc.eq.0) then
         ifail=0
         go to 10000                    ! all done
      endif
c
      do iter=1,itermax
         do i=1,ivec
            if(.not.iok(i)) then
               xm(i)=(xl(i)+xh(i))/2.D0
            endif
         enddo
         call subr(ivec,iok,xm,fm,ivecd,zinput,ninput,zoutput,noutput)
         do i=1,ivec
            if(.not.iok(i)) then
               s(i)  = SQRT(max(0.D0,(fm(i)**2-fl(i)*fh(i))))
               if (s(i) .eq. 0.D0) then
                  isave=i
                  imsg = 2
                  x(i) = -rbig
               else
                  xnew(i)=xm(i)+(xm(i)-xl(i))*fm(i)/s(i) *
     >               SIGN (1.D0,fl(i)-fh(i))
               end if
            endif
         enddo
         call subr(ivec,iok,xnew,fnew,
     >      ivecd,zinput,ninput,zoutput,noutput)
         nc=0
         do i=1,ivec
            if(.not.iok(i)) then
               if (ABS (x(i)-xnew(i)) .le. eps .or.
     >            ABS (fnew(i)) .le. eta) then
                  iok(i)=.TRUE.
               endif
               x(i) = xnew(i)
            end if
            if(.not.iok(i)) then
               if (SIGN (fm(i),fnew(i)) .ne. fm(i)) then
                  xl(i) = xm(i)
                  fl(i) = fm(i)
                  xh(i) = xnew(i)
                  fh(i) = fnew(i)
               else if (SIGN (fl(i),fnew(i)) .ne. fl(i)) then
                  xh(i) = xnew(i)
                  fh(i) = fnew(i)
               else if (SIGN (fh(i),fnew(i)) .ne. fh(i)) then
                  xl(i) = xnew(i)
                  fl(i) = fnew(i)
               else
                  call errmsg_exit(
     >               'subroutine ZRIDDERX: programming error')
               end if
               if (ABS (xl(i)-xh(i)) .le. eps) then
                  iok(i)=.TRUE.
               else
                  nc=nc+1
               endif
            endif
         enddo
         if(nc.eq.0) then
            ifail=0
            go to 10000
         endif
      end do
c
c  get here if for some point, convergence failed
c
      imsg = 3
      do i=1,ivec
         if(.not.iok(i)) isave=i
      enddo
c
 9999 continue
      if((imsg.gt.0).and.(imsg.le.3)) then
         write (stderr, '(a)')  errmsg(imsg)
         write (stderr, '(a)') 'i xl(i), fl(i), xh(i), fh(i):'
         write (stderr,   *  )  isave, xl(isave), fl(isave),
     >      xh(isave), fh(isave)
      else
         call errmsg_exit('subroutine ZRIDDERX: unspecified error')
      endif
c
      ifail = imsg
c
10000 continue
      deallocate(xl,xh,fl,fh,xm,fm,xnew,fnew,s,iok)
      return
c
      end

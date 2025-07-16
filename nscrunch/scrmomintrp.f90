      SUBROUTINE scrmomintrp(rmcx,zmcx,xsize,imomuse,ic0,icc,idim1,immax, &
     &   inpts,ier)
!
      USE mintrp
      IMPLICIT NONE
       INTEGER 	idim1 ! <momintrp.f>
       INTEGER 	itheta ! <momintrp.f>
       REAL*8 	ztest ! <momintrp.f>
       REAL*8 	rtest ! <momintrp.f>
       REAL*8 	scrYceptcon ! <momintrp.f>
       REAL*8 	y ! <momintrp.f>
       REAL*8 	x ! <momintrp.f>
       REAL*8 	zstat ! <momintrp.f>
       REAL*8 	xguessp ! <momintrp.f>
       REAL*8 	zftol ! <momintrp.f>
       REAL*8 	zxtol ! <momintrp.f>
       REAL*8 	delx ! <momintrp.f>
       REAL*8 	xguess ! <momintrp.f>
       REAL*8 	zxmax ! <momintrp.f>
       REAL*8   zero
       INTEGER 	ic2 ! <momintrp.f>
       INTEGER 	ics ! <momintrp.f>
       INTEGER 	ic0 ! <momintrp.f>  ! index to axis
       INTEGER 	im ! <momintrp.f>
       REAL*8   xsize(idim1)
       INTEGER 	imomuse(idim1) ! <momintrp.f>
       REAL*8 	ztheta ! <momintrp.f>
       INTEGER 	ith ! <momintrp.f>
       INTEGER 	ic ! <momintrp.f>
       INTEGER 	i ! <momintrp.f>
       INTEGER 	ier ! <momintrp.f>
       REAL*8 	zpi ! <momintrp.f>
       INTEGER 	ichkn ! <momintrp.f>
       INTEGER 	inpts ! <momintrp.f>
       INTEGER 	immax ! <momintrp.f>
       integer icskip,icx,ica
 
!
!  dmc 26 June 1997
!  use a low order rotated polynomial interpolation to fill theta
!  lines btw scrunched surfaces; use fft to extract moments at
!  interpolated (non-scrunched) surfaces, truncating to use only
!  imomuse(j) moments at the j'th surface
!
      REAL*8 rmcx(idim1,0:immax,2),zmcx(idim1,0:immax,2)
      integer icc(6)

      REAL*8 rmcx1(0:immax,2),zmcx1(0:immax,2),xsw,wsw,xsrch0
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rwork, zwork
!
! !
!  ncmax -- max no. of surfaces (array dimension) in MINTRP
!  immax -- max no. of moments  (array dimension) in MINTRP
!
!  imomuse(...)  -- actual no. of moments to use
!  icc (moption.eq.2)
!           icc(1) gives index of inmost scrunched surface
!           icc(4) gives index of outermost scrunched surface
!           icc(2) gives a scrunched surface at r/a ~ 0.5
!           icc(3) gives index of surface near outside for numerical
!                  diffentiation of theta line direction...
!  icc (moption.eq.3)
!           icc(1:3) gives ordered index of near-axis surfaces
!           icc(6) gives index of outermost scrunched surface
!           icc(4) gives a scrunched surface at r/a ~ 0.5
!           icc(5) gives index of surface near outside for numerical
!                  diffentiation of theta line direction...
!
! input and modified:
!  rmcx(...) and zmcx(...) -- on input, moments at the scrunched surfaces
!           only; on output, moments sets for all surfaces
!
! output:
!  ier -- completion code.  0=success.
!
! workspace:
!  rwork(...) and zwork(...) (R,Y) expansions workspaces
!  inpts -- 1st dimension of workspaces
!
!  an fft is used with rows of rwork and zwork to get the intermediate
!  moments sets.
!
!---------------------------------------
!  external function f(x), solve f(x)=0 for x to find the
!  theta-line | flux-surface intercept
!
      external scrXceptcon
!
!---------------------------------------
!  local:
!
!  for moments expansion evaluation:  sincos ladder
      REAL*8 zsn(100),zcs(100)
!
!  fft output arrays -- dimension should match largest value of ichk!
!
      REAL*8 zs1d(128),zc1d(128)
!
!  fft workspace arrays -- largest value of ichk, +1
!
      REAL*8 zsq,zdel,ztol
 
      parameter(ichkn=4)
      integer ichk(ichkn)
      integer ilim,imoms
!
      integer lunmsg
!
      logical lyes
!------
!
      data zpi/3.1415926535897931D+00/
!
      data ichk/16,32,64,128/
!
!-------------------------------------------------------------------
!
      call scrgetlun(lunmsg)
!
! local stack arrays
      ALLOCATE(rwork(inpts+1,idim1),zwork(inpts+1,idim1)); rwork=0; zwork=0
!
!
!  0.  error check
!
      ier=0
      if(moption.eq.3) then
         if(icc(1).ge.icc(2)) ier=1
         if(icc(2).ge.icc(3)) ier=1
         if(icc(3).ge.icc(4)) ier=1
         if(icc(4).ge.icc(5)) ier=1
         if(icc(5).ge.icc(6)) ier=1
         if(icc(6).gt.ncmax) ier=2
         ilim=0
         icskip=5
      else if(moption.eq.2) then
         if(icc(1).ge.icc(2)) ier=1
         if(icc(2).ge.icc(3)) ier=1
         if(icc(3).ge.icc(4)) ier=1
         if(icc(4).gt.ncmax) ier=2
         ilim=0
         icskip=3
      endif
!
      if(ier.eq.1) then
         write(lunmsg,'('' ?momintrp -- icc argument order error'')')
      else if(ier.eq.2) then
         write(lunmsg,'('' ?momintrp -- icc argument limit error'')')
      endif
      if(ier.gt.0) go to 1000
!
      ier=3
      do i=1,ichkn
         if(inpts.eq.ichk(i)) ier=0
      enddo
!
      if(ier.eq.3) then
         write(lunmsg, &
     &'('' ?momintrp -- 1st work array dimension not of form 2**N'')')
         write(lunmsg, &
     &'('' ?            or too small or too large.'')')
         go to 1000
      endif
!
!  1.  expand the already scrunched surfaces
!
      do i=1,6
         ic=icc(i)
         if(ic.gt.0) then
            do ith=1,inpts
               ztheta=(ith-1)*2.0D0*zpi/inpts
               imoms=max(ilim,imomuse(ic))
               CALL scrsincos(ztheta,imoms,zsn,zcs)
               rwork(ith,ic)=rmcx(ic,0,1)
               zwork(ith,ic)=zmcx(ic,0,1)
               do im=1,imoms
                  rwork(ith,ic)=rwork(ith,ic)+rmcx(ic,im,1)*zcs(im)+ &
     &               rmcx(ic,im,2)*zsn(im)
                  zwork(ith,ic)=zwork(ith,ic)+zmcx(ic,im,1)*zcs(im)+ &
     &               zmcx(ic,im,2)*zsn(im)
               enddo
7701           format(' ith=',i3,' ic=',i3,' R=',1pe11.4,' Y=',1pe11.4)
            enddo
            rwork(inpts+1,ic)=rwork(1,ic)
            zwork(inpts+1,ic)=zwork(1,ic)
         endif
      enddo
!
!  2.  interpolate the theta lines
!
      if(moption.eq.2) then
         ic2=icc(4)
      else if(moption.eq.3) then
         ic2=icc(6)
      endif
!
!  loop over theta lines
!
      do ith=1,inpts
!
!  transform to rotated frame
!  2x2 rotations
!   T = ( zcosa -zsina )      Tinv = ( zcosa  zsina )
!       ( zsina  zcosa )             (-zsina  zcosa )
!
!  ( x ) = Tinv * ( R - R0 )   ( R ) = ( R0 ) + T * ( x )
!  ( y )          ( Z - Z0 )   ( Z )   ( Z0 )       ( y )
!
!  origin is 1st surface point
!
         CALL scrfsetup(ith,icc,rwork,zwork,inpts,idim1,zxmax,ier)
         if(ier.ne.0) then
            write(lunmsg,*) ' ?scrmomintrp: theta line fit failure, ith=',ith
            go to 1000
         endif
!
!  find theta-line | flux surface intercepts
!  first initial guess
         xguess=xmax
         delx=0.25D0*xmax
!
         zxtol=1.0D-10*xmax               ! x error tolerance
         zftol=1.0D-10*r0                 ! f(x) an error tolerance
!
         do ic=ic2,ic0+1,-1
!
!  use numeric intercept finder on each surface
!  scrunched surfaces are also probed, so that treatment is consistent on
!  all surfaces
!    see scrXceptcon(x,zdum) function
!
            xguessp=xguess
            xguess=xguess-delx
!
            ksurf=ic
            zero = 0.0D0
            CALL scrmfalsi(scrXceptcon,xguess,zftol,zxtol, &
                 &            xspl(1),1.05D0*xmax,100,lunmsg,0,zstat)
            if(zstat.ne.0.0D0) then
               write(lunmsg,*) &
                    ' ?momintrp:  MFALSI error (check input data).'
               ier=99
               exit
            endif
            delx=xguessp-xguess
!  intercept in (x,y) space
            x=xguess
            y=scrYceptcon(x,Rtest,Ztest)
!  transform back to (R,Y) and save
            rwork(ith,ic) = Rtest
            zwork(ith,ic) = Ztest
!
            if(ier.ne.0) exit
         enddo                          ! loop over surfaces, out -> in
!
         if(ier.ne.0) exit
      enddo                             ! loop over theta angles

!  close all contours
      rwork(inpts+1,1:idim1)=rwork(1,1:idim1)
      zwork(inpts+1,1:idim1)=zwork(1,1:idim1)
!
      if(ier.ne.0) go to 1000
!
!  OK now information is available to do an fft to get moments back
!
!  3.  compute the unscrunched moments set from the interpolation results
!
      itheta=inpts
      do ic=ic0+1,ic2
!
!  do this on all surfaces, for consistency...
!
         imoms=max(ilim,imomuse(ic))
!  R moments
         call r8fftsc(rwork(1,ic),itheta,zs1d,zc1d,ier)
         if(ier.ne.0) go to 1000
         rmcx(ic,0,1)=zc1d(1)/(2.0D0*itheta)
         do im=1,imoms
            rmcx(ic,im,1)=zc1d(im+1)/(itheta)
            rmcx(ic,im,2)=zs1d(im+1)/(itheta)
         enddo
!  Z moments
         call r8fftsc(zwork(1,ic),itheta,zs1d,zc1d,ier)
         if(ier.ne.0) go to 1000
         zmcx(ic,0,1)=zc1d(1)/(2.0D0*itheta)
         do im=1,imoms
            zmcx(ic,im,1)=zc1d(im+1)/(itheta)
            zmcx(ic,im,2)=zs1d(im+1)/(itheta)
         enddo
      enddo
!
! moption.eq.3 -- cleanup higher order moments' approach to axis
!
      if(moption.eq.3) then
         !  2nd moment:

         call a2smoo(rmcx(1:idim1,2,1))
         call a2smoo(rmcx(1:idim1,2,2))
         call a2smoo(zmcx(1:idim1,2,1))
         call a2smoo(zmcx(1:idim1,2,2))

         !  moments 3 and up
         ztol=rmcx(ic0,0,1)*1.0d-6
         do im=3,imomuse(ic2)
            xsrch0=0.0d0
            if(im.gt.ilim) then
               do i=ic2-1,ic0+1,-1
                  if(imomuse(i).lt.im) then
                     xsrch0=xsize(i)
                     exit
                  endif
               enddo
            endif
            call axsmoo(rmcx(1:idim1,im,1),im)
            call axsmoo(rmcx(1:idim1,im,2),im)
            call axsmoo(zmcx(1:idim1,im,1),im)
            call axsmoo(zmcx(1:idim1,im,2),im)
         enddo
      endif
!
 1000 continue
!
!---------------------------------------------------------------
!
      DEALLOCATE(rwork, zwork)
!
      return
 
      contains
        subroutine a2smoo(zmoms)
          REAL*8, dimension(:) :: zmoms  ! moments on these surfaces

          integer :: icbrk,ic

          real*8 :: x2,x2sum,x2sum2,ssum,scav,wav

          do ic=ic0+1,ic2
             if(imomuse(ic).gt.2) then
                icbrk=ic
                exit
             endif
          enddo

          x2sum=0
          ssum=0
          do ic=ic0+1,icbrk
             x2=xsize(ic)*xsize(ic)
             ssum=ssum+zmoms(ic)   ! /x2 *x2
             x2sum=x2sum+x2
          enddo

          scav=ssum/x2sum  ! avg scaled moment

          x2sum2=0
          do ic=ic0+1,icbrk-1
             x2=xsize(ic)*xsize(ic)
             x2sum2=x2sum2+x2
             wav = 1 - x2sum2/x2sum
             zmoms(ic)=x2*scav*wav + zmoms(ic)*(1-wav)
          enddo

        end subroutine a2smoo

        subroutine axsmoo(zmoms,im)
 
          REAL*8, dimension(:) :: zmoms  ! moments on these surfaces
          integer im          ! order of moments
 
          ! local...
 
          real*8 :: zdelx,xbrk1,xbrk2
          integer i
 
          !-----------------------
          zdelx=0.15d0

          do i=ic0+1,ic2
             if(abs(zmoms(i)).gt.ztol) then
                if(xsrch0.gt.0.0d0) then
                   xbrk2=max(xsrch0+zdelx,xsize(i-1))
                else
                   xbrk2=xsize(i-1)
                endif
                xbrk1=xbrk2-zdelx
                exit
             endif
          enddo

          if(xbrk1.lt.0.d0) then
             if(im.lt.3) then
                return  ! no change
             else
                xbrk1=0.0d0
                xbrk2=zdelx
             endif
          endif

          do i=ic0+1,ic2
             if(xsize(i).le.xbrk1) then
                zmoms(i)=0.0d0
             else if(xsize(i).lt.xbrk2) then
                xsw=(xbrk2-xsize(i))/(xbrk2-xbrk1)
                wsw=xsw*xsw*xsw*(xsw*(15.0d0-6.0d0*xsw)-10.0d0) + 1.0d0
                zmoms(i)=zmoms(i)*wsw
             endif
          enddo
 
        end subroutine axsmoo
      end


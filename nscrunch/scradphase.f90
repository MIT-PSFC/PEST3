      SUBROUTINE scradphase(raxis,yaxis,rmcx1,ymcx1,immax,imom,zthadj)
!
!  adjust phase of theta parametrization on "1" surface
!    want theta=0 to give R.gt.raxis, y.eq.yaxis
!
!  input:  raxis,yaxis -- mag. axis position
!          rmcx1,ymcx1 -- moments of surface (unadjusted)
!          immax -- max no. of moments (array dimension)
!          imom  -- no. of active moments (use 0:imom).
!
! output:  rmcx1,ymcx1 -- adjusted
!
!  the idea here is to align the phases of two surfaces, to reduce as
!  much as possible, on average, the non-orthogonality of the theta
!  lines -- i.e. remove any "swastika-like" skewing of the theta
!  parametrization.
!
      IMPLICIT NONE
       REAL*8 	delth ! <adphase.f>
       REAL*8 	angspan ! <adphase.f>
       REAL*8 	angnorm ! <adphase.f>
       REAL*8 	zthadj ! <adphase.f>
       INTEGER 	imom ! <adphase.f>
       INTEGER 	im ! <adphase.f>
       REAL*8 	ztheta ! <adphase.f>
       INTEGER 	ith,itp ! <adphase.f>
       REAL*8 	zpi ! <adphase.f>
       INTEGER 	inpts ! <adphase.f>
       INTEGER 	immax ! <adphase.f>
!
       integer iy0   ! interval where y=yaxis falls, with r.gt.raxis
!
       real*8 raxis,yaxis
      REAL*8 rmcx1(0:immax,2)
      REAL*8 ymcx1(0:immax,2)
!
      parameter (inpts=128)
      REAL*8 r0(inpts),y0(inpts),r1(inpts),y1(inpts)
      REAL*8 enorm(2),espan(2)
!
      real*8 zthetap,zthp,zth,zyp,zy,zero
!
      REAL*8 zsn(100),zcs(100)
!
      integer lunmsg
!
      data zpi/3.1415926535897931D+00/
!
!---------------------------------------------------------------------
!
      call scrgetlun(lunmsg)
!
!  expand surfaces into (R,Y) sequences
!
      zero=0
      ztheta=zero
      iy0=0
!
      do ith=1,inpts
         ztheta=(ith-1)*2.0d0*zpi/inpts
         CALL scrsincos(ztheta,imom,zsn,zcs)
         r1(ith)=rmcx1(0,1)
         y1(ith)=ymcx1(0,1)
         do im=1,imom
            r1(ith)=r1(ith)+rmcx1(im,1)*zcs(im)+rmcx1(im,2)*zsn(im)
            y1(ith)=y1(ith)+ymcx1(im,1)*zcs(im)+ymcx1(im,2)*zsn(im)
         enddo
      enddo
!
!  look for interval on large major radius side containing y=yaxis
!
      do ith=1,inpts
         zthetap=ztheta
         ztheta=(ith-1)*2.0d0*zpi/inpts
         if(ith.gt.1) then
            itp=ith-1
            if((r1(ith).gt.raxis).and.(r1(itp).gt.raxis)) then
               if((y1(ith)-yaxis)*(y1(itp)-yaxis).le.zero) then
                  if(iy0.gt.0) then
                     if(y1(itp).ne.yaxis) call errmsg_exit( &
                          & '?scradphase:  mult. y=yaxis intervals?')
                  else
                     iy0=itp
                     zthp=zthetap
                     zth=ztheta
                     zyp=y1(itp)
                     zy=y1(ith)
                  endif
               endif
            endif
         endif
!
      enddo
!
      if(iy0.eq.0) then
         zy=y1(1)
         zyp=y1(inpts)
         zth=zero
         zthp=ztheta-2.0d0*zpi
      endif
!
      if((zy-yaxis)*(zyp-yaxis).gt.zero) then
         call errmsg_exit('?scradphase: y=yaxis interval search failure?')
      endif
!
      zthadj=zthp+(yaxis-zyp)*(zth-zthp)/(zy-zyp)  ! theta -> y=yaxis
!
      CALL scrphasead(rmcx1,ymcx1,immax,imom,zthadj)
!
      return
      end

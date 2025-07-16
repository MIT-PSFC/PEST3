      subroutine scriscrunch(rin1d, zin1d, jmax, ns, &
     &                    rmc1d, rms1d, ymc1d, yms1d, &
     &                    ic1,inc,imom,ier)
!
!  see scriscrunchx (below) for description of arguments.
!  estimate x and call scriscrunchx
!
      IMPLICIT NONE
      REAL*8, dimension(*) :: rin1d,zin1d
      REAL*8, dimension(*) :: rmc1d, rms1d, ymc1d, yms1d
      integer jmax,ns
      integer ic1,inc,imom,ier
!
!-------------------------------------------
      integer i,k
      REAL*8 rsize(ns),xsize(ns)
!
!-------------------------------------------
      rsize=0.0D0
      xsize=0.0D0
!
      do i=inc,ic1,-1
         k=(i-1)*jmax
         rsize(i)=maxval(rin1d(k+1:k+jmax))-minval(rin1d(k+1:k+jmax))
         xsize(i)=rsize(i)/rsize(inc)
      enddo
!
      call scriscrunchx(xsize, rin1d, zin1d, jmax, ns, &
     &                    rmc1d, rms1d, ymc1d, yms1d, &
     &                    ic1,inc,imom,ier)
!
      return
      end
!-----------------------------------------------------------------------
      SUBROUTINE scriscrunchx(xrho,rin1d, zin1d, jmax, ns, &
     &                    rmc1d, rms1d, ymc1d, yms1d, &
     &                    ic1,inc,imom,ier)
!
      USE mintrp
      USE scrunch_inc1
!
!             input:   xrho(1:inc) ... estimated or actual normalized rho
!                      (recommended based on toroidal flux sqrt(phi/philim))
!             input:   rin1d, zin1d -- jmax*ns grid vectors. Indexing is
!                      such that the i-th poloidal point on the j-th radial
!                      surface, (rin1d(k), zin1d(k)), is given by
!                      k=(j-1)*jmax+i. Note that the contours must close,
!                      i.e. rin1d(k1)=rin1d(k2) where k1=(j-1)*jmax+1 and
!                      k2 = k1 + jmax - 1.
!
!             input:   jmax     -- number of poloidal points
!
!             input:   ns       -- number of surfaces
!
!             output:  rmc1d, rms1d, ymc1d, yms1d -- COS and SIN
!                      moments stored in 1-d arrays. Indexing is such
!                      that, eg, rmc1d( (i-1)*(imom+1)+ m+1 ) represents the
!                      m-th mode at the i-th surface.
!
!              input:  ic1:  first surface for which moments are needed
!                      this should be the first **non-singular** surface
!                      not the magnetic axis
!                      e.g. LCENTR+1 in TRCOM
!
!                      note:  the magnetic axis position is estimated
!                             if ic1=1 & no axis position is given.
!
!              input:  inc:  last surface for which moments are needed
!                      this should be the boundary surface.  Generally
!                      inc = ns.
!
!              input:  imom:  number of moments to use
!                      moments 0:imom are used.
!
!  scrunch the inner and outermost surfaces, fit a smooth theta
!  contour between them, and from this build a complete moments set.
!
!  ** caution **
!  there should be a one to one correspondence in the indexing of the
!  RIN and ZIN arrays in with the indexing of the rmcx,ymcx
!  output arrays; i.e. ic1 should point to the first non-singular
!  surface in both array pairs!
!
      IMPLICIT NONE
       REAL*8, DIMENSION(*) :: rin1d, zin1d,xrho
       REAL*8, DIMENSION(*) :: rmc1d, rms1d, ymc1d, yms1d
!
!
!
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE ::  rmcx, ymcx
!
!
!
      LOGICAL ioflag
      LOGICAL iasym
      LOGICAL izak
      data ioflag/.false./               ! descur chat flag
      data iasym/.true./                ! asymmetry flag
      data izak/.false./
 
!
!-----------------------------------------------------------------
!
      INTEGER ns
      INTEGER idum
      INTEGER ier
      INTEGER iok
      INTEGER imom,imomuse(ns),jmomuse(0:imom)
      INTEGER ic0
      INTEGER inum
      INTEGER itm
      INTEGER ith
      INTEGER ict
      INTEGER inth
      INTEGER inc
      INTEGER ic
      INTEGER ic1
      INTEGER i
      INTEGER im
      INTEGER is
      INTEGER inpts
      INTEGER immax
      INTEGER k
      INTEGER jmax
      INTEGER icc(6),icdmin
      REAL*8 zdel
      REAL*8 fsq,xtrim
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rmcx1 ! (0:imom,2)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ymcx1 ! (0:imom,2)

!................................................................
!...begin heap allocation

      ALLOCATE(rmcx1(0:imom,2))
      ALLOCATE(ymcx1(0:imom,2))
      ALLOCATE(rmcx(ns,0:imom,2), ymcx(ns,0:imom,2))
!...end heap allocation
!...begin set dynamic arrays to zero
      rmcx1 = 0.0D0
      ymcx1 = 0.0D0
      rmcx  = 0.0D0
      ymcx  = 0.0D0
!...end set dynamic arrays to zero

!-----------------------------------------------------------------
!
      ier = 0

      CALL scrpreset(2*jmax, ns, imom)
 
      if (facmid.le.0.0D0) facmid=0.551D0
      if (facedg.le.0.0D0) facedg=0.851D0
      if (facdx0.le.0.0D0) facdx0=0.100D0
 
      CALL scrmemory
 
      CALL scrclear
      DO is=1, ns
         k = (is-1)*jmax
         rin(1:jmax,is) = rin1d(k+1:k+jmax)
         zin(1:jmax,is) = zin1d(k+1:k+jmax)
      END DO
!
!  arrange for truncation of no. of moments -- zero out the small ones.
!
      xtrim = 1.0D-4
!
      imomuse=imom
      jmomuse=1
      do is=1,inc
         if(is.lt.ic1) then
            imomuse(is)=0
            cycle
         else if(is.eq.ic1) then
            imomuse(is)=2
            cycle
         else
            imomuse(is)=imomuse(is-1)
         endif
         if(xrho(is).le.0.0D0) cycle
         imomuse(is)=max(1,imomuse(is))
         if(imomuse(is).eq.imom) cycle
         do
            if(xrho(is)**(imomuse(is)+1).le.xtrim) exit
            imomuse(is)=imomuse(is)+1
            jmomuse(imomuse(is))=is
            if(imomuse(is).eq.imom) exit
         enddo
      enddo
 
      if(ns.gt.ncmax) then
         write(lunmsg,*) '  ns=',ns,'; ncmax=',ncmax
         CALL scrabwarn(lunmsg,'?iscrunch:  ns exceeds ncmax parameter')
         ier=99
         go to 199
      endif
!
!  get mag. axis from nearest non-singular surfaces
      if(ic1.le.1) then
         raxis = (4*SUM(rin(1:jmax,ic1))-SUM(rin(1:jmax,ic1+1)))/(3*jmax)
         zaxis = (4*SUM(zin(1:jmax,ic1))-SUM(zin(1:jmax,ic1+1)))/(3*jmax)
      else
         raxis = SUM(rin(1:jmax,ic1-1))/jmax
         zaxis = SUM(zin(1:jmax,ic1-1))/jmax
      endif
      if(ic1.gt.1) then
         do k=1,ic1-1
            rmcx(k,0,1)=raxis
            ymcx(k,0,1)=zaxis
         enddo
      endif
!
!  dmc 23 July 1997 -- compress duplicate points out of RIN,ZIN arrays
!
      do ic=ic1,inc
         inth=jmax
         ict=1
         do ith=2,inth
            itm=ith-1
            if((rin(itm,ic).eq.rin(ith,ic)).and. &
     &         (zin(itm,ic).eq.zin(ith,ic)) ) then
               write(lunmsg,6001) ic,itm,ith,rin(itm,ic),zin(itm,ic)
 6001          format(/ &
     &' ?iscrunch -- duplicate contour point, surface #',i3, &
     &' theta indices: ',2(1x,i4)/ &
     &' ...R=',1pe12.5,' Z=',1pe12.5)
               ier=1
               go to 199
            else
               ict=ict+1
               if(ict.lt.ith) then
                  rin(ict,ic)=rin(ith,ic)
                  zin(ict,ic)=zin(ith,ic)
               endif
            endif
         enddo
         ntheta0=ict
         if( abs(rin(1,ic) - rin(ict,ic)) >  1.D-8) then
            write(lunmsg,*) '?iscrunch:  end point does not match starting point.'
            write(lunmsg,*) ' surface ic=',ic,'  nth=',ict
            write(lunmsg,*) ' rin(1,ic) = ',rin(1,ic)
            write(lunmsg,*) ' rin(ict,ic) = ',rin(ict,ic)
            ier = 2
            go to 199
         endif
         if( abs(zin(1,ic) - zin(ict,ic)) > 1.D-8) then
            write(lunmsg,*) '?iscrunch:  end point does not match starting point.'
            write(lunmsg,*) ' surface ic=',ic,'  nth=',ict
            write(lunmsg,*) ' zin(1,ic) = ',zin(1,ic)
            write(lunmsg,*) ' zin(ict,ic) = ',zin(ict,ic)
            ier = 2
            go to 199
         endif
      enddo
!
!  scrunch innermost and outermost surface
      inum=inc-ic1+1
      ic0=ic1-1
!
      icc=0
      icdmin=99
      if(moption.eq.3) then
         !
         !  select surfaces near axis; select that these are
         !  distinct and bound away from axis; if there are too few
         !  radial points they won't be & we revert to a simpler method.
         !
         do i=1,3
            icc(i)=ic0+i*inum*facdx0      ! near axis
            if(i.eq.1) then
               icdmin=min(icdmin,(icc(1)-ic0))
            else
               icdmin=min(icdmin,(icc(i)-icc(i-1)))
            endif
         enddo
         if(icdmin.eq.0) then
            write(lunmsg,*) '%iscrunch: only ', &
                 inum,' surfaces: reverting to moption=2 radial fit.'
            moption=2
         endif
      endif
      if(moption.eq.2) then
         icc(1)=ic1
         icc(2)=ic0+inum*facmid         ! one in the middle
         icc(3)=ic0+inum*facedg         ! & one near the edge, for numeric diff
         icc(4)=inc
      else if(moption.eq.3) then
         !
         !  moption.eq.3 OK -- get rest of surfaces.
         !
         icc(4)=ic0+inum*facmid         ! one in the middle
         icc(5)=ic0+inum*facedg         ! & one near the edge, for numeric diff
         icc(6)=inc
      else
         CALL scrabwarn(lunmsg,'?iscrunch:  illegal moption value.')
         ier=99
         go to 199
      endif
!
      ic=ic1-1
      do 150 i=1,6
!
         if(icc(i).eq.0) go to 150
!
!  for moments copying mechanism
!
         ic=icc(i)
!
         inth=ntheta0
!
!  do the scrunch
!
         inum=min(17,(imomuse(ic)+1))
         iok=1
         idum=0   ! set to (fortran LUN>0) for scrdescur to save a debug file
 
         CALL scrdesreg(rin(1,ic),zin(1,ic),inth,inum,idum,ioflag,fsq, &
              &      iasym,izak,imom,rmcx1,ymcx1,iok)
!
         if(iok.ne.1) then
            CALL scrabwarn(lunmsg,'?iscrunch:  error code on exit from descur')
            ier=max(1,iok)
            go to 199
         endif
!
!  remove any phase twisting -- set theta=0 -->z=zaxis, large radius side
!  uncomment this to add a twist, for testing...
!debug            call phasead(rmcx1,ymcx1,imom,imom,0.1)
         ict=0
105      continue
         ict=ict+1
         CALL scradphase(raxis,zaxis,rmcx1,ymcx1,imom,imom,zdel)
         if((abs(zdel).gt.1.0D-10*raxis).or.(ict.eq.1)) go to 105
!
!  copy the moments
!
         do is=1,2
            do im=0,imom
               rmcx(ic,im,is)=rmcx1(im,is)
               ymcx(ic,im,is)=ymcx1(im,is)
            enddo
         enddo
!
 150  continue
!
!  set up moments for unscrunched surfaces, by interpolation scheme
!
      inpts = 64
      if(imom.gt.32) inpts = 128

      CALL scrmomintrp(rmcx,ymcx,xrho,imomuse,ic0,icc,ns,imom, &
     &      inpts,iok)
!
      DO is=1, ns
         DO im=0, imom
            k = (is-1)*(imom+1) + im + 1
            rmc1d(k) = rmcx(is, im, 1)
            rms1d(k) = rmcx(is, im, 2)
            ymc1d(k) = ymcx(is, im, 1)
            yms1d(k) = ymcx(is, im, 2)
         END DO
      END DO
!
!-------------------------------------------------------------------
!
      ier = iok              ! ier = 0 means ok
 
199   continue
 
!...begin restore memory
      IF( ALLOCATED(rmcx1) ) DEALLOCATE(rmcx1)
      IF( ALLOCATED(ymcx1) ) DEALLOCATE(ymcx1)
      IF( ALLOCATED(rmcx) ) DEALLOCATE(rmcx)
      IF( ALLOCATED(ymcx) ) DEALLOCATE(ymcx)
!...end restore memory
      CALL scrfree
      return
      end
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"

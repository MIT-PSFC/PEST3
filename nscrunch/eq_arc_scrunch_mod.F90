module eq_arc_scrunch_mod

  ! "private" module for nscrunch library

  implicit NONE

  integer :: ntha
  real*8, dimension(:,:), allocatable :: tha_pkg   ! working TH grid

  real*8, dimension(:,:), allocatable :: Rspl  ! R(TH)
  real*8, dimension(:,:), allocatable :: Zspl  ! Z(TH)
  real*8, dimension(:,:), allocatable :: Lspl  ! L(TH) precise arc-length
  real*8, dimension(:), allocatable :: wk,thmap

  real*8 :: thzero,Lzero,Ltotal

  real*8, parameter :: C2PI = 6.2831853071795862D+00
  real*8, parameter :: CSMALL = 1.0D-14
  real*8, parameter :: ZERO = 0.0D0

  contains

    subroutine eq_arc_scrunch1(Raxis,Zaxis, R,Z, rmc,rms,zmc,zms, imom,ier)

      !--------------------------------------
      !  equal arc scrunch a single surface; define theta=0 as being on
      !  the outer half midplane, drawn horizontally outward from the 
      !  magnetic axis; store R[0:2pi],Z[0:2pi] w.r.t. an equal arc length
      !  defined poloidal angle parameter; also compute FFT moments 
      !  representation

      real*8, intent(in) :: Raxis,Zaxis        ! axis location

      real*8, dimension(:), intent(inout) :: R,Z  ! a surface contour
      !  each point must be distinct; 1st and last points must match.
      !  these are assumed to be smoothly spaced, e.g. equispaced in
      !  angle defined by an arc-tangent from the magnetic axis

      integer, intent(in) :: imom  ! #moments not counting 0th moment

      !  the R & Z cos and sin moments:

      real*8, intent(out) :: rmc(0:imom),rms(0:imom),zmc(0:imom),zms(0:imom)

      integer, intent(out) :: ier

      !--------------------------------------
      !  external routines for root finder...

      external eq_arc_fmidplane
      external eq_arc_findarc

      !--------------------------------------
      !  local variables...

      integer :: i,ii,im,idum,izero,ifound,iertmp,inm1
      real*8 :: zdR,zdZ,zdL,zdth,eps,eta,ztest,zdlmin,zdlmax

      integer, dimension(3) :: ict_eval = (/ 1, 0, 0 /)
      integer, dimension(3) :: ict_der1 = (/ 0, 1, 0 /)

      real*8 :: w0(10),th0(10),wvec(10),thvec(10)  ! for integration of L(TH)
      real*8 :: drdth(10),dzdth(10),dldth,thzero,zlsum,zdum

      logical, dimension(:), allocatable :: imask
      real*8, dimension(:), allocatable :: ltarg,thfnd,th1,th2,zdumma

      !  fft arrays -- dimension should match largest value of ichk!

      integer :: nfft
      real*8, dimension(:), allocatable :: thfft,rfft,zfft,zs1d,zc1d

      !--------------------------------------
      !  debugging...
      !
#ifdef __DEBUG

      logical, save :: idebug = .FALSE.
      integer, save :: ictr = 0

      integer, parameter :: neval = 1001
      real*8 theval(neval),dldth1(neval),dldth2(neval),zdra(neval),zdza(neval)
#endif
      real*8, dimension(:,:), allocatable :: thpkg2,rspl2,zspl2

      !--------------------------------------

#ifdef __DEBUG
      ictr = ictr + 1

      if(ictr.eq.51) then
         idebug=.TRUE.
         write(6,*) ' ictr = ',ictr
      else
         idebug=.FALSE.
      endif
#endif
      ier = 0

      if(size(R).le.3) then
         ier=201   ! a non-singular closed contour must be at least 4 pts
         return
      endif

      if(size(R).ne.size(Z)) then
         ier=202   ! R and Z array sizes must match
         return
      endif

      ntha = size(R)

      if((R(1).ne.R(ntha)).or.(Z(1).ne.Z(ntha))) then
         ier=203   ! contour must close
         return
      endif

      if(allocated(tha_pkg)) deallocate(tha_pkg,Rspl,Zspl,Lspl,wk,thmap)

      allocate(tha_pkg(ntha,4),Rspl(4,ntha),Zspl(4,ntha),Lspl(4,ntha))
      allocate(wk(ntha),thmap(ntha))

      do i=1,ntha
         thmap(i)=(i-1)*C2PI/(ntha-1)
      enddo

      !  create TH grid package for 1d spline evaluation (pspline lib call)
      call r8genxpkg(ntha,thmap,tha_pkg,1,0,0,ZERO,3,ier)
      if(ier.ne.0) then
         ier=205
         return
      endif

      !  create R spline & Z spline -- periodic splines

      Rspl(1,1:ntha)=R
      call r8cspline(thmap,ntha,Rspl,-1,ZERO,-1,ZERO,wk,ntha,idum,ier)
      if(ier.ne.0) then
         ier=206
         return
      endif
      
      Zspl(1,1:ntha)=Z
      call r8cspline(thmap,ntha,Zspl,-1,ZERO,-1,ZERO,wk,ntha,idum,ier)
      if(ier.ne.0) then
         ier=207
         return
      endif

      !  compute data for arc-length spline
      !  also find zone on large major radius side crossing Z=Zaxis

      call integ_wts(th0,w0)

      Lspl(1,1)=ZERO

      izero=0
      ifound=0
      if((R(ntha).gt.Raxis).and.(Z(ntha).eq.Zaxis)) then
         izero=ntha-1
         ifound=ntha
      endif

      do i=1,ntha-1
         wvec = (thmap(i+1)-thmap(i))*w0
         thvec = thmap(i) + (thmap(i+1)-thmap(i))*th0

         ! dR/dth evals

         call r8spvec(ict_der1,10,thvec,10,drdth,ntha,tha_pkg,Rspl,idum,ier)
         if(ier.ne.0) then
            ier=211
            return
         endif

         ! dZ/dth evals

         call r8spvec(ict_der1,10,thvec,10,dzdth,ntha,tha_pkg,Zspl,idum,ier)
         if(ier.ne.0) then
            ier=212
            return
         endif

         Lspl(1,i+1)=Lspl(1,i)
         do ii=1,10
            dldth = sqrt(drdth(ii)*drdth(ii) + dzdth(ii)*dzdth(ii))
            Lspl(1,i+1) = Lspl(1,i+1) + dldth*wvec(ii)
         enddo

         if((R(i).gt.Raxis).and.(Z(i).eq.Zaxis).and.(ifound/=ntha.or.i/=1)) then
            if(izero.gt.0) ier=213
            izero=i
            ifound=i
         else if(min(R(i),R(i+1)).gt.Raxis) then
            ztest = (Z(i)-Zaxis)*(Z(i+1)-Zaxis)
            if(ztest.lt.ZERO) then
               if(izero.gt.0) ier=213
               izero=i
            endif
         endif
         if(ier.ne.0) return

      enddo

      if(izero.eq.0) then
         ier=299  ! search failure, unexpected
         return
      endif

      !  setup L(TH)

      call r8cspline(thmap,ntha,Lspl,-1,ZERO,-1,ZERO,wk,ntha,idum,ier)
      if(ier.ne.0) then
         ier=215
         return
      endif

      Ltotal=Lspl(1,ntha)

#ifdef __DEBUG
      if(idebug) then
         do i=1,neval
            theval(i)=(i-1)*C2PI/(neval-1)
         enddo
         call r8spvec(ict_der1,neval,theval,neval,dldth1, &
              ntha,tha_pkg,Lspl,idum,ier)
         call r8spvec(ict_der1,neval,theval,neval,zdra, &
              ntha,tha_pkg,Rspl,idum,ier)
         call r8spvec(ict_der1,neval,theval,neval,zdza, &
              ntha,tha_pkg,Zspl,idum,ier)

         do i=1,neval
            dldth1(i)=dldth1(i)/(Ltotal/C2PI)
            dldth2(i)=sqrt(zdra(i)**2+zdza(i)**2)/(Ltotal/C2PI)
            zdra(i)=zdra(i)/(Ltotal/C2PI)
            zdza(i)=zdza(i)/(Ltotal/C2PI)
         enddo

         call r8_grafx2(theval,dldth1,dldth2,neval,'rad',' ', &
              'eq_arc_scrunch1 debug plot', &
              '[(dl/dth)/(Ltot/2pi)] from L spline', &
              '[(dl/dth)/(Ltot/2pi)] from R & Z splines')

         call r8_grafx3(theval,dldth2,zdra,zdza,neval,'rad',' ', &
              'eq_arc_scrunch1 debug plot', &
              '[(dl/dth)/(Ltot/2pi)] from R & Z splines', &
              'R & Z contributions separately.')
      endif
#endif

      !  find TH satisfying {R > Raxis, Z(TH)=Zaxis}

      eps=CSMALL*C2PI
      eta=CSMALL*Ltotal

      allocate(imask(ntha-2))
      imask=.FALSE.

      ! reuse thmap; the th working grid data is in tha_pkg(:,1)

      if(ifound.gt.0) then
         thzero=thmap(ifound)
      else

         call zridderx(1, imask, tha_pkg(izero,1), tha_pkg(izero+1,1), &
              eps, eta, eq_arc_fmidplane, thzero, ier, &
              1,Zaxis,1,Zdum,0)

      endif
      if(ier.ne.0) then
         ier=217
         return
      endif

      call r8spvec(ict_eval,1,thzero,1,Lzero,ntha,tha_pkg,Lspl,idum,iertmp)
      ier=max(ier,iertmp)
      call r8spvec(ict_eval,1,thzero,1,R(1),ntha,tha_pkg,Rspl,idum,iertmp)
      ier=max(ier,iertmp)
      call r8spvec(ict_eval,1,thzero,1,Z(1),ntha,tha_pkg,Zspl,idum,iertmp)
      ier=max(ier,iertmp)

      if(ier.gt.0) then
         ier=218
         return
      endif

      !  zero point = 2pi point established...

      R(ntha)=R(1)
      Z(ntha)=Z(1)

      imask(1)=.FALSE.

      !  get the rest...

      allocate(Ltarg(ntha-2),thfnd(ntha-2),th1(ntha-2),th2(ntha-2))
      allocate(zdumma(ntha-2))

      zLsum=Lzero
      do i=1,ntha-2
         zLsum = zLsum + Ltotal/(ntha-1)  ! equal arc increments
         if(zLsum.gt.Ltotal) then
            zLsum=zLsum-Ltotal
            izero=1
         endif
         Ltarg(i)=zLsum
         do
            if(Ltarg(i).gt.Lspl(1,izero+1)) then
               izero=izero+1
            else
               exit
            endif
         enddo
         th1(i)=tha_pkg(izero,1)
         th2(i)=tha_pkg(izero+1,1)
      enddo

      call zridderx(ntha-2, imask, th1, th2, eps, eta, eq_arc_findarc, &
           thfnd, ier, ntha-2, Ltarg, 1, zdumma, 1)

      if(ier.ne.0) then
         ier=219
         return
      endif

      call r8spvec(ict_eval,ntha-2,thfnd,ntha-2,R(2:ntha-1), &
           ntha,tha_pkg,Rspl,idum,iertmp)
      ier=max(ier,iertmp)

      call r8spvec(ict_eval,ntha-2,thfnd,ntha-2,Z(2:ntha-1), &
           ntha,tha_pkg,Zspl,idum,iertmp)
      ier=max(ier,iertmp)

      do i=1,ntha
         thmap(i)=(i-1)*C2PI/(ntha-1)
      enddo
      
      allocate(thpkg2(ntha,4),rspl2(4,ntha),zspl2(4,ntha))
      call r8genxpkg(ntha,thmap,thpkg2,1,0,0,ZERO,3,ier)

      rspl2(1,1:ntha) = r(1:ntha)
      zspl2(1,1:ntha) = z(1:ntha)

      call r8cspline(thmap,ntha,Rspl2,-1,ZERO,-1,ZERO,wk,ntha,idum,ier)
      call r8cspline(thmap,ntha,Zspl2,-1,ZERO,-1,ZERO,wk,ntha,idum,ier)

#ifdef __DEBUG

      if(idebug) then

         call r8spvec(ict_der1,neval,theval,neval,zdra, &
              ntha,thpkg2,Rspl2,idum,ier)
         call r8spvec(ict_der1,neval,theval,neval,zdza, &
              ntha,thpkg2,Zspl2,idum,ier)

         do i=1,neval
            dldth2(i)=sqrt(zdra(i)**2+zdza(i)**2)/(Ltotal/C2PI)
            zdra(i)=abs(zdra(i))/(Ltotal/C2PI)
            zdza(i)=abs(zdza(i))/(Ltotal/C2PI)
         enddo

         call r8_grafx3(theval,dldth2,zdra,zdza,neval,'rad',' ', &
              'eq_arc_scrunch1 debug plot', &
              '[(dl/dth)/(Ltot/2pi)] from R & Z spline set#2', &
              'R & Z contributions separately.')

      endif

#endif
      !  FFT

      nfft=8
      do
         if(nfft.ge.4*imom) exit
         nfft=nfft*2
      enddo

      allocate(thfft(nfft),rfft(nfft),zfft(nfft),zs1d(nfft),zc1d(nfft))

      do i=1,nfft
         thfft(i)=(i-1)*C2PI/nfft
      enddo

      call r8spvec(ict_eval,nfft,thfft,nfft,rfft, &
           ntha,thpkg2,rspl2,idum,ier)
      call r8spvec(ict_eval,nfft,thfft,nfft,zfft, &
           ntha,thpkg2,zspl2,idum,ier)

      call r8fftsc(Rfft,nfft,zs1d,zc1d,ier)

      rmc(0)=zc1d(1)/(2*nfft)
      do im=1,imom
         rmc(im)=zc1d(im+1)/nfft
         rms(im)=zs1d(im+1)/nfft
      enddo

      call r8fftsc(Zfft,nfft,zs1d,zc1d,ier)

      zmc(0)=zc1d(1)/(2*nfft)
      do im=1,imom
         zmc(im)=zc1d(im+1)/nfft
         zms(im)=zs1d(im+1)/nfft
      enddo

#ifdef __DEBUG
      if(idebug) then
         do i=1,neval-1
            call r8sincos(theval(i),imom,zs1d,zc1d)
            zdra(i)=rmc(0)
            zdza(i)=zmc(0)
            do im=1,imom
               zdra(i)=zdra(i)+rmc(im)*zc1d(im)+rms(im)*zs1d(im)
               zdza(i)=zdza(i)+zmc(im)*zc1d(im)+zms(im)*zs1d(im)
            enddo
         enddo
         zdra(neval)=zdra(1)
         zdza(neval)=zdza(1)

         call r8_grafc2(r,z,ntha, zdra,zdza,neval, 'm','m', &
              '[R,Z] contour & FFT expansion version', &
              'eq_arc_scrunch_mod debug plot', ' ')
      endif
#endif

      deallocate(thfft,rfft,zfft,zs1d,zc1d)
      deallocate(imask,Ltarg,thfnd,th1,th2,zdumma)
      deallocate(thpkg2,rspl2,zspl2)

    end subroutine eq_arc_scrunch1

end module eq_arc_scrunch_mod

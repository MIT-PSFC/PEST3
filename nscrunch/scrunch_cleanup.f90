subroutine scrunch_cleanup(xi,ns, rmmc,rmms,zmmc,zmms, kmom, &
     R0patch, xmom_xbrk,xmom_rtol, ierr)

  !  clean up moments representation:
  !   (a) modify R0(x) and Z0(x) so that 
  !       lim[x->0]d[R0]/dx = lim[x->0]d[Z0]/dx = 0
  !       by fiting a parabola that satisfies point and 1st derivative
  !       match at xmom_xbrk and 0 derivative on axis -- THIS MOVES THE
  !       MAGNETIC AXIS SLIGHTLY.
  !       **apply** similar patch to all EVEN numbered moments
  !       execute (a) iff R0patch=.TRUE.
  !
  !   (b) zero out moments that are below tolerance R0*xmom_rtol;
  !       save data on flux location inside of which the moments
  !       profile is below tolerance and should be smoothly set to
  !       zero.  Execute independent of value of R0patch
  !
  !  errors:
  !    ierr=2 -- xmom_xbrk value out of range
  !    ierr=3 -- xmom_rtol exceeds 1.0d-5
  !    ierr=1 -- error (unspecified)
  !
  !  inputs:
  !    ns -- no. of surfaces including axis and bdy
  !    kmom -- no. of highest moment, moments are (0:kmom)
  !    xi(1:ns) -- flux coord. grid [0,1]
  !    xmom_xbrk -- break point for parabola; also, transition width
  !       for zeroing of moments
  !    xmom_rtol -- tolerance that defines when moments values are too
  !       small and should be set to zero
  !
  !  output:
  !    ierr -- completion code (0 = OK)
  !
  !  in/out:
  !    rmmc,rmms,zmmc,zmms all dimensioned (0:kmom,1:ns) -- MODIFIED
  !    
  !-------------------------------------------
  implicit NONE

  integer, intent(in) :: ns,kmom
  real*8, intent(in) :: xi(ns)
  real*8, intent(inout) :: rmmc(0:kmom,ns),rmms(0:kmom,ns)
  real*8, intent(inout) :: zmmc(0:kmom,ns),zmms(0:kmom,ns)
  logical, intent(in) :: R0patch
  real*8, intent(in) :: xmom_xbrk,xmom_rtol
  integer, intent(out) :: ierr 
  !-------------------------------------------
  real*8 :: xcut_rc(0:kmom),xcut_rs(0:kmom),xcut_zc(0:kmom),xcut_zs(0:kmom)
  !-------------------------------------------

  call scrunch_cleanup_all(xi,ns, rmmc,rmms,zmmc,zmms, kmom, &
     R0patch, xmom_xbrk,xmom_rtol, xcut_rc, xcut_rs, xcut_zc, xcut_zs, ierr)

end subroutine scrunch_cleanup

subroutine scrunch_cleanup_all(xi,ns, rmmc,rmms,zmmc,zmms, kmom, &
     R0patch, xmom_xbrk,xmom_rtol, &
     xcut_rc, xcut_rs, xcut_zc, xcut_zs, ierr)

  !  clean up moments representation:
  !   (a) modify R0(x) and Z0(x) so that 
  !       lim[x->0]d[R0]/dx = lim[x->0]d[Z0]/dx = 0
  !       by fiting a parabola that satisfies point and 1st derivative
  !       match at xmom_xbrk and 0 derivative on axis -- THIS MOVES THE
  !       MAGNETIC AXIS SLIGHTLY.
  !   (b) zero out moments that are below tolerance R0*xmom_rtol;
  !       save data on flux location inside of which the moments
  !       profile is below tolerance and should be smoothly set to
  !       zero-- this step skipped if xmom_rtol=0
  !
  !  errors:
  !    ierr=2 -- xmom_xbrk value out of range
  !    ierr=3 -- xmom_rtol exceeds 1.0d-5
  !    ierr=1 -- error (unspecified)
  !
  !  inputs:
  !    ns -- no. of surfaces including axis and bdy
  !    kmom -- no. of highest moment, moments are (0:kmom)
  !    xi(1:ns) -- flux coord. grid [0,1]
  !    R0patch -- option to patch R0 (T/F)
  !    xmom_xbrk -- break point for parabola; also, transition width
  !       for zeroing of moments
  !    xmom_rtol -- tolerance that defines when moments values are too
  !       small and should be set to zero
  !
  !  output:
  !    ierr -- completion code (0 = OK)
  !    xcut_rc,xcut_rs,xcut_zc,xcut_zs all dimensioned (0:kmom)
  !       x >= xcut_rc(j) --> j'th moment of rmmc is unchanged;
  !       x <= xcut_rc(j)-xmom_xbrk --> j'th moment is zeroed
  !         (smooth transition in between)
  !
  !  in/out:
  !    rmmc,rmms,zmmc,zmms all dimensioned (0:kmom,1:ns) -- MODIFIED
  implicit NONE

  integer, intent(in) :: ns,kmom
  real*8, intent(in) :: xi(ns)
  real*8, intent(inout) :: rmmc(0:kmom,ns),rmms(0:kmom,ns)
  real*8, intent(inout) :: zmmc(0:kmom,ns),zmms(0:kmom,ns)
  logical, intent(in) :: R0patch
  real*8, intent(in) :: xmom_xbrk,xmom_rtol
  real*8, intent(out) :: xcut_rc(0:kmom),xcut_rs(0:kmom)
  real*8, intent(out) :: xcut_zc(0:kmom),xcut_zs(0:kmom)
  integer, intent(out) :: ierr
  !----------------------------------------------
  real*8 :: zxrel,zfac
  real*8 :: ztola   ! absolute tolerance inferred
  real*8 :: ztolp1
  integer :: i,im,ibrk
  real*8, parameter :: ZERO=0.0d0
  real*8 :: zmindiff,ztest
  !----------------------------------------------

  ierr=0

  if((xmom_xbrk.lt.0.1d0).or.(xmom_xbrk.gt.0.2d0)) then
     ierr=2
  endif

  if((xmom_rtol.lt.0.0d0).or.(xmom_rtol.gt.1.0d-5)) then
     ierr=3
  endif

  if(ierr.ne.0) return

  !  (1) patch the 0'th moments profiles -- borrow 0'th R sin and Z sin
  !      moments slots (reset to zero later)

  if(R0patch) then

     ztolp1 = 1 + xmom_rtol

     rmms(0,1:ns)=rmmc(0,1:ns)
     zmms(0,1:ns)=zmmc(0,1:ns)
     zmindiff = xi(ns)-xi(1)
     do i=1,ns
        ztest=abs(xi(i)-xmom_xbrk)
        if(ztest.lt.zmindiff) then
           zmindiff=ztest
           ibrk=i
        endif
     enddo

     call patch(xi,0,rmms(0,1:ns),rmmc(0,1:ns),ns,ibrk)
     call patch(xi,0,zmms(0,1:ns),zmmc(0,1:ns),ns,ibrk)

     rmms(0,1:ns)=0
     zmms(0,1:ns)=0

     !  (1a) apply similar patch to all even moments
     do im=2,kmom,2
        rmms(0,1:ns)=rmmc(im,1:ns)
        call patch(xi,im,rmms(0,1:ns),rmmc(im,1:ns),ns,ibrk)
        rmms(0,1:ns)=rmms(im,1:ns)
        call patch(xi,im,rmms(0,1:ns),rmms(im,1:ns),ns,ibrk)
        rmms(0,1:ns)=zmmc(im,1:ns)
        call patch(xi,im,rmms(0,1:ns),zmmc(im,1:ns),ns,ibrk)
        rmms(0,1:ns)=zmms(im,1:ns)
        call patch(xi,im,rmms(0,1:ns),zmms(im,1:ns),ns,ibrk)
        rmms(0,1:ns)=0
     enddo
  endif

  !  (2) zero small moments; 0'th moments not eligible

  ztola = xmom_rtol*rmmc(0,1)  ! R0(0)*tolerance
  if(ztola.le.ZERO) then
     xcut_rc=-1.0d0
     xcut_rs=-1.0d0
     xcut_zc=-1.0d0
     xcut_zs=-1.0d0
     return
  endif

  xcut_rc(0)=-1.0d0
  xcut_rs(0)=-1.0d0
  xcut_zc(0)=-1.0d0
  xcut_zs(0)=-1.0d0

  do im=1,kmom
     call zpatch(rmmc(im,1:ns),xcut_rc(im))
     call zpatch(rmms(im,1:ns),xcut_rs(im))
     call zpatch(zmmc(im,1:ns),xcut_zc(im))
     call zpatch(zmms(im,1:ns),xcut_zs(im))
  enddo

  !----------------------------------------------------------
  contains
    subroutine zpatch(zmom,zcut)
      !
      !  zero patch a single moments profile
      !
      real*8, intent(inout) :: zmom(ns)
      real*8, intent(out) :: zcut

      real*8 :: dfdx,zdelx,a,c

      ibrk=0
      do i=2,ns
         if(abs(zmom(i)).gt.ztola) then
            ibrk=i-1
            exit
         endif
      enddo

      if(ibrk.eq.0) then
         zcut=2.0d0  ! zero everywhere
         zmom=0.0d0

      else if(xi(ibrk).le.xmom_xbrk) then
         zcut=-1.0d0
         return  ! no change, non-negligible value too close to axis

      else
         zcut=xi(ibrk)
         dfdx = (zmom(ibrk+1)-zmom(ibrk))/(xi(ibrk+1)-xi(ibrk))
         a = dfdx/(2*xmom_xbrk)
         c = zmom(ibrk)-dfdx*xmom_xbrk/2
         do i=1,ibrk-1
            if(xi(i).le.zcut-xmom_xbrk) then
               zmom(i)=0
            else
               zdelx=(xi(i)-(zcut-xmom_xbrk))
               zxrel=zdelx/xmom_xbrk
               zfac=10 + zxrel*(6*zxrel - 15)
               zfac= zfac*zxrel*zxrel*zxrel  ! f(x)=6x**5 -15x**4 +10x**3
               !  f(0)=f'(0)=f''(0)=f'(1)=f''(1)=0; f(1)=1
               zdelx=zcut-xi(i)
               zmom(i)=zfac*(a*zdelx*zdelx+c)
            endif
         enddo
      endif

    end subroutine zpatch
 
    !-----------------------------------------------------

    subroutine patch(x,imoment,f_orig,f_new,inum,ibrk)

      integer, intent(in) :: inum
      integer, intent(in) :: imoment  ! 0,2,4,...
      real*8, intent(in) :: x(inum)
      real*8, intent(in) :: f_orig(inum)
      real*8, intent(out) :: f_new(inum)
      integer, intent(in) :: ibrk

      !  parabolic extrapolation to axis ==> d/dx[0'th moments] --> 0

      !-----------------------------
      real*8 :: f1,f2,f3,x1,x2,x3,dfdx,a,c
      !-----------------------------

      f_new = f_orig

      !  If 1st two points are identically zero, assume that no further
      !  patching is needed...

      if((f_new(1).eq.ZERO).and.(f_new(2).eq.ZERO)) return

      f1=f_orig(ibrk+1)-f_orig(ibrk)
      x1=x(ibrk+1)-x(ibrk)

      f2=f_orig(ibrk+2)-f_orig(ibrk)
      x2=x(ibrk+2)-x(ibrk)

      !  local derivate from parabol through f(ibrk:ibrk+2)...
      !  f(0)=0, f(x1)=f1, f(x2)=f2

      dfdx = (f1*x2*x2-f2*x1*x1)/(x1*x2*(x2-x1))
      
      f3=f_orig(ibrk)
      x3=x(ibrk)

      if(imoment.eq.0) then

         !  quadratic to satisfy f(x3)=f3, f'(x3)=dfdx, f'(0)=0
         !  f(x) = a*x*x + b*x + c; f'(x) = 2*a*x + b
         !    f'(0) = 0 --> b = 0
         !    f'(x3)=dfdx --> a = dfdx/(2*x3)
         !    f(x3)=f3 --> c = f3 - a*x3*x3 = f3 - dfdx*x3/2

         a = dfdx/(2*x3)
         c = f3 - dfdx*x3/2

         do i=1,ibrk-1
            x1=x(i)
            f_new(i)=a*x1*x1 + c
         enddo

      else

         !  cubic to satisfy f(x3)=f3, f'(x3)=ddfdx, f'(0)=0, f(0)=0
         !  f(x) = a*x*x*x + c*x*x + b*x + d = 0; f'(x) = 3*a*x*x + 2*c*x + b
         !         f(0)=0 --> d = 0    f'(0) = 0 --> b = 0
         !
         !         a*x3**3 + c*x3**2 = f3
         !       3*a*x3**2 + 2*c*x3 = dfdx

         a = (dfdx - 2*f3/x3)/(x3*x3)
         c = (3*f3/x3 - dfdx)/x3

         do i=1,ibrk-1
            x1=x(i)
            f_new(i)=x1*x1*(a*x1 + c)
         enddo

      endif

    end subroutine patch

end subroutine scrunch_cleanup_all

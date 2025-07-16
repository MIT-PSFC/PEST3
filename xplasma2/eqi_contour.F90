subroutine eqi_contour(xp, yp, cbdy, errscale, cbdy_got, psiav, istat, ier)
 
  ! Using the function Psi(x,y) (already set up in bicubic spline pspl
  ! in the eqi_geq_mod module), compute values (xs, ys) where
  ! f(xs, ys) = f(xp, yp).
  !
  ! values are set at points where 
  !    atan2((ys(j)-y0),(xs(j)-x0)) = thadapt(j,isw)
  !    ...where thadapt(1:nthad,isw) is a grid spanning [0,2pi]
  !       declared in eqi_geq_mod
  !
  ! Method: Hamilton's equation dx/dt=dF/dy; dy/dt=-dF/fx are integrated
  ! starting from (xp, yp) until the cummulative angle between (x-x0, y-y0)
  ! and (xnew-x0, ynew-y0) reaches 2*pi, where (xnew, ynew) are new coordinates
  ! along the trajectory.
  ! 
  ! As a second step, the (xnew, ynew) coordinates are interpolated onto
  ! the thadapt angle grid to yield the (xs, ys).
  !
  ! We are looking for a set of closed contours.  If a contour leads
  ! off the grid, set an error code.
 
  ! dmc modifications -- integrand direction:  counterclockwise around
  !   (x0,y0); integrand "speed" normalized s.t. t is O(1) to go all the
  !   way around-- error flag is set if dF/dx=dF/dy=0.
 
  ! ier=123 exit denotes one of a variety of contour features:
  !   (a) direction of contour reversed: d/dt(arctan((y-y0)/(x-x0)))
  !   (b) contour escaped (R,Z) box
  !   (c) contour landed on grad(psi)=0 point
 
  ! ier=123 is **normal** for a contour corresponding to an open field line
 
  !------------------------------
  !  DMC note: this routine has been carefully tested and has a number
  !  of complexities.  For example, the zeroing-in on the final end point
  !  of a closed contour has to be, and is, carefully handled; the end
  !  point treatment can affect the spline contour accuracy test.
  !
  !  If you work in here, Handle with care!
  !------------------------------

  use eqi_geq_mod
 
  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)
  real(r8), intent(in) :: xp, yp ! starting point
  real(r8), intent(in) :: cbdy   ! min. curvature parameter
  real(r8), intent(in) :: errscale ! accuracy test scale
  real(r8), intent(out) :: cbdy_got ! min. curvature observed
  real(r8), intent(out) :: psiav ! avg Psi value along contour
  integer, intent(out) :: istat
  !                        =0 -- closed surface, curvature OK
  !                              (psiav,cbdy_got) computed
  !                        =1 -- closed surface, curvature too sharp
  !                              (psiav,cbdy_got) computed
  !                        =2 -- open surface or grad(Psi)=0 encountered
  !                              (psiav computed; cbdy_got=0)
  !                        =3 -- some other error
  !                              (psiav=0; cbdy_got=0)
  integer, intent(out) :: ier
 
  real(r8) :: x0,y0,xxmin,xxmax,yymin,yymax
  !
  !  foreach i, |f(xs(i),ys(i)) - f(xp,yp)| < 0.0001*errscale
  !     recommend errscale = max(f)-min(f)
  !
  !--------------------------------------------
 
  integer, parameter :: neq=2, mf=10, itmax=10
  integer itol, itask, istate, iopt, lrw, iwork(20), liw, ieq, iter, iout
  integer irevers,igrad0,ilsode
  real(r8) atol, rtol, rwork(20 + 16*neq)
 
  real(r8), dimension(:), allocatable :: psis
 
  real(r8), parameter :: twopi=6.28318530717958623_r8, eps=1.e-10_r8
  integer, parameter :: nsteps0 = 2500   ! no. of steps for ns=32
  integer :: nsteps                      ! actual number of integration steps
  integer i, iok, npoints, npoints_save, ichk, iang, iangp
  real(r8), dimension(:), allocatable :: aa
  real(r8), dimension(:,:), allocatable :: xy
  real(r8) :: dx(neq), dxnew(neq), xx(neq), xx0(neq)
  real(r8) :: lenscale, t, dt, tout, angle, a, const, vecprod, scaprod
 
  real(r8) :: diff, aver, corr, dtest, radius, minrat, radmax
  real(r8) :: rminc,rmaxc
 
  type(ezspline1_r8), save :: xspl, yspl
 
  character*120 msgbuf
 
  external eqi_xydot, eqi_xyjac
 
#ifdef __DEBUG
  integer, parameter :: lun_dbg = 6   ! used iff __DEBUG macro is set
  integer :: icount_calls = 0
  integer :: icount_debug = -1

  !-------------------------------------------
  icount_calls = icount_calls + 1
  if(icount_calls.eq.icount_debug) then
     write(lun_dbg,*) ' eqi_contour: reached icount_debug: ',icount_debug
  endif
#endif
  !-------------------------------------------

  ier = 0
  iopt = 0
 
  istat = 3

  cbdy_got = 0.0_r8
  psiav = 0.0_r8

  if( errscale <= 0.0_r8 ) then
     ier=999
     go to 1001
  endif
 
  x0=r0
  y0=z0

  rminc=xp
  rmaxc=xp

  xxmin=pspl%x1(1)
  xxmax=pspl%x1(pspl%n1)

  yymin=pspl%x2(1)
  yymax=pspl%x2(pspl%n2)
 
!  DMC: dt halved; nsteps doubled: Feb 2008

  dt = twopi/max(64._r8,2.0_r8*nthad0)
  iter=0
!
  nsteps = 2*nsteps0*max(1.0_r8,nthad0/32.0_r8)
  allocate(xy(neq,nsteps),aa(nsteps))
!
!  initialize integrand
!
  lnorm=sqrt((xp-x0)**2+(yp-y0)**2)
  lmax=1.0_r8/errscale
  errflag=0

  radmax = 10*lnorm
!
!  check direction from eqi_xydot.  If the reverse of the expected
!  direction, change signs.
!
!  starting at (xp,yp), we want to "orbit" around (x0,y0) in a
!  counterclockwise direction  (dmc 5 Feb 2001)
!
  isign_xy=1      ! element of the "eqi_geq_mod" module used by eqi_xydot...
  call eqi_xydot(2, 0.0_r8, (/xp, yp/), dxnew)
  if(errflag.ne.0) then
     istat=2   ! grad(Psi) = 0
     call ezspline_interp(pspl,xp,yp,psiav,iok)
     call ezspline_error(iok)
     ier=123
     go to 1000
  endif
  dtest = (xp-x0)*dxnew(2) - (yp-y0)*dxnew(1)
  if(dtest.lt.0.0_r8) isign_xy = -isign_xy
!
!  in this loop, timesteps are shortened until desired accuracy
!  is achieved.  Mod DMC Mar 2008: this is done even for open contours.
!
  do
     istat=0
     ier=0
     iter=iter+1
     if(iter.gt.itmax) then
        istat=3
        ier=999
        go to 1000
     endif
 
     curvad(1:nthad,isw) = radmax
     iang=1

     dt = dt/2.0_r8
 
     ! lsode stuff
     rtol=eps
     atol=eps*errscale
     itask = 1
     istate = 1
     itol = 1
     iopt = 0
     lrw  = 20 + 16*neq
     liw = 20
 
     t = 0.0_r8
     xy = 0.0_r8
     xx0 = (/x0, y0/)
     xx  = (/xp, yp/)
     i = 1

     npoints= 1
     npoints_save = 0  ! used for end pt.

     lenscale = sqrt( (xp-x0)**2 + (yp-y0)**2 )
 
     angle = 0.0_r8
     const = 1.0_r8
     aa = 0.0_r8
     xy(:,1) = xx
 
     iout = 0
     irevers = 0
     ilsode = 0
     igrad0 = 0
 
     minrat = radmax
 
     do while (npoints < nsteps .and. abs(const*dt) > eps)
        dx = xx - xx0
        tout = t + const*dt
        call lsode_r8 (eqi_xydot, neq, xx, t, tout, itol, rtol, atol, itask, &
             & istate, iopt, rwork, lrw, iwork, liw, eqi_xyjac, mf)
        if(istate<=0) then
           ilsode=1
           exit
        endif
        if(errflag.ne.0) then
           igrad0=1
           exit
        endif
        if((xx(1).lt.xxmin).or.(xx(1).gt.xxmax).or. &
             & (xx(2).lt.yymin).or.(xx(2).gt.yymax)) then
           iout=1
           exit
        endif
        dxnew = xx - xx0
        vecprod = dx(1)*dxnew(2)-dx(2)*dxnew(1)
        scaprod = dx(1)*dxnew(1)+dx(2)*dxnew(2)
        a = atan2(vecprod, scaprod)
! (detect but do not prevent angular velocity reversal -- dmc Feb 2001)
!!!        if(a < 0) then
!!!           a = a + twopi
!!!        endif
        angle = angle + a
        if(abs(angle) > twopi) then
           if(const.ge.0.99_R8) then
              !  enter here when first past end point; step has not yet been
              !  reduced to search for exact end point

              !  npoints points to start of current step, which is now known
              !  to be within a const=1 step of the boundary closing
              npoints_save = npoints

              ! end point curvature: evaluate at *start* of zeroing-in to
              ! angle=twopi end point of closed contour
              if(npoints.ge.10) then
                 call r8_trirad(xy(:,npoints-1),xy(:,npoints),xy(:,2), &
                      dxnew,radius,ichk)
                 if(ichk.eq.0) then
                    curvad(1,isw)=min(curvad(1,isw),radius)
                    curvad(nthad,isw)=min(curvad(nthad,isw),radius)
                    if(nthad.gt.(iang+1)) then
                       curvad(iang+1:nthad-1,isw)=curvad(iang,isw)
                    endif
                    minrat=min(minrat,radius)
                 endif
              endif
           endif
           ! back off
           t = t - const*dt
           const = const/2.0_r8
           angle = angle - a
           xx = xy(:,npoints)
           ! reset lsode state
           istate = 1
        else
           if(angle.le.aa(npoints)) then
              irevers=1
              exit
           endif
           npoints = npoints + 1 ! forward counter
           aa(npoints) = angle

           iangp=iang
           do while(angle.gt.thadapt(iang+1,isw))
              iang = iang + 1
           enddo
           if(iangp.lt.(iang-1)) then
              ! deal with skipped index space
              curvad(iangp+1:iang-1,isw)=curvad(iangp,isw)
           endif

           rminc=min(rminc,xx(1))
           rmaxc=max(rmaxc,xx(1))
           xy(:, npoints) = xx
           if((npoints.ge.3).and.(const.ge.0.99_R8)) then
              ! inaccurate to evaluate curvature at end point where const << 1
              call r8_trirad(xy(:,npoints-2),xy(:,npoints-1),xy(:,npoints), &
                   dxnew,radius,ichk)
              if(ichk.eq.0) then
                 curvad(iang,isw)=min(curvad(iang,isw),radius)
                 curvad(iang+1,isw)=min(curvad(iang+1,isw),radius)
                 minrat=min(minrat,radius)
              endif
           endif
        endif
        i = i + 1
     enddo
     if(iout.eq.1) then
#ifdef __DEBUG        
        write(lun_dbg,*) '%eqi_contour:  contour escapes mapped region.'
#endif
        istat=2
     else if(irevers.eq.1) then
#ifdef __DEBUG        
        write(lun_dbg,*) '%eqi_contour:  contour angular velocity reversed.'
#endif
        istat=2
     else if(igrad0.eq.1) then
#ifdef __DEBUG        
        write(lun_dbg,*) '%eqi_contour:  contour hit grad(psi)=0'
#endif
        istat=2
     else if(ilsode.eq.1) then
#ifdef __DEBUG        
        msgbuf=' '
        write(msgbuf,*) '%eqi_contour:  LSODE error', istate
        write(lun_dbg,*) trim(msgbuf)
#endif
        istat=3
        iter=itmax
     else if(npoints < nsteps) then
        npoints = npoints + 1
        aa(npoints) = twopi
        xy(:, npoints) = xx
        lnorm=(rmaxc-rminc)/2
        minrat=minrat/lnorm
        curvad(1:nthad,isw) = curvad(1:nthad,isw)/lnorm
        if(minrat.lt.cbdy) then
           istat=1
#ifdef __DEBUG        
           msgbuf=' '
           write(msgbuf,*) ' *** minimum normalized curvature = ',minrat
           write(lun_dbg,*) trim(msgbuf)
           msgbuf=' '
           write(msgbuf,*) &
                & '%eqi_contour:  contour violates sharp curvature limit.'
           write(lun_dbg,*) trim(msgbuf)
#endif
        endif
     else
#ifdef __DEBUG        
        msgbuf=' '
        write(msgbuf,*) '%eqi_contour:  nsteps limit reached:', nsteps
        write(lun_dbg,*) trim(msgbuf)
#endif
        istat=3
        iter=itmax
     endif

     if(iter.eq.itmax) cycle
 
     if(istat.eq.2) then
        if(allocated(psis)) deallocate(psis)
        allocate(psis(npoints))
        call ezspline_interp(pspl, npoints, &
             xy(1,1:npoints), xy(2,1:npoints), psis, iok)
        call ezspline_error(iok)

        aver = sum(psis)/real(npoints, r8)
        psiav = aver

     else
        ! contour is closed...
        ! interpolate onto equal angle grid
 
        ! avoid closely packed points at end of contour
        ! these can screw up the spline boundary condition

        if(npoints_save.gt.0) then
           npoints = npoints_save
        endif
        xy(:,npoints)=xy(:,1)  ! enforce exact closing
        aa(npoints)=aa(1)+twopi

        call ezspline_init(xspl, npoints, (/-1, -1/), iok)   ! periodic BC
        call ezspline_error(iok)
        call ezspline_init(yspl, npoints, (/-1, -1/), iok)   ! periodic BC
        call ezspline_error(iok)
 
        xspl%x1 = aa(1:npoints)
        yspl%x1 = aa(1:npoints)
 
        xy(:,npoints) = xy(:,1) ! make periodic
 
        call ezspline_setup(xspl, xy(1,1:npoints), iok)
        call ezspline_error(iok)
        call ezspline_setup(yspl, xy(2,1:npoints), iok)
        call ezspline_error(iok)
 
        call ezspline_interp(xspl, nthad, thadapt(1:nthad,isw), &
             radapt(1:nthad,isw), iok)
        call ezspline_error(iok)
        call ezspline_interp(yspl, nthad, thadapt(1:nthad,isw), &
             zadapt(1:nthad,isw), iok)
        call ezspline_error(iok)
 
        call ezspline_free(xspl, iok)
        call ezspline_error(iok)
        call ezspline_free(yspl, iok)
        call ezspline_error(iok)
 
        ! check
 
        if(allocated(psis)) deallocate(psis)
        allocate(psis(nthad))
        call ezspline_interp(pspl, nthad, &
             radapt(1:nthad,isw), zadapt(1:nthad,isw), psis, iok)
        call ezspline_error(iok)

        aver = sum(psis)/real(nthad, r8)
        psiav = aver
     endif

     diff = maxval(abs(psis-aver))
     if(diff/errscale > 0.0001) then
        ier = 124  ! redo: insufficient accuracy: too much Psi variation.
     else
        if(istat.gt.1) ier = 123
        exit
     endif
 
  enddo
 
  ! clean up
 
1000 continue
 
  deallocate(xy,aa)
  if(allocated(psis)) deallocate(psis)

1001 continue
  if(istat.le.1) then
     cbdy_got=minrat
  else
     cbdy_got=0.0
  endif

  return
 
end subroutine eqi_contour
 
 
subroutine eqi_xydot(neq, t, y, ydot)
 
  use eqi_geq_mod
 
  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)
  integer neq, iok
  real(r8) t, y(neq), ydot(neq)
  real(r8) df(2)
  real(r8) v,zero,norm,rr,zz
  integer :: in1,in2
 
  !-----------------
  rr=y(1)
  zz=y(2)

  in1 = pspl%n1
  in2 = pspl%n2

  if(rr.lt.pspl%x1(1)) then
     rr=pspl%x1(1)
     istep_out = istep_out + 1
  else if(rr.gt.pspl%x1(in1)) then
     rr=pspl%x1(in1)
     istep_out = istep_out + 1
  endif

  if(zz.lt.pspl%x2(1)) then
     zz=pspl%x2(1)
     istep_out = istep_out + 1
  else if(zz.gt.pspl%x2(in2)) then
     zz=pspl%x2(in2)
     istep_out = istep_out + 1
  endif

  call ezspline_gradient(pspl, rr, zz, df, iok)
  call ezspline_error(iok)
 
  ydot(1) = -df(2)*isign_xy
  ydot(2) = +df(1)*isign_xy
 
  zero = 0
  v = sqrt(ydot(1)**2+ydot(2)**2)
  if(v.le.zero) then
     errflag=1
  else
     norm=min(lmax,lnorm/v)
     ydot=norm*ydot
  endif
 
end subroutine eqi_xydot
 
subroutine eqi_xyjac
end subroutine eqi_xyjac

subroutine eqi_xtrz(ierr)

  use xplasma_definitions
  use eqi_rzbox_module

  !  given the interior geometry, extrapolate the exterior...
  !  axisymmetry assumed...

  !  this routine creates an extrapolated rho grid, and an extrapolated
  !  set of nested surfaces, based on legacy TRANSP "bloat" software
  !  (D. McCune 17 Jan 2000)

  implicit NONE

  integer ierr                      ! completion code, 0=OK

  !------------------------------------------------

  logical symflag,axisymm,scrapeoff

  integer inbblo                    ! size of temporary extrap. grid
  real*8 zdrho

  real*8 zrmin,zrmax,zzmin,zzmax,zfacx0,zfacx1,zfacx,zfac,rmax_pt
  real*8 zrmin_pl,zrmax_pl,zzmin_pl,zzmax_pl,raxis,zaxis,zlen,zlentot
  real*8 zdlav,zdlmax,zdr,zdz,zfaclim

  real*8, dimension(:), allocatable :: xiblo,xi_bdy,th_bdy,th_arc,zdlen_bdy
  real*8, dimension(:,:), allocatable :: rz_bdy,th_arc_pkg,thspl,rz_ntk
  real*8, dimension(:,:,:), allocatable :: rmc,ymc,rmcx,ymcx
  real*8, dimension(:), allocatable :: th_arc_ntk,th_xpl_ntk,xi_bdy_ntk
  real*8, dimension(:), allocatable :: th_arc_exp,th_xpl_exp
  real*8, dimension(:), allocatable :: th_xpl_fin,th_arc_fin,xi_bdy_fin
  real*8, dimension(:,:), allocatable :: rz_fin

  integer lcentr,nzones,inzp1
  integer, parameter :: mimom = 32  ! no higher than 32

  integer i,ii,ith,ith_rmax,intk,id,id_Rgrid,id_Zgrid,id_R,id_Z,ids(2)
  integer :: ntheta,nthetax,nrhox
  integer :: iaxis_rhox,iaxis_theta,iaxis_thetax,iadd2pi,iertmp
  integer iexpand,idexp,iexp,imark_endpt,idum

  real*8 :: rho_bdy,rho_vac,zdth_mark

  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: c2pi = 6.2831853071795862D+00

  data symflag/.FALSE./

  !------------------------------------------------

  ierr=0
  zfacx1=1.600d0  ! initial grid expansion factor... will be expanded
  !                 if needed to cover SOL (R,Z) rectangle

  rho_bdy = 1.00d0

  call xplasma_global_info(sp,ierr, axisymm=axisymm, scrapeoff=scrapeoff)
  if(ierr.ne.0) return

  if(.not.axisymm) then
     ierr=107
     call xplasma_errmsg_append(sp,' ...error in eqi_xtrz subroutine.')
     return
  endif

  if(.not.scrapeoff) then
     ierr=109
     call xplasma_errmsg_append(sp,' ...error in eqi_xtrz subroutine.')
     return
  endif

  !  just exit now if the extrapolation has been done already...

  call xplasma_find_item(sp,'__R_extrap',id,ierr, nf_noerr=.TRUE.)
  if(ierr.ne.0) return
  if(id.ne.0) return

  !  OK................................................... proceed...
  !  ad hoc extrapolation step

  !  __R_extrap(rho,theta) and __Z_extrap(rho,theta) will be created.
  !  extrapolation grids (__rhox) & (__thetax) will also be defined

  call xplasma_common_ids(sp,ierr, id_R=id_R, id_Z=id_Z)
  if(ierr.ne.0) return

  call xplasma_prof_info(sp,id_R,ierr, gridId1=iaxis_theta)
  if(ierr.ne.0) return

  call xplasma_grid_size(sp,iaxis_theta,ntheta,ierr)
  if(ierr.ne.0) return

  call xplasma_RZminmax_extended(sp, zrmin, zrmax, zzmin, zzmax, ierr, &
       id_Rgrid=id_Rgrid, id_Zgrid=id_Zgrid)
  if(ierr.ne.0) return

  call xplasma_RZminmax_plasma(sp, zrmin_pl,zrmax_pl, zzmin_pl,zzmax_pl, ierr)
  if(ierr.ne.0) return

  call xplasma_mag_axis(sp, ierr, raxis=raxis, zaxis=zaxis)
  if(ierr.ne.0) return

  !-------------------
  !  OK -- will use arc length theta for bloat expansion
  !  error checking is done... error after this point would be a bug

  allocate(xi_bdy(ntheta),th_bdy(ntheta),th_arc(ntheta),zdlen_bdy(ntheta-1))
  allocate(rz_bdy(ntheta,2),th_arc_pkg(ntheta,4),thspl(2,ntheta))
  allocate(th_arc_exp(2*ntheta))

  ids(1)=id_R; ids(2)=id_Z
  xi_bdy = rho_bdy
  call xplasma_grid(sp,iaxis_theta,th_bdy,iertmp)

  call xplasma_eval_prof(sp,ids, &
       xplasma_theta_coord,th_bdy, xplasma_rho_coord,xi_bdy, &
       rz_bdy,iertmp)

  zlentot=ZERO
  zdlmax=ZERO
  ith_rmax=1
  rmax_pt=rz_bdy(1,1)
  do ith=1,ntheta-1
     if(rz_bdy(ith,1).gt.rmax_pt) then
        ith_rmax=ith
        rmax_pt=rz_bdy(ith,1)
     endif
     zdr=rz_bdy(ith+1,1)-rz_bdy(ith,1)
     zdz=rz_bdy(ith+1,2)-rz_bdy(ith,2)
     zdlen_bdy(ith)=sqrt(zdr*zdr+zdz*zdz)
     zdlmax=max(zdlmax,zdlen_bdy(ith))
     zlentot=zlentot+zdlen_bdy(ith)
  enddo
  zdlav=zlentot/(ntheta-1)

  ! Arc-length theta 0 at ith_rmax

  ith=ith_rmax-1
  iadd2pi=0
  th_arc_pkg(1,1)=0
  th_arc_exp(1)=0
  zlen=0.0d0
  iexpand=0
  iexp=1

  do i=1,ntheta-1
     ith=ith+1
     th_arc(ith)=c2pi*(zlen/zlentot)
     thspl(1,i)=th_bdy(ith)+iadd2pi*c2pi
     if(ith.eq.ntheta) then
        ith=1
        iadd2pi=1
        th_arc(1)=th_arc(ntheta)
     endif
     zlen=zlen+zdlen_bdy(ith)
     th_arc_pkg(i+1,1)=c2pi*(zlen/zlentot)
     if(zdlen_bdy(ith).ge.2.500d0*zdlav) then
        idexp = (zdlen_bdy(ith)/zdlav) - 1
        iexpand = iexpand + idexp
        do ii=1,idexp
           iexp=iexp+1
           th_arc_exp(iexp) = th_arc_pkg(i,1) + &
                ii*(th_arc_pkg(i+1,1)-th_arc_pkg(i,1))/(idexp+1)
        enddo
     endif
     iexp=iexp+1
     th_arc_exp(iexp)=th_arc_pkg(i+1,1)
  enddo
  thspl(1,ntheta)=thspl(1,1)+c2pi
  th_arc(ntheta)=th_arc(1)+c2pi

  call r8genxpkg(ntheta,th_arc_pkg(1:ntheta,1),th_arc_pkg,1,0,0,ZERO,3,iertmp)
  call r8mkspline(th_arc_pkg(1:ntheta,1),ntheta,thspl,-1,ZERO,-1,ZERO,idum, &
       iertmp)

  !-------------------
  !  bloat gridsize...

  call r8ntk(intk)
  allocate(th_arc_ntk(intk),th_xpl_ntk(intk),xi_bdy_ntk(intk),rz_ntk(intk,2))

  do ith=1,intk
     th_arc_ntk(ith)=(ith-1)*c2pi/intk  ! 0:2pi-dtheta
     xi_bdy_ntk(ith)=rho_bdy
  enddo

  !  transform back to xplasma splines (w/possible 2pi shift) from arclength
  !  angle parameter...

  call r8vecspline((/1,0,0/),intk,th_arc_ntk,intk, &
       th_xpl_ntk,ntheta,th_arc_pkg,thspl,idum,iertmp)

  !  corresponding (R,Z) locations:

  call xplasma_eval_prof(sp,ids, &
       xplasma_theta_coord,th_xpl_ntk, xplasma_rho_coord,xi_bdy_ntk, &
       rz_ntk,iertmp)

  !  now form final grid & corresponding arclength thetas

  nthetax=ntheta+iexpand
  allocate(th_xpl_exp(nthetax),th_arc_fin(nthetax),th_xpl_fin(nthetax))

  call r8vecspline((/1,0,0/),nthetax,th_arc_exp(1:nthetax),nthetax, &
       th_xpl_exp,ntheta,th_arc_pkg,thspl,idum,iertmp)

  zdth_mark=2*c2pi
  do iexp=1,nthetax
     if(abs(th_xpl_exp(iexp)-th_bdy(ntheta)).lt.zdth_mark) then
        zdth_mark=abs(th_xpl_exp(iexp)-th_bdy(ntheta))
        imark_endpt=iexp
     endif
  enddo

  if(imark_endpt.eq.nthetax) then
     th_arc_fin=th_arc_exp(1:nthetax)
     th_xpl_fin=th_xpl_exp
     th_xpl_fin(1)=th_bdy(1)
     th_xpl_fin(nthetax)=th_bdy(ntheta)
  else
     th_xpl_fin(1)=th_bdy(1)
     th_arc_fin(1)=th_arc_exp(imark_endpt)
     ii=1
     ith=imark_endpt
     do i=2,nthetax-1
        ith=ith+1
        if(ith.eq.nthetax) ith=1
        ii=ii+1
        if(th_xpl_exp(ith).gt.th_bdy(ntheta)) then
           th_xpl_fin(ii)=th_xpl_exp(ith)-c2pi
        else
           th_xpl_fin(ii)=th_xpl_exp(ith)
        endif
        th_arc_fin(ii)=th_arc_exp(ith)
     enddo
     th_xpl_fin(nthetax)=th_bdy(ntheta)
     th_arc_fin(nthetax)=th_arc_fin(1)
  endif

  allocate(rz_fin(nthetax,2),xi_bdy_fin(nthetax))

  xi_bdy_fin=rho_bdy
  if(iexpand.gt.0) then
     call xplasma_eval_prof(sp,ids, &
          xplasma_theta_coord,th_xpl_fin, xplasma_rho_coord,xi_bdy_fin, &
          rz_fin,iertmp)
  else
     rz_fin = rz_bdy
  endif

  !-------------------
  call xplasma_author_set(sp,xplasma_xmhd,iertmp)

  zdrho= 0.04d0   ! extrapolation rho step size...
  zfacx0=zfacx1
  zfaclim=4.0d0

10 continue
  zfacx=zfacx0*max((zrmax-zrmin)/(zrmax_pl-zrmin_pl), &
       (zzmax-zzmin)/(zzmax_pl-zzmin_pl))

  zfacx=min(zfacx,zfaclim)  ! sanity limit

  rho_vac=zfacx

  inbblo=5+(rho_vac-rho_bdy)/zdrho

  lcentr=1
  nzones=2
  inzp1=nzones+1

  allocate(xiblo(inbblo))
  allocate(rmc(inzp1,0:mimom,2))
  allocate(ymc(inzp1,0:mimom,2))
  allocate(rmcx(inbblo,0:mimom,2))
  allocate(ymcx(inbblo,0:mimom,2))

  !  form temporary grid for bloat routine

  xiblo(1)=0.00d0
  xiblo(2)=0.70d0
  xiblo(3)=rho_bdy
  do i=4,inbblo
     xiblo(i)=xiblo(i-1)+zdrho
  enddo
  rho_vac=xiblo(inbblo)
  nrhox=inbblo-2

  rmc=ZERO
  ymc=ZERO
  rmcx=ZERO
  ymcx=ZERO

  rmc(lcentr,0,1)=raxis
  ymc(lcentr,0,1)=zaxis

  !  setup startup grid -- last two surfaces inside the plasma
  !  complete interior moments

  call eqi_xtrzf(symflag,intk,xiblo(1:inzp1),rmc,ymc, &
       inzp1,mimom,rz_ntk)

  !  extrapolate from plasma bdy out to user specified "rho" ~ r/a

  call r8bloata0(symflag,1,nzones,xiblo(1),inbblo,inbblo,inzp1, &
       rmc(1,0,1),ymc(1,0,1),mimom,mimom,6, &
       rmcx(1,0,1),ymcx(1,0,1),ierr)

  if(ierr.ne.0) then
     ierr=9999
     call xplasma_errmsg_append(sp,' ?? eqi_xtrz:  BLOAT failure.')
     call xplasma_errmsg_append(sp,'    r8bloata0 call failed in eqi_xtrz.')
     go to 1000
  endif

  call eqi_xtrz_test(zrmin,zrmax,zzmin,zzmax, &
       xiblo(1),inbblo,rmcx(1,0,1),ymcx(1,0,1),mimom,zfac,ierr)
  if(ierr.ne.0) then

  !  try to increase the expansion range once, then give up

     if(zfacx0.gt.zfacx1) then
        ierr=9999
        call xplasma_errmsg_append(sp,' ?? eqi_xtrz:  BLOAT failure.')
        call xplasma_errmsg_append(sp,'    Failed to cover (R,Z) rectangle.')
        go to 1000
     else
        zfacx0=(1.1d0)*zfac*zfacx0
        zfaclim=(1.1d0)*zfac*zfaclim
        deallocate(xiblo,rmc,ymc,rmcx,ymcx)
        go to 10
     endif
  endif

  !  extrapolate (rho,theta) space

  call xplasma_create_grid(sp,'__RHOX',xplasma_rhox_coord,xiblo(3:inbblo), &
       iaxis_rhox,ierr, label='extrapolated rho grid')

  call xplasma_create_grid(sp,'__THETAX',xplasma_thx_coord,th_xpl_fin, &
       iaxis_thetax,ierr, label='theta grid for extrapolated region')

  call eqi_xtrza(xiblo(1),inbblo,rmcx(1,0,1),ymcx(1,0,1),mimom, &
       inzp1,xiblo(3:inbblo),nrhox,th_xpl_fin,nthetax,th_arc_fin,rz_fin, &
       iaxis_rhox,iaxis_thetax,raxis,zaxis)

1000 continue
  deallocate(xiblo)

  deallocate(rmc)
  deallocate(ymc)

  deallocate(rmcx)
  deallocate(ymcx)

  deallocate(xi_bdy,th_bdy,zdlen_bdy,rz_bdy)
  deallocate(th_arc,th_arc_pkg,thspl)
  deallocate(th_arc_ntk,th_xpl_ntk,xi_bdy_ntk,rz_ntk)
  deallocate(th_arc_exp,th_xpl_exp)
  deallocate(th_arc_fin,th_xpl_fin,xi_bdy_fin,rz_fin)

  call xplasma_author_clear(sp,xplasma_xmhd,iertmp)

end subroutine eqi_xtrz

!------------------------------------------------------
subroutine eqi_xtrzf(symflag,intk,zrho,rmc,ymc,insurf,mimom,rz_ntk)

  use xplasma_definitions
  use eqi_rzbox_module

  !  get data and do fft to define a surface moments expansion

  implicit NONE

  integer intk                      ! no. of theta points
  logical symflag                   ! symmetry flag
  integer insurf                    ! #of surfaces in rmc,ymc
  real*8 zrho(insurf)               ! surface coordinate
  integer mimom                     ! #of moments in rmc,ymc
  real*8 rmc(insurf,0:mimom,2)      ! R moments (surf. i written)
  real*8 ymc(insurf,0:mimom,2)      ! Y moments (surf. i written)

  real*8 rz_ntk(intk,2)             ! (R,Z) contour at boundary

  !------------------------------------------------
  integer i,ith

  !  automatic:
  real*8 zrwk(intk),zzwk(intk),raxis,zaxis
  real*8 zrmc(3,0:mimom,2),zymc(3,0:mimom,2)

  !------------------------------------------------
  !  caller must assure that zrho(insurf) is at the boundary
  !  and zrho(insurf-1) is some distance in from the boundary

  raxis=rmc(1,0,1)
  zaxis=ymc(1,0,1)

  zrmc=0.0d0
  zymc=0.0d0

  zrwk=0.3d0*raxis+0.7d0*rz_ntk(1:intk,1)
  zzwk=0.7d0*zaxis+0.7d0*rz_ntk(1:intk,2)

  call r8tkbmmc(symflag,zrmc,zymc,mimom,mimom,zrwk,zzwk)
  do i=0,mimom
     rmc(insurf-1,i,1)=zrmc(3,i,1)
     rmc(insurf-1,i,2)=zrmc(3,i,2)
     ymc(insurf-1,i,1)=zymc(3,i,1)
     ymc(insurf-1,i,2)=zymc(3,i,2)
  enddo

  zrmc=0.0d0
  zymc=0.0d0

  zrwk=rz_ntk(1:intk,1)
  zzwk=rz_ntk(1:intk,2)

  call r8tkbmmc(symflag,zrmc,zymc,mimom,mimom,zrwk,zzwk)
  do i=0,mimom
     rmc(insurf,i,1)=zrmc(3,i,1)
     rmc(insurf,i,2)=zrmc(3,i,2)
     ymc(insurf,i,1)=zymc(3,i,1)
     ymc(insurf,i,2)=zymc(3,i,2)
  enddo

end subroutine eqi_xtrzf

!-----------------------------------------------------------------
subroutine eqi_xtrz_test(zrmin,zrmax,zzmin,zzmax, &
     xiblo,inbblo,rmcx,zmcx,mimom,zfac,ierr)

  !  test if extrapolation covers the R,Z cartesian grid

  use xplasma_definitions
  use eqi_rzbox_module

  implicit NONE
  !  input:
  real*8 zrmin,zrmax,zzmin,zzmax    ! range that must be covered
  integer inbblo                    ! no. of moments surfaces
  real*8 xiblo(inbblo)              ! surfaces where moments are given
  integer mimom                     ! no. of moments in Fourier Expansion
  real*8 rmcx(inbblo,0:mimom,2)     ! R moments
  real*8 zmcx(inbblo,0:mimom,2)     ! Z moments

  !  output:
  real*8 zfac                       ! estimate of further expansion needed
  integer ierr                      ! completion code, 0=OK

  !  if ierr=1, zfac is a multiplicative factor to apply to the range
  !  of the previous extrapolation, to get the right amount of total
  !  extrapolation

  !---------------------------------------
  integer ntest
  parameter (ntest=200)

  real*8 zr0,zz0
  real*8 zrtest(ntest),zztest(ntest),ztheta
  real*8 zsn(mimom),zcs(mimom),zRsum,zZsum
  real*8 zrminc,zrmaxc,zzminc,zzmaxc

  integer itop,ileft,ibot,iright,ith,itp,imi,im
  integer itest,ilim1,ilim2

  real*8 ztest1,ztest2,zxpmin
  real*8, parameter :: c2pi = 6.2831853071795862D+00
  real*8, parameter :: czero= 0.0d0
  real*8, parameter :: chalf= 0.5d0

  !---------------------------------------
  zr0=chalf*(zrmin+zrmax)
  zz0=chalf*(zzmin+zzmax)

  do ith=1,ntest
     ztheta=(ith-1)*c2pi/(ntest-1)
     call r8sincos(ztheta,mimom,zsn,zcs)

     imi=inbblo
     zRsum=rmcx(imi,0,1)
     zZsum=zmcx(imi,0,1)
     do im=1,mimom
        zRsum=zRsum+zcs(im)*rmcx(imi,im,1)+zsn(im)*rmcx(imi,im,2)
        zZsum=zZsum+zcs(im)*zmcx(imi,im,1)+zsn(im)*zmcx(imi,im,2)
     enddo

     zrtest(ith)=zRsum
     zztest(ith)=zZsum
  enddo

  !  find pts near corners

  itop=0
  ileft=0
  ibot=0
  iright=0

  do ith=1,ntest-1
     itp=ith+1
     if(((zrtest(itp)-zr0).gt.czero).and. &
          ((zztest(itp)-zz0).gt.czero)) then
        ztest1 = (zzmax-zz0)*(zrmax-zrtest(ith))- &
             (zrmax-zr0)*(zzmax-zztest(ith))
        ztest2 = (zzmax-zz0)*(zrmax-zrtest(itp))- &
             (zrmax-zr0)*(zzmax-zztest(itp))
        if(ztest1*ztest2.le.czero) then
           itop=itp
        endif
     else if(((zrtest(itp)-zr0).lt.czero).and. &
          ((zztest(itp)-zz0).gt.czero)) then
        ztest1 = (zzmax-zz0)*(zrmin-zrtest(ith))- &
             (zrmin-zr0)*(zzmax-zztest(ith))
        ztest2 = (zzmax-zz0)*(zrmin-zrtest(itp))- &
             (zrmin-zr0)*(zzmax-zztest(itp))
        if(ztest1*ztest2.le.czero) then
           ileft=itp
        endif
     else if(((zrtest(itp)-zr0).lt.czero).and. &
          ((zztest(itp)-zz0).lt.czero)) then
        ztest1 = (zzmin-zz0)*(zrmin-zrtest(ith))- &
             (zrmin-zr0)*(zzmin-zztest(ith))
        ztest2 = (zzmin-zz0)*(zrmin-zrtest(itp))- &
             (zrmin-zr0)*(zzmin-zztest(itp))
        if(ztest1*ztest2.le.czero) then
           ibot=itp
        endif
     else if(((zrtest(itp)-zr0).gt.czero).and. &
          ((zztest(itp)-zz0).lt.czero)) then
        ztest1 = (zzmin-zz0)*(zrmax-zrtest(ith))- &
             (zrmax-zr0)*(zzmin-zztest(ith))
        ztest2 = (zzmin-zz0)*(zrmax-zrtest(itp))- &
             (zrmax-zr0)*(zzmin-zztest(itp))
        if(ztest1*ztest2.le.czero) then
           iright=itp
        endif
     endif
  enddo

  !  sanity checks...

  if((itop.eq.0).or.(ibot.eq.0).or.(iright.eq.0).or.(ileft.eq.0)) then
     call xplasma_errmsg_append(sp,' ?? eqi_xtrz_test:  algorithm error 1!')
     ierr=9999
     zfac=1.5d0
     return
  endif

  itest=0
  if(ileft.gt.itop) itest=itest+1
  if(ibot.gt.ileft) itest=itest+1
  if(iright.gt.ibot) itest=itest+1
  if(itop.gt.iright) itest=itest+1
  if(itest.ne.3) then
     call xplasma_errmsg_append(sp,' ?? eqi_xtrz_test:  algorithm error 2!')
     ierr=9999
     zfac=1.5d0
     return
  endif

  zrmaxc=zrtest(iright)
  zrminc=zrtest(ileft)
  zzmaxc=zztest(itop)
  zzminc=zztest(ibot)

  !  min Z on top

  ilim1=itop
  ilim2=ileft
  if(ilim2.lt.ilim1) ilim2=ilim2+ntest
  do ith=ilim1,ilim2
     itp=ith
     if(itp.gt.ntest) itp=itp-ntest
     zzmaxc=min(zzmaxc,zztest(itp))
  enddo

  !  max R on left

  ilim1=ileft
  ilim2=ibot
  if(ilim2.lt.ilim1) ilim2=ilim2+ntest
  do ith=ilim1,ilim2
     itp=ith
     if(itp.gt.ntest) itp=itp-ntest
     zrminc=max(zrminc,zrtest(itp))
  enddo

  !  max Z at bottom

  ilim1=ibot
  ilim2=iright
  if(ilim2.lt.ilim1) ilim2=ilim2+ntest
  do ith=ilim1,ilim2
     itp=ith
     if(itp.gt.ntest) itp=itp-ntest
     zzminc=max(zzminc,zztest(itp))
  enddo

  !  min R on right

  ilim1=iright
  ilim2=itop
  if(ilim2.lt.ilim1) ilim2=ilim2+ntest
  do ith=ilim1,ilim2
     itp=ith
     if(itp.gt.ntest) itp=itp-ntest
     zrmaxc=min(zrmaxc,zrtest(itp))
  enddo

  ierr=0
  if((zrminc.le.zrmin).and.(zrmaxc.ge.zrmax).and. &
       (zzminc.le.zzmin).and.(zzmaxc.ge.zzmax)) then
     continue
  else
     ierr=1
     zxpmin=1.1d0
     zfac=max(zxpmin, &
          (zrmax-zr0)/(zrmaxc-zr0),(zzmax-zz0)/(zzmaxc-zz0), &
          (zrmin-zr0)/(zrminc-zr0),(zzmin-zz0)/(zzminc-zz0))
  endif

end subroutine eqi_xtrz_test

!-----------------------------------------------------------------
subroutine eqi_xtrza(xiblo,inbblo,rmcx,zmcx,mimom,ibdy, &
     rhox,nrhox,thetax,nthetax,th_fft,rz_xbdy,iaxis_rhox,iaxis_thetax, &
     raxis,zaxis)

  use xplasma_definitions
  use eqi_rzbox_module

  implicit NONE

  !  given the extrapolated surface data in moments form, convert to
  !  spline representation & store as xplasma profiles

  !  input:
  integer inbblo                    ! no. of moments surfaces
  real*8 xiblo(inbblo)              ! surfaces where moments are given
  integer mimom                     ! no. of moments in Fourier Expansion
  real*8 rmcx(inbblo,0:mimom,2)     ! R moments
  real*8 zmcx(inbblo,0:mimom,2)     ! Z moments
  integer ibdy                      ! index of plasma bdy surf. in xiblo

  integer :: nrhox                  ! extrap. rho grid
  real*8 :: rhox(nrhox)
  integer :: nthetax                ! extrap. theta grid
  real*8 :: thetax(nthetax)
  real*8 :: th_fft(nthetax)
  real*8 :: rz_xbdy(nthetax,2)      ! R & Z bdy values @ rho=1

  integer :: iaxis_rhox,iaxis_thetax  ! extrapolation grid Ids

  real*8, intent(in) :: raxis, zaxis

  !--------------------------------------
  integer :: ispline,iertmp,itheta,irho,inum,id,ishift
  real*8, dimension(:,:), allocatable :: Rx,Zx,distx,atanx
  real*8, dimension(:), allocatable :: atanx0
  real*8 :: zth,zthmin,zatan,zdra,zdza
  character*80 atanx0_label
  real*8,parameter :: C2PI = 6.2831853071795862D+00
  !--------------------------------------
  !  it will be important to use the same fit for the extrapolated R & Z
  !  profiles as was used for the core profiles...

  call xplasma_global_info(sp,iertmp, rzOrder=ispline)

  !------------------
  !  1.  set boundary values to match core region

  allocate(rx(nthetax,nrhox),zx(nthetax,nrhox))
  allocate(distx(nthetax,nrhox),atanx(nthetax,nrhox),atanx0(nrhox))

  Rx(1:nthetax,1) = rz_xbdy(1:nthetax,1)
  Zx(1:nthetax,1) = rz_xbdy(1:nthetax,2)

  !  2.  find the R & Z sequence for each theta, do a 1d spline on this,
  !     and insert results in 2d spline data arrays

  inum=inbblo-ibdy

  do itheta=1,nthetax
     if(itheta.eq.nthetax) then
        do irho=2,nrhox
           !  periodic condition:
           Rx(nthetax,irho)=Rx(1,irho)
           Zx(nthetax,irho)=Zx(1,irho)
        enddo
     else
        zth=th_fft(itheta)
        call eqi_xtrz_fill(zth,Rx(itheta,2:nrhox),Zx(itheta,2:nrhox),inum, &
             rmcx,zmcx,inbblo,mimom,ibdy)
     endif
  enddo

  !  3.  setup the profile splines...

  call xplasma_create_2dprof(sp,'__R_extrap', &
       iaxis_thetax,iaxis_rhox,Rx,id,iertmp, &
       ispline=ispline,label='R(rho,theta) extrapolation',units='m')

  call xplasma_create_2dprof(sp,'__Z_extrap', &
       iaxis_thetax,iaxis_rhox,Zx,id,iertmp, &
       ispline=ispline,label='Z(rho,theta) extrapolation',units='m')

  do irho=1,nrhox
     do itheta=1,nthetax-1
        zdra=Rx(itheta,irho)-raxis
        zdza=Zx(itheta,irho)-zaxis
        distx(itheta,irho)=sqrt(zdra*zdra+zdza*zdza)
        atanx(itheta,irho)=atan2(zdza,zdra)
     enddo
     distx(nthetax,irho)=distx(1,irho)
     !  now make atanx(theta,rho) be relative to 1st theta atan.
     zatan=atanx(1,irho)
     atanx0(irho)=zatan
     atanx(1,irho)=0.0d0
     atanx(nthetax,irho)=C2PI
     ishift=0
     do itheta=2,nthetax-1
        if(ishift.eq.0) then
           if(atanx(itheta,irho).lt.zatan) then
              ishift=1
           else
              zatan=atanx(itheta,irho)
           endif
        endif
        atanx(itheta,irho)=atanx(itheta,irho)-atanx0(irho)+ishift*C2PI
     enddo
  enddo

  call xplasma_create_2dprof(sp,'__distx_axis',iaxis_thetax,iaxis_rhox, &
       distx,id,iertmp, &
       label='distance to axis',units='m')

  call xplasma_create_2dprof(sp,'__atanx_axis',iaxis_thetax,iaxis_rhox, &
       atanx,id,iertmp, &
       label='adjusted arctan((Zx-Zaxis)/(Rx-Raxis))',units='rad')

  call xplasma_grid_info(sp,iaxis_thetax,iertmp, xmin=zthmin)
 
  atanx0_label=' '
  write(atanx0_label, &
       '(" arctan((Zx-Zaxis)/(Rx-Raxis)) vs. rho @theta=",1pd11.4)') zthmin

  call xplasma_create_prof(sp,'__atanx0_axis',iaxis_rhox, &
       atanx0,id,iertmp, &
       label=trim(atanx0_label),units='rad')
 
  deallocate(Rx,Zx,distx,atanx,atanx0)

end subroutine eqi_xtrza

!------------------------------------------------------------------
subroutine eqi_xtrz_fill(ztheta,Rdata,Zdata,inum, &
     rmcx,zmcx,inbblo,mimom,ibdy)

  !  setup 1d spline vs. rho (beyond plasma bdy) along a fixed
  !  theta line.

  !  the moments input defines an axial point, a non-singular surface
  !  just inside the boundary, the boundary surface, and then a set of
  !  extrapolated surfaces...

  implicit NONE

  real*8 ztheta                     ! theta angle (w.r.t. Fourier exp.)

  integer inum                      ! number of pts
  real*8 :: Rdata(inum),Zdata(inum) ! data points to be computed

  integer inbblo                    ! #of surfaces in moments array
  integer mimom                     ! #of moments in moments array

  integer ibdy                      ! bdy surface index in rmcx,zmcx
  real*8 rmcx(inbblo,0:mimom,2)     ! R moments
  real*8 zmcx(inbblo,0:mimom,2)     ! Z moments

  !--------------------------------------

  integer irho,isurf,im
  real*8 zsn(mimom),zcs(mimom)      ! sincos ladder (automatic storage)
  real*8 zRsum,zZsum

  !--------------------------------------

  !  the 1st data points of R1spl, Z1spl match the interior R,Z set and
  !  are presumed to already have been set ***

  call r8sincos(ztheta,mimom,zsn,zcs)

  do isurf=ibdy+1,inbblo
     irho=isurf-ibdy

     zRsum=rmcx(isurf,0,1)
     zZsum=zmcx(isurf,0,1)
     do im=1,mimom
        zRsum=zRsum+zcs(im)*rmcx(isurf,im,1)+zsn(im)*rmcx(isurf,im,2)
        zZsum=zZsum+zcs(im)*zmcx(isurf,im,1)+zsn(im)*zmcx(isurf,im,2)
     enddo
     Rdata(irho)=zRsum
     Zdata(irho)=zZsum
  enddo

end subroutine eqi_xtrz_fill

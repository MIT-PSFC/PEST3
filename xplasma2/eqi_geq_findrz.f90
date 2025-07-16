!! subroutine eqi_geq_findrz(ns,lmask,x,dif,idum,psi,idum2,rzvec,id2)
!!   (original routine removed DMC Mar 2008)

subroutine eqi_find_curv(nsrch, lmask, psivec, result, id1, zinput, id2i, &
     zoutput, id2o)

  use eqi_geq_mod
  implicit NONE
  integer, parameter :: r8=selected_real_kind(12,100)

  ! subroutine called by root finder "zriddery" to find psi-surface with
  ! desired matching curvature.  Although formally a collection of searches
  ! can be done, in practice this routine only supports 1 search at a time.

  ! this routine calls eqi_find_curvature which calls zridderx -> a root
  ! finder within a root finder; 2d root finding operation.

  integer, intent(in) :: nsrch   ! number of searches, must = 1
  logical, intent(in) :: lmask(nsrch) ! mask searches where soln is known
 
  real*8, intent(in) :: psivec(nsrch)    ! psi location to try next
 
  real*8, intent(out) :: result(nsrch)   ! curvature difference at surfaces
 
  integer, intent(in) :: id1,id2i,id2o   ! auxilliary array dimensions
 
  real*8, intent(in) :: zinput(id1,id2i) ! target curvature value here
  real*8, intent(out) :: zoutput(id1,id2o)  ! output: 0=OK, 1= arg error
 
  !----------------------------------
  real*8 :: psi_in,xmin,curv_out,curv_sought
  !----------------------------------

  result = 0.0_r8
  if(nsrch.ne.1) then
     zoutput=1.0_r8
     return
  else
     zoutput=0.0_r8
  endif

  ! find Psi value on all rays (eqi_geq_mod: (r0,0) -> {radapt(:),zadapt(:)})
  ! compute normalized curvature

  psi_in = psivec(1)
  curv_sought = zinput(1,1)
  xmin = 0.4_r8  ! expect to look out towards the edge...

  !  use curvad
  call eqi_find_curvature(nthad,psi_in,xmin,psep,rzwk,curvec,curv_out)

  result(1) = curv_out-curv_sought

end subroutine eqi_find_curv

subroutine eqi_find_curvature(ns,psi,xmin,psi_range,rzvec,curvature,mincurv)

  ! find curvature properties of surface at a target psi value, as
  ! measured from a set of rays {(r0,z0)->(radapt(:),zadapt(:))}
  !   (see declarations in eqi_geq_mod

  use eqi_geq_mod
  implicit NONE
  integer, parameter :: r8=selected_real_kind(12,100)

  integer, intent(in) :: ns  ! number of rays; should be equal to nthad

  real*8, intent(in) :: psi  ! target psi value 

  real*8, intent(in) :: xmin ! minimum normalized distance from (r0,z0) at
  !  which to start search (in range [0:1]; could be 0; less than 1).
  !  0 corresponds to (r0,z0); 1 corresponds to r/zadapt.

  real*8, intent(in) :: psi_range  ! Psi(max)-Psi(min) in domain sought over

  real*8, intent(out) :: rzvec(ns,2)   ! (R,Z) values of psi surface found
  real*8, intent(out) :: curvature(ns) ! curvature along surface
  real*8, intent(out) :: mincurv       ! extrapolated minimum curvature

  !---------------------------------
  real*8 :: psivec(ns),vecmin(ns),vec1(ns),xvec(ns)
  real*8 :: tol
  logical :: lmask(ns)
  external :: eqi_geq_findrzs
  integer :: imin,imax,ichk,ii,ifail
  real*8 :: rmin,rmax,circen(2),radius,radmax,c1,c2,c3,aa,bb,xx
  real*8 :: delx12,delx23,r1,r2,r3
  !---------------------------------

  if(ns.ne.nthad) then
     call errmsg_exit('?internal call error, eqi_find_curvature.')
  endif

  psivec=psi
  lmask = .FALSE.
  lmask(ns)=.TRUE.  ! first & last pts match; don't recalculate
  vecmin = xmin
  vec1 = 1.0_r8

  tol = 1.0e-12_r8

  call zridderx(ns, lmask, vecmin, vec1, tol, tol*psi_range, &
       eqi_geq_findrzs, xvec, ifail, ns, psivec, 1, rzvec, 2)

  if(ifail.ne.0) then
     call errmsg_exit(' ?eqi_find_curvature: unexpected zridderx error.')
  endif

  rzvec(ns,1:2)=rzvec(1,1:2)
  imin=1
  imax=1
  rmin=rzvec(1,1)
  rmax=rmin
  do ii=2,ns-1
     if(rzvec(ii,1).lt.rmin) then
        rmin=rzvec(ii,1)
        imin=ii
     endif
     if(rzvec(ii,1).gt.rmax) then
        rmax=rzvec(ii,1)
        imax=ii
     endif
  enddo

  ! quadratic fit
  if(imin.eq.1) then
     r1=rzvec(ns-1,1)
     r2=rmin
     r3=rzvec(2,1)
     delx12=thadapt(ns,isw)-thadapt(ns-1,isw)
     delx23=thadapt(2,isw)-thadapt(1,isw)
  else
     r1=rzvec(imin-1,1)
     r2=rmin
     r3=rzvec(imin+1,1)
     delx12=thadapt(imin,isw)-thadapt(imin-1,isw)
     delx23=thadapt(imin+1,isw)-thadapt(imin,isw)
  endif

  call quad_minmax(delx12,delx23,r1,r2,r3,rmin)

  if(imax.eq.1) then
     r1=rzvec(ns-1,1)
     r2=rmax
     r3=rzvec(2,1)
     delx12=thadapt(ns,isw)-thadapt(ns-1,isw)
     delx23=thadapt(2,isw)-thadapt(1,isw)
  else
     r1=rzvec(imax-1,1)
     r2=rmax
     r3=rzvec(imax+1,1)
     delx12=thadapt(imax,isw)-thadapt(imax-1,isw)
     delx23=thadapt(imax+1,isw)-thadapt(imax,isw)
  endif

  call quad_minmax(delx12,delx23,r1,r2,r3,rmax)

  lnorm = (rmax-rmin)/2
  radmax=10*lnorm

  !  end pts first

  imin=1
  call r8_trirad(rzvec(ns,:),rzvec(1,:),rzvec(2,:),circen,radius,ichk)
  if(ichk.eq.0) then
     curvature(1)=min(radmax,radius)
     mincurv=curvature(1)
  else
     mincurv=radmax
     curvature(1)=radmax
  endif
  curvature(ns)=curvature(1)

  ! search the rest

  do ii=2,ns-1
     call r8_trirad(rzvec(ii-1,:),rzvec(ii,:),rzvec(ii+1,:),circen,radius,ichk)
     if(ichk.eq.0) then
        if(radius.lt.mincurv) then
           imin=ii
           mincurv=radius
        endif
        curvature(ii)=min(radmax,radius)
     else
        curvature(ii)=radmax
     endif
  enddo

  ! normalize the curvature (mincurv not normalized; recomputed below)

  curvature = curvature/lnorm

  ! quadratic fit
  if(imin.eq.1) then
     c1=curvature(ns-1)
     c2=curvature(1)
     c3=curvature(2)
     delx12=thadapt(ns,isw)-thadapt(ns-1,isw)
     delx23=thadapt(2,isw)-thadapt(1,isw)
  else
     c1=curvature(imin-1)
     c2=curvature(imin)
     c3=curvature(imin+1)
     delx12=thadapt(imin,isw)-thadapt(imin-1,isw)
     delx23=thadapt(imin+1,isw)-thadapt(imin,isw)
  endif

  call quad_minmax(delx12,delx23,c1,c2,c3,mincurv)
  
end subroutine eqi_find_curvature

subroutine eqi_geq_findrzs(ns,lmask,x,dif,idum,psi,idum2,rzvec,id2)
 
  use eqi_geq_mod
  implicit NONE
  integer, parameter :: r8=selected_real_kind(12,100)
 
  !
  ! rootfinder (zridderx) subroutine -- difference btw target and actual
  !    abs(psi) for various points along a collection of line segments
  !    (r0,z0) to {radapt(1:ns,isw),zadapt(1:ns,isw)}.
  !
  ! see declarations in eqi_geq_mod
  !
 
  !------------------------------ passed
 
  integer, intent(in) :: ns        ! no. of segments
  logical, intent(in) :: lmask(ns) ! mask segments where target psi was reached
 
  real(r8), intent(in) :: x(ns)    ! parameter, 0->(r0,z0), 1->(r/zseek*(:))
 
  real(r8), intent(out) :: dif(ns) ! psi difference at surfaces
 
  integer, intent(in) :: idum,idum2,id2
 
  real(r8), intent(in) :: psi(ns)  ! target psi values
 
  real(r8), intent(out) :: rzvec(ns,id2)  ! (R,Z) values
 
  !------------------------------ local
 
  integer jvec,j,ier
  real(r8) rtmp(ns),ztmp(ns),psitmp(ns)
  real(r8), parameter :: ZERO=0.0_r8
 
  !------------------------------ executable code
 
  jvec=0
  do j=1,ns
     if(.not.lmask(j)) then
        jvec=jvec+1
        rtmp(jvec)=r0+x(j)*(radapt(j,isw)-r0)
        ztmp(jvec)=z0+x(j)*(zadapt(j,isw)-z0)
     endif
  enddo
 
  !
  !  R0,Z0 and r/zseek are inside the mapped region
  !  0 < x(j) < 1  --> rtmp,ztmp are inside mapped region
  !  so spline call should be safe
  !
 
  call ezspline_interp(pspl, jvec, rtmp, ztmp, psitmp, ier)
  if(ier.ne.0) then
     call errmsg_exit(' ??eqi_geq_findrzs:  unexpected spline error.')
  endif
 
  jvec=0
  do j=1,ns
     if(.not.lmask(j)) then
        jvec=jvec+1
        dif(j)=abs(psi(j))-psitmp(jvec)
        rzvec(j,1)=rtmp(jvec)
        rzvec(j,2)=ztmp(jvec)
     endif
  enddo
 
  return
end subroutine eqi_geq_findrzs

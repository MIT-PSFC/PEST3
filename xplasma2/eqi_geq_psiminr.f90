subroutine eqi_geq_psiminr(ns,lmask,rr,dpsidr,idum,zerr,idum2,dpsidz,id2)
 
  use eqi_geq_mod
  implicit NONE
  integer, parameter :: r8=selected_real_kind(12,100)
 
  !
  ! rootfinder (zridderx) subroutine -- 
  !    find minimum (dpsi/dR = 0) point along line of fixed Z=z0.
  !
  ! z0 and psi spline, are all in eqi_geq_mod
  !
  !------------------------------ passed
 
  integer, intent(in) :: ns        ! no. of surfaces
  logical, intent(in) :: lmask(ns) ! mask surfaces where soln is known
 
  real(r8), intent(in) :: rr(ns)   ! R value to try now
 
  real(r8), intent(out) :: dpsidr(ns) ! d(psi)/dR returned
 
  integer, intent(in) :: idum,idum2,id2
 
  real(r8), intent(out) :: zerr(ns)   ! =0; set =1 if error occurs
 
  real(r8), intent(out) :: dpsidz(ns) ! d(psi)/dZ returned
 
  !------------------------------ local
 
  integer jvec,j,ier
  real(r8) rtmp(ns),ztmp(ns),psitmp(ns),zgrad(ns,2)
 
  !------------------------------ executable code
 
  jvec=0
  do j=1,ns
     if(.not.lmask(j)) then
        jvec=jvec+1
        rtmp(jvec)=rr(j)
        ztmp(jvec)=z0
     endif
  enddo
 
  !
  !  R0,Z0 and RPP,ZPP are inside the mapped region
  !  0 < x(j) < 1  --> rtmp,ztmp are inside mapped region
  !  so spline call should be safe
  !
 
  call ezspline_gradient(pspl, jvec, rtmp, ztmp, zgrad, ier)
  if(ier.ne.0) then
     zgrad=0
     zerr=1
  else
     zerr=0
  endif
 
  jvec=0
  do j=1,ns
     if(.not.lmask(j)) then
        jvec=jvec+1
        dpsidr(j)=zgrad(jvec,1)
        dpsidz(j)=zgrad(jvec,2)
     endif
  enddo
 
  return
end subroutine eqi_geq_psiminr

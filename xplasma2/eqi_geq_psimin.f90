subroutine eqi_geq_psimin(ns,lmask,zz,dpsidz,idum,rzin,id1,rout,id2)
 
  use eqi_geq_mod
  implicit NONE
  integer, parameter :: r8=selected_real_kind(12,100)
 
  !
  ! rootfinder (zriddery) subroutine -- 
  !    find minimum Psi location
  !    by finding d(psi)/dZ=0 where the derivative is taken
  !    at a location where d(psi)/dR=0 along some line Z=const=(zz)
  !
  ! z0 and psi spline, are all in eqi_geq_mod
  !
  !------------------------------ passed
 
  integer, intent(in) :: ns        ! no. of surfaces
  logical, intent(in) :: lmask(ns) ! mask surfaces where soln is known
 
  real(r8), intent(in) :: zz(ns)   ! Z value to try now
 
  real(r8), intent(out) :: dpsidz(ns) ! d(psi)/dZ returned
 
  integer, intent(in) :: idum,id1,id2
 
  real(r8), intent(in) :: rzin(ns,id1) ! range of R variation to try; tolerance
 
  real(r8), intent(out) :: rout(ns,id2) ! {R,err} returned
 
  !------------------------------ local
 
  real(r8) rr1,rr2,rr_out,tol,tolf,zerr(1),dpsidy(1)
  integer ifail,jfail
  external eqi_geq_psiminr

  !------------------------------ executable code
 
  if(ns.ne.1) then
     rout(:,2)=1
     rout(:,1)=0
     dpsidz=0
     return
  endif

  z0=zz(1)  ! now will look for d(Psi)/dR=0 along Z=z0...
  rr1=rzin(1,1)-rzin(1,2)  ! min R to look at
  rr2=rzin(1,1)+rzin(1,2)  ! max R to look at
  tol=rzin(1,3)
  tolf=rzin(1,4)

  !  dpsidz will also be evaluated...
  
  call zridderx(ns,lmask,rr1,rr2,tol,tolf,eqi_geq_psiminr,rr_out,ifail, &
       1,zerr,1,dpsidy,1)
  jfail=abs(zerr(1))+0.5_R8

  if((ifail.ne.0).or.(jfail.ne.0)) then
     rout(:,2)=1
     rout(:,1)=0
     dpsidz=0
  else
     rout(:,2)=0
     rout(:,1)=rr_out
     dpsidz=dpsidy(1)
  endif
  
  return
end subroutine eqi_geq_psimin

subroutine eqi_tau_bounce(zrho,indx,itrap,inum,zx,ztol,zcut, &
     zbmin,zbmax,zvtaub,ierr)
 
  !  evaluate zero banana width approximation to v*tau_bounce
  !
  !  this is the integral over the limits of the orbit of
  !     dtheta*(dlp/dtheta)*(B/Bp)*(1/sqrt(1-B/Brefl))
  !
  !  where Brefl is a trapped orbit's reflection point (vpll=0,v=vperp)
  !  for passing orbits (Brefl > Bmax) the integrand is over the
  !  entire surface contour; for trapped orbits the integrand is
  !  twice the value from theta_lower to theta_upper where these
  !  limits need to be found by a root finder so that
  !     B(theta_lower)=B(theta_upper)=Brefl
  !
  !  the integrand formula follows from the assumptions that v
  !  and mu=vperp**2/B are constants of motion
 
  !  axisymmetry is assumed...

  !  This is an internal private subroutine; error checking done in
  !  xplasma_ctran.f90:xplasma_vtaub contained subroutine.

  !------------------
 
  use xplasma_definitions
  use eqi_rzbox_module
  use eqi_vtaub_mod

  implicit NONE
 
  real*8, intent(in) :: zrho   ! surface on which to calculate
  !  note -- must not be magnetic axis; non-singular surface required.

  integer, intent(in) :: indx  ! vtaub integrand data element indx
  !  (could be a thread number for parallel implementations)
  !  (see eqi_vtaub_mod)

  integer, intent(in) :: itrap ! no. of eval points of trapped orbits
  integer, intent(in) :: inum  ! no. of eval points
  real*8, intent(in) :: zx(inum) ! evaluation points:
 
  !   zx(j) = (Brefl(j)-Bmin)/(Bmax-Bmin)
  ! where Bmin and Bmax, the min and mod(B) on the zrho surface, will
  ! be found; zx > 0 is required; zx(j+1) > zx(j) required for all j < inum.
 
  real*8, intent(in) :: ztol  ! relative error tolerance for integral
  real*8, intent(in) :: zcut  ! minimum allowed value for (1-B/Brefl)
 
  !---
 
  real*8, intent(out) :: zbmin,zbmax   ! Bmin and Bmax on zrho surface (T)
  real*8, intent(out) :: zvtaub(inum)  ! integral values returned (m)
  !  divide by ion velocity (m/sec) to get bounce time in seconds.
 
  integer, intent(out) :: ierr  ! completion code, 0=OK
 
  !-------------------------------------------
  real*8, external :: f_eq_tau_bounce  ! integrand function
  integer i,ieval,iwarn,ilast,id_R,id_th,inth,ith
  integer :: id_Z,id_Bmod,id_BR,id_BZ
  integer, parameter :: ilimit=1000, ilimit4 = 4*ilimit
  real*8, dimension(:), allocatable :: zbvals,zthbr1,zthbr2,zwork
  integer, dimension(:), allocatable :: iwork
  real*8 zthbmin,zthbmax,zaux(4),zabserr
  character*120 zmsg

  real*8, dimension(:), allocatable :: zwkrho,zmodb,zbpol,zdldth
  real*8, dimension(:,:), allocatable :: zwkfth
  integer :: ilist(5),iderivs(5)

  real*8, parameter :: CPI = 3.1415926535897931D+00
  real*8, parameter :: ZERO = 0.0D0
  !-------------------------------------------
 
  zbmin=0
  zbmax=0
  zvtaub=0
  ierr=0

  !  basic grid info
  call xplasma_common_ids(sp,ierr, &
       id_Bmod=id_Bmod,id_BR=id_BR,id_BZ=id_BZ,id_R=id_R,id_Z=id_Z)
  if(ierr.ne.0) return

  if(min(id_Bmod,id_BR,id_BZ).eq.0) then
     ierr=9999
     call xplasma_errmsg_append(sp,' ?eqi_tau_bounce: field data missing.')
     return
  endif

  ilist(1)=id_Bmod; iderivs(1)=0
  ilist(2)=id_BR;   iderivs(2)=0
  ilist(3)=id_BZ;   iderivs(3)=0
  ilist(4)=id_R;    iderivs(4)=1
  ilist(5)=id_Z;    iderivs(5)=1

  call xplasma_prof_info(sp,id_R,ierr, gridId1=id_th)
  if(ierr.ne.0) return

  call xplasma_grid_info(sp,id_th,ierr, size=inth)
  if(ierr.ne.0) return
 
  !  allocate integration data element

  call elem_alloc(indx,inth,ierr)  ! eqi_vtaub_mod 
  if(ierr.ne.0) then
     call errstr(ierr,zmsg)
     ierr=9999
     call xplasma_errmsg_append(sp,zmsg)
     return
  endif

  !  get data for integration

  allocate(zwkrho(inth),zwkfth(inth,5)); zwkrho = zrho
  allocate(zmodb(inth),zbpol(inth),zdldth(inth))

  call xplasma_eval_prof(sp,ilist, &
       xplasma_theta_coord,vta(indx)%th, &
       xplasma_rho_coord,zwkrho, &
       zwkfth,ierr, &
       ideriv1s=iderivs)

  if(ierr.ne.0) return

  do ith=1,inth
     zmodb(ith)=abs(zwkfth(ith,1))
     zbpol(ith)=sqrt(zwkfth(ith,2)**2+zwkfth(ith,3)**2)
     zdldth(ith)=sqrt(zwkfth(ith,4)**2+zwkfth(ith,5)**2)
  enddo

  call elem_init(indx,zmodb,zbpol,zdldth,ierr)

  deallocate(zwkrho,zwkfth,zmodb,zbpol,zdldth)

  if(ierr.ne.0) then
     call errstr(ierr,zmsg)
     ierr=9999
     call xplasma_errmsg_append(sp,zmsg)
     return
  endif

  !-----------------------------------------------------
  !  OK: field and dl/dth splines ready...

  !  get integration limits for trapped population
 
  allocate(zbvals(inum),zthbr1(inum),zthbr2(inum))
  zthbr1=-CPI
  zthbr2=CPI
 
  if(itrap.eq.0) then
     call eqi_bbox(zrho,zbmin,zthbmin,zbmax,zthbmax,ierr)
  else
     call eqi_bbvec(zrho,itrap,zx(1:itrap), &
          zbmin,zbmax,zbvals(1:itrap),zthbr1(1:itrap),zthbr2(1:itrap),ierr)
  endif
  if(ierr.ne.0) go to 999
 
  do i=itrap+1,inum
     zbvals(i)=zbmin+zx(i)*(zbmax-zbmin)
  enddo
 
  !  calculate integrands...
 
  allocate(iwork(ilimit),zwork(ilimit4))
 
  zaux(2)=zcut
  zaux(3)=zrho
  zaux(4)=indx+0.1d0

  do i=1,inum
     zaux(1)=zbvals(i)
     call tr_r8qags(f_eq_tau_bounce,zaux,zthbr1(i),zthbr2(i),ZERO,ztol, &
          zvtaub(i),zabserr,ieval,iwarn,ilimit,ilimit4,ilast,iwork,zwork)
  enddo
 
  zvtaub(1:itrap)=2*zvtaub(1:itrap)  ! trapped orbit goes back on self
 
  deallocate(iwork,zwork)
 
999 continue
  deallocate(zthbr1,zthbr2,zbvals)
 
end subroutine eqi_tau_bounce
 
!------------------------------------------------------------
real*8 function f_eq_tau_bounce(zchi,aux)
 
  !  integrand routine

  use eqi_vtaub_mod
  implicit NONE

  real*8, intent(in) :: zchi  ! poloidal angle
  real*8 aux(*)               ! auxilliary data
 
  !  aux(1) = Brefl
  !  aux(2) = cutoff for (1-B/Brefl)
  !  aux(3) = rho
  !  aux(4) = data element index
 
  real*8 :: zrho,zcut,zbrefl,zbp,zb,zdldth,zvpv,ztheta,zdth
  real*8, parameter :: ONE = 1.0d0
  integer :: indx,ith,inth
  !--------------------
 
  zrho=aux(3)
  zcut=aux(2)
  zbrefl=aux(1)
  indx=aux(4)

  !  bring argument in range; it is assumed to be nearly in range...

  ztheta=zchi

  if(ztheta.lt.-ZPI) then
     do
        ztheta = ztheta + 2*ZPI
        if(ztheta.ge.-ZPI) exit
     enddo

  else if(ztheta.gt.ZPI) then
     do
        ztheta = ztheta - 2*ZPI
        if(ztheta.le.ZPI) exit
     enddo

  endif

  !  find the zone & evaluate splines

  inth = vta(indx)%nth
  ith = 1 + (ztheta+ZPI)*(inth-1)/(2*ZPI)
  ith = max(1,min((inth-1),ith))
  zdth = ztheta - vta(indx)%th(ith)

  zb  = vta(indx)%fmodb(1,ith) + zdth*(vta(indx)%fmodb(2,ith) + &
       zdth*(vta(indx)%fmodb(3,ith) + zdth*vta(indx)%fmodb(4,ith)))

  zbp = vta(indx)%fbpol(1,ith) + zdth*(vta(indx)%fbpol(2,ith) + &
       zdth*(vta(indx)%fbpol(3,ith) + zdth*vta(indx)%fbpol(4,ith)))

  zdldth = vta(indx)%fdldth(1,ith) + zdth*(vta(indx)%fdldth(2,ith) + &
       zdth*(vta(indx)%fdldth(3,ith) + zdth*vta(indx)%fdldth(4,ith)))

  zvpv = sqrt(max(zcut,(ONE-zb/zbrefl)))  ! vpll/v = sqrt(1-B/Brefl)

  f_eq_tau_bounce = zdldth*zb/(zbp*zvpv)

end function f_eq_tau_bounce

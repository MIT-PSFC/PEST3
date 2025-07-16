subroutine pstCompMetric(nt1, ns, the, psi, ier)
  !
  ! compute metric quantities
  !
  ! pletzer@pppl.gov Mon May  1 15:24:15 EDT 2000
  !
  USE pstcom
  USE newmet
  USE mtrik1

  use i2mex_mod

  implicit none
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)

  integer, intent(in) :: nt1, ns
  real(r8), intent(in) :: the(*), psi(*)
  integer, intent(out) :: ier

  real(r8), dimension(:,:), allocatable :: zdum, zdum2, zdum3
  real(r8), dimension(:), allocatable :: zf
  integer i, iok, j, nt, jextra
  real(i2mex_r8) signj


!!$  include 'CUCCCC.inc'

  ier = 0
  nt = nt1 - 1

  signj = 2._i2mex_r8*(i2mex_o%isThetaClockwise-0.5_i2mex_r8)

  allocate(zdum(nt1, ns), stat=iok)
  allocate(zdum2(nt1, ns), stat=iok)
  allocate(zdum3(nt1, ns), stat=iok)
  allocate(zf(nt1 ), stat=iok)

  call i2mex_getP(ns, psi, pa, ier)
  call i2mex_error(ier)
  call i2mex_getPP(ns, psi, ppa, ier)
  call i2mex_error(ier)
  call i2mex_getQ(ns, psi, qa, ier)
  call i2mex_error(ier)
  call i2mex_getQP(ns, psi, qpa, ier)
  call i2mex_error(ier)
  call i2mex_getQPP(ns, psi, qppa, ier)
  call i2mex_error(ier)
  call i2mex_getG(ns, psi, ga, ier)
  call i2mex_error(ier)
  call i2mex_getGP(ns, psi, gpa, ier)
  call i2mex_error(ier)
  r = 1.0_r8 ! R * Bext 
  fa = r*ga/qa
  fpa = (r*gpa - qpa*fa)/qa


  call i2mex_getX(nt1, ns, the, psi, zdum, ier)
  call i2mex_error(ier)
  xa(1:nt1, 1:ns) = zdum(1:nt1, 1:ns)
  xa(nt1+1:nt1+4,:) =zdum(2:5,:)
  call i2mex_getZ(nt1, ns, the, psi, zdum, ier)
  call i2mex_error(ier)
  za(1:nt1, 1:ns) = zdum(1:nt1, 1:ns)
  za(nt1+1:nt1+4,:) =zdum(2:5,:)

  call i2mex_getGradX(nt1, ns, the, psi, zdum, zdum2, ier)
  call i2mex_error(ier)
  xdth(1:nt1, 1:ns) = zdum(1:nt1, 1:ns)
  xdps(1:nt1, 1:ns) = zdum2(1:nt1, 1:ns)
  xdth(nt1+1:nt1+4,:) =zdum(2:5,:)
  xdps(nt1+1:nt1+4,:) =zdum2(2:5,:)

  call i2mex_getGradZ(nt1, ns, the, psi, zdum, zdum2, ier)
  call i2mex_error(ier)
  zdth(1:nt1, 1:ns) = zdum(1:nt1, 1:ns)
  zdps(1:nt1, 1:ns) = zdum2(1:nt1, 1:ns)
  zdth(nt1+1:nt1+4,:) =zdum(2:5,:)
  zdps(nt1+1:nt1+4,:) =zdum2(2:5,:)

  ! ap: fixed indexing bug 6 Dec 2006
  xinf(1:mth2) = xa(3:nt1+3, ns) !!xa(3:nt1-1+3, ns)
  zinf(1:mth2) = za(3:nt1+3, ns) !!za(3:nt1-1+3, ns)

  xsq = xa**2

  call i2mex_getJ(nt1, ns, the, psi, zdum, ier)
  call i2mex_error(ier)
  xjacob(1:nt1, 1:ns) = zdum(1:nt1, 1:ns)
  xjacob(nt1+1:nt1+4,:) = zdum(2:5,:)

  call i2mex_getJp(nt1, ns, the, psi, zdum, ier)
  call i2mex_error(ier)
  xjprym(1:nt1, 1:ns) = zdum(1:nt1, 1:ns)
  xjprym(nt1+1:nt1+4,:) = zdum(2:5,:)

  call i2mex_getMetric(nt1, ns, the, psi, zdum3, zdum, zdum2, ier)
  call i2mex_error(ier)
  grpssq(1:nt1, 1:ns) = zdum2(1:nt1, 1:ns)
  grpssq(nt1+1:nt1+4,:) = zdum2(2:5,:)
  grpsth(1:nt1, 1:ns) = zdum(1:nt1, 1:ns)
  grpsth(nt1+1:nt1+4,:) = zdum(2:5,:)
  
  xsqdps = 2._r8 * xa * xdps


  !
  ! compute delta = theta_p - theta
  !
  call i2mex_getDelta(nt1, ns, the, psi, zdum, ier)
  delta(1:nt1, 1:ns) = zdum(1:nt1, 1:ns)
  delta(nt1+1:nt1+4,:) = zdum(2:5,:)
  
!!$  do i=1, ns
!!$     zf(1:nt1) = ga(i)*xjacob(1:nt1,i)/ &
!!$          & (qa(i)* xsq(1:nt1,i))
!!$     call i2mex_integratePeriodic1d(nt1, the, zf, nt1, the, delta(1,i))
!!$     delta(1:nt1,i) = delta(1:nt1,i) - the(1:nt1)
!!$     delta(nt1+1:nt1+4,i) = delta(2:5,i)
!!$  enddo
!!$    ! extrapolate to axis
!!$
!!$    do i=1, nt1+4
!!$       delta(i,1) = FCCCC0(delta(i,2), delta(i,3), delta(i,4), delta(i,5), &
!!$            &                 psi(2),  psi(3),  psi(4),  psi(5), &
!!$            & psi(1))
!!$    enddo

  !
  ! compute d(q*delta)/d psi
  !
  call i2mex_getQDeltaP(nt1, ns, the, psi, zdum, ier)
  qdelp(1:nt1, 1:ns) = zdum(1:nt1, 1:ns)
  qdelp(nt1+1:nt1+4,:) = zdum(2:5,:)
!!$
!!$  print *,' qdelp(1,1:10)=',qdelp(1,1:10) 
!!$  print *,' qdelp(2,1:10)=',qdelp(2,1:10) 
!!$  print *,' qdelp(3,1:10)=',qdelp(3,1:10) 

!!$  do i=1, ns
!!$     zf(1:nt1) = gpa(i)*xjacob(1:nt1,i)/xsq(1:nt1,i) &
!!$          & + ga(i)*xjprym(1:nt1,i)/xsq(1:nt1,i) &
!!$          & - ga(i)*xjacob(1:nt1,i)*xsqdps(1:nt1,i)/xsq(1:nt1,i)**2
!!$     call i2mex_integratePeriodic1d(nt1, the, zf, nt1, the, qdelp(1,i))
!!$     qdelp(1:nt1,i) = qdelp(1:nt1,i) - qpa(i)*the(1:nt1)
!!$     qdelp(nt1+1:nt1+4,i) = qdelp(2:5,i)
!!$  enddo
!!$    ! extrapolate to axis
!!$
!!$    do i=1, nt1+4
!!$       qdelp(i,1) = FCCCC0(qdelp(i,2), qdelp(i,3), qdelp(i,4), qdelp(i,5), &
!!$            &                 psi(2),  psi(3),  psi(4),  psi(5), &
!!$            & psi(1))
!!$    enddo
!!$  print *,' qdelp(1,1:10)=',qdelp(1,1:10) 
!!$  print *,' qdelp(2,1:10)=',qdelp(2,1:10) 
!!$  print *,' qdelp(3,1:10)=',qdelp(3,1:10) 



!!$  call i2mex_getGsError(nt1, ns, the, psi, zdum2, iok)
!!$  call i2mex_2Dx(nt1, ns, the, psi, zdum2, 'gserror.dx', iok)


  deallocate(zdum, stat=iok)
  deallocate(zdum2, stat=iok)
  deallocate(zdum3, stat=iok)
  deallocate(zf, stat=iok)


  if(iok/=0) ier=1

end subroutine pstCompMetric

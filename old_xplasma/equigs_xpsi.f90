subroutine equigs_xpsi(id_xphi,id_q,xpsi_name,id_xpsi,ierr)

  !  using an xplasma sqrt(Phi/Philim) grid (id_xphi) and
  !  a q profile (id_q) construct a function sqrt(Psi/Psilim) with
  !  name (xpsi_name) and return its xplasma id (id_xpsi)

  implicit NONE

  integer, intent(in) :: id_xphi    ! xplasma id: x=sqrt(Phi/Philim) grid
  integer, intent(in) :: id_q       ! xplasma id: q(x) profile
  character*(*), intent(in) :: xpsi_name    ! name for sqrt(Psi/Psilim)

  integer, intent(out) :: id_xpsi   ! id of sqrt(Psi/Psilim) profile created.
  integer, intent(out) :: ierr      ! status code, 0=OK

  !----------------------------------------------------
  integer :: ix,ixqmin,inx,id_rho,lunmsg,inrho,idum
  real*8 :: zdelta,zqlim,zqi,zx,zdx,zs1,zs2
  real*8, dimension(:), allocatable :: x,q,qsm,xpsi

  real*8, parameter :: QLIM0=10.0d0 ! limit axial q value
  real*8, parameter :: ZERO=0.0d0
  real*8, parameter :: ONE=1.0d0
  real*8, parameter :: TWO=1.0d0
  real*8, parameter :: ZDMIN=0.05d0 ! min smoothing
  !----------------------------------------------------

  call eq_get_lunerr(lunmsg)

  id_xpsi=0
  ierr=1

  !--------------------------------------
  !  look up equilibrium grid to set smoothing parameter

  call eq_ganum('__RHO',id_rho)
  if(id_rho.eq.0) then
     write(lunmsg,*) ' ?equigs_xpsi: could not find __RHO grid.'
     return
  endif

  call eq_ngrid(id_rho,inrho)

  !  smoothing parameter:
  zdelta = max(ZDMIN,ONE/(inrho-1))

  !--------------------------------------
  !  fetch xphi grid (may or not be same as __RHO grid) and allocate
  !  and fetch q profile

  call eq_ngrid(id_xphi,inx)
  if(inx.eq.0) then
     write(lunmsg,*) ' ?equigs_xpsi: grid id_xphi = ',id_xphi,' shows 0 size.'
     return
  endif

  allocate(x(inx),q(inx),qsm(inx),xpsi(inx))
  
  do
     call eq_grid(id_xphi,x,inx,idum,ierr)
     if(ierr.ne.0) exit

     call eq_rgetf(inx,x,id_q,0,q,ierr)
     exit
  enddo

  if(ierr.ne.0) then
     deallocate(x,q,qsm,xpsi)
     return
  endif

  !--------------------------------------
  !  smooth q profile; impose QLIM0 inbound from a q profile local minimum...

  ixqmin=1

  do ix=inx,2,-1
     if(q(ix-1).gt.q(ix)) then
        ixqmin=ix
        exit
     endif
  enddo

  zqlim=max(QLIM0,2*q(ixqmin))

  do ix=1,ixqmin-1
     q(ix)=min(zqlim,q(ix))
  enddo

  qsm = q

  call r8_qksmooth(inx,x,qsm,zdelta,ZERO,ZERO,ierr) ! fix end pts

  if(ierr.ne.0) then
     deallocate(x,q,qsm,xpsi)
     return
  endif

  !--------------------------------------
  ! compute xpsi = Psi/Psilim = integral[0,x](dx*x/q) / integral[0,1](dx*x/q)

  xpsi(1) = 0
  do ix = 2,inx
     zdx=x(ix)-x(ix-1)
     zx=(x(ix)+x(ix-1))/TWO
     zqi=TWO/(qsm(ix)+qsm(ix-1))

     xpsi(ix) = xpsi(ix-1) + zdx*zx*zqi
  enddo

  ! normalize to [0,1]

  xpsi = xpsi/xpsi(inx)

  ! sqrt

  xpsi(2:inx) = sqrt(xpsi(2:inx))

  !--------------------------------------
  !  create profile

  zs1=(xpsi(2)-xpsi(1))/(x(2)-x(1))
  zs2=(xpsi(inx)-xpsi(inx-1))/(x(inx)-x(inx-1))

  call eqm_rhofun(101,id_xphi,xpsi_name,xpsi, &    ! Akima Hermite, xpsi' BC
       1,zs1,1,zs2,id_xpsi,ierr)

  deallocate(x,q,qsm,xpsi)

  return
end subroutine equigs_xpsi

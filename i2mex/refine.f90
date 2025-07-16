subroutine i2mex_refineESCPAndQ(mpol, ier)

  ! Refine the equilibrium using the ESC equilibrium code in p and q mode

  use i2mex_mod
  use ezspline
  use ezspline_obj
  implicit none
  integer, intent(in) :: mpol ! number of Fourier harmonics
  integer, intent(out) :: ier

  real(i2mex_r8), dimension(:), allocatable :: psi, phi, p, q, g, sstor, the
  real(i2mex_r8), dimension(:,:), allocatable :: x, z
  real(i2mex_r8), dimension(:), allocatable :: t
  real(i2mex_r8), dimension(:,:), allocatable :: cspl
  real(i2mex_r8) phimin, phimax
  integer, parameter :: na1 = 21
  real(i2mex_r8) sesc(na1), pesc(na1), qesc(na1), wk(na1)
  real(i2mex_r8) pesc0(na1), jesc(na1), r0, b0, rbtor
  real(i2mex_r8) det, x0, z0, x1, z1, x2, z2

  integer, parameter :: nPt1 = 33
  real(i2mex_r8), dimension(nPt1) :: xbound, zbound, tbound
  character*10, parameter :: runid = "NSTX.00000"
  type(ezspline1_r8) :: f_spl

  integer nt1, ns, iok, i, niter

  ier = 0

  call i2mex_getOriNt1(nt1, iok)
  call i2mex_getOriNs(ns, iok)
  ALLOCATE(x(nt1, ns))
  ALLOCATE(z(nt1, ns))
  ALLOCATE(t(nt1))
  ALLOCATE(cspl(4, ns))
  ALLOCATE(the(nt1))
  ALLOCATE(psi(ns))
  ALLOCATE(phi(ns))
  ALLOCATE(sstor(ns))
  ALLOCATE(p(ns))
  ALLOCATE(q(ns))
  ALLOCATE(g(ns))
  call i2mex_getOriPsi(ns, psi, iok)
  call i2mex_getOriT(nt1, the, iok)
  call i2mex_getPhi(ns, psi, phi, iok)
  phimin = minval(phi)
  phimax = maxval(phi)
  sstor = sqrt((phi-phimin)/(phimax-phimin))
  call i2mex_getP(ns, psi, p, iok)
  call i2mex_getQ(ns, psi, q, iok)
  
  t = i2mex_twopi_r8* (/ &
       ( real(i-1, i2mex_r8)/real(nt1-1, i2mex_r8), i=1, nt1 ) &
       & /)

  call i2mex_getX(nt1, ns, the, psi, x, iok) 
  call i2mex_getZ(nt1, ns, the, psi, z, iok) 

 if(iok/=0) ier = 113

  tbound = i2mex_twopi_r8* (/ &
       ( real(i-1, i2mex_r8)/real(nPt1-1, i2mex_r8), i=1, nPt1 ) &
       & /)
  call i2mex_interpPeriodic1d(nt1, t, x(1, ns), nPt1, tbound, xbound)
  call i2mex_interpPeriodic1d(nt1, t, z(1, ns), nPt1, tbound, zbound)


  ! compute the profiles on the uniform sstor mesh

  sesc = (/ (real(i-1,i2mex_r8)/real(na1-1,i2mex_r8), i=1, na1) /)

  ! pressure

  call EZspline_init(f_spl, ns, (/0,0/), iok)
  f_spl%x1 = sstor
  call EZspline_setup(f_spl, p, iok) ! coefficient set up
  call EZspline_interp(f_spl, na1, sesc, pesc, iok)
  call EZspline_free(f_spl, ier)

  ! safety factor

  call EZspline_init(f_spl, ns, (/0,0/), iok)
  f_spl%x1 = sstor
  call EZspline_setup(f_spl, q, iok) ! coefficient set up
  call EZspline_interp(f_spl, na1, sesc, qesc, iok)
  call EZspline_free(f_spl, ier)

  if(iok/=0) ier = 114

  ! prepare esc

  call esc0(runid)

  
  ! run in p=0 & j// mode first

  call i2mex_getB0(b0, iok)
  call i2mex_getRmajor(r0, iok)
  rbtor = r0*b0
  jesc = 2._i2mex_r8*b0*(1.0_i2mex_r8 - sesc*sesc)/(r0*q(1))
  niter = 1
  do i = 1, niter
     pesc0 = 0.0 ! real(i-1,i2mex_r8)*pesc/real(niter-1,i2mex_r8)
     call esc_p_jparallel(pesc0, jesc, xbound, zbound, nPt1, &
     & rbtor, r0)
     call esc(iok, 0, mpol, 0.0_i2mex_r8)
     if(iok/=0) then
        print *,'--warning iok=',iok,' when running in p & j// mode'
        iok = 115
     endif
  enddo


  ! run in p and q mode

  niter = 5
  do i = 1, niter
     pesc0 = real(i-1,i2mex_r8)*pesc/real(niter-1,i2mex_r8)
     call esc_p_q(pesc0, qesc, xbound, zbound, nPt1, &
     & rbtor, r0)
     call esc(iok, 0, mpol, 0.0_i2mex_r8)
     if(iok/=0) then
        print *,'--warning iok=',iok,' when running in p & q mode (iter=',i,')'
        ier = 116
     endif
  enddo

  call escGetPsibar(psi, ns)
  call escGetP(p, ns)
  call escGetQ(q, ns)
  call escGetG(g, ns)
  call escGetR(x, nt1, ns, i2mex_o%isThetaClockwise)  
  call escGetZ(z, nt1, ns, i2mex_o%isThetaClockwise) 


  ! clean up

  call esc1


  call i2mex_free(iok)

  call i2mex_init("ESC", nt1, ns, the, psi, p, g, q, x, z, iok)
  
  DEALLOCATE(x)
  DEALLOCATE(z)
  DEALLOCATE(t)
  DEALLOCATE(cspl)
  DEALLOCATE(the)
  DEALLOCATE(psi)
  DEALLOCATE(phi)
  DEALLOCATE(sstor)
  DEALLOCATE(p)
  DEALLOCATE(q)
  DEALLOCATE(g)
  
end subroutine i2mex_refineESCPAndQ
  

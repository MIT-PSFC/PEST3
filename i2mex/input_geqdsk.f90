subroutine i2mex_fromGeqdsk(filename, it_orientation, ier)
 
  ! read G-EQDSK file and and perform inverse map
 
   use i2mex_mod
   use geqdsk_mod
   use ezspline_obj
   use ezspline
   use cont_mod
 
 
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
 
  character*(*), intent(in) :: filename
  integer, intent(in) :: it_orientation
  integer, intent(out) :: ier
 
  type(geqdsk) :: geq
  type(ezspline2_r8) :: pspl
 
  real(r8), dimension(:), allocatable :: p, q, g, psi, psi_geq, psis, &
       & the, xgrid, zgrid, tlike, q_est
  real(r8), dimension(:, :), allocatable :: x, z
  real(r8), dimension(:, :), allocatable :: xin, zin
  real(r8) :: pmax, pmin, xmin, xmax, zmin, zmax, surf, xm, zm, xp, zp, &
       & pmag, psep, ttot
  integer iok, i, j, ioutboard(1)
 
  real(r8), dimension(:,:), allocatable :: integrand
  real(r8), dimension(:), allocatable :: qres
  real(r8) x0, z0, x1, z1, x2, z2, det
  real(r8) ztol, delta, dlim
 
  ier = 0
 
  call geq_init(geq, filename, iok)
  call geq_error(iok)
  if(iok/=0) then
     ier = 117
     return
  endif
 
  xmin = geq%RLEFT_M
  xmax = geq%RLEFT_M + geq%RDIM_M
  zmin = geq%ZMID_M - geq%ZDIM_M/2.0_r8
  zmax = geq%ZMID_M + geq%ZDIM_M/2.0_r8
 
 
  allocate(xgrid(geq%nw))
  allocate(zgrid(geq%nh))
 
  xgrid=xmin + (xmax - xmin)*(/ (real(i-1,r8)/real(geq%NW-1,r8), i=1,geq%NW) /)
  zgrid=zmin + (zmax - zmin)*(/ (real(i-1,r8)/real(geq%NH-1,r8), i=1,geq%NH) /)
 
!!$  i2mex_nt1 = 257
!!$  i2mex_ns = 401
 
  allocate(psi(i2mex_ns), psi_geq(geq%NW), p(i2mex_ns), g(i2mex_ns), q(i2mex_ns), the(i2mex_nt1))
  allocate(tlike(i2mex_ns), q_est(i2mex_ns))
  allocate(x(i2mex_nt1, i2mex_ns), z(i2mex_nt1, i2mex_ns))
  allocate(xin(i2mex_nt1, i2mex_ns), zin(i2mex_nt1, i2mex_ns))
  allocate(psis(i2mex_nt1))
 
  xm = geq%RMAXIS_M
  zm = geq%ZMAXIS_M
  ioutboard = maxloc(geq%RBBBS_M) ! most outboard point
  xp = geq%RBBBS_M(ioutboard(1))
  zp = geq%ZBBBS_M(ioutboard(1))
 
  call ezspline_init(pspl, geq%NW, geq%NH, (/0,0/), (/0,0/), iok)
  call ezspline_error(iok)
 
  ! toroidal current sign
  cont_sign = 1.0_r8
  if(geq%SIBRY_Wb__Rad < geq%SIMAG_Wb__Rad) cont_sign = -1.0_r8
 
  ! We like to have psi = 0 on magnetic axis, and increasing
  ! outwards (this will change the current orientation but that's
  ! ok).
  geq%psirz_Wb__Rad = cont_sign*(geq%psirz_Wb__Rad - geq%SIMAG_Wb__Rad)
  geq%SIBRY_Wb__Rad = cont_sign*(geq%SIBRY_Wb__Rad - geq%SIMAG_Wb__Rad)
  geq%SIMAG_Wb__Rad = 0.0_r8
  geq%FFPRIM_T2M2Rad__Wb = cont_sign*geq%FFPRIM_T2M2Rad__Wb
  geq%pprime_NtRad__M2Wb = cont_sign*geq%pprime_NtRad__M2Wb
  geq%CURRENT_A = cont_sign*geq%CURRENT_A
  ! set back
  cont_sign = 1.0_r8
 
  pspl%x1 = xgrid
  pspl%x2 = zgrid
  !pspl%isHermite = 1 ! Akima Hermite
  call ezspline_setup(pspl, geq%PSIRZ_Wb__Rad, iok)
  call ezspline_error(iok)
 
  ! interpolate to get psi on axis
  call ezspline_interp(pspl, xm, zm, pmag, iok)
  call ezspline_error(iok)
  ! psi on separatrix
  call ezspline_interp(pspl, xp, zp, psep, iok)
  call ezspline_error(iok)
 
  do j = 2, i2mex_ns
     ! rough estimate for q (used to determine angle integration step)
     cont_q = geq%qpsi(1) + &
          & real(j,r8)*(geq%qpsi(geq%NW) - geq%qpsi(1))/real(i2mex_ns,r8)
 
     surf = real(j-1, r8)/real(i2mex_ns-1, r8)
     xp = xm + surf * i2mex_LAST_NORM_SURFACE_IS*(geq%RBBBS_M(ioutboard(1)) - xm)
     zp = zm + surf * i2mex_LAST_NORM_SURFACE_IS*(geq%ZBBBS_M(ioutboard(1)) - zm)
     call i2mex_contour(geq%NW, geq%NH, xgrid, zgrid, geq%PSIRZ_Wb__Rad,  &
       & xm, zm, xp, zp, i2mex_nt1, ttot, xin(1, j), zin(1, j), iok)
     call i2mex_error(iok)
 
     ! get psi on contour
     call ezspline_interp(pspl, i2mex_nt1, xin(:, j), zin(:, j), psis, iok)
     call ezspline_error(iok)
 
     ! average psi on contour
     psi(j) = sum(psis)/real(i2mex_nt1, r8)
 
     ! save tot time
     tlike(j) = ttot
 
  enddo
  if(iok/=0) ier = 120
  psi(1) = pmag
 
  xin(1:i2mex_nt1, 1) = xm
  zin(1:i2mex_nt1, 1) = zm
 
  ! normalized to psi on magnetic axis = 0
  psi = psi - psi(1)
  pmax = MAXVAL(psi)
  pmin = MINVAL(psi)
 
  psi_geq = 0._r8 + ((pmax-pmin)/i2mex_LAST_NORM_SURFACE_IS)* &
       & (/ (real(i-1, r8)/real(geq%nw-1,r8), i=1, geq%nw) /)
 
 
  call ezspline_free(pspl, iok)
 
  the = i2mex_twopi_r8* (/ (real(i-1,r8)/real(i2mex_nt1-1,r8), i=1, i2mex_nt1) /)
 
  ! pressure
  call i2mex_interp1d(geq%NW, psi_geq, &
       & geq%PRES_Nt__M2 * 4.e-7_r8 * i2mex_pi_r8, &
       & i2mex_ns, psi, p)
  ! q
  call i2mex_interp1d(geq%NW, psi_geq, &
       & geq%QPSI                                , &
       & i2mex_ns, psi, q)
  ! g
  call i2mex_interp1d(geq%NW, psi_geq, &
       & abs(geq%FPOL_TM)                        , &
       & i2mex_ns, psi, g)
 
  ! estimate for q
  q_est(2:i2mex_ns) = abs(g(2:i2mex_ns))*tlike(2:i2mex_ns)/(i2mex_twopi_r8)
  q_est(1) = q_est(2) + (q_est(3)-q_est(2))*(psi(1)-psi(2))/(psi(3)-psi(2))
 
  q = q_est ! take our own q
 
  ! initialize i2mex
  ! get theta orientation
  x0 = (maxval(xin(:,1))+minval(xin(:,1)))/2._i2mex_r8
  z0 = (maxval(zin(:,1))+minval(zin(:,1)))/2._i2mex_r8
  x1 = xin(1, i2mex_ns)
  z1 = zin(1, i2mex_ns)
  x2 = xin(2, i2mex_ns)
  z2 = zin(2, i2mex_ns)
  det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)
 
 
  if( it_orientation*det < 0.0_i2mex_r8 ) then
     do i = 1, i2mex_nt1
        x(i, 1:i2mex_ns) = xin(i, 1:i2mex_ns)
        z(i, 1:i2mex_ns) = zin(i, 1:i2mex_ns)
     enddo
  else
     do i = 1, i2mex_nt1
        x(i, 1:i2mex_ns) = REAL(xin(i2mex_nt1-i+1, 1:i2mex_ns), i2mex_r8)
        z(i, 1:i2mex_ns) = REAL(zin(i2mex_nt1-i+1, 1:i2mex_ns), i2mex_r8)
     enddo
  endif
 
  ! to remap poloidally set i2mex_remap=1
  ! Jacobian ~ x^mx /(|grad psi|^npsi B^kb)
  if(i2mex_remap/=0) then
     call i2mex_poloidalRemap(i2mex_mx, i2mex_npsi, i2mex_kb, &
          & i2mex_nt1, i2mex_ns, the, psi, g, x, z, iok)
  endif
  
  call i2mex_init('GEQDSK_'//filename, &
       i2mex_nt1, i2mex_ns, the, psi, p, g, q, x, z, iok)
  call i2mex_error(iok)
 
  if(i2mex_direct>0) then
        ! Build Psi(R,Z) representation, if required...
        ztol  = 1.0e-6*xmax
        dlim  = 0.05_r8 * Xmax ! box boundaries
        delta = 0.0_r8         ! smoothing parameter

        i2mex_o%Rleft = xmin
        i2mex_o%Rrigh = xmax
        i2mex_o%Zbot = zmin
        i2mex_o%Ztop = zmax
        

        call eqm_rzgrid(xgrid, zgrid, -1, -1, geq%NW, geq%NH, ztol, &
             & i2mex_o%id_Rgrid, i2mex_o%id_Zgrid, iok)
        if(iok/=0) then
           ier = 151
           call i2mex_error(ier)
        endif

        call eqm_cbdy(5, &
             & (/xmin, xmax, xmax, xmin, xmin /), &
             & (/zmin, zmin, zmax, zmax, zmin /), &
             & iok)
        if(iok/=0) ier = 154

  
        ! interpolation order is Akima
        call eqm_rzfunda('PSI', i2mex_o%id_psirz, geq%psirz_Wb__Rad, &
             & geq%NW, geq%NH, i2mex_akima, delta, iok)
        if(iok/=0) ier = 153
  endif
 
  deallocate(xgrid)
  deallocate(zgrid)
  deallocate(psi, psi_geq, p, g, q, the)
  deallocate(tlike, q_est)
  deallocate(x, z)
  deallocate(xin, zin)
  deallocate(psis)
 
  call geq_free(geq, iok)
  if(iok/=0) ier = 150
 
 
end subroutine i2mex_fromGeqdsk
 
subroutine i2mex_fromGeqdskThruEscQ(filename, it_orientation, ier)
 
  ! read G-EQDSK file, take boundary coodinates, pressure and q profiles
  ! and compute equilibrium using Esc
 
   use i2mex_mod
   use geqdsk_mod
   use cont_mod
 
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
 
  character*(*), intent(in) :: filename
  integer, intent(in) :: it_orientation
  integer, intent(out) :: ier
 
  type(geqdsk) :: geq
  real(r8), dimension(:), allocatable :: pgeq, qgeq, psigeq, sgeq
  integer, parameter :: na1 = 21 ! esc input profile grid size
  real(r8), dimension(na1) :: pesc, pesc0, jesc, qesc, sesc
  real(r8), dimension(:), allocatable :: psi, p, q, g, the
  real(r8), dimension(:, :), allocatable :: x, z
  real(r8) :: r0btor, r0, b0, pmag, psep
  integer iok, i, mpol, niter, ioutboard(1)
 
  integer, parameter :: nbound1 = 33
  real(r8) :: xm, zm, xp, zp, xmin, xmax, zmin, zmax, ttot
  real(r8), dimension(:), allocatable :: xgrid, zgrid
  real(r8), dimension(nbound1) :: xbound, zbound
 
  ier = 0
 
  call geq_init(geq, filename, iok)
  call geq_error(iok)
  if(iok/=0) then
     ier = 117
     return
  endif
 
  ! toroidal current sign
  cont_sign = 1.0_r8
  if(geq%SIBRY_Wb__Rad < geq%SIMAG_Wb__Rad) cont_sign = -1.0_r8
 
  ! We like to have psi = 0 on magnetic axis, and increasing
  ! outwards (this will change the current orientation but that's
  ! ok).
  geq%psirz_Wb__Rad = cont_sign*(geq%psirz_Wb__Rad - geq%SIMAG_Wb__Rad)
  geq%SIBRY_Wb__Rad = cont_sign*(geq%SIBRY_Wb__Rad - geq%SIMAG_Wb__Rad)
  geq%SIMAG_Wb__Rad = 0.0_r8
  geq%FFPRIM_T2M2Rad__Wb = cont_sign*geq%FFPRIM_T2M2Rad__Wb
  geq%pprime_NtRad__M2Wb = cont_sign*geq%pprime_NtRad__M2Wb
  geq%CURRENT_A = cont_sign*geq%CURRENT_A
  ! set back
  cont_sign = 1.0_r8
 
  ! recompute boundary points
  xmin = geq%RLEFT_M
  xmax = geq%RLEFT_M + geq%RDIM_M
  zmin = geq%ZMID_M - geq%ZDIM_M/2.0_r8
  zmax = geq%ZMID_M + geq%ZDIM_M/2.0_r8
 
  allocate(xgrid(geq%nw))
  allocate(zgrid(geq%nh))
 
  xgrid=xmin + (xmax - xmin)*(/ (real(i-1,r8)/real(geq%NW-1,r8), i=1,geq%NW) /)
  zgrid=zmin + (zmax - zmin)*(/ (real(i-1,r8)/real(geq%NH-1,r8), i=1,geq%NH) /)
 
  xm = geq%RMAXIS_M
  zm = geq%ZMAXIS_M
  ! starting point
  ioutboard  = maxloc(geq%RBBBS_M)
  i = ioutboard(1)
  xp = xm + i2mex_LAST_NORM_SURFACE_IS * (geq%RBBBS_M(i) - xm)
  zp = zm + i2mex_LAST_NORM_SURFACE_IS * (geq%ZBBBS_M(i) - zm)
 
  ! integration direction
  cont_sign = 1.0_r8
  if(geq%SIBRY_Wb__Rad < geq%SIMAG_Wb__Rad) cont_sign = -1.0_r8
 
  cont_q = geq%qpsi(geq%NW)

  call i2mex_contour(geq%NW, geq%NH, xgrid, zgrid, geq%PSIRZ_Wb__Rad,  &
       & xm, zm, xp, zp, nbound1, ttot, xbound, zbound, iok)
  call i2mex_error(iok)
 
  deallocate(xgrid)
  deallocate(zgrid)
 
 
  ! radial profiles are vs poloidal flux
  ! compute toroidal flux
 
  i2mex_ns = geq%nw
 
  ALLOCATE(pgeq(i2mex_ns), qgeq(i2mex_ns), psigeq(i2mex_ns), sgeq(i2mex_ns))
 
  ! Pascals -> mu0*Pa
  pgeq = geq%pres_Nt__M2 * i2mex_fourpi_r8 * 1.0e-7_i2mex_r8
  qgeq = geq%qpsi
 
 
  psigeq = abs(geq%SIBRY_Wb__Rad - geq%SIMAG_Wb__Rad) * &
       & (/ ( real(i-1, r8)/real(i2mex_ns-1, r8), i = 1, i2mex_ns ) /)
 
  ! integrate q d psi to get toroidal flux
  call i2mex_integrate1d(i2mex_ns, psigeq, qgeq, i2mex_ns, psigeq, sgeq)
  sgeq = sgeq/sgeq(i2mex_ns)
  sgeq = sqrt(sgeq) !  = sqrt(norm toroidal flux)
 
  ! esc wants the radial profiles on equidistant sqrt(norm phi) mesh
  sesc = (/ (real(i-1, r8)/real(na1-1, r8), i = 1, na1) /)
 
  ! interpolate the profiles onto the new mesh
  call i2mex_interp1d(i2mex_ns, sgeq, qgeq, na1, sesc, qesc)
  call i2mex_interp1d(i2mex_ns, sgeq, pgeq, na1, sesc, pesc)
 
  !
  ! run esc
  !
  r0 = geq%RCENTR_M
  b0 = abs(geq%BCENTR_T)
  r0btor = r0*b0
 
  call esc0("dum")
 
  ! run in p & j// mode first
 
  mpol = 4
 
  niter = 4
  do i = 1, 1
     pesc0 = real(i-1,i2mex_r8)*pesc/real(niter-1,i2mex_r8)
     jesc = 2._i2mex_r8*b0*(1.0_i2mex_r8 - sesc*sesc)/(r0*qesc(1))
     call esc_p_jparallel(pesc0, jesc, xbound, zbound, nbound1, &
          & r0btor, r0)
     call esc(iok, 0, mpol, 0.0_i2mex_r8)
     if(iok/=0) then
        print *,'--warning error flag ',iok,' raised when running ESC in p & j// mode'
        ier = 118
     endif
  enddo
 
  ! run in p and q mode
 
  niter = 4
  do i = 1, niter
     pesc0 = real(i-1,i2mex_r8)*pesc/real(niter-1,i2mex_r8)
     call esc_p_q(pesc0, qesc, xbound, zbound, nbound1, &
     & r0btor, r0)
     call esc(iok, 0, mpol, 0.0_i2mex_r8)
     if(iok/=0) then
        print *,'--warning error flag ',iok,' raised when running ESC in p & q mode (iter=',i,')'
        ier = 118
     endif
  enddo
 
  i2mex_nt1 = 65
  ALLOCATE(psi(i2mex_ns), p(i2mex_ns), q(i2mex_ns), g(i2mex_ns))
  ALLOCATE(x(i2mex_nt1, i2mex_ns), z(i2mex_nt1, i2mex_ns))
  ALLOCATE(the(i2mex_nt1))
 
  call escGetPsibar(psi, i2mex_ns)
  call escGetP(p, i2mex_ns)
  call escGetQ(q, i2mex_ns)
  call escGetG(g, i2mex_ns)
  call escGetR(x, i2mex_nt1, i2mex_ns, it_orientation)
  call escGetZ(z, i2mex_nt1, i2mex_ns, it_orientation)
 
 
  ! clean up
 
  call esc1
 
  ! is psi monotonic increasing?
  ioutboard = maxloc(psi)
  i = ioutboard(1)
  if(i /= i2mex_ns) then
     ier = 119
     ! back off from plasma edge
     i2mex_ns = i
  endif
 
  !
  ! initialize
  !
 
  the = i2mex_twopi_r8*(/ ( real(i-1, r8)/real(i2mex_nt1-1, r8), i = 1, i2mex_nt1 ) /)
  psi = abs(psi - psi(1))
 
  ! to remap poloidally set i2mex_remap=1
  ! Jacobian ~ x^mx /(|grad psi|^npsi B^kb)
  if(i2mex_remap/=0) then
     call i2mex_poloidalRemap(i2mex_mx, i2mex_npsi, i2mex_kb, &
          & i2mex_nt1, i2mex_ns, the, psi, g, x, z, ier)
  endif
 
  call i2mex_init('GEQDSK_escQ_'//filename, &
       i2mex_nt1, i2mex_ns, the, psi, p, g, q, x, z, iok)
  call i2mex_error(iok)
 
  DEALLOCATE(the)
  DEALLOCATE(pgeq, qgeq, psigeq, sgeq)  
  DEALLOCATE(psi, p, q, g)
  DEALLOCATE(x, z)
 
  call geq_free(geq, iok)
 
end subroutine i2mex_fromGeqdskThruEscQ
 

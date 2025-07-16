program mex2gtc

  use ezcdf
  use i2mex_mod
  use ezspline_obj
  use ezspline
  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)


!!!!!  integer, parameter :: nt1=129, ns=101
  integer nt1, ns  

  real(r8) :: psi_min, psi_max, xmax, xmin, zmax, zmin, rerror
  real(r8) :: eps, area, vol, elong, beta, betap_, betat, betan, ima, li
  real(r8) :: r0, rm, scale
  real(r8) :: torfluxmin, torfluxmax, psidum, polfun
  real(r8), dimension(:), allocatable ::  psii
  real(r8), dimension(:), allocatable ::  t,  psi, p, g, q
  real(r8), dimension(:), allocatable ::  torfluxi, torflux
  real(r8), dimension(:,:), allocatable ::  xa, za,  xj, delta, Bfield, nufunc
  real(r8), dimension(:,:), allocatable ::  gtt, gtp, gpp
  real(r8), dimension(:,:), allocatable ::  gserror
  real(r8), dimension(:,:), allocatable ::  temp
  character :: answ

  type(EZspline1_r8) :: spline_1d! 1-d EZspline object

  integer i, j, ier, inputFormat, nsi, nti
  character*130 inputFile
  character*12 filename, cdum
  real(r8) time
  integer m, n, k, tic, toc, tics_per_second, dims(3), iu

  character*8 today
  character*10 now
  character*80 lon_date

  ! Remap coordinate system
  i2mex_remap=1
  ! set jacobian for Boozer coordinates
!!!  i2mex_mx=0; i2mex_npsi=0; i2mex_kb=2 ! Boozer 
  ! set jacobian for Equal-arc coordinates
  i2mex_mx=1; i2mex_npsi=1; i2mex_kb=0 ! Equal-Arc
!!!  i2mex_direct = 0  ! =1 to allow direct representation to co-exist

  write(*,*)' Choose your theta variable by specifying the Jacobian dependence'
  write(*,*)' J ~ X**m/(|grad psi|**n |B|**k). Popular choices are:'
  write(*,*)' PEST     : m=2, n=0, k=0'
  write(*,*)' Boozer   : m=0, n=0, k=2'
  write(*,*)' Equal-arc: m=1, n=1, k=0'
  write(*,'("Enter m,n,k: ",$)')
  read(*,*) i2mex_mx, i2mex_npsi, i2mex_kb

  write(*,*)' Enter the input equilibrium file format'
  write(*,*)' -1  for CHEASE INP1 format.'
  write(*,*)'  0  for TRANSP data in UFILE format (requires MDSPlus connection)'
  write(*,*)' +1  for CHEASE inp1.cdf format'
  write(*,*)' +2  for JSOLVER eqdsk.cdf format'
  write(*,*)' +3  for EFIT G-EQDSK format'
  write(*,*)' +4  for EFIT G-EQDSK format, rerunning equilibrium through ESC-Q'
  write(*,*)' +5  for Menard''s psipqgRZ netCDF format'
  read(5,*) inputFormat
  if (inputFormat<-1 .or. inputFormat>5) inputFormat=1
  if (inputFormat==0) then
     write(*,*) ' Enter the path or MDSPlus node (requires MDSPlus connection)'
     write(*,*) ' Example:'
     write(*,*) 'MDS+:BIRCH.PPPL.GOV:TRANSP(NSTX.00,11115P07)'
     read(5,'(a80)') inputFile
     write(*,*) ' Enter the time slice in sec (data will be interpolated)'
     read(5,*) time
     call i2mex_fromMDSPlus(inputFile, &
          & time, nt1, ns, i2mex_counterclockwise, ier)
  else
     write(*,*) ' Enter input file name '
     if(inputFormat == 3 .or. inputFormat==4) then
        write(*,*) ' This can be an MDSPlus path as in '
        write(*,*) ' "MDS+/REDUCE:BIRCH.PPPL.GOV:8501:EFIT01(103984;t=0.23)"'
        write(*,*) ' for instance'
     endif
     if(inputFormat==-1) then
        write(*,*) ' Depending on your platform, you have the choice between:'
        write(*,*) ' INP1_ieee-le for binary in little endian format (Alphas, ix86...), or'
        write(*,*) ' INP1_ieee-be for binary in big endian format (Sun, ...). '
     endif
     read(5,*) inputFile
     call i2mex_ReadEquilibrium(inputFormat, inputFile, i2mex_clockwise, ier)

  endif

  call i2mex_error(ier)

! Get dimensions of original grid from the input file
  call i2mex_getOriNs(nsi, ier) 
  call i2mex_getOriNt1(nti, ier) 
  write(*,*) 'Original grid: n_surfaces =',nsi,'  n_theta =',nti
! Set new grid dimensions to original grid
  write(*,'("New grid size ns,nt1 (default=original): ",$)')
  read(*,'(a12)') cdum
  if (cdum(1:1)/="") then
     read(cdum,*) ns,nt1
  else
     ns = nsi
     nt1 = nti
  endif
  write(*,*) 'New grid: ns =',ns,'  nt1 =',nt1

  ALLOCATE(t(nt1))
  ALLOCATE(psii(nsi))
  ALLOCATE(psi(ns), p(ns), q(ns), g(ns))
  ALLOCATE(xa(nt1, ns), za(nt1, ns), xj(nt1, ns), delta(nt1, ns),&
         & Bfield(nt1, ns), nufunc(nt1,ns))
  ALLOCATE(gtt(nt1, ns), gtp(nt1, ns), gpp(nt1, ns))
  ALLOCATE(gserror(nt1, ns))
  ALLOCATE(temp(nt1, ns))

! Define the theta grid as d_theta = constant
  t = i2mex_twopi_r8* (/ (real(i-1, r8)/real(nt1-1, r8), i=1, nt1) /)

! Get the original poloidal flux grid (psii)
  call i2mex_getOriPsi(nsi, psii, ier) 
  psi_min = psii(1); psi_max = psii(nsi)

! Now GTC and Orbit3d use the "toroidal flux" as their flux surface (radial)
! coordinate. We need to convert all quantities from a (theta,pol_flux) grid
! to a (theta,tor_flux) grid. We use the i2mex_getPhi function to do this.
! First we get the toroidal flux on the original psii grid points:
  allocate(torfluxi(nsi))
  call i2mex_getPhi(nsi, psii, torfluxi, ier)
  torfluxmin = torfluxi(1)
  torfluxmax = torfluxi(nsi)

! We want to end up with a uniform grid in the toroidal flux coordinate
! so let's define the torflux array of "ns" elements as:
  allocate(torflux(ns))
  torflux = torfluxmin + (torfluxmax-torfluxmin)* &
       & (/ (real(i-1, r8)/real(ns-1, r8), i=1, ns) /)
!!!  torflux(1) = torflux(1) + 1.d-20
!!!  psi = psi_min + (psi_max-psi_min)* &
!!!       & (/ (real(i-1, r8)/real(ns-1, r8), i=1, ns) /)

  call EZspline_init(spline_1d, nsi, (/0,0/), ier)
  call EZspline_error(ier)

  spline_1d%x1 = torfluxi
  call EZspline_setup(spline_1d, psii, ier)
  call EZspline_error(ier)

  call EZspline_interp(spline_1d, nsi, torflux, psi, ier)
  call EZspline_error(ier)

! We need to find the poloidal flux grid that corresponds to our uniform
! toroidal flux grid. To do this, we use i2mex_getPhi again but we have
! to iterate while changing the value of psi, the poloidal flux.
! The first point is already known from psii(1)

  
  call i2mex_getP(ns, psi, p, ier)
  call i2mex_error(ier)
  call i2mex_getQ(ns, psi, q, ier)
  call i2mex_error(ier)
  call i2mex_getG(ns, psi, g, ier)
  call i2mex_error(ier)
  call i2mex_getX(nt1, ns, t, psi, xa, ier)
  call i2mex_error(ier)
  call i2mex_getZ(nt1, ns, t, psi, za, ier)
  call i2mex_error(ier)
  call i2mex_getJ(nt1, ns, t, psi, xj, ier)
  call i2mex_error(ier)

! The Jacobian must be multiplied by q since we want it in terms of the
! toroidal flux, and not the poloidal flux.
  do i=1,ns
     xj(:,i) = xj(:,i)*q(i)
  enddo

  call i2mex_getMetric(nt1, ns, t, psi, gtt, gtp, gpp, ier)
  call i2mex_error(ier)

! gtp and gpp are given in terms of the poloidal flux (psi_pol) but we need
! them in terms of the toroidal flux (psi_tor). Since we have:
!         Grad(psi_poloidal) = Grad(psi_toroidal)/q 
! then gtp (Grad(theta).Grad(psi_pol)) becomes gtp_new = gtp/q and
! gpp (Grad(psi_pol).Grad(psi_pol) becomes  gpp_new = gpp/(q*q)
  do i=1,ns
     gtp(:,i) = gtp(:,i)/q(i)
     gpp(:,i) = gpp(:,i)/(q(i)*q(i))
  enddo

  call i2mex_getDelta(nt1, ns, t, psi, delta, ier)

! Evaluate nu(theta,psi) = q(psi)*delta(theta,psi). This is needed to
! calculate zeta and its derivatives --> zeta = phi - nu(theta,psi)
! where phi is the cylindrical geometry angle (the PEST1 code uses phi
! as its toroidal angle. See Grimm et al., "Methods in Computational
! Physics, volume 16", p. 253).
  do i=1,ns
     nufunc(:,i) = q(i)*delta(:,i)
  enddo


  ! check GS error
  call i2mex_getGsError(nt1, ns, t, psi, gserror, ier)
  call i2mex_error(ier)
  call i2mex_2Dx(nt1, ns, t, psi, gserror, 'gserror.dx', ier)
  call i2mex_error(ier)
  rerror = sum(sum(gserror, dim=1))/(real(nt1,r8)*real(ns,r8))
  write(*,*)' '
  write(*,'(a,e10.2)')' Average Grad-Shafranov error: ', rerror
  if(abs(rerror)>0.01_r8) write(*,*) ' --Warning--: large rel GS error > 1%!'
  if(abs(rerror)>1._r8  ) then 
     write(*,*) ' --WARNING--: HUGE rel GS error > 100%!'
     write(*,*) ' Unless this is due to a highly localized error'
     write(*,*) ' (typically near the separatrix), this equilibrium'
     write(*,*) ' is worthless. To find out, send file gserror.dx'
     write(*,*) ' to pletzer@pppl.gov.'
     write(*,*) ' Continue at your own risk.'
  endif
     

  xmax = maxval(maxval(xa,dim=1))
  xmin = minval(minval(xa,dim=1))
  zmax = maxval(maxval(za,dim=1))
  zmin = minval(minval(za,dim=1))
  eps = (xmax-xmin)/(xmax+xmin)
  r0 = (xmax+xmin)/2.0_r8
  rm = sum(xa(1:nt1-1,1))/real(nt1-1,r8)
  elong = (zmax-zmin)/(xmax-xmin)
  call i2mex_getSurface(area, ier)
  call i2mex_getVolume(vol, ier)
  call i2mex_getBeta(nt1, ns, t, psi, beta, ier)
  call i2mex_getBetaPoloidal(nt1, ns, t, psi, betap_, ier)
  call i2mex_getBetaToroidal(nt1, ns, t, psi, betat, ier)
  call i2mex_getBetaN(nt1, ns, t, psi, betan, ier)
  call i2mex_getPlasmaCurrent(nt1, ns, t, psi, ima, ier)
  ima = ima/(0.4_r8*i2mex_pi_r8)
  call i2mex_getLi(nt1, ns, t, psi, li, ier)
  write(*,*)' '
  write(*,*)'  a/R0    Area   Volume   Elong   Beta'
  write(*,'(5f8.3)') eps, area, vol, elong, beta
  write(*,*)' '
  write(*,*)' Beta-p  Beta-t  Beta-N    I-MA     li'
  write(*,'(5f8.3)') betap_, betat, betan, ima, li
  write(*,*)' '

  write(*,*)' '
  write(*,*)           'Toroidal magnetic field strength in:        plasma   vacuum'
  write(*,'(a,3f9.4)') ' At geometric centre position R0=',r0,  g(1)/r0,     g(ns)/r0
  write(*,'(a,3f9.4)') ' At magnetic axis position    Rm=',rm,  g(1)/rm,     g(ns)/rm
  write(*,*)'Do you want to rescale B so as to get 1 [T] on the magnetic axis? [Y/y]'
  read(5,'(a1)') answ
  if(answ/='N' .or. answ/='n') then
     scale = rm/g(1)
     write(*,'(a,f10.4)') ' B-field scaling factor ', scale
     ! rescale
     p = scale**2 * p 
     g = scale* g
     psi = scale* psi
     gpp = scale**2 * gpp 
     xj = xj / scale
  endif
  
! Bfield**2 = (|grad_Psi_pol|**2 + g**2)/X**2
  do i=1,ns
     Bfield(1:nt1,i) = sqrt((gpp(1:nt1,i) + g(i)**2)/xa(1:nt1,i)**2)
  enddo


  ! save in map01.cdf file

  filename = 'map01.cdf'
  write(*,*)' '
  write(*,*)' Data will now be saved in file ', filename

  call cdfOpn(iu, filename, 'w')
  if(iu==0) then
     print*,'--Oh oh...-- failed to open ', filename
     ier = 1
  endif
  
  dims=(/1, 1, 1/)
  call cdfDefVar(iu, 'mth', dims, 'INT')
  call cdfDefVar(iu, 'nosurf', dims, 'INT')
  call cdfDefVar(iu, 'mx', dims, 'INT')
  call cdfDefVar(iu, 'npsi', dims, 'INT')
  call cdfDefVar(iu, 'kb', dims, 'INT')
  call cdfDefVar(iu, 'time', dims, 'R8')
!!$  call cdfDefVar(iu, 'r', dims, 'R8')
  call cdfDefVar(iu, 'xma', dims, 'R8')
  dims=(/80, 1, 1/)
  call cdfDefVar(iu, 'title', dims, 'CHAR')
  call cdfDefVar(iu, 'date', dims, 'CHAR')
  
  call cdfDefVar(iu, 'R0_major_axis', (/1,1,1/), 'R8')

  dims=(/nsi, 1, 1/)
  call cdfDefVar(iu, 'psi_pol_ori', dims, 'R8')
  call cdfDefVar(iu, 'psi_tor_ori', dims, 'R8')

  dims=(/ns, 1, 1/)
  call cdfDefVar(iu, 'psibar', dims, 'R8')
  call cdfDefVar(iu, 'psi_tor', dims, 'R8')
!!$  call cdfDefVar(iu, 'psival', dims, 'R8')
  call cdfDefVar(iu, 'p', dims, 'R8')
  call cdfDefVar(iu, 'g', dims, 'R8')
  call cdfDefVar(iu, 'q', dims, 'R8')
  dims=(/nt1, ns, 1/)
  call cdfDefVar(iu, 'x', dims, 'R8')
  call cdfDefVar(iu, 'z', dims, 'R8')
  call cdfDefVar(iu, 'xjacob', dims, 'R8')
  call cdfDefVar(iu, 'Grad_psi_square', dims, 'R8')
  call cdfDefVar(iu, 'Grad_theta_square', dims, 'R8')
  call cdfDefVar(iu, 'Grad_theta_psi', dims, 'R8')
  call cdfDefVar(iu, 'nufunc', dims, 'R8')
  call cdfDefVar(iu, 'Bfield', dims, 'R8')

  call cdfPutVar(iu, 'title', i2mex_o%label, ier)
  call cdfPutVar(iu, 'time', time, ier)
  call date_and_time(date=today, time=now)
  lon_date = today(1:4)//'/'//today(5:6)//'/'// &
       & today(7:8)//' at '//now(1:2)//':'//now(3:4)
  call cdfPutVar(iu, 'date', lon_date, ier)

  call cdfPutVar(iu, 'mth', nt1-1, ier)
  call cdfPutVar(iu, 'nosurf', ns, ier)
  call cdfPutVar(iu, 'mx', i2mex_mx, ier)
  call cdfPutVar(iu, 'npsi', i2mex_npsi, ier)
  call cdfPutVar(iu, 'kb', i2mex_kb, ier)

  call cdfPutVar(iu, 'psi_pol_ori', psii, ier)
  call cdfPutVar(iu, 'psi_tor_ori', torfluxi, ier)

!!$  call cdfPutVar(iu, 'r', 1.0_r8, ier) ! we don't rescale so r=1?
  call cdfPutVar(iu, 'R0_major_axis', rm, ier)
  ! save both: pol flux in Wb/rad=psibar and in Wb=psival to avoid confusion
  call cdfPutVar(iu, 'psibar', psi, ier)       ! in [Wb/rd]
!!$  call cdfPutVar(iu, 'psival', psi*i2mex_twopi_r8, ier) ! in [Wb]
  call cdfPutVar(iu, 'psi_tor', torflux, ier)

  
  call cdfPutVar(iu, 'p', p, ier) 
  call cdfPutVar(iu, 'g', g, ier) 
  call cdfPutVar(iu, 'q', q, ier) 
  temp = 0.0_r8
  temp(1:nt1, 1:ns) = xa
  call cdfPutVar(iu, 'x', temp, ier)
  temp(1:nt1, 1:ns) = za
  call cdfPutVar(iu, 'z', temp, ier)
  temp(1:nt1, 1:ns) = xj
  call cdfPutVar(iu, 'xjacob', temp, ier)
  temp(1:nt1, 1:ns) = gpp
  call cdfPutVar(iu, 'Grad_psi_square', temp, ier)
  temp(1:nt1, 1:ns) = gtt
  call cdfPutVar(iu, 'Grad_theta_square', temp, ier)
  temp(1:nt1, 1:ns) = gtp
  call cdfPutVar(iu, 'Grad_theta_psi', temp, ier)
  temp(1:nt1, 1:ns) = nufunc
  call cdfPutVar(iu, 'nufunc', temp, ier)
  temp(1:nt1, 1:ns) = Bfield
  call cdfPutVar(iu, 'Bfield', temp, ier)

  call cdfCls(iu)
  
  DEALLOCATE(t)
  DEALLOCATE(psii)
  DEALLOCATE(psi, p, q, g, torfluxi, torflux)
  DEALLOCATE(xa, za, xj, delta, Bfield, nufunc)
  DEALLOCATE(gtt, gtp, gpp)
  DEALLOCATE(gserror)
  DEALLOCATE(temp)

  call i2mex_free(ier)

  print*,'All done.'
  stop
end program 

subroutine pstCheckMercier(ier)

  ! Check if Mercier index is > 0 on a rational surface and return
  ! error flag ier > 0 if there is such a surface.

  use i2mex_mod
!!$  use pstcom

  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)
  integer, parameter :: nsg1 = 11 ! 
  integer, parameter :: nmax=3! max toroidal number
  integer, intent(out) :: ier ! error flag

  integer :: ns, iok, mmin, mmax, m, i, j, n, nrats
  real(r8) qmin, qmax, zdi, zqrat, psileft, psiright
  real(r8), dimension(:), allocatable :: psi, q
  real(r8) psis(nsg1), zes(nsg1), zfs(nsg1), zhs(nsg1), psi_padded(nsg1+2)

  ier = 0
  iok = 0
  call i2mex_getOriNs(ns, iok)

  allocate(psi(ns), q(ns))

  call i2mex_getOriPsi(ns, psi, iok)
  call i2mex_getQ(ns, psi, q, iok)
  
  qmin = minval(q)
  qmax = maxval(q)
  mmin = ceiling(nmax * qmin)
  mmax = floor(nmax * qmax)

  i = 1
  psis = psi(ns)
  psis(1) = psi(1)
  ! nsg1 is the max number of singular surfaces in pstcom
     do m = mmin, mmax
        zqrat = real(m,r8)/real(nmax,r8)
        if(qmin < zqrat .and. zqrat < qmax) then
        psileft = psi(1)
        if(i >= 2) psileft = psis(i-1)
        psiright = psi(ns)
        call i2mex_findRationalSurface(m, nmax, &
             & psileft, psiright, psis(i), iok)
!!$        call i2mex_error(iok)
        if(iok == 0) then
           i = i + 1
           if(i > nsg1) exit
        endif
        endif
     enddo

     nrats = i - 1

  ! nrats no of rational surfaces found
  ! rational surfaces are psis(1)...psis(nrats)
  psi_padded(1) = psi(1)
  psi_padded(2:nrats+1)= psis(1:nrats)
  psi_padded(nrats+2) = psi(ns)
  ! first psi must be on axis
  call i2mex_getGlasserGreenJohnsonEFH(nrats+1, psi_padded, &
       & zes, zfs, zhs, iok)

  do j = 1, nrats
     zdi = zes(j+1) + zfs(j+1) + zhs(j+1) - 0.25_r8
     if(zdi > 0._r8) then
        write(*,'(a,e10.2,a,e10.2)') ' *** Mercier unstable D_I =',zdi, &
             & ' @ psi/psi_a = ', psis(j)/psi(ns)
        
        ier = ier + 1
     endif
  enddo

  if(allocated(psi)) deallocate(psi)
  if(allocated(q)) deallocate(q)

end subroutine pstCheckMercier

subroutine eqi_geq_psfind(psi_geq, q_geq, ngeq, xi, ns, fbdy, &
     & psi, phi, q, ier)
  !
  !  EFIT G-EQDSK post-processing routine
  !
  !  find ns psi surfaces which are equally spaced in sqrt toroidal flux,
  !  and which span from psi(mag. axis) to
  !                      fbdy*(psi(separatrix)-psi(mag. axis)), where
  !  `fbdy' is a parameter, slightly less than one, which specifies how
  !  far to pull back from the separatrix (last closed flux surface).
 
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
 
  ! ------------------ arguments ---------------------------------
 
  integer, intent(in) :: ngeq        ! no. of pts in G-EQDSK 1d profiles
  real(r8), intent(in) :: psi_geq(ngeq)  ! G-EQDSK psi grid (strict ascending)
  real(r8), intent(in) :: q_geq(ngeq)    ! G-EQDSK q(psi)
 
  integer, intent(in) :: ns          ! no. of surfaces wanted
  real(r8), intent(in) :: xi(ns)     ! spacing of surfaces in sqrt(tor.flux)
                                     ! strict ascending, xi(1)=0, xi(ns)=1
 
  real(r8), intent(in) :: fbdy       ! back-off parameter for outermost surface
 
                                     ! on output psi(1)=0, psi(ns)=
                                     !   fbdy*(psi_geq(ngeq)-psi_geq(1))
 
  real(r8), intent(out) :: psi(ns)   ! the selected psi surfaces
  real(r8), intent(out) :: phi(ns)   ! the enclosed toroidal flux (Wb)
                                     !   satisfying xi(j)=sqrt(phi(j)/phi(ns))
  real(r8), intent(out) :: q(ns)     ! the q profile at {psi(1),...,psi(ns)}
 
  integer, intent(out) :: ier        ! completion code, 0=OK
 
  ! ------------------ local variables and arrays ----------------
 
  integer :: intmp,ixtra(ngeq)
  real(r8), dimension(:), allocatable :: psi_tmp,q_tmp
  !  psi grid (a) truncated for bdy adjustment & (b) with extra points
  !  to resolve sqrt(flux) near axis.

  real(r8) psi_bdy,r8ns,zdpsi
 
  real(r8) xi2(ns)                   ! xi squared ~ normalized tor. flux
  real(r8) tflux
 
  real(r8) zdum,tol,twopi
 
  integer luntty                     ! for messages
 
  integer i,j,ii,iglb,ibrk
 
  ! ------------------ executable code ---------------------------
 
  ! sanity checks
 
  ier = 0
 
  if((fbdy.lt.0.95_R8).or.(fbdy.gt.1.00_R8)) then
     call eqi_fromg_err( &
          ' ?eqi_geq_psfind:  fbdy not between 0.95 and 1.00')
     ier=120
  endif
 
  if((xi(1).ne.0.0_R8).or.(abs(xi(ns)-1.0_R8).gt.0.000001_R8)) then
     call eqi_fromg_err( &
          ' ?eqi_geq_psfind:  xi does not range from 0.00 to 1.00')
     ier=120
  endif
 
  do i=1,ngeq-1
     if(psi_geq(i+1).le.psi_geq(i)) then
        call eqi_fromg_err( &
             & ' ?eqi_geq_psfind:  G-EQDSK psi not in ascending order.')
        ier=1
        exit
     endif
  enddo
 
  do i=1,ns-1
     if(xi(i+1).le.xi(i)) then
        call eqi_fromg_err( &
             & ' ?eqi_geq_psfind:  xi coordinate not in ascending order.')
        ier=1
        exit
     endif
  enddo
 
  if(ier.ne.0) then
     ier=120
     return
  endif
 
  !
  !  ok......... last G-EQDSK psi value wanted:
 
  psi_bdy=psi_geq(1)+fbdy*(psi_geq(ngeq)-psi_geq(1))
 
  ! find GLB of selected psi_bdy in orig. G-EQDSK grid.
 
  do i=1,ngeq
     if(psi_geq(i).lt.psi_bdy) iglb=i
  enddo
 
  ! from this construct a temporary grid which is the same as the
  ! original, except in the edge region
 
  ibrk=max(1,(iglb-1))
  ixtra(ibrk:ngeq)=0
  intmp=ngeq
  if(ibrk.gt.1) then
     r8ns = ns
     ixtra(1) = sqrt(r8ns) - 1   ! pack extra pts near axis due to sqrt
     intmp = intmp + ixtra(1)    ! singularity in resolution
     do i=2,ibrk-1           
        ixtra(i)=ixtra(i-1)/2
        intmp = intmp + ixtra(i)
     enddo
  endif
 
  allocate(psi_tmp(intmp),q_tmp(intmp))

  ii = 0
  do i=1,ngeq
     if(i.le.ibrk) then
        ii=ii+1
        psi_tmp(ii)=psi_geq(i)
        if(ixtra(i).gt.0) then
           !  pack in extrap points
           zdpsi=(psi_geq(i+1)-psi_geq(i))/(ixtra(i)+1)
           do j=1,ixtra(i)
              ii=ii+1
              psi_tmp(ii)=psi_tmp(ii-1)+zdpsi
           enddo
        endif
     else
        ii=ii+1
        psi_tmp(ii)=psi_geq(ibrk) + &
             (i-ibrk)*(psi_bdy-psi_geq(ibrk))/(ngeq-ibrk)
     endif
  enddo
 
  ! interpolate q to this grid
 
  call eqi_geq_intrp(psi_geq, q_geq, ngeq, psi_tmp, intmp, q_tmp, ier)
  if(ier.ne.0) return
 
  !
  ! using relation dphi = q*dpsi, find the desired surfaces
  !
 
  xi2 = xi*xi
 
  tol=1.0e-12_R8
  zdum=0.0_R8
  call xdistrib(2, intmp,psi_tmp,q_tmp, 0,zdum,0,zdum, &
       & ns,xi2,psi, tflux, tol, ier)
  if(ier.ne.0) then
     ier=120
     return
  endif
 
  deallocate(q_tmp,psi_tmp)

  !
  ! ok get the q profile on this grid
  !
 
  call eqi_geq_intrp(psi_geq, q_geq, ngeq, psi, ns, q, ier)
  if(ier.ne.0) return
 
  !
  ! and the toroidal flux...  Wb, not Wb/rad -> factor of 2pi
  !
 
  twopi=6.2831853071795862_R8
  phi = twopi*tflux*xi2
 
  return
 
end subroutine eqi_geq_psfind

!  this file contains work routines for xplasma_integ_work.f90...
!  mostly f77 style interface, object pointer in eqi_intg_module

subroutine xplasma_integ_wk1(icode,id_R,id_Z,id_Bmod,id_BR,id_BZ, &
     nvec,thvec,rhovec,ansvec,ier)

  use xplasma_definitions
  use eqi_intg_module
  implicit NONE

  integer, intent(in) :: icode                ! code: what to evaluate
  integer, intent(in) :: id_R,id_Z,id_Bmod,id_BR,id_BZ   ! xplasma profile ids

  integer, intent(in) :: nvec                 ! vector size
  real*8, intent(in) :: rhovec(nvec),thvec(nvec) ! input vectors

  real*8, intent(out) :: ansvec(nvec)         ! result...
  integer, intent(out) :: ier                 ! status, 0=OK

  !-------------------------------
  integer :: id,ideriv_rho,ideriv_th
  !-------------------------------

  call xplasma_integ_gid(icode,id_R,id_Z,id_Bmod,id_BR,id_BZ, &
       id,ideriv_th,ideriv_rho,ier)

  if(ier.ne.0) return

  call xplasma_eval_prof(sp,id, &
       xplasma_theta_coord,thvec,xplasma_rho_coord,rhovec, &
       ansvec,ier, ideriv1=ideriv_th,ideriv2=ideriv_rho)

end subroutine xplasma_integ_wk1

subroutine xplasma_integ_wk2(icode,id_R,id_Z,id_Bmod,id_BR,id_BZ, &
     nvec,th_intrp,rho_intrp,ansvec,ier)

  use xplasma_definitions
  use eqi_intg_module
  implicit NONE

  integer, intent(in) :: icode                ! code: what to evaluate
  integer, intent(in) :: id_R,id_Z,id_Bmod,id_BR,id_BZ   ! xplasma profile ids

  integer, intent(in) :: nvec                 ! vector size
  type (xpeval) :: th_intrp,rho_intrp         ! pre-evaluated lookup info

  real*8, intent(out) :: ansvec(nvec)         ! result...
  integer, intent(out) :: ier                 ! status, 0=OK

  !-------------------------------
  integer :: id,ideriv_rho,ideriv_th
  !-------------------------------

  call xplasma_integ_gid(icode,id_R,id_Z,id_Bmod,id_BR,id_BZ, &
       id,ideriv_th,ideriv_rho,ier)

  if(ier.ne.0) return

  call xplasma_eval_prof(sp,id, &
       th_intrp, rho_intrp, &
       ansvec,ier, ideriv1=ideriv_th,ideriv2=ideriv_rho)

end subroutine xplasma_integ_wk2

subroutine xplasma_integ_gid(icode,id_R,id_Z,id_Bmod,id_BR,id_BZ, &
     id,ideriv_th,ideriv_rho,ier)

  !  choose profile id...

  use xplasma_definitions
  use eqi_intg_module
  implicit NONE

  integer, intent(in) :: icode                ! code: what to evaluate
  integer, intent(in) :: id_R,id_Z,id_Bmod,id_BR,id_BZ   ! xplasma profile ids

  integer, intent(out) :: id                  ! id chosen from above set
  integer, intent(out) :: ideriv_th           ! theta derivative control
  integer, intent(out) :: ideriv_rho          ! rho derivative control
  integer, intent(out) :: ier                 ! status, 0=OK

  !------------------------------------------

  id = 0
  ier = 0
  ideriv_th = 0
  ideriv_rho = 0

  if(icode.gt.2*eqi_intg_bign) then

     id = icode - 2*eqi_intg_bign

  else if(icode.gt.eqi_intg_bign) then

     id = icode - eqi_intg_bign

  else if(icode.eq.1) then
     id = id_R

  else if(icode.eq.2) then
     id = id_R
     ideriv_th=1

  else if(icode.eq.3) then
     id = id_R
     ideriv_rho=1

  else if(icode.eq.4) then
     id = id_Z

  else if(icode.eq.5) then
     id = id_Z
     ideriv_th=1

  else if(icode.eq.6) then
     id = id_Z
     ideriv_rho=1

  else if(icode.eq.7) then
     id = id_BMOD

  else if(icode.eq.8) then
     id = id_BR

  else if(icode.eq.9) then
     id = id_BZ

  else

     ier=611
     call xplasma_errmsg_append(sp, &
          ' invalid integrand term id code passed to xplasma_integ_wk1')
  endif

  if(id.eq.0) then
     ier=611
     if(icode.gt.6) then
        call xplasma_errmsg_append(sp, &
             ' integrand term eval failed -- B fields not available.')
     else
        call xplasma_errmsg_append(sp, &
             ' integrand term eval failed -- (R,Z) geometry not available.')
     endif
  endif

end subroutine xplasma_integ_gid

subroutine xplasma_integ_fwk(icode,nth,wk,local_wk,fans,ier)
  
  !  evaluate integrand f -- f*dth to be evaluated

  use xplasma_definitions
  use eqi_intg_module
  implicit NONE

  integer, intent(in) :: icode     ! integration selection code

  integer, intent(in) :: nth       ! number of values to evaluate & array dim.

  real*8, intent(in) :: wk(nth,eqi_intg_nfi)  ! cached data
  real*8, intent(in) :: local_wk(nth)      ! cached data
  !  R, dR/dtheta, ...

  real*8, intent(out) :: fans(nth)   ! integrand vector returned
  integer, intent(out) :: ier      ! set if icode not recognized...

  !--------------------------
  integer :: ith
  real*8 :: zmui
  real*8, parameter :: c2pi = 6.2831853071795862D+00
  real*8, parameter :: cmu0 = 2*c2pi*1.0d-7  ! 4*pi*1.0d-7
  !--------------------------

  ier=0
  if(icode.gt.2*eqi_intg_bign) then
     do ith=1,nth
        fans(ith)=c2pi*local_wk(ith)
     enddo
     
  else if(icode.eq.eqi_intg_vol) then
     do ith=1,nth
        fans(ith)=-c2pi*wk(ith,eqi_intg_r)* &
             wk(ith,eqi_intg_drdth)*wk(ith,eqi_intg_z)
     enddo

  else if(icode.eq.eqi_intg_area) then
     do ith=1,nth
        fans(ith)=-wk(ith,eqi_intg_drdth)*wk(ith,eqi_intg_z)
     enddo

  else if(icode.eq.eqi_intg_dvdrho) then
     do ith=1,nth
        fans(ith)=c2pi*wk(ith,eqi_intg_r)*abs( &
             wk(ith,eqi_intg_drdth)*wk(ith,eqi_intg_dzdrho) - &
             wk(ith,eqi_intg_dzdth)*wk(ith,eqi_intg_drdrho))
     enddo

  else if(icode.eq.eqi_intg_itor) then
     zmui=1.0d0/cmu0
     do ith=1,nth
        fans(ith)=abs(wk(ith,eqi_intg_drdth)*wk(ith,eqi_intg_br) + &
             wk(ith,eqi_intg_dzdth)*wk(ith,eqi_intg_bz))*zmui
        !! fans(ith)=sqrt((wk(ith,eqi_intg_drdth)**2+wk(ith,eqi_intg_dzdth)**2)* &
        !!     (wk(ith,eqi_intg_br)**2+wk(ith,eqi_intg_bz)**2))*zmui
     enddo

  else if(icode.eq.eqi_intg_lpol) then
     do ith=1,nth
        fans(ith)=sqrt(wk(ith,eqi_intg_drdth)**2+wk(ith,eqi_intg_dzdth)**2)
     enddo

  else if(icode.eq.eqi_intg_surf) then
     do ith=1,nth
        fans(ith)=sqrt(wk(ith,eqi_intg_drdth)**2+wk(ith,eqi_intg_dzdth)**2)* &
             c2pi*wk(ith,eqi_intg_r)
     enddo

  else if(icode.eq.eqi_intg_booz_qovg) then
     do ith=1,nth
        fans(ith)=sqrt(wk(ith,eqi_intg_drdth)**2+wk(ith,eqi_intg_dzdth)**2)/ &
             sqrt(wk(ith,eqi_intg_br)**2+wk(ith,eqi_intg_bz)**2)/ &
             (c2pi*wk(ith,eqi_intg_r)**2)
     enddo

  else
     fans=0
     ier=611
     call xplasma_errmsg_append(sp, &
          ' ?invalid integration code passed to xplasma_integ_fwk')
  endif

end subroutine xplasma_integ_fwk


subroutine xplasma_integ_fwkdv(icode,nth_eval,nth_dim,wk,local_wk, &
     bmax,fans,ier)
  
  !  evaluate integrand f*dV/drho -- f*dV/drho*dth to be evaluated

  use xplasma_definitions
  use eqi_intg_module
  implicit NONE

  integer, intent(in) :: icode     ! integration selection code

  integer, intent(in) :: nth_eval  ! number of theta values to evaluate
  integer, intent(in) :: nth_dim   ! number of theta values (array dim.)

  real*8, intent(in) :: wk(nth_dim,eqi_intg_nfi)  ! cached data
  real*8, intent(in) :: local_wk(nth_dim)

  !  R, dR/dtheta, ...

  real*8, intent(inout) :: bmax    ! max(B) on surface, used by some integrands

  real*8, intent(out) :: fans(nth_eval) ! integrand vector returned
  integer, intent(out) :: ier      ! set if icode not recognized...

  !--------------------------
  integer :: ith,idv
  real*8 :: zval,zR,zB,zdl,zd2,zgx,zh,zarg
  real*8, parameter :: c2pi = 6.2831853071795862D+00
  real*8, parameter :: c1 = 1.0d0
  real*8, parameter :: chalf = 0.5d0
  real*8, parameter :: c0 = 0.0d0
  real*8, parameter :: ceps6 = 1.0d-6
  !--------------------------

  idv = eqi_intg_2pirdetj
  ier=0

  if(icode.gt.eqi_intg_bign) then
     do ith=1,nth_eval
        zval=local_wk(ith)
        fans(ith)=wk(ith,idv)*zval
     enddo
  
  else if(icode.eq.eqi_intg_r2i) then
     do ith=1,nth_eval
        zval=wk(ith,eqi_intg_r)
        fans(ith)=wk(ith,idv)/(zval*zval)
     enddo

  else if(icode.eq.eqi_intg_ri) then
     do ith=1,nth_eval
        zval=wk(ith,eqi_intg_r)
        fans(ith)=wk(ith,idv)/zval
     enddo

  else if(icode.eq.eqi_intg_r3i) then
     do ith=1,nth_eval
        zval=wk(ith,eqi_intg_r)
        fans(ith)=wk(ith,idv)/(zval*zval*zval)
     enddo

  else if(icode.eq.eqi_intg_r1) then
     do ith=1,nth_eval
        zval=wk(ith,eqi_intg_r)
        fans(ith)=wk(ith,idv)*zval
     enddo

  else if(icode.eq.eqi_intg_r2) then
     do ith=1,nth_eval
        zval=wk(ith,eqi_intg_r)
        fans(ith)=wk(ith,idv)*zval*zval
     enddo

  else if(icode.eq.eqi_intg_br2) then
     do ith=1,nth_eval
        zval=wk(ith,eqi_intg_r)
        fans(ith)=wk(ith,idv)*wk(ith,eqi_intg_modb)*zval*zval
     enddo

  else if(icode.eq.eqi_intg_b2i) then
     do ith=1,nth_eval
        zval=wk(ith,eqi_intg_modb)
        fans(ith)=wk(ith,idv)/(zval*zval)
     enddo

  else if(icode.eq.eqi_intg_bi) then
     do ith=1,nth_eval
        zval=wk(ith,eqi_intg_modb)
        fans(ith)=wk(ith,idv)/zval
     enddo

  else if(icode.eq.eqi_intg_b1) then
     do ith=1,nth_eval
        zval=wk(ith,eqi_intg_modb)
        fans(ith)=wk(ith,idv)*zval
     enddo

  else if(icode.eq.eqi_intg_b2) then
     do ith=1,nth_eval
        zval=wk(ith,eqi_intg_modb)
        fans(ith)=wk(ith,idv)*zval*zval
     enddo

  else if(icode.eq.eqi_intg_bz2) then
     do ith=1,nth_eval
        zval=wk(ith,eqi_intg_bz)
        fans(ith)=wk(ith,idv)*zval*zval
     enddo

  else if(icode.eq.eqi_intg_g2b2i) then
     do ith=1,nth_eval
        zdl = sqrt(wk(ith,eqi_intg_drdth)**2+wk(ith,eqi_intg_dzdth)**2)
        zd2 = abs(wk(ith,eqi_intg_drdth)*wk(ith,eqi_intg_dzdrho) - &
             wk(ith,eqi_intg_dzdth)*wk(ith,eqi_intg_drdrho))
        zgx = zdl/zd2
        zB = wk(ith,eqi_intg_modb)
        fans(ith)=wk(ith,idv)*(zgx*zgx)/(zB*zB)
     enddo

  else if(icode.eq.eqi_intg_g2r2i) then
     do ith=1,nth_eval
        zdl = sqrt(wk(ith,eqi_intg_drdth)**2+wk(ith,eqi_intg_dzdth)**2)
        zd2 = abs(wk(ith,eqi_intg_drdth)*wk(ith,eqi_intg_dzdrho) - &
             wk(ith,eqi_intg_dzdth)*wk(ith,eqi_intg_drdrho))
        zgx = zdl/zd2
        zR = wk(ith,eqi_intg_r)
        fans(ith)=wk(ith,idv)*(zgx*zgx)/(zR*zR)
     enddo

  else if(icode.eq.eqi_intg_g2r3i) then
     do ith=1,nth_eval
        zdl = sqrt(wk(ith,eqi_intg_drdth)**2+wk(ith,eqi_intg_dzdth)**2)
        zd2 = abs(wk(ith,eqi_intg_drdth)*wk(ith,eqi_intg_dzdrho) - &
             wk(ith,eqi_intg_dzdth)*wk(ith,eqi_intg_drdrho))
        zgx = zdl/zd2
        zR = wk(ith,eqi_intg_r)
        fans(ith)=wk(ith,idv)*(zgx*zgx)/(zR*zR*zR)
     enddo

  else if(icode.eq.eqi_intg_g2r2) then
     do ith=1,nth_eval
        zdl = sqrt(wk(ith,eqi_intg_drdth)**2+wk(ith,eqi_intg_dzdth)**2)
        zd2 = abs(wk(ith,eqi_intg_drdth)*wk(ith,eqi_intg_dzdrho) - &
             wk(ith,eqi_intg_dzdth)*wk(ith,eqi_intg_drdrho))
        zgx = zdl/zd2
        zR = wk(ith,eqi_intg_r)
        fans(ith)=wk(ith,idv)*(zgx*zgx*zR*zR)
     enddo

  else if(icode.eq.eqi_intg_giri) then
     do ith=1,nth_eval
        zdl = sqrt(wk(ith,eqi_intg_drdth)**2+wk(ith,eqi_intg_dzdth)**2)
        zd2 = abs(wk(ith,eqi_intg_drdth)*wk(ith,eqi_intg_dzdrho) - &
             wk(ith,eqi_intg_dzdth)*wk(ith,eqi_intg_drdrho))
        zgx = zdl/zd2
        zR = wk(ith,eqi_intg_r)
        fans(ith)=wk(ith,idv)/(zgx*zR)
     enddo

  else if(icode.eq.eqi_intg_g2) then
     do ith=1,nth_eval
        zdl = sqrt(wk(ith,eqi_intg_drdth)**2+wk(ith,eqi_intg_dzdth)**2)
        zd2 = abs(wk(ith,eqi_intg_drdth)*wk(ith,eqi_intg_dzdrho) - &
             wk(ith,eqi_intg_dzdth)*wk(ith,eqi_intg_drdrho))
        zgx = zdl/zd2
        fans(ith)=wk(ith,idv)*(zgx*zgx)
     enddo

  else if(icode.eq.eqi_intg_g1) then
     do ith=1,nth_eval
        zdl = sqrt(wk(ith,eqi_intg_drdth)**2+wk(ith,eqi_intg_dzdth)**2)
        zd2 = abs(wk(ith,eqi_intg_drdth)*wk(ith,eqi_intg_dzdrho) - &
             wk(ith,eqi_intg_dzdth)*wk(ith,eqi_intg_drdrho))
        zgx = zdl/zd2
        fans(ith)=wk(ith,idv)*(zgx)
     enddo

  else if(icode.eq.eqi_intg_fhnclass) then
     if(bmax.eq.c0) then
        do ith=1,nth_eval
           bmax=max(bmax,wk(ith,eqi_intg_modb))
        enddo
     endif
     do ith=1,nth_eval
        zb=wk(ith,eqi_intg_modb)
        zh=zb/bmax              ! B/Bmax
        zarg=max((c1-zh),ceps6)
        fans(ith)=wk(ith,idv)*(c1-sqrt(zarg)*(c1+chalf*zh))/(zh*zh)
     enddo

  else if(icode.eq.eqi_intg_fhnc_tsc) then
     if(bmax.eq.c0) then
        do ith=1,nth_eval
           bmax=max(bmax,wk(ith,eqi_intg_modb))
        enddo
     endif
     do ith=1,nth_eval
        zb=wk(ith,eqi_intg_modb)
        zh=zb/bmax              ! B/Bmax
        zarg=max((c1-zh),ceps6)
        fans(ith)=wk(ith,idv)*(c1/zb**2)*sqrt(zarg)*(c1-zarg/3)
     enddo

  else if(icode.eq.eqi_intg_fhmmx) then
     if(bmax.eq.c0) then
        do ith=1,nth_eval
           bmax=max(bmax,wk(ith,eqi_intg_modb)) ! Bmax -- max on flux surface
        enddo
     endif
     do ith=1,nth_eval
        zb=wk(ith,eqi_intg_modb)
        zh=zb/bmax              ! B/Bmax
        zarg=max((c1-zh),ceps6)
        fans(ith)=wk(ith,idv)*sqrt(zarg)  ! sqrt(1-B/Bmax)
     enddo

  else
     fans=0
     ier=611
     call xplasma_errmsg_append(sp, &
          ' ?invalid integration code passed to xplasma_integ_fwkdv')
  endif

end subroutine xplasma_integ_fwkdv

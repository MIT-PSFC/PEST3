module xplasma_nclass

  !  code to compute NCLASS Pfirsch-Schlutter moments and related profiles
  !  ported from old_xplasma DMC Feb 2008

  use xplasma_obj
  use xplasma_flxint
  use xplasma_rzgeo
  use xplasma_ctran
  implicit NONE

  private

  public :: xplasma_psmom, xplasma_ftrap

  !--------------------------
  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: ZSMALL = 1.0d-10
  !--------------------------

  ! (originally from eq_psmom_mod...)

  ! -- parameters --
  integer, parameter :: r8 = kind(0.0D0)                   ! primary kind

  real(kind=r8), parameter :: znc_pi    = 3.1415926535897932_r8 ! pi

  real(kind=r8), parameter :: zrho_min  = 1.e-7_r8         ! rho not allowed below this for xplasma
  real(kind=r8), parameter :: zmax_theta_error = 1.e-3_r8  ! ztheta array must span 2pi within this tolerance
                                                           ! the boundary
  ! -- flags for flow control --
  integer, parameter :: itest=0              ! nonzero for testing

  ! -- xplasma signs --
  integer :: nsnccwb                         ! +1 if Btoroidal  is ccw, -1 o.w.
  integer :: nsnccwi                         ! +1 if plasma current  is ccw, -1 o.w.
  integer, parameter :: ichi_ccw = 1         ! +1 if chi increases ccw, -1 o.w. [ update Sept2006, assume ccw ]

  ! -- existing xplasma function ids --
  integer :: id_rho, id_chi            ! xplasma internal grid codes __RHO, __CHI
  integer :: id_br,  id_bz, id_bmod    ! xplasma internal b field codes BR, BZ, BMOD
  integer :: id_psi                    ! PSI

  ! -- rho,theta grids --
  integer :: nc_rhob                   ! size of xplasma boundary grid
  integer :: nc_theta                  ! size of poloidal grid
  integer :: nc_points                 ! number of rho points to evaluate 
                                             ! nclass quantities
  logical :: is_rhob                   ! .true. if the points and rhob arrays are 
                                             ! the same
  real(kind=r8) :: zrho_axis           ! rhos below this value are considered axis values

  real(kind=r8), allocatable, dimension(:) :: zrhob        ! boundary rho
  real(kind=r8), allocatable, dimension(:) :: ztheta       ! poloidal angle
  real(kind=r8), allocatable, dimension(:) :: ztheta_xp    ! xplasma poloidal angle (may be weird because of chi_cwdir change)
  real(kind=r8), allocatable, dimension(:) :: zpoints      ! rho points to use for nclass quantities

  ! -- moments --
  integer :: nc_mom                          ! number of moments desired for Fm

  ! -- new xplasma function ids --
  integer :: id_bigtheta                     ! xplasma code for computed BIGTHETA if is_rhob
  integer :: id_gamma                        ! xplasma code for computed GAMMA if is_rhob


  ! -- computed quantities --
  real(kind=r8), allocatable, dimension(:,:) :: zbigtheta  ! the BIG theta on the chi,points grid
  real(kind=r8), allocatable, dimension(:)   :: zgamma     ! gamma  at the points
  real(kind=r8), allocatable, dimension(:)   :: zdvol      ! dV/dxi at the points, at axis d2V/dxi2
  real(kind=r8), allocatable, dimension(:)   :: zdpsi      ! dpsi/dxi at the points, at axis d2psi/dxi2
  real(kind=r8), allocatable, dimension(:)   :: zngb2      ! <(n.grad(B))**2> at the points

  real(kind=r8), allocatable, dimension(:)   :: zb2        ! <B**2> at the points
  real(kind=r8), allocatable, dimension(:)   :: zbgrth     ! <B.grad(theta)> at the points
  real(kind=r8)                              :: zbmod_zero ! Bmod at axis

  real(kind=r8), allocatable, dimension(:,:) :: zmom_a     ! Am moment
  real(kind=r8), allocatable, dimension(:,:) :: zmom_b     ! Bm moment
  real(kind=r8), allocatable, dimension(:,:) :: zmom_c     ! Cm moment
  real(kind=r8), allocatable, dimension(:,:) :: zmom_d     ! Dm moment
  real(kind=r8), allocatable, dimension(:,:) :: zmom_f     ! Fm Pfirsch-Schluter moment (m,irho)

  !
  ! --- working variables ---
  ! iztotal may not be changed, the rest can be used in any subroutine
  !
  integer :: iztotal                              ! maximum number of points needed in 
                                                  ! xplasma call

  real(kind=r8), dimension(:),     allocatable :: zrho, zchi       ! for xplasma call
  real(kind=r8), dimension(:,:),   allocatable :: zans             ! hold B fields
  real(kind=r8), dimension(:,:,:), allocatable :: zgans            ! hold dB/dR,dB/dZ field
  real(kind=r8), allocatable :: zdet(:), ztens(:,:,:)              ! holds metrics

  real(kind=r8), dimension(:,:),   allocatable :: zintegrand1      ! integrand for defining THETA
  real(kind=r8), dimension(:,:),   allocatable :: zintegrand2      ! integrand for dV/dxi
  real(kind=r8), dimension(:,:),   allocatable :: zintegrand3      ! integrand for <(n.grad(B))**2>
  real(kind=r8), dimension(:,:),   allocatable :: zintegrand4      ! integrand for <B**2>
  real(kind=r8), dimension(:,:),   allocatable :: zresult1,zresult2,zresult3   ! surface integration 
  real(kind=r8), dimension(:,:),   allocatable :: zresult4                     ! surface integration

  real(kind=r8) :: br,bz,bmod                              ! B fields
  real(kind=r8) :: db_drho, db_dtheta                      ! dB/drho, dB/dtheta
  real(kind=r8) :: drho_dR, drho_dZ, dtheta_dR,dtheta_dZ   ! jacobian elements
  real(kind=r8) :: jacob                                   ! 2D jacobian
  real(kind=r8) :: db_dR, db_dZ                            ! dB/dR, dB/dZ

contains
  
  !
  ! ----------------- xplasma_ftrap ----------------
  ! Compute the trapping fraction from the integral
  !
  !      ftrap = 1 - (3/4)*<h^2>*integral_0to1[ lambda/<sqrt(1-lambda*h)>  dlambda]
  !          h = B/Bmax
  !
  ! The integral is done over alpha, where lambda=sin(pi/2 * alpha).  The alpha
  ! integrand has the form,
  !        lambda_integrand = lambda/<sqrt(1-lambda*h)>
  !        alpha_integrand  = (pi/2)*sqrt(1-lambda^2)*lambda_integrand
  !
  subroutine xplasma_ftrap(s, rho, ftrap, ier, rhomin, nsimp, integrand, lambda, b2)

    type (xplasma), pointer :: s   ! object for which to compute profiles...

    real*8,  intent(in)  :: rho(:)     ! radial coordinate
    real*8,  intent(out) :: ftrap(:)   ! trapping fraction
    integer, intent(out) :: ier        ! nonzero on error

    real*8,  intent(in),  optional :: rhomin          ! rho values below this are deemed on axis
    integer, intent(in),  optional :: nsimp           ! number of simpson intervals to use
    real*8,  intent(out), optional :: integrand(:,:)  ! lambda integrand at each simpson point ; nsimp must be defined
                                                      ! and this should have shape [2*nsimp+1, size(rho)]
    real*8,  intent(out), optional :: lambda(:)       ! when integrand is given, this is lambda at each simpson point
                                                      ! and this should have shape [2*nsimp+1]
    real*8,  intent(out), optional :: b2(:)           ! <B2> for checking flux surface integration [size(rho)]
 
    integer, parameter :: ipts = 5    ! number of gaussian points is 2*ipts between theta grid points using the xplasma standard

    integer :: iveca            ! number of output points
    integer :: isize2,isize2_1  ! theta grid size and interval size
    integer :: isinteg          ! number of quadrature points on one flux surface
    integer :: isall            ! total number of quadrature points on each flux surface
    integer :: i,j,k,ixp,ith,ii ! temps
    integer :: ksimp            ! number of simpson intervals for alpha integral
    integer :: ilamb            ! loop over lambda quadrature points
    integer :: maxsimp          ! totalnumber of simpson points

    real(r8) :: x5(ipts),w5(ipts)  ! gaussian points and weights
    real(r8) :: zxcen, zdx2        ! center and half-width of theta interval
    real(r8) :: dalpha             ! interval size in alpha
    real(r8) :: wlamb              ! lambda quadrature weight
    real(r8) :: xlamb              ! lambda value
    real(r8) :: nlamb              ! lambda quadrature numerator
    real(r8) :: zz                 ! temp

    real(r8), allocatable :: ztheta(:)   ! theta at interval boundaries
    real(r8), allocatable :: zth(:)      ! theta at all quadrature points
    real(r8), allocatable :: wth(:)      ! weight at all quadrature points
    real(r8), allocatable :: rhovec(:)   ! rho at all quadrature and radial points
    real(r8), allocatable :: thvec(:)    ! theta at all quadrature and radial points
    real(r8), allocatable :: bmax(:)     ! Bmax at radial points
    real(r8), allocatable :: zdet(:)     ! jacobian at all points
    real(r8), allocatable :: rvec(:)     ! R at all points
    real(r8), allocatable :: h(:)        ! B/Bmax at all points
    real(r8), allocatable :: h2(:)       ! dVol/dxi*<h^2> at radial points

    ier   = 0
    iveca = size(rho)

    call xplasma_global_info(s, ier, init_check=.TRUE.)
    if(ier.ne.0) then
       ier=777 ; return
    endif

    if(iveca<=0 .or. size(ftrap)/=iveca) then
       call xplasma_errmsg_append(s, '?xplasma_ftrap:"rho" dimension inconsistency')
       ier=9999 ;return
    endif

    if (present(nsimp)) then
       if (nsimp<1) then
          call xplasma_errmsg_append(s, '?xplasma_ftrap: bad nsimp argument')
          ier=1 ; return
       end if
       ksimp = nsimp
    else
       ksimp = 64   ! hard coded, this gives a max error around 1.e-6 at innermost flux surfaces
    end if

    if (present(integrand)) then
       if (.not. present(nsimp)) then
          call xplasma_errmsg_append(s, '?xplasma_ftrap: nsimp argument must be present with integrand argument')
          ier=1 ; return
       end if

       if (size(integrand,1)/=(2*nsimp+1) .or. size(integrand,2)/=iveca) then
          call xplasma_errmsg_append(s, '?xplasma_ftrap: integrand  dimension inconsistency')
          ier=9999 ; return
       endif

       if (present(lambda)) then
          if (size(lambda)/=(2*nsimp+1)) then
             call xplasma_errmsg_append(s, '?xplasma_ftrap: lambda dimension inconsistency')
             ier=9999 ;return
          end if
          lambda = 0.d0
       end if

       integrand = 0.d0
    end if          
          
    if (present(b2)) then
       if (size(b2)/=iveca) then
          call xplasma_errmsg_append(s, '?xplasma_ftrap: b2 dimension inconsistency')
          ier=9999 ;return
       end if
    endif

    ier    = 0
    isize2 = 0
    ftrap  = 0._r8

    call xplasma_integ_g10info(x5,w5)    ! grab gaussian quadrature data

    if (present(rhomin)) then
       zrho_axis = rhomin
    else
       zrho_axis = zrho_min
    end if

    zrho_axis = min(0.2_r8,max(zrho_min,zrho_axis))

    if (ier==0) call xnc_get_id(s, "__CHI", 0, id_chi,  ier) 
    if (ier==0) call xnc_get_id(s, "BMOD",  1, id_bmod, ier) 
    if (ier==0) call xplasma_grid_size(s,id_chi,isize2, ier)

    if (ier/=0 .or. isize2<=2) then
       call xplasma_errmsg_append(s, '? xplasma_ftrap: Unable to get some of the xplasma IDs')
       ier=1 ; return
    end if

    ! -- exit through 100 from now on --

    isize2_1 = isize2-1
    isinteg  = isize2_1*2*ipts
    isall    = isinteg*iveca

    allocate(ztheta(isize2),zth(isinteg),wth(isinteg))
    allocate(rhovec(isall), thvec(isall), bmax(iveca), zdet(isall), h(isall))
    allocate(h2(iveca), rvec(isall))

    ! -- get theta grid points --
    do i=1,isize2
       ztheta(i) = (-znc_pi) + real(i-1,r8)*2.0_r8*znc_pi/isize2_1
    end do

    ! -- get theta quadrature points, weights --
    ixp = -2*ipts
    do ith = 1, isize2_1   ! loop over theta grid intervals       
       ixp=ixp+2*ipts
       
       zxcen= (ztheta(ith+1)+ztheta(ith))/2._r8
       zdx2 = (ztheta(ith+1)-ztheta(ith))/2._r8
       
       do ii=1,ipts
          zth(ixp+ii)      = zxcen+x5(ii)*zdx2
          zth(ixp+ipts+ii) = zxcen-x5(ii)*zdx2
          wth(ixp+ii)      = w5(ii)*zdx2
          wth(ixp+ipts+ii) = w5(ii)*zdx2
       enddo
    end do
    
    ! -- get rho,theta at all points and Bmax on flux surface --
    ixp=0
    do k=1, iveca
       call xplasma_Bminmax(s, rho(k), ier, bmax=bmax(k))

       if (ier/=0) then
          call xplasma_errmsg_append(s, '? xplasma_ftrap: error finding Bmax on flux surface')
          ier=1 ; goto 100
       else if (bmax(k)<=0._r8) then
          call xplasma_errmsg_append(s, '? xplasma_ftrap: zero B field on flux surface')
          ier=1 ; goto 100
       end if

       do j = 1, isinteg
          ixp = ixp+1
          rhovec(ixp) = rho(k)
          thvec(ixp)  = zth(j)
       end do
    end do

    ! -- Jacobian at all points --
    call xplasma_rzjac(s, rhovec, thvec, ier, rzdetj=zdet, r=rvec)

    if (ier/=0) then
       call xplasma_errmsg_append(s, '? xplasma_ftrap: bad call to xplasma_rzjac')
       ier=1 ; goto 100
    end if

    zdet = rvec*abs(zdet)

    ! -- form h=B/Bmax at all points and get h2 --
    call xplasma_eval_prof(s, id_bmod, &
         xplasma_rho_coord,   rhovec, &
         xplasma_theta_coord, thvec, &
         h, ier)

    if (ier/=0) then
       call xplasma_errmsg_append(s, '? xplasma_ftrap: bad call to xplasma_eval_prof for Bmod')
       ier=1 ; goto 100
    end if

    ixp=0
    do k=1, iveca
       h2(k) = 0._r8
       do ii = 1, isinteg
          ixp = ixp+1
          h(ixp) = min(1._r8, h(ixp)/bmax(k))

          h2(k) = h2(k) + wth(ii) * zdet(ixp) * h(ixp)**2
       end do
    end do

    if (present(b2)) then
       do k=1, iveca
          if (rho(k)<=zrho_axis) then
             b2(k) = bmax(k)**2
          else
             ixp = (k-1)*isinteg
             zz  = sum(wth(1:isinteg)*zdet(ixp+1:ixp+isinteg))  ! flux surface integration

             if (zz<=0._r8) then
                call xplasma_errmsg_append(s, '? xplasma_ftrap: bad dV/dxi for <B2>')
                ier=1 ; goto 100
             end if

             b2(k) = bmax(k)**2 * h2(k)/zz
          end if
       end do
    end if

    ! -- now do integration over alpha, lambda=sin((pi/2)*alpha) --
    dalpha = 1._r8/(2*ksimp)

    maxsimp = 2*ksimp+1

    do ilamb = 2, maxsimp          ! loop over simpson points, lambda=0 has a zero integrand
       if (mod(ilamb,2)==0) then
          wlamb = 4._r8*dalpha/3._r8
       else 
          wlamb = 2._r8*dalpha/3._r8
       end if

       if (ilamb==maxsimp) then
          if (.not. present(integrand)) cycle  ! alpha integrand is zero in this case

          xlamb = 1._r8            ! because znc_pi isn't quite pi
          nlamb = 0._r8
          wlamb = 0._r8            ! integrand in lambda is nonzero so need to do this step
       else
          xlamb = sin(znc_pi*(ilamb-1)*dalpha/2._r8)
          nlamb = sin(znc_pi*(ilamb-1)*dalpha)
       end if


       if (present(integrand) .and. present(lambda)) lambda(ilamb)=xlamb

       do k=1, iveca               ! loop over each flux surface
          if (rho(k)<=zrho_axis) then
             ftrap(k) = 0._r8
             if (present(integrand)) integrand(ilamb,k)= xlamb/max(.01_r8,sqrt(1._r8-xlamb))
          else
             ixp = (k-1)*isinteg
             zz  = sum(wth(1:isinteg)*zdet(ixp+1:ixp+isinteg)*sqrt(1._r8-xlamb*h(ixp+1:ixp+isinteg)))  ! flux surface integration

             if (zz<=0._r8) then
                call xplasma_errmsg_append(s, '? xplasma_ftrap: bad denominator in lambda integrand')
                ier=1 ; goto 100
             end if
             
             ftrap(k) = ftrap(k) + wlamb*(nlamb/zz)                   ! add integrand in alpha
             
             if (present(integrand)) integrand(ilamb,k)=xlamb/zz      ! integrand in lambda
          end if
       end do
    end do

    where (rho>zrho_axis)
       ftrap = max(1._r8 - 0.75_r8*h2*(znc_pi/4._r8)*ftrap, 0._r8)
    else where
       ftrap = 0._r8
    end where

100 continue
    if (allocated(ztheta)) deallocate(ztheta, zth, wth)
    if (allocated(rhovec)) deallocate(rhovec, thvec, bmax, zdet, h)
    if (allocated(h2))     deallocate(h2, rvec)

  end subroutine xplasma_ftrap

  !------------------------------------------
  !  THIS the only PUBLIC entry point...

  subroutine xplasma_psmom(s, &
       rho, rhomin, psmom, gamma, ngrdb2, b2, zerror, ier)

    type (xplasma), pointer :: s   ! object for which to compute profiles...

    real*8,  intent(in), dimension(:) :: rho        ! radial coordinate
    real*8,  intent(in) :: rhomin  ! rho values below this are deemed on axis

    real*8,  intent(out), dimension(:,:) :: psmom   ! PS moments
    !  FIRST dimension: moment index; SECOND dimension: radial coordinate
    real*8,  intent(out), dimension(:) :: gamma     ! wayne's gamma 
         ! = 2*pi/integral_0_2pi[Bmod/B.grad(theta)]_dtheta
    real*8,  intent(out), dimension(:) :: ngrdb2    ! <(n.grad(B))**2>
    real*8,  intent(out), dimension(:) :: b2        ! <B**2>
    real*8,  intent(out), dimension(:) :: zerror    ! convergence estimate
         ! |sum(Fm) - <(n.grad(B))**2>/<B**2>| / max(<(n.grad(B))**2>/<B**2>)

    integer, intent(out) :: ier                ! nonzero on error
    ! ier=777 means the passed xplasma pointer is invalid or the object
    !         has not been initialized.

    ! (rho) vector must be monotonic incureasing; 
    !       size of all 1d outputs must match size(rho)

    !       size(psmom,2) = size(rho) also expected.

    !----------------------
    real*8 :: yr02                             ! max(<(n.grad(B))**2>/<B**2>)
    integer :: i,iveca,nmom
    !----------------------

    ier = 0
    iveca=size(rho)
    nmom =size(psmom,1)
    if((iveca.le.0).or.(nmom.le.0)) then
       ier=9999
       call xplasma_errmsg_append(s, &
            'xplasma_psmom: "rho" or "psmom" has empty dimension.')
       return
    endif

    call xplasma_global_info(s,ier, init_check=.TRUE.)
    if(ier.ne.0) then
       ier=777
       return
    endif

    if(size(psmom,2).ne.iveca) ier=ier+1
    if(size(gamma).ne.iveca) ier=ier+1
    if(size(ngrdb2).ne.iveca) ier=ier+1
    if(size(b2).ne.iveca) ier=ier+1
    if(size(zerror).ne.iveca) ier=ier+1

    if(ier.ne.0) then
       ier=9999
       call xplasma_errmsg_append(s, &
            ' "rho" dimension inconsistency in xplasma_psmom')
       return
    endif

    ier=0
    do
       call xnc_get_grid(s, rhomin, ier)                  ; if (ier/=0) exit
       call xnc_set_points(iveca, rho, nmom, ier)         ; if (ier/=0) exit
       call xnc_allocate
  
       psmom = ZERO
       gamma = ZERO
       ngrdb2 = ZERO
       b2 = ZERO
       zerror = ZERO

       call xnc_get_bigtheta(s,ier)  ; if (ier/=0) exit
       call xnc_get_abcd(s,ier)      ; if (ier/=0) exit

       psmom = zmom_f
  
       gamma = zgamma
       ngrdb2 = zngb2
       b2 = zb2

       yr02 = max(1.D-20, maxval(zngb2(1:iveca)/zb2(1:iveca)))
       do i = 1, iveca
          zerror(i) = abs(sum(psmom(1:nmom,i)) - zngb2(i)/zb2(i))/yr02
       end do

       exit
    enddo

    if(ier.ne.0) then
       ier=9999
    endif

    call xnc_free

  end subroutine xplasma_psmom

  !========================================================================
  ! all the rest of the code is PRIVATE...
  !========================================================================

  ! ------------------------- xnc_free ------------------------------------
  ! free up the allocated memory
  !
  subroutine xnc_free
    if(allocated(zbigtheta)) &
         deallocate(zbigtheta,zgamma,zdvol,zdpsi,zngb2,zb2,zbgrth)

    if(allocated(zmom_a)) &
         deallocate(zmom_a, zmom_b, zmom_c, zmom_d, zmom_f)

    if(allocated(zrho)) &
         deallocate(zrho, zchi, zdet)

    if(allocated(zans)) &
         deallocate(zans, zgans, ztens)

    if(allocated(zintegrand1)) &
         deallocate(zintegrand1,zintegrand2,zintegrand3,zintegrand4)

    if(allocated(zresult1)) &
         deallocate(zresult1,zresult2,zresult3,zresult4)

    if(allocated(zrhob)) &
         deallocate(zrhob,ztheta,ztheta_xp)

    if(allocated(zpoints)) deallocate(zpoints)

  end subroutine xnc_free

  ! ------------------------- xnc_get_grid ------------------------------------
  ! get size nc_rhob, nc_theta and data zrhob,ztheta for xplasma grids and 
  ! get the basic xplasma profile ids
  !
  subroutine xnc_get_grid(s,rhomin, ier)

    type (xplasma), pointer :: s   ! object for which to compute profiles...

    real*8,  intent(in)   :: rhomin    ! rhos below this should be considered axis values
    integer, intent(out)  :: ier       ! nonzero on error

    !--------------------

    integer :: isize1, isize2   ! temporary array size
    integer :: ierr1, ierr2     ! error codes
    integer :: i                ! loop variable

    !--------------------

    ier=0
    zrho_axis = rhomin

    !
    ! --- toroidal signs ---
    ! 

    call xplasma_global_info(s,ier, bphi_ccw=nsnccwb, jphi_ccw=nsnccwi)
    if (nsnccwb==0 .or. nsnccwi==0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_grid: xplasma not initialized -- nsnccwb=0')
       goto 100
    end if

    !
    ! --- get ids ---
    !
    ierr1=0
    call xnc_get_id(s, "__RHO", 0, id_rho,  ierr2) ; ierr1=ierr1+ierr2
    call xnc_get_id(s, "__CHI", 0, id_chi,  ierr2) ; ierr1=ierr1+ierr2
    call xnc_get_id(s, "BR",    1, id_br,   ierr2) ; ierr1=ierr1+ierr2
    call xnc_get_id(s, "BZ",    1, id_bz,   ierr2) ; ierr1=ierr1+ierr2
    call xnc_get_id(s, "BMOD",  1, id_bmod, ierr2) ; ierr1=ierr1+ierr2
    call xnc_get_id(s, "PSI",   1, id_psi,  ierr2) ; ierr1=ierr1+ierr2

    if (ierr1>0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_grid: Unable to get some of the xplasma IDs')
       goto 100
    end if

    call xplasma_grid_size(s,id_rho,isize1,ierr1)
    call xplasma_grid_size(s,id_chi,isize2,ierr2)

    if( isize1*isize2 > 0 ) then
       nc_rhob=isize1
       nc_theta=isize2

       allocate(zrhob(nc_rhob),ztheta(nc_theta),ztheta_xp(nc_theta))

       call xplasma_grid(s,id_rho,zrhob,ierr1)
       call xplasma_grid(s,id_chi,ztheta_xp,ierr2)

       if (ierr1/=0 .or. ierr2/=0) then
          call xplasma_errmsg_append(s, &
               '? xnc_get_grid: Error fetching grid data')
          goto 100
       end if

       do i=1,isize2
          ztheta(i) = (-znc_pi) + (i-1)*2.0_r8*znc_pi/max(1,isize2-1)
       end do

    else
       call xplasma_errmsg_append(s,'?xnc_get_grid: Zero size grids')
       goto 100
    end if
    return

100 continue
    ier=1
    return
  end subroutine xnc_get_grid

  !
  ! ------------------------- xnc_set_points ------------------------------------
  ! setup points over which the moments will be computed
  !
  subroutine xnc_set_points(npoints, points, nmom, ier)
    integer,                    intent(in) :: npoints            ! number of rho points
    real(kind=r8),dimension(*), intent(in) :: points             ! rho
    integer,                    intent(in) :: nmom               ! number of moments for Fm to calculate
    integer,                    intent(out):: ier                ! nonzero on error

    !-------------------------
    integer        :: i               ! loop variable
    !-------------------------

    ier = 0

    ! -- store --
    nc_mom = nmom

    allocate(zpoints(npoints))
    nc_points = npoints
    do i = 1, npoints
       zpoints(i) = points(i)
       if (points(i)<zrho_axis) zpoints(i)=0._r8
    end do

  end subroutine xnc_set_points
  !
  ! ------------------------- xnc_allocate ------------------------------------
  !
  subroutine xnc_allocate
    !
    ! allocate arrays -- call after xnc_get_grid, xnc_set_points
    !
    integer :: i           ! loop variable
    integer :: inc_points  ! number of allocated points in zpoints

    ! -- points array the same as rhob? --
    is_rhob = (nc_points == nc_rhob)
    if (is_rhob) then
       do i = 1, nc_points
          if (abs(zrhob(i)-zpoints(i))>ZSMALL) is_rhob=.false.
       end do
    end if

    inc_points = size(zpoints)

    if (itest>0) print *, '%xnc_allocate: Allocating xplasma_nclass arrays.'

    allocate(zbigtheta(nc_theta, inc_points))
    allocate(zgamma(inc_points))
    allocate(zdvol(inc_points))
    allocate(zdpsi(inc_points))
    allocate(zngb2(inc_points))
    allocate(zb2(inc_points))
    allocate(zbgrth(inc_points))

    allocate(zmom_a(nc_mom, inc_points), zmom_b(nc_mom, inc_points),&
         zmom_c(nc_mom, inc_points), zmom_d(nc_mom, inc_points),&
         zmom_f(nc_mom, inc_points))

    ! -- allocate working arrays --
    iztotal = max(8*max(64,nc_mom), max(2,nc_theta)*max(nc_rhob,inc_points))
                       ! account for large fft buffer

    allocate(zrho(iztotal), zchi(iztotal), zdet(iztotal))
    allocate(zans(iztotal,3), zgans(iztotal,2,1), ztens(3,3,iztotal))

    allocate(zintegrand1(nc_theta, nc_rhob),&
         zintegrand2(nc_theta,     nc_rhob),&
         zintegrand3(nc_theta,     nc_rhob),&
         zintegrand4(nc_theta,     nc_rhob) &
         )

    allocate(zresult1(nc_theta-1, inc_points),&
         zresult2(1,              inc_points),&
         zresult3(1,              inc_points),&
         zresult4(1,              inc_points) &
         )

  end subroutine xnc_allocate


  !
  ! ------------------------ xnc_get_id ------------------------
  ! returns the xplasma id for the named stored function,
  !
  subroutine xnc_get_id(s, name, itype, id, ier)

    type (xplasma), pointer :: s   ! object for which to compute profiles...

    character(*), intent(in)    :: name   ! name of xplasma id
    integer,      intent(in)    :: itype  ! 0 for axis, 1 for profile
    integer,      intent(inout) :: id     ! resulting id
    integer,      intent(out)   :: ier    ! error code
 
    !----------------------------------

    ier = 0

    if (itype == 0) then
       call xplasma_gridId(s, name, id)
    else
       call xplasma_profId(s, name, id)
    end if

    if (id <= 0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_id: Could not find the xplasma ID for item: '// &
            trim(name))
       ier=1
    end if
  end subroutine xnc_get_id

  !
  ! ------------------------- xnc_get_bigtheta ------------------------------------
  ! finds the function THETA(theta)
  ! computes gamma, dV/dxi, dpsi/dxi, <(n.grad(B))**2>, <B**2>
  !
  subroutine xnc_get_bigtheta(s,ier)

    type (xplasma), pointer :: s   ! object for which to compute profiles...

    integer, intent(out) :: ier   ! error flag

    !------------------------------

    integer :: itotal      ! number of points to evaluate Jacobian*B
    integer :: i,j,k       ! loop variables
    integer :: ierr        ! error return 
    integer :: iertmp,idi  ! tmp error code; 1d integrator ID
    integer :: idi_2d      ! 2d integrator ID
    integer :: iauth       ! author set flag
    integer :: idtmpf      ! f(chi,rho) temporary spline

    real(kind=r8) :: zs             ! scale factor
    real(kind=r8) :: zstart         ! theta where theta and THETA agree
    real(kind=r8) :: ndgb           ! n.grad(B)
    real(kind=r8) :: zadet          ! absolute value of jacobian

    real*8        :: d2v,d2va(1)    ! d2V/dxi2
    real*8 :: zrhominu

    !------------------------------

    ier  = 0    ! ok, bad variable naming
    ierr = 0

    iauth = 0

    ! --- compute integrands on original grid ---
    itotal = nc_theta*nc_rhob   ! number of points for evaluation

    k = 1
    do i = 1, nc_rhob
       do j = 1, nc_theta
          zrho(k) = zrhob(i)
          zchi(k) = ztheta_xp(j)
          k = k+1
       end do
    end do

    ! -- BR, BZ, BMOD --
    call xplasma_eval_prof(s, (/id_br, id_bz, id_bmod/), &
         xplasma_rho_coord, zrho(1:itotal), &
         xplasma_theta_coord, zchi(1:itotal), &
         zans(1:itotal,1:3),ierr)
    if (ierr /= 0) goto 100

    ! -- dBMOD/drho, dBMOD/dtheta --
    call xplasma_eval_prof(s, (/id_bmod, id_bmod/), &
         xplasma_rho_coord, zrho(1:itotal), &
         xplasma_theta_coord, zchi(1:itotal), &
         zgans(1:itotal,1:2,1), ierr, &
         ideriv1s=(/1,0/), ideriv2s=(/0,1/))
    if (ierr /= 0) goto 100

    ! -- J matrix, det=-J --
    call xplasma_rzjac(s, zrho(1:itotal), zchi(1:itotal), ierr, &
         r=ztens(3,3,1:itotal), rzdetj=zdet(1:itotal), &
         drdrho=  ztens(1,1,1:itotal) ,dzdrho=  ztens(1,2,1:itotal), &
         drdtheta=ztens(2,1,1:itotal) ,dzdtheta=ztens(2,2,1:itotal)  )

    zdet(1:itotal)=zdet(1:itotal)*ztens(3,3,1:itotal)
    if (ierr /= 0) goto 100
    
    ! -- fill integrands --
    k = 1
    zbmod_zero = 0._r8
    do i = 1, nc_rhob    ! loop over flux surfaces

       ! -- fill integrands for this flux surface --
       do j = 1, nc_theta
          br   = zans(k,1)
          bz   = zans(k,2)
          bmod = zans(k,3)

          db_drho   = zgans(k,1,1)
          db_dtheta = zgans(k,2,1)
          
          zadet= abs(zdet(k))

          if (zrhob(i)<=zrho_axis) then
             zintegrand3(j,i)=0._r8
             zbmod_zero = bmod
          else
             ! -- inverse jacobian matrix --
             jacob = ztens(1,1,k)*ztens(2,2,k)-ztens(1,2,k)*ztens(2,1,k)    ! 2d jacobian -- sqrt(g)/R, axisymmetry warning
          
             drho_dR   =  ztens(2,2,k)/jacob
             drho_dZ   = -ztens(2,1,k)/jacob  ! note ztens(2,1,k) = dR/dchi
             dtheta_dR = -ztens(1,2,k)/jacob  !      ztens(1,2,k) = dZ/drho
             dtheta_dZ =  ztens(1,1,k)/jacob
          
             db_dR = db_drho*drho_dR + db_dtheta*dtheta_dR
             db_dZ = db_drho*drho_dZ + db_dtheta*dtheta_dZ
          
             ndgb = (br*db_dR+bz*db_dZ)/bmod

             zintegrand3(j,i) = zadet*(ndgb**2)  ! Jacobian * (B.grad(B)/B)**2
          end if

          zintegrand2(j,i) = zadet             ! Jacobian
          zintegrand1(j,i) = zadet * bmod      ! Jacobian * BMOD 
          zintegrand4(j,i) = zadet * bmod**2   ! Jacobian * BMOD**2

          k = k+1
       end do
    end do

    ! --- initialize surface integrals for THETA ---

    ! some temporary items are created, then removed...
    iauth=1
    call xplasma_author_set(s,"xplasma_nclass",iertmp)  

    idi=0
    idi_2d=0
    idtmpf=0

    zrhominu = max(1.0e-11_r8,min(1.0e-6_r8,zrho_min))

    call xplasma_create_integ(s, "xplasma_nclass_1dinteg", &
         zpoints(1:nc_points), idi, ierr, &
         rhomin=zrhominu, cache_enable = .TRUE.)
    if (ierr /= 0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_bigtheta: Error creating 1d integrator.')
       ier=1 ; goto 50
    end if

    call xplasma_augment_integ(s, "xplasma_nclass_2dinteg", idi, &
         ztheta(1:nc_theta), idi_2d, ierr)
    if (ierr /= 0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_bigtheta: Error creating 2d integrator.')
       ier=1 ; goto 50
    end if

    !
    ! integration over surfaces -- integration is over both phi and theta
    !
    ! -- THETA integration --
    zresult1 = 0.0_r8     ! must be set for xplasma
    call xplasma_create_2dprof(s,"__TMPF_NCLASS", &
         id_chi, id_rho, zintegrand1, idtmpf, ierr, &
         ispline=2,ibcx2a=0,ibcx2b=0)
    if (ierr /= 0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_bigtheta: Error creating zintegrand1 tmp fcn.')
       ier=1 ; goto 50
    end if

    call xplasma_2d_zonint(s, idi_2d, "I(__TMPF_NCLASS)S", zresult1, ierr)
    if (ierr /= 0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_bigtheta: Error integrating zresult1.')
       ier=1 ; goto 50
    endif

    ! -- dV/dxi --
    zresult2 = 0.0_r8     ! must be set for xplasma
    call xplasma_create_2dprof(s,"__TMPF_NCLASS", &
         id_chi, id_rho, zintegrand2, idtmpf, ierr, &
         ispline=2,ibcx2a=0,ibcx2b=0)
    if (ierr /= 0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_bigtheta: Error creating zintegrand2 tmp fcn.')
       ier=1 ; goto 50
    end if

    call xplasma_rho_zonint(s, idi, "I(__TMPF_NCLASS)S", zresult2(1,:), ierr)
    if (ierr /= 0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_bigtheta: Error integrating zresult2.')
       ier=1 ; goto 50
    endif

    ! -- (n.grad(B))**2 --
    zresult3 = 0.0_r8     ! must be set for xplasma
    call xplasma_create_2dprof(s,"__TMPF_NCLASS", &
         id_chi, id_rho, zintegrand3, idtmpf, ierr, &
         ispline=2,ibcx2a=0,ibcx2b=0)
    if (ierr /= 0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_bigtheta: Error creating zintegrand3 tmp fcn.')
       ier=1 ; goto 50
    end if

    call xplasma_rho_zonint(s, idi, "I(__TMPF_NCLASS)S", zresult3(1,:), ierr)
    if (ierr /= 0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_bigtheta: Error integrating zresult3.')
       ier=1 ; goto 50
    endif

    ! -- B**2 --
    zresult4 = 0.0_r8     ! must be set for xplasma
    call xplasma_create_2dprof(s,"__TMPF_NCLASS", &
         id_chi, id_rho, zintegrand4, idtmpf, ierr, &
         ispline=2,ibcx2a=0,ibcx2b=0)
    if (ierr /= 0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_bigtheta: Error creating zintegrand4 tmp fcn.')
       ier=1 ; goto 50
    end if

    call xplasma_rho_zonint(s, idi, "I(__TMPF_NCLASS)S", zresult4(1,:), ierr)
    if (ierr /= 0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_bigtheta: Error integrating zresult4.')
       ier=1 ; goto 50
    endif

    ! --- normalize so THETA and theta agree at ztheta(1) ---
    zstart = ztheta(1)
    do k = 1, nc_points
       if (zpoints(k)<zrho_axis) then
          call xplasma_volume(s,(/ 0.0_r8 /),d2va,ierr, ideriv=2)
               ! array args for debugger
          d2v = d2va(1)
          if (ierr/=0) then
             call xplasma_errmsg_append(s, &
                  ' ? xnc_get_bigtheta: error calling xplasma_volume')
             ier=1 ; goto 50
          end if
          if (zbmod_zero<=0._r8) then
             call xplasma_errmsg_append(s, &
                  ' ? xnc_get_bigtheta: rho points do not include the axis???')
             ier=1 ; goto 50
          end if
          zgamma(k) = zbmod_zero*d2v     ! replace <B>dV/dxi with <B>d2V/dxi2
          zdvol(k)  = d2v                ! prevent divide by zero later on
          zngb2(k)  = 0._r8
          zb2(k)    = zbmod_zero*zbmod_zero
          zbigtheta(1:nc_theta,k) = ztheta(1:nc_theta)  ! not quite right but does it matter?
       else
          do j = 2, nc_theta-1
             zresult1(j,k)  = zresult1(j,k)+zresult1(j-1,k)    ! accumulate interval integrations
          end do
       
          ! -- 2pi comes from phi integration -- 
          zgamma(k)= zresult1(nc_theta-1,k)               ! (dV/dxi)*<B>, not gamma yet
          zdvol(k) = zresult2(1,k)                        ! (dV/dxi)
          zngb2(k) = max(0._r8,zresult3(1,k)/zdvol(k))    ! <(n.grad(B))**2>
          zb2(k)   = zresult4(1,k)/zdvol(k)               ! <B**2>
          
          zbigtheta(1,k) = zstart
          zs = 2.0_r8*znc_pi/zresult1(nc_theta-1,k) ! force THETA to have 2pi domain
          do j = 2, nc_theta
             zbigtheta(j,k) = zstart + zs*zresult1(j-1,k)
          end do
       end if
    end do

    ! --- get gamma and dpsi/dxi ---
    call xplasma_eval_prof(s,id_psi,zpoints,zdpsi,ierr, ideriv1=1)
    !print *, 'zdpsi=',zdpsi
    if (ierr /= 0) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_bigtheta: error fetching dpsi/dxi')
       ier=1 ; goto 50
    end if

    !print *, 'ichi_ccw=',ichi_ccw,'  nsnccwi=',nsnccwi
    do i = 1, nc_points
       if (zpoints(i)<zrho_axis) then
          call xplasma_eval_prof(s,id_psi,(/ 0.0_r8 /),zdpsi(i:i),ierr, &
               ideriv1=2)  ! replace dpsi/dxi with d2psi/dxi2, on axis
          if (ierr /= 0) then
             call xplasma_errmsg_append(s, &
                  '? xnc_get_bigtheta: error d2psi/dxi2')
             ier=1 ; goto 50
          end if
       end if
       zdpsi(i)  = (-ichi_ccw*nsnccwi)*zdpsi(i)
       zgamma(i) = 4.0_r8*znc_pi*znc_pi*zdpsi(i)/zgamma(i)
       !print *, 'zgamma(',i,')=',zgamma(i), '   zdpsi=',zdpsi(i)
    end do

    ! --- register curves ---
    if (is_rhob .and. itest>0) then
       call xplasma_create_2dprof(s,"BIGTHETA",id_chi,id_rho,zbigtheta, &
            id_bigtheta,ierr, &
            ispline=2,ibcx1a=0,ibcx1b=0)
       if(ierr.ne.0) then
          call xplasma_errmsg_append(s, &
               '? xnc_get_bigtheta: Error registering BIGTHETA.')
          ier=1; goto 50
       endif

       call xplasma_create_1dprof(s,"XNC_GAMMA",id_rho,zgamma,id_gamma,ierr, &
            ispline=1,ibca=0,ibcb=0)
       if (ierr /= 0) then
          call xplasma_errmsg_append(s, &
               '? xnc_get_bigtheta: Error registering XNC_GAMMA.')
          ier=1 ; goto 50
       end if
    else
       id_bigtheta = 0
       id_gamma    = 0
    end if

50  continue   ! there used to be a deallocation statment here
    if(idtmpf.gt.0) call xplasma_remove_item(s,idtmpf,iertmp)
    if(idi.gt.0) call xplasma_remove_item(s,idi,iertmp)
    if(idi_2d.gt.0) call xplasma_remove_item(s,idi_2d,iertmp)
    if(iauth.gt.0) call xplasma_author_clear(s,"xplasma_nclass",iertmp)  
    return

100 call xplasma_errmsg_append(s,'? xnc_get_bigtheta: Error in setup')
    ier=1
    goto 50
  end subroutine xnc_get_bigtheta


  !
  ! ------------------------- xnc_get_abcd ------------------------------------
  ! compute the Fm moments using an FFT after having computed the THETA(theta) function.
  ! The process is to use a spline of theta(THETA) to evaluate theta at a uniform set of
  ! THETA FFT points.  Find the integrand values at these points then perform the FFT
  ! to get the moments.
  !
  subroutine xnc_get_abcd(s,ier)

    type (xplasma), pointer :: s   ! object for which to compute profiles...

    integer, intent(out) :: ier    ! error flag

    integer, parameter :: ncroll=4 ! increase fft density to account for aliasing

    integer :: i,j,k               ! loop variable
    integer :: ilin                ! linearity return code
    integer :: ierr                ! error code
    integer :: nfft                ! size of fft to compute
    integer :: mwk                 ! size of work and output arrays

    real(kind=r8), allocatable, dimension(:)   :: big_uniform   ! uniformly spaced THETA
    real(kind=r8), allocatable, dimension(:,:) :: little_spline ! spline coefficients for theta(THETA)
    real(kind=r8), allocatable, dimension(:)   :: zwork
    real(kind=r8), allocatable, dimension(:)   :: little   ! theta(THETA)
    real(kind=r8), allocatable, dimension(:)   :: dlittle  ! dtheta(THETA)/dTHETA
    real(kind=r8), allocatable, dimension(:)   :: fft_a,fft_b ! arguments to fft coeficients
    real(kind=r8) :: b               ! THETA
    real(kind=r8) :: blast           ! last THETA in loop
    real(kind=r8) :: pi2             ! 2.*pi
    real(kind=r8) :: zmin,  zmax     ! allowed range for spline theta
    real(kind=r8) :: zbmin, zbmax    ! allowed range for spline THETA
    real(kind=r8) :: db              ! delta THETA in spline interval
    real(kind=r8) :: zs              ! scaling
    real(kind=r8) :: zadet           ! absolute value of jacobian

    real(kind=r8),    allocatable, dimension(:) :: st,ct  ! fft results

    ier = 0
    nfft = 64
    do while((ncroll*nc_mom)>(nfft/2))
       nfft = 2*nfft            ! get space for fft
    end do
    mwk = (nfft/2)+1

    if (nfft>size(zrho,1) .or. nc_theta<2) then
       call xplasma_errmsg_append(s, &
            '? xnc_get_abcd: Not enough points to support fft')
       ier=1 ; return
    end if

    ! --- allocation ---
    allocate(st(mwk), ct(mwk))
    allocate(big_uniform(nfft), little(nfft), dlittle(nfft))
    allocate(fft_a(nfft), fft_b(nfft))
    allocate(little_spline(4,nc_theta), zwork(nc_theta))

    do j = 1, nfft
       big_uniform(j) = (2.0_r8*znc_pi*(j-1))/nfft
    end do

    zmin  = ztheta(1)        - zmax_theta_error
    zmax  = ztheta(nc_theta) + zmax_theta_error
    pi2 = 2.0_r8*znc_pi
    do i = 1, nc_points
       ! -- form a nonuniform spline theta(THETA) --
       little_spline(1,:) = ztheta(:)
       call r8cspline(zbigtheta(:,i), nc_theta, little_spline, -1, 0.0_r8, -1, 0.0_r8, &
            zwork, nc_theta, ilin, ierr)
       if (ierr /= 0) then
          call xplasma_errmsg_append(s, &
               '? xnc_get_abcd:  Error forming spline of theta with BigTheta axis')
          call xplasma_errmsg_append(s, &
               '? xnc_get_abcd:  this could be due to a negative mod(B) or jacobian after interpolation')          
          ier=1 ; goto 50
       end if

       ! --- find values of theta,BR,BZ,BMOD,J at fft point of THETA ---
       zbmin = zbigtheta(1,i)        - zmax_theta_error
       zbmax = zbigtheta(nc_theta,i) + zmax_theta_error
       k = 1
       blast = zbmin-2.0_r8*pi2  ! force computation of new k
       do j = 1, nfft
          ! -- evaluate theta(THETA), dtheta(THETA) --
          b = big_uniform(j)

          ! -- bring b in range of zbigtheta --
          do while (b>zbmax)
             b = b-pi2
          end do
          do while (b<zbmin)
             b = b+pi2
          end do
             
          ! - find k so that zbigtheta(k,i)<b<zbigtheta(k+1,i) -
          !   k restricted to [1,nc_theta-1] interval
          if (abs(b-blast)>znc_pi) then
             k = int(nc_theta*(b-zbmin)/pi2)     ! guess spline interval
             k = max(1,min(nc_theta-1,k))
          end if
          do while (b<zbigtheta(k,i) .and. k>1)
             k = k-1
          end do
          do while (b>zbigtheta(k+1,i) .and. k<(nc_theta-1))
             k = k+1
          end do

          db = b - zbigtheta(k,i)
          little(j)  = little_spline(1,k) + db*(little_spline(2,k) + &
               db*(little_spline(3,k) + db*little_spline(4,k)))
          dlittle(j) = little_spline(2,k) + db*(2.0_r8*little_spline(3,k) + &
               3.0_r8*db*little_spline(4,k))

          ! -- bring little(j) in range --
          do while(little(j)>zmax)
             little(j) = little(j)-pi2
          end do
          do while(little(j)<zmin)
             little(j) = little(j)+pi2
          end do

          zrho(j) = zpoints(i)
          zchi(j) = little(j)
          blast = b
       end do

       ! -- BR, BZ, BMOD --
       call xplasma_eval_prof(s, (/id_br, id_bz, id_bmod/), &
            xplasma_rho_coord, zrho(1:nfft), &
            xplasma_theta_coord, zchi(1:nfft), &
            zans(1:nfft,1:3),ierr)
       if (ierr /= 0) then
          call xplasma_errmsg_append(s, &
               '? xnc_get_abcd:  Error in call to xplasma_eval_prof (1)')
          ier=1 ; goto 50
       end if
 
       ! -- dBMOD/drho, dBMOD/dtheta --
       call xplasma_eval_prof(s, (/id_bmod, id_bmod/), &
            xplasma_rho_coord, zrho(1:nfft), &
            xplasma_theta_coord, zchi(1:nfft), &
            zgans(1:nfft,1:2,1), ierr, &
            ideriv1s=(/1,0/), ideriv2s=(/0,1/))
       if (ierr /= 0) then
          call xplasma_errmsg_append(s, &
               '? xnc_get_abcd:  Error in call to xplasma_eval_prof (2)')
          ier=1 ; goto 50
       end if
       
       ! -- J matrix, -J --
       call xplasma_rzjac(s, zrho(1:nfft), zchi(1:nfft), ierr, &
            r=ztens(3,3,1:nfft), rzdetj=zdet(1:nfft), &
            drdrho=  ztens(1,1,1:nfft) ,dzdrho=  ztens(1,2,1:nfft), &
            drdtheta=ztens(2,1,1:nfft) ,dzdtheta=ztens(2,2,1:nfft)  )

       zdet(1:nfft)=zdet(1:nfft)*ztens(3,3,1:nfft)
       if (ierr /= 0) then
          call xplasma_errmsg_append(s, &
               '? xnc_get_abcd:  Error in call to xplasma_rzjac')
          ier=1 ; goto 50
       end if

       ! -- fill array for argument to fft --
       if (zpoints(i)<zrho_axis) then
          fft_b = 0._r8
          fft_a = 0._r8
       else
          do j=1, nfft
             br   = zans(j,1)
             bz   = zans(j,2)
             bmod = zans(j,3)
             zadet= abs(zdet(j))
             
             db_drho   = zgans(j,1,1)
             db_dtheta = zgans(j,2,1)
             
             ! -- inverse jacobian matrix --
             jacob = ztens(1,1,j)*ztens(2,2,j)-ztens(1,2,j)*ztens(2,1,j)    ! 2d jacobian -- sqrt(g)/R, axisymmetry warning
             
             drho_dR   =  ztens(2,2,j)/jacob
             drho_dZ   = -ztens(2,1,j)/jacob  ! note ztens(2,1,j) = dR/dchi
             dtheta_dR = -ztens(1,2,j)/jacob  !      ztens(1,2,j) = dZ/drho
             dtheta_dZ =  ztens(1,1,j)/jacob
             
             db_dR = db_drho*drho_dR + db_dtheta*dtheta_dR
             db_dZ = db_drho*drho_dZ + db_dtheta*dtheta_dZ
             
             fft_b(j) = (br*db_dR+bz*db_dZ)*zadet*dlittle(j)
             fft_a(j) = fft_b(j)/bmod
          end do
       end if

       ! get Am,Bm,Cm,Dm moments from FFT
       zs = (2.0_r8*znc_pi*znc_pi)/(nfft*zdvol(i))  ! common factor
       call r8fftsc(fft_a, nfft, st, ct, ierr)
       do k = 1, nc_mom
          zmom_a(k,i) = zs*st(k+1)
          zmom_c(k,i) = zs*ct(k+1)
       end do

       zs = zs * zgamma(i)
       call r8fftsc(fft_b, nfft, st, ct, ierr)
       do k = 1, nc_mom
          zmom_b(k,i) = zs*st(k+1)
          zmom_d(k,i) = zs*ct(k+1)
       end do

       ! finish off the Fm moment
       zs = zdvol(i)/(2.0_r8*znc_pi*znc_pi*zdpsi(i)*zb2(i))
       zbgrth(i) = 2.0_r8*znc_pi*znc_pi*zdpsi(i)/zdvol(i)
       !print *, 'zbgrth(',i,')=',zbgrth(i)
       do k = 1, nc_mom
          zmom_f(k,i) = zs*(zmom_a(k,i)*zmom_b(k,i)+zmom_c(k,i)*zmom_d(k,i))
       end do

    end do  ! end loop over points i

    ! --- deallocation ---
50  deallocate(st, ct)
    deallocate(big_uniform, little, dlittle)
    deallocate(fft_a, fft_b)
    deallocate(little_spline, zwork)

  end subroutine xnc_get_abcd

end module xplasma_nclass

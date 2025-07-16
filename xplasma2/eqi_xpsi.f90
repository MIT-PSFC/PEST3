subroutine eqi_xpsi(nR,nZ,zR,zZ,zg_bdy,zpsi_bdy,nsnccwb,nsnccwi, &
     bra,bza,bphia,psia,ierr)

  use xplasma_definitions
  use eqi_rzbox_module

  use ezspline
  use ezspline_obj

  !  AD HOC field and Psi(R,Z) extrapolation
  !  smooth, stable, fast, reasonable but not physical!

  !  mod DMC 19 Apr 2006 -- "imethod" argument deleted.
  !  DMC May 2006 -- adapted for Xplasma Vsn 2...

  implicit NONE

  integer, intent(in) :: nR,nZ         ! R & Z grid sizes
  real*8, intent(in) :: zR(nR),zZ(nZ)  ! R & Z grids
  real*8, intent(in) :: zg_bdy         ! R*|Bphi| at plasma bdy
  real*8, intent(in) :: zpsi_bdy       ! Psi at bdy

  integer, intent(in) :: nsnccwb       ! Bphi field sign
  integer, intent(in) :: nsnccwi       ! Jphi sign
  !  +1:  points counter-clockwise (ccw) in tokamak viewed from above


  !  output fields and Psi -- only at points beyond plasma boundary!
  !      BR         BZ         Bphi         Psi
  real*8 bra(nR,nZ),bza(nR,nZ),bphia(nR,nZ),psia(nR,nZ)


  integer ierr                         ! completion code, 0=normal

  !  compute the internal BR and BZ on a cartesian grid; then use an
  !  extrapolation method to get BR and BZ outside.
  !    also compute psi(poloidal)

  !  earlier versions extrapolated psi(R,Z) (hence "xpsi" name) but
  !  the current algorithm works with B directly -- 2d axisymmetric
  !  geometry.  method seems less error prone (dmc Jan 2000)

  !  the (R,Z) (rho,chi) binlinear inverse map is extended to the
  !  extrapolated pseudo-magnetic coordinate space beyond the plasma
  !  boundary

  !  the points beyond the boundary are sorted into extrapolated rho
  !  bins to ease interpolation from (rho,chi)->(R,Z) of BR and BZ

  !  2d axisymmetric code

  !  mod dmc 20 Aug 2003 -- imethod argument added --
  !  mod dmc 19 Apr 2006 -- imethod deleted -- ONLY the "ad hoc" method
  !                is available now.

  !  DELETED:
  !   imethod=1 -- use del.B = 0 extrapolation, with del x B = 0 also
  !                until current has to be introduced to stabilize and 
  !                smooth the result (traditional method).

  !                the result can be rough, because the extrapolation
  !                through the eddy-current laden vessel, and coil region
  !                of the tokamak, is exponentially unstable, so...

  !  USED IN ALL CASES:
  !   imethod=2 -- ad hoc method:  for rho>rhobdy, an ad hoc nested
  !                (rho,chi) coordinate system has been formed, a 
  !                pseudo-flux-coordinate system.  Lpol(rho) gives
  !                the poloidal path length of each rho surface:

  !                for Bpol:
  !                Bpol(rho,chi) = Bpol(rhobdy,chi)*Lpol(rhobdy)/Lpol(rho)
  !                is assumed, with direction tangent to the extrapolated
  !                rho surface.

  !                for Bphi, the vacuum value nsnccwb*g(rhobdy)/R is taken.

  !                for the flux function psi: psi is not just a function
  !                of rho for rho>rhobdy.  Instead, psi(rho) is integrated
  !                along the contour [R,Z](rho) @ R=Rmax for each rho surface,
  !                using dpsi = dR*r*Bpol (Bpol=Bz at Rmax of each flux
  !                surface, since the surface tangent is vertical there),
  !                forming psirmx(rho)
  !                Also, chirmx(rho) = chi at Rmax of each rho surface is
  !                computed, and, at rhobdy, f(chi) = R*Bpol/(Rmax*Bpol(Rmax))
  !                is computed.  Then,
  !                   psi(rho,chi) = 
  !                      psi(rhobdy) + (psirmx(rho)-psi(rhobdy))*f(chi_adj)
  !                where chi_adj is a phase adjusted chi based on chirmx(rho).

  !                this leads to a quick, well behaved extrapolation that
  !                has "reasonable" properties for modeling purposes, but
  !                of course is merely ad hoc and non-physical: del.B, 
  !                del x B do not satisfy expected properties, nor does
  !                Bpol = grad(psi)/R outside the boundary...

  !                The results can be plotted with xplasma test programs
  !                which start with a fixed flux bdy interior equilibrium,
  !                e.g. trxpl (accessing TRANSP results)-- dmc Aug 2003.
  !---------------------------------------------------------------------
  !  splines defined for region beyond plasma boundary (rho > rhobdy):

  type (ezspline1_r8) :: lpolspl,psispl,chispl,rbvspl

  !  vs. extrapolated rho coordinate:

  !  lpolspl(rho) = poloidal path length vs. rho
  !  psispl(rho) = psi @ R = Rmax, vs. rho
  !  chispl(rho) = chi (poloidal angle) @ R = Rmax, vs. rho

  !  vs. poloidal angle
  !  fbvspl(chi-chi[Rmax]) = 
  !          [R*Bp](chi;rho=rhobdy)/[Rmax*Bp(Rmax;rho=rhobdy)]

  !  psi(rho,chi) evaluation procedure:
  !    1: find chix = chi-chispl(rho); apply 2pi phase shift if needed;
  !       evaluate fbvspl(chix)
  !    2. evaluate psispl(rho)-psispl(rhobdy)
  !    3. psi = psispl(rhobdy)+fbvspl(chix)*[psispl(rho)-psispl(rhobdy)]

  !  Bpol(rho,chi) evaluation procedure:
  !    mod(Bpol)=Bpol(rhobdy,chi)*Lpolspl(rhobdy)/Lpolspl(rho)
  !    direction = tangent to extrapolated rho-surface
  !    sign = consistent with poloidal field induced by plasma current 
  !           (based on nsnccwi)

  integer isignb                    ! poloidal field sign factor 

  integer ix
  real*8 zchioff
  real*8 zlpol_bdy,zlpol,zlpolp,zRav,zdR,zBZav,zdpsi
  real*8 zbpol,znorm

  integer, parameter :: inchib = 201 ! no. of chi pts at bdy
  real*8, dimension(:), allocatable :: zrpsix,zbzpsix,zwk,zchx,zrho1
  real*8, dimension(:), allocatable :: zbra,zbza,zbphia,zpsia,zlpola
  real*8 zrbvals(inchib,3),zrhobdy(inchib),zf(inchib)

  !---------------------------------------------------------------------

  real*8 zrho,zchi,zZtarg,zphdum,zdum

  real*8, dimension(:), allocatable :: zRtmp,zZtmp,zrhotmp,zchitmp,zphitmp
  logical, dimension(:), allocatable :: ioutsid

  integer jvec,ii,jv,iermax

  real*8, dimension(:), allocatable :: zchia,zpsfac
  real*8, dimension(:,:), allocatable :: zrztan,zbpola

  !  the exterior psi needs a boundary condition to make the poloidal
  !  B field continuous at the boundary
  !------------------------------------------
  integer :: nrhox,id_rhox,id_r,id_br,id_bz,ilist(3)
  integer :: id_Rx,id_Zx,id_chix,id_map,iR,iZ,inc

  real*8, dimension(:), allocatable :: rhox

  real*8, dimension(:), pointer :: eqbuf
  integer, dimension(:), pointer :: idata

  integer :: lrhomap,lchimap

  real*8, parameter :: ZERO=0.0d0
  real*8, parameter :: CPI =3.1415926535897931D+00
  real*8, parameter :: C2PI=6.2831853071795862D+00
  !------------------------------------------

  ierr=0

  call xplasma_find_item(sp,'__RHOX',id_rhox,ierr)
  if(ierr.ne.0) return

  call xplasma_find_item(sp,'__THETAX',id_chix,ierr)
  if(ierr.ne.0) return

  call xplasma_find_item(sp,'__R_EXTRAP',id_Rx,ierr)
  if(ierr.ne.0) return

  call xplasma_find_item(sp,'__Z_EXTRAP',id_Zx,ierr)
  if(ierr.ne.0) return

  call xplasma_grid_size(sp,id_rhox,nrhox,ierr)
  if(ierr.ne.0) return

  call xplasma_common_ids(sp,ierr, id_R=id_R, id_BR=id_br, id_BZ=id_bz)
  if(ierr.ne.0) return

  allocate(rhox(nrhox))
  call xplasma_grid(sp,id_rhox,rhox,ierr)
  if(ierr.ne.0) then
     deallocate(rhox)
     return
  endif

  do

     !  create Lpol and psi splines

     call ezspline_init(lpolspl,nrhox,(/0, 0/),ierr)
     if(ierr.eq.0) then
        lpolspl%x1 = rhox
     else
        call ezspline_error(ierr)
        exit
     endif

     call ezspline_init(psispl,nrhox,(/0, 0/),ierr)
     if(ierr.eq.0) then
        psispl%x1 = rhox
     else
        call ezspline_error(ierr)
        exit
     endif

     call ezspline_init(chispl,nrhox,(/0, 0/),ierr)
     if(ierr.eq.0) then
        chispl%x1 = rhox
     else
        call ezspline_error(ierr)
        exit
     endif

     call ezspline_init(rbvspl,inchib,(/-1, -1/),ierr) ! periodic
     if(ierr.eq.0) then
        do ix=1,inchib
           rbvspl%x1(ix) = (ix-1)*c2pi/(inchib-1)
        enddo
     else
        call ezspline_error(ierr)
        exit
     endif

     zrhobdy=rhox(1)     ! vector of bdy rho values [inchib]

     allocate(zwk(nrhox),zchx(nrhox),zrho1(nrhox))
     allocate(zrpsix(nrhox),zbzpsix(nrhox))
     
     zrho1=rhox(1)       ! vector of bdy rho values [nrhox]

     !  poloidal path length Lpol(rho) of extrapolated surfaces...

     call lpol_gen(zwk,rhox)
     call ezspline_setup(lpolspl,zwk,ierr)
     if(ierr.ne.0) then
        call ezspline_error(ierr)
        exit
     endif

     !  extrapolated psi spline
     !    zwk(...) contains Lpol(x) at grid pts
     !    zchx(...) will contain chi value corresponding to Rmax at 
     !              each x surface

     zphdum=ZERO
     do ix=1,nrhox
        call eqi_rzbox(rhox(ix),zphdum, &
             .FALSE.,zdum,zdum,.TRUE.,zrpsix(ix),zchx(ix), &
             .FALSE.,zdum,zdum,.FALSE.,zdum,zdum, ierr)
        if(ierr.ne.0) exit
     enddo
     if(ierr.ne.0) exit

     !  Bz(rhobdy,chi@Rmax(rho)) evaluated here:

     call xplasma_eval_prof(sp,id_bz, &
          xplasma_rho_coord,zrho1,xplasma_theta_coord,zchx,zbzpsix,ierr)
     if(ierr.ne.0) exit

     !  chi@Rmax vs. rho spline; assure no branch cut...

     zchioff=0
     do ix=2,nrhox
        zchx(ix)=zchx(ix)+zchioff
        if(zchx(ix)-zchx(ix-1).gt.CPI) then
           zchioff=zchioff-C2PI
           zchx(ix)=zchx(ix)-C2PI
        else if(zchx(ix)-zchx(ix-1).lt.-CPI) then
           zchioff=zchioff+C2PI
           zchx(ix)=zchx(ix)+C2PI
        endif
     enddo

     call ezspline_setup(chispl,zchx,ierr)
     if(ierr.ne.0) then
        call ezspline_error(ierr)
        exit
     endif

     !  save sign factor...
     if(zbzpsix(1).lt.0) then 
        isignb=-1
     else
        isignb=1
     endif

     !  will use L(rhobdy)/L(rho)*mod(Bpol(rho,Rmax)) for psi spline...

     zlpol_bdy=zwk(1)
     zbzpsix=abs(zbzpsix)*zlpol_bdy     ! abs(Bz)*Lpol(rhobdy)
     zlpol=zlpol_bdy

     zwk(1)=zpsi_bdy  ! psi at bdy, input
     do ix=2,nrhox
        zlpolp=zlpol                ! Lpol at previous extrap. surface.
        zlpol=zwk(ix)               ! Lpol at next extrap. surface
        zRav=(zrpsix(ix)+zrpsix(ix-1))/2
        zBZav=(zbzpsix(ix)/zlpol+zbzpsix(ix-1)/zlpolp)/2
        zdR=(zrpsix(ix)-zrpsix(ix-1))
        zdpsi=zdR*zRav*zBZav
        zwk(ix)=zwk(ix-1)+zdpsi
     enddo

     call ezspline_setup(psispl,zwk,ierr)
     if(ierr.ne.0) then
        call ezspline_error(ierr)
        exit
     endif

     !  rbv spline: evaluate R(chi),BR(chi),BZ(chi) at desired points,
     !  at rhobdy

     ilist(1)=id_R; ilist(2)=id_BR; ilist(3)=id_BZ
     call xplasma_eval_prof(sp,ilist, &
          xplasma_rho_coord,zrhobdy, &
          xplasma_theta_coord,zchx(1)+rbvspl%x1(1:inchib), &
          zrbvals,ierr)
     if(ierr.ne.0) exit

     do ix=1,inchib
        zf(ix)=zrbvals(ix,1)*sqrt(zrbvals(ix,2)**2+zrbvals(ix,3)**2)
     enddo
     zf=zf/zf(1)                    ! normalized to dimensionless variation
     call ezspline_setup(rbvspl,zf,ierr)
     if(ierr.ne.0) then
        call ezspline_error(ierr)
        exit
     endif

     !  make sure rhomap is available

     call xplasma_ctrans(sp,ierr, &
          R_in=zR(nR/2), Z_in=zZ(nZ/2), &
          rho_out=zrho, theta_out=zchi, maptype=3)
     if(ierr.ne.0) exit

     call xplasma_find_item(sp,'__FASTMAP',id_map,ierr)
     if(ierr.ne.0) exit

     call xplasma_blackbox_retrieve(sp, id_map, ierr, &
          ia_ptr = idata, r8a_ptr = eqbuf)
     if(ierr.ne.0) exit

     lrhomap=idata(1)
     lchimap=idata(2)

     exit
  enddo

  if(allocated(zwk)) deallocate(zwk,zchx,zrho1,zrpsix,zbzpsix)

  if(ierr.ne.0) then
     call xplasma_errmsg_append(sp, &
          ' ?eqi_xpsi ad hoc B(R,Z) extrapolation setup failure.')
     call spl_free
     deallocate(rhox)
     nullify(eqbuf)
     nullify(idata)
     return
  endif

  allocate(zRtmp(nR*nZ),zZtmp(nR*nZ),zrhotmp(nR*nZ),zchitmp(nR*nZ))
  allocate(zphitmp(nR*nZ),ioutsid(nR*nZ))

  jvec=0
  ii=0
  do iz=1,nZ
     do ir=1,nR
        ii=ii+1
        inc=(iz-1)*nR+(ir-1)
        zrho=eqbuf(lrhomap+inc)
        zchi=eqbuf(lchimap+inc)
        if(zrho.le.rhox(1)) then
           ioutsid(ii)=.FALSE.

        else
           !  (R,Z) point is outside plasma bdy -- get extrapolated (rho,chi)
           !  coordinates only, for now

           ioutsid(ii)=.TRUE.

           jvec=jvec+1
           zRtmp(jvec)=zR(ir)
           zZtmp(jvec)=zZ(iz)
           zrhotmp(jvec)=zrho
           zchitmp(jvec)=zchi
           zphitmp(jvec)=ZERO

        endif
     enddo
  enddo

  !  evaluate exterior vector -- just rho,chi (external) coordinates for now

  allocate(zbra(jvec),zbza(jvec),zbphia(jvec),zpsia(jvec),zrho1(jvec))
  allocate(zlpola(jvec),zrztan(jvec,2),zbpola(jvec,2))
  allocate(zchia(jvec),zpsfac(jvec))

  !  evaluate ad hoc extrapolated Bpol

  iermax=0

  zrho1=rhox(1)

  ilist(1)=id_Rx; ilist(2)=id_Zx
  call xplasma_eval_prof(sp,ilist(1:2), &
       xplasma_rhox_coord,zrhotmp(1:jvec), &
       xplasma_thx_coord,zchitmp(1:jvec), &
       zrztan,ierr, ideriv2=1)
  iermax=max(iermax,ierr)

  ilist(1)=id_BR; ilist(2)=id_BZ
  call xplasma_eval_prof(sp,ilist(1:2), &
       xplasma_rho_coord,zrho1,xplasma_theta_coord,zchitmp(1:jvec), &
       zbpola,ierr)
  iermax=max(iermax,ierr)

  call ezspline_interp(lpolspl,jvec,zrhotmp,zlpola,ierr)
  iermax=max(iermax,ierr)

  do jv=1,jvec
     zbpol=sqrt(zbpola(jv,1)**2+zbpola(jv,2)**2)*(zlpol_bdy/zlpola(jv))
     znorm=sqrt(zrztan(jv,1)**2+zrztan(jv,2)**2)
     zbra(jv)=isignb*zrztan(jv,1)*zbpol/znorm
     zbza(jv)=isignb*zrztan(jv,2)*zbpol/znorm
     zbphia(jv)=nsnccwb*zg_bdy/zRtmp(jv)
  enddo

  !  evaluate ad hoc psi

  call ezspline_interp(psispl,jvec,zrhotmp,zpsia,ierr) ! psi @Rmax
  iermax=max(iermax,ierr)

  call ezspline_interp(chispl,jvec,zrhotmp,zchia,ierr)
  iermax=max(iermax,ierr)

  do jv=1,jvec
     zchia(jv)=zchitmp(jv)-zchia(jv) ! adjusted chi
     do                 ! put in range [0,2pi]
        if(zchia(jv).lt.ZERO) then
           zchia(jv)=zchia(jv)+C2PI
        else if(zchia(jv).gt.C2PI) then
           zchia(jv)=zchia(jv)-C2PI
        else
           exit
        endif
     enddo
  enddo

  call ezspline_interp(rbvspl,jvec,zchia,zpsfac,ierr)
  iermax=max(iermax,ierr)

  do jv=1,jvec
     zpsia(jv)=zpsi_bdy+(zpsia(jv)-zpsi_bdy)*zpsfac(jv)
  enddo

  ii=0
  jvec=0
  do iz=1,nz
     do ir=1,nR
        ii=ii+1
        if(ioutsid(ii)) then
           jvec=jvec+1
           bra(ir,iz)=zbra(jvec)
           bza(ir,iz)=zbza(jvec)
           bphia(ir,iz)=zbphia(jvec)
           psia(ir,iz)=zpsia(jvec)
        endif
     enddo                    ! iR
  enddo                             ! Z loop

  deallocate(zbra,zbza,zbphia,zpsia,zrho1)
  deallocate(zlpola,zrztan,zbpola)
  deallocate(zchia,zpsfac)

  call spl_free

  nullify(eqbuf)
  nullify(idata)

contains
  subroutine spl_free

    call ezspline_free(lpolspl,ierr)
    call ezspline_free(psispl,ierr)
    call ezspline_free(rbvspl,ierr)
    call ezspline_free(chispl,ierr)

  end subroutine spl_free

  subroutine lpol_gen(lpol_est,rho)
    real*8, dimension(:), intent(out) :: lpol_est
    real*8, dimension(:), intent(in) ::  rho

    integer :: i,ith,isize
    integer, parameter :: inth=500

    !------------------------------
    real*8 :: thvec(inth)
    real*8 :: rz(inth,2),zsum,zdr,zdz
    integer :: ict(10)

    type (xpeval) :: th_intrp, rho_intrp
    !------------------------------

    isize = size(rho)

    do i=1,inth
       thvec(i)=(i-1)*c2pi/(inth-1)
    enddo

    call xplasma_x_lookup(sp,id_chix, thvec, th_intrp,ierr)
    if(ierr.ne.0) then
       call xpeval_free(th_intrp)
       return
    endif

    do i = 1,isize
       call xplasma_x_lookup(sp,id_rhox, rho(i), rho_intrp,ierr)
       if(ierr.ne.0) exit

       call xplasma_eval_prof(sp,id_Rx,th_intrp,rho_intrp,rz(1:inth,1),ierr)
       if(ierr.ne.0) exit

       call xplasma_eval_prof(sp,id_Zx,th_intrp,rho_intrp,rz(1:inth,2),ierr)
       if(ierr.ne.0) exit

       zsum=0
       do ith=2,inth
          zdr = rz(ith,1)-rz(ith-1,1)
          zdz = rz(ith,2)-rz(ith-1,2)
          zsum = zsum + sqrt(zdr*zdr + zdz*zdz)
       enddo
               
       lpol_est(i)=zsum

    enddo

    call xpeval_free(rho_intrp)
    call xpeval_free(th_intrp)

  end subroutine lpol_gen

end subroutine eqi_xpsi

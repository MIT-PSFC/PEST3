subroutine eqi_geqdsk(lun_geqdsk,geqdsk_lbl, &
     Rmin,Rmax,Zmin,Zmax,zcur, &
     id_p,id_q,id_psi_in,nh,nv,nb, &
     ierr)
  
  !  write GEQDSK file from xplasma module

  !  **correction** dmc 9 Nov 2001:
  !  follow Lang Lao sign convention for G-EQDSK files:
  !  -- psi always increasing
  !  -- direction of current specified in pcur scalar *only*

  !-----------------------------------
  !  DMC Sep 2010: add code to set magnetic axis to match Psi(R,Z) min
  !
  !  Note: surely no "thread safety" here.  A shared xplasma2 pointer "sp"
  !  from the module "eqi_rzbox_module" is used.  The call to eqi_geq_axis
  !  could involve an excursion into a root finder which needs a module with
  !  data elements that use the SAVE attribute.
  !
  !  At present in the TRANSP/xplasma2 world this routine is called only
  !  from the xplasma_profs module, and this routine makes the only existing
  !  call to eqi_geq_axis.  So... arguments could be changed to pass down
  !  xplasma2 object references and establish thread safety, but... it has
  !  not been attempted as of today.  (DMC noted Jan 2011).
  !-----------------------------------

  use xplasma_definitions
  use eqi_rzbox_module
  implicit NONE

  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

  !  input:

  integer lun_geqdsk                ! logical unit where to write file

  character*48 geqdsk_lbl           ! 48 character label for GEQDSK file
  real*8 Rmin,Rmax                  ! (Rmin,Rmax) of Psi(R,Z)
  real*8 Zmin,Zmax                  ! (Zmin,Zmax) of Psi(R,Z)
  real*8 zcur                       ! plasma current (amps)

  !  [Rmin,Rmax]x[Zmin,Zmax] must contain the core plasma but not exceed the
  !  available (R,Z) grids.

  integer id_p                      ! xplasma id:  Pressure profile
  integer id_q                      ! xplasma id:  q profile
  integer id_psi_in                 ! xplasma id:  Psi (or 0 to use default)

  integer nh                        ! #of GEQDSK horizontal pts
  integer nv                        ! #of GEQDSK vertical pts

  !  note nh also controls the number of pts in the 1d profiles

  integer nb                        ! #of pts in bdy contour

  !  output:

  integer ierr                      ! completion code (0=OK, file written).

  !  local data arrays...

  real*8 psirz(nh,nv)               ! Psi(R,Z) to be written
  real*8 psi(nh)                    ! local psi grid
  real*8 fpol(nh)                   ! f = R*Bt vs. (psi)
  real*8 ffprime(nh)                ! f*f', ' means d/dpsi, vs. psi
  real*8 pres(nh)                   ! p vs. psi
  real*8 pprime(nh)                 ! p' vs. psi
  real*8 qpsi(nh)                   ! q vs. psi

  real*8 xpsi(nh)                   ! rho grid (for equal steps in psi)
  real*8 rgrid(nh)                  ! R grid for psi(R,Z)
  real*8 zgrid(nv)                  ! Z grid for psi(R,Z)
  real*8 zztmp(nh)                  ! array of Z(k) for vector eval

  real*8 thgrid(nb)
  real*8 xtmp(nb)
  real*8 rbdy(nb),zbdy(nb)          ! contour boundary
  real*8 rlim(nb),zlim(nb)          ! limiter boundary
  integer nblim                     ! #of limiter points retained

  real*8 :: psi0a                   ! Psi0 offset to machine axis (if avail.)

  !----------------------
  logical, parameter :: intrp_flag = .TRUE.
  !----------------------

  !  scalars

  real*8 rdim,zdim,rleft,zmid,rmaxis,zmaxis
  real*8 rmaxis_save,zmaxis_save ! for debug testing
  real*8 zpsimag,zpsibdy
  real*8 rcentr,bcentr
  real*8 pcur

  real*8 zrmin,zrmax,zzmin,zzmax

  !  dummies

  integer idum
  real*8 xdum


  !------------------------------------------------------
  !  work stuff

  integer :: nsnccwi,nsnccwb,id_g,id_psi,id_R,id_Z

  integer i,j,igot
  integer ifcns(4)
  real*8 zbuf(4)

  real*8 zbufa(nh,4)
  real*8 zbufd(nh,4)

  real*8 zbufb(nb,2)

  real*8 xsmall
  real*8 zdpsidx
  real*8 :: zpsi0,zpsia

  real*8, parameter :: ZERO = 0.0_R8
  real*8, parameter :: C2PI = 6.2831853071795862_R8
  !------------------------------------------------------

  idum=0
  xdum=ZERO

  !  computational region, as per GEQDSK specification

  rdim=Rmax-Rmin
  zdim=Zmax-Zmin
  rleft=Rmin
  zmid=0.5_R8*(Zmin+Zmax)

  call xplasma_common_ids(sp,ierr,id_g=id_g,id_psi=id_psi,id_R=id_R,id_Z=id_Z)
  if(ierr.ne.0) return

  if(id_g.eq.0) then
     ierr=9999
     call xplasma_errmsg_append(sp,' ?eqi_geqdsk: g(rho) not found.')
  endif

  if(id_psi_in.gt.0) then
     id_psi = id_psi_in
  endif

  if(id_psi.le.0) then
     ierr=9999
     call xplasma_errmsg_append(sp,' ?eqi_geqdsk: psi(rho) not found.')
  endif

  if(min(id_R,id_Z).eq.0) then
     ierr=9999
     call xplasma_errmsg_append(sp,' ?eqi_geqdsk: [R,Z](rho,theta) not found.')
  endif
  if(ierr.ne.0) return

  !  magnetic axis

  call xplasma_mag_axis(sp,ierr, raxis=rmaxis,zaxis=zmaxis)
  if(ierr.ne.0) return

  rmaxis_save = rmaxis ! for debug testing
  zmaxis_save = zmaxis

  !  pol. flux & axis & @ bdy -- zpsibdy > zpsimag as per xplasma & G-eqdsk
  !  convention

  call xplasma_psi_range(sp,zpsimag,zpsibdy)

  !  signed vacuum field-- at geometric center (in R) of LCFS

  call xplasma_global_info(sp, ierr, bphi_ccw=nsnccwb, jphi_ccw=nsnccwi)

  call xplasma_RZminmax_plasma(sp, zrmin,zrmax,zzmin,zzmax, ierr)
  if(ierr.ne.0) return

  rcentr=0.5*(zrmin+zrmax)
  bcentr=nsnccwb/rcentr   ! will multiply in bdy R*Bphi later

  !  signed total plasma current

  pcur=nsnccwi*abs(zcur)

  !  ok form profiles:  equispaced psi grid

  psi(1)=zpsimag
  psi(nh)=zpsibdy
  do i=2,nh-1
     psi(i)=abs(zpsimag)+(i-1)*(abs(zpsibdy)-abs(zpsimag))/(nh-1)
  enddo

  !  corresponding x grid

  call xplasma_rhopsi_find(sp,psi,xpsi,ierr)
  if(ierr.ne.0) return

  !  and for finite difference evals near the axis:

  xsmall=xpsi(1)+1.0e-4_R8*(xpsi(2)-xpsi(1))

  !  evaluate profiles & derivatives

  ifcns(1)=id_g
  ifcns(2)=id_p
  ifcns(3)=id_psi
  ifcns(4)=id_q

  !  fcn values

  call xplasma_eval_prof(sp,ifcns(1:4),xpsi,zbufa,ierr)
  if(ierr.ne.0) go to 990

  bcentr = zbufa(nh,1)*bcentr  ! multiply in R*(vacuum field) magnitude

  do i=1,nh
     fpol(i)=zbufa(i,1)*nsnccwb
     pres(i)=zbufa(i,2)
     qpsi(i)=zbufa(i,4)
  enddo

  !  derivatives -- off axis

  call xplasma_eval_prof(sp,ifcns(1:4),xpsi,zbufd,ierr, ideriv1=1)
  if(ierr.ne.0) go to 990

  !  df/dpsi = (df/dx)/(dpsi/dx)
  !   ...ok except @mag. axis where dpsi/dx -> 0

  do i=2,nh
     ffprime(i)=fpol(i)*zbufd(i,1)*nsnccwb/zbufd(i,3)
     pprime(i)=zbufd(i,2)/zbufd(i,3)
  enddo

  !  do finite difference estimate for the axis

  call xplasma_eval_prof(sp,ifcns,xsmall,zbuf,ierr)
  if(ierr.ne.0) go to 990

  zdpsidx=max(1.0e-8_R8*(psi(2)-psi(1))/(xpsi(2)-xpsi(1)), &
       (zbuf(3)-psi(1))/(xsmall-xpsi(1)))
  ffprime(1)=fpol(1)*nsnccwb* &
       ((zbuf(1)-zbufa(1,1))/(xsmall-xpsi(1)))/zdpsidx
  pprime(1)=((zbuf(2)-zbufa(1,2))/(xsmall-xpsi(1)))/zdpsidx

  !  OK... now get psi(R,Z)

  do i=1,nh
     rgrid(i)=Rmin+(i-1)*(Rmax-Rmin)/(nh-1)
  enddo

  do j=1,nv
     zgrid(j)=Zmin+(j-1)*(Zmax-Zmin)/(nv-1)
     zztmp=zgrid(j)
     call xplasma_eval_prof(sp,id_psi, &
          xplasma_R_coord,rgrid,xplasma_Z_coord,zztmp, &
          psirz(1:nh,j),ierr)
     if(ierr.ne.0) go to 990
     
     !  leave psi sign unchanged -- G-EQDSK convention is same as XPLASMA
     !  convention
  enddo

  call xplasma_get_psi0(sp,psi0a,igot,ierr)
  if(ierr.ne.0) go to 990

  !  adjust Psi(R,Z); restoring original EFIT offset; reset extrema as well

  psirz = psirz + psi0a
  zpsimag = zpsimag + psi0a
  zpsibdy = zpsibdy + psi0a

  !-------------------------------------------------------
  !  and now get the plasma boundary specification...

  ifcns(1)=id_R
  ifcns(2)=id_Z
  do i=1,nb
     xtmp(i)=1.0_R8
     thgrid(i)=(i-1)*C2PI/(nb-1)
  enddo

  !  evaluate all but last point; force exact periodicity

  call xplasma_eval_prof(sp,ifcns(1:2), &
       xplasma_theta_coord,thgrid(1:nb-1), xplasma_rho_coord,xtmp(1:nb-1), &
       zbufb(1:nb-1,1:2),ierr)
  if(ierr.ne.0) go to 990
  zbufb(nb,1)=zbufb(1,1)
  zbufb(nb,2)=zbufb(1,2)

  !  copy into output arrays

  rbdy(1:nb)=zbufb(1:nb,1)
  zbdy(1:nb)=zbufb(1:nb,2)

  !  and now, finally, the limiter...

  call xplasma_limcon(sp,rlim,zlim,igot,ierr, tol=ZERO)
  if(ierr.ne.0) go to 990

  nblim=nb

  !-------------------------------------------------------
  if(intrp_flag) then

     ! write out and read back from ascii file -- to get the
     ! effect of loss of precision due to ascii I/O in EQDSK format

     rmaxis = rmaxis_save
     zmaxis = zmaxis_save

     call eqi_geq_axis(nh,rgrid,nv,zgrid,psirz,nb,rbdy,zbdy, &
          rmaxis,zmaxis,zPsi0,zPsia,ierr)
     if(ierr.ne.0) then
        call xplasma_errmsg_append(sp, &
             ' ?Newton loop error in eqi_geq_axis, eqi_geqdsk (xplasma).')
        go to 990
     endif

     ! ** adjust Psi(R,Z) so that the true magnetic axis value is zpsimag

     psirz = psirz + (zpsimag - zPsi0)

     ! ** adjust boundary value per result of eqi_geq_axis evaluation
     !    (include the general adjustment to all of psirz(:,:) ...)

     zpsibdy = zPsia + (zpsimag - zPsi0)

  endif

  !  write the GEQDSK file...

2000 format(a48,3i4)
2001 format(5e16.9)
2002 format(2i5)

  !  from notes (General Atomic, Lang Lao 02/21/00)

  idum=0
  xdum=ZERO

  write(lun_geqdsk,2000) geqdsk_lbl,idum,nh,nv
  write(lun_geqdsk,2001) rdim,zdim,rcentr,rleft,zmid
  write(lun_geqdsk,2001) rmaxis,zmaxis,zpsimag,zpsibdy,bcentr
  write(lun_geqdsk,2001) pcur,zpsimag,xdum,rmaxis,xdum
  write(lun_geqdsk,2001) zmaxis,xdum,zpsibdy,xdum,xdum
  write(lun_geqdsk,2001) (fpol(i),i=1,nh)
  write(lun_geqdsk,2001) (pres(i),i=1,nh)
  write(lun_geqdsk,2001) (ffprime(i),i=1,nh)
  write(lun_geqdsk,2001) (pprime(i),i=1,nh)
  write(lun_geqdsk,2001) ((psirz(i,j),i=1,nh),j=1,nv)
  write(lun_geqdsk,2001) (qpsi(i),i=1,nh)
  write(lun_geqdsk,2002) nb,nblim
  write(lun_geqdsk,2001) (rbdy(i),zbdy(i),i=1,nb)
  write(lun_geqdsk,2001) (rlim(i),zlim(i),i=1,nblim)

  !  that is all

  ierr=0

  go to 1000

  !-----------------------------------
  !  error exit...

990 continue

  ierr=9999
  call xplasma_errmsg_append(sp,' ?data lookup error in eqi_geqdsk (xplasma).')

1000 continue
  return

end subroutine eqi_geqdsk

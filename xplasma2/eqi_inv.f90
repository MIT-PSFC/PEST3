subroutine eqi_inv(ivec,zR,zZ,zphi,reltol,rtol,zrho,zth,nregion,ierr)

  !  **vectorized**
  !  mapping:  (R,Z,phi) -> (rho,chi,phi)

  use xplasma_definitions
  use eqi_rzbox_module
  IMPLICIT NONE

  integer ivec
  REAL*8 zR(ivec),zZ(ivec),zphi(ivec) ! input (R,Z,phi) target
  REAL*8 reltol                       ! error tolerance (relative)
  REAL*8 rtol                         ! error tolerance, R,Z units.

  !  on input, (R,Z,phi) contains the
  !  initial guess.  a good initial guess improves performance, but is
  !  generally not required for accurate results.

  REAL*8 zrho(ivec),zth(ivec)       ! output (rho,theta,(phi unchanged))
  integer nregion(ivec)             ! input/output region code
  ! if fast map is in use: nregion(ivec)=3 means pt is outside RxZ rectangle...

  integer ierr                      ! output error code; 0=OK

  !  on output,
  !      nregion = 1 means:  rho_axis.le.rho.le.rho_bdy
  !      nregion = 2 means:  outside plasma but inside mapped region
  !         rho.gt.rho_bdy is returned at these points
  !      nregion = 3 means:  outside plasma and outside mapped region
  !         rho = rho_vac is returned at these points

  !   ierr .ne. zero means a serious error in the code
  !--------------------------------------------------------------

  integer i,jvec

  logical imask(ivec),iskipi,iskipo,iexact
  integer iregion(ivec)
  real*8 zR_use(ivec),zZ_use(ivec)
  real*8 zrho_use(ivec),zth_use(ivec),zphi_use(ivec),zdist(ivec)
  real*8 zth_bdy1(ivec),zphi_bdy1(ivec)
  real*8 :: zrho_bdy,zrho_plas,zrho_plus,zrho_minus
  real*8 :: rmin,rmax,zmin,zmax,zdnorm
  logical :: sol

  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: ONE = 1.0d0
  real*8, parameter :: EPS = 1.0d-15

  !--------------------------------------------------------------
  !  form vector of points that are not clearly outside...

  ierr=0

  call xplasma_RZminmax_extended(sp, rmin,rmax, zmin,zmax, ierr, sol=sol)
  if(ierr.ne.0) return

  zdnorm=max((rmax-rmin),(zmax-zmin))/2

  zrho_plas = ONE
  zrho_plus = zrho_plas+EPS
  zrho_minus = zrho_plas-EPS

  if(sol) then
     call eqi_extrap_rho_bdy(zrho_bdy)
  else
     zrho_bdy = ONE
  endif

  iskipi=.TRUE.
  iskipo=.FALSE.
  iexact=.TRUE.
  call eqi_bdfind(ivec,zR,zZ,zphi,zrho_bdy,iskipi,iskipo,iexact, &
       zth_use,zphi_use,zdist,ierr)

  if(ierr.ne.0) return

  do i=1,ivec

     if(zdist(i).ge.0.0d0) then
        imask(i)=.TRUE.
        zth(i)=zth_use(i)
        if(zdist(i).gt.0.0d0) then
           zrho(i)=zrho_bdy*(zdnorm+zdist(i))/zdnorm
           nregion(i)=3 ! beyond last surface
        else
           zrho(i)=zrho_bdy
           if(sol) then
              nregion(i)=2 ! at SOL bdy
           else
              nregion(i)=1 ! at plasma bdy
           endif
        endif

     else
        nregion(i)=1  ! may get corrected later...
        imask(i)=.FALSE.
        if(zrho(i).gt.zrho_plas) then
           zrho(i)=min(zrho_bdy,max(zrho_plus,zrho(i)))
        else
           zrho(i)=min(zrho_minus,zrho(i))
        endif

     endif
  enddo

  !  OK all with imask(i)=.FALSE. are in the mapped region and Newton
  !  search can proceed...

  jvec=0
  do i=1,ivec
     if(.not.imask(i)) then
        jvec=jvec+1
        zR_use(jvec)=zR(i)
        zZ_use(jvec)=zZ(i)
        zrho_use(jvec)=zrho(i)
        zth_use(jvec)=zth(i)
        zphi_use(jvec)=zphi(i)
     endif
  enddo

  !  But first, if there is a scrape off region, separate points across
  !  the boundary...

  if(sol) then

     iskipi=.TRUE.
     iskipo=.TRUE.
     iexact=.TRUE.
     call eqi_bdfind(jvec,zR_use,zZ_use,zphi_use,zrho_plas,iskipi,iskipo, &
          iexact,zth_bdy1,zphi_bdy1,zdist,ierr)
     if(ierr.ne.0) return

     do i=1,jvec
        if((zrho_use(i).gt.zrho_plas).and.(zdist(i).lt.ZERO)) then
           zrho_use(i)=zrho_minus
           zth_use(i)=zth_bdy1(i)
        else if((zrho_use(i).lt.zrho_plas).and.(zdist(i).gt.ZERO)) then
           zrho_use(i)=zrho_plus
           zth_use(i)=zth_bdy1(i)
        endif
     enddo
           
  endif

  !  inverse map evaluation by 2x2 Newton's method

  if(jvec.gt.0) then
     call eqi_xinv2d(jvec,zrho_use,zth_use,zphi_use,zR_use,zZ_use,iregion, &
          reltol,rtol,zrho_bdy, &
          ierr)

     jvec=0
     do i=1,ivec
        if(.not.imask(i)) then
           jvec=jvec+1
           zrho(i)=zrho_use(jvec)
           zth(i)=zth_use(jvec)
           nregion(i)=iregion(jvec)
        endif
     enddo
  endif

  !-----------------------

end subroutine eqi_inv

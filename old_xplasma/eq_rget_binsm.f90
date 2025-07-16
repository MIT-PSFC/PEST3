subroutine eq_rget_binsm(inrho,zrho,id,zdelta,iint,zval,ierr)
  !
  !  given an integrated profile function id,
  !
  !  return an integrated profile, smoothed onto the user's grid
  !  from the step-function defined by binning the density on the function's
  !  original grid.
  !
  implicit NONE
 
  integer, intent(in) :: inrho           ! target gridsize, must be .gt.1
  real*8, intent(in) :: zrho(inrho)      ! target grid-- strict ascending order
 
  integer, intent(in) :: id              ! function id of f(rho) profile
 
  real*8, intent(in) :: zdelta           ! smoothing parameter
                                         ! if <0, zdelta=delta(rho) assumed.
 
  integer, intent(in) :: iint            ! data type id
  !  iint=1 -- area integral profile being smoothed --
  !            Examples: current density profiles
  !  iint=2 -- volume integral profile being smoothed --
  !            Examples: density, energy density, ptcl/power/momentum source
  !                      functions, ...
  !
  !    for iint > 0, the data is binned on the original grid, but smoothed
  !    on the user grid.  The unsmoothed unintegrated profile (e.g. density
  !    profile) becomes a step function with step widths matching the original
  !    grid of the data.  It is smoothed using a triangular weighting function
  !    whose base half-width (smoothing parameter) is  zdelta.
  !
  !    If the smoothing parameter is too small, the step-
  !    function-like behavior of the binned density profile will remain
  !    evident in the result
  !
  !      zdelta =~ delta(rho)
  !
  !    is recommended, where delta(rho) is the stepsize of the function's
  !    underlying rho grid.  The recommended value is used if the passed
  !    value for zdelta is < 0.  If zdelta=0, no smoothing occurs.
  !    if zdelta > 0, smoothing occurs using the user specified zdelta.
  !
  !    note that the result is independent of the interpolating function
  !    xplasma has for the (integrated) profile function #id.
  !
  !  NOTE: smoothing integrated profiles directly, e.g. with eq_rget_sm,
  !    is hazardous, because, if d/dV of the smoothed result is evaluated
  !    near the axis it will generally have a pole singularity there.
 
  real*8, intent(out) :: zval(inrho)     ! results of smoothed interpolation
 
  !  zval(rho) is a weighted average of the original profile over the
  !            range [rho-zdelta,rho+zdelta] with a triangular weighting
  !            function whose weight goes to 0 at the interval enpoints.
 
  integer ierr  ! completion code 0=OK
 
!$r8real_input:  zrho,zdelta
!$r8real_output:  zval
 
  !--------------------------------------------
  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: bcav = 0.5d0  ! smoothing bc
  real*8 zsm,zrhomin,zrhomax
  integer ix,lunerr
 
  !  for integral-constrained smoothing
  real*8, dimension(:), allocatable :: zrhoi  ! user grid, extended/merged
  integer, dimension(:), allocatable :: imapu ! map back to user grid
  real*8, dimension(:), allocatable :: zvali  ! integrated profile at zrhoi
  real*8, dimension(:), allocatable :: zvol   ! associated total vol or area
 
  integer inumloc,inumrho,iuser,iorig
  real*8, dimension(:), allocatable :: zrholoc  ! grid for density profile
  real*8, dimension(:), allocatable :: zdens  ! density profile
 
  real*8 ztol,zdsum,zdens0
 
  !  for checking type of function for integral-constrained smoothing:
  !  must be f(rho) -- mag. coordinate.
 
  integer magdim,magfit,cyldim,cylfit,iwarn
  integer id_grid,ingrid,idum,ixp
  real*8, dimension(:), allocatable :: zrho_orig  ! orig. rho of f(rho)
  integer, dimension(:), allocatable :: imapo ! map back to original grid
 
  !--------------------------------------------

  call eqi_rget_binck(inrho,zrho,id,iint,ierr)
  if(ierr.ne.0) return

  !-------------------------------------------------------
  !  smoothing with integration...
 
  call eq_rholim(zrhomin,zrhomax)
  zsm=min(zdelta,(zrhomax-zrhomin)/4) ! sanity check on smoothing parameter
 
  !--------------------------------------------
 
  if(zsm.eq.ZERO) then
     call eq_rget_bin(inrho,zrho,id,iint,zval,ierr)
     return
  endif

  !--------------------------------------------

  call eq_fun(id,magdim,magfit,cyldim,cylfit,iwarn,ierr)
  if(ierr.ne.0) return
  if(magdim.ne.1) then
     call eq_get_lunerr(lunerr)
     write(lunerr,*) ' ?eq_rget_binsm: iint=',iint
     write(lunerr,*) '  integral-constrained smoothing available only'
     write(lunerr,*) '  for profiles of form f(rho).'
     ierr=1
     return
  endif
  ztol=1.d-8*(zrhomax-zrhomin)
 
  !  get original grid...
 
  call eq_fun_rhogrid(id,id_grid)
  call eq_ngrid(id_grid,ingrid)
 
  allocate(zrho_orig(ingrid))
  call eq_grid(id_grid,zrho_orig,ingrid,idum,ierr)
  if(ierr.ne.0) then
     deallocate(zrho_orig)
     return
  endif

  !  if delta < 0 reset smoothing to match grid spacing

  if(zsm.lt.0.0d0) then
     zsm=zrho_orig(ingrid)-zrho_orig(1) ! safe upper limit
     do ix=1,inrho-1
        zsm=min(zsm,(zrho_orig(ix+1)-zrho_orig(ix)))
     enddo
  endif

  !  merge user grid and data grid
 
  allocate(zrhoi(inrho+ingrid))
  allocate(imapu(inrho+ingrid),imapo(inrho+ingrid)); imapu=0; imapo=0
 
  inumrho=1
  zrhoi(1)=zrhomin
  imapo(1)=1
  iorig=2         ! 1st pt of original grid is always zrhomin
 
  iuser=1         ! 1st pt of user grid might also be at zrhomin
  if(zrho(1).lt.zrhomin+ztol) then
     imapu(1)=1
     iuser=2
  endif
 
  do
     if((iorig.gt.ingrid).and.(iuser.gt.inrho)) exit
     if(iorig.gt.ingrid) then
        ! iuser.le.inrho -- must be last pt...
        if(iuser.lt.inrho) then
           call eq_errmsg('?eq_rget_binsm: grid merge algorithm failure.')
           ierr=99
           deallocate(zrhoi,imapu,imapo)
           return
        endif
        imapu(inumrho)=inrho
        iuser=inrho+1
        cycle
     endif
     if(iuser.gt.inrho) then
        inumrho=inumrho+1
        imapo(inumrho)=iorig
        zrhoi(inumrho)=zrho_orig(iorig)
        iorig=iorig+1
        cycle
     endif
 
     ! pts from both grids available; chose the minimum.
 
     if(zrho(iuser).lt.zrho_orig(iorig)) then
        inumrho=inumrho+1
        imapu(inumrho)=iuser
        zrhoi(inumrho)=zrho(iuser)
        iuser=iuser+1
        if(zrho_orig(iorig).lt.zrho(iuser-1)+ztol) then
           imapo(inumrho)=iorig
           iorig=iorig+1      ! pts are nearly equal
        endif
     else
        inumrho=inumrho+1
        imapo(inumrho)=iorig
        zrhoi(inumrho)=zrho_orig(iorig)
        iorig=iorig+1
        if(zrho(iuser).lt.zrho_orig(iorig-1)+ztol) then
           imapu(inumrho)=iuser
           iuser=iuser+1
        endif
     endif
  enddo
 
  !  now: zrhoi(1:inumrho) contains ordered merger of
  !   zrho(1:inrho) and zrho_orig(1:inorig) with duplicates
  !   (equal to tolerance ztol) removed.
  !  and:
  !   if(imapo(k).gt.0) then zrhoi(k) = zrho_orig(imapo(k))
  !   if(imapu(k).gt.0) then zrhoi(k) = zrho(imapu(k))
  !     point back into the original two grids.
 
  allocate(zvali(inumrho),                &  ! integrated values &
       zvol(inumrho))                        ! integrated volume/area
 
  inumloc=inumrho-1                          ! grid & density profile
  allocate(zrholoc(inumloc),zdens(inumloc))
 
  do ix=1,inumloc
     zrholoc(ix)=(zrhoi(ix)+zrhoi(ix+1))/2
  enddo
 
  do
     !  fetch integrated profile on merged grid
 
     call eq_rgetf(inumrho,zrhoi,id,0,zvali,ierr)
     if(ierr.ne.0) exit
 
     !  fetching volume (area) elements
     !  fetch integrated volumes (areas)
 
     if(iint.eq.1) call eq_area(inumrho,zrhoi,0,zvol,ierr)
     if(iint.eq.2) call eq_volume(inumrho,zrhoi,0,zvol,ierr)
     if(ierr.ne.0) exit
 
     !  compute density profile -- binned as per original grid
 
     ixp=1
     do ix=1,inumloc
        if(imapo(ix+1).gt.0) then
           zdens0=(zvali(ix+1)-zvali(ixp))/(zvol(ix+1)-zvol(ixp))
           zdens(ixp:ix)=zdens0
           ixp=ix+1
        endif
     enddo
 
     !  apply smooth to density profile
 
     call r8_qksmooth(inumloc,zrholoc,zdens,zsm,ZERO,bcav,ierr)
     if(ierr.ne.0) exit
 
     !  integrate and map back to get answer
     if(imapu(1).eq.1) zval(1)=0
     zdsum=0
     do ix=1,inumloc
        zdsum=zdsum+zdens(ix)*(zvol(ix+1)-zvol(ix))
        if(imapu(ix+1).ne.0) then
           zval(imapu(ix+1))=zdsum
        endif
     enddo
     exit  ! normal exit
  enddo
 
  deallocate(zrholoc,zdens)
  deallocate(zvol,zvali)
  deallocate(zrhoi,zrho_orig)
  deallocate(imapo,imapu)
 
end subroutine eq_rget_binsm

subroutine eq_rget_bin(inrho,zrho,id,iint,zval,ierr)
  !
  !  given an integrated profile function id,
  !  return an integrated profile on the user's grid.
  !
  !  if the interpolation order is >1, just do an interpolation.
  !  if the interpolation order is =1 (piecewise linear), do an integration
  !  of the step function of the inferred density
  !
  implicit NONE

  integer, intent(in) :: inrho           ! target gridsize, must be .gt.1
  real*8, intent(in) :: zrho(inrho)      ! target grid-- strict ascending order

  integer, intent(in) :: id              ! function id of f(rho) profile

  integer, intent(in) :: iint            ! data type id
  !  iint=1 -- area integral profile being smoothed -- 
  !            Examples: current density profiles
  !  iint=2 -- volume integral profile being smoothed -- 
  !            Examples: density, energy density, ptcl/power/momentum source
  !                      functions, ...
  !
  real*8, intent(out) :: zval(inrho)     ! results of interpolation/integration

  integer ierr  ! completion code 0=OK

  !----------------------------------------
  real*8 :: zvolarho(inrho),zp,zvolp
  real*8, parameter :: ZERO = 0.0d0
  !----------------------------------------
  integer magdim,magfit,cyldim,cylfit,iwarn,inx,idum,ix,irho
  integer lunerr,id_grid
  real*8, dimension(:), allocatable :: zx,zvola,zfint
  !----------------------------------------

  call eqi_rget_binck(inrho,zrho,id,iint,ierr)
  if(ierr.ne.0) return

  call eq_fun(id,magdim,magfit,cyldim,cylfit,iwarn,ierr)
  if(ierr.ne.0) return

  if(magdim.ne.1) then
     call eq_get_lunerr(lunerr)
     write(lunerr,*) ' ?eq_rget_bin: iint=',iint
     write(lunerr,*) '  integral-binned lookup available only'
     write(lunerr,*) '  for profiles of form f(rho).'
     ierr=1
     return
  endif

  if(magfit.ge.1) then
     !  just do the interpolation of the underlying function...
     call eq_rgetf(inrho,zrho,id,0,zval,ierr)

  else
     !  get original funtion grid and associated volumes or areas
 
     call eq_fun_rhogrid(id,id_grid)
     call eq_ngrid(id_grid,inx)
 
     allocate(zx(inx),zvola(inx),zfint(inx))
     call eq_grid(id_grid,zx,inx,idum,ierr)
     call eq_rgetf(inx,zx,id,0,zfint,ierr)

     if(iint.eq.1) then
        call eq_area(inx,zx,0,zvola,ierr)
        call eq_area(inrho,zrho,0,zvolarho,ierr)
     else
        call eq_volume(inx,zx,0,zvola,ierr)
        call eq_volume(inrho,zrho,0,zvolarho,ierr)
     endif

     !  do the step function integration

     zval(1:2)=ZERO
     zvolp=ZERO

     ix=2
     irho=2
     if(zrho(1).gt.zx(1)) irho=1

     zp = (zfint(ix)-zfint(ix-1))/(zvola(ix)-zvola(ix-1))
     do
        if((zx(ix).lt.zrho(irho)).and.(ix.lt.inx)) then
           zval(irho)=zval(irho)+zp*(zvola(ix)-zvolp)
           zvolp=zvola(ix)
           ix=ix+1
           zp = (zfint(ix)-zfint(ix-1))/(zvola(ix)-zvola(ix-1))
        else
           zval(irho)=zval(irho)+zp*(zvolarho(irho)-zvolp)
           zvolp=zvolarho(irho)

           irho=irho+1
           if(irho.gt.inrho) exit
           zval(irho)=zval(irho-1)
        endif
     enddo

     deallocate(zx,zvola,zfint)

  endif
 
end subroutine eq_rget_bin

subroutine eqi_rget_binck(inrho,zrho,id,iint,ierr)

  ! shared error checking...

  implicit NONE
  integer,intent(in) :: inrho
  real*8,intent(in) :: zrho(inrho)
  integer,intent(in) :: id
  integer,intent(in) :: iint

  integer, intent(out) :: ierr

  !----------------------------------
  integer :: ix,lunerr
  real*8 :: ztol,zrhomax,zrhomin
  !----------------------------------
  call eq_get_lunerr(lunerr)
  ierr=0
 
  if((iint.lt.1).or.(iint.gt.2)) then
     write(lunerr,*) ' ?eq_rget_binck:  iint=1 or iint=2 expected, got: ', &
          'iint = ',iint
     ierr=1
  endif
  if(inrho.le.1) then
     write(lunerr,*) &
          ' ?eq_rget_binck: output grid must contain at least 2 pts.'
     ierr=1
  endif
  do ix=2,inrho
     if(zrho(ix).le.zrho(ix-1)) then
        write(lunerr,*) ' ?eq_rget_binck: zrho(...) not in ascending order.'
        ierr=1
        exit
     endif
  enddo
  if(ierr.ne.0) return
 
  call eq_rholim(zrhomin,zrhomax)
  ztol=1.d-8*(zrhomax-zrhomin)
  !  need extra sanity checks here...
 
  if(zrho(1).lt.zrhomin-ztol) then
     write(lunerr,*) ' ?eq_rget_binck: user zrho(1) < rho_min:'
     write(lunerr,*) '  zrho(1)=',zrho(1),'; rho_min=',zrhomin
     ierr=1
     return
  endif
  if(zrho(inrho).gt.zrhomax+ztol) then
     write(lunerr,*) ' ?eq_rget_binck: user zrho(inrho) > rho_max:'
     write(lunerr,*) '  zrho(inrho)=',zrho(inrho),'; rho_max=',zrhomax
     ierr=1
     return
  endif
  do ix=1,inrho-1
     if(zrho(ix).gt.zrho(ix+1)-2*ztol) then
        write(lunerr,*) ' ?eq_rget_binck: user grid pts too closely spaced:'
        write(lunerr,*) '  zrho(',ix,')=',zrho(ix)
        write(lunerr,*) '  zrho(',ix+1,')=',zrho(ix+1)
        ierr=1
        return
     endif
  enddo

end subroutine eqi_rget_binck

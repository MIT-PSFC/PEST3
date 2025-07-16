subroutine eq_rget_sm(inrho,zrho,id,zdelta,zval,ierr)
  !
  !  return a profile interpolated and then smoothed.
  !  caution: not recommended for volume- or area-integrated profiles; see
  !    eq_rget_binsm(...)
  !
  implicit NONE
 
  integer, intent(in) :: inrho           ! target gridsize, must be .gt.1
  real*8, intent(in) :: zrho(inrho)      ! target grid-- strict ascending order
 
  integer, intent(in) :: id              ! function id of f(rho) profile
 
  real*8, intent(in) :: zdelta           ! smoothing parameter
 
  !  zdelta defines the base half-width of a triangular weighting function
  !  used for smoothing.  The smoothed output profile at each rho is
  !  weighted average of the original function over the range
  !    (rho-zdelta,rho+zdelta) with weight 1/zdelta at rho linearly
  !  declining to weight 0 at the end points.  (area of weighting
  !  triangle is unity).
 
  real*8, intent(out) :: zval(inrho)     ! results of smoothed interpolation
 
  !  zval(rho) is a weighted average of the original profile over the
  !            range [rho-zdelta,rho+zdelta] with a triangular weighting
  !            function whose weight goes to 0 at the interval enpoints.
 
  integer ierr  ! completion code 0=OK
 
  !--------------------------------------------
  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: bcav = 0.5d0  ! smoothing bc
  real*8 zsm,zrhomin,zrhomax
  integer ix,lunerr
 
  !--------------------------------------------
  call eq_get_lunerr(lunerr)
  ierr=0
  if(inrho.le.1) then
     write(lunerr,*) ' ?eq_rget_sm: output grid must contain at least 2 pts.'
     ierr=1
  endif
  do ix=2,inrho
     if(zrho(ix).le.zrho(ix-1)) then
        write(lunerr,*) ' ?eq_rget_sm: zrho(...) not in ascending order.'
        ierr=1
        exit
     endif
  enddo
  if(ierr.ne.0) return
 
  if(zdelta.le.ZERO) then
     !  interpolate without smoothing
     call eq_rgetf(inrho,zrho,id,0,zval,ierr)
  else
     !  interpolate with smoothing
     call eq_rholim(zrhomin,zrhomax)
     zsm=min(zdelta,(zrhomax-zrhomin)/4) ! sanity check on smoothing parameter
 
     !  just interpolate and smooth the data (no integrals involved)
     call eq_rgetf(inrho,zrho,id,0,zval,ierr)
     if(ierr.ne.0) return
     call r8_qksmooth(inrho,zrho,zval,zsm,bcav,bcav,ierr)
  endif
 
end subroutine eq_rget_sm

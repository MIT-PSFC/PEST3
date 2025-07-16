subroutine f_upwind0(UPWIND0,zdiffus,zdravi,zveloc,zphip1,zphi)

  implicit NONE

  !  evaluate upwind adjustment
  !  on output:  zphip1=0.5, zphi=0.5 => no upwind adjustment
  !  on output:  zphip1=0.0, zphi=1.0 => full upwind adjustment (outflow)
  !  on output:  zphip1=1.0, zphi=0.0 => full upwind adjustment (inflow)
  !    intermediate adjustments possible; zphi+zphip1 = 1.0 always.

  !-------------------------------------------

  real*8, intent(in) :: UPWIND0   ! upwind parameter (dimensionless)

  !  roughly, this limits fluctuations df/f that can occur due to lack
  !  of upwind differencing; the smaller the value set, the stronger
  !  the upwind adjustment that will be made

  !  the following dimensional quanties can be cgs or mks but choose
  !  one or the other

  real*8, intent(in) :: zdiffus   ! diffusivity  (cm2/sec) or (m2/sec)  .ge. 0
  real*8, intent(in) :: zdravi    ! <1/dr> (1/cm) or (1/m) -- surface spacing
  real*8, intent(in) :: zveloc    ! flow velocity (cm/sec) or (m/sec)

  !  output

  real*8, intent(out) :: zphip1   ! weight on outward zone
  real*8, intent(out) :: zphi     ! weight on inward zone

  !---------------------------------------
  ! local:
  real*8 zalph
  !---------------------------------------

  if(zdiffus.eq.0.0d0) then

     !  special case: DIFFUSIVITY = ZERO

     if(zveloc.gt.0.0d0) then
        zphi=1.0d0
        zphip1=0.0

     else if(zveloc.lt.0.0d0) then
        zphi=0.0
        zphip1=1.0d0

     else
        zphi=0.5d0
        zphip1=0.5d0
     endif

     return

  endif

  if(zveloc.eq.0.0d0) then

     !  special case: VELOCITY = ZERO

     zphi=0.5d0
     zphip1=0.5d0

     return

  endif

  !  normal case:  both D and v are non-zero; look at ratio...
  !  normalizing factor UPWIND0 from TRANSP namelist...

  zalph = min(1.0d0, UPWIND0*zdiffus*zdravi/abs(zveloc) )

  if(zveloc.gt.0.0d0) then

     zphip1=zalph*0.5d0
     zphi=1.0d0-zphip1

  else

     zphi=zalph*0.5d0
     zphip1=1.0d0-zphi

  endif

  return
end subroutine f_upwind0

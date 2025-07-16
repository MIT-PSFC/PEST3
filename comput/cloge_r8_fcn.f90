real*8 function cloge_r8_fcn(znecm3i,zteevi,zzeffi,iwarn)

  ! compute the electron coulomb log in REAL*8 precision
  ! formula from TRANSP source/magcor/loglam.for
  ! this code 0d, no "trcom".

  ! warning flag set if any values are out of expected range:
  !   zzeffi < 1, zne13i < 1, zteevi < 1

  real*8, intent(in) :: znecm3i  ! input electron density, cm**-3
  real*8, intent(in) :: zteevi   ! input electron temperature, eV
  real*8, intent(in) :: zzeffi   ! input Zeff = sum(ni*Zi**2)/ne

  integer, intent(out) :: iwarn  ! 0: normal; 1: input range warning

  !-------------------------------------------
  real*8 :: zne13,zteev,zzeff
  real*8 :: zlamde
  real*8, parameter :: ONE = 1.0d0
  !-------------------------------------------

  iwarn = 0

  if(znecm3i.lt.ONE) then
     zne13=ONE
     iwarn=1
  else
     zne13=znecm3i
  endif

  if(zteevi.lt.ONE) then
     zteev=ONE
     iwarn=1
  else
     zteev=zteevi
  endif

  if(zzeffi.lt.ONE) then
     zzeff=ONE
     iwarn=1
  else
     zzeff=zzeffi
  endif

  ! here the actual computation...

  if(zteev.lt.(50.0d0)) then
     ZLAMDE=1.54d10*zteev*SQRT(zteev)/ &
          (zzeff*sqrt(zne13))
  else
     ZLAMDE=1.09d11*zteev/ &
          (zzeff*sqrt(zne13))
  endif

  cloge_r8_fcn = log(zlamde)

end function cloge_r8_fcn

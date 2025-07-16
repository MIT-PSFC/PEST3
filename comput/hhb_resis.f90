subroutine hhb_resis(zneCM3i,zteEVi,zzeffi, &
     deltai,lpolCM,philimWB,rho,dVdrhoCM3,xiotb,vacRBphiTCM, &
     eta_spitzer,eta_nc,iwarn)

  ! 0d CGS implementation of resistivity calculation
  ! This is HHB resistivity as used in TRANSP, magcor/resis.for

  ! Reference:S. P. Hirshman, R. J. Hawryluk, B. Birge,
  ! Nucl. Fusion 17, 611 (1977). ("HHB")).

  ! All information passed through arguments; no modules or COMMON
  ! REAL*8 precision

  ! implemented in TRANSP/NTCC comput library
  ! reference to real*8 function cloge_r8_fcn, also in comput library

  !----------------------------
  implicit NONE

  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  REAL*8, PARAMETER :: ZERO = 0.0_R8
  REAL*8, PARAMETER :: ONE  = 1.0_R8
  REAL*8, PARAMETER :: XMU  = 1.25664E-6_R8
  REAL*8, PARAMETER :: TWOPI= 6.2831853071795862_R8

  !----------------------------
  ! arguments:

  real*8 :: zneCM3i  ! electron density, cm**-3, iwarn=1 set if .le. 1 ptcl/cm3
  real*8 :: zteEVi   ! electron temperature, eV, iwarn=1 set if .le. 1 eV
  real*8 :: zZeffi   ! Zeff; iwarn=1 set if .le.1

  ! if iwarn=1 is set, a non-zero resistivity is still returned, but it is
  ! likely to be rubbish

  ! for the following, iwarn=2 is set and zero resistivity is returned, if
  ! any of these arguments are zero or negative

  real*8 :: deltai   ! inverse aspect ratio r/R
  real*8 :: LpolCM   ! poloidal path length, cm
  real*8 :: philimWb ! enclosed toroidal flux at plasma boundary, Wb
  real*8 :: rho      ! flux surface label: sqrt(PhiWB/PhilimWB) where
                     ! phiWB is the enclosed toroidal flux at the rho surface
  real*8 :: dVdrhoCM3    ! dVol/drho at the flux surface, cm**3
  real*8 :: xiotb    ! iota(bar) = 1/q = dPsi/dPhi
  real*8 :: vacRBphiTCM  ! toroidal field (R*B_phi) at vacuum boundary, T*cm

  ! output...

  real*8, intent(out) :: eta_spitzer  ! Spitzer resistivity, Ohm*cm
  real*8, intent(out) :: eta_nc   ! HHB neoclassical resistivity, Ohm*cm

  integer, intent(out) :: iwarn   ! warning flag, 0=normal

  !----------------------------
  !  local variables (several names taken from TRANSP's magcor/resis.for)

  real*8 :: zneCM3,zteEV,zZeff

  real*8 :: cloge_r8_fcn,cloge

  real*8 :: zgamm,zconst

  real*8 :: zbpav,zcuror

  real*8 :: zvstae,zcr,zxi,zd1m,zft,zfmn,zfnc,zte32i

  integer :: jwarn

  !----------------------------

  iwarn = 0
  eta_spitzer = ZERO
  eta_nc = ZERO

  ! check plasma parameter inputs

  if(znecm3i.lt.ONE) then
     zneCM3=ONE
     iwarn=1
  else
     zneCM3=znecm3i
  endif

  if(zteevi.lt.ONE) then
     zteEV=ONE
     iwarn=1
  else
     zteEV=zteevi
  endif

  if(zzeffi.lt.ONE) then
     zZeff=ONE
     iwarn=1
  else
     zZeff=zzeffi
  endif

  ! check field and geometry inputs

  if(deltai.le.ZERO) iwarn=2
  if(LpolCM.le.ZERO) iwarn=2
  if(PhilimWb.le.ZERO) iwarn=2
  if(rho.le.ZERO) iwarn=2

  if(dVdrhoCM3.le.ZERO) iwarn=2
  if(xiotb.le.ZERO) iwarn=2
  if(vacRBphiTCM.le.ZERO) iwarn=2

  if(iwarn.eq.2) return

  !----------------------------------
  ! OK-- get electron Coulomb log

  cloge = cloge_r8_fcn(zneCM3,zteEV,zZeff,jwarn)

  !----------------------------------
  ! Spitzer resistivity

  !  for traditional (HHB) form of Spitzer resistivity
  ZGAMM=ZZEFF*0.581_R8*(2.67_R8+ZZEFF)/(1.13_R8+ZZEFF)
  zconst=5.22E-3_R8

  zte32i = ONE/(zteEV*sqrt(zteEV))

  eta_spitzer = zconst*cloge*ZGAMM*zte32i

  !----------------------------------
  ! NC resistivity 

  ! use Bpol volume average
  ! Bpol = (1/R)*grad(Psi) = (1/R)*grad(rho)*dPsi/drho
  !    dPsi/drho = xiotb*rho*PhilimWb/pi
  ! <Bpol> = [2pi*int(dl*R/grad(rho)*Bpol]/dVdrho
  !        =  2*Lpol*xiotb*rho*PhilimWb/dVdrho

  zBPAV = PhilimWB*2.E4_R8*rho*xiotb*LpolCM/dVdrhoCM3

  ! this formulation uses
  ! mu0*Ip = 2pi*r*<Bpol>; zcuror = Ip/r = 2pi*<Bpol>/mu0; convert to Amps/cm

  zcuror = TWOPI*zBPAV/(100.0_R8*XMU)

  ! code transcribed from magcor/resis.for

  ZVSTAE=3.46E-9_R8*zneCM3*cloge*vacRBPHITCM/(zcuror*zteEV*zteEV*sqrt(deltai))
  ZCR=0.56_R8*(3.0_R8-ZZEFF)/((3.0_R8+ZZEFF)*ZZEFF)
  ZXI=0.58_R8+0.2_R8*ZZEFF
  ZD1M=ONE-deltai
  ZFT=ONE-ZD1M*ZD1M/(SQRT(ONE-deltai*deltai)*(ONE+1.46_R8*SQRT(deltai)))
  ZFMN=ZFT/(ONE+ZXI*ZVSTAE)
  ZFNC=ONE/((ONE-ZFMN)*(ONE-ZCR*ZFMN))

  eta_nc = ZFNC*eta_spitzer

end subroutine hhb_resis

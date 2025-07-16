!----------------------------------------------------------------------
!      1.4_r8    function initialization.
!.......................................................................
!
!     this SUBROUTINE pstwill set up various arrays of quantities needed
!     in setting up the matrix elements.  
!.......................................................................
!***********************************************************
!       version for numerically calculated equilibrium
!***********************************************************
!
      SUBROUTINE pstfunint
!.......................................................................
 USE i2mex_mod
 USE pstcom

 USE l21com

 USE l22com

 USE mtrik1

 USE temps

 USE comggf

 USE newmet
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 TOP
      real*8 TOPS
      real*8 BOTTOM
      real*8 BOTTM2
      real*8 PSQSUM
      real*8 BP2
      INTEGER NSURF
      real*8 ZHL
      real*8 ZHR
      real*8 ZDPSI
      real*8 B2SUM
      real*8 VTSUM
      real*8 SSUM
      real*8 SJPHI
      real*8 DSUBR
      real*8 DSUBEE
      real*8 BP2SUM
      real*8 RGRPS
      real*8 BETA3
      real*8 ARAT
      real*8 XRAD
      real*8 ELLIP
      real*8 ZZMIN
      real*8 ZZMAX
      real*8 ZXMIN
      real*8 ZXMAX
      real*8 XIN
      real*8 XOUT
      real*8 ZAVERA
      real*8 ZCUR
      real*8 ZFACT
      real*8 ZZ0

      real*8 the(nths0+1)
      real*8 ze(nosurf), zf(nosurf), zh(nosurf)
      integer ier
      INTEGER I


!
!
      psia(1:nosurf) = (psibig(1:nosurf)-psimin)  / twopi
!
!
!      to set the scaling exactly by q-edge adjust scale
      scqdp1(1) =  1.0_r8 

!
!aplet 8/10/96      if( imode >  1 ) go to 41
!
!      compute di and gam for use in the small singular soln.
!      and the volume averaged beta
!

      the = twopi*(/ (real(i-1,r8)/real(nths0,r8), i=1,nths0+1) /)

!      scale factor for b0 on axis
!
      call i2mex_getB0(b0SquareCentre, ier)
      b0SquareCentre = b0SquareCentre**2

!
! changed dpsi(1) to zdpsi: aplet  12.1_r8 .93
! adapted to non-equidistant mesh: aplet  5.4_r8 .93

!
      top =  0.0_r8 
      tops =  0.0_r8 
      bottom =  0.0_r8 
      bottm2 =  0.0_r8 

      psqsum =  0.0_r8 
      bp2    =  0.0_r8 
!
      do surf = 1, nosurf ! changed from 1..nosurf 10/10/99
         nsurf = surf
         zhl =  0._r8 
         if( surf >   1    )  &
              zhl = (psia(surf) - psia(surf-1))/ 2.0_r8 
         zhr =  0._r8 
         if( surf < nosurf)  &
              zhr = (psia(surf+1) - psia(surf))/ 2.0_r8 
         zdpsi = zhl + zhr      
         !
         CALL pstdsubi(nsurf,b2sum,vtsum,ssum,sjphi,dsubr,dsubee,bp2sum,rgrps)
         !
         top = top + pa(surf)*vtsum*zdpsi
         tops = tops + pa(surf) * ssum * zdpsi
         bottom = bottom + b2sum * zdpsi
         bottm2 = bottm2 + vtsum*zdpsi
         ajphi(surf) = - sjphi
!!$         di(surf) = dsubee
!!$         dr(surf) = dsubr
         rgradpsi(surf) = rgrps
         psqsum = psqsum + pa(surf)**2 * vtsum * zdpsi
         bp2 = bp2 + bp2sum * zdpsi
enddo

      call i2mex_getGlasserGreenJohnsonEFH(nosurf, psia, ze, zf, zh, ier)
      dr = ze + zf + zh**2
      di = ze + zf + zh - 0.25_r8


call pstToAxis1(psia, ajphi)
call pstToAxis1(psia, di)
call pstToAxis1(psia, dr)
call pstToAxis1(psia, rgradpsi)

!
! plasma parameters
!
      call i2mex_getSurface(plasmaSurface, ier)
      call i2mex_error(ier)
      call i2mex_getVolume(plasmaVolume, ier)
      call i2mex_error(ier)

      rMagnetic = xma

      call i2mex_getVolumeAveragedPressure(nths0+1, nosurf, the, psia, &
           & volumeAveragePressure, ier)
      call i2mex_error(ier)

      call i2mex_getVolumeAveragedBSquare(nths0+1, nosurf, the, psia, &
           & volumeAverageBSquare, ier)
      call i2mex_error(ier)
      
!
! betas...
!
      volumeAverageBeta =  100._r8  * two * &
           & volumeAveragePressure/volumeAverageBSquare

      beta3 =  100._r8 * two*(xma/r)**2 * top / (bottm2*b0SquareCentre)
!
! inverse aspect ratio...
!
      call pstshape(arat,xrad,ellip,triangularity &
           , zzmin, zzmax, zxmin, zxmax)
!
      xin = zxmin
      xout = zxmax
      minorRadius = (xout-xin)/ 2._r8 
      rCentre = (xout+xin)/ 2._r8 
      majorRadius = rCentre
      inverseAspectRatio = minorRadius / rCentre

!
! averaged radius...
!
       zavera = sqrt( abs(plasmaSurface/pye) )
!
! plasma elongation...
!
      elongation = plasmaSurface/(pye * minorRadius**2 )
!
      call i2mex_getPlasmaCurrent(nths0+1, nosurf, the, psia, &
           & totalToroidalCurrent, ier)
      call i2mex_error(ier)
!
! poloidal beta's...
!
      betaPoloidal =  4._r8  * twopi * tops / totalToroidalCurrent**2
!!$      print*,'old betaPoloidal=',betaPoloidal
!!$      call i2mex_getBetaPoloidal(nths0+1, nosurf, the, psia, &
!!$           & betaPoloidal, ier)
!!$      call i2mex_error(ier)
!!$      print*,'new betaPoloidal=',betaPoloidal
!
! toroidal beta...
!
      call i2mex_getBetaToroidal(nths0+1, nosurf, the, psia, &
           & betaToroidal, ier)
      call i2mex_error(ier)
      betaToroidal = betaToroidal * 100._r8
!
! beta*
!
      betaStar =  100._r8 *  2.0_r8  * sqrt(psqsum/bottm2)*(rCentre/r)**2/b0SquareCentre
!
! troyon's factor
!
! aplet 11/9/96      TroyonG = beta3 / totalToroidalCurrent
      call i2mex_getBetaN(nths0+1, nosurf, the, psia, &
           & TroyonG, ier)
      call i2mex_error(ier)
      
! aplet 11/9/96      gstar = betaStar / totalToroidalCurrent
      gStar = betaStar * 2*twopi* 0.1_r8  / (totalToroidalCurrent/ (minorRadius*b0SquareCentre))
!
! pressure peaking factor...
!
   pressurePeakingFactor =  0.0_r8 
      if(pa(1)  >    0._r8 )  pressurePeakingFactor = pa(1) * bottm2 / top
!
! inductance...
!
!     fourPiInductance2 = twopi * bpsq /(xma*sumj**2*dth*dpsi(1))
      fourPiInductance = 2*twopi*bp2 / (rCentre*totalToroidalCurrent**2)
!!$      print*,'old fourpiinductance=',fourpiinductance
!!$      call i2mex_getLi(nths0+1, nosurf, the, psia, &
!!$           & fourpiinductance, ier)
!!$      call i2mex_error(ier)
!!$      print*,'new fourpiinductance=',fourpiinductance
!

      if (imode == 1) then
      write(itty  ,8300) rCentre, minorRadius, xma, zxmin, zxmax, zzmin, zzmax &
                       , inverseAspectRatio, zavera,  &
                    plasmaSurface, elongation, ellip, triangularity, plasmaVolume
      write(outmod,8300) rCentre, minorRadius, xma, zxmin, zxmax, zzmin, zzmax &
                       , inverseAspectRatio, zavera,   &
                    plasmaSurface, elongation, ellip, triangularity, plasmaVolume
      write(itty  ,8301) totalToroidalCurrent, b0SquareCentre, betaPoloidal, betaToroidal, betaStar, &
                    beta3, TroyonG, gstar, volumeAveragePressure, volumeAverageBSquare, &
                    volumeAverageBeta, pressurePeakingFactor &
                  , fourPiInductance
      write(outmod,8301) totalToroidalCurrent, b0SquareCentre, betaPoloidal, betaToroidal, betaStar, &
                    beta3, TroyonG, gstar, volumeAveragePressure, volumeAverageBSquare, &
                    volumeAverageBeta, pressurePeakingFactor &
                  , fourPiInductance
      end if
!
!
 8300 format(/1x,'figures of merit:'// &
  " major radius              r0 = ",f13.6/ &
  " minor radius               a = ",f13.6/ &
  " magnetic axis             rm = ",f13.6/ &
  "                         rmin = ",f13.6/ &
  "                         rmax = ",f13.6/ &
  "                         zmin = ",f13.6/ &
  "                         zmax = ",f13.6/ &
  " inverse aspect ratio epsilon = ",f13.6/ &
  " averaged minor radius        = ",f13.6/ &
  " plasma surface               = ",f13.6/ &
  " elongation                 k = ",f13.6/ &
  " ellipticity                e = ",f13.6/ &
  " triangularity          delta = ",f13.6/ &
  " plasma volume                = ",f13.6/       )
!
!
 8301 format(   " total toroidal current     i = ",f13.6/ &
  " at geometric centre    b0**2 = ",f13.6/  &
 "                beta-poloidal = ",f13.6/  &
 "                beta-toroidal = ",f13.6/  &
 "                    beta-star = ",f13.6," %"/ &
  "                 beta (b-mid) = ",f13.6," %"/ &
  " troyon factor              g = ",f13.6/&
   "                       g-star = ",f13.6/&
   " vol. averaged pressure < p > = ",f13.6/&
   " vol. averaged b**2  < b**2 > = ",f13.6/&
   " vol. averaged beta    <beta> = ",f13.6," %"/&
   " p(0)/<p>                 ppf = ",f13.6/&
   " 4 pi internal inductance  li = ",f13.6/)
!
 8302 format(   " lundquist number           s = ",f13.2)
!
!       1.4_r8 .3  graph equilibrium properties.
!
      write(outmod,9000)
!
!
 9000 format(1x," original equilibrium profiles passed from",  &
    " mapping",//)
!
      write(outequ,8000) r
      do 10 surf = 1,nosurf
        write(outequ,8001) psia(surf), &
  qa(surf),qpa(surf),pa(surf), ppa(surf),  &
  di(surf), dr(surf), ajphi(surf) &
 ,rgradpsi(surf)
 10   continue
!
!-.1372E
 8000 format(1x,'r = ',f10.4/ &
  '   psia',8x,'q',10x,'qp',11x,'p',10x,'pp',9x,'d_i',9x,'d_r',4x,'<j.grad_phi>', &
1x,'r.grad psi')
!
 8001 format(1x,f9.5,8(1x,e11.4))
!
! aplet  5.5_r8 .94  create new file which can be used as input to
! chease (file expeq for experimental equilibrium profiles...
!
      zcur = totalToroidalCurrent
!
! rescaling of plasma with respect to average radius ~ major radius
!
      zfact = (zxmax + zxmin)/ 2._r8 
      zz0 = (zzmax + zzmin)/ 2._r8 
      zz0 = zz0/zfact
      zcur = zcur/zfact
!
 1001 format(i5)
 1002 format(2e18.8)
 1003 format(1e18.8)
!
      return
 7000 CALL psterrmes(outpst,'funint')
      return
      end



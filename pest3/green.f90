      SUBROUTINE pstgreen

 USE pstcom
 USE l22com
 implicit none 
 integer, parameter :: r8 = selected_real_kind(12,100)

!
       REAL*8   pn ! <green.f>
       REAL*8   pp ! <green.f>
       REAL*8   aleg0 ! <green.f>
       REAL*8   aleg1 ! <green.f>
       REAL*8   sqpi ! <green.f>
       REAL*8   pii ! <green.f>
       REAL*8   gamg ! <green.f>
       INTEGER  nloc ! <green.f>
       REAL*8   xs2 ! <green.f>
       REAL*8   xt2 ! <green.f>
       REAL*8   xp2 ! <green.f>
       REAL*8   xm2 ! <green.f>
       REAL*8   zm2 ! <green.f>
       REAL*8   r14 ! <green.f>
       REAL*8   r1sq ! <green.f>
       REAL*8   r1 ! <green.f>
       REAL*8   s ! <green.f>
       INTEGER  kloc ! <green.f>
       REAL*8   ak ! <green.f>
       REAL*8   ak02 ! <green.f>
       REAL*8   gg ! <green.f>
       REAL*8   aval1 ! <green.f>
       REAL*8   aval2 ! <green.f>
       REAL*8   aval3 ! <green.f>
       REAL*8   aval4 ! <green.f>
       REAL*8   aval5 ! <green.f>
       REAL*8   aval6 ! <green.f>
 REAL*8	 	 	 :: aval0
 COMMON /c2f90aval0/ aval0 
      pm =  0._r8 
      pn =  0._r8 
      pp =  0._r8 
      aleg0 =  0._r8 
      aleg1 =  0._r8 
!
      sqpi  = sqrt(pye)
      pii  = two / pye
      gamg  = sqpi
!     n is floating point in common
      nloc  = n +  0.1E0_r8  
      xs2  = xs * xs
      xt2  = xt * xt
      xp2  = xs2 + xt2
      xm2  = xt2 - xs2
      zm2  = ( zt - zs )**2
      r14  = xm2*xm2 + zm2*zm2 + two*xp2*zm2
      r1sq = sqrt( r14 )
      r1 = sqrt( r1sq )
      s  = (xp2 + zm2 )/r1sq
!
!     use upwards recurrence relations...
!
      CALL pstaleg ( s,nloc, pm,pn,pp, aleg0,aleg1 )
!
!      calculate gamma(half-n)
      kloc=0
      ak=zero
      if ( nloc  ==  0 )  go to 10
    5 kloc = kloc+1
      ak = float(kloc)
      ak02 = half-ak
      gamg = gamg / ak02
      if ( kloc  /=  nloc )  go to 5
!
!     now pp contains p(n+1) and pn contains p(n)
!
   10 gg  = -two * sqpi * gamg / r1
      bval  = -gg*pn
      aval1 = ( n*(xs2+xt2+zm2)*(xt2-xs2-zm2)+xt2*(xm2+zm2))*pn
      aval2 = two*xt*xs*(xm2-zm2)*pp
      aval3 = ztp*(aval1+aval2) / xt
      aval4 = ( two*n+one)*(xp2+zm2)*pn+four*xt*xs*pp
      aval5 = xtp*(zt-zs)*aval4
      aval6 =(aval3-aval5) / ( xt*r14 )
      aval = - xt2*aval6 * gg / twopi
      aval0 = ztp*(two*xs*(zm2-xm2)*aleg1 - xt*(xm2+zm2)*aleg0)
      aval0 = aval0 + xtp*(zt-zs)*(four*xt*xs*aleg1+(xp2+zm2)*aleg0)
      aval0 = -aval0*xt / (r14*r1)
      return
      end

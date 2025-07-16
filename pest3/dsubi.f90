!......................................................................
      SUBROUTINE pstdsubi(nsurf,b2sum,vtsum,ssum,sjphi,dsubr,dsubee,bp2sum,rgrps)
!......................................................................
!
!      purpose: to compute mercier's ideal criteria di and gam for use
!      in the small singular solution on surface - nsurf.
!      jm.  dec. 1982
!
 USE pstcom
 USE l22com
 USE newmet
 USE mtrik1
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER NSURF
      real*8 B2SUM
      real*8 VTSUM
      real*8 SSUM
      real*8 SJPHI
      real*8 DSUBR
      real*8 DSUBEE
      real*8 BP2SUM
      real*8 RGRPS
      real*8 G
      real*8 GPV
      real*8 Q
      real*8 QP
      real*8 PP
      real*8 G2
      real*8 R2G2
      real*8 R2GGP
      real*8 SUM1
      real*8 SUM2
      real*8 SUM3
      real*8 SUM4
      real*8 SUM5
      real*8 SUM6
      real*8 SUM7
      real*8 XV
      real*8 B2
      real*8 BP2
      real*8 ALPHA2
      real*8 WW1
      real*8 WW2
      real*8 WW3
      real*8 WW4
      real*8 VP


!
!      set surface quantities
!
 INTEGER	 	 	 :: t
!
      g = ga(nsurf)
      gpv = gpa(nsurf)
      q = qa(nsurf)
      qp = qpa(nsurf)
      pp = ppa(nsurf)
!
!      set disk address to read in the surface metric
!
      g2 = g * g
      r2g2 = r2 * g2
      r2ggp = r2 * g * gpv
      sum1 =  0.0_r8 
      sum2 =  0.0_r8 
      sum3 =  0.0_r8 
      sum4 =  0.0_r8 
      sum5 =  0.0_r8 
      sum6 =  0.0_r8 
      sum7 =  0.0_r8 
      rgrps =  0.0_r8 
      b2sum =  0.0_r8 
      vtsum =  0.0_r8 
      ssum =  0.0_r8 
      sjphi =  0.0_r8 
      bp2sum =  0.0_r8 
!
      do 10 t = 1, mth
      xv = sqrt(xsq(t,nsurf))
      b2 = (grpssq(t,nsurf) + r2g2)/xsq(t,nsurf)
      bp2 = grpssq(t,nsurf)/xsq(t,nsurf)
      alpha2 = qp*qp*grpssq(t,nsurf) / b2
      b2sum = b2sum + b2*xjacob(t,nsurf)*dth
      vtsum = vtsum + xjacob(t,nsurf)*dth
      ssum = ssum + xjacob(t,nsurf)/xv * dth
      sjphi = sjphi + ( xv*pp + r2ggp / xv) * xjacob(t,nsurf) * dth/xv
      bp2sum = bp2sum + bp2*xjacob(t,nsurf)*dth
!
      sum1 = sum1 + (xjacob(t,nsurf)/alpha2)/mth
      sum2 = sum2 + xjacob(t,nsurf) / mth
      rgrps = rgrps +  xjacob(t,nsurf)*grpssq(t,nsurf) /(mth*xv)
!     if(.not. lscale) go to 10

! aplet  24.4_r8 .97       if(pp  ==   0.0_r8 )go to 10
      sum3 = sum3 + xjacob(t,nsurf) / grpssq(t,nsurf)
      sum4 = sum4 + xjprym(t,nsurf)
      sum5 = sum5 + xjacob(t,nsurf) * xsq(t,nsurf) / grpssq(t,nsurf)
      sum6 = sum6 + xjacob(t,nsurf) / (xsq(t,nsurf)*grpssq(t,nsurf))
      sum7 = sum7 + xjacob(t,nsurf) * grpssq(t,nsurf) / xsq(t,nsurf)
   10 continue
!
!      introduce factor of 1/twopi to fix beta-poloidal.
      sjphi = sjphi / twopi
!
!     if(.not. lscale) return
!!$      di(nsurf) = - 0.25_r8 
! aplet  24.4_r8 .97      if(pp  ==   0.0_r8 ) return
      sum3 = sum3 / mth
      sum4 = sum4 / mth
      sum5 = sum5 / mth
      sum6 = sum6 / mth
      sum7 = sum7 / mth
!
      ww1 = one + sum7 /(r*g*q)
      ww2 = pp/(r*g*qp**2)
      ww3 = one + (r*g)**3 * sum6 / q
      ww4 = sum4 - pp * sum5
      vp = sum2
!
      dsubee= -( 0.5_r8  - pp*r*g*sum3/qp)**2 - ww2*q*ww3*ww4
      dsubr = -ww2*(q*ww3*ww4 - vp*ww3/ww1*(qp-two*r*g*pp*sum3) &
              -pp/(r*g)*(vp*ww3/ww1)**2)
      return
 7000 CALL psterrmes(outpst,'dsubi')
      end



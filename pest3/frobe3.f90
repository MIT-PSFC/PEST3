      SUBROUTINE pstfrobe3
!.....................................................
!
!      this SUBROUTINE pstcomputes matrix elements of
!      the coefficients
!      of the frobenius expansion around rational surfaces.
!
! aplet: last version of  26.1.94
!
!     new version of routine froben. performs the following
!     calculations:
! 1)  frobenius expansion of big, prescribed solution.
! 2)  source term formed by:
!      p^t \chi - k \xi --> alxi1e and alxi1o.
!      q^t [ \chi - g p \xi ] ) --> qdchie and qdchio.
!      ( m' - nq )[ \chi - g p \xi ] --> ddchie and ddchio.
!      - i ( g^{-1} \chi - p \xi ) --> gm1dco and gm1dce.
!      \xi --> x1frbo and x1frbe.
!      d_\psi xi --> x1fpro and x1fpre.
!      \chi --> c1frbo and c1frbe.
! for reference consult j. plasma phys., vol. 43, pp 291- 310._r8 
!
! nb. note that we work here with -( 0._r8 , 1._r8 ) * the q operator and
! -( 0._r8 , 1._r8 ) *chi, both being purely real for up-down symmetric
! plasmas.
!
!
 USE pstcom

 USE l22com

 USE newmet

 USE r33com

 USE comggf

 USE comfft
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER MS1
      INTEGER NRAT
      INTEGER I
      INTEGER JRATL
      INTEGER JRATR
      real*8 ZDPSI
      real*8 ZDPSI2
      real*8 ZDIFF
      INTEGER LM1
      INTEGER LM2
      INTEGER JINC
      INTEGER LS
      INTEGER JS
      INTEGER JJS
      real*8 QP
      real*8 QPN
      real*8 QPP
      real*8 QPPN
      real*8 ZGM1S
      real*8 ZRS
      real*8 ZMUS
      real*8 PB
      real*8 PB1
      real*8 PB2
      real*8 PS
      real*8 PS1
      real*8 PS2
      real*8 ZNORM
      real*8 ZNORM2
      INTEGER JJMAX1
      INTEGER JINC1
      INTEGER JJ1
      INTEGER JJ2
      INTEGER J1
      INTEGER J2
      INTEGER J
      INTEGER JMIN
      INTEGER JMAX
      INTEGER L1
      real*8 ZPBK
      real*8 ZER11
      real*8 ZER12
      real*8 ZER21
      real*8 ZER22
      real*8 ZTOT
      INTEGER J0
      real*8 ZD2
      real*8 ZD1
      real*8 ZD12
      real*8 ZD22
      real*8 ZD110
      real*8 ZD210
      real*8 ZPSING
      INTEGER JSIDE
      INTEGER JSUMIN
      INTEGER JSUMAX
      real*8 ZPARIT
      real*8 ZDLAYB
      real*8 ZDLAY2
      real*8 ZXS
      real*8 ZAXS
      real*8 ZLNX
      real*8 ZD
      real*8 ZARG
      real*8 ZY
      real*8 ZHPLUS
      real*8 ZDHPLS
      real*8 ALNQS
      real*8 ALNQ1
      INTEGER JS1
      INTEGER JS2
      INTEGER JINCR
      real*8 ZDLA10
      real*8 ZNEW
      real*8 ZOLD
      INTEGER JSUR1
      real*8 ZNEWOL
      real*8 ZLN
      real*8 ZPOWER
      real*8 ZEXPAO
      real*8 ZEXPAE
      real*8 ZEXPGO
      real*8 ZEXPGE
      real*8 ZERRAO
      real*8 ZERRAE
      real*8 ZERRGO
      real*8 ZERRGE
      real*8 ZXNEW
      real*8 ZXOLD
      real*8 ZABSX
      real*8 ZDX
      INTEGER J3
      INTEGER JLMIN
      INTEGER JLMAX



!
 COMPLEX*16	 	 	 :: zc
 COMPLEX*16	 	 	 :: zc1
 COMPLEX*16	 	 	 :: zc2
 COMPLEX*16	 	 	 :: zxibig
 COMPLEX*16	 	 	 :: zchbig
 COMPLEX*16	 	 	 :: zxsmal
 COMPLEX*16	 	 	 :: zcsmal
 COMPLEX*16	 	 	 :: zxadjb
 COMPLEX*16	 	 	 :: zcadjb
 COMPLEX*16	 	 	 :: zxadjs
 COMPLEX*16	 	 	 :: zcadjs
 COMPLEX*16	 	 	 :: z11
 COMPLEX*16	 	 	 :: z12
 COMPLEX*16	 	 	 :: z21
 COMPLEX*16	 	 	 :: z22
 COMPLEX*16	 	 	 :: zhchi
 COMPLEX*16	 	 	 :: zhxi
 COMPLEX*16	 	 	 :: zdhchi
 COMPLEX*16	 	 	 :: zdhxi
 COMPLEX*16	 	 	 :: zqs
 COMPLEX*16	 	 	 :: zqss
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: zrpr
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: zqpr
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: zgm1pr
 COMPLEX*16, DIMENSION(2,2) 	 :: zainv
 COMPLEX*16, DIMENSION(2,2) 	 :: zaa
 COMPLEX*16, DIMENSION(2,2) 	 :: zaaa
 COMPLEX*16, DIMENSION(2) 	 :: zay
 COMPLEX*16, DIMENSION(nfn) 	 :: zaoold
 COMPLEX*16, DIMENSION(nfn) 	 :: zaeold
 COMPLEX*16, DIMENSION(nfn) 	 :: zgoold
 COMPLEX*16, DIMENSION(nfn) 	 :: zgeold
 COMPLEX*16, DIMENSION(nfn) 	 :: zaonew
 COMPLEX*16, DIMENSION(nfn) 	 :: zaenew
 COMPLEX*16, DIMENSION(nfn) 	 :: zgonew
 COMPLEX*16, DIMENSION(nfn) 	 :: zgenew
 COMPLEX*16, DIMENSION(2,nfn,nfn) 	 :: qm
 COMPLEX*16, DIMENSION(2,nfn,nfn) 	 :: rm
 COMPLEX*16, DIMENSION(2,nfn,nfn) 	 :: gmm1
 REAL*8, DIMENSION(nfm) 	 :: zenr
 REAL*8, DIMENSION(nfm) 	 :: zeni
!
!
      ms1 = msing(ms+1)
!
! no of rational surface...
!
      nrat = 1
      do i = 1, ms
         nrat = nrat + ( msing(i+1)-msing(i) )*msub(i)
      end do
      nrat = nrat + msub(ms)/2
!
! no of previous rational surface
!
      jratl = 1
      if ( ms >  1 ) then
         do i = 1,ms-1
            jratl = jratl + ( msing(i+1)-msing(i) )*msub(i)
         end do
         jratl = jratl + msub(ms-1)/2
      end if
      !
      ! no of following rational surface
      !
      jratr = nusurf
      if ( ms < nosing ) then
         jratr = 1
         do i = 1,ms+1
            jratr = jratr + ( msing(i+1)-msing(i) )*msub(i)
         end do
         jratr = jratr - msub(ms)/2
      end if
!
      zdpsi   = psinew(nrat+1) - psinew(nrat  )
      zdpsi2  = psinew(nrat  ) - psinew(nrat-1)
      zdiff = abs(zdpsi2-zdpsi)
      if ( zdiff >   1.E-6_r8  ) then
      write(outmod,1001)
      write(itty,1001)
!
!
 1001 format(1x,'**** warning in frobe2: ' &
   ,'adjacent nodes not equidistant')
      end if
!
      lm1 = lmin(1) - 1
      jmax1 = lmax(1) - lmin(1) + 1
      lm2 = -(jmax2-1)/2 - 1
      if ( jmax2  ==  1 )  lm2 = lm1
      jinc = lm1 - lm2
      ls = n*qa(nrat)+ 0.1_r8 
      js = ls - lmin(1) + 1
      jjs = js + jinc
!
! skip frobenius coefficient computation if imode  >    1._r8 ..
!
      if( imode  >   1 ) go to 104
!
!.................................
!
! first run ...
!
      surf = nrat
      qp = qpa(surf)
      qpn = n*qp
!     qpp =  0.5_r8  * ( qpa(nrat+1)-qpa(nrat-1) )/zdpsi
!!$      qpp = (qa(nrat-1) -  2._r8 *qa(nrat) + qa(nrat+1))/zdpsi**2
      qpp = qppa(nrat)
      qppn = n*qpp
!
! operators...
! at left-hand side of rational surface...
!
      surf = nrat - 1
!!$  26.2.97     call thfft2
!call thefft
call pstthfft2
!
! store matrices for derivatives...
!
      rm(1,1:jmax1,1:jmax1) = ttr(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
      qm(1,1:jmax1,1:jmax1) = ttq(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
      gmm1(1,1:jmax1,1:jmax1) = gm1(1:jmax1,1:jmax1)
!
! at right-hand side of rational surface...
!
      surf = nrat + 1
!!$  26.2_r8 .97     call thfft2
!call thefft 
call pstthfft2
!
!
! store matrices for derivatives...
!
      rm(2,1:jmax1,1:jmax1) = ttr(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
      qm(2,1:jmax1,1:jmax1) = ttq(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
      gmm1(2,1:jmax1,1:jmax1) = gm1(1:jmax1,1:jmax1)
!
! at the rational surface....
!
      surf = nrat
!!$   26.2_r8 .97    call thfft2
!call thefft
call pstthfft2
!
      zqpr(1:jmax1,1:jmax1)   = &
 (   qm(2,1:jmax1,1:jmax1) -   qm(1,1:jmax1,1:jmax1) )/( 2.0_r8 *zdpsi)
      zrpr(1:jmax1,1:jmax1)   = &
 (   rm(2,1:jmax1,1:jmax1) -   rm(1,1:jmax1,1:jmax1) )/( 2.0_r8 *zdpsi)
      zgm1pr(1:jmax1,1:jmax1) = &
 ( gmm1(2,1:jmax1,1:jmax1) - gmm1(1,1:jmax1,1:jmax1) )/( 2.0_r8 *zdpsi)
!
 102  continue
!
      zgm1s = real( gm1(js,js)  )
      zqs   =  ttq(jjs,jjs)
	  zqss  = conjg( zqs ) 
      zrs   = real( ttr(jjs,jjs) )

!
      zmus = xmu(ms)
      pb = - 0.5_r8  - zmus
      pb1 = pb+ 1._r8 
      pb2 = pb +  2._r8 
      ps = - 1.0_r8 -pb
      ps1 = ps+ 1._r8 
      ps2 = ps+ 2._r8 
      alphb(ms) = pb
      alphs(ms) = ps
!
! cmatch = 2 mu f_0
! 
      cmatch(ms) =  2.0_r8 *twopi2*qpn**2*zmus/zgm1s
!
! norm factor rs 2 mu f
!
     rsnorm(ms) = cmatch(ms)* rgradpsi(nrat)
!
     beplay(ms) = betapol(nrat)
     epslay(ms) = aiaspect(nrat)
     aqplay(ms) = alqolp(nrat)
!
 101  continue
!
      write(outmod,1000)surf,qa(surf),qp,qpp,pb,ps,ms &
   ,cmatch(ms),nrat,psinew(nrat) &
   ,rsnorm(ms),aiaspect(nrat),betapol(nrat),alqolp(nrat)
      write(itty,1000)surf,qa(surf),qp,qpp,pb,ps,ms &
   ,cmatch(ms),nrat,psinew(nrat) &
   ,rsnorm(ms),aiaspect(nrat),betapol(nrat),alqolp(nrat)
!
!
 1000 format(/   ' rational surface at surface no.  ',i4 &
  /   '                              q = ',e17.9/ &
  ' 1st derivative of q         qp = ',e17.9/ &
  ' 2nd derivative of q        qpp = ',e17.9/ &
  ' big exponent                pb = ',e17.9/ &
  ' small exponent              ps = ',e17.9/ &
  ' matching coef.     cmatch(',i3,') = ', e17.9/ &
  ' rational surface   psinew(',i3,') = ', e17.9/ &
  ' norm factor 2 mu f <r grad psi>     ', e17.9/ &
  ' inverse aspect ratio              = ', f17.9/ &
  ' beta poloidal 2 p/B_theta^2       = ', f17.9/ &
  ' Lq/Lp                             = ', f17.9)
!
! check operators...
!
      if( lchkap ) then
!
      jjmax1= jmax1+ jinc
      jinc1 = 1 + jinc
      write(outmod,1030)
      do jj1 = jinc1,jjmax1
      write(outmod,1010) ( real(ttq(jj1,jj2)), &
                     jj2 = jinc1,jjmax1)
      end do
      write(outmod,1031)
      do jj1 = jinc1,jjmax1
      write(outmod,1010) (aimag(ttq(jj1,jj2)), &
                     jj2 = jinc1,jjmax1)
      end do
!
      write(outmod,1011)
      do jj1 = jinc1,jjmax1
      write(outmod,1010) ( real(ttr(jj1,jj2)), &
                     jj2 = jinc1,jjmax1)
      end do
      write(outmod,1012)
      do jj1 = jinc1,jjmax1
      write(outmod,1010) (aimag(ttr(jj1,jj2)), &
                     jj2 = jinc1,jjmax1)
      end do
!
      write(outmod,1020)
      do j1 = 1,jmax1
      write(outmod,1024) ( real(gm1(j1,j2)), &
                     j2 = 1,jmax1)
      end do
      write(outmod,1026)
      do j1 = 1,jmax1
      write(outmod,1024) (aimag(gm1(j1,j2)), &
                     j2 = 1,jmax1)
      end do
!
      write(outmod,1023)
      do j1 = 1,jmax1
      write(outmod,1024) ( real(geem(j1,j2)), &
                     j2 = 1,jmax1)
      end do
      write(outmod,1027)
      do j1 = 1,jmax1
      write(outmod,1024) (aimag(geem(j1,j2)), &
                     j2 = 1,jmax1)
      end do
!
      write(outmod,1021)
      do j1 = 1,jmax1
      write(outmod,1010) ( real(zqpr(j1,j2)), &
                     j2=1,jmax1 )
      end do
      write(outmod,1028)
      do j1 = 1,jmax1
      write(outmod,1010) (aimag(zqpr(j1,j2)), &
                     j2=1,jmax1 )
      end do
!
      write(outmod,1022)
      do j1 = 1,jmax1
      write(outmod,1010) ( real(zrpr(j1,j2)), &
                     j2=1,jmax1 )
      end do
      write(outmod,1029)
      do j1 = 1,jmax1
      write(outmod,1010) (aimag(zrpr(j1,j2)), &
                     j2=1,jmax1 )
      end do
!
      write(outmod,1025)
      do j1 = 1,jmax1
      write(outmod,1010) ( real(zgm1pr(j1,j2)), &
                     j2=1,jmax1 )
      end do
      write(outmod,1032)
      do j1 = 1,jmax1
      write(outmod,1010) (aimag(zgm1pr(j1,j2)), &
                     j2=1,jmax1 )
      end do
!
 1030 format(1x,'re ttq(jj1,jj2) :')
 1031 format(1x,'im ttq(jj1,jj2) :')
 1010 format(9(1x,e11.4) )
 1024 format(9(f12.6) )
 1011 format(1x,'re ttr(jj1,jj2) :')
 1012 format(1x,'im ttr(jj1,jj2) :')
 1020 format(1x,'re gm1(j1,j2) :')
 1026 format(1x,'im gm1(j1,j2) :')
 1023 format(1x,'re geem(j1,j2) :')
 1027 format(1x,'im geem(j1,j2) :')
 1021 format(1x,'re zqpr(j1,j2) :')
 1028 format(1x,'im zqpr(j1,j2) :')
 1022 format(1x,'re zrpr(j1,j2) :')
 1029 format(1x,'im zrpr(j1,j2) :')
 1025 format(1x,'re zgm1pr(j1,j2) :')
 1032 format(1x,'im zgm1pr(j1,j2) :')
!
      end if
!
! calculate coefficients of frobenius expansion...
!
      do 24 j1 = 1,jmax1
      xix0(j1,ms)  = ( 0._r8 , 0._r8 ) 
      chi0(j1,ms)  = ( 0._r8 , 0._r8 ) 
      xilog(j1,ms) = ( 0._r8 , 0._r8 ) 
      chlog(j1,ms) = ( 0._r8 , 0._r8 ) 
      xicon(j1,ms) = ( 0._r8 , 0._r8 ) 
      chcon(j1,ms) = ( 0._r8 , 0._r8 ) 
      xix1(j1,ms)  = ( 0._r8 , 0._r8 ) 
      chi1(j1,ms)  = ( 0._r8 , 0._r8 ) 
      xix2(j1,ms)  = ( 0._r8 , 0._r8 ) 
      chi2(j1,ms)  = ( 0._r8 , 0._r8 ) 
      xsml0(j1,ms) = ( 0._r8 , 0._r8 ) 
      csml0(j1,ms) = ( 0._r8 , 0._r8 ) 
      xsml1(j1,ms) = ( 0._r8 , 0._r8 ) 
      csml1(j1,ms) = ( 0._r8 , 0._r8 ) 
 24   continue
!
! lowest order big frobenius coefficients...
!
      zxibig = ( 1.0_r8 , 0.0_r8 ) 
      xix0(js,ms) = zxibig/cmatch(ms)
!
! calculate 2 versions of chi0 and average...
!
!  1.4_r8 .97 if( abs(pb+ 1._r8 ) <= tiny( 1._r8 )   .AND.  ltransf ) then
if( abs(pb+ 1._r8 ) <=  1.E-6_r8   .AND.  ltransf ) then
zc1 =  0._r8 
zc2 = zc1
else
      zc1 = ( zqs - qpn*pb )/zgm1s
!!$zc2=zc1
     zc2 = zrs/(zqss + qpn*(pb+ 1._r8 ) )
     if( abs(zc1-zc2)/abs(zc1) >    0.01_r8  ) &
   write(itty,*)' *** inaccurate zchbig: ',zc1,' ',zc2
end if
!
!!$      zchbig = zc1
zchbig = (zc1+zc2)/ 2._r8 

      chi0(js,ms) = zchbig/cmatch(ms)
!
! xix0 and chi0 for the small solution...
!
      zxsmal = ( 1.0_r8 , 0.0_r8 ) 

      zc1 = (zqs - qpn*ps)/zgm1s
!  1.4_r8 .97 if( abs(ps) <= tiny( 1._r8 )   .AND.  ltransf ) then
if( abs(ps) <=  1.E-6_r8   .AND.  ltransf ) then
zc2 = zc1
else
!!$zc2=zc1
     zc2 = zrs/(zqss + qpn*(ps+ 1._r8 ) )
     if( abs(zc1-zc2)/abs(zc1) >    0.01_r8  ) &
   write(itty,*)' *** inaccurate zcsmal: ',zc1,' ',zc2
end if
!
!!$      zcsmal = zc1
zcsmal = (zc1+zc2)/ 2._r8 

      xsml0(js,ms) = zxsmal
      csml0(js,ms) = zcsmal
!
! adjoint vectors...
!
! big...
!
      zxadjb = (zqss + qpn*(pb+ 1._r8 ))/(qpn*( 2.0_r8 *pb+ 1.0_r8 ))
      zc1 = zgm1s/(qpn*( 2.0_r8 *pb+ 1.0_r8 ) )
!  1.4_r8 .97 if( abs(pb+ 1._r8 ) <= tiny( 1._r8 )   .AND.  ltransf ) then
if( abs(pb+ 1._r8 ) <=  1.E-6_r8   .AND.  ltransf ) then
zc2 = zc1
else
!!$zc2=zc1
     zc2 = zxadjb*( zqs - qpn*pb )/zrs
     if( abs(zc1-zc2)/abs(zc1) >    0.01_r8  ) &
   write(itty,*)' *** inaccurate zcadjb: ',zc1,' ',zc2
end if
!
!!$      zcadjb = zc1
zcadjb = (zc1 + zc2)/ 2._r8 

!
! small...
!
      zxadjs = (zqss + qpn*(ps+ 1._r8 ))/(qpn*( 2.0_r8 *ps+ 1.0_r8 ))
      zc1 = zgm1s/(qpn*( 2.0_r8 *ps+ 1.0_r8 ) )
!  1.4_r8 .97 if( abs(ps) <= tiny( 1._r8 )   .AND.  ltransf ) then
if( abs(ps) <=  1.E-6_r8    .AND.  ltransf ) then
zxadjs =  0._r8 
zc2 = zc1
else
!!$zc2=zc1
     zc2 = zxadjs*( zqs - qpn*ps )/zrs
     if( abs(zc1-zc2)/abs(zc1) >    0.01_r8  ) &
   write(itty,*)' *** inaccurate zcadjs: ',zc1,' ',zc2
end if
!
!!$      zcadjs = zc1
zcadjs = (zc1 + zc2)/ 2._r8 

!
!  1.4_r8 .97      if (abs(pb+ 1.0_r8 ) <= tiny( 1._r8 ) ) then
if (abs(pb+ 1.0_r8 ) <=  1.E-6_r8 ) then
!
! zero beta case...
!
      lzerop = .true.
      if( xsmnus(ms) /=  0._r8  ) then
      write(itty   ,*)'frobe3 warning: xsmnus(ms) = ',xsmnus(ms), &
                 ' ne 0'
      write(outmod ,*)'frobe3 warning: xsmnus(ms) = ',xsmnus(ms), &
                 ' ne 0'
      end if
      if( xsplus(ms) /=  0._r8  ) then
      write(itty   ,*)'frobe3 warning: xsplus(ms) = ',xsplus(ms), &
                 ' ne 0'
      write(outmod ,*)'frobe3 warning: xsplus(ms) = ',xsplus(ms), &
                 ' ne 0'
      end if
!
!aplet  22.6_r8 .93      xsplus(ms) =  0.0_r8 
!
! logarithmic singularity: frobenius coefficient...
!
      xilog(js,ms) = ( 1.0_r8 , 0.0_r8 ) 
      chlog(js,ms) = zchbig - qpn/zgm1s
!
      zc = - zgm1s*( zrpr(js,js) - zchbig*( 2._r8 *real( zqpr(js,js) ) &
      + qppn) &
      + zchbig**2 *zgm1pr(js,js) )/qpn**2
!
      xilog(js,ms) = xilog(js,ms) * zc /cmatch(ms)
      chlog(js,ms) = chlog(js,ms) * zc /cmatch(ms)

      XiLogTerm(ms) = xilog(js,ms)

!
! constant term: frobenius coefficient...
!
      xicon(js,ms) = ( 1.0_r8 , 0.0_r8 ) 
      chcon(js,ms) = zchbig
!
      zc2 = - zc + ( &
              zqpr(js,js) - zchbig*zgm1pr(js,js) &
              )/qpn &
     + qpp/( 2.0_r8 *qp)
!
      xicon(js,ms) = xicon(js,ms) * zc2 /cmatch(ms)
      chcon(js,ms) = chcon(js,ms) * zc2 /cmatch(ms)
!
      else
!
! order pb + 1, non-resonant...
!
         do j = 1,2
            jmin = 1
            if( j == 2 ) jmin = js + 1
            jmax = js - 1
            if( j == 2 ) jmax = jmax1
            do j1 = jmin,jmax
               jj1 = j1 + jinc
               l1 = j1 + lm1
               xix1(j1,ms) = (- ttq(jj1,jjs)*zxibig + gm1(j1,js)*zchbig) &
                    /((l1 - ls)*(pb +  1.0_r8 )) /cmatch(ms)
               chi1(j1,ms) = (- ttr(jj1,jjs)*zxibig + conjg(ttq(jjs,jj1)) &
                    *zchbig) &
                    /((l1 - ls)*(pb +  1.0_r8 )) /cmatch(ms)
            end do
         end do
!
! order ps + 1, non-resonant...
!
         do j = 1,2
            jmin = 1
            if( j == 2 ) jmin = js + 1
            jmax = js - 1
            if( j == 2 ) jmax = jmax1
            do j1 = jmin,jmax
               jj1 = j1 + jinc
               l1 = j1 + lm1
               xsml1(j1,ms) = (- ttq(jj1,jjs)*zxsmal + gm1(j1,js)*zcsmal) &
                    /((l1 - ls)*(ps +  1.0_r8 ))
               csml1(j1,ms) = (- ttr(jj1,jjs)*zxsmal + conjg(ttq(jjs,jj1)) &
                    *zcsmal) &
                    /((l1 - ls)*(ps +  1.0_r8 ))
            end do
         end do
!
!
! aplet  22.6_r8 .94      if( zmus >    0.5_r8  ) then
      if( zmus >    0.45_r8  ) then
!
      zpbk = pb +  2.0_r8  
!
! order pb + 1, resonant...
!
!            construct [ a + inq'(pb+1) ]^{-1}...
!
      zainv(1,1) = - ( zxsmal*zxadjs/(pb+ 1._r8 -ps) &
              +   zxibig*zxadjb            )/qpn
      zainv(1,2) =   ( zxsmal*zcadjs/(pb+ 1._r8 -ps) &
              +   zxibig*zcadjb            )/qpn
      zainv(2,1) =   ( zcsmal*zxadjs/(pb+ 1._r8 -ps) &
              +   zchbig*zxadjb            )/qpn
      zainv(2,2) =   ( zcsmal*zcadjs/(pb+ 1._r8 -ps) &
              +   zchbig*zcadjb            )/qpn
!
! check...
!
      z11 = - zainv(1,1)*( - zqs + qpn*(pb+ 1._r8 ) ) &
       + zainv(1,2)*zrs
      z12 = + zainv(1,1)*zgm1s &
       + zainv(1,2)*(  zqss + qpn*(pb+ 2._r8 ) )
      z21 = + zainv(2,1)*( - zqs + qpn*(pb+ 1._r8 ) ) &
       + zainv(2,2)*zrs
      z22 = + zainv(2,1)*zgm1s &
       - zainv(2,2)*(  zqss + qpn*(pb+ 2._r8 ) )
      zer11 = abs(z11- 1._r8 )
      zer12 = abs(z12)
      zer21 = abs(z21)
      zer22 = abs(z22- 1._r8 )
      ztot = zer11 + zer12 + zer21 + zer22
      if( ztot >   0.001_r8  ) &
  write(itty,*)' *** inaccurate inversion: ',z11,' ',z12,' ' &
 ,z21,' ',z22
!
! sum of <m| a |m'> <m'| y_1>...
!
      zay(1) =  0.0_r8 
      zay(2) =  0.0_r8 
      do j = 1,2
         jmin = 1
         jmax = js -1
         if( j == 2) then
            jmin = js + 1
            jmax = jmax1
         end if
         do j1 = jmin,jmax
            jj1 = j1 + jinc
            l1 = j1 + lm1
            zay(1) = zay(1) - ttq(jjs,jj1)*xix1(j1,ms)  &
                 + gm1(js,j1) *chi1(j1,ms)
            zay(2) = zay(2) + ttr(jjs,jj1)*xix1(j1,ms)  &
                 - conjg( ttq(jj1,jjs) )*chi1(j1,ms)
         end do
      end do
!
! add contribution from higher order derivative
!
      zaa(1,1) =  - zqpr(js,js) + qppn*pb/ 2.0_r8 
      zaa(1,2) =  + zgm1pr(js,js)
      zaa(2,1) =  + zrpr(js,js)
      zaa(2,2) =  + conjg( zqpr(js,js) ) + qppn*(pb/ 2._r8  + 1._r8 )
!
      zay(1) = zay(1) + zaa(1,1)*xix0(js,ms) + zaa(1,2)*chi0(js,ms)
      zay(2) = zay(2) + zaa(2,1)*xix0(js,ms) - zaa(2,2)*chi0(js,ms)
!
! resonant frobenius coefficient...
!
      xix1(js,ms) = + zainv(1,1)*zay(1) - zainv(1,2)*zay(2)
      chi1(js,ms) = - zainv(2,1)*zay(1) - zainv(2,2)*zay(2)
!
! non-resonant order pb+ 2._r8 ..
!
      do j = 1,2
         jmin = 1
         if(j ==  2) jmin = js + 1
         jmax = js - 1
         if(j ==  2) jmax = jmax1
         do j1 = jmin,jmax
            jj1 = j1 + jinc
            l1 = j1 + lm1
            xix2(j1,ms) = - zqpr(j1,js)*xix0(js,ms)  &
                 + qpn*(pb+ 1._r8 )*xix1(j1,ms) &
                 + zgm1pr(j1,js)*chi0(js,ms)
            chi2(j1,ms) = - zrpr(j1,js)*xix0(js,ms) &
                 + conjg( zqpr(js,j1) )*chi0(js,ms)  &
                 + qpn*(pb+ 2._r8 )*chi1(j1,ms)
            do j2 = 1,jmax1
               jj2 = j2 + jinc
               xix2(j1,ms) = xix2(j1,ms) &
                    - ttq(jj1,jj2)*xix1(j2,ms) + gm1(j1,j2)*chi1(j2,ms)
               chi2(j1,ms) = chi2(j1,ms) &
                    - ttr(jj1,jj2)*xix1(j2,ms) + conjg(ttq(jj2,jj1)) &
                    *chi1(j2,ms)
            end do
            xix2(j1,ms) = xix2(j1,ms)/( (l1 - ls)*(pb+ 2._r8 ) )
            chi2(j1,ms) = chi2(j1,ms)/( (l1 - ls)*(pb+ 2._r8 ) )
         end do
      end do
!
      end if
!
      end if
!
      if( zmus >=    1.0_r8  ) then
      write(itty,*)  '***frobe2: zmus = ',zmus,' >= 1'
      write(outmod,*)'***frobe2: zmus = ',zmus,' >= 1'
      stop 'ERROR: Mercier index too large in frobe2'
      end if
!-----------------------------------------------------
      j0 = - lm1
      write(itty,1014) ls
      write(itty,1014) ls
      write(itty,1015)  real(zxsmal),  real(zxadjs), &
                   real(zxibig),  real(zxadjb)
      write(itty,1015)  real(zcsmal),  real(zcadjs), &
                   real(zchbig),  real(zcadjb)
      write(itty,1040)
      write(itty,1015) aimag(zxsmal), aimag(zxadjs), &
                  aimag(zxibig), aimag(zxadjb)
      write(itty,1015) aimag(zcsmal), aimag(zcadjs), &
                  aimag(zchbig), aimag(zcadjb)
      write(itty,1008)
      write(itty,1013)  real(xilog(js,ms)), aimag(xilog(js,ms)) &
                ,  real(xicon(js,ms)), aimag(xicon(js,ms))
      write(itty,1013)  real(chlog(js,ms)), aimag(chlog(js,ms)) &
                ,  real(chcon(js,ms)), aimag(chcon(js,ms))
!
      write(itty,1016) (j, j=lm1+j0,jmax1+lm1)
      write(itty,1017) ( real(xix1(j1,ms)), j1 = j0,jmax1)
      write(itty,1017) ( real(chi1(j1,ms)), j1 = j0,jmax1)
      write(itty,1019) (j, j=lm1+j0,jmax1+lm1)
      write(itty,1017) (aimag(xix1(j1,ms)), j1 = j0,jmax1)
      write(itty,1017) (aimag(chi1(j1,ms)), j1 = j0,jmax1)
      write(itty,1018) (j, j=lm1+j0,jmax1+lm1)
      write(itty,1017) ( real(xix2(j1,ms)), j1 = j0,jmax1)
      write(itty,1017) ( real(chi2(j1,ms)), j1 = j0,jmax1)
      write(itty,1009) (j, j=lm1+j0,jmax1+lm1)
      write(itty,1017) (aimag(xix2(j1,ms)), j1 = j0,jmax1)
      write(itty,1017) (aimag(chi2(j1,ms)), j1 = j0,jmax1)
!
      write(outmod,1014) ls
      write(outmod,1015)  real(zxsmal),  real(zxadjs), &
                   real(zxibig),  real(zxadjb)
      write(outmod,1015)  real(zcsmal),  real(zcadjs), &
                   real(zchbig),  real(zcadjb)
      write(outmod,1040)
      write(outmod,1015) aimag(zxsmal), aimag(zxadjs), &
                  aimag(zxibig), aimag(zxadjb)
      write(outmod,1015) aimag(zcsmal), aimag(zcadjs), &
                  aimag(zchbig), aimag(zcadjb)
      write(outmod,1008)
      write(outmod,1013)  real(xilog(js,ms)), aimag(xilog(js,ms)) &
                ,  real(xicon(js,ms)), aimag(xicon(js,ms))
      write(outmod,1013)  real(chlog(js,ms)), aimag(chlog(js,ms)) &
                ,  real(chcon(js,ms)), aimag(chcon(js,ms))
!
      write(outmod,1016) (j, j=lm1+j0,jmax1+lm1)
      write(outmod,1017) ( real(xix1(j1,ms)), j1 = j0,jmax1)
      write(outmod,1017) ( real(chi1(j1,ms)), j1 = j0,jmax1)
      write(outmod,1019) (j, j=lm1+j0,jmax1+lm1)
      write(outmod,1017) (aimag(xix1(j1,ms)), j1 = j0,jmax1)
      write(outmod,1017) (aimag(chi1(j1,ms)), j1 = j0,jmax1)
      write(outmod,1018) (j, j=lm1+j0,jmax1+lm1)
      write(outmod,1017) ( real(xix2(j1,ms)), j1 = j0,jmax1)
      write(outmod,1017) ( real(chi2(j1,ms)), j1 = j0,jmax1)
      write(outmod,1009) (j, j=lm1+j0,jmax1+lm1)
      write(outmod,1017) (aimag(xix2(j1,ms)), j1 = j0,jmax1)
      write(outmod,1017) (aimag(chi2(j1,ms)), j1 = j0,jmax1)
!
! aplet  11.5_r8 .94
!
!aplet  22.6_r8 .94      if( lsmnus   .AND.  xmu(ms) >   0.5_r8  ) then
      if( lsmnus ) then
!
! subtract off small solution from big expansion...
!
      if( xmu(ms) <  1._r8  ) xsmnus(ms) = - xix1(js,ms)
!
      write(outmod,1090) ms,xsmnus(ms)
      write(itty  ,1090) ms,xsmnus(ms)
!
 1090 format(/1x,'coefficient for subtraction of small solution',    &
        ' has been reset to:',           ' xsmnus(',i2,') = ',e15.8,/)
!
      end if
!
! aplet  11.5_r8 .94
!
!
 1008 format(/1x,'zero beta coefficients:' &
 /1x,10x,'re log',11x,'im log',9x,'re const',9x,                                   'im const')
 1013 format(4e17.8)
!
!
 1014 format(/1x,'frobenius eigenvector and their adjoint for '&
  ,'ls = ',i3//  1x,'re vector: small',1x,'re adj.v.: small', &
 1x,'re vector:   big',1x,'re adj.v.:   big')
!
!
 1040 format(/  1x,'im vector: small',1x,'im adj.v.: small',  1x,'im vector:   big',1x,'im adj.v.:   big')
!	 
 1015 format(4e17.8)
 1016 format(/1x,'order pb+1 re coefficients:'/1x,9i12)
 1019 format(/1x,'order pb+1 im coefficients:'/1x,9i12)
 1017 format(1x,9e12.4)
 1018 format(/1x,'order pb+2 re coefficients:'/1x,9i12)
 1009 format(/1x,'order pb+2 im coefficients:'/1x,9i12)
!
!.....
 104  continue
!.....
!
      return
      end



